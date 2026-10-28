#' Run SuSiE fine-mapping across all clumped loci in a GWAS
#'
#' Lead SNPs from plink --clump output are turned into fine-mapping regions by
#' \code{define_finemap_regions()}: a window of \code{window_kb} either side of
#' each lead, with overlapping windows merged into a single region (up to
#' \code{max_region_kb} wide) so that a signal is only fine-mapped once. For each
#' region an LD matrix is computed from the 1000 Genomes reference panel via
#' plink and \code{susieR::susie_rss} is run.
#'
#' Outputs written to \code{output_finemap_dir}:
#' \itemize{
#'   \item \verb{<CHR>_<BP>_finemap.tsv.gz} per region (BP is the most significant
#'     lead SNP in the region): all GWAS variants in the region with SuSiE
#'     Z-scores, posterior inclusion probabilities (\code{PIP}), credible-set
#'     membership (\code{CS}) and per-credible-set log Bayes factors
#'     \code{LBF_1}, \code{LBF_2}, ...
#'   \item \code{finemap_regions.tsv}: one row per region (bounds, leads, status,
#'     number of variants with LD, number of credible sets, SuSiE convergence).
#'   \item \code{finemap_credible_sets.tsv}: one row per credible set (size,
#'     coverage, purity, lead variant and its PIP, log Bayes factor).
#' }
#'
#' @param gwas filename or dataframe of a standardised GWAS
#' @param clumped_file plink --clump output file
#' @param ancestry ancestry code matching 1000 Genomes panel prefix (EUR, EAS, AFR, AMR, SAS)
#' @param default_n GWAS sample size when not inferrable from \code{gwas}; see
#' @param output_finemap_dir directory to write one
#'   \verb{<CHR>_<BP>_finemap.tsv.gz} per fine-mapping region.
#' @param completion_file path to a sentinel file written on successful completion
#'   (one line: number of fine-mapping regions). If NULL, defaults to
#'   \verb{<output_finemap_dir>/finemap_complete_<dir_basename>.txt}.
#' @param window_kb half-width of the fine-mapping window around each lead SNP in kb
#' @param max_causal maximum number of causal signals per locus (SuSiE L, default 10)
#' @param coverage credible set coverage (default 0.95)
#' @param min_abs_corr minimum absolute correlation for credible sets (default 0.5)
#' @param max_region_kb maximum width of a merged fine-mapping region in kb (default 3000)
#' @return invisibly, the combined per-SNP finemap data.table (LD-matched SNPs only)
#' @export
finemap_gwas <- function(gwas,
                         clumped_file,
                         ancestry,
                         default_n,
                         output_finemap_dir,
                         completion_file = NULL,
                         window_kb = 1000,
                         max_causal = 10,
                         coverage = 0.95,
                         min_abs_corr = 0.5,
                         max_region_kb = 3000) {

  gwas <- get_file_or_dataframe(gwas)
  data.table::setDT(gwas)

  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF", "N", "N_CASE", "N_CONTROL")
  for (col in intersect(numeric_cols, colnames(gwas))) {
    data.table::set(gwas, j = col, value = as.numeric(gwas[[col]]))
  }

  data.table::setkey(gwas, CHR, BP)

  has_n <- "N" %in% colnames(gwas)
  if (!has_n && (length(default_n) != 1L || is.na(default_n))) {
    stop("Fine-mapping requires GWAS sample size")
  }

  if (!dir.exists(output_finemap_dir)) {
    dir.create(output_finemap_dir, recursive = TRUE)
  }

  if (is.null(completion_file)) {
    completion_file <- finemap_completion_file(output_finemap_dir)
  }

  lead_snps <- read_clumped_leads(clumped_file)
  regions <- define_finemap_regions(lead_snps, window_bp = window_kb * 1000, max_region_bp = max_region_kb * 1000)
  regions_file <- file.path(output_finemap_dir, "finemap_regions.tsv")
  credible_sets_file <- file.path(output_finemap_dir, "finemap_credible_sets.tsv")

  if (nrow(regions) == 0) {
    message("No clumped SNPs found; no per-locus finemap files written.")
    data.table::fwrite(add_region_status_columns(regions), regions_file, sep = "\t")
    data.table::fwrite(empty_credible_set_summary(), credible_sets_file, sep = "\t")
    write_finemap_complete_marker(completion_file, 0L, output_finemap_dir)
    empty <- data.table::data.table(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), PIP = numeric(), CS = integer()
    )
    return(invisible(empty))
  }

  has_rsid <- "RSID" %in% colnames(gwas)

  n_loci <- nrow(regions)
  message(paste(nrow(lead_snps), "clumped lead SNP(s) grouped into", n_loci, "fine-mapping region(s)"))

  # Pre-extract per-region subsets so workers don't need the full GWAS
  locus_subsets <- vector("list", n_loci)
  window_subsets <- vector("list", n_loci)
  for (i in seq_len(n_loci)) {
    chr_rows <- gwas[.(regions$CHR[i]), nomatch = NULL]
    window_dt <- chr_rows[BP >= regions$START[i] & BP <= regions$END[i]]
    window_subsets[[i]] <- window_dt
    locus_subsets[[i]] <- window_dt[!is.na(BETA) & !is.na(SE) & SE > 0]
  }
  rm(gwas); gc(verbose = FALSE)

  locus_outcome <- function(status, n_snps = NA_integer_, result = NULL, converged = NA) {
    list(status = status, n_snps = n_snps, result = result, converged = converged)
  }

  process_one_locus <- function(i) {
    lead_rsid <- regions$LEAD_SNP[i]
    lead_chr <- regions$CHR[i]
    locus_name <- regions$LOCUS[i]

    tryCatch({
      locus_gwas <- locus_subsets[[i]]
      window_gwas <- window_subsets[[i]]

      if (nrow(locus_gwas) < 2L) {
        message(paste("Skipping locus", locus_name, "- fewer than 2 SNPs in window"))
        return(locus_outcome("skipped: fewer than 2 SNPs in region", nrow(locus_gwas)))
      }

      if (!has_rsid || all(is.na(locus_gwas$RSID))) {
        message(paste("Skipping locus", locus_name, "- no RSIDs available for LD computation"))
        return(locus_outcome("skipped: no RSIDs for LD computation", 0L))
      }
      locus_gwas <- locus_gwas[!is.na(RSID) & nchar(RSID) > 0L]

      ld_result <- compute_ld_matrix(locus_gwas$RSID, lead_chr, ancestry)
      if (is.null(ld_result)) {
        message(paste("Skipping locus", locus_name, "- LD matrix computation failed"))
        return(locus_outcome("skipped: LD matrix computation failed", 0L))
      }

      shared_rsids <- intersect(locus_gwas$RSID, ld_result$snps)
      if (length(shared_rsids) < 2L) {
        rm(ld_result); gc(verbose = FALSE)
        message(paste("Skipping locus", locus_name, "- too few shared SNPs between GWAS and LD panel"))
        return(locus_outcome("skipped: fewer than 2 SNPs in LD panel", length(shared_rsids)))
      }

      locus_gwas <- locus_gwas[match(shared_rsids, RSID)]
      ld_idx <- match(shared_rsids, ld_result$snps)
      R <- ld_result$matrix[ld_idx, ld_idx]
      rm(ld_result); gc(verbose = FALSE)

      z_scores <- locus_gwas$BETA / locus_gwas$SE
      n <- if (has_n && !is.na(as.numeric(locus_gwas$N[1L]))) as.numeric(locus_gwas$N[1L]) else default_n
      if (is.na(n)) {
        stop("Fine-mapping requires GWAS sample size")
      }

      locus_lbf <- run_susie_for_locus(
        z_scores = z_scores,
        ld_matrix = R,
        snp_info = locus_gwas,
        n = n,
        lead_snp = lead_rsid,
        max_causal = max_causal,
        coverage = coverage,
        min_abs_corr = min_abs_corr
      )
      rm(R, locus_gwas); gc(verbose = FALSE)

      if (nrow(locus_lbf) == 0L) {
        return(locus_outcome("skipped: SuSiE failed", length(shared_rsids)))
      }

      cs_summary <- attr(locus_lbf, "cs_summary")
      converged <- attr(locus_lbf, "converged")
      if (!is.null(cs_summary) && nrow(cs_summary) > 0L) {
        cs_summary <- cbind(data.table::data.table(LOCUS = locus_name), cs_summary)
      }

      lbf_nm <- grep("^LBF_[0-9]+$", names(locus_lbf), value = TRUE)
      finemap_cols <- intersect(c("RSID", "Z", "PIP", "CS", lbf_nm), names(locus_lbf))
      finemap_join <- locus_lbf[, ..finemap_cols]

      out_gwas <- merge(window_gwas, finemap_join, by = "RSID", all.x = TRUE, sort = FALSE)
      rm(window_gwas, finemap_join)

      out_file <- file.path(output_finemap_dir, paste0(locus_name, "_finemap.tsv.gz"))
      data.table::fwrite(out_gwas, out_file, sep = "\t", compress = "gzip")
      rm(out_gwas)

      list(status = "fine-mapped", n_snps = length(shared_rsids), result = locus_lbf,
           converged = if (is.null(converged)) NA else converged, cs_summary = cs_summary)
    }, error = function(e) {
      message(paste("Error processing locus", locus_name, ":", conditionMessage(e)))
      return(locus_outcome(paste("error:", conditionMessage(e))))
    })
  }

  finemap_workers <- min(
    calculate_parallelism(memory_per_worker_mb = 6000L),
    n_loci
  )
  locus_results <- parallel::mclapply(
    seq_len(n_loci),
    process_one_locus,
    mc.cores = finemap_workers
  )
  rm(locus_subsets, window_subsets); gc(verbose = FALSE)

  is_error <- vapply(locus_results, function(x) inherits(x, "try-error") || !is.list(x), FUN.VALUE = logical(1))
  if (any(is_error)) {
    message(paste(sum(is_error), "locus worker(s) returned errors"))
    locus_results[is_error] <- lapply(locus_results[is_error], function(x) locus_outcome("error: worker failed"))
  }

  regions <- add_region_status_columns(regions)
  regions$STATUS <- vapply(locus_results, function(x) x$status, FUN.VALUE = character(1))
  regions$N_SNPS <- vapply(locus_results, function(x) as.integer(x$n_snps), FUN.VALUE = integer(1))
  regions$N_CS <- vapply(locus_results, function(x) {
    if (is.null(x$cs_summary)) 0L else nrow(x$cs_summary)
  }, FUN.VALUE = integer(1))
  regions$CONVERGED <- vapply(locus_results, function(x) as.logical(x$converged)[1], FUN.VALUE = logical(1))
  data.table::fwrite(regions, regions_file, sep = "\t")

  cs_tables <- Filter(function(x) !is.null(x) && nrow(x) > 0L, lapply(locus_results, function(x) x$cs_summary))
  credible_sets <- if (length(cs_tables) == 0L) {
    empty_credible_set_summary()
  } else {
    data.table::rbindlist(cs_tables, fill = TRUE)
  }
  data.table::fwrite(credible_sets, credible_sets_file, sep = "\t")

  all_lbf <- Filter(Negate(is.null), lapply(locus_results, function(x) x$result))
  rm(locus_results)

  if (length(all_lbf) == 0L) {
    combined_lbf <- data.table::data.table(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), PIP = numeric(), CS = integer()
    )
  } else {
    combined_lbf <- data.table::rbindlist(all_lbf, fill = TRUE)
  }
  rm(all_lbf); gc(verbose = FALSE)

  if (nrow(combined_lbf) == 0L) {
    message("No loci produced SuSiE results.")
  }

  message(paste("Fine-mapping complete:",
    n_loci, "region(s),",
    sum(regions$STATUS == "fine-mapped"), "fine-mapped,",
    nrow(credible_sets), "credible set(s),",
    "outputs in", output_finemap_dir))

  write_finemap_complete_marker(completion_file, n_loci, output_finemap_dir)

  return(invisible(combined_lbf))
}


#' Read lead SNPs (SNP, CHR, BP and, when present, P) from a plink --clump output file.
#' @keywords internal
read_clumped_leads <- function(clumped_file) {
  header <- names(data.table::fread(clumped_file, nrows = 0))
  columns <- intersect(c("SNP", "CHR", "BP", "P"), header)
  leads <- data.table::fread(clumped_file, select = columns)
  if (!"P" %in% names(leads)) leads$P <- rep(NA_real_, nrow(leads))
  leads
}


#' Define fine-mapping regions from clumped lead SNPs.
#'
#' Each lead gets a window of \code{window_bp} either side. Windows that overlap
#' on the same chromosome are merged, so a signal is fine-mapped in one region
#' rather than once per nearby lead. A merged region never grows beyond
#' \code{max_region_bp}: if adding a lead would exceed it, the region is closed
#' and the next region starts immediately after it, so regions never overlap and
#' every lead lies inside exactly one region. Each region is named after its most
#' significant lead (\verb{<CHR>_<BP>}), or its first lead when P is unavailable.
#'
#' @param lead_snps data.frame with SNP, CHR, BP and optionally P
#' @param window_bp half-width of the window around each lead, in base pairs
#' @param max_region_bp maximum region width in base pairs (raised to 2 * window_bp if smaller)
#' @return data.table with LOCUS, CHR, START, END, LEAD_SNP, LEAD_BP, N_LEADS, LEADS
#' @keywords internal
define_finemap_regions <- function(lead_snps, window_bp, max_region_bp = 3e6) {
  CHR = BP = P = SNP = REGION = NULL # silence R CMD check notes for data.table columns
  empty <- data.table::data.table(
    LOCUS = character(), CHR = numeric(), START = numeric(), END = numeric(),
    LEAD_SNP = character(), LEAD_BP = numeric(), N_LEADS = integer(), LEADS = character()
  )
  if (is.null(lead_snps) || nrow(lead_snps) == 0L) return(empty)

  if (max_region_bp < 2 * window_bp) {
    message(paste0("max_region_kb (", max_region_bp / 1000, ") is smaller than a single window (",
                   2 * window_bp / 1000, " kb); using ", 2 * window_bp / 1000, " kb"))
    max_region_bp <- 2 * window_bp
  }

  leads <- data.table::data.table(
    SNP = as.character(lead_snps$SNP),
    CHR = as.numeric(lead_snps$CHR),
    BP = as.numeric(lead_snps$BP),
    P = if ("P" %in% names(lead_snps)) suppressWarnings(as.numeric(lead_snps$P)) else NA_real_
  )
  leads <- leads[!is.na(CHR) & !is.na(BP)]
  if (nrow(leads) == 0L) return(empty)
  data.table::setorder(leads, CHR, BP)

  region_id <- integer(nrow(leads))
  starts <- numeric(0)
  ends <- numeric(0)
  current_chr <- NA_real_
  for (i in seq_len(nrow(leads))) {
    chr <- leads$CHR[i]
    bp <- leads$BP[i]
    k <- length(starts)
    if (k == 0L || chr != current_chr || bp - window_bp > ends[k]) {
      starts <- c(starts, max(0, bp - window_bp))
      ends <- c(ends, bp + window_bp)
      current_chr <- chr
    } else {
      new_end <- max(ends[k], bp + window_bp)
      if (new_end - starts[k] <= max_region_bp) {
        ends[k] <- new_end
      } else if (bp <= ends[k]) {
        ends[k] <- max(ends[k], starts[k] + max_region_bp)
      } else {
        starts <- c(starts, ends[k] + 1)
        ends <- c(ends, bp + window_bp)
      }
    }
    region_id[i] <- length(starts)
  }
  leads$REGION <- region_id

  regions <- leads[, {
    lead_idx <- if (all(is.na(P))) 1L else which.min(P)
    list(
      CHR = CHR[1],
      LEAD_SNP = SNP[lead_idx],
      LEAD_BP = BP[lead_idx],
      N_LEADS = .N,
      LEADS = paste(SNP, collapse = ",")
    )
  }, by = REGION]
  regions$START <- starts[regions$REGION]
  regions$END <- ends[regions$REGION]
  regions$LOCUS <- paste0(
    format(regions$CHR, scientific = FALSE, trim = TRUE), "_",
    format(regions$LEAD_BP, scientific = FALSE, trim = TRUE, digits = 20)
  )
  regions[, c("LOCUS", "CHR", "START", "END", "LEAD_SNP", "LEAD_BP", "N_LEADS", "LEADS")]
}


#' Add the per-region status columns written to finemap_regions.tsv.
#' @keywords internal
add_region_status_columns <- function(regions) {
  regions <- data.table::copy(regions)
  regions$STATUS <- rep(NA_character_, nrow(regions))
  regions$N_SNPS <- rep(NA_integer_, nrow(regions))
  regions$N_CS <- rep(NA_integer_, nrow(regions))
  regions$CONVERGED <- rep(NA, nrow(regions))
  regions
}


#' Empty credible-set summary with the columns of finemap_credible_sets.tsv.
#' @keywords internal
empty_credible_set_summary <- function() {
  data.table::data.table(
    LOCUS = character(), CS = integer(), CS_SIZE = integer(), COVERAGE = numeric(),
    MIN_ABS_CORR = numeric(), MEAN_ABS_CORR = numeric(), MEDIAN_ABS_CORR = numeric(),
    LEAD_SNP = character(), LEAD_RSID = character(), LEAD_BP = numeric(), LEAD_PIP = numeric(),
    LOG_BF = numeric()
  )
}


#' Write a completion sentinel file (one line: expected lead count) for Snakemake.
#' @keywords internal
write_finemap_complete_marker <- function(completion_file, n_loci, output_finemap_dir = NULL) {
  if (is.null(completion_file)) {
    if (is.null(output_finemap_dir) || length(output_finemap_dir) != 1L || !nzchar(output_finemap_dir)) {
      stop("completion_file is NULL; provide output_finemap_dir or an explicit completion_file path")
    }
    completion_file <- finemap_completion_file(output_finemap_dir)
  }
  if (!is.character(completion_file) || length(completion_file) != 1L || !nzchar(completion_file)) {
    stop("completion_file must be a non-empty character path")
  }
  dir_path <- dirname(completion_file)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  writeLines(as.character(as.integer(n_loci)), completion_file)
}


#' Default finemap completion sentinel path from the output directory name.
#' @keywords internal
finemap_completion_file <- function(output_finemap_dir, completion_file = NULL) {
  if (!is.null(completion_file)) {
    return(completion_file)
  }
  label <- basename(normalizePath(output_finemap_dir, winslash = "/"))
  file.path(output_finemap_dir, paste0("finemap_complete_", label, ".txt"))
}


#' Compute LD correlation matrix from 1000 Genomes via plink
#' @param rsids character vector of RSIDs
#' @param chr chromosome number
#' @param ancestry ancestry code (EUR, EAS, etc.)
#' @return list with components \code{matrix} (numeric LD matrix) and
#'   \code{snps} (character vector of SNP IDs in matrix order), or NULL on failure
compute_ld_matrix <- function(rsids, chr, ancestry) {
  bfile <- file.path(thousand_genomes_dir, ancestry)
  tmpdir <- tempdir()
  snp_file <- tempfile(tmpdir = tmpdir, fileext = ".snps")
  out_prefix <- tempfile(tmpdir = tmpdir, pattern = "ld_")

  writeLines(rsids, snp_file)

  cmd <- paste(
    "plink1.9",
    "--bfile", bfile,
    "--chr", chr,
    "--extract", snp_file,
    "--r", "square",
    "--keep-allele-order",
    "--out", out_prefix
  )
  exit_code <- run_system(cmd, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(out_prefix, ".ld")
  log_file <- paste0(out_prefix, ".log")

  if (exit_code != 0 || !file.exists(ld_file)) {
    if (file.exists(log_file)) {
      plink_log <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
      warning(
        "PLINK LD matrix computation failed for chr ", chr, " (ancestry ", ancestry,
        ", exit code ", exit_code, "). Command:\n", cmd, "\nPLINK log:\n", plink_log,
        call. = FALSE, immediate. = TRUE
      )
    } else {
      warning(
        "PLINK LD matrix computation failed for chr ", chr, " (ancestry ", ancestry,
        ", exit code ", exit_code, "). Command:\n", cmd,
        call. = FALSE, immediate. = TRUE
      )
    }
    unlink(c(snp_file, paste0(out_prefix, c(".ld", ".nosex", ".log", ".bim", ".bed", ".fam"))),
           force = TRUE)
    return(NULL)
  }

  bim_file <- paste0(out_prefix, ".bim")
  if (file.exists(bim_file)) {
    snp_order <- data.table::fread(bim_file, header = FALSE, select = 2)$V2
  } else {
    bim_ref <- paste0(bfile, ".bim")
    all_bim <- data.table::fread(bim_ref, header = FALSE, select = c(1, 2, 4))
    colnames(all_bim) <- c("CHR", "SNP", "BP")
    all_bim <- all_bim[CHR == chr & SNP %in% rsids]
    data.table::setorder(all_bim, BP)
    snp_order <- all_bim$SNP
    rm(all_bim)
  }

  ld_raw <- data.table::fread(ld_file, header = FALSE)
  R <- as.matrix(ld_raw)
  rm(ld_raw)

  if (nrow(R) != length(snp_order) || ncol(R) != length(snp_order)) {
    rm(R)
    unlink(c(snp_file, paste0(out_prefix, c(".ld", ".nosex", ".log", ".bim", ".bed", ".fam"))),
           force = TRUE)
    return(NULL)
  }

  rownames(R) <- snp_order
  colnames(R) <- snp_order

  unlink(c(snp_file, paste0(out_prefix, c(".ld", ".nosex", ".log", ".bim", ".bed", ".fam"))),
         force = TRUE)

  return(list(matrix = R, snps = snp_order))
}


#' Wrapper so tests can mock SuSiE with \code{testthat::local_mocked_bindings()}.
#' @keywords internal
run_susie_rss_impl <- function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
  susieR::susie_rss(
    z = z, R = R, n = n, L = L, coverage = coverage, min_abs_corr = min_abs_corr, verbose = verbose
  )
}


#' Build per-credible-set LBF columns (\code{LBF_1}, \code{LBF_2}, ...) from a SuSiE fit.
#' @keywords internal
susie_lbf_columns <- function(fitted, p) {
  lv <- fitted$lbf_variable
  if (is.null(lv) || !is.matrix(lv) || nrow(lv) < 1L || ncol(lv) != p) {
    return(data.table::data.table())
  }

  n_row <- nrow(lv)
  cs_list <- fitted$sets$cs
  cs_idx <- fitted$sets$cs_index
  col_vecs <- vector("list", 0L)

  if (!is.null(cs_list) && length(cs_list) > 0L) {
    n_cs <- length(cs_list)
    use_idx <- !is.null(cs_idx) && length(cs_idx) == n_cs
    for (j in seq_len(n_cs)) {
      row_i <- if (use_idx) cs_idx[j] else j
      if (is.finite(row_i) && row_i >= 1L && row_i <= n_row) {
        col_vecs[[paste0("LBF_", j)]] <- lv[row_i, ]
      } else {
        col_vecs[[paste0("LBF_", j)]] <- rep(NA_real_, p)
      }
    }
  } else {
    for (j in seq_len(n_row)) {
      col_vecs[[paste0("LBF_", j)]] <- lv[j, ]
    }
  }

  data.table::as.data.table(col_vecs)
}


#' Run SuSiE on a single locus
#' @param z_scores numeric vector of z-scores
#' @param ld_matrix square LD correlation matrix
#' @param snp_info data frame with SNP, CHR, BP, RSID columns (same order as z_scores)
#' @param n sample size
#' @param lead_snp RSID of the lead (clumped) SNP
#' @param max_causal max causal signals (L)
#' @param coverage credible set coverage
#' @param min_abs_corr minimum absolute correlation for CS purity
#' @return data.table with SNP, CHR, BP, RSID, Z, PIP, CS, and \code{LBF_1}, \code{LBF_2}, ... per signal.
#'   Attributes \code{cs_summary} (one row per credible set, see \code{susie_credible_set_summary()})
#'   and \code{converged} (SuSiE convergence flag) are attached.
run_susie_for_locus <- function(z_scores, ld_matrix, snp_info, n, lead_snp,
                                max_causal = 10, coverage = 0.95,
                                min_abs_corr = 0.5) {

  p <- length(z_scores)
  fitted <- tryCatch(
    run_susie_rss_impl(
      z = z_scores,
      R = ld_matrix,
      n = n,
      L = max_causal,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      verbose = FALSE
    ),
    error = function(e) {
      message(paste("SuSiE failed for locus", lead_snp, ":", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(fitted)) {
    return(data.table::data.table(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), PIP = numeric(), CS = integer()
    ))
  }

  lbf_wide <- susie_lbf_columns(fitted, p)
  pip <- if (!is.null(fitted$pip) && length(fitted$pip) == p) as.numeric(fitted$pip) else rep(NA_real_, p)
  cs_summary <- susie_credible_set_summary(fitted, snp_info, pip)
  converged <- if (is.null(fitted$converged)) NA else isTRUE(fitted$converged)

  cs_membership <- rep(NA_integer_, p)
  if (!is.null(fitted$sets) && !is.null(fitted$sets$cs)) {
    for (cs_i in seq_along(fitted$sets$cs)) {
      snp_indices <- fitted$sets$cs[[cs_i]]
      cs_membership[snp_indices] <- cs_i
    }
  }
  rm(fitted); gc(verbose = FALSE)

  result <- data.table::data.table(
    SNP = snp_info$SNP,
    CHR = snp_info$CHR,
    BP = snp_info$BP,
    RSID = snp_info$RSID,
    Z = z_scores,
    PIP = pip,
    CS = cs_membership
  )

  if (ncol(lbf_wide) > 0L) {
    result <- cbind(result, lbf_wide)
  }
  data.table::setattr(result, "cs_summary", cs_summary)
  data.table::setattr(result, "converged", converged)
  result
}


#' Summarise SuSiE credible sets: one row per credible set with its size, coverage,
#' purity (absolute LD correlations within the set), lead variant (highest PIP) and log Bayes factor.
#' @param fitted SuSiE fit (as returned by \code{susieR::susie_rss})
#' @param snp_info data frame with SNP, BP, RSID columns in the same order as the fitted variants
#' @param pip numeric vector of posterior inclusion probabilities
#' @return data.table with CS, CS_SIZE, COVERAGE, MIN_ABS_CORR, MEAN_ABS_CORR, MEDIAN_ABS_CORR,
#'   LEAD_SNP, LEAD_RSID, LEAD_BP, LEAD_PIP, LOG_BF
#' @keywords internal
susie_credible_set_summary <- function(fitted, snp_info, pip) {
  LOCUS = NULL # silence R CMD check note for data.table column
  cs_list <- fitted$sets$cs
  summary <- empty_credible_set_summary()
  summary[, LOCUS := NULL]
  if (is.null(cs_list) || length(cs_list) == 0L) return(summary)

  sets <- fitted$sets
  purity <- sets$purity
  purity_value <- function(j, column) {
    if (is.null(purity) || !column %in% colnames(purity) || nrow(purity) < j) return(NA_real_)
    as.numeric(purity[j, column])
  }
  cs_index <- sets$cs_index
  lbf <- fitted$lbf

  rows <- lapply(seq_along(cs_list), function(j) {
    idx <- as.integer(cs_list[[j]])
    lead <- if (all(is.na(pip[idx]))) idx[1] else idx[which.max(pip[idx])]
    lbf_row <- if (!is.null(cs_index) && length(cs_index) >= j) cs_index[j] else j
    data.table::data.table(
      CS = j,
      CS_SIZE = length(idx),
      COVERAGE = if (!is.null(sets$coverage) && length(sets$coverage) >= j) as.numeric(sets$coverage[j]) else NA_real_,
      MIN_ABS_CORR = purity_value(j, "min.abs.corr"),
      MEAN_ABS_CORR = purity_value(j, "mean.abs.corr"),
      MEDIAN_ABS_CORR = purity_value(j, "median.abs.corr"),
      LEAD_SNP = as.character(snp_info$SNP[lead]),
      LEAD_RSID = as.character(snp_info$RSID[lead]),
      LEAD_BP = as.numeric(snp_info$BP[lead]),
      LEAD_PIP = pip[lead],
      LOG_BF = if (!is.null(lbf) && length(lbf) >= lbf_row) as.numeric(lbf[lbf_row]) else NA_real_
    )
  })
  data.table::rbindlist(rows)
}