#' Run SuSiE fine-mapping across all clumped loci in a GWAS
#'
#' Lead SNPs from plink --clump output are turned into fine-mapping regions by
#' \code{define_finemap_regions()}: a window of \code{window_kb} either side of
#' each lead, with overlapping windows merged into a single region (up to
#' \code{max_region_kb} wide) so that a signal is only fine-mapped once. For each
#' region an LD matrix is taken from the LD reference (\code{ld_reference}): the
#' precomputed Pan-UKB LD pack, matched on SNP keys, or the 1000 Genomes panel via
#' plink, matched on RSIDs; then \code{susieR::susie_rss} is run.
#'
#' Outputs written to \code{output_finemap_dir}:
#' \itemize{
#'   \item \verb{<CHR>_<BP>_finemap.tsv.gz} per region (BP is the most significant
#'     lead SNP in the region): all GWAS variants in the region with SuSiE
#'     Z-scores, posterior inclusion probabilities (\code{PIP}), credible-set
#'     membership (\code{CS}) and per-credible-set log Bayes factors
#'     \code{LBF_1}, \code{LBF_2}, ...
#'   \item \code{finemap_regions.tsv}: one row per region (bounds, leads, status,
#'     number of variants with LD, number of credible sets, SuSiE convergence, LD
#'     reference, its sample size \code{N_REF}, the share of variant pairs beyond the
#'     LD pack's stored radius \code{FRAC_BEYOND_LD_RADIUS} (set to r = 0) and the
#'     smallest eigenvalue of the pack LD matrix, \code{LD_MIN_EIGENVALUE}; when it is below 0
#'     the matrix was made positive semi-definite by clipping its negative eigenvalues).
#'   \item \code{finemap_credible_sets.tsv}: one row per credible set (size,
#'     coverage, purity, lead variant and its PIP, log Bayes factor).
#' }
#'
#' @param gwas filename or dataframe of a standardised GWAS
#' @param clumped_file plink --clump output file
#' @param ancestry ancestry code: EUR, EAS, AFR, AMR or SAS for 1000 Genomes; EUR, AFR, AMR, CSA (or SAS),
#'   EAS or MID for the Pan-UKB LD pack
#' @param default_n GWAS sample size when not inferrable from \code{gwas}; see
#' @param output_finemap_dir directory to write one
#'   \verb{<CHR>_<BP>_finemap.tsv.gz} per fine-mapping region.
#' @param completion_file path to a sentinel file written on successful completion
#'   (one line: number of fine-mapping regions). If NULL, defaults to
#'   \verb{<output_finemap_dir>/finemap_complete_<dir_basename>.txt}.
#' @param window_kb half-width of the fine-mapping window around each lead SNP in kb
#' @param max_causal SuSiE L: "auto" (default) for L = 1 when the GWAS sample size is more than \code{max_n_ratio}
#'   times the LD reference's and 10 otherwise (see \code{choose_max_causal()}), or a number
#' @param max_n_ratio for \code{max_causal = "auto"}: largest GWAS N / LD reference N that still uses L = 10
#' @param coverage credible set coverage (default 0.95)
#' @param min_abs_corr minimum absolute correlation for credible sets (default 0.5)
#' @param max_region_kb maximum width of a merged fine-mapping region, or the tile width, in kb (default 3000)
#' @param region_strategy "merge" (default), "none" or "tile"; see \code{define_finemap_regions()}
#' @param ld_bfile PLINK prefix of in-sample (or other user-supplied) genotypes in \code{build}; when given, LD
#'   comes from these instead of \code{ld_reference}
#' @param remove_ld_outliers if TRUE, variants flagged by the LD-mismatch check (\code{ld_mismatch_qc()}) are
#'   removed and SuSiE is run again; by default they are only reported
#' @param ld_reference LD source, one of \code{ld_reference_options} ("1000genomes" or "pan_ukb")
#' @param build reference build of the GWAS positions (used to match SNP keys in the Pan-UKB LD pack)
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
                         max_region_kb = 3000,
                         ld_reference = ld_reference_options$thousand_genomes,
                         build = reference_builds$GRCh38,
                         region_strategy = "merge",
                         max_n_ratio = 10,
                         ld_bfile = NULL,
                         remove_ld_outliers = FALSE) {
  validate_ld_reference(ld_reference)
  user_ld <- has_user_ld(ld_bfile)
  id_column <- ld_reference_id_column(ld_reference, ld_bfile)
  use_pack <- !user_ld && ld_reference == ld_reference_options$pan_ukb
  ld_source <- if (user_ld) paste0("user:", basename(ld_bfile)) else ld_reference
  n_ref <- ld_reference_n_samples(ld_reference, ancestry, build, ld_bfile)

  gwas <- get_file_or_dataframe(gwas)
  data.table::setDT(gwas)

  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF", "N", "N_EFF", "N_CASES", "N_CONTROLS")
  for (col in intersect(numeric_cols, colnames(gwas))) {
    data.table::set(gwas, j = col, value = as.numeric(gwas[[col]]))
  }

  data.table::setkey(gwas, CHR, BP)

  gwas_n <- analysis_sample_size(gwas, default_n)
  if (length(gwas_n) != 1L || is.na(gwas_n)) {
    stop("Fine-mapping requires GWAS sample size")
  }
  l_choice <- choose_max_causal(max_causal, gwas_n, n_ref, max_n_ratio)
  max_causal <- l_choice$L
  message("SuSiE L = ", max_causal, " (", l_choice$rule, "; GWAS N = ", round(gwas_n), ", LD reference N = ", n_ref, ")")

  if (!dir.exists(output_finemap_dir)) {
    dir.create(output_finemap_dir, recursive = TRUE)
  }

  if (is.null(completion_file)) {
    completion_file <- finemap_completion_file(output_finemap_dir)
  }

  lead_snps <- read_clumped_leads(clumped_file)
  regions <- define_finemap_regions(lead_snps, window_bp = window_kb * 1000, max_region_bp = max_region_kb * 1000,
                                    strategy = region_strategy, build = build)
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
    list(status = status, n_snps = n_snps, result = result, converged = converged,
         beyond_radius = NA_real_, min_eigenvalue = NA_real_, ld_s = NA_real_, n_outliers = NA_integer_,
         n_allele_mismatch = NA_integer_, n = NA_real_)
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

      if (id_column == "RSID" && (!has_rsid || all(is.na(locus_gwas$RSID)))) {
        message(paste("Skipping locus", locus_name, "- no RSIDs available for LD computation"))
        return(locus_outcome("skipped: no RSIDs for LD computation", 0L))
      }
      ids <- locus_gwas[[id_column]]
      locus_gwas <- locus_gwas[!is.na(ids) & nchar(ids) > 0L & !duplicated(ids)]

      ld_result <- get_ld_matrix(locus_gwas[[id_column]], lead_chr, ancestry, build, ld_reference, ld_bfile)
      if (is.null(ld_result)) {
        message(paste("Skipping locus", locus_name, "- LD matrix computation failed"))
        return(locus_outcome("skipped: LD matrix computation failed", 0L))
      }

      shared_ids <- intersect(locus_gwas[[id_column]], ld_result$snps)
      if (length(shared_ids) < 2L) {
        rm(ld_result); gc(verbose = FALSE)
        message(paste("Skipping locus", locus_name, "- too few shared SNPs between GWAS and LD panel"))
        return(locus_outcome("skipped: fewer than 2 SNPs in LD panel", length(shared_ids)))
      }

      locus_gwas <- locus_gwas[match(shared_ids, locus_gwas[[id_column]])]
      ld_idx <- match(shared_ids, ld_result$snps)
      R <- ld_result$matrix[ld_idx, ld_idx, drop = FALSE]
      n_allele_mismatch <- 0L
      if (!is.null(ld_result$a1)) {
        # LD matched on RSID: orient r to the GWAS effect allele, drop variants whose alleles differ
        aligned <- align_ld_to_gwas(R, locus_gwas, ld_result$a1[ld_idx], ld_result$a2[ld_idx])
        R <- aligned$matrix
        locus_gwas <- locus_gwas[aligned$keep]
        shared_ids <- shared_ids[aligned$keep]
        n_allele_mismatch <- aligned$n_mismatch
        if (length(shared_ids) < 2L) {
          return(locus_outcome("skipped: fewer than 2 SNPs with alleles matching the LD panel", length(shared_ids)))
        }
      }
      beyond_radius <- if (is.null(ld_result$beyond_band_fraction)) NA_real_ else ld_result$beyond_band_fraction
      rm(ld_result); gc(verbose = FALSE)
      min_eigenvalue <- NA_real_
      if (use_pack) {
        # pairs beyond the pack's radius are r = 0 and r is quantised; keep R positive semi-definite
        psd <- make_ld_positive_semidefinite(R)
        R <- psd$matrix
        min_eigenvalue <- psd$min_eigenvalue
      }

      z_scores <- locus_gwas$BETA / locus_gwas$SE
      # median per-variant N_EFF (binary traits) or N in the region; the input N where the file has none
      n <- analysis_sample_size(locus_gwas, default_n)
      if (is.na(n)) {
        stop("Fine-mapping requires GWAS sample size")
      }

      mismatch <- ld_mismatch_qc(z_scores, R, n)
      ld_outlier <- mismatch$outlier
      if (isTRUE(remove_ld_outliers) && any(ld_outlier)) {
        message(paste("Locus", locus_name, ": removing", sum(ld_outlier), "variant(s) flagged by the LD-mismatch check"))
        keep <- !ld_outlier
        R <- R[keep, keep, drop = FALSE]
        z_scores <- z_scores[keep]
        locus_gwas <- locus_gwas[keep]
        shared_ids <- shared_ids[keep]
        ld_outlier <- ld_outlier[keep]
        if (length(shared_ids) < 2L) {
          return(locus_outcome("skipped: fewer than 2 SNPs after removing LD outliers", length(shared_ids)))
        }
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
        return(locus_outcome("skipped: SuSiE failed", length(shared_ids)))
      }
      locus_lbf$LD_OUTLIER <- ld_outlier
      locus_lbf$LOCUS <- locus_name
      in_core <- tile_core_flags(locus_lbf$BP, regions, i)
      locus_lbf$TILE_CORE <- in_core

      cs_summary <- attr(locus_lbf, "cs_summary")
      converged <- attr(locus_lbf, "converged")
      if (!is.null(cs_summary) && nrow(cs_summary) > 0L) {
        cs_summary <- cbind(data.table::data.table(LOCUS = locus_name), cs_summary)
        # tiles overlap: a credible set is reported by the tile whose centre is nearest its lead variant
        cs_summary$LEAD_IN_TILE_CORE <- tile_core_flags(cs_summary$LEAD_BP, regions, i)
        cs_summary$LD_MISMATCH_S <- mismatch$s
        cs_summary$N_LD_OUTLIERS_IN_CS <- vapply(cs_summary$CS, function(k) {
          as.integer(sum(locus_lbf$LD_OUTLIER[!is.na(locus_lbf$CS) & locus_lbf$CS == k]))
        }, integer(1))
        cs_summary$LD_FLAG <- ld_mismatch_flag(mismatch$s, sum(mismatch$outlier), cs_summary$N_LD_OUTLIERS_IN_CS,
                                               regions$MHC[i])
      }

      lbf_nm <- grep("^LBF_[0-9]+$", names(locus_lbf), value = TRUE)
      finemap_cols <- intersect(c(id_column, "Z", "PIP", "CS", "LD_OUTLIER", "TILE_CORE", lbf_nm), names(locus_lbf))
      finemap_join <- locus_lbf[, ..finemap_cols]

      out_gwas <- merge(window_gwas, finemap_join, by = id_column, all.x = TRUE, sort = FALSE)
      rm(window_gwas, finemap_join)

      out_file <- file.path(output_finemap_dir, paste0(locus_name, "_finemap.tsv.gz"))
      data.table::fwrite(out_gwas, out_file, sep = "\t", compress = "gzip")
      rm(out_gwas)

      list(status = "fine-mapped", n_snps = length(shared_ids), result = locus_lbf,
           converged = if (is.null(converged)) NA else converged, cs_summary = cs_summary,
           beyond_radius = beyond_radius, min_eigenvalue = min_eigenvalue,
           ld_s = mismatch$s, n_outliers = sum(mismatch$outlier), n_allele_mismatch = n_allele_mismatch, n = n)
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
  regions$LD_REFERENCE <- ld_source
  regions$N_REF <- n_ref
  regions$N_GWAS <- vapply(locus_results, function(x) as.numeric(x$n %||% NA)[1], FUN.VALUE = numeric(1))
  regions$L <- max_causal
  regions$L_RULE <- l_choice$rule
  regions$LD_MISMATCH_S <- vapply(locus_results, function(x) as.numeric(x$ld_s %||% NA)[1], FUN.VALUE = numeric(1))
  regions$N_LD_OUTLIERS <- vapply(locus_results, function(x) as.integer(x$n_outliers %||% NA)[1], FUN.VALUE = integer(1))
  regions$N_ALLELE_MISMATCH <- vapply(locus_results, function(x) as.integer(x$n_allele_mismatch %||% NA)[1],
                                      FUN.VALUE = integer(1))
  regions$LD_OUTLIERS_REMOVED <- isTRUE(remove_ld_outliers)
  regions$FRAC_BEYOND_LD_RADIUS <- vapply(locus_results, function(x) as.numeric(x$beyond_radius %||% NA)[1],
                                          FUN.VALUE = numeric(1))
  regions$LD_MIN_EIGENVALUE <- vapply(locus_results, function(x) as.numeric(x$min_eigenvalue %||% NA)[1],
                                      FUN.VALUE = numeric(1))
  regions$LD_FLAG <- ifelse(regions$STATUS == "fine-mapped",
                            ld_mismatch_flag(regions$LD_MISMATCH_S, regions$N_LD_OUTLIERS, 0L, regions$MHC),
                            NA_character_)
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
  leads <- leads[!is.na(leads$SNP) & nzchar(as.character(leads$SNP))]
  if (!"P" %in% names(leads)) leads$P <- rep(NA_real_, nrow(leads))
  leads
}


#' Define fine-mapping regions from clumped lead SNPs.
#'
#' Each lead gets a window of \code{window_bp} either side.  Three strategies:
#' \describe{
#'   \item{merge (default)}{Windows that overlap on the same chromosome are merged, so a signal is fine-mapped
#'     in one region rather than once per nearby lead.  A merged region never grows beyond
#'     \code{max_region_bp}: if adding a lead would exceed it, the region is closed and the next region starts
#'     immediately after it, so regions never overlap and every lead lies inside exactly one region.}
#'   \item{none}{One window per lead, not merged.  Windows of nearby leads overlap, so one signal can be
#'     fine-mapped (and reported) more than once.}
#'   \item{tile}{Overlapping windows are first grouped as for \code{merge} (with no size limit); each group is
#'     then covered by tiles \code{max_region_bp} wide, starting every \code{max_region_bp / 3} (3 Mb tiles
#'     every 1 Mb by default, as in Weissbrod et al. 2020).  Tiles without a lead are dropped.  Each variant's
#'     result is taken from the tile whose centre is nearest (\code{TILE_CORE} in the outputs).}
#' }
#' Each region is named after its most significant lead (\verb{<CHR>_<BP>}), or its first lead when P is
#' unavailable; tiles sharing a lead get a suffix.  \code{MHC} marks regions overlapping the MHC
#' (chr6:28.5-33.5 Mb), where LD is long-range and complex and fine-mapping is unreliable.
#'
#' @param lead_snps data.frame with SNP, CHR, BP and optionally P
#' @param window_bp half-width of the window around each lead, in base pairs
#' @param max_region_bp maximum region width (merge) or tile width (tile) in base pairs (raised to
#'   2 * window_bp if smaller)
#' @param strategy "merge", "none" or "tile"
#' @param build reference build of the positions, for the MHC flag
#' @return data.table with LOCUS, CHR, START, END, LEAD_SNP, LEAD_BP, N_LEADS, LEADS, REGION_STRATEGY,
#'   CLUSTER (tile only) and MHC
#' @keywords internal
define_finemap_regions <- function(lead_snps, window_bp, max_region_bp = 3e6, strategy = "merge",
                                   build = reference_builds$GRCh38) {
  CHR = BP = P = SNP = REGION = NULL # silence R CMD check notes for data.table columns
  strategy <- match.arg(strategy, c("merge", "none", "tile"))
  empty <- data.table::data.table(
    LOCUS = character(), CHR = numeric(), START = numeric(), END = numeric(),
    LEAD_SNP = character(), LEAD_BP = numeric(), N_LEADS = integer(), LEADS = character(),
    REGION_STRATEGY = character(), CLUSTER = integer(), MHC = logical()
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

  # each region: its bounds and the leads (row numbers of `leads`) inside it
  region_list <- list()
  cluster_of <- integer(0)
  add_region <- function(start, end, rows, cluster = NA_integer_) {
    region_list[[length(region_list) + 1L]] <<- list(start = start, end = end, rows = rows)
    cluster_of <<- c(cluster_of, cluster)
  }

  if (strategy == "none") {
    for (i in seq_len(nrow(leads))) add_region(max(0, leads$BP[i] - window_bp), leads$BP[i] + window_bp, i)
  } else {
    cap <- if (strategy == "merge") max_region_bp else Inf
    starts <- numeric(0); ends <- numeric(0); members <- list(); current_chr <- NA_real_
    for (i in seq_len(nrow(leads))) {
      chr <- leads$CHR[i]; bp <- leads$BP[i]; k <- length(starts)
      if (k == 0L || chr != current_chr || bp - window_bp > ends[k]) {
        starts <- c(starts, max(0, bp - window_bp)); ends <- c(ends, bp + window_bp)
        members[[k + 1L]] <- i; current_chr <- chr
      } else {
        new_end <- max(ends[k], bp + window_bp)
        if (new_end - starts[k] <= cap) {
          ends[k] <- new_end; members[[k]] <- c(members[[k]], i)
        } else if (bp <= ends[k]) {
          ends[k] <- max(ends[k], starts[k] + cap); members[[k]] <- c(members[[k]], i)
        } else {
          starts <- c(starts, ends[k] + 1); ends <- c(ends, bp + window_bp); members[[k + 1L]] <- i
        }
      }
    }
    for (k in seq_along(starts)) {
      if (strategy == "merge" || ends[k] - starts[k] <= max_region_bp) {
        add_region(starts[k], ends[k], members[[k]], if (strategy == "tile") k else NA_integer_)
        next
      }
      step <- max_region_bp / 3
      tile_starts <- seq(starts[k], max(starts[k], ends[k] - max_region_bp + step), by = step)
      for (tile_start in tile_starts) {
        tile_end <- min(tile_start + max_region_bp, ends[k])
        rows <- members[[k]][leads$BP[members[[k]]] >= tile_start & leads$BP[members[[k]]] <= tile_end]
        if (length(rows) > 0) add_region(tile_start, tile_end, rows, k)
      }
    }
  }

  regions <- data.table::rbindlist(lapply(region_list, function(r) {
    region_leads <- leads[r$rows]
    lead_idx <- if (all(is.na(region_leads$P))) 1L else which.min(region_leads$P)
    data.table::data.table(CHR = region_leads$CHR[1], START = r$start, END = r$end,
                           LEAD_SNP = region_leads$SNP[lead_idx], LEAD_BP = region_leads$BP[lead_idx],
                           N_LEADS = nrow(region_leads), LEADS = paste(region_leads$SNP, collapse = ","))
  }))
  regions$REGION_STRATEGY <- strategy
  regions$CLUSTER <- cluster_of
  locus <- paste0(format(regions$CHR, scientific = FALSE, trim = TRUE), "_",
                  format(regions$LEAD_BP, scientific = FALSE, trim = TRUE, digits = 20))
  duplicated_name <- locus %in% locus[duplicated(locus)]
  locus[duplicated_name] <- paste0(locus[duplicated_name], "_",
                                   format(regions$START[duplicated_name], scientific = FALSE, trim = TRUE))
  regions$LOCUS <- locus
  mhc <- mhc_bounds(build)
  regions$MHC <- regions$CHR == 6 & regions$START <= mhc[2] & regions$END >= mhc[1]
  if (any(regions$MHC)) {
    message(sum(regions$MHC), " fine-mapping region(s) overlap the MHC (chr6:", mhc[1], "-", mhc[2],
            "); their results are flagged (MHC column) and should be interpreted with caution")
  }
  regions[, c("LOCUS", "CHR", "START", "END", "LEAD_SNP", "LEAD_BP", "N_LEADS", "LEADS", "REGION_STRATEGY",
              "CLUSTER", "MHC")]
}


#' Bounds of the MHC (the GRC's chr6 MHC region) in a reference build.
#' @keywords internal
mhc_bounds <- function(build = reference_builds$GRCh38) {
  if (identical(build, reference_builds$GRCh37)) c(28477797, 33448354) else c(28510120, 33480577)
}


#' For tiled regions: which positions take their result from tile \code{i} (the tile whose centre is nearest).
#' Always TRUE for the merge and none strategies.
#' @keywords internal
tile_core_flags <- function(bp, regions, i) {
  if (length(bp) == 0 || !"CLUSTER" %in% names(regions) || is.na(regions$CLUSTER[i])) return(rep(TRUE, length(bp)))
  same <- which(regions$CLUSTER == regions$CLUSTER[i])
  if (length(same) <= 1) return(rep(TRUE, length(bp)))
  centres <- (regions$START[same] + regions$END[same]) / 2
  nearest <- same[apply(abs(outer(as.numeric(bp), centres, "-")), 1, which.min)]
  nearest == i
}


#' Add the per-region status columns written to finemap_regions.tsv.
#' @keywords internal
add_region_status_columns <- function(regions) {
  regions <- data.table::copy(regions)
  regions$STATUS <- rep(NA_character_, nrow(regions))
  regions$N_SNPS <- rep(NA_integer_, nrow(regions))
  regions$N_CS <- rep(NA_integer_, nrow(regions))
  regions$CONVERGED <- rep(NA, nrow(regions))
  regions$LD_REFERENCE <- rep(NA_character_, nrow(regions))
  regions$N_REF <- rep(NA_real_, nrow(regions))
  regions$FRAC_BEYOND_LD_RADIUS <- rep(NA_real_, nrow(regions))
  regions$LD_MIN_EIGENVALUE <- rep(NA_real_, nrow(regions))
  regions$N_GWAS <- rep(NA_real_, nrow(regions))
  regions$L <- rep(NA_integer_, nrow(regions))
  regions$L_RULE <- rep(NA_character_, nrow(regions))
  regions$LD_MISMATCH_S <- rep(NA_real_, nrow(regions))
  regions$N_LD_OUTLIERS <- rep(NA_integer_, nrow(regions))
  regions$N_ALLELE_MISMATCH <- rep(NA_integer_, nrow(regions))
  regions$LD_OUTLIERS_REMOVED <- rep(NA, nrow(regions))
  regions$LD_FLAG <- rep(NA_character_, nrow(regions))
  regions
}


#' Threshold on the SuSiE-RSS estimate s above which a region's z-scores and LD are reported as inconsistent.
#' s is 0 when they are consistent; values well above 0 suggest LD mismatch (Zou et al. 2022).
#' @keywords internal
ld_mismatch_s_threshold <- 0.1


#' Flag fine-mapping results that may be affected by LD mismatch or the MHC.
#'
#' \describe{
#'   \item{OUTLIER_IN_CS}{a variant in the credible set is an LD outlier (its z-score is unlikely given the
#'     others' and the LD; a flipped allele is the typical cause). The credible set itself is suspect.}
#'   \item{LD_MISMATCH}{the region's estimated inconsistency s is above \code{ld_mismatch_s_threshold}: the LD
#'     reference fits the z-scores poorly, so credible sets other than the lead one may be spurious.}
#'   \item{OUTLIERS_IN_REGION}{the region has LD outliers, outside this credible set.}
#'   \item{MHC}{the region overlaps the MHC.}
#' }
#' Several flags are joined with ";"; "OK" when none applies.  The flags matter most for multi-signal
#' fine-mapping (L > 1): with L = 1 the posterior inclusion probabilities don't depend on the LD matrix.
#' @param s estimated inconsistency per region (NA if not computed)
#' @param n_outliers_region LD outliers in the region
#' @param n_outliers_cs LD outliers in the credible set (0 for a region-level flag)
#' @param mhc whether the region overlaps the MHC
#' @return character vector of flags
#' @keywords internal
ld_mismatch_flag <- function(s, n_outliers_region, n_outliers_cs = 0L, mhc = FALSE) {
  if (length(s) == 0 || length(n_outliers_region) == 0) return(character())
  n <- max(length(s), length(n_outliers_region), length(n_outliers_cs), length(mhc))
  s <- rep_len(as.numeric(s), n)
  n_outliers_region <- rep_len(as.integer(n_outliers_region), n)
  n_outliers_cs <- rep_len(as.integer(n_outliers_cs), n)
  mhc <- rep_len(as.logical(mhc), n)
  vapply(seq_len(n), function(i) {
    in_cs <- !is.na(n_outliers_cs[i]) && n_outliers_cs[i] > 0
    flags <- c(
      if (in_cs) "OUTLIER_IN_CS",
      if (!is.na(s[i]) && s[i] > ld_mismatch_s_threshold) "LD_MISMATCH",
      if (!in_cs && !is.na(n_outliers_region[i]) && n_outliers_region[i] > 0) "OUTLIERS_IN_REGION",
      if (isTRUE(mhc[i])) "MHC"
    )
    if (length(flags) == 0) "OK" else paste(flags, collapse = ";")
  }, character(1))
}


#' Fine-mapping results and LD diagnostics of one output directory, for the HTML reports.
#'
#' @param finemap_dir directory written by \code{finemap_gwas()}
#' @return list with \code{regions} and \code{credible_sets} (data.tables, empty when missing) and \code{summary},
#'   a one-row data.table: regions, regions fine-mapped, L and the rule that set it, LD reference, credible sets,
#'   flagged credible sets, regions with s above the threshold, regions with LD outliers, MHC regions
#' @export
finemap_diagnostics <- function(finemap_dir) {
  read_or_empty <- function(f) {
    if (file.exists(f) && file.size(f) > 0) data.table::fread(f) else data.table::data.table()
  }
  regions <- read_or_empty(file.path(finemap_dir, "finemap_regions.tsv"))
  credible_sets <- read_or_empty(file.path(finemap_dir, "finemap_credible_sets.tsv"))
  column <- function(dt, name, default = NA) if (name %in% names(dt)) dt[[name]] else rep(default, nrow(dt))
  count <- function(x) as.integer(sum(x, na.rm = TRUE))
  cs_flags <- column(credible_sets, "LD_FLAG")
  summary <- data.table::data.table(
    REGIONS = nrow(regions),
    FINE_MAPPED = count(column(regions, "STATUS") == "fine-mapped"),
    L = paste(unique(stats::na.omit(column(regions, "L"))), collapse = ", "),
    L_RULE = paste(unique(stats::na.omit(column(regions, "L_RULE"))), collapse = "; "),
    LD_REFERENCE = paste(unique(stats::na.omit(column(regions, "LD_REFERENCE"))), collapse = ", "),
    CREDIBLE_SETS = nrow(credible_sets),
    CREDIBLE_SETS_FLAGGED = count(!is.na(cs_flags) & cs_flags != "OK"),
    REGIONS_LD_MISMATCH = count(column(regions, "LD_MISMATCH_S") > ld_mismatch_s_threshold),
    REGIONS_WITH_LD_OUTLIERS = count(column(regions, "N_LD_OUTLIERS") > 0),
    MHC_REGIONS = count(column(regions, "MHC", FALSE))
  )
  list(regions = regions, credible_sets = credible_sets, summary = summary)
}



#' Empty credible-set summary with the columns of finemap_credible_sets.tsv.
#' @keywords internal
empty_credible_set_summary <- function() {
  data.table::data.table(
    LOCUS = character(), CS = integer(), CS_SIZE = integer(), COVERAGE = numeric(),
    MIN_ABS_CORR = numeric(), MEAN_ABS_CORR = numeric(), MEDIAN_ABS_CORR = numeric(),
    LEAD_SNP = character(), LEAD_RSID = character(), LEAD_BP = numeric(), LEAD_PIP = numeric(),
    LOG_BF = numeric(), LEAD_IN_TILE_CORE = logical(), LD_MISMATCH_S = numeric(), N_LD_OUTLIERS_IN_CS = integer(),
    LD_FLAG = character()
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
#' @param build reference build of the GWAS; selects the 1000 Genomes panel (GRCh38 or GRCh37)
#' @return list with components \code{matrix} (numeric LD matrix),
#'   \code{snps} (character vector of SNP IDs in matrix order) and the panel's alleles \code{a1} (the allele r
#'   refers to) and \code{a2}, or NULL on failure
compute_ld_matrix <- function(rsids, chr, ancestry, build = reference_builds$GRCh38) {
  CHR = SNP = NULL # silence R CMD check notes for data.table columns
  bfile <- file.path(thousand_genomes_dir_for(build, check = TRUE), ancestry)
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

  # PLINK writes the matrix in .bim order; keep the panel's alleles to orient r to the GWAS effect allele
  bim_ref <- paste0(bfile, ".bim")
  all_bim <- data.table::fread(bim_ref, header = FALSE, select = c(1, 2, 4, 5, 6),
                               colClasses = list(character = c(1, 2, 5, 6)))
  colnames(all_bim) <- c("CHR", "SNP", "BP", "A1", "A2")
  all_bim <- all_bim[CHR == as.character(chr) & SNP %in% rsids]
  snp_order <- all_bim$SNP
  a1 <- all_bim$A1
  a2 <- all_bim$A2
  rm(all_bim)

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

  return(list(matrix = R, snps = snp_order, a1 = a1, a2 = a2))
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
  summary[, c("LOCUS", "LEAD_IN_TILE_CORE", "LD_MISMATCH_S", "N_LD_OUTLIERS_IN_CS", "LD_FLAG") := NULL]
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
      LEAD_RSID = if (is.null(snp_info$RSID)) NA_character_ else as.character(snp_info$RSID[lead]),
      LEAD_BP = as.numeric(snp_info$BP[lead]),
      LEAD_PIP = pip[lead],
      LOG_BF = if (!is.null(lbf) && length(lbf) >= lbf_row) as.numeric(lbf[lbf_row]) else NA_real_
    )
  })
  data.table::rbindlist(rows)
}