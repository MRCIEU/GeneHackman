#' Run SuSiE fine-mapping across all clumped loci in a GWAS
#'
#' For each lead SNP from plink --clump output, extracts a window, computes an
#' LD matrix from the 1000 Genomes reference panel via plink, and runs
#' susieR::susie_rss.  Writes one TSV per locus: all GWAS variants within the
#' genomic window around the lead SNP, with SuSiE Z-scores, credible-set
#' membership, and per-credible-set (independent signal) columns \code{LBF_1}, \code{LBF_2}, ...
#'
#' @param gwas filename or dataframe of a standardised GWAS
#' @param clumped_file plink --clump output file
#' @param ancestry ancestry code matching 1000 Genomes panel prefix (EUR, EAS, AFR, AMR, SAS)
#' @param default_n GWAS sample size when not inferrable from \code{gwas}; see
#' @param output_finemap_dir directory to write one
#'   \verb{<CHR>_<BP>_finemap.tsv.gz} per clumped locus.
#' @param completion_file path to a sentinel file written on successful completion
#'   (one line: expected lead count). If NULL, defaults to
#'   \verb{<output_finemap_dir>/finemap_complete.txt}.
#' @param window_kb half-width of the fine-mapping window in kb (default 1000 = ±1 Mb)
#' @param max_causal maximum number of causal signals per locus (SuSiE L, default 10)
#' @param coverage credible set coverage (default 0.95)
#' @param min_abs_corr minimum absolute correlation for credible sets (default 0.5)
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
                         min_abs_corr = 0.5) {

  gwas <- get_file_or_dataframe(gwas)
  data.table::setDT(gwas)

  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF", "N", "N_CASE", "N_CONTROL")
  for (col in intersect(numeric_cols, colnames(gwas))) {
    data.table::set(gwas, j = col, value = as.numeric(gwas[[col]]))
  }

  data.table::setkey(gwas, CHR, BP)

  if (!dir.exists(output_finemap_dir)) {
    dir.create(output_finemap_dir, recursive = TRUE)
  }

  lead_snps <- data.table::fread(clumped_file, select = c("SNP", "CHR", "BP"))
  if (nrow(lead_snps) == 0) {
    message("No clumped SNPs found; no per-locus finemap files written.")
    write_finemap_complete_marker(completion_file, 0L)
    empty <- data.table::data.table(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), CS = integer()
    )
    return(invisible(empty))
  }

  window_bp <- window_kb * 1000
  has_rsid <- "RSID" %in% colnames(gwas)
  has_n <- "N" %in% colnames(gwas)

  n_loci <- nrow(lead_snps)

  # Pre-extract per-locus subsets so workers don't need the full GWAS
  locus_subsets <- vector("list", n_loci)
  window_subsets <- vector("list", n_loci)
  for (i in seq_len(n_loci)) {
    lead_chr <- as.numeric(lead_snps$CHR[i])
    lead_bp <- as.numeric(lead_snps$BP[i])
    bp_lo <- lead_bp - window_bp
    bp_hi <- lead_bp + window_bp

    chr_rows <- gwas[.(lead_chr), nomatch = NULL]
    window_dt <- chr_rows[BP >= bp_lo & BP <= bp_hi]
    window_subsets[[i]] <- window_dt
    locus_subsets[[i]] <- window_dt[!is.na(BETA) & !is.na(SE) & SE > 0]
  }
  rm(gwas); gc(verbose = FALSE)

  process_one_locus <- function(i) {
    lead_rsid <- lead_snps$SNP[i]
    lead_chr <- as.numeric(lead_snps$CHR[i])
    lead_bp <- as.numeric(lead_snps$BP[i])

    tryCatch({
      locus_gwas <- locus_subsets[[i]]
      window_gwas <- window_subsets[[i]]

      if (nrow(locus_gwas) < 2L) {
        message(paste("Skipping locus", lead_rsid, "- fewer than 2 SNPs in window"))
        return(NULL)
      }

      if (!has_rsid || all(is.na(locus_gwas$RSID))) {
        message(paste("Skipping locus", lead_rsid, "- no RSIDs available for LD computation"))
        return(NULL)
      }
      locus_gwas <- locus_gwas[!is.na(RSID) & nchar(RSID) > 0L]

      ld_result <- compute_ld_matrix(locus_gwas$RSID, lead_chr, ancestry)
      if (is.null(ld_result)) {
        message(paste("Skipping locus", lead_rsid, "- LD matrix computation failed"))
        return(NULL)
      }

      shared_rsids <- intersect(locus_gwas$RSID, ld_result$snps)
      if (length(shared_rsids) < 2L) {
        rm(ld_result); gc(verbose = FALSE)
        message(paste("Skipping locus", lead_rsid, "- too few shared SNPs between GWAS and LD panel"))
        return(NULL)
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

      if (nrow(locus_lbf) == 0L) return(NULL)

      lbf_nm <- grep("^LBF_[0-9]+$", names(locus_lbf), value = TRUE)
      finemap_cols <- intersect(c("RSID", "Z", "CS", lbf_nm), names(locus_lbf))
      finemap_join <- locus_lbf[, ..finemap_cols]

      out_gwas <- merge(window_gwas, finemap_join, by = "RSID", all.x = TRUE, sort = FALSE)
      rm(window_gwas, finemap_join)

      lead_chr_bp <- paste0(
        lead_chr, "_",
        format(as.numeric(lead_bp), scientific = FALSE, trim = TRUE, digits = 20)
      )
      safe_locus <- gsub("[^A-Za-z0-9._-]+", "_", lead_chr_bp)
      out_file <- file.path(output_finemap_dir, paste0(safe_locus, "_finemap.tsv.gz"))
      data.table::fwrite(out_gwas, out_file, sep = "\t", compress = "gzip")
      rm(out_gwas)

      locus_lbf
    }, error = function(e) {
      message(paste("Error processing locus", lead_rsid, ":", conditionMessage(e)))
      return(NULL)
    })
  }

  number_of_simultaneous_jobs <- 5L

  locus_results <- parallel::mclapply(
    seq_len(n_loci),
    process_one_locus,
    mc.cores = number_of_simultaneous_jobs
  )
  rm(locus_subsets, window_subsets); gc(verbose = FALSE)

  is_valid <- vapply(locus_results, function(x) {
    is.data.frame(x) || data.table::is.data.table(x)
  }, FUN.VALUE = logical(1))

  n_errors <- sum(vapply(locus_results, inherits, "try-error", FUN.VALUE = logical(1)))
  if (n_errors > 0L) {
    message(paste(n_errors, "locus worker(s) returned errors"))
  }

  all_lbf <- locus_results[is_valid]
  rm(locus_results)

  if (length(all_lbf) == 0L) {
    combined_lbf <- data.table::data.table(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), CS = integer()
    )
  } else {
    combined_lbf <- data.table::rbindlist(all_lbf, fill = TRUE)
  }
  rm(all_lbf); gc(verbose = FALSE)

  if (nrow(combined_lbf) == 0L) {
    message("No loci produced SuSiE results.")
  }

  message(paste("Fine-mapping complete:",
    nrow(lead_snps), "clumped loci,",
    nrow(combined_lbf), "LD-matched SNPs with finemap stats,",
    "outputs in", output_finemap_dir))

  write_finemap_complete_marker(completion_file, nrow(lead_snps))

  return(invisible(combined_lbf))
}


#' Write a completion sentinel file (one line: expected lead count) for Snakemake.
#' @keywords internal
write_finemap_complete_marker <- function(completion_file, n_loci) {
  dir_path <- dirname(completion_file)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  writeLines(as.character(as.integer(n_loci)), completion_file)
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
    "--threads", as.character(number_of_cpus_available),
    "--bfile", bfile,
    "--chr", chr,
    "--extract", snp_file,
    "--r square",
    "--out", out_prefix
  )
  exit_code <- run_system(cmd, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(out_prefix, ".ld")
  snp_order_file <- paste0(out_prefix, ".nosex")

  if (exit_code != 0 || !file.exists(ld_file)) {
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
    all_bim <- all_bim[CHR == chr & SNP %chin% rsids]
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
#' @return data.table with SNP, CHR, BP, RSID, Z, CS, and \code{LBF_1}, \code{LBF_2}, ... per signal
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
      Z = numeric(), CS = integer()
    ))
  }

  lbf_wide <- susie_lbf_columns(fitted, p)

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
    CS = cs_membership
  )

  if (ncol(lbf_wide) > 0L) {
    result <- cbind(result, lbf_wide)
  }
  result
}
