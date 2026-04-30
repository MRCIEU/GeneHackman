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
#'   \verb{<CHR>_<BP>_finemap.tsv.gz} per clumped locus
#' @param window_kb half-width of the fine-mapping window in kb (default 1000 = ±1 Mb)
#' @param max_causal maximum number of causal signals per locus (SuSiE L, default 10)
#' @param coverage credible set coverage (default 0.95)
#' @param min_abs_corr minimum absolute correlation for credible sets (default 0.5)
#' @details Loci are fine-mapped in parallel with [parallel::mclapply()] (`mc.cores = 2`).
#'   On Windows this runs sequentially. Each job writes its own output file.
#' @return invisibly, the combined per-SNP finemap tibble (LD-matched SNPs only)




#' @export
finemap_gwas <- function(gwas,
                         clumped_file,
                         ancestry,
                         default_n,
                         output_finemap_dir,
                         window_kb = 1000,
                         max_causal = 10,
                         coverage = 0.95,
                         min_abs_corr = 0.5) {

  gwas <- get_file_or_dataframe(gwas)
  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF", "N", "N_CASE", "N_CONTROL")
  for (col in intersect(numeric_cols, colnames(gwas))) {
    gwas[[col]] <- as.numeric(gwas[[col]])
  }

  if (!dir.exists(output_finemap_dir)) {
    dir.create(output_finemap_dir, recursive = TRUE)
  }

  lead_snps <- data.table::fread(clumped_file, select = c("SNP", "CHR", "BP"))
  if (nrow(lead_snps) == 0) {
    message("No clumped SNPs found; no per-locus finemap files written.")
    empty <- tibble::tibble(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), CS = integer()
    )
    return(invisible(empty))
  }

  window_bp <- window_kb * 1000

  process_one_locus <- function(i) {
    lead_rsid <- lead_snps$SNP[i]
    lead_chr <- as.numeric(lead_snps$CHR[i])
    lead_bp <- as.numeric(lead_snps$BP[i])

    locus_gwas <- dplyr::filter(gwas,
      CHR == lead_chr &
      BP >= (lead_bp - window_bp) &
      BP <= (lead_bp + window_bp) &
      !is.na(BETA) & !is.na(SE) & SE > 0
    )

    if (nrow(locus_gwas) < 2) {
      message(paste("Skipping locus", lead_rsid, "- fewer than 2 SNPs in window"))
      return(NULL)
    }

    if (!"RSID" %in% colnames(locus_gwas) || all(is.na(locus_gwas$RSID))) {
      message(paste("Skipping locus", lead_rsid, "- no RSIDs available for LD computation"))
      return(NULL)
    }
    locus_gwas <- dplyr::filter(locus_gwas, !is.na(RSID) & nchar(RSID) > 0)

    ld_result <- compute_ld_matrix(locus_gwas$RSID, lead_chr, ancestry)
    if (is.null(ld_result)) {
      message(paste("Skipping locus", lead_rsid, "- LD matrix computation failed"))
      return(NULL)
    }

    shared_rsids <- intersect(locus_gwas$RSID, ld_result$snps)
    if (length(shared_rsids) < 2) {
      message(paste("Skipping locus", lead_rsid, "- too few shared SNPs between GWAS and LD panel"))
      return(NULL)
    }

    locus_gwas <- dplyr::filter(locus_gwas, RSID %in% shared_rsids)
    locus_gwas <- locus_gwas[match(shared_rsids, locus_gwas$RSID), ]
    ld_idx <- match(shared_rsids, ld_result$snps)
    R <- ld_result$matrix[ld_idx, ld_idx]

    z_scores <- locus_gwas$BETA / locus_gwas$SE
    n <- `if`("N" %in% colnames(locus_gwas) && !is.na(as.numeric(locus_gwas$N[1])), as.numeric(locus_gwas$N[1]), default_n)
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

    if (nrow(locus_lbf) == 0) {
      return(NULL)
    }

    lbf_nm <- grep("^LBF_[0-9]+$", names(locus_lbf), value = TRUE)
    finemap_join <- dplyr::select(locus_lbf, dplyr::any_of(c(
      "RSID", "Z", "CS", lbf_nm
    )))

    out_gwas <- dplyr::filter(gwas,
      CHR == lead_chr,
      BP >= (lead_bp - window_bp),
      BP <= (lead_bp + window_bp)
    )
    out_gwas <- dplyr::left_join(out_gwas, finemap_join, by = "RSID")

    lead_chr_bp <- paste0(
      lead_chr, "_",
      format(as.numeric(lead_bp), scientific = FALSE, trim = TRUE, digits = 20)
    )
    safe_locus <- gsub("[^A-Za-z0-9._-]+", "_", lead_chr_bp)
    out_file <- file.path(output_finemap_dir, paste0(safe_locus, "_finemap.tsv.gz"))
    vroom::vroom_write(out_gwas, out_file)

    locus_lbf
  }

  n_loci <- nrow(lead_snps)
  locus_results <- parallel::mclapply(
    seq_len(n_loci),
    process_one_locus,
    mc.cores = 2L
  )
  all_lbf <- locus_results[!vapply(locus_results, is.null, FUN.VALUE = logical(1))]
  combined_lbf <- dplyr::bind_rows(all_lbf)

  if (nrow(combined_lbf) == 0) {
    message("No loci produced SuSiE results.")
  }

  message(paste("Fine-mapping complete:",
    nrow(lead_snps), "clumped loci,",
    nrow(combined_lbf), "LD-matched SNPs with finemap stats,",
    "outputs in", output_finemap_dir))

  return(invisible(combined_lbf))
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
    "--r square",
    "--out", out_prefix
  )
  exit_code <- run_system(cmd, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(out_prefix, ".ld")
  snp_order_file <- paste0(out_prefix, ".nosex")

  if (exit_code != 0 || !file.exists(ld_file)) {
    return(NULL)
  }

  bim_file <- paste0(out_prefix, ".bim")
  if (file.exists(bim_file)) {
    snp_order <- data.table::fread(bim_file, header = FALSE, select = 2)$V2
  } else {
    bim_ref <- paste0(bfile, ".bim")
    all_bim <- data.table::fread(bim_ref, header = FALSE, select = c(1, 2, 4))
    colnames(all_bim) <- c("CHR", "SNP", "BP")
    all_bim <- dplyr::filter(all_bim, CHR == chr & SNP %in% rsids)
    all_bim <- dplyr::arrange(all_bim, BP)
    snp_order <- all_bim$SNP
  }

  ld_raw <- data.table::fread(ld_file, header = FALSE)
  R <- as.matrix(ld_raw)

  if (nrow(R) != length(snp_order) || ncol(R) != length(snp_order)) {
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
    return(tibble::tibble())
  }

  n_row <- nrow(lv)
  cs_list <- fitted$sets$cs
  cs_idx <- fitted$sets$cs_index
  col_vecs <- list()

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

  tibble::as_tibble(col_vecs)
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
#' @return tibble with SNP, CHR, BP, RSID, Z, CS, and \code{LBF_1}, \code{LBF_2}, ... per signal
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
    return(tibble::tibble(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), CS = integer()
    ))
  }

  lbf_wide <- susie_lbf_columns(fitted, p)

  cs_membership <- rep(NA_integer_, p)
  if (!is.null(fitted$sets) && !is.null(fitted$sets$cs)) {
    for (cs_idx in seq_along(fitted$sets$cs)) {
      snp_indices <- fitted$sets$cs[[cs_idx]]
      cs_membership[snp_indices] <- cs_idx
    }
  }

  result <- tibble::tibble(
    SNP = snp_info$SNP,
    CHR = snp_info$CHR,
    BP = snp_info$BP,
    RSID = snp_info$RSID,
    Z = z_scores,
    CS = cs_membership
  )
  dplyr::bind_cols(result, lbf_wide)
}
