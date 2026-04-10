#' Run SuSiE fine-mapping across all clumped loci in a GWAS
#'
#' For each lead SNP from plink --clump output, extracts a window, computes an
#' LD matrix from the 1000 Genomes reference panel via plink, and runs
#' susieR::susie_rss.  Produces per-SNP log Bayes factors and a filtered GWAS
#' containing only credible-set SNPs.
#'
#' @param gwas filename or dataframe of a standardised GWAS
#' @param clumped_file plink --clump output file
#' @param ancestry ancestry code matching 1000 Genomes panel prefix (EUR, EAS, AFR, AMR, SAS)
#' @param n GWAS sample size
#' @param output_lbf_file path to write combined per-locus LBF table
#' @param output_credible_set_file path to write credible-set-filtered GWAS
#' @param window_kb fine-mapping window half-width in kb (default 500)
#' @param max_causal maximum number of causal signals per locus (SuSiE L, default 10)
#' @param coverage credible set coverage (default 0.95)
#' @param min_abs_corr minimum absolute correlation for credible sets (default 0.5)
#' @return invisibly, the combined LBF tibble
#' @import data.table
#' @import dplyr
#' @import vroom
#' @import tibble
#' @export
finemap_gwas <- function(gwas,
                         clumped_file,
                         ancestry,
                         n,
                         output_lbf_file,
                         output_credible_set_file,
                         window_kb = 500,
                         max_causal = 10,
                         coverage = 0.95,
                         min_abs_corr = 0.5) {

  gwas <- get_file_or_dataframe(gwas)
  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF")
  for (col in intersect(numeric_cols, colnames(gwas))) {
    gwas[[col]] <- as.numeric(gwas[[col]])
  }

  lead_snps <- data.table::fread(clumped_file, select = c("SNP", "CHR", "BP"))
  if (nrow(lead_snps) == 0) {
    message("No clumped SNPs found; writing empty output files.")
    empty_lbf <- tibble::tibble(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), LBF = numeric(), CS = integer(), LEAD_SNP = character()
    )
    vroom::vroom_write(empty_lbf, output_lbf_file)
    vroom::vroom_write(gwas[0, ], output_credible_set_file)
    return(invisible(empty_lbf))
  }

  window_bp <- window_kb * 1000
  all_lbf <- list()

  for (i in seq_len(nrow(lead_snps))) {
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
      next
    }

    if (!"RSID" %in% colnames(locus_gwas) || all(is.na(locus_gwas$RSID))) {
      message(paste("Skipping locus", lead_rsid, "- no RSIDs available for LD computation"))
      next
    }
    locus_gwas <- dplyr::filter(locus_gwas, !is.na(RSID) & nchar(RSID) > 0)

    ld_result <- compute_ld_matrix(locus_gwas$RSID, lead_chr, ancestry)
    if (is.null(ld_result)) {
      message(paste("Skipping locus", lead_rsid, "- LD matrix computation failed"))
      next
    }

    shared_rsids <- intersect(locus_gwas$RSID, ld_result$snps)
    if (length(shared_rsids) < 2) {
      message(paste("Skipping locus", lead_rsid, "- too few shared SNPs between GWAS and LD panel"))
      next
    }

    locus_gwas <- dplyr::filter(locus_gwas, RSID %in% shared_rsids)
    locus_gwas <- locus_gwas[match(shared_rsids, locus_gwas$RSID), ]
    ld_idx <- match(shared_rsids, ld_result$snps)
    R <- ld_result$matrix[ld_idx, ld_idx]

    z_scores <- locus_gwas$BETA / locus_gwas$SE

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
    all_lbf[[i]] <- locus_lbf
  }

  combined_lbf <- dplyr::bind_rows(all_lbf)

  if (nrow(combined_lbf) == 0) {
    message("No loci produced SuSiE results.")
    combined_lbf <- tibble::tibble(
      SNP = character(), CHR = numeric(), BP = numeric(), RSID = character(),
      Z = numeric(), LBF = numeric(), CS = integer(), LEAD_SNP = character()
    )
  }

  vroom::vroom_write(combined_lbf, output_lbf_file)

  cs_snps <- dplyr::filter(combined_lbf, !is.na(CS))
  credible_set_gwas <- dplyr::filter(gwas, SNP %in% cs_snps$SNP)
  vroom::vroom_write(credible_set_gwas, output_credible_set_file)

  message(paste("Fine-mapping complete:",
    nrow(lead_snps), "loci,",
    nrow(combined_lbf), "SNPs with LBFs,",
    nrow(credible_set_gwas), "SNPs in credible sets"))

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

  # plink --r square writes the SNP order into .ld in the same order as the
  # extracted variants. Read the .bim-like info from the log or use the
  # nosex/fam to get the order. Safer: read the pruned bim.
  bim_file <- paste0(out_prefix, ".bim")
  if (file.exists(bim_file)) {
    snp_order <- data.table::fread(bim_file, header = FALSE, select = 2)$V2
  } else {
    # Fallback: the SNPs in .ld follow input extraction order; plink sorts by
    # genomic position so we cannot assume our input order.  Read the log
    # to find the count and use the filtered set.
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


#' Wrapper for susieR::susie_rss (mockable in tests).
#' @keywords internal
run_susie_rss_impl <- function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
  susieR::susie_rss(
    z = z,
    R = R,
    n = n,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    verbose = verbose
  )
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
#' @return tibble with SNP, CHR, BP, RSID, Z, LBF, CS, LEAD_SNP
run_susie_for_locus <- function(z_scores, ld_matrix, snp_info, n, lead_snp,
                                max_causal = 10, coverage = 0.95,
                                min_abs_corr = 0.5) {

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
      Z = numeric(), LBF = numeric(), CS = integer(), LEAD_SNP = character()
    ))
  }

  # Per-SNP log Bayes factor: sum across all L effects
  if (!is.null(fitted$lbf_variable)) {
    per_snp_lbf <- colSums(fitted$lbf_variable)
  } else {
    per_snp_lbf <- rep(NA_real_, length(z_scores))
  }

  # Map each SNP to its credible set (if any)
  cs_membership <- rep(NA_integer_, length(z_scores))
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
    LBF = per_snp_lbf,
    CS = cs_membership,
    LEAD_SNP = lead_snp
  )

  return(result)
}
