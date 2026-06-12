#' Run pairwise  colocalization on finemapped GWAS results
#'
#' For each pair of GWAS datasets, identifies overlapping finemapped signals
#' (lead SNPs within ±\code{overlap_kb} kb), then runs \code{coloc::coloc.bf_bf}
#' on shared SNPs using the per-signal LBF columns produced by the fine-mapping
#' pipeline.
#'
#' @param finemap_dirs named character vector of finemap output directories,
#'   one per GWAS. Names are used as trait labels in the output.
#' @param overlap_kb distance in kb to define overlapping signals (default 1000 = ±1 Mb)
#' @param p1 prior probability a SNP is associated with trait 1
#' @param p2 prior probability a SNP is associated with trait 2
#' @param p12 prior probability a SNP is associated with both traits
#' @param output_file path to write the combined coloc results TSV
#' @return tibble of coloc results for every overlapping signal pair



#' @export
run_bf_bf_coloc <- function(finemap_dirs,
                            overlap_kb = 1000,
                            p1 = 1e-4,
                            p2 = 1e-4,
                            p12 = 5e-6,
                            output_file = NULL) {

  trait_names <- names(finemap_dirs)
  if (is.null(trait_names) || any(trait_names == "")) {
    trait_names <- basename(finemap_dirs)
    names(finemap_dirs) <- trait_names
  }

  locus_data <- load_all_finemap_loci(finemap_dirs)

  if (length(locus_data) < 2) {
    message("Need at least 2 traits with finemapped loci for colocalization.")
    empty <- tibble::tibble(
      trait1 = character(), trait2 = character(),
      locus1 = character(), locus2 = character(),
      n_snps = integer(),
      PP.H0.abf = numeric(), PP.H1.abf = numeric(),
      PP.H2.abf = numeric(), PP.H3.abf = numeric(), PP.H4.abf = numeric()
    )
    if (!is.null(output_file)) vroom::vroom_write(empty, output_file)
    return(invisible(empty))
  }

  pairs <- utils::combn(trait_names, 2, simplify = FALSE)
  overlap_bp <- overlap_kb * 1000
  all_results <- list()

  for (pair in pairs) {
    t1 <- pair[1]
    t2 <- pair[2]
    loci1 <- locus_data[[t1]]
    loci2 <- locus_data[[t2]]

    for (l1 in loci1) {
      for (l2 in loci2) {
        if (l1$chr != l2$chr) next
        if (abs(l1$lead_bp - l2$lead_bp) > overlap_bp) next

        result <- coloc_bf_bf_for_locus_pair(l1, l2, t1, t2, p1, p2, p12)
        if (!is.null(result)) {
          all_results <- c(all_results, list(result))
        }
      }
    }
  }

  combined <- dplyr::bind_rows(all_results)

  if (nrow(combined) == 0) {
    message("No overlapping finemapped signals found across trait pairs.")
  } else {
    message(paste(" colocalization complete:", nrow(combined), "signal pair(s) tested."))
  }

  if (!is.null(output_file)) {
    vroom::vroom_write(combined, output_file)
  }

  return(invisible(combined))
}


#' Load all finemapped locus files from a set of directories
#' @keywords internal
load_all_finemap_loci <- function(finemap_dirs) {
  locus_data <- list()

  for (trait in names(finemap_dirs)) {
    dir_path <- finemap_dirs[[trait]]
    if (!dir.exists(dir_path)) {
      message(paste("Finemap directory not found for", trait, ":", dir_path))
      next
    }

    files <- list.files(dir_path, pattern = "_finemap\\.tsv\\.gz$", full.names = TRUE)
    if (length(files) == 0) {
      message(paste("No finemap locus files found in", dir_path))
      next
    }

    trait_loci <- list()
    for (f in files) {
      locus <- vroom::vroom(f, show_col_types = FALSE)
      lbf_cols <- grep("^LBF_[0-9]+$", colnames(locus), value = TRUE)
      if (length(lbf_cols) == 0) next

      locus_name <- sub("_finemap\\.tsv\\.gz$", "", basename(f))
      parts <- strsplit(locus_name, "_")[[1]]
      chr <- as.numeric(parts[1])
      bp <- as.numeric(parts[2])

      trait_loci <- c(trait_loci, list(list(
        data = locus,
        lbf_cols = lbf_cols,
        chr = chr,
        lead_bp = bp,
        locus_name = locus_name
      )))
    }

    if (length(trait_loci) > 0) {
      locus_data[[trait]] <- trait_loci
    }
  }

  return(locus_data)
}


#' Run coloc::coloc.bf_bf for one pair of overlapping loci
#' @keywords internal
coloc_bf_bf_for_locus_pair <- function(l1, l2, trait1_name, trait2_name,
                                       p1, p2, p12) {
  d1 <- l1$data
  d2 <- l2$data

  if (!"SNP" %in% colnames(d1) || !"SNP" %in% colnames(d2)) {
    return(NULL)
  }

  shared_snps <- intersect(d1$SNP, d2$SNP)
  if (length(shared_snps) < 2) return(NULL)

  d1 <- dplyr::filter(d1, SNP %in% shared_snps) |>
    dplyr::filter(!duplicated(SNP)) |>
    dplyr::arrange(SNP)
  d2 <- dplyr::filter(d2, SNP %in% shared_snps) |>
    dplyr::filter(!duplicated(SNP)) |>
    dplyr::arrange(SNP)

  shared_snps <- intersect(d1$SNP, d2$SNP)
  d1 <- dplyr::filter(d1, SNP %in% shared_snps)
  d2 <- dplyr::filter(d2, SNP %in% shared_snps)

  bf1 <- build_lbf_matrix(d1, l1$lbf_cols)
  bf2 <- build_lbf_matrix(d2, l2$lbf_cols)

  if (is.null(bf1) || is.null(bf2)) return(NULL)

  result <- tryCatch(
    coloc::coloc.bf_bf(bf1 = bf1, bf2 = bf2, p1 = p1, p2 = p2, p12 = p12),
    error = function(e) {
      message(paste("coloc.bf_bf failed for", trait1_name, l1$locus_name,
                    "vs", trait2_name, l2$locus_name, ":", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(result)) return(NULL)

  parse_pairwise_bf_bf_result(result, trait1_name, trait2_name,
                             l1$locus_name, l2$locus_name, length(shared_snps))
}


#' Build a named BF matrix from locus data
#' @keywords internal
build_lbf_matrix <- function(locus_df, lbf_cols) {
  lbf_cols <- intersect(lbf_cols, colnames(locus_df))
  if (length(lbf_cols) == 0) return(NULL)

  has_cs <- "CS" %in% colnames(locus_df)
  if (has_cs) {
    has_signal <- !all(is.na(locus_df$CS))
  } else {
    has_signal <- TRUE
  }
  if (!has_signal) return(NULL)

  bf_mat <- t(as.matrix(locus_df[, lbf_cols, drop = FALSE]))
  bf_mat[is.na(bf_mat)] <- 0
  colnames(bf_mat) <- locus_df$SNP
  rownames(bf_mat) <- lbf_cols
  return(bf_mat)
}


#' Parse the result of coloc.bf_bf for a pairwise GWAS comparison
#' @keywords internal
parse_pairwise_bf_bf_result <- function(result, trait1, trait2, locus1, locus2, n_snps) {
  if (!is.null(result$summary)) {
    summary_df <- result$summary
    if (is.data.frame(summary_df) || is.matrix(summary_df)) {
      rows <- list()
      for (i in seq_len(nrow(summary_df))) {
        row <- summary_df[i, ]
        rows <- c(rows, list(tibble::tibble(
          trait1 = trait1,
          trait2 = trait2,
          locus1 = locus1,
          locus2 = locus2,
          hit1 = as.character(row[["hit1"]]),
          hit2 = as.character(row[["hit2"]]),
          n_snps = as.integer(row[["nsnps"]]),
          PP.H0.abf = as.numeric(row[["PP.H0.abf"]]),
          PP.H1.abf = as.numeric(row[["PP.H1.abf"]]),
          PP.H2.abf = as.numeric(row[["PP.H2.abf"]]),
          PP.H3.abf = as.numeric(row[["PP.H3.abf"]]),
          PP.H4.abf = as.numeric(row[["PP.H4.abf"]])
        )))
      }
      return(dplyr::bind_rows(rows))
    }
  }

  tibble::tibble(
    trait1 = trait1,
    trait2 = trait2,
    locus1 = locus1,
    locus2 = locus2,
    hit1 = NA_character_,
    hit2 = NA_character_,
    n_snps = n_snps,
    PP.H0.abf = NA_real_,
    PP.H1.abf = NA_real_,
    PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_,
    PP.H4.abf = NA_real_
  )
}
