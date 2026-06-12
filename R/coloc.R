#' run_coloc_on_list_of_datasets: run coloc on a list of datasets
#' @param first_gwas_list list of first GWAS files
#' @param second_gwas_list list of second GWAS files
#' @param exposure_name_list list of exposure names
#' @param chr_list list of chromosomes
#' @param bp_list list of base pair positions
#' @param range range to filter in base pairs
#' @param default_n default sample size
#' @param output_file file to save the results
#' @return tibble of coloc results



#' @export
run_coloc_on_list_of_datasets <- function(first_gwas_list=list(),
                                          second_gwas_list=list(),
                                          exposure_name_list=list(),
                                          chr_list=list(),
                                          bp_list=list(),
                                          range=500000,
                                          default_n=NA,
                                          output_file) {
  input_lengths <- c(length(first_gwas_list), length(second_gwas_list), length(chr_list), length(bp_list))
  if (var(input_lengths) != 0) {
    stop("Error: Input lengths are not equal")
  }

  data_for_coloc <- tibble::tibble(
      first_gwas = first_gwas_list,
      second_gwas = second_gwas_list,
      chr = chr_list,
      bp = bp_list,
      exposure_name = exposure_name_list
  )

  coloc_results <- apply(data_for_coloc,1, function(row) {
    coloc_columns <- c("SNP", "CHR", "BP", "P", "SE", "N", "EAF")
    first_gwas <- get_file_or_dataframe(row[['first_gwas']], columns = coloc_columns)
    second_gwas <- get_file_or_dataframe(row[['second_gwas']], columns = coloc_columns)

    coloc_analysis(first_gwas, second_gwas, row[['exposure_name']], as.numeric(row[['chr']]), as.numeric(row[['bp']]), range, default_n=default_n)
  }) |> dplyr::bind_rows()

  if (!missing(output_file)) {
    vroom::vroom_write(coloc_results, output_file)
  }
  return(coloc_results)
}

#' Run  colocalization on significant QTL MR results using finemapped data
#'
#' For each significant MR result, loads the finemapped GWAS locus that covers
#' the MR hit, computes LBF scores for the QTL data via \code{convert_z_to_lbf},
#' and runs \code{coloc::coloc.bf_bf}.
#'
#' @param mr_results_file MR results file
#' @param finemap_dir directory containing finemapped GWAS locus files
#' @param qtl_dataset QTL dataset name (metabrain or eqtlgen)
#' @param study_type study type for LBF conversion ("continuous" or "categorical")
#' @param exposures optional character vector of exposures to filter
#' @param default_n default sample size for QTL data
#' @param output_file output file path




#' @export
run_coloc_on_qtl_mr_results <- function(mr_results_file,
                                        finemap_dir,
                                        qtl_dataset,
                                        study_type = "continuous",
                                        exposures = c(),
                                        default_n = NA,
                                        output_file) {

  mr_results <- get_file_or_dataframe(mr_results_file) |>
    dplyr::filter(p.adjusted < 0.05)

  if (length(exposures) > 0) {
    mr_results <- subset(mr_results, exposure %in% exposures)
  }

  empty_result <- tibble::tibble(
    exposure = character(),
    locus = character(),
    hit1 = character(),
    hit2 = character(),
    n_snps = integer(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric()
  )

  if (nrow(mr_results) == 0) {
    vroom::vroom_write(empty_result, output_file)
    return(invisible(empty_result))
  }

  coloc_results <- apply(mr_results, 1, function(mr_result) {
    qtl_gwas <- load_qtl_gwas_for_mr_result(mr_result, qtl_dataset)
    if (is.null(qtl_gwas)) return(NULL)

    chr <- as.numeric(mr_result[["CHR"]])
    bp <- as.numeric(mr_result[["BP"]])
    exposure_name <- mr_result[["EXPOSURE"]]

    coloc_bf_bf_qtl_analysis(
      finemap_dir = finemap_dir,
      qtl_gwas = qtl_gwas,
      exposure_name = exposure_name,
      chr = chr,
      bp = bp,
      study_type = study_type,
      default_n = default_n
    )
  }) |> dplyr::bind_rows()

  if (nrow(coloc_results) == 0) coloc_results <- empty_result

  vroom::vroom_write(coloc_results, output_file)
  return(invisible(coloc_results))
}


#' Load the QTL GWAS file for an MR result row
#' @keywords internal
load_qtl_gwas_for_mr_result <- function(mr_result, qtl_dataset) {
  if (qtl_dataset == qtl_datasets$metabrain) {
    outcome <- unlist(strsplit(mr_result[["outcome"]], "_"))
    brain_region <- outcome[1]
    ancestry <- toupper(outcome[2])
    qtl_gwas_file <- paste0(metabrain_gwas_dir, "/", brain_region, "/",
                            mr_result[["EXPOSURE"]], "_", ancestry, ".tsv.gz")
  } else if (qtl_dataset == qtl_datasets$eqtlgen) {
    qtl_gwas_file <- paste0(eqtlgen_gwas_dir, "/", mr_result[["outcome"]], "/",
                            mr_result[["EXPOSURE"]], ".tsv.gz")
  } else {
    stop("Error: qtl dataset not supported for coloc right now")
  }

  if (!file.exists(qtl_gwas_file)) {
    message(paste("QTL GWAS file not found:", qtl_gwas_file))
    return(NULL)
  }
  get_file_or_dataframe(qtl_gwas_file)
}


#' Run  coloc between a finemapped GWAS locus and QTL data
#'
#' Finds the finemapped locus file covering the region around chr:bp,
#' extracts the GWAS LBF matrix, computes QTL LBF via \code{convert_z_to_lbf},
#' and runs \code{coloc::coloc.bf_bf}.
#'
#' @param finemap_dir finemap output directory containing locus TSV files
#' @param qtl_gwas QTL GWAS data frame
#' @param exposure_name exposure label
#' @param chr chromosome of the MR hit
#' @param bp base position of the MR hit
#' @param study_type study type for LBF conversion
#' @param default_n default sample size
#' @param range_bp overlap window in bp (default 1e6 = ±1Mb)
#' @return tibble with coloc results or NULL



#' @export
coloc_bf_bf_qtl_analysis <- function(finemap_dir,
                                     qtl_gwas,
                                     exposure_name,
                                     chr,
                                     bp,
                                     study_type = "continuous",
                                     default_n = NA,
                                     range_bp = 1e6) {

  finemap_locus <- find_finemap_locus_for_region(finemap_dir, chr, bp, range_bp)
  if (is.null(finemap_locus)) {
    message(paste("No finemapped locus found near", chr, ":", bp))
    return(NULL)
  }

  lbf_cols <- grep("^LBF_[0-9]+$", colnames(finemap_locus), value = TRUE)
  if (length(lbf_cols) == 0) {
    message(paste("No LBF columns in finemapped locus for", chr, ":", bp))
    return(NULL)
  }

  numeric_cols <- c("CHR", "BP", "BETA", "SE", "P", "EAF", "N")
  for (col in intersect(numeric_cols, colnames(qtl_gwas))) {
    qtl_gwas[[col]] <- as.numeric(qtl_gwas[[col]])
  }

  qtl_gwas <- dplyr::filter(qtl_gwas,
    !is.na(BETA) & !is.na(SE) & SE > 0 &
    !is.na(EAF) & EAF > 0 & EAF < 1
  )

  qtl_region <- gwas_region(qtl_gwas, chr, bp, range_bp)
  if (nrow(qtl_region) < 2) {
    message(paste("Too few QTL SNPs in region around", chr, ":", bp))
    return(NULL)
  }

  shared_snps <- intersect(finemap_locus$SNP, qtl_region$SNP)
  if (length(shared_snps) < 2) {
    message(paste("Too few shared SNPs for", exposure_name, "at", chr, ":", bp))
    return(NULL)
  }

  gwas_sub <- dplyr::filter(finemap_locus, SNP %in% shared_snps) |>
    dplyr::filter(!duplicated(SNP)) |>
    dplyr::arrange(SNP)
  qtl_sub <- dplyr::filter(qtl_region, SNP %in% shared_snps) |>
    dplyr::filter(!duplicated(SNP)) |>
    dplyr::arrange(SNP)

  shared_snps <- intersect(gwas_sub$SNP, qtl_sub$SNP)
  gwas_sub <- dplyr::filter(gwas_sub, SNP %in% shared_snps)
  qtl_sub <- dplyr::filter(qtl_sub, SNP %in% shared_snps)

  gwas_bf <- build_bf_matrix(gwas_sub, lbf_cols)
  if (is.null(gwas_bf)) return(NULL)

  qtl_n <- `if`("N" %in% colnames(qtl_sub) && !is.na(as.numeric(qtl_sub$N[1])),
                as.numeric(qtl_sub$N[1]), default_n)
  if (is.na(qtl_n)) {
    message(paste("Cannot determine N for QTL data for", exposure_name))
    return(NULL)
  }

  z_scores <- qtl_sub$BETA / qtl_sub$SE
  qtl_lbf <- convert_z_to_lbf(
    z = z_scores,
    se = qtl_sub$SE,
    eaf = qtl_sub$EAF,
    sample_size = qtl_n,
    study_type = study_type
  )
  names(qtl_lbf) <- qtl_sub$SNP

  result <- tryCatch(
    coloc::coloc.bf_bf(bf1 = gwas_bf, bf2 = qtl_lbf),
    error = function(e) {
      message(paste("coloc.bf_bf failed for", exposure_name, "at", chr, ":", bp,
                    ":", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(result)) return(NULL)

  locus_name <- paste0(chr, "_", bp)
  parsed <- parse_bf_bf_result(result, exposure_name, locus_name, length(shared_snps))
  return(parsed)
}


#' Build a named BF matrix from a finemapped locus data frame
#' @keywords internal
build_bf_matrix <- function(locus_df, lbf_cols) {
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


#' Find the finemapped locus file that covers a genomic region
#' @keywords internal
find_finemap_locus_for_region <- function(finemap_dir, chr, bp, range_bp = 1e6) {
  if (!dir.exists(finemap_dir)) return(NULL)

  files <- list.files(finemap_dir, pattern = "_finemap\\.tsv\\.gz$", full.names = TRUE)
  if (length(files) == 0) return(NULL)

  best_file <- NULL
  best_dist <- Inf

  for (f in files) {
    locus_name <- sub("_finemap\\.tsv\\.gz$", "", basename(f))
    parts <- strsplit(locus_name, "_")[[1]]
    locus_chr <- as.numeric(parts[1])
    locus_bp <- as.numeric(parts[2])

    if (locus_chr != chr) next
    dist <- abs(locus_bp - bp)
    if (dist <= range_bp && dist < best_dist) {
      best_dist <- dist
      best_file <- f
    }
  }

  if (is.null(best_file)) return(NULL)
  vroom::vroom(best_file, show_col_types = FALSE)
}


#' Parse coloc.bf_bf result for QTL analysis
#' @keywords internal
parse_bf_bf_result <- function(result, exposure_name, locus_name, n_snps) {
  if (!is.null(result$summary)) {
    summary_df <- result$summary
    if (is.data.frame(summary_df) || is.matrix(summary_df)) {
      rows <- list()
      for (i in seq_len(nrow(summary_df))) {
        row <- summary_df[i, ]
        rows <- c(rows, list(tibble::tibble(
          exposure = exposure_name,
          locus = locus_name,
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
    exposure = exposure_name,
    locus = locus_name,
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

#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas second gwas to be run through coloc
#' @param exposure_name name of exposure to perform coloc on
#' @param chr chromosome to perform coloc on
#' @param bp base pair position to perform coloc on
#' @param range range to filter in base pairs
#' @param default_n default sample size
#' @return tibble of coloc results (h0 - h4)


#' @export
coloc_analysis <- function(first_gwas, second_gwas, exposure_name, chr=NA, bp=NA, range=NA, default_n=NA) {
  numeric_columns <- c("P", "SE", "EAF")
  first_gwas[,numeric_columns] <- lapply(first_gwas[,numeric_columns,drop=FALSE], as.numeric)
  second_gwas[,numeric_columns] <- lapply(second_gwas[,numeric_columns,drop=FALSE], as.numeric)

  first_gwas <- dplyr::filter(first_gwas, EAF > 0 & EAF < 1)
  second_gwas <- dplyr::filter(second_gwas, EAF > 0 & EAF < 1)

  if (!is.na(chr) & !is.na(bp) & !is.na(range)) {
    first_gwas <- gwas_region(first_gwas, chr, bp, range)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
  }

  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  first_gwas <- harmonised_gwases[[1]]
  second_gwas <- harmonised_gwases[[2]]

  first_n <- `if`("N" %in% colnames(first_gwas) && !is.na(as.numeric(first_gwas$N[1])), as.numeric(first_gwas$N[1]), default_n)
  second_n <- `if`("N" %in% colnames(second_gwas) && !is.na(as.numeric(second_gwas$N[1])), as.numeric(second_gwas$N[1]), default_n)

  if (is.na(first_n) || is.na(second_n)) {
    stop("Error: N (sample size) must be present to complete coloc analysis.  Please include N in gwas or populate default_n param")
  }

  first_coloc_dataset <- list(
    pvalues = first_gwas$P,
    N = first_n,
    varbeta = first_gwas$SE^2,
    type = "quant",
    snp = first_gwas$SNP,
    MAF = first_gwas$EAF
  )

  second_coloc_dataset <- list(
    pvalues = second_gwas$P,
    N = second_n,
    varbeta = second_gwas$SE^2,
    type = "quant",
    snp = second_gwas$SNP,
    MAF = second_gwas$EAF
  )

  result <- coloc::coloc.abf(dataset1 = first_coloc_dataset, dataset2 = second_coloc_dataset)
  coloc_results <- tibble::tribble(
    ~exposure, ~h0, ~h1, ~h2, ~h3, ~h4,
    exposure_name, result$summary[2], result$summary[3], result$summary[4], result$summary[5], result$summary[6]
  )

  return(coloc_results)
}
