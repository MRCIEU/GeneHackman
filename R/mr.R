#' @keywords internal
mr_list_files <- function(dir, pattern, full.names = FALSE) {
  list.files(dir, pattern = pattern, full.names = full.names)
}

#' Check that a QTL dataset is one GeneHackman provides
#' @keywords internal
check_qtl_dataset <- function(dataset) {
  if (length(dataset) != 1 || !dataset %in% qtl_datasets) {
    stop("Unknown QTL dataset '", paste(dataset, collapse = ","), "'; options: ", paste(qtl_datasets, collapse = ", "))
  }
  invisible(dataset)
}

#' MR of an outcome GWAS against the cis QTLs of a QTL dataset
#'
#' The exposures are the dataset's genes or proteins, instrumented by their cis QTL top hits
#' (\verb{QTL_DATA_DIR/<dataset>/<build>/top_hits/cis.tsv.gz}).
#' @param gwas_filename standardised outcome GWAS
#' @param dataset one of \code{qtl_datasets}
#' @param exposures optional exposure names to keep
#' @param results_output output file of MR results
#' @param build reference build of the GWAS; the QTL files of the same build are used
perform_mr_on_qtl_dataset <- function(gwas_filename, dataset, exposures=c(), results_output,
                                      build = reference_builds$GRCh38) {
  check_qtl_dataset(dataset)
  top_hits_dir <- qtl_dataset_dir(dataset, "top_hits", build)
  qtl_files <- mr_list_files(top_hits_dir, pattern = paste0("^", qtl_subcategory, "\\.tsv(\\.gz)?$"), full.names = T)
  if (length(qtl_files) == 0) {
    stop("No ", qtl_subcategory, " QTL top hits for ", dataset, " in ", top_hits_dir)
  }
  run_mr_on_qtl_data(gwas_filename, qtl_files = qtl_files, results_output = results_output, exposures = exposures)
}




run_mr_on_qtl_data <- function(gwas_filename, qtl_files, results_output, exposures=c()) {
  all_qtl_mr_results <- lapply(qtl_files, function(qtl_file) {
    qtl_exposure <- TwoSampleMR::read_exposure_data(qtl_file,
                                                    snp_col="SNP",
                                                    effect_allele_col="EA",
                                                    other_allele_col="OA",
                                                    eaf_col="EAF",
                                                    beta_col="BETA",
                                                    se_col="SE",
                                                    pval_col="P",
                                                    phenotype_col="EXPOSURE",
                                                    sep="\t"
    )

    #Don't need to clump because the qtl data sets are already only top hits
    #qtl_exposure <- TwoSampleMR::clump_data(qtl_exposure, pop=ancestry)

    gwas_outcome_data <- TwoSampleMR::read_outcome_data(gwas_filename,
                                                        snp_col = "SNP",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "EA",
                                                        other_allele_col = "OA",
                                                        pval_col = "P",
                                                        eaf_col = "EAF",
                                                        snps = qtl_exposure$SNP,
                                                        sep="\t",
    )
    gwas_outcome_data$outcome <- file_prefix(qtl_file)

    harmonised_data <- TwoSampleMR::harmonise_data(qtl_exposure, gwas_outcome_data)
    mr_results <- TwoSampleMR::mr(harmonised_data, method_list=c("mr_wald_ratio", "mr_ivw"))
    if (length(exposures) > 0) {
      mr_results <- subset(mr_results, exposure %in% exposures)
    }
    harmonised_data$SNP <- toupper(harmonised_data$SNP)

    # One row per exposure: an exposure can have several instruments (one cis lead per LD block), and
    # mr() then gives one IVW row for it.  SNP is the strongest instrument used, F_STAT the weakest one's.
    instruments_used <- harmonised_data |>
      dplyr::filter(mr_keep) |>
      dplyr::distinct(SNP, EXPOSURE = exposure)
    qtl_dataset <- vroom::vroom(qtl_file, show_col_types=F) |>
      calculate_f_statistic()
    relevant_qtl_results <- qtl_dataset |>
      dplyr::mutate(SNP = toupper(SNP)) |>
      dplyr::inner_join(instruments_used, by = c("SNP", "EXPOSURE")) |>
      dplyr::arrange(P) |>
      dplyr::summarise(SNP = dplyr::first(SNP), F_STAT = min(F_STAT), .by = EXPOSURE)

    mr_results <- merge(mr_results, relevant_qtl_results, by.x="exposure", by.y="EXPOSURE")
    return(mr_results)
  }) |> dplyr::bind_rows()


  all_qtl_mr_results <- tidyr::separate(all_qtl_mr_results, col = "SNP", into = c("CHR", "BP", "EA", "OA"), sep = "[:_]", remove = F) |>
    dplyr::rename(BETA="b", SE="se", P="pval", EXPOSURE="exposure") |>
    dplyr::mutate(p.adjusted = p.adjust(P, "fdr"))

  vroom::vroom_write(all_qtl_mr_results, results_output)
}

#'
#'


compare_interesting_mr_results <- function(pqtl_mr_results, forest_plot_output_file) {
  interesting_pqtl_mr_results <- vroom::vroom(pqtl_mr_results, show_col_types=F) |>
    dplyr::filter(p.adjusted < 0.05) |>
    dplyr::select(EXPOSURE, dplyr::everything())

  grouped_forest_plot(interesting_pqtl_mr_results,
                      "Comparing MR Results across pQTL data sets",
                      "outcome",
                      forest_plot_output_file
  )
}
