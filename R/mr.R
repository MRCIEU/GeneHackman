#'
perform_mr_on_metabrain_datasets <- function(gwas_filename, ancestry="EUR", subcategory=NULL, exposures=c(), results_output) {
  file_pattern <- paste0(tolower(subcategory), "_", tolower(ancestry))
  metabrain_top_hits <- list.files(metabrain_top_hits_dir, pattern = file_pattern, full.names = T)
  print(metabrain_top_hits)
  print(exposures)

  run_mr_on_qtl_data(gwas_filename, qtl_files = metabrain_top_hits, results_output = results_output, exposures = exposures)
}

perform_mr_on_eqtlgen_datasets <- function(gwas_filename, subcategory=NULL, exposures=c(), results_output) {
  file_pattern <- tolower(subcategory)
  eqtlgen_top_hits <- list.files(eqtlgen_top_hits_dir, pattern = file_pattern, full.names = T)

  run_mr_on_qtl_data(gwas_filename, results_output = results_output, qtl_files = eqtlgen_top_hits, exposures = exposures)
}


#'
perform_mr_on_pqtl_datasets <- function(gwas_filename, results_output, qtl_dir=pqtl_top_hits_dir, subcategory=None, cis_only = T) {
  pqtl_file_pattern <- if(cis_only) ".*_cis_only.*" else ".*_cis_trans.*"
  pqtl_datasets <- list.files(qtl_dir, pattern = pqtl_file_pattern, full.names = T)

  run_mr_on_qtl_data(gwas_filename, qtl_files = pqtl_datasets, results_output = results_output)
}

#' @import dplyr
#' @import vroom
#' @import TwoSampleMR
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

    qtl_dataset <- vroom::vroom(qtl_file, show_col_types=F) |>
      calculate_f_statistic()
    matching <- match(mr_results$exposure, qtl_dataset$EXPOSURE)
    mr_results$SNP <- qtl_dataset$SNP[matching]
    mr_results$F_STAT <- qtl_dataset$F_STAT[matching]

    return(mr_results)
  }) |> dplyr::bind_rows()

  all_qtl_mr_results <- tidyr::separate(all_qtl_mr_results, col = "SNP", into = c("CHR", "BP", "EA", "OA"), sep = "[:_]", remove = F) |>
    dplyr::rename(BETA="b", SE="se", P="pval", EXPOSURE="exposure") |>
    dplyr::mutate(p.adjusted = p.adjust(P, "fdr"))

  vroom::vroom_write(all_qtl_mr_results, results_output)
}

#'
#'
#' @import dplyr
#' @import vroom
compare_interesting_mr_results <- function(pqtl_mr_results, forest_plot_output_file) {
  pqtl_mr_results <- vroom::vroom(pqtl_mr_results, show_col_types=F)
  significant_pqtl_mr_results <- subset(pqtl_mr_results, p.adjusted < 0.05)

  interesting_pqtl_mr_results <- subset(pqtl_mr_results, exposure %in% significant_pqtl_mr_results$exposure) |>
    dplyr::select(exposure, dplyr::everything())
  interesting_pqtl_mr_results$BETA <- interesting_pqtl_mr_results$b
  interesting_pqtl_mr_results$SE <- interesting_pqtl_mr_results$se

  grouped_forest_plot(interesting_pqtl_mr_results,
                      "Comparing MR Results across pQTL data sets",
                      "outcome",
                      forest_plot_output_file
  )
}
