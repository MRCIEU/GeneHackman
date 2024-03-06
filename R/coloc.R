#'
#' @import dplyr
#' @import vroom
#' @import tibble
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

#' @import dplyr
#' @import tidyr
#' @import vroom
#' @import utils
run_coloc_on_qtl_mr_results <- function(mr_results_file,
                                        gwas_file,
                                        qtl_dataset,
                                        exposures=c(),
                                        default_n=NA,
                                        output_file) {
  coloc_columns <- c("SNP", "CHR", "BP", "P", "SE", "N", "EAF")
  gwas <- get_file_or_dataframe(gwas_file, columns = coloc_columns) |> dplyr::filter(EAF > 0 & EAF < 1)
  range <- 500000

  mr_results <- get_file_or_dataframe(mr_results_file) |>
    dplyr::filter(p.adjusted < 0.05)

  if (length(exposures) > 0) {
    mr_results <- subset(mr_results, exposure %in% exposures)
  }
  coloc_results <- apply(mr_results, 1, function(mr_result) {
    if (qtl_dataset == qtl_datasets$metabrain) {
      outcome <- unlist(strsplit(mr_result[["outcome"]], "_"))
      brain_region <- outcome[1]
      ancestry <- toupper(outcome[2])

      qtl_gwas_file <- paste0(metabrain_gwas_dir, "/", brain_region, "/", mr_result[["EXPOSURE"]], "_", ancestry, ".tsv.gz")
      qtl_gwas <- get_file_or_dataframe(qtl_gwas_file) |> dplyr::filter(EAF > 0 & EAF < 1)
    } else if (qtl_dataset == qtl_datasets$eqtlgen) {
      qtl_gwas_file <- paste0(eqtlgen_gwas_dir, "/", mr_result[["outcome"]], "/", mr_result[["EXPOSURE"]], ".tsv.gz")
      qtl_gwas <- get_file_or_dataframe(qtl_gwas_file) |> dplyr::filter(EAF > 0 & EAF < 1)
    } else {
      stop("Error: qtl dataset not supported for coloc right now")
    }
    chr <- as.numeric(mr_result[["CHR"]])
    bp <- as.numeric(mr_result[["BP"]])
    result <- coloc_analysis(gwas, qtl_gwas, mr_result[["EXPOSURE"]], chr, bp, range, default_n=default_n)

    #TODO: should we create Miami Plots for coloc results?
    #miami_filename <- paste0(user_results_dir, "plots/mr_metabrain_coloc_", file_prefix(gwas_file), "_", mr_result[["EXPOSURE"]], ".png")
    #miami_plot(qtl_gwas_file, gwas_file, miami_filename, paste0("Miami Plot of", mr_result[["EXPOSURE"]]), chr, bp, range)

    return(result)
  }) |> dplyr::bind_rows()

  vroom::vroom_write(coloc_results, output_file)
}

#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns tibble of coloc results (h0 - h4)
#' @import coloc
#' @import tibble
#' @export
coloc_analysis <- function(first_gwas, second_gwas, exposure_name, chr=NA, bp=NA, range=NA, default_n=NA) {
  if (!is.na(chr) & !is.na(bp) & !is.na(range)) {
    first_gwas <- gwas_region(first_gwas, chr, bp, range)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
  }
  numeric_columns <- c("P", "SE", "EAF")
  first_gwas[,numeric_columns] <- lapply(first_gwas[,numeric_columns,drop=FALSE], as.numeric)
  second_gwas[,numeric_columns] <- lapply(second_gwas[,numeric_columns,drop=FALSE], as.numeric)

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
