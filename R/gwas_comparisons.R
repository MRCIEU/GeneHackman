#' compare_replication_across_all_gwas_permutations
#' @description takes a vector of gwasse and associated clumps, then runs the expected vs. observed algorithm
#'
#' @param gwas_filenames list of filenames of gwases
#' @param clumped_filenames list of clumped list filesnames genetred from gwas, in same order as gwases
#' @param result_output_file data frame of all results of gwas comparisons concatenated
#' @param variants_output_file data frame of every SNP comparison from clumped list concatenated

compare_replication_across_all_gwas_permutations <- function(gwas_filenames = c(),
                                                             clumped_filenames = c(),
                                                             result_output_file,
                                                             variants_output_file) {
  if (length(gwas_filenames) < 2) {
    empty_dataframe <- data.frame(SNP=c())
    vroom::vroom_write(empty_dataframe, result_output_file)
    vroom::vroom_write(empty_dataframe, variants_output_file)
    return()
  }

  expected_vs_observed_results <- list()
  expected_vs_observed_variants <- list()

  # every unordered pair once, with the earlier-listed GWAS (and its clumped hits) as the discovery dataset
  for (pair in utils::combn(seq_along(gwas_filenames), 2, simplify = FALSE)) {
    i <- pair[1]
    j <- pair[2]
    result <- compare_two_gwases_from_clumped_hits(gwas_filenames[i], gwas_filenames[j], clumped_filenames[i])
    expected_vs_observed_results <- c(expected_vs_observed_results, list(result$res))
    expected_vs_observed_variants <- c(expected_vs_observed_variants, list(result$variants))
  }
  merged_results <- do.call("rbind", expected_vs_observed_results)
  merged_variants <- do.call("rbind", expected_vs_observed_variants)

  vroom::vroom_write(merged_results, result_output_file)
  vroom::vroom_write(merged_variants, variants_output_file)
}



compare_two_gwases_from_clumped_hits <- function(first_gwas, second_gwas, clumped_snps) {
  comparison_name <- paste0(file_prefix(first_gwas), "_vs_", file_prefix(second_gwas))
  comparison_columns <- c("SNP", "BETA", "SE", "P", "RSID")

  first_gwas <- filter_gwas_by_clumped_results(first_gwas, clumped_snps)
  second_gwas <- get_file_or_dataframe(second_gwas, columns = comparison_columns)

  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  results <- expected_vs_observed_replication(harmonised_gwases[[1]]$BETA,
                                              harmonised_gwases[[1]]$SE,
                                              harmonised_gwases[[2]]$BETA,
                                              harmonised_gwases[[2]]$SE
  )
  results$res$COMPARISON <- comparison_name
  results$variants$COMPARISON <- comparison_name
  results$variants$SNP <- harmonised_gwases[[1]]$SNP
  results$variants <- dplyr::select(results$variants, COMPARISON, SNP, dplyr::everything())
  results$res <- dplyr::select(results$res, COMPARISON, dplyr::everything())

  return(results)
}

#' Expected vs observed replication rates
#'
#' @description For a set of effects that have discovery and replication betas and SEs, this function determines the
#'  extent to which the observed replication rate matches the expected replication rate.
#' The expected replication rate is based on the assumption that the replication dataset has the same effect sizes
#'  but that the power may be different (e.g. due to allele frequencies or sample sizes) and is reflected in the replication standard errors.
#' It assesses replication based on concordance of effect direction across discovery and replication, and p-values
#'  surpassing a user-specified p-value threshold.
#'
#' doi: 10.1038/nature17671
#'
#' @param b_disc vector of clumped incidence hit effects
#' @param se_disc the standard errors for incidence effects
#' @param b_rep corresponding vector of associations in progression
#' @param se_rep standard errors of effects in progression
#' @param alpha p-value threshold to check for replication of incidence hits in progression (e.g. try 0.05 or 1e-5)

expected_vs_observed_replication <- function(b_disc, se_disc, b_rep, se_rep, alpha = 0.05) {
  p_sign <- pnorm(-abs(b_disc) / se_disc) *
    pnorm(-abs(b_disc) / se_rep) + (
      (1 - pnorm(-abs(b_disc) / se_disc)) *
      (1 - pnorm(-abs(b_disc) / se_rep))
  )

  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) +
    (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))

  # two-sided, to match the two-sided threshold used for the expected count (p_sig)
  p_rep <- 2 * pnorm(abs(b_rep) / se_rep, lower.tail = FALSE)

  res <- tibble::tibble(
    N = length(b_disc),
    METRIC = c("Sign", "Sign", "P-value", "P-value"),
    DATUM = c("Expected", "Observed", "Expected", "Observed"),
    VALUE = c(
        sum(p_sign, na.rm=TRUE),
        sum(sign(b_disc) == sign(b_rep)),
        sum(p_sig, na.rm=TRUE),
        sum(p_rep < alpha, na.rm=TRUE)
    ),
    P_DIFF = c(
      NA_real_,
            binom.test(VALUE[2], N[2], VALUE[1]/N[2])$p.value,
            NA_real_,
            binom.test(VALUE[4], N[4], VALUE[3]/N[4])$p.value
    )
  )

  res_per_variant <- tibble::tibble(
    EXPECTED = p_sig,
    OBSERVED = p_rep < alpha,
    REPLICATION_FAIL = EXPECTED > 0.95 & !OBSERVED,
    EXPECTED_SIGN = p_sign,
    OBSERVED_SIGN = sign(b_disc) == sign(b_rep),
    SIGN_FAIL = EXPECTED_SIGN > 0.95 & !OBSERVED_SIGN
  )

  return(list(res = res, variants = res_per_variant))
}





dataset_label_from_gwas_path <- function(path) {
  x <- basename(path)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  tools::file_path_sans_ext(x)
}


compare_heterogeneity_across_ancestries <- function(gwas_filenames = c(),
                                                    clumped_filenames = c(),
                                                    ancestry_list = c(),
                                                    heterogeneity_score_file,
                                                    heterogeneity_plot_file,
                                                    heterogeneity_plots_per_snp_file) {

  if (length(gwas_filenames) < 2) {
    empty_dataframe <- data.frame(SNP=c())
    vroom::vroom_write(empty_dataframe, heterogeneity_score_file)
    file.create(heterogeneity_plot_file)
    file.create(heterogeneity_plots_per_snp_file)
    return()
  }

  comparison_columns <- c("SNP", "BETA", "SE", "RSID")
  clumped_snps <- data.table::rbindlist(lapply(clumped_filenames, data.table::fread))$SNP
  clumped_snps <- unique(clumped_snps)

  dataset_labels <- make.unique(
    vapply(gwas_filenames, dataset_label_from_gwas_path, character(1)),
    sep = "_"
  )

  gwases <- lapply(gwas_filenames, function(gwas_filename) {
    get_file_or_dataframe(gwas_filename, columns = comparison_columns) |>
      dplyr::filter(RSID %in% clumped_snps)
  })
  for (i in seq_along(gwases)) {
    gwases[[i]]$dataset <- dataset_labels[[i]]
  }

  harmonised_gwases <- rlang::inject(harmonise_gwases(!!!gwases))
  gwases_by_dataset <- stats::setNames(harmonised_gwases, dataset_labels)

  heterogeneity_results <- calculate_heterogeneity_scores(gwases_by_dataset, heterogeneity_score_file)
  plot_heritability_contribution_per_ancestry(heterogeneity_results$Qj, heterogeneity_plot_file)
  plot_snps_with_heterogeneity(gwases_by_dataset, heterogeneity_results$Q, heterogeneity_plots_per_snp_file)
}


plot_snps_with_heterogeneity <- function(gwases_by_dataset, heterogeneity_resuts_q, heterogeneity_forest_plot_file) {
  are_heterogeneous <- which(heterogeneity_resuts_q$Qpval < 0.05 / nrow(heterogeneity_resuts_q))

  forest_plot_rows <- data.table::rbindlist(lapply(gwases_by_dataset, function(gwas) {
    gwas$Qpval <- heterogeneity_resuts_q$Qpval
    gwas[are_heterogeneous, ]
  }))

  if ("dataset" %in% colnames(forest_plot_rows)) {
    data.table::setnames(forest_plot_rows, "dataset", "GWAS")
  }

  grouped_forest_plot(forest_plot_rows,
                      title = "Plot of SNPs That Are Heterogeneous",
                      group_column = "GWAS",
                      output_file = heterogeneity_forest_plot_file,
  )
}

#' Test for heterogeneity of effect estimates between GWAS datasets
#'
#' @description For each SNP this function will provide a Cochran's Q test statistic -
#'  a measure of heterogeneity of effect sizes between datasets. A low p-value means high heterogeneity.
#'  In addition, for every SNP it gives a per-dataset p-value —
#'  interpreted as whether a particular study is an outlier at that SNP.
#'
#' @param gwases_by_ancestry Named list of data frames, one per GWAS dataset, with at least BETA, SE and SNP columns
#' @param heterogeneity_score_file file to save the heterogeneity scores
#'
#' @return List
#' - Q = vector of p-values for Cochrane's Q statistic for each SNP
#' - Qj = Data frame of per-population outlier q values for each SNP




calculate_heterogeneity_scores <- function(gwases_by_ancestry, heterogeneity_score_file) {
  beta <- lapply(gwases_by_ancestry, \(x) x$BETA) |> dplyr::bind_cols()
  se <- lapply(gwases_by_ancestry, \(x) x$SE) |> dplyr::bind_cols()
  o <- lapply(seq_len(nrow(beta)), \(i) {
    fixed_effects_meta_analysis(as.numeric(beta[i,]), as.numeric(se[i,]))
  })

  Q <- tibble::tibble(SNP = gwases_by_ancestry[[1]]$SNP,
                      Qpval = sapply(o, \(x) x$Qpval)
  )

  Qj <- lapply(o, \(x) x$Qjpval)
  Qj <-  do.call(rbind, Qj) |>
    tibble::as_tibble() |>
    dplyr::rename(setNames(paste0("V", seq_along(gwases_by_ancestry)), names(gwases_by_ancestry))) |>
    dplyr::mutate(SNP = gwases_by_ancestry[[1]]$SNP)

  joined_heteogeneity_data <- merge(Q, Qj, by="SNP")

  vroom::vroom_write(joined_heteogeneity_data, heterogeneity_score_file)
  return(list(Q=Q, Qj=Qj))
}

fixed_effects_meta_analysis <- function(betas, standard_errors) {
  w <- 1 / standard_errors^2
  beta <- sum(betas * w) / sum(w)
  se <- sqrt(1 / sum(w))

  Qj <- w * (beta-betas)^2
  Q <- sum(Qj)
  Qdf <- length(betas)-1
  Qjpval <- stats::pchisq(Qj, 1, lower.tail=FALSE)
  Qpval <- stats::pchisq(Q, Qdf, lower.tail=FALSE)

  return(list(beta=beta, se=se, Qpval=Qpval, Qj=Qj, Qjpval=Qjpval))
}


#' Parse an LDSC --h2 log.
#'
#' @param log_file LDSC log file
#' @return one-row data.frame: H2_OBS, H2_OBS_SE, H2_LIAB, H2_LIAB_SE (NA unless prevalences were given),
#'   LAMBDA_GC, MEAN_CHISQ, INTERCEPT, INTERCEPT_SE, RATIO, RATIO_SE and RATIO_NOTE (e.g. LDSC's
#'   "Ratio < 0 (usually indicates GC correction)")
#' @export
parse_ldsc_h2_log <- function(log_file) {
  log_lines <- readLines(log_file, warn = FALSE)
  value_and_se <- function(label) {
    line <- grep(paste0("^", label, ":"), log_lines, value = TRUE)
    if (length(line) == 0) return(c(NA_real_, NA_real_))
    after_label <- sub("^[^:]*:", "", line[1])
    numbers <- regmatches(after_label, gregexpr("-?[0-9.]+(e-?[0-9]+)?|nan|NA", after_label, ignore.case = TRUE))[[1]]
    as.numeric(suppressWarnings(c(numbers[1], if (length(numbers) > 1) numbers[2] else NA)))
  }
  observed <- value_and_se("Total Observed scale h2")
  liability <- value_and_se("Total Liability scale h2")
  intercept <- value_and_se("Intercept")
  ratio <- value_and_se("Ratio")
  ratio_line <- grep("^Ratio", log_lines, value = TRUE)
  data.frame(
    H2_OBS = observed[1], H2_OBS_SE = observed[2],
    H2_LIAB = liability[1], H2_LIAB_SE = liability[2],
    LAMBDA_GC = value_and_se("Lambda GC")[1], MEAN_CHISQ = value_and_se("Mean Chi\\^2")[1],
    INTERCEPT = intercept[1], INTERCEPT_SE = intercept[2],
    RATIO = ratio[1], RATIO_SE = ratio[2],
    RATIO_NOTE = if (length(ratio_line) > 0 && !grepl("^Ratio:", ratio_line[1])) trimws(ratio_line[1]) else NA_character_
  )
}

#' Parse the "Summary of Genetic Correlation Results" table of an LDSC --rg log.
#' @param log_file LDSC log file
#' @return data.frame with one row per pair (p1, p2, rg, se, z, p and the h2 / intercept columns LDSC prints)
#' @export
parse_ldsc_rg_log <- function(log_file) {
  log_lines <- readLines(log_file, warn = FALSE)
  start <- grep("^Summary of Genetic Correlation Results", log_lines)
  if (length(start) == 0) return(data.frame())
  table_lines <- log_lines[(start[1] + 1):length(log_lines)]
  end <- which(!nzchar(trimws(table_lines)))
  if (length(end) > 0) table_lines <- table_lines[seq_len(end[1] - 1)]
  if (length(table_lines) < 2) return(data.frame())
  utils::read.table(text = table_lines, header = TRUE, stringsAsFactors = FALSE, na.strings = c("NA", "nan"))
}

#' Tables of LDSC heritability and genetic correlation results, with reliability flags.
#'
#' Heritability: the z-score of the observed-scale estimate (H2_Z) and a flag when it is below 4, the threshold
#' below which genetic correlations are unreliable (Zheng et al. 2017, LD Hub).  Genetic correlation: a flag
#' when |rg| > 1 or when either trait's H2_Z is below 4; RG_RELIABLE is FALSE in either case, and such
#' estimates should not be interpreted.  For binary traits with prevalences, H2_LIAB is the liability-scale
#' heritability (with effective N and a sample prevalence of 0.5; Grotzinger et al. 2023).
#'
#' @param h2_logs LDSC --h2 logs, one per GWAS in the order of \code{sumstats}
#' @param rg_logs LDSC --rg logs
#' @param sumstats munged sumstats file of each GWAS (to name the traits in the --rg tables)
#' @param labels trait labels, in the order of \code{sumstats}
#' @param ancestry ancestry of the LD scores
#' @param population_prevalences population prevalence used for each trait (NA or "nan" for none)
#' @return list with \code{h2} and \code{rg} data.frames
#' @export
ldsc_results_tables <- function(h2_logs, rg_logs, sumstats, labels = basename(sumstats), ancestry = NA_character_,
                                population_prevalences = rep(NA, length(sumstats))) {
  population_prevalences <- suppressWarnings(as.numeric(population_prevalences))
  h2 <- do.call(rbind, lapply(h2_logs, parse_ldsc_h2_log))
  h2 <- cbind(data.frame(ANCESTRY = ancestry, TRAIT = labels, POP_PREV = population_prevalences), h2)
  h2$H2_Z <- h2$H2_OBS / h2$H2_OBS_SE
  h2$FLAG_LOW_H2_Z <- is.na(h2$H2_Z) | h2$H2_Z < 4
  h2$FLAG_RATIO <- !is.na(h2$RATIO_NOTE)

  rg <- do.call(rbind, lapply(rg_logs, parse_ldsc_rg_log))
  if (is.null(rg) || nrow(rg) == 0) {
    return(list(h2 = h2, rg = data.frame()))
  }
  trait_of <- function(path) labels[match(basename(path), basename(sumstats))]
  h2_z <- stats::setNames(h2$H2_Z, labels)
  rg <- data.frame(ANCESTRY = ancestry, TRAIT_1 = trait_of(rg$p1), TRAIT_2 = trait_of(rg$p2),
                   RG = rg$rg, RG_SE = rg$se, RG_Z = rg$z, RG_P = rg$p,
                   GCOV_INTERCEPT = rg$gcov_int %||% NA_real_, GCOV_INTERCEPT_SE = rg$gcov_int_se %||% NA_real_)
  rg$H2_Z_1 <- unname(h2_z[rg$TRAIT_1])
  rg$H2_Z_2 <- unname(h2_z[rg$TRAIT_2])
  rg$FLAG_RG_OUT_OF_BOUNDS <- !is.na(rg$RG) & abs(rg$RG) > 1
  rg$FLAG_LOW_H2_Z <- is.na(rg$H2_Z_1) | is.na(rg$H2_Z_2) | rg$H2_Z_1 < 4 | rg$H2_Z_2 < 4
  rg$RG_RELIABLE <- !is.na(rg$RG) & !rg$FLAG_RG_OUT_OF_BOUNDS & !rg$FLAG_LOW_H2_Z
  list(h2 = h2, rg = rg)
}
