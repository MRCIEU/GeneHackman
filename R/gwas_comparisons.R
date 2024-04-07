#' compare_replication_across_all_gwas_permutations
#' @description takes a vector of gwasse and associated clumps, then runs the expected vs. observed algorithm
#'
#' @param gwas_filenames: list of filenames of gwases
#' @param clumped_filenames: list of clumped list filesnames genetred from gwas, in same order as gwases
#' @param result_output_file: data frame of all results of gwas comparisons concatenated
#' @param variants_output_file: data frame of every SNP comparison from clumped list concatenated
#' @import vroom
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

  for (i in seq_along(gwas_filenames)) {
    for (j in i:length(gwas_filenames)) {
      j <- j + 1
      if (j <= length(gwas_filenames)) {
        result <- compare_two_gwases_from_clumped_hits(gwas_filenames[i], gwas_filenames[j], clumped_filenames[i])
        expected_vs_observed_results[[i+j]] <- result$res
        expected_vs_observed_variants[[i+j]] <- result$variants
      }
    }
  }
  merged_results <- do.call("rbind", expected_vs_observed_results)
  merged_variants <- do.call("rbind", expected_vs_observed_variants)

  vroom::vroom_write(merged_results, result_output_file)
  vroom::vroom_write(merged_variants, variants_output_file)
}

#' @import dplyr
#' @import data.table
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
#' @import tibble
expected_vs_observed_replication <- function(b_disc, se_disc, b_rep, se_rep, alpha = 0.05) {
  p_sign <- pnorm(-abs(b_disc) / se_disc) *
    pnorm(-abs(b_disc) / se_rep) + (
      (1 - pnorm(-abs(b_disc) / se_disc)) *
      (1 - pnorm(-abs(b_disc) / se_rep))
  )

  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) +
    (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))

  p_rep <- pnorm(abs(b_rep) / se_rep, lower.tail = FALSE)

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

#' @import data.table
#' @import dplyr
#' @import rlang
#' @import stats
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

  gwases <- lapply(gwas_filenames, function(gwas_filename) {
    get_file_or_dataframe(gwas_filename, columns = comparison_columns) |>
      dplyr::filter(RSID %in% clumped_snps)
  })
  gwas_ancestry_names <- c()
  for (i in seq_along(gwases)) {
    gwases[[i]]$ancestry <- ancestry_list[[i]]
    gwas_ancestry_names[[i]] <- paste0(i, "_", ancestry_list[[i]])
  }

  harmonised_gwases <- rlang::inject(harmonise_gwases(!!!gwases))
  gwases_by_ancestry <- stats::setNames(harmonised_gwases, gwas_ancestry_names)

  heterogeneity_results <- calculate_heterogeneity_scores(gwases_by_ancestry, heterogeneity_score_file)
  plot_heritability_contribution_per_ancestry(heterogeneity_results$Qj, heterogeneity_plot_file)
  plot_snps_with_heterogeneity(gwases_by_ancestry, heterogeneity_results$Q, heterogeneity_plots_per_snp_file)
}

#' @import data.table
plot_snps_with_heterogeneity <- function(gwases_by_ancestry, heterogeneity_resuts_q, heterogeneity_forest_plot_file) {
  are_heterogeneous <- which(heterogeneity_resuts_q$Qpval < 0.05 / nrow(heterogeneity_resuts_q))

  forest_plot_rows <- data.table::rbindlist(lapply(gwases_by_ancestry, function(gwas) {
    gwas$Qpval <- heterogeneity_resuts_q$Qpval
    gwas[are_heterogeneous, ]
  }))

  grouped_forest_plot(forest_plot_rows,
                      title = "Plot of SNPs That Are Heterogeneous",
                      group_column = "ancestry",
                      output_file = heterogeneity_forest_plot_file,
                      #p_value_column = "Qpval" ??
  )
}

#' Test for heterogeneity of effect estimates between populations
#'
#' @description For each SNP this function will provide a Cochran's Q test statistic -
#'  a measure of heterogeneity of effect sizes between populations. A low p-value means high heterogeneity.
#'  In addition, for every SNP it gives a per population p-value -
#'  this can be interpreted as asking for each SNP is a particular giving an outlier estimate.
#'
#' @param gwases_by_ancestry Named list of data frames, one for each population, with at least bhat, se and snp columns
#'
#' @return List
#' - Q = vector of p-values for Cochrane's Q statistic for each SNP
#' - Qj = Data frame of per-population outlier q values for each SNP
#' @import tibble
#' @import dplyr
#' @import vroom
#' @import stats
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
