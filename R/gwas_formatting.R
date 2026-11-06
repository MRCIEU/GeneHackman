#' standardise_gwas: takes an input gwas, changes headers, standardises allelic input, adds RSID, makes life easier
#' @param gwas filename of gwas to standardise or dataframe of gwas
#' @param output_file file to save standardised gwas
#' @param N sample size of GWAS: fills the N column where it is missing; a per-variant N column in the file
#'   takes precedence
#' @param populate_rsid_option if you want RSID populated or not
#' @param input_reference_build reference build of CHR and BP of data.  Defaults to GRCh37
#' @param output_reference_build reference build of CHR and BP of the output.  Defaults to GRCh38
#' @param input_columns column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @param output_columns column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @param remove_extra_columns logical flag on if to remove extra columns
#' @param populate_eaf If TRUE, fill missing \code{EAF} from the LD reference: the 1000 Genomes
#'   \code{.frq} file via RSID matching (a partial RSID population is performed automatically when
#'   RSIDs are not already present), or the Pan-UKB LD pack via SNP keys.  Requires \code{ancestry}.
#' @param ancestry Ancestry code (EUR, EAS, AFR, AMR, SAS) matching the reference panel prefix; required when \code{populate_eaf} is TRUE.
#' @param flip_alleles If TRUE (default), keep EA/OA in alphabetical order with flipped
#'   BETA/EAF/Z in the output.  When FALSE the original allele order and effect direction
#'   are restored before saving; all internal steps still operate on the canonical
#'   (alphabetically flipped) representation.
#' @param ld_reference reference for \code{populate_eaf}: "1000genomes" or "pan_ukb"
#' @param complement_minus_strand liftover: keep non-palindromic variants that map to the minus strand,
#'   with complemented alleles (by default they are dropped)
#' @param max_liftover_drop_fraction liftover: stop if more than this fraction of variants is dropped
#' @param study_type "continuous" (default) or "binary" (case-control); for binary traits an effective
#'   sample size column N_EFF is added, see \code{populate_effective_n()}
#' @param n_case,n_control number of cases and controls of a binary trait (whole study), used for N_EFF
#'   when the file has no per-variant N_CASES and N_CONTROLS
#' @param check_build if TRUE, check the positions against the Pan-UKB LD pack's GRCh37 and GRCh38 keys
#'   and stop when they match the other build much better than \code{input_reference_build}
#' @return modified gwas: saves new gwas in output_file if present; a liftover QC report is written
#'   beside it (\verb{<output>.liftover_qc.tsv}, \verb{<output>.liftover_dropped.tsv.gz}) when the build changes


#' @export
standardise_gwas <- function(gwas,
                             output_file,
                             N=0,
                             populate_rsid_option=populate_rsid_options$none,
                             input_reference_build=reference_builds$GRCh37,
                             output_reference_build=reference_builds$GRCh38,
                             input_columns="default",
                             output_columns="default",
                             remove_extra_columns=F,
                             populate_eaf=FALSE,
                             ancestry=NULL,
                             flip_alleles=TRUE,
                             ld_reference=ld_reference_options$thousand_genomes,
                             complement_minus_strand=FALSE,
                             max_liftover_drop_fraction=0.05,
                             check_build=TRUE,
                             study_type=study_types$continuous,
                             n_case=NA,
                             n_control=NA) {
  input_gwas_columns <- resolve_column_map(input_columns)
  output_gwas_columns <- resolve_column_map(output_columns)

  if (populate_eaf && (is.null(ancestry) || length(ancestry) == 0 || !nzchar(as.character(ancestry)[1]))) {
    stop("populate_eaf is TRUE but ancestry is missing; set ancestry to EUR, EAS, AFR, AMR, or SAS to match your LD reference panel.")
  }

  gwas <- get_file_or_dataframe(gwas) |>
    change_column_names(input_gwas_columns, remove_extra_columns) |>
    standardise_columns(N) |>
    filter_incomplete_rows() |>
    populate_effective_n(study_type, n_case, n_control)

  if (isTRUE(check_build)) {
    check_declared_build(gwas, input_reference_build)
  }
  report_prefix <- NULL
  if (!missing(output_file) && shiny::isTruthy(output_file)) {
    report_prefix <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", output_file)
  }
  gwas <- convert_reference_build_via_liftover(gwas, input_reference_build, output_reference_build,
                                               complement_minus_strand = complement_minus_strand,
                                               max_drop_fraction = max_liftover_drop_fraction,
                                               report_prefix = report_prefix)

  gwas <- standardise_alleles(gwas)

  rsid_option <- populate_rsid_option
  if (isTRUE(populate_eaf) && ld_reference == ld_reference_options$thousand_genomes &&
      identical(rsid_option, populate_rsid_options$none) && !"RSID" %in% colnames(gwas)) {
    rsid_option <- populate_rsid_options$partial
  }

  gwas <- gwas |>
    health_check() |>
    populate_rsid(rsid_option, build = output_reference_build) |>
    populate_gene_names()

  if (isTRUE(populate_eaf)) {
    gwas <- populate_eaf_from_reference_panel(gwas, ancestry, ld_reference, output_reference_build)
  }

  if (!flip_alleles) {
    gwas <- unflip_alleles(gwas)
  } else {
    gwas$.ALLELES_FLIPPED <- NULL
  }

  gwas <- change_column_names(gwas, output_gwas_columns)

  if (!missing(output_file) && shiny::isTruthy(output_file)) {
    vroom::vroom_write(gwas, output_file)
  }
  return(gwas)
}

#' harmonise_gwases: takes a list of gwases, get the SNPs in common
#' across all datasets arranged to be in the same order
#'
#' @param ... elipses of gwases
#' @return list of harmonised gwases

#' @export
harmonise_gwases <- function(...) {
  gwases <- list(...)

  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
  message(paste("Number of shared SNPs after harmonisation:", length(snpids)))

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, SNP %in% snpids & !duplicated(SNP)) |>
      dplyr::arrange(SNP)
  })

  return(gwases)
}

filter_incomplete_rows <- function(gwas) {
  filtered_gwas <- gwas[!is.na(gwas$OA) & !is.null(gwas$OA) &
                        !is.na(gwas$EA) & !is.null(gwas$EA) &
                        !is.na(gwas$CHR) & !is.null(gwas$CHR) &
                        !is.na(gwas$BP) & !is.null(gwas$BP),
  ]

  filtered_rows <- nrow(gwas) - nrow(filtered_gwas)
  if (nrow(filtered_gwas) == 0) {
    stop("Error: all rows have been filtered from GWAS due to lack of information.  Stopping")
  } else if (filtered_rows > 0) {
    warning(paste("Warning: Filtering out ", filtered_rows, "rows due to NULLs and NAs"))
  }
  return(filtered_gwas)
}

#' Normalise CHR for GWAS rows (chr-prefix, whitespace, leading zeros) so downstream joins match PLINK panels.
#' @keywords internal
normalise_gwas_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("^0+", "", x)
  x
}

standardise_columns <- function(gwas, N) {
  gwas_columns <- colnames(gwas)

  gwas <- fill_sample_size(gwas, N)

  if (!all(c("CHR", "BP") %in% gwas_columns) && "SNP" %in% gwas_columns) {
    gwas <- populate_chr_bp_from_variant_id(gwas)
  }

  if (all(c("OR", "OR_LB", "OR_UB") %in% gwas_columns) && !all(c("BETA", "SE") %in% colnames(gwas))) {
    gwas <- convert_or_to_beta(gwas)
  }

  if ("LOG_P" %in% gwas_columns && ! "P" %in% gwas_columns) {
    gwas <- convert_negative_log_p_to_p(gwas)
  }

  if ("Z" %in% gwas_columns && !"BETA" %in% gwas_columns) {
    gwas <- convert_z_score_to_beta(gwas)
  }

  if ("BP" %in% colnames(gwas)) gwas$BP <- as.numeric(gwas$BP)
  if ("P" %in% colnames(gwas)) {
    gwas$P <- as.numeric(gwas$P)
    n_zero <- sum(gwas$P == 0, na.rm = TRUE)
    if (n_zero > 0) {
      warning(n_zero, " P values of 0 set to ", .Machine$double.xmin,
              " (the smallest positive double); fine-mapping and coloc use BETA and SE, not P")
      gwas$P[!is.na(gwas$P) & gwas$P == 0] <- .Machine$double.xmin
    }
  }
  if ("BETA" %in% colnames(gwas)) {
    gwas$BETA <- as.numeric(gwas$BETA)
  }

  if ("CHR" %in% colnames(gwas)) {
    gwas$CHR <- normalise_gwas_chr(gwas$CHR)
  }

  return(gwas)
}

#' Fill the N column from the sample size given in the input YAML.
#'
#' A per-variant N in the file takes precedence: the YAML value only fills rows where N is missing (or the
#' whole column when the file has none).  In meta-analyses N differs between variants, and a single value
#' would misstate it for variants not measured in every cohort.
#' @keywords internal
fill_sample_size <- function(gwas, N) {
  N <- suppressWarnings(as.numeric(N))
  if (length(N) != 1 || is.na(N) || N <= 0) return(gwas)
  if (!"N" %in% colnames(gwas)) {
    gwas$N <- N
    return(gwas)
  }
  per_variant <- suppressWarnings(as.numeric(gwas$N))
  missing_n <- is.na(per_variant) | per_variant <= 0
  if (any(!missing_n)) {
    message("Using the GWAS file's per-variant N; the N given in the input (", N, ") fills ",
            sum(missing_n), " variant(s) without one")
  }
  per_variant[missing_n] <- N
  gwas$N <- per_variant
  return(gwas)
}

#' Effective sample size of a case-control study: \eqn{4 / (1/N_{cases} + 1/N_{controls})} (Willer et al. 2010).
#' @param n_case number of cases
#' @param n_control number of controls
#' @return effective sample size; NA where either count is missing or not positive
#' @export
effective_sample_size <- function(n_case, n_control) {
  n_case <- suppressWarnings(as.numeric(n_case))
  n_control <- suppressWarnings(as.numeric(n_control))
  n_eff <- 4 / (1 / n_case + 1 / n_control)
  n_eff[is.na(n_case) | is.na(n_control) | n_case <= 0 | n_control <= 0] <- NA_real_
  n_eff
}

#' Normalise a study type: "continuous" (also "quantitative", "quant") or "binary" (also "case_control",
#' "categorical", "cc").
#' @keywords internal
normalise_study_type <- function(study_type) {
  type <- tolower(trimws(as.character(study_type %||% study_types$continuous)))
  if (length(type) == 0 || is.na(type) || !nzchar(type)) return(study_types$continuous)
  if (type %in% c("continuous", "quantitative", "quant")) return(study_types$continuous)
  if (type %in% c("binary", "case_control", "case-control", "categorical", "cc")) return(study_types$binary)
  stop("study_type must be continuous or binary (case-control), not ", study_type)
}

#' Add an effective sample size column (N_EFF) for a binary trait.
#'
#' In order of preference:
#' \enumerate{
#'   \item an N_EFF column already in the file is kept (for a meta-analysis, the sum of the cohorts'
#'     effective sample sizes; Grotzinger et al. 2023);
#'   \item per-variant N_CASES and N_CONTROLS columns: N_EFF = 4 / (1/N_CASES + 1/N_CONTROLS);
#'   \item a per-variant N column and the study's case and control counts: N_EFF = N x 4 phi (1 - phi), with
#'     phi = n_case / (n_case + n_control), which assumes the same case fraction at every variant;
#'   \item the study's case and control counts alone: one N_EFF for every variant.
#' }
#' Continuous traits are returned unchanged.  A binary trait with none of these keeps only N, with a warning.
#' @param gwas GWAS data frame
#' @param study_type "continuous" or "binary"
#' @param n_case,n_control numbers of cases and controls in the whole study
#' @return gwas, with N_EFF for binary traits
#' @export
populate_effective_n <- function(gwas, study_type = study_types$continuous, n_case = NA, n_control = NA) {
  if (normalise_study_type(study_type) != study_types$binary) return(gwas)
  columns <- colnames(gwas)
  if ("N_EFF" %in% columns) {
    gwas$N_EFF <- suppressWarnings(as.numeric(gwas$N_EFF))
    return(gwas)
  }
  if (all(c("N_CASES", "N_CONTROLS") %in% columns)) {
    gwas$N_EFF <- effective_sample_size(gwas$N_CASES, gwas$N_CONTROLS)
    if (!all(is.na(gwas$N_EFF))) return(gwas)
  }
  study_n_eff <- effective_sample_size(n_case, n_control)
  if (is.na(study_n_eff)) {
    warning("Binary trait without case and control counts (N_case and N_control in the input, or N_CASES and ",
            "N_CONTROLS columns): N_EFF can't be computed, so the total N is used, which overstates the ",
            "information in an unbalanced case-control study")
    return(gwas)
  }
  if ("N" %in% columns && any(!is.na(suppressWarnings(as.numeric(gwas$N))))) {
    case_fraction <- as.numeric(n_case) / (as.numeric(n_case) + as.numeric(n_control))
    gwas$N_EFF <- suppressWarnings(as.numeric(gwas$N)) * 4 * case_fraction * (1 - case_fraction)
    gwas$N_EFF[is.na(gwas$N_EFF)] <- study_n_eff
  } else {
    gwas$N_EFF <- study_n_eff
  }
  return(gwas)
}

#' Sample size to use for a set of variants: the median of N_EFF (binary traits) or N, ignoring missing values.
#' @param gwas GWAS rows (e.g. one fine-mapping region)
#' @param default_n used when the rows have no sample size
#' @keywords internal
analysis_sample_size <- function(gwas, default_n = NA_real_) {
  for (column in c("N_EFF", "N")) {
    if (column %in% colnames(gwas)) {
      n <- suppressWarnings(as.numeric(gwas[[column]]))
      n <- n[!is.na(n) & n > 0]
      if (length(n) > 0) return(stats::median(n))
    }
  }
  return(suppressWarnings(as.numeric(default_n)))
}

#' Fill missing CHR / BP columns from \code{chr:pos[:a1:a2]} style variant identifiers in the SNP column.
#' Accepts an optional \code{chr} prefix, \code{:}, \code{_} or \code{-} separators,
#' numeric chromosomes and X, Y, XY, M or MT.
#' Rows whose identifier does not match (e.g. rsIDs) are left with missing CHR / BP and are removed later by
#' \code{filter_incomplete_rows()}. Alleles embedded in the identifier are not used: EA and OA must come from
#' their own columns, because the identifier does not say which allele is the effect allele.
#' @keywords internal
populate_chr_bp_from_variant_id <- function(gwas) {
  variant_id_pattern <- stringr::regex("^(?:chr)?([0-9]{1,2}|X|Y|XY|MT|M)[:_-]([0-9]+)(?:[:_-].*)?$",
                                       ignore_case = TRUE)
  parsed <- stringr::str_match(stringr::str_trim(as.character(gwas$SNP)), variant_id_pattern)
  matched <- !is.na(parsed[, 1])

  if (!any(matched)) {
    return(gwas)
  }
  if (!all(matched)) {
    warning(paste("Could not parse CHR and BP from", sum(!matched), "SNP identifiers; these rows will be removed"))
  }

  if (!"CHR" %in% colnames(gwas)) gwas$CHR <- toupper(parsed[, 2])
  if (!"BP" %in% colnames(gwas)) gwas$BP <- as.numeric(parsed[, 3])
  return(gwas)
}

health_check <- function(gwas) {
  gwas_columns <- colnames(gwas)
  if ("P" %in% gwas_columns && (nrow(gwas[gwas$P <= 0 | gwas$P > 1, ]) > 0)) {
    warning("GWAS has some P values outside accepted range")
  }
  if ("EAF" %in% gwas_columns && (nrow(gwas[gwas$EAF < 0 | gwas$EAF > 1, ]) > 0)) {
    warning("GWAS has some EAF values outside accepted range")
  }
  return(gwas)
}

#' change_column_names: private function that takes a named list of column names
#'  and changes the supplied data frame's column names accordingly
#'
#' @param gwas dataframe of gwas to standardise column names
#' @param columns named list for
#' @param remove_extra_columns logical flag on if we are removing extra columns
#' @return gwas with new column names
change_column_names <- function(gwas, columns = list(), remove_extra_columns = F) {
  for (name in names(columns)) {
    #this deletes an existing column that we're about to rename, so we don't have 2 columns
    already <- name != columns[[name]]
    already <- already %in% TRUE
    if (name %in% names(gwas) & already) {
      gwas <- gwas[ , -which(names(gwas) %in% c(name))]
    }
    names(gwas)[names(gwas) == columns[name]] <- name
  }

  if (remove_extra_columns) {
    columns_to_remove <- setdiff(colnames(gwas), names(columns))
    gwas <- gwas[,-which(colnames(gwas) %in% columns_to_remove)]
  }

  return(gwas)
}

#' Reorder EA/OA alphabetically and flip BETA, EAF, Z to match.
#' Rows that were flipped are marked in a logical \code{.ALLELES_FLIPPED} column
#' so they can be restored later by \code{unflip_alleles()}.
#' @keywords internal
standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  gwas$.ALLELES_FLIPPED <- to_flip

  if (any(to_flip)) {
    if ("EAF" %in% names(gwas)) gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
    if ("Z" %in% names(gwas)) gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  gwas <- dplyr::select(gwas, SNP, CHR, BP, EA, OA, dplyr::everything())

  return(gwas)
}

#' Restore original allele order for rows that were flipped by \code{standardise_alleles()}.
#' Reverses BETA, EAF, Z and swaps EA/OA back, then rebuilds the SNP column.
#' @keywords internal
unflip_alleles <- function(gwas) {
  if (!".ALLELES_FLIPPED" %in% names(gwas)) return(gwas)

  to_unflip <- gwas$.ALLELES_FLIPPED
  if (any(to_unflip)) {
    if ("EAF" %in% names(gwas)) gwas$EAF[to_unflip] <- 1 - gwas$EAF[to_unflip]
    gwas$BETA[to_unflip] <- -1 * gwas$BETA[to_unflip]
    if ("Z" %in% names(gwas)) gwas$Z[to_unflip] <- -1 * gwas$Z[to_unflip]

    temp <- gwas$OA[to_unflip]
    gwas$OA[to_unflip] <- gwas$EA[to_unflip]
    gwas$EA[to_unflip] <- temp
  }

  gwas$.ALLELES_FLIPPED <- NULL
  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  return(gwas)
}

#' convert_or_to_beta: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param gwas dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return gwas with new columns BETA and SE

#' @export
convert_or_to_beta <- function(gwas) {
  gwas <- get_file_or_dataframe(gwas)
  if (!all(c("OR", "OR_LB", "OR_UB") %in% colnames(gwas))) {
    stop("Need OR, OR_LB + OR_UB to complete conversion")
  }

  z_score <- stats::qnorm(.975, mean = 0, sd = 1) #1.96
  gwas$BETA <- log(gwas$OR)
  gwas$SE <- (log(gwas$OR_LB) - gwas$BETA) / -z_score

  return(gwas)
}

#' convert_beta_to_or: Given a BETA and SE, calculates the OR and lower and upper bounds
#' @param gwas dataframe with the following columns: BETA, SE
#' @return gwas with new columns OR, OR_LB, OR_UB

#' @export
convert_beta_to_or <- function(gwas) {
  gwas <- get_file_or_dataframe(gwas)
  z_score <- stats::qnorm(.975, mean = 0, sd = 1) #1.96

  gwas$OR <- exp(gwas$BETA)
  gwas$OR_LB <- exp(gwas$BETA - z_score * gwas$SE)
  gwas$OR_UB <- exp(gwas$BETA + z_score * gwas$SE)
  return(gwas)
}

#' Effect size and standard error from a Z-score, allele frequency and sample size.
#'
#' For files that give only Z (no BETA and SE), the effect is put on the standardised scale (per standard
#' deviation of the trait): \eqn{SE = 1 / \sqrt{2p(1-p)(N + Z^2)}} and \eqn{BETA = Z \times SE}
#' (Zhu et al. 2016).  P is computed from Z when missing.  Needs EAF and N (a column, or \code{sample_size}).
#' @param gwas data frame with Z, EAF and N
#' @param sample_size sample size used where N is missing
#' @return gwas with BETA, SE and P; rows that can't be converted are removed with a warning
#' @keywords internal
convert_z_score_to_beta <- function(gwas, sample_size = NA) {
  if (!"EAF" %in% colnames(gwas)) {
    stop("The GWAS gives Z but no BETA and SE: converting Z needs an EAF column (and N)")
  }
  n <- if ("N" %in% colnames(gwas)) suppressWarnings(as.numeric(gwas$N)) else rep(NA_real_, nrow(gwas))
  n[is.na(n)] <- suppressWarnings(as.numeric(sample_size))
  if (all(is.na(n))) {
    stop("The GWAS gives Z but no BETA and SE: converting Z needs a sample size (an N column or N in the input)")
  }
  z <- as.numeric(gwas$Z)
  eaf <- as.numeric(gwas$EAF)
  gwas$SE <- 1 / sqrt(2 * eaf * (1 - eaf) * (n + z^2))
  gwas$BETA <- z * gwas$SE
  if (!"P" %in% colnames(gwas)) gwas$P <- 2 * stats::pnorm(-abs(z))
  gwas$P[is.na(gwas$P)] <- 2 * stats::pnorm(-abs(z[is.na(gwas$P)]))

  ok <- is.finite(gwas$BETA) & is.finite(gwas$SE) & gwas$SE > 0
  if (any(!ok)) {
    warning(sum(!ok), " variant(s) removed: Z can't be converted to BETA and SE (missing Z, EAF or N, or EAF of 0 or 1)")
  }
  return(gwas[ok, , drop = FALSE])
}

convert_z_to_p <- function(gwas) {
  gwas$P <- 2 * pnorm(-abs(gwas$Z))
  return(gwas)
}

# #TODO: don't think this is right...
# convert_beta_and_se_to_z_score <- function(gwas) {
#   gwas$Z <- gwas$BETA / gwas$SE
#   return(gwas)
# }

convert_negative_log_p_to_p <- function(gwas) {
  gwas$P <- 10^-gwas$LOG_P
  return(gwas)
}

convert_p_to_negative_log_p <- function(gwas) {
  gwas$LOG_P <- -log10(gwas$P)
  return(gwas)
}

#' calculate_f_statistic: calculates the F-statistic from the P-value
#' @param gwas dataframe with P column
#' @return gwas with F_STAT column

#' @export
calculate_f_statistic <- function(gwas) {
  gwas$F_STAT <- stats::qchisq(gwas$P, 1, low=F)
  return(gwas)
}


calculate_lambda_statistic <- function(gwas) {
  lambda <- median(stats::qchisq(1 - gwas$P, 1)) / stats::qchisq(0.5, 1)
  return(lambda)
}

resolve_column_map <- function(column_map) {
  column_map_file <- system.file("extdata", "predefined_column_maps.csv", package = "GeneHackman")
  if (!file.exists(column_map_file)) {
    predefined_column_maps <- vroom::vroom("../inst/extdata/predefined_column_maps.csv")
  } else {
    predefined_column_maps <- vroom::vroom(column_map_file, show_col_types = F)
  }
  predefined_column_maps <- tibble::column_to_rownames(predefined_column_maps, "name")

  if (is.vector(column_map) && length(column_map) > 1) {
    resolved_column_map <- column_map
  } else if (is.character(column_map) && length(column_map) == 1) {
    if (column_map %in% row.names(predefined_column_maps)) {
      predefined_map <- predefined_column_maps[column_map, ]
      resolved_column_map <- as.list(predefined_map)
    } else {
      split_map <- split_string_into_named_list(column_map)
      if (length(split_map) == 0) stop(paste("Error resolving column map for", column_map))
      resolved_column_map <- split_map
    }
  } else {
    stop(paste("Error resolving column map for", column_map))
  }

  return(resolved_column_map)
}
