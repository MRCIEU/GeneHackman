#' standardise_gwas: takes an input gwas, changes headers, standardises allelic input, adds RSID, makes life easier
#' @param gwas: filename of gwas to standardise or dataframe of gwas
#' @param output_file: file to save standardised gwas
#' @param N: sample size of GWAS (if GWAS has defined N column, defaults to that value)
#' @param populate_rsid_option: if you want RSID populated or not
#' @param input_reference_build: reference build of CHR and BP of data.  Defaults to GRCh37
#' @param output_reference_build: reference build of CHR and BP of data.  Defaults to GRCh37
#' @param input_columns: column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @param output_columns: column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @return modified gwas: saves new gwas in {output_file} if present
#' @import vroom
#' @import shiny
#' @export
standardise_gwas <- function(gwas,
                             output_file,
                             N=0,
                             populate_rsid_option=F,
                             input_reference_build=reference_builds$GRCh37,
                             output_reference_build=reference_builds$GRCh37,
                             input_columns="default",
                             output_columns="default") {
  input_gwas_columns <- resolve_column_map(input_columns)
  output_gwas_columns <- resolve_column_map(output_columns)

  #TODO: if we need to add bespoke input format wrangling here, we can

  gwas <- get_file_or_dataframe(gwas) |>
    change_column_names(input_gwas_columns) |>
    standardise_columns(N) |>
    filter_incomplete_rows() |>
    convert_reference_build_via_liftover(input_reference_build, output_reference_build) |>
    standardise_alleles() |>
    health_check() |>
    populate_rsid(populate_rsid_option) |>
    populate_gene_names() |>
    change_column_names(output_gwas_columns)

  if (!missing(output_file) && shiny::isTruthy(output_file)) {
    vroom::vroom_write(gwas, output_file)
  }
  return(gwas)
}

#' harmonise_gwases: takes a list of gwases, get the SNPs in common
#' across all datasets arranged to be in the same order
#'
#' @param: elipses of gwases
#' @return: list of harmonised gwases
#' @import dplyr
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

#' @import tidyr
standardise_columns <- function(gwas, N) {
  gwas_columns <- colnames(gwas)

  if (!"N" %in% gwas_columns && N > 0) {
    gwas$N <- N
  }

  if (!all(c("CHR", "BP") %in% gwas_columns)) {
    if(all(grepl("\\d:\\d", gwas$SNP))) {
      gwas <- tidyr::separate(data = gwas, col = "SNP", into = c("CHR", "BP"), sep = "[:_]", remove = F)
      gwas$BP <- as.numeric(gwas$BP)
    }
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

  if ("BP" %in% gwas_columns) gwas$BP <- as.numeric(gwas$BP)
  if ("P" %in% gwas_columns) {
    gwas$P <- as.numeric(gwas$P)
    gwas$P[gwas$P == 0] <- .Machine$double.xmin
  }

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
#' @param: gwas dataframe of gwas to standardise column names
#' @param: columns named list for
#' @param: opposite_mapping logical flag on if we are mapping from key to value or vice verca
change_column_names <- function(gwas, columns = list()) {
  for (name in names(columns)) {
   names(gwas)[names(gwas) == columns[name]] <- name
  }

  return(gwas)
}

#' @import dplyr
standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
  if (any(to_flip)) {
    gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
    gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]

    temp <- gwas$OA[to_flip]
    gwas$OA[to_flip] <- gwas$EA[to_flip]
    gwas$EA[to_flip] <- temp
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  gwas <- dplyr::select(gwas, SNP, CHR, BP, EA, OA, dplyr::everything())

  return(gwas)
}

#' convert_or_to_beta: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param gwas: dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return gwas with new columns BETA and SE
#' @import stats
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

#' @import stats
#' @export
convert_beta_to_or <- function(gwas) {
  gwas <- get_file_or_dataframe(gwas)
  z_score <- stats::qnorm(.975, mean = 0, sd = 1) #1.96

  gwas$OR <- exp(gwas$BETA)
  gwas$OR_LB <- exp(gwas$BETA - z_score * gwas$SE)
  gwas$OR_UB <- exp(gwas$BETA + z_score * gwas$SE)
  return(gwas)
}

#' convert_z_score_to_beta
#' taken from here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8432599/ under "Correlation of trans-eQTL effects"
#' @export
convert_z_score_to_beta <- function(gwas) {
  gwas$BETA <- gwas$Z / sqrt(
    (2 * gwas$EAF) * (1 - gwas$EAF) * (gwas$N + gwas$Z^2)
  )
  gwas$SE <- 1 / sqrt(
    (2 * gwas$EAF) * (1 - gwas$EAF) * (gwas$N + gwas$Z^2)
  )

  return(gwas)
}

#TODO: don't think this is right...
# convert_beta_to_z_score <- function(gwas) {
#   mean_beta <- mean(gwas$BETA)
#   sd_beta <- sd(gwas$BETA)
#   gwas$Z <- (gwas$BETA - mean_beta) / sd_beta
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

#' @import stats
convert_z_to_p <- function(gwas) {
  gwas$P <- 2 * pnorm(-abs(gwas$Z))
  return(gwas)
}

#' @import stats
#' @export
calculate_f_statistic <- function(gwas) {
  gwas$F_STAT <- stats::qchisq(gwas$P, 1, low=F)
  return(gwas)
}

resolve_column_map <- function(column_map) {
  column_map_file <- system.file("extdata", "predefined_column_maps.csv", package = "gwaspipeline")
  if (!file.exists(column_map_file)) {
    predefined_column_maps <- vroom::vroom("../inst/extdata/predefined_column_maps.csv")
  }
  else {
    predefined_column_maps <- vroom::vroom(column_map_file, show_col_types = F)
  }
  predefined_column_maps <- tibble::column_to_rownames(predefined_column_maps, "name")

  if (is.vector(column_map) && length(column_map) > 1) {
    return(column_map)
  }
  else if (is.character(column_map) && length(column_map) == 1) {
    if (column_map %in% row.names(predefined_column_maps)) {
      predefined_map <- predefined_column_maps[column_map, ]
      return(as.list(predefined_map))
    }
    else {
      split_map <- split_string_into_named_list(column_map)
      if (length(split_map) == 0) stop(paste("Error resolving column map for", column_map))
      return(split_map)
    }
  }
  else stop(paste("Error resolving column map for", column_map))
}
