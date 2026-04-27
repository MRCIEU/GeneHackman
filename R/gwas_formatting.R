#' standardise_gwas: takes an input gwas, changes headers, standardises allelic input, adds RSID, makes life easier
#' @param gwas filename of gwas to standardise or dataframe of gwas
#' @param output_file file to save standardised gwas
#' @param N sample size of GWAS (if GWAS has defined N column, defaults to that value)
#' @param populate_rsid_option if you want RSID populated or not
#' @param input_reference_build reference build of CHR and BP of data.  Defaults to GRCh37
#' @param output_reference_build reference build of CHR and BP of data.  Defaults to GRCh37
#' @param input_columns column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @param output_columns column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @param remove_extra_columns logical flag on if to remove extra columns
#' @param populate_eaf If TRUE, fill missing \code{EAF} from the 1000 Genomes LD panel (\code{ANCESTRY.bim} + \code{ANCESTRY.frq} under \code{THOUSAND_GENOMES_DIR}).
#' @param ancestry Ancestry code (EUR, EAS, AFR, AMR, SAS) matching the reference panel prefix; required when \code{populate_eaf} is TRUE.
#' @param flip_alleles If TRUE (default), reorder EA/OA alphabetically and flip BETA, EAF, Z accordingly.
#'   Set to FALSE to keep alleles as-is (useful when you only need to populate RSIDs, EAF, or change build).
#' @return modified gwas: saves new gwas in output_file if present


#' @export
standardise_gwas <- function(gwas,
                             output_file,
                             N=0,
                             populate_rsid_option=populate_rsid_options$none,
                             input_reference_build=reference_builds$GRCh37,
                             output_reference_build=reference_builds$GRCh37,
                             input_columns="default",
                             output_columns="default",
                             remove_extra_columns=F,
                             populate_eaf=FALSE,
                             ancestry=NULL,
                             flip_alleles=TRUE) {
  input_gwas_columns <- resolve_column_map(input_columns)
  output_gwas_columns <- resolve_column_map(output_columns)

  if (populate_eaf && (is.null(ancestry) || length(ancestry) == 0 || !nzchar(as.character(ancestry)[1]))) {
    stop("populate_eaf is TRUE but ancestry is missing; set ancestry to EUR, EAS, AFR, AMR, or SAS to match your LD reference panel.")
  }
  if (!flip_alleles && identical(populate_rsid_option, populate_rsid_options$partial)) {
    stop("Partial RSID population requires alphabetically ordered alleles (the lookup table stores CHR:BP_EA_OA with EA < OA). ",
         "Either set flip_alleles=TRUE, or use populate_rsid='full' which matches on CHR+BP+alleles directly.")
  }

  gwas <- get_file_or_dataframe(gwas) |>
    change_column_names(input_gwas_columns, remove_extra_columns) |>
    standardise_columns(N) |>
    filter_incomplete_rows() |>
    convert_reference_build_via_liftover(input_reference_build, output_reference_build) |>
    standardise_alleles(flip = flip_alleles) |>
    health_check() |>
    populate_eaf_from_reference_panel(ancestry, populate_eaf) |>
    populate_rsid(populate_rsid_option) |>
    populate_gene_names() |>
    change_column_names(output_gwas_columns)

  if (!missing(output_file) && shiny::isTruthy(output_file)) {
    vroom::vroom_write(gwas, output_file)
  }
  return(gwas)
}

#' change_snp_identifiers: changes the SNP identifiers of a GWAS to the output reference build
#' @param gwas filename of gwas to standardise or dataframe of gwas
#' @param output_file file to save standardised gwas
#' @param populate_rsid_option if you want RSID populated or not
#' @param input_reference_build reference build of CHR and BP of data.  Defaults to GRCh37
#' @param output_reference_build reference build of CHR and BP of data.  Defaults to GRCh37
#' @param input_columns column header map used for renaming GWAS:
#'   can be list() of key/value pairs, string of row in predefined_column_maps, or comma separated list of keys=values
#' @return modified gwas: saves new gwas in output_file if present

#' @export
change_snp_identifiers <- function(gwas,
                                   output_file,
                                   populate_rsid_option=populate_rsid_options$none,
                                   input_reference_build=reference_builds$GRCh37,
                                   output_reference_build=reference_builds$GRCh37,
                                   input_columns="default") {
  input_gwas_columns <- resolve_column_map(input_columns)

  # Load the original GWAS (can be a filename or a dataframe)
  original_gwas <- get_file_or_dataframe(gwas)

  # Track original row order so we can map new annotations back
  original_gwas$row_number <- seq_len(nrow(original_gwas))

  # Create an updated version with new coordinates / RSIDs, keeping row_number
  updated_gwas <- original_gwas |>
    change_column_names(input_gwas_columns) |>
    convert_reference_build_via_liftover(input_reference_build, output_reference_build) |>
    standardise_alleles() |>
    populate_rsid(populate_rsid_option) |>
    populate_gene_names()

  # Map new RSID and BP (if present) back onto the original data using row_number
  annotation_cols <- c("RSID", "BP")
  annotation_cols <- annotation_cols[annotation_cols %in% colnames(updated_gwas)]

  if (length(annotation_cols) > 0) {
    updated_subset <- updated_gwas[, c("row_number", annotation_cols), drop = FALSE]
    for (col in annotation_cols) {
      original_gwas[[col]] <- updated_subset[[col]][match(original_gwas$row_number, updated_subset$row_number)]
    }
  }

  # Remove helper row_number column
  original_gwas$row_number <- NULL

  if (!missing(output_file) && shiny::isTruthy(output_file)) {
    vroom::vroom_write(original_gwas, output_file)
  }

  return(original_gwas)
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


standardise_columns <- function(gwas, N) {
  gwas_columns <- colnames(gwas)

  if (N > 0) {
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
  if ("BETA" %in% gwas_columns) {
    gwas$BETA <- as.numeric(gwas$BETA)
  }

  if ("BETA" %in% gwas_columns) {
    gwas$BETA <- as.numeric(gwas$BETA)
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

#' @param flip If TRUE (default), reorder EA/OA into alphabetical order and flip
#'   BETA, EAF, and Z to match.  When FALSE, alleles are uppercased and the SNP
#'   column is rebuilt, but allele order (and effect direction) is left as-is.

#' @keywords internal
standardise_alleles <- function(gwas, flip = TRUE) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  if (flip) {
    to_flip <- (gwas$EA > gwas$OA) & (!gwas$EA %in% c("D", "I"))
    if (any(to_flip)) {
      gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
      gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]
      if ('Z' %in% names(gwas)) {
        gwas$Z[to_flip] <- -1 * gwas$Z[to_flip]
      }

      temp <- gwas$OA[to_flip]
      gwas$OA[to_flip] <- gwas$EA[to_flip]
      gwas$EA[to_flip] <- temp
    }
  }

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  gwas <- dplyr::select(gwas, SNP, CHR, BP, EA, OA, dplyr::everything())

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

convert_z_score_to_beta <- function(gwas, sample_size) {
  gwas$SE_new <- 1 / sqrt(2 * gwas$EAF * (1 - gwas$EAF) * sample_size)
  gwas$BETA_new <- gwas$Z * gwas$SE_new

  correction_gwas <- dplyr::filter(gwas, !is.na(BETA) & !is.null(BETA))
  correction_gwas$BETA_new[which(!is.finite(correction_gwas$BETA_new))] <- NA
  correction <- lm(correction_gwas$BETA_new ~ correction_gwas$BETA, na.action=na.omit)$coef[2]

  gwas$BETA_new <- gwas$BETA_new / correction
  gwas$SE_new <- gwas$SE_new / correction
  gwas$P_new <- abs(2 * pnorm(abs(gwas$Z), lower.tail = F))

  gwas <- dplyr::mutate(gwas, BETA = dplyr::if_else(is.na(BETA), BETA_new, BETA),
                              SE = dplyr::if_else(is.na(SE), SE_new, SE),
                              P = dplyr::if_else(is.na(P), P_new, P) ) |>
    dplyr::select(-BETA_new, -SE_new, -P_new) |>
    dplyr::filter(!is.na(BETA) & !is.na(SE) & BETA != Inf & SE != Inf)

  return(gwas)
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
