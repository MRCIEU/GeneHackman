#' get_file_or_dataframe: function that takes either dataframe OR file name, and returns a subset of
#'   columns, rows, or both.
#'   accepted input file types: txt, csv, tsv, zip, gz
#' @param input either data frame or string input of file
#' @param columns vector of strings matching column names to filter columns
#' @param snps vector of SNP names matching SNPs to filter rows
#' @return output data frame of GWAS


#' @export
get_file_or_dataframe <- function(input, columns=NULL, snps=NULL) {
  if (is.data.frame(input)) {
    output <- dplyr::select(input, `if`(is.null(columns), dplyr::all_of(dplyr::everything()), dplyr::all_of(columns))) |>
        dplyr::filter(`if`(is.null(snps), T, SNP %in% snps))
  } else {
    if (!file.exists(input)) {
      stop(paste("Error:", input, "can't be found"))
    } else {
      if (!is.null(snps)) {
        output <- vroom_snps(input, snps) |>
            dplyr::select(`if`(is.null(columns), dplyr::all_of(dplyr::everything()), dplyr::all_of(columns)))
      } else {
        if (is.null(columns)) {
          output <- vroom::vroom(input, col_type = vroom::cols(vroom::col_character()), show_col_types=F)
        } else {
          output <- vroom::vroom(input, col_type = vroom::cols(vroom::col_character()), col_select = dplyr::all_of(columns), show_col_types=F)
        }
      }
    }
    output <- dplyr::filter(output , `if`(is.null(snps), T, SNP %in% snps))
  }
  return(output)
}


filter_gwas_by_clumped_results <- function(gwas, clumped_results) {
  #vroom has trouble reading plink --clump output, so using fread
  rsids <- data.table::fread(clumped_results, select = "SNP")$SNP
  gwas <- get_file_or_dataframe(gwas) |>
    dplyr::filter(RSID %in% rsids)
  return(gwas)
}

#' vroom_snps: If you only need to get a handful of SNPs out of a whole GWAS, this saves time and memory
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
#' @param gwas_file file of gwas to get SNPs from
#' @param snps vector of SNP names to get
#' @return dataframe of SNPs

#' @export
vroom_snps <- function(gwas_file, snps=c()){
  snps <- paste(snps, collapse="\t|")

  if (endsWith(gwas_file, ".gz")) {
    if (Sys.info()["sysname"] == "Darwin") {
      grep_command <- paste0("zcat < ", gwas_file, " | head -n 1 && rg -Iz '", snps, "' ", gwas_file)
    } else {
      grep_command <- paste0("zcat ", gwas_file, " | head -n 1 && rg -Iz '", snps, "' ", gwas_file)
    }
  } else {
    grep_command <- paste0("head -n 1", gwas_file, " && rg -I '", snps, "' ", gwas_file)
  }

  snps_in_gwas <- vroom::vroom(pipe(grep_command), col_type = vroom::cols(vroom::col_character()), show_col_types = F)
  return(snps_in_gwas)
}

#' gwas_region: filter a GWAS file to a region
#' @param gwas dataframe with the following columns: CHR, BP
#' @param chr chromosome
#' @param bp base pair position
#' @param range range to filter in base pairs
#' @return gwas with filtered rows

#' @export
gwas_region <- function(gwas, chr, bp, range = 500000) {
  return(dplyr::filter(gwas, CHR == chr &BP > (bp - floor(range/2)) & BP < (bp + floor(range/2))))
}


split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

parse_gwas_input_column_maps <- function(input_column_string) {
  column_map_as_a_string <- unlist(strsplit(input_column_string, '[:]'))
}

extract_numbers_from_string <- function(string) {
  regmatches(string, gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", string))
}

split_string_into_named_list <- function(input_string) {
  split <- unlist(strsplit(input_string, '[=,]'))
  names <- split[c(T, F)]
  values <- split[c(F, T)]

  return(structure(as.list(values), names=names))
}

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub("\\..*", "", file_name)
  file_prefix <- sub("_std", "", file_prefix)
  return(file_prefix)
}


create_dir_for_files <- function(...) {
  filenames <- list(...)

  for (filename in filenames) {
    filepath <- file.path(stringr::str_extract(filename, "^(.*/)"))
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  }
}

map_rsid_list_to_snps <- function(gwas, rsids=c()) {
  gwas <- subset(gwas, RSID %in% rsids)
  return(gwas$SNP)
}




create_html_from_rmd <- function(rmd_file, params = list(), output_file) {
  temp_file <- tempfile(fileext = ".Rmd", tmpdir = "/tmp")
  file.copy(rmd_file, temp_file, overwrite = TRUE)

  rmarkdown::render(temp_file,
                    output_file = output_file,
                    params = params,
                    quiet = T
  )
}

get_docker_image_tag <- function() {
  docker_version <- get_env_var("DOCKER_VERSION")
  if (is.null(docker_version)) {
    return("latest")
  } else {
    return(docker_version)
  }
  #return(packageVersion("GeneHackman"))
}

#' Generate log Bayes Factor from Z-score
#'
#' @param z Z-score
#' @param se Standard error
#' @param eaf Allele frequency
#' @param sample_size Sample size
#' @param study_type Study type
#' @param effect_priors Effect priors
#'
#' @return Log Bayes Factor
#' @internal
convert_z_to_lbf <- function(
  z,
  se,
  eaf,
  sample_size,
  study_type,
  effect_priors = c(continuous = 0.15, categorical = 0.2)
) {
  estimated_sd <- estimate_variance(se, eaf, sample_size)
  if (study_type == study_types$continuous) {
    sd_prior <- effect_priors[study_types$continuous] * estimated_sd
  } else {
    sd_prior <- effect_priors[study_types$categorical]
  }
  r <- sd_prior^2 / (sd_prior^2 + se^2)
  lbf <- (log(1 - r) + (r * z^2)) / 2
  return(lbf)
}

#' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
#'
#' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
#' var(X) = 2*maf*(1-maf)
#' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
#'
#' @title Estimate trait variance, internal function
#' @param SE vector of standard errors
#' @param EAF vector of MAF (same length as SE)
#' @param n sample size
#' 
#' @return estimated standard deviation of Y
#' @internal
estimate_variance <- function(se, eaf, n) {
  oneover <- 1 / se^2
  nvx <- 2 * n * eaf * (1 - eaf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[["oneover"]]
  if (cf < 0) {
    stop(
      "Estimated sdY is negative - this can happen with small datasets, or those with errors. ",
      "A reasonable estimate of sdY is required to continue."
    )
  }
  return(sqrt(cf))
}