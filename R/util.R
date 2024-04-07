#' get_file_or_dataframe: function that takes either dataframe OR file name, and returns a subset of
#'   columns, rows, or both.
#'   accepted input file types: txt, csv, tsv, zip, gz
#' @param input: either data frame or string input of file
#' @param columns: vector of strings matching column names to filter columns
#' @param snps: vector of SNP names matching SNPs to filter rows
#' @return output: data frame of GWAS
#' @import dplyr
#' @import vroom
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
#' @import data.table
#' @import dplyr
filter_gwas_by_clumped_results <- function(gwas, clumped_results) {
  #vroom has trouble reading plink --clump output, so using fread
  rsids <- data.table::fread(clumped_results, select = "SNP")$SNP
  gwas <- get_file_or_dataframe(gwas) |>
    dplyr::filter(RSID %in% rsids)
  return(gwas)
}

#' vroom_snps: If you only need to get a handful of SNPs out of a whole GWAS, this saves time and memory
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
#' @import vroom
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

#' @import dplyr
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

#' @import stringr
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

#' @import rmarkdown
#' @import httr
#' @import prettydoc
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
  return("1.0.0")
  #return(packageVersion("GeneHackman"))
}

#' @import vroom
run_sqlite_command <- function(db_name, query, col_names = c()) {
  query <- paste0('"', query, '"')
  sqlite_command <- paste("sqlite3", db_name, query)
  output <- system(sqlite_command, intern = T, wait = T)
  result <- vroom::vroom(output, col_names = F, col_type = vroom::cols(vroom::col_character()), show_col_types=F)

  if (is.vector(col_names)) colnames(result) <- col_names
  return(result)
}
