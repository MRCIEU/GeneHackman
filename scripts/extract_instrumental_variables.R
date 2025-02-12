source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Extract instrumental variables from GWASes")

parser <- add_argument(parser,
                       "--gwas_filenames",
                       help = "Comma separated list of GWAS filenames",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--clumped_filenames",
                       help = "Comma separated list of clumped SNP filenames",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--result_output",
                       help = "Output file name of results from gwas comparison",
                       type = "character"
)

args <- parse_args(parser)

gwas_filenames <- split_string_into_vector(args$gwas_filenames)
clumped_filenames <- split_string_into_vector(args$clumped_filenames)

if (length(gwas_filenames) - length(clumped_filenames) > 1) {
  stop("There are not enough clumped files (compared to gwases) to successfully run this operation.")
}

# Initialize list to store results
instrumental_variables_list <- list()

# Extract instruments for each GWAS
for (i in seq_along(gwas_filenames)) {
  clumped_file <- if (i <= length(clumped_filenames)) clumped_filenames[i] else NULL
  
  instrumental_variables_list[[i]] <- extract_instrumental_variables(
    gwas_filename = gwas_filenames[i],
    clumped_filename = clumped_file
  )
}
