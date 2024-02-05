source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Standardise a GWAS to be ready for the rest of the pipeline")

parser <- add_argument(parser, "--input_gwas",
                       help = "Comma separated list of filenames of GWASes to standardise",
                       type = "character"
)
parser <- add_argument(parser, "--N",
                       help = "Sample size of GWAS",
                       type = "numeric",
                       default = 0
)
parser <- add_argument(parser, "--input_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character"
)
parser <- add_argument(parser, "--input_format",
                       help = "Input format of the gwas (ie. metal, bolt, plink, default)",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--output_format",
                       help = "Output format of the gwas (ie. metal, bolt, plink, default)",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Comma separated list of filenames of the standardised GWASes",
                       type = "character"
)
parser <- add_argument(parser, "--input_build",
                       help = paste(c("Input reference builds, options:", reference_builds), collapse = " "),
                       type = "character",
                       default = NULL
)
parser <- add_argument(parser, "--output_build",
                       help = paste(c("Output reference builds, options:", reference_builds), collapse = " "),
                       type = "character",
                       default = "GRCh37"
)
parser <- add_argument(parser, "--populate_rsid",
                       help = "Should GWAS Populate RSID (based on 1000 genomes data)",
                       type = "logical",
                       default = F
)
parser <- add_argument(parser, "--output_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character",
                       default = NULL
)

args <- parse_args(parser)
create_dir_for_files(args$output_gwas)

input_column_map <-split_string_into_named_list(args$input_columns)
output_column_map <-split_string_into_named_list(args$output_columns)

standardise_gwas(args$input_gwas,
                 args$output_gwas,
                 N = args$N,
                 input_format = args$input_format,
                 populate_rsid_option = args$populate_rsid,
                 input_reference_build = args$input_build,
                 output_reference_build = args$output_build,
                 input_column_map = input_column_map,
                 output_column_map = output_column_map
)
