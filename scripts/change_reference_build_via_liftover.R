source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Convert reference build (eg. )")

parser <- add_argument(parser, "--input_gwas",
                       help = "Filename of GWASes to change reference build",
                       type = "character"
)
parser <- add_argument(parser, "--input_reference_build",
											 help = paste(c("Input reference builds, options:", reference_builds), collapse = " "),
                       type = "character"
)
parser <- add_argument(parser, "--output_reference_build",
					   help = paste(c("Output reference builds, options:", reference_builds), collapse = " "),
                       type = "character"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Filename of GWASes to save changes to",
                       type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_gwas)

convert_reference_build_via_liftover(args$input_gwas,
                                     input_reference_build = args$input_reference_build,
                                     output_reference_build = args$output_reference_build,
                                     output_file = args$output_gwas
)
