source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Mendelian randomization of an outcome GWAS against the cis QTLs of a QTL dataset")

parser <- add_argument(parser, "--gwas_filename",
                       help = "Standardised outcome GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--dataset",
                       help = paste(c("QTL dataset, options:", qtl_datasets), collapse = " "),
                       type = "character",
                       default = NULL
)
parser <- add_argument(parser, "--exposures",
                       help = "List of exposures to focus on when running MR",
                       type = "character",
                       default = "",
                       nargs = Inf
)
parser <- add_argument(parser, "--build",
                       help = "Reference build of the GWAS; QTL files of the same build are used (GRCh38 or GRCh37)",
                       type = "character",
                       default = "GRCh38"
)
parser <- add_argument(parser, "--output_file",
                       help = "Output file name of results from MR",
                       type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_file)
exposures <- split_string_into_vector(args$exposures)

perform_mr_on_qtl_dataset(args$gwas_filename, args$dataset, exposures, args$output_file, build = args$build)
