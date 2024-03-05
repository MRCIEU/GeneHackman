source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Perform a coloc analysis (and create a miami plot) for the coloc analysis")

parser <- add_argument(parser, "--mr_results_filename", help = "filename of first GWAS", type = "character" )
parser <- add_argument(parser, "--gwas_filename",
                       help = "filename of first GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--N",
                       help = "Sample size of GWAS",
                       type = "numeric",
                       default = 0
)
parser <- add_argument(parser, "--exposures",
                       help = "List of exposures to perform coloc on, if none provided, runs on exposures with lowest p-vals",
                       type = "character",
                       default = "",
                       nargs = Inf
)
parser <- add_argument(parser, "--dataset",
                       help = paste(c("QTL dataset, options:", qtl_datasets), collapse = " "),
                       default = "exposure",
                       type = "character"
)
parser <- add_argument(parser, "--output_file",
                       help = "filename of coloc analysis to save",
                       type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_file, paste0(user_results_dir, "/plots"))
exposures <- split_string_into_vector(args$exposures)

run_coloc_on_qtl_mr_results(args$mr_results_filename, args$gwas_filename, args$dataset, exposures, output_file=args$output_file, default_n=args$N)
