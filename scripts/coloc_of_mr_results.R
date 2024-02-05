source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Perform a coloc analysis (and create a miami plot) for the coloc analysis")

parser <- add_argument(parser, "--mr_results_filename",
					   help = "filename of first GWAS",
					   type = "character"
)
parser <- add_argument(parser, "--gwas_filename",
                       help = "filename of first GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--exposures",
                       help = "List of exposures to perform coloc on, if none provided, runs on exposures with lowest p-vals",
                       type = "character",
					   default = "",
                       nargs = Inf
)
parser <- add_argument(parser, "--qtl_dataset",
                       help = "Name of mr dataset (metabrain, pqtl)",
                       default = "exposure",
                       type = "character"
)
parser <- add_argument(parser, "--output_file",
					   help = "filename of coloc analysis to save",
					   type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_file, paste0(Sys.getenv("RESULTS_DIR"), "/plots"))
exposures <- split_string_into_vector(args$exposures)

run_coloc_on_qtl_mr_results(args$mr_results_filename, args$gwas_filename, args$qtl_dataset, args$ancestry, exposures, args$output_file)
