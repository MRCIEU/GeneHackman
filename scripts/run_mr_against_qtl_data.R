source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Comparison of expected vs observed replication rates")

parser <- add_argument(parser, "--gwas_filename",
                       help = "Comma separated list of GWAS filenames",
                       type = "character"
)
parser <- add_argument(parser, "--ancestry",
                       help = "Comma separated list of ancestries related to GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--dataset",
                       help = paste(c("QTL dataset, options:", qtl_datasets), collapse = " "),
                       type = "character",
                       default = NULL
)
parser <- add_argument(parser, "--subcategory",
                       help = "Subcategory of type of QTL analysis to run (depends on dataset)",
                       type = "character",
                       default = NULL
)
parser <- add_argument(parser, "--exposures",
                       help = "List of exposures to focus on when running MR",
                       type = "character",
                       default = "",
                       nargs = Inf
)
parser <- add_argument(parser, "--output_file",
                       help = "Output file name of results from MR",
                       type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_file)
exposures <- split_string_into_vector(args$exposures)

if (args$dataset == qtl_datasets$metabrain) {
  perform_mr_on_metabrain_datasets(args$gwas_filename, args$ancestry, args$subcategory, exposures, args$output_file)
}
else if (args$dataset == qtl_datasets$eqtlgen) {
  perform_mr_on_eqtlgen_datasets(args$gwas_filename, args$subcategory, exposures, args$output_file)
}
