source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Comparison of expected vs observed replication rates")

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
parser <- add_argument(parser,
                       "--variants_output",
                       help = "Output file name of variants from gas comparison",
                       type = "character"
)

args <- parse_args(parser)

gwas_filenames <- split_string_into_vector(args$gwas_filenames)
clumped_filenames <- split_string_into_vector(args$clumped_filenames)

if (length(gwas_filenames) - length(clumped_filenames) > 1) {
  stop("There are not enough clumped files (compared to gwases) to successfully run this operation.")
}

create_dir_for_files(args$result_output, args$variants_output)
compare_replication_across_all_gwas_permutations(gwas_filenames, clumped_filenames, args$result_output, args$variants_output)
