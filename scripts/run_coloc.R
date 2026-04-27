source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Run BF-BF colocalization on finemapped GWAS results")

parser <- add_argument(parser, "--finemap_dirs",
                       help = "Space-separated list of finemap output directories",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--trait_names",
                       help = "Space-separated trait labels (same order as finemap_dirs)",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--overlap_kb",
                       help = "Distance in kb to define overlapping signals (default 1000 = ±1 Mb)",
                       type = "numeric",
                       default = 1000
)
parser <- add_argument(parser, "--p1",
                       help = "Prior: SNP associated with trait 1",
                       type = "numeric",
                       default = 1e-4
)
parser <- add_argument(parser, "--p2",
                       help = "Prior: SNP associated with trait 2",
                       type = "numeric",
                       default = 1e-4
)
parser <- add_argument(parser, "--p12",
                       help = "Prior: SNP associated with both traits",
                       type = "numeric",
                       default = 5e-6
)
parser <- add_argument(parser, "--output_file",
                       help = "Output file for coloc results",
                       type = "character"
)

args <- parse_args(parser)

finemap_dirs <- split_string_into_vector(args$finemap_dirs)
trait_names <- split_string_into_vector(args$trait_names)

if (length(trait_names) != length(finemap_dirs)) {
  stop("Number of trait names must match number of finemap directories")
}
names(finemap_dirs) <- trait_names

create_dir_for_files(args$output_file)

run_bf_bf_coloc(
  finemap_dirs = finemap_dirs,
  overlap_kb = args$overlap_kb,
  p1 = args$p1,
  p2 = args$p2,
  p12 = args$p12,
  output_file = args$output_file
)
