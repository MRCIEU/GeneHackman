source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Generate locus zoom plots for significant coloc results")

parser <- add_argument(parser, "--coloc_results",
                       help = "Path to coloc_results.tsv",
                       type = "character"
)
parser <- add_argument(parser, "--gwas_files",
                       help = "Space-separated list of standardised GWAS file paths",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--trait_names",
                       help = "Space-separated trait labels (same order as gwas_files)",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--ancestry",
                       help = "Ancestry for LD reference panel (EUR, EAS, etc.)",
                       type = "character",
                       default = "EUR"
)
parser <- add_argument(parser, "--pp_h4_threshold",
                       help = "Minimum PP.H4.abf to plot (default 0.8)",
                       type = "numeric",
                       default = 0.8
)
parser <- add_argument(parser, "--window_kb",
                       help = "Half-width of plotting window in kb (default 500)",
                       type = "numeric",
                       default = 500
)
parser <- add_argument(parser, "--output_dir",
                       help = "Directory to write locus zoom PNGs",
                       type = "character"
)
parser <- add_argument(parser, "--completion_file",
                       help = "Sentinel file written on successful completion",
                       type = "character"
)

args <- parse_args(parser)

gwas_files <- split_string_into_vector(args$gwas_files)
trait_names <- split_string_into_vector(args$trait_names)

if (length(trait_names) != length(gwas_files)) {
  stop("Number of trait names must match number of GWAS files")
}
names(gwas_files) <- trait_names

if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
  stop("EnsDb.Hsapiens.v75 is required for gene annotations. ",
       "Install with: BiocManager::install('EnsDb.Hsapiens.v75')")
}
ens_db <- "EnsDb.Hsapiens.v75"

locus_zoom_coloc(
  coloc_results_file = args$coloc_results,
  gwas_files = gwas_files,
  ancestry = args$ancestry,
  ens_db = ens_db,
  output_dir = args$output_dir,
  pp_h4_threshold = args$pp_h4_threshold,
  window_kb = args$window_kb,
  completion_file = args$completion_file
)
