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
parser <- add_argument(parser, "--ld_reference",
                       help = "LD reference: 1000genomes or pan_ukb",
                       type = "character",
                       default = "1000genomes"
)
parser <- add_argument(parser, "--build",
                       help = "Reference build of the standardised GWAS (GRCh38 or GRCh37); selects the gene annotation",
                       type = "character",
                       default = "GRCh38"
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

ens_db <- ensembl_db_for_build(args$build)
if (!requireNamespace(ens_db, quietly = TRUE)) {
  stop(ens_db, " is required for ", args$build, " gene annotations. ",
       "Install with: BiocManager::install('", ens_db, "')")
}
# locuszoomr requires the EnsDb package to be loaded (not just installed); see ?locuszoomr::locus
suppressPackageStartupMessages(library(ens_db, character.only = TRUE))

locus_zoom_coloc(
  coloc_results_file = args$coloc_results,
  gwas_files = gwas_files,
  ancestry = args$ancestry,
  ens_db = ens_db,
  output_dir = args$output_dir,
  pp_h4_threshold = args$pp_h4_threshold,
  window_kb = args$window_kb,
  completion_file = args$completion_file,
  ld_reference = args$ld_reference,
  build = args$build
)
