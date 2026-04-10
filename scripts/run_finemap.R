source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Run SuSiE fine-mapping on clumped loci from a GWAS")

parser <- add_argument(parser, "--gwas_filename",
                       help = "Standardised GWAS filename",
                       type = "character"
)
parser <- add_argument(parser, "--clumped_filename",
                       help = "plink --clump output file",
                       type = "character"
)
parser <- add_argument(parser, "--ancestry",
                       help = "Ancestry code matching 1000 Genomes panel (EUR, EAS, AFR, AMR, SAS)",
                       type = "character"
)
parser <- add_argument(parser, "--N",
                       help = "GWAS sample size",
                       type = "numeric"
)
parser <- add_argument(parser, "--output_lbf_file",
                       help = "Output file for per-SNP log Bayes factors across all loci",
                       type = "character"
)
parser <- add_argument(parser, "--output_credible_set_file",
                       help = "Output file for credible-set-filtered GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--window_kb",
                       help = "Fine-mapping window half-width in kb",
                       type = "numeric",
                       default = 500
)
parser <- add_argument(parser, "--max_causal",
                       help = "Maximum number of causal signals per locus (SuSiE L)",
                       type = "numeric",
                       default = 10
)
parser <- add_argument(parser, "--coverage",
                       help = "Credible set coverage",
                       type = "numeric",
                       default = 0.95
)
parser <- add_argument(parser, "--min_abs_corr",
                       help = "Minimum absolute correlation for credible set purity",
                       type = "numeric",
                       default = 0.5
)

args <- parse_args(parser)
create_dir_for_files(args$output_lbf_file, args$output_credible_set_file)

finemap_gwas(gwas = args$gwas_filename,
             clumped_file = args$clumped_filename,
             ancestry = args$ancestry,
             n = args$N,
             output_lbf_file = args$output_lbf_file,
             output_credible_set_file = args$output_credible_set_file,
             window_kb = args$window_kb,
             max_causal = args$max_causal,
             coverage = args$coverage,
             min_abs_corr = args$min_abs_corr
)
