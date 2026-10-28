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
                       help = paste("Ancestry code (1000genomes: EUR, EAS, AFR, AMR, SAS;",
                                    "pan_ukb: EUR, AFR, AMR, CSA/SAS, EAS, MID)"),
                       type = "character"
)
parser <- add_argument(parser, "--ld_reference",
                       help = paste("LD reference: 1000genomes (PLINK on the 1000 Genomes panel)",
                                    "or pan_ukb (precomputed LD pack)"),
                       type = "character",
                       default = "1000genomes"
)
parser <- add_argument(parser, "--build",
                       help = "Reference build of the standardised GWAS (GRCh37 or GRCh38)",
                       type = "character",
                       default = "GRCh37"
)
parser <- add_argument(parser, "--N",
                       help = "GWAS sample size",
                       type = "numeric"
)
parser <- add_argument(parser, "--output_finemap_dir",
                       help = "Directory for one <CHR>_<BP>_finemap.tsv.gz per fine-mapping region",
                       type = "character"
)
parser <- add_argument(parser, "--window_kb",
                       help = "Fine-mapping window half-width in kb (1000 = ±1 Mb)",
                       type = "numeric",
                       default = 1000
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
parser <- add_argument(parser, "--max_region_kb",
                       help = "Maximum width (kb) of a region formed by merging overlapping lead-SNP windows",
                       type = "numeric",
                       default = 3000
)
parser <- add_argument(parser, "--completion_file",
                       help = "Sentinel file written on successful completion",
                       type = "character"
)

args <- parse_args(parser)
if (!dir.exists(args$output_finemap_dir)) {
  dir.create(args$output_finemap_dir, recursive = TRUE)
}

invisible(finemap_gwas(gwas = args$gwas_filename,
             clumped_file = args$clumped_filename,
             ancestry = args$ancestry,
             default_n = args$N,
             output_finemap_dir = args$output_finemap_dir,
             completion_file = args$completion_file,
             window_kb = args$window_kb,
             max_causal = args$max_causal,
             coverage = args$coverage,
             min_abs_corr = args$min_abs_corr,
             max_region_kb = args$max_region_kb,
             ld_reference = args$ld_reference,
             build = args$build
))
