source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Clump a standardised GWAS with the precomputed Pan-UKB LD pack (PLINK 1.9 --clump rules)")

parser <- add_argument(parser, "--gwas",
                       help = "Standardised GWAS filename",
                       type = "character"
)
parser <- add_argument(parser, "--ancestry",
                       help = "Ancestry code (EUR, AFR, AMR, CSA, EAS, MID; SAS maps to CSA)",
                       type = "character"
)
parser <- add_argument(parser, "--build",
                       help = "Reference build of the standardised GWAS (GRCh37 or GRCh38)",
                       type = "character",
                       default = "GRCh37"
)
parser <- add_argument(parser, "--clump_arguments",
                       help = "PLINK-style clumping arguments, e.g. '--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000'",
                       type = "character",
                       default = ""
)
parser <- add_argument(parser, "--output",
                       help = "Output .clumped file (PLINK layout)",
                       type = "character"
)

args <- parse_args(parser)
invisible(clump_gwas_with_ld_pack(gwas = args$gwas,
                                  ancestry = args$ancestry,
                                  build = args$build,
                                  clump_arguments = args$clump_arguments,
                                  output_file = args$output))
