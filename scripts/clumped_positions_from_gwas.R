source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Replace the panel positions PLINK writes in a .clumped file with the GWAS's own positions")
parser <- add_argument(parser, "--clumped", help = "PLINK .clumped file (rewritten in place)", type = "character")
parser <- add_argument(parser, "--gwas", help = "Standardised GWAS the clumps came from", type = "character")
args <- parse_args(parser)

invisible(use_gwas_positions_in_clumped(args$clumped, args$gwas))
