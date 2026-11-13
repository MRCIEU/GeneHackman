source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Clump a standardised GWAS with the LD reference of the clumping step")

parser <- add_argument(parser, "--gwas",
                       help = "Standardised GWAS filename",
                       type = "character"
)
parser <- add_argument(parser, "--ancestry",
                       help = paste("Ancestry code (1000genomes: EUR, EAS, AFR, AMR, SAS;",
                                    "pan_ukb: EUR, AFR, AMR, CSA or SAS, EAS, MID)"),
                       type = "character"
)
parser <- add_argument(parser, "--ld_reference",
                       help = paste("pan_ukb (Pan-UKB LD pack), 1000genomes (1000 Genomes panel)",
                                    "or the PLINK prefix of your own genotypes"),
                       type = "character"
)
parser <- add_argument(parser, "--build",
                       help = "Reference build of the standardised GWAS (GRCh37 or GRCh38)",
                       type = "character",
                       default = "GRCh38"
)
parser <- add_argument(parser, "--clump_arguments",
                       help = "PLINK-style clumping arguments, e.g. '--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000'",
                       type = "character",
                       default = ""
)
parser <- add_argument(parser, "--threads",
                       help = "PLINK threads",
                       type = "numeric",
                       default = 1
)
parser <- add_argument(parser, "--output",
                       help = "Output .clumped file (PLINK layout)",
                       type = "character"
)

args <- parse_args(parser)
if (is.na(args$ld_reference) || !nzchar(args$ld_reference)) stop("--ld_reference is required")
invisible(clump_gwas(gwas = args$gwas,
                     ancestry = args$ancestry,
                     ld_reference = args$ld_reference,
                     build = args$build,
                     clump_arguments = trimws(args$clump_arguments),
                     output_file = args$output,
                     threads = args$threads))
