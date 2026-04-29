source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Standardise a GWAS to be ready for the rest of the pipeline")

parser <- add_argument(parser, "--input_gwas",
                       help = "Comma separated list of filenames of GWASes to standardise",
                       type = "character"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Comma separated list of filenames of the standardised GWASes",
                       type = "character"
)
parser <- add_argument(parser, "--N",
                       help = "Sample size of GWAS",
                       type = "numeric",
                       default = 0
)
parser <- add_argument(parser, "--input_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--output_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--input_build",
                       help = paste(c("Input reference builds, options:", reference_builds), collapse = " "),
                       type = "character",
                       default = NULL
)
parser <- add_argument(parser, "--output_build",
                       help = paste(c("Output reference builds, options:", reference_builds), collapse = " "),
                       type = "character",
                       default = reference_builds$GRCh37
)
parser <- add_argument(parser, "--populate_rsid",
											 help = paste(c("Should GWAS Populate RSID, options:", populate_rsid_options), collapse = " "),
                       type = "character",
                       default = "none"
)
parser <- add_argument(parser, "--populate_eaf",
                       help = "If TRUE, fill missing EAF from LD reference panel (.bim + .frq under thousand genomes path)",
                       type = "logical",
                       default = FALSE
)
parser <- add_argument(parser, "--ancestry",
                       help = "Ancestry code for EAF reference (EUR, EAS, AFR, AMR, SAS); required when --populate_eaf is true",
                       type = "character",
                       default = ""
)
parser <- add_argument(parser, "--remove_extra_columns",
											 help = "Remove additional columns in GWAS that are not specified as an input column",
											 type = "logical",
											 default = FALSE
)
parser <- add_argument(parser, "--flip_alleles",
											 help = "If TRUE, reorder EA/OA alphabetically and flip BETA/EAF/Z. Set FALSE to keep alleles as-is.",
											 type = "logical",
											 default = TRUE
)

args <- parse_args(parser)
create_dir_for_files(args$output_gwas)

invisible(standardise_gwas(
                 gwas = args$input_gwas,
                 output_file = args$output_gwas,
                 N = args$N,
                 populate_rsid_option = args$populate_rsid,
                 input_reference_build = args$input_build,
                 output_reference_build = args$output_build,
                 input_columns = args$input_columns,
                 output_columns = args$output_columns,
								 remove_extra_columns = args$remove_extra_columns,
                 populate_eaf = args$populate_eaf,
                 ancestry = args$ancestry,
                 flip_alleles = args$flip_alleles
))
