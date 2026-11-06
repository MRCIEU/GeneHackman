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
                       default = reference_builds$GRCh38
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
                       help = paste("Ancestry code for EAF reference (EUR, EAS, AFR, AMR, SAS; pan_ukb also CSA, MID);",
                                    "required when --populate_eaf is true"),
                       type = "character",
                       default = ""
)
parser <- add_argument(parser, "--ld_reference",
                       help = paste("Reference used to fill EAF: 1000genomes (.frq, matched on RSID)",
                                    "or pan_ukb (LD pack, matched on SNP key)"),
                       type = "character",
                       default = "1000genomes"
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
parser <- add_argument(parser, "--complement_minus_strand",
                       help = "Liftover: keep non-palindromic minus-strand variants with complemented alleles",
                       type = "logical",
                       default = FALSE
)
parser <- add_argument(parser, "--max_liftover_drop_fraction",
                       help = "Liftover: stop if more than this fraction of variants is dropped",
                       type = "numeric",
                       default = 0.05
)
parser <- add_argument(parser, "--check_build",
                       help = "Check positions against the Pan-UKB LD pack's GRCh37/GRCh38 keys before lifting",
                       type = "logical",
                       default = TRUE
)
parser <- add_argument(parser, "--study_type",
                       help = "continuous or binary (case-control); binary traits get an effective sample size column N_EFF",
                       type = "character",
                       default = "continuous"
)
parser <- add_argument(parser, "--N_case",
                       help = "Number of cases (binary traits), used for N_EFF when the file has no N_CASES/N_CONTROLS",
                       type = "numeric",
                       default = NA
)
parser <- add_argument(parser, "--N_control",
                       help = "Number of controls (binary traits)",
                       type = "numeric",
                       default = NA
)
parser <- add_argument(parser, "--strand_check",
                       help = "Check alleles and strand against the LD reference (needs --ancestry) and fix strand flips",
                       type = "logical",
                       default = TRUE
)
parser <- add_argument(parser, "--drop_ambiguous_palindromic",
                       help = "Remove A/T and C/G SNVs whose strand can't be told from allele frequency",
                       type = "logical",
                       default = FALSE
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
                 flip_alleles = args$flip_alleles,
                 ld_reference = args$ld_reference,
                 complement_minus_strand = args$complement_minus_strand,
                 max_liftover_drop_fraction = args$max_liftover_drop_fraction,
                 check_build = args$check_build,
                 study_type = args$study_type,
                 n_case = args$N_case,
                 n_control = args$N_control,
                 strand_check = args$strand_check,
                 drop_ambiguous_palindromic = args$drop_ambiguous_palindromic
))
