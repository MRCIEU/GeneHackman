source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Correct for Collider Bias between incidence and subsequent GWASes")

parser <- add_argument(parser, "--incidence_gwas",
                       help = "Indience GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--subsequent_gwas",
                       help = "Subsequent GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--clumped_file",
                       help = "Clumped SNP list to run collider bias corrections against",
                       type = "character"
)
parser <- add_argument(parser, "--collider_bias_results_output",
                       help = "Output file comparing Collider Bias Slopes",
                       type = "character"
)
parser <- add_argument(parser, "--collider_bias_slopehunter_output",
                       help = "Output file of SlopeHunter Corrected GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--harmonised_effects_output",
                       help = "Output file of the harmonised effects file created by slopehunter, used in corrections",
                       type = "character"
)
parser <- add_argument(parser, "--p_value_thresholds",
                       help = "p value threshold to run corrections",
                       type = "character",
                       default = "0.1 0.01 0.001 1e-05",
                       nargs = Inf
)

args <- parse_args(parser)
p_value_thresholds <- as.numeric(split_string_into_vector(args$p_value_thresholds))

create_dir_for_files(args$collider_bias_results_output, args$harmonised_effects_output, args$collider_bias_slopehunter_output)

conduct_collider_bias_analysis(args$incidence_gwas,
                               args$subsequent_gwas,
                               args$clumped_file,
                               args$collider_bias_results_output,
                               args$harmonised_effects_output,
                               args$collider_bias_slopehunter_output,
                               p_value_thresholds
)
