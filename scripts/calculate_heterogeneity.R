source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Caluclating heterogeneity scores between a list of GWASes")

parser <- add_argument(parser,
                       "--gwas_filenames",
                       help = "Comma separated list of GWAS filenames",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--clumped_filenames",
                       help = "Comma separated list of clumped SNP filenames",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--ancestry_list",
                       help = "Comma separated list of ancestry names",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--heterogeneity_scores_output",
                       help = "Output file name of results from heterogeneity score",
                       type = "character"
)
parser <- add_argument(parser,
                       "--heterogeneity_plot_output",
                       help = "Output file name of heterogenetiy plot",
                       type = "character"
)
parser <- add_argument(parser,
                       "--heterogeneity_plot_per_snp_output",
                       help = "Output file name of heterogeneity plot per SNP",
                       type = "character"
)

args <- parse_args(parser)

gwas_filenames <- split_string_into_vector(args$gwas_filenames)
clumped_filenames <- split_string_into_vector(args$clumped_filenames)
ancestry_list <- split_string_into_vector(args$ancestry_list)

if (length(gwas_filenames) != length(clumped_filenames) && length(gwas_filenames != length(ancestry_list))) {
  stop("Error: size of --gwas_filenames, --clumped_filenames, and --ancestry_list must be equal.")
}

create_dir_for_files(args$heterogeneity_scores_output,
                                        args$heterogeneity_plot_output,
                                        args$heterogeneity_plot_per_snp_output
)

compare_heterogeneity_across_ancestries(gwas_filenames,
                                        clumped_filenames,
                                        ancestry_list,
                                        args$heterogeneity_scores_output,
                                        args$heterogeneity_plot_output,
                                        args$heterogeneity_plot_per_snp_output
)
