source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Run  colocalization on significant MR results using finemapped GWAS data")

parser <- add_argument(parser, "--mr_results_filename", help = "MR results filename", type = "character")
parser <- add_argument(parser, "--finemap_dir",
                       help = "Directory containing finemapped GWAS locus files",
                       type = "character"
)
parser <- add_argument(parser, "--N",
                       help = "Sample size of the QTL dataset, used only when the QTL files have no N column",
                       type = "numeric",
                       default = 0
)
parser <- add_argument(parser, "--study_type",
                       help = "Study type of the QTL trait for LBF conversion (continuous or categorical)",
                       type = "character",
                       default = "continuous"
)
parser <- add_argument(parser, "--exposures",
                       help = "List of exposures to perform coloc on",
                       type = "character",
                       default = "",
                       nargs = Inf
)
parser <- add_argument(parser, "--dataset",
                       help = paste(c("QTL dataset, options:", qtl_datasets), collapse = " "),
                       default = "exposure",
                       type = "character"
)
parser <- add_argument(parser, "--p1",
                       help = "Prior: SNP associated with the GWAS trait",
                       type = "numeric",
                       default = 1e-4
)
parser <- add_argument(parser, "--p2",
                       help = "Prior: SNP associated with the QTL trait",
                       type = "numeric",
                       default = 1e-4
)
parser <- add_argument(parser, "--p12",
                       help = "Prior: SNP associated with both traits",
                       type = "numeric",
                       default = 5e-6
)
parser <- add_argument(parser, "--output_file",
                       help = "Output filename for coloc results",
                       type = "character"
)

args <- parse_args(parser)
create_dir_for_files(args$output_file, paste0(user_results_dir, "/plots"))
exposures <- split_string_into_vector(args$exposures)

run_coloc_on_qtl_mr_results(
  mr_results_file = args$mr_results_filename,
  finemap_dir = args$finemap_dir,
  qtl_dataset = args$dataset,
  study_type = args$study_type,
  exposures = exposures,
  default_n = args$N,
  output_file = args$output_file,
  p1 = args$p1,
  p2 = args$p2,
  p12 = args$p12
)
