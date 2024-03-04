source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Create a volcano plot from results file (typically MR results)")

parser <- add_argument(parser, "--results_filename",
             help = "Results file to make plot",
             type = "character"
)
parser <- add_argument(parser, "--title",
             help = "filename of first GWAS",
             type = "character",
             default = "Volcano Plot of MR Results"
)
parser <- add_argument(parser, "--label",
             help = "column name to be used as label name in plot",
             type = "character",
             default = "exposure"
)
parser <- add_argument(parser, "--num_labels",
             help = "Number of labels to include in volcano plot",
             type = "double",
             default = 30
)
parser <- add_argument(parser, "--output_file",
             help = "filename of plot to save",
             type = "character"
)

args <- parse_args(parser)

create_dir_for_files(args$output_file)
volcano_plot(args$results_file, args$title, label=args$label, num_labels = args$num_labels, args$output_file)
