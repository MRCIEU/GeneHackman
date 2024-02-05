source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Create a manhattan plot from GWAS summary statistics")

parser <- add_argument(parser, "--first_gwas", help = "filename of first GWAS", type = "character")
parser <- add_argument(parser, "--second_gwas", help = "filename of second GWAS", type = "character")
parser <- add_argument(parser, "--miami_filename", help = "filename of plot to save", type = "character")
parser <- add_argument(parser, "--title", help = "title of plot", default = "Comparing GWASes", type = "character")

args <- parse_args(parser)

miami_plot(args$first_gwas, args$second_gwas, args$miami_filename, args$title)
