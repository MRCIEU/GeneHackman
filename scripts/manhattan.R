source("load.R")

library(argparser, quietly = TRUE)

parser <- arg_parser("Create a miami plot, comparing 2 GWASes")

parser <- add_argument(parser, "--gwas_filename", help = "filename of GWAS", type = "character")
parser <- add_argument(parser, "--manhattan_filename", help = "filename of plot", type = "character")
parser <- add_argument(parser, "--qq_filename", help = "filename of qq plot", type = "character")
parser <- add_argument(parser, "--include_qq", help = "include qq plot", default = TRUE, type = "logical")

args <- parse_args(parser)

create_dir_for_files(args$manhattan_filename, args$qq_filename)
manhattan_and_qq(args$gwas_filename, args$manhattan_filename, args$qq_filename, args$include_qq)
