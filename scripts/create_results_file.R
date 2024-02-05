source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Create a markdown file of results")

parser <- add_argument(parser, "--rmd_file",
                       help = "name of rmd file to generate from (choose on in markdown folder)",
                       type = "character"
)
parser <- add_argument(parser, "--params",
                       help = "Params to pass into markdown, a=b,c=d",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--output_file",
                            help = "Name of output file",
                            type = "character"
)

args <- parse_args(parser)

rmd_params <- split_string_into_named_list(args$params)
create_html_from_rmd(args$rmd_file, rmd_params, args$output_file)
