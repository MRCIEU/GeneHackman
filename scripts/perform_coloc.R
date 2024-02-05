source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Perform a coloc analysis (and create a miami plot) for the coloc analysis")

parser <- add_argument(parser, "--first_gwas",
                       help = "filename of first GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--second_gwas",
                       help = "filename of second GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--coloc_filename",
                       help = "filename of coloc analysis to save",
                       type = "character"
)
parser <- add_argument(parser, "--miami_filename",
                       help = "filename of miami plot to save",
                       type = "character"
)
parser <- add_argument(parser, "--exposure_name",
                       help = "Name of exposure to perform coloc on",
                       default = "exposure",
                       type = "character"
)
parser <- add_argument(parser, "--chr",
                       help = "Chromosome to perform coloc on",
                       type = "double"
)
parser <- add_argument(parser, "--bp",
                       help = "bp, or position, of to perform coloc on",
                       type = "double"
)
parser <- add_argument(parser, "--range",
                       help = "bp, or position, of to perform coloc on",
                       default = 500000,
                       type = "double"
)

args <- parse_args(parser)
create_dir_for_files(args$coloc_filename, args$miami_filename)

coloc_analysis(args$first_gwas,
               args$second_gwas,
               args$coloc_filename,
               args$exposure_name,
               args$chr,
               args$bp,
               args$range
)

miami_plot(args$first_gwas,
           args$second_gwas,
           args$miami_filename,
           args$title,
           args$chr,
           args$bp,
           args$range
)
