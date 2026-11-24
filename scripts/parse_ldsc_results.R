source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Tables of LDSC heritability and genetic correlation results, with reliability flags")
parser <- add_argument(parser, "--logs", help = "Comma separated LDSC logs (*_h2_<i>.log and *_rg_<i>.log)", type = "character")
parser <- add_argument(parser, "--sumstats", help = "Comma separated munged sumstats, in GWAS order", type = "character")
parser <- add_argument(parser, "--labels", help = "Comma separated trait labels, in GWAS order", type = "character")
parser <- add_argument(parser, "--population_prevalences", help = "Comma separated population prevalences (nan for none)",
                       type = "character", default = "")
parser <- add_argument(parser, "--ancestry", help = "Ancestry of the LD scores", type = "character")
parser <- add_argument(parser, "--h2_output", help = "Heritability table", type = "character")
parser <- add_argument(parser, "--rg_output", help = "Genetic correlation table", type = "character")
args <- parse_args(parser)

split <- function(x) if (is.na(x) || !nzchar(x)) character() else strsplit(x, ",", fixed = TRUE)[[1]]
logs <- split(args$logs)
sumstats <- split(args$sumstats)
prevalences <- split(args$population_prevalences)
if (length(prevalences) == 0) prevalences <- rep(NA, length(sumstats))

tables <- ldsc_results_tables(h2_logs = grep("_h2_[0-9]+\\.log$", logs, value = TRUE),
                              rg_logs = grep("_rg_[0-9]+\\.log$", logs, value = TRUE),
                              sumstats = sumstats, labels = split(args$labels), ancestry = args$ancestry,
                              population_prevalences = prevalences)
data.table::fwrite(tables$h2, args$h2_output, sep = "\t")
data.table::fwrite(tables$rg, args$rg_output, sep = "\t")
if (any(tables$h2$FLAG_LOW_H2_Z)) {
  message("Heritability z-score below 4 for: ", paste(tables$h2$TRAIT[tables$h2$FLAG_LOW_H2_Z], collapse = ", "),
          "; their genetic correlations are flagged as unreliable")
}
