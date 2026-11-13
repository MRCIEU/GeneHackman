# Concordance of GeneHackman's harmonisation with MungeSumstats (Murphy et al. 2021) on the same GWAS file,
# for the response to Reviewers 1 (Major 3) and 2 (comment 1).  See scripts/validation/README.md.
#
# Both tools standardise the file to --output_build and assign rsIDs; variants are matched on chromosome, position
# and the unordered allele pair (compare_harmonisation() in R/validation.R).  With --chrpos_only the input is first
# rewritten so variants are identified only by chr:pos:a1:a2 (no CHR, BP or RSID columns), which is the case
# Reviewer 1 reported GeneHackman failed on.
#
# Writes <output_prefix>_summary.tsv, <output_prefix>_discordant.tsv.gz and both tools' outputs.
source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Compare GeneHackman and MungeSumstats harmonisation of one GWAS")
parser <- add_argument(parser, "--gwas", help = "GWAS summary statistics file", type = "character")
parser <- add_argument(parser, "--input_build", help = "Build of the file (GRCh37 or GRCh38)", default = "GRCh37")
parser <- add_argument(parser, "--output_build", help = "Build to harmonise to", default = "GRCh38")
parser <- add_argument(parser, "--columns", help = "GeneHackman column map (preset or key=value list)", default = "default")
parser <- add_argument(parser, "--N", help = "Sample size, when the file has none", type = "numeric", default = 0)
parser <- add_argument(parser, "--chrpos_only", help = "Rewrite the input to chr:pos:a1:a2 identifiers only (Reviewer 1's case)",
                       flag = TRUE)
parser <- add_argument(parser, "--output_prefix", help = "Prefix of the output files", type = "character")
parser <- add_argument(parser, "--threads", help = "Threads for MungeSumstats", type = "numeric", default = 4)
args <- parse_args(parser)

if (!requireNamespace("MungeSumstats", quietly = TRUE)) {
  stop("MungeSumstats is not installed; see scripts/validation/README.md")
}
dir.create(dirname(args$output_prefix), recursive = TRUE, showWarnings = FALSE)

# the same input for both tools, in GeneHackman's column names
# Only the mapped columns are kept: MungeSumstats guesses columns from their names, so extra ones (e.g. HERMES's logP,
# N_case, p_het) could be read as something else and make the two tools' inputs differ.
input <- data.table::as.data.table(change_column_names(get_file_or_dataframe(args$gwas), resolve_column_map(args$columns),
                                                       remove_extra_columns = TRUE))
if (isTRUE(args$chrpos_only)) {
  input <- data.table::as.data.table(as_chrpos_identifiers(input))
  message("Input rewritten to chr:pos:a1:a2 identifiers only, e.g. ", input$SNP[1])
}
input_file <- paste0(args$output_prefix, "_input.tsv.gz")
data.table::fwrite(input, input_file, sep = "\t")

# GeneHackman: full rsID lookup in dbSNP of the output build
start <- Sys.time()
genehackman <- standardise_gwas(input_file, output_file = paste0(args$output_prefix, "_genehackman_std.tsv.gz"),
                                N = args$N, populate_rsid_option = populate_rsid_options$full,
                                input_reference_build = args$input_build, output_reference_build = args$output_build)
message("GeneHackman: ", nrow(genehackman), " variants in ", format(Sys.time() - start))

# MungeSumstats: A2 is the effect allele, A1 the other allele
munge_input <- data.table::copy(input)
data.table::setnames(munge_input, c("EA", "OA", "EAF"), c("A2", "A1", "FRQ"), skip_absent = TRUE)
munge_file <- paste0(args$output_prefix, "_mungesumstats_input.tsv.gz")
data.table::fwrite(munge_input, munge_file, sep = "\t")
start <- Sys.time()
munged <- MungeSumstats::format_sumstats(
  munge_file, ref_genome = args$input_build,
  convert_ref_genome = if (args$output_build != args$input_build) args$output_build else NULL,
  return_data = TRUE, return_format = "data.table", nThread = args$threads,
  save_path = paste0(args$output_prefix, "_mungesumstats.tsv.gz"),
  log_folder_ind = TRUE, log_folder = paste0(args$output_prefix, "_mungesumstats_logs"))
munged <- data.table::as.data.table(if (is.list(munged) && !is.data.frame(munged)) munged$sumstats else munged)
message("MungeSumstats: ", nrow(munged), " variants in ", format(Sys.time() - start))

result <- compare_harmonisation(genehackman, munged)
summary <- rbind(data.frame(MEASURE = c("input_variants", "chrpos_only", "input_build", "output_build"),
                            VALUE = c(nrow(input), isTRUE(args$chrpos_only), args$input_build, args$output_build)),
                 result$summary)
data.table::fwrite(summary, paste0(args$output_prefix, "_summary.tsv"), sep = "\t")
data.table::fwrite(result$discordant, paste0(args$output_prefix, "_discordant.tsv.gz"), sep = "\t")
print(summary, row.names = FALSE)
