# Concordance of GeneHackman's harmonisation with MungeSumstats (Murphy et al. 2021) on the same GWAS file.
#
# Run from the scripts/ directory, in an R with GeneHackman's dependencies and MungeSumstats (plus the
# SNPlocs/BSgenome packages MungeSumstats asks for), e.g.:
#   Rscript validation/compare_with_mungesumstats.R --gwas hf_all.tsv.gz --input_build GRCh37 \
#     --output_build GRCh38 --columns metal --N 1665481 --output_prefix /results/validation/hf_all
#
# Both tools standardise the file; variants are matched on chromosome, position and the unordered allele pair,
# and the script reports, for the matched variants: the share with the same rsID, the share with the same
# effect direction after aligning to a common effect allele, the correlation of the aligned effects, and the
# variants kept by only one tool.  Writes <output_prefix>_summary.tsv and <output_prefix>_discordant.tsv.gz.
source("load.R")
library(argparser, quietly = TRUE)

parser <- arg_parser("Compare GeneHackman and MungeSumstats harmonisation of one GWAS")
parser <- add_argument(parser, "--gwas", help = "GWAS summary statistics file", type = "character")
parser <- add_argument(parser, "--input_build", help = "Build of the file (GRCh37 or GRCh38)", default = "GRCh37")
parser <- add_argument(parser, "--output_build", help = "Build to harmonise to", default = "GRCh38")
parser <- add_argument(parser, "--columns", help = "GeneHackman column map (preset or key=value list)", default = "default")
parser <- add_argument(parser, "--N", help = "Sample size, when the file has none", type = "numeric", default = 0)
parser <- add_argument(parser, "--output_prefix", help = "Prefix of the output files", type = "character")
parser <- add_argument(parser, "--threads", help = "Threads for MungeSumstats", type = "numeric", default = 4)
args <- parse_args(parser)

if (!requireNamespace("MungeSumstats", quietly = TRUE)) {
  stop("MungeSumstats is not installed: BiocManager::install('MungeSumstats') (and the SNPlocs/BSgenome ",
       "packages for ", args$input_build, " and ", args$output_build, " it asks for)")
}
dir.create(dirname(args$output_prefix), recursive = TRUE, showWarnings = FALSE)

# GeneHackman: full rsID lookup in dbSNP so both tools assign rsIDs
genehackman <- standardise_gwas(args$gwas, output_file = paste0(args$output_prefix, "_genehackman_std.tsv.gz"),
                                N = args$N, populate_rsid_option = populate_rsid_options$full,
                                input_reference_build = args$input_build, output_reference_build = args$output_build,
                                input_columns = args$columns)
genehackman <- data.table::as.data.table(genehackman)

# MungeSumstats reads the original file; rename the columns GeneHackman was told about to names it recognises
raw <- data.table::as.data.table(change_column_names(get_file_or_dataframe(args$gwas), resolve_column_map(args$columns)))
munge_input <- tempfile(fileext = ".tsv.gz")
data.table::fwrite(raw[, intersect(c("SNP", "RSID", "CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"),
                                   names(raw)), with = FALSE] |>
                     data.table::setnames(c("EA", "OA", "EAF"), c("A2", "A1", "FRQ"), skip_absent = TRUE),
                   munge_input, sep = "\t")
munged <- MungeSumstats::format_sumstats(munge_input, ref_genome = args$input_build,
                                         convert_ref_genome = if (args$output_build != args$input_build) args$output_build,
                                         return_data = TRUE, return_format = "data.table", nThread = args$threads,
                                         save_path = paste0(args$output_prefix, "_mungesumstats.tsv.gz"),
                                         log_folder_ind = TRUE, log_folder = dirname(args$output_prefix))
munged <- data.table::as.data.table(if (is.list(munged) && !is.data.frame(munged)) munged$sumstats else munged)

# MungeSumstats: A2 is the effect allele.  Match on CHR, BP and the unordered allele pair.
pair_key <- function(chr, bp, a, b) paste0(chr, ":", bp, "_", pmin(toupper(a), toupper(b)), "_", pmax(toupper(a), toupper(b)))
genehackman[, KEY := pair_key(CHR, BP, EA, OA)]
munged[, KEY := pair_key(CHR, BP, A1, A2)]
matched <- merge(genehackman[, .(KEY, RSID_GH = RSID, EA_GH = toupper(EA), BETA_GH = BETA)],
                 munged[, .(KEY, RSID_MSS = SNP, EA_MSS = toupper(A2), BETA_MSS = BETA)], by = "KEY")
matched[, BETA_MSS_ALIGNED := ifelse(EA_MSS == EA_GH, BETA_MSS, -BETA_MSS)]

summary <- data.table::data.table(
  MEASURE = c("variants_genehackman", "variants_mungesumstats", "variants_matched", "only_genehackman",
              "only_mungesumstats", "rsid_identical_pct", "effect_direction_concordant_pct", "aligned_beta_correlation"),
  VALUE = c(nrow(genehackman), nrow(munged), nrow(matched), sum(!genehackman$KEY %in% munged$KEY),
            sum(!munged$KEY %in% genehackman$KEY),
            100 * mean(matched$RSID_GH == matched$RSID_MSS, na.rm = TRUE),
            100 * mean(sign(matched$BETA_GH) == sign(matched$BETA_MSS_ALIGNED), na.rm = TRUE),
            stats::cor(matched$BETA_GH, matched$BETA_MSS_ALIGNED, use = "complete.obs")))
data.table::fwrite(summary, paste0(args$output_prefix, "_summary.tsv"), sep = "\t")
discordant <- matched[RSID_GH != RSID_MSS | sign(BETA_GH) != sign(BETA_MSS_ALIGNED)]
data.table::fwrite(discordant, paste0(args$output_prefix, "_discordant.tsv.gz"), sep = "\t")
print(summary)
