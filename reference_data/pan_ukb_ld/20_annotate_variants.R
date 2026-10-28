# Stage 20: give the kept variants GeneHackman keys in GRCh37 and GRCh38.
#
#   Rscript 20_annotate_variants.R --pop EUR
#
# Reads raw/<POP>/variants_raw.tsv.gz (CHR, BP37, REF, ALT, AF, IDX_GLOBAL, RSID_UKB) and writes
# raw/<POP>/variants_annot_chr<C>.tsv.gz with EA, OA, EAF, SIGN, SNP37, BP38, SNP38 and LIFT_STATUS.
#
# Keys are built with GeneHackman's own standardise_alleles(), so they match a standardised GWAS
# character for character: EA is the alphabetically smaller allele, EAF its frequency, and
# SIGN = +1 when EA is the ALT allele the source LD is counted on (-1 otherwise).
#
# GRCh38 positions come from UCSC liftOver with the same chain file and the same coordinate
# convention as convert_reference_build_via_liftover() (BED start = BP, new BP = mapped start);
# the result is cross-checked against that function on a sample of variants.  Variants that do not
# map, map to the minus strand (alleles would be complemented) or to another chromosome, or that
# collide on a GRCh38 key are kept as rows but get BP38/SNP38 = NA and a LIFT_STATUS.

args <- commandArgs(trailingOnly = TRUE)
pop <- args[which(args == "--pop") + 1]
if (length(pop) != 1 || is.na(pop)) stop("usage: Rscript 20_annotate_variants.R --pop EUR")

work_dir <- Sys.getenv("WORK_DIR")
gh_dir <- Sys.getenv("GENEHACKMAN_DIR")
if (!nzchar(work_dir) || !nzchar(gh_dir)) stop("WORK_DIR and GENEHACKMAN_DIR must be set (see config.sh)")

# GeneHackman finds the chain under $PIPELINE_DATA_DIR/genomic_data/liftover and liftOver on PATH.
liftover_data <- file.path(work_dir, "liftover_data")
chain_dir <- file.path(liftover_data, "genomic_data", "liftover")
chain_name <- "hg19ToHg38.over.chain.gz"
dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)
chain_file <- Sys.getenv("CHAIN_FILE")
target_chain <- normalizePath(file.path(chain_dir, chain_name), mustWork = FALSE)
if (nzchar(chain_file) && normalizePath(chain_file) != target_chain) {
  file.copy(chain_file, file.path(chain_dir, chain_name), overwrite = TRUE)
}
if (!file.exists(file.path(chain_dir, chain_name))) stop("Chain file missing: ", file.path(chain_dir, chain_name))
Sys.setenv(PIPELINE_DATA_DIR = liftover_data)
liftover_bin <- Sys.getenv("LIFTOVER_BIN")
if (nzchar(liftover_bin)) Sys.setenv(PATH = paste(dirname(liftover_bin), Sys.getenv("PATH"), sep = ":"))
if (!nzchar(Sys.which("liftOver")) && file.exists(file.path(work_dir, "bin", "liftOver"))) {
  Sys.setenv(PATH = paste(file.path(work_dir, "bin"), Sys.getenv("PATH"), sep = ":"))
}

suppressMessages(pkgload::load_all(gh_dir, quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))
if (!file.exists(liftover_binary)) stop("liftOver not found (set LIFTOVER_BIN or run 01_setup_env.sh)")
message("Using GeneHackman at ", gh_dir, "; liftOver ", liftover_binary, "; chain ",
        file.path(liftover_chain_dir, chain_name))

raw_dir <- file.path(work_dir, "raw", pop)
raw <- data.table::fread(file.path(raw_dir, "variants_raw.tsv.gz"),
                         colClasses = list(character = c("CHR", "REF", "ALT", "RSID_UKB")))
raw[, ROW_ID := .I]
message(pop, ": annotating ", nrow(raw), " variants")

# --- GRCh37 keys and orientation, via GeneHackman's allele standardisation
pseudo <- data.frame(CHR = raw$CHR, BP = raw$BP37, EA = toupper(raw$ALT), OA = toupper(raw$REF),
                     BETA = 1, EAF = raw$AF, ROW_ID = raw$ROW_ID)
std37 <- data.table::as.data.table(standardise_alleles(pseudo))
data.table::setorder(std37, ROW_ID)
stopifnot(identical(std37$ROW_ID, raw$ROW_ID))
annot <- data.table::data.table(
  CHR = raw$CHR, BP37 = raw$BP37, EA = std37$EA, OA = std37$OA, EAF = round(std37$EAF, 6),
  SIGN = as.integer(std37$BETA), SNP37 = std37$SNP, REF37 = raw$REF, ALT37 = raw$ALT,
  RSID_UKB = raw$RSID_UKB, IDX_GLOBAL = raw$IDX_GLOBAL, ROW_ID = raw$ROW_ID
)
stopifnot(all(annot$SIGN %in% c(-1L, 1L)))
stopifnot(all(ifelse(annot$SIGN == 1L, annot$EA == toupper(annot$ALT37), annot$EA == toupper(annot$REF37))))
if (anyDuplicated(annot$SNP37)) stop("Duplicate GRCh37 keys in the source variant list")

# --- GRCh38 positions: one strand-aware liftOver run over every variant
bed_in <- tempfile(fileext = ".bed")
bed_out <- tempfile(fileext = ".bed")
unmapped <- tempfile(fileext = ".unmapped")
data.table::fwrite(
  data.table::data.table(paste0("chr", annot$CHR), format(annot$BP37, scientific = FALSE, trim = TRUE),
                         format(annot$BP37 + 1, scientific = FALSE, trim = TRUE), annot$ROW_ID, 0, "+"),
  bed_in, sep = "\t", col.names = FALSE, quote = FALSE
)
status <- system2(liftover_binary, c(bed_in, file.path(liftover_chain_dir, chain_name), bed_out, unmapped))
if (status != 0) stop("liftOver failed with exit code ", status)
lifted <- data.table::fread(bed_out, header = FALSE, colClasses = list(character = c(1, 6)),
                            col.names = c("CHROM38", "START38", "END38", "ROW_ID", "SCORE", "STRAND38"))
annot <- merge(annot, lifted[, .(ROW_ID, CHROM38 = sub("^chr", "", CHROM38), START38, STRAND38)],
               by = "ROW_ID", all.x = TRUE)
data.table::setorder(annot, ROW_ID)

annot[, LIFT_STATUS := data.table::fcase(
  is.na(START38), "unmapped",
  CHROM38 != CHR, "chrom_change",
  STRAND38 != "+", "strand_flip",
  default = "ok"
)]
annot[, BP38 := data.table::fifelse(LIFT_STATUS == "ok", as.numeric(START38), NA_real_)]

# --- GRCh38 keys, again via standardise_alleles (alleles are already canonical, so this only formats)
ok <- annot$LIFT_STATUS == "ok"
pseudo38 <- data.frame(CHR = annot$CHR[ok], BP = annot$BP38[ok], EA = annot$EA[ok], OA = annot$OA[ok],
                       BETA = 1, ROW_ID = annot$ROW_ID[ok])
std38 <- data.table::as.data.table(standardise_alleles(pseudo38))
stopifnot(all(std38$BETA == 1))
annot[, SNP38 := NA_character_]
annot[match(std38$ROW_ID, ROW_ID), SNP38 := std38$SNP]
duplicated38 <- !is.na(annot$SNP38) & (duplicated(annot$SNP38) | duplicated(annot$SNP38, fromLast = TRUE))
annot[duplicated38, `:=`(LIFT_STATUS = "duplicate_grch38", BP38 = NA_real_, SNP38 = NA_character_)]

# --- cross-check against GeneHackman's own GWAS liftover on a sample
set.seed(20261028)
sample_rows <- sort(sample(which(annot$LIFT_STATUS == "ok"), min(200000L, sum(annot$LIFT_STATUS == "ok"))))
check_in <- data.frame(CHR = annot$CHR[sample_rows], BP = annot$BP37[sample_rows], ROW_ID = annot$ROW_ID[sample_rows])
check_out <- suppressMessages(convert_reference_build_via_liftover(check_in, reference_builds$GRCh37,
                                                                   reference_builds$GRCh38))
check <- merge(data.table::as.data.table(check_out)[, .(ROW_ID, BP_PACKAGE = BP)], annot[, .(ROW_ID, BP38)],
               by = "ROW_ID")
n_disagree <- sum(check$BP_PACKAGE != check$BP38)
n_package_dropped <- length(sample_rows) - nrow(check)
message(pop, ": liftover cross-check on ", length(sample_rows), " variants: ", n_disagree, " disagree, ",
        n_package_dropped, " not returned by convert_reference_build_via_liftover()")
if (n_disagree > 0) stop("GRCh38 positions disagree with convert_reference_build_via_liftover(); not writing keys")

# --- write per chromosome
columns <- c("CHR", "BP37", "BP38", "EA", "OA", "EAF", "SNP37", "SNP38", "LIFT_STATUS", "REF37", "ALT37",
             "RSID_UKB", "SIGN", "IDX_GLOBAL")
for (chrom in unique(annot$CHR)) {
  out_file <- file.path(raw_dir, paste0("variants_annot_chr", chrom, ".tsv.gz"))
  data.table::fwrite(annot[CHR == chrom, ..columns], paste0(out_file, ".tmp.gz"), sep = "\t", na = "NA")
  file.rename(paste0(out_file, ".tmp.gz"), out_file)
}
lift_counts <- annot[, .N, by = LIFT_STATUS]
summary_file <- file.path(raw_dir, "annotate_summary.tsv")
data.table::fwrite(rbind(lift_counts,
                         data.table::data.table(LIFT_STATUS = c("crosscheck_sampled", "crosscheck_disagree",
                                                                "crosscheck_dropped_by_package"),
                                                N = c(length(sample_rows), n_disagree, n_package_dropped))),
                   summary_file, sep = "\t")
print(lift_counts)
message(pop, ": wrote ", length(unique(annot$CHR)), " annotated chromosome files to ", raw_dir)
