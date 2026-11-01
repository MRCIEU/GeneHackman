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
# GRCh38 positions come from GeneHackman's convert_reference_build_via_liftover(), i.e. the same
# liftOver chains, coordinate convention and QC that lift a GWAS.  Variants it drops (unmapped,
# non_primary_contig, chrom_change, strand_flip, round_trip, or colliding on a GRCh38 key:
# duplicate_grch38) are kept as rows but get BP38/SNP38 = NA and that LIFT_STATUS.  Variants that
# collide on a GRCh37 key (LIFT_STATUS duplicate_grch37) get neither key.  The round-trip check needs
# hg38ToHg19.over.chain.gz next to hg19ToHg38.over.chain.gz (01_setup_env.sh fetches both).

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
for (needed in c(chain_name, "hg38ToHg19.over.chain.gz")) {
  if (!file.exists(file.path(chain_dir, needed))) {
    stop("Chain file missing: ", file.path(chain_dir, needed), " (run 01_setup_env.sh, which fetches both chains)")
  }
}
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
# Two source variants can share a GeneHackman key: e.g. a deletion AT>A and an insertion A>AT at one
# position both become CHR:BP_A_AT, because the key does not record which allele is the reference.
# A GWAS row with that key could be either, so neither gets a key; both stay in the matrix.
duplicated37 <- duplicated(annot$SNP37) | duplicated(annot$SNP37, fromLast = TRUE)
if (any(duplicated37)) {
  examples <- annot[duplicated37][order(SNP37)][seq_len(min(6L, sum(duplicated37)))]
  message(pop, ": ", sum(duplicated37), " variants share a GRCh37 key with another variant and get no keys, e.g. ",
          paste0(examples$SNP37, " (", examples$REF37, ">", examples$ALT37, ")", collapse = ", "))
}

# --- GRCh38 positions, with GeneHackman's own liftover (the same function and QC that lift a GWAS)
# convert_reference_build_via_liftover() keeps a variant only if it maps to the same chromosome on the
# plus strand and lifts back to where it started; everything else is dropped with a reason, which
# becomes LIFT_STATUS here (the variant keeps its row and LD, but gets no GRCh38 key).  So a GWAS lifted
# by GeneHackman and the pack agree on which variants have GRCh38 keys.
lift_in <- data.frame(CHR = annot$CHR, BP = annot$BP37, EA = annot$EA, OA = annot$OA, ROW_ID = annot$ROW_ID)
report_prefix <- file.path(raw_dir, "liftover")
lifted <- suppressMessages(convert_reference_build_via_liftover(
  lift_in, reference_builds$GRCh37, reference_builds$GRCh38, max_drop_fraction = 1, report_prefix = report_prefix
))
dropped <- data.table::fread(paste0(report_prefix, ".liftover_dropped.tsv.gz"), select = c("ROW_ID", "LIFTOVER_DROP_REASON"))
annot[, LIFT_STATUS := "ok"]
annot[match(dropped$ROW_ID, ROW_ID), LIFT_STATUS := dropped$LIFTOVER_DROP_REASON]
annot[LIFT_STATUS == "duplicate_target", LIFT_STATUS := "duplicate_grch38"]
annot[, BP38 := NA_real_]
annot[match(lifted$ROW_ID, ROW_ID), BP38 := as.numeric(lifted$BP)]
stopifnot(all(is.na(annot$BP38) == (annot$LIFT_STATUS != "ok")))
annot[duplicated37, LIFT_STATUS := "duplicate_grch37"]
annot[LIFT_STATUS != "ok", BP38 := NA_real_]

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

annot[LIFT_STATUS == "duplicate_grch37", SNP37 := NA_character_]

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
data.table::fwrite(lift_counts, summary_file, sep = "\t")
print(lift_counts)
message(pop, ": wrote ", length(unique(annot$CHR)), " annotated chromosome files to ", raw_dir)
