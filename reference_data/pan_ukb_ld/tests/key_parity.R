# Key parity: the pack's SNP37/SNP38 keys must equal the SNP column GeneHackman writes when it
# standardises a GWAS, in both builds.
#
#   WORK_DIR=<tmp> GENEHACKMAN_DIR=<repo> CHAIN_FILE=<hg19ToHg38.over.chain.gz> LIFTOVER_BIN=<liftOver> \
#     Rscript reference_data/pan_ukb_ld/tests/key_parity.R
#
# Takes the variants of tests/testthat/data/test_data_small.tsv.gz as if they were Pan-UKB variants (REF/ALT
# assigned at random, so both orientations occur), runs 20_annotate_variants.R on them, then standardises
# the GWAS itself with standardise_gwas() to GRCh37 and to GRCh38 and compares keys and EAF.
suppressPackageStartupMessages(library(data.table))
work_dir <- Sys.getenv("WORK_DIR")
gh_dir <- Sys.getenv("GENEHACKMAN_DIR")
here <- file.path(gh_dir, "reference_data", "pan_ukb_ld")
raw_dir <- file.path(work_dir, "raw", "TST")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

gwas <- as.data.frame(fread(file.path(gh_dir, "tests/testthat/data/test_data_small.tsv.gz")))
set.seed(3)
alt_is_ea <- runif(nrow(gwas)) < 0.5
raw <- data.table(CHR = as.character(gwas$CHR), BP37 = gwas$BP,
                  REF = ifelse(alt_is_ea, gwas$OA, gwas$EA), ALT = ifelse(alt_is_ea, gwas$EA, gwas$OA),
                  AF = ifelse(alt_is_ea, gwas$EAF, 1 - gwas$EAF), RSID_UKB = gwas$RSID)
raw <- raw[order(as.integer(CHR), BP37)]
raw[, IDX_GLOBAL := seq_len(.N) - 1]
fwrite(raw, file.path(raw_dir, "variants_raw.tsv.gz"), sep = "\t")

status <- system2(file.path(R.home("bin"), "Rscript"), c(file.path(here, "20_annotate_variants.R"), "--pop", "TST"))
if (status != 0) stop("20_annotate_variants.R failed")
annot <- rbindlist(lapply(list.files(raw_dir, "^variants_annot_chr", full.names = TRUE), fread, na.strings = "NA"))

Sys.setenv(PIPELINE_DATA_DIR = file.path(work_dir, "liftover_data"))
suppressMessages(pkgload::load_all(gh_dir, quiet = TRUE))
std37 <- as.data.table(suppressMessages(standardise_gwas(gwas, output_reference_build = "GRCh37")))
std38 <- as.data.table(suppressMessages(standardise_gwas(gwas, output_reference_build = "GRCh38")))
hg38_fixture <- fread(file.path(gh_dir, "tests/testthat/data/test_data_small_hg38.tsv.gz"))

report <- function(label, ok) {
  cat(sprintf("%-4s %s\n", if (ok) "OK" else "FAIL", label))
  ok
}
ok <- c(
  report(sprintf("GRCh37 keys: %d of %d standardised GWAS keys are pack keys", sum(std37$SNP %in% annot$SNP37),
                 nrow(std37)), all(std37$SNP %in% annot$SNP37)),
  report("GRCh37 EAF of the effect allele matches", isTRUE(all.equal(annot$EAF[match(std37$SNP, annot$SNP37)],
                                                                      std37$EAF, tolerance = 1e-6))),
  report(sprintf("GRCh38 keys: %d of %d standardised (lifted) GWAS keys are pack keys; %d pack variants not ok (%s)",
                 sum(std38$SNP %in% annot$SNP38), nrow(std38), sum(annot$LIFT_STATUS != "ok"),
                 paste(names(table(annot$LIFT_STATUS)), table(annot$LIFT_STATUS), collapse = ", ")),
         all(std38$SNP[std38$SNP %in% annot$SNP38] %in% annot$SNP38) &&
           mean(std38$SNP %in% annot$SNP38) > 0.999),
  report(sprintf("GRCh38 keys agree with the committed hg38 test GWAS: %d of %d",
                 sum(hg38_fixture$SNP %in% annot$SNP38), nrow(hg38_fixture)),
         mean(hg38_fixture$SNP %in% annot$SNP38) > 0.999)
)
missing38 <- setdiff(std38$SNP, annot$SNP38)
if (length(missing38)) cat("GRCh38 GWAS keys not in the pack:", head(missing38, 5), "\n")
extra <- annot[LIFT_STATUS == "ok" & !SNP38 %in% std38$SNP]
if (nrow(extra)) cat(nrow(extra), "pack variants have a GRCh38 key that standardise_gwas() did not produce, e.g.",
                     head(paste(extra$SNP37, "->", extra$SNP38), 3), "\n")
if (!all(ok)) quit(status = 1)
