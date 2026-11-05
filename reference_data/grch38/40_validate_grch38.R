# Checks on the GRCh38 reference data, written to $WORK_DIR/validation_report.md:
#
#   Rscript 40_validate_grch38.R [--pop EUR] [--chrom 22]
#
# 1. Pan-UKB LD pack GRCh38 keys vs the natively GRCh38 1000G 30x panel: the share of pack SNP38 keys found in
#    the panel, the share of common panel variants (MAF >= 5%) found in the pack, and EAF agreement.  This is
#    the first check of the pack's lifted keys against native GRCh38 data.
# 2. GeneHackman's liftover vs native positions: the GRCh37 phase 3 panel lifted to GRCh38 with
#    convert_reference_build_via_liftover(), compared with the 30x panel's positions of the same rsIDs.
# 3. Panel continuity: r between the same rsID pairs in the GRCh37 phase 3 panel and the GRCh38 30x panel.
# 4. Files present in the bundle (panels, .frq, marker map, dbSNP, LD scores, chains).
# Needs PLINK 1.9 ($PLINK19) and GeneHackman (R_CMD); settings from config.sh.
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  args[i + 1]
}
pop <- get_arg("--pop", "EUR")
chrom <- get_arg("--chrom", "22")
bundle <- Sys.getenv("BUNDLE_DIR")
old_data <- Sys.getenv("OLD_PIPELINE_DATA_DIR")
pack_dir <- Sys.getenv("LD_PACK_DIR")
work_dir <- Sys.getenv("WORK_DIR")
plink <- Sys.getenv("PLINK19", "plink")
suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))
options(GeneHackman.liftover_chain_dir = file.path(bundle, "genomic_data", "liftover"))

report <- c("# GRCh38 reference data: validation", "", paste0("Population ", pop, ", chromosome ", chrom, "."), "")
add <- function(...) report <<- c(report, paste0(...))
panel38 <- file.path(bundle, "genomic_data", "1000genomes", "b38_dbsnp156", pop)
panel37 <- file.path(old_data, "genomic_data", "1000genomes", "b37_dbsnp156", pop)

read_panel <- function(prefix) {
  bim <- fread(paste0(prefix, ".bim"), header = FALSE, col.names = c("CHR", "RSID", "CM", "BP", "A1", "A2"),
               colClasses = list(character = c(1, 2, 5, 6)))
  frq <- fread(paste0(prefix, ".frq"), select = c("SNP", "A1", "MAF"), colClasses = list(character = c("SNP", "A1")))
  # .bim A1 is the alphabetical EA; .frq A1 is the minor allele (PLINK 1.9 --freq), so flip MAF where they differ
  i <- match(bim$RSID, frq$SNP)
  bim[, EAF := ifelse(frq$A1[i] == A1, frq$MAF[i], 1 - frq$MAF[i])]
  bim[, KEY := toupper(paste0(CHR, ":", BP, "_", pmin(A1, A2), "_", pmax(A1, A2)))]
  bim
}

# 1. pack GRCh38 keys vs native GRCh38 panel
p38 <- read_panel(panel38)[CHR == chrom]
pack <- fread(file.path(pack_dir, resolve_ld_ancestry(pop, "pan_ukb"), paste0("variants_chr", chrom, ".tsv.gz")),
              select = c("SNP38", "EAF", "LIFT_STATUS"), na.strings = c("", "NA"))
pack_ok <- pack[!is.na(SNP38)]
in_panel <- match(pack_ok$SNP38, p38$KEY)
common <- p38[pmin(EAF, 1 - EAF) >= 0.05]
eaf_diff <- pack_ok$EAF[!is.na(in_panel)] - p38$EAF[in_panel[!is.na(in_panel)]]
add("## 1. Pan-UKB GRCh38 keys vs the native GRCh38 1000G 30x panel", "")
add("- Pack variants with a GRCh38 key: ", nrow(pack_ok), " of ", nrow(pack), " (",
    paste(names(table(pack$LIFT_STATUS)), table(pack$LIFT_STATUS), collapse = ", "), ")")
add(sprintf("- Pack GRCh38 keys found in the panel: %.1f%%", 100 * mean(!is.na(in_panel))))
add(sprintf("- Panel variants with MAF >= 5%% found in the pack: %.1f%%", 100 * mean(common$KEY %in% pack_ok$SNP38)))
add(sprintf("- EAF agreement (shared variants): r = %.4f, mean |difference| = %.4f, %d with |difference| > 0.2",
            cor(pack_ok$EAF[!is.na(in_panel)], p38$EAF[in_panel[!is.na(in_panel)]]), mean(abs(eaf_diff)),
            sum(abs(eaf_diff) > 0.2)))
add("", "Expect most pack keys in the panel (UK Biobank is imputed from other references, so not all), EAF r > 0.95,",
    " and very few |difference| > 0.2 (which would mean swapped alleles).", "")

# 2. GeneHackman liftover vs native positions
if (dir.exists(dirname(panel37))) {
  p37 <- read_panel(panel37)[CHR == chrom]
  lifted <- suppressMessages(convert_reference_build_via_liftover(
    data.frame(CHR = p37$CHR, BP = p37$BP, EA = pmin(p37$A1, p37$A2), OA = pmax(p37$A1, p37$A2), RSID = p37$RSID),
    "GRCh37", "GRCh38", max_drop_fraction = 1))
  shared <- merge(as.data.table(lifted)[, .(RSID, BP_LIFTED = BP)], p38[, .(RSID, BP_NATIVE = BP)], by = "RSID")
  qc <- attr(lifted, "liftover_qc")
  add("## 2. GeneHackman liftover vs native GRCh38 positions", "")
  add("- GRCh37 panel variants lifted: ", nrow(lifted), " of ", nrow(p37), " (dropped: ",
      paste(qc$REASON[qc$N > 0 & !qc$REASON %in% c("input", "kept")], qc$N[qc$N > 0 & !qc$REASON %in% c("input", "kept")],
            collapse = ", "), ")")
  add(sprintf("- Same position as the 30x panel, among %d shared rsIDs: %.3f%% (%d differ)", nrow(shared),
              100 * mean(shared$BP_LIFTED == shared$BP_NATIVE), sum(shared$BP_LIFTED != shared$BP_NATIVE)))
  if (any(shared$BP_LIFTED != shared$BP_NATIVE)) {
    add("- Examples: ", paste(head(shared[BP_LIFTED != BP_NATIVE, paste(RSID, BP_LIFTED, BP_NATIVE)], 5), collapse = "; "))
  }
  add("", "Expect 100% (rsIDs whose positions differ were usually re-mapped or merged in dbSNP).", "")

  # 3. panel continuity: r for the same rsID pairs, for 3,000 shared variants in a 5 Mb window
  window <- range(p38$BP)[1] + c(5e6, 10e6)
  ids <- intersect(p37$RSID, p38[BP >= window[1] & BP <= window[2]]$RSID)
  set.seed(1)
  ids <- sample(ids, min(3000L, length(ids)))
  snps <- tempfile()
  writeLines(ids, snps)
  ld <- function(prefix) {
    out <- tempfile()
    system2(plink, c("--bfile", prefix, "--chr", chrom, "--extract", snps, "--r", "--ld-window-kb", "500",
                     "--ld-window", "99999", "--ld-window-r2", "0", "--keep-allele-order", "--out", out),
            stdout = FALSE, stderr = FALSE)
    fread(paste0(out, ".ld"))[, .(PAIR = paste(SNP_A, SNP_B), R)]
  }
  both <- merge(ld(panel37), ld(panel38), by = "PAIR", suffixes = c("37", "38"))
  add("## 3. Panel continuity (phase 3 GRCh37 vs 30x GRCh38, same rsIDs)", "")
  add(sprintf("- %d variant pairs within 500 kb: r correlation %.4f, mean |difference| %.4f, sign agreement for |r| > 0.3: %.2f%%",
              nrow(both), cor(both$R37, both$R38), mean(abs(both$R37 - both$R38)),
              100 * mean(sign(both$R37[abs(both$R37) > 0.3]) == sign(both$R38[abs(both$R37) > 0.3]))))
  add("", "Mostly the same people, re-sequenced: expect r correlation > 0.95 and sign agreement ~100%.", "")
} else {
  add("## 2-3. Skipped: GRCh37 panel not found at ", dirname(panel37), " (set OLD_PIPELINE_DATA_DIR)", "")
}

# 4. files
add("## 4. Bundle contents", "")
expect <- c(file.path("genomic_data/1000genomes/b38_dbsnp156", paste0(c("EUR", "EAS", "AFR", "AMR", "SAS"), ".bed")),
            file.path("genomic_data/1000genomes/b38_dbsnp156", paste0(c("EUR", "EAS", "AFR", "AMR", "SAS"), ".frq")),
            "genomic_data/1000genomes/b38_dbsnp156/marker_to_rsid.tsv.gz",
            "genomic_data/dbsnp/b38_dbsnp156/chr1.fst",
            file.path("genomic_data/liftover", c("hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz",
                                                 "hg18ToHg38.over.chain.gz", "hg18ToHg19.over.chain.gz",
                                                 "hg19ToHg18.over.chain.gz")),
            "genomic_data/pan_ukb_ld/MANIFEST.json",
            file.path("LDSCORE/pan_ukb", c("w_hm3.snplist", paste0(c("EUR", "AFR", "AMR", "CSA", "EAS", "MID"),
                                                                  "/ld.l2.ldscore.gz"))))
present <- file.exists(file.path(bundle, expect))
add(paste0("- ", ifelse(present, "OK      ", "MISSING "), expect))
writeLines(report, file.path(work_dir, "validation_report.md"))
cat(report, sep = "\n")
