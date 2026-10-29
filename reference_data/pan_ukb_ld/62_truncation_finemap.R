# Phase 0: does storing LD only within RADIUS_KB change fine-mapping?
#
#   Rscript 62_truncation_finemap.R --region phase0/EUR_region --gwas <standardised GWAS> [--radius-kb 500]
#
# Runs SuSiE-RSS in the phase 0 test region twice, with the full dense LD and with LD set to 0 beyond
# the radius (plus the same positive semi-definite fix the pipeline applies), and reports how much the
# PIPs and credible sets differ.  The GWAS must be GRCh37 and standardised by GeneHackman (SNP keys).

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  args[i + 1]
}
region_stem <- get_arg("--region")
gwas_file <- get_arg("--gwas")
radius_bp <- as.numeric(get_arg("--radius-kb", Sys.getenv("RADIUS_KB", "1000"))) * 1000
if (is.null(region_stem) || is.null(gwas_file)) stop("usage: --region <stem> --gwas <file> [--radius-kb 500]")
work_dir <- Sys.getenv("WORK_DIR", ".")
if (!startsWith(region_stem, "/")) region_stem <- file.path(work_dir, region_stem)

suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))

variants <- data.table::fread(paste0(region_stem, ".variants.tsv"))
n <- nrow(variants)
R_full <- matrix(readBin(paste0(region_stem, ".f32"), what = "numeric", size = 4, n = n * n, endian = "little"), n, n)
gwas <- data.table::fread(gwas_file)
gwas <- gwas[SNP %in% variants$SNP37 & !is.na(BETA) & !is.na(SE) & SE > 0]
idx <- match(gwas$SNP, variants$SNP37)
z <- gwas$BETA / gwas$SE
n_gwas <- if ("N" %in% names(gwas)) stats::median(gwas$N, na.rm = TRUE) else stop("GWAS needs an N column")
R_full <- R_full[idx, idx]
distance <- abs(outer(variants$BP37[idx], variants$BP37[idx], "-"))
R_trunc <- R_full
R_trunc[distance > radius_bp] <- 0
psd <- make_ld_positive_semidefinite(R_trunc)
message(length(z), " GWAS variants in the region; ridge added after truncation: ", signif(psd$ridge, 3))

fit <- function(R) susieR::susie_rss(z = z, R = R, n = n_gwas, L = 10, coverage = 0.95, min_abs_corr = 0.5)
full <- fit(make_ld_positive_semidefinite(R_full)$matrix)
trunc <- fit(psd$matrix)
cs_sets <- function(f) vapply(f$sets$cs, function(x) paste(sort(gwas$SNP[x]), collapse = ","), character(1))
result <- list(
  n_variants = length(z), radius_kb = radius_bp / 1000, ridge = psd$ridge,
  pip_correlation = stats::cor(full$pip, trunc$pip), max_abs_pip_difference = max(abs(full$pip - trunc$pip)),
  n_cs_full = length(full$sets$cs), n_cs_truncated = length(trunc$sets$cs),
  identical_credible_sets = setequal(cs_sets(full), cs_sets(trunc))
)
print(result)
out <- paste0(region_stem, "_truncation_finemap.tsv")
data.table::fwrite(data.table::as.data.table(result), out, sep = "\t")
message("wrote ", out)
