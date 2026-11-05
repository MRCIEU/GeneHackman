# Build GeneHackman's QTL data (cis QTLs, GRCh38) from the genotype-phenotype map's standardised study data:
#
#   Rscript 31_qtl_from_gpm.R [--datasets eqtlgen,brainmeta,ukb_pqtl] [--cores 16] [--limit N]
#
# Datasets (all EUR, GRCh38 in the genotype-phenotype map):
#   eqtlgen    <- source "eqtlgen", studies eqtlgen-cis-<ENSG>      (blood cis-eQTLs, Vosa 2021)
#   brainmeta  <- source "brain_meta", studies BrainMeta-cis-eQTL-*  (brain cis-eQTLs, Qi 2022)
#   ukb_pqtl   <- source "ukb_ppp", genome-wide pQTLs per protein    (plasma pQTLs, Sun 2023): cis leads only
#
# For each study (exposure):
#   1. cis leads: rows of extracted_snps.tsv with cis_trans == "cis".  For ukb_pqtl the lead must also lie within
#      1 Mb of the gene (gene_ranges_hg38.tsv); leads outside are dropped and counted.
#   2. gwas/cis/<EXPOSURE>.tsv.gz: the standardised (pre-imputation) LD-block regions of those leads, joined, one
#      row per SNP, with N (the study's sample size).  This is what coloc uses.
#   3. top_hits/cis.tsv.gz: each cis lead with P < 5e-8 (one per LD block) is an MR instrument.
# EXPOSURE is the gene symbol, or <symbol>_<ENSG / OID> when the symbol is not unique in the dataset.
#
# Output: NEW_QTL_DATA_DIR/<dataset>/GRCh38/{top_hits/cis.tsv.gz, gwas/cis/<EXPOSURE>.tsv.gz, exposures.tsv,
# summary.tsv}.  Per-study results are cached in WORK_DIR/qtl_work, so a re-run only does missing studies.
# Settings from the environment (config.sh): NEW_QTL_DATA_DIR, WORK_DIR, DBSNP_B38_DIR, GENEHACKMAN_DIR, and
# GPM_STUDIES / GPM_GENE_RANGES (defaults below).
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  args[i + 1]
}
datasets <- strsplit(get_arg("--datasets", "eqtlgen,brainmeta,ukb_pqtl"), ",")[[1]]
cores <- as.integer(get_arg("--cores", "8"))
limit <- as.integer(get_arg("--limit", "0"))
p_instrument <- 5e-8
cis_window <- 1e6

gpm_studies <- Sys.getenv("GPM_STUDIES",
                          "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/studies_processed.tsv.gz")
gpm_gene_ranges <- Sys.getenv("GPM_GENE_RANGES",
                              "/local-scratch/projects/genotype-phenotype-map/data/gene_ranges_hg38.tsv")
new_dir <- Sys.getenv("NEW_QTL_DATA_DIR")
work_dir <- file.path(Sys.getenv("WORK_DIR"), "qtl_work")
dbsnp_dir <- Sys.getenv("DBSNP_B38_DIR")
if (!nzchar(new_dir) || !nzchar(Sys.getenv("WORK_DIR"))) stop("set NEW_QTL_DATA_DIR and WORK_DIR (see config.sh)")

suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))

dataset_sources <- list(
  eqtlgen = function(s) s[source == "eqtlgen" & startsWith(study_name, "eqtlgen-cis-")],
  brainmeta = function(s) s[source == "brain_meta"],
  ukb_pqtl = function(s) s[source == "ukb_ppp"]
)
unknown <- setdiff(datasets, names(dataset_sources))
if (length(unknown) > 0) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

region_columns <- c("SNP", "CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P")
study_top_hit_columns <- c("SNP", "CHR", "BP", "EA", "OA", "EAF", "BETA", "SE", "P", "N", "id", "ENSEMBL_ID", "GENE")
top_hit_columns <- c(study_top_hit_columns, "EXPOSURE")

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# EXPOSURE: the gene symbol where it is unique in the dataset, else <symbol>_<identifier>
exposure_names <- function(studies, dataset) {
  symbol <- studies$gene
  identifier <- switch(dataset,
    ukb_pqtl = sub("^.*:(OID[0-9]+):.*$", "\\1", studies$probe),
    ifelse(is.na(studies$ensg), sub("\\..*$", "", studies$probe), studies$ensg)
  )
  symbol[is.na(symbol)] <- identifier[is.na(symbol)]
  duplicated_symbol <- symbol %in% symbol[duplicated(symbol)]
  name <- ifelse(duplicated_symbol, paste0(symbol, "_", identifier), symbol)
  if (anyDuplicated(name)) name <- make.unique(name, sep = "_")
  safe_name(name)
}

gene_ranges <- unique(fread(gpm_gene_ranges)[ancestry == "EUR", .(chr = as.character(chr), start, stop, gene_name)])

process_study <- function(study, dataset, out_gwas_dir, cache_dir) {
  cache_file <- file.path(cache_dir, paste0(study$EXPOSURE, ".tsv"))
  if (file.exists(cache_file)) {
    return(fread(cache_file, colClasses = list(character = c("CHR", "SNP")), na.strings = c("", "NA")))
  }
  status <- data.table(EXPOSURE = study$EXPOSURE, STATUS = NA_character_, CIS_LEADS = 0L, OUTSIDE_WINDOW = 0L,
                       INSTRUMENTS = 0L, REGION_SNPS = 0L)
  result <- function(top_hits = NULL) {
    out <- if (is.null(top_hits) || nrow(top_hits) == 0) cbind(status, empty_top_hits) else cbind(status, top_hits)
    fwrite(out, paste0(cache_file, ".tmp"), sep = "\t")
    file.rename(paste0(cache_file, ".tmp"), cache_file)
    out
  }
  study_dir <- sub("/+$", "", study$extracted_location)
  leads_file <- file.path(study_dir, "extracted_snps.tsv")
  if (!file.exists(leads_file)) { status$STATUS <- "no_extracted_snps"; return(result()) }
  leads <- fread(leads_file, colClasses = list(character = "chr"))
  if (!"cis_trans" %in% names(leads)) { status$STATUS <- "no_cis_lead"; return(result()) }  # empty file: no leads
  leads <- leads[cis_trans == "cis"]
  if (dataset == "ukb_pqtl" && nrow(leads) > 0) {
    gene <- gene_ranges[gene_name == study$gene]
    if (nrow(gene) > 0) {
      inside <- vapply(seq_len(nrow(leads)), function(i) {
        any(gene$chr == leads$chr[i] & leads$bp[i] >= gene$start - cis_window & leads$bp[i] <= gene$stop + cis_window)
      }, logical(1))
      status$OUTSIDE_WINDOW <- sum(!inside)
      leads <- leads[inside]
    }
  }
  status$CIS_LEADS <- nrow(leads)
  if (nrow(leads) == 0) { status$STATUS <- "no_cis_lead"; return(result()) }

  region_files <- file.path(study_dir, "standardised", basename(leads$file))
  missing <- !file.exists(region_files)
  if (all(missing)) { status$STATUS <- "no_region_files"; return(result()) }
  leads <- leads[!missing]
  regions <- lapply(region_files[!missing], fread, select = region_columns,
                    colClasses = list(character = c("CHR", "SNP", "EA", "OA")))
  region <- unique(rbindlist(regions), by = "SNP")
  setorder(region, CHR, BP)
  region[, N := study$sample_size]
  status$REGION_SNPS <- nrow(region)
  gwas_file <- file.path(out_gwas_dir, paste0(study$EXPOSURE, ".tsv.gz"))
  fwrite(region, paste0(gwas_file, ".tmp.gz"), sep = "\t")
  file.rename(paste0(gwas_file, ".tmp.gz"), gwas_file)

  # instruments: the lead variant of each cis lead with P < 5e-8 (the strongest row at the lead position);
  # regions[[k]] is the region file of leads[k]
  strong <- which(leads$log_p > -log10(p_instrument))
  top_hits <- rbindlist(lapply(strong, function(k) {
    rows <- regions[[k]][CHR == leads$chr[k] & BP == leads$bp[k]]
    if (nrow(rows) == 0) return(NULL)
    rows[order(P)][1]
  }))
  if (nrow(top_hits) > 0) {
    top_hits <- unique(top_hits, by = "SNP")
    top_hits[, `:=`(N = study$sample_size, id = study$study_name,
                    ENSEMBL_ID = if (is.na(study$ensg)) NA_character_ else study$ensg, GENE = study$gene)]
    top_hits <- top_hits[, ..study_top_hit_columns]
  }
  status$INSTRUMENTS <- nrow(top_hits)
  status$STATUS <- "ok"
  result(top_hits)
}

# RSIDs of the instruments from dbSNP (DBSNP_B38_DIR), one chromosome at a time with GeneHackman's own matcher
# (chrpos_to_rsid() itself runs its workers through future, which needs GeneHackman installed, not load_all()ed)
add_rsids <- function(top_hits) {
  per_chrom <- split(copy(top_hits)[, CHR := as.character(CHR)], by = "CHR")
  out <- rbindlist(lapply(per_chrom, function(d) {
    suppressMessages(process_chromosome(d, "CHR", "BP", build = basename(dbsnp_dir), dbsnp_dir = dirname(dbsnp_dir),
                                        flip = "allow", alt_rsids = FALSE, p = function(...) NULL,
                                        nea_col = "OA", ea_col = "EA"))
  }), fill = TRUE)
  stopifnot(nrow(out) == nrow(top_hits))
  out[, c(top_hit_columns[1:10], "RSID", top_hit_columns[11:14]), with = FALSE]
}

# one all-NA row: a study without instruments still leaves its status row in the cache
empty_top_hits <- data.table(SNP = NA_character_, CHR = NA_character_, BP = NA_integer_, EA = NA_character_,
                             OA = NA_character_, EAF = NA_real_, BETA = NA_real_, SE = NA_real_, P = NA_real_,
                             N = NA_real_, id = NA_character_, ENSEMBL_ID = NA_character_, GENE = NA_character_)

all_studies <- fread(gpm_studies, na.strings = c("", "NA"))
for (dataset in datasets) {
  studies <- dataset_sources[[dataset]](all_studies)
  studies <- studies[order(study_name)]
  studies[, EXPOSURE := exposure_names(studies, dataset)]
  if (limit > 0) studies <- head(studies, limit)
  target <- file.path(new_dir, dataset, reference_builds$GRCh38)
  out_gwas_dir <- file.path(target, "gwas", "cis")
  cache_dir <- file.path(work_dir, dataset)
  dir.create(out_gwas_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(target, "top_hits"), recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  message(dataset, ": ", nrow(studies), " studies, ", cores, " cores")

  results <- parallel::mclapply(seq_len(nrow(studies)), function(i) {
    tryCatch(process_study(studies[i], dataset, out_gwas_dir, cache_dir),
             error = function(e) structure(conditionMessage(e), class = "study_error", study = studies$study_name[i]))
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(results, function(x) inherits(x, "study_error") || inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    first <- results[[which(failed)[1]]]
    stop(dataset, ": ", sum(failed), " studies failed, e.g. ", attr(first, "study"), ": ", first)
  }
  combined <- rbindlist(results, fill = TRUE)
  per_exposure <- unique(combined[, .(EXPOSURE, STATUS, CIS_LEADS, OUTSIDE_WINDOW, INSTRUMENTS, REGION_SNPS)])
  exposures <- merge(studies[, .(EXPOSURE, study_name, gene, ensg, probe, sample_size, tissue)], per_exposure,
                     by = "EXPOSURE", all.x = TRUE)
  fwrite(exposures, file.path(target, "exposures.tsv"), sep = "\t")

  top_hits <- combined[!is.na(SNP) & nzchar(SNP), ..top_hit_columns]
  if (nzchar(dbsnp_dir) && nrow(top_hits) > 0) top_hits <- add_rsids(top_hits)
  setorder(top_hits, EXPOSURE, P)
  fwrite(top_hits, file.path(target, "top_hits", "cis.tsv.gz"), sep = "\t")

  summary <- data.table(
    dataset = dataset, studies = nrow(studies),
    with_cis_lead = sum(exposures$CIS_LEADS > 0, na.rm = TRUE),
    cis_leads = sum(exposures$CIS_LEADS, na.rm = TRUE),
    ukb_leads_outside_1mb = sum(exposures$OUTSIDE_WINDOW, na.rm = TRUE),
    exposures_with_instruments = sum(exposures$INSTRUMENTS > 0, na.rm = TRUE),
    instruments = nrow(top_hits),
    instruments_with_rsid = if ("RSID" %in% names(top_hits)) sum(!is.na(top_hits$RSID)) else NA_integer_,
    gwas_files = length(list.files(out_gwas_dir, pattern = "\\.tsv\\.gz$")),
    no_region_files = sum(exposures$STATUS == "no_region_files", na.rm = TRUE),
    no_extracted_snps = sum(exposures$STATUS == "no_extracted_snps", na.rm = TRUE)
  )
  fwrite(summary, file.path(target, "summary.tsv"), sep = "\t")
  print(summary)
}
