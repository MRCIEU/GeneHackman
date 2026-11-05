# Make a copy of GeneHackman's QTL data in another reference build, with GeneHackman's own liftover and QC:
#
#   Rscript 30_lift_qtl_data.R [--datasets eqtlgen,brainmeta,ukb_pqtl] [--from GRCh38] [--to GRCh37] [--cores 8]
#
# Reads NEW_QTL_DATA_DIR/<dataset>/<from>/{top_hits,gwas}/... (written by 31_qtl_from_gpm.R) and writes
# NEW_QTL_DATA_DIR/<dataset>/<to>/... with the same file names.  Each file is lifted with
# convert_reference_build_via_liftover() (the QC a GWAS gets: unmapped, other-contig, chromosome-change, strand-flip,
# round-trip and duplicate variants are dropped) and its SNP keys are rebuilt in the new build.  RSIDs don't depend
# on the build and are kept.  Finished files are skipped on a re-run.  Counts per file go to
# NEW_QTL_DATA_DIR/<dataset>/<to>/liftover_summary.tsv.  The chains are read from CHAIN_DIR (both directions, for
# the round-trip check).
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  args[i + 1]
}
datasets <- strsplit(get_arg("--datasets", "eqtlgen,brainmeta,ukb_pqtl"), ",")[[1]]
from_build <- get_arg("--from", "GRCh38")
to_build <- get_arg("--to", "GRCh37")
cores <- as.integer(get_arg("--cores", "8"))
qtl_dir <- Sys.getenv("NEW_QTL_DATA_DIR")
if (!nzchar(qtl_dir)) stop("set NEW_QTL_DATA_DIR (see config.sh)")

suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))
if (nzchar(Sys.getenv("CHAIN_DIR"))) options(GeneHackman.liftover_chain_dir = Sys.getenv("CHAIN_DIR"))
stopifnot(from_build %in% reference_builds, to_build %in% reference_builds, from_build != to_build)

lift_file <- function(in_file, out_file) {
  if (file.exists(out_file)) return(NULL)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  qtl <- as.data.frame(fread(in_file, colClasses = list(character = c("CHR", "SNP", "EA", "OA"))))
  n_in <- nrow(qtl)
  lifted <- suppressMessages(convert_reference_build_via_liftover(qtl, from_build, to_build, max_drop_fraction = 1))
  attr(lifted, "liftover_qc") <- NULL
  lifted[[paste0("BP", sub("GRCh", "", from_build))]] <- NULL
  if (nrow(lifted) > 0) {
    lifted$SNP <- paste0(lifted$CHR, ":", format(lifted$BP, scientific = FALSE, trim = TRUE), "_",
                         toupper(lifted$EA), "_", toupper(lifted$OA))
  }
  tmp <- paste0(out_file, ".tmp.gz")
  fwrite(lifted, tmp, sep = "\t")
  file.rename(tmp, out_file)
  return(data.table(FILE = sub(paste0("^", qtl_dir, "/"), "", in_file), INPUT = n_in, KEPT = nrow(lifted)))
}

for (dataset in datasets) {
  source_dir <- file.path(qtl_dir, dataset, from_build)
  target <- file.path(qtl_dir, dataset, to_build)
  files <- list.files(source_dir, pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("^(top_hits|gwas)/", sub(paste0("^", source_dir, "/"), "", files))]
  if (length(files) == 0) stop("No QTL files under ", source_dir, "/{top_hits,gwas}")
  message(dataset, ": lifting ", length(files), " files from ", from_build, " to ", to_build)
  out_files <- file.path(target, sub(paste0("^", source_dir, "/"), "", files))
  results <- parallel::mclapply(seq_along(files), function(i) {
    tryCatch(lift_file(files[i], out_files[i]), error = function(e) structure(conditionMessage(e), class = "lift_error"))
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(results, function(x) inherits(x, "lift_error") || inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sum(failed), " files failed, e.g. ", files[which(failed)[1]], ": ", results[[which(failed)[1]]])
  for (extra in c("exposures.tsv", "summary.tsv")) {
    if (file.exists(file.path(source_dir, extra))) file.copy(file.path(source_dir, extra), file.path(target, extra),
                                                             overwrite = TRUE)
  }
  summary <- rbindlist(results)
  if (nrow(summary) > 0) {
    summary_file <- file.path(target, "liftover_summary.tsv")
    fwrite(summary, summary_file, sep = "\t", append = file.exists(summary_file))
  }
  all_summary <- fread(file.path(target, "liftover_summary.tsv"))
  message(dataset, ": ", nrow(all_summary), " files lifted; variants kept ", sum(all_summary$KEPT), " of ",
          sum(all_summary$INPUT), sprintf(" (%.2f%% dropped)", 100 * (1 - sum(all_summary$KEPT) / sum(all_summary$INPUT))))
}
