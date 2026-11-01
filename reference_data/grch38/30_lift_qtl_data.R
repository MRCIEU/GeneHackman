# Make GRCh38 copies of GeneHackman's QTL data (eQTLGen, MetaBrain), which are used in the output build.
#
#   Rscript 30_lift_qtl_data.R [--datasets eqtlgen,metabrain] [--cores 8]
#
# Reads OLD_QTL_DATA_DIR/<dataset>/{top_hits,gwas}/... (the 1.1.0 layout) and writes
# NEW_QTL_DATA_DIR/<dataset>/GRCh38/{top_hits,gwas}/... with the same file names.  For each dataset it first
# works out which build the files are in (matching a sample of variants against the Pan-UKB LD pack's GRCh37 and
# GRCh38 keys): files already in GRCh38 are copied, GRCh37 files are lifted with GeneHackman's own liftover
# (convert_reference_build_via_liftover(), with the same QC as a GWAS) and their SNP keys rebuilt.
# Finished files are skipped on a re-run.  A summary goes to NEW_QTL_DATA_DIR/<dataset>/GRCh38/liftover_summary.tsv.
# The GRCh37 originals are put at NEW_QTL_DATA_DIR/<dataset>/GRCh37 by RUNBOOK.md (hard links, no copy).
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0) return(default)
  args[i + 1]
}
datasets <- strsplit(get_arg("--datasets", "eqtlgen,metabrain"), ",")[[1]]
cores <- as.integer(get_arg("--cores", "8"))
old_dir <- Sys.getenv("OLD_QTL_DATA_DIR")
new_dir <- Sys.getenv("NEW_QTL_DATA_DIR")
pack_dir <- Sys.getenv("LD_PACK_DIR")
if (!nzchar(old_dir) || !nzchar(new_dir)) stop("set OLD_QTL_DATA_DIR and NEW_QTL_DATA_DIR (see config.sh)")

suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))
if (nzchar(Sys.getenv("CHAIN_DIR"))) options(GeneHackman.liftover_chain_dir = Sys.getenv("CHAIN_DIR"))

detect_build <- function(files) {
  sample <- rbindlist(lapply(head(files[order(-file.size(files))], 5), fread, select = c("CHR", "BP", "EA", "OA"),
                             colClasses = list(character = "CHR")))
  share <- tryCatch(check_declared_build(sample, reference_builds$GRCh37, pack_dir),
                    error = function(e) c(GRCh37 = 0, GRCh38 = 1))
  if (is.null(share)) stop("Could not check the QTL build: is LD_PACK_DIR (", pack_dir, ") a Pan-UKB LD pack?")
  return(if (share[["GRCh38"]] > share[["GRCh37"]]) reference_builds$GRCh38 else reference_builds$GRCh37)
}

lift_file <- function(in_file, out_file, build) {
  if (file.exists(out_file)) return(NULL)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (build == reference_builds$GRCh38) {
    file.copy(in_file, out_file)
    return(data.table(FILE = in_file, INPUT = NA_integer_, KEPT = NA_integer_))
  }
  qtl <- as.data.frame(fread(in_file, colClasses = list(character = c("CHR", "SNP"))))
  n_in <- nrow(qtl)
  lifted <- suppressMessages(convert_reference_build_via_liftover(qtl, reference_builds$GRCh37, reference_builds$GRCh38,
                                                                  max_drop_fraction = 1))
  attr(lifted, "liftover_qc") <- NULL
  lifted$BP37 <- NULL
  if ("SNP" %in% colnames(lifted)) {
    lifted$SNP <- toupper(paste0(lifted$CHR, ":", format(lifted$BP, scientific = FALSE, trim = TRUE), "_",
                                 sub("^[^_]*_", "", lifted$SNP)))
  }
  tmp <- paste0(out_file, ".tmp.gz")
  fwrite(lifted, tmp, sep = "\t")
  file.rename(tmp, out_file)
  return(data.table(FILE = in_file, INPUT = n_in, KEPT = nrow(lifted)))
}

for (dataset in datasets) {
  source_dir <- file.path(old_dir, dataset)
  files <- list.files(source_dir, pattern = "\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("^(top_hits|gwas)/", sub(paste0("^", source_dir, "/"), "", files))]
  if (length(files) == 0) stop("No QTL files under ", source_dir, "/{top_hits,gwas}")
  build <- detect_build(files[grepl("/gwas/", files)])
  message(dataset, ": ", length(files), " files, in ", build, if (build == "GRCh38") " (copied)" else " (lifted)")
  target <- file.path(new_dir, dataset, reference_builds$GRCh38)
  out_files <- file.path(target, sub(paste0("^", source_dir, "/"), "", files))
  results <- parallel::mclapply(seq_along(files), function(i) lift_file(files[i], out_files[i], build),
                                mc.cores = cores)
  failed <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sum(failed), " files failed, e.g. ", files[which(failed)[1]], ": ", results[[which(failed)[1]]])
  summary <- rbindlist(results)
  if (nrow(summary) > 0) {
    fwrite(summary, file.path(target, "liftover_summary.tsv"), sep = "\t", append = file.exists(file.path(target, "liftover_summary.tsv")))
    message(dataset, ": ", nrow(summary), " files written; variants kept ",
            sum(summary$KEPT, na.rm = TRUE), " of ", sum(summary$INPUT, na.rm = TRUE))
  }
}
