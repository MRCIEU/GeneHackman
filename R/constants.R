get_env_var <- function(env_var_name, default_value=NULL) {
  if (Sys.getenv(env_var_name) == "") {
    return(default_value)
  } else {
    return(Sys.getenv(env_var_name))
  }
}

# Slurm allocated memory (MB) when running under Slurm, else total system memory.
available_memory <- function() {
  slurm <- Sys.getenv("SLURM_MEM_PER_NODE", "")
  if (nzchar(slurm)) {
    v <- suppressWarnings(as.numeric(slurm))
    if (!is.na(v) && v > 0) return(v)
  }
  os <- Sys.info()[["sysname"]]
  if (os == "Linux") {
    meminfo <- tryCatch(readLines("/proc/meminfo", n = 1), error = function(e) "")
    kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", meminfo)))
    if (!is.na(kb) && kb > 0) return(kb / 1024)
  } else if (os == "Darwin") {
    raw <- tryCatch(system2("sysctl", "-n hw.memsize", stdout = TRUE, stderr = FALSE),
                    error = function(e) "")
    bytes <- suppressWarnings(as.numeric(raw))
    if (!is.na(bytes) && bytes > 0) return(bytes / (1024^2))
  }
  return(NA)
}

# Slurm node CPUs when SLURM_* env is set, else parallel::detectCores().
available_cpus <- function() {
  slurm <- Sys.getenv("SLURM_CPUS_ON_NODE", "")
  if (nzchar(slurm)) {
    v <- suppressWarnings(as.integer(slurm))
    if (!is.na(v)) {
      return(max(1L, v))
    }
  }
  dc <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(dc) || dc < 1L) {
    return(1L)
  } else {
    return(as.integer(dc))
  }
}

# Shared by RSID population, finemap, etc. (load.R sources all R/*.R; keep only this definition).
calculate_parallelism <- function(max_workers = 10L, memory_per_worker_mb = 8000L) {
  cpus <- available_cpus()
  mem <- available_memory()
  by_mem <- if (!is.na(mem) && mem > 0) {
    round(floor(mem / memory_per_worker_mb))
  } else {
    return(cpus)
  }
  return(max(1L, min(max_workers, by_mem, cpus)))
}

.strip_path <- function(x) {
  return(sub("/+$", "", trimws(x)))
}

project_dir <- .strip_path(get_env_var("PROJECT_DIR", ""))
if (nzchar(project_dir)) {
  user_data_dir <- file.path(project_dir, "data")
  user_results_dir <- file.path(project_dir, "results")
} else {
  user_data_dir <- ""
  user_results_dir <- ""
}

pipeline_data_dir <- sub("/+$", "", get_env_var("PIPELINE_DATA_DIR", ""))
qtl_directory <- sub("/+$", "", get_env_var("QTL_DATA_DIR", ""))
is_on_cluster <- nzchar(Sys.getenv("SLURM_JOB_ID")) || nzchar(Sys.getenv("PBS_JOBID"))

genomic_data_dir <- file.path(pipeline_data_dir, "genomic_data")
# 1000 Genomes genotype panels (PLINK bfiles, .frq, marker_to_rsid) per reference build; must match
# thousand_genomes_dir() in snakemake/util/constants.smk.  GRCh38 is the 1000G 30x release (main reference
# data download); GRCh37 is phase 3 (download_reference_data.sh --grch37, served from release 1.1.0).
thousand_genomes_builds <- list(GRCh38 = "b38_dbsnp156", GRCh37 = "b37_dbsnp156")
# Must match snakemake/util/constants.smk LD_PACK_DIR (precomputed Pan-UKB LD, downloaded with download_reference_data.sh).
ld_pack_dir <- file.path(genomic_data_dir, "pan_ukb_ld")

liftover_chain_dir <- file.path(genomic_data_dir, "liftover")
liftover_binary <- Sys.which("liftOver")
if (!nzchar(liftover_binary)) {
  liftover_binary <- "/usr/local/bin/liftOver"
}
# QTL datasets, cis QTLs only, at QTL_DATA_DIR/<dataset>/<build>/{top_hits,gwas}/cis... (published in GRCh38 only)
# (from the genotype-phenotype map 1.0.0 study data).
# Must match qtl_dataset_options in snakemake/util/constants.smk.
qtl_datasets <- list(eqtlgen="eqtlgen", brainmeta="brainmeta", ukb_pqtl="ukb_pqtl")
# Only cis QTLs are used as instruments: trans QTLs are too often pleiotropic for MR.
qtl_subcategory <- "cis"
populate_rsid_options <- list(none="none", partial="partial", full="full")
populate_chr_bp_options <- list(none="none", partial="partial", full="full")
rsid_builds <- list(GRCh37="b37_dbsnp156", GRCh38="b38_dbsnp156")
study_types <- list(continuous="continuous", categorical="categorical", binary="binary")
#' Directory of the 1000 Genomes panel for a reference build.
#'
#' @param build reference build of the positions the panel is used with (one of \code{reference_builds})
#' @param check if TRUE, stop with download instructions when the panel is not installed
#' @keywords internal
thousand_genomes_dir_for <- function(build = reference_builds$GRCh38, check = FALSE) {
  subdir <- thousand_genomes_builds[[build]]
  if (is.null(subdir)) {
    stop("No 1000 Genomes panel for reference build ", build, "; available: ",
         paste(names(thousand_genomes_builds), collapse = ", "))
  }
  dir <- file.path(genomic_data_dir, "1000genomes", subdir)
  if (check && !dir.exists(dir)) {
    stop("1000 Genomes ", build, " panel not found at ", dir, ". ",
         "Download it with ./download_reference_data.sh --ancestries <ancestries>",
         if (build == reference_builds$GRCh37) " --grch37 (GRCh37 output)" else "", ".")
  }
  return(dir)
}

#' Directories of a QTL dataset's files for a reference build: \verb{QTL_DATA_DIR/<dataset>/<build>/<kind>}.
#' The published QTL data are GRCh38 only.
#'
#' @param dataset one of \code{qtl_datasets}
#' @param kind "top_hits" (MR instruments) or "gwas" (full summary statistics, for coloc)
#' @param build reference build of the GWAS the QTLs are compared with
#' @keywords internal
qtl_dataset_dir <- function(dataset, kind, build = reference_builds$GRCh38) {
  dir <- file.path(qtl_directory, dataset, build, kind)
  if (!dir.exists(dir)) {
    stop("QTL data for ", dataset, " in ", build, " not found at ", dir, ". ",
         if (build != reference_builds$GRCh38) "The QTL data are GRCh38 only: use output build GRCh38."
         else paste0("Download them with ./download_reference_data.sh --qtl ", dataset, " and set QTL_DATA_DIR."))
  }
  return(dir)
}
