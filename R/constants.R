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
    floor(mem / memory_per_worker_mb)
  } else {
    cpus
  }
  max(1L, min(max_workers, by_mem, cpus))
}

.strip_path <- function(x) {
  sub("/+$", "", trimws(x))
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
# Must match snakemake/util/constants.smk THOUSAND_GENOMES_DIR (used by clumping + LD).
thousand_genomes_dir <- file.path(genomic_data_dir, "1000genomes", "b37_dbsnp156")

liftover_chain_dir <- file.path(genomic_data_dir, "liftover")
liftover_binary <- Sys.which("liftOver")
if (!nzchar(liftover_binary)) {
  liftover_binary <- "/usr/local/bin/liftOver"
}
pqtl_top_hits_dir <- file.path(qtl_directory, "pqtl")
metabrain_top_hits_dir <- file.path(qtl_directory, "metabrain", "top_hits")
metabrain_gwas_dir <- file.path(qtl_directory, "metabrain", "gwas")
eqtlgen_top_hits_dir <- file.path(qtl_directory, "eqtlgen", "top_hits")
eqtlgen_gwas_dir <- file.path(qtl_directory, "eqtlgen", "gwas")

qtl_datasets <- list(metabrain="metabrain", eqtlgen="eqtlgen")
populate_rsid_options <- list(none="none", partial="partial", full="full")
populate_chr_bp_options <- list(none="none", partial="partial", full="full")
rsid_builds <- list(GRCh37="b37_dbsnp156", GRCh38="b38_dbsnp156")
study_types <- list(continuous="continuous", categorical="categorical")