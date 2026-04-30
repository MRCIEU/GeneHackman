get_env_var <- function(env_var_name, default_value=NULL) {
  if (Sys.getenv(env_var_name) == "") {
    return(default_value)
  } else {
    return(Sys.getenv(env_var_name))
  }
}

user_data_dir <- get_env_var("DATA_DIR", "")
user_results_dir <- get_env_var("RESULTS_DIR", "")
pipeline_data_dir <- sub("/+$", "", get_env_var("PIPELINE_DATA_DIR", ""))
qtl_directory <- sub("/+$", "", get_env_var("QTL_DATA_DIR", ""))
number_of_cpus_available <- as.numeric(get_env_var("SLURM_CPUS_ON_NODE", 1))

genomic_data_dir <- file.path(pipeline_data_dir, "genomic_data")
# Must match snakemake/util/constants.smk THOUSAND_GENOMES_DIR (used by clumping + LD).
thousand_genomes_dir <- file.path(genomic_data_dir, "1000genomes", "b37_dbsnp156")

liftover_dir <- file.path(genomic_data_dir, "liftover")
pqtl_top_hits_dir <- file.path(qtl_directory, "pqtl")
metabrain_top_hits_dir <- file.path(qtl_directory, "metabrain", "top_hits")
metabrain_gwas_dir <- file.path(qtl_directory, "metabrain", "gwas")
eqtlgen_top_hits_dir <- file.path(qtl_directory, "eqtlgen", "top_hits")
eqtlgen_gwas_dir <- file.path(qtl_directory, "eqtlgen", "gwas")

qtl_datasets <- list(metabrain="metabrain", eqtlgen="eqtlgen")
populate_rsid_options <- list(none="none", partial="partial", full="full")
populate_chr_bp_options <- list(none="none", partial="partial", full="full")
rsid_builds <- list(GRCh37="b37_dbsnp156", GRCh38="b38_dbsnp156")
