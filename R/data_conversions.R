rsid_builds <- list(GRCh37="b37_dbsnp156", GRCh38="b38_dbsnp156")

#' @export
populate_gene_names <- function(gwas) {
  return(gwas)
  if ("ENSEMBL_ID" %in% colnames(gwas) && !"GENE_NAME" %in% colnames(gwas)) {
    return(ensembl_id_to_gene_name(gwas))
  } else if ("GENE_NAME" %in% colnames(gwas) && !"ENSEMBL_ID" %in% colnames(gwas)) {
    return(gene_name_to_ensembl_id(gwas))
  }
  return(gwas)
}

ensembl_id_to_gene_name <- function(gwas) {
  sqlite_db <- paste0(genomic_data_dir, "EnsDb.Hsapiens.v79.sqlite")
  ensembl_id_list <- toString(gwas$ENSEMBL_ID)
  command <- paste0("select gene_id || ',' || gene_name from gene where gene_id in (", ensembl_id_list, ")")

  gene_map <- run_sqlite_command(sqlite_db, command, c("ENSEMBL_ID", "GENE_NAME"))
  gwas$GENE_NAME <- gene_map$GENE_NAME[match(gwas$ENSEMBL_ID, gene_map$ENSEMBL_ID)]
}

gene_name_to_ensembl_id <- function(gwas) {
  sqlite_db <- paste0(genomic_data_dir, "EnsDb.Hsapiens.v79.sqlite")
  ensembl_id_list <- toString(gwas$GENE_NAME)
  command <- paste0("select gene_id || ',' || gene_name from gene where gene_id in (", ensembl_id_list, ")")

  gene_map <- run_sqlite_command(sqlite_db, command, c("ENSEMBL_ID", "GENE_NAME"))
  gwas$ENSEMBL_ID <- gene_map$ENSEMBL_ID[match(gwas$GENE_NAME, gene_map$GENE_NAME)]
}

#' @export
populate_rsid <- function(gwas, option = F) {
  start <- Sys.time()
  if (option == F || predefined_column_maps$default$RSID %in% colnames(gwas)) {
    message("Skipping RSID population for GWAS")
  } else {
    gwas <- populate_full_rsids(gwas)
  }

  message("RSID population time taken:")
  message(Sys.time() - start)
  return (gwas)
}

# populate_partial_rsids <- function(gwas) {
#   message("populating RSIDs based on 1000genomes...")
#   marker_to_rsid_file <- paste0(genomic_data_dir, "1000genomes/marker_to_rsid.tsv.gz")
#   chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select = c("HG37", "RSID"), show_col_types=F)
#   gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]
#
#   return(gwas)
# }

#' @import data.table
#' @import tibble
#' @import genepi.utils
#' @import future
populate_full_rsids <- function(gwas, build = rsid_builds$GRCh37) {
  dbsnp_dir <- paste0(genomic_data_dir, "dbsnp")
  if (!build %in% rsid_builds) stop(paste("Error: invalid rsid build option:", build))
  gc()

  gwas <- data.table::as.data.table(gwas)
  future::plan(future::multisession, workers = number_of_cpus_available)
  gwas <- genepi.utils::chrpos_to_rsid(gwas, "CHR", "BP", "EA", "OA", flip = "allow", dbsnp_dir=dbsnp_dir, build=build, alt_rsids = F)
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}
