#' @export
populate_gene_names <- function(gwas) {
  if ("ENSEMBL_ID" %in% colnames(gwas) && !"GENE_NAME" %in% colnames(gwas)) {
    return(ensembl_id_to_gene_name(gwas))
  } else if ("GENE_NAME" %in% colnames(gwas) && !"ENSEMBL_ID" %in% colnames(gwas)) {
    return(gene_name_to_ensembl_id(gwas))
  }
  return(gwas)
}

ensembl_id_to_gene_name <- function(gwas) {
  gene_map <- vroom::vroom(paste0(genomic_data_dir, "gene_name_map.tsv"))
  gwas$GENE_NAME <- gene_map$GENE_NAME[match(gwas$ENSEMBL_ID, gene_map$ENSEMBL_ID)]
  return(gwas)
}

gene_name_to_ensembl_id <- function(gwas) {
  gene_map <- vroom::vroom(paste0(genomic_data_dir, "gene_name_map.tsv"))
  gwas$ENSEMBL_ID <- gene_map$ENSEMBL_ID[match(gwas$GENE_NAME, gene_map$GENE_NAME)]
  return(gwas)
}

#' @export
populate_rsid <- function(gwas, option = populate_rsid_options$none) {
  gc()
  start_time <- Sys.time()
  if (option == populate_rsid_options$none || "RSID" %in% colnames(gwas)) {
    message("Skipping RSID population for GWAS")
    return(gwas)
  } else if (option == populate_rsid_options$partial) {
    gwas <- populate_partial_rsids(gwas)
  } else if (option == populate_rsid_options$full){
    gwas <- populate_full_rsids(gwas)
  }
  else {
    stop(paste("Invalid RSID population option", option))
  }

  print(paste0("RSID population option: ", option, ". Time taken:"))
  print(Sys.time() - start_time)
  return(gwas)
}

populate_partial_rsids <- function(gwas) {
  message("populating RSIDs based on 1000genomes...")
  marker_to_rsid_file <- paste0(genomic_data_dir, "1000genomes/marker_to_rsid.tsv.gz")
  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select = c("HG37", "RSID"), show_col_types=F)
  gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]

  return(gwas)
}

#' @import data.table
#' @import tibble
#' @import genepi.utils
#' @import future
populate_full_rsids <- function(gwas, build = rsid_builds$GRCh37) {
  dbsnp_dir <- paste0(genomic_data_dir, "dbsnp")
  if (!build %in% rsid_builds) stop(paste("Error: invalid rsid build option:", build))

  gwas <- data.table::as.data.table(gwas)
  gwas <- genepi.utils::chrpos_to_rsid(gwas, "CHR", "BP", "EA", "OA", flip = "allow", dbsnp_dir=dbsnp_dir, build=build, alt_rsids = F, parallel_cores=number_of_cpus_available)
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}
