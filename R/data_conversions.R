#' populate_gene_names: populate the gene names from the ENSEMBL_ID column
#' @param gwas dataframe with the following columns: ENSEMBL_ID
#' @return gwas with new column GENE_NAME


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
  gene_map <- vroom::vroom(file.path(genomic_data_dir, "gene_name_map.tsv"))
  gwas$GENE_NAME <- gene_map$GENE_NAME[match(gwas$ENSEMBL_ID, gene_map$ENSEMBL_ID)]
  return(gwas)
}

gene_name_to_ensembl_id <- function(gwas) {
  gene_map <- vroom::vroom(file.path(genomic_data_dir, "gene_name_map.tsv"))
  gwas$ENSEMBL_ID <- gene_map$ENSEMBL_ID[match(gwas$GENE_NAME, gene_map$GENE_NAME)]
  return(gwas)
}

#' Populate missing EAF from the LD reference panel \code{.frq} file via RSID matching.
#'
#' Expects the GWAS to have an \code{RSID} column (populated by \code{populate_rsid})
#' and alleles in canonical (alphabetically sorted) order from \code{standardise_alleles}.
#' Only reads the lightweight \code{.frq} file; the \code{.bim} is not needed.
#'
#' @param gwas Standardised GWAS tibble with RSID, EA, OA columns.
#' @param ancestry One of EUR, EAS, AFR, AMR, SAS matching the reference panel prefix.
#' @export
populate_eaf_from_reference_panel <- function(gwas, ancestry) {
  if ("EAF" %in% colnames(gwas) && !all(is.na(gwas$EAF))) {
    if (sum(is.na(gwas$EAF)) == 0) {
      message("EAF column already fully populated; skipping reference panel lookup.")
      return(gwas)
    }
  }

  allowed <- c("EUR", "EAS", "AFR", "AMR", "SAS")
  if (!ancestry %in% allowed) {
    stop("ancestry must be one of: ", paste(allowed, collapse = ", "))
  }
  frq_path <- paste0(file.path(thousand_genomes_dir, ancestry), ".frq")
  if (!file.exists(frq_path)) {
    stop("LD reference .frq not found: ", frq_path,
         "\nGenerate with: plink --bfile ",
         file.path(thousand_genomes_dir, ancestry), " --freq --out ",
         file.path(thousand_genomes_dir, ancestry))
  }

  if (!"RSID" %in% colnames(gwas) || all(is.na(gwas$RSID))) {
    warning("No RSIDs available for EAF lookup; skipping EAF population.")
    return(gwas)
  }

  if (!"EAF" %in% names(gwas)) {
    gwas$EAF <- NA_real_
  }
  need <- which(is.na(gwas$EAF) & !is.na(gwas$RSID))
  if (length(need) == 0L) {
    message("EAF population: no rows need filling (all have EAF or lack RSID).")
    return(gwas)
  }

  rsids_needed <- unique(gwas$RSID[need])
  message("Reading reference .frq for ", length(rsids_needed), " RSIDs...")

  ref <- data.table::fread(frq_path, header = TRUE, showProgress = FALSE, data.table = TRUE)
  names(ref) <- stringr::str_trim(names(ref))
  req <- c("SNP", "A1", "A2", "MAF")
  if (!all(req %in% names(ref))) {
    stop(".frq file must contain columns SNP, A1, A2, MAF (plink 1.9 --freq)")
  }
  ref <- ref[, .(SNP, A1 = toupper(A1), A2 = toupper(A2), MAF = as.numeric(MAF))]
  ref <- ref[SNP %in% rsids_needed]
  ref <- unique(ref, by = "SNP")

  if (nrow(ref) == 0L) {
    message("EAF population (", ancestry, " reference): 0 RSIDs matched in .frq; ",
            length(need), " still missing.")
    return(gwas)
  }

  idx <- match(gwas$RSID[need], ref$SNP)
  matched <- !is.na(idx)

  ea <- toupper(gwas$EA[need[matched]])
  oa <- toupper(gwas$OA[need[matched]])
  ra1 <- ref$A1[idx[matched]]
  ra2 <- ref$A2[idx[matched]]
  maf <- ref$MAF[idx[matched]]

  eaf <- data.table::fcase(
    ea == ra1 & oa == ra2, maf,
    ea == ra2 & oa == ra1, 1 - maf,
    default = NA_real_
  )

  gwas$EAF[need[matched]] <- eaf

  n_fill <- sum(!is.na(eaf))
  n_total_need <- sum(is.na(gwas$EAF)) + n_fill
  message(
    "EAF population (", ancestry, " reference): filled ", n_fill,
    " of ", n_total_need, " missing values; ", n_total_need - n_fill, " still missing."
  )
  gwas
}

#' populate_rsid: populate the RSID column from the SNP column
#' @param gwas dataframe with the following columns: SNP
#' @param option option to populate the RSID column
#' @return gwas with new column RSID


#' @export
populate_rsid <- function(gwas, option = populate_rsid_options$none) {
  gc()
  start_time <- Sys.time()
  if (option == populate_rsid_options$none || "RSID" %in% colnames(gwas)) {
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
  marker_to_rsid_file <- file.path(thousand_genomes_dir, "marker_to_rsid.tsv.gz")
  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select = c("HG37", "RSID"), show_col_types=F)
  gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]

  return(gwas)
}




populate_full_rsids <- function(gwas) {
  build <- rsid_builds$GRCh37
  dbsnp_dir <- file.path(genomic_data_dir, "dbsnp")
  if (!build %in% rsid_builds) stop(paste("Error: invalid rsid build option:", build))

  parallel_cores <- ifelse (locally_running, 1, number_of_cpus_available)
  message(paste0("Using ", parallel_cores, " cores for full RSID population"))

  gwas <- data.table::as.data.table(gwas)
  gwas <- chrpos_to_rsid(
    gwas,
    "CHR",
    "BP", 
    "EA", 
    "OA",
    flip = "allow",
    dbsnp_dir = dbsnp_dir,
    build = build,
    alt_rsids = FALSE,
    parallel_cores=parallel_cores
  )
  print('finished chrpos_to_rsid')
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}
