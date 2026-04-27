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

#' Allele complement for strand alignment (SNPs only).
#' @keywords internal
allele_complement_vec <- function(alleles) {
  u <- toupper(alleles)
  r <- c(A = "T", T = "A", G = "C", C = "G", I = "I", D = "D")
  out <- unname(r[u])
  out[is.na(u) | u == ""] <- NA_character_
  out
}

#' Map effect allele frequency using reference A1 allele frequency (PLINK --freq).
#' @keywords internal
map_ref_to_eaf <- function(ea, oa, ra1, ra2, a1_freq) {
  ea <- toupper(ea)
  oa <- toupper(oa)
  ra1 <- toupper(ra1)
  ra2 <- toupper(ra2)
  dplyr::case_when(
    ea == ra1 & oa == ra2 ~ a1_freq,
    ea == ra2 & oa == ra1 ~ 1 - a1_freq,
    ea == allele_complement_vec(ra2) & oa == allele_complement_vec(ra1) ~ 1 - a1_freq,
    ea == allele_complement_vec(ra1) & oa == allele_complement_vec(ra2) ~ a1_freq,
    TRUE ~ NA_real_
  )
}

#' Populate missing EAF from the LD reference panel (1000 Genomes PLINK binary + .frq).
#'
#' Expects \code{plink --bfile <prefix> --freq --out <prefix>} to exist beside the
#' \code{.bed/.bim/.fam} files, i.e. \code{<thousand_genomes_dir>/<ANCESTRY>.frq}.
#'
#' @param gwas Standardised GWAS tibble (CHR, BP, EA, OA; SNP id built after allele ordering).
#' @param ancestry One of EUR, EAS, AFR, AMR, SAS matching the reference panel prefix.
#' @param enabled If FALSE, returns \code{gwas} unchanged.
#' @export
populate_eaf_from_reference_panel <- function(gwas, ancestry, enabled = FALSE) {
  if (!enabled) {
    return(gwas)
  }
  allowed <- c("EUR", "EAS", "AFR", "AMR", "SAS")
  if (!ancestry %in% allowed) {
    stop("ancestry must be one of: ", paste(allowed, collapse = ", "))
  }
  prefix <- file.path(thousand_genomes_dir, ancestry)
  bim_path <- paste0(prefix, ".bim")
  frq_path <- paste0(prefix, ".frq")
  if (!file.exists(bim_path)) {
    stop("LD reference .bim not found: ", bim_path)
  }
  if (!file.exists(frq_path)) {
    stop(
      "LD reference .frq not found: ", frq_path,
      "\nGenerate with: plink --bfile ", prefix, " --freq --out ", prefix
    )
  }

  ref_bim <- vroom::vroom(
    bim_path,
    delim = "\t",
    col_names = c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"),
    show_col_types = FALSE
  )
  ref_frq <- vroom::vroom(frq_path, show_col_types = FALSE)
  names(ref_frq) <- stringr::str_trim(names(ref_frq))
  req <- c("CHR", "SNP", "A1", "A2", "MAF")
  if (!all(req %in% colnames(ref_frq))) {
    stop(".frq file must contain columns CHR, SNP, A1, A2, MAF (plink 1.9 --freq)")
  }

  ref_bim <- dplyr::mutate(ref_bim, CHR = as.character(CHR))
  ref_frq <- dplyr::mutate(dplyr::select(ref_frq, dplyr::all_of(req)), CHR = as.character(CHR))

  ref <- dplyr::inner_join(ref_bim, ref_frq, by = c("CHR", "SNP"))
  ref <- dplyr::distinct(ref, CHR, BP, .keep_all = TRUE)
  ref <- dplyr::transmute(
    ref,
    CHR = as.character(CHR),
    BP = as.integer(BP),
    A1 = toupper(A1),
    A2 = toupper(A2),
    MAF = as.numeric(MAF)
  )

  gwas$CHR <- as.character(gwas$CHR)
  gwas$BP <- as.integer(gwas$BP)

  joined <- dplyr::left_join(gwas, ref, by = c("CHR", "BP"))
  computed_eaf <- map_ref_to_eaf(joined$EA, joined$OA, joined$A1, joined$A2, joined$MAF)

  if (!"EAF" %in% names(gwas)) {
    gwas$EAF <- NA_real_
  }
  need <- is.na(gwas$EAF)
  fill <- need & !is.na(computed_eaf)
  gwas$EAF[fill] <- computed_eaf[fill]

  n_need <- sum(need)
  n_fill <- sum(fill)
  message(
    "EAF population (", ancestry, " reference): filled ", n_fill,
    " of ", n_need, " missing values; ", n_need - n_fill, " still missing."
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
    parallel_cores=number_of_cpus_available
  )
  print('finished chrpos_to_rsid')
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}
