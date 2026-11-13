#' Check and fix strand against a reference panel's alleles and allele frequencies.
#'
#' Run on a GWAS with GeneHackman SNP keys in canonical allele order (after \code{standardise_alleles()}).  Each
#' variant is compared with the reference (the same build):
#' \describe{
#'   \item{match}{the key is in the reference.}
#'   \item{strand_flipped}{a non-palindromic SNV whose key matches only after complementing both alleles
#'     (A<->T, C<->G): the alleles are complemented, and BETA and EAF flipped if that changes the allele order.}
#'   \item{palindromic_flipped}{an A/T or C/G SNV whose EAF is on the other side of 0.5 from the reference's,
#'     with both minor allele frequencies below \code{palindromic_threshold}: reported on the other strand,
#'     so BETA and EAF are flipped.}
#'   \item{palindromic_ambiguous}{an A/T or C/G SNV with a minor allele frequency at or above
#'     \code{palindromic_threshold} (0.42, as TwoSampleMR) in the GWAS or the reference: its strand can't be told
#'     from frequency.  Kept unless \code{drop_ambiguous_palindromic}.}
#'   \item{palindromic_unchecked}{an A/T or C/G SNV without a GWAS EAF to compare.}
#'   \item{eaf_discordant}{alleles match but the EAF differs from the reference's by more than
#'     \code{max_eaf_difference}: reported, not changed (it may be a different population, or an allele-coding
#'     error in the file).}
#'   \item{not_in_reference}{kept, not checked.}
#' }
#' @param gwas standardised GWAS (SNP, EA, OA, BETA and optionally EAF, Z)
#' @param reference data.frame with SNP (key in canonical allele order) and REF_EAF (frequency of the first
#'   allele of the key)
#' @param palindromic_threshold minor allele frequency from which palindromic SNVs are ambiguous
#' @param drop_ambiguous_palindromic remove ambiguous palindromic SNVs
#' @param max_eaf_difference EAF difference reported as discordant
#' @return list with \code{gwas} (fixed, with an ALLELE_QC column) and \code{report} (counts per status)
#' @export
harmonise_against_reference <- function(gwas, reference, palindromic_threshold = 0.42,
                                        drop_ambiguous_palindromic = FALSE, max_eaf_difference = 0.2) {
  ea <- toupper(gwas$EA)
  oa <- toupper(gwas$OA)
  snv <- nchar(ea) == 1 & nchar(oa) == 1 & ea %in% c("A", "C", "G", "T") & oa %in% c("A", "C", "G", "T")
  palindromic <- snv & complement_alleles(ea) == oa
  ref_eaf <- reference$REF_EAF[match(gwas$SNP, reference$SNP)]
  in_reference <- gwas$SNP %in% reference$SNP
  status <- ifelse(in_reference, "match", "not_in_reference")

  # non-palindromic SNVs reported on the other strand
  comp_ea <- complement_alleles(ea)
  comp_oa <- complement_alleles(oa)
  swap <- comp_ea > comp_oa
  comp_key <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = FALSE, trim = TRUE), "_",
                             ifelse(swap, comp_oa, comp_ea), "_", ifelse(swap, comp_ea, comp_oa)))
  strand_flip <- !in_reference & snv & !palindromic & comp_key %in% reference$SNP
  if (any(strand_flip)) {
    flip_sign <- strand_flip & swap
    gwas$EA[strand_flip] <- ifelse(swap, comp_oa, comp_ea)[strand_flip]
    gwas$OA[strand_flip] <- ifelse(swap, comp_ea, comp_oa)[strand_flip]
    gwas <- flip_effects(gwas, flip_sign)
    gwas$SNP[strand_flip] <- comp_key[strand_flip]
    ref_eaf[strand_flip] <- reference$REF_EAF[match(comp_key[strand_flip], reference$SNP)]
    status[strand_flip] <- "strand_flipped"
  }

  has_eaf <- "EAF" %in% colnames(gwas)
  eaf <- if (has_eaf) as.numeric(gwas$EAF) else rep(NA_real_, nrow(gwas))
  checked_palindromic <- palindromic & in_reference
  status[checked_palindromic & is.na(eaf)] <- "palindromic_unchecked"
  comparable <- checked_palindromic & !is.na(eaf) & !is.na(ref_eaf)
  maf <- pmax(pmin(eaf, 1 - eaf), pmin(ref_eaf, 1 - ref_eaf))
  ambiguous <- comparable & maf >= palindromic_threshold
  other_strand <- comparable & !ambiguous & (eaf - 0.5) * (ref_eaf - 0.5) < 0
  status[ambiguous] <- "palindromic_ambiguous"
  status[other_strand] <- "palindromic_flipped"
  gwas <- flip_effects(gwas, other_strand)
  if (has_eaf) eaf <- as.numeric(gwas$EAF)

  discordant <- status %in% c("match", "strand_flipped") & !palindromic & !is.na(eaf) & !is.na(ref_eaf) &
    abs(eaf - ref_eaf) > max_eaf_difference
  status[discordant] <- "eaf_discordant"

  gwas$ALLELE_QC <- status
  if (drop_ambiguous_palindromic) gwas <- gwas[status != "palindromic_ambiguous", , drop = FALSE]
  statuses <- c("match", "strand_flipped", "palindromic_flipped", "palindromic_ambiguous", "palindromic_unchecked",
                "eaf_discordant", "not_in_reference")
  report <- data.frame(STATUS = statuses, N = vapply(statuses, function(x) sum(status == x), numeric(1)),
                       ACTION = c("kept", "alleles complemented", "BETA and EAF flipped",
                                  if (drop_ambiguous_palindromic) "removed" else "kept", "kept", "kept (reported)",
                                  "kept"))
  rownames(report) <- NULL
  list(gwas = gwas, report = report)
}

#' Flip BETA, Z and EAF of the rows in \code{which} (the effect allele is the other allele).
#' @keywords internal
flip_effects <- function(gwas, which) {
  if (!any(which)) return(gwas)
  gwas$BETA[which] <- -gwas$BETA[which]
  if ("Z" %in% colnames(gwas)) gwas$Z[which] <- -gwas$Z[which]
  if ("EAF" %in% colnames(gwas)) gwas$EAF[which] <- 1 - gwas$EAF[which]
  gwas
}

#' Reference alleles and frequencies for \code{harmonise_against_reference()}, from the LD reference.
#'
#' Pan-UKB: the pack's variant table (SNP37/SNP38 keys and EAF).  1000 Genomes: the panel's .bim positions and
#' .frq allele frequencies, as keys in canonical allele order.  The user's own genotypes (a PLINK prefix as
#' \code{ld_reference}): their .bim and frequencies computed with PLINK (\code{genotype_allele_frequencies()}).
#' @return data.frame with SNP and REF_EAF, or NULL when the reference isn't installed
#' @keywords internal
allele_reference <- function(chromosomes, ancestry, ld_reference = ld_reference_options$thousand_genomes,
                             build = reference_builds$GRCh38) {
  tryCatch({
    if (is_ld_genotypes(ld_reference)) return(genotype_allele_frequencies(ld_reference, chromosomes))
    if (validate_ld_reference(ld_reference) == ld_reference_options$pan_ukb) {
      key_column <- ld_pack_key_column(build)
      available <- ld_pack_chromosomes(ancestry)
      tables <- lapply(intersect(as.character(chromosomes), available), function(chr) {
        variants <- load_ld_pack_variants(chr, ancestry)
        data.frame(SNP = variants[[key_column]], REF_EAF = variants$EAF)
      })
      reference <- do.call(rbind, tables)
      return(reference[!is.na(reference$SNP), , drop = FALSE])
    }
    panel <- file.path(thousand_genomes_dir_for(build), resolve_ld_ancestry(ancestry, ld_reference))
    if (!file.exists(paste0(panel, ".frq")) || !file.exists(paste0(panel, ".bim"))) return(NULL)
    bim <- data.table::fread(paste0(panel, ".bim"), header = FALSE, select = c(1, 4), col.names = c("CHR", "BP"),
                             colClasses = list(character = 1))
    frq <- data.table::fread(paste0(panel, ".frq"), header = TRUE, select = c("A1", "A2", "MAF"))
    keep <- bim$CHR %in% as.character(chromosomes)
    bim <- bim[keep]
    frq <- frq[keep]
    first <- pmin(frq$A1, frq$A2)
    data.frame(SNP = toupper(paste0(bim$CHR, ":", format(bim$BP, scientific = FALSE, trim = TRUE), "_", first, "_",
                                    pmax(frq$A1, frq$A2))),
               REF_EAF = ifelse(frq$A1 == first, frq$MAF, 1 - frq$MAF))
  }, error = function(e) {
    message("Allele check against the reference skipped: ", conditionMessage(e))
    NULL
  })
}

.genotype_frequency_cache <- new.env(parent = emptyenv())

#' Allele frequencies of the user's own genotypes, as GeneHackman SNP keys and the frequency of the key's first allele.
#'
#' Frequencies are computed once per PLINK prefix with \code{plink1.9 --freq} and kept for the R session.  Used for
#' the strand check (\code{allele_reference()}) and to fill EAF (\code{populate_eaf_from_reference_panel()}).
#' @param bfile PLINK prefix (.bed/.bim/.fam) in the build of the GWAS
#' @param chromosomes chromosomes to keep (NULL for all)
#' @return data.frame with SNP (key in canonical allele order) and REF_EAF (frequency of its first allele)
#' @keywords internal
genotype_allele_frequencies <- function(bfile, chromosomes = NULL) {
  if (is.null(.genotype_frequency_cache[[bfile]])) {
    out_prefix <- tempfile(pattern = "genotype_freq_")
    on.exit(unlink(paste0(out_prefix, c(".frq", ".log", ".nosex"))), add = TRUE)
    cmd <- paste("plink1.9", "--bfile", bfile, "--freq", "--out", out_prefix)
    exit_code <- run_system(cmd, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (exit_code != 0 || !file.exists(paste0(out_prefix, ".frq"))) {
      stop("PLINK could not compute allele frequencies from ", bfile, " (", cmd, ")")
    }
    bim <- read_bim(bfile)
    frq <- data.table::fread(paste0(out_prefix, ".frq"), header = TRUE, select = c("A1", "A2", "MAF"),
                             colClasses = list(character = c("A1", "A2")))
    chr <- sub("^chr", "", bim$V1, ignore.case = TRUE)
    a1 <- toupper(frq$A1)
    a2 <- toupper(frq$A2)
    first <- pmin(a1, a2)
    .genotype_frequency_cache[[bfile]] <- data.frame(
      CHR = chr,
      SNP = paste0(chr, ":", format(bim$V4, scientific = FALSE, trim = TRUE), "_", first, "_", pmax(a1, a2)),
      REF_EAF = ifelse(a1 == first, frq$MAF, 1 - frq$MAF)
    )
  }
  reference <- .genotype_frequency_cache[[bfile]]
  if (!is.null(chromosomes)) reference <- reference[reference$CHR %in% as.character(chromosomes), , drop = FALSE]
  reference <- reference[!duplicated(reference$SNP), , drop = FALSE]
  return(reference[, c("SNP", "REF_EAF")])
}
