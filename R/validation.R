#' Compare GeneHackman's harmonised GWAS with another tool's (e.g. MungeSumstats) for the same input file.
#'
#' Variants are matched on chromosome, position and the unordered allele pair, so the comparison doesn't depend
#' on which allele each tool chose as the effect allele.  The other tool's effects are then aligned to
#' GeneHackman's effect allele.
#'
#' @param genehackman GeneHackman output: CHR, BP, EA, OA, BETA and RSID
#' @param other the other tool's output: CHR, BP, its effect allele (\code{other_ea}), other allele
#'   (\code{other_oa}), BETA and its rsID column (\code{other_rsid})
#' @param other_ea,other_oa,other_rsid column names in \code{other} (MungeSumstats: A2 is the effect allele, A1 the
#'   other, SNP the rsID)
#' @return list with \code{summary} (MEASURE, VALUE) and \code{discordant} (matched variants whose rsID or effect
#'   direction differ)
#' @export
compare_harmonisation <- function(genehackman, other, other_ea = "A2", other_oa = "A1", other_rsid = "SNP") {
  key <- function(chr, bp, a, b) {
    a <- toupper(a); b <- toupper(b)
    paste0(sub("^chr", "", as.character(chr), ignore.case = TRUE), ":", format(as.numeric(bp), scientific = FALSE,
           trim = TRUE), "_", pmin(a, b), "_", pmax(a, b))
  }
  gh <- data.frame(KEY = key(genehackman$CHR, genehackman$BP, genehackman$EA, genehackman$OA),
                   RSID_GENEHACKMAN = if ("RSID" %in% names(genehackman)) as.character(genehackman$RSID) else NA,
                   EA_GENEHACKMAN = toupper(genehackman$EA), BETA_GENEHACKMAN = as.numeric(genehackman$BETA))
  ot <- data.frame(KEY = key(other$CHR, other$BP, other[[other_ea]], other[[other_oa]]),
                   RSID_OTHER = as.character(other[[other_rsid]]), EA_OTHER = toupper(other[[other_ea]]),
                   BETA_OTHER = as.numeric(other$BETA))
  gh <- gh[!duplicated(gh$KEY), ]
  ot <- ot[!duplicated(ot$KEY), ]
  matched <- merge(gh, ot, by = "KEY")
  matched$BETA_OTHER_ALIGNED <- ifelse(matched$EA_OTHER == matched$EA_GENEHACKMAN, matched$BETA_OTHER,
                                       -matched$BETA_OTHER)
  both_rsid <- !is.na(matched$RSID_GENEHACKMAN) & nzchar(matched$RSID_GENEHACKMAN) &
    !is.na(matched$RSID_OTHER) & grepl("^rs", matched$RSID_OTHER)
  same_rsid <- matched$RSID_GENEHACKMAN == matched$RSID_OTHER
  nonzero <- matched$BETA_GENEHACKMAN != 0 & matched$BETA_OTHER_ALIGNED != 0
  same_sign <- sign(matched$BETA_GENEHACKMAN) == sign(matched$BETA_OTHER_ALIGNED)
  pct <- function(x) if (length(x) == 0) NA_real_ else 100 * mean(x, na.rm = TRUE)
  summary <- data.frame(
    MEASURE = c("variants_genehackman", "variants_other", "variants_matched", "only_genehackman", "only_other",
                "rsid_both_assigned", "rsid_identical_pct", "rsid_only_genehackman", "rsid_only_other",
                "effect_direction_concordant_pct", "aligned_beta_correlation"),
    VALUE = c(nrow(gh), nrow(ot), nrow(matched), sum(!gh$KEY %in% ot$KEY), sum(!ot$KEY %in% gh$KEY),
              sum(both_rsid), pct(same_rsid[both_rsid]),
              sum(!is.na(matched$RSID_GENEHACKMAN) & nzchar(matched$RSID_GENEHACKMAN) & !both_rsid),
              sum(!is.na(matched$RSID_OTHER) & grepl("^rs", matched$RSID_OTHER) & !both_rsid),
              pct(same_sign[nonzero]),
              if (nrow(matched) > 2) stats::cor(matched$BETA_GENEHACKMAN, matched$BETA_OTHER_ALIGNED,
                                                use = "complete.obs") else NA_real_))
  discordant <- matched[(both_rsid & !same_rsid) | (nonzero & !same_sign), ]
  list(summary = summary, discordant = discordant)
}

#' Rewrite a GWAS so variants are identified only by \code{chr:pos:a1:a2} (no CHR, BP or RSID columns).
#'
#' Reproduces the case Reviewer 1 reported (files with only \code{chr:pos:a1:a2} identifiers, to be harmonised
#' to GRCh38 with rsIDs).  The alleles in the identifier are written as \code{OA:EA}, the usual REF:ALT order;
#' EA and OA stay as columns.
#' @param gwas data frame with CHR, BP, EA, OA
#' @return gwas with SNP = \code{chr<CHR>:<BP>:<OA>:<EA>} and without CHR, BP and RSID
#' @export
as_chrpos_identifiers <- function(gwas) {
  gwas$SNP <- paste0("chr", gwas$CHR, ":", format(as.numeric(gwas$BP), scientific = FALSE, trim = TRUE), ":",
                     gwas$OA, ":", gwas$EA)
  keep <- setdiff(names(gwas), c("CHR", "BP", "RSID"))
  if (data.table::is.data.table(gwas)) return(gwas[, keep, with = FALSE])
  gwas[, keep, drop = FALSE]
}
