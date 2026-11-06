test_that("compare_harmonisation matches on the allele pair and aligns effects to GeneHackman's effect allele", {
  genehackman <- data.frame(CHR = 1, BP = c(100, 200, 300, 400), EA = c("A", "C", "G", "A"), OA = c("G", "T", "T", "C"),
                            BETA = c(0.1, -0.2, 0.3, 0.05), RSID = c("rs1", "rs2", "rs3", "rs4"))
  # MungeSumstats-like: A2 is the effect allele; variant 2 reported for the other allele; variant 3 sign wrong;
  # variant 4 missing; one variant only in MungeSumstats
  munged <- data.frame(CHR = 1, BP = c(100, 200, 300, 500), A1 = c("G", "C", "T", "A"), A2 = c("A", "T", "G", "C"),
                       BETA = c(0.1, 0.2, -0.3, 0.1), SNP = c("rs1", "rs2", "rs99", "rs5"))
  result <- compare_harmonisation(genehackman, munged)
  value <- stats::setNames(result$summary$VALUE, result$summary$MEASURE)
  expect_equal(unname(value[c("variants_matched", "only_genehackman", "only_other")]), c(3, 1, 1))
  expect_equal(unname(value["rsid_identical_pct"]), 200 / 3)
  expect_equal(unname(value["effect_direction_concordant_pct"]), 200 / 3)
  expect_equal(sort(result$discordant$KEY), "1:300_G_T")
})

test_that("as_chrpos_identifiers reproduces the reviewer's input, and standardisation recovers CHR, BP and alleles", {
  gwas <- data.frame(CHR = "1", BP = 1000, EA = "G", OA = "A", EAF = 0.8, BETA = -0.1, SE = 0.02, P = 1e-6,
                     RSID = "rs1")
  chrpos <- as_chrpos_identifiers(gwas)
  expect_equal(chrpos$SNP, "chr1:1000:A:G")
  expect_false(any(c("CHR", "BP", "RSID") %in% names(chrpos)))
  out <- suppressWarnings(suppressMessages(standardise_gwas(
    chrpos, N = 1000, input_reference_build = "GRCh38", output_reference_build = "GRCh38",
    populate_rsid_option = populate_rsid_options$none, check_build = FALSE, strand_check = FALSE)))
  expect_equal(c(out$SNP, out$EA, out$BETA), c("1:1000_A_G", "A", "0.1"))
})
