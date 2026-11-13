# Known-answer harmonisation tests: each case has an input row and the exact output expected from
# standardise_gwas() (no liftover, no RSID lookup), so the tests check identifiers, alleles, effect sign,
# allele frequency and P rather than only that the code runs.

standardise_known <- function(gwas, ...) {
  suppressWarnings(suppressMessages(standardise_gwas(
    gwas, N = 1000, input_reference_build = "GRCh38", output_reference_build = "GRCh38",
    populate_rsid_option = populate_rsid_options$none, check_build = FALSE, strand_check = FALSE, ...)))
}

base_row <- function(...) {
  row <- list(CHR = "1", BP = 1000, EA = "A", OA = "G", EAF = 0.2, BETA = 0.1, SE = 0.02, P = 1e-6)
  utils::modifyList(row, list(...)) |> as.data.frame()
}

test_that("alleles already in canonical order are unchanged", {
  out <- standardise_known(base_row())
  expect_equal(out$SNP, "1:1000_A_G")
  expect_equal(c(out$EA, out$OA), c("A", "G"))
  expect_equal(c(out$BETA, out$EAF), c(0.1, 0.2))
})

test_that("swapped effect and other alleles flip BETA and EAF, and give the same key", {
  out <- standardise_known(base_row(EA = "G", OA = "A", BETA = -0.1, EAF = 0.8))
  expect_equal(out$SNP, "1:1000_A_G")
  expect_equal(c(out$EA, out$OA), c("A", "G"))
  expect_equal(c(out$BETA, out$EAF), c(0.1, 0.2))
})

test_that("flip_alleles = FALSE returns the input orientation", {
  out <- standardise_known(base_row(EA = "G", OA = "A", BETA = -0.1, EAF = 0.8), flip_alleles = FALSE)
  expect_equal(c(out$EA, out$OA, out$SNP), c("G", "A", "1:1000_G_A"))
  expect_equal(c(out$BETA, out$EAF), c(-0.1, 0.8))
})

test_that("lower-case alleles and chr-prefixed chromosomes are normalised", {
  out <- standardise_known(base_row(CHR = "chr1", EA = "g", OA = "a", BETA = -0.1, EAF = 0.8))
  expect_equal(out$SNP, "1:1000_A_G")
  expect_equal(out$BETA, 0.1)
  x <- standardise_known(base_row(CHR = "chrX"))
  expect_equal(x$SNP, "X:1000_A_G")
})

test_that("indels: multi-base alleles and I/D coding have one key whichever allele is the effect allele", {
  insertion <- standardise_known(base_row(EA = "AT", OA = "A", BETA = 0.1, EAF = 0.3))
  expect_equal(c(insertion$SNP, insertion$BETA, insertion$EAF), c("1:1000_A_AT", "-0.1", "0.7"))
  i_d <- standardise_known(base_row(EA = "I", OA = "D", BETA = 0.1, EAF = 0.3))
  d_i <- standardise_known(base_row(EA = "D", OA = "I", BETA = -0.1, EAF = 0.7))
  expect_equal(i_d$SNP, d_i$SNP)
  expect_equal(c(i_d$BETA, i_d$EAF), c(d_i$BETA, d_i$EAF))
})

test_that("odds ratios with confidence intervals give log-OR and its SE", {
  gwas <- base_row(OR = 1.5, OR_LB = 1.2, OR_UB = 1.875)
  gwas$BETA <- NULL; gwas$SE <- NULL
  out <- standardise_known(gwas)
  expect_equal(out$BETA, log(1.5))
  expect_equal(out$SE, (log(1.5) - log(1.2)) / qnorm(0.975))
})

test_that("Z-only input gives BETA / SE equal to Z", {
  gwas <- base_row(Z = -4)
  gwas$BETA <- NULL; gwas$SE <- NULL; gwas$P <- NULL
  out <- standardise_known(gwas)
  expect_equal(out$BETA / out$SE, -4)
  expect_equal(out$P, 2 * pnorm(-4))
})

test_that("-log10 P is converted and P = 0 is kept as the smallest double", {
  gwas <- base_row(LOG_P = 10); gwas$P <- NULL
  expect_equal(standardise_known(gwas)$P, 1e-10)
  expect_equal(standardise_known(base_row(P = 0))$P, .Machine$double.xmin)
})

test_that("chr:pos:a1:a2 identifiers fill CHR and BP (alleles come from EA/OA)", {
  gwas <- base_row(SNP = "chr12:345:G:A", EA = "G", OA = "A", BETA = -0.1, EAF = 0.8)
  gwas$CHR <- NULL; gwas$BP <- NULL
  out <- standardise_known(gwas)
  expect_equal(out$SNP, "12:345_A_G")
  expect_equal(out$BETA, 0.1)
})

test_that("the per-variant N in the file is kept over the input N", {
  out <- standardise_known(base_row(N = 750))
  expect_equal(out$N, 750)
})

test_that("presets map effect alleles correctly: OpenGWAS gives effects for ALT, GWAS Catalog for effect_allele", {
  opengwas <- resolve_column_map("opengwas")
  expect_equal(c(opengwas$EA, opengwas$OA), c("ALT", "REF"))
  catalog <- resolve_column_map("gwas_catalog")
  expect_equal(c(catalog$EA, catalog$OA, catalog$BP), c("effect_allele", "other_allele", "base_pair_location"))
  ssf <- data.frame(chromosome = 3, base_pair_location = 500, effect_allele = "T", other_allele = "C",
                    effect_allele_frequency = 0.9, beta = 0.2, standard_error = 0.05, p_value = 1e-4, n = 5000)
  out <- standardise_known(ssf, input_columns = "gwas_catalog")
  expect_equal(c(out$SNP, out$BETA, out$EAF), c("3:500_C_T", "-0.2", "0.1"))
})

test_that("strand: complemented non-palindromic SNVs are fixed, palindromic SNVs resolved by frequency", {
  reference <- data.frame(SNP = c("1:100_A_G", "1:200_A_T", "1:300_C_G", "1:400_A_C", "1:500_A_G"),
                          REF_EAF = c(0.3, 0.2, 0.45, 0.6, 0.1))
  gwas <- data.frame(SNP = c("1:100_C_T", "1:200_A_T", "1:300_C_G", "1:400_A_C", "1:500_A_G", "1:600_A_G"),
                     CHR = 1, BP = c(100, 200, 300, 400, 500, 600),
                     EA = c("C", "A", "C", "A", "A", "A"), OA = c("T", "T", "G", "C", "G", "G"),
                     BETA = rep(0.1, 6), EAF = c(0.7, 0.8, 0.5, 0.58, 0.6, 0.3))
  checked <- harmonise_against_reference(gwas, reference)
  out <- checked$gwas
  # C/T on the minus strand is G/A: EA C -> G, so the canonical key is A_G with EA = A (BETA and EAF flip)
  expect_equal(out$SNP[1], "1:100_A_G")
  expect_equal(c(out$BETA[1], out$EAF[1]), c(-0.1, 0.3))
  expect_equal(out$ALLELE_QC, c("strand_flipped", "palindromic_flipped", "palindromic_ambiguous", "match",
                                "eaf_discordant", "not_in_reference"))
  expect_equal(c(out$BETA[2], out$EAF[2]), c(-0.1, 0.2))
  expect_equal(out$BETA[3:6], rep(0.1, 4))
  expect_equal(checked$report$N[checked$report$STATUS == "palindromic_ambiguous"], 1)
  dropped <- harmonise_against_reference(gwas, reference, drop_ambiguous_palindromic = TRUE)$gwas
  expect_false("1:300_C_G" %in% dropped$SNP)
})
