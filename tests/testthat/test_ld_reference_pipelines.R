# Pipelines using the Pan-UKB LD pack (ld_reference = "pan_ukb"), against the tiny fixture pack.
fixture_pack <- testthat::test_path("data", "ld_pack_tiny")

# Point the package at the fixture pack for the duration of `code`
with_fixture_pack <- function(code) {
  old <- options(GeneHackman.ld_pack_dir = fixture_pack)
  on.exit(options(old))
  force(code)
}

# GWAS on the fixture's chr22 variants with a strong signal at one variant, built from the pack's LD
pack_gwas <- function(build = "GRCh37", rsid = FALSE, causal = 20) {
  variants <- load_ld_pack_variants("22", "EUR", fixture_pack)
  variants <- variants[LIFT_STATUS == "ok"]
  R <- get_ld_matrix_from_pack(variants$SNP37, "22", "EUR", "GRCh37", fixture_pack)$matrix
  z <- R[, causal] * 12
  gwas <- data.table::data.table(
    SNP = if (build == "GRCh37") variants$SNP37 else variants$SNP38,
    CHR = 22, BP = if (build == "GRCh37") variants$BP37 else variants$BP38,
    EA = variants$EA, OA = variants$OA, EAF = variants$EAF,
    BETA = z * 0.01, SE = 0.01, P = 2 * stats::pnorm(-abs(z)), N = 50000
  )
  if (rsid) gwas$RSID <- paste0("rs", seq_len(nrow(gwas)))
  gwas
}

test_that("ld_reference.finemap_gwas fine-maps with the LD pack, matching on SNP keys without RSIDs", {
  with_fixture_pack({
    for (build in c("GRCh37", "GRCh38")) {
      gwas <- pack_gwas(build)
      gwas_file <- tempfile(fileext = ".tsv.gz")
      data.table::fwrite(gwas, gwas_file, sep = "\t")
      lead <- gwas[which.min(P)]
      clump_file <- tempfile(fileext = ".clumped")
      writeLines(c("CHR F SNP BP P", paste(22, 1, lead$SNP, lead$BP, lead$P)), clump_file)
      out_dir <- tempfile("finemap_pack_")

      res <- finemap_gwas(gwas_file, clump_file, ancestry = "EUR", default_n = 50000, output_finemap_dir = out_dir,
                          window_kb = 30, max_causal = 3, ld_reference = "pan_ukb", build = build)

      regions <- data.table::fread(file.path(out_dir, "finemap_regions.tsv"))
      expect_equal(regions$STATUS, "fine-mapped")
      expect_equal(regions$LD_REFERENCE, "pan_ukb")
      expect_equal(regions$N_REF, 400)
      expect_gte(regions$FRAC_BEYOND_LD_RADIUS, 0)
      expect_false(is.na(regions$LD_RIDGE))
      expect_gt(regions$N_SNPS, 10)

      locus <- data.table::fread(list.files(out_dir, pattern = "_finemap.tsv.gz$", full.names = TRUE))
      expect_equal(nrow(locus), sum(abs(gwas$BP - lead$BP) <= 30000))
      expect_true(all(c("PIP", "CS") %in% names(locus)))
      expect_equal(locus$SNP[which.max(locus$PIP)], gwas$SNP[20])
      expect_gt(max(locus$PIP, na.rm = TRUE), 0.5)
    }
  })
})

test_that("ld_reference.finemap_gwas with 1000genomes still requires RSIDs", {
  gwas <- pack_gwas()
  gwas_file <- tempfile(fileext = ".tsv.gz")
  data.table::fwrite(gwas, gwas_file, sep = "\t")
  clump_file <- tempfile(fileext = ".clumped")
  writeLines(c("CHR F SNP BP P", paste(22, 1, gwas$SNP[20], gwas$BP[20], gwas$P[20])), clump_file)
  out_dir <- tempfile("finemap_1kg_")
  finemap_gwas(gwas_file, clump_file, ancestry = "EUR", default_n = 50000, output_finemap_dir = out_dir,
               window_kb = 30)
  regions <- data.table::fread(file.path(out_dir, "finemap_regions.tsv"))
  expect_equal(regions$STATUS, "skipped: no RSIDs for LD computation")
  expect_equal(regions$LD_REFERENCE, "1000genomes")
})

test_that("ld_reference.populate_eaf_from_reference_panel fills EAF from the LD pack by SNP key", {
  with_fixture_pack({
    gwas <- pack_gwas("GRCh38")
    expected <- gwas$EAF
    gwas$EAF[1:10] <- NA
    gwas <- rbind(gwas, data.table::data.table(SNP = "22:1_A_C", CHR = 22, BP = 1, EA = "A", OA = "C", EAF = NA,
                                               BETA = 0, SE = 1, P = 1, N = 1))
    filled <- populate_eaf_from_reference_panel(gwas, "EUR", ld_reference = "pan_ukb", build = "GRCh38")
    expect_equal(filled$EAF[1:10], expected[1:10])
    expect_true(is.na(filled$EAF[nrow(filled)]))
  })
})

test_that("ld_reference.standardise_gwas fills EAF from the LD pack without populating RSIDs", {
  with_fixture_pack({
    gwas <- pack_gwas()[1:15]
    expected <- gwas$EAF
    gwas$EAF <- NULL
    standardised <- standardise_gwas(gwas, populate_eaf = TRUE, ancestry = "EUR", ld_reference = "pan_ukb")
    expect_false("RSID" %in% names(standardised))
    expect_equal(standardised$EAF[match(gwas$SNP, standardised$SNP)], expected)
  })
})

test_that("ld_reference.LocusZoom helpers use pack keys for the index variant and r2", {
  with_fixture_pack({
    gwas <- pack_gwas()
    index <- find_index_snp_from_ld_panel(22, gwas$BP[20] + 10, 5000, "EUR", "pan_ukb", "GRCh37")
    expect_equal(index, gwas$SNP[20])
    r2 <- compute_ld_r2_for_locus(gwas$SNP[15:25], 22, "EUR", index, "pan_ukb", "GRCh37")
    expect_equal(unname(r2[index]), 1)
    expected <- get_ld_matrix_from_pack(gwas$SNP[15:25], "22", "EUR", "GRCh37", fixture_pack)$matrix[6, ]^2
    expect_equal(unname(r2[gwas$SNP[15:25]]), unname(expected))
  })
})
