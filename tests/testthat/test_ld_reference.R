pack_dir <- testthat::test_path("data", "ld_pack_tiny")

read_truth <- function(pop, chr) {
  unname(as.matrix(data.table::fread(file.path(pack_dir, paste0("truth_", pop, "_chr", chr, ".tsv.gz")),
                                     header = FALSE)))
}

# What the pack should return: truth within the stored radius, the dead-zone applied, quantised
expected_pack_matrix <- function(pop, chr) {
  variants <- load_ld_pack_variants(chr, pop, pack_dir)
  truth <- read_truth(pop, chr)
  bp <- variants$BP37
  expected <- truth
  expected[abs(outer(bp, bp, "-")) > 40000] <- 0
  expected[abs(expected) < 0.02] <- 0
  diag(expected) <- 1
  return(expected)
}

test_that("ld_reference.get_ld_matrix_from_pack matches the source LD within quantisation error", {
  for (pop in c("EUR", "CSA")) {
    for (chr in c("21", "22")) {
      variants <- load_ld_pack_variants(chr, pop, pack_dir)
      expected <- expected_pack_matrix(pop, chr)
      snps <- rev(variants$SNP37)
      result <- get_ld_matrix_from_pack(snps, chr, pop, "GRCh37", pack_dir)

      expect_equal(result$snps, snps)
      expect_true(isSymmetric(result$matrix))
      expect_equal(unname(diag(result$matrix)), rep(1, length(snps)))
      idx <- match(result$snps, variants$SNP37)
      # the truth file is rounded to 6 decimals, so skip values sitting on the dead-zone edge
      clear <- abs(abs(read_truth(pop, chr)[idx, idx]) - 0.02) > 1e-5
      expect_lte(max(abs(result$matrix - expected[idx, idx])[clear]), 0.5 / 127 + 1e-6)
    }
  }
})

test_that("ld_reference.get_ld_matrix_from_pack orients r to the canonical effect allele", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  truth <- read_truth("EUR", "22")
  # pick a strongly correlated pair where exactly one variant has EA == REF (sign flipped vs ALT)
  flipped <- variants$EA == variants$REF37
  pairs <- which(abs(truth) > 0.3 & outer(flipped, !flipped) & upper.tri(truth) &
                   abs(outer(variants$BP37, variants$BP37, "-")) <= 40000, arr.ind = TRUE)
  expect_gt(nrow(pairs), 0)
  i <- pairs[1, 1]
  j <- pairs[1, 2]
  result <- get_ld_matrix_from_pack(variants$SNP37[c(i, j)], "22", "EUR", "GRCh37", pack_dir)
  expect_equal(sign(result$matrix[1, 2]), sign(truth[i, j]))
})

test_that("ld_reference.get_ld_matrix_from_pack looks up GRCh38 keys and drops unknown or unmapped variants", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  unmapped <- variants[LIFT_STATUS != "ok"]
  mapped <- variants[LIFT_STATUS == "ok"]
  expect_gt(nrow(unmapped), 0)

  wanted <- c(mapped$SNP38[1:5], "22:1_A_C", unmapped$SNP37[1])
  result <- get_ld_matrix_from_pack(wanted, "22", "EUR", "GRCh38", pack_dir)
  expect_equal(result$snps, mapped$SNP38[1:5])

  result37 <- get_ld_matrix_from_pack(mapped$SNP37[1:5], "22", "EUR", "GRCh37", pack_dir)
  expect_equal(unname(result$matrix), unname(result37$matrix))

  expect_null(get_ld_matrix_from_pack(c("22:1_A_C", mapped$SNP38[1]), "22", "EUR", "GRCh38", pack_dir))
})

test_that("ld_reference.get_ld_matrix_from_pack reports pairs beyond the stored radius", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  near <- get_ld_matrix_from_pack(variants$SNP37[1:3], "22", "EUR", "GRCh37", pack_dir)
  expect_equal(near$beyond_band_fraction, 0)

  far <- c(variants$SNP37[1], variants$SNP37[nrow(variants)])
  far_result <- get_ld_matrix_from_pack(far, "22", "EUR", "GRCh37", pack_dir)
  expect_equal(far_result$beyond_band_fraction, 1)
  expect_equal(far_result$matrix[1, 2], 0)
})

test_that("ld_reference.resolve_ld_ancestry maps SAS to CSA for the pack and validates codes", {
  expect_equal(resolve_ld_ancestry("SAS", "pan_ukb"), "CSA")
  expect_equal(resolve_ld_ancestry("mid", "pan_ukb"), "MID")
  expect_equal(resolve_ld_ancestry("SAS", "1000genomes"), "SAS")
  expect_error(resolve_ld_ancestry("MID", "1000genomes"), "not available")
  expect_error(resolve_ld_ancestry("EUR", "gnomad"), "ld_reference must be one of")

  sas <- get_ld_matrix_from_pack(load_ld_pack_variants("21", "CSA", pack_dir)$SNP37[1:4], "21", "SAS",
                                 "GRCh37", pack_dir)
  expect_equal(nrow(sas$matrix), 4)
  expect_equal(ld_pack_n_samples("SAS", pack_dir), 120)
})

test_that("ld_reference.make_ld_positive_semidefinite clips negative eigenvalues only when needed", {
  ok <- make_ld_positive_semidefinite(diag(3))
  expect_equal(ok$matrix, diag(3))
  expect_equal(ok$min_eigenvalue, 1)

  not_psd <- matrix(c(1, 0.9, 0, 0.9, 1, 0.9, 0, 0.9, 1), 3)
  fixed <- make_ld_positive_semidefinite(not_psd)
  expect_lt(fixed$min_eigenvalue, 0)
  expect_gte(min(eigen(fixed$matrix, symmetric = TRUE, only.values = TRUE)$values), -1e-9)
  expect_equal(diag(fixed$matrix), rep(1, 3))
  # clipping changes R much less than adding |lambda_min| to the diagonal would
  ridge <- not_psd + diag(abs(fixed$min_eigenvalue), 3)
  expect_lt(max(abs(fixed$matrix - not_psd)), max(abs(stats::cov2cor(ridge) - not_psd)))
})

test_that("ld_reference.get_ld_matrix dispatches to PLINK for 1000genomes", {
  local_mocked_bindings(
    compute_ld_matrix = function(rsids, chr, ancestry) list(matrix = diag(length(rsids)), snps = rsids,
                                                             ancestry = ancestry)
  )
  result <- get_ld_matrix(c("rs1", "rs2"), 1, "EUR", ld_reference = "1000genomes")
  expect_equal(result$snps, c("rs1", "rs2"))
  expect_equal(result$ancestry, "EUR")
  expect_error(get_ld_matrix("rs1", 1, "EUR", ld_reference = "other"), "ld_reference must be one of")
})

test_that("ld_reference.read_ld_pack_info gives a clear error when the pack is missing", {
  expect_error(read_ld_pack_info("EUR", tempfile("no_pack_")), "reference_data/pan_ukb_ld/README.md")
})
