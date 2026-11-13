liftover_test_gwas <- function() {
  data.frame(
    SNP = paste0("v", 1:6),
    CHR = "1",
    BP = c(100000, 200000, 300000, 400000, 500000, 600000),
    EA = "A",
    OA = c("G", "G", "G", "T", "G", "G"),
    BETA = 0.1, SE = 0.01, P = 1e-3,
    stringsAsFactors = FALSE
  )
}

test_that("liftover.convert_reference_build_via_liftover drops untrustworthy variants with a reason", {
  local_mocked_bindings(run_liftover = mock_liftover_shift(
    shift = 1000,
    overrides = list("2" = "unmapped", "3" = list(chrom = "chr1_KI270706v1_random"), "4" = list(strand = "-"),
                     "5" = list(chrom = "chr21")),
    reverse_overrides = list("6" = list(end = 123))
  ))
  local_mocked_bindings(liftover_chain_available = function(...) TRUE)
  prefix <- file.path(withr::local_tempdir(), "lifted")

  result <- convert_reference_build_via_liftover(liftover_test_gwas(), reference_builds$GRCh37,
                                                 reference_builds$GRCh38, max_drop_fraction = 1,
                                                 report_prefix = prefix)

  expect_equal(result$SNP, "v1")
  expect_equal(result$BP, 101000)
  expect_equal(result$BP37, 100000)
  expect_false(".ROW_ID" %in% colnames(result))

  qc <- attr(result, "liftover_qc")
  counts <- stats::setNames(qc$N, qc$REASON)
  expect_equal(unname(counts[c("input", "kept", "unmapped", "non_primary_contig", "strand_flip", "chrom_change",
                               "round_trip", "duplicate_target")]), c(6, 1, 1, 1, 1, 1, 1, 0))

  dropped <- data.table::fread(paste0(prefix, ".liftover_dropped.tsv.gz"))
  expect_setequal(dropped$LIFTOVER_DROP_REASON,
                  c("unmapped", "non_primary_contig", "strand_flip", "chrom_change", "round_trip"))
  expect_true(file.exists(paste0(prefix, ".liftover_qc.tsv")))
})

test_that("liftover.convert_reference_build_via_liftover drops different positions that collide on one variant", {
  gwas <- liftover_test_gwas()[1:2, ]
  local_mocked_bindings(run_liftover = mock_liftover_shift(
    shift = 0, overrides = list("1" = list(end = 150000), "2" = list(end = 150000)),
    reverse_overrides = list("1" = list(end = 100000), "2" = list(end = 200000))
  ))
  local_mocked_bindings(liftover_chain_available = function(...) TRUE)
  result <- convert_reference_build_via_liftover(gwas, "GRCh37", "GRCh38", max_drop_fraction = 1)
  expect_equal(nrow(result), 0)
  qc <- attr(result, "liftover_qc")
  expect_equal(qc$N[qc$REASON == "duplicate_target"], 2)
})

test_that("liftover.convert_reference_build_via_liftover complements non-palindromic minus-strand variants on request", {
  gwas <- liftover_test_gwas()[1:2, ]
  gwas$OA[2] <- "T"  # v2 is A/T: palindromic
  flip_all <- list("1" = list(strand = "-"), "2" = list(strand = "-"))
  local_mocked_bindings(run_liftover = mock_liftover_shift(shift = 10, overrides = flip_all))
  local_mocked_bindings(liftover_chain_available = function(...) TRUE)

  result <- convert_reference_build_via_liftover(gwas, "GRCh37", "GRCh38", complement_minus_strand = TRUE,
                                                 max_drop_fraction = 1)
  expect_equal(result$SNP, "v1")
  expect_equal(c(result$EA, result$OA), c("T", "C"))
})

test_that("liftover.convert_reference_build_via_liftover stops when too many variants are dropped", {
  local_mocked_bindings(run_liftover = mock_liftover_shift(overrides = list("1" = "unmapped", "2" = "unmapped")))
  local_mocked_bindings(liftover_chain_available = function(...) TRUE)
  expect_error(convert_reference_build_via_liftover(liftover_test_gwas(), "GRCh37", "GRCh38"),
               "check the declared build")
})

test_that("liftover.convert_reference_build_via_liftover returns the input when builds match", {
  gwas <- liftover_test_gwas()
  expect_identical(convert_reference_build_via_liftover(gwas, "GRCh38", "GRCh38"), gwas)
})

test_that("liftover.convert_reference_build_via_liftover handles real chain edge cases", {
  local_real_liftover()
  # GRCh37 positions with known GRCh38 behaviour under UCSC hg19ToHg38
  gwas <- data.frame(
    SNP = c("round", "plus", "inverted", "alt_contig", "other_chrom", "unstable"),
    CHR = "1",
    BP = c(20000000, 100000000, 144500000, 143400000, 143151070, 1583243),
    EA = "A", OA = "G", BETA = 0.1, SE = 0.01, P = 0.01,
    stringsAsFactors = FALSE
  )
  result <- convert_reference_build_via_liftover(gwas, "GRCh37", "GRCh38", max_drop_fraction = 1)

  expect_equal(result$SNP, c("round", "plus"))
  expect_equal(result$BP, c(19673507, 99534444))
  qc <- attr(result, "liftover_qc")
  counts <- stats::setNames(qc$N, qc$REASON)
  expect_equal(unname(counts[c("strand_flip", "non_primary_contig", "chrom_change", "round_trip")]), c(1, 1, 1, 1))

  back <- convert_reference_build_via_liftover(result, "GRCh38", "GRCh37", max_drop_fraction = 1)
  expect_equal(back$BP, c(20000000, 100000000))
})

test_that("liftover.check_declared_build stops when positions match the other build", {
  pack_dir <- withr::local_tempdir()
  dir.create(file.path(pack_dir, "EUR"))
  bp <- format(1:200 * 1000, scientific = FALSE, trim = TRUE)
  bp38 <- format(1:200 * 1000 + 7, scientific = FALSE, trim = TRUE)
  pack <- data.frame(SNP37 = paste0("1:", bp, "_A_G"), SNP38 = paste0("1:", bp38, "_A_G"))
  data.table::fwrite(pack, file.path(pack_dir, "EUR", "variants_chr1.tsv.gz"), sep = "\t")
  gwas37 <- data.frame(CHR = "1", BP = 1:200 * 1000, EA = "G", OA = "A")

  share <- check_declared_build(gwas37, "GRCh37", pack_dir)
  expect_equal(share[["GRCh37"]], 1)
  expect_error(check_declared_build(gwas37, "GRCh38", pack_dir), "Set build: GRCh37")
  expect_message(check_declared_build(gwas37, "GRCh37", file.path(pack_dir, "missing")), "Build check skipped")
})
