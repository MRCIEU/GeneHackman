# Sample size (N, N_EFF), Z-only input, fine-mapping options (L rule, region strategies, MHC flag),
# user-supplied LD and LD-mismatch diagnostics.

test_that("a per-variant N in the file takes precedence over the input N", {
  gwas <- data.frame(SNP = c("a", "b", "c"), N = c(1000, NA, 0))
  filled <- suppressMessages(fill_sample_size(gwas, 5000))
  expect_equal(filled$N, c(1000, 5000, 5000))
  expect_equal(fill_sample_size(data.frame(SNP = "a"), 5000)$N, 5000)
  expect_equal(fill_sample_size(gwas, 0)$N, gwas$N)
})

test_that("effective sample size follows 4 / (1/cases + 1/controls)", {
  expect_equal(effective_sample_size(1000, 1000), 2000)
  expect_equal(effective_sample_size(3841, 1e6), 4 / (1 / 3841 + 1 / 1e6))
  expect_true(is.na(effective_sample_size(NA, 10)))
  expect_true(is.na(effective_sample_size(0, 10)))
})

test_that("populate_effective_n uses, in order: N_EFF, N_CASES/N_CONTROLS, N scaled by case fraction, study counts", {
  expect_equal(populate_effective_n(data.frame(N = 10), "continuous", 1, 1), data.frame(N = 10))
  kept <- populate_effective_n(data.frame(N_EFF = c(5, 6), N = 100), "binary", 10, 90)
  expect_equal(kept$N_EFF, c(5, 6))
  per_variant <- populate_effective_n(data.frame(N_CASES = c(100, 50), N_CONTROLS = c(900, 950)), "binary")
  expect_equal(per_variant$N_EFF, effective_sample_size(c(100, 50), c(900, 950)))
  scaled <- populate_effective_n(data.frame(N = c(1000, 500)), "binary", n_case = 100, n_control = 900)
  expect_equal(scaled$N_EFF, c(1000, 500) * 4 * 0.1 * 0.9)
  study <- populate_effective_n(data.frame(BETA = 1:2), "case_control", n_case = 100, n_control = 900)
  expect_equal(study$N_EFF, rep(effective_sample_size(100, 900), 2))
  expect_warning(populate_effective_n(data.frame(N = 1), "binary"), "N_EFF can't be computed")
  expect_error(normalise_study_type("ordinal"), "continuous or binary")
})

test_that("populate_effective_n derives N_CONTROLS = N - N_CASES when the file gives cases and the total", {
  gwas <- data.frame(N_CASES = c(100, 50, NA, 10), N = c(1000, 1000, 1000, 5))
  result <- populate_effective_n(gwas, "binary")
  expect_equal(result$N_CONTROLS, c(900, 950, NA, NA))
  expect_equal(result$N_EFF, c(effective_sample_size(c(100, 50), c(900, 950)), NA, NA))

  # variants without per-variant counts get the study's N_EFF
  filled <- populate_effective_n(gwas, "binary", n_case = 200, n_control = 1800)
  expect_equal(filled$N_EFF[3:4], rep(effective_sample_size(200, 1800), 2))

  # an N_CONTROLS column in the file is kept as it is
  given <- populate_effective_n(data.frame(N_CASES = 100, N_CONTROLS = 300, N = 1000), "binary")
  expect_equal(given$N_CONTROLS, 300)
  expect_equal(given$N_EFF, effective_sample_size(100, 300))
})

test_that("standardise_gwas takes the per-variant totals from the file, not the input's N, for N_CONTROLS", {
  # derive_n_controls runs before the input's N fills a missing N column, so no controls are invented
  gwas <- data.frame(N_CASES = c(100, 50))
  expect_false("N_CONTROLS" %in% colnames(derive_n_controls(gwas)))
  expect_equal(derive_n_controls(data.frame(N_CASES = 100, N = "1000"))$N_CONTROLS, 900)
})

test_that("analysis_sample_size prefers the median N_EFF, then the median N, then the default", {
  expect_equal(analysis_sample_size(data.frame(N = c(10, 20, 30), N_EFF = c(4, 5, NA))), 4.5)
  expect_equal(analysis_sample_size(data.frame(N = c(10, 20, 30))), 20)
  expect_equal(analysis_sample_size(data.frame(BETA = 1), default_n = 7), 7)
})

test_that("Z-only input is converted to standardised BETA and SE (known answer)", {
  gwas <- data.frame(Z = c(5, -2), EAF = c(0.3, 0.5), N = c(10000, 10000))
  converted <- convert_z_score_to_beta(gwas)
  se <- 1 / sqrt(2 * c(0.3, 0.5) * c(0.7, 0.5) * (10000 + c(25, 4)))
  expect_equal(converted$SE, se)
  expect_equal(converted$BETA, c(5, -2) * se)
  expect_equal(converted$BETA / converted$SE, c(5, -2))
  expect_equal(converted$P, 2 * pnorm(-abs(c(5, -2))))
  expect_error(convert_z_score_to_beta(data.frame(Z = 1)), "EAF")
  expect_error(convert_z_score_to_beta(data.frame(Z = 1, EAF = 0.2)), "sample size")
})

test_that("max_causal is used as given, and a small LD reference gives a warning, not a different L", {
  expect_equal(resolve_max_causal(10), 10L)
  expect_equal(resolve_max_causal("5"), 5L)
  expect_error(resolve_max_causal("auto"), "positive whole number")
  expect_error(resolve_max_causal(0), "positive whole number")
  expect_error(resolve_max_causal(2.5), "positive whole number")
  expect_match(ld_reference_size_warning(1.7e6, 503), "3,380 times that of the LD reference \\(503 individuals\\)")
  expect_true(is.na(ld_reference_size_warning(4000, 503)))
  expect_true(is.na(ld_reference_size_warning(4000, NA)))
  expect_match(ld_reference_size_warning(1e5, 503, label = "EUR"), "^EUR: ")
  expect_equal(ld_size_warning_markdown(c(NA, "")), "")
  expect_match(ld_size_warning_markdown(c("A.", NA, "A.")), "Warning: LD reference much smaller than the GWAS.\\*\\* A\\.\n")
})

test_that("region strategy none keeps one window per lead", {
  leads <- data.frame(SNP = c("a", "b"), CHR = 1, BP = c(1e6, 1.5e6), P = c(1e-9, 1e-10))
  merged <- define_finemap_regions(leads, window_bp = 5e5, max_region_bp = 3e6)
  expect_equal(nrow(merged), 1)
  separate <- define_finemap_regions(leads, window_bp = 5e5, max_region_bp = 3e6, strategy = "none")
  expect_equal(nrow(separate), 2)
  expect_equal(separate$START, c(5e5, 1e6))
  expect_true(all(separate$REGION_STRATEGY == "none"))
})

test_that("region strategy tile covers a long cluster with overlapping tiles and assigns each position once", {
  leads <- data.frame(SNP = paste0("s", 1:6), CHR = 2, BP = seq(2e6, 7e6, by = 1e6), P = 1e-9)
  tiles <- define_finemap_regions(leads, window_bp = 1e6, max_region_bp = 3e6, strategy = "tile")
  expect_true(nrow(tiles) > 1)
  expect_true(all(tiles$END - tiles$START <= 3e6))
  expect_equal(length(unique(tiles$LOCUS)), nrow(tiles))
  expect_true(all(unlist(strsplit(tiles$LEADS, ",")) %in% leads$SNP))
  positions <- seq(1.5e6, 7.5e6, by = 1e5)
  owners <- sapply(seq_len(nrow(tiles)), function(i) {
    inside <- positions >= tiles$START[i] & positions <= tiles$END[i]
    tile_core_flags(positions, tiles, i) & inside
  })
  covered <- sapply(positions, function(p) any(tiles$START <= p & tiles$END >= p))
  expect_true(all(rowSums(owners)[covered] == 1))
  short <- define_finemap_regions(leads[1:2, ], window_bp = 5e5, max_region_bp = 3e6, strategy = "tile")
  expect_equal(nrow(short), 1)
  expect_true(all(tile_core_flags(c(1, 2), short, 1)))
})

test_that("regions overlapping the MHC are flagged in the output build", {
  leads <- data.frame(SNP = c("mhc", "other"), CHR = c(6, 6), BP = c(32e6, 50e6), P = 1e-9)
  regions <- suppressMessages(define_finemap_regions(leads, window_bp = 5e5, build = "GRCh38"))
  expect_equal(regions$MHC, c(TRUE, FALSE))
})

test_that("align_ld_to_gwas flips r for swapped alleles and drops mismatched ones", {
  R <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3)
  gwas <- data.frame(EA = c("A", "C", "G"), OA = c("G", "T", "T"))
  aligned <- align_ld_to_gwas(R, gwas, a1 = c("A", "T", "C"), a2 = c("G", "C", "A"))
  expect_equal(aligned$keep, c(TRUE, TRUE, FALSE))
  expect_equal(aligned$n_mismatch, 1)
  expect_equal(aligned$matrix, matrix(c(1, -0.5, -0.5, 1), 2))
})

simulate_genotypes <- function(n = 600, p = 12, seed = 1) {
  set.seed(seed)
  haplotype <- function() {
    h <- matrix(0L, n, p)
    h[, 1] <- rbinom(n, 1, 0.4)
    for (j in 2:p) h[, j] <- ifelse(runif(n) < 0.8, h[, j - 1], rbinom(n, 1, 0.4))
    h
  }
  haplotype() + haplotype()
}

test_that("ld_mismatch_qc flags a variant whose allele was flipped (positive control)", {
  genotypes <- simulate_genotypes(p = 40, seed = 2)
  R <- cor(genotypes)
  n <- 50000
  b <- rep(0, ncol(R)); b[4] <- 0.15
  z <- as.vector(R %*% b) * sqrt(n) * sd(genotypes[, 4])
  consistent <- ld_mismatch_qc(z, R, n)
  expect_false(any(consistent$outlier))
  flipped <- z; flipped[5] <- -flipped[5]
  qc <- ld_mismatch_qc(flipped, R, n)
  expect_true(qc$outlier[5])
  expect_equal(which(qc$outlier), 5L)
  expect_lt(consistent$s, 0.1)
  expect_gt(qc$s, 0.5)
})

write_test_bfile <- function(genotypes, prefix, chr = 1, a1 = NULL, a2 = NULL) {
  p <- ncol(genotypes)
  a1 <- a1 %||% rep("A", p)
  a2 <- a2 %||% rep("G", p)
  map <- data.frame(chr, paste0("v", seq_len(p)), 0, seq_len(p) * 1000)
  ped_alleles <- do.call(cbind, lapply(seq_len(p), function(j) {
    g <- genotypes[, j]
    cbind(ifelse(g >= 1, a1[j], a2[j]), ifelse(g == 2, a1[j], a2[j]))
  }))
  ped <- cbind(paste0("f", seq_len(nrow(genotypes))), paste0("i", seq_len(nrow(genotypes))), 0, 0, 0, -9, ped_alleles)
  utils::write.table(map, paste0(prefix, ".map"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(ped, paste0(prefix, ".ped"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  system2("plink1.9", c("--file", prefix, "--make-bed", "--keep-allele-order", "--out", prefix), stdout = FALSE, stderr = FALSE)
  # make A1 the first listed allele for every variant, as given
  bim <- utils::read.table(paste0(prefix, ".bim"), colClasses = "character")
  invisible(bim)
}

test_that("user-supplied PLINK LD is matched on position and alleles and oriented to the GWAS effect allele", {
  skip_if(!nzchar(Sys.which("plink1.9")), "plink1.9 not installed")
  genotypes <- simulate_genotypes(n = 300, p = 4)
  prefix <- file.path(withr::local_tempdir(), "user_ld")
  write_test_bfile(genotypes, prefix, a1 = c("A", "A", "C", "G"), a2 = c("G", "T", "T", "T"))
  bim <- utils::read.table(paste0(prefix, ".bim"), colClasses = "character")
  # GWAS keys in GeneHackman's canonical order (alphabetical EA_OA)
  keys <- paste0("1:", bim$V4, "_", pmin(bim$V5, bim$V6), "_", pmax(bim$V5, bim$V6))
  ld <- compute_ld_matrix_from_bfile(c(keys, "1:999999_A_C"), 1, prefix)
  expect_equal(ld$snps, keys)
  expected <- cor(genotypes)
  # the panel's A1 is the minor (counted) allele; flip where the GWAS EA is the other one
  counted <- bim$V5
  ea <- pmin(bim$V5, bim$V6)
  dosage_ea <- sapply(seq_len(4), function(j) {
    g_a1 <- if (counted[j] == c("A", "A", "C", "G")[j]) genotypes[, j] else 2 - genotypes[, j]
    if (ea[j] == counted[j]) g_a1 else 2 - g_a1
  })
  expect_equal(unname(ld$matrix), unname(cor(dosage_ea)), tolerance = 1e-4)
  expect_null(compute_ld_matrix_from_bfile(keys[1], 1, prefix))
})

test_that("ld_reference_n_samples counts the user-supplied .fam", {
  prefix <- file.path(withr::local_tempdir(), "x")
  writeLines(c("a a 0 0 0 -9", "b b 0 0 0 -9"), paste0(prefix, ".fam"))
  expect_equal(ld_reference_n_samples("1000genomes", "EUR", bfile = prefix), 2)
  expect_equal(ld_reference_id_column("1000genomes", prefix), "SNP")
  expect_equal(ld_reference_id_column("1000genomes"), "RSID")
})

test_that("finemap_gwas runs with user-supplied LD, warns about the small reference and reports LD diagnostics", {
  skip_if(!nzchar(Sys.which("plink1.9")), "plink1.9 not installed")
  genotypes <- simulate_genotypes(n = 400, p = 30, seed = 3)
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "in_sample")
  write_test_bfile(genotypes, prefix)
  bim <- utils::read.table(paste0(prefix, ".bim"), colClasses = "character")
  ea <- pmin(bim$V5, bim$V6); oa <- pmax(bim$V5, bim$V6)
  dosage_ea <- sapply(seq_len(30), function(j) if (bim$V5[j] == ea[j]) genotypes[, j] else 2 - genotypes[, j])
  R <- cor(dosage_ea)
  b <- rep(0, 30); b[10] <- 0.05
  z <- as.vector(R %*% b) * sqrt(40000) * sd(dosage_ea[, 10])
  gwas <- data.table::data.table(SNP = paste0("1:", bim$V4, "_", ea, "_", oa), CHR = 1, BP = as.numeric(bim$V4),
                                 EA = ea, OA = oa, BETA = z / 100, SE = 0.01, P = 2 * pnorm(-abs(z)),
                                 N_EFF = 40000, RSID = paste0("rs", 1:30))
  clumped <- file.path(dir, "clumped.tsv")
  data.table::fwrite(data.frame(CHR = 1, SNP = "rs10", BP = gwas$BP[10], P = gwas$P[10]), clumped, sep = "\t")
  out <- file.path(dir, "finemap")
  suppressWarnings(suppressMessages(finemap_gwas(gwas, clumped, "EUR", default_n = NA, output_finemap_dir = out,
                                                 window_kb = 100, max_causal = 1, ld_bfile = prefix)))
  regions <- data.table::fread(file.path(out, "finemap_regions.tsv"))
  expect_equal(regions$STATUS, "fine-mapped")
  expect_equal(regions$L, 1L)
  expect_match(regions$LD_SIZE_WARNING, "100 times that of the LD reference \\(400 individuals\\)")
  expect_equal(regions$N_GWAS, 40000)
  expect_equal(regions$LD_REFERENCE, "user:in_sample")
  expect_equal(regions$N_LD_OUTLIERS, 0L)
  expect_equal(regions$N_SNPS, 30L)
  expect_false(regions$MHC)
  per_variant <- data.table::fread(list.files(out, pattern = "_finemap.tsv.gz$", full.names = TRUE))
  expect_equal(per_variant$SNP[which.max(per_variant$PIP)], gwas$SNP[10])
  expect_true(all(c("LD_OUTLIER", "TILE_CORE") %in% names(per_variant)))
  expect_equal(regions$LD_FLAG, "OK")
  credible_sets <- data.table::fread(file.path(out, "finemap_credible_sets.tsv"))
  expect_equal(credible_sets$LD_FLAG, "OK")
  expect_equal(credible_sets$N_LD_OUTLIERS_IN_CS, 0L)
})

test_that("ld_mismatch_flag combines outliers, mismatch and MHC", {
  expect_equal(ld_mismatch_flag(0.01, 0L), "OK")
  expect_equal(ld_mismatch_flag(0.3, 2L, 1L, TRUE), "OUTLIER_IN_CS;LD_MISMATCH;MHC")
  expect_equal(ld_mismatch_flag(NA, 1L, 0L), "OUTLIERS_IN_REGION")
  expect_equal(ld_mismatch_flag(c(0.01, 0.5), c(0L, 0L)), c("OK", "LD_MISMATCH"))
  expect_equal(ld_mismatch_flag(numeric(), integer()), character())
  expect_equal(ld_mismatch_flag(0.01, 3L, c(0L, 2L)), c("OUTLIERS_IN_REGION", "OUTLIER_IN_CS"))
})

test_that("finemap_diagnostics summarises regions and flagged credible sets, and coloc rows get the flags", {
  dir <- withr::local_tempdir()
  data.table::fwrite(data.frame(LOCUS = c("1_100", "2_200"), STATUS = "fine-mapped", L = 1L, LD_SIZE_WARNING = c(NA, "big GWAS"),
                                LD_REFERENCE = "1000genomes", LD_MISMATCH_S = c(0.02, 0.4), N_LD_OUTLIERS = c(0L, 3L),
                                MHC = FALSE, LD_FLAG = c("OK", "LD_MISMATCH;OUTLIERS_IN_REGION")),
                     file.path(dir, "finemap_regions.tsv"), sep = "\t")
  data.table::fwrite(data.frame(LOCUS = c("1_100", "2_200"), CS = 1L, LD_FLAG = c("OK", "OUTLIER_IN_CS;LD_MISMATCH")),
                     file.path(dir, "finemap_credible_sets.tsv"), sep = "\t")
  diagnostics <- finemap_diagnostics(dir)
  expect_equal(diagnostics$summary$REGIONS, 2L)
  expect_equal(diagnostics$summary$CREDIBLE_SETS_FLAGGED, 1L)
  expect_equal(diagnostics$summary$REGIONS_LD_MISMATCH, 1L)
  expect_equal(diagnostics$summary$REGIONS_WITH_LD_OUTLIERS, 1L)
  expect_true(diagnostics$summary$LD_SIZE_WARNING)
  expect_equal(finemap_diagnostics(withr::local_tempdir())$summary$REGIONS, 0L)

  coloc <- data.frame(trait1 = "A", trait2 = "B", locus1 = c("1_100", "2_200"), locus2 = c("1_100", "9_9"),
                      PP.H4.abf = c(0.9, 0.95))
  annotated <- annotate_coloc_with_ld_flags(coloc, c(A = dir, B = withr::local_tempdir()))
  expect_equal(annotated$ld_flag1, c("OK", "LD_MISMATCH;OUTLIERS_IN_REGION"))
  expect_true(all(is.na(annotated$ld_flag2)))
  expect_equal(annotated$ld_check, c(FALSE, TRUE))
})

test_that("the fine-mapping report renders with the diagnostics", {
  dir <- withr::local_tempdir()
  data.table::fwrite(data.frame(LOCUS = "1_100", STATUS = "fine-mapped", N_SNPS = 30L, N_CS = 1L, L = 1L,
                                LD_MISMATCH_S = 0.4, N_LD_OUTLIERS = 1L, MHC = FALSE, LD_FLAG = "LD_MISMATCH"),
                     file.path(dir, "finemap_regions.tsv"), sep = "\t")
  data.table::fwrite(data.frame(LOCUS = "1_100", CS = 1L, CS_SIZE = 2L, LEAD_RSID = "rs1", LEAD_PIP = 0.9,
                                LD_FLAG = "OUTLIER_IN_CS"), file.path(dir, "finemap_credible_sets.tsv"), sep = "\t")
  html <- tempfile(fileext = ".html")
  rmd <- if (file.exists("../../R/markdown/finemap.Rmd")) "../../R/markdown/finemap.Rmd" else
    system.file("markdown", "finemap.Rmd", package = "GeneHackman")
  skip_if(!nzchar(rmd) || !file.exists(rmd), "finemap.Rmd not found")
  suppressWarnings(create_html_from_rmd(normalizePath(rmd), list(finemap_dirs = dir, labels = "Trait_A",
                                                                  multi_ancestry_dir = "none"), html))
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(page, "OUTLIER_IN_CS")
  expect_match(page, "Trait_A")
})
