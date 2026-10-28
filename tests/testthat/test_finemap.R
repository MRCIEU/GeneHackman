test_that("finemap.finemap_gwas stops when sample size cannot be resolved", {
  # Mock LD so we reach the N check; otherwise plink failure skips the locus
  # before default_n is validated.
  local_mocked_bindings(
    compute_ld_matrix = function(rsids, chr, ancestry) {
      n <- length(rsids)
      list(matrix = diag(n), snps = rsids)
    }
  )

  gwas <- tibble::tibble(
    SNP = c("1:100000_A_G", "1:100150_T_C"),
    CHR = c(1, 1),
    BP = c(100000, 100150),
    RSID = c("rs1", "rs2"),
    BETA = c(0.1, 0.05),
    SE = c(0.02, 0.02),
    P = c(0.01, 0.05),
    EA = c("A", "T"),
    OA = c("G", "C"),
    EAF = c(0.3, 0.4)
  )
  gwas_file <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(gwas, gwas_file)
  clump_file <- tempfile(fileext = ".clumped")
  writeLines(c("SNP\tCHR\tBP", "rs1\t1\t100000"), clump_file)
  out_dir <- tempfile("finemap_")
  expect_error(
    finemap_gwas(
      gwas_file,
      clump_file,
      ancestry = "EUR",
      default_n = NA,
      output_finemap_dir = out_dir
    ),
    regexp = "Fine-mapping requires GWAS sample size"
  )
})

test_that("finemap.finemap_gwas writes no locus files when clump has no rows", {
  gwas <- tibble::tibble(
    SNP = c("1:100000_A_G", "1:100150_T_C"),
    CHR = c(1, 1),
    BP = c(100000, 100150),
    RSID = c("rs1", "rs2"),
    BETA = c(0.1, 0.05),
    SE = c(0.02, 0.02),
    P = c(0.01, 0.05),
    EA = c("A", "T"),
    OA = c("G", "C"),
    EAF = c(0.3, 0.4)
  )
  gwas_file <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(gwas, gwas_file)

  clump_file <- tempfile(fileext = ".clumped")
  writeLines("SNP\tCHR\tBP", clump_file)

  out_dir <- tempfile("finemap_")
  dir.create(out_dir)
  complete_file <- file.path(out_dir, paste0("finemap_complete_", basename(out_dir), ".txt"))

  res <- finemap_gwas(
    gwas_file,
    clump_file,
    ancestry = "EUR",
    default_n = 10000L,
    output_finemap_dir = out_dir
  )

  expect_equal(nrow(res), 0)
  expect_equal(length(list.files(out_dir, pattern = "_finemap\\.tsv\\.gz$")), 0)
  expect_true(file.exists(complete_file))
  expect_identical(readLines(complete_file), "0")
})

test_that("finemap.compute_ld_matrix returns NULL when plink fails", {
  local_mocked_bindings(
    run_system = function(command, wait = TRUE, intern = FALSE,
                          ignore.stdout = FALSE, ignore.stderr = FALSE) {
      1L
    }
  )
  expect_null(compute_ld_matrix(c("rs1", "rs2"), 1L, "EUR"))
})

test_that("finemap.run_susie_for_locus maps per-CS LBF_k and credible sets", {
  local_mocked_bindings(
    run_susie_rss_impl = function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
      p <- length(z)
      list(
        lbf_variable = matrix(c(0.5, 0.25, 0.5, 0.75), nrow = 2, ncol = p),
        sets = list(cs = list(c(1L, 2L)))
      )
    }
  )

  snp_info <- tibble::tibble(
    SNP = c("1:100000_A_G", "1:100150_T_C"),
    CHR = c(1, 1),
    BP = c(100000, 100150),
    RSID = c("rs1", "rs2")
  )
  z <- c(1, 2)
  R <- diag(2)

  out <- run_susie_for_locus(
    z_scores = z,
    ld_matrix = R,
    snp_info = snp_info,
    n = 1000L,
    lead_snp = "rs1"
  )

  expect_equal(nrow(out), 2)
  expect_equal(out$LBF_1, c(0.5, 0.5))
  expect_equal(out$CS, c(1L, 1L))
})

test_that("finemap.finemap_gwas runs end-to-end with mocked LD and SuSiE", {
  local_mocked_bindings(
    compute_ld_matrix = function(rsids, chr, ancestry) {
      n <- length(rsids)
      list(matrix = diag(n), snps = rsids)
    },
    run_susie_rss_impl = function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
      p <- length(z)
      list(
        lbf_variable = matrix(1, nrow = 2, ncol = p),
        sets = list(cs = list(seq_len(p)))
      )
    }
  )

  gwas <- tibble::tibble(
    SNP = c("1:100000_A_G", "1:100150_T_C"),
    CHR = c(1, 1),
    BP = c(100000, 100150),
    RSID = c("rs1", "rs2"),
    BETA = c(0.1, 0.05),
    SE = c(0.02, 0.02),
    P = c(0.01, 0.05),
    EA = c("A", "T"),
    OA = c("G", "C"),
    EAF = c(0.3, 0.4)
  )
  gwas_file <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(gwas, gwas_file)

  clump_file <- tempfile(fileext = ".clumped")
  writeLines(c("SNP\tCHR\tBP", "rs1\t1\t100000"), clump_file)

  out_dir <- tempfile("finemap_")
  dir.create(out_dir)
  complete_file <- file.path(out_dir, paste0("finemap_complete_", basename(out_dir), ".txt"))

  res <- finemap_gwas(
    gwas_file,
    clump_file,
    ancestry = "EUR",
    default_n = 10000L,
    output_finemap_dir = out_dir,
    window_kb = 500
  )

  expect_equal(nrow(res), 2)
  expect_setequal(res$RSID, c("rs1", "rs2"))
  expect_true(all(res$CS == 1L))

  locus_file <- file.path(out_dir, "1_100000_finemap.tsv.gz")
  expect_true(file.exists(locus_file))
  locus_read <- vroom::vroom(locus_file, show_col_types = FALSE)
  expect_equal(nrow(locus_read), 2)
  expect_true("LBF_1" %in% names(locus_read))
  expect_setequal(locus_read$SNP, gwas$SNP)
  expect_false(any(grepl("^LEAD_", names(locus_read))))
  expect_true("PIP" %in% names(locus_read))
  expect_true(file.exists(complete_file))
  expect_identical(readLines(complete_file), "1")
})

test_that("finemap.run_susie_for_locus returns empty tibble when SuSiE errors", {
  local_mocked_bindings(
    run_susie_rss_impl = function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
      stop("mock failure")
    }
  )

  snp_info <- tibble::tibble(
    SNP = "1:100000_A_G",
    CHR = 1,
    BP = 100000,
    RSID = "rs1"
  )
  out <- run_susie_for_locus(
    z_scores = 1,
    ld_matrix = matrix(1, 1, 1),
    snp_info = snp_info,
    n = 1000L,
    lead_snp = "rs1"
  )
  expect_equal(nrow(out), 0)
})


test_that("finemap.run_susie_for_locus reports PIPs and a credible-set summary", {
  local_mocked_bindings(
    run_susie_rss_impl = function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
      list(
        pip = c(0.9, 0.08, 0.02),
        lbf = c(12.5, 0.1),
        lbf_variable = matrix(c(5, 2, 1, 0, 0, 0), nrow = 2, byrow = TRUE),
        converged = TRUE,
        sets = list(
          cs = list(L1 = c(1L, 2L)),
          cs_index = 1L,
          coverage = 0.98,
          purity = data.frame(min.abs.corr = 0.8, mean.abs.corr = 0.85, median.abs.corr = 0.86)
        )
      )
    }
  )
  snp_info <- tibble::tibble(
    SNP = c("1:100_A_G", "1:200_C_T", "1:300_A_C"),
    CHR = c(1, 1, 1),
    BP = c(100, 200, 300),
    RSID = c("rs1", "rs2", "rs3")
  )

  out <- run_susie_for_locus(z_scores = c(5, 3, 1), ld_matrix = diag(3), snp_info = snp_info,
                             n = 1000L, lead_snp = "rs1")

  expect_equal(out$PIP, c(0.9, 0.08, 0.02))
  expect_equal(out$CS, c(1L, 1L, NA))
  expect_true(attr(out, "converged"))
  cs <- attr(out, "cs_summary")
  expect_equal(nrow(cs), 1)
  expect_equal(cs$CS_SIZE, 2)
  expect_equal(cs$LEAD_RSID, "rs1")
  expect_equal(cs$LEAD_PIP, 0.9)
  expect_equal(cs$COVERAGE, 0.98)
  expect_equal(cs$MIN_ABS_CORR, 0.8)
  expect_equal(cs$LOG_BF, 12.5)
})

test_that("finemap.define_finemap_regions merges overlapping lead windows and names regions by the top lead", {
  leads <- data.table::data.table(
    SNP = c("rs1", "rs2", "rs3", "rs4"),
    CHR = c(1, 1, 1, 2),
    BP = c(1e6, 1.01e6, 5e6, 1e6),
    P = c(1e-8, 1e-12, 1e-9, 1e-10)
  )

  regions <- define_finemap_regions(leads, window_bp = 5e5, max_region_bp = 3e6)

  expect_equal(nrow(regions), 3)
  expect_equal(regions$LOCUS, c("1_1010000", "1_5000000", "2_1000000"))
  expect_equal(regions$N_LEADS, c(2L, 1L, 1L))
  expect_equal(regions$START[1], 5e5)
  expect_equal(regions$END[1], 1.51e6)
})

test_that("finemap.define_finemap_regions caps merged regions without overlapping them", {
  leads <- data.table::data.table(
    SNP = paste0("rs", 1:6),
    CHR = rep(1, 6),
    BP = c(1e6, 1.8e6, 2.6e6, 3.4e6, 4.2e6, 5.0e6)
  )

  regions <- define_finemap_regions(leads, window_bp = 5e5, max_region_bp = 2e6)

  expect_true(all(regions$END - regions$START <= 2e6))
  expect_true(all(regions$START[-1] > regions$END[-nrow(regions)]))
  expect_equal(sum(regions$N_LEADS), 6)
  for (i in seq_len(nrow(leads))) {
    expect_equal(sum(leads$BP[i] >= regions$START & leads$BP[i] <= regions$END), 1)
  }
})

test_that("finemap.finemap_gwas fine-maps nearby leads once and writes region and credible-set tables", {
  local_mocked_bindings(
    compute_ld_matrix = function(rsids, chr, ancestry) {
      n <- length(rsids)
      list(matrix = diag(n), snps = rsids)
    },
    run_susie_rss_impl = function(z, R, n, L, coverage, min_abs_corr, verbose = FALSE) {
      p <- length(z)
      list(
        pip = rep(1 / p, p),
        lbf = 3,
        lbf_variable = matrix(1, nrow = 1, ncol = p),
        converged = TRUE,
        sets = list(cs = list(L1 = seq_len(p)), cs_index = 1L, coverage = 0.95,
                    purity = data.frame(min.abs.corr = 1, mean.abs.corr = 1, median.abs.corr = 1))
      )
    }
  )

  gwas <- tibble::tibble(
    SNP = c("1:100000_A_G", "1:100150_C_T", "1:120000_A_C"),
    CHR = c(1, 1, 1),
    BP = c(100000, 100150, 120000),
    RSID = c("rs1", "rs2", "rs3"),
    BETA = c(0.1, 0.05, 0.2),
    SE = c(0.02, 0.02, 0.02),
    P = c(0.01, 0.05, 1e-20),
    EA = c("A", "C", "A"),
    OA = c("G", "T", "C"),
    EAF = c(0.3, 0.4, 0.2)
  )
  gwas_file <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(gwas, gwas_file)
  clump_file <- tempfile(fileext = ".clumped")
  writeLines(c("CHR\tF\tSNP\tBP\tP", "1\t1\trs1\t100000\t1e-9", "1\t1\trs3\t120000\t1e-20"), clump_file)
  out_dir <- tempfile("finemap_")

  finemap_gwas(gwas_file, clump_file, ancestry = "EUR", default_n = 10000L,
               output_finemap_dir = out_dir, window_kb = 500)

  expect_equal(list.files(out_dir, pattern = "_finemap\\.tsv\\.gz$"), "1_120000_finemap.tsv.gz")
  regions <- data.table::fread(file.path(out_dir, "finemap_regions.tsv"))
  expect_equal(nrow(regions), 1)
  expect_equal(regions$LEAD_SNP, "rs3")
  expect_equal(regions$N_LEADS, 2)
  expect_equal(regions$STATUS, "fine-mapped")
  expect_equal(regions$N_CS, 1)
  credible_sets <- data.table::fread(file.path(out_dir, "finemap_credible_sets.tsv"))
  expect_equal(nrow(credible_sets), 1)
  expect_equal(credible_sets$LOCUS, "1_120000")
  expect_equal(credible_sets$CS_SIZE, 3)
  expect_identical(readLines(file.path(out_dir, paste0("finemap_complete_", basename(out_dir), ".txt"))), "1")
})
