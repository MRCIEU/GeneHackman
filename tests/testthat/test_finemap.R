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

  res <- finemap_gwas(
    gwas_file,
    clump_file,
    ancestry = "EUR",
    default_n = 10000L,
    output_finemap_dir = out_dir
  )

  expect_equal(nrow(res), 0)
  expect_equal(length(list.files(out_dir, pattern = "_finemap\\.tsv\\.gz$")), 0)
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
