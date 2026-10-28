test_that("coloc.coloc_analysis returns a successfull coloc analysis", {
  gwas <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)

  #suppressing warning around coloc not using great data from the GWAS
  result <- suppressWarnings(coloc_analysis(gwas,
                           gwas,
                           exposure_name = "exposure",
                           chr = 1,
                           bp = 1232132,
                           range = 5000000
  ))
  expect_equal(nrow(result), 1)
  expect_true(all(result$h0 > 0.95 ))
  expect_true(all(result$h1 < 0.01 ))
})

test_that("coloc.run_coloc_on_list_of_datasets returns a file out coloc results", {
  gwas <- "data/test_data_small.tsv.gz"

  #suppressing warning around coloc not using great data from the GWAS
  result <- suppressWarnings(run_coloc_on_list_of_datasets(c(gwas),
                                            c(gwas),
                                            exposure_name = c("exposure"),
                                            chr = c(1),
                                            bp = c(1232132),
                                            range = c(5000000)
  ))
  expect_equal(nrow(result), 1)
  expect_true(all(result$h0 > 0.95 ))
  expect_true(all(result$h1 < 0.01 ))
})

write_qtl_finemap_dir <- function(lbf = c(4, 0.2, 0.1, NA)) {
  finemap_dir <- tempfile("qtl_finemap_")
  dir.create(finemap_dir)
  locus <- tibble::tibble(
    SNP = c("1:1000_A_G", "1:2000_C_T", "1:3000_A_C", "1:4000_G_T"),
    CHR = 1,
    BP = c(1000, 2000, 3000, 4000),
    CS = c(1L, NA, NA, NA),
    LBF_1 = lbf
  )
  vroom::vroom_write(locus, file.path(finemap_dir, "1_1000_finemap.tsv.gz"))
  finemap_dir
}

qtl_test_gwas <- function(with_n = TRUE) {
  qtl <- tibble::tibble(
    SNP = c("1:1000_A_G", "1:2000_C_T", "1:3000_A_C", "1:4000_G_T"),
    CHR = "1",
    BP = c("1000", "2000", "3000", "4000"),
    BETA = c("0.12", "0.01", "0.02", "0.01"),
    SE = "0.04",
    EAF = "0.3"
  )
  if (with_n) qtl$N <- "1000"
  qtl
}

test_that("coloc.coloc_bf_bf_qtl_analysis skips the QTL when its sample size is unknown", {
  finemap_dir <- write_qtl_finemap_dir()
  expect_message(
    result <- coloc_bf_bf_qtl_analysis(finemap_dir, qtl_test_gwas(with_n = FALSE), "GENE1", 1, 1000, default_n = 0),
    "Cannot determine N for QTL data"
  )
  expect_null(result)
})

test_that("coloc.coloc_bf_bf_qtl_analysis uses the supplied priors and only fine-mapped SNPs", {
  finemap_dir <- write_qtl_finemap_dir()
  default_prior <- coloc_bf_bf_qtl_analysis(finemap_dir, qtl_test_gwas(), "GENE1", 1, 1000)
  generous_prior <- coloc_bf_bf_qtl_analysis(finemap_dir, qtl_test_gwas(), "GENE1", 1, 1000, p12 = 1e-4)

  expect_equal(default_prior$n_snps, 3L)
  expect_true(generous_prior$PP.H4.abf > default_prior$PP.H4.abf)
})

test_that("coloc.find_finemap_locus_for_region prefers the fine-mapping region that contains the position", {
  finemap_dir <- tempfile("regions_")
  dir.create(finemap_dir)
  vroom::vroom_write(tibble::tibble(SNP = "a", REGION = "A"), file.path(finemap_dir, "1_100000_finemap.tsv.gz"))
  vroom::vroom_write(tibble::tibble(SNP = "b", REGION = "B"), file.path(finemap_dir, "1_2000000_finemap.tsv.gz"))
  data.table::fwrite(data.table::data.table(LOCUS = c("1_100000", "1_2000000"), CHR = c(1, 1),
                                            START = c(0, 1600001), END = c(1600000, 2500000)),
                     file.path(finemap_dir, "finemap_regions.tsv"), sep = "\t")

  # the nearest lead to 1.5 Mb is 1_2000000, but 1.5 Mb lies inside the 1_100000 region
  locus <- find_finemap_locus_for_region(finemap_dir, 1, 1500000)
  expect_equal(locus$REGION, "A")
})
