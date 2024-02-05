test_that("util.vroom_snps gets just a handful of SNPS from a GWAS", {
  gwas <- "data/test_data_small.tsv.gz"
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  result <- vroom_snps(gwas, snps)
  expect_true(all(result$SNP %in% snps))
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  file <- "data/test_data_small.tsv.gz"
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  result <- get_file_or_dataframe(file, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  file <- "data/test_data_small.tsv.gz"
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  columns <- c("SNP")
  result <- get_file_or_dataframe(file, columns=columns, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 1)
})

test_that("util.get_file_or_dataframe returns a data frame with column subset if given a file", {
  file <- "data/test_data_small.tsv.gz"
  columns <- c("SNP")

  result <- get_file_or_dataframe(file, columns=columns)

  df <- vroom::vroom(file, show_col_types=F)
  expect_equal(nrow(result), nrow(df))
  expect_equal(ncol(result), 1)
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  df <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
  snps <- c("19:12436574_A_G", "10:100392738_C_T")

  result <- get_file_or_dataframe(df, snps=snps)

  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
})


test_that("util.get_file_or_dataframe returns a data frame with column subset", {
  df <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
  columns <- c("SNP")

  result <- get_file_or_dataframe(df, columns=columns)

  expect_equal(nrow(result), nrow(df))
  expect_equal(ncol(result), 1)
})


test_that("util.get_file_or_dataframe returns a data frame with column and row subset", {
  df <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  columns <- c("SNP")

  result <- get_file_or_dataframe(df, columns=columns, snps=snps)

  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 1)
})

