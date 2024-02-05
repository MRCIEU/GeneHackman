test_that("gwas_comparison.compare_two_gwases_from_clumped_hits works well", {
  result <- compare_two_gwases_from_clumped_hits("data/test_data_small.tsv.gz",
                                                 "data/test_data_small.tsv.gz",
                                                 "data/clumped_snps.tsv.gz"
  )

  expect_equal(nrow(result$res), 4)
  expect_true(all(result$res$VALUE > 28))
  expect_equal(nrow(result$variants), 29)
})

test_that("gwas_comparison.run_cross_section_of_gwas_comparisons works well", {
  input_gwases <- c("data/test_data_small.tsv.gz", "data/test_data_small.tsv.gz", "data/test_data_small.tsv.gz")
  input_clumps <- c("data/clumped_snps.tsv.gz", "data/clumped_snps.tsv.gz", "data/clumped_snps.tsv.gz")

  output_file <- tempfile(fileext = ".tsv")
  variants_output_file <- tempfile(fileext = ".tsv")

  compare_replication_across_all_gwas_permutations(input_gwases, input_clumps, output_file, variants_output_file)

  expect_true(file.exists(output_file))
  results <- vroom::vroom(output_file, show_col_types=F)
  expect_equal(nrow(results), 4 * length(input_gwases))
  expect_true(all(results$VALUE > 28))

  expect_true(file.exists(variants_output_file))
  variants <- vroom::vroom(variants_output_file, show_col_types=F)
  expect_equal(nrow(variants), 29 * length(input_gwases))
})

test_that("gwas_comparison.compare_heterogeneity_across_ancestries works well", {
  input_gwases <- c("data/test_data_small.tsv.gz",
                    "data/test_data_small.tsv.gz",
                    "data/test_data_small.tsv.gz"
  )
  input_clumps <- c("data/clumped_snps.tsv.gz",
                    "data/clumped_snps.tsv.gz",
                    "data/clumped_snps.tsv.gz"
  )
  input_ancestries <- c("EUR", "AFR", "SAS")

  heterogeneity_file  <- tempfile(fileext = ".tsv")
  first_plot  <- tempfile(fileext = ".png")
  second_plot <- tempfile(fileext = ".png")

  suppressWarnings(
      compare_heterogeneity_across_ancestries(input_gwases, input_clumps, input_ancestries, heterogeneity_file, first_plot, second_plot)
  )

  expect_true(file.exists(heterogeneity_file))
  results <- vroom::vroom(heterogeneity_file, show_col_types=F)
  expect_true(all(results$Qpval == 1))

  expect_true(file.exists(first_plot))
  expect_true(file.exists(second_plot))
})
