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

test_that("gwas_comparison.run_cross_section_of_gwas_comparisons deson't error with only one gwas", {
  input_gwases <- "data/test_data_small.tsv.gz"
  input_clumps <- "data/clumped_snps.tsv.gz"

  output_file <- tempfile(fileext = ".tsv")
  variants_output_file <- tempfile(fileext = ".tsv")

  compare_replication_across_all_gwas_permutations(input_gwases, input_clumps, output_file, variants_output_file)

  expect_true(file.exists(output_file))
  results <- vroom::vroom(output_file, show_col_types=F)
  expect_equal(nrow(results), 0)

  expect_true(file.exists(variants_output_file))
  variants <- vroom::vroom(variants_output_file, show_col_types=F)
  expect_equal(nrow(variants), 0)
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

test_that("gwas_comparison.expected_vs_observed_replication counts observed replications with a two-sided p-value", {
  # replication z-scores of 1.8 (two-sided p = 0.072) and 2.5 (two-sided p = 0.012)
  result <- expected_vs_observed_replication(
    b_disc = rep(1, 4), se_disc = rep(0.1, 4),
    b_rep = c(0.18, 0.18, 0.25, 0.25), se_rep = rep(0.1, 4),
    alpha = 0.05
  )

  observed <- result$res$VALUE[result$res$METRIC == "P-value" & result$res$DATUM == "Observed"]
  expect_equal(observed, 2)
  expect_equal(result$variants$OBSERVED, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("gwas_comparison.compare_replication_across_all_gwas_permutations keeps every pair with four GWAS", {
  input_gwases <- vapply(c("gwas_a", "gwas_b", "gwas_c", "gwas_d"), function(name) {
    path <- file.path(tempdir(), paste0(name, ".tsv.gz"))
    file.copy("data/test_data_small.tsv.gz", path, overwrite = TRUE)
    path
  }, character(1))
  input_clumps <- rep("data/clumped_snps.tsv.gz", 4)

  output_file <- tempfile(fileext = ".tsv")
  variants_output_file <- tempfile(fileext = ".tsv")
  compare_replication_across_all_gwas_permutations(input_gwases, input_clumps, output_file, variants_output_file)

  results <- vroom::vroom(output_file, show_col_types = FALSE)
  expect_equal(nrow(results), 4 * 6)
  expect_setequal(unique(results$COMPARISON), c(
    "gwas_a_vs_gwas_b", "gwas_a_vs_gwas_c", "gwas_a_vs_gwas_d",
    "gwas_b_vs_gwas_c", "gwas_b_vs_gwas_d", "gwas_c_vs_gwas_d"
  ))
})
