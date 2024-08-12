test_that("collider_bias.correct_for_collider_bias works well", {
  collider_bias_results_file <- tempfile(fileext = ".tsv")
  adjusted_file <- tempfile(fileext = ".tsv.gz")
  harmonised_results <- tempfile(fileext = ".tsv.gz")

  suppressWarnings(
      conduct_collider_bias_analysis("data/test_data_small.tsv.gz",
                                   "data/test_data_small.tsv.gz",
                                   "data/clumped_snps.tsv.gz",
                                   "slopehunter",
                                   0.001,
                                   collider_bias_results_file,
                                   harmonised_results,
                                   adjusted_file,
                                   p_value_thresholds = c(0.001)
    )
  )

  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file, show_col_types=F)
  expect_equal(nrow(result), 3)
  expect_true(all(result$BETA >= 1))

  expect_true(file.exists(harmonised_results))
  expect_true(file.exists(adjusted_file))

  adjusted_result <- vroom::vroom(adjusted_file, show_col_types=F)
  expect_true(all(c("P", "BETA", "SE") %in% colnames(adjusted_result)))
})


test_that("collider_bias.correct_for_collider_bias adjusts other adjustment method too", {
  collider_bias_results_file <- tempfile(fileext = ".tsv")
  adjusted_file <- tempfile(fileext = ".tsv.gz")
  harmonised_results <- tempfile(fileext = ".tsv.gz")

  suppressWarnings(
    conduct_collider_bias_analysis("data/test_data_small.tsv.gz",
                                   "data/test_data_small.tsv.gz",
                                   "data/clumped_snps.tsv.gz",
                                   "cwls",
                                   0.001,
                                   collider_bias_results_file,
                                   harmonised_results,
                                   adjusted_file,
                                   p_value_thresholds = c(0.001)
    )
  )

  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file, show_col_types=F)
  expect_equal(nrow(result), 3)
  expect_true(all(result$BETA >= 1))

  expect_true(file.exists(harmonised_results))
  expect_true(file.exists(adjusted_file))

  adjusted_result <- vroom::vroom(adjusted_file, show_col_types=F)
  expect_true(all(c("P", "BETA", "SE") %in% colnames(adjusted_result)))
})
