test_that("collider_bias.correct_for_collider_bias works well", {
  collider_bias_results_file <- tempfile(fileext = ".tsv")
  slopehunter_file <- tempfile(fileext = ".tsv.gz")
  harmonised_results <- tempfile(fileext = ".tsv.gz")

  suppressWarnings(
      conduct_collider_bias_analysis("data/test_data_small.tsv.gz",
                                   "data/test_data_small.tsv.gz",
                                   "data/clumped_snps.tsv.gz",
                                   collider_bias_results_file,
                                   harmonised_results,
                                   slopehunter_file,
                                   p_value_thresholds = c(0.001)
    )
  )

  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file, show_col_types=F)
  expect_equal(nrow(result), 3)
  expect_true(all(result$BETA >= 1))

  expect_true(file.exists(harmonised_results))
  expect_true(file.exists(slopehunter_file))

  slopehunter_result <- vroom::vroom(slopehunter_file, show_col_types=F)
  expect_true(all(c("P", "BETA", "SE") %in% colnames(slopehunter_result)))
})
