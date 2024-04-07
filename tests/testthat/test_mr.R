mr_results <- tempfile(fileext = ".tsv")

test_that("mr.perform_mr_on_pqtl_datasets runs cis only mr against a gwas", {
  local_mock(
    list.files = function(dir, pattern, full.names) {
      return("data/some_eqtlgen_hits.tsv.gz")
    }
  )
  perform_mr_on_eqtlgen_datasets("data/test_data_small.tsv.gz", results_output = mr_results)

  expect_true(file.exists(mr_results))
  results <- vroom::vroom(mr_results, show_col_types=F)
  expect_equal(nrow(results), 4)
})

test_that("mr.compare_interesting_mr_results creates a grouped forest plot", {
  forest_plot_output_file <- tempfile(fileext = ".png")
  compare_interesting_mr_results(mr_results, forest_plot_output_file)

  expect_true(file.exists(forest_plot_output_file))
  expect_gt(file.info(forest_plot_output_file)$size, 10000)
})
