test_that("graphs.miami_plot given a specific range, creates a png file ", {
  miami_plot_file <- tempfile(fileext = ".png")

  result <- miami_plot("data/test_data_small.tsv.gz",
                       "data/test_data_small.tsv.gz",
                       miami_plot_file,
                       title = "ooyooo",
                       chr = 1,
                       bp = 5232132,
                       range = 1000000
  )
  expect_true(file.exists(miami_plot_file))
  print(miami_plot_file)
  expect_gt(file.info(miami_plot_file)$size, 8000)
})
