test_that( "r_markdown test rmd works", {
  input_string <- "test=/tmp/mr_results.tsv,test2=/tmp/mr_forest.png"
  html_output_file  <- tempfile(fileext = ".html")

  input_params <- split_string_into_named_list(input_string)
  create_html_from_rmd("data/test.Rmd", input_params, html_output_file)

  expect_true(file.exists(html_output_file))
})
