mock_liftover_bed_output <- function(bed_file_output, unmapped = NULL) {
  fixture <- testthat::test_path("data/test_data_small_hg38.bed.gz")
  if (!file.exists(fixture)) {
    testthat::skip(paste("missing fixture:", fixture))
  }
  if (grepl("\\.gz$", fixture, ignore.case = TRUE) &&
      !grepl("\\.gz$", bed_file_output, ignore.case = TRUE)) {
    con <- gzfile(fixture, "rt")
    on.exit(close(con), add = TRUE)
    writeLines(readLines(con), bed_file_output)
  } else {
    file.copy(fixture, bed_file_output, overwrite = TRUE)
  }
  if (!is.null(unmapped) && nzchar(unmapped)) {
    file.create(unmapped)
  }
}

test_that("liftover.convert_reference_build_via_liftover returns updated data frame", {
  local_mocked_bindings(
    run_liftover = function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
      mock_liftover_bed_output(bed_file_output)
    }
  )
  result <- convert_reference_build_via_liftover(
    "data/test_data_small.tsv.gz",
    reference_builds$GRCh37,
    reference_builds$GRCh38
  )

  expect_equal(nrow(result), 75048)
  expect_equal(nrow(result[is.na(result$BP), ]), 0)
  expect_false("NEW_BP" %in% colnames(result))
  expect_true("BP37" %in% colnames(result))
})

test_that("liftover.convert_reference_build_via_liftover saved output and unmapped", {
  local_mocked_bindings(
    run_liftover = function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
      mock_liftover_bed_output(bed_file_output, unmapped = unmapped)
    }
  )
  output_file <- tempfile(fileext = ".tsv.gz")
  convert_reference_build_via_liftover(
    "data/test_data_small.tsv.gz",
    reference_builds$GRCh37,
    reference_builds$GRCh38,
    output_file
  )

  expect_true(file.exists(output_file))
  expect_true(file.exists(paste0(dirname(output_file), file_prefix(output_file), ".unmapped")))

  result <- vroom::vroom(output_file, show_col_types = FALSE)
  expect_equal(nrow(result), 75048)
  expect_equal(nrow(result[is.na(result$BP), ]), 0)
  expect_false("NEW_BP" %in% colnames(result))
  expect_true("BP37" %in% colnames(result))
})
