test_that("liftover.convert_reference_build_via_liftover returns updated data frame", {
  local_mock(
      system = function(command, wait, intern) {
        bed_output_file <- unlist(strsplit(command, " "))[4]
        file.copy("data/test_data_small_hg38.bed.gz", bed_output_file)
      }
  )
  result <- convert_reference_build_via_liftover("data/test_data_small.tsv.gz",
                                       reference_builds$GRCh37,
                                       reference_builds$GRCh38
  )

  expect_equal(nrow(result), 75048)
  expect_equal(nrow(result[is.na(result$BP), ]), 0)
  expect_false("NEW_BP" %in% colnames(result))
  expect_true("BP37" %in% colnames(result))
})

test_that("liftover.convert_reference_build_via_liftover saved output and unmapped", {
  local_mock(
    system = function(command, wait) {
      bed_output_file <- unlist(strsplit(command, " "))[4]
      unmapped_file <- unlist(strsplit(command, " "))[5]
      file.copy("data/test_data_small_hg38.bed.gz", bed_output_file)
      file.create(unmapped_file)
    }
  )
  output_file <- tempfile(fileext = ".tsv.gz")
  convert_reference_build_via_liftover("data/test_data_small.tsv.gz",
                                                 reference_builds$GRCh37,
                                                 reference_builds$GRCh38,
                                                 output_file
  )

  expect_true(file.exists(output_file))
  expect_true(file.exists(paste0(dirname(output_file), file_prefix(output_file), ".unmapped")))

  result <- vroom::vroom(output_file, show_col_types=F)
  expect_equal(nrow(result), 75048)
  expect_equal(nrow(result[is.na(result$BP), ]), 0)
  expect_false("NEW_BP" %in% colnames(result))
  expect_true("BP37" %in% colnames(result))
})
