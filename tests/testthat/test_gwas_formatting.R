test_that("gwas_formatting.standardise_gwas standardises a gwas", {
  test_gwas <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
  output_file  <- tempfile(fileext = ".tsv.gz")

  standardise_gwas("data/test_data_small.tsv.gz", output_file)
  result <- vroom::vroom(output_file, show_col_types=F)

  expect_equal(nrow(result), nrow(test_gwas))
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas with bespoke_column_map standardises a gwas", {
  output_file  <- tempfile(fileext = ".tsv.gz")
  map <- "SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID"
  bespoke_column_map <- split_string_into_named_list(map)

  standardise_gwas("data/test_data_tiny.tsv.gz", output_file, input_columns = bespoke_column_map)
  result <- vroom::vroom(output_file, show_col_types=F)

  expect_equal(nrow(result), 12)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas with input_column_map as a string standardises a gwas", {
  output_file  <- tempfile(fileext = ".tsv.gz")
  map <- "SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID"

  standardise_gwas("data/test_data_tiny.tsv.gz", output_file, input_columns = map)
  result <- vroom::vroom(output_file, show_col_types=F)

  expect_equal(nrow(result), 12)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas standardises a METAL output", {
  output_file  <- tempfile(fileext = ".tsv.gz")

  standardise_gwas("data/metal_output.txt", output_file, input_columns = "metal")
  result <- vroom::vroom(output_file, show_col_types=F)

  metal_lines <-  as.numeric(R.utils::countLines("data/metal_output.txt")) - 1
  expect_equal(nrow(result), metal_lines)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas standardises an OpenGWAS PheWAS output", {
  output_file  <- tempfile(fileext = ".tsv.gz")

  standardise_gwas("data/opengwas_phewas_result.tsv", output_file, input_columns = "opengwas_phewas")
  result <- vroom::vroom(output_file, show_col_types=F)

  opengwas_lines <-  as.numeric(R.utils::countLines("data/opengwas_phewas_result.tsv")) - 1
  expect_equal(nrow(result), opengwas_lines)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas standardises an gwama output", {
  output_file  <- tempfile(fileext = ".tsv.gz")

  filename <- "data/gwama_output.txt.gz"
  standardise_gwas(filename, output_file, input_columns = "gwama")
  result <- vroom::vroom(output_file, show_col_types=F)

  number_of_lines <-  as.numeric(R.utils::countLines(filename)) - 1
  expect_equal(nrow(result), number_of_lines)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas with remove_extra_columns removes extra columns", {
  filename <- "data/gwama_output.txt.gz"
  result <- standardise_gwas(filename, output_file=F, input_columns = "gwama", remove_extra_columns = T)

  suppressWarnings(
    expect_true(all(is.null(result$effects), is.null(result$q_statistic), is.null(result$i2), is.null(result$n_studies)))
  )
  expect_true(all(!is.null(result$BETA), !is.null(result$SE), !is.null(result$P), !is.null(result$RSID)))
})


test_that("gwas_formatting.standardise_gwas standardises an IEU ukb pipeline output", {
  output_file  <- tempfile(fileext = ".tsv.gz")

  filename <- "data/ukb_pipeline.txt.gz"
  standardise_gwas(filename, output_file, input_columns = "ieu_ukb")
  result <- vroom::vroom(output_file, show_col_types=F)

  opengwas_lines <-  as.numeric(R.utils::countLines(filename)) - 1
  expect_equal(nrow(result), opengwas_lines)
  floating_point_tolerance <- 1e-10
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.convert_beta_to_or and back returns the same results", {
  original_gwas <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)

  gwas <- convert_beta_to_or(original_gwas)
  gwas <- convert_or_to_beta(gwas)

  floating_point_tolerance <- 1e-10
  expect_true(all(abs(gwas$BETA - original_gwas$BETA) < floating_point_tolerance))
  expect_true(all(abs(gwas$SE - original_gwas$SE) < floating_point_tolerance))
})

# test_that("gwas_formatting.convert_beta_and_se_to_z_score and back returns the same results", {
#   original_gwas <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
#
#   gwas <- convert_beta_and_se_to_z_score(original_gwas)
#   gwas <- convert_z_score_to_beta(gwas)
#   print(abs(gwas$BETA - original_gwas$BETA))
#
#   floating_point_tolerance <- 1e-10
#   expect_true(all(abs(gwas$BETA - original_gwas$BETA) < floating_point_tolerance))
#   expect_true(all(abs(gwas$SE - original_gwas$SE) < floating_point_tolerance))
# })
