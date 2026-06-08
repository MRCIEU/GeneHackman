test_that("gwas_formatting.standardise_gwas standardises a gwas", {
  test_gwas <- vroom::vroom("data/test_data_small.tsv.gz", show_col_types=F)
  output_file  <- tempfile(fileext = ".tsv.gz")

  standardise_gwas("data/test_data_small.tsv.gz", output_file)
  result <- vroom::vroom(output_file, show_col_types=F)

  expect_equal(nrow(result), nrow(test_gwas))
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas deletes clashing columns", {
  test_gwas <- vroom::vroom("data/test_data_tiny_existing_column.tsv.gz", show_col_types=F)
  output_file  <- tempfile(fileext = ".tsv.gz")
  map <- "N=NEW_N,SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID"
  bespoke_column_map <- split_string_into_named_list(map)

  standardise_gwas("data/test_data_tiny_existing_column.tsv.gz", output_file, input_columns = bespoke_column_map)
  result <- vroom::vroom(output_file, show_col_types=F)

  expect_true('N' %in% names(result))
  expect_true(all(result$N == 123))
  expect_false('NEW_N' %in% names(result))

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

test_that("standardise_gwas flip_alleles=FALSE preserves EA/OA order and effect direction", {
  tmp <- tempfile(fileext = ".tsv.gz")
  # Use G/C (not T/A): vroom reads bare "T" in a file as logical TRUE.
  g <- tibble::tibble(
    CHR = 1L,
    BP = 100000L,
    EA = "G",
    OA = "C",
    P = 0.05,
    BETA = 0.15,
    SE = 0.03,
    EAF = 0.42
  )
  vroom::vroom_write(g, tmp)

  out_no_flip <- standardise_gwas(
    tmp,
    tempfile(fileext = ".tsv.gz"),
    input_reference_build = reference_builds$GRCh37,
    output_reference_build = reference_builds$GRCh37,
    flip_alleles = FALSE
  )
  expect_equal(out_no_flip$EA[1], "G")
  expect_equal(out_no_flip$OA[1], "C")
  expect_equal(out_no_flip$BETA[1], 0.15)
  expect_equal(out_no_flip$EAF[1], 0.42)
  expect_match(out_no_flip$SNP[1], "^1:100000_G_C$")

  out_flip <- standardise_gwas(
    tmp,
    tempfile(fileext = ".tsv.gz"),
    input_reference_build = reference_builds$GRCh37,
    output_reference_build = reference_builds$GRCh37,
    flip_alleles = TRUE
  )
  expect_equal(out_flip$EA[1], "C")
  expect_equal(out_flip$OA[1], "G")
  expect_equal(out_flip$BETA[1], -0.15)
  expect_equal(out_flip$EAF[1], 1 - 0.42)
  expect_match(out_flip$SNP[1], "^1:100000_C_G$")
})

test_that("standardise_gwas errors when populate_eaf is TRUE and ancestry is missing", {
  expect_error(
    standardise_gwas(
      "data/test_data_tiny.tsv.gz",
      tempfile(fileext = ".tsv.gz"),
      populate_eaf = TRUE,
      ancestry = NULL
    ),
    "populate_eaf is TRUE but ancestry is missing"
  )
  expect_error(
    standardise_gwas(
      "data/test_data_tiny.tsv.gz",
      tempfile(fileext = ".tsv.gz"),
      populate_eaf = TRUE,
      ancestry = ""
    ),
    "populate_eaf is TRUE but ancestry is missing"
  )
})

test_that("standardise_gwas flip_alleles=FALSE with partial RSID succeeds (internal flip handles it)", {
  out <- standardise_gwas(
    "data/test_data_tiny.tsv.gz",
    tempfile(fileext = ".tsv.gz"),
    input_columns = "SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID",
    flip_alleles = FALSE,
    populate_rsid_option = populate_rsid_options$partial
  )
  expect_equal(nrow(out), 12L)
})

test_that("standardise_gwas allows flip_alleles FALSE with populate_rsid none", {
  out_none <- standardise_gwas(
    "data/test_data_tiny.tsv.gz",
    tempfile(fileext = ".tsv.gz"),
    input_columns = "SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID",
    flip_alleles = FALSE,
    populate_rsid_option = populate_rsid_options$none
  )
  expect_equal(nrow(out_none), 12L)
})

test_that("standardise_alleles always flips and unflip_alleles restores original", {
  g <- tibble::tibble(
    CHR = 1L,
    BP = 100L,
    EA = "T",
    OA = "A",
    BETA = 0.2,
    EAF = 0.1,
    P = 0.5,
    SE = 0.1
  )
  flipped <- GeneHackman:::standardise_alleles(g)
  expect_equal(flipped$BETA, -0.2)
  expect_equal(flipped$EA, "A")
  expect_equal(flipped$OA, "T")
  expect_equal(flipped$EAF, 1 - 0.1)
  expect_true(flipped$.ALLELES_FLIPPED[1])

  restored <- GeneHackman:::unflip_alleles(flipped)
  expect_equal(restored$BETA, 0.2)
  expect_equal(restored$EA, "T")
  expect_equal(restored$OA, "A")
  expect_equal(restored$EAF, 0.1)
  expect_false(".ALLELES_FLIPPED" %in% names(restored))
})

test_that("filter_incomplete_rows stops when all rows are incomplete", {
  bad <- tibble::tibble(
    CHR = 1L,
    BP = 100000L,
    EA = NA_character_,
    OA = "A",
    P = 0.05,
    BETA = 0.1,
    SE = 0.02,
    EAF = 0.3
  )
  tmp <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(bad, tmp)
  expect_error(
    standardise_gwas(
      tmp,
      tempfile(fileext = ".tsv.gz"),
      input_reference_build = reference_builds$GRCh37,
      output_reference_build = reference_builds$GRCh37
    ),
    "all rows have been filtered from GWAS"
  )
})

test_that("resolve_column_map errors on unresolvable column_map value", {
  expect_error(
    GeneHackman:::resolve_column_map(1L),
    "Error resolving column map"
  )
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
