write_marker_to_rsid <- function(columns) {
  marker_file <- tempfile(fileext = ".tsv.gz")
  vroom::vroom_write(tibble::as_tibble(columns), marker_file)
  marker_file
}

test_that("data_conversions.populate_full_rsids uses the dbSNP resource of the requested build", {
  dbsnp_dir <- tempfile("dbsnp_")
  dir.create(file.path(dbsnp_dir, "b38_dbsnp156"), recursive = TRUE)
  captured_build <- NULL
  local_mocked_bindings(
    chrpos_to_rsid = function(dt, chr_col, pos_col, ea_col, nea_col, flip, alt_rsids, build, dbsnp_dir,
                              parallel_cores, ...) {
      captured_build <<- build
      dt$RSID <- "rs1"
      dt
    }
  )
  gwas <- tibble::tibble(SNP = "1:100_A_G", CHR = "1", BP = 100, EA = "A", OA = "G")

  result <- populate_full_rsids(gwas, build = reference_builds$GRCh38, dbsnp_dir = dbsnp_dir)

  expect_equal(captured_build, "b38_dbsnp156")
  expect_equal(result$RSID, "rs1")
})

test_that("data_conversions.populate_full_rsids stops when the build's dbSNP resource is missing", {
  dbsnp_dir <- tempfile("dbsnp_")
  dir.create(file.path(dbsnp_dir, "b37_dbsnp156"), recursive = TRUE)
  gwas <- tibble::tibble(SNP = "1:100_A_G", CHR = "1", BP = 100, EA = "A", OA = "G")

  expect_error(populate_full_rsids(gwas, build = reference_builds$GRCh38, dbsnp_dir = dbsnp_dir),
               "dbSNP reference for GRCh38 not found")
  expect_error(populate_full_rsids(gwas, build = reference_builds$GRCh36, dbsnp_dir = dbsnp_dir),
               "not available for reference build GRCh36")
})

test_that("data_conversions.populate_partial_rsids matches GRCh37 positions on HG37", {
  marker_file <- write_marker_to_rsid(list(HG37 = c("1:100_A_G", "1:200_C_T"), RSID = c("rs1", "rs2")))
  gwas <- tibble::tibble(SNP = c("1:200_C_T", "1:300_A_C"), CHR = "1", BP = c(200, 300),
                         EA = c("C", "A"), OA = c("T", "C"))

  result <- populate_partial_rsids(gwas, build = reference_builds$GRCh37, marker_to_rsid_file = marker_file)

  expect_equal(result$RSID, c("rs2", NA))
})

test_that("data_conversions.populate_partial_rsids uses the GRCh37 position kept by liftover for GRCh38 output", {
  marker_file <- write_marker_to_rsid(list(HG37 = c("1:100_A_G", "1:200_C_T"), RSID = c("rs1", "rs2")))
  # GRCh38 positions; BP37 holds the original GRCh37 position
  gwas <- tibble::tibble(SNP = c("1:1100_A_G", "1:1200_C_T"), CHR = "1", BP = c(1100, 1200), BP37 = c("100", "200"),
                         EA = c("A", "C"), OA = c("G", "T"))

  result <- populate_partial_rsids(gwas, build = reference_builds$GRCh38, marker_to_rsid_file = marker_file)

  expect_equal(result$RSID, c("rs1", "rs2"))
})

test_that("data_conversions.populate_partial_rsids stops for GRCh38 positions it cannot map", {
  marker_file <- write_marker_to_rsid(list(HG37 = "1:100_A_G", RSID = "rs1"))
  gwas <- tibble::tibble(SNP = "1:1100_A_G", CHR = "1", BP = 1100, EA = "A", OA = "G")

  expect_error(populate_partial_rsids(gwas, build = reference_builds$GRCh38, marker_to_rsid_file = marker_file),
               "set populate_rsid: true")
})

test_that("data_conversions.populate_partial_rsids uses an HG38 column when the marker map has one", {
  marker_file <- write_marker_to_rsid(list(HG37 = "1:100_A_G", HG38 = "1:1100_A_G", RSID = "rs1"))
  gwas <- tibble::tibble(SNP = "1:1100_A_G", CHR = "1", BP = 1100, EA = "A", OA = "G")

  result <- populate_partial_rsids(gwas, build = reference_builds$GRCh38, marker_to_rsid_file = marker_file)

  expect_equal(result$RSID, "rs1")
})
