mr_results <- tempfile(fileext = ".tsv")

test_that("mr.perform_mr_on_qtl_dataset runs cis only mr against a gwas, one row per exposure", {
  listed <- NULL
  local_mocked_bindings(
    mr_list_files = function(dir, pattern, full.names) {
      listed <<- c(dir = dir, pattern = pattern)
      return("data/some_eqtlgen_hits.tsv.gz")
    },
    qtl_dataset_dir = function(dataset, kind, build) file.path("qtl", dataset, build, kind)
  )
  perform_mr_on_qtl_dataset("data/test_data_small.tsv.gz", "eqtlgen", results_output = mr_results)

  expect_equal(listed[["dir"]], file.path("qtl", "eqtlgen", "GRCh38", "top_hits"))
  expect_true(grepl(listed[["pattern"]], "cis.tsv.gz") && !grepl(listed[["pattern"]], "trans.tsv.gz"))
  expect_true(file.exists(mr_results))
  results <- vroom::vroom(mr_results, show_col_types=F)
  expect_gt(nrow(results), 0)
  expect_equal(anyDuplicated(results$EXPOSURE), 0)
  hits <- vroom::vroom("data/some_eqtlgen_hits.tsv.gz", show_col_types = F)
  # the reported SNP is one of the exposure's own instruments
  expect_true(all(paste(results$EXPOSURE, results$SNP) %in% paste(hits$EXPOSURE, toupper(hits$SNP))))
})

test_that("mr.perform_mr_on_qtl_dataset rejects datasets GeneHackman does not provide", {
  expect_error(perform_mr_on_qtl_dataset("data/test_data_small.tsv.gz", "metabrain", results_output = tempfile()),
               "Unknown QTL dataset 'metabrain'")
})

test_that("coloc.load_qtl_gwas_for_mr_result reads gwas/cis/<EXPOSURE>.tsv.gz of the dataset and build", {
  qtl_root <- withr::local_tempdir()
  local_mocked_bindings(qtl_dataset_dir = function(dataset, kind, build) file.path(qtl_root, dataset, build, kind))
  dir.create(file.path(qtl_root, "brainmeta", "GRCh37", "gwas", "cis"), recursive = TRUE)
  region <- data.frame(SNP = "1:1000_A_G", CHR = 1, BP = 1000, EA = "A", OA = "G", BETA = 0.1, SE = 0.01, P = 1e-20)
  vroom::vroom_write(region, file.path(qtl_root, "brainmeta", "GRCh37", "gwas", "cis", "GENE1.tsv.gz"))

  mr_result <- c(outcome = "cis", EXPOSURE = "GENE1")
  loaded <- load_qtl_gwas_for_mr_result(mr_result, "brainmeta", build = "GRCh37")
  expect_equal(loaded$SNP, "1:1000_A_G")
  expect_null(suppressMessages(load_qtl_gwas_for_mr_result(c(outcome = "cis", EXPOSURE = "GENE2"), "brainmeta", "GRCh37")))
})

test_that("coloc.run_coloc_on_qtl_mr_results filters MR results on the EXPOSURE column", {
  mr <- tempfile(fileext = ".tsv")
  vroom::vroom_write(data.frame(EXPOSURE = c("GENE1", "GENE2"), outcome = "cis", CHR = 1, BP = c(1000, 2000),
                                p.adjusted = 0.01), mr)
  asked <- character()
  local_mocked_bindings(load_qtl_gwas_for_mr_result = function(mr_result, qtl_dataset, build) {
    asked <<- c(asked, mr_result[["EXPOSURE"]])
    NULL
  })
  out <- tempfile(fileext = ".tsv")
  run_coloc_on_qtl_mr_results(mr, finemap_dir = tempdir(), qtl_dataset = "eqtlgen", exposures = "GENE2", output_file = out)
  expect_equal(asked, "GENE2")
})

test_that("mr.compare_interesting_mr_results creates a grouped forest plot", {
  forest_plot_output_file <- tempfile(fileext = ".png")
  compare_interesting_mr_results(mr_results, forest_plot_output_file)

  expect_true(file.exists(forest_plot_output_file))
  expect_gt(file.info(forest_plot_output_file)$size, 10000)
})
