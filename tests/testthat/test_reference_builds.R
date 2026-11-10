test_that("constants.thousand_genomes_dir_for picks the panel of the build and explains missing ones", {
  expect_match(thousand_genomes_dir_for("GRCh38"), "1000genomes/b38_dbsnp156$")
  expect_match(thousand_genomes_dir_for("GRCh37"), "1000genomes/b37_dbsnp156$")
  expect_error(thousand_genomes_dir_for("GRCh36"), "No 1000 Genomes panel for reference build GRCh36")
  skip_if(dir.exists(thousand_genomes_dir_for("GRCh37")), "the GRCh37 panel is installed")
  expect_error(thousand_genomes_dir_for("GRCh37", check = TRUE), "--grch37")
})

test_that("constants.qtl_dataset_dir uses the build's QTL files and stops when they are missing", {
  # a dataset that is never installed, so the test doesn't depend on which QTL data QTL_DATA_DIR holds
  expect_error(qtl_dataset_dir("not_a_dataset", "top_hits", "GRCh38"), "not_a_dataset/GRCh38/top_hits")
})

test_that("locus_zoom.ensembl_db_for_build picks the gene annotation of the build", {
  expect_equal(ensembl_db_for_build("GRCh38"), "EnsDb.Hsapiens.v86")
  expect_equal(ensembl_db_for_build("GRCh37"), "EnsDb.Hsapiens.v75")
  expect_error(ensembl_db_for_build("GRCh36"), "No gene annotation")
})

test_that("clumping.use_gwas_positions_in_clumped replaces panel positions with the GWAS's", {
  clumped <- tempfile(fileext = ".clumped")
  writeLines(c(" CHR    F         SNP         BP          P    TOTAL   NSIG    S05    S01   S001  S0001    SP2",
               "   1    1       rs100      10000   1.00e-10        1      0      0      0      0      1 rs200(1)",
               "   1    1       rs300      90000   2.00e-09        0      0      0      0      0      0 NONE",
               ""), clumped)
  gwas <- data.frame(RSID = c("rs100", "rs200", "rs300"), BP = c(10500, 20500, 90000), CHR = 1)

  expect_equal(use_gwas_positions_in_clumped(clumped, gwas), 1)
  result <- data.table::fread(clumped)
  expect_equal(result$BP, c(10500, 90000))
  expect_equal(result$SNP, c("rs100", "rs300"))
  expect_equal(result$SP2, c("rs200(1)", "NONE"))
})

test_that("clumping.use_gwas_positions_in_clumped drops the NA row PLINK's trailing blank lines produce", {
  clumped <- tempfile(fileext = ".clumped")
  writeLines(c(" CHR    F         SNP         BP          P    TOTAL   NSIG    S05    S01   S001  S0001    SP2",
               "   1    1       rs100      10000   1.00e-10        1      0      0      0      0      1 rs200(1)",
               "", "", "          "), clumped)
  use_gwas_positions_in_clumped(clumped, data.frame(RSID = "rs100", BP = 10500, CHR = 1))
  result <- data.table::fread(clumped)
  expect_equal(result$SNP, "rs100")
  expect_false(any(grepl("^\\s*NA\\s*$", readLines(clumped))))
  expect_equal(nrow(read_clumped_leads(clumped)), 1)
})

test_that("clumping.use_gwas_positions_in_clumped leaves an empty clump file alone", {
  clumped <- tempfile(fileext = ".clumped")
  writeLines(paste(clumped_file_columns, collapse = " "), clumped)
  expect_equal(use_gwas_positions_in_clumped(clumped, data.frame(RSID = "rs1", BP = 1)), 0)
})
