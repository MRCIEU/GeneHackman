# Tests for R/coloc_bf.R ( colocalization on finemapped loci)

min_locus_tbl <- function(snps, lbf1, lbf2 = NULL, cs = 1L) {
  n <- length(snps)
  if (is.null(lbf2)) {
    lbf2 <- rep(0.1, n)
  }
  tibble::tibble(
    SNP = snps,
    CS = rep(cs, length.out = n),
    LBF_1 = lbf1,
    LBF_2 = lbf2
  )
}

test_that("coloc_bf.build_lbf_matrix is rows=signals, cols=SNPs", {
  d <- min_locus_tbl(
    snps = c("1:10_A_G", "1:20_T_C"),
    lbf1 = c(0.5, 0.25),
    lbf2 = c(0.1, 0.2)
  )
  m <- GeneHackman:::build_lbf_matrix(d, c("LBF_1", "LBF_2"))
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(colnames(m), d$SNP)
  expect_equal(rownames(m), c("LBF_1", "LBF_2"))
  expect_equal(as.numeric(m[1, ]), c(0.5, 0.25))
})

test_that("coloc_bf.build_lbf_matrix returns NULL without LBF columns", {
  d <- tibble::tibble(SNP = "a", CS = 1L)
  expect_null(GeneHackman:::build_lbf_matrix(d, c("LBF_1")))
})

test_that("coloc_bf.build_lbf_matrix returns NULL when all CS are NA", {
  d <- tibble::tibble(
    SNP = c("s1", "s2"),
    CS = c(NA_integer_, NA_integer_),
    LBF_1 = c(1, 1)
  )
  expect_null(GeneHackman:::build_lbf_matrix(d, "LBF_1"))
})

test_that("coloc_bf.parse_pairwise_bf_bf_result reads coloc summary rows", {
  sm <- tibble::tibble(
    nsnps = 3L,
    hit1 = "LBF_1",
    hit2 = "LBF_1",
    PP.H0.abf = 0.1,
    PP.H1.abf = 0.2,
    PP.H2.abf = 0.15,
    PP.H3.abf = 0.15,
    PP.H4.abf = 0.4
  )
  out <- GeneHackman:::parse_pairwise_bf_bf_result(
    list(summary = sm),
    trait1 = "A",
    trait2 = "B",
    locus1 = "1_1",
    locus2 = "1_1",
    n_snps = 3L
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$trait1, "A")
  expect_equal(out$PP.H4.abf, 0.4)
})

test_that("coloc_bf.parse_pairwise_bf_bf_result fallback when no summary", {
  out <- GeneHackman:::parse_pairwise_bf_bf_result(
    list(),
    "A", "B", "1_1", "1_1",
    n_snps = 5L
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_snps, 5L)
  expect_true(is.na(out$PP.H4.abf))
})

test_that("coloc_bf.load_all_finemap_loci loads LBF and chr/bp from filename", {
  d1 <- min_locus_tbl(c("1:1_A_G", "1:2_T_C"), c(0.1, 0.2), c(0, 0))
  dir1 <- tempfile("fm_")
  dir.create(dir1)
  vroom::vroom_write(
    d1,
    file.path(dir1, "1_100000_finemap.tsv.gz")
  )
  fd <- c(trait1 = dir1)
  loci <- GeneHackman:::load_all_finemap_loci(fd)
  expect_length(loci, 1L)
  expect_length(loci$trait1, 1L)
  expect_equal(loci$trait1[[1]]$chr, 1)
  expect_equal(loci$trait1[[1]]$lead_bp, 100000)
  expect_true("LBF_1" %in% loci$trait1[[1]]$lbf_cols)
})

test_that("coloc_bf.load_all_finemap_loci skips files without LBF columns", {
  d <- tibble::tibble(SNP = "x", BETA = 1)
  dir1 <- tempfile("fm_")
  dir.create(dir1)
  vroom::vroom_write(d, file.path(dir1, "1_100000_finemap.tsv.gz"))
  expect_equal(length(GeneHackman:::load_all_finemap_loci(c(t = dir1))$t), 0L)
})

test_that("coloc_bf.coloc_bf_bf_for_locus_pair returns NULL for <2 shared SNPs", {
  l1 <- list(
    data = tibble::tibble(SNP = "a", LBF_1 = 1, CS = 1L),
    lbf_cols = "LBF_1",
    chr = 1L,
    lead_bp = 100L,
    locus_name = "1_100"
  )
  l2 <- list(
    data = tibble::tibble(SNP = "b", LBF_1 = 1, CS = 1L),
    lbf_cols = "LBF_1",
    chr = 1L,
    lead_bp = 100L,
    locus_name = "1_100"
  )
  expect_null(
    GeneHackman:::coloc_bf_bf_for_locus_pair(l1, l2, "A", "B", 1e-4, 1e-4, 5e-6)
  )
})

test_that("coloc_bf.coloc_bf_bf_for_locus_pair runs coloc on two aligned loci", {
  snps <- c("1:1_A_G", "1:2_T_C", "1:3_G_T")
  d <- min_locus_tbl(
    snps = snps,
    lbf1 = c(0.5, 0.3, 0.2),
    lbf2 = c(0.1, 0.4, 0.1)
  )
  l1 <- list(
    data = d,
    lbf_cols = c("LBF_1", "LBF_2"),
    chr = 1L,
    lead_bp = 100000L,
    locus_name = "1_100000"
  )
  l2 <- list(
    data = d,
    lbf_cols = c("LBF_1", "LBF_2"),
    chr = 1L,
    lead_bp = 100500L,
    locus_name = "1_100500"
  )
  res <- GeneHackman:::coloc_bf_bf_for_locus_pair(
    l1, l2, "trait1", "trait2", 1e-4, 1e-4, 5e-6
  )
  expect_true(!is.null(res))
  expect_true(all(c("trait1", "trait2", "locus1", "PP.H4.abf") %in% names(res)))
  expect_true(all(is.finite(res$PP.H0.abf + res$PP.H4.abf)))
})

write_pair_finemap_dirs <- function() {
  snps <- c("1:1_A_G", "1:2_T_C", "1:3_G_T")
  d <- min_locus_tbl(
    snps = snps,
    lbf1 = c(0.4, 0.35, 0.25),
    lbf2 = c(0.2, 0.2, 0.2)
  )
  dir1 <- tempfile("coc_a_")
  dir2 <- tempfile("coc_b_")
  dir.create(dir1)
  dir.create(dir2)
  vroom::vroom_write(d, file.path(dir1, "1_100000_finemap.tsv.gz"))
  vroom::vroom_write(d, file.path(dir2, "1_100200_finemap.tsv.gz"))
  c(a = dir1, b = dir2)
}

test_that("coloc_bf.run_bf_bf_coloc needs two traits with loci", {
  d <- tempfile("one_")
  dir.create(d)
  d2 <- min_locus_tbl(c("1:1_A_G", "1:2_T_C"), c(0.1, 0.2), c(0, 0))
  vroom::vroom_write(d2, file.path(d, "1_1_finemap.tsv.gz"))
  out <- run_bf_bf_coloc(c(only = d), overlap_kb = 1000)
  expect_equal(nrow(out), 0L)
})

test_that("coloc_bf.run_bf_bf_coloc pairs overlapping loci and writes file", {
  dirs <- write_pair_finemap_dirs()
  out_path <- tempfile(fileext = ".tsv")
  res <- run_bf_bf_coloc(
    finemap_dirs = dirs,
    overlap_kb = 1000,
    output_file = out_path
  )
  expect_true(nrow(res) >= 1L)
  expect_true(file.exists(out_path))
  disk <- vroom::vroom(out_path, show_col_types = FALSE)
  expect_equal(nrow(disk), nrow(res))
  expect_equal(disk$trait1[1], "a")
  expect_equal(disk$trait2[1], "b")
})

test_that("coloc_bf.run_bf_bf_coloc uses basename when finemap_dirs unnamed", {
  dirs <- write_pair_finemap_dirs()
  names(dirs) <- NULL
  res <- run_bf_bf_coloc(dirs, overlap_kb = 1000)
  expect_true(all(c(basename(dirs[1]), basename(dirs[2])) %in% c(res$trait1, res$trait2)))
})

test_that("coloc_bf.run_bf_bf_coloc finds no pair when loci on different chromosomes", {
  snps <- c("1:1_A_G", "1:2_T_C", "1:3_G_T")
  d1 <- min_locus_tbl(snps, c(0.4, 0.35, 0.25), c(0.2, 0.2, 0.2))
  d2 <- min_locus_tbl(snps, c(0.4, 0.35, 0.25), c(0.2, 0.2, 0.2))
  dir1 <- tempfile("coc_chr1_")
  dir2 <- tempfile("coc_chr2_")
  dir.create(dir1)
  dir.create(dir2)
  vroom::vroom_write(d1, file.path(dir1, "1_100000_finemap.tsv.gz"))
  vroom::vroom_write(d2, file.path(dir2, "2_100000_finemap.tsv.gz"))
  res <- run_bf_bf_coloc(c(x = dir1, y = dir2), overlap_kb = 1000)
  expect_equal(nrow(res), 0L)
})
