pack_dir <- testthat::test_path("data", "ld_pack_tiny")

# A GWAS on the fixture's EUR chr22 variants, in either build, with chosen P values
fixture_gwas <- function(p, build = "GRCh37", rsid = TRUE) {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  ok <- variants$LIFT_STATUS == "ok"
  key <- if (build == "GRCh37") variants$SNP37 else variants$SNP38
  bp <- if (build == "GRCh37") variants$BP37 else variants$BP38
  gwas <- data.table::data.table(SNP = key, CHR = 22, BP = bp, EA = variants$EA, OA = variants$OA,
                                 P = p, RSID = if (rsid) paste0("rs", seq_len(nrow(variants))) else NA_character_)
  gwas[ok | build == "GRCh37"]
}

test_that("clumping.parse_clump_arguments uses PLINK defaults and rejects unsupported flags", {
  expect_equal(parse_clump_arguments(""), list(p1 = 1e-4, p2 = 0.01, r2 = 0.5, kb = 250))
  expect_equal(parse_clump_arguments(" --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"),
               list(p1 = 5e-8, p2 = 0.01, r2 = 0.01, kb = 1000))
  expect_equal(parse_clump_arguments("--clump-field P --clump-p2 0.05")$p2, 0.05)
  expect_error(parse_clump_arguments("--clump-best"), "flag/value pairs")
  expect_error(parse_clump_arguments("--clump-range genes.txt"), "Unsupported clumping argument")
  expect_error(parse_clump_arguments("--clump-p1 abc"), "Invalid value")
})

test_that("clumping.clump_gwas_with_ld_pack follows PLINK's greedy rules", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  R <- get_ld_matrix_from_pack(variants$SNP37, "22", "EUR", "GRCh37", pack_dir)$matrix
  n <- nrow(variants)
  lead <- 3
  near <- abs(variants$BP37 - variants$BP37[lead]) <= 40000
  linked <- setdiff(which(near & R[lead, ]^2 >= 0.2), lead)
  unlinked_far <- which(!near)[1]
  expect_gt(length(linked), 0)

  p <- rep(0.5, n)
  p[lead] <- 1e-12
  p[linked] <- 1e-6
  p[unlinked_far] <- 1e-9
  gwas <- fixture_gwas(p)
  clumps <- clump_gwas_with_ld_pack(gwas, "EUR", "GRCh37", " --clump-p1 1e-5 --clump-r2 0.2 --clump-kb 40",
                                    pack_dir = pack_dir)

  expect_equal(clumps$SNP[1], paste0("rs", lead))
  expect_equal(clumps$TOTAL[1], length(linked))
  expect_setequal(strsplit(gsub("\\(1\\)", "", clumps$SP2[1]), ",")[[1]], paste0("rs", linked))
  expect_true(paste0("rs", unlinked_far) %in% clumps$SNP)
  expect_equal(sum(clumps$TOTAL) + nrow(clumps), sum(p <= 0.01))
  expect_false(any(clumps$SNP %in% paste0("rs", linked)))
})

test_that("clumping.clump_gwas_with_ld_pack gives the same clumps from GRCh38 keys", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  set.seed(4)
  p <- 10^-runif(nrow(variants), 0, 12)
  p[variants$LIFT_STATUS != "ok"] <- 0.9
  clumps37 <- clump_gwas_with_ld_pack(fixture_gwas(p, "GRCh37"), "EUR", "GRCh37", " --clump-p1 1e-6 --clump-r2 0.1",
                                      pack_dir = pack_dir)
  clumps38 <- clump_gwas_with_ld_pack(fixture_gwas(p, "GRCh38"), "EUR", "GRCh38", " --clump-p1 1e-6 --clump-r2 0.1",
                                      pack_dir = pack_dir)
  expect_gt(nrow(clumps37), 1)
  expect_equal(clumps38$SNP, clumps37$SNP)
  expect_equal(clumps38$SP2, clumps37$SP2)
  expect_equal(clumps38$BP, variants$BP38[match(clumps38$SNP, paste0("rs", seq_len(nrow(variants))))])
})

test_that("clumping.clump_gwas_with_ld_pack skips variants absent from the pack and caps the window", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  gwas <- fixture_gwas(c(1e-10, rep(0.5, nrow(variants) - 1)), rsid = FALSE)
  gwas <- rbind(gwas, data.table::data.table(SNP = "22:16050001_A_C", CHR = 22, BP = 16050001, EA = "A",
                                             OA = "C", P = 1e-20, RSID = NA_character_))
  expect_message(clumps <- clump_gwas_with_ld_pack(gwas, "EUR", "GRCh37", " --clump-p1 1e-8 --clump-kb 1000",
                                                    pack_dir = pack_dir),
                 "wider than the LD pack radius")
  expect_equal(clumps$SNP, variants$SNP37[1])
  expect_equal(clumps$SP2, "NONE")

  out <- tempfile(fileext = ".clumped")
  clump_gwas_with_ld_pack(gwas, "EUR", "GRCh37", " --clump-p1 1e-8", output_file = out, pack_dir = pack_dir)
  leads <- read_clumped_leads(out)
  expect_equal(leads$SNP, variants$SNP37[1])
  expect_equal(leads$BP, variants$BP37[1])
  expect_equal(leads$P, 1e-10)

  empty <- tempfile(fileext = ".clumped")
  none <- clump_gwas_with_ld_pack(gwas[P > 1e-5], "EUR", "GRCh37", " --clump-p1 1e-8", output_file = empty,
                                  pack_dir = pack_dir)
  expect_equal(nrow(none), 0)
  expect_equal(nrow(read_clumped_leads(empty)), 0)
})

test_that("clumping.clump_gwas_with_ld_pack output works with RSID-based downstream filters", {
  variants <- load_ld_pack_variants("22", "EUR", pack_dir)
  gwas <- fixture_gwas(c(1e-10, rep(0.5, nrow(variants) - 1)))
  out <- tempfile(fileext = ".clumped")
  clump_gwas_with_ld_pack(gwas, "EUR", "GRCh37", " --clump-p1 1e-8", output_file = out, pack_dir = pack_dir)
  filtered <- filter_gwas_by_clumped_results(gwas, out)
  expect_equal(filtered$RSID, "rs1")
})
