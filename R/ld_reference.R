#' LD reference options: precomputed Pan-UKB LD pack, or 1000 Genomes genotypes via PLINK.
#' @keywords internal
ld_reference_options <- list(pan_ukb = "pan_ukb", thousand_genomes = "1000genomes")

#' Populations available in the Pan-UKB LD pack, and GeneHackman ancestry codes mapped onto them.
#' @keywords internal
ld_pack_ancestries <- c("EUR", "AFR", "AMR", "CSA", "EAS", "MID")
ld_pack_ancestry_aliases <- c(SAS = "CSA")
thousand_genomes_ancestries <- c("EUR", "EAS", "AFR", "AMR", "SAS")

#' Directory of the Pan-UKB LD pack: \code{PIPELINE_DATA_DIR/genomic_data/pan_ukb_ld}, unless
#' \code{options(GeneHackman.ld_pack_dir = ...)} points elsewhere.
#' @keywords internal
default_ld_pack_dir <- function() {
  getOption("GeneHackman.ld_pack_dir", ld_pack_dir)
}

.ld_pack_cache <- new.env(parent = emptyenv())
.ld_pack_chunk_cache <- new.env(parent = emptyenv())
.ld_pack_chunk_cache$values <- list()
.ld_pack_chunk_cache$order <- character(0)

#' Check an LD reference option, returning it unchanged.
#' @keywords internal
validate_ld_reference <- function(ld_reference) {
  if (is.null(ld_reference) || !ld_reference %in% unlist(ld_reference_options)) {
    stop("ld_reference must be one of: ", paste(unlist(ld_reference_options), collapse = ", "))
  }
  return(ld_reference)
}

#' Map a GeneHackman ancestry code onto the population used by the LD reference.
#'
#' For the Pan-UKB pack, SAS is mapped to CSA (Pan-UKB's Central/South Asian group); MID and CSA are
#' only available there.  For 1000 Genomes the five super-populations are used unchanged.
#' @param ancestry ancestry code, e.g. EUR
#' @param ld_reference one of \code{ld_reference_options}
#' @return population code used by the LD reference
#' @keywords internal
resolve_ld_ancestry <- function(ancestry, ld_reference = ld_reference_options$pan_ukb) {
  validate_ld_reference(ld_reference)
  code <- toupper(trimws(as.character(ancestry)))
  if (ld_reference == ld_reference_options$pan_ukb) {
    if (code %in% names(ld_pack_ancestry_aliases)) code <- unname(ld_pack_ancestry_aliases[code])
    allowed <- ld_pack_ancestries
  } else {
    allowed <- thousand_genomes_ancestries
  }
  if (!code %in% allowed) {
    stop("Ancestry ", ancestry, " is not available for ld_reference '", ld_reference, "'. Use one of: ",
         paste(allowed, collapse = ", "), if (ld_reference == ld_reference_options$pan_ukb) " (SAS maps to CSA)")
  }
  return(code)
}

#' Column of the pack's variant table holding GeneHackman SNP keys for a build.
#' @keywords internal
ld_pack_key_column <- function(build) {
  key <- switch(as.character(build), GRCh37 = "SNP37", GRCh38 = "SNP38", NULL)
  if (is.null(key)) stop("LD pack lookups support GRCh37 and GRCh38 only, not ", build)
  return(key)
}

#' Read \code{<pack>/<POP>/info.tsv} (format version, scale, radius, sample size, ...).
#' @keywords internal
read_ld_pack_info <- function(ancestry, pack_dir = default_ld_pack_dir()) {
  pop <- resolve_ld_ancestry(ancestry, ld_reference_options$pan_ukb)
  cache_key <- paste("info", pack_dir, pop, sep = "|")
  if (!is.null(.ld_pack_cache[[cache_key]])) return(.ld_pack_cache[[cache_key]])

  info_file <- file.path(pack_dir, pop, "info.tsv")
  if (!file.exists(info_file)) {
    stop("LD pack not found for ancestry ", ancestry, " (", pop, "): ", info_file,
         "\nSee reference_data/pan_ukb_ld/README.md to build or download the Pan-UKB LD pack.")
  }
  lines <- readLines(info_file, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  info <- stats::setNames(lapply(parts, function(x) if (length(x) > 1) x[2] else ""),
                          vapply(parts, `[`, character(1), 1))
  if (!identical(as.integer(info$format_version), 1L)) {
    stop("Unsupported LD pack format_version in ", info_file, ": ", info$format_version)
  }
  .ld_pack_cache[[cache_key]] <- info
  return(info)
}

#' Chromosomes present in the pack for one ancestry.
#' @keywords internal
ld_pack_chromosomes <- function(ancestry, pack_dir = default_ld_pack_dir()) {
  pop <- resolve_ld_ancestry(ancestry, ld_reference_options$pan_ukb)
  read_ld_pack_info(pop, pack_dir)
  files <- list.files(file.path(pack_dir, pop), pattern = "^variants_chr[0-9XY]+\\.tsv\\.gz$")
  sub("^variants_chr([0-9XY]+)\\.tsv\\.gz$", "\\1", files)
}

#' Number of individuals the pack's LD was estimated from, for one ancestry.
#' @keywords internal
ld_pack_n_samples <- function(ancestry, pack_dir = default_ld_pack_dir()) {
  as.numeric(read_ld_pack_info(ancestry, pack_dir)$n_samples)
}

#' Chunk index of one chromosome of the pack.
#' @keywords internal
load_ld_pack_chunks <- function(chr, ancestry, pack_dir = default_ld_pack_dir()) {
  pop <- resolve_ld_ancestry(ancestry, ld_reference_options$pan_ukb)
  cache_key <- paste("chunks", pack_dir, pop, chr, sep = "|")
  if (!is.null(.ld_pack_cache[[cache_key]])) return(.ld_pack_cache[[cache_key]])
  chunk_file <- file.path(pack_dir, pop, paste0("ld_chr", chr, ".chunks.tsv"))
  if (!file.exists(chunk_file)) stop("LD pack chunk index not found: ", chunk_file)
  chunks <- data.table::fread(chunk_file, colClasses = list(numeric = c("OFFSET", "NBYTES", "NVALUES")))
  .ld_pack_cache[[cache_key]] <- chunks
  return(chunks)
}

#' Variant table of one chromosome of the pack, with each row's offset inside its chunk.
#' @keywords internal
load_ld_pack_variants <- function(chr, ancestry, pack_dir = default_ld_pack_dir()) {
  ROW_OFFSET = SEG_LEN = CHUNK = IDX = NULL # silence R CMD check notes for data.table columns
  pop <- resolve_ld_ancestry(ancestry, ld_reference_options$pan_ukb)
  cache_key <- paste("variants", pack_dir, pop, chr, sep = "|")
  if (!is.null(.ld_pack_cache[[cache_key]])) return(.ld_pack_cache[[cache_key]])

  read_ld_pack_info(pop, pack_dir)
  variants_file <- file.path(pack_dir, pop, paste0("variants_chr", chr, ".tsv.gz"))
  if (!file.exists(variants_file)) stop("LD pack variants not found: ", variants_file)
  variants <- data.table::fread(variants_file, colClasses = list(character = c("CHR", "SNP37", "SNP38", "EA", "OA")),
                                na.strings = c("", "NA"))
  chunks <- load_ld_pack_chunks(chr, pop, pack_dir)
  variants[, CHUNK := findInterval(IDX, chunks$ROW_START)]
  variants[, ROW_OFFSET := cumsum(SEG_LEN) - SEG_LEN, by = CHUNK]
  .ld_pack_cache[[cache_key]] <- variants
  return(variants)
}

#' Decompress one chunk of the pack into int8 values.
#'
#' The most recently used chunks are cached (\code{options(GeneHackman.ld_chunk_cache = 64)} chunks,
#' a few MB each), because clumping and overlapping regions read the same chunks repeatedly.
#' @keywords internal
read_ld_pack_chunk <- function(chr, ancestry, chunk_row, pack_dir = default_ld_pack_dir()) {
  pop <- resolve_ld_ancestry(ancestry, ld_reference_options$pan_ukb)
  cache_key <- paste("chunk", pack_dir, pop, chr, chunk_row$CHUNK, sep = "|")
  cached <- .ld_pack_chunk_cache$values[[cache_key]]
  if (!is.null(cached)) return(cached)

  values <- decompress_ld_pack_chunk(chr, pop, chunk_row, pack_dir)
  max_entries <- getOption("GeneHackman.ld_chunk_cache", 64)
  .ld_pack_chunk_cache$order <- c(.ld_pack_chunk_cache$order, cache_key)
  .ld_pack_chunk_cache$values[[cache_key]] <- values
  while (length(.ld_pack_chunk_cache$order) > max_entries) {
    .ld_pack_chunk_cache$values[[.ld_pack_chunk_cache$order[1]]] <- NULL
    .ld_pack_chunk_cache$order <- .ld_pack_chunk_cache$order[-1]
  }
  return(values)
}

#' Read and decompress one chunk of the pack's \code{ld_chr<C>.bin} file.
#' @keywords internal
decompress_ld_pack_chunk <- function(chr, pop, chunk_row, pack_dir) {
  bin_file <- file.path(pack_dir, pop, paste0("ld_chr", chr, ".bin"))
  con <- file(bin_file, "rb")
  on.exit(close(con))
  seek(con, where = chunk_row$OFFSET)
  payload <- readBin(con, what = "raw", n = chunk_row$NBYTES)
  values <- readBin(memDecompress(payload, type = "xz"), what = "integer", size = 1, signed = TRUE,
                    n = chunk_row$NVALUES)
  if (length(values) != chunk_row$NVALUES) {
    stop("Corrupt LD pack chunk ", chunk_row$CHUNK, " on chr", chr, " (", pop, ")")
  }
  return(values)
}

#' Symmetric LD matrix for pack rows (0-based \code{IDX}), in the order given.
#' @keywords internal
ld_pack_matrix_for_rows <- function(rows, chr, ancestry, pack_dir = default_ld_pack_dir()) {
  variants <- load_ld_pack_variants(chr, ancestry, pack_dir)
  chunks <- load_ld_pack_chunks(chr, ancestry, pack_dir)
  scale <- as.numeric(read_ld_pack_info(ancestry, pack_dir)$scale)

  n <- length(rows)
  ord <- order(rows)
  sorted_rows <- rows[ord]
  out <- diag(n)
  seg_len <- variants$SEG_LEN[sorted_rows + 1]
  row_offset <- variants$ROW_OFFSET[sorted_rows + 1]
  chunk_of_row <- findInterval(sorted_rows, chunks$ROW_START)

  for (chunk_index in unique(chunk_of_row)) {
    values <- read_ld_pack_chunk(chr, ancestry, chunks[chunk_index], pack_dir)
    for (pos in which(chunk_of_row == chunk_index)) {
      i <- sorted_rows[pos]
      seg <- seg_len[pos]
      if (seg == 0 || pos == n) next
      stop_pos <- findInterval(i + seg, sorted_rows)
      if (stop_pos <= pos) next
      later <- (pos + 1):stop_pos
      r <- values[row_offset[pos] + (sorted_rows[later] - i)] / scale
      out[pos, later] <- r
      out[later, pos] <- r
    }
  }

  inverse <- integer(n)
  inverse[ord] <- seq_len(n)
  return(out[inverse, inverse, drop = FALSE])
}

#' LD matrix from the precomputed LD pack for GeneHackman SNP keys.
#'
#' @param snps character vector of SNP keys (\code{CHR:BP_EA_OA}, canonical allele order) in \code{build}
#' @param chr chromosome
#' @param ancestry ancestry code (SAS maps to CSA)
#' @param build reference build of the SNP keys (GRCh37 or GRCh38)
#' @param pack_dir LD pack directory
#' @return list with \code{matrix} (r oriented to EA) and \code{snps} (the keys found, in input order),
#'   plus \code{beyond_band_fraction}: the share of variant pairs further apart than the stored
#'   radius, which are returned as r = 0; or NULL when fewer than two keys are in the pack
#' @keywords internal
get_ld_matrix_from_pack <- function(snps, chr, ancestry, build = reference_builds$GRCh38,
                                    pack_dir = default_ld_pack_dir()) {
  variants <- load_ld_pack_variants(chr, ancestry, pack_dir)
  key_column <- ld_pack_key_column(build)
  keys <- variants[[key_column]]
  wanted <- unique(toupper(as.character(snps)))
  idx <- match(wanted, keys)
  found <- wanted[!is.na(idx)]
  if (length(found) < 2) return(NULL)
  rows <- variants$IDX[idx[!is.na(idx)]]

  R <- ld_pack_matrix_for_rows(rows, chr, ancestry, pack_dir)
  rownames(R) <- found
  colnames(R) <- found

  radius <- as.numeric(read_ld_pack_info(ancestry, pack_dir)$radius_bp)
  bp <- sort(variants$BP37[rows + 1])
  within <- findInterval(bp + radius, bp) - seq_along(bp)
  total_pairs <- length(bp) * (length(bp) - 1) / 2
  beyond <- if (total_pairs > 0) 1 - sum(within) / total_pairs else 0

  return(list(matrix = R, snps = found, beyond_band_fraction = beyond))
}

#' Make an LD matrix positive semi-definite, changing it as little as possible.
#'
#' Quantisation, the dead-zone and the zero-filled pairs beyond the stored radius leave R with
#' negative eigenvalues, which SuSiE-RSS does not accept.  When the smallest eigenvalue is below
#' \code{tolerance}, negative eigenvalues are set to 0 and the result is rescaled to a unit
#' diagonal (the nearest positive semi-definite correlation matrix, approximately).  This is used
#' instead of adding \code{|lambda_min|} to the diagonal: on 1000 Genomes EUR chr22 LD stored like the
#' pack, that ridge was 0.1-1.0 for regions of 250-1,900 variants, which shrank every correlation
#' enough to lose or invent credible sets in simulations; eigenvalue clipping left SuSiE results
#' almost unchanged.
#' @return list with \code{matrix} and \code{min_eigenvalue} (smallest eigenvalue before any change)
#' @keywords internal
make_ld_positive_semidefinite <- function(R, tolerance = -1e-6) {
  if (nrow(R) < 2) return(list(matrix = R, min_eigenvalue = 1))
  decomposition <- eigen(R, symmetric = TRUE)
  min_eigen <- min(decomposition$values)
  if (min_eigen >= tolerance) return(list(matrix = R, min_eigenvalue = min_eigen))
  vectors <- decomposition$vectors
  clipped <- vectors %*% (pmax(decomposition$values, 0) * t(vectors))
  clipped <- stats::cov2cor((clipped + t(clipped)) / 2)
  dimnames(clipped) <- dimnames(R)
  return(list(matrix = clipped, min_eigenvalue = min_eigen))
}

#' LD matrix for a set of variants from the configured LD reference.
#'
#' @param ids variant identifiers: GeneHackman SNP keys for \code{pan_ukb} and for a user-supplied \code{bfile},
#'   RSIDs for \code{1000genomes}
#' @param chr chromosome
#' @param ancestry ancestry code
#' @param build reference build of the GWAS: the SNP key column for \code{pan_ukb}, the 1000 Genomes panel
#'   for \code{1000genomes}
#' @param ld_reference one of \code{ld_reference_options}
#' @param bfile PLINK prefix of user-supplied genotypes (e.g. in-sample), in the GWAS's build; when given it is
#'   used instead of \code{ld_reference}
#' @return list with \code{matrix} and \code{snps} (identifiers in matrix order), or NULL on failure.  For
#'   1000 Genomes the panel's alleles are returned too (\code{a1}, \code{a2}), for \code{align_ld_to_gwas()}.
#' @keywords internal
get_ld_matrix <- function(ids, chr, ancestry, build = reference_builds$GRCh38,
                          ld_reference = ld_reference_options$thousand_genomes, bfile = NULL) {
  if (has_user_ld(bfile)) {
    return(compute_ld_matrix_from_bfile(ids, chr, bfile))
  }
  validate_ld_reference(ld_reference)
  if (ld_reference == ld_reference_options$pan_ukb) {
    return(get_ld_matrix_from_pack(ids, chr, ancestry, build, default_ld_pack_dir()))
  }
  return(compute_ld_matrix(ids, chr, resolve_ld_ancestry(ancestry, ld_reference), build))
}

#' TRUE when a user-supplied PLINK prefix is given.
#' @keywords internal
has_user_ld <- function(bfile) {
  !is.null(bfile) && length(bfile) == 1 && !is.na(bfile) && nzchar(bfile)
}

#' Identifier column used to match a GWAS to the LD reference.
#' @keywords internal
ld_reference_id_column <- function(ld_reference, bfile = NULL) {
  if (has_user_ld(bfile)) return("SNP")
  if (validate_ld_reference(ld_reference) == ld_reference_options$pan_ukb) "SNP" else "RSID"
}

#' Number of individuals the LD comes from: the user's .fam, the 1000 Genomes panel's .fam or the pack's info.
#' @return a number, or NA when it can't be found
#' @keywords internal
ld_reference_n_samples <- function(ld_reference, ancestry, build = reference_builds$GRCh38, bfile = NULL) {
  count_fam <- function(prefix) {
    fam <- paste0(prefix, ".fam")
    if (!file.exists(fam)) return(NA_real_)
    as.numeric(length(readLines(fam, warn = FALSE)))
  }
  if (has_user_ld(bfile)) return(count_fam(bfile))
  if (validate_ld_reference(ld_reference) == ld_reference_options$pan_ukb) {
    return(tryCatch(ld_pack_n_samples(ancestry), error = function(e) NA_real_))
  }
  tryCatch(count_fam(file.path(thousand_genomes_dir_for(build), resolve_ld_ancestry(ancestry, ld_reference))),
           error = function(e) NA_real_)
}

#' LD matrix from user-supplied PLINK genotypes (e.g. in-sample), matched on position and alleles.
#'
#' Variants are matched to the .bim by chromosome, position and both alleles (in either order), so the
#' variant IDs in the .bim can be anything unique.  Correlations are returned for the GWAS effect allele:
#' PLINK gives r between A1 dosages, and the sign is flipped for every variant whose A1 is the GWAS's other
#' allele.  GWAS variants with no match, or with more than one, are left out.
#' @param snp_keys GeneHackman SNP keys (\verb{CHR:BP_EA_OA}) in the build of \code{bfile}
#' @param chr chromosome
#' @param bfile PLINK prefix (.bed/.bim/.fam)
#' @return list with \code{matrix} and \code{snps} (the SNP keys, in matrix order), or NULL on failure
#' @keywords internal
compute_ld_matrix_from_bfile <- function(snp_keys, chr, bfile) {
  V1 = V2 = V4 = V5 = V6 = ROW = NULL # silence R CMD check notes for data.table columns
  bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE, colClasses = list(character = c(1, 2, 5, 6)))
  bim[, ROW := .I]
  bim <- bim[sub("^chr", "", V1, ignore.case = TRUE) == as.character(chr)]
  if (nrow(bim) == 0) return(NULL)
  bp <- format(bim$V4, scientific = FALSE, trim = TRUE)
  a1 <- toupper(bim$V5)
  a2 <- toupper(bim$V6)
  key_a1_first <- paste0(chr, ":", bp, "_", a1, "_", a2)
  key_a2_first <- paste0(chr, ":", bp, "_", a2, "_", a1)

  keys <- unique(toupper(as.character(snp_keys)))
  hits_a1 <- match(keys, key_a1_first)
  hits_a2 <- match(keys, key_a2_first)
  ambiguous <- !is.na(hits_a1) & !is.na(hits_a2) |
    keys %in% key_a1_first[duplicated(key_a1_first)] | keys %in% key_a2_first[duplicated(key_a2_first)]
  row <- ifelse(is.na(hits_a1), hits_a2, hits_a1)
  keep <- !is.na(row) & !ambiguous
  if (sum(keep) < 2) return(NULL)
  keys <- keys[keep]
  row <- row[keep]
  sign <- ifelse(is.na(hits_a1[keep]), -1, 1)

  ids <- bim$V2[row]
  duplicated_ids <- bim$V2[duplicated(bim$V2)]
  if (any(ids %in% c(".", duplicated_ids))) {
    stop("Variant IDs in ", bfile, ".bim must be unique (and not '.') to extract LD; set them with e.g. ",
         "plink2 --bfile ", bfile, " --set-all-var-ids @:#:\\$r:\\$a --make-bed --out <new prefix>")
  }

  order_in_bim <- order(bim$ROW[row])
  keys <- keys[order_in_bim]
  ids <- ids[order_in_bim]
  sign <- sign[order_in_bim]

  snp_file <- tempfile(fileext = ".snps")
  out_prefix <- tempfile(pattern = "user_ld_")
  on.exit(unlink(c(snp_file, paste0(out_prefix, c(".ld", ".nosex", ".log")))), add = TRUE)
  writeLines(ids, snp_file)
  cmd <- paste("plink1.9", "--bfile", bfile, "--extract", snp_file, "--r", "square",
               "--keep-allele-order", "--out", out_prefix)
  exit_code <- run_system(cmd, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  ld_file <- paste0(out_prefix, ".ld")
  if (exit_code != 0 || !file.exists(ld_file)) {
    log_file <- paste0(out_prefix, ".log")
    warning("PLINK LD computation from ", bfile, " failed for chr ", chr, ". Command:\n", cmd,
            if (file.exists(log_file)) paste0("\nPLINK log:\n", paste(readLines(log_file), collapse = "\n")),
            call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  R <- as.matrix(data.table::fread(ld_file, header = FALSE))
  if (nrow(R) != length(keys)) return(NULL)
  R <- R * outer(sign, sign)
  R[is.nan(R)] <- 0
  dimnames(R) <- list(keys, keys)
  return(list(matrix = R, snps = keys))
}

#' Align an LD matrix to the GWAS effect alleles, using the reference panel's alleles.
#'
#' For LD matched on RSID (the 1000 Genomes panel), the panel's A1/A2 are compared with each variant's EA/OA:
#' the sign of r is flipped where A1 is the GWAS's other allele, and variants whose alleles don't match in
#' either order (e.g. a different allele at a multi-allelic site, or a strand difference) are removed.
#' @param R LD matrix, rows and columns in the order of \code{gwas}
#' @param gwas GWAS rows with EA and OA
#' @param a1,a2 the panel's alleles for the same variants
#' @return list with \code{matrix}, \code{keep} (logical, rows of \code{gwas} kept) and \code{n_mismatch}
#' @keywords internal
align_ld_to_gwas <- function(R, gwas, a1, a2) {
  ea <- toupper(gwas$EA)
  oa <- toupper(gwas$OA)
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  same <- ea == a1 & oa == a2
  swapped <- ea == a2 & oa == a1
  keep <- !is.na(same) & !is.na(swapped) & (same | swapped)
  sign <- ifelse(swapped, -1, 1)[keep]
  R <- R[keep, keep, drop = FALSE] * outer(sign, sign)
  list(matrix = R, keep = keep, n_mismatch = sum(!keep))
}

#' Check consistency between z-scores and an LD matrix (SuSiE-RSS diagnostics, Zou et al. 2022).
#'
#' \code{susieR::estimate_s_rss()} estimates s, the inconsistency between z and R (0 = consistent; values
#' above ~0.1-0.2 suggest LD mismatch), and \code{susieR::kriging_rss()} compares each z-score with its
#' expectation given the others.  Variants with \code{logLR > 2} and \code{|z| > 2} are flagged as outliers,
#' as in the susieR vignette; a flipped allele shows up as one of these.
#' @return list with \code{s} and \code{outlier} (logical, one per variant)
#' @keywords internal
ld_mismatch_qc <- function(z, R, n) {
  failed <- list(s = NA_real_, outlier = rep(FALSE, length(z)))
  if (length(z) < 3) return(failed)
  tryCatch({
    s <- susieR::estimate_s_rss(z, R, n = n)
    kriging <- suppressWarnings(suppressMessages(susieR::kriging_rss(z, R, n = n, s = s)))
    dist <- kriging$conditional_dist
    outlier <- !is.na(dist$logLR) & dist$logLR > 2 & abs(dist$z) > 2
    list(s = as.numeric(s), outlier = outlier)
  }, error = function(e) {
    message("LD mismatch check failed: ", conditionMessage(e))
    failed
  })
}

#' Number of single effects (SuSiE L) for a GWAS fine-mapped with a given LD reference.
#'
#' With \code{max_causal = "auto"}, L = 1 when the GWAS sample size is more than \code{max_n_ratio} times the
#' LD reference's: a reference that much smaller than the GWAS gives miscalibrated multi-signal fine-mapping
#' (Benner et al. 2017; Zou et al. 2022), while with L = 1 the posterior inclusion probabilities don't depend on
#' the LD matrix (Weissbrod et al. 2020).  Otherwise L = 10.  A number is used as given.
#' @return list with \code{L} and \code{rule} (a short explanation recorded in finemap_regions.tsv)
#' @keywords internal
choose_max_causal <- function(max_causal, gwas_n, n_ref, max_n_ratio = 10, default_l = 10) {
  if (!identical(tolower(as.character(max_causal)), "auto")) {
    l <- suppressWarnings(as.integer(max_causal))
    if (is.na(l) || l < 1) stop("max_causal must be auto or a positive whole number, not ", max_causal)
    return(list(L = l, rule = "set by user"))
  }
  if (is.na(gwas_n) || is.na(n_ref) || n_ref <= 0) {
    return(list(L = default_l, rule = "auto: LD reference size unknown"))
  }
  ratio <- gwas_n / n_ref
  if (ratio > max_n_ratio) {
    warning(sprintf(paste("GWAS N (%.0f) is %.0f times the LD reference's (%.0f): fine-mapping with L = 1.",
                          "Use in-sample LD (ld_bfile), or set finemap.max_causal to fit more signals."),
                    gwas_n, ratio, n_ref), call. = FALSE, immediate. = TRUE)
    return(list(L = 1L, rule = sprintf("auto: N / N_REF = %.1f > %s", ratio, max_n_ratio)))
  }
  list(L = as.integer(default_l), rule = sprintf("auto: N / N_REF = %.1f <= %s", ratio, max_n_ratio))
}
