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
get_ld_matrix_from_pack <- function(snps, chr, ancestry, build = reference_builds$GRCh37,
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

#' Make an LD matrix positive semi-definite by adding a small ridge when needed.
#'
#' Quantisation and the zero-filled pairs beyond the stored radius can leave R with slightly
#' negative eigenvalues, which SuSiE-RSS does not accept.  When the smallest eigenvalue is below
#' \code{tolerance}, \code{|lambda_min| + 1e-6} is added to the diagonal.
#' @keywords internal
make_ld_positive_semidefinite <- function(R, tolerance = -1e-6) {
  if (nrow(R) < 2) return(list(matrix = R, ridge = 0))
  min_eigen <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eigen >= tolerance) return(list(matrix = R, ridge = 0))
  ridge <- abs(min_eigen) + 1e-6
  return(list(matrix = R + diag(ridge, nrow(R)), ridge = ridge))
}

#' LD matrix for a set of variants from the configured LD reference.
#'
#' @param ids variant identifiers: GeneHackman SNP keys for \code{pan_ukb}, RSIDs for \code{1000genomes}
#' @param chr chromosome
#' @param ancestry ancestry code
#' @param build reference build of the SNP keys (only used for \code{pan_ukb})
#' @param ld_reference one of \code{ld_reference_options}
#' @return list with \code{matrix} and \code{snps} (identifiers in matrix order), or NULL on failure
#' @keywords internal
get_ld_matrix <- function(ids, chr, ancestry, build = reference_builds$GRCh37,
                          ld_reference = ld_reference_options$thousand_genomes) {
  validate_ld_reference(ld_reference)
  if (ld_reference == ld_reference_options$pan_ukb) {
    return(get_ld_matrix_from_pack(ids, chr, ancestry, build, default_ld_pack_dir()))
  }
  return(compute_ld_matrix(ids, chr, resolve_ld_ancestry(ancestry, ld_reference)))
}

#' Identifier column used to match a GWAS to the LD reference.
#' @keywords internal
ld_reference_id_column <- function(ld_reference) {
  if (validate_ld_reference(ld_reference) == ld_reference_options$pan_ukb) "SNP" else "RSID"
}
