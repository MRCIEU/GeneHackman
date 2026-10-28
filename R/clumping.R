#' Header of a PLINK 1.9 \code{.clumped} file, used for the output of the LD-pack clumper too.
#' @keywords internal
clumped_file_columns <- c("CHR", "F", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")

#' Parse PLINK \code{--clump-*} arguments.
#'
#' Supports \code{--clump-p1}, \code{--clump-p2}, \code{--clump-r2} and \code{--clump-kb} (with PLINK 1.9's
#' defaults 1e-4, 0.01, 0.5 and 250), plus \code{--clump-field P}.  Anything else is an error, so a
#' setting is never silently ignored.
#' @param arguments string such as \code{"--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"}
#' @return list with p1, p2, r2 and kb
#' @keywords internal
parse_clump_arguments <- function(arguments) {
  settings <- list(p1 = 1e-4, p2 = 0.01, r2 = 0.5, kb = 250)
  tokens <- strsplit(trimws(paste(arguments, collapse = " ")), "\\s+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) %% 2 != 0) stop("plink_clump_arguments must be flag/value pairs: ", arguments)
  for (i in seq_len(length(tokens) %/% 2) * 2 - 1) {
    flag <- tokens[i]
    value <- tokens[i + 1]
    if (flag == "--clump-field") {
      if (value != "P") stop("--clump-field other than P is not supported with ld_reference: pan_ukb")
      next
    }
    name <- switch(flag, "--clump-p1" = "p1", "--clump-p2" = "p2", "--clump-r2" = "r2", "--clump-kb" = "kb", NULL)
    if (is.null(name)) {
      stop("Unsupported clumping argument for ld_reference: pan_ukb: ", flag,
           " (supported: --clump-p1, --clump-p2, --clump-r2, --clump-kb)")
    }
    number <- suppressWarnings(as.numeric(value))
    if (is.na(number) || number < 0) stop("Invalid value for ", flag, ": ", value)
    settings[[name]] <- number
  }
  return(settings)
}

#' Clump a GWAS using the precomputed LD pack, following PLINK 1.9's \code{--clump}.
#'
#' Variants are first matched to the pack on the GeneHackman SNP key (\code{CHR:BP_EA_OA}) in
#' \code{build}; variants absent from the pack are dropped, as PLINK drops variants absent from the
#' reference genotypes.  Then, greedily: the most significant unclumped variant with P <= p1 becomes an
#' index variant, and every unclumped variant with P <= p2 within \code{kb} of it and r^2 >= r2 joins its
#' clump.  Each variant joins at most one clump.  The window is capped at the pack's stored radius,
#' because r is not stored for pairs further apart.
#'
#' @param gwas standardised GWAS (file or data frame) with SNP, CHR, BP, P and optionally RSID
#' @param ancestry ancestry code (SAS maps to CSA)
#' @param build reference build of the GWAS (GRCh37 or GRCh38)
#' @param clump_arguments PLINK-style clumping arguments (see \code{parse_clump_arguments})
#' @param output_file optional path for a PLINK-format \code{.clumped} file
#' @param pack_dir LD pack directory
#' @return data.table in PLINK \code{.clumped} layout; the SNP column holds the GWAS RSID when available,
#'   otherwise the SNP key
#' @export
clump_gwas_with_ld_pack <- function(gwas, ancestry, build = reference_builds$GRCh37, clump_arguments = "",
                                    output_file = NULL, pack_dir = default_ld_pack_dir()) {
  SNP = CHR = BP = P = RSID = ROW = NULL # silence R CMD check notes for data.table columns
  settings <- parse_clump_arguments(clump_arguments)
  radius_kb <- as.numeric(read_ld_pack_info(ancestry, pack_dir)$radius_bp) / 1000
  if (settings$kb > radius_kb) {
    message("Clumping window of ", settings$kb, " kb is wider than the LD pack radius (", radius_kb,
            " kb): r is not stored for variants further apart, so clumping uses ", radius_kb, " kb")
  }
  window_bp <- min(settings$kb, radius_kb) * 1000

  gwas <- data.table::as.data.table(get_file_or_dataframe(gwas))
  required <- c("SNP", "CHR", "BP", "P")
  if (!all(required %in% names(gwas))) stop("GWAS must contain columns: ", paste(required, collapse = ", "))
  if (!"RSID" %in% names(gwas)) gwas[, RSID := NA_character_]
  gwas <- gwas[!is.na(P) & P <= settings$p2, list(SNP = toupper(SNP), CHR = as.character(CHR), BP, P, RSID)]

  key_column <- ld_pack_key_column(build)
  clumps <- list()
  n_absent <- 0L
  for (chr in unique(gwas[P <= settings$p1, CHR])) {
    variants <- load_ld_pack_variants(chr, ancestry, pack_dir)
    pool <- gwas[CHR == chr]
    in_pack <- pool$SNP %in% variants[[key_column]]
    n_absent <- n_absent + sum(!in_pack)
    pool <- pool[in_pack]
    pool[, ROW := variants$IDX[match(SNP, variants[[key_column]])]]
    data.table::setorder(pool, P, BP)
    clumped <- rep(FALSE, nrow(pool))

    for (i in which(pool$P <= settings$p1)) {
      if (clumped[i]) next
      clumped[i] <- TRUE
      nearby <- which(!clumped & abs(pool$BP - pool$BP[i]) <= window_bp)
      members <- integer(0)
      if (length(nearby) > 0) {
        R <- ld_pack_matrix_for_rows(c(pool$ROW[i], pool$ROW[nearby]), chr, ancestry, pack_dir)
        r2 <- R[1, -1]^2
        members <- nearby[r2 >= settings$r2]
        clumped[members] <- TRUE
      }
      clumps[[length(clumps) + 1]] <- format_clump_row(pool, i, members)
    }
  }
  if (n_absent > 0) {
    message(n_absent, " GWAS variants with P <= ", settings$p2, " are not in the LD pack and were skipped")
  }

  result <- if (length(clumps) > 0) data.table::rbindlist(clumps) else
    data.table::as.data.table(stats::setNames(replicate(length(clumped_file_columns), character(0), simplify = FALSE),
                                              clumped_file_columns))
  if (nrow(result) > 0) {
    result <- result[order(as.numeric(CHR), as.numeric(P))]
  }
  message(nrow(result), " clumps from ", sum(gwas$P <= settings$p1), " variants with P <= ", settings$p1)
  if (!is.null(output_file)) write_clumped_file(result, output_file)
  return(result)
}

#' One row of a PLINK-format clump table.
#' @keywords internal
format_clump_row <- function(pool, index, members) {
  label <- function(rows) data.table::fifelse(is.na(pool$RSID[rows]) | !nzchar(pool$RSID[rows]),
                                              pool$SNP[rows], pool$RSID[rows])
  p <- pool$P[members]
  data.table::data.table(
    CHR = pool$CHR[index], F = 1L, SNP = label(index), BP = pool$BP[index], P = pool$P[index],
    TOTAL = length(members), NSIG = sum(p > 0.05), S05 = sum(p > 0.01 & p <= 0.05), S01 = sum(p > 0.001 & p <= 0.01),
    S001 = sum(p > 1e-4 & p <= 0.001), S0001 = sum(p <= 1e-4),
    SP2 = if (length(members) == 0) "NONE" else paste0(label(members), "(1)", collapse = ",")
  )
}

#' Write a clump table in PLINK's whitespace-separated \code{.clumped} layout.
#' @keywords internal
write_clumped_file <- function(clumps, output_file) {
  P = NULL # silence R CMD check note for data.table column
  if (nrow(clumps) == 0) {
    writeLines(paste(clumped_file_columns, collapse = " "), output_file)
    return(invisible(output_file))
  }
  out <- data.table::copy(clumps)
  out[, P := format(as.numeric(P), digits = 4, scientific = TRUE)]
  data.table::fwrite(out, output_file, sep = " ", quote = FALSE)
  return(invisible(output_file))
}
