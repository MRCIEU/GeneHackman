reference_builds <- list(GRCh36="GRCh36", GRCh37="GRCh37", GRCh38="GRCh38")
available_liftover_conversions <- list()
available_liftover_conversions[[paste0(reference_builds$GRCh36, reference_builds$GRCh37)]] <- "hg18ToHg19.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh36, reference_builds$GRCh38)]] <- "hg18ToHg38.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh38, reference_builds$GRCh37)]] <- "hg38ToHg19.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh37, reference_builds$GRCh38)]] <- "hg19ToHg38.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh37, reference_builds$GRCh36)]] <- "hg19ToHg18.over.chain.gz"

#' Reasons a variant is dropped during liftover, in the order they are checked.
#' @keywords internal
liftover_drop_reasons <- c("unmapped", "non_primary_contig", "chrom_change", "strand_flip", "round_trip",
                           "duplicate_target")

#' convert_reference_build_via_liftover: Change reference build of BP marker from allow list of liftOver conversions
#'
#' Positions are lifted with UCSC liftOver, and variants whose new position can't be trusted are
#' dropped, each with a reason:
#' \itemize{
#'   \item \code{unmapped}: liftOver found no position;
#'   \item \code{non_primary_contig}: mapped to an alt, random or unplaced contig;
#'   \item \code{chrom_change}: mapped to another chromosome (paralogous or relocated sequence);
#'   \item \code{strand_flip}: mapped to the minus strand (sequence inverted between builds), so the
#'     alleles would need complementing.  With \code{complement_minus_strand = TRUE}, non-palindromic
#'     variants are kept with complemented alleles and only palindromic ones are dropped;
#'   \item \code{round_trip}: lifting the new position back does not return the original one (a
#'     conversion-unstable position, Ormond et al. 2021);
#'   \item \code{duplicate_target}: two different input positions land on the same variant.
#' }
#' A QC table (counts per reason) is attached as \code{attr(result, "liftover_qc")}, and when
#' \code{report_prefix} is given it is written to \verb{<report_prefix>.liftover_qc.tsv} with the dropped
#' rows in \verb{<report_prefix>.liftover_dropped.tsv.gz}.
#'
#' @param gwas GWAS (file or dataframe) of standardised GWAS
#' @param input_reference_build string reference build, found in reference_builds list
#' @param output_reference_build string reference build that GWAS is to change to, found in reference_builds list
#' @param output_file optional output file name to save to (the QC report is written beside it)
#' @param complement_minus_strand keep non-palindromic minus-strand variants, with complemented alleles
#' @param max_drop_fraction stop if more than this fraction of variants is dropped; a large loss
#'   usually means the declared input build is wrong
#' @param report_prefix optional path prefix for the QC report files
#' @return gwas input is altered and returned
#' @export
convert_reference_build_via_liftover <- function(gwas,
                                                 input_reference_build=reference_builds$GRCh37,
                                                 output_reference_build=reference_builds$GRCh38,
                                                 output_file,
                                                 complement_minus_strand=FALSE,
                                                 max_drop_fraction=0.05,
                                                 report_prefix=NULL) {
  gwas <- get_file_or_dataframe(gwas)

  if (input_reference_build == output_reference_build) {
    return(gwas)
  }
  start <- Sys.time()

  liftover_conversion <- available_liftover_conversions[[paste0(input_reference_build, output_reference_build)]]
  if (is.null(liftover_conversion)) {
    stop(paste(c("Error: liftOver combination of", input_reference_build, output_reference_build, "not recocognised.",
                  "Reference builds must be one of:", reference_builds), collapse = " "))
  }

  original_gwas_size <- nrow(gwas)
  gwas[[paste0("BP", sub("GRCh", "", input_reference_build))]] <- gwas$BP
  gwas$.ROW_ID <- seq_len(nrow(gwas))

  forward <- lift_positions(gwas$CHR, gwas$BP, gwas$.ROW_ID, input_reference_build, output_reference_build)
  matched <- match(gwas$.ROW_ID, forward$ROW_ID)
  new_chr <- forward$CHR[matched]
  new_bp <- forward$BP[matched]
  new_strand <- forward$STRAND[matched]

  reason <- rep(NA_character_, nrow(gwas))
  reason[is.na(new_bp)] <- "unmapped"
  undecided <- is.na(reason)
  reason[undecided & !is_primary_chromosome(new_chr)] <- "non_primary_contig"
  undecided <- is.na(reason)
  reason[undecided & chromosome_label(new_chr) != chromosome_label(gwas$CHR)] <- "chrom_change"

  undecided <- is.na(reason)
  minus <- undecided & new_strand == "-"
  if (complement_minus_strand) {
    palindromic <- is_palindromic(gwas$EA, gwas$OA)
    reason[minus & palindromic] <- "strand_flip"
    complement <- minus & !palindromic
    gwas$EA[complement] <- complement_alleles(gwas$EA[complement])
    gwas$OA[complement] <- complement_alleles(gwas$OA[complement])
  } else {
    reason[minus] <- "strand_flip"
  }

  undecided <- which(is.na(reason))
  if (length(undecided) > 0 && liftover_chain_available(output_reference_build, input_reference_build)) {
    back <- lift_positions(new_chr[undecided], new_bp[undecided], gwas$.ROW_ID[undecided],
                           output_reference_build, input_reference_build)
    back_matched <- match(gwas$.ROW_ID[undecided], back$ROW_ID)
    returns <- !is.na(back_matched) &
      chromosome_label(back$CHR[back_matched]) == chromosome_label(gwas$CHR[undecided]) &
      back$BP[back_matched] == gwas$BP[undecided]
    reason[undecided[!returns]] <- "round_trip"
  } else if (length(undecided) > 0) {
    message("No ", output_reference_build, " to ", input_reference_build,
            " chain available: skipping the round-trip check for conversion-unstable positions")
  }

  kept <- which(is.na(reason))
  if (length(kept) > 0) {
    target_key <- paste(chromosome_label(gwas$CHR[kept]), new_bp[kept],
                        pmin(toupper(gwas$EA[kept]), toupper(gwas$OA[kept])),
                        pmax(toupper(gwas$EA[kept]), toupper(gwas$OA[kept])))
    source_position <- paste(gwas$CHR[kept], gwas$BP[kept])
    sources_per_key <- tapply(source_position, target_key, function(x) length(unique(x)))
    reason[kept[sources_per_key[target_key] > 1]] <- "duplicate_target"
  }

  qc <- liftover_qc_table(reason, original_gwas_size, input_reference_build, output_reference_build,
                          liftover_conversion)
  dropped <- gwas[!is.na(reason), , drop = FALSE]
  dropped$LIFTOVER_DROP_REASON <- reason[!is.na(reason)]
  dropped$.ROW_ID <- NULL

  if (is.null(report_prefix) && !missing(output_file) && shiny::isTruthy(output_file)) {
    report_prefix <- file.path(dirname(output_file), file_prefix(output_file))
  }
  if (!is.null(report_prefix)) {
    data.table::fwrite(qc, paste0(report_prefix, ".liftover_qc.tsv"), sep = "\t")
    data.table::fwrite(dropped, paste0(report_prefix, ".liftover_dropped.tsv.gz"), sep = "\t")
  }

  keep <- is.na(reason)
  gwas$BP <- new_bp
  gwas <- gwas[keep, setdiff(colnames(gwas), ".ROW_ID"), drop = FALSE]

  n_dropped <- original_gwas_size - nrow(gwas)
  counts <- qc[qc$N > 0 & qc$REASON %in% liftover_drop_reasons, ]
  message(paste0("Liftover ", input_reference_build, " -> ", output_reference_build, ": kept ", nrow(gwas), " of ",
                 original_gwas_size, " variants; dropped ", n_dropped,
                 if (nrow(counts) > 0) paste0(" (", paste(counts$REASON, counts$N, sep = " ", collapse = ", "), ")")
                 else ""))
  if (n_dropped > max_drop_fraction * original_gwas_size) {
    stop(paste0("Liftover dropped ", n_dropped, " of ", original_gwas_size, " variants (more than ",
                100 * max_drop_fraction, "%). This usually means the GWAS is not in ", input_reference_build,
                ": check the declared build. Counts per reason: ",
                paste(counts$REASON, counts$N, sep = " ", collapse = ", ")))
  }
  attr(gwas, "liftover_qc") <- qc

  if(!missing(output_file) && shiny::isTruthy(output_file)) {
    vroom::vroom_write(gwas, output_file)
  }

  print("Liftover time taken:")
  print(Sys.time() - start)
  return(gwas)
}

#' Counts of kept and dropped variants per reason, as written to \verb{.liftover_qc.tsv}.
#' @keywords internal
liftover_qc_table <- function(reason, n_input, input_build, output_build, chain) {
  counts <- vapply(liftover_drop_reasons, function(r) sum(reason == r, na.rm = TRUE), numeric(1))
  data.frame(
    REASON = c("input", "kept", liftover_drop_reasons),
    N = c(n_input, sum(is.na(reason)), counts),
    INPUT_BUILD = input_build,
    OUTPUT_BUILD = output_build,
    CHAIN = chain,
    stringsAsFactors = FALSE
  )
}

#' Lift 1-based positions with UCSC liftOver.
#'
#' Each variant is written as the one-base BED interval \verb{[BP-1, BP)} with a row identifier and
#' \code{+} strand, so liftOver reports the strand of the target sequence.  The new 1-based position is
#' the end of the mapped interval.
#' @return data.frame with ROW_ID, CHR (without "chr"), BP and STRAND for the variants liftOver mapped
#' @keywords internal
lift_positions <- function(chr, bp, row_id, input_build, output_build) {
  bed_file_input <- tempfile(fileext = ".bed")
  bed_file_output <- tempfile(fileext = ".bed")
  unmapped <- tempfile(fileext = ".unmapped")
  on.exit(unlink(c(bed_file_input, bed_file_output, unmapped)), add = TRUE)

  create_bed_file_from_positions(chr, bp, row_id, bed_file_input)
  run_liftover(bed_file_input, bed_file_output, input_build, output_build, unmapped)
  return(read_lifted_bed_file(bed_file_output))
}

#' Write the liftOver input BED (6 columns: chrom, start, end, name, score, strand).
#' @keywords internal
create_bed_file_from_positions <- function(chr, bp, row_id, output_file) {
  bed <- data.frame(
    CHROM = paste0("chr", ucsc_chromosome_name(chr)),
    START = format(as.numeric(bp) - 1, scientific = FALSE, trim = TRUE),
    END = format(as.numeric(bp), scientific = FALSE, trim = TRUE),
    NAME = format(row_id, scientific = FALSE, trim = TRUE),
    SCORE = 0,
    STRAND = "+",
    stringsAsFactors = FALSE
  )
  data.table::fwrite(bed, output_file, sep = "\t", col.names = FALSE, quote = FALSE)
  return(invisible(bed))
}

#' Read liftOver output written from \code{create_bed_file_from_positions()}.
#' @keywords internal
read_lifted_bed_file <- function(bed_file) {
  if (!file.exists(bed_file) || file.size(bed_file) == 0) {
    return(data.frame(ROW_ID = integer(0), CHR = character(0), BP = numeric(0), STRAND = character(0)))
  }
  lifted <- data.table::fread(bed_file, header = FALSE, sep = "\t", colClasses = "character")
  data.frame(
    ROW_ID = as.integer(lifted[[4]]),
    CHR = sub("^chr", "", lifted[[1]]),
    BP = as.numeric(lifted[[3]]),
    STRAND = if (ncol(lifted) >= 6) lifted[[6]] else "+",
    stringsAsFactors = FALSE
  )
}

#' TRUE when the chain file for converting \code{input_build} to \code{output_build} is installed.
#' @keywords internal
liftover_chain_available <- function(input_build, output_build) {
  conversion <- available_liftover_conversions[[paste0(input_build, output_build)]]
  return(!is.null(conversion) && file.exists(file.path(liftover_chain_directory(), conversion)))
}

#' Directory holding the liftOver chain files: \code{PIPELINE_DATA_DIR/genomic_data/liftover}, unless
#' \code{options(GeneHackman.liftover_chain_dir = ...)} points elsewhere.
#' @keywords internal
liftover_chain_directory <- function() {
  getOption("GeneHackman.liftover_chain_dir", liftover_chain_dir)
}

#' The liftOver binary: on the PATH or \code{/usr/local/bin/liftOver}, unless
#' \code{options(GeneHackman.liftover_binary = ...)} points elsewhere.
#' @keywords internal
liftover_binary_path <- function() {
  getOption("GeneHackman.liftover_binary", liftover_binary)
}

#' Chromosome names as UCSC chain files use them (X, Y, M), from GWAS codes (23, 24, 25, 26, MT, ...).
#' @keywords internal
ucsc_chromosome_name <- function(chr) {
  chr <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  chr[chr %in% c("23", "25")] <- "X"
  chr[chr == "24"] <- "Y"
  chr[chr %in% c("26", "MT")] <- "M"
  return(chr)
}

#' Chromosome label used to compare chromosomes across builds (23 and X are the same chromosome).
#' @keywords internal
chromosome_label <- function(chr) {
  return(ucsc_chromosome_name(chr))
}

#' TRUE for chromosomes 1-22, X, Y and M; FALSE for alt, random and unplaced contigs.
#' @keywords internal
is_primary_chromosome <- function(chr) {
  return(!is.na(chr) & ucsc_chromosome_name(chr) %in% c(as.character(1:22), "X", "Y", "M"))
}

#' TRUE for A/T and C/G variants, whose strand can't be told from the alleles.
#' @keywords internal
is_palindromic <- function(ea, oa) {
  return(toupper(ea) == complement_alleles(oa))
}

#' Complement DNA alleles (A<->T, C<->G), keeping the base order; other characters are left alone.
#' @keywords internal
complement_alleles <- function(alleles) {
  return(chartr("ACGTacgt", "TGCAtgca", toupper(alleles)))
}


#' Thin wrapper so tests can mock via testthat::local_mocked_bindings().
#' @keywords internal
run_system <- function(command, wait = TRUE, intern = FALSE,
                         ignore.stdout = FALSE, ignore.stderr = FALSE) {
  base::system(command, wait = wait, intern = intern,
               ignore.stdout = ignore.stdout, ignore.stderr = ignore.stderr)
}

run_liftover <- function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
  liftover_conversion <- available_liftover_conversions[[paste0(input_build, output_build)]]
  chain_file <- file.path(liftover_chain_directory(), liftover_conversion)
  if (!file.exists(chain_file)) {
    stop("liftOver chain file not found: ", chain_file, ". Download the reference data (see README.md).")
  }
  liftover_command <- paste(shQuote(liftover_binary_path()), shQuote(bed_file_input), shQuote(chain_file),
                            shQuote(bed_file_output), shQuote(unmapped))
  status <- run_system(liftover_command, wait = TRUE, ignore.stderr = TRUE)
  if (!is.null(status) && status != 0) {
    stop("liftOver failed (exit code ", status, "): ", liftover_command)
  }
}


#' Check that a GWAS's positions are in the build it is declared to be in.
#'
#' Builds \code{CHR:BP} keys with alphabetically ordered alleles for one chromosome of the GWAS and
#' counts how many match the Pan-UKB LD pack's GRCh37 keys and its GRCh38 keys (the pack lists about
#' 10 million common variants per ancestry in both builds).  If the other build matches clearly
#' better, the declared build is wrong and standardisation stops, instead of lifting from the wrong
#' build and silently producing wrong positions.  Skipped (with a message) when the pack is not
#' installed or the build is not GRCh37/GRCh38.
#'
#' @param gwas data frame with CHR, BP, EA and OA
#' @param declared_build the build the GWAS is declared to be in
#' @param pack_dir directory of the Pan-UKB LD pack
#' @param max_rows number of rows sampled from the chromosome with the most rows
#' @return invisibly, a named vector with the share of sampled keys found in each build
#' @keywords internal
check_declared_build <- function(gwas, declared_build, pack_dir = default_ld_pack_dir(), max_rows = 20000) {
  builds <- c(reference_builds$GRCh37, reference_builds$GRCh38)
  if (!declared_build %in% builds || nrow(gwas) == 0) return(invisible(NULL))
  population <- intersect(c("EUR", ld_pack_ancestries), list.dirs(pack_dir, full.names = FALSE, recursive = FALSE))
  if (length(population) == 0) {
    message("Build check skipped: no Pan-UKB LD pack at ", pack_dir)
    return(invisible(NULL))
  }
  chr_label <- chromosome_label(gwas$CHR)
  autosomal <- chr_label %in% as.character(1:22)
  if (!any(autosomal)) return(invisible(NULL))
  chrom <- names(which.max(table(chr_label[autosomal])))
  variants_file <- file.path(pack_dir, population[1], paste0("variants_chr", chrom, ".tsv.gz"))
  if (!file.exists(variants_file)) {
    message("Build check skipped: ", variants_file, " not found")
    return(invisible(NULL))
  }

  rows <- which(chr_label == chrom)
  if (length(rows) > max_rows) rows <- rows[round(seq(1, length(rows), length.out = max_rows))]
  ea <- toupper(gwas$EA[rows])
  oa <- toupper(gwas$OA[rows])
  keys <- paste0(chrom, ":", format(as.numeric(gwas$BP[rows]), scientific = FALSE, trim = TRUE), "_",
                 pmin(ea, oa), "_", pmax(ea, oa))
  pack <- data.table::fread(variants_file, select = c("SNP37", "SNP38"), colClasses = "character",
                            na.strings = c("", "NA"))
  share <- c(GRCh37 = mean(keys %in% pack$SNP37), GRCh38 = mean(keys %in% pack$SNP38))
  other_build <- setdiff(builds, declared_build)
  message(sprintf("Build check on chr%s (%d variants): %.1f%% match GRCh37 keys, %.1f%% match GRCh38 keys",
                  chrom, length(rows), 100 * share[["GRCh37"]], 100 * share[["GRCh38"]]))
  if (share[[other_build]] > 0.2 && share[[other_build]] > 3 * share[[declared_build]]) {
    stop(sprintf(paste0("The GWAS is declared as %s, but its positions match %s much better (%.1f%% vs %.1f%% ",
                        "of variants on chr%s match the Pan-UKB LD pack). Set build: %s for this GWAS."),
                 declared_build, other_build, 100 * share[[other_build]], 100 * share[[declared_build]],
                 chrom, other_build))
  }
  return(invisible(share))
}
