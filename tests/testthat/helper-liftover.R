# A stand-in for UCSC liftOver in unit tests: maps every position by a fixed shift (forward builds add
# `shift`, the reverse direction subtracts it), and lets a test override single rows by their BED
# name (the row identifier convert_reference_build_via_liftover() writes).
mock_liftover_shift <- function(shift = 1000, forward = c("GRCh37", "GRCh38"), overrides = list(),
                                reverse_overrides = list()) {
  function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
    bed <- utils::read.table(bed_file_input, sep = "\t", colClasses = "character",
                             col.names = c("chrom", "start", "end", "name", "score", "strand"))
    is_forward <- identical(c(input_build, output_build), forward)
    delta <- if (is_forward) shift else -shift
    rows <- list()
    rules <- if (is_forward) overrides else reverse_overrides
    for (i in seq_len(nrow(bed))) {
      rule <- rules[[bed$name[i]]]
      if (identical(rule, "unmapped")) next
      chrom <- bed$chrom[i]
      end <- as.numeric(bed$end[i]) + delta
      strand <- "+"
      if (is.list(rule)) {
        if (!is.null(rule$chrom)) chrom <- rule$chrom
        if (!is.null(rule$end)) end <- rule$end
        if (!is.null(rule$strand)) strand <- rule$strand
      }
      rows[[length(rows) + 1]] <- data.frame(chrom, format(end - 1, scientific = FALSE),
                                             format(end, scientific = FALSE), bed$name[i], 0, strand)
    }
    out <- if (length(rows)) do.call(rbind, rows) else data.frame()
    utils::write.table(out, bed_file_output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    file.create(unmapped)
  }
}

# Real UCSC liftOver with GeneHackman's chain files, when GH_LIFTOVER_CHAIN_DIR points at them (and
# liftOver is on the PATH or in GH_LIFTOVER_BIN); otherwise the test is skipped.
local_real_liftover <- function(env = parent.frame()) {
  chain_dir <- Sys.getenv("GH_LIFTOVER_CHAIN_DIR")
  binary <- Sys.getenv("GH_LIFTOVER_BIN", Sys.which("liftOver"))
  if (!nzchar(chain_dir) || !file.exists(file.path(chain_dir, "hg19ToHg38.over.chain.gz")) ||
      !file.exists(file.path(chain_dir, "hg38ToHg19.over.chain.gz")) || !nzchar(binary) || !file.exists(binary)) {
    testthat::skip("set GH_LIFTOVER_CHAIN_DIR (and GH_LIFTOVER_BIN) to run liftover tests with real chains")
  }
  withr::local_options(list(GeneHackman.liftover_chain_dir = chain_dir, GeneHackman.liftover_binary = binary),
                       .local_envir = env)
}
