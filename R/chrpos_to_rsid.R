#' @title Chromosome & position data to variant RSID
#'
#' @param dt a data.frame like object, or file path, with at least columns (chrom<chr>, pos<int>, ea<chr>, nea<chr>)
#' @param chr_col a string column name; chromosome position
#' @param pos_col a string column name; base position
#' @param ea_col a string column name; effect allele
#' @param nea_col a string column name; non effect allele
#' @param build a string, options: "b37_dbsnp156", "b38_dbsnp156" (corresponds to the appropriate data directory)
#' @param flip a string, options: "report", "allow", "no_flip"
#' @param alt_rsids a logical, whether to return additional alternate RSIDs
#' @param verbose a logical, runtime reporting
#' @param dbsnp_dir a string file path to the dbSNP .fst file directory - see setup documentation
#' @param parallel_cores an integer, the number of cores/workers to set up the `future::multisession` with
#'
#' @return a tibble with an RSID column, or (if `alt_rsids` is `TRUE`) a list with
#'   tibbles: `data` and `alt_rsids`
#' @export
#' @importFrom data.table setkey setkeyv
#' @importFrom progressr progressor
#' @importFrom furrr future_map
#' @importFrom fst fst read_fst
#' @importFrom tibble as_tibble
#'
chrpos_to_rsid <- function(dt,
                           chr_col,
                           pos_col,
                           ea_col=NULL,
                           nea_col=NULL,
                           flip="allow",
                           alt_rsids=FALSE,
                           build="b37_dbsnp156",
                           dbsnp_dir = NA_character_,
                           parallel_cores = parallel::detectCores(),
                           verbose = TRUE) {

  #-- Checks -----------------------------------------------------------

  # import / convert input
  dt <- data.table::as.data.table(dt)

  stopifnot("`chr_col` must be a column name in `dt`" = chr_col %in% colnames(dt))
  stopifnot("`pos_col` must be a column name in `dt`" = pos_col %in% colnames(dt))
  stopifnot("`pos_col` must be numeric/integer type" = is.numeric(dt[,get(pos_col)]))
  stopifnot("`ea_col` and `nea_col` must be either both provided or both NULL" = sum(c(is.null(nea_col),is.null(ea_col)))%in%c(0,2))
  if(all(!c(is.null(ea_col), is.null(nea_col)))) {
    stopifnot("`ea_col` must be a column name in `dt`" = ea_col %in% colnames(dt))
    stopifnot("`nea_col` must be a column name in `dt`" = nea_col %in% colnames(dt))
    stopifnot("`ea_col` column must be character type" = is.character(dt[,get(ea_col)]))
    stopifnot("`nea_col` column must be character type" = is.character(dt[,get(nea_col)]))
  }
  stopifnot("`chr_col` column must be character, integer, or numeric type" = is.character(dt[,get(chr_col)]) |
                                                                             is.integer(dt[,get(chr_col)]) |
                                                                             is.numeric(dt[,get(chr_col)]))
  if(is.numeric(dt[,get(chr_col)])) {
    stopifnot("`chr_col`, if numeric, must be whole numbers" = all(dt[,get(chr_col)] %% 1 == 0))
  }
  stopifnot("`dbsnp_dir` must be a valid directory path" = dir.exists(dbsnp_dir))
  available_builds <- list.dirs(dbsnp_dir, full.names=FALSE)[-1]
  build <- match.arg(build, available_builds)
  flip  <- match.arg(flip, c("report", "allow", "no_flip"))


  #-- Processing ---------------------------------------------------------

  # data adjustments
  # 1) remove chromosome 'Chr6' --> '6' and recode X/Y
  #    (?i) case-insensitive
  #    ^ at the beginning of the string
  #    (?:chr)? non-capture group, matches "chr" 0-1 times
  #    ([0-9XYMT]+) capture group //1, matches 1 or more occurances of group [0-9XYMT]
  #    \\1 replace with string found in capture group 1
  chr_col_orig = "chr_col_orig"
  dt[, (chr_col_orig) := get(chr_col)]
  dt[, (chr_col) := sub("(?i)^(?:chr)?([0-9XY]+)", toupper("\\1"), toupper(get(chr_col)))]
  # 2) replace X[23,25] and Y[24]
  dt[, (chr_col) := ifelse(get(chr_col) %in% c("X","25") ,"23",ifelse(get(chr_col)=="Y","24",get(chr_col)))]
  # 3) delete the RSID column ahead of re-adding it, if it exists
  if("RSID" %in% names(dt)) {
    dt[, "RSID" := NULL]
  }

  # split the dt into a list of data.tables, by chromosome
  dts <- split(dt, by=chr_col)

  if(verbose) {
    message(paste("RSID mapping...\nChromosomes to process: ", paste0(names(dts), collapse=", ")))
  }

  # set up parallel processing and process each chromosome independently
  future::plan(future::multisession, workers=parallel_cores)
  progressr::with_progress({

    # set up progress bar
    p <- progressr::progressor(steps = 7*length(dts)) # going to update 7 times in the function

    # run each chromosome in parallel
    # out[[1]] is the RSID data
    # out[[2]] is the alternative RSID data (other RSIDs for the same CHR:BP:A1:A2)
    out <- furrr::future_map(.x = dts,
                             .f = process_chromosome,
                             chr_col=chr_col, pos_col=pos_col, nea_col=nea_col, ea_col=ea_col, build=build, dbsnp_dir=dbsnp_dir, flip=flip, alt_rsids=alt_rsids, p=p,
                             .options=furrr::furrr_options(seed=TRUE))

  })

  future::plan(future::sequential)
  rm(dts)
  gc(verbose = FALSE)

  # deal with alternative RSIDs if requested
  if(alt_rsids) {

    out_data <- data.table::rbindlist( lapply(out, `[[`, 1) )
    out_alts <- data.table::rbindlist( lapply(out, `[[`, 2) )

  } else {
    out_data <- data.table::rbindlist(out)
  }
  rm(out)
  gc(verbose = FALSE)

  # put back the original chromosome column and remove the temporary one
  out_data[, (chr_col) := chr_col_orig]
  out_data[, chr_col_orig := NULL]

  # RSID to the front
  data.table::setcolorder(out_data, c("RSID", names(out_data)[names(out_data) != "RSID"]))

  # report
  msg <- paste0("RSID coverage ", round(100*sum(!is.na(out_data[["RSID"]]))/nrow(out_data), digits=2), "% (", sum(!is.na(out_data[["RSID"]])), "/", nrow(out_data), ")")
  if(flip!="no_flip") {
    msg <- paste0(msg, ", of which ", round(100*sum(out_data[["rsid_flip_match"]], na.rm=TRUE)/nrow(out_data), digits=2),"% (", sum(out_data[["rsid_flip_match"]], na.rm=TRUE), "/", nrow(out_data), ") where found flipping alleles")
  }
  message(msg)

  # remove flipping flag from output if not wanted
  if(flip=="allow") {
    out_data[, "rsid_flip_match" := NULL]
  }

  # return as same type as input
  # data.table::setattr(out_data, "orig_type", attr(dt,"orig_type"))
  if (alt_rsids) {
    return(list(
      "data" = tibble::as_tibble(out_data),
      "alt_rsids" = tibble::as_tibble(out_alts)
    ))
  }
  return(tibble::as_tibble(out_data))
}


# this function could be defined in the above function. However I took it out to avoid
# the gwas data.table being captured in each futures environment. See this discussion
# for the details: https://furrr.futureverse.org/articles/gotchas.html#function-environments-and-large-objects
process_chromosome <- function(chrom_dt, chr_col, pos_col, build, dbsnp_dir, flip, alt_rsids, p, nea_col=NULL, ea_col=NULL) {
  message(paste0("Processing chromosome ", chrom_dt[[1, chr_col]]))

  # silence RMDcheck warning
  RSID = i.RSID = baseRSID = rsid_flip_match = REF = ALT = NULL

  # increment progress bar #1
  p()

  # the chromosome to process, a string value e.g. "1"
  chrom <- chrom_dt[[1, chr_col]]

  # the dbSNP `.fst` file for this chromosome
  dbSNP_build_dir  <- file.path(dbsnp_dir, build)
  dbSNP_path <- file.path(dbSNP_build_dir, paste0("chr", chrom, ".fst"))
  stopifnot("Missing .fst reference file - is the dbSNP directory set correctly?" = file.exists(dbSNP_path))

  # read just the dbSNP position data (an integer vector) (reduce amount of data read in)
  dbSNP_pos <- fst::read_fst(dbSNP_path, "BP", as.data.table=TRUE)

  # increment progress bar #2
  p()

  # set as keys
  data.table::setkey(dbSNP_pos, "BP")
  data.table::setkeyv(chrom_dt, pos_col)

  # get the indices of the positions that are needed / also found in the input data
  dbSNP_pos[chrom_dt, "found" := get(pos_col)]
  row_idxs <- which(!is.na(dbSNP_pos[["found"]]))
  rm(dbSNP_pos)

  # define whether to work with alleles, or just CHR:POS
  alleles <- !(is.null(nea_col) & is.null(ea_col))

  # define the key depending on which columns are provided
  if(alleles) {
    dbSNP_key     <- c("CHR","BP","REF","ALT")
    dbSNP_keyflip <- c("CHR","BP","ALT","REF")
    data_key      <- c(chr_col, pos_col, nea_col, ea_col)
  } else {
    dbSNP_key <- c("CHR","BP")
    data_key  <- c(chr_col, pos_col)
  }

  # if matches found
  if(length(row_idxs)>0) {

    # increment progress bar #3
    p()

    # create a fst object which allows row access without reading the whole file
    dbSNP_fst <- fst::fst(dbSNP_path)

    # read the needed rows
    dbSNP_data <- dbSNP_fst[row_idxs, c("RSID", dbSNP_key)] |> data.table::as.data.table()
    rm(dbSNP_fst)

    # increment progress bar #4
    p()

    if(alleles) {
      # split the ALT column which can be a comma separate vector of alternative alleles
      dbSNP_data[, "ALT" := lapply(.SD, strsplit, split=','), .SDcols="ALT"]

      # make data.table longer, one allele combination per row
      dbSNP_data <- dbSNP_data[, lapply(.SD, unlist), by=1:nrow(dbSNP_data)]
      dbSNP_data[, nrow := NULL]

      # recode as D/I if input data has D/I coding
      if(any(grepl("^(D|I)$", chrom_dt[[nea_col]]))) {
        dbSNP_data[, "ALT" := data.table::fcase(nchar(ALT)< nchar(REF), "D",
                                                 nchar(ALT)> nchar(REF), "I",
                                                 nchar(ALT)==nchar(REF), ALT)]
        dbSNP_data[, "REF" := data.table::fcase(ALT=="D", "I",
                                                ALT=="I", "D",
                                                rep(TRUE, nrow(dbSNP_data)), REF)]
      }
    }

    # increment progress bar #5
    p()

    if(alt_rsids) {
      # get the alt rsids (some rsids code for the same position, chr, alt, and ref...)
      alt_rsid_data  <- dbSNP_data[duplicated(dbSNP_data[, dbSNP_key, with = FALSE]),]
    }

    # take unique (first rsid occurance)
    dbSNP_data <- dbSNP_data[!duplicated(dbSNP_data[, dbSNP_key, with = FALSE]),]

    # set the keys to match expected way round
    data.table::setkeyv(dbSNP_data, dbSNP_key)
    data.table::setkeyv(chrom_dt, data_key)
    chrom_dt[dbSNP_data, "RSID" := i.RSID]

    if(alt_rsids) {
      # set key
      data.table::setkeyv(alt_rsid_data, dbSNP_key)

      # map the chosen RSID to the alternative RSIDS
      alt_rsid_data[chrom_dt, "baseRSID" := i.RSID]
    }

    # increment progress bar #6
    p()

    # flip the alleles
    if(flip!="no_flip" & alleles) {

      # set the keys to flipped alleles
      data.table::setkeyv(dbSNP_data, dbSNP_keyflip)

      # add in flipped matches and a logical flag
      chrom_dt[dbSNP_data, c("RSID", "rsid_flip_match") := list(i.RSID, TRUE)]
      chrom_dt[!is.na(RSID) & is.na(rsid_flip_match), rsid_flip_match := FALSE]

      if(alt_rsids) {
        # get any alt rsids that match flipped
        data.table::setkeyv(alt_rsid_data, dbSNP_keyflip)
        alt_rsid_data[chrom_dt, "baseRSID" := i.RSID]
      }

    }

  # no matches found. skip all the processing but add the columns
  } else {
    p()
    p()
    p()
    p()
    chrom_dt[["RSID"]] <- NA_character_
    if(alt_rsids) alt_rsid_data <- data.table::copy(chrom_dt)
    if(flip!="no_flip" & alleles) chrom_dt[, rsid_flip_match := NA]
    if(alt_rsids) {
      alt_rsid_data[, "baseRSID" := NA_character_]
      alt_rsid_data <- alt_rsid_data[!is.na(baseRSID), ]
    }
  }

  # increment progress bar #7
  p()

  # return
  if(alt_rsids) {
    alt_rsid_data <- alt_rsid_data[is.na(baseRSID), ]
    return(list("data"=chrom_dt, "alt_rsids"=alt_rsid_data))
  } else {
    return(chrom_dt)
  }

}