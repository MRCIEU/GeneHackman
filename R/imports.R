# NAMESPACE: consolidated imports with `except` to avoid R CMD check WARNINGs from
# symbol clashes between e.g. dplyr/stats, data.table/rlang, TwoSampleMR/MR, R.utils, …
# (see 00install.out: "Warning: replacing previous import … when loading 'GeneHackman'").
# All @import and @rawNamespace must be in the same roxygen block before the following NULL
# (otherwise roxygen2 drops them from NAMESPACE).
#
#' @keywords internal
#' @rawNamespace import(MendelianRandomization, except = c(mr_ivw, mr_median))
#' @rawNamespace import(SlopeHunter, except = c(format_data))
#' @rawNamespace import(data.table, except = c(between, first, last))
#' @rawNamespace import(dplyr, except = c(between, first, last, filter, lag))
#' @rawNamespace import(stats, except = c(filter, lag))
#' @rawNamespace import(rlang, except = c(`:=`))
#' @rawNamespace import(R.utils, except = c(env, extract, reset, setProgress, timestamp, validate))
#' @rawNamespace import(future, except = c(run))
#' @rawNamespace import(shiny, except = c(setProgress, validate))
#' @rawNamespace import(tidyr, except = c(extract))
#' @rawNamespace import(utils, except = c(timestamp))
#' @import coloc
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
#' @import grDevices
#' @import graphics
#' @import httr
#' @import prettydoc
#' @import qqman
#' @import rmarkdown
#' @import stringr
#' @import tibble
#' @import vroom
#' @import susieR
#' @import TwoSampleMR
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom furrr future_map
#' @importFrom fst fst
#' @importFrom fst read_fst
#' @importFrom progressr progressor
#' @importFrom tibble as_tibble
NULL
