source("../R/util.R")
source("../R/constants.R")
lapply(list.files("../R/", full.names = T, pattern="\\.R$"), source)

# Prevent stray Rplots.pdf from being written in batch mode (all pipelines).
if (!interactive()) pdf(file = nullfile())

options(error = function() traceback(20))
