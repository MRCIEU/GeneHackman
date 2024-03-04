source("../R/util.R")
source("../R/constants.R")
lapply(list.files("../R/", full.names = T, pattern="\\.R$"), source)

options(error = function() traceback(20))
