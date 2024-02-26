#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011
source("../R/util.R")
source("../R/constants.R")
lapply(list.files("../R/", full.names = T, pattern="\\.R$"), source)

options(error = function() traceback(20))
