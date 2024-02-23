install.packages(c("devtools", "argparser", "gmp", "corrplot", "broom", "conflicted", "nloptr", "Cairo"),
                 repos = "http://cran.us.r-project.org")

#also needed: forestplot, data.table, devtools
#install.packages('https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_23.8.tar.gz', repos = NULL, type = 'source')
.libPaths( c( .libPaths(), "/usr/lib/R/site-library") )

devtools::install_github(c("MRCIEU/gwaspipeline", "phenoscanner/phenoscanner", "suchestoncampbelllab/gwasurvivr"))
