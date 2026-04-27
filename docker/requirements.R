options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("/usr/lib/R/site-library", .libPaths()))

install.packages("remotes")
pkgs <- c("argparser", "gmp", "corrplot", "broom", "conflicted", "nloptr", "Cairo", "susieR")
install.packages(pkgs)

install.packages("devtools")

remotes::install_github(c("MRCIEU/GeneHackman"), upgrade = "never")
