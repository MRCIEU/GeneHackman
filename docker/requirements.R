options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("/usr/lib/R/site-library", .libPaths()))

install.packages("remotes")
pkgs <- c("argparser", "gmp", "corrplot", "broom", "conflicted", "nloptr", "Cairo", "susieR")
install.packages(pkgs)

install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v75"), update = FALSE, ask = FALSE)

install.packages(c("locuszoomr", "patchwork"))

remotes::install_github(c("MRCIEU/GeneHackman"), upgrade = "never")
