options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("/usr/lib/R/site-library", .libPaths()))

install.packages("remotes")
pkgs <- c("argparser", "gmp", "corrplot", "broom", "conflicted", "nloptr", "Cairo", "Rfast", "susieR")
install.packages(pkgs)

install.packages("devtools")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# gene annotation for LocusZoom plots: Ensembl 86 for GRCh38 output, Ensembl 75 for GRCh37
BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86", "EnsDb.Hsapiens.v75"), update = FALSE, ask = FALSE)

install.packages(c("locuszoomr", "patchwork"))

remotes::install_deps(
  "docker",
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never",
  repos = BiocManager::repositories()
)
