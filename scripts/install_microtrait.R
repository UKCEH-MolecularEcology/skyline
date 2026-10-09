#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait
# EI notes:
#  - R's write check (file.access) can misreport /ei as read-only -> install into a staging
#    library on local /tmp, then copy into the env library
#  - microtrait's DESCRIPTION declares no dependencies -> install them explicitly (list from its README/NAMESPACE)
#  - GitHub packages are fetched as codeload tarballs (no API, no PAT needed)

############################## LOG
log <- file(snakemake@log[["out"]], open = "wt")
sink(log); sink(log, type = "message")

Sys.setenv(GITHUB_PAT = "", GITHUB_TOKEN = "")
options(repos = BiocManager::repositories(), Ncpus = 4)

env_lib <- .libPaths()[1]
stage   <- file.path("/tmp", paste0("rlib_microtrait_", Sys.getenv("USER"), "_", Sys.getpid()))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(stage, .libPaths()))
message("env lib: ", env_lib, "\nstaging lib: ", stage)

missing_pkgs <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

# 1. CRAN + Bioconductor dependencies (microtrait README + NAMESPACE imports; gRodon imports)
deps <- c("R.utils", "RColorBrewer", "ape", "assertthat", "checkmate", "corrplot", "doParallel",
          "dplyr", "foreach", "fs", "futile.logger", "ggplot2", "gtools", "kmed", "lazyeval",
          "magrittr", "nlme", "pheatmap", "readr", "seqinr", "stringr", "tibble", "tictoc",
          "tidyr", "vegan", "purrr", "matrixStats", "testthat", "ggsci",
          "Biostrings", "coRdon", "ComplexHeatmap")
todo <- missing_pkgs(deps)
message("Installing ", length(todo), " dependencies: ", paste(todo, collapse = ", "))
if (length(todo)) BiocManager::install(todo, lib = stage, update = FALSE, ask = FALSE)
still <- missing_pkgs(deps)
if (length(still)) stop("Dependencies failed to install: ", paste(still, collapse = ", "))

# 2. GitHub packages as tarballs: gRodon (imported by microtrait), then microtrait
gh_install <- function(repo) {
  tb <- file.path(tempdir(), paste0(basename(repo), ".tar.gz"))
  download.file(paste0("https://codeload.github.com/", repo, "/tar.gz/HEAD"), tb, mode = "wb")
  install.packages(tb, repos = NULL, type = "source", lib = stage)
}
if (!requireNamespace("gRodon", quietly = TRUE))     gh_install("jlw-ecoevo/gRodon")
if (!requireNamespace("microtrait", quietly = TRUE)) gh_install("ukaraoz/microtrait")

# 3. copy staged packages into the env library
pkgs <- list.dirs(stage, recursive = FALSE, full.names = TRUE)
message("Copying ", length(pkgs), " packages into the env library: ", paste(basename(pkgs), collapse = ", "))
for (p in pkgs) {
  dest <- file.path(env_lib, basename(p))
  if (dir.exists(dest)) unlink(dest, recursive = TRUE)
  stopifnot(file.copy(p, env_lib, recursive = TRUE))
}
.libPaths(env_lib)

# 4. load + HMM models (downloaded from GitHub releases / dbCAN, written inside the package)
library(Biostrings)
library(microtrait)
microtrait::prep.hmmmodels()

message("Biostrings ", as.character(packageVersion("Biostrings")),
        "; microtrait ", as.character(packageVersion("microtrait")))
unlink(stage, recursive = TRUE)   # only after everything is loaded
file.create(snakemake@output[[1]])
