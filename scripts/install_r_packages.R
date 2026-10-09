#!/usr/bin/env Rscript
# Install extra CRAN packages into the R library of the conda env running this Rscript.
# EI: run on the software node (internet). R's write check misreports /ei as read-only, so packages
# are installed into a staging library on local /tmp and then copied into the env library.
# Usage: <env>/bin/Rscript scripts/install_r_packages.R ggsci [more packages ...]
pkgs <- commandArgs(trailingOnly = TRUE)
if (length(pkgs) == 0) pkgs <- "ggsci"
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
env_lib <- .libPaths()[1]
stage <- file.path("/tmp", paste0("rlib_", Sys.getenv("USER"), "_", Sys.getpid()))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(stage, .libPaths()))
todo <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(todo) > 0) install.packages(todo, lib = stage)
for (p in list.dirs(stage, recursive = FALSE, full.names = TRUE)) {
  dest <- file.path(env_lib, basename(p))
  if (dir.exists(dest)) unlink(dest, recursive = TRUE)
  stopifnot(file.copy(p, env_lib, recursive = TRUE))
}
unlink(stage, recursive = TRUE)
.libPaths(env_lib)
for (p in pkgs) cat(p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED", "\n")
