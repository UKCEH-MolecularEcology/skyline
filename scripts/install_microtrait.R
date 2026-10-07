#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait

############################## LOG
log <- file(snakemake@log[["out"]], open = "wt")
sink(log); sink(log, type = "message")

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("microtrait", quietly = TRUE)) {
  remotes::install_github("ukaraoz/microtrait", dependencies = TRUE, upgrade = "never")
}
library(microtrait)
microtrait::prep.hmmmodels()

message("microtrait installed: ", as.character(packageVersion("microtrait")))
file.create(snakemake@output[[1]])
