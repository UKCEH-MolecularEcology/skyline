#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait -- genome-set matrices

############################## LOG
log <- file(snakemake@log[["out"]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(microtrait))
library(tictoc)

rds_files <- unlist(snakemake@params[["rds"]])
missing <- rds_files[!file.exists(rds_files)]
if (length(missing) > 0) stop(length(missing), " rds files missing, e.g. ", missing[1])
ids <- sub("\\.microtrait\\.rds$", "", basename(rds_files))
message("Number of rds files: ", length(rds_files))

tictoc::tic("make.genomeset.results")
genomeset_results <- make.genomeset.results(rds_files = rds_files, ids = ids,
                                            ncores = snakemake@threads)
tictoc::toc()

# saving genome set results
saveRDS(genomeset_results, file = snakemake@output[["rds"]])

tables <- snakemake@output[["tables"]]
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
for (n in names(genomeset_results)) {
  x <- genomeset_results[[n]]
  if (is.data.frame(x) || is.matrix(x)) {
    write.table(as.data.frame(x), file.path(tables, paste0(n, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote ", n)
  }
}
message("Microtrait analysis completed successfully.")
