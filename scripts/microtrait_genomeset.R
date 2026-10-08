#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait -- genome-set matrices
# NOTE: make.genomeset.results(growthrate = T / optimumT = T) fails when any genome lacks a gRodon
#       estimate (growthrate_d / ogt is empty: "arguments imply differing number of rows: 1, 0").
#       So the matrices are built without them and mingentime / optimumT are added here, NA if missing.

############################## LOG
log <- file(snakemake@log[["out"]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(microtrait))
library(tictoc)

rds_files <- unlist(snakemake@params[["rds"]])
missing <- rds_files[!file.exists(rds_files)]
if (length(missing) > 0) stop(length(missing), " rds files missing, e.g. ", missing[1])
ids <- sub("\\.microtrait\\.rds$", "", basename(rds_files))
ncores <- snakemake@threads
message("Number of rds files: ", length(rds_files))

tictoc::tic("make.genomeset.results")
genomeset_results <- make.genomeset.results(rds_files = rds_files, ids = ids,
                                            growthrate = FALSE, optimumT = FALSE, ncores = ncores)
tictoc::toc()

# mingentime (gRodon; microtrait stores it as growthrate_d, in hours) and optimum growth temperature
get_value <- function(f, field) {
  x <- readRDS(f)[[field]]
  if (is.null(x) || length(x) == 0) NA_real_ else suppressWarnings(as.numeric(x[[1]]))
}
extra <- tibble::tibble(
  id         = ids,
  mingentime = unlist(parallel::mclapply(rds_files, get_value, field = "growthrate_d", mc.cores = ncores)),
  optimumT   = unlist(parallel::mclapply(rds_files, get_value, field = "ogt",          mc.cores = ncores))
)
message("Genomes without a mingentime estimate: ", sum(is.na(extra$mingentime)),
        "; without optimumT: ", sum(is.na(extra$optimumT)))
for (n in names(genomeset_results)) {
  genomeset_results[[n]] <- dplyr::left_join(genomeset_results[[n]], extra, by = "id")
}

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
