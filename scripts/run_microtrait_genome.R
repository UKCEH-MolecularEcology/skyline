#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait -- one genome per job

############################## LOG
log <- file(snakemake@log[["out"]], open = "wt")
sink(log); sink(log, type = "message")

############################## LIBS
suppressPackageStartupMessages(library(microtrait))

# Biostrings set ops (masking workaround, as in rehab_analysis' run_microtrait.R)
for (f in c("setequal", "intersect", "union", "collapse", "setdiff")) {
  if (exists(f, envir = asNamespace("Biostrings"), inherits = FALSE))
    assign(f, get(f, envir = asNamespace("Biostrings")), envir = globalenv())
}

fa      <- snakemake@input[["fa"]]
out_rds <- snakemake@output[["rds"]]
out_dir <- dirname(out_rds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Genome: ", fa)

res <- extract.traits(in_file = fa, out_dir = out_dir)

# microtrait names the rds after the genome; normalise to the expected output
if (!file.exists(out_rds)) {
  cand <- if (!is.null(res$rds_file)) res$rds_file else
    list.files(out_dir, pattern = "\\.microtrait\\.rds$", full.names = TRUE)
  if (length(cand) != 1) stop("Could not find microtrait rds for ", fa)
  file.rename(cand, out_rds)
}
message("Done: ", out_rds)
