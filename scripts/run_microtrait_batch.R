#!/usr/bin/Rscript
# https://github.com/ukaraoz/microtrait -- one batch of genomes per job,
# genomes run in parallel; finished genomes (rds present) are skipped

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

mags   <- unlist(snakemake@params[["mags"]])
fas    <- unlist(snakemake@input[["fa"]])
outdir <- snakemake@params[["outdir"]]
ncores <- snakemake@threads
message("Batch of ", length(mags), " genomes on ", ncores, " cores")

run_one <- function(i) {
  mag <- mags[i]
  od  <- file.path(outdir, mag)
  rds <- file.path(od, paste0(mag, ".microtrait.rds"))
  if (file.exists(rds)) return(NA_character_)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    res <- extract.traits(in_file = fas[i], out_dir = od)
    # microtrait names the rds after the genome; normalise to the expected name
    if (!file.exists(rds)) {
      cand <- if (!is.null(res$rds_file)) res$rds_file else
        list.files(od, pattern = "\\.microtrait\\.rds$", full.names = TRUE)
      if (length(cand) != 1) stop("no rds produced")
      file.rename(cand, rds)
    }
    NA_character_
  }, error = function(e) paste0(mag, ": ", conditionMessage(e)))
}

res  <- parallel::mclapply(seq_along(mags), run_one, mc.cores = ncores, mc.preschedule = FALSE)
errs <- vapply(res, function(x) if (inherits(x, "try-error")) as.character(x) else
               if (is.na(x)) NA_character_ else x, character(1))
errs <- errs[!is.na(errs)]
if (length(errs) > 0) {
  message("FAILED (", length(errs), "):\n", paste(errs, collapse = "\n"))
  stop("microtrait failed for ", length(errs), " genome(s); rerun to retry only those")
}
message("Done: ", length(mags), " genomes")
file.create(snakemake@output[["done"]])
