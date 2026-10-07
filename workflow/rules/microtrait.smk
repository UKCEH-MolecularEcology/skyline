"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-07]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp microtrait_all
Latest modification: batched (annotations.batch_size MAGs per job, per-MAG rds, finished MAGs skipped)
Purpose: To run microtrait on the dereplicated MQ MAGs
"""

import os

MT_DIR = os.path.join(ANNOT_RESULTS_DIR, "microtrait/per_genome")

localrules: microtrait_all, install_microtrait

###################
# RULES
###################
rule microtrait_all:
    input:
        os.path.join(ANNOT_RESULTS_DIR, "microtrait/genomeset_results.rds")
    output:
        touch("status/microtrait.done")


# installing microtrait from GitHub + preparing HMM models (local rule: needs internet)
rule install_microtrait:
    output:
        "status/microtrait_installed.done"
    log:
        out=os.path.join(ANNOT_RESULTS_DIR, "logs/microtrait/install_microtrait.log")
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Installing microtrait"
    script:
        os.path.join(SRC_DIR, "install_microtrait.R")

# microtrait on a batch of MAGs -> per_genome/{mag}/{mag}.microtrait.rds
rule microtrait_batch:
    input:
        fa=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]],
        installed="status/microtrait_installed.done"
    output:
        done=os.path.join(ANNOT_RESULTS_DIR, "microtrait/batches/{batch}.done")
    log:
        out=os.path.join(ANNOT_RESULTS_DIR, "logs/microtrait/{batch}.log")
    params:
        mags=lambda wildcards: ANNOT_BATCHES[wildcards.batch],
        outdir=MT_DIR
    threads:
        config["microtrait"]["threads"]
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Running microtrait on {wildcards.batch}"
    script:
        os.path.join(SRC_DIR, "run_microtrait_batch.R")

# genome-set trait matrices
rule microtrait_genomeset:
    input:
        expand(os.path.join(ANNOT_RESULTS_DIR, "microtrait/batches/{batch}.done"), batch=sorted(ANNOT_BATCHES))
    output:
        rds=os.path.join(ANNOT_RESULTS_DIR, "microtrait/genomeset_results.rds"),
        tables=directory(os.path.join(ANNOT_RESULTS_DIR, "microtrait/tables"))
    log:
        out=os.path.join(ANNOT_RESULTS_DIR, "logs/microtrait/genomeset.log")
    params:
        rds=[os.path.join(MT_DIR, m, m + ".microtrait.rds") for m in ANNOT_MAGS]
    threads:
        config["microtrait"]["genomeset_threads"]
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Creating microtrait genome-set results for {} MAGs".format(len(ANNOT_MAGS))
    script:
        os.path.join(SRC_DIR, "microtrait_genomeset.R")
