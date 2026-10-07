"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-07]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp microtrait_all
Latest modification: per-MAG jobs (restartable) + genome-set aggregation
Purpose: To run microtrait on the dereplicated MQ MAGs
"""

import os

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

# microtrait per MAG
rule microtrait_mag:
    input:
        fa=os.path.join(ANNOT_MAGS_DIR, "{mag}." + ANNOT_MAGS_EXT),
        installed="status/microtrait_installed.done"
    output:
        rds=os.path.join(ANNOT_RESULTS_DIR, "microtrait/per_genome/{mag}/{mag}.microtrait.rds")
    log:
        out=os.path.join(ANNOT_RESULTS_DIR, "logs/microtrait/{mag}.log")
    threads:
        config["microtrait"]["threads"]
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Running microtrait on {wildcards.mag}"
    script:
        os.path.join(SRC_DIR, "run_microtrait_genome.R")

# genome-set trait matrices
rule microtrait_genomeset:
    input:
        rds=expand(os.path.join(ANNOT_RESULTS_DIR, "microtrait/per_genome/{mag}/{mag}.microtrait.rds"), mag=ANNOT_MAGS)
    output:
        rds=os.path.join(ANNOT_RESULTS_DIR, "microtrait/genomeset_results.rds"),
        tables=directory(os.path.join(ANNOT_RESULTS_DIR, "microtrait/tables"))
    log:
        out=os.path.join(ANNOT_RESULTS_DIR, "logs/microtrait/genomeset.log")
    threads:
        config["microtrait"]["genomeset_threads"]
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Creating microtrait genome-set results for {} MAGs".format(len(ANNOT_MAGS))
    script:
        os.path.join(SRC_DIR, "microtrait_genomeset.R")
