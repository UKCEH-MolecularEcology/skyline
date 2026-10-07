"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-07]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 1 -rp coverm_check_all
Latest modification:
Purpose: To QC and collate the existing per-sample CoverM (genome mode) outputs against the
         dereplicated MQ MAG set (MQ_MAGs_all_0.95) and the Skyline sequencing sheet
"""

import os

# parameters
COVERM_DIR = config["coverm_check"]["coverm_dir"]
COVERM_SAMPLES = sorted(glob_wildcards(os.path.join(COVERM_DIR, "{sid,[^/]+}_output_coverm.tsv")).sid)

localrules: coverm_check_all

###################
# RULES
###################
rule coverm_check_all:
    input:
        os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_check_report.txt")
    output:
        touch("status/coverm_check.done")


# collating + QC of the per-sample CoverM tables
rule coverm_check:
    input:
        tsv=expand(os.path.join(COVERM_DIR, "{sid}_output_coverm.tsv"), sid=COVERM_SAMPLES),
        meta=config["coverm_check"]["metadata"]
    output:
        report=os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_check_report.txt"),
        long=os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_long.tsv.gz"),
        samples=os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_sample_summary.tsv"),
        genomes=os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_genome_summary.tsv"),
        ra=os.path.join(ANNOT_RESULTS_DIR, "coverm/coverm_rel_abund_matrix.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/coverm/coverm_check.log")
    params:
        src=os.path.join(SRC_DIR, "coverm_check.py"),
        mags_dir=ANNOT_MAGS_DIR,
        ext=ANNOT_MAGS_EXT,
        min_covfrac=config["coverm_check"]["min_covfrac"],
        max_unmapped=config["coverm_check"]["max_unmapped"]
    conda:
        os.path.join(ENV_DIR, "coverm_check.yaml")
    message:
        "Collating and checking CoverM outputs for {} samples".format(len(COVERM_SAMPLES))
    shell:
        "(date && "
        "python {params.src} --files {input.tsv} --metadata {input.meta} "
        "--genomes {params.mags_dir} --ext {params.ext} --out $(dirname {output.report}) "
        "--min-covfrac {params.min_covfrac} --max-unmapped {params.max_unmapped} && "
        "date) &> >(tee {log})"
