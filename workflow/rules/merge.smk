"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-09]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 1 -rp merge_all
Purpose: To merge all MAG annotations (CoverM, microTrait, porTraits, RGI, BacMet, antiSMASH, GECCO, BiG-SLiCE)
         into one table, one row per MAG (scripts/merge_annotations.py)
"""

import os

MERGE_OUT = os.path.join(ANNOT_RESULTS_DIR, "merged", "skyline_MAG_annotations.tsv.gz")
# wait for every annotation step that is enabled in config["steps"]
MERGE_AFTER = {"coverage_check": "status/coverm_check.done", "traits": "status/microtrait.done",
               "amr": "status/amr.done", "portraits": "status/portraits.done",
               "bgc": "status/bgc.done", "bgc_families": "status/bgc_families.done"}
MERGE_INPUTS = [v for k, v in MERGE_AFTER.items() if k in config["steps"]]

localrules: merge_all

rule merge_all:
    input:
        MERGE_OUT
    output:
        touch("status/merge.done")

rule merge_annotations:
    input:
        MERGE_INPUTS
    output:
        table=MERGE_OUT,
        columns=MERGE_OUT.replace(".tsv.gz", ".columns.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/merge_annotations.log")
    params:
        script=os.path.join(SRC_DIR, "merge_annotations.py"),
        granularity=config.get("merge", {}).get("mt_granularity", 3)
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "Merging all MAG annotations into one table"
    shell:
        "(date && python {params.script} --results-dir {ANNOT_RESULTS_DIR} --mags-dir {ANNOT_MAGS_DIR} "
        "--mags-ext {ANNOT_MAGS_EXT} --mt-granularity {params.granularity} --out {output.table} && date) &> >(tee {log})"
