"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-09]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp community_analysis_all
Purpose: Community-level analysis of MAG annotations (microTrait, porTraits, RGI, BacMet, antiSMASH, BiG-SLiCE GCFs):
         rhizosphere vs bulk soil, biochar vs control, across seasonal (fortnightly), daily and diel (4-hourly)
         sampling (scripts/community_analysis.R, R in the existing microTrait conda env).
Outputs: <results>/analysis_tables/<set>/..., <results>/figures/<set>/...,
         <results>/figures/community_analysis_report.pdf (all figures, tables, notes and the R code)
"""

import os

CA_TABLES = os.path.join(ANNOT_RESULTS_DIR, "analysis_tables")
CA_REPORT = os.path.join(ANNOT_RESULTS_DIR, "figures", "community_analysis_report.pdf")
# wait for every annotation step enabled in config["steps"]
CA_AFTER = {"coverage_check": "status/coverm_check.done", "traits": "status/microtrait.done",
            "amr": "status/amr.done", "portraits": "status/portraits.done",
            "bgc": "status/bgc.done", "bgc_families": "status/bgc_families.done"}

localrules: community_analysis_all

rule community_analysis_all:
    input:
        CA_REPORT
    output:
        touch("status/community_analysis.done")

rule community_analysis:
    input:
        steps=[v for k, v in CA_AFTER.items() if k in config["steps"]],
        metadata=config["coverm_check"]["metadata"]
    output:
        summary=os.path.join(CA_TABLES, "summary_terms.tsv"),
        report=CA_REPORT
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "community_analysis.log")
    params:
        script=os.path.join(SRC_DIR, "community_analysis.R"),
        sets=",".join(config.get("community_analysis", {}).get("sets",
                      ["microtrait", "portraits", "amr", "bacmet", "bgc", "gcf"])),
        nperm=config.get("community_analysis", {}).get("nperm", 999),
        min_mags=config.get("community_analysis", {}).get("min_mags", 20),
        granularity=config.get("community_analysis", {}).get("granularity", 3)
    conda:
        os.path.join(ENV_DIR, "microtrait.yaml")
    message:
        "Community-level analysis (compartment x treatment x time) of: {params.sets}"
    shell:
        "(date && Rscript {params.script} --results-dir {ANNOT_RESULTS_DIR} --metadata {input.metadata} "
        "--sets {params.sets} --nperm {params.nperm} --min-mags {params.min_mags} --granularity {params.granularity} "
        "&& date) &> >(tee {log})"
