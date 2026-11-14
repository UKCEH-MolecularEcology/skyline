"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-11-14]
Run: step "cs_profile" (after "extract", "cluster", "dereplicate"), e.g. with config/config_16S_CS2007.yaml
Purpose: Profile dereplicated MAGs across an amplicon survey (Countryside Survey 2007 16S) through the
         16S OTU -> long-read 16S -> contig -> MAG links: how many survey OTUs link to MAGs, and in which
         habitats / environmental zones the MAGs occur (scripts/cs_mag_profile.py)
"""

import os

CSP = config["cs_profile"]
CSP_IDS = [str(i) for i in CSP.get("cluster_ids", [97])]
CSP_DIR = os.path.join(RESULTS_DIR, "CS_profile")
CSP_OUT = ["otu_mag_links.tsv", "link_summary.tsv", "dmag_sample_abundance.tsv.gz",
           "dmag_habitat.tsv", "dmag_soil.tsv", "dmag_summary.tsv"]

localrules: cs_profile_all

rule cs_profile_all:
    input:
        expand(os.path.join(CSP_DIR, "{id}", "{f}"), id=CSP_IDS, f=CSP_OUT)
    output:
        touch("status/cs_profile.done")

rule cs_mag_profile:
    input:
        links=os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}_with_tax.tsv"),
        otu_table=CSP["otu_table"],
        metadata=CSP["metadata"]
    output:
        [os.path.join(CSP_DIR, "{id}", f) for f in CSP_OUT]
    log:
        os.path.join(RESULTS_DIR, "logs/CS_profile/cs_mag_profile_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "cs_mag_profile.py"),
        taxonomy=CSP.get("taxonomy", ""),
        prefix=config.get("asv_prefix", ""),
        sample_col=CSP.get("sample_column", "PLOT"),
        habitat=" ".join(CSP.get("habitat_columns", [])),
        numeric=" ".join(CSP.get("numeric_columns", [])),
        out_dir=lambda wildcards: os.path.join(CSP_DIR, wildcards.id)
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "Profiling MAGs across the survey via 16S links ({wildcards.id}% clusters)"
    shell:
        "(date && python {params.src} --links {input.links} --otu-table {input.otu_table} "
        "--metadata {input.metadata} --taxonomy '{params.taxonomy}' --prefix '{params.prefix}' "
        "--sample-column {params.sample_col} --habitat-columns {params.habitat} "
        "--numeric-columns {params.numeric} --out-dir {params.out_dir} && date) &> >(tee {log})"
