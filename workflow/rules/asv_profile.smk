"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-11-15]
Run: step "asv_profile", with "extract", "cluster", "dereplicate" in the same run (config_16S_2026.yaml)
Purpose: MAG occurrence across the Skyline amplicon samples, from this run's ASV -> dMAG links and the ASV count
         table (same outputs as the Countryside Survey profile; scripts/cs_mag_profile.py).
         Outputs in <results_dir>/ASV_profile/<id>/.
"""

import os

ASP = config["asv_profile"]
ASP_IDS = [str(i) for i in ASP.get("cluster_ids", CLUSTER_IDS)]
ASP_OUT = ["otu_mag_links.tsv", "link_summary.tsv", "dmag_sample_abundance.tsv.gz",
           "dmag_habitat.tsv", "dmag_soil.tsv", "dmag_summary.tsv"]

localrules: asv_profile_all

rule asv_profile_all:
    input:
        expand(os.path.join(RESULTS_DIR, "ASV_profile", "{aid}", "{f}"), aid=ASP_IDS, f=ASP_OUT),
        expand(os.path.join(RESULTS_DIR, "ASV_profile", "{aid}", "dmag_16S_link_sources.tsv"), aid=ASP_IDS)
    output:
        touch("status/asv_profile.done")

rule asv_mag_profile:
    input:
        links=os.path.join(RESULTS_DIR, "dereplicated_OTU/{aid}_cluster/asv_to_dmag_mapping_{aid}_with_tax.tsv"),
        counts=ASP["feature_table"],
        metadata=ASP["sample_table"]
    output:
        [os.path.join(RESULTS_DIR, "ASV_profile", "{aid}", f) for f in ASP_OUT]
    log:
        os.path.join(RESULTS_DIR, "logs/ASV_profile/asv_mag_profile_{aid}.log")
    wildcard_constraints:
        aid="|".join(ASP_IDS)
    params:
        src=os.path.join(SRC_DIR, "cs_mag_profile.py"),
        taxonomy=ASP.get("taxonomy", ""),
        sample_col=ASP.get("sample_column", "ID"),
        rows="--features-in-rows" if ASP.get("features_in_rows", True) else "",
        replace=ASP.get("sample_name_replace", ""),
        habitat=" ".join(ASP.get("group_columns", [])),
        numeric=" ".join(ASP.get("numeric_columns", [])),
        out_dir=lambda wildcards: os.path.join(RESULTS_DIR, "ASV_profile", wildcards.aid)
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "Skyline: MAG occurrence across the amplicon samples ({wildcards.aid}% clusters)"
    shell:
        "(date && python {params.src} --links {input.links} --otu-table {input.counts} {params.rows} "
        "--sample-name-replace '{params.replace}' --metadata {input.metadata} --taxonomy '{params.taxonomy}' "
        "--sample-column {params.sample_col} --habitat-columns {params.habitat} "
        "--numeric-columns {params.numeric} --out-dir {params.out_dir} && date) &> >(tee {log})"

# per dMAG: which MAGs (winner and/or other dRep cluster members) carried the 16S evidence for its links
rule asv_dmag_link_sources:
    input:
        links=os.path.join(RESULTS_DIR, "dereplicated_OTU/{aid}_cluster/asv_to_dmag_mapping_{aid}_with_tax.tsv"),
        cdb=os.path.join(DREP_DIR, "data_tables", "Cdb.csv"),
        wdb=os.path.join(DREP_DIR, "data_tables", "Wdb.csv")
    output:
        os.path.join(RESULTS_DIR, "ASV_profile", "{aid}", "dmag_16S_link_sources.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/ASV_profile/dmag_16S_link_sources_{aid}.log")
    wildcard_constraints:
        aid="|".join(ASP_IDS)
    params:
        src=os.path.join(SRC_DIR, "dmag_link_sources.py")
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "Skyline: MAGs behind each dMAG's 16S links ({wildcards.aid}%)"
    shell:
        "(date && python {params.src} {input.links} {input.cdb} {input.wdb} {output} && date) &> >(tee {log})"

