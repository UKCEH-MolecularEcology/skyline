"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/bin_info.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To merge bin info with vsearch cluster outputs following LR-ASV clustering with amplicon ASVs
"""

import os

localrules: add_sample_to_clustering_circ, add_sample_to_clustering_consensus, add_sample_to_gtdbtk_bac_summary, add_sample_to_gtdbtk_arc_summary, bin_info_all, combine_all_bin_info

###################
# RULES
###################
rule mag_link_all:
    input:
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_mag_mapping_{id}.tsv"), id=CLUSTER_IDS),
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_summary.tsv"), sample=SAMPLES),
        os.path.join(RESULTS_DIR, "bin_info/all_samples_gtdbtk_summary.tsv"),
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_mag_mapping_{id}_with_tax.tsv"), id=CLUSTER_IDS)
    output:
        touch("status/mag_link.done")

# Linking the ASVs to LR-MAGs
rule map_asvs_to_mags:
    input:
        enriched=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters_with_bin.csv")
    output:
        mapping=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_mag_mapping_{id}.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/asv_to_mag_mapping_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "map_asvs_to_mags.py")
    message:
        "Mapping ASVs to MAGs & sample names for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.enriched} {output.mapping} && date) &> >(tee {log})
        """

# Concatenating the GTDBtk summaries for the bins
rule concat_gtdbtk_summaries:
    input:
        arc=lambda wildcards: os.path.join(RESULTS_DIR, "bin_info", wildcards.sample, f"{wildcards.sample}_gtdbtk_arc_summary.tsv"),
        bac=lambda wildcards: os.path.join(RESULTS_DIR, "bin_info", wildcards.sample, f"{wildcards.sample}_gtdbtk_bac_summary.tsv")
    output:
        combined=os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/{sample}_concat_gtdbtk.log")
    params:
        src=os.path.join(SRC_DIR, "concat_gtdbtk_summaries.py")
    message:
        "Concatenating GTDB-Tk summaries for sample {wildcards.sample}"
    shell:
        """
        (date && python {params.src} results/bin_info/{wildcards.sample} {output.combined} && date) &> >(tee {log})
        """

rule concat_all_samples_gtdbtk_summaries:
    input:
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_summary.tsv"), sample=SAMPLES)
    output:
        combined=os.path.join(RESULTS_DIR, "bin_info/all_samples_gtdbtk_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/concat_all_samples_gtdbtk.log")
    message:
        "Concatenating all per-sample GTDB-Tk summaries into one file"
    shell:
        """
        # Extract header from the first input file
        header=$(head -n1 {input[0]})
        echo "$header" > {output.combined}

        # Loop through all per-sample summary files and append content skipping header
        for file in {input}; do
            tail -n +2 "$file"
        done >> {output.combined}
        """

# Merging ASV-MAG mapping with taxonomy
rule merge_asv_to_mag_with_gtdbtk:
    input:
        mapping=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_mag_mapping_{id}.tsv"),
        gtdbtk=rules.concat_all_samples_gtdbtk_summaries.output.combined
    output:
        merged=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_mag_mapping_{id}_with_tax.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/asv_to_mag_mapping_{id}_with_tax.log")
    params:
        script=os.path.join(SRC_DIR, "merge_asv_to_mag_with_gtdbtk.py")
    message:
        "Merging ASV-MAG mappings with GTDB-Tk taxonomy for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.script} {input.mapping} {input.gtdbtk} {output.merged} && date) &> >(tee {log})
        """
