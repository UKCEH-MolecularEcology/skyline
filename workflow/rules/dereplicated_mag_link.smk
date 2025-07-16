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
rule dereplicated_mag_link_all:
    input:
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}_with_tax.tsv"), id=CLUSTER_IDS)
    output:
        touch("status/dereplicated_mag_link.done")

# gathering the bin info for all dereplicated mags
rule dereplicated_all_bin_info:
    input:
        os.path.join(DREP_DIR, "results/contig_MAGs_dMAGs.tsv")
    output:
        merged=os.path.join(RESULTS_DIR, "bin_info/all_sample_bin_info.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/combine_all_bin_info.log")
    message:
        "Combining bin info from all samples into one file"
    shell:
        "(date && cp -v {input[0]} {output.merged} && date) &> >(tee {log})"

rule merge_dereplicated_bin_into_clusters:
    input:
        bin_info=rules.dereplicated_all_bin_info.output.merged,
        clusters=lambda wildcards: os.path.join(RESULTS_DIR, f"OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_mixed_clusters_split.csv")
    output:
        merged=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters_with_bin.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/merge_with_bin_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "merge_bin_into_clusters_updated.py")
    message:
        "Merging bin info into clusters for cluster ID {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.bin_info} {input.clusters} {output.merged} && date) &> >(tee {log}) 
        """
        # check the params.src to make sure column names match

# Linking the ASVs to dereplicated LR-MAGs
rule map_asvs_to_dreplicated_mags:
    input:
        enriched=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters_with_bin.csv")
    output:
        mapping=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "map_asvs_to_mags.py")
    message:
        "Mapping ASVs to MAGs & sample names for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.enriched} {output.mapping} && date) &> >(tee {log})
        """

# Concatenating the GTDBtk summaries for the dereplicated bins
rule concat_dereplicated_gtdbtk_summaries:
    input:
        arc=os.path.join(DREP_DIR, "gtdb/gtdbtk.ar53.summary.tsv"),
        bac=os.path.join(DREP_DIR, "gtdb/gtdbtk.bac120.summary.tsv")
    output:
        combined=os.path.join(RESULTS_DIR, "bin_info/all_sample_gtdbtk_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/all_sample_dereplicated_concat_gtdbtk.log")
    params:
        src=os.path.join(SRC_DIR, "concat_gtdbtk_summaries.py")
    message:
        "Concatenating GTDB-Tk summaries for all mags"
    shell:
        """
        (date && python {params.src} results/bin_info/ {output.combined} && date) &> >(tee {log})
        """        

# Merging ASV to dereplicated MAG mapping with taxonomy
rule merge_asv_to_dreplicated_mag_with_gtdbtk:
    input:
        mapping=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}.tsv"),
        gtdbtk=rules.concat_dereplicated_gtdbtk_summaries.output.combined
    output:
        merged=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}_with_tax.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/asv_to_dereplicated_mag_mapping_{id}_with_tax.log")
    params:
        script=os.path.join(SRC_DIR, "merge_asv_to_mag_with_gtdbtk.py")
    message:
        "Merging ASV-MAG mappings with GTDB-Tk taxonomy for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.script} {input.mapping} {input.gtdbtk} {output.merged} && date) &> >(tee {log})
        """
