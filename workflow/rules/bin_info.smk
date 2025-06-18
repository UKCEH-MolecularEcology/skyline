"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/bin_info.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To merge bin info with vsearch cluster outputs following LR-ASV clustering with amplicon ASVs
"""

import os

localrules: add_sample_to_clustering_circ, add_sample_to_clustering_consensus, add_sample_to_gtdbtk_bac_summary, add_sample_to_gtdbtk_arc_summary, bin_info_all

###################
# RULES
###################
rule bin_info_all:
    input:
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_clustering_circ.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_clustering_consensus.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_bac_summary.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_arc_summary.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters_split.csv"), id=CLUSTER_IDS)
    output:
        touch("status/bin_info.done")

# Rules to add sample info to existing files
rule add_sample_to_clustering_circ:
    input:
        circ=os.path.join(ASS_DIR, "{sample}/metamdbg/binning/circ/clustering_circ.csv")
    output:
        circ_out=os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_clustering_circ.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/{sample}_clustering_circ.log")
    message:
        "Adding sample name to clustering_circ for {wildcards.sample}"
    shell:
        """
        (date && awk -v sample="{wildcards.sample}" 'BEGIN{{FS=OFS=","}} NR==1{{$3="sample"}} NR>1{{$3=sample}} {{print $1, $2, $3}}' {input.circ} > {output.circ_out} && date) &> >(tee {log})
        """

rule add_sample_to_clustering_consensus:
    input:
        cons=os.path.join(ASS_DIR, "{sample}/metamdbg/binning/consensus_LR/clustering_consensus_LR.csv")
    output:
        cons_out=os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_clustering_consensus.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/{sample}_clustering_consensus.log")
    message:
        "Adding sample name to clustering_consensus_LR for {wildcards.sample}"
    shell:
        """
        (date && awk -v sample="{wildcards.sample}" 'BEGIN{{FS=OFS=","}} NR==1{{$3="sample"}} NR>1{{$3=sample}} {{print $1, $2, $3}}' {input.cons} > {output.cons_out} && date) &> >(tee {log})
        """

rule add_sample_to_gtdbtk_bac_summary:
    input:
        gtdb=os.path.join(ASS_DIR, "{sample}/metamdbg/MAGs/gtdb_LR/classify/gtdbtk.bac120.summary.tsv")
    output:
        gtdb_out=os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_bac_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/{sample}_gtdbtk_bac_summary.log")
    message:
        "Adding sample name to GTDB bacterial summary for {wildcards.sample}"
    shell:
        """
        (date && awk -v sample="{wildcards.sample}" 'BEGIN{{FS=OFS="\\t"}} NR==1{{$17="sample"}} NR>1{{$17=sample}} {{print}}' {input.gtdb} > {output.gtdb_out} && date) &> >(tee {log})
        """

rule add_sample_to_gtdbtk_arc_summary:
    input:
        gtdb=os.path.join(ASS_DIR, "{sample}/metamdbg/MAGs/gtdb_LR/classify/gtdbtk.ar53.summary.tsv")
    output:
        gtdb_out=os.path.join(RESULTS_DIR, "bin_info/{sample}/{sample}_gtdbtk_arc_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bin_info/{sample}_gtdbtk_arc_summary.log")
    message:
        "Adding sample name to GTDB archaeal summary for {wildcards.sample}"
    shell:
        """
        (date && awk -v sample="{wildcards.sample}" 'BEGIN{{FS=OFS="\\t"}} NR==1{{$17="sample"}} NR>1{{$17=sample}} {{print}}' {input.gtdb} > {output.gtdb_out} && date) &> >(tee {log})
        """

rule split_mixed_cluster_members:
    input:
        mixed=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters.csv")
    output:
        split=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters_split.csv")
    log:
        os.path.join(RESULTS_DIR, "OTU/{id}_cluster/split_mixed_clusters_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "split_mixed_cluster_members.py")
    message:
        "Splitting members into individual lines for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.mixed} {output.split} && date) &> >(tee {log})
        """
