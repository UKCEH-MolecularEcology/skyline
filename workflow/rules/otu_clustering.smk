"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/kraken2.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on reads
"""

import os

localrules: cat_asvs, otu_all

CLUSTER_IDS=config["cluster_ids"]

###################
# RULES
###################
rule otu_all:
    input:
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otus_{id}.fasta"), id=CLUSTER_IDS),
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otu_table_{id}.txt"), id=CLUSTER_IDS)
    output:
        touch("status/otu_clustering.done")        

# Concatenating amplicon and LR ASVs
rule cat_asvs:
    input:
        amplicon=config["ASV"],
        lr=os.path.join(RESULTS_DIR, "extracted/extracted_16S.fa")
    output:
        all_cat=os.path.join(RESULTS_DIR, "concat_16S/amplicon_lr_ASVs.fa")
    log:
        os.path.join(RESULTS_DIR, "logs/concat/amplicon_lr_concat.log")
    message:
        "Concatenating amplicon and extracted-LR 16S ASVs"
    shell:
        "(date && cat {input} > {output.all_cat} && date) &> >(tee {log})"

# OTU clustering
rule cluster_otus:
    input:
        rules.cat_asvs.output.all_cat
    output:
        otus=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otus_{id}.fasta"),
        uc=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}.uc")
    log:
        os.path.join(RESULTS_DIR, "logs/cluster/cluster_{id}.log")        
    message:
        "Cluster OTUs at percentage: {wildcards.id}"
    params:
        id=lambda wildcards: float(wildcards.id)/100
    envmodules:
        "vsearch"
#    conda:
#        os.path.join(ENV_DIR, "vsearch.yaml")
    shell:
        "(date && mkdir -p $(dirname {output.otus}) && "
        "vsearch --cluster_fast {input} --id {params.id} --centroids {output.otus} --uc {output.uc} && "
        "date) &> >(tee {log})"

# OTU table
rule make_otu_table:
    input:
        uc=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}.uc"),
        fasta=rules.cat_asvs.output.all_cat,
        centroids=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otus_{id}.fasta")
    output:
        table=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otu_table_{id}.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/cluster/otu_table_{id}.log")
    message:
        "Generating OTU table at {wildcards.id}% identity"
#    conda:
#        os.path.join(ENV_DIR, "vsearch.yaml")
    params:
        id=lambda wildcards: float(wildcards.id)/100
    envmodules:
        "vsearch"
    shell:
        "(date && "
        "vsearch --usearch_global {input.fasta} --db {input.centroids} --id {params.id} --otutabout {output.table} && "
        "date) &> >(tee {log})"
        


