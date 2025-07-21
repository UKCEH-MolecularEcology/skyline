"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/compare_globdb.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To use fastANI and compare against SMAG catalogue and GLOB_DB mags
"""

import os

# parameters
CLUSTER_IDS = config["cluster_ids"]

localrules: 

###################
# RULES
###################
rule compare_globdb_all:
    input:
        expand(os.path.join(RESULTS_DIR, "fastani/skyline_vs_{catalogue}_ani.txt"), catalogue=["smag", "globdb"])
    output:
        touch("status/compare_globdb.done")


# creating the query and reference list for skyline mags and SMAG/GLOB_DB mags
rule create_lists:
    input:
        skyline=os.path.join(DREP_DIR, "drep/dereplicated_genomes/"),
        globdb=os.path.join(DB_DIR, "globdb/globdb_r226_genome_fasta/")
    output:
        sky_list=os.path.join(RESULTS_DIR, "fastani/skyline_mags.txt"),
        glob_list=os.path.join(RESULTS_DIR, "fastani/glob_mags.txt"),
        smag_list=os.path.join(RESULTS_DIR, "fastani/smag_mags.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/fastani/query_ref_lists.log")
    message:
        "Creating query and reference lists for downstream fastANI processing"
    shell:
        "(date && "
        "find {input.skyline}/ -name "*.fa" > {output.sky_list} && "
        "find {input.globdb}/ -name "*.fa" > {output.glob_list} && "
        "find {input.globdb}/ -name "SMAGOTU_*.fa" > {output.smag_list} && "
        "date) &> >(tee {log})"

# Comparing Skyline and SMAG catalogue
rule compare_skyline_smag:
    input:
        query=rules.create_lists.output.smag_list,
        ref=rules.create_lists.output.smag_list
    output:
        sky_smag=os.path.join(RESULTS_DIR, "fastani/skyline_vs_smag_ani.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/fastani/compare_skyline_smag.log")
    conda:
        "anvio-8"
    threads:
        config["fastani"]["threads"]
    message:
        "Comparing Skyline MAGs to SMAG catalogue"
    shell:
        "(date && "
        "fastANI --ql {input.query} --rl {input.ref} -o {output.sky_smag} -t {threads} && "
        "date) &> >(tee {log})"

# Comparing Skyline and GLOB_DB catalogue
rule compare_skyline_globdb:
    input:
        query=rules.create_lists.output.smag_list,
        ref=rules.create_lists.output.glob_list
    output:
        sky_glob=os.path.join(RESULTS_DIR, "fastani/skyline_vs_globdb_ani.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/fastani/compare_skyline_globdb.log")
    conda:
        "anvio-8"
    threads:
        config["fastani"]["threads"]
    message:
        "Comparing Skyline MAGs to GLOB_DB catalogue"
    shell:
        "(date && "
        "fastANI --ql {input.query} --rl {input.ref} -o {output.sky_glob} -t {threads} && "
        "date) &> >(tee {log})"


