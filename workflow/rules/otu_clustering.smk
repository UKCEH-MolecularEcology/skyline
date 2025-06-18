"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/otu_clustering.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To clusters ASVs into OTUs
"""

import os

localrules: cat_asvs, otu_all, parse_cluster_logs, summarize_cluster_counts, combine_cluster_summaries

CLUSTER_IDS=config["cluster_ids"]

###################
# RULES
###################
rule otu_all:
    input:
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otus_{id}.fasta"), id=CLUSTER_IDS),
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/otu_table_{id}.txt"), id=CLUSTER_IDS),
        os.path.join(RESULTS_DIR, "OTU/otu_cluster_log_summary.tsv"),
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_amplicon_only.csv"), id=CLUSTER_IDS),
        os.path.join(RESULTS_DIR, "OTU/all_cluster_summaries.tsv")
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

# Cluster log summary extraction
rule parse_cluster_logs:
    input:
        logs=expand(os.path.join(RESULTS_DIR, "logs/cluster/cluster_{id}.log"), id=CLUSTER_IDS)
    output:
        summary=os.path.join(RESULTS_DIR, "OTU/otu_cluster_log_summary.tsv")
    message:
        "Parsing OTU clustering logs to summarize cluster stats"
    run:
        import re
        import pandas as pd

        summary_data = []
        pattern = re.compile(
            r"Clusters: (?P<clusters>\d+).*?"
            r"Size min (?P<min>\d+), max (?P<max>\d+), avg (?P<avg>[0-9.]+)\n"
            r"Singletons: (?P<singletons>\d+), (?P<singletons_seq_pct>[0-9.]+)% of seqs, (?P<singletons_cluster_pct>[0-9.]+)% of clusters"
        )

        for log_file in input.logs:
            with open(log_file) as f:
                content = f.read()
                match = pattern.search(content)
                if match:
                    id_match = re.search(r"cluster_(\d+)\.log", log_file)
                    cluster_pct = id_match.group(1)
                    stats = match.groupdict()
                    stats["cluster_pct"] = cluster_pct
                    stats["total_seqs"] = 107605  # Hardcoded from your input
                    summary_data.append(stats)

        df = pd.DataFrame(summary_data)
        df = df[["cluster_pct", "total_seqs", "clusters", "min", "max", "avg", 
                 "singletons", "singletons_seq_pct", "singletons_cluster_pct"]]
        df.columns = ["Cluster %", "Total Seqs", "Clusters", "Min Size", "Max Size", "Avg Size", 
                      "Singletons", "% Seqs as Singletons", "% Clusters as Singletons"]

        df.sort_values(by="Cluster %", inplace=True)
        df.to_csv(output.summary, sep="\t", index=False)

# Analyze OTU clusters
rule analyse_otu_clusters:
    input:
        uc=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}.uc")
    output:
        mixed=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_mixed_clusters.csv"),
        amp_only=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_amplicon_only.csv"),
        lr_only=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_longread_only.csv"),
        sing_amp=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_singleton_amplicon.csv"),
        sing_lr=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_singleton_longread.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/cluster/analyse_{id}.log")
    message:
        "Analyse cluster composition for identity {wildcards.id}"
    script:
        os.path.join(SRC_DIR, "analyse_uc_clusters.py")

#rule summarize_cluster_counts:
#    input:
#        amplicon_only=rules.analyse_otu_clusters.output.amp_only,
#        longread_only=rules.analyse_otu_clusters.output.lr_only,
#        mixed=rules.analyse_otu_clusters.output.mixed,
#        singleton_amplicon=rules.analyse_otu_clusters.output.sing_amp,
#        singleton_longread=rules.analyse_otu_clusters.output.sing_lr
#    output:
#        summary=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/cluster_{id}_summary.tsv")
#    log:
#        os.path.join(RESULTS_DIR, "logs/cluster/cluster_{id}_summary.log")
#    message:
#        "Summarizing cluster counts for identity {wildcards.id}"
#    shell:
#        """
#        (
#            echo -e "Cluster ID\tAmplicon Only\tLongread Only\tMixed\tSingleton Amplicon\tSingleton Longread" &&
#            echo -e "{wildcards.id}\t$(wc -l < {input.amplicon_only})\t$(wc -l < {input.longread_only})\t$(wc -l < {input.mixed})\t$(wc -l < {input.singleton_amplicon})\t$(wc -l < {input.singleton_longread})
#        ) > {output.summary}
#        """ + " &> >(tee {log})"

rule summarize_cluster_counts:
    input:
        lambda wildcards: expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/clusters_{id}_{type}.csv"),
                                 id=wildcards.id,
                                 type=["amplicon_only", "longread_only", "mixed_clusters", "singleton_amplicon", "singleton_longread"])
    output:
        summary=os.path.join(RESULTS_DIR, "OTU/{id}_cluster/cluster_{id}_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/cluster/cluster_{id}_summary.log")
    message:
        "Summarizing cluster counts for identity {wildcards.id}"
    shell:
        r"""
        (
            echo -e "Cluster ID\tAmplicon Only\tLongread Only\tMixed\tSingleton Amplicon\tSingleton Longread" &&
            echo -e "{wildcards.id}\t$(wc -l < results/OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_amplicon_only.csv)\t\
$(wc -l < results/OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_longread_only.csv)\t\
$(wc -l < results/OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_mixed_clusters.csv)\t\
$(wc -l < results/OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_singleton_amplicon.csv)\t\
$(wc -l < results/OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_singleton_longread.csv)"
        ) > {output.summary}
        """ + " &> >(tee {log})"

rule combine_cluster_summaries:
    input:
        expand(os.path.join(RESULTS_DIR, "OTU/{id}_cluster/cluster_{id}_summary.tsv"), id=CLUSTER_IDS)
    output:
        combined=os.path.join(RESULTS_DIR, "OTU/all_cluster_summaries.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/cluster/combined_cluster_summary.log")
    message:
        "Combining all cluster summaries into one file"
    shell:
        """
        head -n 1 {input[0]} > {output.combined}
        tail -n +2 -q {input} >> {output.combined}
        """ + " &> >(tee {log})"

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
        


