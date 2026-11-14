"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/bin_info.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To merge bin info with vsearch cluster outputs following LR-ASV clustering with amplicon ASVs
"""

import os

# parameters
CLUSTER_IDS = config["cluster_ids"]

# contig -> MAG -> dMAG table and GTDB-Tk summaries; defaults = 2025 layout of drep_dir
DMAG = config.get("dmag_link") or {}
DMAG_BUILD = bool(DMAG.get("build_contig_table", False))
DMAG_CONTIGS = (os.path.join(RESULTS_DIR, "dereplicated_OTU", "contig_MAGs_dMAGs.tsv") if DMAG_BUILD
                else DMAG.get("contig_table", os.path.join(DREP_DIR, "results/contig_MAGs_dMAGs.tsv")))
DMAG_GTDB_ARC = DMAG.get("gtdb_arc", os.path.join(DREP_DIR, "gtdb/gtdbtk.ar53.summary.tsv"))
DMAG_GTDB_BAC = DMAG.get("gtdb_bac", os.path.join(DREP_DIR, "gtdb/gtdbtk.bac120.summary.tsv"))

localrules: merge_dereplicated_bin_into_clusters, map_asvs_to_dreplicated_mags, concat_dereplicated_gtdbtk_summaries, merge_asv_to_dmag_with_gtdbtk

###################
# RULES
###################
rule dereplicated_mag_link_all:
    input:
        expand(os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/mixed_clusters_{id}_dereplicated_bin.csv"), id=CLUSTER_IDS),
        expand(os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}.tsv"), id=CLUSTER_IDS),
        os.path.join(RESULTS_DIR, "dereplicated_OTU/dmag_gtdbtk_summary.tsv"),
        expand(os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}_with_tax.tsv"), id=CLUSTER_IDS)
    output:
        touch("status/dereplicated_mag_link.done")


# merging the dereplicated bin and contig info with the cluster info
rule merge_dereplicated_bin_into_clusters:
    input:
        bin_info=DMAG_CONTIGS,
        clusters=lambda wildcards: os.path.join(RESULTS_DIR, f"OTU/{wildcards.id}_cluster/clusters_{wildcards.id}_mixed_clusters_split.csv")
    output:
        merged=os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/mixed_clusters_{id}_dereplicated_bin.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/OTU/{id}_cluster/merge_with_bin_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "merge_dereplicated_bins_to_clusters.py")
    message:
        "Merging dereplicated bin info into clusters for cluster ID {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.bin_info} {input.clusters} {output.merged} && date) &> >(tee {log}) 
        """

# Linking the ASVs to dereplicated LR-MAGs
rule map_asvs_to_dreplicated_mags:
    input:
        enriched=rules.merge_dereplicated_bin_into_clusters.output.merged,
        dmag_map=DMAG_CONTIGS
    output:
        mapping=os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}.log")
    params:
        src=os.path.join(SRC_DIR, "map_asvs_to_dmags.py")
    message:
        "Mapping ASVs to dereplicated MAGs for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.src} {input.enriched} {input.dmag_map} {output.mapping} && date) &> >(tee {log})
        """
#rule map_asvs_to_dreplicated_mags:
#    input:
#        enriched=rules.merge_dereplicated_bin_into_clusters.output.merged,
#    output:
#        mapping=os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}.tsv")
#    log:
#        os.path.join(RESULTS_DIR, "logs/dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}.log")
#    params:
#        src=os.path.join(SRC_DIR, "map_asvs_to_dmags.py")
#    message:
#        "Mapping ASVs to MAGs & sample names for cluster {wildcards.id}"
#    shell:
#        """
#        (date && python {params.src} {input.enriched} {output.mapping} && date) &> >(tee {log})
#        """

# Concatenating the GTDBtk summaries for the dereplicated bins
rule concat_dereplicated_gtdbtk_summaries:
    input:
        arc=DMAG_GTDB_ARC,
        bac=DMAG_GTDB_BAC
    output:
        combined=os.path.join(RESULTS_DIR, "dereplicated_OTU/dmag_gtdbtk_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/dereplicated_OTU/dmag_concat_gtdbtk.log")
    params:
        src=os.path.join(SRC_DIR, "concat_gtdbtk_summaries.py")
    message:
        "Concatenating GTDB-Tk summaries for all dmags"
    shell:
        """
        (date && awk 'NR==1 || FNR>1' {input.arc} {input.bac} > {output.combined} && date) &> >(tee {log})
        """

# Merging ASV to dereplicated MAG mapping with taxonomy
rule merge_asv_to_dmag_with_gtdbtk:
    input:
        mapping=rules.map_asvs_to_dreplicated_mags.output.mapping,
        gtdbtk=rules.concat_dereplicated_gtdbtk_summaries.output.combined
    output:
        merged=os.path.join(RESULTS_DIR, "dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}_with_tax.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/dereplicated_OTU/{id}_cluster/asv_to_dmag_mapping_{id}_with_tax.log")
    params:
        script=os.path.join(SRC_DIR, "merge_asv_dmag_with_gtdbtk.py")
    message:
        "Merging ASV-MAG mappings with GTDB-Tk taxonomy for cluster {wildcards.id}"
    shell:
        """
        (date && python {params.script} {input.mapping} {input.gtdbtk} {output.merged} && date) &> >(tee {log})
        """


# contig -> MAG -> dMAG table for a dRep collection (Cdb/Wdb + MAG FASTA headers), checked against the assemblies
if DMAG_BUILD:
    rule build_dmag_contig_table:
        input:
            cdb=os.path.join(DREP_DIR, "data_tables", "Cdb.csv"),
            wdb=os.path.join(DREP_DIR, "data_tables", "Wdb.csv")
        output:
            DMAG_CONTIGS
        log:
            os.path.join(RESULTS_DIR, "logs/dereplicated_OTU/build_contig_table.log")
        params:
            src=os.path.join(SRC_DIR, "build_dmag_contig_table.py"),
            genome_dirs=" ".join(DMAG.get("genome_dirs", [os.path.join(DREP_DIR, "dereplicated_genomes")])),
            prefix_map=" ".join("{}={}".format(k, v) for k, v in (DMAG.get("mag_prefix_to_assembly") or {}).items()),
            # templates passed through literally (lambdas stop Snakemake treating {sample}/{asm} as wildcards)
            default_asm=lambda wildcards: DMAG.get("default_assembly", "hifi_individual_{sample}"),
            asm_contigs=lambda wildcards: DMAG.get("assembly_contigs", "{ass_dir}/{asm}/contigs/contigs.fa").replace("{ass_dir}", ASS_DIR)
        message:
            "Building the contig -> MAG -> dMAG table from the dRep tables"
        shell:
            "(date && python {params.src} --cdb {input.cdb} --wdb {input.wdb} --genome-dirs {params.genome_dirs} "
            "--prefix-map {params.prefix_map} --default-assembly '{params.default_asm}' "
            "--assembly-contigs '{params.asm_contigs}' --out {output} && date) &> >(tee {log})"

