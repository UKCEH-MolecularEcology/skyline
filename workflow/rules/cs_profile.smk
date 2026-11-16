"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-11-14]
Run: step "cs_profile", together with "extract", "cluster", "dereplicate" in the same run (config_16S_2026.yaml)
Purpose: Countryside Survey 2007 16S (Zenodo 10.5281/zenodo.12744255) linked to the dereplicated MAGs, and MAG
         occurrence across the survey. Uses this run's long-read 16S (cat_ass_16S), contig -> MAG -> dMAG table
         and GTDB-Tk summary, so the CS steps start as soon as those exist:
           1. long-read 16S trimmed to the CS region (V3-V4, CS primers)
           2. clustered with the CS OTUs (renamed ASV_CS_<OTU>)
           3. OTU -> long-read 16S -> contig -> MAG -> dMAG links (same scripts as dereplicated_mag_link)
           4. how many CS OTUs link to MAGs, MAG abundance across samples, habitats / soil (cs_mag_profile.py)
         Outputs in cs_profile.results_dir (results/16S_CS2007/).
"""

import os

CSP = config["cs_profile"]
CSP_DIR = os.path.abspath(CSP["results_dir"])
CSP_IDS = [str(i) for i in CSP.get("cluster_ids", [97])]
CSP_PREFIX = CSP.get("otu_prefix", "ASV_CS_")
CSP_OUT = ["otu_mag_links.tsv", "link_summary.tsv", "dmag_sample_abundance.tsv.gz",
           "dmag_habitat.tsv", "dmag_soil.tsv", "dmag_summary.tsv"]

# per-cluster tables can be large at low identity (90%): those steps are submitted jobs (resources in slurm.yaml)
localrules: cs_profile_all, cs_cat_otus

wildcard_constraints:
    csid="|".join(CSP_IDS)

rule cs_profile_all:
    input:
        expand(os.path.join(CSP_DIR, "CS_profile", "{csid}", "{f}"), csid=CSP_IDS, f=CSP_OUT),
        expand(os.path.join(CSP_DIR, "CS_profile", "{csid}", "dmag_16S_link_sources.tsv"), csid=CSP_IDS)
    output:
        touch("status/cs_profile.done")

# 1. long-read 16S (all 14 assemblies, untrimmed, from cat_ass_16S) -> CS region
rule cs_extract_16S:
    input:
        rules.cat_ass_16S.output.cat_fa
    output:
        ext_fa=os.path.join(CSP_DIR, "extracted/extracted_16S.fa"),
        map_fa=os.path.join(CSP_DIR, "extracted/extracted_map.txt")
    log:
        os.path.join(CSP_DIR, "logs/extract_16S.log")
    params:
        src=os.path.join(SRC_DIR, "modified_in_silico_pcr.pl"),
        fwd=CSP["fwd"],
        rev=CSP["rev"]
    message:
        "CS: trimming the long-read 16S to the CS amplicon region"
    shell:
        "(date && mkdir -p $(dirname {output.ext_fa}) && "
        "perl {params.src} -s {input} -a {params.fwd} -b {params.rev} -e -f {output.ext_fa} -n {output.map_fa} && "
        "echo \"long-read 16S in the CS region: $(grep -c '>' {output.ext_fa})\" && date) &> >(tee {log})"

# 2. CS OTUs (prefixed, so they are typed as amplicon) + trimmed long-read 16S
rule cs_cat_otus:
    input:
        otus=CSP["repseqs"],
        lr=rules.cs_extract_16S.output.ext_fa
    output:
        os.path.join(CSP_DIR, "concat_16S/cs_otus_lr_16S.fa")
    params:
        prefix=CSP_PREFIX
    shell:
        "awk -v p='{params.prefix}' '/^>/{{sub(/^>/, \">\" p)}} {{print}}' {input.otus} > {output} && cat {input.lr} >> {output}"

rule cs_cluster:
    input:
        rules.cs_cat_otus.output
    output:
        uc=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}.uc"),
        otus=os.path.join(CSP_DIR, "OTU/{csid}_cluster/otus_{csid}.fasta")
    log:
        os.path.join(CSP_DIR, "logs/cluster_{csid}.log")
    params:
        id=lambda wildcards: float(wildcards.csid) / 100,
        strand=CSP.get("cluster_strand", "both")    # long-read 16S are extracted in contig orientation
    envmodules:
        "vsearch/2.11.1"
    message:
        "CS: clustering CS OTUs with the long-read 16S at {wildcards.csid}%"
    shell:
        "(date && vsearch --cluster_fast {input} --id {params.id} --strand {params.strand} "
        "--centroids {output.otus} --uc {output.uc} && date) &> >(tee {log})"

rule cs_analyse_clusters:
    input:
        uc=rules.cs_cluster.output.uc
    output:
        mixed=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_mixed_clusters.csv"),
        amp_only=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_amplicon_only.csv"),
        lr_only=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_longread_only.csv"),
        sing_amp=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_singleton_amplicon.csv"),
        sing_lr=os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_singleton_longread.csv")
    log:
        os.path.join(CSP_DIR, "logs/analyse_{csid}.log")
    script:
        os.path.join(SRC_DIR, "analyse_uc_clusters.py")

rule cs_split_mixed:
    input:
        rules.cs_analyse_clusters.output.mixed
    output:
        os.path.join(CSP_DIR, "OTU/{csid}_cluster/clusters_{csid}_mixed_clusters_split.csv")
    log:
        os.path.join(CSP_DIR, "logs/split_mixed_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "split_mixed_cluster_members.py")
    message:
        "CS: splitting mixed clusters ({wildcards.csid}%)"
    shell:
        "(date && python {params.src} {input} {output} && date) &> >(tee {log})"

# 3. links, through this run's contig -> MAG -> dMAG table and GTDB-Tk summary
rule cs_merge_bins:
    input:
        bin_info=DMAG_CONTIGS,
        clusters=rules.cs_split_mixed.output
    output:
        os.path.join(CSP_DIR, "dereplicated_OTU/{csid}_cluster/mixed_clusters_{csid}_dereplicated_bin.csv")
    log:
        os.path.join(CSP_DIR, "logs/merge_bins_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "merge_dereplicated_bins_to_clusters.py")
    message:
        "CS: adding contig -> MAG -> dMAG to the clusters ({wildcards.csid}%)"
    shell:
        "(date && mkdir -p $(dirname {output}) && python {params.src} {input.bin_info} {input.clusters} {output} && date) "
        "&> >(tee {log})"

rule cs_map_otus:
    input:
        enriched=rules.cs_merge_bins.output,
        dmag_map=DMAG_CONTIGS
    output:
        os.path.join(CSP_DIR, "dereplicated_OTU/{csid}_cluster/otu_to_dmag_mapping_{csid}.tsv")
    log:
        os.path.join(CSP_DIR, "logs/map_otus_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "map_asvs_to_dmags.py")
    message:
        "CS: OTU -> dMAG links ({wildcards.csid}%)"
    shell:
        "(date && python {params.src} {input.enriched} {input.dmag_map} {output} && date) &> >(tee {log})"

rule cs_merge_tax:
    input:
        mapping=rules.cs_map_otus.output,
        gtdbtk=rules.concat_dereplicated_gtdbtk_summaries.output.combined
    output:
        os.path.join(CSP_DIR, "dereplicated_OTU/{csid}_cluster/otu_to_dmag_mapping_{csid}_with_tax.tsv")
    log:
        os.path.join(CSP_DIR, "logs/merge_tax_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "merge_asv_dmag_with_gtdbtk.py")
    message:
        "CS: adding GTDB-Tk taxonomy ({wildcards.csid}%)"
    shell:
        "(date && python {params.src} {input.mapping} {input.gtdbtk} {output} && date) &> >(tee {log})"

# 4. MAG occurrence across the survey
rule cs_mag_profile:
    input:
        links=rules.cs_merge_tax.output,
        otu_table=CSP["otu_table"],
        metadata=CSP["metadata"]
    output:
        [os.path.join(CSP_DIR, "CS_profile", "{csid}", f) for f in CSP_OUT]
    log:
        os.path.join(CSP_DIR, "logs/cs_mag_profile_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "cs_mag_profile.py"),
        taxonomy=CSP.get("taxonomy", ""),
        prefix=CSP_PREFIX,
        sample_col=CSP.get("sample_column", "PLOT"),
        habitat=" ".join(CSP.get("habitat_columns", [])),
        numeric=" ".join(CSP.get("numeric_columns", [])),
        out_dir=lambda wildcards: os.path.join(CSP_DIR, "CS_profile", wildcards.csid)
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "CS: MAG occurrence across the survey ({wildcards.csid}% clusters)"
    shell:
        "(date && python {params.src} --links {input.links} --otu-table {input.otu_table} "
        "--metadata {input.metadata} --taxonomy '{params.taxonomy}' --prefix '{params.prefix}' "
        "--sample-column {params.sample_col} --habitat-columns {params.habitat} "
        "--numeric-columns {params.numeric} --out-dir {params.out_dir} && date) &> >(tee {log})"

# per dMAG: which MAGs (winner and/or other dRep cluster members) carried the 16S evidence for its links
rule cs_dmag_link_sources:
    input:
        links=rules.cs_merge_tax.output,
        cdb=os.path.join(DREP_DIR, "data_tables", "Cdb.csv"),
        wdb=os.path.join(DREP_DIR, "data_tables", "Wdb.csv")
    output:
        os.path.join(CSP_DIR, "CS_profile", "{csid}", "dmag_16S_link_sources.tsv")
    log:
        os.path.join(CSP_DIR, "logs/dmag_16S_link_sources_{csid}.log")
    params:
        src=os.path.join(SRC_DIR, "dmag_link_sources.py")
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")      # python + pandas
    message:
        "CS: MAGs behind each dMAG's 16S links ({wildcards.csid}%)"
    shell:
        "(date && python {params.src} {input.links} {input.cdb} {input.wdb} {output} && date) &> >(tee {log})"

