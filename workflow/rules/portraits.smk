"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-07]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp portraits_all
Latest modification:
Purpose: To predict phenotypic traits of the dereplicated MQ MAGs with porTraits (grp-bork/porTraits;
         BacDive-AI, GenomeSPOT, Traitar, MICROPHERRET), run per batch of MAGs in a single SLURM job.
Notes (EI): porTraits code, Singularity images and the Pfam clan table are fetched by local rules
         (software node, internet); batches run on compute nodes without internet. GTDB-Tk is skipped
         (params.skip_gtdbtk, patched into main.nf).
"""

import os

PT = config["portraits"]
PT_TOOLS = os.path.join(config["annotations"]["work_dir"], "tools")
PT_REPO = os.path.join(PT_TOOLS, "porTraits_{}".format(PT["commit"][:7]))
PT_SIF_DIR = os.path.join(config["annotations"]["work_dir"], "singularity_cache", "portraits")
PT_SIFS = {k: os.path.join(PT_SIF_DIR, "{}.sif".format(k)) for k in PT["images"]}
PT_PFAM = os.path.join(ANNOT_DBS_DIR, "pfam/Pfam31.0/Pfam-A.clans.tsv.gz")
PT_DIR = os.path.join(ANNOT_RESULTS_DIR, "portraits")

localrules: portraits_all, install_portraits, pull_portraits_image, download_pfam_clans, collate_portraits

###################
# RULES
###################
rule portraits_all:
    input:
        os.path.join(PT_DIR, "portraits_results.tsv.gz")
    output:
        touch("status/portraits.done")


# porTraits code at a pinned commit + GTDB-Tk skip patch (local: needs internet)
rule install_portraits:
    output:
        main=os.path.join(PT_REPO, "main.nf")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/install_portraits.log")
    params:
        repo=PT["repo"],
        commit=PT["commit"],
        dest=PT_REPO,
        patch=srcdir("../portraits/patch_skip_gtdbtk.py")
    message:
        "Installing porTraits ({})".format(PT["commit"][:7])
    shell:
        "(date && rm -rf {params.dest} && mkdir -p $(dirname {params.dest}) && "
        "wget -q -O {params.dest}.tar.gz https://codeload.github.com/{params.repo}/tar.gz/{params.commit} && "
        "mkdir -p {params.dest} && tar -xzf {params.dest}.tar.gz -C {params.dest} --strip-components=1 && "
        "rm {params.dest}.tar.gz && python3 {params.patch} {output.main} && date) &> >(tee {log})"

# Singularity images as local SIFs (local: needs internet).
# One job per image: a failed pull never removes images that already succeeded.
# Pull into <sif>.tmp and rename only when complete; build space on node-local disk
# (/ei does not allow the lchown that unpacking needs); retry once.
rule pull_portraits_image:
    output:
        sif=os.path.join(PT_SIF_DIR, "{img}.sif")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/pull_image_{img}.log")
    params:
        uri=lambda wildcards: "docker://" + PT["images"][wildcards.img],
        build_tmp=lambda wildcards: os.path.join(PT.get("singularity_tmpdir", "/tmp"),
                                                 "{}_sif_build_{}".format(os.environ.get("USER", "user"), wildcards.img))
    wildcard_constraints:
        img="|".join(PT["images"])
    message:
        "Pulling porTraits image: {wildcards.img}"
    shell:
        "(date && mkdir -p $(dirname {output.sif}) && df -h $(dirname {params.build_tmp}) && "
        "for attempt in 1 2; do "
        "rm -rf {params.build_tmp} {output.sif}.tmp && mkdir -p {params.build_tmp} && "
        "if SINGULARITY_TMPDIR={params.build_tmp} APPTAINER_TMPDIR={params.build_tmp} "
        "singularity pull {output.sif}.tmp {params.uri}; then break; fi; "
        "echo \"pull attempt $attempt failed\"; done && "
        "rm -rf {params.build_tmp} && test -s {output.sif}.tmp && mv {output.sif}.tmp {output.sif} && "
        "date) &> >(tee {log})"

# Pfam 31.0 clan table (matches eggNOG 5.0.2; porTraits' emapper2matrix default)
rule download_pfam_clans:
    output:
        PT_PFAM
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/download_pfam_clans.log")
    params:
        url=PT["pfam_clans_url"]
    message:
        "Downloading Pfam 31.0 clan table"
    shell:
        "(date && mkdir -p $(dirname {output}) && wget -q -O {output} {params.url} && date) &> >(tee {log})"

# porTraits on one batch of MAGs (Nextflow local executor inside this job)
rule portraits_batch:
    input:
        fa=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]],
        main=rules.install_portraits.output.main,
        sifs=list(PT_SIFS.values()),
        pfam=PT_PFAM
    output:
        os.path.join(PT_DIR, "batches/{batch}/collated/portraits_results.tsv.gz")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/{batch}.log")
    params:
        rundir=lambda wildcards: os.path.join(PT_DIR, "batches", wildcards.batch),
        config=srcdir("../portraits/ei.config"),
        models=PT["models"],
        markers=PT["recognise_markers"],
        eggnog=PT["eggnog_db"],
        mem=PT["mem"],
        eggnog_cpus=PT["eggnog_cpus"],
        eggnog_mem=PT["eggnog_mem"],
        nxf_home=os.path.join(config["annotations"]["work_dir"], "nextflow_home"),
        sif_portraits=PT_SIFS["portraits"],
        sif_recognise=PT_SIFS["recognise"],
        sif_eggnog=PT_SIFS["eggnog"]
    threads:
        PT["threads"]
    conda:
        os.path.join(ENV_DIR, "nextflow.yaml")
    message:
        "Running porTraits on {wildcards.batch}"
    shell:
        "(date && mkdir -p {params.rundir}/input && cd {params.rundir} && "
        "for f in {input.fa}; do ln -sf $f input/; done && "
        "export NXF_OFFLINE=true NXF_HOME={params.nxf_home} NXF_ANSI_LOG=false "
        "PORTRAITS_CPUS={threads} PORTRAITS_MEM='{params.mem}' "
        "PORTRAITS_EGGNOG_CPUS={params.eggnog_cpus} PORTRAITS_EGGNOG_MEM='{params.eggnog_mem}' "
        "PORTRAITS_SIF_PORTRAITS={params.sif_portraits} PORTRAITS_SIF_RECOGNISE={params.sif_recognise} "
        "PORTRAITS_SIF_EGGNOG={params.sif_eggnog} && "
        "nextflow run {input.main} -c {params.config} -work-dir {params.rundir}/work -resume "
        "--input_dir {params.rundir}/input --output_dir {params.rundir} "
        "--metatraits_models {params.models} --recognise_marker_genes {params.markers} "
        "--eggnog_db {params.eggnog} --pfam_clade_map {input.pfam} && "
        "rm -rf {params.rundir}/work && date) &> >(tee {log})"

# all batches -> one table
rule collate_portraits:
    input:
        expand(os.path.join(PT_DIR, "batches/{batch}/collated/portraits_results.tsv.gz"), batch=sorted(ANNOT_BATCHES))
    output:
        os.path.join(PT_DIR, "portraits_results.tsv.gz")
    message:
        "Collating porTraits results for {} MAGs".format(len(ANNOT_MAGS))
    run:
        import gzip
        header = None
        with gzip.open(output[0], "wt") as out:
            for f in input:
                with gzip.open(f, "rt") as fh:
                    h = fh.readline()
                    if header is None:
                        header = h
                        out.write(h)
                    elif h != header:
                        raise ValueError("Header mismatch in " + f)
                    for line in fh:
                        out.write(line)
