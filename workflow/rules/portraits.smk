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
PT_PFAM = os.path.join(ANNOT_DBS_DIR, "pfam", PT.get("pfam_release", "Pfam31.0"), "Pfam-A.clans.tsv.gz")
PT_DIR = os.path.join(ANNOT_RESULTS_DIR, "portraits")
# eggNOG-mapper (porTraits' own container: emapper 2.1.12 + eggNOG 5.0.2) run once per batch, split per MAG
PT_EMAP_DIR = os.path.join(PT_DIR, "eggnog_mapper")
PT_EGGNOG_MODE = PT.get("eggnog_mode", "batch")      # "batch" (fast) or "per_genome" (porTraits default)
PT_FAA = config["annotations"].get("prodigal_existing", "")

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

# Pfam clan table: maps eggNOG-mapper Pfam names -> stable PF accessions (the models' features)
rule download_pfam_clans:
    output:
        PT_PFAM
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/download_pfam_clans.log")
    params:
        url=PT["pfam_clans_url"]
    message:
        "Downloading Pfam clan table ({})".format(PT.get("pfam_release", "Pfam31.0"))
    shell:
        "(date && mkdir -p $(dirname {output}) && wget -q -O {output} {params.url} && date) &> >(tee {log})"

# eggNOG-mapper on all proteins of a batch at once (same container, version and database porTraits uses),
# DIAMOND database copied to node-local disk; results split into <PT_EMAP_DIR>/<mag>/<mag>.emapper.annotations
rule eggnog_batch:
    input:
        faa=lambda wildcards: [PT_FAA.format(mag=m) for m in ANNOT_BATCHES[wildcards.batch]],
        sif=PT_SIFS["eggnog"]
    output:
        done=os.path.join(PT_EMAP_DIR, "batches", "{batch}.done")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/portraits/eggnog_{batch}.log")
    params:
        mags=lambda wildcards: " ".join(ANNOT_BATCHES[wildcards.batch]),
        db=PT["eggnog_db"],
        outdir=PT_EMAP_DIR,
        rundir=lambda wildcards: os.path.join(PT_EMAP_DIR, "runs", wildcards.batch),
        local=lambda wildcards: os.path.join(PT.get("local_tmpdir", "/tmp"), "{}_eggnog_{}".format(os.environ.get("USER", "user"), wildcards.batch)),
        split=srcdir("../portraits/split_emapper.py"),
        block=PT.get("eggnog_block_size", 8),
        chunks=PT.get("eggnog_index_chunks", 1)
    threads:
        PT.get("eggnog_batch_threads", 32)
    message:
        "eggNOG-mapper (2.1.12 / eggNOG 5.0.2) on {wildcards.batch}"
    shell:
        "(date && rm -rf {params.local} && mkdir -p {params.rundir} {params.local} && "
        "mags=({params.mags}); faas=({input.faa}); "
        "for i in ${{!mags[@]}}; do "
        "awk -v m=${{mags[$i]}} '/^>/{{sub(/^>/, \">\" m \"|\")}} {{print}}' ${{faas[$i]}}; "
        "done > {params.local}/proteins.faa && grep -c '>' {params.local}/proteins.faa && "
        "cp {params.db}/eggnog_proteins.dmnd {params.local}/ && "
        "singularity exec {input.sif} emapper.py -i {params.local}/proteins.faa --data_dir {params.db} "
        "--dmnd_db {params.local}/eggnog_proteins.dmnd -m diamond --dmnd_algo 0 --cpu {threads} "
        "--block_size {params.block} --index_chunks {params.chunks} "
        "--temp_dir {params.local} --output_dir {params.rundir} --output {wildcards.batch} --override && "
        "singularity exec {input.sif} python3 {params.split} --annotations {params.rundir}/{wildcards.batch}.emapper.annotations "
        "--mags {params.mags} --outdir {params.outdir} && "
        "rm -rf {params.local} && touch {output.done} && date) &> >(tee {log})"

# porTraits on one batch of MAGs (Nextflow local executor inside this job)
rule portraits_batch:
    input:
        fa=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]],
        main=rules.install_portraits.output.main,
        sifs=list(PT_SIFS.values()),
        pfam=PT_PFAM,
        # eggNOG-mapper output of this batch (eggnog_mode "batch"); porTraits then skips its per-genome eggNOG step
        emapper=lambda wildcards: [os.path.join(PT_EMAP_DIR, "batches", wildcards.batch + ".done")] if PT_EGGNOG_MODE == "batch" else []
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
        sif_eggnog=PT_SIFS["eggnog"],
        max_skip=PT.get("max_skipped_frac", 0.05),
        emapper_dir=PT_EMAP_DIR if PT_EGGNOG_MODE == "batch" else ""
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
        "PORTRAITS_SIF_EGGNOG={params.sif_eggnog} PORTRAITS_EMAPPER_DIR='{params.emapper_dir}' && "
        "nextflow run {input.main} -c {params.config} -work-dir {params.rundir}/work -resume "
        "--input_dir {params.rundir}/input --output_dir {params.rundir} "
        "--metatraits_models {params.models} --recognise_marker_genes {params.markers} "
        "--eggnog_db {params.eggnog} --pfam_clade_map {input.pfam} && "
        "(grep 'Error is ignored' .nextflow.log > skipped_tasks.txt || true) && "
        "n=$(ls input | wc -l) && s=$(wc -l < skipped_tasks.txt) && "
        "echo \"skipped tasks: $s (genomes in batch: $n)\" && "
        # fail (and keep the work dir for debugging) if skipped tasks exceed max_skipped_frac of genomes
        "if [ $(awk -v s=$s -v n=$n 'BEGIN{{print (s > {params.max_skip} * n) ? 1 : 0}}') -eq 1 ]; then "
        "echo \"ERROR: too many skipped tasks ($s for $n genomes); work dir kept\"; "
        "rm -f {output}; exit 1; fi && "
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
