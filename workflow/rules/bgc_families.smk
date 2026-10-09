"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-08]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp bgc_families_all
Latest modification:
Purpose: To group antiSMASH regions of the dereplicated MQ MAGs into gene cluster families (GCFs) with
         BiG-SLiCE 2, and summarise GCFs per MAG and per sample (CoverM).
Notes (EI): BiG-SLiCE HMM models are downloaded once (software node) into <dbs_dir>/bigslice;
         the clustering runs as one cluster job on all regions.
"""

import os

BS = config["bgc_families"]
BS_DIR = os.path.join(ANNOT_RESULTS_DIR, "bgc", "bigslice")
BS_MODELS = os.path.join(ANNOT_DBS_DIR, "bigslice", "bigslice-models")
BS_REGIONS = os.path.join(ANNOT_RESULTS_DIR, "bgc", "antismash_regions_all.tsv")
BS_ANTISMASH = os.path.join(ANNOT_RESULTS_DIR, "bgc", "antismash", "per_mag")
BS_COVERM = os.path.join(ANNOT_RESULTS_DIR, "coverm",
                         "coverm_rel_abund_masked_cf{}.tsv".format(config.get("coverm_check", {}).get("min_covfrac", 0.1)))

localrules: bgc_families_all, download_bigslice_models, bigslice_input

###################
# RULES
###################
rule bgc_families_all:
    input:
        os.path.join(BS_DIR, "tables", "gcf_summary.tsv")
    output:
        touch("status/bgc_families.done")


# BiG-SLiCE HMM models (local: needs internet; same source/MD5 as download_bigslice_hmmdb)
rule download_bigslice_models:
    output:
        done=os.path.join(ANNOT_DBS_DIR, "bigslice", "bigslice-models.done")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/download_bigslice_models.log")
    params:
        url=BS["models_url"],
        md5=BS["models_md5"],
        models=BS_MODELS
    message:
        "Downloading BiG-SLiCE HMM models"
    shell:
        "(date && d=$(dirname {params.models}) && mkdir -p $d && cd $d && "
        "wget -c -O bigslice-models.tar.gz {params.url} && "
        "echo '{params.md5}  bigslice-models.tar.gz' | md5sum -c - && "
        "rm -rf {params.models} && mkdir -p {params.models} && tar -xzf bigslice-models.tar.gz -C {params.models} && "
        "rm bigslice-models.tar.gz && ls {params.models} && touch {output.done} && date) &> >(tee {log})"

# BiG-SLiCE input folder: one dataset, one sub-folder per MAG with its antiSMASH region GenBank files (symlinks)
rule bigslice_input:
    input:
        regions=BS_REGIONS
    output:
        datasets=os.path.join(BS_DIR, "input", "datasets.tsv")
    params:
        ds_dir=os.path.join(BS_DIR, "input", "skyline"),
        antismash=BS_ANTISMASH
    message:
        "Preparing BiG-SLiCE input from antiSMASH regions"
    run:
        import glob, shutil
        mags = set()
        with open(input.regions) as fh:
            next(fh)
            for line in fh:
                mags.add(line.split("\t", 1)[0])
        shutil.rmtree(params.ds_dir, ignore_errors=True)
        n = 0
        for m in sorted(mags):
            os.makedirs(os.path.join(params.ds_dir, m), exist_ok=True)
            for g in glob.glob(os.path.join(params.antismash, m, "*.region*.gbk")):
                os.symlink(g, os.path.join(params.ds_dir, m, os.path.basename(g))); n += 1
        with open(output.datasets, "w") as out:
            out.write("# dataset_name\tdataset_path\ttaxonomy_path\tdescription\n")
            out.write("skyline\tskyline\t\tSkyline MQ MAGs (antiSMASH regions)\n")
        print("BiG-SLiCE input: {} regions from {} MAGs".format(n, len(mags)))

# GCF clustering
rule bigslice_run:
    input:
        datasets=rules.bigslice_input.output.datasets,
        models=rules.download_bigslice_models.output.done
    output:
        db=os.path.join(BS_DIR, "output", "result", "data.db")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/bigslice_run.log")
    params:
        indir=os.path.join(BS_DIR, "input"),
        outdir=os.path.join(BS_DIR, "output"),
        models=BS_MODELS,
        threshold=BS["threshold"],
        extra=BS.get("extra", "")
    threads:
        BS["threads"]
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")
    message:
        "BiG-SLiCE GCF clustering (threshold {})".format(BS["threshold"])
    shell:
        "(date && pip show bigslice pyhmmer | grep -E '^(Name|Version)' && rm -rf {params.outdir} && "
        "bigslice -i {params.indir} {params.outdir} -t {threads} --threshold {params.threshold} "
        "--program_db_folder {params.models} {params.extra} && date) &> >(tee {log})"

# tables: region -> GCF, GCF summary, GCF x sample abundance (if the CoverM matrix exists)
rule bigslice_tables:
    input:
        db=rules.bigslice_run.output.db,
        regions=BS_REGIONS
    output:
        membership=os.path.join(BS_DIR, "tables", "gcf_membership.tsv"),
        summary=os.path.join(BS_DIR, "tables", "gcf_summary.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/bigslice_tables.log")
    params:
        script=srcdir("../bgc/bigslice_tables.py"),
        coverm=BS_COVERM
    conda:
        os.path.join(ENV_DIR, "bigslice.yaml")
    message:
        "Summarising BiG-SLiCE GCFs"
    shell:
        "(date && python {params.script} --db {input.db} --regions {input.regions} --coverm {params.coverm} "
        "--outdir $(dirname {output.summary}) && date) &> >(tee {log})"
