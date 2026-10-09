"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-08]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp bgc_all
Latest modification:
Purpose: To predict biosynthetic gene clusters (BGCs) in the dereplicated MQ MAGs with antiSMASH and GECCO,
         per batch of MAGs, collated into one table per tool.
Notes (EI): antiSMASH runs from the official antismash/standalone image (databases included), pulled once
         on the software node as a local SIF; GECCO from a conda env. Batch jobs run on compute nodes.
"""

import os

BGC = config["bgc"]
BGC_DIR = os.path.join(ANNOT_RESULTS_DIR, "bgc")
BGC_SIF = os.path.join(config["annotations"]["work_dir"], "singularity_cache", "bgc", "antismash.sif")

localrules: bgc_all, pull_antismash_image, collate_antismash, collate_gecco, bgc_overlap

###################
# RULES
###################
rule bgc_all:
    input:
        os.path.join(BGC_DIR, "antismash_regions_all.tsv"),
        os.path.join(BGC_DIR, "gecco_clusters_all.tsv"),
        os.path.join(BGC_DIR, "bgc_overlap.tsv")
    output:
        touch("status/bgc.done")


################### antiSMASH
# official image incl. databases (local: needs internet); pull into <sif>.tmp, rename when complete
rule pull_antismash_image:
    output:
        sif=BGC_SIF
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/pull_antismash_image.log")
    params:
        uri="docker://" + BGC["antismash"]["image"],
        build_tmp=os.path.join(BGC.get("singularity_tmpdir", "/tmp"), "{}_sif_build_antismash".format(os.environ.get("USER", "user")))
    message:
        "Pulling antiSMASH image ({})".format(BGC["antismash"]["image"])
    shell:
        "(date && mkdir -p $(dirname {output.sif}) && df -h $(dirname {params.build_tmp}) && "
        "for attempt in 1 2; do "
        "rm -rf {params.build_tmp} {output.sif}.tmp && mkdir -p {params.build_tmp} && "
        "if SINGULARITY_TMPDIR={params.build_tmp} APPTAINER_TMPDIR={params.build_tmp} "
        "singularity pull {output.sif}.tmp {params.uri}; then break; fi; "
        "echo \"pull attempt $attempt failed\"; done && "
        "rm -rf {params.build_tmp} && test -s {output.sif}.tmp && mv {output.sif}.tmp {output.sif} && "
        "date) &> >(tee {log})"

# antiSMASH on a batch of MAGs, several MAGs in parallel; finished MAGs are skipped on reruns
rule antismash_batch:
    input:
        fa=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]],
        sif=BGC_SIF
    output:
        done=os.path.join(BGC_DIR, "antismash/batches/{batch}.done")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/antismash_{batch}.log")
    params:
        outdir=os.path.join(BGC_DIR, "antismash/per_mag"),
        pairs=lambda wildcards: " ".join("{} {}".format(m, os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT))
                                          for m in ANNOT_BATCHES[wildcards.batch]),
        mags=lambda wildcards: " ".join(ANNOT_BATCHES[wildcards.batch]),
        parallel=BGC["antismash"]["parallel"],
        cpus=BGC["antismash"]["cpus_per_mag"],
        extra=BGC["antismash"].get("extra", ""),
        runner=srcdir("../bgc/run_antismash_one.sh"),
        regions=srcdir("../bgc/antismash_regions.py"),
        max_fail=BGC.get("max_failed_frac", 0.05)
    threads:
        BGC["antismash"]["parallel"] * BGC["antismash"]["cpus_per_mag"]
    message:
        "antiSMASH on {wildcards.batch}"
    shell:
        "(date && mkdir -p {params.outdir} && "
        "echo {params.pairs} | xargs -n 2 -P {params.parallel} sh -c "
        "'bash {params.runner} \"$0\" \"$1\" {params.outdir} {input.sif} {params.cpus} {params.regions} {params.extra}' && "
        "n=0; f=0; for m in {params.mags}; do n=$((n+1)); grep -q '^ok' {params.outdir}/$m/STATUS || f=$((f+1)); done; "
        "echo \"antiSMASH failed for $f of $n MAGs\"; "
        "if [ $(awk -v f=$f -v n=$n 'BEGIN{{print (f > {params.max_fail} * n) ? 1 : 0}}') -eq 1 ]; then exit 1; fi && "
        "touch {output.done} && date) &> >(tee {log})"

# all MAGs -> one table of regions
rule collate_antismash:
    input:
        expand(os.path.join(BGC_DIR, "antismash/batches/{batch}.done"), batch=sorted(ANNOT_BATCHES))
    output:
        os.path.join(BGC_DIR, "antismash_regions_all.tsv")
    params:
        outdir=os.path.join(BGC_DIR, "antismash/per_mag")
    message:
        "Collating antiSMASH regions for {} MAGs".format(len(ANNOT_MAGS))
    run:
        header, failed = None, []
        with open(output[0], "w") as out:
            for m in ANNOT_MAGS:
                f = os.path.join(params.outdir, m, m + ".regions.tsv")
                if not os.path.exists(f):
                    failed.append(m); continue
                with open(f) as fh:
                    h = fh.readline()
                    if header is None:
                        header = h; out.write(h)
                    for line in fh:
                        out.write(line)
        with open(output[0].replace(".tsv", ".failed_mags.txt"), "w") as fh:
            fh.write("\n".join(failed) + ("\n" if failed else ""))
        print("antiSMASH: {} MAGs without results (see *.failed_mags.txt)".format(len(failed)))


################### GECCO
rule gecco_batch:
    input:
        fa=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]]
    output:
        done=os.path.join(BGC_DIR, "gecco/batches/{batch}.done")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/gecco_{batch}.log")
    params:
        outdir=os.path.join(BGC_DIR, "gecco/per_mag"),
        pairs=lambda wildcards: " ".join("{} {}".format(m, os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT))
                                          for m in ANNOT_BATCHES[wildcards.batch]),
        parallel=BGC["gecco"]["parallel"],
        cpus=BGC["gecco"]["cpus_per_mag"],
        extra=BGC["gecco"].get("extra", ""),
        mags=lambda wildcards: " ".join(ANNOT_BATCHES[wildcards.batch]),
        max_fail=BGC.get("max_failed_frac", 0.05)
    threads:
        BGC["gecco"]["parallel"] * BGC["gecco"]["cpus_per_mag"]
    conda:
        os.path.join(ENV_DIR, "gecco.yaml")
    message:
        "GECCO on {wildcards.batch}"
    shell:
        # one output folder per MAG; a MAG counts as done once its .clusters.tsv exists (GECCO writes it last)
        "(date && gecco --version && mkdir -p {params.outdir} && "
        "echo {params.pairs} | xargs -n 2 -P {params.parallel} sh -c '"
        "d={params.outdir}/$0; if ls $d/*.clusters.tsv >/dev/null 2>&1; then echo \"skip $0\"; exit 0; fi; "
        "rm -rf $d $d.FAILED && gecco run --genome $1 --output-dir $d --jobs {params.cpus} {params.extra} > $d.log 2>&1 "
        "&& rm -f $d.log && echo \"done $0\" || {{ mv $d.log $d.FAILED; echo \"FAILED $0\"; }}' && "
        "n=0; f=0; for m in {params.mags}; do n=$((n+1)); [ -e {params.outdir}/$m.FAILED ] && f=$((f+1)); done; "
        "echo \"GECCO failed for $f of $n MAGs\"; "
        "if [ $(awk -v f=$f -v n=$n 'BEGIN{{print (f > {params.max_fail} * n) ? 1 : 0}}') -eq 1 ]; then exit 1; fi && "
        "touch {output.done} && date) &> >(tee {log})"

# all MAGs -> one table of clusters (Genome column added)
rule collate_gecco:
    input:
        expand(os.path.join(BGC_DIR, "gecco/batches/{batch}.done"), batch=sorted(ANNOT_BATCHES))
    output:
        os.path.join(BGC_DIR, "gecco_clusters_all.tsv")
    params:
        outdir=os.path.join(BGC_DIR, "gecco/per_mag")
    message:
        "Collating GECCO clusters for {} MAGs".format(len(ANNOT_MAGS))
    run:
        import glob
        header = None
        with open(output[0], "w") as out:
            for m in ANNOT_MAGS:
                for f in sorted(glob.glob(os.path.join(params.outdir, m, "*.clusters.tsv"))):
                    with open(f) as fh:
                        h = fh.readline()
                        if header is None:
                            header = h; out.write("Genome\t" + h)
                        for line in fh:
                            out.write(m + "\t" + line)


################### antiSMASH vs GECCO
# regions/clusters matched by MAG + contig + coordinate overlap
rule bgc_overlap:
    input:
        antismash=os.path.join(BGC_DIR, "antismash_regions_all.tsv"),
        gecco=os.path.join(BGC_DIR, "gecco_clusters_all.tsv")
    output:
        pairs=os.path.join(BGC_DIR, "bgc_overlap.tsv"),
        summary=os.path.join(BGC_DIR, "bgc_overlap_summary.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/bgc/bgc_overlap.log")
    params:
        script=srcdir("../bgc/bgc_overlap.py"),
        min_bp=BGC.get("overlap_min_bp", 1)
    message:
        "Matching antiSMASH regions and GECCO clusters"
    shell:
        "(date && python3 {params.script} --antismash {input.antismash} --gecco {input.gecco} "
        "--out {output.pairs} --summary {output.summary} --min_overlap_bp {params.min_bp} && "
        "cat {output.summary} && date) &> >(tee {log})"
