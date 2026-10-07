"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-07]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4 -rp amr_all
Latest modification: batched (annotations.batch_size MAGs per job; finished MAGs skipped)
Purpose: To identify antimicrobial (RGI/CARD) and metal/biocide (BacMet2 EXP) resistance genes
         in the dereplicated MQ MAGs
"""

import os

# reuse existing prodigal proteins where available, otherwise call genes here
PRODIGAL_EXISTING = config["annotations"].get("prodigal_existing", "")

def faa_for(wildcards):
    if PRODIGAL_EXISTING:
        f = PRODIGAL_EXISTING.format(mag=wildcards.mag)
        if os.path.exists(f):
            return f
    return os.path.join(ANNOT_RESULTS_DIR, "prodigal/{}.faa".format(wildcards.mag))

def faa_for_mag(mag):
    class _W: pass
    w = _W(); w.mag = mag
    return faa_for(w)

RGI_DIR = os.path.join(ANNOT_RESULTS_DIR, "amr/rgi")

# NOTE: setup_rgi_db is NOT local -- it runs in the RGI container, which can only see /ei on compute nodes
localrules: amr_all, download_rgi_db, collate_rgi, download_bacmet, collate_bacmet

###################
# RULES
###################
rule amr_all:
    input:
        os.path.join(ANNOT_RESULTS_DIR, "amr/rgi_all.tsv"),
        os.path.join(ANNOT_RESULTS_DIR, "amr/bacmet_all.tsv")
    output:
        touch("status/amr.done")


################### RGI
# downloading CARD into the shared database folder (skipped if card.json exists)
rule download_rgi_db:
    output:
        json=os.path.join(ANNOT_DBS_DIR, "rgi/card.json")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/download_rgi_db.log")
    params:
        url=config["rgi"]["db_url"]
    message:
        "Downloading CARD data for RGI"
    shell:
        "(date && cd $(dirname {output.json}) && "
        "wget -O card-data.tar.bz2 {params.url} --no-check-certificate && "
        "tar -xjf card-data.tar.bz2 && rm card-data.tar.bz2 && "
        "date) &> >(tee {log})"

# loading CARD into <dbs>/rgi/localDB
# NOTE: to make sure that the same DB is used for all targets
rule setup_rgi_db:
    input:
        rules.download_rgi_db.output.json
    output:
        "status/rgi_setup.done"
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/setup_rgi_db.log")
    container:
        config["rgi"]["container"]
    message:
        "Loading the CARD database for RGI"
    shell:
        "(date && "
        "(cd $(dirname {input}) && "
        " (test -d localDB || rgi load --card_json card.json --local) && "
        " rgi database --version --local > card_version.txt) && "
        "cat $(dirname {input})/card_version.txt && "
        "touch {output} && date) &> >(tee {log})"

# running RGI on a batch of MAGs (contigs) -> amr/rgi/{mag}.txt
rule rgi_batch:
    input:
        fna=lambda wildcards: [os.path.join(ANNOT_MAGS_DIR, m + "." + ANNOT_MAGS_EXT) for m in ANNOT_BATCHES[wildcards.batch]],
        setup="status/rgi_setup.done"
    output:
        done=os.path.join(ANNOT_RESULTS_DIR, "amr/rgi/batches/{batch}.done")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/rgi/{batch}.log")
    threads:
        config["rgi"]["threads"]
    params:
        db=os.path.join(ANNOT_DBS_DIR, "rgi"),
        outdir=RGI_DIR,
        aligner=config["rgi"]["alignment_tool"],
        extra=config["rgi"]["extra"]
    container:
        config["rgi"]["container"]
    message:
        "Running RGI on {wildcards.batch}"
    shell:
        # --local reads ./localDB, hence running from the DB folder; MAGs with output are skipped
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run once
        "(date && mkdir -p {params.outdir} && cd {params.db} && failed=0 && "
        "for fna in {input.fna}; do "
        "mag=$(basename $fna .{ANNOT_MAGS_EXT}); out={params.outdir}/$mag; "
        "if [ -s $out.txt ]; then echo \"skip $mag\"; continue; fi; "
        "rgi main --input_sequence $fna --output_file $out --input_type contig --local -a {params.aligner} --clean {params.extra} -n {threads} || "
        "rgi main --input_sequence $fna --output_file $out --input_type contig --local -a {params.aligner} --clean {params.extra} -n {threads} || true; "
        "if [ ! -s $out.txt ]; then echo \"FAILED $mag\"; failed=$((failed+1)); fi; "
        "done && "
        "if [ $failed -gt 0 ]; then echo \"$failed MAG(s) failed; rerun to retry only those\"; exit 1; fi && "
        "touch {output.done} && date) &> >(tee {log})"

# collating RGI tables
rule collate_rgi:
    input:
        expand(os.path.join(ANNOT_RESULTS_DIR, "amr/rgi/batches/{batch}.done"), batch=sorted(ANNOT_BATCHES))
    output:
        os.path.join(ANNOT_RESULTS_DIR, "amr/rgi_all.tsv")
    params:
        txt=[os.path.join(RGI_DIR, m + ".txt") for m in ANNOT_MAGS]
    message:
        "Collating RGI results for {} MAGs".format(len(ANNOT_MAGS))
    run:
        header_written = False
        with open(output[0], "w") as out:
            for f in params.txt:
                mag = os.path.basename(f)[:-len(".txt")]
                with open(f) as fh:
                    header = fh.readline()
                    if not header_written:
                        out.write("Genome\t" + header)
                        header_written = True
                    for line in fh:
                        out.write(mag + "\t" + line)


################### BACMET
# gene calling for MAGs without existing prodigal proteins
rule prodigal_mag:
    input:
        os.path.join(ANNOT_MAGS_DIR, "{mag}." + ANNOT_MAGS_EXT)
    output:
        faa=os.path.join(ANNOT_RESULTS_DIR, "prodigal/{mag}.faa"),
        gff=os.path.join(ANNOT_RESULTS_DIR, "prodigal/{mag}.gff")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/prodigal/{mag}.log")
    conda:
        os.path.join(ENV_DIR, "prodigal.yaml")
    message:
        "Running prodigal on {wildcards.mag}"
    shell:
        "(date && prodigal -i {input} -a {output.faa} -f gff -o {output.gff} -p single -q && date) &> >(tee {log})"

# downloading BacMet2 experimentally confirmed DB (skipped if present)
rule download_bacmet:
    output:
        fasta=os.path.join(ANNOT_DBS_DIR, "bacmet/BacMet2_EXP_database.fasta")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/download_bacmet.log")
    params:
        url=config["bacmet"]["db_url"],
        map_url=config["bacmet"]["map_url"]
    message:
        "Downloading BacMet2 EXP"
    shell:
        "(date && cd $(dirname {output.fasta}) && "
        "wget -O $(basename {output.fasta}) {params.url} && "
        "(wget -N {params.map_url} || echo 'WARNING: BacMet mapping file not downloaded') && "
        "date) &> >(tee {log})"

rule bacmet_makedb:
    input:
        rules.download_bacmet.output.fasta
    output:
        os.path.join(ANNOT_DBS_DIR, "bacmet/bacmet_exp.dmnd")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/bacmet_makedb.log")
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    message:
        "Building DIAMOND database for BacMet2 EXP"
    shell:
        "(date && diamond makedb --in {input} -d $(dirname {output})/bacmet_exp && date) &> >(tee {log})"

# DIAMOND blastp: proteins of a batch of MAGs (headers prefixed "<mag>|") vs BacMet2 EXP
rule bacmet_batch:
    input:
        faa=lambda wildcards: [faa_for_mag(m) for m in ANNOT_BATCHES[wildcards.batch]],
        db=rules.bacmet_makedb.output
    output:
        tsv=os.path.join(ANNOT_RESULTS_DIR, "amr/bacmet/batches/{batch}.bacmet.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs/amr/bacmet/{batch}.log")
    threads:
        config["bacmet"]["threads"]
    params:
        mags=lambda wildcards: " ".join(ANNOT_BATCHES[wildcards.batch]),
        pid=config["bacmet"]["id"],
        qcov=config["bacmet"]["qcov"],
        evalue=config["bacmet"]["evalue"]
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    message:
        "Running BacMet (DIAMOND) on {wildcards.batch}"
    shell:
        "(date && "
        "mags=({params.mags}); faas=({input.faa}); "
        "for i in ${{!mags[@]}}; do "
        "awk -v m=${{mags[$i]}} '/^>/{{sub(/^>/, \">\" m \"|\")}} {{print}}' ${{faas[$i]}}; "
        "done > {output.tsv}.faa && "
        "diamond blastp -q {output.tsv}.faa -d {input.db} -o {output.tsv} "
        "--id {params.pid} --query-cover {params.qcov} --evalue {params.evalue} "
        "--max-target-seqs 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore && "
        "rm {output.tsv}.faa && date) &> >(tee {log})"

# collating BacMet hits (qseqid = "<mag>|<protein>")
rule collate_bacmet:
    input:
        expand(os.path.join(ANNOT_RESULTS_DIR, "amr/bacmet/batches/{batch}.bacmet.tsv"), batch=sorted(ANNOT_BATCHES))
    output:
        os.path.join(ANNOT_RESULTS_DIR, "amr/bacmet_all.tsv")
    message:
        "Collating BacMet results for {} MAGs".format(len(ANNOT_MAGS))
    run:
        cols = ["Genome", "qseqid", "sseqid", "bacmet_id", "bacmet_gene", "pident", "length",
                "qlen", "slen", "qcovhsp", "evalue", "bitscore"]
        with open(output[0], "w") as out:
            out.write("\t".join(cols) + "\n")
            for f in input:
                with open(f) as fh:
                    for line in fh:
                        p = line.rstrip("\n").split("\t")
                        mag, prot = p[0].split("|", 1)
                        sid = p[1].split("|")
                        out.write("\t".join([mag, prot, p[1], sid[0], sid[1] if len(sid) > 1 else ""] + p[2:]) + "\n")
