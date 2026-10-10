"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-09]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --use-singularity -rp gene_catalogue_all
Purpose: Gene catalogue from all short-read (per-sample) and long-read assemblies (Prodigal contigs.faa):
         symlink inputs -> unique gene IDs (<sample>|<gene>) -> exact dereplication (MMseqs2 clusthash) ->
         cascaded clustering 99 / 95 / 90 / 70 / 50 % (MMseqs2 linclust on the previous level's representatives) ->
         gene-level annotation of the representatives at one level (RGI protein mode, BacMet, eggNOG-mapper v3 /
         eggNOG 7), in chunks.
Outputs (<annotations>/gene_catalogue/):
    input/<sample>.faa                      symlinks to the original contigs.faa
    renamed/<sample>.faa.gz                 headers '<sample>|<gene>'
    clusters/clusters_<L>.tsv.gz            representative -> member (members are representatives of the level above;
                                            L=100: exact duplicates of the original genes)
    clusters/summary.tsv                    representatives per level
    reps/reps_<L>.faa.gz                    representative sequences of the annotated level
    annotation/<L>/{rgi,bacmet,eggnog}_genes.tsv   collated per-representative annotations
"""

import os, glob

GC = config["gene_catalogue"]
GC_DIR = os.path.join(ANNOT_RESULTS_DIR, "gene_catalogue")
GC_LEVELS = [100] + sorted([int(l) for l in GC["levels"] if int(l) < 100], reverse=True)   # 100 = exact
GC_PREV = {l: GC_LEVELS[i - 1] for i, l in enumerate(GC_LEVELS) if i > 0}
GC_ANNOT = int(GC["annotate_level"])
assert GC_ANNOT in GC_LEVELS, "gene_catalogue.annotate_level must be one of the clustering levels"

# ---- input assemblies: one entry per (prefix, path pattern with {name})
def gc_find_inputs():
    found = {}
    for src in GC["inputs"]:
        for name in glob_wildcards(src["pattern"].replace("{name}", "{name,[^/]+}")).name:
            sample = "{}_{}".format(src["prefix"], name)
            if sample in GC.get("exclude", []):
                continue
            if sample in found:
                raise ValueError("duplicate gene catalogue sample name: " + sample)
            found[sample] = src["pattern"].format(name=name)
    return found
GC_INPUTS = gc_find_inputs()
GC_SAMPLES = sorted(GC_INPUTS)
GC_DB = os.path.join(GC_DIR, "db")

localrules: gene_catalogue_all, gc_link_inputs, gc_summary, gc_download_eggnog_sif

###################
# RULES
###################
rule gene_catalogue_all:
    input:
        os.path.join(GC_DIR, "clusters", "summary.tsv"),
        expand(os.path.join(GC_DIR, "annotation", str(GC_ANNOT), "{tool}_genes.tsv"),
               tool=[t for t in ["rgi", "bacmet", "eggnog"] if t in GC.get("annotate", ["rgi", "bacmet", "eggnog"])])
    output:
        touch("status/gene_catalogue.done")


################### inputs
# symlinks named by sample (sr_* short-read per-sample, lr_* long-read assemblies) + manifest
rule gc_link_inputs:
    output:
        manifest=os.path.join(GC_DIR, "input", "manifest.tsv")
    message:
        "Gene catalogue: linking {} assemblies".format(len(GC_SAMPLES))
    run:
        d = os.path.dirname(output.manifest)
        os.makedirs(d, exist_ok=True)
        with open(output.manifest, "w") as out:
            out.write("sample\tsource\n")
            for s in GC_SAMPLES:
                link = os.path.join(d, s + ".faa")
                if os.path.lexists(link):
                    os.remove(link)
                os.symlink(GC_INPUTS[s], link)
                out.write("{}\t{}\n".format(s, GC_INPUTS[s]))

# unique gene IDs: '<sample>|<first word of the Prodigal header>', trailing '*' removed, optional min length
rule gc_rename:
    input:
        manifest=rules.gc_link_inputs.output.manifest
    output:
        faa=os.path.join(GC_DIR, "renamed", "{sample}.faa.gz"),
        n_genes=os.path.join(GC_DIR, "renamed", "{sample}.count")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "rename_{sample}.log")
    params:
        src=lambda wildcards: os.path.join(GC_DIR, "input", wildcards.sample + ".faa"),
        min_len=GC.get("min_length", 0)
    threads:
        4
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    wildcard_constraints:
        sample="|".join(GC_SAMPLES) if GC_SAMPLES else "NONE"
    message:
        "Gene catalogue: unique IDs for {wildcards.sample}"
    shell:
        "(date && awk -v s={wildcards.sample} -v m={params.min_len} -v c={output.n_genes} '"
        "function flush() {{ if (id != \"\" && length(seq) >= m) {{ print \">\" s \"|\" id; print seq; n++ }} }} "
        "/^>/ {{ flush(); split(substr($0, 2), a, \" \"); id = a[1]; seq = \"\"; next }} "
        "{{ gsub(/[*[:space:]]/, \"\"); seq = seq $0 }} "
        "END {{ flush(); print n + 0 > c }}' {params.src} | pigz -p {threads} > {output.faa} && "
        "echo \"genes: $(cat {output.n_genes})\" && date) &> >(tee {log})"


################### MMseqs2 database + clustering
rule gc_createdb:
    input:
        expand(os.path.join(GC_DIR, "renamed", "{sample}.faa.gz"), sample=GC_SAMPLES)
    output:
        os.path.join(GC_DB, "genes.dbtype")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "createdb.log")
    params:
        db=os.path.join(GC_DB, "genes")
    threads:
        GC["mmseqs"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue: MMseqs2 database of {} assemblies".format(len(GC_SAMPLES))
    shell:
        "(date && mkdir -p $(dirname {params.db}) && "
        "mmseqs createdb {input} {params.db} --dbtype 1 && "
        "echo \"genes: $(wc -l < {params.db}.index)\" && date) &> >(tee {log})"

# level 100: identical sequences only (hash-based; low memory)
rule gc_exact:
    input:
        rules.gc_createdb.output
    output:
        reps=os.path.join(GC_DB, "reps_100.dbtype"),
        tsv=os.path.join(GC_DIR, "clusters", "clusters_100.tsv.gz")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "cluster_100.log")
    params:
        db=os.path.join(GC_DB, "genes"),
        hash=os.path.join(GC_DB, "hash_100"),
        clu=os.path.join(GC_DB, "clu_100"),
        reps=os.path.join(GC_DB, "reps_100")
    threads:
        GC["mmseqs"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue: exact dereplication (100%)"
    shell:
        "(date && mmseqs clusthash {params.db} {params.hash} --min-seq-id 1.0 --threads {threads} && "
        "mmseqs clust {params.db} {params.hash} {params.clu} --threads {threads} && "
        "mmseqs createsubdb {params.clu} {params.db} {params.reps} --subdb-mode 1 && "
        "mmseqs createtsv {params.db} {params.db} {params.clu} {params.clu}.tsv --threads {threads} && "
        "pigz -p {threads} -c {params.clu}.tsv > {output.tsv} && rm {params.clu}.tsv && "
        "echo \"representatives: $(wc -l < {params.reps}.index)\" && date) &> >(tee {log})"

# levels < 100: linclust on the previous level's representatives
rule gc_linclust:
    input:
        lambda wildcards: os.path.join(GC_DB, "reps_{}.dbtype".format(GC_PREV[int(wildcards.level)]))
    output:
        reps=os.path.join(GC_DB, "reps_{level}.dbtype"),
        tsv=os.path.join(GC_DIR, "clusters", "clusters_{level}.tsv.gz")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "cluster_{level}.log")
    wildcard_constraints:
        level="|".join(str(l) for l in GC_LEVELS if l < 100)
    params:
        prev=lambda wildcards: os.path.join(GC_DB, "reps_{}".format(GC_PREV[int(wildcards.level)])),
        clu=lambda wildcards: os.path.join(GC_DB, "clu_{}".format(wildcards.level)),
        reps=lambda wildcards: os.path.join(GC_DB, "reps_{}".format(wildcards.level)),
        tmp=lambda wildcards: os.path.join(GC_DIR, "tmp", "linclust_{}".format(wildcards.level)),
        ident=lambda wildcards: int(wildcards.level) / 100,
        cov=GC["mmseqs"]["coverage"],
        cov_mode=GC["mmseqs"]["cov_mode"],
        mem=GC["mmseqs"]["split_memory_limit"],
        extra=GC["mmseqs"].get("extra", "")
    threads:
        GC["mmseqs"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue: clustering at {wildcards.level}% identity"
    shell:
        "(date && rm -rf {params.tmp} && mkdir -p {params.tmp} && "
        "mmseqs linclust {params.prev} {params.clu} {params.tmp} --min-seq-id {params.ident} -c {params.cov} "
        "--cov-mode {params.cov_mode} --split-memory-limit {params.mem} --threads {threads} {params.extra} && "
        "mmseqs createsubdb {params.clu} {params.prev} {params.reps} --subdb-mode 1 && "
        "mmseqs createtsv {params.prev} {params.prev} {params.clu} {params.clu}.tsv --threads {threads} && "
        "pigz -p {threads} -c {params.clu}.tsv > {output.tsv} && rm {params.clu}.tsv && rm -rf {params.tmp} && "
        "echo \"representatives: $(wc -l < {params.reps}.index)\" && date) &> >(tee {log})"

rule gc_summary:
    input:
        genes=rules.gc_createdb.output,
        reps=expand(os.path.join(GC_DB, "reps_{level}.dbtype"), level=GC_LEVELS)
    output:
        os.path.join(GC_DIR, "clusters", "summary.tsv")
    message:
        "Gene catalogue: summary"
    run:
        n = lambda p: sum(1 for _ in open(p[:-len(".dbtype")] + ".index"))
        with open(output[0], "w") as out:
            out.write("level\trepresentatives\tclustering_input\n")
            out.write("genes\t{}\t\n".format(n(input.genes[0])))
            for l in GC_LEVELS:
                out.write("{}\t{}\t{}\n".format(l, n(os.path.join(GC_DB, "reps_{}.dbtype".format(l))),
                                                "genes" if l == 100 else "reps_{}".format(GC_PREV[l])))


################### annotation of the representatives (one level), in chunks
rule gc_reps_fasta:
    input:
        os.path.join(GC_DB, "reps_{level}.dbtype")
    output:
        os.path.join(GC_DIR, "reps", "reps_{level}.faa.gz")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "reps_fasta_{level}.log")
    params:
        reps=lambda wildcards: os.path.join(GC_DB, "reps_{}".format(wildcards.level))
    threads:
        8
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    shell:
        "(date && mmseqs convert2fasta {params.reps} {output}.tmp.faa && "
        "pigz -p {threads} -c {output}.tmp.faa > {output} && rm {output}.tmp.faa && date) &> >(tee {log})"

checkpoint gc_split_reps:
    input:
        rules.gc_reps_fasta.output
    output:
        directory(os.path.join(GC_DIR, "annotation", "{level}", "chunks"))
    params:
        n=GC["chunk_size"]
    shell:
        "mkdir -p {output} && zcat {input} | awk -v n={params.n} -v d={output} '"
        "/^>/ {{ if (c % n == 0) {{ if (f) close(f); f = sprintf(\"%s/chunk_%05d.faa\", d, int(c / n)) }} c++ }} "
        "{{ print > f }}'"

def gc_chunks(wildcards):
    d = checkpoints.gc_split_reps.get(level=wildcards.level).output[0]
    return sorted(glob_wildcards(os.path.join(d, "chunk_{chunk}.faa")).chunk)

rule gc_rgi_chunk:
    input:
        faa=os.path.join(GC_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        setup="status/rgi_setup.done"
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "rgi", "chunk_{chunk}.txt")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "rgi_{level}_{chunk}.log")
    params:
        db=os.path.join(ANNOT_DBS_DIR, "rgi"),
        prefix=lambda wildcards, output: output[0][:-len(".txt")]
    threads:
        GC["rgi_threads"]
    container:
        config["rgi"]["container"]
    shell:
        "(date && mkdir -p $(dirname {output}) && cd {params.db} && "
        "rgi main --input_sequence {input.faa} --output_file {params.prefix} --input_type protein --local "
        "-a DIAMOND --clean -n {threads} && date) &> >(tee {log})"

rule gc_bacmet_chunk:
    input:
        faa=os.path.join(GC_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        db=os.path.join(ANNOT_DBS_DIR, "bacmet", "bacmet_exp.dmnd")
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "bacmet", "chunk_{chunk}.tsv")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "bacmet_{level}_{chunk}.log")
    params:
        pid=config["bacmet"]["id"],
        qcov=config["bacmet"]["qcov"],
        evalue=config["bacmet"]["evalue"]
    threads:
        GC["bacmet_threads"]
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    shell:
        "(date && diamond blastp -q {input.faa} -d {input.db} -o {output} "
        "--id {params.pid} --query-cover {params.qcov} --evalue {params.evalue} "
        "--max-target-seqs 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore && date) &> >(tee {log})"

# eggNOG-mapper v3 (eggNOG 7) official Singularity image, downloaded as a SIF (no image build)
rule gc_download_eggnog_sif:
    output:
        GC["eggnog"]["sif"]
    params:
        url=GC["eggnog"]["sif_url"]
    shell:
        "mkdir -p $(dirname {output}) && wget -q -O {output}.tmp {params.url} && mv {output}.tmp {output}"

rule gc_eggnog_chunk:
    input:
        faa=os.path.join(GC_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        sif=GC["eggnog"]["sif"]
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "eggnog", "chunk_{chunk}.emapper.annotations")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue", "eggnog_{level}_{chunk}.log")
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        db=GC["eggnog"]["db"],
        tmp=lambda wildcards: os.path.join(GC.get("local_tmpdir", "/tmp"),
                                           "{}_emapper_{}_{}".format(os.environ.get("USER", "user"), wildcards.level, wildcards.chunk)),
        extra=GC["eggnog"].get("extra", "")
    threads:
        GC["eggnog"]["threads"]
    shell:
        "(date && mkdir -p {params.outdir} {params.tmp} && "
        "singularity exec {input.sif} emapper.py -i {input.faa} --output_dir {params.outdir} "
        "-o chunk_{wildcards.chunk} --data_dir {params.db} --cpu {threads} --temp_dir {params.tmp} --override {params.extra} && "
        "rm -rf {params.tmp} && date) &> >(tee {log})"

# collate per-chunk outputs (one header)
def gc_chunk_files(tool, ext):
    def f(wildcards):
        return [os.path.join(GC_DIR, "annotation", wildcards.level, tool, "chunk_{}.{}".format(c, ext))
                for c in gc_chunks(wildcards)]
    return f

rule gc_collate:
    input:
        rgi=gc_chunk_files("rgi", "txt")
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "rgi_genes.tsv")
    shell:
        "awk 'FNR == 1 && NR != 1 {{ next }} {{ print }}' {input.rgi} > {output}"

rule gc_collate_bacmet:
    input:
        gc_chunk_files("bacmet", "tsv")
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "bacmet_genes.tsv")
    shell:
        "(printf 'gene\\tsseqid\\tpident\\tlength\\tqlen\\tslen\\tqcovhsp\\tevalue\\tbitscore\\n'; cat {input}) > {output}"

rule gc_collate_eggnog:
    input:
        gc_chunk_files("eggnog", "emapper.annotations")
    output:
        os.path.join(GC_DIR, "annotation", "{level}", "eggnog_genes.tsv")
    shell:
        "awk '/^##/ {{ next }} /^#query/ {{ if (h++) next; sub(/^#/, \"\") }} {{ print }}' {input} > {output}"
