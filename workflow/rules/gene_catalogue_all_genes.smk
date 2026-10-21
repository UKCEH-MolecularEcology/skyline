"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-21]
Run: step "gene_catalogue_all_genes" in config.yaml (target: gene_catalogue_all_genes_all)
Purpose: Second gene catalogue that keeps ALL genes (no length filter), built with the 64-bit-key MMseqs2
         binary (https://mmseqs.com/latest/mmseqs-linux-avx2-64bit.tar.gz; MMseqs2 issue #1154): the standard
         build stores database keys as 32-bit integers and silently wraps above 2^32 (~4.29 billion) entries.
         Same workflow as rules/gene_catalogue.smk (exact -> 99/95/90/70/50 linclust -> RGI/BacMet/eggNOG on the
         annotated level -> inherited annotations at the other levels); rule names gca_*, outputs in
         <annotations>/gene_catalogue_all_genes/, logs in logs/gene_catalogue_all_genes/.
Disk: intermediates are removed as soon as nothing needs them (renamed/ after createdb; hash and per-level clustering
      DBs after each level; a level's representative DB once the next level is built; per-chunk annotation outputs after
      collation). db/ (sequence + header DBs, ~1 TB) is removed at the end unless gene_catalogue_all_genes.keep_db: true.
Config: section 'gene_catalogue_all_genes' overrides 'gene_catalogue' key by key (min_length 0, mmseqs_bin, ...).
"""

import os, glob, re
from snakemake.exceptions import WorkflowError

def _gca_merge(base, over):
    out = dict(base)
    for k, v in over.items():
        out[k] = _gca_merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out
GCA = _gca_merge(config["gene_catalogue"], config.get("gene_catalogue_all_genes") or {})
config["gene_catalogue_all_genes"] = GCA            # shell commands read {config[gene_catalogue_all_genes][mmseqs_bin]}
if not os.access(str(GCA.get("mmseqs_bin", "")), os.X_OK):
    raise WorkflowError(
        "gene_catalogue_all_genes.mmseqs_bin ('{}') is not an executable. Download the 64-bit-key MMseqs2 build:\n"
        "  wget https://mmseqs.com/latest/mmseqs-linux-avx2-64bit.tar.gz && tar xzf mmseqs-linux-avx2-64bit.tar.gz\n"
        "and set mmseqs_bin to <dir>/mmseqs/bin/mmseqs".format(GCA.get("mmseqs_bin")))
GCA_DIR = os.path.join(ANNOT_RESULTS_DIR, "gene_catalogue_all_genes")
GCA_LEVELS = [100] + sorted([int(l) for l in GCA["levels"] if int(l) < 100], reverse=True)   # 100 = exact
GCA_PREV = {l: GCA_LEVELS[i - 1] for i, l in enumerate(GCA_LEVELS) if i > 0}
GCA_ANNOT = int(GCA["annotate_level"])
assert GCA_ANNOT in GCA_LEVELS, "gene_catalogue_all_genes.annotate_level must be one of the clustering levels"

# ---- input assemblies: one entry per (prefix, path pattern with {name})
# NOTE: plain glob, not glob_wildcards -- glob_wildcards walks the whole tree below the pattern's
#       fixed prefix (os.walk), which on project scratch took >30 min at parse time.
#       '*' matches within one path component, i.e. the same as {name,[^/]+}.
def gca_find_inputs():
    found = {}
    for src in GCA["inputs"]:
        pat = src["pattern"]
        parts = [re.escape(p) for p in pat.split("{name}")]
        rx = re.compile("^" + parts[0] + "(?P<name>[^/]+)" + parts[1]
                        + "".join("(?P=name)" + p for p in parts[2:]) + "$")
        for path in sorted(glob.glob(pat.replace("{name}", "*"))):
            m = rx.match(path)
            if not m:
                continue
            sample = "{}_{}".format(src["prefix"], m.group("name"))
            if sample in GCA.get("exclude", []):
                continue
            if sample in found:
                raise ValueError("duplicate gene catalogue sample name: " + sample)
            found[sample] = path
    return found
GCA_INPUTS = gca_find_inputs()
GCA_SAMPLES = sorted(GCA_INPUTS)
GCA_DB = os.path.join(GCA_DIR, "db")
GCA_TOOLS = [t for t in ["rgi", "bacmet", "eggnog"] if t in GCA.get("annotate", ["rgi", "bacmet", "eggnog"])]
# levels that inherit the annotations of the annotated level (default: every clustering level except exact)
GCA_PROP = [int(l) for l in GCA.get("propagate_levels", [l for l in GCA_LEVELS if l not in (100, GCA_ANNOT)])]
assert all(l in GCA_LEVELS and l != GCA_ANNOT for l in GCA_PROP), \
    "gene_catalogue_all_genes.propagate_levels must be clustering levels other than annotate_level"

# ---- safety: this file deletes intermediates, so refuse to run if any path could reach inputs or shared data.
#      Deletions are limited to GCA_DIR (renamed/, db/, tmp/, per-chunk annotation files) and node-local
#      <local_tmpdir>/<user>_gca_* scratch; original assemblies are only read (via symlinks in input/).
_real = lambda p: os.path.realpath(os.path.abspath(str(p)))
_inside = lambda p, d: _real(p) == _real(d) or _real(p).startswith(_real(d).rstrip(os.sep) + os.sep)
if os.path.basename(GCA_DIR) != "gene_catalogue_all_genes" or not _inside(GCA_DIR, ANNOT_RESULTS_DIR) \
        or _real(GCA_DIR) == _real(ANNOT_RESULTS_DIR):
    raise WorkflowError("unsafe gene_catalogue_all_genes output directory: " + GCA_DIR)
_protected = list(GCA_INPUTS.values()) + [GCA["eggnog"]["db"], GCA["eggnog"]["sif"],
              os.path.join(ANNOT_DBS_DIR, "rgi"), os.path.join(ANNOT_DBS_DIR, "bacmet")]
for _p in _protected:
    if _inside(_p, GCA_DIR) or _inside(GCA_DIR, _p):
        raise WorkflowError("input or shared database overlaps the all-genes output directory: " + str(_p))
if _real(GCA.get("local_tmpdir", "/tmp")) in ("/", _real(ANNOT_RESULTS_DIR)):
    raise WorkflowError("gene_catalogue_all_genes.local_tmpdir must be a scratch directory, not " + GCA.get("local_tmpdir"))

localrules: gene_catalogue_all_genes_all, gca_link_inputs, gca_summary

###################
# RULES
###################
rule gene_catalogue_all_genes_all:
    input:
        os.path.join(GCA_DIR, "clusters", "summary.tsv"),
        expand(os.path.join(GCA_DIR, "annotation", str(GCA_ANNOT), "{tool}_genes.tsv"), tool=GCA_TOOLS),
        expand(os.path.join(GCA_DIR, "annotation", "{level}", "{tool}_reps.tsv.gz"), level=GCA_PROP, tool=GCA_TOOLS)
    output:
        touch("status/gene_catalogue_all_genes.done")
    params:
        db=GCA_DB,
        tmp=os.path.join(GCA_DIR, "tmp"),
        keep_db=GCA.get("keep_db", False)
    run:
        import shutil
        if not params.keep_db:
            for d in (params.db, params.tmp):
                if not os.path.lexists(d):
                    continue
                if os.path.islink(d) or not _inside(d, GCA_DIR) or _real(d) == _real(GCA_DIR):
                    raise WorkflowError("refusing to remove " + d)
                shutil.rmtree(d)     # does not follow symlinks inside (subdb data links are removed, not their targets)


################### inputs
# symlinks named by sample (sr_* short-read per-sample, lr_* long-read assemblies) + manifest
rule gca_link_inputs:
    output:
        manifest=os.path.join(GCA_DIR, "input", "manifest.tsv")
    message:
        "Gene catalogue (all genes): linking {} assemblies".format(len(GCA_SAMPLES))
    run:
        d = os.path.dirname(output.manifest)
        os.makedirs(d, exist_ok=True)
        with open(output.manifest, "w") as out:
            out.write("sample\tsource\n")
            for s in GCA_SAMPLES:
                link = os.path.join(d, s + ".faa")
                if os.path.islink(link):
                    os.remove(link)          # removes the link only, never the original assembly
                elif os.path.exists(link):
                    raise WorkflowError("not a symlink, refusing to replace: " + link)
                os.symlink(GCA_INPUTS[s], link)
                out.write("{}\t{}\n".format(s, GCA_INPUTS[s]))

# unique gene IDs: '<sample>|<first word of the Prodigal header>', trailing '*' removed, optional min length
rule gca_rename:
    input:
        manifest=rules.gca_link_inputs.output.manifest
    output:
        faa=temp(os.path.join(GCA_DIR, "renamed", "{sample}.faa.gz")),   # deleted once createdb has read it
        n_genes=os.path.join(GCA_DIR, "renamed", "{sample}.count")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "rename_{sample}.log")
    params:
        src=lambda wildcards: os.path.join(GCA_DIR, "input", wildcards.sample + ".faa"),
        min_len=GCA.get("min_length", 0)
    threads:
        2
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    wildcard_constraints:
        sample="|".join(GCA_SAMPLES) if GCA_SAMPLES else "NONE"
    message:
        "Gene catalogue (all genes): unique IDs for {wildcards.sample}"
    shell:
        "(date && awk -v s={wildcards.sample} -v m={params.min_len} -v c={output.n_genes} '"
        "function flush() {{ if (id != \"\" && length(seq) >= m) {{ print \">\" s \"|\" id; print seq; n++ }} }} "
        "/^>/ {{ flush(); split(substr($0, 2), a, \" \"); id = a[1]; seq = \"\"; next }} "
        "{{ gsub(/[*[:space:]]/, \"\"); seq = seq $0 }} "
        "END {{ flush(); print n + 0 > c }}' {params.src} | pigz -p {threads} > {output.faa} && "
        "echo \"genes: $(cat {output.n_genes})\" && date) &> >(tee {log})"


################### MMseqs2 database + clustering
rule gca_createdb:
    input:
        expand(os.path.join(GCA_DIR, "renamed", "{sample}.faa.gz"), sample=GCA_SAMPLES)
    output:
        dbtype=temp(os.path.join(GCA_DB, "genes.dbtype")),
        n_genes=os.path.join(GCA_DIR, "clusters", "genes.count")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "createdb.log")
    params:
        db=os.path.join(GCA_DB, "genes")
    threads:
        4    # createdb is essentially single-threaded (3.5% CPU efficiency at 16)
    # --shuffle 0: skips createRenumberedDB (single-threaded renumbering of every entry)
    # --write-lookup 0: the .lookup is not used downstream (names come from genes_h)
    # 64-bit-key MMseqs2 build (mmseqs_bin): no 2^32 limit on the number of genes
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue (all genes): MMseqs2 database of {} assemblies".format(len(GCA_SAMPLES))
    shell:
        "(date && mkdir -p $(dirname {params.db}) && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} createdb {input} {params.db} --dbtype 1 --shuffle 0 --write-lookup 0 --threads {threads} && "
        "wc -l < {params.db}.index > {output.n_genes} && "
        "echo \"genes: $(cat {output.n_genes})\" && date) &> >(tee {log})"

# level 100: identical sequences only (hash-based)
rule gca_exact:
    input:
        rules.gca_createdb.output.dbtype
    output:
        reps=temp(os.path.join(GCA_DB, "reps_100.dbtype")),
        tsv=os.path.join(GCA_DIR, "clusters", "clusters_100.tsv.gz"),
        n_reps=os.path.join(GCA_DIR, "clusters", "reps_100.count")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "cluster_100.log")
    params:
        db=os.path.join(GCA_DB, "genes"),
        hash=os.path.join(GCA_DB, "hash_100"),
        clu=os.path.join(GCA_DB, "clu_100"),
        reps=os.path.join(GCA_DB, "reps_100")
    threads:
        GCA["mmseqs"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue (all genes): exact dereplication (100%)"
    # clusthash (~3.5 h at 4 billion genes) is skipped on a rerun if it already completed;
    # clust (set cover) holds the whole hash graph in memory (~1 TB at 4 billion genes)
    shell:
        "(date && "
        "if [ -e {params.hash}.dbtype ]; then echo \"reusing existing {params.hash}\"; "
        "else {config[gene_catalogue_all_genes][mmseqs_bin]} clusthash {params.db} {params.hash} --min-seq-id 1.0 --threads {threads}; fi && "
        "rm -f {params.clu} {params.clu}.index {params.clu}.dbtype {params.clu}.[0-9]* && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} clust {params.db} {params.hash} {params.clu} --threads {threads} && "
        "rm -f {params.hash} {params.hash}.index {params.hash}.dbtype {params.hash}.[0-9]* && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} createsubdb {params.clu} {params.db} {params.reps} --subdb-mode 1 && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} createtsv {params.db} {params.db} {params.clu} {params.clu}.tsv --threads {threads} && "
        "pigz -p {threads} -c {params.clu}.tsv > {output.tsv} && rm {params.clu}.tsv && "
        "rm -f {params.clu} {params.clu}.index {params.clu}.dbtype {params.clu}.[0-9]* && "
        "wc -l < {params.reps}.index > {output.n_reps} && "
        "echo \"representatives: $(cat {output.n_reps})\" && date) &> >(tee {log})"

# levels < 100: linclust on the previous level's representatives
rule gca_linclust:
    input:
        lambda wildcards: os.path.join(GCA_DB, "reps_{}.dbtype".format(GCA_PREV[int(wildcards.level)]))
    output:
        reps=temp(os.path.join(GCA_DB, "reps_{level}.dbtype")),
        tsv=os.path.join(GCA_DIR, "clusters", "clusters_{level}.tsv.gz"),
        n_reps=os.path.join(GCA_DIR, "clusters", "reps_{level}.count")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "cluster_{level}.log")
    wildcard_constraints:
        level="|".join(str(l) for l in GCA_LEVELS if l < 100)
    params:
        prev=lambda wildcards: os.path.join(GCA_DB, "reps_{}".format(GCA_PREV[int(wildcards.level)])),
        clu=lambda wildcards: os.path.join(GCA_DB, "clu_{}".format(wildcards.level)),
        reps=lambda wildcards: os.path.join(GCA_DB, "reps_{}".format(wildcards.level)),
        tmp=lambda wildcards: os.path.join(GCA_DIR, "tmp", "linclust_{}".format(wildcards.level)),
        ident=lambda wildcards: int(wildcards.level) / 100,
        cov=GCA["mmseqs"]["coverage"],
        cov_mode=GCA["mmseqs"]["cov_mode"],
        mem=GCA["mmseqs"]["split_memory_limit"],
        extra=GCA["mmseqs"].get("extra", ""),
        # the previous level's representative DB is only needed by this level -- except the annotated level's,
        # which gca_reps_fasta also reads (left in place; small)
        drop_prev=lambda wildcards: "" if GCA_PREV[int(wildcards.level)] == GCA_ANNOT else "yes"
    threads:
        GCA["mmseqs"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    message:
        "Gene catalogue (all genes): clustering at {wildcards.level}% identity"
    shell:
        "(date && rm -rf {params.tmp} && mkdir -p {params.tmp} && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} linclust {params.prev} {params.clu} {params.tmp} --min-seq-id {params.ident} -c {params.cov} "
        "--cov-mode {params.cov_mode} --split-memory-limit {params.mem} --threads {threads} {params.extra} && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} createsubdb {params.clu} {params.prev} {params.reps} --subdb-mode 1 && "
        "{config[gene_catalogue_all_genes][mmseqs_bin]} createtsv {params.prev} {params.prev} {params.clu} {params.clu}.tsv --threads {threads} && "
        "pigz -p {threads} -c {params.clu}.tsv > {output.tsv} && rm {params.clu}.tsv && rm -rf {params.tmp} && "
        "rm -f {params.clu} {params.clu}.index {params.clu}.dbtype {params.clu}.[0-9]* && "
        "wc -l < {params.reps}.index > {output.n_reps} && "
        "if [ -n \"{params.drop_prev}\" ]; then rm -f {params.prev} {params.prev}.index {params.prev}.dbtype "
        "{params.prev}_h {params.prev}_h.index {params.prev}_h.dbtype; fi && "
        "echo \"representatives: $(cat {output.n_reps})\" && date) &> >(tee {log})"

rule gca_summary:
    input:
        genes=rules.gca_createdb.output.n_genes,
        reps=expand(os.path.join(GCA_DIR, "clusters", "reps_{level}.count"), level=GCA_LEVELS)
    output:
        os.path.join(GCA_DIR, "clusters", "summary.tsv")
    message:
        "Gene catalogue (all genes): summary"
    run:
        n = lambda p: int(open(p).read().strip())
        with open(output[0], "w") as out:
            out.write("level\trepresentatives\tclustering_input\n")
            out.write("genes\t{}\t\n".format(n(input.genes)))
            for l in GCA_LEVELS:
                out.write("{}\t{}\t{}\n".format(l, n(os.path.join(GCA_DIR, "clusters", "reps_{}.count".format(l))),
                                                "genes" if l == 100 else "reps_{}".format(GCA_PREV[l])))


################### annotation of the representatives (one level), in chunks
rule gca_reps_fasta:
    input:
        os.path.join(GCA_DB, "reps_{level}.dbtype")
    output:
        os.path.join(GCA_DIR, "reps", "reps_{level}.faa.gz")
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "reps_fasta_{level}.log")
    params:
        reps=lambda wildcards: os.path.join(GCA_DB, "reps_{}".format(wildcards.level))
    threads:
        8
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    shell:
        "(date && {config[gene_catalogue_all_genes][mmseqs_bin]} convert2fasta {params.reps} {output}.tmp.faa && "
        "pigz -p {threads} -c {output}.tmp.faa > {output} && rm {output}.tmp.faa && date) &> >(tee {log})"

checkpoint gca_split_reps:
    input:
        rules.gca_reps_fasta.output
    output:
        directory(os.path.join(GCA_DIR, "annotation", "{level}", "chunks"))
    params:
        n=GCA["chunk_size"]
    shell:
        "mkdir -p {output} && zcat {input} | awk -v n={params.n} -v d={output} '"
        "/^>/ {{ if (c % n == 0) {{ if (f) close(f); f = sprintf(\"%s/chunk_%05d.faa\", d, int(c / n)) }} c++ }} "
        "{{ print > f }}'"

def gca_chunks(wildcards):
    d = checkpoints.gca_split_reps.get(level=wildcards.level).output[0]
    return sorted(glob_wildcards(os.path.join(d, "chunk_{chunk}.faa")).chunk)

rule gca_rgi_chunk:
    input:
        faa=os.path.join(GCA_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        setup="status/rgi_setup.done"
    output:
        temp(os.path.join(GCA_DIR, "annotation", "{level}", "rgi", "chunk_{chunk}.txt"))
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "rgi_{level}_{chunk}.log")
    params:
        db=os.path.join(ANNOT_DBS_DIR, "rgi"),
        prefix=lambda wildcards, output: output[0][:-len(".txt")]
    threads:
        GCA["rgi_threads"]
    container:
        config["rgi"]["container"]
    shell:
        "(date && mkdir -p $(dirname {output}) && cd {params.db} && "
        "rgi main --input_sequence {input.faa} --output_file {params.prefix} --input_type protein --local "
        "-a DIAMOND --clean -n {threads} && rm -f {params.prefix}.json && date) &> >(tee {log})"

rule gca_bacmet_chunk:
    input:
        faa=os.path.join(GCA_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        db=os.path.join(ANNOT_DBS_DIR, "bacmet", "bacmet_exp.dmnd")
    output:
        temp(os.path.join(GCA_DIR, "annotation", "{level}", "bacmet", "chunk_{chunk}.tsv"))
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "bacmet_{level}_{chunk}.log")
    params:
        pid=config["bacmet"]["id"],
        qcov=config["bacmet"]["qcov"],
        evalue=config["bacmet"]["evalue"]
    threads:
        GCA["bacmet_threads"]
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    shell:
        "(date && diamond blastp -q {input.faa} -d {input.db} -o {output} "
        "--id {params.pid} --query-cover {params.qcov} --evalue {params.evalue} "
        "--max-target-seqs 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore && date) &> >(tee {log})"

rule gca_eggnog_chunk:
    input:
        faa=os.path.join(GCA_DIR, "annotation", "{level}", "chunks", "chunk_{chunk}.faa"),
        sif=GCA["eggnog"]["sif"]
    output:
        temp(os.path.join(GCA_DIR, "annotation", "{level}", "eggnog", "chunk_{chunk}.emapper.annotations"))
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "eggnog_{level}_{chunk}.log")
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        db=GCA["eggnog"]["db"],
        tmp=lambda wildcards: os.path.join(GCA.get("local_tmpdir", "/tmp"),
                                           "{}_gca_emapper_{}_{}".format(os.environ.get("USER", "user"), wildcards.level, wildcards.chunk)),
        extra=GCA["eggnog"].get("extra", "")
    threads:
        GCA["eggnog"]["threads"]
    # eggNOG annotation databases are copied to node-local {params.tmp}/db: SQLite lookups over NFS
    # made every annotation worker time out (2 s/seed limit). The DIAMOND db stays on NFS (symlink).
    # If a previous attempt already finished the DIAMOND search, its seed_orthologs are reused
    # (kept as .seed_orthologs.keep so --override cannot remove them) and only annotation runs.
    shell:
        "(set -e; trap 'rm -rf {params.tmp}' EXIT; date; "
        "mkdir -p {params.outdir} {params.tmp}/db {params.tmp}/work; "
        "for f in {params.db}/*; do case $(basename $f) in "
        "mmseqs|*.sif|eggnog_proteins.dmnd|download.*) ;; *) cp -r $f {params.tmp}/db/ ;; esac; done; "
        "ln -s {params.db}/eggnog_proteins.dmnd {params.tmp}/db/; "
        "echo \"local eggNOG db: $(du -sh {params.tmp}/db | cut -f1)\"; date; "
        "S={params.outdir}/chunk_{wildcards.chunk}.emapper.seed_orthologs; "
        "if [ -s $S ]; then mv $S $S.keep; fi; "
        "if [ -s $S.keep ]; then echo \"reusing DIAMOND hits: $S.keep\"; "
        "MODE=\"-m no_search --annotate_hits_table $S.keep\"; "
        "else MODE=\"-i {input.faa}\"; fi; "
        "singularity exec -B /ei {input.sif} emapper.py $MODE --output_dir {params.outdir} "
        "-o chunk_{wildcards.chunk} --data_dir {params.tmp}/db --cpu {threads} --temp_dir {params.tmp}/work "
        "--override {params.extra}; "
        "rm -f {params.outdir}/chunk_{wildcards.chunk}.emapper.hits $S $S.keep; "
        "date) &> >(tee {log})"

# collate per-chunk outputs (one header)
def gca_chunk_files(tool, ext):
    def f(wildcards):
        return [os.path.join(GCA_DIR, "annotation", wildcards.level, tool, "chunk_{}.{}".format(c, ext))
                for c in gca_chunks(wildcards)]
    return f

rule gca_collate:
    input:
        rgi=gca_chunk_files("rgi", "txt")
    output:
        os.path.join(GCA_DIR, "annotation", "{level}", "rgi_genes.tsv")
    shell:
        "awk 'FNR == 1 && NR != 1 {{ next }} {{ print }}' {input.rgi} > {output}"

rule gca_collate_bacmet:
    input:
        gca_chunk_files("bacmet", "tsv")
    output:
        os.path.join(GCA_DIR, "annotation", "{level}", "bacmet_genes.tsv")
    shell:
        "(printf 'gene\\tsseqid\\tpident\\tlength\\tqlen\\tslen\\tqcovhsp\\tevalue\\tbitscore\\n'; cat {input}) > {output}"

rule gca_collate_eggnog:
    input:
        gca_chunk_files("eggnog", "emapper.annotations")
    output:
        os.path.join(GCA_DIR, "annotation", "{level}", "eggnog_genes.tsv")
    shell:
        "awk '/^##/ {{ next }} /^#query/ {{ if (h++) next; sub(/^#/, \"\") }} {{ print }}' {input} > {output}"


################### annotations at the other levels (inherited from the annotated level)
# below annotate_level (70, 50): the representatives are annotated representatives -> their rows
# above annotate_level (95, 99, 100): rep_<level> -> its annotated representative (through the cluster chain)
#   + that representative's annotation; also writes annotation/<level>/reps_to_<annotate_level>.tsv.gz
def gca_prop_clusters(wildcards):
    l = int(wildcards.level)
    lv = [l] if l < GCA_ANNOT else [x for x in GCA_LEVELS if GCA_ANNOT <= x < l]
    return [os.path.join(GCA_DIR, "clusters", "clusters_{}.tsv.gz".format(x)) for x in lv]

rule gca_propagate:
    input:
        clusters=gca_prop_clusters,
        ann=expand(os.path.join(GCA_DIR, "annotation", str(GCA_ANNOT), "{tool}_genes.tsv"), tool=GCA_TOOLS)
    output:
        expand(os.path.join(GCA_DIR, "annotation", "{{level}}", "{tool}_reps.tsv.gz"), tool=GCA_TOOLS)
    log:
        os.path.join(ANNOT_RESULTS_DIR, "logs", "gene_catalogue_all_genes", "propagate_{level}.log")
    params:
        script=os.path.join(SRC_DIR, "gc_propagate.py"),
        levels=",".join(str(l) for l in GCA_LEVELS),
        tools=",".join(GCA_TOOLS),
        clusters_dir=os.path.join(GCA_DIR, "clusters"),
        annot_dir=os.path.join(GCA_DIR, "annotation", str(GCA_ANNOT)),
        out_dir=lambda wildcards: os.path.join(GCA_DIR, "annotation", wildcards.level),
        memory=GCA.get("propagate_memory", "300GB"),
        tmp=lambda wildcards: os.path.join(GCA.get("local_tmpdir", "/tmp"), "{}_gca_propagate_{}".format(os.environ.get("USER", "user"), wildcards.level))
    threads:
        GCA.get("propagate_threads", 32)
    conda:
        os.path.join(ENV_DIR, "duckdb.yaml")
    wildcard_constraints:
        level="|".join(str(l) for l in GCA_PROP) if GCA_PROP else "NONE"
    message:
        "Gene catalogue (all genes): level-{wildcards.level} annotations inherited from level " + str(GCA_ANNOT)
    shell:
        "(trap 'rm -rf {params.tmp}' EXIT; date && "
        "python {params.script} --level {wildcards.level} --annot-level " + str(GCA_ANNOT) + " "
        "--levels {params.levels} --clusters-dir {params.clusters_dir} --annot-dir {params.annot_dir} "
        "--tools {params.tools} --out-dir {params.out_dir} --tmp {params.tmp} "
        "--threads {threads} --memory {params.memory} && date) &> >(tee {log})"

