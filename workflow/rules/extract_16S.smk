"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/kraken2.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on reads
"""

import os

localrules: install_hyperex, extract_16S_all

###################
# RULES
###################
rule extract_16S_all:
    input:
        expand(f"{RESULTS_DIR}/{{sample}}/{{sample}}_16S_seqs.fa", sample=SAMPLES),
        "submodules/hyperex_installed.txt",
        os.path.join(RESULTS_DIR, "hyperex/extracted_16S.fa"),
        os.path.join(RESULTS_DIR, "extracted/extracted_16S.fa")
    output:
        touch("status/extract_16S.done")        

# Renaming
rule rename_16S_headers:
    input:
        lambda wildcards: os.path.join(ASS_DIR, wildcards.sample, "metamdbg", "annotation", "16S_seqs.fa")
    output:
        renamed = f"{RESULTS_DIR}/{{sample}}/{{sample}}_16S_seqs.fa"
    wildcard_constraints:
        sid="|".join(SAMPLES)
    message:
        "Copying and renaming the fasta headers for the 16S from assembly for: {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.renamed})
        awk -v prefix="{wildcards.sample}" '
            /^>/ {{sub(/^>/, ">" prefix "-"); print; next}}
            {{print}}
        ' {input} > {output.renamed}
        echo "✔️ Wrote {output.renamed}"
        """

rule cat_ass_16S:
    input:
        expand(os.path.join(RESULTS_DIR, "{sample}/{sample}_16S_seqs.fa"), sample=SAMPLES)
    output:
        cat_fa=os.path.join(RESULTS_DIR, "concat_16S/concat_ass_16S.fa")
    log:
        os.path.join(RESULTS_DIR, "logs/concat/cat_ass_16S.log")        
    message:
        "Concatenating the 16S sequences from all assemblies"
    shell:
        "(date && cat {input} > {output} && "
        "date) &> >(tee {log})"

rule install_hyperex:
    output:
        touch("submodules/hyperex_installed.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/install/hyperex_install.log")
    conda:
        "rust"
#       os.path.join(ENV_DIR, "rust.yaml")
    shell:
        """
        mkdir -p submodules
        if [ ! -d submodules/hyperex ]; then
            git clone https://github.com/Ebedthan/hyperex.git submodules/hyperex
        fi
        
        cd submodules/hyperex
        
        cargo build --release
        cargo test
        
        # Install hyperex into the submodules directory (local cargo install --root)
        cargo install --path . --root ../hyperex_install
        
        cd ../..
        
        # Marker file to indicate successful installation
        touch {output}
        """

rule hyperex_16S:
    input:
        rules.cat_ass_16S.output.cat_fa
    output:
        ext_fa=os.path.join(RESULTS_DIR, "hyperex/extracted_16S.fa")
    conda:
        "rust"
    log:
        os.path.join(RESULTS_DIR, "logs/hyperex/extracting_16S.log")
    params:
        cargo_path=config["hyperex"]["cargo_path"],
        path=config["hyperex"]["path"],
        fwd=config["hyperex"]["fwd"],
        rev=config["hyperex"]["rev"]
    threads:
        16
    message:
        "hyperex run to trim the 16S sequences from the concatenated assemblies"
    shell:
        "(date && export PATH={params.path}:$PATH && "
        "export PATH={params.cargo_path}:$PATH && "
        "hyperex -p $(basename -s '.fa' {output.ext_fa}) -f {params.fwd} -r {params.rev} {input} && "
        "date) &> >(tee {log})"

rule extract_16S:
    input:
        rules.cat_ass_16S.output.cat_fa
    output:
        ext_fa=os.path.join(RESULTS_DIR, "extracted/extracted_16S.fa"),
        map_fa=os.path.join(RESULTS_DIR, "extracted/extracted_map.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/extraction/perl_extract_16S.log")
    params:
        src=os.path.join(SRC_DIR, "modified_in_silico_pcr.pl"),
        fwd=config["hyperex"]["fwd"],
        rev=config["hyperex"]["rev"]
    message:
        "Extracting 16S using an in-silico PCR perl script"
    shell:
        "(date && "
        "perl {params.src} -s {input} -a {params.fwd} -b {params.rev} -e -f {output.ext_fa} -n {output.map_fa} && "
        "date) &> >(tee {log})"


