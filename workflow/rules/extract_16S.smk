"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/kraken2.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on reads
"""

import os

localrules: install_hyperex

###################
# RULES
###################
rule extract_16S_all:
    input:
        expand(f"{RESULTS_DIR}/{{sample}}/{{sample}}_16S_seqs.fa", sample=SAMPLES),
        "submodules/hyperex_installed.txt"
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
