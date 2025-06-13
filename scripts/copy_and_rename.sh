#!/bin/bash
# Script to copy 16S_seqs.fa files from HiFi assemblies and rename headers to include the sample name
# Usage: ./copy_and_rename.sh
# Absolute path to SSA directory with original 16S files
base_path="/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/SSA"

# Output will be created in your current working directory
output_base="./results"

# List of sample folder names
samples=(
    b10_d3_0200_con
    b10_d3_1400_con
    b20_d3_0200_bio
    b26_t3_con
    r15_d3_0200_con
    r15_d3_1400_con
    r30_d3_0200_bio
    r30_d3_1400_bio
)

# Make results folder
mkdir -p "$output_base"

# Loop through sample names
for sample in "${samples[@]}"; do
    input_file="$base_path/$sample/metamdbg/annotation/16S_seqs.fa"
    sample_output_dir="$output_base/$sample"
    output_file="${sample_output_dir}/${sample}_16S_seqs.fa"

    mkdir -p "$sample_output_dir"

    if [[ -f "$input_file" ]]; then
        awk -v prefix="$sample" '
            /^>/ {sub(/^>/, ">" prefix "-"); print; next}
            {print}
        ' "$input_file" > "$output_file"
        echo "✔️ Wrote $output_file"
    else
        echo "⚠️ Missing: $input_file"
    fi
done
