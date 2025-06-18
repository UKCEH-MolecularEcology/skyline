# skyline
Analysis for Skyline project in collaboration with Earlham Institute on the Decoding Biodiversity project

This repository contains a Snakemake-based workflow for analyzing **16S rRNA genes** from long-read metagenomic assemblies. The pipeline is designed for efficient, reproducible extraction, annotation, and sub-region (e.g., V4) analysis of 16S genes using long-read data.

---

## 🔬 Overview

The workflow supports:

- Assembly-based identification of full-length 16S rRNA genes
- FASTA header standardization and metadata linking
- Sub-region extraction (e.g., V4) using [`hyperex`](https://github.com/Ebedthan/hyperex)
- Integration with MAGs and downstream comparative analyses

---

## 🧩 Submodules

This project uses [`hyperex`](https://github.com/Ebedthan/hyperex) as a **git submodule** to extract V-regions from full-length 16S genes.

To clone with submodules:

```bash
git clone --recurse-submodules git@github.com:UKCEH-MolecularEcology/skyline.git
```

To update submodules later:

```bash
git submodule update --init --recursive
```

---

## 📁 Workflow Structure

Key rules in the Snakemake workflow include:

- `extract_16S_all`: Ensures 16S genes are extracted and renamed per sample
- `install_hyperex`: Clones and installs hyperex from source via Cargo
- `extract_v4`: Uses hyperex to extract V4 (or other) subregions from full-length 16S genes

---

## 🔧 Installation

### Dependencies

- Snakemake ≥ 7
- Python ≥ 3.8
- Rust + Cargo (for installing `hyperex`)
- Conda (for environment management)

### Setup Instructions

1. Clone this repo and its submodules:

```bash
git clone --recurse-submodules git@github.com:UKCEH-MolecularEcology/skyline.git
cd skyline
```

2. Create the Snakemake environment:

```bash
conda env create -f envs/snakemake.yaml
conda activate snakemake
```

3. Install Rust (if not already installed):

```bash
curl https://sh.rustup.rs -sSf | sh
source $HOME/.cargo/env
```

4. Run the workflow:

```bash
snakemake --use-conda --profile workflow/profiles/slurm
```

---

## 🧪 Example Inputs

Sample configuration includes:

- `samples.csv` - Sample table with paths to long-read assemblies and metadata
- `config.yaml` - Pipeline configuration

---

## 📊 Output

- Full-length 16S gene FASTA files per sample
- Renamed and standardized headers
- Extracted V4 regions (`*_v4.fa`)
- Final status indicators in `status/`

---

## 📎 Credits

- 16S extraction and assembly: Skyline team (UKCEH Molecular Ecology Group) and Quince group (Earlham Institute)
- V-region extraction tool: [`hyperex`](https://github.com/Ebedthan/hyperex)
- Maintainer: [@susheelbhanu](https://github.com/susheelbhanu)

---

## 📄 License

MIT License unless otherwise stated.

---

## 🧠 Citation

If you use this pipeline, please cite:

> Busi et al. (in prep). **16S rRNA gene region analysis from long-read metagenomes using the Skyline pipeline**.

