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

## 🧬 MAG annotations (`annotations` branch)

Steps for the dereplicated MQ MAG set (`MAG_collection/MQ_MAGs_all_0.95/dereplicated_genomes`).
Turn them on via `steps:` in `config/config.yaml`; paths are under `annotations:`.

| step | rules | output (`annotations.results_dir`) |
|---|---|---|
| `coverage_check` | `coverm_check.smk`: collates `MAG_collection/Coverm_output/*_output_coverm.tsv` and checks each holds exactly the dereplicated set. Also checks rel. abundance + unmapped ≈ 100 %, flags high-unmapped samples and undetected MAGs, and matches samples to the sequencing sheet (`dna_tube_num`, else `ceh_sample_code`; repairs the fused `dna_plate_numdna_plate_col` header) | `coverm/`: report, long table, MAG × sample matrices, covered-fraction-masked rel. abundance, summaries |
| `traits` | `microtrait.smk`: microTrait per MAG, then `make.genomeset.results` | `microtrait/genomeset_results.rds`, `microtrait/tables/*.tsv` |
| `amr` | `amr.smk`: RGI 6.0.3 (CARD, contigs) and DIAMOND vs BacMet2 EXP (≥80 % id/qcov, e ≤ 1e-5). Reuses `prodigal_annotations/{mag}.faa` if present, else calls genes with Prodigal | `amr/rgi_all.tsv`, `amr/bacmet_all.tsv` |

CARD (`rgi/`, incl. `localDB`) and BacMet (`bacmet/`) are downloaded once into
`annotations.dbs_dir` by local rules, as is the microTrait install. Run Snakemake from a
node with internet access.

```bash
# steps: ["coverage_check", "traits", "amr"] in config/config.yaml
snakemake --profile workflow/profiles/slurm -n
snakemake --profile workflow/profiles/slurm --jobs 100 \
  --conda-prefix /ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/MAG_collection/MQ_MAGs_all_0.95/annotations/susbus/envs
```

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

