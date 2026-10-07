#!/bin/bash -l

##############################
# SLURM
# NOTE: resources for the Snakemake *controller* only -- rule jobs are submitted
#       by Snakemake via workflow/profiles/slurm
#
# Usage (from the repo root):
#   bash sbatch_snakemake.sh setup                      # on the SOFTWARE node (internet): envs, images, CARD/BacMet downloads, microtrait
#   sbatch sbatch_snakemake.sh dryrun                   # dry-run, see the .slurm.out log
#   sbatch sbatch_snakemake.sh run                      # full run
#   sbatch sbatch_snakemake.sh run --until coverm_check # extra args are passed to snakemake

#SBATCH -J skyline_smk
#SBATCH -p ei-long             # USER_INPUT: no time limit on ei-long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o %x_%j.slurm.out     # matches *slurm*.out in .gitignore

set -eo pipefail

##############################
# SETTINGS

MODE="${1:-run}"; shift || true
REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

# working/annotation folder: envs, caches and temp files live here, never $HOME or /tmp
WORK="/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/MAG_collection/MQ_MAGs_all_0.95/annotations/susbus" # USER_INPUT
DBS="/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/data/Databases" # USER_INPUT (= annotations.dbs_dir)
PREFIX="${WORK}/envs"
SMK_ENV="snakemake"            # USER_INPUT: conda env with snakemake 7
SMK_JOBS="${SMK_JOBS:-100}"    # max. concurrent cluster jobs

mkdir -p "${WORK}/tmp" "${WORK}/conda_pkgs" "${WORK}/singularity_cache" "${PREFIX}"
export TMPDIR="${WORK}/tmp"
export CONDA_PKGS_DIRS="${WORK}/conda_pkgs"
export SINGULARITY_CACHEDIR="${WORK}/singularity_cache"
export APPTAINER_CACHEDIR="${WORK}/singularity_cache"
export MPLCONFIGDIR="${WORK}/tmp"   # matplotlib cache (RGI) off read-only $HOME

##############################
# ENV

eval "$(conda shell.bash hook)"
conda activate "${SMK_ENV}"

cd "${REPO_DIR}"
echo "[$(date)] host=$(hostname) mode=${MODE} repo=${REPO_DIR}"
snakemake --version

# conda envs + singularity images (RGI) both kept under ${WORK}
# NOTE (EI): images are pulled on the software node (internet) and run on compute nodes,
#            where /ei is visible inside containers (not on the software node)
SMK_COMMON=(--use-conda --conda-frontend conda --conda-prefix "${PREFIX}"
            --use-singularity --singularity-prefix "${WORK}/singularity_cache")
SMK=(snakemake --profile workflow/profiles/slurm "${SMK_COMMON[@]}")

##############################
# RUN

case "${MODE}" in
  setup)
    # porTraits setup targets only if the step is enabled in config.yaml
    PT_SETUP=()
    if grep -q '^steps:.*"portraits"' config/config.yaml; then
      PT_SETUP=("${WORK}/tools/porTraits_2804b21/main.nf"
                "${WORK}/singularity_cache/portraits/portraits.sif"
                "${WORK}/singularity_cache/portraits/recognise.sif"
                "${WORK}/singularity_cache/portraits/eggnog.sif"
                "${DBS}/pfam/Pfam31.0/Pfam-A.clans.tsv.gz")
    fi
    # run on the software node (needs internet); everything runs locally, no SLURM submission.
    # Step 1 also pulls the container images; CARD is only downloaded here and loaded
    # (rgi load, inside the container) as the first cluster job of the run.
    SMK_LOCAL=(snakemake -s workflow/Snakefile --configfile config/config.yaml "${SMK_COMMON[@]}" \
               --cores "${SMK_CORES:-4}" --rerun-incomplete -rp)
    "${SMK_LOCAL[@]}" --conda-create-envs-only "$@"
    "${SMK_LOCAL[@]}" "$@" \
      "${DBS}/rgi/card.json" \
      status/microtrait_installed.done \
      "${PT_SETUP[@]}" \
      "${DBS}/bacmet/bacmet_exp.dmnd"
    ;;
  dryrun)
    "${SMK[@]}" -n "$@"
    ;;
  run)
    "${SMK[@]}" --jobs "${SMK_JOBS}" "$@"
    ;;
  unlock)
    "${SMK[@]}" --unlock
    ;;
  *)
    echo "Unknown mode '${MODE}' (use: setup | dryrun | run | unlock)"; exit 1
    ;;
esac

echo "[$(date)] done"
