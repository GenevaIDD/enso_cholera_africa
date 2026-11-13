#!/bin/bash
# slurm/submit_kl_array.sh — one brms_kl occurrence fit per array slot
# (04_fit_occ_brms_kl.R via the 04_fit_run.R dispatcher), over the configs of
# one batch dir. Submitted by stages/04_fit_kl.sh.
#
# Usage (normally via the stage driver):
#   sbatch --array=0-23 analysis/bash/slurm/submit_kl_array.sh <batch_dir>
#
# Slot i fits the i-th config of <batch_dir> in sorted order (manifest.yml
# excluded), the convention of slurm/submit_fit.sh. Resources: -c 32 covers
# 6 chains x 5 threads. A config whose fit
# already exists is skipped by the worker (overwrite: no).
# Memory: 24G = ~4x the measured peak (4.5-5.5 GB per fit, 8 chains x 4 threads,
# country-only ENSO slopes).
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=03:00:00
#SBATCH --mem=24G
#SBATCH -c 32

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

# Disable BLAS thread proliferation — Stan runs one process per chain.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIGDIR="$1"
if [[ -z "$CONFIGDIR" || ! -d "$CONFIGDIR" ]]; then
  echo "Error: usage: sbatch --array=... slurm/submit_kl_array.sh <batch_dir>" >&2
  exit 1
fi

CONFIGNAMES=($(ls "$CONFIGDIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort))
THISCONFIG="${CONFIGDIR}/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}"

echo "Batch dir : ${CONFIGDIR}"
echo "Config    : ${THISCONFIG}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"
echo "CPUs      : ${SLURM_CPUS_PER_TASK:-<unset>}"

Rscript analysis/04_fit_run.R -C "${THISCONFIG}" || exit 1

echo "***** End of script *****"
date
