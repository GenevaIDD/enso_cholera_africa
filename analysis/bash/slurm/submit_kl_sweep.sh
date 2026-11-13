#!/bin/bash
# slurm/submit_kl_sweep.sh — compare the fits of a brms_kl sweep
# (05_fit_postprocess_kl_sweep.R) once they are all on disk. Submitted by
# stages/04_fit_kl.sh --postprocess, afterok the fit array.
#
# Usage:
#   sbatch analysis/bash/slurm/submit_kl_sweep.sh <batch_dir> [05 options...]
#   e.g. ... <batch_dir> --mode_tol 1.5 --n_agree 4
#
# Memory: every fit (12-25 MB on disk, several times that in memory) is held at
# once for the cross-fit LOO; measured 0.3 GB peak for a 6-fit sweep, so 16G.
#SBATCH --job-name=kl_sweep_compare
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH -c 4

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

CONFIGDIR="$1"
if [[ -z "$CONFIGDIR" || ! -d "$CONFIGDIR" ]]; then
  echo "Error: usage: sbatch slurm/submit_kl_sweep.sh <batch_dir> [05 options...]" >&2
  exit 1
fi
shift

Rscript analysis/05_fit_postprocess_kl_sweep.R -d "${CONFIGDIR}" "$@" || exit 1

echo "***** End of script *****"
date
