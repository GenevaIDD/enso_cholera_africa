#!/bin/bash
# slurm/submit_kl_compile.sh — compile every Stan program of a brms_kl
# batch once, before its fit array (utils/kl_compile_batch.R: one
# `04_fit_run.R --compile_only` per config, in parallel; configs sharing a
# program meet at its compile lock, the first compiles, the others reuse it).
# Submitted by stages/04_fit_kl.sh; the fit array waits on it (afterok), so no
# array task compiles or execs an executable another task is still writing.
# Same environment as slurm/submit_kl_array.sh (cholera_env), so the executable matches.
#
# Usage:
#   sbatch analysis/bash/slurm/submit_kl_compile.sh <batch_dir>
#
# Memory: up to one R session per config (data load + brms::stancode) at once;
# measured 17.7 GB peak for 6 configs (~3 GB each), so 40G covers the 12
# sessions -c 12 runs at once.
#SBATCH --job-name=kl_compile
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=40G
#SBATCH -c 12

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIGDIR="$1"
if [[ -z "$CONFIGDIR" || ! -d "$CONFIGDIR" ]]; then
  echo "Error: usage: sbatch slurm/submit_kl_compile.sh <batch_dir>" >&2
  exit 1
fi

Rscript analysis/utils/kl_compile_batch.R -d "${CONFIGDIR}" || exit 1

echo "***** End of script *****"
date
