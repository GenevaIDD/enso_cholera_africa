#!/bin/bash
# slurm/submit_preprocess.sh — SBATCH payload for one 00-03 data-leg step.
#
# Generic single-job payload: the stage driver names the script and its
# arguments, and sets --job-name / -c / --mem / --time on the sbatch line
# (the header below is only the fallback).
#
# Args: $1 = analysis/NN*.R script, $2.. = its arguments.
#
# Submitted by analysis/bash/stages/{00_ingest,01_preprocess,02_population,
# 02_crosswalk,02_climate,03_prepare}.sh. Do not invoke directly.
#
# R 4.3.2 + GDAL: the 00-02 steps read and write sf geometry (R/4.2.1 has no sf
# here), as in slurm/submit_population.sh.
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH -c 1

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env GEO

echo "***** Beginning of script *****"
date

# Parallelism is furrr workers (future::availableCores() reads the allocation);
# keep BLAS single-threaded inside each of them.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ $# -lt 1 || ! -f "$1" ]]; then
  echo "Error: usage: slurm/submit_preprocess.sh <analysis/NN*.R> [args...]" >&2
  exit 1
fi

echo "Script : $1"
echo "Args   : ${*:2}"
echo "CPUs   : ${SLURM_CPUS_PER_TASK:-1}"

Rscript "$@" || exit 1

echo "***** End of script *****"
date
