#!/bin/bash
# slurm/submit_smoothing_compile.sh — SBATCH script for the smoothing compile
# (reduction) step: runs 08c_smoothing_compile.R once over all partitions.
#
# Submitted by analysis/bash/stages/08_smoothing.sh (typically with
# --dependency=afterok:<arrayid>). May also be submitted directly once the
# runner array has finished:
#   sbatch analysis/bash/slurm/submit_smoothing_compile.sh <config.yml>
#
# Arg: $1 = config path.
#
#SBATCH --output=logs/%x_%j.log
# Partition: the cluster profile's BIGMEM class (lib/cluster.sh); list every
# partition that can fit the job: SLURM starts it on the first with free
# resources. Each must allow the --time below.
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH -c 1

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="$1"
if [[ -z "$CONFIG" ]]; then
  echo "Error: usage: slurm/submit_smoothing_compile.sh <config.yml> [--no_isweight]" >&2
  exit 1
fi
if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG" >&2; exit 1; fi

shift          # remaining args (e.g. --no_isweight) forwarded to R script
COMPILE_ARGS="$*"

echo "Config       : ${CONFIG}"
echo "Compile args : ${COMPILE_ARGS:-<none>}"

Rscript analysis/08c_smoothing_compile.R -C "${CONFIG}" ${COMPILE_ARGS} || exit 1

echo "***** End of script *****"
date
