#!/bin/bash
# Usage: sbatch [--partition P] [--time T] slurm/submit_fit.sh <batch_name> [SKIP_POSTPROCESS] [PARTITION] [TIME_LIMIT]
#
# PARTITION/TIME_LIMIT positional args (3rd/4th) are informational only — by
# the time this script body runs, SLURM has already allocated resources per
# the #SBATCH pragmas below (or per any --partition/--time override passed to
# sbatch / stages/04_fit.sh at submission time). They default to this
# script's own #SBATCH defaults and are echoed for traceability.
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=3-23:00:00
#SBATCH --mem=20G
#SBATCH -c 16

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
DATE=$(date +%d-%b-%H:%M:%S)
echo $DATE

# Disable BLAS thread proliferation — brms/Stan use their own threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

BATCH="$1"
SKIP_POSTPROCESS="${2:-FALSE}"
PARTITION="${3:-<profile>}"
TIME_LIMIT="${4:-3-23:00:00}"

if [[ -z "$BATCH" ]]; then
  echo "Error: usage: slurm/submit_fit.sh <batch_name> [SKIP_POSTPROCESS] [PARTITION] [TIME_LIMIT]" >&2
  exit 1
fi

CONFIGDIR="analysis/configs/${BATCH}"

if [[ ! -d "$CONFIGDIR" ]]; then
  echo "Error: config directory not found: $CONFIGDIR" >&2
  exit 1
fi

CONFIGNAMES=($(ls "$CONFIGDIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort))
THISCONFIG="${CONFIGDIR}/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}"

echo "Batch     : ${BATCH}"
echo "Config    : ${THISCONFIG}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Partition : ${SLURM_JOB_PARTITION:-$PARTITION} (requested: ${PARTITION})"
echo "Time limit: ${TIME_LIMIT}"

Rscript analysis/04_fit_run.R -C "${THISCONFIG}" || exit 1

if [[ "${SKIP_POSTPROCESS}" != "TRUE" ]]; then
  Rscript analysis/05_fit_postprocess.R -C "${THISCONFIG}" || exit 1
fi

echo "***** End of script *****"
date
