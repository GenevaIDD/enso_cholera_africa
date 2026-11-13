#!/bin/bash
# slurm/submit_loo_compare.sh — SBATCH script for ONE LOO-comparison array task.
#
# Called by SLURM; each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a
# "<batch_dir>::<group>" entry in the ';'-joined job list, then invokes
# 05_fit_postprocess_loo_lags.R for that (track, group). Each task loads every single-lag
# fit in the group, computes/caches per-fit brms::loo, and writes the ranking RDS
# + elpd-by-lag PNG under generated_data/loo_compare/<batch>/ (and, with
# --lag_sidecar, the per-group lag_loo_sidecar_*.rds there too).
#
# Args: $1 = jobs string ("<dir>::<grp>;<dir>::<grp>;..."), $2 = redo flag ("--redo" or ""),
#       $3 = optional extra 05_fit_postprocess_loo_lags.R arguments, word-split
#            (e.g. "--lag_sidecar --sidecar_only"). The wall time is set by the
#            stage driver's --time (sbatch --time overrides the header below).
# Submitted by analysis/bash/stages/05_loo_compare.sh. Do not invoke directly.
#
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=2:00:00
# Memory: dominated by holding one brms_cb fit in RAM at a time plus its loo
# object; the lags are loaded serially, not all at once. 32G leaves comfortable
# margin for the largest grouped occ/size fits. (Measured locally:
# --lag_sidecar over 23 cached-LOO single-lag occurrence fits of ~11 MB peaked
# at 0.6 GB RSS and took 18 s.)
#SBATCH --mem=32G
#SBATCH -c 4

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

# Disable BLAS thread proliferation — loo/brms manage their own threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

JOBS_STR="$1"
REDO="${2:-}"
EXTRA="${3:-}"

if [[ -z "$JOBS_STR" ]]; then
  echo "Error: usage: slurm/submit_loo_compare.sh <jobs_str> [--redo]" >&2
  exit 1
fi

# Split the ';'-joined job list and select this task's entry.
IFS=';' read -r -a JOBS <<< "$JOBS_STR"
ENTRY="${JOBS[$SLURM_ARRAY_TASK_ID]}"
if [[ -z "$ENTRY" ]]; then
  echo "Error: no job entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} (have ${#JOBS[@]} job(s))." >&2
  exit 1
fi

BATCH_DIR="${ENTRY%%::*}"
GROUP="${ENTRY##*::}"

if [[ ! -d "$BATCH_DIR" ]]; then echo "Error: batch dir not found: $BATCH_DIR" >&2; exit 1; fi
if [[ -z "$GROUP" ]]; then echo "Error: empty group parsed from '$ENTRY'." >&2; exit 1; fi

echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Batch dir : ${BATCH_DIR}"
echo "Group     : ${GROUP}"
echo "Redo      : ${REDO:-<off>}"
echo "Extra     : ${EXTRA:-<none>}"

# ${EXTRA} is deliberately unquoted: it carries several flags.
Rscript analysis/05_fit_postprocess_loo_lags.R -d "${BATCH_DIR}" --group "${GROUP}" ${REDO} ${EXTRA} || exit 1

echo "***** End of script *****"
date
