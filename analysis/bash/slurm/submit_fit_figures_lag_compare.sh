#!/bin/bash
# slurm/submit_fit_figures_lag_compare.sh — SBATCH script for ONE
# cross-lag figures array task (06_fit_figures_lag_compare.R).
#
# Each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a "<batch_dir>::<group>"
# entry in the ';'-joined job list (same convention as slurm/submit_loo_compare.sh),
# then draws the lag_compare_erf / lag_compare_erf_country / lag_compare_elpd
# figures for that (track, group). Reads cb_surfaces (written by
# 05_fit_postprocess_cb.R) and the loo_compare_lags ranking (written by
# 05_fit_postprocess_loo_lags.R under generated_data/loo_compare/<batch>/,
# located via utils/loo_lags.R) — both must already exist for every lag in
# the group, or the corresponding panel is silently skipped (see script's own
# "surfaces missing" / "No loo_compare_lags ranking found" messages).
#
# Args: $1 = jobs string ("<dir>::<grp>;<dir>::<grp>;..."), $2 = extra flags
#       forwarded to the R script (e.g. "--format pdf --dpi 300"), optional.
# Submitted by analysis/bash/stages/06_fit_figures_lag_compare.sh. Do not invoke directly.
#
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH -c 4

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

# Disable BLAS thread proliferation — brms/dlnm manage their own threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

JOBS_STR="$1"
EXTRA="${2:-}"

if [[ -z "$JOBS_STR" ]]; then
  echo "Error: usage: slurm/submit_fit_figures_lag_compare.sh <jobs_str> [extra_flags]" >&2
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
echo "Extra     : ${EXTRA:-<none>}"

Rscript analysis/06_fit_figures_lag_compare.R -d "${BATCH_DIR}" --group "${GROUP}" ${EXTRA} || exit 1

echo "***** End of script *****"
date
