#!/bin/bash
# slurm/submit_postprocess_cb.sh — SBATCH script for ONE DLNM crossbasis
# postprocess array task (05_fit_postprocess_cb.R).
#
# Each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to the Nth *.yml in
# the batch dir (same sort/index convention as slurm/submit_fit.sh), loads
# that config's brms_cb fit + cb_meta sidecar, and writes the cb_surfaces
# artefact via dlnm::crosspred() per posterior draw.
#
# Args: $1 = batch dir, $2 = extra flags forwarded verbatim to the R script
#       (e.g. "--cumul --n_at 100"), optional.
# Submitted by analysis/bash/stages/05_postprocess_cb.sh. Do not invoke directly.
#
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=1:00:00
#SBATCH --mem=16G
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

BATCH_DIR="$1"
EXTRA="${2:-}"

if [[ -z "$BATCH_DIR" ]]; then
  echo "Error: usage: slurm/submit_postprocess_cb.sh <batch_dir> [extra_flags]" >&2
  exit 1
fi
if [[ ! -d "$BATCH_DIR" ]]; then
  echo "Error: batch dir not found: $BATCH_DIR" >&2
  exit 1
fi

CONFIGNAMES=($(ls "$BATCH_DIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort))
THISCONFIG="${BATCH_DIR%/}/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]:-}" ]]; then
  echo "Error: no config at index ${SLURM_ARRAY_TASK_ID} (have ${#CONFIGNAMES[@]})." >&2
  exit 1
fi

echo "Batch dir : ${BATCH_DIR}"
echo "Config    : ${THISCONFIG}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Extra     : ${EXTRA:-<none>}"

Rscript analysis/05_fit_postprocess_cb.R -C "${THISCONFIG}" ${EXTRA} || exit 1

echo "***** End of script *****"
date
