#!/bin/bash
# Run DLNM crossbasis postprocessing + figures locally on a config batch.
# Companion to stages/04_fit.sh (which fits + does plain postprocess on
# HPC) — use this afterwards to fill in the cb_surfaces + cb figures that the
# HPC batch doesn't produce.
#
# For each config in the batch (skipping manifest.yml), runs:
#   05_fit_postprocess_cb.R  (unless --skip-postprocess)
#   06_fit_figures_cb.R      (unless --skip-figures)
#
# Usage:
#   bash analysis/bash/tools/run_batch_cb_local.sh --batch <batch_name>
#   bash analysis/bash/tools/run_batch_cb_local.sh --batch <batch_name> --skip-postprocess
#   bash analysis/bash/tools/run_batch_cb_local.sh --batch <batch_name> --skip-figures
#   bash analysis/bash/tools/run_batch_cb_local.sh --batch <batch_name> --cumul --format png
#
# Example:
#   bash analysis/bash/tools/run_batch_cb_local.sh --batch prod_run_occ_grouped_enso_singlelag

set -euo pipefail

BATCH=""
SKIP_POSTPROCESS=FALSE
SKIP_FIGURES=FALSE
CUMUL=FALSE
FORMAT="pdf"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch|-b)          BATCH="$2";          shift 2 ;;
    --skip-postprocess)  SKIP_POSTPROCESS=TRUE; shift  ;;
    --skip-figures)      SKIP_FIGURES=TRUE;     shift  ;;
    --cumul)             CUMUL=TRUE;            shift  ;;
    --format)            FORMAT="$2";          shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --batch <batch_name> [--skip-postprocess] [--skip-figures] [--cumul] [--format pdf|png]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$BATCH" ]]; then
  echo "Error: --batch is required." >&2
  echo "Usage: $0 --batch <batch_name> [--skip-postprocess] [--skip-figures] [--cumul] [--format pdf|png]" >&2
  exit 1
fi

CONFIGDIR="analysis/configs/${BATCH}"

if [[ ! -d "$CONFIGDIR" ]]; then
  echo "Error: config directory '$CONFIGDIR' not found" >&2
  exit 1
fi

CONFIGNAMES=($(ls "$CONFIGDIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort))
N=${#CONFIGNAMES[@]}

if [[ "$N" -eq 0 ]]; then
  echo "Error: no config files found in '$CONFIGDIR'" >&2
  exit 1
fi

CUMUL_FLAG=""
if [[ "$CUMUL" == "TRUE" ]]; then
  CUMUL_FLAG="--cumul"
fi

echo "Batch '${BATCH}': ${N} config(s), skip_postprocess=${SKIP_POSTPROCESS}, skip_figures=${SKIP_FIGURES}, cumul=${CUMUL}, format=${FORMAT}"

for CFG in "${CONFIGNAMES[@]}"; do
  THISCONFIG="${CONFIGDIR}/${CFG}"
  echo "=== ${THISCONFIG} ==="

  if [[ "${SKIP_POSTPROCESS}" != "TRUE" ]]; then
    Rscript analysis/05_fit_postprocess_cb.R -C "${THISCONFIG}" ${CUMUL_FLAG} || exit 1
  fi

  if [[ "${SKIP_FIGURES}" != "TRUE" ]]; then
    Rscript analysis/06_fit_figures_cb.R -C "${THISCONFIG}" --format "${FORMAT}" || exit 1
  fi
done

echo "***** Done: ${N} config(s) processed *****"
