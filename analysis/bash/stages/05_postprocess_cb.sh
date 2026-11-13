#!/bin/bash
# stages/05_postprocess_cb.sh — submit the DLNM crossbasis postprocess
# step (05_fit_postprocess_cb.R) as a SLURM ARRAY job.
#
# Sibling to stages/04_fit.sh: that script only runs 05_fit_postprocess.R
# (RE stats / coef summaries / loo) per config. For brms_cb backends (DLNM /
# ENSO crossbasis sweeps such as the *_grouped_enso_singlelag batches), the
# cb_surfaces artefact consumed by 06_fit_figures_lag_compare.R (and
# 06_fit_figures_cb*.R) is written by 05_fit_postprocess_cb.R instead, and is
# NOT produced by the plain batch array — run this afterward.
#
# Usage:
#   bash analysis/bash/stages/05_postprocess_cb.sh --batch-dir <dir>
#   bash analysis/bash/stages/05_postprocess_cb.sh --batch-dir <dir> --concurrency 10
#   bash analysis/bash/stages/05_postprocess_cb.sh --batch-dir <dir> --extra "--cumul --n_at 100"
#   bash analysis/bash/stages/05_postprocess_cb.sh --batch-dir <dir> --list
#
# Prints the submitted array job id on stdout as "CB_JOB_ID=<id>" (last line)
# so callers (e.g. an orchestrator script) can capture it for
# --dependency=afterok chaining.
#
# Options:
#   --batch-dir DIR     Directory of *.yml configs (required; manifest.yml excluded)
#   --concurrency K     Cap simultaneous array slots   [default: all]
#   --partition P       SLURM partition                [default: the cluster profile's,
#                       CHOLERA_RES_POSTPROCESS_CB, else its SHORT class]
#   --time T            SLURM time limit                [default: 1:00:00]
#   --extra "FLAGS"     Extra flags forwarded to 05_fit_postprocess_cb.R (e.g. "--cumul --n_at 100")
#   --list              Print the config list and exit (no submit)

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"

BATCH_DIR=""
CONCURRENCY=""
PARTITION=""
TIME_LIMIT="1:00:00"
EXTRA=""
LIST_ONLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch-dir|-d)   BATCH_DIR="$2";    shift 2 ;;
    --concurrency)    CONCURRENCY="$2";  shift 2 ;;
    --partition)      PARTITION="$2";    shift 2 ;;
    --time)           TIME_LIMIT="$2";   shift 2 ;;
    --extra)          EXTRA="$2";        shift 2 ;;
    --list)           LIST_ONLY=true;    shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --batch-dir <dir> [--concurrency K] [--partition P] [--time T] [--extra \"FLAGS\"] [--list]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$BATCH_DIR" ]]; then
  echo "Error: --batch-dir is required." >&2
  exit 1
fi
if [[ ! -d "$BATCH_DIR" ]]; then
  echo "Error: batch dir not found: $BATCH_DIR" >&2
  exit 1
fi

mapfile -t CONFIGNAMES < <(ls "$BATCH_DIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort)
N=${#CONFIGNAMES[@]}

if [[ "$N" -eq 0 ]]; then
  echo "Error: no config files found in '$BATCH_DIR'" >&2
  exit 1
fi

echo "=== CB postprocess: ${N} config(s) in ${BATCH_DIR} ===" >&2
for i in "${!CONFIGNAMES[@]}"; do
  printf '  [%d] %s\n' "$i" "${CONFIGNAMES[$i]}" >&2
done

if [[ "$LIST_ONLY" == "true" ]]; then
  echo "(--list: not submitting)" >&2
  exit 0
fi

ARRAY_SPEC="0-$((N - 1))"
[[ -n "$CONCURRENCY" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${CONCURRENCY}"

mkdir -p "$CHOLERA_LOG_DIR"

BATCH_NAME="$(basename "${BATCH_DIR%/}")"
echo "--- Submitting array ${ARRAY_SPEC} as job 'postprocess_cb_${BATCH_NAME}' ..." >&2
CB_JOB_ID=$(cholera_sbatch SHORT --parsable \
  --array="${ARRAY_SPEC}" \
  --job-name="postprocess_cb_${BATCH_NAME}" \
  ${PARTITION:+--partition="${PARTITION}"} \
  --time="${TIME_LIMIT}" \
  analysis/bash/slurm/submit_postprocess_cb.sh \
  "$BATCH_DIR" "$EXTRA")

echo "CB postprocess array job id: ${CB_JOB_ID}" >&2
echo "CB_JOB_ID=${CB_JOB_ID}"
