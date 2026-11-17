#!/bin/bash
# stages/06_fit_figures_lag_compare.sh — submit the cross-lag comparison
# figures (06_fit_figures_lag_compare.R) as a SLURM array, one task per
# (track, group).
#
# Needs, per lag config: cb_surfaces (05_fit_postprocess_cb.R) and the group's
# loo_compare_lags ranking (generated_data/loo_compare/<batch>/): run
# stages/05_postprocess_cb.sh and stages/05_loo_compare.sh first, or chain with
# --dependency afterok:<ids>.
#
# Usage:
#   bash analysis/bash/stages/06_fit_figures_lag_compare.sh --occ-dir <dir>
#   bash analysis/bash/stages/06_fit_figures_lag_compare.sh --track dur --dur-dir <dir>
#   bash analysis/bash/stages/06_fit_figures_lag_compare.sh --occ-dir <dir> --format pdf --dpi 300
#   bash analysis/bash/stages/06_fit_figures_lag_compare.sh --occ-dir <dir> --dependency afterok:12345 --list
#
# Options:
#   --track T          occ | size | dur | mon | both (occ+size) [default: occ]
#   --groups G1,...    ENSO groups [default: enso_east,enso_horn,enso_south,
#                      enso_west]
#   --occ-dir / --size-dir / --dur-dir / --mon-dir DIR   Batch config dirs
#   --concurrency K    Cap simultaneous array slots [default: all]
#   --format {png|pdf} Figure format [default: png]
#   --dpi N            Figure resolution [default: 150]
#   --dependency SPEC  Forwarded to sbatch --dependency
#   --list             Print the (track, group) jobs and exit
# Last stdout line: FIG_JOB_ID=<id>.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"

TRACK="occ"
GROUP_CSV="enso_east,enso_horn,enso_south,enso_west"
OCC_DIR=""
SIZE_DIR=""
DUR_DIR=""
MON_DIR=""
CONCURRENCY=""
FORMAT="png"
DPI="150"
DEPENDENCY=""
LIST_ONLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --track)        TRACK="$2";        shift 2 ;;
    --groups)       GROUP_CSV="$2";    shift 2 ;;
    --occ-dir)      OCC_DIR="$2";      shift 2 ;;
    --size-dir)     SIZE_DIR="$2";     shift 2 ;;
    --dur-dir)      DUR_DIR="$2";      shift 2 ;;
    --mon-dir)      MON_DIR="$2";      shift 2 ;;
    --concurrency)  CONCURRENCY="$2";  shift 2 ;;
    --format)       FORMAT="$2";       shift 2 ;;
    --dpi)          DPI="$2";          shift 2 ;;
    --dependency)   DEPENDENCY="$2";   shift 2 ;;
    --list)         LIST_ONLY=true;    shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--track {occ|size|dur|mon|both}] [--groups G1,G2,...] [--occ-dir DIR] [--size-dir DIR] [--dur-dir DIR] [--mon-dir DIR] [--concurrency K] [--format png|pdf] [--dpi N] [--dependency SPEC] [--list]" >&2
      exit 1 ;;
  esac
done

case "$TRACK" in
  occ|size|dur|mon|both) ;;
  *) echo "Error: --track must be one of {occ|size|dur|mon|both}, got '$TRACK'." >&2; exit 1 ;;
esac

declare -a TRACK_DIRS=()
if [[ "$TRACK" == "occ" || "$TRACK" == "both" ]]; then
  [[ -z "$OCC_DIR" ]] && { echo "Error: --occ-dir required for --track occ/both." >&2; exit 1; }
  TRACK_DIRS+=("$OCC_DIR")
fi
if [[ "$TRACK" == "size" || "$TRACK" == "both" ]]; then
  [[ -z "$SIZE_DIR" ]] && { echo "Error: --size-dir required for --track size/both." >&2; exit 1; }
  TRACK_DIRS+=("$SIZE_DIR")
fi
if [[ "$TRACK" == "dur" ]]; then
  [[ -z "$DUR_DIR" ]] && { echo "Error: --dur-dir required for --track dur." >&2; exit 1; }
  TRACK_DIRS+=("$DUR_DIR")
fi
if [[ "$TRACK" == "mon" ]]; then
  [[ -z "$MON_DIR" ]] && { echo "Error: --mon-dir required for --track mon." >&2; exit 1; }
  TRACK_DIRS+=("$MON_DIR")
fi

IFS=',' read -r -a GROUP_ARR <<< "$GROUP_CSV"

declare -a JOBS=()
for d in "${TRACK_DIRS[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: batch config dir not found: $d" >&2
    exit 1
  fi
  for g in "${GROUP_ARR[@]}"; do
    g_trimmed="${g//[[:space:]]/}"
    [[ -z "$g_trimmed" ]] && continue
    JOBS+=("${d}::${g_trimmed}")
  done
done

N_JOBS=${#JOBS[@]}
if [[ "$N_JOBS" -eq 0 ]]; then
  echo "Error: empty job list (check --track / --groups / --occ-dir / --size-dir / --dur-dir / --mon-dir)." >&2
  exit 1
fi

echo "=== Cross-lag figures: ${N_JOBS} array task(s) ===" >&2
for i in "${!JOBS[@]}"; do
  printf '  [%d] %s\n' "$i" "${JOBS[$i]}" >&2
done

if [[ "$LIST_ONLY" == "true" ]]; then
  echo "(--list: not submitting)" >&2
  exit 0
fi

JOBS_STR=$(IFS=';'; echo "${JOBS[*]}")
EXTRA="--format ${FORMAT} --dpi ${DPI}"

ARRAY_SPEC="0-$((N_JOBS - 1))"
[[ -n "$CONCURRENCY" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${CONCURRENCY}"

mkdir -p "$CHOLERA_LOG_DIR"

DEP_ARGS=()
[[ -n "$DEPENDENCY" ]] && DEP_ARGS=(--dependency="$DEPENDENCY")

echo "--- Submitting array ${ARRAY_SPEC} as job 'fit_figures_lag_compare' ${DEPENDENCY:+(dependency: $DEPENDENCY)} ..." >&2
FIG_JOB_ID=$(cholera_sbatch SHORT --parsable \
  --array="${ARRAY_SPEC}" \
  --job-name="fit_figures_lag_compare" \
  ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} \
  analysis/bash/slurm/submit_fit_figures_lag_compare.sh \
  "$JOBS_STR" "$EXTRA")

echo "Cross-lag figures array job id: ${FIG_JOB_ID}" >&2
echo "FIG_JOB_ID=${FIG_JOB_ID}"
