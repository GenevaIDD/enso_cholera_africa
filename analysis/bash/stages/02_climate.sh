#!/bin/bash
# stages/02_climate.sh — submit 02_preprocess_climate.R as one Slurm job.
#
# --stage all (default) runs the centroid match, the per-gid extraction and the
# database write in one job. With a new outbreak vintage only the match and the
# gids new to the caches cost anything; the per-gid caches
# (generated_data/climate_stats/, copied to the cluster with the inputs)
# carry the rest. The parquet corpus must be on the cluster
# (tools/download_jahn2026_parquet.sh, login node).
#
# Usage:
#   bash analysis/bash/stages/02_climate.sh [--stage all|match|extract|db] [--indir DIR]
#       [--after JOBID[:JOBID...]] [--cpus 16] [--mem 64G] [--time 08:00:00] [--plan] [-- EXTRA R ARGS]
#
# Arguments after `--` go to the R script unchanged (e.g. --force TRUE,
# --stats_db <staged.sqlite>). -s TRUE (rebuild climate_stats) is always passed.
#
# Prints CLIMATE_JOB_ID=.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

STAGE="all"
INDIR=""
AFTER=""
CPUS=16
MEM="64G"
TIME="08:00:00"
PLAN=false
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage) STAGE="$2"; shift 2 ;;
    --indir) INDIR="$2"; shift 2 ;;
    --after) AFTER="$2"; shift 2 ;;
    --cpus)  CPUS="$2";  shift 2 ;;
    --mem)   MEM="$2";   shift 2 ;;
    --time)  TIME="$2";  shift 2 ;;
    --plan)  PLAN=true;  shift   ;;
    --)      shift; EXTRA=("$@"); break ;;
    *) echo "Unknown argument: $1" >&2
       echo "Usage: $0 [--stage all|match|extract|db] [--indir DIR] [--after JOBID] [--cpus N] [--mem M] [--time T] [--plan] [-- R args]" >&2
       exit 1 ;;
  esac
done
case "$STAGE" in all|match|extract|db) ;; *) echo "Error: --stage must be all, match, extract or db." >&2; exit 1 ;; esac

ARGS=(analysis/02_preprocess_climate.R --stage "$STAGE" -s TRUE)
[[ -n "$INDIR" ]] && ARGS+=(--indir "$INDIR")
ARGS+=("${EXTRA[@]}")
DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")
cholera_sbatch_cmd LONG --parsable --job-name="prep_02_climate_${STAGE}" -c "$CPUS" --mem="$MEM" --time="$TIME" \
  ${DEP[@]+"${DEP[@]}"} analysis/bash/slurm/submit_preprocess.sh "${ARGS[@]}"
CMD=("${CHOLERA_SBATCH_CMD[@]}")

if [[ "$PLAN" == "true" ]]; then
  echo "PLAN: ${CMD[*]}" >&2
  echo "CLIMATE_JOB_ID=PLAN_02_climate"
  exit 0
fi

mkdir -p "$CHOLERA_LOG_DIR"
ID=$("${CMD[@]}")
echo "CLIMATE_JOB_ID=${ID}"
