#!/bin/bash
# stages/02_crosswalk.sh — submit 02_preprocess_crosswalk.R as one Slurm job.
#
# Reads 01's outbreak_shapefiles.gpkg (vintage folder) and 02c's compiled
# population crosswalk; writes the LP -> gid crosswalk and u_gids_crosswalk.rds
# (the gid universe of 02_preprocess_climate.R --gid_source crosswalk).
#
# Usage:
#   bash analysis/bash/stages/02_crosswalk.sh --config <config.yml> [--indir DIR] [--after JOBID[:JOBID...]] [--plan]
#
# --after takes a colon-joined id list (afterok on all): campaigns/preprocess.sh
# passes 01's outbreaks job and 02c's compile job.
#
# Prints CROSSWALK_JOB_ID=.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

CONFIG=""
INDIR=""
AFTER=""
PLAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C) CONFIG="$2"; shift 2 ;;
    --indir)     INDIR="$2";  shift 2 ;;
    --after)     AFTER="$2";  shift 2 ;;
    --plan)      PLAN=true;   shift   ;;
    *) echo "Unknown argument: $1" >&2
       echo "Usage: $0 --config <config.yml> [--indir DIR] [--after JOBID[:JOBID...]] [--plan]" >&2
       exit 1 ;;
  esac
done
if [[ -z "$CONFIG" ]]; then
  echo "Error: --config is required (the population block names the crosswalk inputs)." >&2
  exit 1
fi

ARGS=(analysis/02_preprocess_crosswalk.R -C "$CONFIG")
[[ -n "$INDIR" ]] && ARGS+=(--indir "$INDIR")
DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")
# Sizing: CHOLERA_RES_CROSSWALK in the cluster profile overrides these defaults.
cholera_sbatch_cmd LONG --res CROSSWALK --defaults "--cpus-per-task=1 --mem=16G --time=02:00:00" \
  --parsable --job-name=prep_02_crosswalk \
  ${DEP[@]+"${DEP[@]}"} analysis/bash/slurm/submit_preprocess.sh "${ARGS[@]}"
CMD=("${CHOLERA_SBATCH_CMD[@]}")

if [[ "$PLAN" == "true" ]]; then
  echo "PLAN: ${CMD[*]}" >&2
  echo "CROSSWALK_JOB_ID=PLAN_02_crosswalk"
  exit 0
fi

mkdir -p "$CHOLERA_LOG_DIR"
ID=$("${CMD[@]}")
echo "CROSSWALK_JOB_ID=${ID}"
