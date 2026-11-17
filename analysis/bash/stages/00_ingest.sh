#!/bin/bash
# stages/00_ingest.sh — submit 00_ingest_outbreakextractr.R as one Slurm job.
#
# Reads the OutbreakExtractR outputs of the outbreak vintage
# (paths.R::default_outbreak_data_dir(), or --indir) and writes the legacy RDS
# files into the same folder. The vintage folder must be on the cluster first.
#
# Usage:
#   bash analysis/bash/stages/00_ingest.sh [--indir DIR] [--after JOBID] [--cpus 8] [--mem 48G] [--time 04:00:00] [--plan]
#
# --after  submit with --dependency=afterok:JOBID (campaigns/preprocess.sh)
# --plan   print the sbatch line, submit nothing
#
# Prints INGEST_JOB_ID=<id>.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

INDIR=""
AFTER=""
CPUS=8          # 00 caps its furrr workers at 8
MEM="48G"       # 8.8 GB resident in the parent locally, plus the workers' copies
TIME="04:00:00"
PLAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --indir) INDIR="$2"; shift 2 ;;
    --after) AFTER="$2"; shift 2 ;;
    --cpus)  CPUS="$2";  shift 2 ;;
    --mem)   MEM="$2";   shift 2 ;;
    --time)  TIME="$2";  shift 2 ;;
    --plan)  PLAN=true;  shift   ;;
    *) echo "Unknown argument: $1" >&2
       echo "Usage: $0 [--indir DIR] [--after JOBID] [--cpus N] [--mem M] [--time T] [--plan]" >&2
       exit 1 ;;
  esac
done

ARGS=(analysis/00_ingest_outbreakextractr.R)
[[ -n "$INDIR" ]] && ARGS+=(--indir "$INDIR")
DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")
cholera_sbatch_cmd LONG --parsable --job-name=prep_00_ingest -c "$CPUS" --mem="$MEM" --time="$TIME" \
  ${DEP[@]+"${DEP[@]}"} analysis/bash/slurm/submit_preprocess.sh "${ARGS[@]}"
CMD=("${CHOLERA_SBATCH_CMD[@]}")

if [[ "$PLAN" == "true" ]]; then
  echo "PLAN: ${CMD[*]}" >&2
  echo "INGEST_JOB_ID=PLAN_00_ingest"
  exit 0
fi

mkdir -p "$CHOLERA_LOG_DIR"
ID=$("${CMD[@]}")
echo "INGEST_JOB_ID=${ID}"
