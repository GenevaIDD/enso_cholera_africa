#!/bin/bash
# stages/01_preprocess.sh — submit the 01 outbreak steps as Slurm jobs:
#   01_preprocess_outbreaks.R (-f FALSE: no figures on the cluster)
#     -> afterok: 01b_preprocess_outbreaks_monthly.R
#     -> afterok: 01_preprocess_national_zero_months.R (the deploy's rule)
#
# 01b and 01_nz both read 01's outbreak metadata, so they wait for it. All
# three write under the outbreak-data root (paths.R::default_outbreak_data_dir()).
#
# The national-zero months follow the deploy YAML's deploy.episodes block
# (--deploy, read by analysis/utils/preprocess_paths.R), or --nz-rule /
# --nz-threshold given explicitly; one of the two is required.
#
# Usage:
#   bash analysis/bash/stages/01_preprocess.sh --deploy D.yml [--indir DIR] [--after JOBID] [--cpus 8] [--plan]
#   bash analysis/bash/stages/01_preprocess.sh --nz-rule r2 --nz-threshold 60 [...]
#
# Prints OUTBREAKS_JOB_ID=, MONTHLY_JOB_ID=, NZ_JOB_ID=.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

DEPLOY=""
INDIR=""
AFTER=""
NZ_RULE=""
NZ_THRESHOLD=""
CPUS=8
TIME="04:00:00"
PLAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deploy)       DEPLOY="$2";       shift 2 ;;
    --indir)        INDIR="$2";        shift 2 ;;
    --after)        AFTER="$2";        shift 2 ;;
    --nz-rule)      NZ_RULE="$2";      shift 2 ;;
    --nz-threshold) NZ_THRESHOLD="$2"; shift 2 ;;
    --cpus)         CPUS="$2";         shift 2 ;;
    --time)         TIME="$2";         shift 2 ;;
    --plan)         PLAN=true;         shift   ;;
    *) echo "Unknown argument: $1" >&2
       echo "Usage: $0 (--deploy D.yml | --nz-rule R --nz-threshold N) [--indir DIR] [--after JOBID] [--cpus N] [--time T] [--plan]" >&2
       exit 1 ;;
  esac
done

if [[ -n "$DEPLOY" ]]; then
  KV=$(Rscript analysis/utils/preprocess_paths.R --deploy "$DEPLOY")
  [[ -z "$NZ_RULE"      ]] && NZ_RULE=$(echo "$KV" | grep '^NZ_RULE=' | cut -d= -f2-)
  [[ -z "$NZ_THRESHOLD" ]] && NZ_THRESHOLD=$(echo "$KV" | grep '^NZ_THRESHOLD=' | cut -d= -f2-)
fi
if [[ -z "$NZ_THRESHOLD" ]]; then
  echo "Error: pass --deploy D.yml or --nz-rule/--nz-threshold (the national-zero months to build)." >&2
  exit 1
fi

IND=()
[[ -n "$INDIR" ]] && IND=(--indir "$INDIR")
PAYLOAD=analysis/bash/slurm/submit_preprocess.sh

# submit NAME CPUS MEM AFTER_ID -- SCRIPT ARGS...; echoes the job id
submit() {
  local name=$1 cpus=$2 mem=$3 after=$4; shift 5
  local dep=()
  [[ -n "$after" ]] && dep=(--dependency="afterok:${after}")
  cholera_sbatch_cmd LONG --parsable --job-name="$name" -c "$cpus" --mem="$mem" --time="$TIME" \
    ${dep[@]+"${dep[@]}"} "$PAYLOAD" "$@"
  local cmd=("${CHOLERA_SBATCH_CMD[@]}")
  if [[ "$PLAN" == "true" ]]; then
    echo "PLAN: ${cmd[*]}" >&2
    echo "PLAN_${name}"
  else
    "${cmd[@]}"
  fi
}

[[ "$PLAN" == "true" ]] || mkdir -p "$CHOLERA_LOG_DIR"
OB=$(submit prep_01_outbreaks "$CPUS" 32G "$AFTER" -- \
       analysis/01_preprocess_outbreaks.R -r TRUE -f FALSE "${IND[@]}")
echo "OUTBREAKS_JOB_ID=${OB}"
MO=$(submit prep_01b_monthly 1 8G "$OB" -- \
       analysis/01b_preprocess_outbreaks_monthly.R "${IND[@]}")
echo "MONTHLY_JOB_ID=${MO}"
RULE=()
[[ -n "$NZ_RULE" ]] && RULE=(--rule "$NZ_RULE")   # "" = the legacy rule = 01_nz's default
NZ=$(submit prep_01_national_zero 1 8G "$OB" -- \
       analysis/01_preprocess_national_zero_months.R "${IND[@]}" \
       "${RULE[@]}" --zero_threshold "$NZ_THRESHOLD")
echo "NZ_JOB_ID=${NZ}"
