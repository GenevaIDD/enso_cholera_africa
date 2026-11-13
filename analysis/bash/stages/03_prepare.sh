#!/bin/bash
# stages/03_prepare.sh — submit the 03 model-dataset steps as Slurm jobs, after
# the 00-02 data leg (campaigns/preprocess.sh), so a rollout's inputs are built
# on the cluster that fits them.
#
#   [--after IDS] -> 03_prepare_occ.R            default panel (sim panels read it)
#                 -> 03_prepare_occ.R            the deploy's panel (national zeros,
#                                                at-risk start, ADM1 overlap, coverage report)
#                 -> 03_prepare_siz.R -> 03_prepare_dur.R -> 03_prepare_mon.R
#                 -> 03_prepare_sim_panel.R -> 03_prepare_sim_panel_projection.R
#                 -> 03_prepare_sim_panel_extended.R (--extended-end)
#   [--after IDS] -> 03_prepare_who_monthly.R
#   [--after IDS] -> 03_prepare_rain_adm1.R      (reads climate_stats.sqlite)
#
# Both occurrence calls pass --redo TRUE (03_prepare_occ.R defaults to FALSE, 03_prepare_siz.R
# to TRUE): with FALSE the occurrence prep can reuse a stale matched_ob_gids_seasonality_occ.rds,
# whose episode digest then differs from the size / duration / monthly preps and deploy_fits
# refuses the pair (assertEpisodeCoherence).
#
# The model and simulation steps run one after another: they share caches in
# generated_data/ (seasonality matches, contagion distances) and the sim panels
# re-download ENSO into the caches the others read. Back up the ENSO caches
# before submitting (the panels' stamps do not cover ENSO values). Every step
# runs with --force: the panels of the previous vintage are replaced.
#
# The panel spec, the projection's ENSO continuation and its CFSv2 correction
# come from the deploy YAML (analysis/utils/preprocess_paths.R --deploy), so
# bash names no artefact and no run.
#
# Usage:
#   bash analysis/bash/stages/03_prepare.sh --deploy D.yml [--after JOBID[:JOBID...]]
#       [--extended-end YYYY-MM-DD] [--baseline 1990_2019] [--enso-max-lead 11]
#       [--contagion-partition-source enso_grouping,hmm_clusters] [--who-zero-iso3 TZA]
#       [--time 04:00:00] [--plan]
#
# Prints one *_JOB_ID= line per step and PREPARE_JOB_IDS= (all of them,
# colon-joined: run_pipeline.sh --after takes it).

set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/../../.."

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

DEPLOY=""
AFTER=""
EXT_END=""
BASELINE="1990_2019"
ENSO_MAX_LEAD=11
PART_SOURCE="enso_grouping,hmm_clusters"   # the contagion registries the panels carry
WHO_ZERO_ISO3="TZA"                          # unchanged WHO cumulative = zero only here
TIME="04:00:00"
PLAN=false

usage() {
  echo "Usage: $0 --deploy D.yml [--after IDS] [--extended-end YYYY-MM-DD] [--baseline B] [--enso-max-lead N] [--contagion-partition-source S] [--who-zero-iso3 ISO3] [--time T] [--plan]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deploy)                     DEPLOY="$2";        shift 2 ;;
    --after)                      AFTER="$2";         shift 2 ;;
    --extended-end)               EXT_END="$2";       shift 2 ;;
    --baseline)                   BASELINE="$2";      shift 2 ;;
    --enso-max-lead)              ENSO_MAX_LEAD="$2"; shift 2 ;;
    --contagion-partition-source) PART_SOURCE="$2";   shift 2 ;;
    --who-zero-iso3)              WHO_ZERO_ISO3="$2"; shift 2 ;;
    --time)                       TIME="$2";          shift 2 ;;
    --plan)                       PLAN=true;          shift   ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$DEPLOY" && -f "$DEPLOY" ]] || { echo "Error: --deploy D.yml is required and must exist." >&2; usage; exit 1; }
[[ -z "$AFTER" || "$AFTER" =~ ^[A-Za-z0-9_]+(:[A-Za-z0-9_]+)*$ ]] || { echo "Error: --after must be colon-joined job ids." >&2; exit 1; }
[[ -z "$EXT_END" || "$EXT_END" =~ ^[0-9]{4}-[0-9]{2}-01$ ]] || { echo "Error: --extended-end must be a month start, YYYY-MM-01." >&2; exit 1; }

KV=$(Rscript analysis/utils/preprocess_paths.R --deploy "$DEPLOY")
get() { echo "$KV" | grep "^$1=" | head -1 | cut -d= -f2- || true; }
EPISODE_KEY=$(get EPISODE_KEY)
NZ_RULE=$(get NZ_RULE)
NZ_THRESHOLD=$(get NZ_THRESHOLD)
AT_RISK_START=$(get AT_RISK_START)
ADM1_OVERLAP=$(get ADM1_OVERLAP)
read -r -a ENSO_CONT <<< "$(get ENSO_CONT_ARGS)"
read -r -a CFSV2 <<< "$(get CFSV2_ARGS)"
[[ -n "$EPISODE_KEY" ]] || { echo "Error: no EPISODE_KEY from preprocess_paths.R --deploy $DEPLOY." >&2; exit 1; }

PAYLOAD=analysis/bash/slurm/submit_preprocess.sh

# submit NAME CPUS MEM AFTER_IDS -- SCRIPT ARGS...; echoes the job id
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
COMMON=(-b "$BASELINE" --episode_key "$EPISODE_KEY")
IDS=()

OCC0=$(submit prep_03_occ_default 1 32G "$AFTER" -- \
         analysis/03_prepare_occ.R "${COMMON[@]}" --contagion_partition_source "$PART_SOURCE" \
         --enso_max_lead "$ENSO_MAX_LEAD" --redo TRUE --force)
echo "OCC_DEFAULT_JOB_ID=${OCC0}"; IDS+=("$OCC0")

NZ=()
[[ -n "$NZ_THRESHOLD" ]] && NZ=(--national_zero_fill --national_zero_threshold "$NZ_THRESHOLD")
[[ -n "$NZ_RULE" ]] && NZ+=(--national_zero_rule "$NZ_RULE")
OCC1=$(submit prep_03_occ_deploy 1 32G "$OCC0" -- \
         analysis/03_prepare_occ.R "${COMMON[@]}" --contagion_partition_source "$PART_SOURCE" \
         --enso_max_lead "$ENSO_MAX_LEAD" "${NZ[@]}" --at_risk_start "$AT_RISK_START" \
         --adm1_overlap "$ADM1_OVERLAP" --contagion_coverage_report --redo TRUE --force)
echo "OCC_DEPLOY_JOB_ID=${OCC1}"; IDS+=("$OCC1")

SIZ=$(submit prep_03_siz 1 32G "$OCC1" -- analysis/03_prepare_siz.R "${COMMON[@]}" --force)
echo "SIZ_JOB_ID=${SIZ}"; IDS+=("$SIZ")
DUR=$(submit prep_03_dur 1 32G "$SIZ" -- \
        analysis/03_prepare_dur.R "${COMMON[@]}" --enso_max_lead "$ENSO_MAX_LEAD" --force)
echo "DUR_JOB_ID=${DUR}"; IDS+=("$DUR")
MON=$(submit prep_03_mon 1 32G "$DUR" -- analysis/03_prepare_mon.R "${COMMON[@]}" --force)
echo "MON_JOB_ID=${MON}"; IDS+=("$MON")

SIM=$(submit prep_03_sim_panel 1 32G "$MON" -- \
        analysis/03_prepare_sim_panel.R "${COMMON[@]}" --enso_max_lead "$ENSO_MAX_LEAD" --force)
echo "SIM_PANEL_JOB_ID=${SIM}"; IDS+=("$SIM")
PROJ=$(submit prep_03_sim_projection 1 32G "$SIM" -- \
         analysis/03_prepare_sim_panel_projection.R "${COMMON[@]}" --enso_max_lead "$ENSO_MAX_LEAD" \
         "${ENSO_CONT[@]}" "${CFSV2[@]}" --force)
echo "SIM_PROJECTION_JOB_ID=${PROJ}"; IDS+=("$PROJ")
if [[ -n "$EXT_END" ]]; then
  EXT=$(submit prep_03_sim_extended 1 32G "$PROJ" -- \
          analysis/03_prepare_sim_panel_extended.R --episode_key "$EPISODE_KEY" \
          --enso_max_lead "$ENSO_MAX_LEAD" --end "$EXT_END" --force)
  echo "SIM_EXTENDED_JOB_ID=${EXT}"; IDS+=("$EXT")
else
  echo "note: no --extended-end: 03_prepare_sim_panel_extended.R not submitted" >&2
fi

WHO=$(submit prep_03_who_monthly 1 8G "$AFTER" -- \
        analysis/03_prepare_who_monthly.R --zero_iso3 "$WHO_ZERO_ISO3" --force)
echo "WHO_MONTHLY_JOB_ID=${WHO}"; IDS+=("$WHO")
RAIN=$(submit prep_03_rain_adm1 1 64G "$AFTER" -- analysis/03_prepare_rain_adm1.R)
echo "RAIN_ADM1_JOB_ID=${RAIN}"; IDS+=("$RAIN")

echo "PREPARE_JOB_IDS=$(IFS=:; echo "${IDS[*]}")"
