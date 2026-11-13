#!/bin/bash
# stages/07_hindcast.sh — submit the season-by-season forecast-skill test
# (hindcast) as a set of spatPomp SLURM ARRAY JOBS, one per
# (ENSO group x launch year) config produced by make_hindcast_configs.R, with
# one ARRAY TASK per COUNTRY in that config's simulator.scope.countries and
# intra-task parallelism over REPLICATE CHUNKS only.
#
# This mirrors stages/07_forecast_groups.sh: it submits each config through
# analysis/bash/slurm/submit_forecast_groups.sh (the spatPomp array script
# that runs 07b_simulator_spatpomp.R per country) — NOT the v1 task-array path.
# The hindcast configs declare simulator.backend: spatpomp and carry the
# filter.handoff_month / enso.forecast_source: observed seams, which only the
# spatpomp engine reads; routing them through the v1 replay kernel would ignore
# both the backend key and the seams.
#
# Each hindcast_<group>_<year>.yml is an independent array job (distinct
# run_tag/run_id -> distinct parts/merged dirs) and can run fully in parallel;
# within one array job countries are independent array tasks, and within one
# array task replicate chunks are independent background processes.
#
# The config set is REQUIRED: either --deploy <prod_run<N>.yml> (config dir,
# hindcast manifest and partition_by are read from `deploy_prod_run.R --stage
# resolve`) or an explicit --config-dir. Explicit --config-dir / --manifest /
# --partition_by / --obs-panel win over --deploy. There is no run-specific
# default. --deploy also supplies the --postprocess observed panel (OBS_PANEL,
# the survival panel matching the run's episodes.key); without --deploy or
# --obs-panel it falls back to the lp_season panel.
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml):
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D                      # submit all
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --list               # print, don't submit
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --group enso_east    # one ENSO group
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --glob '*_2015.yml'  # one launch year
#   bash analysis/bash/stages/07_hindcast.sh --config-dir <dir>               # an explicit config set instead
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --n_replicates 300 --chunk_size 50 --Np 10000
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --dry-run            # forward --dry-run to the Rscript
#   bash analysis/bash/stages/07_hindcast.sh --deploy $D --postprocess --pp-workers 8 --pp-mem 64G
#                                     # chain 07c/07d --hindcast, 8 runs scored in parallel
#
# --n_replicates/--chunk_size/--Np/--dry-run/--overwrite are forwarded to
# slurm/submit_forecast_groups.sh and then (per replicate chunk) to
# 07b_simulator_spatpomp.R.

set -euo pipefail

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

CONFIG_DIR=""      # required: --config-dir or --deploy (no run-specific default)
DEPLOY=""
GLOB=""            # optional extra shell glob within CONFIG_DIR (e.g. '*_2015.yml')
GROUP=""           # optional ENSO group filter (e.g. enso_east)
SUBMITTER="analysis/bash/slurm/submit_forecast_groups.sh"
PP_SUBMITTER="analysis/bash/slurm/submit_hindcast_postprocess.sh"
LIST_ONLY=0
PARTITION_BY=""    # ""/country = per-country array; enso_region|registry = one whole-unit task per config
POSTPROCESS=0      # --postprocess: chain a whole-manifest 07c/07d job afterok all arrays
PP_WORKERS=""      # --pp-workers N: CPUs (= 07c/07d --workers) of that job (default: payload's -c 8)
PP_MEM=""          # --pp-mem SIZE: its memory (default: the payload's #SBATCH --mem, 120G)
DRYRUN=0           # tracked separately from EXTRA to gate the postprocess chain
MANIFEST=""        # 07c/07d --hindcast manifest; --manifest > --deploy > hindcast_horizons.yml
OBS_PANEL=""       # 07c/07d --obs-panel; --obs-panel > --deploy (matches episodes.key) > lp_season panel
EXTRA=""           # forwarded verbatim to the array submitter
CONCURRENCY=""     # --concurrency K: %K on each array (one-task arrays in registry/enso_region mode: no effect)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config-dir)   CONFIG_DIR="$2";  shift 2 ;;
    --deploy)       DEPLOY="$2";      shift 2 ;;
    --glob)         GLOB="$2";        shift 2 ;;
    --group)        GROUP="$2";       shift 2 ;;
    --partition_by) PARTITION_BY="$2"; EXTRA="$EXTRA --partition_by $2"; shift 2 ;;
    --n_replicates|--chunk_size|--Np) EXTRA="$EXTRA $1 $2"; shift 2 ;;
    --dry-run)      EXTRA="$EXTRA --dry-run"; DRYRUN=1; shift ;;
    --overwrite)    EXTRA="$EXTRA --overwrite"; shift ;;
    --postprocess)  POSTPROCESS=1;    shift ;;
    --pp-workers)   PP_WORKERS="$2";  shift 2 ;;
    --pp-mem)       PP_MEM="$2";      shift 2 ;;
    --manifest)     MANIFEST="$2";    shift 2 ;;
    --obs-panel)    OBS_PANEL="$2";   shift 2 ;;
    --list)         LIST_ONLY=1;      shift ;;
    --concurrency)  CONCURRENCY="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 (--deploy DEPLOY.yml | --config-dir DIR) [--group G] [--glob PAT] [--partition_by enso_region|registry] [--n_replicates N] [--chunk_size N] [--Np N] [--dry-run] [--overwrite] [--postprocess] [--manifest FILE] [--obs-panel FILE] [--concurrency K] [--list]" >&2
      exit 1 ;;
  esac
done

# Run paths come from the deploy config (R owns them; KEY=value contract of
# deploy_prod_run.R --stage resolve). Explicit flags win.
if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  [[ -n "$CONFIG_DIR" ]] || CONFIG_DIR=$(deploy_key HINDCAST_CONFIG_DIR)
  [[ -n "$MANIFEST" ]]   || MANIFEST=$(deploy_key HINDCAST_MANIFEST)
  [[ -n "$OBS_PANEL" ]]  || OBS_PANEL=$(deploy_key OBS_PANEL)
  _dp="$(deploy_key PARTITION_BY)"
  if [[ -z "$PARTITION_BY" && ( "$_dp" == "enso_region" || "$_dp" == "registry" ) ]]; then
    PARTITION_BY="$_dp"; EXTRA="$EXTRA --partition_by $_dp"
  fi
fi
MANIFEST="${MANIFEST:-analysis/configs/hindcast_horizons.yml}"
OBS_PANEL="${OBS_PANEL:-generated_data/outbreaks_survival_data_base-1990_2019.rds}"
if [[ -z "$CONFIG_DIR" ]]; then
  echo "Error: no config set. Pass --deploy analysis/configs/deploy/prod_run<N>.yml or --config-dir DIR." >&2
  exit 1
fi
echo "Config dir: ${CONFIG_DIR}${DEPLOY:+  (deploy: $DEPLOY)}"

if [[ ${POSTPROCESS} -eq 1 && ! -f "$PP_SUBMITTER" ]]; then
  echo "Error: postprocess submitter not found: $PP_SUBMITTER" >&2
  exit 1
fi

if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "Error: config dir not found: $CONFIG_DIR" >&2
  echo "       (run: Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_hindcast first)" >&2
  exit 1
fi
if [[ ! -f "$SUBMITTER" ]]; then
  echo "Error: array submitter not found: $SUBMITTER" >&2
  exit 1
fi

# Build the select pattern: prefer explicit --glob, else --group, else all.
if [[ -n "$GLOB" ]]; then
  PATTERN="$GLOB"
elif [[ -n "$GROUP" ]]; then
  PATTERN="hindcast_${GROUP}_*.yml"
else
  PATTERN="hindcast_*.yml"
fi

shopt -s nullglob
CONFIGS=("$CONFIG_DIR"/$PATTERN)
shopt -u nullglob

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
  echo "Error: no configs matched '$CONFIG_DIR/$PATTERN'." >&2
  exit 1
fi

mkdir -p "$CHOLERA_LOG_DIR"

# Extract simulator.scope.countries (comma-separated ISO3 list) from a config —
# drives both the --array=0-N-1 size and the --countries list each array task
# indexes into via SLURM_ARRAY_TASK_ID.
countries_for_config() {
  sed -n 's/^[[:space:]]*countries:[[:space:]]*\(.*\)$/\1/p' "$1" | tr -d ' \r' | head -1
}

echo "--- Hindcast submission: ${#CONFIGS[@]} config(s) from $CONFIG_DIR"
[[ -n "$EXTRA" ]] && echo "    Forwarding to submitter:${EXTRA}"

submitted=0
ARRAY_IDS=()    # job ids of the submitted array jobs, for the postprocess afterok dep
for cfg in "${CONFIGS[@]}"; do
  stem="$(basename "$cfg" .yml)"        # e.g. hindcast_enso_east_2015
  rest="${stem#hindcast_}"              # -> enso_east_2015
  year="${rest##*_}"                    # -> 2015
  group="${rest%_*}"                    # -> enso_east
  short="${group#enso_}"                # -> east
  jobname="sphc_${short}_${year}"

  countries="$(countries_for_config "$cfg")"
  if [[ -z "$countries" ]]; then
    echo "ERROR: could not extract simulator.scope.countries from ${cfg}" >&2
    exit 1
  fi
  IFS=',' read -ra country_arr <<< "$countries"
  n_countries=${#country_arr[@]}

  # Region mode: the whole ENSO group is ONE forward partition → one array task
  # per config (no per-country fan-out). The submitter ignores --countries under
  # --partition_by enso_region. Country mode: one array task per ISO3.
  if [[ "${PARTITION_BY}" == "enso_region" || "${PARTITION_BY}" == "registry" ]]; then
    array_spec="0-0"
    scope_desc="${PARTITION_BY} unit (1 task; ${n_countries} member countries: ${countries})"
  else
    array_spec="0-$((n_countries - 1))"
    scope_desc="${n_countries} countries: ${countries}"
  fi

  if [[ ${LIST_ONLY} -eq 1 ]]; then
    echo "[list] ${jobname}  --array=${array_spec}  (${scope_desc})  <-  ${cfg}"
    continue
  fi

  echo "--- Submitting ${jobname} --array=${array_spec} (${scope_desc}) (${cfg}) ..."
  jid=$(cholera_sbatch SHORT --parsable \
    --job-name="${jobname}" \
    --array="${array_spec}${CONCURRENCY:+%${CONCURRENCY}}" \
    "$SUBMITTER" \
    --config "${cfg}" --countries "${countries}" ${EXTRA})
  echo "    array job id: ${jid}"
  ARRAY_IDS+=("${jid}")
  submitted=$((submitted + 1))
done

[[ ${LIST_ONLY} -eq 1 ]] && exit 0

echo ""
echo "Submitted ${submitted} array job(s)."

# --- Optional: chain the hindcast postprocess + figures ----------------------
# ONE whole-manifest job (07c --hindcast scores every run in the manifest, 07d
# --hindcast draws the skill figures), gated on ALL submitted array jobs via
# afterok. Skipped under --dry-run (the array tasks only validate+exit, so there
# would be nothing to score) and when nothing was submitted.
if [[ ${POSTPROCESS} -eq 1 ]]; then
  if [[ ${DRYRUN} -eq 1 ]]; then
    echo "--- (postprocess skipped: --dry-run submits no real work to score)"
  elif [[ ${submitted} -eq 0 || ${#ARRAY_IDS[@]} -eq 0 ]]; then
    echo "--- (postprocess skipped: no array jobs were submitted)"
  else
    dep="afterok"
    for id in "${ARRAY_IDS[@]}"; do dep="${dep}:${id}"; done
    echo "--- Submitting hindcast postprocess (dependency=${dep}) ..."
    PP_JOB_ID=$(cholera_sbatch SHORT --parsable --dependency="${dep}" \
      --job-name="hc_postprocess" ${PP_WORKERS:+-c ${PP_WORKERS}} ${PP_MEM:+--mem=${PP_MEM}} \
      "$PP_SUBMITTER" "${MANIFEST}" "${OBS_PANEL}")
    echo "Postprocess job id: ${PP_JOB_ID}  (manifest=${MANIFEST}, obs-panel=${OBS_PANEL})"
  fi
fi

echo "Monitor with: squeue -u \$USER --format='%.10i %.14j %.8T %.10M %R'"
