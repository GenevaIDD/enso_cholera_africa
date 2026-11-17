#!/bin/bash
# stages/07_event_attribution.sh — submit the past-ENSO-event attribution runs
# (factual + counterfactual configs per event x unit, from deploy_prod_run.R
# --stage gen_events) as SLURM arrays through slurm/submit_forecast_groups.sh
# (07b_simulator_spatpomp.R).
#
# Factual configs are hindcasts (observed ENSO, interior smoothing handoff);
# counterfactuals differ only by simulator.enso.counterfactual (Lenssen
# analogs from a per-event cache). partition_by registry / enso_region: one
# task per config; otherwise one per country. No postprocess: 09
# (report_event_attrib) reads the partitions directly.
#
# A config whose run has _SUCCESS in every partition
# (<outdir>/sim_spatpomp/<run_tag>/parts/tu=*/[repoff=*/]_SUCCESS) is skipped
# (07b does not protect finished runs); --list marks it [skip], --overwrite
# submits it anyway.
#
# The config set is REQUIRED: --deploy <prod_run<N>.yml> (EVENT_CONFIG_DIR,
# PARTITION_BY from deploy_prod_run.R --stage resolve) or --config-dir.
# Explicit flags win.
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml):
#   bash analysis/bash/stages/07_event_attribution.sh --deploy $D --list   # print only
#   bash analysis/bash/stages/07_event_attribution.sh --deploy $D          # submit all
#   bash analysis/bash/stages/07_event_attribution.sh --deploy $D --event elnino_2015
#   bash analysis/bash/stages/07_event_attribution.sh --deploy $D --arm factual   # or counterfactual
#   bash analysis/bash/stages/07_event_attribution.sh --deploy $D --n_replicates 100 --chunk_size 50
#
# --n_replicates / --chunk_size / --Np / --dry-run / --overwrite go through
# the payload to 07b.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

CONFIG_DIR=""
DEPLOY=""
EVENT=""           # optional event id filter (e.g. elnino_2015)
ARM=""             # optional: factual | counterfactual
SUBMITTER="analysis/bash/slurm/submit_forecast_groups.sh"
LIST_ONLY=0
PARTITION_BY=""
EXTRA=""
CONCURRENCY=""
OVERWRITE=0        # --overwrite: also resubmit configs whose run is complete

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config-dir)   CONFIG_DIR="$2";  shift 2 ;;
    --deploy)       DEPLOY="$2";      shift 2 ;;
    --event)        EVENT="$2";       shift 2 ;;
    --arm)          ARM="$2";         shift 2 ;;
    --partition_by) PARTITION_BY="$2"; EXTRA="$EXTRA --partition_by $2"; shift 2 ;;
    --n_replicates|--chunk_size|--Np) EXTRA="$EXTRA $1 $2"; shift 2 ;;
    --dry-run)      EXTRA="$EXTRA --dry-run"; shift ;;
    --overwrite)    EXTRA="$EXTRA --overwrite"; OVERWRITE=1; shift ;;
    --list)         LIST_ONLY=1;      shift ;;
    --concurrency)  CONCURRENCY="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 (--deploy DEPLOY.yml | --config-dir DIR) [--event ID] [--arm factual|counterfactual] [--partition_by enso_region|registry] [--n_replicates N] [--chunk_size N] [--Np N] [--dry-run] [--overwrite] [--concurrency K] [--list]" >&2
      exit 1 ;;
  esac
done

if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  [[ -n "$CONFIG_DIR" ]] || CONFIG_DIR=$(deploy_key EVENT_CONFIG_DIR)
  _dp="$(deploy_key PARTITION_BY)"
  if [[ -z "$PARTITION_BY" && ( "$_dp" == "enso_region" || "$_dp" == "registry" ) ]]; then
    PARTITION_BY="$_dp"; EXTRA="$EXTRA --partition_by $_dp"
  fi
fi
if [[ -z "$CONFIG_DIR" ]]; then
  echo "Error: no config set. Pass --deploy analysis/configs/deploy/prod_run<N>.yml (with deploy.event_attribution set) or --config-dir DIR." >&2
  exit 1
fi
if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "Error: config dir not found: $CONFIG_DIR" >&2
  echo "       (run: Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_events first)" >&2
  exit 1
fi
if [[ -n "$ARM" && "$ARM" != "factual" && "$ARM" != "counterfactual" ]]; then
  echo "Error: --arm must be factual or counterfactual, got '$ARM'." >&2
  exit 1
fi
echo "Config dir: ${CONFIG_DIR}${DEPLOY:+  (deploy: $DEPLOY)}"

shopt -s nullglob
CONFIGS=("$CONFIG_DIR"/event_${EVENT:-*}_*.yml)
shopt -u nullglob
if [[ ${#CONFIGS[@]} -eq 0 ]]; then
  echo "Error: no configs matched '$CONFIG_DIR/event_${EVENT:-*}_*.yml'." >&2
  exit 1
fi

mkdir -p "$CHOLERA_LOG_DIR"

countries_for_config() {
  sed -n 's/^[[:space:]]*countries:[[:space:]]*\(.*\)$/\1/p' "$1" | tr -d ' \r' | head -1
}

# The counterfactual config's Lenssen cache must exist before its run starts
# (07b would stop at the ENSO sampler). Read from the config, never built here.
lenssen_cache_for_config() {
  sed -n 's/^[[:space:]]*cache_path:[[:space:]]*\(.*\)$/\1/p' "$1" | tr -d ' \r"' | head -1
}

# The parts dir of a config's run: <outdir>/sim_spatpomp/<simulator.run_tag>/parts
# (07b's spatpomp_parts_dir; outdir defaults to generated_data/). Same reader as
# stages/07_forecast_groups.sh.
parts_dir_for_config() {
  local tag outdir
  tag=$(sed -n 's/^  run_tag:[[:space:]]*//p' "$1" | head -1 | tr -d "\"' \r")
  outdir=$(sed -n 's/^outdir:[[:space:]]*//p' "$1" | head -1 | tr -d "\"' \r")
  echo "${outdir:-generated_data/}"
  echo "${tag}"
}

# 0 when every partition of the run has its _SUCCESS sentinel: each tu=<unit> dir
# (or each of its repoff=<chunk> dirs) and at least $2 tu dirs.
run_complete() {
  local parts="$1" min_tu="$2" n_tu=0 tu c
  [[ -d "$parts" ]] || return 1
  for tu in "$parts"/tu=*; do
    [[ -d "$tu" ]] || continue
    n_tu=$((n_tu + 1))
    if compgen -G "$tu/repoff=*" >/dev/null; then
      for c in "$tu"/repoff=*; do [[ -f "$c/_SUCCESS" ]] || return 1; done
    else
      [[ -f "$tu/_SUCCESS" ]] || return 1
    fi
  done
  (( n_tu >= min_tu ))
}

submitted=0
skipped=0
for cfg in "${CONFIGS[@]}"; do
  stem="$(basename "$cfg" .yml)"          # event_elnino_2015_hmmCE[_cf-neutral]
  if [[ "$stem" == *_cf-* ]]; then arm="counterfactual"; else arm="factual"; fi
  [[ -n "$ARM" && "$ARM" != "$arm" ]] && continue
  rest="${stem#event_}"                   # elnino_2015_hmmCE[_cf-neutral]
  jobname="spev_${rest}"

  countries="$(countries_for_config "$cfg")"
  if [[ -z "$countries" ]]; then
    echo "ERROR: could not extract simulator.scope.countries from ${cfg}" >&2
    exit 1
  fi
  IFS=',' read -ra country_arr <<< "$countries"
  n_countries=${#country_arr[@]}

  # Skip-complete guard (see header); before the cache check: a finished run
  # needs no Lenssen cache.
  if [[ ${OVERWRITE} -eq 0 ]]; then
    mapfile -t _pd < <(parts_dir_for_config "$cfg")
    _parts="${_pd[0]%/}/sim_spatpomp/${_pd[1]}/parts"
    _min_tu=1
    [[ "${PARTITION_BY}" == "enso_region" || "${PARTITION_BY}" == "registry" ]] || _min_tu=${n_countries}
    if [[ -n "${_pd[1]}" ]] && run_complete "${_parts}" "${_min_tu}"; then
      echo "[skip] ${jobname}: run complete (${_parts}); --overwrite resubmits it  <-  ${cfg}"
      skipped=$((skipped + 1))
      continue
    fi
  fi

  if [[ "$arm" == "counterfactual" ]]; then
    cache="$(lenssen_cache_for_config "$cfg")"
    if [[ -z "$cache" || ! -f "$cache" ]]; then
      echo "ERROR: ${cfg}: Lenssen cache '${cache}' not found; build it first (deploy --stage gen_events prints the command)." >&2
      exit 1
    fi
  fi

  if [[ "${PARTITION_BY}" == "enso_region" || "${PARTITION_BY}" == "registry" ]]; then
    array_spec="0-0"
    scope_desc="${PARTITION_BY} unit (1 task; ${n_countries} member countries)"
  else
    array_spec="0-$((n_countries - 1))"
    scope_desc="${n_countries} countries"
  fi

  if [[ ${LIST_ONLY} -eq 1 ]]; then
    echo "[list] ${jobname}  --array=${array_spec}  (${arm}; ${scope_desc})  <-  ${cfg}"
    continue
  fi

  echo "--- Submitting ${jobname} --array=${array_spec} (${arm}; ${scope_desc}) ..."
  jid=$(cholera_sbatch SHORT --parsable \
    --job-name="${jobname}" \
    --array="${array_spec}${CONCURRENCY:+%${CONCURRENCY}}" \
    "$SUBMITTER" \
    --config "${cfg}" --countries "${countries}" ${EXTRA})
  echo "    array job id: ${jid}"
  submitted=$((submitted + 1))
done

[[ ${LIST_ONLY} -eq 1 ]] && { echo "[list] ${skipped} complete run(s) would be skipped."; exit 0; }
echo ""
echo "Submitted ${submitted} array job(s); skipped ${skipped} complete run(s)."
echo "Monitor with: squeue -u \$USER --format='%.10i %.24j %.8T %.10M %R'"
