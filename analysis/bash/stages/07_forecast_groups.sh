#!/bin/bash
# stages/07_forecast_groups.sh — submit the per-ENSO-group forecast and
# counterfactual spatPomp runs: one SLURM array job per (arm, group) config
# (made by utils/make_forecast_configs.R), one array task per country of the
# config's simulator.scope.countries, replicate chunks in parallel inside a
# task (--chunk_size). Jobs are independent; no dependency chain between them.
#
# Run tools/check_forecast_groups_files.sh first. TUs without a historical
# smoothing cache are dropped from their group's run, not an error
# (07b_simulator_spatpomp.R::load_smoothing_distribution()).
#
# The config set is REQUIRED: --deploy <prod_run<N>.yml> (config dir and
# partition_by from deploy_prod_run.R --stage resolve) or --config-dir /
# CONFIG_DIR, which wins. No default, so stale configs are never paired with
# the current forecast_horizons.yml.
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml):
#   bash analysis/bash/stages/07_forecast_groups.sh --deploy $D            # all
#   bash analysis/bash/stages/07_forecast_groups.sh --deploy $D --list     # print only
#   bash analysis/bash/stages/07_forecast_groups.sh --deploy $D --arm all --postprocess
#   bash analysis/bash/stages/07_forecast_groups.sh --deploy $D --groups east,south
#   bash analysis/bash/stages/07_forecast_groups.sh --deploy $D --n_replicates 300 --Np 10000
#
# --arm: factual (forecast_<unit>.yml), attrib (forecast_<unit>_attrib.yml,
# from --stage gen_attrib), counterfactual (counterfactual_<unit>_set-<s>.yml),
# both [default], all (every arm present).
# --n_replicates / --chunk_size / --Np / --dry-run / --overwrite go through
# slurm/submit_forecast_groups.sh to 07b (replicates default to the config's
# simulator.n_replicates; --dry-run validates config + artefacts and exits).
# --postprocess: per-config 07c --aggregate + 07d job afterok its array
# (slurm/submit_forecast_postprocess.sh); skipped under --dry-run. The
# cross-group combine lives in campaigns/forecast_groups_pipeline.sh.
#
# A config whose run has _SUCCESS in every partition
# (generated_data/sim_spatpomp/<run_tag>/parts/tu=*/[repoff=*/]_SUCCESS) is
# skipped (07b does not protect finished runs); --overwrite submits it anyway.
# Last stdout line: FORECAST_JOB_IDS=<id>:<id>:... when anything was submitted.

set -euo pipefail

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

# No default config set: it must come from --config-dir / CONFIG_DIR or
# --deploy (resolved below): a hardcoded default would drift from the manifest.
CONFIG_DIR="${CONFIG_DIR:-}"
DEPLOY=""
ARM="both"
GROUPS_FILTER=""
SETS_FILTER=""
LIST_ONLY=0
CONCURRENCY=""   # --concurrency K: %K on each array (one-task arrays in registry/enso_region mode: no effect)
PARTITION_BY=""    # ""/country = per-country array; enso_region|registry = one whole-unit task per config
PP_SUBMITTER="analysis/bash/slurm/submit_forecast_postprocess.sh"
POSTPROCESS=0      # --postprocess: chain a PER-CONFIG 07c/07d job afterok that config's array
DRYRUN=0           # tracked separately from EXTRA to gate the postprocess chain
OVERWRITE=0        # --overwrite: also resubmit configs whose run is complete
EXTRA=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arm)
      ARM="$2"; shift 2 ;;
    --config-dir)
      CONFIG_DIR="$2"; shift 2 ;;
    --deploy)
      DEPLOY="$2"; shift 2 ;;
    --groups)
      GROUPS_FILTER="$2"; shift 2 ;;
    --analog-sets|--analog_sets)
      SETS_FILTER="$2"; shift 2 ;;
    --partition_by)
      PARTITION_BY="$2"; EXTRA="$EXTRA --partition_by $2"; shift 2 ;;
    --n_replicates|--chunk_size|--Np)
      EXTRA="$EXTRA $1 $2"; shift 2 ;;
    --dry-run)
      EXTRA="$EXTRA --dry-run"; DRYRUN=1; shift ;;
    --overwrite)
      EXTRA="$EXTRA --overwrite"; OVERWRITE=1; shift ;;
    --postprocess)
      POSTPROCESS=1; shift ;;
    --list)
      LIST_ONLY=1; shift ;;
    --concurrency)
      CONCURRENCY="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 (--deploy DEPLOY.yml | --config-dir DIR) [--arm factual|attrib|counterfactual|both|all] [--groups east,south,...] [--analog-sets all,nina,...] [--partition_by enso_region|registry] [--n_replicates N] [--chunk_size N] [--Np N] [--dry-run] [--overwrite] [--postprocess] [--concurrency K] [--list]" >&2
      exit 1 ;;
  esac
done

# Run paths come from the deploy config (R owns them; KEY=value contract of
# deploy_prod_run.R --stage resolve). Explicit --config-dir / --partition_by win.
if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  [[ -n "$CONFIG_DIR" ]] || CONFIG_DIR=$(deploy_key FORECAST_CONFIG_DIR)
  _dp="$(deploy_key PARTITION_BY)"
  if [[ -z "$PARTITION_BY" && ( "$_dp" == "enso_region" || "$_dp" == "registry" ) ]]; then
    PARTITION_BY="$_dp"; EXTRA="$EXTRA --partition_by $_dp"
  fi
fi
if [[ -z "$CONFIG_DIR" ]]; then
  echo "Error: no config set. Pass --deploy analysis/configs/deploy/prod_run<N>.yml or --config-dir DIR (or set CONFIG_DIR)." >&2
  exit 1
fi
echo "Config dir: ${CONFIG_DIR}${DEPLOY:+  (deploy: $DEPLOY)}"

if [[ ${POSTPROCESS} -eq 1 && ! -f "$PP_SUBMITTER" ]]; then
  echo "Error: postprocess submitter not found: $PP_SUBMITTER" >&2
  exit 1
fi

case "${ARM}" in
  factual)        ARMS=(factual) ;;
  attrib)         ARMS=(attrib) ;;
  counterfactual) ARMS=(counterfactual) ;;
  both)           ARMS=(factual counterfactual) ;;
  all)            ARMS=(factual attrib counterfactual) ;;
  *)
    echo "--arm must be one of: factual, attrib, counterfactual, both, all" >&2
    exit 1 ;;
esac

# The configs of one arm (one per unit, sorted by the glob).
configs_for_arm() {
  local f
  case "$1" in
    factual)        for f in "${CONFIG_DIR}"/forecast_*.yml; do
                      [[ -e "$f" && "$f" != *_attrib.yml ]] && echo "$f"; done ;;
    attrib)         for f in "${CONFIG_DIR}"/forecast_*_attrib.yml; do
                      [[ -e "$f" ]] && echo "$f"; done ;;
    counterfactual) for f in "${CONFIG_DIR}"/counterfactual_*.yml; do
                      [[ -e "$f" ]] && echo "$f"; done ;;
  esac
  return 0
}

# The parts dir of a config's run: <outdir>/sim_spatpomp/<simulator.run_tag>/parts
# (07b's spatpomp_parts_dir; outdir defaults to generated_data/).
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

mkdir -p "$CHOLERA_LOG_DIR"

# Extract simulator.scope.countries (comma-separated ISO3 list) from a
# generated config — drives both the --array=0-N-1 size and the --countries
# list each array task indexes into via SLURM_ARRAY_TASK_ID.
countries_for_config() {
  sed -n 's/^[[:space:]]*countries:[[:space:]]*\(.*\)$/\1/p' "$1" | tr -d ' \r' | head -1
}

submitted=0
skipped=0
ARRAY_IDS=()    # array job ids, printed as FORECAST_JOB_IDS= for downstream afterok chains
for arm in "${ARMS[@]}"; do
  prefix="forecast"
  [[ "${arm}" == "counterfactual" ]] && prefix="counterfactual"
  # One config per unit: an ENSO group (enso_<g>) or, under partition_by:
  # registry, a partition-registry label of any name.
  while IFS= read -r cfg; do
    [[ -n "${cfg}" ]] || continue
    stem="$(basename "${cfg}" .yml)"      # e.g. forecast_enso_east
                                          #   or forecast_hmmCE_attrib
                                          #   or counterfactual_enso_east_set-nina
    group="${stem#${prefix}_}"            # -> enso_east[_set-nina] / hmmCE_attrib
    group="${group%_attrib}"              # attribution factual: forecast_<unit>_attrib

    # Counterfactual configs carry a trailing _set-<analog_set> token. It MUST be
    # stripped before the --groups filter: otherwise short= "east_set-nina", nothing
    # ever equals the requested "east", and the loop silently submits zero jobs.
    aset=""
    case "${group}" in
      *_set-*)
        aset="${group##*_set-}"           # -> nina
        group="${group%_set-*}"           # -> enso_east
        ;;
    esac
    short="${group#enso_}"                # -> east

    if [[ -n "${GROUPS_FILTER}" ]]; then
      match=0
      IFS=',' read -ra wanted <<< "${GROUPS_FILTER}"
      for w in "${wanted[@]}"; do
        [[ "${short}" == "${w}" ]] && match=1
      done
      [[ ${match} -eq 1 ]] || continue
    fi

    if [[ -n "${SETS_FILTER}" && -n "${aset}" ]]; then
      match=0
      IFS=',' read -ra wanted_sets <<< "${SETS_FILTER}"
      for w in "${wanted_sets[@]}"; do
        [[ "${aset}" == "${w}" ]] && match=1
      done
      [[ ${match} -eq 1 ]] || continue
    fi

    countries="$(countries_for_config "${cfg}")"
    if [[ -z "${countries}" ]]; then
      echo "ERROR: could not extract simulator.scope.countries from ${cfg}" >&2
      exit 1
    fi
    IFS=',' read -ra country_arr <<< "${countries}"
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

    arm_tag="fo"
    [[ "${arm}" == "attrib" ]] && arm_tag="fa"
    [[ "${arm}" == "counterfactual" ]] && arm_tag="cf"
    jobname="spfc_${arm_tag}_${short}"
    # Keep job names distinguishable across analog sets (abbreviated to stay within
    # SLURM's practical job-name width).
    [[ -n "${aset}" ]] && jobname="${jobname}_${aset:0:3}"

    # Skip-complete guard (see header).
    if [[ ${OVERWRITE} -eq 0 ]]; then
      mapfile -t _pd < <(parts_dir_for_config "${cfg}")
      _parts="${_pd[0]%/}/sim_spatpomp/${_pd[1]}/parts"
      _min_tu=1
      [[ "${PARTITION_BY}" == "enso_region" || "${PARTITION_BY}" == "registry" ]] || _min_tu=${n_countries}
      if [[ -n "${_pd[1]}" ]] && run_complete "${_parts}" "${_min_tu}"; then
        echo "[skip] ${jobname}: run complete (${_parts}); --overwrite resubmits it  <-  ${cfg}"
        skipped=$((skipped + 1))
        continue
      fi
    fi

    if [[ ${LIST_ONLY} -eq 1 ]]; then
      echo "[list] ${jobname}  --array=${array_spec}  (${scope_desc})  <-  ${cfg}"
      continue
    fi

    echo "--- Submitting ${jobname} --array=${array_spec} (${scope_desc}) (${cfg}) ..."
    jid=$(cholera_sbatch SHORT --parsable \
      --job-name="${jobname}" \
      --array="${array_spec}${CONCURRENCY:+%${CONCURRENCY}}" \
      analysis/bash/slurm/submit_forecast_groups.sh \
      --config "${cfg}" --countries "${countries}" ${EXTRA})
    echo "    array job id: ${jid}"
    ARRAY_IDS+=("${jid}")
    submitted=$((submitted + 1))

    # --- Optional: chain THIS config's postprocess + figures -----------------
    # One PER-CONFIG job (07c -C cfg --aggregate --engine spatpomp merges this
    # run's partitions, 07d -C cfg --engine spatpomp draws its per-run figures),
    # gated on just this config's array job via afterok. Skipped under --dry-run
    # (the array tasks only validate+exit, so there would be nothing to merge).
    if [[ ${POSTPROCESS} -eq 1 && ${DRYRUN} -eq 0 ]]; then
      pp_jid=$(cholera_sbatch SHORT --parsable --dependency="afterok:${jid}" \
        --job-name="fc_pp_${jobname}" \
        "$PP_SUBMITTER" "${cfg}")
      echo "    postprocess job id: ${pp_jid}  (afterok:${jid}, cfg=${cfg})"
    fi
  done < <(configs_for_arm "${arm}")
done

if [[ ${POSTPROCESS} -eq 1 && ${DRYRUN} -eq 1 ]]; then
  echo "--- (postprocess skipped: --dry-run submits no real work to merge)"
fi

[[ ${LIST_ONLY} -eq 1 ]] && exit 0

echo ""
echo "Submitted ${submitted} array job(s); skipped ${skipped} complete run(s)."
echo "Monitor with: squeue -u \$USER --format='%.10i %.14j %.8T %.10M %R'"
if [[ ${#ARRAY_IDS[@]} -gt 0 ]]; then
  echo "FORECAST_JOB_IDS=$(IFS=:; echo "${ARRAY_IDS[*]}")"
fi
