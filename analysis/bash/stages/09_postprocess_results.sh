#!/bin/bash
# stages/09_postprocess_results.sh — derive the manuscript-results artefacts
# (phase 09) as SLURM arrays, one array per dependency wave.
#
# 09a_postprocess_results_tasks.R writes one task YAML per step of the step
# registry (utils/results_derive.R) and assigns waves: steps of one wave never
# read each other's artefacts. Wave k + 1 waits afterok on wave k, so no two
# jobs derive the same artefact. Each submission has its own task folder
# (tasks/sub-<timestamp>-<pid>/), so overlapping submissions are safe.
#
# Run-agnostic: the results config names the run (results.deploy_config, set
# by deploy_prod_run.R --stage align_results). No --resume: fresh steps are
# no-ops, so resubmitting re-derives only what is stale; --steps submits a
# subset (waves renumbered from 1).
#
# Usage:
#   bash analysis/bash/stages/09_postprocess_results.sh --config analysis/configs/results/manuscript.yml
#   bash analysis/bash/stages/09_postprocess_results.sh --config <cfg> --steps hindcast,oni
#   bash analysis/bash/stages/09_postprocess_results.sh --config <cfg> --redo     # ignore freshness
#   bash analysis/bash/stages/09_postprocess_results.sh --config <cfg> --mem 96G --cpus 8 --time 12:00:00
#   bash analysis/bash/stages/09_postprocess_results.sh --config <cfg> --list     # steps, waves, status
#   bash analysis/bash/stages/09_postprocess_results.sh --config <cfg> --dry-run  # print sbatch only
#
# --partition P overrides the cluster profile's (CHOLERA_RES_RESULTS, else its
# SHORT class); the attribution steps need ~100 GB (run_pipeline.sh asks for
# 160G), so on a cluster with smaller nodes pass partitions with bigger ones.
# --after IDS (':' or ','): wave 1 waits afterok on upstream jobs, e.g. the
# forecast arrays of stages/07_forecast_groups.sh.
# Prints RESULTS_JOB_IDS=<wave1 id>:<wave2 id>:... (afterok-ready). Then,
# locally: 10_manuscript_figures.R and quarto render of 11.

set -euo pipefail

# The GEO stack (R + sf / GDAL, as slurm/submit_population.sh):
# 09 needs sf (figures_cb.R maps, taxdat basemaps), absent from the R 4.2.1 stack.
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env --purge GEO

CONFIG=""
STEPS=""
REDO=false
LIST=false
DRY_RUN=false
ARRAY_CONCURRENCY=1000   # array cap: 1000 = no practical limit
MEM=""
CPUS=""
TIME=""
PARTITION=""
AFTER=""

usage() {
  echo "Usage: $0 --config <results cfg.yml> [--steps id,id] [--redo] [--concurrency K] [--mem M] [--cpus N] [--time T] [--partition P] [--after ID:ID] [--list] [--dry-run]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)   CONFIG="$2";            shift 2 ;;
    --steps)       STEPS="$2";             shift 2 ;;
    --redo)        REDO=true;              shift   ;;
    --concurrency) ARRAY_CONCURRENCY="$2"; shift 2 ;;
    --mem)         MEM="$2";               shift 2 ;;
    --cpus)        CPUS="$2";              shift 2 ;;
    --time)        TIME="$2";              shift 2 ;;
    --partition)   PARTITION="$2";         shift 2 ;;
    --after)       AFTER="${2//,/:}";      shift 2 ;;
    --list)        LIST=true;              shift   ;;
    --dry-run)     DRY_RUN=true;           shift   ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then echo "Error: --config is required." >&2; usage; exit 1; fi
if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG" >&2; exit 1; fi
AFTER="${AFTER#:}"; AFTER="${AFTER%:}"
if [[ -n "$AFTER" && ! "$AFTER" =~ ^[^:]+(:[^:]+)*$ ]]; then
  echo "Error: --after expects job ids separated by ':' or ',' (got '$AFTER')." >&2; exit 1
fi

TASK_ARGS=()
[[ -n "$STEPS" ]]        && TASK_ARGS+=(--steps "$STEPS")
[[ "$REDO" == "true" ]]  && TASK_ARGS+=(--redo_cache TRUE)

if [[ "$LIST" == "true" ]]; then
  Rscript analysis/09a_postprocess_results_tasks.R -C "$CONFIG" ${TASK_ARGS[@]+"${TASK_ARGS[@]}"} --list
  exit 0
fi

JOB_NAME="results_$(basename "$CONFIG" .yml)"
mkdir -p "$CHOLERA_LOG_DIR"

# ---------------------------------------------------------------------------
# 1. Build task YAMLs (one per step, each with its wave).
# ---------------------------------------------------------------------------
echo "--- Building task YAMLs (09a_postprocess_results_tasks.R) ..." >&2
# `|| true` / out_key(): under `set -euo pipefail` a failing 09a or a missing key
# would otherwise exit here silently, before the error message below.
MK_OUTPUT=$(Rscript analysis/09a_postprocess_results_tasks.R -C "$CONFIG" ${TASK_ARGS[@]+"${TASK_ARGS[@]}"} 2>/dev/null) || true
out_key() { echo "$MK_OUTPUT" | grep "^$1=" | cut -d= -f2- || true; }

TASKS_DIR=$(out_key TASKS_DIR)
N_TASKS=$(out_key N_TASKS)
N_WAVES=$(out_key N_WAVES)

if [[ -z "$TASKS_DIR" || -z "$N_TASKS" || -z "$N_WAVES" ]]; then
  echo "Error: 09a_postprocess_results_tasks.R did not emit TASKS_DIR / N_TASKS / N_WAVES." >&2
  echo "       Re-run it by hand to see why: Rscript analysis/09a_postprocess_results_tasks.R -C $CONFIG" >&2
  exit 1
fi
if [[ "$N_TASKS" -eq 0 ]]; then echo "Error: zero tasks generated." >&2; exit 1; fi
echo "Tasks dir : $TASKS_DIR" >&2
echo "Task count: $N_TASKS in $N_WAVES wave(s)" >&2

# ---------------------------------------------------------------------------
# 2. One array per wave, each chained afterok on the previous one.
# ---------------------------------------------------------------------------
SBATCH_OVERRIDES=()
[[ -n "$MEM"  ]] && SBATCH_OVERRIDES+=(--mem="$MEM")
[[ -n "$CPUS" ]] && SBATCH_OVERRIDES+=(-c "$CPUS")
[[ -n "$TIME" ]] && SBATCH_OVERRIDES+=(--time="$TIME")
[[ -n "$PARTITION" ]] && SBATCH_OVERRIDES+=(--partition="$PARTITION")

# Wave 1 waits on --after (upstream jobs); wave k + 1 on wave k.
PREV_JOB_ID="$AFTER"
JOB_IDS=()
for (( k = 1; k <= N_WAVES; k++ )); do
  IDS=$(out_key "WAVE_${k}")
  if [[ -z "$IDS" ]]; then echo "Error: no WAVE_${k}= line from 09a." >&2; exit 1; fi
  DEP_ARGS=()
  [[ -n "$PREV_JOB_ID" ]] && DEP_ARGS+=(--dependency="afterok:${PREV_JOB_ID}")
  cholera_sbatch_cmd SHORT --parsable ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} ${SBATCH_OVERRIDES[@]+"${SBATCH_OVERRIDES[@]}"} \
    --array="${IDS}%${ARRAY_CONCURRENCY}" --job-name="${JOB_NAME}_w${k}" \
    analysis/bash/slurm/submit_results.sh "$CONFIG" "$TASKS_DIR"
  CMD=("${CHOLERA_SBATCH_CMD[@]}")
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[dry-run] wave ${k}: ${CMD[*]}" >&2
    PREV_JOB_ID="<wave${k}>"
  else
    echo "--- Submitting wave ${k} (tasks ${IDS}) ${DEP_ARGS[*]:-} ..." >&2
    PREV_JOB_ID=$("${CMD[@]}")
    echo "Wave ${k} array job id: ${PREV_JOB_ID}" >&2
  fi
  JOB_IDS+=("$PREV_JOB_ID")
done

echo "RESULTS_JOB_IDS=$(IFS=:; echo "${JOB_IDS[*]}")"
