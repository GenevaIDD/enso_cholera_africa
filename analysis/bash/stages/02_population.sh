#!/bin/bash
# stages/02_population.sh — submit the LP x gid population extraction as a
# SLURM array, or resume by re-submitting only missing tasks.
#
# Usage:
#   bash analysis/bash/stages/02_population.sh --config <config.yml>
#   bash analysis/bash/stages/02_population.sh --config <config.yml> --resume
#   bash analysis/bash/stages/02_population.sh --config <config.yml> --task-unit country_year --countries BDI --years 2020
#   bash analysis/bash/stages/02_population.sh --config <config.yml> --concurrency 10
#   bash analysis/bash/stages/02_population.sh --config <config.yml> --after <jobid> --plan
#
# --task-unit country (default): one array task per country,
#   all years, 02b_preprocess_population_run_country.R (geometry once per
#   country, Q_g cache). country_year: the per-(country, year) array of
#   02b_preprocess_population_run.R. --years applies to country_year only.
# --after JOBID: the array waits for that job (afterok; campaigns/preprocess.sh
#   chains it after 00, which writes the LP shapefiles).
# --plan: print the steps and sbatch lines; builds no tasks, submits nothing.
#
# Steps:
#   0. REFUSE to submit unless the boundary-cache manifest is complete.
#   1. Build task YAMLs via Rscript analysis/02a_preprocess_population_tasks.R.
#   2. Count task files; derive --array spec.
#      --resume: re-submit only tasks with a missing or stale partition
#      (population_core.R::popPartFresh(): _SUCCESS built from the current LP
#      corpus); emit a sparse spec.
#   3. sbatch the array, then 02c_preprocess_population_compile.R afterok on it.
#
# Prints POP_ARRAY_JOB_ID= and POP_COMPILE_JOB_ID=.
#
# The layout of the other array stages (08_smoothing.sh), plus step 0.
#
# NOTE: do NOT run two arrays against the same raster_dir at
# once. The truncated-LZW delete-and-retry path in download_worldpop_constrained()
# (OutbreakExtractR R/get_pop.R:52-63) has no cross-worker lock, so a second
# array can delete a raster a first array is mid-read on. --resume is the
# supported way to fill gaps.

set -euo pipefail

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
# Same stack as the payloads (GEO: R + sf / GDAL): 02a and the --resume check run here,
# and the population modules need sf.
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

CONFIG=""
RESUME=false
COUNTRIES=""
YEARS=""
ARRAY_CONCURRENCY=1000   # matches config_defaults.yml population.slurm.array_concurrency
CACHE_DIR="data/country_boundaries"
SKIP_BOUNDARY_CHECK=false
TASK_UNIT="country"
AFTER=""
PLAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)            CONFIG="$2";            shift 2 ;;
    --resume)               RESUME=true;            shift   ;;
    --countries)            COUNTRIES="$2";         shift 2 ;;
    --years)                YEARS="$2";             shift 2 ;;
    --concurrency)          ARRAY_CONCURRENCY="$2"; shift 2 ;;
    --cache-dir)            CACHE_DIR="$2";         shift 2 ;;
    --skip-boundary-check)  SKIP_BOUNDARY_CHECK=true; shift ;;
    --task-unit)            TASK_UNIT="$2";         shift 2 ;;
    --after)                AFTER="$2";             shift 2 ;;
    --plan)                 PLAN=true;              shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --config <config.yml> [--task-unit country|country_year] [--resume] [--countries A,B] [--years Y1,Y2] [--concurrency K] [--after JOBID] [--plan]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "Error: --config is required." >&2
  exit 1
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "Error: config file not found: $CONFIG" >&2
  exit 1
fi

case "$TASK_UNIT" in
  country)      PAYLOAD=analysis/bash/slurm/submit_population_country.sh ;;
  country_year) PAYLOAD=analysis/bash/slurm/submit_population.sh ;;
  *) echo "Error: --task-unit must be country or country_year, got '$TASK_UNIT'." >&2; exit 1 ;;
esac
if [[ "$TASK_UNIT" == "country" && -n "$YEARS" ]]; then
  echo "Error: --years subsets per-(country, year) tasks; with --task-unit country every year is built." >&2
  exit 1
fi

JOB_NAME="pop_$(basename "$CONFIG" .yml)"
DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")
# 02c sizing: CHOLERA_RES_POP_COMPILE in the cluster profile overrides these defaults.
COMPILE_RES=(--res POP_COMPILE --defaults "--cpus-per-task=1 --mem=32G --time=04:00:00")

if [[ "$PLAN" == "true" ]]; then
  echo "PLAN: boundary-cache gate on ${CACHE_DIR}/_manifest.csv$([[ "$SKIP_BOUNDARY_CHECK" == "true" ]] && echo " (skipped)")" >&2
  echo "PLAN: Rscript analysis/02a_preprocess_population_tasks.R -C $CONFIG --task_unit $TASK_UNIT${COUNTRIES:+ --countries $COUNTRIES}${YEARS:+ --years $YEARS}" >&2
  cholera_sbatch_cmd LONG --parsable --array="<0..N-1>%${ARRAY_CONCURRENCY}" --job-name="${JOB_NAME}" \
    ${DEP[@]+"${DEP[@]}"} "$PAYLOAD" "$CONFIG" "<TASKS_DIR>"
  echo "PLAN: ${CHOLERA_SBATCH_CMD[*]}" >&2
  echo "POP_ARRAY_JOB_ID=PLAN_pop_array"
  cholera_sbatch_cmd LONG "${COMPILE_RES[@]}" --parsable --job-name="${JOB_NAME}_compile" \
    --dependency=afterok:PLAN_pop_array analysis/bash/slurm/submit_preprocess.sh \
    analysis/02c_preprocess_population_compile.R -C "$CONFIG"
  echo "PLAN: ${CHOLERA_SBATCH_CMD[*]}" >&2
  echo "POP_COMPILE_JOB_ID=PLAN_pop_compile"
  exit 0
fi

# ---------------------------------------------------------------------------
# 0. Boundary-cache gate.
# ---------------------------------------------------------------------------
# Compute nodes CAN reach the internet here, but a boundary fetch that fails on
# a worker does not fail the task — it silently disables the UN adjustment. The
# manifest is the only place that distinction is visible before the numbers are
# already wrong.
if [[ "$SKIP_BOUNDARY_CHECK" == "true" ]]; then
  echo "!!! --skip-boundary-check: UN adjustment may be silently skipped."
else
  MANIFEST="${CACHE_DIR}/_manifest.csv"
  if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: boundary manifest not found: $MANIFEST" >&2
    echo "       Build the national-boundary cache (${CACHE_DIR}) on a node with internet first." >&2
    exit 1
  fi
  # Column 3 is `ok`; R writes logicals as TRUE/FALSE.
  N_BAD=$(awk -F, 'NR>1 && $3 !~ /TRUE/ {n++} END {print n+0}' "$MANIFEST")
  if [[ "$N_BAD" -gt 0 ]]; then
    echo "Error: $N_BAD country/ies missing from the boundary cache." >&2
    awk -F, 'NR>1 && $3 !~ /TRUE/ {gsub(/"/,"",$1); print "       - " $1}' "$MANIFEST" >&2
    echo "       Their UN adjustment would be SILENTLY skipped. Rebuild the boundary cache" >&2
    echo "       (${CACHE_DIR}) on a node with internet first." >&2
    exit 1
  fi
  echo "Boundary cache: complete ($(($(wc -l < "$MANIFEST") - 1)) countries)."
fi

# ---------------------------------------------------------------------------
# 1. Build (or rebuild) task YAMLs.
# ---------------------------------------------------------------------------
echo "--- Building task YAMLs (02a_preprocess_population_tasks.R) ..."
MK_ARGS=(-C "$CONFIG" --task_unit "$TASK_UNIT")
[[ -n "$COUNTRIES" ]] && MK_ARGS+=(--countries "$COUNTRIES")
[[ -n "$YEARS"     ]] && MK_ARGS+=(--years "$YEARS")

MK_OUTPUT=$(Rscript analysis/02a_preprocess_population_tasks.R "${MK_ARGS[@]}" 2>/dev/null)

TASKS_DIR=$(echo "$MK_OUTPUT" | grep '^TASKS_DIR=' | cut -d= -f2-)
N_TASKS=$(  echo "$MK_OUTPUT" | grep '^N_TASKS='   | cut -d= -f2-)

if [[ -z "$TASKS_DIR" || -z "$N_TASKS" ]]; then
  echo "Error: 02a_preprocess_population_tasks.R did not emit TASKS_DIR / N_TASKS." >&2
  echo "       (captured stdout: ${MK_OUTPUT})" >&2
  exit 1
fi

echo "Tasks dir : $TASKS_DIR"
echo "Task count: $N_TASKS"

if [[ "$N_TASKS" -eq 0 ]]; then
  echo "Error: zero tasks generated — check config." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 2. Determine array spec.
# ---------------------------------------------------------------------------
if [[ "$RESUME" == "true" ]]; then
  echo "--- Resume mode: scanning for missing or stale partitions ..."
  MISSING_IDS=$(Rscript --vanilla - "$TASKS_DIR" 2>/dev/null << 'RSCRIPT'
suppressPackageStartupMessages(library(yaml))
source("analysis/utils/population_core.R")   # popPartFresh()
args      <- commandArgs(trailingOnly = TRUE)
tasks_dir <- args[1]
files     <- sort(list.files(tasks_dir, pattern = "^task_[0-9]+\\.yml$", full.names = TRUE))
missing   <- integer(0)
md5       <- list()
for (i in seq_along(files)) {
  task <- yaml::read_yaml(files[[i]])
  lp   <- task$lp_shapefiles
  if (is.null(md5[[lp]])) md5[[lp]] <- unname(tools::md5sum(path.expand(lp)))
  # A country task (part_dirs) is done only when every year is fresh.
  dirs <- unlist(task$part_dirs %||% task$part_dir)
  if (!all(vapply(dirs, popPartFresh, logical(1), lp_shapefiles = lp, lp_md5 = md5[[lp]])))
    missing <- c(missing, i - 1L)
}
if (!length(missing)) cat("NONE\n") else cat(paste(missing, collapse = ","), "\n")
RSCRIPT
)

  MISSING_IDS=$(echo "$MISSING_IDS" | tail -1 | tr -d '[:space:]')

  if [[ "$MISSING_IDS" == "NONE" ]]; then
    echo "All $N_TASKS task(s) already completed (fresh _SUCCESS found). Nothing to submit."
    echo "POP_ARRAY_JOB_ID="
    echo "POP_COMPILE_JOB_ID="
    exit 0
  fi

  # Condense the comma-separated list into a SLURM range spec (e.g. "3,7,12-15").
  ARRAY_SPEC=$(python3 -c "
ids = sorted(int(x) for x in '${MISSING_IDS}'.split(','))
ranges, i = [], 0
while i < len(ids):
    j = i
    while j + 1 < len(ids) and ids[j+1] == ids[j] + 1:
        j += 1
    ranges.append(str(ids[i]) if i == j else f'{ids[i]}-{ids[j]}')
    i = j + 1
print(','.join(ranges))
" 2>/dev/null || echo "$MISSING_IDS")

  echo "Missing task IDs (0-indexed): $ARRAY_SPEC"
else
  ARRAY_SPEC="0-$((N_TASKS - 1))"
fi

# ---------------------------------------------------------------------------
# 3. Submit.
# ---------------------------------------------------------------------------
mkdir -p "$CHOLERA_LOG_DIR"
echo "--- Submitting array ${ARRAY_SPEC}%${ARRAY_CONCURRENCY} as job '${JOB_NAME}' (${TASK_UNIT}) ..."
ARRAY_ID=$(cholera_sbatch LONG --parsable \
  --array="${ARRAY_SPEC}%${ARRAY_CONCURRENCY}" \
  --job-name="${JOB_NAME}" \
  ${DEP[@]+"${DEP[@]}"} \
  "$PAYLOAD" \
  "$CONFIG" "$TASKS_DIR")
echo "POP_ARRAY_JOB_ID=${ARRAY_ID}"

# 02c reads every partition: afterok on the whole array.
COMPILE_ID=$(cholera_sbatch LONG "${COMPILE_RES[@]}" --parsable \
  --job-name="${JOB_NAME}_compile" \
  --dependency="afterok:${ARRAY_ID}" \
  analysis/bash/slurm/submit_preprocess.sh \
  analysis/02c_preprocess_population_compile.R -C "$CONFIG")
echo "POP_COMPILE_JOB_ID=${COMPILE_ID}"
