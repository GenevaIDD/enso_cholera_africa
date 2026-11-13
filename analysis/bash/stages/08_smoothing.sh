#!/bin/bash
# stages/08_smoothing.sh — submit the spatPomp historical-smoothing pipeline
# as a SLURM ARRAY job, then chain the compile (reduction) job.
#
# Each array task runs one
# (TU, replicate-chunk) via 08b_smoothing_spatpomp.R --task; 08c_smoothing_compile.R
# recombines the partitions into importance-weighted ribbons + figures.
#
# Usage:
#   bash analysis/bash/stages/08_smoothing.sh --config <cfg.yml>
#   bash analysis/bash/stages/08_smoothing.sh --config <cfg.yml> --n_reps 100 --Np 25000 --chunk_size 10
#   bash analysis/bash/stages/08_smoothing.sh --config <cfg.yml> --resume      # only missing tasks
#   bash analysis/bash/stages/08_smoothing.sh --config <cfg.yml> --no-compile  # skip the chained compile
#
# Steps:
#   1. Build task YAMLs via 08a_smoothing_tasks.R (idempotent).
#   2. Derive --array spec (--resume scans parts/ for missing _SUCCESS sentinels).
#   3. sbatch the runner array; capture its job id.
#   4. sbatch the compile job with --dependency=afterok:<arrayid> (unless --no-compile).

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

CONFIG=""
RESUME=false
COMPILE=true
ARRAY_CONCURRENCY=1000   # array cap: 1000 = no practical limit
EXTRA=""          # forwarded to make_tasks (n_reps/Np/workers/chunk_size/tus)
COMPILE_EXTRA=""  # forwarded to compile step (e.g. --no_isweight)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)     CONFIG="$2";                  shift 2 ;;
    --resume)        RESUME=true;                  shift   ;;
    --no-compile)    COMPILE=false;                shift   ;;
    --concurrency)   ARRAY_CONCURRENCY="$2";       shift 2 ;;
    --no_isweight)   COMPILE_EXTRA="$COMPILE_EXTRA --no_isweight"; shift ;;
    --n_reps|--Np|--workers|--chunk_size|--tus)
                     EXTRA="$EXTRA $1 $2";         shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --config <cfg.yml> [--n_reps N] [--Np N] [--workers N] [--chunk_size N] [--tus L] [--resume] [--no-compile] [--no_isweight] [--concurrency K]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then echo "Error: --config is required." >&2; exit 1; fi
if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG" >&2; exit 1; fi

JOB_NAME="smoothing_$(basename "$CONFIG" .yml)"
mkdir -p "$CHOLERA_LOG_DIR"

# ---------------------------------------------------------------------------
# 1. Build task YAMLs.
# ---------------------------------------------------------------------------
echo "--- Building task YAMLs (08a_smoothing_tasks.R) ..."
MK_OUTPUT=$(Rscript analysis/08a_smoothing_tasks.R -C "$CONFIG" ${EXTRA} 2>/dev/null)

TASKS_DIR=$(echo "$MK_OUTPUT" | grep '^TASKS_DIR=' | cut -d= -f2-)
N_TASKS=$(  echo "$MK_OUTPUT" | grep '^N_TASKS='   | cut -d= -f2-)

if [[ -z "$TASKS_DIR" || -z "$N_TASKS" ]]; then
  echo "Error: 08a_smoothing_tasks.R did not emit TASKS_DIR / N_TASKS." >&2
  echo "       (captured: ${MK_OUTPUT})" >&2
  exit 1
fi
echo "Tasks dir : $TASKS_DIR"
echo "Task count: $N_TASKS"
if [[ "$N_TASKS" -eq 0 ]]; then echo "Error: zero tasks generated." >&2; exit 1; fi

# ---------------------------------------------------------------------------
# 2. Determine array spec.
# ---------------------------------------------------------------------------
if [[ "$RESUME" == "true" ]]; then
  echo "--- Resume mode: scanning for missing _SUCCESS sentinels ..."
  MISSING_IDS=$(Rscript --vanilla - "$TASKS_DIR" 2>/dev/null << 'RSCRIPT'
suppressPackageStartupMessages(library(yaml))
args      <- commandArgs(trailingOnly = TRUE)
tasks_dir <- args[1]
files     <- sort(list.files(tasks_dir, pattern = "^task_[0-9]+\\.yml$", full.names = TRUE))
missing   <- integer(0)
for (i in seq_along(files)) {
  task     <- yaml::read_yaml(files[[i]])
  sentinel <- file.path(task$output_part, "_SUCCESS")
  if (!file.exists(sentinel)) missing <- c(missing, i - 1L)
}
if (!length(missing)) cat("NONE\n") else cat(paste(missing, collapse = ","), "\n")
RSCRIPT
)
  MISSING_IDS=$(echo "$MISSING_IDS" | tail -1 | tr -d '[:space:]')
  if [[ "$MISSING_IDS" == "NONE" ]]; then
    echo "All $N_TASKS task(s) already completed (_SUCCESS found). Nothing to submit."
    ARRAY_SPEC=""
  else
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
  fi
else
  ARRAY_SPEC="0-$((N_TASKS - 1))"
fi

# ---------------------------------------------------------------------------
# 3. Submit runner array (capture job id for compile dependency).
# ---------------------------------------------------------------------------
ARRAY_JOB_ID=""
if [[ -n "$ARRAY_SPEC" ]]; then
  echo "--- Submitting array ${ARRAY_SPEC}%${ARRAY_CONCURRENCY} as job '${JOB_NAME}' ..."
  ARRAY_JOB_ID=$(cholera_sbatch BIGMEM --parsable \
    --array="${ARRAY_SPEC}%${ARRAY_CONCURRENCY}" \
    --job-name="${JOB_NAME}" \
    analysis/bash/slurm/submit_smoothing.sh \
    "$CONFIG" "$TASKS_DIR")
  echo "Runner array job id: ${ARRAY_JOB_ID}"
fi

# ---------------------------------------------------------------------------
# 4. Chain the compile (reduction) job.
# ---------------------------------------------------------------------------
if [[ "$COMPILE" == "true" ]]; then
  DEP_ARG=""
  [[ -n "$ARRAY_JOB_ID" ]] && DEP_ARG="--dependency=afterok:${ARRAY_JOB_ID}"
  echo "--- Submitting compile job (${DEP_ARG:-no dependency}) ..."
  COMPILE_JOB_ID=$(cholera_sbatch BIGMEM --parsable ${DEP_ARG} \
    --job-name="${JOB_NAME}_compile" \
    analysis/bash/slurm/submit_smoothing_compile.sh \
    "$CONFIG" ${COMPILE_EXTRA})
  echo "Compile job id: ${COMPILE_JOB_ID}"
fi
