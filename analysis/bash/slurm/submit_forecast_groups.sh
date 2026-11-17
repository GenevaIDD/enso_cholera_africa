#!/bin/bash
# slurm/submit_forecast_groups.sh — SBATCH array payload for one per-ENSO-group
# forecast or counterfactual spatPomp run (07b_simulator_spatpomp.R): one array
# task per country, replicate chunks run as background Rscript processes.
#
# Country mode (filter.partition_by = country): --only_country restricts each
# task to one TU; --replicate_offset makes replicate ids unique per chunk and
# nests output under parts/tu=<ISO3>/repoff=<offset>/ (07b's
# aggregate_partitions() globs the _SUCCESS sentinels recursively).
# Region mode (--partition_by enso_region): the whole group is one partition,
# a single task (index 0), no --only_country; output under
# parts/tu=<region>/repoff=<offset>/.
#
# Submitted by stages/07_forecast_groups.sh (computes --array and
# --countries from simulator.scope.countries), or directly:
#   sbatch --job-name=spfc_fo_east --array=0-5 \
#          analysis/bash/slurm/submit_forecast_groups.sh \
#          --config analysis/configs/forecast_country_by_group/forecast_enso_east.yml \
#          --countries BDI,COD,KEN,RWA,TZA,UGA \
#          [--n_replicates 1000] [--chunk_size 50] [--Np 10000] [--dry-run] [--overwrite]
#
# Replicates: --n_replicates, else the config's simulator.n_replicates, else
# 100. --chunk_size splits them into ceil(n_replicates / chunk_size) chunks run
# concurrently in the task [default: one serial chunk]; size it with -c / --mem
# (each chunk is a single-threaded Rscript holding one country's model, draws
# and panel plus the smoothing cache).
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH -c 4

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date
hostname
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-<none>}  SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-<none>}"

# Disable BLAS thread proliferation — spatPomp/pomp use their own threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG=""
COUNTRIES=""
N_REPLICATES=""
CHUNK_SIZE=""
PARTITION_BY=""    # ""/country = per-country array task; enso_region|registry = whole-unit single task
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)       CONFIG="$2";       shift 2 ;;
    --countries)    COUNTRIES="$2";    shift 2 ;;
    --n_replicates) N_REPLICATES="$2"; shift 2 ;;
    --chunk_size)   CHUNK_SIZE="$2";   shift 2 ;;
    --partition_by) PARTITION_BY="$2"; shift 2 ;;
    --Np)           EXTRA_ARGS="$EXTRA_ARGS --Np $2";     shift 2 ;;
    --dry-run)      EXTRA_ARGS="$EXTRA_ARGS --dry-run";   shift ;;
    --overwrite)    EXTRA_ARGS="$EXTRA_ARGS --overwrite"; shift ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1 ;;
  esac
done

if [[ -z "${CONFIG}" ]]; then
  echo "ERROR: --config is required" >&2
  exit 1
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: must be submitted as a SLURM array job (sbatch --array=0-N-1 ...) — SLURM_ARRAY_TASK_ID is unset." >&2
  exit 1
fi

# Region mode (filter.partition_by=enso_region): the whole ENSO group is ONE
# forward partition, so this is a single array task (index 0) that runs the
# entire region in one process — no --only_country, no per-country fan-out. 07b
# derives the region TU from the frame and assembles its handoff IC from the
# per-country smoothing caches. Replicate-chunk fan-out still applies.
# Country mode (default): one array task per ISO3, indexed by SLURM_ARRAY_TASK_ID.
# Registry mode (filter.partition_by=registry) is the same single-task path: a
# unit of the partition registry couples its member countries (possibly of
# several fit groups) in one process.
if [[ "${PARTITION_BY}" == "enso_region" || "${PARTITION_BY}" == "registry" ]]; then
  ONLY_COUNTRY_ARG=""
  LABEL="$(basename "${CONFIG}" .yml)"
  # NB: do NOT forward --partition_by to 07b_simulator_spatpomp.R — it has no such
  # CLI flag and reads filter.partition_by from the config's `filter:` block (the
  # generators write `partition_by: enso_region` into every region config). The
  # --partition_by seen here is a bash-side switch only: it gates --only_country
  # (dropped in region mode) and the log label, nothing forwarded to the Rscript.
  echo "config       : ${CONFIG}"
  echo "partition_by : ${PARTITION_BY}  (whole unit, single task; array index ${SLURM_ARRAY_TASK_ID})"
else
  if [[ -z "${COUNTRIES}" ]]; then
    echo "ERROR: --countries is required (comma-separated ISO3 list, indexed by SLURM_ARRAY_TASK_ID)" >&2
    exit 1
  fi
  IFS=',' read -ra COUNTRY_ARR <<< "${COUNTRIES}"
  ISO3="${COUNTRY_ARR[$SLURM_ARRAY_TASK_ID]:-}"
  if [[ -z "${ISO3}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range for --countries '${COUNTRIES}' (${#COUNTRY_ARR[@]} countries)." >&2
    exit 1
  fi
  ONLY_COUNTRY_ARG="--only_country ${ISO3}"
  LABEL="${ISO3}"
  echo "config       : ${CONFIG}"
  echo "country      : ${ISO3}  (array index ${SLURM_ARRAY_TASK_ID} of ${#COUNTRY_ARR[@]}: ${COUNTRIES})"
fi

# Replicate count: the flag, else the config's simulator.n_replicates, else the
# runtime default 100 (see the header).
NREP_SOURCE="--n_replicates"
if [[ -z "${N_REPLICATES}" ]]; then
  N_REPLICATES=$(Rscript -e 'x <- yaml::read_yaml(commandArgs(TRUE)[1])$simulator$n_replicates; if (!is.null(x)) cat(as.integer(x))' "${CONFIG}" 2>/dev/null || true)
  NREP_SOURCE="config simulator.n_replicates"
  if [[ -z "${N_REPLICATES}" ]]; then
    N_REPLICATES=100
    NREP_SOURCE="runtime default (the config has no simulator.n_replicates)"
  fi
fi
if ! [[ "${N_REPLICATES}" =~ ^[0-9]+$ ]] || (( N_REPLICATES < 1 )); then
  echo "ERROR: n_replicates '${N_REPLICATES}' (from ${NREP_SOURCE}) is not a positive integer." >&2
  exit 1
fi
CHUNK_SIZE="${CHUNK_SIZE:-${N_REPLICATES}}"   # default: one chunk == no fan-out

echo "n_replicates : ${N_REPLICATES}   chunk_size: ${CHUNK_SIZE}   (from ${NREP_SOURCE})"
echo "extra        : ${EXTRA_ARGS}"

mkdir -p logs

pids=()
chunk_logs=()
offset=0
chunk_id=0
while (( offset < N_REPLICATES )); do
  remaining=$(( N_REPLICATES - offset ))
  this_chunk=$(( CHUNK_SIZE < remaining ? CHUNK_SIZE : remaining ))
  chunk_id=$(( chunk_id + 1 ))
  chunk_log="${CHOLERA_LOG_DIR}/${SLURM_JOB_NAME:-spfc}_${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID}_${LABEL}_chunk${chunk_id}.log"
  chunk_logs+=("${chunk_log}")
  echo "--- launching chunk ${chunk_id} for ${LABEL}: replicate_offset=${offset} n_replicates=${this_chunk} -> ${chunk_log}"
  Rscript analysis/07b_simulator_spatpomp.R \
    -C "${CONFIG}" \
    ${ONLY_COUNTRY_ARG} \
    --n_replicates "${this_chunk}" \
    --replicate_offset "${offset}" \
    --n_replicates_total "${N_REPLICATES}" \
    ${EXTRA_ARGS} \
    > "${chunk_log}" 2>&1 &
  pids+=($!)
  offset=$(( offset + this_chunk ))
done

echo "--- waiting on ${#pids[@]} replicate-chunk process(es) for ${LABEL}..."
fail=0
for i in "${!pids[@]}"; do
  if ! wait "${pids[$i]}"; then
    echo "ERROR: chunk $((i + 1)) for ${LABEL} failed — see ${chunk_logs[$i]}" >&2
    fail=1
  fi
done

if [[ ${fail} -ne 0 ]]; then
  echo "ERROR: at least one replicate-chunk process failed for ${LABEL}." >&2
  exit 1
fi

echo "***** End of script *****"
date
