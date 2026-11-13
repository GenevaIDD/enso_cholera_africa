#!/bin/bash
# slurm/submit_population.sh — SBATCH script for one population array task.
#
# Called by SLURM; each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a
# task_NNNN.yml via `ls | sort`, then invokes 02b_preprocess_population_run.R --task.
#
# Args: $1 = config path, $2 = tasks_dir (absolute path to task YAML dir).
#
# Submitted by analysis/bash/stages/02_population.sh. Do not invoke directly.
#
# --mem=16G, NOT the house 8G: raster::raster() on a 100 m constrained WorldPop
# raster is 500 MB-1 GB resident, and NGA / COD are the worst cases. Benchmark
# one of those as a single task before trusting this number.
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH -c 1

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env GEO

echo "***** Beginning of script *****"
date

# The work is I/O- and GEOS-bound on one core; stop BLAS spawning threads it
# cannot use.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="$1"
TASKS_DIR="$2"

if [[ -z "$CONFIG" || -z "$TASKS_DIR" ]]; then
  echo "Error: usage: slurm/submit_population.sh <config.yml> <tasks_dir>" >&2
  exit 1
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "Error: config not found: $CONFIG" >&2
  exit 1
fi
if [[ ! -d "$TASKS_DIR" ]]; then
  echo "Error: tasks directory not found: $TASKS_DIR" >&2
  exit 1
fi

# Map SLURM_ARRAY_TASK_ID (0-indexed) to a task YAML via ls | sort.
TASK_FILE=$(ls "${TASKS_DIR}"/task_*.yml | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [[ -z "$TASK_FILE" ]]; then
  echo "Error: no task YAML for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${TASKS_DIR}" >&2
  exit 1
fi

echo "Config    : ${CONFIG}"
echo "Tasks dir : ${TASKS_DIR}"
echo "Task file : ${TASK_FILE}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"

Rscript analysis/02b_preprocess_population_run.R -C "${CONFIG}" --task "${TASK_FILE}" || exit 1

echo "***** End of script *****"
date
