#!/bin/bash
# slurm/submit_population_country.sh — SBATCH script for one PER-COUNTRY
# population array task (sibling of submit_population.sh, kept side by side).
#
# Each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a task_NNNN.yml via
# `ls | sort` and runs 02b_preprocess_population_run_country.R --task: one
# country, every year, geometry built once, the years extracted in parallel on
# the task's CPUs (forked workers, one raster each).
#
# Args: $1 = config path, $2 = tasks_dir (absolute path to task YAML dir).
#
# Submitted by analysis/bash/stages/02_population.sh --task-unit country. Do
# not invoke directly.
#
# -c 8 / --mem=64G: up to 8 years at once. The serial tasks of the Nov_2026 run
# peaked at 6.5 GB RSS with one raster resident, so 8 workers stay
# under 64 GB; COD took ~9.5 min per year, i.e. ~2.5 h serial, ~20 min on 8.
# --time=06:00:00: the slowest serial country (COD) took 2.4 h.
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH -c 8

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env GEO

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="$1"
TASKS_DIR="$2"

if [[ -z "$CONFIG" || -z "$TASKS_DIR" ]]; then
  echo "Error: usage: slurm/submit_population_country.sh <config.yml> <tasks_dir>" >&2
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

TASK_FILE=$(ls "${TASKS_DIR}"/task_*.yml | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [[ -z "$TASK_FILE" ]]; then
  echo "Error: no task YAML for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${TASKS_DIR}" >&2
  exit 1
fi

echo "Config    : ${CONFIG}"
echo "Tasks dir : ${TASKS_DIR}"
echo "Task file : ${TASK_FILE}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"

Rscript analysis/02b_preprocess_population_run_country.R -C "${CONFIG}" --task "${TASK_FILE}" || exit 1

echo "***** End of script *****"
date
