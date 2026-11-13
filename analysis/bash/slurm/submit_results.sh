#!/bin/bash
# slurm/submit_results.sh — SBATCH script for ONE manuscript-results
# step (phase 09b) of a SLURM array.
#
# Called by SLURM; each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a
# task_NNNN.yml written by 09a_postprocess_results_tasks.R via `ls | sort`, then
# runs that one step with 09b_postprocess_results.R --task. The allocated cores
# are passed as --workers (forked per-group fit loads, per-run hindcast scoring;
# artefacts are identical for any value).
#
# Args: $1 = results config path, $2 = tasks_dir.
# Submitted by analysis/bash/stages/09_postprocess_results.sh (one array per
# dependency wave). Do not invoke directly.
#
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --time=06:00:00
# Memory: sized for the heaviest steps (enso_relation / track_params load one
# brms fit per ENSO group; smoothing / report_smoothing read every 08b part.qs;
# hindcast scores every manifest run — the hindcast postprocess job uses 64G).
# The stage driver's --mem / --cpus / --time override these per submission.
#SBATCH --mem=64G
#SBATCH -c 4

# R 4.3.2 + GDAL stack (as OutbreakExtractR): 09 needs sf, absent from R 4.2.1.
source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env --purge GEO

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="$1"
TASKS_DIR="$2"

if [[ -z "$CONFIG" || -z "$TASKS_DIR" ]]; then
  echo "Error: usage: slurm/submit_results.sh <config.yml> <tasks_dir>" >&2
  exit 1
fi
if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG" >&2; exit 1; fi
if [[ ! -d "$TASKS_DIR" ]]; then echo "Error: tasks dir not found: $TASKS_DIR" >&2; exit 1; fi

# Map SLURM_ARRAY_TASK_ID (0-indexed) to a task YAML via ls | sort.
TASK_FILE=$(ls "${TASKS_DIR}"/task_*.yml | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
if [[ -z "$TASK_FILE" ]]; then
  echo "Error: no task YAML for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${TASKS_DIR}" >&2
  exit 1
fi
WORKERS="${SLURM_CPUS_PER_TASK:-1}"

echo "Config    : ${CONFIG}"
echo "Tasks dir : ${TASKS_DIR}"
echo "Task file : ${TASK_FILE}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Workers   : ${WORKERS}"

Rscript analysis/09b_postprocess_results.R -C "${CONFIG}" --task "${TASK_FILE}" --workers "${WORKERS}" || exit 1

echo "***** End of script *****"
date
