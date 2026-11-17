#!/bin/bash
# slurm/submit_smoothing.sh — SBATCH script for ONE smoothing array task.
#
# Called by SLURM; each array slot maps $SLURM_ARRAY_TASK_ID (0-indexed) to a
# task_NNNN.yml via `ls | sort`, then invokes 08b_smoothing_spatpomp.R --task.
# Each task is one (TU, replicate-chunk), so the per-task footprint is much
# smaller than a single-job run — request modest mem and 1 core.
#
# Args: $1 = config path, $2 = tasks_dir (absolute path to task YAML dir).
# Submitted by analysis/bash/stages/08_smoothing.sh. Do not invoke directly.
#
#SBATCH --output=logs/%x_%A_%a.log
# Partition: the cluster profile's BIGMEM class (lib/cluster.sh); SLURM starts
# the job on the first listed partition with free resources. Each must allow
# the --time below.
#SBATCH --time=08:00:00
# Memory is dominated by the particle store bpfilter keeps for the smoothed
# trajectory (filter_traj, spatpomp_core.R::run_one_bpfilter_replicate()):
# about states x U x Np x months x 8 B under traj_storage: full.
# simulator.filter.traj_storage pruned | pruned_float32 (the default) keeps only
# surviving ancestors and cuts it by an order of magnitude or more.
# 75G covers the pruned store on registry (transmission-unit) partitions;
# traj_storage: full needs about 120G. Raise --mem for larger Np.
#SBATCH --mem=75G
#SBATCH -c 1

source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env

echo "***** Beginning of script *****"
date

# Disable BLAS thread proliferation — spatPomp/pomp use their own threading.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="$1"
TASKS_DIR="$2"

if [[ -z "$CONFIG" || -z "$TASKS_DIR" ]]; then
  echo "Error: usage: slurm/submit_smoothing.sh <config.yml> <tasks_dir>" >&2
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

echo "Config    : ${CONFIG}"
echo "Tasks dir : ${TASKS_DIR}"
echo "Task file : ${TASK_FILE}"
echo "Array ID  : ${SLURM_ARRAY_TASK_ID}"

Rscript analysis/08b_smoothing_spatpomp.R -C "${CONFIG}" --task "${TASK_FILE}" || exit 1

echo "***** End of script *****"
date
