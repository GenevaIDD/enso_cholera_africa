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
# Partition: the cluster profile's BIGMEM class (lib/cluster.sh); list every
# partition that can fit the job: SLURM starts it on the first with free
# resources. Each must allow the --time below.
#SBATCH --time=08:00:00
# Memory: dominated by the particle store that spatPomp's bpfilter keeps to
# build the smoothed trajectory (filter_traj=TRUE, spatpomp_core.R::
# run_one_bpfilter_replicate(); save_states is NOT set on this path): every
# resampled particle of every month and block, 11+ states × U × Np × M months
# × 8 B, plus the ancestry vectors. Sized for country-partitioned runs
# (partition_by=country) at Np=10000: the largest country (NGA, U=121, M=192)
# ≈ 17 GB raw, ~28 GB peak once the fits/draws/panel (~6 GB) and extraction
# copies are included → 35G leaves ~20% margin.
# NB: full-TU runs (U≤319) or larger Np need much more (the deployed hmmCE unit
# at Np=10000, ~213 months: ~90 GB peak) — bump --mem accordingly for those.
# simulator.filter.traj_storage: pruned | pruned_float32 (config_defaults.yml)
# keeps only the particles that are still ancestors of the current population
# and cuts this store by an order of magnitude or more; it is the default
# (deployed hmmCE unit at Np 30000, pruned_float32: 32 GB peak). The figures above
# are for traj_storage: full.
# 75G covers the pruned store on registry (transmission-unit) partitions, the
# largest being the hmmCE unit (16 countries; ~40-45 GB at Np 30000),
# and fits ~94 GB nodes. A traj_storage: full run needs 120G.
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
