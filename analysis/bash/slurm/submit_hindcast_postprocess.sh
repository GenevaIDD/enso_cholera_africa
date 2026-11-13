#!/bin/bash
# slurm/submit_hindcast_postprocess.sh — SBATCH script for the hindcast
# postprocess + figures step: once ALL per-(group, year) hindcast array jobs have
# finished, score every run in the manifest and draw the skill figures.
#
# Runs two whole-manifest orchestrators (NOT per-config — each internally
# iterates the per-(group, year) runs listed in the manifest):
#   1. 07c_simulator_postprocess.R --hindcast  — aggregates each run's partitions
#      (do_aggregate=TRUE) and writes the combined launch-year-tagged skill table.
#   2. 07d_simulator_figures.R --hindcast      — draws the hindcast skill figures.
#
# Submitted by analysis/bash/stages/07_hindcast.sh --postprocess (typically
# with --dependency=afterok:<all submitted array ids>). May also be run directly
# once the hindcast array jobs have finished:
#   sbatch analysis/bash/slurm/submit_hindcast_postprocess.sh \
#          analysis/configs/hindcast_horizons.yml \
#          generated_data/outbreaks_survival_data_base-1990_2019.rds
#
# Args: $1 = hindcast manifest (hindcast_horizons.yml)
#       $2 = observed survival panel RDS (--obs-panel)
#
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=10:00:00
# 8 runs scored in parallel by default (07c/07d --workers = the CPUs allocated
# here); each forked worker holds one run + the observed panel, hence 64G.
# Override at submission: sbatch -c N --mem=SIZE (stages/07_hindcast.sh
# --pp-workers N --pp-mem SIZE); -c 1 runs serially.
# 120G / 10 h sizes 10000 replicates: 5000 took 2 h 20 and 29.3 GB, so
# 64G / 4 h is too tight once doubled.
#SBATCH --mem=120G
#SBATCH -c 8

# R 4.3.2 + GDAL stack (as OutbreakExtractR and the 09 results job): the
# 07c/07d figures need sf, absent from the R 4.2.1 stack.
source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env --purge GEO

echo "***** Beginning of script *****"
date

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

MANIFEST="$1"
OBS_PANEL="$2"
if [[ -z "${MANIFEST}" || -z "${OBS_PANEL}" ]]; then
  echo "Error: usage: slurm/submit_hindcast_postprocess.sh <hindcast_horizons.yml> <obs_panel.rds>" >&2
  exit 1
fi
if [[ ! -f "${MANIFEST}"  ]]; then echo "Error: manifest not found: ${MANIFEST}" >&2;  exit 1; fi
if [[ ! -f "${OBS_PANEL}" ]]; then echo "Error: obs panel not found: ${OBS_PANEL}" >&2; exit 1; fi

echo "Manifest  : ${MANIFEST}"
echo "Obs panel : ${OBS_PANEL}"

# One forked worker per allocated CPU (07c/07d --workers; runs are independent).
# 8 with this script's default `-c 8`; stages/07_hindcast.sh --pp-workers N changes it.
WORKERS="${SLURM_CPUS_PER_TASK:-1}"
echo "Workers   : ${WORKERS}"

echo "--- 07c_simulator_postprocess.R --hindcast (aggregate + skill) ..."
Rscript analysis/07c_simulator_postprocess.R \
  --hindcast --manifest "${MANIFEST}" --obs-panel "${OBS_PANEL}" --workers "${WORKERS}" || exit 1

echo "--- 07d_simulator_figures.R --hindcast (skill figures) ..."
Rscript analysis/07d_simulator_figures.R \
  --hindcast --manifest "${MANIFEST}" --obs-panel "${OBS_PANEL}" --workers "${WORKERS}" || exit 1

echo "***** End of script *****"
date
