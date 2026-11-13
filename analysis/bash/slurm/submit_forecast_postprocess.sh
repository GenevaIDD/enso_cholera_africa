#!/bin/bash
# slurm/submit_forecast_postprocess.sh — SBATCH script for the per-config
# forecast/counterfactual postprocess + figures step: once ONE config's forward
# array job has finished, merge its partitions and draw its per-run figures.
#
# Runs, for the single config passed in:
#   1. 07c_simulator_postprocess.R -C <cfg> --aggregate --engine spatpomp
#      — merge all successful partitions (all countries x replicate chunks) into
#        the run's merged/ dir.
#   2. 07d_simulator_figures.R -C <cfg> --engine spatpomp --if-validated
#      — per-run simulator figures for this (arm, group). They are validation
#        figures (PIT, CRPS, ...) read from validation/metrics; a forward run
#        is never validated, so --if-validated makes 07d skip instead of fail.
#
# This is the PER-CONFIG postprocess only. The CROSS-GROUP combine
# (07c --combine-groups) and the cross-group comparison figures are run by
# hand, since they reduce across every group's merged output at once.
#
# Submitted by analysis/bash/stages/07_forecast_groups.sh --postprocess
# (with --dependency=afterok:<that config's array id>). May also be run directly
# once a config's array job has finished:
#   sbatch analysis/bash/slurm/submit_forecast_postprocess.sh \
#          analysis/configs/forecast_country_by_group_<run>/forecast_<unit>.yml
#
# Arg: $1 = config path.
#
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH -c 1

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

CONFIG="$1"
if [[ -z "${CONFIG}" ]]; then
  echo "Error: usage: slurm/submit_forecast_postprocess.sh <config.yml>" >&2
  exit 1
fi
if [[ ! -f "${CONFIG}" ]]; then echo "Error: config not found: ${CONFIG}" >&2; exit 1; fi

echo "Config : ${CONFIG}"

echo "--- 07c_simulator_postprocess.R --aggregate (engine spatpomp) ..."
Rscript analysis/07c_simulator_postprocess.R -C "${CONFIG}" --aggregate --engine spatpomp || exit 1

echo "--- 07d_simulator_figures.R (engine spatpomp) ..."
Rscript analysis/07d_simulator_figures.R -C "${CONFIG}" --engine spatpomp --if-validated || exit 1

echo "***** End of script *****"
date
