#!/bin/bash
# stages/04_fit_kl.sh — fit a brms_kl config batch (e.g. the class sweep) as a
# SLURM array with the kl resources, and optionally compare the sweep once
# every fit has finished.
#
#   first : slurm/submit_kl_compile.sh  (utils/kl_compile_batch.R: every program
#           of the batch compiled once, before any fit; --no-precompile skips it)
#   array : slurm/submit_kl_array.sh  (one 04_fit_run.R per config, -c 32, 24G,
#           afterok the compile job: the tasks reuse the executables, none compiles)
#   then  : slurm/submit_kl_sweep.sh  (05_fit_postprocess_kl_sweep.R, afterok)  [--postprocess]
#
# stages/04_fit.sh cannot be used for brms_kl batches: its payload requests
# -c 16 / 20G and runs 05_fit_postprocess.R, which has nothing to do for a
# kl_fit. Existing fits are skipped by the worker (overwrite: no), so a
# resubmission refits only what is missing.
#
# Usage:
#   bash analysis/bash/stages/04_fit_kl.sh --batch kl_sweep_classes_prodtu_all --postprocess
#   bash analysis/bash/stages/04_fit_kl.sh --deploy <deploy.yml> --postprocess
#   bash analysis/bash/stages/04_fit_kl.sh --batch <name|dir> [--concurrency K] [--partition P]
#        [--time T] [--cores N] [--mem M] [--postprocess] [--mode_tol X] [--n_agree N] [--dry-run] [--list]
#        [--no-precompile]
#
# --deploy reads the batch from the deploy YAML's kl_sweep.batch (KL_SWEEP_BATCH
# of deploy_prod_run.R --stage resolve), and the comparison's chain selection
# from kl_sweep.n_agree / kl_sweep.mode_tol; an explicit --batch / --n_agree /
# --mode_tol wins. This is the
# `classes` phase of bash/tools/run_pipeline.sh, which chains on the
# KL_FIT_JOB_ID= / KL_SWEEP_JOB_ID= lines printed on submission.
#
# The batch is generated with
#   Rscript analysis/utils/make_batch_configs.R --batch <name>
# or, for a deploy's sweep, Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_batches

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"

BATCH=""
DEPLOY=""
CONCURRENCY=""
PARTITION=""              # "" = the cluster profile's SHORT partition(s)
TIME_LIMIT="03:00:00"     # a short-sampler fit takes 35-40 min (200 warm-up / 200 sampling)
CORES=""
MEM=""
POSTPROCESS=FALSE
MODE_TOL=""
N_AGREE=""
DRY_RUN=FALSE
LIST=FALSE
PRECOMPILE=TRUE

usage() {
  echo "Usage: $0 --batch <name|dir> | --deploy <deploy.yml> [--concurrency K] [--partition P] [--time T] [--cores N] [--mem M]" >&2
  echo "          [--postprocess] [--mode_tol X] [--n_agree N] [--dry-run] [--list] [--no-precompile]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch|-b)       BATCH="$2";       shift 2 ;;
    --deploy)         DEPLOY="$2";      shift 2 ;;
    --concurrency|-k) CONCURRENCY="$2"; shift 2 ;;
    --partition|-p)   PARTITION="$2";   shift 2 ;;
    --time|-t)        TIME_LIMIT="$2";  shift 2 ;;
    --cores|-c)       CORES="$2";       shift 2 ;;
    --mem)            MEM="$2";         shift 2 ;;
    --postprocess)    POSTPROCESS=TRUE; shift   ;;
    --mode_tol)       MODE_TOL="$2";    shift 2 ;;
    --n_agree)        N_AGREE="$2";     shift 2 ;;
    --dry-run)        DRY_RUN=TRUE;     shift   ;;
    --list)           LIST=TRUE;        shift   ;;
    --no-precompile)  PRECOMPILE=FALSE; shift   ;;
    -h|--help)        usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  if [[ -z "$BATCH" ]]; then
    BATCH=$(deploy_key KL_SWEEP_BATCH)
    if [[ -z "$BATCH" ]]; then
      echo "Error: $DEPLOY sets no deploy.kl_sweep.batch (no class sweep in this rollout)." >&2; exit 1
    fi
  fi
  [[ -n "$N_AGREE" ]]  || N_AGREE=$(deploy_key KL_SWEEP_N_AGREE)
  [[ -n "$MODE_TOL" ]] || MODE_TOL=$(deploy_key KL_SWEEP_MODE_TOL)
fi
if [[ -z "$BATCH" ]]; then
  echo "Error: --batch or --deploy is required." >&2; usage; exit 1
fi
if [[ -d "$BATCH" ]]; then
  CONFIGDIR="${BATCH%/}"
else
  CONFIGDIR="analysis/configs/${BATCH}"
fi
if [[ ! -d "$CONFIGDIR" ]]; then
  echo "Error: config directory '$CONFIGDIR' not found" >&2; exit 1
fi
NAME="$(basename "$CONFIGDIR")"

CONFIGNAMES=($(ls "$CONFIGDIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | sort))
N=${#CONFIGNAMES[@]}
if [[ "$N" -eq 0 ]]; then
  echo "Error: no config files found in '$CONFIGDIR'" >&2; exit 1
fi

if [[ "$LIST" == TRUE ]]; then
  for i in "${!CONFIGNAMES[@]}"; do printf "%3d  %s\n" "$i" "${CONFIGNAMES[$i]}"; done
  exit 0
fi

ARRAY_SPEC="0-$((N - 1))"
[[ -n "$CONCURRENCY" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${CONCURRENCY}"

FIT_ARGS=(--parsable --array="${ARRAY_SPEC}" --job-name="${NAME}" --time="${TIME_LIMIT}")
[[ -n "$PARTITION" ]] && FIT_ARGS+=(--partition="${PARTITION}")
[[ -n "$CORES" ]] && FIT_ARGS+=(-c "${CORES}")
[[ -n "$MEM" ]]   && FIT_ARGS+=(--mem="${MEM}")
cholera_sbatch_cmd SHORT --parsable --job-name="${NAME}_compile" ${PARTITION:+--partition="${PARTITION}"} \
  analysis/bash/slurm/submit_kl_compile.sh "${CONFIGDIR}"
COMPILE_CMD=("${CHOLERA_SBATCH_CMD[@]}")

PP_OPTS=()
[[ -n "$MODE_TOL" ]] && PP_OPTS+=(--mode_tol "${MODE_TOL}")
[[ -n "$N_AGREE" ]]  && PP_OPTS+=(--n_agree "${N_AGREE}")

echo "Batch '${NAME}': ${N} config(s), array ${ARRAY_SPEC}, partition=${PARTITION:-<profile SHORT>}, time=${TIME_LIMIT}, postprocess=${POSTPROCESS}, precompile=${PRECOMPILE}"
if [[ "$DRY_RUN" == TRUE ]]; then
  if [[ "$PRECOMPILE" == TRUE ]]; then
    echo "[dry-run] ${COMPILE_CMD[*]}"
    cholera_sbatch_cmd SHORT --dependency=afterok:"<compile>" "${FIT_ARGS[@]}" analysis/bash/slurm/submit_kl_array.sh "${CONFIGDIR}"
  else
    cholera_sbatch_cmd SHORT "${FIT_ARGS[@]}" analysis/bash/slurm/submit_kl_array.sh "${CONFIGDIR}"
  fi
  echo "[dry-run] ${CHOLERA_SBATCH_CMD[*]}"
  if [[ "$POSTPROCESS" == TRUE ]]; then
    cholera_sbatch_cmd LONG --parsable --dependency=afterok:"<array>" analysis/bash/slurm/submit_kl_sweep.sh "${CONFIGDIR}" ${PP_OPTS[@]+"${PP_OPTS[@]}"}
    echo "[dry-run] ${CHOLERA_SBATCH_CMD[*]}"
  fi
  exit 0
fi

if [[ "$PRECOMPILE" == TRUE ]]; then
  CJOB=$("${COMPILE_CMD[@]}")
  echo "Compile job: ${CJOB}"
  echo "KL_COMPILE_JOB_ID=${CJOB}"
  FIT_ARGS+=(--dependency=afterok:"${CJOB}")
fi
JOB=$(cholera_sbatch SHORT "${FIT_ARGS[@]}" analysis/bash/slurm/submit_kl_array.sh "${CONFIGDIR}")
echo "Fit array job: ${JOB}"
echo "KL_FIT_JOB_ID=${JOB}"
if [[ "$POSTPROCESS" == TRUE ]]; then
  PP=$(cholera_sbatch LONG --parsable --dependency=afterok:"${JOB}" \
         analysis/bash/slurm/submit_kl_sweep.sh "${CONFIGDIR}" ${PP_OPTS[@]+"${PP_OPTS[@]}"})
  echo "Sweep comparison job: ${PP} (afterok ${JOB})"
  echo "KL_SWEEP_JOB_ID=${PP}"
fi
