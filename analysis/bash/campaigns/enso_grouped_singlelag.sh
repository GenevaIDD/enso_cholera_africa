#!/bin/bash
# Submit the ENSO-grouped single-lag sweep (one config per ENSO group x lag)
# for the occurrence track as a SLURM array job. The batch name is the deploy
# config's fits.occ_dir. Siblings: campaigns/dur_enso_grouped_singlelag.sh and
# campaigns/mon_enso_grouped_singlelag.sh. (The size track is retired downstream
# and is not fit by this campaign.)
#
# Usage:
#   bash analysis/bash/campaigns/enso_grouped_singlelag.sh --deploy DEPLOY.yml [--skip-postprocess] [--concurrency K] [--partition P[,P2...]] [--time T]
#
# Options:
#   --deploy FILE        Deploy config naming the batch (REQUIRED)
#   --skip-postprocess   Pass through to stages/04_fit.sh (skip 05_fit_postprocess.R)
#   --concurrency K      Cap simultaneous array slots (e.g. --concurrency 10)
#   --partition P        Pass through to stages/04_fit.sh [default there: the cluster
#                        profile's LONG partition]. A comma list lets SLURM start each
#                        array task on whichever listed partition frees up first; --time
#                        must fit the shortest listed partition's limit.
#   --time T             Pass through to stages/04_fit.sh [default there: 3-23:00:00]
#
# Example:
#   bash analysis/bash/campaigns/enso_grouped_singlelag.sh --deploy analysis/configs/deploy/prod_run<N>.yml --concurrency 10

set -euo pipefail

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

SKIP_POSTPROCESS=""
CONCURRENCY=""
PARTITION_ARG=""
TIME_ARG=""
DEPLOY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deploy)            DEPLOY="$2";                           shift 2 ;;
    --skip-postprocess)  SKIP_POSTPROCESS="--skip-postprocess"; shift   ;;
    --concurrency|-k)    CONCURRENCY="--concurrency $2";        shift 2 ;;
    --partition|-p)      PARTITION_ARG="--partition $2";        shift 2 ;;
    --time|-t)           TIME_ARG="--time $2";                  shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --deploy DEPLOY.yml [--skip-postprocess] [--concurrency K] [--partition P[,P2...]] [--time T]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$DEPLOY" ]]; then
  echo "Error: --deploy analysis/configs/deploy/prod_run<N>.yml is required (it names the batch)." >&2
  exit 1
fi
# Batch names come from the deploy config (R owns them; KEY=value contract of
# deploy_prod_run.R --stage resolve). Generate the batch dirs first with
#   Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_batches
RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }

SCRIPT="analysis/bash/stages/04_fit.sh"
BATCH=$(deploy_key OCC_BATCH)

echo "=== Submitting ${BATCH} (occurrence) ==="
bash "$SCRIPT" --batch "$BATCH" $SKIP_POSTPROCESS $CONCURRENCY $PARTITION_ARG $TIME_ARG

echo "Done. Array job submitted."
