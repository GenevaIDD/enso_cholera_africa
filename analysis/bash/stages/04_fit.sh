#!/bin/bash
# Submit a config batch as a SLURM array job.
# Each array slot fits one config (04_fit_run.R) then postprocesses it (05_fit_postprocess.R).
#
# Usage:
#   bash analysis/bash/stages/04_fit.sh --batch <batch_name>
#   bash analysis/bash/stages/04_fit.sh --batch <batch_name> --skip-postprocess
#   bash analysis/bash/stages/04_fit.sh --batch <batch_name> --concurrency 5
#   bash analysis/bash/stages/04_fit.sh --batch <batch_name> --partition <partition> --time 3-23:00:00
#
# Example:
#   bash analysis/bash/stages/04_fit.sh --batch prod_run_occ_grouped_enso_singlelag

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"

BATCH=""
SKIP_POSTPROCESS=FALSE
CONCURRENCY=""   # no limit by default; set e.g. --concurrency 10 to cap parallel slots
PARTITION=""                # "" = the cluster profile's LONG partition (lib/cluster.sh)
TIME_LIMIT="3-23:00:00"     # default matches slurm/submit_fit.sh's #SBATCH --time

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch|-b)           BATCH="$2";             shift 2 ;;
    --skip-postprocess)   SKIP_POSTPROCESS=TRUE;  shift   ;;
    --concurrency|-k)     CONCURRENCY="$2";       shift 2 ;;
    --partition|-p)       PARTITION="$2";         shift 2 ;;
    --time|-t)            TIME_LIMIT="$2";        shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --batch <batch_name> [--skip-postprocess] [--concurrency K] [--partition P] [--time T]" >&2
      exit 1 ;;
  esac
done

if [[ -z "$BATCH" ]]; then
  echo "Error: --batch is required." >&2
  echo "Usage: $0 --batch <batch_name> [--skip-postprocess] [--concurrency K] [--partition P] [--time T]" >&2
  exit 1
fi

CONFIGDIR="analysis/configs/${BATCH}"

if [[ ! -d "$CONFIGDIR" ]]; then
  echo "Error: config directory '$CONFIGDIR' not found" >&2
  exit 1
fi

N=$(ls "$CONFIGDIR" | grep '\.yml$' | grep -v '^manifest\.yml$' | wc -l)

if [[ "$N" -eq 0 ]]; then
  echo "Error: no config files found in '$CONFIGDIR'" >&2
  exit 1
fi

# Build --array spec: optionally cap concurrency with %K.
ARRAY_SPEC="0-$((N - 1))"
if [[ -n "$CONCURRENCY" ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${CONCURRENCY}"
fi

echo "Submitting batch '${BATCH}': ${N} config(s), array ${ARRAY_SPEC}, skip_postprocess=${SKIP_POSTPROCESS}, partition=${PARTITION:-<profile LONG>}, time=${TIME_LIMIT}"
cholera_sbatch LONG \
  --array="${ARRAY_SPEC}" \
  --job-name="${BATCH}" \
  ${PARTITION:+--partition="${PARTITION}"} \
  --time="${TIME_LIMIT}" \
  analysis/bash/slurm/submit_fit.sh \
  "${BATCH}" "${SKIP_POSTPROCESS}" "${PARTITION}" "${TIME_LIMIT}"
