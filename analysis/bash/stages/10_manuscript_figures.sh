#!/bin/bash
# stages/10_manuscript_figures.sh — submit the manuscript figures (phase 10)
# and the key-numbers report (phase 11) as one SLURM job.
#
# Login-node driver over slurm/submit_manuscript.sh (tier 0), which
# runs, in one job: the package check (10/11 R packages, taxdat, sf's GDAL
# linkage, quarto), Rscript analysis/10_manuscript_figures.R and
# quarto render analysis/11_manuscript_numbers.qmd. Both read only the
# phase-09 artefacts (results_require()), so run this after
# stages/09_postprocess_results.sh: pass its RESULTS_JOB_IDS as --after.
#
# Run-agnostic: the results config names the run (results.deploy_config), so
# no --deploy (as stages/09_postprocess_results.sh).
#
# Usage:
#   bash analysis/bash/stages/10_manuscript_figures.sh --config analysis/configs/results/manuscript.yml
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --after 123:456        # afterok the 09 waves
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --figures 4,H --supplementary false
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --skip-10              # report only
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --skip-11              # figures only
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --check                # package check here, then submit
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --check-only           # package check here, no submit
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --mem 64G --time 04:00:00
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --dry-run              # print the sbatch command only
#   bash analysis/bash/stages/10_manuscript_figures.sh --config <cfg> --allow-stale-inputs   # read artefacts whose inputs changed
#
# --after IDS (colon- or comma-separated job ids): --dependency=afterok on them.
# --check / --check-only run the payload's package check on this node (the
# payload loads its own modules), so a missing package fails before queueing.
#
# --allow-stale-inputs exports RESULTS_ALLOW_STALE_INPUTS=true to the job:
# results_require() warns instead of stopping on artefacts whose inputs changed
# since 09 derived them (stamp, schema and decision checks still stop). For a
# rerun that must not pick up an upstream re-prep; use deliberately.
#
# Prints MANUSCRIPT_JOB_ID=<id> on stdout.

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"

CONFIG=""
FIGURES=""
SUPPLEMENTARY=""
SKIP_10=false
SKIP_11=false
CHECK=false
CHECK_ONLY=false
DRY_RUN=false
MEM=""
CPUS=""
TIME=""
AFTER=""
JOB_NAME=""
ALLOW_STALE=false

usage() {
  echo "Usage: $0 --config <results cfg.yml> [--figures 1,2,3,4,H] [--supplementary true|false] [--skip-10] [--skip-11] [--after ID:ID] [--mem M] [--cpus N] [--time T] [--job-name NAME] [--allow-stale-inputs] [--check] [--check-only] [--dry-run]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)     CONFIG="$2";        shift 2 ;;
    --figures)       FIGURES="$2";       shift 2 ;;
    --supplementary) SUPPLEMENTARY="$2"; shift 2 ;;
    --skip-10)       SKIP_10=true;       shift   ;;
    --skip-11)       SKIP_11=true;       shift   ;;
    --after)         AFTER="${2//,/:}";  shift 2 ;;
    --mem)           MEM="$2";           shift 2 ;;
    --cpus)          CPUS="$2";          shift 2 ;;
    --time)          TIME="$2";          shift 2 ;;
    --job-name)      JOB_NAME="$2";      shift 2 ;;
    --allow-stale-inputs) ALLOW_STALE=true; shift ;;
    --check)         CHECK=true;         shift   ;;
    --check-only)    CHECK_ONLY=true;    shift   ;;
    --dry-run)       DRY_RUN=true;       shift   ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then echo "Error: --config is required." >&2; usage; exit 1; fi
if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG (run from the repo root)" >&2; exit 1; fi
if [[ "$SKIP_10" == "true" && "$SKIP_11" == "true" ]]; then
  echo "Error: --skip-10 and --skip-11 together leave nothing to run." >&2; exit 1
fi
AFTER="${AFTER#:}"; AFTER="${AFTER%:}"
if [[ -n "$AFTER" && ! "$AFTER" =~ ^[^:]+(:[^:]+)*$ ]]; then
  echo "Error: --after expects job ids separated by ':' or ',' (got '$AFTER')." >&2; exit 1
fi

PAYLOAD=analysis/bash/slurm/submit_manuscript.sh
PAYLOAD_ARGS=(--config "$CONFIG")
[[ -n "$FIGURES" ]]         && PAYLOAD_ARGS+=(--figures "$FIGURES")
[[ -n "$SUPPLEMENTARY" ]]   && PAYLOAD_ARGS+=(--supplementary "$SUPPLEMENTARY")
[[ "$SKIP_10" == "true" ]]  && PAYLOAD_ARGS+=(--skip-10)
[[ "$SKIP_11" == "true" ]]  && PAYLOAD_ARGS+=(--skip-11)

# ---------------------------------------------------------------------------
# 1. Optional package check on this node (same check the job runs first).
# ---------------------------------------------------------------------------
if [[ "$CHECK" == "true" || "$CHECK_ONLY" == "true" ]]; then
  echo "--- Package check (${PAYLOAD} --check-only) ..." >&2
  bash "$PAYLOAD" "${PAYLOAD_ARGS[@]}" --check-only >&2
  [[ "$CHECK_ONLY" == "true" ]] && exit 0
fi

# ---------------------------------------------------------------------------
# 2. Submit the 10 + 11 job.
# ---------------------------------------------------------------------------
mkdir -p "$CHOLERA_LOG_DIR"

# Job name: --job-name, else the payload's #SBATCH default (manuscript_10_11).
SBATCH_ARGS=(--parsable)
[[ -n "$JOB_NAME" ]] && SBATCH_ARGS+=(--job-name="$JOB_NAME")
[[ -n "$AFTER" ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER}")
[[ -n "$MEM"  ]]  && SBATCH_ARGS+=(--mem="$MEM")
[[ -n "$CPUS" ]]  && SBATCH_ARGS+=(-c "$CPUS")
[[ -n "$TIME" ]]  && SBATCH_ARGS+=(--time="$TIME")
[[ "$ALLOW_STALE" == "true" ]] && SBATCH_ARGS+=(--export=ALL,RESULTS_ALLOW_STALE_INPUTS=true)

cholera_sbatch_cmd SHORT "${SBATCH_ARGS[@]}" "$PAYLOAD" "${PAYLOAD_ARGS[@]}"
CMD=("${CHOLERA_SBATCH_CMD[@]}")
if [[ "$DRY_RUN" == "true" ]]; then
  echo "[dry-run] ${CMD[*]}" >&2
  JOB_ID="<manuscript>"
else
  echo "--- Submitting 10 + 11 ${AFTER:+(afterok:${AFTER}) }..." >&2
  JOB_ID=$("${CMD[@]}")
  echo "Manuscript job id: ${JOB_ID}" >&2
fi

echo "MANUSCRIPT_JOB_ID=${JOB_ID}"
