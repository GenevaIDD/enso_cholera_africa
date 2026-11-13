#!/bin/bash
# campaigns/postfit.sh — full postprocessing + figures
# pipeline for grouped-ENSO single-lag batches of any track (brms_cb configs
# sl-1..11 x grp-{east,horn,south,west}; the deploy config's fits.<track>_dir),
# chained as dependent SLURM array jobs. Run this after the
# 04_fit_run.R fitting arrays (campaigns/enso_grouped_singlelag.sh,
# campaigns/dur_enso_grouped_singlelag.sh, campaigns/mon_enso_grouped_singlelag.sh
# with --deploy <deploy.yml>) have completed for every config of the tracks
# requested. The default --tracks occ runs the occurrence track only.
#
# Chain, per track (each step an array job; steps 1 and 2 run in parallel,
# step 3 waits on both via --dependency=afterok). Tracks are independent of
# each other and are submitted back to back:
#   1. CB postprocess   (05_fit_postprocess_cb.R)     — one task per config (44)
#   2. LOO lag compare  (05_fit_postprocess_loo_lags.R) — one task per group (4)
#   3. Lag-compare figs (06_fit_figures_lag_compare.R)  — one task per group (4)
#      needs BOTH step 1 (cb_surfaces) and step 2 (the loo_compare_lags ranking
#      under generated_data/loo_compare/<batch>/)
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml):
#   bash analysis/bash/campaigns/postfit.sh --deploy $D
#   bash analysis/bash/campaigns/postfit.sh --deploy $D --tracks occ,dur,mon
#   bash analysis/bash/campaigns/postfit.sh --occ-dir <dir>
#   bash analysis/bash/campaigns/postfit.sh --tracks dur --dur-dir <dir>
#   bash analysis/bash/campaigns/postfit.sh --deploy $D --no-redo-loo   # reuse cached per-fit loo
#   bash analysis/bash/campaigns/postfit.sh --deploy $D --list          # print job lists, submit nothing
#
# Options:
#   --deploy FILE    Deploy config; supplies the batch dirs + groups via
#                    deploy_prod_run.R --stage resolve
#   --tracks LIST    Comma list of tracks, from {occ,dur,mon} [default: occ]
#                    (size is retired downstream)
#   --occ-dir DIR    Occurrence batch config dir [default: from --deploy]
#   --dur-dir DIR    Duration batch config dir   [default: from --deploy]
#   --mon-dir DIR    Monthly batch config dir    [default: from --deploy]
#                    Each requested track needs its dir, from --deploy or the flag.
#   --groups G,...   ENSO groups      [default: from --deploy, else enso_east,enso_horn,enso_south,enso_west]
#   Explicit --*-dir / --groups win over --deploy.
#   --no-redo-loo    Reuse cached per-fit LOO in step 2 (drop --redo)
#   --format F       Figure format forwarded to step 3 [default: png]
#   --concurrency K  Cap simultaneous array slots in all 3 steps [default: all]
#   --list           Print the job lists for all 3 steps and exit (no submit)

set -euo pipefail

TRACKS="occ"
OCC_DIR=""
DUR_DIR=""
MON_DIR=""
DEPLOY=""
GROUP_CSV=""
REDO_LOO="--redo"
FORMAT="png"
LIST_ONLY=false
CONCURRENCY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tracks)       TRACKS="$2";   shift 2 ;;
    --occ-dir)      OCC_DIR="$2";  shift 2 ;;
    --dur-dir)      DUR_DIR="$2";  shift 2 ;;
    --mon-dir)      MON_DIR="$2";  shift 2 ;;
    --deploy)       DEPLOY="$2";   shift 2 ;;
    --groups)       GROUP_CSV="$2"; shift 2 ;;
    --no-redo-loo)  REDO_LOO="";   shift   ;;
    --format)       FORMAT="$2";   shift 2 ;;
    --list)         LIST_ONLY=true; shift  ;;
    --concurrency)  CONCURRENCY="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--deploy DEPLOY.yml] [--tracks occ,dur,mon] [--occ-dir DIR] [--dur-dir DIR] [--mon-dir DIR] [--groups G1,G2,...] [--no-redo-loo] [--format png|pdf] [--concurrency K] [--list]" >&2
      exit 1 ;;
  esac
done

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env
if [[ -n "$DEPLOY" ]] && ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found on PATH (tried the cluster profile's cholera_env above, if" >&2
  echo "it loads anything). --deploy resolves batch dirs via Rscript. Load" >&2
  echo "an R module (or add Rscript to PATH) and re-run." >&2
  exit 1
fi

# Batch dirs + groups come from the deploy config (R owns them; KEY=value
# contract of deploy_prod_run.R --stage resolve). Explicit flags win.
if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  [[ -n "$OCC_DIR"   ]] || OCC_DIR=$(deploy_key OCC_DIR)
  [[ -n "$DUR_DIR"   ]] || DUR_DIR=$(deploy_key DUR_DIR)
  [[ -n "$MON_DIR"   ]] || MON_DIR=$(deploy_key MON_DIR)
  [[ -n "$GROUP_CSV" ]] || GROUP_CSV=$(deploy_key GROUPS | tr ' ' ',')
fi
GROUP_CSV="${GROUP_CSV:-enso_east,enso_horn,enso_south,enso_west}"

# Validate every requested track up front, so a bad track/dir submits nothing.
IFS=',' read -r -a TRACK_ARR <<< "$TRACKS"
declare -a TRACK_LIST=()
declare -A TRACK_DIR=()
for t in "${TRACK_ARR[@]}"; do
  t="${t//[[:space:]]/}"
  [[ -z "$t" ]] && continue
  case "$t" in
    occ)  d="$OCC_DIR"  ;;
    dur)  d="$DUR_DIR"  ;;
    mon)  d="$MON_DIR"  ;;
    *) echo "Error: --tracks entries must be occ, dur or mon, got '$t'." >&2; exit 1 ;;
  esac
  if [[ -z "$d" ]]; then
    echo "Error: no batch for track '$t'. Pass --deploy analysis/configs/deploy/prod_run<N>.yml or --${t}-dir DIR." >&2
    exit 1
  fi
  if [[ ! -d "$d" ]]; then
    echo "Error: batch dir not found for track '$t': $d" >&2
    exit 1
  fi
  TRACK_LIST+=("$t")
  TRACK_DIR[$t]="$d"
done
if [[ ${#TRACK_LIST[@]} -eq 0 ]]; then
  echo "Error: empty --tracks." >&2
  exit 1
fi

LIST_FLAG=()
[[ "$LIST_ONLY" == "true" ]] && LIST_FLAG=(--list)
CONC_FLAG=()
[[ -n "$CONCURRENCY" ]] && CONC_FLAG=(--concurrency "$CONCURRENCY")

declare -a SUMMARY=()
declare -a BATCHES=()

for TRACK in "${TRACK_LIST[@]}"; do
  BATCH_DIR="${TRACK_DIR[$TRACK]}"
  BATCHES+=("$(basename "$BATCH_DIR")")

  echo ""
  echo "############################################################"
  echo "Track: ${TRACK}  (batch: ${BATCH_DIR})"
  echo "############################################################"

  echo "============================================================"
  echo "[${TRACK}] Step 1/3: CB postprocess (05_fit_postprocess_cb.R), per config"
  echo "============================================================"
  CB_OUTPUT=$(bash analysis/bash/stages/05_postprocess_cb.sh --batch-dir "$BATCH_DIR" "${CONC_FLAG[@]}" "${LIST_FLAG[@]}")
  echo "$CB_OUTPUT"
  CB_JOB_ID=$(echo "$CB_OUTPUT" | grep '^CB_JOB_ID=' | cut -d= -f2- || true)

  echo ""
  echo "============================================================"
  echo "[${TRACK}] Step 2/3: LOO lag comparison (05_fit_postprocess_loo_lags.R), per group"
  echo "============================================================"
  LOO_OUTPUT=$(bash analysis/bash/stages/05_loo_compare.sh \
    --track "$TRACK" "--${TRACK}-dir" "$BATCH_DIR" --groups "$GROUP_CSV" \
    $( [[ -z "$REDO_LOO" ]] && echo "--no-redo" ) \
    "${CONC_FLAG[@]}" "${LIST_FLAG[@]}")
  echo "$LOO_OUTPUT"
  LOO_JOB_ID=$(echo "$LOO_OUTPUT" | grep '^Array job id:' | awk '{print $NF}' || true)

  if [[ "$LIST_ONLY" == "true" ]]; then
    echo ""
    echo "============================================================"
    echo "[${TRACK}] Step 3/3: Lag-compare figures (06_fit_figures_lag_compare.R), per group"
    echo "============================================================"
    bash analysis/bash/stages/06_fit_figures_lag_compare.sh \
      --track "$TRACK" "--${TRACK}-dir" "$BATCH_DIR" --groups "$GROUP_CSV" --format "$FORMAT" "${CONC_FLAG[@]}" --list
    continue
  fi

  if [[ -z "$CB_JOB_ID" || -z "$LOO_JOB_ID" ]]; then
    echo "Error: [${TRACK}] failed to capture one of CB_JOB_ID='${CB_JOB_ID}' / LOO_JOB_ID='${LOO_JOB_ID}'." >&2
    exit 1
  fi

  echo ""
  echo "============================================================"
  echo "[${TRACK}] Step 3/3: Lag-compare figures (06_fit_figures_lag_compare.R), per group"
  echo "  (dependency: afterok:${CB_JOB_ID}:${LOO_JOB_ID})"
  echo "============================================================"
  FIG_OUTPUT=$(bash analysis/bash/stages/06_fit_figures_lag_compare.sh \
    --track "$TRACK" "--${TRACK}-dir" "$BATCH_DIR" --groups "$GROUP_CSV" --format "$FORMAT" \
    "${CONC_FLAG[@]}" --dependency "afterok:${CB_JOB_ID}:${LOO_JOB_ID}")
  echo "$FIG_OUTPUT"
  FIG_JOB_ID=$(echo "$FIG_OUTPUT" | grep '^FIG_JOB_ID=' | cut -d= -f2- || true)

  SUMMARY+=("${TRACK}: CB=${CB_JOB_ID}  LOO=${LOO_JOB_ID}  FIGURES(dep)=${FIG_JOB_ID}")
done

if [[ "$LIST_ONLY" == "true" ]]; then
  echo ""
  echo "(--list: nothing submitted)"
  exit 0
fi


echo ""
echo "============================================================"
echo "Submitted:"
for s in "${SUMMARY[@]}"; do echo "  ${s}"; done
echo "Watch:  squeue -u \$USER"
echo "Outputs: cb_surfaces via make*RunFile (generated_data/), LOO rankings under generated_data/loo_compare/<batch>/, figures dir via makeFigureDir"
echo "============================================================"
