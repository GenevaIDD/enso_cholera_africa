#!/bin/bash
# stages/05_loo_compare.sh — submit the per-group LOO lag re-selection
# (05_fit_postprocess_loo_lags.R) as a SLURM array: one task per
# (track, group), after the grouped ENSO single-lag sweeps finish.
#
# Each task writes its ranking to generated_data/loo_compare/<batch>/
# (paths.R::make_loo_compare_run_file()), which is gitignored: copy it back to
# the local checkout if the fits ran on a cluster.
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml):
#   bash analysis/bash/stages/05_loo_compare.sh --deploy $D
#   bash analysis/bash/stages/05_loo_compare.sh --deploy $D --track occ --concurrency 4
#   bash analysis/bash/stages/05_loo_compare.sh --track mon --mon-dir <dir>
#   bash analysis/bash/stages/05_loo_compare.sh --deploy $D --no-redo --list
#   bash analysis/bash/stages/05_loo_compare.sh --deploy $D --no-redo \
#       --extra "--lag_sidecar --sidecar_only" --time 1:00:00   # sidecar only
#
# Options:
#   --deploy FILE        Deploy config: batch dirs and group list via
#                        deploy_prod_run.R --stage resolve
#   --track T            occ | dur | mon | all [default: all = occ+dur+mon]
#   --groups G1,G2,...   ENSO groups [default: from --deploy, else enso_east,
#                        enso_horn, enso_south, enso_west]
#   --occ-dir / --dur-dir / --mon-dir DIR
#                        Batch config dirs [default: from --deploy; required
#                        otherwise]. Explicit dirs / --groups win over --deploy.
#   --concurrency K      Cap simultaneous array slots [default: all]
#   --no-redo            Reuse cached per-fit LOO (drop --redo)
#   --lag-sidecar        Pass --lag_sidecar (per-(batch, group) sidecar too)
#   --sidecar-only       Pass --sidecar_only (sidecar alone; rankings untouched,
#                        so 09's enso_relation freshness does not move)
#   --extra "ARGS"       Extra 05_fit_postprocess_loo_lags.R arguments, verbatim
#   --time T             Wall time per task [default: the payload's #SBATCH]
#   --dependency SPEC    Forwarded to sbatch --dependency
#   --list               Print the (track, group) jobs and exit

set -euo pipefail

# Env entrypoint: the cluster profile's R environment (lib/cluster.sh,
# cholera_env; a no-op where the profile loads nothing, e.g. locally).
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env

TRACK="all"
# NB: do NOT name this GROUPS — that is a reserved bash special array (user gids).
GROUP_CSV=""
DEPLOY=""
OCC_DIR=""
DUR_DIR=""
MON_DIR=""
CONCURRENCY=""
REDO="--redo"
EXTRA=""
TIME=""
DEPENDENCY=""
LIST_ONLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --track)        TRACK="$2";        shift 2 ;;
    --deploy)       DEPLOY="$2";       shift 2 ;;
    --groups)       GROUP_CSV="$2";    shift 2 ;;
    --occ-dir)      OCC_DIR="$2";      shift 2 ;;
    --dur-dir)      DUR_DIR="$2";      shift 2 ;;
    --mon-dir)      MON_DIR="$2";      shift 2 ;;
    --concurrency)  CONCURRENCY="$2";  shift 2 ;;
    --no-redo)      REDO="";           shift   ;;
    --lag-sidecar)  EXTRA="${EXTRA:+$EXTRA }--lag_sidecar";  shift ;;
    --sidecar-only) EXTRA="${EXTRA:+$EXTRA }--sidecar_only"; shift ;;
    --extra)        EXTRA="${EXTRA:+$EXTRA }$2"; shift 2 ;;
    --time)         TIME="$2";         shift 2 ;;
    --dependency)   DEPENDENCY="$2";   shift 2 ;;
    --list)         LIST_ONLY=true;    shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--deploy DEPLOY.yml] [--track {occ|dur|mon|all}] [--groups G1,G2,...] [--occ-dir DIR] [--dur-dir DIR] [--mon-dir DIR] [--concurrency K] [--no-redo] [--lag-sidecar] [--sidecar-only] [--extra \"ARGS\"] [--time T] [--dependency SPEC] [--list]" >&2
      exit 1 ;;
  esac
done

case "$TRACK" in
  occ|dur|mon|all) ;;
  *) echo "Error: --track must be one of {occ|dur|mon|all}, got '$TRACK'." >&2; exit 1 ;;
esac

# Batch dirs come from the deploy config (R owns them; KEY=value contract of
# deploy_prod_run.R --stage resolve). Explicit flags win.
if [[ -n "$DEPLOY" ]]; then
  RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve) || exit 1
  deploy_key() { echo "$RESOLVED" | grep "^$1=" | cut -d= -f2-; }
  [[ -n "$OCC_DIR"  ]] || OCC_DIR=$(deploy_key OCC_DIR)
  [[ -n "$DUR_DIR"  ]] || DUR_DIR=$(deploy_key DUR_DIR)
  [[ -n "$MON_DIR"  ]] || MON_DIR=$(deploy_key MON_DIR)
  [[ -n "$GROUP_CSV" ]] || GROUP_CSV=$(deploy_key GROUPS | tr ' ' ',')
fi
GROUP_CSV="${GROUP_CSV:-enso_east,enso_horn,enso_south,enso_west}"

# ---------------------------------------------------------------------------
# 1. Build the (batch_dir, group) job list as one entry per array task.
#    Entry format: "<batch_dir>::<group>" (":: " is the field separator the
#    inner submit script splits on).
# ---------------------------------------------------------------------------
declare -a TRACK_DIRS=()
if [[ "$TRACK" == "occ" || "$TRACK" == "all" ]]; then
  [[ -z "$OCC_DIR" ]] && { echo "Error: --occ-dir or --deploy required for --track $TRACK." >&2; exit 1; }
  TRACK_DIRS+=("$OCC_DIR")
fi
if [[ "$TRACK" == "dur" || "$TRACK" == "all" ]]; then
  [[ -z "$DUR_DIR" ]] && { echo "Error: --dur-dir or --deploy required for --track $TRACK." >&2; exit 1; }
  TRACK_DIRS+=("$DUR_DIR")
fi
if [[ "$TRACK" == "mon" || "$TRACK" == "all" ]]; then
  [[ -z "$MON_DIR" ]] && { echo "Error: --mon-dir or --deploy required for --track $TRACK." >&2; exit 1; }
  TRACK_DIRS+=("$MON_DIR")
fi

IFS=',' read -r -a GROUP_ARR <<< "$GROUP_CSV"

declare -a JOBS=()
for d in "${TRACK_DIRS[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: batch config dir not found: $d" >&2
    exit 1
  fi
  for g in "${GROUP_ARR[@]}"; do
    g_trimmed="${g//[[:space:]]/}"
    [[ -z "$g_trimmed" ]] && continue
    JOBS+=("${d}::${g_trimmed}")
  done
done

N_JOBS=${#JOBS[@]}
if [[ "$N_JOBS" -eq 0 ]]; then
  echo "Error: empty job list (check --track / --groups)." >&2
  exit 1
fi

echo "=== LOO lag re-selection: ${N_JOBS} array task(s) ==="
[[ -n "$EXTRA" ]] && echo "    05 extra args: ${EXTRA}"
[[ -n "$TIME"  ]] && echo "    wall time    : ${TIME}"
for i in "${!JOBS[@]}"; do
  printf '  [%d] %s\n' "$i" "${JOBS[$i]}"
done

if [[ "$LIST_ONLY" == "true" ]]; then
  echo "(--list: not submitting)"
  exit 0
fi

# Join the job list with ';' into one arg the inner script indexes by array id.
JOBS_STR=$(IFS=';'; echo "${JOBS[*]}")

# ---------------------------------------------------------------------------
# 2. Submit the array.
# ---------------------------------------------------------------------------
ARRAY_SPEC="0-$((N_JOBS - 1))"
[[ -n "$CONCURRENCY" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${CONCURRENCY}"

mkdir -p "$CHOLERA_LOG_DIR"

declare -a TIME_FLAG=()
[[ -n "$TIME" ]] && TIME_FLAG=(--time="$TIME")
declare -a DEP_ARGS=()
[[ -n "$DEPENDENCY" ]] && DEP_ARGS=(--dependency="$DEPENDENCY")

echo "--- Submitting array ${ARRAY_SPEC} as job 'loo_compare_lags' ${REDO:+(--redo)} ${EXTRA:+($EXTRA)} ${DEPENDENCY:+(dependency: $DEPENDENCY)} ..."
ARRAY_JOB_ID=$(cholera_sbatch SHORT --parsable \
  --array="${ARRAY_SPEC}" \
  --job-name="loo_compare_lags" \
  ${TIME_FLAG[@]+"${TIME_FLAG[@]}"} \
  ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} \
  analysis/bash/slurm/submit_loo_compare.sh \
  "$JOBS_STR" "$REDO" "$EXTRA")

echo "Array job id: ${ARRAY_JOB_ID}"
echo "Done."
