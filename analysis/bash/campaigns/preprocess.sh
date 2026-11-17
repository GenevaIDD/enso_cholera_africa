#!/bin/bash
# campaigns/preprocess.sh — the 00 -> 02 data leg on a SLURM cluster, one command.
#
# Rebuilds every outbreak-vintage-dependent artefact of phases 00-02 as chained
# Slurm jobs (afterok):
#
#   00 ingest                    stages/00_ingest.sh
#    ├─ 01 outbreaks             stages/01_preprocess.sh
#    │   ├─ 01b monthly
#    │   └─ 01 national zero     (the --deploy's deploy.episodes rule)
#    └─ 02 population            stages/02_population.sh (per-country array -> 02c)
#   02 crosswalk                 stages/02_crosswalk.sh  (after 01 outbreaks + 02c)
#   02 climate                   stages/02_climate.sh    (after 02 crosswalk)
#
# Population reads 00's LP shapefiles (not 01's gpkg), so it runs beside 01.
# The order follows analysis/pipeline_steps.yml.
#
# Run-agnostic: it takes no run. --deploy is read for the national-zero rule
# only (or pass --nz-rule / --nz-threshold).
#
# Before the first run on a vintage (login node, from the repo root):
#   git pull                                     # the outbreak vintage and climate caches copied over
#   bash analysis/bash/campaigns/preprocess.sh --check-packages
#   bash analysis/bash/campaigns/preprocess.sh --config C --deploy D --plan
# After it finishes, copy the outputs back to the local checkout
# (utils/preprocess_paths.R --outputs lists them).
#
# Usage:
#   bash analysis/bash/campaigns/preprocess.sh --config <population config.yml> --deploy <deploy.yml>
#       [--from 00|01|02pop|02xw|02clim] [--until 00|01|02pop|02xw|02clim]
#       [--indir DIR] [--task-unit country|country_year] [--plan] [--check-packages]
#
# --from / --until  run a contiguous slice of 00 01 02pop 02xw 02clim; steps
#                   outside it are assumed done (no dependency on them)
# --indir           outbreak vintage folder [default: paths.R::default_outbreak_data_dir()]
# --plan            every stage prints its sbatch lines; nothing is submitted
# --check-packages  check the R packages of the 00-02 steps (and stages/03_prepare.sh) on this node; exit

set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/../../.."

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
cholera_env GEO

STEPS=(00 01 02pop 02xw 02clim)
CONFIG=""
DEPLOY=""
NZ_RULE=""
NZ_THRESHOLD=""
INDIR=""
TASK_UNIT="country"
FROM="00"
UNTIL="02clim"
PLAN=false
CHECK=false

usage() {
  echo "Usage: $0 --config C.yml (--deploy D.yml | --nz-rule R --nz-threshold N) [--from S] [--until S] [--indir DIR] [--task-unit U] [--plan] [--check-packages]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)      CONFIG="$2";       shift 2 ;;
    --deploy)         DEPLOY="$2";       shift 2 ;;
    --nz-rule)        NZ_RULE="$2";      shift 2 ;;
    --nz-threshold)   NZ_THRESHOLD="$2"; shift 2 ;;
    --indir)          INDIR="$2";        shift 2 ;;
    --task-unit)      TASK_UNIT="$2";    shift 2 ;;
    --from)           FROM="$2";         shift 2 ;;
    --until)          UNTIL="$2";        shift 2 ;;
    --plan)           PLAN=true;         shift   ;;
    --check-packages) CHECK=true;        shift   ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ "$CHECK" == "true" ]]; then
  Rscript --vanilla - <<'RSCRIPT'
pkgs <- c("dplyr", "purrr", "stringr", "lubridate", "tidyr", "sf", "lwgeom", "furrr",
          "future", "optparse", "yaml", "qs2", "glue", "readr", "arrow", "RSQLite", "DBI",
          "exactextractr", "raster", "giscoR", "taxdat", "OutbreakExtractR", "digest",
          # stages/03_prepare.sh (same R 4.3.2 stack): survival and tidyverse are
          # not shipped with it; ncdf4 reads the CFSv2 NetCDF (projection panel).
          "survival", "tidyverse", "ncdf4")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
for (p in pkgs) cat(sprintf("  %-4s %s%s\n", if (ok[[p]]) "OK" else "MISS", p,
                            if (ok[[p]]) paste0(" ", utils::packageVersion(p)) else ""))
cat(sprintf("R %s; .libPaths(): %s\n", getRversion(), paste(.libPaths(), collapse = " ")))
if (!all(ok)) { cat("Missing:", paste(pkgs[!ok], collapse = ", "), "\n"); quit(status = 1) }
RSCRIPT
  exit $?
fi

idx() { local i; for i in "${!STEPS[@]}"; do [[ "${STEPS[$i]}" == "$1" ]] && { echo "$i"; return; }; done; }
I_FROM=$(idx "$FROM" || true); I_UNTIL=$(idx "$UNTIL" || true)
if [[ -z "$I_FROM" || -z "$I_UNTIL" || $I_FROM -gt $I_UNTIL ]]; then
  echo "Error: --from/--until must be steps of ${STEPS[*]}, in order (got ${FROM}..${UNTIL})." >&2; exit 1
fi
want() { local i; i=$(idx "$1"); [[ $i -ge $I_FROM && $i -le $I_UNTIL ]]; }

if [[ -z "$CONFIG" ]] && { want 02pop || want 02xw; }; then
  echo "Error: --config is required for 02 population / crosswalk (its population: block)." >&2; usage; exit 1
fi
if want 01 && [[ -z "$DEPLOY" && -z "$NZ_THRESHOLD" ]]; then
  echo "Error: 01 needs --deploy D.yml or --nz-rule/--nz-threshold." >&2; usage; exit 1
fi

COMMON=()
[[ "$PLAN" == "true" ]] && COMMON+=(--plan)
IND=()
[[ -n "$INDIR" ]] && IND=(--indir "$INDIR")

# run_stage KEY... -- CMD...: run a stage driver, echo its output, export KEY=value
run_stage() {
  local keys=() out k
  while [[ "$1" != "--" ]]; do keys+=("$1"); shift; done; shift
  echo "== $*" >&2
  out=$("$@")
  echo "$out"
  for k in "${keys[@]}"; do
    printf -v "$k" '%s' "$(echo "$out" | grep "^${k}=" | tail -1 | cut -d= -f2- || true)"
  done
}
join_ids() { local IFS=:; local a=(); local x; for x in "$@"; do [[ -n "$x" ]] && a+=("$x"); done; echo "${a[*]}"; }

INGEST_JOB_ID=""; OUTBREAKS_JOB_ID=""; POP_COMPILE_JOB_ID=""; CROSSWALK_JOB_ID=""

if want 00; then
  run_stage INGEST_JOB_ID -- bash analysis/bash/stages/00_ingest.sh "${IND[@]}" "${COMMON[@]}"
fi

if want 01; then
  NZ=()
  [[ -n "$DEPLOY" ]] && NZ+=(--deploy "$DEPLOY")
  [[ -n "$NZ_RULE$NZ_THRESHOLD" ]] && NZ+=(--nz-rule "$NZ_RULE" --nz-threshold "$NZ_THRESHOLD")
  AFT=(); [[ -n "$INGEST_JOB_ID" ]] && AFT=(--after "$INGEST_JOB_ID")
  run_stage OUTBREAKS_JOB_ID MONTHLY_JOB_ID NZ_JOB_ID -- \
    bash analysis/bash/stages/01_preprocess.sh "${NZ[@]}" "${IND[@]}" "${AFT[@]}" "${COMMON[@]}"
fi

if want 02pop; then
  AFT=(); [[ -n "$INGEST_JOB_ID" ]] && AFT=(--after "$INGEST_JOB_ID")
  # --resume judges freshness against the LP corpus on disk NOW; with 00 queued
  # in this chain that corpus is about to change, so submit every task (each
  # task still skips the years whose partition is fresh when it runs).
  RES=(); want 00 || RES=(--resume)
  run_stage POP_ARRAY_JOB_ID POP_COMPILE_JOB_ID -- \
    bash analysis/bash/stages/02_population.sh --config "$CONFIG" --task-unit "$TASK_UNIT" \
      "${RES[@]}" "${AFT[@]}" "${COMMON[@]}"
fi

if want 02xw; then
  DEPS=$(join_ids "$OUTBREAKS_JOB_ID" "$POP_COMPILE_JOB_ID")
  AFT=(); [[ -n "$DEPS" ]] && AFT=(--after "$DEPS")
  run_stage CROSSWALK_JOB_ID -- \
    bash analysis/bash/stages/02_crosswalk.sh --config "$CONFIG" "${IND[@]}" "${AFT[@]}" "${COMMON[@]}"
fi

if want 02clim; then
  DEPS=$(join_ids "$CROSSWALK_JOB_ID" "$OUTBREAKS_JOB_ID")
  AFT=(); [[ -n "$DEPS" ]] && AFT=(--after "$DEPS")
  run_stage CLIMATE_JOB_ID -- \
    bash analysis/bash/stages/02_climate.sh "${IND[@]}" "${AFT[@]}" "${COMMON[@]}"
fi

echo ""
if [[ "$PLAN" == "true" ]]; then
  echo "Plan only: nothing submitted (steps ${FROM}..${UNTIL})."
else
  echo "Submitted steps ${FROM}..${UNTIL}. When done, copy the outputs back (Rscript analysis/utils/preprocess_paths.R --outputs lists them)."
fi
