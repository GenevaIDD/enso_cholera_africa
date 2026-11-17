#!/bin/bash
# tools/run_pipeline.sh — run a deploy (configs/deploy/<run>.yml) end to end on
# a SLURM cluster (cluster profile, lib/cluster.sh) or this machine (--local).
#
# A relay: each PHASE submits its jobs with the existing drivers, then submits
# the NEXT phase as a small SLURM job (this script, `--phase <next>`) with
# --dependency=afterok on the jobs it needs; later task lists depend on files
# earlier jobs write. Needs `sbatch` on the compute nodes. --local (or
# CHOLERA_BACKEND=local): no relay, phases run inline through lib/local/sbatch
# ($CHOLERA_CORES tasks at a time). It only chains existing drivers in the
# order of analysis/PIPELINE.md (tiers: analysis/bash/README.md).
#
# Phases (-> what the next phase waits on):
#   classes   kl class sweep (only when deploy sets kl_sweep.batch):
#             gen_batches, stages/04_fit_kl.sh --postprocess -> the sweep
#             comparison; stops here until deploy.kl_sweep.cell is pinned
#   fits      [gen_grouping], gen_batches if a fit dir lacks manifest.yml,
#             campaigns/{enso,dur_enso,mon_enso}_grouped_singlelag.sh
#             -> all fit arrays
#   postfit   campaigns/postfit.sh (05 cb surfaces + LOO, 06 figures)
#             -> the cb + LOO arrays (not the figures job)
#   deploy    deploy_fits, gen_templates, smoothing check,
#             stages/08_smoothing.sh -> the smoothing compile job
#   sim       gen_forecast / gen_attrib / gen_hindcast / gen_events, event
#             preflight (refuses a missing Lenssen cache),
#             07_forecast_groups.sh --arm all, 07_hindcast.sh, and
#             07_event_attribution.sh -> forecast arrays, hindcast
#             postprocess and event arrays (not per-run 07c/07d jobs)
#   attrib    only with --from/--phase attrib: gen_attrib +
#             07_forecast_groups.sh --arm attrib (needs
#             deploy.outputs.forecast_attrib_manifest) -> its arrays
#   results   align_results (points results/manuscript.yml at this deploy),
#             stages/09_postprocess_results.sh -> every 09 wave
#   manuscript  stages/10_manuscript_figures.sh (10 + quarto 11) -> that job
#
# Usage (D=analysis/configs/deploy/<run>.yml), from the repo root:
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D            # whole chain
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --plan     # submit nothing
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from fits   # restart
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --until sim   # stop early
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --local --cores 8
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --after 123:456
#
# Options:
#   --deploy FILE         Deploy config (REQUIRED); every path comes from it.
#   --from PHASE          First phase: classes | fits | postfit | deploy | sim |
#                         attrib | results | manuscript [default: classes when
#                         the deploy sets kl_sweep.batch, else fits]
#   --until PHASE         Last phase [default: manuscript]
#   --phase PHASE         Run exactly this phase, then relay (relay jobs)
#   --results-config F    Results config of the results / manuscript phases
#                         [default: analysis/configs/results/manuscript.yml]
#   --concurrency K       Cap simultaneous array slots in every phase
#                         [default: each driver's own]
#   --partition P         Fits and kl sweep: SLURM partition(s)
#   --time T              Fits and kl sweep: time limit
#   --dep TYPE            Between phases: afterok | afterany [default: afterok];
#                         with afterok one failed task stops the chain
#   --relay-time T        Relay job time [default: $CHOLERA_RELAY_TIME, else
#                         03:00:00]; memory / cores: $CHOLERA_RES_RELAY
#   --event-replicates N  Sim: --n_replicates for past-event runs [default: each
#                         config's simulator.n_replicates; prefer the deploy key]
#   --results-mem M       09 task memory [default: $CHOLERA_RESULTS_MEM, else 160G]
#   --results-partition P 09 partition(s) [default: the profile's
#                         CHOLERA_RES_RESULTS, else its SHORT class]
#   --after IDS           Queue the first phase as a relay on these job ids
#   --force-gen           Pass --force to gen_templates and gen_attrib
#   --plan                Print what each phase would run and exit
#   --local               Run on this machine, every phase inline, no relays
#   --cores N             --local: parallel tasks [default: $CHOLERA_CORES, nproc]
#   Hindcast postprocess: --pp-workers $CHOLERA_HINDCAST_PP_WORKERS (default 8),
#   --pp-mem $CHOLERA_HINDCAST_PP_MEM (default 240G). Cluster profile:
#   ${CHOLERA_CLUSTER_PROFILE:-analysis/bash/cluster.profile}.
#
# Restart with --from <phase> after a failure: drivers skip existing outputs
# (smoothing --resume; forecast / event runs with a _SUCCESS sentinel are not
# resubmitted; this script never passes --overwrite). Logs:
# logs/pipe_<phase>_<deploy tag>_<jobid>.log.

set -euo pipefail
# --local / --cores first: lib/cluster.sh reads CHOLERA_BACKEND / CHOLERA_CORES.
for ((_i = 1; _i <= $#; _i++)); do
  case "${!_i}" in
    --local) export CHOLERA_BACKEND=local ;;
    --cores) _j=$((_i + 1)); export CHOLERA_CORES="${!_j}" ;;
  esac
done
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"
# The project root: CHOLERA_PROJECT_ROOT (a sandbox whose analysis/ is a
# symlink, tools/run_synthetic_pipeline.sh), else this checkout.
cd "${CHOLERA_PROJECT_ROOT:-$(dirname "$(readlink -f "$0")")/../../..}"
cholera_env

PHASES=(classes fits postfit deploy sim attrib results manuscript)
DEPLOY=""; FROM=""; UNTIL="manuscript"; ONLY=""; PLAN=false
# The results config align_results edits (deploy_prod_run.R's MANUSCRIPT_RESULTS).
RESULTS_CONFIG="analysis/configs/results/manuscript.yml"
CONCURRENCY=""; PARTITION=""; TIME_LIMIT=""; DEP="afterok"
RELAY_TIME="${CHOLERA_RELAY_TIME:-03:00:00}"; RESULTS_MEM="${CHOLERA_RESULTS_MEM:-160G}"
RESULTS_PARTITION=""  # "" = the cluster profile's (CHOLERA_RES_RESULTS, else SHORT)
EVENT_REPLICATES=""   # "" = each past-event config's simulator.n_replicates
AFTER=""              # job ids the first phase waits on (relay job); "" = run it now
GEN_FORCE=()          # (--force) with --force-gen: gen_templates / gen_attrib overwrite

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deploy)      DEPLOY="$2"; shift 2 ;;
    --from)        FROM="$2"; shift 2 ;;
    --until)       UNTIL="$2"; shift 2 ;;
    --phase)       ONLY="$2"; shift 2 ;;
    --concurrency) CONCURRENCY="$2"; shift 2 ;;
    --partition)   PARTITION="$2"; shift 2 ;;
    --time)        TIME_LIMIT="$2"; shift 2 ;;
    --dep)         DEP="$2"; shift 2 ;;
    --relay-time)  RELAY_TIME="$2"; shift 2 ;;
    --results-mem) RESULTS_MEM="$2"; shift 2 ;;
    --results-partition) RESULTS_PARTITION="$2"; shift 2 ;;
    --results-config) RESULTS_CONFIG="$2"; shift 2 ;;
    --event-replicates) EVENT_REPLICATES="$2"; shift 2 ;;
    --after)       AFTER="$2"; shift 2 ;;
    --plan)        PLAN=true; shift ;;
    --force-gen)   GEN_FORCE=(--force); shift ;;
    --local)       shift ;;   # read above (CHOLERA_BACKEND=local)
    --cores)       shift 2 ;; # read above (CHOLERA_CORES)
    -h|--help)     awk 'NR > 1 && /^set -euo pipefail/ { exit } NR > 1' "$0"; exit 0 ;;
    *) echo "Error: unknown option '$1' (see --help)" >&2; exit 1 ;;
  esac
done

[[ -n "$DEPLOY" && -f "$DEPLOY" ]] || { echo "Error: --deploy <prod_run<N>.yml> is required and must exist." >&2; exit 1; }
[[ "$DEP" == "afterok" || "$DEP" == "afterany" ]] || { echo "Error: --dep must be afterok or afterany." >&2; exit 1; }
[[ -z "$EVENT_REPLICATES" || "$EVENT_REPLICATES" =~ ^[1-9][0-9]*$ ]] || { echo "Error: --event-replicates must be a positive integer." >&2; exit 1; }
[[ -f "$RESULTS_CONFIG" ]] || { echo "Error: --results-config $RESULTS_CONFIG not found." >&2; exit 1; }
[[ -z "$AFTER" || "$AFTER" =~ ^[0-9]+(:[0-9]+)*$ ]] || { echo "Error: --after must be colon-joined job ids (e.g. 123:456)." >&2; exit 1; }
[[ -z "$AFTER" || -z "$ONLY" ]] || { echo "Error: --after and --phase are exclusive (relay jobs carry their own dependency)." >&2; exit 1; }

phase_index() { local i; for i in "${!PHASES[@]}"; do [[ "${PHASES[$i]}" == "$1" ]] && { echo "$i"; return; }; done
                echo "Error: unknown phase '$1' (one of: ${PHASES[*]})" >&2; exit 1; }

# Run paths from the deploy config (R owns them; deploy_prod_run.R --stage resolve).
RESOLVED=$(Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage resolve 2>/dev/null | grep -E '^[A-Z_]+=' || true)
# A key resolve does not print reads as "" (|| true: under pipefail a grep with
# no match would abort the `X=$(get_key ...)` assignment, and the script).
get_key() { echo "$RESOLVED" | grep "^$1=" | head -1 | cut -d= -f2- || true; }
DEPLOY_TAG=$(get_key DEPLOY_TAG)
HISTORICAL_TEMPLATE=$(get_key HISTORICAL_TEMPLATE)
[[ -n "$DEPLOY_TAG" && -n "$HISTORICAL_TEMPLATE" ]] || { echo "Error: could not resolve DEPLOY_TAG / HISTORICAL_TEMPLATE from $DEPLOY." >&2; exit 1; }
KL_SWEEP_BATCH=$(get_key KL_SWEEP_BATCH)
ATTRIB_MANIFEST=$(get_key FORECAST_ATTRIB_MANIFEST)
# Past-ENSO-event attribution: on when deploy.event_attribution.events is set.
EVENTS=$(get_key EVENTS)
EVENT_CONFIG_DIR=$(get_key EVENT_CONFIG_DIR)
EVENT_MANIFEST=$(get_key EVENT_MANIFEST)

# The whole chain is what the deploy declares: it starts at the class sweep
# when the deploy derives its grouping from one.
[[ -n "$FROM" ]] || { [[ -n "$KL_SWEEP_BATCH" ]] && FROM=classes || FROM=fits; }
for p in "$FROM" "$ONLY"; do
  [[ "$p" == classes && -z "$KL_SWEEP_BATCH" ]] && { echo "Error: phase 'classes' needs deploy.kl_sweep.batch in $DEPLOY." >&2; exit 1; }
  [[ "$p" == attrib && -z "$ATTRIB_MANIFEST" ]] && { echo "Error: phase 'attrib' needs deploy.outputs.forecast_attrib_manifest in $DEPLOY." >&2; exit 1; }
done
I_FROM=$(phase_index "$FROM"); I_UNTIL=$(phase_index "$UNTIL")
(( I_FROM <= I_UNTIL )) || { echo "Error: --from $FROM comes after --until $UNTIL." >&2; exit 1; }

CONC_PASS=()
[[ -n "$CONCURRENCY" ]] && CONC_PASS=(--concurrency "$CONCURRENCY")
FIT_PASS=("${CONC_PASS[@]}")
EVENT_PASS=("${CONC_PASS[@]}")
[[ -n "$EVENT_REPLICATES" ]] && EVENT_PASS+=(--n_replicates "$EVENT_REPLICATES")
[[ -n "$PARTITION"   ]] && FIT_PASS+=(--partition "$PARTITION")
[[ -n "$TIME_LIMIT"  ]] && FIT_PASS+=(--time "$TIME_LIMIT")
RESULTS_PASS=(--mem "$RESULTS_MEM")
[[ -n "$RESULTS_PARTITION" ]] && RESULTS_PASS+=(--partition "$RESULTS_PARTITION")

mkdir -p "$CHOLERA_LOG_DIR"
LOCAL=false
[[ "$CHOLERA_BACKEND" == local ]] && LOCAL=true

# Pull every job id a driver printed (sbatch's "Submitted batch job N" and the
# drivers' KEY=<id> / "... job id: <id>" lines), colon-joined for --dependency.
# Under `set -euo pipefail` a grep that matches nothing fails its pipeline, and a
# failing `ids=$(...)` aborts the script: every id capture below ends in
# `|| true`, so an empty capture reaches its explicit handling (the attrib
# phase's no-array branch; relay()'s "no job ids captured" error elsewhere).
collect_ids() {
  grep -oE '(Submitted batch job [0-9]+|^[A-Z_]*JOB_ID=[0-9]+|[Aa]rray job id: *[0-9]+|[Cc]ompile job id: *[0-9]+)' \
    | grep -oE '[0-9]+$' | paste -sd: -
}

# gen_batches only when a batch dir named by the deploy (resolve keys $@) has no
# manifest.yml yet: the stage loads make_batch_configs.R, whose packages
# (tidyverse) the cluster R library may lack, and existing dirs are skipped anyway.
need_batches() {
  local k d
  for k in "$@"; do d=$(get_key "$k"); [[ -n "$d" && ! -f "$d/manifest.yml" ]] && return 0; done
  return 1
}
gen_batches_if_needed() {
  if need_batches "$@"; then
    Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_batches
  else
    echo "[run_pipeline] batch dirs ($*) already generated: gen_batches skipped."
  fi
}

# Submit the next phase as a relay job depending on $1 (colon-joined ids).
relay() {
  local next="$1" ids="$2"
  local i_next; i_next=$(phase_index "$next")
  if (( i_next > I_UNTIL )); then echo "[run_pipeline] --until ${UNTIL}: chain stops here."; return; fi
  if [[ "$LOCAL" == true ]]; then
    # Local backend: every job of this phase has already run; no relay job.
    run_phase "$next"
    return
  fi
  [[ -n "$ids" ]] || { echo "Error: no job ids captured to chain phase '${next}' on; stopping." >&2; exit 1; }
  local args=(--deploy "$DEPLOY" --phase "$next" --until "$UNTIL" --dep "$DEP" --relay-time "$RELAY_TIME"
              --results-mem "$RESULTS_MEM" --results-config "$RESULTS_CONFIG" "${CONC_PASS[@]}")
  [[ -n "$EVENT_REPLICATES" ]] && args+=(--event-replicates "$EVENT_REPLICATES")
  [[ -n "$RESULTS_PARTITION" ]] && args+=(--results-partition "$RESULTS_PARTITION")
  [[ ${#GEN_FORCE[@]} -gt 0 ]] && args+=(--force-gen)
  local jid
  jid=$(cholera_sbatch SHORT --res RELAY --defaults "--mem=8G --cpus-per-task=1" \
          --parsable --job-name="pipe_${next}_${DEPLOY_TAG}" --time="$RELAY_TIME" \
          --dependency="${DEP}:${ids}" \
          --wrap="bash analysis/bash/tools/run_pipeline.sh $(printf '%q ' "${args[@]}")")
  echo "[run_pipeline] phase '${next}' queued as job ${jid} (${DEP}:${ids})"
  echo "PIPE_${next^^}_JOB_ID=${jid}"
}

run_phase() {
  local phase="$1" out ids
  echo "============================================================"
  echo "[run_pipeline] ${DEPLOY_TAG}: phase '${phase}'  ($(date '+%F %T'))"
  echo "============================================================"
  case "$phase" in
    classes)
      gen_batches_if_needed KL_SWEEP_DIR
      out=$(bash analysis/bash/stages/04_fit_kl.sh --deploy "$DEPLOY" --postprocess "${FIT_PASS[@]}")
      echo "$out"
      if [[ -z "$(get_key KL_SWEEP_CELL)" ]]; then
        echo "[run_pipeline] ${DEPLOY_TAG}: deploy.kl_sweep.cell is not pinned: the chain stops after the sweep."
        echo "  Read $(get_key KL_SWEEP_COMPARE_DIR)/{fits,loo}.csv once the comparison job ends,"
        echo "  set deploy.kl_sweep.cell (and deploy.fits.groups) in $DEPLOY, then rerun with --from fits."
        return
      fi
      # gen_grouping (start of `fits`) reads the comparison tables, itself afterok the fit array.
      ids=$(echo "$out" | grep '^KL_SWEEP_JOB_ID=' | collect_ids || true)
      relay fits "$ids" ;;
    fits)
      if [[ -n "$KL_SWEEP_BATCH" ]]; then
        Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_grouping
      fi
      gen_batches_if_needed OCC_DIR DUR_DIR MON_DIR
      out=""
      for c in enso_grouped_singlelag dur_enso_grouped_singlelag mon_enso_grouped_singlelag; do
        local o; o=$(bash "analysis/bash/campaigns/${c}.sh" --deploy "$DEPLOY" "${FIT_PASS[@]}")
        echo "$o"; out+="$o"$'\n'
      done
      ids=$(echo "$out" | collect_ids || true)
      relay postfit "$ids" ;;
    postfit)
      out=$(bash analysis/bash/campaigns/postfit.sh --deploy "$DEPLOY" --tracks occ,dur,mon "${CONC_PASS[@]}")
      echo "$out"
      # deploy_fits needs the LOO rankings (and the CB surfaces); the lag-compare
      # figures are not waited on, so a figure failure does not stop the chain.
      ids=$(echo "$out" | grep -E '^CB_JOB_ID=|[Aa]rray job id:' | collect_ids || true)
      relay deploy "$ids" ;;
    deploy)
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage deploy_fits
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_templates ${GEN_FORCE[@]+"${GEN_FORCE[@]}"}
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage smoothing
      out=$(bash analysis/bash/stages/08_smoothing.sh --config "$HISTORICAL_TEMPLATE" "${CONC_PASS[@]}")
      echo "$out"
      ids=$(echo "$out" | grep -E '[Cc]ompile job id:' | collect_ids || true)
      relay sim "$ids" ;;
    sim)
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_forecast
      # The attribution factual (no-op when deploy.outputs.forecast_attrib_manifest is "").
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_attrib ${GEN_FORCE[@]+"${GEN_FORCE[@]}"}
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_hindcast
      # Past ENSO events (09 report_event_attrib / report_event_af_model); off
      # when deploy.event_attribution.events is "" (gen_events refuses an empty
      # list, so it is not called then).
      if [[ -n "$EVENTS" ]]; then
        Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_events
        # Preflight before anything is submitted: --list refuses a counterfactual
        # config whose per-event Lenssen cache is missing.
        bash analysis/bash/stages/07_event_attribution.sh --deploy "$DEPLOY" --list
      else
        echo "[run_pipeline] ${DEPLOY_TAG}: deploy.event_attribution.events is unset: no past-event attribution."
      fi
      local fc hc ev=""
      fc=$(bash analysis/bash/stages/07_forecast_groups.sh --deploy "$DEPLOY" --arm all --postprocess "${CONC_PASS[@]}")
      echo "$fc"
      hc=$(bash analysis/bash/stages/07_hindcast.sh --deploy "$DEPLOY" --postprocess \
             --pp-workers "${CHOLERA_HINDCAST_PP_WORKERS:-8}" --pp-mem "${CHOLERA_HINDCAST_PP_MEM:-240G}" "${CONC_PASS[@]}")   # 240G: 64G runs out of memory on a full deploy
      echo "$hc"
      if [[ -n "$EVENTS" ]]; then
        ev=$(bash analysis/bash/stages/07_event_attribution.sh --deploy "$DEPLOY" "${EVENT_PASS[@]}")
        echo "$ev"
      fi
      echo "[run_pipeline] ${DEPLOY_TAG}: forecast / counterfactual, hindcast${EVENTS:+ and past events} submitted."
      # 09 needs the forecast arrays (parts/), the hindcast 07c scores (the
      # postprocess job, itself afterok every hindcast array) and the event
      # arrays (no event postprocess: 09 reads their partitions). The forecast
      # per-config postprocess lines ("postprocess job id:") are not collected.
      ids=$( { echo "$fc" | grep -E '[Aa]rray job id:'
               echo "$hc" | grep -E '[Aa]rray job id:|^Postprocess job id:' \
                 | sed -E 's/^Postprocess job id: *([0-9]+).*/Submitted batch job \1/'
               echo "$ev" | grep -E '[Aa]rray job id:' || true; } | collect_ids || true)
      echo "  Not chained: the cross-group combine (07c --combine-groups) and its figures (run by hand)."
      relay results "$ids" ;;
    attrib)
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage gen_attrib ${GEN_FORCE[@]+"${GEN_FORCE[@]}"}
      out=$(bash analysis/bash/stages/07_forecast_groups.sh --deploy "$DEPLOY" --arm attrib --postprocess "${CONC_PASS[@]}")
      echo "$out"
      echo "[run_pipeline] ${DEPLOY_TAG}: attribution factual stage done."
      ids=$(echo "$out" | grep -E '[Aa]rray job id:' | collect_ids || true)
      if [[ -z "$ids" ]]; then
        # Nothing submitted: the arm is complete (or the launches are aligned);
        # 09 can start at once -- unless --until stops the chain here.
        if (( $(phase_index results) > I_UNTIL )); then
          echo "[run_pipeline] --until ${UNTIL}: chain stops here (no attribution array was needed)."
          return
        fi
        echo "[run_pipeline] ${DEPLOY_TAG}: no attribution array to wait on; running results now."
        run_phase results
        return
      fi
      relay results "$ids" ;;
    results)
      Rscript analysis/utils/deploy_prod_run.R -C "$DEPLOY" --stage align_results --results-config "$RESULTS_CONFIG"
      out=$(bash analysis/bash/stages/09_postprocess_results.sh --config "$RESULTS_CONFIG" "${RESULTS_PASS[@]}" "${CONC_PASS[@]}")
      echo "$out"
      echo "[run_pipeline] ${DEPLOY_TAG}: 09 waves submitted."
      # 10 + 11 read every 09 artefact: wait on all waves (RESULTS_JOB_IDS=w1:w2:...).
      ids=$(echo "$out" | grep '^RESULTS_JOB_IDS=' | tail -1 | cut -d= -f2- || true)
      relay manuscript "$ids" ;;
    manuscript)
      out=$(bash analysis/bash/stages/10_manuscript_figures.sh --config "$RESULTS_CONFIG" \
              --job-name "manuscript_${DEPLOY_TAG}")
      echo "$out"
      id=$(echo "$out" | grep '^MANUSCRIPT_JOB_ID=' | tail -1 | cut -d= -f2- || true)
      echo "[run_pipeline] ${DEPLOY_TAG}: 10 + 11 submitted as job ${id}; end of the chain." ;;
  esac
}

if [[ "$PLAN" == "true" ]]; then
  KL_CELL=$(get_key KL_SWEEP_CELL)
  echo "[run_pipeline] ${DEPLOY_TAG} (${DEPLOY}): phases ${FROM} .. ${UNTIL}, chained ${DEP}"
  [[ -n "$AFTER" ]] && echo "  start   : phase '${FROM}' queued as a relay job, ${DEP}:${AFTER}"
  for ((i = I_FROM; i <= I_UNTIL; i++)); do
    case "${PHASES[$i]}" in
      classes) echo "  classes : deploy_prod_run.R --stage gen_batches (if the sweep dir is missing); stages/04_fit_kl.sh --deploy $DEPLOY --postprocess ${FIT_PASS[*]} (sweep $KL_SWEEP_BATCH, cell ${KL_CELL:-<not pinned: chain stops here>})" ;;
      fits)    echo "  fits    : deploy_prod_run.R --stage ${KL_SWEEP_BATCH:+gen_grouping, }gen_batches (if a fit dir is missing); campaigns/{enso,dur_enso,mon_enso}_grouped_singlelag.sh --deploy $DEPLOY ${FIT_PASS[*]}" ;;
      postfit) echo "  postfit : campaigns/postfit.sh --deploy $DEPLOY --tracks occ,dur,mon ${CONC_PASS[*]}" ;;
      deploy)  echo "  deploy  : deploy_prod_run.R --stage deploy_fits, gen_templates, smoothing; stages/08_smoothing.sh --config $HISTORICAL_TEMPLATE ${CONC_PASS[*]}" ;;
      sim)     echo "  sim     : deploy_prod_run.R --stage gen_forecast, gen_attrib${ATTRIB_MANIFEST:+ (-> $ATTRIB_MANIFEST)}, gen_hindcast${EVENTS:+, gen_events (-> $EVENT_MANIFEST)}; stages/07_forecast_groups.sh --arm all --postprocess ${CONC_PASS[*]}; stages/07_hindcast.sh --postprocess ${CONC_PASS[*]}"
               if [[ -n "$EVENTS" ]]; then
                 echo "            events: stages/07_event_attribution.sh --deploy $DEPLOY --list (preflight), then stages/07_event_attribution.sh --deploy $DEPLOY ${EVENT_PASS[*]} (events $EVENTS; configs $EVENT_CONFIG_DIR; complete runs skipped)"
                 if [[ -n "$EVENT_REPLICATES" ]]; then
                   echo "            event replicates: $EVENT_REPLICATES (--event-replicates)"
                 else
                   echo "            event replicates: each event config's simulator.n_replicates (no --event-replicates; payload fallback 100)"
                 fi
                 echo "            results waits afterok on the forecast arrays, the hindcast postprocess and the event arrays"
               else
                 echo "            events: (skipped: deploy.event_attribution.events is unset)"
               fi ;;
      attrib)  if (( I_FROM == $(phase_index attrib) )); then
                 echo "  attrib  : deploy_prod_run.R --stage gen_attrib (-> $ATTRIB_MANIFEST); stages/07_forecast_groups.sh --arm attrib --postprocess ${CONC_PASS[*]}"
               else
                 echo "  attrib  : (skipped: the sim phase submits the attribution factual)"
               fi ;;
      results) echo "  results : deploy_prod_run.R --stage align_results; stages/09_postprocess_results.sh --config $RESULTS_CONFIG ${RESULTS_PASS[*]} ${CONC_PASS[*]}" ;;
      manuscript) echo "  manuscript: stages/10_manuscript_figures.sh --config $RESULTS_CONFIG (10 + 11, one job)" ;;
    esac
  done
  exit 0
fi

if [[ -n "$ONLY" ]]; then
  phase_index "$ONLY" >/dev/null
  run_phase "$ONLY"
elif [[ -n "$AFTER" ]]; then
  relay "$FROM" "$AFTER"
else
  run_phase "$FROM"
fi
