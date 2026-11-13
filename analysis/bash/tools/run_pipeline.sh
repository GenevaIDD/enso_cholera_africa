#!/bin/bash
# tools/run_pipeline.sh — run a deploy (configs/deploy/<run>.yml) end to end on a SLURM
# cluster (the cluster profile, lib/cluster.sh) or this machine (--local) as
# a chain of jobs: [kl class sweep ->] fits -> postfit -> deploy +
# smoothing -> forecast / counterfactual / attribution factual + hindcast +
# past ENSO events -> manuscript results (09) -> manuscript figures and key
# numbers (10 + 11).
#
# A relay. Each PHASE submits its jobs with the existing drivers, then submits
# the NEXT phase as a small SLURM job (this script again, `--phase <next>`)
# with --dependency=afterok on the jobs the next phase needs. A later phase
# cannot be submitted up front because its drivers build their task lists from
# files the earlier jobs write (LOO rankings, deployed fits, smoothing cache).
# --local (or CHOLERA_BACKEND=local): no relay; every phase runs inline, one
# after the other, its jobs run by lib/local/sbatch on this machine
# ($CHOLERA_CORES array tasks at a time; tools/run_synthetic_pipeline.sh).
# Needs `sbatch` on the compute nodes (the relay jobs submit the next phase).
#
#   classes  (only when the deploy sets kl_sweep.batch) deploy_prod_run.R
#            --stage gen_batches (the sweep dir; the fit dirs wait for the
#            registry), then stages/04_fit_kl.sh --deploy D --postprocess (the
#            brms_kl class sweep + 05_fit_postprocess_kl_sweep.R)
#            -> next phase afterok the sweep comparison job, if
#               deploy.kl_sweep.cell is pinned; otherwise the chain stops here:
#               read generated_data/kl_sweep_compare/<batch>/{fits,loo}.csv,
#               pin the cell in the deploy YAML and rerun with --from fits
#   fits     [deploy_prod_run.R --stage gen_grouping: analysis/groupings/
#            <deploy.grouping>.{csv,yml} from the pinned cell, when kl_sweep.batch
#            is set], deploy_prod_run.R --stage gen_batches (only when a fit dir
#            has no manifest.yml: the batch maker needs tidyverse),
#            then campaigns/{enso,dur_enso,mon_enso}_grouped_singlelag.sh --deploy D
#            -> next phase afterok all fit arrays
#   postfit  campaigns/postfit.sh --deploy D --tracks occ,dur,mon
#            (05 cb surfaces + LOO lag rankings, then 06 lag-compare figures)
#            -> next phase afterok the CB + LOO arrays (the figures job is not waited on)
#   deploy   deploy_prod_run.R --stage deploy_fits, gen_templates, smoothing (check),
#            then stages/08_smoothing.sh --config <historical template>
#            -> next phase afterok the smoothing compile job
#   sim      deploy_prod_run.R --stage gen_forecast, gen_attrib (no-op unless
#            deploy.outputs.forecast_attrib_manifest is set), gen_hindcast,
#            gen_events (only when deploy.event_attribution.events is set), then
#            [stages/07_event_attribution.sh --deploy D --list: a preflight that
#            refuses a missing per-event Lenssen cache before anything is
#            submitted], stages/07_forecast_groups.sh --deploy D --arm all
#            --postprocess (the factual, the attribution factual and the
#            counterfactual arms, in parallel; complete runs are skipped),
#            stages/07_hindcast.sh --deploy D --postprocess and
#            [stages/07_event_attribution.sh --deploy D (the past-ENSO-event
#            factual + counterfactual arrays; complete runs are skipped;
#            --event-replicates N passes --n_replicates N, else each config's
#            simulator.n_replicates, set by deploy.event_attribution.n_replicates)]
#            -> next phase (results; `attrib` is skipped, its arm is submitted
#               here) afterok the forecast arrays, the hindcast postprocess job
#               and the event arrays (09 reads the forecast parts/ directly, the
#               hindcast 07c scores and the event partitions (report_event_attrib,
#               report_event_af_model); the per-run forecast 07c/07d jobs are not
#               waited on, so a per-run figure failure does not stop the chain)
#   attrib   entered only with --from attrib / --phase attrib: adds the
#            attribution factual (the factual launched at the counterfactual
#            handoff)
#            to a rollout whose sim phase already ran. deploy_prod_run.R --stage
#            gen_attrib, then stages/07_forecast_groups.sh --deploy D --arm attrib
#            --postprocess. Needs deploy.outputs.forecast_attrib_manifest.
#            -> next phase afterok the attribution arrays
#   results  deploy_prod_run.R --stage align_results (points
#            analysis/configs/results/manuscript.yml at this deploy: a text edit
#            of results.deploy_config + the deploy-derived inputs.*, comments kept;
#            inputs.attrib_manifest switches the attribution steps to the
#            attribution factual),
#            then stages/09_postprocess_results.sh --config <that file>: one
#            09b array per dependency wave, chained afterok
#            -> next phase afterok every 09 wave
#   manuscript  stages/10_manuscript_figures.sh --config <that file> (slurm/submit_manuscript.sh):
#            package check, 10_manuscript_figures.R, quarto render 11
#            -> next phase afterok the 10 + 11 job
#
# This is a deliberate exception to the "no umbrella entry point" rule for
# campaigns/ (analysis/bash/README.md): it adds no logic of its own, it only
# chains the existing drivers in the order of analysis/PIPELINE.md.
#
# Usage (D=analysis/configs/deploy/prod_run<N>.yml), from the repo root on the login node:
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D                  # whole chain (from classes when kl_sweep.batch is set, else fits)
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --until classes  # fit + compare the kl sweep only
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from fits      # after pinning kl_sweep.cell
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from postfit   # restart at a phase
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --until deploy   # stop after a phase
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --until sim      # everything but 09
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --event-replicates 9950  # past-event runs at an explicit replicate count
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --until results  # everything but 10 + 11
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from attrib    # add the attribution factual to a finished sim phase, then 09 + 10 + 11
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from results   # align_results + 09, then 10 + 11
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --from manuscript  # 10 + 11
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --plan           # print the phases, submit nothing
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --after 123:456  # first phase queued after those jobs (e.g. stages/03_prepare.sh's PREPARE_JOB_IDS)
#   bash analysis/bash/tools/run_pipeline.sh --deploy $D --results-config analysis/configs/results/manuscript_run<N>.yml --plan
#                                         # a run under review: 09 / 10 / 11 on its own results config
#
# Options:
#   --deploy FILE         Deploy config (REQUIRED); every path comes from it.
#   --from PHASE          First phase: classes | fits | postfit | deploy | sim | attrib | results | manuscript
#                         [default: classes when the deploy sets kl_sweep.batch, else fits]
#   --until PHASE         Last phase [default: manuscript]
#   --phase PHASE         Run exactly this phase, then relay (used by the relay jobs)
#   --concurrency K       Cap simultaneous array slots in every phase: fits, postfit
#                         (cb / LOO / figures), smoothing, the forecast / hindcast / event
#                         arrays (whole-unit one-task arrays there: no effect) and
#                         the 09 waves.
#                         Carried through the relay jobs. Unset = each driver's default
#                         (1000 in the drivers that cap: 02 population, 07 v1, 08, 09).
#   --partition P         Fits: SLURM partition(s), passed to the fit campaigns and
#                         the kl sweep (04_fit_kl.sh)
#   --time T              Fits: time limit, passed to the fit campaigns and the kl sweep
#   --dep TYPE            Dependency type between phases: afterok | afterany [default: afterok].
#                         With afterok a single failed array task stops the chain
#                         (SLURM leaves the next phase pending as DependencyNeverSatisfied);
#                         afterany carries on regardless.
#   --relay-time T        Time limit of each relay job [default: $CHOLERA_RELAY_TIME of the
#                         cluster profile, else 03:00:00]. Its memory / cores: $CHOLERA_RES_RELAY,
#                         else --mem=8G --cpus-per-task=1.
#   --event-replicates N  Sim: replicates of every past-event run, passed to
#                         stages/07_event_attribution.sh as --n_replicates N
#                         [default: none: each event config's simulator.n_replicates
#                         (deploy.event_attribution.n_replicates, 0 = the forecast
#                         template's, e.g. 9950), else the payload's 100].
#                         Prefer the deploy key: the configs then record the count.
#   --results-mem M       Memory of every 09 array task [default: $CHOLERA_RESULTS_MEM of the
#                         cluster profile, else 160G]: the
#                         attribution steps load the whole forecast + counterfactual
#                         ADM1-month runs (~10,000 replicates: ~100 GB peak, above
#                         slurm/submit_results.sh's 64G).
#   --results-partition P Partition(s) of every 09 array task [default: the cluster
#                         profile's (CHOLERA_RES_RESULTS, else its SHORT class)].
#                         Carried through the relay jobs.
#   --results-config F    Results config the results / manuscript phases
#                         align and read [default: analysis/configs/results/manuscript.yml].
#                         A run reviewed before manuscript.yml moves to it (e.g.
#                         analysis/configs/results/manuscript_<run>.yml) passes its own.
#   --after IDS           Colon-joined Slurm job ids (e.g. the 03 prep jobs of
#                         stages/03_prepare.sh): the first phase is not run now but
#                         queued as a relay job (--dep type, default afterok) on them.
#   --force-gen           Pass --force to deploy_prod_run.R --stage gen_templates and
#                         gen_attrib: overwrite generated templates that differ (a rerun
#                         whose fits changed, e.g. tools/run_synthetic_pipeline.sh).
#   --plan                Print what each phase would run and exit.
#   --local               Run on this machine (CHOLERA_BACKEND=local): every phase inline,
#                         no relay jobs (--after: the jobs have run, the first phase runs now).
#   --cores N             --local: parallel array tasks [default: $CHOLERA_CORES, else nproc].
#   The hindcast postprocess of the sim phase gets --pp-workers $CHOLERA_HINDCAST_PP_WORKERS
#   (default 8) --pp-mem $CHOLERA_HINDCAST_PP_MEM (default 240G; 64G runs out of memory on a
#   full deploy). Cluster profile: ${CHOLERA_CLUSTER_PROFILE:-analysis/bash/cluster.profile}.
#
# Restarting: after fixing a failure, rerun with --from <phase>. The drivers
# are idempotent where they can be (fits/postprocess skip existing outputs
# under overwrite: no; smoothing has --resume, see stages/08_smoothing.sh;
# the forecast arms and the past-event runs skip complete runs: --from sim
# resubmits only the forecast / event configs whose run lacks a _SUCCESS
# sentinel in some partition; 07_forecast_groups.sh / 07_event_attribution.sh
# --overwrite, run by hand, forces a resubmit; run_pipeline.sh never passes it).
# Every phase logs to logs/pipe_<phase>_<deploy tag>_<jobid>.log.

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
