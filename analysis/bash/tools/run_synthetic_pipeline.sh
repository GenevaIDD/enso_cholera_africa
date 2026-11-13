#!/bin/bash
# tools/run_synthetic_pipeline.sh — the whole pipeline on the synthetic world, locally.
#
# Builds a sandbox (its own data/, generated_data/, figures/ and a fake
# outbreak-data root), writes the synthetic world's raw inputs
# (analysis/00_ingest_synthetic.R: fake countries, analysis/groupings/
# synthetic_countries.csv), builds the panels (the data phase, below), then
# hands over to the production chain on this machine:
#   tools/run_pipeline.sh --deploy configs/deploy/synth_run1.yml --local --from fits
# (the same stage drivers and slurm/ payloads as on a cluster, run by
# lib/local/sbatch; CHOLERA_BACKEND=local). No restricted data, no cluster.
#
# Phases (--from / --until):
#   data       00 synthetic world -> 00 ingest -> 01 -> 02 population, crosswalk,
#              climate -> 03 panels, ENSO continuation library, WHO monthly, rainfall
#   fits .. manuscript   run_pipeline.sh's phases (see its header): fits, postfit,
#              deploy (smoothing), sim (forecast / attribution / counterfactual arms,
#              hindcast, past ENSO events), results (09, results config
#              manuscript_synth.yml), manuscript (10 + 11). The reports phase (12)
#              is not run.
#
# Usage (from the repository root):
#   bash analysis/bash/tools/run_synthetic_pipeline.sh                      # everything
#   bash analysis/bash/tools/run_synthetic_pipeline.sh --until data         # inputs + panels only (minutes)
#   bash analysis/bash/tools/run_synthetic_pipeline.sh --from sim --cores 8
#
# Options:
#   --sandbox DIR   sandbox root [default: synthetic_world]; created, `analysis`
#                   symlinked into it; every step runs with DIR as working dir
#   --from PHASE    first phase [default: data]
#   --until PHASE   last phase  [default: manuscript]
#   --cores N       parallel processes (data phase workers; CHOLERA_CORES, the
#                   array tasks run at once by lib/local/sbatch) [default: 4]
#   --seed N        seed of the synthetic world [default: 1]
#
# Runtime on a 4-core laptop: data ~2 min, fits ~10 min, the rest ~minutes
# per arm. Logs: <sandbox>/logs/synthetic/<step>.log (the data steps and
# run_pipeline.log), one per job in <sandbox>/logs/synthetic/jobs/.

set -euo pipefail

SANDBOX="synthetic_world"; FROM="data"; UNTIL="manuscript"; CORES=4; SEED=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sandbox) SANDBOX="$2"; shift 2 ;;
    --from)    FROM="$2"; shift 2 ;;
    --until)   UNTIL="$2"; shift 2 ;;
    --cores)   CORES="$2"; shift 2 ;;
    --seed)    SEED="$2"; shift 2 ;;
    -h|--help) sed -n '2,37p' "$0"; exit 0 ;;
    *) echo "Error: unknown option $1" >&2; exit 1 ;;
  esac
done

PHASES=(data fits postfit deploy sim results manuscript)
idx() { local i; for i in "${!PHASES[@]}"; do [[ "${PHASES[$i]}" == "$1" ]] && { echo "$i"; return; }; done
        echo "Error: unknown phase '$1' (${PHASES[*]})" >&2; exit 1; }
I_FROM=$(idx "$FROM"); I_UNTIL=$(idx "$UNTIL")
want() { local i; i=$(idx "$1"); (( i >= I_FROM && i <= I_UNTIL )); }

REPO=$(pwd)
[[ -f "$REPO/analysis/00_ingest_synthetic.R" ]] || { echo "Error: run from the repository root." >&2; exit 1; }
mkdir -p "$SANDBOX"
cd "$SANDBOX"
[[ -e analysis ]] || ln -s "$REPO/analysis" analysis
mkdir -p data generated_data figures logs/synthetic secure
touch .synthetic_sandbox   # 00_ingest_synthetic.R refuses to run anywhere else
export CHOLERA_OUTBREAK_DATA_DIR="$PWD/secure/cholera_outbreak_data"
export CHOLERA_PROJECT_ROOT="$PWD"   # 11's root (analysis/ is a symlink into the repository)
LOG=logs/synthetic
D=analysis/configs/deploy/synth_run1.yml
R=analysis/configs/results/manuscript_synth.yml
B=1990_2019
C=(-b "$B" --episode_key gid_season)

step() {   # step NAME CMD...: run, log, stop on failure
  local name=$1; shift
  local t0=$SECONDS
  if "$@" > "$LOG/$name.log" 2>&1; then
    printf "[synthetic] %-28s ok   (%ds)\n" "$name" $(( SECONDS - t0 ))
  else
    printf "[synthetic] %-28s FAIL (%ds) -- tail of %s:\n" "$name" $(( SECONDS - t0 )) "$LOG/$name.log" >&2
    tail -25 "$LOG/$name.log" >&2
    exit 1
  fi
}

# --- data ---------------------------------------------------------------------
if want data; then
  POP=analysis/configs/pop_extract/pop_extract_synth.yml
  step 00_synthetic_world    Rscript analysis/00_ingest_synthetic.R --seed "$SEED" --force
  step 00_ingest             Rscript analysis/00_ingest_outbreakextractr.R
  step 01_outbreaks          Rscript analysis/01_preprocess_outbreaks.R -r TRUE -f FALSE
  step 01b_monthly           Rscript analysis/01b_preprocess_outbreaks_monthly.R
  step 01_national_zero      Rscript analysis/01_preprocess_national_zero_months.R --rule r3 --zero_threshold 60
  step 02a_population_tasks  Rscript analysis/02a_preprocess_population_tasks.R -C "$POP" --task_unit country
  for t in generated_data/population/tasks/*.yml; do
    step "02b_population_$(basename "$t" .yml)" Rscript analysis/02b_preprocess_population_run_country.R -C "$POP" --task "$t" --workers "$CORES"
  done
  step 02c_population        Rscript analysis/02c_preprocess_population_compile.R -C "$POP"
  step 02_crosswalk          Rscript analysis/02_preprocess_crosswalk.R -C "$POP"
  step 02_climate            Rscript analysis/02_preprocess_climate.R --stage all -s TRUE --workers "$CORES"
  P=enso_grouping_synth,hmm_clusters_synth
  step 03_occ_default        Rscript analysis/03_prepare_occ.R "${C[@]}" --contagion_partition_source "$P" \
                               --enso_max_lead 11 --redo TRUE --force
  step 03_occ_deploy         Rscript analysis/03_prepare_occ.R "${C[@]}" --contagion_partition_source "$P" \
                               --enso_max_lead 11 --national_zero_fill --national_zero_threshold 60 --national_zero_rule r3 \
                               --at_risk_start after_episode --adm1_overlap censor_orphan --contagion_coverage_report \
                               --redo TRUE --force
  step 03_siz                Rscript analysis/03_prepare_siz.R "${C[@]}" --force
  step 03_dur                Rscript analysis/03_prepare_dur.R "${C[@]}" --enso_max_lead 11 --force
  step 03_mon                Rscript analysis/03_prepare_mon.R "${C[@]}" --force
  step 03_enso_library       Rscript analysis/03_prepare_enso_continuation_library.R
  step 03_enso_library_skill Rscript analysis/03_prepare_enso_continuation_library.R --validate --deploy "$D"
  step 03_sim_panel          Rscript analysis/03_prepare_sim_panel.R "${C[@]}" --enso_max_lead 11 --no_refresh_enso --force
  step 03_sim_projection     Rscript analysis/03_prepare_sim_panel_projection.R "${C[@]}" --enso_max_lead 11 \
                               --no_refresh_enso --cfsv2_correction sd \
                               --enso_continuation analog --enso_continuation_library generated_data/enso_continuation_library.qs \
                               --enso_continuation_k 30 --enso_continuation_blend_months 3 --enso_continuation_min_sep 24 \
                               --enso_continuation_max_per_model Inf --enso_continuation_seam last_full_oni --force
  step 03_sim_extended       Rscript analysis/03_prepare_sim_panel_extended.R --episode_key gid_season \
                               --enso_max_lead 11 --end 2026-08-01 --force
  step 03_who_monthly        Rscript analysis/03_prepare_who_monthly.R --zero_iso3 "" --force
  step 03_rain_adm1          Rscript analysis/03_prepare_rain_adm1.R
fi

# --- fits .. manuscript: the production chain, on this machine ---------------------
# --force-gen: a rerun whose fits moved a lag regenerates the synthetic templates.
if (( I_UNTIL >= $(idx fits) )); then
  PIPE_FROM="${PHASES[$I_FROM]}"
  [[ "$PIPE_FROM" == data ]] && PIPE_FROM=fits
  export CHOLERA_LOG_DIR="$LOG/jobs"
  step run_pipeline bash analysis/bash/tools/run_pipeline.sh --deploy "$D" --local --cores "$CORES" \
    --from "$PIPE_FROM" --until "$UNTIL" --results-config "$R" --force-gen
  grep -E '^\[run_pipeline\]' "$LOG/run_pipeline.log" | sed 's/^/  /'
fi

echo "[synthetic] done: phases ${PHASES[$I_FROM]} .. ${PHASES[$I_UNTIL]} in $SANDBOX (logs in $SANDBOX/$LOG)"
