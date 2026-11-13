#!/usr/bin/env Rscript
# analysis/09b_postprocess_results.R
# =============================================================================
# Phase 09 — MANUSCRIPT results postprocessing (derivation only).
#
# Derives the aggregated tables analysis/10_manuscript_figures.R plots, and
# NOTHING else: no ggplot object is built here; the derivation/plotting seam
# is explained in analysis/utils/results_core.R's header.
#
# Inputs come from the config-driven `results:` block, caching goes through
# results_core.R's makeResultsRunFile()/results_read_if_fresh()/results_write()
# contract, and each derive_*() returns its table rather than writing a cache
# path itself.
#
# CACHING CONTRACT
# -----------------
#   "cache" = DERIVED aggregations only. This script NEVER re-fits a model
#   and NEVER re-runs the bpfilter / smoothing. It reads raw outputs already
#   on disk (brms fits + cb_meta sidecars, the 08b smoothing part.qs
#   partitions, the 07b forecast/counterfactual monthly_summary.qs
#   partitions, ENSO caches) and writes only the small derived tables to
#   generated_data/results/<tag>/. Each artefact is stamped with the decisions
#   (results_core.R: results_write()) that produced it and is invalidated by
#   either a stale input (mtime) OR a changed decision (the stamp diff) — see
#   results_core.R's freshness-contract section header.
#
# Usage
# -----
#   Rscript analysis/09b_postprocess_results.R -C analysis/configs/results/manuscript.yml
#   Rscript analysis/09b_postprocess_results.R -C ... --redo_cache TRUE
#   Rscript analysis/09b_postprocess_results.R -C ... --list                 # steps, waves, artefact status
#   Rscript analysis/09b_postprocess_results.R -C ... --steps oni,enso_relation:occurrence
#   Rscript analysis/09b_postprocess_results.R -C ... --task <tasks_dir>/task_0003.yml
#
# On the HPC: bash analysis/bash/stages/09_postprocess_results.sh --config <cfg>
# (09a_postprocess_results_tasks.R writes one task per step; each dependency
# wave is one SLURM array, chained afterok).
#
# Structure
# ---------
#   The derivers live in analysis/utils/results_derive.R. This file is the CLI:
#   parse options, build the read-only context (results_context()), derive.
#   With no step selection it runs every step (results_dispatch()); --steps or
#   --task run a subset of the step registry (.RESULTS_STEP_TABLE,
#   results_run_steps()). Every artefact is identical whichever way it is
#   derived.
#   Everything runs behind `sys.nframe() == 0L`, so source()-ing this file has
#   no side effects.
# =============================================================================

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(optparse)
    library(yaml)
    library(tibble)
    library(tidyr)
    library(readr)
    library(purrr)
    library(dplyr)
    library(stringr)
    library(lubridate)
    library(qs2)
    library(sf)
    library(brms)
    # Attach rstan up front: left to load lazily, S4 dispatch on a fit's stanfit
    # attaches it mid-dispatch, which fails on R 4.3.2
    # ("Unable to load module class_model_base: attempt to apply non-function").
    library(rstan)
    library(dlnm)
  })

  # Verified closure (grepped against every function call in the derivation
  # half): results_core carries the track/group registries, the freshness +
  # provenance contract, and the manifest/config accessors. gid/sim_io/
  # simulator_core/figures_cb/loo_lags are each called by the derivers in utils/results_derive.R
  # (build_location_frame, dropOverlappingADM1FromFrame, build_coverage_mask,
  # flattenGid, build_onebasis_grid/erf_from_draws/get_re_vec/get_re_draws_mat,
  # select_best_lag). The Phase-11 "key numbers" report adds population.R
  # (loadGidPopulation/resolveGidPopulation), dur_post.R (durationRatioDraws/
  # durationRatioSummary), and dlnm.R (build_cb_columns_from_Q, for the
  # construction-faithful duration RMD newdata). plots.R adds getCovarType(), the
  # classifier results_core's param_family() wraps for derive_track_params();
  # plots.R's own closure (clim_compile for getPrecipVars/getTempVars, clim_db)
  # comes with it. Deliberately NOT spatpomp_core.R: nothing here may be able to
  # invoke a particle filter. deploy_core.R is sourced so that a non-empty
  # results.deploy_config resolves (make_options_from_config() hard-stops if the
  # pointer is set but resolve_results_inputs_from_deploy() is absent — the
  # "dormant branch / missing source" trap); it is inert when deploy_config is "".
  # obs_sim_compare / sim_metrics / hindcast_skill: the hindcast-skill section
  # (ensure_hindcast(); scoring in utils/hindcast_skill.R). sim_metrics is pure
  # at load; 09 calls only load_observed_from_panel() from it.
  # anchor / diagnostics: the supplementary tables (report_anchor reads the
  # baseline-rate tables; convergence runs fit_diagnostics() on the deployed fits).
  # who_annual: report_who_enso_ratio's WHO annual cases (buildWhoAnnualCases(),
  # whoWindowTotals(); supplement fig. WR).
  # season_calendar: group_labels' per-country season ends (table S30);
  # kl_sweep: report_membership's adjusted_rand_index() / kl_match_classes().
  for (.m in c("results_core", "gid", "sim_io", "coverage", "simulator_core", "figures_cb", "loo_lags",
               "population", "dur_post", "dlnm", "plots", "clim_compile", "clim_db", "rain_enso",
               "deploy_core", "obs_sim_compare", "sim_metrics", "hindcast_skill",
               "anchor", "diagnostics", "season_calendar", "kl_sweep",
               "who_dashboard", "who_annual"))
    source(file.path("analysis/utils", paste0(.m, ".R")))
  source("analysis/utils/results_derive.R")   # results_context() / results_dispatch() + every deriver

  opt <- make_options_from_config(
    option_list = list(
      make_option("--redo_cache", type = "logical", default = FALSE,
                  help = "Re-derive every artefact from raw outputs regardless of freshness [default: %default]"),
      make_option("--workers", type = "integer", default = NULL,
                  help = "Forked workers for the per-group fit loads and per-run hindcast scoring; overrides results.workers. Artefacts are identical for any value."),
      make_option("--steps", type = "character", default = NULL,
                  help = "Comma-separated step ids to run (see --list), e.g. oni,enso_relation:occurrence [default: every step]"),
      make_option("--task", type = "character", default = NULL,
                  help = "A task YAML written by 09a_postprocess_results_tasks.R: run its one step (SLURM array payload)"),
      make_option("--list", action = "store_true", default = FALSE,
                  help = "Print the steps (id, wave, kinds, artefact status) and exit")),
    model = NULL)

  if (!is.null(opt[["task"]]) && !is.null(opt[["steps"]]))
    stop("[09b] pass --steps or --task, not both.", call. = FALSE)
  step_ids <- if (!is.null(opt[["steps"]])) strsplit(opt[["steps"]], ",", fixed = TRUE)[[1]]
  if (!is.null(opt[["task"]])) {
    if (!file.exists(opt[["task"]])) stop("[09b] task file not found: ", opt[["task"]], call. = FALSE)
    task <- yaml::read_yaml(opt[["task"]])
    cfg_path <- results_config_path()
    if (!identical(normalizePath(task$config, mustWork = FALSE),
                   normalizePath(cfg_path, mustWork = FALSE)))
      stop("[09b] task ", opt[["task"]], " was written for config ", task$config,
           ", not ", cfg_path, ". Re-run 09a_postprocess_results_tasks.R.", call. = FALSE)
    step_ids <- task$step
    # The task records 09a's --redo_cache; the 09b flag can only add to it.
    if (isTRUE(task$redo_cache)) opt[["redo_cache"]] <- TRUE
  }

  ctx <- results_context(opt)
  if (isTRUE(opt[["list"]])) {
    options(width = 200)
    print(results_steps_table(ctx, step_ids), row.names = FALSE)
  } else if (is.null(step_ids)) {
    results_dispatch(ctx)
  } else {
    results_run_steps(ctx, step_ids)
    message("[09b] done: ", paste(step_ids, collapse = ", "))
  }
  # Exit here, explicitly. --workers > 1 forks (base.R::mapWorkers() ->
  # parallel::mclapply); a child the parallel package fails to terminate shares
  # this script's open file offset, so on falling off the end Rscript read on
  # from a moved offset and died with `unexpected '[' in "["` -- exit 1 after
  # every artefact was written, which fails a SLURM afterok chain.
  quit(save = "no", status = 0)
}
