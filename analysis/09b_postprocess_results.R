#!/usr/bin/env Rscript
# 09b_postprocess_results.R — derive the manuscript results tables (no plots)
#
# Usage:
#   Rscript analysis/09b_postprocess_results.R \
#     -C analysis/configs/results/manuscript.yml
#   Rscript analysis/09b_postprocess_results.R -C <cfg> --list
#   Rscript analysis/09b_postprocess_results.R -C <cfg> \
#     --steps oni,enso_relation:occurrence
#   Rscript analysis/09b_postprocess_results.R -C <cfg> \
#     --task <tasks_dir>/task_0003.yml
#
# Details:
#   Derives the tables 10_manuscript_figures.R plots and nothing else; inputs
#   come from the config's results: block. Never re-fits a model or re-runs a
#   filter; writes only stamped derived tables to generated_data/results/<tag>/.
#   No selection runs every step; --steps / --task run a subset (identical
#   artefacts). --redo_cache re-derives regardless of freshness. On the HPC:
#   bash analysis/bash/stages/09_postprocess_results.sh (one array per wave).
#   Everything runs behind sys.nframe() == 0L: source() has no side effects.
#
# See: analysis/utils/results_core.R (caching, freshness, derive/plot seam)

# 1. Preamble ------------------------------------------------------------------
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

  # Module closure of the derivation half (who calls what:
  # analysis/utils/README.md). Deliberately NOT spatpomp_core.R: nothing here
  # may run a particle filter. deploy_core.R is sourced so a non-empty
  # results.deploy_config resolves (else make_options_from_config() hard-stops:
  # the dormant-branch trap); it is inert when deploy_config is "".
  for (.m in c(
    "results_core",
    "gid",
    "sim_io",
    "coverage",
    "simulator_core",
    "figures_cb",
    "loo_lags",
    "population",
    "dur_post",
    "dlnm",
    "plots",
    "clim_compile",
    "clim_db",
    "rain_enso",
    "deploy_core",
    "obs_sim_compare",
    "sim_metrics",
    "hindcast_skill",
    "anchor",
    "diagnostics",
    "season_calendar",
    "kl_sweep",
    "who_dashboard",
    "who_annual"
  )) {
    source(file.path("analysis/utils", paste0(.m, ".R")))
  }
  source("analysis/utils/results_derive.R") # results_context() / results_dispatch() + every deriver

  # 2. Options -----------------------------------------------------------------
  opt <- make_options_from_config(
    option_list = list(
      make_option(
        "--redo_cache",
        type = "logical",
        default = FALSE,
        help = "Re-derive every artefact from raw outputs regardless of freshness [default: %default]"
      ),
      make_option(
        "--workers",
        type = "integer",
        default = NULL,
        help = "Forked workers for the per-group fit loads and per-run hindcast scoring; overrides results.workers. Artefacts are identical for any value."
      ),
      make_option(
        "--steps",
        type = "character",
        default = NULL,
        help = "Comma-separated step ids to run (see --list), e.g. oni,enso_relation:occurrence [default: every step]"
      ),
      make_option(
        "--task",
        type = "character",
        default = NULL,
        help = "A task YAML written by 09a_postprocess_results_tasks.R: run its one step (SLURM array payload)"
      ),
      make_option(
        "--list",
        action = "store_true",
        default = FALSE,
        help = "Print the steps (id, wave, kinds, artefact status) and exit"
      )
    ),
    model = NULL
  )

  if (!is.null(opt[["task"]]) && !is.null(opt[["steps"]])) {
    stop("[09b] pass --steps or --task, not both.", call. = FALSE)
  }
  step_ids <- if (!is.null(opt[["steps"]])) {
    strsplit(opt[["steps"]], ",", fixed = TRUE)[[1]]
  }
  if (!is.null(opt[["task"]])) {
    if (!file.exists(opt[["task"]])) {
      stop("[09b] task file not found: ", opt[["task"]], call. = FALSE)
    }
    task <- yaml::read_yaml(opt[["task"]])
    cfg_path <- results_config_path()
    if (
      !identical(
        normalizePath(task$config, mustWork = FALSE),
        normalizePath(cfg_path, mustWork = FALSE)
      )
    ) {
      stop(
        "[09b] task ",
        opt[["task"]],
        " was written for config ",
        task$config,
        ", not ",
        cfg_path,
        ". Re-run 09a_postprocess_results_tasks.R.",
        call. = FALSE
      )
    }
    step_ids <- task$step
    # The task records 09a's --redo_cache; the 09b flag can only add to it.
    if (isTRUE(task$redo_cache)) opt[["redo_cache"]] <- TRUE
  }

  # 3. Derive the results ------------------------------------------------------
  # CACHING CONTRACT: "cache" = DERIVED aggregations only. This script NEVER
  # re-fits a model and NEVER re-runs the bpfilter / smoothing. It reads raw
  # outputs already on disk (brms fits + cb_meta sidecars, the 08b smoothing
  # part.qs partitions, the 07b forecast/counterfactual monthly_summary.qs
  # partitions, ENSO caches) and writes only the small derived tables to
  # generated_data/results/<tag>/. Caching goes through results_core.R's
  # make_results_run_file()/results_read_if_fresh()/results_write(); each
  # derive_*() returns its table rather than writing a cache path itself. Each
  # artefact is stamped with the decisions (results_write()) that produced it
  # and is invalidated by either a stale input (mtime) OR a changed decision
  # (the stamp diff) — see results_core.R's freshness-contract section header.
  # The derivers live in analysis/utils/results_derive.R; with no step
  # selection results_dispatch() runs every step, --steps / --task run a
  # subset of .RESULTS_STEP_TABLE (results_run_steps()).
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
  # Exit here, explicitly. --workers > 1 forks (base.R::map_workers() ->
  # parallel::mclapply); a child the parallel package fails to terminate shares
  # this script's open file offset, so on falling off the end Rscript read on
  # from a moved offset and died with `unexpected '[' in "["` -- exit 1 after
  # every artefact was written, which fails a SLURM afterok chain.
  quit(save = "no", status = 0)
}
