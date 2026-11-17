#!/usr/bin/env Rscript
# 09a_postprocess_results_tasks.R — one task YAML per manuscript-results step,
#   for the SLURM array driver
#
# Usage:
#   Rscript analysis/09a_postprocess_results_tasks.R \
#     -C analysis/configs/results/manuscript.yml
#   Rscript analysis/09a_postprocess_results_tasks.R -C <cfg> \
#     --steps hindcast,oni
#   Rscript analysis/09a_postprocess_results_tasks.R -C <cfg> --redo_cache TRUE
#   Rscript analysis/09a_postprocess_results_tasks.R -C <cfg> --list
#
# Details:
#   Steps come from .RESULTS_STEP_TABLE (utils/results_derive.R), each with a
#   dependency wave. Derives nothing and loads no fit; builds the results
#   context so missing inputs fail on the login node. Writes
#   <outdir>/results/<config stem>/tasks/<submission>/task_NNNN.yml (a fresh
#   sub-<timestamp>-<pid> folder per call) and prints TASKS_DIR, N_TASKS,
#   N_WAVES, WAVE_<k> for analysis/bash/stages/09_postprocess_results.sh.
#   --redo_cache is recorded in each task; --list writes nothing.

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
    library(dlnm)
  })

  # Verified closure (grepped against every function call in the derivation
  # half): results_core carries the track/group registries, the freshness +
  # provenance contract, and the manifest/config accessors. gid/sim_io/
  # simulator_core/figures_cb/loo_lags are each called by the derivers in utils/results_derive.R
  # (build_location_frame, drop_overlapping_adm1_from_frame, build_coverage_mask,
  # flatten_gid, build_onebasis_grid/erf_from_draws/get_re_vec/get_re_draws_mat,
  # select_best_lag). The Phase-11 "key numbers" report adds population.R
  # (load_gid_population/resolve_gid_population), dur_post.R (duration_ratio_draws/
  # duration_ratio_summary), and dlnm.R (build_cb_columns_from_Q, for the
  # construction-faithful duration RMD newdata). plots.R adds get_covar_type(), the
  # classifier results_core's param_family() wraps for derive_track_params();
  # plots.R's own closure (clim_compile for get_precip_vars/get_temp_vars, clim_db)
  # comes with it. Deliberately NOT spatpomp_core.R: nothing here may be able to
  # invoke a particle filter. deploy_core.R is sourced so that a non-empty
  # results.deploy_config resolves (make_options_from_config() hard-stops if the
  # pointer is set but resolve_results_inputs_from_deploy() is absent — the
  # "dormant branch / missing source" trap); it is inert when deploy_config is "".
  # obs_sim_compare / sim_metrics / hindcast_skill: the hindcast-skill section
  # (ensure_hindcast(); scoring in utils/hindcast_skill.R). sim_metrics is pure
  # at load; 09 calls only load_observed_from_panel() from it.
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
    "hindcast_skill"
  )) {
    source(file.path("analysis/utils", paste0(.m, ".R")))
  }
  source("analysis/utils/results_derive.R") # results_context() / results_dispatch() + every deriver

  # 2. Options and results context ---------------------------------------------
  opt <- make_options_from_config(
    option_list = list(
      make_option(
        "--redo_cache",
        type = "logical",
        default = FALSE,
        help = "Record redo_cache in every task: 09b re-derives regardless of freshness [default: %default]"
      ),
      make_option(
        "--steps",
        type = "character",
        default = NULL,
        help = "Comma-separated step ids to materialise [default: every step]"
      ),
      make_option(
        "--list",
        action = "store_true",
        default = FALSE,
        help = "Print the steps (id, wave, kinds, artefact status) and exit without writing"
      )
    ),
    model = NULL
  )

  ctx <- results_context(opt)
  step_ids <- if (!is.null(opt[["steps"]])) {
    strsplit(opt[["steps"]], ",", fixed = TRUE)[[1]]
  }
  if (isTRUE(opt[["list"]])) {
    tbl <- results_steps_table(ctx, step_ids)
    options(width = 200)
    message(paste(
      capture.output(print(tbl, row.names = FALSE)),
      collapse = "\n"
    ))
    quit(status = 0)
  }

  # 3. Select steps and waves --------------------------------------------------
  # The steps are the registry in analysis/utils/results_derive.R
  # (.RESULTS_STEP_TABLE expanded by results_steps(): track-keyed steps once per
  # results.scope.fig1_tracks entry, the hindcast step only when
  # results.scope.hindcast.enable). Each step gets a dependency WAVE: steps of
  # one wave never read each other's artefacts, so the driver submits one array
  # per wave and chains wave k + 1 afterok on wave k.
  # 09b_postprocess_results.R --task <task.yml> runs one step.
  steps <- results_select_steps(ctx, step_ids)
  if (!length(steps)) {
    stop("[09a] no steps selected.", call. = FALSE)
  }
  # A subset keeps each step's own wave; renumber the waves present to 1..K so
  # the driver never submits an empty array.
  waves <- vapply(steps, `[[`, integer(1), "wave")
  waves <- match(waves, sort(unique(waves)))

  # 4. Write task YAMLs --------------------------------------------------------
  # One fresh <submission> folder (sub-<timestamp>-<pid>) per call, so a later
  # submission never overwrites the task files of an earlier one whose waves are
  # still pending (a shared folder would make a pending wave-2 job fail with
  # "no task YAML"). NNNN - 1 is the SLURM array index.
  tasks_dir <- file.path(
    make_results_run_file(
      ctx$opt,
      "results_tasks_dir",
      config_path = ctx$CONFIG_PATH
    ),
    sprintf("sub-%s-%d", format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid())
  )
  if (dir.exists(tasks_dir)) {
    stop("[09a] submission folder already exists: ", tasks_dir, call. = FALSE)
  }
  dir.create(tasks_dir, recursive = TRUE)

  for (i in seq_along(steps)) {
    s <- steps[[i]]
    yaml::write_yaml(
      list(
        task_id = sprintf("task_%04d", i),
        step = s$id,
        wave = waves[[i]],
        config = ctx$CONFIG_PATH,
        kinds = s$kinds,
        after = s$after,
        artefacts = results_step_paths(ctx, s),
        redo_cache = isTRUE(ctx$REDO)
      ),
      file.path(tasks_dir, sprintf("task_%04d.yml", i))
    )
  }
  message(sprintf(
    "[09a] wrote %d task YAML(s) in %d wave(s) to %s",
    length(steps),
    max(waves),
    tasks_dir
  ))

  # 5. Print for bash capture --------------------------------------------------
  # Stdout, for bash capture:
  #   TASKS_DIR=<path>
  #   N_TASKS=<N>
  #   N_WAVES=<K>
  #   WAVE_<k>=<comma-separated 0-based array indices>    (one line per wave)
  cat(sprintf("TASKS_DIR=%s\n", tasks_dir))
  cat(sprintf("N_TASKS=%d\n", length(steps)))
  cat(sprintf("N_WAVES=%d\n", max(waves)))
  for (k in seq_len(max(waves))) {
    cat(sprintf(
      "WAVE_%d=%s\n",
      k,
      paste(which(waves == k) - 1L, collapse = ",")
    ))
  }
}
