# results_derive.R — manuscript-results derivers (phase 09)
#
# Provides:
#   results_context()    resolve the results config into the read-only `ctx`
#   results_dispatch()   run the selected steps of the step registry
#   results_steps(), results_step_waves(), results_select_steps()
#                        the step registry (.RESULTS_STEP_TABLE)
#   ensure_artefact()    derive-or-reuse one artefact under the stamp
#   derive_*(), ensure_*()   one deriver per artefact; each takes `ctx` first
#
# Details:
#   09b_postprocess_results.R is a guarded CLI: results_context(opt) then
#   results_dispatch(ctx); 09a writes one task per step and
#   stages/09_postprocess_results.sh runs them as wave-chained arrays.
#   Sources nothing: callers source the closure 09 sources (results_core.R
#   and the modules 09 lists). Deliberately NOT spatpomp_core.R: nothing here
#   may invoke a particle filter.
#   Gate: tests/test_results_artefacts_equivalence.R (every 09 artefact
#   identical to a pre-refactor baseline).

# ---------------------------------------------------------------------------
# Context and dispatch
# ---------------------------------------------------------------------------

# The closure 09 sources for this module: results_core.R (which sources its
# own closure) plus gid, sim_io, simulator_core, figures_cb, loo_lags,
# population, dur_post, dlnm, plots, clim_compile, clim_db, deploy_core; for
# the hindcast section obs_sim_compare, sim_metrics, hindcast_skill; for the
# supplementary tables anchor (anchor_enabled(), anchor_spec()) and diagnostics
# (fit_diagnostics()).
#' results_context
#'
#' Resolve the results config into the read-only context every deriver reads:
#' the `results:` blocks, manifest, sets, thresholds and every input path, in
#' the order 09 always resolved them (the manifest is loaded and the prep
#' artefacts' episode definitions are checked here, so a missing input fails
#' at the same point with the same message).
#'
#' @param opt options from make_options_from_config() (with `results:`)
#' @return environment holding the names listed in `.RESULTS_CTX_VARS`
results_context <- function(opt) {
  if (is.null(opt$results)) {
    stop(
      "[09] Config must contain a `results:` block. Got: ",
      paste(names(opt), collapse = ", ")
    )
  }

  CONFIG_PATH <- results_config_path()
  RES <- opt$results
  REDO <- isTRUE(opt$redo_cache)
  # Forked workers for the per-group fit loads and the per-run hindcast scoring
  # (map_workers(), base.R). CLI --workers > results.workers > 1. Speed only:
  # every artefact is identical for any value, so it is not stamped.
  WORKERS <- as.integer(opt$workers %||% RES$workers %||% 1L)
  if (length(WORKERS) != 1L || is.na(WORKERS) || WORKERS < 1L) {
    stop("[09] results.workers / --workers must be a single integer >= 1.")
  }

  INPUTS <- RES$inputs
  SCOPE <- RES$scope
  DERIVATION <- RES$derivation
  THRESHOLDS <- RES$thresholds

  MANIFEST_PATH <- INPUTS$manifest
  GROUPING_PATH <- INPUTS$grouping
  SIM_ROOT <- INPUTS$sim_root
  SMOOTHING_PARTS_DIR <- INPUTS$smoothing_parts_dir
  LAG_LOO_BATCH <- INPUTS$lag_loo_batch
  # Run-specific inputs have no default ("" = unset): the results YAML must set
  # them, or derive them from results.deploy_config.
  if (is.null(SMOOTHING_PARTS_DIR) || !nzchar(SMOOTHING_PARTS_DIR)) {
    stop(
      "[09] results.inputs.smoothing_parts_dir is unset. Set it in the results ",
      "YAML or point results.deploy_config at a deploy/prod_run<N>.yml.",
      call. = FALSE
    )
  }
  BEST_LAG_DIR <- INPUTS$best_lag_dir

  MAIN_SET <- as.character(SCOPE$main_set)
  EXCLUDE_GROUPS <- setdiff(
    as.character(SCOPE$exclude_groups %||% character(0)),
    ""
  )
  FIG1_TRACKS <- as.character(SCOPE$fig1_tracks)

  SMOOTH_START <- as.Date(DERIVATION$smoothing_start)
  N_IRR_DRAWS <- as.integer(DERIVATION$n_irr_draws)
  ENFORCE_OCC_FLOOR <- isTRUE(DERIVATION$enforce_occ_floor)
  # Forecast arms launching at different months (results_core.R::load_arm_run()).
  ARM_ALIGNMENT <- as.character(DERIVATION$arm_alignment %||% "stitch")
  if (
    length(ARM_ALIGNMENT) != 1L ||
      !ARM_ALIGNMENT %in% c("stitch", "common_start")
  ) {
    stop(
      "[09] results.derivation.arm_alignment must be 'stitch' or 'common_start'; got '",
      paste(ARM_ALIGNMENT, collapse = ","),
      "'.",
      call. = FALSE
    )
  }

  PROB_THRESHOLD <- THRESHOLDS$prob_threshold
  LOW_SIGNAL <- THRESHOLDS$low_signal

  # --- Phase-11 "key numbers" report decisions -------------------------------
  POP_REF_YEAR <- as.integer(DERIVATION$pop_ref_year)
  # Per-capita rate population: "ref" (POP_REF_YEAR) or
  # "window" (month-weighted mean over the months read; rate_pop_weights()).
  RATE_POP_YEAR <- as.character(DERIVATION[["rate_pop_year"]] %||% "ref")
  if (!RATE_POP_YEAR %in% c("ref", "window")) {
    stop(
      "[09] results.derivation.rate_pop_year must be 'ref' or 'window'; got '",
      RATE_POP_YEAR,
      "'.",
      call. = FALSE
    )
  }
  DURATION_RATIO_HORIZON <- as.integer(DERIVATION$duration_ratio_horizon)
  POP_PROB_CUTOFFS <- as.numeric(unlist(THRESHOLDS$pop_prob_cutoffs))
  POP_MAGNITUDES <- as.numeric(unlist(THRESHOLDS$pop_magnitudes))
  # ONI conditions every composite IRR / population table is evaluated at
  # (named via figures_cb.R::enso_conds_named()), and the ENSO-sensitive
  # population threshold (|log IRR| > log(1 + SENSITIVE_IRR) on any track).
  ENSO_CONDS <- enso_conds_named(DERIVATION$enso_conds)
  SENSITIVE_IRR <- as.numeric(THRESHOLDS$sensitive_irr)
  # Magnitude axis of the ENSO exposure tables (IRR beyond 1 + m; 0 = any effect).
  ENSO_MAGNITUDES <- sort(unique(as.numeric(unlist(
    THRESHOLDS$enso_magnitudes
  ))))
  if (
    !length(ENSO_MAGNITUDES) ||
      anyNA(ENSO_MAGNITUDES) ||
      any(ENSO_MAGNITUDES < 0)
  ) {
    stop(
      "[09] results.thresholds.enso_magnitudes must be non-negative numbers; got ",
      paste(THRESHOLDS$enso_magnitudes, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    length(SENSITIVE_IRR) != 1L ||
      !is.finite(SENSITIVE_IRR) ||
      SENSITIVE_IRR <= 0
  ) {
    stop(
      "[09] results.thresholds.sensitive_irr must be one positive number; got ",
      paste(THRESHOLDS$sensitive_irr, collapse = ", "),
      call. = FALSE
    )
  }
  # Forecast vs counterfactual exposure (report_attrib_exposure): magnitude axis
  # (mean forecast >= (1 + m) x mean counterfactual) and the IPCC robustness
  # shares (>= ROBUST_EXCEED of the pairs beyond (1 + m), >= ROBUST_AGREE f > cf).
  ATTRIB_MAGNITUDES <- sort(unique(as.numeric(unlist(
    THRESHOLDS$attrib_magnitudes
  ))))
  if (
    !length(ATTRIB_MAGNITUDES) ||
      anyNA(ATTRIB_MAGNITUDES) ||
      any(ATTRIB_MAGNITUDES < 0)
  ) {
    stop(
      "[09] results.thresholds.attrib_magnitudes must be non-negative numbers; got ",
      paste(THRESHOLDS$attrib_magnitudes, collapse = ", "),
      call. = FALSE
    )
  }
  ROBUST_EXCEED <- as.numeric(THRESHOLDS$robust_exceed)
  ROBUST_AGREE <- as.numeric(THRESHOLDS$robust_agree)
  for (k in c("robust_exceed", "robust_agree")) {
    v <- as.numeric(THRESHOLDS[[k]])
    if (length(v) != 1L || !is.finite(v) || v <= 0 || v > 1) {
      stop(
        "[09] results.thresholds.",
        k,
        " must be one number in (0, 1]; got ",
        paste(THRESHOLDS[[k]], collapse = ", "),
        call. = FALSE
      )
    }
  }
  # report_attrib_metrics: the cases_ge event threshold (the database's
  # outbreak-size convention), the one-sided exact level of the PR lower bound
  # when the counterfactual never has the event, and the S/N emergence levels.
  PR_CASE_THRESHOLD <- as.numeric(THRESHOLDS$pr_case_threshold)
  if (
    length(PR_CASE_THRESHOLD) != 1L ||
      !is.finite(PR_CASE_THRESHOLD) ||
      PR_CASE_THRESHOLD <= 0
  ) {
    stop(
      "[09] results.thresholds.pr_case_threshold must be one positive number; got ",
      paste(THRESHOLDS$pr_case_threshold, collapse = ", "),
      call. = FALSE
    )
  }
  PR_ALPHA <- as.numeric(THRESHOLDS$pr_alpha)
  if (
    length(PR_ALPHA) != 1L ||
      !is.finite(PR_ALPHA) ||
      PR_ALPHA <= 0 ||
      PR_ALPHA >= 1
  ) {
    stop(
      "[09] results.thresholds.pr_alpha must be one number in (0, 1); got ",
      paste(THRESHOLDS$pr_alpha, collapse = ", "),
      call. = FALSE
    )
  }
  SN_LEVELS <- sort(unique(as.numeric(unlist(THRESHOLDS$sn_levels))))
  if (!length(SN_LEVELS) || anyNA(SN_LEVELS) || any(SN_LEVELS <= 0)) {
    stop(
      "[09] results.thresholds.sn_levels must be positive numbers; got ",
      paste(THRESHOLDS$sn_levels, collapse = ", "),
      call. = FALSE
    )
  }

  # Explosive-replicate guard: the cap of the "drop" mode and
  # the modes every attribution cell is derived under (drop | none).
  EXPLOSIVE_CAP <- as.numeric(unlist(THRESHOLDS[["explosive_case_cap"]]))
  if (
    length(EXPLOSIVE_CAP) != 1L ||
      !is.finite(EXPLOSIVE_CAP) ||
      EXPLOSIVE_CAP <= 0
  ) {
    stop(
      "[09] results.thresholds.explosive_case_cap must be one positive number; got ",
      paste(THRESHOLDS[["explosive_case_cap"]], collapse = ", "),
      call. = FALSE
    )
  }
  # Per-gid monthly cap of the "drop" mode (Inf = none).
  EXPLOSIVE_GID_CAP <- as.numeric(unlist(THRESHOLDS[[
    "explosive_gid_month_cap"
  ]]))
  if (
    length(EXPLOSIVE_GID_CAP) != 1L ||
      is.na(EXPLOSIVE_GID_CAP) ||
      EXPLOSIVE_GID_CAP <= 0
  ) {
    stop(
      "[09] results.thresholds.explosive_gid_month_cap must be one positive number (.inf = none); got ",
      paste(THRESHOLDS[["explosive_gid_month_cap"]], collapse = ", "),
      call. = FALSE
    )
  }
  EXPLOSIVE_MODES <- as.character(unlist(DERIVATION[["explosive_modes"]]))
  if (
    !length(EXPLOSIVE_MODES) ||
      !all(EXPLOSIVE_MODES %in% c("drop", "none")) ||
      anyDuplicated(EXPLOSIVE_MODES)
  ) {
    stop(
      "[09] results.derivation.explosive_modes must be distinct values of {drop, none}; got ",
      paste(EXPLOSIVE_MODES, collapse = ", "),
      call. = FALSE
    )
  }

  # ENSO phase rule of report_smoothing's by-phase table:
  # "threshold" (legacy) or "episode" (the NOAA CPC runs of the oni artefact).
  ENSO_PHASE_SMOOTHING <- as.character(
    DERIVATION[["enso_phase"]][["smoothing"]] %||% "threshold"
  )
  if (
    length(ENSO_PHASE_SMOOTHING) != 1L ||
      !ENSO_PHASE_SMOOTHING %in% c("threshold", "episode")
  ) {
    stop(
      "[09] results.derivation.enso_phase.smoothing must be 'threshold' or 'episode'; got '",
      paste(ENSO_PHASE_SMOOTHING, collapse = ","),
      "'.",
      call. = FALSE
    )
  }

  SETS <- resolve_results_sets(RES)

  MAN <- load_manifest(MANIFEST_PATH)
  # Every excluded name must be a group of the registry.
  .resolve_exclude_groups(EXCLUDE_GROUPS, MAN, GROUPING_PATH)
  # Library-mode counterfactual cache: the counterfactual runs
  # drew their ENSO from the cache the manifest records (no-op for a legacy
  # manifest, which has no enso_counterfactual).
  for (.aset in SETS) {
    assert_counterfactual_cache_identity(
      MAN[["enso_counterfactual"]],
      manifest_runs(MAN, "counterfactual", .aset)
    )
  }
  # Group short names of the manifest's fit registry (enso_grouping_prod_run has
  # enso_nga); equals results_core's GROUP_SHORT under enso_grouping.
  GROUP_SHORT_RUN <- results_group_meta(list(MAN = MAN))$short
  # The smoothing arm (Figure 3) and the stitched factual prefix are read from
  # SMOOTHING_PARTS_DIR, which must be the parts/ of the cache every forward run
  # handed off from -- else two different filter runs are glued together.
  .si <- unique(unlist(lapply(manifest_runs(MAN, "factual"), function(r) {
    yaml::read_yaml(r$config)$simulator$filter$smoothing_input
  })))
  .want <- unique(vapply(.si, smoothing_parts_dir_of, character(1)))
  if (length(.want) && !identical(.want, sub("/+$", "", SMOOTHING_PARTS_DIR))) {
    stop(
      "[09] results.inputs.smoothing_parts_dir is '",
      SMOOTHING_PARTS_DIR,
      "' but the ",
      "manifest's factual runs hand off from ",
      paste(.si, collapse = ", "),
      " (parts: ",
      paste(.want, collapse = ", "),
      "). Align it (deploy_prod_run.R --stage align_results).",
      call. = FALSE
    )
  }

  # Attribution factual (results.inputs.attrib_manifest; deploy_prod_run.R
  # --stage gen_attrib, utils/make_forecast_attrib_configs.R): every attribution
  # step reads its factual arm from MAN_ATTRIB (results_attrib_ctx()), launched at
  # the counterfactual handoff, while MAN's factual stays the displayed forecast.
  # "" = no attribution factual: MAN_ATTRIB is MAN.
  ATTRIB_MANIFEST_PATH <- as.character(INPUTS$attrib_manifest %||% "")
  MAN_DISPLAY <- MAN
  MANIFEST_PATH_DISPLAY <- MANIFEST_PATH
  MAN_ATTRIB <- MAN
  ATTRIB_MANIFEST_FILE <- MANIFEST_PATH
  if (nzchar(ATTRIB_MANIFEST_PATH)) {
    MAN_ATTRIB <- load_manifest(ATTRIB_MANIFEST_PATH)
    ATTRIB_MANIFEST_FILE <- ATTRIB_MANIFEST_PATH
    results_attrib_manifest_check(
      MAN,
      MAN_ATTRIB,
      MANIFEST_PATH,
      ATTRIB_MANIFEST_PATH
    )
    .si_a <- unique(unlist(lapply(
      manifest_runs(MAN_ATTRIB, "factual"),
      function(r) {
        yaml::read_yaml(r$config)$simulator$filter$smoothing_input
      }
    )))
    if (!setequal(.si_a, .si)) {
      stop(
        "[09] the attribution factual hands off from ",
        paste(.si_a, collapse = ", "),
        ", the forecast factual from ",
        paste(.si, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  # "attrib" = the attribution steps read another factual than the display;
  # "forecast" = both read the forecast manifest's factual.
  ATTRIB_FACTUAL <- if (results_attrib_factual_differs(MAN, MAN_ATTRIB)) {
    "attrib"
  } else {
    "forecast"
  }
  # The one definition of Africa of the attribution steps (results.scope.africa;
  # results_africa_window()) and the analog continuation's skill
  # by lead behind table S28 ("" = no CRPSS column).
  AFRICA_WIN <- results_africa_window(
    MAN,
    manifest_groups_inc(MAN, EXCLUDE_GROUPS),
    SCOPE[["africa"]]
  )
  CONT_SKILL_PATH <- as.character(INPUTS[["enso_continuation_skill"]] %||% "")

  # Shared raw-output paths. Episode-key aware (utils/paths.R): the historical
  # literals under the default `episodes.key: lp_season`, the "_ep-gs" artefacts
  # under gid_season. The baseline stays pinned at 1990_2019 as before.
  .ep_opt <- list(baseline = "1990_2019", episodes = opt$episodes)
  # The simulator panel the smoothing scope's location frame is built from: the
  # run's own episode key. For the deployed gid_season run (one
  # episodes block) the gid_season panel gives the same location frame as
  # the lp_season one (599 rows, every column equal); artefacts derived from the
  # lp_season panel name it as an input and re-derive once on a --redo_cache run.
  PANEL_PATH <- as.character(make_sim_panel_file(.ep_opt))
  SURV_PATH <- as.character(make_surv_dat_file(.ep_opt))
  SIZE_DATA_PATH <- as.character(make_size_dat_file(.ep_opt))
  MON_DATA_PATH <- as.character(make_monthly_dat_file(.ep_opt))
  DUR_DATA_PATH <- as.character(make_duration_dat_file(.ep_opt))
  # The four prep artefacts must share one episode definition whenever any of
  # them -- or this config's episodes.key -- is gid_season.
  sim_episode_check(
    preps = list(
      surv = SURV_PATH,
      size_data = SIZE_DATA_PATH,
      mon_data = MON_DATA_PATH,
      dur_data = DUR_DATA_PATH
    ),
    what = "09b_postprocess_results inputs",
    episode_key = episode_key_of(opt)
  )
  # Canonical per-(gid, year) population Q_g table (02c_preprocess_population_compile.R
  # output). The manuscript "key numbers" report joins every population number at
  # POP_REF_YEAR through this ONE table via resolve_gid_population(), so the
  # population-covered / significant-effect-population / threshold-population
  # tables cannot drift on which population source they used.
  GID_POP_PATH <- make_pop_run_file(opt, "pop_gid")
  # WPP national totals: the T2 rescale of gid populations outside the gid
  # table's years, read only by rate_pop_year = "window".
  WPP_PATH <- make_pop_run_file(opt, "pop_wpp")
  HMM_PATH <- "data/hmm_country_clusters.csv"
  # GADM ADM_1 layer (figures_cb.R::load_africa_sf()'s default, the ADM1
  # polygons fig. MU draws): the modelled_units step's ADM1 universe.
  GADM_ADM1_PATH <- "data/geodata/gadm_410-levels.gpkg"
  AGG_TS_PATH <- make_agg_ts_file()
  AGG_ADM_TS_PATH <- make_agg_adm_ts_file()
  OB_META_PATH <- as.character(make_outbreak_metadata_file()) # 01's raw outbreak set (report_dataset)
  ENSO_CACHE_DIR <- "generated_data/"
  LENSSEN_CACHE <- "generated_data/lenssen_nino34_analogs.rds"
  # Lag averaging over the single-lag fits (lag_averaging step).
  LAG_AVG <- lag_avg_cfg(DERIVATION$lag_averaging)
  # Figure 1d's rainfall tables (rain_enso step): the ADM1 x month rainfall
  # series of 03_prepare_rain_adm1.R and the method (results.derivation.rain_enso).
  RAIN_PATH <- make_adm1_precip_file(opt$outdir %||% "generated_data/")
  RAIN_DERIV <- DERIVATION$rain_enso
  # Supplement fig. SE's method (seasonality step, results.derivation.seasonality).
  SEAS_DERIV <- DERIVATION[["seasonality"]]
  # Supplement fig. WR's method (report_who_enso_ratio, results.derivation.who_enso_ratio).
  WER_DERIV <- DERIVATION[["who_enso_ratio"]]
  # Confirmed-subset agreement table (report_confirmed_agreement,
  # results.derivation.confirmed_agreement).
  AGREE_DERIV <- DERIVATION[["confirmed_agreement"]]
  # Table T1a's curated literature (results.inputs.enso_rain_mechanisms; "" = none),
  # copied into the rain_enso artefact.
  RAIN_MECH_PATH <- as.character(INPUTS[["enso_rain_mechanisms"]] %||% "")
  # Occurrence-panel sidecars: each step is on only when its sidecar exists
  # ("" = absent), so a panel built without them (older panels) gains no step and no
  # artefact changes. D32: the national-zero audit
  # of 01_preprocess_national_zero_months.R, beside the panel's months file.
  .if_exists <- function(p) {
    if (length(p) && nzchar(p) && file.exists(p)) as.character(p) else ""
  }
  NZ_AUDIT_PATH <- if (national_zero_fill_of(opt)) {
    .if_exists(national_zero_audit_path(national_zero_months_file_of(opt)))
  } else {
    ""
  }
  # D31: the overlap sidecars of the survival panel (03_prepare_occ.R).
  OCC_EXCL_PATH <- .if_exists(occ_panel_sidecar_paths(SURV_PATH)$adm1_overlap)
  # D33: the contagion-coverage sidecar (03 --contagion_coverage_report).
  CONTAGION_COV_PATH <- .if_exists(contagion_coverage_csv_path(SURV_PATH))

  message(sprintf(
    "[09] config=%s sets=%s fig1_tracks=%s redo_cache=%s workers=%d attribution factual=%s",
    CONFIG_PATH,
    paste(SETS, collapse = ","),
    paste(FIG1_TRACKS, collapse = ","),
    REDO,
    WORKERS,
    if (ATTRIB_FACTUAL == "attrib") {
      ATTRIB_MANIFEST_PATH
    } else {
      "forecast manifest"
    }
  ))

  #' The archive of currently-deployed per-track/per-group best-lag fits
  #' (collected by deploy_prod_run.R --stage deploy_fits from the SAME
  #' select_best_lag() rule this script calls directly below). Unlike
  #' forecast_horizons.yml (which only carries each group's deployed fit files
  #' for the spatPomp forecast simulator), this archive has fit AND cb_meta paths
  #' for every deployed modelling track (occurrence/duration/monthly), keyed by
  #' (track, group, kind) — see the file header for the exact kind strings.
  ENSO_BEST_LAG_DIR <- BEST_LAG_DIR
  ENSO_BEST_LAG_MANIFEST <- file.path(ENSO_BEST_LAG_DIR, "MANIFEST.csv")

  # Hindcast skill (utils/hindcast_skill.R). The manifest is read here so a
  # missing manifest fails at context time, like MAN above. The runs must share
  # this config's episodes.key; ensure_hindcast() refuses a run simulated under
  # another, so they are scored against SURV_PATH / MON_DATA_PATH.
  HINDCAST_MANIFEST_PATH <- INPUTS$hindcast_manifest
  HC_SCOPE <- SCOPE$hindcast
  HC_DERIV <- DERIVATION$hindcast
  HC_MAN <- if (isTRUE(HC_SCOPE$enable)) {
    read_hindcast_manifest(HINDCAST_MANIFEST_PATH, SURV_PATH)
  } else {
    NULL
  }

  # Supplementary tables: processing counts of 00/01/02 (cascade, table S1) and
  # the brms_kl class-sweep comparison dir (05_fit_postprocess_kl_sweep.R output;
  # "" = the report_kl_sweep step is disabled).
  PREP_COUNTS_DIR <- "generated_data/prep_counts"
  KL_SWEEP_DIR <- as.character(RES$kl_sweep_dir %||% "")
  # Its deployed cell (report_membership; "" = that step is off).
  KL_SWEEP_CELL <- as.character(RES$kl_sweep_cell %||% "")
  # Past-ENSO-event attribution manifest (deploy --stage gen_events); "" = the
  # report_event_attrib step is off.
  EVENT_MANIFEST_PATH <- as.character(RES$event_manifest %||% "")
  # Posterior interval of the attributable fraction (af_posterior /
  # event_af_posterior steps; fixed-draw runs only) and the paired-bootstrap
  # resample count of report_attrib_exposure / report_attrib_metrics /
  # report_event_attrib (0 = no bootstrap: interval columns NA, points unchanged).
  AF_POSTERIOR <- isTRUE(RES$af_posterior) # top level, unstamped: selects the two steps
  # Replicates each (unit, draw) needs after the explosive drop in the
  # af_posterior steps (results.af_posterior_min_reps; default AF_POST_MIN_REPS).
  # Top level, unstamped: a change re-derives a present artefact only with --redo_cache.
  AF_MIN_REPS <- as.integer(
    RES[["af_posterior_min_reps"]] %||% AF_POST_MIN_REPS
  )
  if (length(AF_MIN_REPS) != 1L || is.na(AF_MIN_REPS) || AF_MIN_REPS < 2L) {
    stop(
      "[09] results.af_posterior_min_reps must be one integer >= 2; got ",
      paste(RES[["af_posterior_min_reps"]], collapse = ", "),
      call. = FALSE
    )
  }
  # af_posterior steps: the explosive guard
  # scope (stamped derivation.af_guard_scope) and the unidentified threshold
  # on the Monte-Carlo noise share (stamped thresholds.af_noise_share_max).
  AF_GUARD_SCOPE <- as.character(DERIVATION[["af_guard_scope"]] %||% "unit")
  if (length(AF_GUARD_SCOPE) != 1L || !AF_GUARD_SCOPE %in% c("unit", "group")) {
    stop(
      "[09] results.derivation.af_guard_scope must be 'unit' or 'group'; got ",
      paste(DERIVATION[["af_guard_scope"]], collapse = ", "),
      call. = FALSE
    )
  }
  AF_NOISE_SHARE_MAX <- as.numeric(THRESHOLDS[["af_noise_share_max"]] %||% 0.8)
  if (
    length(AF_NOISE_SHARE_MAX) != 1L ||
      is.na(AF_NOISE_SHARE_MAX) ||
      AF_NOISE_SHARE_MAX <= 0
  ) {
    stop(
      "[09] results.thresholds.af_noise_share_max must be one positive number (.inf = never flag); got ",
      paste(THRESHOLDS[["af_noise_share_max"]], collapse = ", "),
      call. = FALSE
    )
  }
  BOOT_B <- as.integer(DERIVATION$attrib_boot_b %||% ATTRIB_BOOT_B)
  if (length(BOOT_B) != 1L || is.na(BOOT_B) || BOOT_B < 0L) {
    stop(
      "[09] results.derivation.attrib_boot_b must be one integer >= 0; got ",
      paste(DERIVATION$attrib_boot_b, collapse = ", "),
      call. = FALSE
    )
  }
  # af_posterior inputs kept from this job's attribution loads (per-replicate
  # window totals; af_posterior_stash_put()). An environment, so the event
  # views of .results_ctx_with() share it.
  AF_POST_STASH <- new.env(parent = emptyenv())
  # The last load_group_arm_runs() result of this job, reused by the next step
  # that asks for the same runs (attribution, report_forecast,
  # report_attrib_exposure / _metrics, af_posterior: one load, not four).
  RUN_LOAD_MEMO <- new.env(parent = emptyenv())
  # The deployed fits' specifications the convergence step reads in its fit
  # load, for report_model_spec (table S34) of the same job.
  MODEL_SPEC_STASH <- new.env(parent = emptyenv())

  list2env(mget(.RESULTS_CTX_VARS, envir = environment()), parent = emptyenv())
}

#' Check results.scope.exclude_groups against the fit registry and the manifest:
#' a name that is no group of `grouping_path` is an error
#' naming the registry (a stale exclusion, e.g. enso_central under the kl
#' registry); a registry group the manifest neither generated nor skipped is a
#' no-op, said in a message. Returns `exclude` unchanged.
.resolve_exclude_groups <- function(exclude, man, grouping_path) {
  exclude <- setdiff(as.character(unlist(exclude)), "")
  if (!length(exclude)) {
    return(invisible(exclude))
  }
  reg_groups <- unique(load_enso_grouping(grouping_path)$group)
  unknown <- setdiff(exclude, reg_groups)
  if (length(unknown)) {
    stop(
      sprintf(
        paste0(
          "[09] results.scope.exclude_groups names %s, not a group of the registry %s ",
          "(its groups: %s). Remove the name."
        ),
        paste(unknown, collapse = ", "),
        grouping_path,
        paste(sort(reg_groups), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  man_groups <- c(names(man$groups), as.character(unlist(man$skipped_groups)))
  idle <- setdiff(exclude, man_groups)
  if (length(idle)) {
    message(sprintf(
      "[09] results.scope.exclude_groups: %s is not in the manifest (no-op).",
      paste(idle, collapse = ", ")
    ))
  }
  invisible(exclude)
}

#' Names results_context() exports (the former 09 script-level variables).
.RESULTS_CTX_VARS <- c(
  "opt",
  "CONFIG_PATH",
  "RES",
  "REDO",
  "WORKERS",
  "INPUTS",
  "SCOPE",
  "DERIVATION",
  "THRESHOLDS",
  "MANIFEST_PATH",
  "GROUPING_PATH",
  "SIM_ROOT",
  "SMOOTHING_PARTS_DIR",
  "LAG_LOO_BATCH",
  "BEST_LAG_DIR",
  "MAIN_SET",
  "EXCLUDE_GROUPS",
  "FIG1_TRACKS",
  "SMOOTH_START",
  "N_IRR_DRAWS",
  "ENFORCE_OCC_FLOOR",
  "ARM_ALIGNMENT",
  "PROB_THRESHOLD",
  "LOW_SIGNAL",
  "POP_REF_YEAR",
  "DURATION_RATIO_HORIZON",
  "POP_PROB_CUTOFFS",
  "POP_MAGNITUDES",
  "ENSO_CONDS",
  "SENSITIVE_IRR",
  "ENSO_MAGNITUDES",
  "ATTRIB_MAGNITUDES",
  "ROBUST_EXCEED",
  "ROBUST_AGREE",
  "PR_CASE_THRESHOLD",
  "PR_ALPHA",
  "SN_LEVELS",
  "EXPLOSIVE_CAP",
  "EXPLOSIVE_GID_CAP",
  "EXPLOSIVE_MODES",
  "ENSO_PHASE_SMOOTHING",
  "SETS",
  "MAN",
  "MAN_DISPLAY",
  "MANIFEST_PATH_DISPLAY",
  "MAN_ATTRIB",
  "ATTRIB_MANIFEST_PATH",
  "ATTRIB_MANIFEST_FILE",
  "ATTRIB_FACTUAL",
  "AFRICA_WIN",
  "CONT_SKILL_PATH",
  "PANEL_PATH",
  ".ep_opt",
  "SURV_PATH",
  "RAIN_MECH_PATH",
  "SIZE_DATA_PATH",
  "MON_DATA_PATH",
  "DUR_DATA_PATH",
  "GID_POP_PATH",
  "HMM_PATH",
  "GADM_ADM1_PATH",
  "AGG_TS_PATH",
  "AGG_ADM_TS_PATH",
  "OB_META_PATH",
  "ENSO_CACHE_DIR",
  "LENSSEN_CACHE",
  "ENSO_BEST_LAG_DIR",
  "ENSO_BEST_LAG_MANIFEST",
  "HINDCAST_MANIFEST_PATH",
  "HC_SCOPE",
  "HC_DERIV",
  "HC_MAN",
  "GROUP_SHORT_RUN",
  "PREP_COUNTS_DIR",
  "KL_SWEEP_DIR",
  "EVENT_MANIFEST_PATH",
  "AF_POSTERIOR",
  "AF_MIN_REPS",
  "AF_GUARD_SCOPE",
  "AF_NOISE_SHARE_MAX",
  "BOOT_B",
  "AF_POST_STASH",
  "RUN_LOAD_MEMO",
  "MODEL_SPEC_STASH",
  "RAIN_PATH",
  "RAIN_DERIV",
  "SEAS_DERIV",
  "WER_DERIV",
  "AGREE_DERIV",
  "LAG_AVG",
  "KL_SWEEP_CELL",
  "RATE_POP_YEAR",
  "WPP_PATH",
  "NZ_AUDIT_PATH",
  "OCC_EXCL_PATH",
  "CONTAGION_COV_PATH"
)

#' ensure_artefact
#'
#' The freshness-read / derive / stamp-and-write tail every ensure_*() wrapper
#' shares: build the artefact path for `kind`
#' (keyed by `...`, e.g. track = or sets =), return the cached artefact when
#' results_read_if_fresh() accepts it, else run `derive()`, optionally drop its
#' `.inputs` element, and write it through results_write() with the kind's
#' RESULTS_SCHEMA token. Each wrapper still assembles its own `inputs`.
#'
#' @param kind RESULTS_SCHEMA key; the artefact kind is "results_<kind>"
#' @param inputs results_inputs() of the artefact
#' @param derive zero-argument function returning the artefact
#' @param ... key arguments forwarded to make_results_run_file()
#' @param strip_inputs drop `out$.inputs` before writing
#' @return the artefact, invisibly
ensure_artefact <- function(
  ctx,
  kind,
  inputs,
  derive,
  ...,
  strip_inputs = FALSE
) {
  path <- make_results_run_file(
    ctx$opt,
    paste0("results_", kind),
    ...,
    schema = RESULTS_SCHEMA[[kind]],
    config_path = ctx$CONFIG_PATH
  )
  cached <- results_read_if_fresh(path, inputs, ctx$RES, ctx$REDO)
  if (!is.null(cached)) {
    return(invisible(cached))
  }
  out <- derive()
  if (strip_inputs) {
    out$.inputs <- NULL
  }
  results_write(
    out,
    path,
    ctx$RES,
    ctx$CONFIG_PATH,
    RESULTS_SCHEMA[[kind]],
    inputs = as.character(inputs),
    globs = attr(inputs, "globs") %||% list(),
    kind = kind
  )
  invisible(out)
}

#' results_dispatch
#'
#' Derive (or reuse when fresh) every manuscript artefact, in 09's order: every
#' step of the registry (results_steps()), in table order.
results_dispatch <- function(ctx) {
  results_run_steps(ctx)
  # Hindcast forecast skill (10 figure H, notebooks/hindcast_skill.qmd, 11):
  # the step is absent from the registry when disabled; say so, as before.
  if (!isTRUE(ctx$HC_SCOPE$enable)) {
    message("[09] hindcast: results.scope.hindcast.enable is false; skipped.")
  }
  message("[09] done.")
  invisible(NULL)
}

# ===========================================================================
# Step registry (09a materialises one SLURM task per step; 09b runs them)
# ===========================================================================

#' The phase-09 postprocessing steps, in derivation order (the order the former
#' hard-coded dispatch called them, so a full run is call-for-call identical).
#'
#' One row per step stem. `key = "track"` rows expand to one step per
#' results.scope.fig1_tracks entry (id `<stem>:<track>`); `key = "none"` rows
#' are a single step (id `<stem>`). `after` names the steps whose artefacts this
#' one reads: `enso_relation:{track}` is the same track's step, `enso_relation:*`
#' every track's. `kinds` are the RESULTS_SCHEMA artefact kinds the step writes
#' (each kind owned by exactly one step). `hindcast` is ONE step on purpose: its
#' three tables and report share one scoring of every hindcast run
#' (ensure_hindcast()'s memoised bundle); `enabled` drops it when
#' results.scope.hindcast.enable is false.
.RESULTS_STEP_TABLE <- list(
  list(
    stem = "enso_relation",
    key = "track",
    after = character(),
    kinds = "enso_relation"
  ),
  list(
    stem = "track_params",
    key = "track",
    after = character(),
    kinds = "track_params"
  ),
  list(
    stem = "onset_rate",
    key = "none",
    after = character(),
    kinds = "onset_rate"
  ),
  list(
    stem = "attribution",
    key = "none",
    after = character(),
    kinds = "attribution"
  ),
  list(stem = "oni", key = "none", after = character(), kinds = "oni"),
  list(
    stem = "smoothing",
    key = "none",
    after = character(),
    kinds = "smoothing"
  ),
  list(
    stem = "report_dataset",
    key = "none",
    after = character(),
    kinds = "report_dataset"
  ),
  list(
    stem = "report_effects",
    key = "track",
    after = "enso_relation:{track}",
    kinds = "report_effects"
  ),
  # after "oni": the episode rule of results.derivation.enso_phase.smoothing reads its runs.
  list(
    stem = "report_smoothing",
    key = "none",
    after = "oni",
    kinds = "report_smoothing"
  ),
  list(
    stem = "report_forecast",
    key = "none",
    after = character(),
    kinds = "report_forecast"
  ),
  list(
    stem = "report_enso_exposure",
    key = "none",
    after = "enso_relation:*",
    kinds = "report_enso_exposure"
  ),
  list(
    stem = "hindcast",
    key = "none",
    after = character(),
    kinds = c(
      "hindcast_runs",
      "hindcast_occ",
      "hindcast_cases",
      "report_hindcast"
    ),
    enabled = function(ctx) isTRUE(ctx$HC_SCOPE$enable)
  ),
  # Supplementary tables. convergence reads every deployed fit of a
  # track; report_convergence (S8) and report_software (S13) read it.
  list(
    stem = "convergence",
    key = "track",
    after = character(),
    kinds = c("convergence", "report_model_spec")
  ), # + the fit-level model spec (S34), same fit load
  list(
    stem = "report_convergence",
    key = "none",
    after = "convergence:*",
    kinds = "report_convergence"
  ),
  list(
    stem = "report_lags",
    key = "none",
    after = "enso_relation:*",
    kinds = "report_lags"
  ),
  list(
    stem = "report_anchor",
    key = "none",
    after = character(),
    kinds = "report_anchor"
  ),
  list(
    stem = "report_software",
    key = "none",
    after = "convergence:*",
    kinds = "report_software"
  ),
  list(
    stem = "report_kl_sweep",
    key = "none",
    after = character(),
    kinds = "report_kl_sweep",
    enabled = function(ctx) isTRUE(nzchar(ctx$KL_SWEEP_DIR %||% ""))
  ),
  # Population in ADM1s with more ENSO-attributable cases, forecast vs
  # counterfactual; reads the runs report_forecast reads.
  list(
    stem = "report_attrib_exposure",
    key = "none",
    after = character(),
    kinds = "report_attrib_exposure"
  ),
  # Attributable fraction / number, probability ratio (+ FAR, PS, PNS) and
  # signal-to-noise of the same paired runs.
  list(
    stem = "report_attrib_metrics",
    key = "none",
    after = character(),
    kinds = "report_attrib_metrics"
  ),
  # Attribution of past strong ENSO events: the paired event runs
  # of results.event_manifest; off when it is empty.
  list(
    stem = "report_event_attrib",
    key = "none",
    after = character(),
    kinds = "report_event_attrib",
    enabled = function(ctx) isTRUE(nzchar(ctx$EVENT_MANIFEST_PATH %||% ""))
  ),
  # Its analytic cross-check: crossbasis-only AF on the observed panel (loads
  # the deployed occurrence and monthly fits).
  list(
    stem = "report_event_af_model",
    key = "none",
    after = character(),
    kinds = "report_event_af_model",
    enabled = function(ctx) isTRUE(nzchar(ctx$EVENT_MANIFEST_PATH %||% ""))
  ),
  # Posterior interval of the attributable fraction / number of
  # the forecast runs and of the event runs: fixed-draw runs only, on when
  # results.af_posterior is true. Run in the same 09 job as the
  # attribution steps, they reuse their load (ctx$AF_POST_STASH).
  list(
    stem = "af_posterior",
    key = "none",
    after = character(),
    kinds = "af_posterior",
    enabled = function(ctx) isTRUE(ctx$AF_POSTERIOR)
  ),
  list(
    stem = "event_af_posterior",
    key = "none",
    after = character(),
    kinds = "event_af_posterior",
    enabled = function(ctx) {
      isTRUE(ctx$AF_POSTERIOR) &&
        isTRUE(nzchar(ctx$EVENT_MANIFEST_PATH %||% ""))
    }
  ),
  # Posterior interval of the Childs population of report_attrib_exposure,
  # from the same fixed-draw runs.
  list(
    stem = "attrib_exposure_posterior",
    key = "none",
    after = character(),
    kinds = "attrib_exposure_posterior",
    enabled = function(ctx) isTRUE(ctx$AF_POSTERIOR)
  ),
  # The information component of the forecast attribution:
  # displayed factual vs attribution factual vs counterfactual, per draw. On
  # with an attribution factual (results.inputs.attrib_manifest) and the fixed
  # draw design (results.af_posterior).
  list(
    stem = "report_attrib_information",
    key = "none",
    after = character(),
    kinds = "report_attrib_information",
    enabled = function(ctx) {
      identical(ctx$ATTRIB_FACTUAL, "attrib") && isTRUE(ctx$AF_POSTERIOR)
    }
  ),
  # Figure 1d's rainfall tables: the rainfall anomaly-ONI
  # correlation tiles and the groups' rainy seasons, from
  # 03_prepare_rain_adm1.R's ADM1 series.
  list(
    stem = "rain_enso",
    key = "none",
    after = character(),
    kinds = "rain_enso"
  ),
  # Group display lags (style.group_labels "lag") and the per-country season
  # ends of table S30: the deployed lags read from
  # the fits, members, counts; ids only.
  list(
    stem = "group_labels",
    key = "none",
    after = character(),
    kinds = "group_labels"
  ),
  # Pareto-k, k-robust elpd differences and stacking / pseudo-BMA+ over the
  # single-lag fits, the stacked population ERF:
  # from the 05 lag sidecars; on with results.derivation.lag_averaging.enable.
  list(
    stem = "lag_averaging",
    key = "track",
    after = character(),
    kinds = "lag_averaging",
    enabled = function(ctx) isTRUE(ctx$LAG_AVG$enable)
  ),
  # Membership uncertainty of the deployed classes:
  # responsibilities, p_max, stability, sensitivity, M1 / M2; on when
  # results.kl_sweep_dir and results.kl_sweep_cell are both set.
  list(
    stem = "report_membership",
    key = "none",
    after = character(),
    kinds = "report_membership",
    enabled = function(ctx) {
      isTRUE(
        nzchar(ctx$KL_SWEEP_DIR %||% "") && nzchar(ctx$KL_SWEEP_CELL %||% "")
      )
    }
  ),
  # Per-country WHO term scale of the smoothing (table D30): 08a's
  # who_scale.csv, or a legacy preview of the implied WHO/DB ratios for a fixed-scale
  # smoothing (older runs).
  list(
    stem = "report_who_scale",
    key = "none",
    after = character(),
    kinds = "report_who_scale",
    enabled = function(ctx) isTRUE(nzchar(ctx$SMOOTHING_PARTS_DIR %||% ""))
  ),
  # Counterfactual analog sets of the forecast and event caches (D29), appended.
  list(
    stem = "report_analog_sets",
    key = "none",
    after = character(),
    kinds = "report_analog_sets"
  ),
  # National-zero audit (table D32): the country-year audit 01
  # writes beside the panel's months file; on when that file exists.
  list(
    stem = "report_national_zero_audit",
    key = "none",
    after = character(),
    kinds = "report_national_zero_audit",
    enabled = function(ctx) isTRUE(nzchar(ctx$NZ_AUDIT_PATH %||% ""))
  ),
  # Occurrence-panel exclusions (table D31): the 03 sidecars of
  # the survival panel; on when they exist.
  list(
    stem = "report_occ_panel_exclusions",
    key = "none",
    after = character(),
    kinds = "report_occ_panel_exclusions",
    enabled = function(ctx) isTRUE(nzchar(ctx$OCC_EXCL_PATH %||% ""))
  ),
  # Contagion coverage (table D33): the neighbour-contagion
  # term's surveillance coverage; on when 03 wrote its sidecar.
  list(
    stem = "report_contagion_coverage",
    key = "none",
    after = character(),
    kinds = "report_contagion_coverage",
    enabled = function(ctx) isTRUE(nzchar(ctx$CONTAGION_COV_PATH %||% ""))
  ),
  # Supplement maps: the modelled units by
  # level, the ADM1s without a recorded episode and the fit membership (figs
  # MU, GM, KC), so 10 only reads.
  list(
    stem = "modelled_units",
    key = "none",
    after = "onset_rate",
    kinds = "modelled_units"
  ),
  # Figure 1c's onset-season column: the deployed occurrence fits'
  # seasonal spline per group and calendar month (loads the occurrence fits).
  list(
    stem = "occ_seasonality",
    key = "none",
    after = character(),
    kinds = "occ_seasonality"
  ),
  # Supplement fig. SE: rainfall climatology and rainy season vs
  # the observed onset, cessation and case rates by calendar month, with the
  # fitted occurrence and monthly seasonality (loads the deployed mon fits).
  list(
    stem = "seasonality",
    key = "none",
    after = c("rain_enso", "occ_seasonality"),
    kinds = "seasonality"
  ),
  # Supplement fig. WR: the WHO-observed ratio of the target year
  # to the ENSO-neutral years against the smoothing / ENSO-counterfactual ratio
  # of means over their common months (reads the smoothing and the
  # counterfactual runs of the units holding the configured groups); after
  # group_labels for the deployed occurrence lags. On when the smoothing and
  # the counterfactual share months and can be paired
  # (.who_enso_ratio_enabled(): off for a common_start counterfactual and for
  # runs without handoff_map.csv).
  list(
    stem = "report_who_enso_ratio",
    key = "none",
    after = "group_labels",
    kinds = "report_who_enso_ratio",
    enabled = function(ctx) .who_enso_ratio_enabled(ctx)
  ),
  # Figure 1c's top row: the composite
  # El Nino / La Nina cycle over the 2-year axis, from the oni artefact.
  list(
    stem = "oni_composite",
    key = "none",
    after = "oni",
    kinds = "oni_composite"
  ),
  # Supplementary table of forecast vs counterfactual incidence by Africa, group
  # and country: window totals of the runs report_attrib_exposure
  # reads, per attribution cell.
  list(
    stem = "report_attrib_incidence",
    key = "none",
    after = character(),
    kinds = "report_attrib_incidence"
  ),
  # Supplementary table of episodes with a confirmed case vs all episodes of
  # the analysis set: ENSO contrasts of occurrence, duration and
  # monthly cases, and the within-country characteristics ratios.
  list(
    stem = "report_confirmed_agreement",
    key = "none",
    after = character(),
    kinds = "report_confirmed_agreement"
  )
)

#' The context of the attribution steps: the manifest whose factual arm is the
#' attribution factual (results.inputs.attrib_manifest), with MANIFEST_PATH its
#' file. Unchanged when the attribution reads the forecast manifest's factual.
#' MAN_DISPLAY / MANIFEST_PATH_DISPLAY keep the forecast manifest.
results_attrib_ctx <- function(ctx) {
  if (!identical(ctx$ATTRIB_FACTUAL, "attrib")) {
    return(ctx)
  }
  .results_ctx_with(
    ctx,
    MAN = ctx$MAN_ATTRIB,
    MANIFEST_PATH = ctx$ATTRIB_MANIFEST_FILE
  )
}

#' The manifest files an attribution artefact depends on: its own (the
#' attribution manifest under results_attrib_ctx()) and the forecast manifest,
#' whose factual launch starts the post_fork window. One file when they agree.
results_attrib_manifest_files <- function(ctx) {
  unique(c(ctx$MANIFEST_PATH, ctx$MANIFEST_PATH_DISPLAY %||% ctx$MANIFEST_PATH))
}

#' The displayed factual's launch over `groups` / `sets`: the start of every
#' post_fork window (the months both forecast arms simulate freely). Equals the
#' window's own factual launch when the attribution reads the forecast manifest.
display_factual_start <- function(ctx, groups, sets) {
  manifest_arm_starts(ctx$MAN_DISPLAY %||% ctx$MAN, groups, sets)$factual
}

#' The zero-argument closure that runs one step: the step's ensure_*() call.
.results_step_runner <- function(ctx, stem, track = NULL) {
  force(ctx)
  force(track)
  switch(
    stem,
    enso_relation = function() ensure_enso_relation(ctx, track),
    track_params = function() ensure_track_params(ctx, track),
    onset_rate = function() ensure_onset_rate(ctx),
    attribution = function() ensure_attribution(ctx, ctx$SETS),
    oni = function() ensure_oni(ctx, ctx$SETS),
    smoothing = function() ensure_smoothing(ctx),
    report_dataset = function() ensure_report_dataset(ctx),
    report_effects = function() ensure_report_effects(ctx, track),
    report_smoothing = function() ensure_report_smoothing(ctx),
    report_forecast = function() ensure_report_forecast(ctx, ctx$SETS),
    report_enso_exposure = function() ensure_report_enso_exposure(ctx),
    report_attrib_exposure = function() {
      ensure_report_attrib_exposure(ctx, ctx$SETS)
    },
    report_attrib_incidence = function() {
      ensure_report_attrib_incidence(ctx, ctx$SETS)
    },
    report_attrib_metrics = function() {
      ensure_report_attrib_metrics(ctx, ctx$SETS)
    },
    report_event_attrib = function() ensure_report_event_attrib(ctx),
    report_event_af_model = function() ensure_report_event_af_model(ctx),
    af_posterior = function() ensure_af_posterior(ctx, ctx$SETS),
    report_attrib_information = function() {
      ensure_report_attrib_information(ctx, ctx$SETS)
    },
    event_af_posterior = function() ensure_event_af_posterior(ctx),
    attrib_exposure_posterior = function() {
      ensure_attrib_exposure_posterior(ctx, ctx$SETS)
    },
    hindcast = function() ensure_hindcast(ctx),
    convergence = function() ensure_convergence(ctx, track),
    report_convergence = function() ensure_report_convergence(ctx),
    report_lags = function() ensure_report_lags(ctx),
    report_anchor = function() ensure_report_anchor(ctx),
    report_software = function() ensure_report_software(ctx),
    report_kl_sweep = function() ensure_report_kl_sweep(ctx),
    rain_enso = function() ensure_rain_enso(ctx),
    group_labels = function() ensure_group_labels(ctx),
    lag_averaging = function() ensure_lag_averaging(ctx, track),
    report_membership = function() ensure_report_membership(ctx),
    report_who_scale = function() ensure_report_who_scale(ctx),
    report_analog_sets = function() ensure_report_analog_sets(ctx),
    report_national_zero_audit = function() {
      ensure_report_national_zero_audit(ctx)
    },
    report_occ_panel_exclusions = function() {
      ensure_report_occ_panel_exclusions(ctx)
    },
    report_contagion_coverage = function() {
      ensure_report_contagion_coverage(ctx)
    },
    modelled_units = function() ensure_modelled_units(ctx),
    occ_seasonality = function() ensure_occ_seasonality(ctx),
    seasonality = function() ensure_seasonality(ctx),
    oni_composite = function() ensure_oni_composite(ctx),
    report_who_enso_ratio = function() ensure_report_who_enso_ratio(ctx),
    report_confirmed_agreement = function() {
      ensure_report_confirmed_agreement(ctx)
    },
    stop("[09] no runner for step stem '", stem, "'.", call. = FALSE)
  )
}

#' results_steps
#'
#' Expand .RESULTS_STEP_TABLE against the context: track rows over
#' ctx$FIG1_TRACKS, `after` placeholders resolved, disabled steps dropped.
#'
#' @return list of steps, each list(id, stem, track, after, kinds, wave, run),
#'   in derivation order
results_steps <- function(ctx) {
  tracks <- ctx$FIG1_TRACKS
  steps <- list()
  for (row in .RESULTS_STEP_TABLE) {
    if (is.function(row$enabled) && !isTRUE(row$enabled(ctx))) {
      next
    }
    keys <- if (identical(row$key, "track")) tracks else list(NULL)
    for (trk in keys) {
      after <- unlist(
        lapply(row$after, function(a) {
          if (endsWith(a, ":*")) {
            paste0(sub(":\\*$", "", a), ":", tracks)
          } else if (!is.null(trk)) {
            gsub("{track}", trk, a, fixed = TRUE)
          } else {
            a
          }
        }),
        use.names = FALSE
      )
      id <- if (is.null(trk)) row$stem else paste0(row$stem, ":", trk)
      steps[[id]] <- list(
        id = id,
        stem = row$stem,
        track = trk,
        after = as.character(after),
        kinds = row$kinds,
        run = .results_step_runner(ctx, row$stem, trk)
      )
    }
  }
  waves <- results_step_waves(steps)
  for (id in names(steps)) {
    steps[[id]]$wave <- waves[[id]]
  }
  steps
}

#' results_step_waves
#'
#' Dependency wave of each step: 1 with no `after`, else 1 + the largest wave
#' among its `after` steps. Steps of one wave never read each other's artefacts,
#' so they can run concurrently; wave k + 1 starts once wave k has finished.
#'
#' @param steps named list of steps (each with `after`)
#' @return named integer vector, one wave per step
results_step_waves <- function(steps) {
  ids <- names(steps)
  unknown <- setdiff(unlist(lapply(steps, `[[`, "after")), ids)
  if (length(unknown)) {
    stop(
      "[09] step dependencies name unknown steps: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  wave <- setNames(rep(NA_integer_, length(ids)), ids)
  for (pass in seq_along(ids)) {
    for (id in ids) {
      dep <- steps[[id]]$after
      if (!length(dep)) {
        wave[[id]] <- 1L
        next
      }
      if (all(!is.na(wave[dep]))) wave[[id]] <- 1L + max(wave[dep])
    }
    if (!anyNA(wave)) break
  }
  if (anyNA(wave)) {
    stop(
      "[09] step dependencies are cyclic: ",
      paste(ids[is.na(wave)], collapse = ", "),
      call. = FALSE
    )
  }
  wave
}

#' results_step_paths -- the artefact file(s) a step writes (the paths
#' ensure_artefact() builds for its kinds).
results_step_paths <- function(ctx, step) {
  vapply(
    step$kinds,
    function(k) {
      make_results_run_file(
        ctx$opt,
        paste0("results_", k),
        track = step$track,
        sets = ctx$SETS,
        schema = RESULTS_SCHEMA[[k]],
        config_path = ctx$CONFIG_PATH
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' results_select_steps -- the registry steps named by `ids` (NULL = all), in
#' registry order; an unknown id stops with the valid ones.
results_select_steps <- function(ctx, ids = NULL) {
  steps <- results_steps(ctx)
  if (is.null(ids)) {
    return(steps)
  }
  ids <- unique(trimws(ids[nzchar(trimws(ids))]))
  bad <- setdiff(ids, names(steps))
  if (length(bad)) {
    stop(
      "[09] unknown step(s): ",
      paste(bad, collapse = ", "),
      ". Valid steps: ",
      paste(names(steps), collapse = ", "),
      call. = FALSE
    )
  }
  steps[names(steps) %in% ids]
}

#' results_run_steps
#'
#' Run the selected steps (NULL = every step) in registry order. Upstream steps
#' are NOT added: a step that reads another's artefact (report_effects reads
#' enso_relation) reuses it when fresh and derives it in-process otherwise, as
#' the ensure_*() calls always did. The SLURM driver runs the steps one wave
#' after another so no two jobs derive the same artefact.
results_run_steps <- function(ctx, ids = NULL) {
  steps <- results_select_steps(ctx, ids)
  for (s in steps) {
    message(sprintf("[09] step %s (wave %d)", s$id, s$wave))
    s$run()
  }
  invisible(names(steps))
}

#' results_steps_table
#'
#' One row per step: id, wave, the kinds it writes and whether every artefact
#' file exists ("present"/"missing"; freshness is decided by each ensure_*()
#' against its own inputs at run time).
results_steps_table <- function(ctx, ids = NULL) {
  steps <- results_select_steps(ctx, ids)
  data.frame(
    task = seq_along(steps) - 1L,
    id = names(steps),
    wave = vapply(steps, `[[`, integer(1), "wave"),
    kinds = vapply(
      steps,
      function(s) paste(s$kinds, collapse = ","),
      character(1)
    ),
    status = vapply(
      steps,
      function(s) {
        if (all(file.exists(results_step_paths(ctx, s)))) {
          "present"
        } else {
          "missing"
        }
      },
      character(1)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


# ===========================================================================
# Hindcast forecast skill (utils/hindcast_skill.R does the scoring)
# ===========================================================================

#' hindcast_inputs
#'
#' Freshness inputs of every hindcast artefact: the manifest, each run's
#' config / validation metrics / merged monthly summary, the observed panels,
#' the population table, the grouping and (when set) the deploy config.
hindcast_inputs <- function(ctx, idx, panels) {
  dep <- ctx$RES$deploy_config %||% ""
  results_inputs(
    files = c(
      ctx$HINDCAST_MANIFEST_PATH,
      idx$config,
      idx$metrics_path,
      idx$monthly_summary_path,
      panels$surv_path,
      panels$mon_path,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH,
      if (nzchar(dep)) dep,
      hindcast_skill_input_files(ctx)
    )
  )
}

#' hindcast_skill_input_files -- the extra inputs of the monthly / season case
#' skill (results.derivation.hindcast.case_skill_supports non-empty): the WPP
#' table under case_pop_scale "wpp", the coverage-mask sources under case_mask
#' "coverage". None while the skill is off (existing artefacts stay fresh).
hindcast_skill_input_files <- function(ctx) {
  cfg <- ctx$HC_DERIV
  if (!length(unlist(cfg[["case_skill_supports"]]))) {
    return(character(0))
  }
  c(
    if (identical(cfg[["case_pop_scale"]], "wpp")) ctx$WPP_PATH,
    if (identical(cfg[["case_mask"]], "coverage")) {
      lpm <- occ_panel_lp_members_path(ctx$SURV_PATH) # coverage_lp_bridge()'s sidecar
      c(
        ctx$AGG_TS_PATH,
        if (file.exists(ctx$AGG_ADM_TS_PATH)) ctx$AGG_ADM_TS_PATH,
        as.character(national_zero_months_file_of(ctx$opt)),
        if (file.exists(lpm)) lpm
      )
    }
  )
}

#' results_coverage_mask_raw
#'
#' build_coverage_mask() of this config's survival panel (TL-TR surveillance
#' windows plus national-zero months), as tibble(gid, month); NULL when the
#' aggregated pre-outbreak series are absent. Shared by the smoothing (which
#' further cuts it to the panel's at-risk cells) and the hindcast case skill.
results_coverage_mask_raw <- function(ctx, surv_dat) {
  if (!file.exists(ctx$AGG_TS_PATH)) {
    return(NULL)
  }
  build_coverage_mask(
    coverage_lp_bridge(surv_dat, ctx$SURV_PATH),
    ctx$AGG_TS_PATH,
    if (file.exists(ctx$AGG_ADM_TS_PATH)) ctx$AGG_ADM_TS_PATH else NULL,
    national_zero_path = national_zero_months_file_of(ctx$opt)
  ) %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))
}

#' derive_hindcast_bundle -- score every run (one read of each run).
derive_hindcast_bundle <- function(ctx, idx, panels) {
  cfg <- ctx$HC_DERIV
  skill_on <- length(unlist(cfg[["case_skill_supports"]])) > 0
  surv <- readRDS(panels$surv_path)
  coverage <- if (skill_on && identical(cfg[["case_mask"]], "coverage")) {
    results_coverage_mask_raw(ctx, surv)
  }
  wpp <- if (skill_on && identical(cfg[["case_pop_scale"]], "wpp")) {
    load_wpp_table(.pop_artefact_or_null(ctx$WPP_PATH))
  }
  hindcast_skill_bundle(
    idx,
    obs_panel = hindcast_obs_panel(surv),
    mon_all = build_mon_obs_frame(panels$mon_path),
    gid_pop = load_gid_population(ctx$GID_POP_PATH),
    cfg = cfg,
    tag = "[09 hindcast]",
    workers = ctx$WORKERS,
    coverage = coverage,
    wpp = wpp
  )
}

#' ensure_hindcast
#'
#' The four hindcast artefacts. The 44 runs are scored at most once per call
#' (a lazily cached bundle shared by the three derive closures), and not at all
#' when every artefact is fresh. report_hindcast is a pure function of the other
#' three.
ensure_hindcast <- function(ctx) {
  if (!isTRUE(ctx$HC_SCOPE$enable)) {
    message("[09] hindcast: results.scope.hindcast.enable is false; skipped.")
    return(invisible(NULL))
  }
  # Every run must share this config's episodes.key (results.inputs.episode_key,
  # deploy-derived): hindcast_run_index() stops on any run simulated under
  # another key, so the runs' panels are exactly SURV_PATH / MON_DATA_PATH.
  idx <- hindcast_run_index(
    ctx$HC_MAN,
    episode_key = episode_key_of(ctx$opt),
    tag = "[09 hindcast]",
    national_zero_fill = national_zero_fill_of(ctx$opt),
    national_zero_threshold = national_zero_threshold_of(ctx$opt),
    occ_panel = occ_panel_spec_of(ctx$opt)
  )
  panels <- hindcast_panel_paths(
    idx,
    tag = "[09 hindcast]",
    occ_panel = occ_panel_spec_of(ctx$opt)
  )
  if (
    !identical(
      c(panels$surv_path, panels$mon_path),
      c(ctx$SURV_PATH, ctx$MON_DATA_PATH)
    )
  ) {
    stop(
      "[09 hindcast] the runs' observed panels (",
      panels$surv_path,
      ", ",
      panels$mon_path,
      ") differ from this config's (",
      ctx$SURV_PATH,
      ", ",
      ctx$MON_DATA_PATH,
      ").",
      call. = FALSE
    )
  }
  inputs <- hindcast_inputs(ctx, idx, panels)
  bundle <- local({
    cache <- NULL
    function() {
      if (is.null(cache)) {
        cache <<- derive_hindcast_bundle(ctx, idx, panels)
      }
      cache
    }
  })
  runs <- ensure_artefact(ctx, "hindcast_runs", inputs, function() {
    bundle()$runs
  })
  occ <- ensure_artefact(ctx, "hindcast_occ", inputs, function() bundle()$occ)
  cases <- ensure_artefact(ctx, "hindcast_cases", inputs, function() {
    bundle()$cases
  })
  ensure_artefact(ctx, "report_hindcast", inputs, function() {
    hindcast_report_tables(
      runs,
      occ,
      cases,
      headline_occ_reference = ctx$HC_SCOPE$headline_occ_reference,
      headline_case_support = ctx$HC_SCOPE$headline_case_support,
      headline_case_skill_reference = ctx$HC_SCOPE[[
        "headline_case_skill_reference"
      ]] %||%
        "clim_pooled"
    )
  })
}

# ===========================================================================
# Section 1: ENSO exposure-response derivation (all tracks; feeds Figures 1 & 2)
# ===========================================================================

#' Batch dir holding one track's single-lag LOO ranking artefacts
#' (analysis/05_fit_postprocess_loo_lags.R output) — the same rankings
#' generated_data/enso_best_lag_fits/ was built from via select_best_lag().
#'
#' `lag_loo_batch` is a per-track map with keys `default`, `occurrence`,
#' `duration` and `monthly` (the tracks of `results.scope.fig1_tracks`; the
#' size track is retired downstream and has no key) rather than
#' one shared string, because the tracks are not refit in lockstep: a later
#' prod run can refit occurrence alone (e.g. to drop the leaked
#' ongoing_outbreak_country_lag1/tu_lag1 spillover covariates) while duration
#' and monthly keep resolving against the earlier batch they were fit in. The
#' deploy tool emits the same four keys (`.deploy_derived_results_inputs()`,
#' `align_results_to_deploy()`). A track with no explicit entry (or "", the
#' unset sentinel) falls back to `default`. Still accepts a bare string for
#' back-compat with any config that has not adopted the per-track form.
LAG_LOO_DIR <- function(ctx, track) {
  # "" is the unset sentinel (the run-specific default), so it falls through to
  # `default` exactly like a missing entry.
  unset <- function(b) is.null(b) || !nzchar(b)
  batch <- if (is.list(ctx$LAG_LOO_BATCH)) {
    b <- ctx$LAG_LOO_BATCH[[track]]
    if (unset(b)) ctx$LAG_LOO_BATCH$default else b
  } else {
    ctx$LAG_LOO_BATCH
  }
  if (unset(batch)) {
    stop(sprintf(
      "[09] results.inputs.lag_loo_batch has no entry for track '%s' and no 'default' fallback.",
      track
    ))
  }
  sprintf(
    "analysis/configs/%s_%s_grouped_enso_singlelag",
    batch,
    TRACK_LOO_DIR_STUB[[track]]
  )
}

# load_best_lag_manifest(ENSO_BEST_LAG_DIR) / track_fit_files(): utils/results_core.R (shared with
# the grouped figure scripts; the manifest carries its archive dir).

#' Most-recent loo_compare_lags_<group>_<YYYYMMDD>.rds `$ranking` tibble per
#' group, stacked with an `enso_group` column. No selection happens here —
#' see selected_lags_by_group() below — this just assembles the raw ranking
#' table plotted as Figure 1's best-lag ELPD-diff tiles.
read_lag_ranking <- function(ctx, track, groups) {
  batch <- LAG_LOO_DIR(ctx, track)
  dplyr::bind_rows(lapply(groups, function(g) {
    # utils/loo_lags.R: <outdir>/loo_compare/<batch>/ first, legacy
    # in-config-dir fallback; latest by filename date, then mtime.
    rk <- read_latest_ranking(batch, group = g, opt = ctx$opt, required = FALSE)
    if (is.null(rk)) {
      warning(sprintf(
        "[09] ENSO exposure-response (%s): no loo_compare_lags rds found for group %s",
        track,
        g
      ))
      return(NULL)
    }
    rk %>% dplyr::mutate(enso_group = g)
  }))
}

# Every ranking cache for a track's batch (all dates, all groups) -- the
# freshness-input list for the artefacts derived from read_lag_ranking().
lag_loo_input_files <- function(ctx, track) {
  loo_compare_ranking_files(
    LAG_LOO_DIR(ctx, track),
    opt = ctx$opt,
    latest_only = FALSE
  )$file
}

#' results_lag_selection -- the deploy.lag_selection block of the results
#' config's deploy (results.deploy_config), normalised by lag_selection_cfg();
#' the current longest_within_se rule when there is no deploy config or block.
results_lag_selection <- function(ctx) {
  dep <- ctx$RES$deploy_config %||% ""
  blk <- if (nzchar(dep) && file.exists(dep)) {
    yaml::read_yaml(dep)$deploy$lag_selection
  } else {
    NULL
  }
  lag_selection_cfg(blk)
}

#' Apply select_best_lag() (analysis/utils/loo_lags.R) per group to a stacked
#' ranking table, returning a named (by enso_group) integer vector of picks.
#' `enforce_occ_floor` defaults to FALSE (results.derivation.enforce_occ_floor):
#' the ENSO exposure-response per-track ERF/ADM1-IRR panels are, by default,
#' each track's own unconstrained best-lag pick, not one floored at the
#' occurrence lag.
selected_lags_by_group <- function(
  ctx,
  rk,
  track,
  occ_lag_by_group = NULL,
  enforce_occ_floor = ctx$ENFORCE_OCC_FLOOR
) {
  groups <- unique(rk$enso_group)
  # The deploy's selection rule (deploy.lag_selection, loo_lags.R::lag_selection_cfg()),
  # so the outlined lag is the deployed fit's under either rule.
  ls <- results_lag_selection(ctx)
  sel <- vapply(
    groups,
    function(g) {
      select_best_lag(
        dplyr::filter(rk, enso_group == g),
        track = track,
        occ_lag = occ_lag_by_group[[g]],
        enforce_occ_floor = enforce_occ_floor,
        se_mult = ls$se_mult,
        min_lag = ls$min_lag,
        rule = ls$rule,
        occ_floor_offset = ls$occ_floor_offset
      )$lag
    },
    integer(1)
  )
  stats::setNames(sel, groups)
}

#' lag_rank_with_selection
#'
#' The stacked LOO lag ranking of one track over `groups`, with the lag the
#' deploy's selection rule picks marked `selected` (shared by
#' derive_enso_relation() and derive_report_lags(), so Figure 1's outlined
#' lag and tables S9/S10 cannot disagree). Warns when the rule's pick differs
#' from the deployed best-lag fit's lag.
#'
#' select_best_lag()'s pick is rule-dependent (analysis/utils/loo_lags.R): the
#' longest lag within se_mult SEs of the rank-1 fit (longest_within_se), or the
#' best-ranked lag at or above the occurrence floor (best_elpd) -- not
#' necessarily rank 1 itself. Marking it explicitly lets the ELPD-diff tiles
#' (Figure 1, 10_manuscript_figures.R) outline the fit actually used to build
#' every other ENSO exposure-response panel.
#'
#' @param mnf load_best_lag_manifest() of the deployed archive
#' @return list(lag_rank = ranking + `selected`, selected = named integer
#'   (group -> rule pick), occ_lag = named list or NULL (the occurrence floor
#'   inputs, when the rule uses them))
lag_rank_with_selection <- function(ctx, track, groups, mnf) {
  lag_rank <- read_lag_ranking(ctx, track, groups)
  occ_lag_by_group <- NULL
  # The occurrence lags are needed for the downstream floor: under
  # enforce_occ_floor, and under the best_elpd rule with an occ_floor_offset.
  ls_rule <- results_lag_selection(ctx)
  if (
    !identical(track, "occurrence") &&
      (isTRUE(ctx$ENFORCE_OCC_FLOOR) ||
        (identical(ls_rule$rule, "best_elpd") &&
          !is.null(ls_rule$occ_floor_offset)))
  ) {
    occ_lag_by_group <- as.list(selected_lags_by_group(
      ctx,
      read_lag_ranking(ctx, "occurrence", groups),
      "occurrence"
    ))
  }
  sel_lag_by_group <- selected_lags_by_group(
    ctx,
    lag_rank,
    track,
    occ_lag_by_group,
    ctx$ENFORCE_OCC_FLOOR
  )
  for (g in groups) {
    # No ranking for this group: read_lag_ranking() already warned.
    if (!g %in% names(sel_lag_by_group)) {
      next
    }
    deployed_lag <- track_fit_files(mnf, track, g)$lag
    if (
      !is.na(sel_lag_by_group[[g]]) &&
        !identical(as.integer(sel_lag_by_group[[g]]), as.integer(deployed_lag))
    ) {
      warning(sprintf(
        paste(
          "[09] ENSO exposure-response (%s): group %s: select_best_lag() pick (lag=%s,",
          "enforce_occ_floor=%s) != deployed enso_best_lag_fits lag (%s)"
        ),
        track,
        g,
        sel_lag_by_group[[g]],
        ctx$ENFORCE_OCC_FLOOR,
        deployed_lag
      ))
    }
  }
  lag_rank <- lag_rank %>%
    dplyr::mutate(selected = single_lag == unname(sel_lag_by_group[enso_group]))
  list(
    lag_rank = lag_rank,
    selected = sel_lag_by_group,
    occ_lag = occ_lag_by_group
  )
}

#' 1/IRR inversion so duration's cessation-hazard estimand reads on the same
#' "high = epidemiologically worse" convention as occurrence/size (see
#' build_fig2()'s shared ADM1-map colour scale, 10_manuscript_figures.R).
#'
#' Not called; kept side by side per the no-deprecation rule. The
#' artefacts written by this script hold the RAW cessation-hazard IRR and
#' any orientation flip belongs to the figure layer. Note what the inversion
#' does and does not claim: 1/HR is a monotone reorientation of the hazard
#' ratio, exactly a ratio of expected durations only under a constant baseline
#' hazard (duration.basehaz_duration: none) -- not under the default `factor`
#' baseline, whose hazard spans a ~492x range across mso bands (dur_post.R:27).
invert_duration_irr <- function(df) {
  if ("IRR" %in% names(df)) {
    df$IRR <- 1 / df$IRR
  }
  if (all(c("IRR_lo", "IRR_hi") %in% names(df))) {
    old_lo <- df$IRR_lo
    old_hi <- df$IRR_hi
    df$IRR_lo <- 1 / old_hi
    df$IRR_hi <- 1 / old_lo
  }
  if ("logRR" %in% names(df)) {
    df$logRR <- -df$logRR
  }
  if (all(c("lo", "hi") %in% names(df))) {
    old_lo <- df$lo
    old_hi <- df$hi
    df$lo <- -old_hi
    df$hi <- -old_lo
  }
  df
}

# extract_adm1_irr() / extract_country_erf(): utils/figures_cb.R (the same extraction
# the grouped figure scripts use).

#' Map a version-less ADM1 gid (e.g. "BDI_1", the flattened form the ENSO fits'
#' grouping levels carry) to the canonical gid_population key (e.g. "BDI_1_1",
#' ISO3_ADM1n_suffix). The suffix segment is stripped from every ADM1 row of
#' gid_pop to build the crosswalk, so it works whatever the GADM version token
#' is; ids that already match (or have no ADM1 row) fall through unchanged.
#' Without this, resolve_gid_population() returns all-NA for the effects-track
#' ADM1 population and the significant-effect population collapses to zero.
canonical_adm1_gid <- function(gids, gid_pop) {
  adm1 <- gid_pop[gid_pop$gadm_level == "ADM1", , drop = FALSE]
  key <- sub("_[0-9]+$", "", adm1$gid)
  xwalk <- stats::setNames(adm1$gid, key)
  xwalk <- xwalk[!duplicated(names(xwalk))]
  mapped <- unname(xwalk[as.character(gids)])
  ifelse(is.na(mapped), as.character(gids), mapped)
}

#' Summarise a C(conditions) x D(draws) population-mass matrix across draws:
#' one row per condition with pop_mean, pop_lo (q025), pop_med, pop_hi (q975).
#' The shared posterior-quantile contract of every per-draw population table
#' (dir_pop_ci() and its Africa rollup, the ENSO-sensitive population).
.pop_summary_rows <- function(mat, enso_conds) {
  dplyr::bind_cols(enso_cond_cols(enso_conds), .pop_summary_cols(mat))
}

#' The condition columns of a population row: one row per ONI condition.
enso_cond_cols <- function(conds) {
  tibble::tibble(enso_label = names(conds), enso_val = unname(conds))
}

#' Mean, q025, median, q975 across draws (columns) of a C x D population-mass
#' matrix, one row per condition (row).
.pop_summary_cols <- function(mat) {
  qs <- apply(mat, 1L, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE) # 3 x C
  tibble::tibble(
    pop_mean = rowMeans(mat),
    pop_lo = qs[1L, ],
    pop_med = qs[2L, ],
    pop_hi = qs[3L, ]
  )
}

#' Median + 95% CrI rows of a C x D population-mass matrix, tagged with a
#' direction (the significant-effect population CrI of derive_report_effects()).
.pop_ci_rows <- function(mat, direction, enso_conds) {
  .pop_summary_rows(mat, enso_conds) %>%
    dplyr::transmute(
      enso_label,
      enso_val,
      direction = direction,
      pop_lo,
      pop_med,
      pop_hi
    )
}

#' Population-weighted per-draw mass: sum over ADM1 a of pop[a] * mask[a, c, d].
#' Accumulated sequentially in ADM1 order (not %*%), so the doubles equal the
#' former per-ADM1 loop in dir_pop_ci() exactly.
#' @param mask logical (or 0/1) array [A x C x D]
#' @param pop numeric[A]; NA counts as 0
#' @return numeric matrix [C x D]
pop_mass_by_draw <- function(mask, pop) {
  d <- dim(mask)
  pop[is.na(pop)] <- 0
  acc <- matrix(0, d[2L], d[3L])
  for (a in seq_len(d[1L])) {
    acc <- acc + matrix(mask[a, , ], d[2L], d[3L]) * pop[[a]]
  }
  acc
}

#' Bind per-group extract_adm1_logrr_draws() results into one ADM1 x C x D array
#' (groups stacked along the ADM1 dimension, in the order given).
bind_adm1_draws <- function(parts) {
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    return(NULL)
  }
  dims <- lapply(parts, function(p) dim(p$logrr))
  C <- dims[[1L]][2L]
  D <- dims[[1L]][3L]
  if (any(vapply(dims, function(x) x[2L] != C || x[3L] != D, logical(1)))) {
    stop(
      "bind_adm1_draws: groups disagree on the condition/draw dimensions",
      call. = FALSE
    )
  }
  A <- sum(vapply(dims, `[`, numeric(1), 1L))
  logrr <- array(
    NA_real_,
    dim = c(A, C, D),
    dimnames = list(NULL, dimnames(parts[[1L]]$logrr)[[2L]], NULL)
  )
  off <- 0L
  for (p in parts) {
    n <- dim(p$logrr)[1L]
    logrr[off + seq_len(n), , ] <- p$logrr
    off <- off + n
  }
  adm1 <- dplyr::bind_rows(lapply(parts, `[[`, "adm1"))
  dimnames(logrr)[[1L]] <- adm1$adm1_gid
  list(logrr = logrr, adm1 = adm1)
}

#' ADM1 x C x D array aligned onto a shared ADM1 universe (rows in `universe`
#' order); ADM1s absent from `adm1_gid` get `fill`. Used to OR per-track masks
#' whose fits carry different ADM1 random-effect levels.
align_adm1_array <- function(arr, adm1_gid, universe, fill = FALSE) {
  d <- dim(arr)
  out <- array(
    fill,
    dim = c(length(universe), d[2L], d[3L]),
    dimnames = list(universe, dimnames(arr)[[2L]], NULL)
  )
  hit <- match(universe, adm1_gid)
  ok <- which(!is.na(hit))
  if (length(ok)) {
    out[ok, , ] <- arr[hit[ok], , , drop = FALSE]
  }
  out
}

#' ENSO-sensitive mask: |log IRR| beyond log(1 + thr), either direction.
sensitive_mask <- function(logrr, thr) abs(logrr) > log1p(thr)

#' Directions of the ENSO exposure tables (derive_report_enso_exposure()).
ENSO_EXPOSURE_DIRECTIONS <- c("adverse", "protective", "either")

#' Exposure mask of one track's log-IRR array: beyond log(1 + m) in the adverse
#' direction (`adverse_sign` * x, TRACK_ADVERSE_SIGN), the protective one, or
#' either (sensitive_mask(), the sign ignored).
exposure_mask <- function(logrr, adverse_sign, direction, m) {
  switch(
    direction,
    adverse = adverse_sign * logrr > log1p(m),
    protective = -adverse_sign * logrr > log1p(m),
    either = sensitive_mask(logrr, m),
    stop("[09] unknown exposure direction '", direction, "'.", call. = FALSE)
  )
}

#' Population shares: for every population column present in `d` (pop,
#' pop_point, pop_mean, pop_lo, pop_med, pop_hi) add its share of the scope's
#' population `pop_universe` (pct, pct_point, pct_mean, ...; a fraction in
#' [0, 1]). The denominator is a constant per scope, so the share of a
#' posterior quantile is the quantile of the per-draw share.
add_pop_shares <- function(d, denom = "pop_universe") {
  for (col in intersect(
    c("pop", "pop_point", "pop_mean", "pop_lo", "pop_med", "pop_hi"),
    names(d)
  )) {
    d[[sub("^pop", "pct", col)]] <- ifelse(
      d[[denom]] > 0,
      d[[col]] / d[[denom]],
      NA_real_
    )
  }
  d
}

#' Observed sub-ADM1 population map: for every ADM1 unit, the summed population
#' of the sub-units ACTUALLY observed in the data (its ADM2 districts, or the
#' ADM1 itself when it is an ADM1-only unit), NOT the whole-province ADM1
#' population. This is the containment map the manuscript population tables join
#' on, so a significant/forecast effect in one district is credited only that
#' district's population, never the entire state (which double-counts geography
#' and inflates every population total above the actual coverage — the mismatch
#' vs. the R1 "total population covered" descriptor).
#'
#' `frame` supplies the observed gids (+ country/AFRO_region); the overlap filter
#' (drop_overlapping_adm1_from_frame) removes coarse ADM1 gids that also have nested
#' ADM2 children, mirroring R1 (Table 4) and R4. `recover_dropped = TRUE` re-adds
#' the handful of ADM1 units whose single gid STRING carries both an ADM1 and an
#' ADM2 row (so the string-level overlap filter removed the legitimate ADM2 obs
#' too, collapsing that unit's observed population to zero); their population is
#' recovered from their observed ADM2-level sub-units in surv_dat. R4 keeps
#' recover_dropped = FALSE for byte-identical output; R1/R2 set it TRUE so their
#' two population universes reconcile (Table-4 total rises to match).
#'
#' @return one row per ADM1 unit: adm1_unit, country, AFRO_region, pop,
#'   n_gid (contributing gids), n_pop_missing.
build_adm1_obs_pop <- function(
  ctx,
  frame,
  surv_dat,
  gid_pop,
  tag = "adm1-obs-pop",
  pop_ref_year = ctx$POP_REF_YEAR,
  recover_dropped = FALSE,
  wpp = NULL
) {
  fc <- drop_overlapping_adm1_from_frame(
    frame %>% dplyr::mutate(gid = as.character(gid)),
    surv_dat,
    tag
  ) %>%
    dplyr::mutate(adm1_unit = gid_to_adm1(gid, tag))
  prf <- resolve_gid_population(
    fc$gid,
    pop_ref_year,
    gid_pop,
    wpp = wpp,
    on_missing = "na"
  )
  fc$pop <- prf$pop
  fc$pop_missing <- prf$pop_tier == "missing"
  adm1_pop <- fc %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(
      country = dplyr::first(country),
      AFRO_region = dplyr::first(AFRO_region),
      pop = sum(pop, na.rm = TRUE),
      n_gid = dplyr::n(),
      n_pop_missing = sum(pop_missing),
      .groups = "drop"
    )
  if (isTRUE(recover_dropped)) {
    su <- surv_dat %>%
      dplyr::transmute(
        adm1_unit = as.character(ADM1),
        admin_level = as.character(admin_level),
        gid = as.character(gid),
        country = as.character(country),
        AFRO_region = as.character(AFRO_region)
      ) %>%
      dplyr::filter(!is.na(adm1_unit))
    dropped <- setdiff(unique(su$adm1_unit), adm1_pop$adm1_unit)
    if (length(dropped)) {
      rec_src <- su %>%
        dplyr::filter(adm1_unit %in% dropped, admin_level == "ADM2") %>%
        dplyr::distinct(adm1_unit, gid, country, AFRO_region)
      rp <- resolve_gid_population(
        rec_src$gid,
        pop_ref_year,
        gid_pop,
        wpp = NULL,
        on_missing = "na"
      )
      rec_src$pop <- rp$pop
      rec_src$pop_missing <- rp$pop_tier == "missing"
      rec <- rec_src %>%
        dplyr::group_by(adm1_unit) %>%
        dplyr::summarise(
          country = dplyr::first(country),
          AFRO_region = dplyr::first(AFRO_region),
          pop = sum(pop, na.rm = TRUE),
          n_gid = dplyr::n(),
          n_pop_missing = sum(pop_missing),
          .groups = "drop"
        )
      adm1_pop <- dplyr::bind_rows(adm1_pop, rec)
      message(sprintf(
        paste0(
          "[%s] recovered %d over-dropped ADM1 unit(s) from ",
          "observed ADM2 sub-units (+%.2fM)."
        ),
        tag,
        nrow(rec),
        sum(rec$pop, na.rm = TRUE) / 1e6
      ))
    }
  }
  adm1_pop
}

#' Named observed-sub-unit population per ADM1 unit (names = adm1_unit) over
#' the occurrence survival dataset's containment universe, recover_dropped =
#' TRUE: the basis every ENSO-effect population table (significant-effect,
#' ENSO-sensitive) is credited on, so the two can never drift apart.
obs_adm1_pop_map <- function(ctx, surv_dat, gid_pop, tag) {
  adm1_obs <- build_adm1_obs_pop(
    ctx,
    surv_dat %>% dplyr::distinct(gid, country, AFRO_region),
    surv_dat,
    gid_pop,
    tag,
    ctx$POP_REF_YEAR,
    recover_dropped = TRUE
  )
  stats::setNames(adm1_obs$pop, adm1_obs$adm1_unit)
}

#' Observed sub-ADM1 population of the simulated forecast universe
#' (report_forecast, report_attrib_exposure): build_adm1_obs_pop() with
#' recover_dropped = FALSE (the pre-helper forecast universe), joined to each
#' simulated ADM1's ENSO group and its short region name.
build_forecast_adm1_pop <- function(
  ctx,
  frame,
  surv_dat,
  gid_pop,
  adm1_month,
  tag = "report-forecast"
) {
  out <- build_adm1_obs_pop(
    ctx,
    frame,
    surv_dat,
    gid_pop,
    tag,
    ctx$POP_REF_YEAR,
    recover_dropped = FALSE
  ) %>%
    dplyr::left_join(
      dplyr::distinct(adm1_month, adm1_unit, enso_group),
      by = "adm1_unit"
    ) %>%
    dplyr::mutate(region = unname(ctx$GROUP_SHORT_RUN[enso_group]))
  # rate_pop_year = "window": the per-capita rates' population
  # `pop_rate`, the month-weighted mean over the months of `adm1_month`; `pop`
  # (people living in) stays at POP_REF_YEAR. Absent under "ref".
  if (identical(ctx[["RATE_POP_YEAR"]], "window")) {
    w <- rate_pop_weights(adm1_month$month)
    wpp <- load_wpp_table(.pop_artefact_or_null(ctx$WPP_PATH))
    pr <- Reduce(
      `+`,
      lapply(names(w), function(y) {
        p <- build_adm1_obs_pop(
          ctx,
          frame,
          surv_dat,
          gid_pop,
          tag,
          as.integer(y),
          recover_dropped = FALSE,
          wpp = wpp
        )
        w[[y]] * p$pop[match(out$adm1_unit, p$adm1_unit)]
      })
    )
    out$pop_rate <- pr
  }
  out
}

#' .sim_pop_year_fields
#'
#' Additive provenance fields of a simulation-derived payload (no schema
#' bump): `sim_pop_year`, the simulator's
#' simulator.location_frame.pop_year as the manifests record it
#' (make_*_configs.R write it only when not the default 2020), and
#' `rate_pop_year` under "window". Both absent at the defaults, so a default
#' payload is unchanged. The manifests given (e.g. the displayed forecast's and
#' the attribution factual's) must agree.
#' @return a list to c() onto the payload.
.sim_pop_year_fields <- function(ctx, mans) {
  vals <- vapply(
    Filter(Negate(is.null), mans),
    function(m) as.character(m[["pop_year"]] %||% .SIM_POP_YEAR_DEFAULT),
    character(1)
  )
  if (length(unique(vals)) > 1L) {
    stop(
      "[09] the manifests read by one step record different simulator population years (",
      paste(vals, collapse = " vs "),
      "); the arms must share one pop_year.",
      call. = FALSE
    )
  }
  out <- list()
  if (
    length(vals) && !identical(vals[[1L]], as.character(.SIM_POP_YEAR_DEFAULT))
  ) {
    out$sim_pop_year <- vals[[1L]]
  }
  if (identical(ctx[["RATE_POP_YEAR"]], "window")) {
    out$rate_pop_year <- "window"
  }
  out
}

#' rate_pop_weights
#'
#' Calendar-year weights of the months `months` reads (each distinct month
#' counts once): named numeric (names = years) summing to 1, the weights of
#' the month-weighted mean population of rate_pop_year = "window". An empty / all-NA
#' `months` is an error.
rate_pop_weights <- function(months) {
  m <- unique(lubridate::floor_date(as.Date(months[!is.na(months)]), "month"))
  if (!length(m)) {
    stop(
      "[09] rate_pop_weights: no month to weight the population over.",
      call. = FALSE
    )
  }
  tab <- table(format(m, "%Y"))
  stats::setNames(as.numeric(tab) / sum(tab), names(tab))
}

#' The window totals behind pair_adm1_cum() for every analog set and both
#' metrics at once: cases / outbreaks per (run, analog_set, adm1_unit,
#' replicate_id) of `adm1_month` (one .sum_by_keys() pass per cell).
pair_adm1_cum_table <- function(adm1_month) {
  .sum_by_keys(
    adm1_month,
    c("run", "analog_set", "adm1_unit", "replicate_id"),
    c(cases = "cases", outbreaks = "outbreaks")
  )
}

#' Forecast vs counterfactual cumulative totals per ADM1 and replicate, paired by
#' replicate_id: the forecast arm against counterfactual set `aset`, `value_col`
#' summed over every row of the (already windowed) adm1_month frame.
#' @return tibble(adm1_unit, replicate_id, f, cf)
pair_adm1_cum <- function(adm1_month, value_col, aset, cum = NULL) {
  # `cum` (pair_adm1_cum_table() of the same rows): the per-(run, set, ADM1,
  # replicate) totals summed once per cell for every metric and set; the ~10^8-row
  # table is then not filtered and summed again per (metric, set).
  cum_adm1 <- if (is.null(cum)) {
    .sum_by_keys(
      dplyr::filter(adm1_month, run == "forecast" | analog_set == aset),
      c("run", "adm1_unit", "replicate_id"),
      c(value = value_col)
    )
  } else {
    .sum_by_keys(
      cum[cum$run == "forecast" | cum$analog_set %in% aset, ],
      c("run", "adm1_unit", "replicate_id"),
      c(value = value_col)
    )
  }
  dplyr::inner_join(
    dplyr::filter(cum_adm1, run == "forecast") %>%
      dplyr::transmute(adm1_unit, replicate_id, f = value),
    dplyr::filter(cum_adm1, run == "counterfactual") %>%
      dplyr::transmute(adm1_unit, replicate_id, cf = value),
    by = c("adm1_unit", "replicate_id")
  )
}

#' Per-draw directional-population 95% CrI, the posterior companion to the
#' significant-effect population sums in derive_report_effects(). For every
#' retained posterior draw and each ONI condition we sum the population of every
#' ADM1 whose combined log-IRR is > 0 (increase) or < 0 (decrease) on that draw,
#' then take median + 95% quantiles ACROSS draws.
#'
#' Reads the per-draw log-IRR array cached in the enso_relation artefact
#' (`adm1_draws`, extract_adm1_logrr_draws(): the same seed, subsample and
#' crossbasis rows as extract_adm1_irr()), so no fit is reloaded and the point
#' estimate and this CrI come off one shared posterior.
#'
#' @param draws enso_relation artefact's `adm1_draws` (list(logrr, adm1))
#' @param pop_map named numeric: observed-sub-unit population per adm1_unit
dir_pop_ci <- function(
  draws,
  group_label,
  pop_map,
  enso_conds,
  group_short = GROUP_SHORT
) {
  sel <- which(draws$adm1$enso_group == group_label)
  if (!length(sel)) {
    return(NULL)
  }
  arr <- draws$logrr[sel, , , drop = FALSE]
  adm1 <- draws$adm1[sel, , drop = FALSE]
  pops <- unname(pop_map[adm1$adm1_gid])

  region_inc <- pop_mass_by_draw(arr > 0, pops)
  region_dec <- pop_mass_by_draw(arr < 0, pops)
  ctry <- unique(adm1$country)
  ctry_inc <- lapply(stats::setNames(ctry, ctry), function(cc) {
    k <- adm1$country == cc
    pop_mass_by_draw(arr[k, , , drop = FALSE] > 0, pops[k])
  })
  ctry_dec <- lapply(stats::setNames(ctry, ctry), function(cc) {
    k <- adm1$country == cc
    pop_mass_by_draw(arr[k, , , drop = FALSE] < 0, pops[k])
  })

  q_rows <- function(mat, direction) .pop_ci_rows(mat, direction, enso_conds)
  region_ci <- dplyr::bind_rows(
    q_rows(region_inc, "increase"),
    q_rows(region_dec, "decrease")
  ) %>%
    dplyr::mutate(
      enso_group = group_label,
      region = unname(group_short[group_label])
    )
  country_ci <- dplyr::bind_rows(lapply(ctry, function(cc) {
    dplyr::bind_rows(
      q_rows(ctry_inc[[cc]], "increase"),
      q_rows(ctry_dec[[cc]], "decrease")
    ) %>%
      dplyr::mutate(country = cc)
  })) %>%
    dplyr::mutate(
      enso_group = group_label,
      region = unname(group_short[group_label])
    )
  # region_inc/region_dec (C x D) are the group's per-draw population mass; the
  # caller sums them element-wise across groups (pairing draw index across the
  # independent group posteriors) to build the Africa-wide per-draw total.
  list(
    region_ci = region_ci,
    country_ci = country_ci,
    region_inc = region_inc,
    region_dec = region_dec
  )
}

derive_enso_relation <- function(ctx, track) {
  man <- ctx$MAN
  groups <- manifest_groups_all(man)
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  message(sprintf(
    "[09] ENSO exposure-response (%s): deriving ERF + ADM1 IRR for groups: %s",
    track,
    paste(groups, collapse = ", ")
  ))

  fit_files <- character(0)
  meta_files <- character(0)
  present <- character(0)
  for (g in groups) {
    ff <- track_fit_files(mnf, track, g)
    fit_files <- c(fit_files, ff$fit)
    meta_files <- c(meta_files, ff$meta)
    if (!file.exists(ff$fit) || !file.exists(ff$meta)) {
      warning(sprintf(
        "[09] ENSO exposure-response (%s): missing fit/meta for %s: %s",
        track,
        g,
        ff$fit
      ))
      next
    }
    present <- c(present, g)
  }
  # One fit per group, loaded and reduced on its own worker (ctx$WORKERS;
  # serial when 1). Every extractor seeds its own draw subsample, so the result
  # does not depend on the worker count.
  per_group <- map_workers(
    present,
    function(g) {
      ff <- track_fit_files(mnf, track, g)
      message(sprintf("[09]   loading %s | %s", track, g))
      brms_fit <- readRDS(ff$fit)
      meta <- readRDS(ff$meta)
      vm <- meta[["enso"]]
      cb_names <- vm$cb_col_names
      grid <- build_onebasis_grid(vm, n_at = 200L)

      b_pop_draws <- brms::as_draws_matrix(
        brms_fit,
        variable = sprintf("b_%s", cb_names)
      ) |>
        unclass()
      colnames(b_pop_draws) <- cb_names
      out <- list(
        erf = erf_from_draws(b_pop_draws, grid$Z_cen, grid$at_grid) %>%
          dplyr::mutate(
            enso_group = g,
            range_lo = vm$range[1],
            range_hi = vm$range[2]
          ),
        adm1 = extract_adm1_irr(
          brms_fit,
          meta,
          enso_conds = ctx$ENSO_CONDS,
          group_label = g,
          n_draws = ctx$N_IRR_DRAWS
        ),
        # Per-draw log-IRR array (same seed/subsample as the summaries above),
        # cached so dir_pop_ci() and the ENSO-sensitive population never reload a fit.
        draws = extract_adm1_logrr_draws(
          brms_fit,
          meta,
          enso_conds = ctx$ENSO_CONDS,
          group_label = g,
          n_draws = ctx$N_IRR_DRAWS
        ),
        country_erf = extract_country_erf(
          brms_fit,
          meta,
          grid,
          group_label = g,
          n_draws = ctx$N_IRR_DRAWS
        )
      )
      rm(brms_fit)
      gc(verbose = FALSE)
      out
    },
    workers = ctx$WORKERS
  )
  names(per_group) <- present
  pick <- function(k) lapply(per_group, `[[`, k)
  erf_list <- pick("erf")
  adm1_list <- pick("adm1")
  draws_list <- pick("draws")
  country_erf_list <- pick("country_erf")
  erf_all <- dplyr::bind_rows(erf_list)
  adm1_all <- dplyr::bind_rows(adm1_list)
  country_erf_all <- dplyr::bind_rows(country_erf_list)
  # NOTE the duration track is NO LONGER inverted here: the ENSO exposure-response
  # artefacts now carry the duration model's own estimand, the CESSATION-HAZARD
  # IRR (high = ends sooner = SHORTER outbreak), the same scale 06_fit_figures.R
  # plots it on. Any orientation flip (invert_duration_irr(), above) is the
  # figure layer's business.
  if (!nrow(erf_all)) {
    stop(sprintf(
      "[09] ENSO exposure-response (%s): no ERF rows extracted",
      track
    ))
  }
  if (!nrow(adm1_all)) {
    stop(sprintf(
      "[09] ENSO exposure-response (%s): no ADM1 IRR rows extracted",
      track
    ))
  }
  if (!nrow(country_erf_all)) {
    warning(sprintf(
      "[09] ENSO exposure-response (%s): no country ERF rows extracted",
      track
    ))
  }

  lag_rank <- lag_rank_with_selection(ctx, track, groups, mnf)$lag_rank

  list(
    erf = erf_all,
    adm1_irr = adm1_all,
    lag_rank = lag_rank,
    enso_grp = load_enso_grouping(ctx$GROUPING_PATH),
    groups = groups,
    country_erf = country_erf_all,
    adm1_draws = bind_adm1_draws(draws_list),
    enso_conds = ctx$ENSO_CONDS,
    n_draws = ctx$N_IRR_DRAWS,
    .inputs = results_inputs(files = c(fit_files, meta_files))
  )
}

#' Deployed best-lag fit + cb_meta paths of one track over every manifest
#' group, for an ensure_*() wrapper's freshness inputs. Lenient: an unreadable
#' archive or a missing (track, group) row contributes nothing (the deriver
#' itself fails loudly).
#' @return list(fit = character, meta = character)
track_fit_input_files <- function(ctx, track) {
  groups <- manifest_groups_all(ctx$MAN)
  mnf <- tryCatch(
    load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR),
    error = function(e) NULL
  )
  ff <- if (!is.null(mnf)) {
    lapply(groups, function(g) {
      tryCatch(track_fit_files(mnf, track, g), error = function(e) NULL)
    })
  } else {
    list()
  }
  ff <- Filter(Negate(is.null), ff)
  list(
    fit = vapply(ff, `[[`, character(1), "fit"),
    meta = vapply(ff, `[[`, character(1), "meta")
  )
}

#' Freshness-guarded ENSO exposure-response artefact: derive only if
#' stale/absent/decision-diff.
ensure_enso_relation <- function(ctx, track) {
  tf <- track_fit_input_files(ctx, track)
  fit_files <- tf$fit
  meta_files <- tf$meta
  lag_loo_files <- lag_loo_input_files(ctx, track)
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$GROUPING_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      fit_files,
      meta_files,
      lag_loo_files
    )
  )

  ensure_artefact(
    ctx,
    "enso_relation",
    inputs,
    function() derive_enso_relation(ctx, track),
    track = track,
    strip_inputs = TRUE
  )
}

# ===========================================================================
# Section 1b: region/country parameter summaries (supplementary plots)
# ===========================================================================

#' Tidy population + country parameter summaries from each ENSO-group best-lag
#' fit, for every model term EXCEPT ENSO (that is Figures 1/2). Reuses the ENSO
#' path's per-group fit-resolution loop and the brms_group_re_scaffold()
#' draw-loading contract; the only difference from the ENSO extractors is that
#' the parameters are scalar (identity projection, no centred onebasis) and the
#' tier accumulation runs over fixef rownames instead of a crossbasis.
#'
#' Country tier: for every term that carries a country-level RE (intercept + the
#' configured climate slopes) accumulate b_pop + region-RE + country-RE per draw
#' and summarise exp(.) to a mean + 95% CrI. Terms with no country RE (spillover,
#' ttlo, seasonal, month) appear only at the population tier (brms::fixef() with
#' 95% CrI, exponentiated). Under grouping the region tier is NULL (the ENSO
#' group is the region facet), so re_draws_mat("region", .) is a zero matrix.
derive_track_params <- function(ctx, track) {
  groups <- manifest_groups_all(ctx$MAN)
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  message(sprintf(
    "[09] track params (%s): extracting region/country parameters for groups: %s",
    track,
    paste(groups, collapse = ", ")
  ))

  fit_files <- character(0)
  present <- character(0)
  for (g in groups) {
    ff <- track_fit_files(mnf, track, g)
    fit_files <- c(fit_files, ff$fit)
    if (!file.exists(ff$fit)) {
      warning(sprintf(
        "[09] track params (%s): missing fit for %s: %s",
        track,
        g,
        ff$fit
      ))
      next
    }
    present <- c(present, g)
  }
  # One fit per group on its own worker (ctx$WORKERS; serial when 1). Each
  # returns its list of row tibbles; the draw subsample is seeded per fit.
  per_group <- map_workers(
    present,
    function(g) {
      ff <- track_fit_files(mnf, track, g)
      rows <- list()
      message(sprintf("[09]   loading %s | %s", track, g))
      fit <- readRDS(ff$fit)

      all_terms <- rownames(brms::fixef(fit))
      fam <- vapply(all_terms, param_family, character(1))
      variables <- all_terms[!is.na(fam)]
      if (!length(variables)) {
        warning(sprintf(
          "[09] track params (%s): no non-ENSO terms for %s",
          track,
          g
        ))
        return(rows)
      }

      sc <- brms_group_re_scaffold(fit, variables, n_draws = ctx$N_IRR_DRAWS)

      # --- population tier: fixef mean + 95% CrI, exponentiated ----------------
      fx <- brms::fixef(fit, probs = c(0.025, 0.975))[variables, , drop = FALSE]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        track = track,
        enso_group = g,
        tier = "population",
        level = g,
        term = variables,
        family = unname(fam[variables]),
        mean = exp(fx[, "Estimate"]),
        q025 = exp(fx[, "Q2.5"]),
        q975 = exp(fx[, "Q97.5"])
      )

      # --- country tier: b_pop + region-RE + country-RE, terms with a country RE
      re_ctry_arr <- sc$re_mean[[sc$grp$country]]
      if (is.null(re_ctry_arr)) {
        warning(sprintf(
          "[09] track params (%s): country RE not found for %s",
          track,
          g
        ))
        rm(fit)
        gc(verbose = FALSE)
        return(rows)
      }
      ctry_levels <- dimnames(re_ctry_arr)[[1L]]
      ctry_vars <- intersect(variables, dimnames(re_ctry_arr)[[3L]])
      if (length(ctry_vars)) {
        for (lbl in ctry_levels) {
          parsed <- sc$parse_level(lbl, "country")
          comb_d <- sc$b_pop_draws +
            sc$re_draws_mat("region", parsed$region) +
            sc$re_draws_mat("country", lbl) # n_draws x K
          ratio_d <- exp(comb_d[, ctry_vars, drop = FALSE])
          rows[[length(rows) + 1L]] <- tibble::tibble(
            track = track,
            enso_group = g,
            tier = "country",
            level = parsed$country,
            term = ctry_vars,
            family = unname(fam[ctry_vars]),
            mean = colMeans(ratio_d),
            q025 = apply(ratio_d, 2L, quantile, 0.025),
            q975 = apply(ratio_d, 2L, quantile, 0.975)
          )
        }
      }
      rm(fit)
      gc(verbose = FALSE)
      rows
    },
    workers = ctx$WORKERS
  )
  rows <- do.call(c, per_group)

  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    stop(sprintf("[09] track params (%s): no parameter rows extracted", track))
  }
  out
}

#' Freshness-guarded region/country parameter artefact: derive only if
#' stale/absent/decision-diff. Inputs mirror ensure_enso_relation (same
#' best-lag archive), minus the lag-LOO rankings (not used here).
ensure_track_params <- function(ctx, track) {
  fit_files <- track_fit_input_files(ctx, track)$fit
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$GROUPING_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      fit_files
    )
  )

  ensure_artefact(
    ctx,
    "track_params",
    inputs,
    function() derive_track_params(ctx, track),
    track = track
  )
}

# ===========================================================================
# Section 1c: per-gid outbreak-onset rate (Figure 1b inset map)
# ===========================================================================

#' Onsets per surveillance-covered month for every gid in the survival panel.
#'
#' Denominator: 12 x the gid's covered season-years. Every covered season-year
#' has rows in the panel (under gid_season, gid_at_risk_years() in
#' 03_prepare_occ.R unions the member LPs' coverage in whole season-years;
#' under lp_season, per LP), and a season-year starts at t_start - t_start_seas,
#' rounded to the month. Covered months therefore include months during an
#' outbreak, unlike the at-risk row count.
#' Numerator: distinct onset months per gid -- the gid x month union of the
#' descriptive notebook's Table 1b. Under gid_season this equals the onset
#' count (one history per gid); under lp_season it collapses the onsets that
#' several location-period histories put on the same gid-month.
derive_onset_rate <- function(ctx) {
  surv <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  grp <- load_enso_grouping(ctx$GROUPING_PATH)
  out <- surv %>%
    dplyr::mutate(
      month = lubridate::floor_date(as.Date(t_start), "month"),
      season_year = lubridate::round_date(
        as.Date(t_start) - t_start_seas,
        "month"
      ),
      outbreak = as.logical(outbreak)
    ) %>%
    dplyr::group_by(gid, country) %>%
    dplyr::summarise(
      onsets = dplyr::n_distinct(month[outbreak %in% TRUE]),
      covered_months = 12L * dplyr::n_distinct(season_year),
      .groups = "drop"
    ) %>%
    dplyr::left_join(grp, by = c("country" = "iso3")) %>%
    dplyr::rename(enso_group = group) %>%
    dplyr::mutate(rate = onsets / covered_months)
  if (anyDuplicated(out$gid)) {
    stop("[09] onset rate: a gid maps to more than one country")
  }
  if (any(out$rate > 1)) {
    stop(
      "[09] onset rate: ",
      sum(out$rate > 1),
      " gid(s) with more onset months than covered months"
    )
  }
  message(sprintf(
    "[09] onset rate: %d gids, median %.4f, max %.4f onsets per covered month",
    nrow(out),
    stats::median(out$rate),
    max(out$rate)
  ))
  out
}

#' Freshness-guarded onset-rate artefact.
ensure_onset_rate <- function(ctx) {
  inputs <- results_inputs(files = c(ctx$SURV_PATH, ctx$GROUPING_PATH))
  ensure_artefact(ctx, "onset_rate", inputs, function() derive_onset_rate(ctx))
}

# ===========================================================================
# Section 1c': the modelled units and the fit membership (supplement maps MU,
# GM, KC)
# ===========================================================================

#' Status of every modelled unit and of every ADM1 of a modelled country
#' without one (fig. MU).
#' A country is modelled when it has a fit group in `rate`; its units are the
#' gids of the survival panel. Admin level: gid_level() of the flat gid with
#' Ghana's "GHA5_2" read as "GHA_5_2"; ADM1 key: gid.R::gid_to_adm1().
#'
#' @param rate onset_rate artefact (gid, country, onsets, enso_group)
#' @param adm1_gids GID_1 values of the GADM ADM_1 layer (any form gid_to_adm1() parses)
#' @return list(units = tibble(gid, country, enso_group, level, adm1, status in
#'   adm1 | adm2 | no_onset), gaps = tibble(adm1, country, status = "no_episode"))
modelled_units_tables <- function(rate, adm1_gids) {
  r <- rate[!is.na(rate$enso_group), , drop = FALSE]
  flat <- sub("^([A-Z]{3})([0-9])", "\\1_\\2", flatten_gid(as.character(r$gid)))
  lvl <- gid_level(flat)
  units <- tibble::tibble(
    gid = as.character(r$gid),
    country = as.character(r$country),
    enso_group = as.character(r$enso_group),
    level = lvl,
    adm1 = gid_to_adm1(r$gid, tag = "modelled_units", strict = FALSE),
    status = ifelse(
      r$onsets > 0,
      ifelse(lvl == "ADM2", "adm2", "adm1"),
      "no_onset"
    )
  )
  ids <- unique(stats::na.omit(gid_to_adm1(
    adm1_gids,
    tag = "modelled_units",
    strict = FALSE
  )))
  ids <- sort(ids[
    substr(ids, 1, 3) %in% unique(units$country) & !ids %in% units$adm1
  ])
  list(
    units = units,
    gaps = tibble::tibble(
      adm1 = ids,
      country = substr(ids, 1, 3),
      status = rep("no_episode", length(ids))
    )
  )
}

#' The GID_1 column of the GADM ADM_1 layer (attributes only).
.gadm_adm1_gids <- function(path) {
  d <- suppressWarnings(sf::st_read(
    require_input(path, "GADM geopackage"),
    quiet = TRUE,
    query = 'SELECT GID_1 FROM "ADM_1"'
  ))
  as.character(d$GID_1)
}

#' modelled_units (v1, unkeyed): the analysis set by unit (fig. MU) and the fit
#' registry's membership (figs GM, KC). `members` is the deployed partition
#' (ctx$GROUPING_PATH): 10 takes country -> group from here, never from the
#' registry CSV.
#' @return list(units, gaps (modelled_units_tables()), members = tibble(iso3, group))
derive_modelled_units <- function(ctx) {
  rate <- ensure_onset_rate(ctx)
  out <- modelled_units_tables(rate, .gadm_adm1_gids(ctx$GADM_ADM1_PATH))
  grp <- load_enso_grouping(ctx$GROUPING_PATH)
  out$members <- tibble::tibble(
    iso3 = as.character(grp$iso3),
    group = as.character(grp$group)
  )
  message(sprintf(
    "[09] modelled units: %d units (%s), %d ADM1s without a recorded episode, %d members",
    nrow(out$units),
    paste(
      sprintf(
        "%s %d",
        names(table(out$units$status)),
        as.integer(table(out$units$status))
      ),
      collapse = ", "
    ),
    nrow(out$gaps),
    nrow(out$members)
  ))
  out
}

#' Freshness-guarded modelled-units artefact.
ensure_modelled_units <- function(ctx) {
  inputs <- results_inputs(
    files = c(ctx$SURV_PATH, ctx$GROUPING_PATH, ctx$GADM_ADM1_PATH)
  )
  ensure_artefact(ctx, "modelled_units", inputs, function() {
    derive_modelled_units(ctx)
  })
}

# --- Figure 1c's onset-season column ----------------------------------------
#
# occ_seasonality (v1, unkeyed): the occurrence track's seasonal baseline per
# ENSO group and calendar month. Per group, the deployed occurrence fit's
# cyclic spline of t_start_seas (utils/formula.R `cyclic`) at the 12
# season-month starts, exp(s - mean(s)) per draw (relative to the annual
# mean), posterior mean; each gid of the group's countries rotates it to
# calendar months by its season start min_seas; the mean over the gids.
# v2 (Figure 1c seasonality rows): + the
# per-draw curves behind that mean (`draws`, per group a draws x 12 matrix in
# calendar months, the gid mean taken per draw), so seasonality can carry an
# occurrence band at its own `interval`; `rel` unchanged.

# t_start_seas (days since the season year's start, occ_data.R) at the start of
# season months 1..12 in a non-leap year; findInterval() on it gives the season
# month of any row (leap-year offsets never cross a month boundary).
.OCC_SEAS_MONTH_STARTS <- c(
  0L,
  31L,
  59L,
  90L,
  120L,
  151L,
  181L,
  212L,
  243L,
  273L,
  304L,
  334L
)

#' Season start of every gid of the survival panel: the calendar month at which
#' its cholera season year starts (occ_data.R: t_start_seas = yday(t_start
#' rotated by -(min_seas - 1)) - 1), modal over its rows.
#' @return tibble(gid, country, min_seas)
occ_gid_min_seas <- function(surv) {
  surv %>%
    dplyr::transmute(
      gid = as.character(gid),
      country = as.character(country),
      cal_month = lubridate::month(as.Date(t_start)),
      seas_month = findInterval(
        as.numeric(t_start_seas),
        .OCC_SEAS_MONTH_STARTS
      ),
      min_seas = ((cal_month - seas_month) %% 12L) + 1L
    ) %>%
    dplyr::count(gid, country, min_seas) %>%
    dplyr::group_by(gid, country) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(gid, country, min_seas)
}

#' The brms label of a fit's cyclic t_start_seas smooth (the formula term with
#' its whitespace removed, as posterior_smooths() expects it).
occ_seas_smooth_label <- function(fit) {
  lab <- attr(stats::terms(stats::formula(fit$formula$formula)), "term.labels")
  lab <- gsub("[[:space:]]+", "", lab)
  hit <- lab[startsWith(lab, "s(t_start_seas")]
  if (length(hit) != 1L) {
    stop(
      "[09] occ_seasonality: expected one s(t_start_seas, ...) term in the occurrence formula; found ",
      length(hit),
      ".",
      call. = FALSE
    )
  }
  hit
}

#' Each group's deployed occurrence fit file (read_deployed_lag() on the
#' best-lag dir); an error when one is missing.
occ_seasonality_fit_files <- function(ctx, groups) {
  grouping <- as.character(
    ctx$MAN$grouping %||% sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  )
  stats::setNames(
    vapply(
      groups,
      function(g) {
        dl <- read_deployed_lag(
          ctx$ENSO_BEST_LAG_DIR,
          g,
          grouping,
          "occ",
          episode_key = episode_key_of(ctx$opt),
          national_zero_fill = national_zero_fill_of(ctx$opt),
          national_zero_threshold = national_zero_threshold_of(ctx$opt),
          occ_panel = occ_panel_spec_of(ctx$opt)
        )
        if (is.null(dl)) {
          stop(
            sprintf(
              "[09] occ_seasonality: no deployed occurrence fit for group %s in %s (grouping %s).",
              g,
              ctx$ENSO_BEST_LAG_DIR,
              grouping
            ),
            call. = FALSE
          )
        }
        dl$file
      },
      character(1)
    ),
    groups
  )
}

#' A group's season-month curves (draws x 12) rotated to calendar months per
#' gid by its min_seas and averaged over the gids, per draw (pure).
#' @param r matrix draws x 12, season months 1..12
#' @param min_seas integer, one per gid
#' @return matrix draws x 12, calendar months 1..12
occ_seas_calendar_draws <- function(r, min_seas) {
  m <- matrix(0, nrow(r), 12L)
  for (ms in min_seas) {
    idx <- ((ms - 1L) + (seq_len(12L) - 1L)) %% 12L + 1L # season month j -> calendar month
    m[, idx] <- m[, idx] + r
  }
  m / length(min_seas)
}

#' occ_seasonality v2.
#' @return list(seasonality = tibble(enso_group, month, rel, n_gid), draws =
#'   named list per group of draws x 12 matrices (calendar months), fits =
#'   tibble(enso_group, file, smooth), best_lag_dir, grouping)
derive_occ_seasonality <- function(ctx) {
  groups <- as.character(results_group_meta(ctx)$order)
  files <- occ_seasonality_fit_files(ctx, groups)
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  cal <- occ_gid_min_seas(readRDS(ctx$SURV_PATH))
  smooths <- character(0)
  draws <- list()
  rows <- lapply(groups, function(g) {
    fit <- readRDS(files[[g]])
    sm <- occ_seas_smooth_label(fit)
    smooths[[g]] <<- sm
    s <- brms::posterior_smooths(
      fit,
      smooth = sm,
      newdata = data.frame(t_start_seas = .OCC_SEAS_MONTH_STARTS)
    ) # draws x 12
    gids <- dplyr::filter(
      cal,
      country %in% as.character(reg$iso3[reg$group == g])
    )
    if (!nrow(gids)) {
      stop(
        "[09] occ_seasonality: no gid of group ",
        g,
        " in the survival panel ",
        ctx$SURV_PATH,
        call. = FALSE
      )
    }
    m <- occ_seas_calendar_draws(exp(sweep(s, 1, rowMeans(s))), gids$min_seas)
    draws[[g]] <<- m
    tibble::tibble(
      enso_group = g,
      month = seq_len(12L),
      rel = colMeans(m),
      n_gid = dplyr::n_distinct(gids$gid)
    )
  })
  out <- dplyr::bind_rows(rows)
  message(sprintf(
    "[09] occ_seasonality: %d groups, rel %.2f..%.2f (%s)",
    length(groups),
    min(out$rel),
    max(out$rel),
    paste(basename(files), collapse = ", ")
  ))
  list(
    seasonality = out,
    draws = draws,
    fits = tibble::tibble(
      enso_group = groups,
      file = unname(files),
      smooth = unname(smooths[groups])
    ),
    best_lag_dir = ctx$ENSO_BEST_LAG_DIR,
    grouping = ctx$GROUPING_PATH
  )
}

#' Freshness-guarded occurrence-seasonality artefact.
ensure_occ_seasonality <- function(ctx) {
  groups <- as.character(results_group_meta(ctx)$order)
  inputs <- results_inputs(
    files = c(
      occ_seasonality_fit_files(ctx, groups),
      ctx$SURV_PATH,
      ctx$GROUPING_PATH
    )
  )
  ensure_artefact(ctx, "occ_seasonality", inputs, function() {
    derive_occ_seasonality(ctx)
  })
}

# --- Supplement fig. SE: seasonality against the rainy season ---
#
# seasonality (v1, unkeyed), the reviewer's "does cholera peak in the rainy
# season of every region, and do case counts, outbreak frequency and cessation
# share one seasonal pattern". Per ENSO group and calendar month:
#   rain     the group's population-weighted rainfall climatology (mm) and
#            whether the month is in its main rainy season: rain_enso$month,
#            weighting "pop" (the Fig. 1d outlines; read, never recomputed);
#   rates    the observed rate of each track on the rows the fits keep
#            (results.derivation.seasonality.admin_level through
#            gid.R::apply_admin_level_filter()), by calendar month:
#              occurrence  onsets per 1,000 at-risk unit-months (survival
#                          panel: `outbreak` rows over exposure_days at risk,
#                          month of t_start);
#              cessation   episode ends per 100 active episode-months
#                          (duration panel: `ended` over exposure days);
#              cases       median suspected cases per active episode-month and
#                          its IQR (monthly panel; cases_gap months left out,
#                          their 0 is a reporting gap, not an observation).
#            Occurrence and cessation carry exact Poisson intervals at
#            `interval`. A month with fewer than min_n events (cases:
#            episode-months) is `sparse` and takes no part in the peak;
#   fitted   the occurrence model's seasonal spline (09's occ_seasonality; v2:
#            with its interval, from occ_seasonality's draws) and the monthly
#            model's calendar harmonics
#            (deployed mon fit per group: exp(h - mean h) per draw, posterior
#            mean and interval), each also `scaled` to the observed mean rate
#            of its track so 10 draws them without arithmetic. The duration
#            model has no seasonal term, so cessation has no fitted curve;
#   peaks    per group x track (observed and fitted) the peak month, whether
#            it is a rainy-season month, and the signed circular lag in months
#            from the centre of the nearest rainy season (positive: after it).
# The occurrence panel's load-time masks of the fits (national zero fill,
# censor_orphan) are not re-applied: the rates describe the panel's rows.

# Days per calendar month: converts exposure days to unit-months.
.SEAS_DAYS_PER_MONTH <- 365.25 / 12

#' Exact Poisson interval of a rate k / e at `level` (0 lower bound at k = 0).
seasonality_pois_ci <- function(k, e, level = 0.95) {
  a <- (1 - level) / 2
  list(
    lo = ifelse(k > 0, stats::qgamma(a, k), 0) / e,
    hi = stats::qgamma(1 - a, k + 1) / e
  )
}

#' Observed monthly rates of the three tracks per ENSO group (pure).
#' @param surv survival panel rows (country, admin_level, ADM1, t_start, outbreak, exposure_days)
#' @param dur duration panel rows (country, admin_level, ADM1, month, ended, exposure)
#' @param mon monthly panel rows (country, admin_level, ADM1, month, cases_month, cases_gap)
#' @param country_group named character, ISO3 -> group id
#' @param groups group ids to keep, in display order
#' @param admin_level apply_admin_level_filter() option
#' @param min_n events (cases: episode-months) below which a month is sparse
#' @param interval interval level
#' @return tibble(enso_group, track, month, rate, lo, hi, n, exposure, sparse);
#'   occurrence per 1,000 unit-months, cessation per 100 episode-months, cases
#'   per active episode-month (rate = median, lo / hi = quartiles, n = months).
seasonality_track_rates <- function(
  surv,
  dur,
  mon,
  country_group,
  groups,
  admin_level = "ADM1-ADM2",
  min_n = 5L,
  interval = 0.95
) {
  keep <- function(d) {
    d <- apply_admin_level_filter(d, admin_level)
    d$enso_group <- unname(country_group[as.character(d$country)])
    d[!is.na(d$enso_group) & d$enso_group %in% groups, , drop = FALSE]
  }
  pois <- function(d, track, scale) {
    d %>%
      dplyr::group_by(enso_group, month) %>%
      dplyr::summarise(
        n = sum(k),
        exposure = sum(e) / .SEAS_DAYS_PER_MONTH,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        track = track,
        rate = scale * n / exposure,
        lo = scale * seasonality_pois_ci(n, exposure, interval)$lo,
        hi = scale * seasonality_pois_ci(n, exposure, interval)$hi
      )
  }
  occ <- keep(surv) %>%
    dplyr::transmute(
      enso_group,
      month = lubridate::month(as.Date(t_start)),
      k = as.integer(as.logical(outbreak) %in% TRUE),
      e = as.numeric(exposure_days)
    ) %>%
    pois("occurrence", 1000)
  ces <- keep(dur) %>%
    dplyr::transmute(
      enso_group,
      month = lubridate::month(as.Date(month)),
      k = as.integer(ended),
      e = as.numeric(exposure)
    ) %>%
    pois("cessation", 100)
  cas <- keep(mon) %>%
    dplyr::filter(!(as.logical(cases_gap) %in% TRUE)) %>%
    dplyr::transmute(
      enso_group,
      month = lubridate::month(as.Date(month)),
      cases = as.numeric(cases_month)
    ) %>%
    dplyr::group_by(enso_group, month) %>%
    dplyr::summarise(
      rate = stats::median(cases),
      lo = stats::quantile(cases, 0.25, names = FALSE),
      hi = stats::quantile(cases, 0.75, names = FALSE),
      n = dplyr::n(),
      exposure = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(track = "cases")
  out <- dplyr::bind_rows(occ, ces, cas)
  # every group x track x month present (a month without rows: rate NA, n 0)
  out <- tidyr::expand_grid(
    enso_group = groups,
    track = c("occurrence", "cessation", "cases"),
    month = 1:12
  ) %>%
    dplyr::left_join(out, by = c("enso_group", "track", "month")) %>%
    dplyr::mutate(
      n = dplyr::coalesce(as.numeric(n), 0),
      exposure = dplyr::coalesce(exposure, 0),
      sparse = n < min_n
    )
  dplyr::select(
    out,
    enso_group,
    track,
    month,
    rate,
    lo,
    hi,
    n,
    exposure,
    sparse
  )
}

#' The monthly model's calendar seasonality from its harmonic coefficient
#' draws (pure): exp(h - mean_m h) per draw, h(m) = b1 sin(2 pi m / 12) +
#' b2 cos(2 pi m / 12) + b3 sin(4 pi m / 12) + b4 cos(4 pi m / 12)
#' (mon_data.R's harmonic2 columns on month_num = calendar month). Also the
#' duration model's, under duration.seasonal_term harmonic2 (dur_data.R builds
#' the same columns): a relative cessation hazard.
#' @param draws matrix draws x 4, columns b_sin_month, b_cos_month, b_sin2_month, b_cos2_month
#' @return tibble(month, rel, lo, hi)
seasonality_mon_harmonic <- function(draws, level = 0.95) {
  m <- 1:12
  X <- cbind(
    sin(2 * pi * m / 12),
    cos(2 * pi * m / 12),
    sin(4 * pi * m / 12),
    cos(4 * pi * m / 12)
  )
  h <- draws[,
    c("b_sin_month", "b_cos_month", "b_sin2_month", "b_cos2_month"),
    drop = FALSE
  ] %*%
    t(X)
  r <- exp(sweep(h, 1, rowMeans(h)))
  a <- (1 - level) / 2
  tibble::tibble(
    month = m,
    rel = colMeans(r),
    lo = apply(r, 2, stats::quantile, a, names = FALSE),
    hi = apply(r, 2, stats::quantile, 1 - a, names = FALSE)
  )
}

#' Peak month of each group x track x source and its position against the
#' group's rainy season (pure).
#' @param rates seasonality_track_rates()
#' @param fitted tibble(enso_group, track, month, rel)
#' @param rain tibble(enso_group, month, in_season)
#' @return tibble(enso_group, track, source ("observed" | "fitted"), peak_month,
#'   in_rainy, lag_from_rain, rainy_months: the rainy season(s) as runs,
#'   "Nov-Apr" or "Mar-May, Oct-Nov")
seasonality_peaks <- function(rates, fitted, rain) {
  pk <- dplyr::bind_rows(
    rates %>%
      dplyr::filter(!sparse, is.finite(rate)) %>%
      dplyr::mutate(source = "observed", v = rate),
    fitted %>% dplyr::mutate(source = "fitted", v = rel)
  ) %>%
    dplyr::group_by(enso_group, track, source) %>%
    dplyr::slice_max(v, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(enso_group, track, source, peak_month = month)
  runs <- lapply(split(rain, rain$enso_group), function(r) {
    on <- sort(r$month[r$in_season])
    if (!length(on)) {
      return(list(months = on, centres = numeric(0), label = "none"))
    }
    rr <- circular_month_runs(seq_len(12L) %in% on)
    list(
      months = on,
      centres = vapply(
        rr,
        function(x) ((x[1] - 1) + (length(x) - 1) / 2) %% 12 + 1,
        numeric(1)
      ),
      # each run in season order, e.g. "Nov-Apr" or "Mar-May, Nov-Dec"
      label = paste(
        vapply(
          rr,
          function(x) {
            if (length(x) == 1L) {
              month.abb[x]
            } else {
              paste(month.abb[x[1]], month.abb[x[length(x)]], sep = "-")
            }
          },
          ""
        ),
        collapse = ", "
      )
    )
  })
  lag_of <- function(g, p) {
    cs <- runs[[g]]$centres
    if (!length(cs) || is.na(p)) {
      return(NA_real_)
    }
    d <- ((p - cs + 6) %% 12) - 6
    d[which.min(abs(d))]
  }
  pk %>%
    dplyr::mutate(
      in_rainy = mapply(
        function(g, p) p %in% runs[[g]]$months,
        enso_group,
        peak_month
      ),
      lag_from_rain = mapply(lag_of, enso_group, peak_month),
      rainy_months = vapply(enso_group, function(g) runs[[g]]$label, "")
    )
}

#' Each group's deployed monthly-case fit file (read_deployed_lag(), kind
#' "mon"); an error when one is missing.
seasonality_mon_fit_files <- function(ctx, groups) {
  grouping <- as.character(
    ctx$MAN$grouping %||% sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  )
  stats::setNames(
    vapply(
      groups,
      function(g) {
        dl <- read_deployed_lag(
          ctx$ENSO_BEST_LAG_DIR,
          g,
          grouping,
          "mon",
          episode_key = episode_key_of(ctx$opt)
        )
        if (is.null(dl)) {
          stop(
            sprintf(
              "[09] seasonality: no deployed monthly-case fit for group %s in %s (grouping %s).",
              g,
              ctx$ENSO_BEST_LAG_DIR,
              grouping
            ),
            call. = FALSE
          )
        }
        dl$file
      },
      character(1)
    ),
    groups
  )
}

#' Each group's deployed duration fit file (read_deployed_lag(), kind "dur"),
#' NA where a group has none: the cessation curve is optional, drawn only for
#' the fits that carry duration.seasonal_term's harmonics.
seasonality_dur_fit_files <- function(ctx, groups) {
  grouping <- as.character(
    ctx$MAN$grouping %||% sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  )
  stats::setNames(
    vapply(
      groups,
      function(g) {
        dl <- read_deployed_lag(
          ctx$ENSO_BEST_LAG_DIR,
          g,
          grouping,
          "dur",
          episode_key = episode_key_of(ctx$opt)
        )
        if (is.null(dl)) NA_character_ else dl$file
      },
      character(1)
    ),
    groups
  )
}

#' seasonality v2 (v1 + the occurrence band).
#' @return list(rain, rates, fitted, peaks, meta). `fitted` holds a cessation
#'   curve (with its 95% band) only for groups whose deployed duration fit
#'   carries the harmonic2 terms (duration.seasonal_term); meta$dur_fits then
#'   names those fits, and is absent when no fit carries them.
derive_seasonality <- function(ctx) {
  d <- ctx$SEAS_DERIV
  groups <- as.character(results_group_meta(ctx)$order)
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  cg <- stats::setNames(as.character(reg$group), as.character(reg$iso3))
  rain_art <- ensure_rain_enso(ctx)
  occ_art <- ensure_occ_seasonality(ctx)
  rain <- tibble::as_tibble(rain_art$month) %>%
    dplyr::filter(weighting == "pop", group %in% groups) %>%
    dplyr::transmute(
      enso_group = group,
      month = as.integer(mo),
      mm,
      in_season,
      share
    )
  miss <- setdiff(groups, rain$enso_group)
  if (length(miss)) {
    stop(
      "[09] seasonality: rain_enso has no pop-weighted months for group(s) ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  rates <- seasonality_track_rates(
    readRDS(require_input(ctx$SURV_PATH, "survival data")),
    readRDS(require_input(ctx$DUR_DATA_PATH, "duration data")),
    readRDS(require_input(ctx$MON_DATA_PATH, "monthly case data")),
    cg,
    groups,
    admin_level = as.character(d$admin_level),
    min_n = as.integer(d$min_n),
    interval = as.numeric(d$interval)
  )
  files <- seasonality_mon_fit_files(ctx, groups)
  vars <- c("b_sin_month", "b_cos_month", "b_sin2_month", "b_cos2_month")
  mon_fit <- dplyr::bind_rows(lapply(groups, function(g) {
    fit <- readRDS(files[[g]])
    if (!all(vars %in% brms::variables(fit))) {
      stop(
        "[09] seasonality: the monthly-case fit of group ",
        g,
        " has no harmonic2 seasonal terms (",
        basename(files[[g]]),
        ").",
        call. = FALSE
      )
    }
    dplyr::mutate(
      seasonality_mon_harmonic(
        as.matrix(posterior::as_draws_matrix(fit, variable = vars)),
        as.numeric(d$interval)
      ),
      enso_group = g,
      track = "cases",
      .before = 1
    )
  }))
  # Cessation (duration.seasonal_term harmonic2): only the fits that carry the
  # four harmonics; a fit without them contributes no row (the deployed fits).
  dur_files <- seasonality_dur_fit_files(ctx, groups)
  dur_fit <- dplyr::bind_rows(lapply(groups[!is.na(dur_files)], function(g) {
    fit <- readRDS(dur_files[[g]])
    if (!all(vars %in% brms::variables(fit))) {
      return(NULL)
    }
    dplyr::mutate(
      seasonality_mon_harmonic(
        as.matrix(posterior::as_draws_matrix(fit, variable = vars)),
        as.numeric(d$interval)
      ),
      enso_group = g,
      track = "cessation",
      .before = 1
    )
  }))
  # v2: the occurrence band at `interval`, from occ_seasonality's per-draw curves.
  a <- (1 - as.numeric(d$interval)) / 2
  occ_band <- dplyr::bind_rows(lapply(groups, function(g) {
    m <- occ_art$draws[[g]]
    if (is.null(m)) {
      stop(
        "[09] seasonality: occ_seasonality has no draws for group ",
        g,
        call. = FALSE
      )
    }
    tibble::tibble(
      enso_group = g,
      month = seq_len(12L),
      lo = apply(m, 2, stats::quantile, a, names = FALSE),
      hi = apply(m, 2, stats::quantile, 1 - a, names = FALSE)
    )
  }))
  occ_fit <- tibble::as_tibble(occ_art$seasonality) %>%
    dplyr::filter(enso_group %in% groups) %>%
    dplyr::transmute(
      enso_group,
      track = "occurrence",
      month = as.integer(month),
      rel
    ) %>%
    dplyr::left_join(occ_band, by = c("enso_group", "month"))
  obs_mean <- rates %>%
    dplyr::group_by(enso_group, track) %>%
    dplyr::summarise(obs_mean = mean(rate, na.rm = TRUE), .groups = "drop")
  fitted <- dplyr::bind_rows(occ_fit, dur_fit, mon_fit) %>%
    dplyr::left_join(obs_mean, by = c("enso_group", "track")) %>%
    dplyr::mutate(
      scaled = rel * obs_mean,
      scaled_lo = lo * obs_mean,
      scaled_hi = hi * obs_mean
    ) %>%
    dplyr::select(-obs_mean)
  peaks <- seasonality_peaks(
    rates,
    dplyr::select(fitted, enso_group, track, month, rel),
    dplyr::select(rain, enso_group, month, in_season)
  )
  message(sprintf(
    "[09] seasonality: %d groups; observed peaks in the rainy season: %s",
    length(groups),
    paste(
      sprintf(
        "%s %d/%d",
        c("occurrence", "cessation", "cases"),
        vapply(
          c("occurrence", "cessation", "cases"),
          function(t) {
            sum(peaks$in_rainy[peaks$source == "observed" & peaks$track == t])
          },
          0L
        ),
        length(groups)
      ),
      collapse = ", "
    )
  ))
  meta <- list(
    groups = groups,
    mon_fits = tibble::tibble(enso_group = groups, file = unname(files)),
    admin_level = as.character(d$admin_level),
    min_n = as.integer(d$min_n),
    interval = as.numeric(d$interval),
    units = c(
      occurrence = "onsets per 1,000 at-risk unit-months",
      cessation = "episode ends per 100 active episode-months",
      cases = "median suspected cases per active episode-month"
    )
  )
  if (nrow(dur_fit)) {
    g_dur <- unique(dur_fit$enso_group)
    meta$dur_fits <- tibble::tibble(
      enso_group = g_dur,
      file = unname(dur_files[g_dur])
    )
  }
  list(rain = rain, rates = rates, fitted = fitted, peaks = peaks, meta = meta)
}

#' Freshness-guarded seasonality artefact (supplement fig. SE).
ensure_seasonality <- function(ctx) {
  d <- ctx$SEAS_DERIV
  need <- c("admin_level", "min_n", "interval")
  if (is.null(d) || !all(need %in% names(d))) {
    stop(
      "[09] seasonality: results.derivation.seasonality must set ",
      paste(need, collapse = ", "),
      call. = FALSE
    )
  }
  groups <- as.character(results_group_meta(ctx)$order)
  inputs <- results_inputs(
    files = c(
      ctx$SURV_PATH,
      ctx$DUR_DATA_PATH,
      ctx$MON_DATA_PATH,
      ctx$GROUPING_PATH,
      ctx$RAIN_PATH,
      ctx$GID_POP_PATH,
      occ_seasonality_fit_files(ctx, groups),
      seasonality_mon_fit_files(ctx, groups),
      seasonality_dur_fit_files(ctx, groups)
    )
  ) # NA (no dur fit) dropped by results_inputs()
  ensure_artefact(ctx, "seasonality", inputs, function() {
    derive_seasonality(ctx)
  })
}

# --- Figure 1c's composite ONI row ------------------------------------------
#
# oni_composite (v1, unkeyed): the observed ONI over a 2-year cycle (x = 1..24:
# the months of Y and Y+1) per ENSO episode kind, averaged over the episodes
# of the oni artefact's `phases` (El Nino, La Nina runs). Each episode is
# aligned on its peak (El Nino: max ONI; La Nina: min ONI) inside the run:
# Y = the peak's year if the peak falls in Jul-Dec, else the year before. A
# cycle whose 24 months are not all observed is left out; two episodes of one
# kind in the same cycle count once.

#' Composite ONI cycle (pure).
#' @param obs tibble(date, enso_index): the observed ONI
#' @param phases tibble(phase, TL, TR): the episode runs
#' @return tibble(kind ("elnino" | "nina"), x, mean, n = cycles)
oni_composite_cycle <- function(obs, phases) {
  obs <- dplyr::mutate(
    obs,
    y = lubridate::year(date),
    m = lubridate::month(date)
  )
  ep <- dplyr::filter(phases, phase %in% c("elnino", "nina"))
  rows <- lapply(seq_len(nrow(ep)), function(i) {
    w <- dplyr::filter(obs, date >= ep$TL[i], date <= ep$TR[i])
    if (!nrow(w)) {
      return(NULL)
    }
    pk <- if (ep$phase[i] == "elnino") {
      w[which.max(w$enso_index), ]
    } else {
      w[which.min(w$enso_index), ]
    }
    y0 <- if (pk$m >= 7) pk$y else pk$y - 1L
    cyc <- dplyr::filter(obs, y %in% c(y0, y0 + 1L)) %>%
      dplyr::mutate(
        x = as.integer((y - y0) * 12L + m),
        kind = ep$phase[i],
        y0 = y0
      )
    if (nrow(cyc) != 24L) {
      return(NULL)
    }
    cyc
  })
  dplyr::bind_rows(rows) %>%
    dplyr::distinct(kind, y0, x, .keep_all = TRUE) %>%
    dplyr::group_by(kind, x) %>%
    dplyr::summarise(mean = mean(enso_index), n = dplyr::n(), .groups = "drop")
}

#' oni_composite v1.
#' @return list(cycle = oni_composite_cycle(), episodes = tibble(kind, n))
derive_oni_composite <- function(ctx) {
  oni <- ensure_oni(ctx, ctx$SETS)
  cyc <- oni_composite_cycle(oni$obs_full, oni$phases)
  n <- dplyr::distinct(cyc, kind, n)
  message(sprintf(
    "[09] oni_composite: %s",
    paste(sprintf("%s %d cycles", n$kind, n$n), collapse = ", ")
  ))
  list(cycle = cyc, episodes = n)
}

#' Freshness-guarded composite-ONI artefact: its input is the oni artefact.
ensure_oni_composite <- function(ctx) {
  ensure_oni(ctx, ctx$SETS)
  inputs <- results_inputs(
    files = make_results_run_file(
      ctx$opt,
      "results_oni",
      sets = ctx$SETS,
      schema = RESULTS_SCHEMA[["oni"]],
      config_path = ctx$CONFIG_PATH
    )
  )
  ensure_artefact(ctx, "oni_composite", inputs, function() {
    derive_oni_composite(ctx)
  })
}

# --- Supplement fig. WR: observed vs modelled ENSO ratio ---------
#
# report_who_enso_ratio (v1, unkeyed). For the countries of
# results.derivation.who_enso_ratio.groups (West only), two ratios of the target year,
# the year of the attribution
# window (2026):
#   modelled  ratio of means of the smoothing (the factual) to the ENSO
#             counterfactual of each analogue set over their common months, the
#             attribution window's stitch months (counterfactual launch to
#             factual launch, e.g. May-August 2026), cases summed over the
#             window per replicate. Spread (model_lo .. model_hi): each smoothing
#             replicate's total over the mean of the
#             counterfactual replicates launched from it (handoff_map.csv), at
#             `interval` and 50%; the paired bootstrap of the ratio of means
#             over smoothing replicates is kept as rom_lo / rom_hi.
#   observed  WHO suspected cases of the target year over the mean of the
#             ENSO-neutral years (each of `year_definitions`;
#             enso.R::enso_year_states(), the lag_occ window at the group's
#             deployed occurrence lag from group_labels), population-
#             standardised to the target year (WPP), two estimators:
#               same_window  the window's calendar months in the WHO monthly
#                            series (03_prepare_who_monthly.R; complete months
#                            only, so dashboard years only);
#               ytd_annual   the target year to the WHO cut-off against full
#                            years (who_annual.R::build_who_annual_cases(): annual
#                            report + dashboard weekly, absence NA, censored
#                            target-year series NA) -- a lower bound.
#             Interval: bootstrap over the neutral years (none with one year).
# Both ratios add `pseudo_count` cases to numerator and denominator. A group
# row pools the countries with an observed ratio in its estimator x definition,
# over neutral years where every one of them has a value; its modelled ratio
# is over the same countries. Every analogue set of ctx$SETS is derived; 10
# draws results.style.who_enso_ratio's.

#' Smoothing window totals per replicate and country, zero-filled (pure).
#' @param events_mon smoothing events (rep_id, gid, t_start, cases_month[, country])
#' @param rep_ids every smoothing replicate
#' @return tibble(rep_id, country, cases)
who_enso_smoothing_totals <- function(events_mon, rep_ids, countries, months) {
  ev <- tibble::as_tibble(events_mon)
  if (!"country" %in% names(ev)) {
    ev$country <- substr(as.character(ev$gid), 1L, 3L)
  }
  tot <- ev %>%
    dplyr::filter(
      as.Date(t_start) %in% as.Date(months),
      country %in% countries
    ) %>%
    dplyr::summarise(
      cases = sum(as.numeric(cases_month)),
      .by = c(rep_id, country)
    )
  tidyr::expand_grid(
    rep_id = sort(unique(as.integer(rep_ids))),
    country = sort(unique(countries))
  ) %>%
    dplyr::left_join(tot, by = c("rep_id", "country")) %>%
    dplyr::mutate(cases = dplyr::coalesce(cases, 0))
}

#' Counterfactual window totals averaged over the replicates launched from each
#' smoothing replicate (pure). Every replicate must cover every window month and
#' every smoothing replicate must have the same number of children (else the
#' plain counterfactual mean is not the mean of the parent means).
#' @param country_month load_sim_run_monthly()$country_month (replicate_id, month, n_cases, country)
#' @param handoff read_handoff_map() (replicate_id, smoothing_rep_id)
#' @return tibble(smoothing_rep_id, country, cases, n_children)
who_enso_cf_parent_means <- function(
  country_month,
  handoff,
  countries,
  months
) {
  if (is.null(handoff) || !nrow(handoff)) {
    stop(
      "[09] report_who_enso_ratio: a counterfactual run has no handoff_map.csv (runs before 2026-10-20 ",
      "cannot be paired with the smoothing).",
      call. = FALSE
    )
  }
  cm <- tibble::as_tibble(country_month) %>%
    dplyr::filter(
      as.Date(month) %in% as.Date(months),
      country %in% countries
    ) %>%
    dplyr::summarise(
      cases = sum(as.numeric(n_cases)),
      n_months = dplyr::n_distinct(month),
      .by = c(replicate_id, country)
    )
  if (any(cm$n_months != length(months))) {
    stop(
      "[09] report_who_enso_ratio: a counterfactual replicate does not cover every window month.",
      call. = FALSE
    )
  }
  hm <- dplyr::distinct(
    tibble::as_tibble(handoff),
    replicate_id = as.integer(replicate_id),
    smoothing_rep_id = as.integer(smoothing_rep_id)
  )
  cm <- dplyr::left_join(cm, hm, by = "replicate_id")
  if (anyNA(cm$smoothing_rep_id)) {
    stop(
      "[09] report_who_enso_ratio: counterfactual replicates missing from handoff_map.csv.",
      call. = FALSE
    )
  }
  out <- cm %>%
    dplyr::summarise(
      cases = mean(cases),
      n_children = dplyr::n(),
      .by = c(smoothing_rep_id, country)
    )
  if (length(unique(out$n_children)) != 1L) {
    stop(
      "[09] report_who_enso_ratio: unequal counterfactual replicates per smoothing replicate (",
      paste(sort(unique(out$n_children)), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  out
}

#' Modelled ratio of one country set (pure): ratio of means, per-smoothing-
#' replicate spread, paired bootstrap of the ratio of means.
#' @param S who_enso_smoothing_totals(); C who_enso_cf_parent_means() of one set
#' @return one-row tibble
who_enso_model_ratio <- function(
  S,
  C,
  countries,
  pseudo = 1,
  interval = 0.95,
  n_boot = 2000L,
  seed = 1L
) {
  sv <- S %>%
    dplyr::filter(country %in% countries) %>%
    dplyr::summarise(x = sum(cases), .by = rep_id) %>%
    dplyr::arrange(rep_id)
  cs <- C %>%
    dplyr::filter(country %in% countries) %>%
    dplyr::summarise(x = sum(cases), .by = smoothing_rep_id)
  cv <- cs$x[match(sv$rep_id, cs$smoothing_rep_id)]
  if (anyNA(cv)) {
    stop(
      "[09] report_who_enso_ratio: smoothing replicate(s) with no counterfactual child.",
      call. = FALSE
    )
  }
  a <- (1 - interval) / 2
  rr <- (sv$x + pseudo) / (cv + pseudo)
  n <- length(sv$x)
  rb <- .with_seed(seed, {
    idx <- matrix(sample.int(n, n_boot * n, replace = TRUE), nrow = n_boot)
    (rowMeans(matrix(sv$x[idx], nrow = n_boot)) + pseudo) /
      (rowMeans(matrix(cv[idx], nrow = n_boot)) + pseudo)
  })
  q <- function(x, p) unname(stats::quantile(x, p))
  tibble::tibble(
    mean_smoothing = mean(sv$x),
    mean_counterfactual = mean(cv),
    model_ratio = (mean(sv$x) + pseudo) / (mean(cv) + pseudo),
    model_rep_median = stats::median(rr),
    model_lo = q(rr, a),
    model_q25 = q(rr, 0.25),
    model_q75 = q(rr, 0.75),
    model_hi = q(rr, 1 - a),
    prob_rep_ratio_gt1 = mean(rr > 1),
    rom_lo = q(rb, a),
    rom_hi = q(rb, 1 - a),
    n_smoothing_reps = n
  )
}

#' Observed ratio of one country set (pure): the target year over the mean of
#' the neutral years in which every country has a value, bootstrap over those
#' years, and the ratios to the single highest / lowest year.
#' @param obs (country, year, cases_std)
#' @return one-row tibble
who_enso_observed_ratio <- function(
  obs,
  countries,
  neutral_years,
  target,
  pseudo = 1,
  interval = 0.95,
  n_boot = 2000L,
  seed = 1L
) {
  o <- obs %>%
    dplyr::filter(country %in% countries, !is.na(cases_std)) %>%
    dplyr::summarise(x = sum(cases_std), n = dplyr::n(), .by = year)
  full <- o$n == length(unique(countries))
  yrs <- sort(o$year[full & o$year %in% neutral_years & o$year < target])
  xt <- o$x[full & o$year == target]
  out <- tibble::tibble(
    obs_target_std = NA_real_,
    obs_neutral_mean_std = NA_real_,
    obs_ratio = NA_real_,
    obs_lo = NA_real_,
    obs_hi = NA_real_,
    obs_single_year_min = NA_real_,
    obs_single_year_max = NA_real_,
    n_neutral = length(yrs),
    neutral_years = paste(yrs, collapse = " ")
  )
  if (!length(xt) || !length(yrs)) {
    return(out)
  }
  xn <- o$x[match(yrs, o$year)]
  a <- (1 - interval) / 2
  out$obs_target_std <- xt
  out$obs_neutral_mean_std <- mean(xn)
  out$obs_ratio <- (xt + pseudo) / (mean(xn) + pseudo)
  out$obs_single_year_min <- (xt + pseudo) / (max(xn) + pseudo)
  out$obs_single_year_max <- (xt + pseudo) / (min(xn) + pseudo)
  if (length(xn) > 1L) {
    rb <- .with_seed(
      seed,
      (xt + pseudo) /
        (rowMeans(matrix(
          xn[sample.int(length(xn), n_boot * length(xn), replace = TRUE)],
          nrow = n_boot
        )) +
          pseudo)
    )
    out$obs_lo <- unname(stats::quantile(rb, a))
    out$obs_hi <- unname(stats::quantile(rb, 1 - a))
  }
  out
}

#' Cases standardised to the target year's population (pure).
#' @param d (country, year, cases); pop (country, year, pop)
who_enso_standardise <- function(d, pop, target) {
  pt <- pop %>%
    dplyr::filter(year == target) %>%
    dplyr::select(country, pop_target = pop)
  out <- d %>%
    dplyr::left_join(pop, by = c("country", "year")) %>%
    dplyr::left_join(pt, by = "country") %>%
    dplyr::mutate(cases_std = cases * pop_target / pop)
  bad <- out$country[
    !is.na(out$cases) & (is.na(out$pop) | is.na(out$pop_target))
  ]
  if (length(bad)) {
    stop(
      "[09] report_who_enso_ratio: no WPP population for ",
      paste(unique(bad), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

#' The model, observed and comparison tables (pure).
#' @param annual build_who_annual_cases()$annual; window_tot who_window_totals()
#' @param pop (country, year, pop); states named list group -> enso_year_states()
#' @param country_group named character ISO3 -> group; model_countries the simulated ones
#' @param S smoothing totals; C_by_set named list set -> parent means
#' @param cfg list(groups, year_definitions, pseudo_count, interval, n_boot, seed)
#' @return list(model, comparison, observed = list(ytd_annual, same_window))
who_enso_ratio_tables <- function(
  annual,
  window_tot,
  pop,
  states,
  country_group,
  model_countries,
  S,
  C_by_set,
  target,
  cfg
) {
  pr <- function(f, ...) {
    f(
      ...,
      pseudo = cfg$pseudo_count,
      interval = cfg$interval,
      n_boot = cfg$n_boot,
      seed = cfg$seed
    )
  }
  obs <- list(
    ytd_annual = who_enso_standardise(
      dplyr::filter(annual, !is.na(cases)) %>%
        dplyr::select(country, year, cases),
      pop,
      target
    ),
    same_window = who_enso_standardise(window_tot, pop, target)
  )
  model <- dplyr::bind_rows(lapply(cfg$groups, function(g) {
    mc <- sort(intersect(
      names(country_group)[country_group == g],
      model_countries
    ))
    dplyr::bind_rows(lapply(names(C_by_set), function(s) {
      dplyr::bind_rows(
        lapply(mc, function(cc) {
          pr(who_enso_model_ratio, S, C_by_set[[s]], cc) %>%
            dplyr::mutate(scope = cc)
        }),
        pr(who_enso_model_ratio, S, C_by_set[[s]], mc) %>%
          dplyr::mutate(scope = g)
      ) %>%
        dplyr::mutate(
          enso_group = g,
          analog_set = s,
          countries = c(mc, paste(mc, collapse = " ")),
          level = c(rep("country", length(mc)), "group")
        )
    }))
  })) %>%
    dplyr::relocate(enso_group, analog_set, scope, level, countries)
  comparison <- dplyr::bind_rows(lapply(cfg$groups, function(g) {
    mc <- sort(intersect(
      names(country_group)[country_group == g],
      model_countries
    ))
    st <- states[[g]]
    dplyr::bind_rows(lapply(cfg$year_definitions, function(d) {
      ny <- st$year[st$year < target & st[[paste0("state_", d)]] %in% "neutral"]
      dplyr::bind_rows(lapply(names(obs), function(e) {
        per <- dplyr::bind_rows(lapply(mc, function(cc) {
          pr(who_enso_observed_ratio, obs[[e]], cc, ny, target) %>%
            dplyr::mutate(scope = cc, countries = cc)
        }))
        with_obs <- per$scope[is.finite(per$obs_ratio)]
        grp <- if (length(with_obs)) {
          pr(who_enso_observed_ratio, obs[[e]], with_obs, ny, target) %>%
            dplyr::mutate(
              scope = g,
              countries = paste(with_obs, collapse = " ")
            )
        }
        dplyr::bind_rows(per, grp) %>%
          dplyr::mutate(
            level = ifelse(scope == g, "group", "country"),
            enso_group = g,
            year_definition = d,
            estimator = e,
            state_target = st[[paste0("state_", d)]][st$year == target] %||%
              NA_character_
          )
      }))
    }))
  }))
  comparison <- tidyr::crossing(comparison, analog_set = names(C_by_set)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      model = list(pr(
        who_enso_model_ratio,
        S,
        C_by_set[[analog_set]],
        strsplit(countries, " ")[[1]]
      ))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(model) %>%
    dplyr::relocate(
      enso_group,
      year_definition,
      estimator,
      analog_set,
      scope,
      level,
      countries
    )
  list(model = model, comparison = comparison, observed = obs)
}

#' The step's config, validated (stops on a missing key or an unknown value).
.who_enso_ratio_cfg <- function(ctx) {
  d <- ctx$WER_DERIV
  need <- c(
    "groups",
    "year_definitions",
    "first_year",
    "censor_weeks",
    "pseudo_count",
    "interval",
    "n_boot",
    "seed"
  )
  if (is.null(d) || !all(need %in% names(d))) {
    stop(
      "[09] report_who_enso_ratio: results.derivation.who_enso_ratio must set ",
      paste(need, collapse = ", "),
      call. = FALSE
    )
  }
  groups <- setdiff(as.character(unlist(d$groups)), "")
  defs <- as.character(unlist(d$year_definitions))
  bad <- setdiff(defs, names(ENSO_YEAR_DEFINITIONS))
  if (length(bad) || !length(defs)) {
    stop(
      "[09] report_who_enso_ratio: year_definitions must be among ",
      paste(names(ENSO_YEAR_DEFINITIONS), collapse = ", "),
      "; got ",
      paste(defs, collapse = ", "),
      call. = FALSE
    )
  }
  known <- as.character(results_group_meta(ctx)$order)
  if (!length(groups) || length(setdiff(groups, known))) {
    stop(
      "[09] report_who_enso_ratio: groups must be fit groups of the registry (",
      paste(known, collapse = ", "),
      "); got ",
      paste(groups, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    groups = groups,
    year_definitions = defs,
    first_year = as.integer(d$first_year),
    censor_weeks = as.integer(d$censor_weeks),
    pseudo_count = as.numeric(d$pseudo_count),
    interval = as.numeric(d$interval),
    n_boot = as.integer(d$n_boot),
    seed = as.integer(d$seed)
  )
}

#' The counterfactual runs of the units holding `groups`, one row per unit x set.
.who_enso_ratio_runs <- function(ctx, groups, sets) {
  units <- manifest_units(ctx$MAN)
  units <- names(units)[vapply(
    units,
    function(u) any(groups %in% u$groups),
    logical(1)
  )]
  if (!length(units)) {
    stop(
      "[09] report_who_enso_ratio: no manifest unit holds ",
      paste(groups, collapse = ", "),
      call. = FALSE
    )
  }
  dplyr::bind_rows(lapply(sets, function(s) {
    runs <- manifest_runs(ctx$MAN, "counterfactual", s)
    dplyr::bind_rows(lapply(units, function(un) {
      if (is.null(runs[[un]]$config)) {
        stop(
          "[09] report_who_enso_ratio: unit '",
          un,
          "' has no counterfactual run for set '",
          s,
          "'.",
          call. = FALSE
        )
      }
      tibble::tibble(
        unit = un,
        analog_set = s,
        config = runs[[un]]$config,
        tag = yaml::read_yaml(runs[[un]]$config)$simulator$run_tag
      )
    }))
  }))
}

#' Whether report_who_enso_ratio can run on this config (the step's `enabled`):
#' a smoothing run, at least one group, stitched arms (results.derivation.
#' arm_alignment "stitch" and a counterfactual launch before the factual one),
#' and counterfactual runs that record their handoff (handoff_map.csv). A run whose
#' partitions are absent here counts as pairable, so a
#' missing output still fails loudly in the step instead of disabling it.
#' A context without a manifest (the registry tests) is judged on the config.
.who_enso_ratio_enabled <- function(ctx) {
  if (!isTRUE(nzchar(ctx$SMOOTHING_PARTS_DIR %||% ""))) {
    return(FALSE)
  }
  if (!identical(as.character(ctx$ARM_ALIGNMENT %||% "stitch"), "stitch")) {
    return(FALSE)
  }
  groups <- setdiff(as.character(unlist(ctx$WER_DERIV[["groups"]])), "")
  if (!length(groups)) {
    return(FALSE)
  }
  if (is.null(ctx$MAN)) {
    return(TRUE)
  }
  sets <- as.character(ctx$SETS)
  starts <- tryCatch(
    manifest_arm_starts(ctx$MAN, groups, sets),
    error = function(e) NULL
  )
  if (
    is.null(starts) ||
      !isTRUE(as.Date(starts$counterfactual) < as.Date(starts$factual))
  ) {
    return(FALSE)
  }
  runs <- tryCatch(
    .who_enso_ratio_runs(ctx, groups, sets),
    error = function(e) NULL
  )
  if (is.null(runs)) {
    return(FALSE)
  }
  for (t in unique(runs$tag)) {
    done <- dirname(partition_sentinels(file.path(ctx$SIM_ROOT, t, "parts")))
    if (length(done) && !any(file.exists(file.path(done, "handoff_map.csv")))) {
      return(FALSE)
    }
  }
  TRUE
}

#' report_who_enso_ratio v1.
#' @return list(model, comparison, observed, annual, window_totals, year_states, meta)
derive_report_who_enso_ratio <- function(ctx) {
  cfg <- .who_enso_ratio_cfg(ctx)
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  cg <- stats::setNames(as.character(reg$group), as.character(reg$iso3))
  countries <- sort(names(cg)[cg %in% cfg$groups])
  sets <- as.character(ctx$SETS)
  win <- results_arm_window(ctx, cfg$groups, sets)
  if (!isTRUE(win$stitched) || !length(win$stitch_months)) {
    stop(
      "[09] report_who_enso_ratio: the counterfactual and factual launch together, so the smoothing and the ",
      "counterfactual share no month (needs results.derivation.arm_alignment 'stitch' and an earlier ",
      "counterfactual launch).",
      call. = FALSE
    )
  }
  months <- as.Date(win$stitch_months)
  target <- unique(as.integer(format(months, "%Y")))
  if (length(target) != 1L) {
    stop(
      "[09] report_who_enso_ratio: the common months span two years (",
      format(min(months)),
      " to ",
      format(max(months)),
      ").",
      call. = FALSE
    )
  }
  S <- who_enso_smoothing_totals(
    win$smoothing_events$events_mon,
    win$smoothing_events$rep_ids,
    countries,
    months
  )
  runs <- .who_enso_ratio_runs(ctx, cfg$groups, sets)
  series <- as.character(ctx$RES$series %||% "latent")
  C_by_set <- stats::setNames(
    lapply(sets, function(s) {
      rs <- runs[runs$analog_set == s, ]
      dplyr::bind_rows(lapply(seq_len(nrow(rs)), function(i) {
        r <- load_sim_run_monthly(
          ctx$SIM_ROOT,
          rs$tag[i],
          min(months),
          max(months),
          series = series,
          tag = "[09]"
        )
        hm <- read_handoff_map(ctx$SIM_ROOT, rs$tag[i])
        if (!is.null(hm) && any(as.Date(hm$handoff_month) != min(months))) {
          stop(
            "[09] report_who_enso_ratio: run ",
            rs$tag[i],
            " was not handed off at ",
            format(min(months)),
            call. = FALSE
          )
        }
        who_enso_cf_parent_means(r$country_month, hm, countries, months)
      }))
    }),
    sets
  )
  model_countries <- sort(unique(unlist(lapply(C_by_set, `[[`, "country"))))
  rep_off <- setdiff(
    unique(unlist(lapply(C_by_set, `[[`, "smoothing_rep_id"))),
    unique(S$rep_id)
  )
  if (length(rep_off)) {
    stop(
      "[09] report_who_enso_ratio: counterfactual runs handed off from smoothing replicates absent from the ",
      "smoothing parts (",
      paste(utils::head(rep_off, 5L), collapse = ", "),
      ").",
      call. = FALSE
    )
  }

  wr <- build_who_annual_cases(
    readRDS(require_input(make_who_annual_file("report"), "WHO annual report")),
    readRDS(require_input(
      make_who_annual_file("dashboard_weekly"),
      "WHO dashboard weekly series"
    )),
    countries,
    cfg$first_year:target,
    target,
    cfg$censor_weeks
  )
  window_tot <- who_window_totals(
    readRDS(require_input(
      make_who_monthly_file(),
      "WHO monthly table (03_prepare_who_monthly.R)"
    )),
    countries,
    as.integer(format(months, "%m"))
  )
  wpp <- load_wpp_table(.pop_artefact_or_null(ctx$WPP_PATH))
  if (is.null(wpp)) {
    stop(
      "[09] report_who_enso_ratio: WPP national totals not found: ",
      ctx$WPP_PATH,
      call. = FALSE
    )
  }
  pop <- tibble::as_tibble(wpp) %>%
    dplyr::transmute(
      country = as.character(iso3),
      year = as.integer(year),
      pop = as.numeric(pop_total)
    )
  oni <- fetch_enso_index(
    source = "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  )
  lags <- ensure_group_labels(ctx)$groups
  lag_of <- stats::setNames(as.integer(lags$lag_occ), lags$enso_group)
  if (anyNA(lag_of[cfg$groups])) {
    stop(
      "[09] report_who_enso_ratio: no deployed occurrence lag for ",
      paste(cfg$groups[is.na(lag_of[cfg$groups])], collapse = ", "),
      call. = FALSE
    )
  }
  states <- stats::setNames(
    lapply(cfg$groups, function(g) {
      enso_year_states(oni, cfg$first_year:target, lag_of[[g]])
    }),
    cfg$groups
  )
  tabs <- who_enso_ratio_tables(
    wr$annual,
    window_tot,
    pop,
    states,
    cg[cg %in% cfg$groups],
    model_countries,
    S,
    C_by_set,
    target,
    cfg
  )
  n_children <- unique(unlist(lapply(C_by_set, `[[`, "n_children")))
  message(sprintf(
    "[09] report_who_enso_ratio: %s, %s; %d smoothing x %s counterfactual replicates per set; %s censored",
    paste(cfg$groups, collapse = ","),
    paste(format(range(months), "%Y-%m"), collapse = " to "),
    length(unique(S$rep_id)),
    paste(n_children, collapse = ","),
    if (length(wr$censored)) paste(wr$censored, collapse = ",") else "none"
  ))
  list(
    model = tabs$model,
    comparison = tabs$comparison,
    observed = tabs$observed,
    annual = who_enso_standardise(wr$annual, pop, target),
    window_totals = window_tot,
    year_states = dplyr::bind_rows(lapply(cfg$groups, function(g) {
      dplyr::mutate(
        states[[g]],
        enso_group = g,
        lag_occ = lag_of[[g]],
        .before = 1
      )
    })),
    meta = c(
      cfg,
      list(
        countries = countries,
        model_countries = model_countries,
        censored = wr$censored,
        cutoff = wr$cutoff,
        window = months,
        target_year = target,
        analog_sets = sets,
        lag_occ = lag_of[cfg$groups],
        n_smoothing_reps = length(unique(S$rep_id)),
        n_children_per_smoothing_rep = n_children,
        series = series,
        runs = runs,
        smoothing_parts_dir = ctx$SMOOTHING_PARTS_DIR,
        enso_rule = enso_phase_definition("threshold", "year_mean_oni")$label,
        year_definition_labels = ENSO_YEAR_DEFINITIONS
      )
    )
  )
}

#' Freshness-guarded report_who_enso_ratio artefact (supplement fig. WR).
ensure_report_who_enso_ratio <- function(ctx) {
  cfg <- .who_enso_ratio_cfg(ctx)
  runs <- .who_enso_ratio_runs(ctx, cfg$groups, as.character(ctx$SETS))
  wpp <- as.character(ctx$WPP_PATH %||% "")
  inputs <- results_inputs(
    files = c(
      make_who_annual_file("report"),
      make_who_annual_file("dashboard_weekly"),
      make_who_monthly_file(),
      if (nzchar(wpp) && file.exists(wpp)) wpp,
      .oni_cache_file(ctx),
      ctx$GROUPING_PATH,
      ctx$MANIFEST_PATH,
      unique(runs$config),
      # the deployed occurrence lags of lag_occ (group_labels reads them from here)
      if (
        nzchar(ctx$ENSO_BEST_LAG_MANIFEST %||% "") &&
          file.exists(ctx$ENSO_BEST_LAG_MANIFEST)
      ) {
        ctx$ENSO_BEST_LAG_MANIFEST
      }
    ),
    globs = c(
      smoothing_parts_glob(ctx),
      lapply(unique(runs$tag), function(t) {
        list(
          dir = file.path(ctx$SIM_ROOT, t, "parts"),
          pattern = "_SUCCESS$",
          recursive = TRUE
        )
      })
    )
  )
  ensure_artefact(ctx, "report_who_enso_ratio", inputs, function() {
    derive_report_who_enso_ratio(ctx)
  })
}

# ===========================================================================
# Section 2: forecast / counterfactual aggregate derivation (Figures 3 & 4)
# ===========================================================================

#' This run's binding of utils/results_core.R::load_sim_run_monthly(): `cases` is
#' the mon track's `cases_month` (one draw per ACTIVE month; the size track's
#' onset-lumped `Z` is not used). Requires the
#' run to have been fit with grouping.monthly_enable: true.
#' `window_start`/`window_end` bound the ADM1-level table only (Figure 4's
#' attribution window).
#' `arm` / `win` (results_arm_window()): the factual ("forecast") arm is
#' stitched from the smoothing when the window starts before its launch
#' (utils/results_core.R::load_arm_run()); without `win` the read is plain.
load_run <- function(
  ctx,
  run_tag,
  window_start,
  window_end,
  arm = "forecast",
  win = NULL,
  ...
) {
  load_arm_run(
    ctx$SIM_ROOT,
    run_tag,
    arm,
    window_start,
    window_end,
    arm_starts = win$arm_starts,
    smoothing_events = win$smoothing_events,
    factual_source = win$factual_source %||% "simulation",
    series = ctx$RES$series,
    tag = "[09]",
    ...
  )
}

#' results_arm_window
#'
#' The forecast / counterfactual window of Figures 3-4 and the forecast report:
#' each arm's launch (utils/sim_io.R::manifest_arm_starts()), the attribution
#' start under results.derivation.arm_alignment, the ratio label, and -- only
#' when the window starts before the factual launch -- the smoothing events the
#' factual prefix is stitched from (read once here, shared by every group).
results_arm_window <- function(ctx, groups, sets, read_events = TRUE) {
  arm_starts <- manifest_arm_starts(ctx$MAN, groups, sets)
  attrib_start <- resolve_attrib_start(
    ctx$MAN,
    groups,
    sets,
    alignment = ctx$ARM_ALIGNMENT
  )
  stitched <- attrib_start < arm_starts$factual
  list(
    attrib_start = attrib_start,
    factual_start = arm_starts$factual,
    arm_starts = arm_starts,
    stitched = stitched,
    stitch_months = if (stitched) {
      m <- seq(attrib_start, arm_starts$factual, by = "month")
      m[m < arm_starts$factual]
    } else {
      as.Date(character(0))
    },
    ratio_label = attribution_ratio_label(
      arm_starts$factual,
      arm_starts$counterfactual
    ),
    smoothing_events = if (stitched && read_events) {
      load_smoothing_events(ctx, ctx$SMOOTHING_PARTS_DIR, attrib_start)
    } else {
      NULL
    }
  )
}

#' The smoothing partitions' _SUCCESS glob, a freshness input of every artefact
#' read from the smoothing (and of the forecast artefacts when stitched).
smoothing_parts_glob <- function(ctx) {
  list(list(
    dir = ctx$SMOOTHING_PARTS_DIR,
    pattern = "_SUCCESS$",
    recursive = TRUE
  ))
}

#' The per-arm frames of a load_group_arm_runs() result bound once per 09 job
#' (kept in ctx$RUN_LOAD_MEMO with the load they come from, so attribution
#' steps do not bind_rows() the same ~10^8 rows again):
#' list(adm1_month, country_month).
.bound_run_frames <- function(ctx, runs) {
  memo <- ctx$RUN_LOAD_MEMO
  same <- is.environment(memo) && identical(memo$value, runs) # same object: a pointer check
  if (same && !is.null(memo$bound)) {
    return(memo$bound)
  }
  bind1 <- function(l) if (length(l) == 1L) l[[1L]] else dplyr::bind_rows(l) # bound at load: no copy
  b <- list(
    adm1_month = bind1(runs$adm1_month),
    country_month = bind1(runs$country_month)
  )
  if (same) {
    memo$bound <- b
  }
  b
}

#' dplyr::left_join(df, lut, by = key) for a lookup table with unique keys: the
#' columns of `lut` matched onto `df` by `key`, in df's row order (missing keys
#' match missing keys, as dplyr). A lut with repeated keys, or a column df
#' already has, takes left_join().
.lookup_join <- function(df, lut, key) {
  if (
    anyDuplicated(lut[[key]]) || any(setdiff(names(lut), key) %in% names(df))
  ) {
    return(dplyr::left_join(df, lut, by = key))
  }
  i <- match(df[[key]], lut[[key]])
  for (v in setdiff(names(lut), key)) {
    df[[v]] <- lut[[v]][i]
  }
  df
}

#' The per-run frames of a load bound once, each as a one-element list (same
#' shape: bind_rows() of it is the bound frame): the per-run copies are freed
#' as soon as each table is bound, instead of living beside the bound copy for
#' the whole 09 job (~2x the run frames at 10^4 replicates).
.bind_run_lists <- function(country_month, adm1_month) {
  am <- dplyr::bind_rows(adm1_month)
  rm(adm1_month)
  cm <- dplyr::bind_rows(country_month)
  rm(country_month)
  list(country_month = list(all = cm), adm1_month = list(all = am))
}

#' load_group_arm_runs
#'
#' Load every included group's forecast run and its counterfactual run per
#' analog set over [attrib_start, group horizon], tagged with run / analog_set /
#' enso_group / group_horizon_end (shared by derive_attribution_aggregates()
#' and derive_report_forecast()).
#' `verbose` prints the per-run progress lines the attribution deriver logged.
#' List element names are bookkeeping only (bind_rows() drops them).
#'
#' `win` (results_arm_window(), optional): the factual arm is then read through
#' load_arm_run(), stitched from the smoothing when attrib_start precedes its
#' launch; NULL keeps the plain read.
#' @return list(country_month, adm1_month: lists of tagged frames,
#'   group_horizon: named character, cfgs: the run configs read,
#'   gids: one tibble (gid, country, enso_group) of each group's forecast-arm
#'   simulated gids; a counterfactual arm simulating a different gid set warns;
#'   stitch: per forecast run, its load's `stitch` element when stitched)
load_group_arm_runs <- function(
  ctx,
  groups,
  sets,
  attrib_start,
  verbose = FALSE,
  win = NULL
) {
  # One load per 09 job: the last result is kept in ctx$RUN_LOAD_MEMO under a
  # hash of everything the loaders read (manifest, results block, sim root,
  # smoothing parts, arm alignment, groups, sets, start, window), so another step
  # asking for the same runs reuses it and a different request (another event
  # view, other sets) never does. The previous load is dropped before a new one
  # (one load in memory at a time).
  memo <- ctx$RUN_LOAD_MEMO
  if (!is.environment(memo)) {
    return(.load_group_arm_runs_impl(
      ctx,
      groups,
      sets,
      attrib_start,
      verbose,
      win
    ))
  }
  key <- rlang::hash(list(
    ctx$MAN,
    ctx$RES,
    ctx$SIM_ROOT,
    ctx$SMOOTHING_PARTS_DIR,
    ctx$ARM_ALIGNMENT,
    groups,
    sets,
    as.character(attrib_start),
    win,
    .africa_load_ends(ctx, ctx$MAN, groups)
  ))
  if (identical(memo$key, key)) {
    message("[09]   reusing this job's run load (same runs, window and sets)")
    return(memo$value)
  }
  memo$key <- NULL
  memo$value <- NULL
  memo$bound <- NULL
  out <- .load_group_arm_runs_impl(
    ctx,
    groups,
    sets,
    attrib_start,
    verbose,
    win
  )
  memo$key <- key
  memo$value <- out
  out
}

.load_group_arm_runs_impl <- function(
  ctx,
  groups,
  sets,
  attrib_start,
  verbose = FALSE,
  win = NULL
) {
  man <- ctx$MAN
  if (length(man$units)) {
    return(.load_unit_arm_runs(ctx, groups, sets, attrib_start, verbose, win))
  }
  group_horizon <- setNames(
    vapply(
      groups,
      function(g) as.character(man$groups[[g]]$horizon_end),
      character(1)
    ),
    groups
  )
  load_end <- .africa_load_ends(ctx, man, groups) # past the firm horizon: the Africa window's continuation
  country_month <- list()
  adm1_month <- list()
  cfgs <- character(0)
  gids <- list()
  stitch <- list()
  for (g in groups) {
    horizon_end_g <- as.Date(group_horizon[[g]])
    read_end_g <- if (is.null(load_end)) {
      horizon_end_g
    } else {
      max(horizon_end_g, load_end[[g]])
    }
    tag_f <- group_run_tag(man, g, "forecast")
    cfgs <- c(cfgs, group_config_path(man, g, "forecast"))
    if (verbose) {
      message(sprintf(
        "[09]   group %-12s forecast=%s (-> %s)",
        g,
        tag_f,
        read_end_g
      ))
    }
    run_f <- load_run(
      ctx,
      tag_f,
      attrib_start,
      read_end_g,
      arm = "forecast",
      win = win
    )
    stitch[[g]] <- run_f$stitch
    gids[[g]] <- run_f$gids %>% dplyr::transmute(gid, country, enso_group = g)
    country_month[[paste0(g, "|forecast")]] <- run_f$country_month %>%
      dplyr::mutate(
        run = "forecast",
        analog_set = NA_character_,
        enso_group = g,
        group_horizon_end = horizon_end_g
      )
    adm1_month[[paste0(g, "|forecast")]] <- run_f$adm1_month %>%
      dplyr::mutate(
        run = "forecast",
        analog_set = NA_character_,
        enso_group = g,
        group_horizon_end = horizon_end_g
      )
    for (aset in sets) {
      tag_c <- group_run_tag(man, g, "counterfactual", analog_set = aset)
      cfgs <- c(cfgs, group_config_path(man, g, "counterfactual", aset))
      if (verbose) {
        message(sprintf(
          "[09]   group %-12s counterfactual[%s]=%s",
          g,
          aset,
          tag_c
        ))
      }
      run_c <- load_run(
        ctx,
        tag_c,
        attrib_start,
        read_end_g,
        arm = "counterfactual",
        win = win
      )
      if (!setequal(run_c$gids$gid, gids[[g]]$gid)) {
        warning(sprintf(
          "[09] group %s: counterfactual[%s] simulates a different gid set than the forecast arm",
          g,
          aset
        ))
      }
      key <- paste0(g, "|cf|", aset)
      country_month[[key]] <- run_c$country_month %>%
        dplyr::mutate(
          run = "counterfactual",
          analog_set = aset,
          enso_group = g,
          group_horizon_end = horizon_end_g
        )
      adm1_month[[key]] <- run_c$adm1_month %>%
        dplyr::mutate(
          run = "counterfactual",
          analog_set = aset,
          enso_group = g,
          group_horizon_end = horizon_end_g
        )
    }
  }
  c(
    .bind_run_lists(country_month, adm1_month),
    list(
      group_horizon = group_horizon,
      cfgs = cfgs,
      gids = dplyr::bind_rows(gids),
      stitch = stitch
    )
  )
}

#' .load_unit_arm_runs
#'
#' v2 manifests (partition_by: registry): runs are per partition unit and a
#' unit may mix fit groups. Each unit run is loaded once over [attrib_start,
#' unit horizon], every country censored at its own fit group's firm horizon
#' (results.scope.beyond_firm_horizon: "drop" (default) or "flag", which keeps
#' the later rows flagged), then every row is relabelled with its country's
#' FIT group -- Figures 3 and 4 stay keyed by the fit grouping. Same return
#' shape as load_group_arm_runs(), plus a `unit` column. A group the Africa
#' window reads past its firm horizon (results.scope.africa) is censored at
#' its end in that window instead; `group_horizon_end` stays the firm horizon.
.load_unit_arm_runs <- function(
  ctx,
  groups,
  sets,
  attrib_start,
  verbose = FALSE,
  win = NULL
) {
  man <- ctx$MAN
  cg <- manifest_country_groups(man)
  firm <- manifest_country_horizons(man, "firm")
  load_end <- .africa_load_ends(ctx, man, groups)
  if (!is.null(load_end)) {
    ext <- names(firm)[cg[names(firm)] %in% names(load_end)]
    firm[ext] <- format(pmax(
      as.Date(firm[ext]),
      as.Date(unname(load_end[cg[ext]]))
    ))
  }
  keep <- identical(ctx$RES$scope$beyond_firm_horizon %||% "drop", "flag")
  group_horizon <- setNames(
    vapply(
      groups,
      function(g) as.character(man$groups[[g]]$horizon_end),
      character(1)
    ),
    groups
  )
  units <- manifest_units(man)
  units <- units[vapply(
    units,
    function(u) any(groups %in% u$groups),
    logical(1)
  )]
  country_month <- list()
  adm1_month <- list()
  cfgs <- character(0)
  gids <- list()
  stitch <- list()
  tag_of <- function(cfg) yaml::read_yaml(cfg)$simulator$run_tag
  # Group and horizon looked up per distinct ISO3 / group, not per row (~10^8
  # rows at 10^4 replicates: as.Date() on a character per row was ~15 min).
  gh_date <- stats::setNames(
    as.Date(unname(group_horizon)),
    names(group_horizon)
  )
  relabel <- function(df, iso, arm_run, arm_set, un) {
    u_iso <- unique(iso)
    g <- unname(cg[u_iso])[match(iso, u_iso)]
    df %>%
      dplyr::mutate(
        run = .env$arm_run,
        analog_set = .env$arm_set,
        enso_group = g,
        group_horizon_end = unname(gh_date[g]),
        unit = .env$un
      ) %>%
      dplyr::filter(enso_group %in% groups)
  }
  for (un in names(units)) {
    u <- units[[un]]
    h_u <- as.Date(u$horizon_end)
    arms <- c(
      list(list(
        run = "forecast",
        aset = NA_character_,
        cfg = u$factual_config
      )),
      lapply(sets, function(a) {
        primary <- as.character(
          u$primary_analog_set %||% man$primary_analog_set %||% "all"
        )
        cfg <- u$counterfactual_sets[[a]]$config %||%
          (if (identical(a, primary)) u$counterfactual_config else NULL)
        if (is.null(cfg)) {
          stop(
            sprintf(
              "[09] unit '%s' has no counterfactual config for analog_set '%s'",
              un,
              a
            ),
            call. = FALSE
          )
        }
        list(run = "counterfactual", aset = a, cfg = cfg)
      })
    )
    for (a in arms) {
      tag <- tag_of(a$cfg)
      cfgs <- c(cfgs, a$cfg)
      if (verbose) {
        message(sprintf(
          "[09]   unit %-12s %s%s=%s (-> %s; groups %s)",
          un,
          a$run,
          if (is.na(a$aset)) "" else paste0("[", a$aset, "]"),
          tag,
          h_u,
          paste(u$groups, collapse = ",")
        ))
      }
      r <- load_run(
        ctx,
        tag,
        attrib_start,
        h_u,
        arm = a$run,
        win = win,
        horizon_by_country = firm,
        keep_beyond = keep
      )
      if (!is.null(r$stitch)) {
        stitch[[un]] <- r$stitch
      }
      key <- paste(un, a$run, a$aset, sep = "|")
      country_month[[key]] <- relabel(
        r$country_month,
        r$country_month$country,
        a$run,
        a$aset,
        un
      )
      adm1_month[[key]] <- relabel(
        r$adm1_month,
        substr(r$adm1_month$adm1_unit, 1L, 3L),
        a$run,
        a$aset,
        un
      )
      if (identical(a$run, "forecast")) {
        gids[[un]] <- r$gids %>%
          dplyr::transmute(
            gid,
            country,
            enso_group = unname(cg[country]),
            unit = un
          ) %>%
          dplyr::filter(enso_group %in% groups)
      }
    }
  }
  c(
    .bind_run_lists(country_month, adm1_month),
    list(
      group_horizon = group_horizon,
      cfgs = cfgs,
      gids = dplyr::bind_rows(gids),
      stitch = stitch
    )
  )
}

#' Replicate ids present in EVERY element of `dfs` (list of data frames each
#' with a replicate_id column). The Africa-total sum/ratio is only valid over
#' a replicate id that exists in every included group's run -- silently
#' summing over whatever happens to be present would pool mismatched draws.
common_replicate_ids <- function(dfs) {
  ids <- lapply(dfs, function(d) unique(d$replicate_id))
  common <- Reduce(intersect, ids)
  ns <- vapply(ids, length, integer(1))
  if (length(unique(ns)) > 1L || length(common) < max(ns)) {
    message(sprintf(
      "[09]   replicate-id pairing: group counts = %s -> %d common ids used",
      paste(ns, collapse = ","),
      length(common)
    ))
  }
  common
}

# Explosive-replicate cap: a replicate whose
# group's cumulative cases over a window reach results.thresholds.explosive_case_cap
# (ctx$EXPLOSIVE_CAP) in either
# arm is dropped from that group's pair under explosive mode "drop"
# (results.derivation.explosive_modes, ctx$EXPLOSIVE_MODES; "none" keeps every
# replicate). Every helper takes the cap as a required argument. The cap is
# c(case = <group window cap>, gid_month = <per-gid
# monthly cap>) (.ctx_explosive_cap(), .cell_cap()): a replicate is explosive
# when its group window total reaches `case` OR any gid of the group reaches
# `gid_month` in any window month (results.thresholds.explosive_gid_month_cap,
# ctx$EXPLOSIVE_GID_CAP; Inf = none, the former rule). A bare number is the
# group cap alone.

#' The group window cap / the per-gid monthly cap of an explosive cap spec.
.cap_case <- function(cap) {
  unname(if ("case" %in% names(cap)) cap[["case"]] else cap[[1L]])
}
.cap_gid <- function(cap) {
  unname(if ("gid_month" %in% names(cap)) cap[["gid_month"]] else Inf)
}

#' Whether replicates are explosive under `cap`: the group window total
#' (`cum_cases`) at or above the group cap, or the largest gid-month
#' (`max_gid_cases`, window_cum_per_rep()) at or above the per-gid cap. A
#' finite per-gid cap without max_gid_cases is an error, never a silent pass.
.is_explosive <- function(cum_cases, max_gid_cases, cap) {
  e <- !(cum_cases < .cap_case(cap))
  g <- .cap_gid(cap)
  if (is.finite(g)) {
    if (is.null(max_gid_cases)) {
      stop(
        "[09] the per-gid monthly explosive cap needs max_gid_cases (load_sim_run_monthly()) in the window totals.",
        call. = FALSE
      )
    }
    e <- e | !(max_gid_cases < g)
  }
  e
}

#' dplyr::anti_join(x, y, by) for the ~10^8-row ADM1 tables of a 10^4-replicate
#' run: the rows of x whose `by` keys are not in y, in x's order, missing keys
#' matching missing keys (dplyr's na_matches = "na"). data.table on-join; an
#' empty x or y takes the dplyr path. Identical output:
#' tests/test_results_derive_fast.R.
.anti_join_rows <- function(x, y, by) {
  if (!nrow(x)) {
    return(dplyr::anti_join(x, y, by = by))
  }
  # Nothing to drop: x itself (what anti_join() returns, without the copy).
  if (!nrow(y)) {
    return(x)
  }
  hit <- data.table::as.data.table(x[, by, drop = FALSE])[
    data.table::as.data.table(unique(y[, by, drop = FALSE])),
    on = by,
    which = TRUE,
    nomatch = 0L
  ]
  if (!length(hit)) {
    return(x)
  }
  x[-sort(unique(hit)), , drop = FALSE]
}

#' The distinct (enso_group, run, analog_set, replicate_id) keys of `kept`: the
#' paired-replicate lookups of attrib_exposure_tables() / attrib_metrics_tables()
#' read this instead of filtering the whole ADM1 table once per (set, group).
.replicate_keys <- function(kept) {
  k <- c("enso_group", "run", "analog_set", "replicate_id")
  as.data.frame(unique(data.table::as.data.table(kept[, k, drop = FALSE])))
}

#' Per-replicate cumulative cases / onsets over each group's attribution window
#' [attrib_start, group_horizon_end] of a country_month frame, grouped by `keys`
#' (e.g. "enso_group", or c("country", "enso_group")) then run, analog_set and
#' replicate_id.
window_cum_per_rep <- function(country_month, attrib_start, keys) {
  w <- dplyr::filter(
    country_month,
    month >= attrib_start,
    month <= group_horizon_end
  )
  by <- c(keys, "run", "analog_set", "replicate_id")
  out <- .sum_by_keys(
    w,
    by,
    c(cum_cases = "n_cases", cum_outbreaks = "n_outbreaks")
  )
  # The largest gid-month of the window (per-gid explosive cap); same keys,
  # same sort as the sums.
  if ("max_gid_cases" %in% names(w)) {
    out$max_gid_cases <- .sum_by_keys(
      w,
      by,
      c(max_gid_cases = "max_gid_cases"),
      op = "max"
    )$max_gid_cases
  }
  out
}

#' Rows of a window_cum_per_rep() table below the explosive-replicate cap.
drop_explosive_replicates <- function(cum, cap) {
  cum[
    !.is_explosive(cum$cum_cases, cum[["max_gid_cases"]], cap),
    ,
    drop = FALSE
  ]
}

#' The group-level replicates drop_explosive_replicates() removes: the exact
#' complement, as (enso_group, run, analog_set, replicate_id) keys plus their
#' cum_cases.
explosive_replicate_keys <- function(country_month, attrib_start, cap) {
  cum <- window_cum_per_rep(country_month, attrib_start, "enso_group")
  cum[
    .is_explosive(cum$cum_cases, cum[["max_gid_cases"]], cap),
    ,
    drop = FALSE
  ] %>%
    dplyr::select(enso_group, run, analog_set, replicate_id, cum_cases)
}

#' Per-country window totals (window_cum_per_rep(..., c("country", "enso_group")))
#' without the GROUP-level explosive replicates (explosive_replicate_keys()):
#' Figure 4a's country ratios drop the same replicates as their group's.
country_cum_kept <- function(country_month, start, cap) {
  cc <- window_cum_per_rep(country_month, start, c("country", "enso_group"))
  if (!is.finite(.cap_case(cap)) && !is.finite(.cap_gid(cap))) {
    return(cc)
  }
  .anti_join_rows(
    cc,
    explosive_replicate_keys(country_month, start, cap),
    c("enso_group", "run", "analog_set", "replicate_id")
  )
}

#' explosive_pair_summary (table S29)
#'
#' What the explosive guard ("drop" mode) removes from the pairs of one window,
#' per analog set and scope (each group, and Africa: the replicate ids common
#' to every group on both arms, dropped when explosive in any group or arm).
#' Dropping a pair when either arm reaches the cap conditions on the outcome:
#' where ENSO raises incidence the factual arm reaches it more often, so the
#' guard pulls the AF towards 0.
#'
#' @param cum per-replicate group window totals: tibble(enso_group, run
#'   ("forecast" | "counterfactual"), analog_set, replicate_id, cum_cases)
#'   (window_cum_per_rep(..., "enso_group"))
#' @param cap the explosive cap (a bare group cap, or c(case, gid_month)); with
#'   a per-gid cap, `cum` needs max_gid_cases; @param sets the counterfactual sets
#' @return one row per set x scope: level, scope (group id or "Africa"),
#'   n_rep_f, n_rep_cf (replicates per arm), n_expl_f, n_expl_cf (replicates at
#'   the cap per arm), n_pairs, n_drop_f_only, n_drop_cf_only, n_drop_both,
#'   n_pairs_kept, share_dropped, sign_p (exact two-sided binomial test of
#'   F-only vs CF-only drops: is the guard one-sided?), an_dropped (sum of
#'   F - CF cases over the dropped pairs), af_none, af_drop (AF of the mean
#'   cases with every pair / the kept pairs), af_none_minus_drop
explosive_pair_summary <- function(cum, cap, sets) {
  # Both caps first: inside tibble() a column named `cap` would mask the argument.
  cap_case <- .cap_case(cap)
  cap_gid <- .cap_gid(cap)
  cum <- dplyr::mutate(cum, run = as.character(run))
  groups <- sort(unique(cum$enso_group))
  af_of <- function(f, cf) {
    mf <- mean(f)
    if (length(f) && isTRUE(mf > 0)) 1 - mean(cf) / mf else NA_real_
  }
  one <- function(level, scope, aset, f, cf, fx, cfx) {
    # f, cf: named numeric (replicate id -> window cases) of one scope;
    # fx, cfx: named logical (replicate id -> explosive)
    ids <- intersect(names(f), names(cf))
    fe <- fx[ids]
    ce <- cfx[ids]
    drop <- fe | ce
    n_fo <- sum(fe & !ce)
    n_co <- sum(ce & !fe)
    tibble::tibble(
      level = level,
      scope = scope,
      analog_set = aset,
      cap = cap_case,
      gid_month_cap = cap_gid,
      n_rep_f = length(f),
      n_rep_cf = length(cf),
      n_expl_f = sum(fx),
      n_expl_cf = sum(cfx),
      n_pairs = length(ids),
      n_drop_f_only = n_fo,
      n_drop_cf_only = n_co,
      n_drop_both = sum(fe & ce),
      n_pairs_kept = sum(!drop),
      share_dropped = if (length(ids)) mean(drop) else NA_real_,
      sign_p = if (n_fo + n_co > 0) {
        stats::binom.test(n_fo, n_fo + n_co, 0.5)$p.value
      } else {
        NA_real_
      },
      an_dropped = sum(f[ids][drop] - cf[ids][drop]),
      af_none = af_of(f[ids], cf[ids]),
      af_drop = af_of(f[ids][!drop], cf[ids][!drop]),
      af_none_minus_drop = af_none - af_drop
    )
  }
  cum$.expl <- .is_explosive(cum$cum_cases, cum[["max_gid_cases"]], cap)
  vec <- function(d) stats::setNames(d$cum_cases, as.character(d$replicate_id))
  flg <- function(d) stats::setNames(d$.expl, as.character(d$replicate_id))
  dplyr::bind_rows(lapply(sets, function(aset) {
    fr <- cum[cum$run == "forecast", , drop = FALSE]
    cr <- cum[
      cum$run == "counterfactual" & cum$analog_set %in% aset,
      ,
      drop = FALSE
    ]
    rows <- lapply(groups, function(g) {
      fg <- fr[fr$enso_group == g, ]
      cg <- cr[cr$enso_group == g, ]
      one("region", g, aset, vec(fg), vec(cg), flg(fg), flg(cg))
    })
    # Africa: ids present in every group on both arms; a pair is dropped when
    # explosive in any group of either arm.
    ids <- Reduce(
      intersect,
      lapply(groups, function(g) {
        intersect(
          as.character(fr$replicate_id[fr$enso_group == g]),
          as.character(cr$replicate_id[cr$enso_group == g])
        )
      })
    )
    tot <- function(d) {
      d <- d[as.character(d$replicate_id) %in% ids, ]
      list(
        sum = tapply(d$cum_cases, as.character(d$replicate_id), sum),
        any = tapply(d$.expl, as.character(d$replicate_id), any)
      )
    }
    tf <- tot(fr)
    tc <- tot(cr)
    fe <- tf$any[ids]
    ce <- tc$any[ids]
    drop <- fe | ce
    f <- tf$sum[ids]
    cf <- tc$sum[ids]
    n_fo <- sum(fe & !ce)
    n_co <- sum(ce & !fe)
    afr <- tibble::tibble(
      level = "continent",
      scope = "Africa",
      analog_set = aset,
      cap = cap_case,
      gid_month_cap = cap_gid,
      n_rep_f = length(ids),
      n_rep_cf = length(ids),
      n_expl_f = sum(fe),
      n_expl_cf = sum(ce),
      n_pairs = length(ids),
      n_drop_f_only = n_fo,
      n_drop_cf_only = n_co,
      n_drop_both = sum(fe & ce),
      n_pairs_kept = sum(!drop),
      share_dropped = if (length(ids)) mean(drop) else NA_real_,
      sign_p = if (n_fo + n_co > 0) {
        stats::binom.test(n_fo, n_fo + n_co, 0.5)$p.value
      } else {
        NA_real_
      },
      an_dropped = sum(f[drop] - cf[drop]),
      af_none = af_of(f, cf),
      af_drop = af_of(f[!drop], cf[!drop]),
      af_none_minus_drop = af_none - af_drop
    )
    dplyr::bind_rows(afr, rows)
  }))
}

#' explosive_pair_summary() of every window of a cell table on one load's
#' country_month frame (each window cut by .cell_frame()), rows tagged
#' `window` and the factual label (which factual the pairs use).
explosive_summary_cells <- function(country_month, cells, sets, factual) {
  w <- cells[!duplicated(cells$window), , drop = FALSE]
  cd <- cells[cells$explosive_mode == "drop", , drop = FALSE]
  if (!nrow(cd)) {
    return(tibble::tibble())
  }
  cap <- .cell_cap(cd[1L, ])
  cm <- dplyr::mutate(country_month, run = as.character(run))
  dplyr::bind_rows(lapply(seq_len(nrow(w)), function(i) {
    explosive_pair_summary(
      window_cum_per_rep(.cell_frame(cm, w[i, ]), w$start[i], "enso_group"),
      cap,
      sets
    ) %>%
      dplyr::mutate(window = w$window[i], factual = factual, .before = 1L)
  }))
}

# ---------------------------------------------------------------------------
# Attribution cells: ONE table of the window x
# explosive-mode cells every attribution step derives (report_attrib_exposure,
# report_attrib_metrics, af_posterior, attrib_exposure_posterior,
# report_attrib_information, the event attribution), instead of a windows /
# guards loop per step.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# One definition of Africa. Each ENSO group is forecast to its own firm horizon
# (the last ONI month the firm forecast covers + the group's minimum deployed
# lag). With
# results.scope.africa.window_end a date, Africa is every included group over
# [attribution start, window_end]; a group whose firm horizon ends earlier is
# read on into the simulated analog continuation of ENSO, up to
# continuation_max_lead months past its firm horizon (0 = firm only, "to_end" =
# to window_end) and, with max_counterfactual_winter ("YYYY/YY"), no further than
# its last month whose ONI reads stay before the following DJF (every track:
# t - lag <= Nov of the winter's second year). "legacy" keeps the two former
# definitions (and derives no common cell).
# ---------------------------------------------------------------------------

#' Add `k` months to first-of-month Dates (k may be Inf: returns NA).
.add_months <- function(d, k) {
  out <- .month_from_index(
    .month_index(as.Date(d)) + ifelse(is.finite(k), k, 0)
  )
  out[!is.finite(k)] <- NA
  out
}

#' Months from `a` to `b` (first-of-month Dates): b - a.
.months_between <- function(a, b) {
  .month_index(as.Date(b)) - .month_index(as.Date(a))
}

#' Parse results.scope.africa.continuation_max_lead: a non-negative integer, or
#' "to_end" (Inf).
.parse_continuation_max_lead <- function(x) {
  x <- unlist(x %||% 0L)
  if (length(x) == 1L && identical(as.character(x), "to_end")) {
    return(Inf)
  }
  v <- suppressWarnings(as.numeric(as.character(x)))
  if (length(v) != 1L || is.na(v) || v < 0 || v != round(v)) {
    stop(
      "[09] results.scope.africa.continuation_max_lead must be a non-negative integer or \"to_end\"; got ",
      paste(x, collapse = ","),
      call. = FALSE
    )
  }
  v
}

#' Parse results.scope.africa.max_counterfactual_winter ("" = no cap, or
#' "YYYY/YY": the last DJF winter a track may read): the last ONI month allowed
#' (November of the winter's second year), or NA.
.parse_max_cf_winter <- function(x) {
  x <- as.character(unlist(x %||% ""))
  if (!length(x) || !nzchar(x[1])) {
    return(as.Date(NA))
  }
  m <- regmatches(x, regexec("^([0-9]{4})/([0-9]{2})$", x))[[1]]
  if (length(m) != 3L || (as.integer(m[2]) + 1L) %% 100L != as.integer(m[3])) {
    stop(
      "[09] results.scope.africa.max_counterfactual_winter must be \"\" or \"YYYY/YY\" (e.g. \"2026/27\"); got '",
      x,
      "'.",
      call. = FALSE
    )
  }
  as.Date(sprintf("%d-11-01", as.integer(m[2]) + 1L))
}

#' The deployed lags of a manifest group (occurrence / duration / monthly).
.manifest_group_lags <- function(man, g) {
  x <- man$groups[[g]]
  stats::setNames(
    vapply(
      c("occ", "dur", "mon"),
      function(t) {
        as.integer(x[[paste0("best_lag_", t)]] %||% NA_integer_)
      },
      integer(1)
    ),
    c("occurrence", "duration", "monthly")
  )
}

#' results_africa_window
#'
#' Africa of the attribution steps (R2: the only definition): its groups, its
#' window end and each group's last month in it.
#'
#' @param man    forecast manifest (the displayed one; the attribution manifest
#'   is a copy with the same groups and horizons)
#' @param groups included groups (manifest_groups_inc())
#' @param africa results.scope.africa (window_end, continuation_max_lead,
#'   max_counterfactual_winter); NULL = legacy
#' @return list(mode = "legacy" | "common", groups, end, firm (named Date: firm
#'   horizons), group_ends (named Date), lead (named int: months each group is
#'   read past its firm horizon, >= 0), censored (groups ending before `end`),
#'   extended (some group is read past its firm horizon), firm_oni_end (the last
#'   month of the firm ENSO forecast, manifest_last_firm_oni()), lags (named
#'   list of per-track lags), min_lag, max_lead, max_cf_winter)
results_africa_window <- function(man, groups, africa = NULL) {
  africa <- africa %||% list()
  we <- as.character(unlist(africa[["window_end"]] %||% "legacy"))
  if (length(we) != 1L) {
    stop(
      "[09] results.scope.africa.window_end must be one value.",
      call. = FALSE
    )
  }
  firm <- stats::setNames(
    as.Date(vapply(
      groups,
      function(g) as.character(man$groups[[g]]$horizon_end %||% NA),
      character(1)
    )),
    groups
  )
  if (anyNA(firm)) {
    stop(
      "[09] Africa window: no firm horizon (groups.<g>.horizon_end) in the manifest for ",
      paste(groups[is.na(firm)], collapse = ", "),
      call. = FALSE
    )
  }
  lags <- lapply(stats::setNames(groups, groups), function(g) {
    .manifest_group_lags(man, g)
  })
  min_lag <- vapply(
    groups,
    function(g) {
      as.integer(man$groups[[g]]$min_lag %||% min(lags[[g]], na.rm = TRUE))
    },
    integer(1)
  )
  seam <- manifest_last_firm_oni(man)
  out <- list(
    groups = groups,
    firm = firm,
    lags = lags,
    min_lag = min_lag,
    firm_oni_end = if (is.null(seam)) {
      as.Date(NA)
    } else {
      as.Date(as.character(seam))
    }
  )
  if (identical(we, "legacy")) {
    end <- as.Date(man$africa_min_horizon %||% min(firm))
    ge <- stats::setNames(rep(end, length(groups)), groups)
    return(c(
      out,
      list(
        mode = "legacy",
        end = end,
        group_ends = ge,
        lead = stats::setNames(rep(0L, length(groups)), groups),
        censored = character(0),
        extended = FALSE,
        max_lead = 0,
        max_cf_winter = ""
      )
    ))
  }
  end <- suppressWarnings(as.Date(we))
  if (is.na(end) || format(end, "%d") != "01") {
    stop(
      "[09] results.scope.africa.window_end must be \"legacy\" or a first-of-month date \"YYYY-MM-01\"; got '",
      we,
      "'.",
      call. = FALSE
    )
  }
  max_lead <- .parse_continuation_max_lead(africa[["continuation_max_lead"]])
  wint <- as.character(unlist(africa[["max_counterfactual_winter"]] %||% ""))
  oni_cap <- .parse_max_cf_winter(wint)
  ge <- pmin(
    rep(end, length(groups)),
    .add_months(firm, max_lead),
    na.rm = TRUE
  )
  if (!is.na(oni_cap)) {
    ge <- pmin(ge, .add_months(rep(oni_cap, length(groups)), min_lag))
  }
  ge <- stats::setNames(as.Date(ge), groups)
  a0 <- suppressWarnings(as.Date(as.character(man$attribution_start %||% NA)))
  if (!is.na(a0) && any(ge < a0)) {
    stop(
      "[09] Africa window: group(s) ",
      paste(groups[ge < a0], collapse = ", "),
      " would end (",
      paste(format(ge[ge < a0]), collapse = ", "),
      ") before the attribution start ",
      format(a0),
      "; raise continuation_max_lead / max_counterfactual_winter or exclude them.",
      call. = FALSE
    )
  }
  lead <- pmax(0L, as.integer(.months_between(firm, ge)))
  c(
    out,
    list(
      mode = "common",
      end = end,
      group_ends = ge,
      lead = stats::setNames(lead, groups),
      censored = groups[ge < end],
      extended = any(ge > firm),
      max_lead = max_lead,
      max_cf_winter = wint
    )
  )
}

#' The Africa window of a context: ctx$AFRICA_WIN, else (a context built
#' without one, e.g. a test's) the legacy definition of its manifest.
.africa_win <- function(ctx, groups) {
  ctx$AFRICA_WIN %||% results_africa_window(ctx$MAN, groups, NULL)
}

#' Per-group last month the runs are loaded to: the firm horizon, or the
#' group's end in the Africa window when it reads past it. NULL = the firm
#' horizons (Africa window not extended, or not this manifest's).
.africa_load_ends <- function(ctx, man, groups) {
  aw <- ctx$AFRICA_WIN
  if (is.null(aw) || !isTRUE(aw$extended)) {
    return(NULL)
  }
  firm <- stats::setNames(
    as.Date(vapply(
      names(aw$firm),
      function(g) as.character(man$groups[[g]]$horizon_end %||% NA),
      character(1)
    )),
    names(aw$firm)
  )
  if (anyNA(firm) || any(firm != aw$firm)) {
    stop(
      "[09] the Africa window was resolved on another manifest (firm horizons differ).",
      call. = FALSE
    )
  }
  e <- pmax(aw$firm, aw$group_ends)
  e[intersect(names(e), groups)]
}

#' The Africa window as stored in the attribution artefacts (`africa_window`):
#' list(mode, start, end, n_months, groups, group_ends, firm, lead, censored,
#' max_lead, max_cf_winter, label).
africa_window_payload <- function(aw, attrib_start) {
  st <- as.Date(attrib_start)
  list(
    mode = aw$mode,
    start = st,
    end = aw$end,
    n_months = window_n_months(st, aw$end),
    groups = aw$groups,
    group_ends = aw$group_ends,
    firm = aw$firm,
    lead = aw$lead,
    censored = aw$censored,
    max_lead = aw$max_lead,
    max_cf_winter = aw$max_cf_winter,
    label = sprintf("%s–%s", format(st, "%b %Y"), format(aw$end, "%b %Y"))
  )
}

#' DJF winters ("YYYY/YY") spanned by ONI months [a, b] (first-of-month Dates).
djf_winters <- function(a, b) {
  if (is.na(a) || is.na(b) || b < a) {
    return(character(0))
  }
  m <- seq(as.Date(a), as.Date(b), by = "month")
  mo <- as.integer(format(m, "%m"))
  y <- as.integer(format(m, "%Y"))
  y1 <- ifelse(mo == 12L, y, ifelse(mo <= 2L, y - 1L, NA_integer_))
  y1 <- sort(unique(stats::na.omit(y1)))
  sprintf("%d/%02d", y1, (y1 + 1L) %% 100L)
}

#' The analog continuation's skill by lead (03_prepare_enso_continuation_library.R
#' --validate; results.inputs.enso_continuation_skill): tibble(lead,
#' crpss_vs_clim) of method "analog"; NULL when the path is "".
read_continuation_skill <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  s <- utils::read.csv(
    require_input(
      path,
      "ENSO continuation skill (03_prepare_enso_continuation_library.R --validate)"
    ),
    stringsAsFactors = FALSE
  )
  if (!all(c("method", "lead", "crpss_vs_clim") %in% names(s))) {
    stop(
      "[09] ",
      path,
      ": expected columns method, lead, crpss_vs_clim.",
      call. = FALSE
    )
  }
  s <- s[s$method == "analog", c("lead", "crpss_vs_clim")]
  tibble::tibble(
    lead = as.integer(s$lead),
    crpss_vs_clim = as.numeric(s$crpss_vs_clim)
  )
}

#' The DJF winter(s) the counterfactual ENSO analogs are constrained on, read
#' from the primary counterfactual config's Lenssen cache: a legacy cache's
#' set_definition$peak_season_year (one winter, unchanged); a library cache's
#' set_definition$winters (primary + the other read winters when held neutral,
#' several "YYYY/YY" values). NA when unreadable.
cf_constrained_winter <- function(man) {
  u <- manifest_units(man)
  if (!length(u)) {
    return(NA_character_)
  }
  cfg <- u[[1]]$counterfactual_config %||%
    u[[1]]$counterfactual_sets[[1]]$config
  if (is.null(cfg) || !file.exists(cfg)) {
    return(NA_character_)
  }
  cp <- yaml::read_yaml(cfg)$simulator$enso$counterfactual$cache_path
  if (is.null(cp) || !file.exists(cp)) {
    return(NA_character_)
  }
  sd <- tryCatch(readRDS(cp)$meta$set_definition, error = function(e) NULL)
  .djf <- function(y) {
    sprintf("%d/%02d", as.integer(y), (as.integer(y) + 1L) %% 100L)
  }
  # A library-mode cache records every classified winter:
  # the primary, plus the other read winters when they are held neutral
  # (other_winters "neutral"; "any" constrains the primary only).
  w <- sd$winters
  if (is.list(w) && length(w$primary)) {
    yrs <- as.integer(unlist(w$primary))
    if (identical(as.character(w$other_winters %||% "neutral"), "neutral")) {
      yrs <- c(yrs, as.integer(unlist(w$other)))
    }
    return(.djf(sort(unique(yrs))))
  }
  # Legacy cache: one primary winter (unchanged).
  y <- sd$peak_season_year
  if (is.null(y)) {
    return(NA_character_)
  }
  .djf(y)
}

# The constrained winter(s) of cf_constrained_winter() as one table cell: the
# single legacy winter as before, several library winters ", "-joined; NA kept.
.cf_winter_cell <- function(cf_winter) {
  w <- stats::na.omit(as.character(cf_winter))
  if (!length(w)) NA_character_ else paste(w, collapse = ", ")
}

#' The counterfactual-config Lenssen cache behind cf_constrained_winter() (a
#' freshness input of the attribution artefact), or character(0).
.cf_cache_path <- function(man) {
  u <- manifest_units(man)
  cfg <- if (length(u)) {
    u[[1]]$counterfactual_config %||% u[[1]]$counterfactual_sets[[1]]$config
  } else {
    NULL
  }
  if (is.null(cfg) || !file.exists(cfg)) {
    return(character(0))
  }
  cp <- yaml::read_yaml(cfg)$simulator$enso$counterfactual$cache_path
  if (is.null(cp) || !file.exists(cp)) character(0) else cp
}

#' africa_scope_table (table S28)
#'
#' One row per group of the Africa window: countries, firm horizon, end in the
#' window, months read past the firm horizon, and per track (occurrence /
#' duration / monthly) its lag and its longest continuation lead in the window
#' (the ONI month read at the group's end, month - lag, past the last firm ONI
#' month manifest_last_firm_oni(); 0 = no continuation) with the analog
#' continuation's CRPSS vs climatology at that lead (NA at lead 0: CFSv2 /
#' observed); the longest of them; and the DJF winters the window reads (the
#' ONI months [start - lag, end - lag] of every track) -- those other than the
#' counterfactual's constrained winter compare the factual's continuation with
#' an unconstrained counterfactual winter.
#' @param aw results_africa_window(); @param skill read_continuation_skill() or NULL
#' @param cf_winter cf_constrained_winter()
africa_scope_table <- function(
  aw,
  man,
  attrib_start,
  skill = NULL,
  cf_winter = NA_character_
) {
  st <- as.Date(attrib_start)
  crpss_at <- function(l) {
    if (is.null(skill) || is.na(l) || l <= 0L) {
      NA_real_
    } else {
      skill$crpss_vs_clim[match(l, skill$lead)]
    }
  }
  dplyr::bind_rows(lapply(aw$groups, function(g) {
    lg <- aw$lags[[g]]
    e <- aw$group_ends[[g]]
    lead_t <- vapply(
      lg,
      function(L) {
        if (is.na(L) || is.na(aw$firm_oni_end)) {
          NA_integer_
        } else {
          max(
            0L,
            as.integer(.months_between(aw$firm_oni_end, .add_months(e, -L)))
          )
        }
      },
      integer(1)
    )
    winters <- sort(unique(unlist(lapply(lg[!is.na(lg)], function(L) {
      djf_winters(.add_months(st, -L), .add_months(e, -L))
    }))))
    lmax <- suppressWarnings(max(lead_t, na.rm = TRUE))
    tibble::tibble(
      enso_group = g,
      n_countries = length(unlist(man$groups[[g]]$iso3)),
      firm_horizon = aw$firm[[g]],
      group_end = e,
      window_end = aw$end,
      months_past_firm = aw$lead[[g]],
      censored = e < aw$end,
      lag_occurrence = lg[["occurrence"]],
      lag_duration = lg[["duration"]],
      lag_monthly = lg[["monthly"]],
      lead_occurrence = lead_t[["occurrence"]],
      lead_duration = lead_t[["duration"]],
      lead_monthly = lead_t[["monthly"]],
      crpss_occurrence = crpss_at(lead_t[["occurrence"]]),
      crpss_duration = crpss_at(lead_t[["duration"]]),
      crpss_monthly = crpss_at(lead_t[["monthly"]]),
      lead_max = if (is.finite(lmax)) as.integer(lmax) else NA_integer_,
      crpss_at_lead_max = crpss_at(
        if (is.finite(lmax)) as.integer(lmax) else NA_integer_
      ),
      cf_winters_read = paste(winters, collapse = ", "),
      cf_winter_constrained = .cf_winter_cell(cf_winter),
      cf_winters_unconstrained = paste(
        setdiff(winters, stats::na.omit(cf_winter)),
        collapse = ", "
      ),
      firm_oni_end = aw$firm_oni_end
    )
  }))
}

#' attrib_cell_table
#'
#' The cells of the attribution steps (pure): one row per window x explosive
#' mode, windows outer, modes inner (the loop order of the former
#' windows x guards loops, so the bootstrap draws and row order are unchanged).
#' With a common Africa window (results_africa_window(), mode "common") the
#' legacy cells (attribution, post_fork: each group over its own firm window)
#' are followed by attribution_common and post_fork_common (every group to its
#' end in the Africa window) and attribution_firm (every group cut at
#' min(firm horizon, window end): no continuation; `censored_groups` lists the
#' groups ending before the window end).
#'
#' @param attrib_start    first month of the attribution window
#' @param post_fork_start the displayed factual's launch (display_factual_start());
#'   NULL = the attribution window only (event attribution: both arms launch
#'   at the event launch)
#' @param africa_win results_africa_window() (NULL or mode "legacy": the legacy
#'   cells only)
#' @param cap   explosive-replicate cap of the "drop" mode
#' @param modes explosive modes, in order: "drop" (a replicate whose group
#'   window cases reach `cap` in either arm is dropped from the pair) | "none"
#' @return tibble(window, start, end (NA: each group to its own firm horizon),
#'   group_ends (list: NULL, or named Date = each group's last month in the
#'   cell), censored_groups, explosive_mode, explosive_guard (mode == "drop"),
#'   cap (Inf: none))
attrib_cell_table <- function(
  attrib_start,
  post_fork_start = NULL,
  africa_win = NULL,
  cap,
  modes = c("drop", "none")
) {
  if (
    !length(modes) || !all(modes %in% c("drop", "none")) || anyDuplicated(modes)
  ) {
    stop(
      "[09] explosive modes must be distinct values of {drop, none}; got ",
      paste(modes, collapse = ","),
      call. = FALSE
    )
  }
  a <- as.Date(attrib_start)
  pf <- if (is.null(post_fork_start)) NULL else max(a, as.Date(post_fork_start))
  common <- !is.null(africa_win) && identical(africa_win$mode, "common")
  # Legacy cells: cut at the firm horizons when the load reaches past them.
  firm_ge <- if (common && isTRUE(africa_win$extended)) {
    africa_win$firm
  } else {
    NULL
  }
  w <- function(window, start, end = as.Date(NA), ge = NULL) {
    list(
      window = window,
      start = start,
      end = end,
      group_ends = ge,
      censored = if (is.null(ge) || is.na(end)) {
        ""
      } else {
        paste(names(ge)[ge < end], collapse = ",")
      }
    )
  }
  win <- list(w("attribution", a, ge = firm_ge))
  if (!is.null(pf)) {
    win <- c(win, list(w("post_fork", pf, ge = firm_ge)))
  }
  if (common) {
    e <- africa_win$end
    win <- c(win, list(w("attribution_common", a, e, africa_win$group_ends)))
    if (!is.null(pf)) {
      win <- c(win, list(w("post_fork_common", pf, e, africa_win$group_ends)))
    }
    win <- c(
      win,
      list(w(
        "attribution_firm",
        a,
        e,
        stats::setNames(
          as.Date(pmin(africa_win$firm, e)),
          names(africa_win$firm)
        )
      ))
    )
  }
  # Both caps first: inside tibble() the column `cap` would mask the argument.
  cap_case <- .cap_case(cap)
  cap_gid <- .cap_gid(cap)
  rows <- list()
  for (x in win) {
    for (md in modes) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        window = x$window,
        start = x$start,
        end = x$end,
        group_ends = list(x$group_ends),
        censored_groups = x$censored,
        explosive_mode = md,
        explosive_guard = identical(md, "drop"),
        cap = if (identical(md, "drop")) cap_case else Inf,
        gid_cap = if (identical(md, "drop")) cap_gid else Inf
      )
    }
  }
  dplyr::bind_rows(rows)
}

#' attrib_cells
#'
#' The attribution cells of this results context (attrib_cell_table() on
#' ctx$AFRICA_WIN).
attrib_cells <- function(ctx, attrib_start, post_fork_start = NULL) {
  attrib_cell_table(
    attrib_start,
    post_fork_start,
    africa_win = ctx$AFRICA_WIN,
    cap = .ctx_explosive_cap(ctx),
    modes = .ctx_explosive_modes(ctx)
  )
}

#' The explosive cap / modes of a context (results.thresholds.explosive_case_cap,
#' results.derivation.explosive_modes); a context without them is an error.
.ctx_explosive_cap <- function(ctx) {
  cap <- ctx$EXPLOSIVE_CAP
  if (is.null(cap)) {
    stop(
      "[09] the results context has no EXPLOSIVE_CAP (results.thresholds.explosive_case_cap).",
      call. = FALSE
    )
  }
  c(case = cap, gid_month = ctx$EXPLOSIVE_GID_CAP %||% Inf)
}

#' The explosive cap of one attribution cell (a row of attrib_cell_table()):
#' c(case, gid_month); a cell table without gid_cap has no per-gid cap.
.cell_cap <- function(cl) {
  c(
    case = cl$cap[[1L]],
    gid_month = if ("gid_cap" %in% names(cl)) cl$gid_cap[[1L]] else Inf
  )
}
.ctx_explosive_modes <- function(ctx) ctx$EXPLOSIVE_MODES %||% c("drop", "none")

#' The cells of one past ENSO event (event attribution, event_af_posterior):
#' the whole run from the launch, named by the event, x the explosive modes.
#' The Africa window is the forecast's: never applied to an event.
event_cells <- function(ctx, launch, event) {
  ec <- attrib_cell_table(
    launch,
    NULL,
    africa_win = NULL,
    cap = .ctx_explosive_cap(ctx),
    modes = .ctx_explosive_modes(ctx)
  )
  ec$window <- event
  ec
}

#' The cells of one per-winter sub-window of an event: the
#' event cells (one per explosive mode) named by the sub-window, starting at
#' its earliest bound and cut per key by .cell_frame() through `bounds`
#' (`sw_rows`: the sub-window's rows of sim_io.R::event_sub_windows_table()).
event_sub_window_cells <- function(ctx, launch, sw_rows) {
  ec <- event_cells(ctx, launch, sw_rows$sub_window[1L])
  ec$start <- min(sw_rows$start)
  ec$bounds <- rep(
    list(dplyr::select(sw_rows, level, key, start, end)),
    nrow(ec)
  )
  ec
}

#' Every cell of an event: the whole run (event_cells()) and, when the event
#' has sub-windows (`swt`, event_sub_windows_table()), each sub-window's.
#' Without sub-windows it is event_cells() itself.
event_cells_all <- function(ctx, launch, event, swt = NULL) {
  ec <- event_cells(ctx, launch, event)
  if (is.null(swt)) {
    return(ec)
  }
  dplyr::bind_rows(c(
    list(ec),
    lapply(split(swt, swt$k), function(r) {
      event_sub_window_cells(ctx, launch, r)
    })
  ))
}

#' The cells as stored in the attribution artefacts (`cells`): one row per
#' window x mode, group ends flattened ("enso_east=2027-10-01,...").
cells_payload <- function(cells) {
  dplyr::mutate(
    cells,
    group_ends = vapply(
      group_ends,
      function(ge) {
        if (is.null(ge)) {
          ""
        } else {
          paste(
            sprintf("%s=%s", names(ge), format(as.Date(ge))),
            collapse = ","
          )
        }
      },
      character(1)
    )
  )
}

#' The window starts of a cell table, one per window: named Date (window ->
#' start), the `windows` element of the attribution artefacts.
.cell_starts <- function(cells) {
  w <- cells[!duplicated(cells$window), , drop = FALSE]
  stats::setNames(w$start, w$window)
}

#' The rows of a run frame (adm1_month / country_month) in one cell: months
#' from the cell's start and, when the cell has group ends, up to its group's
#' end (group_horizon_end then set to it, so window_cum_per_rep() and the
#' explosive keys cut there too).
.cell_frame <- function(df, cell) {
  # A filter that keeps every row returns the frame itself, not a copy (the
  # usual case: a cell starting at the load's own start; ~10^8 rows at 10^4
  # replicates).
  keep0 <- df$month >= cell$start[[1]]
  out <- if (!anyNA(keep0) && all(keep0)) {
    df
  } else {
    dplyr::filter(df, month >= cell$start[[1]])
  }
  rm(keep0)
  ge <- cell$group_ends[[1]]
  if (!is.null(ge)) {
    e <- unname(as.Date(ge)[match(out$enso_group, names(ge))])
    keep <- !is.na(e) & out$month <= e
    if (!all(keep)) {
      out <- out[keep, , drop = FALSE]
    }
    out$group_horizon_end <- e[keep]
  }
  # Per-key bounds (an event's per-winter sub-window): only
  # cells built by event_sub_window_cells() have the column.
  if ("bounds" %in% names(cell) && !is.null(cell$bounds[[1]])) {
    out <- .cell_bounds_frame(out, cell$bounds[[1]])
  }
  out
}

#' Rows of a run frame inside per-key [start, end] bounds (event_sub_windows_table()
#' rows of one sub-window: level "country" keyed by the row's country, from
#' `country` or the ADM1 gid's ISO3 prefix, or level "group" keyed by
#' enso_group); group_horizon_end is set to the row's bound end.
.cell_bounds_frame <- function(df, bounds) {
  lev <- unique(bounds$level)
  if (length(lev) != 1L) {
    stop("[09] sub-window bounds mix levels.", call. = FALSE)
  }
  key <- if (identical(lev, "group")) {
    df$enso_group
  } else if ("country" %in% names(df)) {
    df$country
  } else {
    substr(df$adm1_unit, 1L, 3L)
  }
  i <- match(as.character(key), bounds$key)
  st <- bounds$start[i]
  en <- bounds$end[i]
  keep <- !is.na(i) & df$month >= st & df$month <= en
  out <- df[keep, , drop = FALSE]
  out$group_horizon_end <- en[keep]
  out
}

#' Tag a cell's table with its window and explosive guard (the leading columns
#' of every attribution table); `window = FALSE` for the event cells.
.cell_tag <- function(d, cell, window = TRUE) {
  if (window) {
    dplyr::mutate(
      d,
      window = cell$window[[1]],
      explosive_guard = cell$explosive_guard[[1]],
      .before = 1L
    )
  } else {
    dplyr::mutate(d, explosive_guard = cell$explosive_guard[[1]], .before = 1L)
  }
}

#' Monthly ribbon summary of replicate totals: quantiles of `n_cases` and
#' `n_outbreaks` per group of `df` (cases_q025 .. cases_q975, onsets_q025 ..
#' onsets_q975), the column set of every attribution ribbon.
.case_onset_ribbon <- function(df) {
  # One quantile() call per column and group (one sort) instead of seven; the
  # columns keep quantile()'s names ("2.5%", ...) as the per-call form did.
  P <- c(
    q025 = .025,
    q10 = .10,
    q25 = .25,
    q50 = .50,
    q75 = .75,
    q90 = .90,
    q975 = .975
  )
  out <- df %>%
    dplyr::summarise(
      .qc = list(stats::quantile(n_cases, P)),
      .qo = list(stats::quantile(n_outbreaks, P)),
      .groups = "drop"
    )
  qc <- matrix(as.numeric(unlist(out$.qc)), ncol = length(P), byrow = TRUE)
  qo <- matrix(as.numeric(unlist(out$.qo)), ncol = length(P), byrow = TRUE)
  lab <- names(stats::quantile(0, P))
  nm <- function(x, k) stats::setNames(x, rep(lab[k], length(x)))
  out$.qc <- NULL
  out$.qo <- NULL
  for (k in seq_along(P)) {
    out[[paste0("cases_", names(P)[k])]] <- nm(qc[, k], k)
  }
  for (k in seq_along(P)) {
    out[[paste0("onsets_", names(P)[k])]] <- nm(qo[, k], k)
  }
  out
}

derive_attribution_aggregates <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  win <- results_arm_window(ctx, groups, sets)
  attrib_start <- win$attrib_start
  africa_min_horizon <- as.Date(man$africa_min_horizon)
  aw <- .africa_win(ctx, groups) # Africa (results.scope.africa): groups, end, per-group ends
  message(sprintf(
    "[09] Figures 2-3: deriving aggregates | groups=%s | sets=%s | attrib_start=%s (%s; factual launch %s, counterfactual %s) | africa_min=%s | africa window %s (%s)",
    paste(groups, collapse = ","),
    paste(sets, collapse = ","),
    attrib_start,
    ctx$ARM_ALIGNMENT,
    win$arm_starts$factual,
    win$arm_starts$counterfactual,
    africa_min_horizon,
    aw$mode,
    format(aw$end)
  ))

  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  frame <- build_location_frame(
    survival_data = surv_dat,
    hmm_clusters_path = ctx$HMM_PATH,
    outbreak_metadata_path = NULL,
    size_data_path = ctx$SIZE_DATA_PATH
  )
  country_region <- frame %>% dplyr::distinct(country, AFRO_region)

  # Forecast + counterfactual frames of one manifest: the forecast manifest
  # (ctx; the displayed forecast, Figure 3/4 ribbons) or the attribution one
  # (results_attrib_ctx(); Figure 4's paired statistics).
  arm_frames <- function(cx, w) {
    runs <- load_group_arm_runs(
      cx,
      groups,
      sets,
      w$attrib_start,
      verbose = TRUE,
      win = w
    )
    bf <- .bound_run_frames(cx, runs)
    list(
      runs = runs,
      country_month = bf$country_month %>%
        dplyr::mutate(
          run = factor(run, levels = c("forecast", "counterfactual"))
        ) %>%
        .lookup_join(country_region, "country"),
      adm1_month = bf$adm1_month %>%
        dplyr::mutate(
          run = factor(run, levels = c("forecast", "counterfactual"))
        )
    )
  }
  fr <- arm_frames(ctx, win)
  runs <- fr$runs
  group_horizon <- runs$group_horizon
  cfgs <- runs$cfgs
  country_month_both <- fr$country_month
  adm1_month_both <- fr$adm1_month
  fr <- NULL

  # --- Figure 3 monthly ribbons -------------------------------------------
  # Per group, full own-horizon monthly ribbon.
  group_reps <- .sum_by_keys(
    country_month_both,
    c("run", "analog_set", "replicate_id", "month", "enso_group"),
    c(n_cases = "n_cases", n_outbreaks = "n_outbreaks")
  )
  # `continuation`: months past the group's firm horizon, read from the analog
  # continuation of ENSO (only when the Africa window reaches past it).
  cont_flag <- function(d) {
    dplyr::mutate(d, continuation = month > unname(aw$firm[enso_group]))
  }
  region_ribbon <- group_reps %>%
    dplyr::group_by(run, analog_set, month, enso_group) %>%
    .case_onset_ribbon() %>%
    cont_flag()

  # Per-country monthly ribbon (supplementary forecast-by-country figure; the
  # 07d_simulator_figures_counterfactual.R fig6/fig7 series), each country to its
  # own firm horizon (load_group_arm_runs() censors) and stitched like the groups.
  country_ribbon <- country_month_both %>%
    dplyr::group_by(run, analog_set, month, enso_group, country) %>%
    .case_onset_ribbon() %>%
    cont_flag()

  # Africa ribbon: paired replicate ids across every included group, each group
  # up to its end in the Africa window (legacy: africa_min_horizon for all);
  # n_groups_continuation = the groups read past their firm horizon that month.
  afr_end <- aw$group_ends
  africa_ribbon_of <- function(group_reps) {
    dplyr::bind_rows(lapply(unique(group_reps$run), function(rn) {
      asets_here <- if (rn == "forecast") NA_character_ else sets
      dplyr::bind_rows(lapply(asets_here, function(aset) {
        sub <- group_reps %>%
          dplyr::filter(
            run == rn,
            if (rn == "forecast") TRUE else analog_set == aset,
            month <= unname(afr_end[enso_group])
          )
        by_group <- split(sub, sub$enso_group)
        common <- common_replicate_ids(by_group)
        sub %>%
          dplyr::filter(replicate_id %in% common) %>%
          dplyr::group_by(replicate_id, month) %>%
          dplyr::summarise(
            n_cases = sum(n_cases),
            n_outbreaks = sum(n_outbreaks),
            .groups = "drop"
          ) %>%
          dplyr::mutate(run = rn, analog_set = aset)
      }))
    })) %>%
      dplyr::group_by(run, analog_set, month) %>%
      .case_onset_ribbon() %>%
      dplyr::mutate(
        run = factor(run, levels = c("forecast", "counterfactual")),
        n_groups_continuation = vapply(
          month,
          function(m) {
            sum(m > aw$firm & m <= afr_end)
          },
          integer(1)
        )
      )
  }
  africa_ribbon <- africa_ribbon_of(group_reps)

  # --- The attribution pair (results.inputs.attrib_manifest) ---------------
  # Everything below is a paired factual-vs-counterfactual statistic, so it reads
  # the attribution factual (launched at the counterfactual handoff) when there
  # is one; the ribbons above stay the displayed forecast. The display frames
  # are released first: one run load in memory at a time. attrib_ribbon /
  # attrib_africa_ribbon: the attribution factual's own ribbons (Figure 4).
  attrib_ribbon <- NULL
  attrib_africa_ribbon <- NULL
  win_a <- win
  if (identical(ctx$ATTRIB_FACTUAL, "attrib")) {
    country_month_both <- NULL
    adm1_month_both <- NULL
    group_reps <- NULL
    runs <- NULL
    invisible(gc(verbose = FALSE))
    cx_a <- results_attrib_ctx(ctx)
    win_a <- results_arm_window(cx_a, groups, sets)
    if (win_a$attrib_start != attrib_start) {
      stop(
        sprintf(
          "[09] attribution: the attribution window starts %s, the display window %s.",
          win_a$attrib_start,
          attrib_start
        ),
        call. = FALSE
      )
    }
    fa <- arm_frames(cx_a, win_a)
    cfgs <- unique(c(cfgs, fa$runs$cfgs))
    country_month_both <- fa$country_month
    adm1_month_both <- fa$adm1_month
    fa <- NULL
    reps_a <- .sum_by_keys(
      dplyr::filter(country_month_both, run == "forecast"),
      c("run", "analog_set", "replicate_id", "month", "enso_group"),
      c(n_cases = "n_cases", n_outbreaks = "n_outbreaks")
    )
    attrib_ribbon <- reps_a %>%
      dplyr::group_by(run, analog_set, month, enso_group) %>%
      .case_onset_ribbon() %>%
      cont_flag()
    attrib_africa_ribbon <- africa_ribbon_of(reps_a)
  }
  # The ADM1 frame of Figure 4b's probability map: each group's own firm window
  # (cut there when the Africa window loaded continuation months).
  adm1_firm <- if (isTRUE(aw$extended)) {
    dplyr::filter(adm1_month_both, month <= group_horizon_end)
  } else {
    adm1_month_both
  }
  adm1_group_horizon <- adm1_firm %>%
    dplyr::distinct(adm1_unit, enso_group, group_horizon_end)

  # --- Figure 4a: per-replicate PAIRED cumulative-case ratio --------------
  # Explosive replicates excluded (results.thresholds.explosive_case_cap).
  cap <- .ctx_explosive_cap(ctx)
  cum_per_rep_group <- drop_explosive_replicates(
    window_cum_per_rep(country_month_both, attrib_start, "enso_group"),
    cap
  )

  # Per-country counterpart of cum_per_rep_group, same attribution window --
  # feeds Figure 4a's bottom (country) panel, one ratio per country. The
  # replicates dropped are the GROUP-level explosive ones, as in every other
  # table.
  cum_per_rep_country <- country_cum_kept(country_month_both, attrib_start, cap)

  ratio_one <- function(df_f, df_c, value_col) {
    m <- dplyr::inner_join(
      df_f,
      df_c,
      by = "replicate_id",
      suffix = c("_f", "_c")
    )
    m$ratio <- m[[paste0(value_col, "_f")]] / m[[paste0(value_col, "_c")]]
    m$ratio[!is.finite(m$ratio)] <- NA_real_
    m
  }
  METRIC_COL_RATIO <- c(cases = "cum_cases", onsets = "cum_outbreaks")
  ratio_stats_of <- function(cum_per_rep_group, cum_per_rep_country) {
    dplyr::bind_rows(lapply(names(METRIC_COL_RATIO), function(metric) {
      val_col <- METRIC_COL_RATIO[[metric]]
      dplyr::bind_rows(lapply(sets, function(aset) {
        per_group <- dplyr::bind_rows(lapply(groups, function(g) {
          df_f <- dplyr::filter(
            cum_per_rep_group,
            enso_group == g,
            run == "forecast"
          )
          df_c <- dplyr::filter(
            cum_per_rep_group,
            enso_group == g,
            run == "counterfactual",
            analog_set == aset
          )
          m <- ratio_one(df_f, df_c, val_col)
          tibble::tibble(
            scope = ctx$GROUP_SHORT_RUN[[g]],
            analog_set = aset,
            ratio = m$ratio,
            fc_val = m[[paste0(val_col, "_f")]],
            level = "region",
            region = g
          )
        }))
        # Per-country ratio, nested within its own ENSO group -- `region` keeps
        # the group KEY (not the short label) so a plot can key straight into
        # GROUP_COLOURS/GROUP_LABELS, the convention Figure 1 uses for its
        # enso_group facets/palette.
        per_country <- dplyr::bind_rows(lapply(groups, function(g) {
          countries_g <- unique(cum_per_rep_country$country[
            cum_per_rep_country$enso_group == g
          ])
          dplyr::bind_rows(lapply(countries_g, function(ctry) {
            df_f <- dplyr::filter(
              cum_per_rep_country,
              country == ctry,
              enso_group == g,
              run == "forecast"
            )
            df_c <- dplyr::filter(
              cum_per_rep_country,
              country == ctry,
              enso_group == g,
              run == "counterfactual",
              analog_set == aset
            )
            m <- ratio_one(df_f, df_c, val_col)
            tibble::tibble(
              scope = ctry,
              analog_set = aset,
              ratio = m$ratio,
              fc_val = m[[paste0(val_col, "_f")]],
              level = "country",
              region = g
            )
          }))
        }))
        # Africa total: pair replicate ids common to every group, sum the metric
        # across groups on each side, THEN take the ratio.
        f_list <- lapply(groups, function(g) {
          dplyr::filter(cum_per_rep_group, enso_group == g, run == "forecast")
        })
        c_list <- lapply(groups, function(g) {
          dplyr::filter(
            cum_per_rep_group,
            enso_group == g,
            run == "counterfactual",
            analog_set == aset
          )
        })
        names(f_list) <- groups
        names(c_list) <- groups
        common_f <- common_replicate_ids(f_list)
        common_c <- common_replicate_ids(c_list)
        common_ids <- intersect(common_f, common_c)
        afr_f <- dplyr::bind_rows(f_list) %>%
          dplyr::filter(replicate_id %in% common_ids) %>%
          dplyr::group_by(replicate_id) %>%
          dplyr::summarise(
            dplyr::across(dplyr::all_of(val_col), sum),
            .groups = "drop"
          )
        afr_c <- dplyr::bind_rows(c_list) %>%
          dplyr::filter(replicate_id %in% common_ids) %>%
          dplyr::group_by(replicate_id) %>%
          dplyr::summarise(
            dplyr::across(dplyr::all_of(val_col), sum),
            .groups = "drop"
          )
        m <- ratio_one(afr_f, afr_c, val_col)
        dplyr::bind_rows(
          per_group,
          per_country,
          tibble::tibble(
            scope = "Africa",
            analog_set = aset,
            ratio = m$ratio,
            fc_val = m[[paste0(val_col, "_f")]],
            level = "continent",
            region = NA_character_
          )
        )
      })) %>%
        dplyr::mutate(metric = metric)
    })) %>%
      dplyr::group_by(metric, scope, analog_set, level, region) %>%
      dplyr::summarise(
        med = median(ratio, na.rm = TRUE),
        lo = quantile(ratio, .025, na.rm = TRUE, names = FALSE),
        q25 = quantile(ratio, .25, na.rm = TRUE, names = FALSE),
        q75 = quantile(ratio, .75, na.rm = TRUE, names = FALSE),
        hi = quantile(ratio, .975, na.rm = TRUE, names = FALSE),
        # Median forecast (factual-arm) case count across the
        # paired replicates, for Figure 4a's dot-size encoding
        # (bubble size = absolute scale, independent of the
        # ratio itself, which only encodes relative shift).
        med_forecast = median(fc_val, na.rm = TRUE),
        # Same PROB_THRESHOLD flag (results.thresholds.prob_threshold,
        # default 0.50) Figure 4b's ADM1 map already uses in BOTH
        # roles: the fractional shift ("> 50% more/fewer") AND
        # the significance probability ("> 50% chance of").
        # Paired per-replicate `ratio` (vs. fig4b's unpaired
        # outer() over raw sums) makes this the more appropriate
        # estimator here, but the threshold semantics match.
        prob_larger = mean(ratio > (1 + ctx$PROB_THRESHOLD), na.rm = TRUE),
        prob_smaller = mean(ratio < (1 - ctx$PROB_THRESHOLD), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        direction = factor(
          dplyr::case_when(
            prob_larger >= ctx$PROB_THRESHOLD ~ "increase",
            prob_smaller >= ctx$PROB_THRESHOLD ~ "decrease",
            TRUE ~ "no change"
          ),
          levels = c("increase", "no change", "decrease")
        )
      )
  }
  ratio_stats <- ratio_stats_of(cum_per_rep_group, cum_per_rep_country)

  # --- Figure 4b: ADM1 signed probability of shift, same window -----------
  METRIC_COL_ADM1 <- c(cases = "cases", onsets = "outbreaks")
  prob_df <- dplyr::bind_rows(lapply(names(METRIC_COL_ADM1), function(metric) {
    dplyr::bind_rows(lapply(sets, function(aset) {
      adm1_firm %>%
        dplyr::filter(run == "forecast" | analog_set == aset) %>%
        adm1_prob_stats(METRIC_COL_ADM1[[metric]], ctx$PROB_THRESHOLD) %>%
        finish_adm1_prob(
          low_signal_thresh = ctx$LOW_SIGNAL,
          sig_threshold = ctx$PROB_THRESHOLD,
          country_region = country_region,
          adm1_group_horizon = adm1_group_horizon
        ) %>%
        dplyr::mutate(analog_set = aset, metric = metric)
    }))
  }))

  # --- Figure 3c: paired replicate totals over the same window ------------
  # The cum_per_rep_group sums behind Figure 4a's ratio, summarised as absolute
  # totals per arm (both metrics from the same replicates, so the explosive-
  # replicate guard on cases applies to onsets too). Per analog set: Africa sums
  # the replicate ids present in every group on BOTH arms of that set, so its
  # forecast row also depends on the set.
  qv <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
  cum_totals_of <- function(cum_per_rep_group) {
    dplyr::bind_rows(lapply(sets, function(aset) {
      per_group <- cum_per_rep_group %>%
        dplyr::filter(run == "forecast" | analog_set == aset) %>%
        dplyr::transmute(
          enso_group,
          run = as.character(run),
          replicate_id,
          cases = cum_cases,
          onsets = cum_outbreaks
        )
      common <- Reduce(
        intersect,
        lapply(split(per_group, per_group$run), function(a) {
          Reduce(
            intersect,
            lapply(split(a, a$enso_group), function(d) unique(d$replicate_id))
          )
        })
      )
      africa <- per_group %>%
        dplyr::filter(replicate_id %in% common) %>%
        dplyr::group_by(run, replicate_id) %>%
        dplyr::summarise(
          cases = sum(cases),
          onsets = sum(onsets),
          .groups = "drop"
        ) %>%
        dplyr::mutate(enso_group = "Africa")
      dplyr::bind_rows(per_group, africa) %>%
        tidyr::pivot_longer(
          c(cases, onsets),
          names_to = "metric",
          values_to = "total"
        ) %>%
        dplyr::group_by(scope = enso_group, run, metric) %>%
        dplyr::summarise(
          n_rep = dplyr::n(),
          lo = qv(total, .025),
          q25 = qv(total, .25),
          med = qv(total, .5),
          q75 = qv(total, .75),
          hi = qv(total, .975),
          .groups = "drop"
        ) %>%
        dplyr::mutate(analog_set = aset)
    }))
  }
  cum_totals <- cum_totals_of(cum_per_rep_group)

  # --- Every attribution cell (attrib_cells(): window x explosive mode) ------
  # cum_totals / ratio_stats above are the attribution x drop cell, each group
  # over its own firm window (the legacy readers'); the *_cells tables repeat
  # them per cell, tagged window / explosive_guard / explosive_mode -- with a
  # common Africa window (results.scope.africa) its attribution_common /
  # post_fork_common / attribution_firm cells too, Africa summing the groups
  # over the same calendar window. Figure 4 reads the headline cell
  # (results.style.attrib_headline).
  cells <- attrib_cells(
    ctx,
    attrib_start,
    display_factual_start(ctx, groups, sets)
  )
  cell_rows <- lapply(seq_len(nrow(cells)), function(i) {
    cl <- cells[i, ]
    cm <- .cell_frame(country_month_both, cl)
    cg <- drop_explosive_replicates(
      window_cum_per_rep(cm, cl$start, "enso_group"),
      .cell_cap(cl)
    )
    cc <- country_cum_kept(cm, cl$start, .cell_cap(cl))
    tag <- function(d) {
      dplyr::mutate(
        d,
        window = cl$window,
        explosive_guard = cl$explosive_guard,
        explosive_mode = cl$explosive_mode,
        .before = 1L
      )
    }
    list(
      cum_totals = tag(cum_totals_of(cg)),
      ratio_stats = tag(ratio_stats_of(cg, cc))
    )
  })
  cum_totals_cells <- dplyr::bind_rows(lapply(cell_rows, `[[`, "cum_totals"))
  ratio_stats_cells <- dplyr::bind_rows(lapply(cell_rows, `[[`, "ratio_stats"))
  cell_rows <- NULL

  # Table S28 (africa_scope): per group, the window, continuation leads per
  # track, the analog continuation's CRPSS at them and the counterfactual
  # winters read.
  africa_scope <- africa_scope_table(
    aw,
    man,
    attrib_start,
    read_continuation_skill(ctx$CONT_SKILL_PATH),
    cf_constrained_winter(man)
  )

  list(
    region_ribbon = region_ribbon,
    africa_ribbon = africa_ribbon,
    country_ribbon = country_ribbon,
    ratio_stats = ratio_stats,
    prob_df = prob_df,
    cum_totals = cum_totals,
    cum_totals_cells = cum_totals_cells,
    ratio_stats_cells = ratio_stats_cells,
    cells = cells_payload(cells),
    africa_window = africa_window_payload(aw, attrib_start),
    africa_scope = africa_scope,
    explosive_cap = .cap_case(cap),
    explosive_gid_month_cap = .cap_gid(cap),
    explosive_modes = .ctx_explosive_modes(ctx),
    attrib_start = attrib_start,
    africa_min_horizon = africa_min_horizon,
    group_horizon = group_horizon,
    groups = groups,
    analog_sets = sets,
    # Per-arm launches (attribution v8): the Forecast series is the smoothing
    # output over stitch_months (empty when the arms launch together).
    factual_start = win$factual_start,
    arm_starts = win$arm_starts,
    stitch_months = win$stitch_months,
    ratio_label = win$ratio_label,
    # The factual of the paired statistics (ratio_stats, prob_df, cum_totals):
    # "attrib" = the attribution factual (attrib_arm_starts, its ribbons
    # attrib_ribbon / attrib_africa_ribbon), "forecast" = the displayed one.
    attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
    attrib_arm_starts = win_a$arm_starts,
    attrib_ratio_label = win_a$ratio_label,
    attrib_ribbon = attrib_ribbon,
    attrib_africa_ribbon = attrib_africa_ribbon,
    .inputs = results_inputs(files = c(cfgs, ctx$SURV_PATH, ctx$SIZE_DATA_PATH))
  )
}

#' Forecast-manifest + attribution-manifest run configs of the steps that read
#' both (attribution, report_forecast); the forecast manifest's alone without an
#' attribution factual.
.display_attrib_configs <- function(ctx, groups, sets) {
  cfgs <- results_run_configs(ctx$MAN, groups, sets)
  if (identical(ctx$ATTRIB_FACTUAL, "attrib")) {
    cfgs <- unique(c(cfgs, results_run_configs(ctx$MAN_ATTRIB, groups, sets)))
  }
  cfgs
}

ensure_attribution <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- .display_attrib_configs(ctx, groups, sets)
  inputs <- results_inputs(
    files = c(
      unique(c(
        ctx$MANIFEST_PATH,
        ctx$ATTRIB_MANIFEST_FILE %||% ctx$MANIFEST_PATH
      )),
      cfgs,
      # table S28's inputs: continuation skill, counterfactual analog cache
      if (nzchar(ctx$CONT_SKILL_PATH %||% "")) ctx$CONT_SKILL_PATH,
      .cf_cache_path(man)
    ),
    globs = if (
      results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched
    ) {
      smoothing_parts_glob(ctx)
    } else {
      list()
    }
  )

  ensure_artefact(
    ctx,
    "attribution",
    inputs,
    function() derive_attribution_aggregates(ctx, sets),
    sets = sets,
    strip_inputs = TRUE
  )
}

# ===========================================================================
# Section 3: historical smoothing derivation (multi-scope, multi-metric)
# ===========================================================================

#' load_smoothing_scope
#'
#' Shared set-up of derive_smoothing() and derive_report_smoothing():
#' require the raw inputs, build the location frame
#' (dropping ADM1 units with nested ADM2 children; `label` names the caller in
#' that filter's log line), load the ENSO grouping and the smoothing events,
#' and keep the included groups' onset (`events_all`) and mon-track case
#' (`events_mon_all`) events.
#' @return list(groups, surv_dat, frame, enso_grp, events_all, events_mon_all)
load_smoothing_scope <- function(ctx, label) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  for (p in c(
    ctx$PANEL_PATH,
    ctx$SURV_PATH,
    ctx$SIZE_DATA_PATH,
    ctx$MON_DATA_PATH
  )) {
    require_input(p, "raw output")
  }

  surv_dat <- readRDS(ctx$SURV_PATH)
  frame <- build_location_frame(
    readRDS(ctx$PANEL_PATH),
    hmm_clusters_path = if (file.exists(ctx$HMM_PATH)) ctx$HMM_PATH else NULL,
    size_data_path = ctx$SIZE_DATA_PATH
  )
  frame <- drop_overlapping_adm1_from_frame(frame, surv_dat, label)

  enso_grp <- load_enso_grouping(ctx$GROUPING_PATH)
  ev_raw <- load_smoothing_events(
    ctx,
    ctx$SMOOTHING_PARTS_DIR,
    ctx$SMOOTH_START
  )
  .attach_group <- function(ev) {
    ev %>%
      dplyr::left_join(enso_grp, by = c("country" = "iso3")) %>%
      dplyr::filter(group %in% groups)
  }
  events_all <- .attach_group(ev_raw$events) # onsets metric
  events_mon_all <- .attach_group(ev_raw$events_mon) # cases metric (mon track)
  list(
    groups = groups,
    surv_dat = surv_dat,
    frame = frame,
    enso_grp = enso_grp,
    events_all = events_all,
    events_mon_all = events_mon_all
  )
}

#' smoothing_scope_gids
#'
#' The gid sets the smoothing ribbons are summed over: "Africa" (every frame
#' gid in a country with smoothing events) plus one per included ENSO group,
#' named by GROUP_SHORT.
smoothing_scope_gids <- function(
  frame,
  events_all,
  enso_grp,
  groups,
  group_short = GROUP_SHORT
) {
  scope_gids <- list(
    Africa = unique(as.character(frame$gid[
      frame$country %in% events_all$country
    ]))
  )
  for (g in groups) {
    isos_g <- enso_grp$iso3[enso_grp$group == g]
    scope_gids[[group_short[[g]]]] <- unique(as.character(frame$gid[
      frame$country %in% isos_g
    ]))
  }
  scope_gids
}

#' Read the sparse event tables straight from the 08b smoothing partitions, the
#' same source analysis/08c_smoothing_compile.R itself binds together. Returns
#' both channels:
#'   $events     (part.qs -> $events)     one row per Yout==1 ONSET -- drives the
#'                                        onsets metric.
#'   $events_mon (part.qs -> $events_mon)  one row per (rep, gid, ACTIVE month)
#'                                        carrying `cases_month` -- drives the
#'                                        cases metric (see
#'                                        smoothing_events_mon_long()'s docstring
#'                                        in utils/smoothing_core.R).
#' Both come from the event tables rather than being re-derived from the raw
#' per-country filter state (the mon mean/dispersion state slots `mu`/`size_phi`
#' are filter state, not the per-month case draw).
load_smoothing_events <- function(ctx, parts_dir, filt_start) {
  read_smoothing_parts_events(
    parts_dir,
    filt_start,
    ctx$RES$series,
    tag = "[09]"
  )
} # utils/results_core.R

#' Weighted-quantile ribbon of `metric` (cases: sum(cases_month) over the
#' events_mon table; onsets: count of onset-event rows) over a scope's event
#' subset, zero-filled onto `all_months_vec` for every rep_id in `rep_ids_all`
#' (not just rep_ids present in this scope's subset -- a scope/month with zero
#' events must still contribute a zero).
ribbon_from_events <- function(
  events_scope,
  rep_ids_all,
  rep_weights,
  all_months_vec,
  metric
) {
  # zero-filled (rep_id, t_start) totals: rep_month_totals()
  agg_full <- rep_month_totals(
    events_scope,
    rep_ids_all,
    all_months_vec,
    metric
  )
  agg_full <- dplyr::left_join(agg_full, rep_weights, by = "rep_id")
  # A scope with no events (no months) has an empty ribbon; dplyr would still
  # evaluate wquantile() once on zero-length input, which has no quantile.
  if (!nrow(agg_full)) {
    return(tibble::tibble(
      t_start = as.Date(character(0)),
      q025 = numeric(0),
      q10 = numeric(0),
      q25 = numeric(0),
      q50 = numeric(0),
      q75 = numeric(0),
      q90 = numeric(0),
      q975 = numeric(0)
    ))
  }
  agg_full %>%
    dplyr::group_by(t_start) %>%
    dplyr::summarise(
      q025 = wquantile(value, w, 0.025),
      q10 = wquantile(value, w, 0.10),
      q25 = wquantile(value, w, 0.25),
      q50 = wquantile(value, w, 0.50),
      q75 = wquantile(value, w, 0.75),
      q90 = wquantile(value, w, 0.90),
      q975 = wquantile(value, w, 0.975),
      .groups = "drop"
    )
}

# Weighted quantiles: base.R::wquantile() (shared with smoothing_core.R).

.mask_to_coverage <- function(df, coverage_mask) {
  if (is.null(coverage_mask)) {
    return(df)
  }
  df %>%
    dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::semi_join(coverage_mask, by = c("gid", "t_start" = "month"))
}

derive_smoothing <- function(ctx) {
  message(
    "[09] smoothing: loading cases/onsets directly from part.qs (same source as the pipeline's own smoothing figures)"
  )
  man <- ctx$MAN
  sm <- load_smoothing_scope(ctx, "fig2-smoothing") # included groups, frame, grouping, event tables
  groups <- sm$groups
  surv_dat <- sm$surv_dat
  frame <- sm$frame
  enso_grp <- sm$enso_grp
  events_all <- sm$events_all
  events_mon_all <- sm$events_mon_all
  message(sprintf(
    "[09] smoothing: %d countries across %d included groups (%d onset rows, %d mon rows)",
    length(unique(events_all$country)),
    length(groups),
    nrow(events_all),
    nrow(events_mon_all)
  ))

  rep_ids_seen <- sort(unique(c(events_all$rep_id, events_mon_all$rep_id)))
  # Uniform replicate weights: importance weighting on this dataset provably
  # collapses to ESS ~= 1 over a 14-year, 44-country joint filter run (verified
  # against the pipeline's own merged ribbons); 08c_smoothing_compile.R itself
  # exposes --no_isweight for exactly this failure mode.
  rep_weights <- tibble::tibble(rep_id = rep_ids_seen) %>%
    dplyr::mutate(w = 1 / dplyr::n())

  surv_panel_obs <- build_obs_panel(surv_dat)

  # coverage_raw: every surveillance-covered unit-month (the validation
  # tables' denominator); coverage_mask: its at-risk cells (the ribbons').
  coverage_raw <- results_coverage_mask_raw(ctx, surv_dat)
  coverage_mask <- coverage_raw
  if (!is.null(coverage_mask)) {
    coverage_mask <- coverage_mask %>%
      dplyr::inner_join(
        dplyr::transmute(surv_panel_obs, gid, month = t_start),
        by = c("gid", "month")
      ) %>%
      dplyr::distinct(gid, month)
  }

  events_masked <- .mask_to_coverage(events_all, coverage_mask)
  events_mon_masked <- .mask_to_coverage(events_mon_all, coverage_mask)
  # Observed cases on the mon track's own estimand: cases_month at `month`, one
  # row per active month -- NOT a whole-episode total lumped onto floor_date(TL).
  # Only the rows of the admin level the filter counts per (gid, month)
  # (footprint.R::finest_level_rows()): summing every level double-counted nested
  # reports under lp_season; one level per gid (gid_season) is unchanged.
  # mon_all keeps the active months whose cases are unknown (the validation
  # tables' observed activity); mon_dat is the case series, as before.
  mon_all <- readRDS(ctx$MON_DATA_PATH) %>%
    dplyr::mutate(
      gid = as.character(gid),
      t_start = as.Date(month),
      obs_cases = as.numeric(cases_month)
    ) %>%
    dplyr::filter(!is.na(t_start), t_start >= ctx$SMOOTH_START)
  mon_dat <- mon_all %>%
    dplyr::filter(!is.na(obs_cases)) %>%
    finest_level_rows("obs_cases", c("gid", "t_start"))
  mon_dat_masked <- .mask_to_coverage(mon_dat, coverage_mask)
  surv_obs_masked <- .mask_to_coverage(
    dplyr::rename(surv_panel_obs, t_start = t_start),
    coverage_mask
  )

  scope_gids <- smoothing_scope_gids(
    frame,
    events_all,
    enso_grp,
    groups,
    ctx$GROUP_SHORT_RUN
  )

  # Each metric reads its own event table: onsets from the Yout==1 table, cases
  # from the mon table (active months, so a strictly wider month grid).
  ev_by_metric <- list(
    cases = list(all = events_mon_all, masked = events_mon_masked),
    onsets = list(all = events_all, masked = events_masked)
  )

  # One scope's ribbons (coverage-masked and unmasked) + observed series, per
  # metric. Month grids: `full_grid = FALSE` (Africa / group scopes, as before)
  # uses the months with at least one simulated event; `full_grid = TRUE`
  # (country scopes, which are sparse) uses every month from SMOOTH_START for the
  # unmasked series and every surveillance-covered month of the scope's gids for
  # the masked one, zero-filled, with the observed series zero-filled there too.
  all_months <- seq(
    ctx$SMOOTH_START,
    max(c(events_all$t_start, events_mon_all$t_start), na.rm = TRUE),
    by = "month"
  )
  scope_ribbons <- function(scope, gids, full_grid = FALSE) {
    out <- list()
    covered <- if (full_grid) {
      if (is.null(coverage_mask)) {
        all_months
      } else {
        sort(unique(coverage_mask$month[
          coverage_mask$gid %in% gids & coverage_mask$month >= ctx$SMOOTH_START
        ]))
      }
    }
    for (metric in c("cases", "onsets")) {
      ev_all_s <- dplyr::filter(ev_by_metric[[metric]]$all, gid %in% gids)
      ev_masked_s <- dplyr::filter(ev_by_metric[[metric]]$masked, gid %in% gids)
      months_unmasked <- if (full_grid) {
        all_months
      } else {
        sort(unique(ev_all_s$t_start))
      }
      months_masked <- if (full_grid) {
        covered
      } else {
        sort(unique(ev_masked_s$t_start))
      }
      obs_s <- if (metric == "cases") {
        mon_dat_masked %>%
          dplyr::filter(gid %in% gids, t_start %in% months_masked) %>%
          dplyr::group_by(t_start) %>%
          dplyr::summarise(
            value = sum(obs_cases, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        surv_obs_masked %>%
          dplyr::filter(
            gid %in% gids,
            t_start %in% months_masked,
            obs_Y == 1
          ) %>%
          dplyr::group_by(t_start) %>%
          dplyr::summarise(value = dplyr::n(), .groups = "drop")
      }
      if (full_grid) {
        obs_s <- dplyr::left_join(
          tibble::tibble(t_start = as.Date(months_masked)),
          obs_s,
          by = "t_start"
        ) %>%
          dplyr::mutate(value = dplyr::coalesce(as.numeric(value), 0))
      }
      out[[paste(scope, metric, sep = "|")]] <- list(
        scope = scope,
        metric = metric,
        ribbon_unmasked = ribbon_from_events(
          ev_all_s,
          rep_ids_seen,
          rep_weights,
          months_unmasked,
          metric
        ),
        ribbon_masked = ribbon_from_events(
          ev_masked_s,
          rep_ids_seen,
          rep_weights,
          months_masked,
          metric
        ),
        obs = obs_s
      )
    }
    out
  }

  results <- list()
  for (scope in names(scope_gids)) {
    results <- c(results, scope_ribbons(scope, scope_gids[[scope]]))
  }

  # Per-country ribbons (supplementary smoothing-by-country figure; the
  # 08c_smoothing_compile.R fig3 series): every country of an included group
  # with smoothing events, scope = ISO3, keyed "<ISO3>|<metric>".
  countries <- sort(unique(as.character(events_all$country)))
  by_country <- list()
  for (iso in countries) {
    by_country <- c(
      by_country,
      scope_ribbons(
        iso,
        unique(as.character(frame$gid[frame$country == iso])),
        full_grid = TRUE
      )
    )
  }
  country_group <- stats::setNames(
    enso_grp$group[match(countries, enso_grp$iso3)],
    countries
  )

  # --- v4: checks of the reconstruction against the observations (fig. SV) ---
  # Denominator: every surveillance-covered unit-month of the scope's gids from
  # SMOOTH_START to the last smoothed month (without a coverage mask: the
  # panel's rows plus the observed active months); at-risk = a survival-panel row.
  last_month <- max(c(events_all$t_start, events_mon_all$t_start), na.rm = TRUE)
  obs_mon <- dplyr::distinct(mon_all, gid, t_start) %>%
    dplyr::left_join(
      mon_dat %>%
        dplyr::group_by(gid, t_start) %>%
        dplyr::summarise(obs_cases = sum(obs_cases), .groups = "drop"),
      by = c("gid", "t_start")
    ) %>%
    dplyr::filter(t_start <= last_month)
  cells <- if (is.null(coverage_raw)) {
    dplyr::bind_rows(
      dplyr::select(surv_panel_obs, gid, t_start),
      dplyr::select(obs_mon, gid, t_start)
    )
  } else {
    dplyr::transmute(coverage_raw, gid, t_start = month)
  }
  cells <- cells %>%
    dplyr::filter(
      gid %in% scope_gids$Africa,
      t_start >= ctx$SMOOTH_START,
      t_start <= last_month
    ) %>%
    dplyr::distinct(gid, t_start) %>%
    dplyr::left_join(
      dplyr::select(surv_panel_obs, gid, t_start, obs_onset = obs_Y),
      by = c("gid", "t_start")
    )
  validation <- smoothing_validation_tables(
    dplyr::select(events_all, rep_id, gid, t_start),
    dplyr::select(events_mon_all, rep_id, gid, t_start, cases_month),
    rep_ids_seen,
    frame,
    cells,
    obs_mon,
    scope_gids
  )

  list(
    by_scope_metric = results,
    by_country = by_country,
    country_group = country_group,
    groups = groups,
    validation = validation
  )
}

#' smoothing_validation_tables
#'
#' Checks of the historical reconstruction against the observations (smoothing
#' v4; supplementary fig. SV), every replicate at uniform weight (as the
#' ribbons), from the smoothing event tables and the observed panel:
#'   pit_hist   randomised PIT of the observed monthly onsets (binary, at-risk
#'              unit-months: p = share of replicates with an onset) and cases
#'              (covered unit-months with known cases: the replicates' case
#'              draws, 0 when inactive), counted in `n_bins` equal bins per
#'              scope x metric (a storage resolution; 10 re-bins for display);
#'   seasonal   onset frequency per calendar month (onsets / at-risk
#'              unit-months), observed and across replicates;
#'   duration   share of episodes (runs of consecutive active months in a
#'              unit) by duration in months, observed and across replicates;
#'   active     fraction of covered unit-months in an active episode;
#'   contagion  onset frequency among at-risk unit-months not active themselves
#'              in the previous month, by whether an outbreak was ongoing in the
#'              previous month in the same ADM1 (`adm1`), in another ADM1 of the
#'              country (`country`) or in another country of the transmission
#'              unit (`tu`; units without a TU left out) -- the definitions of
#'              the kernel's update_lags(); state "ongoing" / "none".
#'   season_onsets  share of unit-seasons (gid x cholera year from min_seas,
#'              those with a covered cell) by number of onsets (0, 1, 2, 3 =
#'              three or more), observed and across replicates; with
#'              max_season_onsets (the most onsets in one unit-season of any
#'              replicate, every simulated onset) and n_unit_seasons. NULL / NA
#'              when the frame has no min_seas.
#' Replicate summaries: q025 / q25 / q50 / q75 / q975 of the per-replicate
#' statistic; `obs` the observed one. Scopes as the ribbons (Africa + the
#' groups' stored short names).
#' @param ev_onset tibble(rep_id, gid, t_start): one row per simulated onset
#' @param ev_mon tibble(rep_id, gid, t_start, cases_month): simulated active months
#' @param rep_ids every replicate id (one without events counts as zeros)
#' @param frame location frame (gid, adm1_unit, country, tu)
#' @param cells tibble(gid, t_start, obs_onset): covered unit-months, obs_onset
#'   0 / 1 on at-risk rows, NA elsewhere
#' @param obs_mon tibble(gid, t_start, obs_cases): observed active months (any
#'   coverage) with their cases (NA = unknown)
#' @param scopes named list scope -> gids (smoothing_scope_gids())
#' @param seed seed of the randomised PIT
#' @param n_bins PIT bins stored
smoothing_validation_tables <- function(
  ev_onset,
  ev_mon,
  rep_ids,
  frame,
  cells,
  obs_mon,
  scopes,
  seed = 42L,
  n_bins = 100L
) {
  fr <- frame %>%
    dplyr::transmute(
      gid = as.character(gid),
      adm1_unit = as.character(adm1_unit),
      country = as.character(country),
      tu = as.character(tu)
    ) %>%
    dplyr::distinct(gid, .keep_all = TRUE)
  gids <- fr$gid
  N <- length(gids)
  keep <- function(d) {
    dplyr::mutate(d, gid = as.character(gid), t_start = as.Date(t_start)) %>%
      dplyr::filter(gid %in% gids)
  }
  ev_onset <- keep(ev_onset)
  ev_mon <- keep(ev_mon)
  obs_mon <- keep(obs_mon)
  cells <- keep(cells) %>% dplyr::distinct(gid, t_start, .keep_all = TRUE)
  if (!nrow(cells)) {
    stop(
      "[09] smoothing validation: no covered unit-month in the location frame.",
      call. = FALSE
    )
  }
  m0 <- min(.month_index(c(
    cells$t_start,
    ev_onset$t_start,
    ev_mon$t_start,
    obs_mon$t_start
  )))
  M <- max(.month_index(c(
    cells$t_start,
    ev_onset$t_start,
    ev_mon$t_start,
    obs_mon$t_start
  ))) -
    m0 +
    1L
  lin <- function(g, d) match(g, gids) + (.month_index(d) - m0) * N
  qrow <- function(x) {
    q <- if (all(is.na(x))) {
      rep(NA_real_, 5)
    } else {
      stats::quantile(
        x,
        c(.025, .25, .5, .75, .975),
        na.rm = TRUE,
        names = FALSE
      )
    }
    tibble::tibble(q025 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q975 = q[5])
  }
  R <- length(rep_ids)
  S <- length(scopes)
  sc <- names(scopes)

  # Cells: linear index, the previous month's, scopes, calendar month.
  c_lin <- lin(cells$gid, cells$t_start)
  c_row <- match(cells$gid, gids)
  c_col <- .month_index(cells$t_start) - m0 + 1L
  risk <- which(!is.na(cells$obs_onset))
  in_scope <- lapply(scopes, function(g) cells$gid %in% g)
  in_risk <- lapply(in_scope, function(x) x[risk])
  moy_r <- as.integer(format(cells$t_start[risk], "%m"))
  # Contagion keys of the at-risk cells: ADM1 / country / TU rows and the
  # previous month's column (none in the first month).
  adm1_lv <- unique(fr$adm1_unit)
  ct_lv <- unique(fr$country)
  tu_lv <- unique(stats::na.omit(fr$tu))
  adm1_ct <- fr$country[match(adm1_lv, fr$adm1_unit)]
  ct_tu <- fr$tu[match(ct_lv, fr$country)]
  r_a <- match(fr$adm1_unit[c_row[risk]], adm1_lv)
  r_c <- match(fr$country[c_row[risk]], ct_lv)
  r_t <- match(fr$tu[c_row[risk]], tu_lv)
  r_ok <- c_col[risk] > 1L
  r_p <- pmax(c_col[risk] - 1L, 1L)
  r_prev <- c_row[risk] + (r_p - 1L) * N
  # Observed cases of the covered cells: 0 when not active, NA when unknown.
  om_lin <- lin(obs_mon$gid, obs_mon$t_start)
  o_case <- obs_mon$obs_cases[match(c_lin, om_lin)]
  o_case[!(c_lin %in% om_lin)] <- 0
  case_sel <- which(!is.na(o_case))

  # Run lengths of consecutive active months per unit, counted per scope (S x M).
  run_lengths <- function(A) {
    out <- matrix(0, S, M)
    w <- which(A, arr.ind = TRUE)
    if (!nrow(w)) {
      return(out)
    }
    w <- w[order(w[, 1], w[, 2]), , drop = FALSE]
    brk <- c(TRUE, w[-1, 1] != w[-nrow(w), 1] | w[-1, 2] != w[-nrow(w), 2] + 1L)
    len <- tabulate(cumsum(brk))
    g_run <- gids[w[brk, 1]]
    for (k in seq_len(S)) {
      out[k, ] <- tabulate(len[g_run %in% scopes[[k]]], M)
    }
    out
  }
  # One replicate's (or the observation's) statistics from its active matrix A
  # (N x M) and its onset indicator on the cells.
  stats_of <- function(A, on_c) {
    Ai <- A
    storage.mode(Ai) <- "integer"
    ad <- rowsum(Ai, fr$adm1_unit, reorder = FALSE)[adm1_lv, , drop = FALSE] > 0
    ctn <- rowsum(ad * 1L, adm1_ct, reorder = FALSE)[ct_lv, , drop = FALSE]
    cta <- ctn > 0
    has_tu <- !is.na(ct_tu)
    tu_on <- if (length(tu_lv)) {
      tun <- rowsum(
        cta[has_tu, , drop = FALSE] * 1L,
        ct_tu[has_tu],
        reorder = FALSE
      )[tu_lv, , drop = FALSE]
      tun[cbind(ifelse(is.na(r_t), 1L, r_t), r_p)] - cta[cbind(r_c, r_p)] >= 1L
    } else {
      rep(NA, length(risk))
    }
    ind <- list(
      adm1 = ad[cbind(r_a, r_p)],
      country = (ctn[cbind(r_c, r_p)] - ad[cbind(r_a, r_p)]) >= 1L,
      tu = ifelse(is.na(r_t), NA, tu_on)
    )
    free <- r_ok & !A[r_prev]
    o_r <- on_c[risk]
    cont <- array(NA_real_, c(3L, S, 2L, 2L)) # indicator x scope x state (ongoing, none) x (onsets, cells)
    for (k in 1:3) {
      for (s in seq_len(S)) {
        x <- ind[[k]]
        base <- free & in_risk[[s]]
        on_n <- base & x %in% TRUE
        on_f <- base & x %in% FALSE
        cont[k, s, , ] <- c(
          sum(o_r[on_n]),
          sum(o_r[on_f]),
          sum(on_n),
          sum(on_f)
        )
      }
    }
    seas <- vapply(
      seq_len(S),
      function(s) {
        n <- tabulate(moy_r[in_risk[[s]]], 12L)
        ifelse(
          n > 0,
          tabulate(moy_r[in_risk[[s]] & o_r == 1L], 12L) / pmax(n, 1L),
          NA_real_
        )
      },
      numeric(12)
    )
    act_c <- A[c_lin]
    list(
      seasonal = seas,
      active = vapply(in_scope, function(x) mean(act_c[x]), numeric(1)),
      contagion = cont,
      runs = run_lengths(A)
    )
  }

  on_split <- split(
    lin(ev_onset$gid, ev_onset$t_start),
    factor(ev_onset$rep_id, levels = rep_ids)
  )
  act_split <- split(
    seq_len(nrow(ev_mon)),
    factor(ev_mon$rep_id, levels = rep_ids)
  )
  em_lin <- lin(ev_mon$gid, ev_mon$t_start)
  p_sum <- numeric(length(risk))
  n_lt <- n_eq <- numeric(length(case_sel))
  seas <- array(NA_real_, c(R, 12L, S))
  act <- matrix(NA_real_, R, S)
  runs <- array(0, c(R, S, M))
  cont <- array(NA_real_, c(R, 3L, S, 2L))
  y_case <- o_case[case_sel]
  for (r in seq_len(R)) {
    ia <- act_split[[r]]
    A <- matrix(FALSE, N, M)
    A[em_lin[ia]] <- TRUE
    on_c <- as.integer(c_lin %in% on_split[[r]])
    p_sum <- p_sum + on_c[risk]
    hit <- match(c_lin[case_sel], em_lin[ia])
    x <- ifelse(is.na(hit), 0, ev_mon$cases_month[ia][hit])
    x[is.na(x)] <- 0
    n_lt <- n_lt + (x < y_case)
    n_eq <- n_eq + (x == y_case)
    st <- stats_of(A, on_c)
    seas[r, , ] <- st$seasonal
    act[r, ] <- st$active
    runs[r, , ] <- st$runs
    cont[r, , , ] <- st$contagion[,,, 1L] / st$contagion[,,, 2L]
  }
  A_obs <- matrix(FALSE, N, M)
  A_obs[om_lin] <- TRUE
  on_obs <- as.integer(cells$obs_onset %in% 1L)
  so <- stats_of(A_obs, on_obs)

  # Randomised PIT (binary onsets; cases on the replicates' empirical
  # distribution), binned per scope.
  ph <- p_sum / max(R, 1L)
  y_on <- on_obs[risk]
  pits <- .with_seed(
    seed,
    list(
      on = ifelse(
        y_on == 0L,
        stats::runif(length(risk), 0, 1 - ph),
        stats::runif(length(risk), 1 - ph, 1)
      ),
      ca = stats::runif(
        length(case_sel),
        n_lt / max(R, 1L),
        (n_lt + n_eq) / max(R, 1L)
      )
    )
  )
  pit_on <- pits$on
  pit_ca <- pits$ca
  brks <- seq(0, 1, length.out = n_bins + 1L)
  hist_of <- function(pit, scope, metric) {
    tibble::tibble(
      scope = scope,
      metric = metric,
      bin = seq_len(n_bins),
      lo = brks[-length(brks)],
      hi = brks[-1L],
      n = tabulate(
        findInterval(pit, brks, rightmost.closed = TRUE, all.inside = TRUE),
        n_bins
      )
    )
  }
  pit_hist <- dplyr::bind_rows(lapply(seq_len(S), function(s) {
    dplyr::bind_rows(
      hist_of(pit_on[in_risk[[s]]], sc[s], "onsets"),
      hist_of(pit_ca[in_scope[[s]][case_sel]], sc[s], "cases")
    )
  }))

  seasonal <- dplyr::bind_rows(lapply(seq_len(S), function(s) {
    dplyr::bind_rows(lapply(1:12, function(m) {
      dplyr::bind_cols(
        tibble::tibble(
          scope = sc[s],
          month = m,
          obs = so$seasonal[m, s],
          n_cells = sum(in_risk[[s]] & moy_r == m)
        ),
        qrow(seas[, m, s])
      )
    }))
  }))
  active <- dplyr::bind_rows(lapply(seq_len(S), function(s) {
    dplyr::bind_cols(
      tibble::tibble(
        scope = sc[s],
        obs = unname(so$active[s]),
        n_cells = sum(in_scope[[s]])
      ),
      qrow(act[, s])
    )
  }))
  tot <- apply(runs, 3, sum)
  d_max <- max(c(which(tot > 0), which(colSums(so$runs) > 0), 1L))
  duration <- dplyr::bind_rows(lapply(seq_len(S), function(s) {
    m <- matrix(runs[, s, ], nrow = R)
    rs <- rowSums(m)
    sh <- m / ifelse(rs > 0, rs, NA_real_)
    ob <- so$runs[s, ]
    obs_sh <- if (sum(ob) > 0) ob / sum(ob) else rep(NA_real_, M)
    dplyr::bind_rows(lapply(seq_len(d_max), function(d) {
      dplyr::bind_cols(
        tibble::tibble(
          scope = sc[s],
          duration = d,
          obs = obs_sh[d],
          n_obs = ob[d]
        ),
        qrow(sh[, d])
      )
    }))
  }))
  ind_lv <- c("adm1", "country", "tu")
  st_lv <- c("ongoing", "none")
  contagion <- dplyr::bind_rows(lapply(1:3, function(k) {
    dplyr::bind_rows(lapply(seq_len(S), function(s) {
      dplyr::bind_rows(lapply(1:2, function(j) {
        dplyr::bind_cols(
          tibble::tibble(
            scope = sc[s],
            indicator = ind_lv[k],
            state = st_lv[j],
            obs = so$contagion[k, s, j, 1L] / so$contagion[k, s, j, 2L],
            n_obs = so$contagion[k, s, j, 2L]
          ),
          qrow(cont[, k, s, j])
        )
      }))
    }))
  }))

  # Onsets per unit-season: the occurrence fit keeps
  # one onset per unit and cholera year, which simulator.onset_eligibility:
  # season imposes on the onsets that reset u_days (the detected ones under
  # covariate_feed: detected, so read it with results.series: detected).
  # Unit-seasons: the gid x cholera year (from min_seas) of the covered cells;
  # a replicate's onsets anywhere in such a season count. NULL without min_seas.
  season_onsets <- NULL
  max_season_onsets <- NA_integer_
  n_unit_seasons <- NA_integer_
  if ("min_seas" %in% names(frame)) {
    ms <- as.integer(frame$min_seas[match(gids, as.character(frame$gid))])
    us_of <- function(g, d) {
      i <- match(g, gids)
      sprintf("%d_%d", i, (.month_index(d) - (ms[i] - 1L)) %/% 12L)
    }
    us_cells <- us_of(cells$gid, cells$t_start)
    U <- unique(us_cells)
    n_unit_seasons <- length(U)
    us_gid <- gids[as.integer(sub("_.*", "", U))]
    us_in <- lapply(scopes, function(g) us_gid %in% g)
    ev_us <- us_of(ev_onset$gid, ev_onset$t_start)
    max_season_onsets <- if (nrow(ev_onset)) {
      max(tabulate(match(
        paste(ev_onset$rep_id, ev_us),
        unique(paste(ev_onset$rep_id, ev_us))
      )))
    } else {
      0L
    }
    n_cap <- 3L # bins 0, 1, 2, 3 (= three or more)
    bin_counts <- function(cnt) {
      vapply(
        seq_len(S),
        function(s) {
          tabulate(pmin(cnt[us_in[[s]]], n_cap) + 1L, n_cap + 1L)
        },
        numeric(n_cap + 1L)
      )
    }
    ev_us_split <- split(
      match(ev_us, U),
      factor(ev_onset$rep_id, levels = rep_ids)
    )
    sh <- array(NA_real_, c(R, n_cap + 1L, S))
    n_us <- unname(vapply(us_in, sum, numeric(1)))
    for (r in seq_len(R)) {
      b <- bin_counts(tabulate(stats::na.omit(ev_us_split[[r]]), length(U)))
      sh[r, , ] <- sweep(b, 2L, ifelse(n_us > 0, n_us, NA_real_), "/")
    }
    ob <- bin_counts(tabulate(
      match(us_cells[cells$obs_onset %in% 1L], U),
      length(U)
    ))
    season_onsets <- dplyr::bind_rows(lapply(seq_len(S), function(s) {
      dplyr::bind_rows(lapply(0:n_cap, function(k) {
        dplyr::bind_cols(
          tibble::tibble(
            scope = sc[s],
            n_onsets = k,
            obs = if (n_us[s] > 0) ob[k + 1L, s] / n_us[s] else NA_real_,
            n_obs = ob[k + 1L, s],
            n_unit_seasons = n_us[s]
          ),
          qrow(sh[, k + 1L, s])
        )
      }))
    }))
  }
  list(
    pit_hist = pit_hist,
    seasonal = seasonal,
    duration = duration,
    active = active,
    contagion = contagion,
    season_onsets = season_onsets,
    max_season_onsets = max_season_onsets,
    n_unit_seasons = n_unit_seasons,
    n_reps = R,
    n_cells = length(c_lin),
    n_at_risk = length(risk),
    n_case_cells = length(case_sel),
    pit_bins = n_bins,
    pit_seed = seed
  )
}

ensure_smoothing <- function(ctx) {
  # make_national_zero_months_file(): read by derive_smoothing()'s build_coverage_mask().
  # The panel's LP member bridge (coverage.R::coverage_lp_bridge()) is an input
  # only when it exists, so a panel built before the sidecar keeps its stamp.
  lpm <- occ_panel_lp_members_path(ctx$SURV_PATH)
  # v4 (validation tables): the coverage mask's aggregated series (its
  # denominator) and the HMM cluster file (the frame's TU, the contagion
  # check's third level), when they exist.
  inputs <- results_inputs(
    files = c(
      ctx$PANEL_PATH,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$MON_DATA_PATH,
      ctx$GROUPING_PATH,
      ctx$MANIFEST_PATH,
      as.character(national_zero_months_file_of(ctx$opt)),
      if (file.exists(lpm)) lpm,
      Filter(file.exists, c(ctx$AGG_TS_PATH, ctx$AGG_ADM_TS_PATH, ctx$HMM_PATH))
    ),
    globs = smoothing_parts_glob(ctx)
  )

  ensure_artefact(ctx, "smoothing", inputs, function() derive_smoothing(ctx))
}

# ===========================================================================
# Section 4: config-faithful ENSO envelope (Figure 3 ONI inset)
# ===========================================================================
# This reconstruction of what the simulator actually consumed reads sd_correction /
# rebase / target_peak / n_replicates
# / seed straight from the deployed factual config, and derives sim_months the
# same way .load_artefacts_grouped() (utils/simulations.R) does: the later of
# the horizon start and the smoothing-handoff pre-scan, through the end of the
# projection climate panel actually pointed at by the config.

#' Earliest smoothing-handoff month across the per-country hist_filter_tu=*.qs
#' caches (mirrors analysis/07b_simulator_spatpomp.R:227-255): take the last
#' 8-digit token of each basename as %Y%m%d, then min().
earliest_smoothing_handoff <- function(smoothing_input_dir) {
  fs <- Sys.glob(file.path(smoothing_input_dir, "hist_filter_tu=*_Np=*.qs"))
  if (!length(fs)) {
    stop("[09] no hist_filter_tu=*.qs caches under ", smoothing_input_dir)
  }
  toks <- regmatches(
    basename(fs),
    regexpr("[0-9]{8}(?=\\.qs$)", basename(fs), perl = TRUE)
  )
  min(as.Date(toks, format = "%Y%m%d"))
}

derive_oni_envelope <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  g0 <- groups[1]
  fac_cfg <- yaml::read_yaml(group_config_path(man, g0, "forecast"))$simulator
  # Each arm's ENSO matrix starts at its own launch (or the earlier smoothing-
  # handoff pre-scan), as 07b ran it; the arms may launch at different months
  # (results_arm_window()). `handoff` (Figure 3's vline) is the attribution start.
  win <- results_arm_window(ctx, groups, sets, read_events = FALSE)
  smoothing_input <- fac_cfg$filter$smoothing_input
  sm_handoff <- if (!is.null(smoothing_input) && dir.exists(smoothing_input)) {
    earliest_smoothing_handoff(smoothing_input)
  } else {
    NULL
  }
  panel <- readRDS(require_input(
    fac_cfg$climate_panel,
    "projection climate panel"
  ))
  months_from <- function(launch) {
    enso_start <- if (is.null(sm_handoff)) launch else min(launch, sm_handoff)
    sort(unique(panel$t_start[panel$t_start >= enso_start]))
  }
  horizon_start <- win$arm_starts$factual
  sim_months <- months_from(horizon_start)
  cf_months <- months_from(win$arm_starts$counterfactual)

  n_replicates <- as.integer(fac_cfg$n_replicates %||% 100L)
  seed <- as.integer(fac_cfg$seed %||% 20260603L)
  rebase <- fac_cfg$enso$cfsv2_rebase %||% "seam"
  sd_correction <- isTRUE(fac_cfg$enso$cfsv2_sd_correction)
  target_peak <- fac_cfg$enso$cfsv2_target_peak %||% NULL
  calibration <- if (nzchar(fac_cfg$enso$cfsv2_calibration %||% "")) {
    fac_cfg$enso$cfsv2_calibration
  } else {
    NULL
  }
  calibration_target <- fac_cfg$enso$cfsv2_calibration_target %||% "ersst"
  # ENSO past the CFSv2 span as the forecast ran it (simulator.enso.continuation;
  # absent = neutral, the generator writes the block only when not neutral).
  cont <- fac_cfg$enso$continuation %||% list()
  cont_method <- cont[["method"]] %||% "neutral"
  message(sprintf(
    "[09] ONI inset (config-faithful): n_rep=%d seed=%d rebase=%s sd_correction=%s calibration=%s continuation=%s sim_months=%s..%s",
    n_replicates,
    seed,
    rebase,
    sd_correction,
    if (is.null(calibration)) {
      "none"
    } else {
      paste0(calibration, " (", calibration_target, ")")
    },
    cont_method,
    min(sim_months),
    max(sim_months)
  ))

  historical_oni <- fetch_enso_index(
    source = "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  )
  historical_n34 <- tryCatch(
    fetch_enso_index(
      source = "nino34",
      cache_dir = ctx$ENSO_CACHE_DIR,
      refresh = FALSE
    ),
    error = function(e) NULL
  )
  cfsv2_data <- fetch_enso_cfsv2_ensemble(
    cache_dir = ctx$ENSO_CACHE_DIR,
    init = "latest"
  )

  summ_arm <- function(mat, arm, aset, months = sim_months) {
    tibble::tibble(
      date = months,
      arm = arm,
      analog_set = aset,
      mean = rowMeans(mat),
      lo = apply(mat, 1, quantile, 0.05),
      hi = apply(mat, 1, quantile, 0.95),
      # 95% envelope (oni v2): the manuscript figures' uniform 50% + 95% bands.
      lo95 = apply(mat, 1, quantile, 0.025),
      hi95 = apply(mat, 1, quantile, 0.975)
    )
  }

  fac <- sample_enso_cfsv2_trajectories(
    cfsv2_data = cfsv2_data,
    n_rep = n_replicates,
    historical_oni = historical_oni,
    historical_nino34 = historical_n34,
    sim_months = sim_months,
    max_lag = 12L,
    rebase = rebase,
    sd_correction = sd_correction,
    target_peak = target_peak,
    calibration = calibration,
    calibration_target = calibration_target,
    seed = seed,
    beyond_cfsv2 = cont_method,
    ar1_phi = cont[["phi"]],
    ar1_fit_window = cont[["fit_window"]] %||% "1991_2020",
    legacy_seam = isTRUE(cont[["legacy_seam"]]),
    analog_library = cont[["library"]] %||%
      "generated_data/enso_continuation_library.qs",
    analog_k = cont[["k"]] %||% 30L,
    analog_blend_months = cont[["blend_months"]] %||% 3,
    analog_min_sep = cont[["min_sep"]] %||% 24L,
    analog_max_per_model = cont[["max_per_model"]] %||% Inf,
    analog_exclude_models = unlist(cont[["exclude_models"]]),
    analog_seam = cont[["seam"]] %||% "last_cfsv2"
  )
  ctf_env <- dplyr::bind_rows(lapply(sets, function(aset) {
    cf_cfg <- yaml::read_yaml(group_config_path(
      man,
      g0,
      "counterfactual",
      aset
    ))$simulator
    # The counterfactual block, resolved as the simulator resolves it
    # (utils/enso.R::lenssen_counterfactual_args(); anchor default = the
    # counterfactual's own horizon.start).
    lcf <- lenssen_counterfactual_args(
      cf_cfg$enso$counterfactual,
      win$arm_starts$counterfactual
    )
    # The config's own cache (a per-run cache_path), not ctx$LENSSEN_CACHE.
    ctf <- sample_enso_lenssen_scenarios(
      cache_path = lcf$cache_path,
      sim_months = cf_months,
      n_rep = n_replicates,
      historical_oni = historical_oni,
      max_lag = 12L,
      horizon_start = lcf$horizon_start,
      mean_offset = lcf$mean_offset,
      log_var_scaling = lcf$log_var_scaling,
      seed = seed,
      analog_set = aset
    )
    summ_arm(ctf$enso_lag0, "counterfactual", aset, months = cf_months)
  }))
  env <- dplyr::bind_rows(
    summ_arm(fac$enso_lag0, "forecast", NA_character_),
    ctf_env
  )
  obs_full <- historical_oni %>%
    dplyr::transmute(
      date = as.Date(sprintf("%d-%02d-01", year, month)),
      enso_index
    )
  obs <- obs_full %>%
    dplyr::filter(
      date >= min(sim_months, cf_months),
      date <= max(sim_months, cf_months)
    )
  # ENSO phases (oni v3) for Figures 1 and 3: NOAA episodes
  # (.classify_enso_episodes) on the observed ONI extended by the CFSv2 ensemble
  # mean up to the last month carrying CFSv2 information (the forecast arm's
  # mean is exactly 0 beyond the CFSv2 span), so an episode still running at
  # the end of the observed record is not right-censored to "neutral". One row
  # per run of one phase, TL / TR = first / last month -/+ 15 days; 10 clips
  # them to each panel's window.
  fc_mean <- env %>%
    dplyr::filter(
      arm == "forecast",
      date > max(obs_full$date),
      abs(mean) > 1e-9
    ) %>%
    dplyr::transmute(date, enso_index = mean)
  ser <- dplyr::bind_rows(obs_full, fc_mean) %>% dplyr::arrange(date)
  ser$phase <- classify_enso_phase(ser$enso_index, "episode")
  phases <- ser %>%
    dplyr::mutate(
      run = cumsum(phase != dplyr::lag(phase, default = phase[1]))
    ) %>%
    dplyr::group_by(run) %>%
    dplyr::summarise(
      phase = phase[1],
      TL = min(date) - 15,
      TR = max(date) + 15,
      .groups = "drop"
    ) %>%
    dplyr::select(-run)
  list(
    env = env,
    obs = obs,
    obs_full = obs_full,
    phases = phases,
    # How `phases` was made (additive field; 10 falls back
    # to its legacy caption when absent).
    phase_rule = enso_phase_definition("episode", "monthly_oni_cfsv2_ext"),
    handoff = win$attrib_start,
    analog_sets = sets,
    # Per-arm launches (oni v4); equal to `handoff` when the arms launch together.
    handoffs = list(
      forecast = win$arm_starts$factual,
      counterfactual = win$arm_starts$counterfactual
    ),
    .inputs = results_inputs(files = c(fac_cfg$climate_panel))
  )
}

ensure_oni <- function(ctx, sets) {
  man <- ctx$MAN
  g0 <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)[1]
  cfgs <- c(
    group_config_path(man, g0, "forecast"),
    vapply(
      sets,
      function(a) group_config_path(man, g0, "counterfactual", a),
      character(1)
    )
  )
  # The counterfactual Lenssen caches are inputs too: a rebuilt
  # cache re-derives the envelope.
  caches <- unique(vapply(
    cfgs[-1L],
    function(p) {
      lenssen_counterfactual_args(
        yaml::read_yaml(p)$simulator$enso$counterfactual,
        Sys.Date()
      )$cache_path
    },
    character(1)
  ))
  inputs <- results_inputs(files = c(ctx$MANIFEST_PATH, cfgs, caches))

  ensure_artefact(
    ctx,
    "oni",
    inputs,
    function() derive_oni_envelope(ctx, sets),
    sets = sets,
    strip_inputs = TRUE
  )
}

# ===========================================================================
# Section 4b: Phase-11 "key numbers" report derivations (R1-R4)
# ===========================================================================
# Four report_* artefacts, each a bundle of tidy, kable-ready tibbles read ONLY
# by analysis/11_manuscript_numbers.qmd. Derivation lives here (phase 09) under
# the same results_read_if_fresh -> derive -> results_write freshness contract
# as Figures 1-3; the notebook only displays. Numbers are stored NUMERIC (never
# pre-formatted): the notebook's formatter helpers enforce "reasonable digits"
# centrally, and modelling tables carry mean/median + 95% CrI columns while the
# simulation tables carry median + 95% (q025/q975) columns.

#' gid -> GADM tiling level, ported verbatim from
#' notebooks/descriptive_occurrence_size.Rmd: surv$gid/size$gid are ALREADY
#' flattened (dots -> "_" upstream in 03_prepare_*), so this counts underscores
#' after stripping a trailing "_<n>" suffix, NOT dots.
gid_level <- function(g) {
  n_sep <- stringr::str_count(stringr::str_remove(g, "_\\d+$"), "_")
  dplyr::case_when(
    n_sep == 0 ~ "ADM0",
    n_sep == 1 ~ "ADM1",
    n_sep == 2 ~ "ADM2",
    n_sep >= 3 ~ "ADM3",
    TRUE ~ NA_character_
  )
}

#' Scalar quantile with a finite-only guard and rounding, for the dataset
#' descriptive tables (returns NA_real_ on an all-missing vector).
.qm <- function(x, p, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  round(as.numeric(stats::quantile(x, p, names = FALSE, na.rm = TRUE)), digits)
}

#' Mean, median and 95% interval of a draw vector (one-row tibble).
.q_mean95 <- function(x) {
  tibble::tibble(
    mean = mean(x, na.rm = TRUE),
    med = stats::median(x, na.rm = TRUE),
    lo = stats::quantile(x, .025, na.rm = TRUE, names = FALSE),
    hi = stats::quantile(x, .975, na.rm = TRUE, names = FALSE)
  )
}

#' Calendar months in [start, end], both ends inclusive (2026-06 .. 2026-06 = 1).
window_n_months <- function(start, end) {
  start <- as.POSIXlt(as.Date(start))
  end <- as.POSIXlt(as.Date(end))
  as.integer((end$year - start$year) * 12L + (end$mon - start$mon) + 1L)
}

#' Population denominators of the SIMULATED gids per scope, at POP_REF_YEAR,
#' from the canonical gid population table (resolve_gid_population(), flattened
#' keys). The simulator treats an ADM1 gid and the ADM2 gids nested in it as
#' separate units producing separate cases, so their populations are summed
#' exactly as simulated; n_gid_adm1_overlap counts the ADM1 gids that also have
#' an ADM2 child in the set (for transparency, not corrected).
#'
#' @param gids tibble(gid, country, enso_group) (load_group_arm_runs()$gids)
#' @return tibble(level in continent/region/country, scope, enso_group, region,
#'   pop, n_gid, n_gid_pop_missing, n_gid_adm1_overlap)
build_sim_scope_pop <- function(
  gids,
  gid_pop,
  pop_ref_year,
  group_short = GROUP_SHORT,
  pop_year_weights = NULL,
  wpp = NULL
) {
  g <- gids %>% dplyr::mutate(gid_flat = flatten_gid(gid))
  if (is.null(pop_year_weights)) {
    prf <- resolve_gid_population(
      g$gid_flat,
      pop_ref_year,
      gid_pop,
      wpp = NULL,
      on_missing = "na"
    )
    g$pop <- prf$pop
    g$pop_missing <- prf$pop_tier == "missing"
  } else {
    # rate_pop_year = "window": month-weighted mean over the
    # years of rate_pop_weights(); missing when any year is.
    prs <- lapply(names(pop_year_weights), function(y) {
      resolve_gid_population(
        g$gid_flat,
        as.integer(y),
        gid_pop,
        wpp = wpp,
        on_missing = "na"
      )
    })
    g$pop <- Reduce(`+`, Map(function(p, w) w * p$pop, prs, pop_year_weights))
    g$pop_missing <- Reduce(
      `|`,
      lapply(prs, function(p) p$pop_tier == "missing")
    )
  }
  n_us <- stringr::str_count(g$gid_flat, "_")
  g$adm1_unit <- stringr::str_extract(g$gid_flat, "^[A-Z]{3}_[0-9]+")
  adm2_units <- unique(g$adm1_unit[n_us >= 3L])
  g$adm1_overlap <- n_us == 2L & g$adm1_unit %in% adm2_units
  summ <- function(df) {
    dplyr::summarise(
      df,
      pop = sum(pop, na.rm = TRUE),
      n_gid = dplyr::n(),
      n_gid_pop_missing = sum(pop_missing),
      n_gid_adm1_overlap = sum(adm1_overlap),
      .groups = "drop"
    )
  }
  dplyr::bind_rows(
    summ(g) %>%
      dplyr::mutate(
        level = "continent",
        scope = "Africa",
        enso_group = NA_character_
      ),
    summ(dplyr::group_by(g, enso_group)) %>%
      dplyr::mutate(level = "region", scope = unname(group_short[enso_group])),
    summ(dplyr::group_by(g, enso_group, country)) %>%
      dplyr::mutate(level = "country", scope = country)
  ) %>%
    dplyr::mutate(region = unname(group_short[enso_group])) %>%
    dplyr::select(
      level,
      scope,
      enso_group,
      region,
      pop,
      n_gid,
      n_gid_pop_missing,
      n_gid_adm1_overlap
    )
}

#' Per-replicate horizon totals of a (replicate_id, month, n_cases, n_outbreaks)
#' frame grouped by `key_cols` + arm, with the scope's simulated-month count
#' n_months and its population `pop` already joined as columns.
#' inc_ann_100k annualises the window total: tot_cases / n_months * 12 / pop * 1e5
#' (NA when the population is missing or zero).
per_rep_totals <- function(df, key_cols) {
  .sum_by_keys(
    df,
    c(key_cols, "run", "analog_set", "replicate_id", "n_months", "pop"),
    c(tot_cases = "n_cases", tot_onsets = "n_outbreaks")
  ) %>%
    dplyr::mutate(
      mean_cases = tot_cases / n_months,
      mean_onsets = tot_onsets / n_months,
      inc_ann_100k = dplyr::if_else(
        is.na(pop) | pop <= 0,
        NA_real_,
        tot_cases / n_months * 12 / pop * 1e5
      )
    )
}

#' median + 95% (q025/q975) across a vector of per-replicate values -- the
#' simulation-output reporting convention, one row.
.q_med95 <- function(x) {
  tibble::tibble(
    med = stats::median(x, na.rm = TRUE),
    lo = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    hi = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

# --- R1: dataset descriptors + Table 1 -------------------------------------

derive_report_dataset <- function(ctx) {
  for (p in c(
    ctx$SURV_PATH,
    ctx$SIZE_DATA_PATH,
    ctx$DUR_DATA_PATH,
    ctx$GID_POP_PATH
  )) {
    require_input(p, "raw output / population table")
  }

  surv <- readRDS(ctx$SURV_PATH) %>% dplyr::mutate(gid = as.character(gid))
  size <- readRDS(ctx$SIZE_DATA_PATH) %>%
    dplyr::filter(as.character(admin_level) != "ADM0") %>%
    dplyr::mutate(gid = as.character(gid))
  dur <- readRDS(ctx$DUR_DATA_PATH) %>%
    dplyr::filter(as.character(admin_level) != "ADM0")

  # --- headline counts (from the size dataset, the analysis outbreak set) ---
  counts <- tibble::tibble(
    metric = c("n_outbreaks", "n_countries", "n_gids", "n_suspected_cases"),
    value = c(
      nrow(size),
      dplyr::n_distinct(size$country),
      dplyr::n_distinct(size$gid),
      sum(size$suspected_cases, na.rm = TRUE)
    )
  )

  # --- unique gids by GADM tiling level ------------------------------------
  gids_by_level <- size %>%
    dplyr::distinct(gid) %>%
    dplyr::mutate(gid_level = gid_level(gid)) %>%
    dplyr::count(gid_level, name = "n_gids") %>%
    dplyr::arrange(gid_level)

  # --- Table 1 by (corpus) admin level, Zheng et al. layout ----------------
  size_lvl <- size %>%
    dplyr::mutate(
      admin_level = as.character(admin_level),
      duration_weeks = as.numeric(TR - TL) / 7,
      attack_rate_1k = ifelse(
        !is.na(pop) & pop > 0,
        1000 * suspected_cases / pop,
        NA_real_
      )
    )
  table1 <- size_lvl %>%
    dplyr::group_by(admin_level) %>%
    dplyr::summarise(
      n_outbreaks = dplyr::n(),
      n_gids = dplyr::n_distinct(gid),
      n_years = dplyr::n_distinct(year_seas),
      total_cases = sum(suspected_cases, na.rm = TRUE),
      cases_q25 = .qm(suspected_cases, .25, 0),
      cases_med = .qm(suspected_cases, .5, 0),
      cases_q75 = .qm(suspected_cases, .75, 0),
      dur_wk_q25 = .qm(duration_weeks, .25, 1),
      dur_wk_med = .qm(duration_weeks, .5, 1),
      dur_wk_q75 = .qm(duration_weeks, .75, 1),
      pop_q25 = .qm(pop, .25, 0),
      pop_med = .qm(pop, .5, 0),
      pop_q75 = .qm(pop, .75, 0),
      ar_q25 = .qm(attack_rate_1k, .25, 2),
      ar_med = .qm(attack_rate_1k, .5, 2),
      ar_q75 = .qm(attack_rate_1k, .75, 2),
      .groups = "drop"
    )
  # Episode-level duration (months) per admin level from the duration dataset.
  dur_lvl <- dur %>%
    dplyr::distinct(episode_id, admin_level, duration_months) %>%
    dplyr::mutate(admin_level = as.character(admin_level)) %>%
    dplyr::group_by(admin_level) %>%
    dplyr::summarise(
      n_episodes = dplyr::n_distinct(episode_id),
      dur_mo_q25 = .qm(duration_months, .25, 1),
      dur_mo_med = .qm(duration_months, .5, 1),
      dur_mo_q75 = .qm(duration_months, .75, 1),
      .groups = "drop"
    )
  table1 <- dplyr::left_join(table1, dur_lvl, by = "admin_level") %>%
    dplyr::arrange(admin_level)

  # --- temporal coverage + total population covered ------------------------
  # Total covered population is the sum of the OBSERVED sub-ADM1 populations
  # (build_adm1_obs_pop, the same map R2/R4 join on), with recover_dropped = TRUE
  # so the few ADM1 units the string-level overlap filter over-drops are re-added
  # from their observed ADM2 sub-units. This keeps the R1 "total population
  # covered" descriptor reconciled with the R2 effects-population ceiling.
  gid_pop <- load_gid_population(ctx$GID_POP_PATH)
  adm1_obs <- build_adm1_obs_pop(
    ctx,
    surv %>% dplyr::distinct(gid, country, AFRO_region),
    surv,
    gid_pop,
    "report-dataset",
    ctx$POP_REF_YEAR,
    recover_dropped = TRUE
  )
  coverage <- tibble::tibble(
    date_min = min(as.Date(surv$t_start), na.rm = TRUE),
    date_max = max(as.Date(surv$t_start), na.rm = TRUE),
    total_location_months = nrow(surv),
    n_gids_covered = sum(adm1_obs$n_gid),
    n_adm1_covered = nrow(adm1_obs),
    pop_ref_year = ctx$POP_REF_YEAR,
    total_pop_covered = sum(adm1_obs$pop, na.rm = TRUE),
    n_gids_pop_missing = sum(adm1_obs$n_pop_missing)
  )

  # --- table S2: canonical episodes by ENSO group x country x admin level ---
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  episodes_by_country <- size %>%
    dplyr::mutate(
      admin_level = as.character(admin_level),
      enso_group = reg$group[match(country, reg$iso3)]
    ) %>%
    dplyr::group_by(enso_group, country, admin_level) %>%
    dplyr::summarise(
      n_episodes = dplyr::n(),
      n_gids = dplyr::n_distinct(gid),
      total_cases = sum(suspected_cases, na.rm = TRUE),
      first_start = min(as.Date(TL), na.rm = TRUE),
      last_end = max(as.Date(TR), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(enso_group, country, admin_level)

  # --- v4: analysis panel counts, ADM0 excluded ----------------------------
  surv_sub <- if ("admin_level" %in% names(surv)) {
    dplyr::filter(surv, as.character(admin_level) != "ADM0")
  } else {
    surv
  }
  coverage <- dplyr::mutate(
    coverage,
    n_gids_panel = dplyr::n_distinct(surv_sub$gid),
    n_countries_panel = dplyr::n_distinct(surv_sub$country),
    n_location_months_panel = nrow(surv_sub)
  )

  # --- v4: raw outbreak set (01's outbreak metadata, ADM0 excluded) ---------
  raw <- readRDS(require_input(ctx$OB_META_PATH, "outbreak metadata (01)")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(admin_level = as.character(admin_level)) %>%
    dplyr::filter(admin_level != "ADM0")
  raw_counts <- tibble::tibble(
    n_outbreaks = nrow(raw),
    n_locations = dplyr::n_distinct(raw$location),
    n_countries = dplyr::n_distinct(raw$country),
    n_suspected_cases = sum(raw$suspected_cases, na.rm = TRUE)
  )
  raw_table1 <- .zheng_table(
    raw %>%
      dplyr::mutate(
        duration_weeks = as.numeric(as.Date(TR) - as.Date(TL)) / 7,
        attack_rate_1k = ifelse(
          !is.na(pop) & pop > 0,
          1000 * suspected_cases / pop,
          NA_real_
        )
      )
  )

  # --- v4: one row per analysis outbreak, with its ENSO state ---------------
  # ENSO state of the outbreak's cholera season-year: the ONI of the DJF whose
  # December lies in that season-year (a season-year y starts in month
  # min_seas of year y, so it always contains December y; DJF is centred on
  # January y + 1), threshold phase at +/-0.5 (seasonal, no persistence: not
  # the CPC episode rule). Rows feed 11's gtsummary tables.
  djf <- fetch_enso_index(
    source = "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  ) %>%
    dplyr::filter(as.integer(month) == 1L) %>%
    dplyr::transmute(year_seas = as.integer(year) - 1L, oni_djf = enso_index)
  dur_ep <- dur %>%
    dplyr::distinct(uid, duration_months) %>%
    dplyr::distinct(uid, .keep_all = TRUE)
  outbreaks <- size_lvl %>%
    dplyr::transmute(
      uid,
      gid,
      country,
      afro_subregion = afro_subregion_of(country)$afro_subregion,
      enso_group = reg$group[match(country, reg$iso3)],
      region = unname(ctx$GROUP_SHORT_RUN[enso_group]),
      admin_level,
      year_seas = as.integer(as.character(year_seas)),
      cases = suspected_cases,
      duration_weeks,
      pop,
      attack_rate_1k
    ) %>%
    dplyr::left_join(dur_ep, by = "uid") %>%
    dplyr::left_join(djf, by = "year_seas") %>%
    dplyr::mutate(
      enso_state = unname(ENSO_PHASE_LABELS[classify_enso_phase(
        oni_djf,
        "threshold"
      )])
    ) %>%
    dplyr::select(-uid)
  stopifnot(nrow(outbreaks) == nrow(size))

  # --- v5: blocks Africa + WHO-AFRO subregion (tables S2 and S3) -------------
  # table1_region: table1's descriptors by admin level and "All" levels, per
  # block. enso_region: outbreaks per season-year ENSO state (and "All") per
  # block, their share of the block, cholera years, cases, duration (weeks) and
  # attack rate, and outbreaks per cholera year of that state (n_state_years:
  # the record's cholera years in the state, the same for every block, so a
  # subregion's quiet years count); enso_region_p: Kruskal-Wallis p across the
  # three states.
  if (anyNA(outbreaks$afro_subregion)) {
    stop(
      "[09] report_dataset: no WHO-AFRO subregion for ",
      paste(
        unique(outbreaks$country[is.na(outbreaks$afro_subregion)]),
        collapse = ", "
      )
    )
  }
  blocks <- c(
    "Africa",
    intersect(AFRO_SUBREGION_ORDER, outbreaks$afro_subregion)
  )
  in_block <- function(b) {
    if (b == "Africa") outbreaks else outbreaks[outbreaks$afro_subregion == b, ]
  }
  set_summ <- function(x) {
    dplyr::summarise(
      x,
      n_outbreaks = dplyr::n(),
      n_gids = dplyr::n_distinct(gid),
      n_years = dplyr::n_distinct(year_seas),
      total_cases = sum(cases, na.rm = TRUE),
      cases_q25 = .qm(cases, .25, 0),
      cases_med = .qm(cases, .5, 0),
      cases_q75 = .qm(cases, .75, 0),
      dur_wk_q25 = .qm(duration_weeks, .25, 1),
      dur_wk_med = .qm(duration_weeks, .5, 1),
      dur_wk_q75 = .qm(duration_weeks, .75, 1),
      pop_q25 = .qm(pop, .25, 0),
      pop_med = .qm(pop, .5, 0),
      pop_q75 = .qm(pop, .75, 0),
      ar_q25 = .qm(attack_rate_1k, .25, 2),
      ar_med = .qm(attack_rate_1k, .5, 2),
      ar_q75 = .qm(attack_rate_1k, .75, 2)
    )
  }
  table1_region <- dplyr::bind_rows(lapply(blocks, function(b) {
    x <- in_block(b)
    dplyr::bind_rows(lapply(c(sort(unique(x$admin_level)), "All"), function(l) {
      dplyr::mutate(
        set_summ(if (l == "All") x else x[x$admin_level == l, ]),
        region = b,
        admin_level = l,
        .before = 1
      )
    }))
  }))
  states <- unname(ENSO_PHASE_LABELS[c("nina", "neutral", "elnino")])
  state_years <- c(
    vapply(
      states,
      function(st) {
        dplyr::n_distinct(outbreaks$year_seas[outbreaks$enso_state %in% st])
      },
      integer(1)
    ),
    All = dplyr::n_distinct(outbreaks$year_seas)
  )
  enso_region <- dplyr::bind_rows(lapply(blocks, function(b) {
    x <- in_block(b)
    dplyr::bind_rows(lapply(c(states, "All"), function(st) {
      y <- if (st == "All") x else x[x$enso_state %in% st, ]
      set_summ(y) %>%
        dplyr::select(-n_gids, -total_cases, -dplyr::starts_with("pop_")) %>%
        dplyr::mutate(
          region = b,
          enso_state = st,
          share = n_outbreaks / nrow(x),
          n_state_years = unname(state_years[st]),
          per_year = n_outbreaks / n_state_years,
          .before = 1
        )
    }))
  }))
  enso_region_p <- dplyr::bind_rows(lapply(blocks, function(b) {
    x <- in_block(b)
    dplyr::bind_rows(lapply(
      c("cases", "duration_weeks", "attack_rate_1k"),
      function(v) {
        tibble::tibble(
          region = b,
          variable = v,
          p_value = tryCatch(
            stats::kruskal.test(x[[v]], factor(x$enso_state, states))$p.value,
            error = function(e) NA_real_
          )
        )
      }
    ))
  }))

  list(
    counts = counts,
    gids_by_level = gids_by_level,
    table1 = table1,
    coverage = coverage,
    episodes_by_country = episodes_by_country,
    cascade = report_cascade(ctx),
    raw_counts = raw_counts,
    raw_table1 = raw_table1,
    outbreaks = outbreaks,
    table1_region = table1_region,
    enso_region = enso_region,
    enso_region_p = enso_region_p,
    # How outbreaks$enso_state was made (additive field).
    enso_state_rule = enso_phase_definition("threshold", "djf_oni"),
    .inputs = results_inputs(
      files = c(
        ctx$SURV_PATH,
        ctx$SIZE_DATA_PATH,
        ctx$DUR_DATA_PATH,
        ctx$GID_POP_PATH,
        ctx$OB_META_PATH,
        .oni_cache_file(ctx)
      )
    )
  )
}

#' The ONI cache fetch_enso_index("oni") reads (an input of report_dataset).
.oni_cache_file <- function(ctx) {
  file.path(ctx$ENSO_CACHE_DIR, "enso_oni_cache.rds")
}

#' Outbreak descriptors in the Zheng et al. (2022) layout, by admin level plus
#' an "All" row: counts, total cases, and median [IQR] of cases, duration
#' (weeks), population and attack rate (/1,000). `d` has admin_level, location
#' or gid, country, suspected_cases, duration_weeks, pop, attack_rate_1k.
.zheng_table <- function(d) {
  loc <- if ("location" %in% names(d)) "location" else "gid"
  summ <- function(x) {
    x %>%
      dplyr::summarise(
        n_outbreaks = dplyr::n(),
        n_locations = dplyr::n_distinct(.data[[loc]]),
        n_countries = dplyr::n_distinct(country),
        total_cases = sum(suspected_cases, na.rm = TRUE),
        cases_q25 = .qm(suspected_cases, .25, 0),
        cases_med = .qm(suspected_cases, .5, 0),
        cases_q75 = .qm(suspected_cases, .75, 0),
        dur_wk_q25 = .qm(duration_weeks, .25, 1),
        dur_wk_med = .qm(duration_weeks, .5, 1),
        dur_wk_q75 = .qm(duration_weeks, .75, 1),
        pop_q25 = .qm(pop, .25, 0),
        pop_med = .qm(pop, .5, 0),
        pop_q75 = .qm(pop, .75, 0),
        ar_q25 = .qm(attack_rate_1k, .25, 2),
        ar_med = .qm(attack_rate_1k, .5, 2),
        ar_q75 = .qm(attack_rate_1k, .75, 2),
        .groups = "drop"
      )
  }
  dplyr::bind_rows(
    summ(dplyr::group_by(d, admin_level)) %>% dplyr::arrange(admin_level),
    dplyr::mutate(summ(d), admin_level = "All")
  )
}

#' The processing-step count files and prep stamps the cascade reads.
#' @noRd
.CASCADE_STEPS <- c(
  "00_ingest_outbreakextractr",
  "01_preprocess_outbreaks",
  "02_preprocess_crosswalk"
)
.cascade_files <- function(ctx) {
  c(
    step_counts_path(.CASCADE_STEPS, ctx$PREP_COUNTS_DIR),
    vapply(
      c(
        ctx$SIZE_DATA_PATH,
        ctx$SURV_PATH,
        ctx$DUR_DATA_PATH,
        ctx$MON_DATA_PATH
      ),
      prep_policy_stamp_path,
      character(1)
    )
  )
}

#' report_cascade
#'
#' Table S1: counts at each processing step, as a long tibble (step, metric,
#' value) in pipeline order. The 00 / 01 / 02 counts come from their
#' generated_data/prep_counts/<step>.yml files (a missing file is reported with
#' a warning), the 03 counts from
#' the size / occurrence / duration / monthly prep stamps (numeric scalars
#' n_*, cases_*, median_*, span_*; span_ carries the duration panel's episode
#' span IQR).
report_cascade <- function(ctx) {
  rows <- list()
  for (st in .CASCADE_STEPS) {
    cnt <- read_step_counts(st, ctx$PREP_COUNTS_DIR)
    if (is.null(cnt)) {
      warning(sprintf(
        "[09] report_dataset: no processing counts for %s (%s); re-run that step to record them.",
        st,
        step_counts_path(st, ctx$PREP_COUNTS_DIR)
      ))
      next
    }
    cnt <- cnt[setdiff(names(cnt), c("step", "written"))]
    rows[[st]] <- tibble::tibble(
      step = st,
      metric = names(cnt),
      value = as.numeric(unlist(cnt))
    )
  }
  stamps <- c(
    `03_prepare_siz` = ctx$SIZE_DATA_PATH,
    `03_prepare_occ` = ctx$SURV_PATH,
    `03_prepare_dur` = ctx$DUR_DATA_PATH,
    `03_prepare_mon` = ctx$MON_DATA_PATH
  )
  for (st in names(stamps)) {
    s <- read_prep_policy_stamp(stamps[[st]])
    if (is.null(s)) {
      warning("[09] report_dataset: no prep stamp for ", stamps[[st]])
      next
    }
    keep <- names(s)[
      grepl("^(n_|cases_|median_|span_)", names(s)) &
        vapply(s, function(x) is.numeric(x) && length(x) == 1L, logical(1))
    ]
    rows[[st]] <- tibble::tibble(
      step = st,
      metric = keep,
      value = as.numeric(unlist(s[keep]))
    )
  }
  dplyr::bind_rows(rows)
}

ensure_report_dataset <- function(ctx) {
  inputs <- results_inputs(
    files = c(
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$DUR_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH,
      ctx$OB_META_PATH,
      .oni_cache_file(ctx),
      .cascade_files(ctx)
    )
  )
  ensure_artefact(
    ctx,
    "report_dataset",
    inputs,
    function() derive_report_dataset(ctx),
    strip_inputs = TRUE
  )
}

# --- R2: BRMS effect sizes per track ---------------------------------------

#' Construction-faithful duration RMD newdata: one row per months-since-onset
#' 1..horizon, taken from the fit's OWN design frame (so mso_band / the mso &
#' season M-spline columns are the actual fitted values, never re-derived), with
#' only the ENSO cross-basis columns overwritten to the onebasis row at a
#' constant ONI value via build_cb_columns_from_Q() (bit-identical to fit time).
#' Rows are keyed by `mso` when the fit frame carries it, else by the mso->band
#' map from the duration dataset.
build_dur_newdata <- function(fit, vm, oni_value, horizon, mso_to_band) {
  d <- fit$data
  cb_names <- vm$cb_col_names
  pick_row <- function(m) {
    if ("mso" %in% names(d)) {
      r <- d[d$mso == m, , drop = FALSE]
      if (!nrow(r)) {
        r <- d[order(abs(d$mso - m)), , drop = FALSE][1, , drop = FALSE]
      }
    } else if ("mso_band" %in% names(d)) {
      band <- mso_to_band[[as.character(m)]]
      r <- if (!is.null(band)) {
        d[as.character(d$mso_band) == band, , drop = FALSE]
      } else {
        d[0, , drop = FALSE]
      }
      if (!nrow(r)) r <- d[1, , drop = FALSE]
    } else {
      r <- d[1, , drop = FALSE]
    }
    r[1, , drop = FALSE]
  }
  nd <- do.call(rbind, lapply(seq_len(horizon), pick_row))
  z <- build_cb_columns_from_Q(oni_value, vm) # 1 x K, raw onebasis
  for (j in seq_along(cb_names)) {
    nd[[cb_names[j]]] <- as.numeric(z[1, j])
  }
  nd
}

DUR_HARMONIC_COLS <- c("sin_month", "cos_month", "sin2_month", "cos2_month")

#' Does a duration fit carry duration.seasonal_term's calendar harmonics?
dur_fit_has_harmonics <- function(fit) {
  all(DUR_HARMONIC_COLS %in% names(fit$data))
}

#' Overwrite the harmonic columns of a build_dur_newdata() frame so row m (mso m)
#' sits in calendar month ((onset_month - 1) + (m - 1)) mod 12 + 1: an episode
#' starting in `onset_month`, the month advancing with months since onset (the
#' construction load_duration_data() uses on the panel's calendar month).
set_dur_onset_month <- function(nd, onset_month) {
  mo <- ((onset_month - 1L) + seq_len(nrow(nd)) - 1L) %% 12L + 1L
  nd$sin_month <- sin(2 * pi * mo / 12)
  nd$cos_month <- cos(2 * pi * mo / 12)
  nd$sin2_month <- sin(4 * pi * mo / 12)
  nd$cos2_month <- cos(4 * pi * mo / 12)
  nd
}

#' duration_ratio_draws() for one ONI contrast. A fit without the harmonics is
#' scored once on build_dur_newdata() (unchanged). A fit with them is scored for
#' each of the 12 onset months; the restricted mean durations are averaged with
#' equal weight per onset month (a typical episode, whatever month it starts) and the
#' ratio is taken per draw on those averages.
#' Returns the duration_ratio_draws() list; tail_mass is the largest over onset
#' months, warned about once.
dur_rmd_ratio_draws <- function(
  fit,
  vm,
  oni_value,
  horizon,
  mso_to_band,
  re_formula = NA,
  tail_tol = 0.05
) {
  nd_hi <- build_dur_newdata(fit, vm, oni_value, horizon, mso_to_band)
  nd_lo <- build_dur_newdata(fit, vm, 0, horizon, mso_to_band)
  if (!dur_fit_has_harmonics(fit)) {
    return(duration_ratio_draws(
      fit,
      newdata_hi = nd_hi,
      newdata_lo = nd_lo,
      horizon = horizon,
      re_formula = re_formula,
      tail_tol = tail_tol
    ))
  }
  per <- lapply(1:12, function(s) {
    withCallingHandlers(
      duration_ratio_draws(
        fit,
        newdata_hi = set_dur_onset_month(nd_hi, s),
        newdata_lo = set_dur_onset_month(nd_lo, s),
        horizon = horizon,
        re_formula = re_formula,
        tail_tol = tail_tol
      ),
      warning = function(w) {
        if (startsWith(conditionMessage(w), "duration_ratio_draws: median S")) {
          invokeRestart("muffleWarning")
        }
      }
    )
  })
  rmd_hi <- Reduce(`+`, lapply(per, `[[`, "rmd_hi")) / 12
  rmd_lo <- Reduce(`+`, lapply(per, `[[`, "rmd_lo")) / 12
  tail_mass <- max(vapply(per, `[[`, numeric(1), "tail_mass"))
  if (tail_mass > tail_tol) {
    warning(sprintf(
      paste(
        "dur_rmd_ratio_draws: median S(%d) = %.3f > %.3f for at least one onset month --",
        "report this as a %d-month RESTRICTED mean duration ratio."
      ),
      horizon,
      tail_mass,
      tail_tol,
      horizon
    ))
  }
  list(
    ratio = rmd_hi / rmd_lo,
    rmd_hi = rmd_hi,
    rmd_lo = rmd_lo,
    horizon = horizon,
    tail_mass = tail_mass
  )
}

derive_report_effects <- function(ctx, track) {
  enso_rel <- ensure_enso_relation(ctx, track) # reuse the ENSO exposure-response extraction as compute source
  groups <- enso_rel$groups
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  gid_pop <- load_gid_population(ctx$GID_POP_PATH)

  # --- observed sub-ADM1 population map (shared with R1/R4) -----------------
  # Population credited to each ADM1 unit is the sum of the sub-units ACTUALLY
  # observed within it, NOT the whole province (which double-counts geography and
  # inflated every effects-population total to ~769M, above R1's ~547M coverage).
  # The occurrence survival dataset is the single containment universe for every
  # track; recover_dropped = TRUE re-adds the few units the string-level
  # overlap filter over-drops, matching the recovered R1 (Table-4) total.
  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  pop_by_adm1 <- obs_adm1_pop_map(ctx, surv_dat, gid_pop, "report-effects")

  # --- best lag deployed per group -----------------------------------------
  best_lag <- dplyr::bind_rows(lapply(groups, function(g) {
    tibble::tibble(
      enso_group = g,
      region = unname(ctx$GROUP_SHORT_RUN[g]),
      best_lag = as.integer(track_fit_files(mnf, track, g)$lag)
    )
  }))

  # nearest-grid ERF row to each of the 4 composite ONI conditions
  conds <- ctx$ENSO_CONDS
  pick_cond <- function(df) {
    dplyr::bind_rows(lapply(seq_along(conds), function(i) {
      j <- which.min(abs(df$at - conds[[i]]))
      dplyr::mutate(
        df[j, , drop = FALSE],
        enso_label = names(conds)[i],
        enso_val = unname(conds[[i]])
      )
    }))
  }

  effects_region <- enso_rel$erf %>%
    dplyr::group_by(enso_group) %>%
    dplyr::group_modify(~ pick_cond(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      enso_group,
      region = unname(ctx$GROUP_SHORT_RUN[enso_group]),
      enso_label,
      enso_val,
      IRR,
      IRR_lo,
      IRR_hi,
      # v6: the fitted ONI range of the group's crossbasis
      # (enso_relation; boundary knots), so 11 can flag a
      # condition outside it.
      dplyr::across(dplyr::any_of(c("range_lo", "range_hi")))
    )

  effects_country <- if (
    !is.null(enso_rel$country_erf) && nrow(enso_rel$country_erf)
  ) {
    enso_rel$country_erf %>%
      dplyr::group_by(enso_group, country) %>%
      dplyr::group_modify(~ pick_cond(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(
        enso_group,
        region = unname(ctx$GROUP_SHORT_RUN[enso_group]),
        country,
        enso_label,
        enso_val,
        IRR,
        IRR_lo,
        IRR_hi,
        dplyr::across(dplyr::any_of(c("range_lo", "range_hi")))
      )
  } else {
    tibble::tibble()
  }

  # --- population with a significant per-ADM1 ENSO effect -------------------
  # enso_rel$adm1_irr$lo/hi are on the LOG scale: significant increase = lo > 0,
  # decrease = hi < 0. adm1_gid is the flattened gid == the adm1_unit key of the
  # observed-sub-unit map, so each significant ADM1 is credited only the
  # population of the districts observed within it (never the whole province).
  adm1 <- enso_rel$adm1_irr %>% dplyr::mutate(adm1_gid = as.character(adm1_gid))
  adm1$pop <- unname(pop_by_adm1[adm1$adm1_gid])
  adm1 <- adm1 %>%
    dplyr::mutate(
      region = unname(ctx$GROUP_SHORT_RUN[enso_group]),
      direction = dplyr::case_when(
        significant & lo > 0 ~ "increase",
        significant & hi < 0 ~ "decrease",
        TRUE ~ "none"
      )
    )
  sig <- dplyr::filter(adm1, direction != "none")
  # n_pop_missing BEFORE pop: dplyr::summarise evaluates left-to-right, so the
  # missing count must read the raw per-ADM1 pop column, not the summed scalar.
  sig_pop_region <- sig %>%
    dplyr::group_by(region, enso_group, enso_label, enso_val, direction) %>%
    dplyr::summarise(
      n_adm1 = dplyr::n(),
      n_pop_missing = sum(is.na(pop)),
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    )
  sig_pop_country <- sig %>%
    dplyr::group_by(
      region,
      enso_group,
      country,
      enso_label,
      enso_val,
      direction
    ) %>%
    dplyr::summarise(
      n_adm1 = dplyr::n(),
      n_pop_missing = sum(is.na(pop)),
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    )
  # Africa-wide total: the same significant-effect population summed across every
  # region (no region/country grouping).
  sig_pop_africa <- sig %>%
    dplyr::group_by(enso_label, enso_val, direction) %>%
    dplyr::summarise(
      n_adm1 = dplyr::n(),
      n_pop_missing = sum(is.na(pop)),
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    )

  # --- 95% CrI on the significant-effect population -------------------------
  # The point-estimate sums above use each ADM1's posterior-summary lo/hi to
  # decide significance, then sum populations deterministically. dir_pop_ci()
  # instead sums population per posterior draw (an ADM1 counts on a draw iff its
  # combined log-IRR crosses zero on that draw), giving a genuine median + 95%
  # CrI on "population with an increased/decreased ENSO effect", read off the
  # per-draw array cached in enso_relation (no fit reload). Best-effort: a
  # failure leaves the CrI columns NA and never breaks the rest.
  dir_ci <- tryCatch(
    {
      parts <- lapply(groups, function(g) {
        dir_pop_ci(
          enso_rel$adm1_draws,
          g,
          pop_map = pop_by_adm1,
          enso_conds = conds,
          group_short = ctx$GROUP_SHORT_RUN
        )
      })
      parts <- parts[!vapply(parts, is.null, logical(1L))]
      # Africa-wide CrI: the per-draw region_inc/region_dec matrices are independent
      # group posteriors sharing one draw index d, so summing them element-wise
      # (Reduce(`+`)) yields valid Monte Carlo samples of the Africa-wide population
      # sum. Quantiles are NOT additive, so this must be done from the raw matrices
      # here, never by summing per-group CrIs in the notebook.
      inc_mats <- lapply(parts, `[[`, "region_inc")
      dec_mats <- lapply(parts, `[[`, "region_dec")
      africa <- if (length(inc_mats)) {
        dplyr::bind_rows(
          .pop_ci_rows(Reduce(`+`, inc_mats), "increase", conds),
          .pop_ci_rows(Reduce(`+`, dec_mats), "decrease", conds)
        )
      } else {
        tibble::tibble()
      }
      list(
        region = dplyr::bind_rows(lapply(parts, `[[`, "region_ci")),
        country = dplyr::bind_rows(lapply(parts, `[[`, "country_ci")),
        africa = africa
      )
    },
    error = function(e) {
      warning("[09] dir_pop_ci table failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(dir_ci) && nrow(dir_ci$region)) {
    sig_pop_region <- dplyr::left_join(
      sig_pop_region,
      dir_ci$region,
      by = c("region", "enso_group", "enso_label", "enso_val", "direction")
    )
    sig_pop_country <- dplyr::left_join(
      sig_pop_country,
      dir_ci$country,
      by = c(
        "region",
        "enso_group",
        "country",
        "enso_label",
        "enso_val",
        "direction"
      )
    )
  } else {
    sig_pop_region <- dplyr::mutate(
      sig_pop_region,
      pop_lo = NA_real_,
      pop_med = NA_real_,
      pop_hi = NA_real_
    )
    sig_pop_country <- dplyr::mutate(
      sig_pop_country,
      pop_lo = NA_real_,
      pop_med = NA_real_,
      pop_hi = NA_real_
    )
  }
  if (!is.null(dir_ci) && nrow(dir_ci$africa)) {
    sig_pop_africa <- dplyr::left_join(
      sig_pop_africa,
      dir_ci$africa,
      by = c("enso_label", "enso_val", "direction")
    )
  } else {
    sig_pop_africa <- dplyr::mutate(
      sig_pop_africa,
      pop_lo = NA_real_,
      pop_med = NA_real_,
      pop_hi = NA_real_
    )
  }

  # --- shares of the population at risk -------------------------------------
  # Denominator per scope = observed-sub-unit population of every ADM1 this
  # track's fits carry (the universe the significance test ran over), so
  # "x% of the population" is read off the same basis as the numerator.
  univ <- adm1 %>% dplyr::distinct(adm1_gid, enso_group, region, country, pop)
  den_region <- univ %>%
    dplyr::group_by(region, enso_group) %>%
    dplyr::summarise(pop_universe = sum(pop, na.rm = TRUE), .groups = "drop")
  den_country <- univ %>%
    dplyr::group_by(region, enso_group, country) %>%
    dplyr::summarise(pop_universe = sum(pop, na.rm = TRUE), .groups = "drop")
  sig_pop_region <- add_pop_shares(dplyr::left_join(
    sig_pop_region,
    den_region,
    by = c("region", "enso_group")
  ))
  sig_pop_country <- add_pop_shares(dplyr::left_join(
    sig_pop_country,
    den_country,
    by = c("region", "enso_group", "country")
  ))
  sig_pop_africa <- add_pop_shares(dplyr::mutate(
    sig_pop_africa,
    pop_universe = sum(univ$pop, na.rm = TRUE)
  ))

  # --- duration only: restricted-mean-duration (RMD) ratio -----------------
  # The REAL "outbreaks last XX times longer" estimand (median + 95% CrI),
  # contrasting each composite ONI condition vs neutral (0) over
  # DURATION_RATIO_HORIZON months. Best-effort: any failure (e.g. a fit frame
  # lacking the columns the newdata needs) leaves rmd_ratio NULL and never
  # breaks the rest of the artefact.
  rmd_ratio <- NULL
  if (identical(track, "duration")) {
    rmd_ratio <- tryCatch(
      {
        dur_dat <- readRDS(require_input(ctx$DUR_DATA_PATH, "duration dataset"))
        mmap <- dur_dat %>% dplyr::distinct(mso, mso_band)
        mso_to_band <- stats::setNames(
          as.character(mmap$mso_band),
          as.character(mmap$mso)
        )
        H <- ctx$DURATION_RATIO_HORIZON
        dplyr::bind_rows(lapply(groups, function(g) {
          ff <- track_fit_files(mnf, track, g)
          if (!file.exists(ff$fit) || !file.exists(ff$meta)) {
            return(NULL)
          }
          fit <- readRDS(ff$fit)
          meta <- readRDS(ff$meta)
          vm <- meta[["enso"]]
          res_g <- dplyr::bind_rows(lapply(seq_along(conds), function(i) {
            v <- unname(conds[[i]])
            dr <- dur_rmd_ratio_draws(
              fit,
              vm,
              v,
              H,
              mso_to_band,
              re_formula = NA
            )
            dplyr::mutate(
              duration_ratio_summary(dr),
              enso_label = names(conds)[i],
              enso_val = v
            )
          }))
          rm(fit)
          gc(verbose = FALSE)
          dplyr::mutate(
            res_g,
            enso_group = g,
            region = unname(ctx$GROUP_SHORT_RUN[g])
          )
        }))
      },
      error = function(e) {
        warning("[09] duration RMD ratio table failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  list(
    track = track,
    estimand = if (identical(track, "duration")) {
      "cessation hazard ratio (high = ends sooner = shorter)"
    } else {
      "incidence rate ratio"
    },
    pop_ref_year = ctx$POP_REF_YEAR,
    best_lag = best_lag,
    effects_region = effects_region,
    effects_country = effects_country,
    sig_pop_region = sig_pop_region,
    sig_pop_country = sig_pop_country,
    sig_pop_africa = sig_pop_africa,
    rmd_ratio = rmd_ratio
  )
}

ensure_report_effects <- function(ctx, track) {
  tf <- track_fit_input_files(ctx, track)
  fit_files <- tf$fit
  meta_files <- tf$meta
  lag_loo_files <- lag_loo_input_files(ctx, track)
  extra <- if (identical(track, "duration")) ctx$DUR_DATA_PATH else character(0)
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$GROUPING_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      ctx$GID_POP_PATH,
      extra,
      fit_files,
      meta_files,
      lag_loo_files
    )
  )
  ensure_artefact(
    ctx,
    "report_effects",
    inputs,
    function() derive_report_effects(ctx, track),
    track = track
  )
}

# --- R3: historical smoothing summaries ------------------------------------

#' Per-(rep, month) metric total over a scope's event subset, zero-filled onto
#' the full rep x month grid -- a scope/month with no events must still
#' contribute a zero to the monthly mean, exactly as ribbon_from_events() does
#' for the Figure-2 ribbons. `metric`: "cases" sums the mon-track cases_month,
#' "onsets" counts onset-event rows.
rep_month_totals <- function(
  events_scope,
  rep_ids_all,
  all_months_vec,
  metric
) {
  agg <- if (metric == "cases") {
    events_scope %>%
      dplyr::group_by(rep_id, t_start) %>%
      dplyr::summarise(value = sum(cases_month, na.rm = TRUE), .groups = "drop")
  } else {
    events_scope %>%
      dplyr::group_by(rep_id, t_start) %>%
      dplyr::summarise(value = dplyr::n(), .groups = "drop")
  }
  grid <- expand.grid(
    rep_id = rep_ids_all,
    t_start = all_months_vec,
    stringsAsFactors = FALSE
  )
  grid$t_start <- as.Date(grid$t_start, origin = "1970-01-01")
  out <- merge(grid, agg, by = c("rep_id", "t_start"), all.x = TRUE)
  out$value[is.na(out$value)] <- 0
  tibble::as_tibble(out)
}

derive_report_smoothing <- function(ctx) {
  man <- ctx$MAN
  sm <- load_smoothing_scope(ctx, "report-smoothing") # included groups, frame, grouping, event tables
  groups <- sm$groups
  surv_dat <- sm$surv_dat
  frame <- sm$frame
  enso_grp <- sm$enso_grp
  events_all <- sm$events_all
  events_mon_all <- sm$events_mon_all
  rep_ids_seen <- sort(unique(c(events_all$rep_id, events_mon_all$rep_id)))

  # ENSO phase per calendar month (results.derivation.enso_phase.smoothing):
  # "threshold" = the threshold phase of the monthly ONI (-/+0.5, no
  # persistence; the cache derive_oni_envelope() reads), the legacy table;
  # "episode" = the NOAA CPC episodes of the `oni` artefact, the runs Figures 1
  # and 3 shade, so the table and the shading agree.
  phase_rule <- smoothing_phase_rule(ctx)
  oni <- smoothing_month_phases(ctx)

  scope_gids <- smoothing_scope_gids(
    frame,
    events_all,
    enso_grp,
    groups,
    ctx$GROUP_SHORT_RUN
  )
  ev_by_metric <- list(cases = events_mon_all, onsets = events_all)

  monthly_overall <- list()
  monthly_by_year <- list()
  monthly_by_phase <- list()
  for (scope in names(scope_gids)) {
    gids <- scope_gids[[scope]]
    for (metric in c("cases", "onsets")) {
      ev_s <- dplyr::filter(ev_by_metric[[metric]], gid %in% gids)
      months_vec <- sort(unique(ev_s$t_start))
      if (!length(months_vec)) {
        next
      }
      # Per-rep per-month totals (zero-filled), tagged with year + ONI phase.
      rmt <- rep_month_totals(ev_s, rep_ids_seen, months_vec, metric) %>%
        dplyr::mutate(
          year = lubridate::year(t_start),
          mon = lubridate::month(t_start)
        ) %>%
        dplyr::left_join(oni, by = c("year", "mon" = "month")) %>%
        dplyr::mutate(phase = tidyr::replace_na(phase, "neutral"))

      # Mean monthly value = per-rep mean over months, then median + 95% across
      # replicates (the simulation-output reporting convention).
      ov <- rmt %>%
        dplyr::group_by(rep_id) %>%
        dplyr::summarise(mean_monthly = mean(value), .groups = "drop")
      monthly_overall[[paste(scope, metric)]] <-
        dplyr::mutate(.q_med95(ov$mean_monthly), scope = scope, metric = metric)

      by_yr <- rmt %>%
        dplyr::group_by(rep_id, year) %>%
        dplyr::summarise(mean_monthly = mean(value), .groups = "drop") %>%
        dplyr::group_by(year) %>%
        dplyr::group_modify(~ .q_med95(.x$mean_monthly)) %>%
        dplyr::ungroup()
      monthly_by_year[[paste(scope, metric)]] <-
        dplyr::mutate(by_yr, scope = scope, metric = metric)

      by_ph <- rmt %>%
        dplyr::group_by(rep_id, phase) %>%
        dplyr::summarise(mean_monthly = mean(value), .groups = "drop") %>%
        dplyr::group_by(phase) %>%
        dplyr::group_modify(~ .q_med95(.x$mean_monthly)) %>%
        dplyr::ungroup()
      monthly_by_phase[[paste(scope, metric)]] <-
        dplyr::mutate(by_ph, scope = scope, metric = metric)
    }
  }

  list(
    monthly_overall = dplyr::bind_rows(monthly_overall),
    monthly_by_year = dplyr::bind_rows(monthly_by_year),
    monthly_by_phase = dplyr::bind_rows(monthly_by_phase),
    # How monthly_by_phase's phases were made (additive field).
    phase_rule = phase_rule,
    groups = groups
  )
}

#' The phase definition report_smoothing uses (results.derivation.enso_phase.smoothing).
smoothing_phase_rule <- function(ctx) {
  if (identical(ctx$ENSO_PHASE_SMOOTHING, "episode")) {
    enso_phase_definition("episode", "monthly_oni_cfsv2_ext")
  } else {
    enso_phase_definition("threshold", "monthly_oni")
  }
}

#' Phase per calendar month for report_smoothing: tibble(year, month, phase).
#' "threshold" reads the ONI cache (the legacy table, bit-identical); "episode"
#' expands the `oni` artefact's episode runs (one row per run, TL / TR = first /
#' last month -/+ 15 days) to their months.
smoothing_month_phases <- function(ctx) {
  if (!identical(ctx$ENSO_PHASE_SMOOTHING, "episode")) {
    return(
      fetch_enso_index(
        source = "oni",
        cache_dir = ctx$ENSO_CACHE_DIR,
        refresh = FALSE
      ) %>%
        dplyr::transmute(
          year = as.integer(year),
          month = as.integer(month),
          phase = tidyr::replace_na(
            classify_enso_phase(enso_index, "threshold"),
            "neutral"
          )
        )
    )
  }
  ph <- ensure_oni(ctx, ctx$SETS)$phases
  smoothing_phases_from_runs(ph)
}

#' Months of each episode run of the `oni` artefact's `phases` table.
smoothing_phases_from_runs <- function(ph) {
  rows <- lapply(seq_len(nrow(ph)), function(i) {
    d <- seq(as.Date(ph$TL[i]) + 15, as.Date(ph$TR[i]) - 15, by = "1 month")
    tibble::tibble(
      year = as.integer(format(d, "%Y")),
      month = as.integer(format(d, "%m")),
      phase = ph$phase[i]
    )
  })
  dplyr::bind_rows(rows)
}

#' The `oni` artefact file (report_smoothing's input under the episode rule).
.oni_artefact_file <- function(ctx) {
  make_results_run_file(
    ctx$opt,
    "results_oni",
    sets = ctx$SETS,
    schema = RESULTS_SCHEMA[["oni"]],
    config_path = ctx$CONFIG_PATH
  )
}

ensure_report_smoothing <- function(ctx) {
  # The ONI cache is an input under either rule (a refreshed ONI re-derives the
  # table); the `oni` artefact under episode.
  extra <- c(
    .oni_cache_file(ctx),
    if (identical(ctx$ENSO_PHASE_SMOOTHING, "episode")) .oni_artefact_file(ctx)
  )
  inputs <- results_inputs(
    files = c(
      ctx$PANEL_PATH,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$MON_DATA_PATH,
      ctx$GROUPING_PATH,
      ctx$MANIFEST_PATH,
      extra
    ),
    globs = smoothing_parts_glob(ctx)
  )
  ensure_artefact(ctx, "report_smoothing", inputs, function() {
    derive_report_smoothing(ctx)
  })
}

# --- R4: forecast / counterfactual horizon numbers -------------------------

derive_report_forecast <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  win <- results_arm_window(ctx, groups, sets)
  attrib_start <- win$attrib_start
  africa_min_horizon <- as.Date(man$africa_min_horizon)
  # Africa (results.scope.africa): every group up to its end in the Africa
  # window (legacy: africa_min_horizon for all, the former shared window).
  aw <- .africa_win(ctx, groups)
  afr_end <- aw$group_ends
  gid_pop <- load_gid_population(ctx$GID_POP_PATH)

  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  frame <- build_location_frame(
    survival_data = surv_dat,
    hmm_clusters_path = ctx$HMM_PATH,
    outbreak_metadata_path = NULL,
    size_data_path = ctx$SIZE_DATA_PATH
  )
  country_region <- frame %>% dplyr::distinct(country, AFRO_region)

  # --- load forecast + counterfactual monthly partitions (as Figures 2-3) ---
  # Of the forecast manifest (horizon totals: the displayed forecast) and, for
  # the paired ratios and population tables below, of the attribution manifest
  # when there is an attribution factual (results_attrib_ctx()).
  arm_frames <- function(cx, w) {
    runs <- load_group_arm_runs(cx, groups, sets, w$attrib_start, win = w)
    bf <- .bound_run_frames(cx, runs)
    country_month_both <- .lookup_join(
      bf$country_month,
      country_region,
      "country"
    )
    # adm1_month: each group's own firm window (pop_grid / pop_shift, the mixed
    # Africa), cut there when the Africa window loaded continuation months.
    list(
      runs = runs,
      country_month = country_month_both,
      adm1_month = if (isTRUE(aw$extended)) {
        dplyr::filter(bf$adm1_month, month <= group_horizon_end)
      } else {
        bf$adm1_month
      },
      cm_win = country_month_both %>%
        dplyr::filter(month >= w$attrib_start, month <= group_horizon_end),
      afr_win = country_month_both %>%
        dplyr::filter(
          month >= w$attrib_start,
          month <= unname(afr_end[enso_group])
        )
    )
  }
  fr <- arm_frames(ctx, win)
  runs <- fr$runs
  group_horizon <- runs$group_horizon
  cfgs <- runs$cfgs
  adm1_month_both <- fr$adm1_month
  cm_win <- fr$cm_win
  afr_win <- fr$afr_win
  fr <- NULL

  # ---------------------- horizon totals (median + 95%) --------------------
  # Population of the simulated gids per scope (canonical table at POP_REF_YEAR)
  # and the months each scope's total spans: each group's own [attrib_start,
  # horizon_end], Africa's shared [attrib_start, africa_min_horizon]. Totals are
  # over the window; inc_ann_100k annualises them so windows of different
  # length compare. The annualising month count is the number of months the
  # runs actually SIMULATED inside the window (monthly_summary is a full
  # replicate x gid x month grid, so every simulated month is present), not the
  # calendar length: a run that stops short of its configured horizon must not
  # be divided by months it never simulated. A shortfall is warned about and
  # both counts are kept in the `n_months` table.
  # rate_pop_year = "window": the incidence denominators are the
  # month-weighted mean population over the months the totals read (cm_win).
  scope_pop <- if (identical(ctx[["RATE_POP_YEAR"]], "window")) {
    build_sim_scope_pop(
      runs$gids,
      gid_pop,
      ctx$POP_REF_YEAR,
      ctx$GROUP_SHORT_RUN,
      pop_year_weights = rate_pop_weights(cm_win$month),
      wpp = load_wpp_table(.pop_artefact_or_null(ctx$WPP_PATH))
    )
  } else {
    build_sim_scope_pop(
      runs$gids,
      gid_pop,
      ctx$POP_REF_YEAR,
      ctx$GROUP_SHORT_RUN
    )
  }
  sim_months <- function(df) dplyr::n_distinct(df$month)
  n_window_grp <- vapply(
    groups,
    function(g) window_n_months(attrib_start, group_horizon[[g]]),
    integer(1)
  )
  n_months_grp <- vapply(
    groups,
    function(g) {
      # The forecast arm's analog_set is NA, which split() would drop as a group.
      cm_g <- cm_win[cm_win$enso_group == g, ]
      by_arm <- vapply(
        split(
          cm_g,
          list(cm_g$run, dplyr::coalesce(as.character(cm_g$analog_set), "-")),
          drop = TRUE
        ),
        sim_months,
        integer(1)
      )
      if (length(unique(by_arm)) > 1L) {
        warning(sprintf(
          "[09] report_forecast: group %s arms simulate different month counts (%s); using the minimum",
          g,
          paste(sprintf("%s=%d", names(by_arm), by_arm), collapse = ", ")
        ))
      }
      min(by_arm)
    },
    integer(1)
  )
  n_window_afr <- window_n_months(attrib_start, aw$end)
  n_months_afr <- sim_months(afr_win)
  for (g in groups[n_months_grp != n_window_grp]) {
    warning(sprintf(
      paste(
        "[09] report_forecast: group %s simulates %d of the %d months in [%s, %s];",
        "totals cover the simulated months and incidence is annualised over them"
      ),
      g,
      n_months_grp[[g]],
      n_window_grp[[g]],
      attrib_start,
      group_horizon[[g]]
    ))
  }
  pop_region <- dplyr::filter(scope_pop, level == "region") %>%
    dplyr::select(enso_group, pop)
  pop_country <- dplyr::filter(scope_pop, level == "country") %>%
    dplyr::transmute(enso_group, country = scope, pop)
  pop_africa <- scope_pop$pop[scope_pop$level == "continent"]
  cm_win_n <- cm_win %>%
    dplyr::mutate(n_months = unname(n_months_grp[enso_group]))

  QTY <- c(
    "tot_cases",
    "tot_onsets",
    "mean_cases",
    "mean_onsets",
    "inc_ann_100k"
  )
  summ_totals <- function(per_rep, key_cols) {
    tidyr::pivot_longer(
      per_rep,
      dplyr::all_of(QTY),
      names_to = "quantity",
      values_to = "v"
    ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(
        key_cols,
        "run",
        "analog_set",
        "quantity"
      )))) %>%
      dplyr::summarise(
        med = median(v, na.rm = TRUE),
        lo = quantile(v, .025, na.rm = TRUE, names = FALSE),
        hi = quantile(v, .975, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      )
  }

  totals_region <- summ_totals(
    per_rep_totals(
      dplyr::left_join(cm_win_n, pop_region, by = "enso_group"),
      "enso_group"
    ),
    "enso_group"
  ) %>%
    dplyr::mutate(
      level = "region",
      scope = unname(ctx$GROUP_SHORT_RUN[enso_group])
    )
  totals_country <- summ_totals(
    per_rep_totals(
      dplyr::left_join(cm_win_n, pop_country, by = c("enso_group", "country")),
      c("enso_group", "country")
    ),
    c("enso_group", "country")
  ) %>%
    dplyr::mutate(level = "country", scope = country)
  # Africa: paired replicate ids common to every group, summed across groups,
  # within the Africa window (aw: every group to its end in it).
  scenarios <- dplyr::bind_rows(
    tibble::tibble(run = "forecast", analog_set = NA_character_),
    tibble::tibble(run = "counterfactual", analog_set = sets)
  )
  totals_africa <- dplyr::bind_rows(lapply(
    seq_len(nrow(scenarios)),
    function(i) {
      rn <- scenarios$run[i]
      as_ <- scenarios$analog_set[i]
      sub <- afr_win %>%
        dplyr::filter(run == rn, if (is.na(as_)) TRUE else analog_set == as_)
      common <- common_replicate_ids(split(sub, sub$enso_group))
      pr <- sub %>%
        dplyr::filter(replicate_id %in% common) %>%
        dplyr::group_by(replicate_id, month) %>%
        dplyr::summarise(
          n_cases = sum(n_cases),
          n_outbreaks = sum(n_outbreaks),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          run = rn,
          analog_set = as_,
          n_months = n_months_afr,
          pop = pop_africa
        )
      summ_totals(per_rep_totals(pr, character(0)), character(0))
    }
  )) %>%
    dplyr::mutate(level = "continent", scope = "Africa")
  horizon_totals <- dplyr::bind_rows(
    totals_region,
    totals_country,
    totals_africa
  ) %>%
    dplyr::mutate(region = unname(ctx$GROUP_SHORT_RUN[enso_group])) %>%
    dplyr::select(
      level,
      scope,
      enso_group,
      region,
      run,
      analog_set,
      quantity,
      med,
      lo,
      hi
    )
  n_months <- dplyr::bind_rows(
    tibble::tibble(
      level = "continent",
      scope = "Africa",
      enso_group = NA_character_,
      n_months = n_months_afr,
      n_months_window = n_window_afr
    ),
    tibble::tibble(
      level = "region",
      scope = unname(ctx$GROUP_SHORT_RUN[groups]),
      enso_group = groups,
      n_months = unname(n_months_grp),
      n_months_window = unname(n_window_grp)
    )
  )

  # --- the attribution pair (results.inputs.attrib_manifest) ----------------
  # The ratios and population tables below are paired factual-vs-counterfactual
  # statistics: they read the attribution factual when there is one (display
  # frames released first, one load in memory at a time).
  win_a <- win
  if (identical(ctx$ATTRIB_FACTUAL, "attrib")) {
    adm1_month_both <- NULL
    cm_win <- NULL
    afr_win <- NULL
    runs <- NULL
    invisible(gc(verbose = FALSE))
    cx_a <- results_attrib_ctx(ctx)
    win_a <- results_arm_window(cx_a, groups, sets)
    if (win_a$attrib_start != attrib_start) {
      stop(
        sprintf(
          "[09] report_forecast: the attribution window starts %s, the display window %s.",
          win_a$attrib_start,
          attrib_start
        ),
        call. = FALSE
      )
    }
    fa <- arm_frames(cx_a, win_a)
    adm1_month_both <- fa$adm1_month
    cm_win <- fa$cm_win
    afr_win <- fa$afr_win
    fa <- NULL
  }

  # ---------------------- paired forecast/cf ratios ------------------------
  ratio_one <- function(df_f, df_c, value_col) {
    m <- dplyr::inner_join(
      df_f,
      df_c,
      by = "replicate_id",
      suffix = c("_f", "_c")
    )
    r <- m[[paste0(value_col, "_f")]] / m[[paste0(value_col, "_c")]]
    r[!is.finite(r)] <- NA_real_
    r
  }
  summ_ratio <- function(r) {
    tibble::tibble(
      med = median(r, na.rm = TRUE),
      lo = quantile(r, .025, na.rm = TRUE, names = FALSE),
      hi = quantile(r, .975, na.rm = TRUE, names = FALSE),
      prob_increase = mean(r > 1, na.rm = TRUE)
    )
  }

  # Cumulative-total ratio and peak (per-rep max monthly) ratio share the same
  # region/country/Africa pairing skeleton, so build it once and swap the
  # per-(scope, arm, rep) reduction.
  CO <- c(cases = "n_cases", onsets = "n_outbreaks")
  cum_group <- .sum_by_keys(
    cm_win,
    c("enso_group", "run", "analog_set", "replicate_id"),
    CO
  )
  cum_country <- .sum_by_keys(
    cm_win,
    c("enso_group", "country", "run", "analog_set", "replicate_id"),
    CO
  )
  # Region peak = max over months of the REGION's monthly total (cm_win rows are
  # country-months since load_sim_run_monthly() keys countries by gid ISO3;
  # before that fix each region was a single pseudo-country, so this matches
  # the pre-fix values).
  peak_group <- .sum_by_keys(
    .sum_by_keys(
      cm_win,
      c("enso_group", "run", "analog_set", "replicate_id", "month"),
      c(n_cases = "n_cases", n_outbreaks = "n_outbreaks")
    ),
    c("enso_group", "run", "analog_set", "replicate_id"),
    CO,
    op = "max"
  )
  peak_country <- .sum_by_keys(
    cm_win,
    c("enso_group", "country", "run", "analog_set", "replicate_id"),
    CO,
    op = "max"
  )

  # Africa arm reduction: sum across groups per (rep, month) over common ids,
  # then per-rep either sum (cumulative) or max (peak) over the window.
  afr_arm <- function(rn, as_, value_col, reduce) {
    sub <- afr_win %>%
      dplyr::filter(run == rn, if (is.na(as_)) TRUE else analog_set == as_)
    common <- common_replicate_ids(split(sub, sub$enso_group))
    per_month <- sub %>%
      dplyr::filter(replicate_id %in% common) %>%
      dplyr::group_by(replicate_id, month) %>%
      dplyr::summarise(v = sum(.data[[value_col]]), .groups = "drop")
    per_month %>%
      dplyr::group_by(replicate_id) %>%
      dplyr::summarise(val = reduce(v), .groups = "drop")
  }
  METRIC_COLS <- c(cases = "cases", onsets = "onsets")
  AFRICA_MON_COL <- c(cases = "n_cases", onsets = "n_outbreaks")
  build_ratio <- function(grp_df, cty_df, reduce) {
    dplyr::bind_rows(lapply(names(METRIC_COLS), function(metric) {
      vc <- METRIC_COLS[[metric]]
      dplyr::bind_rows(lapply(sets, function(aset) {
        per_group <- dplyr::bind_rows(lapply(groups, function(g) {
          df_f <- dplyr::filter(grp_df, enso_group == g, run == "forecast")
          df_c <- dplyr::filter(
            grp_df,
            enso_group == g,
            run == "counterfactual",
            analog_set == aset
          )
          dplyr::mutate(
            summ_ratio(ratio_one(df_f, df_c, vc)),
            level = "region",
            scope = unname(ctx$GROUP_SHORT_RUN[g])
          )
        }))
        per_country <- dplyr::bind_rows(lapply(groups, function(g) {
          ctys <- unique(cty_df$country[cty_df$enso_group == g])
          dplyr::bind_rows(lapply(ctys, function(ct) {
            df_f <- dplyr::filter(
              cty_df,
              country == ct,
              enso_group == g,
              run == "forecast"
            )
            df_c <- dplyr::filter(
              cty_df,
              country == ct,
              enso_group == g,
              run == "counterfactual",
              analog_set == aset
            )
            dplyr::mutate(
              summ_ratio(ratio_one(df_f, df_c, vc)),
              level = "country",
              scope = ct
            )
          }))
        }))
        af_f <- afr_arm(
          "forecast",
          NA_character_,
          AFRICA_MON_COL[[metric]],
          reduce
        )
        af_c <- afr_arm(
          "counterfactual",
          aset,
          AFRICA_MON_COL[[metric]],
          reduce
        )
        afr <- dplyr::mutate(
          summ_ratio(ratio_one(af_f, af_c, "val")),
          level = "continent",
          scope = "Africa"
        )
        dplyr::bind_rows(per_group, per_country, afr) %>%
          dplyr::mutate(analog_set = aset, metric = metric)
      }))
    }))
  }
  cum_ratio <- build_ratio(cum_group, cum_country, sum)
  peak_ratio <- build_ratio(peak_group, peak_country, max)

  # ---------------------- threshold-population grid ------------------------
  # Population in admin units with >= p posterior chance of a >= m fractional
  # increase/decrease over the horizon, for every (cutoff x magnitude) pair.
  # Observed sub-ADM1 population map (shared with R1/R2 via build_adm1_obs_pop).
  # recover_dropped = FALSE keeps this forecast universe byte-identical to the
  # pre-helper behaviour; R4's own enso_group/region join follows.
  adm1_pop <- build_forecast_adm1_pop(
    ctx,
    frame,
    surv_dat,
    gid_pop,
    adm1_month_both
  )

  ADM1_METRIC <- c(cases = "cases", onsets = "outbreaks")
  grid_rows <- list()
  for (metric in names(ADM1_METRIC)) {
    for (aset in sets) {
      base_cum <- adm1_prob_cum(
        dplyr::filter(adm1_month_both, run == "forecast" | analog_set == aset),
        ADM1_METRIC[[metric]]
      ) # once per (metric, set), not per magnitude
      for (m in ctx$POP_MAGNITUDES) {
        st <- adm1_prob_stats_cum(base_cum, m, keep_means = FALSE) %>%
          dplyr::left_join(adm1_pop, by = "adm1_unit")
        for (p in ctx$POP_PROB_CUTOFFS) {
          for (dir in c("increase", "decrease")) {
            flagged <- if (dir == "increase") {
              dplyr::filter(st, !is.na(prob_larger), prob_larger >= p)
            } else {
              dplyr::filter(st, !is.na(prob_smaller), prob_smaller >= p)
            }
            mk <- function(level, scope, sub) {
              tibble::tibble(
                metric = metric,
                analog_set = aset,
                magnitude = m,
                cutoff = p,
                direction = dir,
                level = level,
                scope = scope,
                n_adm1 = nrow(sub),
                pop = sum(sub$pop, na.rm = TRUE),
                n_pop_missing = sum(sub$n_pop_missing, na.rm = TRUE)
              )
            }
            rows <- list(mk("continent", "Africa", flagged))
            for (rg in unique(stats::na.omit(flagged$region))) {
              rows <- c(
                rows,
                list(mk("region", rg, dplyr::filter(flagged, region == rg)))
              )
            }
            for (ct in unique(stats::na.omit(flagged$country))) {
              rows <- c(
                rows,
                list(mk("country", ct, dplyr::filter(flagged, country == ct)))
              )
            }
            grid_rows <- c(grid_rows, rows)
          }
        }
      }
    }
  }
  pop_grid <- dplyr::bind_rows(grid_rows) %>%
    dplyr::mutate(africa_def = "mixed") # each group over its firm window

  # ---------------- per-replicate population under shift (95% CrI) ----------
  # pop_grid above is a probability-threshold COUNT (population in units with
  # >= p chance of a >= m shift) and is deterministic given the replicate
  # distribution. pop_shift is its simulation-output CrI companion: for each
  # replicate we pair forecast vs counterfactual cumulative totals per ADM1
  # (by replicate_id, exactly as ratio_one() pairs the ratio draws), sum the
  # population of every ADM1 whose forecast exceeds (1+m)x [< (1-m)x] its paired
  # counterfactual, then take median + 95% across replicates. No cutoff axis:
  # the spread over replicates IS the uncertainty.
  pop_shift_rows <- list()
  cum_all <- pair_adm1_cum_table(adm1_month_both) # once, every metric and set
  for (metric in names(ADM1_METRIC)) {
    for (aset in sets) {
      vc <- ADM1_METRIC[[metric]]
      paired <- pair_adm1_cum(adm1_month_both, vc, aset, cum = cum_all) %>%
        dplyr::left_join(adm1_pop, by = "adm1_unit")
      for (m in ctx$POP_MAGNITUDES) {
        for (dir in c("increase", "decrease")) {
          shifted <- if (dir == "increase") {
            dplyr::mutate(paired, hit = f > (1 + m) * cf)
          } else {
            dplyr::mutate(paired, hit = f < (1 - m) * cf)
          }
          # Per replicate: sum(pop[hit], na.rm = TRUE), by tapply() (the same sum()
          # over the same elements in row order as the former group_by(replicate_id);
          # ~10^4 replicates x ~40 scopes per cell).
          per_rep_pop <- function(sub) {
            unname(tapply(
              ifelse(sub$hit, sub$pop, NA_real_),
              sub$replicate_id,
              sum,
              na.rm = TRUE
            ))
          }
          qrow <- function(level, scope, sub) {
            pv <- per_rep_pop(sub)
            tibble::tibble(
              metric = metric,
              analog_set = aset,
              magnitude = m,
              direction = dir,
              level = level,
              scope = scope,
              pop_lo = quantile(pv, .025, names = FALSE, na.rm = TRUE),
              pop_med = median(pv, na.rm = TRUE),
              pop_hi = quantile(pv, .975, names = FALSE, na.rm = TRUE)
            )
          }
          # continent: restrict to replicate ids common to every group so the
          # per-replicate sum spans the same Africa in each draw (mirrors afr_arm).
          common <- Reduce(
            intersect,
            split(shifted$replicate_id, shifted$enso_group)
          )
          cont_sub <- dplyr::filter(shifted, replicate_id %in% common)
          rows <- list(qrow("continent", "Africa", cont_sub))
          for (rg in unique(stats::na.omit(shifted$region))) {
            rows <- c(
              rows,
              list(qrow("region", rg, dplyr::filter(shifted, region == rg)))
            )
          }
          for (ct in unique(stats::na.omit(shifted$country))) {
            rows <- c(
              rows,
              list(qrow("country", ct, dplyr::filter(shifted, country == ct)))
            )
          }
          pop_shift_rows <- c(pop_shift_rows, rows)
        }
      }
    }
  }
  pop_shift <- dplyr::bind_rows(pop_shift_rows) %>%
    dplyr::mutate(africa_def = "mixed")

  list(
    horizon_totals = horizon_totals,
    cum_ratio = cum_ratio,
    peak_ratio = peak_ratio,
    pop_grid = pop_grid,
    pop_shift = pop_shift,
    scope_pop = scope_pop,
    n_months = n_months,
    attrib_start = attrib_start,
    africa_min_horizon = africa_min_horizon,
    pop_ref_year = ctx$POP_REF_YEAR,
    analog_sets = sets,
    groups = groups,
    factual_start = win$factual_start,
    arm_starts = win$arm_starts,
    # horizon_totals: the displayed forecast; cum_ratio / peak_ratio /
    # pop_grid / pop_shift: the factual of attrib_factual ("attrib" = the
    # attribution factual, launching at attrib_arm_starts$factual).
    attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
    attrib_arm_starts = win_a$arm_starts,
    # The Africa of horizon_totals / cum_ratio / peak_ratio / n_months
    # (results.scope.africa); pop_grid / pop_shift stay the mixed sum
    # (africa_def "mixed": each group over its own firm window).
    africa_window = africa_window_payload(aw, attrib_start)
  ) %>%
    c(.sim_pop_year_fields(ctx, list(ctx$MAN_DISPLAY, ctx$MAN_ATTRIB))) # report_forecast v5
}

ensure_report_forecast <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- .display_attrib_configs(ctx, groups, sets)
  inputs <- results_inputs(
    files = c(
      unique(c(
        ctx$MANIFEST_PATH,
        ctx$ATTRIB_MANIFEST_FILE %||% ctx$MANIFEST_PATH
      )),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = if (
      results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched
    ) {
      smoothing_parts_glob(ctx)
    } else {
      list()
    }
  )
  ensure_artefact(
    ctx,
    "report_forecast",
    inputs,
    function() derive_report_forecast(ctx, sets),
    sets = sets
  )
}

# --- R5: population in ENSO-sensitive areas ---------------------------------
#
# Shared machinery for every per-draw ENSO population table: the ADM1 draw
# universe (enso_draw_universe), the Africa / region / country scopes
# (enso_scopes), masks (U x C x D logical arrays, one per track plus
# cross-track combinations) reduced to population per draw (mask_pop_rows) or
# to a probability-threshold count (prob_grid_rows). report_enso_exposure is
# thin over these (as is the unregistered derive_report_sensitive_pop()).

#' The ADM1 draw universe of the FIG1 tracks: validates that every track's
#' enso_relation artefact carries per-draw arrays at this config's conditions
#' with one shared draw count, joins the observed-sub-unit population, and
#' returns the per-track log-IRR arrays aligned onto the union of ADM1s
#' (absent ADM1 -> NA; mask it with enso_track_masks()).
#'
#' Draws are paired by index across tracks and across groups: the tracks are
#' independent fits, each subsampled with the same seed and draw count, so draw
#' d of a cross-track statistic combines draw d of each fit -- the same pairing
#' the Africa rollup of dir_pop_ci() makes across groups.
#'
#' @return list(tracks, conds, n_draws, universe = tibble(adm1_gid, country,
#'   enso_group, region, pop, pop_missing), U = adm1_gid vector, pop0 = pop
#'   with NA -> 0, logrr = named list of [U x C x D] arrays, rel = the
#'   enso_relation artefacts)
enso_draw_universe <- function(ctx, tag) {
  tracks <- ctx$FIG1_TRACKS
  conds <- ctx$ENSO_CONDS
  rel <- stats::setNames(
    lapply(tracks, function(t) ensure_enso_relation(ctx, t)),
    tracks
  )
  for (t in tracks) {
    if (is.null(rel[[t]]$adm1_draws)) {
      stop(
        "[09] ",
        tag,
        ": enso_relation (",
        t,
        ") carries no adm1_draws.",
        call. = FALSE
      )
    }
    if (!identical(dimnames(rel[[t]]$adm1_draws$logrr)[[2L]], names(conds))) {
      stop(
        "[09] ",
        tag,
        ": enso_relation (",
        t,
        ") was derived at other ONI conditions.",
        call. = FALSE
      )
    }
  }
  n_d <- vapply(rel, function(r) dim(r$adm1_draws$logrr)[3L], numeric(1))
  if (length(unique(n_d)) != 1L) {
    stop(
      "[09] ",
      tag,
      ": tracks disagree on the draw count (",
      paste(sprintf("%s=%d", tracks, n_d), collapse = ", "),
      "); draws cannot be paired.",
      call. = FALSE
    )
  }

  gid_pop <- load_gid_population(ctx$GID_POP_PATH)
  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  pop_map <- obs_adm1_pop_map(ctx, surv_dat, gid_pop, tag)

  # ADM1 universe: union over tracks, metadata from the first track carrying it.
  universe <- dplyr::bind_rows(lapply(tracks, function(t) {
    dplyr::select(rel[[t]]$adm1_draws$adm1, adm1_gid, country, enso_group)
  })) %>%
    dplyr::distinct(adm1_gid, .keep_all = TRUE) %>%
    dplyr::mutate(
      region = unname(ctx$GROUP_SHORT_RUN[enso_group]),
      pop = unname(pop_map[adm1_gid]),
      pop_missing = is.na(pop)
    )
  U <- universe$adm1_gid
  logrr <- lapply(rel, function(r) {
    align_adm1_array(
      r$adm1_draws$logrr,
      r$adm1_draws$adm1$adm1_gid,
      U,
      fill = NA_real_
    )
  })
  list(
    tracks = tracks,
    conds = conds,
    n_draws = unname(n_d[[1L]]),
    universe = universe,
    U = U,
    pop0 = ifelse(is.na(universe$pop), 0, universe$pop),
    logrr = logrr,
    rel = rel
  )
}

#' Per-track logical masks over the draw universe: `test(logrr, track)` on each
#' aligned array, with the NA of an absent ADM1 counted as FALSE.
enso_track_masks <- function(dr, test) {
  lapply(stats::setNames(dr$tracks, dr$tracks), function(t) {
    m <- test(dr$logrr[[t]], t)
    m[is.na(m)] <- FALSE
    m
  })
}

#' Cross-track combinations of per-track masks: `any` (union) and
#' `concordant2` (at least two tracks agree on the draw).
enso_combine_masks <- function(masks) {
  masks$any <- Reduce(`|`, masks)
  masks$concordant2 <- Reduce(
    `+`,
    lapply(masks[setdiff(names(masks), "any")], function(x) x * 1L)
  ) >=
    2L
  masks
}

#' Africa / each ENSO region / each country of a draw universe, with the ADM1
#' row indices of each scope.
enso_scopes <- function(universe) {
  # group -> short name, as enso_draw_universe() set it from the run's registry
  short <- dplyr::distinct(universe, enso_group, region)
  short <- stats::setNames(short$region, short$enso_group)
  sc <- dplyr::bind_rows(
    tibble::tibble(
      level = "continent",
      scope = "Africa",
      enso_group = NA_character_,
      country = NA_character_
    ),
    universe %>%
      dplyr::distinct(enso_group) %>%
      dplyr::transmute(
        level = "region",
        scope = unname(short[enso_group]),
        enso_group,
        country = NA_character_
      ),
    universe %>%
      dplyr::distinct(enso_group, country) %>%
      dplyr::transmute(level = "country", scope = country, enso_group, country)
  )
  sc$region <- unname(short[sc$enso_group])
  sc$idx <- lapply(seq_len(nrow(sc)), function(k) {
    switch(
      sc$level[k],
      continent = seq_len(nrow(universe)),
      region = which(universe$enso_group == sc$enso_group[k]),
      country = which(
        universe$enso_group == sc$enso_group[k] &
          universe$country == sc$country[k]
      )
    )
  })
  sc
}

# Scope columns + universe counts shared by every per-draw population row.
.scope_cols <- function(dr, sc, k) {
  idx <- sc$idx[[k]]
  tibble::tibble(
    level = sc$level[k],
    scope = sc$scope[k],
    enso_group = sc$enso_group[k],
    region = sc$region[k],
    country = sc$country[k],
    n_adm1_universe = length(idx),
    n_pop_missing = sum(dr$universe$pop_missing[idx]),
    pop_universe = sum(dr$pop0[idx])
  )
}

#' Population exposed per draw: for every mask x scope, the population of the
#' ADM1s whose mask is TRUE on each draw, summarised across draws
#' (.pop_summary_rows: mean, q025, median, q975) plus the mean ADM1 count, with
#' shares of the scope population. `point`, when given, is a parallel list of
#' U x C logical point-estimate matrices (adds n_adm1_point / pop_point).
#' `cond_cols` labels the C axis, one row per column of the masks: the ONI
#' conditions by default; report_attrib_exposure passes tibble(analog_set) and
#' bootstrap resamples as the D axis.
mask_pop_rows <- function(
  masks,
  dr,
  point = NULL,
  cond_cols = enso_cond_cols(dr$conds)
) {
  sc <- enso_scopes(dr$universe)
  rows <- list()
  for (t in names(masks)) {
    for (k in seq_len(nrow(sc))) {
      idx <- sc$idx[[k]]
      mk <- masks[[t]][idx, , , drop = FALSE]
      r <- dplyr::bind_cols(
        .scope_cols(dr, sc, k),
        cond_cols,
        .pop_summary_cols(pop_mass_by_draw(mk, dr$pop0[idx]))
      ) %>%
        dplyr::mutate(track = t, n_adm1_mean = rowMeans(colSums(mk, dims = 1L)))
      if (!is.null(point)) {
        pt <- point[[t]][idx, , drop = FALSE]
        r <- dplyr::mutate(
          r,
          n_adm1_point = colSums(pt),
          pop_point = colSums(pt * dr$pop0[idx])
        )
      }
      rows[[length(rows) + 1L]] <- r
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::relocate(
      track,
      level,
      scope,
      enso_group,
      region,
      country,
      dplyr::all_of(names(cond_cols))
    ) %>%
    add_pop_shares()
}

#' One scope's point-population row from a U x C logical matrix `hit`: the ADM1
#' count and population of the scope's ADM1s with hit TRUE, per condition column
#' (labelled by `cond_cols`), after the scope columns; `...` adds constant
#' columns between the condition columns and the counts (e.g. cutoff = p).
.point_scope_row <- function(hit, dr, sc, k, cond_cols, ...) {
  idx <- sc$idx[[k]]
  h <- hit[idx, , drop = FALSE]
  extra <- list(...)
  dplyr::bind_cols(
    .scope_cols(dr, sc, k),
    cond_cols,
    if (length(extra)) tibble::as_tibble(extra),
    tibble::tibble(
      n_adm1 = unname(colSums(h)),
      pop = unname(colSums(h * dr$pop0[idx]))
    )
  )
}

#' Probability-threshold population: for every mask x scope x cutoff p, the
#' population (and ADM1 count) of the ADM1s whose posterior probability of the
#' masked event -- the mean of the mask over draws -- is at least p. No
#' interval: the probability already carries the posterior uncertainty (the
#' ENSO analogue of the forecast threshold-population grid).
#'
#' `credible_union` (names of masks, e.g. the per-track ones) adds the rows
#' track "any_credible": the TRACK-WISE union, an ADM1 counts when SOME single
#' track has P >= p on its own. Unlike the mask `any` (the certainty of the
#' union), it does not inflate the certainty of several weak tracks: with
#' independent tracks P(any) = 1 - prod(1 - p_t), so three tracks at 0.65 each
#' pass 0.95.
prob_grid_rows <- function(masks, dr, cutoffs, credible_union = NULL) {
  sc <- enso_scopes(dr$universe)
  rows <- list()
  prob <- lapply(masks, function(m) apply(m, c(1L, 2L), mean)) # U x C each
  hits <- lapply(stats::setNames(names(masks), names(masks)), function(t) {
    function(p) prob[[t]] >= p
  })
  if (length(credible_union)) {
    stopifnot(all(credible_union %in% names(masks)))
    hits$any_credible <- function(p) {
      Reduce(`|`, lapply(credible_union, function(t) prob[[t]] >= p))
    }
  }
  for (t in names(hits)) {
    for (k in seq_len(nrow(sc))) {
      for (p in cutoffs) {
        rows[[length(rows) + 1L]] <- .point_scope_row(
          hits[[t]](p),
          dr,
          sc,
          k,
          enso_cond_cols(dr$conds),
          cutoff = p
        ) %>%
          dplyr::mutate(track = t)
      }
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::relocate(
      track,
      level,
      scope,
      enso_group,
      region,
      country,
      enso_label,
      enso_val,
      cutoff
    ) %>%
    add_pop_shares()
}

#' SUPERSEDED by report_enso_exposure direction `either` at
#' m = sensitive_irr, which reproduces its per-draw rows; kept side by side, no
#' longer a registered step (no ensure_*, no artefact kind).
#'
#' Population living in "ENSO-sensitive" ADM1 units: an ADM1 is sensitive at an
#' ONI condition on a posterior draw when |log IRR| > log(1 + SENSITIVE_IRR) for
#' a track (either direction: a protective effect counts as well as an adverse
#' one; the duration cessation-hazard ratio is used unflipped, which is
#' equivalent under the symmetric rule). Track "any" is the union over the
#' FIG1_TRACKS. Companion with a direction and a magnitude axis:
#' derive_report_enso_exposure().
#'
#' @return list(sensitive_pop = tibble(track, level, scope, enso_group, region,
#'   country, enso_label, enso_val, n_adm1_universe, n_pop_missing, pop_universe,
#'   n_adm1_point, pop_point, n_adm1_mean, pop_mean, pop_lo, pop_med, pop_hi, and
#'   the shares pct_point, pct_mean, pct_lo, pct_med, pct_hi of pop_universe),
#'   threshold, rule, tracks, enso_conds, n_draws, pop_ref_year, universe)
derive_report_sensitive_pop <- function(ctx) {
  thr <- ctx$SENSITIVE_IRR
  dr <- enso_draw_universe(ctx, "report-sensitive")
  conds <- dr$conds
  U <- dr$U

  masks <- enso_track_masks(dr, function(x, t) sensitive_mask(x, thr))
  masks$any <- Reduce(`|`, masks)
  # Point estimate from each track's posterior-mean ADM1 IRR (U x C).
  point <- lapply(dr$rel, function(r) {
    a <- r$adm1_irr %>% dplyr::mutate(adm1_gid = as.character(adm1_gid))
    m <- matrix(
      FALSE,
      length(U),
      length(conds),
      dimnames = list(U, names(conds))
    )
    ri <- match(a$adm1_gid, U)
    ci <- match(a$enso_label, names(conds))
    ok <- !is.na(ri) & !is.na(ci)
    m[cbind(ri[ok], ci[ok])] <- abs(log(a$IRR[ok])) > log1p(thr)
    m
  })
  point$any <- Reduce(`|`, point)

  sensitive_pop <- mask_pop_rows(masks, dr, point) %>%
    dplyr::select(
      track,
      level,
      scope,
      enso_group,
      region,
      country,
      enso_label,
      enso_val,
      n_adm1_universe,
      n_pop_missing,
      pop_universe,
      n_adm1_point,
      pop_point,
      n_adm1_mean,
      pop_mean,
      pop_lo,
      pop_med,
      pop_hi,
      pct_point,
      pct_mean,
      pct_lo,
      pct_med,
      pct_hi
    )

  list(
    sensitive_pop = sensitive_pop,
    threshold = thr,
    rule = sprintf(
      "|log IRR| > log(1 + %g) on any of %s (either direction)",
      thr,
      paste(dr$tracks, collapse = ", ")
    ),
    tracks = dr$tracks,
    enso_conds = conds,
    n_draws = dr$n_draws,
    pop_ref_year = ctx$POP_REF_YEAR,
    universe = dr$universe
  )
}

#' Freshness inputs shared by the per-draw ENSO population artefacts: manifest,
#' grouping, best-lag archive, survival panel, gid population and every FIG1
#' track's fit + cb_meta files.
enso_pop_inputs <- function(ctx) {
  tf <- lapply(ctx$FIG1_TRACKS, function(t) track_fit_input_files(ctx, t))
  results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$GROUPING_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      ctx$SURV_PATH,
      ctx$GID_POP_PATH,
      unlist(lapply(tf, `[[`, "fit")),
      unlist(lapply(tf, `[[`, "meta"))
    )
  )
}

# --- R6: ENSO exposure by direction, magnitude and certainty -----------------

#' Population exposed to an ENSO effect, by DIRECTION (adverse: more / longer
#' cholera; protective: the opposite; either: |log IRR| beyond the magnitude in
#' either direction, the duration ratio unflipped -- the former
#' report_sensitive_pop rule), MAGNITUDE (IRR beyond 1 + m, m in
#' ENSO_MAGNITUDES; m = 0 is "any effect in that direction") and, for the grid,
#' CERTAINTY (posterior probability >= p, p in POP_PROB_CUTOFFS). Each track's
#' adverse direction is TRACK_ADVERSE_SIGN (utils/results_core.R): occurrence
#' and monthly up, the duration cessation hazard DOWN (= longer outbreaks).
#' Per-track rows plus `any` (union across tracks, within the direction) and
#' `concordant2` (at least two tracks meet the rule on the draw; under `either`
#' each may meet it in its own direction); the grid also has `any_credible`
#' (the track-wise union of the per-track certain sets, prob_grid_rows()).
#'
#' Two readings of the same masks:
#'   * exposure -- the per-draw companion with a 95% CrI: on every draw sum
#'     the population of the ADM1s whose drawn effect exceeds m (mask_pop_rows);
#'     over the m axis this is the population exceedance curve. m = 0 is the
#'     directional counterpart of report_effects' per-draw dir_pop_ci() CrI;
#'     direction `either` at m = sensitive_irr reproduces the former
#'     report_sensitive_pop rows.
#'   * grid -- the probability-threshold count with no interval
#'     (prob_grid_rows): population where P(effect >= m) >= p. p = 0.95 at
#'     m = 0 is "significant" in the usual sense.
#'
#' @return list(exposure = tibble(direction, magnitude, track, level, scope,
#'   enso_group, region, country, enso_label, enso_val, n_adm1_universe,
#'   n_pop_missing, pop_universe, n_adm1_mean, pop_mean, pop_lo, pop_med, pop_hi,
#'   pct_*), grid = tibble(direction, magnitude, track, ..., cutoff, n_adm1, pop,
#'   pct), magnitudes, cutoffs, adverse_sign, tracks, enso_conds, n_draws,
#'   pop_ref_year, universe)
derive_report_enso_exposure <- function(ctx) {
  dr <- enso_draw_universe(ctx, "report-enso-exposure")
  mags <- ctx$ENSO_MAGNITUDES
  cutoffs <- ctx$POP_PROB_CUTOFFS
  sign_of <- TRACK_ADVERSE_SIGN[dr$tracks]
  if (anyNA(sign_of)) {
    stop(
      "[09] report_enso_exposure: no adverse direction registered for track(s) ",
      paste(dr$tracks[is.na(sign_of)], collapse = ", "),
      " (TRACK_ADVERSE_SIGN).",
      call. = FALSE
    )
  }
  exposure <- list()
  grid <- list()
  for (direction in ENSO_EXPOSURE_DIRECTIONS) {
    for (m in mags) {
      masks <- enso_combine_masks(enso_track_masks(dr, function(x, t) {
        exposure_mask(x, sign_of[[t]], direction, m)
      }))
      exposure[[length(exposure) + 1L]] <- mask_pop_rows(masks, dr) %>%
        dplyr::mutate(direction = direction, magnitude = m)
      grid[[length(grid) + 1L]] <- prob_grid_rows(
        masks,
        dr,
        cutoffs,
        credible_union = dr$tracks
      ) %>%
        dplyr::mutate(direction = direction, magnitude = m)
    }
  }
  list(
    exposure = dplyr::bind_rows(exposure) %>%
      dplyr::relocate(direction, magnitude),
    grid = dplyr::bind_rows(grid) %>% dplyr::relocate(direction, magnitude),
    magnitudes = mags,
    cutoffs = cutoffs,
    adverse_sign = sign_of,
    tracks = dr$tracks,
    enso_conds = dr$conds,
    n_draws = dr$n_draws,
    pop_ref_year = ctx$POP_REF_YEAR,
    universe = dr$universe
  )
}

ensure_report_enso_exposure <- function(ctx) {
  ensure_artefact(
    ctx,
    "report_enso_exposure",
    enso_pop_inputs(ctx),
    function() derive_report_enso_exposure(ctx)
  )
}

# --- R7: population in ADM1s with more ENSO-attributable cases ---------------
#
# Forecast arm vs counterfactual arm over the attribution window (the window and
# runs of report_forecast / Figure 4), per ADM1, paired by replicate_id (the
# arms share each replicate's smoothing handoff state and posterior draw). Two
# estimands per metric x analog set x magnitude m (thresholds.attrib_magnitudes):
#   * childs: the ADM1 counts when mean forecast / mean counterfactual >= 1 + m
#     (the impact-attribution convention of Childs et al. 2025, PNAS: a threshold
#     on the expected-incidence ratio), with a 95% interval on the population
#     from a paired bootstrap of the replicate ids within each group;
#   * childs_median: the same rule on the ratio of medians (insensitive to the
#     heavy upper tail of the replicate totals), also with a bootstrap interval;
#   * IPCC AR6 Atlas robustness classes from the paired replicates (the
#     magnitude m is the variability threshold): robust when >= robust_exceed of
#     the pairs have f > (1 + m) cf AND >= robust_agree have f > cf (ties
#     disagree); conflicting when the first holds but not the second; none
#     otherwise.
# Replicates excluded from Figure 4a's group ratio (results.thresholds.explosive_case_cap) are
# dropped from every ADM1 of their group and arm (explosive_guard = TRUE; the
# FALSE rows keep them, a sensitivity). Two windows: `attribution` (the
# inherited window, from the earlier arm launch, the factual series stitched
# from the smoothing before its own launch) and `post_fork` (from the factual
# launch: both arms free simulations; equal to `attribution` when the arms
# launch together). The headline slice is attribution x explosive_guard. Scope rows reuse the R5
# machinery (enso_scopes, mask_pop_rows, .point_scope_row) with the analog set
# as the condition axis and bootstrap resamples as the draw axis.

#' Bootstrap resamples drawn by attrib_exposure_tables(), and their fixed seed
#' (both stored in the artefact).
ATTRIB_BOOT_B <- 500L
#' The context's bootstrap resample count (results.derivation.attrib_boot_b), or
#' the default for a context without it (e.g. a test's).
.ctx_boot_b <- function(ctx) as.integer(ctx$BOOT_B %||% ATTRIB_BOOT_B)
ATTRIB_BOOT_SEED <- 20261023L

#' Evaluate `expr` under set.seed(seed), restoring the caller's RNG state.
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) {
    old <- get(".Random.seed", envir = globalenv())
  }
  on.exit(
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  set.seed(seed)
  expr
}

#' The Childs rule, shared by the point estimate and every bootstrap resample,
#' for both statistics (means, medians): the forecast statistic is positive and
#' at least (1 + m) times the counterfactual's. A zero counterfactual with a
#' positive forecast is an infinite ratio and counts; only 0 vs 0 does not
#' (dropping zero counterfactuals would remove most ADM1s from the ratio of
#' medians -- the typical counterfactual replicate has no outbreak there -- and
#' with them the strongest signals).
childs_hit <- function(stat_f, stat_cf, m) {
  stat_f > 0 & stat_f >= (1 + m) * stat_cf
}

#' A ratio of two non-negative statistics: Inf when only the denominator is 0
#' (the ADM1 counts), NA when both are 0 (no ratio).
.attrib_ratio <- function(num, den) {
  ifelse(den > 0, num / den, ifelse(num > 0, Inf, NA_real_))
}

#' Per-ADM1 summary of a pair_adm1_cum() table at magnitude m: n_pairs,
#' mean_f, mean_cf, ratio, the medians med_f, med_cf and ratio_median (ratios:
#' Inf when only the counterfactual is 0, NA when both are), p_exceed = share of pairs with
#' f > (1 + m) cf, p_sign = share with f > cf (strict: ties disagree), the
#' Childs hits (ratio of means, ratio of medians) and the IPCC class.
attrib_exposure_classes <- function(paired, m, robust_exceed, robust_agree) {
  paired %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_f = mean(f),
      mean_cf = mean(cf),
      med_f = stats::median(f),
      med_cf = stats::median(cf),
      p_exceed = mean(f > (1 + m) * cf),
      p_sign = mean(f > cf),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ratio = .attrib_ratio(mean_f, mean_cf),
      ratio_median = .attrib_ratio(med_f, med_cf),
      childs = childs_hit(mean_f, mean_cf, m),
      childs_median = childs_hit(med_f, med_cf, m),
      class = ipcc_class(p_exceed, p_sign, robust_exceed, robust_agree)
    )
}

#' The IPCC AR6 Atlas robustness class from the pair shares: "robust" when
#' >= robust_exceed of the pairs exceed (1 + m) and >= robust_agree agree in
#' sign, "conflicting" when only the first holds, else "none".
ipcc_class <- function(p_exceed, p_sign, robust_exceed, robust_agree) {
  dplyr::case_when(
    p_exceed >= robust_exceed & p_sign >= robust_agree ~ "robust",
    p_exceed >= robust_exceed ~ "conflicting",
    TRUE ~ "none"
  )
}

#' IPCC class of one scope's paired replicate totals at magnitude m (as
#' attrib_exposure_classes() per ADM1): 1-row tibble(n_pairs, p_exceed, p_sign,
#' class). Ties disagree (f > cf strictly).
ipcc_scope_stats <- function(f, cf, m, robust_exceed, robust_agree) {
  n <- length(f)
  pe <- if (n) mean(f > (1 + m) * cf) else NA_real_
  ps <- if (n) mean(f > cf) else NA_real_
  tibble::tibble(
    n_pairs = n,
    p_exceed = pe,
    p_sign = ps,
    class = if (n) {
      ipcc_class(pe, ps, robust_exceed, robust_agree)
    } else {
      NA_character_
    }
  )
}

#' The ADM1 universe of the forecast arm as the dr-shaped list the R5 scope
#' helpers read (universe, U, pop0). Population from build_forecast_adm1_pop()
#' (NA -> pop_missing, 0 in pop0); country = the ADM1's ISO3 prefix.
attrib_adm1_universe <- function(adm1_month, adm1_pop, group_short) {
  u <- adm1_month %>%
    dplyr::filter(run == "forecast") %>%
    dplyr::distinct(adm1_unit, enso_group) %>%
    dplyr::arrange(enso_group, adm1_unit)
  if (anyDuplicated(u$adm1_unit)) {
    stop(
      "[09] report_attrib_exposure: ADM1 unit(s) in more than one ENSO group: ",
      paste(
        head(unique(u$adm1_unit[duplicated(u$adm1_unit)]), 5L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  pop <- adm1_pop$pop[match(u$adm1_unit, adm1_pop$adm1_unit)]
  universe <- tibble::tibble(
    adm1_gid = u$adm1_unit,
    country = stringr::str_extract(u$adm1_unit, "^[A-Z]{3}"),
    enso_group = u$enso_group,
    region = unname(group_short[u$enso_group]),
    pop = pop,
    pop_missing = is.na(pop)
  )
  # pop_rate (rate_pop_year = "window"): the per-capita rates'
  # population; NULL (rates on `pop`) when adm1_pop has no such column.
  pop_rate <- if ("pop_rate" %in% names(adm1_pop)) {
    adm1_pop$pop_rate[match(u$adm1_unit, adm1_pop$adm1_unit)]
  } else {
    NULL
  }
  list(
    universe = universe,
    U = universe$adm1_gid,
    pop0 = ifelse(is.na(pop), 0, pop),
    pop_rate0 = if (is.null(pop_rate)) {
      NULL
    } else {
      ifelse(is.na(pop_rate), 0, pop_rate)
    }
  )
}

#' n x B multinomial resample counts (column b: how many times each of n paired
#' replicates is drawn in resample b); n = 0 gives a 0 x B matrix.
boot_weights <- function(n, B) {
  if (n == 0L) {
    return(matrix(0L, 0L, B))
  }
  matrix(
    vapply(
      seq_len(B),
      function(b) tabulate(sample.int(n, n, replace = TRUE), n),
      integer(n)
    ),
    n,
    B
  )
}

#' Paired cumulative totals as ADM1 x replicate matrices over the replicate ids
#' `ids` (rows `units`). A pair missing from `paired` is an error: the
#' simulator writes a full replicate x gid x month grid.
pair_adm1_matrices <- function(paired, units, ids) {
  i <- match(paired$adm1_unit, units)
  j <- match(paired$replicate_id, ids)
  keep <- !is.na(i) & !is.na(j)
  F <- CF <- matrix(NA_real_, length(units), length(ids))
  F[cbind(i[keep], j[keep])] <- paired$f[keep]
  CF[cbind(i[keep], j[keep])] <- paired$cf[keep]
  if (anyNA(F) || anyNA(CF)) {
    stop(
      "[09] report_attrib_exposure: incomplete forecast/counterfactual pairing for ADM1(s) ",
      paste(
        head(units[rowSums(is.na(F) | is.na(CF)) > 0], 5L),
        collapse = ", "
      ),
      " -- expected every replicate id on every ADM1 of the group.",
      call. = FALSE
    )
  }
  list(F = F, CF = CF)
}

#' Bootstrap statistics of the paired matrices, ADM1 x B each: the resampled
#' means (F %*% W / n) and medians (resample b repeats replicate j W[j, b]
#' times) of both arms. Thresholded per magnitude by childs_hit(); computed
#' once per (metric, set, group).
boot_adm1_stats <- function(mats, W) {
  A <- nrow(mats$F)
  n <- ncol(mats$F)
  B <- ncol(W)
  if (n == 0L) {
    z <- matrix(0, A, B)
    return(list(mean_f = z, mean_cf = z, med_f = z, med_cf = z))
  }
  rowmed <- function(M) {
    if (ncol(M) == 1L) M[, 1L] else apply(M, 1L, stats::median)
  }
  med_f <- med_cf <- matrix(NA_real_, A, B)
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), W[, b])
    med_f[, b] <- rowmed(mats$F[, idx, drop = FALSE])
    med_cf[, b] <- rowmed(mats$CF[, idx, drop = FALSE])
  }
  list(
    mean_f = mats$F %*% W / n,
    mean_cf = mats$CF %*% W / n,
    med_f = med_f,
    med_cf = med_cf
  )
}

#' U x S logical matrices of the per-ADM1 class table `cls` (columns adm1_unit,
#' analog_set, childs, class, means, medians, n_pairs): childs, robust,
#' conflicting, none, undefined / undefined_median (no ratio: forecast and
#' counterfactual statistic both 0) and paired (n_pairs > 0). An ADM1 absent
#' from `cls` for a set is FALSE everywhere.
attrib_class_hits <- function(cls, U, sets) {
  mk <- function(v) {
    out <- matrix(FALSE, length(U), length(sets), dimnames = list(U, sets))
    i <- match(cls$adm1_unit, U)
    j <- match(cls$analog_set, sets)
    keep <- !is.na(i) & !is.na(j)
    out[cbind(i[keep], j[keep])] <- v[keep]
    out
  }
  list(
    childs = mk(cls$childs),
    childs_median = mk(cls$childs_median),
    robust = mk(cls$class == "robust"),
    conflicting = mk(cls$class == "conflicting"),
    none = mk(cls$class == "none"),
    undefined = mk(cls$mean_f == 0 & cls$mean_cf == 0),
    undefined_median = mk(cls$med_f == 0 & cls$med_cf == 0),
    paired = mk(cls$n_pairs > 0)
  )
}

#' The pure core of report_attrib_exposure (no IO) for ONE window and one
#' explosive-replicate cap: the per-ADM1 table and the Africa / region / country
#' population rows of every estimand (childs, childs_median, robust,
#' conflicting, none).
#'
#' @param adm1_month forecast + counterfactual adm1_month rows over the window
#'   (load_group_arm_runs(); run, analog_set, enso_group, adm1_unit,
#'   replicate_id, month, cases, outbreaks)
#' @param country_month the matching country_month rows (for the explosive guard)
#' @param window_start first month of the window (the explosive guard sums from it)
#' @param cap explosive-replicate cap (Inf: keep every replicate)
#' @param adm1_pop build_forecast_adm1_pop()
#' @param group_short enso_group -> short region name (ctx$GROUP_SHORT_RUN)
#' @return list(exposure, adm1, universe, explosive)
attrib_exposure_tables <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  window_start,
  mags,
  robust_exceed,
  robust_agree,
  n_boot = ATTRIB_BOOT_B,
  boot_seed = ATTRIB_BOOT_SEED,
  cap
) {
  adm1_month <- dplyr::mutate(adm1_month, run = as.character(run))
  country_month <- dplyr::mutate(country_month, run = as.character(run))
  expl <- explosive_replicate_keys(country_month, window_start, cap)
  if (nrow(expl)) {
    message(sprintf(
      "[09] report_attrib_exposure: %d explosive replicate(s) (group window cases >= %g or a gid-month >= %g) dropped from their group's ADM1s",
      nrow(expl),
      .cap_case(cap),
      .cap_gid(cap)
    ))
  }
  kept <- .anti_join_rows(
    adm1_month,
    expl,
    by = c("enso_group", "run", "analog_set", "replicate_id")
  )
  rk <- .replicate_keys(kept)
  cum_all <- pair_adm1_cum_table(kept) # once per cell, every metric and set

  dr <- attrib_adm1_universe(adm1_month, adm1_pop, group_short)
  sc <- enso_scopes(dr$universe)
  groups <- unique(dr$universe$enso_group)
  units_g <- lapply(stats::setNames(groups, groups), function(g) {
    dr$U[dr$universe$enso_group == g]
  })
  cond <- tibble::tibble(analog_set = sets)

  # Paired replicate ids per (set, group) -- the explosive guard is cases-based,
  # so the ids are shared by both metrics -- and the bootstrap weights, drawn
  # once (loop order: set, then group) and shared by every metric and magnitude.
  ids <- W <- list()
  .with_seed(
    boot_seed,
    for (aset in sets) {
      for (g in groups) {
        sub <- rk[rk$enso_group == g, , drop = FALSE]
        id_g <- sort(intersect(
          unique(sub$replicate_id[sub$run == "forecast"]),
          unique(sub$replicate_id[
            sub$run == "counterfactual" & sub$analog_set %in% aset
          ])
        ))
        if (!length(id_g)) {
          warning(sprintf(
            "[09] report_attrib_exposure: group %s, set %s has no paired replicate left",
            g,
            aset
          ))
        }
        ids[[aset]][[g]] <- id_g
        W[[aset]][[g]] <- boot_weights(length(id_g), n_boot)
      }
    }
  )

  ADM1_METRIC <- c(cases = "cases", onsets = "outbreaks")
  BOOT_CLASSES <- c(childs = "mean", childs_median = "med")
  exposure <- list()
  adm1_rows <- list()
  for (metric in names(ADM1_METRIC)) {
    vc <- ADM1_METRIC[[metric]]
    paired <- lapply(stats::setNames(sets, sets), function(aset) {
      pair_adm1_cum(kept, vc, aset, cum = cum_all)
    })
    # Bootstrap means / medians once per (set, group); thresholded per m below.
    bst <- lapply(stats::setNames(sets, sets), function(aset) {
      lapply(groups, function(g) {
        boot_adm1_stats(
          pair_adm1_matrices(
            dplyr::filter(paired[[aset]], adm1_unit %in% units_g[[g]]),
            units_g[[g]],
            ids[[aset]][[g]]
          ),
          W[[aset]][[g]]
        )
      })
    })
    for (m in mags) {
      cls <- dplyr::bind_rows(lapply(sets, function(aset) {
        attrib_exposure_classes(
          paired[[aset]],
          m,
          robust_exceed,
          robust_agree
        ) %>%
          dplyr::mutate(analog_set = aset)
      }))
      hits <- attrib_class_hits(cls, dr$U, sets)

      boot <- lapply(BOOT_CLASSES, function(st) {
        arr <- array(FALSE, c(length(dr$U), length(sets), n_boot))
        for (s in seq_along(sets)) {
          for (gi in seq_along(groups)) {
            b <- bst[[sets[s]]][[gi]]
            arr[match(units_g[[gi]], dr$U), s, ] <- childs_hit(
              b[[paste0(st, "_f")]],
              b[[paste0(st, "_cf")]],
              m
            )
          }
        }
        arr
      })
      # unname(): the point matrices' dimnames would otherwise name the
      # colSums() of n_adm1_point / pop_point.
      childs_rows <- mask_pop_rows(
        boot,
        dr,
        point = lapply(hits[names(BOOT_CLASSES)], unname),
        cond_cols = cond
      ) %>%
        dplyr::rename(class = track) %>%
        dplyr::mutate(n_adm1 = n_adm1_point, pop = pop_point, pct = pct_point)
      class_rows <- dplyr::bind_rows(lapply(
        c("robust", "conflicting", "none"),
        function(cl) {
          dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
            .point_scope_row(hits[[cl]], dr, sc, k, cond)
          })) %>%
            dplyr::mutate(class = cl)
        }
      )) %>%
        add_pop_shares()
      counts <- dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
        idx <- sc$idx[[k]]
        tibble::tibble(
          level = sc$level[k],
          scope = sc$scope[k],
          analog_set = sets,
          und_mean = unname(colSums(hits$undefined[idx, , drop = FALSE])),
          und_med = unname(colSums(hits$undefined_median[idx, , drop = FALSE])),
          n_adm1_paired = unname(colSums(hits$paired[idx, , drop = FALSE]))
        )
      }))
      exposure[[length(exposure) + 1L]] <- dplyr::bind_rows(
        childs_rows,
        class_rows
      ) %>%
        dplyr::left_join(counts, by = c("level", "scope", "analog_set")) %>%
        dplyr::mutate(
          n_adm1_undefined = ifelse(
            class == "childs_median",
            und_med,
            und_mean
          ),
          metric = metric,
          magnitude = m
        ) %>%
        dplyr::select(-und_mean, -und_med)
      adm1_rows[[length(adm1_rows) + 1L]] <- cls %>%
        dplyr::left_join(
          dplyr::select(
            dr$universe,
            adm1_unit = adm1_gid,
            country,
            enso_group,
            region,
            pop
          ),
          by = "adm1_unit"
        ) %>%
        dplyr::mutate(metric = metric, magnitude = m)
    }
  }
  exposure <- dplyr::bind_rows(exposure) %>%
    dplyr::relocate(
      metric,
      magnitude,
      class,
      level,
      scope,
      enso_group,
      region,
      country,
      analog_set,
      n_adm1,
      pop,
      pct
    )
  adm1 <- dplyr::bind_rows(adm1_rows) %>%
    dplyr::relocate(
      adm1_unit,
      country,
      enso_group,
      region,
      metric,
      analog_set,
      magnitude,
      n_pairs,
      mean_f,
      mean_cf,
      ratio,
      med_f,
      med_cf,
      ratio_median,
      p_exceed,
      p_sign,
      childs,
      childs_median,
      class,
      pop
    )
  list(
    exposure = exposure,
    adm1 = adm1,
    universe = dr$universe,
    explosive = expl
  )
}

#' attrib_exposure_tables() over the attribution cells (attrib_cells(): the
#' attribution window from the attribution start, post_fork from the displayed
#' factual's launch, each crossed with the explosive modes), rows tagged
#' `window` and `explosive_guard`. The frames are cut to each cell's months
#' (.cell_frame()); the bootstrap seed is the same in every cell.
#' @param cells attrib_cell_table()
#' @return list(exposure, adm1, universe, explosive (the dropped replicates per
#'   window, guarded rows only), windows = named Date starts)
attrib_exposure_windows <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  cells,
  mags,
  robust_exceed,
  robust_agree,
  n_boot = ATTRIB_BOOT_B,
  boot_seed = ATTRIB_BOOT_SEED
) {
  out <- list()
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    o <- attrib_exposure_tables(
      .cell_frame(adm1_month, cl),
      .cell_frame(country_month, cl),
      adm1_pop,
      group_short,
      sets,
      cl$start,
      mags,
      robust_exceed,
      robust_agree,
      n_boot = n_boot,
      boot_seed = boot_seed,
      cap = .cell_cap(cl)
    )
    tag <- function(d) .cell_tag(d, cl)
    out[[length(out) + 1L]] <- list(
      exposure = tag(o$exposure),
      adm1 = tag(o$adm1),
      explosive = tag(o$explosive),
      universe = o$universe
    )
  }
  list(
    exposure = dplyr::bind_rows(lapply(out, `[[`, "exposure")),
    adm1 = dplyr::bind_rows(lapply(out, `[[`, "adm1")),
    explosive = dplyr::bind_rows(lapply(out, `[[`, "explosive")) %>%
      dplyr::filter(explosive_guard),
    universe = out[[1L]]$universe,
    windows = .cell_starts(cells)
  )
}

#' Inputs shared by report_attrib_exposure and report_attrib_metrics: the
#' included groups, the arm window, the forecast + counterfactual frames of every
#' (group, set) run (adm1_month, country_month) and the observed ADM1 population.
#' @param tag label of build_adm1_obs_pop()'s messages
#' @param factual_source "simulation" (default: the factual run, stitched from
#'   the smoothing before its launch as before) or "smoothing" (event
#'   attribution: every factual month read from the smoothing replicate each
#'   forecast replicate was handed off from; needs `smoothing_events`).
#' @param smoothing_events read_smoothing_parts_events() output for
#'   factual_source = "smoothing" (read once by the caller across events).
#' @param base NULL, or list(gid_pop, surv_dat, frame) already loaded (the event
#'   attribution calls this once per event and factual source).
attrib_arm_inputs <- function(
  ctx,
  sets,
  tag,
  factual_source = "simulation",
  smoothing_events = NULL,
  base = NULL
) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  smooth_fact <- identical(factual_source, "smoothing")
  win <- results_arm_window(ctx, groups, sets, read_events = !smooth_fact)
  if (smooth_fact) {
    if (is.null(smoothing_events)) {
      stop(
        "[09] attrib_arm_inputs: factual_source = 'smoothing' needs smoothing_events.",
        call. = FALSE
      )
    }
    win$factual_source <- "smoothing"
    win$smoothing_events <- smoothing_events
  }
  attrib_start <- win$attrib_start
  base <- base %||% attrib_base_inputs(ctx)
  runs <- load_group_arm_runs(ctx, groups, sets, attrib_start, win = win)
  bf <- .bound_run_frames(ctx, runs)
  adm1_month_both <- bf$adm1_month
  country_month_both <- bf$country_month
  adm1_pop <- build_forecast_adm1_pop(
    ctx,
    base$frame,
    base$surv_dat,
    base$gid_pop,
    adm1_month_both,
    tag
  )
  list(
    man = man,
    groups = groups,
    win = win,
    attrib_start = attrib_start,
    adm1_month = adm1_month_both,
    country_month = country_month_both,
    adm1_pop = adm1_pop,
    gids = runs$gids
  )
}

#' The inputs attrib_arm_inputs() reads besides the runs: gid population,
#' survival panel and location frame.
attrib_base_inputs <- function(ctx) {
  gid_pop <- load_gid_population(ctx$GID_POP_PATH)
  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  frame <- build_location_frame(
    survival_data = surv_dat,
    hmm_clusters_path = ctx$HMM_PATH,
    outbreak_metadata_path = NULL,
    size_data_path = ctx$SIZE_DATA_PATH
  )
  list(gid_pop = gid_pop, surv_dat = surv_dat, frame = frame)
}

derive_report_attrib_exposure <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  message(sprintf(
    "[09] report_attrib_exposure: groups=%s | sets=%s | m=%s | robust %g/%g | B=%d",
    paste(groups, collapse = ","),
    paste(sets, collapse = ","),
    paste(ctx$ATTRIB_MAGNITUDES, collapse = ","),
    ctx$ROBUST_EXCEED,
    ctx$ROBUST_AGREE,
    .ctx_boot_b(ctx)
  ))
  inp <- attrib_arm_inputs(ctx, sets, "report-attrib-exposure")
  win <- inp$win
  attrib_start <- inp$attrib_start
  post_fork_start <- display_factual_start(ctx, groups, sets)
  cells <- attrib_cells(ctx, attrib_start, post_fork_start)
  af_posterior_stash_put(
    ctx,
    paste("forecast", paste(sort(sets), collapse = "-")),
    inp,
    cells
  )

  out <- attrib_exposure_windows(
    inp$adm1_month,
    inp$country_month,
    inp$adm1_pop,
    ctx$GROUP_SHORT_RUN,
    sets,
    cells,
    ctx$ATTRIB_MAGNITUDES,
    ctx$ROBUST_EXCEED,
    ctx$ROBUST_AGREE,
    n_boot = .ctx_boot_b(ctx)
  )
  c(
    out,
    list(
      magnitudes = ctx$ATTRIB_MAGNITUDES,
      robust_exceed = ctx$ROBUST_EXCEED,
      robust_agree = ctx$ROBUST_AGREE,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      n_boot = .ctx_boot_b(ctx),
      boot_seed = ATTRIB_BOOT_SEED,
      attrib_start = attrib_start,
      africa_min_horizon = as.Date(man$africa_min_horizon),
      pop_ref_year = ctx$POP_REF_YEAR,
      analog_sets = sets,
      groups = groups,
      factual_start = win$factual_start,
      arm_starts = win$arm_starts,
      post_fork_start = post_fork_start,
      attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
      cells = cells_payload(cells),
      africa_window = africa_window_payload(
        .africa_win(ctx, groups),
        attrib_start
      )
    )
  ) # report_attrib_exposure v4
}

ensure_report_attrib_exposure <- function(ctx, sets) {
  ctx <- results_attrib_ctx(ctx)
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- results_run_configs(man, groups, sets)
  inputs <- results_inputs(
    files = c(
      results_attrib_manifest_files(ctx),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = if (
      results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched
    ) {
      smoothing_parts_glob(ctx)
    } else {
      list()
    }
  )
  ensure_artefact(
    ctx,
    "report_attrib_exposure",
    inputs,
    function() derive_report_attrib_exposure(ctx, sets),
    sets = sets
  )
}

# ===========================================================================
# report_attrib_incidence: window totals and annualised incidence
# of the paired forecast / ENSO-counterfactual runs of report_attrib_exposure,
# for Africa, each group AND each country, in every attribution cell (window x
# explosive mode). Supplementary table of forecast vs counterfactual incidence
# by country (the attribution artefact's cum_totals_cells stop at the groups).
# Africa and group rows are the cum_totals_cells rows of the attribution
# artefact (same replicates, same explosive guard); country rows drop the
# GROUP-level explosive replicates (country_cum_kept(), as Figure 4a).
# Incidence annualises the window total as report_forecast does:
# tot / n_months * 12 / pop * 1e5, n_months the months the scope's group
# simulates in the cell (Africa: the cell's months), pop the simulated gids'
# population (build_sim_scope_pop(); rate_pop_year "window": month-weighted
# over the cell's months).
# ===========================================================================

#' Per-replicate window totals of one cell, long: (level, scope, enso_group,
#' country, run, analog_set, replicate_id, cases, onsets). Group and Africa rows
#' as the attribution step's cum_totals_of(): per set, Africa sums the replicate
#' ids present in every group on both arms.
attrib_incidence_per_rep <- function(cg, cc, sets) {
  dplyr::bind_rows(lapply(sets, function(aset) {
    grp <- cg %>%
      dplyr::filter(run == "forecast" | analog_set == aset) %>%
      dplyr::transmute(
        enso_group,
        run = as.character(run),
        replicate_id,
        cases = cum_cases,
        onsets = cum_outbreaks
      )
    common <- Reduce(
      intersect,
      lapply(split(grp, grp$run), function(a) {
        Reduce(
          intersect,
          lapply(split(a, a$enso_group), function(d) unique(d$replicate_id))
        )
      })
    )
    afr <- grp %>%
      dplyr::filter(replicate_id %in% common) %>%
      dplyr::group_by(run, replicate_id) %>%
      dplyr::summarise(
        cases = sum(cases),
        onsets = sum(onsets),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        level = "continent",
        scope = "Africa",
        enso_group = NA_character_,
        country = NA_character_
      )
    ctry <- cc %>%
      dplyr::filter(run == "forecast" | analog_set == aset) %>%
      dplyr::transmute(
        level = "country",
        scope = country,
        enso_group,
        country,
        run = as.character(run),
        replicate_id,
        cases = cum_cases,
        onsets = cum_outbreaks
      )
    dplyr::bind_rows(
      afr,
      dplyr::mutate(
        grp,
        level = "region",
        scope = enso_group,
        country = NA_character_
      ),
      ctry
    ) %>%
      dplyr::mutate(analog_set = aset)
  }))
}

#' Window totals and annualised incidence of one cell: median, q025 / q975
#' (and q25 / q75) across replicates per scope x arm x set x metric.
#' @param per_rep attrib_incidence_per_rep() rows
#' @param scope_pop build_sim_scope_pop() rows (level, scope, enso_group, pop)
#' @param n_months named integer: months simulated in the cell per group, and
#'   "Africa"
attrib_incidence_table <- function(per_rep, scope_pop, n_months) {
  pop_key <- scope_pop %>%
    dplyr::transmute(
      level,
      key = dplyr::case_when(
        level == "continent" ~ "Africa",
        level == "region" ~ enso_group,
        TRUE ~ scope
      ),
      pop
    )
  qv <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
  per_rep %>%
    tidyr::pivot_longer(
      c(cases, onsets),
      names_to = "metric",
      values_to = "total"
    ) %>%
    dplyr::group_by(
      level,
      scope,
      enso_group,
      country,
      run,
      analog_set,
      metric
    ) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      tot_lo = qv(total, .025),
      tot_q25 = qv(total, .25),
      tot_med = qv(total, .5),
      tot_q75 = qv(total, .75),
      tot_hi = qv(total, .975),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      key = ifelse(level == "continent", "Africa", scope),
      n_months = unname(n_months[ifelse(
        level == "continent",
        "Africa",
        enso_group
      )])
    ) %>%
    dplyr::left_join(pop_key, by = c("level", "key")) %>%
    # Quantiles of the annualised rate are the scaled quantiles of the total
    # (a fixed positive factor per scope); NA without a population.
    dplyr::mutate(
      f = dplyr::if_else(
        is.na(pop) | pop <= 0,
        NA_real_,
        12 / n_months / pop * 1e5
      ),
      inc_lo = tot_lo * f,
      inc_q25 = tot_q25 * f,
      inc_med = tot_med * f,
      inc_q75 = tot_q75 * f,
      inc_hi = tot_hi * f
    ) %>%
    dplyr::select(-key, -f)
}

derive_report_attrib_incidence <- function(ctx, sets) {
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  message(sprintf(
    "[09] report_attrib_incidence: groups=%s | sets=%s | rate_pop_year=%s",
    paste(groups, collapse = ","),
    paste(sets, collapse = ","),
    ctx[["RATE_POP_YEAR"]] %||% "ref"
  ))
  inp <- attrib_arm_inputs(ctx, sets, "report-attrib-incidence")
  win <- inp$win
  attrib_start <- inp$attrib_start
  post_fork_start <- display_factual_start(ctx, groups, sets)
  cells <- attrib_cells(ctx, attrib_start, post_fork_start)
  gid_pop <- load_gid_population(ctx$GID_POP_PATH)
  wpp <- if (identical(ctx[["RATE_POP_YEAR"]], "window")) {
    load_wpp_table(.pop_artefact_or_null(ctx$WPP_PATH))
  } else {
    NULL
  }

  rows <- lapply(seq_len(nrow(cells)), function(i) {
    cl <- cells[i, ]
    cm <- .cell_frame(inp$country_month, cl)
    cg <- drop_explosive_replicates(
      window_cum_per_rep(cm, cl$start, "enso_group"),
      .cell_cap(cl)
    )
    cc <- country_cum_kept(cm, cl$start, .cell_cap(cl))
    n_months <- c(
      vapply(
        groups,
        function(g) dplyr::n_distinct(cm$month[cm$enso_group == g]),
        integer(1)
      ),
      Africa = dplyr::n_distinct(cm$month)
    )
    scope_pop <- if (identical(ctx[["RATE_POP_YEAR"]], "window")) {
      build_sim_scope_pop(
        inp$gids,
        gid_pop,
        ctx$POP_REF_YEAR,
        ctx$GROUP_SHORT_RUN,
        pop_year_weights = rate_pop_weights(cm$month),
        wpp = wpp
      )
    } else {
      build_sim_scope_pop(
        inp$gids,
        gid_pop,
        ctx$POP_REF_YEAR,
        ctx$GROUP_SHORT_RUN
      )
    }
    tag <- function(d) {
      .cell_tag(d, cl) %>%
        dplyr::mutate(
          explosive_mode = cl$explosive_mode[[1]],
          .after = explosive_guard
        )
    }
    list(
      incidence = tag(attrib_incidence_table(
        attrib_incidence_per_rep(cg, cc, sets),
        scope_pop,
        n_months
      )),
      scope_pop = tag(scope_pop)
    )
  })
  list(
    incidence = dplyr::bind_rows(lapply(rows, `[[`, "incidence")),
    scope_pop = dplyr::bind_rows(lapply(rows, `[[`, "scope_pop")),
    cells = cells_payload(cells),
    windows = .cell_starts(cells),
    explosive_cap = ctx$EXPLOSIVE_CAP,
    explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
    explosive_modes = .ctx_explosive_modes(ctx),
    attrib_start = attrib_start,
    pop_ref_year = ctx$POP_REF_YEAR,
    rate_pop_year = ctx[["RATE_POP_YEAR"]] %||% "ref",
    analog_sets = sets,
    groups = groups,
    factual_start = win$factual_start,
    arm_starts = win$arm_starts,
    post_fork_start = post_fork_start,
    attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
    africa_window = africa_window_payload(
      .africa_win(ctx, groups),
      attrib_start
    )
  ) # report_attrib_incidence v1
}

ensure_report_attrib_incidence <- function(ctx, sets) {
  ctx <- results_attrib_ctx(ctx)
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- results_run_configs(man, groups, sets)
  inputs <- results_inputs(
    files = c(
      results_attrib_manifest_files(ctx),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = if (
      results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched
    ) {
      smoothing_parts_glob(ctx)
    } else {
      list()
    }
  )
  ensure_artefact(
    ctx,
    "report_attrib_incidence",
    inputs,
    function() derive_report_attrib_incidence(ctx, sets),
    sets = sets
  )
}

# ===========================================================================
# report_attrib_metrics: attribution metrics of the paired
# forecast / ENSO-counterfactual runs of report_attrib_exposure, beside its
# Childs and IPCC-robustness classes. Three families, so the manuscript can pick
# one later:
#   * attributable number / fraction (AN = F - CF, AF = AN / F; Gasparrini &
#     Leone 2014, Vicedo-Cabrera et al. 2021, Childs et al. 2025) of cases and
#     onsets, per replicate and of the replicate means;
#   * probability ratio PR = p1 / p0 of an outbreak event (World Weather
#     Attribution) with FAR = 1 - 1/PR (= the probability of necessary
#     causation PN), PS = (p1 - p0) / (1 - p0) and PNS = p1 - p0 (Hannart &
#     Naveau), for two events: onset_any (an onset in the window) and cases_ge
#     (window cases >= thresholds.pr_case_threshold); when the counterfactual
#     never has the event PR is infinite and pr_lo is the exact one-sided bound
#     p1 / (1 - alpha^(1/n)) (Paciorek et al.);
#   * signal-to-noise of the paired difference, mean / sd across replicates
#     (Hawkins & Sutton), and the population of the ADM1s where |S/N| reaches
#     thresholds.sn_levels.
# Same axes as report_attrib_exposure: analog sets x window x explosive guard,
# ADM1 + country + region + Africa. Scope totals are sums of the ADM1 rows per
# replicate (country / region within one group's paired ids; Africa over the
# ids common to every group), so ADM1 and scope rows are consistent. Bootstrap
# weights: the same draws as attrib_exposure_tables() per (set, group), plus one
# resample of the common ids for Africa.
# ===========================================================================

#' Type-1 (no interpolation) quantiles, safe with Inf (type 7 gives NaN on
#' 0 * Inf); NA when nothing is finite or infinite to rank.
.boot_q <- function(x, p) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(rep(NA_real_, length(p)))
  }
  stats::quantile(x, p, type = 1, names = FALSE)
}

#' Replicate quantiles q025 / med / q975 (NA-safe).
.q3 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(q025 = NA_real_, med = NA_real_, q975 = NA_real_))
  }
  stats::setNames(
    stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE),
    c("q025", "med", "q975")
  )
}

#' Attributable number / fraction of paired per-replicate totals.
#' @param f,cf forecast / counterfactual totals of the n paired replicates
#' @param W n x B bootstrap resample counts (boot_weights())
#' @param pop the scope's population (for an_per_100k)
#' @return 1-row tibble
af_stats <- function(f, cf, W, pop = NA_real_) {
  n <- length(f)
  an_r <- f - cf
  af_r <- ifelse(f > 0, an_r / f, NA_real_)
  mean_f <- if (n) mean(f) else NA_real_
  mean_cf <- if (n) mean(cf) else NA_real_
  an_mean <- mean_f - mean_cf
  af_m <- ifelse(isTRUE(mean_f > 0), an_mean / mean_f, NA_real_)
  if (n && ncol(W)) {
    mf_b <- as.numeric(crossprod(W, f)) / n
    mcf_b <- as.numeric(crossprod(W, cf)) / n
    an_b <- mf_b - mcf_b
    af_b <- ifelse(mf_b > 0, an_b / mf_b, NA_real_)
  } else {
    an_b <- af_b <- numeric(0)
  }
  qa <- .q3(an_r)
  qf <- .q3(af_r)
  ba <- .boot_q(an_b, c(.025, .975))
  bf <- .boot_q(af_b, c(.025, .975))
  tibble::tibble(
    n_pairs = n,
    n_f_zero = sum(f == 0),
    mean_f = mean_f,
    mean_cf = mean_cf,
    an_mean = an_mean,
    an_lo = ba[1],
    an_hi = ba[2],
    an_q025 = qa[["q025"]],
    an_med = qa[["med"]],
    an_q975 = qa[["q975"]],
    af_of_means = af_m,
    af_lo = bf[1],
    af_hi = bf[2],
    af_q025 = qf[["q025"]],
    af_med = qf[["med"]],
    af_q975 = qf[["q975"]],
    an_per_100k = ifelse(
      is.finite(pop) & pop > 0,
      an_mean / pop * 1e5,
      NA_real_
    )
  )
}

# PR family from event probabilities (vectorised over resamples).
.pr_family <- function(p1, p0) {
  pr <- .attrib_ratio(p1, p0)
  list(
    pr = pr,
    far = ifelse(is.na(pr), NA_real_, ifelse(is.infinite(pr), 1, 1 - 1 / pr)),
    ps = ifelse(p0 < 1, (p1 - p0) / (1 - p0), NA_real_),
    pns = p1 - p0
  )
}

#' Probability ratio and probabilities of causation of a paired event.
#' @param ef,ecf 0/1 event indicators of the n paired replicates (forecast, counterfactual)
#' @param W n x B resample counts
#' @param alpha one-sided level of the exact bound on p0 when no counterfactual
#'   replicate has the event
#' @return 1-row tibble: n_pairs, p1, p0, p0_zero, p0_upper, pr (+ lo/hi), far,
#'   ps, pns (+ lo/hi each)
pr_stats <- function(ef, ecf, W, alpha) {
  n <- length(ef)
  if (!n) {
    return(tibble::tibble(
      n_pairs = 0L,
      p1 = NA_real_,
      p0 = NA_real_,
      p0_zero = NA,
      p0_upper = NA_real_,
      pr = NA_real_,
      pr_lo = NA_real_,
      pr_hi = NA_real_,
      far = NA_real_,
      far_lo = NA_real_,
      far_hi = NA_real_,
      ps = NA_real_,
      ps_lo = NA_real_,
      ps_hi = NA_real_,
      pns = NA_real_,
      pns_lo = NA_real_,
      pns_hi = NA_real_
    ))
  }
  p1 <- mean(ef)
  p0 <- mean(ecf)
  pt <- .pr_family(p1, p0)
  p0_upper <- 1 - alpha^(1 / n)
  p1_b <- as.numeric(crossprod(W, ef)) / n
  p0_b <- as.numeric(crossprod(W, ecf)) / n
  bt <- .pr_family(p1_b, p0_b)
  ci <- lapply(bt, .boot_q, p = c(.025, .975))
  p0_zero <- p0 == 0 && p1 > 0
  # p0 = 0: the resampled PR is Inf in every resample with p1_b > 0; its lower
  # end is the exact bound (p0 <= p0_upper with confidence 1 - alpha).
  pr_lo <- if (p0_zero) p1 / p0_upper else ci$pr[1]
  tibble::tibble(
    n_pairs = n,
    p1 = p1,
    p0 = p0,
    p0_zero = p0_zero,
    p0_upper = p0_upper,
    pr = pt$pr,
    pr_lo = pr_lo,
    pr_hi = if (is.na(pt$pr)) NA_real_ else ci$pr[2],
    far = pt$far,
    far_lo = ci$far[1],
    far_hi = ci$far[2],
    ps = pt$ps,
    ps_lo = ci$ps[1],
    ps_hi = ci$ps[2],
    pns = pt$pns,
    pns_lo = ci$pns[1],
    pns_hi = ci$pns[2]
  )
}

#' report_event_attrib: stamp a cell table with its past event (and sub-window,
#' factual source). The probability-ratio tables already carry `event`, the
#' outbreak event of the ratio (onset_any / cases_ge, results.headline.pr_event);
#' it is renamed `pr_event` first, as the past event's name would overwrite it
#' and leave two indistinguishable rows per cell (report_event_attrib v4).
#' @noRd
.event_attrib_stamp <- function(d, ...) {
  if ("event" %in% names(d)) {
    d <- dplyr::rename(d, pr_event = "event")
  }
  dplyr::mutate(d, ..., .before = 1L)
}

#' Signal-to-noise of the paired difference d = f - cf: mean(d) / sd(d) (NA
#' when sd is 0 or fewer than 2 pairs), with p_sign = share of pairs with f > cf.
sn_stats <- function(f, cf) {
  d <- f - cf
  n <- length(d)
  sdd <- if (n >= 2L) stats::sd(d) else NA_real_
  tibble::tibble(
    n_pairs = n,
    mean_diff = if (n) mean(d) else NA_real_,
    sd_diff = sdd,
    sn = ifelse(is.finite(sdd) & sdd > 0, mean(d) / sdd, NA_real_),
    p_sign = if (n) mean(d > 0) else NA_real_
  )
}

#' Per-replicate F / CF totals of a scope (row k of enso_scopes()): sums of the
#' scope's ADM1 rows in the group matrices. A country / region lies within one
#' group and uses its paired ids; Africa sums every group over `africa_ids`.
#' @param mats named list (by group) of pair_adm1_matrices() results
#' @param ids named list (by group) of the paired replicate ids (matrix columns)
#' @param units_g named list (by group) of the ADM1 units (matrix rows)
#' @return list(f, cf, group) (group NA for Africa)
scope_pair_vectors <- function(mats, ids, units_g, dr, sc, k, africa_ids) {
  rows_u <- dr$U[sc$idx[[k]]]
  if (sc$level[k] == "continent") {
    f <- cf <- numeric(length(africa_ids))
    for (g in names(mats)) {
      j <- match(africa_ids, ids[[g]])
      f <- f + colSums(mats[[g]]$F[, j, drop = FALSE])
      cf <- cf + colSums(mats[[g]]$CF[, j, drop = FALSE])
    }
    return(list(f = unname(f), cf = unname(cf), group = NA_character_))
  }
  g <- sc$enso_group[k]
  r <- match(rows_u, units_g[[g]])
  list(
    f = unname(colSums(mats[[g]]$F[r, , drop = FALSE])),
    cf = unname(colSums(mats[[g]]$CF[r, , drop = FALSE])),
    group = g
  )
}

#' Pure core of report_attrib_metrics for one window and one explosive cap
#' (mirrors attrib_exposure_tables(); same inputs). `ipcc_mags` (v2): the
#' magnitudes of the scope-level IPCC class table `ipcc_scope` (Africa, ENSO
#' regions, countries: the class of each scope's summed paired totals, as
#' attrib_exposure_classes() per ADM1); NULL derives none (the event cores).
#' @return list(af_adm1, af_scope, pr_adm1, pr_scope, pr_pop, sn_adm1, sn_scope,
#'   sn_pop, [ipcc_scope,] universe, explosive)
attrib_metrics_tables <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  window_start,
  pr_case_threshold,
  pr_alpha,
  sn_levels,
  n_boot = ATTRIB_BOOT_B,
  boot_seed = ATTRIB_BOOT_SEED,
  cap,
  ipcc_mags = NULL,
  robust_exceed = NA_real_,
  robust_agree = NA_real_
) {
  adm1_month <- dplyr::mutate(adm1_month, run = as.character(run))
  country_month <- dplyr::mutate(country_month, run = as.character(run))
  expl <- explosive_replicate_keys(country_month, window_start, cap)
  kept <- .anti_join_rows(
    adm1_month,
    expl,
    by = c("enso_group", "run", "analog_set", "replicate_id")
  )
  rk <- .replicate_keys(kept)
  cum_all <- pair_adm1_cum_table(kept) # once per cell, every metric and set

  dr <- attrib_adm1_universe(adm1_month, adm1_pop, group_short)
  sc <- enso_scopes(dr$universe)
  groups <- unique(dr$universe$enso_group)
  units_g <- lapply(stats::setNames(groups, groups), function(g) {
    dr$U[dr$universe$enso_group == g]
  })
  cond <- tibble::tibble(analog_set = sets)
  adm1_cols <- dplyr::select(
    dr$universe,
    adm1_unit = adm1_gid,
    country,
    enso_group,
    region,
    pop
  )

  # Paired ids and bootstrap weights: per (set, group) exactly as
  # attrib_exposure_tables() draws them, then one resample of the ids common to
  # every group for Africa.
  ids <- W <- afr_ids <- W_afr <- list()
  .with_seed(boot_seed, {
    for (aset in sets) {
      for (g in groups) {
        sub <- rk[rk$enso_group == g, , drop = FALSE]
        id_g <- sort(intersect(
          unique(sub$replicate_id[sub$run == "forecast"]),
          unique(sub$replicate_id[
            sub$run == "counterfactual" & sub$analog_set %in% aset
          ])
        ))
        ids[[aset]][[g]] <- id_g
        W[[aset]][[g]] <- boot_weights(length(id_g), n_boot)
      }
    }
    for (aset in sets) {
      afr_ids[[aset]] <- sort(Reduce(intersect, ids[[aset]]))
      W_afr[[aset]] <- boot_weights(length(afr_ids[[aset]]), n_boot)
    }
  })

  # Matrices per (metric, set, group): ADM1 x paired replicate.
  VC <- c(cases = "cases", onsets = "outbreaks")
  mats <- lapply(VC, function(vc) {
    lapply(stats::setNames(sets, sets), function(aset) {
      paired <- pair_adm1_cum(kept, vc, aset, cum = cum_all)
      lapply(stats::setNames(groups, groups), function(g) {
        pair_adm1_matrices(
          dplyr::filter(paired, adm1_unit %in% units_g[[g]]),
          units_g[[g]],
          ids[[aset]][[g]]
        )
      })
    })
  })

  af_a <- af_s <- sn_a <- sn_s <- pr_a <- pr_s <- ip_s <- list()
  sn_hit <- list() # metric -> U x S matrix of sn
  pr_hit <- list() # event -> list(point, lower) U x S matrices
  EVENTS <- list(
    onset_any = list(metric = "onsets", fn = function(x) x > 0),
    cases_ge = list(metric = "cases", fn = function(x) x >= pr_case_threshold)
  )
  for (metric in names(VC)) {
    sn_hit[[metric]] <- matrix(NA_real_, length(dr$U), length(sets))
  }
  for (ev in names(EVENTS)) {
    pr_hit[[ev]] <- list(
      point = matrix(FALSE, length(dr$U), length(sets)),
      lower = matrix(FALSE, length(dr$U), length(sets))
    )
  }

  for (s in seq_along(sets)) {
    aset <- sets[s]
    for (g in groups) {
      Wg <- W[[aset]][[g]]
      for (i in seq_along(units_g[[g]])) {
        u <- units_g[[g]][i]
        ui <- match(u, dr$U)
        key <- dplyr::filter(adm1_cols, adm1_unit == u)
        for (metric in names(VC)) {
          m <- mats[[metric]][[aset]][[g]]
          f <- m$F[i, ]
          cf <- m$CF[i, ]
          af_a[[length(af_a) + 1L]] <- dplyr::bind_cols(
            key,
            tibble::tibble(metric = metric, analog_set = aset),
            af_stats(
              f,
              cf,
              Wg,
              if (is.null(dr$pop_rate0)) {
                key$pop
              } else {
                dr$pop_rate0[ui]
              }
            )
          )
          sn <- sn_stats(f, cf)
          sn_hit[[metric]][ui, s] <- sn$sn
          sn_a[[length(sn_a) + 1L]] <- dplyr::bind_cols(
            key,
            tibble::tibble(metric = metric, analog_set = aset),
            sn
          )
        }
        for (ev in names(EVENTS)) {
          m <- mats[[EVENTS[[ev]]$metric]][[aset]][[g]]
          pr <- pr_stats(
            as.numeric(EVENTS[[ev]]$fn(m$F[i, ])),
            as.numeric(EVENTS[[ev]]$fn(m$CF[i, ])),
            Wg,
            pr_alpha
          )
          pr_hit[[ev]]$point[ui, s] <- isTRUE(pr$pr > 1)
          pr_hit[[ev]]$lower[ui, s] <- isTRUE(pr$pr_lo > 1)
          pr_a[[length(pr_a) + 1L]] <- dplyr::bind_cols(
            key,
            tibble::tibble(event = ev, analog_set = aset),
            pr
          )
        }
      }
    }
    # Scope rows from the summed series.
    for (k in seq_len(nrow(sc))) {
      scol <- .scope_cols(dr, sc, k)
      for (metric in names(VC)) {
        v <- scope_pair_vectors(
          mats[[metric]][[aset]],
          ids[[aset]],
          units_g,
          dr,
          sc,
          k,
          afr_ids[[aset]]
        )
        Wk <- if (is.na(v$group)) W_afr[[aset]] else W[[aset]][[v$group]]
        af_s[[length(af_s) + 1L]] <- dplyr::bind_cols(
          scol,
          tibble::tibble(metric = metric, analog_set = aset),
          af_stats(
            v$f,
            v$cf,
            Wk,
            if (is.null(dr$pop_rate0)) {
              scol$pop_universe
            } else {
              sum(dr$pop_rate0[sc$idx[[k]]])
            }
          )
        )
        sn_s[[length(sn_s) + 1L]] <- dplyr::bind_cols(
          scol,
          tibble::tibble(metric = metric, analog_set = aset),
          sn_stats(v$f, v$cf)
        )
        for (m in ipcc_mags) {
          ip_s[[length(ip_s) + 1L]] <- dplyr::bind_cols(
            scol,
            tibble::tibble(metric = metric, analog_set = aset, magnitude = m),
            ipcc_scope_stats(v$f, v$cf, m, robust_exceed, robust_agree)
          )
        }
      }
      for (ev in names(EVENTS)) {
        v <- scope_pair_vectors(
          mats[[EVENTS[[ev]]$metric]][[aset]],
          ids[[aset]],
          units_g,
          dr,
          sc,
          k,
          afr_ids[[aset]]
        )
        Wk <- if (is.na(v$group)) W_afr[[aset]] else W[[aset]][[v$group]]
        pr_s[[length(pr_s) + 1L]] <- dplyr::bind_cols(
          scol,
          tibble::tibble(event = ev, analog_set = aset),
          pr_stats(
            as.numeric(EVENTS[[ev]]$fn(v$f)),
            as.numeric(EVENTS[[ev]]$fn(v$cf)),
            Wk,
            pr_alpha
          )
        )
      }
    }
  }

  # Population of the ADM1s whose PR exceeds 1 (point / bootstrap lower end).
  pr_pop <- dplyr::bind_rows(lapply(names(EVENTS), function(ev) {
    dplyr::bind_rows(lapply(c("point", "lower"), function(rule) {
      dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
        .point_scope_row(
          pr_hit[[ev]][[rule]],
          dr,
          sc,
          k,
          cond,
          event = ev,
          rule = rule
        )
      }))
    }))
  })) %>%
    add_pop_shares()
  # Population of the ADM1s where the signal has emerged at level L.
  DIRS <- list(
    increase = function(x, L) x >= L,
    decrease = function(x, L) x <= -L,
    either = function(x, L) abs(x) >= L
  )
  sn_pop <- dplyr::bind_rows(lapply(names(VC), function(metric) {
    dplyr::bind_rows(lapply(sn_levels, function(L) {
      dplyr::bind_rows(lapply(names(DIRS), function(d) {
        hit <- DIRS[[d]](sn_hit[[metric]], L)
        hit[is.na(hit)] <- FALSE
        dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
          .point_scope_row(
            hit,
            dr,
            sc,
            k,
            cond,
            metric = metric,
            sn_level = L,
            direction = d
          )
        }))
      }))
    }))
  })) %>%
    add_pop_shares()

  out <- list(
    af_adm1 = dplyr::bind_rows(af_a),
    af_scope = dplyr::bind_rows(af_s),
    pr_adm1 = dplyr::bind_rows(pr_a),
    pr_scope = dplyr::bind_rows(pr_s),
    pr_pop = pr_pop,
    sn_adm1 = dplyr::bind_rows(sn_a),
    sn_scope = dplyr::bind_rows(sn_s),
    sn_pop = sn_pop,
    universe = dr$universe,
    explosive = expl
  )
  if (!is.null(ipcc_mags)) {
    out$ipcc_scope <- dplyr::bind_rows(ip_s)
  }
  out
}

#' attrib_metrics_tables() over the attribution cells (as
#' attrib_exposure_windows()); rows tagged `window`, `explosive_guard`.
#' @param cells attrib_cell_table()
attrib_metrics_windows <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  cells,
  pr_case_threshold,
  pr_alpha,
  sn_levels,
  n_boot = ATTRIB_BOOT_B,
  boot_seed = ATTRIB_BOOT_SEED,
  ipcc_mags = NULL,
  robust_exceed = NA_real_,
  robust_agree = NA_real_
) {
  TABLES <- c(
    "af_adm1",
    "af_scope",
    "pr_adm1",
    "pr_scope",
    "pr_pop",
    "sn_adm1",
    "sn_scope",
    "sn_pop",
    "explosive",
    if (!is.null(ipcc_mags)) "ipcc_scope"
  )
  res <- list()
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    o <- attrib_metrics_tables(
      .cell_frame(adm1_month, cl),
      .cell_frame(country_month, cl),
      adm1_pop,
      group_short,
      sets,
      cl$start,
      pr_case_threshold,
      pr_alpha,
      sn_levels,
      n_boot = n_boot,
      boot_seed = boot_seed,
      cap = .cell_cap(cl),
      ipcc_mags = ipcc_mags,
      robust_exceed = robust_exceed,
      robust_agree = robust_agree
    )
    tag <- function(d) .cell_tag(d, cl)
    res[[length(res) + 1L]] <- c(
      lapply(o[TABLES], tag),
      list(universe = o$universe)
    )
  }
  out <- lapply(stats::setNames(TABLES, TABLES), function(t) {
    dplyr::bind_rows(lapply(res, `[[`, t))
  })
  out$explosive <- dplyr::filter(out$explosive, explosive_guard)
  c(out, list(universe = res[[1L]]$universe, windows = .cell_starts(cells)))
}

derive_report_attrib_metrics <- function(ctx, sets) {
  inp <- attrib_arm_inputs(ctx, sets, "report-attrib-metrics")
  win <- inp$win
  attrib_start <- inp$attrib_start
  message(sprintf(
    "[09] report_attrib_metrics: groups=%s | sets=%s | attrib_start=%s | pr_case_threshold=%g alpha=%g | sn_levels=%s | B=%d",
    paste(inp$groups, collapse = ","),
    paste(sets, collapse = ","),
    attrib_start,
    ctx$PR_CASE_THRESHOLD,
    ctx$PR_ALPHA,
    paste(ctx$SN_LEVELS, collapse = ","),
    .ctx_boot_b(ctx)
  ))
  post_fork_start <- display_factual_start(ctx, inp$groups, sets)
  cells <- attrib_cells(ctx, attrib_start, post_fork_start)
  af_posterior_stash_put(
    ctx,
    paste("forecast", paste(sort(sets), collapse = "-")),
    inp,
    cells
  )
  out <- attrib_metrics_windows(
    inp$adm1_month,
    inp$country_month,
    inp$adm1_pop,
    ctx$GROUP_SHORT_RUN,
    sets,
    cells,
    ctx$PR_CASE_THRESHOLD,
    ctx$PR_ALPHA,
    ctx$SN_LEVELS,
    n_boot = .ctx_boot_b(ctx),
    ipcc_mags = ctx$ATTRIB_MAGNITUDES,
    robust_exceed = ctx$ROBUST_EXCEED,
    robust_agree = ctx$ROBUST_AGREE
  )
  # Table S29: what the explosive guard drops, per window x set x scope, for
  # this factual (forecast = the displayed one, attrib = the attribution factual).
  out$explosive_summary <- explosive_summary_cells(
    inp$country_month,
    cells,
    sets,
    ctx$ATTRIB_FACTUAL %||% "forecast"
  )
  c(
    out,
    list(
      pr_case_threshold = ctx$PR_CASE_THRESHOLD,
      pr_alpha = ctx$PR_ALPHA,
      sn_levels = ctx$SN_LEVELS,
      magnitudes = ctx$ATTRIB_MAGNITUDES,
      robust_exceed = ctx$ROBUST_EXCEED,
      robust_agree = ctx$ROBUST_AGREE,
      n_boot = .ctx_boot_b(ctx),
      boot_seed = ATTRIB_BOOT_SEED,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      attrib_start = attrib_start,
      factual_start = win$factual_start,
      arm_starts = win$arm_starts,
      africa_min_horizon = as.Date(inp$man$africa_min_horizon),
      pop_ref_year = ctx$POP_REF_YEAR,
      analog_sets = sets,
      groups = inp$groups,
      post_fork_start = post_fork_start,
      attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
      cells = cells_payload(cells),
      africa_window = africa_window_payload(
        .africa_win(ctx, inp$groups),
        attrib_start
      )
    ),
    .sim_pop_year_fields(ctx, list(inp$man))
  ) # report_attrib_metrics v3
}

ensure_report_attrib_metrics <- function(ctx, sets) {
  ctx <- results_attrib_ctx(ctx)
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- results_run_configs(man, groups, sets)
  inputs <- results_inputs(
    files = c(
      results_attrib_manifest_files(ctx),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = if (
      results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched
    ) {
      smoothing_parts_glob(ctx)
    } else {
      list()
    }
  )
  ensure_artefact(
    ctx,
    "report_attrib_metrics",
    inputs,
    function() derive_report_attrib_metrics(ctx, sets),
    sets = sets
  )
}

# ===========================================================================
# report_event_attrib: attribution of PAST strong ENSO events.
# For each event of the event manifest (utils/make_event_attribution_configs.R,
# deploy --stage gen_events), the paired factual / counterfactual runs of every
# simulation unit, both handed off from the smoothing at the event launch:
# factual = observed ENSO, counterfactual = Lenssen analogs of one set (default
# neutral: the event did not develop). The report_attrib_metrics families (AF /
# AN, PR + FAR / PS / PNS, S/N) and the report_attrib_exposure classes (Childs,
# IPCC robustness) are derived by the same pure cores, per event x factual
# source x explosive guard. Two factual sources, both headline candidates:
#   simulation  the free factual run (observed ENSO);
#   smoothing   the smoothing replicate each factual replicate was handed off
#               from (data-conditioned; only the smoothing's distinct replicates).
# Plus the observed totals of the window (survival-panel onsets, monthly
# cases) and the observed-anchored attributable number AN_obs = AF x observed.
# Off unless results.event_manifest is set (top level, unstamped).
# ===========================================================================

#' A copy of the results context with fields replaced (an event's manifest view).
.results_ctx_with <- function(ctx, ...) {
  e <- list2env(as.list(ctx, all.names = TRUE), parent = emptyenv())
  vals <- list(...)
  for (k in names(vals)) {
    assign(k, vals[[k]], envir = e)
  }
  e
}

#' Observed totals of one event per ADM1 over [launch, the country's last
#' simulated month]: onsets (survival-panel onset months) and cases (monthly
#' case counts; reported-but-unquantified months count 0), with the gid-months
#' the panels cover (n_cells_onset / n_cells_cases) out of every simulated
#' gid-month (n_cells).
#' @param obs          event_observed_panel(surv_dat, mon_dat): gid, t_start, obs_Y, cases_month
#' @param gids         tibble(gid, country) of the simulated gids
#' @param country_last named Date (ISO3 -> last simulated month)
#' @param country_first NULL (every country from `launch`) or named Date (ISO3
#'   -> first month: an event's per-winter sub-window)
event_observed_adm1 <- function(
  obs,
  gids,
  country_last,
  launch,
  country_first = NULL
) {
  launch <- as.Date(launch)
  g <- dplyr::distinct(gids, gid, country) %>%
    dplyr::mutate(
      gid = as.character(gid),
      adm1_unit = gid_to_adm1(gid, "event_observed_adm1"),
      last = as.Date(unname(country_last[country]))
    )
  if (anyNA(g$last)) {
    stop(
      "[09] event_observed_adm1: no simulated month for country(ies) ",
      paste(unique(g$country[is.na(g$last)]), collapse = ", "),
      call. = FALSE
    )
  }
  g$first <- if (is.null(country_first)) {
    rep(launch, nrow(g))
  } else {
    as.Date(unname(country_first[g$country]))
  }
  if (anyNA(g$first)) {
    stop(
      "[09] event_observed_adm1: no first month for country(ies) ",
      paste(unique(g$country[is.na(g$first)]), collapse = ", "),
      call. = FALSE
    )
  }
  g$n_months <- vapply(
    seq_len(nrow(g)),
    function(i) {
      if (g$last[i] < g$first[i]) {
        0L
      } else {
        length(seq(g$first[i], g$last[i], by = "month"))
      }
    },
    integer(1)
  )
  cells <- g %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(n_cells = sum(n_months), .groups = "drop")
  ob <- obs %>%
    dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::inner_join(
      dplyr::select(g, gid, adm1_unit, first, last),
      by = "gid"
    ) %>%
    dplyr::filter(t_start >= first, t_start <= last) %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(
      obs_onsets = sum(obs_Y, na.rm = TRUE),
      obs_cases = sum(cases_month, na.rm = TRUE),
      n_cells_onset = sum(!is.na(obs_Y)),
      n_cells_cases = sum(!is.na(cases_month)),
      .groups = "drop"
    )
  dplyr::left_join(cells, ob, by = "adm1_unit") %>%
    dplyr::mutate(dplyr::across(
      c(obs_onsets, obs_cases, n_cells_onset, n_cells_cases),
      ~ dplyr::coalesce(as.numeric(.x), 0)
    ))
}

#' event_observed_panel
#'
#' Observed (gid, t_start) panel of the event attribution: obs_Y from the
#' survival panel (build_obs_panel()) and cases_month from the monthly case
#' frame, each gid-month's footprint total = the largest per-admin-level sum
#' (footprint.R::finest_level_totals(), the base-only twin of the filter's
#' .aggregate_footprint_counts(), which lives in spatpomp_core.R and needs
#' spatPomp, absent from 09's R stack). Reported-but-unquantified months
#' (cases_gap) stay NA.
event_observed_panel <- function(surv_dat, mon_dat) {
  obs <- build_obs_panel(surv_dat)
  m <- mon_dat %>%
    dplyr::mutate(
      gid = as.character(gid),
      t_start = as.Date(month),
      cases_month = as.numeric(cases_month)
    )
  if ("cases_gap" %in% names(m)) {
    m <- dplyr::filter(m, !(cases_gap %in% TRUE))
  }
  m <- dplyr::filter(m, !is.na(cases_month))
  cm <- if ("admin_level" %in% names(m)) {
    finest_level_totals(m, "cases_month", keys = c("gid", "t_start")) %>%
      dplyr::transmute(gid, t_start, cases_month = total)
  } else {
    m %>%
      dplyr::group_by(gid, t_start) %>%
      dplyr::summarise(cases_month = sum(cases_month), .groups = "drop")
  }
  dplyr::full_join(obs, cm, by = c("gid", "t_start")) %>%
    dplyr::arrange(gid, t_start)
}

#' Observed totals per scope (Africa, region, country) of the ADM1 universe.
event_observed_scopes <- function(obs_adm1, universe) {
  sc <- enso_scopes(universe)
  o <- obs_adm1[match(universe$adm1_gid, obs_adm1$adm1_unit), , drop = FALSE]
  dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
    idx <- sc$idx[[k]]
    n <- sum(o$n_cells[idx], na.rm = TRUE)
    tibble::tibble(
      level = sc$level[k],
      scope = sc$scope[k],
      enso_group = sc$enso_group[k],
      region = sc$region[k],
      country = sc$country[k],
      obs_onsets = sum(o$obs_onsets[idx], na.rm = TRUE),
      obs_cases = sum(o$obs_cases[idx], na.rm = TRUE),
      cover_onsets = if (n > 0) {
        sum(o$n_cells_onset[idx], na.rm = TRUE) / n
      } else {
        NA_real_
      },
      cover_cases = if (n > 0) {
        sum(o$n_cells_cases[idx], na.rm = TRUE) / n
      } else {
        NA_real_
      }
    )
  }))
}

#' Observed-anchored attributable number: the scope's AF (of the replicate
#' means, with its bootstrap interval) times the observed window total of the
#' same metric; the observed total is held fixed. Plus the factual mean against
#' the observed total (obs_over_f), the skill check of the factual arm.
#' `cover` is the surveillance coverage of the simulated gid-months for both
#' metrics: the share with a survival-panel row (cover_onsets). Monthly case
#' rows exist only in episode months, so cover_cases is not a coverage.
event_an_observed <- function(af_scope, obs_scope) {
  ob <- dplyr::bind_rows(
    dplyr::transmute(
      obs_scope,
      level,
      scope,
      enso_group,
      metric = "cases",
      observed = obs_cases,
      cover = cover_onsets
    ),
    dplyr::transmute(
      obs_scope,
      level,
      scope,
      enso_group,
      metric = "onsets",
      observed = obs_onsets,
      cover = cover_onsets
    )
  )
  af_scope %>%
    dplyr::left_join(ob, by = c("level", "scope", "enso_group", "metric")) %>%
    dplyr::mutate(
      an_obs = af_of_means * observed,
      an_obs_lo = af_lo * observed,
      an_obs_hi = af_hi * observed,
      obs_over_f = ifelse(mean_f > 0, observed / mean_f, NA_real_)
    )
}

#' Per-event cells: attrib_metrics_tables() + attrib_exposure_tables() for both
#' explosive guards, rows tagged with `explosive_guard`. Both arms launch at the
#' event launch, so the window is the whole run (no stitched prefix).
#' `cells`: one sub-window's cells (event_sub_window_cells());
#' the frames are then cut to it (.cell_frame()) and its start is the window
#' start. NULL (default) = the whole run, as before.
event_attrib_cells <- function(inp, ctx, sets, cells = NULL) {
  launch <- inp$attrib_start
  sub <- !is.null(cells)
  ec <- if (sub) cells else event_cells(ctx, launch, "event") # one window (the whole run) x the explosive modes
  cells <- list()
  for (i in seq_len(nrow(ec))) {
    cl <- ec[i, ]
    am <- if (sub) .cell_frame(inp$adm1_month, cl) else inp$adm1_month
    cm <- if (sub) .cell_frame(inp$country_month, cl) else inp$country_month
    st <- if (sub) cl$start[[1]] else launch
    m <- attrib_metrics_tables(
      am,
      cm,
      inp$adm1_pop,
      ctx$GROUP_SHORT_RUN,
      sets,
      st,
      ctx$PR_CASE_THRESHOLD,
      ctx$PR_ALPHA,
      ctx$SN_LEVELS,
      cap = .cell_cap(cl),
      n_boot = .ctx_boot_b(ctx)
    )
    x <- attrib_exposure_tables(
      am,
      cm,
      inp$adm1_pop,
      ctx$GROUP_SHORT_RUN,
      sets,
      st,
      ctx$ATTRIB_MAGNITUDES,
      ctx$ROBUST_EXCEED,
      ctx$ROBUST_AGREE,
      cap = .cell_cap(cl),
      n_boot = .ctx_boot_b(ctx)
    )
    tag <- function(d) .cell_tag(d, cl, window = FALSE)
    cells[[as.character(cl$explosive_guard)]] <- c(
      lapply(
        m[c(
          "af_adm1",
          "af_scope",
          "pr_adm1",
          "pr_scope",
          "pr_pop",
          "sn_adm1",
          "sn_scope",
          "sn_pop",
          "explosive"
        )],
        tag
      ),
      list(
        exposure = tag(x$exposure),
        exposure_adm1 = tag(x$adm1),
        universe = m$universe
      )
    )
  }
  tabs <- setdiff(names(cells[[1L]]), "universe")
  out <- lapply(stats::setNames(tabs, tabs), function(t) {
    dplyr::bind_rows(lapply(cells, `[[`, t))
  })
  out$explosive <- dplyr::filter(out$explosive, explosive_guard)
  c(out, list(universe = cells[[1L]]$universe))
}

EVENT_FACTUAL_SOURCES <- c("simulation", "smoothing")

derive_report_event_attrib <- function(ctx) {
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  want_parts <- sub("/+$", "", smoothing_parts_dir_of(em$smoothing_input))
  if (!identical(want_parts, sub("/+$", "", ctx$SMOOTHING_PARTS_DIR))) {
    stop(
      sprintf(
        paste0(
          "[09] report_event_attrib: the event runs hand off from %s but ",
          "results.inputs.smoothing_parts_dir is %s; the smoothing-as-factual reading ",
          "would pair them with another smoothing run."
        ),
        want_parts,
        ctx$SMOOTHING_PARTS_DIR
      ),
      call. = FALSE
    )
  }
  sets <- as.character(em$analog_set)
  events <- names(em$events)
  message(sprintf(
    "[09] report_event_attrib: %d event(s) %s | set %s | factual sources %s | B=%d",
    length(events),
    paste(events, collapse = ","),
    sets,
    paste(EVENT_FACTUAL_SOURCES, collapse = "/"),
    .ctx_boot_b(ctx)
  ))
  base <- attrib_base_inputs(ctx)
  mon_dat <- readRDS(require_input(ctx$MON_DATA_PATH, "monthly case data"))
  obs <- event_observed_panel(base$surv_dat, mon_dat)
  first_launch <- min(as.Date(vapply(em$events, `[[`, "", "launch")))
  sm_ev <- load_smoothing_events(ctx, ctx$SMOOTHING_PARTS_DIR, first_launch)

  acc <- list()
  meta <- list()
  obs_rows <- list()
  an_rows <- list()
  # Per-winter sub-windows (schema-3 manifests): one cell per
  # sub-window beside the whole event, in `sub_windows` (absent otherwise).
  acc_sw <- list()
  obs_sw <- list()
  an_sw <- list()
  win_sw <- list()
  for (ev in events) {
    view <- event_manifest_view(em, ev)
    ctx_e <- .results_ctx_with(
      ctx,
      MAN = view,
      MANIFEST_PATH = ctx$EVENT_MANIFEST_PATH,
      AFRICA_WIN = NULL
    ) # the forecast Africa window never applies to an event
    e <- em$events[[ev]]
    swt <- event_sub_windows_table(em, ev)
    if (!is.null(swt)) {
      win_sw[[ev]] <- dplyr::mutate(swt, event = ev, .before = 1L)
    }
    obs_ev_sw <- NULL
    # Library-mode event caches: the runs drew their ENSO from
    # the recorded cache. Legacy events carry no lenssen_digest: nothing checked.
    if (!is.null(e$lenssen_digest)) {
      assert_counterfactual_cache_identity(
        list(cache_path = e$lenssen_cache, digest = e$lenssen_digest),
        manifest_runs(view, "counterfactual", sets),
        label = sprintf("[09] report_event_attrib %s:", ev)
      )
    }
    obs_ev <- NULL
    for (fs in EVENT_FACTUAL_SOURCES) {
      message(sprintf(
        "[09]   event %s (launch %s) factual = %s",
        ev,
        e$launch,
        fs
      ))
      inp <- attrib_arm_inputs(
        ctx_e,
        sets,
        "report-event-attrib",
        factual_source = fs,
        smoothing_events = if (identical(fs, "smoothing")) sm_ev else NULL,
        base = base
      )
      if (identical(fs, "simulation")) {
        af_posterior_stash_put(
          ctx,
          paste("event", ev, paste(sort(sets), collapse = "-")),
          inp,
          event_cells_all(ctx, inp$attrib_start, ev, swt)
        )
      }
      cl <- event_attrib_cells(inp, ctx, sets)
      # Table S29 for the event: what the explosive guard drops (window = the event).
      cl$explosive_summary <- explosive_summary_cells(
        inp$country_month,
        event_cells(ctx, inp$attrib_start, ev),
        sets,
        fs
      )
      if (is.null(obs_ev)) {
        cm_f <- dplyr::filter(inp$country_month, run == "forecast")
        country_last <- tapply(cm_f$month, cm_f$country, max)
        country_last <- stats::setNames(
          as.Date(country_last, origin = "1970-01-01"),
          names(country_last)
        )
        obs_adm1 <- event_observed_adm1(obs, inp$gids, country_last, e$launch)
        obs_ev <- list(
          adm1 = obs_adm1,
          scope = event_observed_scopes(obs_adm1, cl$universe)
        )
        obs_rows[[ev]] <- list(
          adm1 = dplyr::mutate(obs_adm1, event = ev, .before = 1L),
          scope = dplyr::mutate(obs_ev$scope, event = ev, .before = 1L)
        )
      }
      an_rows[[paste(ev, fs)]] <- dplyr::mutate(
        event_an_observed(cl$af_scope, obs_ev$scope),
        event = ev,
        factual_source = fs,
        .before = 1L
      )
      for (t in setdiff(names(cl), "universe")) {
        acc[[t]][[paste(ev, fs)]] <- .event_attrib_stamp(
          cl[[t]],
          event = ev,
          factual_source = fs
        )
      }
      if (identical(fs, EVENT_FACTUAL_SOURCES[1L])) {
        acc$universe[[ev]] <- dplyr::mutate(
          cl$universe,
          event = ev,
          .before = 1L
        )
      }
      if (!is.null(swt)) {
        # Observed totals per sub-window: each country over its own sub-window
        # bounds, cut at its last simulated month.
        if (is.null(obs_ev_sw)) {
          obs_ev_sw <- event_sub_window_observed(
            obs,
            inp$gids,
            country_last,
            swt,
            view,
            cl$universe,
            ev
          )
        }
        for (k in sort(unique(swt$k))) {
          r <- swt[swt$k == k, , drop = FALSE]
          id <- r$sub_window[1L]
          sc <- event_sub_window_cells(ctx, inp$attrib_start, r)
          cl_k <- event_attrib_cells(inp, ctx, sets, sc)
          cl_k$explosive_summary <- explosive_summary_cells(
            inp$country_month,
            sc,
            sets,
            fs
          )
          tag_k <- function(d) {
            .event_attrib_stamp(
              d,
              event = ev,
              sub_window = id,
              factual_source = fs
            )
          }
          for (t in setdiff(names(cl_k), "universe")) {
            acc_sw[[t]][[paste(ev, fs, id)]] <- tag_k(cl_k[[t]])
          }
          an_sw[[paste(ev, fs, id)]] <- tag_k(event_an_observed(
            cl_k$af_scope,
            obs_ev_sw[[id]]$scope
          ))
        }
        if (identical(fs, EVENT_FACTUAL_SOURCES[1L])) {
          obs_sw[[ev]] <- obs_ev_sw
        }
      }
    }
    meta[[ev]] <- tibble::tibble(
      event = ev,
      phase = e$phase,
      label = e$label,
      launch = as.Date(e$launch),
      peak = as.Date(e$peak),
      peak_season_year = as.integer(e$peak_season_year),
      window_end = as.Date(e$africa_max_horizon),
      lenssen_n_set = as.integer(e$lenssen_n_set %||% NA)
    )
  }
  out <- lapply(acc, dplyr::bind_rows)
  res <- c(
    out,
    list(
      observed_adm1 = dplyr::bind_rows(lapply(obs_rows, `[[`, "adm1")),
      observed_scope = dplyr::bind_rows(lapply(obs_rows, `[[`, "scope")),
      an_observed = dplyr::bind_rows(an_rows),
      events = dplyr::bind_rows(meta),
      group_windows = dplyr::bind_rows(lapply(events, function(ev) {
        g <- em$events[[ev]]$groups
        tibble::tibble(
          event = ev,
          enso_group = names(g),
          season_end = as.Date(vapply(g, `[[`, "", "season_end")),
          window_end = as.Date(vapply(g, `[[`, "", "horizon_end")),
          max_lag = vapply(g, function(x) as.integer(x$max_lag), integer(1))
        )
      })),
      # Per-country windows: each country's season end
      # and firm window end (its own under per-country season ends, else
      # its group's); group_windows$window_end is the latest member's.
      country_windows = dplyr::select(
        season_windows_table(em, "event"),
        -window_type
      ) %>%
        dplyr::rename(event = window),
      analog_set = sets,
      factual_sources = EVENT_FACTUAL_SOURCES,
      pr_case_threshold = ctx$PR_CASE_THRESHOLD,
      pr_alpha = ctx$PR_ALPHA,
      sn_levels = ctx$SN_LEVELS,
      magnitudes = ctx$ATTRIB_MAGNITUDES,
      robust_exceed = ctx$ROBUST_EXCEED,
      robust_agree = ctx$ROBUST_AGREE,
      n_boot = .ctx_boot_b(ctx),
      boot_seed = ATTRIB_BOOT_SEED,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      pop_ref_year = ctx$POP_REF_YEAR,
      event_manifest = ctx$EVENT_MANIFEST_PATH
    ),
    .sim_pop_year_fields(ctx, list(em))
  )
  # v3: per-winter sub-windows of schema-3 events, the whole
  # event's tables unchanged above. Only when an event has sub-windows, so a
  # schema-2 manifest derives the v2 payload exactly.
  if (length(win_sw)) {
    res$sub_windows <- c(
      lapply(acc_sw, dplyr::bind_rows),
      list(
        observed_adm1 = dplyr::bind_rows(lapply(names(obs_sw), function(ev) {
          dplyr::bind_rows(lapply(names(obs_sw[[ev]]), function(id) {
            dplyr::mutate(
              obs_sw[[ev]][[id]]$adm1,
              event = ev,
              sub_window = id,
              .before = 1L
            )
          }))
        })),
        observed_scope = dplyr::bind_rows(lapply(names(obs_sw), function(ev) {
          dplyr::bind_rows(lapply(names(obs_sw[[ev]]), function(id) {
            dplyr::mutate(
              obs_sw[[ev]][[id]]$scope,
              event = ev,
              sub_window = id,
              .before = 1L
            )
          }))
        })),
        an_observed = dplyr::bind_rows(an_sw),
        windows = dplyr::bind_rows(win_sw)
      )
    )
  }
  res # report_event_attrib v4
}

#' Rows (country, month) inside one sub-window's bounds (`rk`: its rows of
#' event_sub_windows_table()): per country at the country level, else the
#' bounds of `group` (the fit group of every row).
event_sub_window_keep <- function(country, month, rk, group = NA_character_) {
  key <- if (identical(rk$level[1L], "country")) {
    as.character(country)
  } else {
    rep(as.character(group), length(country))
  }
  i <- match(key, rk$key)
  !is.na(i) & as.Date(month) >= rk$start[i] & as.Date(month) <= rk$end[i]
}

#' Observed totals of each per-winter sub-window of one event:
#' event_observed_adm1() with every country from its own sub-window start
#' to min(its sub-window end, its last simulated month), and the scope sums on
#' the whole event's ADM1 universe. Group-level bounds (a by_group manifest)
#' apply to the group's countries (`view$countries`).
#' @return named list (sub-window id) of list(adm1, scope)
event_sub_window_observed <- function(
  obs,
  gids,
  country_last,
  swt,
  view,
  universe,
  event = ""
) {
  isos <- names(country_last)
  grp_of <- vapply(
    isos,
    function(i) {
      as.character((view$countries[[i]] %||% list())$group %||% NA_character_)
    },
    ""
  )
  out <- list()
  for (k in sort(unique(swt$k))) {
    r <- swt[swt$k == k, , drop = FALSE]
    key <- if (identical(r$level[1L], "country")) isos else grp_of
    i <- match(key, r$key)
    if (anyNA(i)) {
      stop(
        sprintf(
          "[09] report_event_attrib %s / %s: no sub-window bounds for country(ies) %s.",
          event,
          r$sub_window[1L],
          paste(isos[is.na(i)], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    first <- stats::setNames(r$start[i], isos)
    last <- stats::setNames(
      as.Date(
        pmin(as.numeric(country_last[isos]), as.numeric(r$end[i])),
        origin = "1970-01-01"
      ),
      isos
    )
    adm1 <- event_observed_adm1(
      obs,
      gids,
      last,
      min(r$start),
      country_first = first
    )
    out[[r$sub_window[1L]]] <- list(
      adm1 = adm1,
      scope = event_observed_scopes(adm1, universe)
    )
  }
  out
}

ensure_report_event_attrib <- function(ctx) {
  if (!nzchar(ctx$EVENT_MANIFEST_PATH %||% "")) {
    message(
      "[09] report_event_attrib: results.event_manifest is empty; skipped."
    )
    return(invisible(NULL))
  }
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  cfgs <- unlist(lapply(em$events, function(e) {
    unlist(lapply(e$units, function(u) {
      c(u$factual_config, vapply(u$counterfactual_sets, `[[`, "", "config"))
    }))
  }))
  inputs <- results_inputs(
    files = c(
      ctx$EVENT_MANIFEST_PATH,
      cfgs,
      ctx$SURV_PATH,
      ctx$MON_DATA_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = smoothing_parts_glob(ctx)
  )
  ensure_artefact(ctx, "report_event_attrib", inputs, function() {
    derive_report_event_attrib(ctx)
  })
}

# ===========================================================================
# report_event_af_model: ANALYTIC (crossbasis-only) attributable
# fraction of the same past events, a cross-check of report_event_attrib with
# no epidemic dynamics. On the OBSERVED rows of each event window (occurrence:
# survival-panel onset months; monthly: observed monthly cases), per posterior
# draw d of the deployed fit and a counterfactual ENSO trajectory drawn for it
# from the event's Lenssen analog set:
#   delta_eta(t, d) = [Z(ONI_obs(t - L)) - Z(ONI_cf,d(t - L))] . b_d(ADM1)
#   AN_d = sum_t y_t (1 - exp(-delta_eta(t, d))),  AF_d = AN_d / sum_t y_t
# (backward attributable number of a distributed-lag model, Gasparrini & Leone
# 2014). b_d = population + region + country + ADM1 slope draws (country-level
# coefficients for an ADM1 without a random effect). Z is the fit's ENSO basis
# at its single lag L, clamped as the forward runs clamp it. Months before the
# launch keep the observed ONI in both arms, as in the simulations.
# ===========================================================================

#' Month index (year * 12 + month - 1) of first-of-month Dates.
.month_index <- function(d) {
  as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m")) - 1L
}
.month_from_index <- function(i) {
  as.Date(sprintf("%04d-%02d-01", i %/% 12L, i %% 12L + 1L))
}

#' event_model_af_rows
#'
#' Per-draw attributable outcome of observed rows under one fit (pure; no brms).
#' @param rows    tibble(adm1_gid, month (Date), y): observed rows with y > 0
#' @param coef_of function(adm1_gid) -> D x K coefficient draws of that ADM1
#' @param vm      the fit's cb_meta ENSO entry (single_lag basis: onebasis_attrs,
#'   cb_col_names, range)
#' @param lag     the fit's ENSO lag (exposure month = month - lag; negative = lead)
#' @param oni     named numeric: "YYYY-MM-01" -> observed ONI
#' @param cf      D x M counterfactual ONI from the launch (column 1 = launch)
#' @param launch  first counterfactual month
#' @param clamp_margin NULL (no clamp) or the forward runs' clamp margin
#' @return list(units = ADM1 gids, observed = per-unit sum of y,
#'   an = units x D matrix of attributable y)
event_model_af_rows <- function(
  rows,
  coef_of,
  vm,
  lag,
  oni,
  cf,
  launch,
  clamp_margin = NULL
) {
  D <- nrow(cf)
  units <- sort(unique(rows$adm1_gid))
  if (!nrow(rows)) {
    return(list(
      units = character(0),
      observed = numeric(0),
      an = matrix(0, 0L, D)
    ))
  }
  ex <- .month_index(as.Date(rows$month)) - as.integer(lag)
  x_obs <- unname(oni[format(.month_from_index(ex))])
  if (anyNA(x_obs)) {
    stop(
      "[09] event_model_af_rows: no observed ONI for exposure month(s) ",
      paste(
        head(unique(format(.month_from_index(ex[is.na(x_obs)]))), 5L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  col <- ex - .month_index(as.Date(launch)) + 1L
  is_cf <- col >= 1L
  if (any(col[is_cf] > ncol(cf))) {
    stop(
      "[09] event_model_af_rows: counterfactual trajectories end before exposure month ",
      format(.month_from_index(max(ex))),
      call. = FALSE
    )
  }
  basis <- function(x) {
    if (!is.null(clamp_margin)) {
      x <- .clamp_cb_range(x, vm$range, clamp_margin)
    }
    build_cb_columns_from_Q(x, vm)
  }
  n <- nrow(rows)
  Z_obs <- basis(x_obs)
  X_cf <- matrix(x_obs, n, D)
  if (any(is_cf)) {
    X_cf[is_cf, ] <- t(cf[, col[is_cf], drop = FALSE])
  }
  Z_cf <- basis(as.vector(X_cf)) # (n * D) x K, column-major in d
  an <- matrix(0, length(units), D, dimnames = list(units, NULL))
  obs <- stats::setNames(numeric(length(units)), units)
  for (u in units) {
    r <- which(rows$adm1_gid == u)
    nr <- length(r)
    B <- coef_of(u) # D x K
    eta_obs <- Z_obs[r, , drop = FALSE] %*% t(B) # nr x D
    ii <- rep(r, D) + rep((seq_len(D) - 1L) * n, each = nr)
    eta_cf <- matrix(
      rowSums(
        Z_cf[ii, , drop = FALSE] * B[rep(seq_len(D), each = nr), , drop = FALSE]
      ),
      nr,
      D
    )
    an[u, ] <- colSums(rows$y[r] * (1 - exp(-(eta_obs - eta_cf))))
    obs[[u]] <- sum(rows$y[r])
  }
  list(units = units, observed = obs, an = an)
}

#' Scope summaries of event_model_af_rows() results bound over groups: AN and
#' AF = AN / observed per draw, mean and 2.5 / 97.5% quantiles, P(AN > 0).
#' @param universe tibble(adm1_gid, country, enso_group, region)
#' @param observed named numeric over universe$adm1_gid (0 where no y)
#' @param an       universe x D matrix (rows in universe order, 0 where no y)
event_model_af_scopes <- function(universe, observed, an) {
  sc <- enso_scopes(universe)
  summ <- function(v) {
    c(
      mean = mean(v),
      lo = stats::quantile(v, .025, names = FALSE),
      hi = stats::quantile(v, .975, names = FALSE)
    )
  }
  one <- function(lev, scope, eg, rg, ctry, idx) {
    o <- sum(observed[idx])
    a <- colSums(an[idx, , drop = FALSE])
    sa <- summ(a)
    sf <- if (o > 0) {
      summ(a / o)
    } else {
      c(mean = NA_real_, lo = NA_real_, hi = NA_real_)
    }
    tibble::tibble(
      level = lev,
      scope = scope,
      enso_group = eg,
      region = rg,
      country = ctry,
      n_adm1 = length(idx),
      observed = o,
      an_mean = sa[["mean"]],
      an_lo = sa[["lo"]],
      an_hi = sa[["hi"]],
      af_mean = sf[["mean"]],
      af_lo = sf[["lo"]],
      af_hi = sf[["hi"]],
      p_an_pos = mean(a > 0)
    )
  }
  dplyr::bind_rows(
    dplyr::bind_rows(lapply(seq_len(nrow(sc)), function(k) {
      one(
        sc$level[k],
        sc$scope[k],
        sc$enso_group[k],
        sc$region[k],
        sc$country[k],
        sc$idx[[k]]
      )
    })),
    dplyr::bind_rows(lapply(seq_len(nrow(universe)), function(i) {
      one(
        "adm1",
        universe$adm1_gid[i],
        universe$enso_group[i],
        universe$region[i],
        universe$country[i],
        i
      )
    }))
  )
}

#' Per-ADM1 coefficient-draw lookup of one fit: population + region + country
#' + ADM1 slope draws, the country-level combination for an ADM1 the fit has
#' no random effect for, population draws alone for an unseen country.
event_model_coef_lookup <- function(sc) {
  adm1_lv <- dimnames(sc$re_mean[[sc$grp$adm1]] %||% array(dim = c(0, 0, 0)))[[
    1L
  ]] %||%
    character(0)
  ctry_lv <- dimnames(
    sc$re_mean[[sc$grp$country]] %||% array(dim = c(0, 0, 0))
  )[[1L]] %||%
    character(0)
  by_adm1 <- list()
  by_ctry <- list()
  for (l in ctry_lv) {
    p <- sc$parse_level(l, "country")
    by_ctry[[p$country]] <- sc$b_pop_draws +
      sc$re_draws_mat("region", p$region) +
      sc$re_draws_mat("country", l)
  }
  for (l in adm1_lv) {
    p <- sc$parse_level(l, "adm1")
    by_adm1[[p$adm1_gid]] <- (by_ctry[[p$country]] %||% sc$b_pop_draws) +
      sc$re_draws_mat("adm1", l)
  }
  function(adm1_gid) {
    by_adm1[[adm1_gid]] %||%
      by_ctry[[substr(adm1_gid, 1L, 3L)]] %||%
      sc$b_pop_draws
  }
}

#' Counterfactual ONI per draw for one event: D rows drawn (seeded, with
#' replacement) from the event cache's analog set; columns from the launch.
event_cf_trajectories <- function(cache_path, analog_set, D, seed) {
  cache <- readRDS(cache_path)
  rows <- enso_analog_set_rows(cache, analog_set)
  if (!length(rows)) {
    stop(
      "[09] event_cf_trajectories: analog set '",
      analog_set,
      "' is empty in ",
      cache_path,
      call. = FALSE
    )
  }
  pick <- .with_seed(seed, sample(rows, D, replace = TRUE))
  cache$oni_trajectories[pick, , drop = FALSE]
}

EVENT_MODEL_TRACKS <- c(occurrence = "onsets", monthly = "cases")

derive_report_event_af_model <- function(ctx) {
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  events <- names(em$events)
  groups <- names(em$events[[1L]]$groups)
  D <- as.integer(ctx$N_IRR_DRAWS)
  message(sprintf(
    "[09] report_event_af_model: %d event(s) x tracks %s x %d group(s); D = %d draws",
    length(events),
    paste(names(EVENT_MODEL_TRACKS), collapse = "/"),
    length(groups),
    D
  ))
  oni_df <- fetch_enso_index(
    "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  )
  oni <- stats::setNames(
    oni_df$enso_index,
    sprintf("%04d-%02d-01", oni_df$year, oni_df$month)
  )
  surv_dat <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  mon_dat <- readRDS(require_input(ctx$MON_DATA_PATH, "monthly case data"))
  obs <- event_observed_panel(surv_dat, mon_dat) %>%
    dplyr::mutate(gid = as.character(gid), country = gid_iso3(gid))
  rm(surv_dat, mon_dat)
  # Clamp margin of the forward runs (both arms of every event share the template).
  sim1 <- yaml::read_yaml(
    em$events[[1L]]$units[[1L]]$factual_config
  )$simulator$enso
  clamp_margin <- if (isTRUE(sim1$clamp_to_knots)) {
    as.numeric(sim1$clamp_margin %||% 0.25)
  } else {
    NULL
  }
  cf <- lapply(stats::setNames(events, events), function(ev) {
    event_cf_trajectories(
      em$events[[ev]]$lenssen_cache,
      em$analog_set,
      D,
      ATTRIB_BOOT_SEED + match(ev, events)
    )
  })
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  fit_key <- c(occurrence = "occ_fit_file", monthly = "mon_fit_file")
  # Each country's own window end: countries.<iso>.
  # firm_horizon_end when the manifest has per-country season ends, else its
  # group's horizon_end (identical for a by_group manifest).
  ctry_end <- lapply(stats::setNames(events, events), function(ev) {
    as.Date(manifest_country_horizons(event_manifest_view(em, ev), "firm"))
  })
  # Per-winter sub-windows of schema-3 events (NULL per event otherwise).
  swt_all <- lapply(stats::setNames(events, events), function(ev) {
    event_sub_windows_table(em, ev)
  })

  jobs <- expand.grid(
    track = names(EVENT_MODEL_TRACKS),
    group = groups,
    stringsAsFactors = FALSE
  )
  res <- map_workers(
    seq_len(nrow(jobs)),
    function(j) {
      tr <- jobs$track[j]
      g <- jobs$group[j]
      ff <- track_fit_files(mnf, tr, g)
      want <- em$events[[1L]]$groups[[g]][[fit_key[[tr]]]]
      if (!identical(basename(ff$fit), basename(want))) {
        stop(
          sprintf(
            "[09] report_event_af_model: %s / %s: the best-lag archive deploys %s but the event runs used %s.",
            tr,
            g,
            basename(ff$fit),
            basename(want)
          ),
          call. = FALSE
        )
      }
      message(sprintf("[09]   loading %s | %s", tr, g))
      fit <- readRDS(ff$fit)
      vm <- readRDS(ff$meta)[["enso"]]
      sc <- brms_group_re_scaffold(
        fit,
        vm$cb_col_names,
        n_draws = D,
        seed = 42L
      )
      rm(fit)
      gc(verbose = FALSE)
      coef_of <- event_model_coef_lookup(sc)
      lag <- as.integer(vm$single_lag)
      isos <- unlist(em$events[[1L]]$groups[[g]]$iso3)
      lapply(stats::setNames(events, events), function(ev) {
        e <- em$events[[ev]]
        ce <- ctry_end[[ev]]
        if (any(!isos %in% names(ce))) {
          stop(
            sprintf(
              "[09] report_event_af_model: %s / %s: no window end for %s in the event manifest.",
              ev,
              g,
              paste(setdiff(isos, names(ce)), collapse = ", ")
            ),
            call. = FALSE
          )
        }
        win <- obs %>%
          dplyr::filter(country %in% isos, t_start >= as.Date(e$launch)) %>%
          dplyr::filter(t_start <= ce[country]) %>%
          dplyr::mutate(adm1_gid = gid_to_adm1(gid, "report_event_af_model"))
        y <- if (identical(tr, "occurrence")) win$obs_Y else win$cases_month
        rows <- tibble::tibble(
          adm1_gid = win$adm1_gid,
          month = win$t_start,
          y = as.numeric(y)
        ) %>%
          dplyr::filter(!is.na(y), y > 0)
        r <- event_model_af_rows(
          rows,
          coef_of,
          vm,
          lag,
          oni,
          cf[[ev]],
          as.Date(e$launch),
          clamp_margin
        )
        # Per-winter sub-windows: the same rows cut to each
        # sub-window's bounds (per country, or the group's), same draws and
        # counterfactual trajectories (they start at the event launch).
        swt <- swt_all[[ev]]
        sw <- if (!is.null(swt)) {
          lapply(split(swt, swt$sub_window), function(rk) {
            wk <- win[
              event_sub_window_keep(win$country, win$t_start, rk, g),
              ,
              drop = FALSE
            ]
            yk <- if (identical(tr, "occurrence")) wk$obs_Y else wk$cases_month
            rows_k <- tibble::tibble(
              adm1_gid = wk$adm1_gid,
              month = wk$t_start,
              y = as.numeric(yk)
            ) %>%
              dplyr::filter(!is.na(y), y > 0)
            list(
              r = event_model_af_rows(
                rows_k,
                coef_of,
                vm,
                lag,
                oni,
                cf[[ev]],
                as.Date(e$launch),
                clamp_margin
              )
            )
          })
        }
        out <- list(
          units = unique(c(r$units, win$adm1_gid)),
          r = r,
          lag = lag,
          enso_group = g,
          track = tr
        )
        if (!is.null(sw)) {
          out$sw <- sw
        }
        out
      })
    },
    workers = ctx$WORKERS
  )

  adm1 <- list()
  scopes <- list()
  sw_adm1 <- list()
  sw_scopes <- list()
  for (ev in events) {
    for (tr in names(EVENT_MODEL_TRACKS)) {
      parts <- lapply(res[jobs$track == tr], `[[`, ev)
      universe <- dplyr::bind_rows(lapply(parts, function(p) {
        tibble::tibble(
          adm1_gid = sort(unique(p$units)),
          enso_group = p$enso_group
        )
      })) %>%
        dplyr::distinct(adm1_gid, .keep_all = TRUE) %>%
        dplyr::mutate(
          country = substr(adm1_gid, 1L, 3L),
          region = unname(ctx$GROUP_SHORT_RUN[enso_group])
        )
      observed <- stats::setNames(numeric(nrow(universe)), universe$adm1_gid)
      an <- matrix(
        0,
        nrow(universe),
        D,
        dimnames = list(universe$adm1_gid, NULL)
      )
      for (p in parts) {
        if (length(p$r$units)) {
          observed[p$r$units] <- p$r$observed
          an[p$r$units, ] <- p$r$an
        }
      }
      tab <- event_model_af_scopes(universe, observed, an) %>%
        dplyr::mutate(
          event = ev,
          track = tr,
          metric = EVENT_MODEL_TRACKS[[tr]],
          .before = 1L
        )
      scopes[[paste(ev, tr)]] <- dplyr::filter(tab, level != "adm1")
      adm1[[paste(ev, tr)]] <- dplyr::filter(tab, level == "adm1")
      # Sub-window cells on the whole event's ADM1 universe.
      for (id in unique(swt_all[[ev]]$sub_window)) {
        obs_k <- stats::setNames(numeric(nrow(universe)), universe$adm1_gid)
        an_k <- matrix(
          0,
          nrow(universe),
          D,
          dimnames = list(universe$adm1_gid, NULL)
        )
        for (p in parts) {
          q <- p$sw[[id]]$r
          if (length(q$units)) {
            obs_k[q$units] <- q$observed
            an_k[q$units, ] <- q$an
          }
        }
        tab_k <- event_model_af_scopes(universe, obs_k, an_k) %>%
          dplyr::mutate(
            event = ev,
            sub_window = id,
            track = tr,
            metric = EVENT_MODEL_TRACKS[[tr]],
            .before = 1L
          )
        sw_scopes[[paste(ev, tr, id)]] <- dplyr::filter(tab_k, level != "adm1")
        sw_adm1[[paste(ev, tr, id)]] <- dplyr::filter(tab_k, level == "adm1")
      }
    }
  }
  lags <- dplyr::bind_rows(lapply(seq_len(nrow(jobs)), function(j) {
    tibble::tibble(
      track = jobs$track[j],
      enso_group = jobs$group[j],
      lag = res[[j]][[1L]]$lag
    )
  }))
  out <- list(
    scope = dplyr::bind_rows(scopes),
    adm1 = dplyr::bind_rows(adm1),
    lags = lags,
    n_draws = D,
    analog_set = em$analog_set,
    clamp_margin = clamp_margin %||% NA_real_,
    event_manifest = ctx$EVENT_MANIFEST_PATH
  )
  # v2: per-winter sub-window cells of schema-3 events; absent
  # otherwise, so a schema-2 manifest derives the v1 payload exactly.
  if (length(sw_scopes)) {
    out$sub_windows <- list(
      scope = dplyr::bind_rows(sw_scopes),
      adm1 = dplyr::bind_rows(sw_adm1),
      windows = dplyr::bind_rows(lapply(names(swt_all), function(ev) {
        if (!is.null(swt_all[[ev]])) {
          dplyr::mutate(swt_all[[ev]], event = ev, .before = 1L)
        }
      }))
    )
  }
  out # report_event_af_model v2
}

ensure_report_event_af_model <- function(ctx) {
  if (!nzchar(ctx$EVENT_MANIFEST_PATH %||% "")) {
    message(
      "[09] report_event_af_model: results.event_manifest is empty; skipped."
    )
    return(invisible(NULL))
  }
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  fits <- unlist(lapply(c("occurrence", "monthly"), function(tr) {
    unlist(track_fit_input_files(ctx, tr))
  }))
  inputs <- results_inputs(
    files = c(
      ctx$EVENT_MANIFEST_PATH,
      fits,
      ctx$SURV_PATH,
      ctx$MON_DATA_PATH,
      vapply(em$events, `[[`, "", "lenssen_cache"),
      file.path(ctx$ENSO_CACHE_DIR, "enso_oni_cache.rds")
    )
  )
  ensure_artefact(ctx, "report_event_af_model", inputs, function() {
    derive_report_event_af_model(ctx)
  })
}

# ===========================================================================
# Supplementary tables: convergence (S8, fig. S5), lag selection
# and rankings (S9, S10, fig. S6), baseline-rate anchoring (S7), software (S13)
# and the brms_kl class-sweep comparison. The dataset tables S1 / S2 extend
# report_dataset above.
# ===========================================================================

# --- MCMC convergence of the deployed fits -----------------------------------

# --- Fit-level model specification -----------
# report_model_spec (track-keyed, table S34): what each DEPLOYED fit actually
# contains -- formula terms (fixed / random / offsets / dispersion), the RE
# grouping factors, family and link, the fixed-effect coefficients and the
# priors -- read from the fit objects in the convergence step's own loop (no
# extra fit load), and checked against analysis/model_terms_deployed.yml (the
# config-level statement gated by tests/test_deployed_model_terms.R).

MODEL_TERMS_DEPLOYED <- "analysis/model_terms_deployed.yml"

#' Normalised term string: no whitespace, single quotes.
.norm_term <- function(x) gsub("\"", "'", gsub("[[:space:]]+", "", x))

#' RE grouping factors of RE terms ("(a || g1/g2)" -> g1, g1:g2).
.re_groups_of <- function(re_terms) {
  rhs <- sub("^\\(?.*\\|+\\s*", "", sub("\\)$", "", re_terms))
  unique(unlist(lapply(gsub("[[:space:]]", "", rhs), function(g) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    vapply(
      seq_along(parts),
      function(i) paste(parts[seq_len(i)], collapse = ":"),
      character(1)
    )
  })))
}

#' model_spec_of
#'
#' The specification of one fitted model (a brmsfit, or a list with the same
#' fields: formula = list(formula, pforms), family = list(family, link),
#' prior = data.frame, and the coefficient names in `variables`).
#' @return list(terms = tibble(component (mu | <dpar>), kind (fixed | re |
#'   offset), term), re_groups, family, link, response, fixef, priors (tibble
#'   prior, class, coef, group, dpar, source))
model_spec_of <- function(fit) {
  bf <- fit$formula
  f <- if (inherits(bf, "formula")) bf else bf$formula
  split_terms <- function(fo, component) {
    tt <- stats::terms(fo, keep.order = TRUE)
    lab <- attr(tt, "term.labels")
    offs <- if (is.null(attr(tt, "offset"))) {
      character(0)
    } else {
      vapply(
        attr(tt, "offset"),
        function(i) deparse1(attr(tt, "variables")[[i + 1L]]),
        character(1)
      )
    }
    is_re <- grepl("\\|", lab)
    tibble::tibble(
      component = component,
      kind = c(ifelse(is_re, "re", "fixed"), rep("offset", length(offs))),
      term = .norm_term(c(ifelse(is_re, paste0("(", lab, ")"), lab), offs))
    )
  }
  terms <- split_terms(f, "mu")
  for (dp in names(bf$pforms %||% list())) {
    terms <- dplyr::bind_rows(terms, split_terms(bf$pforms[[dp]], dp))
  }
  fam <- fit$family %||% bf$family
  vars <- if (inherits(fit, "brmsfit")) {
    brms::variables(fit)
  } else {
    as.character(fit$variables %||% character(0))
  }
  pr <- if (inherits(fit, "brmsfit")) {
    as.data.frame(brms::prior_summary(fit))
  } else {
    as.data.frame(fit$prior %||% data.frame())
  }
  keep <- intersect(
    c("prior", "class", "coef", "group", "dpar", "source"),
    names(pr)
  )
  list(
    terms = terms,
    re_groups = .re_groups_of(terms$term[
      terms$kind == "re" & terms$component == "mu"
    ]),
    family = as.character(fam$family %||% NA_character_),
    link = as.character(fam$link %||% NA_character_),
    response = deparse1(f[[2]]),
    fixef = sub("^b_", "", grep("^b_", vars, value = TRUE)),
    priors = tibble::as_tibble(pr[, keep, drop = FALSE])
  )
}

#' One track's entries of model_terms_deployed.yml, expanded ({k} over
#' 1..count; {nbr} matched by pattern), normalised: list(fixed (fixed_effects +
#' baseline), re, re_groups, offsets, shape, family, link).
.deployed_terms_expected <- function(spec_track) {
  ex <- function(entries) {
    unlist(lapply(entries %||% list(), function(e) {
      if (is.list(e)) {
        vapply(
          seq_len(as.integer(e$count)),
          function(k) gsub("{k}", k, e$term, fixed = TRUE),
          character(1)
        )
      } else {
        as.character(e)
      }
    }))
  }
  list(
    fixed = .norm_term(c(
      ex(spec_track$fixed_effects),
      ex(spec_track$baseline)
    )),
    re = .norm_term(ex(spec_track$re_terms)),
    re_groups = as.character(unlist(spec_track$re_groups)),
    offsets = .norm_term(ex(spec_track$offsets)),
    shape = .norm_term(ex(spec_track$shape_terms)),
    family = as.character(spec_track$family %||% NA),
    link = as.character(spec_track$link %||% NA)
  )
}

#' The spec of one fit against its track's yml entry: one row per item
#' (fixed, re, re_groups, offsets, shape, family, link) with the expected and
#' found sets (";"-joined), the missing / unexpected members and `ok`. The
#' neighbour-contagion column ({nbr}) matches any found `nbr_lag1_*` term.
model_spec_check <- function(spec, expected) {
  found <- list(
    fixed = spec$terms$term[
      spec$terms$component == "mu" & spec$terms$kind == "fixed"
    ],
    re = spec$terms$term[
      spec$terms$component == "mu" & spec$terms$kind == "re"
    ],
    re_groups = spec$re_groups,
    offsets = spec$terms$term[
      spec$terms$component == "mu" & spec$terms$kind == "offset"
    ],
    shape = spec$terms$term[
      spec$terms$component != "mu" & spec$terms$kind != "offset"
    ],
    family = spec$family,
    link = spec$link
  )
  nbr <- unique(regmatches(
    unlist(found),
    regexpr("nbr_lag1_[A-Za-z0-9_]+", unlist(found))
  ))
  dplyr::bind_rows(lapply(names(found), function(k) {
    want <- expected[[k]] %||% character(0)
    if (length(nbr) == 1L) {
      want <- gsub("{nbr}", nbr, want, fixed = TRUE)
    }
    got <- found[[k]] %||% character(0)
    tibble::tibble(
      item = k,
      expected = paste(sort(want), collapse = "; "),
      found = paste(sort(got), collapse = "; "),
      missing = paste(sort(setdiff(want, got)), collapse = "; "),
      unexpected = paste(sort(setdiff(got, want)), collapse = "; "),
      ok = setequal(want, got)
    )
  }))
}

#' The report_model_spec tables of one track from the per-group specs
#' (list group -> model_spec_of()): terms, priors, fixef, check (vs the yml)
#' and one summary row per fit.
model_spec_tables <- function(specs, track, lags, yml = NULL) {
  exp_t <- if (is.null(yml) || is.null(yml$tracks[[track]])) {
    NULL
  } else {
    .deployed_terms_expected(yml$tracks[[track]])
  }
  key <- function(d, g) {
    dplyr::mutate(
      tibble::as_tibble(d),
      track = track,
      enso_group = g,
      lag = as.integer(lags[[g]] %||% NA),
      .before = 1L
    )
  }
  gs <- names(specs)
  check <- if (is.null(exp_t)) {
    tibble::tibble()
  } else {
    dplyr::bind_rows(lapply(gs, function(g) {
      key(model_spec_check(specs[[g]], exp_t), g)
    }))
  }
  fits <- dplyr::bind_rows(lapply(gs, function(g) {
    s <- specs[[g]]
    t <- s$terms
    key(
      tibble::tibble(
        response = s$response,
        family = s$family,
        link = s$link,
        n_fixed = sum(t$kind == "fixed" & t$component == "mu"),
        n_re = sum(t$kind == "re" & t$component == "mu"),
        re_groups = paste(s$re_groups, collapse = "; "),
        offsets = paste(t$term[t$kind == "offset"], collapse = "; "),
        n_fixef = length(s$fixef),
        n_priors = nrow(s$priors),
        matches_deployed_terms = if (nrow(check)) {
          all(check$ok[check$enso_group == g])
        } else {
          NA
        }
      ),
      g
    )
  }))
  list(
    fits = fits,
    terms = dplyr::bind_rows(lapply(gs, function(g) key(specs[[g]]$terms, g))),
    priors = dplyr::bind_rows(lapply(gs, function(g) {
      key(specs[[g]]$priors, g)
    })),
    fixef = dplyr::bind_rows(lapply(gs, function(g) {
      key(tibble::tibble(coef = specs[[g]]$fixef), g)
    })),
    check = check
  )
}

#' Convergence diagnostics of every deployed best-lag fit of one track
#' (utils/diagnostics.R::fit_diagnostics()): per-parameter R-hat / ESS, the
#' per-chain sampler table, sampler settings and the package versions stored
#' in each fit, plus one scalar row per fit.
#' @return list(scalars, summary, sampler, settings, versions, thresholds)
derive_convergence <- function(ctx, track) {
  groups <- manifest_groups_all(ctx$MAN)
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  present <- Filter(
    function(g) {
      ff <- track_fit_files(mnf, track, g)
      ok <- file.exists(ff$fit)
      if (!ok) {
        warning(sprintf(
          "[09] convergence (%s): missing fit for %s: %s",
          track,
          g,
          ff$fit
        ))
      }
      ok
    },
    groups
  )
  message(sprintf(
    "[09] convergence (%s): groups %s",
    track,
    paste(present, collapse = ", ")
  ))
  per_group <- map_workers(
    present,
    function(g) {
      ff <- track_fit_files(mnf, track, g)
      message(sprintf("[09]   loading %s | %s", track, g))
      fit <- readRDS(ff$fit)
      d <- fit_diagnostics(fit)
      spec <- model_spec_of(fit) # report_model_spec, from the same load
      fit <- NULL
      key <- function(x) {
        dplyr::mutate(
          x,
          track = track,
          enso_group = g,
          lag = as.integer(ff$lag),
          .before = 1
        )
      }
      gc(verbose = FALSE)
      list(
        scalars = key(fit_diagnostics_row(d)),
        summary = key(d$summary),
        sampler = key(d$sampler),
        settings = key(d$settings),
        versions = key(d$versions),
        thresholds = d$thresholds,
        spec = spec,
        lag = as.integer(ff$lag)
      )
    },
    workers = ctx$WORKERS
  )
  if (!length(per_group)) {
    stop(sprintf("[09] convergence (%s): no deployed fit found", track))
  }
  # The fits' specifications for ensure_report_model_spec() of this job.
  if (is.environment(ctx$MODEL_SPEC_STASH)) {
    assign(
      track,
      list(
        specs = stats::setNames(lapply(per_group, `[[`, "spec"), present),
        lags = stats::setNames(lapply(per_group, `[[`, "lag"), present)
      ),
      envir = ctx$MODEL_SPEC_STASH
    )
  }
  pick <- function(k) dplyr::bind_rows(lapply(per_group, `[[`, k))
  list(
    scalars = pick("scalars"),
    summary = pick("summary"),
    sampler = pick("sampler"),
    settings = pick("settings"),
    versions = pick("versions"),
    thresholds = per_group[[1]]$thresholds
  )
}

ensure_convergence <- function(ctx, track) {
  fit_files <- track_fit_input_files(ctx, track)$fit
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$GROUPING_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      fit_files
    )
  )
  out <- ensure_artefact(
    ctx,
    "convergence",
    inputs,
    function() derive_convergence(ctx, track),
    track = track
  )
  # The same step writes the fit-level model specification (table S34), from
  # the fits this derivation just loaded.
  ensure_report_model_spec(ctx, track)
  invisible(out)
}

#' report_model_spec v1 (table S34) of one track: model_spec_tables() of the
#' deployed fits, from this job's convergence load (ctx$MODEL_SPEC_STASH) or,
#' when convergence was fresh, a load of its own.
derive_report_model_spec <- function(ctx, track) {
  st <- if (is.environment(ctx$MODEL_SPEC_STASH)) {
    ctx$MODEL_SPEC_STASH[[track]]
  } else {
    NULL
  }
  if (is.null(st)) {
    message(sprintf(
      "[09] report_model_spec (%s): convergence was not derived in this job; loading the fits",
      track
    ))
    mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
    gs <- Filter(
      function(g) file.exists(track_fit_files(mnf, track, g)$fit),
      manifest_groups_all(ctx$MAN)
    )
    one <- map_workers(
      gs,
      function(g) {
        ff <- track_fit_files(mnf, track, g)
        list(spec = model_spec_of(readRDS(ff$fit)), lag = as.integer(ff$lag))
      },
      workers = ctx$WORKERS
    )
    st <- list(
      specs = stats::setNames(lapply(one, `[[`, "spec"), gs),
      lags = stats::setNames(lapply(one, `[[`, "lag"), gs)
    )
  } else {
    message(sprintf(
      "[09] report_model_spec (%s): from the convergence step's fit load",
      track
    ))
  }
  yml <- if (file.exists(MODEL_TERMS_DEPLOYED)) {
    yaml::read_yaml(MODEL_TERMS_DEPLOYED)
  } else {
    NULL
  }
  out <- model_spec_tables(st$specs, track, st$lags, yml)
  bad <- dplyr::filter(out$check, !ok)
  if (nrow(bad)) {
    warning(
      sprintf(
        "[09] report_model_spec (%s): deployed fit(s) differ from %s: %s",
        track,
        MODEL_TERMS_DEPLOYED,
        paste(
          unique(sprintf("%s %s", bad$enso_group, bad$item)),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  c(out, list(model_terms_file = MODEL_TERMS_DEPLOYED, track = track)) # report_model_spec v1
}

ensure_report_model_spec <- function(ctx, track) {
  fit_files <- track_fit_input_files(ctx, track)$fit
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      fit_files,
      if (file.exists(MODEL_TERMS_DEPLOYED)) MODEL_TERMS_DEPLOYED
    )
  )
  ensure_artefact(
    ctx,
    "report_model_spec",
    inputs,
    function() derive_report_model_spec(ctx, track),
    track = track
  )
}

#' Paths of the convergence artefacts of every track (report freshness inputs).
#' @noRd
.convergence_paths <- function(ctx) {
  vapply(
    ctx$FIG1_TRACKS,
    function(t) {
      make_results_run_file(
        ctx$opt,
        "results_convergence",
        track = t,
        schema = RESULTS_SCHEMA[["convergence"]],
        config_path = ctx$CONFIG_PATH
      )
    },
    character(1)
  )
}

#' The deploy convergence gate of the results config's deploy
#' (deploy_core.R::deploy_convergence_cfg(); `allow` ids not checked against a
#' registry), or the in-function defaults (action off, scope population, the
#' DIAG_THRESHOLDS) when results.deploy_config is unset; plus the deploy block
#' (NULL without one) for the short group ids of `allow`.
.convergence_gate_cfg <- function(ctx) {
  dep_path <- ctx$RES$deploy_config %||% ""
  dep <- if (nzchar(dep_path)) load_deploy_config(dep_path) else NULL
  list(
    cfg = deploy_convergence_cfg(dep %||% list(), validate_allow = FALSE),
    dep = dep
  )
}

#' Per deployed fit the deploy-gate verdict (diagnostics.R::fit_diagnostics_gate()),
#' from the per-parameter summary and per-chain sampler table the convergence
#' artefact stores for each fit (the inputs fit_diagnostics_gate() reads from
#' fit_diagnostics(), so no fit is reloaded). gate_verdict: pass | allowed
#' (fails but its <track>:<short group> id is in deploy.convergence.allow) | fail.
#' @param conv a convergence artefact (one track); @param gc .convergence_gate_cfg()
.convergence_gate_rows <- function(conv, gc) {
  cfg <- gc$cfg
  keys <- dplyr::distinct(conv$scalars, track, enso_group)
  dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    tr <- keys$track[[i]]
    g <- keys$enso_group[[i]]
    pick <- function(x) {
      as.data.frame(x[x$track == tr & x$enso_group == g, , drop = FALSE])
    }
    v <- fit_diagnostics_gate(
      list(summary = pick(conv$summary), sampler = pick(conv$sampler)),
      scope = cfg$scope,
      rhat_max = cfg$rhat_max,
      ess_min = cfg$ess_min,
      max_divergent = cfg$max_divergent,
      rhat_any = cfg$rhat_any,
      ebfmi_min = cfg$ebfmi_min
    )
    short <- if (is.null(gc$dep)) {
      sub("^enso_", "", g)
    } else {
      .deploy_group_short(gc$dep, g)
    }
    id <- paste0(TRACK_CODE[[tr]], ":", short)
    tibble::tibble(
      track = tr,
      enso_group = g,
      gate_scope = v$scope,
      gate_max_rhat_scope = v$max_rhat_scope,
      gate_max_rhat_variable = v$max_rhat_variable,
      gate_reasons = v$reasons,
      gate_verdict = if (isTRUE(v$pass)) {
        "pass"
      } else if (id %in% cfg$allow) {
        "allowed"
      } else {
        "fail"
      }
    )
  }))
}

#' Table S8: one row per deployed fit (track x group) with the scalar
#' diagnostics, the sampler settings and (v2) the deploy-gate
#' verdict of each fit (.convergence_gate_rows(): gate_verdict pass | allowed
#' | fail, gate_reasons, the scope's worst R-hat and its parameter) under the
#' deploy's convergence block; `gate` = that block (action, scope, thresholds,
#' allow).
derive_report_convergence <- function(ctx) {
  conv <- lapply(ctx$FIG1_TRACKS, function(t) ensure_convergence(ctx, t))
  gc <- .convergence_gate_cfg(ctx)
  fits <- dplyr::bind_rows(lapply(conv, function(x) {
    dplyr::left_join(
      x$scalars,
      dplyr::select(
        x$settings,
        track,
        enso_group,
        iter,
        warmup,
        adapt_delta,
        max_treedepth
      ),
      by = c("track", "enso_group")
    )
  }))
  gate <- dplyr::bind_rows(lapply(conv, .convergence_gate_rows, gc = gc))
  fits <- dplyr::left_join(fits, gate, by = c("track", "enso_group"))
  list(
    fits = fits,
    thresholds = conv[[1]]$thresholds,
    tracks = ctx$FIG1_TRACKS,
    gate = gc$cfg[c(
      "action",
      "scope",
      "rhat_max",
      "ess_min",
      "max_divergent",
      "allow",
      "rhat_any",
      "ebfmi_min"
    )]
  ) # report_convergence v2
}

ensure_report_convergence <- function(ctx) {
  dep <- ctx$RES$deploy_config %||% ""
  ensure_artefact(
    ctx,
    "report_convergence",
    results_inputs(files = c(.convergence_paths(ctx), if (nzchar(dep)) dep)),
    function() derive_report_convergence(ctx)
  )
}

# --- Lag selection and LOO rankings ------------------------------------------

#' Tables S9 / S10 and fig. S6: the stacked LOO lag ranking of every track
#' with the deploy rule's pick marked, read from the enso_relation artefacts
#' (lag_rank, built by lag_rank_with_selection() -- the ranking Figure 1
#' outlines, so the tables cannot disagree with it); one row per track x group
#' with the deployed lag, the rule's pick, the rank-1 lag and the deployed
#' lag's elpd difference; and the firm forecast horizon of every group and
#' country from the forecast manifest.
#' @return list(ranking, selected, horizons, countries, rule)
# --- Lag floor, report_lags v2 rows ---------

#' The smallest lag the deploy's selection rule can pick for one track / group
#' (the eligible set is single_lag >= floor). The rule's min_lag; for duration
#' / monthly, under rule best_elpd with an occ_floor_offset, also >= occ_lag -
#' offset (select_best_lag_best_elpd()), and under longest_within_se with
#' enforce_occ_floor, >= occ_lag (select_best_lag_downstream()'s eligible
#' fallback). -Inf = unbounded.
#' @param rule lag_selection_cfg() list; @param occ_lag the group's occurrence lag
lag_floor_of <- function(
  rule,
  track,
  occ_lag = NULL,
  enforce_occ_floor = FALSE
) {
  floor <- as.numeric(rule$min_lag)
  down <- !identical(track, "occurrence")
  if (down && !is.null(occ_lag) && !is.na(occ_lag)) {
    if (identical(rule$rule, "best_elpd") && !is.null(rule$occ_floor_offset)) {
      floor <- max(
        floor,
        as.numeric(occ_lag) - as.numeric(rule$occ_floor_offset)
      )
    }
    if (
      identical(rule$rule, "longest_within_se") && isTRUE(enforce_occ_floor)
    ) {
      floor <- max(floor, as.numeric(occ_lag))
    }
  }
  floor
}

#' results_lag_floor -- lag_floor_of() under the results config's deploy rule
#' (results_lag_selection()) and results.derivation.enforce_occ_floor.
results_lag_floor <- function(ctx, track, occ_lag = NULL) {
  lag_floor_of(
    results_lag_selection(ctx),
    track,
    occ_lag,
    isTRUE(ctx$ENFORCE_OCC_FLOOR)
  )
}

#' One report_lags `selected` row (v2) from one group's ranking (rank,
#' single_lag, elpd_diff, se_diff vs the rank-1 lag, `selected`; Pareto-k
#' columns optional): v1's columns plus the deploy floor, whether it binds
#' (rank-1 lag below it) and its paired elpd cost, the rank-2 lag and the best
#' lag >= 0 ("no lead") with their elpd difference to rank 1, and the deployed
#' lag's Pareto-k counts when the ranking has them.
report_lags_selected_row <- function(r, track, g, deployed_lag, floor) {
  r <- r[order(r$rank), , drop = FALSE]
  dl <- as.integer(deployed_lag)
  dep <- r[r$single_lag == dl, , drop = FALSE]
  sel <- r$single_lag[r$selected %in% TRUE]
  r1 <- as.integer(r$single_lag[1L])
  r2 <- if (nrow(r) >= 2L) r[2L, , drop = FALSE] else r[0L, , drop = FALSE]
  nn <- r[r$single_lag >= 0, , drop = FALSE][1L, , drop = FALSE]
  fin <- is.finite(floor)
  active <- fin && r1 < floor
  pk <- function(col) {
    if (col %in% names(dep) && nrow(dep)) dep[[col]][1L] else NA
  }
  tibble::tibble(
    track = track,
    enso_group = g,
    deployed_lag = dl,
    selected_lag = if (length(sel)) as.integer(sel[1]) else NA_integer_,
    rank1_lag = r1,
    elpd_diff_deployed = if (nrow(dep)) dep$elpd_diff[1] else NA_real_,
    se_diff_deployed = if (nrow(dep)) dep$se_diff[1] else NA_real_,
    n_lags = nrow(r),
    floor_lag = if (fin) as.integer(floor) else NA_integer_,
    bound_active = active,
    bound_cost = if (active && nrow(dep)) dep$elpd_diff[1] else NA_real_,
    bound_cost_se = if (active && nrow(dep)) dep$se_diff[1] else NA_real_,
    rank2_lag = if (nrow(r2)) as.integer(r2$single_lag) else NA_integer_,
    elpd_diff_rank2 = if (nrow(r2)) r2$elpd_diff else NA_real_,
    se_diff_rank2 = if (nrow(r2)) r2$se_diff else NA_real_,
    best_nonneg_lag = if (nrow(nn) && !is.na(nn$single_lag)) {
      as.integer(nn$single_lag)
    } else {
      NA_integer_
    },
    elpd_diff_nonneg = if (nrow(nn)) nn$elpd_diff else NA_real_,
    se_diff_nonneg = if (nrow(nn)) nn$se_diff else NA_real_,
    n_k_gt07_deployed = as.integer(pk("n_k_gt07")),
    k_max_deployed = as.numeric(pk("k_max"))
  )
}

derive_report_lags <- function(ctx) {
  groups <- manifest_groups_all(ctx$MAN)
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  rule <- results_lag_selection(ctx)
  occ_dl <- vapply(
    groups,
    function(g) {
      tryCatch(
        as.integer(track_fit_files(mnf, "occurrence", g)$lag),
        error = function(e) NA_integer_
      )
    },
    integer(1)
  )
  ranking <- list()
  selected <- list()
  psrc <- character(0)
  for (trk in ctx$FIG1_TRACKS) {
    rk <- ensure_enso_relation(ctx, trk)$lag_rank
    if (is.null(rk) || !nrow(rk)) {
      stop(
        sprintf(
          "[09] report_lags: the %s enso_relation artefact has no lag ranking.",
          trk
        ),
        call. = FALSE
      )
    }
    rk <- dplyr::mutate(rk, track = trk, .before = 1)
    # Pareto-k columns (v2): passed through when the ranking carries them
    # (05), else NA (utils/loo_lags.R::ranking_pareto()).
    rk <- dplyr::bind_rows(lapply(unique(rk$enso_group), function(g) {
      x <- rk[rk$enso_group == g, , drop = FALSE]
      y <- ranking_pareto(as.data.frame(x))
      psrc[[paste(trk, x$enso_group[1])]] <<- attr(y, "pareto_source")
      tibble::as_tibble(y)
    }))
    ranking[[trk]] <- rk
    selected[[trk]] <- dplyr::bind_rows(lapply(
      intersect(groups, unique(rk$enso_group)),
      function(g) {
        r <- dplyr::filter(rk, enso_group == g)
        dl <- as.integer(track_fit_files(mnf, trk, g)$lag)
        report_lags_selected_row(
          r,
          trk,
          g,
          dl,
          results_lag_floor(ctx, trk, occ_dl[[g]])
        )
      }
    ))
  }
  man <- ctx$MAN
  horizons <- dplyr::bind_rows(lapply(groups, function(g) {
    x <- man$groups[[g]]
    if (is.null(x)) {
      return(NULL)
    }
    tibble::tibble(
      enso_group = g,
      best_lag_occ = as.integer(x$best_lag_occ %||% NA),
      best_lag_dur = as.integer(x$best_lag_dur %||% NA),
      best_lag_mon = as.integer(x$best_lag_mon %||% NA),
      min_lag = as.integer(x$min_lag %||% NA),
      horizon_end = as.character(x$horizon_end %||% NA),
      n_countries = as.integer(x$n_countries %||% NA)
    )
  }))
  horizons$last_cfsv2_month <- as.character(man$last_cfsv2_month %||% NA)
  # The firm horizons' ONI end (additive): = last_cfsv2_month
  # under the legacy rule, the month before under seam last_full_oni.
  horizons$firm_oni_end <- as.character(manifest_last_firm_oni(man) %||% NA)
  horizons$firm_horizon_rule <- as.character(
    man[["firm_horizon_rule"]] %||% "last_cfsv2 + min lag"
  )
  fh <- manifest_country_horizons(man, "firm")
  cg <- manifest_country_groups(man)
  countries <- tibble::tibble(
    iso3 = names(fh),
    enso_group = unname(cg[names(fh)]),
    firm_horizon_end = unname(fh)
  )
  list(
    ranking = dplyr::bind_rows(ranking),
    selected = dplyr::bind_rows(selected),
    horizons = horizons,
    countries = countries,
    rule = rule[c("rule", "se_mult", "min_lag", "occ_floor_offset")],
    pareto_source = psrc,
    # Table T2 (v3): the per-cell LOO summary at the selection rule's se_mult.
    loo_summary = report_lags_loo_summary(
      dplyr::bind_rows(ranking),
      dplyr::bind_rows(selected),
      as.numeric(rule$se_mult)
    )
  ) # report_lags v3
}

#' report_lags_loo_summary
#'
#' Table T2 (supplement), per track x group, from the ranking (elpd_diff /
#' se_diff relative to the rank-1 lag of the cell): rank-1 and deployed lag,
#' the deployed lag's cost, rank 2, lag 0, the worst lag, the number of lags
#' with elpd_diff >= -se_mult x se_diff (rank 1 included), and the rank-1 fit's
#' elpd_loo (se_elpd) and p_loo. No no-ENSO column.
report_lags_loo_summary <- function(ranking, selected, se_mult = 2) {
  keys <- unique(ranking[, c("track", "enso_group")])
  col <- function(r, k) if (k %in% names(r) && nrow(r)) r[[k]][1] else NA_real_
  dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    r <- ranking[
      ranking$track == keys$track[i] & ranking$enso_group == keys$enso_group[i],
      ,
      drop = FALSE
    ]
    r <- r[order(r$rank), , drop = FALSE]
    dlg <- selected$deployed_lag[
      selected$track == keys$track[i] &
        selected$enso_group == keys$enso_group[i]
    ]
    dlg <- if (length(dlg)) as.integer(dlg[1]) else NA_integer_
    dep <- r[r$single_lag %in% dlg, , drop = FALSE]
    r2 <- if (nrow(r) >= 2L) r[2L, , drop = FALSE] else r[0L, , drop = FALSE]
    l0 <- r[r$single_lag == 0L, , drop = FALSE]
    w <- r[which.min(r$elpd_diff), , drop = FALSE]
    tibble::tibble(
      track = keys$track[i],
      enso_group = keys$enso_group[i],
      n_lags = nrow(r),
      rank1_lag = as.integer(r$single_lag[1]),
      deployed_lag = dlg,
      deployed_elpd_diff = col(dep, "elpd_diff"),
      deployed_se_diff = col(dep, "se_diff"),
      rank2_lag = as.integer(col(r2, "single_lag")),
      rank2_elpd_diff = col(r2, "elpd_diff"),
      rank2_se_diff = col(r2, "se_diff"),
      lag0_elpd_diff = col(l0, "elpd_diff"),
      lag0_se_diff = col(l0, "se_diff"),
      worst_lag = as.integer(col(w, "single_lag")),
      worst_elpd_diff = col(w, "elpd_diff"),
      worst_se_diff = col(w, "se_diff"),
      n_within_se = sum(r$elpd_diff >= -se_mult * r$se_diff, na.rm = TRUE),
      se_mult = se_mult,
      elpd_loo_rank1 = col(r, "elpd_loo"),
      se_elpd_rank1 = col(r, "se_elpd"),
      p_loo_rank1 = col(r, "p_loo")
    )
  }))
}

ensure_report_lags <- function(ctx) {
  rel <- vapply(
    ctx$FIG1_TRACKS,
    function(t) {
      make_results_run_file(
        ctx$opt,
        "results_enso_relation",
        track = t,
        schema = RESULTS_SCHEMA[["enso_relation"]],
        config_path = ctx$CONFIG_PATH
      )
    },
    character(1)
  )
  dep <- ctx$RES$deploy_config %||% ""
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      rel,
      if (nzchar(dep)) dep
    )
  )
  ensure_artefact(ctx, "report_lags", inputs, function() {
    derive_report_lags(ctx)
  })
}

# --- Lag averaging over the single-lag fits --
#
# lag_averaging (v1, track-keyed, wave 1; on with results.derivation.
# lag_averaging.enable): from the 05 lag sidecars (05_fit_postprocess_loo_lags.R
# --lag_sidecar: per (batch, group) the N x L pointwise elpd and Pareto-k of
# every single-lag fit and each lag's population crossbasis draws), per group:
#   lags      per lag: elpd_loo, Pareto-k counts, the paired elpd difference
#             (and SE) to the rank-1 lag and to the deployed lag, the k-robust
#             difference (rows with k > k_threshold at ANY lag dropped, rank 1
#             re-found on the rest), and the weight of every method x lag set
#             (stacking = loo::stacking_weights(), pseudobma_plus =
#             loo::pseudobma_weights(BB = TRUE), seeded; lag sets `all` and
#             `eligible` = single_lag >= the deploy floor, results_lag_floor())
#   erf       the population ERF on one grid (the intersection of every lag's
#             ONI range, 200 points incl. 0, centred exactly at ONI 0), per lag
#             and for each method x lag set's mixture: n_irr_draws draws whose
#             lag labels are systematic samples by weight (seeded); the
#             deployed / rank-2 / unbounded rank-1 lags are rows of the per-lag
#             table (roles in `roles`)
#   irr       IRR at enso_conds of every lag and mixture, median, 95% CrI, P(IRR > 1)
# Only the population ERF and the IRR at the ENSO conditions are averaged; no
# ADM1 map, no attribution at an averaged lag. A track with no sidecar is
# skipped with a message.

LAG_AVG_SEED <- 20261027L
LAG_AVG_GRID_N <- 200L
LAG_AVG_METHODS <- c("stacking", "pseudobma_plus")
LAG_AVG_SETS <- c("all", "eligible")

#' The validated results.derivation.lag_averaging block.
lag_avg_cfg <- function(x) {
  x <- x %||% list()
  m <- as.character(unlist(x$methods %||% LAG_AVG_METHODS))
  ls <- as.character(unlist(x$lag_sets %||% LAG_AVG_SETS))
  if (!length(m) || !all(m %in% LAG_AVG_METHODS)) {
    stop(
      "[09] results.derivation.lag_averaging.methods must be among ",
      paste(LAG_AVG_METHODS, collapse = ", "),
      call. = FALSE
    )
  }
  if (!length(ls) || !all(ls %in% LAG_AVG_SETS)) {
    stop(
      "[09] results.derivation.lag_averaging.lag_sets must be among ",
      paste(LAG_AVG_SETS, collapse = ", "),
      call. = FALSE
    )
  }
  kt <- as.numeric(x$k_threshold %||% 0.7)
  if (length(kt) != 1L || !is.finite(kt) || kt <= 0) {
    stop(
      "[09] results.derivation.lag_averaging.k_threshold must be one positive number.",
      call. = FALSE
    )
  }
  list(
    enable = isTRUE(x$enable),
    methods = m,
    lag_sets = ls,
    k_threshold = kt,
    report_min_weight = as.numeric(x$report_min_weight %||% 0.01)
  )
}

#' Paired elpd difference of every column of E (N x L pointwise elpd) to
#' column `ref`: sum and SE = sqrt(N var(diff)) (loo::loo_compare()'s SE).
lag_elpd_diff <- function(E, ref) {
  N <- nrow(E)
  d <- E - E[, ref]
  list(diff = colSums(d), se = sqrt(N * apply(d, 2, stats::var)))
}

#' k-robust elpd differences: rows with Pareto k > thr at any lag dropped,
#' rank 1 re-found on the remaining rows, paired differences to it.
lag_krobust <- function(E, K, thr) {
  keep <- rowSums(K > thr, na.rm = TRUE) == 0L
  if (!any(keep)) {
    return(list(
      n_dropped = nrow(E),
      rank1 = NA_integer_,
      diff = rep(NA_real_, ncol(E)),
      se = rep(NA_real_, ncol(E))
    ))
  }
  Ek <- E[keep, , drop = FALSE]
  r1 <- which.max(colSums(Ek))
  d <- lag_elpd_diff(Ek, r1)
  list(n_dropped = sum(!keep), rank1 = r1, diff = d$diff, se = d$se)
}

#' Lag weights of one method on the columns `cols` of E (0 elsewhere).
lag_weights <- function(E, cols, method, seed = LAG_AVG_SEED) {
  w <- numeric(ncol(E))
  if (length(cols) == 1L) {
    w[cols] <- 1
    return(w)
  }
  ww <- switch(
    method,
    stacking = as.numeric(loo::stacking_weights(E[, cols, drop = FALSE])),
    pseudobma_plus = .with_seed(
      seed,
      as.numeric(loo::pseudobma_weights(E[, cols, drop = FALSE], BB = TRUE))
    ),
    stop("lag_weights: unknown method ", method, call. = FALSE)
  )
  w[cols] <- ww / sum(ww)
  w
}

#' Population log-IRR draws of one single-lag fit at exposures `x`: the fit's
#' ENSO onebasis (cb_meta enso entry) at x minus at 0, times the S x K
#' population coefficient draws -> length(x) x S (exactly 0 at x = 0).
lag_erf_logrr <- function(vm, b_draws, x) {
  Zx <- build_cb_columns_from_Q(as.numeric(x), vm)
  Z0 <- build_cb_columns_from_Q(0, vm)
  Zc <- Zx - matrix(Z0[1L, ], nrow(Zx), ncol(Zx), byrow = TRUE)
  Zc %*% t(b_draws[, colnames(Zx), drop = FALSE])
}

#' The common ERF grid: the intersection of the lags' ONI ranges, n points,
#' plus 0 when inside (so the centring point is on the grid).
lag_common_grid <- function(ranges, n = LAG_AVG_GRID_N) {
  lo <- max(vapply(ranges, `[`, 0, 1L))
  hi <- min(vapply(ranges, `[`, 0, 2L))
  if (!(hi > lo)) {
    stop("lag_common_grid: the lags' ONI ranges do not overlap.", call. = FALSE)
  }
  g <- seq(lo, hi, length.out = if (lo < 0 && hi > 0) n - 1L else n)
  sort(unique(c(g, if (lo < 0 && hi > 0) 0)))
}

#' Evenly spaced thinning of S draws to n (all when S <= n).
.thin_idx <- function(S, n) {
  if (S <= n) seq_len(S) else unique(round(seq(1, S, length.out = n)))
}

#' Mixture draws over lags: M draws whose lag labels are systematic samples by
#' weight w (one uniform offset, seeded), each taking column perm[j] of its
#' lag's draw matrix (a seeded permutation of 1..M; every matrix has >= M
#' columns). With w = (1, 0, ...) the mixture is lag 1's draws, permuted.
#' @param mats list of (n_x x >= M) draw matrices, one per lag (w's order)
#' @return n_x x M matrix, attributes lag_index and draw_index
lag_mixture_draws <- function(mats, w, M, seed = LAG_AVG_SEED) {
  w <- as.numeric(w)
  stopifnot(length(mats) == length(w), all(w >= 0), sum(w) > 0)
  cw <- cumsum(w) / sum(w)
  u <- .with_seed(seed, list(u0 = stats::runif(1L) / M, perm = sample.int(M)))
  pos <- u$u0 + (seq_len(M) - 1L) / M
  a <- pmin(findInterval(pos, c(0, cw), rightmost.closed = TRUE), length(w))
  a <- vapply(a, function(i) if (w[i] > 0) i else which(w > 0)[1L], integer(1)) # guard ties at 0
  out <- matrix(NA_real_, nrow(mats[[1L]]), M)
  for (l in unique(a)) {
    j <- which(a == l)
    if (ncol(mats[[l]]) < M) {
      stop(
        "lag_mixture_draws: lag ",
        l,
        " has fewer than M draws.",
        call. = FALSE
      )
    }
    out[, j] <- mats[[l]][, u$perm[j], drop = FALSE]
  }
  structure(out, lag_index = a, draw_index = u$perm)
}

#' Pointwise ERF summary of a log-IRR draw matrix (rows = grid points).
.lag_erf_summary <- function(mat, at) {
  q <- function(p) apply(mat, 1L, stats::quantile, p, names = FALSE)
  tibble::tibble(
    at = at,
    logRR = apply(mat, 1L, stats::median),
    lo = q(0.025),
    hi = q(0.975),
    lo50 = q(0.25),
    hi50 = q(0.75)
  )
}

#' IRR at the ENSO conditions from a log-IRR draw matrix (rows = conds).
.lag_irr_rows <- function(mat, conds) {
  tibble::tibble(
    enso_label = names(conds),
    oni = unname(conds),
    irr = exp(apply(mat, 1L, stats::median)),
    irr_lo = exp(apply(mat, 1L, stats::quantile, 0.025, names = FALSE)),
    irr_hi = exp(apply(mat, 1L, stats::quantile, 0.975, names = FALSE)),
    p_gt1 = rowMeans(mat > 0)
  )
}

#' Lag averaging of one (track, group) sidecar (pure; no file access).
#' @param sc a read_lag_loo_sidecar() object
#' @param deployed_lag the group's deployed lag of the track
#' @param floor its deploy floor (results_lag_floor(); -Inf = unbounded)
#' @param cfg lag_avg_cfg(); @param conds named ONI conditions
#' @param n_draws draws of each per-lag ERF and of each mixture
#' @return list(lags, erf, irr, roles, grid, n_obs, n_dropped_k)
lag_averaging_group <- function(
  sc,
  deployed_lag,
  floor,
  cfg,
  conds,
  n_draws,
  seed = LAG_AVG_SEED
) {
  check_lag_sidecar_rows(sc)
  E <- sc$elpd_pointwise
  K <- sc$pareto_k
  lags <- as.integer(sc$lags)
  labs <- names(sc$lags) %||% colnames(E)
  if (!identical(colnames(E), labs)) {
    stop(
      "lag_averaging_group: sidecar columns and lags disagree.",
      call. = FALSE
    )
  }
  dl <- as.integer(deployed_lag)
  dep <- match(dl, lags)
  if (is.na(dep)) {
    stop(
      sprintf(
        "lag_averaging_group: the deployed lag %d is not in the sidecar (%s / %s).",
        dl,
        sc$track,
        sc$group
      ),
      call. = FALSE
    )
  }
  elpd <- colSums(E)
  r1 <- which.max(elpd)
  ord <- order(-elpd)
  d1 <- lag_elpd_diff(E, r1)
  dd <- lag_elpd_diff(E, dep)
  kr <- lag_krobust(E, K, cfg$k_threshold)
  pk <- sc$pareto[match(lags, sc$pareto$single_lag), , drop = FALSE]
  elig <- which(lags >= floor)
  tab <- tibble::tibble(
    single_lag = lags,
    label = labs,
    elpd_loo = unname(elpd),
    rank = match(seq_along(lags), ord),
    elpd_diff = unname(d1$diff),
    se_diff = unname(d1$se),
    elpd_diff_deployed = unname(dd$diff),
    se_diff_deployed = unname(dd$se),
    elpd_diff_krobust = unname(kr$diff),
    se_diff_krobust = unname(kr$se),
    n_k_gt07 = as.integer(pk$n_k_gt07),
    n_k_gt1 = as.integer(pk$n_k_gt1),
    k_max = as.numeric(pk$k_max),
    share_k_bad = as.numeric(pk$share_k_bad),
    k_flag = as.integer(colSums(K > cfg$k_threshold, na.rm = TRUE)),
    eligible = seq_along(lags) %in% elig,
    deployed = seq_along(lags) == dep,
    rank1 = seq_along(lags) == r1,
    rank2 = seq_along(lags) == ord[2L]
  )
  W <- list()
  for (m in cfg$methods) {
    for (ls in cfg$lag_sets) {
      cols <- if (identical(ls, "all")) seq_along(lags) else elig
      if (!length(cols)) {
        next
      }
      w <- lag_weights(E, cols, m, seed)
      tab[[sprintf("w_%s_%s", m, ls)]] <- w
      W[[paste(m, ls, sep = ":")]] <- w
    }
  }
  # Per-lag draws on the common grid + the ENSO conditions.
  grid <- lag_common_grid(lapply(sc$cb, function(x) {
    as.numeric(x$cb_meta$range)
  }))
  x <- c(grid, unname(conds))
  ng <- length(grid)
  M <- as.integer(n_draws)
  mats <- lapply(sc$cb[labs], function(cb) {
    b <- cb$b_draws
    lag_erf_logrr(cb$cb_meta, b[.thin_idx(nrow(b), M), , drop = FALSE], x)
  })
  M <- min(M, vapply(mats, ncol, 0L))
  erf <- list()
  irr <- list()
  for (i in seq_along(lags)) {
    mt <- mats[[i]]
    erf[[length(erf) + 1L]] <- dplyr::mutate(
      .lag_erf_summary(mt[seq_len(ng), , drop = FALSE], grid),
      variant = "lag",
      method = NA_character_,
      lag_set = NA_character_,
      single_lag = lags[i],
      .before = 1L
    )
    irr[[length(irr) + 1L]] <- dplyr::mutate(
      .lag_irr_rows(mt[ng + seq_along(conds), , drop = FALSE], conds),
      variant = "lag",
      method = NA_character_,
      lag_set = NA_character_,
      single_lag = lags[i],
      .before = 1L
    )
  }
  for (k in names(W)) {
    mx <- lag_mixture_draws(
      lapply(mats, function(m) m[, seq_len(M), drop = FALSE]),
      W[[k]],
      M,
      seed
    )
    p <- strsplit(k, ":", fixed = TRUE)[[1]]
    erf[[length(erf) + 1L]] <- dplyr::mutate(
      .lag_erf_summary(mx[seq_len(ng), , drop = FALSE], grid),
      variant = "mixture",
      method = p[1],
      lag_set = p[2],
      single_lag = NA_integer_,
      .before = 1L
    )
    irr[[length(irr) + 1L]] <- dplyr::mutate(
      .lag_irr_rows(mx[ng + seq_along(conds), , drop = FALSE], conds),
      variant = "mixture",
      method = p[1],
      lag_set = p[2],
      single_lag = NA_integer_,
      .before = 1L
    )
  }
  list(
    lags = tab,
    erf = dplyr::bind_rows(erf),
    irr = dplyr::bind_rows(irr),
    roles = tibble::tibble(
      role = c("deployed", "rank1", "rank2", "rank1_krobust"),
      single_lag = c(
        dl,
        lags[r1],
        lags[ord[2L]],
        if (is.na(kr$rank1)) NA_integer_ else lags[kr$rank1]
      )
    ),
    grid = grid,
    n_obs = nrow(E),
    n_dropped_k = kr$n_dropped,
    n_draws = M,
    floor_lag = if (is.finite(floor)) as.integer(floor) else NA_integer_
  )
}

#' The latest lag sidecar file of every manifest group for a track (NA where none).
.lag_sidecar_files <- function(ctx, track) {
  batch <- LAG_LOO_DIR(ctx, track)
  vapply(
    manifest_groups_all(ctx$MAN),
    function(g) {
      latest_lag_sidecar_file(batch, group = g, opt = ctx$opt) %||%
        NA_character_
    },
    character(1)
  )
}

derive_lag_averaging <- function(ctx, track) {
  cfg <- ctx$LAG_AVG
  batch <- LAG_LOO_DIR(ctx, track)
  files <- .lag_sidecar_files(ctx, track)
  groups <- names(files)[!is.na(files)]
  code <- TRACK_CODE[[track]]
  mnf <- tryCatch(
    load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR),
    error = function(e) NULL
  )
  dep_lag <- function(g, trk) {
    v <- ctx$MAN$groups[[g]][[paste0("best_lag_", TRACK_CODE[[trk]])]]
    if (is.null(v) && !is.null(mnf)) {
      v <- track_fit_files(mnf, trk, g)$lag
    }
    as.integer(v)
  }
  out <- list()
  for (g in groups) {
    sc <- read_lag_loo_sidecar(batch, group = g, opt = ctx$opt)
    if (!identical(sc$track, track) || !identical(as.character(sc$group), g)) {
      stop(
        sprintf(
          "[09] lag_averaging: sidecar %s is for %s / %s, not %s / %s.",
          attr(sc, "file"),
          sc$track,
          sc$group,
          track,
          g
        ),
        call. = FALSE
      )
    }
    occ <- if (identical(track, "occurrence")) {
      NULL
    } else {
      dep_lag(g, "occurrence")
    }
    message(sprintf(
      "[09] lag_averaging %s | %s: %d lags, N = %d",
      track,
      g,
      length(sc$lags),
      sc$n_obs
    ))
    r <- lag_averaging_group(
      sc,
      dep_lag(g, track),
      results_lag_floor(ctx, track, occ),
      cfg,
      ctx$ENSO_CONDS,
      ctx$N_IRR_DRAWS
    )
    tag <- function(d) {
      dplyr::mutate(d, track = track, enso_group = g, .before = 1L)
    }
    out[[g]] <- list(
      lags = tag(r$lags),
      erf = tag(r$erf),
      irr = tag(r$irr),
      roles = tag(r$roles),
      meta = tibble::tibble(
        track = track,
        enso_group = g,
        n_obs = r$n_obs,
        n_dropped_k = r$n_dropped_k,
        n_draws = r$n_draws,
        floor_lag = r$floor_lag,
        grid_lo = min(r$grid),
        grid_hi = max(r$grid),
        sidecar = attr(sc, "file")
      )
    )
  }
  pick <- function(f) dplyr::bind_rows(lapply(out, `[[`, f))
  list(
    lags = pick("lags"),
    erf = pick("erf"),
    irr = pick("irr"),
    roles = pick("roles"),
    meta = pick("meta"),
    missing_groups = setdiff(names(files), groups),
    settings = cfg,
    seed = LAG_AVG_SEED,
    track = track,
    batch = batch
  ) # lag_averaging v1
}

ensure_lag_averaging <- function(ctx, track) {
  files <- .lag_sidecar_files(ctx, track)
  if (all(is.na(files))) {
    message(sprintf(
      paste0(
        "[09] lag_averaging (%s): no lag sidecar in %s (written by ",
        "05_fit_postprocess_loo_lags.R --lag_sidecar, HPC run H1); skipped."
      ),
      track,
      make_loo_compare_run_file(
        ctx$opt,
        "loo_compare_dir",
        batch = LAG_LOO_DIR(ctx, track)
      )
    ))
    return(invisible(NULL))
  }
  if (any(is.na(files))) {
    message(sprintf(
      "[09] lag_averaging (%s): no sidecar for %s; those groups are left out.",
      track,
      paste(names(files)[is.na(files)], collapse = ", ")
    ))
  }
  dep <- ctx$RES$deploy_config %||% ""
  dir <- make_loo_compare_run_file(
    ctx$opt,
    "loo_compare_dir",
    batch = LAG_LOO_DIR(ctx, track)
  )
  inputs <- results_inputs(
    files = c(ctx$MANIFEST_PATH, files[!is.na(files)], if (nzchar(dep)) dep),
    globs = list(list(dir = dir, pattern = LAG_LOO_SIDECAR_RE))
  )
  ensure_artefact(
    ctx,
    "lag_averaging",
    inputs,
    function() derive_lag_averaging(ctx, track),
    track = track
  )
}

# --- Group display lags and per-country season ends ---------
#
# group_labels (v1, unkeyed, wave 1): the numbers 10 / 11 build the lag-class
# display labels from (results_core.R::group_display_meta(), style.group_labels
# "lag") -- each group's deployed occ / dur / mon lag READ FROM THE FITS
# (paths.R::read_deployed_lag() on the best-lag archive, asserted equal to the
# archive MANIFEST and to the forecast manifest's best_lag_*), its members and
# counts -- plus the per-country season ends of table S30 (season_ends,
# season_windows). Ids only: no display string is stored.

#' Each group's deployed lag per track from the deployed fit files, checked
#' against the archive MANIFEST.csv (track_fit_files()) and the forecast
#' manifest (groups.<g>.best_lag_<code>). A missing fit is NA for duration /
#' monthly and an error for occurrence (the display label needs it).
#' @return tibble(enso_group, lag_occ, lag_dur, lag_mon)
group_deployed_lags <- function(ctx, groups) {
  grouping <- as.character(
    ctx$MAN$grouping %||% sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  )
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  code <- c(occurrence = "occ", duration = "dur", monthly = "mon")
  rows <- lapply(groups, function(g) {
    v <- vapply(
      names(code),
      function(trk) {
        dl <- read_deployed_lag(
          ctx$ENSO_BEST_LAG_DIR,
          g,
          grouping,
          code[[trk]],
          episode_key = episode_key_of(ctx$opt),
          national_zero_fill = national_zero_fill_of(ctx$opt),
          national_zero_threshold = national_zero_threshold_of(ctx$opt),
          occ_panel = occ_panel_spec_of(ctx$opt)
        )
        if (is.null(dl)) {
          if (identical(trk, "occurrence")) {
            stop(
              sprintf(
                "[09] group_labels: no deployed occurrence fit for group %s in %s (grouping %s).",
                g,
                ctx$ENSO_BEST_LAG_DIR,
                grouping
              ),
              call. = FALSE
            )
          }
          return(NA_integer_)
        }
        lag <- as.integer(dl$lag)
        arch <- tryCatch(
          as.integer(track_fit_files(mnf, trk, g)$lag),
          error = function(e) NA_integer_
        )
        if (!is.na(arch) && !identical(arch, lag)) {
          stop(
            sprintf(
              "[09] group_labels: %s / %s: the fit file says lag %d, the archive MANIFEST %d.",
              trk,
              g,
              lag,
              arch
            ),
            call. = FALSE
          )
        }
        man <- ctx$MAN$groups[[g]][[paste0("best_lag_", code[[trk]])]]
        if (!is.null(man) && !identical(as.integer(man), lag)) {
          stop(
            sprintf(
              "[09] group_labels: %s / %s: the fit file says lag %d, the forecast manifest %d.",
              trk,
              g,
              lag,
              as.integer(man)
            ),
            call. = FALSE
          )
        }
        lag
      },
      integer(1)
    )
    tibble::tibble(
      enso_group = g,
      lag_occ = v[["occurrence"]],
      lag_dur = v[["duration"]],
      lag_mon = v[["monthly"]]
    )
  })
  dplyr::bind_rows(rows)
}

#' Per-country season-end table (S30) from a gid calendar.
#' @param se_tbl season_calendar.R::season_end_by_country() output, or NULL
#'   (no panel: the season columns are NA)
#' @param country_group named character ISO3 -> fit group
#' @param group_month named integer group -> the former group season-end month
#'   (deploy hindcast.season_end_by_group), or NULL
#' @param gid_source the survival panel's .min_seas_source.csv rows (gid,
#'   min_seas_source), or NULL: calendar_default stays NA until 03 writes it
#' @param gid_country named character gid -> ISO3 (for gid_source)
#' @return tibble(iso3, enso_group, season_end_month, min_seas_mode,
#'   modal_share, n_gid, n_distinct, group_season_end_month, n_gid_default,
#'   calendar_default), one row per country of country_group
season_ends_country_table <- function(
  se_tbl,
  country_group,
  group_month = NULL,
  gid_source = NULL,
  gid_country = NULL
) {
  iso <- sort(names(country_group))
  i <- if (is.null(se_tbl)) {
    rep(NA_integer_, length(iso))
  } else {
    match(iso, se_tbl$iso3)
  }
  col <- function(f, na) {
    if (is.null(se_tbl)) rep(na, length(iso)) else se_tbl[[f]][i]
  }
  out <- tibble::tibble(
    iso3 = iso,
    enso_group = unname(country_group[iso]),
    season_end_month = as.integer(col("season_end_month", NA_integer_)),
    min_seas_mode = as.integer(col("min_seas_mode", NA_integer_)),
    modal_share = as.numeric(col("modal_share", NA_real_)),
    n_gid = as.integer(col("n_gid", NA_integer_)),
    n_distinct = as.integer(col("n_distinct", NA_integer_)),
    group_season_end_month = if (is.null(group_month)) {
      NA_integer_
    } else {
      as.integer(unname(group_month[unname(country_group[iso])]))
    }
  )
  if (is.null(gid_source) || is.null(gid_country)) {
    out$n_gid_default <- NA_integer_
    out$calendar_default <- NA
  } else {
    src <- data.frame(
      gid = as.character(gid_source$gid),
      dflt = as.character(gid_source$min_seas_source) == "default_january",
      stringsAsFactors = FALSE
    )
    src$iso3 <- unname(gid_country[src$gid])
    src <- src[!is.na(src$iso3), , drop = FALSE]
    n_def <- tapply(src$dflt, src$iso3, sum)
    n_all <- tapply(src$dflt, src$iso3, length)
    out$n_gid_default <- unname(as.integer(n_def[iso]))
    out$calendar_default <- unname(
      as.integer(n_def[iso]) == as.integer(n_all[iso])
    )
  }
  out
}

#' Per-country windows of a hindcast or event manifest: per window (launch
#' year / event) x country, the country's season end and firm / unit horizons
#' (manifest_country_horizons(): a country's own entry first, else its
#' group's / unit's). One row per country of the window's units.
#' @param man hindcast manifest (read_hindcast_manifest() / yaml) or event
#'   manifest (read_event_manifest())
#' @param type "hindcast" | "event"
#' @return tibble(window_type, window, iso3, enso_group, unit, launch,
#'   season_end, firm_horizon_end, unit_horizon_end)
season_windows_table <- function(man, type = c("hindcast", "event")) {
  type <- match.arg(type)
  one <- function(m, window, launch_of, se_country, se_group, year = NULL) {
    cg <- manifest_country_groups(m)
    if (!length(cg)) {
      return(NULL)
    }
    fh <- manifest_country_horizons(m, "firm", year = year)
    uh <- manifest_country_horizons(m, "unit", year = year)
    iso <- names(cg)
    un <- vapply(
      iso,
      function(i) as.character(m$countries[[i]]$unit %||% cg[[i]]),
      character(1)
    )
    la <- vapply(
      iso,
      function(i) as.character(launch_of(cg[[i]]) %||% NA),
      character(1)
    )
    se <- vapply(
      iso,
      function(i) as.character(se_country(i) %||% se_group(cg[[i]]) %||% NA),
      character(1)
    )
    tibble::tibble(
      window_type = type,
      window = as.character(window),
      iso3 = iso,
      enso_group = unname(cg),
      unit = unname(un),
      launch = as.Date(unname(la)),
      season_end = as.Date(unname(se)),
      firm_horizon_end = as.Date(unname(fh[iso])),
      unit_horizon_end = as.Date(unname(uh[iso]))
    )
  }
  rows <- if (identical(type, "event")) {
    lapply(names(man$events), function(ev) {
      e <- man$events[[ev]]
      one(
        event_manifest_view(man, ev),
        ev,
        function(g) e$launch,
        function(i) e$countries[[i]]$season_end,
        function(g) e$groups[[g]]$season_end
      )
    })
  } else {
    yrs <- sort(unique(unlist(lapply(man$groups, function(g) names(g$years)))))
    lapply(yrs, function(y) {
      one(
        man,
        y,
        function(g) man$groups[[g]]$years[[y]]$launch,
        function(i) {
          ((man$countries[[i]][["years"]] %||% list())[[y]] %||%
            list())$season_end
        },
        function(g) man$groups[[g]]$years[[y]]$season_end,
        year = y
      )
    })
  }
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    return(out)
  }
  dplyr::arrange(out, window_type, window, iso3)
}

#' The former group season-end months (deploy hindcast.season_end_by_group):
#' the hindcast manifest's groups.<g>.season_end_month when it carries them
#' (by_group), else the results deploy config's hindcast.season_end_by_group.
.group_season_months <- function(ctx, hc_man) {
  m <- if (!is.null(hc_man)) {
    unlist(lapply(hc_man$groups, function(x) x$season_end_month))
  } else {
    NULL
  }
  if (length(m)) {
    return(stats::setNames(as.integer(m), names(m)))
  }
  dep <- ctx$RES$deploy_config %||% ""
  if (!nzchar(dep) || !file.exists(dep)) {
    return(NULL)
  }
  s <- as.character(
    yaml::read_yaml(dep)$deploy$hindcast$season_end_by_group %||% ""
  )
  if (!nzchar(s)) {
    return(NULL)
  }
  kv <- strsplit(strsplit(s, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
  stats::setNames(
    as.integer(vapply(kv, `[`, "", 2L)),
    trimws(vapply(kv, `[`, "", 1L))
  )
}

#' Inputs of group_labels' season-end tables: the hindcast manifest (when the
#' hindcast is on or its file exists), the event manifest, the gid calendar's
#' panel (the hindcast manifest's season_end_source panel, else its
#' climate_panel, else the event manifest's) and the survival panel's
#' .min_seas_source.csv (03_prepare_occ.R; absent until its next run).
.group_labels_season_inputs <- function(ctx) {
  hp <- as.character(ctx$HINDCAST_MANIFEST_PATH %||% "")
  hc <- if (nzchar(hp) && file.exists(hp)) yaml::read_yaml(hp) else NULL
  ep <- as.character(ctx$EVENT_MANIFEST_PATH %||% "")
  em <- if (nzchar(ep) && file.exists(ep)) read_event_manifest(ep) else NULL
  panel <- as.character(
    hc$season_end_source$panel %||%
      hc$climate_panel %||%
      em$climate_panel %||%
      ""
  )
  list(
    hc = hc,
    hc_path = if (!is.null(hc)) hp else "",
    em = em,
    em_path = if (!is.null(em)) ep else "",
    panel = panel,
    src_csv = sub("\\.rds$", ".min_seas_source.csv", ctx$SURV_PATH)
  )
}

derive_group_labels <- function(ctx) {
  meta <- results_group_meta(ctx)
  groups <- as.character(meta$order)
  lags <- group_deployed_lags(ctx, groups)
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  members <- tibble::tibble(
    enso_group = as.character(reg$group),
    iso3 = as.character(reg$iso3)
  ) %>%
    dplyr::filter(enso_group %in% groups) %>%
    dplyr::arrange(match(enso_group, groups), iso3)
  cnt <- dplyr::count(members, enso_group, name = "n_countries")
  grp <- lags %>%
    dplyr::mutate(
      sidecar_order = match(enso_group, groups),
      .after = enso_group
    ) %>%
    dplyr::left_join(cnt, by = "enso_group") %>%
    dplyr::mutate(
      n_countries = dplyr::coalesce(n_countries, 0L),
      in_manifest = enso_group %in% names(ctx$MAN$groups)
    )

  si <- .group_labels_season_inputs(ctx)
  cal <- if (nzchar(si$panel) && file.exists(si$panel)) {
    panel_season_calendar(readRDS(si$panel), si$panel)
  } else {
    message(
      "[09] group_labels: no climate panel (",
      if (nzchar(si$panel)) si$panel else "none named",
      "): season-end columns NA."
    )
    NULL
  }
  se_tbl <- if (!is.null(cal)) season_end_by_country(cal) else NULL
  src <- if (file.exists(si$src_csv)) {
    utils::read.csv(si$src_csv, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  if (is.null(src)) {
    message(
      "[09] group_labels: ",
      si$src_csv,
      " absent (03_prepare_occ.R writes it at its next run): ",
      "calendar_default NA."
    )
  }
  gid_country <- if (!is.null(cal)) {
    stats::setNames(cal$country, cal$gid)
  } else {
    NULL
  }
  cg <- stats::setNames(members$enso_group, members$iso3)
  season_ends <- season_ends_country_table(
    se_tbl,
    cg,
    .group_season_months(ctx, si$hc),
    src,
    gid_country
  )
  windows <- dplyr::bind_rows(
    if (!is.null(si$hc)) season_windows_table(si$hc, "hindcast"),
    if (!is.null(si$em)) season_windows_table(si$em, "event")
  )
  list(
    groups = grp,
    members = members,
    season_ends = season_ends,
    season_windows = windows,
    season_end_mode = list(
      hindcast = as.character(si$hc$season_end %||% "by_group"),
      event = as.character(si$em$season_end %||% "by_group")
    ),
    season_end_source = list(
      panel = si$panel,
      min_seas_digest = if (!is.null(cal)) {
        as.character(attr(cal, "digest"))
      } else {
        NA_character_
      },
      min_seas_source_csv = if (!is.null(src)) si$src_csv else NA_character_,
      rule = paste(
        "season_end_month = min_seas - 1; gid -> country: modal min_seas,",
        "ties -> the later end (the larger month 1-12)"
      )
    ),
    grouping = as.character(meta$source),
    best_lag_dir = ctx$ENSO_BEST_LAG_DIR
  ) # group_labels v1
}

ensure_group_labels <- function(ctx) {
  si <- .group_labels_season_inputs(ctx)
  side <- sub("\\.csv$", ".yml", ctx$GROUPING_PATH)
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$ENSO_BEST_LAG_MANIFEST,
      ctx$GROUPING_PATH,
      if (file.exists(side)) side,
      si$hc_path,
      si$em_path,
      if (nzchar(si$panel) && file.exists(si$panel)) si$panel,
      if (file.exists(si$src_csv)) si$src_csv
    )
  )
  ensure_artefact(ctx, "group_labels", inputs, function() {
    derive_group_labels(ctx)
  })
}

# --- Baseline-rate anchoring -------------------------------------------------

#' Resolved options of the per-lag config a deployed fit of `track` / `group`
#' was fitted from (results.deploy_config's fits.<track>_dir batch), without
#' loading the fit.
deployed_track_opt <- function(ctx, track, group) {
  dep_path <- ctx$RES$deploy_config %||% ""
  if (!nzchar(dep_path) || !file.exists(dep_path)) {
    stop(
      "[09] report_anchor: results.deploy_config is unset or missing; the anchor tables ",
      "are resolved from the deployed fit configs.",
      call. = FALSE
    )
  }
  dep <- yaml::read_yaml(dep_path)$deploy
  spec <- .deploy_track_specs()[[TRACK_CODE[[track]]]]
  mnf <- load_best_lag_manifest(ctx$ENSO_BEST_LAG_DIR)
  lag <- track_fit_files(mnf, track, group)$lag
  cfg <- .deploy_cfg_path(
    file.path("analysis/configs", dep$fits[[spec$batch_key]]),
    group,
    lag
  )
  if (!file.exists(cfg)) {
    stop(
      "[09] report_anchor: deployed fit config not found: ",
      cfg,
      call. = FALSE
    )
  }
  o <- NULL
  invisible(utils::capture.output(
    o <- suppressMessages(make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", cfg),
      model = spec$model
    ))
  ))
  o
}

#' Base-rate table path of each track (one continent-wide table per track; the
#' first manifest group's config resolves it). NA when the track is not anchored.
#' @noRd
.anchor_table_paths <- function(ctx) {
  g <- manifest_groups_all(ctx$MAN)[1]
  vapply(
    ctx$FIG1_TRACKS,
    function(t) {
      o <- deployed_track_opt(ctx, t, g)
      code <- TRACK_CODE[[t]]
      if (!anchor_enabled(o, code)) {
        return(NA_character_)
      }
      make_run_file(o, anchor_spec(code)$kind)
    },
    character(1)
  )
}

#' Table S7: the empirical baseline-rate tables of the three tracks
#' (utils/anchor.R::build_base_rate_table(), as persisted at fit time): every
#' country x admin-level cell with its tier (own estimate / zero events /
#' borrowed), the continental rates, and the between-cell variance tau2 and
#' level contrasts.
derive_report_anchor <- function(ctx) {
  paths <- .anchor_table_paths(ctx)
  scal <- function(x) {
    if (is.null(x) || !length(x)) {
      NA
    } else if (length(x) == 1L) {
      x
    } else {
      paste(format(x), collapse = " / ")
    }
  }
  out <- list(cell = list(), global = list(), meta = list())
  for (t in names(paths)) {
    p <- paths[[t]]
    if (is.na(p)) {
      message("[09] report_anchor: ", t, " is not anchored; skipped.")
      next
    }
    require_input(p, "baseline-rate table")
    brt <- readRDS(p)
    out$cell[[t]] <- dplyr::mutate(
      tibble::as_tibble(brt$cell),
      track = t,
      .before = 1
    )
    out$global[[t]] <- dplyr::mutate(
      tibble::as_tibble(brt$global),
      track = t,
      .before = 1
    )
    m <- brt$meta
    out$meta[[t]] <- tibble::tibble(
      track = t,
      file = basename(p),
      tau2 = as.numeric(scal(m$tau2)),
      estimator = as.character(scal(m$estimator)),
      n_cells = as.integer(scal(m$n_cells)),
      n_countries = as.integer(scal(m$n_countries)),
      level_factor = list(m$level_factor),
      level_factor_ci = list(m$level_factor_ci),
      n_own = sum(brt$cell$mr_tier == "T1", na.rm = TRUE),
      n_zero = sum(brt$cell$mr_tier == "T1-zero", na.rm = TRUE),
      n_borrowed = sum(!brt$cell$mr_tier %in% c("T1", "T1-zero"), na.rm = TRUE)
    )
  }
  list(
    cell = dplyr::bind_rows(out$cell),
    global = dplyr::bind_rows(out$global),
    meta = dplyr::bind_rows(out$meta),
    files = paths
  )
}

ensure_report_anchor <- function(ctx) {
  paths <- .anchor_table_paths(ctx)
  dep <- ctx$RES$deploy_config %||% ""
  inputs <- results_inputs(
    files = c(
      stats::na.omit(paths),
      ctx$ENSO_BEST_LAG_MANIFEST,
      if (nzchar(dep)) dep
    )
  )
  ensure_artefact(ctx, "report_anchor", inputs, function() {
    derive_report_anchor(ctx)
  })
}

# --- Software ----------------------------------------------------------------

#' The forward runs whose simulation environment table S13 reports: the
#' displayed forecast factual, the attribution factual when it is a run of its
#' own, and the counterfactual of every analog set; named list arm -> list of
#' manifest_runs() entries (list(config, merged_dir), named by unit).
.report_software_sim_runs <- function(ctx) {
  cf <- lapply(ctx$SETS, function(a) {
    manifest_runs(ctx$MAN, "counterfactual", a)
  })
  c(
    list(forecast = manifest_runs(ctx$MAN_DISPLAY %||% ctx$MAN, "factual")),
    if (identical(ctx$ATTRIB_FACTUAL, "attrib")) {
      list(attribution = manifest_runs(ctx$MAN_ATTRIB, "factual"))
    },
    list(counterfactual = do.call(c, unname(cf)))
  )
}

#' report_software_sim_versions
#'
#' The simulation environment of the forward runs: the
#' `software` 07b spatPomp records in each run's partition_meta.yml
#' (provenance.R::sim_software_versions()), one row per arm x unit x software,
#' distinct. A run without the field (older runs) or without a partition_meta.yml here
#' gives one NA row.
#' @return tibble(arm, unit, software, version)
report_software_sim_versions <- function(runs) {
  rows <- list()
  for (arm in names(runs)) {
    for (i in seq_along(runs[[arm]])) {
      # by index: units repeat across sets
      un <- names(runs[[arm]])[i] %||% NA_character_
      pm <- partition_meta_path(runs[[arm]][[i]]$merged_dir)
      sw <- if (!is.na(pm)) yaml::read_yaml(pm)$software else NULL
      rows[[length(rows) + 1L]] <- if (length(sw)) {
        tibble::tibble(
          arm = arm,
          unit = un,
          software = names(sw),
          version = unname(vapply(
            sw,
            function(v) as.character(v %||% NA_character_)[1L],
            character(1)
          ))
        )
      } else {
        tibble::tibble(
          arm = arm,
          unit = un,
          software = NA_character_,
          version = NA_character_
        )
      }
    }
  }
  if (!length(rows)) {
    return(tibble::tibble(
      arm = character(0),
      unit = character(0),
      software = character(0),
      version = character(0)
    ))
  }
  dplyr::distinct(dplyr::bind_rows(rows))
}

#' Table S13: R and package versions of the session deriving the results
#' (provenance.R::software_versions()), the versions each deployed fit
#' recorded at fit time (the computing cluster's) and the simulation
#' environment of the forward runs (`sim_versions`; additive, no schema bump:
#' an artefact may lack it).
derive_report_software <- function(ctx) {
  conv <- lapply(ctx$FIG1_TRACKS, function(t) ensure_convergence(ctx, t))
  fv <- dplyr::bind_rows(lapply(conv, function(x) {
    dplyr::left_join(
      x$versions,
      dplyr::select(x$settings, track, enso_group, backend),
      by = c("track", "enso_group")
    )
  }))
  list(
    session = software_versions(),
    fit_versions = dplyr::distinct(dplyr::select(fv, -enso_group, -lag)),
    fit_versions_by_fit = fv,
    sim_versions = report_software_sim_versions(.report_software_sim_runs(ctx))
  )
}

ensure_report_software <- function(ctx) {
  pm <- vapply(
    unlist(.report_software_sim_runs(ctx), recursive = FALSE),
    function(r) partition_meta_path(r$merged_dir),
    character(1)
  )
  ensure_artefact(
    ctx,
    "report_software",
    results_inputs(
      files = c(.convergence_paths(ctx), unique(stats::na.omit(pm)))
    ),
    function() derive_report_software(ctx)
  )
}

# --- brms_kl class-sweep comparison ------------------------------------------

#' Files of the class-sweep comparison (05_fit_postprocess_kl_sweep.R output).
#' @noRd
.kl_sweep_files <- function(ctx) {
  d <- ctx$KL_SWEEP_DIR
  file.path(
    d,
    c(
      "fits.csv",
      "classes.csv",
      "partition.csv",
      "loo.csv",
      "diagnostics_chains.csv",
      "loo.rds"
    )
  )
}

#' The grouping-model comparison: per sweep fit its modes, convergence on all /
#' agreeing chains, country-block LOO (elpd difference, SE, stacking weight,
#' Pareto k) and MAP partition with its agreement with the deployed fit
#' registry (results.inputs.grouping).
derive_report_kl_sweep <- function(ctx) {
  f <- .kl_sweep_files(ctx)
  for (p in f) {
    require_input(p, "kl sweep comparison (05_fit_postprocess_kl_sweep.R)")
  }
  rd <- function(p) {
    tibble::as_tibble(utils::read.csv(p, stringsAsFactors = FALSE))
  }
  fits <- rd(f[1])
  classes <- rd(f[2])
  part <- rd(f[3])
  loo <- rd(f[4])
  settings <- readRDS(f[6])$settings
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  reg_short <- stats::setNames(sub("^enso_", "", reg$group), reg$iso3)
  part <- dplyr::mutate(
    part,
    deployed_group = unname(reg_short[iso3]),
    agrees = !is.na(deployed_group) & group == deployed_group
  )
  agreement <- dplyr::summarise(
    dplyr::group_by(part, id),
    n_countries = dplyr::n(),
    n_agree = sum(agrees),
    share_agree = mean(agrees),
    .groups = "drop"
  )
  src <- sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  side <- sub("\\.csv$", ".yml", ctx$GROUPING_PATH)
  list(
    fits = dplyr::left_join(fits, agreement, by = "id"),
    loo = loo,
    classes = classes,
    partition = part,
    selected = list(
      source = src,
      digest = tryCatch(registry_digest(src), error = function(e) {
        NA_character_
      }),
      sidecar_sha256 = if (file.exists(side)) {
        digest::digest(file = side, algo = "sha256")
      } else {
        NA_character_
      }
    ),
    settings = settings,
    dir = ctx$KL_SWEEP_DIR
  )
}

# --- Membership uncertainty of the deployed classes ----------
#
# report_membership (v1, unkeyed, wave 1; on when results.kl_sweep_dir AND the
# top-level results.kl_sweep_cell are set): from the extended 05 kl-sweep
# comparison (responsibilities.csv, country_mixture.csv, stability.csv,
# loo_country.csv, loo.csv, fits.csv, loo.rds), for the deployed cell:
#   country      each country's deployed group, p_max (largest mean
#                responsibility, agreeing chains), the full responsibility row
#                (r_<group>), its prior registry group, the M1 contrast / lag
#                (mixture vs hard MAP) and `carried_over` (countries the fit
#                lacks, added by fill_from: NA)
#   group        per deployed group: hard members, soft size sum_c r_ck, median /
#                min p_max (the count under style.membership_low is display)
#   stability    stability.csv (sampled partitions vs the MAP, ARI)
#   sensitivity  per sweep cell x chain set: ARI vs the deployed partition,
#                countries agreeing by label and after the optimal class
#                matching (kl_sweep.R::kl_match_classes()), the countries that
#                move, the pairwise country-block elpd difference vs the
#                deployed cell with SE = sqrt(C var(diff)) from loo_country.csv,
#                Pareto k > 0.7 counts of both cells and the kl stacking weight
#                -- descriptive (country-block PSIS-LOO; high k expected)
#   mixture_m1   country_mixture.csv rows of the cell (agreeing chains)
#   mixture_m2   per track and country: sum_k r_ck ERF_k at group k's deployed
#                lag, from the 05 lag sidecars' population draws (the reader of
#                lag_averaging), vs the hard (deployed group) ERF: IRR at
#                enso_conds and the ERF on a common grid; skipped per track with
#                a message while the sidecars are absent
# 09 stops unless the cell's agreeing-chain MAP partition is the deployed
# registry (results.inputs.grouping), naming the countries that differ.

MEMBERSHIP_SEED <- 20261027L

#' Files of report_membership (05_fit_postprocess_kl_sweep.R outputs).
.membership_files <- function(ctx) {
  file.path(
    ctx$KL_SWEEP_DIR,
    c(
      "responsibilities.csv",
      "country_mixture.csv",
      "stability.csv",
      "loo_country.csv",
      "loo.csv",
      "fits.csv",
      "partition.csv",
      "loo.rds"
    )
  )
}

#' MAP partition per (id, chains) from responsibility rows: each country's
#' class of largest mean responsibility (first on ties) and its p_max.
#' @param resp responsibilities.csv rows (id, chains, iso3, class, label, mean)
membership_map_partition <- function(resp) {
  resp <- resp[
    order(resp$id, resp$chains, resp$iso3, -resp$mean, resp$class),
    ,
    drop = FALSE
  ]
  top <- resp[!duplicated(resp[, c("id", "chains", "iso3")]), , drop = FALSE]
  tibble::tibble(
    id = top$id,
    chains = top$chains,
    iso3 = top$iso3,
    label = top$label,
    p_max = top$mean
  )
}

#' Stop unless the cell's MAP partition (agreeing chains) is the deployed
#' registry on every fitted country.
membership_check_map <- function(map_cell, reg_group) {
  dep <- unname(reg_group[map_cell$iso3])
  bad <- which(is.na(dep) | dep != map_cell$label)
  if (length(bad)) {
    stop(
      sprintf(
        paste0(
          "[09] report_membership: the agreeing-chain MAP partition of the sweep cell is not the ",
          "deployed registry for %d countr%s: %s. Point results.kl_sweep_cell at the cell the ",
          "registry was derived from (deploy.kl_sweep.cell)."
        ),
        length(bad),
        if (length(bad) == 1L) "y" else "ies",
        paste(
          sprintf(
            "%s (cell %s, registry %s)",
            map_cell$iso3[bad],
            map_cell$label[bad],
            ifelse(is.na(dep[bad]), "absent", dep[bad])
          ),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Pairwise country-block elpd difference of two sweep cells on one chain set
#' (sum of the per-country differences, SE = sqrt(C var(diff))) and their
#' Pareto k > 0.7 counts. Descriptive: country-block PSIS-LOO.
membership_elpd_diff <- function(loo_country, id, ref, chains) {
  a <- loo_country[
    loo_country$id == id & loo_country$chains == chains,
    ,
    drop = FALSE
  ]
  b <- loo_country[
    loo_country$id == ref & loo_country$chains == chains,
    ,
    drop = FALSE
  ]
  if (!nrow(a) || !nrow(b)) {
    return(list(
      diff = NA_real_,
      se = NA_real_,
      n_k07 = NA_integer_,
      n_k07_ref = NA_integer_
    ))
  }
  ab <- merge(
    a[, c("iso3", "elpd", "pareto_k")],
    b[, c("iso3", "elpd", "pareto_k")],
    by = "iso3",
    suffixes = c("", "_ref")
  )
  d <- ab$elpd - ab$elpd_ref
  list(
    diff = sum(d),
    se = sqrt(length(d) * stats::var(d)),
    n_k07 = sum(ab$pareto_k > 0.7, na.rm = TRUE),
    n_k07_ref = sum(ab$pareto_k_ref > 0.7, na.rm = TRUE)
  )
}

#' Sensitivity of the deployed partition to the sweep cell and chain set.
#' @param maps membership_map_partition() of every cell; @param dep named
#'   character iso3 -> deployed group (fitted countries); @param cell deployed id
membership_sensitivity <- function(maps, dep, loo_country, loo_tab, cell) {
  keys <- unique(maps[, c("id", "chains")])
  dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    m <- maps[
      maps$id == keys$id[i] & maps$chains == keys$chains[i],
      ,
      drop = FALSE
    ]
    m <- m[m$iso3 %in% names(dep), , drop = FALSE]
    d <- unname(dep[m$iso3])
    mt <- kl_match_classes(m$label, d)
    to <- unname(mt$map[m$label])
    moved <- which(is.na(to) | to != d)
    e <- membership_elpd_diff(loo_country, keys$id[i], cell, keys$chains[i])
    sw <- loo_tab$stacking_weight[
      loo_tab$id == keys$id[i] & loo_tab$chains == keys$chains[i]
    ]
    tibble::tibble(
      id = keys$id[i],
      chains = keys$chains[i],
      deployed_cell = identical(keys$id[i], cell),
      n_countries = nrow(m),
      ari = adjusted_rand_index(m$label, d),
      n_agree_label = sum(m$label == d),
      n_agree_matched = mt$n_agree,
      moved = if (length(moved)) {
        paste(
          sprintf(
            "%s %s->%s",
            m$iso3[moved],
            d[moved],
            ifelse(is.na(to[moved]), "unmatched", to[moved])
          ),
          collapse = ", "
        )
      } else {
        ""
      },
      elpd_diff = e$diff,
      se_diff = e$se,
      n_k07 = e$n_k07,
      n_k07_deployed = e$n_k07_ref,
      stacking_weight = if (length(sw)) sw[1] else NA_real_,
      note = "descriptive (country-block PSIS-LOO, high Pareto k expected)"
    )
  }))
}

#' Each cell's MAP class per fitted country, matched to the deployed labels
#' (kl_match_classes()): the per-country rows behind `sensitivity` (fig. MB).
membership_cell_partitions <- function(maps, dep) {
  keys <- unique(maps[, c("id", "chains")])
  dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    m <- maps[
      maps$id == keys$id[i] & maps$chains == keys$chains[i],
      ,
      drop = FALSE
    ]
    m <- m[m$iso3 %in% names(dep), , drop = FALSE]
    d <- unname(dep[m$iso3])
    mt <- kl_match_classes(m$label, d)
    tibble::tibble(
      id = keys$id[i],
      chains = keys$chains[i],
      iso3 = m$iso3,
      label = m$label,
      matched = unname(mt$map[m$label]),
      deployed = d,
      p_max = m$p_max
    )
  }))
}

#' Country table: deployed group, p_max, responsibility row, prior group, M1.
#' @param resp_cell responsibilities of the cell, agreeing chains
#' @param m1 country_mixture rows of the cell; @param reg tibble(iso3, group)
membership_country_table <- function(resp_cell, m1, reg) {
  w <- tidyr::pivot_wider(
    dplyr::select(resp_cell, iso3, label, mean),
    names_from = label,
    values_from = mean,
    names_prefix = "r_"
  )
  map <- membership_map_partition(resp_cell)
  pg <- dplyr::distinct(dplyr::select(resp_cell, iso3, prior_group))
  mix <- function(q, col, nm) {
    x <- m1[m1$quantity == q, c("iso3", col), drop = FALSE]
    names(x)[2] <- nm
    x
  }
  out <- tibble::tibble(
    iso3 = as.character(reg$iso3),
    group = as.character(reg$group)
  ) %>%
    dplyr::left_join(
      dplyr::select(map, iso3, map_label = label, p_max),
      by = "iso3"
    ) %>%
    dplyr::left_join(pg, by = "iso3") %>%
    dplyr::left_join(w, by = "iso3") %>%
    dplyr::mutate(carried_over = is.na(p_max))
  for (a in list(
    c("kcon", "mix_mean", "m1_kcon_mix"),
    c("kcon", "mix_q2.5", "m1_kcon_mix_lo"),
    c("kcon", "mix_q97.5", "m1_kcon_mix_hi"),
    c("kcon", "map_mean", "m1_kcon_map"),
    c("kcon", "diff_mean", "m1_kcon_diff"),
    c("kmu", "mix_mean", "m1_kmu_mix"),
    c("kmu", "map_mean", "m1_kmu_map")
  )) {
    out <- dplyr::left_join(out, mix(a[1], a[2], a[3]), by = "iso3")
  }
  dplyr::arrange(out, group, iso3)
}

#' Group summary: hard members (fitted), carried-over members, soft size
#' sum_c r_ck over the fitted countries, median / min p_max of the hard members.
membership_group_table <- function(country) {
  rcols <- grep("^r_", names(country), value = TRUE)
  groups <- unique(country$group)
  dplyr::bind_rows(lapply(groups, function(g) {
    x <- country[country$group == g, , drop = FALSE]
    fit <- x[!x$carried_over, , drop = FALSE]
    rc <- paste0("r_", g)
    tibble::tibble(
      group = g,
      n_hard = nrow(fit),
      n_carried_over = sum(x$carried_over),
      soft_size = if (rc %in% rcols) {
        sum(country[[rc]], na.rm = TRUE)
      } else {
        NA_real_
      },
      p_max_median = if (nrow(fit)) stats::median(fit$p_max) else NA_real_,
      p_max_min = if (nrow(fit)) min(fit$p_max) else NA_real_
    )
  }))
}

#' M2: per country the membership-weighted log-IRR sum_k r_ck * logRR_k (draws
#' paired by index across the groups' independent fits) vs the hard log-IRR of
#' its deployed group. Pure.
#' @param R countries x groups matrix of mean responsibilities (rows sum to 1)
#' @param logrr named list group -> (n_x x M) draw matrix, same x and M
#' @param hard named character country -> deployed group
#' @param x evaluation points; @param conds named ONI conditions (a subset of x)
#' @return list(irr = per country x cond: mix / hard IRR median + CrI, P(>1),
#'   erf = per country on the non-condition points: mix median + CrI)
membership_m2 <- function(R, logrr, hard, x, conds) {
  g <- colnames(R)
  if (!all(g %in% names(logrr))) {
    stop(
      "membership_m2: no ERF draws for group(s) ",
      paste(setdiff(g, names(logrr)), collapse = ", "),
      call. = FALSE
    )
  }
  ic <- match(unname(conds), x)
  irr <- list()
  erf <- list()
  q <- function(m, p) apply(m, 1L, stats::quantile, p, names = FALSE)
  for (c in rownames(R)) {
    mix <- Reduce(`+`, lapply(g, function(k) R[c, k] * logrr[[k]]))
    hd <- logrr[[hard[[c]]]]
    for (what in c("mixture", "hard")) {
      m <- if (what == "mixture") mix else hd
      irr[[length(irr) + 1L]] <- tibble::tibble(
        iso3 = c,
        erf = what,
        enso_label = names(conds),
        oni = unname(conds),
        irr = exp(apply(m[ic, , drop = FALSE], 1L, stats::median)),
        irr_lo = exp(q(m[ic, , drop = FALSE], 0.025)),
        irr_hi = exp(q(m[ic, , drop = FALSE], 0.975)),
        p_gt1 = rowMeans(m[ic, , drop = FALSE] > 0)
      )
    }
    gi <- setdiff(seq_along(x), ic)
    erf[[c]] <- tibble::tibble(
      iso3 = c,
      at = x[gi],
      logRR = apply(mix[gi, , drop = FALSE], 1L, stats::median),
      lo = q(mix[gi, , drop = FALSE], 0.025),
      hi = q(mix[gi, , drop = FALSE], 0.975),
      logRR_hard = apply(hd[gi, , drop = FALSE], 1L, stats::median)
    )
  }
  list(irr = dplyr::bind_rows(irr), erf = dplyr::bind_rows(erf))
}

#' M2 inputs of one track: each deployed group's population draws at its
#' deployed lag from the 05 lag sidecar, on a common grid (intersection of the
#' groups' ONI ranges) plus the ENSO conditions. NULL (with a message) when a
#' group's sidecar or its deployed lag is missing.
membership_m2_track <- function(ctx, track, groups) {
  batch <- LAG_LOO_DIR(ctx, track)
  code <- TRACK_CODE[[track]]
  ent <- list()
  for (g in groups) {
    sc <- read_lag_loo_sidecar(
      batch,
      group = g,
      opt = ctx$opt,
      required = FALSE
    )
    lag <- ctx$MAN$groups[[g]][[paste0("best_lag_", code)]]
    if (is.null(sc) || is.null(lag)) {
      message(sprintf(
        "[09] report_membership M2 (%s): no lag sidecar or deployed lag for %s; M2 skipped for this track.",
        track,
        g
      ))
      return(NULL)
    }
    i <- match(as.integer(lag), as.integer(sc$lags))
    if (is.na(i)) {
      stop(
        sprintf(
          "[09] report_membership M2: the deployed %s lag %d of %s is not in its sidecar.",
          track,
          as.integer(lag),
          g
        ),
        call. = FALSE
      )
    }
    ent[[g]] <- sc$cb[[names(sc$lags)[i]]]
  }
  grid <- lag_common_grid(
    lapply(ent, function(e) as.numeric(e$cb_meta$range)),
    100L
  )
  x <- c(grid, unname(ctx$ENSO_CONDS))
  idx <- lapply(ent, function(e) {
    .thin_idx(nrow(e$b_draws), as.integer(ctx$N_IRR_DRAWS))
  })
  M <- min(lengths(idx))
  list(
    x = x,
    logrr = lapply(stats::setNames(names(ent), names(ent)), function(g) {
      lag_erf_logrr(
        ent[[g]]$cb_meta,
        ent[[g]]$b_draws[idx[[g]][seq_len(M)], , drop = FALSE],
        x
      )
    }),
    lags = vapply(ent, function(e) as.integer(e$single_lag), integer(1))
  )
}

derive_report_membership <- function(ctx) {
  f <- .membership_files(ctx)
  for (p in f) {
    require_input(
      p,
      "kl sweep comparison (05_fit_postprocess_kl_sweep.R, membership outputs)"
    )
  }
  rd <- function(p) {
    tibble::as_tibble(utils::read.csv(p, stringsAsFactors = FALSE))
  }
  resp <- rd(f[1])
  m1 <- rd(f[2])
  stab <- rd(f[3])
  lc <- rd(f[4])
  lt <- rd(f[5])
  fits <- rd(f[6])
  settings <- readRDS(f[8])$settings
  cell <- ctx$KL_SWEEP_CELL
  if (!cell %in% fits$id) {
    stop(
      sprintf(
        "[09] report_membership: results.kl_sweep_cell '%s' is not a fit of %s (fits: %s).",
        cell,
        ctx$KL_SWEEP_DIR,
        paste(fits$id, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  reg <- load_enso_grouping(ctx$GROUPING_PATH)
  reg_group <- stats::setNames(as.character(reg$group), reg$iso3)
  maps <- membership_map_partition(resp)
  map_cell <- maps[maps$id == cell & maps$chains == "agree", , drop = FALSE]
  membership_check_map(map_cell, reg_group)
  resp_cell <- resp[resp$id == cell & resp$chains == "agree", , drop = FALSE]
  country <- membership_country_table(
    resp_cell,
    m1[m1$id == cell, , drop = FALSE],
    reg
  )
  dep <- reg_group[map_cell$iso3]
  sens <- membership_sensitivity(maps, dep, lc, lt, cell)
  # M2 per track (membership-weighted deployed ERFs); needs the lag sidecars.
  groups <- intersect(names(ctx$MAN$groups), unique(resp_cell$label))
  fitted <- country[!country$carried_over, , drop = FALSE]
  R <- as.matrix(fitted[, paste0("r_", groups), drop = FALSE])
  dimnames(R) <- list(fitted$iso3, groups)
  m2 <- list()
  m2_erf <- list()
  m2_skip <- character(0)
  for (trk in ctx$FIG1_TRACKS) {
    inp <- membership_m2_track(ctx, trk, groups)
    if (is.null(inp)) {
      m2_skip <- c(m2_skip, trk)
      next
    }
    r <- membership_m2(
      R,
      inp$logrr,
      stats::setNames(fitted$group, fitted$iso3),
      inp$x,
      ctx$ENSO_CONDS
    )
    m2[[trk]] <- dplyr::mutate(r$irr, track = trk, .before = 1L)
    m2_erf[[trk]] <- dplyr::mutate(r$erf, track = trk, .before = 1L)
  }
  list(
    country = country,
    group = membership_group_table(country),
    stability = stab,
    sensitivity = sens,
    cell_partition = membership_cell_partitions(maps, dep),
    mixture_m1 = m1[m1$id == cell, , drop = FALSE],
    mixture_m2 = if (length(m2)) dplyr::bind_rows(m2) else NULL,
    mixture_m2_erf = if (length(m2_erf)) dplyr::bind_rows(m2_erf) else NULL,
    m2_skipped_tracks = m2_skip,
    cell = cell,
    chains = "agree",
    agreeing = fits$agreeing[fits$id == cell],
    settings = list(
      n_agree = settings$n_agree,
      mode_tol = settings$mode_tol,
      k_threshold = 0.7,
      m1_chains = "agree",
      seed = MEMBERSHIP_SEED
    ),
    dir = ctx$KL_SWEEP_DIR,
    grouping = sub("\\.csv$", "", basename(ctx$GROUPING_PATH))
  ) # report_membership v1
}

ensure_report_membership <- function(ctx) {
  if (!nzchar(ctx$KL_SWEEP_DIR %||% "") || !nzchar(ctx$KL_SWEEP_CELL %||% "")) {
    message(
      "[09] report_membership: results.kl_sweep_dir / results.kl_sweep_cell unset; skipped."
    )
    return(invisible(NULL))
  }
  side <- sub("\\.csv$", ".yml", ctx$GROUPING_PATH)
  sidecar_globs <- lapply(ctx$FIG1_TRACKS, function(t) {
    list(
      dir = make_loo_compare_run_file(
        ctx$opt,
        "loo_compare_dir",
        batch = LAG_LOO_DIR(ctx, t)
      ),
      pattern = LAG_LOO_SIDECAR_RE
    )
  })
  inputs <- results_inputs(
    files = c(
      .membership_files(ctx),
      ctx$GROUPING_PATH,
      if (file.exists(side)) side,
      ctx$MANIFEST_PATH
    ),
    globs = sidecar_globs
  )
  # The cell is an unstamped key: an artefact derived for another cell is re-derived.
  path <- make_results_run_file(
    ctx$opt,
    "results_report_membership",
    schema = RESULTS_SCHEMA[["report_membership"]],
    config_path = ctx$CONFIG_PATH
  )
  if (file.exists(path)) {
    old <- tryCatch(qs2::qs_read(path), error = function(e) NULL)
    if (!is.null(old) && !identical(old$cell, ctx$KL_SWEEP_CELL)) unlink(path)
  }
  ensure_artefact(ctx, "report_membership", inputs, function() {
    derive_report_membership(ctx)
  })
}

# ===========================================================================
# Occurrence-panel tables
# ===========================================================================

#' report_national_zero_audit (v1, unkeyed, wave 1; table D32): the national-zero
#' audit of the panel's months file (01_preprocess_national_zero_months.R,
#' utils/national_zero.R::national_zero_audit()), one row per WHO country-year,
#' tagged with the rule and threshold, plus the rule parameters and totals.
derive_report_national_zero_audit <- function(ctx) {
  p <- require_input(ctx$NZ_AUDIT_PATH, "national-zero audit")
  months <- as.character(national_zero_months_file_of(ctx$opt))
  st <- read_prep_policy_stamp(months) %||% list()
  spec <- occ_panel_spec_of(ctx$opt)
  au <- tibble::as_tibble(utils::read.csv(p, stringsAsFactors = FALSE))
  au <- dplyr::mutate(
    au,
    rule = if (nzchar(spec$national_zero_rule)) {
      spec$national_zero_rule
    } else {
      "legacy"
    },
    zero_threshold = spec$national_zero_threshold,
    .before = 1
  )
  params <- tibble::tibble(
    rule = au$rule[1] %||% "legacy",
    zero_threshold = spec$national_zero_threshold,
    zero_day_frac = st$zero_day_frac %||% NA_real_,
    window = st$rule_window %||% st$window %||% NA_integer_,
    head_window = st$head_window %||% 0L,
    db_tolerance = st$rule_db_tolerance %||% st$db_tolerance %||% NA_real_,
    consistency_scope = st$consistency_scope %||% "threshold",
    start = st$start %||% NA_character_,
    end = st$end %||% NA_character_,
    exclude = st$exclude %||% NA_character_,
    rule_params_digest = st$rule_params_digest %||% NA_character_
  )
  totals <- au %>%
    dplyr::summarise(
      country_years = dplyr::n(),
      months_claimed = sum(months_claimed),
      months_who_below_threshold = sum(months_who_below_threshold),
      head_unclaimed = sum(head_unclaimed),
      tail_unclaimed = sum(tail_unclaimed),
      auto_excluded_exact = sum(auto_excluded_exact),
      auto_excluded_threshold = sum(auto_excluded_threshold),
      country_years_db_exceeds_who = sum(db_exceeds_who %in% TRUE),
      country_years_static_excluded = sum(static_excluded %in% TRUE)
    )
  list(
    audit = au,
    params = params,
    totals = totals,
    months_path = months,
    audit_path = p
  )
}

ensure_report_national_zero_audit <- function(ctx) {
  if (!nzchar(ctx$NZ_AUDIT_PATH %||% "")) {
    message(
      "[09] report_national_zero_audit: no national-zero audit beside the panel's months file; skipped."
    )
    return(invisible(NULL))
  }
  months <- as.character(national_zero_months_file_of(ctx$opt))
  inputs <- results_inputs(
    files = c(ctx$NZ_AUDIT_PATH, months, prep_policy_stamp_path(months))
  )
  ensure_artefact(ctx, "report_national_zero_audit", inputs, function() {
    derive_report_national_zero_audit(ctx)
  })
}

# D31 columns, in order (the supplement reads them by name).
OCC_PANEL_EXCLUSION_COLS <- c(
  "country",
  "adm1_overlap_parent_episodes",
  "adm1_overlap_orphan",
  "adm1_overlap_covered",
  "child_rows_censored",
  "rows_masked_running_episode",
  "onsets_masked_running_episode",
  "nz_rows_masked",
  "onsets_discarded_lp_second_onset",
  "onsets_discarded_after_episode_end",
  "onsets_discarded_lp_calendar_only",
  "cases_discarded_lp_second_onset",
  "onsets_redate_duplicate",
  "members_absorbed",
  "members_absorbed_after_gap",
  "adm0_rows_dropped",
  "cases_adm0_dropped"
)

#' Per-country table D31 from the overlap sidecars (pure; one row per country
#' with any exclusion, zeros where a sidecar has none).
occ_panel_exclusions_table <- function(
  overlap,
  discarded,
  adm0,
  masked = NULL
) {
  cnt <- function(d, nm) {
    if (is.null(d) || !nrow(d)) {
      tibble::tibble(country = character(), !!nm := integer())
    } else {
      dplyr::count(d, country = as.character(country), name = nm)
    }
  }
  sm <- function(d, v, nm) {
    if (is.null(d) || !nrow(d)) {
      tibble::tibble(country = character(), !!nm := numeric())
    } else {
      d %>%
        dplyr::group_by(country = as.character(country)) %>%
        dplyr::summarise(
          !!nm := sum(as.numeric(.data[[v]]), na.rm = TRUE),
          .groups = "drop"
        )
    }
  }
  ov <- overlap %||% tibble::tibble(country = character(), class = character())
  dc <- discarded %||%
    tibble::tibble(
      country = character(),
      reason = character(),
      after_episode_end = logical(),
      lp_calendar_only = logical(),
      cases = numeric()
    )
  mk <- masked %||%
    tibble::tibble(
      country = character(),
      reason = character(),
      outbreak = logical(),
      national_zero_fill = logical()
    )
  lp2 <- dc[dc$reason == "lp_second_onset", , drop = FALSE]
  run <- mk[mk$reason == "running_episode", , drop = FALSE]
  parts <- list(
    cnt(ov, "adm1_overlap_parent_episodes"),
    cnt(ov[ov$class == "orphan", , drop = FALSE], "adm1_overlap_orphan"),
    cnt(ov[ov$class == "covered", , drop = FALSE], "adm1_overlap_covered"),
    cnt(
      mk[mk$reason == "adm1_orphan_child", , drop = FALSE],
      "child_rows_censored"
    ),
    cnt(run, "rows_masked_running_episode"),
    cnt(
      run[run$outbreak %in% TRUE, , drop = FALSE],
      "onsets_masked_running_episode"
    ),
    cnt(mk[mk$national_zero_fill %in% TRUE, , drop = FALSE], "nz_rows_masked"),
    cnt(lp2, "onsets_discarded_lp_second_onset"),
    cnt(
      lp2[lp2$after_episode_end %in% TRUE, , drop = FALSE],
      "onsets_discarded_after_episode_end"
    ),
    cnt(
      lp2[lp2$lp_calendar_only %in% TRUE, , drop = FALSE],
      "onsets_discarded_lp_calendar_only"
    ),
    sm(lp2, "cases", "cases_discarded_lp_second_onset"),
    cnt(
      dc[dc$reason == "redate_duplicate", , drop = FALSE],
      "onsets_redate_duplicate"
    ),
    cnt(dc[dc$reason == "absorbed_member", , drop = FALSE], "members_absorbed"),
    cnt(
      dc[dc$reason == "absorbed_after_gap", , drop = FALSE],
      "members_absorbed_after_gap"
    ),
    if (is.null(adm0) || !nrow(adm0)) {
      tibble::tibble(
        country = character(),
        adm0_rows_dropped = integer(),
        cases_adm0_dropped = numeric()
      )
    } else {
      dplyr::transmute(
        adm0,
        country = as.character(country),
        adm0_rows_dropped = as.integer(n_adm0_episodes),
        cases_adm0_dropped = as.numeric(cases_adm0)
      )
    }
  )
  cs <- sort(unique(unlist(lapply(parts, `[[`, "country"))))
  out <- tibble::tibble(country = cs)
  for (p in parts) {
    out <- dplyr::left_join(out, p, by = "country")
  }
  for (nm in setdiff(OCC_PANEL_EXCLUSION_COLS, names(out))) {
    out[[nm]] <- 0L
  }
  out <- out[, OCC_PANEL_EXCLUSION_COLS]
  dplyr::mutate(out, dplyr::across(-country, ~ dplyr::coalesce(.x, 0)))
}

#' report_occ_panel_exclusions (v1, unkeyed, wave 1; table D31): what the
#' occurrence panel leaves out, by country, from the 03 sidecars of the survival
#' panel: ADM1-ADM2 overlap parent episodes (orphan / covered),
#' masked rows, discarded onsets by reason, ADM0 rows.
derive_report_occ_panel_exclusions <- function(ctx) {
  sc <- occ_panel_sidecar_paths(ctx$SURV_PATH)
  rd <- function(p) {
    if (file.exists(p)) {
      tibble::as_tibble(utils::read.csv(p, stringsAsFactors = FALSE))
    } else {
      NULL
    }
  }
  spec <- occ_panel_spec_of(ctx$opt)
  tab <- occ_panel_exclusions_table(
    rd(sc$adm1_overlap),
    rd(sc$discarded_onsets),
    rd(sc$adm0_dropped),
    rd(sc$masked_rows)
  )
  totals <- dplyr::summarise(tab, dplyr::across(-country, sum))
  list(
    by_country = tab,
    totals = totals,
    spec = tibble::tibble(
      at_risk_start = spec$at_risk_start,
      adm1_overlap = spec$adm1_overlap,
      panel = basename(ctx$SURV_PATH)
    ),
    sidecars = unlist(sc)
  )
}

ensure_report_occ_panel_exclusions <- function(ctx) {
  if (!nzchar(ctx$OCC_EXCL_PATH %||% "")) {
    message(
      "[09] report_occ_panel_exclusions: no item 11 sidecars beside the survival panel; skipped."
    )
    return(invisible(NULL))
  }
  sc <- unlist(occ_panel_sidecar_paths(ctx$SURV_PATH))
  inputs <- results_inputs(files = c(sc[file.exists(sc)], ctx$SURV_PATH))
  ensure_artefact(ctx, "report_occ_panel_exclusions", inputs, function() {
    derive_report_occ_panel_exclusions(ctx)
  })
}

# D33 columns, in order (long over the share columns; the supplement reads them by name).
CONTAGION_COVERAGE_COLS <- c(
  "country",
  "group",
  "year",
  "share_column",
  "panel_rows",
  "nbr_cov_share_mean",
  "nbr_cov_share_p10",
  "nbr_cov_share_p50",
  "rows_zero_cover",
  "own_uncovered_rows",
  "adm1_months",
  "adm1_months_covered",
  "adm1_coverage_share"
)

#' Table D33 (pure): per country x calendar year x share column, the panel rows'
#' kernel-weighted neighbour coverage at t-1 (the fit's implicit coverage; the
#' simulator's is 1, every universe unit being simulated), own-ADM1 coverage,
#' and the ADM1-month coverage of the sidecar.
#' @param panel panel rows with country, t_start, ADM1, own_cov_lag1 and the
#'   share columns `share_cols`.
#' @param side the 03 sidecar (ADM1, t_start, month, covered, in_universe).
#' @param grouping data.frame(iso3, group) or NULL.
contagion_coverage_table <- function(panel, side, share_cols, grouping = NULL) {
  grp <- if (is.null(grouping)) {
    character()
  } else {
    stats::setNames(as.character(grouping$group), as.character(grouping$iso3))
  }
  base <- tibble::tibble(
    country = as.character(panel$country),
    year = as.integer(lubridate::year(as.Date(panel$t_start))),
    own = as.numeric(panel$own_cov_lag1)
  )
  long <- dplyr::bind_rows(lapply(share_cols, function(cc) {
    dplyr::mutate(base, share_column = cc, share = as.numeric(panel[[cc]]))
  }))
  rows <- long %>%
    dplyr::group_by(country, year, share_column) %>%
    dplyr::summarise(
      panel_rows = dplyr::n(),
      nbr_cov_share_mean = mean(share),
      nbr_cov_share_p10 = as.numeric(stats::quantile(share, 0.1)),
      nbr_cov_share_p50 = as.numeric(stats::quantile(share, 0.5)),
      rows_zero_cover = sum(share == 0),
      own_uncovered_rows = sum(own == 0),
      .groups = "drop"
    )
  adm <- side %>%
    dplyr::filter(in_universe %in% TRUE) %>%
    dplyr::mutate(
      country = substr(as.character(ADM1), 1L, 3L),
      year = as.integer(lubridate::year(as.Date(month)))
    ) %>%
    dplyr::distinct(country, ADM1, month, covered, year) %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(
      adm1_months = dplyr::n(),
      adm1_months_covered = sum(covered %in% TRUE),
      .groups = "drop"
    )
  rows %>%
    dplyr::left_join(adm, by = c("country", "year")) %>%
    dplyr::mutate(
      adm1_months = dplyr::coalesce(adm1_months, 0L),
      adm1_months_covered = dplyr::coalesce(adm1_months_covered, 0L),
      adm1_coverage_share = ifelse(
        adm1_months > 0,
        adm1_months_covered / adm1_months,
        NA_real_
      ),
      group = unname(grp[country])
    ) %>%
    dplyr::select(dplyr::all_of(CONTAGION_COVERAGE_COLS)) %>%
    dplyr::arrange(share_column, country, year)
}

#' report_contagion_coverage (v1, unkeyed, wave 1; table D33): the contagion
#' coverage of the survival panel built with 03 --contagion_coverage_report
#' (stamp only -- the covariate is unchanged).
derive_report_contagion_coverage <- function(ctx) {
  side <- tibble::as_tibble(utils::read.csv(
    require_input(ctx$CONTAGION_COV_PATH, "contagion coverage"),
    stringsAsFactors = FALSE
  ))
  st <- read_prep_policy_stamp(ctx$SURV_PATH) %||% list()
  share_cols <- strsplit(as.character(st$contagion_share_cols %||% ""), ",")[[
    1L
  ]]
  surv <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
  if (!length(share_cols)) {
    share_cols <- grep("^nbr_cov_lag1_", names(surv), value = TRUE)
  }
  miss <- setdiff(c(share_cols, "own_cov_lag1"), names(surv))
  if (length(miss)) {
    stop(
      "[09] report_contagion_coverage: the survival panel lacks ",
      paste(miss, collapse = ", "),
      " (built without --contagion_coverage_report?).",
      call. = FALSE
    )
  }
  if (!"ADM1" %in% names(surv)) {
    surv$ADM1 <- gid_to_adm1(surv$gid, "report_contagion_coverage")
  }
  grouping <- if (file.exists(ctx$GROUPING_PATH %||% "")) {
    utils::read.csv(ctx$GROUPING_PATH)
  } else {
    NULL
  }
  tab <- contagion_coverage_table(surv, side, share_cols, grouping)
  by_group <- tab %>%
    dplyr::group_by(share_column, group) %>%
    dplyr::summarise(
      panel_rows = sum(panel_rows),
      nbr_cov_share_mean = sum(nbr_cov_share_mean * panel_rows) /
        sum(panel_rows),
      rows_zero_cover = sum(rows_zero_cover),
      own_uncovered_rows = sum(own_uncovered_rows),
      .groups = "drop"
    )
  stamp <- st[grep("^contagion_", names(st))]
  list(
    by_country_year = tab,
    by_group = by_group,
    stamp = stamp,
    share_cols = share_cols,
    sidecar = ctx$CONTAGION_COV_PATH
  )
}

ensure_report_contagion_coverage <- function(ctx) {
  if (!nzchar(ctx$CONTAGION_COV_PATH %||% "")) {
    message(
      "[09] report_contagion_coverage: no contagion-coverage sidecar beside the survival panel; skipped."
    )
    return(invisible(NULL))
  }
  inputs <- results_inputs(
    files = c(
      ctx$CONTAGION_COV_PATH,
      ctx$SURV_PATH,
      prep_policy_stamp_path(ctx$SURV_PATH),
      ctx$GROUPING_PATH
    )
  )
  ensure_artefact(ctx, "report_contagion_coverage", inputs, function() {
    derive_report_contagion_coverage(ctx)
  })
}

ensure_report_kl_sweep <- function(ctx) {
  if (!nzchar(ctx$KL_SWEEP_DIR %||% "")) {
    message("[09] report_kl_sweep: results.kl_sweep_dir is empty; skipped.")
    return(invisible(NULL))
  }
  inputs <- results_inputs(files = c(.kl_sweep_files(ctx), ctx$GROUPING_PATH))
  ensure_artefact(ctx, "report_kl_sweep", inputs, function() {
    derive_report_kl_sweep(ctx)
  })
}

# ===========================================================================
# af_posterior / event_af_posterior: the POSTERIOR interval of the
# attributable fraction (AF = 1 - CF / F) and number (AN = F - CF) of paired
# factual / counterfactual runs, for the forecast (Figure 4) and the past ENSO
# events (Figure 3). The paired bootstrap of report_attrib_metrics /
# report_event_attrib is the Monte-Carlo precision of the replicate means
# (width ~ 1 / sqrt(R)); this is the uncertainty of the ENSO effect: the spread
# of the expected AF over posterior draws.
#
# Needs runs with a fixed draw design (simulator.forecast.param_draw_mode:
# cycle, replicate r runs draw ((r - 1) mod max_draws) + 1, e.g. 199 draws x 50
# replicates); a run with random draws is refused. Per cell (forecast window or
# event) x analog set x explosive guard x metric x scope (Africa, ENSO group,
# country):
#   F_d, CF_d  expected window totals under draw d: the mean over the unit's
#              replicates that ran d, summed over the scope's units AT THE SAME
#              d (d indexes the group's own fit in every unit: one parameter
#              vector); AF_d = 1 - CF_d / F_d, AN_d = F_d - CF_d.
#   Africa     the groups are separate fits (independent posteriors): AF_MC
#              combinations of one draw per group.
#   noise      F_d / CF_d still carry the within-draw spread (epidemic noise,
#              the handoff state, the counterfactual ENSO trajectory). Its
#              variance is estimated directly by the delta method on the ratio
#              of means, var(AF_d) ~ sum_u var_u(cf - R_d f) / n_u / F_d^2
#              (R_d = CF_d / F_d; AN: sum_u var_u(f - cf) / n_u), and AF_d is
#              shrunk towards the pooled value by sqrt(1 - mean noise var /
#              var(AF_d)) before its quantiles are taken ("adj"; "raw" =
#              unshrunk). Validated on a synthetic known posterior
#              (tests/test_results_af_posterior.R).
# The per-replicate window totals come from the same run load as
# the other attribution steps when they run in one 09 job (ctx$AF_POST_STASH).
# ===========================================================================

AF_POST_MC <- 4000L # Africa: combinations of one draw per group
AF_POST_SEED <- 20261026L
AF_POST_MIN_REPS <- 20L # default of results.af_posterior_min_reps (replicates per (unit, draw) after the explosive drop)

#' Per-replicate window totals of one cell from a load (attrib_arm_inputs()):
#' cases and onsets per (cell, enso_group, unit, country, run, analog_set,
#' replicate_id) over months >= `window_start`, and the explosive-replicate
#' keys of that window (09's guard). A manifest without units (v1) keys the
#' unit by the group.
af_posterior_window_reps <- function(
  adm1_month,
  country_month,
  cell,
  window_start,
  cap
) {
  st <- as.Date(window_start)
  am <- dplyr::filter(adm1_month, month >= st)
  if (!"unit" %in% names(am)) {
    am$unit <- am$enso_group
  }
  cm <- dplyr::filter(
    dplyr::mutate(country_month, run = as.character(run)),
    month >= st
  )
  reps <- .sum_by_keys(
    dplyr::mutate(
      am,
      run = as.character(run),
      country = substr(adm1_unit, 1L, 3L)
    ),
    c("enso_group", "unit", "country", "run", "analog_set", "replicate_id"),
    c(cases = "cases", onsets = "outbreaks")
  ) %>%
    dplyr::mutate(cell = cell, .before = 1L)
  # Per-replicate group window maxima of the largest gid-month: the per-gid
  # explosive cap of attrib_information_explosive(), which rebuilds the group
  # totals from `reps` (no gid detail).
  gcum <- window_cum_per_rep(cm, st, "enso_group")
  gmax <- if ("max_gid_cases" %in% names(gcum)) {
    dplyr::mutate(
      dplyr::select(
        gcum,
        enso_group,
        run,
        analog_set,
        replicate_id,
        max_gid_cases
      ),
      cell = cell,
      .before = 1L
    )
  } else {
    NULL
  }
  list(
    reps = reps,
    explosive = dplyr::mutate(
      explosive_replicate_keys(cm, st, cap),
      cell = cell,
      .before = 1L
    ),
    gmax = gmax
  )
}

#' The af_posterior inputs of a load: groups and the per-window totals of the
#' attribution cells (`cells`: attrib_cell_table(); one window total per
#' window, the explosive keys at the "drop" mode's cap). Returns
#' list(groups, cells = named Date window starts, cell_table, reps, explosive).
af_posterior_stash_build <- function(inp, cells) {
  w <- cells[!duplicated(cells$window), , drop = FALSE]
  cd <- cells[cells$explosive_mode == "drop", , drop = FALSE]
  cap <- if (nrow(cd)) .cell_cap(cd[1L, ]) else Inf
  parts <- lapply(seq_len(nrow(w)), function(i) {
    af_posterior_window_reps(
      .cell_frame(inp$adm1_month, w[i, ]),
      .cell_frame(inp$country_month, w[i, ]),
      w$window[i],
      w$start[i],
      cap = cap
    )
  })
  list(
    groups = inp$groups,
    cells = .cell_starts(cells),
    cell_table = cells,
    reps = dplyr::bind_rows(lapply(parts, `[[`, "reps")),
    explosive = dplyr::bind_rows(lapply(parts, `[[`, "explosive")),
    gmax = dplyr::bind_rows(lapply(parts, `[[`, "gmax"))
  )
}

#' The explosive guards of a cell table (af_posterior_tables()'s `guards`),
#' in the cells' mode order.
.cell_guards <- function(cells) unique(cells$explosive_guard)

#' Keep a load's af_posterior inputs for a later af_posterior step of the same
#' 09 job (only when results.af_posterior is on; small: per-replicate
#' country totals, not the ADM1 x month load). The stash fills only when an
#' attribution step actually DERIVES (a fresh artefact loads nothing), and only
#' within one process: under 09a's one-task-per-step SLURM arrays every step
#' loads on its own. Building it is a group_by over the whole load per cell
#' (~10^8 rows at 9,950 replicates: minutes, not seconds).
af_posterior_stash_put <- function(ctx, key, inp, cells) {
  if (!isTRUE(ctx$AF_POSTERIOR) || is.null(ctx$AF_POST_STASH)) {
    return(invisible(NULL))
  }
  if (is.null(ctx$AF_POST_STASH[[key]])) {
    assign(key, af_posterior_stash_build(inp, cells), envir = ctx$AF_POST_STASH)
  }
  invisible(NULL)
}

#' Posterior draw of every replicate of one simulator run, from its partitions'
#' assignment.qs (simulator_core.R): tibble(replicate_id, draw).
af_posterior_run_draws <- function(sim_root, run_tag) {
  fs <- list.files(
    file.path(sim_root, run_tag, "parts"),
    pattern = "^assignment\\.qs$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(fs)) {
    stop(
      "[09] af_posterior: no assignment.qs under ",
      file.path(sim_root, run_tag, "parts"),
      call. = FALSE
    )
  }
  d <- dplyr::bind_rows(lapply(fs, function(f) {
    a <- qs2::qs_read(f)
    tibble::tibble(
      replicate_id = as.integer(a$replicate_ids),
      draw = as.integer(a$posterior_idx_occ)
    )
  }))
  if (anyDuplicated(d$replicate_id)) {
    stop(
      "[09] af_posterior: replicate ids repeated across the partitions of ",
      run_tag,
      call. = FALSE
    )
  }
  d
}

#' The parameter-draw mode a simulator config runs, as 07b resolves it
#' (spatpomp_core.R::.forecast_param_indices()): unset = the default "random"
#' (get_default_options()); blank = spread_param_draws ? "cycle" : "inherited";
#' otherwise the value as set ("random" / "cycle" / "inherited" / "joint").
#' @param sim a run config's `simulator` block (raw YAML)
.sim_param_draw_mode <- function(sim) {
  fc <- sim$forecast %||% list()
  mode <- as.character(fc[["param_draw_mode"]] %||% "random")
  if (!length(mode) || is.na(mode[1]) || !nzchar(mode[1])) {
    mode <- if (isTRUE(fc[["spread_param_draws"]] %||% FALSE)) {
      "cycle"
    } else {
      "inherited"
    }
  }
  mode[1]
}

#' af_posterior_check_design
#'
#' The fixed-draw design of one unit's arms, or an error naming what is wrong:
#' every arm config runs param_draw_mode cycle (or every arm joint) with the
#' same max_draws S, every arm's replicates ran draw ((r - 1) mod S) + 1, and
#' the arms share their replicate -> draw map. `arm_draws`: named list arm ->
#' tibble(replicate_id, draw); `arm_cfgs`: the parsed configs, same names.
#' Under joint, replicate r runs the draw of its smoothing
#' replicate ((r - 1) mod S_s) + 1, which a shared_cycle smoothing ran as draw
#' j = its replicate id: the same cycle index with S = S_s, so S must equal
#' `n_smoothing` when it is known.
#' @return tibble(replicate_id, draw); attr "max_draws" = S
af_posterior_check_design <- function(
  unit,
  arm_cfgs,
  arm_draws,
  n_smoothing = NULL
) {
  modes <- vapply(
    arm_cfgs,
    function(cf) .sim_param_draw_mode(cf$simulator %||% list()),
    character(1)
  )
  fixed <- modes %in% c("cycle", "joint")
  if (!all(fixed)) {
    stop(
      sprintf(
        paste0(
          "[09] af_posterior: unit %s arm(s) %s run param_draw_mode '%s', not cycle or ",
          "joint. The posterior interval needs a fixed draw design (e.g. 199 draws x 50 ",
          "replicates)."
        ),
        unit,
        paste(names(modes)[!fixed], collapse = ", "),
        paste(unique(modes[!fixed]), collapse = "/")
      ),
      call. = FALSE
    )
  }
  if (length(unique(modes)) != 1L) {
    stop(
      sprintf(
        "[09] af_posterior: unit %s arms mix param_draw_mode cycle and joint (%s).",
        unit,
        paste(sprintf("%s: %s", names(modes), modes), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  joint <- identical(modes[[1L]], "joint")
  # Unset = the simulator's default 200 (07b / simulations.R / smoothing_core.R; run_design_row()).
  S <- unique(vapply(
    arm_cfgs,
    function(cf) as.integer(cf$simulator$max_draws %||% 200L),
    integer(1)
  ))
  if (length(S) != 1L || is.na(S)) {
    stop(
      sprintf(
        "[09] af_posterior: unit %s arms disagree on simulator.max_draws (unset = 200).",
        unit
      ),
      call. = FALSE
    )
  }
  if (
    joint &&
      !is.na(as.integer(n_smoothing %||% NA)) &&
      S != as.integer(n_smoothing)
  ) {
    stop(
      sprintf(
        paste0(
          "[09] af_posterior: unit %s runs param_draw_mode joint with max_draws %d but %d ",
          "smoothing replicates: the joint design has one draw per smoothing replicate."
        ),
        unit,
        S,
        as.integer(n_smoothing)
      ),
      call. = FALSE
    )
  }
  # Uniform handoff starts replicate r from smoothing replicate
  # ((r - 1) mod n_smoothing) + 1: a max_draws sharing a factor with
  # n_smoothing ties each draw to a few smoothing states (200 vs 100: one), so
  # draw and handoff state are confounded. "shuffled" handoff (a seeded
  # permutation per block of n_smoothing ids) does not, so only arms on
  # uniform handoff are checked.
  handoff <- vapply(
    arm_cfgs,
    function(cf) {
      as.character(cf$simulator$forecast[["handoff_weighting"]] %||% "uniform")
    },
    character(1)
  )
  n_sm <- as.integer(n_smoothing %||% NA)
  # joint: draw and state are paired by design (S = n_sm), not confounded.
  if (!joint && !is.na(n_sm) && any(handoff == "uniform")) {
    gcd <- function(a, b) if (b == 0L) a else gcd(b, a %% b)
    if (gcd(S, n_sm) != 1L) {
      warning(
        sprintf(
          paste0(
            "[09] af_posterior: unit %s: max_draws %d and the %d smoothing replicates share ",
            "a factor (gcd %d): each draw starts from only %d smoothing states."
          ),
          unit,
          S,
          n_sm,
          gcd(S, n_sm),
          n_sm %/% gcd(S, n_sm)
        ),
        call. = FALSE
      )
    }
  }
  for (a in names(arm_draws)) {
    d <- arm_draws[[a]]
    if (any(d$draw != ((d$replicate_id - 1L) %% S) + 1L)) {
      stop(
        sprintf(
          "[09] af_posterior: unit %s arm %s: draws are not the cycle index ((r - 1) mod %d) + 1.",
          unit,
          a,
          S
        ),
        call. = FALSE
      )
    }
  }
  ref <- arm_draws[[1L]]
  for (a in names(arm_draws)[-1L]) {
    j <- dplyr::inner_join(
      ref,
      arm_draws[[a]],
      by = "replicate_id",
      suffix = c("", ".b")
    )
    if (nrow(j) != nrow(ref) || any(j$draw != j$draw.b)) {
      stop(
        sprintf(
          "[09] af_posterior: unit %s arm %s does not share the factual arm's replicate -> draw map.",
          unit,
          a
        ),
        call. = FALSE
      )
    }
  }
  structure(ref, max_draws = S)
}

#' The replicate -> draw map of every unit holding one of `groups` (manifest
#' `man`, its arms for `sets`), checked by af_posterior_check_design().
#' tibble(unit, replicate_id, draw); attr "max_draws" (named by unit).
af_posterior_draws <- function(ctx, man, groups, sets) {
  arms <- results_unit_arm_configs(man, groups, sets)
  S <- integer(0)
  out <- dplyr::bind_rows(lapply(names(arms), function(un) {
    paths <- arms[[un]]
    names(paths) <- c("factual", sets)
    cfgs <- lapply(paths, yaml::read_yaml)
    d <- lapply(cfgs, function(cf) {
      af_posterior_run_draws(ctx$SIM_ROOT, cf$simulator$run_tag)
    })
    meta <- file.path(
      cfgs[[1L]]$simulator$filter$smoothing_input %||% "",
      "meta.yml"
    )
    n_sm <- if (file.exists(meta)) yaml::read_yaml(meta)$nrep else NULL
    chk <- af_posterior_check_design(un, cfgs, d, n_smoothing = n_sm)
    S[[un]] <<- attr(chk, "max_draws")
    dplyr::mutate(chk, unit = un, .before = 1L)
  }))
  structure(out, max_draws = S)
}

#' Paired factual / counterfactual window totals per replicate and country of
#' one analog set, with the replicate's draw: tibble(cell, enso_group, unit,
#' country, replicate_id, draw, f, cf). `guard`: a replicate explosive in either
#' arm is dropped. `guard_scope` "unit" (the default): from every group and country of
#' its simulation unit, so the groups
#' of one unit keep a common replicate set and the replicate-level Africa total
#' is never missing a group; "group": from every country of
#' its group only.
af_posterior_pairs <- function(
  reps,
  explosive,
  draws,
  aset,
  metric,
  guard = TRUE,
  guard_scope = c("unit", "group")
) {
  guard_scope <- match.arg(guard_scope)
  r <- dplyr::filter(reps, run == "forecast" | analog_set == aset)
  if (guard) {
    ex <- dplyr::filter(explosive, run == "forecast" | analog_set == aset)
    if (guard_scope == "unit") {
      ex <- dplyr::inner_join(
        ex,
        dplyr::distinct(reps, cell, enso_group, unit),
        by = c("cell", "enso_group"),
        relationship = "many-to-many"
      ) %>%
        dplyr::distinct(cell, unit, run, analog_set, replicate_id)
      r <- dplyr::anti_join(
        r,
        ex,
        by = c("cell", "unit", "run", "analog_set", "replicate_id")
      )
    } else {
      r <- dplyr::anti_join(
        r,
        ex,
        by = c("cell", "enso_group", "run", "analog_set", "replicate_id")
      )
    }
  }
  keys <- c("cell", "enso_group", "unit", "country", "replicate_id")
  arm <- function(rn, col) {
    dplyr::filter(r, run == rn) %>%
      dplyr::select(dplyr::all_of(keys), !!col := dplyr::all_of(metric))
  }
  p <- dplyr::inner_join(
    arm("forecast", "f"),
    arm("counterfactual", "cf"),
    by = keys
  ) %>%
    dplyr::left_join(draws, by = c("cell", "unit", "replicate_id"))
  if (anyNA(p$draw)) {
    stop(
      "[09] af_posterior: replicates without a posterior draw.",
      call. = FALSE
    )
  }
  p
}

#' Replicates per (cell, group, unit, draw) after the guard, against the
#' design: stops when a draw of the design is missing from a group's unit or has
#' fewer than `min_reps` replicates. Returns the smallest count.
af_posterior_check_reps <- function(pairs, draws, min_reps = AF_POST_MIN_REPS) {
  n <- pairs %>%
    dplyr::distinct(cell, enso_group, unit, draw, replicate_id) %>%
    dplyr::count(cell, enso_group, unit, draw, name = "n")
  design <- dplyr::count(
    dplyr::distinct(draws, cell, unit, draw),
    cell,
    unit,
    name = "n_design"
  )
  have <- dplyr::count(n, cell, enso_group, unit, name = "n_have") %>%
    dplyr::left_join(design, by = c("cell", "unit"))
  bad <- dplyr::filter(have, n_have < n_design)
  if (nrow(bad)) {
    stop(
      sprintf(
        "[09] af_posterior: %s / %s / unit %s has %d of the design's %d draws after the explosive drop.",
        bad$cell[1],
        bad$enso_group[1],
        bad$unit[1],
        bad$n_have[1],
        bad$n_design[1]
      ),
      call. = FALSE
    )
  }
  low <- dplyr::filter(n, n < min_reps)
  if (nrow(low)) {
    stop(
      sprintf(
        paste0(
          "[09] af_posterior: %s / %s / unit %s / draw %d has %d replicates (< %d) after the ",
          "explosive drop: too few to estimate the within-draw noise."
        ),
        low$cell[1],
        low$enso_group[1],
        low$unit[1],
        low$draw[1],
        low$n[1],
        min_reps
      ),
      call. = FALSE
    )
  }
  min(n$n)
}

#' Expected totals per draw of each scope `by` ("enso_group", "country" or
#' "africa"), with the within-draw (co)variance terms of the unit means summed
#' over units (a = sum var(cf) / n, b = sum cov(cf, f) / n, c = sum var(f) / n);
#' only draws present in every unit of the scope.
#' "africa": one scope over every row, summed per
#' (unit, replicate) BEFORE the per-draw moments, so the groups simulated
#' together in a unit keep their within-replicate dependence in both the mean and
#' the noise terms; across units the simulations are independent and draw d is
#' the same index in every unit (each group at its own posterior draw d), so
#' summing the unit terms at d is exact. enso_group is NA on those rows.
#' tibble(cell, scope, enso_group, draw, f, cf, a, b, c, n_min).
af_posterior_draw_table <- function(pairs, by = "enso_group") {
  if (identical(by, "africa")) {
    pairs <- dplyr::mutate(pairs, africa = "Africa", enso_group = NA_character_)
  }
  # data.table for the two replicate-level passes (~10^6 / ~10^4 groups at
  # 10^4 replicates); the same R functions per group, so identical output.
  per_rep <- .sum_by_keys(
    dplyr::mutate(pairs, scope = .data[[by]]),
    c("cell", "scope", "enso_group", "unit", "replicate_id", "draw"),
    c(f = "f", cf = "cf")
  )
  n_units <- dplyr::count(
    dplyr::distinct(per_rep, cell, scope, unit),
    cell,
    scope,
    name = "n_units"
  )
  v0 <- function(x) if (length(x) > 1L) stats::var(x) else 0
  c0 <- function(x, y) if (length(x) > 1L) stats::cov(x, y) else 0
  by_ud <- c("cell", "scope", "enso_group", "unit", "draw")
  ud <- data.table::as.data.table(per_rep)[,
    {
      n <- .N
      list(
        n = n,
        a = v0(cf) / n,
        b = c0(cf, f) / n,
        c = v0(f) / n,
        f = mean(f),
        cf = mean(cf)
      )
    },
    keyby = by_ud
  ]
  data.table::setkey(ud, NULL)
  if (anyNA(ud[, by_ud, with = FALSE])) {
    data.table::setorderv(ud, by_ud, na.last = TRUE)
  }
  tibble::as_tibble(as.data.frame(ud)) %>%
    dplyr::group_by(cell, scope, enso_group, draw) %>%
    dplyr::summarise(
      k = dplyr::n(),
      n_min = min(n),
      f = sum(f),
      cf = sum(cf),
      a = sum(a),
      b = sum(b),
      c = sum(c),
      .groups = "drop"
    ) %>%
    dplyr::left_join(n_units, by = c("cell", "scope")) %>%
    dplyr::filter(k == n_units) %>%
    dplyr::select(-k, -n_units)
}

#' AF_d and AN_d with their delta-method noise variances from the per-draw
#' terms: data.frame(af, af_noise, an, an_noise).
.afp_stats <- function(f, cf, a, b, c) {
  r <- cf / f
  data.frame(
    af = ifelse(f > 0, 1 - r, NA_real_),
    af_noise = ifelse(f > 0, (a - 2 * r * b + r^2 * c) / f^2, NA_real_),
    an = f - cf,
    an_noise = a + c - 2 * b
  )
}

#' Quantiles of one per-draw statistic, raw and shrunk towards `center` by the
#' noise share mean(noise) / var(x): one-row tibble with prefix `p` (af / an).
#'
#' `noise_share_max` (results.thresholds.af_noise_share_max): a statistic whose noise
#' share is NOT below it (or not finite)
#' has an unidentified posterior spread: `<p>_identified` is FALSE and its
#' adjusted quantiles are NA (the pooled value and the raw quantiles stay), so no
#' reader can take the collapsed interval for a certain one. Inf = the rule
#' before: never flagged, a share >= 1 collapses every draw to the pooled value.
#' `<p>_var_between` = max(0, var(x) - mean(noise)), the between-draw variance
#' the adjustment keeps.
.afp_summary <- function(x, noise, center, p, noise_share_max = Inf) {
  q <- function(v, pr) unname(stats::quantile(v, pr, na.rm = TRUE))
  PR <- c(lo = .025, q25 = .25, med = .5, q75 = .75, hi = .975)
  v_raw <- stats::var(x, na.rm = TRUE)
  share <- mean(noise, na.rm = TRUE) / v_raw
  s <- if (is.finite(share)) sqrt(max(0, 1 - share)) else NA_real_
  adj <- center + s * (x - center)
  ident <- if (is.finite(noise_share_max)) {
    is.finite(share) && share < noise_share_max
  } else {
    is.finite(share)
  }
  adj_q <- if (ident || !is.finite(noise_share_max)) {
    vapply(PR, function(pr) q(adj, pr), numeric(1))
  } else {
    stats::setNames(rep(NA_real_, length(PR)), names(PR))
  }
  # Diagnostic: the share on robust scales (median noise / IQR variance), less
  # sensitive to the tail draws of a scope with near-zero F_d.
  share_iqr <- stats::median(noise, na.rm = TRUE) /
    (stats::IQR(x, na.rm = TRUE) / 1.349)^2
  out <- c(
    stats::setNames(center, paste0(p, "_pool")),
    stats::setNames(
      vapply(PR, function(pr) q(x, pr), numeric(1)),
      paste0(p, "_raw_", names(PR))
    ),
    stats::setNames(
      c(share, s, share_iqr),
      paste0(p, c("_noise_share", "_shrink", "_noise_share_iqr"))
    ),
    stats::setNames(adj_q, paste0(p, "_adj_", names(PR)))
  )
  out <- tibble::as_tibble(as.list(out))
  if (is.finite(noise_share_max)) {
    out[[paste0(p, "_identified")]] <- ident
    out[[paste0(p, "_var_between")]] <- max(
      0,
      v_raw - mean(noise, na.rm = TRUE)
    )
  }
  out
}

#' af_posterior_tables
#'
#' Pure core of af_posterior / event_af_posterior: the per-draw AF / AN of
#' every cell x analog set x explosive guard x metric x scope (Africa, ENSO
#' groups, countries) and their summaries.
#'
#' @param reps per-replicate window totals (af_posterior_window_reps()), any
#'   number of cells
#' @param explosive the explosive keys of the same cells
#' @param draws tibble(cell, unit, replicate_id, draw): the design per cell
#' @param sets analog sets; @param group_short group -> short label
#' @param guards explosive guards to derive (TRUE = dropped, 09's headline)
#' @param africa "joint" (the default): Africa per
#'   draw is the replicate-level sum of every unit (af_posterior_draw_table(by =
#'   "africa")), S draw rows; "independent": one draw
#'   per group at random, `n_mc` combinations under `seed`
#' @param guard_scope explosive guard scope, af_posterior_pairs()
#' @param noise_share_max unidentified threshold, .afp_summary() (Inf = never)
#' @return list(scope, draws, n_min): scope = one row per cell x set x guard x
#'   metric x scope with af_* / an_* pool, raw and adjusted quantiles, noise
#'   share, shrink and robust share (diagnostic), the identified flag and
#'   between-draw variance (finite noise_share_max), n_draws, n_mc (Africa
#'   under "independent" only), n_min; draws = AF_d / AN_d per draw (Africa
#'   "independent": the combinations, draw = combination index)
af_posterior_tables <- function(
  reps,
  explosive,
  draws,
  sets,
  group_short,
  guards = c(TRUE, FALSE),
  min_reps = AF_POST_MIN_REPS,
  n_mc = AF_POST_MC,
  seed = AF_POST_SEED,
  africa = c("joint", "independent"),
  guard_scope = c("unit", "group"),
  noise_share_max = Inf
) {
  africa <- match.arg(africa)
  guard_scope <- match.arg(guard_scope)
  scope_rows <- list()
  draw_rows <- list()
  n_min_all <- Inf
  for (aset in sets) {
    for (guard in guards) {
      for (metric in c("cases", "onsets")) {
        pairs <- af_posterior_pairs(
          reps,
          explosive,
          draws,
          aset,
          metric,
          guard,
          guard_scope = guard_scope
        )
        n_min_all <- min(
          n_min_all,
          af_posterior_check_reps(pairs, draws, min_reps)
        )
        tg <- af_posterior_draw_table(pairs, "enso_group")
        tc <- af_posterior_draw_table(pairs, "country")
        ta <- if (africa == "joint") {
          af_posterior_draw_table(pairs, "africa")
        } else {
          NULL
        }
        pool <- function(t) 1 - sum(t$cf) / sum(t$f)
        one <- function(
          cell,
          level,
          key,
          grp,
          t,
          st,
          af_center,
          an_center,
          n_draws = sum(is.finite(st$af)),
          n_mc = NA_integer_
        ) {
          lab <- dplyr::case_when(
            level == "continent" ~ "Africa",
            level == "region" ~ unname(group_short[key]),
            TRUE ~ key
          )
          scope_rows[[length(scope_rows) + 1L]] <<- dplyr::bind_cols(
            tibble::tibble(
              cell = cell,
              analog_set = aset,
              explosive_guard = guard,
              metric = metric,
              level = level,
              key = key,
              scope = lab,
              enso_group = grp,
              region = if (is.na(grp)) {
                NA_character_
              } else {
                unname(group_short[grp])
              },
              country = if (level == "country") key else NA_character_,
              n_draws = n_draws,
              n_mc = n_mc,
              n_min = if (is.null(t$n_min)) NA_integer_ else min(t$n_min)
            ),
            .afp_summary(st$af, st$af_noise, af_center, "af", noise_share_max),
            .afp_summary(st$an, st$an_noise, an_center, "an", noise_share_max)
          )
          draw_rows[[length(draw_rows) + 1L]] <<- tibble::tibble(
            cell = cell,
            analog_set = aset,
            explosive_guard = guard,
            metric = metric,
            level = level,
            key = key,
            draw = if (is.null(t$draw)) seq_len(nrow(st)) else t$draw,
            af = st$af,
            an = st$an
          )
        }
        for (cell in unique(tg$cell)) {
          t_c <- dplyr::filter(tg, cell == !!cell)
          for (g in unique(t_c$scope)) {
            t <- dplyr::filter(t_c, scope == g)
            one(
              cell,
              "region",
              g,
              g,
              t,
              .afp_stats(t$f, t$cf, t$a, t$b, t$c),
              pool(t),
              mean(t$f - t$cf)
            )
          }
          if (africa == "joint") {
            # Africa: the replicate-level total of every unit per draw.
            t <- dplyr::filter(ta, cell == !!cell)
            one(
              cell,
              "continent",
              "Africa",
              NA_character_,
              t,
              .afp_stats(t$f, t$cf, t$a, t$b, t$c),
              pool(t),
              mean(t$f - t$cf)
            )
          } else {
            # Africa "independent" (legacy): one draw per group, n_mc
            # combinations, as if the groups' outcomes were independent.
            pick <- .with_seed(
              seed,
              lapply(split(t_c, t_c$scope), function(t_g) {
                j <- sample.int(nrow(t_g), n_mc, replace = TRUE)
                as.matrix(t_g[j, c("f", "cf", "a", "b", "c")])
              })
            )
            tot <- Reduce(`+`, pick)
            gm <- dplyr::summarise(
              dplyr::group_by(t_c, scope),
              f = mean(f),
              cf = mean(cf),
              .groups = "drop"
            )
            one(
              cell,
              "continent",
              "Africa",
              NA_character_,
              list(n_min = t_c$n_min),
              .afp_stats(
                tot[, "f"],
                tot[, "cf"],
                tot[, "a"],
                tot[, "b"],
                tot[, "c"]
              ),
              pool(gm),
              sum(gm$f - gm$cf),
              n_draws = min(table(t_c$scope)),
              n_mc = n_mc
            )
          }
          t_cc <- dplyr::filter(tc, cell == !!cell)
          for (k in unique(t_cc$scope)) {
            t <- dplyr::filter(t_cc, scope == k)
            one(
              cell,
              "country",
              k,
              t$enso_group[1],
              t,
              .afp_stats(t$f, t$cf, t$a, t$b, t$c),
              pool(t),
              mean(t$f - t$cf)
            )
          }
        }
      }
    }
  }
  list(
    scope = dplyr::bind_rows(scope_rows),
    draws = dplyr::bind_rows(draw_rows),
    n_min = n_min_all
  )
}

derive_af_posterior <- function(ctx, sets) {
  key <- paste("forecast", paste(sort(sets), collapse = "-"))
  st <- ctx$AF_POST_STASH[[key]]
  if (is.null(st)) {
    inp <- attrib_arm_inputs(ctx, sets, "af-posterior")
    st <- af_posterior_stash_build(
      inp,
      attrib_cells(
        ctx,
        inp$attrib_start,
        display_factual_start(ctx, inp$groups, sets)
      )
    )
  } else {
    message(
      "[09] af_posterior: reusing the window totals of this job's attribution load"
    )
  }
  d <- af_posterior_draws(ctx, ctx$MAN, st$groups, sets)
  draws <- dplyr::bind_rows(lapply(names(st$cells), function(k) {
    dplyr::mutate(d, cell = k, .before = 1L)
  }))
  message(sprintf(
    "[09] af_posterior: groups=%s | sets=%s | cells=%s | max_draws=%s",
    paste(st$groups, collapse = ","),
    paste(sets, collapse = ","),
    paste(names(st$cells), collapse = ","),
    paste(unique(attr(d, "max_draws")), collapse = ",")
  ))
  tabs <- af_posterior_tables(
    st$reps,
    st$explosive,
    draws,
    sets,
    ctx$GROUP_SHORT_RUN,
    guards = .cell_guards(st$cell_table),
    min_reps = ctx$AF_MIN_REPS,
    guard_scope = ctx$AF_GUARD_SCOPE,
    noise_share_max = ctx$AF_NOISE_SHARE_MAX
  )
  c(
    tabs,
    list(
      windows = st$cells,
      analog_sets = sets,
      groups = st$groups,
      max_draws = attr(d, "max_draws"),
      min_reps = ctx$AF_MIN_REPS,
      africa = "joint",
      guard_scope = ctx$AF_GUARD_SCOPE,
      noise_share_max = ctx$AF_NOISE_SHARE_MAX,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
      cells = cells_payload(st$cell_table),
      africa_window = africa_window_payload(
        .africa_win(ctx, st$groups),
        st$cells[["attribution"]]
      )
    )
  ) # af_posterior v3
}

#' Freshness globs of the runs' draw maps: every arm's partition assignment.qs
#' (a rerun with the same configs rewrites them).
.afp_assignment_globs <- function(ctx, cfgs) {
  lapply(unname(cfgs), function(cf) {
    list(
      dir = file.path(
        ctx$SIM_ROOT,
        yaml::read_yaml(cf)$simulator$run_tag,
        "parts"
      ),
      pattern = "^assignment\\.qs$",
      recursive = TRUE
    )
  })
}

ensure_af_posterior <- function(ctx, sets) {
  ctx <- results_attrib_ctx(ctx)
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- results_run_configs(man, groups, sets)
  inputs <- results_inputs(
    files = c(
      results_attrib_manifest_files(ctx),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = c(
      if (results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched) {
        smoothing_parts_glob(ctx)
      } else {
        list()
      },
      .afp_assignment_globs(ctx, cfgs)
    )
  )
  ensure_artefact(
    ctx,
    "af_posterior",
    inputs,
    function() derive_af_posterior(ctx, sets),
    sets = sets
  )
}

derive_event_af_posterior <- function(ctx) {
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  sets <- as.character(em$analog_set)
  base <- NULL
  parts <- list()
  meta <- list()
  parts_sw <- list()
  for (ev in names(em$events)) {
    view <- event_manifest_view(em, ev)
    swt <- event_sub_windows_table(em, ev)
    key <- paste("event", ev, paste(sort(sets), collapse = "-"))
    st <- ctx$AF_POST_STASH[[key]]
    if (is.null(st)) {
      ctx_e <- .results_ctx_with(
        ctx,
        MAN = view,
        MANIFEST_PATH = ctx$EVENT_MANIFEST_PATH,
        AFRICA_WIN = NULL
      )
      base <- base %||% attrib_base_inputs(ctx)
      inp <- attrib_arm_inputs(
        ctx_e,
        sets,
        "event-af-posterior",
        factual_source = "simulation",
        base = base
      )
      st <- af_posterior_stash_build(
        inp,
        event_cells_all(ctx, inp$attrib_start, ev, swt)
      )
    } else {
      message(sprintf(
        "[09] event_af_posterior: %s reuses report_event_attrib's load",
        ev
      ))
    }
    d <- af_posterior_draws(ctx, view, st$groups, sets)
    message(sprintf(
      "[09] event_af_posterior: %s groups=%s max_draws=%s",
      ev,
      paste(st$groups, collapse = ","),
      paste(unique(attr(d, "max_draws")), collapse = ",")
    ))
    # One cell per window: the event and, schema-3 manifests,
    # each per-winter sub-window (the draws repeated per cell, as af_posterior).
    cells_ev <- if (is.null(swt)) ev else names(st$cells)
    d_cells <- if (is.null(swt)) {
      dplyr::mutate(d, cell = ev, .before = 1L)
    } else {
      dplyr::bind_rows(lapply(cells_ev, function(k) {
        dplyr::mutate(d, cell = k, .before = 1L)
      }))
    }
    tabs <- af_posterior_tables(
      st$reps,
      st$explosive,
      d_cells,
      sets,
      ctx$GROUP_SHORT_RUN,
      guards = .cell_guards(st$cell_table),
      min_reps = ctx$AF_MIN_REPS,
      guard_scope = ctx$AF_GUARD_SCOPE,
      noise_share_max = ctx$AF_NOISE_SHARE_MAX
    )
    if (!is.null(swt)) {
      parts_sw[[ev]] <- list(
        scope = dplyr::rename(
          dplyr::mutate(
            dplyr::filter(tabs$scope, cell != ev),
            event = ev,
            .before = 1L
          ),
          sub_window = cell
        ),
        draws = dplyr::rename(
          dplyr::mutate(
            dplyr::filter(tabs$draws, cell != ev),
            event = ev,
            .before = 1L
          ),
          sub_window = cell
        )
      )
      tabs$scope <- dplyr::filter(tabs$scope, cell == ev)
      tabs$draws <- dplyr::filter(tabs$draws, cell == ev)
    }
    parts[[ev]] <- tabs
    e <- em$events[[ev]]
    meta[[ev]] <- tibble::tibble(
      event = ev,
      phase = e$phase,
      label = e$label,
      launch = as.Date(e$launch),
      peak = as.Date(e$peak),
      window_end = as.Date(e$africa_max_horizon),
      max_draws = max(attr(d, "max_draws"))
    )
  }
  rn <- function(x) dplyr::rename(x, event = cell)
  res <- list(
    scope = rn(dplyr::bind_rows(lapply(parts, `[[`, "scope"))),
    draws = rn(dplyr::bind_rows(lapply(parts, `[[`, "draws"))),
    n_min = min(vapply(parts, `[[`, numeric(1), "n_min")),
    events = dplyr::bind_rows(meta),
    analog_set = sets,
    factual_source = "simulation",
    min_reps = ctx$AF_MIN_REPS,
    africa = "joint",
    guard_scope = ctx$AF_GUARD_SCOPE,
    noise_share_max = ctx$AF_NOISE_SHARE_MAX,
    explosive_cap = ctx$EXPLOSIVE_CAP,
    explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
    explosive_modes = .ctx_explosive_modes(ctx),
    event_manifest = ctx$EVENT_MANIFEST_PATH
  )
  # v3: the per-winter sub-window cells of schema-3 events
  # (n_min is then over every cell); absent otherwise, so a schema-2 manifest
  # derives the v2 payload exactly.
  if (length(parts_sw)) {
    res$sub_windows <- list(
      scope = dplyr::bind_rows(lapply(parts_sw, `[[`, "scope")),
      draws = dplyr::bind_rows(lapply(parts_sw, `[[`, "draws"))
    )
  }
  res # event_af_posterior v4
}

ensure_event_af_posterior <- function(ctx) {
  if (!nzchar(ctx$EVENT_MANIFEST_PATH %||% "")) {
    message(
      "[09] event_af_posterior: results.event_manifest is empty; skipped."
    )
    return(invisible(NULL))
  }
  em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
  cfgs <- unlist(lapply(em$events, function(e) {
    unlist(lapply(e$units, function(u) {
      c(u$factual_config, vapply(u$counterfactual_sets, `[[`, "", "config"))
    }))
  }))
  inputs <- results_inputs(
    files = c(
      ctx$EVENT_MANIFEST_PATH,
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = .afp_assignment_globs(ctx, cfgs)
  )
  ensure_artefact(ctx, "event_af_posterior", inputs, function() {
    derive_event_af_posterior(ctx)
  })
}

# ===========================================================================
# attrib_exposure_posterior: the posterior interval of the Childs
# population of report_attrib_exposure (the population of the ADM1s whose mean
# forecast is at least (1 + m) x the mean counterfactual), from the fixed-draw
# runs of af_posterior. Per ADM1 and draw d, F_d / CF_d are the means over the
# replicates that ran d, AF_d = 1 - CF_d / F_d with its delta-method noise
# variance (.afp_stats()), and the Childs rule F >= (1 + m) CF is
# AF >= m / (1 + m). As in af_posterior, AF_d is shrunk towards the ADM1's
# pooled AF by sqrt(1 - noise share) before the rule is applied (variant "adj",
# the headline); "raw" applies the rule to the unshrunk per-draw means, so it
# counts the epidemic noise of ~50 replicates as posterior spread (an upper
# bound on the width). An ADM1 whose noise share is >= 1 keeps its pooled class
# on every draw: its posterior spread is not identified. A draw with F_d = 0 is
# no hit unless CF_d = 0 too (then the pooled class). Population per draw: the
# population of the ADM1s that pass; every scope from AF_MC combinations of one
# draw per group (independent fits; a region's quantiles are those of its own
# draws, resampled). Same windows (attribution / post_fork), explosive guards
# and magnitudes as report_attrib_exposure; pop_point is its Childs population.
# Gate: tests/test_results_attrib_exposure_posterior.R.
# ===========================================================================

#' Per-ADM1 x draw means of a paired table (adm1_unit, draw, f, cf) with the
#' within-draw (co)variance terms of the means (a = var(cf) / n, b = cov / n,
#' c = var(f) / n), as af_posterior_draw_table():
#' tibble(adm1_unit, draw, n, f, cf, a, b, c).
attrib_adm1_draw_table <- function(paired) {
  v0 <- function(x) if (length(x) > 1L) stats::var(x) else 0
  c0 <- function(x, y) if (length(x) > 1L) stats::cov(x, y) else 0
  ud <- data.table::as.data.table(paired)[,
    {
      n <- .N
      list(
        n = n,
        f = mean(f),
        cf = mean(cf),
        a = v0(cf) / n,
        b = c0(cf, f) / n,
        c = v0(f) / n
      )
    },
    keyby = c("adm1_unit", "draw")
  ]
  tibble::as_tibble(as.data.frame(ud))
}

#' Childs hits per ADM1 and draw of an attrib_adm1_draw_table() at magnitude m:
#' the table with af, af_noise, af_pool (1 - sum CF_d / sum F_d), noise_share,
#' shrink, pool_hit and the hits `adj` (shrunk AF_d >= m / (1 + m)) and `raw`
#' (childs_hit(F_d, CF_d, m)).
#'
#' `noise_share_max`: an ADM1 whose noise share
#' is not below it (or not finite) is unidentified (`identified` FALSE) and its
#' `adj` hit is NA on every draw instead of its pooled class; Inf = no
#' check (never NA: a share >= 1 keeps the pooled class).
childs_draw_hits <- function(dt, m, noise_share_max = Inf) {
  st <- .afp_stats(dt$f, dt$cf, dt$a, dt$b, dt$c)
  x <- dplyr::mutate(dt, af = st$af, af_noise = st$af_noise) %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::mutate(
      F_pool = sum(f),
      CF_pool = sum(cf),
      noise_share = mean(af_noise, na.rm = TRUE) / stats::var(af, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      af_pool = ifelse(F_pool > 0, 1 - CF_pool / F_pool, NA_real_),
      shrink = ifelse(
        is.finite(noise_share),
        sqrt(pmax(0, 1 - noise_share)),
        0
      ),
      pool_hit = childs_hit(F_pool, CF_pool, m),
      identified = is.finite(noise_share) & noise_share < noise_share_max
    )
  x$adj <- ifelse(
    x$shrink == 0,
    x$pool_hit,
    ifelse(
      is.na(x$af),
      x$cf == 0 & x$pool_hit,
      x$af_pool + x$shrink * (x$af - x$af_pool) >= m / (1 + m)
    )
  )
  if (is.finite(noise_share_max)) {
    x$adj[!x$identified] <- NA
  }
  x$raw <- childs_hit(x$f, x$cf, m)
  dplyr::select(x, -F_pool, -CF_pool)
}

#' The pure core of attrib_exposure_posterior (no IO) for ONE window and one
#' explosive-replicate cap.
#'
#' @param adm1_month,country_month,adm1_pop,group_short,sets,window_start,mags,cap
#'   as attrib_exposure_tables()
#' @param draws tibble(unit, replicate_id, draw): the fixed-draw design
#'   (af_posterior_draws()); the ADM1's unit is read from adm1_month (`unit`, or
#'   the group for a v1 manifest)
#' @return list(exposure: one row per metric x magnitude x variant (adj | raw)
#'   x scope x analog set with pop_mean / pop_lo / pop_med / pop_hi (+ pct_*),
#'   n_adm1_mean and the point n_adm1_point / pop_point; adm1: per ADM1 x
#'   metric x set x magnitude, af_pool, noise_share, shrink, the point hit and
#'   the share of draws that pass (p_adj, p_raw); n_min: smallest replicate
#'   count per (unit, draw))
attrib_exposure_posterior_tables <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  window_start,
  mags,
  draws,
  cap,
  min_reps = AF_POST_MIN_REPS,
  n_mc = AF_POST_MC,
  seed = AF_POST_SEED,
  africa = c("joint", "independent"),
  guard_scope = c("unit", "group"),
  noise_share_max = Inf
) {
  africa <- match.arg(africa)
  guard_scope <- match.arg(guard_scope)
  adm1_month <- dplyr::mutate(adm1_month, run = as.character(run))
  country_month <- dplyr::mutate(country_month, run = as.character(run))
  if (!"unit" %in% names(adm1_month)) {
    adm1_month$unit <- adm1_month$enso_group
  }
  expl <- explosive_replicate_keys(country_month, window_start, cap)
  kept <- if (guard_scope == "unit") {
    # Issue 7: a replicate explosive in any group of a unit leaves the whole unit.
    ex_u <- dplyr::inner_join(
      expl,
      dplyr::distinct(adm1_month, enso_group, unit),
      by = "enso_group",
      relationship = "many-to-many"
    ) %>%
      dplyr::distinct(unit, run, analog_set, replicate_id)
    .anti_join_rows(
      adm1_month,
      ex_u,
      by = c("unit", "run", "analog_set", "replicate_id")
    )
  } else {
    .anti_join_rows(
      adm1_month,
      expl,
      by = c("enso_group", "run", "analog_set", "replicate_id")
    )
  }
  cum_all <- pair_adm1_cum_table(kept)

  dr <- attrib_adm1_universe(adm1_month, adm1_pop, group_short)
  groups <- unique(dr$universe$enso_group)
  rows_g <- lapply(stats::setNames(groups, groups), function(g) {
    which(dr$universe$enso_group == g)
  })
  a2u <- dplyr::distinct(
    dplyr::filter(adm1_month, run == "forecast"),
    adm1_unit,
    enso_group,
    unit
  )
  if (anyDuplicated(a2u$adm1_unit)) {
    stop(
      "[09] attrib_exposure_posterior: ADM1 unit(s) in more than one simulation unit.",
      call. = FALSE
    )
  }
  # "independent" (legacy): one draw per group for each of the n_mc
  # combinations, shared by every metric, set and magnitude (and, at a fixed
  # seed, by every window and guard). "joint": the draw axis is the draw index
  # itself, every ADM1 at the same d (the draws present in every unit).
  dg <- dplyr::left_join(
    dplyr::distinct(a2u, enso_group, unit),
    dplyr::distinct(draws, unit, draw),
    by = "unit",
    relationship = "many-to-many"
  )
  d_g <- lapply(stats::setNames(groups, groups), function(g) {
    sort(unique(dg$draw[dg$enso_group == g]))
  })
  if (africa == "independent") {
    pick <- .with_seed(
      seed,
      lapply(d_g, function(v) sample.int(length(v), n_mc, replace = TRUE))
    )
    n_ax <- n_mc
  } else {
    d_all <- sort(Reduce(intersect, split(dg$draw, dg$unit)))
    n_ax <- length(d_all)
  }
  flag <- is.finite(noise_share_max)

  ADM1_METRIC <- c(cases = "cases", onsets = "outbreaks")
  S <- length(sets)
  U <- length(dr$U)
  exposure <- list()
  adm1_rows <- list()
  n_min <- Inf
  for (metric in names(ADM1_METRIC)) {
    tabs <- lapply(stats::setNames(sets, sets), function(aset) {
      p <- pair_adm1_cum(kept, ADM1_METRIC[[metric]], aset, cum = cum_all) %>%
        dplyr::left_join(a2u, by = "adm1_unit") %>%
        dplyr::left_join(draws, by = c("unit", "replicate_id"))
      if (anyNA(p$draw)) {
        stop(
          "[09] attrib_exposure_posterior: replicates without a posterior draw.",
          call. = FALSE
        )
      }
      chk <- dplyr::mutate(
        dplyr::distinct(p, enso_group, unit, replicate_id, draw),
        cell = "w"
      )
      n_min <<- min(
        n_min,
        af_posterior_check_reps(chk, dplyr::mutate(draws, cell = "w"), min_reps)
      )
      list(
        dt = attrib_adm1_draw_table(p),
        mean = dplyr::summarise(
          dplyr::group_by(p, adm1_unit),
          mean_f = mean(f),
          mean_cf = mean(cf),
          .groups = "drop"
        )
      )
    })
    for (m in mags) {
      arr <- list(
        adj = array(FALSE, c(U, S, n_ax)),
        raw = array(FALSE, c(U, S, n_ax))
      )
      # Issue 9: unidentified ADM1s leave the adj population (their adj hit is
      # NA) and are reported apart: their population (unidentified) and the
      # part of it whose pooled estimate passes the rule (unidentified_pass).
      if (flag) {
        arr$unidentified <- arr$unidentified_pass <- array(FALSE, c(U, S, n_ax))
      }
      point <- matrix(FALSE, U, S)
      pt_un <- pt_unp <- matrix(FALSE, U, S)
      for (s in seq_len(S)) {
        h <- childs_draw_hits(tabs[[s]]$dt, m, noise_share_max)
        mn <- tabs[[s]]$mean
        point[match(mn$adm1_unit, dr$U), s] <- childs_hit(
          mn$mean_f,
          mn$mean_cf,
          m
        )
        if (flag) {
          h$unidentified <- !h$identified
          h$unidentified_pass <- !h$identified & h$pool_hit %in% TRUE
          hu <- dplyr::distinct(h, adm1_unit, unidentified, unidentified_pass)
          iu <- match(hu$adm1_unit, dr$U)
          pt_un[iu[!is.na(iu)], s] <- hu$unidentified[!is.na(iu)]
          pt_unp[iu[!is.na(iu)], s] <- hu$unidentified_pass[!is.na(iu)]
        }
        # The population masks count an unidentified ADM1 (adj NA) as no hit.
        h_m <- h
        if (flag) {
          h_m$adj[is.na(h_m$adj)] <- FALSE
        }
        if (africa == "joint") {
          i <- match(h_m$adm1_unit, dr$U)
          j <- match(h_m$draw, d_all)
          keep <- !is.na(i) & !is.na(j)
          for (v in names(arr)) {
            arr[[v]][cbind(i[keep], s, j[keep])] <- h_m[[v]][keep]
          }
        } else {
          for (g in groups) {
            i <- match(h_m$adm1_unit, dr$U)
            j <- match(h_m$draw, d_g[[g]])
            keep <- !is.na(i) & dr$universe$enso_group[i] == g & !is.na(j)
            for (v in names(arr)) {
              M <- matrix(FALSE, U, length(d_g[[g]]))
              M[cbind(i[keep], j[keep])] <- h_m[[v]][keep]
              arr[[v]][rows_g[[g]], s, ] <- M[
                rows_g[[g]],
                pick[[g]],
                drop = FALSE
              ]
            }
          }
        }
        adm1_rows[[length(adm1_rows) + 1L]] <- h %>%
          dplyr::group_by(adm1_unit) %>%
          dplyr::summarise(
            n_draws = dplyr::n(),
            af_pool = af_pool[1],
            noise_share = noise_share[1],
            shrink = shrink[1],
            identified = identified[1],
            p_adj = mean(adj),
            p_raw = mean(raw),
            .groups = "drop"
          ) %>%
          {
            if (flag) . else dplyr::select(., -identified)
          } %>%
          dplyr::left_join(
            dplyr::transmute(
              mn,
              adm1_unit,
              mean_f,
              mean_cf,
              childs = childs_hit(mean_f, mean_cf, m)
            ),
            by = "adm1_unit"
          ) %>%
          dplyr::mutate(metric = metric, analog_set = sets[s], magnitude = m)
      }
      pts <- list(adj = point, raw = point)
      # The adj point counts identified ADM1s only, as its interval does, so
      # adj + unidentified_pass = the Childs point population of every ADM1.
      if (flag) {
        pts$adj <- point & !pt_un
        pts$unidentified <- pt_un
        pts$unidentified_pass <- pt_unp
      }
      exposure[[length(exposure) + 1L]] <- mask_pop_rows(
        arr,
        dr,
        point = pts[names(arr)],
        cond_cols = tibble::tibble(analog_set = sets)
      ) %>%
        dplyr::rename(variant = track) %>%
        dplyr::mutate(metric = metric, magnitude = m, .before = 1L)
    }
  }
  adm1 <- dplyr::bind_rows(adm1_rows) %>%
    dplyr::left_join(
      dplyr::select(
        dr$universe,
        adm1_unit = adm1_gid,
        country,
        enso_group,
        region,
        pop
      ),
      by = "adm1_unit"
    ) %>%
    dplyr::relocate(
      adm1_unit,
      country,
      enso_group,
      region,
      metric,
      analog_set,
      magnitude
    )
  list(
    exposure = dplyr::bind_rows(exposure),
    adm1 = adm1,
    n_min = n_min,
    explosive = expl
  )
}

#' attrib_exposure_posterior_tables() over the attribution cells
#' (attrib_cell_table()), rows tagged `window` and `explosive_guard` (as
#' attrib_exposure_windows()).
attrib_exposure_posterior_windows <- function(
  adm1_month,
  country_month,
  adm1_pop,
  group_short,
  sets,
  cells,
  mags,
  draws,
  ...
) {
  out <- list()
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    o <- attrib_exposure_posterior_tables(
      .cell_frame(adm1_month, cl),
      .cell_frame(country_month, cl),
      adm1_pop,
      group_short,
      sets,
      cl$start,
      mags,
      draws,
      cap = .cell_cap(cl),
      ...
    )
    tag <- function(x) .cell_tag(x, cl)
    out[[length(out) + 1L]] <- list(
      exposure = tag(o$exposure),
      adm1 = tag(o$adm1),
      n_min = o$n_min
    )
  }
  list(
    exposure = dplyr::bind_rows(lapply(out, `[[`, "exposure")),
    adm1 = dplyr::bind_rows(lapply(out, `[[`, "adm1")),
    n_min = min(vapply(out, `[[`, numeric(1), "n_min"))
  )
}

derive_attrib_exposure_posterior <- function(ctx, sets) {
  inp <- attrib_arm_inputs(ctx, sets, "attrib-exposure-posterior")
  d <- af_posterior_draws(ctx, ctx$MAN, inp$groups, sets)
  cell_tab <- attrib_cells(
    ctx,
    inp$attrib_start,
    display_factual_start(ctx, inp$groups, sets)
  )
  cells <- .cell_starts(cell_tab)
  message(sprintf(
    "[09] attrib_exposure_posterior: groups=%s | sets=%s | m=%s | cells=%s | max_draws=%s",
    paste(inp$groups, collapse = ","),
    paste(sets, collapse = ","),
    paste(ctx$ATTRIB_MAGNITUDES, collapse = ","),
    paste(names(cells), collapse = ","),
    paste(unique(attr(d, "max_draws")), collapse = ",")
  ))
  out <- attrib_exposure_posterior_windows(
    inp$adm1_month,
    inp$country_month,
    inp$adm1_pop,
    ctx$GROUP_SHORT_RUN,
    sets,
    cell_tab,
    ctx$ATTRIB_MAGNITUDES,
    d,
    min_reps = ctx$AF_MIN_REPS,
    guard_scope = ctx$AF_GUARD_SCOPE,
    noise_share_max = ctx$AF_NOISE_SHARE_MAX
  )
  c(
    out,
    list(
      windows = cells,
      magnitudes = ctx$ATTRIB_MAGNITUDES,
      analog_sets = sets,
      groups = inp$groups,
      max_draws = attr(d, "max_draws"),
      min_reps = ctx$AF_MIN_REPS,
      africa = "joint",
      guard_scope = ctx$AF_GUARD_SCOPE,
      noise_share_max = ctx$AF_NOISE_SHARE_MAX,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      attrib_start = inp$attrib_start,
      factual_start = inp$win$factual_start,
      pop_ref_year = ctx$POP_REF_YEAR,
      post_fork_start = cells[["post_fork"]],
      attrib_factual = ctx$ATTRIB_FACTUAL %||% "forecast",
      cells = cells_payload(cell_tab),
      africa_window = africa_window_payload(
        .africa_win(ctx, inp$groups),
        inp$attrib_start
      )
    )
  ) # attrib_exposure_posterior v3
}

ensure_attrib_exposure_posterior <- function(ctx, sets) {
  ctx <- results_attrib_ctx(ctx)
  man <- ctx$MAN
  groups <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)
  cfgs <- results_run_configs(man, groups, sets)
  inputs <- results_inputs(
    files = c(
      results_attrib_manifest_files(ctx),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = c(
      if (results_arm_window(ctx, groups, sets, read_events = FALSE)$stitched) {
        smoothing_parts_glob(ctx)
      } else {
        list()
      },
      .afp_assignment_globs(ctx, cfgs)
    )
  )
  ensure_artefact(
    ctx,
    "attrib_exposure_posterior",
    inputs,
    function() derive_attrib_exposure_posterior(ctx, sets),
    sets = sets
  )
}

# ===========================================================================
# report_attrib_information: the information component of the
# forecast attribution. With an
# attribution factual (results.inputs.attrib_manifest) three factual /
# counterfactual pairs share every replicate's smoothing replicate and
# posterior draw:
#   total        S - C   the displayed (stitched) factual vs the counterfactual
#   information  S - M   the displayed factual vs the attribution factual: what
#                        the months between the two launches add (the smoothed
#                        trajectory over them, and observed vs forecast ENSO
#                        inputs there) beyond re-simulating them
#   enso         M - C   the attribution factual vs the counterfactual: the ENSO
#                        trajectory only (af_posterior's headline under the
#                        attribution factual)
# so AN_total = AN_information + AN_enso per replicate (AF does not add). Each
# pair goes through af_posterior_tables() (per-draw AF / AN, noise-adjusted
# quantiles), relabelled as forecast / counterfactual. Needs the fixed draw
# design of af_posterior (results.af_posterior). Both factuals start from the
# same smoothing replicate of each replicate id (checked on the handoff maps):
# the attribution factual from its April state, which the smoother traced back
# from the August data, so the information component is "realised vs
# re-simulated May-Aug", not the value of the May-Aug data alone.
# ===========================================================================

#' The three pairs of report_attrib_information as af_posterior_tables() inputs:
#' list(total, information, enso), each list(reps, explosive) with the pair's
#' factual as run "forecast" (analog_set NA) and its reference as run
#' "counterfactual" (the information pair: analog_set "information").
#' @param reps_s,expl_s window totals / explosive keys of the displayed load
#'   (af_posterior_window_reps(): forecast = the stitched factual)
#' @param reps_a,expl_a the same of the attribution load (forecast = the
#'   attribution factual; counterfactual = the analog sets)
attrib_information_pairs <- function(reps_s, reps_a, expl_s, expl_a, sets) {
  pick <- function(x, rn) {
    out <- dplyr::filter(x, run == rn)
    if (rn == "counterfactual") {
      out <- dplyr::filter(out, analog_set %in% sets)
    }
    out
  }
  as_ref <- function(x) {
    dplyr::mutate(x, run = "counterfactual", analog_set = "information")
  }
  pair <- function(f_reps, c_reps, f_expl, c_expl) {
    list(
      reps = dplyr::bind_rows(f_reps, c_reps),
      explosive = dplyr::bind_rows(f_expl, c_expl)
    )
  }
  list(
    total = pair(
      pick(reps_s, "forecast"),
      pick(reps_a, "counterfactual"),
      pick(expl_s, "forecast"),
      pick(expl_a, "counterfactual")
    ),
    information = pair(
      pick(reps_s, "forecast"),
      as_ref(pick(reps_a, "forecast")),
      pick(expl_s, "forecast"),
      as_ref(pick(expl_a, "forecast"))
    ),
    enso = pair(
      pick(reps_a, "forecast"),
      pick(reps_a, "counterfactual"),
      pick(expl_a, "forecast"),
      pick(expl_a, "counterfactual")
    )
  )
}

#' Pure core of report_attrib_information: af_posterior_tables() of the three
#' pairs, rows tagged `component`, plus the decomposition table (an_pool of each
#' component per cell x analog set x guard x metric x scope, and the gap
#' total - (information + enso), 0 up to the replicates the explosive guard drops
#' differently per pair).
#' @param draws tibble(cell, unit, replicate_id, draw): the shared design
#' @return list(scope, draws, decomposition, n_min)
attrib_information_tables <- function(
  reps_s,
  reps_a,
  expl_s,
  expl_a,
  draws,
  sets,
  group_short,
  ...
) {
  pr <- attrib_information_pairs(reps_s, reps_a, expl_s, expl_a, sets)
  comp_sets <- list(total = sets, information = "information", enso = sets)
  tabs <- lapply(names(pr), function(k) {
    t <- af_posterior_tables(
      pr[[k]]$reps,
      pr[[k]]$explosive,
      draws,
      comp_sets[[k]],
      group_short,
      ...
    )
    t$scope <- dplyr::mutate(t$scope, component = k, .before = 1L)
    t$draws <- dplyr::mutate(t$draws, component = k, .before = 1L)
    t
  })
  scope <- dplyr::bind_rows(lapply(tabs, `[[`, "scope"))
  keys <- c(
    "cell",
    "explosive_guard",
    "metric",
    "level",
    "key",
    "scope",
    "enso_group"
  )
  info <- dplyr::filter(scope, component == "information") %>%
    dplyr::select(dplyr::all_of(keys), information = an_pool)
  decomposition <- dplyr::filter(scope, component != "information") %>%
    dplyr::select(
      dplyr::all_of(c(keys, "analog_set", "component")),
      an_pool
    ) %>%
    tidyr::pivot_wider(names_from = component, values_from = an_pool) %>%
    dplyr::left_join(info, by = keys) %>%
    dplyr::mutate(gap = total - (information + enso)) %>%
    dplyr::relocate(
      cell,
      analog_set,
      explosive_guard,
      metric,
      level,
      key,
      scope,
      enso_group,
      total,
      information,
      enso,
      gap
    )
  list(
    scope = scope,
    draws = dplyr::bind_rows(lapply(tabs, `[[`, "draws")),
    decomposition = decomposition,
    n_min = min(vapply(tabs, `[[`, numeric(1), "n_min"))
  )
}

#' Table S29 for the information decomposition: explosive_pair_summary()
#' of the stitched (displayed) factual against the counterfactual
#' (`stitched_vs_cf`) and against the attribution factual launched in May
#' (`stitched_vs_may`, set "information"), per window, from the per-replicate
#' country window totals of both loads (af_posterior_window_reps()).
attrib_information_explosive <- function(
  reps_s,
  reps_a,
  cap,
  sets,
  gmax_s = NULL,
  gmax_a = NULL
) {
  keys <- c("cell", "enso_group", "run", "analog_set", "replicate_id")
  grp <- function(r, gm) {
    out <- .sum_by_keys(
      dplyr::mutate(r, run = as.character(run)),
      keys,
      c(cum_cases = "cases")
    )
    # The per-gid cap reads each replicate's largest gid-month
    # (af_posterior_window_reps()$gmax).
    if (!is.null(gm) && nrow(gm)) {
      out <- dplyr::left_join(
        out,
        dplyr::mutate(gm, run = as.character(run)),
        by = keys
      )
    }
    out
  }
  gs <- grp(reps_s, gmax_s)
  ga <- grp(reps_a, gmax_a)
  dplyr::bind_rows(lapply(unique(ga$cell), function(k) {
    s <- dplyr::filter(gs, cell == k, run == "forecast")
    a <- dplyr::filter(ga, cell == k)
    dplyr::bind_rows(
      explosive_pair_summary(
        dplyr::bind_rows(s, dplyr::filter(a, run == "counterfactual")),
        cap,
        sets
      ) %>%
        dplyr::mutate(comparison = "stitched_vs_cf", .before = 1L),
      explosive_pair_summary(
        dplyr::bind_rows(
          s,
          dplyr::mutate(
            dplyr::filter(a, run == "forecast"),
            run = "counterfactual",
            analog_set = "information"
          )
        ),
        cap,
        "information"
      ) %>%
        dplyr::mutate(comparison = "stitched_vs_may", .before = 1L)
    ) %>%
      dplyr::mutate(window = k, .before = 1L)
  }))
}

#' Refuse two factual runs per unit that do not start each replicate from the
#' same smoothing replicate (their handoff_map.csv, 07b): the information
#' component would then mix in a change of smoothing particle.
attrib_information_check_handoff <- function(
  sim_root,
  man_s,
  man_a,
  groups,
  sets
) {
  cs <- results_unit_arm_configs(man_s, groups, sets)
  ca <- results_unit_arm_configs(man_a, groups, sets)
  for (un in names(ca)) {
    tags <- c(
      stitched = yaml::read_yaml(cs[[un]][[1L]])$simulator$run_tag,
      attrib = yaml::read_yaml(ca[[un]][[1L]])$simulator$run_tag
    )
    maps <- lapply(tags, function(tg) read_handoff_map(sim_root, tg))
    if (any(vapply(maps, is.null, logical(1)))) {
      stop(
        sprintf(
          "[09] report_attrib_information: unit %s: no handoff_map.csv under %s.",
          un,
          paste(tags[vapply(maps, is.null, logical(1))], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    j <- dplyr::inner_join(
      maps$stitched[, c("replicate_id", "smoothing_rep_id")],
      maps$attrib[, c("replicate_id", "smoothing_rep_id")],
      by = "replicate_id",
      suffix = c("_s", "_a")
    )
    if (
      nrow(j) != nrow(maps$attrib) ||
        any(j$smoothing_rep_id_s != j$smoothing_rep_id_a)
    ) {
      stop(
        sprintf(
          paste0(
            "[09] report_attrib_information: unit %s: the two factuals do not start every ",
            "replicate from the same smoothing replicate (%d of %d replicates differ)."
          ),
          un,
          sum(j$smoothing_rep_id_s != j$smoothing_rep_id_a) +
            nrow(maps$attrib) -
            nrow(j),
          nrow(maps$attrib)
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

derive_report_attrib_information <- function(ctx, sets) {
  cx_a <- results_attrib_ctx(ctx)
  groups <- manifest_groups_inc(ctx$MAN, ctx$EXCLUDE_GROUPS)
  pf <- display_factual_start(ctx, groups, sets)
  # The attribution load: this job's attribution steps' window totals when they
  # ran here (ctx$AF_POST_STASH), else loaded now; then the displayed load.
  st_a <- ctx$AF_POST_STASH[[paste(
    "forecast",
    paste(sort(sets), collapse = "-")
  )]]
  if (is.null(st_a)) {
    inp <- attrib_arm_inputs(cx_a, sets, "attrib-information")
    st_a <- af_posterior_stash_build(
      inp,
      attrib_cells(ctx, inp$attrib_start, pf)
    )
    inp <- NULL
  } else {
    message(
      "[09] report_attrib_information: reusing the attribution load's window totals"
    )
  }
  cells <- st_a$cells
  inp <- attrib_arm_inputs(ctx, sets, "attrib-information")
  if (as.Date(inp$attrib_start) != as.Date(cells[["attribution"]])) {
    stop(
      "[09] report_attrib_information: the displayed and attribution windows start differently.",
      call. = FALSE
    )
  }
  st_s <- af_posterior_stash_build(inp, st_a$cell_table)
  inp <- NULL

  d_a <- af_posterior_draws(cx_a, cx_a$MAN, groups, sets)
  d_s <- af_posterior_draws(ctx, ctx$MAN_DISPLAY %||% ctx$MAN, groups, sets)
  j <- dplyr::inner_join(
    d_a,
    d_s,
    by = c("unit", "replicate_id"),
    suffix = c("_a", "_s")
  )
  if (nrow(j) != nrow(d_a) || any(j$draw_a != j$draw_s)) {
    stop(
      "[09] report_attrib_information: the displayed and attribution factuals do not share the ",
      "replicate -> draw map.",
      call. = FALSE
    )
  }
  attrib_information_check_handoff(
    ctx$SIM_ROOT,
    ctx$MAN_DISPLAY %||% ctx$MAN,
    cx_a$MAN,
    groups,
    sets
  )
  draws <- dplyr::bind_rows(lapply(names(cells), function(k) {
    dplyr::mutate(d_a, cell = k, .before = 1L)
  }))
  message(sprintf(
    "[09] report_attrib_information: groups=%s | sets=%s | cells=%s | factual launches %s (display) / %s (attribution)",
    paste(groups, collapse = ","),
    paste(sets, collapse = ","),
    paste(names(cells), collapse = ","),
    format(pf),
    format(manifest_arm_starts(cx_a$MAN, groups, sets)$factual)
  ))
  tabs <- attrib_information_tables(
    st_s$reps,
    st_a$reps,
    st_s$explosive,
    st_a$explosive,
    draws,
    sets,
    ctx$GROUP_SHORT_RUN,
    guards = .cell_guards(st_a$cell_table),
    min_reps = ctx$AF_MIN_REPS,
    guard_scope = ctx$AF_GUARD_SCOPE,
    noise_share_max = ctx$AF_NOISE_SHARE_MAX
  )
  tabs$explosive_summary <- attrib_information_explosive(
    st_s$reps,
    st_a$reps,
    .ctx_explosive_cap(ctx),
    sets,
    st_s$gmax,
    st_a$gmax
  )
  c(
    tabs,
    list(
      windows = cells,
      analog_sets = sets,
      groups = groups,
      max_draws = attr(d_a, "max_draws"),
      display_factual_start = pf,
      attrib_factual_start = manifest_arm_starts(
        cx_a$MAN,
        groups,
        sets
      )$factual,
      min_reps = ctx$AF_MIN_REPS,
      africa = "joint",
      guard_scope = ctx$AF_GUARD_SCOPE,
      noise_share_max = ctx$AF_NOISE_SHARE_MAX,
      explosive_cap = ctx$EXPLOSIVE_CAP,
      explosive_gid_month_cap = ctx$EXPLOSIVE_GID_CAP %||% Inf,
      explosive_modes = .ctx_explosive_modes(ctx),
      cells = cells_payload(st_a$cell_table),
      africa_window = africa_window_payload(
        .africa_win(ctx, groups),
        cells[["attribution"]]
      )
    )
  ) # report_attrib_information v2
}

ensure_report_attrib_information <- function(ctx, sets) {
  groups <- manifest_groups_inc(ctx$MAN, ctx$EXCLUDE_GROUPS)
  cfgs <- .display_attrib_configs(ctx, groups, sets)
  inputs <- results_inputs(
    files = c(
      unique(c(ctx$MANIFEST_PATH, ctx$ATTRIB_MANIFEST_FILE)),
      cfgs,
      ctx$SURV_PATH,
      ctx$SIZE_DATA_PATH,
      ctx$GID_POP_PATH,
      ctx$GROUPING_PATH
    ),
    globs = c(smoothing_parts_glob(ctx), .afp_assignment_globs(ctx, cfgs))
  )
  ensure_artefact(
    ctx,
    "report_attrib_information",
    inputs,
    function() derive_report_attrib_information(ctx, sets),
    sets = sets
  )
}

# --- Figure 1d: rainfall anomaly vs ONI, rainy seasons ----------------------

#' Figure 1d's rainfall tables (rain_enso v1), per group of the results'
#' grouping (inputs.grouping), from 03_prepare_rain_adm1.R's ADM1 series
#' (utils/clim_compile.R):
#'   corr     Pearson r of the group-mean 3-month rainfall deviation (baseline
#'            sd units) with ONI at each lag, per centre month, with p_r, the
#'            Fisher-z interval and agree_r (share of the group's countries
#'            significant with the group's sign; n_ctry of them);
#'   seasons  each group's main rainy season(s), population-weighted (`pop`,
#'            ADM1 population in season_pop_year) and per ADM1 (`adm1`);
#'   month    the per-month share of the group in season behind `seasons`.
derive_rain_enso <- function(ctx) {
  d <- ctx$RAIN_DERIV
  bl <- as.integer(strsplit(as.character(d$baseline), "_", fixed = TRUE)[[1]])
  lags <- seq.int(as.integer(d$lags[[1]]), as.integer(d$lags[[2]]))
  grp <- load_enso_grouping(ctx$GROUPING_PATH)
  G <- setNames(grp$group, grp$iso3)
  adm1_m <- qs2::qs_read(ctx$RAIN_PATH)
  oni <- readRDS(.oni_cache_file(ctx)) %>%
    dplyr::transmute(
      year = as.integer(year),
      month = as.integer(month),
      oni = enso_index
    )
  dev <- precip_dev_series(adm1_m, bl, as.numeric(d$dry_mm_month))
  corr <- precip_oni_group_corr(dev, G, oni, lags, as.integer(d$min_n))
  pop <- adm1_population_from_gid(
    qs2::qs_read(ctx$GID_POP_PATH),
    as.integer(d$season_pop_year)
  )
  rs <- adm1_rain_seasons(
    adm1_m,
    G,
    pop,
    bl,
    as.numeric(d$season_share),
    as.numeric(d$season_min_annual_mm)
  )
  message(sprintf(
    "[09] rain_enso: %d ADM1 series (dry share %.2f), %d groups, lags %d..%d; seasons on %d ADM1s (%d with population)",
    length(dev$units),
    dev$dry_share,
    length(unique(G)),
    min(lags),
    max(lags),
    rs$n_units,
    rs$n_pop_matched
  ))
  out <- list(
    corr = corr,
    seasons = rs$seasons,
    month = rs$month,
    baseline = bl,
    lags = lags,
    months = range(dev$months),
    n_adm1 = length(dev$units)
  )
  # Table T1b's rainfall columns (v2), per season weighting and group.
  out$group_lags <- dplyr::bind_rows(lapply(
    unique(as.character(rs$seasons$weighting)),
    function(w) {
      rain_enso_lag_summary(out, sort(unique(G)), w)
    }
  ))
  out$mechanisms <- rain_enso_mechanisms(ctx$RAIN_MECH_PATH, sort(unique(G)))
  out # rain_enso v2
}

#' Table T1a's curated literature (results.inputs.enso_rain_mechanisms), read
#' and checked against the registry's groups: list(path, groups = <the YAML's
#' groups>). NULL when the input is "".
rain_enso_mechanisms <- function(path, groups) {
  if (!nzchar(path %||% "")) {
    return(NULL)
  }
  y <- yaml::read_yaml(require_input(
    path,
    "curated ENSO-rainfall mechanisms (table T1a)"
  ))
  miss <- setdiff(groups, names(y$groups))
  if (length(miss)) {
    stop(
      "[09] rain_enso: ",
      path,
      " has no entry for group(s) ",
      paste(miss, collapse = ", "),
      " (results.inputs.enso_rain_mechanisms).",
      call. = FALSE
    )
  }
  list(path = path, groups = y$groups)
}

#' rain_enso_lag_summary
#'
#' Table T1b's rainfall columns (supplement), per group: the rainy season(s)
#' (`seasons`, `weighting`), the in-season centre month and lag with the
#' largest |r| (`corr`), and the lag with the largest |mean r| over the
#' in-season centre months. The deployed cholera lags are joined in 11.
#' @param rain list(corr, seasons, month) of rain_enso
#' @return tibble(weighting, enso_group, rainy_season, best_centre_month,
#'   best_lag, best_r, best_p, best_n, best_agree, season_mean_peak_lag,
#'   season_mean_peak_r)
rain_enso_lag_summary <- function(rain, groups, weighting = "pop") {
  mon <- rain$month[
    rain$month$weighting == weighting & rain$month$in_season %in% TRUE,
    ,
    drop = FALSE
  ]
  sea <- rain$seasons[rain$seasons$weighting == weighting, , drop = FALSE]
  dplyr::bind_rows(lapply(groups, function(g) {
    ins <- mon$mo[mon$group == g]
    cr <- rain$corr[
      rain$corr$series == g &
        rain$corr$centre %in% ins &
        is.finite(rain$corr$r),
      ,
      drop = FALSE
    ]
    best <- if (nrow(cr)) cr[which.max(abs(cr$r)), , drop = FALSE] else NULL
    prof <- if (nrow(cr)) {
      dplyr::summarise(dplyr::group_by(cr, lag), r = mean(r), .groups = "drop")
    } else {
      NULL
    }
    pk <- if (!is.null(prof) && nrow(prof)) {
      prof[which.max(abs(prof$r)), , drop = FALSE]
    } else {
      NULL
    }
    tibble::tibble(
      weighting = weighting,
      enso_group = g,
      rainy_season = paste(sea$months[sea$group == g], collapse = "; "),
      best_centre_month = if (!is.null(best)) {
        month.abb[best$centre]
      } else {
        NA_character_
      },
      best_lag = if (!is.null(best)) as.integer(best$lag) else NA_integer_,
      best_r = if (!is.null(best)) best$r else NA_real_,
      best_p = if (!is.null(best)) best$p_r else NA_real_,
      best_n = if (!is.null(best)) as.integer(best$n) else NA_integer_,
      best_agree = if (!is.null(best) && "agree_r" %in% names(best)) {
        best$agree_r
      } else {
        NA_real_
      },
      season_mean_peak_lag = if (!is.null(pk)) {
        as.integer(pk$lag)
      } else {
        NA_integer_
      },
      season_mean_peak_r = if (!is.null(pk)) pk$r else NA_real_
    )
  }))
}

ensure_rain_enso <- function(ctx) {
  d <- ctx$RAIN_DERIV
  need <- c(
    "baseline",
    "lags",
    "dry_mm_month",
    "min_n",
    "season_share",
    "season_min_annual_mm",
    "season_pop_year"
  )
  if (is.null(d) || !all(need %in% names(d))) {
    stop(
      "[09] rain_enso: results.derivation.rain_enso must set ",
      paste(need, collapse = ", "),
      call. = FALSE
    )
  }
  require_input(
    ctx$RAIN_PATH,
    "ADM1 monthly rainfall series (Rscript analysis/03_prepare_rain_adm1.R)"
  )
  inputs <- results_inputs(
    files = c(
      ctx$RAIN_PATH,
      ctx$GROUPING_PATH,
      ctx$GID_POP_PATH,
      .oni_cache_file(ctx),
      ctx$RAIN_MECH_PATH
    )
  )
  ensure_artefact(ctx, "rain_enso", inputs, function() derive_rain_enso(ctx))
}

# ===========================================================================
# report_who_scale (v1, unkeyed, wave 1): the
# per-country WHO term scale of the smoothing, table D30.
# ===========================================================================

#' Files of report_who_scale: the smoothing run's who_scale.csv / meta (08a,
#' estimated mode) next to its parts/, and the merged ribbons (08c).
.who_scale_files <- function(ctx) {
  run_dir <- dirname(sub("/+$", "", ctx$SMOOTHING_PARTS_DIR))
  list(
    dir = run_dir,
    csv = file.path(run_dir, "who_scale.csv"),
    meta = file.path(run_dir, "who_scale_meta.yml"),
    ribbons = file.path(run_dir, "merged", "smoothing_ribbons.rds")
  )
}

#' The optional WHO-vs-database consistency sidecar of the national zero-fill:
#' the db_who_consistency() table persisted next
#' to the national-zero months file. NULL when absent (read only if present).
.who_db_consistency_path <- function(ctx) {
  thr <- national_zero_threshold_of(ctx$opt)
  sub(
    "\\.rds$",
    "_db_who_consistency.csv",
    as.character(make_national_zero_months_file(thr))
  )
}

#' One partition per TU: the first part.qs of each tu=<TU>/ directory (frame
#' metadata is the same in every chunk of a TU).
.who_scale_frame_meta <- function(ctx) {
  tu_dirs <- list.dirs(ctx$SMOOTHING_PARTS_DIR, recursive = FALSE)
  tu_dirs <- tu_dirs[grepl("^tu=", basename(tu_dirs))]
  if (!length(tu_dirs)) {
    stop(
      "[09] report_who_scale: no tu=* partitions under ",
      ctx$SMOOTHING_PARTS_DIR,
      call. = FALSE
    )
  }
  fm <- lapply(tu_dirs, function(d) {
    p <- sort(Sys.glob(file.path(d, "chunk_*", "part.qs")))
    if (!length(p)) {
      return(NULL)
    }
    qs2::qs_read(p[1])$frame_meta
  })
  dplyr::distinct(dplyr::bind_rows(fm), gid, .keep_all = TRUE)
}

#' derive_report_who_scale
#'
#' mode "estimated": the smoothing ran with who_obs.scale_mode estimated; the
#' table is 08a's who_scale.csv (the one the filter scored) and its meta.
#' mode "legacy_preview": a fixed-scale smoothing (older runs; who_scale.csv
#' absent). The implied ratios are recomputed with estimate_who_scale() on the
#' run's own frame (part.qs frame_meta: gids after the overlap drop, who_f at
#' the run's population year) and the filter's observation rule
#' (build_obs_panel(), monthly_case_obs zero_filled) over the default window
#' (simulator.filter.who_obs.scale_est defaults); the filter scored scale 1,
#' i.e. multiplier_scored = 1 / who_f. `scale` / `multiplier` there are what
#' the default estimator (kappa 0) would give, `scale_k500` with kappa = 500
#' -- a preview, not what the run used.
#' Both modes add the posterior check of 08c (WHO months inside the 95 % band)
#' when the merged ribbons carry it, and the db_who_consistency sidecar
#' when present.
derive_report_who_scale <- function(ctx) {
  f <- .who_scale_files(ctx)
  if (file.exists(f$csv)) {
    mode <- "estimated"
    tbl <- tibble::as_tibble(utils::read.csv(f$csv, stringsAsFactors = FALSE))
    meta <- if (file.exists(f$meta)) yaml::read_yaml(f$meta) else list()
    tbl$multiplier_scored <- tbl$multiplier
  } else {
    mode <- "legacy_preview"
    fm <- .who_scale_frame_meta(ctx)
    if (!"who_f" %in% names(fm)) {
      stop(
        "[09] report_who_scale: the smoothing partitions carry no who_f (WHO term off); ",
        "nothing to report.",
        call. = FALSE
      )
    }
    se <- get_default_options()$simulator$filter$who_obs$scale_est
    surv <- readRDS(require_input(ctx$SURV_PATH, "survival data"))
    mon <- readRDS(require_input(ctx$MON_DATA_PATH, "monthly case frame"))
    obs <- build_obs_panel(surv, mon, zero_fill = TRUE)
    who_path <- require_input(
      make_who_monthly_file(),
      "WHO dashboard monthly table (03_prepare_who_monthly.R)"
    )
    who <- who_obs_table(readRDS(who_path), se$start, se$end)
    who$cases_who <- round(who$cases_who)
    db <- who_db_observed_months(obs, fm, se$start, se$end)
    cov <- dplyr::summarise(
      dplyr::group_by(fm, country),
      who_f = dplyr::first(who_f),
      .groups = "drop"
    ) %>%
      dplyr::mutate(pop_covered = NA_real_, pop_national = NA_real_)
    est <- function(k) {
      estimate_who_scale(
        db,
        who,
        cov,
        se$start,
        se$end,
        min_db_cases = se$min_db_cases,
        min_months = se$min_months,
        shrink_cases = k,
        fallback = se$fallback,
        fixed_scale = 1,
        clamp = se$clamp
      )
    }
    tbl <- est(0)
    tbl$scale_k500 <- est(500)$scale[match(tbl$country, est(500)$country)]
    tbl$multiplier_scored <- who_multiplier(tbl$who_f, 1)
    meta <- list(
      scale_mode = "fixed",
      fixed_scale = 1,
      scale_est = list(
        start = format(as.Date(se$start)),
        end = format(as.Date(se$end)),
        min_db_cases = se$min_db_cases,
        min_months = se$min_months,
        fallback = se$fallback,
        clamp = se$clamp
      ),
      obs_rule = "build_obs_panel(zero_fill = TRUE)"
    )
  }
  post <- NULL
  if (file.exists(f$ribbons)) {
    rb <- tryCatch(readRDS(f$ribbons), error = function(e) NULL)
    post <- rb[["who_country"]][["summary"]] %||% rb[["who_cmp"]][["summary"]]
  }
  cons_path <- .who_db_consistency_path(ctx)
  cons <- if (file.exists(cons_path)) {
    tibble::as_tibble(utils::read.csv(cons_path, stringsAsFactors = FALSE))
  } else {
    NULL
  }
  message(sprintf(
    "[09] report_who_scale: %s, %d countries (%s)%s%s",
    mode,
    nrow(tbl),
    paste(
      sprintf("%s %d", names(table(tbl$source)), as.integer(table(tbl$source))),
      collapse = ", "
    ),
    if (is.null(post)) "" else "; posterior check from the merged ribbons",
    if (is.null(cons)) "" else "; db_who_consistency sidecar joined"
  ))
  # The national monthly series itself (03_prepare_who_monthly.R: WHO dashboard
  # weekly counts and sitrep differences per country-month), released as D36.
  wm_path <- make_who_monthly_file()
  national_monthly <- if (file.exists(wm_path)) {
    wm <- tibble::as_tibble(readRDS(wm_path))
    wm[
      order(wm$iso3, wm$t_start),
      intersect(
        c("iso3", "t_start", "cases", "n_days", "partial", "source"),
        names(wm)
      )
    ]
  } else {
    NULL
  }
  list(
    mode = mode,
    table = tbl,
    meta = meta,
    posterior = post,
    consistency = cons,
    national_monthly = national_monthly,
    smoothing_run_dir = f$dir
  ) # report_who_scale v2
}

ensure_report_who_scale <- function(ctx) {
  f <- .who_scale_files(ctx)
  cons <- .who_db_consistency_path(ctx)
  inputs <- results_inputs(
    files = c(
      if (file.exists(f$csv)) {
        c(f$csv, if (file.exists(f$meta)) f$meta)
      } else {
        c(ctx$SURV_PATH, ctx$MON_DATA_PATH)
      },
      if (file.exists(make_who_monthly_file())) make_who_monthly_file(),
      if (file.exists(f$ribbons)) f$ribbons,
      if (file.exists(cons)) cons
    )
  )
  ensure_artefact(ctx, "report_who_scale", inputs, function() {
    derive_report_who_scale(ctx)
  })
}

# ===========================================================================
# Counterfactual analog sets: table D29
# ===========================================================================
# What the counterfactual runs sampled from: per arm (the 2026 forecast, each
# past event) x analog set, the pool size, its model mix, the seam error (mean
# launch-month ONI of the set minus the observed ONI there) and the mean of each
# classified DJF winter, flagged under min_set_size. Read from the Lenssen
# caches themselves, so it covers legacy and library-mode caches.

#' Pure core: the D29 rows of one cache. `oni_lkp` = named published ONI
#' ("YYYY-MM-DD"). Returns list(sets, model_shares, winters).
analog_set_tables <- function(cache, arm, cache_path, oni_lkp) {
  m <- cache$meta %||% list()
  labels <- cache$set_labels %||%
    rep(NA_character_, nrow(cache$oni_trajectories))
  model <- cache$model %||% rep(NA_character_, length(labels))
  launch <- as.Date(m$horizon_start %||% NA)
  obs <- if (!is.na(launch)) {
    suppressWarnings(as.numeric(oni_lkp[format(launch)]))
  } else {
    NA_real_
  }
  min_size <- as.integer(
    m$spec$min_set_size %||%
      m$min_set_size %||%
      ENSO_COUNTERFACTUAL_DEFAULTS$min_set_size
  )
  match <- as.character(m$match %||% "legacy")
  def <- m$set_definition %||% list()
  prim <- as.integer(unlist(
    def$winters$primary %||% def$peak_season_year %||% NA_integer_
  ))
  stat_mat <- cache$stat_matrix %||%
    matrix(
      cache$set_stat %||% rep(NA_real_, length(labels)),
      ncol = 1L,
      dimnames = list(NULL, as.character(prim[1L]))
    )
  # Winter roles: primary; constrained (another read winter that must be
  # neutral, other_winters "neutral"); or
  # read_unconstrained (other_winters "any", e.g. an event registry's
  # other_winters: any): classified and recorded, constrains
  # no set. `constrained_winters` lists the primary + constrained ones, which
  # under "neutral" is every column, as before.
  yrs_all <- as.integer(colnames(stat_mat))
  unconstr <- identical(
    as.character(def$winters$other_winters %||% "neutral"),
    "any"
  )
  roles <- ifelse(
    yrs_all %in% prim,
    "primary",
    if (unconstr) "read_unconstrained" else "constrained"
  )
  yrs_con <- yrs_all[roles != "read_unconstrained"]
  sets <- c("all", setdiff(c(ENSO_ANALOG_SETS, "mixed"), "all"))
  rows <- lapply(sets, function(s) {
    idx <- if (identical(s, "all")) seq_along(labels) else which(labels == s)
    sh <- if (length(idx)) prop.table(table(model[idx])) else numeric(0)
    tibble::tibble(
      arm = arm,
      cache_path = cache_path,
      match = match,
      # every primary winter (";"-joined: a multi-winter event;
      # one winter otherwise, as before)
      launch = format(launch),
      primary_winter = paste(
        sprintf("%d/%02d", prim, (prim + 1L) %% 100L),
        collapse = ";"
      ),
      constrained_winters = paste(
        sprintf("%d/%02d", yrs_con, (yrs_con + 1L) %% 100L),
        collapse = ";"
      ),
      set = s,
      n = length(idx),
      share_of_pool = length(idx) / max(1L, length(labels)),
      n_models = length(sh),
      max_model_share = if (length(sh)) max(sh) else NA_real_,
      top_model = if (length(sh)) names(sh)[which.max(sh)] else NA_character_,
      seam_error = if (length(idx) && is.finite(obs)) {
        mean(cache$oni_trajectories[idx, 1L]) - obs
      } else {
        NA_real_
      },
      min_set_size = min_size,
      below_min_set_size = !identical(s, "all") && length(idx) < min_size,
      cache_digest = as.character(m$digest %||% NA_character_)
    )
  })
  shares <- lapply(sets, function(s) {
    idx <- if (identical(s, "all")) seq_along(labels) else which(labels == s)
    if (!length(idx)) {
      return(NULL)
    }
    sh <- prop.table(table(model[idx]))
    tibble::tibble(
      arm = arm,
      set = s,
      model = names(sh),
      share = as.numeric(sh)
    )
  })
  wint <- lapply(sets, function(s) {
    idx <- if (identical(s, "all")) seq_along(labels) else which(labels == s)
    if (!length(idx)) {
      return(NULL)
    }
    yrs <- yrs_all
    tibble::tibble(
      arm = arm,
      set = s,
      winter = sprintf("%d/%02d", yrs, (yrs + 1L) %% 100L),
      role = roles,
      mean_djf_oni = as.numeric(colMeans(stat_mat[idx, , drop = FALSE]))
    )
  })
  list(
    sets = dplyr::bind_rows(rows),
    model_shares = dplyr::bind_rows(shares),
    winters = dplyr::bind_rows(wint)
  )
}

#' The (arm, cache path) pairs of the context: the forecast counterfactual
#' configs' caches and, with an event manifest, each event's cache.
analog_set_caches <- function(ctx) {
  man <- ctx$MAN
  g0 <- manifest_groups_inc(man, ctx$EXCLUDE_GROUPS)[1]
  fc <- unique(vapply(
    ctx$SETS,
    function(a) {
      lenssen_counterfactual_args(
        yaml::read_yaml(group_config_path(
          man,
          g0,
          "counterfactual",
          a
        ))$simulator$enso$counterfactual,
        Sys.Date()
      )$cache_path
    },
    character(1)
  ))
  out <- tibble::tibble(arm = "forecast", cache_path = fc)
  if (nzchar(ctx$EVENT_MANIFEST_PATH %||% "")) {
    em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
    out <- dplyr::bind_rows(
      out,
      tibble::tibble(
        arm = names(em$events),
        cache_path = vapply(
          em$events,
          function(e) as.character(e$lenssen_cache %||% NA),
          character(1)
        )
      )
    )
  }
  out[!is.na(out$cache_path), , drop = FALSE]
}

derive_report_analog_sets <- function(ctx) {
  oni <- fetch_enso_index(
    source = "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  )
  oni_lkp <- stats::setNames(
    as.numeric(oni$enso_index),
    sprintf("%d-%02d-01", as.integer(oni$year), as.integer(oni$month))
  )
  cc <- analog_set_caches(ctx)
  parts <- lapply(seq_len(nrow(cc)), function(i) {
    analog_set_tables(
      readRDS(require_input(cc$cache_path[i], "Lenssen analog cache")),
      cc$arm[i],
      cc$cache_path[i],
      oni_lkp
    )
  })
  rd <- run_design_configs(ctx)
  list(
    sets = dplyr::bind_rows(lapply(parts, `[[`, "sets")),
    model_shares = dplyr::bind_rows(lapply(parts, `[[`, "model_shares")),
    winters = dplyr::bind_rows(lapply(parts, `[[`, "winters")),
    # v2: table D34, the replicate / draw / analog
    # design of every simulation arm, read from one run config per arm.
    run_design = run_design_table(rd)
  )
}

ensure_report_analog_sets <- function(ctx) {
  cc <- analog_set_caches(ctx)
  rd <- run_design_configs(ctx)
  inputs <- results_inputs(
    files = c(
      ctx$MANIFEST_PATH,
      ctx$EVENT_MANIFEST_PATH,
      cc$cache_path,
      .oni_cache_file(ctx),
      ctx$ATTRIB_MANIFEST_PATH,
      ctx$HINDCAST_MANIFEST_PATH,
      ctx$RES$deploy_config,
      stats::na.omit(c(rd$config, rd$cache_path))
    )
  )
  ensure_artefact(ctx, "report_analog_sets", inputs, function() {
    derive_report_analog_sets(ctx)
  })
}

# ===========================================================================
# Simulation run design: table D34, a slot of
# report_analog_sets v2
# ===========================================================================
# One row per simulation arm (forecast, counterfactual, attribution factual,
# past events, hindcast, smoothing): replicates, posterior draws, the draw
# design, the handoff weighting, the analog pool (per model, models) and the
# ENSO clamp margin -- read from one run config of the arm (every unit of an
# arm shares the template) and its Lenssen cache. The %||% defaults are those
# of get_default_options()$simulator and of the simulator's own readers
# (simulations.R / smoothing_core.R: max_draws 200).

RUN_DESIGN_ARMS <- c(
  "forecast",
  "counterfactual",
  "attribution_factual",
  "events",
  "hindcast",
  "smoothing"
)

#' Pure core: the D34 row of one arm.
#' @param cfg the run config (list), or NULL (arm absent: an all-NA row)
#' @param cache the Lenssen cache the arm's ENSO counterfactual samples from
#'   (list with $meta and $model), or NULL (no analog pool)
#' @param handoff FALSE for an arm that is not launched from smoothing states
#'   (the smoothing itself)
run_design_row <- function(
  arm,
  cfg,
  cache = NULL,
  handoff = TRUE,
  config = NA_character_,
  cache_path = NA_character_
) {
  sim <- cfg$simulator %||% list()
  have <- !is.null(cfg)
  n <- if (have) as.integer(sim$n_replicates %||% 100L) else NA_integer_
  S <- if (have) as.integer(sim$max_draws %||% 200L) else NA_integer_
  mode <- if (have) .sim_param_draw_mode(sim) else NA_character_
  enso <- sim$enso %||% list()
  clamp <- if (have && isTRUE(enso[["clamp_to_knots"]] %||% FALSE)) {
    as.numeric(enso[["clamp_margin"]] %||% 0.25)
  } else {
    NA_real_
  }
  m <- cache$meta %||% list()
  n_models <- m$n_models %||%
    (if (!is.null(cache$model)) length(unique(cache$model)) else NULL)
  tibble::tibble(
    arm = arm,
    n_replicates = n,
    max_draws = S,
    reps_per_draw = if (
      isTRUE(mode %in% c("cycle", "joint")) && !is.na(n) && !is.na(S) && S > 0L
    ) {
      n / S
    } else {
      NA_real_
    },
    draw_mode = mode,
    handoff = if (have && isTRUE(handoff)) {
      as.character(sim$forecast[["handoff_weighting"]] %||% "uniform")
    } else {
      NA_character_
    },
    n_analogs_per_model = as.integer(m$n_analogs_per_model %||% NA_integer_),
    n_models = as.integer(n_models %||% NA_integer_),
    n_analogs_total = if (!is.null(cache)) {
      as.integer(m$n_total_analogs %||% NROW(cache$oni_trajectories))
    } else {
      NA_integer_
    },
    analog_match = if (!is.null(cache)) {
      as.character(m$match %||% "legacy")
    } else {
      NA_character_
    },
    clamp_margin = clamp,
    clamp_to_knots = if (have) {
      isTRUE(enso[["clamp_to_knots"]] %||% FALSE)
    } else {
      NA
    },
    config = as.character(config),
    cache_path = as.character(cache_path)
  )
}

#' The first `config` leaf of a manifest list (hindcast manifests nest
#' units / groups x years x config).
.first_config_leaf <- function(x) {
  v <- unlist(x)
  hit <- v[grepl("(^|\\.)(factual_)?config$", names(v)) & nzchar(v)]
  if (length(hit)) as.character(hit[[1L]]) else NA_character_
}

#' One run config (and its Lenssen cache) per arm of the context: the first
#' partition unit of the forecast manifest (factual; the main-set
#' counterfactual), of the attribution manifest, of the first past event, of
#' the hindcast manifest, and the deploy's smoothing template. Arms the
#' context does not have keep config NA.
#' @return tibble(arm, config, cache_path, handoff)
run_design_configs <- function(ctx) {
  out <- tibble::tibble(
    arm = RUN_DESIGN_ARMS,
    config = NA_character_,
    cache_path = NA_character_,
    handoff = RUN_DESIGN_ARMS != "smoothing"
  )
  set1 <- function(a, cfg, cache = NA_character_) {
    i <- match(a, out$arm)
    out$config[i] <<- as.character(cfg %||% NA_character_)
    out$cache_path[i] <<- as.character(cache %||% NA_character_)
  }
  cf_cache <- function(path) {
    tryCatch(
      lenssen_counterfactual_args(
        yaml::read_yaml(path)$simulator$enso$counterfactual,
        Sys.Date()
      )$cache_path,
      error = function(e) NA_character_
    )
  }
  groups <- manifest_groups_inc(ctx$MAN, ctx$EXCLUDE_GROUPS)
  ua <- results_unit_arm_configs(ctx$MAN, groups, ctx$MAIN_SET)
  if (length(ua)) {
    set1("forecast", ua[[1L]][[1L]])
    cf <- ua[[1L]][[ctx$MAIN_SET]]
    if (!is.null(cf) && !is.na(cf) && file.exists(cf)) {
      set1("counterfactual", cf, cf_cache(cf))
    }
  }
  if (!is.null(ctx$MAN_ATTRIB)) {
    ua2 <- results_unit_arm_configs(ctx$MAN_ATTRIB, groups, character(0))
    if (length(ua2)) set1("attribution_factual", ua2[[1L]][[1L]])
  }
  if (
    nzchar(ctx$EVENT_MANIFEST_PATH %||% "") &&
      file.exists(ctx$EVENT_MANIFEST_PATH)
  ) {
    em <- read_event_manifest(ctx$EVENT_MANIFEST_PATH)
    if (length(em$events)) {
      set1(
        "events",
        em$events[[1L]]$units[[1L]]$factual_config %||%
          .first_config_leaf(em$events[[1L]]),
        em$events[[1L]]$lenssen_cache %||% NA_character_
      )
    }
  }
  hp <- as.character(ctx$HINDCAST_MANIFEST_PATH %||% "")
  if (length(hp) == 1L && nzchar(hp) && file.exists(hp)) {
    set1("hindcast", .first_config_leaf(yaml::read_yaml(hp)))
  }
  dp <- as.character(ctx$RES$deploy_config %||% "")
  if (length(dp) == 1L && nzchar(dp) && file.exists(dp)) {
    dy <- yaml::read_yaml(dp)
    set1(
      "smoothing",
      dy[["deploy"]][["templates"]][["historical"]] %||%
        dy[["templates"]][["historical"]]
    )
  }
  out$config[!is.na(out$config) & !file.exists(out$config)] <- NA_character_
  out$cache_path[
    !is.na(out$cache_path) & !file.exists(out$cache_path)
  ] <- NA_character_
  out
}

#' D34 from run_design_configs(): one run_design_row() per arm.
run_design_table <- function(rd) {
  dplyr::bind_rows(lapply(seq_len(nrow(rd)), function(i) {
    cfg <- if (!is.na(rd$config[i])) yaml::read_yaml(rd$config[i]) else NULL
    cache <- if (!is.na(rd$cache_path[i])) readRDS(rd$cache_path[i]) else NULL
    run_design_row(
      rd$arm[i],
      cfg,
      cache,
      rd$handoff[i],
      rd$config[i],
      rd$cache_path[i]
    )
  }))
}

# ---------------------------------------------------------------------------
# Supplementary table: episodes with a confirmed case vs all episodes of the
# analysis set (Table 9 and the within-country characteristics comparison)
# ---------------------------------------------------------------------------

# Tracks and contrasts of the table, in column order (classify_enso_phase() codes).
CONFIRMED_AGREEMENT_TRACKS <- c(
  occ = "Occurrence",
  dur = "Duration",
  mon = "Monthly cases"
)
CONFIRMED_AGREEMENT_CONTRASTS <- c("elnino", "nina")

#' Cluster bootstrap of `f(rows)` over `cluster` (location periods: one can
#' hold several outbreaks, which are not independent): c(est, lo, hi) at
#' `level`, from `B` resamples of whole clusters.
.agree_boot_cluster <- function(cluster, f, B, level) {
  idx <- split(seq_along(cluster), cluster)
  k <- length(idx)
  est <- f(seq_along(cluster))
  if (k < 2) {
    return(c(est = est, lo = NA_real_, hi = NA_real_))
  }
  bs <- replicate(
    B,
    f(unlist(idx[sample.int(k, k, replace = TRUE)], use.names = FALSE))
  )
  q <- stats::quantile(
    bs,
    c((1 - level) / 2, (1 + level) / 2),
    na.rm = TRUE,
    names = FALSE
  )
  c(est = est, lo = q[1], hi = q[2])
}

#' Ratio of the medians of x in `g` over the rest, with its cluster-bootstrap
#' interval; NA with fewer than 2 positive values on either side.
.agree_boot_ratio <- function(x, g, cluster, B, level) {
  ok <- is.finite(x) & x > 0
  x <- x[ok]
  g <- g[ok]
  cluster <- cluster[ok]
  if (sum(g) < 2 || sum(!g) < 2) {
    return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  .agree_boot_cluster(
    cluster,
    function(i) {
      a <- x[i][g[i]]
      b <- x[i][!g[i]]
      if (!length(a) || !length(b)) {
        NA_real_
      } else {
        stats::median(a) / stats::median(b)
      }
    },
    B,
    level
  )
}

#' Onset rate ratio of phase `ph` vs neutral from a quasi-Poisson GLM of the
#' monthly onset counts on phase (dispersion floored at 1), Wald interval at
#' `level`. `months` = the window's months (mi, phase); `onset_mi` = the onset
#' months of the outbreaks.
.agree_rate_ratio <- function(months, onset_mi, ph, level) {
  na <- c(est = NA_real_, lo = NA_real_, hi = NA_real_)
  u <- months
  u$n <- as.integer(table(factor(onset_mi, levels = u$mi)))
  tot <- tapply(u$n, u$phase, sum)
  keep <- names(tot)[tot > 0]
  if (!all(c("neutral", ph) %in% keep)) {
    return(na)
  }
  u <- u[u$phase %in% keep, , drop = FALSE]
  u$phase <- stats::relevel(factor(u$phase), "neutral")
  fit <- stats::glm(n ~ phase, family = stats::quasipoisson(), data = u)
  disp <- summary(fit)$dispersion
  V <- stats::vcov(fit) * max(1, disp) / disp
  k <- paste0("phase", ph)
  z <- stats::qnorm((1 + level) / 2)
  b <- stats::coef(fit)[[k]]
  se <- sqrt(V[k, k])
  c(est = exp(b), lo = exp(b - z * se), hi = exp(b + z * se))
}

#' confirmed_agreement_table
#'
#' ENSO contrasts against neutral (El Niño, La Niña) of three descriptive
#' tracks, for all outbreaks and for the subset with at least `min_confirmed`
#' confirmed cases, both restricted to the same onset window: occurrence (onset
#' rate ratio, quasi-Poisson over the window's months), duration (ratio of
#' median duration in weeks) and monthly cases (ratio of median suspected cases
#' per active month), the last two with a cluster bootstrap over location
#' periods. Phase: the threshold phase (+/-0.5) of the onset month's ONI.
#' Rows: Africa, each group, and the group's countries with confirmed
#' outbreaks; a confirmed contrast needs `min_neutral` neutral and `min_phase`
#' phase onsets, and a country row is kept when one of its confirmed contrasts
#' does. Pure (no file reads).
#'
#' @param ob episodes: country, location_period_id, TL, TR, suspected_cases,
#'   confirmed_cases (NA = not reported); duration (TR - TL) / 7 weeks, as
#'   report_dataset
#' @param members registry (iso3, group); @param group_order group ids in order
#' @param oni monthly ONI (year, month, enso_index)
#' @param cfg list(min_confirmed, min_neutral, min_phase, window ("confirmed" |
#'   "full"), boot_b, seed, level)
#' @return list(long = one row per row x set x track x contrast, counts = onsets
#'   per row x set x phase, characteristics_by_block = the within-country ratio,
#'   confirmed / other, of suspected cases and duration per row, window =
#'   c(start, end) "YYYY-MM", cfg)
confirmed_agreement_table <- function(ob, members, group_order, oni, cfg) {
  need <- c(
    "country",
    "location_period_id",
    "TL",
    "TR",
    "suspected_cases",
    "confirmed_cases"
  )
  if (!all(need %in% names(ob))) {
    stop(
      "confirmed_agreement_table: the outbreak table lacks ",
      paste(setdiff(need, names(ob)), collapse = ", "),
      call. = FALSE
    )
  }
  B <- as.integer(cfg$boot_b)
  level <- as.numeric(cfg$level)
  mi_of <- function(d) {
    as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m")) - 1L
  }
  months <- data.frame(
    mi = as.integer(oni$year) * 12L + as.integer(oni$month) - 1L,
    phase = classify_enso_phase(oni$enso_index, "threshold"),
    stringsAsFactors = FALSE
  )
  months <- months[!is.na(oni$enso_index), , drop = FALSE]
  d <- data.frame(
    country = as.character(ob$country),
    lp = as.character(ob$location_period_id),
    onset_mi = mi_of(as.Date(ob$TL)),
    duration_wk = as.numeric(as.Date(ob$TR) - as.Date(ob$TL)) / 7,
    suspected = as.numeric(ob$suspected_cases),
    confirmed = !is.na(ob$confirmed_cases) &
      ob$confirmed_cases >= cfg$min_confirmed,
    stringsAsFactors = FALSE
  )
  d$mon_cases <- d$suspected / (d$duration_wk / (365.25 / 12 / 7))
  d$group <- members$group[match(d$country, members$iso3)]
  d$phase <- months$phase[match(d$onset_mi, months$mi)]
  if (anyNA(d$phase)) {
    stop(
      "confirmed_agreement_table: no ONI for some onset months.",
      call. = FALSE
    )
  }
  win <- if (identical(cfg$window, "full")) {
    range(d$onset_mi)
  } else {
    range(d$onset_mi[d$confirmed])
  }
  d <- d[d$onset_mi >= win[1] & d$onset_mi <= win[2], , drop = FALSE]
  wm <- months[months$mi >= win[1] & months$mi <= win[2], , drop = FALSE]
  set.seed(as.integer(cfg$seed))

  rows <- list(list(block = "Africa", level = "africa", d = d))
  for (g in group_order) {
    dg <- d[!is.na(d$group) & d$group == g, , drop = FALSE]
    nc <- table(dg$country[dg$confirmed])
    ctry <- if (length(nc)) {
      names(nc)[order(-as.integer(nc), names(nc))]
    } else {
      character(0)
    }
    rows <- c(
      rows,
      list(list(block = g, level = "group", d = dg)),
      lapply(ctry, function(ct) {
        list(
          block = ct,
          level = "country",
          d = dg[dg$country == ct, , drop = FALSE]
        )
      })
    )
  }
  sel <- function(x, s) if (s == "conf") x[x$confirmed, , drop = FALSE] else x
  stat <- function(x, track, ph) {
    if (track == "occ") {
      return(.agree_rate_ratio(wm, x$onset_mi, ph, level))
    }
    y <- x[x$phase %in% c(ph, "neutral"), , drop = FALSE]
    v <- if (track == "dur") y$duration_wk else y$mon_cases
    .agree_boot_ratio(v, y$phase == ph, y$lp, B, level)
  }
  long <- dplyr::bind_rows(lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    dplyr::bind_rows(lapply(c("all", "conf"), function(s) {
      x <- sel(r$d, s)
      n_neu <- sum(x$phase == "neutral")
      dplyr::bind_rows(lapply(names(CONFIRMED_AGREEMENT_TRACKS), function(tr) {
        dplyr::bind_rows(lapply(CONFIRMED_AGREEMENT_CONTRASTS, function(ph) {
          n_ph <- sum(x$phase == ph)
          ok <- if (s == "conf") {
            n_neu >= cfg$min_neutral && n_ph >= cfg$min_phase
          } else {
            n_neu >= 1 && n_ph >= 1
          }
          ci <- if (ok) {
            stat(x, tr, ph)
          } else {
            c(est = NA_real_, lo = NA_real_, hi = NA_real_)
          }
          tibble::tibble(
            row = i,
            block = r$block,
            level = r$level,
            set = s,
            track = tr,
            contrast = ph,
            n_phase = n_ph,
            n_neutral = n_neu,
            computed = ok,
            est = unname(ci["est"]),
            lo = unname(ci["lo"]),
            hi = unname(ci["hi"])
          )
        }))
      }))
    }))
  }))
  keep <- unique(long$row[
    long$level != "country" | (long$set == "conf" & long$computed)
  ])
  long <- long[long$row %in% keep, , drop = FALSE]
  # A pair differs in sign when its two estimates fall strictly on opposite
  # sides of 1 (a ratio of exactly 1 sides with neither); sign_differs_ci
  # (v2) also needs at least one of the two intervals to exclude 1.
  excl1 <- function(lo, hi) !is.na(lo) & !is.na(hi) & (lo > 1 | hi < 1)
  long <- long %>%
    dplyr::group_by(row, track, contrast) %>%
    dplyr::mutate(
      sign_differs = sum(!is.na(est)) == 2 &&
        log(est[set == "all"]) * log(est[set == "conf"]) < 0,
      sign_differs_ci = sign_differs & any(excl1(lo, hi))
    ) %>%
    dplyr::ungroup()
  counts <- dplyr::bind_rows(lapply(keep, function(i) {
    dplyr::bind_rows(lapply(c("all", "conf"), function(s) {
      x <- sel(rows[[i]]$d, s)
      tibble::tibble(
        row = i,
        block = rows[[i]]$block,
        level = rows[[i]]$level,
        set = s,
        n_nina = sum(x$phase == "nina"),
        n_neutral = sum(x$phase == "neutral"),
        n_elnino = sum(x$phase == "elnino")
      )
    }))
  }))
  # Within-country ratio, confirmed / other, of suspected cases and duration
  # per kept row: the row's countries with a confirmed episode, onsets in the
  # window; lm(log y ~ confirmed + country) (no country term for one country),
  # cluster bootstrap over location periods. After the contrasts, so their
  # bootstrap draws are unchanged.
  fe_ratio <- function(dd, v) {
    dd <- dd[is.finite(dd[[v]]) & dd[[v]] > 0, , drop = FALSE]
    if (length(unique(dd$confirmed)) < 2) {
      return(NA_real_)
    }
    dd$y <- log(dd[[v]])
    f <- if (length(unique(dd$country)) > 1) {
      y ~ confirmed + country
    } else {
      y ~ confirmed
    }
    exp(stats::coef(stats::lm(f, data = dd))[["confirmedTRUE"]])
  }
  zq <- c((1 - level) / 2, (1 + level) / 2)
  characteristics <- dplyr::bind_rows(lapply(keep, function(i) {
    x <- rows[[i]]$d
    x <- x[x$country %in% unique(x$country[x$confirmed]), , drop = FALSE]
    idx <- split(seq_len(nrow(x)), x$lp)
    k <- length(idx)
    dplyr::bind_rows(
      lapply(
        c(suspected_cases = "suspected", duration_weeks = "duration_wk"),
        function(v) {
          est <- if (nrow(x)) fe_ratio(x, v) else NA_real_
          q <- if (is.na(est) || k < 2) {
            c(NA_real_, NA_real_)
          } else {
            bs <- replicate(
              B,
              fe_ratio(
                x[
                  unlist(
                    idx[sample.int(k, k, replace = TRUE)],
                    use.names = FALSE
                  ),
                  ,
                  drop = FALSE
                ],
                v
              )
            )
            stats::quantile(bs, zq, na.rm = TRUE, names = FALSE)
          }
          tibble::tibble(
            row = i,
            block = rows[[i]]$block,
            level = rows[[i]]$level,
            est = est,
            lo = q[1],
            hi = q[2],
            n_episodes = nrow(x),
            n_confirmed = sum(x$confirmed)
          )
        }
      ),
      .id = "variable"
    )
  }))
  fmt <- function(mi) sprintf("%d-%02d", mi %/% 12L, mi %% 12L + 1L)
  list(
    long = long,
    counts = counts,
    characteristics_by_block = characteristics,
    window = c(start = fmt(win[1]), end = fmt(win[2])),
    cfg = cfg
  )
}

#' confirmed_characteristics_table
#'
#' Episodes with a confirmed case against the other episodes of the countries
#' that report one, onsets in the confirmed window: the geometric-mean ratio
#' (confirmed / other) of suspected cases, duration (weeks, (TR - TL) / 7),
#' attack rate and population from lm(log y ~ confirmed + country), with a
#' cluster-bootstrap interval over location periods. Pure. Call after
#' confirmed_agreement_table() in the same RNG stream.
#' @param ob as confirmed_agreement_table(), plus pop
#' @return list(ratios = one row per variable (variable, est, lo, hi, n),
#'   coverage = episodes / confirmed / countries, countries = confirmed
#'   episodes per country)
confirmed_characteristics_table <- function(ob, cfg) {
  d <- data.frame(
    country = as.character(ob$country),
    lp = as.character(ob$location_period_id),
    TL = as.Date(ob$TL),
    suspected_cases = as.numeric(ob$suspected_cases),
    duration_weeks = as.numeric(as.Date(ob$TR) - as.Date(ob$TL)) / 7,
    pop = as.numeric(ob$pop),
    confirmed = !is.na(ob$confirmed_cases) &
      ob$confirmed_cases >= cfg$min_confirmed,
    stringsAsFactors = FALSE
  )
  d$attack_rate_1k <- ifelse(
    !is.na(d$pop) & d$pop > 0,
    1000 * d$suspected_cases / d$pop,
    NA_real_
  )
  nc <- table(d$country[d$confirmed])
  coverage <- tibble::tibble(
    n_episodes = nrow(d),
    n_countries = dplyr::n_distinct(d$country),
    n_confirmed = sum(d$confirmed),
    n_countries_confirmed = length(nc)
  )
  mi <- as.integer(format(d$TL, "%Y")) *
    12L +
    as.integer(format(d$TL, "%m")) -
    1L
  win <- if (identical(cfg$window, "full")) {
    range(mi)
  } else {
    range(mi[d$confirmed])
  }
  w <- d[d$country %in% names(nc) & mi >= win[1] & mi <= win[2], , drop = FALSE]
  idx <- split(seq_len(nrow(w)), w$lp)
  k <- length(idx)
  fe_ratio <- function(dd, v) {
    dd <- dd[is.finite(dd[[v]]) & dd[[v]] > 0, , drop = FALSE]
    if (dplyr::n_distinct(dd$confirmed) < 2) {
      return(NA_real_)
    }
    dd$y <- log(dd[[v]])
    f <- if (dplyr::n_distinct(dd$country) > 1) {
      y ~ confirmed + country
    } else {
      y ~ confirmed
    }
    exp(stats::coef(stats::lm(f, data = dd))[["confirmedTRUE"]])
  }
  z <- c((1 - cfg$level) / 2, (1 + cfg$level) / 2)
  ratios <- dplyr::bind_rows(lapply(
    c("suspected_cases", "duration_weeks", "attack_rate_1k", "pop"),
    function(v) {
      bs <- replicate(
        as.integer(cfg$boot_b),
        fe_ratio(
          w[unlist(idx[sample.int(k, k, replace = TRUE)], use.names = FALSE), ],
          v
        )
      )
      q <- stats::quantile(bs, z, na.rm = TRUE, names = FALSE)
      tibble::tibble(
        variable = v,
        est = fe_ratio(w, v),
        lo = q[1],
        hi = q[2],
        n_episodes = nrow(w),
        n_confirmed = sum(w$confirmed),
        n_countries = dplyr::n_distinct(w$country)
      )
    }
  ))
  list(
    ratios = ratios,
    coverage = coverage,
    countries = tibble::tibble(
      country = names(nc),
      n_confirmed = as.integer(nc)
    ) %>%
      dplyr::arrange(dplyr::desc(n_confirmed), country)
  )
}

#' report_confirmed_agreement (v2, unkeyed; v2 adds sign_differs_ci): confirmed_agreement_table() and
#' confirmed_characteristics_table() on the analysis set (the size dataset,
#' ADM0 excluded, as report_dataset), the results' grouping registry and the
#' ONI cache, under results.derivation.confirmed_agreement.
derive_report_confirmed_agreement <- function(ctx) {
  ob <- readRDS(require_input(ctx$SIZE_DATA_PATH, "size dataset (03)")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(as.character(admin_level) != "ADM0")
  members <- utils::read.csv(
    require_input(ctx$GROUPING_PATH, "grouping registry"),
    stringsAsFactors = FALSE
  )
  ord <- intersect(names(ctx$GROUP_SHORT_RUN), unique(members$group))
  ord <- c(ord, sort(setdiff(unique(members$group), ord)))
  oni <- fetch_enso_index(
    source = "oni",
    cache_dir = ctx$ENSO_CACHE_DIR,
    refresh = FALSE
  )
  out <- confirmed_agreement_table(ob, members, ord, oni, ctx$AGREE_DERIV)
  c(
    out,
    characteristics = list(confirmed_characteristics_table(ob, ctx$AGREE_DERIV))
  )
}

ensure_report_confirmed_agreement <- function(ctx) {
  inputs <- results_inputs(
    files = c(ctx$SIZE_DATA_PATH, ctx$GROUPING_PATH, .oni_cache_file(ctx))
  )
  ensure_artefact(ctx, "report_confirmed_agreement", inputs, function() {
    derive_report_confirmed_agreement(ctx)
  })
}
