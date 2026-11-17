# simulations.R — shared artefact loading and frame validation of the
#   outbreak simulator
#
# Provides:
#   load_artefacts()          fits, draws, panels and frame of one sim config
#   check_frame_vs_draws()    the location frame agrees with the fit's levels
#   .fit_opt_from_attr()      the options a fit was trained under
#
# Details:
#   Sourced by 07b_simulator_spatpomp.R, 03_prepare_sim_panel_projection.R and
#   the simulator tests; 07b_simulator_v1.R keeps its own copies (see
#   load_artefacts()).

# ---------------------------------------------------------------------------
# .fit_opt_from_attr
# ---------------------------------------------------------------------------

.fit_opt_from_attr <- function(fit) {
  # The fit-time opt is recorded as an attribute by 04*/05* scripts (best effort);
  # fall back to a minimal opt with the expected backend slots.
  opt_fit <- attr(fit, "opt")
  if (is.null(opt_fit)) {
    opt_fit <- list()
  }
  opt_fit
}

# ---------------------------------------------------------------------------
# .check_fit_pair_consistency
# ---------------------------------------------------------------------------
# Only the occurrence fit's baseline is read downstream (it resolves the
# survival panel and so the climate normalisation). A duration/monthly fit
# trained on another baseline would get covariates on the wrong climatology,
# silently. Warns rather than stops. opt_dur / opt_mon NULL skip that track.
.check_fit_pair_consistency <- function(
  opt_occ,
  label = NULL,
  opt_dur = NULL,
  opt_mon = NULL
) {
  b_occ <- opt_occ$baseline
  for (.tr in c("duration", "monthly")) {
    opt_tr <- if (identical(.tr, "duration")) opt_dur else opt_mon
    if (is.null(opt_tr)) {
      next
    }
    b_tr <- opt_tr$baseline
    if (!is.null(b_occ) && !is.null(b_tr) && !identical(b_occ, b_tr)) {
      warning(
        "[utils_simulations] ",
        if (!is.null(label)) paste0(label, ": ") else "",
        "occurrence fit was trained with baseline '",
        b_occ,
        "' but ",
        .tr,
        " fit with baseline '",
        b_tr,
        "'. Only the occurrence baseline is used to build the climate panel, ",
        "so ",
        .tr,
        " draws will be applied to covariates standardised against a ",
        "different climatology. Refit one track or set both to the same baseline.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# load_artefacts
# ---------------------------------------------------------------------------
#
# Loads the fits (occurrence, optional duration / monthly) and the climate
# panel, builds the location frame, draws bundles, per-replicate ENSO
# scenarios and, optionally, hybrid warmup panels. In bpfilter mode
# (simulator.filter.enable) it also attaches the survival panel as
# art$surv_panel; the caller joins its (gid, t_start, outbreak) onto
# panel_tu_warmup before build_obs_for_tu().

# ---------------------------------------------------------------------------
# .discover_group_fit — locate a per-group fit RDS in fits_dir by token glob.
#
# Per-group fits are written by the make_*_run_file family with the token
# `_grp-<source>-<group>`. Globbing by that token is more robust
# than reconstructing the full filename (which varies by backend/admin/basehaz).
# `which` is "occ" (prefix outbreaks_survival), "dur" (prefix
# outbreak_duration), or "mon" (prefix outbreak_monthly,
# make_monthly_run_file()); anything else is an error (the size track is no
# longer simulated). Only `_fit_` files match (re_stats/coef_summary/loo/
# cb_meta are excluded).
# ---------------------------------------------------------------------------
.discover_group_fit <- function(
  fits_dir,
  source,
  group,
  which,
  episode_key = "lp_season",
  national_zero_fill = FALSE,
  national_zero_threshold = 0L,
  occ_panel = NULL
) {
  # Fit-file stems: utils/paths.R::.FIT_FILE_STEM (shared with read_deployed_lag()).
  if (!(length(which) == 1L && which %in% names(.FIT_FILE_STEM))) {
    stop(
      sprintf(
        paste0(
          "[load_artefacts] .discover_group_fit(): unknown track ",
          "'%s' -- expected one of 'occ', 'dur', 'mon'."
        ),
        which
      ),
      call. = FALSE
    )
  }
  prefix <- .FIT_FILE_STEM[[which]]
  pattern <- file.path(
    fits_dir,
    sprintf("%s*_fit_*_grp-%s-%s_*.rds", prefix, source, group)
  )
  hits <- Sys.glob(pattern)
  # Episode key (utils/paths.R::.episode_key_path_token()): a gid_season fit
  # carries "_ep-gs" in its name and is invisible to an lp_season run, and
  # vice versa -- both may sit in one fits_dir without a "more than one
  # match" error. No "_ep-gs" file exists unless one was fitted, so the
  # default behaviour is unchanged.
  hits <- hits[
    grepl("_ep-gs", basename(hits), fixed = TRUE) ==
      identical(episode_key, "gid_season")
  ]
  # Occurrence-panel variant ("_nz[<thr>][-<rule>]", "_ae", "_a1c"): an exact
  # whole-suffix match (utils/paths.R::.occ_panel_variant_file_match()). `occ_panel`
  # (a spec) replaces the national_zero_* pair when given.
  if (identical(which, "occ")) {
    hits <- hits[
      if (is.null(occ_panel)) {
        .national_zero_fill_file_match(
          hits,
          national_zero_fill,
          national_zero_threshold
        )
      } else {
        .occ_panel_variant_file_match(hits, occ_panel)
      }
    ]
  }
  if (length(hits) == 0L) {
    stop(
      sprintf(
        "[load_artefacts] no %s fit for group '%s' (source '%s') in %s ",
        which,
        group,
        source,
        fits_dir
      ),
      sprintf("(glob: %s)", basename(pattern))
    )
  }
  if (length(hits) > 1L) {
    stop(
      sprintf(
        "[load_artefacts] %d %s fits match group '%s' in %s — ",
        length(hits),
        which,
        group,
        fits_dir
      ),
      "expected exactly one. Matches: ",
      paste(basename(hits), collapse = ", ")
    )
  }
  hits[[1L]]
}

# ---------------------------------------------------------------------------
# .lags_required_from_draws — collect the panel lag columns a fit consumes.
#
# Returns the `enso_lag<n>` FE coefficient names present in an occurrence /
# duration / monthly draws bundle plus any `<var>_lag<k>` columns implied by a cross-basis
# meta sidecar (max_lag per variable). Used to validate the simulator climate
# panel carries the union of every group's lags.
# ---------------------------------------------------------------------------
.lags_required_from_draws <- function(draws, cb_meta) {
  enso_fe <- grep("^enso_(lag|lead)[0-9]+$", colnames(draws$fe), value = TRUE)
  cb_cols <- character(0)
  if (!is.null(cb_meta)) {
    for (v in names(cb_meta)) {
      # Latent cross-basis variables (utils/dlnm.R) have no panel column by
      # construction -- `mso` is built from the duration model's state inside the
      # kernel month loop. Without this skip the `sl == 0L` branch below would
      # demand a panel column literally named "mso" and stop() the run.
      if (v %in% .LATENT_CB_VARS) {
        next
      }
      vm <- cb_meta[[v]]
      if (!is.null(vm$single_lag)) {
        sl <- as.integer(vm$single_lag)
        cb_cols <- c(
          cb_cols,
          if (identical(v, "enso")) {
            # negative = enso_lead<k>
            enso_lag_col_name(sl)
          } else if (sl == 0L) {
            v
          } else {
            sprintf("%s_lag%d", v, sl)
          }
        )
      } else {
        max_lag <- as.integer(vm$lag[2])
        if (identical(v, "enso")) {
          cb_cols <- c(cb_cols, sprintf("enso_lag%d", seq(0L, max_lag)))
        } else {
          cb_cols <- c(cb_cols, v, sprintf("%s_lag%d", v, seq_len(max_lag)))
        }
      }
    }
  }
  unique(c(enso_fe, cb_cols))
}

#' Replicate chunking of the ENSO scenarios (07b --replicate_offset /
#' --n_replicates_total, set by slurm/submit_forecast_groups.sh).
#' A chunk runs replicates offset + 1 .. offset + n of a TU, but the ENSO
#' samplers draw columns 1 .. n_rep from a fixed seed, so a chunk sampling only
#' its own n would give replicate offset + j the trajectory of replicate j, and
#' every chunk would repeat the first chunk's trajectories. With the TU's total
#' known, the loader samples all of them and keeps the chunk's columns, so
#' replicate r gets the same trajectory however the replicates are chunked.
#' NULL (no slicing, n_rep = n) when the run is not chunked, or when the total is
#' unknown (opt$replicate_total unset: earlier callers, byte-identical).
#' @return NULL or list(n_build = total, cols = offset + seq_len(n))
.enso_chunk_cols <- function(opt, n_replicates) {
  off <- as.integer(opt$replicate_offset %||% 0L)
  tot <- as.integer(opt$replicate_total %||% NA_integer_)
  if (is.na(tot) || (off == 0L && tot == n_replicates)) {
    return(NULL)
  }
  if (off < 0L || off + n_replicates > tot) {
    stop(
      sprintf(
        "[utils_simulations] replicate chunk %d..%d exceeds --n_replicates_total %d",
        off + 1L,
        off + n_replicates,
        tot
      ),
      call. = FALSE
    )
  }
  list(n_build = tot, cols = off + seq_len(n_replicates))
}

#' This chunk's columns of every ENSO scenario matrix (identity when `chunk` is NULL).
.enso_chunk_slice <- function(enso_scenarios, chunk) {
  if (is.null(chunk) || is.null(enso_scenarios)) {
    return(enso_scenarios)
  }
  enso_scenarios[] <- lapply(enso_scenarios, function(m) {
    m[, chunk$cols, drop = FALSE]
  })
  enso_scenarios
}

# 07b_simulator_v1.R keeps its own copies of these functions (side-by-side
# rule). If the logic here diverges from 07b_simulator_v1.R, the divergence is
# intentional (e.g. the surv_panel extension that 07b_simulator_v1.R does not
# need).
load_artefacts <- function(opt, enso_sim_months_start = NULL) {
  sim <- opt$simulator

  # -------------------------------------------------------------------------
  # Country-grouping mode: load one occ (+ optional dur/mon) fit per group,
  # route each country to its group's draws. Mutually exclusive with the
  # single occurrence_fit/duration_fit/monthly_fit keys.
  # -------------------------------------------------------------------------
  grp_cfg <- sim$grouping %||% list(enable = FALSE)
  is_grouped <- isTRUE(grp_cfg$enable)
  if (is_grouped) {
    return(.load_artefacts_grouped(
      opt,
      enso_sim_months_start = enso_sim_months_start
    ))
  }

  message("[utils_simulations] reading occurrence fit: ", sim$occurrence_fit)
  occ_fit <- readRDS(sim$occurrence_fit)
  message("[utils_simulations] reading climate panel:  ", sim$climate_panel)
  panel <- readRDS(sim$climate_panel)

  opt_occ <- .fit_opt_from_attr(occ_fit)

  if (is.null(opt_occ$occurrence)) {
    # Every supported occurrence fit is a brmsfit; the extractor dispatches on it.
    be <- "brms"
    opt_occ$occurrence <- list(backend = be)
  }

  # ---------------------------------------------------------------------------
  # Duration fit (optional). When simulator.duration_fit is unset, draws_dur /
  # dur_fit stay NULL and simulate_one_tu() takes the has_dur=FALSE path.
  # ---------------------------------------------------------------------------
  has_duration_fit <- !is.null(sim$duration_fit) && nzchar(sim$duration_fit)
  dur_fit <- NULL
  opt_dur <- NULL
  if (has_duration_fit) {
    message("[utils_simulations] reading duration fit:   ", sim$duration_fit)
    dur_fit <- readRDS(sim$duration_fit)
    opt_dur <- .fit_opt_from_attr(dur_fit)
    if (is.null(opt_dur$duration)) {
      be <- if (inherits(dur_fit, "brmsfit")) "brms" else "brms"
      opt_dur$duration <- list(backend = be)
    }
  }

  # ---------------------------------------------------------------------------
  # Monthly case-count fit (`mon` track, optional). When set, the kernels
  # emit a case count at every ACTIVE month with mu driven by months-since-onset
  # instead of once at onset. `active`/`mso` are duration-model latent states,
  # so a duration fit is a hard requirement -- there is no invented fallback.
  # ---------------------------------------------------------------------------
  has_monthly_fit <- !is.null(sim$monthly_fit) && nzchar(sim$monthly_fit)
  mon_fit <- NULL
  opt_mon <- NULL
  if (has_monthly_fit) {
    if (!has_duration_fit) {
      stop(
        "simulator.monthly_fit requires simulator.duration_fit -- the ",
        "monthly case-count kernel reads `active`/`mso`, which are duration-model ",
        "latent states."
      )
    }
    message("[utils_simulations] reading monthly fit:    ", sim$monthly_fit)
    mon_fit <- readRDS(sim$monthly_fit)
    opt_mon <- .fit_opt_from_attr(mon_fit)
    if (is.null(opt_mon$monthly)) {
      opt_mon$monthly <- list(
        backend = if (inherits(mon_fit, "brmsfit")) "brms" else "brms"
      )
    }
  }

  .check_fit_pair_consistency(opt_occ, opt_dur = opt_dur, opt_mon = opt_mon)

  # Layer 4 (artefact chain, see pop_stamp_sidecar_path() / fit_pop_stamp_path() in
  # utils.R): does the monthly fit's own companion stamp (if any) match the
  # population that is CURRENTLY on disk for opt_mon? The monthly track is the
  # population-offset case-count model. require_stamp = FALSE (the default) --
  # a fit predating the stamp must keep loading; a mismatch (both stamps
  # present but disagreeing) still stops the run.
  if (has_monthly_fit) {
    assert_fit_pop_provenance(
      opt_mon,
      sim$monthly_fit,
      what = "the monthly fit"
    )
  }

  # Episode definition (utils/provenance.R::sim_episode_check()): every fit and
  # every prep artefact this run pairs must share one episode digest whenever
  # any of them -- or the config's episodes.key -- is gid_season. Silent for a
  # historical (lp_season) run.
  sim_episode_check(
    fits = list(occ = occ_fit, dur = dur_fit, mon = mon_fit),
    preps = list(
      surv = make_surv_dat_file(list(
        baseline = opt_occ$baseline %||% opt$baseline %||% "1990_2019",
        episodes = opt$episodes
      )),
      size_data = (opt$simulator$location_frame %||% list())$size_data_path %||%
        make_size_dat_file(opt),
      climate_panel = sim$climate_panel
    ),
    what = "load_artefacts",
    episode_key = episode_key_of(opt)
  )

  # ---------------------------------------------------------------------------
  # DLNM cross-basis metadata sidecars (optional; mirrors 09_simulator.R).
  # Derived from the explicit fit paths by string transform; config override
  # honoured if present. NULL when the fit carries no cross-basis.
  # ---------------------------------------------------------------------------
  .load_cb_meta <- function(fit_path, override) {
    cb_path <- override %||%
      sub("_fit_base-", "_cb_meta_base-", fit_path, fixed = TRUE)
    if (!is.null(cb_path) && file.exists(cb_path)) {
      message("[utils_simulations] reading cross-basis meta: ", cb_path)
      readRDS(cb_path)
    } else {
      NULL
    }
  }
  cb_meta_occ <- .load_cb_meta(sim$occurrence_fit, sim$occurrence_cb_meta)
  cb_meta_mon <- if (has_monthly_fit) {
    .load_cb_meta(sim$monthly_fit, sim$monthly_cb_meta)
  } else {
    NULL
  }
  # Duration cross-basis. `04_fit_dur_brms_cb.R` exists and the deployed
  # per-group best-lag duration fits are `_cb` fits carrying `enso_cb<k>` fixed
  # effects, so this sidecar is NOT optional. Without it the duration
  # LP's enso_cb coefficients multiply nothing (evaluate_lp() skips a FE whose
  # covariate is absent) and the track is silently ENSO-blind.
  cb_meta_dur <- if (has_duration_fit) {
    .load_cb_meta(sim$duration_fit, sim$duration_cb_meta)
  } else {
    NULL
  }

  # ADM1-ADM2 overlap filter's authoritative admin_level/ADM1 columns (loaded
  # here, ahead of build_location_frame(), so the same read also supplies
  # `admin_level_data`).
  .surv_overlap <- tryCatch(
    readRDS(make_surv_dat_file(list(
      baseline = opt_occ$baseline %||% opt$baseline %||% "1990_2019",
      episodes = opt$episodes
    ))),
    error = function(e) NULL
  )

  lf_opts <- opt$simulator$location_frame %||% list()
  .fpop <- sim_frame_pop_args(opt$simulator, opt_occ)
  .fpop$pop_window <- sim_pop_window(opt$simulator)
  frame <- build_location_frame(
    panel,
    hmm_clusters_path = lf_opts$hmm_clusters_path,
    outbreak_metadata_path = lf_opts$outbreak_metadata_path,
    size_data_path = lf_opts$size_data_path,
    # NULL (hence fully inert) unless one of the three fits was actually made
    # with `<track>.offset_baserate: true` -- see resolve_base_rate_paths().
    base_rate_table_paths = resolve_base_rate_paths(
      list(occ = opt_occ, dur = opt_dur, mon = opt_mon),
      lf_opts
    ),
    admin_level_data = .surv_overlap,
    # make_pop_run_file() reads only `outdir`/`population` from the opt, so the
    # occurrence fit's opt resolves the same global population artefacts.
    # Population arguments in one place (simulator_core.R::sim_frame_pop_args():
    # location_frame.gid_pop_path / wpp_path / pop_year).
    gid_pop_path = .fpop$gid_pop_path,
    wpp_path = .fpop$wpp_path,
    pop_year = .fpop$pop_year,
    pop_window = .fpop$pop_window
  )
  # gids whose min_seas was defaulted, checked on the final frame.
  .ms_defaulted <- attr(frame, "min_seas_defaulted")
  message(sprintf(
    "[utils_simulations] location frame: %d locations, %d ADM1, %d countries, %d TU(s) (incl. NA)",
    nrow(frame),
    length(unique(frame$adm1_unit)),
    length(unique(frame$country)),
    length(unique(frame$tu))
  ))

  # Apply scope.countries filter (ISO3 subset; severs TU-lag → warn)
  ctry_scope <- sim$scope$countries %||% "all"
  if (!identical(ctry_scope, "all")) {
    iso3_req <- trimws(strsplit(paste(ctry_scope, collapse = ","), ",")[[1L]])
    n_before <- nrow(frame)
    frame <- frame %>%
      dplyr::filter(stringr::str_extract(adm1_unit, "^[A-Z]{3}") %in% iso3_req)
    if (nrow(frame) == 0L) {
      stop(sprintf(
        "scope.countries filter ('%s') left 0 locations — check ISO3 codes.",
        paste(iso3_req, collapse = ", ")
      ))
    }
    message(sprintf(
      "[utils_simulations] scope.countries='%s': %d → %d locations",
      paste(iso3_req, collapse = ", "),
      n_before,
      nrow(frame)
    ))
    warning(paste0(
      "scope.countries severs the TU-lag: o_tu reflects only the requested ",
      "countries, not the full TU at fit time. Use scope.countries='all' for ",
      "unbiased TU-lag evaluation."
    ))
  }

  max_draws <- as.integer(sim$max_draws %||% 200L)
  draws_occ <- extract_occurrence_draws(occ_fit, opt_occ, max_draws = max_draws)
  draws_dur <- if (has_duration_fit) {
    extract_duration_draws(dur_fit, opt_dur, max_draws = max_draws)
  } else {
    NULL
  }
  draws_mon <- if (has_monthly_fit) {
    extract_monthly_draws(mon_fit, opt_mon, max_draws = max_draws)
  } else {
    NULL
  }

  # Zero-filled climate gate: a fit with a local-climate term must not be
  # replayed on the panel's zero-filled post-DB months. Inert without one.
  .assert_no_zero_filled_climate(
    .draws_fe_names(draws_occ, draws_dur, draws_mon),
    panel,
    window = .sim_panel_window(sim, enso_sim_months_start),
    label = "load_artefacts"
  )

  # ADM1-ADM2 overlap filter — mirror the smoothing pipeline (smoothing_core.R)
  # so the forecast frame matches the smoothed frame it consumes. Uses the survival
  # data's authoritative admin_level/ADM1 columns; `.surv_overlap` was already read
  # above (also supplied build_location_frame()'s `admin_level_data`), so no second
  # copy is loaded here.
  frame <- drop_overlapping_adm1_from_frame(
    frame,
    .surv_overlap,
    "utils_simulations"
  )
  # TTLO level "none": each gid's first recorded onset, from
  # which the kernels seed no-prior-onset units; only when the fit has the level.
  if (.draws_has_ttlo_none(draws_occ)) {
    frame$first_onset <- frame_first_onset(.surv_overlap, frame$gid)
  }
  rm(.surv_overlap)

  frame <- check_frame_vs_draws(
    frame,
    draws_occ,
    sim$scope,
    draws_mon = draws_mon
  )
  # The season onset gate needs every simulated gid's own season start.
  check_onset_eligibility_frame(
    frame,
    onset_eligibility_of(opt$simulator),
    "load_artefacts",
    defaulted = .ms_defaulted
  )

  n_replicates <- as.integer(sim$n_replicates %||% 100L)
  enso_chunk <- .enso_chunk_cols(opt, n_replicates) # chunked: sample the TU's total, keep ours
  n_enso <- enso_chunk$n_build %||% n_replicates

  # ---------------------------------------------------------------------------
  # Per-replicate ENSO forecast scenarios
  # Three branches (mutually exclusive; counterfactual key takes precedence):
  #   enso.counterfactual.scenario = "analog"  → sample_enso_analog_scenarios()
  #   enso.counterfactual.scenario = "lenssen" → sample_enso_lenssen_scenarios()
  #   (default / factual, enso.c3s_plume set)  → sample_enso_forecast_trajectories()
  # ---------------------------------------------------------------------------
  enso_scenarios <- NULL
  ctf_scenario <- sim$enso$counterfactual$scenario # NULL for factual runs

  forecast_source <- sim$enso$forecast_source %||% "cfsv2"

  # Observed-ENSO ("perfect foresight") hindcast mode. The entire forecast
  # window is historical, so 03_prepare_sim_panel.R's panel already carries the
  # observed enso_lag<N> columns for those months. Leave enso_scenarios = NULL
  # so simulate_one_tu() reads those panel columns unmodified — no CFSv2/C3S
  # override. Factual-only (a counterfactual scenario still takes precedence).
  observed_enso <- is.null(ctf_scenario) &&
    identical(forecast_source, "observed")
  if (observed_enso) {
    message(
      "[utils_simulations] ENSO forecast_source=observed (perfect foresight): ",
      "using panel enso_lag columns unmodified (no scenario override)."
    )
  } else if (
    !is.null(ctf_scenario) ||
      !is.null(sim$enso$c3s_plume) ||
      identical(forecast_source, "cfsv2")
  ) {
    tryCatch(
      {
        enso_source <- sim$enso$source %||% "oni"
        historical_oni <- fetch_enso_index(
          source = enso_source,
          cache_dir = "generated_data/",
          refresh = FALSE
        )
        # Derive simulation months from the panel (t_start >= horizon start).
        # In smoothing-input mode the effective per-TU forecast window can start
        # earlier than sim$horizon$start (at the auto-derived handoff). When the
        # caller supplies that earlier start (enso_sim_months_start), extend the
        # ENSO scenario window back to cover it so the per-replicate override
        # matrices span the full simulated horizon; otherwise they end up shorter
        # than the panel and get silently recycled downstream (see the ENSO
        # override in simulate_one_tu / simulator_core.R).
        horizon_start <- as.Date(sim$horizon$start)
        enso_start <- horizon_start
        if (!is.null(enso_sim_months_start)) {
          enso_start <- min(horizon_start, as.Date(enso_sim_months_start))
          if (enso_start < horizon_start) {
            message(sprintf(
              "[utils_simulations] ENSO scenario window extended back to %s (< sim horizon %s) to cover the smoothing handoff.",
              enso_start,
              horizon_start
            ))
          }
        }
        sim_months <- sort(unique(panel$t_start[panel$t_start >= enso_start]))

        if (!is.null(ctf_scenario)) {
          # ── Counterfactual branches ─────────────────────────────────────────
          message(
            "[utils_simulations] ENSO counterfactual scenario: ",
            ctf_scenario
          )

          if (identical(ctf_scenario, "analog")) {
            analog_start <- as.Date(
              sim$enso$counterfactual$analog_start %||%
                stop(
                  "enso.counterfactual.analog_start is required for the analog scenario."
                )
            )
            analog_end <- as.Date(
              sim$enso$counterfactual$analog_end %||%
                stop(
                  "enso.counterfactual.analog_end is required for the analog scenario."
                )
            )
            enso_scenarios <- sample_enso_analog_scenarios(
              historical_oni = historical_oni,
              sim_months = sim_months,
              analog_start = analog_start,
              analog_end = analog_end,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              n_rep = n_enso
            )
            message(sprintf(
              "[utils_simulations] analog ENSO: %d lags, %d months, %d replicates (deterministic)",
              length(enso_scenarios),
              if (length(enso_scenarios) > 0L) {
                nrow(enso_scenarios[[1L]])
              } else {
                0L
              },
              if (length(enso_scenarios) > 0L) {
                ncol(enso_scenarios[[1L]])
              } else {
                0L
              }
            ))
          } else if (identical(ctf_scenario, "lenssen")) {
            # One resolution of the counterfactual block (utils/enso.R::
            # lenssen_counterfactual_args()): cache path, analog
            # anchor (lenssen_horizon_start, NULL / "" = horizon_start), offsets
            # and analog set (`%||% "all"` there is the effective default).
            lcf <- lenssen_counterfactual_args(
              sim$enso$counterfactual,
              horizon_start
            )
            enso_horizon <- lcf$horizon_start
            if (enso_horizon != horizon_start) {
              message(sprintf(
                "[utils_simulations] Lenssen ENSO horizon: %s (overrides sim horizon %s)",
                enso_horizon,
                horizon_start
              ))
            }
            enso_scenarios <- sample_enso_lenssen_scenarios(
              cache_path = lcf$cache_path,
              sim_months = sim_months,
              n_rep = n_enso,
              historical_oni = historical_oni,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              horizon_start = enso_horizon,
              mean_offset = lcf$mean_offset,
              log_var_scaling = lcf$log_var_scaling,
              seed = as.integer(sim$seed %||% 20260603L),
              analog_set = lcf$analog_set
            )
            message(sprintf(
              "[utils_simulations] Lenssen ENSO: %d lags, %d months × %d replicates (analog_set=%s)",
              length(enso_scenarios),
              if (length(enso_scenarios) > 0L) {
                nrow(enso_scenarios[[1L]])
              } else {
                0L
              },
              if (length(enso_scenarios) > 0L) {
                ncol(enso_scenarios[[1L]])
              } else {
                0L
              },
              as.character(sim$enso$counterfactual$analog_set %||% "all")
            ))
          } else {
            stop(
              "[utils_simulations] Unknown enso.counterfactual.scenario: '",
              ctf_scenario,
              "'. Must be 'analog' or 'lenssen'."
            )
          }
        } else {
          # ── Factual branch: forecast_source = cfsv2 (default) | c3s ───────────
          # Raw Niño 3.4 for the rolling-mean conversion (shared by both sources).
          historical_nino34 <- tryCatch(
            fetch_enso_index(
              source = "nino34",
              cache_dir = "generated_data/",
              refresh = FALSE
            ),
            error = function(e) {
              message(
                "[utils_simulations] nino34 cache not found — using ONI as Niño 3.4 proxy."
              )
              NULL
            }
          )

          if (identical(forecast_source, "cfsv2")) {
            message(
              "[utils_simulations] building per-replicate ENSO trajectories from CFSv2 ensemble"
            )
            cfsv2_data <- fetch_enso_cfsv2_ensemble(
              url = sim$enso$cfsv2_url %||%
                "https://www.cpc.ncep.noaa.gov/products/CFSv2/dataInd3/nino34Mon.nc",
              cache_dir = sim$enso$cache_dir %||% "generated_data/",
              init = sim$enso$cfsv2_init %||% "latest"
            )
            message(sprintf(
              "[utils_simulations] CFSv2 ensemble: %d rows, target dates %s – %s, %d members",
              nrow(cfsv2_data),
              min(cfsv2_data$target_date),
              max(cfsv2_data$target_date),
              length(unique(cfsv2_data$member))
            ))
            enso_scenarios <- sample_enso_cfsv2_trajectories(
              cfsv2_data = cfsv2_data,
              n_rep = n_enso,
              historical_oni = historical_oni,
              historical_nino34 = historical_nino34,
              sim_months = sim_months,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              rebase = sim$enso$cfsv2_rebase %||% "seam",
              sd_correction = isTRUE(sim$enso$cfsv2_sd_correction),
              target_peak = sim$enso$cfsv2_target_peak %||% NULL,
              calibration = if (nzchar(sim$enso$cfsv2_calibration %||% "")) {
                sim$enso$cfsv2_calibration
              } else {
                NULL
              },
              calibration_target = sim$enso$cfsv2_calibration_target %||%
                "ersst",
              seed = as.integer(sim$seed %||% 20260603L),
              # past the CFSv2 span: neutral (default) | ar1 | analog (simulator.enso.continuation)
              beyond_cfsv2 = (sim$enso$continuation %||% list())$method %||%
                "neutral",
              ar1_phi = (sim$enso$continuation %||% list())$phi,
              ar1_fit_window = (sim$enso$continuation %||%
                list())$fit_window %||%
                "1991_2020",
              legacy_seam = isTRUE(
                (sim$enso$continuation %||% list())$legacy_seam
              ),
              analog_library = (sim$enso$continuation %||% list())[[
                "library"
              ]] %||%
                "generated_data/enso_continuation_library.qs",
              analog_k = (sim$enso$continuation %||% list())[["k"]] %||% 30L,
              analog_blend_months = (sim$enso$continuation %||% list())[[
                "blend_months"
              ]] %||%
                3,
              analog_min_sep = (sim$enso$continuation %||% list())[[
                "min_sep"
              ]] %||%
                24L,
              analog_max_per_model = (sim$enso$continuation %||% list())[[
                "max_per_model"
              ]] %||%
                Inf,
              analog_exclude_models = unlist((sim$enso$continuation %||%
                list())[["exclude_models"]]),
              # analog seam: last_cfsv2 (default) | last_full_oni
              analog_seam = (sim$enso$continuation %||% list())[["seam"]] %||%
                "last_cfsv2"
            )
          } else {
            message(
              "[utils_simulations] building per-replicate ENSO forecast trajectories (C3S plume)"
            )
            plume_data <- read_enso_c3s_plume(sim$enso$c3s_plume)
            message(sprintf(
              "[utils_simulations] C3S plume: %d rows, dates %s – %s, quantile levels %s",
              nrow(plume_data),
              min(plume_data$date),
              max(plume_data$date),
              paste(sort(unique(plume_data$q_level)), collapse = ", ")
            ))
            enso_scenarios <- sample_enso_forecast_trajectories(
              plume_data = plume_data,
              n_rep = n_enso,
              historical_oni = historical_oni,
              historical_nino34 = historical_nino34,
              sim_months = sim_months,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              seed = as.integer(sim$seed %||% 20260603L)
            )
          }
          message(sprintf(
            "[utils_simulations] ENSO scenarios built: %d lags, %d months × %d replicates",
            length(enso_scenarios),
            if (length(enso_scenarios) > 0L) nrow(enso_scenarios[[1]]) else 0L,
            if (length(enso_scenarios) > 0L) ncol(enso_scenarios[[1]]) else 0L
          ))
        }
      },
      error = function(e) {
        # Was: warn + enso_scenarios <<- NULL (silent fallback to the shared
        # panel's ENSO columns, which are unpopulated/NA in the forecast horizon
        # for projection panels). That silently swallowed a real config bug
        # (as.Date(NULL) + is.na() on lenssen_horizon_start) and produced a
        # spurious outbreak spike at the counterfactual handoff. Fail loudly
        # instead — a broken ENSO scenario must not be allowed to silently
        # degrade into unoverridden panel ENSO.
        stop(
          "[utils_simulations] ENSO scenario build failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  } else if (is.null(ctf_scenario) && is.null(sim$enso$c3s_plume)) {
    # Reached only when forecast_source is a factual value we do not recognise.
    # Fail loudly rather than silently fall through to unoverridden panel ENSO
    # (a typo like "cfsv2x" must not quietly change the driver).
    stop(
      sprintf(
        "[utils_simulations] Unknown enso.forecast_source '%s' (expected 'cfsv2', 'c3s', or 'observed').",
        forecast_source
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Per-replicate warmup panels (optional, from warmup config key)
  # ---------------------------------------------------------------------------
  warmup_panels_list <- NULL
  if (isTRUE(sim$warmup$use_sim_trajectories)) {
    message("[utils_simulations] building per-replicate hybrid warmup panels")
    tryCatch(
      {
        sim_output_path <- sim$warmup$sim_output_path
        if (is.null(sim_output_path) || !nzchar(sim_output_path)) {
          stop(
            "simulator.warmup.sim_output_path is required when use_sim_trajectories=true."
          )
        }
        if (!file.exists(sim_output_path)) {
          stop("simulator.warmup.sim_output_path not found: ", sim_output_path)
        }

        surv_path <- tryCatch(make_surv_dat_file(opt_occ), error = function(e) {
          NULL
        })
        if (is.null(surv_path) || !file.exists(surv_path %||% "")) {
          surv_panel_warmup <- panel
          message(
            "[utils_simulations] warmup: using projection panel as surrogate surv_panel ",
            "(survival panel path not resolved from fit opt)."
          )
        } else {
          message(
            "[utils_simulations] warmup: reading survival panel: ",
            surv_path
          )
          surv_panel_warmup <- readRDS(surv_path)
        }

        horizon_start <- as.Date(sim$horizon$start)
        warmup_months <- as.integer(sim$warmup$warmup_months %||% 12L)
        warmup_start <- horizon_start %m-% months(warmup_months)
        warmup_end <- horizon_start

        scoring_start <- tryCatch(
          as.Date(
            sim$warmup$scoring_start %||%
              format(warmup_end %m-% months(21L), "%Y-%m-%d")
          ),
          error = function(e) warmup_end %m-% months(21L)
        )
        scoring_end <- tryCatch(
          as.Date(
            sim$warmup$scoring_end %||%
              format(warmup_end %m-% months(1L), "%Y-%m-%d")
          ),
          error = function(e) warmup_end %m-% months(1L)
        )
        score_temperature <- as.numeric(sim$warmup$score_temperature %||% 1.0)

        warmup_panels_list <- build_hybrid_warmup_panels(
          sim_monthly_path = sim_output_path,
          surv_panel = surv_panel_warmup,
          warmup_start = warmup_start,
          warmup_end = warmup_end,
          n_rep = n_replicates,
          scoring_start = scoring_start,
          scoring_end = scoring_end,
          score_temperature = score_temperature,
          tu_map = frame %>%
            dplyr::select(gid, country, tu) %>%
            dplyr::distinct()
        )
        message(sprintf(
          "[utils_simulations] warmup panels built: %d replicates, each covering %s – %s",
          length(warmup_panels_list),
          warmup_start,
          warmup_end
        ))
      },
      error = function(e) {
        warning(
          "[utils_simulations] warmup panel build failed: ",
          conditionMessage(e),
          " — falling back to shared panel-based warmup for all replicates."
        )
        warmup_panels_list <<- NULL
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Survival panel for bpfilter warmup (filter.enable = TRUE)
  #
  # The projection climate panel does not carry an `outbreak` column, but
  # build_obs_for_tu requires one to construct the Bernoulli observation table
  # for the block particle filter. When bpfilter is enabled we load the
  # survival dataset (which does have `outbreak`) and attach it to `art` so
  # the caller can left-join it onto panel_tu_warmup before calling
  # build_spatpomp_tu_c.
  # ---------------------------------------------------------------------------
  surv_panel <- NULL
  if (isTRUE(sim$filter$enable)) {
    baseline <- opt_occ$baseline %||% opt$baseline %||% "1990_2019"
    surv_path <- tryCatch(
      make_surv_dat_file(list(baseline = baseline, episodes = opt$episodes)),
      error = function(e) NULL
    )
    if (!is.null(surv_path) && file.exists(surv_path)) {
      message(
        "[utils_simulations] bpfilter: reading survival panel for warmup obs: ",
        surv_path
      )
      surv_panel <- readRDS(surv_path)
    } else {
      warning(
        "[utils_simulations] bpfilter: survival panel not found at '",
        surv_path %||% "NULL",
        "' — 'outbreak' column will be absent from warmup panel."
      )
    }
  }

  list(
    occ_fit = occ_fit,
    dur_fit = dur_fit,
    mon_fit = mon_fit,
    panel = panel,
    frame = frame,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    enso_scenarios = .enso_chunk_slice(enso_scenarios, enso_chunk),
    warmup_panels_list = warmup_panels_list,
    surv_panel = surv_panel,
    cb_meta_occ = cb_meta_occ,
    cb_meta_dur = cb_meta_dur,
    cb_meta_mon = cb_meta_mon
  )
}

# ---------------------------------------------------------------------------
# .load_artefacts_grouped — country-grouping artefact loader
#
# Loads one occ (+ optional dur/mon) fit per group from
# `simulator.grouping.fits_dir`, routes
# each simulated country to its group's draws bundle, and returns per-group
# bundles plus a country→group map. The location frame gains a `group` column.
# Random effects in grouped fits nest country/ADM1 only (region dropped), so
# downstream LP evaluation must pass drop_region = TRUE.
#
# Phased scope: the advanced ENSO-scenario (counterfactual / C3S plume),
# hybrid-warmup, and bpfilter paths are not yet wired for grouped mode; they
# error here with a clear message rather than silently mis-routing.
# ---------------------------------------------------------------------------
.load_artefacts_grouped <- function(opt, enso_sim_months_start = NULL) {
  sim <- opt$simulator
  grp <- sim$grouping

  if (isTRUE(sim$warmup$use_sim_trajectories)) {
    stop(
      "[load_artefacts] grouped mode does not yet support hybrid warmup ",
      "(simulator.warmup.use_sim_trajectories)."
    )
  }
  surv_panel <- NULL
  if (isTRUE(sim$filter$enable)) {
    baseline <- opt$baseline
    surv_path <- tryCatch(
      make_surv_dat_file(list(baseline = baseline, episodes = opt$episodes)),
      error = function(e) NULL
    )
    if (!is.null(surv_path) && file.exists(surv_path)) {
      message(
        "[utils_simulations] bpfilter (grouped): reading survival panel: ",
        surv_path
      )
      surv_panel <- readRDS(surv_path)
    } else {
      warning(
        "[utils_simulations] bpfilter (grouped): survival panel not found at '",
        surv_path %||% "NULL",
        "' — 'outbreak' column will be absent from warmup panel."
      )
    }
  }

  source <- grp$source %||% stop("simulator.grouping.source is required")
  fits_dir <- grp$fits_dir %||% stop("simulator.grouping.fits_dir is required")
  if (!dir.exists(fits_dir)) {
    stop("[load_artefacts] simulator.grouping.fits_dir not found: ", fits_dir)
  }

  message("[utils_simulations] reading climate panel:  ", sim$climate_panel)
  panel <- readRDS(sim$climate_panel)

  # ---------------------------------------------------------------------------
  # Per-replicate ENSO forecast scenarios (group-agnostic; same block as the
  # non-grouped loader). The ENSO ensemble depends only on the panel, config,
  # and n_replicates — it is independent of group routing. Each group's
  # cross-basis meta is applied to these shared trajectories at forecast time.
  # Three branches: counterfactual analog, counterfactual lenssen, factual C3S.
  # ---------------------------------------------------------------------------
  n_replicates <- as.integer(sim$n_replicates %||% 100L)
  enso_chunk <- .enso_chunk_cols(opt, n_replicates) # chunked: sample the TU's total, keep ours
  n_enso <- enso_chunk$n_build %||% n_replicates
  enso_scenarios <- NULL
  ctf_scenario <- sim$enso$counterfactual$scenario
  forecast_source <- sim$enso$forecast_source %||% "cfsv2"

  if (
    !is.null(ctf_scenario) ||
      !is.null(sim$enso$c3s_plume) ||
      identical(forecast_source, "cfsv2")
  ) {
    tryCatch(
      {
        enso_source <- sim$enso$source %||% "oni"
        historical_oni <- fetch_enso_index(
          source = enso_source,
          cache_dir = "generated_data/",
          refresh = FALSE
        )
        # See the non-grouped loader for the rationale: extend the ENSO scenario
        # window back to the smoothing handoff (enso_sim_months_start) when
        # supplied, so the per-replicate override matrices cover the full
        # simulated horizon rather than being silently recycled downstream.
        horizon_start <- as.Date(sim$horizon$start)
        enso_start <- horizon_start
        if (!is.null(enso_sim_months_start)) {
          enso_start <- min(horizon_start, as.Date(enso_sim_months_start))
          if (enso_start < horizon_start) {
            message(sprintf(
              "[utils_simulations] ENSO scenario window extended back to %s (< sim horizon %s) to cover the smoothing handoff (grouped).",
              enso_start,
              horizon_start
            ))
          }
        }
        sim_months <- sort(unique(panel$t_start[panel$t_start >= enso_start]))

        if (!is.null(ctf_scenario)) {
          message(
            "[utils_simulations] ENSO counterfactual scenario (grouped): ",
            ctf_scenario
          )

          if (identical(ctf_scenario, "analog")) {
            analog_start <- as.Date(
              sim$enso$counterfactual$analog_start %||%
                stop(
                  "enso.counterfactual.analog_start is required for the analog scenario."
                )
            )
            analog_end <- as.Date(
              sim$enso$counterfactual$analog_end %||%
                stop(
                  "enso.counterfactual.analog_end is required for the analog scenario."
                )
            )
            enso_scenarios <- sample_enso_analog_scenarios(
              historical_oni = historical_oni,
              sim_months = sim_months,
              analog_start = analog_start,
              analog_end = analog_end,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              n_rep = n_enso
            )
            message(sprintf(
              "[utils_simulations] analog ENSO (grouped): %d lags, %d months, %d replicates",
              length(enso_scenarios),
              if (length(enso_scenarios) > 0L) {
                nrow(enso_scenarios[[1L]])
              } else {
                0L
              },
              if (length(enso_scenarios) > 0L) {
                ncol(enso_scenarios[[1L]])
              } else {
                0L
              }
            ))
          } else if (identical(ctf_scenario, "lenssen")) {
            # One resolution of the counterfactual block (utils/enso.R::
            # lenssen_counterfactual_args()): cache path, analog
            # anchor (lenssen_horizon_start, NULL / "" = horizon_start), offsets
            # and analog set (`%||% "all"` there is the effective default).
            lcf <- lenssen_counterfactual_args(
              sim$enso$counterfactual,
              horizon_start
            )
            enso_horizon <- lcf$horizon_start
            if (enso_horizon != horizon_start) {
              message(sprintf(
                "[utils_simulations] Lenssen ENSO horizon (grouped): %s (overrides sim horizon %s)",
                enso_horizon,
                horizon_start
              ))
            }
            enso_scenarios <- sample_enso_lenssen_scenarios(
              cache_path = lcf$cache_path,
              sim_months = sim_months,
              n_rep = n_enso,
              historical_oni = historical_oni,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              horizon_start = enso_horizon,
              mean_offset = lcf$mean_offset,
              log_var_scaling = lcf$log_var_scaling,
              seed = as.integer(sim$seed %||% 20260603L),
              analog_set = lcf$analog_set
            )
            message(sprintf(
              "[utils_simulations] Lenssen ENSO (grouped): %d lags, %d months × %d replicates (analog_set=%s)",
              length(enso_scenarios),
              if (length(enso_scenarios) > 0L) {
                nrow(enso_scenarios[[1L]])
              } else {
                0L
              },
              if (length(enso_scenarios) > 0L) {
                ncol(enso_scenarios[[1L]])
              } else {
                0L
              },
              as.character(sim$enso$counterfactual$analog_set %||% "all")
            ))
          } else {
            stop(
              "[utils_simulations] Unknown enso.counterfactual.scenario: '",
              ctf_scenario,
              "'. Must be 'analog' or 'lenssen'."
            )
          }
        } else {
          # Factual branch: forecast_source = cfsv2 (default) | c3s
          historical_nino34 <- tryCatch(
            fetch_enso_index(
              source = "nino34",
              cache_dir = "generated_data/",
              refresh = FALSE
            ),
            error = function(e) {
              message(
                "[utils_simulations] nino34 cache not found — using ONI as Niño 3.4 proxy."
              )
              NULL
            }
          )

          if (identical(forecast_source, "cfsv2")) {
            message(
              "[utils_simulations] building per-replicate ENSO trajectories from CFSv2 ensemble (grouped)"
            )
            cfsv2_data <- fetch_enso_cfsv2_ensemble(
              url = sim$enso$cfsv2_url %||%
                "https://www.cpc.ncep.noaa.gov/products/CFSv2/dataInd3/nino34Mon.nc",
              cache_dir = sim$enso$cache_dir %||% "generated_data/",
              init = sim$enso$cfsv2_init %||% "latest"
            )
            message(sprintf(
              "[utils_simulations] CFSv2 ensemble: %d rows, target dates %s – %s, %d members",
              nrow(cfsv2_data),
              min(cfsv2_data$target_date),
              max(cfsv2_data$target_date),
              length(unique(cfsv2_data$member))
            ))
            enso_scenarios <- sample_enso_cfsv2_trajectories(
              cfsv2_data = cfsv2_data,
              n_rep = n_enso,
              historical_oni = historical_oni,
              historical_nino34 = historical_nino34,
              sim_months = sim_months,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              rebase = sim$enso$cfsv2_rebase %||% "seam",
              sd_correction = isTRUE(sim$enso$cfsv2_sd_correction),
              target_peak = sim$enso$cfsv2_target_peak %||% NULL,
              calibration = if (nzchar(sim$enso$cfsv2_calibration %||% "")) {
                sim$enso$cfsv2_calibration
              } else {
                NULL
              },
              calibration_target = sim$enso$cfsv2_calibration_target %||%
                "ersst",
              seed = as.integer(sim$seed %||% 20260603L),
              # past the CFSv2 span: neutral (default) | ar1 | analog (simulator.enso.continuation)
              beyond_cfsv2 = (sim$enso$continuation %||% list())$method %||%
                "neutral",
              ar1_phi = (sim$enso$continuation %||% list())$phi,
              ar1_fit_window = (sim$enso$continuation %||%
                list())$fit_window %||%
                "1991_2020",
              legacy_seam = isTRUE(
                (sim$enso$continuation %||% list())$legacy_seam
              ),
              analog_library = (sim$enso$continuation %||% list())[[
                "library"
              ]] %||%
                "generated_data/enso_continuation_library.qs",
              analog_k = (sim$enso$continuation %||% list())[["k"]] %||% 30L,
              analog_blend_months = (sim$enso$continuation %||% list())[[
                "blend_months"
              ]] %||%
                3,
              analog_min_sep = (sim$enso$continuation %||% list())[[
                "min_sep"
              ]] %||%
                24L,
              analog_max_per_model = (sim$enso$continuation %||% list())[[
                "max_per_model"
              ]] %||%
                Inf,
              analog_exclude_models = unlist((sim$enso$continuation %||%
                list())[["exclude_models"]]),
              # analog seam: last_cfsv2 (default) | last_full_oni
              analog_seam = (sim$enso$continuation %||% list())[["seam"]] %||%
                "last_cfsv2"
            )
          } else {
            message(
              "[utils_simulations] building per-replicate ENSO forecast trajectories (grouped, C3S plume)"
            )
            plume_data <- read_enso_c3s_plume(sim$enso$c3s_plume)
            message(sprintf(
              "[utils_simulations] C3S plume: %d rows, dates %s – %s, quantile levels %s",
              nrow(plume_data),
              min(plume_data$date),
              max(plume_data$date),
              paste(sort(unique(plume_data$q_level)), collapse = ", ")
            ))
            enso_scenarios <- sample_enso_forecast_trajectories(
              plume_data = plume_data,
              n_rep = n_enso,
              historical_oni = historical_oni,
              historical_nino34 = historical_nino34,
              sim_months = sim_months,
              max_lag = 12L,
              max_lead = panel_max_enso_lead(panel),
              seed = as.integer(sim$seed %||% 20260603L)
            )
          }
          message(sprintf(
            "[utils_simulations] ENSO scenarios built (grouped): %d lags, %d months × %d replicates",
            length(enso_scenarios),
            if (length(enso_scenarios) > 0L) nrow(enso_scenarios[[1]]) else 0L,
            if (length(enso_scenarios) > 0L) ncol(enso_scenarios[[1]]) else 0L
          ))
        }
      },
      error = function(e) {
        # See the non-grouped loader above: fail loudly instead of silently
        # falling back to unoverridden (possibly NA) panel ENSO.
        stop(
          "[utils_simulations] ENSO scenario build failed (grouped): ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  reg <- load_grouping(source)
  country_to_group <- stats::setNames(reg$group, reg$iso3)

  # scope.countries (set per generated forecast/counterfactual config -- see
  # make_forecast_configs.R) determines which group(s) this run actually
  # touches; resolved once here and reused below both to pick a representative
  # deployed fit (for pop / base-rate path resolution, next) and to filter the
  # frame (Apply scope.countries filter, below).
  iso3_req <- if (!identical(sim$scope$countries %||% "all", "all")) {
    trimws(strsplit(paste(sim$scope$countries, collapse = ","), ",")[[1L]])
  } else {
    character(0)
  }

  # Per-track base_rate anchor tables (utils/anchor.R) and the population
  # sidecars must be resolved from a deployed FIT's own training opt, not the
  # simulator's top-level `opt`. Each table is global per track (pooled across
  # groups), so any one deployed fit resolves the path -- EXCEPT that not every
  # group necessarily has a fit deployed. Prefer the group(s) implied by
  # scope.countries so we peek a fit that's actually relevant to this run;
  # fall back to every group in the registry, in order, until one resolves.
  #
  # NB: a candidate only counts as resolved if attr(fit, "opt") is non-NULL.
  # .fit_opt_from_attr() (used elsewhere) papers over a missing attribute
  # with list() -- fine for its callers, which just want a best-effort opt,
  # but here a missing opt would silently resolve default paths. Read attr()
  # directly so a fit with no recorded opt is treated as a failed candidate.
  candidate_groups <- if (length(iso3_req) > 0L) {
    sort(unique(stats::na.omit(country_to_group[iso3_req])))
  } else {
    character(0)
  }
  if (length(candidate_groups) == 0L) {
    candidate_groups <- sort(unique(reg$group))
  }

  # duration/monthly are opt-in (grouping.duration_enable / monthly_enable,
  # resolved further below), and an occ fit predating the anchor simply has
  # `offset_baserate: false`, which resolve_base_rate_paths() skips -- so for the
  # anchors this is best-effort: a track with no discoverable fit contributes
  # NULL and is dropped. The hard anchor error lives in resolve_base_rate_paths(),
  # and only fires for a fit that WAS trained with the anchor but whose sidecar
  # is missing. The occurrence representative is required (below): it also
  # resolves the population sidecars.
  repr_anchor_fit_opts <- list()
  for (.tr in c("occ", "dur", "mon")) {
    for (.g in candidate_groups) {
      .o <- tryCatch(
        attr(
          readRDS(.discover_group_fit(
            fits_dir,
            source,
            .g,
            .tr,
            episode_key_of(opt),
            national_zero_fill_of(opt),
            national_zero_threshold_of(opt),
            occ_panel = occ_panel_spec_of(opt)
          )),
          "opt"
        ),
        error = function(e) NULL
      )
      if (!is.null(.o)) {
        repr_anchor_fit_opts[[.tr]] <- .o
        break
      }
    }
  }
  repr_occ_fit_opt <- repr_anchor_fit_opts$occ
  if (is.null(repr_occ_fit_opt)) {
    stop(
      "[load_artefacts] could not find any deployed occurrence fit with a ",
      "recorded 'opt' attribute in '",
      fits_dir,
      "' among candidate group(s): ",
      paste(candidate_groups, collapse = ", "),
      " -- needed to resolve the population sidecar paths. Pass ",
      "simulator.location_frame.gid_pop_path / wpp_path explicitly, or ",
      "check that the fits_dir holds this run's per-group occurrence fits."
    )
  }

  # ADM1-ADM2 overlap filter's authoritative admin_level/ADM1 columns — reuse
  # surv_panel if the bpfilter warmup already loaded it (lines above);
  # otherwise read it here. Loaded ahead of build_location_frame() so the same
  # object also supplies `admin_level_data`,
  # and is reused again below for drop_overlapping_adm1_from_frame().
  .surv_overlap <- surv_panel
  if (is.null(.surv_overlap)) {
    .surv_overlap <- tryCatch(
      readRDS(make_surv_dat_file(list(
        baseline = opt$baseline,
        episodes = opt$episodes
      ))),
      error = function(e) NULL
    )
  }

  lf_opts <- sim$location_frame %||% list()
  .fpop <- sim_frame_pop_args(sim, repr_occ_fit_opt)
  .fpop$pop_window <- sim_pop_window(sim)
  frame <- build_location_frame(
    panel,
    hmm_clusters_path = lf_opts$hmm_clusters_path,
    outbreak_metadata_path = lf_opts$outbreak_metadata_path,
    size_data_path = lf_opts$size_data_path,
    # NULL (hence fully inert) unless a deployed fit was actually made with
    # `<track>.offset_baserate: true` -- see resolve_base_rate_paths().
    base_rate_table_paths = resolve_base_rate_paths(
      repr_anchor_fit_opts,
      lf_opts
    ),
    admin_level_data = .surv_overlap,
    # Population arguments in one place (simulator_core.R::sim_frame_pop_args():
    # location_frame.gid_pop_path / wpp_path / pop_year).
    gid_pop_path = .fpop$gid_pop_path,
    wpp_path = .fpop$wpp_path,
    pop_year = .fpop$pop_year,
    pop_window = .fpop$pop_window
  )
  # gids whose min_seas was defaulted, checked on the final frame.
  .ms_defaulted <- attr(frame, "min_seas_defaulted")

  # Apply scope.countries filter (same as non-grouped path; iso3_req resolved
  # once, above).
  if (length(iso3_req) > 0L) {
    frame <- frame %>%
      dplyr::filter(stringr::str_extract(adm1_unit, "^[A-Z]{3}") %in% iso3_req)
    if (nrow(frame) == 0L) {
      stop(sprintf(
        "scope.countries filter ('%s') left 0 locations.",
        paste(iso3_req, collapse = ", ")
      ))
    }
  }

  # Attach group column; validate every simulated country maps to a group.
  frame$group <- unname(country_to_group[frame$country])
  unmapped <- unique(frame$country[is.na(frame$group)])
  if (length(unmapped) > 0L) {
    stop(sprintf(
      paste0(
        "[load_artefacts] %d simulated country(ies) absent from grouping registry ",
        "'%s': %s. Add them to analysis/groupings/%s.csv or restrict scope.countries."
      ),
      length(unmapped),
      source,
      paste(unmapped, collapse = ", "),
      source
    ))
  }

  groups <- sort(unique(frame$group))
  message(sprintf(
    "[utils_simulations] grouped mode: %d group(s) [%s], %d locations",
    length(groups),
    paste(groups, collapse = ", "),
    nrow(frame)
  ))

  max_draws <- as.integer(sim$max_draws %||% 200L)

  # Duration track (optional; opt-in via simulator.grouping.duration_enable
  # so batches with no per-group duration fits deployed yet keep loading
  # unchanged). When enabled, every group in the registry MUST have a
  # discoverable duration fit -- same hard-error convention as occ, since
  # this is an explicit opt-in rather than a best-effort sidecar like cb_meta.
  duration_enabled <- isTRUE(grp$duration_enable)

  # Monthly case-count track: same opt-in convention, plus the hard duration
  # requirement (`active`/`mso` are duration-model latent states).
  monthly_enabled <- isTRUE(grp$monthly_enable)
  if (monthly_enabled && !duration_enabled) {
    stop(
      "simulator.grouping.monthly_enable requires ",
      "simulator.grouping.duration_enable -- the monthly case-count kernel reads ",
      "`active`/`mso`, which are duration-model latent states."
    )
  }

  draws_occ_by_group <- list()
  draws_dur_by_group <- if (duration_enabled) list() else NULL
  draws_mon_by_group <- if (monthly_enabled) list() else NULL
  cb_meta_occ_by_group <- list()
  cb_meta_dur_by_group <- if (duration_enabled) list() else NULL
  cb_meta_mon_by_group <- if (monthly_enabled) list() else NULL
  required_lags <- character(0)

  for (g in groups) {
    occ_path <- .discover_group_fit(
      fits_dir,
      source,
      g,
      "occ",
      episode_key_of(opt),
      national_zero_fill_of(opt),
      national_zero_threshold_of(opt),
      occ_panel = occ_panel_spec_of(opt)
    )
    message(sprintf(
      "[utils_simulations] group '%s' occ:  %s",
      g,
      basename(occ_path)
    ))

    occ_fit <- readRDS(occ_path)

    opt_occ <- .fit_opt_from_attr(occ_fit)
    if (is.null(opt_occ$occurrence)) {
      be <- "brms" # every supported occurrence fit is a brmsfit
      opt_occ$occurrence <- list(backend = be)
    }

    dur_fit <- NULL
    opt_dur <- NULL
    if (duration_enabled) {
      dur_path <- .discover_group_fit(
        fits_dir,
        source,
        g,
        "dur",
        episode_key_of(opt)
      )
      message(sprintf(
        "[utils_simulations] group '%s' dur:  %s",
        g,
        basename(dur_path)
      ))
      dur_fit <- readRDS(dur_path)
      opt_dur <- .fit_opt_from_attr(dur_fit)
      if (is.null(opt_dur$duration)) opt_dur$duration <- list(backend = "brms")
    }

    mon_fit <- NULL
    opt_mon <- NULL
    if (monthly_enabled) {
      mon_path <- .discover_group_fit(
        fits_dir,
        source,
        g,
        "mon",
        episode_key_of(opt)
      )
      message(sprintf(
        "[utils_simulations] group '%s' mon:  %s",
        g,
        basename(mon_path)
      ))
      mon_fit <- readRDS(mon_path)
      opt_mon <- .fit_opt_from_attr(mon_fit)
      if (is.null(opt_mon$monthly)) opt_mon$monthly <- list(backend = "brms_cb")
    }

    .check_fit_pair_consistency(
      opt_occ,
      label = sprintf("group '%s'", g),
      opt_dur = opt_dur,
      opt_mon = opt_mon
    )

    # Layer 4 (artefact chain, see fit_pop_stamp_path() in utils.R) -- same check
    # as the single-fit path above, per group (monthly fits only).
    if (monthly_enabled) {
      assert_fit_pop_provenance(
        opt_mon,
        mon_path,
        what = sprintf("the monthly fit for group '%s'", g)
      )
    }

    # Episode definition, per group (see the ungrouped loader above).
    sim_episode_check(
      fits = list(occ = occ_fit, dur = dur_fit, mon = mon_fit),
      preps = list(
        surv = make_surv_dat_file(list(
          baseline = opt$baseline,
          episodes = opt$episodes
        )),
        size_data = (opt$simulator$location_frame %||%
          list())$size_data_path %||%
          make_size_dat_file(opt),
        climate_panel = sim$climate_panel
      ),
      what = sprintf("load_artefacts, group '%s'", g),
      episode_key = episode_key_of(opt)
    )

    # Per-group cross-basis meta sidecars (derived from each group's fit path).
    # Duration IS included: `04_fit_dur_brms_cb.R` produces `_cb` duration fits
    # whose FE carry `enso_cb<k>`, and the deployed per-group best-lag duration
    # fits are exactly those. Omitting `cb_dur` here left the duration LP with
    # no `enso_cb*` covariates, and evaluate_lp() silently skips a FE whose
    # covariate is absent -- so the duration track's whole ENSO response was
    # dropped at replay while the fit reported a significant effect.
    cb_occ <- if (
      file.exists(sub("_fit_base-", "_cb_meta_base-", occ_path, fixed = TRUE))
    ) {
      readRDS(sub("_fit_base-", "_cb_meta_base-", occ_path, fixed = TRUE))
    } else {
      NULL
    }
    cb_dur <- if (
      duration_enabled &&
        file.exists(sub("_fit_base-", "_cb_meta_base-", dur_path, fixed = TRUE))
    ) {
      readRDS(sub("_fit_base-", "_cb_meta_base-", dur_path, fixed = TRUE))
    } else {
      NULL
    }
    cb_mon <- if (
      monthly_enabled &&
        file.exists(sub("_fit_base-", "_cb_meta_base-", mon_path, fixed = TRUE))
    ) {
      readRDS(sub("_fit_base-", "_cb_meta_base-", mon_path, fixed = TRUE))
    } else {
      NULL
    }

    d_occ <- extract_occurrence_draws(occ_fit, opt_occ, max_draws = max_draws)
    d_dur <- if (duration_enabled) {
      extract_duration_draws(dur_fit, opt_dur, max_draws = max_draws)
    } else {
      NULL
    }
    d_mon <- if (monthly_enabled) {
      extract_monthly_draws(mon_fit, opt_mon, max_draws = max_draws)
    } else {
      NULL
    }

    draws_occ_by_group[[g]] <- d_occ
    if (duration_enabled) {
      draws_dur_by_group[[g]] <- d_dur
    }
    if (monthly_enabled) {
      draws_mon_by_group[[g]] <- d_mon
    }
    cb_meta_occ_by_group[[g]] <- cb_occ
    if (duration_enabled) {
      cb_meta_dur_by_group[[g]] <- cb_dur
    }
    if (monthly_enabled) {
      cb_meta_mon_by_group[[g]] <- cb_mon
    }
    required_lags <- c(required_lags, .lags_required_from_draws(d_occ, cb_occ))
    # Duration contributes to the panel lag union too: its best lag is selected
    # independently of occ/mon (e.g. Horn is occ=11, mon=11, dur=8), so a panel
    # built from the other tracks' lags alone can be missing
    # the column the duration cross-basis needs.
    if (duration_enabled) {
      required_lags <- c(
        required_lags,
        .lags_required_from_draws(d_dur, cb_dur)
      )
    }
    if (monthly_enabled) {
      required_lags <- c(
        required_lags,
        .lags_required_from_draws(d_mon, cb_mon)
      )
    }
  }

  # Panel lag-union validation: the climate panel must carry every lag column
  # any group consumes (different groups may use different enso/climate lags).
  required_lags <- unique(required_lags)
  missing_lags <- setdiff(required_lags, names(panel))
  if (length(missing_lags) > 0L) {
    stop(sprintf(
      paste0(
        "[load_artefacts] climate panel lacks %d lag column(s) required across ",
        "groups: %s. Re-run 03_prepare_sim_panel.R so the panel carries the union ",
        "of every group's lags."
      ),
      length(missing_lags),
      paste(missing_lags, collapse = ", ")
    ))
  }

  message(sprintf(
    "[utils_simulations] grouped lag union (%d cols) present in panel: %s",
    length(required_lags),
    if (length(required_lags)) {
      paste(required_lags, collapse = ", ")
    } else {
      "(none)"
    }
  ))

  # Zero-filled climate gate, over every group's occ/dur/mon fit.
  .assert_no_zero_filled_climate(
    .draws_fe_names(draws_occ_by_group, draws_dur_by_group, draws_mon_by_group),
    panel,
    window = .sim_panel_window(sim, enso_sim_months_start),
    label = "load_artefacts (grouped)"
  )

  # ADM1-ADM2 overlap filter — mirror the smoothing pipeline (smoothing_core.R)
  # so the forecast frame matches the smoothed frame it consumes. `.surv_overlap`
  # was already loaded above (also supplied build_location_frame()'s
  # `admin_level_data`), so no second read here.
  frame <- drop_overlapping_adm1_from_frame(
    frame,
    .surv_overlap,
    "utils_simulations grouped"
  )
  # TTLO level "none": see the ungrouped loader.
  if (.draws_has_ttlo_none(draws_occ_by_group)) {
    frame$first_onset <- frame_first_onset(.surv_overlap, frame$gid)
  }
  if (is.null(surv_panel)) {
    rm(.surv_overlap)
  }

  # Per-group frame-vs-draws check (region-dropped keys). Each group's frame
  # subset is validated against that group's occurrence draws; restrict policy
  # may drop unseen locations. Reassemble the checked subsets into one frame.
  frame_checked <- do.call(
    rbind,
    lapply(groups, function(g) {
      fg <- frame[frame$group == g, , drop = FALSE]
      check_frame_vs_draws(
        fg,
        draws_occ_by_group[[g]],
        sim$scope,
        drop_region = TRUE,
        draws_mon = if (monthly_enabled) draws_mon_by_group[[g]] else NULL
      )
    })
  )
  frame_checked <- tibble::as_tibble(frame_checked)
  # The season onset gate needs every simulated gid's own season start.
  check_onset_eligibility_frame(
    frame_checked,
    onset_eligibility_of(opt$simulator),
    "load_artefacts",
    defaulted = .ms_defaulted
  )

  list(
    grouped = TRUE,
    occ_fit = NULL,
    panel = panel,
    frame = frame_checked,
    draws_occ = draws_occ_by_group[[groups[[1L]]]], # representative
    draws_dur = if (duration_enabled) {
      draws_dur_by_group[[groups[[1L]]]]
    } else {
      NULL
    },
    draws_mon = if (monthly_enabled) {
      draws_mon_by_group[[groups[[1L]]]]
    } else {
      NULL
    },
    draws_occ_by_group = draws_occ_by_group,
    draws_dur_by_group = draws_dur_by_group, # NULL unless grouping.duration_enable
    draws_mon_by_group = draws_mon_by_group, # NULL unless grouping.monthly_enable
    cb_meta_occ_by_group = cb_meta_occ_by_group,
    cb_meta_dur_by_group = cb_meta_dur_by_group, # NULL unless grouping.duration_enable
    cb_meta_mon_by_group = cb_meta_mon_by_group,
    country_to_group = country_to_group,
    enso_scenarios = .enso_chunk_slice(enso_scenarios, enso_chunk),
    warmup_panels_list = NULL,
    surv_panel = surv_panel,
    cb_meta_occ = cb_meta_occ_by_group[[groups[[1L]]]],
    cb_meta_dur = if (duration_enabled) {
      cb_meta_dur_by_group[[groups[[1L]]]]
    } else {
      NULL
    },
    cb_meta_mon = if (monthly_enabled) {
      cb_meta_mon_by_group[[groups[[1L]]]]
    } else {
      NULL
    }
  )
}

# ---------------------------------------------------------------------------
# check_frame_vs_draws — admin-level alignment guards
# ---------------------------------------------------------------------------

#' Validate and (if needed) filter the location frame against the occurrence
#' draws bundle. Three checks:
#'
#' 1. ADM1-level gids (`gid == adm1_unit`) under `scope.admin_level = "ADM2"`:
#'    warn that they are simulated at ADM1 resolution.
#' 2. Gids with no column in `draws$re_intercept$adm2`: the kernel falls back
#'    to the parent REs with adm2 RE = 0. `scope.unseen_adm1 = "restrict"`
#'    (default) drops them (error if none remain); `"marginalise"` warns and
#'    keeps them.
#' 3. ADM1 units absent from the adm1 RE matrix: same zero fallback, warning
#'    only.
#'
#' @param frame tibble from build_location_frame().
#' @param draws_occ draws bundle from extract_occurrence_draws().
#' @param sim_scope list from `opt$simulator$scope`.
#' @param drop_region logical — TRUE for grouped (region-dropped) fits, so the
#'   RE-matrix lookups use `country/ADM1`-only composite keys. Default FALSE.
#' @param draws_mon optional monthly draws bundle; its adm2 (leaf) and adm1
#'   intercept keys are checked against the frame too, warn only, never
#'   restrict. NULL (the default) skips it.
#' @return frame, possibly with rows removed (restrict policy).
check_frame_vs_draws <- function(
  frame,
  draws_occ,
  sim_scope,
  drop_region = FALSE,
  draws_mon = NULL
) {
  admin_level_scope <- sim_scope$admin_level %||% "ADM2"
  unseen_policy <- sim_scope$unseen_adm1 %||% "restrict"

  # ── 1. ADM1-level gids when scope declares ADM2 ───────────────────────────
  # ADM1-level gids are 3-part `ISO_X_1` (e.g. NGA_19_1); adm1_unit is the 2-part
  # `ISO_X` (NGA_19). The prior `gid == adm1_unit` test compared 3-part vs 2-part
  # so it never matched and the guard silently no-opped. Recognise the ADM1 form
  # as adm1_unit + "_1" (matching get_adm1_gid()).
  is_adm1_gid <- !is.na(frame$gid) &
    !is.na(frame$adm1_unit) &
    (as.character(frame$gid) == paste0(frame$adm1_unit, "_1"))
  n_adm1_gids <- sum(is_adm1_gid)
  if (admin_level_scope == "ADM2" && n_adm1_gids > 0L) {
    shown <- head(frame$gid[is_adm1_gid], 5)
    warning(sprintf(
      paste0(
        "[check_frame_vs_draws] scope.admin_level=ADM2 but %d location(s) in the ",
        "frame are ADM1-level (gid == adm1_unit). These are from countries where ",
        "the fit used ADM1 data; they will be simulated at ADM1 resolution. ",
        "Example gids: %s%s"
      ),
      n_adm1_gids,
      paste(shown, collapse = ", "),
      if (n_adm1_gids > 5L) sprintf(" ... (+%d more)", n_adm1_gids - 5L) else ""
    ))
  }

  # ── 2. Gids absent from adm2 RE matrix ───────────────────────────────────
  adm2_cols <- colnames(draws_occ$re_intercept$adm2)
  if (!is.null(adm2_cols) && length(adm2_cols) > 0L) {
    unseen_mask <- !(re_level_key(frame, "adm2", drop_region) %in% adm2_cols)
    n_unseen <- sum(unseen_mask)
    if (n_unseen > 0L) {
      unseen_ctries <- unique(stringr::str_extract(
        frame$adm1_unit[unseen_mask],
        "^[A-Z]{3}"
      ))
      base_msg <- sprintf(
        paste0(
          "[check_frame_vs_draws] %d location(s) in the frame have no adm2 ",
          "intercept RE in the draws (gid not seen at fit time). The kernel ",
          "falls back to parent-level REs only (adm2 RE = 0). Affected ",
          "countries: %s%s. unseen_adm1='%s'"
        ),
        n_unseen,
        paste(head(unseen_ctries, 5), collapse = ", "),
        if (length(unseen_ctries) > 5L) {
          sprintf(" ... (+%d more)", length(unseen_ctries) - 5L)
        } else {
          ""
        },
        unseen_policy
      )
      if (unseen_policy == "restrict") {
        message(base_msg, " — dropping from frame.")
        frame <- frame[!unseen_mask, , drop = FALSE]
        if (nrow(frame) == 0L) {
          stop(paste0(
            "[check_frame_vs_draws] unseen_adm1='restrict' removed ",
            "all locations from the frame. Check that scope.admin_level ",
            "and the fit's training data are compatible."
          ))
        }
        message(sprintf(
          "[check_frame_vs_draws] frame reduced to %d location(s) after restrict.",
          nrow(frame)
        ))
      } else {
        warning(base_msg, " — proceeding with marginalised adm2 RE.")
      }
    }
  }

  # ── 3. ADM1 units absent from adm1 RE matrix ─────────────────────────────
  adm1_cols <- colnames(draws_occ$re_intercept$adm1)
  if (!is.null(adm1_cols) && length(adm1_cols) > 0L) {
    unseen_adm1_mask <- !(re_level_key(frame, "adm1", drop_region) %in%
      adm1_cols)
    n_unseen_adm1 <- sum(unseen_adm1_mask)
    if (n_unseen_adm1 > 0L) {
      unseen_units <- unique(frame$adm1_unit[unseen_adm1_mask])
      warning(sprintf(
        paste0(
          "[check_frame_vs_draws] %d location(s) span %d ADM1 unit(s) absent from ",
          "the adm1 intercept RE matrix — adm1 RE will be zero for these (falls ",
          "back to country + region REs only). ADM1 units: %s%s"
        ),
        n_unseen_adm1,
        length(unseen_units),
        paste(head(unseen_units, 5), collapse = ", "),
        if (length(unseen_units) > 5L) {
          sprintf(" ... (+%d more)", length(unseen_units) - 5L)
        } else {
          ""
        }
      ))
    }
  }

  # ── 4. Magnitude track (monthly): warn-only leaf/ADM1 coverage ───────────────
  # The occurrence restrict above defines the frame; here we only REPORT how
  # much of that frame the monthly fit never saw at the leaf (adm2) and ADM1
  # tiers, i.e. where its RE contribution is a silent 0.
  for (.tr in "monthly") {
    bnd <- draws_mon
    if (is.null(bnd) || is.null(bnd$re_intercept)) {
      next
    }
    for (.lv in c("adm2", "adm1")) {
      cn <- colnames(bnd$re_intercept[[.lv]])
      if (is.null(cn) || !length(cn)) {
        next
      }
      keys <- re_level_key(frame, .lv, drop_region)
      alt <- re_level_key(frame, .lv, !drop_region)
      miss <- !(keys %in% cn) & !(alt %in% cn)
      if (any(miss)) {
        warning(
          sprintf(
            paste0(
              "[check_frame_vs_draws] %s: %d of %d frame location(s) have no %s ",
              "intercept RE in the %s draws -- their %s RE is 0 at replay. e.g. %s"
            ),
            .tr,
            sum(miss),
            nrow(frame),
            .lv,
            .tr,
            .lv,
            paste(utils::head(frame$gid[miss], 5), collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }

  frame
}
