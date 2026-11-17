# 07b_simulator_v1.R — v1 pure-R cholera outbreak simulator
#
# Usage:
#   Rscript analysis/07_simulator_run.R \
#     -C analysis/configs/sim_smoketest/sim_smoke.yml
#   Rscript analysis/07b_simulator_v1.R -C <cfg> --dry-run
#   Rscript analysis/07b_simulator_v1.R -C <cfg> \
#     --task generated_data/sim/<tag>/tasks/task_0001.yml
#
# Details:
#   Backend "v1" of 07_simulator_run.R. Parses YAML + CLI, validates and loads
#   the fitted artefacts, extracts posterior draws (simulator_core.R) and runs
#   the per-TU AR kernel simulate_one_tu(); outputs under outdir/sim/<run_tag>/.
#   Full-run mode loops over replicate x transmission unit; --task runs one
#   SLURM array slot.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/dlnm.R") # build_cb_columns_from_Q() -- required by simulator_core.R's simulate_one_tu() whenever the loaded fit is a DLNM crossbasis (brms_cb)
source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatial_contagion.R") # build_partition_nbr_kernel()/nbr_from_adm1_active()/resolve_use_nbr() -- required by simulator_core.R's simulate_one_tu() whenever simulator.nbr_contagion.enable is set (modules never source each other)
source("analysis/utils/anchor.R")

# 2. CLI options ---------------------------------------------------------------

option_list <- c(
  list(
    make_option(
      c("--task"),
      default = NULL,
      type = "character",
      help = "Path to a single-task YAML (SLURM array slot mode). When set, runs only the assignment in the file."
    )
  ),
  simulator_common_options() # --dry-run / --overwrite / --n_replicates / --workers (utils/config.R)
)

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = c("overwrite"),
  model = NULL # simulator has no track block
)

# Late CLI overrides for nested simulator keys (they live under opt$simulator,
# not at the top level, so make_options_from_config routed them to top-level
# flags). Shared applier lives in utils/config.R.
opt <- apply_simulator_common_cli(opt)
require_simulator_block(opt)

print_options(opt)

# 3. Validation ----------------------------------------------------------------

#' Validate the resolved simulator config: required paths, scope, horizon.
#' Errors loudly on any condition that would make the run unrunnable;
#' warns on softer issues (e.g. workers > 1 on parallel.backend slurm).
validate_simulator_options <- function(opt) {
  sim <- opt$simulator
  stopifnot(!is.null(sim))

  required_paths <- c(
    occurrence_fit = sim$occurrence_fit,
    climate_panel = sim$climate_panel
  )
  for (nm in names(required_paths)) {
    p <- required_paths[[nm]]
    if (is.null(p) || !nzchar(p)) {
      stop(sprintf("simulator.%s is required (got NULL/empty).", nm))
    }
    if (!file.exists(p)) {
      stop(sprintf("simulator.%s not found on disk: %s", nm, p))
    }
  }

  # Horizon: parse YYYY-MM-DD; end must be strictly after start.
  start <- tryCatch(as.Date(sim$horizon$start), error = function(e) NA)
  end <- tryCatch(as.Date(sim$horizon$end), error = function(e) NA)
  if (is.na(start) || is.na(end)) {
    stop(sprintf(
      "simulator.horizon must be YYYY-MM-DD; got start=%s end=%s",
      sim$horizon$start,
      sim$horizon$end
    ))
  }
  if (end <= start) {
    stop(sprintf(
      "simulator.horizon: end (%s) must be after start (%s).",
      end,
      start
    ))
  }

  # Only 'paired' exists; the size track is retired from the simulator.
  if (!identical(sim$posterior_strategy %||% "paired", "paired")) {
    stop(sprintf(
      "simulator.posterior_strategy must be 'paired'; got '%s'.",
      sim$posterior_strategy
    ))
  }
  if (!sim$scope$admin_level %in% c("ADM1", "ADM2")) {
    stop(sprintf(
      "simulator.scope.admin_level must be 'ADM1' or 'ADM2'; got '%s'.",
      sim$scope$admin_level
    ))
  }
  if (!sim$scope$unseen_adm1 %in% c("restrict", "marginalise")) {
    stop(sprintf(
      "simulator.scope.unseen_adm1 must be 'restrict' or 'marginalise'; got '%s'.",
      sim$scope$unseen_adm1
    ))
  }
  if (!sim$parallel$backend %in% c("furrr", "slurm")) {
    stop(sprintf(
      "simulator.parallel.backend must be 'furrr' or 'slurm'; got '%s'.",
      sim$parallel$backend
    ))
  }
  if (!as.integer(sim$n_replicates) >= 1) {
    stop(sprintf(
      "simulator.n_replicates must be >= 1; got %s.",
      sim$n_replicates
    ))
  }
  if (!as.integer(sim$parallel$chunk_size) >= 1) {
    stop(sprintf(
      "simulator.parallel.chunk_size must be >= 1; got %s.",
      sim$parallel$chunk_size
    ))
  }

  # v1 scope guard — per the design note's "v1 scope" section.
  if (!identical(sim$scenario$climate, "factual")) {
    warning(
      sprintf(
        "simulator.scenario.climate = '%s' is a v2 feature; v1 only supports 'factual'. ",
        sim$scenario$climate
      ),
      "Proceeding with v1 historical replay; the panel adapter will error on non-factual data."
    )
  }

  invisible(opt)
}

validate_simulator_options(opt)

# 4. Output paths + run directory ----------------------------------------------

# Derive run_tag from the config filename if not explicitly set.
if (is.null(opt$simulator$run_tag)) {
  opt$simulator$run_tag <- derive_sim_run_tag("sim")
  message("[07_simulate_v1] derived run_tag = ", opt$simulator$run_tag)
}

run_dir <- make_run_file(opt, kind = "sim")
tasks_dir <- make_run_file(opt, kind = "sim_tasks_dir")
parts_dir <- make_run_file(opt, kind = "sim_parts_dir")

if (!isTRUE(opt$dry_run)) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
}

message("[07_simulate_v1] run_dir   = ", run_dir)
message("[07_simulate_v1] tasks_dir = ", tasks_dir)
message("[07_simulate_v1] parts_dir = ", parts_dir)

# 5. Artefact load + draws extraction ------------------------------------------
# (step 2)
#
# In dry-run mode we skip the (potentially expensive) fit reads + extractions.
# Otherwise we load: occurrence fit (+ optional duration / monthly fits) + the
# climate / survival panel, then build the location frame and the draws
# bundles. The result is held in a single environment that both modes
# (full-run and --task) consume.

.fit_opt_from_attr <- function(fit) {
  # The fit-time opt is recorded as an attribute by 04*/05* scripts (best effort);
  # fall back to a minimal opt with the expected backend slots.
  opt_fit <- attr(fit, "opt")
  if (is.null(opt_fit)) {
    opt_fit <- list()
  }
  opt_fit
}

load_artefacts <- function(opt) {
  sim <- opt$simulator

  # Country-grouping mode: delegate entirely to the shared grouped loader from
  # utils/simulations.R.  Sourced into a local env to avoid shadowing 07's own
  # private helpers (check_frame_vs_draws etc.).  The local env's parent is the
  # globalenv so it inherits utils/config.R / utils.R / simulator_core.R symbols.
  if (isTRUE(sim$grouping$enable)) {
    .usim_env <- new.env(parent = globalenv())
    source("analysis/utils/simulations.R", local = .usim_env)
    return(.usim_env$load_artefacts(opt))
  }

  message("[07_simulate_v1] reading occurrence fit: ", sim$occurrence_fit)
  occ_fit <- readRDS(sim$occurrence_fit)
  message("[07_simulate_v1] reading climate panel:  ", sim$climate_panel)
  panel <- readRDS(sim$climate_panel)

  opt_occ <- .fit_opt_from_attr(occ_fit)

  if (is.null(opt_occ$occurrence)) {
    # Every supported occurrence fit is a brmsfit; the extractor dispatches on it.
    be <- "brms"
    opt_occ$occurrence <- list(backend = be)
  }

  # ---------------------------------------------------------------------------
  # Duration fit (optional). simulator.duration_fit is unset for every
  # config without the duration track -- draws_dur / dur_fit stay NULL,
  # and downstream simulate_one_tu() take the has_dur=FALSE legacy path.
  # ---------------------------------------------------------------------------
  has_duration_fit <- !is.null(sim$duration_fit) && nzchar(sim$duration_fit)
  dur_fit <- NULL
  opt_dur <- NULL
  if (has_duration_fit) {
    message("[07_simulate_v1] reading duration fit:   ", sim$duration_fit)
    dur_fit <- readRDS(sim$duration_fit)
    opt_dur <- .fit_opt_from_attr(dur_fit)
    if (is.null(opt_dur$duration)) opt_dur$duration <- list(backend = "brms")
  }

  # ---------------------------------------------------------------------------
  # Monthly case-count fit (`mon` track, optional; the sole case estimand
  # since the size track was retired). When set, the kernel emits a case count
  # at every ACTIVE month with mu driven by months-since-onset. `active`/`mso`
  # are duration-model latent states,
  # so a duration fit is a hard requirement -- there is no invented fallback.
  # (Mirrors utils/simulations.R::load_artefacts(); both copies must stay in
  # step or this runner silently diverges from the shared loader.)
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
    message("[07_simulate_v1] reading monthly fit:    ", sim$monthly_fit)
    mon_fit <- readRDS(sim$monthly_fit)
    opt_mon <- .fit_opt_from_attr(mon_fit)
    if (is.null(opt_mon$monthly)) {
      opt_mon$monthly <- list(backend = "brms_cb")
    }
    # Layer 4 (artefact chain; mirrors utils/simulations.R::load_artefacts()):
    # the monthly fit is the population-offset case-count model, so its
    # companion pop stamp (if any) must match the population currently on disk.
    # require_stamp = FALSE -- a fit predating the stamp keeps loading.
    assert_fit_pop_provenance(
      opt_mon,
      sim$monthly_fit,
      what = "the monthly fit"
    )
  }

  # ---------------------------------------------------------------------------
  # DLNM cross-basis metadata sidecars (optional).
  # Derived from the explicit fit paths by string transform (the fit's opt
  # attribute infers backend as "brms", not "brms_cb", so make_run_file path
  # reconstruction is unreliable). A config override is honoured if present.
  # NULL when the fit carries no cross-basis (zero behaviour change).
  # ---------------------------------------------------------------------------
  .load_cb_meta <- function(fit_path, override) {
    cb_path <- override %||%
      sub("_fit_base-", "_cb_meta_base-", fit_path, fixed = TRUE)
    if (!is.null(cb_path) && file.exists(cb_path)) {
      message("[07_simulate_v1] reading cross-basis meta: ", cb_path)
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
  # best-lag duration fits are `_cb` fits carrying `enso_cb<k>` fixed effects, so
  # this sidecar is not optional-by-construction. Omitting it makes the duration
  # LP silently ENSO-blind (evaluate_lp() skips a FE whose covariate is absent).
  cb_meta_dur <- if (has_duration_fit) {
    .load_cb_meta(sim$duration_fit, sim$duration_cb_meta)
  } else {
    NULL
  }

  # admin_level source for anchor resolution — this loader (unlike
  # utils/simulations.R / smoothing_core.R) does not otherwise read the
  # occurrence survival panel, so this read exists solely to supply
  # `admin_level_data`; freed immediately after build_location_frame().
  .surv_overlap <- tryCatch(
    readRDS(make_surv_dat_file(list(
      baseline = opt_occ$baseline %||% opt$baseline %||% "1990_2019"
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
    # Each anchor resolves off its OWN fit's training opt (the sidecar name
    # carries the estimator + that track's own tokens), never off `opt`.
    # opt_dur/opt_mon are NULL unless those optional fits were configured, and
    # resolve_base_rate_paths() drops a NULL track -- so this is fully inert for
    # every fit made without `<track>.offset_baserate: true`.
    base_rate_table_paths = resolve_base_rate_paths(
      list(occ = opt_occ, dur = opt_dur, mon = opt_mon),
      lf_opts
    ),
    admin_level_data = .surv_overlap,
    # make_pop_run_file() reads only `outdir`/`population` off the opt, so the
    # occurrence fit's training opt is as good a key as any track's.
    # Population arguments in one place (simulator_core.R::sim_frame_pop_args():
    # location_frame.gid_pop_path / wpp_path / pop_year).
    gid_pop_path = .fpop$gid_pop_path,
    wpp_path = .fpop$wpp_path,
    pop_year = .fpop$pop_year,
    pop_window = .fpop$pop_window
  )
  # gids whose min_seas was defaulted, checked on the final frame.
  .ms_defaulted <- attr(frame, "min_seas_defaulted")
  # TTLO level "none": first recorded onset per frame gid,
  # attached below only when the occurrence fit carries the level.
  .first_onset <- if (!is.null(.surv_overlap)) {
    tryCatch(frame_first_onset(.surv_overlap, frame$gid), error = function(e) {
      NULL
    })
  } else {
    NULL
  }
  rm(.surv_overlap)
  message(sprintf(
    "[07_simulate_v1] location frame: %d locations, %d ADM1, %d countries, %d TU(s) (incl. NA)",
    nrow(frame),
    length(unique(frame$adm1_unit)),
    length(unique(frame$country)),
    length(unique(frame$tu))
  ))

  # Apply scope.countries filter (ISO3 subset; severs TU-lag → warn)
  ctry_scope <- sim$scope$countries %||% "all"
  if (!identical(ctry_scope, "all")) {
    iso3_req <- as.character(ctry_scope)
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
      "[07_simulate_v1] scope.countries='%s': %d → %d locations",
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

  # Zero-filled climate gate (mirrors utils/simulations.R::load_artefacts()).
  .assert_no_zero_filled_climate(
    .draws_fe_names(draws_occ, draws_dur, draws_mon),
    panel,
    window = .sim_panel_window(sim),
    label = "07_simulate_v1 load_artefacts"
  )

  if (.draws_has_ttlo_none(draws_occ)) {
    if (is.null(.first_onset)) {
      stop(
        "[07_simulate_v1] the occurrence fit carries the TTLO level 'none' but the survival ",
        "panel (first recorded onsets) could not be read.",
        call. = FALSE
      )
    }
    frame$first_onset <- .first_onset
  }
  frame <- check_frame_vs_draws(frame, draws_occ, sim$scope)
  # The season onset gate needs every simulated gid's own season start.
  check_onset_eligibility_frame(
    frame,
    onset_eligibility_of(opt$simulator),
    "07b_simulator_v1",
    defaulted = .ms_defaulted
  )

  n_replicates <- as.integer(sim$n_replicates %||% 100L)

  # ---------------------------------------------------------------------------
  # Per-replicate ENSO forecast scenarios
  # Three branches (mutually exclusive; counterfactual key takes precedence):
  #   enso.counterfactual.scenario = "analog"  → sample_enso_analog_scenarios()
  #   enso.counterfactual.scenario = "lenssen" → sample_enso_lenssen_scenarios()
  #   (default / factual, enso.c3s_plume set)  → sample_enso_forecast_trajectories()
  # ---------------------------------------------------------------------------
  enso_scenarios <- NULL
  ctf_scenario <- sim$enso$counterfactual$scenario # NULL for factual runs

  if (!is.null(ctf_scenario) || !is.null(sim$enso$c3s_plume)) {
    tryCatch(
      {
        enso_source <- sim$enso$source %||% "oni"
        historical_oni <- fetch_enso_index(
          source = enso_source,
          cache_dir = "generated_data/",
          refresh = FALSE
        )
        # Derive simulation months from the panel (t_start >= horizon start).
        horizon_start <- as.Date(sim$horizon$start)
        sim_months <- sort(unique(panel$t_start[
          panel$t_start >= horizon_start
        ]))

        if (!is.null(ctf_scenario)) {
          # ── Counterfactual branches ─────────────────────────────────────────
          message(
            "[07_simulate_v1] ENSO counterfactual scenario: ",
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
              n_rep = n_replicates
            )
            message(sprintf(
              "[07_simulate_v1] analog ENSO: %d lags, %d months, %d replicates (deterministic)",
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
            # lenssen_horizon_start: date at which piControl analogs take over from
            # observed ONI; NULL / "" = sim$horizon$start (default).
            # One resolution of the block, analog_set included
            # (utils/enso.R::lenssen_counterfactual_args()).
            lcf <- lenssen_counterfactual_args(
              sim$enso$counterfactual,
              horizon_start
            )
            enso_horizon <- lcf$horizon_start
            if (enso_horizon != horizon_start) {
              message(sprintf(
                "[07_simulate_v1] Lenssen ENSO horizon: %s (overrides sim horizon %s)",
                enso_horizon,
                horizon_start
              ))
            }
            enso_scenarios <- sample_enso_lenssen_scenarios(
              cache_path = lcf$cache_path,
              sim_months = sim_months,
              n_rep = n_replicates,
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
              "[07_simulate_v1] Lenssen ENSO: %d lags, %d months × %d replicates",
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
          } else {
            stop(
              "[07_simulate_v1] Unknown enso.counterfactual.scenario: '",
              ctf_scenario,
              "'. Must be 'analog' or 'lenssen'."
            )
          }
        } else {
          # ── Factual branch (C3S plume) ────────────────────────────────────────
          message(
            "[07_simulate_v1] building per-replicate ENSO forecast trajectories"
          )
          plume_data <- read_enso_c3s_plume(sim$enso$c3s_plume)
          message(sprintf(
            "[07_simulate_v1] C3S plume: %d rows, dates %s – %s, quantile levels %s",
            nrow(plume_data),
            min(plume_data$date),
            max(plume_data$date),
            paste(sort(unique(plume_data$q_level)), collapse = ", ")
          ))

          # Raw Niño 3.4 for the rolling-mean conversion (preferred source).
          historical_nino34 <- tryCatch(
            fetch_enso_index(
              source = "nino34",
              cache_dir = "generated_data/",
              refresh = FALSE
            ),
            error = function(e) {
              message(
                "[07_simulate_v1] nino34 cache not found — using ONI as Niño 3.4 proxy."
              )
              NULL
            }
          )

          enso_scenarios <- sample_enso_forecast_trajectories(
            plume_data = plume_data,
            n_rep = n_replicates,
            historical_oni = historical_oni,
            historical_nino34 = historical_nino34,
            sim_months = sim_months,
            max_lag = 12L,
            max_lead = panel_max_enso_lead(panel),
            seed = as.integer(sim$seed %||% 20260603L)
          )
          message(sprintf(
            "[07_simulate_v1] ENSO scenarios built: %d lags, %d months × %d replicates",
            length(enso_scenarios),
            if (length(enso_scenarios) > 0L) nrow(enso_scenarios[[1]]) else 0L,
            if (length(enso_scenarios) > 0L) ncol(enso_scenarios[[1]]) else 0L
          ))
        }
      },
      error = function(e) {
        warning(
          "[07_simulate_v1] ENSO scenario build failed: ",
          conditionMessage(e),
          " — falling back to shared panel ENSO values."
        )
        enso_scenarios <<- NULL
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Per-replicate warmup panels (optional, from warmup config key)
  # ---------------------------------------------------------------------------
  warmup_panels_list <- NULL
  if (isTRUE(sim$warmup$use_sim_trajectories)) {
    message("[07_simulate_v1] building per-replicate hybrid warmup panels")
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

        # Load the survival panel for observed lag columns (used for observed gids
        # in the warmup window). Prefer opt_occ$baseline (from fit attribute); fall
        # back to opt$baseline (top-level config key) when the fit carries no opt.
        surv_path <- tryCatch(
          make_surv_dat_file(list(
            baseline = opt_occ$baseline %||% opt$baseline
          )),
          error = function(e) NULL
        )
        if (is.null(surv_path) || !file.exists(surv_path %||% "")) {
          # Fall back to panel itself (has lag columns from 03c/03d).
          surv_panel_warmup <- panel
          message(
            "[07_simulate_v1] warmup: using projection panel as surrogate surv_panel ",
            "(survival panel path not resolved from fit opt)."
          )
        } else {
          message(
            "[07_simulate_v1] warmup: reading survival panel: ",
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
          "[07_simulate_v1] warmup panels built: %d replicates, each covering %s – %s",
          length(warmup_panels_list),
          warmup_start,
          warmup_end
        ))
      },
      error = function(e) {
        warning(
          "[07_simulate_v1] warmup panel build failed: ",
          conditionMessage(e),
          " — falling back to shared panel-based warmup for all replicates."
        )
        warmup_panels_list <<- NULL
      }
    )
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
    enso_scenarios = enso_scenarios,
    warmup_panels_list = warmup_panels_list,
    cb_meta_occ = cb_meta_occ,
    cb_meta_dur = cb_meta_dur,
    cb_meta_mon = cb_meta_mon
  )
}

# 6. Frame vs draws checks -----------------------------------------------------
# check_frame_vs_draws — admin-level alignment guards.

#' Validate and (if needed) filter the location frame against the occurrence
#' draws bundle. Three checks:
#'
#' 1. **ADM1-level gids under ADM2 scope.** When `scope.admin_level = "ADM2"`
#'    but some gids satisfy `gid == adm1_unit`, those locations come from
#'    ADM1-level fit data and will be simulated at ADM1 resolution. A warning
#'    is issued so the caller knows the declared admin level is only partially
#'    achieved.
#'
#' 2. **Gids absent from the adm2 RE matrix (unseen_adm1 policy).** If a gid
#'    is in the frame but has no column in `draws$re_intercept$adm2`, the
#'    kernel silently falls back to the parent-level REs (adm1 + country +
#'    region) with adm2 RE = 0. Behaviour is controlled by `scope.unseen_adm1`:
#'    - `"restrict"` (default): drop those gids from the frame and emit a
#'      message. Errors if all locations are removed.
#'    - `"marginalise"`: warn and proceed; the kernel's silent-zero fallback is
#'      already the correct marginalisation.
#'
#' 3. **ADM1 units absent from the adm1 RE matrix.** Same silent-zero fallback
#'    applies at the adm1 bucket. Always a warning (no filtering policy).
#'
#' @param frame tibble from build_location_frame().
#' @param draws_occ draws bundle from extract_occurrence_draws().
#' @param sim_scope list from `opt$simulator$scope`.
#' @return frame, possibly with rows removed (restrict policy).
check_frame_vs_draws <- function(frame, draws_occ, sim_scope) {
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
    unseen_mask <- !(re_level_key(frame, "adm2") %in% adm2_cols)
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
    unseen_adm1_mask <- !(re_level_key(frame, "adm1") %in% adm1_cols)
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

  frame
}

# Kernel wrapper around simulate_one_tu(): slices artefacts down to one TU.
run_one_task <- function(
  tu,
  replicate_ids,
  posterior_idx_occ,
  seed,
  output_part,
  art,
  opt,
  dry_run = FALSE,
  posterior_idx_dur = NULL,
  posterior_idx_mon = NULL
) {
  # "_NA_TU_" is the synthetic bucket for locations whose country is absent
  # from the HMM cluster file (tu = NA at fit time).
  frame_tu <- if (identical(tu, "_NA_TU_")) {
    art$frame %>% dplyr::filter(is.na(.data$tu))
  } else {
    # `.env$tu` is required: under data-masking a bare `tu` on the RHS binds to
    # the column, making the predicate `tu == tu` (always TRUE) and simulating
    # every TU per task. `.env$tu` forces the function argument.
    art$frame %>% dplyr::filter(.data$tu == .env$tu)
  }
  if (nrow(frame_tu) == 0) {
    warning(sprintf(
      "[run_one_task] TU '%s' has no locations in the frame — skipping.",
      tu
    ))
    return(invisible(NULL))
  }
  panel_tu <- build_climate_panel_for_tu(
    art$panel,
    frame_tu,
    horizon = opt$simulator$horizon,
    scenario = opt$simulator$scenario$climate %||% "factual",
    end_mode = "closed" # forward horizon [start, end] (utils/base.R)
  )

  # Grouped mode: build the per-group routing list from art's by_group fields.
  # draws_dur_by_group is NULL unless grouping.duration_enable was set at load
  # time (see .load_artefacts_grouped()); simulate_one_tu() only consults it
  # when has_dur (draws_dur != NULL) is TRUE. draws_mon_by_group /
  # cb_meta_mon_by_group follow the same opt-in rule via
  # grouping.monthly_enable and has_mon.
  grp_arg <- if (isTRUE(art$grouped)) {
    list(
      draws_occ_by_group = art$draws_occ_by_group,
      draws_dur_by_group = art$draws_dur_by_group,
      draws_mon_by_group = art$draws_mon_by_group,
      cb_meta_occ_by_group = art$cb_meta_occ_by_group,
      cb_meta_dur_by_group = art$cb_meta_dur_by_group,
      cb_meta_mon_by_group = art$cb_meta_mon_by_group
    )
  } else {
    NULL
  }

  # Neighbour-contagion build-time guard: on only when the config flag is set
  # AND the (per-group) occurrence fit carries an nbr_lag1_r* FE. FALSE ⇒ every
  # nbr branch in simulate_one_tu() is inert, byte-identical to no contagion.
  nbr_cfg <- nbr_cfg_from_sim(opt$simulator)
  nbr_spec <- resolve_nbr_spec(nbr_cfg, art$draws_occ, grp_arg)
  use_nbr <- nbr_spec$use

  simulate_one_tu(
    tu = tu,
    frame_tu = frame_tu,
    panel_tu = panel_tu,
    draws_occ = art$draws_occ,
    draws_dur = art$draws_dur,
    draws_mon = art$draws_mon,
    replicate_ids = replicate_ids,
    posterior_idx_occ = posterior_idx_occ,
    posterior_idx_dur = posterior_idx_dur,
    posterior_idx_mon = posterior_idx_mon,
    seed = seed,
    output_part = output_part,
    horizon = opt$simulator$horizon,
    out_fmt = opt$simulator$output$format %||% "qs",
    dry_run = dry_run,
    enso_scenarios = art$enso_scenarios,
    warmup_panels_list = art$warmup_panels_list,
    cb_meta_occ = art$cb_meta_occ,
    cb_meta_dur = art$cb_meta_dur,
    cb_meta_mon = art$cb_meta_mon,
    grouping = grp_arg,
    drop_country_tu_spillover = opt$simulator$drop_country_tu_spillover %||%
      FALSE,
    use_nbr = use_nbr,
    nbr_rho_km = nbr_cfg$rho_km %||% 500L,
    nbr_gadm_path = nbr_cfg$gadm_path %||% "data/geodata/gadm_valid.gpkg",
    nbr_cache_dir = nbr_cfg$cache_dir %||% "generated_data/",
    nbr_spec = nbr_spec,
    clamp_enso_to_knots = isTRUE(opt$simulator$enso$clamp_to_knots),
    clamp_enso_margin = as.numeric(opt$simulator$enso$clamp_margin %||% 0.25),
    # Monthly cases per location <= incidence_cap x pop (simulator.incidence_cap).
    incidence_cap = incidence_cap_of(opt$simulator),
    # Onset gate (simulator.onset_eligibility).
    onset_eligibility = onset_eligibility_of(opt$simulator),
    # Cessation-month exposure (simulator.cessation_exposure).
    cessation_exposure = cessation_exposure_of(opt$simulator)
  )
}

#' Dry-run placeholder — prints the planned assignment without invoking the
#' kernel. Used by --dry-run so the validator can be exercised cheaply without
#' paying the fit-load + draws-extraction cost.
simulate_one_tu_stub <- function(
  tu,
  replicate_ids,
  posterior_idx,
  seed,
  output_part,
  dry_run = TRUE
) {
  message(sprintf(
    "[DRY] tu=%s  replicates=[%s]  posterior_idx=[%s]  seed=%s  out=%s",
    tu,
    paste(range(replicate_ids), collapse = ".."),
    paste(range(posterior_idx), collapse = ".."),
    seed,
    output_part
  ))
  invisible(NULL)
}

# 7. Task plan -----------------------------------------------------------------
# TU x replicate-chunk assignments.
# build_task_plan() and write_task_yamls() live in analysis/utils.R (sourced
# above) as the single source of truth shared with 07a_simulator_tasks.R.

# 8. Dispatch ------------------------------------------------------------------

art <- if (isTRUE(opt$dry_run)) NULL else load_artefacts(opt)

# Fallback: when *every* location still has NA TU after the location-frame
# rehydration step (HMM cluster file truly absent, not just missing for some
# countries), collapse them into a single synthetic TU `_ALL_` so the kernel
# can still run with the TU spillover lag effectively disabled. Step-3
# rehydration normally fills `tu` from `data/hmm_country_clusters.csv`.
if (!is.null(art) && all(is.na(art$frame$tu))) {
  message(
    "[07_simulate_v1] every location has NA TU after rehydration — collapsing ",
    "to a single synthetic TU `_ALL_` (TU spillover lag will be inactive)."
  )
  art$frame$tu <- "_ALL_"
}

if (!is.null(opt$task)) {
  # ----- single-task mode (SLURM array slot) -----
  if (!file.exists(opt$task)) {
    stop("--task file not found: ", opt$task)
  }
  task <- yaml::read_yaml(opt$task)
  required_task_keys <- c(
    "tu",
    "replicate_ids",
    "posterior_idx",
    "seed",
    "output_part"
  )
  missing_keys <- setdiff(required_task_keys, names(task))
  if (length(missing_keys) > 0) {
    stop(
      "Task YAML missing required keys: ",
      paste(missing_keys, collapse = ", ")
    )
  }

  message("[07_simulate_v1] single-task mode: ", opt$task)
  if (isTRUE(opt$dry_run)) {
    simulate_one_tu_stub(
      tu = task$tu,
      replicate_ids = task$replicate_ids,
      posterior_idx = task$posterior_idx,
      seed = task$seed,
      output_part = task$output_part
    )
  } else {
    run_one_task(
      tu = task$tu,
      replicate_ids = task$replicate_ids,
      posterior_idx_occ = task$posterior_idx_occ %||% task$posterior_idx,
      posterior_idx_dur = task$posterior_idx_dur %||% task$posterior_idx,
      posterior_idx_mon = task$posterior_idx_mon %||% task$posterior_idx,
      seed = task$seed,
      output_part = task$output_part,
      art = art,
      opt = opt
    )
  }
} else {
  # ----- full-run mode -----
  plan <- build_task_plan(opt, frame = if (!is.null(art)) art$frame else NULL)
  message(sprintf(
    "[07_simulate_v1] full-run plan: %d TU x %d chunk = %d tasks",
    length(unique(vapply(plan, `[[`, character(1), "tu"))),
    length(unique(vapply(plan, `[[`, character(1), "chunk_id"))),
    length(plan)
  ))

  if (!isTRUE(opt$dry_run)) {
    task_paths <- write_task_yamls(plan, tasks_dir)
    message(
      "[07_simulate_v1] wrote ",
      length(task_paths),
      " task YAML(s) under ",
      tasks_dir
    )
  }

  backend <- opt$simulator$parallel$backend %||% "furrr"

  if (identical(backend, "slurm")) {
    # SLURM path: task YAMLs are now on disk; execution is delegated to the
    # SLURM array.  Emit instructions and exit.
    if (!isTRUE(opt$dry_run)) {
      message(sprintf(
        "[07_simulate_v1] parallel.backend=slurm — %d task YAML(s) ready in %s",
        length(plan),
        tasks_dir
      ))
      message("[07_simulate_v1] Submit via:")
      message(sprintf(
        "[07_simulate_v1]   bash analysis/bash/stages/07_simulator_v1.sh --config %s",
        .flag_val("-C") %||% .flag_val("--config") %||% "<config.yml>"
      ))
    } else {
      message(sprintf(
        "[07_simulate_v1] [DRY] slurm: %d task(s) would be written to %s",
        length(plan),
        tasks_dir
      ))
    }
  } else {
    # furrr path: run tasks in parallel locally (workers=1 = sequential).
    if (isTRUE(opt$dry_run)) {
      for (task in plan) {
        simulate_one_tu_stub(
          tu = task$tu,
          replicate_ids = task$replicate_ids,
          posterior_idx = task$posterior_idx,
          seed = task$seed,
          output_part = task$output_part
        )
      }
    } else {
      suppressPackageStartupMessages({
        library(furrr)
        library(future)
      })
      workers <- as.integer(opt$simulator$parallel$workers %||% 1L)
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential), add = TRUE)
      message(sprintf(
        "[07_simulate_v1] furrr: %d worker(s), %d task(s).",
        workers,
        length(plan)
      ))

      # RNG: each task uses set.seed(task$seed + rep_id) inside simulate_one_tu(),
      # so results are independent of which worker processes the task (S7). Task
      # base seeds are strided by n_replicates (see build_task_plan) so the
      # (seed + rep_id) ranges never collide across (task, replicate) pairs.
      furrr::future_map(
        plan,
        function(task) {
          run_one_task(
            tu = task$tu,
            replicate_ids = task$replicate_ids,
            posterior_idx_occ = task$posterior_idx_occ %||% task$posterior_idx,
            posterior_idx_dur = task$posterior_idx_dur %||% task$posterior_idx,
            posterior_idx_mon = task$posterior_idx_mon %||% task$posterior_idx,
            seed = task$seed,
            output_part = task$output_part,
            art = art,
            opt = opt
          )
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    }
  }
}

message("[07_simulate_v1] done.")
