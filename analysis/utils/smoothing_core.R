# smoothing_core.R — shared kernel of the chunked spatPomp historical
#   smoothing (phase 08)
#
# Provides:
#   smoothing_make_config(), smoothing_load_artefacts()   08a: config, inputs
#   smoothing_filter_chunk()      08b Phase A: particle filter of one chunk
#   extract_smoothing_one_tu()    08b Phase B: smoothed trajectories of a TU
#   agg_and_ribbon(), agg_and_ribbon_mon(), agg_and_ribbon_active(),
#   ribbon_plot(), panel_pair()   08c: compile
#
# Details:
#   A `_core` kernel: it sources its own dependency closure (one line per
#   module, with why it is needed), so the 08 scripts source this file.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R") # drop_overlapping_adm1_from_frame() -- required by smoothing_load_artefacts()'s ADM1/ADM2 overlap filter
source("analysis/utils/paths.R") # make_pop_run_file()/make_agg_ts_file() -- required by smoothing_default_config()
source("analysis/utils/provenance.R") # sim_episode_check() -- required by smoothing_load_artefacts()
source("analysis/utils/dlnm.R") # build_cb_design()/build_mso_basis() -- required by extract_smoothing_one_tu() for brms_cb mon fits
source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population()/load_wpp_table() -- required by build_location_frame()'s gid_pop_path/wpp_path branches
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/sim_io.R") # build_obs_panel() & co.
source("analysis/utils/coverage.R") # build_coverage_mask() -- required by 08c_smoothing_compile.R
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatpomp_core_c.R") # Csnippet variant; sources spatpomp_core.R internally
source("analysis/utils/spatial_contagion.R") # build_partition_nbr_kernel()/resolve_use_nbr() -- required by build_spatpomp_tu_c() and smoothing_filter_one() whenever simulator.nbr_contagion.enable is set
source("analysis/utils/anchor.R") # resolve_base_rate_paths()/anchor_spec() -- transitive dependency of simulations.R's grouped build_location_frame() call (utils/README.md)
source("analysis/utils/simulations.R")
source("analysis/utils/who_dashboard.R") # who_obs_table()/who_coverage_fraction()/who_block_share() -- the `_who` measurement (simulator.filter.who_obs)
source("analysis/utils/footprint.R") # finest_level_rows() -- 08c's observed monthly cases (smoothing_obs_cases_mon())

# ---------------------------------------------------------------------------
# Config: defaults + YAML (-C) merge + CLI overrides
# ---------------------------------------------------------------------------
# Priority: CLI flag > YAML config > default.
# `args` is a named list that may contain: config, n_reps, Np, workers, tus
# (NULLs are ignored). Returns a config list `S` consumed by the loaders below.

smoothing_default_config <- function() {
  list(
    OCC_FIT_PATH = "generated_data/outbreaks_survival_poisson_fit_base-1990_2019_bh-seasonal_ctry-all_adm-adm1-adm2_uncorr_init-0.2_enso-onilag5-slope_noloc_adm2i_bhslp.rds",
    PANEL_PATH = "generated_data/sim_climate_panel_base-1990_2019.rds",
    SURV_PATH = "generated_data/outbreaks_survival_data_base-1990_2019.rds",
    # 03_prepare_siz.R's episode table. Not a fit input: build_location_frame()
    # reads it for the population / admin_level fill, and 08c_smoothing_compile.R
    # for the observed ongoing-outbreak series.
    SIZE_DATA_PATH = "generated_data/outbreak_size_analysis_data.rds",
    # True per-(episode, month) case panel (03_prepare_mon.R's `cases_month`),
    # i.e. what the mon track is actually fit to reproduce. Used by
    # 08c_smoothing_compile.R to build the mon-track's observed-cases series
    # (see obs_cases_mon_* there).
    MON_DATA_PATH = "generated_data/outbreaks_monthly_data_base-1990_2019.rds",
    # Raw `simulator.location_frame` config block (NULL when no -C config was
    # given). Kept verbatim -- not pre-resolved into per-track paths -- because
    # resolve_base_rate_paths() (utils/anchor.R) needs the per-track fit `opt`s,
    # which for the ungrouped path aren't known until after the fits are
    # loaded, and for the grouped path aren't known until a representative
    # per-group fit is discovered (see smoothing_load_artefacts()). NULL/empty
    # keeps every existing smoothing run byte-for-byte: with no anchored track,
    # resolve_base_rate_paths() returns NULL and no anchor column is resolved
    # onto the frame.
    LF_CFG = NULL,
    # Canonical gid population Q_g (02c_preprocess_population_compile.R). NULL-if-absent
    # so an existing tree keeps the legacy route; when present, the simulator
    # forecasts from a real population instead of the retired `pop <- 1.0`.
    GID_POP_PATH = make_pop_run_file(
      list(outdir = "generated_data"),
      "pop_gid"
    ),
    WPP_PATH = make_pop_run_file(list(outdir = "generated_data"), "pop_wpp"),
    # simulator.location_frame.pop_year (simulator_core.R::sim_pop_year_spec();
    # threaded from the config by smoothing_make_config()).
    POP_YEAR = .SIM_POP_YEAR_DEFAULT,
    HMM_PATH = "data/hmm_country_clusters.csv",
    AGG_TS_PATH = make_agg_ts_file(),
    AGG_ADM_TS_PATH = make_agg_adm_ts_file(),
    FILTER_START = as.Date("2010-01-01"),
    FILTER_END = as.Date("2023-12-01"),
    BURNIN_MONTHS = 24L,
    BURNIN_START = as.Date("2008-01-01"),
    N_REPS = 100L,
    NP = 25000L,
    WORKERS = 4L,
    SEED_ROOT = 20260609L,
    MAX_DRAWS_OCC = 200L,
    MAX_DRAWS_MON = 200L,
    MAX_DRAWS_DUR = 200L,
    MEASUREMENT = "occurrence_bernoulli",
    DELTA_FN = .SP_STATE_DELTA_FN,
    DELTA_FP = .SP_STATE_DELTA_FP,
    # simulator.filter.case_fp_rate / case_fp_mean (`_actn` measurements).
    CASE_FP_RATE = .SP_CASE_FP_RATE,
    CASE_FP_MEAN = .SP_MON_MU_FLOOR,
    COVARIATE_FEED = "latent",
    MONTHLY_CASE_OBS = "none",
    # simulator.filter.traj_storage / traj_prune_every (resolve_filter_cfg());
    # "pruned" is the default (identical() to "full").
    TRAJ_STORAGE = "pruned",
    TRAJ_PRUNE_EVERY = 1L,
    WHO_OBS = list(enable = FALSE),
    BLOCK_GRANULARITY = "country",
    PARTITION_BY = "transmission_unit",
    PARTITION_SOURCE = "",
    OUT_DIR = "figures/scratchpad/s10_historical_smoothing",
    outdir = "generated_data",
    run_tag = "s10_historical_smoothing",
    GROUPED = FALSE,
    GRP_SOURCE = NULL,
    GRP_FITS_DIR = NULL,
    # Duration track (optional, mirrors simulator.duration_fit on the forecast
    # side -- utils/simulations.R::load_artefacts()). NULL keeps the
    # pre-duration smoothing behaviour byte-for-byte: draws_dur stays NULL and
    # extract_smoothing_one_tu() takes its legacy onset-only path.
    DUR_FIT_PATH = NULL,
    GRP_DUR_ENABLE = FALSE,
    # Monthly case-count track (optional, mirrors simulator.monthly_fit on
    # the forecast side). Requires the duration track -- `active`/`mso` are
    # duration-model latent states. NULL keeps pre-mon behaviour byte-for-byte.
    MON_FIT_PATH = NULL,
    MON_CB_META_PATH = NULL,
    GRP_MON_ENABLE = FALSE,
    # Neighbour-contagion (distance-weighted) -- mirrors simulator.nbr_contagion
    # on the forecast side. FALSE keeps pre-change smoothing byte-for-byte: no
    # kernel built, no WNBR matvec emitted in the filter Csnippet. The on/off
    # decision additionally depends on the occurrence fit carrying an nbr_* FE
    # (resolved per-TU in smoothing_filter_one via resolve_use_nbr()).
    NBR_ENABLE = FALSE,
    NBR_RHO_KM = 500L,
    NBR_GADM_PATH = "data/geodata/gadm_valid.gpkg",
    NBR_CACHE_DIR = "generated_data/",
    # Monthly incidence cap (simulator.incidence_cap): monthly cases per unit
    # <= cap x pop, in the C filter's mon mean and in the post-hoc case draw.
    INCIDENCE_CAP = incidence_cap_of(NULL),
    # Onset eligibility (simulator.onset_eligibility): the C
    # filter's onset gate; "inactive" = the default.
    ONSET_ELIGIBILITY = onset_eligibility_of(NULL),
    # Cessation-month exposure (simulator.cessation_exposure):
    # the C filter's mon mean and the post-hoc case draw; "full" = the default.
    CESSATION_EXPOSURE = cessation_exposure_of(NULL),
    # Draw design (simulator.filter.draw_design):
    # "per_tu_random" = per-unit random draws (the default).
    DRAW_DESIGN = smoothing_draw_design_of(NULL)
  )
}

#' smoothing_who_scale_paths
#'
#' Where 08a writes the estimated WHO scale of a smoothing run:
#' <outdir>/sim_spatpomp/<run_tag>/who_scale.csv and who_scale_meta.yml,
#' next to tasks/ and parts/ (09 finds them as dirname(smoothing_parts_dir)).
smoothing_who_scale_paths <- function(S) {
  run_dir <- dirname(make_smooth_run_file(
    list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
    kind = "smooth_parts_dir"
  ))
  list(
    dir = run_dir,
    csv = file.path(run_dir, "who_scale.csv"),
    meta = file.path(run_dir, "who_scale_meta.yml")
  )
}

#' smoothing_write_who_scale
#'
#' Write the estimated WHO scale table (art$who_scale) and its meta (digest,
#' estimator settings, db_overlap, WHO rows skipped, and -- recorded, not
#' calibrated -- the detection rates delta_fn / delta_fp) for a smoothing run.
#' Called by 08a (and by 08b's local-loop mode). NULL (nothing written) in
#' fixed mode.
#' @return the table digest (digest_who_scale()), or NULL.
smoothing_write_who_scale <- function(S, art) {
  if (is.null(art$who_scale)) {
    return(invisible(NULL))
  }
  p <- smoothing_who_scale_paths(S)
  dir.create(p$dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(as.data.frame(art$who_scale), p$csv, row.names = FALSE)
  dg <- digest_who_scale(art$who_scale)
  se <- S$WHO_OBS$scale_est
  yaml::write_yaml(
    list(
      who_scale_digest = dg,
      scale_mode = S$WHO_OBS$scale_mode,
      scale_est = list(
        start = format(se$start),
        end = format(se$end),
        min_db_cases = se$min_db_cases,
        min_months = se$min_months,
        shrink_cases = se$shrink_cases,
        fallback = se$fallback,
        clamp = se$clamp
      ),
      fixed_scale = S$WHO_OBS$scale,
      db_overlap = S$WHO_OBS$db_overlap,
      who_window = list(
        start = format(S$WHO_OBS$start),
        end = format(S$WHO_OBS$end)
      ),
      who_rows_skipped = if (is.null(art$who_skipped)) 0L else art$who_skipped,
      delta_fn = S$DELTA_FN,
      delta_fp = S$DELTA_FP,
      run_tag = S$run_tag,
      written_at = format(Sys.time())
    ),
    p$meta
  )
  message("[smoothing] who_obs scale table -> ", p$csv, " (digest ", dg, ")")
  invisible(dg)
}

#' smoothing_frame_pop_args
#'
#' The population arguments of the smoothing's build_location_frame() calls
#' (simulator_core.R::sim_frame_pop_args() on the smoothing state): the
#' configured location_frame paths or the S defaults (NULL when absent), the
#' threaded S$POP_YEAR, and the months the smoothing reads (burn-in start to
#' filter end) as the pop_year-mode window.
smoothing_frame_pop_args <- function(S) {
  pa <- sim_frame_pop_args(
    list(location_frame = S[["LF_CFG"]]),
    gid_pop_path = S$GID_POP_PATH,
    wpp_path = S$WPP_PATH
  )
  pa$pop_year <- S[["POP_YEAR"]] %||% .SIM_POP_YEAR_DEFAULT
  pa$pop_window <- c(as.Date(S$BURNIN_START), as.Date(S$FILTER_END))
  pa
}

smoothing_make_config <- function(args = list()) {
  S <- smoothing_default_config()
  cfg_n_reps <- NULL
  cfg_Np <- NULL
  cfg_workers <- NULL

  if (!is.null(args$config)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("[smoothing] yaml package required for -C config support.")
    }
    if (!file.exists(args$config)) {
      stop("[smoothing] config not found: ", args$config)
    }

    # Resolve through make_options_from_config() -- the same entry point
    # 07b_simulator_v1.R / 07b_simulator_spatpomp.R use -- rather than a bare
    # yaml::read_yaml(), so the get_default_options()$simulator defaults
    # (utils/config.R; e.g. location_frame.*, grouping.*,
    # filter.delta_fn/delta_fp) reach 08 smoothing exactly as they reach the
    # 07 forecast for the same config file. A raw read would silently miss
    # every key the config does not spell out.
    cfg <- yaml::read_yaml(args$config)
    opt <- make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", args$config),
      model = NULL
    )
    sim <- opt$simulator %||% list()
    flt <- sim$filter %||% list()

    # n_replicates / filter.Np / parallel.workers are the one deliberate
    # exception: smoothing and forecast runs share the same `simulator:`
    # config schema but intentionally use DIFFERENT numeric defaults when a
    # config omits them (e.g. Np 25000 here vs. the forecast side's 500;
    # workers 4 here vs. 1) -- see smoothing_default_config() above. Reading
    # these three from the merged `opt` would silently overwrite that
    # deliberate divergence with get_default_options()'s forecast-oriented
    # values (parallel.workers defaults to 1L there), so they stay sourced
    # from the raw, unmerged config exactly as before.
    raw_sim <- cfg$simulator %||% list()
    raw_flt <- raw_sim$filter %||% list()
    raw_par <- raw_sim$parallel %||% list()

    if (!is.null(sim$occurrence_fit)) {
      S$OCC_FIT_PATH <- sim$occurrence_fit
    }
    if (!is.null(sim$climate_panel)) {
      S$PANEL_PATH <- sim$climate_panel
    }
    if (!is.null(sim$duration_fit) && nzchar(sim$duration_fit)) {
      S$DUR_FIT_PATH <- sim$duration_fit
    }
    if (!is.null(sim$monthly_fit) && nzchar(sim$monthly_fit)) {
      S$MON_FIT_PATH <- sim$monthly_fit
    }
    if (!is.null(sim$monthly_cb_meta) && nzchar(sim$monthly_cb_meta)) {
      S$MON_CB_META_PATH <- sim$monthly_cb_meta
    }

    lf <- sim$location_frame %||% list()
    if (!is.null(lf$hmm_clusters_path)) {
      S$HMM_PATH <- lf$hmm_clusters_path
    }
    if (!is.null(lf$size_data_path)) {
      S$SIZE_DATA_PATH <- lf$size_data_path
    }
    # Population keys: a template's location_frame.pop_year / gid_pop_path /
    # wpp_path must reach 08 as well as 07. One reader for the year; the paths as the
    # forecast loaders
    # take them (simulator_core.R::sim_frame_pop_args()).
    S$POP_YEAR <- sim_pop_year_spec(sim)
    if (!is.null(lf[["gid_pop_path"]])) {
      S$GID_POP_PATH <- lf[["gid_pop_path"]]
    }
    if (!is.null(lf[["wpp_path"]])) {
      S$WPP_PATH <- lf[["wpp_path"]]
    }
    # Stashed raw; resolve_base_rate_paths() reads lf$base_rate_table_paths
    # (per-track overrides) itself at the point of use in
    # smoothing_load_artefacts().
    S$LF_CFG <- lf

    baseline <- cfg$baseline %||% "1990_2019"
    # Episode key (top-level `episodes:` block): selects the "_ep-gs" prep
    # artefacts under gid_season and is checked against every fit and panel in
    # smoothing_load_artefacts(). The default resolves the historical paths.
    S$EPISODE_KEY <- episode_key_of(opt)
    # National zero-fill occurrence panel/fits ("_nz", episodes.national_zero_fill):
    # the occurrence fit lookups below must ask for the matching fit.
    S$NATIONAL_ZERO_FILL <- national_zero_fill_of(opt)
    S$NATIONAL_ZERO_THRESHOLD <- national_zero_threshold_of(opt)
    # The whole occurrence-panel spec (utils/paths.R::occ_panel_spec_of()): the occurrence
    # fit lookups match its variant tokens exactly.
    S$OCC_PANEL <- occ_panel_spec_of(opt)
    ep_opt <- list(baseline = baseline, episodes = opt$episodes)
    surv_candidate <- as.character(make_surv_dat_file(ep_opt))
    if (file.exists(surv_candidate)) {
      S$SURV_PATH <- surv_candidate
    }
    mon_candidate <- as.character(make_monthly_dat_file(ep_opt))
    if (file.exists(mon_candidate)) {
      S$MON_DATA_PATH <- mon_candidate
    }
    if (is.null(lf$size_data_path) && identical(S$EPISODE_KEY, "gid_season")) {
      S$SIZE_DATA_PATH <- as.character(make_size_dat_file(ep_opt))
    }

    # get_default_options() fills filter.warmup_window with the forecast side's
    # fallback (2021-07-01 .. 2023-05-01, 07b_simulator_spatpomp.R), so this
    # `%||% list()` (= keep BURNIN_START / BURNIN_MONTHS above) is reached only
    # by a hand-built opt; every smoothing config sets warmup_window itself.
    # forecast_window defaults to NULL, so its `%||% list()` stays live.
    ww <- flt$warmup_window %||% list()
    fw <- flt$forecast_window %||% list()
    if (!is.null(fw$start)) {
      S$FILTER_START <- as.Date(fw$start)
    }
    if (!is.null(fw$end)) {
      S$FILTER_END <- as.Date(fw$end)
    }
    if (!is.null(ww$start)) {
      S$BURNIN_START <- as.Date(ww$start)
    }
    if (!is.null(ww$start) && !is.null(fw$start)) {
      S$BURNIN_MONTHS <- as.integer(round(
        as.numeric(difftime(
          as.Date(fw$start),
          as.Date(ww$start),
          units = "days"
        )) /
          30.44
      ))
    }

    if (!is.null(cfg$seed)) {
      S$SEED_ROOT <- as.integer(cfg$seed)
    }

    S$outdir <- cfg$outdir %||% "generated_data"
    S$run_tag <- sim$run_tag %||% "s10_historical_smoothing"
    S$OUT_DIR <- file.path("figures", "scratchpad", S$run_tag)

    cfg_n_reps <- as.integer(raw_sim$n_replicates %||% S$N_REPS)
    cfg_Np <- as.integer(raw_flt$Np %||% S$NP)
    cfg_workers <- as.integer(raw_par$workers %||% S$WORKERS)

    cfg_draws <- as.integer(sim$max_draws %||% 200L)
    S$MAX_DRAWS_OCC <- cfg_draws
    S$MAX_DRAWS_MON <- cfg_draws
    S$MAX_DRAWS_DUR <- cfg_draws
    # Config-level default: the monthly case-count measurement. A config that omits
    # `filter.measurement` now selects the mon path (and therefore hard-requires
    # simulator.monthly_fit + simulator.duration_fit) rather than silently
    # running the Bernoulli filter. `smoothing_default_config()` keeps
    # occurrence_bernoulli for the no-config path, which cannot supply either
    # artefact.
    # Measurement-side keys resolved by the shared resolve_filter_cfg()
    # (spatpomp_core.R) -- the same function the forecast driver uses.
    fc <- resolve_filter_cfg(flt)
    S$MEASUREMENT <- fc$measurement
    S$DELTA_FN <- fc$delta_fn
    S$DELTA_FP <- fc$delta_fp
    S$CASE_FP_RATE <- fc$case_fp_rate
    S$CASE_FP_MEAN <- fc$case_fp_mean
    S$COVARIATE_FEED <- fc$covariate_feed
    S$MONTHLY_CASE_OBS <- fc$monthly_case_obs
    S$MON_CASE_COL <- fc$mon_case_col
    S$TRAJ_STORAGE <- fc$traj_storage
    S$TRAJ_PRUNE_EVERY <- fc$traj_prune_every
    # WHO dashboard national monthly counts (`_who` measurement); resolved and
    # cross-checked against the measurement by resolve_filter_cfg().
    S$WHO_OBS <- fc$who_obs
    S$BLOCK_GRANULARITY <- flt$block_granularity %||% "country"
    S$PARTITION_BY <- flt$partition_by %||% "transmission_unit"
    S$PARTITION_SOURCE <- as.character(flt$partition_source %||% "")

    # Neighbour-contagion config (mirrors the forecast side's
    # simulator.nbr_contagion block). Pulled from the merged `sim`, so a
    # config_defaults.yml default reaches smoothing too. The on/off guard is
    # finalised per-TU in smoothing_filter_one() (needs the occurrence fit's FEs).
    nbr <- sim$nbr_contagion %||% list()
    S$NBR_ENABLE <- isTRUE(nbr$enable)
    S$NBR_RHO_KM <- nbr$rho_km %||% 500L
    S$NBR_GADM_PATH <- nbr$gadm_path %||% "data/geodata/gadm_valid.gpkg"
    S$NBR_CACHE_DIR <- nbr$cache_dir %||% "generated_data/"
    # Contagion registry for a registry-partitioned (`_part<token>`) term;
    # defaults to the partition registry (resolved again in nbr_cfg_from_sim()).
    S$NBR_PARTITION_SOURCE <- as.character(nbr$partition_source %||% "")
    S$INCIDENCE_CAP <- incidence_cap_of(sim)
    S$ONSET_ELIGIBILITY <- onset_eligibility_of(sim)
    S$CESSATION_EXPOSURE <- cessation_exposure_of(sim)
    S$DRAW_DESIGN <- smoothing_draw_design_of(flt)

    grp <- sim$grouping %||% list(enable = FALSE)
    if (isTRUE(grp$enable)) {
      S$GROUPED <- TRUE
      S$GRP_SOURCE <- grp$source %||%
        stop("[smoothing] simulator.grouping.source required")
      S$GRP_FITS_DIR <- grp$fits_dir %||%
        stop("[smoothing] simulator.grouping.fits_dir required")
      if (!dir.exists(S$GRP_FITS_DIR)) {
        stop(
          "[smoothing] simulator.grouping.fits_dir not found: ",
          S$GRP_FITS_DIR
        )
      }
      # Opt-in per-group duration fits, same key and same rationale as the
      # forecast loader (utils/simulations.R::load_artefacts()): batches with no
      # per-group duration fits deployed yet keep loading.
      S$GRP_DUR_ENABLE <- isTRUE(grp$duration_enable)
      S$GRP_MON_ENABLE <- isTRUE(grp$monthly_enable)
      message(sprintf(
        "[smoothing] grouped mode: source='%s' fits_dir='%s'",
        S$GRP_SOURCE,
        S$GRP_FITS_DIR
      ))
    }
    message("[smoothing] config loaded: ", args$config)
  }

  # CLI flags win over config; config wins over default.
  S$N_REPS <- args$n_reps %||% cfg_n_reps %||% S$N_REPS
  S$NP <- args$Np %||% cfg_Np %||% S$NP
  S$WORKERS <- args$workers %||% cfg_workers %||% S$WORKERS
  S$BLOCK_GRANULARITY <- args$block_granularity %||%
    S$BLOCK_GRANULARITY %||%
    "country"
  S$PARTITION_BY <- args$partition_by %||%
    S$PARTITION_BY %||%
    "transmission_unit"
  S$PARTITION_SOURCE <- as.character(
    args$partition_source %||% S$PARTITION_SOURCE %||% ""
  )

  S$hist_window <- list(
    start = format(S$BURNIN_START),
    end = format(S$FILTER_END)
  )
  S
}

# ---------------------------------------------------------------------------
# Artefact loading: fits/draws, frame, observed panels
# ---------------------------------------------------------------------------

.load_cb_sidecar <- function(fit_path, override = NULL) {
  cb_path <- override %||%
    sub("_fit_base-", "_cb_meta_base-", fit_path, fixed = TRUE)
  if (identical(cb_path, fit_path)) {
    return(NULL)
  }
  if (!file.exists(cb_path)) {
    return(NULL)
  }
  message("[smoothing] reading cross-basis meta: ", cb_path)
  readRDS(cb_path)
}

smoothing_load_artefacts <- function(S) {
  for (p in c(S$PANEL_PATH, S$SURV_PATH, S$SIZE_DATA_PATH)) {
    if (!file.exists(p)) stop("[smoothing] missing: ", p)
  }

  message("[smoothing] reading climate panel: ", S$PANEL_PATH)
  panel <- readRDS(S$PANEL_PATH)
  message("[smoothing] reading survival data: ", S$SURV_PATH)
  surv_dat <- readRDS(S$SURV_PATH)

  draws_occ_by_group <- NULL
  cb_meta_occ_by_group <- NULL
  country_to_group <- NULL
  draws_dur <- NULL
  draws_dur_by_group <- NULL
  draws_mon <- NULL
  draws_mon_by_group <- NULL
  cb_meta_mon <- NULL
  cb_meta_mon_by_group <- NULL
  # The duration fit has its own cross-basis sidecar for exactly the same reason
  # occ/mon do: `04_fit_dur_brms_cb.R` exists and the deployed per-group
  # best-lag duration fits are `_cb` fits carrying `enso_cb<k>` fixed effects at
  # a lag chosen independently of the other tracks (Horn: dur 8 vs occ 11).
  # Without it the duration LP is evaluated against the occurrence basis.
  cb_meta_dur <- NULL
  cb_meta_dur_by_group <- NULL
  # Every fit this run pairs, for the episode-definition check below.
  ep_fits <- list()

  if (!S$GROUPED) {
    if (!file.exists(S$OCC_FIT_PATH)) {
      stop("[smoothing] missing: ", S$OCC_FIT_PATH)
    }
    message("[smoothing] reading occurrence fit: ", S$OCC_FIT_PATH)
    occ_fit <- readRDS(S$OCC_FIT_PATH)
    ep_fits[["occ"]] <- occ_fit

    cb_meta_occ <- .load_cb_sidecar(S$OCC_FIT_PATH)

    opt_occ <- .fit_opt_from_attr(occ_fit)
    if (is.null(opt_occ$occurrence)) {
      opt_occ$occurrence <- list(backend = "brms")
    }

    draws_occ <- extract_occurrence_draws(
      occ_fit,
      opt_occ,
      max_draws = S$MAX_DRAWS_OCC
    )
    message(sprintf("[smoothing] draws_occ S=%d", draws_occ$S))

    # NULL unless the optional dur/mon blocks below run -- read by
    # resolve_base_rate_paths() ahead of build_location_frame() regardless of
    # whether either track is part of this run.
    opt_dur <- NULL
    opt_mon <- NULL

    # Duration fit (optional). Mirrors utils/simulations.R::load_artefacts(),
    # including the opt$duration backfill for fits whose attached opt predates
    # the key. NULL DUR_FIT_PATH => draws_dur stays NULL and the Phase-B
    # extractor keeps its legacy onset-only behaviour.
    if (!is.null(S$DUR_FIT_PATH)) {
      if (!file.exists(S$DUR_FIT_PATH)) {
        stop("[smoothing] missing duration fit: ", S$DUR_FIT_PATH)
      }
      message("[smoothing] reading duration fit: ", S$DUR_FIT_PATH)
      dur_fit <- readRDS(S$DUR_FIT_PATH)
      ep_fits[["dur"]] <- dur_fit
      opt_dur <- .fit_opt_from_attr(dur_fit)
      if (is.null(opt_dur$duration)) {
        opt_dur$duration <- list(backend = "brms")
      }
      draws_dur <- extract_duration_draws(
        dur_fit,
        opt_dur,
        max_draws = S$MAX_DRAWS_DUR
      )
      cb_meta_dur <- .load_cb_sidecar(S$DUR_FIT_PATH)
      message(sprintf("[smoothing] draws_dur S=%d", draws_dur$S))
    }

    # Monthly case-count fit (optional). Hard-requires the duration track: the mon mu
    # block is driven by `mso`, and `active`/`mso` are duration-model latent
    # states -- there is no invented fallback.
    if (!is.null(S$MON_FIT_PATH)) {
      if (is.null(draws_dur)) {
        stop(
          "[smoothing] simulator.monthly_fit requires simulator.duration_fit ",
          "-- the monthly case-count kernel reads `active`/`mso`, which are ",
          "duration-model latent states."
        )
      }
      if (!file.exists(S$MON_FIT_PATH)) {
        stop("[smoothing] missing monthly fit: ", S$MON_FIT_PATH)
      }
      message("[smoothing] reading monthly fit: ", S$MON_FIT_PATH)
      mon_fit <- readRDS(S$MON_FIT_PATH)
      ep_fits[["mon"]] <- mon_fit
      opt_mon <- .fit_opt_from_attr(mon_fit)
      if (is.null(opt_mon$monthly)) {
        opt_mon$monthly <- list(backend = "brms_cb")
      }
      draws_mon <- extract_monthly_draws(
        mon_fit,
        opt_mon,
        max_draws = S$MAX_DRAWS_MON
      )
      cb_meta_mon <- .load_cb_sidecar(S$MON_FIT_PATH, S$MON_CB_META_PATH)
      message(sprintf("[smoothing] draws_mon S=%d", draws_mon$S))
    }

    # Baseline check of the optional duration / monthly fits against the
    # occurrence fit (warn-only; a no-op when neither track is configured).
    .check_fit_pair_consistency(
      opt_occ,
      label = "smoothing",
      opt_dur = opt_dur,
      opt_mon = opt_mon
    )

    # NULL (hence fully inert) unless one of occ/dur/mon was actually fit with
    # `<track>.offset_baserate: true` -- see resolve_base_rate_paths().
    base_rate_table_paths <- resolve_base_rate_paths(
      list(occ = opt_occ, dur = opt_dur, mon = opt_mon),
      S$LF_CFG
    )
    sm_pop <- smoothing_frame_pop_args(S)
    frame <- build_location_frame(
      panel,
      hmm_clusters_path = if (file.exists(S$HMM_PATH)) S$HMM_PATH else NULL,
      size_data_path = S$SIZE_DATA_PATH,
      base_rate_table_paths = base_rate_table_paths,
      admin_level_data = surv_dat,
      gid_pop_path = sm_pop$gid_pop_path,
      wpp_path = sm_pop$wpp_path,
      pop_year = sm_pop$pop_year,
      pop_window = sm_pop$pop_window
    )
    ms_defaulted <- attr(frame, "min_seas_defaulted") # checked on the final frame
    max_draw_s <- draws_occ$S
  } else {
    reg <- load_grouping(S$GRP_SOURCE)
    country_to_group <- stats::setNames(reg$group, reg$iso3)

    # Representative-fit trick, same as utils/simulations.R's grouped forecast
    # loader: `frame$group` (hence the real `groups` list) doesn't exist until
    # AFTER build_location_frame() runs below, so probe every group in the
    # registry and take the first with a discoverable fit for that track. The
    # base_rate table is global per track (run-invariant filename -- see
    # persist_base_rate_table()), so any one deployed fit resolves the same path.
    # Only probe dur/mon when that track is actually part of this run.
    candidate_groups <- sort(unique(reg$group))
    repr_anchor_fit_opts <- list()
    for (.tr in c(
      "occ",
      if (isTRUE(S$GRP_DUR_ENABLE)) "dur",
      if (isTRUE(S$GRP_MON_ENABLE)) "mon"
    )) {
      for (.g in candidate_groups) {
        .o <- tryCatch(
          attr(
            readRDS(.discover_group_fit(
              S$GRP_FITS_DIR,
              S$GRP_SOURCE,
              .g,
              .tr,
              S$EPISODE_KEY %||% "lp_season",
              isTRUE(S$NATIONAL_ZERO_FILL),
              S$NATIONAL_ZERO_THRESHOLD %||% 0L,
              occ_panel = S$OCC_PANEL
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
    base_rate_table_paths <- resolve_base_rate_paths(
      repr_anchor_fit_opts,
      S$LF_CFG
    )
    sm_pop <- smoothing_frame_pop_args(S)

    frame <- build_location_frame(
      panel,
      hmm_clusters_path = if (file.exists(S$HMM_PATH)) S$HMM_PATH else NULL,
      size_data_path = S$SIZE_DATA_PATH,
      base_rate_table_paths = base_rate_table_paths,
      admin_level_data = surv_dat,
      gid_pop_path = sm_pop$gid_pop_path,
      wpp_path = sm_pop$wpp_path,
      pop_year = sm_pop$pop_year,
      pop_window = sm_pop$pop_window
    )
    ms_defaulted <- attr(frame, "min_seas_defaulted") # checked on the final frame
    frame$group <- unname(country_to_group[frame$country])
    unmapped <- unique(frame$country[is.na(frame$group)])
    if (length(unmapped) > 0L) {
      stop(sprintf(
        paste0(
          "[smoothing] %d country(ies) not in grouping registry '%s': %s. ",
          "Add them to analysis/groupings/%s.csv."
        ),
        length(unmapped),
        S$GRP_SOURCE,
        paste(unmapped, collapse = ", "),
        S$GRP_SOURCE
      ))
    }

    groups <- sort(unique(frame$group))
    message(sprintf(
      "[smoothing] grouped mode: %d group(s) [%s], %d locations",
      length(groups),
      paste(groups, collapse = ", "),
      nrow(frame)
    ))

    draws_occ_by_group <- list()
    cb_meta_occ_by_group <- list()
    if (isTRUE(S$GRP_DUR_ENABLE)) {
      draws_dur_by_group <- list()
      cb_meta_dur_by_group <- list()
    }
    if (isTRUE(S$GRP_MON_ENABLE)) {
      if (!isTRUE(S$GRP_DUR_ENABLE)) {
        stop(
          "[smoothing] simulator.grouping.monthly_enable requires ",
          "simulator.grouping.duration_enable -- the monthly case-count kernel ",
          "reads `active`/`mso`, which are duration-model latent states."
        )
      }
      draws_mon_by_group <- list()
      cb_meta_mon_by_group <- list()
    }
    for (g in groups) {
      occ_path <- .discover_group_fit(
        S$GRP_FITS_DIR,
        S$GRP_SOURCE,
        g,
        "occ",
        S$EPISODE_KEY %||% "lp_season",
        isTRUE(S$NATIONAL_ZERO_FILL),
        S$NATIONAL_ZERO_THRESHOLD %||% 0L,
        occ_panel = S$OCC_PANEL
      )
      message(sprintf("[smoothing] group '%s' occ:  %s", g, basename(occ_path)))

      occ_fit_g <- readRDS(occ_path)
      ep_fits[[paste0("occ:", g)]] <- occ_fit_g
      opt_occ_g <- .fit_opt_from_attr(occ_fit_g)
      if (is.null(opt_occ_g$occurrence)) {
        opt_occ_g$occurrence <- list(backend = "brms")
      }

      d_occ_g <- extract_occurrence_draws(
        occ_fit_g,
        opt_occ_g,
        max_draws = S$MAX_DRAWS_OCC
      )
      message(sprintf("[smoothing] group '%s': draws_occ S=%d", g, d_occ_g$S))

      draws_occ_by_group[[g]] <- d_occ_g
      cb_meta_occ_by_group[[g]] <- .load_cb_sidecar(occ_path)

      # Per-group duration fit (opt-in). Same discovery helper and same
      # hard-error-on-missing convention as occ above.
      if (isTRUE(S$GRP_DUR_ENABLE)) {
        dur_path <- .discover_group_fit(
          S$GRP_FITS_DIR,
          S$GRP_SOURCE,
          g,
          "dur",
          S$EPISODE_KEY %||% "lp_season"
        )
        message(sprintf(
          "[smoothing] group '%s' dur:  %s",
          g,
          basename(dur_path)
        ))
        dur_fit_g <- readRDS(dur_path)
        ep_fits[[paste0("dur:", g)]] <- dur_fit_g
        opt_dur_g <- .fit_opt_from_attr(dur_fit_g)
        if (is.null(opt_dur_g$duration)) {
          opt_dur_g$duration <- list(backend = "brms")
        }
        draws_dur_by_group[[g]] <- extract_duration_draws(
          dur_fit_g,
          opt_dur_g,
          max_draws = S$MAX_DRAWS_DUR
        )
        cb_meta_dur_by_group[[g]] <- .load_cb_sidecar(dur_path)
      }

      # Per-group monthly case-count fit (opt-in), same conventions again.
      if (isTRUE(S$GRP_MON_ENABLE)) {
        mon_path <- .discover_group_fit(
          S$GRP_FITS_DIR,
          S$GRP_SOURCE,
          g,
          "mon",
          S$EPISODE_KEY %||% "lp_season"
        )
        message(sprintf(
          "[smoothing] group '%s' mon:  %s",
          g,
          basename(mon_path)
        ))
        mon_fit_g <- readRDS(mon_path)
        ep_fits[[paste0("mon:", g)]] <- mon_fit_g
        opt_mon_g <- .fit_opt_from_attr(mon_fit_g)
        if (is.null(opt_mon_g$monthly)) {
          opt_mon_g$monthly <- list(backend = "brms_cb")
        }
        draws_mon_by_group[[g]] <- extract_monthly_draws(
          mon_fit_g,
          opt_mon_g,
          max_draws = S$MAX_DRAWS_MON
        )
        cb_meta_mon_by_group[[g]] <- .load_cb_sidecar(mon_path)
      }
    }

    draws_occ <- draws_occ_by_group[[groups[1]]]
    if (isTRUE(S$GRP_DUR_ENABLE)) {
      draws_dur <- draws_dur_by_group[[groups[1]]]
      cb_meta_dur <- cb_meta_dur_by_group[[groups[1]]]
    }
    if (isTRUE(S$GRP_MON_ENABLE)) {
      draws_mon <- draws_mon_by_group[[groups[1]]]
      cb_meta_mon <- cb_meta_mon_by_group[[groups[1]]]
    }
    max_draw_s <- min(vapply(draws_occ_by_group, function(d) d$S, integer(1)))
    message(sprintf(
      "[smoothing] grouped representative: occ S=%d  max_draw_s=%d",
      draws_occ$S,
      max_draw_s
    ))
  }

  # Zero-filled climate gate (simulator_core.R): no fit with a local-climate
  # term may be filtered over the panel's zero-filled post-DB months. Inert
  # without one. Window = burn-in start .. filter end (S$hist_window).
  .assert_no_zero_filled_climate(
    .draws_fe_names(
      draws_occ,
      draws_dur,
      draws_mon,
      draws_occ_by_group,
      draws_dur_by_group,
      draws_mon_by_group
    ),
    panel,
    window = c(S$BURNIN_START, S$FILTER_END),
    label = "smoothing"
  )

  # ── ADM1-ADM2 overlap filter ───────────────────────────────────────────────
  # build_location_frame() ingests every distinct panel gid, including coarse
  # ADM1-level units that also have nested ADM2 children (mixed-admin fits).
  # Keeping both double-counts geography AND injects whole-state ADM1 populations
  # that inflate the population-offset case means (mu = pop x exp(lp)).
  # Delegate to the shared drop_overlapping_adm1_from_frame() (utils.R) — the single
  # source of truth also used by the forecast loader (utils/simulations.R) — so
  # the two pipelines cannot drift on this drop.
  # Episode definition: the fits, the survival panel (obs_Y), the size frame
  # (population / admin-level fallback) and the monthly panel (obs_CM) must
  # share one episode digest whenever any of them -- or episodes.key -- is
  # gid_season (utils/provenance.R::sim_episode_check()).
  sim_episode_check(
    fits = ep_fits,
    preps = list(
      surv = S$SURV_PATH,
      size_data = S$SIZE_DATA_PATH,
      mon_data = S$MON_DATA_PATH,
      climate_panel = S$PANEL_PATH
    ),
    what = "smoothing artefacts",
    episode_key = S$EPISODE_KEY %||% "lp_season"
  )

  frame <- drop_overlapping_adm1_from_frame(frame, surv_dat, "smoothing")
  # The season onset gate needs every smoothed gid's own season start.
  check_onset_eligibility_frame(
    frame,
    S$ONSET_ELIGIBILITY %||% "inactive",
    "smoothing",
    defaulted = ms_defaulted
  )
  # TTLO level "none": each gid's first recorded onset, only
  # when the occurrence fit carries the level (the kernels seed from it).
  if (.draws_has_ttlo_none(if (S$GROUPED) draws_occ_by_group else draws_occ)) {
    frame$first_onset <- frame_first_onset(surv_dat, frame$gid)
  }

  if (all(is.na(frame$tu))) {
    message("[smoothing] all TU = NA — collapsing to _ALL_")
    frame$tu <- "_ALL_"
  }
  message(sprintf(
    "[smoothing] frame: %d locations, %d countries, TUs: %s",
    nrow(frame),
    length(unique(frame$country)),
    paste(sort(unique(stats::na.omit(frame$tu))), collapse = ", ")
  ))

  surv_panel_obs <- obs_panel_for_filter(
    surv_dat,
    list(monthly_case_obs = S$MONTHLY_CASE_OBS %||% "none"),
    S$MON_DATA_PATH,
    surv_path = S$SURV_PATH
  )

  # WHO dashboard term (simulator.filter.who_obs): the scored national monthly
  # counts plus, per frame row, the covered population share of its country
  # (who_f, after the ADM1/ADM2 overlap drop above so nothing is counted twice)
  # and its share of the country term (who_w).
  who_mon <- NULL
  who_scale <- NULL
  who_skipped <- NULL
  if (isTRUE(S$WHO_OBS$enable)) {
    if (!file.exists(S$WHO_OBS$path)) {
      stop(
        "[smoothing] who_obs.path not found: ",
        S$WHO_OBS$path,
        " -- run analysis/03_prepare_who_monthly.R."
      )
    }
    who_mon <- who_obs_table(
      readRDS(S$WHO_OBS$path),
      S$WHO_OBS$start,
      S$WHO_OBS$end
    )
    who_mon$cases_who <- round(who_mon$cases_who)
    # who_f's population year: the integer pop_year, or under a pop_year mode
    # the calendar year of the WHO window start, with the frame's population of
    # that year.
    who_py <- S[["POP_YEAR"]] %||% .SIM_POP_YEAR_DEFAULT
    fr_who <- frame
    if (sim_pop_year_is_mode(who_py)) {
      who_py <- as.integer(format(as.Date(S$WHO_OBS$start), "%Y"))
      fr_who$pop <- fr_who[[sprintf("pop_y%04d", who_py)]]
    }
    cov <- who_coverage_fraction(
      fr_who,
      load_wpp_table(.pop_artefact_or_null(S$WPP_PATH)),
      who_py
    )
    frame$who_f <- cov$who_f[match(as.character(frame$country), cov$country)]
    frame$who_w <- who_block_share(frame, S$WHO_OBS$block_share)
    # Per-country scale and DB-overlap rule: both read the
    # database side on the filter's own observation (surv_panel_obs, frame
    # gids after the overlap drop). Inert at scale_mode fixed + score_both.
    est_on <- identical(S$WHO_OBS[["scale_mode"]], "estimated")
    skip_on <- identical(S$WHO_OBS[["db_overlap"]], "skip_who")
    if (est_on || skip_on) {
      se <- S$WHO_OBS[["scale_est"]]
      db_months <- who_db_observed_months(
        surv_panel_obs,
        frame,
        min(S$WHO_OBS$start, if (est_on) se$start),
        max(S$WHO_OBS$end, if (est_on) se$end)
      )
      if (est_on) {
        who_est <- who_obs_table(readRDS(S$WHO_OBS$path), se$start, se$end)
        who_est$cases_who <- round(who_est$cases_who)
        who_scale <- estimate_who_scale(
          db_months,
          who_est,
          cov,
          se$start,
          se$end,
          min_db_cases = se$min_db_cases,
          min_months = se$min_months,
          shrink_cases = se$shrink_cases,
          fallback = se$fallback,
          fixed_scale = S$WHO_OBS$scale,
          clamp = se$clamp
        )
        frame$who_scale <- who_scale$scale[match(
          as.character(frame$country),
          who_scale$country
        )]
        message(sprintf(
          "[smoothing] who_obs scale (estimated, %s to %s): pooled %.3g; %s",
          format(se$start),
          format(se$end),
          who_scale$scale_pooled[1],
          paste(
            table(who_scale$source),
            names(table(who_scale$source)),
            sep = " ",
            collapse = ", "
          )
        ))
        print(
          as.data.frame(who_scale[, c(
            "country",
            "who_f",
            "est_months",
            "who_cases",
            "db_cases",
            "ratio_raw",
            "scale",
            "multiplier",
            "source"
          )]),
          row.names = FALSE,
          digits = 3
        )
      }
      if (skip_on) {
        n0 <- nrow(who_mon)
        who_mon <- who_drop_db_overlap(who_mon, db_months)
        who_skipped <- n0 - nrow(who_mon)
        message(sprintf(
          "[smoothing] who_obs db_overlap = skip_who: %d of %d WHO country-months dropped (database observes them)",
          who_skipped,
          n0
        ))
      }
    }
    no_who <- setdiff(unique(as.character(frame$country)), who_mon$country)
    message(sprintf(
      paste0(
        "[smoothing] who_obs: %d country-months scored (%s to %s), ",
        "%d countries; no WHO data (unscored): %s"
      ),
      nrow(who_mon),
      format(min(who_mon$t_start)),
      format(max(who_mon$t_start)),
      length(unique(who_mon$country)),
      if (length(no_who)) paste(sort(no_who), collapse = ", ") else "none"
    ))
    print(as.data.frame(cov), row.names = FALSE, digits = 3)
  }

  list(
    panel = panel,
    surv_dat = surv_dat,
    surv_panel_obs = surv_panel_obs,
    frame = frame,
    draws_occ = draws_occ,
    draws_occ_by_group = draws_occ_by_group,
    cb_meta_occ_by_group = cb_meta_occ_by_group,
    # NULL unless simulator.duration_fit (ungrouped) /
    # simulator.grouping.duration_enable (grouped) is set.
    draws_dur = draws_dur,
    draws_dur_by_group = draws_dur_by_group,
    cb_meta_dur_by_group = cb_meta_dur_by_group,
    cb_meta_dur = if (S$GROUPED) NULL else cb_meta_dur,
    # NULL unless simulator.monthly_fit (ungrouped) /
    # simulator.grouping.monthly_enable (grouped) is set. Both hard-require the
    # duration track above.
    draws_mon = draws_mon,
    draws_mon_by_group = draws_mon_by_group,
    cb_meta_mon_by_group = cb_meta_mon_by_group,
    cb_meta_mon = if (S$GROUPED) NULL else cb_meta_mon,
    country_to_group = country_to_group,
    max_draw_s = max_draw_s,
    # NULL unless simulator.filter.who_obs.enable.
    who_mon = who_mon,
    # NULL unless who_obs.scale_mode = estimated: the
    # estimate_who_scale() table; who_skipped = WHO rows dropped by skip_who.
    who_scale = who_scale,
    who_skipped = who_skipped
  )
}

# Stamp a `part` (partition) column on the frame and, in country mode, collapse
# the internal `tu` indicator to a single constant so the cross-country o_tu
# spillover term resolves to 0 (instead of NA) for every single-country
# partition. Local (o, by ADM1) and intra-country (o_ct) spillover are
# untouched. Must run AFTER S$PARTITION_BY is fully resolved (config + CLI +
# task override) and BEFORE smoothing_tu_list() / smoothing_build_precomp().
smoothing_apply_partition <- function(art, S) {
  fr <- art$frame
  if (identical(S$PARTITION_BY %||% "transmission_unit", "country")) {
    fr$part <- as.character(fr$country)
    # Preserve the true TU label before collapsing so the exogenous observed
    # o_tu term can be folded back into the occurrence LP (smoothing_build_precomp
    # -> .fold_observed_otu_into_lp). The simulated o_tu still resolves to 0.
    fr$tu_true <- as.character(fr$tu)
    # single constant TU collapses leave-one-country-out o_tu to 0; match type.
    fr$tu <- if (is.numeric(fr$tu)) 1L else "1"
    message(sprintf(
      "[smoothing] partition_by=country: %d country partition(s), o_tu collapsed",
      length(unique(stats::na.omit(fr$part)))
    ))
  } else if (identical(S$PARTITION_BY, "enso_region")) {
    # Region mode: one partition per ENSO group spanning all member-country ADM1
    # units, so build_spatpomp_tu_c()'s neighbour-contagion kernel (built over
    # frame_tu's ADM1 centroids) couples across country borders inside the
    # particle filter. Unlike country mode, `tu` is left INTACT: the
    # leave-one-country-out o_tu term resolves naturally for countries that share
    # an HMM TU and now sit in the same partition, so no collapse / tu_true /
    # observed-fold is needed (force_o_tu and otu_obs both stay inert).
    col <- assert_partition_column(fr, "enso_region") # "group"; errors w/ grouping.enable hint
    fr$part <- as.character(fr[[col]])
    message(sprintf(
      "[smoothing] partition_by=enso_region: %d region partition(s) (member countries coupled)",
      length(unique(stats::na.omit(fr$part)))
    ))
  } else if (identical(S$PARTITION_BY, "registry")) {
    # Registry mode: one partition per group of an arbitrary country-level
    # registry (filter.partition_source), independent of the fit grouping, so a
    # partition may mix fit groups -- each unit keeps its own `group` (draws,
    # cb_meta, lag) and the per-TU grouping subset routes them. `tu` is left
    # intact exactly as in enso_region mode (no o_tu collapse).
    fr <- stamp_partition_column(fr, "registry", S$PARTITION_SOURCE)
    message(sprintf(
      "[smoothing] partition_by=registry: source '%s', %d partition(s)",
      S$PARTITION_SOURCE,
      length(unique(stats::na.omit(fr$part)))
    ))
  } else {
    fr$part <- if (is.character(fr$tu)) fr$tu else as.character(fr$tu) # NA preserved
  }
  art$frame <- fr
  smoothing_preflight_partition(art, S)
  art
}

# Contagion registry of the smoothing run: simulator.nbr_contagion.
# partition_source, else the partition registry (as nbr_cfg_from_sim()).
.smoothing_nbr_partition_source <- function(S) {
  x <- S$NBR_PARTITION_SOURCE %||% ""
  if (nzchar(x)) x else (S$PARTITION_SOURCE %||% "")
}

# Partition pre-flight for the smoothing phase: the same consolidated check the
# forecast runs (spatpomp_core.R::check_partition_vs_grouping()), on a copy of
# the frame keyed by the smoothing partition (`part`). Inert pieces skip
# themselves: no grouping -> no fit-group checks; contagion off or no nbr FE ->
# no refinement check. Called at the end of smoothing_apply_partition().
smoothing_preflight_partition <- function(art, S) {
  if (is.null(art$frame) || is.null(art$frame$part)) {
    return(invisible(NULL))
  }
  fr <- art$frame
  fr$tu <- fr$part
  grp <- if (isTRUE(S$GROUPED) && !is.null(art$draws_occ_by_group)) {
    list(
      draws_occ_by_group = art$draws_occ_by_group,
      draws_dur_by_group = art$draws_dur_by_group
    )
  } else {
    NULL
  }
  nbr <- if (
    isTRUE(S$NBR_ENABLE) && (!is.null(art$draws_occ) || !is.null(grp))
  ) {
    resolve_nbr_spec(
      list(
        enable = TRUE,
        rho_km = S$NBR_RHO_KM,
        gadm_path = S$NBR_GADM_PATH,
        cache_dir = S$NBR_CACHE_DIR,
        grouping_source = S$GRP_SOURCE,
        partition_source = .smoothing_nbr_partition_source(S)
      ),
      art$draws_occ,
      grp
    )
  } else {
    NULL
  }
  check_partition_vs_grouping(
    fr,
    S$PARTITION_BY %||% "transmission_unit",
    grp,
    nbr
  )
}

# Full partition list (real partitions + synthetic NA bucket), honouring an
# optional subset. Keys off `frame$part` (= tu in TU mode, = country in country
# mode). The `_NA_TU_` bucket is only reachable in TU mode (NA tu); in country
# mode every row carries a country so no NA bucket arises.
smoothing_tu_list <- function(frame, tu_subset = NULL) {
  part <- if (!is.null(frame$part)) frame$part else frame$tu
  tus_real <- sort(unique(stats::na.omit(part)))
  tus <- if (any(is.na(part))) c(tus_real, "_NA_TU_") else tus_real
  if (!is.null(tu_subset)) {
    tus <- intersect(tus, tu_subset)
  }
  if (length(tus) == 0L) {
    stop("[smoothing] no partitions remaining after --tus filter")
  }
  tus
}

# ---------------------------------------------------------------------------
# Draw-index / seed assignment
# ---------------------------------------------------------------------------
# Reproduces the monolithic per-TU seeding so chunked runs match it exactly.
# Returns, for one TU, an integer vector of length N_REPS of occurrence draw
# indices (1..max_draw_s). `tus` MUST be the full sorted TU universe.
# S$DRAW_DESIGN (simulator.filter.draw_design): "per_tu_random" = that seeded
# per-unit sample; "shared_cycle" = replicate j runs draw j in every unit
# (spatpomp_core.R::smoothing_draw_idx(), shared with 07b's native filter).

smoothing_draw_idx_for_tu <- function(S, tu, tus, max_draw_s) {
  smoothing_draw_idx(
    S[["DRAW_DESIGN"]] %||% "per_tu_random",
    S$SEED_ROOT,
    tu,
    tus,
    S$N_REPS,
    max_draw_s
  )
}

# Per-(TU, rep) task seed, matching the monolithic flat task table
# (SEED_ROOT + global row index over uncached_tus × N_REPS).
smoothing_task_seeds <- function(S, tus) {
  rows <- do.call(
    rbind,
    lapply(tus, function(tu) {
      data.frame(tu = tu, rep_id = seq_len(S$N_REPS), stringsAsFactors = FALSE)
    })
  )
  rows$task_seed <- S$SEED_ROOT + seq_len(nrow(rows))
  rows
}

# ---------------------------------------------------------------------------
# Per-TU input pre-computation
# ---------------------------------------------------------------------------

smoothing_build_precomp <- function(S, art, tus) {
  lapply(stats::setNames(tus, tus), function(tu) {
    frame <- art$frame
    part <- if (!is.null(frame$part)) frame$part else frame$tu
    frame_tu <- if (identical(tu, "_NA_TU_")) {
      frame[is.na(part), , drop = FALSE]
    } else {
      frame[!is.na(part) & part == tu, , drop = FALSE]
    }
    if (nrow(frame_tu) == 0L) {
      return(NULL)
    }

    panel_tu <- build_climate_panel_for_tu(
      art$panel,
      frame_tu,
      horizon = S$hist_window,
      scenario = "factual",
      end_mode = "open"
    ) # half-open filter window [start, handoff) (utils/base.R)

    grouping_tu <- if (S$GROUPED) {
      grp_present <- sort(unique(as.character(frame_tu$group)))
      list(
        draws_occ_by_group = art$draws_occ_by_group[grp_present],
        cb_meta_occ_by_group = art$cb_meta_occ_by_group[grp_present],
        # NULL unless simulator.grouping.duration_enable; build_spatpomp_tu_c()
        # keys use_dur off this being non-NULL (spatpomp_core_c.R).
        draws_dur_by_group = if (is.null(art$draws_dur_by_group)) {
          NULL
        } else {
          art$draws_dur_by_group[grp_present]
        },
        # The duration track's own cross-basis. Same role as the occ/mon
        # entries above: each track's best ENSO lag is selected independently,
        # so without this the duration LP reads the occurrence basis.
        cb_meta_dur_by_group = if (is.null(art$cb_meta_dur_by_group)) {
          NULL
        } else {
          art$cb_meta_dur_by_group[grp_present]
        },
        # NULL unless simulator.grouping.monthly_enable; build_spatpomp_tu_c()
        # keys the mon path off the measurement, and errors if these are absent.
        draws_mon_by_group = if (is.null(art$draws_mon_by_group)) {
          NULL
        } else {
          art$draws_mon_by_group[grp_present]
        },
        cb_meta_mon_by_group = if (is.null(art$cb_meta_mon_by_group)) {
          NULL
        } else {
          art$cb_meta_mon_by_group[grp_present]
        }
      )
    } else {
      NULL
    }

    # Country-partition mode: observed cross-country o_tu activity by calendar
    # month, computed from the OTHER countries in this partition's true TU. Fed
    # to build_spatpomp_tu_c so .fold_observed_otu_into_lp can restore the
    # cross-country contagion term that the collapsed tu zeroes out. NULL (no-op)
    # in TU mode or when the partition's TU has no other countries.
    otu_obs <- NULL
    if (!is.null(frame$tu_true)) {
      this_tu <- unique(stats::na.omit(frame_tu$tu_true))
      this_country <- unique(as.character(frame_tu$country))
      if (length(this_tu) == 1L) {
        other_gids <- as.character(frame$gid[
          !is.na(frame$tu_true) &
            frame$tu_true == this_tu &
            !(as.character(frame$country) %in% this_country)
        ])
        if (length(other_gids) > 0L) {
          otu_obs <- art$surv_panel_obs %>%
            dplyr::filter(gid %in% other_gids) %>%
            dplyr::group_by(t_start) %>%
            dplyr::summarise(
              o_tu_obs = as.integer(any(obs_Y == 1L, na.rm = TRUE)),
              .groups = "drop"
            )
        }
      }
    }

    list(
      frame_tu = frame_tu,
      panel_tu = panel_tu,
      obs_data = art$surv_panel_obs,
      grouping = grouping_tu,
      # WHO term inputs for this partition (NULL unless who_obs.enable).
      who_obs = if (isTRUE(S$WHO_OBS$enable)) {
        c(
          list(
            tbl = art$who_mon[
              art$who_mon$country %in% as.character(frame_tu$country),
              ,
              drop = FALSE
            ],
            w = frame_tu$who_w,
            f = frame_tu$who_f
          ),
          S$WHO_OBS[c("phi", "mu_floor", "scale")],
          # eps of the `_actn_who` country mixture; inert otherwise.
          list(fp_rate = S$WHO_OBS[["fp_rate"]] %||% .SP_CASE_FP_RATE),
          # per-unit estimated scale; absent in fixed mode.
          if (!is.null(frame_tu[["who_scale"]])) {
            list(scale_u = frame_tu[["who_scale"]])
          }
        )
      } else {
        NULL
      },
      otu_obs = otu_obs,
      cb_meta_mon = if (S$GROUPED) NULL else art$cb_meta_mon,
      cb_meta_dur = if (S$GROUPED) NULL else art$cb_meta_dur
    )
  })
}

# ---------------------------------------------------------------------------
# Phase A: one bpfilter replicate
# ---------------------------------------------------------------------------
# cb_meta_occ IS passed (via pc$grouping$cb_meta_occ_by_group) so the ENSO
# cross-basis enters lp_static. For these brms_cb fits ENSO is encoded ONLY as
# enso_cb1..K (there is no plain enso_lag* fixed effect), so dropping cb_meta
# removed the ENSO term entirely and inflated the static hazard ~3x (the
# intercept is calibrated assuming the mean-negative cross-basis terms present).
# Burn-in/coverage NAs in the lag columns are imputed to 0 by make_lag_matrix() /
# build_cb_design(), so no NaN leaks into lp_static. Under grouping,
# build_spatpomp_tu_c routes each country to its group's lp components
# (drop_region=TRUE) and builds each unit's cross-basis from its group's cb_meta.

smoothing_filter_one <- function(
  S,
  pc,
  draws_occ,
  tu,
  rep_id,
  draw_d,
  seed,
  draws_dur = NULL,
  draws_mon = NULL
) {
  if (is.null(pc)) {
    return(NULL)
  }
  # Neighbour-contagion guard: config flag AND the (per-group) occurrence fit
  # carries an nbr_* FE. FALSE ⇒ no WNBR matvec in the generated Csnippet,
  # byte-identical to pre-change. Resolved here (not in smoothing_make_config)
  # because it needs the fit FEs, available only after artefact load.
  nbr_spec_i <- resolve_nbr_spec(
    list(
      enable = S$NBR_ENABLE,
      rho_km = S$NBR_RHO_KM,
      gadm_path = S$NBR_GADM_PATH,
      cache_dir = S$NBR_CACHE_DIR,
      grouping_source = S$GRP_SOURCE,
      partition_source = .smoothing_nbr_partition_source(S)
    ),
    draws_occ,
    pc$grouping
  )
  use_nbr_i <- nbr_spec_i$use
  # Outside the tryCatch: a WHO term on ADM1 blocks is a configuration error,
  # not a per-replicate filter failure to be logged and skipped.
  .check_who_obs_blocks(
    S$WHO_OBS %||% list(enable = FALSE),
    S$BLOCK_GRANULARITY %||% "country",
    pc$frame_tu
  )
  tryCatch(
    {
      obj_i <- build_spatpomp_tu_c(
        pc$frame_tu,
        pc$panel_tu,
        draws_occ,
        window = S$hist_window,
        n_draws_in_filter = 1L,
        fixed_draw_idx = draw_d,
        obs_data = pc$obs_data,
        grouping = pc$grouping,
        measurement = S$MEASUREMENT,
        # NULL for every pre-duration config, which keeps use_dur FALSE and the
        # `active`/`mso` state slots pinned at 0 -- byte-identical to the previous
        # behaviour of this function.
        draws_dur = if (S$GROUPED) NULL else draws_dur,
        # NULL unless the mon track is configured; the mon measurement enum is
        # what actually switches the kernel, and it hard-requires draws_dur.
        draws_mon = if (S$GROUPED) NULL else draws_mon,
        cb_meta_mon = if (S$GROUPED) NULL else pc$cb_meta_mon,
        # The duration cross-basis sidecar (without it the duration LP would
        # read the OCCURRENCE basis). NULL under grouping, where the per-group
        # maps travel inside `pc$grouping` instead.
        cb_meta_dur = if (S$GROUPED) NULL else pc$cb_meta_dur,
        otu_obs = pc$otu_obs,
        # Under partition_by=country the collapsed tu makes update_lags()'s o_tu
        # always 0; force it to 1 instead of relying on the observed-months fold
        # (build_spatpomp_tu_c skips .fold_observed_otu_into_lp when force_o_tu
        # is set, so this does not double-count with pc$otu_obs).
        force_o_tu = identical(
          S$PARTITION_BY %||% "transmission_unit",
          "country"
        ),
        # Neighbour-contagion: build the in-kernel WNBR matvec. Gated ONLY by
        # use_nbr_i (config flag AND fit carries the FE), independent of force_o_tu.
        use_nbr = use_nbr_i,
        nbr_rho_km = S$NBR_RHO_KM %||% 500L,
        nbr_gadm_path = S$NBR_GADM_PATH %||% "data/geodata/gadm_valid.gpkg",
        nbr_cache_dir = S$NBR_CACHE_DIR %||% "generated_data/",
        nbr_spec = nbr_spec_i,
        # Observed monthly cases (obs_CM): NULL unless monthly_case_obs != none.
        mon_case_col = S$MON_CASE_COL,
        covariate_feed = S$COVARIATE_FEED %||% "latent",
        # Detection-error rates: consumed only by `*_state_detect` measurements.
        delta_fn = S$DELTA_FN %||% .SP_STATE_DELTA_FN,
        delta_fp = S$DELTA_FP %||% .SP_STATE_DELTA_FP,
        # Inactive-state case mixture: consumed only by the `_actn` measurements.
        case_fp_rate = S$CASE_FP_RATE %||% .SP_CASE_FP_RATE,
        case_fp_mean = S$CASE_FP_MEAN %||% .SP_MON_MU_FLOOR,
        # Monthly incidence cap on the filter's mon mean (mon_mu_cap covariate).
        incidence_cap = S$INCIDENCE_CAP,
        # Onset gate (simulator.onset_eligibility).
        onset_eligibility = S$ONSET_ELIGIBILITY %||% "inactive",
        # Cessation-month exposure (simulator.cessation_exposure).
        cessation_exposure = S$CESSATION_EXPOSURE %||% "full",
        # WHO dashboard term (`_who` measurement only; NULL otherwise).
        who_obs = pc$who_obs
      )
      bl_i <- make_block_list(
        obj_i,
        pc$frame_tu,
        granularity = S$BLOCK_GRANULARITY %||% "country"
      )
      res_i <- run_one_bpfilter_replicate(
        obj_i,
        Np = S$NP,
        block_list = bl_i,
        seed = seed,
        traj_storage = S$TRAJ_STORAGE,
        prune_every = S$TRAJ_PRUNE_EVERY
      )
      list(
        tu = tu,
        rep_id = rep_id,
        draw_d = draw_d,
        filt_months = attr(obj_i, "filt_months"),
        filter_traj = res_i$filter_traj,
        logLik = res_i$logLik,
        # Canonical spatPomp unit ordering used for this smoothing run. Stored
        # so the forecaster can do gid-name-based state extraction when its
        # frame has fewer units (unseen_adm1 restriction differs between
        # smoothing and forecast artefact loaders).
        # method = "radix" pins this to byte-order (C-locale-equivalent)
        # collation, matching the dplyr::arrange(..., .locale = "C") that
        # built the pomp's unit axis in spatpomp_core.R. Under a system ICU
        # locale (e.g. en_US.UTF-8), plain sort() collates "BDI_1_1" before
        # "BDI_11_1"; radix/C collates "BDI_11_1" first. A mismatch here
        # silently swaps two units' state trajectories under their labels.
        gids = sort(as.character(pc$frame_tu$gid), method = "radix")
      )
    },
    error = function(e) {
      warning(sprintf(
        "[smoothing] TU '%s' rep %d: bpfilter failed: %s",
        tu,
        rep_id,
        conditionMessage(e)
      ))
      list(
        tu = tu,
        rep_id = rep_id,
        draw_d = draw_d,
        filt_months = NULL,
        filter_traj = NULL,
        logLik = NA_real_
      )
    }
  )
}

# Run Phase A for a chunk of replicates within one TU. `rep_ids`, `draw_idx`,
# `seeds` are parallel vectors. Uses furrr when workers > 1, else sequential.
smoothing_filter_chunk <- function(
  S,
  pc,
  draws_occ,
  tu,
  rep_ids,
  draw_idx,
  seeds,
  workers = 1L,
  draws_dur = NULL,
  draws_mon = NULL
) {
  stopifnot(
    length(rep_ids) == length(draw_idx),
    length(rep_ids) == length(seeds)
  )
  run_i <- function(i) {
    smoothing_filter_one(
      S,
      pc,
      draws_occ,
      tu,
      rep_ids[i],
      draw_idx[i],
      seeds[i],
      draws_dur = draws_dur,
      draws_mon = draws_mon
    )
  }

  if (workers > 1L && requireNamespace("furrr", quietly = TRUE)) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    furrr::future_map(
      seq_along(rep_ids),
      run_i,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    lapply(seq_along(rep_ids), run_i)
  }
}

# ---------------------------------------------------------------------------
# Phase B: smoothing-trajectory extraction + mon-track NegBin monthly
# case-count regeneration. FILTER_START/END passed
# explicitly.
# ---------------------------------------------------------------------------

#' Pull a U-row block of per-unit state rows out of a filter_traj slice,
#' falling back to a zero matrix when the slot is absent.
#'
#' The `active`/`mso` slots are part of the state vector
#' (spatpomp_core.R .SP_UNIT_STATENAMES), and `check_smoothing_cache_schema()`
#' rejects caches without them. This fallback therefore only fires for
#' hand-built test fixtures; it exists so a missing slot degrades to "duration
#' never fired" rather than a subscript error deep inside the month loop.
.ft_rows_or_zero <- function(ft_slice, stem, U) {
  rn <- paste0(stem, seq_len(U))
  if (!all(rn %in% rownames(ft_slice))) {
    return(matrix(0, nrow = U, ncol = ncol(ft_slice)))
  }
  ft_slice[rn, , drop = FALSE]
}

#' Reconstruct "active DURING month m" (onset OR carried over from a prior
#' month) directly from one raw states.qs record, for the ongoing-outbreak
#' ribbon (08c_smoothing_compile.R). Unlike the sparse Yout onset table, an
#' "ongoing outbreaks" count needs every month an episode is running, not just
#' its onset month -- so this reads the ACTIVE state trajectory directly out
#' of the persisted filter_traj rather than reusing smoothing_events_long().
#'
#' Mirrors the act_prev_mat/act_m reconstruction inside
#' extract_smoothing_one_tu() (same "entering month m" offset-by-one
#' rationale documented there: the cessation hazard can fire and zero `active`
#' in the same rprocess step that writes month m's state, so a genuinely
#' active final month reads as inactive unless combined with that month's
#' Yout). Operates on a raw states.qs record because 08c never rebuilds the
#' full per-TU smoothing object extract_smoothing_one_tu() returns -- states.qs
#' already carries everything needed (filter_traj, filt_months, gids).
#'
#' @param rec one element of a states.qs list: list(tu, rep_id, draw_d,
#'   filt_months, filter_traj, logLik, gids).
#' @return matrix (U x M), 0/1 "active during month m", or NULL when
#'   filter_traj is absent (mirrors extract_smoothing_one_tu()'s NULL-skip).
#'   All-zero when the duration track wasn't configured for this run (the
#'   `active` state slot is pinned at 0 -- see .ft_rows_or_zero()), so the
#'   ongoing-outbreak ribbon degrades gracefully to an all-zero band rather
#'   than erroring.
#' The per-month slice (U x M) of one state stem from a states.qs record, or
#' NULL when the record has no trajectory or the stem is absent (e.g. `det` in
#' a latent-feed run). Column m is the state WRITTEN at filter month m.
.smoothing_state_matrix <- function(rec, stem) {
  if (is.null(rec$filter_traj)) {
    return(NULL)
  }
  ft <- rec$filter_traj
  U <- length(rec$gids)
  M <- length(rec$filt_months)
  rn <- paste0(stem, seq_len(U))
  if (!all(rn %in% dimnames(ft)[[1]])) {
    return(NULL)
  }
  toff <- dim(ft)[3] - M
  matrix(ft[rn, 1L, (toff + 1L):dim(ft)[3]], nrow = U, ncol = M)
}

#' Detected activity during month m (covariate_feed: detected): the `det`
#' state written at month m IS the detected-activity field for that month
#' (update_detection() runs before cessation, so an episode's final month
#' keeps its flag). NULL for latent-feed runs.
.smoothing_detected_matrix <- function(rec) {
  d <- .smoothing_state_matrix(rec, .SP_DET_STATENAME)
  if (is.null(d)) NULL else matrix(as.integer(round(d)), nrow = nrow(d))
}

.smoothing_active_matrix <- function(rec) {
  if (is.null(rec$filter_traj)) {
    return(NULL)
  }
  ft <- rec$filter_traj
  U <- length(rec$gids)
  M <- length(rec$filt_months)

  toff <- dim(ft)[3] - M
  ft_slice <- ft[, 1L, (toff + 1L):dim(ft)[3]]
  if (!is.matrix(ft_slice)) {
    ft_slice <- matrix(
      ft_slice,
      ncol = 1L,
      dimnames = list(names(ft_slice), NULL)
    )
  }
  ft_wide <- ft[, 1L, toff:dim(ft)[3]]
  if (!is.matrix(ft_wide)) {
    ft_wide <- matrix(ft_wide, ncol = 1L, dimnames = list(names(ft_wide), NULL))
  }

  Yout_mat <- ft_slice[paste0("Yout", seq_len(U)), , drop = FALSE]
  act_prev_mat <- .ft_rows_or_zero(ft_wide, "active", U)[,
    seq_len(M),
    drop = FALSE
  ]

  matrix(
    as.integer((Yout_mat == 1L) | (act_prev_mat == 1L)),
    nrow = U,
    ncol = M
  )
}

extract_smoothing_one_tu <- function(
  tu_results,
  frame_tu,
  panel_tu,
  FILTER_START,
  FILTER_END,
  grouping = NULL,
  draws_mon = NULL,
  cb_meta_mon = NULL,
  incidence_cap = NULL,
  cessation_exposure = "full"
) {
  if (is.null(tu_results) || length(tu_results) == 0L) {
    return(NULL)
  }
  # simulator.cessation_exposure: half halves the post-hoc mean of a unit that
  # ceased in the month, as the C filter's mon block does.
  cess_half <- identical(
    check_cessation_exposure(
      cessation_exposure,
      "extract_smoothing_one_tu: cessation_exposure"
    ),
    "half"
  )
  U <- nrow(frame_tu)

  filt_months <- NULL
  for (r in tu_results) {
    if (!is.null(r$filt_months)) {
      filt_months <- r$filt_months
      break
    }
  }
  if (is.null(filt_months)) {
    all_m <- sort(unique(panel_tu$t_start))
    filt_months <- all_m[all_m >= FILTER_START & all_m < FILTER_END]
  }
  M <- length(filt_months)
  if (M == 0L) {
    return(NULL)
  }

  # method = "radix" — must match the C-locale unit axis built in
  # spatpomp_core.R (see comment there); see the .write_partition `gids`
  # comment above for the failure mode.
  sorted_gids <- sort(as.character(frame_tu$gid), method = "radix")
  frame_sorted <- frame_tu[match(sorted_gids, as.character(frame_tu$gid)), ]
  # Per-unit monthly incidence ceiling, frame_sorted order (Inf when NULL).
  mu_cap <- monthly_mu_cap(frame_sorted$pop, incidence_cap)
  # pop_year mode: per-month population and ceiling, NULL for a
  # frame without pop_y<YYYY> columns (the default keeps the static ones).
  pop_mat_s <- frame_pop_matrix(frame_sorted, filt_months)
  mu_cap_mat_s <- if (!is.null(pop_mat_s)) {
    matrix(
      monthly_mu_cap(as.vector(pop_mat_s), incidence_cap),
      nrow = nrow(pop_mat_s)
    )
  } else {
    NULL
  }

  # ── Monthly case-count (`mon`) draws bundles ──────────────────────────────
  # The only track regenerated in Phase B. `has_mon` FALSE (no monthly fit
  # configured) leaves every case cell NA; the onset table is still emitted.
  mon_draws_by_group <- if (!is.null(grouping)) {
    grouping$draws_mon_by_group %||% list()
  } else if (!is.null(draws_mon)) {
    list(`_ALL_` = draws_mon)
  } else {
    list()
  }
  has_mon <- length(mon_draws_by_group) > 0L

  # Plain `enso_lag<N>` fixed effects consumed by the mon LP (non-cross-basis
  # mon fits). Sourced from the mon bundles' own active terms, so a plain-ENSO
  # mon fit is never silently ENSO-blind here. Deployed mon fits
  # are brms_cb (ENSO enters as enso_cb*, rebuilt from cb_meta below), so this
  # set is empty for them.
  # enso_lag<N> and (deployed ENSO leads) enso_lead<k>.
  panel_enso_cols <- grep(
    "^enso_(lag|lead)[0-9]+$",
    names(panel_tu),
    value = TRUE
  )
  mon_enso_cols <- intersect(
    panel_enso_cols,
    unique(unlist(lapply(mon_draws_by_group, function(d) d$active_terms)))
  )

  panel_indexed <- panel_tu %>%
    dplyr::filter(t_start %in% filt_months) %>%
    dplyr::mutate(
      .gid_idx = match(as.character(gid), sorted_gids),
      .m_idx = match(as.Date(t_start), filt_months)
    ) %>%
    dplyr::filter(!is.na(.gid_idx), !is.na(.m_idx))

  build_mat <- function(col, default = 0) {
    M_mat <- matrix(default, nrow = U, ncol = M)
    if (col %in% names(panel_indexed) && nrow(panel_indexed) > 0L) {
      M_mat[cbind(
        panel_indexed$.gid_idx,
        panel_indexed$.m_idx
      )] <- panel_indexed[[col]]
    }
    M_mat
  }

  clim_mats <- list(
    mean_precip = build_mat("mean_precip"),
    dev_precip = build_mat("dev_precip"),
    mean_temp_daymean = build_mat("mean_temp_daymean"),
    dev_temp_daymean = build_mat("dev_temp_daymean")
  )
  for (ec in mon_enso_cols) {
    clim_mats[[ec]] <- build_mat(ec)
  }

  # ENSO is location-invariant: broadcast one month-representative value (the
  # largest-magnitude, i.e. non-NA/non-imputed, cell) across units.
  for (ec in mon_enso_cols) {
    mat <- clim_mats[[ec]]
    month_vals <- apply(mat, 2L, function(cv) cv[which.max(abs(cv))])
    clim_mats[[ec]] <- matrix(month_vals, nrow = U, ncol = M, byrow = TRUE)
  }

  if (nrow(panel_indexed) > 0L) {
    has_data <- matrix(FALSE, nrow = U, ncol = M)
    has_data[cbind(panel_indexed$.gid_idx, panel_indexed$.m_idx)] <- TRUE
    for (cc in c(
      "mean_precip",
      "dev_precip",
      "mean_temp_daymean",
      "dev_temp_daymean"
    )) {
      if (!(cc %in% names(clim_mats))) {
        next
      }
      mat <- clim_mats[[cc]]
      for (i in seq_len(U)) {
        real_ms <- which(has_data[i, ])
        if (length(real_ms) == 0L || max(real_ms) == M) {
          next
        }
        last_m <- max(real_ms)
        mat[i, (last_m + 1L):M] <- mat[i, last_m]
      }
      clim_mats[[cc]] <- mat
    }
  }

  # ── Mon DLNM cross-basis (brms_cb fits) ───────────────────────────────────
  # brms_cb mon fits take ENSO only as enso_cb1..K columns, so rebuild each
  # group's cross-basis from its own cb_meta sidecar (build_cb_design ->
  # make_lag_matrix -> build_cb_columns_from_Q, as the forward simulator does)
  # and inject it into xrow. ENSO is location-invariant: one value per lag
  # column, broadcast across units, keyed by group. `mso` uses build_mso_basis()
  # instead. ONI is clamped to cbm$enso$range before evaluation, so the spline
  # is never extrapolated past the training range.
  enso_month_vec <- function(col, bk = NULL) {
    mat <- build_mat(col)
    vals <- apply(mat, 2L, function(cv) cv[which.max(abs(cv))])
    if (!is.null(bk)) {
      vals <- .clamp_cb_range(vals, bk)
    }
    vals
  }
  .enso_lag_cols_for_meta <- function(cbm) {
    if (is.null(cbm) || is.null(cbm[["enso"]])) {
      return(character(0))
    }
    meta <- cbm[["enso"]]
    # base.R::enso_lag_col_name(): a negative single_lag is a LEAD (enso_lead<k>),
    # the column build_cb_design() reads for it (dlnm.R::.cb_lag_exposure()).
    if (!is.null(meta$single_lag)) {
      enso_lag_col_name(as.integer(meta$single_lag))
    } else if (!is.null(meta$lags)) {
      vapply(as.integer(meta$lags), enso_lag_col_name, character(1))
    } else {
      sprintf("enso_lag%d", seq(0L, as.integer(meta$lag[2])))
    }
  }

  mon_cb_meta_by_group <- if (!is.null(grouping)) {
    grouping$cb_meta_mon_by_group %||% list()
  } else if (!is.null(cb_meta_mon)) {
    list(`_ALL_` = cb_meta_mon)
  } else {
    list()
  }

  mon_cb_by_group <- list() # group -> M × ncb ENSO design matrix
  mso_basis_by_grp <- list() # group -> list(B, mso_max)
  if (has_mon) {
    for (g in names(mon_cb_meta_by_group)) {
      cbm <- mon_cb_meta_by_group[[g]]
      if (is.null(cbm) || length(cbm) == 0L) {
        next
      }
      mso_basis_by_grp[[g]] <- build_mso_basis(cbm)
      if (is.null(cbm[["enso"]])) {
        next
      }
      lag_cols <- .enso_lag_cols_for_meta(cbm)
      bk <- cbm[["enso"]]$range
      df_months <- as.data.frame(
        stats::setNames(
          lapply(lag_cols, function(lc) enso_month_vec(lc, bk = bk)),
          lag_cols
        ),
        check.names = FALSE
      )
      cbdes <- build_cb_design(df_months, cbm, vars = "enso")
      if (!is.null(cbdes)) mon_cb_by_group[[g]] <- as.matrix(cbdes)
    }
  }

  # Monthly exposure offset. The fitted panel carries `log_offset` =
  # log(exposure days); in smoothing there is no observed exposure column, so it
  # is synthesised as log(days_in_month(t)) -- the same approximation the v1 and
  # spatPomp kernels make (simulator_core.R). That is exact for an episode's
  # INTERIOR months. The ONSET month is corrected at the consumer below via
  # .onset_log_offset_adj(); the CESSATION month is still given full-month
  # exposure, matching the forward kernels. See .ONSET_EXPOSURE_FRAC
  # (simulator_core.R) for the defect's measured size and why onset only.
  log_offset_vec <- log(as.numeric(lubridate::days_in_month(filt_months)))
  mon_month_num <- as.numeric(lubridate::month(filt_months))
  # Canonical key, matching the P3-re-keyed factor the mon fit was levelled on
  # (.frame_admin_level_canonical(), simulator_core.R).
  adm2_dummy <- if ("admin_level" %in% names(frame_sorted)) {
    as.numeric(.frame_admin_level_canonical(frame_sorted) == "ADM2")
  } else {
    NULL
  }

  t_seas_mat <- do.call(
    rbind,
    lapply(as.integer(frame_sorted$min_seas), function(ms) {
      as.integer(
        lubridate::yday(
          filt_months %m+% lubridate::period(-(ms - 1L), "months")
        ) -
          1L
      )
    })
  )

  sentinel <- .SP_OTU_SENTINEL
  R <- length(tu_results)
  out_list <- vector("list", R)

  for (r_idx in seq_len(R)) {
    res_r <- tu_results[[r_idx]]
    if (is.null(res_r$filter_traj)) {
      out_list[[r_idx]] <- NULL
      next
    }

    ft <- res_r$filter_traj
    rep_id <- res_r$rep_id
    draw_d <- res_r$draw_d

    # bpfilter(filter_traj=TRUE) prepends the t0 initial state as the first time
    # column, so dim(ft)[3] == M + 1.  Drop that column so ft_slice[, m] aligns
    # with filt_months[m]; reading columns 1..M without this skip shifts every
    # occurrence event one month early and drops the final filtered month
    # (identical logLik but collapsed obs1 recovery — see diag_reproduce_one_task.R).
    toff <- dim(ft)[3] - M
    ft_slice <- ft[, 1L, (toff + 1L):dim(ft)[3]]
    if (!is.matrix(ft_slice)) {
      ft_slice <- matrix(
        ft_slice,
        ncol = 1L,
        dimnames = list(names(ft_slice), NULL)
      )
    }

    Yout_mat <- ft_slice[paste0("Yout", seq_len(U)), , drop = FALSE]
    ud_mat <- ft_slice[paste0("u_days", seq_len(U)), , drop = FALSE]
    o_mat <- ft_slice[paste0("o", seq_len(U)), , drop = FALSE]
    oct_mat <- ft_slice[paste0("o_ct", seq_len(U)), , drop = FALSE]
    otu_mat <- ft_slice[paste0("o_tu", seq_len(U)), , drop = FALSE]
    # Duration latent states: always in filter_traj, all-zero unless the filter
    # ran with draws_dur, so reading them is safe.
    # The recorded active/mso for month m are the state heading into m+1 (the
    # kernel resets them in the same step a cessation fires). So take the state
    # ENTERING month m (ft_wide: one column wider, t0 prepended) and combine it
    # with the onset flag Yout_mat[, m] to recover "active during month m", the
    # same av_new predicate mon_mu_block gates mu on.
    ft_wide <- ft[, 1L, toff:dim(ft)[3]]
    if (!is.matrix(ft_wide)) {
      ft_wide <- matrix(
        ft_wide,
        ncol = 1L,
        dimnames = list(names(ft_wide), NULL)
      )
    }
    act_prev_mat <- .ft_rows_or_zero(ft_wide, "active", U)[,
      seq_len(M),
      drop = FALSE
    ]
    mso_prev_mat <- .ft_rows_or_zero(ft_wide, "mso", U)[,
      seq_len(M),
      drop = FALSE
    ]
    # Recorded active for month m itself (heading into m + 1). A unit active
    # DURING month m whose recorded active is 0 ceased in month m -- the
    # same predicate (av_act = 1, av_new = 0) the C filter's mon block halves
    # under cessation_exposure: half. Read only under half.
    act_rec_mat <- if (cess_half) {
      .ft_rows_or_zero(ft_wide, "active", U)[, 1L + seq_len(M), drop = FALSE]
    } else {
      NULL
    }

    out_Yout <- matrix(0L, nrow = U, ncol = M)
    # Monthly case-count channel: NA everywhere unless the mon track is
    # configured, in which case it carries a draw at every ACTIVE month.
    out_cases_mon <- matrix(NA_real_, nrow = U, ncol = M)
    # `_who` measurement: the filter itself drew this month's cases per unit
    # (state `cm`, the C twin of draw_monthly_cases()) and scored their country
    # sums against the WHO counts. Those are the conditioned draws, so they ARE
    # the case channel; redrawing from mu here would discard the conditioning.
    # NULL for every other measurement (the redraw below runs as before).
    cm_mat <- .smoothing_state_matrix(res_r, "cm")

    for (m in seq_len(M)) {
      yout_m <- as.integer(round(Yout_mat[, m]))
      u_days_m <- ud_mat[, m]
      o_m <- as.integer(round(o_mat[, m]))
      o_ct_m <- as.integer(round(oct_mat[, m]))
      o_tu_raw <- otu_mat[, m]
      o_tu_m <- ifelse(
        o_tu_raw <= sentinel + 0.5,
        NA_integer_,
        as.integer(round(o_tu_raw))
      )

      out_Yout[, m] <- yout_m

      # ── Monthly case-count draw (`mon` track) ────────────────────────────
      # Gated on activity DURING this month (reconstructed below), not on
      # Yout alone, so every month of an ongoing episode emits a case count:
      # a month with no onset anywhere can still have active episodes carried
      # over from an earlier month.
      if (has_mon) {
        act_prev_m <- as.integer(round(act_prev_mat[, m]))
        mso_prev_m <- as.integer(round(mso_prev_mat[, m]))
        # Active DURING month m (pre-cessation): onset this month, or already
        # active entering this month -- see the act_prev_mat comment above.
        act_m <- as.integer((yout_m == 1L) | (act_prev_m == 1L))
        mso_m <- ifelse(yout_m == 1L, 1L, mso_prev_m + 1L)
        ceased_m <- if (cess_half) {
          act_m == 1L & as.integer(round(act_rec_mat[, m])) == 0L
        } else {
          NULL
        }
        act_i <- which(act_m == 1L)
        if (length(act_i) > 0L && !is.null(cm_mat)) {
          out_cases_mon[act_i, m] <- cm_mat[act_i, m]
        }
        if (length(act_i) > 0L && is.null(cm_mat)) {
          xrow_mon <- lapply(clim_mats, function(mat) mat[, m])
          xrow_mon[["t_start_seas"]] <- t_seas_mat[, m]
          xrow_mon <- lapply(xrow_mon, function(v) {
            v[is.na(v)] <- 0
            v
          })
          # monthly.seasonal_term: "harmonic2" (the default) names these four
          # columns; load_monthly_data() builds them the same way at fit time.
          # Harmless under seasonal_term: none -- evaluate_lp() drops any term
          # absent from `fe`.
          mnum <- mon_month_num[m]
          xrow_mon[["sin_month"]] <- rep(sin(2 * pi * mnum / 12), U)
          xrow_mon[["cos_month"]] <- rep(cos(2 * pi * mnum / 12), U)
          xrow_mon[["sin2_month"]] <- rep(sin(4 * pi * mnum / 12), U)
          xrow_mon[["cos2_month"]] <- rep(cos(4 * pi * mnum / 12), U)
          if (!is.null(adm2_dummy)) {
            xrow_mon[["admin_levelADM2"]] <- adm2_dummy
          }

          log_mu_mon <- rep(NA_real_, U)
          phi_mon <- rep(NA_real_, U)
          grp_sorted <- if (is.null(grouping)) {
            rep("_ALL_", U)
          } else {
            as.character(frame_sorted$group)
          }
          for (g in unique(grp_sorted)) {
            d_g <- mon_draws_by_group[[g]]
            if (is.null(d_g)) {
              next
            }
            idx_g <- which(grp_sorted == g)
            s_g <- ((as.integer(draw_d) - 1L) %% d_g$S) + 1L
            xm <- lapply(xrow_mon, function(v) v[idx_g])
            cb_g <- mon_cb_by_group[[g]]
            if (!is.null(cb_g)) {
              for (cn in colnames(cb_g)) {
                xm[[cn]] <- rep(cb_g[m, cn], length(idx_g))
              }
            }
            fs_g <- frame_sorted[idx_g, , drop = FALSE]
            # Bundle-driven offsets (monthly_log_mu()). For a mon bundle these
            # resolve to log(pop) + log_offset (+ the base_rate_mon anchor).
            # The onset-month exposure correction reads `mso_m`, reconstructed
            # above as ifelse(yout_m == 1, 1, mso_prev + 1); on the ACTIVE rows
            # this block scores, mso_m == 1 holds exactly for units onsetting
            # this month (a unit active entering month m has mso_prev >= 1,
            # hence mso_m >= 2), so it is the same predicate as yout_m == 1 and
            # matches the forward kernels' `mso == 1` test. Only this kernel
            # floors pop at machine epsilon before the log.
            lm_g <- monthly_log_mu(
              d_g,
              s_g,
              xm,
              o_m[idx_g],
              o_ct_m[idx_g],
              o_tu_m[idx_g],
              u_days_m[idx_g],
              fs_g,
              log_offset = log_offset_vec[m],
              mso = mso_m[idx_g],
              mso_tab = mso_basis_by_grp[[g]],
              what = "monthly kernel (smoothing)",
              drop_region = !is.null(grouping),
              pop_floor = .Machine$double.eps,
              pop = if (!is.null(pop_mat_s)) pop_mat_s[idx_g, m] else NULL,
              ceased = if (cess_half) ceased_m[idx_g] else NULL
            )
            log_mu_mon[idx_g] <- lm_g
            phi_mon[idx_g] <- compute_phi_vec(
              d_g,
              s_g,
              fs_g,
              drop_region = !is.null(grouping)
            )
          }
          # Monthly incidence cap (simulator.incidence_cap) on the mean and
          # the drawn count, exactly as simulate_one_tu(); Inf when NULL.
          if (!is.null(mu_cap_mat_s)) {
            mu_cap <- mu_cap_mat_s[, m]
          }
          mu_a <- exp(log_mu_mon[act_i])
          if (any(is.finite(mu_cap[act_i]))) {
            mu_a <- pmin(mu_a, mu_cap[act_i])
          }
          cm <- draw_monthly_cases(mu_a, phi_mon[act_i], mu_cap[act_i])
          ok <- !is.na(cm)
          if (any(ok)) out_cases_mon[act_i[ok], m] <- cm[ok]
        }
      }
    }

    out_list[[r_idx]] <- list(
      rep_id = rep_id,
      Yout_mat = out_Yout,
      # All-NA unless the mon track is configured.
      cases_mon_mat = out_cases_mon,
      # Detected-activity field (covariate_feed:
      # detected); NULL for latent-feed runs, which keep
      # their event tables byte-identical.
      det_mat = .smoothing_detected_matrix(res_r)
    )
  }

  list(
    filt_months = filt_months,
    sorted_gids = sorted_gids,
    frame_sorted = frame_sorted,
    reps = out_list
  )
}

#' .smoothing_events_long_impl
#'
#' Shared body of smoothing_events_long() / smoothing_events_mon_long().
#' `kind = "onset"`: one row per (rep, gid, month) with
#' Yout == 1 and a constant Yout = 1L column. `kind = "mon"`: one row per active
#' month with a case draw (non-NA cases_mon_mat) carrying cases_month;
#' replicates without that matrix are skipped and an empty result is NULL.
#' Both carry `detected` only under the detected covariate feed.
#' @noRd
.smoothing_events_long_impl <- function(sm, kind = c("onset", "mon")) {
  kind <- match.arg(kind)
  if (is.null(sm) || length(sm$reps) == 0L) {
    return(NULL)
  }
  filt_months <- sm$filt_months
  sorted_gids <- sm$sorted_gids
  frame_sorted <- sm$frame_sorted
  is_mon <- identical(kind, "mon")
  out <- purrr::map_dfr(sm$reps, function(rep_r) {
    if (is.null(rep_r) || (is_mon && is.null(rep_r$cases_mon_mat))) {
      return(NULL)
    }
    idx <- if (is_mon) {
      which(!is.na(rep_r$cases_mon_mat), arr.ind = TRUE)
    } else {
      which(rep_r$Yout_mat == 1L, arr.ind = TRUE)
    }
    if (nrow(idx) == 0L) {
      return(NULL)
    }
    tb <- tibble::tibble(
      rep_id = rep_r$rep_id,
      t_start = filt_months[idx[, 2L]],
      gid = sorted_gids[idx[, 1L]],
      AFRO_region = frame_sorted$AFRO_region[idx[, 1L]],
      country = frame_sorted$country[idx[, 1L]]
    )
    if (is_mon) {
      tb$cases_month <- rep_r$cases_mon_mat[idx]
    } else {
      tb$Yout <- 1L
    }
    # Only under covariate_feed: detected (no column otherwise).
    if (!is.null(rep_r$det_mat)) {
      tb$detected <- rep_r$det_mat[idx]
    }
    tb
  })
  if (is_mon && (is.null(out) || nrow(out) == 0L)) NULL else out
}

# Assemble a sparse long event table (only Yout==1 cells) from one TU's
# extracted smoothing object. Columns: rep_id, t_start, gid, AFRO_region,
# country, Yout (+ detected under covariate_feed: detected). Case counts live
# in the separate mon-track table, smoothing_events_mon_long().
# (Scratchpad Section 7, per-TU slice.)
smoothing_events_long <- function(sm) .smoothing_events_long_impl(sm, "onset")


# Sparse long table for the monthly case-count (`mon`) channel: one row per
# (rep, gid, month) cell that carries a drawn case count, i.e. every ACTIVE
# month rather than only onset months.
#
# This is a SECOND table rather than a relaxed filter on smoothing_events_long()
# because the two have different row sets: onsets there, active months here.
# Returns NULL when the mon track is not configured (every cell NA), so callers
# can treat it as opt-in.
smoothing_events_mon_long <- function(sm) .smoothing_events_long_impl(sm, "mon")


# ---------------------------------------------------------------------------
# Compile-side helpers
# ---------------------------------------------------------------------------

# wquantile(): utils/base.R (one copy, shared with utils/results_derive.R).

# Softmax importance weights from per-(TU, rep) logLik, summed across TUs by
# rep_id. uniform=TRUE assigns 1/n. Returns tibble(rep_id, w) + prints ESS.
smoothing_rep_weights <- function(rep_logliks, uniform = FALSE) {
  # A replicate whose bpfilter failed in some TU (logLik NA, no trajectory
  # there) gets weight 0: summing the remaining TUs with na.rm would drop a
  # large negative term and hand it nearly all the weight.
  failed <- unique(rep_logliks$rep_id[is.na(rep_logliks$logLik)])
  if (length(failed)) {
    warning(
      sprintf(
        "[smoothing] %d replicate(s) with a failed bpfilter in some TU get weight 0: %s",
        length(failed),
        paste(head(sort(failed), 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (uniform) {
    rw <- rep_logliks %>%
      dplyr::distinct(rep_id) %>%
      dplyr::mutate(w = ifelse(rep_id %in% failed, 0, 1), w = w / sum(w))
    message(sprintf(
      "[smoothing] importance weighting DISABLED: uniform 1/%d.",
      sum(rw$w > 0)
    ))
    return(rw)
  }
  rw <- rep_logliks %>%
    dplyr::group_by(rep_id) %>%
    dplyr::summarise(logLik_total = sum(logLik), .groups = "drop") %>%
    dplyr::mutate(
      log_w = logLik_total - max(logLik_total, na.rm = TRUE),
      w = ifelse(is.na(log_w), 0, exp(log_w)),
      w = w / sum(w)
    ) %>%
    dplyr::select(rep_id, w)
  ess_val <- 1 / sum(rw$w^2)
  message(sprintf(
    "[smoothing] importance weights: ESS = %.1f / %d (%.0f%%)",
    ess_val,
    nrow(rw),
    100 * ess_val / nrow(rw)
  ))
  rw
}

# Surveillance-coverage mask (07b parity). Canonical definition lives in
# utils/sim_io.R — sourced directly in this file's source block
# above (not via the utils.R facade, which this module does not source).

#' .agg_and_ribbon_impl
#'
#' Shared body of agg_and_ribbon() (onsets: Yout -> n_outbreaks),
#' agg_and_ribbon_mon() (cases_month -> n_cases) and agg_and_ribbon_active()
#' (n_active -> n_active). Sums `src_col` per (rep,
#' month, group), zero-fills the full rep x month x group grid, and returns
#' per (month, group) the importance-weighted quantiles q025..q975 and the
#' weighted MEAN (for skewed counts the median q50 sits below the mean, so the
#' mean is the honest central summary for calibration), named <out_col>_q*.
#' @noRd
.agg_and_ribbon_impl <- function(
  long_tbl,
  group_var,
  rep_weights,
  all_months_vec,
  src_col,
  out_col
) {
  agg <- long_tbl %>%
    dplyr::group_by(rep_id, t_start, .data[[group_var]]) %>%
    dplyr::summarise(
      !!out_col := sum(.data[[src_col]], na.rm = TRUE),
      .groups = "drop"
    )

  all_reps <- sort(unique(long_tbl$rep_id))
  all_groups <- sort(unique(long_tbl[[group_var]]))
  grid <- expand.grid(
    rep_id = all_reps,
    t_start = all_months_vec,
    grp = all_groups,
    stringsAsFactors = FALSE
  )
  grid$t_start <- as.Date(grid$t_start, origin = "1970-01-01")
  names(grid)[3] <- group_var

  agg_full <- merge(
    grid,
    agg,
    by = c("rep_id", "t_start", group_var),
    all.x = TRUE
  )
  agg_full[[out_col]][is.na(agg_full[[out_col]])] <- 0L
  agg_full <- dplyr::left_join(agg_full, rep_weights, by = "rep_id")

  q <- function(p) paste0(out_col, "_q", p)
  agg_full %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("t_start", group_var)))) %>%
    dplyr::summarise(
      !!q("025") := wquantile(.data[[out_col]], w, 0.025),
      !!q("10") := wquantile(.data[[out_col]], w, 0.10),
      !!q("25") := wquantile(.data[[out_col]], w, 0.25),
      !!q("50") := wquantile(.data[[out_col]], w, 0.50),
      !!q("75") := wquantile(.data[[out_col]], w, 0.75),
      !!q("90") := wquantile(.data[[out_col]], w, 0.90),
      !!q("975") := wquantile(.data[[out_col]], w, 0.975),
      !!paste0(out_col, "_mean") := stats::weighted.mean(
        .data[[out_col]],
        w,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

# Aggregate a sparse event table into weighted quantile ribbons over
# (t_start, group_var). `all_months_vec` is the month axis to fill zeros on.
# (Scratchpad Section 8.)
agg_and_ribbon <- function(long_tbl, group_var, rep_weights, all_months_vec) {
  .agg_and_ribbon_impl(
    long_tbl,
    group_var,
    rep_weights,
    all_months_vec,
    src_col = "Yout",
    out_col = "n_outbreaks"
  )
}


# Cases-only ribbon aggregator for the mon-track sparse table
# (smoothing_events_mon_long(): one row per (rep, gid, ACTIVE month) cell,
# column `cases_month`, no Yout — mon has no occurrence component of its own,
# it rides on the duration model's active/mso states). Same weighted-quantile
# summary as agg_and_ribbon()'s n_outbreaks_* block, sourced from
# `cases_month`, emitting the `n_cases_*` column stem consumed by
# ribbon_plot()/panel_pair().
agg_and_ribbon_mon <- function(
  long_tbl,
  group_var,
  rep_weights,
  all_months_vec
) {
  .agg_and_ribbon_impl(
    long_tbl,
    group_var,
    rep_weights,
    all_months_vec,
    src_col = "cases_month",
    out_col = "n_cases"
  )
}


# WHO dashboard comparison: national-scale smoothed monthly cases per country
# vs the WHO dashboard counts (08c_smoothing_compile.R, run comparisons).
# `long_mon` is the UNMASKED mon event table (every simulated gid counts
# toward the national total). The ribbon spans the full [start, end) month
# grid, zero-filled, and each band is divided by the covered-population share
# `who_f` (1 where unknown). `who_scale`: NULL = count / who_f; else
# who_multiplier(who_f, scale) applies (scalar or per-country) and the outputs
# gain a `who_scale` column.
#
# @return list(ribbon, obs = WHO complete-month counts, summary = per-country
#   months / totals / months inside, above or below the 95% band).
who_vs_smoothing <- function(
  long_mon,
  who_monthly,
  rep_weights,
  start,
  end,
  who_f = NULL,
  countries = NULL,
  who_scale = NULL
) {
  months <- seq(as.Date(start), as.Date(end) - 1, by = "month")
  obs <- who_obs_table(who_monthly, start, end)
  countries <- countries %||% sort(unique(long_mon$country))
  countries <- intersect(countries, unique(obs$country))
  obs <- obs[obs$country %in% countries, , drop = FALSE]
  if (!length(countries)) {
    return(NULL)
  }
  lm <- long_mon[
    long_mon$country %in% countries & long_mon$t_start %in% months,
    c("rep_id", "t_start", "country", "cases_month")
  ]
  # One zero row per (missing country, first rep, first month) makes the
  # aggregator's group axis cover every requested country; zeros change no sum.
  miss <- setdiff(countries, unique(lm$country))
  if (length(miss)) {
    lm <- dplyr::bind_rows(
      lm,
      data.frame(
        rep_id = rep_weights$rep_id[1],
        t_start = months[1],
        country = miss,
        cases_month = 0
      )
    )
  }
  rw <- rep_weights[
    rep_weights$rep_id %in% unique(long_mon$rep_id),
    ,
    drop = FALSE
  ]
  rib <- .agg_and_ribbon_impl(
    lm,
    "country",
    rw,
    months,
    src_col = "cases_month",
    out_col = "n_cases"
  )
  f <- if (is.null(who_f)) {
    rep(1, length(countries))
  } else {
    unname(who_f[countries]) %||% rep(1, length(countries))
  }
  f[is.na(f)] <- 1
  fc <- stats::setNames(f, countries)
  qcols <- grep("^n_cases_", names(rib), value = TRUE)
  if (is.null(who_scale)) {
    rib[qcols] <- lapply(rib[qcols], function(v) v / fc[rib$country])
  } else {
    sc <- if (length(who_scale) == 1L && is.null(names(who_scale))) {
      stats::setNames(rep(as.numeric(who_scale), length(countries)), countries)
    } else {
      stats::setNames(as.numeric(who_scale[countries]), countries)
    }
    sc[is.na(sc)] <- 1
    mult <- who_multiplier(fc, sc)
    rib[qcols] <- lapply(rib[qcols], function(v) v * mult[rib$country])
    rib$who_scale <- unname(sc[rib$country])
  }
  rib$who_f <- unname(fc[rib$country])
  cmp <- dplyr::inner_join(obs, rib, by = c("country", "t_start"))
  summary <- cmp %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      months = dplyr::n(),
      who_f = dplyr::first(who_f),
      who_total = sum(cases_who),
      sm_mean_total = sum(n_cases_mean),
      sm_median_total = sum(n_cases_q50),
      in95 = sum(cases_who >= n_cases_q025 & cases_who <= n_cases_q975),
      above95 = sum(cases_who > n_cases_q975),
      below95 = sum(cases_who < n_cases_q025),
      .groups = "drop"
    ) %>%
    dplyr::mutate(ratio_mean = sm_mean_total / pmax(who_total, 1))
  if (!is.null(who_scale)) {
    summary$who_scale <- rib$who_scale[match(summary$country, rib$country)]
  }
  list(ribbon = rib, obs = obs, summary = summary)
}

#' One facet per country: smoothed national-scale band(s) vs WHO counts.
#' `ribbons` is a named list of who_vs_smoothing()$ribbon tables (one per run;
#' the names label the runs), `obs` the WHO counts.
plot_who_vs_smoothing <- function(ribbons, obs, title = NULL, ncol = 4L) {
  rb <- dplyr::bind_rows(ribbons, .id = "run")
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = rb,
      ggplot2::aes(
        t_start,
        ymin = n_cases_q025,
        ymax = n_cases_q975,
        fill = run
      ),
      alpha = 0.18
    ) +
    ggplot2::geom_ribbon(
      data = rb,
      ggplot2::aes(t_start, ymin = n_cases_q25, ymax = n_cases_q75, fill = run),
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = rb,
      ggplot2::aes(t_start, n_cases_q50, colour = run)
    ) +
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(t_start, cases_who),
      size = 1.1
    ) +
    ggplot2::facet_wrap(~country, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      x = NULL,
      y = "Suspected cases per month (national scale)",
      fill = NULL,
      colour = NULL,
      title = title,
      subtitle = "Points: WHO dashboard + monthly updates (complete months). Bands: smoothing 50% / 95%, line: median; simulated gids scaled by 1 / covered population share."
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Ongoing-outbreak ("currently active episode") ribbon aggregator. Mirrors
# agg_and_ribbon_mon() exactly (same "presummed count per (rep, t_start,
# group)" input shape, same fill-zero-on-the-month-grid step) but sources
# `n_active` -- a per-(rep, country, month) COUNT already collapsed over units
# by the caller (08c_smoothing_compile.R, from .smoothing_active_matrix()) --
# instead of a per-cell indicator, and emits the `n_active_*` column stem so
# ribbon_plot()/panel_pair() work unchanged against this table too.
agg_and_ribbon_active <- function(
  long_tbl,
  group_var,
  rep_weights,
  all_months_vec
) {
  .agg_and_ribbon_impl(
    long_tbl,
    group_var,
    rep_weights,
    all_months_vec,
    src_col = "n_active",
    out_col = "n_active"
  )
}

# Observed monthly cases (mon track) per month x {Africa, AFRO_region,
# country}, the observed side of 08c's cases ribbons. Each (gid, month) keeps
# only its admin level with the largest total (ties to the finer), the
# filter's own footprint rule (spatpomp_core.R::.aggregate_footprint_counts()),
# so nested reports are not double-counted.
#
# @param mon_dat  03_prepare_mon.R's monthly frame.
# @param gid_meta gid -> AFRO_region / country map.
# @param start,end filter window, [start, end).
# @param months   months the simulation produced.
# @param mask     function(df, label) restricting rows to the coverage
#   footprint; identity by default.
# @return list(africa, region, country) of tibbles with column n_cases.
smoothing_obs_cases_mon <- function(
  mon_dat,
  gid_meta,
  start,
  end,
  months,
  mask = function(df, label) df
) {
  raw <- mon_dat %>%
    dplyr::mutate(
      gid = as.character(gid),
      t_start = as.Date(month),
      cases = as.numeric(cases_month)
    ) %>%
    dplyr::filter(
      !is.na(t_start),
      !is.na(cases),
      t_start >= start,
      t_start < end
    )
  raw <- finest_level_rows(raw, "cases", c("gid", "t_start")) %>% # footprint.R
    mask("obs_cases_mon") %>%
    dplyr::select(-dplyr::any_of(c("AFRO_region", "country"))) %>%
    dplyr::filter(t_start %in% months) %>%
    dplyr::left_join(gid_meta, by = "gid") %>%
    dplyr::filter(!is.na(AFRO_region))

  list(
    africa = raw %>%
      dplyr::group_by(t_start) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(africa = "Africa"),
    region = raw %>%
      dplyr::group_by(t_start, AFRO_region) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop"),
    country = raw %>%
      dplyr::group_by(t_start, country) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop")
  )
}


# Single ribbon panel (observed points over simulated quantile ribbons).
#
# `ylim` (optional, default NULL = unchanged for every pre-existing call site)
# applies a VIEW clamp via coord_cartesian(): ribbons and lines that exceed it
# run off the panel edge rather than being dropped, so no data is silently
# removed and ggplot emits no removed-rows warning. Use it when one outlier
# unit-month compresses the rest of the series into the axis floor.
# Masked series is always grey, unmasked series is always green -- a fixed
# two-colour scheme identifying coverage-mask status, independent of which
# sub-panel (occurrence/cases/mon-cases/active) is being drawn. `fill_col` is
# kept as a formal parameter (existing panel_pair() call sites still pass a
# per-panel colour) but is not used for ribbon/line colour.
.RIBBON_MASKED_COL <- "grey45"
.RIBBON_UNMASKED_COL <- "#2ca02c"

ribbon_plot <- function(
  ribbon,
  obs,
  group_col,
  group_val,
  y_col_stem,
  y_lab,
  title = NULL,
  fill_col = "#377eb8",
  ylim = NULL,
  ribbon_unmasked = NULL
) {
  d <- ribbon[ribbon[[group_col]] == group_val, , drop = FALSE]
  ob <- obs[obs[[group_col]] == group_val, , drop = FALSE]
  q025 <- paste0(y_col_stem, "_q025")
  q10 <- paste0(y_col_stem, "_q10")
  q25 <- paste0(y_col_stem, "_q25")
  q50 <- paste0(y_col_stem, "_q50")
  q75 <- paste0(y_col_stem, "_q75")
  q90 <- paste0(y_col_stem, "_q90")
  q975 <- paste0(y_col_stem, "_q975")
  y_obs <- y_col_stem
  mean_col <- paste0(y_col_stem, "_mean")

  du <- NULL
  if (!is.null(ribbon_unmasked)) {
    du <- ribbon_unmasked[
      ribbon_unmasked[[group_col]] == group_val,
      ,
      drop = FALSE
    ]
  }

  fill_values <- c(masked = .RIBBON_MASKED_COL)
  colour_values <- c(masked = .RIBBON_MASKED_COL)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = t_start)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[q025]], ymax = .data[[q975]], fill = "masked"),
      alpha = 0.12
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[q10]], ymax = .data[[q90]], fill = "masked"),
      alpha = 0.22
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[q25]], ymax = .data[[q75]], fill = "masked"),
      alpha = 0.38
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[[q50]], colour = "masked", linetype = "median"),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[[mean_col]], colour = "masked", linetype = "mean"),
      linewidth = 0.8
    )

  if (!is.null(du)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = du,
        ggplot2::aes(
          ymin = .data[[q025]],
          ymax = .data[[q975]],
          fill = "unmasked"
        ),
        alpha = 0.12
      ) +
      ggplot2::geom_ribbon(
        data = du,
        ggplot2::aes(
          ymin = .data[[q10]],
          ymax = .data[[q90]],
          fill = "unmasked"
        ),
        alpha = 0.22
      ) +
      ggplot2::geom_ribbon(
        data = du,
        ggplot2::aes(
          ymin = .data[[q25]],
          ymax = .data[[q75]],
          fill = "unmasked"
        ),
        alpha = 0.38
      ) +
      ggplot2::geom_line(
        data = du,
        ggplot2::aes(
          y = .data[[q50]],
          colour = "unmasked",
          linetype = "median"
        ),
        linewidth = 0.8
      ) +
      ggplot2::geom_line(
        data = du,
        ggplot2::aes(
          y = .data[[mean_col]],
          colour = "unmasked",
          linetype = "mean"
        ),
        linewidth = 0.8
      )
    fill_values <- c(fill_values, unmasked = .RIBBON_UNMASKED_COL)
    colour_values <- c(colour_values, unmasked = .RIBBON_UNMASKED_COL)
  }

  p <- p +
    ggplot2::geom_point(
      data = ob,
      ggplot2::aes(y = .data[[y_obs]]),
      colour = "black",
      size = 0.9,
      alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(name = NULL, values = fill_values) +
    ggplot2::scale_colour_manual(name = NULL, values = colour_values) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c(median = "solid", mean = "22"),
      breaks = c("median", "mean")
    ) +
    ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ggplot2::labs(
      x = NULL,
      y = y_lab,
      title = if (!is.null(title)) title else group_val
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  p
}

# Stacked occurrence (+ mon-track cases + ongoing-outbreak) panel for one
# group level. `ribbon_mon`/`obs_cases_mon` and `ribbon_active`/`obs_active`
# are optional (default NULL): supply both members of a pair to add that
# panel, omit either to skip it. This is how fig1/fig2/fig3
# (08c_smoothing_compile.R) fold the mon-track cases series and the
# ongoing-outbreak series into the main figures when a monthly/duration fit
# is configured.
#
# `ylim_cases` (optional, default NULL = unchanged) is forwarded to
# ribbon_plot()'s view clamp for the mon-track CASES panel only. The
# occurrence and ongoing-outbreak panels are counts of outbreaks, not cases,
# and are deliberately left on their own free scale.
panel_pair <- function(
  ribbon,
  obs_occ,
  group_col,
  group_val,
  ribbon_mon = NULL,
  obs_cases_mon = NULL,
  ribbon_active = NULL,
  obs_active = NULL,
  ylim_cases = NULL,
  ribbon_unmasked = NULL,
  ribbon_mon_unmasked = NULL,
  ribbon_active_unmasked = NULL
) {
  p_occ <- ribbon_plot(
    ribbon,
    obs_occ,
    group_col,
    group_val,
    "n_outbreaks",
    "Outbreaks / month",
    fill_col = "#377eb8",
    ribbon_unmasked = ribbon_unmasked
  )
  plots <- list(p_occ)
  if (!is.null(ribbon_mon) && !is.null(obs_cases_mon)) {
    plots[[length(plots) + 1L]] <- ribbon_plot(
      ribbon_mon,
      obs_cases_mon,
      group_col,
      group_val,
      "n_cases",
      "Suspected cases / month",
      title = paste0(group_val, " — cases (mon track)"),
      fill_col = "#4daf4a",
      ylim = ylim_cases,
      ribbon_unmasked = ribbon_mon_unmasked
    )
  }
  if (!is.null(ribbon_active) && !is.null(obs_active)) {
    plots[[length(plots) + 1L]] <- ribbon_plot(
      ribbon_active,
      obs_active,
      group_col,
      group_val,
      "n_active",
      "Ongoing outbreaks / month",
      title = paste0(group_val, " — ongoing"),
      fill_col = "#984ea3",
      ribbon_unmasked = ribbon_active_unmasked
    )
  }
  cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")
}
