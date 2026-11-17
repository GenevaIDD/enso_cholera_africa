#!/usr/bin/env Rscript
# make_hindcast_configs.R — per-(unit x launch) hindcast configs for the
#   forecast-skill test: launch from the final smoothing, observed ENSO
#
# Usage:
#   Rscript analysis/utils/deploy_prod_run.R -C <D.yml> --stage gen_hindcast
#   Rscript analysis/utils/make_hindcast_configs.R --template <forecast tpl> \
#     --run_tag_base <base> --config_outdir <dir> [--launch_month 05-01] \
#     [--season_end_by_group "enso_east=7,..."] [--launch_years 2012:2022]
#
# Details:
#   The in-sample counterpart of make_forecast_configs.R (same deployed-lag,
#   anchor and merged-dir machinery): each run hands off from an interior month
#   of the final smoothing cache, reads OBSERVED ENSO (perfect foresight) and
#   runs to the end of the following cholera season, so its output can be
#   scored against observations (07c_simulator_postprocess.R --hindcast).
#   Writes the configs and a manifest. patch_hindcast_config() is shared with
#   make_event_attribution_configs.R (its factual arm).
#
# See: the "Pass 1" section for the horizon, season-end and season-window rules

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(lubridate)
})

source("analysis/utils/config.R") # load_grouping
source("analysis/utils/base.R") # %||%
source("analysis/utils/paths.R") # make_run_file, read_deployed_lag()
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R") # resolve_base_rate_paths()
source("analysis/utils/season_calendar.R") # per-country season ends (--season_end)
source("analysis/utils/simulator_core.R") # sim_pop_year_spec() -- population year (manifest provenance)

# The forecast factual template, base run_tag and config_outdir have NO
# default: they name a prod_run<N>, which is recorded only in
# analysis/configs/deploy/prod_run<N>.yml (deploy.templates.forecast,
# run_tags.hindcast, outputs.hindcast_config_outdir). The standalone CLI below
# refuses to run without them.

option_list <- list(
  make_option(
    "--template",
    type = "character",
    default = NULL,
    help = "Forecast factual template to inherit smoothing_input/grouping/output from"
  ),
  make_option(
    "--fits_dir",
    type = "character",
    default = "generated_data/enso_best_lag_fits",
    help = "Deployed best-fit artifact dir (== grouping.fits_dir); lag read from sl<N> tag"
  ),
  make_option(
    "--grouping",
    type = "character",
    default = "enso_grouping",
    help = "Grouping registry name (analysis/groupings/<name>.csv)"
  ),
  make_option(
    "--climate_panel",
    type = "character",
    default = "generated_data/sim_climate_panel_base-1990_2019.rds",
    help = "HISTORICAL continuous sim panel (observed climate + observed enso_lag)"
  ),
  make_option(
    "--launch_month",
    type = "character",
    default = "05-01",
    help = "Boreal-spring launch as MM-DD (default 05-01). NOT hardcoded."
  ),
  make_option(
    "--launch_years",
    type = "character",
    default = "",
    help = "Colon range 'A:B' or comma list; empty = auto from smoothing coverage"
  ),
  make_option(
    "--season_end",
    type = "character",
    default = "by_group",
    help = paste(
      "Season ends: by_group (--season_end_by_group; default, unchanged output)",
      "| by_country_from_panel (each country's min_seas - 1 from --climate_panel)"
    )
  ),
  make_option(
    "--season_end_by_group",
    type = "character",
    default = "",
    help = "Per-group following-season-end month, 'grp=MM,...' (1-12); groups absent use --horizon_months. Read only under --season_end by_group"
  ),
  make_option(
    "--season_min_window_months",
    type = "integer",
    default = 6L,
    help = "by_country_from_panel: the season end is the first one at least this many months after the launch [default: %default]"
  ),
  make_option(
    "--season_window",
    action = "store_true",
    default = FALSE,
    help = paste(
      "One run per (unit, launch month, season year): each country's season from the",
      "month after its trough (floored at --launch_floor_month) to the next trough;",
      "replaces --launch_month / --season_end. Requires --partition_by registry"
    )
  ),
  make_option(
    "--launch_floor_month",
    type = "integer",
    default = 5L,
    help = "--season_window: earliest launch month of a season (ONI spring barrier) [default: %default]"
  ),
  make_option(
    "--trough_rule",
    type = "character",
    default = "onset_rate",
    help = "--season_window: onset_rate (min_seas fallback) | min_seas [default: %default]"
  ),
  make_option(
    "--trough_min_onsets",
    type = "integer",
    default = 30L,
    help = "--season_window onset_rate: onsets a country needs for its own trough [default: %default]"
  ),
  make_option(
    "--trough_min_stability",
    type = "double",
    default = 0.8,
    help = "--season_window onset_rate: bootstrap share within +-1 month needed [default: %default]"
  ),
  make_option(
    "--trough_years",
    type = "character",
    default = "",
    help = "--season_window: calendar years of the onset panel read, 'A:B' or comma list ('' = all)"
  ),
  make_option(
    "--trough_smooth_months",
    type = "integer",
    default = 3L,
    help = "--season_window: circular smoothing of the onset rate, odd months [default: %default]"
  ),
  make_option(
    "--trough_bootstrap",
    type = "integer",
    default = 500L,
    help = "--season_window: bootstrap resamples of years [default: %default]"
  ),
  make_option(
    "--trough_seed",
    type = "integer",
    default = 1L,
    help = "--season_window: bootstrap seed (per country) [default: %default]"
  ),
  make_option(
    "--survival_panel",
    type = "character",
    default = "",
    help = "--season_window: occurrence panel for the onset rate ('' = the template episodes' panel)"
  ),
  make_option(
    "--horizon_months",
    type = "integer",
    default = 12L,
    help = "Fallback forecast length in months from launch when a group has no season-end entry"
  ),
  make_option(
    "--add_best_lag",
    action = "store_true",
    default = TRUE,
    help = "Extend horizon by min(best_lag_occ, best_lag_dur, best_lag_mon) so ENSO covariates stay in-window"
  ),
  make_option(
    "--no_add_best_lag",
    action = "store_true",
    default = FALSE,
    help = "Disable the best-lag horizon extension"
  ),
  make_option(
    "--smoothing_input",
    type = "character",
    default = "",
    help = "Override smoothing cache dir (default: inherit from template)"
  ),
  make_option(
    "--episode_key",
    type = "character",
    default = "",
    help = paste(
      "Episode definition written into every generated config",
      "(lp_season | gid_season). '' (default) inherits the",
      "template untouched; deploy_prod_run.R passes",
      "deploy.episodes.key. See apply_episode_key_to_sim_config()."
    )
  ),
  make_option(
    "--partition_by",
    type = "character",
    default = "",
    help = paste(
      "Forward-run top-level partition: '' (default, inherit",
      "template; country-partitioned) or 'enso_region' to run",
      "each ENSO group as one region TU (ADM1 contagion couples",
      "across borders; handoff IC assembled from the per-country",
      "smoothing caches). scope.countries/smoothing_input are",
      "unchanged — only filter.partition_by is written."
    )
  ),
  make_option(
    "--partition_source",
    type = "character",
    default = "",
    help = "Partition registry for --partition_by registry (one config per unit x year)"
  ),
  make_option(
    "--nbr_partition_source",
    type = "character",
    default = NULL,
    help = "Neighbour-contagion partition registry [default: --partition_source]"
  ),
  make_option(
    "--run_tag_base",
    type = "character",
    default = NULL,
    help = "Base run_tag (REQUIRED); per-config tag is <base>_<group>_<year>"
  ),
  make_option(
    "--config_outdir",
    type = "character",
    default = NULL,
    help = "Directory to write generated per-(group,year) configs (REQUIRED)"
  ),
  make_option(
    "--manifest",
    type = "character",
    default = "analysis/configs/hindcast_horizons.yml",
    help = "Path to write the hindcast horizons/merged-dir manifest"
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "generated_data/",
    help = "Simulator output root (for computing per-run merged dirs)"
  ),
  make_option(
    "--param_draw_mode",
    type = "character",
    default = "",
    help = "Posterior-draw mode of the hindcasts: '' (patch_hindcast_config's random) | random | joint"
  ),
  make_option(
    "--n_replicates",
    type = "integer",
    default = 0L,
    help = "Replicates per hindcast run; 0 = the template's simulator.n_replicates [default: %default]"
  )
)
# The three seams that make a run a hindcast:
#   1. filter.handoff_month  -> pins the smoothing handoff to <year>-<launch>
#      (an INTERIOR month of the single final smoothing cache) instead of the
#      end of coverage. 07b_simulator_spatpomp.R slices filter_traj at that
#      month (see load_smoothing_distribution() / filter_traj_to_warmup()).
#   2. enso.forecast_source  -> "observed": leave the panel's observed enso_lag
#      columns unmodified (no CFSv2/C3S override). See utils/simulations.R.
#   3. climate_panel         -> the HISTORICAL continuous panel (observed
#      climate + observed enso_lag), NOT the projection panel (zero-filled
#      future climate) the forecast template ships with.
# The smoothing_input, grouping, measurement, delta and output blocks are
# inherited VERBATIM from the forecast template, so the same final smoothing
# cache is reused for every (group, year) launch.
#
# patch_hindcast_config(): one (unit, launch) hindcast config off the inherited
# forecast template -- the three hindcast seams (interior handoff, observed ENSO,
# historical climate panel), the partition and episode keys and the anchors.
# Top level so the
# event attribution (make_event_attribution_configs.R) builds its factual arm
# with the same code; `opt` carries climate_panel, partition_by,
# partition_source, nbr_partition_source and episode_key.
patch_hindcast_config <- function(
  cfg,
  opt,
  run_tag,
  run_id,
  isos,
  launch_chr,
  horizon_chr,
  base_rate_table_paths
) {
  cfg$simulator$run_tag <- run_tag
  cfg$simulator$run_id <- run_id
  cfg$simulator$scope$countries <- paste(isos, collapse = ",")

  # Historical panel (observed climate + observed enso_lag), NOT the projection panel.
  cfg$simulator$climate_panel <- opt$climate_panel

  # Window: launch .. horizon_end (both first-of-month). horizon.end must equal
  # filter.forecast_window.end (07b invariant).
  cfg$simulator$horizon$start <- launch_chr
  cfg$simulator$horizon$end <- horizon_chr
  cfg$simulator$filter$forecast_window$start <- launch_chr
  cfg$simulator$filter$forecast_window$end <- horizon_chr

  # Seam 1: interior-month handoff (07b load_smoothing_distribution override).
  cfg$simulator$filter$handoff_month <- launch_chr

  # Optional region partition. Only written when --partition_by is given, so the
  # default (country-partitioned) output stays byte-for-byte identical. Under
  # enso_region the forward run couples member-country ADM1 units across borders
  # via the contagion kernel; scope.countries (the region member list) and
  # smoothing_input (per-country caches) are left untouched — the region handoff
  # IC is assembled from those same per-country caches by 07b.
  # %||% "" keeps the library entry point NULL-safe for any caller whose opt
  # carries no partition_by key (deploy_prod_run.R's .hindcast_opt() always
  # sets one): nzchar(NULL) is logical(0) and `if (logical(0))` errors.
  if (nzchar(opt$partition_by %||% "")) {
    cfg$simulator$filter$partition_by <- opt$partition_by
  }
  if (identical(opt$partition_by %||% "", "registry")) {
    cfg$simulator$filter$partition_source <- opt$partition_source
    cfg$simulator$nbr_contagion$partition_source <-
      opt$nbr_partition_source %||% opt$partition_source
  }

  # Episode definition: the deploy's key, written into the config and checked
  # against the prep artefacts it will read (the historical climate panel set
  # above, the size frame -- the location frame's population / admin_level
  # source -- the survival panel). See make_forecast_configs.R.
  cfg <- apply_episode_key_to_sim_config(
    cfg,
    opt$episode_key %||% "",
    label = run_tag
  )
  preflight_sim_episode_inputs(cfg, label = run_tag)

  # Seam 3: draw a FRESH random posterior-parameter set across the forecast
  # ensemble. The smoothing cache these hindcasts hand off from was built with
  # n_draws_in_filter=1, so every particle carries draw_idx=1 and the ensemble
  # would otherwise collapse to a single posterior draw, flattening the mon
  # tail.
  # "random" also means a region's member countries share the same per-replicate
  # draw (see .forecast_param_indices()); it is the global default, set here
  # explicitly for clarity/reproducibility of generated hindcast configs.
  cfg$simulator$forecast$param_draw_mode <- "random"

  # Seam 2: observed ENSO (perfect foresight); drop the CFSv2/C3S forecast keys so
  # nothing tries to fetch a future ENSO product for a historical window.
  cfg$simulator$enso$forecast_source <- "observed"
  cfg$simulator$enso$c3s_plume <- NULL
  cfg$simulator$enso$cfsv2_sd_correction <- NULL
  cfg$simulator$enso$cfsv2_target_peak <- NULL
  cfg$simulator$enso$cfsv2_calibration <- NULL
  cfg$simulator$enso$cfsv2_calibration_target <- NULL

  if (!is.null(base_rate_table_paths) && length(base_rate_table_paths)) {
    cfg$simulator$location_frame$base_rate_table_paths <- base_rate_table_paths
  }
  cfg
}

# hindcast_param_draw_mode() / hindcast_n_replicates(): the hindcasts' posterior-draw
# mode and replicate count,
# applied after patch_hindcast_config() (which writes "random" and keeps the
# template's count) by generate_hindcast_configs() only, so the event
# attribution, which reuses patch_hindcast_config(), is unchanged. "" / 0 leave
# the patched config as it was. "joint": replicate r starts
# from smoothing replicate ((r - 1) mod n_smoothing) + 1 and runs its draw (07b
# checks the cache is shared_cycle), so the handoff must be uniform and
# max_draws equal the smoothing's replicate count.
hindcast_param_draw_mode <- function(cfg, mode, n_smoothing = NULL) {
  mode <- as.character(mode %||% "")
  if (!nzchar(mode)) {
    return(cfg)
  }
  if (!mode %in% c("random", "joint")) {
    stop(
      "[make_hindcast_configs] --param_draw_mode must be random or joint; got '",
      mode,
      "'.",
      call. = FALSE
    )
  }
  if (identical(mode, "joint")) {
    handoff <- as.character(
      cfg$simulator$forecast[["handoff_weighting"]] %||% "uniform"
    )
    if (!identical(handoff, "uniform")) {
      stop(
        "[make_hindcast_configs] param_draw_mode joint pairs each replicate with its smoothing replicate: ",
        "the template's simulator.forecast.handoff_weighting must be uniform, not '",
        handoff,
        "'.",
        call. = FALSE
      )
    }
    S <- as.integer(cfg$simulator$max_draws %||% 200L)
    n <- as.integer(n_smoothing %||% NA)
    if (is.na(n) || n < 1L) {
      stop(
        "[make_hindcast_configs] param_draw_mode joint needs the smoothing cache's replicate count ",
        "(meta.yml nrep) to check max_draws against it.",
        call. = FALSE
      )
    }
    if (S != n) {
      stop(
        sprintf(
          paste0(
            "[make_hindcast_configs] param_draw_mode joint: max_draws %d but %d smoothing ",
            "replicates; the joint design runs one draw per smoothing replicate."
          ),
          S,
          n
        ),
        call. = FALSE
      )
    }
  }
  cfg$simulator$forecast$param_draw_mode <- mode
  cfg
}

hindcast_n_replicates <- function(cfg, n) {
  n <- as.integer(n %||% 0L)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop(
      "[make_hindcast_configs] --n_replicates must be an integer >= 0; got '",
      n,
      "'.",
      call. = FALSE
    )
  }
  if (n > 0L) {
    cfg$simulator$n_replicates <- n
  }
  cfg
}

# ---------------------------------------------------------------------------
# generate_hindcast_configs(opt) — library entry point
# ---------------------------------------------------------------------------
# Body wrapped in a function so the deploy orchestrator
# (analysis/utils/deploy_prod_run.R) can call it IN-PROCESS with an `opt` built from
# the deploy block, instead of shelling out to Rscript. The thin CLI wrapper at
# the bottom (guarded by sys.nframe() == 0L) preserves the standalone
# `Rscript analysis/utils/make_hindcast_configs.R ...` entrypoint byte-for-byte.
# The per-run helpers below close over `opt` exactly as before.
generate_hindcast_configs <- function(opt) {
  add_best_lag <- isTRUE(opt$add_best_lag) && !isTRUE(opt$no_add_best_lag)

  # ---------------------------------------------------------------------------
  # Helpers (mirrors of make_forecast_configs.R; the deployed lag of each track
  # comes from utils/paths.R::read_deployed_lag())
  # ---------------------------------------------------------------------------

  merged_dir_for <- function(run_tag, run_id, fmt = "qs") {
    o <- list(
      outdir = opt$outdir,
      simulator = list(
        run_tag = run_tag,
        run_id = run_id,
        output = list(format = fmt)
      )
    )
    make_run_file(o, kind = "spatpomp_merged_dir")
  }

  # Configs are written with provenance.R::write_yaml_with_header().
  write_config <- write_yaml_with_header

  # Smoothing cache coverage from filenames: utils/paths.R::smoothing_cache_coverage()
  # (the intersection across TUs, so every TU has each candidate IC).
  smoothing_coverage <- smoothing_cache_coverage

  parse_kv_int <- function(s) {
    s <- trimws(s)
    if (!nzchar(s)) {
      return(setNames(integer(0), character(0)))
    }
    parts <- trimws(strsplit(s, ",")[[1]])
    parts <- parts[nzchar(parts)]
    out <- integer(0)
    for (p in parts) {
      kv <- strsplit(p, "=")[[1]]
      if (length(kv) != 2L) {
        stop("Bad key=value in '--season_end_by_group': ", p)
      }
      out[trimws(kv[1])] <- as.integer(trimws(kv[2]))
    }
    out
  }

  # ---------------------------------------------------------------------------
  # Resolve inputs
  # ---------------------------------------------------------------------------

  reg <- load_grouping(opt$grouping)
  groups <- sort(unique(reg$group))

  tpl <- yaml::read_yaml(opt$template)
  # Episode definition the deployed fits are looked up under (the "_ep-gs"
  # filename filter of read_deployed_lag()): the deploy's key when given, else
  # the template's own episodes.key.
  ep_key <- if (nzchar(opt$episode_key %||% "")) {
    opt$episode_key
  } else {
    (tpl$episodes %||% list())$key %||% "lp_season"
  }
  # national zero-fill occurrence fits ("_nz"): the template's own declaration
  # (gen_templates writes deploy.episodes.national_zero_fill into it; the deploy
  # re-checks it in deploy_check_sim_template()).
  ep_nz <- isTRUE((tpl$episodes %||% list())$national_zero_fill)
  ep_nzt <- .episodes_nz_threshold(tpl$episodes)
  # The whole occurrence-panel spec (utils/paths.R; panel variant tokens).
  ep_spec <- .occ_panel_spec_from_episodes(utils::modifyList(
    tpl$episodes %||% list(),
    list(key = ep_key)
  ))
  smoothing_input <- if (nzchar(opt$smoothing_input)) {
    opt$smoothing_input
  } else {
    tpl$simulator$filter$smoothing_input
  }
  if (is.null(smoothing_input) || !nzchar(smoothing_input)) {
    stop("No smoothing_input: not in template and --smoothing_input not given.")
  }

  season_end_map <- parse_kv_int(opt$season_end_by_group)
  bad_groups <- setdiff(names(season_end_map), groups)
  if (length(bad_groups)) {
    stop(
      "--season_end_by_group names unknown group(s): ",
      paste(bad_groups, collapse = ", ")
    )
  }

  # Season-end source (utils/season_calendar.R). by_group keeps the group map
  # above (output unchanged) and reads the panel's calendar only for the
  # hemisphere-mix diagnostic; by_country_from_panel derives each country's end.
  season_mode <- as.character(opt[["season_end"]] %||% "by_group")
  if (!season_mode %in% c("by_group", "by_country_from_panel")) {
    stop(
      "[make_hindcast_configs] --season_end must be by_group or by_country_from_panel; got '",
      season_mode,
      "'.",
      call. = FALSE
    )
  }
  by_country <- identical(season_mode, "by_country_from_panel")
  min_win <- as.integer(opt[["season_min_window_months"]] %||% 6L)
  if (by_country && (length(min_win) != 1L || is.na(min_win) || min_win < 1L)) {
    stop(
      "[make_hindcast_configs] --season_min_window_months must be >= 1.",
      call. = FALSE
    )
  }
  se_cal <- NULL
  se_tbl <- NULL
  if (by_country) {
    se_cal <- read_panel_season_calendar(opt$climate_panel)
    se_tbl <- season_end_by_country(se_cal)
  } else if (
    nzchar(opt$climate_panel %||% "") && file.exists(opt$climate_panel)
  ) {
    se_tbl <- tryCatch(
      season_end_by_country(read_panel_season_calendar(opt$climate_panel)),
      error = function(e) {
        message(
          "[make_hindcast_configs] season-end diagnostic skipped: ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }

  parse_years <- function(s) {
    s <- trimws(s %||% "")
    if (!nzchar(s)) {
      return(NULL)
    }
    y <- if (grepl(":", s)) {
      ab <- as.integer(strsplit(s, ":")[[1]])
      seq(ab[1], ab[2])
    } else {
      as.integer(trimws(strsplit(s, ",")[[1]]))
    }
    sort(unique(y[!is.na(y)]))
  }

  # --season_window: per-country troughs (onset rate where supported, else min_seas).
  season_win <- isTRUE(opt[["season_window"]])
  sw_tbl <- NULL
  sw_src <- NULL
  floor_m <- NULL
  if (season_win) {
    if (!identical(opt$partition_by %||% "", "registry")) {
      stop(
        "[make_hindcast_configs] --season_window needs --partition_by registry.",
        call. = FALSE
      )
    }
    floor_m <- as.integer(opt[["launch_floor_month"]] %||% 5L)
    t_rule <- as.character(opt[["trough_rule"]] %||% "onset_rate")
    if (!t_rule %in% c("onset_rate", "min_seas")) {
      stop(
        "[make_hindcast_configs] --trough_rule must be onset_rate or min_seas; got '",
        t_rule,
        "'.",
        call. = FALSE
      )
    }
    se_cal <- read_panel_season_calendar(opt$climate_panel)
    t_years <- parse_years(opt[["trough_years"]])
    onsets <- NULL
    surv_path <- NA_character_
    if (identical(t_rule, "onset_rate")) {
      surv_path <- as.character(opt[["survival_panel"]] %||% "")
      if (!nzchar(surv_path)) {
        cfg0 <- apply_episode_key_to_sim_config(
          tpl,
          opt$episode_key %||% "",
          label = "season_window"
        )
        surv_path <- as.character(make_surv_dat_file(list(
          baseline = cfg0$baseline %||% "1990_2019",
          episodes = cfg0$episodes
        )))
      }
      if (!file.exists(surv_path)) {
        stop(
          "[make_hindcast_configs] --season_window: occurrence panel not found: ",
          surv_path,
          call. = FALSE
        )
      }
      sv <- readRDS(surv_path)
      if (!is.data.frame(sv)) {
        sv <- sv[[1L]]
      }
      onsets <- onset_rate_troughs(
        sv,
        years = t_years,
        smooth_months = as.integer(opt[["trough_smooth_months"]] %||% 3L),
        bootstrap = as.integer(opt[["trough_bootstrap"]] %||% 500L),
        seed = as.integer(opt[["trough_seed"]] %||% 1L)
      )
      rm(sv)
    }
    by_country <- FALSE # the season windows replace --season_end
    sw_tbl <- season_trough_by_country(
      season_end_by_country(se_cal),
      onsets,
      t_rule,
      min_onsets = as.integer(opt[["trough_min_onsets"]] %||% 30L),
      min_stability = as.numeric(opt[["trough_min_stability"]] %||% 0.8)
    )
    onm <- attr(sw_tbl, "onsets_by_month")
    sw_src <- list(
      rule = paste(
        "season = month after the trough .. next trough; launch = season start floored at",
        "launch_floor_month; firm horizon = season end + max(0, min lag); trough = onset-rate",
        "minimum where n_onsets >= min_onsets and stability >= min_stability, else modal min_seas"
      ),
      launch_floor_month = floor_m,
      trough_rule = t_rule,
      min_onsets = as.integer(opt[["trough_min_onsets"]] %||% 30L),
      min_stability = as.numeric(opt[["trough_min_stability"]] %||% 0.8),
      trough_years = if (is.null(t_years)) {
        ""
      } else {
        sprintf("%d:%d", min(t_years), max(t_years))
      },
      smooth_months = as.integer(opt[["trough_smooth_months"]] %||% 3L),
      bootstrap = as.integer(opt[["trough_bootstrap"]] %||% 500L),
      seed = as.integer(opt[["trough_seed"]] %||% 1L),
      survival_panel = surv_path,
      onsets_digest = if (
        !is.null(onm) && requireNamespace("digest", quietly = TRUE)
      ) {
        digest::digest(onm, algo = "xxhash64")
      } else {
        ""
      },
      climate_panel = as.character(attr(se_cal, "path")),
      min_seas_digest = as.character(attr(se_cal, "digest"))
    )
    for (i in seq_len(nrow(sw_tbl))) {
      message(sprintf(
        "[make_hindcast_configs] trough %s: %s (%s; min_seas %d, onsets %s, stability %s)",
        sw_tbl$iso3[i],
        month.abb[sw_tbl$trough_month[i]],
        sw_tbl$trough_source[i],
        sw_tbl$min_seas_mode[i],
        format(sw_tbl$n_onsets[i]),
        format(round(sw_tbl$stability[i], 2))
      ))
    }
  }

  launch_mmdd <- trimws(opt$launch_month)
  if (!grepl("^[0-9]{2}-[0-9]{2}$", launch_mmdd)) {
    stop("--launch_month must be MM-DD (e.g. 05-01), got: ", launch_mmdd)
  }

  # Launch years: explicit (--launch_years) or, by default, those whose IC month
  # (launch - 1 month) lies inside the smoothing cache coverage, parsed from the
  # cache filenames' date tokens.
  cov <- smoothing_coverage(smoothing_input)
  if (nzchar(opt$launch_years)) {
    ly <- trimws(opt$launch_years)
    years <- if (grepl(":", ly)) {
      ab <- as.integer(strsplit(ly, ":")[[1]])
      seq(ab[1], ab[2])
    } else {
      as.integer(trimws(strsplit(ly, ",")[[1]]))
    }
    years <- sort(unique(years[!is.na(years)]))
    if (!length(years)) stop("--launch_years resolved to nothing.")
  } else {
    if (is.null(cov)) {
      stop(
        "Cannot auto-derive launch years: no hist_filter_tu=*.qs in ",
        smoothing_input,
        ". Pass --launch_years explicitly."
      )
    }
    # IC month = launch - 1 must be in [cov$start, cov$end - 1 month].
    y0 <- year(cov$start %m+% months(1L))
    y1 <- year(cov$end)
    cand <- seq(y0, y1)
    ic_ok <- vapply(
      cand,
      function(y) {
        ic <- as.Date(sprintf("%d-%s", y, launch_mmdd)) %m-% months(1L)
        ic >= cov$start && ic <= (cov$end %m-% months(1L))
      },
      logical(1)
    )
    years <- cand[ic_ok]
    if (!length(years)) {
      stop(
        "No launch year has its IC month within smoothing coverage [",
        format(cov$start),
        ", ",
        format(cov$end),
        ")."
      )
    }
  }
  message(sprintf(
    "[make_hindcast_configs] launch %s; years: %s%s",
    launch_mmdd,
    paste(years, collapse = ", "),
    if (!is.null(cov)) {
      sprintf("  (coverage [%s, %s))", format(cov$start), format(cov$end))
    } else {
      ""
    }
  ))

  date_str <- format(Sys.Date(), "%Y%m%d")
  out_fmt <- tpl$simulator$output$format %||% "qs"
  dir.create(opt$config_outdir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Patch one (group, year) config off the inherited template
  # ---------------------------------------------------------------------------

  # The hindcasts' draw mode and replicate count (--param_draw_mode /
  # --n_replicates; hindcast_param_draw_mode()); the smoothing's replicate count
  # is read only when joint needs it.
  hc_n_smoothing <- if (identical(opt$param_draw_mode %||% "", "joint")) {
    mp <- file.path(smoothing_input, "meta.yml")
    if (file.exists(mp)) {
      as.integer(yaml::read_yaml(mp)$nrep %||% NA)
    } else {
      NA_integer_
    }
  }
  patch_hindcast <- function(
    cfg,
    run_tag,
    run_id,
    isos,
    launch_chr,
    horizon_chr,
    base_rate_table_paths
  ) {
    cfg <- patch_hindcast_config(
      cfg,
      opt,
      run_tag,
      run_id,
      isos,
      launch_chr,
      horizon_chr,
      base_rate_table_paths
    )
    cfg <- hindcast_param_draw_mode(cfg, opt$param_draw_mode, hc_n_smoothing)
    hindcast_n_replicates(cfg, opt$n_replicates)
  }

  manifest_groups <- list()
  skipped <- c()
  n_written <- 0L

  # ---- Pass 1: per FIT GROUP -- deployed lags, anchors, per-year firm windows --
  # Horizon per (group, year):
  #   launch      = <year>-<launch_month>
  #   season_end  = first day of --season_end_by_group[group] month strictly after
  #                 launch (rolled +1yr if <= launch); else launch + --horizon_months
  #   horizon_end = season_end + max(0, min(best_lag_occ, best_lag_dur,
  #                 best_lag_mon)) when --add_best_lag.
  #                 A deployed ENSO lead (negative min lag) does not shorten the
  #                 window: the hindcast reads OBSERVED ENSO (perfect foresight),
  #                 so no lag limits which months can be simulated. (The forecast
  #                 keeps last CFSv2 month + min lag, where a lead does limit it.)
  #
  # Season ends (--season_end):
  #   by_group (default)     the --season_end_by_group month of the country's fit
  #                          group, as above; output unchanged. Warns when a
  #                          group's member countries' cholera-year ends (from the
  #                          panel) span more than 4 months.
  #   by_country_from_panel  each country's own cholera-year end, min_seas - 1 of
  #                          its modal gid calendar in --climate_panel
  #                          (utils/season_calendar.R::season_end_by_country());
  #                          season_end = the first such month at least
  #                          --season_min_window_months (default 6) after the
  #                          launch; firm horizon = that + max(0, min lag) per
  #                          country. A group's horizon_end is the max over its
  #                          members and a unit runs to the max over its member
  #                          countries. The manifest gains countries.<iso>
  #                          (season_end_month, min_seas_mode, modal_share,
  #                          years.<y>.firm_horizon_end / unit_horizon_end) and
  #                          season_end_source (panel path + min_seas digest), and
  #                          season_ends_<manifest stem>.csv is written next to it.
  #
  # Season windows (--season_window; replaces --launch_month and --season_end;
  # registry units only):
  #   trough       each country's onset-rate trough in its occurrence panel
  #                (--trough_rule onset_rate, where it has >= --trough_min_onsets
  #                onsets and a trough stable under a bootstrap over years,
  #                >= --trough_min_stability), else its modal min_seas of
  #                --climate_panel (season_calendar.R::season_trough_by_country())
  #   season       the month after the trough to the next trough; launch = the
  #                season start, floored at --launch_floor_month (May: after the
  #                ONI spring predictability barrier); firm horizon = season end
  #                + max(0, min lag) (season_calendar.R::season_window_of())
  #   runs         one per (unit, launch month, season year), keyed "<year>-<MM>":
  #                the unit is simulated whole, only the countries whose season
  #                starts that month are scored in it (units.<u>.years.<key>.scored;
  #                countries.<iso>.years.<key> only under that key)
  #   reporting    countries.<iso> (trough source, floor cut and its onset share),
  #                manifest season_window (rule + inputs) and the CSV.
  group_info <- list()
  for (g in groups) {
    bl <- lapply(stats::setNames(nm = c("occ", "dur", "mon")), function(tr) {
      read_deployed_lag(
        opt$fits_dir,
        g,
        opt$grouping,
        tr,
        episode_key = ep_key,
        national_zero_fill = ep_nz,
        national_zero_threshold = ep_nzt,
        occ_panel = ep_spec
      )
    })
    missing_tr <- names(bl)[vapply(bl, is.null, logical(1))]
    if (length(missing_tr)) {
      message(sprintf(
        "[make_hindcast_configs] SKIP group '%s' (no deployed %s fit)",
        g,
        paste(missing_tr, collapse = "/")
      ))
      skipped <- c(skipped, g)
      next
    }
    lags <- vapply(bl, `[[`, integer(1), "lag")
    min_lag <- min(lags)
    isos <- reg$iso3[reg$group == g]
    if (min_lag < 0L) {
      message(sprintf(
        paste0(
          "[make_hindcast_configs] group '%s': deployed ENSO lead (lags %s) -- the hindcast ",
          "reads observed ONI up to %d month(s) after each forecast month (perfect foresight)."
        ),
        g,
        paste(lags, collapse = "/"),
        -min_lag
      ))
    }

    base_rate_table_paths <- resolve_base_rate_paths(
      lapply(bl, function(b) attr(readRDS(b$file), "opt"))
    )

    # --season_window: each member country's season per season year
    # (season_calendar.R::season_window_of()); the group's season-year entry
    # spans its earliest launch to its latest firm horizon.
    if (season_win) {
      miss_c <- setdiff(isos, sw_tbl$iso3)
      if (length(miss_c)) {
        stop(
          sprintf(
            paste0(
              "[make_hindcast_configs] group '%s': no gid of %s in the climate panel %s, ",
              "so no trough can be derived (--season_window)."
            ),
            g,
            paste(miss_c, collapse = ", "),
            opt$climate_panel
          ),
          call. = FALSE
        )
      }
      tr_c <- stats::setNames(
        sw_tbl$trough_month[match(isos, sw_tbl$iso3)],
        isos
      )
      win <- lapply(stats::setNames(nm = as.character(years)), function(yc) {
        cw <- lapply(stats::setNames(nm = isos), function(iso) {
          season_window_of(
            tr_c[[iso]],
            as.integer(yc),
            floor_m,
            min_lag,
            add_best_lag
          )
        })
        list(
          country = cw,
          launch_chr = min(vapply(
            cw,
            function(w) format(w$launch),
            character(1)
          )),
          horizon_chr = max(vapply(
            cw,
            function(w) format(w$horizon_end),
            character(1)
          ))
        )
      })
      group_info[[g]] <- list(
        bl = bl,
        lags = lags,
        min_lag = min_lag,
        isos = isos,
        base_rate_table_paths = base_rate_table_paths,
        win = win
      )
      next
    }

    # by_country_from_panel: each member country's own window per launch year
    # (season_calendar.R::hindcast_window_end()); the group's horizon is the max
    # over its members, the country horizons ride along for the manifest.
    if (by_country) {
      miss_c <- setdiff(isos, se_tbl$iso3)
      if (length(miss_c)) {
        stop(
          sprintf(
            paste0(
              "[make_hindcast_configs] group '%s': no gid of %s in the climate panel %s, ",
              "so no season end can be derived (--season_end by_country_from_panel)."
            ),
            g,
            paste(miss_c, collapse = ", "),
            opt$climate_panel
          ),
          call. = FALSE
        )
      }
      sem_c <- stats::setNames(
        se_tbl$season_end_month[match(isos, se_tbl$iso3)],
        isos
      )
      win <- lapply(stats::setNames(nm = as.character(years)), function(yc) {
        launch <- as.Date(sprintf("%d-%s", as.integer(yc), launch_mmdd))
        cw <- lapply(isos, function(iso) {
          hindcast_window_end(
            launch,
            sem_c[[iso]],
            min_lag,
            min_win,
            add_best_lag
          )
        })
        ch <- stats::setNames(
          vapply(cw, function(w) format(w$horizon_end), character(1)),
          isos
        )
        cs <- stats::setNames(
          vapply(cw, function(w) format(w$season_end), character(1)),
          isos
        )
        bad <- isos[as.Date(ch) <= launch]
        if (length(bad)) {
          stop(
            sprintf(
              "[make_hindcast_configs] group '%s', launch %s: firm window of %s is not after the launch.",
              g,
              format(launch),
              paste(bad, collapse = ", ")
            ),
            call. = FALSE
          )
        }
        list(
          launch_chr = format(launch),
          horizon_chr = max(ch),
          country_h = ch,
          country_se = cs
        )
      })
      group_info[[g]] <- list(
        bl = bl,
        lags = lags,
        min_lag = min_lag,
        isos = isos,
        base_rate_table_paths = base_rate_table_paths,
        win = win
      )
      next
    }

    # Firm (lag-defined) window per launch year: launch .. season end (or launch +
    # horizon_months) + min_lag.
    win <- lapply(stats::setNames(nm = as.character(years)), function(yc) {
      y <- as.integer(yc)
      launch <- as.Date(sprintf("%d-%s", y, launch_mmdd))
      if (g %in% names(season_end_map)) {
        sem <- season_end_map[[g]]
        end0 <- as.Date(sprintf("%04d-%02d-01", y, sem))
        if (end0 <= launch) {
          end0 <- end0 %m+% years(1L)
        }
        base_end <- end0
      } else {
        base_end <- launch %m+% months(opt$horizon_months)
      }
      # A deployed ENSO lead (negative min_lag) does NOT shorten the window: the
      # hindcast reads OBSERVED ONI from after the launch -- in the forecast window
      # and in the smoothing initial condition (launch - k .. launch - 1) -- the
      # hindcast's perfect-foresight ENSO, so no lag limits the window (adding a
      # lead would cut a May-Nov season with min lag -6 to nothing). Positive lags
      # still extend it.
      horizon_end <- if (add_best_lag) {
        base_end %m+% months(max(min_lag, 0L))
      } else {
        base_end
      }
      if (horizon_end <= launch) {
        stop(
          sprintf(
            paste0(
              "[make_hindcast_configs] group '%s', launch %s: firm window end %s ",
              "(season end %s + min lag %d) is not after the launch."
            ),
            g,
            format(launch),
            format(horizon_end),
            format(base_end),
            min_lag
          ),
          call. = FALSE
        )
      }
      list(launch_chr = format(launch), horizon_chr = format(horizon_end))
    })
    group_info[[g]] <- list(
      bl = bl,
      lags = lags,
      min_lag = min_lag,
      isos = isos,
      base_rate_table_paths = base_rate_table_paths,
      win = win
    )
  }
  if (!by_country && !season_win && !is.null(se_tbl)) {
    for (w in season_end_spread_warnings(
      lapply(group_info, `[[`, "isos"),
      se_tbl
    )) {
      warning("[make_hindcast_configs] ", w, call. = FALSE)
    }
  }

  # ---- Units (see make_forecast_configs.R): a unit is a fit group unless
  # --partition_by registry, where it is a group of --partition_source run to the
  # max member firm horizon of each launch year.
  registry_mode <- identical(opt$partition_by %||% "", "registry")
  units <- if (!registry_mode) {
    lapply(stats::setNames(nm = names(group_info)), function(g) {
      list(members = group_info[[g]]$isos, groups = g)
    })
  } else {
    if (!nzchar(opt$partition_source %||% "")) {
      stop(
        "[make_hindcast_configs] --partition_by registry requires --partition_source.",
        call. = FALSE
      )
    }
    preg <- load_grouping(opt$partition_source)
    assert_partition_labels(preg$group, opt$partition_source)
    fit_of <- stats::setNames(reg$group, reg$iso3)
    uncov <- setdiff(reg$iso3[reg$group %in% names(group_info)], preg$iso3)
    if (length(uncov)) {
      stop(
        sprintf(
          paste0(
            "[make_hindcast_configs] partition registry '%s' does not cover ",
            "%d country(ies) with deployed fits: %s"
          ),
          opt$partition_source,
          length(uncov),
          paste(uncov, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    u <- list()
    for (lab in sort(unique(preg$group))) {
      mem <- preg$iso3[preg$group == lab]
      mem <- mem[mem %in% names(fit_of) & fit_of[mem] %in% names(group_info)]
      if (!length(mem)) {
        next
      }
      u[[lab]] <- list(
        members = mem,
        groups = sort(unique(unname(fit_of[mem])))
      )
    }
    u
  }

  manifest_units <- list()
  country_years <- list() # by_country_from_panel: iso -> year -> horizons
  for (un in names(units)) {
    U <- units[[un]]
    isos <- U$members
    gi <- group_info[U$groups]
    brt <- unique(lapply(gi, `[[`, "base_rate_table_paths"))
    if (length(brt) > 1L) {
      stop(
        sprintf(
          paste0(
            "[make_hindcast_configs] unit '%s' mixes fit groups (%s) with ",
            "different anchor sidecars."
          ),
          un,
          paste(U$groups, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    base_rate_table_paths <- brt[[1L]]

    year_entries <- list()
    if (season_win) {
      # One run per launch month of the season year: the whole unit is simulated,
      # only the countries whose season starts then are scored (to their own firm
      # horizon); the run goes to the latest of those.
      for (y in years) {
        yc <- as.character(y)
        cw <- unlist(
          unname(lapply(gi, function(x) x$win[[yc]]$country)),
          recursive = FALSE
        )[isos]
        launches <- vapply(cw, function(w) format(w$launch), character(1))
        for (lc in sort(unique(launches))) {
          scored <- isos[launches == lc]
          horizon_chr <- max(vapply(
            cw[scored],
            function(w) format(w$horizon_end),
            character(1)
          ))
          ic <- as.Date(lc) %m-% months(1L)
          if (
            !is.null(cov) && (ic < cov$start || ic > (cov$end %m-% months(1L)))
          ) {
            stop(
              sprintf(
                paste0(
                  "[make_hindcast_configs] unit '%s', launch %s: IC month %s outside the ",
                  "smoothing coverage [%s, %s)."
                ),
                un,
                lc,
                format(ic),
                format(cov$start),
                format(cov$end)
              ),
              call. = FALSE
            )
          }
          key <- substr(lc, 1L, 7L)
          k6 <- sub("-", "", key)
          tag <- sprintf("%s_%s_%s", opt$run_tag_base, un, k6)
          rid <- sprintf("%s_%s_%s_hindcast", date_str, un, k6)
          path <- file.path(
            opt$config_outdir,
            sprintf("hindcast_%s_%s.yml", un, key)
          )
          assert_basename_length(c(path, tag, rid), "make_hindcast_configs")
          cfg <- patch_hindcast(
            tpl,
            tag,
            rid,
            isos,
            lc,
            horizon_chr,
            base_rate_table_paths
          )
          tr_of <- sw_tbl[match(scored, sw_tbl$iso3), , drop = FALSE]
          hdr <- generated_yaml_header(
            "analysis/utils/make_hindcast_configs.R",
            c(
              sprintf(
                paste0(
                  "# Hindcast (season windows): partition unit '%s' of registry '%s' (%d countries; ",
                  "fit groups %s), season %d, launch %s."
                ),
                un,
                opt$partition_source,
                length(isos),
                paste(U$groups, collapse = ", "),
                y,
                lc
              ),
              sprintf(
                "# handoff_month=%s (interior smoothing slice); observed ENSO; horizon.end=%s.",
                lc,
                horizon_chr
              ),
              sprintf(
                "# scored here: %s (the other unit members are scored in their own launch's run).",
                paste(scored, collapse = ", ")
              ),
              sprintf(
                "# [%s trough %s (%s) -> season %s .. %s, firm end %s%s]",
                scored,
                month.abb[tr_of$trough_month],
                tr_of$trough_source,
                vapply(cw[scored], function(w) format(w$launch), character(1)),
                vapply(
                  cw[scored],
                  function(w) format(w$season_end),
                  character(1)
                ),
                vapply(
                  cw[scored],
                  function(w) format(w$horizon_end),
                  character(1)
                ),
                vapply(
                  cw[scored],
                  function(w) {
                    if (isTRUE(w$floored)) {
                      sprintf(
                        "; start %s floored to %s",
                        month.abb[w$start_month],
                        month.abb[w$launch_month]
                      )
                    } else {
                      ""
                    }
                  },
                  character(1)
                )
              ),
              sprintf(
                "# smoothing_input (final historical smoothing): %s",
                smoothing_input
              )
            )
          )
          write_config(cfg, path, hdr)
          n_written <- n_written + 1L
          year_entries[[key]] <- list(
            launch = lc,
            horizon_end = horizon_chr,
            launch_year = y,
            run_tag = tag,
            run_id = rid,
            config = path,
            merged_dir = merged_dir_for(tag, rid, out_fmt),
            scored = as.list(scored)
          )
          for (iso in scored) {
            country_years[[iso]][[key]] <- list(
              launch_year = y,
              season_start = lc,
              season_end = format(cw[[iso]]$season_end),
              firm_horizon_end = format(cw[[iso]]$horizon_end),
              unit_horizon_end = horizon_chr
            )
          }
        }
      }
      manifest_units[[un]] <- list(
        members = as.list(isos),
        groups = as.list(U$groups),
        years = year_entries
      )
      message(sprintf(
        "[make_hindcast_configs] unit %-12s groups=%s members=%s  %d run(s)",
        un,
        paste(U$groups, collapse = ","),
        paste(isos, collapse = ","),
        length(year_entries)
      ))
      next
    }
    for (y in years) {
      yc <- as.character(y)
      launch_chr <- gi[[1L]]$win[[yc]]$launch_chr
      if (by_country) {
        # The unit runs to its latest MEMBER country (a group's max may sit in
        # another unit); the censor in every reader only truncates, so the unit
        # must cover each member's firm window.
        ch_all <- unlist(unname(lapply(gi, function(x) x$win[[yc]]$country_h)))
        cs_all <- unlist(unname(lapply(gi, function(x) x$win[[yc]]$country_se)))
        ch_u <- ch_all[isos]
        horizon_chr <- max(ch_u)
        if (anyNA(ch_u) || any(as.Date(ch_u) > as.Date(horizon_chr))) {
          stop(
            sprintf(
              "[make_hindcast_configs] unit '%s' %s: unit horizon %s does not cover every member's firm horizon.",
              un,
              yc,
              horizon_chr
            ),
            call. = FALSE
          )
        }
        for (iso in isos) {
          country_years[[iso]][[yc]] <- list(
            season_end = cs_all[[iso]],
            firm_horizon_end = ch_u[[iso]],
            unit_horizon_end = horizon_chr
          )
        }
      } else {
        horizon_chr <- max(vapply(
          gi,
          function(x) x$win[[yc]]$horizon_chr,
          character(1)
        ))
      }

      tag <- sprintf("%s_%s_%d", opt$run_tag_base, un, y)
      rid <- sprintf("%s_%s_%d_hindcast", date_str, un, y)
      path <- file.path(opt$config_outdir, sprintf("hindcast_%s_%d.yml", un, y))
      assert_basename_length(c(path, tag, rid), "make_hindcast_configs")

      cfg <- patch_hindcast(
        tpl,
        tag,
        rid,
        isos,
        launch_chr,
        horizon_chr,
        base_rate_table_paths
      )

      lag_lines <- vapply(
        names(gi),
        function(g) {
          x <- gi[[g]]
          sprintf(
            "# best_lag_occ=%d best_lag_dur=%d best_lag_mon=%d min_lag=%d; add_best_lag=%s.",
            x$lags[["occ"]],
            x$lags[["dur"]],
            x$lags[["mon"]],
            x$min_lag,
            add_best_lag
          )
        },
        character(1)
      )
      hdr <- generated_yaml_header(
        "analysis/utils/make_hindcast_configs.R",
        c(
          if (!registry_mode) {
            sprintf(
              "# Hindcast: ENSO group '%s' (%d countries), launch %s.",
              un,
              length(isos),
              launch_chr
            )
          } else {
            sprintf(
              "# Hindcast: partition unit '%s' of registry '%s' (%d countries; fit groups %s), launch %s.",
              un,
              opt$partition_source,
              length(isos),
              paste(U$groups, collapse = ", "),
              launch_chr
            )
          },
          sprintf(
            "# handoff_month=%s (interior smoothing slice); observed ENSO; horizon.end=%s.",
            launch_chr,
            horizon_chr
          ),
          if (!registry_mode) {
            lag_lines
          } else {
            paste0(
              sprintf(
                "# [%s firm end %s] ",
                names(gi),
                vapply(gi, function(x) x$win[[yc]]$horizon_chr, character(1))
              ),
              sub("^# ", "", lag_lines)
            )
          },
          if (by_country) {
            sprintf(
              "# [%s season end %s (cholera-year end, >= %d months after launch) -> firm end %s]",
              isos,
              vapply(
                isos,
                function(i) country_years[[i]][[yc]]$season_end,
                character(1)
              ),
              min_win,
              vapply(
                isos,
                function(i) country_years[[i]][[yc]]$firm_horizon_end,
                character(1)
              )
            )
          },
          sprintf(
            "# smoothing_input (final historical smoothing): %s",
            smoothing_input
          )
        )
      )
      write_config(cfg, path, hdr)
      n_written <- n_written + 1L

      year_entries[[yc]] <- list(
        launch = launch_chr,
        horizon_end = horizon_chr,
        run_tag = tag,
        run_id = rid,
        config = path,
        merged_dir = merged_dir_for(tag, rid, out_fmt)
      )
    }

    if (!registry_mode) {
      g <- un
      x <- gi[[g]]
      manifest_groups[[g]] <- list(
        best_lag_occ = x$lags[["occ"]],
        best_lag_dur = x$lags[["dur"]],
        best_lag_mon = x$lags[["mon"]],
        min_lag = x$min_lag,
        n_countries = length(isos),
        iso3 = as.list(isos),
        occ_fit_file = x$bl$occ$file,
        dur_fit_file = x$bl$dur$file,
        mon_fit_file = x$bl$mon$file,
        season_end_month = if (g %in% names(season_end_map)) {
          season_end_map[[g]]
        } else {
          NULL
        },
        years = year_entries
      )
      message(sprintf(
        "[make_hindcast_configs] %-12s occ=%d dur=%d mon=%d min=%d  %d year(s)",
        g,
        x$lags[["occ"]],
        x$lags[["dur"]],
        x$lags[["mon"]],
        x$min_lag,
        length(year_entries)
      ))
    } else {
      manifest_units[[un]] <- list(
        members = as.list(isos),
        groups = as.list(U$groups),
        years = year_entries
      )
      message(sprintf(
        "[make_hindcast_configs] unit %-12s groups=%s members=%s  %d year(s)",
        un,
        paste(U$groups, collapse = ","),
        paste(isos, collapse = ","),
        length(year_entries)
      ))
    }
  }

  # Registry mode: fit groups carry their per-year FIRM windows (no run fields).
  if (registry_mode) {
    for (g in names(group_info)) {
      x <- group_info[[g]]
      manifest_groups[[g]] <- list(
        best_lag_occ = x$lags[["occ"]],
        best_lag_dur = x$lags[["dur"]],
        best_lag_mon = x$lags[["mon"]],
        min_lag = x$min_lag,
        n_countries = length(x$isos),
        iso3 = as.list(x$isos),
        occ_fit_file = x$bl$occ$file,
        dur_fit_file = x$bl$dur$file,
        mon_fit_file = x$bl$mon$file,
        season_end_month = if (g %in% names(season_end_map)) {
          season_end_map[[g]]
        } else {
          NULL
        },
        units = as.list(names(units)[vapply(
          units,
          function(U) g %in% U$groups,
          logical(1)
        )]),
        years = lapply(x$win, function(w) {
          list(launch = w$launch_chr, horizon_end = w$horizon_chr)
        })
      )
    }
  }

  # by_country_from_panel: a group has no single season end (its members' are
  # under countries.<iso>); its horizon_end per year is the max over its members.
  if (by_country) {
    for (g in names(manifest_groups)) {
      manifest_groups[[g]]$season_end_month <- NULL
    }
  }

  if (!length(manifest_groups)) {
    stop(
      "No group has a deployed occurrence, duration AND monthly fit in ",
      opt$fits_dir,
      " (episodes.key = ",
      ep_key,
      ") -- nothing generated."
    )
  }

  manifest <- list(
    generated_at = format(Sys.time()),
    grouping = opt$grouping,
    launch_month = launch_mmdd,
    launch_years = as.list(years),
    add_best_lag = add_best_lag,
    horizon_months = opt$horizon_months,
    smoothing_input = smoothing_input,
    climate_panel = opt$climate_panel,
    skipped_groups = as.list(skipped),
    groups = manifest_groups
  )
  if (registry_mode) {
    manifest <- c(list(schema = 2L), manifest)
    manifest$grouping_digest <- registry_digest(opt$grouping)
    manifest$partition_by <- "registry"
    manifest$partition_source <- opt$partition_source
    manifest$partition_digest <- registry_digest(opt$partition_source)
    manifest$contagion_partition_source <- opt$nbr_partition_source %||%
      opt$partition_source
    manifest$units <- manifest_units
    manifest$countries <- lapply(
      stats::setNames(nm = unlist(lapply(units, `[[`, "members"))),
      function(iso) {
        list(
          group = reg$group[reg$iso3 == iso],
          unit = names(units)[vapply(
            units,
            function(U) iso %in% U$members,
            logical(1)
          )]
        )
      }
    )
  }
  season_rows <- NULL
  if (by_country) {
    # Per-country season ends and firm windows (sim_io.R::manifest_country_horizons()
    # reads countries.<iso>.years.<y>.firm_horizon_end before the group's).
    fit_of <- stats::setNames(reg$group, reg$iso3)
    cty <- list()
    for (un in names(units)) {
      for (iso in units[[un]]$members) {
        r <- se_tbl[se_tbl$iso3 == iso, , drop = FALSE]
        cty[[iso]] <- list(
          group = unname(fit_of[[iso]]),
          unit = un,
          season_end_month = r$season_end_month,
          min_seas_mode = r$min_seas_mode,
          modal_share = r$modal_share,
          years = country_years[[iso]]
        )
        season_rows <- rbind(
          season_rows,
          data.frame(
            iso3 = iso,
            group = unname(fit_of[[iso]]),
            unit = un,
            launch_year = as.integer(names(country_years[[iso]])),
            launch = vapply(
              names(country_years[[iso]]),
              function(yc) {
                group_info[[fit_of[[iso]]]]$win[[yc]]$launch_chr
              },
              character(1)
            ),
            season_end_month = r$season_end_month,
            min_seas_mode = r$min_seas_mode,
            modal_share = r$modal_share,
            n_gid = r$n_gid,
            n_distinct = r$n_distinct,
            group_season_end_month = if (
              fit_of[[iso]] %in% names(season_end_map)
            ) {
              season_end_map[[fit_of[[iso]]]]
            } else {
              NA_integer_
            },
            season_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "season_end"
            ),
            firm_horizon_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "firm_horizon_end"
            ),
            unit_horizon_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "unit_horizon_end"
            ),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        )
      }
    }
    manifest$countries <- cty
    manifest$season_end <- season_mode
    manifest$season_min_window_months <- min_win
    manifest$season_end_source <- season_end_source(se_cal)
  }
  if (season_win) {
    # Per-country troughs, the floor cut and its onset share, and the one season
    # entry per season year under the key of the run that scores it
    # (sim_io.R::manifest_country_horizons() reads only those).
    fit_of <- stats::setNames(reg$group, reg$iso3)
    onm <- attr(sw_tbl, "onsets_by_month")
    num_or_null <- function(v) if (length(v) && !is.na(v)) v else NULL
    cty <- list()
    for (un in names(units)) {
      for (iso in units[[un]]$members) {
        r <- sw_tbl[sw_tbl$iso3 == iso, , drop = FALSE]
        w0 <- season_window_of(r$trough_month, years[1L], floor_m, 0L, FALSE)
        cut_share <- if (!is.null(onm) && iso %in% rownames(onm)) {
          season_cut_share(onm[iso, ], w0$start_month, floor_m)
        } else {
          NA_real_
        }
        cty[[iso]] <- list(
          group = unname(fit_of[[iso]]),
          unit = un,
          trough_month = r$trough_month,
          trough_source = r$trough_source,
          min_seas_mode = r$min_seas_mode,
          modal_share = r$modal_share,
          n_onsets = num_or_null(r$n_onsets),
          onset_trough = num_or_null(r$onset_trough),
          stability = num_or_null(r$stability),
          season_start_month = w0$start_month,
          launch_month = w0$launch_month,
          floored = w0$floored,
          months_cut = w0$months_cut,
          onset_share_cut = num_or_null(cut_share),
          years = country_years[[iso]]
        )
        ky <- names(country_years[[iso]])
        season_rows <- rbind(
          season_rows,
          data.frame(
            iso3 = iso,
            group = unname(fit_of[[iso]]),
            unit = un,
            run_key = ky,
            launch_year = vapply(
              country_years[[iso]],
              `[[`,
              integer(1),
              "launch_year"
            ),
            trough_month = r$trough_month,
            trough_source = r$trough_source,
            min_seas_mode = r$min_seas_mode,
            modal_share = r$modal_share,
            n_onsets = r$n_onsets,
            onset_trough = r$onset_trough,
            stability = r$stability,
            season_start_month = w0$start_month,
            launch_month = w0$launch_month,
            floored = w0$floored,
            months_cut = w0$months_cut,
            onset_share_cut = cut_share,
            season_start = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "season_start"
            ),
            season_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "season_end"
            ),
            firm_horizon_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "firm_horizon_end"
            ),
            unit_horizon_end = vapply(
              country_years[[iso]],
              `[[`,
              character(1),
              "unit_horizon_end"
            ),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        )
      }
    }
    manifest$countries <- cty
    manifest$season_end <- "season_window"
    manifest$launch_month <- NULL # per country (countries.<iso>.launch_month)
    manifest$season_window <- sw_src
    fl <- vapply(cty, function(x) isTRUE(x$floored), logical(1))
    message(sprintf(
      paste0(
        "[make_hindcast_configs] season windows: %d countries, %d with the start ",
        "floored at %s (%s)"
      ),
      length(cty),
      sum(fl),
      month.abb[floor_m],
      if (any(fl)) {
        paste(
          sprintf(
            "%s -%d mo, %s of onsets",
            names(cty)[fl],
            vapply(cty[fl], `[[`, integer(1), "months_cut"),
            vapply(
              cty[fl],
              function(x) {
                if (is.null(x$onset_share_cut)) {
                  "n/a"
                } else {
                  sprintf("%.0f%%", 100 * x$onset_share_cut)
                }
              },
              character(1)
            )
          ),
          collapse = "; "
        )
      } else {
        "none"
      }
    ))
  }
  # Population year: only when not the default.
  manifest <- manifest_add_provenance(
    manifest,
    "pop_year",
    sim_pop_year_spec(tpl$simulator),
    .SIM_POP_YEAR_DEFAULT
  )
  dir.create(dirname(opt$manifest), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(manifest, opt$manifest)
  if (!is.null(season_rows)) {
    message(sprintf(
      "[make_hindcast_configs] season ends -> %s",
      write_season_ends_csv(season_rows, opt$manifest)
    ))
  }

  message(
    "\n[make_hindcast_configs] ----------------------------------------------"
  )
  message(sprintf(
    "  groups generated : %d  (%s)",
    length(manifest_groups),
    paste(names(manifest_groups), collapse = ", ")
  ))
  message(sprintf(
    "  configs written  : %d  (%d %s x %d years%s)",
    n_written,
    length(units),
    if (registry_mode) "units" else "groups",
    length(years),
    if (season_win) " x launch months" else ""
  ))
  if (length(skipped)) {
    message(sprintf("  groups skipped   : %s", paste(skipped, collapse = ", ")))
  }
  message(sprintf("  configs -> %s", opt$config_outdir))
  message(sprintf("  manifest -> %s", opt$manifest))
  message("[make_hindcast_configs] done.")

  invisible(manifest)
}

# ---------------------------------------------------------------------------
# Thin CLI wrapper — standalone `Rscript analysis/utils/make_hindcast_configs.R`
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  opt <- parse_args(OptionParser(option_list = option_list))
  req <- c(
    template = "--template",
    run_tag_base = "--run_tag_base",
    config_outdir = "--config_outdir"
  )
  missing_req <- req[vapply(
    names(req),
    function(k) is.null(opt[[k]]),
    logical(1)
  )]
  if (length(missing_req)) {
    stop(
      "[make_hindcast_configs] ",
      paste(missing_req, collapse = ", "),
      " required (no run-specific default). Normally run via: Rscript ",
      "analysis/utils/deploy_prod_run.R -C analysis/configs/deploy/",
      "prod_run<N>.yml --stage gen_hindcast",
      call. = FALSE
    )
  }
  generate_hindcast_configs(opt)
}
