# 03_prepare_sim_panel_projection.R — build a continuous per-(gid, month)
#   climate and ENSO panel for the simulator forecast horizon.
#
# Usage:
#   Rscript analysis/03_prepare_sim_panel_projection.R \
#     --cal_start 2023-05-01 --cal_end 2026-12-01
#
# Details:
#   Companion to 03_prepare_sim_panel.R (historical replay); --cal_end may
#   pass the climate DB max, where climate columns are zero-filled (the
#   deployed fits are ENSO-only). ENSO: observed ONI up to the cache max,
#   then a CFSv2 point estimate (--forecast_source cfsv2, default; c3s kept),
#   then --enso_continuation past the last CFSv2 target. Writes
#   make_sim_panel_projection_file(opt) + its .policy.rds stamp.

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(DBI)
library(RSQLite)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/simulator_core.R") # derive_min_seas_from_panel
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------

option_list <- list(
  make_option(
    c("-b", "--baseline"),
    default = "1990_2019",
    type = "character",
    help = "Baseline period; must match the occurrence fit [default: %default]"
  ),
  make_option(
    "--cal_start",
    default = "2025-01-01",
    type = "character",
    help = "First calendar month of the simulation horizon [default: %default]"
  ),
  make_option(
    "--cal_end",
    default = "2027-12-01",
    type = "character",
    help = "Last calendar month of the simulation horizon [default: %default]"
  ),
  make_option(
    "--warmup_months",
    default = 12L,
    type = "integer",
    help = "Months prepended before cal_start for warm-up state [default: %default]"
  ),
  make_option(
    "--enso_source",
    default = "oni",
    type = "character",
    help = "ENSO index source for historical months: 'oni' [default: %default]"
  ),
  make_option(
    "--forecast_source",
    default = "cfsv2",
    type = "character",
    help = "Post-cache ENSO point-estimate source: 'cfsv2' (NOAA CPC ensemble, matches simulator default) | 'c3s' (PNG plume median) [default: %default]"
  ),
  make_option(
    "--c3s_plume",
    default = "data/precalc_plume_quantiles_ENSO_NINO34_202605_branded.png",
    type = "character",
    help = "Path to C3S plume PNG (or companion .rds/.csv); only used when --forecast_source c3s [default: %default]"
  ),
  make_option(
    "--cfsv2_url",
    default = "https://www.cpc.ncep.noaa.gov/products/CFSv2/dataInd3/nino34Mon.nc",
    type = "character",
    help = "CFSv2 Nino3.4 monthly NetCDF URL [default: %default]"
  ),
  make_option(
    "--cfsv2_init",
    default = "latest",
    type = "character",
    help = "CFSv2 issuance to use: 'latest' or 'YYYY-MM' [default: %default]"
  ),
  make_option(
    "--cfsv2_rebase",
    default = "seam",
    type = "character",
    help = "CFSv2 seam-anchor rebase: 'seam' | 'none' | numeric offset [default: %default]"
  ),
  make_option(
    "--cfsv2_correction",
    default = "auto",
    type = "character",
    help = paste(
      "CFSv2 correction: auto | none | sd | target_peak | calibration.",
      "'auto' = target_peak when --cfsv2_target_peak is set, else sd",
      "(the historical behaviour) [default: %default]"
    )
  ),
  make_option(
    "--cfsv2_sd_correction",
    default = TRUE,
    action = "store_true",
    help = "Legacy flag (always TRUE); the mode is chosen by --cfsv2_correction [default: %default]"
  ),
  make_option(
    "--cfsv2_target_peak",
    default = NA,
    type = "double",
    help = "Peak-target rescale of the raw CFSv2 anomaly; NA = unset [default: %default]"
  ),
  make_option(
    "--cfsv2_calibration",
    default = "generated_data/enso_cfsv2_calibration.rds",
    type = "character",
    help = "Calibration RDS for --cfsv2_correction calibration (03_prepare_enso_cfsv2_calibration.R) [default: %default]"
  ),
  make_option(
    "--cfsv2_calibration_target",
    default = "ersst",
    type = "character",
    help = "Calibration target: ersst | oisst [default: %default]"
  ),
  make_option(
    "--cfsv2_n_rep",
    default = 2000L,
    type = "integer",
    help = "Number of Monte-Carlo member draws averaged into the CFSv2 point estimate [default: %default]"
  ),
  make_option(
    "--seed",
    default = 20260701L,
    type = "integer",
    help = "Seed for the CFSv2 point-estimate Monte-Carlo average [default: %default]"
  ),
  make_option(
    "--refresh_enso",
    default = TRUE,
    action = "store_true",
    help = "Re-download the ENSO index even if a cache exists [default: %default]"
  ),
  make_option(
    "--no_refresh_enso",
    action = "store_false",
    dest = "refresh_enso",
    help = "Reuse the cached ONI / Nino3.4 / CFSv2 files (generated_data/enso_*_cache.rds) instead of re-downloading them, e.g. to extend --cal_end without moving the ENSO inputs"
  ),
  make_option(
    "--enso_max_lead",
    default = 0L,
    type = "integer",
    help = paste(
      "Also write enso_lead1..K for deployed ENSO-lead fits, from the panel's own",
      "filled enso_lag0 series (observed ONI, then the forecast fill above):",
      "enso_lead<k>(t) = enso_lag0(t + k). 0 = none [default: %default]"
    )
  ),
  # ENSO past the last CFSv2 target (enso.R::sample_enso_cfsv2_trajectories(beyond_cfsv2 =)):
  # set these to the deploy's enso_continuation so the panel's point estimate and
  # its lead columns past the CFSv2 span match the forward run.
  make_option(
    "--enso_continuation",
    default = "neutral",
    type = "character",
    help = "ENSO past the CFSv2 span: neutral | ar1 | analog [default: %default]"
  ),
  make_option(
    "--enso_continuation_library",
    default = "generated_data/enso_continuation_library.qs",
    type = "character",
    help = "analog: library from 03_prepare_enso_continuation_library.R [default: %default]"
  ),
  make_option(
    "--enso_continuation_k",
    default = 30L,
    type = "integer",
    help = "analog: nearest analogs per seam state [default: %default]"
  ),
  make_option(
    "--enso_continuation_blend_months",
    default = 3,
    type = "double",
    help = "analog: e-folding (months) of the seam offset [default: %default]"
  ),
  make_option(
    "--enso_continuation_min_sep",
    default = 24L,
    type = "integer",
    help = "analog: months between analogs within a model [default: %default]"
  ),
  make_option(
    "--enso_continuation_max_per_model",
    default = Inf,
    type = "double",
    help = "analog: cap on analogs per piControl model (Inf = none) [default: %default]"
  ),
  make_option(
    "--enso_continuation_exclude_models",
    default = "CM2.1-NMME",
    type = "character",
    help = "analog: comma-separated piControl models left out ('' = none) [default: %default]"
  ),
  make_option(
    "--enso_continuation_seam",
    default = "last_cfsv2",
    type = "character",
    help = paste(
      "analog: seam of the continuation, last_cfsv2 (a 2-month ONI mean) |",
      "last_full_oni (the month before; review item 8) [default: %default]"
    )
  )
)

option_list <- c(
  option_list,
  list(
    make_option(
      "--episode_key",
      default = "lp_season",
      type = "character",
      help = paste0(
        "lp_season | gid_season: which survival panel defines ",
        "the gid scope and per-gid season calendar; the ",
        "output carries the matching _ep-gs token [default: %default]"
      )
    ),
    # Overwrite guard: a panel of another prep stamp is not replaced
    # without --force.
    prep_force_option("projection panel")
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$episode_key %in% EPISODE_KEYS) {
  stop("--episode_key must be one of ", paste(EPISODE_KEYS, collapse = " | "))
}
opt$episodes <- list(key = opt$episode_key)
print(opt)

cal_start <- as.Date(opt$cal_start)
cal_end <- as.Date(opt$cal_end)
warmup_start <- cal_start %m-% months(as.integer(opt$warmup_months))

CLIMATE_DB_MAX <- as.Date("2023-12-01") # canonical DB upper boundary

message("-- Building continuous sim projection panel")
message("   baseline:      ", opt$baseline)
message("   window:        ", warmup_start, " – ", cal_end)
message("   warmup:        ", warmup_start, " – ", cal_start - 1L)
message("   climate DB max:", CLIMATE_DB_MAX)
message("   post-DB cols:  zero-filled (ENSO-only model)")

if (cal_end < cal_start) {
  stop("--cal_end must be >= --cal_start.")
}

# 3. Scope + min_seas from survival panel --------------------------------------

message("-- Loading survival panel for gid scope and min_seas")
surv_panel <- readRDS(make_surv_dat_file(opt))
gids_meta <- surv_panel %>% distinct(gid, country, AFRO_region)
min_seas_tbl <- derive_min_seas_from_panel(surv_panel)

message("   ", nrow(gids_meta), " gids in scope")

# 4. Continuous month grid (warmup_start → cal_end) ----------------------------

message("-- Building month grid")
months_seq <- seq.Date(warmup_start, cal_end, by = "1 month")

panel <- expand_grid(gid = gids_meta$gid, t_start = months_seq) %>%
  left_join(gids_meta, by = "gid") %>%
  left_join(min_seas_tbl, by = "gid") %>%
  mutate(
    t_start_seas = yday(t_start %m+% months(-(min_seas - 1L))) - 1L
  )

message(
  "   ",
  nrow(panel),
  " rows  (",
  nrow(gids_meta),
  " gids × ",
  length(months_seq),
  " months)"
)

# 5. Climate: join from DB for months ≤ CLIMATE_DB_MAX; zero-fill beyond -------
# Climate columns (mean_precip, dev_precip, mean_temp_daymean,
# dev_temp_daymean) are zero-filled for months beyond the DB max (2023-12-01).
# The occurrence model is ENSO-only, so zero is explicit and safe.

message(
  "-- Querying climate_stats and climate_extremes for months ≤ ",
  CLIMATE_DB_MAX
)
all_gids <- unique(panel$gid)
clim_stats <- get_stats_gid(gids = all_gids)
clim_extrms <- get_extremes_gid(gids = all_gids, baseline = opt$baseline)

clim_dev <- clim_stats %>%
  pivot_longer(cols = c(get_precip_vars(), get_temp_vars())) %>%
  inner_join(
    clim_extrms %>% select(gid, month, name, mean, sd),
    by = c("gid", "month", "name")
  ) %>%
  mutate(dev = value - mean, dev_sd = (value - mean) / sd) %>%
  select(gid, month, name, mean, dev, dev_sd) %>%
  pivot_wider(names_from = "name", values_from = c("mean", "dev", "dev_sd"))

panel <- panel %>%
  left_join(clim_stats, by = c("gid", "t_start" = "month")) %>%
  left_join(clim_dev, by = c("gid", "t_start" = "month"))

# Zero-fill climate columns for months beyond the DB max.
clim_req_cols <- c(
  "mean_precip",
  "dev_precip",
  "mean_temp_daymean",
  "dev_temp_daymean"
)
panel <- panel %>%
  mutate(across(
    all_of(clim_req_cols),
    ~ if_else(t_start > CLIMATE_DB_MAX, 0.0, .)
  ))

message(
  "   Climate joined + zero-filled for ",
  sum(panel$t_start > CLIMATE_DB_MAX),
  " post-DB rows (",
  length(unique(panel$gid[panel$t_start > CLIMATE_DB_MAX])),
  " gids)"
)

# 6. NA assertion: climate only within [cal_start, CLIMATE_DB_MAX] -------------
# NA assertions: climate only within [cal_start, CLIMATE_DB_MAX]; every
# enso_lag* column asserted non-NA for the full [cal_start, cal_end] window
# (section 8).

in_clim_window <- panel %>%
  filter(t_start >= cal_start, t_start <= CLIMATE_DB_MAX)
if (nrow(in_clim_window) > 0) {
  na_counts <- in_clim_window %>%
    summarise(across(
      all_of(clim_req_cols),
      ~ sum(is.na(.)),
      .names = "na_{.col}"
    ))
  if (any(na_counts > 0)) {
    print(na_counts)
    stop(
      "NA climate values in required columns within [cal_start, CLIMATE_DB_MAX] — ",
      "check DB coverage (climate_stats.sqlite / climate_extremes.sqlite)"
    )
  }
  message(
    "   Climate NA check passed for ",
    nrow(in_clim_window),
    " in-climate-window rows"
  )
}

# 7. ENSO: historical ONI + post-cache point estimate (CFSv2 default) ----------
# Beyond the ONI cache max: a deterministic point estimate from NOAA CPC
# CFSv2 (default, `forecast_source: "cfsv2"`, matching the live simulator
# default in config_defaults.yml) — built by Monte-Carlo-averaging
# `sample_enso_cfsv2_trajectories()` (the same tested function
# `utils/simulations.R` uses for the stochastic per-replicate forecast) over
# many member draws, so the panel's point estimate and the simulator's
# per-replicate mean share identical rescale/rebase treatment. The C3S plume
# path (`--forecast_source c3s`) is kept side-by-side (not the default).
# Per-replicate stochastic ENSO is still injected at simulation time (the
# simulator kernel), not baked into the panel.
# Past the last CFSv2 target month the fill follows --enso_continuation
# (neutral | ar1 | analog; set it to the deploy's enso_continuation). The
# extended smoothing panel (03_prepare_sim_panel_extended.R) reads its lead
# columns past the observed ONI from this panel, so smoothing sees the
# continuation only when smoothing end + max lead > last CFSv2 target.

message("-- Fetching historical ENSO index (", opt$enso_source, ")")
enso_dat <- fetch_enso_index(
  source = opt$enso_source,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_enso)
)

oni_max_date <- with(
  enso_dat,
  as.Date(sprintf(
    "%d-%02d-01",
    year[which.max(year * 12 + month)],
    month[which.max(year * 12 + month)]
  ))
)
message("   ONI cache max date: ", oni_max_date)

# Join ENSO lags from historical ONI (baseline fill for every lag column;
# the branches below overwrite the post-cache months).
panel <- join_enso_lags(panel, enso_dat, date_col = "t_start")

.cfsv2_mode <- NA_character_ # resolved CFSv2 correction, stamped below (NA: no CFSv2 used)
if (identical(opt$forecast_source, "cfsv2")) {
  # --- CFSv2 point estimate ------------------------------------------------
  # Reuses sample_enso_cfsv2_trajectories() unchanged (the exact function the
  # live simulator calls in utils/simulations.R, same defaults as
  # config_defaults.yml simulator.enso.*), so the panel's point estimate and
  # the simulator's per-replicate stochastic forecast share identical
  # rescale/rebase treatment and identical historical/lag stitching — this
  # also fixes a prior bug where the plume-median overwrite clobbered
  # already-valid historical ONI whenever the *forecast month* (not the
  # *lagged* month) was beyond the ONI cache. The point estimate itself is a
  # Monte-Carlo average of `cfsv2_n_rep` member draws (member selection is
  # the only randomness in that function; the rescale/rebase transform is
  # identical across members, so the average converges to the corrected
  # ensemble-mean trajectory).
  message("-- Fetching raw Nino3.4 index (for CFSv2 rolling-mean base series)")
  nino34_dat <- tryCatch(
    {
      fetch_enso_index(
        source = "nino34",
        cache_dir = "generated_data/",
        refresh = isTRUE(opt$refresh_enso)
      )
    },
    error = function(e) {
      warning(
        "[03d] fetch_enso_index(nino34) failed (",
        conditionMessage(e),
        ") — falling back to ONI as the raw-anomaly proxy."
      )
      NULL
    }
  )

  cfsv2_dat <- tryCatch(
    {
      message("-- Fetching CFSv2 ensemble forecast: ", opt$cfsv2_url)
      fetch_enso_cfsv2_ensemble(
        url = opt$cfsv2_url,
        cache_dir = "generated_data/",
        init = opt$cfsv2_init,
        refresh = isTRUE(opt$refresh_enso)
      )
    },
    error = function(e) {
      warning(
        "[03d] fetch_enso_cfsv2_ensemble failed (",
        conditionMessage(e),
        ") — post-cache ENSO lags remain historical-ONI-only / NA beyond that."
      )
      NULL
    }
  )

  if (!is.null(cfsv2_dat)) {
    sim_months_all <- sort(unique(panel$t_start))
    mode <- opt$cfsv2_correction
    if (identical(mode, "auto")) {
      mode <- if (is.na(opt$cfsv2_target_peak)) "sd" else "target_peak"
    }
    mode <- match.arg(mode, c("none", "sd", "target_peak", "calibration"))
    .cfsv2_mode <- mode # stamped below (enso.R::cfsv2_correction_provenance())
    if (mode == "target_peak" && is.na(opt$cfsv2_target_peak)) {
      stop("[03d] --cfsv2_correction target_peak needs --cfsv2_target_peak")
    }
    message("-- CFSv2 correction: ", mode)
    .excl_models <- trimws(strsplit(
      opt$enso_continuation_exclude_models %||% "",
      ","
    )[[1]])
    .excl_models <- .excl_models[nzchar(.excl_models)]
    enso_pt <- sample_enso_cfsv2_trajectories(
      cfsv2_data = cfsv2_dat,
      n_rep = as.integer(opt$cfsv2_n_rep),
      historical_oni = enso_dat,
      historical_nino34 = nino34_dat,
      sim_months = sim_months_all,
      max_lag = 12L,
      rebase = opt$cfsv2_rebase,
      sd_correction = mode == "sd",
      target_peak = if (mode == "target_peak") opt$cfsv2_target_peak else NULL,
      calibration = if (mode == "calibration") opt$cfsv2_calibration else NULL,
      calibration_target = opt$cfsv2_calibration_target,
      seed = as.integer(opt$seed),
      beyond_cfsv2 = opt$enso_continuation,
      analog_library = opt$enso_continuation_library,
      analog_k = opt$enso_continuation_k,
      analog_blend_months = opt$enso_continuation_blend_months,
      analog_min_sep = opt$enso_continuation_min_sep,
      analog_max_per_model = opt$enso_continuation_max_per_model,
      analog_exclude_models = .excl_models,
      analog_seam = opt$enso_continuation_seam
    )
    enso_lag_cols <- grep("^enso_lag[0-9]+$", names(panel), value = TRUE)
    for (ec in enso_lag_cols) {
      if (is.null(enso_pt[[ec]])) {
        next
      }
      pt_est <- rowMeans(enso_pt[[ec]])
      lkp <- stats::setNames(pt_est, as.character(sim_months_all))
      panel[[ec]] <- as.numeric(lkp[as.character(panel$t_start)])
    }
    message(
      "   Filled ",
      length(enso_lag_cols),
      " enso_lag* columns from CFSv2 (n_rep=",
      opt$cfsv2_n_rep,
      " Monte-Carlo point estimate; historical months reconstructed via the ",
      "same rolling-mean pipeline as the live simulator)."
    )
  } else {
    warning("[03d] No CFSv2 data — ENSO lags beyond ONI cache remain NA.")
  }
} else {
  # --- C3S plume median (side-by-side legacy path; not the default) -------
  plume_dat <- tryCatch(
    {
      message("-- Reading C3S plume: ", opt$c3s_plume)
      read_enso_c3s_plume(opt$c3s_plume)
    },
    error = function(e) {
      warning(
        "[03d] read_enso_c3s_plume failed (",
        conditionMessage(e),
        ") — post-cache ENSO will be NA (check plume file path)."
      )
      NULL
    }
  )

  # Build a Niño 3.4 median look-up for post-cache months from the plume.
  plume_median_lkp <- if (!is.null(plume_dat)) {
    plume_dat %>%
      group_by(date) %>%
      summarise(
        nino34_median = {
          qs <- q_level
          vals <- nino34_anomaly
          idx50 <- which.min(abs(qs - 0.50))
          if (length(idx50) > 0) {
            vals[idx50[1]]
          } else {
            stats::median(vals, na.rm = TRUE)
          }
        },
        .groups = "drop"
      ) %>%
      {
        stats::setNames(.$nino34_median, as.character(.$date))
      }
  } else {
    NULL
  }

  # Overwrite enso_lag* for months beyond ONI cache using C3S plume median.
  post_cache_rows <- panel$t_start > oni_max_date

  if (any(post_cache_rows) && !is.null(plume_median_lkp)) {
    message(
      "   Overwriting ENSO lags for ",
      sum(post_cache_rows),
      " post-cache rows from C3S plume median"
    )
    enso_lag_cols <- grep("^enso_lag[0-9]+$", names(panel), value = TRUE)
    for (ec in enso_lag_cols) {
      lag_n <- as.integer(sub("enso_lag", "", ec))
      lagged_dates <- panel$t_start[post_cache_rows] %m-% months(lag_n)
      plume_vals <- plume_median_lkp[as.character(lagged_dates)]
      panel[[ec]][post_cache_rows] <- as.numeric(plume_vals)
    }
  } else if (any(post_cache_rows)) {
    warning(
      "[03d] No C3S plume data — ENSO lags for post-cache months remain NA."
    )
  }
}

# E2. ENSO lead columns (deployed lead fits; enso.R::add_enso_leads_from_lag0()) --
# From the FILLED enso_lag0 above, so a lead month past the last ONI month
# reads the same forecast point estimate the panel uses for lag 0. The last k
# months of the panel have no t + k month and stay NA (warned below).
if (isTRUE(opt$enso_max_lead > 0L)) {
  panel <- add_enso_leads_from_lag0(panel, opt$enso_max_lead)
  message(
    "   Wrote enso_lead1..",
    opt$enso_max_lead,
    " from the filled enso_lag0 series"
  )
}

# 8. NA assertion: every enso_lag* column, full [cal_start, cal_end] -----------

in_full_window <- panel %>% filter(t_start >= cal_start, t_start <= cal_end)
enso_lag_cols_chk <- grep(
  "^enso_(lag|lead)[0-9]+$",
  names(in_full_window),
  value = TRUE
)
for (lag_nm in enso_lag_cols_chk) {
  n_na <- sum(is.na(in_full_window[[lag_nm]]))
  if (n_na > 0) {
    warning(
      sprintf("[03d] %s has %d NA(s) in [cal_start, cal_end]. ", lag_nm, n_na),
      "Simulation months requiring these lags will use 0. ",
      "Check ENSO cache coverage and ",
      opt$forecast_source,
      " forecast file."
    )
  } else {
    message(sprintf(
      "   %s NA check passed in [%s, %s]",
      lag_nm,
      cal_start,
      cal_end
    ))
  }
}

# 9. Warm-up lag columns from survival panel -----------------------------------

message("-- Joining warm-up lag columns from survival panel")
lag_cols <- c(
  "ongoing_outbreak_lag1",
  "ongoing_outbreak_country_lag1",
  "ongoing_outbreak_tu_lag1",
  "time_to_last_outbreak"
)
surv_lags <- surv_panel %>%
  select(gid, t_start, any_of(lag_cols)) %>%
  group_by(gid, t_start) %>%
  summarise(
    across(
      any_of(c("ongoing_outbreak_lag1", "ongoing_outbreak_country_lag1")),
      ~ as.integer(max(coalesce(., 0L))),
      .names = "{.col}"
    ),
    across(
      any_of("ongoing_outbreak_tu_lag1"),
      ~ {
        v <- .[!is.na(.)]
        if (length(v) == 0) NA_integer_ else as.integer(max(v))
      },
      .names = "{.col}"
    ),
    across(
      any_of("time_to_last_outbreak"),
      ~ suppressWarnings(min(., na.rm = TRUE)),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

panel <- panel %>%
  left_join(surv_lags, by = c("gid", "t_start")) %>%
  make_time_to_last_outbreak_cat()

# 10. Climate lag columns for DLNM cross-basis reconstruction at sim time ------
# Mirrors 03c: continuous calendar lags per gid, no season-year reset.
# For months beyond CLIMATE_DB_MAX the climate columns are zero-filled,
# so lag columns in that region propagate zeros — safe for the ENSO-only
# occurrence model. Lag k = NA for the first k rows of each gid.

message("-- Computing climate lag columns (lag1..lag11)")
clim_lag_vars <- grep("^(mean|dev)_", names(panel), value = TRUE)
clim_lag_vars <- clim_lag_vars[!grepl("^dev_sd_", clim_lag_vars)]
clim_lag_vars <- clim_lag_vars[!grepl("_lag\\d+$", clim_lag_vars)]

max_lag_clim <- 11L
lag_fns <- setNames(
  lapply(seq_len(max_lag_clim), function(k) {
    force(k)
    function(x) dplyr::lag(x, n = k)
  }),
  paste0("lag", seq_len(max_lag_clim))
)

panel <- panel %>%
  dplyr::arrange(gid, t_start) %>%
  dplyr::group_by(gid) %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(clim_lag_vars),
    lag_fns,
    .names = "{.col}_{.fn}"
  )) %>%
  dplyr::ungroup()

message(
  "   Added ",
  max_lag_clim * length(clim_lag_vars),
  " lag columns (",
  length(clim_lag_vars),
  " vars × lag1..",
  max_lag_clim,
  ")"
)

# 11. Save ---------------------------------------------------------------------
# Output file: make_sim_panel_projection_file(opt) ->
#   generated_data/sim_climate_panel_projection_base-<baseline>.rds

out_path <- make_sim_panel_projection_file(opt)
message("-- Saving to ", out_path)
# Episode record of the survival panel this panel was scoped from, so the
# simulator can refuse to pair it with fits/panels of another episode
# definition (utils/provenance.R::sim_episode_check()). Built before the panel is
# written (same fields and values).
.up <- prep_stamp_inherit(make_surv_dat_file(opt))
.stamp <- prep_policy_stamp(
  pop_source = .up$pop_source,
  composite_policy = .up$composite_policy,
  cover_threshold = .up$cover_threshold,
  composite_tiling = .up$composite_tiling,
  episode_key = .up$episode_key,
  episode_digest = .up$episode_digest,
  # cfsv2_correction / cfsv2_rebase_effective are appended only when not the
  # default (sd + seam; enso.R::cfsv2_correction_provenance()), so a default
  # stamp keeps its field set.
  counts = c(
    list(
      track = "sim_panel_projection",
      baseline = opt$baseline,
      enso_max_lead = as.integer(opt$enso_max_lead %||% 0L),
      enso_continuation = enso_continuation_provenance(list(
        method = opt$enso_continuation,
        library = opt$enso_continuation_library,
        k = opt$enso_continuation_k,
        blend_months = opt$enso_continuation_blend_months,
        min_sep = opt$enso_continuation_min_sep,
        max_per_model = opt$enso_continuation_max_per_model,
        exclude_models = trimws(strsplit(
          opt$enso_continuation_exclude_models %||% "",
          ","
        )[[1]]),
        seam = opt$enso_continuation_seam
      )),
      source_surv_panel = as.character(make_surv_dat_file(opt))
    ),
    cfsv2_correction_provenance(.cfsv2_mode, opt$cfsv2_rebase)
  )
)
# Overwrite guard (the rule of 03_prepare_occ.R): a panel of
# another stamp (another vintage's survival panel, another continuation or
# CFSv2 correction) or with none is not replaced without --force. ENSO values
# are not in the stamp: a re-download under the same stamp passes (back up and
# diff the panels).
assert_prep_stamp_overwrite(
  out_path,
  .stamp,
  force = isTRUE(opt$force),
  rerun_hint = "Rscript analysis/03_prepare_sim_panel_projection.R ... --force"
)
saveRDS(panel, file = out_path)
persist_prep_policy_stamp(out_path, .stamp)
message("   Done. ", nrow(panel), " rows written.")
message("   Date range: ", min(panel$t_start), " – ", max(panel$t_start))
message("   Gids: ", length(unique(panel$gid)))
