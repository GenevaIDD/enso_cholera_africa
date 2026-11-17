# occ_data.R — occurrence (onset-hazard) survival panel: build at prep time,
#   load and guard at fit time
#
# Provides:
#   make_surv_data()            survival rows of one location period / gid
#   gid_at_risk_years(), gid_season_surv_input()   gid-season builder inputs
#   national_zero_surv_rows()   national-zero months as at-risk rows
#   mask_occ_panel_rows()       occurrence-panel spec masks
#   guard_missing_climate_data(), guard_outbreak_episodes(),
#   assert_nz_rows_outside_episodes()   prep-time guards
#   load_surv_data()            read + preprocess the panel for a fit
#   prepare_occ_brms_data()     shared data preamble of the occurrence workers
#   assert_contagion_partition()   fit-time contagion-registry guard
#   make_climate_lags(), scale_precip(), scale_temp()   covariate helpers
#
# Details:
#   Builders used by 03_prepare_occ.R (and the other 03 preps, 01b); loaders
#   and guards by the 04 workers, 05 and 07. Callers source the dependencies
#   (config.R, paths.R, footprint.R, anchor.R as needed).

#' Load and preprocess the survival dataset for modelling.
#'
#' Applies the standard preprocessing pipeline (category construction, scaling,
#' factor coercion, NA filtering) and the country subset when opt$country != "all".
#' Also normalises opt$country to a filename-safe form (commas → dashes) when a
#' subset is requested.
#'
#' Called by 04a_, 04b_, 04c_, and 05_postprocess_results.R to eliminate the
#' identical ~20-line load + preprocess + country-filter block in each script.
#'
#' @param opt  resolved options list (must contain $baseline and $country)
#' @return list with two elements:
#'   $dat  pre-processed tibble ready for modelling
#'   $opt  updated options (opt$country may be mutated to filename-safe form)
#' Fit-time guard for the partition-restricted neighbour covariate
#' (`contagion.partition: group`, columns `nbr_lag1_r<rho>[_rn]_pg`). Those
#' columns depend on WHICH grouping registry prep partitioned by, which the
#' column name does not say; the prep policy stamp (03_prepare_occ.R,
#' `contagion_partition_source`) does. A fit whose `grouping.source` differs
#' from the stamp -- or a survival file prepped without the flag -- stops here
#' with the re-run command, instead of fitting on the wrong neighbourhoods.
#' No-op for every config without `contagion.partition: group`.
assert_contagion_partition <- function(opt, surv_path) {
  if (
    isTRUE(opt$occurrence$covariates$nbr_contagion_lag1) &&
      identical(opt$contagion$partition, "registry")
  ) {
    return(.assert_contagion_registry(opt, surv_path))
  }
  if (
    !isTRUE(opt$occurrence$covariates$nbr_contagion_lag1) ||
      !identical(opt$contagion$partition, "group")
  ) {
    return(invisible(TRUE))
  }
  want <- opt$grouping$source %||% ""
  got <- read_prep_policy_stamp(surv_path)$contagion_partition_source %||% ""
  if (!nzchar(want) || !identical(got, want)) {
    stop(
      sprintf(
        paste0(
          "contagion.partition = 'group' needs %s prepped with ",
          "--contagion_partition_source %s (stamp says '%s'). Re-run:\n",
          "  Rscript analysis/03_prepare_occ.R -b %s --contagion_partition_source %s"
        ),
        surv_path,
        if (nzchar(want)) want else "<grouping.source>",
        got,
        opt$baseline,
        if (nzchar(want)) want else "<grouping.source>"
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# contagion.partition: registry -- the survival file must have been prepped
# with contagion.partition_source among --contagion_partition_source, and the
# registry must not have changed since (content digest in the prep stamp).
.assert_contagion_registry <- function(opt, surv_path) {
  want <- opt$contagion$partition_source %||% ""
  if (!nzchar(want)) {
    stop(
      "contagion.partition = 'registry' requires contagion.partition_source.",
      call. = FALSE
    )
  }
  st <- read_prep_policy_stamp(surv_path)
  srcs <- strsplit(st$contagion_partition_sources %||% "", ",")[[1L]]
  digs <- strsplit(st$contagion_partition_digests %||% "", ",")[[1L]]
  rerun <- sprintf(
    "  Rscript analysis/03_prepare_occ.R -b %s --contagion_partition_source <...,>%s",
    opt$baseline,
    want
  )
  if (!want %in% srcs) {
    stop(
      sprintf(
        paste0(
          "contagion.partition = 'registry' needs %s prepped with registry '%s' ",
          "in --contagion_partition_source (stamp lists: '%s'). Re-run:\n%s"
        ),
        surv_path,
        want,
        paste(srcs, collapse = ","),
        rerun
      ),
      call. = FALSE
    )
  }
  got <- digs[match(want, srcs)]
  now <- registry_digest(want)
  if (
    !is.na(got) &&
      nzchar(got) &&
      got != "NA" &&
      !is.na(now) &&
      !identical(got, now)
  ) {
    stop(
      sprintf(
        paste0(
          "contagion registry '%s' changed since %s was prepped (digest %s, now %s): ",
          "the _part columns describe the old memberships. Re-run:\n%s"
        ),
        want,
        surv_path,
        got,
        now,
        rerun
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_surv_data <- function(opt) {
  # NOTE: local variable must be named `surv_dat_all` — scale_precip() and
  # scale_temp() reference that name from the calling environment via
  # colnames(surv_dat_all). Renaming it would cause those helpers to fail.
  surv_dat_all <- readRDS(make_surv_dat_file(opt = opt))
  assert_contagion_partition(opt, make_surv_dat_file(opt = opt))
  # Episode definition: the panel must be the one the config's episodes.key
  # names (utils/provenance.R). The record rides on the returned data so the
  # fit script can stamp it onto the fit.
  ep_rec <- assert_episode_key(opt, make_surv_dat_file(opt = opt), "occurrence")
  assert_occ_panel_variant(opt, make_surv_dat_file(opt = opt)) # national zero-fill + panel-mask fields

  # P3 (geometry-derived admin level), BEFORE the admin-level filter -- the
  # same position read_size_data()/load_monthly_data() use. Gated inline on the
  # occurrence key (default off) so the default path needs nothing from
  # utils/population.R, which the 04_fit_occ_* scripts source for this.
  if (isTRUE((opt$occurrence %||% list())$canonical_admin_level)) {
    surv_dat_all <- .apply_p3(
      surv_dat_all,
      opt,
      "occurrence",
      tag = "load_surv_data/P3"
    )
  }

  surv_dat_all <- surv_dat_all %>%
    # occurrence.ttlo_no_prior_level: + level "none" (load-time recode)
    make_time_to_last_outbreak_cat(
      no_prior_level = isTRUE((opt$occurrence %||% list())[[
        "ttlo_no_prior_level"
      ]])
    ) %>%
    dplyr::mutate(
      months = factor(months),
      quarter = factor(lubridate::quarter(t_start))
    ) %>%
    scale_precip() %>%
    scale_temp() %>%
    make_factors_anomalies() %>%
    apply_admin_level_filter(admin_level_opt = opt$occurrence$admin_level) %>%
    dplyr::mutate(country = factor(country), ADM1 = factor(ADM1)) %>%
    dplyr::filter(!is.na(t_start), !is.na(t_end))

  if (opt$country != "all") {
    countries_filter <- trimws(strsplit(opt$country, ",")[[1]])
    opt$country <- paste(countries_filter, collapse = "-") # filename-safe
    surv_dat_all <- surv_dat_all %>%
      dplyr::filter(as.character(country) %in% countries_filter) %>%
      dplyr::mutate(country = droplevels(country), ADM1 = droplevels(ADM1))
  }

  # Empirical baseline-rate anchor (utils/anchor.R), gated on
  # occurrence.offset_baserate -- OFF by default, so every config that exists
  # today is bit-inert here. The table is built on the FULL corpus, ahead of the
  # opt$country subset applied above, so a per-country or per-group fit still
  # resolves against the global cell/group/global fallback tiers rather than
  # against a shrunken sample of itself.
  if (anchor_enabled(opt, "occ")) {
    surv_dat_all <- attach_base_rate(surv_dat_all, opt, "occ")
  }

  attr(surv_dat_all, "episode_stamp") <- ep_rec
  list(dat = surv_dat_all, opt = opt)
}


#' Poisson-GLMM view of the occurrence panel (response, exposure, offset).
#'
#' The occurrence track is the only one of the four whose response and exposure
#' do not exist on disk: `03_prepare_occ.R` writes the survival panel, and the
#' Poisson recoding (`y`, `exposure`, `log_offset`) is derived at fit time. A
#' copy of that derivation in each consumer would let a baseline-rate anchor
#' silently disagree with the offset it is supposed to cancel: the anchor's
#' denominator MUST be the same object as `exp(log_offset)`, not a
#' re-derivation of it. Centralising it here makes that guarantee mechanical.
#'
#' @param dat occurrence panel (either the raw `make_surv_dat_file()` frame or the
#'   `load_surv_data()` output -- the three derived columns depend only on
#'   `outbreak` / `t_start_seas` / `t_end_seas`, which both carry).
#' @return `dat` with `y`, `exposure` and `log_offset` added.
.occ_poisson_view <- function(dat) {
  dplyr::mutate(
    dat,
    y = as.integer(outbreak),
    # Real calendar days of the row (make_surv_data()'s exposure_days). The
    # rotated-day span t_end_seas - t_start_seas is the length of the ROTATED
    # month, up to 3 days off the calendar month the row's climate is joined
    # on whenever min_seas != 1; it stays the fallback for panels
    # without exposure_days.
    exposure = pmax(
      dplyr::coalesce(
        if ("exposure_days" %in% names(dat)) {
          as.integer(dat$exposure_days)
        } else {
          NA_integer_
        },
        as.integer(t_end_seas) - as.integer(t_start_seas)
      ),
      1L
    ),
    log_offset = log(exposure)
  )
}

#' prepare_occ_brms_data
#'
#' Shared data preamble of the occurrence brms workers (04_fit_occ_brms.R,
#' 04_fit_occ_brms_cb.R, 04_fit_occ_brms_kl.R): load_surv_data() -> the episode stamp
#' -> the base-rate
#' anchor sidecar when occurrence.offset_baserate is on -> .occ_poisson_view()
#' -> the optional M-spline baseline columns. Same calls in the same order.
#'
#' Call-time dependencies (callers source them): anchor.R
#' (anchor_enabled(), persist_base_rate_table()), splines2 when basehaz == "ms".
#'
#' @param opt     resolved options
#' @param basehaz resolved baseline hazard (.occ_basehaz())
#' @return list(opt = <opt after load_surv_data()>, pois_dat, ms_basis (NULL
#'   unless "ms"), ep_rec = <episode stamp of the prep artefact>)
prepare_occ_brms_data <- function(opt, basehaz) {
  tmp <- load_surv_data(opt)
  # Episode record of the prep artefact the fit trains on (assert_episode_key()
  # in the loader); the worker stamps it onto the fit so the simulator can
  # refuse a fit whose episode definition differs from its paired artefacts.
  ep_rec <- attr(tmp$dat, "episode_stamp")
  opt <- tmp$opt

  # Anchor sidecar (utils/anchor.R); gated on occurrence.offset_baserate.
  if (anchor_enabled(opt, "occ")) {
    persist_base_rate_table(opt, "occ", tmp$dat)
  }

  pois_dat <- .occ_poisson_view(tmp$dat)

  ms_basis <- NULL
  if (identical(basehaz, "ms")) {
    ms_basis <- splines2::mSpline(
      pois_dat$t_start_seas,
      degree = 3,
      df = 5,
      intercept = FALSE
    )
    colnames(ms_basis) <- paste0("ms", seq_len(ncol(ms_basis)))
    pois_dat <- dplyr::bind_cols(pois_dat, tibble::as_tibble(ms_basis))
  }
  list(opt = opt, pois_dat = pois_dat, ms_basis = ms_basis, ep_rec = ep_rec)
}


# Data wrangling ----------------------------------------------------------

apply_offset_months <- function(date, offset) {
  format(date %m+% months(offset), "%Y-%m-%d")
}

#' .clim_dev_wide
#' Wide climatology / anomaly columns for a set of climate variables:
#' `mean_<v>` (baseline mean of the calendar month), `dev_<v>` (value minus
#' that mean) and `dev_sd_<v>` (the anomaly over the baseline sd). The pivot
#' make_surv_data() always ran inline, lifted out unchanged so the climate panel
#' builds the same columns.
#'
#' @param clim_stats wide climate stats (gid, month, <vars>), get_stats_gid()
#' @param clim_extrms long baseline stats (gid, month, name, mean, sd, ...),
#'   get_extremes_gid()
#' @param vars variable columns of `clim_stats` to pivot
#' @param by join keys between the long stats and `clim_extrms`: make_surv_data()
#'   reads one gid and joins on (month, name); multi-gid callers add gid
#' @return tibble (gid, month, mean_<v>, dev_<v>, dev_sd_<v> ...)
.clim_dev_wide <- function(
  clim_stats,
  clim_extrms,
  vars,
  by = c("gid", "month", "name")
) {
  clim_stats %>%
    tidyr::pivot_longer(cols = dplyr::all_of(vars)) %>%
    dplyr::inner_join(
      clim_extrms %>%
        dplyr::select(dplyr::all_of(setdiff(by, "name")), name, mean, sd),
      by = by
    ) %>%
    dplyr::mutate(dev = value - mean, dev_sd = (value - mean) / sd) %>%
    dplyr::select(gid, month, name, mean, dev, dev_sd) %>%
    tidyr::pivot_wider(
      names_from = "name",
      values_from = c("mean", "dev", "dev_sd")
    )
}

# Function to make data for analysis in survival
make_surv_data <- function(x, baseline) {
  # Get clim data
  clim_stats <- get_stats_gid(gid = x$gid[1])
  clim_extrms <- get_extremes_gid(gid = x$gid[1], baseline = baseline)

  clim_extrms_cat <- clim_extrms %>%
    select(-hilo_sd, -mean, -sd) %>%
    mutate(name = stringr::str_c("hlq_", name)) %>%
    tidyr::pivot_wider(names_from = "name", values_from = "hilo_q")

  clim_dev <- .clim_dev_wide(
    clim_stats,
    clim_extrms,
    c(get_precip_vars(), get_temp_vars()),
    by = c("month", "name")
  )

  tibble(
    days_seas = seq.Date(
      min(lubridate::floor_date(x$TL_seas, unit = "year")),
      min(x$TL_seas),
      by = "1 days"
    ),
    days = days_seas %m+% months(x$min_seas[1] - 1),
    months = lubridate::month(days)
  ) %>%
    group_by(months) %>%
    summarise(
      uid = x$uid[1],
      gid = x$gid[1],
      location_period_id = x$location_period_id[1],
      t_start = min(days),
      t_end = max(days),
      t_start_seas = min(lubridate::yday(days_seas)) - 1, # Remove 1 to havd open left intervals

      t_end_seas = max(lubridate::yday(days_seas))
    ) %>%
    mutate(
      outbreak = t_end == max(.$t_end),
      # Real days at risk in the calendar month: the whole month, except
      # the season's last row (onset or coverage end), cut at t_end.
      exposure_days = ifelse(
        outbreak,
        as.integer(t_end - t_start) + 1L,
        unname(as.integer(lubridate::days_in_month(t_start)))
      )
    ) %>%
    arrange(t_start) %>%
    left_join(clim_stats, by = c("gid" = "gid", "t_start" = "month")) %>%
    left_join(clim_extrms_cat, by = c("gid" = "gid", "t_start" = "month")) %>%
    left_join(clim_dev, by = c("gid" = "gid", "t_start" = "month"))
}

#' guard_missing_climate_data
#' Flags survival-data rows with missing climate covariates. make_surv_data()
#' left_joins outbreak-track rows onto climate_stats/climate_extremes, so any
#' row whose month falls outside the raw CHIRPS/ERA5 coverage window
#' (currently 1981-01 to 2024-12, climate_db_provenance.yml) is kept but has NA `precip`/`temp_daymean`
#' (and all derived dev_*/hlq_*/lag columns). This adds an explicit
#' `climate_missing` flag and warns, so 04_fit_* scripts don't silently pass
#' NA covariates into brms.
#'
#' @param dat survival data; must contain `precip`, `temp_daymean`, `t_start`,
#'   `location_period_id`
#'
guard_missing_climate_data <- function(dat) {
  dat <- dat %>%
    mutate(climate_missing = is.na(precip) | is.na(temp_daymean))

  n_missing <- sum(dat$climate_missing)

  if (n_missing > 0) {
    date_range <- range(dat$t_start[dat$climate_missing])
    n_lp <- dplyr::n_distinct(dat$location_period_id[dat$climate_missing])

    warning(
      n_missing,
      " / ",
      nrow(dat),
      " rows (",
      n_lp,
      " location periods) have ",
      "missing climate data (t_start ",
      format(date_range[1]),
      " to ",
      format(date_range[2]),
      "), outside CHIRPS/ERA5 coverage. ",
      "Rows are kept (flagged via `climate_missing`) but contain NA climate ",
      "covariates -- filter explicitly before model fitting.",
      call. = FALSE
    )
  }

  dat
}


#' guard_outbreak_episodes
#' Drops episode-level outbreak rows that fail plausibility checks (QA of the
#' Aug 2026 OutbreakExtractR refresh): non-positive/implausibly small case
#' counts and multi-year spans that are very likely merge artefacts. Written to
#' be called on `outbreak_metadata.rds` right after it is read, so every track
#' describes the identical corpus. Always LOUD: reports what it drops via
#' warning(), never silently.
#'
#' NOT CALLED by any pipeline script (QA only,
#' pinned by tests/test_duration_span_stamp.R [4]). On the duration track its
#' `max_duration_days` drop would cut exactly the long-episode tail the
#' cessation hazard estimates, so dur_data.R::make_duration_panel() warns about
#' episodes longer than 730 days and keeps them instead. What the guard WOULD
#' drop is counted by dur_data.R::outbreak_episode_guard_flags() (same thresholds
#' and classes) and stamped on the duration panel as `n_guard_cases_le0`,
#' `n_guard_cases_lt50`, `n_guard_gt730d` and `n_guard_any` (03_prepare_dur.R).
#'
#' @param dat episode-level data (one row per outbreak episode), as read from
#'   ~/secure/cholera_outbreak_data/outbreak_metadata.rds. Must contain `suspected_cases`,
#'   `TL`, `TR`, `location_period_id`.
#' @param min_cases drop episodes with `suspected_cases` below this (the
#'   documented ingest-time size filter is >=50; some slip through it via
#'   post-ingest aggregation in 01_preprocess_outbreaks.R). Default 50.
#' @param max_duration_days drop episodes spanning more than this many days
#'   (TR - TL), a proxy for stitching/merge artefacts. Default 730 (~2 years).
#' @param tag short prefix for the warning/message text (e.g. "occ" or "siz").
#'
guard_outbreak_episodes <- function(
  dat,
  min_cases = 50,
  max_duration_days = 730,
  tag = "prep"
) {
  duration_days <- as.numeric(dat$TR - dat$TL)

  bad_cases <- !is.na(dat$suspected_cases) & dat$suspected_cases <= 0
  bad_min_size <- !is.na(dat$suspected_cases) &
    dat$suspected_cases > 0 &
    dat$suspected_cases < min_cases
  bad_duration <- !is.na(duration_days) & duration_days > max_duration_days

  drop <- bad_cases | bad_min_size | bad_duration
  n_drop <- sum(drop)

  if (n_drop > 0) {
    warning(
      sprintf(
        paste0(
          "[%s/guard] dropping %d / %d outbreak episode(s): ",
          "%d suspected_cases<=0, %d suspected_cases<%d, %d duration>%dd ",
          "(e.g. location_period_id: %s)"
        ),
        tag,
        n_drop,
        nrow(dat),
        sum(bad_cases),
        sum(bad_min_size),
        min_cases,
        sum(bad_duration),
        max_duration_days,
        paste(
          utils::head(unique(dat$location_period_id[drop]), 10),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  } else {
    message(sprintf(
      "[%s/guard] outbreak-episode plausibility checks: 0 / %d dropped (cases<=0 / cases<%d / duration>%dd)",
      tag,
      nrow(dat),
      min_cases,
      max_duration_days
    ))
  }

  dat[!drop, , drop = FALSE]
}


#' lp_coverage_season_years
#'
#' The season-years a location period's surveillance record covers, on a given
#' season calendar (`min_seas`): one row per (LP, season-year) with that
#' year's coverage window. Shared with add_no_outbreak_years() so the gid-season occurrence builder
#' can compute an LP's coverage on its GID's canonical calendar.
#'
#' @param lp_id one location_period_id
#' @param min_seas the season-calendar start month to read the record on
#' @param all_ts_dat all time series data available
lp_coverage_season_years <- function(lp_id, min_seas, all_ts_dat) {
  # Local names that cannot collide with all_ts_dat columns under data masking.
  .lp_id <- lp_id
  .ms_lp <- min_seas

  ts_dat <- all_ts_dat %>%
    # Get all time series data for this location period id
    filter(location_period_id == .lp_id) %>%
    # Define cholera-season adjusted months
    mutate(
      TL_seas = TL %m+% months(-(.ms_lp - 1)),
      TR_seas = TR %m+% months(-(.ms_lp - 1)),
      yearL_seas = lubridate::year(TL_seas),
      yearR_seas = lubridate::year(TR_seas)
    )

  # Expand to all years covered by the time series data
  all_years <- ts_dat %>%
    group_by(location_period_id, TL_seas, TR_seas) %>%
    # Define all years covered by the data
    group_modify(
      function(.x, .y) {
        tibble(yearL_seas = seq(.x$yearL_seas[1], .x$yearR_seas[1]))
      }
    ) %>%
    group_by(location_period_id, TL_seas, TR_seas) %>%
    # Assign years covered by data
    mutate(
      TL_seas_new0 = lubridate::floor_date(
        `year<-`(TL_seas, yearL_seas),
        "year"
      ),
      TL_seas_new = pmax(TL_seas, TL_seas_new0),
      TR_seas_new0 = lubridate::ceiling_date(
        `year<-`(TR_seas, yearL_seas),
        "year"
      ) -
        1,
      TR_seas_new = pmin(TR_seas, TR_seas_new0),
      TL = TL_seas_new %m+% months((.ms_lp - 1)),
      TR = TR_seas_new %m+% months((.ms_lp - 1)),
      # Also keep the calendar years for plotting
      TL_min_cal = lubridate::floor_date(TL, "month"),
      TR_max_cal = lubridate::ceiling_date(TR, "month") - 1,
      # Also keep the calendar years for plotting
      TL_seas_min_cal = lubridate::floor_date(TL_seas_new, "year"),
      TR_seas_max_cal = lubridate::ceiling_date(TR_seas_new, "year") - 1,
      year_seas_TL = lubridate::year(TL_seas_new)
    ) %>%
    distinct(
      location_period_id,
      TL,
      TR,
      TL_seas = TL_seas_new,
      TR_seas = TR_seas_new,
      year_seas_TL,
      TL_min_cal,
      TR_max_cal,
      TL_seas_min_cal,
      TR_seas_max_cal
    ) %>%
    group_by(year_seas_TL) %>%
    slice_min(TL_seas) %>%
    ungroup()

  all_years
}

#' add_no_outbreak_years
#' Add fake data for all season years without data that are covered by observations
#'
#' @param x this tibble of data to process
#' @param all_ts_dat all time series data available
#'
add_no_outbreak_years <- function(x, all_ts_dat) {
  all_years <- lp_coverage_season_years(
    x$location_period_id[1],
    x$min_seas[1],
    all_ts_dat
  )

  all_years %>%
    # Keep all padded data
    left_join(
      x %>% select(location_period_id, year_seas, TL, TR, TL_seas, TR_seas),
      by = c("location_period_id", "year_seas_TL" = "year_seas"),
      suffix = c(".padded", ".ob")
    ) %>%
    # Set bounds of padded data to seasonal years outside of months with outbreaks
    mutate(
      TL = pmin(TL.padded, TL.ob, na.rm = T),
      TR = case_when(
        is.na(TL.ob) ~ TR.padded,
        T ~ pmax(TL.padded, TL.ob - 1, na.rm = T)
      ),
      TL_seas = pmin(TL_seas.padded, TL_seas.ob, na.rm = T),
      TR_seas = case_when(
        is.na(TL_seas.ob) ~ TR_seas.padded,
        T ~ pmax(TL_seas.padded, TL_seas.ob, na.rm = T)
      )
    ) %>%
    rename(year_seas = year_seas_TL) %>%
    select(-contains(".")) %>%
    select(location_period_id, TL, TR, TL_seas, TR_seas, year_seas)
}


# Levels and the NA level from base.R::.TTLO_LEVELS: the explicit
# labels are exactly the ones cut() derives from these breaks, so the factor is
# identical (tests/test_t2_shared_infra.R).
make_time_to_last_outbreak_cat <- function(df, no_prior_level = FALSE) {
  # The level "none" -- no recorded prior onset of the
  # location -- appended AFTER "[6,Inf]" so "[0,1)" stays the reference. A
  # LOAD-TIME recode (ttlo_no_prior_onset() reads it off any panel); the level
  # contract is shared with the simulator side: label "none",
  # sentinel .TTLO_NONE_DAYS (utils/base.R).
  none <- if (isTRUE(no_prior_level)) {
    ttlo_no_prior_onset(df$t_start, df$time_to_last_outbreak)
  } else {
    NULL
  }
  out <- df %>%
    # Define categories of time to last outbreak
    mutate(
      time_to_last_outbreak = cut(
        (pmax(0, time_to_last_outbreak) + 1) / 365,
        c(0, 1, 2, 3, 4, 5, 6, Inf),
        labels = .TTLO_LEVELS,
        include.lowest = T,
        right = F
      )
    ) %>%
    tidyr::replace_na(list(
      time_to_last_outbreak = .TTLO_LEVELS[length(.TTLO_LEVELS)]
    ))
  if (is.null(none)) {
    return(out)
  }
  lv <- c(.TTLO_LEVELS, .OCC_TTLO_NONE_LABEL)
  out$time_to_last_outbreak <- factor(
    ifelse(none, .OCC_TTLO_NONE_LABEL, as.character(out$time_to_last_outbreak)),
    levels = lv
  )
  out
}

# The TTLO level of "no recorded prior onset" (the simulator side uses the
# same label).
.OCC_TTLO_NONE_LABEL <- "none"

#' active_last_month
#'
#' Active at any time in the calendar month before t_start's month (a span
#' overlap, not the one-day test at t_start - 30, which skipped February for
#' March rows and missed outbreaks inside the previous month). One value per
#' group of span rows; NA spans (no outbreak) count as inactive. Lives here
#' so the contagion-coverage activity grid (spatial_contagion.R::adm1_activity_grid()) reuses it.
active_last_month <- function(TL, TR, t_start) {
  m0 <- lubridate::floor_date(t_start, "month")
  any(coalesce(TL <= m0 - 1 & TR >= (m0 %m-% months(1)), FALSE))
}

#' ttlo_no_prior_onset
#'
#' TRUE where a panel row has no recorded prior onset of its location: the
#' first row of a history (lag() gives NA) or a row whose "last outbreak" is
#' still the TTLO origin (03_prepare_occ.R: .TTLO_ORIGIN, 2000-01-01).
#' @param t_start row start (Date).
#' @param ttlo time_to_last_outbreak in days (numeric, before categorising).
ttlo_no_prior_onset <- function(t_start, ttlo) {
  ttlo <- as.numeric(ttlo)
  is.na(ttlo) |
    (as.numeric(as.Date(t_start)) - ttlo) == as.numeric(.TTLO_ORIGIN)
}

#' mask_occ_panel_rows
#'
#' Masks of the gid_season occurrence panel, applied AFTER the
#' national zero-fill (masking first would let
#' national_zero_surv_rows() re-add the months):
#'   * at_risk_start = "after_episode": a row of gid g inside (TL, TR]
#'     of an EARLIER-season episode of g is not at risk (g is in an outbreak);
#'     an onset there is masked as an occurrence event (the episode stays in
#'     the episode table, so duration / monthly do not change)
#'     (.own_episode_span_hits(rule = "running")).
#'   * adm1_overlap = "censor_orphan": the rows of the ADM2 children of
#'     an ORPHAN ADM1 parent episode (adm1_overlap_parents()) inside its month
#'     span; by construction never an onset (asserted).
#' At the default spec (season, drop) `surv` is returned unchanged (identical,
#' no attributes added).
#'
#' @param surv panel rows (gid, t_start, t_end, outbreak[, national_zero_fill]).
#' @param ep canonical episodes (gid, year_seas, TL, TR, min_seas).
#' @param spec occurrence-panel spec (utils/paths.R::occ_panel_spec_of()).
#' @param overlap adm1_overlap_parents() output (needed for censor_orphan).
#' @return the kept rows; attributes n_rows_masked_running_episode,
#'   n_onsets_masked_running_episode, n_nz_rows_masked,
#'   n_child_rows_censored_adm1 and `masked_rows` (gid, t_start, outbreak,
#'   national_zero_fill, reason) when a mask is on.
mask_occ_panel_rows <- function(surv, ep, spec, overlap = NULL) {
  run_on <- identical(spec$at_risk_start, "after_episode")
  a1_on <- identical(spec$adm1_overlap, "censor_orphan")
  if (!run_on && !a1_on) {
    return(surv)
  }
  nzf <- if ("national_zero_fill" %in% names(surv)) {
    surv$national_zero_fill %in% TRUE
  } else {
    rep(FALSE, nrow(surv))
  }
  ev <- as.logical(surv$outbreak) %in% TRUE
  run_hit <- if (run_on) {
    .own_episode_span_hits(surv, ep, rule = "running")
  } else {
    rep(FALSE, nrow(surv))
  }
  a1_hit <- rep(FALSE, nrow(surv))
  if (a1_on) {
    if (is.null(overlap)) {
      stop(
        "mask_occ_panel_rows: adm1_overlap = censor_orphan needs adm1_overlap_parents().",
        call. = FALSE
      )
    }
    orph <- overlap[overlap$class == "orphan", , drop = FALSE]
    if (nrow(orph)) {
      m <- function(d) lubridate::floor_date(as.Date(d), "month")
      rows <- tibble::tibble(
        .i = seq_len(nrow(surv)),
        gid = as.character(surv$gid),
        month = m(surv$t_start),
        adm1 = gid_to_adm1(surv$gid, "mask_occ_panel_rows"),
        child = as.character(surv$admin_level) == "ADM2"
      )
      j <- rows[rows$child, , drop = FALSE] %>%
        dplyr::inner_join(
          dplyr::transmute(orph, adm1 = ADM1, p_TL = m(TL), p_TR = m(TR)),
          by = "adm1",
          relationship = "many-to-many"
        ) %>%
        dplyr::filter(month >= p_TL, month <= p_TR)
      a1_hit[unique(j$.i)] <- TRUE
    }
    a1_hit <- a1_hit & !run_hit # counted once, under the running-episode mask first
    if (any(a1_hit & ev)) {
      stop(
        "mask_occ_panel_rows: ",
        sum(a1_hit & ev),
        " onset row(s) of an ADM2 child fall inside an ",
        "orphan ADM1 parent episode -- an orphan has no overlapping child episode by definition.",
        call. = FALSE
      )
    }
  }
  drop <- run_hit | a1_hit
  out <- surv[!drop, , drop = FALSE]
  attr(out, "n_rows_masked_running_episode") <- sum(run_hit)
  attr(out, "n_onsets_masked_running_episode") <- sum(run_hit & ev)
  attr(out, "n_nz_rows_masked") <- sum(drop & nzf)
  attr(out, "n_child_rows_censored_adm1") <- sum(a1_hit)
  attr(out, "masked_rows") <- tibble::tibble(
    gid = as.character(surv$gid[drop]),
    t_start = as.Date(surv$t_start[drop]),
    outbreak = ev[drop],
    national_zero_fill = nzf[drop],
    reason = ifelse(run_hit[drop], "running_episode", "adm1_orphan_child")
  )
  out
}

scale_precip <- function(df) {
  df %>%
    group_by(uid) %>%
    mutate(
      # Precipitation scale to by 100mm
      across(
        any_of(
          colnames(df) %>%
            stringr::str_subset("precip|cum7") %>%
            stringr::str_subset("hlq", negate = T) %>%
            # standardised anomalies are unitless: never rescaled
            stringr::str_subset("^dev_sd_", negate = T)
        ),
        ~ . / 100
      )
    ) %>%
    ungroup()
}

scale_temp <- function(df) {
  df %>%
    group_by(uid) %>%
    mutate(
      # Temperature scale to by 10C
      across(
        any_of(
          colnames(df) %>%
            stringr::str_subset("temp") %>%
            stringr::str_subset("hlq", negate = T) %>%
            stringr::str_subset("dev", negate = T)
        ),
        ~ . / 10
      ),
      # Temperature deviations scale to by 1C
      across(
        any_of(
          colnames(df) %>%
            stringr::str_subset("temp") %>%
            stringr::str_subset("hlq", negate = T) %>%
            stringr::str_subset("dev")
        ),
        ~ . / 1
      )
    ) %>%
    ungroup()
}

make_factors_anomalies <- function(df) {
  df %>%
    group_by(uid) %>%
    mutate(
      # Factors for anomalies
      across(contains("hlq"), ~ . == "high")
    ) %>%
    ungroup()
}

#' make_climate_lags
#'
#' Add one-month-lagged climate columns (<col>_lag1) to a survival dataset.
#' Lags are computed within each uid × year_seas group (one cholera season-year)
#' so season boundaries produce NAs rather than spurious carry-over values.
#' The resulting _lag* columns are stored at the raw (unscaled) level; scaling
#' by scale_precip() / scale_temp() at load time (load_surv_data / read_size_data)
#' picks them up automatically via colnames pattern matching.
#'
#' Affected columns: all ^(mean|dev)_ columns except dev_sd_* and existing _lag\d+.
#'
#' @param df       the survival / size dataset
#' @param max_lag  maximum lag to generate (default 1L for backward-compat)
#' @param mode     "panel" (default, occurrence-side season-aware dplyr::lag) or
#'                 "event" (size-side SQLite lookup, one row per outbreak)
#' @param baseline required for mode = "event" (e.g. "1990_2019")
#' @return df with additional <col>_lag1 … <col>_lag{max_lag} columns appended
make_climate_lags <- function(
  df,
  max_lag = 1L,
  mode = c("panel", "event"),
  baseline = NULL
) {
  max_lag <- as.integer(max_lag)
  mode <- match.arg(mode)

  if (mode == "event") {
    # ---- Event mode: outbreak-size dataset (one row per outbreak) ----------
    # Uses gid + TL (outbreak start date) to look up pre-existing monthly
    # climate stats from SQLite and vectorises left_join for each lag k.
    # Boundary months with no DB record become NA → replaced by 0 below.
    if (is.null(baseline)) {
      stop(
        "make_climate_lags: mode='event' requires a non-NULL baseline ",
        "(e.g. '1990_2019')"
      )
    }
    if (!all(c("gid", "TL") %in% colnames(df))) {
      stop(
        "make_climate_lags: mode='event' requires 'gid' and 'TL' columns in df"
      )
    }

    clim_lookup <- .make_size_clim_stats(unique(df$gid), baseline)
    val_cols <- setdiff(colnames(clim_lookup), c("gid", "month"))
    tl_months <- lubridate::floor_date(as.Date(df$TL), "month") # onset month, as make_size_data()

    for (k in seq(1L, max_lag)) {
      target_months <- lubridate::add_with_rollback(tl_months, months(-k))
      lag_join_df <- dplyr::tibble(gid = df$gid, month = target_months) %>%
        dplyr::left_join(clim_lookup, by = c("gid", "month"))
      new_cols <- as.data.frame(lag_join_df[, val_cols, drop = FALSE])
      colnames(new_cols) <- paste0(val_cols, "_lag", k)
      new_cols[is.na(new_cols)] <- 0
      df <- dplyr::bind_cols(df, dplyr::as_tibble(new_cols))
    }
    return(df)
  }

  # ---- Panel mode (default, backward-compatible): occurrence dataset ------
  # Season-aware dplyr::lag() within (uid × season-year) groups.
  # Generalised to k = 1..max_lag; backward-compat: max_lag = 1 adds only _lag1.
  clim_cols <- colnames(df) %>%
    stringr::str_subset("^(mean|dev)_") %>%
    stringr::str_subset("^dev_sd_", negate = TRUE) %>%
    stringr::str_subset("_lag\\d+$", negate = TRUE)

  # Derive season-year grouping from t_start_seas resets (a decrease in
  # t_start_seas signals a new cholera season-year within the same uid).
  # Data must be sorted by uid, t_start beforehand — 03_prepare_occ.R's
  # arrange(location_period_id, t_start) already guarantees this.
  df_out <- df %>%
    dplyr::group_by(uid) %>%
    dplyr::mutate(
      # <= 0, not < 0: a one-row season at t_start_seas 0 followed by the next
      # season's first row (also 0) is a reset too. Within a season-year
      # t_start_seas strictly increases.
      .season_grp = cumsum(c(0L, diff(as.integer(t_start_seas)) <= 0L))
    ) %>%
    dplyr::group_by(uid, .season_grp)

  for (k in seq(1L, max_lag)) {
    lag_k_name <- paste0("{.col}_lag", k)
    df_out <- dplyr::mutate(
      df_out,
      dplyr::across(
        dplyr::all_of(clim_cols),
        ~ dplyr::lag(., k),
        .names = lag_k_name
      )
    )
  }

  df_out %>%
    dplyr::ungroup() %>%
    dplyr::select(-.season_grp)
}

# Data checks -------------------------------------------------------------

check_outbreak_data <- function(dat, ref) {
  testthat::expect_setequal(unique(dat$uid), unique(ref$uid))
}


# ---------------------------------------------------------------------------
# Gid-season occurrence builder (`episodes.key: gid_season`)
# ---------------------------------------------------------------------------
# One survival history per (gid, season-year), built from the canonical
# episode table (utils/footprint.R::build_canonical_episodes()) instead of per
# location period: ONE onset per gid-season at the earliest member onset, the
# gid censored after it, and the at-risk window of a no-episode year = the
# UNION of the gid's member LPs' surveillance coverage, read on the gid's
# canonical season calendar. make_surv_data() itself is reused unchanged.

#' gid_at_risk_years
#'
#' @param lp_gid tibble(location_period_id, gid, min_seas) -- one row per
#'   (LP, gid) pair, `min_seas` the gid's CANONICAL calendar.
#' @param all_ts_dat all time series data available
#' @return tibble(gid, year_seas, TR_seas): per gid-season, the latest coverage
#'   end over the gid's member LPs.
gid_at_risk_years <- function(lp_gid, all_ts_dat) {
  lp_gid <- dplyr::distinct(lp_gid, location_period_id, gid, min_seas)
  cov <- purrr::pmap_dfr(lp_gid, function(location_period_id, gid, min_seas) {
    lp_coverage_season_years(location_period_id, min_seas, all_ts_dat) %>%
      dplyr::transmute(gid = gid, year_seas = year_seas_TL, TR_seas = TR_seas)
  })
  if (!nrow(cov)) {
    return(tibble::tibble(
      gid = character(),
      year_seas = integer(),
      TR_seas = as.Date(character())
    ))
  }
  cov %>%
    dplyr::group_by(gid, year_seas) %>%
    dplyr::summarise(TR_seas = max(TR_seas), .groups = "drop")
}

#' gid_season_surv_input
#'
#' The `x` rows make_surv_data() is mapped over, one per (gid, year_seas):
#' episode years end at the canonical onset (`.has_onset = TRUE`), covered
#' years without an episode end at the union coverage end (`.has_onset =
#' FALSE`, the caller sets `outbreak = FALSE` on those rows exactly as the
#' lp_season path does for its zero-padded years). Episode years are always
#' present, covered or not, as before.
#'
#' @param episodes canonical episodes (gid, year_seas, uid, location_period_id,
#'   min_seas, TL_seas).
#' @param at_risk gid_at_risk_years() output.
gid_season_surv_input <- function(episodes, at_risk) {
  ep <- episodes %>%
    dplyr::transmute(
      gid,
      year_seas = as.integer(year_seas),
      uid,
      location_period_id,
      min_seas,
      TL_seas = as.Date(TL_seas),
      .has_onset = TRUE
    )
  ids <- dplyr::distinct(ep, gid, uid, location_period_id, min_seas)
  if (anyDuplicated(ids$gid)) {
    stop(
      "gid_season_surv_input: a gid carries more than one representative ",
      "uid/location_period_id/min_seas.",
      call. = FALSE
    )
  }
  pad <- at_risk %>%
    dplyr::mutate(year_seas = as.integer(year_seas)) %>%
    dplyr::anti_join(ep, by = c("gid", "year_seas")) %>%
    dplyr::inner_join(ids, by = "gid") %>%
    dplyr::transmute(
      gid,
      year_seas,
      uid,
      location_period_id,
      min_seas,
      TL_seas = as.Date(TR_seas),
      .has_onset = FALSE
    )
  dplyr::bind_rows(ep, pad) %>% dplyr::arrange(gid, year_seas)
}

#' national_zero_surv_rows
#'
#' National zero-fill (`episodes.national_zero_fill`, gid_season only): the
#' no-onset rows to add to the gid-season occurrence panel for every
#' national-zero month (01_preprocess_national_zero_months.R) of a panel gid's
#' country that is not already a panel row. A month at or after an onset in its
#' season-year is never added (the gid is censored after its onset). A new
#' season-year may start mid-season (left-truncated): each row is an exposure
#' interval at its own season month.
#'
#' Rows are built by make_surv_data() over the whole season-year (TL_seas = the
#' season-year's last day) and then filtered to the national-zero months, so
#' they carry exactly the columns, t_start_seas / t_end_seas and climate joins
#' of every other panel row.
#'
#' @param surv_dat the gid-season panel built so far (gid, t_start, outbreak).
#' @param nz_months national-zero months (country, month).
#' @param ids one row per gid: gid, uid, location_period_id, min_seas, country.
#' @param baseline climate baseline, passed to `build`.
#' @param build,map_fn make_surv_data() / purrr::map_dfr by default (the prep
#'   passes furrr::future_map_dfr; tests pass a stub builder).
#' @param episodes canonical episodes (gid, TL, TR) or NULL (legacy). When
#'   given (a named national-zero rule), a candidate month that
#'   lies inside one of its gid's own episode spans -- typically the previous
#'   season's episode still running -- is not added either
#'   (.own_episode_span_hits(rule = "span")); attribute n_in_episode counts them.
#' @return the new rows with outbreak = FALSE and national_zero_fill = TRUE.
#'   Attributes: n_candidates, n_after_onset (dropped), n_season_years,
#'   n_new_season_years (no previous panel row in that season-year), and
#'   n_in_episode only when `episodes` is given.
national_zero_surv_rows <- function(
  surv_dat,
  nz_months,
  ids,
  baseline,
  build = make_surv_data,
  map_fn = purrr::map_dfr,
  episodes = NULL
) {
  ids <- dplyr::distinct(ids, gid, uid, location_period_id, min_seas, country)
  if (anyDuplicated(ids$gid)) {
    stop(
      "national_zero_surv_rows: a gid carries more than one representative ",
      "uid/location_period_id/min_seas/country.",
      call. = FALSE
    )
  }
  seas_of <- function(month, ms) {
    as.integer(lubridate::year(month %m+% months(-(ms - 1L))))
  }

  have <- surv_dat %>%
    dplyr::transmute(
      gid,
      month = lubridate::floor_date(as.Date(t_start), "month")
    ) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(dplyr::select(ids, gid, min_seas), by = "gid") %>%
    dplyr::mutate(year_seas = seas_of(month, min_seas))
  onset <- surv_dat %>%
    dplyr::filter(as.logical(outbreak)) %>%
    dplyr::transmute(
      gid,
      onset_month = lubridate::floor_date(as.Date(t_start), "month")
    ) %>%
    dplyr::inner_join(dplyr::select(ids, gid, min_seas), by = "gid") %>%
    dplyr::mutate(year_seas = seas_of(onset_month, min_seas)) %>%
    dplyr::group_by(gid, year_seas) %>%
    dplyr::summarise(onset_month = min(onset_month), .groups = "drop")

  cand <- ids %>%
    dplyr::inner_join(
      nz_months %>%
        dplyr::transmute(
          country = as.character(country),
          month = as.Date(month)
        ),
      by = "country",
      relationship = "many-to-many"
    ) %>%
    dplyr::anti_join(have, by = c("gid", "month")) %>%
    dplyr::mutate(year_seas = seas_of(month, min_seas)) %>%
    dplyr::left_join(onset, by = c("gid", "year_seas"))
  n_candidates <- nrow(cand)
  after_onset <- !is.na(cand$onset_month) & cand$month >= cand$onset_month
  n_after_onset <- sum(after_onset)
  cand <- cand[!after_onset, ]
  if (!is.null(episodes)) {
    in_ep <- .own_episode_span_hits(
      dplyr::transmute(
        cand,
        gid,
        t_start = month,
        t_end = lubridate::ceiling_date(month, "month") - 1
      ),
      episodes,
      rule = "span"
    )
    n_in_episode <- sum(in_ep)
    cand <- cand[!in_ep, ]
  }

  seasons <- cand %>%
    dplyr::distinct(gid, uid, location_period_id, min_seas, year_seas)
  rows <- if (!nrow(seasons)) {
    NULL
  } else {
    map_fn(split(seasons, seq_len(nrow(seasons))), function(s) {
      x <- s %>%
        dplyr::mutate(TL_seas = as.Date(sprintf("%d-12-31", year_seas)))
      build(x, baseline = baseline)
    })
  }
  if (is.null(rows) || !nrow(rows)) {
    out <- surv_dat[0, , drop = FALSE]
    out$national_zero_fill <- logical()
  } else {
    out <- rows %>%
      dplyr::mutate(
        .month = lubridate::floor_date(as.Date(t_start), "month")
      ) %>%
      dplyr::semi_join(
        cand %>% dplyr::select(gid, .month = month),
        by = c("gid", ".month")
      ) %>%
      dplyr::select(-.month) %>%
      dplyr::mutate(outbreak = FALSE, national_zero_fill = TRUE)
  }
  attr(out, "n_candidates") <- n_candidates
  attr(out, "n_after_onset") <- n_after_onset
  attr(out, "n_season_years") <- nrow(dplyr::distinct(seasons, gid, year_seas))
  attr(out, "n_new_season_years") <- nrow(dplyr::anti_join(
    dplyr::distinct(seasons, gid, year_seas),
    have,
    by = c("gid", "year_seas")
  ))
  if (!is.null(episodes)) {
    attr(out, "n_in_episode") <- n_in_episode
  }
  out
}

#' .own_episode_span_hits
#'
#' Which panel rows lie inside a canonical episode of their OWN gid -- the one
#' predicate the panel masks and the national-zero rule share.
#'
#' @param surv rows with gid, t_start, t_end (Dates).
#' @param ep canonical episodes with gid, TL, TR (Dates); `rule = "running"`
#'   also needs year_seas and min_seas (the gid's canonical calendar).
#' @param rule "span": the row's MONTH lies in [month(TL), month(TR)] of an own
#'   episode (the national-zero gate: a no-onset month inside an outbreak
#'   the database recorded). "running": the row overlaps (TL, TR] of an own
#'   episode of an EARLIER season-year -- `t_start <= TR` and `t_end >= TL + 1`
#'   -- i.e. the gid is still in an outbreak that began in a previous season
#'   (at_risk_start = after_episode); the episode's own onset row
#'   (t_end = TL) is never hit.
#' @return logical, one per row of `surv`.
.own_episode_span_hits <- function(surv, ep, rule = c("span", "running")) {
  rule <- match.arg(rule)
  n <- nrow(surv)
  if (!n || is.null(ep) || !nrow(ep)) {
    return(rep(FALSE, n))
  }
  s <- tibble::tibble(
    .i = seq_len(n),
    gid = as.character(surv$gid),
    t_start = as.Date(surv$t_start),
    t_end = as.Date(surv$t_end)
  )
  e <- tibble::tibble(
    gid = as.character(ep$gid),
    e_TL = as.Date(ep$TL),
    e_TR = as.Date(ep$TR)
  )
  if (rule == "running") {
    e$e_year <- as.integer(ep$year_seas)
    ms <- dplyr::distinct(tibble::tibble(
      gid = as.character(ep$gid),
      min_seas = as.integer(ep$min_seas)
    ))
    if (anyDuplicated(ms$gid)) {
      stop(
        ".own_episode_span_hits: a gid carries more than one min_seas.",
        call. = FALSE
      )
    }
    s <- dplyr::left_join(s, ms, by = "gid") %>%
      dplyr::mutate(
        row_year = as.integer(lubridate::year(
          t_start %m+% months(-(min_seas - 1L))
        ))
      )
  }
  j <- dplyr::inner_join(s, e, by = "gid", relationship = "many-to-many")
  hit <- if (rule == "span") {
    m <- lubridate::floor_date(j$t_start, "month")
    m >= lubridate::floor_date(j$e_TL, "month") &
      m <= lubridate::floor_date(j$e_TR, "month")
  } else {
    !is.na(j$row_year) &
      j$e_year < j$row_year &
      j$t_start <= j$e_TR &
      j$t_end >= j$e_TL + 1
  }
  out <- rep(FALSE, n)
  out[unique(j$.i[hit %in% TRUE])] <- TRUE
  out
}

#' assert_nz_rows_outside_episodes
#'
#' Post-mask gate: no national zero-fill
#' row lies inside its own gid's canonical episode span. Asserted under a named
#' national-zero rule (the fill then excludes those months,
#' national_zero_surv_rows(episodes =)); report-only under the legacy rule, whose
#' fill does not exclude a previous season's running episode.
#'
#' @param surv the panel (gid, t_start, t_end, national_zero_fill).
#' @param ep canonical episodes (gid, TL, TR).
#' @param assert TRUE = stop on any hit; FALSE = message only.
#' @return invisibly, the number of national zero-fill rows inside an own episode.
assert_nz_rows_outside_episodes <- function(surv, ep, assert = TRUE) {
  nz <- surv[surv$national_zero_fill %in% TRUE, , drop = FALSE]
  hits <- .own_episode_span_hits(nz, ep, rule = "span")
  n <- sum(hits)
  if (n > 0L) {
    msg <- sprintf(
      "[occ/national-zero] %d national zero-fill row(s) lie inside their gid's own episode span (%s)",
      n,
      paste(
        utils::head(
          unique(paste(nz$gid[hits], format(nz$t_start[hits], "%Y-%m"))),
          5
        ),
        collapse = ", "
      )
    )
    if (isTRUE(assert)) {
      stop(msg, ".", call. = FALSE)
    }
    message(msg, " -- report only under the legacy rule.")
  }
  invisible(n)
}
