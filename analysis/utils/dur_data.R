# dur_data.R — duration (cessation-hazard) track data preparation: one row per
#   ACTIVE MONTH of an episode.
#
# Provides:
#   make_duration_panel()      episode table -> active-month panel
#   bin_mso()                  months-since-onset bands (mso_band)
#   derive_series_end(), flag_right_censored()   censoring of episodes
#   duration_span_summary()    onset-to-last-active span (stamped by 03)
#   outbreak_episode_guard_flags()   what the occurrence guard WOULD drop
#   duration_episode_rows(), load_duration_data()   read the panel for a fit
#
# Details:
#   Poisson with log offset, like occurrence:
#     y          = as.integer(ended)   1 on the last active month iff uncensored
#     exposure   = days of calendar month m intersected with [TL, TR]
#     log_offset = log(pmax(exposure, 1))
#   Episodes come from the size track's episode table, not a rebuild (see
#   make_duration_panel()).

# ---------------------------------------------------------------------------
# months-since-onset banding
# ---------------------------------------------------------------------------

#' Band edges for the `mso` (months-since-onset) baseline hazard.
#'
#' Left-closed, right-open on `mso` (which is 1-based: the onset month is 1).
#'
#' These edges were FITTED against the realised panel, not assumed, and kept
#' unchanged for the outbreak_data_Nov_2026 vintage, whose outbreaks end at the
#' last week before their six-week closing run. Empirical monthly cessation
#' hazard on that panel (outbreaks_duration_data_base-1990_2019_ep-gs.rds, the
#' deployed prep: 7137 active months, 1404 episodes, 1300 observed cessations):
#'
#'   band   at-risk  events  haz/month
#'   1        1404      44     0.0313
#'   2        1352     239     0.1768
#'   3        1096     257     0.2345
#'   4-5      1439     333     0.2314
#'   6-8      1084     258     0.2380
#'   9+        762     169     0.2218
#'
#' Two things this says. First, the hazard is NOT flat -- a 7.6x range across
#' bands -- so a constant-hazard (geometric) duration model would be wrong and
#' `basehaz_duration: none` is not a serious default. Second, the structure is
#' entirely in the first two months: 3, 4-5, 6-8 and 9+ sit within noise of
#' each other. The tail bands are kept separate anyway (the model can shrink
#' them together; merging them cannot be undone at fit time), but a
#' `c(1, 2, 3, 4, Inf)` variant is the natural parsimony check and is reachable
#' via `duration.mso_breaks` in config without re-running 03_prepare_dur.R.
#'
#' The mso = 1 band is informed by data: the minimum outbreak length is two
#' weeks (the shortest episode runs 13 days from TL to TR), 52 episodes are
#' active in a single calendar month and 44 of them end in it. Before the
#' Nov_2026 vintage the closing run was part of the outbreak, no episode ended
#' in its first month (0 of 1405 on that vintage's panel) and this band's
#' coefficient was prior-dominated.
.DUR_MSO_BREAKS <- c(1, 2, 3, 4, 6, 9, Inf)

#' bin_mso
#'
#' @param mso integer vector, months since onset (1 = onset month).
#' @param breaks numeric, left edges plus a terminal Inf.
#' @return ordered factor with labels "1", "2", "3", "4-5", "6-8", "9+".
bin_mso <- function(mso, breaks = .DUR_MSO_BREAKS) {
  if (length(breaks) < 2L || !is.infinite(breaks[length(breaks)])) {
    stop("bin_mso: `breaks` must end in Inf.", call. = FALSE)
  }

  labs <- character(length(breaks) - 1L)
  for (i in seq_len(length(breaks) - 1L)) {
    lo <- breaks[i]
    hi <- breaks[i + 1L]
    labs[i] <- if (is.infinite(hi)) {
      sprintf("%d+", as.integer(lo))
    } else if (hi - lo == 1) {
      sprintf("%d", as.integer(lo))
    } else {
      sprintf("%d-%d", as.integer(lo), as.integer(hi) - 1L)
    }
  }

  cut(
    as.numeric(mso),
    breaks = breaks,
    labels = labs,
    right = FALSE,
    include.lowest = TRUE,
    ordered_result = FALSE
  )
}


# ---------------------------------------------------------------------------
# right-censoring
# ---------------------------------------------------------------------------

#' derive_series_end
#'
#' Last date on which each location's surveillance record carries data, from
#' the pre-outbreak-extraction aggregated time series written by
#' 01_preprocess_outbreaks.R:107. An episode whose TR sits at that boundary was
#' still ongoing when observation stopped, i.e. right-censored.
#'
#' KEYED ON `location_period_id`, NOT `gid`. build_coverage_mask() (sim_io.R),
#' which bridges location_period_id -> gid, is lossy for this purpose:
#'   * build_coverage_mask() filters the series to admin_level ADM2/ADM3
#'     (sim_io.R:520), discarding the 32 ADM0 and 252 ADM1 series. 31 of 1582
#'     episodes then had no coverage window at all.
#'   * Worse, it silently substitutes a country x ADM2 window for gids with no
#'     series of their own, which put NGA's series end at 2020-01 while NGA
#'     episodes run to 2022-08 -- 30 spurious "censored" flags, none real.
#'   * Keying directly on location_period_id matches all 1582 episodes with no
#'     fallback, and leaves only 2 episodes with TR beyond their own series end.
#' The gid bridge exists in sim_io.R because the SIMULATOR works in gid space.
#' The censoring question is asked of the observation record, which is in LP
#' space, so no bridge is needed here.
#'
#' @param episodes data frame carrying `location_period_id` and `country`.
#' @param aggregated_ts_path path to
#'   ~/secure/cholera_outbreak_data/aggregated_ts_preoutbreak_extraction.rds.
#' @param aggregated_adm_ts_path path to
#'   ~/secure/cholera_outbreak_data/aggregated_admin_level_ts_preoutbreak_extraction.rds. Used
#'   only for LPs absent from the LP-level series -- empty in practice on the
#'   current corpus, kept so a future corpus refresh degrades rather than fails.
#' @return tibble(location_period_id, series_end [Date], series_end_source).
derive_series_end <- function(
  episodes,
  aggregated_ts_path,
  aggregated_adm_ts_path = NULL
) {
  if (!file.exists(aggregated_ts_path)) {
    stop(
      "derive_series_end: aggregated_ts not found: ",
      aggregated_ts_path,
      call. = FALSE
    )
  }

  need <- c("location_period_id", "country")
  miss <- setdiff(need, names(episodes))
  if (length(miss)) {
    stop(
      "derive_series_end: `episodes` is missing: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  lp_end <- readRDS(aggregated_ts_path) %>%
    dplyr::mutate(
      location_period_id = as.character(location_period_id),
      TR = as.Date(TR)
    ) %>%
    dplyr::filter(!is.na(TR)) %>%
    dplyr::group_by(location_period_id) %>%
    dplyr::summarise(series_end = max(TR), .groups = "drop") %>%
    dplyr::mutate(series_end_source = "location_period")

  out <- episodes %>%
    dplyr::distinct(
      location_period_id = as.character(location_period_id),
      country = as.character(country)
    ) %>%
    dplyr::left_join(lp_end, by = "location_period_id")

  # Fallback: country x admin_level window, for LPs with no series of their own.
  missing_lp <- out$location_period_id[is.na(out$series_end)]
  if (
    length(missing_lp) > 0 &&
      !is.null(aggregated_adm_ts_path) &&
      file.exists(aggregated_adm_ts_path)
  ) {
    ctry_end <- readRDS(aggregated_adm_ts_path) %>%
      dplyr::mutate(country = as.character(country), TR = as.Date(TR)) %>%
      dplyr::filter(!is.na(TR)) %>%
      dplyr::group_by(country) %>%
      dplyr::summarise(ctry_end = max(TR), .groups = "drop")

    out <- out %>%
      dplyr::left_join(ctry_end, by = "country") %>%
      dplyr::mutate(
        series_end_source = dplyr::if_else(
          is.na(series_end) & !is.na(ctry_end),
          "country_fallback",
          series_end_source
        ),
        series_end = dplyr::coalesce(series_end, ctry_end)
      ) %>%
      dplyr::select(-ctry_end)

    message(sprintf(
      "[dur/censor] %d location period(s) used the country-level series-end fallback.",
      length(missing_lp)
    ))
  } else if (length(missing_lp) > 0) {
    message(sprintf(
      "[dur/censor] %d location period(s) have no series end -- treated as uncensored.",
      length(missing_lp)
    ))
  }

  out %>%
    dplyr::mutate(
      series_end_source = dplyr::coalesce(series_end_source, "none")
    ) %>%
    dplyr::select(location_period_id, series_end, series_end_source)
}

#' flag_right_censored
#'
#' An episode is right-censored when its recorded end sits at (or beyond) the
#' end of its location's observation record: we saw it stop being reported, not
#' stop happening. TZA is the documented reference case -- its database ends in
#' 2019 -- and 03_prepare_dur.R hard-stops if TZA comes back unflagged.
#'
#' `censored = TRUE` means the cessation event was NOT observed, so the episode
#' contributes at-risk months with `ended = 0` throughout.
#'
#' @param episodes data frame with `location_period_id`, `country`, `TR`.
#' @param series_end tibble from derive_series_end().
#' @param censor_tol_days integer. TR within this many days of the series end
#'   counts as censored. Default 73 = 31 + 42: the corpus is (sub-)weekly and
#'   an outbreak whose last week falls in the final month of the record cannot
#'   be distinguished from one that was still running; and since the
#'   outbreak_data_Nov_2026 vintage an outbreak's TR is the last week before its
#'   six-week closing run (OutbreakExtractR `exclude_closing_tail`), so an
#'   outbreak still running when the record stops ends 42 days before it.
#' @return `episodes` plus series_end, series_end_source, censored.
flag_right_censored <- function(episodes, series_end, censor_tol_days = 73L) {
  episodes %>%
    dplyr::mutate(location_period_id = as.character(location_period_id)) %>%
    dplyr::left_join(series_end, by = "location_period_id") %>%
    dplyr::mutate(
      censored = !is.na(series_end) &
        as.Date(TR) >= (series_end - as.integer(censor_tol_days))
    )
}


# ---------------------------------------------------------------------------
# the active-month panel
# ---------------------------------------------------------------------------

#' .expand_active_months
#'
#' One row per calendar month intersected by [TL, TR]. `exposure` is the number
#' of days of that month the episode was actually active, so the first and last
#' months are partial -- this is what the Poisson log-offset makes the hazard
#' per-day rather than per-row, and it is why a 3-day episode straddling a month
#' boundary does not contribute two full at-risk months.
#'
#' @param ep data frame with episode_id, TL, TR (one row per episode).
#' @return tibble(episode_id, month [Date, first of month], mso, exposure,
#'   n_months).
.expand_active_months <- function(ep) {
  tl <- as.Date(ep$TL)
  tr <- as.Date(ep$TR)
  if (any(is.na(tl) | is.na(tr))) {
    stop(".expand_active_months: NA in TL/TR.", call. = FALSE)
  }
  if (any(tr < tl)) {
    stop(
      ".expand_active_months: ",
      sum(tr < tl),
      " episode(s) with TR < TL.",
      call. = FALSE
    )
  }

  m0 <- lubridate::floor_date(tl, "month")
  m1 <- lubridate::floor_date(tr, "month")
  n <- (lubridate::year(m1) - lubridate::year(m0)) *
    12L +
    (lubridate::month(m1) - lubridate::month(m0)) +
    1L

  idx <- rep(seq_len(nrow(ep)), times = n)
  mso <- sequence(n) # 1..n within each episode
  month <- m0[idx] %m+% months(mso - 1L)

  month_end <- lubridate::ceiling_date(month, "month") - 1L
  # +1 so a single-day overlap counts as one day of exposure, not zero.
  exposure <- as.numeric(pmin(tr[idx], month_end) - pmax(tl[idx], month)) + 1

  tibble::tibble(
    episode_id = ep$episode_id[idx],
    month = month,
    mso = as.integer(mso),
    exposure = exposure,
    n_months = n[idx]
  )
}

# Why this track exists. At fit time the occurrence contagion covariates are ON
# for every month an outbreak spans TL..TR (03_prepare_occ.R::
# active_last_month()): a whole-month span (onset month to last active month)
# of median 4 [IQR 3-7] under gid_season (duration_span_summary(), stamped by
# 03_prepare_dur.R as median_span_months). A simulator fed last month's ONSET
# indicator is on for exactly one month -- roughly a quarter of the fit-time
# duty cycle. The missing ingredient is a model of when a location stops being
# active.
#
# Source of the episode table. The panel is expanded from the SIZE track's
# episode table (generated_data/outbreak_size_analysis_data.rds), not rebuilt
# from outbreak_metadata.rds. Two reasons, both load-bearing:
#   1. It already carries gid / ADM1 / min_seas / AFRO_region resolved through
#      the seasonality-match cascade (03_prepare_siz.R sections A-B). Rebuilding
#      that cascade here would duplicate ~250 lines and could silently drift.
#   2. Phase 05 (dependence check) pairs size and duration residuals PER
#      EPISODE. Deriving both from one episode table makes the keys align by
#      construction rather than by a join that might not be 1:1.
#' make_duration_panel
#'
#' Expand an episode table into the active-month panel the cessation hazard is
#' fitted on, joining month-of-exposure climate, ENSO lags and seasonal time.
#'
#' The climate join is the point of the whole design: an AFT or count model of
#' duration can only use covariates at onset, whereas this panel carries the
#' climate of EACH active month, so within-outbreak weather can drive when the
#' outbreak ends.
#'
#' @param episodes episode-level table -- the size track's
#'   outbreak_size_analysis_data.rds, or any frame carrying uid,
#'   location_period_id, gid, country, ADM1, AFRO_region, admin_level,
#'   min_seas, TL, TR, suspected_cases, pop.
#' @param baseline character, climate baseline period, e.g. "1990_2019".
#' @param enso_dat ENSO index table from fetch_enso_index(); NULL skips the join.
#' @param enso_max_lead integer, ENSO lead columns enso_lead1..k to add
#'   (join_enso_lags(max_lead =); 0 = none, the prior behaviour).
#' @param max_lag integer, climate lag columns to build (DLNM use at fit time).
#' @param censor_tol_days passed to flag_right_censored().
#' @param aggregated_ts_path,aggregated_adm_ts_path passed to derive_series_end().
#' @param long_episode_days integer. Episodes longer than this are WARNED about,
#'   never dropped -- see the guard_outbreak_episodes() note in this file's header.
#' @return the active-month panel.
make_duration_panel <- function(
  episodes,
  baseline,
  enso_dat = NULL,
  enso_max_lead = 0L,
  max_lag = 11L,
  censor_tol_days = 73L,
  aggregated_ts_path = make_agg_ts_file(),
  aggregated_adm_ts_path = make_agg_adm_ts_file(),
  long_episode_days = 730L,
  censor_col = "location_period_id"
) {
  need <- c(
    "uid",
    "location_period_id",
    "gid",
    "country",
    "TL",
    "TR",
    "min_seas",
    censor_col
  )
  miss <- setdiff(need, names(episodes))
  if (length(miss)) {
    stop(
      "make_duration_panel: `episodes` is missing: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  # --- 1. episode identity ------------------------------------------------
  # (uid, TL, TR) is unique on the current corpus (1582/1582). Asserted rather
  # than assumed: a duplicate would silently double every at-risk month for the
  # episode and inflate the effective sample size of the hazard.
  ep <- episodes %>%
    dplyr::mutate(TL = as.Date(TL), TR = as.Date(TR)) %>%
    dplyr::mutate(episode_id = sprintf("%s|%s|%s", uid, format(TL), format(TR)))

  if (anyDuplicated(ep$episode_id)) {
    stop(
      "make_duration_panel: (uid, TL, TR) is not unique -- ",
      sum(duplicated(ep$episode_id)),
      " duplicate episode_id(s).",
      call. = FALSE
    )
  }

  # --- 2. censoring -------------------------------------------------------
  # `censor_col`: which location period's record decides censoring. Under
  # `episodes.key: gid_season` the episode's location_period_id is a
  # per-gid REPRESENTATIVE, not the member whose report ends the episode, so
  # 03_prepare_dur.R passes `censor_lp_id` (utils/footprint.R) instead. The
  # default is the identity.
  lp_orig <- ep$location_period_id
  if (!identical(censor_col, "location_period_id")) {
    ep$location_period_id <- as.character(ep[[censor_col]])
  }
  ser <- derive_series_end(ep, aggregated_ts_path, aggregated_adm_ts_path)
  ep <- flag_right_censored(ep, ser, censor_tol_days = censor_tol_days)
  if (!identical(censor_col, "location_period_id")) {
    ep$location_period_id <- as.character(lp_orig)
  }

  message(sprintf(
    "[dur/censor] %d / %d episode(s) right-censored (%.1f%%), tol = %d days.",
    sum(ep$censored),
    nrow(ep),
    100 * mean(ep$censored),
    censor_tol_days
  ))

  # --- 3. long-episode report (warn, never drop) --------------------------
  dur_days <- as.numeric(ep$TR - ep$TL)
  long <- which(dur_days > long_episode_days)
  if (length(long)) {
    warning(
      sprintf(
        "make_duration_panel: %d episode(s) exceed %d days and are KEPT (max %.0f d). %s",
        length(long),
        long_episode_days,
        max(dur_days),
        paste(
          sprintf("%s (%.0f d)", ep$uid[long], dur_days[long]),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }

  # --- 4. expand to active months -----------------------------------------
  panel <- .expand_active_months(ep)

  # `ended` is the response: 1 on the last active month of an OBSERVED
  # cessation. Censored episodes are all-zero, which is precisely the
  # discrete-time-hazard likelihood contribution of a censored subject -- they
  # inform the hazard through their at-risk months and never through an event.
  panel <- panel %>%
    dplyr::left_join(
      ep %>% dplyr::select(episode_id, censored, series_end, series_end_source),
      by = "episode_id"
    ) %>%
    dplyr::mutate(
      ended = as.integer(mso == n_months & !censored),
      # pmax mirrors 04_fit_occ_brms.R:124 -- log(0) would be -Inf and drop the
      # row's contribution entirely.
      log_offset = log(pmax(exposure, 1))
    )

  # --- 5. carry episode-level metadata ------------------------------------
  # `duration_days` / `duration_months` are episode constants carried onto every
  # row so the panel is self-sufficient for Phase 5 and Phase 7 without a
  # re-join to the size table.
  # `pop_footprint`/`pop_gid`/`pop_cover` are the P1 denominator triple: `pop`
  # alone is the LP reporting footprint H, which is what the *size* track used
  # to offset by while the simulator replays at the whole GADM unit Q, inflating
  # every mu by 1/coverage. The monthly track inherits this panel wholesale
  # (03_prepare_mon.R left-joins monthly cases onto `dur_dat`), so unless the
  # triple rides along here mon can never declare `denominator: gid`. All three
  # go through `intersect()`, so a panel built before P1 (or under
  # `--pop_source legacy`, which emits none of them) is unaffected.
  meta_cols <- intersect(
    c(
      "episode_id",
      "uid",
      "location",
      "location_period_id",
      "gid",
      "seas_gid",
      "country",
      "ADM1",
      "AFRO_region",
      "admin_level",
      "min_seas",
      "suspected_cases",
      "pop",
      "pop_footprint",
      "pop_gid",
      "pop_cover",
      "year_seas",
      "time_to_last_outbreak",
      "TL",
      "TR"
    ),
    names(ep)
  )

  panel <- panel %>%
    dplyr::left_join(
      ep %>% dplyr::select(dplyr::all_of(meta_cols)),
      by = "episode_id"
    ) %>%
    dplyr::mutate(
      duration_days = as.numeric(TR - TL),
      duration_months = n_months
    )

  # --- 6. seasonal time of the ACTIVE MONTH -------------------------------
  # Same rotation as the size/occurrence tracks (03_prepare_siz.R:286-293), but
  # applied to the active month rather than to TL, so a `basehaz_poisson` on
  # t_start_seas tests "outbreaks end when the season ends" rather than
  # re-expressing the onset date.
  panel <- panel %>%
    dplyr::mutate(
      t_start = month,
      t_end = lubridate::ceiling_date(month, "month") - 1L,
      t_start_seas = as.Date(apply_offset_months(month, -(min_seas - 1L))),
      month_seas = lubridate::month(t_start_seas),
      mso_band = bin_mso(mso)
    )

  # --- 7. climate at the active month -------------------------------------
  # .make_size_clim_stats() (siz_data.R:42) is the batched multi-gid form; the
  # per-gid make_size_data()/make_surv_data() helpers would issue one SQLite query
  # per episode.
  #
  # Raw values (`precip`, `temp_daymean`, ...) are joined ALONGSIDE the derived
  # mean_*/dev_* columns, not instead of them. Two reasons: guard_missing_climate_data()
  # (occ_data.R:112) tests `precip` and `temp_daymean` by name, and the other two
  # tracks both carry raw and derived, so a formula term means the same thing in
  # all three. The hlq_* extreme categories are NOT carried -- they belong to the
  # extreme-index sensitivity variants, which lost the CPO/WAIC comparison
  # (CLAUDE.md, "Preferred climate covariates"); add them here if that is revisited.
  gids <- unique(as.character(panel$gid))
  clim_raw <- get_stats_gid(gids) %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = lubridate::floor_date(as.Date(month), "month")
    )
  clim_der <- .make_size_clim_stats(gids, baseline = baseline) %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = lubridate::floor_date(as.Date(month), "month")
    )
  clim <- dplyr::full_join(clim_raw, clim_der, by = c("gid", "month"))

  n_before <- nrow(panel)
  panel <- panel %>%
    dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::left_join(clim, by = c("gid", "month"))
  stopifnot(nrow(panel) == n_before) # climate must be 1 row per (gid, month)

  # --- 8. climate lags ----------------------------------------------------
  # make_climate_lags(mode = "panel") lags within `uid` x season group and works
  # on any monthly panel (occ_data.R:322-387) -- no new mode is needed.
  #
  # It is handed `episode_id` AS `uid` on purpose. Its season grouping keys on a
  # DECREASE in t_start_seas, which does not necessarily occur between two
  # consecutive episodes at the same location (episode A ending in season month
  # 10 followed by episode B starting in season month 11 shows no decrease), and
  # the lag would then reach across the gap into a different outbreak. Grouping
  # on the episode makes "climate one month ago" mean "one month ago WITHIN THIS
  # OUTBREAK", which is the only reading that makes sense for a cessation
  # hazard. Rows at mso <= k keep NA lags, as in the occurrence track.
  if (max_lag > 0L) {
    panel <- panel %>%
      dplyr::arrange(episode_id, mso) %>%
      dplyr::mutate(.uid_orig = uid, uid = episode_id)
    panel <- make_climate_lags(panel, max_lag = max_lag, mode = "panel")
    panel <- panel %>%
      dplyr::mutate(uid = .uid_orig) %>%
      dplyr::select(-.uid_orig)
  }

  # --- 9. ENSO at the active month ----------------------------------------
  if (!is.null(enso_dat)) {
    panel <- join_enso_lags(
      panel,
      enso_dat,
      date_col = "t_start",
      max_lead = enso_max_lead
    )
  }

  # --- 10. missing-climate guard ------------------------------------------
  panel <- guard_missing_climate_data(panel)

  dplyr::arrange(panel, episode_id, mso)
}


# ---------------------------------------------------------------------------
# episode-level summaries for the policy stamp
# ---------------------------------------------------------------------------

#' duration_span_summary
#'
#' Episode length summaries stamped by 03_prepare_dur.R. The
#' reference length is the whole-month SPAN: whole calendar months from the
#' onset month to the last active month, `n_months - 1` (an episode inside one
#' month has span 0). It is the unit of the deploy lag floor
#' `deploy.lag_selection.occ_floor_offset` (loo_lags.R::check_occ_floor_offset()
#' compares the two). `n_months` itself counts the calendar months TOUCHED
#' (span + 1), which is what the older stamp field `median_duration_d` holds
#' despite its name. The elapsed length (TR - TL) / 30.44 is given beside it.
#'
#' @param ep_lvl one row per episode with `n_months` (calendar months touched,
#'   the panel's `duration_months`), `TL`, `TR`.
#' @return list(median_span_months, span_months_q25, span_months_q75,
#'   median_elapsed_months), numeric.
duration_span_summary <- function(ep_lvl) {
  miss <- setdiff(c("n_months", "TL", "TR"), names(ep_lvl))
  if (length(miss)) {
    stop(
      "duration_span_summary: `ep_lvl` is missing: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  span <- as.numeric(ep_lvl$n_months) - 1
  elapsed <- as.numeric(as.Date(ep_lvl$TR) - as.Date(ep_lvl$TL)) / 30.44
  q <- stats::quantile(span, c(0.25, 0.75), names = FALSE)
  list(
    median_span_months = stats::median(span),
    span_months_q25 = q[1],
    span_months_q75 = q[2],
    median_elapsed_months = stats::median(elapsed)
  )
}

# Deliberately NOT called: guard_outbreak_episodes() (occ_data.R). It has no
# pipeline call site and its max_duration_days = 730 drop would truncate exactly
# the tail the cessation hazard is meant to estimate. Long episodes are warned
# about, never dropped; outbreak_episode_guard_flags() counts what the guard
# WOULD drop, and 03_prepare_dur.R stamps those counts (n_guard_*).
#' outbreak_episode_guard_flags
#'
#' Pure, per-episode version of the plausibility checks of
#' occ_data.R::guard_outbreak_episodes(), with the same thresholds and the same
#' (disjoint) case classes: `cases_le0` = suspected_cases <= 0, `cases_lt50` =
#' 0 < suspected_cases < min_cases, `gt730d` = TR - TL > max_duration_days
#' (days), `any` = the rows the guard would drop. NA values never flag. Nothing
#' is dropped: the duration track deliberately keeps every episode (header), so
#' 03_prepare_dur.R only stamps the counts (n_guard_*). The guard is not
#' rewritten to call this (that would add a cross-module call into occ_data.R);
#' tests/test_duration_span_stamp.R checks instead that `any` selects exactly
#' the rows the guard drops.
#'
#' @param dat episode-level data with `suspected_cases`, `TL`, `TR`.
#' @param min_cases,max_duration_days as guard_outbreak_episodes() (50, 730).
#' @return tibble(cases_le0, cases_lt50, gt730d, any), one row per row of `dat`.
outbreak_episode_guard_flags <- function(
  dat,
  min_cases = 50,
  max_duration_days = 730
) {
  miss <- setdiff(c("suspected_cases", "TL", "TR"), names(dat))
  if (length(miss)) {
    stop(
      "outbreak_episode_guard_flags: `dat` is missing: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  cases <- dat$suspected_cases
  duration_days <- as.numeric(as.Date(dat$TR) - as.Date(dat$TL))
  le0 <- !is.na(cases) & cases <= 0
  lt50 <- !is.na(cases) & cases > 0 & cases < min_cases
  gt730 <- !is.na(duration_days) & duration_days > max_duration_days
  tibble::tibble(
    cases_le0 = le0,
    cases_lt50 = lt50,
    gt730d = gt730,
    any = le0 | lt50 | gt730
  )
}


# ---------------------------------------------------------------------------
# fit-time loader
# ---------------------------------------------------------------------------

#' duration_episode_rows
#'
#' The duration track's row rule from `duration.admin_level`. Under
#' "ADM1-ADM2" the panel goes through apply_admin_level_filter() (gid.R), which
#' drops every coarse ADM1 episode whose ADM1 also carries ADM2 rows: the
#' episode set of load_monthly_data(), the occurrence fit's load-time drop and
#' the simulator frame (drop_overlapping_adm1_from_frame()). Orphan parents leave
#' too; occurrence's `episodes.adm1_overlap: censor_orphan` has nothing to
#' censor here, since an orphan has no overlapping child episode.
#'
#' "ADM1" / "ADM2" keep their historical meaning on this track -- the depth of
#' the intercept RE, no row filter -- and are NOT passed to
#' apply_admin_level_filter(), whose "ADM1" would keep ADM1 rows only.
#'
#' Shared by load_duration_data() and the observed-duration readers (09 R1,
#' 07c), so a raw panel and the fitted panel hold the same episodes.
#'
#' @param d duration panel rows (admin_level, ADM1).
#' @param dur the duration config block (opt$duration); NULL = "ADM1".
#' @return d, possibly with the overlapping ADM1 parent episodes removed.
duration_episode_rows <- function(d, dur) {
  lvl <- (dur %||% list())$admin_level %||% "ADM1"
  if (!lvl %in% c("ADM1", "ADM2", "ADM1-ADM2")) {
    stop(
      "duration.admin_level must be one of ADM1, ADM2, ADM1-ADM2; got '",
      lvl,
      "'.",
      call. = FALSE
    )
  }
  if (identical(lvl, "ADM1-ADM2")) {
    apply_admin_level_filter(d, "ADM1-ADM2")
  } else {
    d
  }
}

#' load_duration_data
#'
#' Read + filter + scale the duration panel for 04_fit_dur_brms.R. Mirrors
#' read_size_data()'s column-class scaling (siz_data.R:114-125) so a coefficient
#' on `dev_precip` means the same thing in all three tracks.
#'
#' The size track's population guards are deliberately absent: the cessation
#' hazard carries a TIME offset, not a population offset, so a broken WorldPop
#' value cannot distort it.
#'
#' @param opt config list.
#' @param path duration panel RDS; defaults to the configured baseline's file.
load_duration_data <- function(opt, path = NULL) {
  if (is.null(path)) {
    path <- make_duration_dat_file(opt)
  }
  if (!file.exists(path)) {
    stop(
      "load_duration_data: missing input ",
      path,
      " -- run analysis/03_prepare_dur.R first.",
      call. = FALSE
    )
  }

  d <- readRDS(path)
  ep_rec <- assert_episode_key(opt, path, "duration")

  # P3 (geometry-derived admin level); default off -- see load_surv_data().
  if (isTRUE((opt$duration %||% list())$canonical_admin_level)) {
    d <- .apply_p3(d, opt, "duration", tag = "load_duration_data/P3")
  }

  cols <- colnames(d)
  precip_cols <- setdiff(
    grep("precip|cum7", cols, value = TRUE),
    grep("hlq|^dev_sd_", cols, value = TRUE)
  )
  temp_all <- setdiff(
    grep("temp", cols, value = TRUE),
    grep("hlq", cols, value = TRUE)
  )
  temp_cols <- setdiff(temp_all, grep("dev", temp_all, value = TRUE))
  temp_dev_cols <- intersect(temp_all, grep("dev", temp_all, value = TRUE))

  if (length(precip_cols)) {
    d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(precip_cols), ~ . / 100))
  }
  if (length(temp_cols)) {
    d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(temp_cols), ~ . / 10))
  }
  if (length(temp_dev_cols)) {
    d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(temp_dev_cols), ~ . / 1))
  }

  ctry <- opt$country %||% "all"
  if (!identical(ctry, "all")) {
    keep <- trimws(strsplit(ctry, ",")[[1]])
    d <- dplyr::filter(d, country %in% keep)
  }

  # Episode set and ADM2 leaf columns (same order as load_monthly_data(): after
  # P3 and the country filter).
  d <- duration_episode_rows(d, opt$duration)
  d <- .leaf_columns_for(d, opt$duration, ep_rec, tag = "load_duration_data")

  # Rebuild the band from config so a basehaz_duration sensitivity run does not
  # require re-running 03_prepare_dur.R.
  brk <- opt$duration$mso_breaks %||% .DUR_MSO_BREAKS
  d$mso_band <- bin_mso(d$mso, breaks = as.numeric(brk))

  # Calendar-month seasonality (duration.seasonal_term: harmonic2), the same
  # four columns load_monthly_data() builds, on the active month's calendar
  # month. Off by default: nothing is added.
  if (identical(opt$duration[["seasonal_term"]] %||% "none", "harmonic2")) {
    month_num <- lubridate::month(d$month)
    d$sin_month <- sin(2 * pi * month_num / 12)
    d$cos_month <- cos(2 * pi * month_num / 12)
    d$sin2_month <- sin(4 * pi * month_num / 12)
    d$cos2_month <- cos(4 * pi * month_num / 12)
  }

  d$y <- as.integer(d$ended)

  # Empirical baseline-rate anchor (utils/anchor.R), gated on
  # duration.offset_baserate -- OFF by default, so every config that exists
  # today is bit-inert here. Built on the FULL corpus ahead of the country
  # subset above, so a per-group fit resolves against the global fallback tiers.
  if (anchor_enabled(opt, "dur")) {
    d <- attach_base_rate(d, opt, "dur")
  }

  attr(d, "episode_stamp") <- ep_rec
  d
}
