# mon_data.R — fit-time loader of the monthly case-count panel (mon track)
#
# Provides:
#   load_monthly_data()   read + filter the panel for 04_fit_mon_brms(_cb).R
#
# Details:
#   The panel (one row per active month of an episode, response `cases_month`)
#   is built by 03_prepare_mon.R; this module only filters it and adds the
#   column families the formula builder references but nothing upstream
#   builds: the harmonic seasonal columns and cases_prev. No column-class
#   scaling here: the panel inherits already-scaled precip/temp columns from
#   the duration panel it is joined onto (unlike dur_data.R's loader).

#' load_monthly_data
#'
#' Read + filter the monthly-size panel for 04_fit_mon_brms(_cb).R.
#'
#' Order of operations is load-bearing:
#'   1. read the FULL, unfiltered panel;
#'   2. derive `cases_prev` / `log_cases_prev` from THAT full panel (a
#'      calendar-month join within episode, not a row-order lag) -- so the
#'      lookup is immune to any row this function later drops for a country or
#'      admin_level not in the current run;
#'   3. THEN filter by country and admin_level;
#'   4. derive `month_num` and, when `seasonal_term == "harmonic2"`, the four
#'      harmonic columns -- named at config.R's `.seasonal_terms()` but built
#'      nowhere upstream (the size track's non-"none" seasonal terms have the
#'      same gap; this loader is the first place they are actually constructed);
#'   5. `cases_gap` rows are NEVER filtered out -- they are real reporting
#'      gaps within an episode's window (03_prepare_mon.R), already coalesced
#'      to `cases_month = 0`, and dropping them would silently delete the rows
#'      documented there;
#'   5b. P1/P3/P4, the estimand-coherence apparatus ported from read_size_data()
#'      (siz_data.R) and interleaved at the SAME positions, because the ordering
#'      is what makes it correct: `pop` is repointed at the declared denominator
#'      and `admin_level` re-keyed to the canonical GADM level BEFORE
#'      apply_admin_level_filter() (so the row filter, the anchor's pooling key and
#'      the elasticity's centring key are one keying, not three), while the
#'      `log_pop_c` deviation column is built AFTER every row filter (so its
#'      centres are the means of the rows actually fitted, which is what the
#'      simulator re-centres frame units on);
#'   6. round `cases_month` to the nearest integer just before returning --
#'      `03_prepare_mon.R`'s footprint-collapse apportionment (`utils/footprint.R`)
#'      leaves it fractional, and `family = negbinomial()` requires an integer
#'      response (`brms:::data_response.brmsframe()` errors otherwise). Same
#'      rationale as the reference GAM's `cases = round(cases_month)` and
#'      as this same function's own `cases_prev_raw = round(cases_month)` above --
#'      applied last so it never perturbs the cases_prev calendar-lag lookup,
#'      which already rounds independently off the full unfiltered panel.
#'
#' @param opt config list (post make_options_from_config(); `opt$country` is
#'   already the group-resolved subset when grouping.enable is TRUE -- see
#'   resolve_group_countries(), config.R:963).
#' @param path monthly panel RDS; defaults to the configured baseline's file.
#' @return the filtered monthly panel, with month_num / cases_prev / harmonic
#'   columns attached as configured.
load_monthly_data <- function(opt, path = NULL) {
  if (is.null(path)) {
    path <- make_monthly_dat_file(opt)
  }
  if (!file.exists(path)) {
    stop(
      "load_monthly_data: missing input ",
      path,
      " -- run analysis/03_prepare_mon.R first.",
      call. = FALSE
    )
  }

  panel <- readRDS(path)
  ep_rec <- assert_episode_key(opt, path, "monthly")

  # P0: resolve and validate the estimand ONCE, up front, so a mis-declared
  # config fails before any data work rather than producing a plausible-looking
  # fit on mixed units (read_size_data()'s rule, same reason).
  estimand <- resolve_monthly_estimand(opt)

  # --- cases_prev, from the FULL unfiltered panel (the reference GAM's rule) --
  mon <- opt$monthly %||% list()
  if (isTRUE(mon$prev_cases)) {
    prev_lookup <- panel %>%
      dplyr::transmute(episode_id, month, cases_prev_raw = round(cases_month))

    panel <- panel %>%
      dplyr::mutate(month_prev = month %m-% months(1)) %>%
      dplyr::left_join(
        prev_lookup,
        by = c("episode_id", "month_prev" = "month")
      ) %>%
      dplyr::mutate(
        cases_prev = tidyr::replace_na(cases_prev_raw, 0),
        log_cases_prev = log1p(cases_prev)
      ) %>%
      dplyr::select(-month_prev, -cases_prev_raw)
  }

  d <- panel

  # --- country filter (mirrors load_duration_data(); opt$country is already
  # group-resolved by the time this runs) --------------------------------
  ctry <- opt$country %||% "all"
  if (!identical(ctry, "all")) {
    keep <- trimws(strsplit(ctry, ",")[[1]])
    d <- dplyr::filter(d, country %in% keep)
  }

  # --- P1: point `pop` at the declared denominator ------------------------
  # The monthly mean is mu = pop x exp(eta) with `offset(log(pop))` FIXED, so
  # `pop`'s meaning is the whole absolute scale and nothing in the likelihood
  # can correct it. The panel `pop` is the LP reporting FOOTPRINT H, while
  # the simulator replays at the whole GADM unit Q -- every replayed mu inflated
  # by 1/coverage. `footprint` keeps exactly that behaviour; `gid` is a hard
  # error until 03_prepare_dur.R/03_prepare_mon.R have been re-run to carry the
  # P1 triple through (dur_data.R meta_cols).
  d <- resolve_size_pop_column(d, estimand, tag = "load_monthly_data")

  # --- P3: canonical admin level, BEFORE the admin_level row filter --------
  # Load-bearing position, exactly as in read_size_data(): `admin_level` is the
  # row filter's partition key, the mon anchor's pooling key AND the elasticity
  # centring key. Re-keying after the filter would leave the filter working on
  # corpus labels while the other two pooled on canonical ones -- three keyings
  # of one dataset, which is the defect. Gated on the estimand being opted in,
  # so pre-P0 configs are untouched.
  d <- .apply_p3(d, opt, "monthly", estimand, tag = "load_monthly_data/P3")

  # --- admin_level filter (shared helper; "ADM1-ADM2" keeps both, dropping
  # only the coarse ADM1 units that have an ADM2 sibling) -----------------
  d <- apply_admin_level_filter(
    d,
    admin_level_opt = mon$admin_level %||% "ADM1-ADM2"
  )

  # re_adm2_leaf_scope: adm2_only -> is_adm2 / gid_leaf (utils/paths.R).
  d <- .leaf_columns_for(d, mon, ep_rec, tag = "load_monthly_data")

  # --- calendar-month seasonality ------------------------------------------
  d <- d %>% dplyr::mutate(month_num = lubridate::month(month))

  if (identical(mon$seasonal_term %||% "none", "harmonic2")) {
    d <- d %>%
      dplyr::mutate(
        sin_month = sin(2 * pi * month_num / 12),
        cos_month = cos(2 * pi * month_num / 12),
        sin2_month = sin(4 * pi * month_num / 12),
        cos2_month = cos(4 * pi * month_num / 12)
      )
  }

  # --- integer response for family = negbinomial() (point 6 above) ---------
  d <- d %>% dplyr::mutate(cases_month = round(cases_month))

  # Empirical baseline-rate anchor (utils/anchor.R), gated on
  # monthly.offset_baserate -- OFF by default, so every config that exists today
  # is bit-inert here. Built on the FULL corpus ahead of the country and
  # admin-level subsets above, so a per-group fit resolves against the global
  # fallback tiers. `.mon_anchor_view()` rounds `cases_month` the same way the
  # line above does, so the anchor is a rate for the response actually fitted.
  if (anchor_enabled(opt, "mon")) {
    d <- attach_base_rate(d, opt, "mon")
  }

  # P4: the elasticity deviation column. Built AFTER all row filters so the
  # centres are the means of the rows actually fitted -- centring on the
  # unfiltered panel would put the pivot outside the training set, and the
  # simulator (which reuses these centres via the stamp below) would extrapolate
  # around a population no row ever had.
  if (!identical(estimand$pop_elasticity, "fixed")) {
    d <- .add_log_pop_centred(d, tag = "load_monthly_data/P4")
  }

  attr(d, "mon_estimand") <- monthly_estimand_stamp(
    opt,
    extra = list(
      n_rows = nrow(d),
      cases_total = sum(d$cases_month, na.rm = TRUE),
      admin_level = mon$admin_level %||% NA_character_,
      log_pop_centres = attr(d, "log_pop_centres")
    )
  )
  attr(d, "episode_stamp") <- ep_rec
  d
}
