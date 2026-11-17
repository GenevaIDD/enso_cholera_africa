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
#' Read + filter the monthly-size panel for 04_fit_mon_brms(_cb).R. The order
#' of operations is load-bearing:
#'   1. read the full, unfiltered panel;
#'   2. derive `cases_prev` / `log_cases_prev` from that full panel (a
#'      calendar-month join within episode, immune to later row drops);
#'   3. filter by country and admin_level;
#'   4. derive `month_num` and, for `seasonal_term == "harmonic2"`, the four
#'      harmonic columns (built here, nowhere upstream);
#'   5. never drop `cases_gap` rows: they are real reporting gaps, already 0;
#'   5b. the estimand-coherence steps from read_size_data() at the same
#'      positions: repoint `pop` and re-key `admin_level` before
#'      apply_admin_level_filter(); build `log_pop_c` after every row filter;
#'   6. round `cases_month` last (the footprint apportionment leaves it
#'      fractional; negbinomial needs integers).
#'
#' @param opt config list (post make_options_from_config(); `opt$country` is
#'   already the group-resolved subset when grouping.enable is TRUE -- see
#'   config.R::resolve_group_countries()).
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

  # P0: resolve and validate the estimand once, up front, so a mis-declared
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
