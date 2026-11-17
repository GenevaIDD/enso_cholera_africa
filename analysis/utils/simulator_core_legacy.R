# simulator_core_legacy.R — frozen reference bodies of the forward-kernel
#   functions, before the bit-identical speed-up of simulator_core.R
#
# Provides:
#   simulate_one_tu_legacy(), build_warmup_state_legacy(),
#   evaluate_lp_legacy(), update_lags_legacy()   and the other `_legacy`
#                                  copies they call
#
# Details:
#   Never modified. Each body is renamed `<name>_legacy` and calls the other
#   `_legacy` copies; everything else (re_level_key(), compute_phi_vec(),
#   build_cb_mats inputs, ...) is the unchanged simulator_core.R function.
#   The speed-up changed: evaluate_lp() lookup plan and per-draw coefficient
#   cache, integer-coded update_lags(), index-based warm-up rows, vector event
#   bookkeeping, .bincode() binning, no ifelse() in the step. Only
#   tests/test_simulator_kernel_golden.R calls it, to compare the two kernels
#   live. Callers source simulator_core.R and its closure.

build_warmup_state_legacy <- function(
  panel_tu,
  frame_tu,
  horizon,
  has_dur = FALSE,
  no_prior_onset = NULL
) {
  start <- as.Date(horizon$start)

  pre <- panel_tu %>% dplyr::filter(t_start < start)
  warm <- if (nrow(pre) > 0) {
    pre %>%
      dplyr::group_by(gid) %>%
      dplyr::filter(t_start == max(t_start)) %>%
      dplyr::ungroup()
  } else {
    pre[0, , drop = FALSE]
  }

  # If we have no pre-horizon rows for some location, fall back to the first
  # post-horizon row (so the kernel can still launch in tests where no warm-up
  # tail was provided).
  if (nrow(warm) < nrow(frame_tu)) {
    missing_gids <- setdiff(frame_tu$gid, warm$gid)
    cand <- panel_tu %>% dplyr::filter(gid %in% missing_gids)
    if (nrow(cand) > 0) {
      first_post <- cand %>%
        dplyr::group_by(gid) %>%
        dplyr::filter(t_start == min(t_start)) %>%
        dplyr::ungroup()
      warm <- dplyr::bind_rows(warm, first_post)
    }
  }

  # Align rows with frame_tu order
  warm <- warm[match(frame_tu$gid, warm$gid), , drop = FALSE]

  pick <- function(col, default = 0) {
    if (col %in% names(warm)) {
      v <- warm[[col]]
      v[is.na(v)] <- default
      v
    } else {
      rep(default, nrow(warm))
    }
  }

  # u: step 4 — initialise as days from the panel's factor-level TTLO column.
  ttlo_cat <- if ("time_to_last_outbreak_cat" %in% names(warm)) {
    warm[["time_to_last_outbreak_cat"]]
  } else {
    rep(NA_character_, nrow(warm))
  }
  u_days <- ttlo_init_from_panel(ttlo_cat)
  # "none" level: a unit with no recorded onset before the
  # launch (no_prior_onset, frame_tu order; ttlo_no_prior_onset_at()) whose panel
  # category is NA or [6,Inf] starts at the sentinel. NULL (default, and every
  # fit without the level) leaves u_days as above.
  if (!is.null(no_prior_onset)) {
    if (length(no_prior_onset) != nrow(frame_tu)) {
      stop(
        "build_warmup_state: no_prior_onset must have one value per frame_tu row."
      )
    }
    cat_chr <- as.character(ttlo_cat)
    seed <- no_prior_onset %in%
      TRUE &
      (is.na(cat_chr) | cat_chr == .TTLO_LEVELS[length(.TTLO_LEVELS)])
    u_days[seed] <- .TTLO_NONE_DAYS
  }

  # o_tu: locations with NA TU membership have ongoing_outbreak_tu_lag1 = NA
  # at fit time — preserve that NA so the spillover term is absent from the LP.
  o_tu_val <- as.integer(pick("ongoing_outbreak_tu_lag1", default = 0L))
  o_tu_val[is.na(frame_tu$tu)] <- NA_integer_

  # active/mso: duration-aware state. No historical source populates
  # these into panel_tu (the continuous climate panel) -- that handoff is
  # the smoothing cache's job. Without it,
  # every location warm-starts as inactive; flag it loudly (not silently) so
  # a duration-aware run without a real warm-start is a visible choice, not
  # an accident.
  active0 <- as.integer(pick("active", default = 0L))
  mso0 <- as.integer(pick("mso", default = 0L))
  if (
    isTRUE(has_dur) &&
      !("active" %in% names(warm)) &&
      !("mso" %in% names(warm))
  ) {
    warning(
      "build_warmup_state: duration draws supplied (has_dur = TRUE) but ",
      "panel_tu carries no 'active'/'mso' warm-up columns -- every location ",
      "starts the duration-aware simulation as inactive (active = 0, ",
      "mso = 0). Expected for v1 runs with no smoothing-cache warm start ",
      "(Phase 9); if a warm start was expected, check that panel_tu / ",
      "warmup_panels_list actually carries 'active'/'mso'.",
      call. = FALSE
    )
  }

  list(
    u = u_days,
    o = as.integer(pick("ongoing_outbreak_lag1", 0)),
    o_ct = as.integer(pick("ongoing_outbreak_country_lag1", 0)),
    o_tu = o_tu_val,
    active = active0,
    mso = mso0,
    # Detection flag (covariate_feed: detected) -- only when the warm start
    # carries it (a detected-feed smoothing handoff). Consumers default it to
    # `active`: an observed-panel warm start is detected by construction.
    det = if ("det" %in% names(warm)) as.integer(pick("det", 0L)) else NULL
  )
}

.onset_log_offset_adj_legacy <- function(mso, use_logoffset = TRUE) {
  if (!isTRUE(use_logoffset)) {
    return(0)
  }
  ifelse(as.integer(mso) == 1L, log(.ONSET_EXPOSURE_FRAC), 0)
}

evaluate_lp_legacy <- function(
  draws,
  s,
  xrow,
  o,
  o_ct,
  o_tu,
  u_days,
  frame_tu,
  basehaz_val = 0,
  drop_region = FALSE,
  drop_ct_tu = FALSE,
  nbr = NULL
) {
  N <- nrow(frame_tu)
  fe <- draws$fe[s, , drop = TRUE]

  lp <- rep(0, N)

  # Intercept (population-level)
  if ("(Intercept)" %in% names(fe)) {
    lp <- lp + fe[["(Intercept)"]]
  }
  if ("Intercept" %in% names(fe)) {
    lp <- lp + fe[["Intercept"]]
  }

  # Dual-key lookup helpers. Grouped brms_cb fits have mixed nesting: intercept
  # and spillover REs use bare country keys ("KEN") while ENSO RE slopes use
  # AFRO-prefixed keys ("Central.Africa_KEN"). Neither single drop_region value
  # matches both, so we try the primary format first and fall back to the
  # alternate format for any location that doesn't match.
  re_keys <- function(bucket) {
    list(
      primary = re_level_key(frame_tu, bucket, drop_region),
      alt = re_level_key(frame_tu, bucket, !drop_region)
    )
  }
  match_keys <- function(keys, cn) {
    idx <- match(keys$primary, cn)
    miss <- is.na(idx)
    if (any(miss)) {
      idx[miss] <- match(keys$alt[miss], cn)
    }
    idx
  }

  # Nested intercept REs (region -> country -> adm1 -> adm2)
  add_intercept_re <- function(lp, mat, bucket) {
    if (is.null(mat) || ncol(mat) == 0) {
      return(lp)
    }
    idx <- match_keys(re_keys(bucket), colnames(mat))
    have <- !is.na(idx)
    if (any(have)) {
      lp[have] <- lp[have] + mat[s, idx[have]]
    }
    lp
  }
  lp <- add_intercept_re(lp, draws$re_intercept$region, "region")
  lp <- add_intercept_re(lp, draws$re_intercept$country, "country")
  lp <- add_intercept_re(lp, draws$re_intercept$adm1, "adm1")
  lp <- add_intercept_re(lp, draws$re_intercept$adm2, "adm2")

  # Sum the per-location RE slope contribution for one term across every nested
  # grouping level present in the bundle. A `/`-nested brms slope such as
  # `(enso_lag5 - 1 | AFRO_region/country/ADM1)` expands to region,
  # region:country and region:country:ADM1 group factors — the kernel must add
  # all three (the previous code only read the country, or only the adm1,
  # bucket, silently dropping the others).
  nested_slope <- function(term) {
    rst <- draws$re_slope[[term]]
    out <- rep(0, N)
    if (is.null(rst)) {
      return(out)
    }
    lk <- function(mat, bucket) {
      if (is.null(mat) || ncol(mat) == 0) {
        return(rep(0, N))
      }
      idx <- match_keys(re_keys(bucket), colnames(mat))
      ifelse(is.na(idx), 0, mat[s, idx])
    }
    out +
      lk(rst$region, "region") +
      lk(rst$country, "country") +
      lk(rst$adm1, "adm1") +
      lk(rst$adm2, "adm2")
  }

  # Linear covariate FE + nested RE slopes. ENSO (`enso_lag<n>`) is detected
  # dynamically from the FE names so any configured lag is consumed — this is
  # the headline covariate of the elag5 model.
  # DLNM cross-basis columns (`<var>_cb<k>`, e.g. enso_cb1, mean_precip_cb3) are
  # consumed the same way: each is a flat FE coefficient × its expanded basis
  # value, plus any nested RE slope (re_adm1_slope). The basis values are built
  # by build_cb_design() upstream and threaded in via xrow; here they are just
  # additional linear terms. Empty when the fit carries no cross-basis.
  cb_terms <- grep("_cb[0-9]+$", names(fe), value = TRUE)
  linear_terms <- c(
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    grep("^enso_(lag|lead)[0-9]+$", names(fe), value = TRUE),
    # Distance-weighted neighbour contagion (nbr_lag1_r<rho>[_rn]),
    # detected generically like ENSO so whatever rho the fit carries
    # is scored. The per-step value is reconstructed dynamically by
    # simulate_one_tu_legacy() / the spatPomp R rprocess and threaded in via
    # the `nbr` arg (written onto xrow just below); it is a PLAIN
    # linear term — (fe + nested_slope) * nbr — NOT a FALSE-factor
    # spillover. A no-op when `nbr` is NULL (xrow lacks the column).
    grep(nbr_covariate_pattern(), names(fe), value = TRUE),
    cb_terms,
    # (`logsize_cases`, the duration track's E1 size covariate,
    # is no longer scored: the size track is retired from the
    # simulator and .extract_dur_draws_brms() refuses such fits.)
    # Monthly case-count (`mon`) track terms. Absent from every
    # occurrence/duration `fe`, so adding them here is a
    # no-op for those tracks (the `%in% names(fe)` and
    # `is.null(xv)` guards below make that unconditional).
    #   * the four harmonics of monthly.seasonal_term:
    #     "harmonic2" -- the shipped default
    #     (config_defaults.yml:485), built by load_monthly_data()
    #     (utils/mon_data.R:92-99) at fit time and by the kernel
    #     from the calendar month at simulation time;
    #   * the ADM2 dummy of monthly.admin_level: "ADM1-ADM2",
    #     which is what unlocks the gid-leaf intercept RE
    #     (.intercept_re_term(), utils/config.R:1677).
    # mon's `mso_cb*` and `enso_cb*` need no entry: `cb_terms`
    # above already matches them generically.
    "sin_month",
    "cos_month",
    "sin2_month",
    "cos2_month",
    "admin_levelADM2",
    # P4 population-elasticity deviation term
    # (`size.pop_elasticity: free` / `free_by_level`,
    # utils/config.R; built at fit time by
    # `.add_log_pop_centred()`, utils/siz_data.R:75). Unlike every
    # other entry above, no caller puts this on `xrow` -- it is
    # a pure function of the frame's own population, so it is
    # derived from `frame_tu` just below. It was missing from
    # this whitelist entirely, which meant `b_log_pop_c` was
    # silently dropped at replay: the kernel applied the
    # hard-coded elasticity 1 (via offset(log(pop))) while the
    # fit had calibrated eta against an estimated beta.
    "log_pop_c"
  )

  # P4: rebuild `log_pop_c` from the frame, around the TRAINING centres carried
  # on the draws bundle (`.recover_log_pop_centres()`). Centring on the frame's
  # own mean would move the pivot and rescale every prediction by
  # exp((beta - 1) * shift), so an absent centre set is a hard error, not a 0.
  if ("log_pop_c" %in% names(fe) && is.null(xrow[["log_pop_c"]])) {
    ctr <- draws$log_pop_centres
    if (is.null(ctr) || length(ctr) == 0L) {
      stop(
        "evaluate_lp: the fit has a `log_pop_c` fixed effect but the draws ",
        "bundle carries no `log_pop_centres`. See ",
        ".recover_log_pop_centres() -- there is no safe default (a 0 pivot ",
        "would rescale every prediction by exp((beta - 1) * mean(log pop))).",
        call. = FALSE
      )
    }
    if (!("pop" %in% names(frame_tu))) {
      stop(
        "evaluate_lp: `log_pop_c` needs `pop` on the location frame.",
        call. = FALSE
      )
    }
    # The centres are named by the TRAINING `admin_level`, which
    # `apply_canonical_admin_level()` (P3) has already re-keyed to GADM depth by the
    # time `.add_log_pop_centred()` runs -- so the replay side must key on the
    # frame's canonical column, not the corpus label (see
    # `.frame_admin_level_canonical()`). `admin_level` is only guaranteed on the
    # frame when build_location_frame() had to resolve a mean_rate/anchor cell; a
    # frame that arrived with its rate column already attached can lack it
    # entirely. Fall back to the SAME gid depth rule
    # .recover_log_pop_centres() uses on the training side, so the two sides key
    # the centre lookup identically instead of hard-failing.
    lvl <- .frame_admin_level_canonical(frame_tu)
    if (length(lvl) != N) {
      lvl <- rep(lvl, length.out = N)
    }
    if (anyNA(lvl) && "gid" %in% names(frame_tu)) {
      na_lvl <- is.na(lvl)
      lvl[na_lvl] <- .admin_level_from_gid(frame_tu$gid[na_lvl])
    }
    lvl[is.na(lvl)] <- "ALL"
    # Same "resolve at the coarsest available match" rule the anchor lookup
    # uses (.resolve_frame_anchor()): a frame unit whose admin_level was never in
    # the training corpus borrows its coarsest available parent's centre rather
    # than dropping to NA.
    lvl <- .coarsen_admin_level_to_available(lvl, names(ctr))
    c_i <- unname(ctr[lvl])
    if (length(ctr) == 1L) {
      c_i[is.na(c_i)] <- unname(ctr[[1L]])
    }
    if (anyNA(c_i)) {
      stop(
        "evaluate_lp: no training `log_pop_c` centre for admin_level(s) ",
        paste(unique(lvl[is.na(c_i)]), collapse = ", "),
        " (fit centres: ",
        paste(names(ctr), collapse = ", "),
        ").",
        call. = FALSE
      )
    }
    xrow[["log_pop_c"]] <- log(as.numeric(frame_tu$pop)) - c_i
  }

  # Neighbour-contagion value is supplied dynamically per step (NOT folded into
  # an exogenous xrow upstream): write the recomputed pressure onto the detected
  # nbr column so the generic linear-term loop below scores it exactly like any
  # other continuous covariate. NULL `nbr` ⇒ no write ⇒ term skipped (is.null xv).
  if (!is.null(nbr)) {
    for (ncol_ in grep(nbr_covariate_pattern(), names(fe), value = TRUE)) {
      xrow[[ncol_]] <- nbr
    }
  }

  for (term in linear_terms) {
    if (!(term %in% names(fe))) {
      next
    }
    xv <- xrow[[term]]
    if (is.null(xv)) {
      next
    }
    xv0 <- ifelse(is.na(xv), 0, xv)
    beta_s <- fe[[term]]
    lp <- lp + (beta_s + nested_slope(term)) * xv0
  }

  # Spillover lag FE + nested RE slopes (region + country at fit time).
  # NA lag_vec (NA-TU locations in o_tu) → 0 contribution, matching fit-time
  # behaviour where rows with NA ongoing_outbreak_tu_lag1 are excluded.
  add_spillover <- function(lp, term, lag_vec) {
    if (!(term %in% names(fe))) {
      return(lp)
    }
    alpha_s <- fe[[term]]
    contrib_true <- (alpha_s + nested_slope(term)) * lag_vec
    contrib_true[is.na(lag_vec)] <- 0 # NA lag → term absent at fit time → 0
    # brms stores RE slopes for ALL factor levels, including FALSE.
    # Apply the FALSE-level slope to rows where lag_vec == 0 (NA treated as 0).
    false_key <- paste0(term, "FALSE")
    contrib_false <- if (!is.null(draws$re_slope[[false_key]])) {
      ns_false <- nested_slope(false_key)
      ns_false * ifelse(is.na(lag_vec), 0L, 1L - as.integer(lag_vec))
    } else {
      rep(0, N)
    }
    lp + contrib_true + contrib_false
  }
  lp <- add_spillover(lp, "ongoing_outbreak_lag1", o)
  if (!isTRUE(drop_ct_tu)) {
    lp <- add_spillover(lp, "ongoing_outbreak_country_lag1", o_ct)
    lp <- add_spillover(lp, "ongoing_outbreak_tu_lag1", o_tu)
  }

  # TTLO factor-level FE + nested RE slopes (step 4).
  # Each non-reference bin [1,2), …, [6,Inf] maps to a brms FE column like
  # `time_to_last_outbreak12`. The reference bin [0,1) has no FE contribution
  # (absorbed into the intercept) but brms still stores RE slopes for it under
  # the key `time_to_last_outbreak01` — these must be applied.
  # NA u_days (warm-up not available) → 0 contribution.
  ttlo_cols <- grep("^time_to_last_outbreak", names(fe), value = TRUE)
  if (length(ttlo_cols) > 0) {
    # "none" only when this bundle's fit carries the level;
    # otherwise the sentinel bins to [6,Inf] as the fit coded it.
    bins <- bin_ttlo_legacy(u_days, has_none = .TTLO_NONE_COL %in% ttlo_cols) # NA for NA input
    for (col_name in ttlo_cols) {
      level <- ttlo_col_to_level(col_name)
      in_bin <- !is.na(bins) & (bins == level)
      if (!any(in_bin)) {
        next
      }
      fe_val <- fe[[col_name]]
      re_val <- nested_slope(col_name)
      lp[in_bin] <- lp[in_bin] + fe_val + re_val[in_bin]
    }
    # Reference bin [0,1): no FE term but RE slope is stored and must be applied.
    ref_key <- "time_to_last_outbreak01"
    if (!is.null(draws$re_slope[[ref_key]])) {
      in_ref <- !is.na(bins) & (bins == .TTLO_LEVELS[1L]) # "[0,1)"
      if (any(in_ref)) {
        re_ref <- nested_slope(ref_key)
        lp[in_ref] <- lp[in_ref] + re_ref[in_ref]
      }
    }
  }

  # Seasonal spline (e.g. s(t_start_seas, bs="cc")) — coefficients stored in
  # draws$smooth_coef / draws$smooth_pred_fn, NOT in b_* FE columns.
  if (length(draws$smooth_pred_fn) > 0L && length(draws$smooth_coef) > 0L) {
    t_ss <- xrow[["t_start_seas"]]
    if (!is.null(t_ss)) {
      for (sm_nm in names(draws$smooth_pred_fn)) {
        pred_fn <- draws$smooth_pred_fn[[sm_nm]]
        coef_mat <- draws$smooth_coef[[sm_nm]]
        if (is.null(pred_fn) || is.null(coef_mat)) {
          next
        }
        X_s <- pred_fn(t_ss) # N × K
        lp <- lp + as.numeric(X_s %*% coef_mat[s, ])
      }
    }
  }

  # Add basehaz (per-location constant for this step; see basehaz_eval below)
  lp + basehaz_val
}

.dur_bin_mso_legacy <- function(mso, breaks = .SIM_DUR_MSO_BREAKS) {
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
  as.character(cut(
    as.numeric(mso),
    breaks = breaks,
    labels = labs,
    right = FALSE,
    include.lowest = TRUE
  ))
}

basehaz_dur_eval_legacy <- function(draws_dur, s, mso, N) {
  bh <- draws_dur$basehaz
  switch(
    bh$kind %||% "none",
    factor = {
      out <- rep(0, N)
      lvl <- .dur_bin_mso_legacy(mso, bh$breaks %||% .SIM_DUR_MSO_BREAKS)
      ok <- !is.na(lvl)
      if (any(ok)) {
        cols <- .dur_mso_level_to_col(lvl[ok])
        idx <- match(cols, colnames(draws_dur$fe))
        got <- !is.na(idx)
        if (any(got)) out[which(ok)[got]] <- draws_dur$fe[s, idx[got]]
      }
      out
    },
    none = rep(0, N),
    rep(0, N)
  )
}

bin_ttlo_legacy <- function(days, has_none = FALSE) {
  years <- (pmax(0, days) + 1) / 365
  out <- as.character(
    cut(
      years,
      .TTLO_BREAKS,
      labels = .TTLO_LEVELS,
      include.lowest = TRUE,
      right = FALSE
    )
  )
  if (isTRUE(has_none)) {
    out[!is.na(days) & days >= .TTLO_NONE_DAYS] <- .TTLO_NONE_LEVEL
  }
  out
}

update_lags_legacy <- function(Y, frame_tu) {
  N <- length(Y)
  # ADM1 active state: any location in ADM1 j active at m
  A_adm1 <- tapply(Y, frame_tu$adm1_unit, max, default = 0L)
  # Country active state: any ADM1 of country c active
  adm1_to_country <- frame_tu %>%
    dplyr::distinct(adm1_unit, country) %>%
    tibble::deframe()
  # Order A_adm1 vector by names; build per-country sum + max
  adm1_names <- names(A_adm1)
  ctry_of_adm1 <- adm1_to_country[adm1_names]

  # ADM1-by-country aggregates
  A_ct_max <- tapply(as.integer(A_adm1), ctry_of_adm1, max, default = 0L)
  A_ct_sum <- tapply(as.integer(A_adm1), ctry_of_adm1, sum, default = 0L)

  # TU active state per country: any other country of TU active
  ctry_to_tu <- frame_tu %>%
    dplyr::distinct(country, tu) %>%
    tibble::deframe()
  ctry_names <- names(A_ct_max)
  tu_of_ctry <- ctry_to_tu[ctry_names]

  # tapply silently drops NA keys → NA-TU countries excluded from TU aggregates
  A_tu_sum <- tapply(as.integer(A_ct_max), tu_of_ctry, sum, default = 0L)

  # Broadcast back to per-location
  o_next <- as.integer(A_adm1[frame_tu$adm1_unit])
  # Country leave-one-ADM1-out: (sum over ADM1 of country c — self) >= 1
  ct_for_loc <- frame_tu$country
  ct_sum_loc <- as.integer(A_ct_sum[ct_for_loc])
  self_adm1 <- as.integer(A_adm1[frame_tu$adm1_unit])
  o_ct_next <- as.integer((ct_sum_loc - self_adm1) >= 1L)

  # TU leave-one-country-out: (sum over countries of TU — self_country) >= 1
  # Locations with NA TU → NA (never coerced to 0 — they match fit-time NA rows).
  tu_for_loc <- frame_tu$tu
  tu_sum_loc <- as.integer(ifelse(
    is.na(tu_for_loc),
    NA_integer_,
    A_tu_sum[tu_for_loc]
  ))
  self_country <- as.integer(A_ct_max[ct_for_loc])
  o_tu_next <- as.integer((tu_sum_loc - self_country) >= 1L)
  # NA preserved: do NOT zero-fill NA-TU rows.

  # A_adm1 (named by ADM1 unit, the lag-1 ADM1-active 0/1 grid) is returned so
  # the neighbour-contagion caller can form `K %*% Aa1` without recomputing the
  # same tapply. Existing consumers ignore the extra field.
  list(
    o = o_next,
    o_ct = o_ct_next,
    o_tu = o_tu_next,
    A_adm1 = as.integer(A_adm1) |> stats::setNames(names(A_adm1))
  )
}

advance_ttlo_legacy <- function(u_days, Y, step_days = 365.25 / 12) {
  ifelse(as.integer(Y) == 1L, 0.0, u_days + step_days)
}

monthlyLogMu_legacy <- function(
  draws,
  s,
  xm,
  o,
  o_ct,
  o_tu,
  u_days,
  frame_g,
  log_offset,
  mso,
  mso_tab = NULL,
  what = "monthly kernel",
  drop_region = FALSE,
  drop_ct_tu = FALSE,
  pop_floor = NULL,
  pop = NULL
) {
  if (!is.null(pop)) {
    frame_g$pop <- pop
  }
  xm <- .inject_mso(xm, mso_tab, mso)
  lm <- evaluate_lp_legacy(
    draws,
    s,
    xm,
    o,
    o_ct,
    o_tu,
    u_days,
    frame_g,
    basehaz_val = 0,
    drop_region = drop_region,
    drop_ct_tu = drop_ct_tu
  )
  if (isTRUE(draws$offset_pop %||% TRUE)) {
    lm <- lm +
      log(if (is.null(pop_floor)) frame_g$pop else pmax(frame_g$pop, pop_floor))
  }
  # Onset-month exposure correction: `log_offset` is log(days_in_month(t)),
  # exact for interior months but ~2x too long for the onset month the mon
  # panel carries. See .ONSET_EXPOSURE_FRAC.
  if (isTRUE(draws$offset_logoffset %||% FALSE)) {
    lm <- lm + log_offset + .onset_log_offset_adj_legacy(mso)
  }
  # offset(log(base_rate_mon)): its denominator is pop x exposure days, the
  # two offsets above. 0 for an un-anchored fit.
  lm + .anchor_log_offset(draws, frame_g, NULL, what)
}

simulate_one_tu_legacy <- function(
  tu,
  frame_tu,
  panel_tu,
  draws_occ,
  replicate_ids,
  posterior_idx_occ = NULL,
  seed,
  output_part,
  horizon,
  dt_months = 1L,
  out_fmt = "qs",
  dry_run = FALSE,
  enso_scenarios = NULL,
  warmup_panels_list = NULL,
  cb_meta_occ = NULL,
  grouping = NULL,
  force_o_tu = FALSE,
  draws_dur = NULL,
  cb_meta_dur = NULL,
  posterior_idx_dur = NULL,
  draws_mon = NULL,
  cb_meta_mon = NULL,
  posterior_idx_mon = NULL,
  drop_country_tu_spillover = FALSE,
  use_nbr = FALSE,
  nbr_rho_km = 500L,
  nbr_gadm_path = "data/geodata/gadm_valid.gpkg",
  nbr_cache_dir = "generated_data/",
  nbr_spec = NULL,
  covariate_feed = "latent",
  delta_fn = 0.10,
  delta_fp = 0.05,
  clamp_enso_to_knots = FALSE,
  clamp_enso_margin = 0,
  incidence_cap = NULL
) {
  stopifnot(nrow(frame_tu) > 0, !is.null(draws_occ))
  R <- length(replicate_ids)
  N <- nrow(frame_tu)
  # The duration-aware path is opt-in via a non-NULL draws_dur bundle.
  # See the @param doc above for the hard-branch rationale.
  has_dur <- !is.null(draws_dur)
  # Detected-activity covariate feed (update_detection()). Needs the duration
  # track: only `active` defines an episode whose detection persists.
  if (!covariate_feed %in% c("latent", "detected")) {
    stop("simulate_one_tu: covariate_feed must be 'latent' or 'detected'.")
  }
  use_det <- identical(covariate_feed, "detected")
  if (use_det && !has_dur) {
    stop("simulate_one_tu: covariate_feed = 'detected' requires draws_dur.")
  }
  # Monthly case-count (`mon`) path: opt-in via a non-NULL draws_mon bundle.
  # The sole case estimand (the size track is retired).
  has_mon <- !is.null(draws_mon)
  if (has_mon && !has_dur) {
    stop(
      "simulate_one_tu: the monthly case-count track requires duration draws. ",
      "`mso` and `active` are latent states produced by the duration model, ",
      "with no defensible fallback -- set `simulator.duration_fit` (or ",
      "`simulator.grouping.duration_enable`) alongside `simulator.monthly_fit`."
    )
  }

  # Country-grouping routing. When `grouping` is supplied, each
  # country in this TU is fit by its group's own occ (+ dur/mon) model (region
  # dropped from the RE nesting). The LP evaluation below is partitioned by
  # group: every group's rows are scored with that group's draws bundle,
  # cross-basis metadata and region-dropped RE key space, then scattered back
  # into the length-N vectors. `grouping` carries `draws_occ_by_group`,
  # `cb_meta_occ_by_group` and, per enabled track, `draws_{dur,mon}_by_group` /
  # `cb_meta_{dur,mon}_by_group` (named by group); the per-row
  # group label comes from `frame_tu$group`. Lag bookkeeping (update_lags) stays
  # TU-level so cross-group spillover within a TU is preserved. When NULL the
  # single-bundle path below is byte-identical to the pre-grouping kernel.
  is_grouped <- !is.null(grouping)
  if (is_grouped) {
    if (!("group" %in% names(frame_tu))) {
      stop(
        "simulate_one_tu: grouped routing requires a `group` column on frame_tu."
      )
    }
    grp_vec <- as.character(frame_tu$group)
    ctx_groups <- unique(grp_vec)
    if (any(is.na(grp_vec))) {
      stop(
        "simulate_one_tu: grouped routing found NA group labels in frame_tu."
      )
    }
  } else {
    grp_vec <- rep("__ALL__", N) # non-empty sentinel (valid list name)
    ctx_groups <- "__ALL__"
  }

  # Raw (pre-clamp) per-replicate draw index. Under grouping each group bundle
  # has its own S, so the clamp is applied per group inside the loop; the raw
  # index keeps the same replicate→draw pairing across all groups.
  raw_idx_occ <- if (is.null(posterior_idx_occ)) {
    replicate_ids
  } else {
    as.integer(posterior_idx_occ)
  }
  if (is.null(posterior_idx_occ)) {
    posterior_idx_occ <- ((replicate_ids - 1L) %% draws_occ$S) + 1L
  }
  # Clamp supplied indices to the available number of draws — the caller may
  # not know how many draws each backend exposes.
  posterior_idx_occ <- ((as.integer(posterior_idx_occ) - 1L) %% draws_occ$S) +
    1L

  # Duration draw index — same paired-strategy clamp, only meaningful when
  # has_dur. raw_idx_dur/posterior_idx_dur stay NULL on the legacy path so
  # nothing downstream can accidentally reference them.
  raw_idx_dur <- NULL
  if (has_dur) {
    raw_idx_dur <- if (is.null(posterior_idx_dur)) {
      replicate_ids
    } else {
      as.integer(posterior_idx_dur)
    }
    if (is.null(posterior_idx_dur)) {
      posterior_idx_dur <- ((replicate_ids - 1L) %% draws_dur$S) + 1L
    }
    posterior_idx_dur <- ((as.integer(posterior_idx_dur) - 1L) %% draws_dur$S) +
      1L
  }

  # Monthly draw index — same paired-strategy clamp as occ/dur.
  raw_idx_mon <- NULL
  if (has_mon) {
    raw_idx_mon <- if (is.null(posterior_idx_mon)) {
      replicate_ids
    } else {
      as.integer(posterior_idx_mon)
    }
    if (is.null(posterior_idx_mon)) {
      posterior_idx_mon <- ((replicate_ids - 1L) %% draws_mon$S) + 1L
    }
    posterior_idx_mon <- ((as.integer(posterior_idx_mon) - 1L) %% draws_mon$S) +
      1L
  }

  # Months sequence (one row per gid per month). The panel may be irregular;
  # we drive the simulation by the unique sorted t_start values per gid.
  # The warm-up tail (t_start < horizon start) is read by build_warmup_state
  # to seed the AR state but must NOT be simulated/emitted: restrict the
  # simulated grid to the horizon proper so events fall only inside [start, end].
  months <- sort(unique(panel_tu$t_start))
  months <- months[months >= as.Date(horizon$start)]
  M <- length(months)
  if (M == 0) {
    warning(
      sprintf(
        "[simulate_one_tu] TU '%s': no panel months at or after horizon start %s — ",
        tu,
        horizon$start
      ),
      "all panel rows fall within the warmup window. Skipping this TU."
    )
    return(invisible(list(
      event_table = tibble::tibble(),
      monthly_summary = tibble::tibble()
    )))
  }

  # Wide panel: rows = months, columns = climate vars stacked per gid.
  # Cheaper: index per (gid, month) directly via a key. Warm-up rows get
  # `.m_idx = NA` (their t_start is below the trimmed grid) and are dropped.
  panel_key <- panel_tu %>%
    dplyr::mutate(
      .gid_idx = match(as.character(gid), frame_tu$gid),
      .m_idx = match(t_start, .env$months)
    ) %>%
    dplyr::filter(!is.na(.gid_idx), !is.na(.m_idx))

  build_clim <- function(col, default = 0) {
    M_mat <- matrix(default, nrow = N, ncol = M)
    M_mat[cbind(panel_key$.gid_idx, panel_key$.m_idx)] <- panel_key[[col]]
    M_mat
  }
  climate <- list(
    mean_precip = build_clim("mean_precip"),
    dev_precip = build_clim("dev_precip"),
    mean_temp_daymean = build_clim("mean_temp_daymean"),
    dev_temp_daymean = build_clim("dev_temp_daymean")
  )
  # ENSO (`enso_lag<n>`) covariate(s), if present in the panel — fed to the
  # hazard LP under the same name as the brms FE column.
  # enso_lag<n> and enso_lead<k> (deployed ENSO-lead fits; base.R::enso_lag_col_name()).
  enso_cols <- grep("^enso_(lag|lead)[0-9]+$", names(panel_tu), value = TRUE)
  for (ec in enso_cols) {
    climate[[ec]] <- build_clim(ec)
  }
  # ENSO is a global (location-invariant) index: all gids share the same value
  # in a given month. build_clim() defaults absent gid-months to 0, so gids
  # whose panel data ends before the simulation horizon (e.g. countries with no
  # recent observed outbreaks) silently receive ENSO = 0. With a multi-year
  # contpanel warmup this accumulates into systematically wrong TTLO trajectories
  # for those countries. Fix: for each month, propagate the ENSO value from any
  # present gid to all N rows (take the element with largest |value| so the
  # correct non-zero value wins over the 0 default from absent gids).
  for (ec in enso_cols) {
    mat <- climate[[ec]]
    month_vals <- apply(mat, 2L, function(col) col[which.max(abs(col))])
    climate[[ec]] <- matrix(month_vals, nrow = N, ncol = M, byrow = TRUE)
  }
  # Location-specific climate variables (mean_precip, dev_precip,
  # mean_temp_daymean, dev_temp_daymean): gids whose panel coverage ends before
  # the simulation horizon receive 0 from build_clim(). Forward-fill each gid's
  # last observed value across trailing missing months so absent gids carry a
  # plausible climatological state rather than a spurious zero.
  # Note: current elag occurrence configs use only enso_lag<n>; this fix is
  # relevant for the duration/monthly models and future occurrence configs including these terms.
  {
    has_data <- matrix(FALSE, nrow = N, ncol = M)
    if (nrow(panel_key) > 0L) {
      has_data[cbind(panel_key$.gid_idx, panel_key$.m_idx)] <- TRUE
    }
    loc_clim_cols <- c(
      "mean_precip",
      "dev_precip",
      "mean_temp_daymean",
      "dev_temp_daymean"
    )
    for (lcc in loc_clim_cols) {
      if (!(lcc %in% names(climate))) {
        next
      }
      mat <- climate[[lcc]]
      for (i in seq_len(N)) {
        real_ms <- which(has_data[i, ])
        if (length(real_ms) == 0L || max(real_ms) == M) {
          next
        }
        last_m <- max(real_ms)
        mat[i, (last_m + 1L):M] <- mat[i, last_m]
      }
      climate[[lcc]] <- mat
    }
  }

  # DLNM cross-basis design --------------------------------------------------
  # Reconstruct `<var>_cb<k>` covariate matrices for any cross-basis fit, driven
  # entirely by the cb_meta sidecars (self-describing: argvar/arglag/lag/
  # cb_col_names). Each track (occ / dur / mon) is namespaced separately because
  # the fits may carry different bases for the same variable name. Two regimes:
  #   * climate vars  — replicate-invariant; built once here (below).
  #   * "enso"        — replicate-variant (the enso_lag* matrices may be swapped
  #                     per replicate by enso_scenarios); rebuilt in the loop.
  # `cb_meta` entry "enso" routes make_lag_matrix() to the enso_lag* branch; any
  # other name reads the climate `<var>_lag*` columns (emitted by 03c/03d).
  split_cb <- function(cb_meta) {
    if (is.null(cb_meta) || length(cb_meta) == 0) {
      return(list(enso = NULL, clim = NULL))
    }
    enso <- cb_meta[intersect("enso", names(cb_meta))]
    clim <- cb_meta[setdiff(names(cb_meta), "enso")]
    list(
      enso = if (length(enso)) enso else NULL,
      clim = if (length(clim)) clim else NULL
    )
  }
  # Build cross-basis N×M matrices from a set of lag matrices (one per lag,
  # lag0..max_lag) flattened column-major, then reshaped back. Returns a named
  # list keyed by cb_col_names. `lag_getter(col)` returns an N×M matrix for the
  # named lag column (build_clim for climate, climate[[..]] for enso).
  build_cb_mats <- function(
    meta_list,
    lag_getter,
    loc_invariant = FALSE,
    clamp = FALSE,
    clamp_margin = 0
  ) {
    mats <- list()
    if (is.null(meta_list)) {
      return(mats)
    }
    for (v in names(meta_list)) {
      # Latent cross-basis variables (utils/dlnm.R) are not panel columns and so
      # have no N×M matrix: `mso` differs per particle and is only known inside
      # the month loop. The kernel builds their basis with a direct
      # build_cb_columns_from_Q() lookup table instead (see B_mso below).
      if (v %in% .LATENT_CB_VARS) {
        next
      }
      meta <- meta_list[[v]]
      if (!is.null(meta$single_lag)) {
        sl <- as.integer(meta$single_lag)
        lag_cols <- if (identical(v, "enso")) {
          # negative = enso_lead<k>
          enso_lag_col_name(sl)
        } else if (sl == 0L) {
          v
        } else {
          sprintf("%s_lag%d", v, sl)
        }
        n_lag_cols <- 1L
      } else if (!is.null(meta$lags)) {
        stop(
          "build_cb_mats: multi-lag cross-basis meta (`lags`, brms_kl) is not supported ",
          "by the simulator kernel (variable '",
          v,
          "').",
          call. = FALSE
        )
      } else {
        max_lag <- as.integer(meta$lag[2])
        lag_cols <- if (identical(v, "enso")) {
          sprintf("enso_lag%d", seq(0L, max_lag))
        } else {
          c(v, sprintf("%s_lag%d", v, seq_len(max_lag)))
        }
        n_lag_cols <- max_lag + 1L
      }
      lag_mats <- lapply(lag_cols, lag_getter)
      if (loc_invariant) {
        # ENSO is location-invariant: transform M rows then broadcast.
        Q <- matrix(0, nrow = M, ncol = n_lag_cols)
        for (k in seq_along(lag_mats)) {
          Q[, k] <- lag_mats[[k]][1L, ]
        }
        Q[is.na(Q)] <- 0
        # simulator.enso.clamp_to_knots: hold ENSO beyond the fitted boundary
        # knots widened by simulator.enso.clamp_margin (.clamp_cb_range(), dlnm.R).
        if (clamp) {
          Q <- .clamp_cb_range(Q, meta$range, clamp_margin)
        }
        cb <- build_cb_columns_from_Q(Q, meta) # M × ncb
        for (j in seq_len(ncol(cb))) {
          mats[[meta$cb_col_names[j]]] <-
            matrix(cb[, j], nrow = N, ncol = M, byrow = TRUE)
        }
      } else {
        Q <- matrix(0, nrow = N * M, ncol = n_lag_cols)
        for (k in seq_along(lag_mats)) {
          Q[, k] <- as.numeric(lag_mats[[k]])
        }
        Q[is.na(Q)] <- 0
        cb <- build_cb_columns_from_Q(Q, meta) # (N*M) × ncb
        for (j in seq_len(ncol(cb))) {
          mats[[meta$cb_col_names[j]]] <- matrix(cb[, j], nrow = N, ncol = M)
        }
      }
    }
    mats
  }

  # LP contexts. One per group under grouping; a single context
  # spanning all rows otherwise. Each context owns its draws bundles, cross-basis
  # split metadata, and replicate-invariant climate cross-basis matrices (built
  # once here, indexed by the context's rows at LP time). ENSO cross-basis is
  # replicate-variant and rebuilt inside the loop per context. The single-context
  # case reproduces the pre-grouping `clim_cb_occ` + `occ_split`
  # exactly, so the math is unchanged when `grouping` is NULL.
  make_context <- function(gname) {
    idx <- if (is_grouped) which(grp_vec == gname) else seq_len(N)
    if (is_grouped) {
      d_occ <- grouping$draws_occ_by_group[[gname]]
      cbm_occ <- grouping$cb_meta_occ_by_group[[gname]]
      if (is.null(d_occ)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing occ draws for group '%s'.",
          gname
        ))
      }
    } else {
      d_occ <- draws_occ
      cbm_occ <- cb_meta_occ
    }
    o_split <- split_cb(cbm_occ)
    # Duration draws AND cross-basis, routed per group exactly like occ.
    # `04_fit_dur_brms_cb.R` fits carry `enso_cb<k>` fixed effects, and the
    # duration track's best lag is chosen independently of the other tracks
    # (Horn: occ=11, siz=9, mon=11, dur=8), so the duration basis must be built
    # from its OWN cb_meta -- reusing the occurrence columns would evaluate the
    # duration coefficients on the wrong lag and the wrong spline basis.
    # `d_split` is NULL when the duration fit carries no cross-basis, which
    # makes build_cb_mats() a no-op.
    d_dur <- NULL
    d_split <- NULL
    if (has_dur) {
      d_dur <- if (is_grouped) {
        grouping$draws_dur_by_group[[gname]]
      } else {
        draws_dur
      }
      if (is_grouped && is.null(d_dur)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing duration draws for group '%s'.",
          gname
        ))
      }
      cbm_dur <- if (is_grouped) {
        grouping$cb_meta_dur_by_group[[gname]]
      } else {
        cb_meta_dur
      }
      d_split <- split_cb(cbm_dur)
    }
    # Monthly bundle + its cross-basis, routed per group exactly like
    # occ. `mso` is stripped out of `m_split$clim` by build_cb_mats()'s
    # .LATENT_CB_VARS skip and handled by the B_mso lookup table instead.
    d_mon <- NULL
    m_split <- NULL
    mso_tab <- NULL
    if (has_mon) {
      d_mon <- if (is_grouped) {
        grouping$draws_mon_by_group[[gname]]
      } else {
        draws_mon
      }
      cbm_mon <- if (is_grouped) {
        grouping$cb_meta_mon_by_group[[gname]]
      } else {
        cb_meta_mon
      }
      if (is.null(d_mon)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing monthly draws for group '%s'.",
          gname
        ))
      }
      m_split <- split_cb(cbm_mon)
      # NULL when the fit has no mso basis; else $B is (mso_max + 1) x K, rows
      # indexed by mso = 0..mso_max (read by .inject_mso()).
      mso_tab <- build_mso_basis(cbm_mon)
    }
    list(
      gname = gname,
      idx = idx,
      frame_g = frame_tu[idx, , drop = FALSE],
      draws_occ = d_occ,
      draws_dur = d_dur,
      draws_mon = d_mon,
      drop_region = is_grouped,
      occ_split = o_split,
      dur_split = d_split,
      mon_split = m_split,
      mso_tab = mso_tab,
      clim_cb_occ = build_cb_mats(
        o_split$clim,
        build_clim,
        loc_invariant = FALSE
      ),
      clim_cb_dur = if (has_dur) {
        build_cb_mats(d_split$clim, build_clim, loc_invariant = FALSE)
      } else {
        NULL
      },
      clim_cb_mon = if (has_mon) {
        build_cb_mats(m_split$clim, build_clim, loc_invariant = FALSE)
      } else {
        NULL
      }
    )
  }
  contexts <- lapply(ctx_groups, make_context)
  names(contexts) <- ctx_groups

  # t_start_seas is purely deterministic: yday(t_start shifted by min_seas) - 1.
  # Reading it from the panel (build_clim) defaults to 0 for gid-months that
  # lack a panel row (e.g. gids whose data starts mid-horizon), which forces the
  # seasonal spline to evaluate at the season-start value for all those months
  # and washes out the seasonal signal in the aggregate. Instead, compute it
  # from frame_tu$min_seas and the months grid — always correct.
  # Using lapply over gids (rows) to avoid outer() coercing the Date vector to
  # numeric, which would require an origin argument for as.Date().
  t_start_seas_mat <- do.call(
    rbind,
    lapply(as.integer(frame_tu$min_seas), function(ms) {
      as.integer(
        lubridate::yday(months %m+% lubridate::period(-(ms - 1L), "months")) -
          1L
      )
    })
  ) # N × M

  # Exposure offset — mirrors `offset(log_offset)` in 04b_run_occ_brms.R:114:
  #   exposure = pmax(t_end_seas - t_start_seas, 1L)
  #   log_offset = log(exposure)
  # brms parameters encode log-rate per day; the probability kernel needs
  # lambda = exp(eta) * exposure_days, not exp(eta) * dt_months.
  # Exposure = days_in_month(t_start) — the same for all gids in a given month.
  log_offset_mat <- matrix(
    log(as.integer(lubridate::days_in_month(months))),
    nrow = N,
    ncol = M,
    byrow = TRUE
  )

  # Shared warmup state: computed once when warmup_panels_list is NULL
  # (backward-compatible default); otherwise deferred to per-replicate calls inside
  # the loop so each replicate can start from a distinct importance-sampled
  # simulation trajectory.
  # TTLO level "none": seeded only when the occurrence fit
  # carries it; NULL otherwise.
  npo <- .ttlo_no_prior_for_frame(
    frame_tu,
    horizon$start,
    .draws_has_ttlo_none(
      if (!is.null(grouping)) grouping$draws_occ_by_group else draws_occ
    )
  )
  shared_warm <- if (is.null(warmup_panels_list)) {
    build_warmup_state_legacy(
      panel_tu,
      frame_tu,
      horizon,
      has_dur = has_dur,
      no_prior_onset = npo
    )
  } else {
    NULL
  }

  # Output accumulators
  event_rows <- vector("list", R)
  monthly_Y <- array(0L, dim = c(R, N, M))
  # Per-(replicate, location, month) "is this location in an active
  # episode this month" indicator. NULL on the legacy path (has_dur = FALSE)
  # so it costs nothing and never enters monthly_summary's column set.
  monthly_A <- if (has_dur) array(0L, dim = c(R, N, M)) else NULL
  # covariate_feed = "detected": the detected-activity field per month (the
  # `det` state after update_detection()); NULL otherwise, so latent runs keep
  # their monthly_summary column set.
  monthly_D <- if (use_det) array(0L, dim = c(R, N, M)) else NULL
  # Monthly case-count (mon) track: mu and the NegBin case draw at every ACTIVE
  # month. NULL on the legacy path for the
  # same reason as monthly_A. Left NA (not 0) where a location is inactive, so
  # "no episode" is distinguishable from "active episode that drew zero cases".
  monthly_MU <- if (has_mon) array(NA_real_, dim = c(R, N, M)) else NULL
  # Per-location monthly incidence ceiling (cap x pop), frame_tu row order;
  # Inf everywhere when incidence_cap is NULL (direct callers, tests).
  mu_cap <- if (has_mon) monthly_mu_cap(frame_tu$pop, incidence_cap) else NULL
  # pop_year mode: the U x M per-month population and the
  # per-month ceiling. NULL for a frame without pop_y<YYYY> columns, which
  # keeps the static `pop` / `mu_cap` above (the default, unchanged).
  pop_mat <- if (has_mon) frame_pop_matrix(frame_tu, months) else NULL
  mu_cap_mat <- if (!is.null(pop_mat)) {
    matrix(
      monthly_mu_cap(as.vector(pop_mat), incidence_cap),
      nrow = nrow(pop_mat)
    )
  } else {
    NULL
  }
  monthly_C <- if (has_mon) array(NA_integer_, dim = c(R, N, M)) else NULL

  # Use L'Ecuyer-CMRG streams seeded from (seed, replicate_id) — independent of
  # which array slot a replicate lands in (S7).
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1]), add = TRUE)

  # Snapshot the baseline ENSO covariate matrices BEFORE the replicate loop.
  # The per-replicate override overwrites `climate[[ec]]` in place;
  # without this snapshot the mutation would persist into the next replicate,
  # so a replicate that does not receive an override (e.g. rep_id beyond the
  # available scenario columns) would silently reuse the previous replicate's
  # ENSO trajectory instead of the baseline panel.
  enso_base <- if (!is.null(enso_scenarios)) {
    climate[intersect(names(enso_scenarios), names(climate))]
  } else {
    NULL
  }

  # Neighbour-contagion kernel (within-partition / within-country). Built ONCE
  # per TU, replicate-invariant, only when `use_nbr` (flag on AND the occurrence
  # fit carries an nbr_* FE, resolved by the caller). `K_adm1_nbr` is a dense
  # n_adm1 x n_adm1 matrix in `unique(frame_tu$adm1_unit)` order with zero
  # diagonal; each step reconstructs the lag-1 ADM1-active vector from the `o`
  # state and forms `K %*% Aa1`. use_nbr=FALSE ⇒ NULL ⇒ every nbr block below is
  # a no-op (nbr stays NULL into evaluate_lp).
  K_adm1_nbr <- if (isTRUE(use_nbr)) {
    # Spec resolved by the driver (resolve_nbr_spec()); legacy callers that pass
    # only the scalar nbr_* args get the same resolution here, so the kernel's
    # scale/partition always follow the fit's own nbr FE name.
    nbr_spec <- nbr_spec %||%
      resolve_nbr_spec(
        list(
          enable = TRUE,
          rho_km = nbr_rho_km,
          gadm_path = nbr_gadm_path,
          cache_dir = nbr_cache_dir
        ),
        draws_occ,
        grouping
      )
    nbr_kernel_for_frame(frame_tu, nbr_spec)
  } else {
    NULL
  }

  for (r_idx in seq_len(R)) {
    rep_id <- replicate_ids[r_idx]
    set.seed(seed + rep_id)

    # Per-context posterior draw indices (each group bundle has its own S; clamp
    # the shared raw index to that group's draw count). Single-context case
    # reproduces posterior_idx_occ[r_idx] exactly. Also scatter to
    # length-N per-row vectors for the event-table `posterior_draw_*` columns.
    s_occ_byctx <- vapply(
      contexts,
      function(ctx) {
        as.integer((raw_idx_occ[r_idx] - 1L) %% ctx$draws_occ$S) + 1L
      },
      integer(1)
    )
    s_occ_byrow <- integer(N)
    for (ctx in contexts) {
      s_occ_byrow[ctx$idx] <- s_occ_byctx[[ctx$gname]]
    }
    # Duration draw index per context/row — only meaningful when has_dur.
    s_dur_byctx <- NULL
    s_dur_byrow <- integer(N)
    if (has_dur) {
      s_dur_byctx <- vapply(
        contexts,
        function(ctx) {
          as.integer((raw_idx_dur[r_idx] - 1L) %% ctx$draws_dur$S) + 1L
        },
        integer(1)
      )
      for (ctx in contexts) {
        s_dur_byrow[ctx$idx] <- s_dur_byctx[[ctx$gname]]
      }
    }
    # Monthly draw index per context — only meaningful when has_mon.
    s_mon_byctx <- NULL
    if (has_mon) {
      s_mon_byctx <- vapply(
        contexts,
        function(ctx) {
          as.integer((raw_idx_mon[r_idx] - 1L) %% ctx$draws_mon$S) + 1L
        },
        integer(1)
      )
    }

    # Per-replicate warmup state.
    # When warmup_panels_list is provided, each replicate draws its initial AR
    # state from the importance-sampled mini-panel for that global replicate ID.
    # Indexed by rep_id (global), not r_idx (within-chunk position), so chunks
    # always draw from the correct element of the full n_rep list.
    warm <- if (!is.null(warmup_panels_list)) {
      build_warmup_state_legacy(
        warmup_panels_list[[rep_id]],
        frame_tu,
        horizon,
        has_dur = has_dur,
        no_prior_onset = npo
      )
    } else {
      shared_warm
    }

    o <- warm$o
    o_ct <- warm$o_ct
    o_tu <- warm$o_tu
    if (isTRUE(drop_country_tu_spillover)) {
      o_ct[!is.na(o_ct)] <- 0L
      o_tu[!is.na(o_tu)] <- 0L
    }
    u_days <- warm$u # numeric days vector (step 4)

    # Duration-aware state (NULL/unused on the legacy path). `active`/
    # `mso` seed from build_warmup_state_legacy() (currently always 0 — see that
    # function's doc comment on the smoothing-cache handoff).
    # `open_row` carries the row index
    # (within this replicate's event_rows[[r_idx]]) of each active location's
    # not-yet-closed event row, so cessation/censoring can patch it in place.
    active <- if (has_dur) as.integer(warm$active) else NULL
    mso <- if (has_dur) as.integer(warm$mso) else NULL
    det <- if (use_det) as.integer(warm$det %||% warm$active) else NULL
    open_row <- if (has_dur) rep(NA_integer_, N) else NULL

    # Per-replicate ENSO override.
    # Replace the shared panel's ENSO covariate matrices with the replicate-
    # specific forecast trajectory from `enso_scenarios`. ENSO is location-
    # invariant, so the M-length vector is broadcast to all N rows. Indexing
    # uses rep_id (global replicate ID) so each chunk accesses the correct column
    # of the [M × n_rep] matrix regardless of its within-chunk position.
    if (!is.null(enso_scenarios)) {
      for (ec in names(enso_scenarios)) {
        if (!(ec %in% names(climate))) {
          next
        }
        # Restore the baseline ENSO matrix first so a previous replicate's
        # override never leaks into this one.
        climate[[ec]] <- enso_base[[ec]]
        if (rep_id > ncol(enso_scenarios[[ec]])) {
          stop(sprintf(
            "simulate_one_tu: replicate id %d exceeds available ENSO scenario columns (%d) for '%s' — n_replicates and the ENSO scenario matrix are out of sync.",
            rep_id,
            ncol(enso_scenarios[[ec]]),
            ec
          ))
        }
        scenario_vec <- enso_scenarios[[ec]][, rep_id]
        # Row 1 is the horizon start; the scenario runs to the end of the shared
        # panel, which can be past this TU's horizon. Keep the first M months:
        # recycling a longer vector row-wise would give each unit a different
        # (time-shifted) slice.
        if (length(scenario_vec) < M) {
          stop(sprintf(
            "simulate_one_tu: ENSO scenario '%s' has %d months but TU '%s' simulates %d.",
            ec,
            length(scenario_vec),
            tu,
            M
          ))
        }
        scenario_vec <- scenario_vec[seq_len(M)]
        climate[[ec]] <- matrix(scenario_vec, nrow = N, ncol = M, byrow = TRUE)
      }
    }

    # ENSO cross-basis: rebuilt per replicate from the (possibly scenario-swapped)
    # enso_lag* matrices. Location-invariant, so transform once per month then
    # broadcast. Combined with each context's replicate-invariant climate
    # cross-basis. Under grouping, every group may carry its own basis/max_lag,
    # so the cb columns are namespaced per context.
    enso_getter <- function(col) climate[[col]]
    cb_mats_occ_byctx <- lapply(contexts, function(ctx) {
      c(
        ctx$clim_cb_occ,
        build_cb_mats(
          ctx$occ_split$enso,
          enso_getter,
          loc_invariant = TRUE,
          clamp = clamp_enso_to_knots,
          clamp_margin = clamp_enso_margin
        )
      )
    })
    cb_mats_mon_byctx <- if (!has_mon) {
      NULL
    } else {
      lapply(contexts, function(ctx) {
        c(
          ctx$clim_cb_mon,
          build_cb_mats(
            ctx$mon_split$enso,
            enso_getter,
            loc_invariant = TRUE,
            clamp = clamp_enso_to_knots,
            clamp_margin = clamp_enso_margin
          )
        )
      })
    }
    cb_mats_dur_byctx <- if (!has_dur) {
      NULL
    } else {
      lapply(contexts, function(ctx) {
        c(
          ctx$clim_cb_dur,
          build_cb_mats(
            ctx$dur_split$enso,
            enso_getter,
            loc_invariant = TRUE,
            clamp = clamp_enso_to_knots,
            clamp_margin = clamp_enso_margin
          )
        )
      })
    }

    for (m in seq_len(M)) {
      # One named covariate vector per climate/enso column for this month
      # (full-N; shared across contexts since the climate panel is TU-wide).
      xrow <- lapply(climate, function(M_mat) M_mat[, m])
      xrow[["t_start_seas"]] <- t_start_seas_mat[, m]

      # Neighbour-contagion pressure for this step. `o` entering the step IS the
      # lag-1 ADM1-active indicator broadcast to each location (seeded from the
      # warmup on step 1, carried from the previous step's update_lags otherwise,
      # and already keyed on `active` not `Y` on the duration path). Collapse it
      # back to one value per ADM1 (in kernel order), apply the kernel, broadcast
      # to locations by adm1_unit (ADM2 inherits its parent). NULL when !use_nbr.
      nbr_loc <- NULL
      if (isTRUE(use_nbr)) {
        Aa1 <- tapply(as.integer(o), frame_tu$adm1_unit, max, default = 0L)
        Aa1 <- as.numeric(Aa1[rownames(K_adm1_nbr)])
        Aa1[is.na(Aa1)] <- 0
        nbr_adm1 <- nbr_from_adm1_active(Aa1, K_adm1_nbr)
        names(nbr_adm1) <- rownames(K_adm1_nbr)
        nbr_loc <- as.numeric(nbr_adm1[as.character(frame_tu$adm1_unit)])
        nbr_loc[is.na(nbr_loc)] <- 0
      }

      # Occurrence LP, partitioned by context: each group's rows are scored with
      # that group's draws bundle, cross-basis columns, posterior draw and
      # region-dropped RE key space, then scattered into the length-N log_lambda.
      log_lambda <- numeric(N)
      for (ctx in contexts) {
        idx <- ctx$idx
        xo <- lapply(xrow, function(v) v[idx])
        cbo <- cb_mats_occ_byctx[[ctx$gname]]
        for (cn in names(cbo)) {
          xo[[cn]] <- cbo[[cn]][idx, m]
        }
        bh_g <- basehaz_eval(
          ctx$draws_occ,
          t_start_seas_mat[idx, m],
          length(idx)
        )
        log_lambda[idx] <- evaluate_lp_legacy(
          ctx$draws_occ,
          s_occ_byctx[[ctx$gname]],
          xo,
          o[idx],
          o_ct[idx],
          o_tu[idx],
          u_days[idx],
          ctx$frame_g,
          basehaz_val = bh_g,
          drop_region = ctx$drop_region,
          drop_ct_tu = drop_country_tu_spillover,
          nbr = if (is.null(nbr_loc)) NULL else nbr_loc[idx]
        ) +
          # Empirical baseline-hazard anchor, offset(log(base_rate_occ)). Its
          # denominator is exposure days -- exactly the log_offset_mat term
          # added below -- so with the anchor on, exp(Intercept) is a residual
          # near 1 rather than a whole-corpus hazard scale. Identically 0 for
          # any fit made without it.
          .anchor_log_offset(ctx$draws_occ, frame_tu, idx, "occurrence kernel")
      }
      # Per-step probability: 1 - exp(-λ) via -expm1(-λ).
      # λ = rate_per_day × exposure_days; the offset (log exposure) is added to
      # the LP rather than multiplying exp(eta) by dt_months, because the brms
      # parameters are on the per-day rate scale (offset(log_offset) in the fit).
      lambda <- exp(log_lambda + log_offset_mat[, m])
      p <- -expm1(-lambda)

      # Onset draw. Legacy path: unconditional per-step Bernoulli
      # (byte-identical RNG stream). Duration-aware path: gated on
      # not-already-active (`at_risk`) so a location that is mid-episode
      # cannot "re-onset" every month — see the @param draws_dur doc above.
      # Hard branch,
      # not a masked single call, so the legacy RNG stream is untouched.
      active_before <- if (has_dur) active else NULL
      if (has_dur) {
        at_risk <- active_before == 0L
        Y <- integer(N)
        if (any(at_risk)) {
          Y[at_risk] <- as.integer(stats::rbinom(sum(at_risk), 1L, p[at_risk]))
        }
      } else {
        Y <- as.integer(stats::rbinom(N, 1L, p))
      }

      monthly_Y[r_idx, , m] <- Y

      if (has_dur) {
        # active/mso bookkeeping. new_onset locations start mso = 1 this
        # month; already-active locations (active_before == 1, necessarily
        # Y == 0 there since at_risk excluded them) advance mso by one.
        new_onset <- which(Y == 1L)
        still_active <- which(active_before == 1L)
        mso[still_active] <- mso[still_active] + 1L
        mso[new_onset] <- 1L
        active[new_onset] <- 1L
        monthly_A[r_idx, , m] <- active
      }

      if (any(Y == 1L)) {
        # Event rows
        idx_trigger <- which(Y == 1L)
        if (has_dur) {
          # duration_months/end_month/censored are placeholders here, patched
          # in place (by row index, via `open_row`) at cessation or at the
          # end-of-horizon censoring pass below — the row is created at onset,
          # before the episode's length is known.
          nrow_before <- if (is.null(event_rows[[r_idx]])) {
            0L
          } else {
            nrow(event_rows[[r_idx]])
          }
          event_rows[[r_idx]] <- rbind(
            event_rows[[r_idx]],
            data.frame(
              run_tag = NA_character_, # filled at write-time by 09_simulator
              replicate_id = rep_id,
              gid = frame_tu$gid[idx_trigger],
              country = frame_tu$country[idx_trigger],
              tu = frame_tu$tu[idx_trigger],
              start_month = as.character(months[m]),
              duration_months = NA_integer_, # filled at cessation/censoring
              end_month = NA_character_, # filled at cessation/censoring
              censored = NA_integer_, # filled at cessation/censoring
              posterior_draw_o = s_occ_byrow[idx_trigger],
              posterior_draw_d = s_dur_byrow[idx_trigger],
              stringsAsFactors = FALSE
            )
          )
          open_row[idx_trigger] <- nrow_before + seq_along(idx_trigger)
        } else {
          event_rows[[r_idx]] <- rbind(
            event_rows[[r_idx]],
            data.frame(
              run_tag = NA_character_, # filled at write-time by 09_simulator
              replicate_id = rep_id,
              gid = frame_tu$gid[idx_trigger],
              country = frame_tu$country[idx_trigger],
              tu = frame_tu$tu[idx_trigger],
              start_month = as.character(months[m]),
              duration_months = 1L, # legacy path: no duration model
              posterior_draw_o = s_occ_byrow[idx_trigger],
              stringsAsFactors = FALSE
            )
          )
        }
      }

      # Monthly case-count draw (`mon` track).
      #
      # Placement is load-bearing on both sides:
      #   * AFTER the promotion block, so a location onsetting this month is
      #     already active with mso = 1 and emits a first month of cases;
      #   * BEFORE update_lags_legacy() and the cessation draw, so (a) the spillover
      #     covariates carry the same lag-1 values the occurrence LP used
      #     this month, matching how the fit sees them, and (b) a location
      #     ceasing this month still emits its final month.
      #
      # The offsets are read off the bundle (`offset_pop`, `offset_logoffset =
      # TRUE` pinned by extract_monthly_draws(), and the optional base_rate_mon
      # anchor). mu is gated on `active`, not on `Y`.
      if (has_mon && any(active == 1L)) {
        log_mu_mon <- numeric(N)
        phi_mon_vec <- numeric(N)

        # Latent mso cross-basis: one row of the per-TU lookup table per
        # location, clamped into the fitted boundary knots (build_mso_basis()).
        # Inactive locations (mso = 0) are scored too -- the result is simply
        # discarded by the `active == 1` gate below -- so the vectorised index
        # never sees an out-of-range value.
        xrow_mon <- xrow
        xrow_mon[["month_num"]] <- rep(
          as.integer(lubridate::month(months[m])),
          N
        )
        # monthly.seasonal_term: "harmonic2" (the default) names these four
        # columns; load_monthly_data() (utils/mon_data.R:92-99) builds them the
        # same way at fit time. Harmless when the fit used seasonal_term: none
        # -- evaluate_lp_legacy() drops any term absent from `fe`.
        mnum <- as.numeric(lubridate::month(months[m]))
        xrow_mon[["sin_month"]] <- rep(sin(2 * pi * mnum / 12), N)
        xrow_mon[["cos_month"]] <- rep(cos(2 * pi * mnum / 12), N)
        xrow_mon[["sin2_month"]] <- rep(sin(4 * pi * mnum / 12), N)
        xrow_mon[["cos2_month"]] <- rep(cos(4 * pi * mnum / 12), N)
        # monthly.admin_level: "ADM1-ADM2" makes admin_level a 2-level factor at
        # fit time; brms names the dummy "admin_levelADM2". Absent from the
        # frame (or from `fe`) => contributes nothing. Keyed on the CANONICAL
        # level, because the factor the fit levelled was the P3-re-keyed one
        # (.frame_admin_level_canonical()).
        if ("admin_level" %in% names(frame_tu)) {
          xrow_mon[["admin_levelADM2"]] <-
            as.numeric(.frame_admin_level_canonical(frame_tu) == "ADM2")
        }

        for (ctx in contexts) {
          idx <- ctx$idx
          xm <- lapply(xrow_mon, function(v) v[idx])
          cbm <- cb_mats_mon_byctx[[ctx$gname]]
          for (cn in names(cbm)) {
            xm[[cn]] <- cbm[[cn]][idx, m]
          }
          # Offsets are read off the bundle (monthlyLogMu_legacy()). In forecast mode
          # there is no observed `exposure` column, so the monthly exposure is
          # synthesised as log(days_in_month(t)) -- the same log_offset_mat the
          # occurrence and duration kernels already use, and the same quantity
          # the fit's panel carries for an INTERIOR month. The ONSET month is
          # partial in reality and is corrected to half a month (mso was
          # promoted to 1 above, at the new_onset block, so this reads the right
          # value). The CESSATION month is still given full-month exposure: it
          # is not knowable here, since the cessation hazard is drawn below
          # rather than above. See .ONSET_EXPOSURE_FRAC for the measured residual.
          lm_g <- monthlyLogMu_legacy(
            ctx$draws_mon,
            s_mon_byctx[[ctx$gname]],
            xm,
            o[idx],
            o_ct[idx],
            o_tu[idx],
            u_days[idx],
            ctx$frame_g,
            log_offset = log_offset_mat[idx, m],
            mso = mso[idx],
            mso_tab = ctx$mso_tab,
            what = "monthly kernel",
            drop_region = ctx$drop_region,
            drop_ct_tu = drop_country_tu_spillover,
            pop = if (!is.null(pop_mat)) pop_mat[idx, m] else NULL
          )
          log_mu_mon[idx] <- lm_g
          phi_mon_vec[idx] <- compute_phi_vec(
            ctx$draws_mon,
            s_mon_byctx[[ctx$gname]],
            ctx$frame_g,
            drop_region = ctx$drop_region
          )
        }

        # Monthly incidence cap (simulator.incidence_cap): the mean, and below
        # the drawn count, are held at cap x pop. Inf when incidence_cap is NULL.
        if (!is.null(mu_cap_mat)) {
          mu_cap <- mu_cap_mat[, m]
        }
        mu_mon <- exp(log_mu_mon)
        if (any(is.finite(mu_cap))) {
          mu_mon <- pmin(mu_mon, mu_cap)
        }
        act_i <- which(active == 1L)
        monthly_MU[r_idx, act_i, m] <- mu_mon[act_i]
        cm <- draw_monthly_cases(
          mu_mon[act_i],
          phi_mon_vec[act_i],
          mu_cap[act_i]
        )
        ok <- !is.na(cm)
        if (any(ok)) monthly_C[r_idx, act_i[ok], m] <- cm[ok]
      }

      # Update state for next step. Duration-aware path feeds `active`
      # (ongoing-episode indicator) instead of `Y` (onset-only indicator) to
      # update_lags_legacy(), so the contagion covariates reflect the true duty
      # cycle of an outbreak rather than firing for exactly the onset month
      # — see the @param draws_dur doc above ("the contagion-covariate
      # duty-cycle bug"). advance_ttlo_legacy() stays keyed on Y on both paths: TTLO
      # resets at onset regardless of how long the episode then runs.
      # covariate_feed = "detected": spillover covariates and TTLO follow the
      # detected-activity field (update_detection()), not the latent one.
      if (use_det) {
        det <- update_detection(Y, active, det, delta_fn, delta_fp)
        monthly_D[r_idx, , m] <- det
      }
      lags <- update_lags_legacy(
        if (use_det) {
          det
        } else if (has_dur) {
          active
        } else {
          Y
        },
        frame_tu
      )
      if (isTRUE(force_o_tu)) {
        lags$o_tu <- ifelse(is.na(lags$o_tu), NA_integer_, 1L)
      }
      # drop_country_tu_spillover takes PRECEDENCE over force_o_tu: o_ct/o_tu
      # are forced to 0 (skipping the leave-one-out aggregation entirely),
      # not 1, when both flags are set.
      if (isTRUE(drop_country_tu_spillover)) {
        lags$o_ct <- ifelse(is.na(lags$o_ct), NA_integer_, 0L)
        lags$o_tu <- ifelse(is.na(lags$o_tu), NA_integer_, 0L)
      }
      o <- lags$o
      o_ct <- lags$o_ct
      o_tu <- lags$o_tu
      u_days <- advance_ttlo_legacy(u_days, if (use_det) Y * det else Y)

      if (has_dur) {
        # Cessation hazard, evaluated for every context-partitioned row (like
        # the occurrence LP above) but only drawn/applied for currently
        # active locations. xrow_dur reuses this month's climate covariates;
        # mso feeds the primary months-since-onset baseline hazard via
        # basehaz_dur_eval_legacy().
        log_lambda_dur <- numeric(N)
        xrow_dur <- xrow
        for (ctx in contexts) {
          idx <- ctx$idx
          xd <- lapply(xrow_dur, function(v) v[idx])
          # Duration cross-basis columns, exactly as occ/mon above. The
          # deployed per-group best-lag duration fits are `_cb` fits carrying
          # `enso_cb<k>` fixed effects; without this injection evaluate_lp_legacy()
          # silently skips them (absent covariate) and the track is ENSO-blind.
          # NULL/empty for a duration fit with no cross-basis.
          cbd <- cb_mats_dur_byctx[[ctx$gname]]
          for (cn in names(cbd)) {
            xd[[cn]] <- cbd[[cn]][idx, m]
          }
          bh_g <- basehaz_dur_eval_legacy(
            ctx$draws_dur,
            s_dur_byctx[[ctx$gname]],
            mso[idx],
            length(idx)
          )
          log_lambda_dur[idx] <- evaluate_lp_legacy(
            ctx$draws_dur,
            s_dur_byctx[[ctx$gname]],
            xd,
            o[idx],
            o_ct[idx],
            o_tu[idx],
            u_days[idx],
            ctx$frame_g,
            basehaz_val = bh_g,
            drop_region = ctx$drop_region,
            drop_ct_tu = drop_country_tu_spillover
          ) +
            # offset(log(base_rate_dur)) -- denominator is exposure days, the
            # same log_offset_mat term added below. 0 for un-anchored fits.
            .anchor_log_offset(ctx$draws_dur, frame_tu, idx, "duration kernel")
        }
        # Same Poisson-with-log-offset link as occurrence (04_fit_dur_brms.R /
        # build_duration_formula_parts(): offset(log_offset) on days-in-month
        # exposure, never optional for this track -- hence no bundle gate on the
        # onset correction, unlike the mon block above). Duration shares the mon
        # track's episode-overlap exposure (both panels come from
        # .expand_active_months()), so it carries the same onset defect: without
        # the correction the cessation hazard for an episode's FIRST month is
        # applied to ~31 days instead of ~15, making episodes die too early.
        lambda_dur <- exp(
          log_lambda_dur +
            log_offset_mat[, m] +
            .onset_log_offset_adj_legacy(mso)
        )
        p_dur <- -expm1(-lambda_dur)

        act_idx <- which(active == 1L)
        ended <- integer(N)
        if (length(act_idx) > 0L) {
          ended[act_idx] <- as.integer(stats::rbinom(
            length(act_idx),
            1L,
            p_dur[act_idx]
          ))
        }

        end_idx <- which(ended == 1L)
        if (length(end_idx) > 0L) {
          for (i in end_idx) {
            if (!is.na(open_row[i])) {
              rn <- open_row[i]
              event_rows[[r_idx]]$duration_months[rn] <- mso[i]
              event_rows[[r_idx]]$end_month[rn] <- as.character(months[m])
              event_rows[[r_idx]]$censored[rn] <- 0L
            }
          }
          active[end_idx] <- 0L
          mso[end_idx] <- 0L
          open_row[end_idx] <- NA_integer_
        }
      }
    }

    if (has_dur) {
      # End-of-horizon administrative censoring: any location still active
      # when the simulated grid runs out gets its still-open event row
      # patched as censored at the last simulated month, mirroring the
      # fit-time panel's own administrative-censoring convention
      # (analysis/utils/dur_data.R).
      still_open <- which(active == 1L & !is.na(open_row))
      if (length(still_open) > 0L) {
        for (i in still_open) {
          rn <- open_row[i]
          event_rows[[r_idx]]$duration_months[rn] <- mso[i]
          event_rows[[r_idx]]$end_month[rn] <- as.character(months[M])
          event_rows[[r_idx]]$censored[rn] <- 1L
        }
      }
    }
  }

  event_table <- dplyr::bind_rows(event_rows)

  # Monthly summary: collapse the arrays to long form
  rid <- rep(replicate_ids, each = N * M)
  gid_v <- rep(rep(frame_tu$gid, each = M), times = R)
  m_v <- rep(rep(as.character(months), times = N), times = R)
  # aperm reorders dim c(R, N, M) → c(M, N, R) so as.vector varies month fastest,
  # matching the rid/gid_v/m_v index construction (replicate slowest, month fastest).
  # Without this, as.vector varies R fastest (column-major), misaligning every
  # (replicate_id, gid, month) label with its Y value and flattening all seasonal curves.
  monthly_summary <- tibble::tibble(
    replicate_id = rid,
    gid = gid_v,
    month = m_v,
    Y = as.vector(aperm(monthly_Y, c(3L, 2L, 1L)))
  )
  # Optional channels appended in flag order, so the legacy (has_dur = FALSE)
  # and duration-only column sets are unchanged from the two-branch form this
  # replaced. Each uses the same aperm() reordering as the columns above.
  if (has_dur) {
    monthly_summary$active <- as.vector(aperm(monthly_A, c(3L, 2L, 1L)))
  }
  if (has_mon) {
    monthly_summary$mu_mon <- as.vector(aperm(monthly_MU, c(3L, 2L, 1L)))
    monthly_summary$cases_month <- as.vector(aperm(monthly_C, c(3L, 2L, 1L)))
  }
  if (use_det) {
    monthly_summary$detected <- as.vector(aperm(monthly_D, c(3L, 2L, 1L)))
  }

  if (!dry_run) {
    dir.create(output_part, recursive = TRUE, showWarnings = FALSE)
    ext <- sim_io_ext(out_fmt)
    sim_io_save(
      event_table,
      file.path(output_part, paste0("event_table", ext)),
      out_fmt
    )
    # Format 2 (utils/sim_io.R::sim_monthly_save()): the same table, stored
    # sparsely; every reader gets it back identical through sim_io_read().
    sim_monthly_save(
      monthly_summary,
      file.path(output_part, paste0("monthly_summary", ext)),
      out_fmt
    )
    assignment_obj <- list(
      tu = tu,
      replicate_ids = replicate_ids,
      posterior_idx_occ = posterior_idx_occ,
      seed = seed,
      N_locations = N,
      M_steps = M,
      generated_at = Sys.time()
    )
    if (has_dur) {
      assignment_obj$posterior_idx_dur <- posterior_idx_dur
    }
    if (has_mon) {
      assignment_obj$posterior_idx_mon <- posterior_idx_mon
    }
    qs_save(assignment_obj, file.path(output_part, "assignment.qs"))
    writeLines(format(Sys.time()), file.path(output_part, "_SUCCESS"))
  }

  invisible(list(event_table = event_table, monthly_summary = monthly_summary))
}
