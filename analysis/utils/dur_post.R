# dur_post.R — duration-track postprocessing: the posterior predictive check
#   of the cessation hazard and duration-ratio draws.
#
# Provides:
#   make_duration_pp_check_tbl()   hazard calibration + implied durations
#   duration_ratio_draws()         per-draw duration ratios
#   duration_ratio_summary()       their summary
#
# Details:
#   RE stats reuse the occurrence adapter (compute_nested_surv_re_stats_brms;
#   same Poisson cessation hazard), so only the PPC lives here.
#   The caller passes the duration panel filtered to the same rows the fit saw
#   (under grouping mode: resolve the group's countries into opt$country first).

# Posterior predictive check on the scale the model is used at (the response
# is a 0/1 cessation indicator per active month, so a size-style scatter means
# nothing):
#   1. hazard calibration: observed vs replicated monthly cessation hazard by
#      months-since-onset (checks the mso_band baseline);
#   2. implied durations: per draw, each episode ends at the first month its
#      replicated indicator fires, compared with observed durations.
# Panel 2 uses the observed at-risk window, so replicated cessation totals fall
# short of observed even when calibrated (meta$rep_nofire_median); read the
# shape, not the totals.
# Row alignment: fit$data lacks episode_id / mso, so the panel is re-screened by
# complete cases over the model variables and matched by position; row count and
# response are asserted. Under grouping mode the caller must resolve the group's
# countries into opt$country first (config.R::resolve_group_countries()).
#' Posterior-predictive check table for the duration cessation hazard.
#'
#' @param fit      brmsfit from 04_fit_dur_brms.R / 04_fit_dur_brms_cb.R.
#' @param dur_dat  duration panel filtered to the fit's rows (see note above);
#'                 needs columns episode_id, mso, y, mso_band.
#' @param ndraws   posterior draws to replicate (default 200 -- the summaries
#'                 below are quantiles, so more draws buy little).
#' @param seed     RNG seed for draw subsetting + Poisson replication.
#'
#' @return list(hazard, duration, meta). Both tibbles are pre-summarised to
#'   median + 95% intervals so the saved artefact stays small, mirroring the
#'   compactness of the size pp_check table.
make_duration_pp_check_tbl <- function(
  fit,
  dur_dat,
  ndraws = 200L,
  seed = 1234L
) {
  req <- c("episode_id", "mso", "y", "mso_band")
  miss <- setdiff(req, names(dur_dat))
  if (length(miss)) {
    stop(
      "make_duration_pp_check_tbl: duration panel missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  # --- align the panel to the rows brms actually fitted -------------------
  # brms drops rows with NA in any model variable (e.g. a missing ADM1
  # label), so reproduce that filter: complete cases over the model variables
  # the panel carries (enso_cb* are built in the fitter, absent here, and
  # NA-free). Order is kept, so rows align by position; the assertions below
  # verify it.
  model_vars <- intersect(all.vars(fit$formula$formula), names(dur_dat))
  keep <- stats::complete.cases(dur_dat[, model_vars, drop = FALSE])
  if (any(!keep)) {
    message(sprintf(
      paste0(
        "make_duration_pp_check_tbl: dropping %d of %d panel row(s) ",
        "with NA in a model variable (%s) -- brms dropped the same."
      ),
      sum(!keep),
      nrow(dur_dat),
      paste(
        model_vars[colSums(is.na(dur_dat[, model_vars, drop = FALSE])) > 0],
        collapse = ", "
      )
    ))
    dur_dat <- dur_dat[keep, , drop = FALSE]
  }

  if (nrow(dur_dat) != nrow(fit$data)) {
    stop(
      sprintf(
        paste0(
          "make_duration_pp_check_tbl: row mismatch -- panel has %d rows ",
          "after NA screening, fit has %d. Pass the panel filtered to ",
          "the fit's countries (resolve_group_countries() under ",
          "grouping mode)."
        ),
        nrow(dur_dat),
        nrow(fit$data)
      ),
      call. = FALSE
    )
  }
  if (!identical(as.integer(dur_dat$y), as.integer(fit$data$y))) {
    stop(
      "make_duration_pp_check_tbl: response vectors differ between panel and fit ",
      "-- positional alignment is not valid, refusing to guess.",
      call. = FALSE
    )
  }

  set.seed(seed)
  ndraws <- min(as.integer(ndraws), brms::ndraws(fit))

  # Poisson counts; the cessation EVENT is the indicator that the count fired.
  # (The Poisson-with-log-offset parameterisation is an approximation to a
  # discrete-time Bernoulli hazard, so counts >1 are possible and all mean
  # "ceased this month".)
  yrep <- brms::posterior_predict(fit, ndraws = ndraws)
  ev <- yrep >= 1

  epi <- as.character(dur_dat$episode_id)
  mso <- as.integer(dur_dat$mso)
  obs <- as.integer(dur_dat$y) >= 1L

  # --- 1. cessation hazard by months-since-onset -------------------------

  mso_lev <- sort(unique(mso))
  mso_f <- factor(mso, levels = mso_lev)

  at_risk <- as.integer(table(mso_f))
  obs_ev <- as.integer(tapply(obs, mso_f, sum))

  # (draws x mso-levels) replicated event counts
  rep_ev <- t(vapply(
    seq_len(nrow(ev)),
    function(i) as.numeric(tapply(ev[i, ], mso_f, sum)),
    numeric(length(mso_lev))
  ))
  rep_haz <- sweep(rep_ev, 2, at_risk, "/")
  haz_q <- apply(
    rep_haz,
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )

  # One band label per mso level (the band is a function of mso).
  band_by_mso <- as.character(dur_dat$mso_band)[match(mso_lev, mso)]

  hazard <- tibble::tibble(
    mso = mso_lev,
    mso_band = band_by_mso,
    at_risk = at_risk,
    obs_events = obs_ev,
    obs_hazard = obs_ev / at_risk,
    rep_lo = haz_q[1, ],
    rep_med = haz_q[2, ],
    rep_hi = haz_q[3, ]
  )

  # --- 2. implied duration distribution ----------------------------------
  #
  # Observed duration = the month the (single) observed cessation fired.
  # Right-censored episodes never fire and so contribute no observed duration;
  # they are counted in meta rather than silently folded into the histogram.

  grid <- seq_len(max(mso_lev))
  obs_dur <- tapply(mso[obs], epi[obs], min)
  obs_cnt <- as.integer(table(factor(as.integer(obs_dur), levels = grid)))

  rep_cnt <- matrix(0L, nrow = nrow(ev), ncol = length(grid))
  n_nofire <- integer(nrow(ev))
  n_epi <- length(unique(epi))
  for (i in seq_len(nrow(ev))) {
    hit <- ev[i, ]
    if (!any(hit)) {
      n_nofire[i] <- n_epi
      next
    }
    dur_i <- tapply(mso[hit], epi[hit], min)
    rep_cnt[i, ] <- as.integer(table(factor(as.integer(dur_i), levels = grid)))
    n_nofire[i] <- n_epi - length(dur_i)
  }
  dur_q <- apply(
    rep_cnt,
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )

  duration <- tibble::tibble(
    duration_months = grid,
    obs_n = obs_cnt,
    rep_lo = dur_q[1, ],
    rep_med = dur_q[2, ],
    rep_hi = dur_q[3, ]
  )

  list(
    hazard = hazard,
    duration = duration,
    meta = list(
      ndraws = nrow(ev),
      n_rows = nrow(dur_dat),
      n_episodes = n_epi,
      n_obs_cessations = sum(obs),
      n_censored = n_epi - sum(obs),
      # episodes with no replicated cessation inside follow-up, per draw
      rep_nofire_median = stats::median(n_nofire)
    )
  )
}


# ===========================================================================
# Duration RATIO (restricted-mean-duration) -- the estimand behind
# "outbreaks last XX times longer"
# ===========================================================================
# 1/HR from the cessation hazard equals a ratio of mean durations only under
# a constant baseline hazard; with the mso_band baseline it depends on S0(t).
# So "outbreaks last XX times longer" is the ratio of RESTRICTED mean durations
# (RMD) over a horizon K, per posterior draw.
# For an ENSO contrast this is a scenario ("if ONI stayed at +X.X
# throughout"). The caller builds newdata_hi / newdata_lo as one row per mso in
# 1..horizon, identical except the contrasted covariate, with the basis columns
# the fit's formula names.

#' Posterior duration ratio for a covariate contrast, duration track.
#'
#' Returns draws of RMD(x_hi) / RMD(x_lo), where RMD is the restricted mean
#' duration over `horizon` months. THIS is the quantity that licenses
#' "outbreaks last XX times longer"; 1/HR is not.
#'
#' @param fit        brmsfit from 04_fit_dur_brms.R / 04_fit_dur_brms_cb.R.
#' @param newdata_hi,newdata_lo  one row per mso in 1..horizon, identical in
#'   every column EXCEPT the contrasted covariate; must already carry whatever
#'   basis columns the formula names (mso_band, msoms*, seasms*, enso_cb*),
#'   built with the same basis as the fit -- never re-fit.
#' @param horizon    K months; truncation is explicit, not incidental.
#' @param days_ref   exposure offset per month, held FIXED (default 30.4375).
#'   Letting it vary would mix exposure differences into the contrast.
#' @param re_formula NA = population-level (a typical-group statement);
#'   NULL = includes REs (a specific-country statement). Different claims.
#' @param tail_tol   if median S(K) > tail_tol the UNrestricted mean is not
#'   identified from the observed window -- warn and keep it restricted.
#'
#' @return list(ratio, rmd_hi, rmd_lo, horizon, tail_mass). `ratio` is a vector
#'   of per-draw RMD ratios; summarise it with duration_ratio_summary().
duration_ratio_draws <- function(
  fit,
  newdata_hi,
  newdata_lo,
  horizon = 24L,
  days_ref = 30.4375,
  re_formula = NA,
  tail_tol = 0.05
) {
  stopifnot(nrow(newdata_hi) == horizon, nrow(newdata_lo) == horizon)

  # --- per-draw monthly cessation RATE on the fitted scale -----------------
  # Force the offset to the reference month length instead of whatever the
  # panel row carried, so exposure differences never leak into the contrast.
  rate <- function(nd) {
    nd$log_offset <- log(days_ref)
    exp(brms::posterior_linpred(fit, newdata = nd, re_formula = re_formula))
  } # draws x horizon

  mu_hi <- rate(newdata_hi)
  mu_lo <- rate(newdata_lo)

  # --- rate -> survival ----------------------------------------------------
  # Poisson-with-offset is the cloglog discrete hazard (Holford 1980):
  #   h(k) = 1 - exp(-mu(k)),  S(t) = prod_{k<=t}(1-h(k)) = exp(-cumsum(mu))
  # so the survival curve needs no product at all.
  surv <- function(mu) exp(-t(apply(mu, 1, cumsum))) # draws x horizon, S(1..K)

  S_hi <- surv(mu_hi)
  S_lo <- surv(mu_lo)

  # --- restricted mean duration: E[T ^ K] = sum_{t=0}^{K-1} S(t), S(0)=1 ---
  rmd <- function(S) 1 + rowSums(S[, -ncol(S), drop = FALSE])

  rmd_hi <- rmd(S_hi)
  rmd_lo <- rmd(S_lo)

  tail_mass <- stats::median(pmax(S_hi[, horizon], S_lo[, horizon]))
  if (tail_mass > tail_tol) {
    warning(sprintf(
      paste(
        "duration_ratio_draws: median S(%d) = %.3f > %.3f --",
        "report this as a %d-month RESTRICTED mean duration ratio;",
        "the unrestricted mean is not identified over this window."
      ),
      horizon,
      tail_mass,
      tail_tol,
      horizon
    ))
  }

  # Ratio PER DRAW, then summarise -- never a ratio of posterior summaries.
  list(
    ratio = rmd_hi / rmd_lo,
    rmd_hi = rmd_hi,
    rmd_lo = rmd_lo,
    horizon = horizon,
    tail_mass = tail_mass
  )
}

#' Median + CrI summary of a duration_ratio_draws() result, ready for a sentence.
#'
#' @param dr    list returned by duration_ratio_draws().
#' @param probs quantile probs (default 2.5% / 50% / 97.5%).
#' @return one-row tibble with ratio_lo/med/hi, median RMDs, horizon, tail_mass,
#'   and prob_longer = P(ratio > 1).
duration_ratio_summary <- function(dr, probs = c(0.025, 0.5, 0.975)) {
  q <- stats::quantile(dr$ratio, probs, names = FALSE)
  tibble::tibble(
    ratio_lo = q[1],
    ratio_med = q[2],
    ratio_hi = q[3],
    rmd_hi_med = stats::median(dr$rmd_hi),
    rmd_lo_med = stats::median(dr$rmd_lo),
    horizon = dr$horizon,
    tail_mass = dr$tail_mass,
    prob_longer = mean(dr$ratio > 1)
  )
}
