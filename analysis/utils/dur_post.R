# ===========================================================================
# dur_post.R -- outbreak-DURATION track postprocessing adapters
# ===========================================================================
#
# Sibling of occ_post.R / siz_post.R. The duration track reuses the OCCURRENCE
# RE-stats adapter unchanged (computeNestedSurvReStatsBrms -- same Poisson
# cessation hazard, same AFRO_region/country/ADM1 nesting), so nothing here
# duplicates that. What DOES need its own implementation is the posterior
# predictive check, because the duration response is not a quantity you can
# scatter against a prediction.
#
# WHY THE SIZE PP_CHECK DOES NOT TRANSFER. makeSizePpCheckTbl() (siz_post.R:3)
# is a predicted-vs-observed scatter of `suspected_cases` -- meaningful because
# the size response IS the endpoint of interest. The duration response is
# `y = as.integer(ended)`, a 0/1 cessation indicator on one row per ACTIVE
# MONTH (dur_data.R:9-13). Scattering predicted against observed 0/1 per
# month answers nothing: the endpoint the model is fitted FOR -- how long an
# outbreak lasts -- is a DERIVED quantity, obtained by pushing the per-month
# hazards through to a duration.
#
# So this PPC checks the two things that actually matter, both on the scale the
# model is used at:
#
#   1. HAZARD CALIBRATION -- observed monthly cessation hazard (events /
#      at-risk) against the posterior-predictive hazard, by months-since-onset.
#      This is the direct check on the `mso_band` baseline hazard, which spans
#      a ~492x range across bands (dur_data.R:50-57) and is therefore the part
#      of the model most able to be wrong in a visible way.
#
#   2. IMPLIED DURATION DISTRIBUTION -- the endpoint. For every posterior
#      draw, each episode's implied duration is the first month its replicated
#      cessation indicator fires; the resulting distribution is compared to the
#      observed duration distribution. A model can match the marginal hazard
#      per band and still get the duration distribution wrong (that is exactly
#      what unmodelled episode-level frailty does), so panel 2 is not implied
#      by panel 1.
#
# HOW TO READ PANEL 2 (do not over-read the totals). The panel is a FIXED
# at-risk window: an episode observed to last 4 months contributes exactly 4
# rows, and for an uncensored episode the last of those rows carries y = 1 BY
# CONSTRUCTION. The replicate gets no such conditioning -- it simply asks
# whether the indicator fires somewhere in those 4 months, which at a hazard
# of ~0.2 happens only 1 - 0.8^4 ~= 59% of the time. So replicated cessation
# TOTALS are expected to fall short of observed totals even for a
# well-calibrated hazard, and the shortfall is reported explicitly as
# meta$rep_nofire_median rather than hidden. The informative comparison is the
# SHAPE of the duration distribution and the per-month hazard calibration in
# panel 1 -- not the ratio of replicated to observed cessation counts.
# Simulating past the observed window would need covariates (climate, ENSO) at
# months the panel does not contain, so it is deliberately not attempted here.
#
# ROW ALIGNMENT (load-bearing). `fit$data` carries only formula variables --
# y, the crossbasis columns, mso_band, country, ADM1, log_offset -- and NOT
# episode_id / mso, so durations cannot be derived from the fit alone. The
# panel must be re-loaded and matched positionally, and because brms drops any
# row with an NA in a model variable, the panel is first screened by complete
# cases over those same variables. This is NOT hypothetical: 41 of the 2710
# enso_west rows carry a missing ADM1 label, which the (1 | country/ADM1)
# grouping term removes. brms preserves row ORDER, so the survivors align
# positionally. Both the row count AND the response vector are still asserted
# after that screen, so a formula change which breaks the correspondence
# becomes a hard error here instead of a silently mislabelled figure.
#
# The caller is responsible for passing a panel filtered to the SAME rows the
# fit saw. Under grouping mode that means resolving the group's countries into
# opt$country first (config.R:695 resolve_group_countries) -- the raw config
# leaves opt$country = "all", which would hand over the full 9062-row panel.
# ===========================================================================


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
makeDurationPpCheckTbl <- function(fit, dur_dat, ndraws = 200L, seed = 1234L) {

  req <- c("episode_id", "mso", "y", "mso_band")
  miss <- setdiff(req, names(dur_dat))
  if (length(miss))
    stop("makeDurationPpCheckTbl: duration panel missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

  # --- align the panel to the rows brms actually fitted -------------------
  #
  # brms drops any row with NA in ANY model variable, so the group panel is not
  # in general the same height as fit$data. The enso_west panel is the live
  # example: 41 of its 2710 rows carry a missing ADM1 label, so the
  # (1 | country/ADM1) grouping term removes them and the raw panel is 41 rows
  # too long. Excluding precip/temp lags does not protect against a missing RE
  # grouping label.
  #
  # Rather than assume alignment, reproduce brms's own filter: complete cases
  # over the model variables the panel carries. Row ORDER is untouched by
  # brms, so the surviving rows line up positionally. Intersecting with
  # names(dur_dat) is what makes this safe for the crossbasis backend, whose
  # enso_cb* columns are built inside the fitter and are absent from the panel
  # (they are NA-free functions of enso_lag<N>, which is itself a model
  # variable and so already screened).
  #
  # The two assertions below VERIFY this filter -- if a future formula change breaks the
  # correspondence, this still fails loudly rather than mislabelling a figure.
  model_vars <- intersect(all.vars(fit$formula$formula), names(dur_dat))
  keep       <- stats::complete.cases(dur_dat[, model_vars, drop = FALSE])
  if (any(!keep)) {
    message(sprintf(paste0("makeDurationPpCheckTbl: dropping %d of %d panel row(s) ",
                           "with NA in a model variable (%s) -- brms dropped the same."),
                    sum(!keep), nrow(dur_dat),
                    paste(model_vars[colSums(is.na(dur_dat[, model_vars, drop = FALSE])) > 0],
                          collapse = ", ")))
    dur_dat <- dur_dat[keep, , drop = FALSE]
  }

  if (nrow(dur_dat) != nrow(fit$data))
    stop(sprintf(paste0("makeDurationPpCheckTbl: row mismatch -- panel has %d rows ",
                        "after NA screening, fit has %d. Pass the panel filtered to ",
                        "the fit's countries (resolve_group_countries() under ",
                        "grouping mode)."),
                 nrow(dur_dat), nrow(fit$data)), call. = FALSE)
  if (!identical(as.integer(dur_dat$y), as.integer(fit$data$y)))
    stop("makeDurationPpCheckTbl: response vectors differ between panel and fit ",
         "-- positional alignment is not valid, refusing to guess.", call. = FALSE)

  set.seed(seed)
  ndraws <- min(as.integer(ndraws), brms::ndraws(fit))

  # Poisson counts; the cessation EVENT is the indicator that the count fired.
  # (The Poisson-with-log-offset parameterisation is an approximation to a
  # discrete-time Bernoulli hazard, so counts >1 are possible and all mean
  # "ceased this month".)
  yrep <- brms::posterior_predict(fit, ndraws = ndraws)
  ev   <- yrep >= 1

  epi <- as.character(dur_dat$episode_id)
  mso <- as.integer(dur_dat$mso)
  obs <- as.integer(dur_dat$y) >= 1L

  # --- 1. cessation hazard by months-since-onset -------------------------

  mso_lev <- sort(unique(mso))
  mso_f   <- factor(mso, levels = mso_lev)

  at_risk <- as.integer(table(mso_f))
  obs_ev  <- as.integer(tapply(obs, mso_f, sum))

  # (draws x mso-levels) replicated event counts
  rep_ev <- t(vapply(seq_len(nrow(ev)),
                     function(i) as.numeric(tapply(ev[i, ], mso_f, sum)),
                     numeric(length(mso_lev))))
  rep_haz <- sweep(rep_ev, 2, at_risk, "/")
  haz_q   <- apply(rep_haz, 2, stats::quantile,
                   probs = c(0.025, 0.5, 0.975), names = FALSE)

  # One band label per mso level (the band is a function of mso).
  band_by_mso <- as.character(dur_dat$mso_band)[match(mso_lev, mso)]

  hazard <- tibble::tibble(
    mso         = mso_lev,
    mso_band    = band_by_mso,
    at_risk     = at_risk,
    obs_events  = obs_ev,
    obs_hazard  = obs_ev / at_risk,
    rep_lo      = haz_q[1, ],
    rep_med     = haz_q[2, ],
    rep_hi      = haz_q[3, ]
  )

  # --- 2. implied duration distribution ----------------------------------
  #
  # Observed duration = the month the (single) observed cessation fired.
  # Right-censored episodes never fire and so contribute no observed duration;
  # they are counted in meta rather than silently folded into the histogram.

  grid    <- seq_len(max(mso_lev))
  obs_dur <- tapply(mso[obs], epi[obs], min)
  obs_cnt <- as.integer(table(factor(as.integer(obs_dur), levels = grid)))

  rep_cnt <- matrix(0L, nrow = nrow(ev), ncol = length(grid))
  n_nofire <- integer(nrow(ev))
  n_epi    <- length(unique(epi))
  for (i in seq_len(nrow(ev))) {
    hit <- ev[i, ]
    if (!any(hit)) { n_nofire[i] <- n_epi; next }
    dur_i <- tapply(mso[hit], epi[hit], min)
    rep_cnt[i, ] <- as.integer(table(factor(as.integer(dur_i), levels = grid)))
    n_nofire[i]  <- n_epi - length(dur_i)
  }
  dur_q <- apply(rep_cnt, 2, stats::quantile,
                 probs = c(0.025, 0.5, 0.975), names = FALSE)

  duration <- tibble::tibble(
    duration_months = grid,
    obs_n           = obs_cnt,
    rep_lo          = dur_q[1, ],
    rep_med         = dur_q[2, ],
    rep_hi          = dur_q[3, ]
  )

  list(
    hazard   = hazard,
    duration = duration,
    meta     = list(
      ndraws            = nrow(ev),
      n_rows            = nrow(dur_dat),
      n_episodes        = n_epi,
      n_obs_cessations  = sum(obs),
      n_censored        = n_epi - sum(obs),
      # episodes with no replicated cessation inside follow-up, per draw
      rep_nofire_median = stats::median(n_nofire)
    )
  )
}


# ===========================================================================
# Duration RATIO (restricted-mean-duration) -- the estimand behind
# "outbreaks last XX times longer"
# ===========================================================================
#
# WHY THIS EXISTS. The duration track is a discrete-time cessation hazard
# (Poisson-with-log-offset, Holford form), so exp(beta) is a hazard ratio for
# *ending*, and 1/HR is a MONOTONE reorientation of it -- nothing more. 1/HR
# equals a ratio of expected durations ONLY under a constant baseline hazard
# (basehaz_duration: none): constant hazard => T ~ Exponential(h) =>
# E[T|x] = 1/h(x) => E[T|x1]/E[T|x0] = 1/HR exactly. Under the default factor/ms
# baseline the survival curve is S(t|x) = S0(t)^HR, whose mean-duration ratio
# depends on the SHAPE of S0(t) -- and the mso_band hazards span a ~492x range,
# so that shape carries the statistic. The damage concentrates in months 1-3,
# where almost nothing ends and which dominate the mean. So the sentence
# "outbreaks last XX times longer" CANNOT be sourced from 1/HR.
#
# The honest estimand is the ratio of RESTRICTED mean durations (RMD) over an
# explicit horizon K, computed per posterior draw and then summarised. This is
# the closed-form cousin of the simulation-based implied-duration machinery in
# makeDurationPpCheckTbl() above.
#
# CONTRAST SEMANTICS. For an ENSO contrast this is a *scenario* statement: the
# covariate is held at one value for the whole episode (which no real episode
# does), so any reported number must be worded "if ONI stayed at +X.X
# throughout". The caller builds newdata_hi / newdata_lo as one row per mso in
# 1..horizon, identical in every column EXCEPT the contrasted covariate, using
# the SAME basis columns the fit's formula names (mso_band, enso_cb*, etc.).

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
#'   built with the SAME basis as the fit -- never re-fit.
#' @param horizon    K months; truncation is explicit, not incidental.
#' @param days_ref   exposure offset per month, held FIXED (default 30.4375).
#'   Letting it vary would mix exposure differences into the contrast.
#' @param re_formula NA = population-level (a typical-group statement);
#'   NULL = includes REs (a specific-country statement). Different claims.
#' @param tail_tol   if median S(K) > tail_tol the UNrestricted mean is not
#'   identified from the observed window -- warn and keep it restricted.
#'
#' @return list(ratio, rmd_hi, rmd_lo, horizon, tail_mass). `ratio` is a vector
#'   of per-draw RMD ratios; summarise it with durationRatioSummary().
durationRatioDraws <- function(fit, newdata_hi, newdata_lo,
                               horizon    = 24L,
                               days_ref   = 30.4375,
                               re_formula = NA,
                               tail_tol   = 0.05) {

  stopifnot(nrow(newdata_hi) == horizon, nrow(newdata_lo) == horizon)

  # --- per-draw monthly cessation RATE on the fitted scale -----------------
  # Force the offset to the reference month length instead of whatever the
  # panel row carried, so exposure differences never leak into the contrast.
  rate <- function(nd) {
    nd$log_offset <- log(days_ref)
    exp(brms::posterior_linpred(fit, newdata = nd, re_formula = re_formula))
  }                                      # draws x horizon

  mu_hi <- rate(newdata_hi)
  mu_lo <- rate(newdata_lo)

  # --- rate -> survival ----------------------------------------------------
  # Poisson-with-offset is the cloglog discrete hazard (Holford 1980):
  #   h(k) = 1 - exp(-mu(k)),  S(t) = prod_{k<=t}(1-h(k)) = exp(-cumsum(mu))
  # so the survival curve needs no product at all.
  surv <- function(mu) exp(-t(apply(mu, 1, cumsum)))     # draws x horizon, S(1..K)

  S_hi <- surv(mu_hi)
  S_lo <- surv(mu_lo)

  # --- restricted mean duration: E[T ^ K] = sum_{t=0}^{K-1} S(t), S(0)=1 ---
  rmd <- function(S) 1 + rowSums(S[, -ncol(S), drop = FALSE])

  rmd_hi <- rmd(S_hi)
  rmd_lo <- rmd(S_lo)

  tail_mass <- stats::median(pmax(S_hi[, horizon], S_lo[, horizon]))
  if (tail_mass > tail_tol)
    warning(sprintf(paste("durationRatioDraws: median S(%d) = %.3f > %.3f --",
                          "report this as a %d-month RESTRICTED mean duration ratio;",
                          "the unrestricted mean is not identified over this window."),
                    horizon, tail_mass, tail_tol, horizon))

  # Ratio PER DRAW, then summarise -- never a ratio of posterior summaries.
  list(ratio     = rmd_hi / rmd_lo,
       rmd_hi    = rmd_hi,
       rmd_lo    = rmd_lo,
       horizon   = horizon,
       tail_mass = tail_mass)
}

#' Median + CrI summary of a durationRatioDraws() result, ready for a sentence.
#'
#' @param dr    list returned by durationRatioDraws().
#' @param probs quantile probs (default 2.5% / 50% / 97.5%).
#' @return one-row tibble with ratio_lo/med/hi, median RMDs, horizon, tail_mass,
#'   and prob_longer = P(ratio > 1).
durationRatioSummary <- function(dr, probs = c(0.025, 0.5, 0.975)) {
  q <- stats::quantile(dr$ratio, probs, names = FALSE)
  tibble::tibble(ratio_lo = q[1], ratio_med = q[2], ratio_hi = q[3],
                 rmd_hi_med = stats::median(dr$rmd_hi),
                 rmd_lo_med = stats::median(dr$rmd_lo),
                 horizon = dr$horizon, tail_mass = dr$tail_mass,
                 prob_longer = mean(dr$ratio > 1))
}
