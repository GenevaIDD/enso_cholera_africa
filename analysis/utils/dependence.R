# dependence.R — size <-> duration posterior-dependence helpers for
#   05_fit_postprocess_dependence.R.
#
# Provides:
#   pit_size_nb()                randomized (Dunn-Smyth) PIT of the size NB
#   pit_duration_episode()       censoring-aware chained-hazard PIT, per episode
#   paired_kendall_tau()         Kendall's tau of paired PITs, per draw
#   re_intercept_correlation()   ADM1-level intercept correlation
#
# Details:
#   One PIT per episode (duration: from its LAST active-month row). Tau is
#   computed PER POSTERIOR DRAW (not on posterior means), so the result is a
#   distribution over tau. The ADM1 intercept correlation is the decomposition
#   step: do episode PITs correlate because ADM1s that run big also run long?
#   Design decisions are documented inline where they are made.

# ---------------------------------------------------------------------------
# 1. Size-side PIT (negative-binomial response, one row per episode)
# ---------------------------------------------------------------------------

#' pit_size_nb
#'
#' Randomized quantile (Dunn-Smyth) PIT for the size model's NB response.
#' u_i ~ Uniform(F(y_i - 1), F(y_i)) where F is the fitted NB CDF at draw d's
#' (mu, phi). Vectorized over (draw, row); one independent uniform draw per
#' (draw, row) cell, not one per row shared across draws — sharing would
#' correlate the PIT noise across draws and bias the propagated tau toward 0.
#'
#' @param fit brms fit (size / brms / negbinomial backend).
#' @param dat the exact data frame `fit` was trained on (e.g. `read_size_data(opt)`
#'   re-run against the same config) — supplies `episode_id`, which is not a
#'   covariate and so is never retained in `fit$data`. Row order must match
#'   `fit$data` 1:1 (asserted below); size fitting drops no rows via NA in
#'   current practice, so this normally just means "the same object".
#' @param response_var column in `dat` holding the observed outbreak size.
#' @return list(u = matrix [n_draws x n_episodes], episode_id = chr vector).
pit_size_nb <- function(fit, dat, response_var = "suspected_cases") {
  n_obs_dat <- nrow(dat)
  n_obs_fit <- nrow(fit$data)
  if (n_obs_dat != n_obs_fit) {
    stop(
      sprintf(
        "pit_size_nb: nrow(dat) = %d != nrow(fit$data) = %d -- `dat` is not the ",
        n_obs_dat,
        n_obs_fit
      ),
      "exact training data (brms silently dropped NA rows?). Re-derive `dat` ",
      "from the same config used to produce `fit` and re-check for NAs.",
      call. = FALSE
    )
  }
  if (!"episode_id" %in% names(dat)) {
    stop(
      "pit_size_nb: `dat` has no `episode_id` column -- read_size_data() ",
      "must be run on a repo state that includes the Phase-5 episode_id add.",
      call. = FALSE
    )
  }

  mu <- brms::posterior_epred(fit) # n_draws x n_obs, response scale (offsets included)
  phi <- brms::posterior_epred(fit, dpar = "shape") # n_draws x n_obs (replicated across cols if phi_re == "constant")
  y <- dat[[response_var]]
  if (is.null(y)) {
    stop("pit_size_nb: response_var '", response_var, "' not found in `dat`.")
  }

  list(u = .pit_nb_from_draws(mu, phi, y), episode_id = dat$episode_id)
}

#' .pit_nb_from_draws
#'
#' Pure-math core of pit_size_nb(): randomized NB quantile PIT from an
#' already-extracted (mu, phi, y) triple. Split out from pit_size_nb() so it
#' is testable on synthetic draws without constructing a real brms fit (see
#' analysis/tests/test_duration_dependence.R).
#'
#' @param mu matrix [n_draws x n_obs], NB mean per draw/obs (response scale).
#' @param phi matrix [n_draws x n_obs] (or a length-n_draws vector, recycled
#'   across columns), NB dispersion ("size" in stats::pnbinom's parameterisation).
#' @param y numeric vector [n_obs], observed counts.
#' @return matrix [n_draws x n_obs] of PIT values in (0, 1).
.pit_nb_from_draws <- function(mu, phi, y) {
  n_draws <- nrow(mu)
  n_obs <- ncol(mu)
  if (is.null(dim(phi))) {
    phi <- matrix(phi, nrow = n_draws, ncol = n_obs)
  }
  y_mat <- matrix(y, nrow = n_draws, ncol = n_obs, byrow = TRUE)

  lower <- stats::pnbinom(y_mat - 1, size = phi, mu = mu)
  upper <- stats::pnbinom(y_mat, size = phi, mu = mu)

  lower +
    matrix(stats::runif(n_draws * n_obs), n_draws, n_obs) * (upper - lower)
}

# ---------------------------------------------------------------------------
# 2. Duration-side PIT (chained discrete-time hazard, one row per episode)
# ---------------------------------------------------------------------------

#' pit_duration_episode
#'
#' Randomized, right-censoring-aware PIT for the duration cessation-hazard
#' model, one value per episode. With the Poisson mean lambda_i per row, the
#' discrete hazard gives log(1 - q_i) = -lambda_i, so the survival after month
#' j is S(j) = exp(-sum_{i<=j} lambda_i). Both S(k) and S(k-1) at the episode's
#' last row k come from one global cumulative sum over the (episode_id,
#' mso)-sorted panel. With F = 1 - S:
#'   event at month k:          u ~ Uniform(F(k-1), F(k))
#'   censored after month k:    u ~ Uniform(F(k), 1)
#'
#' @param fit brms fit (duration / brms backend; Poisson, response `y`).
#' @param dat the exact panel `fit` was trained on (e.g. `load_duration_data(opt)`
#'   re-run against the same config). Must carry `episode_id`, `mso`, `ended`.
#'   If `dat` has a `climate_missing` column, rows flagged TRUE are dropped
#'   first (mirrors the silent NA-drop brms performs at fit time -- see
#'   04_fit_dur_brms.R).
#' @return list(u = matrix [n_draws x n_episodes], episode_id = chr vector,
#'   ended = integer vector [1 = uncensored event, 0 = censored], aligned to
#'   the episode_id vector).
pit_duration_episode <- function(fit, dat) {
  if ("climate_missing" %in% names(dat)) {
    dat <- dat[which(!dat$climate_missing), , drop = FALSE] # keep iff FALSE (NA treated as "not kept")
  }
  n_obs_dat <- nrow(dat)
  n_obs_fit <- nrow(fit$data)
  if (n_obs_dat != n_obs_fit) {
    stop(
      sprintf(
        "pit_duration_episode: nrow(dat) = %d != nrow(fit$data) = %d after ",
        n_obs_dat,
        n_obs_fit
      ),
      "dropping climate_missing rows -- `dat` is not the exact training data. ",
      "Re-derive `dat` from the same config used to produce `fit`.",
      call. = FALSE
    )
  }
  need <- c("episode_id", "mso", "ended")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "pit_duration_episode: `dat` missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  lambda <- brms::posterior_epred(fit) # n_draws x n_rows, response scale (offset included)

  .pit_duration_from_draws(lambda, dat$episode_id, dat$mso, dat$ended)
}

#' .pit_duration_from_draws
#'
#' Pure-math core of pit_duration_episode(): the chained-hazard censored PIT,
#' from an already-extracted lambda matrix. Split out so it is testable on a
#' synthetic discrete-time hazard process without constructing a real brms fit
#' (see analysis/tests/test_duration_dependence.R) -- the generative process
#' for the test is exactly: draw lambda -> q = 1 - exp(-lambda) -> simulate
#' Bernoulli(q) per active month -> feed the same lambda back in here and check
#' the resulting PIT is Uniform(0,1).
#'
#' @param lambda matrix [n_draws x n_rows], columns in the SAME row order as
#'   episode_id/mso/ended (no internal re-alignment to `dat`'s original order
#'   is needed by the caller beyond what's passed in here).
#' @param episode_id, mso, ended vectors of length n_rows (panel rows, one row
#'   per active month; NOT required to be pre-sorted -- sorted internally).
#' @return list(u = matrix [n_draws x n_episodes], episode_id = chr vector,
#'   ended = integer vector), aligned to the (sorted, deduplicated) episode order.
.pit_duration_from_draws <- function(lambda, episode_id, mso, ended) {
  ord <- order(episode_id, mso)
  lambda_ord <- lambda[, ord, drop = FALSE]
  episode_id_ord <- episode_id[ord]
  ended_ord <- ended[ord]

  n_draws <- nrow(lambda_ord)

  # Global cumsum of -lambda over the whole sorted panel, per draw, with a
  # leading zero column so G0[, j+1] == cumsum through row j (G0[, 1] == 0 ==
  # "cumsum through row 0"). This lets any group's within-episode cumsum be
  # read off as a difference of two G0 lookups -- no per-episode loop, no
  # per-row re-basing.
  G <- t(apply(-lambda_ord, 1, cumsum)) # n_draws x n_rows
  G0 <- cbind(0, G) # n_draws x (n_rows + 1)

  rl <- rle(episode_id_ord)
  ep_ids <- rl$values
  grp_sizes <- rl$lengths
  last_idx <- cumsum(grp_sizes) # position of each episode's last row
  first_idx <- last_idx - grp_sizes + 1L # position of each episode's first row

  W_end <- G0[, last_idx + 1L, drop = FALSE] - G0[, first_idx, drop = FALSE] # within-episode cumsum through last row
  lambda_last <- lambda_ord[, last_idx, drop = FALSE]
  W_start <- W_end + lambda_last # within-episode cumsum through second-to-last row

  S_end <- exp(W_end)
  S_start <- exp(W_start)
  F_start <- 1 - S_start # F(k-1)
  F_end <- 1 - S_end # F(k)

  ended_last <- ended_ord[last_idx] # 1 = uncensored event at last row, 0 = censored
  n_ep <- length(ep_ids)
  ended_mat <- matrix(ended_last, nrow = n_draws, ncol = n_ep, byrow = TRUE)

  lo <- ifelse(ended_mat == 1, F_start, F_end)
  hi <- ifelse(ended_mat == 1, F_end, 1)

  u <- lo + matrix(stats::runif(n_draws * n_ep), n_draws, n_ep) * (hi - lo)

  list(u = u, episode_id = ep_ids, ended = as.integer(ended_last))
}

# ---------------------------------------------------------------------------
# 3. Paired Kendall's tau, propagated through draws
# ---------------------------------------------------------------------------

#' paired_kendall_tau
#'
#' Kendall's tau between two per-episode PIT sequences, computed
#' independently PER POSTERIOR DRAW and returned as a distribution (not a
#' point estimate). No monotonic pre-transform (e.g. qnorm()) is applied --
#' Kendall's tau is rank-based and invariant to any strictly monotonic
#' transform of either margin, so raw PIT/uniform values are used directly.
#'
#' If `pit_a` and `pit_b` have different draw counts (independently-run
#' fits, e.g. different chains/iter), draws are paired up to
#' min(n_draws_a, n_draws_b) in storage order. Post-warmup posterior draws
#' are approximately exchangeable within a fit, so this pairing is valid for
#' propagating uncertainty even though the two fits' draws are not jointly
#' sampled; it is NOT valid as a claim that draw d of fit A and draw d of fit
#' B represent "the same" posterior state.
#'
#' @param pit_a list(u = matrix, episode_id = chr) — e.g. from pit_size_nb().
#' @param pit_b list(u = matrix, episode_id = chr) — e.g. from pit_duration_episode().
#' @param group optional named chr vector (names = episode_id, values = group
#'   label, e.g. ENSO teleconnection group) to additionally break results out
#'   by group. When supplied, the returned list gains a `by_group` element.
#' @param min_n minimum number of paired episodes required to report a group
#'   (guards against a near-degenerate tau on tiny groups); groups below this
#'   are dropped from `by_group` with a message.
#' @return list(
#'   overall = list(tau_draws, mean, sd, q025, q975, n_episodes),
#'   by_group = named list of the same structure per group (or NULL)
#' )
paired_kendall_tau <- function(pit_a, pit_b, group = NULL, min_n = 5L) {
  common <- intersect(pit_a$episode_id, pit_b$episode_id)
  if (length(common) < min_n) {
    stop(
      "paired_kendall_tau: only ",
      length(common),
      " episode(s) in common between the two PIT sets (need >= ",
      min_n,
      ").",
      call. = FALSE
    )
  }
  idx_a <- match(common, pit_a$episode_id)
  idx_b <- match(common, pit_b$episode_id)

  n_draws <- min(nrow(pit_a$u), nrow(pit_b$u))

  .tau_summary <- function(cols_a, cols_b) {
    ua <- pit_a$u[seq_len(n_draws), cols_a, drop = FALSE]
    ub <- pit_b$u[seq_len(n_draws), cols_b, drop = FALSE]
    tau_draws <- vapply(
      seq_len(n_draws),
      function(d) {
        stats::cor(ua[d, ], ub[d, ], method = "kendall")
      },
      numeric(1)
    )
    list(
      tau_draws = tau_draws,
      mean = mean(tau_draws, na.rm = TRUE),
      sd = stats::sd(tau_draws, na.rm = TRUE),
      q025 = unname(stats::quantile(tau_draws, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(tau_draws, 0.975, na.rm = TRUE)),
      n_episodes = length(cols_a)
    )
  }

  overall <- .tau_summary(idx_a, idx_b)

  by_group <- NULL
  if (!is.null(group)) {
    grp_labels <- unname(group[common])
    if (any(is.na(grp_labels))) {
      message(sprintf(
        "paired_kendall_tau: %d/%d common episode(s) have no group label; excluded from by_group.",
        sum(is.na(grp_labels)),
        length(common)
      ))
    }
    by_group <- list()
    for (g in sort(unique(stats::na.omit(grp_labels)))) {
      sel <- which(grp_labels == g)
      if (length(sel) < min_n) {
        message(sprintf(
          "paired_kendall_tau: group '%s' has only %d paired episode(s) (< min_n = %d); skipped.",
          g,
          length(sel),
          min_n
        ))
        next
      }
      by_group[[g]] <- .tau_summary(idx_a[sel], idx_b[sel])
    }
  }

  list(overall = overall, by_group = by_group, n_common = length(common))
}

# ---------------------------------------------------------------------------
# 4. ADM1-level intercept correlation (decomposition step)
# ---------------------------------------------------------------------------

#' re_intercept_correlation
#'
#' Generic, schema-agnostic correlation between two ADM1-level point
#' estimates (typically the intercept RE, response/multiplicative scale)
#' from a size fit's and a duration fit's re_stats artefacts. Deliberately
#' thin: the two re_stats adapters use different column-naming conventions
#' (size's `adm1` table keys on `var`; duration's `re_stats_adm1` keys on
#' `term`), so the caller is expected to have already filtered each table
#' down to the intercept row per (country, adm1_unit) and to pass in plain
#' two-column-plus-key frames here.
#'
#' @param adm1_a data.frame with columns country, adm1_unit, value (e.g. the
#'   size fit's intercept `mean`, response scale).
#' @param adm1_b data.frame with columns country, adm1_unit, value (e.g. the
#'   duration fit's intercept `mean`, response scale).
#' @param method passed to stats::cor(); "pearson" default since these are
#'   point estimates (posterior means), not draws -- Kendall's tau is used
#'   above for the draw-level PIT correlation, where rank-invariance mattered.
#' @return list(r, n, data) where `data` is the inner-joined comparison frame
#'   (for scatter-plotting downstream), columns country, adm1_unit,
#'   value_a, value_b.
re_intercept_correlation <- function(adm1_a, adm1_b, method = "pearson") {
  need <- c("country", "adm1_unit", "value")
  miss_a <- setdiff(need, names(adm1_a))
  miss_b <- setdiff(need, names(adm1_b))
  if (length(miss_a)) {
    stop(
      "re_intercept_correlation: adm1_a missing: ",
      paste(miss_a, collapse = ", ")
    )
  }
  if (length(miss_b)) {
    stop(
      "re_intercept_correlation: adm1_b missing: ",
      paste(miss_b, collapse = ", ")
    )
  }

  joined <- merge(
    adm1_a[, need],
    adm1_b[, need],
    by = c("country", "adm1_unit"),
    suffixes = c("_a", "_b")
  )

  if (nrow(joined) < 3L) {
    warning(
      "re_intercept_correlation: only ",
      nrow(joined),
      " ADM1 unit(s) in common -- correlation is unreliable."
    )
  }

  r <- if (nrow(joined) >= 2L) {
    stats::cor(joined$value_a, joined$value_b, method = method)
  } else {
    NA_real_
  }

  list(r = r, n = nrow(joined), data = joined)
}
