# helpers_synth.R — synthetic occurrence panels with known truth, shared by the
#   MCMC tests (test_synthetic_equivalence.R, test_kl_*.R).
#
# Provides:
#   make_synthetic_occ_panel()   region / country / ADM1 intercept panel
#   make_synthetic_kl_panel()    + ENSO series, lag-kernel effect, adjacency
#   adjusted_rand_index()        partition agreement
#   kl_synthetic_opt()           brms_kl options for the synthetic panel
#   make_fake_kl_fit()           Stan-free fake kl_fit
#
# Details:
#   Sourced by tests only. Requires dplyr/tidyr/tibble (tidyverse) and, for the
#   kernel panel, analysis/utils/base.R + enso.R (join_enso_lags,
#   enso_lag_col_name) and lubridate.

# make_synthetic_occ_panel() keeps test_synthetic_equivalence.R's RNG call
# order exactly, so that test's draws are unchanged. make_synthetic_kl_panel()
# extends it with an ENSO-like series joined through the real
# enso.R::join_enso_lags(), a Gaussian lag-kernel ENSO effect per latent
# country class, and a line-shaped country adjacency.
#' make_synthetic_occ_panel
#'
#' Region -> country -> ADM1 random intercepts, a monthly panel with
#' sinusoidal climatology + noisy anomalies, a seasonal baseline log-hazard,
#' 30-day exposures and Poisson counts. With the defaults (countries = NULL,
#' adm1_per_country = NULL, augment = NULL) the RNG calls and their order are
#' exactly those of the former inline builder of test_synthetic_equivalence.R.
#'
#' @param countries optional tibble(country, AFRO_region) fixing the country
#'   layer (skips the random region assignment)
#' @param adm1_per_country optional integer: deterministic ADM1 allocation
#'   (that many units per country) instead of random sampling
#' @param augment optional function(panel) -> panel run after the covariates
#'   and before the Poisson draw; it must add a numeric `.eta_extra` column,
#'   which is added to the linear predictor
#' @return tibble panel (truth columns re_r / re_c / re_a / eta / lambda kept)
make_synthetic_occ_panel <- function(
  n_region = 2,
  n_country = 5,
  n_adm1 = 50,
  n_year = 10,
  beta_true = c(
    mean_precip = 0.30,
    dev_precip = -0.20,
    mean_temp_daymean = 0.15,
    dev_temp_daymean = 0.05
  ),
  sd_region = 0.40,
  sd_country = 0.20,
  sd_adm1 = 0.10,
  log_h0_intercept = -5,
  countries = NULL,
  adm1_per_country = NULL,
  augment = NULL
) {
  regions <- sprintf("R%d", seq_len(n_region))
  if (is.null(countries)) {
    countries <- tibble::tibble(
      country = sprintf("C%02d", seq_len(n_country)),
      AFRO_region = sample(regions, n_country, replace = TRUE)
    )
  }
  adm1_country <- if (is.null(adm1_per_country)) {
    sample(countries$country, n_adm1, replace = TRUE)
  } else {
    rep(countries$country, each = adm1_per_country)
  }
  adm1s <- tibble::tibble(
    ADM1 = sprintf("%s_%d", adm1_country, seq_along(adm1_country))
  ) %>%
    dplyr::mutate(country = stringr::str_extract(ADM1, "^[^_]+")) %>%
    dplyr::inner_join(countries, by = "country")

  re_region <- tibble::tibble(
    AFRO_region = regions,
    re_r = stats::rnorm(length(regions), 0, sd_region)
  )
  re_country <- countries %>%
    dplyr::mutate(re_c = stats::rnorm(dplyr::n(), 0, sd_country))
  re_adm1 <- adm1s %>%
    dplyr::mutate(re_a = stats::rnorm(dplyr::n(), 0, sd_adm1))

  panel <- tidyr::expand_grid(
    ADM1 = adm1s$ADM1,
    year = seq_len(n_year),
    month = 1:12
  ) %>%
    dplyr::inner_join(adm1s, by = "ADM1") %>%
    dplyr::inner_join(re_region, by = "AFRO_region") %>%
    dplyr::inner_join(re_country, by = c("country", "AFRO_region")) %>%
    dplyr::inner_join(re_adm1, by = c("ADM1", "country", "AFRO_region"))

  panel <- panel %>%
    dplyr::mutate(
      mean_precip = sin(2 * pi * month / 12),
      dev_precip = stats::rnorm(dplyr::n(), 0, 1),
      mean_temp_daymean = cos(2 * pi * month / 12),
      dev_temp_daymean = stats::rnorm(dplyr::n(), 0, 1),
      log_h0 = log_h0_intercept + 0.6 * sin(2 * pi * (month - 2) / 12)
    )
  if (!is.null(augment)) {
    panel <- augment(panel)
  }
  extra <- if (".eta_extra" %in% names(panel)) panel$.eta_extra else 0

  # Same addition order as the former inline builder (x + 0 == x exactly, so
  # the default path's eta is bit-identical).
  panel %>%
    dplyr::mutate(
      eta = log_h0 +
        beta_true["mean_precip"] * mean_precip +
        beta_true["dev_precip"] * dev_precip +
        beta_true["mean_temp_daymean"] * mean_temp_daymean +
        beta_true["dev_temp_daymean"] * dev_temp_daymean +
        re_r +
        re_c +
        re_a +
        extra,
      lambda = exp(eta),
      exposure = 30L,
      y = stats::rpois(dplyr::n(), lambda * exposure),
      outbreak = y > 0,
      y_surv = as.integer(outbreak),
      t_start_seas = (month - 1) * 30,
      t_end_seas = month * 30,
      months = factor(month, levels = 1:12)
    )
}

#' make_synthetic_kl_panel
#'
#' Synthetic occurrence panel with a Gaussian lag-kernel ENSO effect that
#' differs by latent country class. Countries C01..Cn sit on a line (C01-C02,
#' C02-C03, ...), which gives the spatial tests an adjacency; by default the
#' classes are contiguous blocks along the line. ENSO is a monthly AR(1) series
#' shared by all units, joined with enso.R::join_enso_lags() so the lag columns
#' are exactly the pipeline's enso_lead<k> / enso_lag<k>.
#'
#' ENSO effect of a row in country c (class k = z[c]):
#'   beta[k] * sum_l w_l(mu[k], sigma[k]) * ONI(t - l),  l = lag_min..lag_max
#' (linear in ONI -- representable by a natural-spline onebasis).
#'
#' @param mu_true,sigma_true,beta_true per-class kernel truth (length K)
#' @param z_true integer class per country (default: contiguous blocks)
#' @param enso_scale multiplies a country's ENSO effect (vector per country;
#'   default 1) -- lets a test weaken one country's signal
#' @return panel with columns the brms_kl data path reads (y, log_offset,
#'   exposure, months, t_start, t_start_seas, t_end_seas, country (factor),
#'   AFRO_region, ADM1, gid, enso_lead*/enso_lag*), attributes
#'   truth = list(mu, sigma, beta, z = tibble(iso3, class), lags) and
#'   adjacency = data.frame(iso3_a, iso3_b)
make_synthetic_kl_panel <- function(
  n_country = 6L,
  adm1_per_country = 4L,
  n_year = 10L,
  K = 1L,
  mu_true = 3,
  sigma_true = 1.5,
  beta_true = 0.8,
  z_true = NULL,
  enso_scale = NULL,
  lag_min = -11L,
  lag_max = 11L,
  rho_enso = 0.7,
  log_h0_intercept = -4,
  sd_country = 0.2,
  sd_adm1 = 0.1,
  seed = 1L
) {
  set.seed(seed)
  stopifnot(
    length(mu_true) == K,
    length(sigma_true) == K,
    length(beta_true) == K
  )
  iso <- sprintf("C%02d", seq_len(n_country))
  if (is.null(z_true)) {
    z_true <- sort(rep_len(seq_len(K), n_country))
  }
  if (is.null(enso_scale)) {
    enso_scale <- rep(1, n_country)
  }
  countries <- tibble::tibble(
    country = iso,
    AFRO_region = sprintf("R%d", (seq_len(n_country) - 1L) %% 2L + 1L)
  )

  # Monthly ENSO-like AR(1) series covering every lag/lead of the panel window
  start <- as.Date("2010-01-01")
  months_all <- seq(
    start - 31 * (abs(lag_max) + 14),
    by = "month",
    length.out = 12L * n_year + abs(lag_max) + abs(lag_min) + 30L
  )
  months_all <- seq(
    as.Date(format(months_all[1], "%Y-%m-01")),
    by = "month",
    length.out = length(months_all)
  )
  x <- numeric(length(months_all))
  x[1] <- stats::rnorm(1)
  for (t in 2:length(x)) {
    x[t] <- rho_enso * x[t - 1] + stats::rnorm(1, 0, sqrt(1 - rho_enso^2))
  }
  enso_dat <- tibble::tibble(
    year = as.integer(format(months_all, "%Y")),
    month = as.integer(format(months_all, "%m")),
    enso_index = x
  )
  lags <- seq.int(lag_min, lag_max)
  W <- exp(-0.5 * outer(-mu_true, lags, `+`)^2 / sigma_true^2)
  W <- W / rowSums(W)

  augment <- function(panel) {
    panel <- panel %>%
      dplyr::mutate(
        t_start = as.Date(sprintf("%d-%02d-01", 2009L + year, month))
      )
    panel <- join_enso_lags(
      panel,
      enso_dat,
      date_col = "t_start",
      max_lag = max(lag_max, 0L),
      max_lead = max(-lag_min, 0L)
    )
    X <- as.matrix(panel[, enso_lag_col_name(lags)])
    ci <- match(panel$country, iso)
    k <- z_true[ci]
    panel$.eta_extra <- enso_scale[ci] *
      beta_true[k] *
      rowSums(X * W[k, , drop = FALSE])
    panel
  }

  panel <- make_synthetic_occ_panel(
    n_region = 2,
    n_country = n_country,
    n_year = n_year,
    beta_true = c(
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    ),
    sd_region = 0,
    sd_country = sd_country,
    sd_adm1 = sd_adm1,
    log_h0_intercept = log_h0_intercept,
    countries = countries,
    adm1_per_country = adm1_per_country,
    augment = augment
  )

  panel <- panel %>%
    dplyr::mutate(
      country = factor(country, levels = iso),
      AFRO_region = factor(AFRO_region),
      gid = ADM1,
      log_offset = log(exposure),
      outbreak = as.integer(outbreak)
    ) %>%
    dplyr::select(-.eta_extra)
  attr(panel, "truth") <- list(
    mu = mu_true,
    sigma = sigma_true,
    beta = beta_true,
    z = tibble::tibble(iso3 = iso, class = z_true),
    lags = lags
  )
  attr(panel, "adjacency") <- data.frame(
    iso3_a = iso[-n_country],
    iso3_b = iso[-1L],
    stringsAsFactors = FALSE
  )
  panel
}

#' adjusted_rand_index -- agreement of two partitions (1 = identical up to labels).
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(n) n * (n - 1) / 2
  s_ij <- sum(comb2(tab))
  s_a <- sum(comb2(rowSums(tab)))
  s_b <- sum(comb2(colSums(tab)))
  n <- comb2(sum(tab))
  expected <- s_a * s_b / n
  maxv <- (s_a + s_b) / 2
  if (maxv == expected) {
    return(1)
  }
  (s_ij - expected) / (maxv - expected)
}

#' kl_synthetic_opt
#'
#' Minimal resolved options for fitting a make_synthetic_kl_panel() panel with
#' the brms_kl backend's module functions: every covariate off, nested
#' region/country/ADM1 intercepts only, ENSO lag kernel over lag_min..lag_max
#' on a natural-spline basis, and the given kernel overrides (e.g.
#' list(classes = 2L)). Requires analysis/utils/config.R.
kl_synthetic_opt <- function(
  lag_min = -11L,
  lag_max = 11L,
  df = 3L,
  kernel = list()
) {
  opt <- get_default_options()
  opt$occurrence$backend <- "brms_kl"
  opt$occurrence$covariates[] <- lapply(opt$occurrence$covariates, function(x) {
    if (is.logical(x)) FALSE else x
  })
  re <- opt$occurrence$random_effects
  for (nm in names(re)) {
    if (is.logical(re[[nm]])) re[[nm]] <- FALSE
  }
  re$re_region_country_adm1_intercept <- TRUE
  opt$occurrence$random_effects <- re
  opt$occurrence$dlnm <- list(
    enable = TRUE,
    defaults = list(argvar = list(fun = "ns", df = df)),
    variables = list(enso = list(lags = list(lag_min, lag_max)))
  )
  for (k in names(kernel)) {
    opt$occurrence$kernel[[k]] <- kernel[[k]]
  }
  opt
}

#' make_fake_kl_fit
#'
#' A Stan-free object of class "kl_fit" with the fields the summary helpers
#' read (utils/dlnm_kernel.R, utils/kl_sweep.R, utils/diagnostics.R): draws
#' (posterior::draws_df with kmu[k], kcon[k], kpi[k], kkappa, lp__,
#' kl_log_resp[k,c], kl_log_lik_ctry[c]), diagnostics (per-chain
#' num_divergent / num_max_treedepth / ebfmi), K, class_labels,
#' country_levels, data_dims, chains.
#'
#' @param centres chains x K matrix: each chain's kmu centre per class (chains
#'   with similar rows form one mode).
#' @param class_of integer per country: the class it belongs to (high
#'   responsibility); default round-robin.
#' @param ll_shift added to every kl_log_lik_ctry draw (a better fit = higher).
#' @param lp per-chain mean of lp__.
make_fake_kl_fit <- function(
  centres,
  countries = sprintf("C%02d", 1:6),
  class_of = NULL,
  n_iter = 50L,
  ll_shift = 0,
  lp = NULL,
  seed = 1L,
  labels = NULL,
  divergent = 0L
) {
  set.seed(seed)
  centres <- as.matrix(centres)
  n_ch <- nrow(centres)
  K <- ncol(centres)
  C <- length(countries)
  class_of <- class_of %||% ((seq_len(C) - 1L) %% K + 1L)
  lp <- lp %||% rep(-1000, n_ch)
  rows <- lapply(seq_len(n_ch), function(ch) {
    d <- list()
    for (k in seq_len(K)) {
      d[[sprintf("kmu[%d]", k)]] <- rnorm(n_iter, centres[ch, k], 0.1)
      d[[sprintf("kcon[%d]", k)]] <- rnorm(n_iter, 0.5 * k, 0.05)
      d[[sprintf("kpi[%d]", k)]] <- rep(1 / K, n_iter)
    }
    d[["kkappa"]] <- abs(rnorm(n_iter, 1, 0.1))
    for (k in seq_len(K)) {
      for (c in seq_len(C)) {
        p <- ifelse(class_of[c] == k, 0.9, 0.1 / max(K - 1L, 1L))
        d[[sprintf("kl_log_resp[%d,%d]", k, c)]] <- rep(log(p), n_iter)
      }
    }
    for (c in seq_len(C)) {
      d[[sprintf("kl_log_lik_ctry[%d]", c)]] <- rnorm(
        n_iter,
        -50 + ll_shift,
        0.5
      )
    }
    d[["lp__"]] <- rnorm(n_iter, lp[ch], 1)
    df <- as.data.frame(d, check.names = FALSE)
    df$.chain <- ch
    df$.iteration <- seq_len(n_iter)
    df
  })
  draws <- posterior::as_draws_df(do.call(rbind, rows))
  structure(
    list(
      draws = draws,
      country_levels = countries,
      class_labels = labels %||% sprintf("klc%d", seq_len(K)),
      K = K,
      chains = n_ch,
      diagnostics = list(
        num_divergent = rep(as.integer(divergent), n_ch),
        num_max_treedepth = rep(0L, n_ch),
        ebfmi = rep(0.9, n_ch)
      ),
      data_dims = list(N = 100L, C = C, J = 3L, NL = 23L)
    ),
    class = "kl_fit"
  )
}
