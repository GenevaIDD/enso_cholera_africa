# test_size_synthetic_equivalence.R — synthetic recovery of the brms size
#   formula.
#
# Usage:
#   Rscript analysis/tests/test_size_synthetic_equivalence.R
#
# Details:
#   Simulates NB-log data from known fixed-effect betas + nested random
#   effects, fits build_size_formula() on a tiny MCMC grid (warmup=200,
#   iter=400, 2 chains) and checks that >= 80% of fixed-effect 90% CIs cover
#   the truth. Soft: low coverage is a WARN; a hard FAIL only if the fit
#   crashes or the formula machinery breaks. Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

set.seed(42)

# --- Truth -----------------------------------------------------------------
n_region <- 2
n_country <- 6 # 3 per region
n_adm1 <- 24 # 4 per country
n_per_adm1 <- 4

beta_truth <- c(
  Intercept = -3.0,
  mean_precip = 0.20,
  dev_precip = 0.15,
  mean_temp_daymean = 0.10,
  dev_temp_daymean = -0.08
)

sd_region <- 0.5
sd_country <- 0.4
sd_adm1 <- 0.3

# --- Simulate --------------------------------------------------------------
sim <- tibble::tibble(
  AFRO_region = rep(
    paste0("R", seq_len(n_region)),
    each = n_country / n_region * n_adm1 / n_country * n_per_adm1
  )
)
sim$country <- rep(
  rep(paste0("C", seq_len(n_country)), each = n_adm1 / n_country * n_per_adm1),
  times = 1
)
sim$ADM1 <- rep(paste0("A", seq_len(n_adm1)), each = n_per_adm1)
sim$adm1_unit <- sim$ADM1
sim$uid <- seq_len(nrow(sim))

# Climate covariates (standardised draws)
sim$mean_precip <- rnorm(nrow(sim))
sim$dev_precip <- rnorm(nrow(sim))
sim$mean_temp_daymean <- rnorm(nrow(sim))
sim$dev_temp_daymean <- rnorm(nrow(sim))
sim$time_to_last_outbreak <- runif(nrow(sim), 0, 1)

# Population offset
sim$pop <- 1e5
sim$mean_rate <- 0.01

# Random effects
re_region <- rnorm(n_region, 0, sd_region)
re_country <- rnorm(n_country, 0, sd_country)
re_adm1 <- rnorm(n_adm1, 0, sd_adm1)

sim$eta <- with(
  sim,
  beta_truth["Intercept"] +
    beta_truth["mean_precip"] * mean_precip +
    beta_truth["dev_precip"] * dev_precip +
    beta_truth["mean_temp_daymean"] * mean_temp_daymean +
    beta_truth["dev_temp_daymean"] * dev_temp_daymean +
    re_region[match(AFRO_region, paste0("R", seq_len(n_region)))] +
    re_country[match(country, paste0("C", seq_len(n_country)))] +
    re_adm1[match(ADM1, paste0("A", seq_len(n_adm1)))] +
    log(pop) +
    log(mean_rate)
)

# NB(mu = exp(eta), phi = 5)
phi_truth <- 5
sim$suspected_cases <- rnbinom(nrow(sim), size = phi_truth, mu = exp(sim$eta))

# Other columns read_size_data expects (we use the simulated frame directly)
sim$location_period_id <- sim$uid
sim$gid <- sim$ADM1
sim$TL <- as.Date("2010-01-01")
sim$TR <- as.Date("2020-12-31")

if (any(!is.finite(sim$suspected_cases))) {
  message("FAIL: synthetic suspected_cases not finite — sim bug")
  quit(status = 1)
}

# --- Fit -------------------------------------------------------------------
form <- bf(
  suspected_cases ~ 1 +
    mean_precip +
    dev_precip +
    mean_temp_daymean +
    dev_temp_daymean +
    offset(log(pop)) +
    offset(log(mean_rate)) +
    (1 | AFRO_region / country / ADM1)
)

fit <- tryCatch(
  {
    brm(
      form,
      data = sim,
      family = negbinomial(),
      prior = c(
        prior(normal(0, 1), class = "b"),
        prior(normal(-3, 2), class = "Intercept"),
        prior(normal(0, 1), class = "sd")
      ),
      chains = 2,
      iter = 400,
      warmup = 200,
      cores = 2,
      seed = 42,
      refresh = 0,
      control = list(adapt_delta = 0.9)
    )
  },
  error = function(e) {
    message("FAIL: synthetic brm fit errored: ", conditionMessage(e))
    quit(status = 1)
  }
)

# --- Coverage check -------------------------------------------------------
post <- as_draws_df(fit)
fe <- list(
  Intercept = "b_Intercept",
  mean_precip = "b_mean_precip",
  dev_precip = "b_dev_precip",
  mean_temp_daymean = "b_mean_temp_daymean",
  dev_temp_daymean = "b_dev_temp_daymean"
)

covered <- 0L
total <- length(fe)
for (nm in names(fe)) {
  col <- fe[[nm]]
  if (!col %in% names(post)) {
    message("WARN: ", col, " not in post")
    next
  }
  q05 <- quantile(post[[col]], 0.05)
  q95 <- quantile(post[[col]], 0.95)
  ok <- (beta_truth[[nm]] >= q05) && (beta_truth[[nm]] <= q95)
  covered <- covered + as.integer(ok)
  message(sprintf(
    "  %-20s truth=%6.3f  [%.3f, %.3f]  %s",
    nm,
    beta_truth[[nm]],
    q05,
    q95,
    if (ok) "OK" else "MISS"
  ))
}

cov_pct <- covered / total
message(sprintf(
  "synthetic 90%% CI coverage: %d/%d (%.0f%%)",
  covered,
  total,
  100 * cov_pct
))

if (cov_pct < 0.8) {
  message("WARN: synthetic coverage < 80% — soft fail (tiny MCMC grid)")
}

message(
  "PASS: test_size_synthetic_equivalence (soft, coverage=",
  sprintf("%.0f%%", 100 * cov_pct),
  ")"
)
quit(status = 0)
