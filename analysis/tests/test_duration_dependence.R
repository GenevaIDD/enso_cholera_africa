# test_duration_dependence.R — HARD: the pure-math cores of
#   analysis/utils/dependence.R (size <-> duration posterior dependence).
#
# Usage:
#   Rscript analysis/tests/test_duration_dependence.R
#
# Details:
#   Used by analysis/05_fit_postprocess_dependence.R. A wrong PIT transform does
#   not error: it produces a plausible but wrong Kendall's tau, which feeds the
#   pre-registered decision rule (|tau| < 0.1 and CrI contains 0 => independent
#   draws are adequate for the simulator). A wrong decision either couples
#   draws that don't need it or ships an independent-draws simulator when the
#   tracks co-vary. Fully synthetic: no data/ access, no fitted model.

source("analysis/utils/dependence.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

set.seed(20260831)

cat("== test_duration_dependence ==\n")

# ---------------------------------------------------------------------------
# 1. Size-side PIT (.pit_nb_from_draws) is Uniform(0,1) when the assumed NB
#    (mu, phi) is exactly the process y was simulated from.
# ---------------------------------------------------------------------------
# .pit_nb_from_draws() / .pit_duration_from_draws() are split out of
# pit_size_nb() / pit_duration_episode() so this file can feed them known
# (mu, phi, y) / (lambda, episode_id, mso, ended) triples straight from a
# synthetic generator that shares the exact same discrete-hazard formula (see
# .pit_duration_from_draws()'s own doc comment) and check the textbook result:
# a randomized (censoring-aware) PIT is Uniform(0,1) when the model IS the
# data-generating process.
cat("\n-- NB PIT uniformity\n")

n_obs <- 4000
mu_true <- exp(stats::rnorm(n_obs, log(20), 0.8)) # heterogeneous means, response scale
phi_true <- 5 # constant dispersion (phi_re = "constant" case)
y_sim <- stats::rnbinom(n_obs, size = phi_true, mu = mu_true)

u_nb <- as.vector(.pit_nb_from_draws(
  matrix(mu_true, nrow = 1),
  matrix(phi_true, nrow = 1),
  y_sim
))

expect(all(u_nb >= 0 & u_nb <= 1), "NB PIT values lie in [0, 1]")
ks_nb <- suppressWarnings(stats::ks.test(u_nb, "punif"))
expect(
  ks_nb$p.value > 0.01,
  sprintf("NB PIT is consistent with Uniform(0,1) (KS p = %.3f)", ks_nb$p.value)
)
expect(
  abs(mean(u_nb) - 0.5) < 0.03,
  sprintf("NB PIT mean close to 0.5 (got %.3f)", mean(u_nb))
)

# ---------------------------------------------------------------------------
# 2. Censored duration PIT (.pit_duration_from_draws) is Uniform(0,1),
#    combining uncensored (event) and administratively-censored episodes, when
#    the model's lambda exactly matches the discrete-hazard process the
#    synthetic panel was simulated from (h_k = 1 - exp(-lambda_k) per active
#    month, chained Bernoulli trials, truncated at max_mso).
# ---------------------------------------------------------------------------
cat("\n-- Censored duration PIT uniformity\n")

max_mso <- 12L
lambda_vec <- seq(0.18, 0.05, length.out = max_mso) # per-row TRUE hazard rate (declining)
h_vec <- 1 - exp(-lambda_vec) # implied discrete hazard, q_i

n_ep <- 4000
rows <- vector("list", n_ep)
for (i in seq_len(n_ep)) {
  ended_row <- NA_integer_
  for (k in seq_len(max_mso)) {
    if (stats::runif(1) < h_vec[k]) {
      ended_row <- k
      break
    }
  }
  censored <- is.na(ended_row)
  last_row <- if (censored) max_mso else ended_row
  ended_flag <- rep(0L, last_row)
  if (!censored) {
    ended_flag[last_row] <- 1L
  }
  rows[[i]] <- data.frame(
    episode_id = i,
    mso = seq_len(last_row),
    lambda = lambda_vec[seq_len(last_row)],
    ended = ended_flag
  )
}
panel <- do.call(rbind, rows)

pit_dur <- .pit_duration_from_draws(
  matrix(panel$lambda, nrow = 1),
  panel$episode_id,
  panel$mso,
  panel$ended
)
u_dur <- as.vector(pit_dur$u)

expect(length(u_dur) == n_ep, "one PIT value per episode")
expect(all(u_dur >= 0 & u_dur <= 1), "duration PIT values lie in [0, 1]")
ks_dur <- suppressWarnings(stats::ks.test(u_dur, "punif"))
expect(
  ks_dur$p.value > 0.01,
  sprintf(
    "censored duration PIT is consistent with Uniform(0,1) (KS p = %.3f)",
    ks_dur$p.value
  )
)
expect(
  abs(mean(u_dur) - 0.5) < 0.03,
  sprintf("duration PIT mean close to 0.5 (got %.3f)", mean(u_dur))
)
frac_censored <- mean(pit_dur$ended == 0)
expect(
  frac_censored > 0.05 && frac_censored < 0.95,
  sprintf(
    "synthetic panel mixes censored and uncensored episodes (censored frac = %.2f)",
    frac_censored
  )
)

# ---------------------------------------------------------------------------
# 3. paired_kendall_tau recovers an injected Kendall's tau. Gaussian-copula
#    construction: z2 = rho*z1 + sqrt(1-rho^2)*noise, u = pnorm(z); the
#    resulting Kendall's tau between u1, u2 has the closed form
#    tau = (2/pi) asin(rho), independent of the (here: uniform) margins --
#    exactly the invariance paired_kendall_tau()'s doc comment relies on.
# ---------------------------------------------------------------------------
# paired_kendall_tau() is checked by injecting a known dependence via a
# Gaussian copula (tau = (2/pi) asin(rho) for the copula's latent correlation
# rho) and confirming the recovered posterior mean tau matches, including the
# by_group breakdown and the small-group drop guard (`min_n`).
cat("\n-- tau recovery\n")

rho_true <- 0.6
tau_true <- (2 / pi) * asin(rho_true)

n_pair <- 3000
z1 <- stats::rnorm(n_pair)
z2 <- rho_true * z1 + sqrt(1 - rho_true^2) * stats::rnorm(n_pair)
u1 <- stats::pnorm(z1)
u2 <- stats::pnorm(z2)

ep_ids <- paste0("ep", seq_len(n_pair))
pit_a <- list(u = matrix(u1, nrow = 1), episode_id = ep_ids)
pit_b <- list(u = matrix(u2, nrow = 1), episode_id = ep_ids)

tau_res <- paired_kendall_tau(pit_a, pit_b, min_n = 5L)

expect(tau_res$n_common == n_pair, "all synthetic episodes are paired")
expect(
  abs(tau_res$overall$mean - tau_true) < 0.05,
  sprintf(
    "recovered tau (%.3f) close to injected tau (%.3f) from rho = %.2f",
    tau_res$overall$mean,
    tau_true,
    rho_true
  )
)

# Null case: independent margins should recover tau close to 0.
z1b <- stats::rnorm(n_pair)
z2b <- stats::rnorm(n_pair)
pit_a0 <- list(u = matrix(stats::pnorm(z1b), nrow = 1), episode_id = ep_ids)
pit_b0 <- list(u = matrix(stats::pnorm(z2b), nrow = 1), episode_id = ep_ids)
tau0 <- paired_kendall_tau(pit_a0, pit_b0, min_n = 5L)
expect(
  abs(tau0$overall$mean) < 0.05,
  sprintf(
    "independent margins give tau close to 0 (got %.3f)",
    tau0$overall$mean
  )
)

# By-group breakdown: split into two groups sharing the same injected tau, and
# confirm a group below min_n is dropped rather than erroring.
grp <- setNames(rep(c("A", "B"), length.out = n_pair), ep_ids)
tau_grp <- paired_kendall_tau(pit_a, pit_b, group = grp, min_n = 5L)
expect(
  !is.null(tau_grp$by_group) && all(c("A", "B") %in% names(tau_grp$by_group)),
  "by_group breakdown reports both injected groups"
)
expect(
  all(vapply(
    tau_grp$by_group,
    function(g) abs(g$mean - tau_true) < 0.08,
    logical(1)
  )),
  "each group recovers approximately the same injected tau"
)

grp_small <- setNames(rep(NA_character_, n_pair), ep_ids)
grp_small[1:3] <- "tiny"
tau_small <- suppressMessages(paired_kendall_tau(
  pit_a,
  pit_b,
  group = grp_small,
  min_n = 5L
))
expect(
  length(tau_small$by_group) == 0L,
  "a group below min_n is dropped from by_group rather than erroring"
)

expect(
  inherits(
    try(
      paired_kendall_tau(
        pit_a,
        list(u = pit_b$u, episode_id = ep_ids[1:2]),
        min_n = 5L
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  "fewer common episodes than min_n errors rather than silently reporting a degenerate tau"
)

# ---------------------------------------------------------------------------
# 4. re_intercept_correlation: exact recovery on a hand-built fixture, and the
#    join keeps only (country, adm1_unit) pairs present in both tables.
# ---------------------------------------------------------------------------
# re_intercept_correlation() gets a small hand-built fixture with an exact
# linear relationship.
cat("\n-- RE-intercept correlation\n")

adm1_a <- data.frame(
  country = c("ETH", "ETH", "KEN", "KEN", "COD"),
  adm1_unit = c("A1", "A2", "A1", "A2", "A1"),
  value = c(1, 2, 3, 4, 5)
)
adm1_b <- data.frame(
  country = c("ETH", "ETH", "KEN", "KEN", "COD"),
  adm1_unit = c("A1", "A2", "A1", "A2", "A1"),
  value = c(2, 4, 6, 8, 10)
) # value_b = 2 * value_a exactly

re_corr <- re_intercept_correlation(adm1_a, adm1_b)
expect(re_corr$n == 5L, "all 5 ADM1 units join")
expect(abs(re_corr$r - 1) < 1e-9, "perfectly linear intercepts give r = 1")

re_corr2 <- re_intercept_correlation(adm1_a, adm1_b[1:3, ])
expect(
  re_corr2$n == 3L,
  "join keeps only the (country, adm1_unit) pairs present in both tables"
)

expect(
  inherits(
    try(
      re_intercept_correlation(adm1_a[, c("country", "value")], adm1_b),
      silent = TRUE
    ),
    "try-error"
  ),
  "missing a required column (adm1_unit) errors rather than silently mis-joining"
)

cat(sprintf(
  "\n== test_duration_dependence: %d passed, %d failed ==\n",
  PASS,
  FAIL
))
if (FAIL > 0) {
  quit(status = 1)
}
