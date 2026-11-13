# test_duration_ratio_estimand.R -- the RMD-ratio estimand vs the 1/HR shortcut.
#
# durationRatioDraws()/durationRatioSummary() (analysis/utils/dur_post.R) exist
# because the manuscript claim "outbreaks last XX times longer" is a ratio of
# RESTRICTED MEAN DURATIONS, not the inverse cessation hazard ratio 1/HR. The
# two agree only in the limit; reporting 1/HR at a finite horizon is the bug this
# machinery avoids.
#
# This test pins the one closed-form identity that makes the estimand
# falsifiable: under a CONSTANT cessation hazard (no mso_band / time dependence),
#
#   S(t) = exp(-mu t),  h = 1 - exp(-mu),  RMD_K = 1 + sum_{t=1}^{K-1} exp(-mu t)
#   RMD_K -> 1/h  as K -> inf   =>   RMD_hi / RMD_lo -> h_lo/h_hi = 1/HR.
#
# So a correct implementation must (a) reproduce the finite-K RMD closed form
# exactly, (b) converge to 1/HR only as the horizon grows and the tail mass
# vanishes, and (c) DISAGREE with 1/HR at a short horizon -- and honestly warn
# when the unrestricted mean is not identified over the window.
#
# Fully synthetic: a stub "fit" whose posterior_linpred returns a fixed per-month
# log-rate, so mu is controlled directly and the whole surv->rmd->ratio path is
# exercised with a known answer. No data/ access, no fitted model.
#
# Run: Rscript analysis/tests/test_duration_ratio_estimand.R

suppressPackageStartupMessages({
  library(tibble)
})

source("analysis/utils/dur_post.R")   # durationRatioDraws(), durationRatioSummary()

# expect()/fail() use PASS/FAIL; expect_equal()->expect_true() use .passes/.fails.
PASS <- 0L; FAIL <- 0L
.passes <- 0L; .fails <- 0L
source("analysis/tests/helpers.R")

cat("# duration RMD-ratio estimand\n")

# dur_post.R calls brms::posterior_linpred(); the generic must be resolvable for
# S3 dispatch onto our stub. Skip cleanly (status 0) if brms is not installed --
# this is a pure-math identity check, not a modelling-stack gate.
if (!requireNamespace("brms", quietly = TRUE)) {
  cat("SKIP: brms not installed; durationRatioDraws() needs the ",
      "brms::posterior_linpred generic for S3 dispatch\n", sep = "")
  quit(status = 0)
}

# --- stub fit: constant per-month cessation rate, read straight off newdata ---
# One draw, deterministic. durationRatioDraws() forces nd$log_offset to
# log(days_ref) before calling this; we ignore the offset so `mu` is exactly the
# per-month rate we asked for, giving a known closed-form answer.
posterior_linpred.fakedur <- function(object, newdata, re_formula = NA, ...) {
  matrix(log(newdata$mu), nrow = 1L, ncol = nrow(newdata), byrow = TRUE)
}
fit <- structure(list(), class = "fakedur")

nd <- function(mu, K) tibble::tibble(mu = rep(mu, K))

# Independent closed form for the restricted mean duration over K months.
rmd_closed <- function(mu, K) 1 + sum(exp(-mu * seq_len(K - 1L)))

mu_lo <- 0.05     # slow-cessation arm  (long outbreaks)
mu_hi <- 0.20     # fast-cessation arm  (short outbreaks)
h_lo  <- 1 - exp(-mu_lo)
h_hi  <- 1 - exp(-mu_hi)
HR    <- h_hi / h_lo          # cessation hazard ratio, hi vs lo
inv_HR <- 1 / HR              # the shortcut the estimand must NOT be confused with

# --- 1. finite-K RMD matches the closed form exactly -------------------------
cat("\n[1] surv -> rmd path reproduces the closed-form RMD at K = 24\n")
K1 <- 24L
dr1 <- suppressWarnings(
  durationRatioDraws(fit, nd(mu_hi, K1), nd(mu_lo, K1), horizon = K1))

expect_equal(dr1$rmd_hi[1], rmd_closed(mu_hi, K1),
             "RMD(hi) equals 1 + sum exp(-mu t) over 1..K-1", tol = 1e-9)
expect_equal(dr1$rmd_lo[1], rmd_closed(mu_lo, K1),
             "RMD(lo) equals 1 + sum exp(-mu t) over 1..K-1", tol = 1e-9)
expect_equal(dr1$ratio[1], rmd_closed(mu_hi, K1) / rmd_closed(mu_lo, K1),
             "ratio is the per-draw RMD ratio, not a ratio of summaries",
             tol = 1e-9)

# --- 2. RMD ratio -> 1/HR only as the horizon grows --------------------------
cat("\n[2] RMD ratio converges to 1/HR as K -> inf (tail mass -> 0)\n")
Kbig <- 600L
drB  <- durationRatioDraws(fit, nd(mu_hi, Kbig), nd(mu_lo, Kbig), horizon = Kbig)
sumB <- durationRatioSummary(drB)

expect(drB$tail_mass < 1e-3,
       sprintf("tail mass negligible at K = %d (%.2e)", Kbig, drB$tail_mass))
expect_equal(sumB$ratio_med, inv_HR,
             "RMD ratio at long horizon equals 1/HR", tol = 1e-3)

# --- 3. the estimand DISAGREES with 1/HR at a short horizon -------------------
# This is the whole reason the RMD ratio exists: at a finite, truncated horizon
# it is NOT 1/HR. A regression that quietly returned 1/HR would pass [2] and fail
# here.
cat("\n[3] RMD ratio is NOT 1/HR at a short horizon (falsifies the shortcut)\n")
K3  <- 3L
dr3 <- suppressWarnings(
  durationRatioDraws(fit, nd(mu_hi, K3), nd(mu_lo, K3), horizon = K3))
expect(abs(dr3$ratio[1] - inv_HR) > 1e-2,
       sprintf("short-horizon RMD ratio (%.4f) differs from 1/HR (%.4f)",
               dr3$ratio[1], inv_HR))

# --- 4. honest warning when the unrestricted mean is not identified ----------
# At K = 24 with mu_lo = 0.05, median S(K) = exp(-1.2) = 0.30 >> tail_tol, so the
# function MUST warn that this is a restricted mean over the window.
cat("\n[4] warns when median S(K) exceeds tail_tol (restricted, not unrestricted)\n")
warned <- FALSE
invisible(withCallingHandlers(
  durationRatioDraws(fit, nd(mu_hi, K1), nd(mu_lo, K1), horizon = K1),
  warning = function(w) {
    if (grepl("RESTRICTED", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }))
expect(warned, "durationRatioDraws warns that the mean is only restricted here")
expect(dr1$tail_mass > 0.05,
       sprintf("tail mass at K = %d exceeds tail_tol (%.3f > 0.05)",
               K1, dr1$tail_mass))

# --- 5. summary orientation is self-consistent -------------------------------
# mu_hi > mu_lo => the hi arm cesses faster => shorter outbreaks => RMD ratio < 1
# and P(ratio > 1) = 0. A flipped sign convention would trip here.
cat("\n[5] summary orientation: faster hazard => shorter, ratio < 1\n")
expect(sumB$ratio_med < 1,
       "hi-hazard arm yields a shorter restricted mean duration (ratio < 1)")
expect(sumB$prob_longer == 0,
       "prob_longer = P(RMD ratio > 1) is 0 when the hi arm cesses faster")

# --- 6. onset-month average for a fit with duration.seasonal_term harmonics ---
# results_derive.R's dur_rmd_ratio_draws(): a fit without the four harmonic
# columns is scored once, exactly as before; a fit with them is scored for each
# of the 12 onset months (calendar month advancing with months since onset) and
# the RMDs are averaged with equal weight before the per-draw ratio.
# build_dur_newdata() is stubbed in the module's environment so the rows carry
# the known per-month rate `mu` and the harmonic columns.
cat("\n[6] dur_rmd_ratio_draws: unchanged without harmonics, onset-month mean with them\n")
rd <- new.env()
sys.source("analysis/utils/results_derive.R", envir = rd)
a_sin <- 0.8
posterior_linpred.fakeseas <- function(object, newdata, re_formula = NA, ...) {
  matrix(log(newdata$mu) + a_sin * newdata$sin_month, nrow = 1L, ncol = nrow(newdata), byrow = TRUE)
}
rd$build_dur_newdata <- function(fit, vm, oni_value, horizon, mso_to_band) {
  tibble::tibble(mu = rep(if (oni_value == 0) mu_lo else mu_hi, horizon),
                 sin_month = 0, cos_month = 0, sin2_month = 0, cos2_month = 0)
}
K6 <- 12L
plain <- structure(list(data = data.frame(mso_band = "1")), class = "fakedur")
d_new <- suppressWarnings(rd$dur_rmd_ratio_draws(plain, NULL, 1, K6, NULL))
d_old <- suppressWarnings(durationRatioDraws(plain, rd$build_dur_newdata(plain, NULL, 1, K6, NULL),
                                             rd$build_dur_newdata(plain, NULL, 0, K6, NULL), horizon = K6))
expect(identical(d_new, d_old), "a fit without harmonics gives exactly durationRatioDraws() on build_dur_newdata()")

nd11 <- rd$set_dur_onset_month(rd$build_dur_newdata(plain, NULL, 0, 3L, NULL), 11L)
expect(isTRUE(all.equal(nd11$sin_month, sin(2 * pi * c(11, 12, 1) / 12))) &&
       isTRUE(all.equal(nd11$cos2_month, cos(4 * pi * c(11, 12, 1) / 12))),
       "onset month 11: rows sit in calendar months 11, 12, 1")

seas <- structure(list(data = data.frame(sin_month = 0, cos_month = 0, sin2_month = 0, cos2_month = 0)),
                  class = "fakeseas")
rmd_onset <- function(mu, s) {
  mo <- ((s - 1L) + seq_len(K6) - 1L) %% 12L + 1L
  S <- exp(-cumsum(mu * exp(a_sin * sin(2 * pi * mo / 12))))
  1 + sum(S[-K6])
}
want_hi <- mean(vapply(1:12, function(s) rmd_onset(mu_hi, s), numeric(1)))
want_lo <- mean(vapply(1:12, function(s) rmd_onset(mu_lo, s), numeric(1)))
n_warn <- 0L
d_seas <- withCallingHandlers(rd$dur_rmd_ratio_draws(seas, NULL, 1, K6, NULL),
                              warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
expect(abs(d_seas$rmd_hi - want_hi) < 1e-12 && abs(d_seas$rmd_lo - want_lo) < 1e-12,
       sprintf("RMDs are the equal-weight onset-month means (hi %.4f, lo %.4f)", want_hi, want_lo))
expect(abs(d_seas$ratio - want_hi / want_lo) < 1e-12, "ratio is taken on the averaged RMDs")
expect(n_warn == 1L, sprintf("one restricted-mean warning, not one per onset month (got %d)", n_warn))

tot_pass <- PASS + .passes
tot_fail <- FAIL + .fails
cat(sprintf("\n%d pass, %d fail\n", tot_pass, tot_fail))
if (tot_fail > 0) quit(status = 1)
