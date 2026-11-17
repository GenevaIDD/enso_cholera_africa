# test_duration_kernel.R — HARD: the duration-aware simulator kernel
#   (analysis/utils/simulator_core.R).
#
# Usage:
#   Rscript analysis/tests/test_duration_kernel.R
#
# Details:
#   Covers extract_duration_draws()/.extract_dur_draws_brms() dispatch,
#   basehaz_dur_eval(), the mso_band bin/rename helpers and the has_dur = TRUE
#   branch of simulate_one_tu(): onset gating (no re-onset while active),
#   active/mso bookkeeping, cessation, end-of-horizon censoring, update_lags()
#   fed `active` (not `Y`), and the refusal of a duration fit carrying the
#   retired log(size_cases) covariate. test_duration_optional_backcompat.R pins
#   the has_dur = FALSE path.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# Phase 6 duration-kernel tests\n\n")

# --------------------------------------------------------------------------
# Fixtures (same shape as test_simulator_kernel.R's fixture_frame()).
# --------------------------------------------------------------------------

fixture_frame <- function() {
  expand.grid(
    tu = c("hmmW", "hmmCE"),
    country = c("AAA", "BBB"),
    adm1n = c(1, 2),
    adm2n = c(1, 2),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      country = paste0(country, substr(tu, 4, 4)), # keep countries unique
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", substr(tu, 4, 5)),
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
}

make_panel <- function(frame_tu, months) {
  tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA
    )
}

# ==========================================================================
cat("## 1) .dur_bin_mso: mirrors dur_data.R::bin_mso() boundary semantics\n")
{
  mso <- c(0, 1, 2, 3, 4, 5, 6, 8, 9, 100, NA)
  lv <- .dur_bin_mso(mso)
  expect_true(
    is.na(lv[mso == 0 & !is.na(mso)]),
    "mso=0 (below first break) -> NA"
  )
  expect_equal(lv[!is.na(mso) & mso == 1], "1", "mso=1 -> band '1'")
  expect_equal(lv[!is.na(mso) & mso == 2], "2", "mso=2 -> band '2'")
  expect_equal(lv[!is.na(mso) & mso == 3], "3", "mso=3 -> band '3'")
  expect_equal(lv[!is.na(mso) & mso == 4], "4-5", "mso=4 -> band '4-5'")
  expect_equal(lv[!is.na(mso) & mso == 5], "4-5", "mso=5 -> band '4-5'")
  expect_equal(lv[!is.na(mso) & mso == 6], "6-8", "mso=6 -> band '6-8'")
  expect_equal(lv[!is.na(mso) & mso == 8], "6-8", "mso=8 -> band '6-8'")
  expect_equal(lv[!is.na(mso) & mso == 9], "9+", "mso=9 -> band '9+'")
  expect_equal(lv[!is.na(mso) & mso == 100], "9+", "mso=100 -> band '9+'")
  expect_true(is.na(lv[is.na(mso)]), "NA mso -> NA band")
}

cat("\n## 2) .dur_mso_level_to_col: brms rename() suffix sanitisation\n")
{
  expect_equal(.dur_mso_level_to_col("2"), "mso_band2", "'2' -> mso_band2")
  expect_equal(.dur_mso_level_to_col("3"), "mso_band3", "'3' -> mso_band3")
  expect_equal(
    .dur_mso_level_to_col("4-5"),
    "mso_band4M5",
    "'4-5' -> mso_band4M5"
  )
  expect_equal(
    .dur_mso_level_to_col("6-8"),
    "mso_band6M8",
    "'6-8' -> mso_band6M8"
  )
  expect_equal(.dur_mso_level_to_col("9+"), "mso_band9P", "'9+' -> mso_band9P")
}

cat("\n## 3) basehaz_dur_eval: factor-kind FE lookup by mso band\n")
{
  fe <- matrix(
    c(0, 1.0, 2.0, 3.0, 4.0, 5.0),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "(Intercept)",
        "mso_band2",
        "mso_band3",
        "mso_band4M5",
        "mso_band6M8",
        "mso_band9P"
      )
    )
  )
  draws_dur <- list(
    fe = fe,
    basehaz = list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
  )
  mso <- c(1, 2, 3, 4, 6, 9, NA, 0)
  out <- basehaz_dur_eval(draws_dur, s = 1L, mso = mso, N = length(mso))
  expect_equal(
    out,
    c(0, 1.0, 2.0, 3.0, 4.0, 5.0, 0, 0),
    "reference band '1' + NA + out-of-range mso all contribute 0; others match fe"
  )

  # "none" kind -> always 0, regardless of mso
  draws_none <- list(fe = fe, basehaz = list(kind = "none"))
  out_none <- basehaz_dur_eval(draws_none, s = 1L, mso = mso, N = length(mso))
  expect_true(all(out_none == 0), "kind='none' -> zero contribution")

  # unrecognised kind -> falls back to zero (never errors)
  draws_unk <- list(fe = fe, basehaz = list(kind = "bogus"))
  out_unk <- basehaz_dur_eval(draws_unk, s = 1L, mso = mso, N = length(mso))
  expect_true(
    all(out_unk == 0),
    "unrecognised basehaz kind -> zero fallback, no error"
  )
}

cat(
  "\n## 4) logsize_cases (retired E1 size covariate): duration fits carrying it are refused\n"
)
{
  # The size track is retired downstream: no episode size is
  # drawn, so a duration fit conditioned on log(size_cases) cannot be scored.
  # .extract_dur_draws_brms() must hard-error -- with a message naming the
  # term -- rather than silently scoring it at 0.
  if (requireNamespace("posterior", quietly = TRUE)) {
    dur_df <- data.frame(
      check.names = FALSE,
      b_Intercept = c(-1, -1.1, -0.9),
      b_logsize_cases = c(0.5, 0.4, 0.6)
    )
    err <- tryCatch(
      {
        .extract_dur_draws_brms(dur_df, list(duration = list()))
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    expect_true(
      !is.null(err) && grepl("logsize_cases", err, fixed = TRUE),
      ".extract_dur_draws_brms() refuses a fit with a logsize_cases FE, naming the term"
    )
    expect_true(
      !is.null(err) && grepl("size", err) && grepl("retired", err),
      "the refusal says why (size track retired)"
    )
  } else {
    cat("  [SKIP] posterior not installed -- extractor refusal not exercised\n")
  }

  # The coverage gate does not list it either, so no track can smuggle it in.
  err2 <- tryCatch(
    {
      .assert_lp_terms_consumed(c("(Intercept)", "logsize_cases"), "duration")
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err2) && grepl("logsize_cases", err2, fixed = TRUE),
    ".assert_lp_terms_consumed() rejects logsize_cases (no longer a consumed term)"
  )

  # evaluate_lp() does not score it even if a bundle carries the FE column.
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -1, logsize_cases = 0.5),
    active_terms = c("(Intercept)", "logsize_cases")
  )
  lp <- evaluate_lp(
    draws_dur,
    1L,
    list(logsize_cases = rep(2, N)),
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  expect_equal(lp, rep(-1, N), "evaluate_lp() no longer consumes logsize_cases")
}

cat(
  "\n## 5) extract_duration_draws: unknown backend errors loudly (no silent default)\n"
)
{
  err <- tryCatch(
    extract_duration_draws(
      fit = NULL,
      opt = list(duration = list(backend = "bogus"))
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    grepl("unknown duration backend", err),
    "unrecognised duration.backend raises a clear error before touching `fit`"
  )
}

# ==========================================================================
# simulate_one_tu(): duration-aware path (has_dur = TRUE via non-NULL draws_dur)
# ==========================================================================

cat(
  "\n## 6) simulate_one_tu (duration-aware): forced onset + forced immediate cessation\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date("2024-01-01") # single-month horizon
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-01-31")
  outdir <- tempfile("dur_kernel_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30),
    active_terms = c("(Intercept)")
  )

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "same_month"),
    horizon = horizon,
    dry_run = TRUE
  )

  expect_equal(
    nrow(res$event_table),
    N,
    "every location onsets (forced) -> N event rows"
  )
  expect_true(
    all(res$event_table$duration_months == 1L),
    "cessation drawn (forced) within the onset month itself -> duration_months = 1"
  )
  expect_true(
    all(res$event_table$censored == 0L),
    "episode closed within the horizon -> censored = 0"
  )
  expect_true(
    all(res$event_table$end_month == "2024-01-01"),
    "end_month = the (only) simulated month"
  )
  expect_true(
    all(!is.na(res$event_table$posterior_draw_d)),
    "posterior_draw_d populated for every event row"
  )
  expect_true(
    "active" %in% names(res$monthly_summary),
    "monthly_summary carries the active column when has_dur = TRUE"
  )
  expect_true(
    all(res$monthly_summary$active == 1L),
    "active = 1 for every (loc, month) cell (onset month, all locations)"
  )
}

cat(
  "\n## 7) simulate_one_tu (duration-aware): forced onset, cessation hazard ~0 -> end-of-horizon censoring\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date("2024-01-01")
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-01-31")
  outdir <- tempfile("dur_kernel_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "censored"),
    horizon = horizon,
    dry_run = TRUE
  )

  expect_equal(nrow(res$event_table), N, "every location onsets (forced)")
  expect_true(
    all(res$event_table$censored == 1L),
    "cessation hazard ~0 -> still open at horizon end -> censored = 1"
  )
  expect_true(
    all(res$event_table$duration_months == 1L),
    "1 simulated month elapsed (mso = 1) at administrative censoring"
  )
  expect_true(
    all(res$event_table$end_month == "2024-01-01"),
    "end_month = last simulated month"
  )
}

cat(
  "\n## 8) simulate_one_tu (duration-aware): onset gate -- active location cannot re-onset\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-03-31")
  outdir <- tempfile("dur_kernel_test_")

  # occ hazard forced ~1 every month for every location; dur cessation hazard ~0
  # (episodes never close). Without the onset gate, every location would
  # "re-onset" (and emit a fresh event row) every month regardless.
  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "gate"),
    horizon = horizon,
    dry_run = TRUE
  )

  expect_equal(
    nrow(res$event_table),
    N,
    "exactly one event row per location over the whole horizon (no re-onset)"
  )
  ms <- res$monthly_summary
  expect_true(
    all(ms$Y[ms$month == "2024-01-01"] == 1L),
    "month 1: onset fires for every location"
  )
  expect_true(
    all(ms$Y[ms$month == "2024-02-01"] == 0L),
    "month 2: Y = 0 for every location -- at_risk excluded them (onset gate), not a low-prob miss"
  )
  expect_true(all(ms$Y[ms$month == "2024-03-01"] == 0L), "month 3: same")
  expect_true(
    all(ms$active[ms$month == "2024-01-01"] == 1L),
    "month 1: active = 1"
  )
  expect_true(
    all(ms$active[ms$month == "2024-02-01"] == 1L),
    "month 2: still active (never ceased)"
  )
  expect_true(
    all(ms$active[ms$month == "2024-03-01"] == 1L),
    "month 3: still active"
  )
}

cat(
  "\n## 9) simulate_one_tu (duration-aware): update_lags() is fed `active`, not `Y` (duty-cycle fix)\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-03-31")
  outdir <- tempfile("dur_kernel_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  # Trace every update_lags() call's first argument (the "state fed to
  # contagion" vector) by shadowing it in .GlobalEnv for this test only. Both
  # simulate_one_tu() and this test script are top-level-sourced into
  # .GlobalEnv, so the reassignment is visible to simulate_one_tu()'s free
  # variable lookup at call time.
  captured <- list()
  real_update_lags <- update_lags
  # `...` forwards update_lags()' optional `maps` (kernel speed-up).
  traced_update_lags <- function(state, frame_tu, ...) {
    captured[[length(captured) + 1L]] <<- state
    real_update_lags(state, frame_tu, ...)
  }
  assign("update_lags", traced_update_lags, envir = .GlobalEnv)
  res <- tryCatch(
    simulate_one_tu(
      tu = "hmmW",
      frame_tu = frame_tu,
      panel_tu = panel,
      draws_occ = draws_occ,
      draws_dur = draws_dur,
      replicate_ids = 1L,
      seed = 42L,
      output_part = file.path(outdir, "duty_cycle"),
      horizon = horizon,
      dry_run = TRUE
    ),
    finally = assign("update_lags", real_update_lags, envir = .GlobalEnv)
  )

  expect_equal(
    length(captured),
    3L,
    "update_lags() called once per simulated month"
  )
  # Month 1: every location onsets (forced) -> Y_m1 = active_m1 = all 1s (both
  # candidate feeds agree here, so this month alone can't distinguish the fix).
  expect_true(
    all(captured[[1]] == 1L),
    "month 1 state fed to update_lags: all active/onset"
  )
  # Month 2/3: the onset gate blocks re-onset (Y = 0 everywhere from month 2 on),
  # but every location's episode is still open (dur cessation hazard forced ~0)
  # -> active stays 1. The duty-cycle fix feeds `active` (still all 1s); feeding
  # `Y` would give all 0s from month 2 on, silently switching every contagion
  # lag off the month after onset -- see the @param draws_dur doc on
  # simulate_one_tu().
  expect_true(
    all(captured[[2]] == 1L),
    "month 2: update_lags() fed `active` (=1), not `Y` (=0) -- duty-cycle fix"
  )
  expect_true(
    all(captured[[3]] == 1L),
    "month 3: contagion state remains active for as long as the episode stays open"
  )
}

cat(
  "\n## 10) simulate_one_tu (duration-aware): duration.seasonal_term harmonic2 enters the cessation hazard\n"
)
{
  # Intercept -30 (never ceases) plus 60 * sin_month: cessation certain in the
  # calendar months where sin(2 pi m / 12) > 0 (Jan-May), impossible in Jun-Dec.
  # With forced onsets that is five one-month episodes (Jan..May) and a sixth
  # opened in June, censored at the December horizon after 7 months. The same
  # fit without the harmonic never ceases: one 12-month censored episode.
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- seq(as.Date("2024-01-01"), by = "month", length.out = 12)
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-12-31")
  outdir <- tempfile("dur_kernel_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 30)
  )
  run_h <- function(draws_dur, tag) {
    simulate_one_tu(
      tu = "hmmW",
      frame_tu = frame_tu,
      panel_tu = panel,
      draws_occ = draws_occ,
      draws_dur = draws_dur,
      replicate_ids = 1L,
      seed = 42L,
      output_part = file.path(outdir, tag),
      horizon = horizon,
      dry_run = TRUE
    )
  }

  res_h <- run_h(
    make_draws_bundle_synthetic(
      frame_tu,
      S = 2,
      fe_overrides = c(`(Intercept)` = -30, sin_month = 60),
      active_terms = c("(Intercept)", "sin_month")
    ),
    "harm"
  )
  et <- res_h$event_table
  expect_equal(nrow(et), 6L * N, "harmonic fit: six episodes per location")
  expect_true(
    all(et$duration_months[et$censored == 0L] == 1L) &&
      sum(et$censored == 0L) == 5L * N,
    "harmonic fit: Jan-May episodes cease in their onset month"
  )
  expect_true(
    all(et$duration_months[et$censored == 1L] == 7L) &&
      sum(et$censored == 1L) == N,
    "harmonic fit: the June episode runs to the horizon (7 months, censored)"
  )

  res_0 <- run_h(
    make_draws_bundle_synthetic(
      frame_tu,
      S = 2,
      fe_overrides = c(`(Intercept)` = -30),
      active_terms = "(Intercept)"
    ),
    "noharm"
  )
  expect_true(
    nrow(res_0$event_table) == N &&
      all(res_0$event_table$duration_months == 12L),
    "no harmonic: one 12-month censored episode per location"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
