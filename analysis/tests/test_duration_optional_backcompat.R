# test_duration_optional_backcompat.R — HARD: the duration-aware extension of
#   analysis/utils/simulator_core.R is optional (has_dur = FALSE unchanged).
#
# Usage:
#   Rscript analysis/tests/test_duration_optional_backcompat.R
#
# Details:
#   The onset-only (has_dur = FALSE) code path stays byte-identical: same RNG
#   stream, same event_table/monthly_summary shape. Sections a)-e) below.
#   test_duration_kernel.R covers the has_dur = TRUE branch itself.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/simulations.R")

PASS <- 0L
FAIL <- 0L
.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# Phase 6 optionality / backward-compatibility tests\n\n")

# --------------------------------------------------------------------------
# Fixtures (same shape as test_duration_kernel.R / test_simulator_kernel.R).
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
# build_warmup_state()'s default (has_dur = FALSE) never warns, even when
# panel_tu carries no active/mso columns, and still returns $active/$mso
# (harmless zero-filled additions; no existing caller reads them).
cat(
  "## a) build_warmup_state(has_dur = FALSE, default): no warning, still returns active/mso\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date("2024-01-01")
  panel <- make_panel(frame_tu, months) # no active/mso columns at all
  horizon <- list(start = "2024-01-01", end = "2024-01-31")

  warned <- FALSE
  warm <- withCallingHandlers(
    build_warmup_state(panel, frame_tu, horizon), # has_dur omitted -> default FALSE
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(
    !warned,
    "default has_dur = FALSE: no warning even though panel lacks active/mso"
  )
  expect_true(
    all(c("active", "mso") %in% names(warm)),
    "return value still carries $active/$mso keys (harmless additions)"
  )
  expect_true(
    all(warm$active == 0L) && all(warm$mso == 0L),
    "active/mso default to 0 when panel carries no such columns"
  )
}

# build_warmup_state(has_dur = TRUE) DOES emit the documented warning under the
# same missing-columns condition (a) is only meaningful in contrast to this).
cat(
  "\n## b) build_warmup_state(has_dur = TRUE): DOES warn under the same missing-columns condition\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date("2024-01-01")
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-01-31")

  w <- tryCatch(
    {
      build_warmup_state(panel, frame_tu, horizon, has_dur = TRUE)
      NULL
    },
    warning = function(w) conditionMessage(w)
  )
  expect_true(
    !is.null(w) && grepl("no 'active'/'mso' warm-up columns", w),
    "has_dur = TRUE: warns that panel carries no active/mso warm-up columns"
  )
}

# ==========================================================================
# simulate_one_tu() called with draws_dur = NULL explicitly vs. draws_dur
# omitted entirely produce identical event_table/monthly_summary.
cat(
  "\n## c) simulate_one_tu(): draws_dur = NULL (explicit) vs. omitted -> identical output\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-03-31")
  outdir <- tempfile("dur_backcompat_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -2)
  )

  res_explicit_null <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = NULL,
    replicate_ids = 1L,
    seed = 777L,
    output_part = file.path(outdir, "explicit_null"),
    horizon = horizon,
    dry_run = TRUE
  )
  res_omitted <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    replicate_ids = 1L,
    seed = 777L,
    output_part = file.path(outdir, "omitted"),
    horizon = horizon,
    dry_run = TRUE
  )

  expect_true(
    identical(res_explicit_null$event_table, res_omitted$event_table),
    "event_table identical: draws_dur = NULL vs. omitted"
  )
  expect_true(
    identical(res_explicit_null$monthly_summary, res_omitted$monthly_summary),
    "monthly_summary identical: draws_dur = NULL vs. omitted"
  )
  expect_true(
    !("active" %in% names(res_omitted$monthly_summary)),
    "legacy path monthly_summary never grows an `active` column"
  )
}

# ==========================================================================
# The legacy per-step Bernoulli draw (`Y <- rbinom(N, 1L, p)`, unconditional,
# no onset gate) is untouched: replaying the same L'Ecuyer-CMRG stream and seed
# by hand (RNGkind + set.seed(seed + rep_id), then rbinom() with the
# probabilities computed via the same evaluate_lp()/basehaz_eval() calls the
# kernel itself uses) reproduces monthly_summary$Y exactly, month by month,
# for a multi-month horizon.
cat(
  "\n## d) legacy RNG stream is untouched: manual rbinom() replay matches monthly_summary$Y exactly\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  panel <- make_panel(frame_tu, months)
  horizon <- list(start = "2024-01-01", end = "2024-03-31")
  outdir <- tempfile("dur_backcompat_test_")

  # Moderate (non-forced) intercept so the Bernoulli draw is actually
  # exercised (not saturated at 0 or 1), and S = 1 so the posterior draw
  # index is deterministic regardless of internal clamping arithmetic.
  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -1.5)
  )
  seed <- 2024L
  rep_id <- 1L

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    replicate_ids = rep_id,
    seed = seed,
    output_part = file.path(outdir, "rng_replay"),
    horizon = horizon,
    dry_run = TRUE
  )

  # Manual replay using the kernel's own building blocks (evaluate_lp,
  # basehaz_eval), so this test is robust to internal formula details while
  # still pinning the RNG call sequence itself (RNGkind + set.seed once per
  # replicate, one unconditional rbinom(N, 1, p) per month, no other draws
  # interleaved on the has_dur = FALSE path).
  warm <- build_warmup_state(panel, frame_tu, horizon) # has_dur = FALSE default
  o <- warm$o
  o_ct <- warm$o_ct
  o_tu <- warm$o_tu
  u_days <- warm$u

  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1]), add = TRUE)
  set.seed(seed + rep_id)

  y_expected <- vector("list", length(months))
  log_offset <- log(as.integer(lubridate::days_in_month(months)))
  for (m in seq_along(months)) {
    xrow <- list(
      mean_precip = rep(0, N),
      dev_precip = rep(0, N),
      mean_temp_daymean = rep(0, N),
      dev_temp_daymean = rep(0, N),
      t_start_seas = rep(lubridate::yday(months[m]), N)
    )
    bh_g <- basehaz_eval(draws_occ, xrow$t_start_seas, N)
    log_lambda <- evaluate_lp(
      draws_occ,
      1L,
      xrow,
      o,
      o_ct,
      o_tu,
      u_days,
      frame_tu,
      basehaz_val = bh_g
    )
    lambda <- exp(log_lambda + log_offset[m])
    p <- -expm1(-lambda)
    y_expected[[m]] <- as.integer(stats::rbinom(N, 1L, p))
  }

  ms <- res$monthly_summary
  ok <- TRUE
  for (m in seq_along(months)) {
    y_actual <- ms$Y[ms$month == as.character(months[m])]
    # monthly_summary orders rows by gid (each gid block spans all months in
    # month-order) -- match frame_tu$gid order before comparing.
    y_actual <- y_actual[match(
      frame_tu$gid,
      ms$gid[ms$month == as.character(months[m])]
    )]
    if (!identical(y_actual, y_expected[[m]])) ok <- FALSE
  }
  expect_true(
    ok,
    "manual rbinom() replay (same RNG stream/seed) reproduces monthly_summary$Y for every month"
  )
}

# ==========================================================================
# .check_fit_pair_consistency() (analysis/utils/simulations.R): the
# `opt_dur = NULL` / `opt_mon = NULL` defaults are a no-op (occurrence-only and
# labelled calls never warn); a mismatched duration baseline warns.
cat(
  "\n## e) .check_fit_pair_consistency(): opt_dur/opt_mon = NULL defaults are a no-op\n"
)
{
  opt_occ <- list(baseline = "1990_2019")

  w <- tryCatch(
    {
      .check_fit_pair_consistency(opt_occ) # occurrence-only call
      NULL
    },
    warning = function(w) conditionMessage(w)
  )
  expect_true(
    is.null(w),
    "occurrence-only call (opt_dur/opt_mon omitted): no warning"
  )

  w2 <- tryCatch(
    {
      .check_fit_pair_consistency(opt_occ, label = "grp") # labelled call
      NULL
    },
    warning = function(w) conditionMessage(w)
  )
  expect_true(
    is.null(w2),
    "labelled call (opt_dur/opt_mon omitted): no warning"
  )

  opt_dur_mismatch <- list(baseline = "1980_2009")
  w3 <- tryCatch(
    {
      .check_fit_pair_consistency(opt_occ, opt_dur = opt_dur_mismatch)
      NULL
    },
    warning = function(w) conditionMessage(w)
  )
  expect_true(
    !is.null(w3) && grepl("duration fit", w3),
    "opt_dur supplied with mismatched baseline -> warns (new behaviour, opt-in)"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
