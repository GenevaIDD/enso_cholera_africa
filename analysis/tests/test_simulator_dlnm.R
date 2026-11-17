# test_simulator_dlnm.R — HARD: DLNM cross-basis support in the simulator
#   kernel (utils/simulator_core.R).
#
# Usage:
#   Rscript analysis/tests/test_simulator_dlnm.R
#
# Details:
#   Tests 29-33: evaluate_lp() consumes _cb* FE columns and nested RE slopes on
#   them, a NULL cb_meta keeps simulate_one_tu() output identical,
#   build_cb_columns_from_Q() and cb_meta_mon. Tests 32-33 require the `dlnm`
#   package and are skipped gracefully when absent.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R") # transitive dependency of dlnm.R (utils/README.md)
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/dlnm.R") # build_cb_columns_from_Q() -- used directly
source("analysis/utils/simulator_core.R")

# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# test_simulator_dlnm.R — DLNM cross-basis support\n\n")

# ---------------------------------------------------------------------------
# Shared fixture: 8-location TU (2 countries × 2 adm1 × 2 adm2)
# ---------------------------------------------------------------------------

fixture_frame <- function() {
  expand.grid(
    country = c("AAA", "BBB"),
    adm1n = c(1, 2),
    adm2n = c(1, 2),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = "REG_W",
      tu = "hmmW",
      min_seas = 1L,
      pop = 1e5
    ) %>%
    tibble::as_tibble()
}

frame_tu <- fixture_frame()
N <- nrow(frame_tu)

# ---------------------------------------------------------------------------
# 29) evaluate_lp: _cb* FE columns consumed identically to flat FE (no dlnm).
cat("\n## 29) evaluate_lp: _cb* FE columns consumed as flat linear terms\n")
# ---------------------------------------------------------------------------
# Verify that cb columns in the FE matrix are detected via the _cb[0-9]+$
# pattern and added to linear_terms, then consumed as (FE + nested_slope) * xv.
# No dlnm needed — xrow cb values are supplied directly.
{
  cb_terms <- c("(Intercept)", "mean_precip_cb1", "mean_precip_cb2")
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(
      `(Intercept)` = 0,
      mean_precip_cb1 = 1.0,
      mean_precip_cb2 = 0.5
    ),
    active_terms = cb_terms
  )

  # Confirm the FE matrix carries the cb columns
  expect_true(
    "mean_precip_cb1" %in% colnames(draws$fe),
    "FE matrix has mean_precip_cb1 column"
  )
  expect_true(
    "mean_precip_cb2" %in% colnames(draws$fe),
    "FE matrix has mean_precip_cb2 column"
  )

  # LP = 0 + 1.0 * 2.0 + 0.5 * 3.0 = 3.5 for all locations
  xrow <- list(mean_precip_cb1 = rep(2.0, N), mean_precip_cb2 = rep(3.0, N))
  lp <- evaluate_lp(
    draws,
    s = 1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  expect_equal(lp, rep(3.5, N), "LP = fe_cb1 * xv_cb1 + fe_cb2 * xv_cb2")

  # Zero xrow → LP = intercept = 0
  xrow0 <- list(mean_precip_cb1 = rep(0, N), mean_precip_cb2 = rep(0, N))
  lp0 <- evaluate_lp(
    draws,
    1L,
    xrow0,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  expect_equal(lp0, rep(0, N), "LP = 0 when cb xrow = 0")

  # Missing cb term in xrow (NULL) → silently skipped, contributes 0
  xrow_partial <- list(mean_precip_cb1 = rep(2.0, N)) # cb2 absent
  lp_partial <- evaluate_lp(
    draws,
    1L,
    xrow_partial,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  expect_equal(
    lp_partial,
    rep(2.0, N),
    "absent cb2 in xrow → contributes 0 (skipped)"
  )
}

# ---------------------------------------------------------------------------
# 30) evaluate_lp: nested RE slope on a cb FE term enters LP (no dlnm).
cat("\n## 30) evaluate_lp: nested ADM1 RE slope on a cb FE term\n")
# ---------------------------------------------------------------------------
# Verify that re_slope on a cb term adds the per-location slope to the LP.
{
  adm1_keys <- re_level_key(frame_tu, "adm1")
  adm1_lv <- unique(adm1_keys)
  aaa1_key <- unique(adm1_keys[frame_tu$adm1_unit == "AAA_1"])

  # ADM1 slope for mean_precip_cb1: +1.5 for AAA_1, 0 elsewhere
  slope_mat <- matrix(
    0,
    nrow = 2,
    ncol = length(adm1_lv),
    dimnames = list(NULL, adm1_lv)
  )
  slope_mat[, aaa1_key] <- 1.5

  cb_terms <- c("(Intercept)", "mean_precip_cb1")
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 0, mean_precip_cb1 = 1.0),
    re_slope_overrides = list(mean_precip_cb1 = list(adm1 = slope_mat)),
    active_terms = cb_terms
  )

  xrow <- list(mean_precip_cb1 = rep(2.0, N))
  lp <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )

  # AAA_1 locations: (1.0 + 1.5) * 2 = 5.0; all others: 1.0 * 2 = 2.0
  exp_lp <- ifelse(frame_tu$adm1_unit == "AAA_1", 5.0, 2.0)
  expect_equal(lp, exp_lp, "RE slope on cb term lifts AAA_1 LP: (fe+slope)*xv")
}

# ---------------------------------------------------------------------------
# 31) simulate_one_tu: backward compat, NULL cb_meta produces identical output.
cat("\n## 31) simulate_one_tu: NULL cb_meta produces identical output\n")
# ---------------------------------------------------------------------------
# Backward compat: passing cb_meta_occ = NULL, cb_meta_mon = NULL explicitly
# must give bit-for-bit identical results to not passing them (default NULL).
{
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0.5,
      dev_precip = 0.1,
      mean_temp_daymean = 0.3,
      dev_temp_daymean = 0.05,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA_character_
    )
  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 4,
    fe_overrides = c(`(Intercept)` = -1)
  )

  # Run A: defaults (cb_meta_occ / cb_meta_mon not passed → NULL)
  resA <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:3L,
    seed = 42L,
    output_part = tempfile("dlnm_bkcompat_a_"),
    horizon = horizon,
    dry_run = TRUE
  )

  # Run B: explicit NULL
  resB <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:3L,
    seed = 42L,
    output_part = tempfile("dlnm_bkcompat_b_"),
    horizon = horizon,
    dry_run = TRUE,
    cb_meta_occ = NULL,
    cb_meta_mon = NULL
  )

  expect_equal(
    resA$event_table,
    resB$event_table,
    "NULL cb_meta: event_table identical to default"
  )
  expect_equal(
    resA$monthly_summary,
    resB$monthly_summary,
    "NULL cb_meta: monthly_summary identical to default"
  )
}

# ---------------------------------------------------------------------------
# dlnm-dependent tests (graceful skip if package absent)
# ---------------------------------------------------------------------------

if (!requireNamespace("dlnm", quietly = TRUE)) {
  cat("\n[SKIP] dlnm package not available — skipping tests 32–33\n")
  cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
  if (.fails > 0L) {
    quit(save = "no", status = 1L)
  }
  quit(save = "no", status = 0L)
}

# ---------------------------------------------------------------------------
# 32) [dlnm] build_cb_columns_from_Q: correct shape and bit-for-bit
#     reproducible.
cat(
  "\n## 32) build_cb_columns_from_Q: correct shape and bit-for-bit reproducible\n"
)
# ---------------------------------------------------------------------------
{
  # Minimal cb_meta for mean_precip: linear variable basis × linear lag basis
  # over lags 0..2.  dlnm::crossbasis with argvar/arglag = list(fun="lin") gives
  # 1 variable df × 1 lag df = 1 cross-basis column.
  meta <- list(
    lag = c(0, 2),
    argvar = list(fun = "lin"),
    arglag = list(fun = "lin"),
    cb_col_names = NA_character_ # placeholder — we'll verify ncol, then set it
  )

  # Build a 5-row Q matrix (5 observations × 3 lags)
  set.seed(7L)
  Q_raw <- matrix(rnorm(15), nrow = 5, ncol = 3)

  # Call once to discover cb_col_names from ncol
  cb_raw <- dlnm::crossbasis(
    Q_raw,
    lag = meta$lag,
    argvar = meta$argvar,
    arglag = meta$arglag
  )
  n_cb <- ncol(as.matrix(cb_raw))
  meta$cb_col_names <- paste0("mean_precip_cb", seq_len(n_cb))

  # Call build_cb_columns_from_Q twice; verify shape and bit-for-bit equality
  out1 <- build_cb_columns_from_Q(Q_raw, meta)
  out2 <- build_cb_columns_from_Q(Q_raw, meta)

  expect_equal(nrow(out1), nrow(Q_raw), "output nrow matches Q nrow")
  expect_equal(ncol(out1), n_cb, "output ncol matches cb_col_names length")
  expect_equal(
    colnames(out1),
    meta$cb_col_names,
    "output colnames match cb_col_names"
  )
  expect_equal(
    out1,
    out2,
    "build_cb_columns_from_Q is bit-for-bit reproducible"
  )

  # Q with different values gives different output (non-trivial transform)
  Q_alt <- Q_raw + 1
  out_alt <- build_cb_columns_from_Q(Q_alt, meta)
  expect_true(
    !isTRUE(all.equal(out1, out_alt, tolerance = 1e-12)),
    "different Q input → different cb output (non-trivial transform)"
  )
}

# ---------------------------------------------------------------------------
# 33) [dlnm] simulate_one_tu: cb_meta_mon activates cross-basis in mon LP.
cat(
  "\n## 33) simulate_one_tu: cb_meta_mon activates cross-basis in the monthly (mon) LP\n"
)
# ---------------------------------------------------------------------------
# Occurrence with a very high intercept (onset in month 1) + a duration bundle
# that never ceases keeps every location active, and the monthly count model
# carries a cross-basis term. With cb_meta_mon provided and FE > 0, the mon LP
# (mu_mon) should be elevated compared to a run without cb_meta_mon.
{
  # Build a minimal cb_meta for mean_precip with max_lag = 2
  meta_precip <- list(
    lag = c(0, 2),
    argvar = list(fun = "lin"),
    arglag = list(fun = "lin"),
    cb_col_names = NA_character_
  )
  # Discover cb_col_names
  Q_probe <- matrix(1, nrow = 3, ncol = 3)
  cb_probe <- dlnm::crossbasis(
    Q_probe,
    lag = meta_precip$lag,
    argvar = meta_precip$argvar,
    arglag = meta_precip$arglag
  )
  n_cb_precip <- ncol(as.matrix(cb_probe))
  meta_precip$cb_col_names <- paste0("mean_precip_cb", seq_len(n_cb_precip))

  cb_meta_mon <- list(mean_precip = meta_precip)

  # Panel: 6 months, all mean_precip = 1.0 (and lag columns = 1.0 too)
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  base_panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::arrange(gid, t_start) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 1.0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA_character_
    )

  # Add lag columns: mean_precip_lag1, mean_precip_lag2 (all 1.0; first k rows NA)
  panel <- base_panel %>%
    dplyr::group_by(gid) %>%
    dplyr::mutate(
      mean_precip_lag1 = dplyr::lag(mean_precip, 1L),
      mean_precip_lag2 = dplyr::lag(mean_precip, 2L)
    ) %>%
    dplyr::ungroup()

  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  # Occurrence: huge positive intercept → onset in month 1 everywhere
  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 20),
    basehaz_const = 0
  )
  # Duration: cessation hazard ~0 → every location stays active
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  # Monthly model A: has mean_precip_cb* FE = 1.0 (elevated LP)
  cb_active_terms <- c("(Intercept)", meta_precip$cb_col_names)
  draws_mon_with_cb <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    phi = 5,
    fe_overrides = setNames(
      c(0, rep(1.0, n_cb_precip)),
      c("(Intercept)", meta_precip$cb_col_names)
    ),
    active_terms = cb_active_terms
  )

  # Monthly model B: no cb terms (FE for them implicitly 0)
  draws_mon_no_cb <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    phi = 5,
    fe_overrides = c(`(Intercept)` = 0)
  )
  for (nm in c("draws_mon_with_cb", "draws_mon_no_cb")) {
    b <- get(nm)
    b$offset_pop <- FALSE
    b$offset_logoffset <- FALSE
    assign(nm, b)
  }

  out_with <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon_with_cb,
    replicate_ids = 1L:2L,
    seed = 77L,
    output_part = tempfile("dlnm_mon_with_"),
    horizon = horizon,
    dry_run = TRUE,
    cb_meta_occ = NULL,
    cb_meta_mon = cb_meta_mon
  )

  out_without <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon_no_cb,
    replicate_ids = 1L:2L,
    seed = 77L,
    output_part = tempfile("dlnm_mon_without_"),
    horizon = horizon,
    dry_run = TRUE,
    cb_meta_occ = NULL,
    cb_meta_mon = NULL
  )

  # Both runs should have events (very high occurrence LP)
  expect_true(
    nrow(out_with$event_table) > 0,
    "cb_meta_mon run: events generated with very high occurrence LP"
  )
  expect_true(
    nrow(out_without$event_table) > 0,
    "no-cb run: events generated with very high occurrence LP"
  )

  # mu_mon is populated on every active month. With mean_precip = 1.0
  # (positive) and FE_cb = 1.0 (positive), the cb contribution is positive →
  # mon LP higher than without cb.
  ms_with <- out_with$monthly_summary
  ms_without <- out_without$monthly_summary

  mu_with <- mean(ms_with$mu_mon[!is.na(ms_with$mu_mon)])
  mu_without <- mean(ms_without$mu_mon[!is.na(ms_without$mu_mon)])

  expect_true(!is.na(mu_with), "cb_meta_mon run: mu_mon populated")
  expect_true(!is.na(mu_without), "no-cb run: mu_mon populated")
  expect_true(
    mu_with > mu_without,
    "cb_meta_mon run: mean mon LP elevated relative to no-cb run"
  )
  # No-cb run: intercept 0, no offsets -> mu_mon = exp(0) = 1 exactly.
  expect_true(
    all(abs(ms_without$mu_mon[!is.na(ms_without$mu_mon)] - 1) < 1e-12),
    "no-cb run: mu_mon = exp(Intercept = 0) exactly (cross-basis did not bleed in)"
  )
}

# ---------------------------------------------------------------------------
cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
