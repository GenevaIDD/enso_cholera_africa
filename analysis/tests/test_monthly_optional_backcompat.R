# test_monthly_optional_backcompat.R — HARD: the monthly (`mon`) extension of
#   simulate_one_tu() is optional and leaves has_mon = FALSE byte-identical.
#
# Usage:
#   Rscript analysis/tests/test_monthly_optional_backcompat.R
#
# Details:
#   draws_mon = NULL vs omitted gives identical output with and without
#   duration draws (the mon gate is independent of draws_dur); mu_mon /
#   cases_month columns are opt-in; a plain-intercept draws_mon needs no
#   `dlnm` package. test_monthly_offset_wellformed.R covers the has_mon =
#   TRUE wiring; this file covers has_mon = FALSE non-interference.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

source("analysis/utils/dlnm.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")

.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# monthly case-count (mon) optionality / backward-compatibility tests\n\n")

# --------------------------------------------------------------------------
# Fixtures (same shape as test_duration_optional_backcompat.R).
# --------------------------------------------------------------------------

fixture_frame <- function() {
  tibble::tibble(
    gid = "A_1_1",
    country = "AAA",
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = 1L,
    pop = 1e5,
    adm1_unit = "A_1"
  )
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

months <- as.Date(sprintf("2024-%02d-01", 1:3))
horizon <- list(start = "2024-01-01", end = "2024-03-31")

# ==========================================================================
# a) simulate_one_tu() called with draws_mon = NULL explicit vs. draws_mon
#    omitted entirely produce identical event_table/monthly_summary, with
#    has_dur = FALSE (no duration draws at all -- the onset-only shape).
cat(
  "## a) has_dur = FALSE: simulate_one_tu(draws_mon = NULL) vs. draws_mon omitted -> identical\n"
)
{
  frame_tu <- fixture_frame()
  panel <- make_panel(frame_tu, months)
  outdir <- tempfile("mon_backcompat_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -1.5)
  )

  res_explicit_null <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_mon = NULL,
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
    "event_table identical: draws_mon = NULL vs. omitted (has_dur = FALSE)"
  )
  expect_true(
    identical(res_explicit_null$monthly_summary, res_omitted$monthly_summary),
    "monthly_summary identical: draws_mon = NULL vs. omitted (has_dur = FALSE)"
  )
}

# ==========================================================================
# b) Same identity check with has_dur = TRUE (duration draws present,
#    onset forced, cessation suppressed so `active` is actually populated)
#    -- proves the mon block's zero-effect isn't just "dur was already off
#    so nothing ran"; it's specifically gated on draws_mon, independent of
#    draws_dur.
cat(
  "\n## b) has_dur = TRUE, draws_mon = NULL vs. omitted -> still identical (mon gate is independent of dur)\n"
)
{
  frame_tu <- fixture_frame()
  panel <- make_panel(frame_tu, months)
  outdir <- tempfile("mon_backcompat_test_")

  # Forced onset, cessation hazard ~0 -> `active` is populated every month,
  # so this exercises the has_dur = TRUE branch for real (not a degenerate
  # all-inactive run) while still holding draws_mon at NULL/omitted.
  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  res_explicit_null <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = NULL,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "explicit_null"),
    horizon = horizon,
    dry_run = TRUE
  )
  res_omitted <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "omitted"),
    horizon = horizon,
    dry_run = TRUE
  )

  expect_true(
    identical(res_explicit_null$event_table, res_omitted$event_table),
    "event_table identical: draws_mon = NULL vs. omitted (has_dur = TRUE)"
  )
  expect_true(
    identical(res_explicit_null$monthly_summary, res_omitted$monthly_summary),
    "monthly_summary identical: draws_mon = NULL vs. omitted (has_dur = TRUE)"
  )
  expect_true(
    "active" %in%
      names(res_omitted$monthly_summary) &&
      any(res_omitted$monthly_summary$active == 1L, na.rm = TRUE),
    "sanity: duration IS active in this fixture (not a degenerate all-inactive run)"
  )
}

# ==========================================================================
# c) monthly_summary never grows `mu_mon`/`cases_month` columns unless
#    draws_mon is supplied -- these are opt-in additions, not always-present
#    NA-filled columns (which would silently change every existing
#    caller's data.frame shape / column count).
cat(
  "\n## c) mu_mon/cases_month columns are opt-in: absent unless draws_mon is supplied\n"
)
{
  frame_tu <- fixture_frame()
  panel <- make_panel(frame_tu, months)
  outdir <- tempfile("mon_backcompat_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )

  res_no_mon <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "no_mon"),
    horizon = horizon,
    dry_run = TRUE
  )
  expect_true(
    !any(c("mu_mon", "cases_month") %in% names(res_no_mon$monthly_summary)),
    "monthly_summary has no mu_mon/cases_month columns when draws_mon is absent"
  )

  draws_mon <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = log(20)),
    active_terms = c("(Intercept)")
  )
  draws_mon$offset_pop <- TRUE
  draws_mon$offset_logoffset <- TRUE

  res_mon <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    replicate_ids = 1L,
    seed = 42L,
    output_part = file.path(outdir, "with_mon"),
    horizon = horizon,
    dry_run = TRUE
  )
  expect_true(
    all(c("mu_mon", "cases_month") %in% names(res_mon$monthly_summary)),
    "monthly_summary DOES gain mu_mon/cases_month columns once draws_mon is supplied"
  )
  # Size track retired: no `mu` / `Z` columns on any path.
  expect_equal(
    names(res_no_mon$monthly_summary),
    c("replicate_id", "gid", "month", "Y", "active"),
    "dur-only monthly_summary columns are exactly replicate_id/gid/month/Y/active"
  )
  expect_equal(
    names(res_mon$monthly_summary),
    c("replicate_id", "gid", "month", "Y", "active", "mu_mon", "cases_month"),
    "dur+mon monthly_summary columns are exactly ... + mu_mon/cases_month (no mu/Z)"
  )
}

# ==========================================================================
# d) has_mon = TRUE with a plain-intercept draws_mon bundle (no
#    cb_meta_mon) runs without requiring the `dlnm` package -- confirms
#    build_mso_basis()'s NULL short-circuit (utils/dlnm.R) so a config
#    with no monthly cross-basis term never pulls in a dependency the rest
#    of the run does not otherwise need.
cat(
  "\n## d) has_mon = TRUE with a plain-intercept draws_mon (no cb_meta_mon) needs no `dlnm` package\n"
)
{
  # build_mso_basis(cb_meta_mon = NULL) short-circuits to NULL (utils/dlnm.R)
  # before ever touching dlnm::onebasis/crossbasis, so a config with no
  # monthly cross-basis term should run even in an environment where `dlnm`
  # is not installed. We can't uninstall the package for this test process,
  # but we CAN confirm the call path never reaches a dlnm:: symbol by running
  # end-to-end successfully with cb_meta_mon simply never supplied.
  frame_tu <- fixture_frame()
  panel <- make_panel(frame_tu, months)
  outdir <- tempfile("mon_backcompat_test_")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )
  draws_mon <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = log(20)),
    active_terms = c("(Intercept)")
  )
  draws_mon$offset_pop <- TRUE
  draws_mon$offset_logoffset <- TRUE

  out <- tryCatch(
    {
      simulate_one_tu(
        tu = "hmmW",
        frame_tu = frame_tu,
        panel_tu = panel,
        draws_occ = draws_occ,
        draws_dur = draws_dur,
        draws_mon = draws_mon,
        replicate_ids = 1L,
        seed = 42L,
        output_part = file.path(outdir, "plain_mon"),
        horizon = horizon,
        dry_run = TRUE
      )
    },
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )

  expect_true(
    !is.null(out),
    "plain-intercept draws_mon (cb_meta_mon absent) runs to completion"
  )
  expect_true(
    !is.null(out) &&
      all(!is.na(out$monthly_summary$mu_mon[out$monthly_summary$active == 1L])),
    "mu_mon is populated on active months with no cross-basis term at all"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
