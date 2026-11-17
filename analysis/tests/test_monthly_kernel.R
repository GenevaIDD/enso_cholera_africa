# test_monthly_kernel.R — HARD: the mon track's latent months-since-onset
#   cross-basis and the v1 kernel's active-gated monthly emission.
#
# Usage:
#   Rscript analysis/tests/test_monthly_kernel.R
#
# Details:
#   Three checks: dlnm.R::build_mso_basis() bit-identity with a direct
#   onebasis(), spatpomp_core.R::.inject_mso() clamping at the training
#   boundary, and simulator_core.R::simulate_one_tu() filling mu_mon /
#   cases_month on active cells only. Tests 1-2 need the `dlnm` package and
#   skip gracefully without it; test 3 needs no dlnm (plain-intercept
#   draws_mon, no cb_meta_mon).
#
# See: test_monthly_offset_wellformed.R, test_monthly_optional_backcompat.R

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
source("analysis/utils/spatpomp_core.R") # .inject_mso()

.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# test_monthly_kernel.R — mso cross-basis + active-gated mon emission\n\n")

# Strip class/attribute baggage (dlnm's "onebasis" objects carry fun/degree/
# knots/range attrs that differ in storage type -- e.g. numeric vs integer --
# between two otherwise numerically-identical calls) down to a plain numeric
# matrix, so `identical()` compares values only.
.plain_matrix <- function(m) {
  matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m))
}

# ===========================================================================
# 1-2) dlnm-dependent: B_mso bit-identity + clamping (graceful skip)
# ===========================================================================
if (!requireNamespace("dlnm", quietly = TRUE)) {
  cat("[SKIP] dlnm package not available — skipping tests 1-2\n")
} else {
  # 1) B_mso bit-identity: build_mso_basis() reconstructs the exact same
  #    cross-basis table dlnm::onebasis() would produce evaluated directly on
  #    0:mso_max with the same resolved (fit-time) attributes -- proving
  #    build_cb_columns_from_Q()'s replay-verbatim contract holds for the
  #    shipped production spec (bs, degree = 3, df = 4; config_defaults.yml),
  #    not just a toy linear basis.
  cat(
    "## 1) build_mso_basis(): bit-identical to a direct onebasis() re-evaluation\n"
  )
  {
    # Realistic training range (episodes lasting up to ~12 months) and the
    # shipped production spec (config_defaults.yml:515): bs, degree = 3, df = 4.
    x_train <- seq(0, 12, by = 1)
    argvar <- list(fun = "bs", degree = 3L, df = 4L)
    ob <- do.call(dlnm::onebasis, c(list(x = x_train), argvar))

    cb_names <- sprintf("mso_cb%d", seq_len(ncol(ob)))
    colnames(ob) <- cb_names
    ob_attrs <- attributes(ob)
    ob_attrs[c("dim", "dimnames", "class")] <- NULL

    cb_meta_mon <- list(
      mso = list(
        single_lag = 0L,
        onebasis_attrs = ob_attrs,
        argvar = argvar,
        range = range(x_train),
        df = ncol(ob),
        cb_col_names = cb_names
      )
    )

    tab <- build_mso_basis(cb_meta_mon)
    expect_true(
      !is.null(tab),
      "build_mso_basis() returns a table when cb_meta_mon$mso is present"
    )

    mso_max_expected <- max(1L, as.integer(ceiling(max(range(x_train)))))
    expect_equal(
      tab$mso_max,
      mso_max_expected,
      "mso_max = ceiling(max(range)) of the fit-time training data"
    )
    expect_equal(
      dim(tab$B),
      c(mso_max_expected + 1L, ncol(ob)),
      "B has one row per mso in [0, mso_max] and one column per basis df"
    )

    # Direct re-evaluation: replay the exact same onebasis() call
    # build_cb_columns_from_Q() makes internally (resolved attrs minus range).
    ob_args <- ob_attrs[setdiff(names(ob_attrs), "range")]
    direct <- as.matrix(do.call(
      dlnm::onebasis,
      c(list(x = 0:mso_max_expected), ob_args)
    ))

    expect_true(
      identical(.plain_matrix(tab$B), .plain_matrix(direct)),
      "build_mso_basis()'s B is bit-identical to a direct onebasis(0:mso_max) call"
    )
  }

  # 2) Clamping actually clamps: .inject_mso() with an out-of-range mso index
  #    returns the SAME row as the training-boundary index (not a fresh
  #    extrapolation) -- and a positive control shows what an *unclamped*
  #    direct onebasis() evaluation at that same out-of-range x would have
  #    given, to prove the clamp is suppressing a real divergence, not a
  #    no-op. Mirrors the documented CFSv2 ENSO 26x-inflation failure mode
  #    (config_defaults.yml).
  cat(
    "\n## 2) .inject_mso(): out-of-range mso clamps to the training boundary, not a fresh extrapolation\n"
  )
  {
    # Narrow training range so a modest mso (e.g. a 50-month-old episode) is
    # unambiguously out-of-range -- exercises the clamp for real.
    x_train <- seq(0, 3, by = 1)
    argvar <- list(fun = "bs", degree = 3L, df = 4L)
    ob <- do.call(dlnm::onebasis, c(list(x = x_train), argvar))
    cb_names <- sprintf("mso_cb%d", seq_len(ncol(ob)))
    colnames(ob) <- cb_names
    ob_attrs <- attributes(ob)
    ob_attrs[c("dim", "dimnames", "class")] <- NULL

    cb_meta_mon <- list(
      mso = list(
        single_lag = 0L,
        onebasis_attrs = ob_attrs,
        argvar = argvar,
        range = range(x_train),
        df = ncol(ob),
        cb_col_names = cb_names
      )
    )
    tab <- build_mso_basis(cb_meta_mon)
    expect_equal(tab$mso_max, 3L, "narrow fixture: mso_max = 3")

    # In-range sanity: mso_idx = 2 reads row 3 (mso=0 -> row 1) directly, no clamp.
    xm_in <- .inject_mso(list(), tab, 2L)
    expect_true(
      identical(
        unlist(xm_in[cb_names], use.names = FALSE),
        unname(tab$B[3L, ])
      ),
      "in-range mso_idx reads the corresponding row unmodified"
    )

    # Negative mso floors to 0 (row 1).
    xm_neg <- .inject_mso(list(), tab, -5L)
    expect_true(
      identical(
        unlist(xm_neg[cb_names], use.names = FALSE),
        unname(tab$B[1L, ])
      ),
      "negative mso_idx clamps to 0 (row 1)"
    )

    # Out-of-range high mso clamps to mso_max's row (row 4), NOT a fresh
    # extrapolation past the boundary knot.
    xm_out <- .inject_mso(list(), tab, 50L)
    boundary_row <- tab$B[tab$mso_max + 1L, ]
    expect_true(
      identical(
        unlist(xm_out[cb_names], use.names = FALSE),
        unname(boundary_row)
      ),
      "mso_idx = 50 (far out-of-range) reads the SAME row as mso_idx = mso_max (clamped)"
    )

    # Positive control: what would an UNCLAMPED direct onebasis() evaluation
    # at x = 50 have given, using the identical resolved attrs? bs(degree=3)
    # is an unbounded cubic past the boundary knot, so this diverges sharply
    # from the clamped boundary value -- proving the clamp suppresses a real
    # effect, not a no-op. Mirrors the documented CFSv2 26x-inflation case
    # (config_defaults.yml:640-644).
    ob_args <- ob_attrs[setdiff(names(ob_attrs), "range")]
    unclamped <- as.numeric(do.call(dlnm::onebasis, c(list(x = 50), ob_args)))
    divergence <- max(abs(unclamped - as.numeric(boundary_row)))
    cat(sprintf(
      "     unclamped onebasis(50) vs. clamped boundary row: max abs diff = %.2f\n",
      divergence
    ))
    expect_true(
      divergence > 1,
      "unclamped extrapolation at mso=50 diverges substantially from the clamped value (positive control)"
    )
  }
}

# ===========================================================================
# 3) simulate_one_tu(): mu_mon/cases_month populated iff active == 1
# ===========================================================================
# 3) simulate_one_tu(): mu_mon/cases_month populated on every active
#    (loc, month) cell and NA everywhere else -- a two-location fixture
#    (one forced onset + never-cease, one forced never-onset) makes both
#    active = 1 and active = 0 cells deterministic in a single run.
cat(
  "\n## 3) simulate_one_tu(): mu_mon/cases_month present on active months, NA elsewhere\n"
)
{
  # Single-location fixture, run twice with the occurrence intercept forced
  # to opposite extremes (fe_overrides applies to every location in a given
  # call, so two separate calls -- not per-location RE overrides -- is the
  # deterministic way to get one all-active and one all-inactive trajectory;
  # mirrors test_duration_kernel.R's forced-onset/forced-no-onset pattern).
  frame_tu <- tibble::tibble(
    gid = "A_1_1",
    country = "AAA",
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = 1L,
    pop = 1e5,
    adm1_unit = "A_1"
  )

  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  horizon <- list(start = "2024-01-01", end = "2024-03-31")
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
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
  outdir <- tempfile("mon_kernel_test_")

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
  draws_mon$offset_pop <- FALSE
  draws_mon$offset_logoffset <- FALSE

  # Run A: onset forced every eligible month, duration cessation hazard ~0
  # -> onsets month 1 and stays active every month of the horizon.
  draws_occ_on <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = 30)
  )
  res_on <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ_on,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    replicate_ids = 1L,
    seed = 99L,
    output_part = file.path(outdir, "on"),
    horizon = horizon,
    dry_run = TRUE
  )

  # Run B: onset hazard forced ~0 for the entire horizon -> never active.
  draws_occ_off <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -30)
  )
  res_off <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ_off,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    replicate_ids = 1L,
    seed = 99L,
    output_part = file.path(outdir, "off"),
    horizon = horizon,
    dry_run = TRUE
  )

  ms_on <- res_on$monthly_summary
  ms_off <- res_off$monthly_summary
  expect_true(
    all(c("mu_mon", "cases_month") %in% names(ms_on)),
    "monthly_summary carries mu_mon/cases_month once draws_mon is supplied"
  )

  expect_equal(nrow(ms_on), 3L, "3 monthly rows for the forced-onset run")
  expect_equal(nrow(ms_off), 3L, "3 monthly rows for the forced-inactive run")

  expect_true(
    all(ms_on$active == 1L),
    "forced onset + never-cease: active = 1 every month of the horizon"
  )
  expect_true(
    all(ms_off$active == 0L),
    "forced no-onset: active = 0 every month of the horizon"
  )

  expect_true(
    all(!is.na(ms_on$mu_mon)) && all(ms_on$mu_mon > 0),
    "mu_mon is finite and positive on every active (onset run) month"
  )
  expect_true(
    all(!is.na(ms_on$cases_month)),
    "cases_month is populated (non-NA) on every active (onset run) month"
  )

  expect_true(
    all(is.na(ms_off$mu_mon)),
    "mu_mon is NA on every inactive (never-onset run) month"
  )
  expect_true(
    all(is.na(ms_off$cases_month)),
    "cases_month is NA on every inactive (never-onset run) month"
  )

  ms_both <- dplyr::bind_rows(ms_on, ms_off)
  expect_true(
    all((!is.na(ms_both$mu_mon)) == (ms_both$active == 1L)),
    "mu_mon non-NA exactly matches active == 1, across both runs"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
