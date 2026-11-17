# test_spatpomp_dunit_detect.R — HARD: the detection-error (`*_state_detect`)
#   dunit measurement family.
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_dunit_detect.R
#
# Details:
#   Run from the project root. Fast and fully deterministic (stochastic
#   recovery assertions stay in test_spatpomp_recovery*.R). Checks enum
#   wiring, frozen emitted-C golden strings
#   (analysis/tests/golden/spatpomp_dunit_golden.rds), builder dispatch, the
#   config plumbing of delta_fn / delta_fp and a bpfilter smoke run.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
  library(tidyr)
})

source("analysis/utils.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
.passes <- 0L
.fails <- 0L

cat("# test_spatpomp_dunit_detect — detection-error dunit family\n\n")

# Why the detect family exists: the eps "state-conditioning" dunit Csnippets
# score the onset match with a single near-delta epsilon (.SP_MEAS_EPS =
# 1e-6). That eps has to be that extreme only because it is silently
# overriding the PROCESS-level mismatch rate implied by the fitted occurrence
# hazard `p` itself -- it conflates aleatory Bernoulli disagreement with
# genuine surveillance detection error, and does almost none of the latter's
# nominal job. The `*_state_detect` siblings use instead an explicit
# two-parameter confusion-matrix likelihood
#
#     P(obs_Y = 1 | Yout = 1) = 1 - delta_fn      P(obs_Y = 1 | Yout = 0) = delta_fp
#
# decoupled from `p`. Because delta_fn/delta_fp are per-run CONFIG parameters,
# the detect Csnippets cannot be pre-built at module load: they are constructed
# at dispatch time inside build_spatpomp_tu_c() from
# .sp_onset_detect_c(delta_fn, delta_fp) spliced into the family's template.
# The mon track's detect sibling (built on the _ACTPEN template) is the only
# count-scoring one.
GOLDEN_PATH <- "analysis/tests/golden/spatpomp_dunit_golden.rds"

# ---------------------------------------------------------------------------
# Fixture — deliberately tiny (2 units x 6 months). Nothing here is filtered
# except in Part 5; Parts 3/4 only need a *built* object so dunit_measure()
# can be evaluated point-wise.
# ---------------------------------------------------------------------------
make_fixture <- function(n_months = 6L, base_date = as.Date("2022-01-01")) {
  frame <- tibble::tibble(country = c("AAA", "AAA"), adm1n = 1:2) %>%
    dplyr::mutate(
      tu = "hmmW",
      gid = sprintf("%s_%d_1", country, adm1n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = "REG_hmmW",
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit)

  months <- base_date %m+% months(0:(n_months - 1L))

  panel <- tidyr::expand_grid(gid = frame$gid, t_start = months) %>%
    dplyr::mutate(
      min_seas = 1L,
      country = sub("_[0-9]+_[0-9]+$", "", gid),
      t_start_seas = as.integer(lubridate::yday(t_start) - 1L),
      mean_precip = 0.0,
      dev_precip = 0.0,
      mean_temp_daymean = 0.0,
      dev_temp_daymean = 0.0,
      outbreak = 0L,
      time_to_last_outbreak_cat = "[6,Inf]",
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L
    ) %>%
    tibble::as_tibble()

  obs <- tibble::tibble(
    gid = rep(frame$gid, each = n_months),
    t_start = rep(months, times = nrow(frame)),
    obs_Y = rep(c(1L, 0L, 0L, 1L, 0L, 0L)[seq_len(n_months)], nrow(frame))
  )

  list(
    frame = frame,
    panel = panel,
    months = months,
    obs = obs,
    window = list(
      start = as.character(months[1]),
      end = as.character(months[n_months] %m+% months(1L))
    )
  )
}

fix <- make_fixture()
draws_occ <- make_draws_bundle_synthetic(
  fix$frame,
  S = 1L,
  fe_overrides = c("(Intercept)" = -4.65),
  basehaz_const = 0.0
)

#' Build a fixture spatPomp object under `measurement` (C builder).
build_fix <- function(measurement, builder = build_spatpomp_tu_c, ...) {
  args <- list(
    fix$frame,
    fix$panel,
    draws_occ,
    fix$window,
    fixed_draw_idx = 1L,
    measurement = measurement,
    obs_data = fix$obs,
    ...
  )
  do.call(builder, args)
}

#' Evaluate dunit_measure() at unit 1 for an explicit (state, obs) point.
#' Everything not named is left at 0, which is what the Csnippets see for the
#' slots they do not read.
eval_dunit <- function(
  obj,
  Yout,
  obs_Y,
  mu = 0,
  size_phi = 5,
  obs_CM = NULL,
  log = TRUE
) {
  U <- length(obj@unit_names)
  sn <- as.vector(t(outer(obj@unit_statenames, seq_len(U), paste0)))
  on <- rownames(obj@data)
  x <- matrix(0, nrow = length(sn), ncol = 1L, dimnames = list(sn, NULL))
  x["Yout1", 1L] <- Yout
  x["mu1", 1L] <- mu
  x["size_phi1", 1L] <- size_phi
  y <- matrix(0, nrow = length(on), ncol = 1L, dimnames = list(on, NULL))
  y["obs_Y1", 1L] <- obs_Y
  if ("obs_CM1" %in% on) {
    y["obs_CM1", 1L] <- obs_CM %||% NA_real_
  }
  # dunit_measure() returns a 1x1 (rep x time) matrix; strip to a bare scalar.
  as.numeric(spatPomp::dunit_measure(
    obj,
    y = y,
    x = x,
    unit = 1L,
    time = 1,
    params = pomp::coef(obj),
    log = log
  ))
}

# ===========================================================================
# 1. Enum wiring
# ===========================================================================
cat("## 1. measurement enum wiring for the detect strings\n")
{
  new_bern <- "occurrence_bernoulli_state_detect"
  old_mon <- "occurrence_duration_mon_NB_state_detect"

  expect_true(
    all(c(new_bern, old_mon) %in% .SP_ALL_MEAS),
    "both detect strings are in .SP_ALL_MEAS"
  )
  expect_true(
    all(c(new_bern, old_mon) %in% .SP_STATE_MEAS),
    "both detect strings are in .SP_STATE_MEAS"
  )
  expect_true(
    all(c(new_bern, old_mon) %in% .SP_DETECT_MEAS),
    "both detect strings are in .SP_DETECT_MEAS"
  )
  expect_true(
    !(new_bern %in% .SP_MON_MEAS) && (old_mon %in% .SP_MON_MEAS),
    "occurrence_bernoulli_state_detect is not a mon measurement; the mon detect string is"
  )

  ok <- tryCatch(
    {
      .check_measurement(new_bern)
      .check_measurement(old_mon)
      TRUE
    },
    error = function(e) FALSE
  )
  expect_true(ok, ".check_measurement() accepts both detect strings")

  expect_true(
    all(vapply(c(new_bern, old_mon), .is_detect_meas, logical(1))),
    ".is_detect_meas() is TRUE for both detect strings"
  )
  expect_true(
    !any(vapply(
      c(
        "occurrence_bernoulli_state",
        "occurrence_duration_mon_NB_state_actpen"
      ),
      .is_detect_meas,
      logical(1)
    )),
    ".is_detect_meas() is FALSE for the corresponding eps siblings"
  )
}

# ===========================================================================
# 2. Golden emitted-C strings
# ===========================================================================
# Regenerating the golden file (analysis/tests/golden/spatpomp_dunit_golden.rds),
# only after you have independently justified the emitted-C change:
#   Rscript -e '
#     source("analysis/utils/spatpomp_core.R")
#     g <- list(
#       eps_state        = as.character(.SP_DUNIT_STATE),
#       eps_mon_nb_state = as.character(.SP_DUNIT_MON_NB_STATE),
#       eps_mon_actpen   = as.character(.SP_DUNIT_MON_NB_STATE_ACTPEN),
#       det_state        = as.character(.make_dunit_state(.sp_onset_detect_c())),
#       det_mon_actpen   = as.character(
#                            .make_dunit_mon_nb_state_actpen(.sp_onset_detect_c())))
#     saveRDS(g, "analysis/tests/golden/spatpomp_dunit_golden.rds")'
cat("\n## 2. emitted C matches the frozen golden strings\n")
{
  current <- list(
    eps_state = as.character(.SP_DUNIT_STATE),
    eps_mon_nb_state = as.character(.SP_DUNIT_MON_NB_STATE),
    eps_mon_actpen = as.character(.SP_DUNIT_MON_NB_STATE_ACTPEN),
    det_state = as.character(.make_dunit_state(.sp_onset_detect_c())),
    det_mon_actpen = as.character(
      .make_dunit_mon_nb_state_actpen(.sp_onset_detect_c())
    )
  )
  if (!file.exists(GOLDEN_PATH)) {
    expect_true(
      FALSE,
      sprintf(
        "golden file missing: %s (see the header for the regeneration snippet)",
        GOLDEN_PATH
      )
    )
  } else {
    golden <- readRDS(GOLDEN_PATH)
    expect_true(
      setequal(names(golden), names(current)),
      "golden file covers exactly the 5 expected dunit variants"
    )
    for (nm in intersect(names(golden), names(current))) {
      expect_true(
        identical(golden[[nm]], current[[nm]]),
        sprintf("emitted C unchanged: %s", nm)
      )
    }
  }
}

# ===========================================================================
# 3. Builder dispatch — functional (C builder)
# ===========================================================================
# Checked by VALUE, not by string matching: at (Yout = 0, obs_Y = 1) the eps
# sibling returns log(eps) = log(1e-6) and the detect sibling log(delta_fp).
cat("\n## 3. build_spatpomp_tu_c() dispatch selects the right dunit\n")
{
  pairs <- list(
    list(
      eps = "occurrence_bernoulli_state",
      det = "occurrence_bernoulli_state_detect"
    )
  )
  for (p in pairs) {
    o_eps <- build_fix(p$eps)
    o_det <- build_fix(p$det)
    v_eps <- eval_dunit(o_eps, Yout = 0, obs_Y = 1)
    v_det <- eval_dunit(o_det, Yout = 0, obs_Y = 1)
    expect_equal(
      v_eps,
      log(.SP_MEAS_EPS),
      sprintf("[C builder] %s scores (Yout=0, obs_Y=1) as log(eps)", p$eps),
      tol = 1e-9
    )
    expect_equal(
      v_det,
      log(.SP_STATE_DELTA_FP),
      sprintf(
        "[C builder] %s scores (Yout=0, obs_Y=1) as log(delta_fp)",
        p$det
      ),
      tol = 1e-9
    )
    expect_equal(
      eval_dunit(o_det, Yout = 1, obs_Y = 0),
      log(.SP_STATE_DELTA_FN),
      sprintf(
        "[C builder] %s scores (Yout=1, obs_Y=0) as log(delta_fn)",
        p$det
      ),
      tol = 1e-9
    )
  }
}

# ===========================================================================
# 4. Config plumbing of delta_fn / delta_fp
# ===========================================================================
cat("\n## 4. delta_fn / delta_fp reach the Csnippet; validation; inertness\n")
{
  obj_d <- build_fix(
    "occurrence_bernoulli_state_detect",
    build_spatpomp_tu_c,
    delta_fn = 0.3,
    delta_fp = 0.2
  )
  expect_equal(
    eval_dunit(obj_d, Yout = 1, obs_Y = 0),
    log(0.3),
    "delta_fn = 0.3 passed to build_spatpomp_tu_c() reaches the compiled Csnippet",
    tol = 1e-9
  )
  expect_equal(
    eval_dunit(obj_d, Yout = 0, obs_Y = 1),
    log(0.2),
    "delta_fp = 0.2 passed to build_spatpomp_tu_c() reaches the compiled Csnippet",
    tol = 1e-9
  )
  expect_true(
    grepl(
      "0.3000",
      as.character(.make_dunit_state(.sp_onset_detect_c(0.3, 0.2))),
      fixed = TRUE
    ),
    "the emitted C literal carries the configured delta_fn at %.4f precision"
  )

  # .check_deltas() rejects everything outside the open interval (0, 1):
  # log(0) / log(<0) would otherwise silently produce -Inf / NaN weights.
  bad <- list(0, 1, -0.1, 1.5, NA_real_, NaN, Inf, c(0.1, 0.2), "0.1", NULL)
  n_rejected <- sum(vapply(
    bad,
    function(v) {
      inherits(tryCatch(.check_deltas(v, 0.05), error = function(e) e), "error")
    },
    logical(1)
  ))
  expect_true(
    n_rejected == length(bad),
    sprintf(
      ".check_deltas() rejects all %d invalid delta_fn values",
      length(bad)
    )
  )
  n_rejected_fp <- sum(vapply(
    bad,
    function(v) {
      inherits(tryCatch(.check_deltas(0.1, v), error = function(e) e), "error")
    },
    logical(1)
  ))
  expect_true(
    n_rejected_fp == length(bad),
    sprintf(
      ".check_deltas() rejects all %d invalid delta_fp values",
      length(bad)
    )
  )
  expect_true(
    isTRUE(.check_deltas(0.1, 0.05)),
    ".check_deltas() accepts (0.10, 0.05)"
  )

  # Inertness: bogus deltas must not affect a non-detect measurement, because
  # .sp_onset_detect_c() (and therefore .check_deltas()) is never reached.
  ok <- tryCatch(
    {
      o <- build_fix(
        "occurrence_bernoulli_state",
        build_spatpomp_tu_c,
        delta_fn = 99,
        delta_fp = -3
      )
      isTRUE(all.equal(eval_dunit(o, Yout = 0, obs_Y = 1), log(.SP_MEAS_EPS)))
    },
    error = function(e) FALSE
  )
  expect_true(
    ok,
    "non-detect measurements ignore delta_fn/delta_fp entirely (no validation, no behaviour change)"
  )
}

# ===========================================================================
# 5. Smoke: build + short bpfilter under the Bernoulli detect measurement
# ===========================================================================
cat("\n## 5. bpfilter smoke run under the Bernoulli detect measurement\n")
{
  for (m in c("occurrence_bernoulli_state_detect")) {
    ll <- tryCatch(
      {
        o <- build_fix(m)
        bl <- make_block_list(o, fix$frame)
        bf <- run_bpfilter_tu(o, Np = 100L, block_list = bl, seed = 20260910L)
        pomp::logLik(bf)
      },
      error = function(e) {
        cat("     ERROR:", conditionMessage(e), "\n")
        NA_real_
      }
    )
    expect_true(
      is.finite(ll),
      sprintf("[C builder] %s: bpfilter logLik finite (%.3f)", m, ll %||% NaN)
    )
  }
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
