# test_mon_measurement_enum.R — HARD: config-validation guards for the mon
#   track's `measurement` enum and its hard-require-duration rule.
#
# Usage:
#   Rscript analysis/tests/test_mon_measurement_enum.R
#
# Details:
#   Every check is a "fail loudly at build/setup time" contract: an
#   unguarded predicate chain would silently degrade a mis-wired mon config
#   to the Bernoulli path instead of erroring, a much worse failure than a
#   stop() (a run that "looks superficially fine"). Three independent
#   guards, each at its own natural call site so the fixtures stay minimal
#   (no full spatPomp object build, no C compilation).

suppressPackageStartupMessages({
  library(tibble)
})

source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R

.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# test_mon_measurement_enum — mon-track config-validation guards\n\n")

# ===========================================================================
# 1. .check_measurement() -- unknown value fails loudly, naming the valid set
# ===========================================================================
# 1. .check_measurement() (spatpomp_core.R, used by the C builder): an
#    unknown/typo'd measurement string, and the retired episode-total NegBin
#    (size track) names.
cat("## 1. .check_measurement() rejects an unrecognised measurement string\n")
{
  err <- tryCatch(
    .check_measurement("occurrence_bogus", where = "test"),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.character(err) && grepl("unknown measurement", err),
    "unknown measurement string raises 'unknown measurement'"
  )
  expect_true(
    grepl("occurrence_duration_mon_NB", err),
    "the error names the mon enum value among the supported set"
  )

  ok <- tryCatch(
    {
      .check_measurement("occurrence_duration_mon_NB")
      TRUE
    },
    error = function(e) FALSE
  )
  expect_true(
    ok,
    "the real mon measurement value passes .check_measurement() unchanged"
  )

  # The episode-total NegBin (size track) family is retired with the size
  # track; a stale config naming one must fail here, not degrade silently.
  for (m in c(
    "occurrence_negbin",
    "occurrence_negbin_state",
    "occurrence_negbin_state_detect"
  )) {
    err_m <- tryCatch(
      .check_measurement(m, where = "test"),
      error = function(e) conditionMessage(e)
    )
    expect_true(
      is.character(err_m) && grepl("unknown measurement", err_m),
      sprintf("retired size-track measurement '%s' is rejected", m)
    )
  }
  expect_true(
    !any(grepl("negbin", .SP_ALL_MEAS)),
    ".SP_ALL_MEAS carries no NegBin (size-track) measurement"
  )

  # Issue 8: the normalized `_actn` siblings are valid measurements.
  for (m in c(
    "occurrence_duration_mon_NB_state_detect_actn",
    "occurrence_duration_mon_NB_state_detect_actn_who"
  )) {
    ok_m <- tryCatch(
      {
        .check_measurement(m, where = "test")
        TRUE
      },
      error = function(e) FALSE
    )
    expect_true(
      ok_m && .is_mon_meas(m) && .is_act_meas(m) && .is_actn_meas(m),
      sprintf(
        "'%s' passes .check_measurement() and is a mon / act / actn measurement",
        m
      )
    )
  }
}

# ===========================================================================
# 2. generate_rprocess_csnippet() (spatPomp C generator) -- mon path
#    hard-requires duration, at the Csnippet-generation level so a mis-wired caller fails
#    before emitting C that would not compile / would read uninitialised state.
# ===========================================================================
# 2. generate_rprocess_csnippet() (spatPomp C kernel generator): mon path
#    used without a duration bundle.
cat("\n## 2. generate_rprocess_csnippet(): mon path without duration stops\n")
{
  err <- tryCatch(
    generate_rprocess_csnippet(
      U = 1L,
      memb = NULL,
      cov_bases = NULL,
      measurement = "occurrence_duration_mon_NB",
      use_dur = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.character(err) && grepl("requires use_dur=TRUE", err),
    "generate_rprocess_csnippet: mon measurement + use_dur=FALSE raises 'requires use_dur=TRUE'"
  )
}

# ===========================================================================
# 3. extract_monthly_draws() -- monthly.prev_cases: TRUE is not yet
#    implemented; error before `fit` is ever touched.
# ===========================================================================
# 3. extract_monthly_draws(): monthly.prev_cases: TRUE (not yet
#    implemented; also pinned in test_monthly_offset_wellformed.R -- repeated
#    here so this file alone documents the full mon guard surface).
cat("\n## 3. extract_monthly_draws(): monthly.prev_cases = TRUE stops\n")
{
  err <- tryCatch(
    extract_monthly_draws(
      fit = NULL,
      opt = list(monthly = list(prev_cases = TRUE))
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.character(err) && grepl("prev_cases", err),
    "monthly.prev_cases = TRUE raises a clear error before touching `fit`"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
