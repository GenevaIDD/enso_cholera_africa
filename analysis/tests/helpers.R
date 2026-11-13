# analysis/tests/helpers.R
#
# Shared test-assertion micro-harness helpers: the single definition of
# `expect()` / `fail()` / `expect_true()` / `expect_equal()` for the test files.
#
# Usage: source("analysis/tests/helpers.R") near the top of a test file, after
# (or before — order doesn't matter, see below) initializing its counters:
#   - expect()      needs  PASS <- 0L; FAIL <- 0L
#   - expect_true() needs  .passes <- 0L; .fails <- 0L
# Counters only need to exist by the time the assertion function is actually
# *called*, not by the time it's defined/sourced (R resolves the `<<-` target
# at call time), so source() and counter-init can appear in either order.
#
# A handful of test files define their own local override of one of these
# functions AFTER sourcing this file, because their assertion semantics
# genuinely differ (verified by content diff during C5, not assumed):
#   - test_size_nested_slopes.R: `expect()` also accumulates a `FAILS` label
#     vector for a end-of-run failure summary.
#   - test_spatpomp_c_variant.R, test_spatpomp_filter_forecast.R,
#     test_spatpomp_grouped_c.R: `expect_equal()` uses looser floating-point
#     tolerance (1e-6 instead of this file's 1e-9 default) and/or different
#     failure-diagnostic verbosity, for numerical checks where exact-bitwise
#     agreement isn't expected.
# R's late binding means a later definition in the same top-level environment
# simply shadows this one for that script's run. Do not "fix" those local
# overrides to match this file — they are intentional, verified-different
# behavior, not accidental duplication.

# --- expect(cond, label) — PASS/FAIL counter style ------------------------
expect <- function(cond, label) {
  if (isTRUE(cond)) {
    PASS <<- PASS + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  FAIL  %s\n", label))
  }
}

# --- fail(msg) — hard-fail gate style (message + quit) --------------------
fail <- function(msg) { message("FAIL: ", msg); quit(status = 1) }

# --- expect_true(cond, label) — .passes/.fails counter style --------------
expect_true <- function(cond, label) {
  if (isTRUE(cond)) {
    .passes <<- .passes + 1L
    cat("  PASS  ", label, "\n")
  } else {
    .fails <<- .fails + 1L
    cat("  FAIL  ", label, "\n")
  }
}

# --- expect_equal(a, b, label, tol) — deparse-diagnostic wrapper ----------
# Delegates to expect_true(); default tol = 1e-9. See file header for the 5
# test files that intentionally keep a local override instead of this one.
expect_equal <- function(a, b, label, tol = 1e-9) {
  ok <- isTRUE(all.equal(a, b, tolerance = tol))
  expect_true(ok, paste0(label,
                          if (!ok) sprintf("  [got=%s want=%s]",
                                            paste(deparse(a), collapse = ""),
                                            paste(deparse(b), collapse = ""))
                          else ""))
}

# --- run_rscript(script, args) — run an R CLI, capture output + status -----
# For tests that drive a CLI end to end (e.g. deploy_prod_run.R --stage ...).
# Returns list(out = stdout+stderr lines, status = exit code).
run_rscript <- function(script, args = character()) {
  out <- suppressWarnings(system2("Rscript", c(script, args),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = out, status = if (is.null(status)) 0L else status)
}
