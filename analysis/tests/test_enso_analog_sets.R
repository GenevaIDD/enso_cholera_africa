#!/usr/bin/env Rscript
# test_enso_analog_sets.R — HARD: the Lenssen analog-set machinery.
#
# Usage:
#   Rscript analysis/tests/test_enso_analog_sets.R
#
# Details:
#   Covers classify_enso_analog_sets(), enso_analog_set_rows() and the
#   `analog_set` row filter in sample_enso_lenssen_scenarios(). Cases 7 and 10
#   carry most of the weight. Runs without the fitted models or the real
#   Lenssen cache: everything is synthetic except the optional real-cache
#   check at the end (skipped if absent).

suppressPackageStartupMessages({
  library(lubridate) # %m-% used by sample_enso_lenssen_scenarios
})

source("analysis/utils/base.R") # enso_lag_col_name() (call-time dependency of the samplers)
source("analysis/utils/enso.R") # classify_enso_analog_sets()/enso_analog_set_rows()/sample_enso_lenssen_scenarios() -- used directly

pass <- 0L
fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) {
    pass <<- pass + 1L
    cat("  [PASS]", label, "\n")
  } else {
    fail <<- fail + 1L
    cat("  [FAIL]", label, "\n")
  }
}
errs <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

HS <- as.Date("2026-04-01") # handoff; DJF 2026/27 then lands in columns 9:11

# Build a [n x 61] trajectory matrix whose DJF columns carry the given values.
mk_traj <- function(djf_vals, n_months = 61L) {
  m <- matrix(0, nrow = length(djf_vals), ncol = n_months)
  m[, 9:11] <- djf_vals
  m
}
mk_cache <- function(djf_vals, labels = NULL, n_months = 61L) {
  m <- mk_traj(djf_vals, n_months)
  cl <- if (is.null(labels)) {
    classify_enso_analog_sets(m, HS)$labels
  } else {
    labels
  }
  list(
    oni_trajectories = m,
    n_months = ncol(m),
    n_analogs = nrow(m),
    set_labels = cl,
    meta = list(horizon_start = format(HS))
  )
}

# ---------------------------------------------------------------------------
cat("\n== 1. DJF window resolution and set assignment ==\n")
cl <- classify_enso_analog_sets(mk_traj(c(-9, 0, 9)), HS)
ok(
  "labels are nina/neutral/elnino",
  identical(cl$labels, c("nina", "neutral", "elnino"))
)
ok("DJF 2026/27 resolves to cols 9:11", identical(cl$columns, 9:11))
ok("stat_name derived as peak_djf", identical(cl$stat_name, "peak_djf"))
ok(
  "window is Dec-Jan-Feb",
  identical(format(cl$window), c("2026-12-01", "2027-01-01", "2027-02-01"))
)
ok("counts sum to n rows", sum(cl$counts) == 3L)

cat("\n== 2. CPC thresholds are inclusive on both sides ==\n")
cb <- classify_enso_analog_sets(mk_traj(c(-0.5, -0.49, 0.49, 0.5)), HS)
ok("exactly -0.5 -> nina", identical(cb$labels[1L], "nina"))
ok("-0.49 -> neutral", identical(cb$labels[2L], "neutral"))
ok("+0.49 -> neutral", identical(cb$labels[3L], "neutral"))
ok("exactly +0.5 -> elnino", identical(cb$labels[4L], "elnino"))

cat("\n== 3. Truncated / absent peak season must stop(), never fall back ==\n")
e <- errs(classify_enso_analog_sets(matrix(0, 3, 6), HS)) # Apr..Sep only
ok("6-month horizon -> stop", !is.na(e) && grepl("absent or truncated", e))
e <- errs(classify_enso_analog_sets(matrix(0, 3, 10), HS)) # Dec, Jan, no Feb
ok(
  "partial DJF (Dec,Jan only) -> stop",
  !is.na(e) && grepl("absent or truncated", e)
)
e <- errs(classify_enso_analog_sets(
  mk_traj(c(0, 0)),
  HS,
  peak_season_year = 2031L
))
ok(
  "unavailable peak_season_year -> stop",
  !is.na(e) && grepl("absent or truncated", e)
)
c2 <- classify_enso_analog_sets(mk_traj(c(-9, 9)), HS, peak_season_year = 2026L)
ok("explicit peak_season_year 2026 == auto pick", identical(c2$columns, 9:11))

cat("\n== 4. na_action ==\n")
mna <- mk_traj(c(0, 0))
mna[2L, 10L] <- NA_real_
e <- errs(classify_enso_analog_sets(mna, HS))
ok("na_action='error' (default) stops", !is.na(e) && grepl("NA value", e))
cn <- classify_enso_analog_sets(mna, HS, na_action = "na")
ok("na_action='na' propagates NA label", is.na(cn$labels[2L]))
ok("na_action='na' keeps other rows", identical(cn$labels[1L], "neutral"))
ok("NA rows excluded from counts", sum(cn$counts) == 1L)

cat("\n== 5. Argument validation ==\n")
ok(
  "non-matrix traj_mat rejected",
  !is.na(errs(classify_enso_analog_sets(1:10, HS)))
)
ok(
  "nina >= elnino rejected",
  !is.na(errs(classify_enso_analog_sets(
    mk_traj(0),
    HS,
    thresholds = c(nina = 0.5, elnino = -0.5)
  )))
)
ok(
  "unnamed thresholds rejected",
  !is.na(errs(classify_enso_analog_sets(
    mk_traj(0),
    HS,
    thresholds = c(-0.5, 0.5)
  )))
)
ok(
  "peak_months out of 1..12 rejected",
  !is.na(errs(classify_enso_analog_sets(mk_traj(0), HS, peak_months = c(13L))))
)

cat("\n== 6. enso_analog_set_rows ==\n")
cache3 <- mk_cache(c(-9, 0, 9))
ok(
  "'all' returns every row",
  identical(enso_analog_set_rows(cache3, "all"), 1:3)
)
ok("'nina' returns row 1", identical(enso_analog_set_rows(cache3, "nina"), 1L))
ok(
  "'elnino' returns row 3",
  identical(enso_analog_set_rows(cache3, "elnino"), 3L)
)
e <- errs(enso_analog_set_rows(cache3, "nino"))
ok("unknown set name -> stop", !is.na(e) && grepl("unknown analog_set", e))
nolab <- list(oni_trajectories = mk_traj(c(-9, 0, 9)))
ok(
  "'all' works without set_labels",
  identical(enso_analog_set_rows(nolab, "all"), 1:3)
)
e <- errs(enso_analog_set_rows(nolab, "nina"))
ok(
  "missing set_labels -> actionable stop",
  !is.na(e) && grepl("set_labels", e) && grepl("03_prepare_lenssen_analogs", e)
)
e <- errs(enso_analog_set_rows(mk_cache(c(0, 0)), "nina")) # no nina rows exist
ok("empty set -> stop", !is.na(e) && grepl("selects 0 of", e))
bad <- cache3
bad$set_labels <- c("nina", "neutral") # length mismatch
ok(
  "label/row length mismatch -> stop",
  !is.na(errs(enso_analog_set_rows(bad, "nina")))
)

# ---------------------------------------------------------------------------
# Shared fixture for the sampler cases.
hist_dates <- seq(HS - 400, HS, by = "1 month")
HIST <- data.frame(
  year = as.integer(format(hist_dates, "%Y")),
  month = as.integer(format(hist_dates, "%m")),
  enso_index = 0.0
)
SIM_MONTHS <- seq(HS, by = "1 month", length.out = 14L) # covers cols 1..14

draw <- function(cache_path, analog_set = NULL, n_rep = 200L, seed = 42L) {
  a <- list(
    cache_path = cache_path,
    sim_months = SIM_MONTHS,
    n_rep = n_rep,
    historical_oni = HIST,
    max_lag = 12L,
    horizon_start = HS,
    seed = seed
  )
  if (!is.null(analog_set)) {
    a$analog_set <- analog_set
  }
  suppressMessages(do.call(sample_enso_lenssen_scenarios, a))
}

# A synthetic 3-row cache with constant rows at -9 / 0 / +9 proves the row
# filter actually reaches sample.int(). A filter applied in the wrong place
# (e.g. after sampling) would still return plausible numbers, so "the code
# runs" is not evidence the set took effect.
cat("\n== 7. Row filter reaches sample.int (constant-row separation) ==\n")
# Rows are CONSTANT within a set (-9 / 0 / +9 across every column), so any leakage
# from outside the requested set is detectable exactly, not statistically.
sep <- matrix(rep(c(-9, 0, 9), times = 61L), nrow = 3L, ncol = 61L)
sep_cache <- list(
  oni_trajectories = sep,
  n_months = 61L,
  n_analogs = 3L,
  set_labels = c("nina", "neutral", "elnino"),
  meta = list(horizon_start = format(HS))
)
p_sep <- file.path(tempdir(), "sep_cache.rds")
saveRDS(sep_cache, p_sep)

for (spec in list(c("nina", "-9"), c("neutral", "0"), c("elnino", "9"))) {
  s <- spec[1L]
  v <- as.numeric(spec[2L])
  r <- draw(p_sep, s)
  vals <- unique(as.vector(r[["enso_lag0"]][1L, ])) # lead-0, all replicates
  ok(
    sprintf("analog_set='%s' draws only %+g", s, v),
    length(vals) == 1L && isTRUE(all.equal(vals, v))
  )
}
r_all <- draw(p_sep, "all")
ok(
  "analog_set='all' draws from more than one set",
  length(unique(as.vector(r_all[["enso_lag0"]][1L, ]))) > 1L
)

cat("\n== 8. Sampler guards ==\n")
e <- errs(draw(p_sep, "nino"))
ok("unknown analog_set -> stop", !is.na(e) && grepl("unknown analog_set", e))
nolab_cache <- sep_cache
nolab_cache$set_labels <- NULL
p_nolab <- file.path(tempdir(), "nolab_cache.rds")
saveRDS(nolab_cache, p_nolab)
e <- errs(draw(p_nolab, "nina"))
ok(
  "cache without set_labels -> actionable stop",
  !is.na(e) && grepl("set_labels", e) && grepl("03_prepare_lenssen_analogs", e)
)
ok(
  "cache without set_labels still serves 'all'",
  !is.na(draw(p_nolab, "all")[[1L]][1L, 1L])
)
incons <- sep_cache
incons$n_analogs <- 99L
p_inc <- file.path(tempdir(), "inc_cache.rds")
saveRDS(incons, p_inc)
e <- errs(draw(p_inc, "all"))
ok("n_analogs / nrow mismatch -> stop", !is.na(e) && grepl("inconsistent", e))

cat(
  "\n== 9. Set membership survives the affine transform and lag shifting ==\n"
)
# mean_offset must shift the drawn values but must NOT change which rows are eligible.
r_off <- suppressMessages(sample_enso_lenssen_scenarios(
  cache_path = p_sep,
  sim_months = SIM_MONTHS,
  n_rep = 50L,
  historical_oni = HIST,
  max_lag = 12L,
  horizon_start = HS,
  mean_offset = 2.0,
  seed = 42L,
  analog_set = "nina"
))
vals <- unique(as.vector(r_off[["enso_lag0"]][1L, ]))
ok(
  "mean_offset shifts values but keeps the set",
  length(vals) == 1L &&
    isTRUE(all.equal(vals, -9 + 2))
)
# Pre-horizon lag dates come from HIST for every replicate, independent of the set.
r_n <- draw(p_sep, "nina")
ok(
  "pre-horizon lag is set-independent (from historical ONI)",
  all(r_n[["enso_lag12"]][1L, ] == 0.0)
)

# analog_set = "all" is identical() to sampling the unfiltered pool: this is
# what licenses naming the unfiltered counterfactual output trees `_set-all`
# instead of recomputing them.
cat("\n== 10. analog_set='all' is identical to the unfiltered sampler ==\n")
# Same seed, same everything: passing "all", omitting the argument, and the
# pre-change behaviour (unfiltered pool) must agree bit for bit.
pool <- matrix(rnorm(40 * 61), nrow = 40L, ncol = 61L)
lab <- classify_enso_analog_sets(pool, HS, na_action = "na")$labels
p_pool <- file.path(tempdir(), "pool_cache.rds")
saveRDS(
  list(
    oni_trajectories = pool,
    n_months = 61L,
    n_analogs = 40L,
    set_labels = lab,
    meta = list(horizon_start = format(HS))
  ),
  p_pool
)
a_explicit <- draw(p_pool, "all")
a_default <- draw(p_pool, NULL)
ok("'all' == omitting analog_set", identical(a_explicit, a_default))
# A cache with NO set_labels at all is the literal pre-change cache format.
p_pool_nolab <- file.path(tempdir(), "pool_nolab.rds")
saveRDS(
  list(
    oni_trajectories = pool,
    n_months = 61L,
    n_analogs = 40L,
    meta = list(horizon_start = format(HS))
  ),
  p_pool_nolab
)
ok(
  "'all' on a pre-change (unlabelled) cache == labelled cache",
  identical(a_explicit, draw(p_pool_nolab, "all"))
)
# And the sets must actually differ from 'all'.
sub <- draw(p_pool, names(sort(table(lab), decreasing = TRUE))[1L])
ok("a set differs from 'all'", !identical(sub, a_explicit))

cat("\n== 11. Real cache (skipped if absent) ==\n")
REAL <- "generated_data/lenssen_nino34_analogs.rds"
if (!file.exists(REAL)) {
  cat("  [SKIP] no", REAL, "\n")
} else {
  rc <- readRDS(REAL)
  rcl <- classify_enso_analog_sets(rc$oni_trajectories, rc$meta$horizon_start)
  cat(sprintf(
    "  pool: %d analogs; counts: %s\n",
    nrow(rc$oni_trajectories),
    paste(sprintf("%s=%d", names(rcl$counts), rcl$counts), collapse = " ")
  ))
  ok("every analog is classified", sum(rcl$counts) == nrow(rc$oni_trajectories))
  ok("no set is empty", all(rcl$counts > 0L))
  ok(
    "peak_djf discriminates far better than lead 0",
    sd(rcl$stat) > 5 * sd(rc$oni_trajectories[, 1L])
  )
  if (!is.null(rc$set_labels)) {
    ok(
      "cached set_labels match a fresh classification",
      identical(rc$set_labels, rcl$labels)
    )
    ok(
      "cached model provenance is row-aligned",
      is.null(rc$model) || length(rc$model) == nrow(rc$oni_trajectories)
    )
  } else {
    cat(
      "  [SKIP] cache predates set_labels — rebuild with 03_prepare_lenssen_analogs.R\n"
    )
  }
}

cat(sprintf(
  "\n== test_enso_analog_sets: %d passed, %d failed ==\n",
  pass,
  fail
))
quit(status = if (fail > 0L) 1L else 0L)
