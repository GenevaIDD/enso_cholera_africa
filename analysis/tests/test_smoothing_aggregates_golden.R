# test_smoothing_aggregates_golden.R — HARD: golden characterization of the
#   smoothing-compile aggregators in utils/smoothing_core.R.
#
# Usage:
#   Rscript analysis/tests/test_smoothing_aggregates_golden.R
#   Rscript analysis/tests/test_smoothing_aggregates_golden.R --capture
#
# Details:
#   Synthetic data only, no filter run, no data files. smoothing_events_long()
#   / smoothing_events_mon_long() (event tables from smoothed trajectories) and
#   agg_and_ribbon() / _mon() / _active() (weighted ribbons per month x group)
#   feed 08c_smoothing_compile.R and the phase-09 smoothing figures; the golden
#   pins their output on a synthetic smoothing result (with and without the
#   detected feed, empty replicates included). --capture writes the golden.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(purrr)
})
cli <- parse_args(OptionParser(
  option_list = list(
    make_option("--capture", action = "store_true", default = FALSE)
  )
))
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
# smoothing_core.R source()s its whole kernel closure; only its definitions are
# needed here, so evaluate just the aggregator definitions (+ wquantile).
ex <- parse("analysis/utils/smoothing_core.R", keep.source = FALSE)
want <- c(
  "smoothing_events_long",
  "smoothing_events_mon_long",
  "agg_and_ribbon",
  "agg_and_ribbon_mon",
  "agg_and_ribbon_active",
  "wquantile",
  ".smoothing_events_long_impl",
  ".agg_and_ribbon_impl"
)
for (e in ex) {
  if (
    is.call(e) &&
      as.character(e[[1]]) %in% c("<-", "=") &&
      is.name(e[[2]]) &&
      as.character(e[[2]]) %in% want
  ) {
    eval(e)
  }
}

GOLDEN_PATH <- "analysis/tests/golden/smoothing_aggregates_golden.rds"
cat(
  "# smoothing aggregators golden (mode:",
  if (cli$capture) "capture" else "compare",
  ")\n"
)

set.seed(20261006)
U <- 5
M <- 8
months <- seq(as.Date("2019-01-01"), by = "month", length.out = M)
frame_sorted <- tibble(
  gid = sprintf("AAA_%d_1", 1:U),
  AFRO_region = rep(c("WEST", "EAST"), length.out = U),
  country = rep(c("AAA", "BBB", "CCC"), length.out = U)
)
mk_sm <- function(det) {
  reps <- lapply(1:4, function(r) {
    if (r == 3) {
      return(NULL)
    } # a missing replicate
    Y <- matrix(rbinom(U * M, 1, 0.2), U, M)
    if (r == 4) {
      Y[] <- 0L
    } # a replicate with no onset
    cm <- matrix(ifelse(runif(U * M) < 0.3, rpois(U * M, 30), NA), U, M)
    list(
      rep_id = r,
      Yout_mat = Y,
      cases_mon_mat = cm,
      det_mat = if (det) matrix(rbinom(U * M, 1, 0.5), U, M) else NULL
    )
  })
  list(
    reps = reps,
    filt_months = months,
    sorted_gids = frame_sorted$gid,
    frame_sorted = frame_sorted
  )
}
out <- list()
for (det in c(FALSE, TRUE)) {
  sm <- mk_sm(det)
  ev <- smoothing_events_long(sm)
  evm <- smoothing_events_mon_long(sm)
  out[[paste0("events_det", det)]] <- ev
  out[[paste0("events_mon_det", det)]] <- evm
  w <- tibble(rep_id = 1:4, w = c(.4, .3, .2, .1))
  for (gv in c("country", "AFRO_region")) {
    out[[paste("ribbon", det, gv)]] <- agg_and_ribbon(ev, gv, w, months)
    out[[paste("ribbon_mon", det, gv)]] <- agg_and_ribbon_mon(
      evm,
      gv,
      w,
      months
    )
    act <- ev %>% mutate(n_active = 1L + (Yout %% 2L))
    out[[paste("ribbon_active", det, gv)]] <- agg_and_ribbon_active(
      act,
      gv,
      w,
      months
    )
  }
}
out$events_empty <- smoothing_events_long(list(reps = list()))
out$events_mon_empty <- smoothing_events_mon_long(list(
  reps = list(list(rep_id = 1, cases_mon_mat = NULL))
))

if (cli$capture) {
  saveRDS(out, GOLDEN_PATH)
  cat("  captured", length(out), "entries ->", GOLDEN_PATH, "\n")
} else {
  gold <- readRDS(GOLDEN_PATH)
  expect(
    identical(sort(names(gold)), sort(names(out))),
    "same entries as the golden"
  )
  for (nm in names(gold)) {
    ok <- identical(gold[[nm]], out[[nm]])
    expect(ok, sprintf("identical to golden: %s", nm))
    if (!ok) print(all.equal(gold[[nm]], out[[nm]]))
  }
}
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
