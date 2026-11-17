# test_results_stitched_load_golden.R — HARD: golden characterization of the 09
#   run loader with a smoothing prefix (results_core.R).
#
# Usage:
#   Rscript analysis/tests/test_results_stitched_load_golden.R
#   Rscript analysis/tests/test_results_stitched_load_golden.R --capture
#
# Details:
#   Synthetic temp files only. Covers load_sim_run_monthly() under
#   load_sim_run_stitched() / load_arm_run(factual_source = "smoothing"). Every
#   load must be identical() to analysis/tests/golden/
#   results_stitched_load_golden.rds. --capture writes the golden from the
#   current code (refuses to overwrite).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(qs2)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "gid", "paths", "sim_io", "results_core")) {
  source(sprintf("analysis/utils/%s.R", m))
}

# Guards the fast path of the smoothing prefix: the smoothing rows are
# aggregated per smoothing replicate and the aggregates handed to the forecast
# replicates through handoff_map.csv, instead of a forecast-replicate x gid x
# month grid joined to the events.
# World: smoothing partitions (two TUs, 6 replicates, Jan 2026 - Jun 2027, NA
# and fractional cases); one forecast run of two partitions (KEN, ETH; 20
# replicates; Jan - Jun 2027; handoff maps reusing smoothing replicates). Loads:
# plain, stitched, stitched with a narrower ADM1 window, stitched under
# per-country horizons (drop and flag), smoothing-as-factual (plain and flag).
GOLDEN <- "analysis/tests/golden/results_stitched_load_golden.rds"
capture <- "--capture" %in% commandArgs(trailingOnly = TRUE)

#' stitched_load_world -- the synthetic smoothing + forecast run, under `td`.
stitched_load_world <- function(td, seed = 20261104L) {
  set.seed(seed)
  gids <- list(
    k = c("KEN_1_1", "KEN_1_2", "KEN_2_1"),
    e = c("ETH_1_1", "ETH_2_1")
  )
  sm_months <- seq(as.Date("2026-01-01"), as.Date("2027-06-01"), by = "month")
  sm_parts <- file.path(td, "sm", "parts")
  for (tu in names(gids)) {
    pd <- file.path(sm_parts, paste0("tu=", tu), "chunk_0001")
    dir.create(pd, recursive = TRUE)
    g <- tidyr::expand_grid(rep_id = 1:6, gid = gids[[tu]], t_start = sm_months)
    ev <- g[runif(nrow(g)) < .15, ] |>
      mutate(country = substr(gid, 1, 3), Yout = 1L)
    ev_mon <- g[runif(nrow(g)) < .3, ] |>
      mutate(
        country = substr(gid, 1, 3),
        cases_month = round(rgamma(n(), .6, .02), 1)
      )
    ev_mon$cases_month[sample(nrow(ev_mon), 3)] <- NA
    qs_save(list(events = ev, events_mon = ev_mon), file.path(pd, "part.qs"))
    invisible(file.create(file.path(pd, "_SUCCESS")))
  }
  fmonths <- seq(as.Date("2027-01-01"), as.Date("2027-06-01"), by = "month")
  for (tu in names(gids)) {
    pd <- file.path(td, "sim", "fc", "parts", paste0("tu=", tu))
    dir.create(pd, recursive = TRUE)
    ms <- tidyr::expand_grid(
      replicate_id = 1:20,
      gid = gids[[tu]],
      month = format(fmonths)
    ) |>
      mutate(
        Y = rbinom(n(), 1, .2),
        cases_month = ifelse(
          runif(n()) < .4,
          round(rgamma(n(), .6, .02), 1),
          NA
        )
      )
    qs_save(ms, file.path(pd, "monthly_summary.qs"))
    invisible(file.create(file.path(pd, "_SUCCESS")))
    utils::write.csv(
      data.frame(
        tu = tu,
        replicate_id = 1:20,
        smoothing_rep_id = sample(1:6, 20, TRUE),
        handoff_month = "2027-01-01",
        handoff_weighting = "uniform"
      ),
      file.path(pd, "handoff_map.csv"),
      row.names = FALSE
    )
  }
  list(
    sim = file.path(td, "sim"),
    sev = read_smoothing_parts_events(sm_parts, as.Date("2026-01-01"))
  )
}

#' stitched_loads -- every load the golden pins (stitch$parts$part is the
#' temp-dir dependent partition path, dropped).
stitched_loads <- function(W) {
  hz <- list(KEN = "2027-03-01", ETH = "2027-06-01")
  sm <- function(...) {
    suppressMessages(load_arm_run(
      W$sim,
      "fc",
      "forecast",
      ...,
      smoothing_events = W$sev,
      factual_source = "smoothing"
    ))
  }
  st <- function(...) {
    suppressMessages(load_sim_run_stitched(
      W$sim,
      "fc",
      ...,
      factual_start = "2027-01-01",
      smoothing_events = W$sev
    ))
  }
  out <- list(
    plain = load_sim_run_monthly(W$sim, "fc", "2027-01-01", "2027-06-01"),
    stitched = st("2026-09-01", "2027-06-01"),
    stitched_win = st("2026-09-01", "2027-04-01"),
    stitched_hz = st("2026-09-01", "2027-06-01", horizon_by_country = hz),
    stitched_fl = st(
      "2026-09-01",
      "2027-06-01",
      horizon_by_country = hz,
      keep_beyond = TRUE
    ),
    smoothing = sm("2027-01-01", "2027-06-01"),
    smoothing_fl = sm(
      "2027-01-01",
      "2027-06-01",
      horizon_by_country = hz,
      keep_beyond = TRUE
    )
  )
  lapply(out, function(x) {
    if (!is.null(x$stitch$parts)) {
      x$stitch$parts$part <- NULL
    }
    x
  })
}

if (sys.nframe() == 0L) {
  td <- tempfile("stload_")
  dir.create(td)
  L <- stitched_loads(stitched_load_world(td))
  if (capture) {
    if (file.exists(GOLDEN)) {
      stop("golden exists, not overwritten: ", GOLDEN, call. = FALSE)
    }
    saveRDS(L, GOLDEN)
    cat("captured", GOLDEN, "\n")
    quit(status = 0)
  }
  cat("# stitched / smoothing-as-factual loads vs the pre-speed-up capture\n")
  gold <- readRDS(GOLDEN)
  expect(identical(names(gold), names(L)), "same set of loads")
  for (nm in names(gold)) {
    ok <- identical(gold[[nm]], L[[nm]])
    expect(ok, sprintf("%s: identical to the golden", nm))
    if (!ok) print(all.equal(gold[[nm]], L[[nm]], tolerance = 0))
  }
  expect(
    sum(L$stitched$country_month$month < as.Date("2027-01-01")) ==
      20L * 2L * 4L,
    "prefix is a full replicate x country x month grid (zero-filled)"
  )
  unlink(td, recursive = TRUE)
  cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
  if (FAIL > 0) quit(status = 1)
}
