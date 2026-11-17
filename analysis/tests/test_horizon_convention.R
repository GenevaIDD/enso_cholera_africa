# test_horizon_convention.R — HARD: one horizon-end convention.
#
# Usage:
#   Rscript analysis/tests/test_horizon_convention.R
#
# Details:
#   Quick, no data. Forward horizons are CLOSED [start, end], filter windows
#   HALF-OPEN [start, end) (utils/base.R); pins the month helpers, the
#   climate-panel end modes, the spatPomp filter window, censoring, the arm
#   launch and panel-coverage asserts against that convention.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/simulator_core.R")
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

# [1] horizon_month_seq() / filter_window_months() semantics
cat("# [1] month helpers\n")
expect(
  identical(
    horizon_month_seq("2026-09-01", "2026-11-01"),
    as.Date(c("2026-09-01", "2026-10-01", "2026-11-01"))
  ),
  "closed horizon includes its end"
)
expect(
  identical(
    horizon_month_seq("2026-09-01", "2026-09-01"),
    as.Date("2026-09-01")
  ),
  "start == end: one month"
)
expect(
  !length(horizon_month_seq("2026-09-01", "2026-08-01")),
  "end < start: empty"
)
expect(
  identical(
    filter_window_months("2026-09-01", "2026-11-01"),
    as.Date(c("2026-09-01", "2026-10-01"))
  ),
  "half-open window excludes its end (the handoff)"
)
expect(
  !length(filter_window_months("2026-09-01", "2026-09-01")),
  "half-open, end == start: empty"
)
expect(
  identical(
    horizon_month_seq("2026-09-15", "2026-10-31"),
    as.Date(c("2026-09-01", "2026-10-01"))
  ),
  "mid-month dates read as their months"
)

# [2] build_climate_panel_for_tu(end_mode =): legacy identical to the frozen
#     legacy slice; closed keeps the end month and warns on a one-month
#     shortfall; open accepts a panel ending at end - 1; the > 31-day gap
#     stays a hard error in every mode
cat("# [2] build_climate_panel_for_tu(end_mode =)\n")
old_slice <- function(panel, gids, start, end, warmup_months = 12L) {
  start_warmup <- as.Date(start) - lubridate::days(round(warmup_months * 30.44))
  panel %>%
    dplyr::filter(as.character(gid) %in% gids) %>%
    dplyr::filter(t_start >= start_warmup, t_start <= as.Date(end))
}
mk_panel <- function(last) {
  m <- seq(as.Date("2024-01-01"), as.Date(last), by = "1 month")
  expand.grid(
    gid = c("A_1", "A_2", "B_1"),
    t_start = m,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      t_start_seas = 1,
      mean_precip = 1,
      dev_precip = 0,
      mean_temp_daymean = 20,
      dev_temp_daymean = 0
    )
}
frame <- data.frame(gid = c("A_1", "A_2"))
hz <- list(start = "2026-01-01", end = "2026-08-01")
full <- mk_panel("2026-12-01")
short <- mk_panel("2026-07-01")
gap <- mk_panel("2026-06-01")
for (mode in c("legacy", "closed", "open")) {
  out <- suppressWarnings(build_climate_panel_for_tu(
    full,
    frame,
    hz,
    end_mode = mode
  ))
  expect(
    identical(out, old_slice(full, frame$gid, hz$start, hz$end)),
    sprintf("%s: slice identical to the frozen pre-item-16 slice", mode)
  )
  out_s <- suppressWarnings(build_climate_panel_for_tu(
    short,
    frame,
    hz,
    end_mode = mode
  ))
  expect(
    identical(out_s, old_slice(short, frame$gid, hz$start, hz$end)),
    sprintf("%s: one-month shortfall, same slice as before", mode)
  )
  expect(
    !is.null(errmsg(build_climate_panel_for_tu(
      gap,
      frame,
      hz,
      end_mode = mode
    ))),
    sprintf("%s: a gap > 1 month stays a hard error", mode)
  )
}
expect(
  identical(
    build_climate_panel_for_tu(full, frame, hz),
    old_slice(full, frame$gid, hz$start, hz$end)
  ),
  "default end_mode is legacy"
)
w <- NULL
withCallingHandlers(
  build_climate_panel_for_tu(short, frame, hz, end_mode = "closed"),
  warning = function(x) {
    w <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  }
)
expect(
  !is.null(w) && grepl("2026-08-01", w, fixed = TRUE),
  "closed: one-month shortfall warns, naming the lost month"
)
w <- NULL
withCallingHandlers(
  build_climate_panel_for_tu(short, frame, hz, end_mode = "legacy"),
  warning = function(x) {
    w <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  }
)
expect(is.null(w), "legacy: no warning (unchanged behaviour)")
w <- NULL
withCallingHandlers(
  build_climate_panel_for_tu(short, frame, hz, end_mode = "open"),
  warning = function(x) {
    w <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  }
)
expect(is.null(w), "open: a window ending at the month after the panel is fine")
expect(
  max(
    build_climate_panel_for_tu(full, frame, hz, end_mode = "closed")$t_start
  ) ==
    as.Date("2026-08-01"),
  "closed: the end month is in the slice"
)

# [3] spatpomp_core.R::.sp_filter_window_months() excludes the window end
cat("# [3] .sp_filter_window_months() is half-open\n")
exprs <- parse("analysis/utils/spatpomp_core.R", keep.source = FALSE)
def <- Filter(
  function(e) {
    is.call(e) &&
      identical(as.character(e[[1]]), "<-") &&
      identical(as.character(e[[2]]), ".sp_filter_window_months")
  },
  as.list(exprs)
)
expect(length(def) == 1L, ".sp_filter_window_months found in spatpomp_core.R")
eval(def[[1]])
fm <- .sp_filter_window_months(
  full,
  list(start = "2026-01-01", end = "2026-05-01"),
  data.frame(tu = "x"),
  "test"
)
expect(
  identical(fm, filter_window_months("2026-01-01", "2026-05-01")),
  "filter-window months == filter_window_months() (end = handoff excluded)"
)

# [4] censor_sim_by_country() keeps the month equal to the horizon
cat("# [4] censor_sim_by_country keeps m == h\n")
df <- data.frame(
  gid = c("AAA.1_1", "AAA.1_1", "AAA.1_1"),
  month = as.Date(c("2026-07-01", "2026-08-01", "2026-09-01")),
  cases = 1:3,
  iso = "AAA"
)
cs <- censor_sim_by_country(df, list(AAA = "2026-08-01"), country_col = "iso")
expect(
  identical(as.Date(cs$month), as.Date(c("2026-07-01", "2026-08-01"))),
  "the horizon month is kept, the month after it censored"
)

# [5] assert_arm_launch() refuses horizon.end != filter.forecast_window.end
cat("# [5] assert_arm_launch: both ends of the closed horizon agree\n")
sim <- list(
  horizon = list(start = "2026-09-01", end = "2027-12-01"),
  filter = list(
    forecast_window = list(start = "2026-09-01", end = "2027-12-01")
  )
)
expect(
  identical(assert_arm_launch(sim, "ok"), as.Date("2026-09-01")),
  "equal ends accepted"
)
sim_bad <- sim
sim_bad$filter$forecast_window$end <- "2028-01-01"
e <- errmsg(assert_arm_launch(sim_bad, "bad"))
expect(
  !is.null(e) && grepl("forecast_window.end", e, fixed = TRUE),
  "differing ends refused"
)
sim_noend <- sim
sim_noend$filter$forecast_window$end <- NULL
expect(
  identical(assert_arm_launch(sim_noend, "noend"), as.Date("2026-09-01")),
  "no forecast_window.end: nothing to compare (unchanged)"
)

# [6] assert_panel_covers_horizon(); enso_seam_month()
cat("# [6] assert_panel_covers_horizon / enso_seam_month\n")
expect(
  isTRUE(assert_panel_covers_horizon(full, "2026-12-01", "t")),
  "panel to its end month: covered"
)
e <- errmsg(assert_panel_covers_horizon(short, "2026-08-01", "t"))
expect(
  !is.null(e) && grepl("2026-08-01", e, fixed = TRUE),
  "one-month shortfall: hard error at config time"
)
expect(
  is.na(suppressMessages(assert_panel_covers_horizon(
    tempfile(fileext = ".rds"),
    "2026-08-01",
    "t"
  ))),
  "absent panel path: not checked"
)
expect(
  identical(enso_seam_month(as.Date("2027-06-01")), as.Date("2027-06-01")) &&
    identical(
      enso_seam_month(as.Date("2027-06-01"), "last_cfsv2"),
      as.Date("2027-06-01")
    ),
  "seam last_cfsv2 (default): the last CFSv2 month, unchanged"
)
expect(
  identical(
    enso_seam_month("2027-06-01", "last_full_oni"),
    as.Date("2027-05-01")
  ) &&
    identical(
      enso_seam_month("2027-01-01", "last_full_oni"),
      as.Date("2026-12-01")
    ),
  "seam last_full_oni: the month before"
)
expect(
  !is.null(errmsg(enso_seam_month("2027-06-01", "bogus"))),
  "unknown seam refused"
)

cat(sprintf("\n[test_horizon_convention] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
