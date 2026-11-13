# analysis/tests/test_climate_zero_fill_guard.R
#
# HARD gate on the zero-filled climate guard:
# simulator_core.R::.assert_no_zero_filled_climate(), with its helpers
# .draws_fe_names() and .sim_panel_window().
#
# 03_prepare_sim_panel_projection.R sets the four local-climate columns
# (mean_precip, dev_precip, mean_temp_daymean, dev_temp_daymean) to exactly 0
# past the climate DB's last month. The kernels' forward-fills cover only
# months with no panel row, so a fit with a climate fixed effect would be
# forecast at zero rainfall / temperature without any error. The guard refuses
# that case and is inert for fits without a climate term (all of the deployed fits).
# It is called by utils/simulations.R::load_artefacts() (both loaders),
# 07b_simulator_v1.R's own loader and utils/smoothing_core.R.
#
# Synthetic panels and fake draws bundles only: no fit, no data/ access.
#
# Run: Rscript analysis/tests/test_climate_zero_fill_guard.R

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/base.R")
source("analysis/utils/simulator_core.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# zero-filled climate guard\n")

refuses <- function(expr) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(msg)
}

# Panel: 2023-07 .. 2024-06 for two gids; observed climate up to 2023-12,
# zero-filled after it (the projection script's convention).
months <- seq(as.Date("2023-07-01"), as.Date("2024-06-01"), by = "month")
panel  <- expand.grid(gid = c("KEN.1_1", "KEN.2_1"), t_start = months,
                      stringsAsFactors = FALSE)
post   <- panel$t_start > as.Date("2023-12-01")
set.seed(1)
panel$mean_precip       <- ifelse(post, 0, runif(nrow(panel), 10, 200))
panel$dev_precip        <- ifelse(post, 0, rnorm(nrow(panel)))
panel$mean_temp_daymean <- ifelse(post, 0, runif(nrow(panel), 18, 30))
panel$dev_temp_daymean  <- ifelse(post, 0, rnorm(nrow(panel)))
panel$enso_lag0         <- rnorm(nrow(panel))

fe_bundle <- function(cols) list(fe = matrix(0, nrow = 2, ncol = length(cols),
                                             dimnames = list(NULL, cols)))

deployed <- c("Intercept", "time_to_last_outbreak1", "ongoing_outbreak_lag1",
              "nbr_lag1_r500_rn_parthmmclusters", paste0("enso_cb", 1:5))
with_clim <- c(deployed, "mean_temp_daymean")

# --- [1] inert without a climate fixed effect ---------------------------------
cat("\n[1] inert for the deployed (climate-free) term set\n")
expect(isTRUE(.assert_no_zero_filled_climate(deployed, panel, window = NULL)),
       "no climate FE: passes over the whole zero-filled panel")
expect(isTRUE(.assert_no_zero_filled_climate(c("mso_band2", "log_pop_c", "shape_Intercept"),
                                             panel[, c("gid", "t_start")])),
       "no climate FE: passes even when the panel has no climate columns")

# --- [2] refuses a climate FE over zero-filled rows ---------------------------
cat("\n[2] refuses a climate FE over zero-filled rows\n")
for (fe in c("mean_precip", "dev_precip", "mean_temp_daymean", "dev_temp_daymean",
             "mean_precip_lag2", "dev_temp_daymean_cb3")) {
  expect(refuses(.assert_no_zero_filled_climate(c(deployed, fe), panel, window = NULL)),
         sprintf("climate FE '%s' + zero-filled rows in the window: refused", fe))
}
expect(refuses(.assert_no_zero_filled_climate(
         with_clim, panel, window = as.Date(c("2023-10-01", "2024-03-01")))),
       "window straddling the DB end: refused")
expect(refuses(.assert_no_zero_filled_climate(
         with_clim, panel, window = as.Date(c("2024-01-01", NA)))),
       "open-ended window past the DB end: refused")
msg <- tryCatch(.assert_no_zero_filled_climate(with_clim, panel, label = "unit-test"),
                error = function(e) conditionMessage(e))
expect(grepl("unit-test", msg) && grepl("mean_temp_daymean", msg) &&
         grepl("2024-01-01 to 2024-06-01", msg, fixed = TRUE),
       "error names the label, the climate term and the zero-filled months")

# --- [3] passes when the window holds only observed climate -------------------
cat("\n[3] passes on observed-climate windows\n")
expect(isTRUE(.assert_no_zero_filled_climate(
         with_clim, panel, window = as.Date(c("2023-07-01", "2023-12-01")))),
       "window ending at the DB end: passes")
expect(isTRUE(.assert_no_zero_filled_climate(
         with_clim, panel, window = as.Date(c(NA, "2023-11-01")))),
       "open-start window before the DB end: passes")

# A real month with SOME climate columns at 0 (e.g. dry month, zero anomaly) is
# not a filled row: all four must be exactly 0.
p2 <- panel[!post, ]
p2$mean_precip[1] <- 0; p2$dev_precip[1] <- 0
expect(isTRUE(.assert_no_zero_filled_climate(with_clim, p2)),
       "a month with only some climate columns at 0 is not a filled row")
p3 <- panel[!post, ]
p3[1, c("mean_precip", "dev_precip", "mean_temp_daymean", "dev_temp_daymean")] <- NA
expect(isTRUE(.assert_no_zero_filled_climate(with_clim, p3)),
       "NA climate (not zero) is left to the other gates")

# --- [4] a climate FE with no climate columns in the panel --------------------
cat("\n[4] climate FE but no climate columns\n")
expect(refuses(.assert_no_zero_filled_climate(with_clim, panel[, c("gid", "t_start")])),
       "climate FE + panel without the climate columns: refused")

# --- [5] .draws_fe_names() over single bundles and *_by_group lists -----------
cat("\n[5] .draws_fe_names()\n")
by_group <- list(enso_east = fe_bundle(deployed), fe_group = fe_bundle(c("Intercept", "mean_precip")))
got <- .draws_fe_names(fe_bundle(c("Intercept", "mso_band2")), NULL, by_group, NULL)
expect(setequal(got, c("Intercept", "mso_band2", deployed, "mean_precip")),
       "union over single bundles, NULLs and group lists (incl. a group named 'fe...')")
expect(identical(.draws_fe_names(NULL, NULL), character(0)),
       "all-NULL input gives character(0)")
expect(refuses(.assert_no_zero_filled_climate(.draws_fe_names(by_group), panel)),
       "a climate FE in ONE group of a grouped run is enough to refuse")
expect(isTRUE(.assert_no_zero_filled_climate(
         .draws_fe_names(list(enso_east = fe_bundle(deployed),
                              enso_west = fe_bundle(deployed))), panel)),
       "grouped run with the deployed term set in every group: passes")

# --- [6] .sim_panel_window() --------------------------------------------------
cat("\n[6] .sim_panel_window()\n")
sim <- list(horizon = list(start = "2026-06-01", end = "2027-03-01"))
w <- .sim_panel_window(sim)
expect(identical(w, as.Date(c("2026-06-01", "2027-03-01"))),
       "window = horizon start .. horizon end")
w <- .sim_panel_window(sim, enso_sim_months_start = as.Date("2026-04-01"))
expect(identical(w[1], as.Date("2026-04-01")),
       "an earlier smoothing handoff extends the window back")
sim_wu <- c(sim, list(warmup = list(use_sim_trajectories = TRUE, warmup_months = 6L)))
expect(identical(.sim_panel_window(sim_wu)[1], as.Date("2025-12-01")),
       "hybrid warmup extends the window back by warmup_months")
w <- .sim_panel_window(list(horizon = list(start = "2026-06-01")))
expect(is.na(w[2]) && identical(w[1], as.Date("2026-06-01")),
       "no horizon end: open-ended window")

# --- [7] every panel loader calls the guard ----------------------------------
cat("\n[7] call sites\n")
n_calls <- function(path) {
  src <- readLines(path, warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  sum(grepl(".assert_no_zero_filled_climate(", src, fixed = TRUE))
}
for (cs in list(c("analysis/utils/simulations.R", 2L),      # ungrouped + grouped loader
                c("analysis/utils/smoothing_core.R", 1L),   # smoothing_load_artefacts()
                c("analysis/07b_simulator_v1.R", 1L))) {    # v1's own ungrouped loader
  expect(n_calls(cs[1]) == as.integer(cs[2]),
         sprintf("%s calls the guard %s time(s)", cs[1], cs[2]))
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
