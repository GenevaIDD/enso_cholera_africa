#!/usr/bin/env Rscript
# analysis/tests/test_enso_cfsv2_calibration.R
#
# CFSv2 hindcast calibration (CPC "PDF + spread correction" replica, utils/enso.R):
# fitCFSv2Calibration(), applyCFSv2Calibration(), .cfsv2AnchorAnalysisRows()
# and the calibration mode of sampleEnsoCFSv2Trajectories(); plus the
# crossbasis boundary-knot clamp helper .clamp_cb_range() (utils/dlnm.R) used by
# the forecast-phase simulator (simulator.enso.clamp_to_knots).
#
# Synthetic data only, no MCMC. Registered as a HARD gate in analysis/tests/run_all.R.

suppressPackageStartupMessages({ library(dplyr); library(lubridate) })
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
source("analysis/utils/dlnm.R")

pass <- 0L; fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("  [PASS]", label, "\n") }
  else              { fail <<- fail + 1L; cat("  [FAIL]", label, "\n") }
}

# --- synthetic hindcast: model = 1.5 x observed + 0.4 + member noise ---------
# Observed anomaly y ~ N(0, 1) per (start, lead) target; the model is too
# strong (x1.5) and too warm (+0.4) in the 26 degC climatology, so a correct
# percentile map must undo both, and the spread factor must match the error.
set.seed(1)
starts <- seq(as.Date("1991-01-01"), as.Date("2020-12-01"), by = "month")
leads  <- 0:2
n_mem  <- 20L
targets <- sort(unique(do.call(c, lapply(leads, function(l) starts %m+% months(l)))))
obs_tbl <- tibble(date = targets, enso_index = rnorm(length(targets)))
hind <- bind_rows(lapply(leads, function(l) {
  td <- starts %m+% months(l)
  y  <- obs_tbl$enso_index[match(td, obs_tbl$date)]
  bind_rows(lapply(seq_len(n_mem), function(m)
    tibble(init_date = starts, lead = l, target_date = td, member = m,
           sst = 26 + 0.4 + 1.5 * y + rnorm(length(y), sd = 0.3))))
}))

cat("\n## fitCFSv2Calibration()\n")
cal <- fitCFSv2Calibration(hind, list(ersst = obs_tbl), window = c(1991L, 2020L))
ok("fits 12 start months x 3 leads", nrow(cal$clim) == 36L)
ok("reference samples per cell: 30 obs years, 600 member values",
   all(cal$targets$ersst$ref$n_years == 30L) &&
   all(lengths(cal$targets$ersst$ref$model) == 30L * n_mem))
ok("model climatology ~ 26.4 (warm bias kept in clim)", abs(mean(cal$clim$clim) - 26.4) < 0.1)

cat("\n## .cfsv2QuantileMap()\n")
r <- cal$targets$ersst$ref[1, ]
x_mid <- stats::quantile(r$model[[1]], c(0.25, 0.5, 0.75), names = FALSE)
y_mid <- stats::quantile(r$obs[[1]],   c(0.25, 0.5, 0.75), names = FALSE)
ok("quartiles map onto observed quartiles (within 0.15)",
   max(abs(.cfsv2QuantileMap(x_mid, r$model[[1]], r$obs[[1]]) - y_mid)) < 0.15)
cap <- .cfsv2QuantileMap(c(-99, 99), r$model[[1]], r$obs[[1]])
ok("values beyond the reference range are held at the reference extremes",
   identical(cap, .cfsv2QuantileMap(range(r$model[[1]]), r$model[[1]], r$obs[[1]])) &&
   cap[1] >= min(r$obs[[1]]) && cap[2] <= max(r$obs[[1]]))
ok("map is monotone", !is.unsorted(.cfsv2QuantileMap(seq(-5, 5, 0.1), r$model[[1]], r$obs[[1]])))

cat("\n## applyCFSv2Calibration()\n")
# Real-time matrix: 2 analysis rows (identical members) then 3 forecast rows.
dates <- as.Date(c("2026-08-01", "2026-09-01", "2026-10-01", "2026-11-01", "2026-12-01"))
X <- rbind(rep(1.0, 8), rep(1.2, 8),
           1.5 * c(-1, -.5, -.2, 0, .1, .3, .6, 1) + 0.5,
           1.5 * c(-1, -.5, -.2, 0, .1, .3, .6, 1) + 0.8,
           1.5 * c(-1, -.5, -.2, 0, .1, .3, .6, 1) + 0.2)
rownames(X) <- as.character(dates)
Y <- applyCFSv2Calibration(X, cal, "ersst")
ok("analysis rows left unchanged", identical(Y[1:2, ], X[1:2, ]))
ok("first forecast row detected (row 3)", identical(attr(Y, "first_forecast_row"), 3L))
g1 <- cal$targets$ersst$spread$g[cal$targets$ersst$spread$sm == 10 & cal$targets$ersst$spread$lead == 0]
ok("spread = g x raw spread (lead 0)", abs(sd(Y[3, ]) - g1 * sd(X[3, ])) < 1e-12)
ok("member order preserved (coherent replicates)", identical(order(Y[4, ]), order(X[4, ])))
ok("too-strong model anomaly is shrunk toward obs scale",
   abs(mean(Y[4, ])) < abs(mean(X[4, ])))
ok("unknown target stops", inherits(try(applyCFSv2Calibration(X, cal, "oisst"), silent = TRUE), "try-error"))

cat("\n## .cfsv2AnchorAnalysisRows()\n")
# ONI published to Sep (centre) i.e. ASO not yet; JAS centre = Aug.
oni_lkp  <- stats::setNames(c(0.9, 1.1), c("2026-07-01", "2026-08-01"))
hist_lkp <- stats::setNames(c(0.8, 0.9), c("2026-06-01", "2026-07-01"))
Xa <- X; Xa[1, ] <- 1.6; Xa[2, ] <- 1.9          # Aug, Sep analysis rows (unobserved)
an <- .cfsv2AnchorAnalysisRows(Xa, 3L, hist_lkp, oni_lkp)
# Latest centre with a full window: Aug (Jul obs 0.9, Aug/Sep analysis):
# shift = (3*1.1 - (0.9 + 1.6 + 1.9)) / 2 = -0.55
ok("anchors to latest ONI season (centre Aug)", identical(an$season_centre, as.Date("2026-08-01")))
ok("shift reproduces the ONI 3-month mean", abs((0.9 + mean(an$cfs_wide[1, ]) + mean(an$cfs_wide[2, ])) / 3 - 1.1) < 1e-12)
ok("forecast rows untouched", identical(an$cfs_wide[3:5, ], Xa[3:5, ]))
an0 <- .cfsv2AnchorAnalysisRows(Xa, 3L, stats::setNames(c(0.8, 0.9, 1.0, 1.2), c("2026-06-01", "2026-07-01", "2026-08-01", "2026-09-01")), oni_lkp)
ok("no-op when the observed series covers the analysis rows", an0$delta == 0 && identical(an0$cfs_wide, Xa))

cat("\n## sampleEnsoCFSv2Trajectories(calibration = ...)\n")
cf_long <- expand.grid(i = seq_len(nrow(X)), member = seq_len(ncol(X)))
cfsv2_data <- tibble(init_date = dates[1], target_date = dates[cf_long$i], lead = cf_long$i - 1L,
                     member = cf_long$member, nino34_anomaly = X[cbind(cf_long$i, cf_long$member)])
oni_hist <- tibble(year = c(2026L, 2026L), month = c(7L, 8L), enso_index = c(0.9, 1.1))
n34_hist <- tibble(year = c(2026L, 2026L), month = c(6L, 7L), enso_index = c(0.8, 0.9))
sm <- seq(as.Date("2026-09-01"), as.Date("2026-11-01"), by = "month")
res <- suppressMessages(sampleEnsoCFSv2Trajectories(cfsv2_data, n_rep = 50, historical_oni = oni_hist,
         historical_nino34 = n34_hist, sim_months = sm, max_lag = 1L, seed = 1,
         calibration = cal, calibration_target = "ersst"))
ok("returns enso_lag0/1 matrices [M x n_rep]", all(dim(res$enso_lag0) == c(3L, 50L)))
ok("no NA in forecast ENSO", !anyNA(res$enso_lag0))
# Pin: under calibration the seam rebase is skipped -- the
# sampler sets rebase <- "none" after anchoring the analysis rows
# (sampleEnsoCFSv2Trajectories(), utils/enso.R) -- so simulator.enso.cfsv2_rebase
# is inert there: "seam" (the default), "none" and a numeric shift give the same
# trajectories. A change that makes the rebase act under calibration must
# update this check deliberately.
run_rb <- function(rb, cal_arg = cal)
  suppressMessages(sampleEnsoCFSv2Trajectories(cfsv2_data, n_rep = 50, historical_oni = oni_hist,
    historical_nino34 = n34_hist, sim_months = sm, max_lag = 1L, seed = 1, rebase = rb,
    calibration = cal_arg, calibration_target = "ersst"))
res_rb <- lapply(list(seam = "seam", none = "none", shift = 0.3), run_rb)
ok("calibration on: the default call is rebase = 'seam'", identical(res_rb$seam, res))
ok("calibration on: identical output under rebase 'seam', 'none' and 0.3",
   identical(res_rb$seam, res_rb$none) && identical(res_rb$seam, res_rb$shift))
# Control: without calibration the same numeric shift does move the output,
# so the identity above is the skip, not an inert argument.
ok("calibration off: rebase 0.3 differs from 'none' (control)",
   !identical(run_rb(0.3, NULL)$enso_lag0, run_rb("none", NULL)$enso_lag0))
for (args in list(list(sd_correction = TRUE), list(target_peak = 2.5)))
  ok(paste("mutually exclusive with", names(args)),
     inherits(try(do.call(sampleEnsoCFSv2Trajectories,
                          c(list(cfsv2_data, n_rep = 5, historical_oni = oni_hist, calibration = cal), args)),
                  silent = TRUE), "try-error"))

cat("\n## .clamp_cb_range() (simulator.enso.clamp_to_knots)\n")
ok("no-op inside the range", identical(.clamp_cb_range(c(-1, 0, 2.5), c(-1.6, 2.59)), c(-1, 0, 2.5)))
ok("held exactly at the edges", identical(.clamp_cb_range(c(-3, 3), c(-1.6, 2.59)), c(-1.6, 2.59)))
ok("NULL range is a no-op", identical(.clamp_cb_range(c(-3, 3), NULL), c(-3, 3)))
ok("margin widens each bound by its own magnitude (25%: -1.6 -> -2, 2.4 -> 3)",
   isTRUE(all.equal(.clamp_cb_range(c(-9, -1.8, 2.9, 9), c(-1.6, 2.4), 0.25), c(-2, -1.8, 2.9, 3))))
ok("matrix shape preserved", identical(dim(.clamp_cb_range(matrix(c(-3, 0, 3, 1), 2), c(-1, 2))), c(2L, 2L)))

cat(sprintf("\n[test_enso_cfsv2_calibration] %d passed, %d failed\n", pass, fail))
if (fail > 0L) quit(status = 1L)
