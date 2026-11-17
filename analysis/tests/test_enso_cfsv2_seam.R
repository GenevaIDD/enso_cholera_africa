#!/usr/bin/env Rscript
# test_enso_cfsv2_seam.R — HARD: sample_enso_cfsv2_trajectories() when the ONI
#   record ends after the raw detrended Niño 3.4 record.
#
# Usage:
#   Rscript analysis/tests/test_enso_cfsv2_seam.R
#
# Details:
#   That is CPC's normal state (e.g. Niño 3.4 to Jun, ONI to Aug). A seam at
#   the ONI end would leave NA raw values in between, an NA 3-month rolling
#   mean around them, and the first month past the ONI end on the 0.0
#   fallback: a one-month neutral notch in the forecast ENSO input, in the
#   manuscript ONI envelope AND in the simulator's forecast arm. Synthetic
#   data only.

suppressPackageStartupMessages(library(lubridate))
source("analysis/utils/base.R") # enso_lag_col_name() (call-time dependency of the samplers)
source("analysis/utils/enso.R")

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

# --- synthetic records -------------------------------------------------------
# ONI to T (Aug-2026), Niño 3.4 to T-2 (Jun-2026), CFSv2 init T-10 with targets
# T-10 .. T+10, sim months T-19 .. T+20 (so the horizon outruns CFSv2 by 10).
T_oni <- as.Date("2026-08-01")
T_n34 <- T_oni %m-% months(2L)
months_ <- function(a, b) seq.Date(a, b, by = "1 month")
ym <- function(d) {
  list(year = as.integer(format(d, "%Y")), month = as.integer(format(d, "%m")))
}

oni_dates <- months_(as.Date("2020-01-01"), T_oni)
oni_vals <- 2 * sin(seq_along(oni_dates) / 6) # smooth, non-zero, |x| up to 2
historical_oni <- data.frame(
  year = ym(oni_dates)$year,
  month = ym(oni_dates)$month,
  enso_index = oni_vals
)
n34_dates <- months_(as.Date("2020-01-01"), T_n34)
historical_n34 <- data.frame(
  year = ym(n34_dates)$year,
  month = ym(n34_dates)$month,
  enso_index = oni_vals[seq_along(n34_dates)] + 0.1
)

cfs_targets <- months_(T_oni %m-% months(10L), T_oni %m+% months(10L))
cfsv2_data <- expand.grid(target_date = cfs_targets, member = 1:5)
cfsv2_data$init_date <- T_oni %m-% months(10L)
cfsv2_data$lead <- as.integer(
  interval(cfsv2_data$init_date, cfsv2_data$target_date) %/% months(1)
)
# Strong, non-zero ensemble: 2.3 +/- small member offset
cfsv2_data$nino34_anomaly <- 2.3 + 0.05 * (cfsv2_data$member - 3)

sim_months <- months_(T_oni %m-% months(19L), T_oni %m+% months(20L))

res <- sample_enso_cfsv2_trajectories(
  cfsv2_data = cfsv2_data,
  n_rep = 20L,
  historical_oni = historical_oni,
  historical_nino34 = historical_n34,
  sim_months = sim_months,
  max_lag = 12L,
  rebase = "seam",
  sd_correction = FALSE,
  seed = 1L
)
lag0 <- res$enso_lag0
rownames(lag0) <- as.character(sim_months)
at <- function(d) lag0[as.character(d), ]

cat("[test_enso_cfsv2_seam] lag0 around the seam:\n")
for (d in as.list(months_(T_n34 %m-% months(1L), T_oni %m+% months(3L)))) {
  cat(sprintf(
    "   %s  mean=%7.3f  min=%7.3f\n",
    format(d),
    mean(at(d)),
    min(at(d))
  ))
}

# 1. The month after the ONI end is not a neutral notch: it carries the
#    seam-rebased CFSv2 ensemble (2.3 rebased onto the Jun-2026 raw value),
#    i.e. a value of the same order as its neighbours, not 0.
nxt <- at(T_oni %m+% months(1L))
ok("month after ONI end is not 0.0", all(nxt != 0) && all(abs(nxt) > 0.5))
# 2. No exact zero anywhere inside the CFSv2 span.
in_span <- sim_months[sim_months <= max(cfs_targets) %m-% months(1L)]
ok(
  "no exact 0.0 inside the CFSv2 span",
  !any(lag0[as.character(in_span), ] == 0)
)
# 3. Months with published ONI reproduce it exactly, including the two months
#    whose raw Niño 3.4 is missing (T-1, T).
for (d in as.list(months_(T_n34 %m-% months(2L), T_oni))) {
  ok(
    sprintf("published ONI preferred at %s", format(d)),
    isTRUE(all.equal(unname(at(d)), rep(oni_vals[match(d, oni_dates)], 20L)))
  )
}
# 4. Replicates differ past the ONI end (CFSv2 members are sampled, not a scalar).
ok(
  "replicates spread past the ONI end",
  length(unique(at(T_oni %m+% months(2L)))) > 1
)
# 5. Beyond the CFSv2 span the documented neutral fallback still applies
#    (first month whose whole 3-month window is past the last CFSv2 target).
beyond <- at(max(cfs_targets) %m+% months(2L))
cat(sprintf(
  "   beyond-span month %s: range %.3f..%.3f\n",
  format(max(cfs_targets) %m+% months(2L)),
  min(beyond),
  max(beyond)
))
ok("neutral fallback beyond the CFSv2 span", all(abs(beyond) < 1e-9)) # rollmean is cumsum-based: ~1e-17 residue
# 6. Lagged views are consistent with lag0 (lag 3 at m == lag0 at m-3).
m <- T_oni %m+% months(4L)
ok(
  "enso_lag3 is lag0 shifted by three months",
  isTRUE(all.equal(
    res$enso_lag3[match(m, sim_months), ],
    lag0[as.character(m %m-% months(3L)), ]
  ))
)

# --- numeric rebase given as a string (the CLI's --cfsv2_rebase 0.3) ----------
# A string offset must not be skipped by an is.numeric() check (delta 0, silently).
samp <- function(rb) {
  sample_enso_cfsv2_trajectories(
    cfsv2_data = cfsv2_data,
    n_rep = 20L,
    historical_oni = historical_oni,
    historical_nino34 = historical_n34,
    sim_months = sim_months,
    max_lag = 12L,
    rebase = rb,
    sd_correction = FALSE,
    seed = 1L
  )$enso_lag0
}
r_num <- samp(0.3)
r_str <- samp("0.3")
r_none <- samp("none")
r_zero <- samp(0)
ok("rebase '0.3' (string) == rebase 0.3 (numeric)", identical(r_str, r_num))
ok(
  "rebase 0.3 shifts the forecast (differs from 'none')",
  !identical(r_num, r_none)
)
ok("rebase 0 == 'none'", identical(r_zero, r_none))
ok(
  ".rebase_offset(): numbers and numeric strings, NULL otherwise",
  identical(.rebase_offset("0.3"), 0.3) &&
    identical(.rebase_offset(-1L), -1) &&
    is.null(.rebase_offset("seam")) &&
    is.null(.rebase_offset("none"))
)
ok(
  "cfsv2_correction_provenance(): default sd + seam adds nothing",
  identical(cfsv2_correction_provenance("sd", "seam"), list()) &&
    identical(cfsv2_correction_provenance(NA_character_, "seam"), list())
)
ok(
  "cfsv2_correction_provenance(): calibration records the effective rebase 'none'",
  identical(
    cfsv2_correction_provenance("calibration", "seam"),
    list(cfsv2_correction = "calibration", cfsv2_rebase_effective = "none")
  )
)
ok(
  "cfsv2_correction_provenance(): a string offset is recorded as a number",
  identical(
    cfsv2_correction_provenance("sd", "0.3")$cfsv2_rebase_effective,
    0.3
  )
)

cat(sprintf("\n[test_enso_cfsv2_seam] %d passed, %d failed\n", pass, fail))
if (fail > 0L) {
  quit(status = 1L)
}
