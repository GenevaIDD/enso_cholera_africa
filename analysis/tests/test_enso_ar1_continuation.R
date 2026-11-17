# test_enso_ar1_continuation.R — HARD: neutral and AR(1) continuation of ENSO
#   past the CFSv2 span.
#
# Usage:
#   Rscript analysis/tests/test_enso_ar1_continuation.R
#
# Details:
#   Checks legacy_seam [1], neutral (default) [2], ar1 [3] and the AR(1)
#   helpers [4]. Synthetic inputs only.

suppressPackageStartupMessages({
  library(lubridate)
  library(dplyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/enso.R")
err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

# Previous sampler, for the legacy comparison: the newest commit whose enso.R
# has no `legacy_seam` argument.
old_env <- new.env()
old_src <- NULL
for (rev in c(
  "HEAD",
  "HEAD~1",
  "HEAD~2",
  "HEAD~3",
  "HEAD~4",
  "HEAD~5",
  "HEAD~6"
)) {
  txt <- tryCatch(
    system2(
      "git",
      c("show", paste0(rev, ":analysis/utils/enso.R")),
      stdout = TRUE,
      stderr = FALSE
    ),
    warning = function(w) NULL,
    error = function(e) NULL
  )
  if (length(txt) && !any(grepl("legacy_seam", txt, fixed = TRUE))) {
    old_src <- txt
    break
  }
}
have_old <- !is.null(old_src)
if (have_old) {
  tf <- tempfile(fileext = ".R")
  writeLines(old_src, tf)
  sys.source(tf, envir = old_env)
  legacy_names(old_env)
  environment(old_env$sample_enso_cfsv2_trajectories) <- old_env
}

# --- synthetic inputs ----------------------------------------------------------
hist_d <- seq(as.Date("1985-01-01"), as.Date("2026-06-01"), by = "month")
nino <- 1.2 * sin(seq_along(hist_d) / 7)
hist_n <- tibble::tibble(
  year = year(hist_d),
  month = month(hist_d),
  enso_index = nino
)
oni <- as.numeric(stats::filter(nino, rep(1 / 3, 3), sides = 2))
hist_o <- tibble::tibble(
  year = year(hist_d),
  month = month(hist_d),
  enso_index = oni
) %>%
  dplyr::filter(
    !is.na(enso_index),
    hist_d[seq_along(oni)] <= as.Date("2026-05-01")
  )
cfs_d <- seq(as.Date("2026-05-01"), as.Date("2026-12-01"), by = "month")
cfs <- tidyr::expand_grid(target_date = cfs_d, member = 1:3) %>%
  dplyr::mutate(
    nino34_anomaly = 1.5 + 0.1 * member + 0.05 * as.numeric(factor(target_date))
  )
sim_m <- seq(as.Date("2026-07-01"), as.Date("2027-06-01"), by = "month")
run <- function(fun = sample_enso_cfsv2_trajectories, ...) {
  suppressMessages(fun(
    cfs,
    n_rep = 6L,
    historical_oni = hist_o,
    historical_nino34 = hist_n,
    sim_months = sim_m,
    max_lag = 12L,
    rebase = "none",
    seed = 7L,
    ...
  ))
}

# legacy_seam = TRUE reproduces the sampler as of the parent commit
# (git show HEAD~ ... enso.R, sourced into its own environment) exactly.
cat("\n[1] legacy_seam reproduces the previous sampler\n")
leg <- run(legacy_seam = TRUE)
if (have_old) {
  old <- run(old_env$sample_enso_cfsv2_trajectories)
  expect(
    identical(leg, old),
    "legacy_seam = TRUE == previous sampler, every lag matrix"
  )
} else {
  cat("  SKIP  previous enso.R not found in git history\n")
}

# ONI past the last CFSv2 month is 0; the last span month's ONI averages its
# two available months; everything earlier is identical to legacy; the horizon
# month of a lag-L track reads the unbiased last-span value.
cat("\n[2] neutral (default)\n")
neu <- run()
last <- which(sim_m == as.Date("2026-12-01"))
after <- which(sim_m > as.Date("2026-12-01"))
expect(all(neu$enso_lag0[after, ] == 0), "lag0 past the span = 0")
expect(
  all(leg$enso_lag0[after[1], ] != 0),
  "legacy lag0 one month past the span was not 0 (mixed mean)"
)
cfs_w <- tidyr::pivot_wider(
  cfs,
  names_from = member,
  values_from = nino34_anomaly
)
# member drawn per replicate is the same in both runs (same seed, same draw)
mem_of <- match(
  round(leg$enso_lag0[last - 1L, ], 8),
  round(leg$enso_lag0[last - 1L, ], 8)
)
expect(
  all(abs(neu$enso_lag0[last, ] - leg$enso_lag0[last, ]) > 0),
  "last span month differs from legacy"
)
expect(
  identical(
    neu$enso_lag0[seq_len(last - 1L), ],
    leg$enso_lag0[seq_len(last - 1L), ]
  ),
  "months before the last span month identical to legacy"
)
# lag-L track at the month L after the span end reads the last-span ONI
expect(
  identical(neu$enso_lag3[last + 3L, ], neu$enso_lag0[last, ]),
  "lag3 at span end + 3 reads the last-span ONI"
)

# The first continued raw month equals phi x the member's last CFSv2 value; the
# continuation decays monotonically to 0; deterministic under the seed; months
# before the span end unchanged.
cat("\n[3] ar1\n")
phi <- 0.8
ar <- run(beyond_cfsv2 = "ar1", ar1_phi = phi)
expect(
  identical(
    ar$enso_lag0[seq_len(last - 1L), ],
    neu$enso_lag0[seq_len(last - 1L), ]
  ),
  "months before the last span month identical to neutral"
)
v <- ar$enso_lag0[after, , drop = FALSE]
expect(
  all(apply(v, 2, function(x) all(diff(x) < 0) && all(x > 0))),
  "continuation decays monotonically towards 0"
)
ar2 <- run(beyond_cfsv2 = "ar1", ar1_phi = phi)
expect(identical(ar, ar2), "deterministic under the seed")
# raw continuity: first continued raw month = phi * last raw value; recover from ONI
# ONI(after1) = (raw[last] + phi raw[last] + phi^2 raw[last]) / 3
raw_last <- 3 * ar$enso_lag0[after[1], ] / (1 + phi + phi^2)
expect(
  isTRUE(all.equal(
    ar$enso_lag0[after[2], ],
    raw_last * (phi + phi^2 + phi^3) / 3
  )),
  "ONI past the span follows phi^k x the member's last raw value"
)
raw_prev <- 2 * neu$enso_lag0[last, ] - raw_last # neutral: (prev + last) / 2
expect(
  isTRUE(all.equal(
    ar$enso_lag0[last, ],
    (raw_prev + raw_last + phi * raw_last) / 3
  )),
  "last span month: centred mean over (prev, last, phi x last)"
)
expect(
  !is.null(err_of(run(beyond_cfsv2 = "lenssen"))),
  "unknown continuation method refused"
)

# fit_oni_ar1() recovers phi from a simulated AR(1) series; continue_enso_ar1()
# shape and guard.
cat("\n[4] fit_oni_ar1 / continue_enso_ar1\n")
set.seed(1)
x <- numeric(600)
for (i in 2:600) {
  x[i] <- 0.9 * x[i - 1] + rnorm(1, sd = 0.2)
}
xs <- stats::setNames(
  x,
  format(seq(as.Date("1971-01-01"), by = "month", length.out = 600))
)
f <- fit_oni_ar1(xs, "1971_2020")
expect(
  abs(f$phi - 0.9) < 0.05 && f$n == 600L,
  sprintf("phi recovered (%.3f)", f$phi)
)
expect(!is.null(err_of(fit_oni_ar1(xs, "2030_2040"))), "empty window refused")
m <- continue_enso_ar1(c(1, -2), 3, 0.5)
expect(
  identical(dim(m), c(3L, 2L)) && isTRUE(all.equal(m[, 2], c(-1, -0.5, -0.25))),
  "shape and values"
)
expect(!is.null(err_of(continue_enso_ar1(1, 2, 1.0))), "|phi| >= 1 refused")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
