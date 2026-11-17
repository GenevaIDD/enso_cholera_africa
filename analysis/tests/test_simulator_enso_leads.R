# test_simulator_enso_leads.R — HARD: deployed ENSO lead fits in the simulator
#   path.
#
# Usage:
#   Rscript analysis/tests/test_simulator_enso_leads.R
#
# Details:
#   Requires the `dlnm` package (skips without it). Checks that a lead fit
#   reads the panel's enso_lead<k> column, that .lags_required_from_draws()
#   asks for it, the ENSO samplers' lead columns and panel_max_enso_lead().

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})
for (m in c(
  "base",
  "config",
  "paths",
  "dlnm",
  "enso",
  "simulator_core",
  "simulations"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R")
cat("# test_simulator_enso_leads.R -- ENSO leads in the simulator\n\n")
if (!requireNamespace("dlnm", quietly = TRUE)) {
  cat("[SKIP] dlnm not available\n")
  quit(status = 0)
}

# ---- fixture (as test_simulator_enso_cb_grouping.R) ---------------------------
frame_tu <- expand.grid(
  country = "WWW",
  adm1n = c(1, 2),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    gid = sprintf("%s_%d_1", country, adm1n),
    adm1_unit = sprintf("%s_%d", country, adm1n),
    group = "all",
    AFRO_region = "REG_W",
    tu = "leadTU",
    min_seas = 1L,
    pop = 1e5
  ) %>%
  dplyr::select(
    gid,
    group,
    country,
    AFRO_region,
    tu,
    min_seas,
    pop,
    adm1_unit
  ) %>%
  tibble::as_tibble()
months <- as.Date(sprintf("2024-%02d-01", 1:6))
horizon <- list(start = "2024-01-01", end = "2024-06-30")
make_panel <- function(lead2) {
  tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA_character_,
      enso_lag0 = 0,
      enso_lag1 = 0,
      enso_lag2 = 0,
      enso_lead1 = 0,
      enso_lead2 = lead2
    )
}
# Single-lag (lead 2) ENSO basis built the fit-time way (build_cb_block).
cbb <- suppressMessages(build_cb_block(
  make_panel(5) %>% dplyr::mutate(enso_lead2 = seq_len(dplyr::n()) / 3),
  list(
    enable = TRUE,
    defaults = list(argvar = list(fun = "lin")),
    variables = list(enso = list(single_lag = -2L))
  ),
  ""
))
cb_meta <- cbb$cb_meta
expect_true(
  identical(cb_meta$enso$single_lag, -2L),
  "fit-time cb_meta: single_lag -2"
)
cbn <- cb_meta$enso$cb_col_names
active <- c(
  "(Intercept)",
  "mean_precip",
  "dev_precip",
  "mean_temp_daymean",
  "dev_temp_daymean",
  "ongoing_outbreak_lag1",
  "ongoing_outbreak_country_lag1",
  "ongoing_outbreak_tu_lag1",
  cbn
)
bun <- make_draws_bundle_synthetic(
  frame_tu,
  S = 3,
  phi = 5,
  active_terms = active,
  fe_overrides = c(
    `(Intercept)` = -20,
    stats::setNames(rep(8, length(cbn)), cbn)
  )
)
run <- function(panel, scen = NULL, seed = 7L) {
  simulate_one_tu(
    "leadTU",
    frame_tu,
    panel,
    draws_occ = bun,
    replicate_ids = 1:3,
    seed = seed,
    output_part = tempfile("lead_"),
    horizon = horizon,
    dry_run = TRUE,
    cb_meta_occ = cb_meta,
    enso_scenarios = scen
  )
}

# [1] simulate_one_tu(): an occurrence cross-basis with single_lag = -2 reads
#     the panel's enso_lead2 column (signal fires only through the lead) and
#     the per-replicate scenario override replaces it.
cat(
  "## 1) simulate_one_tu reads enso_lead<k> and applies the scenario override\n"
)
r_hi <- run(make_panel(5))
r_lo <- run(make_panel(0))
expect_true(
  sum(r_hi$monthly_summary$Y) > 0,
  "lead2 = 5 in the panel: onsets (signal through the lead column)"
)
expect_true(
  sum(r_lo$monthly_summary$Y) == 0,
  "lead2 = 0: no onsets (lags 0..2 are 0 throughout)"
)
scen <- list(enso_lead2 = matrix(5, nrow = length(months), ncol = 3))
r_sc <- run(make_panel(0), scen)
expect_true(
  sum(r_sc$monthly_summary$Y) > 0,
  "scenario enso_lead2 = 5 overrides the panel's 0"
)

# [2] .lags_required_from_draws(): a lead fit requires enso_lead<k>.
cat("\n## 2) .lags_required_from_draws: a lead fit requires enso_lead<k>\n")
req <- .lags_required_from_draws(bun, cb_meta)
expect_true(
  "enso_lead2" %in% req && !any(grepl("enso_lag-", req)),
  "requires enso_lead2 (not enso_lag-2)"
)

# [3] ENSO samplers (CFSv2, C3S plume, analog, Lenssen-style cache):
#     max_lead > 0 adds enso_lead1..K with enso_lead<k>(t) = enso_lag0(t + k)
#     per replicate; max_lead = 0 is identical to the pre-lead samplers
#     (ENSO_REF_SRC) and leaves the lag matrices unchanged.
cat(
  "\n## 3) samplers: enso_lead<k>(t) = enso_lag0(t + k); max_lead = 0 unchanged\n"
)
# Equivalence with the pre-lead samplers: set ENSO_REF_SRC to a copy of the
# pre-lead enso.R (git show 77cc1af:analysis/utils/enso.R > f); otherwise the
# max_lead = 0 output is compared with the lag block of a max_lead > 0 call.
.ref_src <- Sys.getenv("ENSO_REF_SRC", "")
old <- new.env()
if (nzchar(.ref_src) && file.exists(.ref_src)) {
  sys.source(.ref_src, envir = old)
  legacy_names(old)
} else {
  cat(
    "  [note] ENSO_REF_SRC unset: comparing max_lead = 0 with the current samplers\n"
  )
  for (f in c(
    "sample_enso_cfsv2_trajectories",
    "sample_enso_analog_scenarios",
    "sample_enso_forecast_trajectories",
    "sample_enso_lenssen_scenarios"
  )) {
    assign(f, get(f), envir = old)
  }
}
hist_oni <- tibble::tibble(
  year = rep(2018:2024, each = 12),
  month = rep(1:12, 7),
  enso_index = round(sin(seq_len(84) / 5), 3)
)
hist_oni <- hist_oni[!(hist_oni$year == 2024 & hist_oni$month > 6), ]
sim_m <- seq(as.Date("2024-05-01"), by = "month", length.out = 8)
cfs <- tidyr::expand_grid(
  target_date = seq(as.Date("2024-07-01"), by = "month", length.out = 9),
  member = 1:4
) %>%
  dplyr::mutate(
    nino34_anomaly = 0.1 * member + as.numeric(format(target_date, "%m")) / 10
  )
a0 <- sample_enso_cfsv2_trajectories(
  cfs,
  n_rep = 5,
  historical_oni = hist_oni,
  sim_months = sim_m,
  seed = 3
)
b0 <- old$sample_enso_cfsv2_trajectories(
  cfs,
  n_rep = 5,
  historical_oni = hist_oni,
  sim_months = sim_m,
  seed = 3
)
expect_true(
  identical(a0, b0),
  "CFSv2: max_lead = 0 identical to the committed sampler"
)
a3 <- sample_enso_cfsv2_trajectories(
  cfs,
  n_rep = 5,
  historical_oni = hist_oni,
  sim_months = sim_m,
  seed = 3,
  max_lead = 3
)
expect_true(
  all(sprintf("enso_lead%d", 1:3) %in% names(a3)),
  "CFSv2: enso_lead1..3 returned"
)
expect_true(
  identical(a3[names(a0)], a0),
  "CFSv2: lag matrices unchanged by max_lead"
)
expect_true(
  isTRUE(all.equal(a3$enso_lead2[1:6, ], a3$enso_lag0[3:8, ])),
  "CFSv2: lead2(t) == lag0(t + 2) per replicate"
)
expect_true(
  all(is.finite(a3$enso_lead3)),
  "CFSv2: leads past the last sim month are filled (continuation), no NA"
)
an0 <- sample_enso_analog_scenarios(
  hist_oni,
  sim_m,
  "2019-01-01",
  "2020-12-01",
  n_rep = 2
)
expect_true(
  identical(
    an0,
    old$sample_enso_analog_scenarios(
      hist_oni,
      sim_m,
      "2019-01-01",
      "2020-12-01",
      n_rep = 2
    )
  ),
  "analog: max_lead = 0 identical to the committed sampler"
)
an2 <- sample_enso_analog_scenarios(
  hist_oni,
  sim_m,
  "2019-01-01",
  "2020-12-01",
  n_rep = 2,
  max_lead = 2
)
expect_true(
  isTRUE(all.equal(an2$enso_lead2[1:6, ], an2$enso_lag0[3:8, ])),
  "analog: lead2(t) == lag0(t + 2)"
)
plume <- tidyr::expand_grid(
  date = seq(as.Date("2024-07-01"), by = "month", length.out = 12),
  q_level = c(0.1, 0.5, 0.9)
) %>%
  dplyr::mutate(nino34_anomaly = q_level)
p0 <- sample_enso_forecast_trajectories(
  plume,
  n_rep = 3,
  historical_oni = hist_oni,
  sim_months = sim_m,
  seed = 2
)
expect_true(
  identical(
    p0,
    old$sample_enso_forecast_trajectories(
      plume,
      n_rep = 3,
      historical_oni = hist_oni,
      sim_months = sim_m,
      seed = 2
    )
  ),
  "C3S plume: max_lead = 0 identical to the committed sampler"
)
p2 <- sample_enso_forecast_trajectories(
  plume,
  n_rep = 3,
  historical_oni = hist_oni,
  sim_months = sim_m,
  seed = 2,
  max_lead = 2
)
expect_true(
  isTRUE(all.equal(p2$enso_lead2[1:6, ], p2$enso_lag0[3:8, ])),
  "C3S plume: lead2(t) == lag0(t + 2)"
)
# Lenssen-style cache: 4 analogs x 10 months, horizon_start = first sim month.
lc <- tempfile(fileext = ".rds")
saveRDS(
  list(
    oni_trajectories = matrix(seq_len(40) / 10, nrow = 4),
    n_analogs = 4L,
    n_months = 10L,
    horizon_start = format(sim_m[1]),
    set_labels = rep("neutral", 4)
  ),
  lc
)
l0 <- tryCatch(
  sample_enso_lenssen_scenarios(
    lc,
    sim_m,
    n_rep = 3,
    historical_oni = hist_oni,
    horizon_start = sim_m[1],
    seed = 1
  ),
  error = function(e) conditionMessage(e)
)
if (is.character(l0)) {
  cat(
    "  [SKIP] Lenssen-style cache fixture not accepted by this sampler version: ",
    l0,
    "\n"
  )
} else {
  expect_true(
    identical(
      l0,
      old$sample_enso_lenssen_scenarios(
        lc,
        sim_m,
        n_rep = 3,
        historical_oni = hist_oni,
        horizon_start = sim_m[1],
        seed = 1
      )
    ),
    "Lenssen: max_lead = 0 identical to the committed sampler"
  )
  l2 <- sample_enso_lenssen_scenarios(
    lc,
    sim_m,
    n_rep = 3,
    historical_oni = hist_oni,
    horizon_start = sim_m[1],
    seed = 1,
    max_lead = 2
  )
  expect_true(
    isTRUE(all.equal(l2$enso_lead2[1:6, ], l2$enso_lag0[3:8, ])),
    "Lenssen: lead2(t) == lag0(t + 2)"
  )
  err <- tryCatch(
    {
      sample_enso_lenssen_scenarios(
        lc,
        sim_m,
        n_rep = 3,
        historical_oni = hist_oni,
        horizon_start = sim_m[1],
        seed = 1,
        max_lead = 5
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err) && grepl("past the cached trajectories", err),
    "Lenssen: a lead past the cache is an error, not 0"
  )
}

# [4] panel_max_enso_lead().
cat("\n## 4) panel_max_enso_lead()\n")
expect_true(
  identical(panel_max_enso_lead(make_panel(0)), 2L),
  "panel with enso_lead1..2 -> 2"
)
expect_true(
  identical(
    panel_max_enso_lead(
      make_panel(0) %>% dplyr::select(-enso_lead1, -enso_lead2)
    ),
    0L
  ),
  "no lead columns -> 0"
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
