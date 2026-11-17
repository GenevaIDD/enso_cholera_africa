# test_simulator_kernel.R — HARD: unit tests for the forward kernel
#   (utils/simulator_core.R) and its frame / panel helpers.
#
# Usage:
#   Rscript analysis/tests/test_simulator_kernel.R
#
# Details:
#   Synthetic fixtures, no real fit object: the AR kernel, contagion lags,
#   TTLO, evaluate_lp(), simulate_one_tu(), the panel adapter and
#   build_location_frame(). Each test prints PASS/FAIL; the script exits
#   non-zero on failure.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R") # transitive dependency of provenance.R, sim_io.R (utils/README.md)
source("analysis/utils/gid.R") # transitive dependency of clim_compile.R, clim_db.R, occ_data.R, siz_data.R (utils/README.md)
source("analysis/utils/paths.R") # transitive dependency of provenance.R, sim_io.R (utils/README.md)
source("analysis/utils/clim_compile.R") # transitive dependency of occ_data.R, siz_data.R (utils/README.md)
source("analysis/utils/clim_db.R") # transitive dependency of occ_data.R, siz_data.R (utils/README.md)
source("analysis/utils/occ_data.R") # transitive dependency of sim_io.R (utils/README.md)
source("analysis/utils/siz_data.R") # transitive dependency of occ_data.R (utils/README.md)
source("analysis/utils/provenance.R") # transitive dependency of sim_io.R / anchor.R (utils/README.md)
source("analysis/utils/sim_io.R") # build_task_plan() -- used directly
source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/simulator_core.R")
source("analysis/utils/anchor.R")

# --------------------------------------------------------------------------
# Test harness ------------------------------------------------------------
# --------------------------------------------------------------------------

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# 09a_simulator_core kernel tests\n\n")

# --------------------------------------------------------------------------
# Fixture: 2 TUs, 2 countries per TU, 2 ADM1 per country, 2 ADM2 per ADM1
# = 16 locations total. We'll restrict to one TU = 8 locations.
# --------------------------------------------------------------------------

fixture_frame <- function() {
  expand.grid(
    tu = c("hmmW", "hmmCE"),
    country = c("AAA", "BBB"),
    adm1n = c(1, 2),
    adm2n = c(1, 2),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      country = paste0(country, substr(tu, 4, 4)), # keep countries unique
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", substr(tu, 4, 5)),
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
}

cat("## 1) build_location_frame derives adm1_unit + TU\n")
{
  frame <- fixture_frame()
  expect_equal(length(unique(frame$adm1_unit)), 8, "8 adm1 units")
  expect_equal(length(unique(frame$country)), 4, "4 countries")
  expect_equal(length(unique(frame$tu)), 2, "2 TUs")
  expect_true(
    all(grepl("^[A-Z]{3}[A-Z]_[0-9]+$", frame$adm1_unit)),
    "adm1_unit matches expected regex"
  )
}

# --------------------------------------------------------------------------
# Sub-frame: restrict to TU = hmmW so we can target the contagion math.
# --------------------------------------------------------------------------

frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
cat(sprintf(
  "\n## 2) update_lags: ADM1 / country / TU OR with leave-one-out (%d locs)\n",
  nrow(frame_tu)
))
{
  N <- nrow(frame_tu)
  # All zero -> all zero
  Y0 <- rep(0L, N)
  L0 <- update_lags(Y0, frame_tu)
  expect_true(all(L0$o == 0), "o all 0 when Y all 0")
  expect_true(all(L0$o_ct == 0), "o_ct all 0 when Y all 0")
  expect_true(all(L0$o_tu == 0), "o_tu all 0 when Y all 0")

  # Trigger one location in country AAAW, ADM1 #1
  Y1 <- Y0
  trig <- which(frame_tu$country == "AAAW" & frame_tu$adm1_unit == "AAAW_1")[1]
  Y1[trig] <- 1L
  L1 <- update_lags(Y1, frame_tu)
  # o = 1 for both ADM2 in AAAW_1, 0 elsewhere
  exp_o <- as.integer(frame_tu$adm1_unit == "AAAW_1")
  expect_equal(L1$o, exp_o, "o = 1 for both ADM2 in triggered ADM1")
  # o_ct: leave-one-ADM1-out. Locations in AAAW_1 see other-ADM1 in AAAW (none active)
  #       => o_ct = 0. Locations in AAAW_2 see AAAW_1 active => o_ct = 1.
  #       Locations in BBBW see no AAAW-active in their own country => 0.
  exp_o_ct <- as.integer(
    frame_tu$country == "AAAW" & frame_tu$adm1_unit != "AAAW_1"
  )
  expect_equal(L1$o_ct, exp_o_ct, "o_ct leave-one-ADM1-out within country")
  # o_tu: leave-one-country-out. AAAW is active in TU, so BBBW locations see o_tu=1.
  #       AAAW locations see "other countries of TU active" — BBBW not active => 0.
  exp_o_tu <- as.integer(frame_tu$country == "BBBW")
  expect_equal(L1$o_tu, exp_o_tu, "o_tu leave-one-country-out within TU")
}

# --------------------------------------------------------------------------
cat("\n## 3) advance_ttlo: reset on Y=1, advance by step_days otherwise\n")
{
  step <- 365.25 / 12
  u <- c(0.0, 547.5, 1095.0) # 0 d, ~18 mo, ~36 mo
  Y <- c(0L, 1L, 0L)
  u_next <- advance_ttlo(u, Y, step_days = step)
  expect_equal(
    u_next,
    c(0.0 + step, 0.0, 1095.0 + step),
    "TTLO resets to 0 on Y=1, advances by step_days otherwise"
  )
  # Default step_days = 365.25/12
  u2 <- advance_ttlo(c(100, 200), c(1L, 0L))
  expect_equal(u2, c(0.0, 200 + 365.25 / 12), "default step_days = 365.25/12")
}

# --------------------------------------------------------------------------
cat("\n## 4) evaluate_lp with zero coefs -> zero LP\n")
{
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 3,
    fe_overrides = NULL,
    basehaz_const = 0
  )
  N <- nrow(frame_tu)
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )
  lp <- evaluate_lp(
    draws,
    s = 1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u = rep(NA_integer_, N),
    frame_tu = frame_tu,
    basehaz_val = basehaz_eval(draws, 0, N)
  )
  expect_true(all(lp == 0), "LP = 0 with zero coefs and zero covars")
}

# --------------------------------------------------------------------------
cat("\n## 5) evaluate_lp: spillover slope amplifies o\n")
{
  N <- nrow(frame_tu)
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -3, ongoing_outbreak_lag1 = 2),
    basehaz_const = 0
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )
  lp_off <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u = rep(NA_integer_, N),
    frame_tu = frame_tu
  )
  lp_on <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(1L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u = rep(NA_integer_, N),
    frame_tu = frame_tu
  )
  expect_equal(lp_off, rep(-3, N), "LP = intercept when o=0")
  expect_equal(lp_on - lp_off, rep(2, N), "LP gap = ongoing slope when o=1")
}

# --------------------------------------------------------------------------
cat(
  "\n## 6) simulate_one_tu: zero hazard -> zero events; high hazard -> all events\n"
)
{
  # Build a fake 12-month panel for the TU
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  months <- as.Date(sprintf("2024-%02d-01", 1:12))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA
    )

  horizon <- list(start = "2024-01-01", end = "2024-12-31")
  outdir <- tempfile("sim_kernel_test_")
  dir.create(outdir, recursive = TRUE)

  ## 6a) zero hazard (LP = -1e6) -> never trigger
  draws_occ_lo <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -1e6)
  )
  res_lo <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ_lo,
    replicate_ids = 1L:2L,
    seed = 42L,
    output_part = file.path(outdir, "lo"),
    horizon = horizon
  )
  expect_equal(
    nrow(res_lo$event_table),
    0L,
    "no events fire when LP very negative"
  )
  expect_true(
    file.exists(file.path(outdir, "lo", "_SUCCESS")),
    "_SUCCESS sentinel written"
  )
  expect_true(
    file.exists(file.path(outdir, "lo", "event_table.qs")),
    "event_table.qs written"
  )
  expect_true(
    file.exists(file.path(outdir, "lo", "monthly_summary.qs")),
    "monthly_summary.qs written"
  )

  ## 6b) huge hazard (LP = +20) -> nearly every (loc, month) fires
  draws_occ_hi <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 20)
  )
  res_hi <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ_hi,
    replicate_ids = 1L:2L,
    seed = 42L,
    output_part = file.path(outdir, "hi"),
    horizon = horizon
  )
  # 2 replicates x 8 locs x 12 months = 192 cells
  expect_equal(
    nrow(res_hi$event_table),
    192L,
    "events fire in every cell with very high hazard"
  )
  # Size track retired: no episode-level size draw, no size
  # columns on the event table or the monthly summary.
  expect_true(
    !any(c("size_cases", "posterior_draw_s") %in% names(res_hi$event_table)),
    "event_table carries no size_cases / posterior_draw_s"
  )
  expect_equal(
    names(res_hi$monthly_summary),
    c("replicate_id", "gid", "month", "Y"),
    "legacy monthly_summary columns are exactly replicate_id/gid/month/Y"
  )
  ass <- qs2::qs_read(file.path(outdir, "hi", "assignment.qs"))
  expect_true(
    is.null(ass$posterior_idx_size) && !is.null(ass$posterior_idx_occ),
    "assignment.qs carries posterior_idx_occ and no posterior_idx_size"
  )
}

# --------------------------------------------------------------------------
cat("\n## 7) simulate_one_tu reproducibility: same seed -> identical outputs\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 4,
    fe_overrides = c(`(Intercept)` = -1)
  ) # moderate

  out1 <- tempfile("sim_repro_a_")
  out2 <- tempfile("sim_repro_b_")
  r1 <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:3L,
    seed = 7L,
    output_part = out1,
    horizon = horizon,
    dry_run = TRUE
  )
  r2 <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:3L,
    seed = 7L,
    output_part = out2,
    horizon = horizon,
    dry_run = TRUE
  )
  expect_equal(
    r1$event_table,
    r2$event_table,
    "identical event_table for same seed"
  )
  expect_equal(
    r1$monthly_summary,
    r2$monthly_summary,
    "identical monthly_summary for same seed"
  )
}

# --------------------------------------------------------------------------
cat("\n## 8) step-3 panel adapter: counterfactual scenario errors out\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  err <- tryCatch(
    build_climate_panel_for_tu(
      panel,
      frame_tu,
      horizon = list(start = "2024-01-01", end = "2024-06-30"),
      scenario = "counterfactual"
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    grepl("v2 feature", err),
    "counterfactual scenario raises a v2-deferral error"
  )
}

cat("\n## 9) step-3 panel adapter: horizon past panel coverage errors\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  err <- tryCatch(
    build_climate_panel_for_tu(
      panel,
      frame_tu,
      horizon = list(start = "2024-01-01", end = "2025-12-31")
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    grepl("panel ends|gap > 1 month", err),
    "horizon end past panel coverage raises a clear error"
  )
}

cat(
  "\n## 10) step-3 panel adapter: warm-up tail kept inside the returned slice\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(c(
    sprintf("2023-%02d-01", 1:12),
    sprintf("2024-%02d-01", 1:6)
  ))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  out <- build_climate_panel_for_tu(
    panel,
    frame_tu,
    horizon = list(start = "2024-01-01", end = "2024-06-30"),
    warmup_months = 3L
  )
  expect_true(
    min(out$t_start) < as.Date("2024-01-01"),
    "warm-up months land in the returned slice"
  )
  expect_true(
    min(out$t_start) >= as.Date("2023-09-01"),
    "warm-up trim respects warmup_months"
  )
}

cat("\n## 11) derive_min_seas_from_panel inverts the seasonal offset\n")
{
  # Three gids with min_seas = 1, 4, 7. Panel rows span one full season-year.
  build_loc <- function(gid, ms) {
    # t_start = first of each month, calendar; t_start_seas mirrors make_surv_data
    months_idx <- 1:12
    cal_first <- as.Date(sprintf("2024-%02d-01", months_idx))
    seas_first <- cal_first %m+% months(-(ms - 1)) # days_seas
    yd <- as.integer(lubridate::yday(seas_first)) - 1L
    tibble::tibble(gid = gid, t_start = cal_first, t_start_seas = yd)
  }
  panel <- dplyr::bind_rows(
    build_loc("AAA_1", 1L),
    build_loc("BBB_1", 4L),
    build_loc("CCC_1", 7L)
  )
  ms <- derive_min_seas_from_panel(panel)
  expect_equal(
    sort(ms$gid),
    c("AAA_1", "BBB_1", "CCC_1"),
    "all 3 gids returned"
  )
  expect_equal(ms$min_seas[ms$gid == "AAA_1"], 1L, "min_seas = 1 recovered")
  expect_equal(ms$min_seas[ms$gid == "BBB_1"], 4L, "min_seas = 4 recovered")
  expect_equal(ms$min_seas[ms$gid == "CCC_1"], 7L, "min_seas = 7 recovered")
}

cat("\n## 12) load_hmm_tu_mapping reads a CSV with iso3 + cluster columns\n")
{
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c(
      "iso3,cluster,extra",
      "AAA,hmmW,x",
      "BBB,hmmW,y",
      "CCC,hmmCE,z",
      "DDD,hmmS,w"
    ),
    tmp
  )
  m <- load_hmm_tu_mapping(tmp)
  expect_equal(sort(m$iso3), c("AAA", "BBB", "CCC", "DDD"), "iso3 list")
  expect_equal(m$tu[m$iso3 == "AAA"], "hmmW", "AAA → hmmW")
  expect_equal(m$tu[m$iso3 == "CCC"], "hmmCE", "CCC → hmmCE")
  # missing file -> informative error
  err <- tryCatch(load_hmm_tu_mapping("/nonexistent.csv"), error = function(e) {
    conditionMessage(e)
  })
  expect_true(
    grepl("HMM clusters file not found", err),
    "missing HMM file raises a clear error"
  )
}

cat(
  "\n## 13) load_population_per_location on a synthetic size RDS (load_mean_rate_frozen() is gone)\n"
)
{
  size_dat <- tibble::tibble(
    gid = c("AAA_1_1", "AAA_1_1", "BBB_2_1", "CCC_3_1"),
    pop = c(1e5, 1.1e5, 5e4, 2e5),
    suspected_cases = c(50, 60, 20, 400),
    admin_level = factor(
      c("ADM2", "ADM2", "ADM2", "ADM2"),
      levels = c("ADM0", "ADM1", "ADM2")
    ),
    TL = as.Date(c("2010-01-01", "2011-01-01", "2012-01-01", "2013-01-01")),
    TR = as.Date(c("2010-02-01", "2011-02-01", "2012-02-01", "2013-02-01"))
  )
  tmp <- tempfile(fileext = ".rds")
  saveRDS(size_dat, tmp)

  pop_map <- load_population_per_location(tmp)
  expect_equal(sort(pop_map$gid), c("AAA_1_1", "BBB_2_1", "CCC_3_1"), "3 gids")
  expect_equal(
    pop_map$pop[pop_map$gid == "AAA_1_1"],
    mean(c(1e5, 1.1e5)),
    "AAA_1_1 pop = mean over two rows"
  )

  # Size track retired downstream: the frozen mean_rate loader is not defined.
  expect_true(
    !exists("load_mean_rate_frozen"),
    "load_mean_rate_frozen() no longer defined"
  )

  # negative path
  err <- tryCatch(
    load_population_per_location("/nonexistent.rds"),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    grepl("size data file not found", err),
    "missing size data file raises a clear error"
  )
}

cat(
  "\n## 14) build_location_frame rehydrates min_seas / tu / pop / admin_level from paths (no mean_rate)\n"
)
{
  # Synthetic panel that LACKS min_seas/tu/pop/admin_level but
  # carries t_start + t_start_seas so min_seas can be derived.
  panel <- tibble::tibble(
    gid = rep(c("AAA_1_1", "BBB_2_1"), each = 12),
    country = rep(c("AAA", "BBB"), each = 12),
    AFRO_region = "REG_A",
    t_start = rep(as.Date(sprintf("2024-%02d-01", 1:12)), times = 2),
    t_start_seas = rep(
      as.integer(lubridate::yday(as.Date(sprintf("2024-%02d-01", 1:12)))) - 1L,
      times = 2
    )
  )

  hmm_csv <- tempfile(fileext = ".csv")
  writeLines(c("iso3,cluster", "AAA,hmmW", "BBB,hmmCE"), hmm_csv)

  size_dat <- tibble::tibble(
    gid = c("AAA_1_1", "BBB_2_1"),
    pop = c(2e5, 3e5),
    suspected_cases = c(40, 90),
    admin_level = factor("ADM2", levels = c("ADM0", "ADM1", "ADM2")),
    TL = as.Date("2018-01-01"),
    TR = as.Date("2018-06-01")
  )
  size_rds <- tempfile(fileext = ".rds")
  saveRDS(size_dat, size_rds)

  frame <- build_location_frame(
    panel,
    hmm_clusters_path = hmm_csv,
    size_data_path = size_rds
  )
  expect_equal(nrow(frame), 2L, "2 locations in frame")
  expect_true(all(!is.na(frame$min_seas)), "min_seas joined for every location")
  expect_true(
    all(frame$min_seas == 1L),
    "min_seas derived as 1 (calendar season)"
  )
  expect_equal(sort(frame$tu), c("hmmCE", "hmmW"), "tu joined from HMM CSV")
  expect_equal(
    frame$pop[frame$gid == "AAA_1_1"],
    2e5,
    "pop joined from size RDS"
  )
  expect_equal(
    unique(frame$admin_level),
    "ADM2",
    "admin_level joined from size RDS (mode per gid)"
  )

  expect_true(
    !any(c("mean_rate", "mr_tier") %in% names(frame)),
    "no mean_rate / mr_tier column on the frame (size track retired)"
  )
  err <- tryCatch(
    {
      build_location_frame(
        panel,
        size_data_path = size_rds,
        mean_rate_table_path = "x.rds"
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err) && grepl("unused argument", err),
    "mean_rate_table_path is no longer an argument"
  )
}

cat(
  "\n## 15) build_location_frame: missing TU rows logged but TU lag still inactive\n"
)
{
  # CCC is absent from the HMM csv → tu should be NA for CCC locations.
  panel <- tibble::tibble(
    gid = c("AAA_1_1", "CCC_9_1"),
    country = c("AAA", "CCC"),
    AFRO_region = c("REG_A", "REG_C"),
    t_start = as.Date(c("2024-01-01", "2024-01-01")),
    t_start_seas = c(0L, 0L)
  )
  hmm_csv <- tempfile(fileext = ".csv")
  writeLines(c("iso3,cluster", "AAA,hmmW"), hmm_csv) # CCC missing
  size_rds <- tempfile(fileext = ".rds")
  saveRDS(
    tibble::tibble(
      gid = c("AAA_1_1", "CCC_9_1"),
      pop = c(1e5, 1e5),
      suspected_cases = c(10, 10),
      admin_level = factor("ADM2", levels = c("ADM0", "ADM1", "ADM2")),
      TL = as.Date("2018-01-01"),
      TR = as.Date("2018-06-01")
    ),
    size_rds
  )

  frame <- build_location_frame(
    panel,
    hmm_clusters_path = hmm_csv,
    size_data_path = size_rds
  )
  expect_equal(frame$tu[frame$gid == "AAA_1_1"], "hmmW", "AAA mapped")
  expect_true(
    is.na(frame$tu[frame$gid == "CCC_9_1"]),
    "CCC NA (absent from HMM file)"
  )
}

cat(
  "\n## 16) build_location_frame: no population -> hard error; no admin_level / no anchors is legal\n"
)
{
  # Population has no default: a placeholder pop of 1 would put units into the
  # forecast with a population of one person, and against an imposed population
  # elasticity of 1.0 that IS the unit's scale. The assertions below encode the
  # contract, in the order the requirements are checked.
  panel <- tibble::tibble(
    gid = "AAA_1_1",
    country = "AAA",
    AFRO_region = "REG_A",
    t_start = as.Date("2024-01-01"),
    t_start_seas = 0L
  )

  # (a) Nothing supplied at all. Population is the first unmet requirement, so
  #     it is the one named.
  err <- tryCatch(
    {
      build_location_frame(panel)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err),
    "no population source -> hard error, not a silent pop = 1.0 default (D4)"
  )
  expect_true(
    !is.null(err) && grepl("gid_pop_path", err, fixed = TRUE),
    "error message names gid_pop_path -- it says what to pass, not just what is missing"
  )
  expect_true(
    !is.null(err) && grepl("one person", err, fixed = TRUE),
    "error message states the consequence the default used to have"
  )

  # (b) A canonical Q_g table, so the population requirement is met. With the
  #     size track's mean_rate gone, a frame with no admin_level and no anchor
  #     tables is LEGAL (the mon `admin_levelADM2` term is already guarded on
  #     the column's presence).
  gid_pop_qs <- tempfile(fileext = ".qs")
  qs2::qs_save(
    tibble::tibble(
      gid = "AAA_1_1",
      year = 2020L,
      Q = 250000,
      iso3 = "AAA",
      gadm_level = "ADM2"
    ),
    gid_pop_qs
  )
  err2 <- tryCatch(
    {
      suppressMessages(build_location_frame(panel, gid_pop_path = gid_pop_qs))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.null(err2),
    "with population resolved, no admin_level and no anchors is not an error"
  )
  # ...but an anchor table with no admin_level source still hard-errors.
  err2b <- tryCatch(
    {
      suppressMessages(build_location_frame(
        panel,
        gid_pop_path = gid_pop_qs,
        base_rate_table_paths = list(mon = "/nonexistent/brt.rds")
      ))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err2b) && grepl("admin_level", err2b),
    "an anchor table without any admin_level source still names admin_level"
  )

  # (c) An EXPLICITLY configured path that does not exist is fatal -- a config
  #     key is an instruction, unlike a defaulted path which may degrade.
  err3 <- tryCatch(
    {
      build_location_frame(panel, gid_pop_path = "/nonexistent/q.qs")
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err3) && grepl("not found", err3),
    "an explicit but missing gid_pop_path is fatal, never silently skipped"
  )

  # (d) Everything supplied: population comes from Q_g, tagged with its tier.
  al_data <- tibble::tibble(gid = "AAA_1_1", admin_level = "ADM2")
  out <- capture.output(
    frame2 <- build_location_frame(
      panel,
      admin_level_data = al_data,
      gid_pop_path = gid_pop_qs,
      pop_year = 2020L
    ),
    type = "message"
  )
  expect_equal(nrow(frame2), 1L, "single-row frame returned")
  expect_true(is.na(frame2$tu), "tu default NA when no HMM path")
  expect_equal(
    frame2$pop,
    250000,
    "pop resolved from the canonical Q_g table, not defaulted"
  )
  expect_equal(frame2$pop_tier, "direct", "pop_tier recorded on the frame")
  expect_equal(
    frame2$admin_level,
    "ADM2",
    "admin_level joined from admin_level_data"
  )
  expect_true(
    !("mean_rate" %in% names(frame2)),
    "no mean_rate column (size track retired)"
  )
  expect_true(
    any(grepl("rehydration", out)),
    "single-summary rehydration log emitted"
  )
  expect_true(
    any(grepl("tiers: direct=1", out, fixed = TRUE)),
    "the pop tier table is emitted into the rehydration log"
  )
}

# ==========================================================================
# Tests: TTLO factor-level roundtrip + NA-TU bucket routing
# ==========================================================================

cat("\n## 17) bin_ttlo: boundary correctness\n")
{
  # 0 days → (0+1)/365 = 0.00274 → [0,1)
  expect_equal(bin_ttlo(0), "[0,1)", "0 days → [0,1)")
  # 363 days → (363+1)/365 = 0.9973 → [0,1); 364 days = 365/365 = 1.0 → [1,2)
  expect_equal(bin_ttlo(363), "[0,1)", "363 days → [0,1)")
  # 365 days → (365+1)/365 = 1.0027 → [1,2)
  expect_equal(bin_ttlo(365), "[1,2)", "365 days → [1,2)")
  # 6*365=2190 days → (2190+1)/365 = 6.003 → [6,Inf]
  expect_equal(bin_ttlo(2190), "[6,Inf]", "2190 days → [6,Inf]")
  # 1000 years → [6,Inf]
  expect_equal(bin_ttlo(1e6), "[6,Inf]", "very large days → [6,Inf]")
  # NA input → NA output
  expect_true(is.na(bin_ttlo(NA_real_)), "NA days → NA bin")
}

cat("\n## 18) ttlo_init_from_panel: factor level → days roundtrip\n")
{
  levels <- c("[0,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,Inf]")
  days <- ttlo_init_from_panel(levels)
  # Every level must round-trip back to itself (lower bound initialises inside its bin)
  recovered <- bin_ttlo(days)
  expect_equal(
    recovered,
    levels,
    "all 7 levels round-trip bin_ttlo ∘ ttlo_init"
  )
  # NA input → top bin (6*365.25 days)
  na_days <- ttlo_init_from_panel(NA_character_)
  expect_equal(na_days, 6 * 365.25, "NA label → 6*365.25 days (top bin)")
}

cat("\n## 19) ttlo_col_to_level: brms column name → factor level\n")
{
  expect_equal(
    ttlo_col_to_level("time_to_last_outbreak12"),
    "[1,2)",
    "12 → [1,2)"
  )
  expect_equal(
    ttlo_col_to_level("time_to_last_outbreak23"),
    "[2,3)",
    "23 → [2,3)"
  )
  expect_equal(
    ttlo_col_to_level("time_to_last_outbreak6Inf"),
    "[6,Inf]",
    "6Inf → [6,Inf]"
  )
}

cat("\n## 20) evaluate_lp: TTLO factor-level FE contributes correctly\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )

  # Bundle with TTLO FE for levels [1,2) and [6,Inf]: fe = 1.5 and 3.0
  ttlo_terms <- c(
    "(Intercept)",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1",
    "time_to_last_outbreak12",
    "time_to_last_outbreak6Inf"
  )
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(
      `(Intercept)` = 0,
      time_to_last_outbreak12 = 1.5,
      time_to_last_outbreak6Inf = 3.0
    ),
    active_terms = ttlo_terms
  )

  # All locations in reference bin [0,1) → LP = 0
  u_ref <- rep(ttlo_init_from_panel("[0,1)"), N)
  lp_ref <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = u_ref,
    frame_tu = frame_tu
  )
  expect_equal(lp_ref, rep(0, N), "reference bin [0,1) → LP = 0")

  # All locations in [1,2) → LP = 1.5
  u_12 <- rep(ttlo_init_from_panel("[1,2)"), N)
  lp_12 <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = u_12,
    frame_tu = frame_tu
  )
  expect_equal(lp_12, rep(1.5, N), "[1,2) bin → LP = fe[12] = 1.5")

  # All locations in [6,Inf] → LP = 3.0
  u_top <- rep(ttlo_init_from_panel("[6,Inf]"), N)
  lp_top <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = u_top,
    frame_tu = frame_tu
  )
  expect_equal(lp_top, rep(3.0, N), "[6,Inf] bin → LP = fe[6Inf] = 3.0")

  # NA u_days → 0 contribution (warm-up not available)
  lp_na <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  expect_equal(lp_na, rep(0, N), "NA u_days → 0 TTLO contribution")
}

cat("\n## 21) update_lags: NA-TU locations get o_tu = NA\n")
{
  # Build a frame with two countries: one with TU, one without
  frame_mix <- tibble::tibble(
    gid = c("AAA_1_1", "AAA_1_2", "ZZZ_9_1", "ZZZ_9_2"),
    country = c("AAA", "AAA", "ZZZ", "ZZZ"),
    AFRO_region = "REG_A",
    adm1_unit = c("AAA_1", "AAA_1", "ZZZ_9", "ZZZ_9"),
    tu = c("hmmW", "hmmW", NA_character_, NA_character_),
    min_seas = 1L,
    pop = 1e5
  )

  Y0 <- rep(0L, 4)
  L0 <- update_lags(Y0, frame_mix)
  expect_true(all(L0$o == 0), "o all 0 when Y all 0")
  expect_true(all(L0$o_ct == 0), "o_ct all 0 when Y all 0")
  # NA-TU locations have NA o_tu; TU-assigned locations have 0 o_tu
  expect_true(all(L0$o_tu[1:2] == 0), "TU locations get o_tu = 0")
  expect_true(all(is.na(L0$o_tu[3:4])), "NA-TU locations get o_tu = NA")

  # Trigger AAA → TU location BBBW would see o_tu=1 (leave-one-out)
  # But in this frame, AAA is in hmmW and ZZZ is NA-TU, so no leave-one-out
  Y1 <- c(1L, 0L, 0L, 0L)
  L1 <- update_lags(Y1, frame_mix)
  # AAA active → o_tu for AAA = 0 (no other country in hmmW)
  expect_equal(L1$o_tu[1:2], c(0L, 0L), "AAA o_tu = 0 (only country in TU)")
  # ZZZ is NA-TU → o_tu stays NA
  expect_true(all(is.na(L1$o_tu[3:4])), "ZZZ NA-TU preserved after trigger")
}

cat("\n## 22) evaluate_lp: NA o_tu contributes 0 not NA to LP\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -2.0, ongoing_outbreak_tu_lag1 = 5.0)
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )
  # o_tu = NA for all locations → should not produce NA in LP
  lp_na_tu <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(NA_integer_, N),
    u_days = rep(0, N),
    frame_tu = frame_tu
  )
  expect_true(all(!is.na(lp_na_tu)), "NA o_tu produces no NA in LP")
  expect_equal(lp_na_tu, rep(-2.0, N), "NA o_tu contributes 0 (only intercept)")
}

cat(
  "\n## 23) build_task_plan: TU buckets, seed disjointness, posterior strategy\n"
)
{
  # Exercises the REAL shared planner (analysis/utils.R), not an inline copy, so
  # two failure modes are pinned: (A) the base seed must stride
  # by R = n_replicates so set.seed(seed + rep_id) yields disjoint ranges across
  # (TU, chunk) tasks; (B) per-track posterior indices must actually be emitted.
  # With the size track retired only posterior_strategy =
  # "paired" exists; anything else is refused.
  mk_opt <- function(strategy = "paired", R = 4L, cs = 2L) {
    list(
      simulator = list(
        n_replicates = R,
        parallel = list(chunk_size = cs),
        scope = list(transmission_units = "all"),
        seed = 42L,
        posterior_strategy = strategy
      )
    )
  }

  frame_with_na <- tibble::tibble(
    gid = c("AAA_1_1", "BBB_2_1", "CCC_9_1"),
    tu = c("hmmW", "hmmCE", NA_character_)
  )

  # --- TU bucketing (NA-TU → _NA_TU_ appended after sorted real TUs) ---------
  plan_na <- build_task_plan(mk_opt(), frame = frame_with_na)
  tus_na <- unique(vapply(plan_na, function(r) r$tu, character(1)))
  expect_equal(
    tus_na,
    c("hmmCE", "hmmW", "_NA_TU_"),
    "NA-TU adds _NA_TU_ bucket after sorted real TUs"
  )

  plan_no_na <- build_task_plan(mk_opt(), frame = frame_with_na[1:2, ])
  tus_no_na <- unique(vapply(plan_no_na, function(r) r$tu, character(1)))
  expect_equal(tus_no_na, c("hmmCE", "hmmW"), "no NA-TU → no _NA_TU_ bucket")

  # --- Bug A: seed disjointness across (TU, chunk) tasks ---------------------
  # The kernel runs set.seed(task$seed + rep_id) for rep_id in replicate_ids;
  # collect every such value across all tasks — they must be pairwise unique.
  eff_seeds <- unlist(lapply(plan_na, function(r) r$seed + r$replicate_ids))
  expect_equal(
    length(eff_seeds),
    length(unique(eff_seeds)),
    "effective per-replicate seeds are disjoint across tasks (stride = R)"
  )
  # Sanity: base seeds stride by exactly R = n_replicates per task.
  base_seeds <- vapply(plan_na, function(r) r$seed, numeric(1))
  expect_equal(
    diff(base_seeds),
    rep(4, length(base_seeds) - 1L),
    "base seed strides by R (=n_replicates) per task"
  )

  # --- Bug B: posterior_strategy ---------------------------------------------
  # paired: occ index == dur index == replicate ids for every task; no size
  # index is emitted any more.
  plan_paired <- build_task_plan(mk_opt("paired"), frame = frame_with_na)
  ok_paired <- all(vapply(
    plan_paired,
    function(r) {
      identical(r$posterior_idx_occ, r$replicate_ids) &&
        identical(r$posterior_idx_dur, r$posterior_idx_occ)
    },
    logical(1)
  ))
  expect_true(ok_paired, "paired: occ and dur indices equal the replicate ids")
  expect_true(
    all(vapply(
      plan_paired,
      function(r) is.null(r$posterior_idx_size),
      logical(1)
    )),
    "no task carries posterior_idx_size"
  )

  # independent: refused (it permuted the retired size track's draws).
  err <- tryCatch(
    {
      build_task_plan(mk_opt("independent"), frame = frame_with_na)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err) && grepl("paired", err),
    "posterior_strategy = 'independent' is refused"
  )
}

cat("\n## 24) simulate_one_tu: NA-TU bucket runs without NA LP contamination\n")
{
  # Frame with all NA-TU (simulates _NA_TU_ bucket)
  frame_na_tu <- fixture_frame() %>%
    dplyr::filter(tu == "hmmW") %>%
    dplyr::mutate(tu = NA_character_)

  N <- nrow(frame_na_tu)
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = frame_na_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_tu_lag1 = NA_integer_,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      time_to_last_outbreak_cat = NA_character_
    )
  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  draws_occ <- make_draws_bundle_synthetic(
    frame_na_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = -1, ongoing_outbreak_tu_lag1 = 10) # large tu slope → would blow up if NA not handled
  )
  res <- simulate_one_tu(
    tu = "_NA_TU_",
    frame_tu = frame_na_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:2L,
    seed = 99L,
    output_part = tempfile("na_tu_test_"),
    horizon = horizon,
    dry_run = TRUE
  )
  expect_true(!is.null(res), "_NA_TU_ bucket returns a result")
  expect_true(
    !is.null(res$monthly_summary),
    "monthly_summary present for NA-TU run"
  )
  # No NA in the Y column — LP should have been finite
  expect_true(
    all(!is.na(res$monthly_summary$Y)),
    "no NA Y values with NA-TU bucket (LP was finite)"
  )
}

# ==========================================================================
# Real-fit naming tests: the synthetic fixtures above all use clean covariate
# names, so the brms factor-level suffix (`...TRUE`) and the extra ENSO FE were
# never exercised — the exact gap that let the simulator silently drop the
# spillover lags and ENSO. These tests feed a hand-built `as_draws_df`-shaped
# frame with real brms column names through the extractor's parser.
# ==========================================================================

cat(
  "\n## 25) .occ_bundle_from_draws_df: brms TRUE-suffix + enso names parse correctly\n"
)
{
  # Column labels must match what re_level_key() generates: brms pastes all
  # ancestor values with "_" (spaces → "."), so nested labels are composite.
  # region="REG_W", country="REG_W_AAAW", adm1="REG_W_AAAW_AAAW_1",
  # adm2/gid="REG_W_AAAW_AAAW_1_AAAW_1_1".
  draws_df <- data.frame(
    check.names = FALSE,
    `.draw` = 1:3,
    `b_Intercept` = c(-3, -3, -3),
    `b_ongoing_outbreak_lag1TRUE` = c(2.0, 2.1, 1.9),
    `b_ongoing_outbreak_country_lag1TRUE` = c(0.5, 0.4, 0.6),
    `b_enso_lag5` = c(0.3, 0.2, 0.4),
    # enso slope nested at region / region:country / region:country:ADM1
    `r_AFRO_region[REG_W,enso_lag5]` = c(0.01, 0.02, 0.03),
    `r_AFRO_region:country[REG_W_AAAW,enso_lag5]` = c(0.04, 0.05, 0.06),
    `r_AFRO_region:country:ADM1[REG_W_AAAW_AAAW_1,enso_lag5]` = c(
      0.07,
      0.08,
      0.09
    ),
    # spillover slope nested at region:country, carrying the TRUE suffix
    `r_AFRO_region:country[REG_W_AAAW,ongoing_outbreak_lag1TRUE]` = c(
      0.1,
      0.2,
      0.3
    ),
    # nested intercept RE at the ADM2 (gid) level
    `r_AFRO_region:country:ADM1:gid[REG_W_AAAW_AAAW_1_AAAW_1_1,Intercept]` = c(
      -0.2,
      -0.1,
      0.0
    )
  )
  bundle <- .occ_bundle_from_draws_df(draws_df)

  expect_true(
    "ongoing_outbreak_lag1" %in% colnames(bundle$fe),
    "FE TRUE-suffix stripped: ongoing_outbreak_lag1"
  )
  expect_true(
    "ongoing_outbreak_country_lag1" %in% colnames(bundle$fe),
    "FE TRUE-suffix stripped: ongoing_outbreak_country_lag1"
  )
  expect_true("enso_lag5" %in% colnames(bundle$fe), "enso_lag5 FE retained")
  expect_true(
    !("ongoing_outbreak_lag1TRUE" %in% colnames(bundle$fe)),
    "no residual TRUE-suffixed FE column"
  )
  # ENSO slope present at all three nested buckets
  expect_true(
    !is.null(bundle$re_slope$enso_lag5$region),
    "enso slope: region bucket"
  )
  expect_true(
    !is.null(bundle$re_slope$enso_lag5$country),
    "enso slope: country bucket"
  )
  expect_true(
    !is.null(bundle$re_slope$enso_lag5$adm1),
    "enso slope: adm1 bucket"
  )
  # Spillover slope key normalised (TRUE stripped) in re_slope
  expect_true(
    !is.null(bundle$re_slope$ongoing_outbreak_lag1$country),
    "spillover slope key normalised in re_slope"
  )
  expect_true(
    !is.null(bundle$re_intercept$adm2),
    "ADM2 (gid) intercept RE parsed"
  )
}

cat("\n## 26) evaluate_lp: enso FE + nested adm1 slope enters the LP\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  # Column names must be composite re_level_key keys, not bare adm1_unit values,
  # because nested_slope() uses re_level_key(frame_tu, "adm1") for the lookup.
  adm1_keys <- re_level_key(frame_tu, "adm1")
  adm1_lv <- unique(adm1_keys)
  aaaw1_key <- unique(adm1_keys[frame_tu$adm1_unit == "AAAW_1"])

  enso_terms <- c(
    "(Intercept)",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1",
    "enso_lag5"
  )
  # adm1 slope: AAAW_1 gets +1.0, all other adm1 units 0
  slope_mat <- matrix(
    0,
    nrow = 2,
    ncol = length(adm1_lv),
    dimnames = list(NULL, adm1_lv)
  )
  slope_mat[, aaaw1_key] <- 1.0
  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(`(Intercept)` = 0, enso_lag5 = 0.5),
    re_slope_overrides = list(enso_lag5 = list(adm1 = slope_mat)),
    active_terms = enso_terms
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N),
    enso_lag5 = rep(2, N)
  )
  lp <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  # AAAW_1 locations: (0.5 + 1.0) * 2 = 3.0; all others: 0.5 * 2 = 1.0
  exp_lp <- ifelse(frame_tu$adm1_unit == "AAAW_1", 3.0, 1.0)
  expect_equal(
    lp,
    exp_lp,
    "enso FE + adm1 slope contribute (3.0 in AAAW_1, 1.0 elsewhere)"
  )
}

cat(
  "\n## 27) .parse_brms_re_draws: mon-shaped (NB count) ADM2 nesting maps to the adm2 bucket\n"
)
{
  # The NB count models (monthly `mon`; formerly also the retired size track)
  # nest by AFRO_region/country/ADM1/ADM2 (ADM2 == gid), so their deepest r_*
  # group is `...:ADM1:ADM2`, unlike the occurrence track's `...:ADM1:gid`. The
  # composite level labels are identical to re_level_key() output. Regression
  # for the dropped-RE bug: the count extractor once hard-coded
  # re_intercept = list(region=NULL,...), so the ADM2 group was never parsed and
  # the predicted count collapsed to the FE-only mean.
  draws_df <- data.frame(
    check.names = FALSE,
    `.draw` = 1:3,
    `b_Intercept` = c(1, 1, 1),
    `b_mean_precip` = c(0.2, 0.3, 0.4),
    `r_AFRO_region[REG_W,Intercept]` = c(0.1, 0.1, 0.1),
    `r_AFRO_region:country[REG_W_AAAW,Intercept]` = c(0.2, 0.2, 0.2),
    `r_AFRO_region:country:ADM1[REG_W_AAAW_AAAW_1,Intercept]` = c(
      0.3,
      0.3,
      0.3
    ),
    `r_AFRO_region:country:ADM1:ADM2[REG_W_AAAW_AAAW_1_AAAW_1_1,Intercept]` = c(
      0.5,
      0.6,
      0.7
    ),
    `shape` = c(4, 5, 6)
  )
  parsed <- .parse_brms_re_draws(draws_df)
  expect_true(
    !is.null(parsed$re_intercept$region),
    "mon: region intercept parsed"
  )
  expect_true(
    !is.null(parsed$re_intercept$country),
    "mon: country intercept parsed"
  )
  expect_true(!is.null(parsed$re_intercept$adm1), "mon: adm1 intercept parsed")
  expect_true(
    !is.null(parsed$re_intercept$adm2),
    "mon: ADM2-group intercept parsed into adm2 bucket (was dropped before)"
  )
  expect_equal(
    colnames(parsed$re_intercept$adm2),
    "REG_W_AAAW_AAAW_1_AAAW_1_1",
    "mon: adm2 colname is the composite re_level_key key"
  )
  expect_true("mean_precip" %in% colnames(parsed$fe), "mon: FE parsed")
}

cat(
  "\n## 28) mon-shaped REs enter the LP: positive ADM2 intercept lifts the predicted count\n"
)
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  adm2_keys <- re_level_key(frame_tu, "adm2")
  # +2.0 intercept on a single gid (AAAW_1_1), 0 elsewhere
  target_gid <- "AAAW_1_1"
  target_key <- adm2_keys[frame_tu$gid == target_gid]
  re_adm2 <- matrix(
    0,
    nrow = 2,
    ncol = length(adm2_keys),
    dimnames = list(NULL, adm2_keys)
  )
  re_adm2[, target_key] <- 2.0
  draws_mon <- list(
    S = 2L,
    fe = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "Intercept")),
    re_intercept = list(
      region = NULL,
      country = NULL,
      adm1 = NULL,
      adm2 = re_adm2
    ),
    re_slope = list(),
    basehaz = list(kind = "none"),
    phi = c(5, 5),
    active_terms = "Intercept"
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )
  lp <- evaluate_lp(
    draws_mon,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu
  )
  exp_lp <- ifelse(frame_tu$gid == target_gid, 2.0, 0.0)
  expect_equal(
    lp,
    exp_lp,
    "ADM2 intercept RE adds to the mon LP for its gid only"
  )
}

# ==========================================================================
# Tests: drop_country_tu_spillover leakage-fix flag
# ==========================================================================

cat("\n## 29) evaluate_lp: drop_ct_tu suppresses BOTH TRUE- and FALSE-branch\n")
cat("       o_ct/o_tu contributions (own-o unaffected)\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  country_keys <- re_level_key(frame_tu, "country")
  country_lv <- unique(country_keys)
  # FALSE-branch RE slope for o_ct: brms stores a nested slope even for the
  # FALSE factor level, so a naive "zero the input" gate would still leak
  # this term into the LP -- the flag must omit the whole add_spillover() call.
  false_mat <- matrix(
    5.0,
    nrow = 2,
    ncol = length(country_lv),
    dimnames = list(NULL, country_lv)
  )

  draws <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(
      `(Intercept)` = -2,
      ongoing_outbreak_lag1 = 100, # own-o: must NOT be gated
      ongoing_outbreak_country_lag1 = 3,
      ongoing_outbreak_tu_lag1 = 7
    ),
    re_slope_overrides = list(
      ongoing_outbreak_country_lag1FALSE = list(country = false_mat)
    )
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )

  # o = 1 (own, TRUE branch, must always contribute), o_ct = 0 (FALSE branch,
  # RE-slope-driven), o_tu = 1 (TRUE branch, FE-driven).
  lp_keep <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(1L, N),
    o_ct = rep(0L, N),
    o_tu = rep(1L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_ct_tu = FALSE
  )
  expect_equal(
    lp_keep,
    rep(-2 + 100 + 5.0 + 7, N),
    "drop_ct_tu=FALSE: own-o, o_ct FALSE-branch RE, and o_tu TRUE-branch FE all contribute"
  )

  lp_drop <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(1L, N),
    o_ct = rep(0L, N),
    o_tu = rep(1L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_ct_tu = TRUE
  )
  expect_equal(
    lp_drop,
    rep(-2 + 100, N),
    "drop_ct_tu=TRUE: o_ct/o_tu contribution fully suppressed (own-o untouched)"
  )
}

cat(
  "\n## 30) simulate_one_tu: drop_country_tu_spillover zeroes warm o_ct/o_tu\n"
)
cat("       AND every step's LP contribution (end-to-end kernel test)\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  N <- nrow(frame_tu)
  # One warm-up month (pre-horizon) seeding o_ct/o_tu = 1 for every unit, plus
  # one in-horizon month to observe the effect on. Single-month horizon keeps
  # this deterministic (no leave-one-out feedback from realised events).
  panel <- tidyr::expand_grid(
    gid = frame_tu$gid,
    t_start = as.Date(c("2023-12-01", "2024-01-01"))
  ) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = ifelse(
        t_start == as.Date("2023-12-01"),
        1L,
        0L
      ),
      ongoing_outbreak_tu_lag1 = ifelse(
        t_start == as.Date("2023-12-01"),
        1L,
        0L
      ),
      time_to_last_outbreak_cat = NA_character_
    )
  horizon <- list(start = "2024-01-01", end = "2024-01-31")

  # Intercept alone must never fire; country/tu FE alone must fire near-certainly
  # when their (warm-seeded) lag is 1 -- mirrors test 6's "-1e6 never / +20 always"
  # deterministic-hazard pattern.
  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 2,
    fe_overrides = c(
      `(Intercept)` = -1e6,
      ongoing_outbreak_country_lag1 = 1e6 + 20,
      ongoing_outbreak_tu_lag1 = 1e6 + 20
    )
  )

  out_keep <- tempfile("sim_dropct_keep_")
  res_keep <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:2L,
    seed = 11L,
    output_part = out_keep,
    horizon = horizon,
    drop_country_tu_spillover = FALSE
  )
  expect_equal(
    nrow(res_keep$event_table),
    N * 2L,
    "drop_country_tu_spillover=FALSE: warm o_ct/o_tu=1 fires every cell"
  )

  out_drop <- tempfile("sim_dropct_drop_")
  res_drop <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:2L,
    seed = 11L,
    output_part = out_drop,
    horizon = horizon,
    drop_country_tu_spillover = TRUE
  )
  expect_equal(
    nrow(res_drop$event_table),
    0L,
    "drop_country_tu_spillover=TRUE: warm o_ct/o_tu zeroed -> intercept-only LP never fires"
  )
}

cat(
  "\n## 31) simulate_one_tu: drop_country_tu_spillover default FALSE preserves\n"
)
cat("       byte-identical output to the pre-existing call signature\n")
{
  frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")
  months <- as.Date(sprintf("2024-%02d-01", 1:3))
  panel <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  horizon <- list(start = "2024-01-01", end = "2024-03-31")

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 3,
    fe_overrides = c(`(Intercept)` = -1)
  )

  r_implicit <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:2L,
    seed = 5L,
    output_part = tempfile("dropct_default_a_"),
    horizon = horizon,
    dry_run = TRUE
  )
  r_explicit <- simulate_one_tu(
    "hmmW",
    frame_tu,
    panel,
    draws_occ = draws_occ,
    replicate_ids = 1L:2L,
    seed = 5L,
    output_part = tempfile("dropct_default_b_"),
    horizon = horizon,
    dry_run = TRUE,
    drop_country_tu_spillover = FALSE
  )
  expect_equal(
    r_implicit$event_table,
    r_explicit$event_table,
    "omitting drop_country_tu_spillover == passing FALSE explicitly"
  )
  expect_equal(
    r_implicit$monthly_summary,
    r_explicit$monthly_summary,
    "omitting drop_country_tu_spillover == passing FALSE explicitly (monthly_summary)"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
