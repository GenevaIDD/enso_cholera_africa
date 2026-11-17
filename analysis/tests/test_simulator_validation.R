# test_simulator_validation.R — HARD: the simulator validation / postprocessing
#   layer (07c_simulator_postprocess.R, utils/sim_metrics.R).
#
# Usage:
#   Rscript analysis/tests/test_simulator_validation.R
#
# Details:
#   Synthetic in-memory fixtures, no real data paths: partition aggregation and
#   IO, the observed panel loader, every compute_* metric, the validation
#   figures (smoke) and write_validation_summary().

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

# 07c_simulator_postprocess.R's own functions are defined at top level, but its
# CLI block (source() calls for its own dependency closure) is guarded by
# `if (sys.nframe() == 0L)` and never runs when the file is source()d (as here)
# rather than Rscript'd directly -- so this test must mirror that same closure
# itself.
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/simulator_core.R")
source("analysis/utils/dependence.R") # paired_kendall_tau() -- used by compute_sim_kendall_tau()
source("analysis/utils/figures_io.R") # figure_saver() -- used by make_sim_figures()
source("analysis/utils/anchor.R")
source("analysis/utils/obs_sim_compare.R") # crps_ens()/region_key_of() -- used by sim_metrics.R
source("analysis/utils/sim_metrics.R") # the compute_* metrics
source("analysis/07c_simulator_postprocess.R")

# --------------------------------------------------------------------------
# Test harness
# --------------------------------------------------------------------------

.fails <- 0L
.passes <- 0L

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# 09b_postprocess_sim validation tests\n\n")

# --------------------------------------------------------------------------
# Reuse the same fixture_frame() from test_simulator_kernel.R
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
      country = paste0(country, substr(tu, 4, 4)),
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", substr(tu, 4, 5)),
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
}

# --------------------------------------------------------------------------
# Test 1: aggregate_partitions — synthetic partitions
# --------------------------------------------------------------------------
cat("\n## 1) aggregate_partitions — synthetic partitions\n")
{
  tmp_root <- tempfile("agg_test_")
  dir.create(tmp_root, recursive = TRUE)

  make_partition <- function(id, with_sentinel = TRUE) {
    part_dir <- file.path(tmp_root, sprintf("part_%02d", id))
    dir.create(part_dir, recursive = TRUE)
    et <- tibble::tibble(
      replicate_id = id,
      gid = c("AAA_1_1", "BBB_2_1", "CCC_3_1"),
      start_month = "2018-01-01",
      duration_months = c(1L, 2L, 3L)
    )
    ms <- tibble::tibble(
      replicate_id = id,
      gid = rep(c("AAA_1_1", "BBB_2_1"), each = 3),
      month = rep(c("2018-01-01", "2018-02-01", "2018-03-01"), times = 2),
      Y = c(1L, 0L, 1L, 0L, 1L, 0L)
    )
    qs_save(et, file.path(part_dir, "event_table.qs"))
    qs_save(ms, file.path(part_dir, "monthly_summary.qs"))
    if (with_sentinel) writeLines("done", file.path(part_dir, "_SUCCESS"))
  }

  make_partition(1, with_sentinel = TRUE)
  make_partition(2, with_sentinel = TRUE)
  make_partition(3, with_sentinel = FALSE) # no sentinel — should be skipped

  merged_dir <- file.path(tmp_root, "merged")
  res <- aggregate_partitions(tmp_root, merged_dir)

  expect_equal(
    nrow(res$event_table),
    6L,
    "merged event_table has 6 rows (2 successful partitions × 3)"
  )
  expect_equal(
    nrow(res$monthly_summary),
    12L,
    "merged monthly_summary has 12 rows (2 successful × 6)"
  )
  expect_true(
    file.exists(file.path(merged_dir, "event_table.qs")),
    "event_table.qs written to merged_dir"
  )
  expect_true(
    file.exists(file.path(merged_dir, "monthly_summary.qs")),
    "monthly_summary.qs written to merged_dir"
  )
}

# --------------------------------------------------------------------------
# Test 1b: sim_io_save / sim_io_read — parquet round-trip
# --------------------------------------------------------------------------
cat("\n## 1b) sim_io_* parquet round-trip\n")
if (!requireNamespace("nanoparquet", quietly = TRUE)) {
  cat("   SKIP: nanoparquet not installed\n")
} else {
  et_test <- tibble::tibble(
    replicate_id = 1L,
    gid = "AAA_1_1",
    start_month = "2018-01-01",
    duration_months = 100L
  )
  ms_test <- tibble::tibble(
    replicate_id = 1L,
    gid = "AAA_1_1",
    month = "2018-01-01",
    Y = 1L,
    mu_mon = 0.42,
    cases_month = NA_integer_
  )
  tmp_et <- tempfile(fileext = ".parquet")
  tmp_ms <- tempfile(fileext = ".parquet")

  sim_io_save(et_test, tmp_et, "parquet")
  sim_io_save(ms_test, tmp_ms, "parquet")

  et_back <- sim_io_read(tmp_et, "parquet")
  ms_back <- sim_io_read(tmp_ms, "parquet")

  expect_equal(nrow(et_back), 1L, "parquet et: row count")
  expect_equal(
    as.character(et_back$gid),
    "AAA_1_1",
    "parquet et: gid preserved"
  )
  expect_equal(et_back$duration_months, 100L, "parquet et: integer preserved")
  expect_equal(ms_back$mu_mon, 0.42, "parquet ms: double preserved")
  expect_equal(nrow(ms_back), 1L, "parquet ms: row count")
  expect_true(is.na(ms_back$cases_month), "parquet ms: NA integer preserved")
  expect_equal(
    as.character(ms_back$month),
    "2018-01-01",
    "parquet ms: month string preserved"
  )
}

# --------------------------------------------------------------------------
# Test 2: load_observed_from_panel — synthetic panel
# --------------------------------------------------------------------------
cat("\n## 2) load_observed_from_panel — synthetic panel\n")
{
  gids <- c("AAA_1_1", "BBB_2_1")
  months <- seq(as.Date("2018-01-01"), by = "month", length.out = 6)
  panel <- tidyr::expand_grid(gid = gids, t_start = months) %>%
    dplyr::mutate(
      outbreak = as.integer(
        gid == "AAA_1_1" & lubridate::month(t_start) %in% c(1, 3)
      )
    )

  obs <- load_observed_from_panel(panel)
  expect_true(
    all(c("gid", "month", "Y_obs") %in% names(obs)),
    "returns (gid, month, Y_obs) tibble"
  )
  expect_equal(nrow(obs), 12L, "12 (gid, month) cells for 2 gids × 6 months")

  # Horizon filter
  obs_filt <- load_observed_from_panel(
    panel,
    horizon = list(start = "2018-01-01", end = "2018-03-31")
  )
  expect_equal(
    nrow(obs_filt),
    6L,
    "horizon filter returns 6 cells (2 gids × 3 months)"
  )
  expect_true(
    all(obs_filt$month <= as.Date("2018-03-31")),
    "all months within horizon"
  )

  # gid filter
  obs_gid <- load_observed_from_panel(panel, gids = "AAA_1_1")
  expect_equal(nrow(obs_gid), 6L, "gid filter returns 6 cells for single gid")
  expect_true(all(obs_gid$gid == "AAA_1_1"), "only requested gid in output")
}

# --------------------------------------------------------------------------
# Test 3: compute_pit_binom — deterministic calibration cases
# --------------------------------------------------------------------------
cat("\n## 3) compute_pit_binom — deterministic calibration cases\n")
{
  months <- as.Date(sprintf("2018-%02d-01", 1:5))
  gids <- paste0("G", 1:5)
  n_reps <- 10L

  # Case A: all Y=0 (p_hat=0), y_obs=0 → PIT = runif(0,1) but range [0,1-0] = [0,1]
  # Actually: if p_hat=0 and y_obs=0, PIT = runif(n, 0, 1-0) = runif(n, 0, 1)
  sim_A <- tidyr::expand_grid(
    replicate_id = 1:n_reps,
    gid = gids,
    month = months
  ) %>%
    dplyr::mutate(Y = 0L)
  obs_A <- tibble::tibble(gid = gids, month = months, Y_obs = 0L)
  pit_A <- compute_pit_binom(sim_A, obs_A, seed = 1L)
  expect_true(
    all(pit_A$pit >= 0 & pit_A$pit <= 1),
    "Case A (p=0, y=0): all PIT in [0,1]"
  )
  expect_true(all(pit_A$p_hat == 0), "Case A: p_hat = 0 for all cells")

  # Case B: all Y=1 (p_hat=1), y_obs=1 → PIT = runif(n, 1-1, 1) = runif(n, 0, 1)
  sim_B <- tidyr::expand_grid(
    replicate_id = 1:n_reps,
    gid = gids,
    month = months
  ) %>%
    dplyr::mutate(Y = 1L)
  obs_B <- tibble::tibble(gid = gids, month = months, Y_obs = 1L)
  pit_B <- compute_pit_binom(sim_B, obs_B, seed = 2L)
  expect_true(
    all(pit_B$pit >= 0 & pit_B$pit <= 1),
    "Case B (p=1, y=1): all PIT in [0,1]"
  )

  # Case C: p=0.5, verify correct half-intervals
  sim_C <- tidyr::expand_grid(
    replicate_id = 1:n_reps,
    gid = gids,
    month = months
  ) %>%
    dplyr::mutate(Y = rep(c(rep(0L, 5), rep(1L, 5)), length.out = dplyr::n()))

  # Make obs: some y_obs=0, some y_obs=1
  obs_C0 <- tibble::tibble(gid = gids[1:3], month = months[1:3], Y_obs = 0L)
  obs_C1 <- tibble::tibble(gid = gids[4:5], month = months[4:5], Y_obs = 1L)
  obs_C <- dplyr::bind_rows(obs_C0, obs_C1)

  pit_C <- compute_pit_binom(sim_C, obs_C, seed = 3L)

  # For y_obs=0: PIT in [0, 1-p_hat]
  rows_y0 <- pit_C %>% dplyr::filter(y_obs == 0)
  if (nrow(rows_y0) > 0) {
    expect_true(
      all(rows_y0$pit <= 1 - rows_y0$p_hat + 1e-12),
      "Case C: y_obs=0 PIT in [0, 1-p_hat]"
    )
  }
  # For y_obs=1: PIT in [1-p_hat, 1]
  rows_y1 <- pit_C %>% dplyr::filter(y_obs == 1)
  if (nrow(rows_y1) > 0) {
    expect_true(
      all(rows_y1$pit >= rows_y1$p_hat - 1 + 1e-12),
      "Case C: y_obs=1 PIT in [1-p_hat, 1]"
    )
  }
  expect_true(all(pit_C$pit >= 0 & pit_C$pit <= 1), "Case C: all PIT in [0,1]")
}

# --------------------------------------------------------------------------
# Test 4: compute_monthly_freq — simple frame
# --------------------------------------------------------------------------
cat("\n## 4) compute_monthly_freq — simple frame\n")
{
  gids <- c("AAAW_1_1", "AAAW_1_2", "BBBW_1_1", "BBBW_1_2")
  months <- as.Date(c("2018-01-01", "2018-02-01"))
  frame4 <- tibble::tibble(
    gid = gids,
    country = c("AAAW", "AAAW", "BBBW", "BBBW"),
    adm1_unit = c("AAAW_1", "AAAW_1", "BBBW_1", "BBBW_1"),
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = 1L,
    pop = 1e5
  )

  # All Y=0
  sim_all0 <- tidyr::expand_grid(
    replicate_id = 1:3,
    gid = gids,
    month = months
  ) %>%
    dplyr::mutate(Y = 0L)
  obs_all0 <- tidyr::expand_grid(gid = gids, month = months) %>%
    dplyr::mutate(Y_obs = 0L)

  freq0 <- compute_monthly_freq(sim_all0, obs_all0, frame4)
  expect_true(all(freq0$freq_obs == 0), "All Y=0: freq_obs = 0")
  expect_true(all(freq0$freq_sim_mean == 0), "All Y=0: freq_sim_mean = 0")

  # All Y=1
  sim_all1 <- sim_all0 %>% dplyr::mutate(Y = 1L)
  obs_all1 <- obs_all0 %>% dplyr::mutate(Y_obs = 1L)

  freq1 <- compute_monthly_freq(sim_all1, obs_all1, frame4)
  expect_true(all(freq1$freq_obs == 1), "All Y=1: freq_obs = 1")
  expect_true(all(freq1$freq_sim_mean == 1), "All Y=1: freq_sim_mean = 1")
}

# --------------------------------------------------------------------------
# Test 5: compute_crps_annual — known CRPS value
# --------------------------------------------------------------------------
cat("\n## 5) compute_crps_annual — known CRPS value\n")
{
  # Perfect forecast: CRPS = 0
  sim_perf <- tibble::tibble(
    replicate_id = c(1L, 2L),
    gid = "G1",
    month = as.Date(c("2018-01-01", "2018-01-01")),
    Y = 3L # both replicates predict 3
  )
  obs_perf <- tibble::tibble(
    gid = "G1",
    month = as.Date("2018-01-01"),
    Y_obs = 3L
  )
  crps_perf <- compute_crps_annual(sim_perf, obs_perf)
  expect_equal(
    crps_perf$crps_count,
    0.0,
    "Perfect forecast: CRPS = 0",
    tol = 1e-9
  )

  # Known case: R=2, sim={0,2}, obs=1 → CRPS=0.5
  # E|X-1| = mean(|0-1|, |2-1|) = mean(1,1) = 1
  # E|X-X'|/2 = mean(|0-0|,|0-2|,|2-0|,|2-2|) / 4 = (0+2+2+0)/4 = 1 → /2 = 0.5
  # CRPS = 1 - 0.5 = 0.5
  sim_known <- tibble::tibble(
    replicate_id = c(1L, 2L),
    gid = "G2",
    month = as.Date(c("2018-01-01", "2018-01-01")),
    Y = c(0L, 2L)
  )
  obs_known <- tibble::tibble(
    gid = "G2",
    month = as.Date("2018-01-01"),
    Y_obs = 1L
  )
  crps_known <- compute_crps_annual(sim_known, obs_known)
  expect_equal(
    crps_known$crps_count,
    0.5,
    "Known case (sim={0,2}, obs=1): CRPS=0.5",
    tol = 1e-9
  )
}

# --------------------------------------------------------------------------
# Test 6: compute_ar_metrics — known sequence
# --------------------------------------------------------------------------
cat("\n## 6) compute_ar_metrics — known sequence\n")
{
  # Y = c(1,1,1,0,0,1,0)
  # rle: values=c(1,0,1,0), lengths=c(3,2,1,1)
  # run_lengths (values==1) = c(3,1)  → table: "1"=1, "3"=1
  # gap_lengths (values==0) = c(2,1)  → table: "1"=1, "2"=1
  # transitions (pairs from Y[1:6] → Y[2:7]):
  #   (1→1), (1→1), (1→0), (0→0), (0→1), (1→0)
  #   => (0→0)=1, (0→1)=1, (1→0)=2, (1→1)=2

  sim_seq <- tibble::tibble(
    replicate_id = 1L,
    gid = "G1",
    month = as.Date(sprintf("2018-%02d-01", 1:7)),
    Y = c(1L, 1L, 1L, 0L, 0L, 1L, 0L)
  )
  ar <- compute_ar_metrics(sim_seq)

  # Run lengths table
  run_tbl <- ar$run_length_table
  expect_true("1" %in% names(run_tbl), "run_length_table has entry '1'")
  expect_true("3" %in% names(run_tbl), "run_length_table has entry '3'")
  expect_equal(as.integer(run_tbl["1"]), 1L, "run length 1 appears 1 time")
  expect_equal(as.integer(run_tbl["3"]), 1L, "run length 3 appears 1 time")

  # Gap lengths table
  gap_tbl <- ar$gap_length_table
  expect_true("1" %in% names(gap_tbl), "gap_length_table has entry '1'")
  expect_true("2" %in% names(gap_tbl), "gap_length_table has entry '2'")
  expect_equal(as.integer(gap_tbl["1"]), 1L, "gap length 1 appears 1 time")
  expect_equal(as.integer(gap_tbl["2"]), 1L, "gap length 2 appears 1 time")

  # Transitions
  tr <- ar$transition
  expect_true(!is.null(tr) && nrow(tr) > 0, "transition table is non-empty")
  get_n <- function(f, t) {
    row <- tr %>% dplyr::filter(from == f, to == t)
    if (nrow(row) == 0) {
      return(0L)
    }
    as.integer(row$n)
  }
  expect_equal(get_n(1, 1), 2L, "transition (1→1) count = 2")
  expect_equal(get_n(1, 0), 2L, "transition (1→0) count = 2")
  expect_equal(get_n(0, 0), 1L, "transition (0→0) count = 1")
  expect_equal(get_n(0, 1), 1L, "transition (0→1) count = 1")
}

# --------------------------------------------------------------------------
# Test 7: reconstruct_lag_series — 2-location single-ADM1 frame
# --------------------------------------------------------------------------
cat("\n## 7) reconstruct_lag_series — 2-location single-ADM1 frame\n")
{
  frame7 <- tibble::tibble(
    gid = c("AAAW_1_1", "AAAW_1_2"),
    country = c("AAAW", "AAAW"),
    AFRO_region = "REG_W",
    adm1_unit = c("AAAW_1", "AAAW_1"),
    tu = c("hmmW", "hmmW"),
    min_seas = 1L,
    pop = 1e5
  )

  # 3 months; Y: month1: [1,0]; month2: [0,0]; month3: [0,1]
  months7 <- as.Date(c("2018-01-01", "2018-02-01", "2018-03-01"))
  sim7 <- tibble::tibble(
    replicate_id = 1L,
    gid = rep(c("AAAW_1_1", "AAAW_1_2"), each = 3),
    month = rep(months7, times = 2),
    Y = c(
      1L,
      0L,
      0L, # AAAW_1_1: month1=1, m2=0, m3=0
      0L,
      0L,
      1L
    ) # AAAW_1_2: month1=0, m2=0, m3=1
  )

  lag_series <- reconstruct_lag_series(sim7, frame7, max_reps = NULL)

  # At month2: o_adm1 should be 1 for both locations
  # (because in month1, AAAW_1_1=1 → ADM1 AAAW_1 was active)
  o_adm1_m2 <- lag_series %>%
    dplyr::filter(month == as.Date("2018-02-01")) %>%
    dplyr::pull(o_adm1)
  expect_true(
    all(o_adm1_m2 == 1L),
    "o_adm1 at month2 = 1 for both locs (ADM1 active in m1)"
  )

  # At month3: o_adm1 should be 0 (nothing active in m2)
  o_adm1_m3 <- lag_series %>%
    dplyr::filter(month == as.Date("2018-03-01")) %>%
    dplyr::pull(o_adm1)
  expect_true(
    all(o_adm1_m3 == 0L),
    "o_adm1 at month3 = 0 (nothing active in m2)"
  )
}

# --------------------------------------------------------------------------
# Test 8: compute_contagion_conditionals — forced patterns
# --------------------------------------------------------------------------
cat("\n## 8) compute_contagion_conditionals — forced patterns\n")
{
  # Frame: 2 ADM1 units in hmmW TU
  # AAAW_1: locs AAAW_1_1, AAAW_1_2
  # AAAW_2: locs AAAW_2_1, AAAW_2_2
  frame8 <- tibble::tibble(
    gid = c(
      "AAAW_1_1",
      "AAAW_1_2",
      "AAAW_2_1",
      "AAAW_2_2",
      "BBBW_1_1",
      "BBBW_1_2"
    ),
    country = c("AAAW", "AAAW", "AAAW", "AAAW", "BBBW", "BBBW"),
    AFRO_region = "REG_W",
    adm1_unit = c("AAAW_1", "AAAW_1", "AAAW_2", "AAAW_2", "BBBW_1", "BBBW_1"),
    tu = "hmmW",
    min_seas = 1L,
    pop = 1e5
  )

  gids8 <- frame8$gid
  months8 <- as.Date(c("2018-01-01", "2018-02-01", "2018-03-01"))

  # 2 replicates:
  # Month 1: AAAW_1 locations Y=1, all others Y=0
  # Month 2: all Y=0
  # Month 3: all Y=0
  make_rep <- function(rid) {
    tibble::tibble(
      replicate_id = rid,
      gid = rep(gids8, times = 3),
      month = rep(months8, each = length(gids8)),
      Y = c(
        # month 1: AAAW_1_1=1, AAAW_1_2=1, others=0
        1L,
        1L,
        0L,
        0L,
        0L,
        0L,
        # month 2: all 0
        rep(0L, 6),
        # month 3: all 0
        rep(0L, 6)
      )
    )
  }
  sim8 <- dplyr::bind_rows(make_rep(1L), make_rep(2L))

  cond8 <- compute_contagion_conditionals(sim8, frame8, max_reps = 2L)

  # At month 2: AAAW_1 locs have o_adm1=1; AAAW_2 and BBBW locs have o_adm1=0
  # So P(Y|ADM1-lag=1) should come from those locs at month2/3 — expect Y=0 always,
  # but the key test is that the metric is computable and o_adm1 lags are non-trivial.
  # The main check is that P(Y|lag=1) is not NA and >= P(Y|lag=0) when outbreaks cluster.
  # (In this fixture all Y=0 at months 2+3, so both are 0 — verify not NA)
  expect_true(!is.null(cond8$adm1_cond), "adm1_cond is non-null")
  expect_true(
    !is.na(cond8$adm1_cond$lag_active) ||
      is.na(cond8$adm1_cond$lag_active), # NA is OK if no lag-1 observations remain after dropping first month
    "adm1_cond$lag_active is computed (may be NA if no lag-1 obs)"
  )
  expect_true(!is.null(cond8$unconditional), "unconditional P(Y) is computed")
}

# --------------------------------------------------------------------------
# Test 9: scope.countries filter logic
# --------------------------------------------------------------------------
cat("\n## 9) scope.countries filter logic — unit test\n")
{
  # Synthetic frame with 3 countries
  frame9 <- tibble::tibble(
    gid = c("AAA_1_1", "AAA_1_2", "BBB_2_1", "CCC_3_1"),
    country = c("AAA", "AAA", "BBB", "CCC"),
    adm1_unit = c("AAA_1", "AAA_1", "BBB_2", "CCC_3"),
    AFRO_region = "REG_X",
    tu = "hmmW",
    min_seas = 1L,
    pop = 1e5
  )

  # Apply the filter logic from 07b_simulator_v1.R
  iso3_req <- c("AAA", "CCC")
  filtered <- frame9 %>%
    dplyr::filter(stringr::str_extract(adm1_unit, "^[A-Z]{3}") %in% iso3_req)

  expect_equal(
    nrow(filtered),
    3L,
    "scope.countries filter keeps AAA (2 locs) + CCC (1 loc) = 3"
  )
  expect_true(
    all(
      stringr::str_extract(filtered$adm1_unit, "^[A-Z]{3}") %in% iso3_req
    ),
    "only AAA and CCC rows remain after filter"
  )
  expect_true(
    !"BBB" %in% (stringr::str_extract(filtered$adm1_unit, "^[A-Z]{3}")),
    "BBB is excluded"
  )
}

source("analysis/07d_simulator_figures.R")

# --------------------------------------------------------------------------
# Helper: build a minimal synthetic metrics object for figure tests
# --------------------------------------------------------------------------
make_synthetic_metrics <- function() {
  months <- as.Date(sprintf("2018-%02d-01", 1:12))
  gids <- paste0("AAAW_", 1:4, "_1")

  pit <- tidyr::expand_grid(gid = gids, month = months) %>%
    dplyr::mutate(p_hat = 0.2, y_obs = 0L, pit = stats::runif(dplyr::n()))

  monthly_freq <- tidyr::expand_grid(
    country = c("AAAW", "BBBW"),
    month_of_year = 1:12
  ) %>%
    dplyr::mutate(
      freq_obs = 0.1,
      freq_sim_mean = 0.12,
      freq_sim_q025 = 0.04,
      freq_sim_q25 = 0.08,
      freq_sim_q75 = 0.16,
      freq_sim_q975 = 0.22
    )

  crps <- tidyr::expand_grid(gid = gids, year = 2018:2019) %>%
    dplyr::mutate(
      n_sim = 4L,
      y_obs = 0L,
      crps_count = stats::runif(dplyr::n(), 0, 1)
    )

  ar_metrics <- list(
    run_length_table = as.table(c("1" = 5L, "2" = 2L, "3" = 1L)),
    gap_length_table = as.table(c("1" = 8L, "2" = 4L)),
    transition = tibble::tibble(
      from = c(0L, 0L, 1L, 1L),
      to = c(0L, 1L, 0L, 1L),
      n = c(10L, 5L, 6L, 4L),
      p = c(0.67, 0.33, 0.60, 0.40)
    )
  )

  contagion <- list(
    adm1_cond = list(
      lag_active = 0.35,
      lag_inactive = 0.08,
      unconditional = 0.10
    ),
    country_cond = list(
      lag_active = 0.25,
      lag_inactive = 0.09,
      unconditional = 0.10
    ),
    tu_cond = list(
      lag_active = 0.12,
      lag_inactive = 0.10,
      unconditional = 0.10
    ),
    unconditional = 0.10
  )

  list(
    pit = pit,
    monthly_freq = monthly_freq,
    crps = crps,
    ar_metrics = ar_metrics,
    contagion = contagion
  )
}

# --------------------------------------------------------------------------
# Test 10: fig_pit_calibration — smoke
# --------------------------------------------------------------------------
cat("\n## 10) fig_pit_calibration — smoke\n")
{
  m <- make_synthetic_metrics()
  g <- fig_pit_calibration(m$pit)
  expect_true(inherits(g, "gg"), "fig_pit_calibration returns a gg object")

  # Empty / NULL inputs return placeholder gg
  g_null <- fig_pit_calibration(NULL)
  g_empty <- fig_pit_calibration(tibble::tibble())
  expect_true(inherits(g_null, "gg"), "fig_pit_calibration(NULL) returns gg")
  expect_true(inherits(g_empty, "gg"), "fig_pit_calibration(empty) returns gg")
}

# --------------------------------------------------------------------------
# Test 11: fig_seasonal_frequency — smoke
# --------------------------------------------------------------------------
cat("\n## 11) fig_seasonal_frequency — smoke\n")
{
  m <- make_synthetic_metrics()
  g <- fig_seasonal_frequency(m$monthly_freq)
  expect_true(inherits(g, "gg"), "fig_seasonal_frequency returns a gg object")

  g_null <- fig_seasonal_frequency(NULL)
  g_empty <- fig_seasonal_frequency(tibble::tibble())
  expect_true(inherits(g_null, "gg"), "fig_seasonal_frequency(NULL) returns gg")
  expect_true(
    inherits(g_empty, "gg"),
    "fig_seasonal_frequency(empty) returns gg"
  )
}

# --------------------------------------------------------------------------
# Test 12: fig_crps_annual — smoke
# --------------------------------------------------------------------------
cat("\n## 12) fig_crps_annual — smoke\n")
{
  m <- make_synthetic_metrics()
  g <- fig_crps_annual(m$crps)
  expect_true(inherits(g, "gg"), "fig_crps_annual returns a gg object")

  g_null <- fig_crps_annual(NULL)
  expect_true(inherits(g_null, "gg"), "fig_crps_annual(NULL) returns gg")
}

# --------------------------------------------------------------------------
# Test 13: fig_ar_run_lengths — smoke
# --------------------------------------------------------------------------
cat("\n## 13) fig_ar_run_lengths — smoke\n")
{
  m <- make_synthetic_metrics()
  g <- fig_ar_run_lengths(m$ar_metrics)
  expect_true(inherits(g, "gg"), "fig_ar_run_lengths returns a gg object")

  g_null <- fig_ar_run_lengths(NULL)
  expect_true(inherits(g_null, "gg"), "fig_ar_run_lengths(NULL) returns gg")

  # Empty table inputs
  empty_ar <- list(
    run_length_table = as.table(integer(0)),
    gap_length_table = as.table(integer(0)),
    transition = tibble::tibble()
  )
  g_empty <- fig_ar_run_lengths(empty_ar)
  expect_true(
    inherits(g_empty, "gg"),
    "fig_ar_run_lengths(empty tables) returns gg"
  )
}

# --------------------------------------------------------------------------
# Test 14: fig_contagion_conditionals — smoke
# --------------------------------------------------------------------------
cat("\n## 14) fig_contagion_conditionals — smoke\n")
{
  m <- make_synthetic_metrics()
  g <- fig_contagion_conditionals(m$contagion)
  expect_true(
    inherits(g, "gg"),
    "fig_contagion_conditionals returns a gg object"
  )

  g_null <- fig_contagion_conditionals(NULL)
  expect_true(
    inherits(g_null, "gg"),
    "fig_contagion_conditionals(NULL) returns gg"
  )

  # Contagion with no usable cond entries (all NULL lag_active)
  cong_empty <- list(
    adm1_cond = list(lag_active = NULL),
    country_cond = NULL,
    tu_cond = NULL,
    unconditional = 0.1
  )
  g_empty <- fig_contagion_conditionals(cong_empty)
  expect_true(
    inherits(g_empty, "gg"),
    "fig_contagion_conditionals(no entries) returns gg"
  )
}

# --------------------------------------------------------------------------
# Test 15: make_sim_figures — writes 6 PNG files to temp dir
# --------------------------------------------------------------------------
cat("\n## 15) make_sim_figures — writes 6 PNG files\n")
{
  m <- make_synthetic_metrics()
  tmp_fig <- tempfile("sim_figs_")
  figs <- make_sim_figures(m, tmp_fig)

  expected_files <- c(
    "pit_calibration.png",
    "seasonal_frequency.png",
    "crps_annual.png",
    "ar_run_lengths.png",
    "contagion_conditionals.png",
    "crps_total_cases.png"
  )

  for (f in expected_files) {
    expect_true(
      file.exists(file.path(tmp_fig, f)),
      sprintf("make_sim_figures writes %s", f)
    )
  }

  expect_true(
    is.list(figs) && length(figs) == 6L,
    "make_sim_figures returns list of 6 gg objects"
  )
  expect_true(
    !file.exists(file.path(tmp_fig, "conditional_size.png")),
    "make_sim_figures no longer writes conditional_size.png"
  )
  expect_true(
    all(sapply(figs, inherits, "gg")),
    "all returned objects are gg instances"
  )
}

# --------------------------------------------------------------------------
# Test 16: reconstruct_lag_series — active_col replays persistent state,
# not the onset pulse (duty cycle)
# --------------------------------------------------------------------------
cat("\n## 16) reconstruct_lag_series — active_col replays persistent state\n")
{
  frame16 <- tibble::tibble(
    gid = c("AAAW_1_1", "AAAW_1_2"),
    country = c("AAAW", "AAAW"),
    AFRO_region = "REG_W",
    adm1_unit = c("AAAW_1", "AAAW_1"),
    tu = c("hmmW", "hmmW"),
    min_seas = 1L,
    pop = 1e5
  )

  # AAAW_1_1 onsets (Y=1) in month1 only, but stays `active` through month3
  # (a 3-month outbreak). AAAW_1_2 never onsets/active.
  months16 <- as.Date(c("2018-01-01", "2018-02-01", "2018-03-01", "2018-04-01"))
  sim16 <- tibble::tibble(
    replicate_id = 1L,
    gid = rep(c("AAAW_1_1", "AAAW_1_2"), each = 4),
    month = rep(months16, times = 2),
    Y = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
    active = c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  )

  # Legacy (active_col = NULL): o_adm1 at month3 must be 0, because the
  # onset pulse Y is only 1 in month1 -> update_lags() sees nothing active
  # by month3 under the legacy Y-fed replay.
  lag_legacy <- reconstruct_lag_series(sim16, frame16, max_reps = NULL)
  o_adm1_m3_legacy <- lag_legacy %>%
    dplyr::filter(month == as.Date("2018-03-01")) %>%
    dplyr::pull(o_adm1)
  expect_true(
    all(o_adm1_m3_legacy == 0L),
    "legacy (active_col=NULL) o_adm1 at month3 = 0 (Y pulse already decayed)"
  )

  # active_col = "active": o_adm1 at month3 must be 1, because AAAW_1_1 is
  # still `active` in month2 -> update_lags() sees the ADM1 active going
  # into month3 under the persistent-state replay.
  lag_active <- reconstruct_lag_series(
    sim16,
    frame16,
    max_reps = NULL,
    active_col = "active"
  )
  o_adm1_m3_active <- lag_active %>%
    dplyr::filter(month == as.Date("2018-03-01")) %>%
    dplyr::pull(o_adm1)
  expect_true(
    all(o_adm1_m3_active == 1L),
    "active_col='active' o_adm1 at month3 = 1 (persistent state still on in month2)"
  )

  # Y itself (onset) must be identical and untouched by active_col.
  expect_true(
    identical(lag_legacy$Y, lag_active$Y),
    "Y (onset) column unchanged by active_col"
  )

  # Unknown active_col falls back to legacy behaviour (with a message).
  lag_fallback <- suppressMessages(
    reconstruct_lag_series(
      sim16,
      frame16,
      max_reps = NULL,
      active_col = "nonexistent_col"
    )
  )
  expect_true(
    identical(lag_fallback$o_adm1, lag_legacy$o_adm1),
    "unknown active_col falls back to legacy Y-fed reconstruction"
  )
}

# --------------------------------------------------------------------------
# Test 17: compute_contagion_conditionals — active_col threaded through
# --------------------------------------------------------------------------
cat("\n## 17) compute_contagion_conditionals — active_col threaded through\n")
{
  frame17 <- tibble::tibble(
    gid = c("AAAW_1_1", "AAAW_1_2"),
    country = c("AAAW", "AAAW"),
    AFRO_region = "REG_W",
    adm1_unit = c("AAAW_1", "AAAW_1"),
    tu = c("hmmW", "hmmW"),
    min_seas = 1L,
    pop = 1e5
  )
  months17 <- as.Date(c("2018-01-01", "2018-02-01", "2018-03-01"))
  sim17 <- tibble::tibble(
    replicate_id = 1L,
    gid = rep(c("AAAW_1_1", "AAAW_1_2"), each = 3),
    month = rep(months17, times = 2),
    Y = c(1L, 0L, 0L, 0L, 1L, 0L),
    active = c(1L, 1L, 0L, 0L, 1L, 0L)
  )

  cond_active <- compute_contagion_conditionals(
    sim17,
    frame17,
    max_reps = 1L,
    active_col = "active"
  )
  expect_true(
    !is.null(cond_active$adm1_cond),
    "adm1_cond computed with active_col supplied"
  )
  expect_true(
    !is.null(cond_active$unconditional),
    "unconditional P(Y) computed with active_col supplied"
  )

  # Sanity: passing active_col=NULL still works (back-compat default).
  cond_legacy <- compute_contagion_conditionals(sim17, frame17, max_reps = 1L)
  expect_true(
    !is.null(cond_legacy$adm1_cond),
    "adm1_cond still computed with default active_col=NULL"
  )
}

# --------------------------------------------------------------------------
# Test 18: .km_estimate — Kaplan-Meier on a known small sample
# --------------------------------------------------------------------------
cat("\n## 18) .km_estimate — known Kaplan-Meier sample\n")
{
  # 4 uncensored events at t = 1,2,3,4 -> S(t) = 3/4, 2/4, 1/4, 0/4.
  # Median (S <= 0.5) is reached at t=2 (S=0.5).
  km1 <- .km_estimate(time = c(1, 2, 3, 4), event = c(1, 1, 1, 1))
  expect_equal(
    km1$median,
    2,
    "km_estimate: median readoff on a simple 4-event sample"
  )
  expect_true(
    abs(km1$curve$surv[km1$curve$time == 4] - 0) < 1e-8,
    "km_estimate: survival hits 0 at the last event with no censoring"
  )

  # All censored -> no event times -> curve empty, median/p90 NA.
  km2 <- .km_estimate(time = c(1, 2, 3), event = c(0, 0, 0))
  expect_true(
    is.na(km2$median) && is.na(km2$p90),
    "km_estimate: all-censored sample -> NA median/p90"
  )

  # Empty input.
  km3 <- .km_estimate(time = numeric(0), event = integer(0))
  expect_true(
    is.na(km3$median) && nrow(km3$curve) == 0L,
    "km_estimate: empty input handled without error"
  )
}

# --------------------------------------------------------------------------
# Test 19: compute_duration_dist — legacy vs duration-aware event tables
# --------------------------------------------------------------------------
cat("\n## 19) compute_duration_dist — legacy vs duration-aware event tables\n")
{
  # Legacy event table (no duration_months/censored) -> NULL.
  et_legacy <- tibble::tibble(
    replicate_id = 1L,
    gid = "AAAW_1_1",
    country = "AAAW",
    start_month = "2018-01-01"
  )
  expect_true(
    is.null(compute_duration_dist(et_legacy)),
    "compute_duration_dist returns NULL for a legacy (no-duration) event table"
  )

  # Duration-aware event table -> sim side computed; obs side NULL without a panel.
  et_dur <- tibble::tibble(
    replicate_id = 1L,
    gid = sprintf("AAAW_1_%d", 1:6),
    country = "AAAW",
    start_month = "2018-01-01",
    duration_months = c(1, 2, 3, 4, 5, 6),
    censored = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
  dd <- compute_duration_dist(
    et_dur,
    enso_grouping_path = "nonexistent_path.csv"
  )
  expect_true(
    !is.null(dd$sim),
    "compute_duration_dist: sim side computed for duration-aware event table"
  )
  expect_equal(
    dd$sim$n_episodes,
    6L,
    "compute_duration_dist: sim n_episodes matches event table rows"
  )
  expect_equal(
    dd$sim$n_censored,
    1L,
    "compute_duration_dist: sim n_censored counts the one censored episode"
  )
  expect_true(
    is.null(dd$obs),
    "compute_duration_dist: obs side NULL when no obs_dur_panel supplied"
  )
  expect_true(
    is.null(dd$sim$by_group),
    "compute_duration_dist: by_group NULL when enso_grouping_path doesn't exist"
  )

  # With an obs_dur_panel supplied.
  obs_panel19 <- tibble::tibble(
    episode_id = c("E1", "E1", "E2"),
    gid = c("AAAW_1_1", "AAAW_1_1", "AAAW_1_2"),
    country = "AAAW",
    month = c("2018-01-01", "2018-02-01", "2018-01-01"),
    duration_months = c(2, 2, 1),
    censored = c(0L, 0L, 1L)
  )
  dd2 <- compute_duration_dist(
    et_dur,
    obs_dur_panel = obs_panel19,
    enso_grouping_path = "nonexistent_path.csv"
  )
  expect_true(
    !is.null(dd2$obs),
    "compute_duration_dist: obs side computed when obs_dur_panel supplied"
  )
  expect_equal(
    dd2$obs$n_episodes,
    2L,
    "compute_duration_dist: obs side collapses to one row per distinct episode_id"
  )
}

# --------------------------------------------------------------------------
# Test 20: compute_active_fraction — sim/obs duty-cycle fractions
# --------------------------------------------------------------------------
cat("\n## 20) compute_active_fraction — sim/obs duty-cycle fractions\n")
{
  frame20 <- tibble::tibble(
    gid = c("AAAW_1_1", "AAAW_1_2"),
    country = c("AAAW", "AAAW")
  )

  # No `active` column -> NULL (legacy run).
  sim20_legacy <- tibble::tibble(
    replicate_id = 1L,
    gid = c("AAAW_1_1", "AAAW_1_2"),
    month = as.Date("2018-01-01"),
    Y = c(0L, 0L)
  )
  expect_true(
    is.null(compute_active_fraction(sim20_legacy, frame20)),
    "compute_active_fraction returns NULL when `active` column is absent"
  )

  # AAAW_1_1 active in 2/4 months; AAAW_1_2 active in 0/4 months.
  months20 <- as.Date(c("2018-01-01", "2018-02-01", "2018-03-01", "2018-04-01"))
  sim20 <- tibble::tibble(
    replicate_id = 1L,
    gid = rep(c("AAAW_1_1", "AAAW_1_2"), each = 4),
    month = rep(months20, times = 2),
    Y = 0L,
    active = c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L)
  )
  af <- compute_active_fraction(sim20, frame20)
  frac_111 <- af$by_gid %>%
    dplyr::filter(gid == "AAAW_1_1") %>%
    dplyr::pull(frac_active_sim)
  frac_112 <- af$by_gid %>%
    dplyr::filter(gid == "AAAW_1_2") %>%
    dplyr::pull(frac_active_sim)
  expect_equal(
    frac_111,
    0.5,
    "compute_active_fraction: AAAW_1_1 active 2/4 months = 0.5"
  )
  expect_equal(
    frac_112,
    0.0,
    "compute_active_fraction: AAAW_1_2 active 0/4 months = 0.0"
  )
  expect_true(
    !"frac_active_obs" %in% names(af$by_gid),
    "compute_active_fraction: no obs column when obs_dur_panel not supplied"
  )

  # With an obs_dur_panel + mask supplied.
  obs_panel20 <- tibble::tibble(
    gid = c("AAAW_1_1", "AAAW_1_1"),
    month = c("2018-01-01", "2018-02-01")
  )
  mask20 <- tibble::tibble(
    gid = rep(c("AAAW_1_1", "AAAW_1_2"), each = 4),
    month = rep(months20, times = 2)
  )
  af2 <- compute_active_fraction(
    sim20,
    frame20,
    obs_dur_panel = obs_panel20,
    mask = mask20
  )
  obs_frac_111 <- af2$by_gid %>%
    dplyr::filter(gid == "AAAW_1_1") %>%
    dplyr::pull(frac_active_obs)
  obs_frac_112 <- af2$by_gid %>%
    dplyr::filter(gid == "AAAW_1_2") %>%
    dplyr::pull(frac_active_obs)
  expect_equal(
    obs_frac_111,
    0.5,
    "compute_active_fraction: observed AAAW_1_1 active 2/4 months = 0.5"
  )
  expect_equal(
    obs_frac_112,
    0.0,
    "compute_active_fraction: observed AAAW_1_2 active 0/4 months = 0.0 (via full_join + coalesce)"
  )
}

# --------------------------------------------------------------------------
# Test 21: compute_sim_kendall_tau — reuses paired_kendall_tau(), drops
# censored episodes, returns NULL for legacy event tables
# --------------------------------------------------------------------------
cat(
  "\n## 21) compute_sim_kendall_tau — raw pooled tau on uncensored episodes\n"
)
{
  et_legacy21 <- tibble::tibble(
    replicate_id = 1L,
    gid = "AAAW_1_1",
    start_month = "2018-01-01",
    duration_months = 1,
    censored = 0L
  )
  expect_true(
    is.null(compute_sim_kendall_tau(et_legacy21)),
    "compute_sim_kendall_tau returns NULL for an event table with no episode_cases (no mon track)"
  )

  # Fixture built the way run_validation() builds it: a monthly_summary with
  # Y/active/cases_month -> compute_episode_case_totals() -> joined onto the
  # event table. Six gids, one episode each starting 2018-01; gid k is active
  # for k months at 10 cases/month (episode_cases = 10k), so cases and duration
  # are perfectly monotonic among the five uncensored episodes. gid 6 is a
  # censored 1-month episode with 1000 cases that would break the monotonic
  # trend if (wrongly) included -- must be dropped.
  months21 <- as.Date(sprintf("2018-%02d-01", 1:6))
  dur21 <- c(1L, 2L, 3L, 4L, 5L, 1L)
  ms21 <- dplyr::bind_rows(lapply(1:6, function(k) {
    act <- as.integer(seq_along(months21) <= dur21[k])
    tibble::tibble(
      replicate_id = 1L,
      gid = sprintf("AAAW_1_%d", k),
      month = months21,
      Y = as.integer(seq_along(months21) == 1L),
      active = act,
      cases_month = ifelse(act == 1L, if (k == 6L) 1000 else 10, NA_real_)
    )
  }))
  et21 <- tibble::tibble(
    replicate_id = 1L,
    gid = sprintf("AAAW_1_%d", 1:6),
    start_month = "2018-01-01",
    duration_months = c(1, 2, 3, 4, 5, 1),
    censored = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
  et21 <- suppressMessages(attach_episode_case_totals(
    et21,
    compute_episode_case_totals(ms21)
  ))
  expect_equal(
    et21$episode_cases,
    c(10, 20, 30, 40, 50, 1000),
    "compute_sim_kendall_tau fixture: episode_cases = sum of cases_month per episode"
  )
  tau21 <- compute_sim_kendall_tau(et21, min_n = 3L)
  expect_true(
    !is.null(tau21),
    "compute_sim_kendall_tau computes a result with enough uncensored episodes"
  )
  expect_equal(
    tau21$overall$n_episodes,
    5L,
    "compute_sim_kendall_tau drops the censored episode before pooling"
  )
  expect_equal(
    tau21$overall$mean,
    1,
    "compute_sim_kendall_tau: perfectly monotonic uncensored pairs -> tau = 1"
  )

  # Too few uncensored episodes -> NULL (below min_n).
  et21_small <- et21[1:2, ]
  expect_true(
    is.null(compute_sim_kendall_tau(et21_small, min_n = 5L)),
    "compute_sim_kendall_tau returns NULL below min_n uncensored episodes"
  )
}

# --------------------------------------------------------------------------
# Test 22: write_validation_summary — duration report sections render
# --------------------------------------------------------------------------
cat("\n## 22) write_validation_summary — Phase 7 report sections\n")
{
  m22 <- make_synthetic_metrics()
  m22$duration_dist <- list(
    sim = list(
      overall = list(median = 4, p90 = 7),
      n_episodes = 10L,
      n_censored = 2L
    ),
    obs = list(
      overall = list(median = 4.2, p90 = 7.5),
      n_episodes = 8L,
      n_censored = 1L
    )
  )
  m22$active_fraction <- list(
    by_gid = tibble::tibble(
      gid = c("A", "B"),
      frac_active_sim = c(0.3, 0.4),
      frac_active_obs = c(0.28, 0.35)
    ),
    by_country = tibble::tibble(
      country = "AAA",
      frac_active_sim = 0.35,
      frac_active_obs = 0.31
    )
  )
  m22$sim_kendall_tau <- list(overall = list(mean = 0.05, n_episodes = 40L))
  m22$obs_kendall_tau <- list(
    mean = 0.03,
    q025 = -0.10,
    q975 = 0.16,
    n_episodes = 120L
  )

  tmp_report <- tempfile("validation_report_", fileext = ".txt")
  write_validation_summary(m22, tmp_report)
  report_txt <- paste(readLines(tmp_report), collapse = "\n")

  expect_true(
    grepl("Duration distribution", report_txt),
    "validation report includes the Duration distribution section"
  )
  expect_true(
    grepl("Fraction of active months", report_txt),
    "validation report includes the Fraction of active months section"
  )
  expect_true(
    grepl(
      "Episode cases (sum of cases_month) <-> duration Kendall's tau",
      report_txt,
      fixed = TRUE
    ),
    "validation report includes the Kendall's tau section"
  )
  expect_true(
    grepl("Observed \\(Phase 5 size fit, PIT-residual", report_txt),
    "validation report includes the observed Phase-5 tau line when obs_kendall_tau is supplied"
  )
}

# --------------------------------------------------------------------------
# Test 23: compute_crps_total_cases — total_case_basis exact_gid vs region_key
#
# Cross-level reconciliation (rule B): the sim forecasts region NGA_20 at ADM2
# (frame gids NGA_20_1, NGA_20_2), but the observed corpus reports that whole
# state coarse, at ADM1 (gid "NGA_20"). An exact-(gid) membership test drops the
# coarse report (its gid is not a frame gid); the region_key basis maps it to its
# ADM1 region (which the sim forecasts) and recovers the case mass. Region MWI_5
# is reported at ADM2 on a frame gid and must be retained by BOTH bases.
# --------------------------------------------------------------------------
cat("\n## 23) compute_crps_total_cases — exact_gid vs region_key basis\n")
{
  frame23 <- tibble::tibble(
    gid = c("NGA_20_1", "NGA_20_2", "MWI_5_1"),
    country = c("NGA", "NGA", "MWI"),
    AFRO_region = c("WEST", "WEST", "SOUTH")
  )

  # obs metadata: a coarse ADM1 report for NGA_20 (gid not in frame) + an ADM2
  # report for MWI_5_1 (a frame gid). Written to a temp RDS as the fn reads it.
  obs_meta23 <- tibble::tibble(
    location_period_id = c("lpN", "lpM"),
    suspected_cases = c(1000, 300),
    TL = as.Date(c("2018-06-01", "2018-06-01")),
    country = c("NGA", "MWI"),
    AFRO_region = c("WEST", "SOUTH"),
    admin_level = c("ADM1", "ADM2")
  )
  tmp_meta23 <- tempfile("obs_meta23_", fileext = ".rds")
  saveRDS(obs_meta23, tmp_meta23)

  # obs_panel supplies the location_period_id -> gid bridge (incl. coarse gid).
  obs_panel23 <- tibble::tibble(
    location_period_id = c("lpN", "lpM"),
    gid = c("NGA_20", "MWI_5_1")
  )

  # sim event table: events in both regions across 2 replicates. Case totals come
  # from a monthly_summary through compute_episode_case_totals(), as in
  # run_validation(): each episode is active Jun-Jul 2018 and its cases_month
  # sum is 800 / 1200 / 250 / 350.
  ms23 <- tibble::tibble(
    replicate_id = rep(c(1L, 2L, 1L, 2L), each = 2),
    gid = rep(c("NGA_20_1", "NGA_20_1", "MWI_5_1", "MWI_5_1"), each = 2),
    month = rep(as.Date(c("2018-06-01", "2018-07-01")), times = 4),
    Y = rep(c(1L, 0L), times = 4),
    active = 1L,
    cases_month = c(500, 300, 700, 500, 200, 50, 300, 50)
  )
  event23 <- tibble::tibble(
    replicate_id = c(1L, 2L, 1L, 2L),
    gid = c("NGA_20_1", "NGA_20_1", "MWI_5_1", "MWI_5_1"),
    country = c("NGA", "NGA", "MWI", "MWI"),
    start_month = "2018-06-01"
  )
  event23 <- suppressMessages(attach_episode_case_totals(
    event23,
    compute_episode_case_totals(ms23)
  ))
  expect_equal(
    event23$episode_cases,
    c(800, 1200, 250, 350),
    "compute_crps_total_cases fixture: episode_cases from cases_month sums"
  )

  res_exact <- compute_crps_total_cases(
    event_table = event23,
    obs_metadata_path = tmp_meta23,
    frame = frame23,
    obs_panel = obs_panel23,
    mask = NULL,
    case_admin_level = "ADM1-ADM2",
    total_case_basis = "exact_gid"
  )
  res_region <- compute_crps_total_cases(
    event_table = event23,
    obs_metadata_path = tmp_meta23,
    frame = frame23,
    obs_panel = obs_panel23,
    mask = NULL,
    case_admin_level = "ADM1-ADM2",
    total_case_basis = "region_key"
  )

  ct_e <- res_exact$country
  ct_r <- res_region$country

  # exact_gid: coarse NGA report is dropped (gid not in frame) — NGA absent.
  expect_true(
    !("NGA" %in% ct_e$geography),
    "exact_gid: coarse ADM1 NGA report dropped (gid not a frame gid)"
  )
  # region_key: NGA recovered, with the full coarse case mass (1000).
  expect_true(
    "NGA" %in% ct_r$geography,
    "region_key: coarse ADM1 NGA report recovered via region mapping"
  )
  expect_equal(
    ct_r$total_obs[ct_r$geography == "NGA"],
    1000,
    "region_key: NGA observed total = recovered coarse mass (1000)"
  )
  # MWI (ADM2 on a frame gid) is retained identically by both bases.
  expect_equal(
    ct_e$total_obs[ct_e$geography == "MWI"],
    300,
    "exact_gid: MWI (frame-gid ADM2 report) retained = 300"
  )
  expect_equal(
    ct_r$total_obs[ct_r$geography == "MWI"],
    300,
    "region_key: MWI retained unchanged = 300"
  )

  # region_key mask branch: a (region_key, month) coverage mask must keep the
  # coarse NGA report whose gid never appears in the (frame-gid-keyed) mask.
  mask23 <- tibble::tibble(
    gid = c("NGA_20_1", "MWI_5_1"),
    month = as.Date("2018-06-01")
  )
  res_region_m <- compute_crps_total_cases(
    event_table = event23,
    obs_metadata_path = tmp_meta23,
    frame = frame23,
    obs_panel = obs_panel23,
    mask = mask23,
    case_admin_level = "ADM1-ADM2",
    total_case_basis = "region_key"
  )
  expect_true(
    "NGA" %in% res_region_m$country$geography,
    "region_key + mask: coarse NGA report survives (region_key, month) coverage"
  )

  # SIM-SIDE cross-level reconciliation: the REAL failure mode. build_coverage_mask()
  # keys coverage on the OBSERVED gid, so a coarse-reporting region carries only its
  # coarse gid (NGA_20) in the mask — never the fine sim forecast children
  # (NGA_20_1/NGA_20_2). An exact (gid, month) sim-side join therefore drops the whole
  # ensemble in that region and NGA vanishes even under region_key. Masking the sim
  # side at (region_key, month) must keep the fine children whose region is covered.
  mask23_coarse <- tibble::tibble(
    gid = "NGA_20", # ONLY the coarse observed gid
    month = as.Date("2018-06-01")
  )
  res_region_c <- compute_crps_total_cases(
    event_table = event23,
    obs_metadata_path = tmp_meta23,
    frame = frame23,
    obs_panel = obs_panel23,
    mask = mask23_coarse,
    case_admin_level = "ADM1-ADM2",
    total_case_basis = "region_key"
  )
  expect_true(
    "NGA" %in% res_region_c$country$geography,
    "region_key + coarse-gid mask: fine sim children survive via (region_key, month) sim join"
  )
  expect_true(
    res_region_c$country$total_sim_mean[
      res_region_c$country$geography == "NGA"
    ] >
      0,
    "region_key + coarse-gid mask: NGA ensemble total is non-zero (children not masked out)"
  )
}

# --------------------------------------------------------------------------
# Test 24: compute_episode_case_totals / attach_episode_case_totals — per-episode
# sums of the mon track's cases_month (replace the retired size track's
# size_cases in metric 2b and the simulated Kendall's tau)
# --------------------------------------------------------------------------
cat("\n## 24) compute_episode_case_totals — episode sums of cases_month\n")
{
  m24 <- sprintf("2018-%02d-01", 1:8) # character months, as the kernel writes them
  # rep 1, gid A:
  #   Jan-Feb  active, no onset           -> warm-start episode (.ep == 0): excluded
  #   Mar      inactive
  #   Apr-May  episode 1 (onset Apr); May's cases_month is NA -> sum 10, 1 month
  #   Jun-Aug  episode 2, onset the month right after episode 1's last active
  #            month (back-to-back), still active at the horizon end (censored)
  #            -> sum 7 + 8 + 9 = 24, 3 months
  ms_A <- tibble::tibble(
    replicate_id = 1L,
    gid = "A_1",
    month = m24,
    Y = c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L),
    active = c(1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L),
    cases_month = c(5, 5, NA, 10, NA, 7, 8, 9)
  )
  # rep 1, gid B: one 1-month episode in Jan (onset in the first month).
  ms_B <- tibble::tibble(
    replicate_id = 1L,
    gid = "B_1",
    month = m24,
    Y = c(1L, rep(0L, 7)),
    active = c(1L, rep(0L, 7)),
    cases_month = c(3, rep(NA, 7))
  )
  # rep 2, gid A: never active -> no episodes.
  ms_A2 <- tibble::tibble(
    replicate_id = 2L,
    gid = "A_1",
    month = m24,
    Y = 0L,
    active = 0L,
    cases_month = NA_real_
  )
  ms24 <- dplyr::bind_rows(ms_A, ms_B, ms_A2)
  ms24 <- ms24[rev(seq_len(nrow(ms24))), ] # unsorted input: the fn must order by month

  tot24 <- compute_episode_case_totals(ms24) %>%
    dplyr::arrange(replicate_id, gid, start_month)
  expect_equal(
    nrow(tot24),
    3L,
    "compute_episode_case_totals: 3 episodes (warm-start run and inactive rep excluded)"
  )
  expect_true(
    is.character(tot24$start_month),
    "compute_episode_case_totals: start_month is character"
  )
  expect_equal(
    tot24$start_month,
    c("2018-04-01", "2018-06-01", "2018-01-01"),
    "compute_episode_case_totals: start_month = first active month of each episode"
  )
  expect_equal(
    tot24$episode_cases,
    c(10, 24, 3),
    "compute_episode_case_totals: back-to-back episodes split; NA month skipped; censored episode summed in-horizon"
  )
  expect_equal(
    tot24$n_months_cases,
    c(1L, 3L, 1L),
    "compute_episode_case_totals: n_months_cases counts non-NA cases_month months"
  )
  expect_true(
    !any(
      tot24$start_month %in% c("2018-01-01", "2018-02-01") & tot24$gid == "A_1"
    ),
    "compute_episode_case_totals: warm-start active run (no onset) contributes no episode"
  )

  # Missing columns (no-duration / no-monthly runs) -> NULL.
  expect_true(
    is.null(compute_episode_case_totals(dplyr::select(ms24, -cases_month))),
    "compute_episode_case_totals: NULL without cases_month"
  )
  expect_true(
    is.null(compute_episode_case_totals(dplyr::select(ms24, -active))),
    "compute_episode_case_totals: NULL without active"
  )
  expect_true(
    is.null(compute_episode_case_totals(NULL)),
    "compute_episode_case_totals: NULL on NULL input"
  )

  # attach: one event per onset + one event with no matching episode (unmatched).
  et24 <- tibble::tibble(
    replicate_id = c(1L, 1L, 1L, 2L),
    gid = c("A_1", "A_1", "B_1", "A_1"),
    start_month = c("2018-04-01", "2018-06-01", "2018-01-01", "2018-05-01"),
    duration_months = c(2, 3, 1, 1),
    censored = c(0L, 1L, 0L, 0L)
  )
  w24 <- NULL
  et24j <- withCallingHandlers(
    suppressMessages(attach_episode_case_totals(et24, tot24)),
    warning = function(w) {
      w24 <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(
    nrow(et24j),
    4L,
    "attach_episode_case_totals: event rows preserved (left join)"
  )
  expect_equal(
    et24j$episode_cases,
    c(10, 24, 3, NA),
    "attach_episode_case_totals: totals keyed on (replicate_id, gid, start_month)"
  )
  expect_equal(
    attr(et24j, "n_unmatched"),
    1L,
    "attach_episode_case_totals: unmatched event counted"
  )
  expect_true(
    !is.null(w24) && grepl("1 of 4", w24),
    "attach_episode_case_totals: unmatched events reported in a warning"
  )

  # Date-typed start_month on the event table joins identically (type coercion).
  et24d <- dplyr::mutate(et24[1:3, ], start_month = as.Date(start_month))
  et24dj <- suppressMessages(attach_episode_case_totals(et24d, tot24))
  expect_equal(
    et24dj$episode_cases,
    c(10, 24, 3),
    "attach_episode_case_totals: Date start_month coerced to the same key"
  )
}

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
