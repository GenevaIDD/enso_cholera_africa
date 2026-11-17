# test_population_gates.R — HARD: the population-denominator QC gates 1-5b
#   (analysis/utils/population_qc.R).
#
# Usage:
#   Rscript analysis/tests/test_population_gates.R
#
# Details:
#   Synthetic fixtures only -- no real extraction data is read, so the test
#   runs identically on a corpus with or without the provenance columns.

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

source("analysis/utils/population_qc.R")

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# population QC gate tests\n\n")

# --------------------------------------------------------------------------
# Shared fixture: 4 location periods exercising pop==0 (gate 2), pop < floor
# but positive (the known "near-zero raster" failure mode), an all-NA pop
# LP with no matching geometry (gate 1 / flag_no_geom), and one clean LP.
# --------------------------------------------------------------------------
make_ob_shp <- function() {
  ob <- tibble::tibble(
    location = c(
      "AFR::AAA::Zero",
      "AFR::BBB::Low",
      "AFR::CCC::NoGeom",
      "AFR::DDD::Good"
    ),
    location_period_id = c("lp1", "lp2", "lp3", "lp4"),
    country = c("AAA", "BBB", "CCC", "DDD"),
    outbreak_number = c(1, 1, 1, 1),
    TL = as.Date("2015-01-01"),
    TR = as.Date("2015-02-01"),
    sCh = c(100, 50, 10, 20),
    pop = c(0, 0.5, NA, 10000),
    threshold = 1,
    adj_factors = c(1.6, 0.4, 1.0, 1.0)
  )
  shp <- st_as_sf(tibble::tibble(
    lctn_pr = c("lp1", "lp2", "lp3", "lp4"),
    area_per_1km2 = c(10, 10, NA, 50),
    geometry = st_sfc(
      st_point(c(0, 0)),
      st_point(c(1, 1)),
      st_point(c(2, 2)),
      st_point(c(3, 3))
    )
  ))
  list(ob = ob, shp = shp)
}

cat("## A) check_population_integrity: gates 1-2\n")
{
  fx <- make_ob_shp()
  lp_pop <- check_population_integrity(fx$ob, fx$shp)
  by_lp <- function(id) lp_pop %>% filter(location_period_id == id)

  expect_true(
    isTRUE(by_lp("lp1")$flag_low_pop),
    "pop == 0 flagged low_pop (below the person-count floor)"
  )
  expect_true(
    isTRUE(by_lp("lp2")$flag_low_pop),
    "pop == 0.5 (near-zero raster failure mode) flagged low_pop"
  )
  expect_true(
    !isTRUE(by_lp("lp4")$flag_low_pop),
    "pop == 10000 not flagged low_pop"
  )
  expect_true(
    isTRUE(by_lp("lp3")$flag_na_pop),
    "all-NA pop flagged na_pop (gate 1: must be NA, never left unflagged)"
  )
  expect_true(
    isTRUE(by_lp("lp3")$flag_no_geom),
    "LP with no matching geometry flagged flag_no_geom"
  )
  expect_true(
    !isTRUE(by_lp("lp1")$flag_no_geom),
    "LP with a matching geometry not flagged flag_no_geom"
  )
}

cat("\n## B) check_attack_rates: warn / impossible / infinite-rate flags\n")
{
  fx <- make_ob_shp()
  evt <- check_attack_rates(fx$ob, fx$shp, ar_warn = 100, ar_impossible = 1000)
  by_lp <- function(id) evt %>% filter(location_period_id == id)

  expect_true(
    isTRUE(by_lp("lp1")$flag_inf_rate),
    "cases > 0 with pop <= 0 -> flag_inf_rate (the zero-routing trap)"
  )
  expect_true(
    !isTRUE(by_lp("lp1")$flag_warn),
    "an infinite-rate LP is not double counted under flag_warn (rate is NA, not > threshold)"
  )
  expect_true(
    isTRUE(by_lp("lp2")$flag_warn) && isTRUE(by_lp("lp2")$flag_impossible),
    "50 cases / 0.5 pop = 1e5 per 1k -> both warn and impossible"
  )
  expect_true(
    !isTRUE(by_lp("lp4")$flag_warn),
    "20 cases / 10000 pop = 2 per 1k -> not flagged"
  )
}

cat("\n## D) check_duplicate_geometry_provenance: gates 4/5a/5b\n")
{
  cat("  D1) absent columns -> skip, not fail (pre-patch corpus)\n")
  fx <- make_ob_shp()
  expect_true(
    is.null(check_duplicate_geometry_provenance(fx$ob)),
    "pop_geom_dup_n/_class absent -> NULL (graceful skip, no error)"
  )

  cat("  D2) present columns -> mismatch + unknown-class detection\n")
  ob_dup <- tibble::tibble(
    location = rep("AFR::GIN::X", 5),
    location_period_id = c("lp1", "lp2", "lp3", "lp4", "lp5"),
    country = "GIN",
    pop_geom_dup_n = c(1, 7, 7, 2, NA),
    pop_geom_dup_class = c(
      "unique",
      "parent_inherited",
      "unique",
      "cross_unit",
      NA
    )
  )
  dup <- check_duplicate_geometry_provenance(ob_dup)
  by_lp <- function(id) dup %>% filter(location_period_id == id)

  expect_true(
    !isTRUE(by_lp("lp1")$flag_n_class_mismatch),
    "dup_n=1 & class=unique -> consistent, not flagged"
  )
  expect_true(
    !isTRUE(by_lp("lp2")$flag_n_class_mismatch),
    "dup_n=7 & class=parent_inherited -> consistent, not flagged"
  )
  expect_true(
    isTRUE(by_lp("lp3")$flag_n_class_mismatch),
    "dup_n=7 & class=unique -> inconsistent (the Guinea/Conakry signature), flagged"
  )
  expect_true(
    !isTRUE(by_lp("lp4")$flag_n_class_mismatch),
    "dup_n=2 & class=cross_unit -> consistent, not flagged"
  )
  expect_true(
    is.na(by_lp("lp5")$flag_n_class_mismatch),
    "both dup_n and class NA -> mismatch flag is NA, not TRUE/FALSE"
  )

  ob_bad_class <- tibble::tibble(
    location = "AFR::X::Y",
    location_period_id = "lp1",
    country = "X",
    pop_geom_dup_n = 3,
    pop_geom_dup_class = "not_a_real_class"
  )
  dup_bad <- check_duplicate_geometry_provenance(ob_bad_class)
  expect_true(
    isTRUE(dup_bad$flag_unknown_class),
    "pop_geom_dup_class outside the taxonomy -> flag_unknown_class"
  )
}

cat("\n## E) check_adj_factor_banding: gate 3 (refined bands)\n")
{
  ob_adj <- tibble::tibble(
    location = paste0("AFR::X::L", 1:6),
    location_period_id = paste0("lp", 1:6),
    country = "X",
    adj_factors = c(1.0, 0.68, 1.5, 0.55, 1.8, 0.3)
    # lp1: dead-center accept | lp2: just inside accept | lp3: accept boundary
    # lp4: accept_flag band [0.5,0.67) | lp5: accept_flag band (1.5,2.0]
    # lp6: out_of_band (below 0.5)
  )
  adj <- check_adj_factor_banding(ob_adj)
  by_lp <- function(id) adj %>% filter(location_period_id == id)

  expect_true(by_lp("lp1")$band == "accept", "adj_factor 1.0 -> accept")
  expect_true(
    by_lp("lp3")$band == "accept",
    "adj_factor 1.5 (upper accept boundary, inclusive) -> accept"
  )
  expect_true(
    by_lp("lp4")$band == "accept_flag",
    "adj_factor 0.55 -> accept_flag (audit's observed range never gets here)"
  )
  expect_true(
    by_lp("lp5")$band == "accept_flag",
    "adj_factor 1.8 -> accept_flag"
  )
  expect_true(
    by_lp("lp6")$band == "out_of_band" && isTRUE(by_lp("lp6")$flag_out_of_band),
    "adj_factor 0.3 (outside [0.5,2.0]) -> out_of_band, flagged"
  )

  ob_missing <- tibble::tibble(
    location = "AFR::X::L",
    location_period_id = "lp1",
    country = "X",
    adj_factors = NA_real_
  )
  expect_true(
    check_adj_factor_banding(ob_missing)$band == "missing",
    "NA adj_factor -> band 'missing', not silently dropped"
  )

  ob_no_col <- tibble::tibble(
    location = "AFR::X::L",
    location_period_id = "lp1",
    country = "X"
  )
  expect_true(
    is.null(check_adj_factor_banding(ob_no_col)),
    "adj_factors column absent -> NULL (graceful skip)"
  )
}

cat(
  "\n## F) check_crosswalk_reconciliation: P_l = sum_g H vs. independent P_direct\n"
)
{
  cat("  F1) absent columns -> skip, not fail (pre-Stage-1 corpus)\n")
  expect_true(
    is.null(check_crosswalk_reconciliation(
      tibble::tibble(lctn_pr = "lp1"),
      tibble::tibble(lctn_pr = "lp1")
    )),
    "xwalk/lp_pop missing required columns -> NULL (graceful skip)"
  )

  cat("  F2) reconciling vs. failing LPs, at two gadm_level tilings\n")
  # lp_ok:   sums to 1000 at both ADM1 and ADM2 tilings, P_direct = 1000 -> pass
  # lp_bad:  ADM1 tiling sums to 800 vs P_direct 1000 (20% off) -> fail
  # lp_noPD: no matching P_direct row at all -> flag_no_p_direct, not a reconcile failure
  xwalk <- tibble::tibble(
    lctn_pr = c("lp_ok", "lp_ok", "lp_ok", "lp_ok", "lp_bad", "lp_noPD"),
    gid = c("A.1_1", "A.2_1", "A.1.1_1", "A.1.2_1", "A.3_1", "A.4_1"),
    gadm_level = c("ADM1", "ADM1", "ADM2", "ADM2", "ADM1", "ADM1"),
    year = 2020L,
    H = c(600, 400, 600, 400, 800, 500)
  )
  lp_pop_direct <- tibble::tibble(
    lctn_pr = c("lp_ok", "lp_bad"),
    year = 2020L,
    P_direct = c(1000, 1000)
  )
  rec <- check_crosswalk_reconciliation(xwalk, lp_pop_direct, tol = 0.01)

  by_lp_lvl <- function(id, lvl) {
    rec %>% filter(lctn_pr == id, gadm_level == lvl)
  }

  expect_true(
    nrow(rec) == 4L,
    "one row per (lctn_pr, year, gadm_level) group (lp_ok x2, lp_bad, lp_noPD)"
  )
  expect_true(
    !isTRUE(by_lp_lvl("lp_ok", "ADM1")$flag_reconcile_fail),
    "lp_ok ADM1 tiling (600+400=1000 vs P_direct 1000) reconciles"
  )
  expect_true(
    !isTRUE(by_lp_lvl("lp_ok", "ADM2")$flag_reconcile_fail),
    "lp_ok ADM2 tiling (600+400=1000 vs P_direct 1000) reconciles"
  )
  expect_true(
    isTRUE(by_lp_lvl("lp_bad", "ADM1")$flag_reconcile_fail),
    "lp_bad (800 vs P_direct 1000, 20% off) fails at tol=0.01"
  )
  expect_true(
    isTRUE(rec$flag_no_p_direct[rec$lctn_pr == "lp_noPD"]),
    "lp_noPD has no matching P_direct row -> flag_no_p_direct, not a silent pass"
  )
  expect_true(
    is.na(rec$flag_reconcile_fail[rec$lctn_pr == "lp_noPD"]),
    "lp_noPD reconcile-fail flag is NA (nothing to compare), not FALSE"
  )
}

cat("\n## G) run_population_gates: orchestration returns all five tables\n")
{
  fx <- make_ob_shp()
  gates <- run_population_gates(fx$ob, fx$shp)
  expect_true(
    all(
      c("lp_pop", "evt", "geom_dup", "adj_band", "reconcile") %in% names(gates)
    ),
    "returns list(lp_pop, evt, geom_dup, adj_band, reconcile)"
  )
  expect_true(
    nrow(gates$lp_pop) == 4L,
    "lp_pop has one row per location period"
  )
  expect_true(
    is.null(gates$geom_dup),
    "geom_dup NULL when the fixture lacks pop_geom_dup_* columns"
  )
  expect_true(
    is.data.frame(gates$adj_band),
    "adj_band populated when adj_factors is present"
  )
  expect_true(
    is.null(gates$reconcile),
    "reconcile NULL when xwalk/lp_pop_direct not supplied (default args)"
  )

  cat("  G2) reconcile populated when xwalk/lp_pop_direct are supplied\n")
  xwalk2 <- tibble::tibble(
    lctn_pr = "lp1",
    gid = "A.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 1000
  )
  lp_pop_direct2 <- tibble::tibble(
    lctn_pr = "lp1",
    year = 2020L,
    P_direct = 1000
  )
  gates2 <- run_population_gates(
    fx$ob,
    fx$shp,
    xwalk = xwalk2,
    lp_pop_direct = lp_pop_direct2
  )
  expect_true(
    is.data.frame(gates2$reconcile),
    "reconcile populated when both inputs supplied"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0) {
  quit(status = 1)
}
