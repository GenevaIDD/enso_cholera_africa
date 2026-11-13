# test_simulator_enso_cb_grouping.R
#
# Unit tests for grouped ENSO cross-basis (cb_meta) routing in simulate_one_tu,
# targeting the ENSO-teleconnection grouping where each group uses its own best lag
# (enso_east=lag1, enso_horn=lag8, enso_south=lag3, enso_west=lag9 per LOO
# comparison).
#
# Fills the gap between:
#   test_simulator_grouping.R — grouped routing with cb_meta_*_by_group = NULL
#   test_simulator_dlnm.R    — DLNM cross-basis in the ungrouped path only
#
# This file tests the INTERSECTION: grouped routing WITH non-null per-group
# ENSO cb_meta (different max_lag per group). This is the configuration
# produced by the occ_grouped_enso_singlelag batch runs when ingested by the
# simulator via load_artefacts(grouped=TRUE).
#
# Tests (numbered to follow test_simulator_dlnm.R series):
#   34) .lags_required_from_draws: per-group cb_meta with distinct max_lag
#       → distinct enso_lag column requirements
#   35) simulate_one_tu grouped routing: lag-2..3 ENSO signal fires 'horn' only
#       (group 'east' sees lag0/lag1=0 → near-zero LP; group 'horn' sees
#        lag2/lag3=5 → large LP through lag-3 cb)
#   36) simulate_one_tu equivalence: when both groups share identical cb_meta
#       and draws, output matches the ungrouped single-bundle run
#
# All tests require the `dlnm` package and are skipped gracefully when absent.
#
# Run: Rscript analysis/tests/test_simulator_enso_cb_grouping.R

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/dlnm.R")   # .LATENT_CB_VARS -- required by .lags_required_from_draws()
source("analysis/utils/simulator_core.R")
source("analysis/utils/simulations.R")

# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------

.fails  <- 0L
.passes <- 0L

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

cat("# test_simulator_enso_cb_grouping.R — grouped ENSO cross-basis routing\n\n")

# ---------------------------------------------------------------------------
# dlnm availability gate
# ---------------------------------------------------------------------------

if (!requireNamespace("dlnm", quietly = TRUE)) {
  cat("[SKIP] dlnm not available — skipping all tests\n")
  cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
  if (.fails > 0L) quit(save = "no", status = 1L)
  quit(save = "no", status = 0L)
}

# ---------------------------------------------------------------------------
# Shared fixture: two ENSO groups, 2 gids each (2 adm1, 1 adm2 per adm1).
#   'east' = enso_east analogue, max_lag = 1
#   'horn' = enso_horn analogue, max_lag = 3  (real best = 8; kept small here)
# ---------------------------------------------------------------------------

fixture_frame <- function() {
  dplyr::bind_rows(
    expand.grid(group = "east", country = "EEE",
                adm1n = c(1, 2), stringsAsFactors = FALSE),
    expand.grid(group = "horn", country = "HHH",
                adm1n = c(1, 2), stringsAsFactors = FALSE)
  ) %>%
    dplyr::mutate(
      gid         = sprintf("%s_%d_1", country, adm1n),
      adm1_unit   = sprintf("%s_%d",   country, adm1n),
      AFRO_region = dplyr::if_else(group == "east", "REG_E", "REG_H"),
      tu          = "ensoTU",
      min_seas    = 1L,
      pop         = 1e5
    ) %>%
    dplyr::select(gid, group, country, AFRO_region, tu,
                  min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
}

frame_tu <- fixture_frame()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Minimal ENSO cb_meta for a crossbasis over lags 0..max_lag.
# argvar = lin (1 df, no intercept for variable space by default) ×
# arglag = lin with intercept=FALSE (1 df; arglag defaults to intercept=TRUE
# which would give 2 df and 2 cb columns — intercept=FALSE reduces to 1).
# Result: exactly 1 cb column named "enso_cb1".
make_enso_cb_meta <- function(max_lag) {
  list(enso = list(
    lag          = c(0L, max_lag),
    argvar       = list(fun = "lin"),
    arglag       = list(fun = "lin", intercept = FALSE),
    cb_col_names = "enso_cb1"
  ))
}

# Panel: 6 months with the supplied named enso_lag* values (location-invariant,
# same value across all gids for each month). Climate covariates = 0 throughout.
make_panel <- function(enso_vals) {
  # enso_vals: named numeric vector, e.g. c(enso_lag0=0, enso_lag1=5, ...)
  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  base   <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas              = lubridate::yday(t_start),
      mean_precip               = 0, dev_precip = 0,
      mean_temp_daymean         = 0, dev_temp_daymean = 0,
      ongoing_outbreak_lag1         = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1      = 0L,
      time_to_last_outbreak_cat = NA_character_
    )
  for (nm in names(enso_vals)) base[[nm]] <- enso_vals[[nm]]
  base
}

# Draws bundle for a group sub-frame with ENSO cb1 FE column.
# Uses make_draws_bundle_synthetic; all RE values = 0, so region-rooted vs.
# region-dropped key mismatch only produces a 0 contribution either way.
make_bundle <- function(frame_g, intercept = 0, fe_enso_cb = 0, phi = 5, S = 3) {
  active <- c("(Intercept)", "mean_precip", "dev_precip",
              "mean_temp_daymean", "dev_temp_daymean",
              "ongoing_outbreak_lag1", "ongoing_outbreak_country_lag1",
              "ongoing_outbreak_tu_lag1", "enso_cb1")
  make_draws_bundle_synthetic(
    frame_g,
    S            = S,
    phi          = phi,
    active_terms = active,
    fe_overrides = c(`(Intercept)` = intercept, enso_cb1 = fe_enso_cb)
  )
}

horizon <- list(start = "2024-01-01", end = "2024-06-30")

# ---------------------------------------------------------------------------
cat("## 34) .lags_required_from_draws: per-group cb_meta → distinct lag sets\n")
# ---------------------------------------------------------------------------
# 'east' cb_meta max_lag=1 → requires enso_lag0, enso_lag1 only.
# 'horn' cb_meta max_lag=3 → requires enso_lag0 .. enso_lag3.
# The two sets must differ.
{
  frame_east <- frame_tu %>% dplyr::filter(group == "east")
  frame_horn <- frame_tu %>% dplyr::filter(group == "horn")

  bun_east <- make_bundle(frame_east, fe_enso_cb = 1)
  bun_horn <- make_bundle(frame_horn, fe_enso_cb = 1)

  cb_east <- make_enso_cb_meta(max_lag = 1L)
  cb_horn <- make_enso_cb_meta(max_lag = 3L)

  req_east <- .lags_required_from_draws(bun_east, cb_east)
  req_horn <- .lags_required_from_draws(bun_horn, cb_horn)

  expect_true(all(c("enso_lag0", "enso_lag1") %in% req_east),
              "east (max_lag=1): enso_lag0 and enso_lag1 required")
  expect_true(!("enso_lag2" %in% req_east),
              "east (max_lag=1): enso_lag2 NOT required")
  expect_true(all(sprintf("enso_lag%d", 0:3) %in% req_horn),
              "horn (max_lag=3): enso_lag0..enso_lag3 required")
  expect_true(!("enso_lag4" %in% req_horn),
              "horn (max_lag=3): enso_lag4 NOT required")
  expect_true(!setequal(sort(req_east), sort(req_horn)),
              "per-group cb_meta → distinct lag column requirement sets")
}

# ---------------------------------------------------------------------------
cat("\n## 35) simulate_one_tu: grouped ENSO cb — lag-2/3 signal fires horn only\n")
# ---------------------------------------------------------------------------
# Panel: enso_lag0=0, enso_lag1=0, enso_lag2=5, enso_lag3=5.
#
# Group 'east' (max_lag=1): build_cb_mats reads enso_lag0..enso_lag1 → Q=[0,0]
#   → cb_east ≈ 0 → LP = intercept (-20) → p ≈ 0 → Y = 0.
#
# Group 'horn' (max_lag=3): build_cb_mats reads enso_lag0..enso_lag3 → Q=[0,0,5,5]
#   → cb_horn > 0 (lin×lin crossbasis of [0,0,5,5] with lag weights ∝ [0,1/3,2/3,1])
#   → LP = -5 + 5 * cb_horn >> 0 → p ≈ 1 → Y = 1.
#
# With seed=55 and 4 replicates × 2 gids × 6 months:
#   east: 48 Bernoulli trials at p ≈ 4.5e-5 → P(any event) < 0.003 → assert Y=0
#   horn: all trials at p ≈ 1 → assert Y=1
{
  frame_east <- frame_tu %>% dplyr::filter(group == "east")
  frame_horn <- frame_tu %>% dplyr::filter(group == "horn")

  # east: low intercept (dominates when cb≈0), horn: moderate intercept + high FE_cb
  bun_east_occ  <- make_bundle(frame_east, intercept = -20, fe_enso_cb = 5,  S = 4)
  bun_horn_occ  <- make_bundle(frame_horn, intercept = -5,  fe_enso_cb = 5,  S = 4)

  cb_east <- make_enso_cb_meta(max_lag = 1L)   # only sees lag 0..1 (both 0)
  cb_horn <- make_enso_cb_meta(max_lag = 3L)   # sees lag 0..3 (lag2=lag3=5)

  panel_tu <- make_panel(c(enso_lag0 = 0, enso_lag1 = 0,
                            enso_lag2 = 5, enso_lag3 = 5))

  res <- simulate_one_tu(
    tu            = "ensoTU",
    frame_tu      = frame_tu,
    panel_tu      = panel_tu,
    draws_occ     = bun_east_occ,    # positional (overridden by grouping$ below)
    replicate_ids = 1L:4L,
    seed          = 55L,
    output_part   = tempfile("enso_cb_grp35_"),
    horizon       = horizon,
    dry_run       = TRUE,
    grouping      = list(
      draws_occ_by_group    = list(east = bun_east_occ, horn = bun_horn_occ),
      cb_meta_occ_by_group  = list(east = cb_east, horn = cb_horn)
    )
  )

  east_gids <- frame_east$gid
  horn_gids <- frame_horn$gid
  ms        <- res$monthly_summary

  # horn: enso_lag2/3 non-zero with max_lag=3 → large cb → p≈1 → Y=1 everywhere
  horn_ms <- ms %>% dplyr::filter(gid %in% horn_gids)
  expect_true(all(horn_ms$Y == 1L),
              "horn group (lag-2/3 ENSO=5, max_lag=3): Y=1 in all months/reps")

  # east: enso_lag0/1 both 0 with max_lag=1 → cb=0 → LP=-20 → p≈0 → Y=0 always
  east_ms <- ms %>% dplyr::filter(gid %in% east_gids)
  expect_true(all(east_ms$Y == 0L),
              "east group (lag-0/1 ENSO=0, max_lag=1, intercept=-20): Y=0 in all months/reps")

  # event_table: all rows from horn (east is silent)
  if (nrow(res$event_table) > 0L) {
    expect_true(all(res$event_table$gid %in% horn_gids),
                "event_table: all events in horn group (east produces none)")
  } else {
    expect_true(FALSE, "event_table: expected horn events but table is empty")
  }
}

# ---------------------------------------------------------------------------
cat("\n## 36) simulate_one_tu: grouped ENSO cb equivalence (single group == ungrouped)\n")
# ---------------------------------------------------------------------------
# When both the ungrouped path (cb_meta_occ=cb) and a single-group grouped path
# (grouping$cb_meta_occ_by_group=list(all=cb)) receive identical cb_meta and
# draws, output must be bit-for-bit identical.  All RE values are 0 so the
# region-rooted vs. region-dropped key-space mismatch contributes 0 either way,
# isolating the cb routing as the only variable.
{
  frame_all <- frame_tu %>% dplyr::mutate(group = "all")

  cb_meta  <- make_enso_cb_meta(max_lag = 1L)
  panel_tu <- make_panel(c(enso_lag0 = 1, enso_lag1 = 2))

  bun_occ  <- make_bundle(frame_all, intercept = -2, fe_enso_cb = 0.5, S = 4)

  common <- list(
    tu            = "ensoTU",
    frame_tu      = frame_all,
    panel_tu      = panel_tu,
    draws_occ     = bun_occ,
    replicate_ids = 1L:4L,
    seed          = 99L,
    horizon       = horizon,
    dry_run       = TRUE
  )

  # Ungrouped path: cb_meta_occ supplied directly
  res_ug <- do.call(simulate_one_tu, c(common, list(
    output_part  = tempfile("enso_equiv36_ug_"),
    cb_meta_occ  = cb_meta
  )))

  # Grouped path: single group "all" wrapping same cb_meta and draws
  res_g <- do.call(simulate_one_tu, c(common, list(
    output_part = tempfile("enso_equiv36_g_"),
    grouping    = list(
      draws_occ_by_group    = list(all = bun_occ),
      cb_meta_occ_by_group  = list(all = cb_meta)
    )
  )))

  expect_equal(res_ug$monthly_summary$Y, res_g$monthly_summary$Y,
               "grouped ENSO cb equiv: Y identical — ungrouped vs single-group grouped")
  expect_equal(nrow(res_ug$event_table), nrow(res_g$event_table),
               "grouped ENSO cb equiv: event_table row count identical")
}

# ---------------------------------------------------------------------------
cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) quit(save = "no", status = 1L)
