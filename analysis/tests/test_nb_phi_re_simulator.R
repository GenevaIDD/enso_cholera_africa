# test_nb_phi_re_simulator.R — unit tests for compute_phi_vec().
#
# compute_phi_vec() is track-agnostic and serves the monthly (mon) NegBin
# dispersion. Fixtures carry no mean_rate column.
#
# Hard-fail gate. Tests the location-specific phi broadcast logic without
# requiring a real brms fit object.  .extract_phi_re_draws() is NOT tested
# here (it calls brms::ranef() and needs a real fit); its correctness is
# validated in the end-to-end smoke-fit step.
#
# Run: Rscript analysis/tests/test_nb_phi_re_simulator.R

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R")   # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

# ---------------------------------------------------------------------------
# Test harness (mirrors test_simulator_kernel.R)
# ---------------------------------------------------------------------------

.fails <- 0L; .passes <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
expect_equal <- function(a, b, label, tol = 1e-9) {
  ok <- isTRUE(all.equal(a, b, tolerance = tol))
  expect_true(ok, paste0(label,
                          if (!ok) sprintf("  [got=%s want=%s]",
                                           paste(a, collapse=","),
                                           paste(b, collapse=",")) else ""))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Minimal frame_tu: 4 locations across 2 regions and 2 countries.
frame_tu <- tibble::tibble(
  gid         = c("ETH_1_1", "ETH_2_1", "MOZ_1_1", "ZAF_1_1"),
  AFRO_region = c("Eastern.Africa", "Eastern.Africa", "Southern.Africa", "Southern.Africa"),
  country     = c("ETH",            "ETH",            "MOZ",             "ZAF"),
  adm1_unit   = c("ETH_1",          "ETH_2",          "MOZ_1",           "ZAF_1"),
  pop         = c(1e6, 2e6, 1.5e6, 3e6)
)
N <- nrow(frame_tu)
S <- 3L   # small posterior for tests

# Helper: build a NB draws bundle with phi constant
mk_draws_constant <- function(phi_vals = c(2.0, 3.5, 1.0)) {
  list(phi = phi_vals, phi_re = "constant")
}

# Helper: build an NB draws bundle with phi as a list (non-constant)
mk_draws_phi_list <- function(phi_re, fe, region = NULL, country = NULL, adm1 = NULL) {
  list(
    phi    = list(fe = fe, region = region, country = country, adm1 = adm1),
    phi_re = phi_re
  )
}

# ---------------------------------------------------------------------------
# 1. constant phi_re: returns rep(scalar, N)
# ---------------------------------------------------------------------------
cat("\n[1] phi_re = constant: scalar broadcast\n")
ds <- mk_draws_constant(c(2.0, 3.5, 1.0))
v <- compute_phi_vec(ds, s = 1L, frame_tu)
expect_true(length(v) == N,            "constant: length == N")
expect_true(all(v == 2.0),             "constant: all values = phi[1] = 2.0")
v2 <- compute_phi_vec(ds, s = 2L, frame_tu)
expect_true(all(v2 == 3.5),            "constant: draw 2 = 3.5")

# ---------------------------------------------------------------------------
# 2. constant phi_re with list-format phi (backward compat branch)
# ---------------------------------------------------------------------------
cat("\n[2] phi_re = constant, phi stored as list (backward compat)\n")
ds_list <- mk_draws_phi_list("constant", fe = c(2.0, 3.5, 1.0))
v <- compute_phi_vec(ds_list, s = 1L, frame_tu)
expect_true(all(v == 2.0),             "constant list: phi$fe used as scalar")

# ---------------------------------------------------------------------------
# 3. region-level phi RE: FE + region RE, exp()
# ---------------------------------------------------------------------------
cat("\n[3] phi_re = region\n")
# Region RE matrix: S × 2 (Eastern.Africa, Southern.Africa)
re_region <- matrix(
  c( 0.5, -0.5,   # draw 1
     0.2,  0.3,   # draw 2
    -0.1,  0.1),  # draw 3
  nrow = S, byrow = TRUE,
  dimnames = list(NULL, c("Eastern.Africa", "Southern.Africa"))
)
fe <- c(1.0, 1.0, 1.0)  # log scale FE intercept (same for all draws)
ds <- mk_draws_phi_list("region", fe = fe, region = re_region)

v1 <- compute_phi_vec(ds, s = 1L, frame_tu)
# Draw 1: FE=1.0
#   Eastern.Africa locations (1,2): log_phi = 1.0 + 0.5 = 1.5  → exp(1.5)
#   Southern.Africa locations (3,4): log_phi = 1.0 + (-0.5) = 0.5 → exp(0.5)
expect_equal(v1[1], exp(1.5), "region draw1: Eastern.Africa loc1")
expect_equal(v1[2], exp(1.5), "region draw1: Eastern.Africa loc2")
expect_equal(v1[3], exp(0.5), "region draw1: Southern.Africa loc3")
expect_equal(v1[4], exp(0.5), "region draw1: Southern.Africa loc4")
expect_true(all(v1 > 0),               "region: all phi > 0")

# ---------------------------------------------------------------------------
# 4. region_country: accumulates region + country contributions
# ---------------------------------------------------------------------------
cat("\n[4] phi_re = region_country\n")
re_country <- matrix(
  c( 0.1, -0.2,  0.0,   # draw 1: ETH, MOZ, ZAF (region-prefixed keys)
     0.0,  0.1, -0.1,
     0.2,  0.0,  0.0),
  nrow = S, byrow = TRUE,
  dimnames = list(NULL, c("Eastern.Africa_ETH", "Southern.Africa_MOZ", "Southern.Africa_ZAF"))
)
ds <- mk_draws_phi_list("region_country", fe = fe,
                         region = re_region, country = re_country)

v1 <- compute_phi_vec(ds, s = 1L, frame_tu)
# Draw 1 / loc 1 (Eastern.Africa, ETH, ETH_1):
#   log_phi = 1.0 + 0.5 (region) + 0.1 (country ETH) = 1.6 → exp(1.6)
expect_equal(v1[1], exp(1.6), "region_country draw1: loc1 ETH")
# Draw 1 / loc 3 (Southern.Africa, MOZ, MOZ_1):
#   log_phi = 1.0 + (-0.5) + (-0.2) = 0.3 → exp(0.3)
expect_equal(v1[3], exp(0.3), "region_country draw1: loc3 MOZ")
# Draw 1 / loc 4 (Southern.Africa, ZAF, ZAF_1):
#   log_phi = 1.0 + (-0.5) + 0.0 = 0.5 → exp(0.5)
expect_equal(v1[4], exp(0.5), "region_country draw1: loc4 ZAF")

# ---------------------------------------------------------------------------
# 5. region_country_adm1: accumulates all three levels
# ---------------------------------------------------------------------------
cat("\n[5] phi_re = region_country_adm1\n")
re_adm1 <- matrix(
  c(0.05, -0.05, 0.10, -0.10,  # draw 1
    0.00,  0.00, 0.00,  0.00,  # draw 2
    0.01, -0.01, 0.02, -0.02), # draw 3
  nrow = S, byrow = TRUE,
  dimnames = list(NULL, c("Eastern.Africa_ETH_ETH_1",
                           "Eastern.Africa_ETH_ETH_2",
                           "Southern.Africa_MOZ_MOZ_1",
                           "Southern.Africa_ZAF_ZAF_1"))
)
ds <- mk_draws_phi_list("region_country_adm1", fe = fe,
                         region = re_region, country = re_country, adm1 = re_adm1)

v1 <- compute_phi_vec(ds, s = 1L, frame_tu)
# Draw 1 / loc 1 (Eastern.Africa, ETH, ETH_1):
#   log_phi = 1.0 + 0.5 + 0.1 + 0.05 = 1.65 → exp(1.65)
expect_equal(v1[1], exp(1.65), "adm1 draw1: loc1 ETH_1")
# Draw 1 / loc 2 (Eastern.Africa, ETH, ETH_2):
#   log_phi = 1.0 + 0.5 + 0.1 + (-0.05) = 1.55 → exp(1.55)
expect_equal(v1[2], exp(1.55), "adm1 draw1: loc2 ETH_2")
expect_true(all(v1 > 0),               "adm1: all phi > 0")

# ---------------------------------------------------------------------------
# 6. Unmatched group key: location not in RE matrix → FE + higher-level RE only
# ---------------------------------------------------------------------------
cat("\n[6] Unmatched ADM1 key: falls back to FE + region + country only\n")
# Use a frame with a location whose adm1_unit is not in the re_adm1 matrix
frame_unseen <- frame_tu[1, ]  # ETH_1 — IS in re_adm1 matrix; swap out
frame_unseen$adm1_unit <- "ETH_99"  # unknown ADM1
ds_partial <- mk_draws_phi_list("region_country_adm1", fe = fe,
                                  region = re_region, country = re_country, adm1 = re_adm1)
v_unseen <- compute_phi_vec(ds_partial, s = 1L, frame_unseen)
# Falls back to FE + region + country (no adm1 contribution):
# Draw 1: 1.0 + 0.5 + 0.1 = 1.6 → exp(1.6)
expect_equal(v_unseen[1], exp(1.6),    "unseen adm1: falls back to region+country")
expect_true(!is.na(v_unseen[1]),       "unseen adm1: no NA produced")

# ---------------------------------------------------------------------------
# 7. Scale guarantee: exp always produces > 0 (even with negative RE)
# ---------------------------------------------------------------------------
cat("\n[7] Scale: result always > 0\n")
re_extreme <- matrix(c(-10, -10), nrow = 1,
                     dimnames = list(NULL, c("Eastern.Africa", "Southern.Africa")))
ds_extreme <- mk_draws_phi_list("region", fe = c(-5), region = re_extreme)
v_extreme <- compute_phi_vec(ds_extreme, s = 1L, frame_tu)
expect_true(all(v_extreme > 0),        "extreme negative RE: phi still > 0")
expect_true(all(is.finite(v_extreme)), "extreme negative RE: phi finite")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n=== %d passed, %d failed ===\n", .passes, .fails))
if (.fails > 0) quit(save = "no", status = 1)
message("PASS  test_nb_phi_re_simulator")
