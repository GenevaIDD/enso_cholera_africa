# test_lp_observation_operator.R — HARD: the LP observation operator
#   c[g<-l] = H[l,g] / Q_g, and its application to simulated gid cases.
#
# Usage:
#   Rscript analysis/tests/test_lp_observation_operator.R
#
# Details:
#   The simulator emits absolute cases per gid (Q_g * r_g), so mapping into
#   LP space is a COVERAGE-weighted sum, not a share-weighted one:
#     E[y_l] = sum_g H[l,g] * r_g = sum_g c[g<-l] * sim_cases_g
#   which stays correct when an LP covers only PART of a gid (c < 1).
#   Exercises build_lp_observation_operator() and
#   expected_lp_from_gid_rates() on small synthetic crosswalk tables -- no
#   sf, no raster, no network, no data/ access.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/population.R") # build_lp_observation_operator()/expected_lp_from_gid_rates() -- used directly
source("analysis/utils/composite_lp.R") # transitive dependency of population.R (utils/README.md)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-8) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}

# --- A. Exact recovery: LP tiles a whole gid (c == 1) -----------------------
cat("## A) LP tiles a whole gid exactly (c == 1) -- exact recovery\n")
{
  xwalk <- tibble(
    lctn_pr = "lp_whole",
    gid = "AAA.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 1000,
    Q = 1000,
    c = 1.0
  )
  cases_by_gid <- tibble(gid = "AAA.1_1", cases = 250)
  out <- expected_lp_from_gid_rates(xwalk, cases_by_gid, year = 2020L)

  expect(nrow(out) == 1L, "one LP in, one LP out")
  expect(
    near(out$expected_cases[out$lctn_pr == "lp_whole"], 250),
    "c == 1 -> expected_cases == the gid's simulated cases exactly"
  )
}

# --- B. Half-coverage halves ------------------------------------------------
cat("\n## B) half-coverage halves the attributed cases\n")
{
  xwalk <- tibble(
    lctn_pr = "lp_half",
    gid = "AAA.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 500,
    Q = 1000,
    c = 0.5
  )
  cases_by_gid <- tibble(gid = "AAA.1_1", cases = 250)
  out <- expected_lp_from_gid_rates(xwalk, cases_by_gid, year = 2020L)

  expect(
    near(out$expected_cases[out$lctn_pr == "lp_half"], 125),
    "c == 0.5 -> expected_cases == half the gid's simulated cases"
  )
}

# --- C. A gid outside the LP's crosswalk rows contributes 0 -----------------
cat("\n## C) gid outside the LP contributes 0, does not leak in\n")
{
  xwalk <- tibble(
    lctn_pr = "lp_a",
    gid = "AAA.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 1000,
    Q = 1000,
    c = 1.0
  )
  # cases_by_gid carries a SECOND gid the crosswalk never maps lp_a to.
  cases_by_gid <- tibble(gid = c("AAA.1_1", "AAA.9_1"), cases = c(250, 99999))
  out <- expected_lp_from_gid_rates(xwalk, cases_by_gid, year = 2020L)

  expect(
    nrow(out) == 1L,
    "only lp_a appears -- the unmapped gid's cases produce no LP row"
  )
  expect(
    near(out$expected_cases[out$lctn_pr == "lp_a"], 250),
    "the huge unmapped-gid case count (99999) does not leak into lp_a's total"
  )
}

# --- D. Composite LP: coverage-weighted sum over multiple gids -------------
cat("\n## D) composite LP: sum_g c[g<-l] * cases_g\n")
{
  xwalk <- tibble(
    lctn_pr = "lp_comp",
    gid = c("AAA.2_1", "AAA.3_1", "AAA.4_1"),
    gadm_level = "ADM1",
    year = 2020L,
    H = c(400, 350, 250),
    Q = c(500, 400, 300)
  ) %>%
    mutate(c = H / Q)
  cases_by_gid <- tibble(
    gid = c("AAA.2_1", "AAA.3_1", "AAA.4_1"),
    cases = c(40, 35, 25)
  )
  out <- expected_lp_from_gid_rates(xwalk, cases_by_gid, year = 2020L)

  expected <- (400 / 500) * 40 + (350 / 400) * 35 + (250 / 300) * 25
  expect(
    near(out$expected_cases[out$lctn_pr == "lp_comp"], expected),
    "composite LP recovers the exact coverage-weighted sum across its 3 gids"
  )
}

# --- E. `by` grouping columns carried through ------------------------------
cat("\n## E) extra `by` grouping columns (e.g. rep_id) are preserved\n")
{
  xwalk <- tibble(
    lctn_pr = "lp_whole",
    gid = "AAA.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 1000,
    Q = 1000,
    c = 1.0
  )
  cases_by_gid <- tibble(
    gid = rep("AAA.1_1", 4),
    rep_id = c(1, 1, 2, 2),
    month = c(1, 2, 1, 2),
    cases = c(10, 20, 30, 40)
  )
  out <- expected_lp_from_gid_rates(
    xwalk,
    cases_by_gid,
    year = 2020L,
    by = c("rep_id", "month")
  )

  expect(
    nrow(out) == 4L,
    "one row per (lctn_pr, rep_id, month) -- by-columns not collapsed"
  )
  expect(
    near(sort(out$expected_cases), sort(cases_by_gid$cases)),
    "with c == 1 and no aggregation collision, each (rep_id, month) recovers its own case count exactly"
  )
}

# --- F. Zero / NA coverage rows are excluded from the operator -------------
cat("\n## F) c == 0 or c == NA rows are excluded from the operator\n")
{
  xwalk <- tibble(
    lctn_pr = c("lp_ok", "lp_zero_c", "lp_na_c"),
    gid = c("AAA.1_1", "AAA.2_1", "AAA.3_1"),
    gadm_level = "ADM1",
    year = 2020L,
    H = c(1000, 0, NA_real_),
    Q = c(1000, 500, 300),
    c = c(1.0, 0, NA_real_)
  )
  op <- build_lp_observation_operator(xwalk, year = 2020L, gadm_level = "ADM1")

  expect(
    nrow(op) == 1L,
    "only the c > 0, non-NA row survives into the operator"
  )
  expect(
    op$lctn_pr == "lp_ok",
    "the surviving row is the clean LP, not the zero/NA-c ones"
  )
}

# --- G. gadm_level filtering -------------------------------------------------
cat("\n## G) gadm_level filtering picks the right tiling\n")
{
  xwalk <- tibble(
    lctn_pr = c("lp_x", "lp_x"),
    gid = c("AAA.1_1", "AAA.1.1_1"),
    gadm_level = c("ADM1", "ADM2"),
    year = 2020L,
    H = c(1000, 1000),
    Q = c(1000, 800),
    c = H / Q
  )
  op1 <- build_lp_observation_operator(xwalk, year = 2020L, gadm_level = "ADM1")
  op2 <- build_lp_observation_operator(xwalk, year = 2020L, gadm_level = "ADM2")

  expect(
    nrow(op1) == 1L && op1$gid == "AAA.1_1",
    "ADM1 tiling returns the ADM1 gid only"
  )
  expect(
    nrow(op2) == 1L && op2$gid == "AAA.1.1_1",
    "ADM2 tiling returns the ADM2 gid only"
  )
  expect(
    near(op2$c, 1.25),
    "ADM2 tiling's own c is used, not mixed with ADM1's"
  )
}

# --- H. no coverage rows for a year -> hard stop() --------------------------
cat(
  "\n## H) no coverage rows for the requested year -> stop(), not an empty result\n"
)
{
  xwalk <- tibble(
    lctn_pr = "lp_x",
    gid = "AAA.1_1",
    gadm_level = "ADM1",
    year = 2020L,
    H = 1000,
    Q = 1000,
    c = 1.0
  )
  ok <- tryCatch(
    {
      build_lp_observation_operator(xwalk, year = 2099L, gadm_level = "ADM1")
      FALSE
    },
    error = function(e) TRUE
  )
  expect(
    ok,
    "requesting a year absent from the crosswalk stops rather than returning zero rows silently"
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
