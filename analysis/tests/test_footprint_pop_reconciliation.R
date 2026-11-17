# test_footprint_pop_reconciliation.R — HARD: the (gid, season) population
#   denominator, reconcile_footprint_pop().
#
# Usage:
#   Rscript analysis/tests/test_footprint_pop_reconciliation.R
#
# Details:
#   `pop = sum(pop, na.rm = TRUE)` in analysis/03_prepare_siz.R aggregated a
#   STOCK with a flow rule: two location periods both reporting the whole gid
#   ADDED their populations (the denominator doubled), and sum() of an EMPTY
#   group returned 0, which travelled downstream as a population.
#   reconcile_footprint_pop() replaces it with pop = min( sum_l H[l,g] , Q_g ):
#   disjoint pieces still sum, overlapping reports SATURATE at the canonical
#   gid population, zero is never a population. Pure synthetic tables: no sf,
#   raster, SQLite, network or data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/population.R") # reconcile_footprint_pop() -- used directly
source("analysis/utils/composite_lp.R") # transitive dependency of population.R (utils/README.md)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-9) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}

# ---------------------------------------------------------------------------
# Fixture — the canonical gid table Q_g. Deliberately NOT ADM2-shaped ids for
# the missing-Q cases, so resolve_gid_population()'s T3 parent_share tier cannot
# quietly rescue them and mask the branch under test.
# ---------------------------------------------------------------------------
gid_pop <- tibble(
  gid = c("AAA.1_1", "AAA.2_1", "AAA.3_1", "AAA.4_1"),
  year = 2020L,
  Q = c(100000, 50000, 20000, 8000),
  iso3 = "AAA",
  gadm_level = "ADM1"
)

cat("# footprint population reconciliation\n")

# --- 1. disjoint pieces sum -------------------------------------------------
cat("\n[1] disjoint location periods tiling a gid still SUM\n")
r <- reconcile_footprint_pop(
  gid = "AAA.1_1",
  year = 2020L,
  H_sum = 30000 + 45000,
  gid_pop = gid_pop
)
expect(
  near(r$pop, 75000) && r$pop_tier == "footprint",
  "two ADM3 children of 30000 + 45000 give 75000, below Q = 100000"
)
expect(
  nrow(r) == 1L && identical(names(r), c("gid", "year", "pop", "pop_tier")),
  "returns one row per bucket with the documented columns"
)

# --- 2. overlapping reports saturate ---------------------------------------
cat("\n[2] overlapping location periods SATURATE, they do not double\n")
r2 <- reconcile_footprint_pop("AAA.2_1", 2020L, H_sum = 2 * 50000, gid_pop)
expect(
  near(r2$pop, 50000) && r2$pop_tier == "gid_capped",
  "two whole-gid reports give Q = 50000, not 100000 — the NGA_19_1 defect"
)
expect(
  sum(c(50000, 50000)) == 100000,
  "(for contrast: the sum() rule this replaces would have returned 100000)"
)
r2b <- reconcile_footprint_pop("AAA.2_1", 2020L, H_sum = 6 * 50000, gid_pop)
expect(
  near(r2b$pop, 50000),
  "the cap is absolute — six overlapping reports still give Q, not 6Q"
)

# --- 3. exact equality is not a cap ----------------------------------------
cat("\n[3] H_sum == Q is the normal case, not an anomaly\n")
r3 <- reconcile_footprint_pop("AAA.3_1", 2020L, H_sum = 20000, gid_pop)
expect(
  near(r3$pop, 20000) && r3$pop_tier == "footprint",
  "a single LP that IS the gid reports 'footprint', not 'gid_capped'"
)
r3b <- reconcile_footprint_pop(
  "AAA.3_1",
  2020L,
  H_sum = 20000 * (1 + 1e-12),
  gid_pop
)
expect(
  r3b$pop_tier == "footprint",
  "floating-point dust above Q does not trip a spurious cap"
)
r3c <- reconcile_footprint_pop("AAA.3_1", 2020L, H_sum = 20000 * 1.001, gid_pop)
expect(
  r3c$pop_tier == "gid_capped",
  "a real 0.1% excess DOES trip the cap — the slack is dust-sized, not a tolerance"
)

# --- 4. the COD_4_6_1 zero --------------------------------------------------
cat("\n[4] zero is never a population\n")
r4 <- reconcile_footprint_pop("AAA.4_1", 2020L, H_sum = 0, gid_pop)
expect(
  near(r4$pop, 8000) && r4$pop_tier == "gid_fallback",
  "H_sum = 0 (sum() of an empty group) falls back to Q — the COD_4_6_1 fix"
)
expect(
  near(reconcile_footprint_pop("AAA.4_1", 2020L, NA_real_, gid_pop)$pop, 8000),
  "H_sum = NA falls back the same way"
)
expect(
  near(reconcile_footprint_pop("AAA.4_1", 2020L, -5, gid_pop)$pop, 8000),
  "a negative H_sum is treated as absent, never propagated"
)

# --- 5. no canonical Q ------------------------------------------------------
cat("\n[5] a gid absent from the canonical table\n")
r5 <- reconcile_footprint_pop("QQQ.9_1", 2020L, H_sum = 1234, gid_pop)
expect(
  near(r5$pop, 1234) && r5$pop_tier == "lp_mass_only",
  "with no Q_g the LP mass is used, and TAGGED so the fallback is countable"
)
expect(
  inherits(
    try(reconcile_footprint_pop("QQQ.9_1", 2020L, 0, gid_pop), silent = TRUE),
    "try-error"
  ),
  "no Q and no mass is a hard stop — it refuses to invent a denominator"
)
r5b <- reconcile_footprint_pop("QQQ.9_1", 2020L, 0, gid_pop, on_missing = "na")
expect(
  is.na(r5b$pop) && r5b$pop_tier == "missing",
  "on_missing = 'na' yields NA + tier 'missing', for QC reporting only"
)

# --- 6. vectorised, order-preserving ---------------------------------------
cat("\n[6] vectorised over buckets, in input order\n")
g <- c("AAA.4_1", "AAA.1_1", "QQQ.9_1", "AAA.2_1")
rv <- reconcile_footprint_pop(g, 2020L, c(0, 75000, 1234, 999999), gid_pop)
expect(
  identical(rv$gid, g),
  "output rows are in input order — safe to bind back onto a grouped frame"
)
expect(
  near(rv$pop, c(8000, 75000, 1234, 50000)),
  "each bucket takes its own branch independently"
)
expect(
  identical(
    rv$pop_tier,
    c("gid_fallback", "footprint", "lp_mass_only", "gid_capped")
  ),
  "all four tiers appear in one vectorised call"
)
expect(
  !any(is.na(rv$pop)) && all(rv$pop > 0),
  "no NA and no zero survives the default on_missing = 'error' path"
)

# --- 7. the year is honoured ------------------------------------------------
cat("\n[7] Q_g is resolved per year, not per gid\n")
gp2 <- bind_rows(
  gid_pop,
  gid_pop %>% filter(gid == "AAA.2_1") %>% mutate(year = 2025L, Q = 60000)
)
r7 <- reconcile_footprint_pop(
  c("AAA.2_1", "AAA.2_1"),
  c(2020L, 2025L),
  c(999999, 999999),
  gp2
)
expect(
  near(r7$pop, c(50000, 60000)),
  "the same gid caps at 50000 in 2020 and 60000 in 2025"
)

wpp <- tibble(iso3 = "AAA", year = c(2020L, 2040L), pop_total = c(1e6, 2e6))
r7b <- reconcile_footprint_pop("AAA.2_1", 2040L, 999999, gid_pop, wpp = wpp)
expect(
  near(r7b$pop, 100000),
  "an out-of-range year rescales Q by the WPP national ratio (2x) before capping"
)
expect(
  inherits(
    try(
      reconcile_footprint_pop("AAA.2_1", 2040L, 0, gid_pop, wpp = NULL),
      silent = TRUE
    ),
    "try-error"
  ),
  "without a WPP table that same year has no Q and stops, rather than guessing"
)

# A MALFORMED wpp table is not an ABSENT one. NULL means "T2 unavailable" and is
# a legitimate caller choice; a table whose population column is misnamed is a
# broken cache. Before hardening, the latter surfaced as `argument is of length
# zero` from inside resolve_gid_population()'s T2 loop -- four frames from the
# cause. This test exists because it is the failure this fixture originally hit.
bad_wpp <- tibble(iso3 = "AAA", year = c(2020L, 2040L), pop = c(1e6, 2e6))
err <- try(
  reconcile_footprint_pop("AAA.2_1", 2040L, 999999, gid_pop, wpp = bad_wpp),
  silent = TRUE
)
expect(
  inherits(err, "try-error") &&
    grepl("pop_total", conditionMessage(attr(err, "condition")), fixed = TRUE),
  "a WPP table missing 'pop_total' errors BY NAME, not with a length-zero crash"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
