# test_case_apportionment.R — HARD: composite LP -> member gid case
#   apportionment.
#
# Usage:
#   Rscript analysis/tests/test_case_apportionment.R
#
# Details:
#   The centroid collapse (analysis/02_preprocess_climate.R) reduces every
#   location period to ONE st_centroid(), so a composite LP (a union of K admin
#   units) dumps all K units' cases AND population onto one arbitrary gid
#   while the other K-1 are recorded as never having had an outbreak.
#   apportion_composite_cases() replaces that with a population-share split
#   w[l->g] = H[l,g] / P_l; sections [2]-[4] pin the three properties that
#   make it defensible. Pure synthetic tables: no sf, raster, SQLite, network
#   or data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/population.R") # .nearest_crosswalk_year()/apportion_composite_cases()/reconcile_footprint_pop() -- used directly
source("analysis/utils/composite_lp.R") # transitive dependency of population.R (utils/README.md)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-9) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}

# ---------------------------------------------------------------------------
# Fixture — four LPs, deliberately covering every branch
#
#   lp_single  1 gid   H = 10000                     -> must be a NO-OP
#   lp_comp7   7 gids  H = 10000 each  (P_l = 70000) -> the Sudan signature
#   lp_uneven  3 gids  H = 1000/2000/7000 (P_l=1e4)  -> unequal shares
#   lp_absent  not in the crosswalk at all           -> passthrough
#
# Crosswalk years are 2015 and 2020 only (WorldPop constrained coverage), while
# the epi corpus runs 2010-2023 — so year snapping is exercised, not assumed.
# ---------------------------------------------------------------------------
mk_xw <- function(lp, gids, H, year, level = "ADM1") {
  tibble(
    lctn_pr = lp,
    gid = gids,
    gadm_level = level,
    year = as.integer(year),
    H = as.numeric(H)
  ) %>%
    group_by(lctn_pr, gadm_level, year) %>%
    mutate(P_l = sum(H), w = H / P_l, Q = H, c = 1) %>%
    ungroup()
}

xwalk <- bind_rows(
  mk_xw("lp_single", "AAA.1_1", 10000, 2015),
  mk_xw("lp_single", "AAA.1_1", 11000, 2020),
  mk_xw("lp_comp7", sprintf("AAA.%d_1", 2:8), rep(10000, 7), 2015),
  mk_xw("lp_comp7", sprintf("AAA.%d_1", 2:8), rep(11000, 7), 2020),
  mk_xw("lp_uneven", sprintf("BBB.%d_1", 1:3), c(1000, 2000, 7000), 2020),
  # the SAME LPs at ADM2 — a finer tiling of lp_single into two halves
  mk_xw("lp_single", c("AAA.1.1_1", "AAA.1.2_1"), c(6000, 4000), 2015, "ADM2"),
  mk_xw("lp_single", c("AAA.1.1_1", "AAA.1.2_1"), c(6600, 4400), 2020, "ADM2")
)

dat <- tibble(
  location_period_id = c("lp_single", "lp_comp7", "lp_uneven", "lp_absent"),
  gid = c("AAA.1_1", "AAA.5_1", "BBB.3_1", "ZZZ.9_1"),
  admin_level = c("ADM1", "ADM1", "ADM1", "ADM1"),
  year_seas = c(2015L, 2015L, 2020L, 2020L),
  suspected_cases = c(500, 7000, 1000, 42),
  pop = c(10000, 70000, 10000, 999)
)

cat("# composite case apportionment\n")

out <- apportion_composite_cases(dat, xwalk)
rep <- attr(out, "apportion_report")

# --- 1. shape ---------------------------------------------------------------
cat("\n[1] one row per (location period x member gid)\n")
expect(
  nrow(out) == 1L + 7L + 3L + 1L,
  "4 input rows -> 12 output rows (1 + 7 + 3 members, + 1 unmatched)"
)
expect(
  setequal(
    out$gid[out$location_period_id == "lp_comp7"],
    sprintf("AAA.%d_1", 2:8)
  ),
  "all seven member gids of the composite are present"
)
expect(
  !("AAA.5_1" %in% out$gid[out$location_period_id == "lp_comp7"]) == FALSE,
  "the original centroid gid is among them, not replaced by a synthetic id"
)
expect(
  all(out$gid_centroid[out$location_period_id == "lp_comp7"] == "AAA.5_1"),
  "the pre-apportionment centroid match is preserved as gid_centroid"
)

# --- 2. conservation of cases ----------------------------------------------
# CASES ARE CONSERVED: sum_g cases * w[l->g] == cases, exactly.
cat("\n[2] cases are conserved exactly\n")
csum <- out %>%
  group_by(location_period_id) %>%
  summarise(s = sum(suspected_cases), .groups = "drop")
ref <- dat %>% select(location_period_id, ref = suspected_cases)
chk <- left_join(csum, ref, by = "location_period_id")
expect(
  near(chk$s, chk$ref),
  sprintf("sum of apportioned cases == original, for all %d LPs", nrow(chk))
)
expect(
  near(sum(out$suspected_cases), sum(dat$suspected_cases)),
  "corpus total is unchanged (8542 cases in, 8542 out)"
)

# --- 3. conservation of population -----------------------------------------
# POPULATION IS CONSERVED: sum_g H[l,g] == P_l, exactly.
cat("\n[3] population is conserved exactly\n")
psum <- out %>%
  filter(pop_tier == "crosswalk_H") %>%
  group_by(location_period_id) %>%
  summarise(s = sum(pop), .groups = "drop")
expect(
  near(psum$s[psum$location_period_id == "lp_comp7"], 70000),
  "sum_g H == P_l for the 7-way composite (70000)"
)
expect(
  near(psum$s[psum$location_period_id == "lp_uneven"], 10000),
  "sum_g H == P_l for the uneven 3-way split (10000)"
)
expect(
  all(
    near(sum(out$apportion_w[out$location_period_id == "lp_comp7"]), 1),
    near(sum(out$apportion_w[out$location_period_id == "lp_uneven"]), 1)
  ),
  "w sums to 1 within every apportioned LP"
)

# --- 4. the attack rate is invariant ---------------------------------------
# THE ATTACK RATE IS INVARIANT: cases and pop are divided by the SAME w, so
# cases/pop per member equals the LP's. This is what makes the split invisible
# to read_size_data()'s max_attack_rate guard; an apportionment that broke
# this would silently manufacture implausible incidence.
cat("\n[4] attack rate is invariant under the split\n")
ar_in <- dat$suspected_cases / dat$pop
names(ar_in) <- dat$location_period_id
ar_out <- out %>%
  filter(pop_tier == "crosswalk_H") %>%
  mutate(ar = suspected_cases / pop)
expect(
  near(
    ar_out$ar[ar_out$location_period_id == "lp_comp7"],
    rep(ar_in[["lp_comp7"]], 7)
  ),
  "every member of the composite carries the LP's attack rate (0.1)"
)
expect(
  near(
    ar_out$ar[ar_out$location_period_id == "lp_uneven"],
    rep(ar_in[["lp_uneven"]], 3)
  ),
  "unequal shares still give a constant attack rate — cases and pop scale together"
)

# --- 5. uneven shares are proportional to H, not equal ----------------------
cat("\n[5] the split follows population mass, not member count\n")
u <- out %>% filter(location_period_id == "lp_uneven") %>% arrange(gid)
expect(
  near(u$apportion_w, c(0.1, 0.2, 0.7)),
  "w = H / P_l = 1000/2000/7000 over 10000"
)
expect(
  near(u$suspected_cases, c(100, 200, 700)),
  "1000 cases split 100/200/700 — NOT 333.3 each"
)

# --- 6. a single-gid LP is a no-op -----------------------------------------
cat("\n[6] non-composite LPs pass through untouched\n")
s <- out %>% filter(location_period_id == "lp_single")
expect(
  nrow(s) == 1L && near(s$suspected_cases, 500) && near(s$apportion_w, 1),
  "single-member LP keeps all its cases and w == 1"
)
expect(
  isFALSE(s$apportioned) && s$n_apportion_members == 1L,
  "single-member LP is not flagged apportioned"
)
expect(
  all(out$apportioned[out$location_period_id == "lp_comp7"]) &&
    all(out$n_apportion_members[out$location_period_id == "lp_comp7"] == 7L),
  "the composite IS flagged, with n_apportion_members == 7"
)

# --- 7. unmatched LPs are loud, never silently dropped ----------------------
cat("\n[7] an LP absent from the crosswalk is tagged, not vanished\n")
a <- out %>% filter(location_period_id == "lp_absent")
expect(
  nrow(a) == 1L && a$pop_tier == "unmatched_lp",
  "unmatched LP survives with pop_tier == 'unmatched_lp'"
)
expect(
  near(a$suspected_cases, 42) && near(a$pop, 999),
  "its cases and upstream pop are untouched"
)
expect(
  rep$n_lp_unmatched == 1L && rep$n_rows_composite == 10L,
  "the attached report counts 1 unmatched LP and 10 composite rows"
)
expect(
  nrow(apportion_composite_cases(dat, xwalk, keep_unmatched = FALSE)) == 11L,
  "keep_unmatched = FALSE drops it explicitly (11 rows, not 12)"
)

# --- 8. year snapping -------------------------------------------------------
cat("\n[8] rows snap to the nearest crosswalk year\n")
dat_yr <- dat %>%
  filter(location_period_id == "lp_single") %>%
  slice(rep(1, 4)) %>%
  mutate(year_seas = c(2010L, 2016L, 2019L, 2023L))
oy <- apportion_composite_cases(dat_yr, xwalk)
expect(
  identical(oy$pop_year, c(2015L, 2015L, 2020L, 2020L)),
  "2010->2015, 2016->2015, 2019->2020, 2023->2020 (crosswalk has 2015, 2020)"
)
expect(
  near(oy$pop, c(10000, 10000, 11000, 11000)),
  "and the population follows the snapped year, not the requested one"
)
expect(
  identical(
    .nearest_crosswalk_year(c(2012L, 2017L, 2018L), c(2015L, 2020L)),
    c(2015L, 2015L, 2020L)
  ),
  "nearest means nearest: 2017 -> 2015 (gap 2 vs 3), 2018 -> 2020 (gap 2 vs 3)"
)
expect(
  identical(.nearest_crosswalk_year(2017L, c(2015L, 2019L)), 2015L),
  "an exact tie breaks DOWNWARD, deterministically — never platform-dependent"
)
expect(
  identical(.nearest_crosswalk_year(NA_integer_, c(2015L, 2020L)), 2015L) &&
    identical(.nearest_crosswalk_year(2020L, integer(0)), NA_integer_),
  "a missing year takes the earliest available; an empty crosswalk yields NA"
)

# --- 9. admin level selects the tiling -------------------------------------
cat("\n[9] the row's admin level selects the crosswalk tiling\n")
dat_l2 <- dat %>%
  filter(location_period_id == "lp_single") %>%
  mutate(admin_level = "ADM2")
o2 <- apportion_composite_cases(dat_l2, xwalk)
expect(
  nrow(o2) == 2L && setequal(o2$gid, c("AAA.1.1_1", "AAA.1.2_1")),
  "the same LP is a 2-way composite at ADM2 — cover is level-relative"
)
expect(
  near(sum(o2$suspected_cases), 500) && near(sum(o2$pop), 10000),
  "conservation holds at ADM2 too"
)
expect(
  near(sort(o2$apportion_w), c(0.4, 0.6)),
  "ADM2 shares are 6000/4000, not the ADM1 single-piece w = 1"
)
expect(
  nrow(apportion_composite_cases(
    dat_l2,
    xwalk,
    level_col = NULL,
    gadm_level = "ADM1"
  )) ==
    1L,
  "level_col = NULL falls back to the scalar gadm_level argument"
)

# --- 10. the occurrence contract: propagate, do not scale -------------------
cat("\n[10] cases_col = NA propagates the row unchanged to every member\n")
oc <- apportion_composite_cases(dat, xwalk, cases_col = NA_character_)
expect(
  near(oc$suspected_cases[oc$location_period_id == "lp_comp7"], rep(7000, 7)),
  "with no cases_col the value is REPLICATED, not divided — the 0/1 outbreak indicator contract"
)
expect(
  sum(oc$suspected_cases[oc$location_period_id == "lp_comp7"]) == 7 * 7000,
  "which deliberately does NOT conserve a sum: `any()` semantics, not mass"
)

# --- 11. failure modes are loud --------------------------------------------
cat("\n[11] silent-corruption guards fire\n")
expect(
  inherits(
    try(apportion_composite_cases(dat, xwalk %>% select(-H), ), silent = TRUE),
    "try-error"
  ),
  "a crosswalk missing H is an error, not a silent no-op"
)
expect(
  inherits(
    try(apportion_composite_cases(dat, xwalk, lp_col = "nope"), silent = TRUE),
    "try-error"
  ),
  "an absent lp_col is an error, not zero matches"
)
expect(
  nrow(apportion_composite_cases(dat[0, ], xwalk)) == 0L,
  "an empty input returns empty, without erroring"
)

zero_H <- xwalk %>% mutate(H = ifelse(gid == "AAA.2_1", 0, H))
oz <- apportion_composite_cases(dat, zero_H)
expect(
  sum(oz$location_period_id == "lp_comp7") == 6L &&
    near(sum(oz$suspected_cases[oz$location_period_id == "lp_comp7"]), 7000),
  "a zero-mass piece is dropped AND w renormalised — cases still conserved over 6"
)


# --- 12. the SEAM: apportion -> collapse -> reconcile ------------------------
# Sections 1-11 test apportion_composite_cases() in isolation, and
# test_footprint_pop_reconciliation.R tests reconcile_footprint_pop() in
# isolation. Neither covers the COMPOSITION, which is what
# analysis/03_prepare_siz.R actually executes: apportion at the ob_ts3 seam,
# then collapse per (gid, year_seas, admin_level) with sum() over BOTH columns,
# then reconcile the population. The interesting claim is a joint one --
# `pop = sum(pop)` inside that collapse must land on exactly sum_l H[l,g], the
# quantity reconcile_footprint_pop() expects as H_sum. If apportionment ever
# emitted P_l instead of H, or the collapse grouped differently, both functions
# would still pass their own tests and the pipeline would silently double-count.
cat("\n[12] the pipeline seam: apportion, collapse, reconcile\n")

# One gid reported TWICE in whole by two different LPs (the NGA_19_1 defect),
# plus one composite LP spanning two other gids.
seam_xw <- bind_rows(
  mk_xw("whole1", "ZZZ.1_1", 10000, 2020),
  mk_xw("whole2", "ZZZ.1_1", 10000, 2020),
  mk_xw("comp", c("ZZZ.2_1", "ZZZ.3_1"), c(3000, 7000), 2020)
)
seam_gp <- tibble(
  gid = c("ZZZ.1_1", "ZZZ.2_1", "ZZZ.3_1"),
  year = 2020L,
  Q = c(10000, 3000, 7000),
  iso3 = "ZZZ",
  gadm_level = "ADM1"
)
seam_in <- tibble(
  location_period_id = c("whole1", "whole2", "comp"),
  gid = c("ZZZ.1_1", "ZZZ.1_1", "ZZZ.2_1"), # centroid guesses
  admin_level = "ADM1",
  year_seas = 2020L,
  suspected_cases = c(100, 50, 500),
  pop = c(10000, 10000, 10000)
)

seam_ap <- apportion_composite_cases(seam_in, seam_xw)

# The collapse, verbatim from 03_prepare_siz.R.
seam_col <- seam_ap %>%
  group_by(gid, year_seas, admin_level) %>%
  summarise(
    suspected_cases = sum(suspected_cases, na.rm = TRUE),
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )

seam_rc <- reconcile_footprint_pop(
  seam_col$gid,
  seam_col$year_seas,
  H_sum = seam_col$pop,
  gid_pop = seam_gp
)

expect(
  near(sum(seam_col$suspected_cases), 650),
  "cases survive apportion + collapse unchanged (100 + 50 + 500)"
)
expect(
  near(seam_col$pop[seam_col$gid == "ZZZ.1_1"], 20000),
  "the raw collapse DOES double the duplicated gid — 10000 + 10000"
)
expect(
  near(seam_rc$pop[seam_rc$gid == "ZZZ.1_1"], 10000) &&
    seam_rc$pop_tier[seam_rc$gid == "ZZZ.1_1"] == "gid_capped",
  "...and reconciliation saturates it back to Q = 10000, tagged 'gid_capped'"
)
expect(
  near(sort(seam_rc$pop[seam_rc$gid != "ZZZ.1_1"]), c(3000, 7000)) &&
    all(seam_rc$pop_tier[seam_rc$gid != "ZZZ.1_1"] == "footprint"),
  "the composite's members keep their own H and are NOT capped"
)

# The attack rate the model actually sees, after the whole chain.
seam_ar <- seam_col$suspected_cases / seam_rc$pop
expect(
  all(is.finite(seam_ar)) && all(seam_ar > 0) && all(seam_ar < 1),
  "every post-chain attack rate is finite, positive and below 1"
)
expect(
  near(seam_ar[seam_rc$gid == "ZZZ.2_1"], seam_ar[seam_rc$gid == "ZZZ.3_1"]),
  "the composite's two members share one attack rate — w cancels exactly"
)
expect(
  near(seam_ar[seam_rc$gid == "ZZZ.1_1"], 150 / 10000),
  "the duplicated gid's rate uses the SATURATED denominator, not 150/20000"
)

# The failure this whole exercise exists to prevent, stated as a number.
expect(
  150 / 20000 < seam_ar[seam_rc$gid == "ZZZ.1_1"],
  "(for contrast: sum(pop) would have halved this unit's apparent incidence)"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
