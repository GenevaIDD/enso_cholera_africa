# test_composite_policy.R — HARD: composite location-period handling policy
#   (`--composite_policy`).
#
# Usage:
#   Rscript analysis/tests/test_composite_policy.R
#
# Details:
#   apportion_composite_cases() expands one location-period observation into K
#   member-gid rows, which the likelihood treats as K independent data points.
#   `--composite_policy` is the fix; this file pins its four policies (listed
#   above the first section) and the dominance-vs-multiplicity definition of
#   "composite" (fixture note). Pure synthetic tables: no sf, raster, SQLite,
#   network or data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/population.R") # apportion_composite_cases() -- used directly
source("analysis/utils/composite_lp.R") # collapse_composite_lps()/lp_composite_status() -- used directly

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-9) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}

# THE ONE THING MOST LIKELY TO GO WRONG, and the reason for section [2]: two
# non-equivalent definitions of "composite" coexist:
#
#   is_composite  DOMINANCE     gid_cover < 0.90        315/1409 LPs   (22.4%)
#   apportioned   MULTIPLICITY  n_members > 1          8779/9036 rows  (97.0%)
#
# `drop` keys on DOMINANCE. Keying it on MULTIPLICITY instead would discard ~97%
# of event rows and would not error.
# ---------------------------------------------------------------------------
# Fixture — the three cases that separate DOMINANCE from MULTIPLICITY
#
#   lp_true_comp  w = .40/.35/.25  gid_cover = .40  COMPOSITE, multi-member
#   lp_dominant   w = .97/.03      gid_cover = .97  not composite, MULTI-MEMBER
#   lp_single     w = 1            gid_cover = 1    not composite, single
#   lp_absent     not in crosswalk                  passthrough
#
# lp_dominant is the whole point: it is the ~88% of row inflation that `drop`
# alone leaves behind and only `collapse` removes.
#
# lp_boundary sits at gid_cover EXACTLY 0.90 — `<` means it is NOT composite.
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
  # deliberately NOT in descending H order, so "argmax H" cannot be faked by
  # taking the first matched row
  mk_xw(
    "lp_true_comp",
    c("AAA.1_1", "AAA.2_1", "AAA.3_1"),
    c(3500, 4000, 2500),
    2015
  ), # cover .40 on AAA.2_1
  mk_xw("lp_dominant", c("BBB.1_1", "BBB.2_1"), c(300, 9700), 2015), # cover .97 on BBB.2_1
  mk_xw("lp_single", "CCC.1_1", 10000, 2015),
  mk_xw("lp_boundary", c("DDD.1_1", "DDD.2_1"), c(9000, 1000), 2015) # cover EXACTLY .90
)

dat <- tibble(
  location_period_id = c(
    "lp_true_comp",
    "lp_dominant",
    "lp_single",
    "lp_boundary",
    "lp_absent"
  ),
  gid = c("AAA.3_1", "BBB.1_1", "CCC.1_1", "DDD.2_1", "ZZZ.9_1"),
  admin_level = "ADM1",
  year_seas = 2015L,
  suspected_cases = c(1000, 2000, 500, 300, 42),
  outbreak = 1L,
  pop = c(10000, 10000, 10000, 10000, 999)
)

# Why a policy: after the expansion the likelihood treats the K member rows as
# K independent data points.
#   occurrence: the 0/1 indicator is PROPAGATED, so "at least one member had an
#     outbreak" is fitted as "all K members had one simultaneously";
#   size: cases and pop are scaled by the same w, so all K siblings carry an
#     IDENTICAL implied attack rate Y_l/P_l and row-wise PSIS-LOO holds out a
#     row whose siblings already encode it. ELPD optimism scales with K.
# K is geographically clustered (KEN/COD/MWI/SDN), so both distortions are
# region-dependent, which is why this reaches the ENSO regional best-lag
# results, selected by LOO. The four policies:
#   keep          expand to K rows (regression baseline)
#   collapse      one row per LP-event at the dominant member, nothing dropped
#   drop          remove TRUE composites, expand the rest
#   drop_collapse remove TRUE composites, collapse the rest        [the default]
cat("# composite location-period policy\n")

# ===========================================================================
cat("\n[1] lp_composite_status: one row per input row, dominance criterion\n")

st <- lp_composite_status(dat, xwalk)

expect(
  nrow(st) == nrow(dat) && identical(st$.row, seq_len(nrow(dat))),
  "one status row per input row, in input order"
)
expect(
  near(st$gid_cover[st$lctn_pr == "lp_true_comp"], 0.40) &&
    near(st$gid_cover[st$lctn_pr == "lp_dominant"], 0.97) &&
    near(st$gid_cover[st$lctn_pr == "lp_single"], 1.00),
  "gid_cover == max_g w, not 1/K and not the input row's own share"
)
expect(
  identical(st$is_composite, c(TRUE, FALSE, FALSE, FALSE, FALSE)),
  "ONLY lp_true_comp is composite — lp_dominant is multi-member but NOT"
)
expect(
  identical(st$n_member, c(3L, 2L, 1L, 2L, 1L)),
  "n_member is multiplicity and is reported separately from dominance"
)
expect(
  near(
    st$n_gid_eff[st$lctn_pr == "lp_true_comp"],
    1 / (0.40^2 + 0.35^2 + 0.25^2)
  ),
  "n_gid_eff is inverse-Simpson over the surviving pieces"
)
expect(
  isFALSE(st$matched[st$lctn_pr == "lp_absent"]) &&
    isFALSE(st$is_composite[st$lctn_pr == "lp_absent"]),
  "an LP the crosswalk never saw is NOT declared composite — never drop on an absent measurement"
)

# ===========================================================================
cat("\n[2] the threshold is `<`, evaluated at run time\n")

expect(
  isFALSE(lp_composite_status(dat, xwalk, cover_threshold = 0.90)$is_composite[
    st$lctn_pr == "lp_boundary"
  ]),
  "gid_cover EXACTLY 0.90 is NOT composite — the test is strict `<`"
)
expect(
  isTRUE(lp_composite_status(dat, xwalk, cover_threshold = 0.95)$is_composite[
    st$lctn_pr == "lp_boundary"
  ]),
  "...and raising the threshold to 0.95 flips it, so the knob is live"
)
expect(
  sum(lp_composite_status(dat, xwalk, cover_threshold = 0.99)$is_composite) ==
    3L,
  "threshold 0.99 catches true_comp + dominant + boundary, never single"
)
expect(
  sum(lp_composite_status(dat, xwalk, cover_threshold = 0.0)$is_composite) ==
    0L,
  "threshold 0 catches nothing — `<` cannot fire on a share"
)

# ===========================================================================
cat("\n[3] collapse: at most one row per input row, NEVER more\n")

col_all <- collapse_composite_lps(dat, xwalk, drop_composites = FALSE)
rep_all <- attr(col_all, "collapse_report")

expect(
  nrow(col_all) == nrow(dat),
  "drop_composites = FALSE keeps every input row — 5 in, 5 out"
)
expect(
  nrow(col_all) <= nrow(dat),
  "the invariant: output rows never exceed input rows"
)
expect(
  rep_all$n_rows_if_keep == 3L + 2L + 1L + 2L + 1L,
  "the report carries the `keep` counterfactual (9 rows) without a second run"
)
expect(
  all(table(col_all$location_period_id) == 1L),
  "exactly one row per retained LP-event"
)

# ===========================================================================
cat("\n[4] the dominant member is chosen by max H, not by input order\n")

expect(
  col_all$gid[col_all$location_period_id == "lp_true_comp"] == "AAA.2_1",
  "lp_true_comp collapses onto AAA.2_1 (H = 4000), not the first-listed AAA.1_1"
)
expect(
  col_all$gid[col_all$location_period_id == "lp_dominant"] == "BBB.2_1",
  "lp_dominant collapses onto BBB.2_1 (H = 9700), not the incoming centroid BBB.1_1"
)
expect(
  col_all$gid_centroid[col_all$location_period_id == "lp_dominant"] ==
    "BBB.1_1",
  "...and the incoming centroid gid is preserved as gid_centroid"
)
tie_xw <- mk_xw("lp_tie", c("TTT.2_1", "TTT.1_1"), c(5000, 5000), 2015)
tie_in <- tibble(
  location_period_id = "lp_tie",
  gid = "TTT.2_1",
  admin_level = "ADM1",
  year_seas = 2015L,
  suspected_cases = 100,
  pop = 10000
)
expect(
  collapse_composite_lps(tie_in, tie_xw, drop_composites = FALSE)$gid ==
    "TTT.1_1",
  "an exact H tie breaks on gid string order — deterministic, not join-order"
)

# ===========================================================================
cat("\n[5] NOTHING is scaled — the whole point of collapse vs apportion\n")

expect(
  near(sort(col_all$suspected_cases), sort(dat$suspected_cases)),
  "case counts pass through untouched — no multiplication by w"
)
expect(
  all(col_all$outbreak == 1L),
  "the occurrence indicator passes through untouched"
)
expect(
  near(sum(col_all$suspected_cases), sum(dat$suspected_cases)),
  "so total case mass is conserved trivially, not by a rounding accident"
)

# ===========================================================================
cat("\n[6] the denominator is P_l, so the attack rate stays Y_l / P_l\n")

expect(
  near(col_all$pop[col_all$location_period_id == "lp_true_comp"], 10000),
  "pop == P_l = sum_g H (10000), NOT the dominant member's H = 4000"
)
expect(
  near(col_all$pop[col_all$location_period_id == "lp_dominant"], 10000),
  "...same for lp_dominant: 10000, not 9700"
)
ap <- apportion_composite_cases(dat, xwalk)
ar_ap <- ap$suspected_cases[ap$location_period_id == "lp_true_comp"] /
  ap$pop[ap$location_period_id == "lp_true_comp"]
ar_col <- col_all$suspected_cases[
  col_all$location_period_id == "lp_true_comp"
] /
  col_all$pop[col_all$location_period_id == "lp_true_comp"]
expect(
  near(ar_col, unique(round(ar_ap, 12))) &&
    length(unique(round(ar_ap, 12))) == 1L,
  "the collapsed row's attack rate EQUALS every apportioned sibling's — w cancels"
)
expect(
  near(col_all$pop_lp[col_all$location_period_id == "lp_true_comp"], 10000),
  "the incoming pop is preserved as pop_lp, same convention as apportion"
)

# ===========================================================================
cat("\n[7] schema parity with the expansion path\n")

expect(
  all(
    c(
      "gid_centroid",
      "pop_lp",
      "pop_tier",
      "pop_year",
      "apportion_w",
      "n_apportion_members",
      "apportioned"
    ) %in%
      names(col_all)
  ),
  "every column apportion_composite_cases() emits is present"
)
expect(
  near(col_all$apportion_w, rep(1, nrow(col_all))) &&
    all(col_all$n_apportion_members == 1L) &&
    !any(col_all$apportioned),
  "apportion_w == 1, n_apportion_members == 1, apportioned == FALSE"
)
expect(
  all(
    col_all$pop_tier[col_all$location_period_id != "lp_absent"] ==
      "lp_collapsed"
  ) &&
    col_all$pop_tier[col_all$location_period_id == "lp_absent"] ==
      "unmatched_lp",
  "pop_tier distinguishes a collapsed row from an unmatched passthrough"
)
expect(
  identical(col_all$n_member_orig, c(3L, 2L, 1L, 2L, 1L)),
  "n_member_orig records the multiplicity that WOULD have been created"
)
expect(
  identical(col_all$collapsed, c(TRUE, TRUE, TRUE, TRUE, FALSE)),
  "`collapsed` flags the reduced rows; the unmatched passthrough is not one"
)

# ===========================================================================
cat("\n[8] drop_collapse: true composites removed, the rest collapsed\n")

dc <- collapse_composite_lps(dat, xwalk, drop_composites = TRUE)
rdc <- attr(dc, "collapse_report")

expect(
  nrow(dc) == 4L && !("lp_true_comp" %in% dc$location_period_id),
  "lp_true_comp is dropped; the other four survive"
)
expect(
  "lp_dominant" %in% dc$location_period_id,
  "lp_dominant SURVIVES — dropping keys on dominance, not multiplicity"
)
expect(
  identical(attr(dc, "collapse_dropped_lp"), "lp_true_comp"),
  "the dropped LP ids are returned, so callers can restate their gates"
)
expect(
  rdc$n_lp_dropped == 1L &&
    rdc$n_rows_dropped == 1L &&
    rdc$n_rows_in - rdc$n_rows_dropped == rdc$n_rows_out,
  "rows_out + rows_dropped == rows_in — the reporting block reconciles"
)
expect(
  near(sum(dc$suspected_cases), sum(dat$suspected_cases) - 1000),
  "dropped case mass equals exactly the dropped LP's cases (1000)"
)
expect(
  identical(
    attr(
      collapse_composite_lps(dat, xwalk, drop_composites = FALSE),
      "collapse_dropped_lp"
    ),
    character(0)
  ),
  "drop_composites = FALSE reports nothing dropped"
)

# ===========================================================================
cat("\n[9] drop (without collapse) = filter on dominance, then expand\n")

# This is what 03_*.R does for --composite_policy drop.
keep_rows <- !lp_composite_status(dat, xwalk)$is_composite
dr <- apportion_composite_cases(dat[keep_rows, , drop = FALSE], xwalk)

expect(
  nrow(dr) == 2L + 1L + 2L + 1L,
  "drop leaves 6 rows: the 3-way composite is gone, the 2-way fringe is NOT"
)
expect(
  nrow(dr) > nrow(dc),
  "drop still carries the multiplicity that drop_collapse removes"
)
expect(
  near(sum(dr$suspected_cases), sum(dat$suspected_cases) - 1000),
  "drop conserves cases among RETAINED LPs — same mass as drop_collapse"
)

# ===========================================================================
cat(
  "\n[10] keep is untouched — the .match_lp_members() extraction is behaviour-free\n"
)

expect(nrow(ap) == 3L + 2L + 1L + 2L + 1L, "keep still expands to 9 rows")
expect(
  near(sum(ap$suspected_cases), sum(dat$suspected_cases)),
  "keep still conserves cases across the expansion"
)
expect(
  all(
    c(
      "apportion_w",
      "n_apportion_members",
      "apportioned",
      "gid_centroid",
      "pop_lp",
      "pop_tier",
      "pop_year"
    ) %in%
      names(ap)
  ) &&
    !any(c("gid_cover", "collapsed", "n_member_orig") %in% names(ap)),
  "keep's schema is unchanged — collapse-only columns did not leak into it"
)

# ===========================================================================
cat("\n[11] silent-corruption guards fire\n")

expect(
  inherits(
    try(collapse_composite_lps(dat, xwalk %>% select(-H)), silent = TRUE),
    "try-error"
  ),
  "a crosswalk missing H is an error, not a silent no-op"
)
expect(
  inherits(
    try(collapse_composite_lps(dat, xwalk, lp_col = "nope"), silent = TRUE),
    "try-error"
  ),
  "an absent lp_col is an error, not zero matches"
)
expect(
  inherits(
    try(
      collapse_composite_lps(dat, xwalk, cover_threshold = 1.5),
      silent = TRUE
    ),
    "try-error"
  ),
  "a cover_threshold outside [0, 1] is an error — it is a share, not a count"
)
expect(
  nrow(collapse_composite_lps(dat[0, ], xwalk)) == 0L,
  "an empty input returns empty, without erroring"
)
expect(
  nrow(collapse_composite_lps(
    dat,
    xwalk,
    drop_composites = FALSE,
    keep_unmatched = FALSE
  )) ==
    4L,
  "keep_unmatched = FALSE drops the crosswalk-absent LP"
)
zero_xw <- xwalk %>% mutate(H = if_else(gid == "AAA.2_1", 0, H))
z <- collapse_composite_lps(dat, zero_xw, drop_composites = FALSE)
expect(
  z$gid[z$location_period_id == "lp_true_comp"] == "AAA.1_1" &&
    near(z$pop[z$location_period_id == "lp_true_comp"], 6000),
  "a zero-mass piece is excluded from BOTH the argmax and P_l"
)

# ===========================================================================
cat(
  "\n[12] tiling: the level an LP is READ at comes from geometry, not the corpus label\n"
)
#
# The bug all four policies inherited. The corpus's admin depth is NOT GADM's:
# Malawi's corpus "ADM2" is a district, which GADM calls ADM1 (GADM's ADM2 is
# the Traditional Authority). Reading that district at GADM ADM2 shatters it
# into its TAs, gid_cover reads the biggest TA's share of its own district, and
# `drop_collapse` deletes a perfectly ordinary district. 256 of 1410 corpus LPs
# were destroyed this way. The crosswalk now carries the geometry-derived
# `lp_level`; `.match_lp_members()` reads it, so ONE change covers both tracks and
# all four policies.
{
  # One district. GADM ADM1 = the district itself; GADM ADM2 = its 25 TAs.
  xw_mwi <- bind_rows(
    mk_xw("lp_mwi", "MWI.11_1", 2585480, 2015, level = "ADM1"),
    mk_xw(
      "lp_mwi",
      sprintf("MWI.11.%d_1", 1:25),
      c(0.394, rep((1 - 0.394) / 24, 24)) * 2585480,
      2015,
      level = "ADM2"
    )
  ) %>%
    mutate(lp_level = "ADM1") # what derive_lp_tiling_level() froze at 02c

  # The corpus calls it "ADM2" — the label that used to pick the tiling.
  dat_mwi <- tibble(
    location_period_id = "lp_mwi",
    gid = "MWI.11.1_1",
    admin_level = "ADM2",
    year_seas = 2015L,
    suspected_cases = 900,
    outbreak = 1L,
    pop = 2585480
  )

  st_derived <- lp_composite_status(dat_mwi, xw_mwi) # default
  st_label <- lp_composite_status(dat_mwi, xw_mwi, tiling = "corpus_label")

  expect(
    !isTRUE(st_derived$is_composite) && near(st_derived$gid_cover, 1),
    "derived tiling reads the district at GADM ADM1: cover 1, NOT composite"
  )
  expect(
    st_derived$n_member == 1L,
    "derived tiling sees ONE member, not the 25 TAs inside the district"
  )
  expect(
    isTRUE(st_label$is_composite) && near(st_label$gid_cover, 0.394, 1e-3),
    "corpus_label reproduces the bug exactly (cover 0.394, composite) -- side-by-side, not deleted"
  )

  # The consequence the whole fix exists for.
  keep_derived <- collapse_composite_lps(
    dat_mwi,
    xw_mwi,
    drop_composites = TRUE
  )
  keep_label <- collapse_composite_lps(
    dat_mwi,
    xw_mwi,
    drop_composites = TRUE,
    tiling = "corpus_label"
  )
  expect(
    nrow(keep_derived) == 1L && keep_derived$gid == "MWI.11_1",
    "drop_collapse RETAINS the district and puts it on its own GADM ADM1 gid"
  )
  expect(
    nrow(keep_label) == 0L,
    "under corpus_label the same district is dropped -- the 256-LP data loss"
  )
  expect(
    near(keep_derived$pop, 2585480),
    "the retained row carries the WHOLE district population, not one TA's"
  )

  # An LP the crosswalk has no lp_level for must not become unreadable.
  xw_na <- xw_mwi %>% mutate(lp_level = NA_character_)
  st_na <- lp_composite_status(dat_mwi, xw_na)
  expect(
    isTRUE(st_na$matched),
    "no derived level -> fall back to the corpus label, never zero matches"
  )

  # A crosswalk predating the fix has no column at all; behaviour is unchanged.
  st_old <- lp_composite_status(dat_mwi, xw_mwi %>% select(-lp_level))
  expect(
    isTRUE(st_old$is_composite),
    "a pre-fix crosswalk (no lp_level column) keeps the old behaviour exactly"
  )

  expect(
    inherits(
      try(lp_composite_status(dat_mwi, xw_mwi, tiling = "nope"), silent = TRUE),
      "try-error"
    ),
    "an unknown tiling value is an error, not a silent default"
  )
  expect(
    inherits(
      try(
        lp_composite_status(
          dat_mwi,
          xw_mwi %>%
            mutate(lp_level = if_else(gadm_level == "ADM1", "ADM1", "ADM2"))
        ),
        silent = TRUE
      ),
      "try-error"
    ),
    "two lp_level values for one LP is an error -- it must be one decision per LP"
  )
}

# ===========================================================================
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
