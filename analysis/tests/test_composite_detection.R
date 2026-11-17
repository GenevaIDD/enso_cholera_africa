# test_composite_detection.R — HARD: composite-LP detection algebra, isolated
#   from geometry.
#
# Usage:
#   Rscript analysis/tests/test_composite_detection.R
#
# Details:
#   Composite status is a property of the crosswalk (mass spread over >1 gid),
#   never of a name string (the upstream `composite_loc` column is backfilled
#   FALSE everywhere and detects nothing). Exercises compute_gid_cover() on
#   small synthetic (lctn_pr, gadm_level, year, gid, H) tables: no sf, raster,
#   network or data/. test_population_crosswalk.R drives the same derivation
#   through the full geometry pipeline.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/population_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-8) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}
row_for <- function(tbl, lp) tbl[tbl$lctn_pr == lp, , drop = FALSE]

# --- A. Singleton LP: unaffiliated with any other gid, maps entirely onto ----
#        itself (gid_cover == 1, n_gid_eff == 1, never composite) -----------
cat("## A) singleton LP -- unaffiliated gid, maps to itself\n")
{
  xw <- tibble(
    lctn_pr = "lp_single",
    gadm_level = "ADM1",
    year = 2020L,
    gid = "AAA.1_1",
    H = 1000
  )
  out <- compute_gid_cover(xw)
  r <- row_for(out, "lp_single")

  expect(
    nrow(out) == 1L,
    "one summary row for one (lctn_pr, gadm_level, year) group"
  )
  expect(
    near(r$gid_cover, 1),
    "singleton gid_cover == 1 (all mass in the one gid it maps to)"
  )
  expect(
    near(r$n_gid_eff, 1),
    "singleton n_gid_eff == 1 (inverse-Simpson of a point mass)"
  )
  expect(r$n_gid == 1L, "singleton n_gid == 1")
  expect(!isTRUE(r$is_composite), "singleton is never flagged composite")
}

# --- B. K-way composites, equal split ---------------------------------------
cat("\n## B) K-way composites (equal split)\n")
{
  mk_equal <- function(lp, k) {
    tibble(
      lctn_pr = lp,
      gadm_level = "ADM1",
      year = 2020L,
      gid = sprintf("AAA.%d_1", seq_len(k)),
      H = 100
    )
  }
  xw <- bind_rows(mk_equal("lp_2way", 2), mk_equal("lp_7way", 7))
  out <- compute_gid_cover(xw)

  r2 <- row_for(out, "lp_2way")
  r7 <- row_for(out, "lp_7way")

  expect(near(r2$gid_cover, 0.5), "2-way equal split -> gid_cover == 1/2")
  expect(near(r2$n_gid_eff, 2), "2-way equal split -> n_gid_eff == 2")
  expect(
    isTRUE(r2$is_composite),
    "2-way equal split flagged composite at default threshold 0.90"
  )

  expect(near(r7$gid_cover, 1 / 7), "7-way equal split -> gid_cover ~= 1/7")
  expect(
    near(r7$n_gid_eff, 7),
    "7-way equal split -> n_gid_eff ~= 7 (matches SDN_3_1-style composite)"
  )
  expect(isTRUE(r7$is_composite), "7-way equal split flagged composite")
}

# --- C. cover_threshold boundary behaviour ----------------------------------
cat("\n## C) cover_threshold boundary\n")
{
  xw <- bind_rows(
    tibble(
      lctn_pr = "lp_at_thr",
      gadm_level = "ADM1",
      year = 2020L,
      gid = c("AAA.1_1", "AAA.2_1"),
      H = c(90, 10)
    ), # gid_cover == 0.90 exactly
    tibble(
      lctn_pr = "lp_below_thr",
      gadm_level = "ADM1",
      year = 2020L,
      gid = c("AAA.1_1", "AAA.2_1"),
      H = c(89, 11)
    ) # gid_cover == 0.89
  )
  out <- compute_gid_cover(xw, cover_threshold = 0.90)
  r_at <- row_for(out, "lp_at_thr")
  r_below <- row_for(out, "lp_below_thr")

  expect(near(r_at$gid_cover, 0.90), "boundary LP gid_cover == 0.90 exactly")
  expect(
    !isTRUE(r_at$is_composite),
    "gid_cover == cover_threshold is NOT composite (strict '<', boundary inclusive on the clean side)"
  )
  expect(
    near(r_below$gid_cover, 0.89),
    "just-below-threshold LP gid_cover == 0.89"
  )
  expect(
    isTRUE(r_below$is_composite),
    "gid_cover just under cover_threshold IS flagged composite"
  )

  # A stricter caller-supplied threshold reclassifies the same LP.
  out_strict <- compute_gid_cover(xw, cover_threshold = 0.95)
  expect(
    isTRUE(row_for(out_strict, "lp_at_thr")$is_composite),
    "the same 0.90-cover LP becomes composite under a stricter cover_threshold = 0.95"
  )
}

# --- D. zero-mass group: gid_cover NA, never composite by default -----------
cat("\n## D) zero-mass group\n")
{
  xw <- tibble(
    lctn_pr = "lp_zero",
    gadm_level = "ADM1",
    year = 2020L,
    gid = c("AAA.1_1", "AAA.2_1"),
    H = c(0, 0)
  )
  out <- compute_gid_cover(xw)
  r <- row_for(out, "lp_zero")

  expect(
    is.na(r$gid_cover),
    "zero total mass -> gid_cover is NA, not 0/undefined-composite"
  )
  expect(is.na(r$n_gid_eff), "zero total mass -> n_gid_eff is NA")
  expect(
    !isTRUE(r$is_composite),
    "NA gid_cover -> is_composite FALSE, not NA (guarded by is.na())"
  )
}

# --- E. grouping correctness: independent (lctn_pr, gadm_level, year) -------
cat("\n## E) grouping is independent per (lctn_pr, gadm_level, year)\n")
{
  xw <- bind_rows(
    # same lctn_pr, two gadm_level tilings of the same underlying LP: ADM1
    # sees it as clean (cover=1), ADM2 sees the finer split (composite).
    tibble(
      lctn_pr = "lp_mixed",
      gadm_level = "ADM1",
      year = 2020L,
      gid = "AAA.1_1",
      H = 1000
    ),
    tibble(
      lctn_pr = "lp_mixed",
      gadm_level = "ADM2",
      year = 2020L,
      gid = c("AAA.1.1_1", "AAA.1.2_1"),
      H = c(600, 400)
    ),
    # same lctn_pr + gadm_level, a later year where it becomes a clean singleton
    tibble(
      lctn_pr = "lp_mixed",
      gadm_level = "ADM1",
      year = 2021L,
      gid = "AAA.1_1",
      H = 1200
    )
  )
  out <- compute_gid_cover(xw)

  expect(
    nrow(out) == 3L,
    "three independent groups: ADM1/2020, ADM2/2020, ADM1/2021"
  )
  r_adm1_2020 <- out %>%
    filter(lctn_pr == "lp_mixed", gadm_level == "ADM1", year == 2020L)
  r_adm2_2020 <- out %>%
    filter(lctn_pr == "lp_mixed", gadm_level == "ADM2", year == 2020L)
  expect(
    near(r_adm1_2020$gid_cover, 1) && !isTRUE(r_adm1_2020$is_composite),
    "ADM1 tiling of lp_mixed/2020 is clean (cover=1)"
  )
  expect(
    near(r_adm2_2020$gid_cover, 0.6) && isTRUE(r_adm2_2020$is_composite),
    "ADM2 tiling of the SAME lp_mixed/2020 is composite (cover=0.6) -- tiling-dependent, correctly independent"
  )
}

# --- F. empty crosswalk: correct skeleton, no error -------------------------
cat("\n## F) empty crosswalk\n")
{
  xw_empty <- tibble(
    lctn_pr = character(0),
    gadm_level = character(0),
    year = integer(0),
    H = numeric(0),
    gid = character(0)
  )
  out <- compute_gid_cover(xw_empty)
  expect(nrow(out) == 0L, "empty crosswalk -> zero-row result, not an error")
  expect(
    all(
      c(
        "lctn_pr",
        "gadm_level",
        "year",
        "P_l",
        "n_gid",
        "gid_cover",
        "n_gid_eff",
        "is_composite"
      ) %in%
        names(out)
    ),
    "empty result still has the full expected column set"
  )
}

# --- G. required-columns guard ----------------------------------------------
cat("\n## G) required-columns guard\n")
{
  ok <- tryCatch(
    {
      compute_gid_cover(tibble(lctn_pr = "x")) # missing gadm_level, H
      FALSE
    },
    error = function(e) TRUE
  )
  expect(
    ok,
    "missing required columns (gadm_level, H) -> stopifnot() error, not a silent NA result"
  )
}

# --- H. derive_lp_tiling_level(): the tiling an LP is READ at -------------------
# Section E above is the bug in miniature: the SAME LP is clean at ADM1 and
# composite at ADM2, so "is this LP composite?" is unanswerable until the
# tiling is fixed. It used to be fixed from the outbreak corpus's `admin_level`
# string, whose admin depth is not GADM's, which manufactured 256 false
# composites out of 1410 corpus LPs. derive_lp_tiling_level() fixes it from the
# LP's own geometry instead.
cat(
  "\n## H) derive_lp_tiling_level -- tiling depth from geometry, not the corpus label\n"
)

lvl_for <- function(tbl, lp) tbl[tbl$lctn_pr == lp, , drop = FALSE]

# H1. The Malawi pattern — THE regression test for this bug.
# lctn_pr 12753 (AFR::MWI::Central::Lilongwe) is labelled "ADM2" by the corpus,
# but it is a district = GADM ADM1. Read at GADM ADM2 it decomposes into its 25
# Traditional Authorities and gid_cover reads 0.394 — Lilongwe City TA's share
# of its own district, not a statement about compositeness at all.
{
  xw <- bind_rows(
    tibble(
      lctn_pr = "lp_mwi_district",
      gadm_level = "ADM1",
      year = 2020L,
      gid = "MWI.11_1",
      H = 2585480,
      piece_area_km2 = 6159
    ),
    tibble(
      lctn_pr = "lp_mwi_district",
      gadm_level = "ADM2",
      year = 2020L,
      gid = sprintf("MWI.11.%d_1", 1:25),
      # one dominant TA at 39.4%, the rest sharing the remainder
      H = c(0.394, rep((1 - 0.394) / 24, 24)) * 2585480,
      piece_area_km2 = c(456, rep((6159 - 456) / 24, 24))
    )
  )
  out <- derive_lp_tiling_level(xw, cover_threshold = 0.90)
  r <- lvl_for(out, "lp_mwi_district")

  expect(nrow(out) == 1L, "one decision row per LP, not one per (level, year)")
  expect(
    identical(r$lp_level, "ADM1"),
    "a district shattered into 25 TAs at GADM ADM2 is read at GADM ADM1"
  )
  expect(
    near(r$pcov_ADM2, 0.394, tol = 1e-3),
    "the ADM2 reading reproduces the reported gid_cover = 0.394"
  )
  expect(
    near(r$pcov_ADM1, 1),
    "the ADM1 reading owns the whole LP (pcov = 1) -- NOT composite"
  )
  expect(
    identical(r$lp_level_reason, "only_adm1_admissible"),
    "reason records that ADM2 failed the dominance test, ADM1 passed"
  )
}

# H2. Population concentration must NOT be read as a finer tiling.
# A province whose capital district holds 92% of its people clears the
# dominance threshold at both levels. Refining to ADM2 would be wrong: it is
# one province, and 92% is a demographic fact, not a boundary fact. The area
# tie-break decides — this is the fixture that rejects "finest admissible".
{
  xw <- bind_rows(
    tibble(
      lctn_pr = "lp_concentrated",
      gadm_level = "ADM1",
      year = 2020L,
      gid = "XXX.1_1",
      H = 1000,
      piece_area_km2 = 100
    ),
    tibble(
      lctn_pr = "lp_concentrated",
      gadm_level = "ADM2",
      year = 2020L,
      gid = c("XXX.1.1_1", "XXX.1.2_1"),
      H = c(920, 80),
      piece_area_km2 = c(20, 80)
    )
  )
  r <- lvl_for(
    derive_lp_tiling_level(xw, cover_threshold = 0.90),
    "lp_concentrated"
  )

  expect(
    identical(r$lp_level, "ADM1"),
    "both levels admissible -> area tie-break keeps the province at ADM1"
  )
  expect(
    identical(r$lp_level_reason, "both_admissible_area_tiebreak"),
    "reason records that the area tie-break, not dominance, made the call"
  )
}

# H3. A genuine sub-ADM1 LP (a city) stays at the finer tiling.
# Its geometry sits inside one unit at both levels, so both acov == 1 and the
# tie goes to ADM2 — the guard against over-coarsening every city to its
# province, which is what "coarsest covering level" would have done.
{
  xw <- bind_rows(
    tibble(
      lctn_pr = "lp_city",
      gadm_level = "ADM1",
      year = 2020L,
      gid = "XXX.2_1",
      H = 500000,
      piece_area_km2 = 40
    ),
    tibble(
      lctn_pr = "lp_city",
      gadm_level = "ADM2",
      year = 2020L,
      gid = "XXX.2.3_1",
      H = 500000,
      piece_area_km2 = 40
    )
  )
  r <- lvl_for(derive_lp_tiling_level(xw), "lp_city")
  expect(
    identical(r$lp_level, "ADM2"),
    "an LP wholly inside one unit at both levels keeps the FINER tiling"
  )
}

# H4. A genuine multi-ADM1 aggregate is still composite (the SDN 7-state
# signature). No level is admissible, so it is read at the coarsest and the
# dominance test downstream correctly flags it.
{
  xw <- bind_rows(
    tibble(
      lctn_pr = "lp_sdn7",
      gadm_level = "ADM1",
      year = 2020L,
      gid = sprintf("SDN.%d_1", 1:7),
      H = 100,
      piece_area_km2 = 500
    ),
    tibble(
      lctn_pr = "lp_sdn7",
      gadm_level = "ADM2",
      year = 2020L,
      gid = sprintf("SDN.1.%d_1", 1:30),
      H = 700 / 30,
      piece_area_km2 = 3500 / 30
    )
  )
  r <- lvl_for(derive_lp_tiling_level(xw), "lp_sdn7")
  expect(
    identical(r$lp_level, "ADM1"),
    "a true multi-unit aggregate is read at ADM1, where dominance is maximal"
  )
  expect(
    identical(r$lp_level_reason, "none_admissible_composite"),
    "no level owns the LP -> flagged as a genuine composite, not rescued"
  )
  expect(
    near(r$pcov_ADM1, 1 / 7),
    "its ADM1 cover is 1/7 -- far under any sane threshold"
  )
}

# H5. Year stability. Population is redistributed across the 2015-2030 span, so
# a per-year rule flips the level for some LPs. Reading one LP at two tilings
# in two years would double-count its pieces, so the decision is taken once per
# LP from the MEDIAN cover.
{
  mk_year <- function(y, cov1) {
    bind_rows(
      tibble(
        lctn_pr = "lp_straddle",
        gadm_level = "ADM1",
        year = y,
        gid = c("YYY.1_1", "YYY.2_1"),
        H = c(cov1, 1 - cov1) * 1000,
        piece_area_km2 = c(80, 20)
      ),
      tibble(
        lctn_pr = "lp_straddle",
        gadm_level = "ADM2",
        year = y,
        gid = c("YYY.1.1_1", "YYY.1.2_1"),
        H = c(0.95, 0.05) * 1000,
        piece_area_km2 = c(60, 40)
      )
    )
  }
  # ADM1 cover straddles the 0.90 threshold year to year; its MEDIAN (0.89) does not.
  xw <- bind_rows(
    mk_year(2019L, 0.88),
    mk_year(2020L, 0.89),
    mk_year(2021L, 0.93)
  )
  out <- derive_lp_tiling_level(xw, cover_threshold = 0.90)

  expect(
    nrow(out) == 1L,
    "an LP whose per-year cover straddles the threshold still gets ONE level"
  )
  expect(
    near(out$pcov_ADM1, 0.89),
    "the level is decided on the MEDIAN year, not the first or the best one"
  )
  expect(
    identical(out$lp_level, "ADM2"),
    "median ADM1 cover 0.89 < 0.90 is inadmissible; ADM2 at 0.95 is chosen"
  )
}

# H6. Guards.
{
  ok <- tryCatch(
    {
      derive_lp_tiling_level(tibble(lctn_pr = "x"))
      FALSE
    },
    error = function(e) TRUE
  )
  expect(ok, "missing required columns -> error, not a silent NA level")
  ok2 <- tryCatch(
    {
      derive_lp_tiling_level(
        tibble(lctn_pr = "x", gid = "g", gadm_level = "ADM1", H = 1),
        cover_threshold = 1.5
      )
      FALSE
    },
    error = function(e) TRUE
  )
  expect(ok2, "cover_threshold outside [0, 1] is an error -- it is a share")
  empty <- derive_lp_tiling_level(tibble(
    lctn_pr = character(0),
    gid = character(0),
    gadm_level = character(0),
    H = numeric(0)
  ))
  expect(
    nrow(empty) == 0L && "lp_level" %in% names(empty),
    "empty crosswalk -> zero-row result with the full column set"
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
