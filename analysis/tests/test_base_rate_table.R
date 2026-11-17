# test_base_rate_table.R — HARD: the shared four-track baseline-rate anchor
#   (analysis/utils/anchor.R).
#
# Usage:
#   Rscript analysis/tests/test_base_rate_table.R
#
# Details:
#   Covers build_base_rate_table(), resolve_base_rate(), attach_base_rate() and
#   resolve_base_rate_paths(): the POOLED estimator occ/dur/mon require (their
#   responses are 0/1 or carry structural zero months, so a per-row log(rate)
#   is -Inf), the builder on all four tracks, and the two invariants that make
#   the anchor safe to add to a fitted model (sections 6 and 9).
#   test_mean_rate_table.R pins the size-only GEOMETRIC estimator; the two are
#   complementary. Synthetic fixtures only, no real outbreak or climate data.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R") # transitive dependency of anchor.R (utils/README.md)
source("analysis/utils/config.R") # .OFFSET_SPEC / .offset_terms() -- gate-consistency section
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/paths.R") # transitive dependency of anchor.R (utils/README.md)
source("analysis/utils/provenance.R") # resolve_mean_rate()/build_mean_rate_table() -- delegate + geometric cross-check
source("analysis/utils/population.R") # transitive dependency of occ_data.R/siz_data.R (utils/README.md)
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/occ_data.R") # .occ_poisson_view() -- the occ view's numerator/denominator
source("analysis/utils/siz_data.R") # resolve_size_estimand()/resolve_size_pop_column() -- the siz view
source("analysis/utils/simulator_core.R") # .BUNDLE_ANCHOR_COL / .anchor_log_offset()
source("analysis/utils/anchor.R") # the unit under test

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R") # expect_true/expect_equal

ENSO_CSV <- "analysis/groupings/enso_grouping.csv"

cat(
  "# base_rate anchor table (build_base_rate_table / resolve_base_rate) tests\n\n"
)

# --------------------------------------------------------------------------
# Shared synthetic corpus: one (country, admin_level, n_rows, events) design
# for all four tracks; only column names differ (the builder is view-agnostic).
# Real ENSO-registry iso3 (KEN/UGA/TZA/RWA are enso_east) exercise T3:
#   KEN, UGA : ADM1 and ADM2 rows with events -> own T1 cells; the two
#              countries spanning >1 level that level_factor needs
#   TZA      : ADM1 rows only -> its ADM2 query resolves at T2
#   BIGN     : a large ADM1 cell -> the w -> 1 shrinkage check
#   ZERO     : exposure but zero events -> T1-zero (pooled estimator only)
#   RWA      : not in training, in the registry -> T3
#   ZZZ      : in neither -> T4
# --------------------------------------------------------------------------
CELLS <- tibble::tribble(
  ~country , ~admin_level , ~n_rows , ~rate ,
  "KEN"    , "ADM1"       ,   40L   , 0.010 ,
  "KEN"    , "ADM2"       ,   40L   , 0.024 ,
  "UGA"    , "ADM1"       ,   40L   , 0.008 ,
  "UGA"    , "ADM2"       ,   40L   , 0.020 ,
  "TZA"    , "ADM1"       ,   12L   , 0.012 ,
  "BIGN"   , "ADM1"       , 4000L   , 0.010 ,
  "ZERO"   , "ADM1"       ,   25L   , 0.000
)

# Bernoulli-ish event draws with a fixed seed, so every track sees the same
# (num, den) pattern and any cross-track difference is the builder's, not the
# fixture's.
make_corpus <- function(seed = 7L) {
  set.seed(seed)
  purrr_rows <- lapply(seq_len(nrow(CELLS)), function(i) {
    ci <- CELLS[i, ]
    den <- rep(30, ci$n_rows) # 30 exposure-days per row
    num <- stats::rbinom(ci$n_rows, size = 1L, prob = min(1, ci$rate * 30))
    tibble::tibble(
      country = ci$country,
      admin_level = ci$admin_level,
      num = as.numeric(num),
      den = as.numeric(den)
    )
  })
  dplyr::bind_rows(purrr_rows)
}
corpus <- make_corpus()

# Sanity: the fixture must actually produce the tiers the tests below name.
stopifnot(
  sum(corpus$num[corpus$country == "ZERO"]) == 0,
  sum(corpus$num[corpus$country == "TZA"]) > 0
)

# --- per-track materialisation of the same corpus --------------------------
#
# `den` is, for every track, exactly the product of that track's other offsets
# (that identity is the entire mechanism: mu = <other offsets> x anchor x
# exp(lp), so exp(Intercept) ~ 1). Each writer below therefore has to produce a
# panel whose view() reconstructs `den` from the panel's own columns.
write_track_panel <- function(track, d) {
  out <- switch(
    track,
    occ = d %>%
      dplyr::transmute(
        country,
        admin_level,
        outbreak = num,
        # .occ_poisson_view(): exposure = pmax(t_end_seas - t_start_seas, 1)
        t_start_seas = 0L,
        t_end_seas = as.integer(den),
        t_start = as.Date("2015-01-01"),
        t_end = as.Date("2015-02-01")
      ),
    dur = d %>%
      dplyr::transmute(
        country,
        admin_level,
        ended = num,
        exposure = den
      ),
    mon = d %>%
      dplyr::transmute(
        country,
        admin_level,
        cases_month = num,
        # den = pop * pmax(exposure, 1). `pop` carries the whole denominator only
        # if exposure clears the pmax() clamp, so put `den` in `exposure` (>= 1 by
        # construction) and leave pop at 1 -- the other split, exposure = den/1000,
        # is silently clamped to 1 and the view then reconstructs 1000, not `den`.
        pop = 1,
        exposure = den
      ),
    siz = d %>%
      dplyr::transmute(
        country,
        admin_level,
        # den = pop. num is 0/1 and den is 30, so the attack rate is 0.033 --
        # inside the 0.10 max_attack_rate guard and above the min_pop = 10 floor.
        suspected_cases = num,
        pop = den,
        TL = as.Date("2015-01-01"),
        TR = as.Date("2015-06-01")
      ),
    stop("write_track_panel: unknown track ", track)
  )
  p <- tempfile(fileext = ".rds")
  saveRDS(out, p)
  p
}

# `opt` per track: anchor ON, size guards relaxed so the fixture survives them,
# and size pinned to its NON-canonical admin level (the canonical-level branch
# needs a real gid population artefact, which a synthetic fixture has none of).
track_opt <- function(track) {
  base <- list(
    occurrence = list(offset_baserate = TRUE),
    duration = list(offset_baserate = TRUE),
    monthly = list(offset_baserate = TRUE),
    size = list(
      offset_meanrate = TRUE,
      offset_baserate = TRUE,
      min_pop = 10,
      max_attack_rate = 0.10,
      canonical_admin_level = FALSE
    )
  )
  base
}

TRACKS <- c("occ", "dur", "mon", "siz")
panels <- stats::setNames(
  lapply(TRACKS, function(tr) {
    write_track_panel(tr, corpus)
  }),
  TRACKS
)

build <- function(track, estimator = "pooled", opt = track_opt(track)) {
  spec <- anchor_spec(track)
  opt[[spec$opt_block]][[spec$cfg_key]] <-
    c(
      opt[[spec$opt_block]][[spec$cfg_key]] %||% list(),
      list(estimator = estimator)
    )
  suppressWarnings(build_base_rate_table(
    opt,
    track,
    path = panels[[track]],
    enso_grouping_path = ENSO_CSV
  ))
}

brt <- stats::setNames(lapply(TRACKS, build), TRACKS)


cat("## 1) build_base_rate_table: one schema, all four tracks\n")
{
  for (tr in TRACKS) {
    b <- brt[[tr]]
    expect_true(
      all(c("cell", "group", "global", "meta") %in% names(b)),
      sprintf("%s: returns list(cell, group, global, meta)", tr)
    )
    expect_true(
      all(
        c(
          "country",
          "admin_level",
          "mean_rate",
          "mr_tier",
          "n_cell",
          "w_shrink",
          "mr_raw",
          "mr_parent",
          "estimator",
          "k",
          "level_factor"
        ) %in%
          names(b$cell)
      ),
      sprintf(
        "%s: cell table has build_mean_rate_table()'s schema verbatim",
        tr
      )
    )
    expect_equal(
      b$meta$track,
      tr,
      sprintf("%s: meta records its own track", tr)
    )
    expect_equal(
      b$meta$column,
      anchor_spec(tr)$column,
      sprintf("%s: meta$column matches .ANCHOR_SPEC", tr)
    )
    expect_true(
      all(is.finite(b$cell$mean_rate) & b$cell$mean_rate > 0),
      sprintf("%s: every cell rate is finite and strictly positive", tr)
    )
  }

  # The four tracks are the same corpus expressed four ways, and each panel is
  # written so that the track's own view() reconstructs the same (num, den). The
  # four views differ only in which columns they read and how they recombine them
  # (occ derives exposure from t_start/t_end, mon multiplies pop x exposure, siz
  # divides by pop alone), so agreeing cell-for-cell is the check that the shared
  # builder is genuinely track-agnostic: everything track-specific lives in the
  # view, nothing in the estimator or the ladder.
  o <- brt$occ$cell %>% dplyr::arrange(country, admin_level)
  for (tr in setdiff(TRACKS, "occ")) {
    x <- brt[[tr]]$cell %>% dplyr::arrange(country, admin_level)
    expect_equal(
      nrow(x),
      nrow(o),
      sprintf("%s and occ cover the same cells", tr)
    )
    expect_equal(
      x$mr_tier,
      o$mr_tier,
      sprintf("%s and occ assign identical tiers", tr)
    )
    expect_true(
      max(abs(x$mean_rate / o$mean_rate - 1)) < 1e-12,
      sprintf(
        "%s and occ agree cell-for-cell on an identical (num, den) corpus",
        tr
      )
    )
  }
}

cat("\n## 2) Tier assignment: T1 / T1-zero / T2 on the cell table\n")
{
  for (tr in TRACKS) {
    cl <- brt[[tr]]$cell
    g <- function(c_, l_) dplyr::filter(cl, country == c_, admin_level == l_)
    expect_equal(
      g("KEN", "ADM1")$mr_tier,
      "T1",
      sprintf("%s: KEN/ADM1 (own cell, events observed) -> T1", tr)
    )
    expect_equal(
      g("TZA", "ADM2")$mr_tier,
      "T2",
      sprintf("%s: TZA/ADM2 (country seen, level not) -> T2", tr)
    )
    expect_true(
      is.na(g("TZA", "ADM2")$n_cell),
      sprintf("%s: TZA/ADM2 has no own-cell rows (n_cell NA)", tr)
    )
    expect_equal(
      g("ZERO", "ADM1")$mr_tier,
      "T1-zero",
      sprintf("%s: ZERO/ADM1 (cell observed, zero events) -> T1-zero", tr)
    )
  }

  # T1-zero is not cosmetic: `w = 0` must make the cell fall back to its country
  # parent exactly, with no continuity correction leaking in.
  z <- dplyr::filter(brt$occ$cell, country == "ZERO", admin_level == "ADM1")
  expect_equal(
    z$w_shrink,
    0,
    "T1-zero: EB weight is exactly 0 (Y = 0 => w = Y/(Y + 1/tau2) = 0)"
  )
  expect_equal(
    z$mean_rate,
    z$mr_parent,
    "T1-zero: the rate IS the country parent, exactly (no +0.5 fudge)"
  )
  expect_true(
    !identical(z$mr_tier, "T2"),
    "T1-zero is a DISTINCT label from T2 -- 'observed, saw nothing' is not 'never observed'"
  )
}

cat("\n## 3) Pooled EB shrinkage: w = Y / (Y + 1/tau2), exactly\n")
{
  b <- brt$occ
  k <- b$meta$k
  tau2 <- b$meta$tau2
  cl <- b$cell

  expect_true(
    is.finite(tau2) && tau2 > 0,
    "tau2_between is finite and positive"
  )
  expect_equal(
    k,
    1 / tau2,
    "pooled k is exactly 1/tau2 (the delta-method reparameterisation)"
  )

  # n_cell counts EVENTS under `pooled`, not rows -- the substitution that lets
  # the n/(n+k) shrinkage code be shared with the geometric branch verbatim.
  ken <- dplyr::filter(cl, country == "KEN", admin_level == "ADM1")
  bign <- dplyr::filter(cl, country == "BIGN", admin_level == "ADM1")
  tza <- dplyr::filter(cl, country == "TZA", admin_level == "ADM1")

  y_ken <- sum(corpus$num[
    corpus$country == "KEN" & corpus$admin_level == "ADM1"
  ])
  expect_equal(
    ken$n_cell,
    y_ken,
    "n_cell is the EVENT count Y_c under pooled, not the row count"
  )
  expect_true(
    ken$n_cell != sum(corpus$country == "KEN" & corpus$admin_level == "ADM1"),
    "...and that is genuinely different from the row count in this fixture"
  )

  for (nm in c("KEN", "BIGN", "TZA")) {
    r <- dplyr::filter(cl, country == nm, admin_level == "ADM1")
    expect_equal(
      r$w_shrink,
      r$n_cell / (r$n_cell + 1 / tau2),
      sprintf("%s/ADM1 w_shrink == Y/(Y + 1/tau2) exactly", nm)
    )
  }
  expect_true(
    bign$w_shrink > tza$w_shrink,
    "more events shrink less (w closer to 1)"
  )
  expect_true(
    bign$w_shrink > 0.98,
    "w -> 1 as Y -> infinity (the large cell is essentially untouched)"
  )

  # The shrunk estimate must sit BETWEEN the raw cell rate and its parent -- a
  # convex combination on the log scale, which is what w is for.
  lo <- pmin(ken$mr_raw, ken$mr_parent)
  hi <- pmax(ken$mr_raw, ken$mr_parent)
  expect_true(
    ken$mean_rate >= lo - 1e-12 && ken$mean_rate <= hi + 1e-12,
    "the shrunk rate lies between the raw cell rate and its parent"
  )
}

cat("\n## 4) resolve_base_rate: the full ladder T1 / T1-zero / T2 / T3 / T4\n")
{
  for (tr in TRACKS) {
    q <- resolve_base_rate(
      country = c("KEN", "ZERO", "TZA", "RWA", "RWA", "ZZZ"),
      admin_level = c("ADM1", "ADM1", "ADM2", "ADM1", "ADM3", "ADM1"),
      brt[[tr]]
    )
    expect_equal(
      q$br_tier,
      c("T1", "T1-zero", "T2", "T3-own", "T3-parent", "T4"),
      sprintf("%s: the six resolvable tiers, in one query", tr)
    )
    expect_true(
      all(is.finite(q$base_rate) & q$base_rate > 0),
      sprintf("%s: every resolved anchor is finite and strictly positive", tr)
    )
  }

  # Shape contract: resolve_base_rate() renames resolve_mean_rate()'s columns but
  # must be a pure delegate otherwise -- one lookup ladder, not two.
  q <- resolve_base_rate("KEN", "ADM1", brt$occ)
  qm <- resolve_mean_rate("KEN", "ADM1", brt$occ)
  expect_true(
    identical(names(q), c("base_rate", "br_tier")),
    "resolve_base_rate() returns tibble(base_rate, br_tier)"
  )
  expect_equal(
    q$base_rate,
    qm$mean_rate,
    "resolve_base_rate() is a pure delegate to resolve_mean_rate() (same value)"
  )
  expect_equal(
    q$br_tier,
    qm$mr_tier,
    "resolve_base_rate() is a pure delegate to resolve_mean_rate() (same tier)"
  )
}

cat(
  "\n## 5) T5: an unresolvable admin_level is a hard error, never a silent default\n"
)
{
  for (tr in TRACKS) {
    err <- tryCatch(
      {
        resolve_base_rate("RWA", "ADM4", brt[[tr]])
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    expect_true(
      !is.null(err) && grepl("T5", err),
      sprintf("%s: unrecognised admin_level -> T5 hard error", tr)
    )
  }
  err2 <- tryCatch(
    {
      resolve_base_rate("ZZZ", "ADM4", brt$occ)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err2) && grepl("T5", err2),
    "country AND admin_level both unresolvable -> T5 hard error"
  )
}

# mu = <that track's OTHER offsets> x anchor x exp(lp) is invariant under a
# global rescale of the table plus a compensating -log(c) intercept shift (the
# "global level is unidentifiable, only contrasts matter" property): the
# single strongest check.
cat(
  "\n## 6) Global rescaling leaves mu invariant, per track (intercept absorption)\n"
)
{
  # THE load-bearing check, generalised from test_mean_rate_table.R to
  # all four tracks.
  #
  #   mu = <that track's OTHER offsets> x anchor x exp(lp)
  #
  # Scaling every rate in the table by a constant c is exactly cancelled by
  # shifting lp by -log(c): the anchor's LEVEL is unidentifiable, only its
  # per-unit CONTRASTS carry information. This is what licenses re-centring the
  # Intercept prior at normal(0, 1) -- and it is also what would break, loudly,
  # if a track's anchor denominator ever stopped matching its offset product.
  c_scale <- 3.7
  rescale <- function(b) {
    b$cell$mean_rate <- b$cell$mean_rate * c_scale
    b$cell$mr_raw <- b$cell$mr_raw * c_scale
    b$cell$mr_parent <- b$cell$mr_parent * c_scale
    b$group$mean_rate <- b$group$mean_rate * c_scale
    b$global$mean_rate <- b$global$mean_rate * c_scale
    b
  }

  qc <- c("KEN", "ZERO", "TZA", "RWA", "ZZZ")
  ql <- c("ADM1", "ADM1", "ADM2", "ADM1", "ADM1")
  lp <- c(0.10, -0.20, 0.05, 0.30, -0.45)
  # The per-track "other offsets" product, i.e. everything in exp(sum(offsets))
  # EXCEPT the anchor: exposure days for occ/dur, pop x exposure for mon, pop
  # for siz. Values are arbitrary and deliberately unequal across rows.
  other_off <- list(
    occ = c(31, 28, 30, 31, 30),
    dur = c(31, 28, 30, 31, 30),
    mon = c(31, 28, 30, 31, 30) * c(5e4, 3e4, 2e4, 6e4, 1e5),
    siz = c(5e4, 3e4, 2e4, 6e4, 1e5)
  )

  for (tr in TRACKS) {
    before <- resolve_base_rate(qc, ql, brt[[tr]])
    after <- resolve_base_rate(qc, ql, rescale(brt[[tr]]))

    expect_true(
      all(abs(after$base_rate / before$base_rate - c_scale) < 1e-9),
      sprintf("%s: every resolved anchor scaled by exactly c", tr)
    )
    expect_equal(
      after$br_tier,
      before$br_tier,
      sprintf("%s: tier assignment is unaffected by a global rescale", tr)
    )

    mu_before <- other_off[[tr]] * before$base_rate * exp(lp)
    mu_after <- other_off[[tr]] * after$base_rate * exp(lp - log(c_scale))
    expect_true(
      all(abs(mu_before - mu_after) < 1e-8 * abs(mu_before)),
      sprintf(
        paste(
          "%s: mu = <other offsets> x anchor x exp(lp) is invariant",
          "under a global rescale + compensating intercept shift"
        ),
        tr
      )
    )
  }
}

cat(
  "\n## 7) The anchor cancels its offsets: exp(Intercept) ~ 1 at the prior mean\n"
)
{
  # The blowup this whole mechanism exists to fix was an intercept prior on the
  # wrong SCALE: normal(0, 1.5) on a linear predictor that had to span the entire
  # exposure/population range. With the anchor in, the implied rate at lp = 0 is
  # the empirical rate itself, so the prior-predictive total lands within an
  # order of magnitude of the observed total instead of ~1,000x above it.
  #
  # Checked here on the corpus the table was built from, per track.
  for (tr in TRACKS) {
    b <- brt[[tr]]
    res <- resolve_base_rate(corpus$country, corpus$admin_level, b)
    mu <- corpus$den * res$base_rate # lp = 0, i.e. Intercept = 0
    ratio <- sum(mu) / sum(corpus$num)
    expect_true(
      is.finite(ratio) && ratio > 0.1 && ratio < 10,
      sprintf(
        "%s: sum(mu) at Intercept = 0 is within 10x of the observed total (ratio %.3f)",
        tr,
        ratio
      )
    )
  }
}

cat("\n## 8) attach_base_rate: joins the column the FORMULA will name\n")
{
  d <- tibble::tibble(
    country = c("KEN", "RWA", "ZZZ"),
    admin_level = c("ADM1", "ADM1", "ADM1")
  )
  for (tr in TRACKS) {
    spec <- anchor_spec(tr)
    out <- attach_base_rate(d, track_opt(tr), tr, brt = brt[[tr]])
    expect_true(
      spec$column %in% names(out),
      sprintf("%s: adds the `%s` column", tr, spec$column)
    )
    expect_true(
      spec$tier_col %in% names(out),
      sprintf("%s: adds the `%s` tier column", tr, spec$tier_col)
    )
    expect_true(
      !is.null(attr(out, spec$attr_name)),
      sprintf("%s: stashes the table as attr(d, '%s')", tr, spec$attr_name)
    )
    expect_true(
      all(is.finite(out[[spec$column]]) & out[[spec$column]] > 0),
      sprintf("%s: joined anchor is finite and strictly positive", tr)
    )
  }
}

# The formula gate (.OFFSET_SPEC, utils/config.R) and the draws-bundle flag
# (.BUNDLE_ANCHOR_COL, utils/simulator_core.R) name the same column for every
# track. Neither module may source the other (config.R is sourced by
# everything; anchor.R depends on paths/provenance/occ_data), so that agreement
# is a convention, and a convention with no test is a bug waiting for a rename.
cat("\n## 9) Formula gate, bundle flag and table column name the SAME thing\n")
{
  # utils/config.R (.OFFSET_SPEC) must not source utils/anchor.R (.ANCHOR_SPEC),
  # and utils/simulator_core.R (.BUNDLE_ANCHOR_COL) must not source either --
  # the dependency closure would be circular / far too wide. So the three
  # registries agree by convention, and this is the test that keeps them honest.
  for (tr in TRACKS) {
    osp <- .offset_spec(tr)
    asp <- anchor_spec(tr)
    expect_equal(
      osp$anchor_col,
      asp$column,
      sprintf("%s: .OFFSET_SPEC$anchor_col == .ANCHOR_SPEC$column", tr)
    )
    expect_equal(
      osp$block,
      asp$opt_block,
      sprintf("%s: .OFFSET_SPEC$block == .ANCHOR_SPEC$opt_block", tr)
    )
  }

  # .BUNDLE_ANCHOR_COL covers occ/dur/mon ONLY, deliberately: size's anchor is
  # `mean_rate` and already travels as the older `offset_meanrate` flag, so a
  # second flag would multiply log(mean_rate) in twice.
  expect_equal(
    sort(names(.BUNDLE_ANCHOR_COL)),
    sort(c("occurrence", "duration", "monthly")),
    ".BUNDLE_ANCHOR_COL covers exactly the three non-size tracks"
  )
  for (tr in c("occ", "dur", "mon")) {
    asp <- anchor_spec(tr)
    expect_equal(
      .BUNDLE_ANCHOR_COL[[asp$opt_block]],
      asp$column,
      sprintf("%s: .BUNDLE_ANCHOR_COL == .ANCHOR_SPEC$column", tr)
    )
  }
  expect_true(
    is.null(.BUNDLE_ANCHOR_COL[["size"]]),
    "size is absent from .BUNDLE_ANCHOR_COL (its anchor rides offset_meanrate)"
  )

  # And the emitted formula offsets: the anchor term must appear iff the gate is
  # on, and must name the same column.
  for (tr in TRACKS) {
    asp <- anchor_spec(tr)
    want <- sprintf("offset(log(%s))", asp$column)
    on <- .offset_terms(track_opt(tr), tr)
    expect_true(
      want %in% on,
      sprintf("%s: gate on -> formula emits %s", tr, want)
    )
  }
  off_opt <- list(
    occurrence = list(),
    duration = list(),
    monthly = list(),
    size = list(offset_meanrate = FALSE)
  )
  for (tr in TRACKS) {
    asp <- anchor_spec(tr)
    want <- sprintf("offset(log(%s))", asp$column)
    expect_true(
      !(want %in% .offset_terms(off_opt, tr)),
      sprintf(
        "%s: gate off -> no anchor offset (every pre-anchor config stays inert)",
        tr
      )
    )
  }
}

cat("\n## 10) .anchor_log_offset(): the kernel side of the same flag\n")
{
  frame <- tibble::tibble(
    base_rate_occ = c(1e-3, 4e-3),
    base_rate_dur = c(2e-3, 2e-3)
  )

  expect_equal(
    .anchor_log_offset(list(offset_baserate = FALSE), frame),
    0,
    "flag off -> the additive identity 0 (every mu block can add it unconditionally)"
  )
  expect_equal(
    .anchor_log_offset(list(), frame),
    0,
    "flag absent (a pre-anchor bundle) -> 0, not an error"
  )
  expect_equal(
    .anchor_log_offset(
      list(offset_baserate = TRUE, anchor_col = "base_rate_occ"),
      frame
    ),
    log(frame$base_rate_occ),
    "flag on -> log(anchor) for every row"
  )
  expect_equal(
    .anchor_log_offset(
      list(offset_baserate = TRUE, anchor_col = "base_rate_occ"),
      frame,
      idx = 2L
    ),
    log(frame$base_rate_occ[2]),
    "idx subsets the frame, matching the mu blocks' per-location indexing"
  )

  # The three guards. A missing column would otherwise silently rescale every
  # simulated prediction -- the exact failure mode the flag protocol exists to
  # prevent.
  e1 <- tryCatch(
    {
      .anchor_log_offset(list(offset_baserate = TRUE), frame)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(e1) && grepl("anchor_col", e1),
    "offset_baserate = TRUE with no anchor_col -> hard error naming anchor_col"
  )
  e2 <- tryCatch(
    {
      .anchor_log_offset(
        list(offset_baserate = TRUE, anchor_col = "base_rate_mon"),
        frame
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(e2) && grepl("base_rate_table_paths", e2),
    "anchored fit + un-anchored frame -> hard error naming build_location_frame()"
  )
  e3 <- tryCatch(
    {
      .anchor_log_offset(
        list(offset_baserate = TRUE, anchor_col = "base_rate_occ"),
        tibble::tibble(base_rate_occ = c(1e-3, 0))
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(e3) && grepl("positive", e3),
    "a non-positive anchor (rate table did not resolve) -> hard error"
  )
}

cat("\n## 11) resolve_base_rate_paths: opt-in, per track, no silent NULL\n")
{
  spec_occ <- anchor_spec("occ")
  sidecar <- tempfile(fileext = ".rds")
  saveRDS(brt$occ, sidecar)

  expect_true(
    is.null(resolve_base_rate_paths(list())),
    "no tracks supplied -> NULL (build_location_frame() takes its inert path)"
  )
  expect_true(
    is.null(resolve_base_rate_paths(list(occ = list(occurrence = list())))),
    "track supplied but anchor OFF -> NULL (every deployed fit today)"
  )

  got <- resolve_base_rate_paths(
    list(occ = list(occurrence = list(offset_baserate = TRUE))),
    list(base_rate_table_paths = list(occ = sidecar))
  )
  expect_equal(
    got,
    list(occ = sidecar),
    "explicit config path wins, exactly as mean_rate_table_path does"
  )

  got2 <- resolve_base_rate_paths(
    list(
      occ = list(occurrence = list(offset_baserate = TRUE)),
      dur = NULL,
      mon = list(monthly = list(offset_baserate = FALSE))
    ),
    list(base_rate_table_paths = list(occ = sidecar))
  )
  expect_equal(
    names(got2),
    "occ",
    "a NULL track and an un-anchored track are both dropped"
  )

  err <- tryCatch(
    {
      resolve_base_rate_paths(
        list(occ = list(occurrence = list(offset_baserate = TRUE))),
        list(base_rate_table_paths = list(occ = "/nonexistent/brt.rds"))
      )
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  expect_true(
    !is.null(err) && grepl("offset_baserate", err),
    "anchored fit + missing sidecar -> hard error, never a silently un-anchored frame"
  )
}

cat(
  "\n## 12) Size: `geometric` still reproduces the original build_mean_rate_table()\n"
)
{
  # The size track's default estimator is pooled (four-track parity).
  # `geometric` is retained, not deprecated, so a
  # like-for-like refit of the deployed archive needs only
  # `size.mean_rate.estimator: geometric` -- which is only true if the shared
  # builder still agrees with provenance.R's original cell-for-cell.
  opt_siz <- track_opt("siz")
  opt_siz$size$analysis_data_path <- panels$siz
  opt_siz$size$mean_rate <- list(estimator = "geometric")

  shared <- suppressWarnings(build_base_rate_table(
    opt_siz,
    "siz",
    path = panels$siz,
    enso_grouping_path = ENSO_CSV
  ))
  orig <- suppressWarnings(build_mean_rate_table(
    opt_siz,
    path = panels$siz,
    enso_grouping_path = ENSO_CSV
  ))

  a <- shared$cell %>% dplyr::arrange(country, admin_level)
  b <- orig$cell %>% dplyr::arrange(country, admin_level)
  expect_equal(nrow(a), nrow(b), "same number of cells")
  expect_equal(a$mr_tier, b$mr_tier, "identical tier assignment")
  expect_true(
    max(abs(a$mean_rate / b$mean_rate - 1)) < 1e-12,
    "identical cell rates (relative difference < 1e-12)"
  )
  expect_equal(shared$meta$k, orig$meta$k, "identical EB shrinkage constant k")

  # ...and pooled genuinely MOVES the numbers (otherwise section 12 would be
  # vacuous). Joined on the cell key rather than compared positionally:
  # `geometric` drops zero-event rows before the estimator sees them, so the all-zero ZERO country has no
  # cells at all in `b` while `pooled` keeps all three of them.
  cmp <- dplyr::inner_join(
    dplyr::select(brt$siz$cell, country, admin_level, mr_pooled = mean_rate),
    dplyr::select(b, country, admin_level, mr_geom = mean_rate),
    by = c("country", "admin_level")
  )
  expect_true(
    nrow(cmp) == nrow(b) && nrow(cmp) > 0,
    "every geometric cell has a pooled counterpart"
  )
  expect_true(
    max(abs(cmp$mr_pooled / cmp$mr_geom - 1)) > 1e-6,
    "pooled differs from geometric (the default switch is a real change, measured in s31)"
  )
  expect_true(
    !any(brt$siz$cell$country == "ZERO" & is.na(brt$siz$cell$mean_rate)) &&
      all(brt$siz$cell$mean_rate[brt$siz$cell$country == "ZERO"] > 0),
    "pooled keeps the all-zero country and gives it a strictly positive parent rate"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0) {
  quit(status = 1)
}
