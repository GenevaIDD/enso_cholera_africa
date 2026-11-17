# test_results_sensitive_pop.R — HARD: the phase-09 per-draw population
#   helpers and the annualised-incidence helpers (utils/results_derive.R).
#
# Usage:
#   Rscript analysis/tests/test_results_sensitive_pop.R
#
# Details:
#   Quick-safe: no fit, no data, no MCMC. Covers the helpers behind
#   report_effects (dir_pop_ci), report_enso_exposure (exposure_mask(); its
#   `either` direction reproduces the former report_sensitive_pop rows) and
#   report_forecast (incidence), plus enso_conds_named() in utils/figures_cb.R.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/results_core.R") # base, figures_cb (enso_conds_named), registries
source("analysis/utils/gid.R") # flatten_gid
source("analysis/utils/population.R") # resolve_gid_population
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
expect_equal <- function(a, b, label, tol = 1e-9) {
  expect(isTRUE(all.equal(a, b, tolerance = tol)), label)
}

cat("# phase-09 per-draw population + incidence helpers\n")

# enso_conds_named(c(-2,-1,1,2)) reproduces the former COMPOSITE_CONDS literal.
# ---- condition labels -------------------------------------------------------
old_lit <- c(
  "La Niña (−2.0)" = -2,
  "La Niña (−1.0)" = -1,
  "El Niño (+1.0)" = 1,
  "El Niño (+2.0)" = 2
)
expect(
  identical(enso_conds_named(c(-2, -1, 1, 2)), old_lit),
  "enso_conds_named(c(-2,-1,1,2)) equals the former COMPOSITE_CONDS literal"
)
expect(identical(COMPOSITE_CONDS, old_lit), "COMPOSITE_CONDS unchanged")
expect(
  identical(
    names(enso_conds_named(list(0.5, 0))),
    c("El Niño (+0.5)", "Neutral (0.0)")
  ),
  "enso_conds_named() accepts a YAML list, labels 0 as neutral, keeps order"
)
expect(
  inherits(try(enso_conds_named(c(1, 1)), silent = TRUE), "try-error"),
  "duplicate conditions are refused"
)

# sensitive_mask() is symmetric and strict at log(1 + thr).
# ---- sensitive_mask ---------------------------------------------------------
thr <- 0.10
b <- log1p(thr)
x <- c(b, -b, b + 1e-9, -b - 1e-9, 0)
expect(
  identical(sensitive_mask(x, thr), c(FALSE, FALSE, TRUE, TRUE, FALSE)),
  "sensitive_mask: strict at +/- log(1+thr), symmetric"
)

# pop_mass_by_draw() equals the former per-ADM1 dir_pop_ci() loop exactly.
# ---- synthetic draw array: 3 ADM1 x 2 conditions x 4 draws ------------------
conds <- enso_conds_named(c(-1, 1))
set.seed(1)
arr <- array(
  rnorm(3 * 2 * 4, sd = 0.2),
  dim = c(3, 2, 4),
  dimnames = list(c("KEN_1", "KEN_2", "UGA_1"), names(conds), NULL)
)
pops <- c(KEN_1 = 100, KEN_2 = NA, UGA_1 = 50)

# former dir_pop_ci() arithmetic, reproduced literally
old_inc <- matrix(0, 2, 4)
for (i in 1:3) {
  p_i <- unname(pops[i])
  if (is.na(p_i)) {
    p_i <- 0
  }
  old_inc <- old_inc + (arr[i, , ] > 0) * p_i
}
expect(
  identical(pop_mass_by_draw(arr > 0, unname(pops)), unname(old_inc)),
  "pop_mass_by_draw() == former per-ADM1 loop (identical doubles, dimnames aside), NA pop -> 0"
)

# align_adm1_array() re-indexes onto a universe, absent ADM1 -> FALSE.
# ---- align_adm1_array -------------------------------------------------------
al <- align_adm1_array(
  arr > 0,
  dimnames(arr)[[1]],
  c("UGA_1", "TZA_9", "KEN_1")
)
expect(
  identical(dim(al), c(3L, 2L, 4L)) &&
    !any(al["TZA_9", , ]) &&
    identical(al["UGA_1", , ], arr["UGA_1", , ] > 0) &&
    identical(al["KEN_1", , ], arr["KEN_1", , ] > 0),
  "align_adm1_array(): universe order, absent ADM1 all FALSE"
)
m1 <- align_adm1_array(
  sensitive_mask(arr, thr),
  dimnames(arr)[[1]],
  c("KEN_1", "UGA_1")
)
m2 <- align_adm1_array(
  sensitive_mask(-arr[1, , , drop = FALSE] * 0, thr),
  "KEN_1",
  c("KEN_1", "UGA_1")
)
expect(
  identical(Reduce(`|`, list(m1, m2)), m1),
  "union with an all-FALSE track is a no-op"
)

# dir_pop_ci() on a cached draw array gives region/country sums by hand;
# bind_adm1_draws() stacks groups.
# ---- bind_adm1_draws + dir_pop_ci ------------------------------------------
part <- function(rows, g) {
  list(
    logrr = arr[rows, , , drop = FALSE],
    adm1 = tibble::tibble(
      adm1_label = rows,
      adm1_gid = rows,
      country = substr(rows, 1, 3),
      enso_group = g
    )
  )
}
draws <- bind_adm1_draws(list(
  part(c("KEN_1", "KEN_2"), "enso_east"),
  NULL,
  part("UGA_1", "enso_horn")
))
expect(
  identical(unname(draws$logrr), unname(arr)) &&
    identical(dimnames(draws$logrr)[[1]], rownames(arr)),
  "bind_adm1_draws(): groups stacked along the ADM1 dimension, NULL parts dropped"
)
dp <- dir_pop_ci(draws, "enso_east", pops, conds)
expect(
  identical(
    dp$region_inc,
    pop_mass_by_draw(arr[1:2, , , drop = FALSE] > 0, unname(pops[1:2]))
  ),
  "dir_pop_ci(): region_inc is the group's population mass per draw"
)
expect(
  nrow(dp$region_ci) == 4L &&
    all(dp$region_ci$region == "East") &&
    identical(sort(unique(dp$country_ci$country)), "KEN"),
  "dir_pop_ci(): one row per condition x direction, country rows"
)
expect(
  is.null(dir_pop_ci(draws, "enso_west", pops, conds)),
  "dir_pop_ci(): absent group -> NULL"
)

# .pop_summary_rows() / .q_mean95().
# ---- summaries --------------------------------------------------------------
mat <- matrix(c(1, 2, 3, 4, 10, 20, 30, 40), 2, byrow = TRUE)
ps <- .pop_summary_rows(mat, conds)
expect_equal(ps$pop_mean, c(2.5, 25), ".pop_summary_rows(): mean across draws")
expect_equal(
  ps$pop_hi,
  unname(apply(mat, 1, quantile, 0.975)),
  ".pop_summary_rows(): q975"
)
expect(
  identical(
    names(.pop_ci_rows(mat, "increase", conds)),
    c("enso_label", "enso_val", "direction", "pop_lo", "pop_med", "pop_hi")
  ),
  ".pop_ci_rows() keeps the report_effects column contract"
)
qm <- .q_mean95(c(1, 2, 3, NA))
expect_equal(
  c(qm$mean, qm$med),
  c(2, 2),
  ".q_mean95(): mean + median, NA dropped"
)

sh <- add_pop_shares(tibble::tibble(
  pop = c(50, 0),
  pop_lo = c(10, 0),
  pop_hi = c(90, 0),
  pop_universe = c(100, 0)
))
expect(
  identical(sh$pct, c(0.5, NA)) &&
    identical(sh$pct_lo, c(0.1, NA)) &&
    identical(sh$pct_hi, c(0.9, NA)) &&
    !"pct_mean" %in% names(sh),
  "add_pop_shares(): pct* = pop* / pop_universe for the columns present, NA on an empty scope"
)

# ---- shared per-draw machinery: universe scopes, mask reductions, grid -------
expect(
  identical(names(TRACK_ADVERSE_SIGN), names(TRACK_CODE)) &&
    TRACK_ADVERSE_SIGN[["duration"]] == -1 &&
    all(TRACK_ADVERSE_SIGN[c("occurrence", "monthly")] == 1),
  "TRACK_ADVERSE_SIGN covers every track; duration's adverse direction is negative"
)
# A 3-ADM1 universe over two tracks: UGA_1 is missing from track b.
uni <- tibble::tibble(
  adm1_gid = c("KEN_1", "KEN_2", "UGA_1"),
  country = c("KEN", "KEN", "UGA"),
  enso_group = c("enso_east", "enso_east", "enso_east"),
  region = "East",
  pop = c(100, NA, 50),
  pop_missing = c(FALSE, TRUE, FALSE)
)
la <- align_adm1_array(arr, rownames(arr), uni$adm1_gid, fill = NA_real_)
lb <- align_adm1_array(
  arr[1:2, , , drop = FALSE] * -1,
  rownames(arr)[1:2],
  uni$adm1_gid,
  fill = NA_real_
)
dr <- list(
  tracks = c("a", "b"),
  conds = conds,
  n_draws = 4L,
  universe = uni,
  U = uni$adm1_gid,
  pop0 = c(100, 0, 50),
  logrr = list(a = la, b = lb)
)
expect(
  all(is.na(lb["UGA_1", , ])),
  "align_adm1_array(fill = NA) marks an absent ADM1 NA"
)
sc <- enso_scopes(uni)
expect(
  identical(sc$level, c("continent", "region", "country", "country")) &&
    identical(sc$idx[[4]], 3L) &&
    identical(sc$idx[[1]], 1:3),
  "enso_scopes(): Africa, one region, two countries with their ADM1 indices"
)
mk <- enso_track_masks(dr, function(x, t) x > 0)
expect(
  identical(dim(mk$a), c(3L, 2L, 4L)) &&
    !any(mk$b["UGA_1", , ]) &&
    identical(mk$a["KEN_1", , ], arr["KEN_1", , ] > 0),
  "enso_track_masks(): applies the test per track, absent ADM1 -> FALSE"
)
cm <- enso_combine_masks(mk)
expect(
  identical(cm$any, mk$a | mk$b) &&
    identical(cm$concordant2, mk$a & mk$b) &&
    identical(names(cm), c("a", "b", "any", "concordant2")),
  "enso_combine_masks(): any = union, concordant2 = both of two tracks"
)
rows <- mask_pop_rows(cm, dr)
afr_a <- dplyr::filter(rows, track == "a", level == "continent")
expect(
  nrow(rows) == 4L * 4L * 2L &&
    identical(afr_a$pop_mean, rowMeans(pop_mass_by_draw(mk$a, dr$pop0))) &&
    all(afr_a$pop_universe == 150) &&
    all(afr_a$n_pop_missing == 1L) &&
    isTRUE(all.equal(afr_a$pct_mean, afr_a$pop_mean / 150)),
  "mask_pop_rows(): one row per mask x scope x condition; population per draw, shares of the scope"
)
uga <- dplyr::filter(rows, track == "a", scope == "UGA")
expect(
  isTRUE(all.equal(
    unname(uga$pop_mean),
    unname(50 * rowMeans(mk$a["UGA_1", , , drop = FALSE][1, , ]))
  )),
  "mask_pop_rows(): a one-ADM1 country scope is that ADM1's population times its hit rate"
)
g <- prob_grid_rows(cm, dr, cutoffs = c(0.5, 1))
pa <- apply(mk$a, c(1, 2), mean)
g1 <- dplyr::filter(g, track == "a", level == "continent", cutoff == 1)
expect(
  nrow(g) == 4L * 4L * 2L * 2L &&
    identical(g1$n_adm1, unname(colSums(pa >= 1))) &&
    identical(g1$pop, unname(colSums((pa >= 1) * dr$pop0))) &&
    isTRUE(all.equal(g1$pct, g1$pop / 150)),
  "prob_grid_rows(): population of ADM1s whose hit probability >= cutoff, with shares"
)

# ---- report_enso_exposure directions (exposure_mask) ------------------------
# `either` at m = thr is the former report_sensitive_pop rule: same masks, same
# per-draw rows (that table also carried a point estimate, not compared here).
old_m <- enso_track_masks(dr, function(x, t) sensitive_mask(x, thr))
old_m$any <- Reduce(`|`, old_m)
new_m <- enso_combine_masks(enso_track_masks(dr, function(x, t) {
  exposure_mask(x, 1, "either", thr)
}))
pc <- c(
  "track",
  "level",
  "scope",
  "enso_label",
  "pop_universe",
  "pop_mean",
  "pop_lo",
  "pop_med",
  "pop_hi",
  "n_adm1_mean",
  "pct_mean",
  "pct_lo",
  "pct_med",
  "pct_hi"
)
expect(
  identical(
    mask_pop_rows(old_m, dr)[, pc],
    dplyr::filter(mask_pop_rows(new_m, dr), track != "concordant2")[, pc]
  ),
  "exposure_mask(either, m = sensitive_irr) reproduces the report_sensitive_pop rows"
)
expect(
  identical(
    exposure_mask(arr, 1, "either", thr),
    exposure_mask(arr, -1, "either", thr)
  ),
  "exposure_mask(either) ignores the adverse sign"
)
# adverse with the duration flip (sign -1) = a negative log IRR, i.e. dir_pop_ci()'s decrease mass
expect(
  identical(
    exposure_mask(arr, TRACK_ADVERSE_SIGN[["duration"]], "adverse", 0),
    arr < 0
  ) &&
    identical(
      exposure_mask(arr, TRACK_ADVERSE_SIGN[["occurrence"]], "adverse", 0),
      arr > 0
    ) &&
    identical(exposure_mask(arr, 1, "protective", 0), arr < 0),
  "exposure_mask(adverse/protective, m = 0) follows TRACK_ADVERSE_SIGN"
)
dpc <- dir_pop_ci(
  list(
    logrr = arr,
    adm1 = tibble::tibble(
      adm1_gid = rownames(arr),
      country = c("KEN", "KEN", "UGA"),
      enso_group = "enso_east"
    )
  ),
  "enso_east",
  pop_map = pops,
  enso_conds = conds,
  group_short = c(enso_east = "East")
)
expect(
  identical(
    pop_mass_by_draw(exposure_mask(arr, -1, "adverse", 0), unname(pops)),
    dpc$region_dec
  ),
  "adverse duration mass at m = 0 == dir_pop_ci() decrease mass"
)
expect(
  inherits(
    try(exposure_mask(arr, 1, "sideways", 0), silent = TRUE),
    "try-error"
  ),
  "exposure_mask(): an unknown direction stops"
)
# concordant2 under `either`: two tracks beyond m in opposite directions count
opp <- list(
  tracks = c("a", "b"),
  logrr = list(a = array(0.5, c(1, 1, 1)), b = array(-0.5, c(1, 1, 1)))
)
co_e <- enso_combine_masks(enso_track_masks(opp, function(x, t) {
  exposure_mask(x, 1, "either", 0.1)
}))
co_a <- enso_combine_masks(enso_track_masks(opp, function(x, t) {
  exposure_mask(x, 1, "adverse", 0.1)
}))
expect(
  isTRUE(co_e$concordant2[1, 1, 1]) && !isTRUE(co_a$concordant2[1, 1, 1]),
  "concordant2 under `either` counts opposite-direction tracks; under `adverse` it does not"
)
# any_credible: the track-wise union of the per-track certain sets
gc <- prob_grid_rows(cm, dr, cutoffs = c(0.5, 1), credible_union = c("a", "b"))
pb <- apply(mk$b, c(1, 2), mean)
gca <- dplyr::filter(
  gc,
  track == "any_credible",
  level == "continent",
  cutoff == 0.5
)
expect(
  identical(dplyr::filter(gc, track != "any_credible"), g) &&
    identical(gca$n_adm1, unname(colSums(pa >= 0.5 | pb >= 0.5))) &&
    identical(gca$pop, unname(colSums((pa >= 0.5 | pb >= 0.5) * dr$pop0))),
  "prob_grid_rows(credible_union): any_credible = some single track has P >= p; other rows unchanged"
)
# three independent tracks at P = 0.65 each: the union passes 0.95, no track does
w <- array(rep(c(rep(TRUE, 13), rep(FALSE, 7)), 1), c(1, 1, 20))
tri <- list(
  tracks = c("a", "b", "c"),
  conds = conds[1],
  universe = uni[1, ],
  U = uni$adm1_gid[1],
  pop0 = 100,
  logrr = NULL
)
tm <- list(
  a = w,
  b = w[,, c(8:20, 1:7), drop = FALSE],
  c = w[,, c(15:20, 1:14), drop = FALSE]
)
tm$any <- Reduce(`|`, tm)
gt <- prob_grid_rows(tm, tri, cutoffs = 0.95, credible_union = c("a", "b", "c"))
gtc <- dplyr::filter(gt, level == "continent")
expect(
  gtc$n_adm1[gtc$track == "any"] == 1L &&
    gtc$n_adm1[gtc$track == "any_credible"] == 0L,
  "union certainty inflates (three tracks at 0.65 -> any >= 0.95); any_credible does not"
)
expect(
  identical(ENSO_EXPOSURE_DIRECTIONS, c("adverse", "protective", "either")),
  "ENSO_EXPOSURE_DIRECTIONS: adverse, protective, either"
)

# window_n_months(), per_rep_totals() annualisation, build_sim_scope_pop().
# ---- incidence helpers ------------------------------------------------------
expect(
  window_n_months("2026-06-01", "2026-06-01") == 1L &&
    window_n_months("2026-06-01", "2027-10-01") == 17L &&
    window_n_months("2026-12-15", "2027-01-01") == 2L,
  "window_n_months(): inclusive, year wrap"
)
toy <- tibble::tibble(
  replicate_id = c(1L, 1L, 2L),
  month = as.Date(c("2026-06-01", "2026-07-01", "2026-06-01")),
  n_cases = c(10, 20, 6),
  n_outbreaks = c(1L, 0L, 1L),
  run = "forecast",
  analog_set = NA_character_,
  n_months = 2L,
  pop = 1e6
)
pr <- per_rep_totals(toy, character(0))
expect_equal(
  pr$inc_ann_100k,
  c(30, 6) / 2 * 12 / 1e6 * 1e5,
  "per_rep_totals(): annualised incidence over the calendar window, not the rows present"
)
expect(
  all(is.na(
    per_rep_totals(dplyr::mutate(toy, pop = 0), character(0))$inc_ann_100k
  )),
  "per_rep_totals(): zero population -> NA incidence"
)

gid_pop <- tibble::tibble(
  gid = c("KEN_1_1", "KEN_1_2_1", "UGA_1_1"),
  gadm_level = c("ADM1", "ADM2", "ADM1"),
  iso3 = c("KEN", "KEN", "UGA"),
  year = 2025L,
  Q = c(1000, 400, 500)
)
gids <- tibble::tibble(
  gid = c("KEN.1_1", "KEN.1.2_1", "UGA.1_1", "UGA.9_1"),
  country = c("KEN", "KEN", "UGA", "UGA"),
  enso_group = c("enso_east", "enso_east", "enso_east", "enso_east")
)
sp <- build_sim_scope_pop(gids, gid_pop, 2025L)
afr <- sp[sp$level == "continent", ]
expect(
  afr$pop == 1900 &&
    afr$n_gid == 4L &&
    afr$n_gid_pop_missing == 1L &&
    afr$n_gid_adm1_overlap == 1L,
  "build_sim_scope_pop(): dotted gids flattened, missing pop counted, ADM1/ADM2 overlap counted"
)
expect(
  identical(sp$pop[sp$level == "country" & sp$scope == "KEN"], 1400) &&
    identical(sp$region[sp$level == "region"], "East"),
  "build_sim_scope_pop(): country + region rows"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
