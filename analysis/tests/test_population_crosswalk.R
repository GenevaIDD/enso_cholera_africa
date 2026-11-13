# test_population_crosswalk.R — the LP x gid crosswalk algebra, end to end.
#
# HARD gate (see analysis/tests/run_all.R).
#
# WHAT THIS COVERS
# ----------------
# runPopulationTask() with `extract_fn` replaced by a deterministic
# area-proportional stub. No raster, no network, no data/ access, no
# OutbreakExtractR install — the whole geometry pipeline (LP selection, vector
# intersection, GEOMETRYCOLLECTION sanitising, sliver dropping, share algebra,
# cross-border binding) runs on a synthetic grid whose answers are known in
# closed form.
#
# THE STUB IS NOT A NO-OP. It reproduces the two behaviours of
# estimate_pop_for_geometries() that the kernel is written around:
#   * the raster is COUNTRY-CLIPPED — mass outside the task's country is 0, which
#     is what makes summing P_direct across tasks correct for a cross-border LP;
#   * `extracted[extracted <= 0] <- NA_real_` (R/add_population.R:508-510) — so
#     the NA-coalescing path for zero-mass pieces is genuinely exercised.
#
# THE FIXTURE (all coordinates in metres, EPSG:3857 so st_area is planar)
#
#     y=1000 +----+----+----+----+----+----+----+----+   +----+----+
#            | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 8  |   | 1  | 2  |
#     y=0    +----+----+----+----+----+----+----+----+   +----+----+
#            0   1k                                 8k  8k        10k
#                          country AAA                   country BBB
#
#   ADM1 = the 8 (resp. 2) cells; ADM2 = each cell split at y = 500.
#   Density: AAA 10000 people/km^2, BBB 20000 people/km^2. A cell is 1 km^2.
#
#   lp_single  x[0, 1000]      = AAA cell 1 exactly
#   lp_comp7   x[1000, 8000]   = AAA cells 2..8  (7-way composite)
#   lp_sliver  x[0, 1001]      = cell 1 + a 1 m strip of cell 2 (sub-threshold)
#   lp_cross   x[7500, 8500]   = half of AAA cell 8 + half of BBB cell 1
#
# Run: Rscript analysis/tests/test_population_crosswalk.R

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

sf::sf_use_s2(FALSE)

source("analysis/utils/population_core.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
near <- function(a, b, tol = 1e-8) {
  isTRUE(length(a) == length(b)) && all(is.finite(a)) && all(abs(a - b) <= tol)
}

CRS <- 3857
DENS <- c(AAA = 10000, BBB = 20000)   # people per km^2

box <- function(x0, x1, y0 = 0, y1 = 1000) {
  sf::st_polygon(list(cbind(c(x0, x1, x1, x0, x0), c(y0, y0, y1, y1, y0))))
}

# --- GADM fixture ------------------------------------------------------------
mk_gadm <- function(iso3, x0, n) {
  a1 <- lapply(seq_len(n), function(i) box(x0 + (i - 1) * 1000, x0 + i * 1000))
  g1 <- sf::st_sf(
    gid        = sprintf("%s.%d_1", iso3, seq_len(n)),
    gadm_level = "ADM1",
    gid_iso3   = iso3,
    geometry   = sf::st_sfc(a1, crs = CRS)
  )
  a2 <- list(); nm <- character(0)
  for (i in seq_len(n)) for (j in 1:2) {
    a2 <- c(a2, list(box(x0 + (i - 1) * 1000, x0 + i * 1000,
                         (j - 1) * 500, j * 500)))
    nm <- c(nm, sprintf("%s.%d.%d_1", iso3, i, j))
  }
  g2 <- sf::st_sf(gid = nm, gadm_level = "ADM2", gid_iso3 = iso3,
                  geometry = sf::st_sfc(a2, crs = CRS))
  rbind(g1, g2)
}

gadm <- rbind(mk_gadm("AAA", 0, 8), mk_gadm("BBB", 8000, 2))

lp_sf <- sf::st_sf(
  lctn_pr  = c("lp_single", "lp_comp7", "lp_sliver", "lp_cross"),
  geometry = sf::st_sfc(box(0, 1000), box(1000, 8000),
                        box(0, 1001), box(7500, 8500), crs = CRS)
)

# --- the stub ----------------------------------------------------------------
country_geom <- lapply(c("AAA", "BBB"), function(k)
  sf::st_union(sf::st_geometry(gadm[gadm$gid_iso3 == k & gadm$gadm_level == "ADM1", ])))
names(country_geom) <- c("AAA", "BBB")

N_CALLS <- new.env(parent = emptyenv()); N_CALLS$n <- 0L

stub_extract <- function(geom_sf, country_iso3, year, raster_dir,
                         boundary_cache_dir, ...) {
  N_CALLS$n <- N_CALLS$n + 1L
  cg <- country_geom[[country_iso3]]
  g  <- sf::st_geometry(geom_sf)
  out <- vapply(seq_along(g), function(i) {
    p <- suppressWarnings(sf::st_intersection(g[i], cg))
    if (!length(p)) return(0)
    sum(suppressWarnings(as.numeric(sf::st_area(p))), na.rm = TRUE) / 1e6 *
      DENS[[country_iso3]]
  }, numeric(1))
  out[out <= 0] <- NA_real_        # mirrors R/add_population.R:508-510
  out
}

run_task <- function(iso3) {
  runPopulationTask(
    iso3 = iso3, year = 2020L, gadm_sf = gadm, lp_sf = lp_sf,
    extract_fn = stub_extract,
    sliver_mass_frac = 1e-2,       # forces the FRAC branch, not sliver_mass_min
    sliver_mass_min  = 1.0,
    wpp_total        = 8 * 10000,
    require_boundary = FALSE,
    quiet            = TRUE
  )
}

cat("# LP x gid population crosswalk\n")

A <- run_task("AAA")
B <- run_task("BBB")

# --- 1. one raster load per task --------------------------------------------
cat("\n[1] a task extracts every geometry role in a single call\n")
expect(N_CALLS$n == 2L,
       "exactly one extract_fn call per (country, year) task")
expect(setequal(unique(A$crosswalk$gadm_level), c("ADM1", "ADM2")),
       "both GADM levels came out of that one call")

# --- 2. LP membership is by positive-area overlap ----------------------------
cat("\n[2] LP selection: positive-area overlap, not edge touching\n")
expect(setequal(A$lp_pop$lctn_pr,
                c("lp_single", "lp_comp7", "lp_sliver", "lp_cross")),
       "AAA's task sees all four LPs")
expect(setequal(B$lp_pop$lctn_pr, "lp_cross"),
       "BBB's task sees only lp_cross — lp_comp7 merely abuts the border")
expect(isTRUE(A$lp_pop$cross_border[A$lp_pop$lctn_pr == "lp_cross"]) &&
         isTRUE(B$lp_pop$cross_border[B$lp_pop$lctn_pr == "lp_cross"]),
       "lp_cross is flagged cross_border in BOTH tasks")
expect(!any(A$lp_pop$cross_border[A$lp_pop$lctn_pr != "lp_cross"]),
       "no single-country LP is flagged cross_border")

# --- 3. gid populations Q_g --------------------------------------------------
cat("\n[3] Q_g and the ADM2-tiles-ADM1 identity\n")
qA <- A$gid_pop
expect(near(qA$Q[qA$gid == "AAA.1_1"], 10000),
       "Q(ADM1 cell) == area x density")
expect(near(sum(qA$Q[qA$gadm_level == "ADM2" & grepl("^AAA\\.1\\.", qA$gid)]),
            qA$Q[qA$gid == "AAA.1_1"]),
       "sum of ADM2 children == parent ADM1")
expect(near(sum(qA$Q[qA$gadm_level == "ADM1"]), 8 * 10000),
       "ADM1 units tile the country")
expect(near(A$manifest$natl_ratio, 1),
       "manifest natl_ratio == 1 against the WPP total")

# --- 4. w sums to 1; c <= 1 --------------------------------------------------
cat("\n[4] share algebra\n")
xwA <- deriveCrosswalkShares(A$crosswalk, A$gid_pop)

wsum <- xwA %>% group_by(lctn_pr, gadm_level) %>%
  summarise(s = sum(w), .groups = "drop")
expect(near(wsum$s, rep(1, nrow(wsum)), tol = 1e-10),
       sprintf("w sums to 1 for every (LP, level) — %d groups", nrow(wsum)))
expect(all(xwA$c <= 1 + 1e-10),
       "c[g<-l] <= 1 everywhere (an LP cannot hold more than the gid does)")
expect(near(xwA$c[xwA$lctn_pr == "lp_single" & xwA$gadm_level == "ADM1"], 1),
       "an LP that IS a gid has c == 1")

# --- 5. composites, with no string parsing -----------------------------------
cat("\n[5] composite detection is geometric\n")
cov1 <- computeGidCover(A$crosswalk) %>% filter(gadm_level == "ADM1")
g <- function(lp, col) cov1[[col]][cov1$lctn_pr == lp]

expect(g("lp_single", "n_gid") == 1L && near(g("lp_single", "gid_cover"), 1) &&
         near(g("lp_single", "n_gid_eff"), 1) && !g("lp_single", "is_composite"),
       "singleton LP: n_gid 1, gid_cover 1, n_gid_eff 1, not composite")
expect(g("lp_comp7", "n_gid") == 7L &&
         near(g("lp_comp7", "gid_cover"), 1 / 7, tol = 1e-9) &&
         near(g("lp_comp7", "n_gid_eff"), 7, tol = 1e-7) &&
         g("lp_comp7", "is_composite"),
       "7-way composite: gid_cover ~ 1/7, n_gid_eff ~ 7, flagged composite")
cov2 <- computeGidCover(A$crosswalk) %>% filter(gadm_level == "ADM2")
expect(near(cov2$gid_cover[cov2$lctn_pr == "lp_single"], 0.5) &&
         near(cov2$n_gid_eff[cov2$lctn_pr == "lp_single"], 2),
       "the same LP is a 2-way composite at ADM2 — cover is level-relative")
expect(!any(computeGidCover(A$crosswalk, cover_threshold = 0.05)$is_composite) &&
         all(computeGidCover(A$crosswalk, cover_threshold = 0.99)$is_composite[
           computeGidCover(A$crosswalk)$n_gid > 1L]),
       "cover_threshold is live: 0.05 declares nothing composite, 0.99 declares every multi-gid LP composite")

# --- 6. slivers are dropped by mass and w is renormalised --------------------
cat("\n[6] sliver handling\n")
sl <- A$crosswalk %>% filter(lctn_pr == "lp_sliver", gadm_level == "ADM1")
expect(nrow(sl) == 1L && sl$gid == "AAA.1_1",
       "the 1 m strip into cell 2 (H = 10 vs threshold 100.1) is dropped")
expect(near(sl$w_dropped, 10 / 10010, tol = 1e-9),
       "w_dropped records the discarded mass fraction (~1e-3)")
sl_w <- xwA %>% filter(lctn_pr == "lp_sliver", gadm_level == "ADM1")
expect(near(sum(sl_w$w), 1),
       "w is renormalised over the surviving pieces")

A_nofilter <- runPopulationTask("AAA", 2020L, gadm, lp_sf,
                                extract_fn = stub_extract,
                                sliver_mass_frac = 0, sliver_mass_min = 0,
                                require_boundary = FALSE, quiet = TRUE)
sl_all <- A_nofilter$crosswalk %>%
  filter(lctn_pr == "lp_sliver", gadm_level == "ADM1") %>% arrange(gid)
expect(nrow(sl_all) == 2L && near(sl_all$H, c(10000, 10)) &&
         all(sl_all$w_dropped == 0),
       "with the threshold off the strip survives — the filter is what dropped it")

# --- 7. reconciliation against the INDEPENDENT whole-LP extraction -----------
cat("\n[7] sum_g H reconciles with P_direct\n")
rec <- reconcileCrosswalk(A$crosswalk, A$lp_pop, tol = 0.01)
r <- function(lp, lvl, col) rec[[col]][rec$lctn_pr == lp & rec$gadm_level == lvl]

expect(near(r("lp_single", "ADM1", "P_l"), r("lp_single", "ADM1", "P_direct")),
       "clean LP: sum_g H == P_direct exactly")
expect(near(r("lp_comp7", "ADM1", "P_l"), 7 * 10000) &&
         near(r("lp_comp7", "ADM1", "P_direct"), 7 * 10000),
       "composite LP: both routes give 70000")
expect(near(r("lp_sliver", "ADM1", "rel_diff"), 10 / 10010, tol = 1e-9) &&
         r("lp_sliver", "ADM1", "ok"),
       "sliver LP misses by exactly the dropped mass, still inside 1% tol")
expect(all(rec$ok),
       "every LP in the AAA task reconciles within 1%")
expect(all(near(
  (A$crosswalk %>% filter(gadm_level == "ADM1") %>% group_by(lctn_pr) %>%
     summarise(s = sum(H), .groups = "drop"))$s,
  (A$crosswalk %>% filter(gadm_level == "ADM2") %>% group_by(lctn_pr) %>%
     summarise(s = sum(H), .groups = "drop"))$s, tol = 1e-6)),
       "ADM1 and ADM2 routes give the same P_l (levels tile independently)")

# --- 8. a cross-border LP is only whole once the tasks are bound -------------
cat("\n[8] cross-border LP binds across two tasks\n")
xw_bound  <- bind_rows(A$crosswalk, B$crosswalk)
gid_bound <- bind_rows(A$gid_pop, B$gid_pop)
lp_bound  <- bind_rows(A$lp_pop, B$lp_pop)

xc <- deriveCrosswalkShares(xw_bound, gid_bound) %>%
  filter(lctn_pr == "lp_cross", gadm_level == "ADM1") %>%
  arrange(gid)

expect(nrow(xc) == 2L && identical(xc$gid, c("AAA.8_1", "BBB.1_1")),
       "pieces from both countries are present after binding")
expect(near(xc$H, c(5000, 10000)),
       "each piece carries its OWN country's density (0.5 km2 x 10000 / x 20000)")
expect(near(unique(xc$P_l), 15000),
       "P_l = 15000 only once both tasks are bound")
expect(near(xc$w, c(1/3, 2/3), tol = 1e-9) && near(sum(xc$w), 1),
       "w is taken over the bound total, not per task")
expect(near(xc$c, c(0.5, 0.5)),
       "c = 0.5 in each gid — the LP covers half of each")
expect(near(unique(xc$gid_cover), 2/3, tol = 1e-9) &&
         near(unique(xc$n_gid_eff), 1.8, tol = 1e-9) &&
         all(xc$is_composite),
       "gid_cover 2/3, n_gid_eff 1.8, composite")

rec_b <- reconcileCrosswalk(xw_bound, lp_bound, tol = 0.01)
expect(near(rec_b$P_direct[rec_b$lctn_pr == "lp_cross" &
                             rec_b$gadm_level == "ADM1"], 15000),
       "P_direct sums the two country-clipped extractions to 15000")
expect(rec_b$ok[rec_b$lctn_pr == "lp_cross" & rec_b$gadm_level == "ADM1"],
       "the cross-border LP reconciles once bound")

# --- 9. failure modes are loud, not silent -----------------------------------
cat("\n[9] silent-corruption guards fire\n")

bad_len <- function(geom_sf, ...) rep(1, nrow(geom_sf) - 1L)
expect(inherits(try(runPopulationTask("AAA", 2020L, gadm, lp_sf,
                                      extract_fn = bad_len,
                                      require_boundary = FALSE, quiet = TRUE),
                    silent = TRUE), "try-error"),
       "a length mismatch from extract_fn is an error, not a recycled vector")

all_na <- function(geom_sf, ...) rep(NA_real_, nrow(geom_sf))
expect(inherits(try(runPopulationTask("AAA", 2020L, gadm, lp_sf,
                                      extract_fn = all_na,
                                      require_boundary = FALSE, quiet = TRUE),
                    silent = TRUE), "try-error"),
       "an all-NA extraction (raster never loaded) is an error, not zeros")

expect(inherits(try(assertBoundaryCached("AAA", tempfile()), silent = TRUE),
                "try-error"),
       "a boundary-cache miss is an error — the UN adjustment cannot be skipped silently")

expect(inherits(try(runPopulationTask("ZZZ", 2020L, gadm, lp_sf,
                                      extract_fn = stub_extract,
                                      require_boundary = FALSE, quiet = TRUE),
                    silent = TRUE), "try-error"),
       "a country with no GADM units is an error, not an empty crosswalk")

# --- 10. GEOMETRYCOLLECTION sanitising ---------------------------------------
cat("\n[10] non-polygon intersection debris is removed before extraction\n")
touch <- sf::st_sfc(list(box(8000, 9000)), crs = CRS)   # abuts AAA at x = 8000
mixed <- suppressWarnings(sf::st_intersection(
  sf::st_geometry(gadm[gadm$gid == "AAA.8_1", ]), touch))
san <- sanitizePolygons(mixed)
expect(all(san$empty),
       "an edge-only intersection sanitises to empty (it would extract as NA)")
expect(!nrow(buildCrosswalkPieces(
         sf::st_sf(lctn_pr = "t", geometry = touch),
         gadm[gadm$gid == "AAA.8_1", ])),
       "buildCrosswalkPieces() emits no row for a zero-area overlap")

pure <- sanitizePolygons(sf::st_geometry(gadm[1:2, ]))
expect(all(pure$repair == "none") && !any(pure$empty),
       "clean polygons are passed through untouched")

# --- 11. duplicate/parent-inherited geometry detection ------------------------
# Standalone tibble fixture, deliberately not built from the sf grid above: the
# defect this targets is
# purely an H-signature coincidence across (gid, year) — the geometry itself
# is irrelevant to detectGeomDupLps(), only the extracted numbers are.
cat("\n[11] duplicate-geometry LPs (detectGeomDupLps)\n")

mk_xw_dup <- function() {
  tibble::tribble(
    ~lctn_pr,     ~gadm_level, ~year, ~gid,      ~H,
    # lp_a and lp_b: byte-identical H at BOTH gids, BOTH years -> a genuine
    # duplicate cluster (mirrors the Conakry case: same polygon, two ids).
    "lp_a",       "ADM1",      2020L, "X.1_1",   1000,
    "lp_a",       "ADM1",      2021L, "X.1_1",   1010,
    "lp_b",       "ADM1",      2020L, "X.1_1",   1000,
    "lp_b",       "ADM1",      2021L, "X.1_1",   1010,
    # lp_c: same TOTAL as lp_a/lp_b but split across two gids -> not a match.
    "lp_c",       "ADM1",      2020L, "X.1_1",    600,
    "lp_c",       "ADM1",      2020L, "X.2_1",    400,
    "lp_c",       "ADM1",      2021L, "X.1_1",    606,
    "lp_c",       "ADM1",      2021L, "X.2_1",    404,
    # lp_d: genuinely unique, single year/gid.
    "lp_d",       "ADM1",      2020L, "X.3_1",   5000,
    # lp_a again at ADM2, but with DIFFERENT numbers than lp_b at ADM2 ->
    # the ADM1 match must not leak into a false ADM2 match.
    "lp_a",       "ADM2",      2020L, "X.1.1_1",  700,
    "lp_a",       "ADM2",      2020L, "X.1.2_1",  300,
    "lp_b",       "ADM2",      2020L, "X.1.1_1",  650,
    "lp_b",       "ADM2",      2020L, "X.1.2_1",  350
  )
}

dup <- detectGeomDupLps(mk_xw_dup())
d <- function(lp, lvl, col) dup[[col]][dup$lctn_pr == lp & dup$gadm_level == lvl]

expect(d("lp_a", "ADM1", "geom_dup_n") == 2L && d("lp_b", "ADM1", "geom_dup_n") == 2L,
       "lp_a and lp_b are flagged as a 2-member cluster at ADM1")
expect(identical(d("lp_a", "ADM1", "geom_dup_group"),
                 d("lp_b", "ADM1", "geom_dup_group")) &&
         !is.na(d("lp_a", "ADM1", "geom_dup_group")),
       "they share the same non-NA geom_dup_group label")
expect(d("lp_a", "ADM1", "geom_dup_group") == "ADM1:lp_a",
       "the group label picks the lexicographically lowest lctn_pr, per house convention")
expect(d("lp_c", "ADM1", "geom_dup_n") == 1L && is.na(d("lp_c", "ADM1", "geom_dup_group")),
       "lp_c has the same P_l as lp_a/lp_b but a different gid split -> not flagged")
expect(d("lp_d", "ADM1", "geom_dup_n") == 1L && is.na(d("lp_d", "ADM1", "geom_dup_group")),
       "a genuinely unique LP is not flagged")
expect(d("lp_a", "ADM2", "geom_dup_n") == 1L && is.na(d("lp_a", "ADM2", "geom_dup_group")),
       "the same LP pair diverges at ADM2 and is correctly NOT flagged at that level")

empty_dup <- detectGeomDupLps(mk_xw_dup()[0, ])
expect(nrow(empty_dup) == 0L,
       "an empty crosswalk returns an empty (not erroring) result")

# --- 12. per-country runner == per-(country, year) tasks ---------------------
# runPopulationCountryTask() builds the geometry once and extracts per year; a
# year's tables must be the per-year task's. With Q_g cached, only the LP roles
# reach the extractor and the tables are unchanged.
cat("\n[12] runPopulationCountryTask() matches runPopulationTask() per year (serial and on 2 workers)\n")
strip_rt <- function(r) {
  r$manifest[c("runtime_sec", "built_at", "gid_pop_cached")] <- NULL
  r
}
YEARS <- c(2020L, 2021L)
seen_roles <- new.env(parent = emptyenv()); seen_roles$n <- integer(0)
count_stub <- function(geom_sf, ...) {
  seen_roles$n <- c(seen_roles$n, nrow(geom_sf))
  stub_extract(geom_sf, ...)
}
for (k in c("AAA", "BBB")) {
  per_year <- lapply(YEARS, function(y) runPopulationTask(
    k, y, gadm, lp_sf, extract_fn = stub_extract, sliver_mass_frac = 1e-2,
    wpp_total = 8e4, require_boundary = FALSE, quiet = TRUE))
  seen_roles$n <- integer(0)
  cold <- runPopulationCountryTask(
    k, YEARS, gadm, lp_sf, extract_fn = count_stub, sliver_mass_frac = 1e-2,
    wpp_totals = rep(8e4, 2), require_boundary = FALSE, quiet = TRUE)
  n_cold <- seen_roles$n
  expect(identical(names(cold), as.character(YEARS)),
         sprintf("%s: one result per year, named by year", k))
  expect(identical(lapply(cold, strip_rt), setNames(lapply(per_year, strip_rt), YEARS)),
         sprintf("%s: cold cache gives the per-year tables exactly", k))
  expect(!any(vapply(cold, function(r) r$manifest$gid_pop_cached, logical(1))),
         sprintf("%s: cold cache is recorded as such", k))

  seen_roles$n <- integer(0)
  warm <- runPopulationCountryTask(
    k, YEARS, gadm, lp_sf, extract_fn = count_stub, sliver_mass_frac = 1e-2,
    wpp_totals = rep(8e4, 2), require_boundary = FALSE, quiet = TRUE,
    gid_pop_cache = lapply(cold, `[[`, "gid_pop"))
  n_gadm <- sum(gadm$gid_iso3 == k)
  expect(identical(lapply(warm, strip_rt), lapply(cold, strip_rt)),
         sprintf("%s: warm Q_g cache gives the same tables", k))
  expect(identical(seen_roles$n, n_cold - n_gadm),
         sprintf("%s: warm cache extracts the LP roles only", k))
  expect(all(vapply(warm, function(r) r$manifest$gid_pop_cached, logical(1))),
         sprintf("%s: warm cache is recorded in the manifest", k))
}

# workers > 1 extracts the years on forked processes; the tables must not move.
par2 <- runPopulationCountryTask(
  "AAA", YEARS, gadm, lp_sf, extract_fn = stub_extract, sliver_mass_frac = 1e-2,
  wpp_totals = rep(8e4, 2), require_boundary = FALSE, quiet = TRUE, workers = 2L)
ser1 <- runPopulationCountryTask(
  "AAA", YEARS, gadm, lp_sf, extract_fn = stub_extract, sliver_mass_frac = 1e-2,
  wpp_totals = rep(8e4, 2), require_boundary = FALSE, quiet = TRUE, workers = 1L)
expect(identical(names(par2), as.character(YEARS)) &&
         identical(lapply(par2, strip_rt), lapply(ser1, strip_rt)),
       "workers = 2 gives the serial tables, named by year")
expect(identical(future::plan("list")[[1]], future::plan("list")[[1]]) &&
         inherits(future::plan(), "sequential"),
       "the future plan is restored after a parallel country task")

seen_roles$n <- integer(0)
lone <- sf::st_sf(lctn_pr = "far", geometry = sf::st_sfc(box(50000, 51000), crs = CRS))
no_lp <- runPopulationCountryTask(
  "BBB", 2020L, gadm, lone, extract_fn = count_stub,
  require_boundary = FALSE, quiet = TRUE,
  gid_pop_cache = list(`2020` = cold[["2020"]]$gid_pop))
expect(!length(seen_roles$n) && !nrow(no_lp[["2020"]]$lp_pop),
       "no LP in the country and Q_g cached: no extraction call at all")

# An LP-only call that is all NA (only unpopulated pieces) is not a raster
# failure when Q_g is cached: the year is redone with every role.
na_lp_stub <- function(geom_sf, ...) {
  out <- stub_extract(geom_sf, ...)
  if (!any(geom_sf$role %in% c("gadm1", "gadm2"))) out[] <- NA_real_
  out
}
redo <- runPopulationCountryTask(
  "AAA", 2020L, gadm, lp_sf, extract_fn = na_lp_stub, sliver_mass_frac = 1e-2,
  wpp_totals = 8e4, require_boundary = FALSE, quiet = TRUE,
  gid_pop_cache = list(`2020` = runPopulationTask("AAA", 2020L, gadm, lp_sf,
    extract_fn = stub_extract, sliver_mass_frac = 1e-2, require_boundary = FALSE,
    quiet = TRUE)$gid_pop))
expect(!isTRUE(redo[["2020"]]$manifest$gid_pop_cached),
       "an all-NA LP-only extraction with Q_g cached falls back to the full call")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
