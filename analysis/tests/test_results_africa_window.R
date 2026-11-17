# test_results_africa_window.R — HARD: one definition of Africa for the
#   attribution steps.
#
# Usage:
#   Rscript analysis/tests/test_results_africa_window.R
#
# Details:
#   Quick-safe: no fit, no simulation output (reads only a committed
#   forecast manifest, a config file). Pins results_derive.R's
#   results_africa_window(), attrib_cell_table(), africa_scope_table() and
#   sim_io.R::manifest_last_firm_oni(), the cells on a synthetic run with
#   continuation months, loader ends and the stamp back-compat.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
D <- as.Date
err <- function(e) {
  inherits(try(suppressMessages(e), silent = TRUE), "try-error")
}
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(a, b, tolerance = tol))

# [1] manifest_last_firm_oni(): firm_oni_end wins, else last_cfsv2_month.
cat("# [1] manifest_last_firm_oni()\n")
expect(
  identical(
    manifest_last_firm_oni(list(last_cfsv2_month = "2027-06-01")),
    "2027-06-01"
  ),
  "without firm_oni_end: last_cfsv2_month"
)
expect(
  identical(
    manifest_last_firm_oni(list(
      last_cfsv2_month = "2027-06-01",
      firm_oni_end = "2027-05-01"
    )),
    "2027-05-01"
  ),
  "firm_oni_end wins (tier-2 item 16)"
)
expect(is.null(manifest_last_firm_oni(list())), "neither: NULL")

# [2] results_africa_window() on the committed manifest: legacy (=
#     africa_min_horizon), the common window (2027-10, to_end), firm
#     only (0), an integer lead, max_counterfactual_winter, and
#     a manifest with firm_oni_end (South one month short of the end
#     gets one continuation month at lead 1). Refusals.
cat("# [2] results_africa_window()\n")
man <- yaml::read_yaml("analysis/tests/fixtures/forecast_horizons_fixture.yml")
groups <- manifest_groups_inc(man, "enso_central")
G5 <- c("enso_east", "enso_horn", "enso_nga", "enso_south", "enso_west")
expect(setequal(groups, G5), "fixture: five groups")
lg <- results_africa_window(
  man,
  groups,
  list(window_end = "legacy", continuation_max_lead = 0L)
)
expect(
  lg$mode == "legacy" &&
    lg$end == D("2026-10-01") &&
    all(lg$group_ends == D("2026-10-01")) &&
    !lg$extended &&
    !length(lg$censored),
  "legacy: every group to africa_min_horizon, nothing extended"
)
expect(
  identical(results_africa_window(man, groups, NULL)$group_ends, lg$group_ends),
  "NULL block = legacy"
)
af <- list(
  window_end = "2027-10-01",
  continuation_max_lead = "to_end",
  max_counterfactual_winter = ""
)
cw <- results_africa_window(man, groups, af)
expect(
  cw$mode == "common" &&
    all(cw$group_ends == D("2027-10-01")) &&
    cw$extended &&
    !length(cw$censored),
  "2027-10 to_end: every group to 2027-10, none censored"
)
expect(
  identical(
    unname(cw$lead[c(
      "enso_west",
      "enso_east",
      "enso_south",
      "enso_horn",
      "enso_nga"
    )]),
    c(12L, 5L, 0L, 0L, 0L)
  ),
  "months past the firm horizon: West 12, East 5, others 0"
)
expect(
  cw$firm_oni_end == D("2027-06-01"),
  "last firm ONI month from the manifest (last_cfsv2_month)"
)
f0 <- results_africa_window(
  man,
  groups,
  modifyList(af, list(continuation_max_lead = 0L))
)
expect(
  f0$group_ends[["enso_west"]] == D("2026-10-01") &&
    f0$group_ends[["enso_east"]] == D("2027-05-01") &&
    setequal(f0$censored, c("enso_west", "enso_east")) &&
    !f0$extended,
  "continuation_max_lead 0: firm horizons, West and East censored"
)
f2 <- results_africa_window(
  man,
  groups,
  modifyList(af, list(continuation_max_lead = 2L))
)
expect(
  f2$group_ends[["enso_west"]] == D("2026-12-01") &&
    f2$group_ends[["enso_east"]] == D("2027-07-01") &&
    f2$lead[["enso_west"]] == 2L,
  "integer lead: two months past each firm horizon"
)
d27 <- results_africa_window(
  man,
  groups,
  modifyList(af, list(max_counterfactual_winter = "2026/27"))
)
expect(
  d27$group_ends[["enso_west"]] == D("2027-03-01") &&
    d27$group_ends[["enso_east"]] == D("2027-10-01") &&
    d27$group_ends[["enso_south"]] == D("2027-10-01") &&
    identical(d27$censored, "enso_west"),
  "max_counterfactual_winter 2026/27 (D2.7): West capped at 2027-03 (mon reads t + 8 <= 2027-11), others unchanged"
)
expect(
  err(results_africa_window(
    man,
    groups,
    modifyList(af, list(window_end = "2027-10-15"))
  )),
  "refused: not a first of month"
)
expect(
  err(results_africa_window(
    man,
    groups,
    modifyList(af, list(continuation_max_lead = -1))
  )),
  "refused: negative lead"
)
expect(
  err(results_africa_window(
    man,
    groups,
    modifyList(af, list(continuation_max_lead = "all"))
  )),
  "refused: unknown lead word"
)
expect(
  err(results_africa_window(
    man,
    groups,
    modifyList(af, list(max_counterfactual_winter = "2026/28"))
  )),
  "refused: malformed winter"
)
expect(
  err(results_africa_window(
    man,
    groups,
    modifyList(af, list(window_end = "2026-03-01"))
  )),
  "refused: the window ends before the attribution start"
)
# Manifest with firm_oni_end: firm horizons tied to the last full ONI
# month; South's firm horizon one month short of the end.
m19 <- man
m19$firm_oni_end <- "2027-05-01"
m19$groups$enso_south$horizon_end <- "2027-09-01"
c19 <- results_africa_window(m19, groups, af)
expect(
  c19$firm_oni_end == D("2027-05-01") &&
    c19$lead[["enso_south"]] == 1L &&
    c19$group_ends[["enso_south"]] == D("2027-10-01"),
  "firm_oni_end: South (firm 2027-09) gets one continuation month to the end"
)

# [3] africa_scope_table() (S28): longest continuation lead per track from
#     the manifest's lags and last firm ONI month, the CRPSS there, the
#     counterfactual winters read.
cat("# [3] africa_scope_table() (S28)\n")
skill <- tibble(lead = 1:24, crpss_vs_clim = round(1 / (1:24), 4))
s <- africa_scope_table(cw, man, D("2026-05-01"), skill, "2026/27")
w <- filter(s, enso_group == "enso_west")
e <- filter(s, enso_group == "enso_east")
so <- filter(s, enso_group == "enso_south")
expect(
  w$lead_monthly == 12L &&
    w$lead_occurrence == 10L &&
    w$lead_duration == 8L &&
    w$lead_max == 12L,
  "West: continuation leads mon 12 / occ 10 / dur 8 (computed from the lags)"
)
expect(
  e$lead_duration == 5L && e$lead_occurrence == 4L && e$lead_monthly == 2L,
  "East: dur 5 / occ 4 / mon 2"
)
expect(
  so$lead_max == 0L && is.na(so$crpss_at_lead_max),
  "South: no continuation, no CRPSS"
)
expect(
  near(w$crpss_at_lead_max, skill$crpss_vs_clim[12]) &&
    near(e$crpss_duration, skill$crpss_vs_clim[5]),
  "CRPSS read at each lead"
)
expect(
  w$cf_winters_read == "2026/27, 2027/28" &&
    w$cf_winters_unconstrained == "2027/28",
  "West reads DJF 2027/28, unconstrained in the counterfactual"
)
expect(
  e$cf_winters_unconstrained == "",
  "East reads only the constrained DJF 2026/27"
)
s19 <- africa_scope_table(c19, m19, D("2026-05-01"), NULL, NA_character_)
expect(
  filter(s19, enso_group == "enso_west")$lead_monthly == 13L &&
    filter(s19, enso_group == "enso_south")$lead_occurrence == 1L &&
    is.na(filter(s19, enso_group == "enso_west")$crpss_monthly),
  "firm_oni_end: leads counted from firm_oni_end (one more month); no skill file -> no CRPSS"
)
expect(
  identical(djf_winters(D("2026-11-01"), D("2027-01-01")), "2026/27") &&
    identical(djf_winters(D("2027-03-01"), D("2027-11-01")), character(0)),
  "djf_winters()"
)

# T-int 7(a): cf_constrained_winter() reads every classified winter of a
# library cache (set_definition$winters); a legacy cache keeps its one primary.
cat("# [3b] cf_constrained_winter(): legacy vs library caches\n")
td <- tempfile("cfw_")
dir.create(td)
cf_man <- function(cache_meta) {
  cp <- file.path(td, paste0(digest::digest(cache_meta), ".rds"))
  saveRDS(list(meta = cache_meta), cp)
  cfg <- file.path(td, paste0(basename(cp), ".yml"))
  yaml::write_yaml(
    list(simulator = list(enso = list(counterfactual = list(cache_path = cp)))),
    cfg
  )
  list(
    groups = list(
      enso_west = list(iso3 = list("NER"), counterfactual_config = cfg)
    )
  )
}
leg <- cf_man(list(set_definition = list(peak_season_year = 2026L)))
lib <- cf_man(list(
  set_definition = list(
    peak_season_year = 2026L,
    winters = list(
      primary = list(2026L),
      other = list(2027L),
      other_winters = "neutral"
    )
  )
))
lib_any <- cf_man(list(
  set_definition = list(
    peak_season_year = 2026L,
    winters = list(
      primary = list(2026L),
      other = list(2027L),
      other_winters = "any"
    )
  )
))
expect(
  identical(cf_constrained_winter(leg), "2026/27"),
  "legacy cache: the one primary winter (unchanged)"
)
expect(
  identical(cf_constrained_winter(lib), c("2026/27", "2027/28")),
  "library cache (other winters neutral): primary + other read winters"
)
expect(
  identical(cf_constrained_winter(lib_any), "2026/27"),
  "library cache, other_winters any: primary only"
)
expect(is.na(cf_constrained_winter(list(groups = list()))), "no unit: NA")
sl <- africa_scope_table(
  cw,
  man,
  D("2026-05-01"),
  skill,
  cf_constrained_winter(lib)
)
wl <- filter(sl, enso_group == "enso_west")
expect(
  wl$cf_winter_constrained == "2026/27, 2027/28" &&
    wl$cf_winters_unconstrained == "",
  "S28: with DJF 2027/28 constrained West reads no unconstrained winter"
)
expect(
  identical(
    s,
    africa_scope_table(
      cw,
      man,
      D("2026-05-01"),
      skill,
      cf_constrained_winter(leg)
    )
  ),
  "S28 from a legacy cache = the single-winter table"
)
expect(
  is.na(filter(s19, enso_group == "enso_west")$cf_winter_constrained),
  "NA winter stays NA"
)
unlink(td, recursive = TRUE)

# [4] attrib_cell_table(): legacy = the two former windows; common = + the
#     attribution_common / post_fork_common / attribution_firm cells, group
#     ends, censored groups.
cat("# [4] attrib_cell_table()\n")
a0 <- D("2026-05-01")
pf <- D("2026-09-01")
ct0 <- attrib_cell_table(a0, pf, africa_win = lg, cap = 2e5)
expect(
  identical(ct0, attrib_cell_table(a0, pf, cap = 2e5)),
  "legacy window: the two former windows only"
)
ct <- attrib_cell_table(a0, pf, africa_win = cw, cap = 2e5)
expect(
  identical(
    unique(ct$window),
    c(
      "attribution",
      "post_fork",
      "attribution_common",
      "post_fork_common",
      "attribution_firm"
    )
  ) &&
    nrow(ct) == 10L,
  "common: five windows x two modes"
)
ge <- function(wn) ct$group_ends[[match(wn, ct$window)]]
expect(
  identical(ge("attribution"), cw$firm) &&
    all(ge("attribution_common") == D("2027-10-01")) &&
    ge("attribution_firm")[["enso_horn"]] == D("2027-10-01") &&
    ge("attribution_firm")[["enso_west"]] == D("2026-10-01"),
  "group ends: legacy cells cut at the firm horizons, common at 2027-10, firm at min(firm, end)"
)
expect(
  ct$censored_groups[match("attribution_firm", ct$window)] ==
    "enso_east,enso_west" &&
    ct$censored_groups[match("attribution_common", ct$window)] == "",
  "censored groups of attribution_firm: East and West"
)
expect(
  ct$start[match("post_fork_common", ct$window)] == pf,
  "post_fork_common starts at the factual launch"
)
expect(
  nrow(attrib_cell_table(a0, NULL, africa_win = cw, cap = 2e5)) == 6L,
  "no post_fork start: no post_fork cells"
)
expect(
  identical(
    event_cells(list(AFRICA_WIN = cw, EXPLOSIVE_CAP = 2e5), a0, "ev")$window,
    c("ev", "ev")
  ),
  "event cells ignore the Africa window"
)

# [5] the cells on a synthetic 2-group run with continuation months: the
#     common cell's continent = the sum of the groups over the common window,
#     region rows sum to it, attribution_firm drops the continuation, the
#     legacy cells equal those of a load without continuation months.
cat("# [5] the cells on a synthetic run with continuation months\n")
mm <- list(
  groups = list(
    A = list(
      horizon_end = "2026-06-01",
      best_lag_occ = -2L,
      best_lag_dur = -1L,
      best_lag_mon = -2L,
      min_lag = -2L,
      iso3 = "AAA"
    ),
    B = list(
      horizon_end = "2026-08-01",
      best_lag_occ = 0L,
      best_lag_dur = 0L,
      best_lag_mon = 0L,
      min_lag = 0L,
      iso3 = "BBB"
    )
  ),
  africa_min_horizon = "2026-06-01",
  last_cfsv2_month = "2026-08-01",
  attribution_start = "2026-05-01"
)
awm <- results_africa_window(
  mm,
  c("A", "B"),
  af <- list(window_end = "2026-08-01", continuation_max_lead = "to_end")
)
expect(
  awm$group_ends[["A"]] == D("2026-08-01") &&
    awm$lead[["A"]] == 2L &&
    awm$extended,
  "A read 2 months past its firm horizon"
)
set.seed(3)
months <- seq(a0, by = "month", length.out = 4)
am <- tidyr::expand_grid(
  tibble(run = c("forecast", "counterfactual"), analog_set = c(NA, "neutral")),
  adm1_unit = c("AAA_1", "AAA_2", "BBB_1"),
  replicate_id = 1:30,
  month = months
) %>%
  mutate(
    enso_group = ifelse(substr(adm1_unit, 1, 3) == "AAA", "A", "B"),
    group_horizon_end = D(ifelse(
      enso_group == "A",
      "2026-06-01",
      "2026-08-01"
    )),
    cases = as.numeric(rpois(n(), ifelse(run == "forecast", 20, 12))),
    outbreaks = rbinom(n(), 1, .3)
  )
mk_cm <- function(a) {
  a %>%
    mutate(country = substr(adm1_unit, 1, 3)) %>%
    group_by(
      replicate_id,
      month,
      country,
      run,
      analog_set,
      enso_group,
      group_horizon_end
    ) %>%
    summarise(
      n_cases = sum(cases),
      n_outbreaks = sum(outbreaks),
      .groups = "drop"
    )
}
cm <- mk_cm(am)
pop <- tibble(adm1_unit = c("AAA_1", "AAA_2", "BBB_1"), pop = c(100, 200, 300))
gs <- c(A = "Ay", B = "Be")
cells <- attrib_cell_table(a0, D("2026-06-01"), africa_win = awm, cap = 2e5)
mw <- suppressMessages(attrib_metrics_windows(
  am,
  cm,
  pop,
  gs,
  "neutral",
  cells,
  50,
  .05,
  c(1, 2),
  n_boot = 0L
))
sc <- function(wn, lv) {
  filter(
    mw$af_scope,
    window == wn,
    explosive_guard,
    metric == "cases",
    level == lv
  )
}
tot <- function(d) {
  d %>%
    group_by(run, replicate_id) %>%
    summarise(v = sum(cases), .groups = "drop") %>%
    group_by(run) %>%
    summarise(m = mean(v), .groups = "drop")
}
ref_common <- tot(am) # every group to 2026-08
ref_firm <- tot(filter(am, month <= group_horizon_end)) # A cut at its firm horizon
afr <- sc("attribution_common", "continent")
expect(
  near(afr$mean_f, ref_common$m[ref_common$run == "forecast"]) &&
    near(afr$mean_cf, ref_common$m[ref_common$run == "counterfactual"]),
  "attribution_common: Africa = every group summed over the common window (A's continuation included)"
)
expect(
  near(sum(sc("attribution_common", "region")$mean_f), afr$mean_f),
  "region rows sum to the continent row"
)
fi <- sc("attribution_firm", "continent")
expect(
  near(fi$mean_f, ref_firm$m[ref_firm$run == "forecast"]),
  "attribution_firm: A cut at its firm horizon"
)
expect(
  near(sc("attribution", "continent")$mean_f, fi$mean_f),
  "legacy attribution = each group over its firm window"
)
mw0 <- suppressMessages(attrib_metrics_windows(
  filter(am, month <= group_horizon_end),
  mk_cm(filter(am, month <= group_horizon_end)),
  pop,
  gs,
  "neutral",
  attrib_cell_table(a0, D("2026-06-01"), cap = 2e5),
  50,
  .05,
  c(1, 2),
  n_boot = 0L
))
expect(
  identical(
    filter(mw$af_scope, window %in% c("attribution", "post_fork")),
    mw0$af_scope
  ) &&
    identical(
      filter(mw$pr_adm1, window %in% c("attribution", "post_fork")),
      mw0$pr_adm1
    ),
  "legacy cells on a load with continuation months == on a firm-only load"
)
ae <- suppressMessages(attrib_exposure_windows(
  am,
  cm,
  pop,
  gs,
  "neutral",
  cells,
  0.5,
  .66,
  .8,
  n_boot = 0L
))
expect(
  setequal(unique(ae$exposure$window), unique(cells$window)),
  "attrib_exposure_windows: every cell"
)

# [6] loader ends and the stamp back-compat of the new keys.
cat("# [6] loader ends, stamp back-compat\n")
ctx <- list2env(list(AFRICA_WIN = awm, MAN = mm))
expect(
  identical(
    .africa_load_ends(ctx, mm, c("A", "B")),
    c(A = D("2026-08-01"), B = D("2026-08-01"))
  ),
  "load ends: A to its end in the window, B to its firm horizon"
)
expect(
  is.null(.africa_load_ends(list2env(list(AFRICA_WIN = lg)), man, groups)),
  "legacy: no extension"
)
mm2 <- mm
mm2$groups$A$horizon_end <- "2026-07-01"
expect(
  err(.africa_load_ends(ctx, mm2, c("A", "B"))),
  "refused: Africa window of another manifest"
)
defs <- get_default_options()$results
st <- list(results = defs[.RESULTS_STAMPED_BLOCKS])
st$results$scope$africa <- NULL
st$results$inputs$enso_continuation_skill <- NULL
expect(
  !length(.results_decision_diff(st, defs)),
  "an artefact stamped before the keys is fresh at their defaults"
)
res_c <- defs
res_c$scope$africa$window_end <- "2027-10-01"
expect(
  length(.results_decision_diff(st, res_c)) == 1L,
  "a common window re-derives it"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
