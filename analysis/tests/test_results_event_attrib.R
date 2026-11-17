# test_results_event_attrib.R — HARD: report_event_attrib
#   (utils/results_derive.R), attribution of past ENSO events.
#
# Usage:
#   Rscript analysis/tests/test_results_event_attrib.R
#
# Details:
#   Quick-safe: synthetic temp files only, no fit, no simulation. Pins the
#   smoothing-as-factual reading, the event manifest view, the per-event
#   cells, observed totals and AN_obs, the analytic crossbasis AF core and
#   the past-event stamp.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(qs2)
})
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
near <- function(a, b, tol = 1e-9) {
  isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
}
cat("# report_event_attrib\n")

td <- tempfile("evattrib_")
dir.create(td)
SIM <- file.path(td, "sim")

# --- [1] smoothing-as-factual ------------------------------------------------
# [1] smoothing-as-factual reading (load_arm_run(factual_source =
#     "smoothing")): every month of the factual run comes from the smoothing
#     replicate each replicate was handed off from (handoff_map.csv), on the
#     run's own month grid, zero-filled; the counterfactual arm is read
#     plainly; the default reading is unchanged.
cat("\n[1] smoothing-as-factual reading\n")
months <- seq(as.Date("2015-01-01"), as.Date("2015-03-01"), by = "month")
sm_pd <- file.path(td, "sm", "parts", "tu=eh", "chunk_0001")
dir.create(sm_pd, recursive = TRUE)
ev <- tibble(
  rep_id = c(1L, 2L, 2L),
  t_start = as.Date(c("2015-01-01", "2015-02-01", "2014-12-01")),
  gid = c("KEN_1_1", "ETH_2_1", "KEN_1_1"),
  country = c("KEN", "ETH", "KEN"),
  Yout = 1L
)
ev_mon <- tibble(
  rep_id = c(1L, 1L, 2L, 2L),
  t_start = as.Date(c("2015-01-01", "2015-02-01", "2015-02-01", "2014-12-01")),
  gid = c("KEN_1_1", "KEN_1_1", "ETH_2_1", "KEN_1_1"),
  country = c("KEN", "KEN", "ETH", "KEN"),
  cases_month = c(40, 25, 7, 999)
)
qs_save(list(events = ev, events_mon = ev_mon), file.path(sm_pd, "part.qs"))
invisible(file.create(file.path(sm_pd, "_SUCCESS")))
fpart <- function(tag, map = TRUE) {
  pd <- file.path(SIM, tag, "parts", "tu=eh")
  dir.create(pd, recursive = TRUE)
  ms <- tidyr::expand_grid(
    replicate_id = 1:3,
    gid = c("KEN_1_1", "ETH_2_1"),
    month = format(months)
  ) %>%
    mutate(Y = 0L, cases_month = 1)
  qs_save(ms, file.path(pd, "monthly_summary.qs"))
  invisible(file.create(file.path(pd, "_SUCCESS")))
  if (map) {
    utils::write.csv(
      data.frame(
        tu = "eh",
        replicate_id = 1:3,
        smoothing_rep_id = c(1L, 2L, 1L),
        handoff_month = "2015-01-01",
        handoff_weighting = "uniform"
      ),
      file.path(pd, "handoff_map.csv"),
      row.names = FALSE
    )
  }
}
fpart("fc")
fpart("cf")
sev <- read_smoothing_parts_events(
  file.path(td, "sm", "parts"),
  as.Date("2011-01-01")
)
starts <- list(factual = months[1], counterfactual = months[1])
sm <- suppressMessages(load_arm_run(
  SIM,
  "fc",
  "forecast",
  months[1],
  months[3],
  arm_starts = starts,
  smoothing_events = sev,
  factual_source = "smoothing"
))
plain <- load_arm_run(
  SIM,
  "fc",
  "forecast",
  months[1],
  months[3],
  arm_starts = starts
)
cm <- sm$country_month
cell <- function(rid, ctry, m, col) {
  cm[[col]][
    cm$replicate_id == rid & cm$country == ctry & cm$month == as.Date(m)
  ]
}
expect(
  identical(cell(1L, "KEN", "2015-01-01", "n_cases"), 40) &&
    identical(cell(1L, "KEN", "2015-02-01", "n_cases"), 25) &&
    identical(cell(3L, "KEN", "2015-01-01", "n_cases"), 40) &&
    identical(cell(2L, "ETH", "2015-02-01", "n_cases"), 7) &&
    identical(cell(1L, "KEN", "2015-01-01", "n_outbreaks"), 1L) &&
    identical(cell(2L, "ETH", "2015-02-01", "n_outbreaks"), 1L),
  "every factual month = the paired smoothing replicate (one smoothing replicate may feed two)"
)
expect(
  sum(cm$n_cases) == 40 + 25 + 40 + 25 + 7,
  "no simulated forecast month survives; out-of-window smoothing events excluded"
)
expect(
  nrow(cm) == 3L * 2L * 3L && identical(sort(unique(cm$month)), months),
  "full replicate x country x month grid on the run's own months"
)
expect(
  identical(names(cm), names(plain$country_month)) &&
    identical(names(sm$adm1_month), names(plain$adm1_month)),
  "same columns as the plain read"
)
expect(
  isTRUE(sm$stitch$replace) && identical(sm$stitch$months, months),
  "stitch element marks the replacement"
)
expect(
  identical(
    load_arm_run(
      SIM,
      "cf",
      "counterfactual",
      months[1],
      months[3],
      arm_starts = starts,
      smoothing_events = sev,
      factual_source = "smoothing"
    ),
    load_sim_run_monthly(SIM, "cf", months[1], months[3])
  ),
  "the counterfactual arm is read plainly whatever the factual source"
)
expect(
  sum(plain$country_month$n_cases) == 3 * 2 * 3,
  "default factual_source: the simulated run"
)
e <- tryCatch(
  {
    load_arm_run(
      SIM,
      "fc",
      "forecast",
      months[1],
      months[3],
      factual_source = "smoothing"
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e) && grepl("smoothing events", e),
  "smoothing source without events is refused"
)

# --- [2] manifest view ------------------------------------------------------
# [2] event_manifest_view(): a v2-shaped manifest, both arms launching at the
#     event launch, firm horizons per country; read_event_manifest() refuses
#     another manifest kind.
cat("\n[2] event_manifest_view()\n")
em <- list(
  schema = 2L,
  kind = "event_attribution",
  analog_set = "neutral",
  grouping = "g",
  smoothing_input = "x/cache",
  events = list(
    elnino_2015 = list(
      launch = "2014-09-01",
      peak = "2015-12-01",
      africa_max_horizon = "2017-07-01",
      groups = list(
        enso_east = list(horizon_end = "2016-09-01", iso3 = list("KEN")),
        enso_horn = list(horizon_end = "2017-07-01", iso3 = list("ETH"))
      ),
      units = list(
        eh = list(
          members = list("KEN", "ETH"),
          groups = list("enso_east", "enso_horn"),
          horizon_end = "2017-07-01",
          factual_config = "f.yml",
          factual_launch = "2014-09-01",
          analog_sets = list("neutral"),
          primary_analog_set = "neutral",
          counterfactual_sets = list(
            neutral = list(config = "c.yml", launch = "2014-09-01")
          )
        )
      ),
      countries = list(
        KEN = list(
          group = "enso_east",
          unit = "eh",
          firm_horizon_end = "2016-09-01"
        ),
        ETH = list(
          group = "enso_horn",
          unit = "eh",
          firm_horizon_end = "2017-07-01"
        )
      )
    )
  )
)
v <- event_manifest_view(em, "elnino_2015")
st <- manifest_arm_starts(v, c("enso_east", "enso_horn"), "neutral")
expect(
  st$factual == as.Date("2014-09-01") &&
    st$counterfactual == as.Date("2014-09-01"),
  "both arms launch at the event launch"
)
expect(
  manifest_attrib_start(v, NULL, "neutral", "stitch") == as.Date("2014-09-01"),
  "attribution window starts at the launch"
)
expect(
  identical(
    unname(manifest_country_horizons(v, "firm")[c("KEN", "ETH")]),
    c("2016-09-01", "2017-07-01")
  ),
  "firm horizons = each group's event window end"
)
expect(
  identical(v$africa_min_horizon, "2016-09-01"),
  "africa_min_horizon = earliest group window end"
)
mp <- file.path(td, "not_event.yml")
yaml::write_yaml(list(schema = 2L, groups = list()), mp)
expect(
  inherits(try(read_event_manifest(mp), silent = TRUE), "try-error"),
  "a forecast manifest is refused"
)
expect(
  inherits(
    try(event_manifest_view(em, "lanina_2020"), silent = TRUE),
    "try-error"
  ),
  "an unknown event is refused"
)

# --- [3] per-event cells ----------------------------------------------------
# [3] event_attrib_cells(): both explosive guards, the attrib_metrics_tables()
#     numbers unchanged.
cat("\n[3] event_attrib_cells()\n")
a0 <- as.Date("2026-05-01")
am <- expand.grid(
  replicate_id = 1:4,
  adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
  run = c("forecast", "counterfactual"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    month = a0,
    analog_set = ifelse(run == "forecast", NA_character_, "neutral"),
    enso_group = ifelse(
      substr(adm1_unit, 1, 3) == "KEN",
      "enso_east",
      "enso_horn"
    ),
    group_horizon_end = as.Date("2026-06-01"),
    cases = case_when(
      run == "counterfactual" ~ 10,
      adm1_unit == "KEN_1" ~ 30,
      TRUE ~ 12
    ),
    outbreaks = as.integer(cases > 10)
  ) %>%
  as_tibble()
am$cases[
  am$run == "forecast" & am$replicate_id == 4 & am$adm1_unit == "KEN_1"
] <- 2e5
cmf <- am %>%
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
inp <- list(
  adm1_month = am,
  country_month = cmf,
  attrib_start = a0,
  adm1_pop = tibble(
    adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
    pop = c(100, 40, 60)
  )
)
ctx <- new.env()
ctx$GROUP_SHORT_RUN <- c(enso_east = "East", enso_horn = "Horn")
ctx$PR_CASE_THRESHOLD <- 20
ctx$PR_ALPHA <- 0.05
ctx$SN_LEVELS <- c(1, 2)
ctx$ATTRIB_MAGNITUDES <- c(0.25, 0.5, 1.0)
ctx$ROBUST_EXCEED <- 0.66
ctx$ROBUST_AGREE <- 0.80
ctx$EXPLOSIVE_CAP <- 2e5
ctx$EXPLOSIVE_MODES <- c("drop", "none")
cl <- suppressMessages(event_attrib_cells(inp, ctx, "neutral"))
af <- filter(cl$af_adm1, metric == "cases", adm1_unit == "KEN_1")
expect(setequal(af$explosive_guard, c(TRUE, FALSE)), "both explosive guards")
expect(
  near(af$af_of_means[af$explosive_guard], 2 / 3) &&
    af$n_pairs[af$explosive_guard] == 3L &&
    af$n_pairs[!af$explosive_guard] == 4L,
  "guarded: replicate 4 dropped, AF 2/3; unguarded keeps it"
)
ref <- suppressMessages(attrib_metrics_tables(
  am,
  cmf,
  inp$adm1_pop,
  ctx$GROUP_SHORT_RUN,
  "neutral",
  a0,
  20,
  0.05,
  c(1, 2),
  cap = 2e5
))
expect(
  identical(
    select(filter(cl$af_scope, explosive_guard), -explosive_guard),
    ref$af_scope
  ),
  "guarded cell == attrib_metrics_tables() at the cap"
)
expect(
  nrow(cl$exposure) > 0 &&
    all(c("explosive_guard") %in% names(cl$exposure_adm1)),
  "Childs / IPCC exposure tables carried per guard"
)
expect(
  all(cl$explosive$explosive_guard),
  "explosive rows only from the guarded cell"
)

# --- [4] observed totals ----------------------------------------------------
# [4] observed totals per ADM1 / scope over each country's simulated months,
#     and AN_obs = AF x observed (interval from the AF interval).
cat("\n[4] observed totals and AN_obs\n")
obs <- tibble(
  gid = c("KEN_1_1", "KEN_1_1", "KEN_1_1", "KEN_1_2", "UGA_1_1", "UGA_1_1"),
  t_start = as.Date(c(
    "2014-08-01",
    "2014-10-01",
    "2016-12-01",
    "2015-01-01",
    "2015-03-01",
    "2015-04-01"
  )),
  obs_Y = c(1L, 1L, 1L, 0L, 1L, NA),
  cases_month = c(500, 60, 70, NA, 80, 90)
)
gids <- tibble(
  gid = c("KEN_1_1", "KEN_1_2", "UGA_1_1"),
  country = c("KEN", "KEN", "UGA")
)
last <- c(KEN = as.Date("2016-09-01"), UGA = as.Date("2015-03-01"))
oa <- event_observed_adm1(obs, gids, last, "2014-09-01")
ken <- filter(oa, adm1_unit == "KEN_1")
uga <- filter(oa, adm1_unit == "UGA_1")
expect(
  ken$obs_onsets == 1 && ken$obs_cases == 60,
  "KEN: before-launch and after-last months excluded"
)
expect(
  ken$n_cells == 2 * 25 && ken$n_cells_onset == 2 && ken$n_cells_cases == 1,
  "cells: 2 gids x 25 months; panel and case cells counted"
)
expect(
  uga$obs_onsets == 1 && uga$obs_cases == 80 && uga$n_cells == 7,
  "UGA: own last month (2015-03)"
)
univ <- tibble(
  adm1_gid = c("KEN_1", "UGA_1"),
  country = c("KEN", "UGA"),
  enso_group = c("enso_east", "enso_horn"),
  region = c("East", "Horn"),
  pop = c(1, 1),
  pop_missing = FALSE
)
os <- event_observed_scopes(oa, univ)
afr <- filter(os, level == "continent")
expect(
  afr$obs_onsets == 2 && afr$obs_cases == 140 && near(afr$cover_onsets, 3 / 57),
  "Africa sums its ADM1s"
)
afs <- tibble(
  level = c("continent", "continent"),
  scope = "Africa",
  enso_group = NA_character_,
  metric = c("cases", "onsets"),
  mean_f = c(200, 4),
  af_of_means = c(0.5, 0.25),
  af_lo = c(0.4, 0.1),
  af_hi = c(0.6, 0.4)
)
an <- event_an_observed(afs, os)
ac <- filter(an, metric == "cases")
ao <- filter(an, metric == "onsets")
expect(
  near(ac$an_obs, 70) &&
    near(ac$an_obs_lo, 56) &&
    near(ac$an_obs_hi, 84) &&
    near(ac$obs_over_f, 0.7),
  "cases: AN_obs = 0.5 x 140, interval 0.4 / 0.6 x 140, observed / factual mean 0.7"
)
expect(
  near(ac$cover, 3 / 57),
  "cases: coverage = the survival panel's (case rows exist only in episodes)"
)
expect(
  near(ao$an_obs, 0.5) && near(ao$cover, 3 / 57),
  "onsets: AN_obs = 0.25 x 2, coverage carried"
)

cat("\n[4b] event_observed_panel()\n")
surv_t <- tibble(
  gid = c("KEN_1_1", "KEN_1_1", "UGA_1_1"),
  t_start = as.Date(c("2015-01-01", "2015-02-01", "2015-01-01")),
  outbreak = c(1L, 0L, 0L)
)
mon_t <- tibble(
  gid = c("KEN_1_1", "KEN_1_1", "KEN_1_1", "KEN_1_1", "UGA_1_1"),
  month = as.Date(c(
    "2015-01-01",
    "2015-01-01",
    "2015-01-01",
    "2015-02-01",
    "2015-03-01"
  )),
  admin_level = c("ADM2", "ADM2", "ADM1", "ADM2", "ADM2"),
  cases_month = c(30, 20, 40, NA, 7),
  cases_gap = c(FALSE, FALSE, FALSE, TRUE, FALSE)
)
op <- event_observed_panel(surv_t, mon_t)
expect(
  op$cases_month[op$gid == "KEN_1_1" & op$t_start == as.Date("2015-01-01")] ==
    50,
  "footprint total = the largest per-level sum (ADM2 30 + 20 over ADM1 40)"
)
expect(
  is.na(op$cases_month[
    op$gid == "KEN_1_1" & op$t_start == as.Date("2015-02-01")
  ]) &&
    op$obs_Y[op$gid == "KEN_1_1" & op$t_start == as.Date("2015-02-01")] == 0L,
  "a cases_gap month stays NA; the onset panel row is kept"
)
expect(
  op$cases_month[op$gid == "UGA_1_1" & op$t_start == as.Date("2015-03-01")] ==
    7 &&
    is.na(op$obs_Y[op$gid == "UGA_1_1" & op$t_start == as.Date("2015-03-01")]),
  "a case month without a survival row is kept (obs_Y NA)"
)

# --- [5] analytic crossbasis AF --------------------------------------------
# [5] report_event_af_model's pure core: closed-form AN on a linear basis
#     (pre-launch exposure, clamp, lead), scope sums, coefficient lookup.
cat("\n[5] analytic AF (report_event_af_model core)\n")
source("analysis/utils/dlnm.R")
vm <- list(
  single_lag = 1L,
  range = c(-1, 1),
  cb_col_names = "enso_cb1",
  onebasis_attrs = list(fun = "lin", intercept = FALSE)
)
oni <- c(`2014-12-01` = 1.0, `2015-01-01` = 0.2, `2015-02-01` = 2.0)
rows <- tibble(
  adm1_gid = c("KEN_1", "KEN_1"),
  month = as.Date(c("2015-01-01", "2015-03-01")),
  y = c(10, 5)
)
cfm <- rbind(c(0, 0, 0), c(1, 1, 1)) # D = 2, columns from 2015-02
coef_of <- function(u) matrix(0.5, 2, 1)
r <- event_model_af_rows(rows, coef_of, vm, 1L, oni, cfm, as.Date("2015-02-01"))
expect(
  near(r$an["KEN_1", 1], 5 * (1 - exp(-1))) &&
    near(r$an["KEN_1", 2], 5 * (1 - exp(-0.5))) &&
    r$observed[["KEN_1"]] == 15,
  "lin basis: pre-launch exposure contributes 0; AN_d = y (1 - exp(-b (x_obs - x_cf,d)))"
)
rc <- event_model_af_rows(
  rows,
  coef_of,
  vm,
  1L,
  oni,
  cfm,
  as.Date("2015-02-01"),
  clamp_margin = 0
)
expect(
  near(rc$an["KEN_1", 1], 5 * (1 - exp(-0.5))) && near(rc$an["KEN_1", 2], 0),
  "clamp at the knots: observed 2.0 held at 1.0"
)
rl <- event_model_af_rows(
  mutate(rows, month = as.Date(c("2014-11-01", "2015-01-01"))),
  coef_of,
  vm,
  -1L,
  oni,
  cfm,
  as.Date("2015-02-01")
)
expect(
  near(rl$an["KEN_1", 1], 5 * (1 - exp(-1))) &&
    near(rl$an["KEN_1", 2], 5 * (1 - exp(-0.5))),
  "a lead (lag -1): 2015-01 cases read the 2015-02 ONI"
)
expect(
  inherits(
    try(
      event_model_af_rows(
        mutate(rows, month = as.Date("2016-01-01")),
        coef_of,
        vm,
        1L,
        oni,
        cfm,
        as.Date("2015-02-01")
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  "a missing observed ONI month is refused"
)
univ <- tibble(
  adm1_gid = c("KEN_1", "UGA_1"),
  country = c("KEN", "UGA"),
  enso_group = c("enso_east", "enso_horn"),
  region = c("East", "Horn")
)
an <- rbind(r$an["KEN_1", ], c(0, 0))
obs_v <- c(KEN_1 = 15, UGA_1 = 0)
ts <- event_model_af_scopes(univ, obs_v, an)
afr <- filter(ts, level == "continent")
horn <- filter(ts, level == "region", enso_group == "enso_horn")
expect(
  near(afr$an_mean, mean(r$an["KEN_1", ])) &&
    near(afr$af_mean, mean(r$an["KEN_1", ] / 15)) &&
    afr$observed == 15,
  "Africa: AN mean over draws, AF = AN / observed per draw"
)
expect(
  is.na(horn$af_mean) && horn$an_mean == 0,
  "a scope with nothing observed has no AF"
)
expect(sum(ts$level == "adm1") == 2L, "per-ADM1 rows kept")
lv <- c("KEN_KEN_1", "KEN_KEN_2")
fake_sc <- list(
  grp = list(region = NA, country = "country", adm1 = "country:ADM1"),
  re_mean = list(
    country = array(
      0,
      c(1, 1, 1),
      dimnames = list("KEN", "Estimate", "enso_cb1")
    ),
    `country:ADM1` = array(
      0,
      c(1, 1, 1),
      dimnames = list(lv[1], "Estimate", "enso_cb1")
    )
  ),
  b_pop_draws = matrix(1, 2, 1),
  re_draws_mat = function(tier, lbl) {
    matrix(switch(tier, region = 0, country = 10, adm1 = 100), 2, 1)
  },
  parse_level = function(lbl, tier) {
    pp <- strsplit(lbl, "_")[[1]]
    if (tier == "country") {
      list(region = NA, country = pp[1], country_lbl = pp[1])
    } else {
      list(
        region = NA,
        country = pp[1],
        country_lbl = pp[1],
        adm1_gid = paste(pp[2:3], collapse = "_")
      )
    }
  }
)
co <- event_model_coef_lookup(fake_sc)
expect(
  all(co("KEN_1") == 111) && all(co("KEN_2") == 11) && all(co("UGA_1") == 1),
  "coefficients: ADM1 RE when fitted, else country level, else population"
)

# [6] .event_attrib_stamp() (report_event_attrib v4): a probability-ratio
#     table keeps its outbreak event as pr_event under the past event's
#     name; tables without `event` are stamped unchanged.
cat("\n[6] past-event stamp keeps the probability-ratio event\n")
prt <- tibble(event = c("onset_any", "cases_ge"), pr = c(1.2, 1.5))
st <- .event_attrib_stamp(
  prt,
  event = "elnino_2015",
  factual_source = "simulation"
)
expect(
  identical(names(st)[1:3], c("event", "factual_source", "pr_event")) &&
    identical(st$pr_event, c("onset_any", "cases_ge")) &&
    all(st$event == "elnino_2015") &&
    !anyDuplicated(st[c("event", "factual_source", "pr_event")]),
  "pr tables: pr_event kept, event = the past event, one row per (event, pr_event)"
)
st2 <- .event_attrib_stamp(
  tibble(metric = "cases", af = 0.1),
  event = "e",
  sub_window = "w1",
  factual_source = "s"
)
expect(
  identical(
    names(st2),
    c("event", "sub_window", "factual_source", "metric", "af")
  ),
  "tables without event: stamped in front, unchanged otherwise"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
