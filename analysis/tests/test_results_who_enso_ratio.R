# test_results_who_enso_ratio.R -- supplementary fig. WR (09 step
# report_who_enso_ratio; helpers in utils/who_annual.R, utils/enso.R and
# utils/results_derive.R; builder in utils/figures_supp.R).
#
#   [1] buildWhoAnnualCases(): an absent country-year is NA (never zero), the
#       annual report's duplicate rows take the maximum, report + dashboard
#       take the mean, the weekly series is summed by epiweek year, a target-
#       year series stopping before the cut-off is censored (NA, kept as
#       cases_censored), duplicate weeks stop.
#   [2] whoWindowTotals(): a year with a partial or missing window month has
#       no total.
#   [3] ensoYearStates(): annual mean, the lag-shifted window, DJF starting
#       and ending the year, threshold phases, NA past the record.
#   [4] who_enso_observed_ratio(): ratio of means with the pseudo-count, only
#       neutral years where every country has a value, no interval from one
#       year, a deterministic bootstrap interval bracketing the point.
#   [5] who_enso_smoothing_totals() zero-fills; who_enso_cf_parent_means()
#       averages the children of each smoothing replicate and stops on unequal
#       children, a missing window month or an unmapped replicate.
#   [6] who_enso_model_ratio(): ratio of means and the per-replicate spread by
#       hand; P(replicate ratio > 1).
#   [7] who_enso_ratio_tables(): the group row pools the countries with an
#       observed ratio; population standardisation; every analogue set.
#   [8] who_enso_ratio_frame() / build_figS_who_enso_ratio() draw.
#
# Synthetic: no artefact, no data/. HARD, quick-safe.
# Run: Rscript analysis/tests/test_results_who_enso_ratio.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
.md <- new.env()
for (m in c("base", "enso", "who_annual")) sys.source(sprintf("analysis/utils/%s.R", m), envir = .md)
.rd <- parse("analysis/utils/results_derive.R", keep.source = FALSE)
for (e in .rd) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                   as.character(e[[2]]) %in% c(".with_seed", "who_enso_smoothing_totals", "who_enso_cf_parent_means",
                                               "who_enso_model_ratio", "who_enso_observed_ratio",
                                               "who_enso_standardise", "who_enso_ratio_tables")) eval(e, .md)
sys.source("analysis/utils/figures_supp.R", envir = .md)
attach(.md, name = "who_enso_ratio_helpers", warn.conflicts = FALSE)
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

cat("# [1] buildWhoAnnualCases\n")
report <- tibble(iso3 = c("AAA", "AAA", "AAA", "BBB", "AAA"),
                 TL = as.Date(c("2020-01-01", "2021-01-01", "2021-01-01", "2021-01-01", "2020-01-01")),
                 sCh = c(10, 5, 3, 8, 12),
                 spatial_scale = factor(c("country", "country", "country", "country", "admin1")))
wk <- function(iso, dates, cases, epiyr) tibble(iso_3_code = iso, date_wk = as.Date(dates), cases = cases, epiyr = epiyr)
weekly <- bind_rows(
  wk("AAA", c("2021-01-04", "2021-01-11"), c(3L, 4L), 2021L),
  wk("BBB", c("2025-12-29", "2026-01-05", "2026-05-04", "2026-09-14"), c(2L, 1L, 6L, 3L), c(2025L, 2026L, 2026L, 2026L)),
  wk("CCC", c("2026-01-05", "2026-01-12"), c(9L, 9L), 2026L))
w <- buildWhoAnnualCases(report, weekly, c("AAA", "BBB", "CCC"), 2020:2026, 2026L, censor_weeks = 8L)
a <- w$annual
cell <- function(c, y, col = "cases") a[[col]][a$country == c & a$year == y]
expect(near(cell("AAA", 2020), 10), "AAA 2020: the country rows' maximum (the admin1 row is not counted)")
expect(near(cell("AAA", 2021), 6), "AAA 2021: mean of the report (max 5) and the dashboard (3 + 4)")
expect(is.na(cell("AAA", 2022)) && cell("AAA", 2022, "source") == "not_in_report", "AAA 2022: absent is NA")
expect(is.na(cell("AAA", 2026)) && cell("AAA", 2026, "source") == "not_on_dashboard",
       "AAA 2026: off the dashboard is NA, not zero")
expect(near(cell("BBB", 2025), 2) && near(cell("BBB", 2026), 10),
       "BBB: the 2025-12-29 week counts in epiweek year 2025")
expect(identical(w$cutoff, as.Date("2026-09-14")), "cut-off = the last 2026 week")
expect(identical(w$censored, "CCC") && is.na(cell("CCC", 2026)) && near(cell("CCC", 2026, "cases_censored"), 18),
       "CCC stops in January: censored, NA, kept as cases_censored")
expect(errs(buildWhoAnnualCases(report, bind_rows(weekly, weekly[1, ]), "AAA", 2020:2026, 2026L)),
       "a duplicate week stops")

cat("# [2] whoWindowTotals\n")
mon <- tibble(iso3 = "AAA", t_start = as.Date(c("2024-05-01", "2024-06-01", "2025-05-01", "2025-06-01", "2023-05-01")),
              cases = c(1, 2, 3, 4, 5), partial = c(FALSE, FALSE, FALSE, TRUE, FALSE))
wt <- whoWindowTotals(mon, "AAA", 5:6)
expect(identical(wt$year, 2024L) && near(wt$cases, 3), "only the complete window (2024) has a total")

cat("# [3] ensoYearStates\n")
oni <- tibble(year = rep(2019:2021, each = 12), month = rep(1:12, 3),
              enso_index = c(rep(-1, 12), rep(0, 12), rep(1, 12)))
st <- ensoYearStates(oni, 2019:2022, lag_months = -4L)
expect(near(st$oni_annual_mean[1:3], c(-1, 0, 1)) && is.na(st$oni_annual_mean[4]), "annual means, NA past the record")
expect(near(st$oni_lag_occ[2], (-4 + 0 * 8) / 12) && st$n_months_lag_occ[2] == 12L,
       "lag -4: September of the previous year to August")
expect(near(st$oni_djf_start[2], 0) && near(st$oni_djf_end[2], 1) && is.na(st$oni_djf_end[3]),
       "DJF starting / ending the year (January values)")
expect(identical(st$state_annual_mean[1:3], c("nina", "neutral", "elnino")) &&
         identical(st$state_lag_occ[2], "neutral") && is.na(st$state_djf_end[3]),
       "threshold phases; NA without data")

cat("# [4] who_enso_observed_ratio\n")
obs <- tibble(country = c("A", "A", "A", "B", "B", "A", "B"), year = c(2018L, 2019L, 2020L, 2018L, 2020L, 2026L, 2026L),
              cases_std = c(10, 30, 50, 0, 20, 99, 1))
r1 <- who_enso_observed_ratio(obs, c("A", "B"), c(2018L, 2019L), 2026L, pseudo = 1, n_boot = 200L, seed = 7L)
expect(r1$neutral_years == "2018" && near(r1$obs_ratio, (100 + 1) / (10 + 1)) && is.na(r1$obs_lo),
       "pooled: 2019 dropped (B has no value); one neutral year, no interval")
r2 <- who_enso_observed_ratio(obs, "A", c(2018L, 2019L, 2020L), 2026L, pseudo = 1, n_boot = 500L, seed = 7L)
r2b <- who_enso_observed_ratio(obs, "A", c(2018L, 2019L, 2020L), 2026L, pseudo = 1, n_boot = 500L, seed = 7L)
expect(near(r2$obs_ratio, 100 / 31) && r2$obs_lo <= r2$obs_ratio && r2$obs_hi >= r2$obs_ratio &&
         identical(r2$obs_lo, r2b$obs_lo),
       "three years: ratio of means, a seeded interval bracketing it")
expect(near(r2$obs_single_year_min, 100 / 51) && near(r2$obs_single_year_max, 100 / 11), "single-year range")
expect(is.na(who_enso_observed_ratio(obs, "A", 2021L, 2026L)$obs_ratio), "no neutral year with data: NA")

cat("# [5] smoothing totals and counterfactual parent means\n")
months <- as.Date(c("2026-05-01", "2026-06-01"))
# No `country` column: the ISO3 prefix of the gid (as AAA_1_2).
ev <- tibble(rep_id = c(1L, 1L, 2L, 2L), gid = c("AAA_1", "AAA_2", "AAA_1", "BBB_1"),
             t_start = as.Date(c("2026-05-01", "2026-06-01", "2026-07-01", "2026-05-01")), cases_month = c(4, 6, 50, 7))
S <- who_enso_smoothing_totals(ev, 1:3, c("AAA", "BBB"), months)
expect(nrow(S) == 6L && near(S$cases[S$rep_id == 1 & S$country == "AAA"], 10) &&
         near(S$cases[S$rep_id == 2 & S$country == "AAA"], 0) && near(S$cases[S$rep_id == 2 & S$country == "BBB"], 7) &&
         near(S$cases[S$rep_id == 3 & S$country == "BBB"], 0),
       "window months only, zero-filled over replicates x countries")
cm <- tidyr::expand_grid(replicate_id = 1:6, month = months, country = c("A", "B")) %>%
  mutate(n_cases = ifelse(country == "A", replicate_id, 0))
hm <- tibble(replicate_id = 1:6, smoothing_rep_id = rep(1:3, each = 2))
C <- who_enso_cf_parent_means(cm, hm, c("A", "B"), months)
expect(near(C$cases[C$smoothing_rep_id == 2 & C$country == "A"], (3 * 2 + 4 * 2) / 2) && all(C$n_children == 2L),
       "parent mean over its two children, each summed over the window")
expect(errs(who_enso_cf_parent_means(cm, tibble(replicate_id = 1:6, smoothing_rep_id = c(1, 1, 1, 2, 2, 3)),
                                     c("A", "B"), months)), "unequal children stop")
expect(errs(who_enso_cf_parent_means(cm[!(cm$replicate_id == 1 & cm$month == months[2]), ], hm, c("A", "B"), months)),
       "a replicate missing a window month stops")
expect(errs(who_enso_cf_parent_means(cm, hm[-1, ], c("A", "B"), months)), "an unmapped replicate stops")
expect(errs(who_enso_cf_parent_means(cm, NULL, c("A", "B"), months)), "no handoff map stops")

cat("# [6] who_enso_model_ratio\n")
S6 <- tibble(rep_id = 1:2, country = "A", cases = c(10, 30))
C6 <- tibble(smoothing_rep_id = 1:2, country = "A", cases = c(5, 10), n_children = 2L)
m6 <- who_enso_model_ratio(S6, C6, "A", pseudo = 1, interval = 0.5, n_boot = 100L, seed = 3L)
rr <- c(11 / 6, 31 / 11)
expect(near(m6$model_ratio, 21 / 8.5) && near(m6$model_rep_median, mean(rr)) &&
         near(m6$model_lo, unname(quantile(rr, 0.25))) && near(m6$model_hi, unname(quantile(rr, 0.75))) &&
         near(m6$prob_rep_ratio_gt1, 1) && m6$n_smoothing_reps == 2L,
       "ratio of means and the per-smoothing-replicate spread")
expect(errs(who_enso_model_ratio(S6, C6[1, ], "A")), "a smoothing replicate without children stops")

cat("# [7] who_enso_ratio_tables\n")
annual <- tibble(country = c("A", "A", "A", "B", "B", "B"), year = c(2018L, 2019L, 2026L, 2018L, 2019L, 2026L),
                 cases = c(10, 20, 60, 5, NA, NA))
pop <- tidyr::expand_grid(country = c("A", "B"), year = 2018:2026) %>% mutate(pop = ifelse(year == 2026L, 200, 100))
states <- list(g1 = tibble(year = c(2018L, 2019L, 2026L), state_annual_mean = c("neutral", "neutral", "elnino")))
Sx <- tidyr::expand_grid(rep_id = 1:2, country = c("A", "B")) %>% mutate(cases = c(10, 1, 30, 1))
Cx <- tidyr::expand_grid(smoothing_rep_id = 1:2, country = c("A", "B")) %>% mutate(cases = c(5, 1, 10, 1), n_children = 2L)
cfg <- list(groups = "g1", year_definitions = "annual_mean", pseudo_count = 1, interval = 0.95, n_boot = 50L, seed = 1L)
tb <- who_enso_ratio_tables(annual, tibble(country = character(), year = integer(), cases = numeric()), pop, states,
                            c(A = "g1", B = "g1"), c("A", "B"), Sx, list(neutral = Cx, all = Cx), 2026L, cfg)
cp <- tb$comparison %>% filter(estimator == "ytd_annual", analog_set == "neutral")
expect(near(cp$obs_ratio[cp$scope == "A"], (60 + 1) / (2 * 15 + 1)),
       "population-standardised (2026 pop / year pop = 2) ratio of means")
expect(is.na(cp$obs_ratio[cp$scope == "B"]) && identical(cp$countries[cp$scope == "g1"], "A"),
       "B has no 2026 value; the group row pools A only")
expect(near(cp$model_ratio[cp$scope == "g1"], cp$model_ratio[cp$scope == "A"]), "group model over the same countries")
expect(setequal(unique(tb$comparison$analog_set), c("neutral", "all")) &&
         nrow(tb$model) == 2L * 3L && identical(unique(cp$state_target), "elnino"),
       "every analogue set; model rows per country + group; the target year's state")

cat("# [8] who_enso_ratio_frame / build_figS_who_enso_ratio\n")
fr <- who_enso_ratio_frame(tb$comparison, "g1", "annual_mean", "neutral")
expect(identical(sort(unique(fr$kind)), c("country", "group")) && all(is.finite(fr$obs_ratio)),
       "frame: observed points only, country and group rows")
expect(errs(who_enso_ratio_frame(tb$comparison, "g1", "djf_end", "neutral")), "an underived definition stops")
fr2 <- bind_rows(fr, mutate(fr, estimator = "same_window"))
p <- tryCatch(build_figS_who_enso_ratio(fr2, "West", "#66a61e", ggplot2::theme_bw(base_size = 7),
                                        c(same_window = "a", ytd_annual = "b"), model_lab = "model"),
              error = function(e) e)
expect(inherits(p, "patchwork"), "a two-panel patchwork")
ok <- tryCatch({ grDevices::pdf(NULL); print(p); grDevices::dev.off(); TRUE }, error = function(e) conditionMessage(e))
expect(isTRUE(ok), paste("the figure draws:", ok))

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
