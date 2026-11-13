# test_sim_censor_helper.R -- utils/sim_io.R::censor_sim_by_country() and the
# forecast/hindcast manifest readers.
#
#   [1] monthly: rows past each country's firm horizon dropped, or kept and
#       flagged with keep_beyond; unmapped ISO3 refused; country_col path.
#   [2] events: starts past the horizon dropped; open episodes flagged
#       right_censored with end_month untouched.
#   [3] manifest readers on a v1 manifest: units = groups, runs = the groups'
#       run fields, firm == unit horizon, country -> group.
#   [4] the same readers on a v2 manifest: units, per-country firm vs unit
#       horizons, counterfactual sets incl. the primary alias.
#   [5] combine_group_runs(): horizon_by_country = NULL is the scalar behaviour;
#       with a map, each country stops at its own horizon.
#   [6] results_core::load_sim_run_monthly(): optional per-country censor and
#       the keep_beyond flag column.
#   [7] manifest_country_horizons() prefers countries.<iso>.firm_horizon_end /
#       unit_horizon_end (hindcast: countries.<iso>.years.<y>.*), falls back to
#       the group / unit otherwise; manifest_has_country_horizons(); the event
#       manifest view keeps them (per-country season ends).
#
# Synthetic temp files only. Run: Rscript analysis/tests/test_sim_censor_helper.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(qs2) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "config", "paths", "provenance", "gid", "sim_io")) source(sprintf("analysis/utils/%s.R", m))
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

hz <- c(KEN = "2027-01-01", ETH = "2027-06-01")
months <- seq(as.Date("2026-11-01"), as.Date("2027-08-01"), by = "month")
ms <- tidyr::expand_grid(replicate_id = 1:2, gid = c("KEN_1_1", "ETH_2_1"), month = months) %>%
  mutate(Y = 1L, cases_month = 2)

cat("\n[1] monthly\n")
c1 <- censor_sim_by_country(ms, hz)
expect(max(c1$month[startsWith(c1$gid, "KEN")]) == as.Date("2027-01-01") &&
         max(c1$month[startsWith(c1$gid, "ETH")]) == as.Date("2027-06-01"), "each country cut at its own horizon")
expect(nrow(c1) == 2L * (3L + 8L), "row count = reps x (KEN 3 + ETH 8 months)")
c2 <- censor_sim_by_country(ms, hz, keep_beyond = TRUE)
expect(nrow(c2) == nrow(ms) && sum(c2$beyond_firm_horizon) == nrow(ms) - nrow(c1), "keep_beyond flags instead of dropping")
ms_chr <- ms %>% mutate(month = format(month))
expect(identical(nrow(censor_sim_by_country(ms_chr, hz)), nrow(c1)), "character month column handled")
expect(!is.null(err_of(censor_sim_by_country(ms, hz["KEN"]))), "ISO3 without a horizon refused")
cc <- ms %>% mutate(country = substr(gid, 1, 3)) %>% select(-gid)
expect(nrow(censor_sim_by_country(cc, hz, country_col = "country")) == nrow(c1), "country_col path")
expect(nrow(censor_sim_by_country(ms[0, ], hz)) == 0L, "empty frame passes through")

cat("\n[2] events\n")
et <- tibble(replicate_id = 1L, gid = c("KEN_1_1", "KEN_1_1", "ETH_2_1"),
             start_month = as.Date(c("2026-12-01", "2027-03-01", "2027-05-01")),
             end_month   = as.Date(c("2027-04-01", "2027-04-01", "2027-05-01")))
e1 <- censor_sim_by_country(et, hz, month_col = "start_month", kind = "events")
expect(nrow(e1) == 2L && identical(e1$right_censored, c(TRUE, FALSE)), "late start dropped; open episode right_censored")
expect(identical(e1$end_month, et$end_month[c(1, 3)]), "end_month untouched")

cat("\n[3] manifest readers, v1\n")
man1 <- list(primary_analog_set = "neutral", groups = list(
  enso_east = list(horizon_end = "2027-01-01", iso3 = list("KEN", "TZA"),
                   factual_config = "f_e.yml", factual_merged_dir = "md_f_e",
                   primary_analog_set = "neutral",
                   counterfactual_sets = list(neutral = list(config = "c_e_n.yml", merged_dir = "md_c_e_n")),
                   counterfactual_config = "c_e_n.yml", counterfactual_merged_dir = "md_c_e_n"),
  enso_horn = list(horizon_end = "2027-06-01", iso3 = list("ETH"),
                   factual_config = "f_h.yml", factual_merged_dir = "md_f_h",
                   counterfactual_config = "c_h.yml", counterfactual_merged_dir = "md_c_h")))
u1 <- manifest_units(man1)
expect(identical(names(u1), c("enso_east", "enso_horn")) && identical(u1$enso_east$members, c("KEN", "TZA")) &&
         identical(u1$enso_east$groups, "enso_east"), "v1 units = groups")
r1 <- manifest_runs(man1, "factual")
expect(identical(r1$enso_horn$merged_dir, "md_f_h"), "v1 factual runs")
rc <- manifest_runs(man1, "counterfactual")
expect(identical(rc$enso_east$merged_dir, "md_c_e_n") && identical(rc$enso_horn$merged_dir, "md_c_h"),
       "v1 counterfactual: set entry, else the flat primary alias")
expect(length(manifest_runs(man1, "counterfactual", "elnino")) == 0L, "absent set -> no run")
expect(identical(manifest_country_horizons(man1, "firm"), manifest_country_horizons(man1, "unit")) &&
         identical(manifest_country_horizons(man1)[["TZA"]], "2027-01-01"), "v1: firm == unit horizon")
expect(identical(unname(manifest_country_groups(man1)[c("KEN", "ETH")]), c("enso_east", "enso_horn")),
       "v1 country -> group")

cat("\n[4] manifest readers, v2\n")
man2 <- c(man1, list(schema = 2L, primary_analog_set = "neutral",
  units = list(eh = list(members = list("KEN", "TZA", "ETH"), groups = list("enso_east", "enso_horn"),
                         horizon_end = "2027-06-01", factual_config = "f_eh.yml",
                         factual_merged_dir = "md_f_eh", primary_analog_set = "neutral",
                         counterfactual_sets = list(neutral = list(config = "c_eh.yml", merged_dir = "md_c_eh")))),
  countries = list(KEN = list(group = "enso_east", unit = "eh", firm_horizon_end = "2027-01-01"),
                   TZA = list(group = "enso_east", unit = "eh", firm_horizon_end = "2027-01-01"),
                   ETH = list(group = "enso_horn", unit = "eh", firm_horizon_end = "2027-06-01"))))
expect(identical(names(manifest_units(man2)), "eh"), "v2 units")
expect(identical(manifest_runs(man2)$eh$merged_dir, "md_f_eh") &&
         identical(manifest_runs(man2, "counterfactual")$eh$merged_dir, "md_c_eh"), "v2 runs per unit")
expect(identical(manifest_country_horizons(man2, "firm")[["KEN"]], "2027-01-01") &&
         identical(manifest_country_horizons(man2, "unit")[["KEN"]], "2027-06-01"), "v2 firm vs unit horizon")

cat("\n[5] combine_group_runs()\n")
td <- tempfile("comb_"); dir.create(td)
md <- file.path(td, "u"); dir.create(md)
qs_save(ms, file.path(md, "monthly_summary.qs")); qs_save(et, file.path(md, "event_table.qs"))
a <- suppressMessages(combine_group_runs(list(u = md), file.path(td, "c1"), "2027-08-01"))
b <- suppressMessages(combine_group_runs(list(u = md), file.path(td, "c2"), "2027-08-01",
                                         horizon_by_country = hz))
expect(nrow(a$monthly_summary) == nrow(ms), "NULL map: scalar cutoff only")
expect(nrow(b$monthly_summary) == nrow(c1) && nrow(b$event_table) == 2L, "map: per-country censor")

cat("\n[6] load_sim_run_monthly()\n")
rc_env <- new.env(); sys.source("analysis/utils/results_core.R", envir = rc_env)
pd <- file.path(td, "sim", "run", "parts", "tu=eh"); dir.create(pd, recursive = TRUE)
qs_save(ms %>% mutate(month = format(month)), file.path(pd, "monthly_summary.qs"))
file.create(file.path(pd, "_SUCCESS"))
L0 <- rc_env$load_sim_run_monthly(file.path(td, "sim"), "run", as.Date("2026-11-01"), as.Date("2027-08-01"))
L1 <- rc_env$load_sim_run_monthly(file.path(td, "sim"), "run", as.Date("2026-11-01"), as.Date("2027-08-01"),
                                  horizon_by_country = hz)
L2 <- rc_env$load_sim_run_monthly(file.path(td, "sim"), "run", as.Date("2026-11-01"), as.Date("2027-08-01"),
                                  horizon_by_country = hz, keep_beyond = TRUE)
expect(nrow(L0$country_month) == 2L * 2L * length(months), "no map: every month")
expect(max(L1$country_month$month[L1$country_month$country == "KEN"]) == as.Date("2027-01-01"),
       "map: KEN stops at its firm horizon")
expect("beyond_firm_horizon" %in% names(L2$country_month) && nrow(L2$country_month) == nrow(L0$country_month),
       "keep_beyond: all months, flag column")
# max_gid_cases (per-gid explosive cap) comes after country.
expect(identical(names(L0$country_month), c("replicate_id", "month", "n_cases", "n_outbreaks", "country", "max_gid_cases", "run")),
       "no map: column layout unchanged")

cat("\n[7] per-country horizons in countries.<iso> (season ends, 2026-10-27)\n")
# event / forecast layout: a country's own firm horizon wins over its group's
man3 <- man2
man3$countries$TZA$firm_horizon_end <- "2026-12-01"
expect(identical(manifest_country_horizons(man3, "firm")[["TZA"]], "2026-12-01") &&
         identical(manifest_country_horizons(man3, "firm")[["KEN"]], "2027-01-01"),
       "countries.<iso>.firm_horizon_end wins; others unchanged")
expect(identical(manifest_country_horizons(man3, "unit")[["TZA"]], "2027-06-01"), "unit horizon from the unit")
man3$countries$TZA$unit_horizon_end <- "2027-03-01"
expect(identical(manifest_country_horizons(man3, "unit")[["TZA"]], "2027-03-01"),
       "countries.<iso>.unit_horizon_end wins over the unit's")
# a v2 manifest whose country entries carry no horizon (registry hindcast today)
man4 <- man2; man4$countries <- lapply(man4$countries, function(x) x[c("group", "unit")])
expect(identical(manifest_country_horizons(man4, "firm"), manifest_country_horizons(man2, "firm")),
       "country entries without horizons: the group's (unchanged)")
expect(!manifest_has_country_horizons(man4) && manifest_has_country_horizons(man2),
       "manifest_has_country_horizons()")
# hindcast layout: groups.<g>.years.<y>, countries.<iso>.years.<y>
hman <- list(groups = list(
  enso_east = list(iso3 = list("KEN", "TZA"), years = list(`2015` = list(horizon_end = "2016-03-01"),
                                                           `2016` = list(horizon_end = "2017-03-01")))),
  countries = list(KEN = list(group = "enso_east", unit = "enso_east",
                              years = list(`2015` = list(firm_horizon_end = "2015-12-01",
                                                         unit_horizon_end = "2016-03-01"))),
                   TZA = list(group = "enso_east", unit = "enso_east")))
hz15 <- manifest_country_horizons(hman, "firm", year = 2015)
expect(identical(hz15[["KEN"]], "2015-12-01") && identical(hz15[["TZA"]], "2016-03-01"),
       "hindcast: countries.<iso>.years.<y>, else the group's year horizon")
expect(identical(manifest_country_horizons(hman, "firm", year = 2016)[["KEN"]], "2017-03-01"),
       "a year the country entry lacks falls back to the group")
expect(identical(manifest_country_horizons(hman, "unit", year = 2015)[["KEN"]], "2016-03-01") &&
         identical(manifest_country_horizons(hman, "unit", year = 2015)[["TZA"]], "2016-03-01"),
       "v1 layout: unit horizon = the country's unit_horizon_end, else the group's")
expect(manifest_has_country_horizons(hman, 2015) && !manifest_has_country_horizons(hman, 2016),
       "manifest_has_country_horizons() per year")
# an event manifest view keeps the per-country horizons
em <- list(kind = "event_attribution", grouping = "g", analog_set = "neutral",
           events = list(e1 = list(launch = "2015-01-01", africa_max_horizon = "2016-06-01",
                                   groups = man2$groups, units = man2$units, countries = man3$countries)))
expect(identical(manifest_country_horizons(event_manifest_view(em, "e1"), "firm")[["TZA"]], "2026-12-01"),
       "event_manifest_view() carries the countries' own horizons")

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
