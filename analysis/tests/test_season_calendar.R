# test_season_calendar.R -- utils/season_calendar.R: per-country cholera-season
# ends from the gid calendar (min_seas) and the window rules of the hindcast and
# past-event generators under --season_end by_country_from_panel.
#
# HARD, quick-safe: pure helpers on synthetic tables, no data, no fit.
#   [1] season_end_month(): min_seas - 1, wrapping 1 -> 12; refuses 0 / 13 / NA
#   [2] season_end_by_country(): mode over gids (one row per gid however many
#       panel rows), ties -> the later end, modal share / n_gid / n_distinct,
#       wrap-around values, refusals (a gid with two calendars, missing columns)
#   [3] season_end_span(): circular spread (Nov-Apr = 5, not 7)
#   [4] first_season_end(): on / strictly after / at least k months after
#   [5] hindcast rule: first end >= launch + 6 months, + max(min_lag, 0)
#   [6] event rule: first end on or after the peak + max positive lag; equals
#       make_event_attribution_configs.R::event_group_window()
#   [7] season_end_spread_warnings(): a hemisphere-mixing group is named
#   [8] gid_min_seas_source(): canonical value (mode, ties -> smallest month,
#       as footprint.R::canonicalSeasonCalendar()) and its best source
#   [9] panel_season_calendar(): one row per gid, digest attribute; the CSV path
#
# Run: Rscript analysis/tests/test_season_calendar.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/season_calendar.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
d <- as.Date

cat("# season calendar\n")

cat("\n[1] season_end_month()\n")
expect(identical(season_end_month(1:12), c(12L, 1:11)), "min_seas 1..12 -> ends 12, 1..11")
expect(identical(season_end_month(1L), 12L), "min_seas 1 (January start) ends in December")
expect(!is.null(err_of(season_end_month(0L))) && !is.null(err_of(season_end_month(13L))) &&
         !is.null(err_of(season_end_month(NA))), "0, 13 and NA are refused")

cat("\n[2] season_end_by_country()\n")
pan <- data.frame(
  gid      = c("A_1", "A_1", "A_2", "A_3", "B_1", "B_2", "C_1", "C_2", "C_3", "C_4"),
  country  = c("AAA", "AAA", "AAA", "AAA", "BBB", "BBB", "CCC", "CCC", "CCC", "CCC"),
  min_seas = c(5L,    5L,    5L,    9L,    5L,    9L,    11L,   11L,   1L,    1L))
st <- season_end_by_country(pan)
expect(identical(st$iso3, c("AAA", "BBB", "CCC")), "one row per country, sorted")
a <- st[st$iso3 == "AAA", ]
expect(a$min_seas_mode == 5L && a$season_end_month == 4L && a$n_gid == 3L &&
         isTRUE(all.equal(a$modal_share, 2 / 3)) && a$n_distinct == 2L,
       "AAA: mode over GIDS (A_1's two rows count once): min_seas 5 -> end April, share 2/3")
b <- st[st$iso3 == "BBB", ]
expect(b$min_seas_mode == 9L && b$season_end_month == 8L && b$modal_share == 0.5,
       "BBB tie min_seas 5 (end Apr) vs 9 (end Aug): the later end, August")
c_ <- st[st$iso3 == "CCC", ]
expect(c_$season_end_month == 12L && c_$min_seas_mode == 1L,
       "CCC tie 11 (end 10) vs 1 (end 12): end 12 (wrap-around min_seas 1) wins")
expect(identical(names(st), c("iso3", "season_end_month", "min_seas_mode", "modal_share", "n_gid",
                              "n_distinct")), "column layout")
bad <- rbind(pan, data.frame(gid = "A_2", country = "AAA", min_seas = 6L))
expect(!is.null(err_of(season_end_by_country(bad))), "a gid with two min_seas is refused")
expect(!is.null(err_of(season_end_by_country(pan[, c("gid", "min_seas")]))), "missing country column refused")
expect(!is.null(err_of(season_end_by_country(pan, aggregate = "latest"))), "only aggregate = mode")
cc <- pan; names(cc)[2] <- "iso"
expect(identical(season_end_by_country(cc, country_col = "iso")$season_end_month, st$season_end_month),
       "country_col")

cat("\n[3] season_end_span()\n")
expect(season_end_span(7L) == 0L, "one month: 0")
expect(season_end_span(c(4L, 11L)) == 5L, "{Apr, Nov}: the arc Nov-Apr, 5 months")
expect(season_end_span(c(7L, 9L, 11L)) == 4L, "{Jul, Sep, Nov}: 4")
expect(season_end_span(c(12L, 1L, 2L)) == 2L, "{Dec, Jan, Feb}: 2 across the year end")

cat("\n[4] first_season_end()\n")
expect(first_season_end(d("2015-05-01"), 5L, 0L) == d("2015-05-01"), "on: the month itself")
expect(first_season_end(d("2015-05-01"), 5L, 1L) == d("2016-05-01"), "strictly after: next year")
expect(first_season_end(d("2015-05-01"), 4L, 0L) == d("2016-04-01"), "earlier month: next year")
expect(first_season_end(d("2015-11-01"), 2L, 0L) == d("2016-02-01"), "wrap across the year end")

cat("\n[5] hindcast_window_end(): first end >= launch + 6 months\n")
L <- d("2015-05-01")
w <- hindcast_window_end(L, 11L, 0L)
expect(w$season_end == d("2015-11-01") && w$horizon_end == d("2015-11-01"), "end Nov: exactly 6 months, counts")
w <- hindcast_window_end(L, 10L, 0L)
expect(w$season_end == d("2016-10-01"), "end Oct: 5 months is too short -> next year (17 months)")
w <- hindcast_window_end(L, 7L, 0L)
expect(w$season_end == d("2016-07-01"), "end Jul: the July 2015 stub (2 months) is skipped")
w <- hindcast_window_end(L, 12L, 3L)
expect(w$season_end == d("2015-12-01") && w$horizon_end == d("2016-03-01"), "+ max(min_lag, 0) = 3")
w <- hindcast_window_end(L, 12L, -6L)
expect(w$horizon_end == d("2015-12-01"), "a lead (negative min lag) does not shorten the window")
w <- hindcast_window_end(L, 12L, 3L, add_best_lag = FALSE)
expect(w$horizon_end == d("2015-12-01"), "add_best_lag FALSE: no extension")
w <- hindcast_window_end(L, 6L, 0L, min_window_months = 1L)
expect(w$season_end == d("2015-06-01"), "min_window_months 1 = the old 'strictly after' rule")

cat("\n[6] event_window_end() and event_group_window()\n")
w <- event_window_end(d("2015-12-01"), 7L, c(0L, -1L, 2L))
expect(w$season_end == d("2016-07-01") && w$horizon_end == d("2016-09-01") && w$max_lag == 2L,
       "peak Dec 2015, end Jul: 2016-07 + 2 = 2016-09")
w <- event_window_end(d("2020-11-01"), 11L, c(-6L, -4L, -8L))
expect(w$season_end == d("2020-11-01") && w$horizon_end == d("2020-11-01") && w$max_lag == 0L,
       "end on the peak month counts; leads add nothing")
w <- event_window_end(d("2015-12-01"), 12L, 10L)
expect(w$horizon_end == d("2016-10-01"), "end December on a December peak: 2015-12 + 10 = 2016-10")
# The group rule as make_event_attribution_configs.R::event_group_window()
# computed it before delegating to event_window_end(), verbatim.
suppressPackageStartupMessages(library(lubridate))
old_rule <- function(peak, sem, lags) {
  peak <- as.Date(peak)
  end0 <- as.Date(sprintf("%04d-%02d-01", as.integer(format(peak, "%Y")), sem))
  if (end0 < peak) end0 <- end0 %m+% years(1L)
  max_lag <- max(c(as.integer(lags), 0L))
  list(season_end = end0, horizon_end = end0 %m+% months(max_lag), max_lag = max_lag)
}
grid <- expand.grid(peak = c("2015-12-01", "2020-11-01", "2023-01-01", "2016-02-01"), sem = 1:12,
                    stringsAsFactors = FALSE)
same <- vapply(seq_len(nrow(grid)), function(i)
  identical(old_rule(grid$peak[i], grid$sem[i], c(3L, -2L)),
            event_window_end(d(grid$peak[i]), grid$sem[i], c(3L, -2L))), logical(1))
expect(all(same), "event_window_end() == the former event_group_window() rule (4 peaks x 12 end months)")

cat("\n[7] season_end_spread_warnings()\n")
se <- data.frame(iso3 = c("CIV", "AGO", "ZMB", "KEN", "TZA"), season_end_month = c(11L, 4L, 5L, 7L, 8L))
wn <- season_end_spread_warnings(list(enso_south = c("CIV", "AGO", "ZMB"), enso_east = c("KEN", "TZA")), se)
expect(length(wn) == 1L && grepl("enso_south", wn) && grepl("6 months", wn),
       "the hemisphere-mixing group (Nov vs Apr/May: 6 months) is named, the other not")
expect(length(season_end_spread_warnings(list(g = "KEN"), se)) == 0L, "one country: nothing")

cat("\n[8] gid_min_seas_source()\n")
mem <- data.frame(gid = c("G1", "G1", "G1", "G2", "G2", "G3"),
                  location_period_id = c("1", "2", "3", "4", "5", "6"),
                  min_seas = c(3L, 3L, 1L, 1L, 4L, 1L))
src <- data.frame(location_period_id = as.character(1:6),
                  min_seas_source = c("default_january", "estimate_adm1", "estimate_gid",
                                      "default_january", "estimate_gid", "default_january"))
gs <- gid_min_seas_source(mem, src)
expect(identical(gs$gid, c("G1", "G2", "G3")), "one row per gid")
expect(gs$min_seas[1] == 3L && gs$min_seas_source[1] == "estimate_adm1" && gs$n_rows_default[1] == 1L,
       "G1: mode 3, best source among the rows carrying it (adm1 beats default)")
expect(gs$min_seas[2] == 1L && gs$min_seas_source[2] == "default_january",
       "G2: tie 1 vs 4 -> smallest month 1 (canonicalSeasonCalendar's rule), defaulted")
expect(gs$min_seas_source[3] == "default_january", "G3: January default")
expect(!is.null(err_of(gid_min_seas_source(mem, src[-1, ]))), "an LP without a source is refused")

cat("\n[9] panel_season_calendar() / season_ends_csv_path()\n")
pc <- panel_season_calendar(pan, "x.rds")
expect(nrow(pc) == 9L && identical(attr(pc, "path"), "x.rds"), "one row per gid, path kept")
expect(is.character(attr(pc, "digest")), "digest attribute")
pc2 <- panel_season_calendar(pan[rev(seq_len(nrow(pan))), ], "x.rds")
expect(identical(attr(pc, "digest"), attr(pc2, "digest")), "digest does not depend on row order")
expect(identical(season_end_source(pc)$panel, "x.rds"), "season_end_source() records the panel")
expect(identical(season_ends_csv_path("analysis/configs/hindcast_horizons_prod_run.yml"),
                 "analysis/configs/season_ends_hindcast_horizons_prod_run.csv"),
       "season_ends_<manifest stem>.csv next to the manifest")

cat("\n[10] season windows: onset-rate troughs, trough choice, window rule (2026-11-01)\n")
expect(isTRUE(all.equal(.sc_circ_smooth(c(3, rep(0, 10), 3), 3L)[c(1, 2, 12)], c(2, 1, 2))) &&
         identical(.sc_circ_smooth(1:12, 1L), as.numeric(1:12)),
       ".sc_circ_smooth(): circular centred mean; k = 1 is the identity")
expect(identical(.sc_circ_dist(c(1L, 12L, 7L), c(12L, 2L, 1L)), c(1L, 2L, 6L)), ".sc_circ_dist() wraps")
# AAA: onsets every year in Jun-Sep, none Dec-Feb (trough Jan); BBB: 2 onsets
# in March only (trough elsewhere, but too few); CCC: at risk, no onset.
sv <- do.call(rbind, lapply(2010:2019, function(y) {
  ms <- seq(as.Date(sprintf("%d-01-01", y)), by = "month", length.out = 12L)
  rbind(data.frame(country = "AAA", t_start = ms,
                   outbreak = as.integer(as.integer(format(ms, "%m")) %in% 6:9)),
        data.frame(country = "BBB", t_start = ms,
                   outbreak = as.integer(y %in% 2010:2011 & format(ms, "%m") == "03")),
        data.frame(country = "CCC", t_start = ms, outbreak = 0L))
}))
ot <- onset_rate_troughs(sv, bootstrap = 50L)
expect(identical(ot$iso3, c("AAA", "BBB", "CCC")) && identical(ot$n_onsets, c(40L, 2L, 0L)) &&
         identical(ot$n_atrisk, c(120L, 120L, 120L)),
       "onset_rate_troughs(): one row per country, onset and at-risk counts")
expect(identical(ot$onset_trough[1], 1L) && identical(ot$stability[1], 1),
       "AAA: trough Jan (smoothed Dec-Feb zero, first minimum), stable under every resample")
expect(is.na(ot$onset_trough[3]) && is.na(ot$stability[3]), "no onsets: trough and stability NA")
expect(identical(unname(attr(ot, "onsets_by_month")["BBB", ]), c(0, 0, 2, rep(0, 9))),
       "onsets_by_month attribute (raw counts)")
expect(identical(onset_rate_troughs(sv, years = 2012:2019, bootstrap = 0L)$n_onsets, c(32L, 0L, 0L)),
       "years = filters the calendar years read")
expect(identical(unlist(onset_rate_troughs(sv, bootstrap = 50L, seed = 7L)[1, -1]),
                 unlist(onset_rate_troughs(sv[sv$country == "AAA", ], bootstrap = 50L, seed = 7L)[1, -1])),
       "a country's stability does not depend on the other countries (seeded per country)")
cal10 <- data.frame(iso3 = c("AAA", "BBB", "CCC"), min_seas_mode = c(4L, 8L, 11L),
                    modal_share = 1, stringsAsFactors = FALSE)
st <- season_trough_by_country(cal10, ot, "onset_rate", min_onsets = 30L, min_stability = 0.8)
expect(identical(st$trough_month, c(1L, 8L, 11L)) &&
         identical(st$trough_source, c("onset_rate", "min_seas", "min_seas")),
       "onset_rate: AAA uses its onset trough; BBB (2 onsets) and CCC (none) fall back to min_seas")
expect(identical(season_trough_by_country(cal10, ot, "onset_rate", min_onsets = 30L,
                                          min_stability = 1.01)$trough_source[1], "min_seas"),
       "an unstable trough (stability below the threshold) falls back to min_seas")
expect(identical(season_trough_by_country(cal10, NULL, "min_seas")$trough_month, c(4L, 8L, 11L)),
       "rule min_seas: always the modal min_seas")
expect(!is.null(err_of(season_trough_by_country(cal10, NULL, "onset_rate"))),
       "rule onset_rate without onsets is refused")
w <- season_window_of(7L, 2015L, 5L, min_lag = 3L)
expect(identical(w$launch, as.Date("2015-08-01")) && identical(w$season_end, as.Date("2016-07-01")) &&
         identical(w$horizon_end, as.Date("2016-10-01")) && !w$floored && identical(w$months_cut, 0L),
       "trough Jul: season Aug 2015 .. Jul 2016, firm end + min lag 3")
w <- season_window_of(1L, 2015L, 5L, min_lag = -6L)
expect(identical(w$start_month, 2L) && identical(w$launch, as.Date("2015-05-01")) && w$floored &&
         identical(w$months_cut, 3L) && identical(w$season_end, as.Date("2016-01-01")) &&
         identical(w$horizon_end, as.Date("2016-01-01")),
       "trough Jan: start Feb floored to May (3 months cut), ends Jan 2016; a lead adds nothing")
w <- season_window_of(12L, 2015L, 5L)
expect(identical(w$launch, as.Date("2015-05-01")) && identical(w$months_cut, 4L) &&
         identical(w$season_end, as.Date("2015-12-01")),
       "trough Dec: start Jan floored to May (4 cut), ends Dec of the same year")
w <- season_window_of(11L, 2015L, 5L)
expect(identical(w$launch, as.Date("2015-12-01")) && identical(w$season_end, as.Date("2016-11-01")),
       "trough Nov: launch Dec, a full 12-month season")
w <- season_window_of(4L, 2015L, 5L, min_lag = 8L, add_best_lag = FALSE)
expect(identical(w$launch, as.Date("2015-05-01")) && !w$floored &&
         identical(w$season_end, as.Date("2016-04-01")) && identical(w$horizon_end, w$season_end),
       "trough Apr: start May (not floored); add_best_lag FALSE ends at the season end")
expect(isTRUE(all.equal(season_cut_share(c(1, 2, 3, 4, rep(0, 8)), 2L, 5L), 0.9)) &&
         identical(season_cut_share(c(1, rep(0, 11)), 6L, 5L), 0) &&
         is.na(season_cut_share(rep(0, 12), 2L, 5L)),
       "season_cut_share(): share of onsets in start .. floor - 1; 0 when not floored; NA without onsets")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
