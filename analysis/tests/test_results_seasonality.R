# test_results_seasonality.R -- supplementary fig. SE (09 step seasonality,
# helpers in utils/results_derive.R; builder in utils/figures_supp.R).
#
#   [1] seasonality_track_rates() on synthetic panels: onsets per 1,000
#       unit-months at risk and ends per 100 episode-months against a direct
#       computation; the ADM1-ADM2 rule drops an ADM1 parent with ADM2 rows;
#       cases_gap months are left out of the case median and quartiles;
#       countries outside the groups are dropped; every group x track x month
#       is present (n 0, rate NA when empty); sparse = n < min_n; the exact
#       Poisson interval has lower bound 0 at no event.
#   [2] seasonality_mon_harmonic(): zero draws give a flat curve at 1; a
#       cosine coefficient c gives exp(c cos(2 pi m / 12)) (mean of h is 0),
#       peaking in December; the interval brackets the mean.
#   [3] seasonality_peaks(): sparse months take no part in the observed peak;
#       in_rainy and the signed circular lag from the centre of the nearest
#       rainy season, for a season crossing the year end and a bimodal group.
#   [4] seasonality_frames() / build_figS_seasonality() build a 4-row
#       patchwork in display order without error.
#   [5] a fitted cessation curve is drawn in row c.
#   [6] occ_seas_calendar_draws() (occ_seasonality v2): each gid's rotation
#       by min_seas, the gid mean per draw; its column means equal the v1
#       rule (mean over gids of the rotated posterior mean).
#   [7] oni_composite_cycle() (Figure 1c's composite ONI): peak alignment
#       (Jul-Dec peak = year Y, Jan-Jun peak = the year before), incomplete
#       cycles left out, two episodes in one cycle counted once.
#
# HARD, quick-safe (no fit, a few seconds).
# Run: Rscript analysis/tests/test_results_seasonality.R

suppressPackageStartupMessages({ library(dplyr) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
# The helpers are 09's (results_derive.R); sourced into a throwaway env with
# their call-time deps (gid.R: applyAdminLevelFilter; clim_compile.R:
# circularMonthRuns; base.R: %||%).
.md <- new.env()
for (m in c("base", "gid", "clim_compile", "rain_enso")) sys.source(sprintf("analysis/utils/%s.R", m), envir = .md)
.rd <- parse("analysis/utils/results_derive.R", keep.source = FALSE)
for (e in .rd) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                   as.character(e[[2]]) %in% c(".SEAS_DAYS_PER_MONTH", "seasonality_pois_ci",
                                               "seasonality_track_rates", "seasonality_mon_harmonic",
                                               "seasonality_peaks", "occ_seas_calendar_draws",
                                               "oni_composite_cycle")) eval(e, .md)
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
DPM <- 365.25 / 12

cat("# [1] seasonality_track_rates\n")
cg <- c(AAA = "g1", BBB = "g1", CCC = "g2", ZZZ = "other")
# survival panel: unit-month rows at risk; AAA_1 is an ADM1 parent with an
# ADM2 child (AAA_1_1) and leaves under ADM1-ADM2; ZZZ is outside the groups.
surv <- tibble(
  country = c("AAA", "AAA", "AAA", "AAA", "BBB", "BBB", "CCC", "AAA", "ZZZ"),
  admin_level = c("ADM2", "ADM2", "ADM2", "ADM2", "ADM1", "ADM1", "ADM1", "ADM1", "ADM1"),
  ADM1 = c("AAA_1", "AAA_1", "AAA_1", "AAA_1", "BBB_1", "BBB_1", "CCC_1", "AAA_1", "ZZZ_1"),
  t_start = as.Date(c("2015-01-01", "2016-01-01", "2015-02-01", "2016-02-01", "2015-01-01", "2015-02-01",
                      "2015-01-01", "2015-01-01", "2015-01-01")),
  outbreak = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE),
  exposure_days = c(10, 31, 28, 28, 15, 28, 31, 31, 31))
dur <- tibble(country = c("AAA", "AAA", "BBB", "CCC"), admin_level = "ADM2",
              ADM1 = c("AAA_1", "AAA_1", "BBB_1", "CCC_1"),
              month = as.Date(c("2015-01-01", "2015-02-01", "2015-01-01", "2015-03-01")),
              ended = c(0L, 1L, 1L, 0L), exposure = c(31, 14, 20, 31))
mon <- tibble(country = c("AAA", "AAA", "AAA", "BBB", "CCC"), admin_level = "ADM2",
              ADM1 = c("AAA_1", "AAA_1", "AAA_1", "BBB_1", "CCC_1"),
              month = as.Date(c("2015-01-01", "2016-01-01", "2017-01-01", "2015-01-01", "2015-01-01")),
              cases_month = c(10, 0, 30, 50, 7), cases_gap = c(FALSE, TRUE, FALSE, FALSE, FALSE))
rt <- .md$seasonality_track_rates(surv, dur, mon, cg, c("g1", "g2"), "ADM1-ADM2", min_n = 2L, interval = 0.95)
g <- function(grp, tr, m) rt[rt$enso_group == grp & rt$track == tr & rt$month == m, , drop = FALSE]
expect(nrow(rt) == 2 * 3 * 12 && !any(rt$enso_group == "other"),
       "every group x track x month present, countries outside the groups dropped")
o1 <- g("g1", "occurrence", 1)
expect(o1$n == 2 && near(o1$exposure, (10 + 31 + 15) / DPM) && near(o1$rate, 1000 * 2 / ((10 + 31 + 15) / DPM)),
       "g1 January onsets: the ADM1 parent AAA_1 dropped, rate per 1,000 unit-months at risk")
o2 <- g("g1", "occurrence", 2)
expect(o2$n == 0 && o2$lo == 0 && near(o2$hi, 1000 * stats::qgamma(0.975, 1) / ((28 + 28 + 28) / DPM)),
       "no onset: exact Poisson interval with lower bound 0")
expect(o1$sparse == FALSE && o2$sparse == TRUE, "sparse = n < min_n")
expect(g("g1", "occurrence", 3)$n == 0 && is.na(g("g1", "occurrence", 3)$rate),
       "a month without rows: n 0, rate NA")
c2 <- g("g1", "cessation", 2)
expect(c2$n == 1 && near(c2$rate, 100 / (14 / DPM)), "g1 February cessation per 100 episode-months")
k1 <- g("g1", "cases", 1)
expect(k1$n == 3 && k1$rate == 30 && near(k1$lo, 20) && near(k1$hi, 40),
       "cases: the cases_gap month left out of the median and quartiles")
expect(g("g2", "cases", 1)$rate == 7 && g("g2", "cessation", 3)$n == 0, "g2 kept apart from g1")

cat("# [2] seasonality_mon_harmonic\n")
v <- c("b_sin_month", "b_cos_month", "b_sin2_month", "b_cos2_month")
z <- matrix(0, 4, 4, dimnames = list(NULL, v))
h0 <- .md$seasonality_mon_harmonic(z)
expect(nrow(h0) == 12 && near(h0$rel, rep(1, 12)) && near(h0$lo, rep(1, 12)) && near(h0$hi, rep(1, 12)),
       "zero draws: flat at 1")
cc <- matrix(c(0, 0.5, 0, 0), 1, 4, dimnames = list(NULL, v))
h1 <- .md$seasonality_mon_harmonic(cc)
expect(near(h1$rel, exp(0.5 * cos(2 * pi * (1:12) / 12))) && which.max(h1$rel) == 12L,
       "cosine coefficient: exp(c cos(2 pi m / 12)), peak in December")
set.seed(3)
dr <- matrix(rnorm(400 * 4, c(0.2, -0.3, 0.1, 0), 0.05), 400, 4, byrow = TRUE, dimnames = list(NULL, v))
h2 <- .md$seasonality_mon_harmonic(dr)
expect(all(h2$lo <= h2$rel & h2$rel <= h2$hi), "the interval brackets the posterior mean")

cat("# [3] seasonality_peaks\n")
rates <- tibble(enso_group = rep(c("s", "b"), each = 12), track = "occurrence", month = rep(1:12, 2),
                rate = c(1, 5, 2, 1, 1, 1, 1, 1, 9, 1, 1, 1,     # s: Sep sparse, so Feb wins
                         1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1, 1),    # b: peak in August
                sparse = c(rep(FALSE, 8), TRUE, rep(FALSE, 3), rep(FALSE, 12)))
fitted <- tibble(enso_group = c("s", "b"), track = "occurrence", month = c(7L, 11L), rel = 2)
rain <- tibble(enso_group = rep(c("s", "b"), each = 12), month = rep(1:12, 2),
               in_season = c(1:12 %in% c(11, 12, 1, 2, 3),      # s: Nov-Mar, centre January
                             1:12 %in% c(3, 4, 5, 10, 11)))      # b: Mar-May (centre Apr) + Oct-Nov (centre 10.5)
pk <- .md$seasonality_peaks(rates, fitted, rain)
p <- function(grp, src) pk[pk$enso_group == grp & pk$source == src, , drop = FALSE]
expect(p("s", "observed")$peak_month == 2 && p("s", "observed")$in_rainy && near(p("s", "observed")$lag_from_rain, 1),
       "a sparse month takes no part in the peak; Feb is +1 month from a Nov-Mar season's centre")
expect(p("s", "fitted")$peak_month == 7 && !p("s", "fitted")$in_rainy && near(abs(p("s", "fitted")$lag_from_rain), 6),
       "July is 6 months from January, outside the season")
expect(p("b", "observed")$peak_month == 8 && near(p("b", "observed")$lag_from_rain, -2.5) &&
         !p("b", "observed")$in_rainy,
       "bimodal group: the lag is from the nearest season's centre (Aug is 2.5 months before 10.5)")
expect(p("b", "fitted")$in_rainy && near(p("b", "fitted")$lag_from_rain, 0.5) &&
         p("b", "fitted")$rainy_months == "Mar-May, Oct-Nov" &&
         p("s", "observed")$rainy_months == "Nov-Mar",
       "a peak inside the second season, and the rainy months listed")

cat("# [4] seasonality_frames / build_figS_seasonality\n")
source("analysis/utils/figures_supp.R")
`%||%` <- .md$`%||%`
art <- list(
  rain = tibble(enso_group = rep(c("g1", "g2"), each = 12), month = rep(1:12, 2), mm = 50 + 10 * (1:24 %% 12),
                in_season = rep(1:12 %in% 5:9, 2), share = 0.6),
  rates = rt,
  fitted = bind_rows(
    tibble(enso_group = rep(c("g1", "g2"), each = 12), track = "occurrence", month = rep(1:12, 2), rel = 1,
           lo = NA_real_, hi = NA_real_, scaled = 20, scaled_lo = NA_real_, scaled_hi = NA_real_),
    tibble(enso_group = rep(c("g1", "g2"), each = 12), track = "cases", month = rep(1:12, 2), rel = 1,
           lo = 0.9, hi = 1.1, scaled = 20, scaled_lo = 18, scaled_hi = 22)))
fr <- seasonality_frames(art, c("g2", "g1", "zz"))
expect(identical(levels(fr$rain$enso_group), c("g2", "g1")) && nrow(fr$shade) == 10,
       "frames in display order, the rainy months as the band")
pw <- tryCatch(build_figS_seasonality(fr, c(g1 = "G1", g2 = "G2"), c(g1 = "red", g2 = "blue"),
                                      ggplot2::theme_bw(), min_n = 2L),
               error = function(e) e)
expect(inherits(pw, "patchwork") && length(pw$patches$plots) + 1L == 4L, "a 4-row patchwork")
ok <- tryCatch({ grDevices::pdf(NULL); print(pw); grDevices::dev.off(); TRUE }, error = function(e) conditionMessage(e))
expect(isTRUE(ok), paste("the figure draws:", ok))

cat("# [5] a fitted cessation curve (duration.seasonal_term harmonic2) is drawn in row c\n")
art_c <- art
art_c$fitted <- bind_rows(art$fitted,
  tibble(enso_group = rep(c("g1", "g2"), each = 12), track = "cessation", month = rep(1:12, 2), rel = 1,
         lo = 0.8, hi = 1.2, scaled = 10, scaled_lo = 8, scaled_hi = 12))
pw_c <- build_figS_seasonality(seasonality_frames(art_c, c("g1", "g2")), c(g1 = "G1", g2 = "G2"),
                               c(g1 = "red", g2 = "blue"), ggplot2::theme_bw(), min_n = 2L)
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
row_c0 <- geoms(pw[[3]]); row_c1 <- geoms(pw_c[[3]])
expect(!any(c("GeomLine", "GeomRibbon") %in% row_c0) && all(c("GeomLine", "GeomRibbon") %in% row_c1),
       "row c: no fitted curve without cessation rows, line + 95% band with them")
ok_c <- tryCatch({ grDevices::pdf(NULL); print(pw_c); grDevices::dev.off(); TRUE }, error = function(e) conditionMessage(e))
expect(isTRUE(ok_c), paste("the figure with a cessation curve draws:", ok_c))

cat("# [6] occ_seas_calendar_draws\n")
set.seed(1)
r <- matrix(exp(rnorm(40 * 12, sd = 0.3)), 40, 12)
ms <- c(1L, 4L, 4L, 11L)
m6 <- .md$occ_seas_calendar_draws(r, ms)
rot <- function(v, s) { o <- numeric(12); o[((s - 1L) + 0:11) %% 12L + 1L] <- v; o }
expect(near(m6[1, ], rowMeans(sapply(ms, function(s) rot(r[1, ], s)))), "draw 1: the mean of the gids' rotations")
v1 <- rowMeans(sapply(ms, function(s) rot(colMeans(r), s)))
expect(near(colMeans(m6), v1), "column means = v1's mean over gids of the rotated posterior mean")
expect(near(.md$occ_seas_calendar_draws(r, 1L), r), "min_seas 1: season month = calendar month")

cat("# [7] oni_composite_cycle\n")
dates <- seq(as.Date("2000-01-01"), as.Date("2004-12-01"), by = "1 month")
obs <- tibble(date = dates, enso_index = 0)
obs$enso_index[obs$date == as.Date("2000-12-01")] <- 2      # El Nino peak Dec 2000 -> Y = 2000
obs$enso_index[obs$date == as.Date("2003-02-01")] <- -1.5   # La Nina peak Feb 2003 -> Y = 2002
obs$enso_index[obs$date == as.Date("2004-11-01")] <- 1      # El Nino peak Nov 2004: Y+1 not observed
ph <- tibble(phase = c("elnino", "nina", "elnino", "elnino"),
             TL = as.Date(c("2000-09-15", "2002-11-15", "2004-08-15", "2000-10-15")),
             TR = as.Date(c("2001-02-15", "2003-04-15", "2004-12-15", "2001-01-15")))
cyc <- .md$oni_composite_cycle(obs, ph)
en <- dplyr::filter(cyc, kind == "elnino"); ln <- dplyr::filter(cyc, kind == "nina")
expect(nrow(en) == 24 && all(en$n == 1) && en$mean[en$x == 12] == 2 && sum(en$mean != 0) == 1,
       "El Nino: one cycle (Dec peak at x = 12; the 2004 cycle incomplete, the duplicate counted once)")
expect(nrow(ln) == 24 && ln$mean[ln$x == 14] == -1.5, "La Nina: Feb peak of Y+1 at x = 14 (Y = the year before)")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
