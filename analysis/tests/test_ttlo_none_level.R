# test_ttlo_none_level.R -- the TTLO level "none" ("no recorded prior onset"),
# simulator side.
#
# HARD gate (see analysis/tests/run_all.R). Synthetic, no MCMC, no data.
#
#   B  bin_ttlo(): default = a frozen copy of the function without the level (the
#      sentinel bins to [6,Inf]); has_none maps the sentinel to "none".
#      ttlo_init_from_panel("none") = sentinel; ttlo_col_to_level().
#   D  .draws_has_ttlo_none(): bundle, grouped list, mixed set refused.
#   W  build_warmup_state(no_prior_onset =): NULL = unchanged; TRUE seeds the
#      sentinel only where the category is NA or [6,Inf]; frameFirstOnset() /
#      ttloNoPriorOnsetAt(); advance_ttlo() keeps the sentinel until an onset.
#   C  .ttlo_bin_c(): without the level the three emitted chains equal the
#      frozen literals character for character; precompute_ttlo_coefs() adds
#      ttlo_coef_b7 only with the level.
#   E  v1 and C with the level on synthetic draws: a large "none" coefficient
#      makes every no-prior unit onset in month 1 and no other unit, in the
#      v1 kernel and the Csnippet kernel alike; a
#      fit without the level is unaffected by first_onset.
#
# Run: Rscript analysis/tests/test_ttlo_none_level.R   (from the project root)

suppressPackageStartupMessages({ library(tibble); library(dplyr); library(lubridate); library(tidyr) })
source("analysis/utils/smoothing_core.R")   # base, simulator_core, spatpomp_core{,_c}
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
cat("# test_ttlo_none_level\n\n")

# ===========================================================================
cat("## B: bins\n")
bin_old <- function(days) {          # frozen copy, simulator_core.R at c032a12
  years <- (pmax(0, days) + 1) / 365
  as.character(cut(years, .TTLO_BREAKS, labels = .TTLO_LEVELS, include.lowest = TRUE, right = FALSE))
}
set.seed(1); d <- c(NA, -5, 0, 364, 365, 3000, runif(50, 0, 4000), .TTLO_NONE_DAYS, .TTLO_NONE_DAYS + 900)
expect(identical(bin_ttlo(d), bin_old(d)), "B1 default bin_ttlo() = the pre-item-11 function (sentinel -> [6,Inf])")
bn <- bin_ttlo(d, has_none = TRUE)
expect(identical(bn[1:56], bin_old(d)[1:56]) && all(tail(bn, 2) == "none"),
       "B2 has_none: sentinel (and beyond) -> none, everything else unchanged")
expect(identical(ttlo_init_from_panel(c("none", "[2,3)", NA)), c(.TTLO_NONE_DAYS, 2 * 365.25, 6 * 365.25)),
       "B3 ttlo_init_from_panel: none -> sentinel; other labels and NA unchanged")
expect(identical(ttlo_col_to_level("time_to_last_outbreaknone"), "none") &&
       identical(ttlo_col_to_level("time_to_last_outbreak6Inf"), "[6,Inf]") &&
       identical(ttlo_col_to_level("time_to_last_outbreak12"), "[1,2)"), "B4 ttlo_col_to_level")
expect(advance_ttlo(.TTLO_NONE_DAYS, 0L) >= .TTLO_NONE_DAYS && advance_ttlo(.TTLO_NONE_DAYS, 1L) == 0,
       "B5 advance_ttlo keeps the sentinel until an onset resets it to 0")

# ===========================================================================
cat("\n## D: .draws_has_ttlo_none()\n")
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2) %>%
  mutate(gid = sprintf("%s_%d_1", country, adm1n), adm1_unit = sprintf("%s_%d", country, adm1n),
         AFRO_region = "REG_W", tu = "hmmW", min_seas = 1L, pop = 1e5) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>% as_tibble()
lvl_terms <- c("(Intercept)", paste0("time_to_last_outbreak", c("12", "23", "34", "45", "56", "6Inf")))
d_lvl  <- make_draws_bundle_synthetic(frame_tu, S = 2L, active_terms = c(lvl_terms, "time_to_last_outbreaknone"),
                                      fe_overrides = c("(Intercept)" = -40, "time_to_last_outbreaknone" = 80))
d_nolv <- make_draws_bundle_synthetic(frame_tu, S = 2L, active_terms = lvl_terms,
                                      fe_overrides = c("(Intercept)" = -40))
expect(.draws_has_ttlo_none(d_lvl) && !.draws_has_ttlo_none(d_nolv) && !.draws_has_ttlo_none(NULL),
       "D1 single bundle")
expect(.draws_has_ttlo_none(list(a = d_lvl, b = d_lvl)) && !.draws_has_ttlo_none(list(a = d_nolv)),
       "D2 grouped list")
expect(!is.null(err_of(.draws_has_ttlo_none(list(a = d_lvl, b = d_nolv)))), "D3 mixed grouped set refused")

# ===========================================================================
cat("\n## W: warm start\n")
months <- as.Date("2024-01-01") %m+% months(0:3)
pre    <- as.Date("2023-12-01")
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = c(pre, months)) %>%
  left_join(select(frame_tu, gid, min_seas), by = "gid") %>%
  mutate(country = sub("_[0-9]+_[0-9]+$", "", gid), t_start_seas = as.integer(yday(t_start) - 1L),
         mean_precip = 0, dev_precip = 0, mean_temp_daymean = 0, dev_temp_daymean = 0,
         outbreak = 0L, ongoing_outbreak_lag1 = 0L, ongoing_outbreak_country_lag1 = 0L,
         ongoing_outbreak_tu_lag1 = 0L,
         time_to_last_outbreak_cat = c("AAA_1_1" = "[6,Inf]", "AAA_2_1" = NA, "BBB_1_1" = "[2,3)",
                                       "BBB_2_1" = "[6,Inf]")[gid]) %>% as_tibble()
horizon <- list(start = "2024-01-01", end = "2024-04-30")
w0 <- build_warmup_state(panel_tu, frame_tu, horizon)
expect(identical(build_warmup_state(panel_tu, frame_tu, horizon, no_prior_onset = NULL), w0),
       "W1 no_prior_onset NULL = the pre-item-11 warm start")
w1 <- build_warmup_state(panel_tu, frame_tu, horizon, no_prior_onset = c(TRUE, TRUE, TRUE, FALSE))
expect(identical(w1$u, c(.TTLO_NONE_DAYS, .TTLO_NONE_DAYS, w0$u[3], w0$u[4])),
       "W2 sentinel only where no prior onset AND category NA / [6,Inf]")
surv <- tibble(gid = c("AAA_1_1", "AAA_1_1", "BBB_1_1", "BBB_2_1"),
               t_start = as.Date(c("2015-03-01", "2018-01-01", "2024-02-01", "2012-06-01")),
               outbreak = c(1L, 1L, 1L, 0L))
fo <- frameFirstOnset(surv, frame_tu$gid)
expect(identical(fo, as.Date(c("2015-03-01", NA, "2024-02-01", NA))), "W3 frameFirstOnset (NA = never)")
expect(identical(ttloNoPriorOnsetAt(surv, frame_tu$gid, "2024-01-01"), c(FALSE, TRUE, TRUE, TRUE)),
       "W4 ttloNoPriorOnsetAt: an onset after the launch does not count")
expect(!is.null(err_of(.ttlo_no_prior_for_frame(frame_tu, "2024-01-01", TRUE))) &&
       is.null(.ttlo_no_prior_for_frame(frame_tu, "2024-01-01", FALSE)),
       "W5 the level without frame$first_onset is refused; no level -> NULL")

# ===========================================================================
cat("\n## C: Csnippet chains\n")
old1 <- paste0("    double yrs = (ud[u] < 0.0 ? 1.0/365.0 : (ud[u] + 1.0)/365.0);\n",
               "    int bin = (yrs<1.0)?0:(yrs<2.0)?1:(yrs<3.0)?2:(yrs<4.0)?3:(yrs<5.0)?4:(yrs<6.0)?5:6;\n")
old2 <- paste0("      double yrs = (ud[u] < 0.0 ? 1.0/365.0 : (ud[u]+1.0)/365.0);\n",
               "      int bin = (yrs<1.0)?0:(yrs<2.0)?1:(yrs<3.0)?2:(yrs<4.0)?3:(yrs<5.0)?4:(yrs<6.0)?5:6;\n")
old3 <- paste0("      double yrs = (ud_n < 0.0 ? 1.0/365.0 : (ud_n + 1.0)/365.0);\n",
               "      int bin = (yrs<1.0)?0:(yrs<2.0)?1:(yrs<3.0)?2:(yrs<4.0)?3:(yrs<5.0)?4:(yrs<6.0)?5:6;\n")
expect(identical(.ttlo_bin_c("ud[u]", "    ", TRUE), old1) &&
       identical(.ttlo_bin_c("ud[u]", "      ", FALSE), old2) &&
       identical(.ttlo_bin_c("ud_n", "      ", TRUE), old3),
       "C1 without the level: the three chains equal the former literals")
expect(grepl("(ud[u] >= 1000000000.0) ? 7 :", .ttlo_bin_c("ud[u]", has_none = TRUE), fixed = TRUE),
       "C2 with the level: sentinel branch to bin 7")
pc0 <- precompute_ttlo_coefs(.subset_draws(d_nolv, 1L), frame_tu)
pc1 <- precompute_ttlo_coefs(.subset_draws(d_lvl, 1L), frame_tu)
expect(identical(names(pc0), paste0("ttlo_coef_b", 0:6)) &&
       identical(names(pc1), paste0("ttlo_coef_b", 0:7)) && all(pc1$ttlo_coef_b7 == 80),
       "C3 precompute_ttlo_coefs: b7 only with the level")

# ===========================================================================
cat("\n## E: kernels with the level (v1, R spatPomp, C spatPomp)\n")
fr_lvl <- frame_tu; fr_lvl$first_onset <- fo           # AAA_2_1, BBB_2_1 never; BBB_1_1 after launch
no_prior <- ttloNoPriorOnsetAt(surv, frame_tu$gid, horizon$start)
# A unit seeded "none" has lp = -40 + 80 > 0 (onset ~ certain); any other unit
# lp = -40 (never). Category [2,3) (BBB_1_1) is not seeded.
pre_cat <- panel_tu$time_to_last_outbreak_cat[panel_tu$t_start == pre][match(frame_tu$gid,
                                                                             panel_tu$gid[panel_tu$t_start == pre])]
want <- unname(no_prior & (is.na(pre_cat) | pre_cat %in% "[6,Inf]"))
v1 <- suppressWarnings(simulate_one_tu(tu = "hmmW", frame_tu = fr_lvl, panel_tu = panel_tu,
                                       draws_occ = d_lvl, replicate_ids = 1:3, seed = 5L,
                                       output_part = tempfile("ttn_"), horizon = horizon,
                                       dry_run = TRUE)$monthly_summary)
y1 <- v1[as.Date(v1$month) == months[1], ]
y1 <- as.vector(tapply(y1$Y, y1$gid, max)[frame_tu$gid])
expect(identical(unname(y1 == 1L), want), "E1 v1 kernel: month-1 onsets = exactly the seeded no-prior units")
window <- list(start = "2024-01-01", end = "2024-05-01")
yout1 <- function(obj) {
  set.seed(9L)
  sm <- spatPomp::simulate(obj, format = "data.frame")
  r <- sm[sm$time == 1, ]
  stats::setNames(r$Yout, r$unitname)[frame_tu$gid]
}
o_c <- build_spatpomp_tu_c(fr_lvl, panel_tu, d_lvl, window, n_draws_in_filter = 1L, fixed_draw_idx = 1L)
yc <- yout1(o_c)
expect(identical(unname(yc == 1), want), "E3 C kernel: month-1 onsets = the seeded units (with the level)")
expect(any(grepl("ttlo_coef_b7", rownames(o_c@covar@table))), "E4 C build carries the b7 covariate")
# Without the level the frame's first_onset is inert: identical to the plain frame.
o_c0  <- build_spatpomp_tu_c(frame_tu, panel_tu, d_nolv, window, n_draws_in_filter = 1L, fixed_draw_idx = 1L)
o_c0f <- build_spatpomp_tu_c(fr_lvl, panel_tu, d_nolv, window, n_draws_in_filter = 1L, fixed_draw_idx = 1L)
expect(identical(o_c0@covar@table, o_c0f@covar@table) && identical(yout1(o_c0), yout1(o_c0f)) &&
       !any(grepl("ttlo_coef_b7", rownames(o_c0@covar@table))),
       "E5 a fit without the level: first_onset inert, no b7")

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
