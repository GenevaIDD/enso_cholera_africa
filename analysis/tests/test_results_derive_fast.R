# test_results_derive_fast.R -- the 09 derivation speed-ups give identical
# output (HARD gate, quick).
#
# At 10^4 replicates per arm 09's replicate-level passes run over ~10^8-row
# tables, so they go through data.table (.sum_by_keys(),
# .anti_join_rows(), .replicate_keys()), the ribbon quantiles take one
# quantile() call per column, report_forecast sums its window totals once per
# metric x set and per-replicate populations by tapply(), and 09 loads the runs
# once per job (ctx$RUN_LOAD_MEMO).
#
#   [1] .sum_by_keys() (sum, max) == dplyr, with missing / factor keys;
#       .anti_join_rows() == dplyr::anti_join (missing keys, no match, empty).
#   [2] every changed derivation helper == its version at commit cecdea9 (the
#       last before these changes) on a synthetic two-group, two-set run:
#       window_cum_per_rep, explosive_replicate_keys, pair_adm1_cum,
#       per_rep_totals, .case_onset_ribbon, af_posterior_window_reps,
#       attrib_exposure_tables and attrib_metrics_tables (bootstrap on and off).
#   [3] load_group_arm_runs() memo: same request -> the stored load, a
#       different request -> a new load.
#   [4] the second round (pair_adm1_cum() fed one pair_adm1_cum_table() per
#       cell, the bound run frames memoised with the load, the country -> region
#       lookup by match(), af_posterior_draw_table()'s replicate passes in
#       data.table) == commit 30eb7b7 (the last before it): pair_adm1_cum with
#       and without the table (several analog sets on the forecast rows too),
#       .lookup_join() == left_join, .bound_run_frames() memo, and
#       af_posterior_draw_table() / af_posterior_tables() on a synthetic
#       two-group, cycled-draw design.
#   [5] the attribution cell table (attrib_cells()):
#       every window x guard loop it replaced == commit 0008799.
#
# Synthetic inputs only. Run: Rscript analysis/tests/test_results_derive_fast.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/results_core.R")
source("analysis/utils/results_derive.R")

cat("\n[1] .sum_by_keys() / .anti_join_rows() == dplyr\n")
set.seed(7)
n <- 20000L
df <- tibble(run = sample(c("forecast", "counterfactual"), n, TRUE),
             analog_set = sample(c(NA, "all", "neutral"), n, TRUE),
             enso_group = sample(c("enso_east", "enso_west"), n, TRUE),
             fac = factor(sample(c("b", "a", NA), n, TRUE), levels = c("b", "a")),
             replicate_id = sample(1:60, n, TRUE), x = as.numeric(rpois(n, 4)), y = rpois(n, 1))
for (op in c("sum", "max")) for (keys in list(c("enso_group", "run", "analog_set", "replicate_id"), c("fac", "run"))) {
  f <- if (op == "sum") sum else max
  ref <- df %>% group_by(across(all_of(keys))) %>% summarise(a = f(x), b = f(y), .groups = "drop")
  expect(identical(.sum_by_keys(df, keys, c(a = "x", b = "y"), op), ref),
         sprintf(".sum_by_keys %s by %s", op, paste(keys, collapse = " x ")))
}
by <- c("enso_group", "run", "analog_set", "replicate_id")
y1 <- distinct(df[sample(n, 300), by])
expect(identical(.anti_join_rows(df, y1, by), anti_join(df, y1, by = by)), ".anti_join_rows (missing keys match)")
y0 <- tibble(enso_group = "none", run = "x", analog_set = NA_character_, replicate_id = 0L)
expect(identical(.anti_join_rows(df, y0, by), anti_join(df, y0, by = by)), ".anti_join_rows, no match")
expect(identical(.anti_join_rows(df, y1[0, ], by), anti_join(df, y1[0, ], by = by)), ".anti_join_rows, empty y")

cat("\n[2] derivation helpers == cecdea9\n")
old_src <- tryCatch(system2("git", c("show", "cecdea9:analysis/utils/results_derive.R"), stdout = TRUE, stderr = FALSE),
                    warning = function(w) NULL, error = function(e) NULL)
if (is.null(old_src)) {
  cat("  SKIP  cecdea9's results_derive.R not readable (git unavailable?)\n")
} else {
  old <- new.env(parent = globalenv())
  for (e in parse(text = old_src)) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, envir = old)
  # Synthetic run: 2 groups x 2 countries x 2 ADM1s, 2 analog sets, 40 replicates,
  # 5 months; a full replicate x ADM1 x month grid per arm.
  adm1 <- tibble(adm1_unit = c("KEN_1", "KEN_2", "UGA_1", "UGA_2", "NGA_1", "NGA_2", "GHA_1", "GHA_2"),
                 enso_group = rep(c("enso_east", "enso_nga"), each = 4))
  months <- seq(as.Date("2026-05-01"), by = "month", length.out = 5)
  arms <- tibble(run = c("forecast", "counterfactual", "counterfactual"), analog_set = c(NA, "all", "neutral"))
  am <- tidyr::expand_grid(arms, adm1, replicate_id = 1:40, month = months) %>%
    mutate(unit = ifelse(enso_group == "enso_east", "hmmCE", "hmmW"),
           cases = as.numeric(rnbinom(n(), mu = ifelse(run == "forecast", 60, 35), size = .3)),
           outbreaks = rbinom(n(), 1, .25), group_horizon_end = as.Date("2026-09-01"))
  am$cases[am$replicate_id == 3 & am$enso_group == "enso_east" & am$run == "forecast"] <- 1e5  # an explosive replicate
  cm <- am %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
    group_by(run, analog_set, enso_group, country, replicate_id, month, group_horizon_end) %>%
    summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
  st <- as.Date("2026-06-01")
  expect(identical(window_cum_per_rep(cm, st, "enso_group"), old$window_cum_per_rep(cm, st, "enso_group")),
         "window_cum_per_rep (analog_set NA on the forecast rows)")
  expect(identical(explosive_replicate_keys(cm, st, 2e5), old$explosive_replicate_keys(cm, st)), "explosive_replicate_keys")
  expect(identical(pair_adm1_cum(am, "cases", "neutral"), old$pair_adm1_cum(am, "cases", "neutral")), "pair_adm1_cum")
  prt <- mutate(cm, n_months = 5L, pop = ifelse(country == "KEN", NA, 1e6 + replicate_id))
  expect(identical(per_rep_totals(prt, "enso_group"), old$per_rep_totals(prt, "enso_group")), "per_rep_totals (NA pop key)")
  rib <- cm %>% mutate(run = factor(run, levels = c("forecast", "counterfactual"))) %>%
    group_by(run, analog_set, month, enso_group)
  expect(identical(.case_onset_ribbon(rib), old$.case_onset_ribbon(rib)), ".case_onset_ribbon")
  # gmax (per-gid explosive cap) is the one added element.
  expect(identical(af_posterior_window_reps(am, cm, "attribution", st, cap = 2e5)[c("reps", "explosive")],
                   old$af_posterior_window_reps(am, cm, "attribution", st)), "af_posterior_window_reps")
  am_nounit <- select(am, -unit)
  expect(identical(af_posterior_window_reps(am_nounit, cm, "post_fork", st, cap = 2e5)[c("reps", "explosive")],
                   old$af_posterior_window_reps(am_nounit, cm, "post_fork", st)), "af_posterior_window_reps, no unit column")
  pop <- tibble(adm1_unit = adm1$adm1_unit, pop = c(1e5, 2e5, NA, 3e5, 4e5, 5e5, 6e5, 7e5))
  gs <- c(enso_east = "East", enso_nga = "NGA+GHA")
  for (B in c(20L, 0L)) {
    a_new <- suppressMessages(attrib_exposure_tables(am, cm, pop, gs, c("all", "neutral"), st, c(.25, .5), .66, .8, cap = 2e5,
                                                     n_boot = B))
    a_old <- suppressMessages(old$attrib_exposure_tables(am, cm, pop, gs, c("all", "neutral"), st, c(.25, .5), .66, .8,
                                                         n_boot = B))
    expect(identical(a_new, a_old), sprintf("attrib_exposure_tables, B = %d", B))
    m_new <- suppressMessages(attrib_metrics_tables(am, cm, pop, gs, c("all", "neutral"), st, 50, .05, c(1, 2), cap = 2e5,
                                                    n_boot = B))
    m_old <- suppressMessages(old$attrib_metrics_tables(am, cm, pop, gs, c("all", "neutral"), st, 50, .05, c(1, 2),
                                                        n_boot = B))
    expect(identical(m_new, m_old), sprintf("attrib_metrics_tables, B = %d", B))
  }
}

cat("\n[3] load_group_arm_runs() memo\n")
calls <- 0L
local({
  .load_group_arm_runs_impl <<- function(ctx, groups, sets, attrib_start, verbose = FALSE, win = NULL) {
    calls <<- calls + 1L; list(groups = groups, sets = sets, n = calls)
  }
})
ctx <- list2env(list(MAN = list(a = 1), RES = list(), SIM_ROOT = "x", SMOOTHING_PARTS_DIR = "y", ARM_ALIGNMENT = "stitch",
                     RUN_LOAD_MEMO = new.env(parent = emptyenv())))
r1 <- suppressMessages(load_group_arm_runs(ctx, c("g1", "g2"), c("all", "neutral"), as.Date("2026-05-01")))
r2 <- suppressMessages(load_group_arm_runs(ctx, c("g1", "g2"), c("all", "neutral"), as.Date("2026-05-01"), verbose = TRUE))
expect(calls == 1L && identical(r1, r2), "same request: one load, the stored result returned")
r3 <- suppressMessages(load_group_arm_runs(ctx, c("g1", "g2"), "neutral", as.Date("2026-05-01")))
expect(calls == 2L && identical(r3$sets, "neutral"), "other sets: a new load")
ctx2 <- list2env(c(as.list(ctx), list(MAN = list(a = 2))))
r4 <- suppressMessages(load_group_arm_runs(ctx2, c("g1", "g2"), "neutral", as.Date("2026-05-01")))
expect(calls == 3L, "another manifest (an event view): a new load")
ctx$RUN_LOAD_MEMO <- NULL
r5 <- suppressMessages(load_group_arm_runs(ctx, c("g1", "g2"), "neutral", as.Date("2026-05-01")))
expect(calls == 4L, "no memo in the context: loads directly")

cat("\n[4] second round == 30eb7b7\n")
old_src2 <- tryCatch(system2("git", c("show", "30eb7b7:analysis/utils/results_derive.R"), stdout = TRUE, stderr = FALSE),
                     warning = function(w) NULL, error = function(e) NULL)
if (is.null(old_src2)) {
  cat("  SKIP  30eb7b7's results_derive.R not readable (git unavailable?)\n")
} else {
  old2 <- new.env(parent = globalenv())
  for (e in parse(text = old_src2)) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, envir = old2)
  set.seed(11)
  adm1 <- tibble(adm1_unit = c("KEN_1", "KEN_2", "UGA_1", "NGA_1", "NGA_2", "GHA_1"),
                 enso_group = c(rep("enso_east", 3), rep("enso_nga", 3)))
  months <- seq(as.Date("2026-05-01"), by = "month", length.out = 4)
  arms <- tibble(run = c("forecast", "counterfactual", "counterfactual", "counterfactual"),
                 analog_set = c(NA, "all", "neutral", NA))
  am <- tidyr::expand_grid(arms, adm1, replicate_id = 1:60, month = months) %>%
    mutate(unit = ifelse(enso_group == "enso_east", "hmmCE", ifelse(adm1_unit == "GHA_1", "hmmG", "hmmW")),
           cases = as.numeric(rnbinom(n(), mu = ifelse(run == "forecast", 50, 30), size = .4)),
           outbreaks = rbinom(n(), 1, .3), group_horizon_end = as.Date("2026-08-01"))
  am$cases[sample(nrow(am), 40)] <- NA
  cum <- pair_adm1_cum_table(am)
  for (vc in c("cases", "outbreaks")) for (aset in c("all", "neutral"))
    expect(identical(pair_adm1_cum(am, vc, aset, cum = cum), old2$pair_adm1_cum(am, vc, aset)),
           sprintf("pair_adm1_cum(cum = table) %s / %s", vc, aset))
  am2 <- mutate(am, analog_set = ifelse(run == "forecast", sample(c("x", "y", NA), n(), TRUE), analog_set))
  expect(identical(pair_adm1_cum(am2, "outbreaks", "all", cum = pair_adm1_cum_table(am2)),
                   old2$pair_adm1_cum(am2, "outbreaks", "all")), "pair_adm1_cum(cum), several sets on the forecast rows")

  cm <- am %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
    group_by(run, analog_set, enso_group, country, replicate_id, month, group_horizon_end) %>%
    summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
  cr <- tibble(country = c("GHA", "KEN", "NGA", "ZZZ"), AFRO_region = c("W", "E", "W", "X"))
  cm_na <- bind_rows(cm, mutate(cm[1:3, ], country = NA_character_))
  expect(identical(.lookup_join(cm_na, cr, "country"), left_join(cm_na, cr, by = "country")),
         ".lookup_join == left_join (unmatched and missing keys)")
  cr_dup <- bind_rows(cr, tibble(country = "KEN", AFRO_region = "E2"))
  expect(identical(suppressWarnings(.lookup_join(cm, cr_dup, "country")),
                   suppressWarnings(left_join(cm, cr_dup, by = "country"))), ".lookup_join, repeated keys -> left_join")

  runs <- list(adm1_month = split(am, am$enso_group), country_month = split(cm, cm$enso_group))
  ctx <- list2env(list(RUN_LOAD_MEMO = new.env(parent = emptyenv())))
  ctx$RUN_LOAD_MEMO$value <- runs
  b1 <- .bound_run_frames(ctx, runs)
  expect(identical(b1, list(adm1_month = bind_rows(runs$adm1_month), country_month = bind_rows(runs$country_month))),
         ".bound_run_frames == bind_rows")
  expect(identical(.bound_run_frames(ctx, runs), b1) && !is.null(ctx$RUN_LOAD_MEMO$bound), ".bound_run_frames: memoised with the load")
  other <- list(adm1_month = runs$adm1_month[1], country_month = runs$country_month[1])
  b2 <- .bound_run_frames(ctx, other)
  expect(nrow(b2$adm1_month) < nrow(b1$adm1_month) && identical(ctx$RUN_LOAD_MEMO$bound, b1),
         ".bound_run_frames: another load is bound, the memo kept")

  # af_posterior: 3 units, 7 cycled draws over 60 replicates (unit hmmG has
  # one country), a replicate with a single member per draw cell for v0().
  am_af <- filter(am, run == "forecast" | analog_set %in% c("all", "neutral"))
  cm_af <- filter(cm, run == "forecast" | analog_set %in% c("all", "neutral"))
  wr <- af_posterior_window_reps(am_af, cm_af, "attribution", as.Date("2026-05-01"), cap = 2e5)
  draws <- distinct(wr$reps, cell, unit, replicate_id) %>% mutate(draw = ((replicate_id - 1L) %% 7L) + 1L)
  for (aset in c("all", "neutral")) for (metric in c("cases", "onsets")) {
    pairs <- af_posterior_pairs(wr$reps, wr$explosive, draws, aset, metric, TRUE)
    for (by in c("enso_group", "country"))
      expect(identical(af_posterior_draw_table(pairs, by), old2$af_posterior_draw_table(pairs, by)),
             sprintf("af_posterior_draw_table %s / %s by %s", aset, metric, by))
  }
  p1 <- af_posterior_pairs(wr$reps, wr$explosive, draws, "all", "cases", TRUE)
  p1 <- p1[!duplicated(p1[, c("cell", "unit", "draw")]) | p1$unit != "hmmG", ]
  expect(identical(af_posterior_draw_table(p1, "country"), old2$af_posterior_draw_table(p1, "country")),
         "af_posterior_draw_table, one replicate per draw in a unit")
  gs <- c(enso_east = "East", enso_nga = "NGA+GHA")
  # The legacy path (kept side by side): the speed-up equivalence is against
  # the older code.
  t_new <- suppressMessages(af_posterior_tables(wr$reps, wr$explosive, draws, c("all", "neutral"), gs,
                                                min_reps = 5L, n_mc = 200L, africa = "independent", guard_scope = "group"))
  t_old <- suppressMessages(old2$af_posterior_tables(wr$reps, wr$explosive, draws, c("all", "neutral"), gs,
                                                     min_reps = 5L, n_mc = 200L))
  expect(identical(t_new, t_old), "af_posterior_tables")
}

cat("\n[5] attribution cells (attrib_cells(), review R1) == 0008799\n")
# The window x guard loops of attrib_exposure_windows,
# attrib_metrics_windows, attrib_exposure_posterior_windows, event_attrib_cells
# and the af_posterior stash (.afp_forecast_cells) read one cell table. At the
# back-compat settings every output equals the loops' at 0008799.
old_src3 <- tryCatch(system2("git", c("show", "0008799:analysis/utils/results_derive.R"), stdout = TRUE, stderr = FALSE),
                     warning = function(w) NULL, error = function(e) NULL)
if (is.null(old_src3)) {
  cat("  SKIP  0008799's results_derive.R not readable (git unavailable?)\n")
} else {
  old3 <- new.env(parent = globalenv())
  for (e in parse(text = old_src3)) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, envir = old3)
  set.seed(5)
  adm1 <- tibble(adm1_unit = c("KEN_1", "KEN_2", "UGA_1", "NGA_1", "NGA_2", "GHA_1"),
                 enso_group = c(rep("enso_east", 3), rep("enso_nga", 3)))
  months <- seq(as.Date("2026-05-01"), by = "month", length.out = 4)
  arms <- tibble(run = c("forecast", "counterfactual", "counterfactual"), analog_set = c(NA, "all", "neutral"))
  am <- tidyr::expand_grid(arms, adm1, replicate_id = 1:42, month = months) %>%
    mutate(unit = ifelse(enso_group == "enso_east", "hmmCE", "hmmW"),
           cases = as.numeric(rnbinom(n(), mu = ifelse(run == "forecast", 50, 30), size = .4)),
           outbreaks = rbinom(n(), 1, .3), group_horizon_end = as.Date("2026-08-01"))
  am$cases[am$replicate_id == 4 & am$enso_group == "enso_east" & am$run == "forecast"] <- 1e5    # explosive, forecast
  am$cases[am$replicate_id == 9 & am$enso_group == "enso_nga" & am$analog_set %in% "all"] <- 1e5 # explosive, cf
  cm <- am %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
    group_by(run, analog_set, enso_group, country, replicate_id, month, group_horizon_end) %>%
    summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
  pop <- tibble(adm1_unit = adm1$adm1_unit, pop = c(1e5, 2e5, NA, 3e5, 4e5, 5e5))
  gs <- c(enso_east = "East", enso_nga = "NGA+GHA")
  a0 <- as.Date("2026-05-01"); pf <- as.Date("2026-07-01")
  cells <- attrib_cell_table(a0, pf, cap = 2e5)
  expect(identical(cells$window, c("attribution", "attribution", "post_fork", "post_fork")) &&
           identical(cells$explosive_guard, c(TRUE, FALSE, TRUE, FALSE)) && identical(cells$cap, c(2e5, Inf, 2e5, Inf)) &&
           identical(.cell_starts(cells), old3$.afp_forecast_cells(list(attrib_start = a0), pf)),
         "cell table: windows x guards in the former loop order, starts == .afp_forecast_cells()")
  sets <- c("all", "neutral")
  ae_new <- suppressMessages(attrib_exposure_windows(am, cm, pop, gs, sets, cells, c(.25, .5), .66, .8, n_boot = 15L))
  ae_old <- suppressMessages(old3$attrib_exposure_windows(am, cm, pop, gs, sets, a0, pf, c(.25, .5), .66, .8, n_boot = 15L))
  expect(identical(ae_new, ae_old), "attrib_exposure_windows")
  mw_new <- suppressMessages(attrib_metrics_windows(am, cm, pop, gs, sets, cells, 50, .05, c(1, 2), n_boot = 15L,
                                                    ipcc_mags = c(.25, .5), robust_exceed = .66, robust_agree = .8))
  mw_old <- suppressMessages(old3$attrib_metrics_windows(am, cm, pop, gs, sets, a0, pf, 50, .05, c(1, 2), n_boot = 15L,
                                                         ipcc_mags = c(.25, .5), robust_exceed = .66, robust_agree = .8))
  expect(identical(mw_new, mw_old), "attrib_metrics_windows")
  draws <- distinct(am, unit, replicate_id) %>% mutate(draw = ((replicate_id - 1L) %% 3L) + 1L)
  pw_new <- suppressMessages(attrib_exposure_posterior_windows(am, cm, pop, gs, sets, cells, c(.25, .5), draws,
                                                               min_reps = 5L, n_mc = 100L, africa = "independent", guard_scope = "group"))
  pw_old <- suppressMessages(old3$attrib_exposure_posterior_windows(am, cm, pop, gs, sets, .cell_starts(cells),
                                                                    c(.25, .5), draws, min_reps = 5L, n_mc = 100L))
  expect(identical(pw_new, pw_old), "attrib_exposure_posterior_windows")
  inp <- list(adm1_month = am, country_month = cm, attrib_start = a0, groups = c("enso_east", "enso_nga"),
              adm1_pop = pop, win = list(factual_start = pf))
  sb_new <- af_posterior_stash_build(inp, cells)
  sb_old <- old3$af_posterior_stash_build(inp, old3$.afp_forecast_cells(inp, pf))
  expect(identical(sb_new[c("groups", "cells", "reps", "explosive")], sb_old), "af_posterior_stash_build")
  dd <- bind_rows(lapply(names(sb_new$cells), function(k) mutate(draws, cell = k, .before = 1L)))
  expect(identical(af_posterior_tables(sb_new$reps, sb_new$explosive, dd, sets, gs, guards = .cell_guards(cells),
                                       min_reps = 5L, n_mc = 100L, africa = "independent", guard_scope = "group"),
                   old3$af_posterior_tables(sb_old$reps, sb_old$explosive, dd, sets, gs, min_reps = 5L, n_mc = 100L)),
         "af_posterior_tables on the stash (guards from the cells)")
  ctx5 <- list2env(list(GROUP_SHORT_RUN = gs, PR_CASE_THRESHOLD = 50, PR_ALPHA = .05, SN_LEVELS = c(1, 2),
                        ATTRIB_MAGNITUDES = c(.25, .5), ROBUST_EXCEED = .66, ROBUST_AGREE = .8, BOOT_B = 15L,
                        EXPLOSIVE_CAP = 2e5, EXPLOSIVE_MODES = c("drop", "none")))
  ev_new <- suppressMessages(event_attrib_cells(inp, ctx5, sets))
  ev_old <- suppressMessages(old3$event_attrib_cells(inp, ctx5, sets))
  expect(identical(ev_new, ev_old), "event_attrib_cells")
  es_new <- af_posterior_stash_build(inp, event_cells(ctx5, a0, "ev1"))
  es_old <- old3$af_posterior_stash_build(inp, stats::setNames(a0, "ev1"))
  expect(identical(es_new[c("groups", "cells", "reps", "explosive")], es_old), "event af_posterior stash")
}

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
