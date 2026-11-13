# test_event_sub_windows.R -- the one-event
# La Nina 2020-23 with per-winter sub-windows.
#
# HARD gate, quick-safe: synthetic data only (the ONI-cache crossing check is
# skipped when generated_data/enso_oni_cache.rds is absent). No fit, no simulation.
#   [1] registry: enso_events_v2.yml parses; its single-winter events load with
#       peaks = peak; lanina_2020_23 (launch 2020-08, three peaks / winters,
#       per_winter); malformed multi-winter entries are refused
#   [2] lanina_2020_23's crossing launch on the first peak (ONI cache) and on a
#       synthetic series
#   [3] event_sub_windows(): consecutive per-winter sub-windows, the last ending
#       at the whole-event window; two peaks in one season refused
#   [4] sim_io.R::event_sub_windows_table(): NULL for a schema-2 event; country
#       or group bounds of a schema-3 event
#   [5] .cell_frame() bounds; cells without bounds cut exactly as before
#   [6] 09 cells on a synthetic event: event_attrib_cells(cells = NULL) equals
#       the frozen v2 function; a sub-window cell equals the whole-event cell of
#       the frames cut to it; event_cells_all() without sub-windows is
#       event_cells(); the af_posterior stash of the whole cell is unchanged
#   [7] observed totals per sub-window (event_observed_adm1(country_first =),
#       frozen equivalence at NULL; sub-windows partition the whole window)
#   [8] event_sub_window_keep() (report_event_af_model's cut)
#
# Run: Rscript analysis/tests/test_event_sub_windows.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(yaml); library(lubridate) })
source("analysis/utils/make_event_attribution_configs.R")
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")
source("analysis/utils/sim_io.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
cat("# event sub-windows\n")

# --- frozen copies of results_derive.R at 3eb4be4 (before the sub-windows) ----
old_event_attrib_cells <- function(inp, ctx, sets) {
  launch <- inp$attrib_start
  ec <- event_cells(ctx, launch, "event")   # one window (the whole run) x the explosive modes
  cells <- list()
  for (i in seq_len(nrow(ec))) {
    cl <- ec[i, ]
    m <- attrib_metrics_tables(inp$adm1_month, inp$country_month, inp$adm1_pop, ctx$GROUP_SHORT_RUN, sets,
                               launch, ctx$PR_CASE_THRESHOLD, ctx$PR_ALPHA, ctx$SN_LEVELS, cap = cl$cap,
                               n_boot = .ctx_boot_b(ctx))
    x <- attrib_exposure_tables(inp$adm1_month, inp$country_month, inp$adm1_pop, ctx$GROUP_SHORT_RUN, sets,
                                launch, ctx$ATTRIB_MAGNITUDES, ctx$ROBUST_EXCEED, ctx$ROBUST_AGREE,
                                cap = cl$cap, n_boot = .ctx_boot_b(ctx))
    tag <- function(d) .cell_tag(d, cl, window = FALSE)
    cells[[as.character(cl$explosive_guard)]] <- c(lapply(m[c("af_adm1", "af_scope", "pr_adm1", "pr_scope", "pr_pop", "sn_adm1", "sn_scope",
                                "sn_pop", "explosive")], tag),
                     list(exposure = tag(x$exposure), exposure_adm1 = tag(x$adm1), universe = m$universe))
  }
  tabs <- setdiff(names(cells[[1L]]), "universe")
  out <- lapply(stats::setNames(tabs, tabs), function(t) dplyr::bind_rows(lapply(cells, `[[`, t)))
  out$explosive <- dplyr::filter(out$explosive, explosive_guard)
  c(out, list(universe = cells[[1L]]$universe))
}

old_event_observed_adm1 <- function(obs, gids, country_last, launch) {
  launch <- as.Date(launch)
  g <- dplyr::distinct(gids, gid, country) %>%
    dplyr::mutate(gid = as.character(gid), adm1_unit = gidAdm1(gid, "event_observed_adm1"),
                  last = as.Date(unname(country_last[country])))
  if (anyNA(g$last))
    stop("[09] event_observed_adm1: no simulated month for country(ies) ",
         paste(unique(g$country[is.na(g$last)]), collapse = ", "), call. = FALSE)
  g$n_months <- vapply(seq_len(nrow(g)), function(i)
    if (g$last[i] < launch) 0L else length(seq(launch, g$last[i], by = "month")), integer(1))
  cells <- g %>% dplyr::group_by(adm1_unit) %>% dplyr::summarise(n_cells = sum(n_months), .groups = "drop")
  ob <- obs %>% dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::inner_join(dplyr::select(g, gid, adm1_unit, last), by = "gid") %>%
    dplyr::filter(t_start >= launch, t_start <= last) %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(obs_onsets = sum(obs_Y, na.rm = TRUE),
                     obs_cases = sum(cases_month, na.rm = TRUE),
                     n_cells_onset = sum(!is.na(obs_Y)),
                     n_cells_cases = sum(!is.na(cases_month)), .groups = "drop")
  dplyr::left_join(cells, ob, by = "adm1_unit") %>%
    dplyr::mutate(dplyr::across(c(obs_onsets, obs_cases, n_cells_onset, n_cells_cases),
                                ~ dplyr::coalesce(as.numeric(.x), 0)))
}

old_cell_frame <- function(df, cell) {
  out <- dplyr::filter(df, month >= cell$start[[1]])
  ge <- cell$group_ends[[1]]
  if (!is.null(ge)) {
    e <- unname(as.Date(ge)[match(out$enso_group, names(ge))])
    keep <- !is.na(e) & out$month <= e
    out <- out[keep, , drop = FALSE]
    out$group_horizon_end <- e[keep]
  }
  out
}


REG2 <- "analysis/configs/enso_events_v2.yml"

# --- [1] registry ------------------------------------------------------------
cat("\n[1] enso_events_v2.yml\n")
ids2 <- names(yaml::read_yaml(REG2)$events)
ids1 <- setdiff(ids2, "lanina_2020_23")
expect(length(ids1) == 5L && identical(tail(ids2, 1L), "lanina_2020_23"),
       "five single-winter events, then lanina_2020_23")
e1 <- load_enso_events(REG2, ids1)
expect(all(vapply(e1, function(e) identical(e$sub_windows, "none") && identical(e$peaks, e$peak) &&
                    identical(e$peak_season_years, e$peak_season_year), logical(1))),
       "single-peak events: peaks = peak, peak_season_years = peak_season_year, sub_windows none")
ln <- load_enso_events(REG2, "lanina_2020_23")$lanina_2020_23
expect(ln$launch == as.Date("2020-08-01") && identical(ln$launch_rule, "crossing") &&
         identical(ln$peaks, as.Date(c("2020-11-01", "2021-11-01", "2022-10-01"))) &&
         identical(ln$peak_season_years, 2020:2022) && identical(ln$sub_windows, "per_winter"),
       "lanina_2020_23: launch 2020-08 (crossing), peaks 2020-11 / 2021-11 / 2022-10, winters 2020-22, per_winter")
expect(ln$peak == as.Date("2022-10-01") && identical(ln$peak_season_year, 2020L),
       "peak = the last peak (whole window), peak_season_year = the first winter")
expect(identical(event_sub_window_ids(ln), c("winter_2020", "winter_2021", "winter_2022")), "sub-window ids")
# Event-level constraint scope: lanina_2020_23 constrains only its
# own winters; every other event has no key and the deploy's rule.
expect(identical(ln$other_winters, "any"), "lanina_2020_23: other_winters any (its own winters only)")
expect(all(vapply(e1, function(e) is.null(e$other_winters), logical(1))), "v1 events: no other_winters element")
cf_dep <- list(match = "library", n_analogs = 100L, winters = "all_read", other_winters = "neutral")
expect(identical(event_counterfactual_cfg(cf_dep, e1$elnino_2015), ensoCounterfactualCfg(cf_dep)),
       "event without the key: the deploy's matcher settings, unchanged")
cf_ln <- event_counterfactual_cfg(cf_dep, ln)
expect(identical(cf_ln$other_winters, "any") &&
         identical(cf_ln[setdiff(names(cf_ln), "other_winters")],
                   ensoCounterfactualCfg(cf_dep)[setdiff(names(cf_ln), "other_winters")]),
       "lanina_2020_23: only other_winters differs from the deploy's settings")
expect(identical(event_counterfactual_cfg(NULL, e1$elnino_2015), ensoCounterfactualCfg(NULL)),
       "no deploy block: the defaults")
bad <- function(entry) {
  f <- tempfile(fileext = ".yml"); yaml::write_yaml(list(events = list(x = entry)), f)
  !is.null(errmsg(load_enso_events(f, "x")))
}
base_e <- list(phase = "lanina", peaks = list("2020-11-01", "2021-11-01"), peak_season_years = list(2020L, 2021L),
               launch = "2020-08-01", sub_windows = "per_winter")
# Field replacement (modifyList() would merge the unnamed peak lists by position).
set_ <- function(e, ...) { v <- list(...); for (k in names(v)) e[[k]] <- v[[k]]; e }
expect(!bad(base_e), "a well-formed two-winter entry loads")
expect(bad(set_(base_e, peak = "2020-11-01")), "peak and peaks together: refused")
expect(bad(set_(base_e, peak_season_years = list(2020L))), "peaks / winters length mismatch: refused")
expect(bad(set_(base_e, peaks = list("2021-11-01", "2020-11-01"))), "non-increasing peaks: refused")
expect(bad(set_(base_e, peak_season_years = list(2021L, 2020L))), "non-increasing winters: refused")
expect(bad(set_(base_e, sub_windows = "per_month")), "unknown sub_windows: refused")
expect(bad(list(phase = "lanina", peak = "2020-11-01", peak_season_year = 2020L, launch = "2020-08-01",
                sub_windows = "per_winter")), "per_winter with one peak: refused")
expect(bad(set_(base_e, launch = "2020-12-01")), "launch not before the first peak: refused")
expect(bad(set_(base_e, peak_season_year = 2020L)), "peaks with a scalar peak_season_year: refused")
expect(bad(set_(base_e, other_winters = "none")), "unknown other_winters: refused")
expect(!bad(set_(base_e, other_winters = "neutral")), "other_winters neutral: loads")

# --- [2] crossing launch ------------------------------------------------------
cat("\n[2] crossing launch on the first peak\n")
if (file.exists("generated_data/enso_oni_cache.rds")) {
  oni <- as.data.frame(fetchEnsoIndex("oni", cache_dir = "generated_data/", refresh = FALSE))
  oni$date <- as.Date(oni$date)
  expect(ensoEventLaunch(oni, ln$peaks[1L], "lanina") == ln$launch,
         "lanina_2020_23: ensoEventLaunch() on the first peak = 2020-08")
} else cat("  (ONI cache check skipped: no generated_data/enso_oni_cache.rds)\n")
syn <- data.frame(date = seq(as.Date("2020-05-01"), by = "month", length.out = 9),
                  enso_index = c(0.04, -0.18, -0.29, -0.43, -0.77, -1.00, -1.11, -1.06, -0.99))
expect(ensoEventLaunch(syn, "2020-11-01", "lanina") == as.Date("2020-08-01"),
       "the registry's 2020 ONI values: crossing at 2020-09, launch 2020-08")

# --- [3] sub-windows ------------------------------------------------------------
cat("\n[3] event_sub_windows()\n")
sw <- event_sub_windows("2020-08-01", ln$peaks, 7L, c(occ = 0L, dur = -1L, mon = 2L))
expect(identical(sw$start, as.Date(c("2020-08-01", "2021-10-01", "2022-10-01"))) &&
         identical(sw$end, as.Date(c("2021-09-01", "2022-09-01", "2023-09-01"))),
       "season end July, max lag 2: 2021-07+2, 2022-07+2, 2023-07+2; consecutive from the launch")
w_whole <- event_group_window("2020-08-01", ln$peak, 7L, c(0L, -1L, 2L))
expect(sw$end[3] == w_whole$horizon_end, "the last sub-window ends at the whole-event window end")
sw_w <- event_sub_windows("2020-08-01", ln$peaks, 11L, c(-6L, -4L, -8L))
expect(identical(sw_w$end, as.Date(c("2020-11-01", "2021-11-01", "2022-11-01"))) &&
         sw_w$start[2] == as.Date("2020-12-01"),
       "season end November, all-lead group: ends on each season end, next starts the month after")
expect(!is.null(errmsg(event_sub_windows("2020-08-01", as.Date(c("2020-09-01", "2020-12-01")), 3L, 0L))),
       "two peaks in one cholera season (both end 2021-03): refused")

# --- [4] manifest table -----------------------------------------------------------
cat("\n[4] event_sub_windows_table()\n")
em <- list(schema = 3L, kind = "event_attribution", analog_set = "neutral", events = list(
  e2 = list(launch = "2014-09-01"),
  e3 = list(launch = "2020-08-01", sub_windows = list(
    list(id = "winter_2020", peak = "2020-11-01", peak_season_year = 2020L,
         groups = list(enso_east = list(start = "2020-08-01", end = "2021-09-01")),
         countries = list(KEN = list(start = "2020-08-01", end = "2021-06-01"),
                          UGA = list(start = "2020-08-01", end = "2021-09-01"))),
    list(id = "winter_2021", peak = "2021-11-01", peak_season_year = 2021L,
         groups = list(enso_east = list(start = "2021-10-01", end = "2022-09-01")),
         countries = list(KEN = list(start = "2021-07-01", end = "2022-06-01"),
                          UGA = list(start = "2021-10-01", end = "2022-09-01")))))))
expect(is.null(event_sub_windows_table(em, "e2")), "an event without sub-windows: NULL")
tb <- event_sub_windows_table(em, "e3")
expect(nrow(tb) == 4L && all(tb$level == "country") && identical(unique(tb$sub_window), c("winter_2020", "winter_2021")) &&
         tb$end[tb$key == "KEN" & tb$k == 2L] == as.Date("2022-06-01"),
       "per-country bounds when the manifest has them")
em_g <- em; for (k in 1:2) em_g$events$e3$sub_windows[[k]]$countries <- NULL
tg <- event_sub_windows_table(em_g, "e3")
expect(nrow(tg) == 2L && all(tg$level == "group") && identical(tg$key, c("enso_east", "enso_east")),
       "group bounds otherwise")

# --- [5] .cell_frame bounds ---------------------------------------------------------
cat("\n[5] .cell_frame()\n")
mths <- seq(as.Date("2020-08-01"), as.Date("2022-09-01"), by = "month")
cm <- expand.grid(month = mths, country = c("KEN", "UGA"), stringsAsFactors = FALSE) %>%
  mutate(enso_group = "enso_east", group_horizon_end = as.Date("2022-09-01"), n_cases = 1) %>% as_tibble()
c_plain <- tibble(window = "x", start = as.Date("2021-01-01"), group_ends = list(c(enso_east = as.Date("2021-12-01"))))
expect(identical(.cell_frame(cm, c_plain), old_cell_frame(cm, c_plain)), "a cell without bounds: cut as before")
c_none <- tibble(window = "x", start = as.Date("2021-01-01"), group_ends = list(NULL))
expect(identical(.cell_frame(cm, c_none), old_cell_frame(cm, c_none)), "no group ends, no bounds: as before")
bw <- tb[tb$k == 2L, ]
c_sw <- tibble(window = "winter_2021", start = min(bw$start), group_ends = list(NULL),
               bounds = list(select(bw, level, key, start, end)))
f2 <- .cell_frame(cm, c_sw)
expect(min(f2$month[f2$country == "KEN"]) == as.Date("2021-07-01") && max(f2$month[f2$country == "KEN"]) == as.Date("2022-06-01") &&
         min(f2$month[f2$country == "UGA"]) == as.Date("2021-10-01") && all(f2$group_horizon_end[f2$country == "UGA"] == as.Date("2022-09-01")),
       "country bounds: each country over its own sub-window, group_horizon_end = its end")
am5 <- tibble(adm1_unit = c("KEN_1", "UGA_1"), month = as.Date("2021-08-01"), enso_group = "enso_east")
expect(identical(.cell_frame(am5, c_sw)$adm1_unit, "KEN_1"), "an ADM1 frame keys the country by its gid prefix")

# --- [6] 09 cells ----------------------------------------------------------------
cat("\n[6] event_attrib_cells() / event_cells_all() on a synthetic event\n")
a0 <- as.Date("2020-08-01")
m6 <- seq(a0, by = "month", length.out = 4)
am <- expand.grid(replicate_id = 1:5, adm1_unit = c("KEN_1", "KEN_2", "ETH_1"), month = m6,
                  run = c("forecast", "counterfactual"), stringsAsFactors = FALSE) %>%
  mutate(analog_set = ifelse(run == "forecast", NA_character_, "neutral"),
         enso_group = ifelse(substr(adm1_unit, 1, 3) == "KEN", "enso_east", "enso_horn"),
         group_horizon_end = max(m6)) %>% as_tibble()
set.seed(86)
am$cases <- round(stats::rexp(nrow(am), 1 / 20) * ifelse(am$run == "forecast", 1.4, 1))
am$outbreaks <- as.integer(am$cases > 15)
mkcm <- function(d) d %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
  group_by(replicate_id, month, country, run, analog_set, enso_group, group_horizon_end) %>%
  summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
inp <- list(adm1_month = am, country_month = mkcm(am), attrib_start = a0,
            adm1_pop = tibble(adm1_unit = c("KEN_1", "KEN_2", "ETH_1"), pop = c(100, 40, 60)))
ctx <- new.env()
ctx$GROUP_SHORT_RUN <- c(enso_east = "East", enso_horn = "Horn")
ctx$PR_CASE_THRESHOLD <- 20; ctx$PR_ALPHA <- 0.05; ctx$SN_LEVELS <- c(1, 2)
ctx$ATTRIB_MAGNITUDES <- c(0.25, 0.5, 1.0); ctx$ROBUST_EXCEED <- 0.66; ctx$ROBUST_AGREE <- 0.80
ctx$EXPLOSIVE_CAP <- 150; ctx$EXPLOSIVE_MODES <- c("drop", "none")
whole <- suppressMessages(event_attrib_cells(inp, ctx, "neutral"))
frozen <- suppressMessages(old_event_attrib_cells(inp, ctx, "neutral"))
expect(identical(whole, frozen), "cells = NULL: identical to the frozen v2 event_attrib_cells()")
# Group-level sub-windows: east [08, 09] | [10, 11]; horn [08, 10] | [11, 11].
swt6 <- tibble(sub_window = rep(c("winter_a", "winter_b"), each = 2), k = rep(1:2, each = 2),
               peak = as.Date(NA), peak_season_year = NA_integer_, level = "group",
               key = rep(c("enso_east", "enso_horn"), 2),
               start = as.Date(c("2020-08-01", "2020-08-01", "2020-10-01", "2020-11-01")),
               end = as.Date(c("2020-09-01", "2020-10-01", "2020-11-01", "2020-11-01")))
cut_to <- function(d, r) {
  i <- match(d$enso_group, r$key)
  out <- d[d$month >= r$start[i] & d$month <= r$end[i], , drop = FALSE]
  out$group_horizon_end <- r$end[match(out$enso_group, r$key)]
  out
}
for (kk in 1:2) {
  r <- swt6[swt6$k == kk, ]
  sc <- event_sub_window_cells(ctx, a0, r)
  got <- suppressMessages(event_attrib_cells(inp, ctx, "neutral", sc))
  am_k <- cut_to(am, r)
  ref <- suppressMessages(event_attrib_cells(list(adm1_month = am_k, country_month = mkcm(am_k),
                                                  attrib_start = min(r$start), adm1_pop = inp$adm1_pop), ctx, "neutral"))
  expect(identical(got[names(got) != "universe"], ref[names(ref) != "universe"]) &&
           identical(got$universe, ref$universe),
         sprintf("sub-window %s cell == the whole-event cell of the frames cut to it", r$sub_window[1]))
}
expect(identical(event_cells_all(ctx, a0, "lanina_2020_23"), event_cells(ctx, a0, "lanina_2020_23")),
       "event_cells_all() without sub-windows is event_cells()")
ca <- event_cells_all(ctx, a0, "lanina_2020_23", swt6)
expect(nrow(ca) == 6L && identical(unique(ca$window), c("lanina_2020_23", "winter_a", "winter_b")) &&
         is.null(ca$bounds[[1]]) && nrow(ca$bounds[[3]]) == 2L,
       "with sub-windows: whole event + one cell per sub-window, x the explosive modes")
s_whole <- af_posterior_stash_build(c(inp, list(groups = c("enso_east", "enso_horn"))), event_cells(ctx, a0, "ev"))
s_all <- af_posterior_stash_build(c(inp, list(groups = c("enso_east", "enso_horn"))),
                                  event_cells_all(ctx, a0, "ev", swt6))
expect(identical(dplyr::filter(s_all$reps, cell == "ev"), s_whole$reps) &&
         identical(dplyr::filter(s_all$explosive, cell == "ev"), s_whole$explosive) &&
         setequal(unique(s_all$reps$cell), c("ev", "winter_a", "winter_b")),
       "af_posterior stash: the whole-event totals unchanged, sub-window totals added")

# --- [7] observed per sub-window --------------------------------------------------
cat("\n[7] observed totals per sub-window\n")
obs <- tibble(gid = c("KEN_1_1", "KEN_1_1", "KEN_1_2", "ETH_1_1", "ETH_1_1"),
              t_start = as.Date(c("2020-08-01", "2020-10-01", "2020-11-01", "2020-09-01", "2020-11-01")),
              obs_Y = c(1L, 1L, 0L, 1L, 1L), cases_month = c(10, 20, NA, 5, 7))
gids <- tibble(gid = c("KEN_1_1", "KEN_1_2", "ETH_1_1"), country = c("KEN", "KEN", "ETH"))
last <- c(KEN = as.Date("2020-11-01"), ETH = as.Date("2020-11-01"))
expect(identical(event_observed_adm1(obs, gids, last, a0), old_event_observed_adm1(obs, gids, last, a0)),
       "event_observed_adm1(country_first = NULL) identical to the frozen v2 function")
view <- list(countries = list(KEN = list(group = "enso_east"), ETH = list(group = "enso_horn")))
univ <- tibble(adm1_gid = c("KEN_1", "ETH_1"), country = c("KEN", "ETH"), enso_group = c("enso_east", "enso_horn"),
               region = c("East", "Horn"), pop = c(1, 1), pop_missing = FALSE)
osw <- event_sub_window_observed(obs, gids, last, swt6, view, univ, "ev")
whole_o <- event_observed_adm1(obs, gids, last, a0)
sum_o <- osw$winter_a$adm1 %>% bind_rows(osw$winter_b$adm1) %>% group_by(adm1_unit) %>%
  summarise(across(c(n_cells, obs_onsets, obs_cases, n_cells_onset, n_cells_cases), sum), .groups = "drop")
expect(isTRUE(all.equal(as.data.frame(sum_o), as.data.frame(arrange(whole_o, adm1_unit)), check.attributes = FALSE)),
       "the sub-windows partition the whole window: their observed totals and cells sum to the whole event's")
ke <- filter(osw$winter_a$adm1, adm1_unit == "KEN_1")
expect(ke$obs_onsets == 1 && ke$obs_cases == 10 && ke$n_cells == 4,
       "east winter_a [2020-08, 2020-09]: one onset, 10 cases, 2 gids x 2 months")
last_short <- c(KEN = as.Date("2020-10-01"), ETH = as.Date("2020-11-01"))
os2 <- event_sub_window_observed(obs, gids, last_short, swt6, view, univ, "ev")
expect(filter(os2$winter_b$adm1, adm1_unit == "KEN_1")$n_cells == 2,
       "a sub-window is cut at the country's last simulated month")

# --- [8] analytic AF cut -----------------------------------------------------------
cat("\n[8] event_sub_window_keep()\n")
kp <- event_sub_window_keep(c("KEN", "UGA", "KEN", "TZA"), as.Date(c("2021-07-01", "2021-07-01", "2022-07-01", "2021-12-01")),
                            tb[tb$k == 2L, ])
expect(identical(kp, c(TRUE, FALSE, FALSE, FALSE)), "country bounds: KEN 2021-07 in, UGA before its start, no bounds out")
kg <- event_sub_window_keep(c("KEN", "UGA"), as.Date(c("2021-09-01", "2021-10-01")), tg[tg$k == 2L, ], "enso_east")
expect(identical(kg, c(FALSE, TRUE)), "group bounds apply to every country of the group")

# --- [9] 10's reading of the sub-windows ---------------------------------------------
cat("\n[9] figures_supp.R::event_sub_window_view()\n")
suppressPackageStartupMessages(library(ggplot2))
source("analysis/utils/figures_supp.R")
ea9 <- list(events = tibble(event = c("elnino_2015", "lanina_2020_23"), label = c("El Nino 2015-16", "La Nina 2020-23"),
                            launch = as.Date(c("2014-09-01", "2020-08-01"))),
            af_scope = tibble(event = "elnino_2015", factual_source = "simulation", explosive_guard = TRUE,
                              level = "continent", scope = "Africa", region = NA_character_, metric = "cases",
                              af_of_means = 0.2, af_lo = 0.1, af_hi = 0.3))
eam9 <- list(scope = tibble(event = "elnino_2015", level = "continent", scope = "Africa", region = NA_character_,
                            metric = "cases", af_mean = 0.1, af_lo = 0, af_hi = 0.2))
expect(is.null(event_sub_window_view(ea9, eam9)), "no sub_windows: NULL (10 draws nothing more)")
ea9$sub_windows <- list(
  windows = tibble(event = "lanina_2020_23", sub_window = c("winter_2020", "winter_2021"), k = 1:2,
                   peak = as.Date(c("2020-11-01", "2021-11-01")), peak_season_year = 2020:2021,
                   level = "group", key = "enso_east", start = as.Date(c("2020-08-01", "2021-10-01")),
                   end = as.Date(c("2021-09-01", "2022-09-01"))),
  af_scope = tibble(event = "lanina_2020_23", sub_window = c("winter_2020", "winter_2021"),
                    factual_source = "simulation", explosive_guard = TRUE, level = "continent", scope = "Africa",
                    region = NA_character_, metric = "cases", af_of_means = c(0.3, -0.1), af_lo = 0, af_hi = 0.4))
sv <- event_sub_window_view(ea9, eam9)
expect(identical(sv$ea$events$event, c("lanina_2020_23:winter_2020", "lanina_2020_23:winter_2021")) &&
         identical(sv$ea$events$label, c("La Nina 2020-23, DJF 2020/21", "La Nina 2020-23, DJF 2021/22")) &&
         identical(unique(sv$ea$af_scope$event), sv$ea$events$event) && nrow(sv$eam$scope) == 0L,
       "one pseudo-event per (event, winter); no analytic sub-windows -> an empty analytic table")
p9 <- build_figS_event_attrib(sv$ea, sv$eam, "cases", "Africa", ggplot2::theme_minimal())
expect(inherits(p9, "ggplot") && nrow(ggplot2::layer_data(p9, 3L)) == 2L,
       "build_figS_event_attrib() draws the sub-window panels (two points)")

cat(sprintf("\n[test_event_sub_windows] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
