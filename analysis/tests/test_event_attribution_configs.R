# test_event_attribution_configs.R — the past-ENSO-event attribution generator
# (analysis/utils/make_event_attribution_configs.R) and its registry
# (analysis/configs/enso_events_v2.yml).
#
# HARD, quick-safe: pure helpers, no fit, no simulation.
#   1. every registry event loads and validates (annual relaunch = +12 months)
#   2. the crossing-rule launches equal utils/enso.R::ensoEventLaunch() on the
#      cached ONI (skipped when generated_data/enso_oni_cache.rds is absent)
#   3. event_group_window(): season containing the peak, + max positive lag
#   4. event_lenssen_lead_months() covers every panel month + the leads
#   5. event_lenssen_status(): absent / horizon / length / season-year / ok
#   6. the shipped deploy blocks name registry events; patch_hindcast_config()
#      is reachable at top level (the factual arm's builder)
#   7. event_param_draw_mode(): '' = unchanged; cycle needs max_draws coprime
#      to the smoothing's replicate count under uniform handoff (shuffled: any)
#   7b. event_n_replicates(): 0 / NULL = the template's count; N sets both
#      arms; negative or NA refused; CLI default 0
#   8. --season_end options: defaults (by_group, 6 months), no partial match on
#      season_end_by_group, event_group_window() = event_window_end()
#      (the by_country_from_panel generator run: test_deployed_lag_horizon.R [9])
#
# Run: Rscript analysis/tests/test_event_attribution_configs.R

suppressPackageStartupMessages({ library(yaml); library(lubridate) })
source("analysis/utils/make_event_attribution_configs.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

REG <- "analysis/configs/enso_events_v2.yml"
cat("# event attribution configs\n")

# --- 1. registry -------------------------------------------------------------
cat("\n[1] registry\n")
ids <- setdiff(names(yaml::read_yaml(REG)$events), "lanina_2020_23")   # the single-winter events
ev <- load_enso_events(REG, ids)
expect(length(ev) == 5L, sprintf("five events load (%s)", paste(ids, collapse = ", ")))
expect(all(vapply(ev, function(e) e$launch < e$peak, logical(1))), "every launch precedes its peak")
expect(ev$lanina_2022$launch == ev$lanina_2021$launch %m+% months(12L),
       "lanina_2022 relaunches 12 months after lanina_2021")
bad_reg <- tempfile(fileext = ".yml")
yaml::write_yaml(list(events = list(x = list(phase = "lanina", peak = "2022-10-01", launch = "2022-07-01",
                                             peak_season_year = 2022L, launch_rule = "annual",
                                             relaunch_of = "y"),
                                    y = list(phase = "lanina", peak = "2021-11-01", launch = "2021-08-01",
                                             peak_season_year = 2021L))), bad_reg)
expect(inherits(try(load_enso_events(bad_reg, "x"), silent = TRUE), "try-error"),
       "an annual launch that is not relaunch_of + 12 months is refused")
expect(inherits(try(load_enso_events(REG, "elnino_1997"), silent = TRUE), "try-error"),
       "an unknown event id is refused")
expect(inherits(try(load_enso_events(REG, ""), silent = TRUE), "try-error"),
       "an empty selection is refused")

# --- 2. crossing rule on the cached ONI -------------------------------------
cat("\n[2] crossing-rule launches vs the ONI cache\n")
oni_cache <- "generated_data/enso_oni_cache.rds"
if (file.exists(oni_cache)) {
  oni <- as.data.frame(fetchEnsoIndex("oni", cache_dir = "generated_data/", refresh = FALSE))
  oni$date <- as.Date(oni$date)
  for (e in ev) if (identical(e$launch_rule, "crossing")) {
    got <- ensoEventLaunch(oni, e$peak, e$phase)
    expect(got == e$launch, sprintf("%s: registry launch %s == ensoEventLaunch() %s",
                                    e$id, e$launch, got))
  }
  # The annual event: the crossing rule names the earlier year's launch.
  expect(ensoEventLaunch(oni, ev$lanina_2022$peak, "lanina") == ev$lanina_2021$launch,
         "lanina_2022: the crossing rule falls back to lanina_2021's launch (why it is annual)")
} else cat("  (skipped: no ONI cache)\n")
syn <- data.frame(date = seq(as.Date("2020-01-01"), by = "month", length.out = 8),
                  enso_index = c(0.1, 0.6, 0.2, 0.4, 0.55, 0.9, 1.2, 0.7))
expect(ensoEventLaunch(syn, "2020-07-01", "elnino") == as.Date("2020-04-01"),
       "ensoEventLaunch(): the month before the last unbroken run into the peak")
expect(inherits(try(ensoEventLaunch(syn, "2020-03-01", "elnino"), silent = TRUE), "try-error"),
       "ensoEventLaunch(): a peak outside the phase is refused")

# --- 3. windows --------------------------------------------------------------
cat("\n[3] event_group_window()\n")
w <- event_group_window("2014-09-01", "2015-12-01", 7L, c(occ = 0L, dur = -1L, mon = 2L))
expect(w$season_end == as.Date("2016-07-01") && w$horizon_end == as.Date("2016-09-01") && w$max_lag == 2L,
       "El Nino peak Dec 2015, season end July: 2016-07 + max lag 2 = 2016-09")
w <- event_group_window("2020-08-01", "2020-11-01", 11L, c(-6L, -4L, -8L))
expect(w$season_end == as.Date("2020-11-01") && w$horizon_end == as.Date("2020-11-01") && w$max_lag == 0L,
       "season end on the peak month counts; all-lead groups add nothing")
w <- event_group_window("2014-09-01", "2015-12-01", 9L, c(10L, 7L, 10L))
expect(w$horizon_end == as.Date("2017-07-01"), "Horn: 2016-09 + 10 = 2017-07")
expect(inherits(try(event_group_window("2020-08-01", "2020-11-01", 7L, -3L), silent = TRUE), "try-error") ==
         FALSE, "a season end after the launch is accepted")
expect(inherits(try(event_group_window("2021-08-01", "2021-01-01", 4L, 0L), silent = TRUE), "try-error"),
       "a window ending before the launch is refused")

# --- 4. Lenssen trajectory length ------------------------------------------
cat("\n[4] event_lenssen_lead_months()\n")
expect(event_lenssen_lead_months("2014-09-01", "2026-08-01", 11L) == 155L,
       "2014-09 .. 2026-08 is 144 panel months; + 11 leads = 155 (n_months 156 >= 144 + 11)")
expect(event_lenssen_lead_months("2026-08-01", "2026-08-01", 0L) == 1L, "one month, no lead")

# --- 5. cache status ---------------------------------------------------------
cat("\n[5] event_lenssen_status()\n")
root <- tempfile("lenssen_events_"); e <- ev$elnino_2015
st <- event_lenssen_status(e, root, 155L, "neutral")
expect(!st$ok && grepl("absent", st$reason), "absent cache")
d <- makeEventLenssenDir(root, e$id); dir.create(d, recursive = TRUE)
saveRDS(list(), makeLenssenAnalogFile(d))
meta <- list(horizon_start = "2014-09-01", lead_months = 155L, n_total_analogs = 550L,
             set_definition = list(peak_season_year = 2015L, counts = list(nina = 0L, neutral = 200L, elnino = 350L)))
mp <- file.path(d, "lenssen_nino34_analogs_meta.yml")
yaml::write_yaml(meta, mp)
st <- event_lenssen_status(e, root, 155L, "neutral")
expect(st$ok && st$n_set == 200L, "matching cache is ready (neutral = 200)")
expect(!event_lenssen_status(e, root, 155L, "nina")$ok, "an empty analog set is not ready")
expect(!event_lenssen_status(e, root, 160L, "neutral")$ok, "too-short trajectories are not ready")
yaml::write_yaml(utils::modifyList(meta, list(horizon_start = "2015-05-01")), mp)
expect(!event_lenssen_status(e, root, 155L, "neutral")$ok, "another horizon_start is not ready")
yaml::write_yaml(utils::modifyList(meta, list(set_definition = list(peak_season_year = 2014L))), mp)
expect(!event_lenssen_status(e, root, 155L, "neutral")$ok, "sets classified at another DJF are not ready")
cmd <- event_lenssen_command(e, root, 155L)
expect(grepl("--horizon_start 2014-09-01", cmd) && grepl("--set_peak_season_year 2015", cmd) &&
         grepl("--lead_months 155", cmd) && grepl("--oni_cache_dir generated_data/", cmd),
       "the build command pins launch, peak season year, length and the shared ONI cache")

# --- 6. deploy blocks + factual builder ------------------------------------
cat("\n[6] deploy blocks\n")
for (f in sort(list.files("analysis/configs/deploy", "\\.yml$", full.names = TRUE))) {
  ea <- yaml::read_yaml(f)$deploy$event_attribution
  evs <- trimws(strsplit(as.character(ea$events %||% ""), ",")[[1L]]); evs <- evs[nzchar(evs)]
  # Each deploy's own registry (enso_events_v2.yml by default).
  reg_f <- as.character(ea[["registry"]] %||% REG)
  if (length(evs))
    expect(file.exists(reg_f) && all(evs %in% names(yaml::read_yaml(reg_f)$events)) &&
             nzchar(ea$run_tag_base) && nzchar(ea$config_outdir) && nzchar(ea$manifest),
           sprintf("%s: events in its registry (%s) and run-specific outputs set", basename(f), basename(reg_f)))
}
expect(is.function(patch_hindcast_config), "patch_hindcast_config() is top level")

# --- 7. event_param_draw_mode() ---------------------------------------------
cat("\n[7] event_param_draw_mode()\n")
cfg0 <- list(simulator = list(max_draws = 199L, forecast = list(param_draw_mode = "random")))
expect(identical(event_param_draw_mode(cfg0, ""), cfg0), "'' leaves the config unchanged (the default)")
expect(identical(event_param_draw_mode(cfg0, NULL), cfg0), "NULL leaves the config unchanged")
expect(identical(event_param_draw_mode(cfg0, "cycle", 100L)$simulator$forecast$param_draw_mode, "cycle"),
       "cycle with max_draws 199 vs 100 smoothing replicates is set")
cfg200 <- cfg0; cfg200$simulator$max_draws <- 200L
e200 <- try(event_param_draw_mode(cfg200, "cycle", 100L), silent = TRUE)
expect(inherits(e200, "try-error") && grepl("e.g. 199", e200),
       "cycle with max_draws 200 vs 100 is refused, suggesting 199")
expect(inherits(try(event_param_draw_mode(cfg0, "cycle", NA), silent = TRUE), "try-error"),
       "cycle without the smoothing replicate count is refused")
cfg200s <- cfg200; cfg200s$simulator$forecast$handoff_weighting <- "shuffled"
expect(identical(event_param_draw_mode(cfg200s, "cycle", 100L)$simulator$forecast$param_draw_mode, "cycle") &&
         identical(event_param_draw_mode(cfg200s, "cycle", NA)$simulator$forecast$param_draw_mode, "cycle"),
       "shuffled handoff: cycle with max_draws 200 vs 100 is set (no coprime check, no count needed)")
cfg200u <- cfg200; cfg200u$simulator$forecast$handoff_weighting <- "uniform"
expect(inherits(try(event_param_draw_mode(cfg200u, "cycle", 100L), silent = TRUE), "try-error"),
       "explicit uniform handoff: max_draws 200 vs 100 still refused")
expect(inherits(try(event_param_draw_mode(cfg0, "sequential"), silent = TRUE), "try-error"),
       "an unknown mode is refused")
expect(identical(event_param_draw_mode(cfg0, "random")$simulator$forecast$param_draw_mode, "random"),
       "random is accepted without the replicate count")

# --- 7b. event_n_replicates() -------------------------------------------------
cat("\n[7b] event_n_replicates()\n")
cfg_n <- list(simulator = list(n_replicates = 9950L))
expect(identical(event_n_replicates(cfg_n, 0L), cfg_n), "0 leaves the template's n_replicates (the default)")
expect(identical(event_n_replicates(cfg_n, NULL), cfg_n), "NULL leaves the template's n_replicates")
expect(identical(event_n_replicates(cfg_n, 5000L)$simulator$n_replicates, 5000L), "N sets simulator.n_replicates")
expect(inherits(try(event_n_replicates(cfg_n, -1L), silent = TRUE), "try-error"), "a negative count is refused")
expect(inherits(try(event_n_replicates(cfg_n, NA), silent = TRUE), "try-error"), "NA is refused")
po7 <- optparse::parse_args(optparse::OptionParser(option_list = option_list), args = character(0))
expect(identical(po7$n_replicates, 0L) && identical(get_default_options()$deploy$event_attribution$n_replicates, 0L),
       "CLI --n_replicates default 0 = deploy.event_attribution.n_replicates default")

# --- 8. season-end options (per-country season ends) -------------
cat("\n[8] --season_end options\n")
po <- optparse::parse_args(optparse::OptionParser(option_list = option_list), args = character(0))
expect(identical(po$season_end, "by_group") && identical(po$season_min_window_months, 6L),
       "CLI defaults: --season_end by_group, --season_min_window_months 6")
expect(is.null(list(season_end_by_group = "enso_east=7")[["season_end"]]),
       "a deploy opt without season_end does not partial-match season_end_by_group ([[ ]] read)")
w <- event_group_window("2014-09-01", "2015-12-01", 7L, c(occ = 0L, dur = -1L, mon = 2L))
expect(identical(w, event_window_end(as.Date("2015-12-01"), 7L, c(0L, -1L, 2L))),
       "event_group_window() is season_calendar.R::event_window_end() + the launch check")
expect(grepl("by_country_from_panel", paste(readLines("analysis/utils/make_event_attribution_configs.R"),
                                            collapse = "\n")),
       "the generator documents by_country_from_panel")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
