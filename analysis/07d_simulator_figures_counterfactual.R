#!/usr/bin/env Rscript
# 07d_simulator_figures_counterfactual.R — forecast vs ENSO-counterfactual figures
#
# Compares the grouped-ENSO-mode spatPomp FACTUAL FORECAST against its
# ENSO-COUNTERFACTUAL companion (Lenssen piControl-analog ENSO), isolating the
# excess cholera risk attributable to the observed/forecast ENSO anomaly
# (FAR-style attribution). Sibling to 07d_simulator_figures.R (validation figures).
# The manifest/run-tag resolution, partition loader, replicate ribbons and ADM1
# paired / unpaired statistics live in utils/results_core.R, shared with
# 09b_postprocess_results.R.
#
# Run tags are resolved per ENSO teleconnection GROUP from the manifest written by
# analysis/utils/make_forecast_configs.R (analysis/configs/forecast_horizons.yml):
# each group's own forecast horizon (`horizon_end`), the common cutoff shared by all
# groups (`africa_min_horizon`), and the factual / counterfactual run tag read from
# that group's own generated config (`simulator.run_tag`).
#
# Design notes
# ------------
#   * Within each group the factual + counterfactual run share posterior_strategy:
#     paired + seed + n_replicates, so replicate_id is comparable 1:1 across arms.
#     The ADM1 ratio map colours the ratio of MEANS (stable at ADM1 granularity, where
#     per-replicate ratios are mostly 0/0 or Inf); significance (map border) is tested
#     on the paired per-replicate difference forecast_i - counterfactual_i.
#   * "Region" = ENSO teleconnection sub-group (the grouping driving the fits), not
#     AFRO_region. enso_central is in the manifest's `skipped_groups` and is dropped.
#   * CASES: the mon track's `cases_month` = one draw per ACTIVE month, the
#     manuscript estimand (needs runs with grouping.monthly_enable: true).
#   * No importance re-weighting of replicates: forward draws from the paired
#     posterior/particle ensemble handed off from the historical smoothing, all
#     equally weighted.
#
# Outputs (figures/sim/<manifest-stem>/counterfactual/ unless --out_dir):
#   fig1_africa.png        Africa-wide monthly cases + outbreaks, truncated to africa_min_horizon
#   fig2_region.png        same by ENSO group, each to its OWN horizon_end
#   fig6_country_outbreaks.png / fig7_country_cases.png   same by country
#   fig3_adm1_ratio_map_{common,full}.png          ADM1 cumulative-incidence ratio (paired test outline)
#   fig4_adm1_outbreak_ratio_map_{common,full}.png same for outbreak-months
#   fig5_adm1_prob_shift_map_{common,full}.png     unpaired nrep x nrep probability-of-shift
#   adm1_cumulative_{incidence,outbreak}_ratio_{common,full}.csv,
#   adm1_cumulative_incidence_prob_{common,full}.csv   the underlying ADM1 tables
# "common" = every ADM1 truncated to africa_min_horizon; "full" = each ADM1 to its own
# group's horizon_end. The derived-aggregate cache lives under
# <sim_root>/_figure_cache/, keyed on the manifest AND every input that changes it.
#
# Usage:
#   Rscript analysis/07d_simulator_figures_counterfactual.R
#   Rscript analysis/07d_simulator_figures_counterfactual.R --manifest analysis/configs/forecast_horizons.yml \
#     --ratio_start 2026-09-01 --ratio_end 2026-12-01 --prob_threshold 0.25
#
# --ratio_start "auto" = the latest observed ONI month (the ENSO handoff);
# --ratio_end "auto" = per-group horizon_end / africa_min_horizon from --manifest (an
# explicit value overrides BOTH windows uniformly).
# --horizon firm (default) = each country to its fit group's firm horizon (last CFSv2
# month + the group's min deployed lag); unit = to the horizon its simulation unit was
# run to (v2 manifests; the forecast's own end, past the CFSv2 span on the ENSO
# continuation). Group horizon = the latest of its countries', common = the earliest.
# unit writes to <figdir>/sim/<manifest-stem>/counterfactual_unit/ (firm figures kept).
# --arm_alignment stitch (default) | common_start: when the two arms launch at
# different months (manifest launch keys, utils/sim_io.R::manifest_arm_starts()),
# "stitch" starts the window at the earlier (counterfactual) launch -- it
# replaces --ratio_start auto -- and takes the factual series from the smoothing
# up to its own launch (--smoothing_parts_dir, default the parts/ beside the
# factual configs' smoothing_input), paired by replicate through 07b's
# handoff_map.csv; "common_start" keeps the first month both arms simulate.
# Equal launches: no-op either way.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(qs2)
  library(cowplot)
  library(sf)
  library(scales)
})

# results_core.R sources base/config/paths/provenance/enso/gid/figures_cb itself.
source("analysis/utils/results_core.R")
source("analysis/utils/population.R")       # build_location_frame() call-time closure
source("analysis/utils/composite_lp.R")     #   (population <-> composite_lp are mutual)
source("analysis/utils/anchor.R")
source("analysis/utils/simulator_core.R")   # build_location_frame()
source("analysis/utils/figures_io.R")

# ---------------------------------------------------------------------------
# Section 0: CLI
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--manifest", type = "character", default = "analysis/configs/forecast_horizons.yml",
              help = "Per-ENSO-group forecast-horizon manifest (make_forecast_configs.R) [default: %default]"),
  make_option("--sim_root", type = "character", default = "generated_data/sim_spatpomp"),
  make_option("--ratio_start", type = "character", default = "auto",
              help = "Start of the cumulative window; 'auto' = latest observed ONI month [default: %default]"),
  make_option("--ratio_end", type = "character", default = "auto",
              help = "End of the window (inclusive); overrides BOTH the common and full windows; 'auto' = from --manifest [default: %default]"),
  make_option("--horizon", type = "character", default = "firm",
              help = "firm | unit: censor each country at its firm horizon or at its simulation unit's horizon (v2 manifests) [default: %default]"),
  make_option("--series", type = "character", default = "latent",
              help = "latent | detected (needs simulator.filter.covariate_feed: detected) [default: %default]"),
  make_option("--out_dir", type = "character", default = NULL,
              help = "Figure directory [default: <figdir>/sim/<manifest-stem>/counterfactual]"),
  make_option("--figdir", type = "character", default = "figures", help = "Figure root [default: %default]"),
  make_option("--cache_dir", type = "character", default = NULL,
              help = "Aggregate cache dir [default: <sim_root>/_figure_cache]"),
  make_option("--ci_level", type = "double", default = 0.50,
              help = "CI level of the paired difference test [default: %default]"),
  make_option("--prob_threshold", type = "double", default = 0.50,
              help = paste("Fractional shift threshold of the probability-of-shift maps (fig5): P that",
                           "forecast cumulative cases exceed counterfactual by more than this fraction",
                           "('larger') or fall short by more ('smaller'), over ALL forecast x counterfactual",
                           "replicate pairs. 0.50 = +/-50%% [default: %default]")),
  make_option("--arm_alignment", type = "character", default = "stitch",
              help = "stitch | common_start: arms launching at different months (see header) [default: %default]"),
  make_option("--smoothing_parts_dir", type = "character", default = NULL,
              help = "Smoothing parts/ the factual series is stitched from [default: beside the factual configs' smoothing_input]"),
  make_option("--redo_cache", type = "logical", default = FALSE, help = "Rebuild the aggregate cache")
)
args <- parse_args(OptionParser(option_list = option_list))

SERIES         <- args$series
SIM_ROOT       <- args$sim_root
CI_LEVEL       <- args$ci_level
PROB_THRESHOLD <- args$prob_threshold
MANIFEST_PATH  <- args$manifest
HORIZON        <- match.arg(args$horizon, c("firm", "unit"))
OUT_DIR   <- args$out_dir %||% makeFigureDir(list(figdir = args$figdir), track = "sim",
                                             config_path = MANIFEST_PATH,
                                             subdir = if (HORIZON == "unit") "counterfactual_unit" else "counterfactual")
CACHE_DIR <- args$cache_dir %||% file.path(SIM_ROOT, "_figure_cache")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
save_fig <- figureSaver(OUT_DIR, width = 10, height = 9, dpi = 150)   # + bg = "white" per call

RATIO_START <- if (identical(args$ratio_start, "auto")) {
  d <- latest_oni_month()
  message(sprintf("[07d] --ratio_start auto -> latest observed ONI month: %s", d))
  d
} else as.Date(args$ratio_start)
RATIO_END_OVERRIDE <- if (identical(args$ratio_end, "auto")) NULL else as.Date(args$ratio_end)

# ---------------------------------------------------------------------------
# Manifest: group -> {run_tag, horizon_end}, plus africa_min_horizon
# ---------------------------------------------------------------------------
man    <- load_manifest(MANIFEST_PATH)
GROUPS <- names(man$groups)
if (HORIZON == "unit" && !length(man$units))
  stop("[07d] --horizon unit needs a v2 manifest (partition_by: registry); this one has no units.")
# ISO3 -> censoring horizon (firm or unit); NULL for v1 manifests (firm = group horizon_end)
COUNTRY_H <- if (length(man$units)) manifest_country_horizons(man, HORIZON) else NULL
GROUP_H <- if (HORIZON == "unit") {
  cg_h <- manifest_country_groups(man)
  vapply(GROUPS, function(g) max(COUNTRY_H[names(cg_h)[cg_h == g]]), character(1))
} else vapply(GROUPS, function(g) as.character(man$groups[[g]]$horizon_end), character(1))
group_horizon_end <- function(g) as.Date(GROUP_H[[g]])

# Per-arm launches (primary analog set, the only one 07d reads) and the stitch.
ARM_ALIGNMENT <- match.arg(args$arm_alignment, c("stitch", "common_start"))
ARM_STARTS <- manifest_arm_starts(man, sets = man$primary_analog_set %||% "all")
STITCH <- identical(ARM_ALIGNMENT, "stitch") && ARM_STARTS$counterfactual < ARM_STARTS$factual
if (STITCH && identical(args$ratio_start, "auto")) {
  RATIO_START <- ARM_STARTS$counterfactual
  message(sprintf("[07d] arms launch at %s (factual) / %s (counterfactual): window from %s, factual series stitched from the smoothing",
                  ARM_STARTS$factual, ARM_STARTS$counterfactual, RATIO_START))
}
STITCH <- STITCH && RATIO_START < ARM_STARTS$factual
SMOOTHING_PARTS_DIR <- if (STITCH) args$smoothing_parts_dir %||% smoothing_parts_dir_of(
  yaml::read_yaml(manifest_runs(man, "factual")[[1L]]$config)$simulator$filter$smoothing_input) else NULL
RATIO_NOTE <- if (attribution_ratio_label(ARM_STARTS$factual, ARM_STARTS$counterfactual) == "ENSO") "" else
  sprintf("\nratio = %s (factual %s)", attribution_ratio_label(ARM_STARTS$factual, ARM_STARTS$counterfactual),
          if (STITCH) sprintf("from the smoothing to %s", format(ARM_STARTS$factual, "%b %Y"))
          else sprintf("launched %s", format(ARM_STARTS$factual, "%b %Y")))
STITCH_VLINE <- if (STITCH) ARM_STARTS$factual else NULL

AFRICA_MIN_HORIZON <- RATIO_END_OVERRIDE %||%
  (if (HORIZON == "unit") as.Date(min(GROUP_H)) else as.Date(man$africa_min_horizon))
FULL_WINDOW_LABEL <- paste(
  sprintf("%s->%s", GROUPS,
          vapply(GROUPS, function(g) format(RATIO_END_OVERRIDE %||% group_horizon_end(g), "%b %Y"), character(1))),
  collapse = "; ")

SURV_PATH      <- "generated_data/outbreaks_survival_data_base-1990_2019.rds"
SIZE_DATA_PATH <- "generated_data/outbreak_size_analysis_data.rds"
HMM_PATH       <- "data/hmm_country_clusters.csv"

RUN_COLORS           <- c(forecast = "#2166ac", counterfactual = "#e08214")
LOW_SIGNAL_CASES     <- 2     # ADM1 is low-signal if forecast+counterfactual expected cases < this
LOW_SIGNAL_OUTBREAKS <- 1     # ... expected outbreak-months < this
RATIO_EPS            <- 0.25  # pseudocount so the ratio is finite when one arm is exactly 0
MAP_LIM              <- c(-2, 2)   # log2 ratio -> 0.25x .. 4x

message(sprintf("[07d] manifest: %s | groups: %s (skipped: %s)", MANIFEST_PATH,
                paste(GROUPS, collapse = ", "), paste(man$skipped_groups %||% character(0), collapse = ", ")))
message(sprintf("[07d] horizon: %s", HORIZON))
message(sprintf("[07d] cases = cases_month (%s) | window start %s | common end (africa_min_horizon) %s | full: %s",
                SERIES, RATIO_START, AFRICA_MIN_HORIZON, FULL_WINDOW_LABEL))

# ---------------------------------------------------------------------------
# Section 1: location frame -> country / ADM1-population crosswalks
# ---------------------------------------------------------------------------
message("[07d] building location frame (country / ADM1-pop crosswalks)...")
frame <- build_location_frame(
  survival_data          = readRDS(SURV_PATH),
  hmm_clusters_path      = HMM_PATH,
  outbreak_metadata_path = NULL,
  size_data_path         = SIZE_DATA_PATH)
country_region <- frame %>% dplyr::distinct(country, AFRO_region)
adm1_pop <- frame %>% dplyr::group_by(adm1_unit) %>%
  dplyr::summarise(pop_adm1 = sum(pop, na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------------
# Section 2: load every group's two runs (cached, keyed on every input that changes them)
# ---------------------------------------------------------------------------
# v2 manifests (partition_by: registry): runs are per partition unit; each unit
# is loaded once and its rows relabelled with each country's FIT group (the
# figures stay keyed by the fit grouping), every country censored at its own
# firm horizon (utils/sim_io.R).
UNIT_MODE <- length(man$units) > 0L
run_tags <- if (UNIT_MODE) {
  unlist(lapply(manifest_units(man), function(u)
    c(yaml::read_yaml(u$factual_config)$simulator$run_tag,
      yaml::read_yaml(u$counterfactual_config %||%
                        u$counterfactual_sets[[u$primary_analog_set %||% man$primary_analog_set]]$config)$simulator$run_tag)))
} else unlist(lapply(GROUPS, function(g)
  c(group_run_tag(man, g, "forecast"), group_run_tag(man, g, "counterfactual"))))
cache_key  <- substr(digest::digest(list(man = man, horizon = HORIZON, ratio_start = RATIO_START,
                                         ratio_end = RATIO_END_OVERRIDE, sim_root = SIM_ROOT,
                                         series = SERIES, run_tags = run_tags,
                                         arm_alignment = ARM_ALIGNMENT, arm_starts = ARM_STARTS,
                                         stitch = STITCH, smoothing_parts_dir = SMOOTHING_PARTS_DIR)), 1, 12)
cache_file <- file.path(CACHE_DIR, sprintf("counterfactual_%s_%s.qs",
                                           tools::file_path_sans_ext(basename(MANIFEST_PATH)), cache_key))

if (args$redo_cache || !file.exists(cache_file)) {
  message("[07d] building aggregate cache -> ", cache_file)
  country_month_list <- list(); adm1_month_list <- list()
  # Read once, shared by every group / unit (results_core.R::load_arm_run()).
  SEV <- if (STITCH) read_smoothing_parts_events(SMOOTHING_PARTS_DIR, RATIO_START, SERIES, tag = "[07d]") else NULL
  arm_starts_arg <- if (STITCH) ARM_STARTS else NULL
  if (UNIT_MODE) {
    cg   <- manifest_country_groups(man)
    firm <- COUNTRY_H   # --horizon firm: firm horizons; unit: the unit's
    for (un in names(manifest_units(man))) {
      u <- manifest_units(man)[[un]]
      tag_fcst <- yaml::read_yaml(u$factual_config)$simulator$run_tag
      tag_ctrf <- yaml::read_yaml(u$counterfactual_config %||%
        u$counterfactual_sets[[u$primary_analog_set %||% man$primary_analog_set]]$config)$simulator$run_tag
      ratio_end_u <- RATIO_END_OVERRIDE %||% as.Date(u$horizon_end)
      message(sprintf("[07d] unit %-10s forecast=%s counterfactual=%s horizon_end=%s (groups %s)",
                      un, tag_fcst, tag_ctrf, ratio_end_u, paste(u$groups, collapse = ",")))
      load_one <- function(tag, arm) load_arm_run(SIM_ROOT, tag, arm, RATIO_START, ratio_end_u,
                                                  arm_starts = arm_starts_arg, smoothing_events = SEV,
                                                  series = SERIES, tag = "[07d]",
                                                  horizon_by_country = firm)
      run_fcst <- load_one(tag_fcst, "forecast"); run_ctrf <- load_one(tag_ctrf, "counterfactual")
      relabel <- function(df, iso) df %>% dplyr::mutate(
        run = dplyr::case_when(run == tag_fcst ~ "forecast", run == tag_ctrf ~ "counterfactual", TRUE ~ run),
        enso_group = unname(cg[iso]),
        group_horizon_end = RATIO_END_OVERRIDE %||% as.Date(unname(GROUP_H[cg[iso]])))
      cm <- dplyr::bind_rows(run_fcst$country_month, run_ctrf$country_month)
      am <- dplyr::bind_rows(run_fcst$adm1_month, run_ctrf$adm1_month)
      country_month_list[[un]] <- relabel(cm, cm$country)
      adm1_month_list[[un]]    <- relabel(am, substr(am$adm1_unit, 1L, 3L))
    }
  }
  for (g in if (UNIT_MODE) character(0) else GROUPS) {
    tag_fcst <- group_run_tag(man, g, "forecast"); tag_ctrf <- group_run_tag(man, g, "counterfactual")
    ratio_end_g <- RATIO_END_OVERRIDE %||% group_horizon_end(g)
    message(sprintf("[07d] group %-10s forecast=%s counterfactual=%s horizon_end=%s",
                    g, tag_fcst, tag_ctrf, ratio_end_g))
    load_one <- function(tag, arm) load_arm_run(SIM_ROOT, tag, arm, RATIO_START, ratio_end_g,
                                                arm_starts = arm_starts_arg, smoothing_events = SEV,
                                                series = SERIES, tag = "[07d]")
    run_fcst <- load_one(tag_fcst, "forecast"); run_ctrf <- load_one(tag_ctrf, "counterfactual")
    relabel_run <- function(df) df %>% dplyr::mutate(run = dplyr::case_when(
      run == tag_fcst ~ "forecast", run == tag_ctrf ~ "counterfactual", TRUE ~ run))
    country_month_list[[g]] <- dplyr::bind_rows(run_fcst$country_month, run_ctrf$country_month) %>%
      relabel_run() %>% dplyr::mutate(enso_group = g, group_horizon_end = ratio_end_g)
    adm1_month_list[[g]] <- dplyr::bind_rows(run_fcst$adm1_month, run_ctrf$adm1_month) %>%
      relabel_run() %>% dplyr::mutate(enso_group = g, group_horizon_end = ratio_end_g)
  }
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(list(country = country_month_list, adm1 = adm1_month_list), file = cache_file)
} else {
  message("[07d] loading aggregate cache: ", cache_file)
  cache <- qs2::qs_read(file = cache_file)
  country_month_list <- cache$country; adm1_month_list <- cache$adm1
  rm(cache)
}

country_month_both <- dplyr::bind_rows(country_month_list) %>%
  dplyr::mutate(run = factor(run, levels = c("forecast", "counterfactual"))) %>%
  dplyr::left_join(country_region, by = "country")
adm1_month_both <- dplyr::bind_rows(adm1_month_list) %>%
  dplyr::mutate(run = factor(run, levels = c("forecast", "counterfactual")))
# The two arms may hand off at different months (e.g. the counterfactual at the
# Lenssen analog anchor, the forecast at the latest smoothed month + 1). Under
# --arm_alignment stitch the factual series already starts at
# the window start (smoothing before its launch), so this is a no-op. Otherwise
# (common_start) the paired ratio / probability maps cumulate from the first
# loaded month, so they are restricted to the months BOTH arms simulate, per
# group; otherwise the earlier-starting arm accumulates extra months. The ribbon
# figures keep each arm's full span.
.arm_start <- adm1_month_both %>%
  dplyr::group_by(enso_group, run) %>% dplyr::summarise(m0 = min(month), .groups = "drop") %>%
  dplyr::group_by(enso_group) %>% dplyr::summarise(common_start = max(m0), .groups = "drop")
if (any(.arm_start$common_start > RATIO_START)) {
  message(sprintf("[07d] ratio/probability maps start at the first month both arms simulate: %s",
                  paste(sprintf("%s=%s", .arm_start$enso_group, format(.arm_start$common_start)),
                        collapse = ", ")))
  adm1_month_both <- adm1_month_both %>%
    dplyr::inner_join(.arm_start, by = "enso_group") %>%
    dplyr::filter(month >= common_start) %>% dplyr::select(-common_start)
  RATIO_START <- min(.arm_start$common_start)   # map labels ("from <month>")
}
# adm1_unit -> enso_group / group_horizon_end (group provenance after cumulating over replicates)
adm1_group_horizon <- adm1_month_both %>% dplyr::distinct(adm1_unit, enso_group, group_horizon_end)

message(sprintf("[07d] countries loaded: %d across %d group(s)",
                dplyr::n_distinct(country_month_both$country), length(GROUPS)))

# ---------------------------------------------------------------------------
# Section 3: Africa-wide + by-group replicate aggregation. Africa totals are truncated to
# africa_min_horizon (past it the shortest-horizon group's run has stopped, which would
# show as an artificial drop); per-group series run to their own horizon_end.
# ---------------------------------------------------------------------------
africa_reps <- country_month_both %>% dplyr::filter(month <= AFRICA_MIN_HORIZON) %>%
  dplyr::group_by(run, replicate_id, month) %>%
  dplyr::summarise(n_cases = sum(n_cases), n_outbreaks = sum(n_outbreaks), .groups = "drop")
region_reps <- country_month_both %>%
  dplyr::group_by(run, replicate_id, month, enso_group) %>%
  dplyr::summarise(n_cases = sum(n_cases), n_outbreaks = sum(n_outbreaks), .groups = "drop")
africa_ribbon  <- summarise_reps(africa_reps, c("run", "month"))
region_ribbon  <- summarise_reps(region_reps, c("run", "month", "enso_group"))
country_ribbon <- summarise_reps(country_month_both, c("run", "month", "country"))

# ---------------------------------------------------------------------------
# Section 4: ribbon figures (forecast vs counterfactual overlay)
# ---------------------------------------------------------------------------
ribbon_plot2 <- function(ribbon_df, y_stem, y_lab, title = NULL, vline = NULL) {
  q <- function(s) paste0(y_stem, "_", s)
  p <- ggplot2::ggplot(ribbon_df, ggplot2::aes(x = month, fill = run, colour = run)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[q("q10")]], ymax = .data[[q("q90")]]), alpha = 0.18, colour = NA) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[q("q25")]], ymax = .data[[q("q75")]]), alpha = 0.35, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[q("q50")]]), linewidth = 0.9)
  if (!is.null(vline))   # africa_min_horizon (common window) [+ the factual launch when stitched]
    p <- p + ggplot2::geom_vline(xintercept = vline, linetype = "dashed", colour = "grey45", linewidth = 0.4)
  p +
    ggplot2::scale_fill_manual(values = RUN_COLORS, name = NULL) +
    ggplot2::scale_colour_manual(values = RUN_COLORS, name = NULL) +
    ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
    ggplot2::labs(x = NULL, y = y_lab, title = title) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
}
CASES_LAB <- "Cases / month"

message("[07d] fig1: Africa-wide forecast vs counterfactual...")
fig1 <- cowplot::plot_grid(
  ribbon_plot2(africa_ribbon, "n_outbreaks", "Outbreaks / month", "Africa — outbreaks per month", vline = c(AFRICA_MIN_HORIZON, STITCH_VLINE)),
  ribbon_plot2(africa_ribbon, "n_cases", CASES_LAB, "Africa — suspected cases per month", vline = c(AFRICA_MIN_HORIZON, STITCH_VLINE)),
  ncol = 1, align = "v")
cowplot::save_plot(file.path(OUT_DIR, "fig1_africa.png"), fig1, base_width = 10, base_height = 7, bg = "white")

message("[07d] fig2: by ENSO teleconnection group, each to its own horizon_end...")
groups <- sort(unique(region_ribbon$enso_group))
fig2 <- cowplot::plot_grid(plotlist = lapply(groups, function(g) {
  d <- dplyr::filter(region_ribbon, enso_group == g)
  cowplot::plot_grid(
    ribbon_plot2(d, "n_outbreaks", "Outbreaks / month",
                 sprintf("%s (horizon: %s)", g, format(group_horizon_end(g), "%b %Y")), vline = c(AFRICA_MIN_HORIZON, STITCH_VLINE)),
    ribbon_plot2(d, "n_cases", CASES_LAB, paste0(g, " — cases"), vline = c(AFRICA_MIN_HORIZON, STITCH_VLINE)),
    ncol = 1, align = "v")
}), ncol = length(groups))
cowplot::save_plot(file.path(OUT_DIR, "fig2_region.png"), fig2,
                   base_width = 4.5 * length(groups), base_height = 7, bg = "white")

n_countries <- dplyr::n_distinct(country_ribbon$country)
facet_ncol  <- min(5L, n_countries); facet_nrow <- ceiling(n_countries / facet_ncol)
country_fig <- function(y_stem, y_lab, title, stem) {
  p <- ribbon_plot2(country_ribbon, y_stem, y_lab, title, vline = c(AFRICA_MIN_HORIZON, STITCH_VLINE)) +
    ggplot2::facet_wrap(~ country, scales = "free_y", ncol = facet_ncol)
  save_fig(p, stem, width = 3.2 * facet_ncol, height = 2.6 * facet_nrow, bg = "white", limitsize = FALSE)
}
message("[07d] fig6/fig7: per country...")
country_fig("n_outbreaks", "Outbreaks / month", "Outbreaks per month by country (forecast vs counterfactual)", "fig6_country_outbreaks")
country_fig("n_cases", CASES_LAB, "Cases per month by country (forecast vs counterfactual)", "fig7_country_cases")

# ---------------------------------------------------------------------------
# Section 5-6: ADM1 ratio (paired) and probability-of-shift (unpaired) maps
# ---------------------------------------------------------------------------
geo     <- load_africa_sf()
afro_sf <- geo$afro_sf
adm1_sf <- geo$adm1_sf
if (is.null(adm1_sf)) stop("[07d] ADM1 shapefile failed to load — see load_africa_sf() warning above.")

#' One ADM1 diverging map (ratio and probability-of-shift variants): black outline where
#' `significant`, grey where the fill is NA (low signal).
plot_adm1_diverging_map <- function(df, fill_col, scale_args, fill_name, title, subtitle) {
  map_sf <- adm1_sf %>% dplyr::inner_join(df, by = c("GID_1" = "adm1_unit"))
  sig_sf <- dplyr::filter(map_sf, significant)
  message(sprintf("[07d] %s: polygons matched %d / %d, outlined significant %d",
                  title, nrow(map_sf), nrow(df), nrow(sig_sf)))
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = afro_sf, fill = "grey70", colour = "grey60", linewidth = 0.15) +
    ggplot2::geom_sf(data = map_sf, ggplot2::aes(fill = .data[[fill_col]]), colour = NA) +
    ggplot2::geom_sf(data = sig_sf, fill = NA, colour = "black", linewidth = 0.45) +
    ggplot2::geom_sf(data = afro_sf, fill = NA, colour = "black", linewidth = 0.12) +
    do.call(ggplot2::scale_fill_gradient2,
            c(list(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                   oob = scales::squish, na.value = "grey80", name = fill_name), scale_args)) +
    ggplot2::labs(title = title, subtitle = paste0(subtitle, RATIO_NOTE)) +
    ggplot2::theme_void(base_size = 11) + ggplot2::theme(legend.position = "right")
}
RATIO_SCALE <- list(limits = MAP_LIM, breaks = c(-2, -1, 0, 1, 2), labels = c("0.25x", "0.5x", "1x", "2x", "4x"))
PROB_SCALE  <- list(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                    labels = c("P=1.0\n(smaller)", "0.5", "0", "0.5", "P=1.0\n(larger)"))

window_df <- function(cutoff) if (is.null(cutoff)) adm1_month_both else dplyr::filter(adm1_month_both, month <= cutoff)

build_adm1_window <- function(cutoff, window_label, file_slug) {
  df <- window_df(cutoff)
  ratio_cases <- paired_adm1_stats(df, "cases", CI_LEVEL) %>%
    finish_adm1_ratio(LOW_SIGNAL_CASES, country_region, adm1_group_horizon, adm1_pop = adm1_pop, ratio_eps = RATIO_EPS)
  ratio_ob <- paired_adm1_stats(df, "outbreaks", CI_LEVEL) %>%
    finish_adm1_ratio(LOW_SIGNAL_OUTBREAKS, country_region, adm1_group_horizon, adm1_pop = NULL, ratio_eps = RATIO_EPS)
  readr::write_csv(ratio_cases, file.path(OUT_DIR, sprintf("adm1_cumulative_incidence_ratio_%s.csv", file_slug)))
  readr::write_csv(ratio_ob,    file.path(OUT_DIR, sprintf("adm1_cumulative_outbreak_ratio_%s.csv", file_slug)))
  message(sprintf("[07d] [%s] cases: %d units, %d low-signal, %d significant | outbreaks: %d units, %d low-signal, %d significant (at %.0f%% CI)",
                  window_label, nrow(ratio_cases), sum(ratio_cases$low_signal), sum(ratio_cases$significant, na.rm = TRUE),
                  nrow(ratio_ob), sum(ratio_ob$low_signal), sum(ratio_ob$significant, na.rm = TRUE), 100 * CI_LEVEL))
  start_txt <- format(RATIO_START, "%b %Y")
  save_fig(plot_adm1_diverging_map(
    ratio_cases, "log2_ratio", RATIO_SCALE, "Forecast /\ncounterfactual\ncum. cases",
    "Cumulative cholera cases: factual forecast vs ENSO counterfactual",
    sprintf("%s window, from %s (grey = <%d expected cases across both arms; black outline = significant paired diff. at %.0f%% CI)",
            window_label, start_txt, LOW_SIGNAL_CASES, 100 * CI_LEVEL)),
    sprintf("fig3_adm1_ratio_map_%s", file_slug), bg = "white")
  save_fig(plot_adm1_diverging_map(
    ratio_ob, "log2_ratio", RATIO_SCALE, "Forecast /\ncounterfactual\ncum. outbreaks",
    "Cumulative cholera outbreaks: factual forecast vs ENSO counterfactual",
    sprintf("%s window, from %s (grey = <%g expected outbreaks across both arms; black outline = significant paired diff. at %.0f%% CI)",
            window_label, start_txt, LOW_SIGNAL_OUTBREAKS, 100 * CI_LEVEL)),
    sprintf("fig4_adm1_outbreak_ratio_map_%s", file_slug), bg = "white")
  invisible(list(cases = ratio_cases, outbreaks = ratio_ob))
}

build_adm1_prob_window <- function(cutoff, window_label, file_slug) {
  prob_df <- adm1_prob_stats(window_df(cutoff), "cases", PROB_THRESHOLD) %>%
    finish_adm1_prob(LOW_SIGNAL_CASES, PROB_THRESHOLD, country_region, adm1_group_horizon)
  readr::write_csv(prob_df, file.path(OUT_DIR, sprintf("adm1_cumulative_incidence_prob_%s.csv", file_slug)))
  message(sprintf("[07d] [%s] prob-of-shift (+/-%.0f%%): %d units, %d low-signal, %d significant",
                  window_label, 100 * PROB_THRESHOLD, nrow(prob_df), sum(prob_df$low_signal, na.rm = TRUE),
                  sum(prob_df$significant, na.rm = TRUE)))
  save_fig(plot_adm1_diverging_map(
    prob_df, "signed_prob", PROB_SCALE, "P(shift)",
    sprintf("Probability cholera cases are >=%.0f%% larger/smaller: forecast vs ENSO counterfactual", 100 * PROB_THRESHOLD),
    sprintf("%s window, from %s (unpaired nrep x nrep replicate comparison; grey = <%d expected cases across both arms; black outline = P(>=%.0f%% shift) >= %.0f%%)",
            window_label, format(RATIO_START, "%b %Y"), LOW_SIGNAL_CASES, 100 * PROB_THRESHOLD, 100 * PROB_THRESHOLD)),
    sprintf("fig5_adm1_prob_shift_map_%s", file_slug), bg = "white")
  invisible(prob_df)
}

common_label <- sprintf("common, through %s (africa_min_horizon)", format(AFRICA_MIN_HORIZON, "%b %Y"))
full_label   <- sprintf("full per-group (%s)", FULL_WINDOW_LABEL)
message("[07d] ADM1 ratio maps (common + full windows)...")
build_adm1_window(AFRICA_MIN_HORIZON, common_label, "common")
build_adm1_window(NULL, full_label, "full")
message(sprintf("[07d] ADM1 probability-of-shift maps (threshold +/-%.0f%%)...", 100 * PROB_THRESHOLD))
build_adm1_prob_window(AFRICA_MIN_HORIZON, common_label, "common")
build_adm1_prob_window(NULL, full_label, "full")

message("[07d] done. Figures + CSV in: ", OUT_DIR)
