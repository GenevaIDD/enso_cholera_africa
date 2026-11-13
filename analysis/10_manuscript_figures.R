#!/usr/bin/env Rscript
# analysis/10_manuscript_figures.R
# =============================================================================
# Phase 10 — MANUSCRIPT (publication) figures. Plotting only.
#
# Figures
# -------
#   1  ENSO exposure-response overview (stacked)
#        a  observed-ONI strip over the data-coverage timeline (rows = full
#           country names, ordered by last case bubble; ENSO-phase background)
#        b  inset map: onsets per surveillance-covered month per gid
#        c  per ENSO group (columns), on a 2-year calendar axis: the
#           composite El Nino / La Nina ONI (09's oni_composite), the mean
#           monthly rainfall (rain_enso), and per track the ELPD lag weights
#           (deployed lag marked) over its seasonality relative to the annual
#           mean (09's seasonality)
#        d  ERF grid (rows = track, columns = group; 50% + 95% CrI;
#           direction arrows in the first column)
#   2  ENSO effects by ONI condition (stacked)
#        a  ADM1 IRR maps, rows = track, columns = ONI condition
#           (results.derivation.enso_conds), a direction key per row
#        b  under each ONI column, the share of the population in
#           ENSO-sensitive ADM1s (report_enso_exposure at
#           results.style.enso_headline), Africa + each ENSO group
#   3  Attribution of the past strong ENSO events (needs
#      results.event_manifest; the smoothing -> forecast Figure 3 is in
#      10_manuscript_figures_legacy.R)
#        a  event brackets over the observed-ONI strip over Africa's smoothing
#           (masked / unmasked) and observed series
#        b  per event (columns): ENSO-attributable fraction for Africa and each
#           ENSO group, posterior median + 50% / 95% (noise-adjusted)
#        c  ENSO-attributable fraction per ADM1 under each event's bars
#      Main: cases; supplementary: onsets.
#   4  The 2026 ENSO forecast against the neutral-ENSO counterfactual
#      (the paired-ratio attribution Figure 4 is in
#      10_manuscript_figures_legacy.R)
#        a  the ONI forecast strip over Africa + the ENSO regions from the
#           analogue window (forecast vs counterfactual), their cumulative
#           totals and the ratio of means 1 / (1 - AF) filled by the scope's
#           IPCC class (report_attrib_metrics ipcc_scope)
#        b  ADM1 ratio of means, IPCC hatching (the fig. AE map)
#        c  ENSO-attributable fraction, Africa + each ENSO region (af_posterior:
#           posterior median + 50% / 95%, noise-adjusted)
#      Main: cases; supplementary: onsets.
#   H  Hindcast forecast skill (results.scope.hindcast; builders in
#      utils/figures_hindcast.R, shared with notebooks/hindcast_skill.qmd)
#        a  outbreak-count CRPSS vs the climatology references, all leads and
#           the 0-7-month class (headline reference solid)
#        b  CRPSS against the headline reference by lead class
#        c  case totals by support: sim/obs ratio, CRPSS vs LOO climatology,
#           95% coverage (headline support solid)
#        d  reliability diagram
#
# Supplementary: Figure 1 ERF at country granularity; Figures 3 and 4 for onsets;
# 4-set ONI-envelope comparison (with Figure 4); hindcast
# calibration (PIT, rank histogram, coverage, ROC), per-country case skill and
# legacy-basis comparison.
#
# DATA CONTRACT — never derives, re-fits, or re-runs the bpfilter. Reads the
# phase-09 aggregates from generated_data/results/<tag>/ via results_require(),
# which refuses a missing, stale, or differently-decided artefact. Figure 1a
# (and the bubble check that places 1b) reads the survival / coverage / monthly
# data directly -- a narrow exception scoped to that panel; Figures 1b, 2 and 4
# read GADM / giscoR polygons as basemaps. The ONI strips and ENSO-phase
# shading of Figures 1 and 3 come from the phase-09 ONI artefact, Figure 1d's
# correlation tiles and rainy seasons from the phase-09 rain_enso artefact.
#
# Usage
# -----
#   Rscript analysis/10_manuscript_figures.R -C analysis/configs/results/manuscript.yml
#   Rscript analysis/10_manuscript_figures.R -C ... --figures 4,H --supplementary false
#
# --figures (comma-separated subset of 1,2,3,4,H, or the supplementary
# S5,S6,FC,SC,AE,EV,MB,MU,GM,KC,DC,SV,FA,AA,PC,SE,WR) and --supplementary
# (true|false) override results.scope.build_figures / .supplementary for this
# run only: they pick which figures are drawn, never which artefacts are read
# or how they were derived (RES and the analog sets stay as the config says).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  # The tidyverse sub-packages 10 and its sourced modules call bare, in
  # tidyverse's core-attach order (stringr, purrr, tibble and readr are not
  # called bare; readr is called as readr::).
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(lubridate)
  library(qs2)
  library(sf)
  library(scales)
  library(patchwork)
  library(ragg)
  library(taxdat)
  library(systemfonts)
})

# results_core.R: track/group registries, schema tokens, results_require(); it
# sources its own base/config/paths/provenance/enso closure. figures_cb.R:
# load_africa_sf(). figures_io.R: saveFigure(). deploy_core.R: so a non-empty
# results.deploy_config resolves (make_options_from_config() hard-stops if the
# pointer is set but resolve_results_inputs_from_deploy() is absent); inert when
# deploy_config is "". figures_manuscript_fig1 / _fig2 / _fig34.R: the Figure
# 1-4 builders (definitions only; they read this script's globals when called,
# so they are sourced here and nowhere else).
for (.m in c("results_core", "gid", "figures_cb", "figures_io", "deploy_core", "figures_hindcast",
             "figures_supp", "figures_manuscript_fig1", "figures_manuscript_fig2",
             "figures_manuscript_fig34"))
  source(file.path("analysis/utils", paste0(.m, ".R")))

# ---------------------------------------------------------------------------
# Section 0: CLI
# ---------------------------------------------------------------------------
# S5 / S6 / FC / SC / AE: supplementary (figures_supp.R); FC = forecast vs
# counterfactual by country, SC = historical reconstruction (smoothing) by
# country, AE = population in ADM1s with more ENSO-attributable cases, EV =
# attributable fraction of the past strong ENSO events (needs results.event_manifest).
# Keys are matched upper-cased.
# MB = membership uncertainty of the deployed classes (report_membership;
# needs results.kl_sweep_dir and results.kl_sweep_cell).
# Supplement maps: MU = modelled units by admin level
# (figS_modelled_units; 09 modelled_units + GADM), GM = the deployed partition
# (figS_group_map), KC = MAP class and probability per country for K = 2..5 at
# the pinned basis (figS_kl_classes_by_K; report_kl_sweep, needs
# results.kl_sweep_dir and results.kl_sweep_cell).
# Supplement figures (figures_supp.R): DC = outbreak data
# coverage by country and season year + by admin level (figS_data_coverage;
# report_dataset$outbreaks); SV = the historical reconstruction against the
# observations: PIT, seasonality, durations, active fraction, contagion
# conditionals (figS_smoothing_validation; smoothing v4 $validation); FA =
# monthly onsets and cases under the attribution factual and the counterfactual
# of every analogue set, Africa + each group (figS_forecast_analogues;
# attribution ribbons); AA = ADM1 maps of the ratio of means and of the
# probability of an increase, one row per analogue set
# (figS_attribution_analogues; report_attrib_exposure$adm1).
# SE = seasonality against the rainy season:
# rainfall climatology with the Fig. 1d rainy season, observed onset, cessation
# and case rates by calendar month with the fitted occurrence and monthly
# seasonality, and the fitted cessation seasonality where the duration fit has
# one (duration.seasonal_term), per ENSO group (figS_seasonality; 09 seasonality).
# WR = observed vs modelled ENSO ratio of the target year: the WHO
# ratio of 2026 to the ENSO-neutral years against the smoothing / ENSO-
# counterfactual ratio of means over their common months, per country of each
# derived group (figS_who_enso_ratio; 09 report_who_enso_ratio;
# style.who_enso_ratio picks the year definition and the analogue set).
# EG = ENSO teleconnection groups (GM, KC and MB panel a in one figure):
# map, membership probability per class and most probable class at K = 2..4
# (figS_enso_groups), and the class paths of MB panel b on their own
# (figS_membership_paths); needs results.kl_sweep_dir and results.kl_sweep_cell.
# One step per figure key (Section 10 defines the fig_step_*() functions, Section
# 11 runs the requested ones in this order); FIGURE_KEYS is the CLI contract.
FIGURE_STEPS <- c("1" = "fig_step_1", "2" = "fig_step_2", "3" = "fig_step_3", "4" = "fig_step_4",
                  H = "fig_step_H", S5 = "fig_step_S5", S6 = "fig_step_S6", FC = "fig_step_FC",
                  SC = "fig_step_SC", AE = "fig_step_AE", EV = "fig_step_EV", MB = "fig_step_MB",
                  MU = "fig_step_MU", GM = "fig_step_GM", KC = "fig_step_KC", DC = "fig_step_DC",
                  PC = "fig_step_PC", SV = "fig_step_SV", FA = "fig_step_FA", AA = "fig_step_AA",
                  SE = "fig_step_SE", WR = "fig_step_WR", EG = "fig_step_EG")
FIGURE_KEYS <- names(FIGURE_STEPS)
opt <- make_options_from_config(option_list = list(
  make_option("--figures", type = "character", default = NULL,
              help = paste0("Comma-separated figures to draw (subset of ",
                            paste(FIGURE_KEYS, collapse = ","),
                            ") [default: results.scope.build_figures]")),
  make_option("--supplementary", type = "character", default = NULL,
              help = "true|false: draw the supplementary variants [default: results.scope.supplementary]")
), model = NULL)

if (is.null(opt$results))
  stop("[10] Config must contain a `results:` block. Got: ",
       paste(names(opt), collapse = ", "))

CONFIG_PATH <- results_config_path()
RES   <- opt$results
SCOPE <- RES$scope
STYLE <- RES$style

# Group labels / order / colours of the results' fit registry
# (results.inputs.grouping), from its sidecar (utils/config.R::
# load_grouping_meta()); for enso_grouping these equal results_core.R's tables.
GROUPING_CSV <- RES$inputs$grouping %||% "analysis/groupings/enso_grouping.csv"
.grp_meta <- load_grouping_meta(sub("\\.csv$", "", basename(GROUPING_CSV)),
                                groupings_dir = dirname(GROUPING_CSV))
GROUP_LABELS    <- .grp_meta$label
GROUP_SHORT     <- .grp_meta$short
GROUP_ORDER     <- .grp_meta$order
GROUP_SHORT_FIG <- .grp_meta$short_fig
GROUP_COLOURS   <- groupColours(.grp_meta)

FIGURES  <- as.character(SCOPE$build_figures)
SUPP     <- isTRUE(SCOPE$supplementary)
# CLI overrides (selection only; SETS below still follows the config, so the
# artefacts read are the ones the config derived).
if (!is.null(opt$figures)) {
  FIGURES <- toupper(trimws(strsplit(opt$figures, ",", fixed = TRUE)[[1]]))
  FIGURES <- FIGURES[nzchar(FIGURES)]
  bad <- setdiff(FIGURES, FIGURE_KEYS)
  if (length(bad) || !length(FIGURES))
    stop("[10] --figures: unknown figure(s) '", paste(bad, collapse = ","),
         "'; choose from ", paste(FIGURE_KEYS, collapse = ","), call. = FALSE)
}
if (!is.null(opt$supplementary)) {
  sv <- tolower(trimws(opt$supplementary))
  if (!sv %in% c("true", "false", "yes", "no"))
    stop("[10] --supplementary must be true or false, got '", opt$supplementary, "'", call. = FALSE)
  if (sv %in% c("true", "yes") && !isTRUE(SCOPE$supplementary))
    stop("[10] --supplementary true needs results.scope.supplementary: true in the config ",
         "(the supplementary analog-set artefacts are only derived then)", call. = FALSE)
  SUPP <- sv %in% c("true", "yes")
}
MAIN_SET <- as.character(SCOPE$main_set)
SUPP_SETS <- if (SUPP) setdiff(as.character(SCOPE$supp_sets %||% character(0)), "") else character(0)
FIG1_TRACKS <- as.character(SCOPE$fig1_tracks)

# Tracks for the stacked Figure 1 / Figure 2 variants. With the size track
# retired downstream `results.scope.fig1_tracks` is exactly the stacked set
# (occurrence/duration/monthly), so the stack follows the config.
STACKED_TRACKS <- FIG1_TRACKS

# Artefact filenames are keyed on SETS: one artefact per set list.
SETS <- resolve_results_sets(RES)

FORMATS <- tolower(as.character(STYLE$formats))
DPI     <- as.integer(STYLE$dpi)
OUT_DIR <- STYLE$out_dir
ERF_XLIM <- as.numeric(unlist(STYLE$erf_xlim))

# Figure 2 presentation constants (not config surface).
#   FIG2_MAP_LO/HI      — IRR cap, symmetric on the log scale; also the ratio
#     colourbar of Figures 4b, AE and AA (ps_ratio_fill_scale()). Figure 1d's
#     ERF axis has its own range (FIG1_ERF_YLIM, up to 10).
#   FIG2_RULE_LWD / FIG2_RULE_TRIM — booktabs \cmidrule-style header rules
#     (ggplot linewidth, mm) under each ENSO phase and each ONI column, trimmed
#     at both ends (the phase-header fills are trimmed to the same width).
#   Figure height: fitted to the layout (fig2_fit_height()), so the maps fill
#     the WIDTH_MM width with no space left around them.
#   FIG2_PHASE_FILL     — La Nina / El Nino header fills: Figure 1's ENSO-phase
#     colours at PS_ENSO_PHASE_ALPHA.
#   FIG2_ROW_KEY_*      — the direction key right of each track row (labels and
#     column width).
#   Map line widths: white ADM1 borders FIG2_ADM1_LWD, significant-ADM1
#     outlines FIG2_SIG_LWD, black country borders FIG2_COUNTRY_BORDER_LWD.
#   Map layout: track titles FIG2_STRIP_GAP_MM from the
#     maps, maps FIG2_SPACING_MM apart (x, y) with no coord_sf() expansion,
#     outer margin FIG2_MARGIN_MM, the legend row FIG2_LEGEND_GAP_MM under the
#     maps (no legend margin).
#   Panel b (population share in ENSO-sensitive areas, under each ONI column):
#     a FIG2_BAR_GAP_MM spacer row (tag "b" + FIG2_BAR_TITLE), the bars
#     (FIG2_BAR_ROW_MM; bar thickness FIG2_BAR_WIDTH of a row) and their axis
#     (FIG2_BAR_AXIS_MM); shares above FIG2_LABEL_INSIDE are printed inside
#     the bar.
FIG2_MAP_LO <- 1 / 3
FIG2_MAP_HI <- 3
FIG2_COUNTRY_BORDER_LWD <- 0.1
FIG2_ADM1_LWD  <- 0.002
FIG2_SIG_LWD   <- 0.05
FIG2_RULE_LWD  <- 0.15
FIG2_RULE_TRIM <- grid::unit(3, "mm")
FIG2_STRIP_GAP_MM  <- 0.5
FIG2_SPACING_MM    <- c(3, 2)
FIG2_MARGIN_MM     <- 1
FIG2_LEGEND_GAP_MM <- 1
FIG2_BAR_GAP_MM    <- 4.5
FIG2_BAR_ROW_MM    <- 20
FIG2_BAR_AXIS_MM   <- 3.5
FIG2_BAR_WIDTH     <- 0.72
FIG2_LABEL_INSIDE  <- 0.85
FIG2_BAR_TITLE     <- "Population in ENSO-sensitive areas (%)"
# FIG2_PHASE_FILL is set in Section 2b, from PS_ENSO_PHASE_COL / _ALPHA.
# What red / blue means on each track row (Figure 1d's arrow labels; the
# cessation IRR is inverted in the maps, so red = lower cessation hazard =
# longer outbreaks).
FIG2_ROW_KEY_LABELS <- list(
  "Outbreak occurrence" = c(up = "more\nonsets",      down = "fewer\nonsets"),
  "Outbreak cessation"  = c(up = "longer\noutbreaks",  down = "shorter\noutbreaks"),
  "Monthly incidence"   = c(up = "more\ncases",       down = "fewer\ncases"))
FIG2_ROW_KEY_MM <- 15
# The group-map inset: fig. GM's partition
# map, unlabelled, in the row-key column right of panel b's last bar panel,
# its bottom FIG2_INSET_RULE_GAP_MM above the Africa rule, growing upward into
# the legend row; FIG2_INSET_GAP_MM off the bars at the column's left,
# FIG2_INSET_SHIFT_MM left of the column's right edge (it keeps its width, so
# it overhangs into the last bar column's empty top-right corner); white
# country borders FIG2_INSET_BORDER_LWD. NULL FIG2_INSET_SHIFT_MM = no inset.
FIG2_INSET_SHIFT_MM    <- 6
FIG2_INSET_GAP_MM      <- 1.5
FIG2_INSET_RULE_GAP_MM <- 0.5
FIG2_INSET_BORDER_LWD  <- 0.08
# The La Niña / El Niño header row: TRUE fills it with FIG2_PHASE_FILL; FALSE
# leaves it white, the rules kept.
FIG2_PHASE_HEADER_FILL <- FALSE


# Figure 1 presentation constants, not config surface, like the FIG2_* ones.
# Figure 1b inset map (panel a's top-right corner):
#   FIG1_INSET_MAX_ROWS — the inset reaches down at most this many panel-a rows,
#     fewer if a deeper box (wider, by the map's fixed aspect) would cover a
#     case bubble.
#   FIG1_INSET_MARGIN   — gap to panel a's top/right border (panel-width fraction).
#   FIG1_INSET_ASPECT   — map height / width of the inset's lon/lat box.
#   FIG1_PANEL_A_ASPECT — panel a's width / height in the 183 x 140 mm layout,
#     measured from the render (1307 x 2632 px at 600 dpi);
#     re-measure if the figure size, panel a's share of the width, the ONI strip
#     or panel a's tick / legend / country-label text changes.
# Panel a's rows: one per country, by date of its last case bubble, latest at
# the bottom (no group blocks; the group is the bubble colour), labelled with
# FIG1_COUNTRY_NAMES (every country of the ENSO groupings; "\n" splits the
# longest name so it does not set panel a's left edge).
# Layout: FIG1_HEIGHT_MM tall at WIDTH_MM; left (a) : right (c / d / e) columns
# FIG1_WIDTHS; FIG1_ONI_STRIP_MM the ONI strip above a; FIG1_A_RIGHT_MM the
# margin right of a (room for its last year label; no other margin between the
# columns); FIG1_A_LEFT_MM / FIG1_A_TOP_MM the left column's outer left / top
# margins (the top one puts the "a" tag, drawn inside the ONI strip, level with
# "c"); FIG1_A_XLIM panel a's window (every coverage period and case bubble,
# with a little room). The right column's own constants (heights, e's bottom
# margin aligning e's panel bottom with a's, ...) are in Section 9a.
FIG1_GADM_PATH      <- "data/geodata/gadm_valid.gpkg"
FIG1_INSET_BBOX     <- c(xmin = -19, xmax = 52.5, ymin = -36.5, ymax = 37.5)
FIG1_INSET_MAX_ROWS <- 16L
FIG1_INSET_MARGIN   <- 0.006
FIG1_INSET_ASPECT   <- unname((FIG1_INSET_BBOX["ymax"] - FIG1_INSET_BBOX["ymin"]) /
                              (FIG1_INSET_BBOX["xmax"] - FIG1_INSET_BBOX["xmin"]))
FIG1_PANEL_A_ASPECT <- 1307 / 2632
FIG1_HEIGHT_MM      <- 140
FIG1_WIDTHS         <- c(0.42, 0.58)
FIG1_ONI_STRIP_MM   <- 11
FIG1_A_RIGHT_MM     <- 2.5
FIG1_A_LEFT_MM      <- 0.5
FIG1_A_TOP_MM       <- 2.96
FIG1_A_XLIM         <- as.Date(c("2009-10-01", "2025-01-31"))
FIG1_COUNTRY_NAMES  <- c(
  AGO = "Angola", BDI = "Burundi", BEN = "Benin", BFA = "Burkina Faso", BWA = "Botswana",
  CAF = "Central African\nRepublic", CIV = "Côte d'Ivoire", CMR = "Cameroon",
  COD = "DR Congo", COG = "Congo", COM = "Comoros", DJI = "Djibouti", ETH = "Ethiopia",
  GAB = "Gabon", GHA = "Ghana", GIN = "Guinea", GNB = "Guinea-Bissau",
  GNQ = "Equatorial Guinea", KEN = "Kenya", LBR = "Liberia", LSO = "Lesotho",
  MDG = "Madagascar", MLI = "Mali", MOZ = "Mozambique", MRT = "Mauritania", MWI = "Malawi",
  NAM = "Namibia", NER = "Niger", NGA = "Nigeria", RWA = "Rwanda", SDN = "Sudan",
  SEN = "Senegal", SLE = "Sierra Leone", SOM = "Somalia", SSD = "South Sudan",
  STP = "São Tomé and Príncipe", SWZ = "Eswatini", TCD = "Chad", TGO = "Togo",
  TZA = "Tanzania", UGA = "Uganda", ZAF = "South Africa", ZMB = "Zambia", ZWE = "Zimbabwe")
# + the synthetic world's fake countries (gid.R::syntheticCountries(); no real row carries their XA? codes)
FIG1_COUNTRY_NAMES  <- c(FIG1_COUNTRY_NAMES, local({
  syn <- syntheticCountries(); stats::setNames(as.character(syn$name), syn$iso3) }))

for (d in c(OUT_DIR, file.path(OUT_DIR, "supplementary")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[10] figures=%s supplementary=%s sets=%s out=%s",
                paste(FIGURES, collapse = ","), SUPP,
                paste(SETS, collapse = ","), OUT_DIR))

# ---------------------------------------------------------------------------
# Section 1: publication style
# ---------------------------------------------------------------------------
WIDTH_MM  <- 183      # double-column
BASE_SIZE <- 7        # pt

#' First available family from <preferred> -> Arial -> Helvetica -> Nimbus Sans
#' -> DejaVu Sans -> "sans". Any substitution is reported.
resolve_font <- function(preferred) {
  chain <- unique(c(preferred, "Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"))
  chain <- chain[!is.na(chain) & nzchar(chain)]
  avail <- tryCatch(unique(systemfonts::system_fonts()$family),
                    error = function(e) character(0))
  hit <- chain[chain %in% avail]
  if (!length(hit)) {
    message("[10] font: none of ", paste(chain, collapse = " / "),
            " installed; falling back to the device default 'sans'")
    return("sans")
  }
  if (!identical(hit[1], chain[1]))
    message("[10] font: '", chain[1], "' not installed; using '", hit[1], "'")
  hit[1]
}
FONT <- resolve_font(STYLE$font)

theme_ms <- function(base_size = BASE_SIZE) {
  ggplot2::theme_bw(base_size = base_size, base_family = FONT) +
    ggplot2::theme(
      panel.background  = element_blank(),
      panel.border      = element_rect(colour = "grey25", fill = NA, linewidth = 0.3),
      panel.grid.major  = element_line(colour = "grey93", linewidth = 0.2),
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, 
                                       margin = margin(b = 1.2, unit = "mm")),
      axis.ticks        = element_line(colour = "grey25", linewidth = 0.25),
      axis.ticks.length = unit(0.8, "mm"),
      axis.text         = element_text(colour = "grey20", size = base_size - 0.5),
      axis.title        = element_text(size = base_size),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(4.5, "mm"),
      legend.text       = element_text(size = base_size - 0.5),
      legend.title      = element_text(size = base_size - 0.5),
      plot.title        = element_text(size = base_size, hjust = 0,
                                       margin = margin(b = 1, unit = "mm")),
      plot.margin       = margin(1.5, 1.5, 1.5, 1.5, unit = "mm"),
      plot.tag          = element_text(size = base_size + 2, face = "bold",
                                       family = FONT, hjust = 0, vjust = 1)
    )
}

TAG_THEME <- ggplot2::theme(
  plot.tag = element_text(size = BASE_SIZE + 2, face = "bold", family = FONT,
                          hjust = 0, vjust = 1))

#' Write one figure as vector PDF (cairo) and/or high-DPI PNG (ragg); sizes in
#' millimetres. One saveFigure() call per format, since the two formats need
#' different devices.
save_figure <- function(plot, stem, height_mm, width_mm = WIDTH_MM, subdir = NULL) {
  dir <- if (is.null(subdir)) OUT_DIR else file.path(OUT_DIR, subdir)
  if ("pdf" %in% FORMATS)
    saveFigure(plot, stem, dir, width = width_mm, height = height_mm,
               units = "mm", formats = "pdf", device = grDevices::cairo_pdf,
               bg = "white", limitsize = FALSE)
  if ("png" %in% FORMATS)
    saveFigure(plot, stem, dir, width = width_mm, height = height_mm,
               units = "mm", dpi = DPI, formats = "png",
               device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

# ---------------------------------------------------------------------------
# Section 2: palettes and labels
# ---------------------------------------------------------------------------
# Forecast colour of the supplementary ONI-envelope comparison (black, as the
# forecast ONI of Figure 4a).
FC_COLOR <- "black"

# The phase sets in the ENSO-phase colours (Section 2b's PS_ENSO_PHASE_COL,
# solid; restated here, Section 2b comes later; neutral a mid grey, as a line
# needs one); every analogue in crimson, apart from the sand of El Nino.
SET_COLORS <- c(all = "#b2182b", nina = "#2f3e4c", neutral = "#878787",
                elnino = "#c49a5a")
SET_LABELS <- c(all = "All analogues", nina = "La Ni\u00f1a analogues",
                neutral = "Neutral analogues", elnino = "El Ni\u00f1o analogues")

# TRACK_* / GROUP_* / COMPOSITE_CONDS / METRIC_* live in utils/results_core.R;
# ERF_XLIM comes from results.style in Section 0.

set_lab <- function(s) { l <- unname(SET_LABELS[as.character(s)]); ifelse(is.na(l), as.character(s), l) }

# ---------------------------------------------------------------------------
# Section 2b: cross-figure conventions (Figures 1-4)
# ---------------------------------------------------------------------------
#   width        — one double-column width for every figure (WIDTH_MM, the
#                  save_figure() default).
#   typeface     — FONT everywhere, maps included (taxdat::map_theme() leaves the
#                  family to the device default): ps_map_font().
#   ENSO phase   — sand (El Nino) / no fill / slate (La Nina) at alpha 0.14,
#                  labels El Nino / Neutral / La Nina (PS_ENSO_*: none of the
#                  group colours (Dark2) nor the maps' red / blue; Neutral
#                  unshaded, white, so La Nina stands apart). Shading bands at PS_ENSO_PHASE_ALPHA,
#                  solid where a line or text needs the colour (Figure 3's
#                  event brackets, the analogue-set lines).
#   legends      — title BASE_SIZE - 1, labels BASE_SIZE - 1.5 (ps_legend_theme()).
#   counterfactual wording — "neutral-ENSO counterfactual" (ps_cf_label()).
#   uncertainty  — median + 50% + 95% bands (phase-09 q025 / q975 columns).
#   map greys    — two states with a legend key: "Not modelled" (outside the
#                  fitted ENSO groups, or no outbreak in the survival panel) and
#                  "Modelled, no ADM1 estimate" (outbreak countries / ADM1s
#                  without an estimate, incl. low-signal ADM1s in Figure 4b).

PS_ENSO_PHASE_COL   <- c(elnino = "#c49a5a", neutral = "white", nina = "#2f3e4c")
PS_ENSO_PHASE_LAB   <- c(elnino = "El Niño", neutral = "Neutral", nina = "La Niña")
PS_ENSO_PHASE_ALPHA <- 0.14
# Figure 2's La Nina / El Nino header fills: the phase colours at the band alpha.
FIG2_PHASE_FILL <- stats::setNames(scales::alpha(unname(PS_ENSO_PHASE_COL[c("nina", "elnino", "neutral")]),
                                                 PS_ENSO_PHASE_ALPHA), c("nina", "elnino", "neutral"))

PS_GREY <- c(not_modelled = "grey85", no_estimate = "grey93")
PS_GREY_LAB <- c(not_modelled = "Not modelled", no_estimate = "Modelled, no ADM1 estimate")
# Same labels wrapped, for narrow legends (Figure 1b inset).
PS_GREY_LAB_WRAP <- c(not_modelled = "Not modelled", no_estimate = "Modelled,\nno ADM1 estimate")

# K for thousands everywhere.
# ONI strip cell height (mm), both Figure 3 strips.
PS_ONI_STRIP_MM <- 14
PS_LABEL_K     <- scales::label_number(scale_cut = scales::cut_short_scale())
# One ONI y range for every ONI strip (Figures 1a, 3a, 4a): 2015/16 peak 2.6,
# CFSv2 high 2.56, headroom on top.
ONI_YLIM <- c(-2.7, 3.2)

#' Legend text sizes shared by every figure.
ps_legend_theme <- function()
  ggplot2::theme(legend.title = ggplot2::element_text(size = BASE_SIZE - 1),
                 legend.text  = ggplot2::element_text(size = BASE_SIZE - 1.5))

#' FONT for maps built on taxdat::map_theme() / theme_void().
ps_map_font <- function() ggplot2::theme(text = ggplot2::element_text(family = FONT))

#' Counterfactual label: "neutral-ENSO counterfactual" for the main set,
#' "<La Nina / El Nino / all>-analogue counterfactual" otherwise.
ps_cf_label <- function(aset, capital = FALSE) {
  lab <- if (identical(aset, "neutral")) "neutral-ENSO counterfactual"
         else sprintf("%s-analogue counterfactual", sub(" analogues$", "", set_lab(aset)))
  if (capital) paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab))) else lab
}

#' Country outlines with a `grey_status` column (not_modelled / no_estimate).
#' `modelled_countries`: ISO3 codes in the fitted ENSO groups with >= 1 onset
#' in the survival panel (phase-09 onset_rate artefact).
ps_country_status <- function(afro_sf, modelled_countries) {
  afro_sf$grey_status <- factor(ifelse(afro_sf$ISO3_CODE %in% modelled_countries,
                                       "no_estimate", "not_modelled"),
                                levels = names(PS_GREY))
  afro_sf
}

#' Base country layer coloured by grey_status, with its own legend, followed by
#' ggnewscale::new_scale_fill() so the map's estimate fill scale stays separate.
#' `legend_order` places the key among the map's legends.
#' `show_legend = FALSE`: same greys, no key.
ps_country_grey_layer <- function(country_sf, legend_order = 2, labels = PS_GREY_LAB, show_legend = TRUE) {
  list(
    ggplot2::geom_sf(data = country_sf, ggplot2::aes(fill = grey_status), colour = NA),
    ggplot2::scale_fill_manual(values = PS_GREY, labels = labels, name = NULL, drop = FALSE,
                               guide = if (!show_legend) "none" else ggplot2::guide_legend(order = legend_order,
                                                             override.aes = list(colour = "grey60",
                                                                                 linewidth = 0.1))),
    ggnewscale::new_scale_fill())
}

#' ENSO-phase background rectangles (one layer per phase, constant fill, so it
#' never touches a panel's own fill scale and survives faceting). `phases`:
#' TL / TR / phase columns; NULL-safe.
ps_phase_layer <- function(phases) {
  if (is.null(phases) || !nrow(phases)) return(NULL)
  lapply(split(phases, as.character(phases$phase), drop = TRUE), function(ph)
    ggplot2::geom_rect(data = ph, ggplot2::aes(xmin = TL, xmax = TR, ymin = -Inf, ymax = Inf),
                       fill = PS_ENSO_PHASE_COL[[as.character(ph$phase[1])]],
                       alpha = PS_ENSO_PHASE_ALPHA, inherit.aes = FALSE))
}

# --- Shared cross-figure helpers, Figures 1-4 ------------------------------
# One set used by every figure that shows these elements:
#   PS_AFRICA_COL       — Africa's colour wherever Africa is a series or a row
#                         (F3_AFRICA_COL is an alias).
#   PS_AFRICA_RULE / PS_AFRICA_GAP_ROWS / ps_africa_rows() — the Africa row of a
#                         row-per-scope panel sits PS_AFRICA_GAP_ROWS rows above
#                         the ENSO groups, with a thin rule centred in the gap.
#   PS_ERRORBAR / ps_interval_layers() — one interval format: the 95% interval
#                         as a thin capped error bar, the 50% interval as a
#                         thicker capped error bar on top, the point estimate as
#                         a point; optional arrows where the 95% interval runs
#                         past the axis limits. Linewidths in ggplot units, caps
#                         in category units (one row / column = 1).
#   PS_ENSO_STATE_TITLE / PS_ENSO_KEY_OUTLINE / PS_KEY_SPACING_MM /
#   ps_enso_phase_guide() — the ENSO-state legend: swatches at the bands' own
#                         alpha (PS_ENSO_PHASE_ALPHA) with a thin outline (a true
#                         0.2 mm; one ggplot linewidth unit = .pt * 25.4 / 96 mm),
#                         1 mm between a label and the next key.
#   PS_TAG_GP           — panel-tag text as grid gpar (the values of TAG_THEME),
#                         for tags drawn as grobs.
PS_AFRICA_COL      <- "#7b1e3c"
PS_AFRICA_RULE     <- list(colour = "grey55", linewidth = 0.3)
PS_AFRICA_GAP_ROWS <- 0.5
PS_ERRORBAR <- list(lwd95 = 0.25, cap95 = 0.45, lwd50 = 0.6, cap50 = 0.22, colour = "grey15",
                    point_size = 1.1, point_shape = 21, point_fill = "white", point_stroke = 0.4,
                    arrow_mm = 0.9)
PS_ENSO_STATE_TITLE <- "ENSO state"
PS_ENSO_KEY_OUTLINE <- list(colour = "grey60", linewidth = 0.2 / (ggplot2::.pt * 25.4 / 96))
PS_KEY_SPACING_MM   <- 1
PS_TAG_GP <- grid::gpar(fontsize = BASE_SIZE + 2, fontface = "bold", fontfamily = FONT)

#' Rows of a row-per-scope panel with Africa set apart.
#' @param scope factor whose levels are the rows TOP to bottom (Africa first)
#' @param africa the Africa level
#' @return list(y = numeric row of each element of `scope`, breaks, labels (for
#'   scale_y_continuous()), rule_y (the rule's position), rule = the rule layer)
ps_africa_rows <- function(scope, africa = "Africa", gap = PS_AFRICA_GAP_ROWS) {
  lv <- levels(scope); n <- length(lv)
  pos <- stats::setNames(rev(seq_len(n)), lv)          # top row = n
  if (africa %in% lv) pos[[africa]] <- pos[[africa]] + gap
  rule_y <- if (africa %in% lv) n - 0.5 + gap / 2 else NA_real_
  list(y = unname(pos[as.character(scope)]), breaks = unname(pos), labels = lv, rule_y = rule_y,
       rule = if (!is.na(rule_y))
         ggplot2::geom_hline(yintercept = rule_y, colour = PS_AFRICA_RULE$colour,
                             linewidth = PS_AFRICA_RULE$linewidth))
}

#' The 95% / 50% interval layers and the point estimate, from column names.
#' @param est,lo,q25,q75,hi columns of the point estimate and the 2.5 / 25 / 75
#'   / 97.5% bounds (q25 = NULL: 95% interval only)
#' @param orientation "y" = values on x, categories on y (horizontal intervals);
#'   "x" = values on y
#' @param pos the category-axis column (for the clip arrows)
#' @param clip value-axis limits; a 95% interval past one gets an arrow at it
#' @param colour the bars' colour; NULL = the plot's colour mapping (e.g. one
#'   colour per series; Figure 4a's cumulative column)
#' @param point FALSE: the bars only (the caller keeps its own point layer)
#' @param data NULL = the plot's data; else a layer data (a data frame or a
#'   function of the plot's data), e.g. one subset per line type
#' @param linetype the bars' line type ("solid"; "dashed" marks an interval
#'   that is not a posterior interval, ps_af_posterior_cols())
#' @return a list of layers (inherits the plot's category mapping)
ps_interval_layers <- function(est, lo, hi, q25 = NULL, q75 = NULL, orientation = c("y", "x"),
                               pos = NULL, clip = NULL, spec = PS_ERRORBAR, colour = spec$colour,
                               point = TRUE, data = NULL, linetype = "solid") {
  orientation <- match.arg(orientation)
  horiz <- orientation == "y"
  bar <- function(a, b, width, lwd) {
    args <- list(if (horiz) ggplot2::aes(xmin = .data[[a]], xmax = .data[[b]])
                 else ggplot2::aes(ymin = .data[[a]], ymax = .data[[b]]),
                 width = width, linewidth = lwd, orientation = orientation)
    if (!is.null(colour)) args$colour <- colour
    if (!is.null(data)) args$data <- data
    if (!identical(linetype, "solid")) args$linetype <- linetype
    do.call(ggplot2::geom_errorbar, args)
  }
  pt_args <- list(if (horiz) ggplot2::aes(x = .data[[est]]) else ggplot2::aes(y = .data[[est]]),
                  shape = spec$point_shape, size = spec$point_size, fill = spec$point_fill,
                  colour = spec$colour, stroke = spec$point_stroke)
  if (!is.null(data)) pt_args$data <- data
  out <- list(bar(lo, hi, spec$cap95, spec$lwd95),
              if (!is.null(q25)) bar(q25, q75, spec$cap50, spec$lwd50),
              if (point) do.call(ggplot2::geom_point, pt_args))
  if (!is.null(clip)) {
    if (is.null(pos)) stop("ps_interval_layers(): `clip` needs `pos`, the category-axis column.", call. = FALSE)
    w <- diff(clip); arr <- grid::arrow(length = grid::unit(spec$arrow_mm, "mm"), type = "closed")
    arrow_at <- function(keep, from, to) {
      a <- if (horiz)
        ggplot2::aes(x = from, xend = to, y = .data[[pos]], yend = .data[[pos]])
      else ggplot2::aes(y = from, yend = to, x = .data[[pos]], xend = .data[[pos]])
      ggplot2::geom_segment(data = keep, mapping = a, inherit.aes = FALSE, linewidth = spec$lwd95,
                            colour = spec$colour, arrow = arr)
    }
    out <- c(out, list(
      arrow_at(function(d) d[!is.na(d[[lo]]) & d[[lo]] < clip[1], , drop = FALSE], clip[1] + 0.02 * w, clip[1] - 0.04 * w),
      arrow_at(function(d) d[!is.na(d[[hi]]) & d[[hi]] > clip[2], , drop = FALSE], clip[2] - 0.02 * w, clip[2] + 0.04 * w)))
  }
  Filter(Negate(is.null), out)
}

#' The bar and interval columns of a posterior AF table (09's af_posterior /
#' event_af_posterior `scope`), in %: an identified scope takes the
#' noise-adjusted median and 50% / 95% bounds; a scope whose posterior spread
#' is unidentified (af_identified FALSE: noise share >=
#' results.thresholds.af_noise_share_max) takes its pooled value and the RAW
#' quantiles, flagged `unidentified` so the panel draws that interval dashed
#' (ps_af_interval_layers()): it carries the replicates' epidemic noise and is
#' not a posterior interval. A table without the flag (older artefacts) is
#' read as identified.
ps_af_posterior_cols <- function(d, scale = 100) {
  ni <- if ("af_identified" %in% names(d)) d$af_identified %in% FALSE else rep(FALSE, nrow(d))
  # Computed outside the data mask: inside mutate() a helper named like a
  # dplyr verb (pick) would be masked by it.
  col_of <- function(adj, raw) scale * ifelse(ni, d[[raw]], d[[adj]])
  d$unidentified <- ni
  d$af_bar <- col_of("af_adj_med", "af_pool")
  d$af_lo  <- col_of("af_adj_lo", "af_raw_lo");  d$af_q25 <- col_of("af_adj_q25", "af_raw_q25")
  d$af_q75 <- col_of("af_adj_q75", "af_raw_q75"); d$af_hi  <- col_of("af_adj_hi", "af_raw_hi")
  d
}

#' ps_interval_layers() of an AF panel split by ps_af_posterior_cols()'s flag:
#' solid posterior intervals, dashed raw intervals of the unidentified scopes.
ps_af_interval_layers <- function(..., unidentified = "unidentified") {
  keep <- function(flag) function(d) d[(d[[unidentified]] %in% TRUE) == flag, , drop = FALSE]
  c(ps_interval_layers(..., data = keep(FALSE)),
    ps_interval_layers(..., data = keep(TRUE), linetype = "22"))
}

#' The note of a panel with unidentified scopes, "" when there is none; `sep`
#' between its two halves ("\n" for a narrow legend cell).
ps_af_unidentified_note <- function(d, unidentified = "unidentified", sep = " ")
  if (any(d[[unidentified]] %in% TRUE))
    paste0("Dashed: posterior spread not identified", sep, "(pooled value, replicate range)") else ""

#' The ENSO-state legend guide (for a fill scale on PS_ENSO_PHASE_COL): title
#' PS_ENSO_STATE_TITLE, swatches at the bands' alpha with PS_ENSO_KEY_OUTLINE,
#' PS_KEY_SPACING_MM between a label and the next key (applied here, per
#' legend; ps_legend_theme() is unchanged). `order`: the legend's place among a
#' plot's legends (0 = ggplot's default). `...` goes to guide_legend().
ps_enso_phase_guide <- function(order = 0L, ...) {
  args <- utils::modifyList(list(title = PS_ENSO_STATE_TITLE, nrow = 1L, order = order,
                                 override.aes = list(alpha = PS_ENSO_PHASE_ALPHA,
                                                     colour = PS_ENSO_KEY_OUTLINE$colour,
                                                     linewidth = PS_ENSO_KEY_OUTLINE$linewidth),
                                 theme = ggplot2::theme(legend.key.spacing.x = grid::unit(PS_KEY_SPACING_MM, "mm"))),
                            list(...))
  do.call(ggplot2::guide_legend, args)
}

# ---------------------------------------------------------------------------
# Section 2c: shared data access and map helpers
# ---------------------------------------------------------------------------
#' One phase-09 artefact of `kind` (a RESULTS_SCHEMA key) via results_require():
#' refuses a missing, stale or differently-decided artefact. `...` are the
#' artefact's keys (track = , sets = ). Read once per run: the artefact is
#' memoised on its resolved path (.results_cache), so the builders that share
#' one (oni, attribution, smoothing, onset_rate, ...) pay one qs read and one
#' provenance check. Callers edit copies (R's copy-on-modify), never the cache.
.results_cache <- new.env(parent = emptyenv())
read_results <- function(kind, ...) {
  path <- makeResultsRunFile(opt, paste0("results_", kind), ...,
                             schema = RESULTS_SCHEMA[[kind]], config_path = CONFIG_PATH)
  if (is.null(.results_cache[[path]]))
    .results_cache[[path]] <- results_require(path, CONFIG_PATH, RES)
  .results_cache[[path]]
}

# Group display labels (results.style.group_labels):
# "registry" keeps the sidecar's names (every table above unchanged); "lag"
# relabels the groups as classes lettered by the deployed occurrence lag,
# longest lead first ("Class A (ENSO lead 6)", short "A -6"), read from the
# phase-09 group_labels artefact (the lags come from the fits); "centroid"
# prints the sidecar's centroid names ("West [GW]", short "GW", display order
# GW, GG, GH, GE, GS). GROUP_SHORT stays the stored data key; GROUP_COLOURS
# stay keyed by id.
GROUP_LABEL_MODE <- as.character(STYLE[["group_labels"]] %||% "registry")
GROUP_DISP <- group_display_meta(
  .grp_meta, if (identical(GROUP_LABEL_MODE, "lag")) read_results("group_labels")$groups else NULL,
  GROUP_LABEL_MODE)
if (!identical(GROUP_LABEL_MODE, "registry")) {
  GROUP_LABELS    <- GROUP_DISP$label
  GROUP_SHORT_FIG <- GROUP_DISP$short_fig
  GROUP_ORDER     <- GROUP_DISP$order
}
# Legend / header title of the group key.
GROUP_KEY_TITLE <- if (identical(GROUP_LABEL_MODE, "lag")) "ENSO lag class" else "ENSO group"
#' Stored group strings (ids; the stored short names in region / scope
#' columns) -> display labels; identity under "registry".
relab_groups <- function(df, cols = NULL) relabel_groups(df, GROUP_DISP, cols)
#' Africa then the groups in GROUP_ORDER, as display labels: the scope order of
#' the supplementary figures that read stored short names (AE, EV, SV).
display_scope_levels <- function() relab_groups(c("Africa", unname(GROUP_SHORT[GROUP_ORDER])))

# Membership certainty of the deployed classes (report_membership): each country's
# p_max (largest posterior class probability in the
# deployed sweep cell); countries below style.membership_low are marked in
# Figures 1-2 and the supplementary country ERFs. Read once, on first use
# (read_results() memoises); NULL (no markers) when results.kl_sweep_dir /
# kl_sweep_cell are unset.
MEMBERSHIP_LOW <- as.numeric(STYLE[["membership_low"]] %||% 0.6)
membership <- function() {
  if (!nzchar(RES$kl_sweep_dir %||% "") || !nzchar(RES[["kl_sweep_cell"]] %||% "")) return(NULL)
  read_results("report_membership")
}
#' ISO3 -> p_max (fitted countries), or NULL.
membership_pmax <- function() {
  m <- membership(); if (is.null(m)) return(NULL)
  x <- m$country[!m$country$carried_over, , drop = FALSE]
  stats::setNames(x$p_max, x$iso3)
}
#' Countries whose class is uncertain (p_max < style.membership_low).
membership_low_isos <- function() {
  p <- membership_pmax(); if (is.null(p)) character(0) else names(p)[p < MEMBERSHIP_LOW]
}
MEMBERSHIP_MARK <- "\u2020"

#' The ENSO exposure-response artefact of every track in `tracks`, named by track.
read_enso_relation <- function(tracks)
  stats::setNames(lapply(tracks, function(trk) read_results("enso_relation", track = trk)), tracks)

#' Modelled countries for the map greys: >= 1 onset in the survival panel
#' (phase-09 onset_rate; Figure 1b's definition). `rate`: that artefact, if
#' already read.
modelled_countries <- function(rate = read_results("onset_rate"))
  unique(rate$country[rate$onsets > 0])

#' The Africa basemap polygons Figures 2 and 4 draw on: ADM1 units, and the
#' country outlines without the small islands (drop_small_islands()). Loaded
#' once per run (.map_cache; seven figures draw on them); a caller that
#' simplifies or joins its copy leaves the cache untouched.
.map_cache <- new.env(parent = emptyenv())
load_africa_maps <- function(what) {
  if (is.null(.map_cache$maps)) {
    geo <- load_africa_sf()
    if (is.null(geo$adm1_sf)) stop("[10] ", what, ": ADM1 shapefile failed to load")
    .map_cache$maps <- list(adm1_sf = geo$adm1_sf, afro_sf = drop_small_islands(geo$afro_sf))
  }
  .map_cache$maps
}

#' utils/figures_cb.R::load_africa_basemap() (lakes + rivers, st_intersects
#' clip to `afro_sf`), read once per run: every caller clips to
#' load_africa_maps()$afro_sf, so one basemap serves them all.
africa_basemap <- function(afro_sf) {
  if (is.null(.map_cache$basemap)) .map_cache$basemap <- load_africa_basemap(afro_sf)
  .map_cache$basemap
}

#' Full-detail GADM polygons simplified for a supplement map, whose PDF size
#' they otherwise set. Planar GEOS (rmapshaper's V8 crashed on the whole ADM1
#' layer, s2 rejects some GADM rings); `tol` 0.02 degrees (~2 km) is far below a
#' map panel's ~100 km per mm, so there is no visible change.
simplify_map_sf <- function(x, tol = 0.02) {
  s2_was <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(s2_was)))
  suppressWarnings(sf::st_simplify(x, preserveTopology = TRUE, dTolerance = tol)) %>% sf::st_make_valid()
}

#' The GADM ADM_1 + ADM_2 polygons of the countries `isos` (FIG1_GADM_PATH), one
#' row per unit with its flattened `gid` (Figure 1b, figure MU); validity and
#' simplification are left to the caller.
read_gadm_units <- function(isos) {
  iso_sql  <- paste0("'", paste(unique(isos), collapse = "','"), "'")
  read_lvl <- function(layer, col)
    sf::st_read(FIG1_GADM_PATH, quiet = TRUE,
                query = sprintf('SELECT * FROM "%s" WHERE GID_0 IN (%s)', layer, iso_sql)) %>%
      dplyr::transmute(gid = flattenGid(as.character(.data[[col]])))
  dplyr::bind_rows(read_lvl("ADM_1", "GID_1"), read_lvl("ADM_2", "GID_2"))
}

#' The ratio colourbar of Figures 2a, 4b and AA: log10 scale capped at
#' FIG2_MAP_LO / FIG2_MAP_HI, ticks 1/3 .. 3. `...` goes to
#' ps_diverging_fill_scale() (key width, title position, legend order).
ps_ratio_fill_scale <- function(name, ...)
  ps_diverging_fill_scale(limits = log10(c(FIG2_MAP_LO, FIG2_MAP_HI)), breaks = log10(c(1 / 3, 1 / 2, 1, 2, 3)),
                          labels = c("≤ 1/3", "1/2", "1", "2", "≥ 3"), name = name, ...)

#' The IPCC-class hatch fills of the attribution ratio maps (Figure 4b, figure
#' AE): single hatch = no robust signal, cross-hatch = conflicting. NULL, with
#' a message naming `what`, when ggplot2 has no pattern fills (< 3.5).
ipcc_hatch_fills <- function(what) {
  if (!pattern_fill_supported()) {
    message("[10] ", what, ": ggplot2 < 3.5, no pattern fills -- robust ADM1s outlined instead")
    return(NULL)
  }
  list(none = mask_hatch(spacing_mm = 0.8, lwd = 0.3, col = "grey35"),
       conflicting = mask_hatch(spacing_mm = 0.8, lwd = 0.3, col = "grey35", cross = TRUE))
}

#' The layers drawn over a map's fill (as taxdat::output_plot_map()): rivers,
#' lakes, then the black country borders. `basemap`: load_africa_basemap().
map_overlay_layers <- function(basemap, afro_sf)
  list(geom_sf(data = basemap$rivers_sf, colour = taxdat:::color_rivers(), linewidth = 0.1),
       geom_sf(data = basemap$lakes_sf, fill = taxdat:::color_lake_fill(),
               colour = taxdat:::color_lake_border(), linewidth = 0.06),
       geom_sf(data = afro_sf, fill = NA, colour = "black", linewidth = FIG2_COUNTRY_BORDER_LWD))

#' Figure 2's 1000 km scale bar and north arrow, bottom-left. `data`: the facet
#' the pair is drawn in (NULL = every panel); `arrow_pad_x_mm`: the arrow's
#' left pad (0.65 mm centres it over the bar).
map_scale_layers <- function(data = NULL, arrow_pad_x_mm = 0.65)
  list(ggspatial::annotation_scale(data = data, location = "bl", width_hint = 0.15,
                                   height = unit(0.6, "mm"), text_cex = 0.35,
                                   line_width = 0.3, text_family = FONT,
                                   pad_x = unit(1, "mm"), pad_y = unit(3, "mm")),
       ggspatial::annotation_north_arrow(data = data, location = "bl",
                                         height = unit(5, "mm"), width = unit(4, "mm"),
                                         pad_x = unit(arrow_pad_x_mm, "mm"), pad_y = unit(6.5, "mm"),
                                         style = ggspatial::north_arrow_minimal(text_size = 5,
                                                                                text_family = FONT)))

#' Track names as shown in Figures 1 and 2 ("Outbreak cessation": the duration
#' track models the cessation hazard).
FIG_TRACK_NAMES <- c(occurrence = "Outbreak occurrence", duration = "Outbreak cessation",
                     monthly = "Monthly incidence")
#' `trk` as a factor of FIG_TRACK_NAMES. `erf`: "... outbreak occurrence" (the
#' end of Figure 1d's "Effect size on rate of..." axis title); `wrap`: last
#' word on its own line.
track_factor <- function(trk, erf = FALSE, wrap = FALSE) {
  lab <- unname(FIG_TRACK_NAMES)
  if (erf)  lab <- paste("...", tolower(lab))
  if (wrap) lab <- sub(" (\\S+)$", "\n\\1", lab)
  factor(trk, levels = names(FIG_TRACK_NAMES), labels = lab)
}

#' The artefact's ENSO-phase runs (NOAA episodes on the observed ONI extended
#' by the CFSv2 ensemble mean) clipped to `x_limits`; phase as a factor in
#' PS_ENSO_PHASE_COL order.
enso_phase_rects <- function(oni, x_limits) {
  # The shading's definition (caption text) from the artefact (oni v4+ field
  # `phase_rule`); artefacts without it fall back to the
  # legacy wording. The legend title stays PS_ENSO_STATE_TITLE.
  message("[10] ENSO-state shading: ", enso_shading_caption(oni))
  oni$phases %>%
    dplyr::filter(TR >= x_limits[1], TL <= x_limits[2]) %>%
    dplyr::mutate(TL = pmax(TL, x_limits[1]), TR = pmin(TR, x_limits[2]),
                  phase = factor(phase, levels = names(PS_ENSO_PHASE_COL)))
}

#' Caption text of the ENSO-state shading: the oni artefact's `phase_rule`
#' label, else the legacy wording.
enso_shading_caption <- function(oni)
  oni$phase_rule$label %||%
    "NOAA CPC episodes on the observed ONI extended by the CFSv2 ensemble mean"

#' A blank left facet strip (outside the axis, transparent text), so patchwork
#' aligns this plot's strip column with a neighbour's left strips instead of
#' inserting an empty column between its axis and panel.
blank_left_strip <- function(p, size = BASE_SIZE) {
  p + facet_grid(rows = vars("x"), switch = "y") +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 90, face = "bold", size = size, colour = "transparent"))
}

# ===========================================================================
# Section 9a: Figure 1 builders (in utils/figures_manuscript_fig1.R) and
# constants — ENSO exposure-response overview (stacked)
# ===========================================================================

# Figure 1 right column (c / d).
# Columns = ENSO teleconnection group (c carries the group header).
#   c  one tile row per group: the occurrence track's ELPD-implied lag weight
#      (pseudo-BMA over the lag grid, from 09's lag_rank) as the group colour's
#      opacity on a log scale (FIG1_W_*), weights below FIG1_W_LIM[1] white;
#      over the rainfall-correlation tiles of
#   c' (untagged, part of c) the rainfall anomaly-ONI correlation per rainfall
#      month and lag (FIG1_CORR_*), the rainy seasons outlined, and right of
#      each group's tiles the occurrence seasonal baseline per calendar month
#      (09's occ_seasonality; viridis on a log scale log-symmetric round 1,
#      FIG1_ONSET_*). The selected occurrence lag is a white dot-dashed line in
#      both (FIG1_BEST_*); the dashed rule on lag 0 (FIG1_LEAD_RULE) splits
#      leads from lags, named under the leftmost group's axis (FIG1_KEY_*).
#   d  the ERF grid (rows = track) on the ENSO-phase shading.
#   FIG1_A_TICK_PT       panel a's tick labels (countries, years, ONI).
#   FIG1_TAG_POS / FIG1_TAG_GUTTER_MM   c / d tags at the far left of each
#     sub-plot, in a gutter (shared by patchwork's alignment) left of c's labels.
#   FIG1_LAG_BREAKS / FIG1_LEAD_RULE    lag axis; the dashed rule on lag 0.
#   FIG1_X_PAD / FIG1_ONSET_GAP / FIG1_ONSET_W / FIG1_ONSET_LABEL_PAD   lag
#     units: the frame beyond the tiles, the onset column's centre right of the
#     last tile, its width, the room right of it for its rotated label.
#   FIG1_FRAME / FIG1_W_FRAME_LWD       the frames, drawn as rectangles in data
#     coordinates (c's weight strip thinner).
#   FIG1_HEADER*         c's nested header strip over the group names, with a
#     rule; FIG1_STRIP_GAP_MM above the strip text.
#   FIG1_CORR_*          the correlation tiles: 09's rain_enso corr table
#     (Pearson r of the group-mean 3-month rainfall deviation with ONI per
#     centre month and lag; method in results.derivation.rain_enso); `*`
#     p < 0.05 and >= half the group's countries significant with the same
#     sign, `.` p < 0.05 for the group mean only.
#   FIG1_ERF_*           d's y axis (Figure 2's fraction labels, up to 10).
#   FIG1_RIGHT_HEIGHTS   the weight strip, the correlation tiles, the legend
#     block (mm; replaced by its measured height), d.
#   FIG1_E_TOP_MM / FIG1_E_BOTTOM_MM   d's margins; the bottom one aligns d's
#     panel bottom with a's (measured on the render; re-measure when the
#     layout changes).
FIG1_A_TICK_PT      <- BASE_SIZE - 2
FIG1_TAG_POS        <- c(0, 1)
FIG1_TAG_GUTTER_MM  <- 2.2
FIG1_LAG_BREAKS     <- c(-6, 0, 6)
FIG1_LEAD_RULE      <- 0
FIG1_BAR_MM         <- 16
FIG1_HEADER         <- GROUP_KEY_TITLE
FIG1_HEADER_PT      <- BASE_SIZE
FIG1_STRIP_GAP_MM   <- 1.2
FIG1_X_PAD          <- 1
FIG1_FRAME          <- list(colour = "grey25", linewidth = 0.3)
FIG1_BEST_COL       <- "white"
FIG1_BEST_LWD       <- 0.2
FIG1_BEST_LTY       <- "dotdash"
# c's weight strip: opacity linear in log10(weight) over FIG1_W_LIM, squished.
FIG1_W_LIM          <- c(0.05, 1)
FIG1_W_BREAKS       <- c(0.05, 0.1, 0.25, 0.5, 1)
FIG1_W_TITLE        <- "ENSO lag selection\n(ELPD weight, log scale)"
FIG1_W_KEY_MM       <- 3
FIG1_W_FRAME_LWD    <- 0.15
FIG1_LEGEND_GAP_MM  <- 1            # between the weight legend row and the tiles' row
# The onset-season column.
FIG1_ONSET_GAP      <- 3
FIG1_ONSET_W        <- 2
FIG1_ONSET_LABEL_PAD <- 3.2
FIG1_ONSET_LABEL    <- "Onset season"
FIG1_ONSET_LABEL_PT <- BASE_SIZE - 2.5
FIG1_ONSET_TITLE    <- "Seasonal onset rate\n(× annual mean)"
FIG1_ONSET_OPTION   <- "viridis"
FIG1_ONSET_LIM      <- c(0.75, 1 / 0.75)   # log-symmetric: 1 at the palette's middle
FIG1_ONSET_BREAKS   <- c(0.75, 1, 1.33)
FIG1_ONSET_LWD      <- 0.15
# The lead / lag key under the leftmost group: "cholera leads ENSO | ENSO
# leads cholera" on one line, the bar on the lead rule, arrows below (y in
# month rows, the bottom row = 1; x in lag units).
FIG1_KEY_TEXT       <- c(left = "cholera leads ENSO", right = "ENSO leads cholera")
FIG1_KEY_PT         <- BASE_SIZE - 2.5
FIG1_KEY_Y          <- -1.35
FIG1_KEY_LINE_SKIP  <- 0.75
FIG1_KEY_BAR_GAP    <- 0.9
FIG1_KEY_ARROW_LEN  <- 12
FIG1_KEY_ARROW_LWD  <- 0.15
FIG1_KEY_HEAD_MM    <- 0.6
FIG1_CORR_COL       <- c(low = "#8c510a", mid = "white", high = "#01665e")
FIG1_CORR_LIM       <- c(-0.7, 0.7)
FIG1_CORR_TITLE     <- "Rainfall anomaly–ONI\ncorrelation (r)"
FIG1_MARK_PT        <- BASE_SIZE - 3.5
FIG1_XTITLE_GAP_MM  <- 0.4
FIG1_X_TITLE        <- "ENSO lag (months)"
# Figure 1c: main rainy season(s) per ENSO group, outlined over the rainfall
# months, from 09's rain_enso seasons (utils/rain_enso.R::rainSeasonsByGroup():
# a month is in a group's season when >= season_share of the group's population
# lives in ADM1s where that month is above the ADM1's average month); 10 reads
# the population-weighted rows.
FIG1_RAIN_SEASON_WEIGHTING <- "pop"
FIG1_RAIN_SEASON_COL <- "black"
FIG1_RAIN_SEASON_LWD <- 0.35
FIG1_RAIN_SEASON_LABEL <- "Rainy season"
FIG1_ERF_YLIM       <- c(1/3, 10)
FIG1_ERF_BREAKS     <- c(1/3, 1/2, 1, 2, 5, 10)
FIG1_ERF_LABELS     <- c("1/3", "1/2", "1", "2", "5", "10")
FIG1_ONI_BREAKS     <- c(-2, 0, 2)
FIG1_E_TOP_MM       <- 0.3
FIG1_E_BOTTOM_MM    <- 2.24
FIG1_RIGHT_HEIGHTS  <- grid::unit(c(2.2, 28.3, 6.6, 74.6), c("null", "null", "mm", "null"))

# Axis / strip text of Figure 1c and 1d.
FIG1_CD_THEME <- theme(axis.text = element_text(size = BASE_SIZE - 1.5, colour = "grey20"),
                       axis.title = element_text(size = BASE_SIZE - 1),
                       strip.text.x = element_text(size = BASE_SIZE - 1.5, margin = margin(b = 0.6, unit = "mm")))

# Strips, legends and tags shared by Figure 1c and 1d.
FIG1_PANEL_THEME <- theme(strip.placement = "outside", strip.background = element_blank(),
                          strip.text.x = element_text(size = BASE_SIZE - 1, face = "plain",
                                                      margin = margin(b = 0.6, unit = "mm")),
                          strip.text.y.left = element_text(angle = 90, size = BASE_SIZE - 1.5),
                          strip.switch.pad.grid = unit(0.5, "mm"),
                          legend.position = "bottom",
                          legend.title = element_text(size = BASE_SIZE - 1.5),
                          legend.text = element_text(size = BASE_SIZE - 2),
                          legend.key.height = unit(1.6, "mm"),
                          plot.tag.position = FIG1_TAG_POS)

# Figure 1c, seasonality design (builders fig1c_*() and build_fig1() in
# utils/figures_manuscript_fig1.R; the c above stays for build_fig1_v91()). Inputs:
# 09's oni_composite, rain_enso,
# seasonality (fitted curves, v2: with the occurrence band) and lag_rank.
#   FIG1C_TRACKS / _SEAS_TRACK / _TRACK_LABELS  the track blocks (10's ids),
#     the seasonality artefact's track names, the names left of each block.
#   FIG1C_PEAK_X / _PHASES / _X_BREAKS  the 2-year axis (x 1..24 = Jan Y ..
#     Dec Y+1): the lag anchor (Dec Y, lag L at x = 12 + L), the Rising / Peak /
#     Decay boxes, the month ticks.
#   FIG1C_REF            d's dashed reference style (fig1_erf_panel()): the
#     peak rule, ONI 0, x annual mean 1.
#   FIG1C_FONT           c's text at d's sizes.
#   FIG1C_KEY_*          the same / following-year key (first track only).
#   FIG1C_VKEY_*         d's vertical "more / fewer" key, leftmost column.
#   FIG1C_HEIGHTS        row heights (relative): the ONI row at 0.5, the lag
#     strips shortened (keyed x0.85, key-less x0.7) in favour of the
#     seasonality rows; c's legend row and d.
FIG1C_TRACKS        <- c("occurrence", "duration", "monthly")
FIG1C_SEAS_TRACK    <- c(occurrence = "occurrence", duration = "cessation", monthly = "cases")
FIG1C_TRACK_LABELS  <- stats::setNames(sub(" (\\S+)$", "\n\\1", unname(FIG_TRACK_NAMES[FIG1C_TRACKS])), FIG1C_TRACKS)
FIG1C_PEAK_X        <- 12
FIG1C_PHASES        <- tibble::tribble(
  ~phase,   ~from, ~to,   # x, inclusive
  "Rising",     6,  10,   # Jun-Oct Y
  "Peak",      11,  13,   # Nov Y - Jan Y+1
  "Decay",     14,  18)   # Feb-Jun Y+1
FIG1C_PHASE_FILL    <- c(Rising = "grey93", Peak = "grey78", Decay = "grey93")
FIG1C_PHASE_RULE    <- list(colour = "grey45", linewidth = 0.12, linetype = "22")
FIG1C_REF           <- list(colour = "grey40", linewidth = 0.3, linetype = "dashed")
FIG1C_LWD           <- list(oni = 0.3, seas = 0.3, frame = 0.1, border = 0.2)
FIG1C_X_BREAKS      <- c(1, 7, 13, 19)
FIG1C_X_TITLE       <- "Month (El Niño year, following year)"
FIG1C_YEAR_LABS     <- c("Y", "Y+1")
FIG1C_ONI_COL       <- c(elnino = "#a87b3a", nina = "#2f3e4c")
FIG1C_ONI_LAB       <- c(elnino = "El Niño", nina = "La Niña")
FIG1C_RAIN_PALETTE  <- "Blues"
FIG1C_RAIN_TITLE    <- "Mean monthly rainfall (mm)"
FIG1C_FONT          <- list(axis_text = BASE_SIZE - 1.5, row_title = BASE_SIZE - 1.5, track = BASE_SIZE - 1,
                            annot = BASE_SIZE - 2.5, legend_title = BASE_SIZE - 1.5, legend_text = BASE_SIZE - 2,
                            annot_col = "grey25", lineheight = 0.85)
FIG1C_SEAS_YLAB     <- "× annual\nmean"
FIG1C_SEAS_LIM      <- list(occurrence = c(1/3, 3), duration = c(1/4, 4), monthly = c(1/3, 3))
FIG1C_SEAS_BREAKS   <- list(occurrence = c(1/3, 1, 3), duration = c(1/4, 1, 4), monthly = c(1/3, 1, 3))
FIG1C_SEAS_LABELS   <- list(occurrence = c("1/3", "1", "3"), duration = c("1/4", "1", "4"),
                            monthly = c("1/3", "1", "3"))
FIG1C_SEAS_BAND_ALPHA <- 0.25
FIG1C_KEY_ARROW     <- list(len = 4, head_mm = 0.9, col = "grey25", lwd = 0.25, y = 2.0, text_y = 2.25,
                            gap = 1.5, start = 1.5)
FIG1C_KEY_TEXT      <- c(left = "same\nyear", right = "following\nyear")
FIG1C_W_YTOP        <- c(key = 3.1, nokey = 1.8)
FIG1C_VKEY_TEXT     <- list(occurrence = c(up = "more\nonsets", down = "fewer\nonsets"),
                            duration   = c(up = "shorter\noutbreaks", down = "longer\noutbreaks"),
                            monthly    = c(up = "more\ncases", down = "fewer\ncases"))
FIG1C_VKEY_SPAN     <- c(from = 0.07, to = 0.9)
FIG1C_VKEY_X        <- 1.1
FIG1C_VKEY_DX       <- 1.25
FIG1C_SPACING_MM    <- 1.2
FIG1C_LABEL_W       <- 0.075
FIG1C_TRACK_SHIFT   <- 0.5
FIG1C_TAG_X         <- -0.085
FIG1C_HEIGHT_MM     <- 174
FIG1C_HEIGHTS <- local({
  h <- list(oni = 0.7, rain = 0.08, w = 0.38, w_nokey = 0.38 * (1.8 - 0.5) / (3.1 - 0.5), seas = 0.33)
  rest <- h$rain + h$w + 2 * h$w_nokey + 3 * h$seas
  k <- (rest + h$oni - 0.5) / rest
  h <- list(oni = 0.5, rain = h$rain * k, w = h$w * k, w_nokey = h$w_nokey * k, seas = h$seas * k)
  freed <- h$w * (1 - 0.85) + 2 * h$w_nokey * (1 - 0.7)
  c(modifyList(h, list(w = h$w * 0.85, w_nokey = h$w_nokey * 0.7, seas = h$seas + freed / 3)),
    list(legend = 0.4, d = 1.9))
})

# ===========================================================================
# Section 9b: Figure 2 — ADM1 IRR maps (stacked) + the group-map inset;
# constants in Section 0, builders in utils/figures_manuscript_fig2.R
# ===========================================================================
# ===========================================================================
# Section 9c: Figure 3 constants — attribution of the past strong ENSO events;
# builders, incl. the forecast panels Figure 4 reuses (fig3_oni_strip(),
# fig3_panel_b() / fig3_panel_c(), f3_*), in utils/figures_manuscript_fig34.R
# ===========================================================================
# Colours, line types and shapes. The Africa-wide arm is bordeaux, regions use
# GROUP_COLOURS, every counterfactual one neutral dark grey (dashed); the
# free-running simulation is dot-dashed; masked smoothing dark goldenrod (so it
# does not share the counterfactual's grey); the counterfactual's point is a
# square, since the open circle is (a)'s observed symbol.
F3_AFRICA_COL <- PS_AFRICA_COL   # alias (Section 2b)
F3_CF_COL     <- "grey35"
F3_MASK_COL   <- "#B8860B"
F3_ONI_COL    <- "grey10"         # CFSv2 forecast ONI (dot-dashed, as every free-running series)
F3_ARM_LTY    <- c(forecast = "dotdash", counterfactual = "22")
F3_ARM_SHAPE  <- c(forecast = 16, counterfactual = 15)
# Figure 4's attribution factual (09's attribution$attrib_factual == "attrib":
# the factual launched at the counterfactual handoff, results.inputs.attrib_manifest):
# its median as a dotted line in b, and the factual of every paired statistic.
F4_ATTRIB_LTY <- "dotted"
F3_LWD        <- list(a = 0.35, b = 0.4, c = c(0.25, 0.6))
# y caps (cases only; bands squished at the cap): Figure 3a (`a`), Figure 4's Africa
# row, its region rows (shared scale).
F3_CASES_YMAX <- c(a = 30000, africa = 45000, region = 30000)
# c's arm padding (larger = forecast and counterfactual closer).
F3_C_ARM_EXPAND  <- 1.2
# Compact legend theme for the cell above c.
F3_LEGEND_THEME <- theme(legend.position = "bottom", legend.justification = "left",
                         legend.box.just = "left",
                         legend.text = element_text(size = BASE_SIZE - 1),
                         legend.title = element_text(size = BASE_SIZE - 1),
                         legend.key.height = unit(2.6, "mm"), legend.key.width = unit(5.5, "mm"),
                         legend.spacing.y = unit(0.2, "mm"), legend.spacing.x = unit(2.5, "mm"),
                         legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0))

# --- Figure 3: attribution of the past strong ENSO events ---
# Event table, phase colours and readings. The event panels read the free
# factual simulation under the headlined explosive mode
# (results.style.attrib_headline.explosive: "drop" = explosive pairs dropped,
# "none" = every pair kept); b's interval is the noise-adjusted
# posterior one (09's event_af_posterior, its median as the bar), c's maps the
# AF of the replicate means.
HL_EXPLOSIVE <- as.character(RES$style$attrib_headline[["explosive"]] %||% "drop")
if (length(HL_EXPLOSIVE) != 1L || !HL_EXPLOSIVE %in% c("drop", "none") ||
    !HL_EXPLOSIVE %in% as.character(unlist(RES$derivation$explosive_modes %||% c("drop", "none"))))
  stop("[10] results.style.attrib_headline$explosive must be one of results.derivation.explosive_modes {",
       paste(unlist(RES$derivation$explosive_modes), collapse = ", "), "}; got '", HL_EXPLOSIVE, "'.", call. = FALSE)
F3_EV_PHASE_COL <- c(elnino = PS_ENSO_PHASE_COL[["elnino"]], lanina = PS_ENSO_PHASE_COL[["nina"]])   # bracket / column-title colours
F3_EV_FACTUAL   <- "simulation"
F3_EV_GUARD     <- identical(HL_EXPLOSIVE, "drop")   # the tables' explosive_guard of the headlined mode
F3_EV_AF_LIM    <- c(-100, 100)   # c's fill limits (%), squished beyond
F3_EV_B_XLIM    <- c(-50, 100)    # b's x range (%), bars clamped and intervals arrowed beyond
F3_EV_MASK_COL  <- "grey45"       # coverage-masked smoothing in a (neutral, not the NGA+GHA yellow)
F3_EV_MAP_PAD_DEG <- 0.5          # c's crop: the event's estimated ADM1s plus this margin (degrees)
# Layout: the event-bracket and ONI strips in mm; panel a, a's legend row, b, c
# and c's legend row as relative heights of the rest of F3_EV_HEIGHT_MM (fixed
# mm rows exclude axes / labels, so a page sized from them clips). The first
# column's left margin; the legend rows are centred on the plot area, i.e.
# WIDTH_MM minus the left axis gutter F3_EV_GUTTER_MM, their legends
# F3_EV_LEG_GAP_MM apart.
F3_EV_STRIP_MM       <- c(strip = 11, oni = PS_ONI_STRIP_MM)
F3_EV_HEIGHTS        <- c(a = 1.35, leg_a = 0.14, b = 1.15, c = 0.78, leg_bc = 0.16)
F3_EV_HEIGHT_MM      <- 182
F3_EV_LEFT_MARGIN_MM <- 0.3
F3_EV_GUTTER_MM      <- 14.5
F3_EV_LEG_GAP_MM     <- 6
# The group-map inset (fig. S10a's partition without labels or dagger, as
# Figure 2's: fig2_group_map_inset()) in b's F3_EV_INSET_COL-th event column,
# in the box F3_EV_INSET_BOX (npc of that column's panel: its bottom-left
# corner, clear of the bars). NULL = no inset.
F3_EV_INSET_COL <- 1L
F3_EV_INSET_BOX <- c(left = 0.01, bottom = 0.01, right = 0.30, top = 0.45)

# ===========================================================================
# Section 9d: Figure 4 constants — the 2026 ENSO forecast: forecast vs the
# neutral-ENSO counterfactual and its ENSO-attributable fraction (the paired
# ratio + ADM1 probability maps Figure 4 is in 10_manuscript_figures_legacy.R);
# builders in utils/figures_manuscript_fig34.R
# ===========================================================================
#   a  Figure 3's forecast panels: the ONI forecast strip over Africa + the
#      ENSO regions from the analogue window (forecast vs counterfactual,
#      fig3_panel_b()), the cumulative totals (fig3_panel_c()'s column for the
#      metric) and the forecast / counterfactual ratio of means, 1 / (1 - AF)
#      from af_posterior, filled by the scope's IPCC class at the headline m
#      (report_attrib_metrics v2 ipcc_scope).
#   b  the AE map (attrib_ratio_map_layers(), utils/figures_supp.R): ADM1 ratio
#      of means, IPCC hatching at the headline m, Figure 2's colour bar.
#   c  ENSO-attributable fraction AF = 1 - CF / F of Africa and each ENSO region
#      (af_posterior: noise-adjusted posterior median, 50% / 95%).
# Main analog set, the headlined attribution cell and explosive mode (09's
# headline slice). Cases in the main figure, outbreak onsets supplementary.
# Layout: one flat grid -- columns = a's time series | cumulative | ratio |
# b / c; rows = the ONI strip (mm), a's upper / lower part (b ends, c starts at
# their border), the legend row (arm keys under a, group keys under c).
FIG4_WIDTHS    <- c(ts = 1.15, cum = 0.6, ratio = 0.6, right = 2.6)
FIG4_HEIGHTS   <- c(upper = 1.32, lower = 0.88, legend = 0.24)   # c taller
FIG4_HEIGHT_MM <- 165
# The group-map inset in c (as Figure 3's, F3_EV_INSET_BOX): npc box of c's
# panel, its lower-left, left of the GH-GS intervals. NULL = no inset.
FIG4_INSET_BOX <- c(left = 0.02, bottom = 0.03, right = 0.27, top = 0.47)
# Ratio axis (no 1.5: the dashed 1 + m line marks it); cumulative totals one
# break per decade; c's AF axis range (%), rounded out to 25 within it.
FIG4_RATIO_BRK <- c(0.25, 0.5, 1, 2, 4, 8, 16)
FIG4_RATIO_LAB <- c("1/4", "1/2", "1", "2", "4", "8", "16")
FIG4_CUM_BRK   <- 10^(0:6)
FIG4_AF_XLIM   <- F3_EV_B_XLIM   # Figure 3b's cut
# IPCC classes of a's ratio bars: robust in dark grey (red would read as
# Africa's colour).
FIG4_CLASS_COLS <- c(robust = "grey30", conflicting = "#f4a582", none = "grey78")
FIG4_CLASS_LAB  <- c(robust = "Robust", conflicting = "Conflicting", none = "No robust signal")
# The headlined attribution cell (results.style.attrib_headline.window; 09's
# attrib_cells()): "attribution" = each group over its own firm window (the
# legacy mixed Africa), "attribution_common" = every group over the common
# Africa window of results.scope.africa.
FIG4_WINDOW     <- as.character(RES$style$attrib_headline[["window"]] %||% "attribution")
# Continuation months (a group read past its firm horizon on the analog ENSO
# continuation): faded under a white veil of this opacity, labelled with the
# group's continuation lead.
FIG4_CONT_VEIL  <- 0.55
# b's margins: top = its tag level with a's (the map is a wrapped element, not
# aligned by the grid); bottom = the white space between its legends and c.
FIG4_MAP_TOP_MM <- 2
FIG4_C_GAP_MM   <- 5

# ===========================================================================
# Section 9e: supplementary figures
# ===========================================================================
#' Supplementary — Figure 1 ERF at country granularity: population fixed effect
#' plus that country's own ENSO-slope RE, one facet per country, coloured by
#' ENSO group. Same deployed best-lag fit and x-axis as Figure 1's ERF grid.
# y-axis label per track: the ratio each track's ERF is on (duration: the
# cessation hazard ratio, < 1 = longer outbreaks).
SUPP_ERF_YLAB <- c(occurrence = "Onset hazard ratio", duration = "Cessation hazard ratio",
                   monthly = "Monthly case incidence rate ratio")
build_fig1_supp_country_erf <- function(track = "occurrence") {
  ce <- read_results("enso_relation", track = track)$country_erf %>%
    dplyr::mutate(enso_group = factor(enso_group, levels = GROUP_ORDER)) %>%
    dplyr::arrange(enso_group, country) %>%
    dplyr::mutate(country = factor(country, levels = unique(country)))
  n_countries <- length(levels(ce$country))
  fig <- ggplot(ce, aes(x = at, y = IRR, ymin = IRR_lo, ymax = IRR_hi,
                        colour = enso_group, fill = enso_group)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_ribbon(alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.5) +
    # "CIV (0.32)": the country's p_max in the deployed sweep cell (report_membership).
    facet_wrap(~ country, ncol = 6, labeller = ggplot2::as_labeller(function(x) {
      p <- membership_pmax()
      if (is.null(p)) return(x)
      ifelse(x %in% names(p), sprintf("%s (%.2f)%s", x, p[x], ifelse(p[x] < MEMBERSHIP_LOW, MEMBERSHIP_MARK, "")), x)
    })) +
    scale_colour_manual(values = GROUP_COLOURS, name = GROUP_KEY_TITLE, labels = GROUP_LABELS) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_y_log10() +
    # integer ONI ticks, horizontal (the 45-degree labels were clipped at the
    # facet edges); no title / subtitle, the caption carries the description
    scale_x_continuous(breaks = -2:2) +
    coord_cartesian(xlim = ERF_XLIM) +
    labs(x = "Oceanic Ni\u00f1o Index (ONI)", y = SUPP_ERF_YLAB[[track]]) +
    theme_ms() +
    theme(legend.position = "bottom",
          strip.text  = element_text(size = BASE_SIZE - 1),
          plot.margin = margin(1.5, 3, 1.5, 1.5, unit = "mm"))
  stem <- if (identical(track, "occurrence")) "fig1_supp_country_erf" else sprintf("fig1_supp_country_erf_%s", track)
  save_figure(fig, stem,
              height_mm = ceiling(n_countries / 6L) * 26 + 25,
              width_mm = 210, subdir = "supplementary")
}

#' Supplementary — region/country parameter caterpillar for one track. Reads the
#' derive_track_params() artefact (09) and plots, for every NON-ENSO parameter
#' family in SUPP_PARAM_FAMILIES[[track]], the country-level posterior ratio
#' (exp scale, 95% CrI) with the per-group population estimate overlaid as a
#' distinct marker. One panel per term (SUPP_PARAM_TERM_LAB, wrapped six to a
#' row); rows are the countries under their ENSO group's pooled ("population")
#' estimate, groups in GROUP_ORDER; colour encodes the group, shape the tier.
#' Terms with a country RE are panel a (every country under its group's pooled
#' row); terms without one (duration bands, monthly seasonal harmonics and
#' population elasticity) are panel b, one row per group, so they do not take
#' a page of empty country rows. One float per track; thin rules separate the
#' groups.
#' Readable labels of the model terms (raw brms names -> display); a term not
#' listed prints its raw name; panels follow this order. Reference levels: time since the last outbreak
#' < 1 year; month 1 of an outbreak.
SUPP_PARAM_TERM_LAB <- c(
  Intercept = "Intercept",
  ongoing_outbreak_lag1TRUE = "Outbreak ongoing in the previous month",
  nbr_lag1_r500_rn_parthmmclusters = "Neighbour pressure (previous month)",
  time_to_last_outbreak12 = "Last outbreak 1-2 years ago",
  time_to_last_outbreak23 = "Last outbreak 2-3 years ago",
  time_to_last_outbreak34 = "Last outbreak 3-4 years ago",
  time_to_last_outbreak45 = "Last outbreak 4-5 years ago",
  time_to_last_outbreak56 = "Last outbreak 5-6 years ago",
  time_to_last_outbreak6Inf = "Last outbreak \u2265 6 years ago or none",
  mso_band2 = "Month 2 of the outbreak", mso_band3 = "Month 3 of the outbreak",
  mso_band4M5 = "Months 4-5 of the outbreak", mso_band6M8 = "Months 6-8 of the outbreak",
  mso_band9P = "Month 9 or later of the outbreak",
  sin_month = "Season, sine (annual)", cos_month = "Season, cosine (annual)",
  sin2_month = "Season, sine (semi-annual)", cos2_month = "Season, cosine (semi-annual)",
  mso_cb1 = "Months since onset, basis 1", mso_cb2 = "Months since onset, basis 2",
  mso_cb3 = "Months since onset, basis 3", mso_cb4 = "Months since onset, basis 4",
  log_pop_c = "Population elasticity (log population, centred)",
  shape_Intercept = "Negative-binomial shape (intercept)")
build_fig_supp_track_params <- function(track = "occurrence") {
  fams <- SUPP_PARAM_FAMILIES[[track]]
  if (is.null(fams)) { warning(sprintf("[10] no SUPP_PARAM_FAMILIES for track %s", track)); return(invisible(NULL)) }
  d <- read_results("track_params", track = track) %>% dplyr::filter(family %in% fams)
  if (!nrow(d)) { warning(sprintf("[10] track params (%s): nothing to plot after family filter", track)); return(invisible(NULL)) }

  term_tbl <- d %>% dplyr::distinct(family, term) %>%
    dplyr::mutate(family = factor(family, levels = fams)) %>%
    dplyr::arrange(match(term, names(SUPP_PARAM_TERM_LAB)), family, term)
  term_lab <- ifelse(term_tbl$term %in% names(SUPP_PARAM_TERM_LAB),
                     unname(SUPP_PARAM_TERM_LAB[term_tbl$term]), term_tbl$term)
  pd <- d %>%
    dplyr::mutate(ylab = ifelse(tier == "population", sprintf("%s pooled", GROUP_SHORT_FIG[enso_group]), level),
                  enso_group = factor(enso_group, levels = GROUP_ORDER),
                  term = factor(term, levels = term_tbl$term, labels = term_lab))
  varying <- unique(as.character(pd$term[pd$tier == "country"]))
  pd_a <- dplyr::filter(pd, term %in% varying)
  pd_b <- dplyr::filter(pd, !term %in% varying, tier == "population")

  # One block of facets; rows top to bottom: per group (GROUP_ORDER) its pooled
  # estimate, then its countries alphabetically.
  param_panel <- function(pdx, ncol, tag) {
    rows <- pdx %>% dplyr::distinct(enso_group, tier, ylab) %>%
      dplyr::mutate(pop = tier == "population") %>%
      dplyr::arrange(enso_group, dplyr::desc(pop), ylab)
    lev <- rev(unique(rows$ylab))
    pdx$ylab <- factor(pdx$ylab, levels = lev)
    # Rules between groups: below the last row of every group but the last.
    last <- rows %>% dplyr::group_by(enso_group) %>% dplyr::slice_tail(n = 1) %>% dplyr::ungroup()
    rules <- match(utils::head(last$ylab, -1L), lev) - 0.5
    ggplot(pdx, aes(x = mean, y = ylab, colour = enso_group, shape = tier)) +
      geom_hline(yintercept = rules, colour = "grey80", linewidth = 0.25) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
      geom_errorbar(aes(xmin = q025, xmax = q975), width = 0, linewidth = 0.4, orientation = "y") +
      geom_point(size = 1.3) +
      facet_wrap(~ term, ncol = ncol, scales = "free_x", drop = TRUE,
                 labeller = ggplot2::label_wrap_gen(width = 22)) +
      scale_x_log10(breaks = scales::breaks_log(4), labels = scales::label_number(drop0trailing = TRUE)) +
      scale_y_discrete(drop = TRUE) +
      # Explicit legend order (group key, then tier): unset, ggplot2's tie-break
      # depended on the figures drawn before this one in the same session.
      scale_colour_manual(values = GROUP_COLOURS, breaks = GROUP_ORDER, labels = GROUP_LABELS[GROUP_ORDER],
                          name = GROUP_KEY_TITLE, drop = FALSE, guide = guide_legend(order = 1)) +
      scale_shape_manual(values = c(country = 16, population = 18), name = NULL, drop = FALSE,
                         limits = c("country", "population"),
                         labels = c(country = "Country", population = "ENSO group (pooled)"),
                         guide = guide_legend(order = 2)) +
      labs(x = "Posterior ratio (exp scale) and 95% CrI", y = NULL, tag = tag) +
      theme_ms() + TAG_THEME +
      theme(legend.position = "bottom",
            strip.text   = element_text(size = BASE_SIZE - 0.5),
            axis.text.y  = element_text(size = BASE_SIZE - 1.5),
            axis.text.x  = element_text(size = BASE_SIZE - 1),
            panel.grid.major.y = element_blank(),
            panel.spacing.x = unit(3, "mm"))
  }
  # Panel heights in mm: strip + axis per facet row, a pitch per data row.
  block_mm <- function(pdx, ncol, pitch) {
    if (!nrow(pdx)) return(0)
    ceiling(dplyr::n_distinct(pdx$term) / ncol) *
      (14 + dplyr::n_distinct(paste(pdx$enso_group, pdx$ylab)) * pitch)
  }
  # Up to six terms on one facet row, else five per row. With one or two
  # panel-a terms (duration), b sits beside a, at the top of a's column.
  n_a    <- dplyr::n_distinct(pd_a$term)
  ncol_a <- if (n_a <= 6L) n_a else 5L
  has_b  <- nrow(pd_b) > 0
  beside <- has_b && n_a <= 2L
  ncol_b <- if (beside) 3L else 5L
  pa <- param_panel(pd_a, ncol_a, if (has_b) "a" else NULL)
  # b has no country rows: a's legend is the figure's only one.
  pb <- if (has_b) param_panel(pd_b, ncol_b, "b") + guides(colour = "none", shape = "none")
  h_a <- block_mm(pd_a, ncol_a, 2.6)
  h_b <- block_mm(pd_b, ncol_b, 4)
  fig <- if (beside) {
    patchwork::wrap_plots(A = pa, B = pb, design = "AB\nA#",
                          widths = c(n_a + 0.6, ncol_b), heights = c(h_b, max(h_a - h_b, 1)),
                          guides = "collect") & theme(legend.position = "bottom")
  } else if (has_b) {
    patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(h_a, h_b),
                          guides = "collect") & theme(legend.position = "bottom")
  } else pa
  stem <- sprintf("figS_params_%s", track)
  save_figure(fig, stem, height_mm = 14 + if (beside) max(h_a, h_b) else h_a + h_b,
              width_mm = WIDTH_MM, subdir = "supplementary")
}

#' All requested analog sets' ENSO forecast input on one panel (vs. the
#' single-set ONI strip of Figure 4a). Only meaningful when >1 set is on hand.
build_fig_oni_comparison <- function() {
  message("[10] Supplementary: ONI envelope comparison (all analog sets)")
  oni <- read_results("oni", sets = SETS)
  key_levels <- c("Observed", "Forecast", set_lab(SETS))
  env <- oni$env %>%
    dplyr::mutate(key = factor(ifelse(arm == "forecast", "Forecast", set_lab(analog_set)),
                               levels = key_levels))
  # Observed ONI in the legend too: one colour key (the ribbons share their
  # line's colour, so they carry no key of their own).
  pal <- c(Observed = "black", Forecast = FC_COLOR, stats::setNames(unname(SET_COLORS[SETS]), set_lab(SETS)))
  fig <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(data = env, aes(date, ymin = lo, ymax = hi, fill = key), alpha = 0.18) +
    geom_line(data = env, aes(date, mean, colour = key), linewidth = 0.6) +
    geom_line(data = oni$obs, aes(date, enso_index, colour = "Observed"), linewidth = 0.45) +
    # Both arms' launches when they differ (oni v4 `handoffs`), else one line.
    geom_vline(xintercept = unique(as.Date(c(oni$handoff, unlist(oni$handoffs)), origin = "1970-01-01")),
               linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    scale_colour_manual(values = pal, breaks = key_levels, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
    labs(x = NULL, y = "Oceanic Ni\u00f1o Index (ONI)") +
    theme_ms()
  save_figure(fig, "fig4_supp_oni_envelope_comparison", height_mm = 80, width_mm = 200, subdir = "supplementary")
}

# ===========================================================================
# Section 9f: Figure H — hindcast forecast skill
# ===========================================================================
#' Reads the four phase-09 hindcast artefacts; every panel is a shared builder
#' from utils/figures_hindcast.R (the notebook draws the same plots), except
#' the supplementary skill-vs-climatology figure (utils/figures_supp.R).
# Display floor of that figure's CRPSS axis (intervals below are clipped and marked).
HS_YLIM_LO <- -1
# Colours of the climatology references and the case supports (keys as the
# hindcast artefacts store them): blues / reddish purple / greys, none of them a
# group colour (Dark2) nor an ENSO-phase colour; the headline ones in blue.
HC_REF_COLOURS <- c(loo_gid = "#0072B2", regional = "#56B4E9", persistence = "#CC79A7",
                    insample_gid = "grey60", loo_gid_fullhist = "#00446b", zero = "grey30",
                    loo_climatology = "grey30")
HC_SUPPORT_COLOURS <- c(atrisk_gids = "#0072B2", obs_positive_keys = "#56B4E9", footprint = "#CC79A7",
                        exact_gid = "grey60", atrisk_keys = "grey35")
build_fig_hindcast <- function() {
  message("[10] Figure H: hindcast forecast skill")
  runs <- read_results("hindcast_runs"); occ <- read_results("hindcast_occ")
  cas  <- read_results("hindcast_cases"); rep <- read_results("report_hindcast")
  hc <- SCOPE$hindcast
  groups <- intersect(GROUP_ORDER, unique(runs$group))
  # Africa first, then the groups in GROUP_ORDER (the order of every other
  # figure; hc_scope_labels() puts Africa last)
  lab <- hc_scope_labels(GROUP_SHORT_FIG, groups)[c("africa", groups)]
  pal <- hc_scope_palette(GROUP_COLOURS, groups)[c("africa", groups)]
  th  <- theme_ms() + ggplot2::theme(legend.position = "bottom")
  refs <- intersect(c(hc$headline_occ_reference, "insample_gid", "regional", "persistence"),
                    occ$meta$references)

  pa <- hc_plot_skill_bars(occ$skill, refs, hc$headline_occ_reference, "window", "crps_fair",
                           c("all", rep$headline$short_lead), lab, th, pal = HC_REF_COLOURS, legend_nrow = 2) +
    ggplot2::labs(title = "Outbreak-count skill vs climatology references")
  pb <- hc_plot_skill_by_lead(occ$skill, hc$headline_occ_reference, "window", "crps_fair",
                              occ$meta$lead_levels, lab, pal, th) +
    ggplot2::labs(title = "Skill by lead") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(.98, .98),
                   legend.justification = c(1, 1))
  pc <- hc_plot_case_support(cas$calibration, cas$skill,
                             intersect(c(hc$headline_case_support, "obs_positive_keys", "footprint"),
                                       cas$meta$supports),
                             hc$headline_case_support, labels = lab, theme = th, pal = HC_SUPPORT_COLOURS) +
    ggplot2::labs(title = "Case totals by support")
  pd <- hc_plot_reliability(occ$reliability, lab, th) + ggplot2::labs(title = "Reliability")
  fig <- (pa | pb) / (pc | pd) +
    patchwork::plot_layout(heights = c(1, 1.4)) +
    patchwork::plot_annotation(tag_levels = "a") & TAG_THEME
  save_figure(fig, "figH_hindcast_skill", height_mm = 230, subdir = "hindcast")

  if (SUPP) {
    supp <- function(p, stem, h) save_figure(p, stem, height_mm = h, subdir = "hindcast")
    supp(hc_plot_skill_bars(occ$skill, occ$meta$references, hc$headline_occ_reference,
                            "window", "crps_fair", occ$meta$lead_levels, lab, th, pal = HC_REF_COLOURS),
         "figH_supp_all_references", 110)
    supp(hc_plot_skill_bars(occ$skill, c("insample_gid", hc$headline_occ_reference),
                            hc$headline_occ_reference, "calendar_year", "crps_emp", "all",
                            lab, th, pal = HC_REF_COLOURS) + ggplot2::labs(title = "Calendar-year target"),
         "figH_supp_legacy_basis", 70)
    # Histograms: ticks at 0 / 0.5 / 1 only (the facets' end labels collided);
    # each axis names its scored quantity (a and c: the randomised PIT of the
    # monthly onset indicator per unit; b: the outbreak count per unit and window).
    unit_x <- ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
    supp((hc_plot_pit_hist(occ$pit_hist, lab, th) + unit_x +
            ggplot2::labs(x = "PIT of the monthly onset indicator per unit") |
            hc_plot_rank_hist(occ$rank_hist, "window", "all", lab, th) + unit_x +
            ggplot2::labs(x = "Normalised rank of the observed outbreak count per unit and window")) /
           (hc_plot_coverage_curve(occ$coverage_curve, lab, pal, th) +
              ggplot2::labs(x = "Nominal central interval width",
                            y = "Empirical coverage of the monthly onset indicator\n(per unit, from its PIT)") |
              hc_plot_roc(occ$roc, dplyr::filter(occ$cell_scores, lead_class == "all",
                                                 reference == hc$headline_occ_reference),
                          lab, pal, th)) +
           patchwork::plot_annotation(tag_levels = "a") & TAG_THEME,
         "figH_supp_calibration", 200)
    supp(hc_plot_country_lollipop(cas$skill, hc$headline_case_support, theme = th),
         "figH_supp_country_cases", 120)
    supp(hc_plot_cases_interval(cas$run_totals, hc$headline_case_support, "region", "all", lab, th),
         "figH_supp_case_intervals", 120)
    # CRPSS against the headline climatologies (occurrence: headline_occ_reference;
    # monthly cases: headline_case_skill_reference on the headline support) by
    # group and lead class, with the season-year bootstrap intervals
    # (figures_supp.R; display floor HS_YLIM_LO).
    supp(build_figS_hindcast_clim_skill(
           hindcast_clim_skill_frame(occ$skill, cas$month_skill, hc$headline_occ_reference,
                                     hc$headline_case_skill_reference, hc$headline_case_support,
                                     lead_levels = c(setdiff(occ$meta$lead_levels, "all"), "all")),
           lab[c("africa", groups)], pal[c("africa", groups)], theme_ms(), ylim_lo = HS_YLIM_LO,
           base_size = BASE_SIZE),
         "figH_supp_skill_vs_climatology", 170)
    # Skill against climatology (a) over case totals with prediction intervals
    # (b) in one figure (the two above combined).
    supp(patchwork::wrap_plots(
           build_figS_hindcast_clim_skill(
             hindcast_clim_skill_frame(occ$skill, cas$month_skill, hc$headline_occ_reference,
                                       hc$headline_case_skill_reference, hc$headline_case_support,
                                       lead_levels = c(setdiff(occ$meta$lead_levels, "all"), "all")),
             lab[c("africa", groups)], pal[c("africa", groups)], theme_ms(), ylim_lo = HS_YLIM_LO,
             base_size = BASE_SIZE),
           hc_plot_cases_interval(cas$run_totals, hc$headline_case_support, "region", "all", lab, th),
           ncol = 1, heights = c(1.3, 1)) +
           patchwork::plot_annotation(tag_levels = "a") & TAG_THEME,
         "figH_supp_hindcast_skill", 230)
  }
}

# ===========================================================================
# Section 10: figure steps, one fig_step_<key>() per FIGURE_STEPS entry
# (Section 11 runs the requested ones; nothing here runs at definition)
# ===========================================================================
fig_step_1 <- function() {
  build_fig1(STACKED_TRACKS)
  if (SUPP) {
    for (trk in FIG1_TRACKS) build_fig1_supp_country_erf(trk)
    # The region/country parameter caterpillars of every non-ENSO term (figs
    # S7-S9): the deployed fits' parameter figures go with Figure 1's.
    for (trk in FIG1_TRACKS) build_fig_supp_track_params(trk)
  }
}

fig_step_2 <- function() {
  build_fig2(STACKED_TRACKS)
}

fig_step_3 <- function() {
  # Cases (main) and onsets (supplementary); the event runs have one
  # counterfactual set (the event manifest's), so no per-set variants.
  build_fig3(c("cases", "onsets"))
}

fig_step_4 <- function() {
  # The 2026 forecast: cases (main), outbreak onsets (supplementary), main
  # analog set; the forecast ONI envelopes of every analog set (supplementary).
  build_fig4("cases", MAIN_SET)
  build_fig4("onsets", MAIN_SET)
  if (SUPP && length(SETS) > 1L) build_fig_oni_comparison()
}

fig_step_H <- function() {
  if (isTRUE(SCOPE$hindcast$enable)) build_fig_hindcast()
  else message("[10] figure H requested but results.scope.hindcast.enable is false; skipped")
}

# Fig. S5: sampler diagnostics of the deployed models (convergence artefacts).
fig_step_S5 <- function() {
  rc   <- read_results("report_convergence")
  summ <- dplyr::bind_rows(lapply(FIG1_TRACKS, function(t) read_results("convergence", track = t)$summary))
  # Threshold lines at the deploy convergence gate's values (rc$gate: R-hat
  # 1.05, ESS 400, the check the deployed fits passed), not the stricter
  # DIAG_THRESHOLDS screening values (rc$thresholds, R-hat 1.01).
  thr <- c(rhat = rc$gate$rhat_max %||% rc$thresholds[["rhat"]], ess = rc$gate$ess_min %||% rc$thresholds[["ess"]])
  save_figure(build_figS5_convergence(rc$fits, summ, thr, GROUP_COLOURS,
                                      GROUP_SHORT_FIG, theme_ms()),
              "figS5_convergence", height_mm = 150, subdir = "supplementary")
}

# The lag_averaging artefacts present for FIG1_TRACKS, or NULL:
# the step is optional (results.derivation.lag_averaging.enable) and needs the
# 05 lag sidecars.
read_lag_averaging <- function() {
  arts <- lapply(stats::setNames(FIG1_TRACKS, FIG1_TRACKS), function(t) {
    p <- makeResultsRunFile(opt, "results_lag_averaging", track = t, schema = RESULTS_SCHEMA[["lag_averaging"]],
                            config_path = CONFIG_PATH)
    if (file.exists(p)) read_results("lag_averaging", track = t) else NULL
  })
  arts <- Filter(Negate(is.null), arts)
  if (!length(arts)) return(NULL)
  lapply(stats::setNames(c("lags", "erf", "irr", "roles"), c("lags", "erf", "irr", "roles")),
         function(f) dplyr::bind_rows(lapply(arts, `[[`, f)))
}

# Fig. S6: LOO lag rankings (report_lags; the ranking Figure 1c outlines), with
# the deploy floor, the rank-2 lag and (lag_averaging) the high Pareto-k lags.
fig_step_S6 <- function() {
  rl <- read_results("report_lags")
  la <- read_lag_averaging()
  save_figure(build_figS6_lag_loo(rl$ranking, rl$selected, GROUP_COLOURS, GROUP_SHORT_FIG,
                                  theme_ms(), k_flags = la$lags),
              "figS6_lag_loo", height_mm = 120, subdir = "supplementary")
}

# Forecast vs counterfactual by country (attribution$country_ribbon; the
# 07d_simulator_figures_counterfactual.R fig6/fig7 series), main analog set.
fig_step_FC <- function() {
  f23 <- read_results("attribution", sets = SETS)
  n_ct <- length(unique(f23$country_ribbon$country))
  vl <- unique(c(as.Date(f23$attrib_start), as.Date(f23$factual_start)))
  for (m in c("cases", "onsets"))
    save_figure(build_figS_forecast_country(f23$country_ribbon, m, MAIN_SET, GROUP_ORDER, GROUP_COLOURS,
                                            GROUP_SHORT_FIG, F3_CF_COL, theme_ms(), vlines = vl),
                sprintf("figS_forecast_country_%s", m), height_mm = 22 + 32 * ceiling(n_ct / 5),
                subdir = "supplementary")
}

# Historical reconstruction by country (smoothing$by_country; the
# 08c_smoothing_compile.R fig3 country panels).
fig_step_SC <- function() {
  sm <- read_results("smoothing")
  n_ct <- length(sm$country_group)
  for (m in c("cases", "onsets"))
    save_figure(build_figS_smoothing_country(sm$by_country, sm$country_group, m, GROUP_ORDER,
                                             GROUP_SHORT_FIG, GROUP_COLOURS, theme_ms(),
                                             start = as.Date(RES$derivation$smoothing_start)),
                sprintf("figS_smoothing_country_%s", m), height_mm = 22 + 32 * ceiling(n_ct / 5),
                subdir = "supplementary")
}

# Population in ADM1s with more ENSO-attributable cases (report_attrib_exposure):
# ratio-of-means map with the IPCC robustness hatching at the headline m
# (style.attrib_headline) | population shares by class, every m. Main analog set.
fig_step_AE <- function() {
  ae <- read_results("report_attrib_exposure", sets = SETS)
  hl_m <- as.numeric(RES$style$attrib_headline$magnitude)
  if (length(hl_m) != 1L || !any(abs(ae$magnitudes - hl_m) < 1e-9))
    stop("[10] figure AE: results.style.attrib_headline$magnitude must be one of thresholds.attrib_magnitudes {",
         paste(ae$magnitudes, collapse = ", "), "}.", call. = FALSE)
  maps <- load_africa_maps("Figure AE")
  basemap <- africa_basemap(maps$afro_sf)
  hatch <- ipcc_hatch_fills("figure AE")
  grey_sf <- ps_country_status(maps$afro_sf, modelled_countries())
  for (m in c("cases", "onsets")) {
    # A fresh grey layer per figure: ggnewscale renames a layer's fill in place
    # on first use, so a reused new_scale_fill() layer breaks the second plot.
    grey <- ps_country_grey_layer(grey_sf, legend_order = 3, labels = PS_GREY_LAB_WRAP)
    # No scale bar: the map's legends take its bottom-left corner (as Figure 4).
    fig <- build_figS_attrib_exposure(
             ae$adm1, relab_groups(ae$exposure, c("scope", "region")), m, MAIN_SET, hl_m, maps$adm1_sf,
             sf::st_bbox(maps$afro_sf),
             base_layers = grey, overlay_layers = map_overlay_layers(basemap, maps$afro_sf),
             limits = c(FIG2_MAP_LO, FIG2_MAP_HI),
             scope_order = display_scope_levels(),
             class_cols = c(robust = "#b2182b", conflicting = "#f4a582", childs = "black"),
             hatch_fills = hatch, map_theme = taxdat::map_theme() + ps_map_font(),
             theme = theme_ms(), base_size = BASE_SIZE, window = FIG4_WINDOW, guard = F3_EV_GUARD) +
      patchwork::plot_annotation(tag_levels = "a") & TAG_THEME
    save_figure(fig, sprintf("figS_attrib_exposure_%s", m), height_mm = 88, subdir = "supplementary")
  }
}

# Attributable fraction of the past strong ENSO events (report_event_attrib +
# report_event_af_model): simulated (both factual readings) vs analytic, per
# event and ENSO group.
fig_step_EV <- function() {
  if (!nzchar(RES$event_manifest %||% ""))
    stop("[10] figure EV: results.event_manifest is empty (no event attribution derived).", call. = FALSE)
  ea  <- read_results("report_event_attrib")
  eam <- read_results("report_event_af_model")
  # The simulation reading's interval is event_af_posterior's (noise-adjusted
  # 2.5% / 97.5% over draws), as in 11's tables; report_event_attrib carries
  # the paired bootstrap only (empty when attrib_boot_b = 0). The smoothing
  # reading has no posterior and keeps its own (bootstrap) columns. A scope
  # whose posterior spread is unidentified has NA adjusted
  # bounds, so it shows its point without an interval, never a zero-width one.
  post <- read_results("event_af_posterior")$scope %>%
    dplyr::transmute(event, explosive_guard, level, scope, metric, factual_source = "simulation",
                     p_lo = af_adj_lo, p_hi = af_adj_hi)
  ea$af_scope <- ea$af_scope %>%
    dplyr::left_join(post, by = c("event", "explosive_guard", "level", "scope", "metric", "factual_source")) %>%
    dplyr::mutate(af_lo = dplyr::coalesce(p_lo, af_lo), af_hi = dplyr::coalesce(p_hi, af_hi)) %>%
    dplyr::select(-p_lo, -p_hi)
  # The stored region short names -> display labels (style.group_labels).
  ea$af_scope <- relab_groups(ea$af_scope, "region"); eam$scope <- relab_groups(eam$scope, "region")
  for (m in c("cases", "onsets"))
    save_figure(build_figS_event_attrib(ea, eam, m, display_scope_levels(),
                                        theme_ms(), base_size = BASE_SIZE,
                                        tag = c(cases = "a", onsets = "b")[[m]], guard = F3_EV_GUARD) + TAG_THEME,
                sprintf("figS_event_attrib_%s", m), height_mm = 75, subdir = "supplementary")
  # Per-winter sub-windows of a multi-winter event (report_event_attrib v3):
  # the same panels, one per (event, winter).
  sv <- event_sub_window_view(ea, eam)
  if (!is.null(sv)) {
    sv$ea$af_scope <- relab_groups(sv$ea$af_scope, "region"); sv$eam$scope <- relab_groups(sv$eam$scope, "region")
    for (m in c("cases", "onsets"))
      save_figure(build_figS_event_attrib(sv$ea, sv$eam, m, display_scope_levels(),
                                          theme_ms(), base_size = BASE_SIZE,
                                          tag = c(cases = "a", onsets = "b")[[m]], guard = F3_EV_GUARD) + TAG_THEME,
                  sprintf("figS_event_attrib_sub_windows_%s", m), height_mm = 75, subdir = "supplementary")
  }
}

# Membership uncertainty of the deployed classes (report_membership): the
# responsibility matrix and the class paths across the sweep cells.
fig_step_MB <- function() {
  mb <- membership()
  if (is.null(mb))
    stop("[10] figure MB: results.kl_sweep_dir / results.kl_sweep_cell are unset (no report_membership).",
         call. = FALSE)
  save_figure(build_figS_membership(mb, GROUP_COLOURS, GROUP_SHORT_FIG, theme_ms(), low = MEMBERSHIP_LOW,
                                    base_size = BASE_SIZE),
              "figS_membership", height_mm = 140, subdir = "supplementary")
}

# Supplement maps of the analysis set and the grouping.
# MU: every unit of the survival panel by admin level, its units with no
# onset, and the ADM1s of the modelled countries holding no modelled unit (no
# recorded episode), all classified in 09 (modelled_units); polygons from the
# GADM ADM_1 / ADM_2 layers as Figure 1b.
fig_step_MU <- function() {
  st <- read_results("modelled_units")
  rate <- read_results("onset_rate")
  poly <- read_gadm_units(st$units$country) %>%
    sf::st_make_valid() %>% rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE)
  maps <- load_africa_maps("Figure MU")
  unit_sf <- dplyr::inner_join(poly, st$units, by = "gid")
  miss <- setdiff(st$units$gid, unit_sf$gid)
  if (length(miss))
    message(sprintf("[10] figure MU: %d unit(s) without a GADM polygon: %s", length(miss),
                    paste(utils::head(miss, 10), collapse = ",")))
  gap_sf <- dplyr::inner_join(dplyr::select(maps$adm1_sf, adm1 = GID_1), st$gaps, by = "adm1") %>%
    simplify_map_sf()
  grey <- ps_country_grey_layer(ps_country_status(maps$afro_sf, modelled_countries(rate)),
                                legend_order = 2, labels = PS_GREY_LAB_WRAP, show_legend = FALSE)
  basemap <- africa_basemap(maps$afro_sf)
  save_figure(build_figS_modelled_units(unit_sf, gap_sf, maps$afro_sf, sf::st_bbox(maps$afro_sf),
                                        base_layers = grey, overlay_layers = map_overlay_layers(basemap, maps$afro_sf),
                                        map_theme = taxdat::map_theme() + ps_map_font(), base_size = BASE_SIZE),
              "figS_modelled_units", height_mm = 170, subdir = "supplementary")
  message(sprintf("[10] figure MU: %d units (%s), %d ADM1s without a recorded episode",
                  nrow(st$units), paste(sprintf("%s %d", names(table(st$units$status)), as.integer(table(st$units$status))),
                                        collapse = ", "), nrow(st$gaps)))
}

# GM: the deployed partition (the fit registry's countries by group), labelled
# by the display short names (style.group_labels; GW..GS under "centroid");
# countries with p_max below style.membership_low marked with Figure 1a's
# dagger (MEMBERSHIP_MARK), no outline (as Figure 1b).
fig_step_GM <- function() {
  members <- read_results("modelled_units")$members   # the fit registry's partition (09)
  maps <- load_africa_maps("Figure GM")
  grey <- ps_country_grey_layer(ps_country_status(maps$afro_sf, members$iso3), legend_order = 2,
                                labels = PS_GREY_LAB_WRAP, show_legend = FALSE)
  save_figure(build_figS_group_map(maps$afro_sf, members, GROUP_COLOURS, GROUP_LABELS[GROUP_ORDER],
                                   GROUP_SHORT_FIG[GROUP_ORDER], sf::st_bbox(maps$afro_sf),
                                   low_isos = membership_low_isos(), low = MEMBERSHIP_LOW, mark = MEMBERSHIP_MARK,
                                   base_layers = grey,
                                   overlay_layers = list(geom_sf(data = maps$afro_sf, fill = NA, colour = "white",
                                                                 linewidth = 0.1)),
                                   map_theme = taxdat::map_theme() + ps_map_font(), base_size = BASE_SIZE),
              "figS_group_map", height_mm = 150, subdir = "supplementary")
}

# KC: the class sweep's MAP partitions for K = 2..5 at the pinned cell's basis
# and prior (report_kl_sweep$partition; agreeing chains), countries ordered by
# their deployed group.
fig_step_KC <- function() {
  cell <- as.character(RES[["kl_sweep_cell"]] %||% "")
  if (!nzchar(RES$kl_sweep_dir %||% "") || !nzchar(cell))
    stop("[10] figure KC: results.kl_sweep_dir / results.kl_sweep_cell are unset (no report_kl_sweep).", call. = FALSE)
  kl <- read_results("report_kl_sweep")
  d <- kl_classes_by_K_frame(kl$partition, kl$fits, cell, 2:5)
  mem <- read_results("modelled_units")$members      # deployed group per country (09)
  dep <- stats::setNames(mem$group, mem$iso3)
  save_figure(build_figS_kl_classes_by_K(d, dep, GROUP_COLOURS, GROUP_SHORT_FIG[GROUP_ORDER], theme_ms(),
                                         base_size = BASE_SIZE),
              "figS_kl_classes_by_K", height_mm = 30 + 4 * dplyr::n_distinct(d$iso3), width_mm = 120,
              subdir = "supplementary")
}

# EG: figs. GM, KC and MB panel a in one figure, and MB panel b on its own.
fig_step_EG <- function() {
  mb <- membership()
  cell <- as.character(RES[["kl_sweep_cell"]] %||% "")
  if (is.null(mb) || !nzchar(RES$kl_sweep_dir %||% "") || !nzchar(cell))
    stop("[10] figure EG: results.kl_sweep_dir / results.kl_sweep_cell are unset (no report_membership / report_kl_sweep).",
         call. = FALSE)
  members <- read_results("modelled_units")$members   # the fit registry's partition (09)
  maps <- load_africa_maps("Figure EG")
  grey <- ps_country_grey_layer(ps_country_status(maps$afro_sf, members$iso3), legend_order = 2,
                                labels = PS_GREY_LAB_WRAP, show_legend = FALSE)
  map_plot <- build_figS_group_map(maps$afro_sf, members, GROUP_COLOURS, GROUP_LABELS[GROUP_ORDER],
                                   GROUP_SHORT_FIG[GROUP_ORDER], sf::st_bbox(maps$afro_sf),
                                   low_isos = membership_low_isos(), low = MEMBERSHIP_LOW, mark = MEMBERSHIP_MARK,
                                   base_layers = grey,
                                   overlay_layers = list(geom_sf(data = maps$afro_sf, fill = NA, colour = "white",
                                                                 linewidth = 0.1)),
                                   map_theme = taxdat::map_theme() + ps_map_font(), base_size = BASE_SIZE)
  kl <- read_results("report_kl_sweep")
  d <- kl_classes_by_K_frame(kl$partition, kl$fits, cell, 2:4)
  dep <- stats::setNames(members$group, members$iso3)
  save_figure(build_figS_enso_groups(map_plot, mb, d, dep, GROUP_COLOURS, GROUP_SHORT_FIG[GROUP_ORDER], theme_ms(),
                                     low = MEMBERSHIP_LOW, base_size = BASE_SIZE) & TAG_THEME,
              "figS_enso_groups", height_mm = 150, subdir = "supplementary")
  # The same with the class paths under the map (panel d).
  save_figure(build_figS_enso_groups(map_plot, mb, d, dep, GROUP_COLOURS, GROUP_SHORT_FIG[GROUP_ORDER], theme_ms(),
                                     low = MEMBERSHIP_LOW, base_size = BASE_SIZE, paths = TRUE) & TAG_THEME,
              "figS_enso_groups_paths", height_mm = 210, subdir = "supplementary")
  save_figure(build_figS_membership_paths(mb, GROUP_COLOURS, GROUP_SHORT_FIG, theme_ms(), base_size = BASE_SIZE),
              "figS_membership_paths", height_mm = 140, subdir = "supplementary")
}

# DC: outbreak data coverage (report_dataset$outbreaks, one row per analysis
# outbreak): episodes per country and cholera season year, countries in ENSO
# group blocks (display labels), and episodes per season year by admin level.
fig_step_DC <- function() {
  fr <- data_coverage_frames(read_results("report_dataset")$outbreaks)
  save_figure(build_figS_data_coverage(fr, GROUP_ORDER, GROUP_SHORT_FIG[GROUP_ORDER], theme_ms(),
                                       country_lab = FIG1_COUNTRY_NAMES, base_size = BASE_SIZE),
              "figS_data_coverage", height_mm = 45 + 3.6 * dplyr::n_distinct(fr$country$country),
              subdir = "supplementary")
}

# SE: seasonality against the rainy season (09 seasonality): columns = ENSO
# groups, rows = rainfall climatology, onset rate, cessation rate, cases per
# active episode-month; the rainy-season band is Fig. 1d's months.
SE_RAIN_FILL <- "#dfe7ee"
SE_RAIN_COL  <- "#6b7f91"
fig_step_SE <- function() {
  art <- read_results("seasonality")
  save_figure(build_figS_seasonality(seasonality_frames(art, GROUP_ORDER), GROUP_SHORT_FIG[GROUP_ORDER],
                                     GROUP_COLOURS, theme_ms(), rain_fill = SE_RAIN_FILL, rain_col = SE_RAIN_COL,
                                     min_n = art$meta$min_n, base_size = BASE_SIZE),
              "figS_seasonality", height_mm = 150, subdir = "supplementary")
}

# WR: observed vs modelled ENSO ratio (09 report_who_enso_ratio), one figure
# per derived group (a suffix only when there are several): panel a the same
# calendar months as the model (WHO monthly), panel b the target year to the
# WHO cut-off against full neutral years (lower bound).
fig_step_WR <- function() {
  art <- read_results("report_who_enso_ratio")
  sty <- STYLE[["who_enso_ratio"]]   # resolved with get_default_options()
  ydef <- as.character(sty[["year_definition"]])
  aset <- as.character(sty[["analog_set"]])
  if (length(ydef) != 1L || length(aset) != 1L)
    stop("[10] figure WR: results.style.who_enso_ratio must set year_definition and analog_set.", call. = FALSE)
  m <- art$meta
  win <- as.Date(m$window)
  win_lab <- sprintf("%s-%s %d", format(min(win), "%b"), format(max(win), "%b"), as.integer(m$target_year))
  for (g in m$groups) {
    d <- who_enso_ratio_frame(art$comparison, g, ydef, aset)
    ny <- function(e) gsub(" ", ", ", d$neutral_years[d$estimator == e & d$kind == "group"][1] %||% "")
    titles <- c(same_window = sprintf("%s, neutral years %s", win_lab, ny("same_window")),
                ytd_annual = sprintf("%d to %s (lower bound) vs full neutral years", as.integer(m$target_year),
                                     format(as.Date(m$cutoff) + 6L, "%d %b")))
    p <- build_figS_who_enso_ratio(d, unname(GROUP_SHORT_FIG[[g]]), unname(GROUP_COLOURS[[g]]), theme_ms(), titles,
                                   model_lab = sprintf("Modelled: factual / %s, %s", ps_cf_label(aset), win_lab),
                                   obs_lab = sprintf("Observed ratio: %d / mean of neutral years", as.integer(m$target_year)),
                                   base_size = BASE_SIZE, family = FONT)
    save_figure(p, if (length(m$groups) > 1L) paste0("figS_who_enso_ratio_", g) else "figS_who_enso_ratio",
                height_mm = 80, subdir = "supplementary")
  }
}

# PC: the processing cascade (09 report_dataset$cascade, also Data S1 D1) as a
# flow diagram, the descriptive notebook's layout; replaces the long table S1.
fig_step_PC <- function() {
  save_figure(build_figS_processing_cascade(processing_cascade_frames(read_results("report_dataset")$cascade),
                                            base_size = BASE_SIZE + 2, family = FONT),
              "figS_processing_cascade", height_mm = 170, subdir = "supplementary")
}

# SV: the historical reconstruction against the observations (smoothing v4
# $validation, 09): PIT histograms re-binned to SV_PIT_BINS, onset seasonality,
# episode durations up to SV_MAX_DURATION months, active fraction per scope,
# onset frequency given an ongoing outbreak nearby. Display constants only.
SV_PIT_BINS     <- 20L
SV_MAX_DURATION <- 24L
fig_step_SV <- function() {
  v <- read_results("smoothing")$validation
  if (is.null(v)) stop("[10] figure SV: the smoothing artefact has no validation tables (smoothing v4).", call. = FALSE)
  for (k in c("pit_hist", "seasonal", "duration", "active", "contagion")) v[[k]] <- relab_groups(v[[k]], "scope")
  save_figure(build_figS_smoothing_validation(v, display_scope_levels(),
                                              theme_ms(), scope = "Africa", pit_bins = SV_PIT_BINS,
                                              max_duration = SV_MAX_DURATION, base_size = BASE_SIZE),
              "figS_smoothing_validation", height_mm = 200, subdir = "supplementary")
}

# FA: the attribution factual against the counterfactual of each analogue set
# (attribution artefact ribbons: attrib_ribbon / attrib_africa_ribbon, else the
# forecast rows, and the counterfactual rows of region_ribbon / africa_ribbon),
# monthly onsets and cases, Africa + each group, from the attribution start to
# the end of the headlined cell (the common Africa window when it is one).
# Every set the config derived (results.scope.supp_sets); all four when the
# supplementary sets are on.
fig_step_FA <- function() {
  f23 <- read_results("attribution", sets = SETS)
  if (length(SETS) < length(ENSO_ANALOG_SETS))
    message("[10] figure FA: only the analogue set(s) ", paste(SETS, collapse = ", "),
            " were derived (results.scope.supplementary / supp_sets).")
  groups <- GROUP_ORDER[GROUP_ORDER %in% f23$groups]
  fa_end <- if (grepl("_common$", FIG4_WINDOW) && identical(f23$africa_window$mode, "common"))
    as.Date(f23$africa_window$end) else NULL
  d_fa <- forecast_analogue_frame(f23, groups, SETS, end = fa_end)
  attrib <- identical(f23$attrib_factual %||% "forecast", "attrib") && !is.null(f23$attrib_ribbon)
  arm_lab <- c(factual = if (attrib) sprintf("Factual launched %s (attribution)",
                                             format(as.Date(f23$attrib_arm_starts$factual %||% f23$attrib_start), "%b %Y"))
                         else "ENSO-informed forecast",
               counterfactual = "Counterfactual (analogue set of the column)")
  save_figure(build_figS_forecast_analogues(d_fa, c(Africa = "Africa", GROUP_SHORT_FIG[groups]), SET_LABELS[SETS],
                                            theme_ms(), cf_col = F3_CF_COL,
                                            scope_col = c(Africa = F3_AFRICA_COL, GROUP_COLOURS[groups]), arm_lab = arm_lab,
                                            base_size = BASE_SIZE),
              "figS_forecast_analogues", height_mm = 60 + 2 * 14 * (length(groups) + 1), subdir = "supplementary")
  # Overlay variant: the factual once and the counterfactual
  # median of every set on one axis, onsets (a) beside cases (b).
  save_figure(build_figS_forecast_analogues_overlay(d_fa, c(Africa = "Africa", GROUP_SHORT_FIG[groups]),
                                                    SET_LABELS[SETS], SET_COLORS[SETS], theme_ms(),
                                                    fac_lab = arm_lab[["factual"]], base_size = BASE_SIZE),
              "figS_forecast_analogues_overlay", height_mm = 40 + 26 * (length(groups) + 1),
              subdir = "supplementary")
}

# AA: attribution maps by analogue set (report_attrib_exposure$adm1 at the
# headlined cell and explosive mode, as Figure 4b): per metric, the ADM1 ratio
# of mean factual to mean counterfactual totals (Figure 4b's log scale and cap,
# FIG2_MAP_LO / HI) and the probability that the factual exceeds the
# counterfactual over the paired replicates (p_sign), one row per set.
fig_step_AA <- function() {
  ae <- read_results("report_attrib_exposure", sets = SETS)
  if (length(SETS) < length(ENSO_ANALOG_SETS))
    message("[10] figure AA: only the analogue set(s) ", paste(SETS, collapse = ", "),
            " were derived (results.scope.supplementary / supp_sets).")
  maps <- load_africa_maps("Figure AA")
  # Full-detail ADM1 in 2 x length(SETS) vector panels made an 18.7 MB PDF.
  maps$adm1_sf <- simplify_map_sf(maps$adm1_sf)
  grey_sf <- ps_country_status(maps$afro_sf, modelled_countries())
  aa_maps <- lapply(c(cases = "cases", onsets = "onsets"), function(m)
    attribution_analogue_map_sf(ae$adm1, m, SETS, maps$adm1_sf, window = FIG4_WINDOW, guard = F3_EV_GUARD))
  ratio_fill <- ps_ratio_fill_scale("Factual ÷ counterfactual (mean)", key_width_mm = 28, title_position = "top",
                                    order = 1)
  prob_fill <- ps_diverging_fill_scale(
    limits = c(0, 1), midpoint = 0.5, breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "25%", "50%", "75%", "100%"),
    name = "P(factual > counterfactual)", key_width_mm = 28, title_position = "top", order = 2)
  save_figure(build_figS_attribution_analogues(
                aa_maps, SET_LABELS[SETS], sf::st_bbox(maps$afro_sf), ratio_fill, prob_fill,
                limits = c(FIG2_MAP_LO, FIG2_MAP_HI),
                base_layers_fn = function() ps_country_grey_layer(grey_sf, legend_order = 3, labels = PS_GREY_LAB_WRAP),
                overlay_layers = list(geom_sf(data = maps$afro_sf, fill = NA, colour = "black",
                                              linewidth = FIG2_COUNTRY_BORDER_LWD)),
                map_theme = taxdat::map_theme() + ps_map_font(), base_size = BASE_SIZE),
              "figS_attribution_analogues", height_mm = 45 + 44 * length(SETS), subdir = "supplementary")
}

# ===========================================================================
# Section 11: dispatch (the requested keys, in FIGURE_STEPS order)
# ===========================================================================
message(sprintf("[10] figures=%s supplementary=%s main_set=%s supp_sets=%s",
                paste(FIGURES, collapse = ","), SUPP, MAIN_SET,
                paste(SUPP_SETS, collapse = ",")))
for (k in FIGURE_KEYS[FIGURE_KEYS %in% FIGURES]) match.fun(FIGURE_STEPS[[k]])()
message("[10] done.")
