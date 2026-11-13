#!/usr/bin/env Rscript
# analysis/10_manuscript_figures_legacy.R
# =============================================================================
# LEGACY: frozen copy of an earlier analysis/10_manuscript_figures.R, kept
# under the side-by-side rule. The current figures are in 10; do not edit this
# file. Its only changes are this note, the file name above and the output
# root below (LEGACY_OUT_DIR): it always
# writes to figures/manuscript_legacy/ (results.style.out_dir is ignored),
# or to $LEGACY_FIGURES_OUT_DIR when set, and refuses any path inside
# figures/manuscript/, so it can never overwrite the current figures.
# Everything below this note is that version of 10.
# =============================================================================
# Phase 10 — MANUSCRIPT (publication) figures. Plotting only.
#
# Figures
# -------
#   1  ENSO exposure-response overview (stacked)
#        a  observed-ONI strip over the data-coverage timeline (rows ordered
#           by last case bubble per group; ENSO-phase background)
#        b  inset map: onsets per surveillance-covered month per gid
#        c  ENSO lag selection: dELPD vs the best lag per track (rows) and
#           ENSO teleconnection group (columns), selected lag outlined
#        d  rainfall anomaly-ONI correlation per rainfall month and lag, per
#           group (see the data contract)
#        e  ERF grid (rows = track, columns = group; 50% + 95% CrI;
#           direction arrows in the first column)
#   2  ADM1 IRR maps (stacked), rows = track, columns = ONI
#      condition (results.derivation.enso_conds), a direction key per row
#   3  Historical smoothing -> forecast vs neutral-ENSO counterfactual
#        a  observed-ONI strip over Africa's smoothing + free-running forecast
#        b  ONI forecast strip over Africa + the four regions from the
#           analogue handoff, forecast vs counterfactual
#        c  cumulative suspected cases and outbreak onsets over the analogue
#           window, forecast vs counterfactual, per scope
#      Main: cases; supplementary: onsets and the other analog sets.
#   4  Attribution over the ENSO-divergence window
#        a  paired forecast / counterfactual ratio forest stack, two columns
#           (suspected cases | outbreak onsets)
#        b  ADM1 signed probability of a >=50% shift, cases above, onsets below
#   H  Hindcast forecast skill (results.scope.hindcast; builders in
#      utils/figures_hindcast.R, shared with notebooks/hindcast_skill.qmd)
#        a  outbreak-count CRPSS vs the climatology references, all leads and
#           the 0-7-month class (headline reference solid)
#        b  CRPSS against the headline reference by lead class
#        c  case totals by support: sim/obs ratio, CRPSS vs LOO climatology,
#           95% coverage (headline support solid)
#        d  reliability diagram
#
# Supplementary: Figure 1 ERF at country granularity; Figures 3 / 4 for the
# other Lenssen analog sets; 4-set ONI-envelope comparison; hindcast
# calibration (PIT, rank histogram, coverage, ROC), per-country case skill and
# legacy-basis comparison.
#
# DATA CONTRACT — never derives, re-fits, or re-runs the bpfilter. Reads the
# phase-09 aggregates from generated_data/results/<tag>/ via results_require(),
# which refuses a missing, stale, or differently-decided artefact. Figure 1a
# (and the bubble check that places 1b) reads the survival / coverage / monthly
# data directly -- a narrow exception scoped to that panel; Figures 1b, 2 and 4
# read GADM / giscoR polygons as basemaps. The ONI strips and ENSO-phase
# shading of Figures 1 and 3 come from the phase-09 ONI artefact.
# Figure 1d reads a rainfall anomaly-ONI correlation table written outside the
# pipeline (FIG1_RAIN_CORR_CSV) -- a second narrow exception, with no
# freshness check.
#
# Usage
# -----
#   Rscript analysis/10_manuscript_figures.R -C analysis/configs/results/manuscript.yml
#   Rscript analysis/10_manuscript_figures.R -C ... --figures 4,H --supplementary false
#
# --figures (comma-separated subset of 1,2,3,4,H, or the supplementary S5,S6,FC,SC,AE) and --supplementary
# (true|false) override results.scope.build_figures / .supplementary for this
# run only: they pick which figures are drawn, never which artefacts are read
# or how they were derived (RES and the analog sets stay as the config says).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  # Explicit tidyverse sub-packages, in tidyverse's own core-attach order.
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(readr)
  library(purrr)
  library(dplyr)
  library(stringr)
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
# deploy_config is "".
for (.m in c("results_core", "figures_cb", "figures_io", "deploy_core", "figures_hindcast",
             "figures_supp"))
  source(file.path("analysis/utils", paste0(.m, ".R")))

# ---------------------------------------------------------------------------
# Section 0: CLI
# ---------------------------------------------------------------------------
# S5 / S6 / FC / SC / AE: supplementary (figures_supp.R); FC = forecast vs
# counterfactual by country, SC = historical reconstruction (smoothing) by
# country, AE = population in ADM1s with more ENSO-attributable cases, EV =
# attributable fraction of the past strong ENSO events (needs results.event_manifest).
FIGURE_KEYS <- c("1", "2", "3", "4", "H", "S5", "S6", "FC", "SC", "AE", "EV")
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
THRESHOLDS <- RES$thresholds

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

PROB_THRESHOLD <- THRESHOLDS$prob_threshold

FORMATS <- tolower(as.character(STYLE$formats))
DPI     <- as.integer(STYLE$dpi)
# LEGACY: the output root is fixed in code (see the header note).
LEGACY_OUT_DIR <- Sys.getenv("LEGACY_FIGURES_OUT_DIR", "figures/manuscript_legacy")
local({
  # Absolute, symlink-resolved path, also for a path that does not exist yet
  # (normalizePath() returns those unchanged): resolve the deepest existing
  # ancestor and re-append the rest.
  abs_path <- function(p) {
    p <- path.expand(p)
    if (!startsWith(p, "/")) p <- file.path(getwd(), p)
    rest <- character(0)
    while (!file.exists(p)) { rest <- c(basename(p), rest); p <- dirname(p) }
    if (any(rest %in% c(".", ".."))) stop("[10-legacy] LEGACY_FIGURES_OUT_DIR: no '.' / '..' in the part that does not exist yet")
    do.call(file.path, as.list(c(normalizePath(p), rest)))
  }
  cur <- abs_path("figures/manuscript")
  out <- abs_path(LEGACY_OUT_DIR)
  if (identical(out, cur) || startsWith(out, paste0(cur, "/")))
    stop("[10-legacy] refusing to write into figures/manuscript/ (", LEGACY_OUT_DIR, ")")
})
STYLE$out_dir <- LEGACY_OUT_DIR
OUT_DIR <- STYLE$out_dir
ERF_XLIM <- as.numeric(unlist(STYLE$erf_xlim))

# Figure 2 presentation constants (not config surface).
#   FIG2_MAP_LO/HI      — IRR cap, symmetric on the log scale; also Figure 1d's
#     ERF y-axis limits, so the two figures share one effect-size range.
#   FIG2_RULE_LWD / FIG2_RULE_TRIM — booktabs \cmidrule-style header rules
#     (ggplot linewidth, mm) under each ENSO phase and each ONI column, trimmed
#     at both ends (the phase-header fills are trimmed to the same width).
#   FIG2_ROW_H_MM       — figure height per track row (+ FIG2_LEGEND_MM for the
#     legend row under the maps).
#   FIG2_PHASE_FILL     — La Nina / El Nino header fills: Figure 1's ENSO-phase
#     colours at PS_ENSO_PHASE_ALPHA.
#   FIG2_ROW_KEY_*      — the direction key right of each track row (labels and
#     column width).
#   Map line widths: white ADM1 borders FIG2_ADM1_LWD, significant-ADM1
#     outlines FIG2_SIG_LWD, black country borders FIG2_COUNTRY_BORDER_LWD.
FIG2_MAP_LO <- 1 / 3
FIG2_MAP_HI <- 3
FIG2_COUNTRY_BORDER_LWD <- 0.1
FIG2_ADM1_LWD  <- 0.002
FIG2_SIG_LWD   <- 0.05
FIG2_RULE_LWD  <- 0.15
FIG2_RULE_TRIM <- grid::unit(3, "mm")
FIG2_ROW_H_MM  <- 38
FIG2_LEGEND_MM <- 8
FIG2_PHASE_FILL <- c(nina = scales::alpha("blue", 0.075), elnino = scales::alpha("red", 0.075),
                     neutral = scales::alpha("grey50", 0.075))
# What red / blue means on each track row (Figure 1d's arrow labels; the
# cessation IRR is inverted in the maps, so red = lower cessation hazard =
# longer outbreaks).
FIG2_ROW_KEY_LABELS <- list(
  "Outbreak occurrence" = c(up = "more\nonsets",      down = "fewer\nonsets"),
  "Outbreak cessation"  = c(up = "longer\noutbreaks",  down = "shorter\noutbreaks"),
  "Monthly incidence"   = c(up = "more\ncases",       down = "fewer\ncases"))
FIG2_ROW_KEY_MM <- 15


# Figure 1 presentation constants, not config surface, like the FIG2_* ones.
# Figure 1b inset map (panel a's top-right corner):
#   FIG1_INSET_MAX_ROWS — the inset reaches down at most this many panel-a rows,
#     fewer if a deeper box (wider, by the map's fixed aspect) would cover a
#     case bubble.
#   FIG1_INSET_MARGIN   — gap to panel a's top/right border (panel-width fraction).
#   FIG1_INSET_ASPECT   — map height / width of the inset's lon/lat box.
#   FIG1_PANEL_A_ASPECT — panel a's width / height in the 183 x 140 mm layout,
#     measured from the render (1412 x 2584 px at 600 dpi);
#     re-measure if the figure size, panel a's share of the width, the ONI strip
#     or panel a's tick / legend text changes.
# Layout: FIG1_HEIGHT_MM tall at WIDTH_MM; left (a) : right (c / d / e) columns
# FIG1_WIDTHS; FIG1_ONI_STRIP_MM the ONI strip above a; FIG1_A_RIGHT_MM the
# margin right of a (room for its last year label; no other margin between the
# columns); FIG1_A_XLIM panel a's window (every coverage period and case
# bubble, with a little room). The right column's own constants (heights,
# e's bottom margin aligning e's panel bottom with a's, ...) are in Section 9a.
FIG1_GADM_PATH      <- "data/geodata/gadm_valid.gpkg"
FIG1_INSET_BBOX     <- c(xmin = -19, xmax = 52.5, ymin = -36.5, ymax = 37.5)
FIG1_INSET_MAX_ROWS <- 9L
FIG1_INSET_MARGIN   <- 0.006
FIG1_INSET_ASPECT   <- unname((FIG1_INSET_BBOX["ymax"] - FIG1_INSET_BBOX["ymin"]) /
                              (FIG1_INSET_BBOX["xmax"] - FIG1_INSET_BBOX["xmin"]))
FIG1_PANEL_A_ASPECT <- 1412 / 2584
FIG1_HEIGHT_MM      <- 140
FIG1_WIDTHS         <- c(0.42, 0.58)
FIG1_ONI_STRIP_MM   <- 11
FIG1_A_RIGHT_MM     <- 2.5
FIG1_A_XLIM         <- as.Date(c("2009-10-01", "2025-01-31"))

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
# Forecast colour of the supplementary ONI-envelope comparison.
FC_COLOR <- "#2166ac"

# Figure 4a significance colouring -- same red/blue as Figure 4b's map scale.
DIRECTION_COLOURS <- c(increase = "#d73027", "no change" = "grey60", decrease = "#2166ac")
# Figure 4 layout: forest-stack dot size range (mm, area-scaled), columns
# (forest stack : maps), legend-row split (a's legend : the map legend) and
# height (relative to the panel row), forest row height (mm; sets the figure
# height), and the outbreak-onset column's ratio range (cases: 0.02-50).
FIG4_DOT_RANGE  <- c(0.6, 3.5)
FIG4_WIDTHS     <- c(0.58, 0.42)
FIG4_LEG_WIDTHS <- c(0.44, 0.56)
FIG4_LEG_ROW    <- 0.11
FIG4_ROW_MM     <- 3.8
FIG4_ONSET_LIM  <- c(0.05, 20)

SET_COLORS <- c(all = "#e08214", nina = "#35978f", neutral = "#878787",
                elnino = "#b2182b")
SET_LABELS <- c(all = "all analogues", nina = "La Ni\u00f1a analogues",
                neutral = "neutral analogues", elnino = "El Ni\u00f1o analogues")

# TRACK_* / GROUP_* / COMPOSITE_CONDS / METRIC_* live in utils/results_core.R;
# ERF_XLIM comes from results.style in Section 0.

set_lab <- function(s) { l <- unname(SET_LABELS[as.character(s)]); ifelse(is.na(l), as.character(s), l) }

# Display-only short group names: GROUP_SHORT_FIG lives in utils/results_core.R
# (shared with notebooks/hindcast_skill.qmd).
group_short_fig <- function(x) unname(GROUP_SHORT_FIG[match(x, GROUP_SHORT)])

# ---------------------------------------------------------------------------
# Section 2b: cross-figure conventions (Figures 1-4)
# ---------------------------------------------------------------------------
# Shared conventions:
#   width        — one double-column width for every figure (WIDTH_MM, the
#                  save_figure() default).
#   typeface     — FONT everywhere, maps included (taxdat::map_theme() leaves the
#                  family to the device default): ps_map_font().
#   ENSO phase   — Figure 1's red / grey / blue at alpha 0.075, labels
#                  El Nino / Neutral / La Nina (PS_ENSO_*).
#   legends      — title BASE_SIZE - 1, labels BASE_SIZE - 1.5 (ps_legend_theme()).
#   counterfactual wording — "neutral-ENSO counterfactual" (ps_cf_label()).
#   uncertainty  — median + 50% + 95% bands (phase-09 q025 / q975 columns).
#   map greys    — two states with a legend key: "Not modelled" (outside the
#                  fitted ENSO groups, or no outbreak in the survival panel) and
#                  "Modelled, no ADM1 estimate" (outbreak countries / ADM1s
#                  without an estimate, incl. low-signal ADM1s in Figure 4b).

PS_ENSO_PHASE_COL   <- c(elnino = "red", neutral = "gray", nina = "blue")
PS_ENSO_PHASE_LAB   <- c(elnino = "El Niño", neutral = "Neutral", nina = "La Niña")
PS_ENSO_PHASE_ALPHA <- 0.075

PS_GREY <- c(not_modelled = "grey85", no_estimate = "grey93")
PS_GREY_LAB <- c(not_modelled = "Not modelled", no_estimate = "Modelled, no ADM1 estimate")
# Same labels wrapped, for narrow legends (Figure 1b inset).
PS_GREY_LAB_WRAP <- c(not_modelled = "Not modelled", no_estimate = "Modelled,\nno ADM1 estimate")

# K for thousands everywhere.
# ONI strip cell height (mm), both Figure 3 strips.
PS_ONI_STRIP_MM <- 14
PS_LABEL_K     <- scales::label_number(scale_cut = scales::cut_short_scale())
# One ONI y range for every ONI strip (Figures 1a, 3a, 3b): 2015/16 peak 2.6,
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

# ---------------------------------------------------------------------------
# Section 2c: shared data access and map helpers
# ---------------------------------------------------------------------------
#' One phase-09 artefact of `kind` (a RESULTS_SCHEMA key) via results_require():
#' refuses a missing, stale or differently-decided artefact. `...` are the
#' artefact's keys (track = , sets = ).
read_results <- function(kind, ...)
  results_require(makeResultsRunFile(opt, paste0("results_", kind), ...,
                                     schema = RESULTS_SCHEMA[[kind]], config_path = CONFIG_PATH),
                  CONFIG_PATH, RES)

#' The ENSO exposure-response artefact of every track in `tracks`, named by track.
read_enso_relation <- function(tracks)
  setNames(lapply(tracks, function(trk) read_results("enso_relation", track = trk)), tracks)

#' Modelled countries for the map greys: >= 1 onset in the survival panel
#' (phase-09 onset_rate; Figure 1b's definition). `rate`: that artefact, if
#' already read.
modelled_countries <- function(rate = read_results("onset_rate"))
  unique(rate$country[rate$onsets > 0])

#' The Africa basemap polygons Figures 2 and 4 draw on: ADM1 units, and the
#' country outlines without the small islands (drop_small_islands()).
load_africa_maps <- function(what) {
  geo <- load_africa_sf()
  if (is.null(geo$adm1_sf)) stop("[10] ", what, ": ADM1 shapefile failed to load")
  list(adm1_sf = geo$adm1_sf, afro_sf = drop_small_islands(geo$afro_sf))
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
  oni$phases %>%
    dplyr::filter(TR >= x_limits[1], TL <= x_limits[2]) %>%
    dplyr::mutate(TL = pmax(TL, x_limits[1]), TR = pmin(TR, x_limits[2]),
                  phase = factor(phase, levels = names(PS_ENSO_PHASE_COL)))
}

#' A blank left facet strip (outside the axis, transparent text), so patchwork
#' aligns this plot's strip column with a neighbour's left strips instead of
#' inserting an empty column between its axis and panel.
blank_left_strip <- function(p, size = BASE_SIZE) {
  p + facet_grid(rows = vars("x"), switch = "y") +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 90, face = "bold", size = size, colour = "transparent"))
}

# ===========================================================================
# Section 9a: Figure 1 builder — ENSO exposure-response overview (stacked)
# ===========================================================================

#' Figure 1b: Africa map of onsets per surveillance-covered month per gid
#' (phase-09 `onset_rate` artefact). Fill and study-country borders in the ENSO
#' group colour, alpha = rate; country greys as every map (PS_GREY: "Not
#' modelled" / "Modelled, no ADM1 estimate", modelled = >= 1 onset in the
#' survival panel), keyed in the caption rather than in the inset. gid tiles
#' from the GADM ADM_1 + ADM_2 layers, GIDs flattened as in
#' descriptive_occurrence_size.Rmd section 3.2. `groups` are the ENSO groups
#' Figure 1 shows (their countries get a coloured border).
build_fig1_inset_map <- function(groups) {
  rate <- read_results("onset_rate")

  iso_sql  <- paste0("'", paste(unique(rate$country), collapse = "','"), "'")
  read_lvl <- function(layer, col) {
    sf::st_read(FIG1_GADM_PATH, quiet = TRUE,
                query = sprintf('SELECT * FROM "%s" WHERE GID_0 IN (%s)', layer, iso_sql)) %>%
      dplyr::transmute(gid = flattenGid(as.character(.data[[col]])))
  }
  gid_sf <- dplyr::bind_rows(read_lvl("ADM_1", "GID_1"), read_lvl("ADM_2", "GID_2")) %>%
    dplyr::filter(gid %in% rate$gid) %>%
    sf::st_make_valid() %>%
    rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE) %>%
    dplyr::inner_join(rate, by = "gid")
  missing <- setdiff(rate$gid, gid_sf$gid)
  if (length(missing))
    message(sprintf("[10] Figure 1b: %d gid(s) without a GADM polygon: %s",
                    length(missing), paste(utils::head(missing, 10), collapse = ",")))

  afro_sf <- ps_country_status(drop_small_islands(load_africa_sf(gadm_path = FIG1_GADM_PATH)$afro_sf),
                               modelled_countries(rate))
  grp_sf  <- afro_sf %>%
    dplyr::inner_join(dplyr::distinct(rate, country, enso_group), by = c("ISO3_CODE" = "country")) %>%
    dplyr::filter(enso_group %in% groups)

  ggplot() +
    ps_country_grey_layer(afro_sf, legend_order = 2, labels = PS_GREY_LAB_WRAP, show_legend = FALSE) +
    geom_sf(data = gid_sf, aes(fill = enso_group, alpha = rate), colour = NA) +
    geom_sf(data = afro_sf, fill = NA, colour = "grey55", linewidth = 0.1) +
    geom_sf(data = grp_sf, aes(colour = enso_group), fill = NA, linewidth = 0.12) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_alpha_continuous(range = c(0.1, 1), trans = "sqrt",
                           name = "Onsets per\ncovered month",
                           labels = scales::label_percent(accuracy = 1),
                           guide = guide_legend(order = 1, nrow = 2)) +
    coord_sf(xlim = FIG1_INSET_BBOX[c("xmin", "xmax")], ylim = FIG1_INSET_BBOX[c("ymin", "ymax")],
             expand = FALSE, datum = NA) +
    labs(tag = "b") +
    theme_void(base_size = 5, base_family = FONT) +
    # Compact legend (title above a two-row rate key) in the South Atlantic
    # corner of the inset.
    theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02),
          legend.justification = c(0, 0), legend.key.size = unit(1.3, "mm"),
          legend.spacing.y = unit(0.4, "mm"), legend.text.position = "right",
          legend.box.just = "left",
          legend.title.position = "top",
          legend.title = element_text(size = 4.5, family = FONT, vjust = 0.5),
          legend.text = element_text(size = 4, family = FONT),
          plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2),
          plot.tag = element_text(size = BASE_SIZE + 2, face = "bold", family = FONT),
          plot.tag.position = c(0.06, 0.95),
          plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm"))
}

#' Place the Figure 1b inset in panel a's top-right corner, FIG1_INSET_MARGIN
#' in from the border, its bottom halfway into the gap below the k-th row from
#' the top: the largest k <= FIG1_INSET_MAX_ROWS whose box (width set by the
#' map's fixed aspect, so the map fills it) covers no case bubble. Fractions
#' come from panel a's own built x/y ranges. `country_levels` is panel a's
#' factor order (bottom to top), `last_bubble` the per-country last case-bubble
#' date as numeric (NA = none).
add_fig1_inset <- function(p_cov, inset, country_levels, last_bubble) {
  pp <- ggplot_build(p_cov)$layout$panel_params[[1]]
  xr <- pp$x.range; yr <- pp$y.range
  n  <- length(country_levels)
  right <- 1 - FIG1_INSET_MARGIN
  top   <- 1 - FIG1_INSET_MARGIN * FIG1_PANEL_A_ASPECT
  box_for <- function(k) {
    y_bottom <- n - k + 0.5
    bottom   <- (y_bottom - yr[1]) / diff(yr)
    left     <- right - (top - bottom) / (FIG1_INSET_ASPECT * FIG1_PANEL_A_ASPECT)
    covered  <- country_levels[seq_len(n) > y_bottom]
    lb <- last_bubble[as.character(covered)]
    list(left = left, bottom = bottom, covered = covered,
         ok = all(is.na(lb) | lb < xr[1] + left * diff(xr)))
  }
  ks <- rev(seq_len(min(FIG1_INSET_MAX_ROWS, n)))
  k  <- ks[vapply(ks, function(k) box_for(k)$ok, logical(1))][1]
  if (is.na(k)) {
    warning("[10] Figure 1b: no inset size avoids panel a's case bubbles; inset omitted")
    return(p_cov)
  }
  bx <- box_for(k)
  message(sprintf("[10] Figure 1b inset: %d row(s) (%s), panel fractions left=%.3f bottom=%.3f",
                  k, paste(rev(as.character(bx$covered)), collapse = ","), bx$left, bx$bottom))
  p_cov + patchwork::inset_element(inset, left = bx$left, bottom = bx$bottom,
                                   right = right, top = top, align_to = "panel", clip = TRUE)
}

#' Figure 1a's ONI strip: Figure 3a's observed-ONI strip (fig3_oni_strip())
#' over panel a's window, carrying the "a" tag, with tighter margins and panel
#' a's tick size.
fig1_oni_strip <- function(oni, x_limits, phases) {
  fig3_oni_strip(oni, NULL, x_limits, mode = "observed", phases = phases) +
    labs(tag = "a") +
    theme(axis.title.y = element_text(margin = margin(r = 0.5, unit = "mm")),
          plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "mm")) +
    theme(axis.text.y = element_text(size = FIG1_A_TICK_PT))
}

# Figure 1 right column (c / d / e). Columns = ENSO
# teleconnection group (c carries the group header), rows of c and e = track.
#   FIG1_A_TICK_PT       panel a's tick labels (countries, years, ONI).
#   FIG1_TAG_POS / FIG1_TAG_GUTTER_MM   c / d / e tags at the far left of each
#     sub-plot, in a gutter (shared by patchwork's alignment) left of c's labels.
#   FIG1_LAG_BREAKS / FIG1_LEAD_RULE    c / d lag axis; dashed lead | lag split.
#   FIG1_ELPD_*          c's fill: dELPD vs the best lag, capped at FIG1_ELPD_CAP
#     (lower values squished), grey floor blending into reversed magma (cut at
#     FIG1_ELPD_MAGMA_END) towards 0 = best; direction arrows on the bar's end
#     labels; the bar FIG1_ELPD_BAR_MM long (the arrow labels need it).
#   FIG1_HEADER*         c's nested header strip over the group names, with a
#     rule; FIG1_STRIP_GAP_MM above the strip text.
#   FIG1_RAIN_CORR_CSV   d's input: the group table (Pearson r
#     of the group-mean 3-month rainfall deviation, baseline 1990-2019, with
#     ONI at lag -11..11, 1981-2024; `*` p < 0.05 and >= half the group's
#     countries significant with the same sign, `.` p < 0.05 for the group
#     mean only). See the data contract in the header.
#   FIG1_ERF_*           e's y axis (Figure 2's fraction labels, up to 10).
#   FIG1_RIGHT_HEIGHTS   c, d, the colourbar row (mm), e.
#   FIG1_E_TOP_MM / FIG1_E_BOTTOM_MM   e's margins; the bottom one aligns e's
#     panel bottom with a's (measured on the render; re-measure when the
#     layout changes).
FIG1_A_TICK_PT      <- BASE_SIZE - 2
FIG1_TAG_POS        <- c(0, 1)
FIG1_TAG_GUTTER_MM  <- 2.2
FIG1_LAG_BREAKS     <- c(-6, 0, 6)
FIG1_LEAD_RULE      <- -0.5
FIG1_ELPD_CAP       <- -25
FIG1_ELPD_GREY      <- "grey90"
FIG1_ELPD_MAGMA_END <- 0.85
FIG1_ELPD_TITLE     <- "ΔELPD vs.\nbest lag"
FIG1_ELPD_ARROWS    <- c(worse = "← worse fit", better = "better fit →")
FIG1_ELPD_BAR_MM    <- 22
FIG1_BAR_MM         <- 16
FIG1_HEADER         <- "ENSO teleconnection group"
FIG1_HEADER_PT      <- BASE_SIZE
FIG1_STRIP_GAP_MM   <- 1.2
FIG1_RAIN_CORR_CSV  <- "generated_data/scratchpad/s68_precip_dev_enso_phase_kl_groups/group_stats.csv"
FIG1_CORR_COL       <- c(low = "#8c510a", mid = "white", high = "#01665e")
FIG1_CORR_LIM       <- c(-0.7, 0.7)
FIG1_CORR_TITLE     <- "Rainfall anomaly–ONI\ncorrelation (r)"
FIG1_MARK_PT        <- BASE_SIZE - 3.5
FIG1_XTITLE_GAP_MM  <- 0.4
# Figure 1d: main rainy season(s) per ENSO group, outlined over the rainfall
# months. Derived outside the pipeline
# (utils/rain_enso.R::rainSeasonsByGroup(): a month is in a group's season
# when >= 50% of the group's population lives in ADM1s where that month is above
# the ADM1's average month, 1990-2019 climatology); 10 reads the population-
# weighted rows.
FIG1_RAIN_SEASON_CSV <- "generated_data/scratchpad/s80_fig1_rain_seasons/rain_seasons_boxes.csv"
FIG1_RAIN_SEASON_WEIGHTING <- "pop"
FIG1_RAIN_SEASON_COL <- "black"
FIG1_RAIN_SEASON_LWD <- 0.35
FIG1_ERF_YLIM       <- c(1/3, 10)
FIG1_ERF_BREAKS     <- c(1/3, 1/2, 1, 2, 5, 10)
FIG1_ERF_LABELS     <- c("1/3", "1/2", "1", "2", "5", "10")
FIG1_ONI_BREAKS     <- c(-2, 0, 2)
FIG1_E_TOP_MM       <- 0.3
FIG1_E_BOTTOM_MM    <- 2.24
FIG1_RIGHT_HEIGHTS  <- grid::unit(c(18.7, 28.3, 6.6, 74.6), c("null", "null", "mm", "null"))

# Axis / strip text of Figure 1c, 1d and 1e.
FIG1_CD_THEME <- theme(axis.text = element_text(size = BASE_SIZE - 1.5, colour = "grey20"),
                       axis.title = element_text(size = BASE_SIZE - 1),
                       strip.text.x = element_text(size = BASE_SIZE - 1.5, margin = margin(b = 0.6, unit = "mm")))

# Strips, legends and tags shared by Figure 1c, 1d and 1e.
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

#' Group short labels for Figure 1's group strips.
fig1_group_labeller <- function() ggplot2::as_labeller(GROUP_SHORT_FIG)

#' Horizontal colourbar for Figure 1's c / d legend row, title left of the bar.
#' `text`: an optional legend.text element (left unset when NULL -- an explicit
#' NULL in theme() would blank the element and drop the labels to the default size).
fig1_bar <- function(bar_mm = FIG1_BAR_MM, text = NULL) {
  th <- theme(legend.title.position = "left",
              legend.key.width = unit(bar_mm, "mm"),
              legend.key.height = unit(1.6, "mm"),
              legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.9))
  if (!is.null(text)) th <- th + theme(legend.text = text)
  guide_colourbar(direction = "horizontal", theme = th)
}

#' Figure 1c's fill colours, light to dark: the grey floor, then magma from
#' FIG1_ELPD_MAGMA_END down to black.
fig1_elpd_colours <- function(n = 9)
  c(FIG1_ELPD_GREY, rev(viridisLite::magma(n, end = FIG1_ELPD_MAGMA_END)))

#' Figure 1c: ENSO lag selection. dELPD vs the best lag (phase-09 lag_rank of
#' each track), rows = track, columns = ENSO teleconnection group under a nested
#' header; the selected lag outlined in white. `lag_rank`: the tracks' lag_rank
#' tables bound with a `track` column.
fig1_elpd_panel <- function(lag_rank) {
  cap <- FIG1_ELPD_CAP
  elpd_data <- lag_rank %>%
    dplyr::mutate(enso_group = factor(enso_group, levels = GROUP_ORDER),
                  track = track_factor(track, wrap = TRUE),
                  header = FIG1_HEADER)
  best <- dplyr::filter(elpd_data, selected)
  # Explicit width: the sparse `best` outline layer would otherwise resolve a
  # different tile width from the full lag grid.
  tile_width <- ggplot2::resolution(elpd_data$single_lag)
  strip_margin <- margin(t = FIG1_STRIP_GAP_MM, b = 0.6, unit = "mm")
  ggplot(elpd_data, aes(x = single_lag, y = track, fill = elpd_diff)) +
    geom_tile(width = tile_width) +
    geom_tile(data = best, fill = NA, colour = "white", linewidth = 0.3, width = tile_width) +
    geom_vline(xintercept = FIG1_LEAD_RULE, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    # Nested header over the group columns, a rule under it; header and group
    # names sized per strip layer.
    ggh4x::facet_nested(. ~ header + enso_group,
                        labeller = labeller(enso_group = fig1_group_labeller(), header = label_value),
                        nest_line = element_line(colour = "black", linewidth = 0.25),
                        strip = ggh4x::strip_nested(
                          by_layer_x = TRUE,
                          text_x = ggh4x::elem_list_text(size = c(FIG1_HEADER_PT, BASE_SIZE - 1),
                                                         margin = list(strip_margin, strip_margin)))) +
    scale_fill_gradientn(
      colours = fig1_elpd_colours(), name = FIG1_ELPD_TITLE,
      limits = c(cap, 0), breaks = c(cap, cap / 2, 0),
      labels = c(paste0("≤ ", cap, "\n", FIG1_ELPD_ARROWS[["worse"]]), cap / 2,
                 paste0("0\n", FIG1_ELPD_ARROWS[["better"]])),
      oob = scales::squish,
      # Per-label justification (cap | mid | 0): the end labels flush with the
      # bar's ends, so the arrows stay within its width. Vectorised hjust is not
      # officially supported by ggplot2 (warns), but guide labels honour it.
      guide = fig1_bar(FIG1_ELPD_BAR_MM,
                       text = suppressWarnings(element_text(size = BASE_SIZE - 2, vjust = 1, lineheight = 0.9,
                                                            hjust = c(0, 0.5, 1))))) +
    scale_x_continuous(breaks = FIG1_LAG_BREAKS, expand = expansion(add = 0.9)) +
    scale_y_discrete(limits = rev(levels(elpd_data$track)), expand = expansion(add = 0.6)) +
    labs(x = NULL, y = NULL, tag = "c") +
    theme_ms() + ps_legend_theme() + FIG1_CD_THEME + FIG1_PANEL_THEME +
    theme(axis.text.y = element_text(face = "plain", colour = "black", size = BASE_SIZE - 2),
          strip.text.x = element_text(size = BASE_SIZE - 1, face = "plain", margin = strip_margin),
          # No x axis (d's serves both); the tag gutter on the left.
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length.x = unit(0, "mm"),
          plot.margin = margin(1.5, 1.5, 0, FIG1_TAG_GUTTER_MM, unit = "mm"),
          legend.position = "bottom")
}

#' Figure 1d: rainfall anomaly-ONI correlation per rainfall month (rows) and
#' lag, per ENSO teleconnection group (columns), from FIG1_RAIN_CORR_CSV,
#' cropped to `lag_range` (c's lags). `groups`: the groups Figure 1 shows.
fig1_corr_panel <- function(groups, lag_range) {
  # Written outside the pipeline; 10 only names its producer in the message.
  if (!file.exists(FIG1_RAIN_CORR_CSV))
    stop("[10] Figure 1d: ", FIG1_RAIN_CORR_CSV, " not found -- run the s68 scratchpad ",
         "(s68_precip_dev_enso_phase_kl_groups.R)", call. = FALSE)
  st <- readr::read_csv(FIG1_RAIN_CORR_CSV, show_col_types = FALSE) %>%
    dplyr::filter(series %in% groups, lag >= lag_range[1], lag <= lag_range[2]) %>%
    dplyr::mutate(enso_group = factor(series, levels = GROUP_ORDER),
                  centre = factor(month.abb[centre], levels = rev(month.abb)),
                  mark = dplyr::case_when(p_r < .05 & agree_r >= .5 ~ "*", p_r < .05 ~ "·", TRUE ~ ""))
  miss <- setdiff(groups, unique(st$series))
  if (length(miss))
    stop("[10] Figure 1d: ", FIG1_RAIN_CORR_CSV, " has no series for group(s) ",
         paste(miss, collapse = ","), call. = FALSE)
  # Rainy seasons: written outside the pipeline, as above. One row per
  # season (first .. last calendar month, possibly wrapping the year); one
  # rectangle per run of consecutive months on the y positions of
  # rev(month.abb), so a Dec -> Jan season gives two.
  if (!file.exists(FIG1_RAIN_SEASON_CSV))
    stop("[10] Figure 1d: ", FIG1_RAIN_SEASON_CSV, " not found -- run the s80 scratchpad ",
         "(s80_fig1_rain_seasons.R)", call. = FALSE)
  seasons <- readr::read_csv(FIG1_RAIN_SEASON_CSV, show_col_types = FALSE) %>%
    dplyr::filter(weighting == FIG1_RAIN_SEASON_WEIGHTING, group %in% groups)
  miss <- setdiff(groups, unique(seasons$group))
  if (length(miss))
    stop("[10] Figure 1d: ", FIG1_RAIN_SEASON_CSV, " has no ", FIG1_RAIN_SEASON_WEIGHTING,
         "-weighted season for group(s) ", paste(miss, collapse = ","),
         " (re-run s80 with Figure 1's grouping)", call. = FALSE)
  seasons <- seasons %>%
    dplyr::rowwise() %>%
    dplyr::reframe(series = group, season = season,
                   month = ((first - 1) + 0:((last - first) %% 12)) %% 12 + 1) %>%
    dplyr::group_by(series, season) %>%
    dplyr::mutate(run = cumsum(c(1, diff(month) != 1))) %>%
    dplyr::mutate(y = 13 - month) %>%
    dplyr::group_by(series, season, run) %>%
    dplyr::summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop") %>%
    dplyr::mutate(enso_group = factor(series, levels = GROUP_ORDER),
                  xmin = lag_range[1] - 0.5, xmax = lag_range[2] + 0.5)
  ggplot(st, aes(x = lag, y = centre, fill = r)) +
    geom_tile() +
    geom_rect(data = seasons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, colour = FIG1_RAIN_SEASON_COL,
              linewidth = FIG1_RAIN_SEASON_LWD) +
    geom_text(aes(label = mark), size = FIG1_MARK_PT / ggplot2::.pt, family = FONT, vjust = 0.75) +
    geom_vline(xintercept = FIG1_LEAD_RULE, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    facet_grid(. ~ enso_group, labeller = fig1_group_labeller()) +
    scale_fill_gradient2(low = FIG1_CORR_COL[["low"]], mid = FIG1_CORR_COL[["mid"]],
                         high = FIG1_CORR_COL[["high"]], midpoint = 0, limits = FIG1_CORR_LIM,
                         oob = scales::squish, name = FIG1_CORR_TITLE, breaks = c(-0.5, 0, 0.5),
                         guide = fig1_bar()) +
    scale_x_continuous(breaks = FIG1_LAG_BREAKS, expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    labs(x = "ENSO lag (months; < 0 = lead)", y = "Rainfall month", tag = "d") +
    theme_ms() + ps_legend_theme() + FIG1_CD_THEME + FIG1_PANEL_THEME +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = BASE_SIZE - 2.5),
          strip.text.x = element_blank(),
          plot.margin = margin(0, 1.5, 0, 0, unit = "mm"),
          axis.title.x = element_text(size = BASE_SIZE - 1,
                                      margin = margin(t = FIG1_XTITLE_GAP_MM, unit = "mm")),
          legend.position = "bottom")
}

#' Figure 1e: ERF grid, rows = track, columns = ENSO teleconnection group (50%
#' + 95% CrI) on the ENSO-phase shading; arrows in the first column say what
#' an effect above / below 1 means per track (e plots the cessation rate
#' itself: above 1 = ends sooner). `erf`: the tracks' erf tables bound with a
#' `track` column.
fig1_erf_panel <- function(erf) {
  erf_data <- erf %>%
    dplyr::mutate(enso_group = factor(enso_group, levels = GROUP_ORDER),
                  track = track_factor(track, erf = TRUE, wrap = TRUE))
  ann <- data.frame(track = factor(levels(erf_data$track), levels = levels(erf_data$track)),
                    enso_group = factor(levels(erf_data$enso_group)[1], levels = levels(erf_data$enso_group)),
                    x = ERF_XLIM[1] + 0.12,
                    up   = c("more\nonsets", "shorter\noutbreaks", "more\ncases"),
                    down = c("fewer\nonsets", "longer\noutbreaks", "fewer\ncases"))
  arrow_layer <- function(y, yend)
    geom_segment(data = ann, aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE,
                 arrow = grid::arrow(length = unit(0.9, "mm"), type = "closed"), linewidth = 0.25,
                 colour = "grey25")
  label_layer <- function(y, label, vjust)
    geom_text(data = ann, aes(x = x + 0.28, y = y, label = .data[[label]]), inherit.aes = FALSE,
              hjust = 0, vjust = vjust, lineheight = 0.85,
              size = (BASE_SIZE - 2.5) / ggplot2::.pt, family = FONT, colour = "grey25")
  ggplot(erf_data, aes(x = at, y = IRR, ymin = IRR_lo, ymax = IRR_hi,
                       colour = enso_group, fill = enso_group)) +
    annotate("rect", xmin = -3, xmax = -.5, ymin = 0, ymax = 100, fill = "blue", alpha = PS_ENSO_PHASE_ALPHA) +
    annotate("rect", xmin = -.5, xmax = .5, ymin = 0, ymax = 100, fill = "gray", alpha = PS_ENSO_PHASE_ALPHA) +
    annotate("rect", xmin = .5, xmax = 3, ymin = 0, ymax = 100, fill = "red", alpha = PS_ENSO_PHASE_ALPHA) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = IRR_lo50, ymax = IRR_hi50), alpha = 0.35, colour = NA) +
    geom_line(linewidth = 0.6) +
    arrow_layer(1.15, 6) + arrow_layer(0.87, 0.36) +
    label_layer(6, "up", 1) + label_layer(0.36, "down", 0) +
    # Track strips on the left, outside the axis; no group strips (c's serve).
    facet_grid(track ~ enso_group, switch = "y", labeller = labeller(enso_group = fig1_group_labeller())) +
    # One "ENSO group" key (colour and fill merged) under e.
    scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_SHORT_FIG, name = "ENSO group") +
    scale_fill_manual(values = GROUP_COLOURS, labels = GROUP_SHORT_FIG, name = "ENSO group") +
    scale_y_log10(breaks = FIG1_ERF_BREAKS, labels = FIG1_ERF_LABELS) +
    scale_x_continuous(breaks = FIG1_ONI_BREAKS) +
    coord_cartesian(xlim = ERF_XLIM, ylim = FIG1_ERF_YLIM) +
    labs(x = "Oceanic Niño Index (ONI)", y = "Effect size on rate of...", tag = "e") +
    theme_ms() + ps_legend_theme() + FIG1_CD_THEME + FIG1_PANEL_THEME +
    theme(strip.text.x = element_blank(),
          # y title hugging the track strips; bottom margin aligns e's panel
          # bottom with a's.
          axis.title.y = element_text(size = BASE_SIZE - 1, margin = margin(r = 0.5, unit = "mm")),
          plot.margin = margin(FIG1_E_TOP_MM, 1.5, FIG1_E_BOTTOM_MM, 0, unit = "mm"),
          axis.title.x = element_text(margin = margin(t = FIG1_XTITLE_GAP_MM, unit = "mm")))
}

#' Figure 1c and 1d's colourbars side by side, cut from the two plots' own
#' gtables (nothing is collected: e keeps its legend in place).
fig1_legend_row <- function(plots, gap_mm = 4) {
  bars <- lapply(plots, function(p) {
    g <- ggplot2::ggplotGrob(p)
    g$grobs[[which(g$layout$name == "guide-box-bottom")]]
  })
  gtable::gtable_matrix("fig1_cd_legends", matrix(list(bars[[1]], grid::nullGrob(), bars[[2]]), nrow = 1),
                        widths = grid::unit.c(sum(bars[[1]]$widths), grid::unit(gap_mm, "mm"),
                                              sum(bars[[2]]$widths)),
                        heights = grid::unit(1, "null"))
}

#' Figure 1a's data-coverage panel: per-country outbreak-period timeline on the
#' ENSO-phase background, monthly suspected cases as bubbles, the inset map (b)
#' in its top-right corner. Reads the raw survival / coverage / monthly data
#' directly -- a narrow exception to the data contract above, scoped to this
#' one panel. `groups`: the ENSO groups Figure 1 shows (inset borders).
fig1_coverage_panel <- function(phases, groups) {
  enso_groups <- read_csv(GROUPING_CSV, show_col_types = FALSE)
  simEpisodeCheck(preps = list(surv = makeSurvDatFile(opt), mon_data = makeMonthlyDatFile(opt)),
                  what = "10_manuscript_figures data-coverage panel",
                  episode_key = episodeKeyOf(opt))
  surv <- readRDS(makeSurvDatFile(opt))
  cov <- readRDS(makeAggTsFile()) %>%
    filter(country %in% unique(surv$country)) %>%
    inner_join(enso_groups, by = c("country" = "iso3"))

  # Monthly case bubbles (non-ADM0 rows).
  mon_pts <- readRDS(makeMonthlyDatFile(opt)) %>%
    filter(country %in% unique(cov$country), admin_level != "ADM0")
  last_bubble <- mon_pts %>%
    group_by(country) %>%
    summarise(last_bubble = max(as.Date(month)), .groups = "drop")

  # Country axis order: grouped by ENSO group; within a group by date of the
  # last case bubble, latest at the bottom (none -> top), so each block's top
  # rows are the countries whose data stop earliest -- room for the inset.
  cov_ctry <- cov %>%
    distinct(country, group) %>%
    left_join(last_bubble, by = "country") %>%
    arrange(group, desc(dplyr::coalesce(as.numeric(last_bubble), -Inf)), country)
  ctry_levels <- unique(cov_ctry$country)

  ob_pts <- mon_pts %>%
    group_by(country, month) %>%
    summarise(tot_cases = sum(cases_month), .groups = "drop") %>%
    inner_join(enso_groups, by = c("country" = "iso3")) %>%
    mutate(country = factor(country, levels = ctry_levels))
  cov <- mutate(cov, country = factor(country, levels = ctry_levels))

  p_cov <- ggplot(cov) +
    geom_rect(data = phases, aes(xmin = TL, xmax = TR, ymin = -Inf, ymax = Inf, fill = phase),
              alpha = PS_ENSO_PHASE_ALPHA, inherit.aes = FALSE) +
    geom_segment(aes(x = TL - 5, xend = TR + 5, y = country, yend = country),
                 linewidth = 3.2, colour = "gray40") +
    geom_segment(aes(x = TL, xend = TR, y = country, yend = country),
                 linewidth = 3, colour = "white") +
    geom_point(data = ob_pts, aes(x = month, y = country, size = tot_cases,
                                  colour = group), alpha = 0.45) +
    scale_size_continuous(trans = "sqrt", range = c(0.1, 3.5),
                          breaks = c(5000, 10000, 20000), labels = PS_LABEL_K) +
    # Window set on the COORD, not the scale: scale limits drop (not clip)
    # every coverage segment starting before the window, which cut bars short
    # of their bubbles.
    scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = FIG1_A_XLIM) +
    labs(x = NULL, y = NULL, size = "Cases/month") +
    scale_color_manual(values = GROUP_COLOURS) +
    scale_fill_manual(values = PS_ENSO_PHASE_COL, labels = PS_ENSO_PHASE_LAB, drop = FALSE) +
    guides(color = "none", size = guide_legend(order = 1),
           fill = guide_legend("ENSO", order = 2,
                               override.aes = list(alpha = 3 * PS_ENSO_PHASE_ALPHA))) +
    theme_ms() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.2, "mm"),
          legend.box.spacing = unit(1, "mm")) +
    ps_legend_theme() +
    # Both legends on one line, title left of its keys, at e's legend text
    # sizes; compact keys and spacing so the line fits panel a's width.
    theme(legend.title.position = "left", legend.box = "horizontal",
          legend.box.just = "left", legend.justification = "left",
          legend.spacing.x = unit(1.2, "mm"), legend.spacing.y = unit(0.3, "mm"),
          legend.key.height = unit(2.6, "mm"), legend.key.spacing.x = unit(0.2, "mm"),
          legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.5),
          legend.text = element_text(size = BASE_SIZE - 2),
          legend.box.margin = margin(0, 0, 0, -4, unit = "mm"),
          plot.margin = margin(0.3, 0.5, 0.5, 0.5, unit = "mm")) +
    theme(legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.5),
          legend.text = element_text(size = BASE_SIZE - 2, margin = margin(l = 0.4, unit = "mm")),
          legend.spacing.x = unit(2, "mm"),
          legend.key.spacing.x = unit(0.1, "mm"),
          legend.key.width = unit(1.8, "mm"),
          legend.title.position = "left",
          legend.box.margin = margin(0, 0, 0, -5, unit = "mm")) +
    theme(axis.text.y = element_text(size = FIG1_A_TICK_PT)) +
    theme(axis.text.x = element_text(size = FIG1_A_TICK_PT))

  add_fig1_inset(p_cov, build_fig1_inset_map(groups),
                 country_levels = ctry_levels,
                 last_bubble = stats::setNames(as.numeric(cov_ctry$last_bubble), cov_ctry$country))
}

#' Figure 1: left column = observed-ONI strip over the data-coverage timeline
#' (a) with the inset map of onsets per covered month per gid (b); right column
#' = ENSO lag selection (c) and the rainfall anomaly-ONI correlation (d) with
#' their colourbars in one row, then the ERF grid (e). ENSO phases (a's
#' background) as Figure 3: NOAA episodes from the phase-09 ONI artefact.
build_fig1 <- function(tracks = STACKED_TRACKS) {
  message(sprintf("[10] Figure 1: tracks=%s", paste(tracks, collapse = ",")))
  rel <- read_enso_relation(tracks)
  bind_track <- function(el) dplyr::bind_rows(lapply(tracks, function(trk) dplyr::mutate(rel[[trk]][[el]], track = trk)))
  erf <- bind_track("erf")
  lag_rank <- bind_track("lag_rank")
  groups <- intersect(GROUP_ORDER, unique(as.character(erf$enso_group)))

  p_c <- fig1_elpd_panel(lag_rank)
  p_d <- fig1_corr_panel(groups, range(lag_rank$single_lag, na.rm = TRUE))
  # One flat patchwork: a free() plot inside a nested patchwork does not render
  # (patchwork 1.3.2). d's labels are freed so "Rainfall month" hugs the month
  # names instead of aligning with c's wider track labels; e's likewise.
  right <- patchwork::wrap_plots(
    p_c + theme(legend.position = "none"),
    patchwork::free(p_d + theme(legend.position = "none"), type = "label", side = "l"),
    patchwork::wrap_elements(full = fig1_legend_row(list(p_c, p_d))),
    patchwork::free(fig1_erf_panel(erf), type = "label", side = "l"),
    ncol = 1, heights = FIG1_RIGHT_HEIGHTS) +
    # No margin on the side facing a.
    patchwork::plot_annotation(theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 0, unit = "pt")))
  right_block <- patchwork::wrap_elements(full = right) + theme(plot.margin = margin(0, 0, 0, 0))

  oni <- read_results("oni", sets = SETS)
  phases <- enso_phase_rects(oni, FIG1_A_XLIM)
  # ONI strip on top of panel a (tag "a" on the strip); "ONI" next to its own
  # tick labels (free label: not pushed out to the width of a's country
  # labels). Left column laid out on its own (wrap_elements): nested, patchwork
  # would align the strip's row with c's and open a gap above a. On the side
  # facing the right column only FIG1_A_RIGHT_MM (5.5 pt = 1.94 mm elsewhere).
  left <- patchwork::wrap_plots(
    patchwork::free(fig1_oni_strip(oni, FIG1_A_XLIM, phases), type = "label", side = "l"),
    fig1_coverage_panel(phases, groups),
    ncol = 1, heights = grid::unit(c(FIG1_ONI_STRIP_MM, 1), c("mm", "null"))) +
    patchwork::plot_annotation(theme = theme(plot.margin = margin(1.94, FIG1_A_RIGHT_MM, 1.94, 1.94, unit = "mm")))
  left_block <- patchwork::wrap_elements(full = left) + theme(plot.margin = margin(0, 0, 0, 0))

  fig <- (left_block | right_block) + patchwork::plot_layout(widths = FIG1_WIDTHS)
  save_figure(fig, "fig1_stacked", height_mm = FIG1_HEIGHT_MM)
}

# ===========================================================================
# Section 9b: Figure 2 builder — ADM1 IRR maps (stacked)
# ===========================================================================
#' Per-column rule under Figure 2's ONI header row: one segment per panel
#' column at the bottom of the last top-strip row, trimmed by FIG2_RULE_TRIM at
#' each end and as wide as the phase rules above it (facet_nested()'s
#' nest_line). Returns a gtable.
add_fig2_header_rules <- function(p) {
  g <- ggplotGrob(p)
  lay <- g$layout
  st  <- lay[grepl("^strip-t", lay$name), ]
  panel_cols <- sort(unique(lay$l[grepl("^panel", lay$name)]))
  rule <- grid::segmentsGrob(x0 = FIG2_RULE_TRIM, x1 = grid::unit(1, "npc") - FIG2_RULE_TRIM,
                             y0 = grid::unit(0, "npc"), y1 = grid::unit(0, "npc"),
                             gp = grid::gpar(col = "black", lwd = FIG2_RULE_LWD * ggplot2::.pt,
                                             lineend = "butt"))
  for (k in seq_along(panel_cols))
    g <- gtable::gtable_add_grob(g, rule, t = max(st$b), b = max(st$b), l = panel_cols[k],
                                 r = panel_cols[k], z = Inf, clip = "off",
                                 name = paste0("oni-rule-", k))
  g
}

#' Narrow each phase header's fill (the strips spanning several panel columns)
#' to the rule line under it: FIG2_RULE_TRIM in from both ends. Works on the
#' gtable from add_fig2_header_rules().
trim_fig2_phase_fill <- function(g) {
  for (i in which(grepl("^strip-t", g$layout$name) & g$layout$l != g$layout$r)) {
    s <- g$grobs[[i]]
    for (k in seq_along(s$grobs)) {
      ch <- s$grobs[[k]]$children
      if (is.null(ch)) next
      for (j in seq_along(ch)) if (inherits(ch[[j]], "rect"))
        ch[[j]]$width <- grid::unit(1, "npc") - 2 * FIG2_RULE_TRIM
      s$grobs[[k]]$children <- ch
    }
    g$grobs[[i]] <- s
  }
  g
}

#' Append a column right of the maps holding, per track row, an unlabelled
#' vertical copy of the map colourbar (same gradient2, Lab-interpolated) with
#' an up / down arrow from its midpoint (IRR = 1) and the direction labels
#' (FIG2_ROW_KEY_LABELS) at the arrow tips, as in Figure 1d.
add_fig2_row_keys <- function(g, track_labels, width_mm = FIG2_ROW_KEY_MM,
                              bar_y = c(0.22, 0.78)) {
  rows <- sort(unique(g$layout$t[grepl("^panel", g$layout$name)]))
  stopifnot(length(rows) == length(track_labels))
  pal  <- scales::div_gradient_pal("#2166ac", "white", "#d73027", "Lab")
  ramp <- rev(pal(seq(0, 1, length.out = 101)))   # top (red) -> bottom (blue)
  txt  <- grid::gpar(fontsize = BASE_SIZE - 1.5, fontfamily = FONT, col = "grey25", lineheight = 0.85)
  arw  <- grid::arrow(length = grid::unit(0.9, "mm"), type = "closed")
  agp  <- grid::gpar(col = "grey25", fill = "grey25", lwd = 0.25 * ggplot2::.pt)
  mid  <- mean(bar_y); gap <- 0.015
  g <- gtable::gtable_add_cols(g, grid::unit(width_mm, "mm"), pos = -1)
  for (k in seq_along(rows)) {
    lab <- FIG2_ROW_KEY_LABELS[[track_labels[k]]]
    if (is.null(lab)) stop("[10] Figure 2: no row-key labels for track '", track_labels[k], "'")
    key <- grid::gTree(children = grid::gList(
      grid::rasterGrob(matrix(ramp, ncol = 1), x = grid::unit(1.5, "mm"), y = grid::unit(mid, "npc"),
                       width = grid::unit(1.2, "mm"), height = grid::unit(diff(bar_y), "npc"),
                       hjust = 0, interpolate = TRUE),
      grid::segmentsGrob(x0 = grid::unit(4.4, "mm"), x1 = grid::unit(4.4, "mm"),
                         y0 = grid::unit(mid + gap, "npc"), y1 = grid::unit(bar_y[2], "npc"),
                         arrow = arw, gp = agp),
      grid::segmentsGrob(x0 = grid::unit(4.4, "mm"), x1 = grid::unit(4.4, "mm"),
                         y0 = grid::unit(mid - gap, "npc"), y1 = grid::unit(bar_y[1], "npc"),
                         arrow = arw, gp = agp),
      grid::textGrob(lab[["up"]], x = grid::unit(5.4, "mm"), y = grid::unit(bar_y[2], "npc"),
                     hjust = 0, vjust = 1, gp = txt),
      grid::textGrob(lab[["down"]], x = grid::unit(5.4, "mm"), y = grid::unit(bar_y[1], "npc"),
                     hjust = 0, vjust = 0, gp = txt)))
    g <- gtable::gtable_add_grob(g, key, t = rows[k], l = ncol(g), clip = "off",
                                 name = paste0("row-key-", k))
  }
  g
}

#' Figure 2: ADM1 IRR maps faceted track (rows, strips on the left) x ONI
#' condition (columns = results.derivation.enso_conds, under an ENSO-phase
#' header filled in Figure 1's phase colours, with booktabs-style rules).
#' Basemap layers as taxdat::output_plot_map() (country fill -> ADM1 fill ->
#' rivers -> lakes -> country borders), built here so the two country greys
#' (PS_GREY; modelled = >= 1 onset, phase-09 `onset_rate`, as Figure 1b) are
#' keyed in the legend. IRR capped to FIG2_MAP_LO/HI; colourbar and grey keys
#' in one row under the maps, a direction key right of each row; scale bar +
#' north arrow on the top-left map. Duration IRR is inverted here so "high =
#' worse" holds on every row.
build_fig2 <- function(tracks = STACKED_TRACKS) {
  message(sprintf("[10] Figure 2 (ADM1 maps, stacked): tracks=%s", paste(tracks, collapse = ",")))
  maps <- load_africa_maps("Figure 2")
  basemap <- load_africa_basemap(maps$afro_sf)   # utils/figures_cb.R (lakes + rivers, st_intersects clip)

  rel <- read_enso_relation(tracks)
  d <- dplyr::bind_rows(lapply(tracks, function(trk) dplyr::mutate(rel[[trk]]$adm1_irr, track = trk))) %>%
    dplyr::mutate(IRR = dplyr::if_else(track == "duration", 1 / IRR, IRR),
                  log_ratio = log10(pmin(pmax(IRR, FIG2_MAP_LO), FIG2_MAP_HI)))
  # enso_label "La Niña (−2.0)" -> phase "La Niña" (header row) + cond "ONI −2.0".
  cond_levels <- unique(d$enso_label)
  d <- d %>%
    dplyr::mutate(phase = factor(sub(" \\(.*$", "", enso_label),
                                 levels = unique(sub(" \\(.*$", "", cond_levels))),
                  cond  = factor(sub("^.*\\((.*)\\)$", "ONI \\1", enso_label),
                                 levels = sub("^.*\\((.*)\\)$", "ONI \\1", cond_levels)),
                  track = track_factor(track))

  map_sf <- maps$adm1_sf %>% dplyr::inner_join(d, by = c("GID_1" = "adm1_gid"))
  sig_sf <- dplyr::filter(map_sf, significant)

  # Scale bar + north arrow on the top-left panel only.
  tl <- data.frame(track = factor(levels(droplevels(d$track))[1], levels = levels(d$track)),
                   phase = factor(levels(d$phase)[1], levels = levels(d$phase)),
                   cond  = factor(levels(d$cond)[1],  levels = levels(d$cond)))

  # Header strips: ggh4x orders them outer layer first -- the phase headers in
  # level order (filled), then one ONI-condition strip per column (unfilled).
  n_ph <- nlevels(d$phase); n_cd <- nlevels(d$cond)
  ph_fill <- unname(FIG2_PHASE_FILL[dplyr::case_when(grepl("Niña", levels(d$phase)) ~ "nina",
                                                     grepl("Niño", levels(d$phase)) ~ "elnino",
                                                     TRUE ~ "neutral")])

  p <- ggplot() +
    ps_country_grey_layer(ps_country_status(maps$afro_sf, modelled_countries()), legend_order = 2) +
    geom_sf(data = map_sf, aes(fill = log_ratio), colour = "white", linewidth = FIG2_ADM1_LWD) +
    map_overlay_layers(basemap, maps$afro_sf) +
    geom_sf(data = sig_sf, fill = NA, colour = "black", linewidth = FIG2_SIG_LWD) +
    # Scale bar 1000 km; arrow centred over it.
    map_scale_layers(data = tl) +
    scale_fill_gradient2(
      midpoint = 0, low = "#2166ac", mid = "white", high = "#d73027",
      limits = log10(c(FIG2_MAP_LO, FIG2_MAP_HI)), oob = scales::squish,
      breaks = log10(c(1/3, 1/2, 1, 2, 3)),
      labels = c("≤ 1/3", "1/2", "1", "2", "≥ 3"),
      name = "Effect size",
      guide = guide_colourbar(direction = "horizontal", title.position = "left", order = 1,
                              theme = theme(legend.key.width = unit(32, "mm"),
                                            legend.key.height = unit(2, "mm")))) +
    ggh4x::facet_nested(track ~ phase + cond, switch = "y",
                        nest_line = element_line(colour = "black", linewidth = FIG2_RULE_LWD),
                        resect = FIG2_RULE_TRIM,
                        strip = ggh4x::strip_nested(by_layer_x = FALSE,
                          # colour/linewidth per strip: elem_list_rect() does not
                          # inherit the theme's blank strip.background border.
                          background_x = ggh4x::elem_list_rect(fill = c(ph_fill, rep(NA, n_cd)),
                                                               colour = rep(NA, n_ph + n_cd),
                                                               linewidth = rep(0, n_ph + n_cd)),
                          text_x = ggh4x::elem_list_text(
                            size = c(rep(BASE_SIZE + 1, n_ph), rep(BASE_SIZE - 0.5, n_cd))))) +
    taxdat::map_theme() +
    ps_map_font() +
    theme(panel.border = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size = BASE_SIZE - 1),
          strip.text = element_text(size = BASE_SIZE + 2),
          strip.text.y.left = element_text(size = BASE_SIZE, angle = 90,
                                           margin = margin(r = 3, unit = "mm")),
          legend.justification = c(0, 0),
          legend.background = element_blank()) +
    ps_legend_theme() +
    # The maps are too small to host the legends: one row below the figure
    # (colourbar + the two greys), top-aligned.
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.box.just = "top", legend.spacing.x = unit(6, "mm"),
          legend.title = element_text(size = BASE_SIZE - 1, vjust = 0.8),
          legend.key.height = unit(2.4, "mm"), legend.key.width = unit(3, "mm"))

  # Passed unevaluated: the gtable edits call ggplotGrob(), whose layout depends
  # on the device's font metrics, so they must run inside save_figure() with
  # the PNG/PDF device open (building the grob first gives a different layout).
  save_figure({ g <- add_fig2_header_rules(p)
                g <- trim_fig2_phase_fill(g)
                add_fig2_row_keys(g, levels(droplevels(d$track))) },
              "fig2_adm1_stacked",
              height_mm = FIG2_ROW_H_MM * length(tracks) + FIG2_LEGEND_MM)
}

# ===========================================================================
# Section 9c: Figure 3 builder — smoothing -> forecast vs counterfactual
# ===========================================================================
# Colours, line types and shapes. The Africa-wide arm is bordeaux, regions use
# GROUP_COLOURS, every counterfactual one neutral dark grey (dashed); the
# free-running simulation is dot-dashed; masked smoothing dark goldenrod (so it
# does not share the counterfactual's grey); the counterfactual's point is a
# square, since the open circle is (a)'s observed symbol.
F3_AFRICA_COL <- "#7b1e3c"
F3_CF_COL     <- "grey35"
F3_MASK_COL   <- "#B8860B"
F3_ONI_COL    <- "grey10"         # CFSv2 forecast ONI (dot-dashed, as every free-running series)
F3_ARM_LTY    <- c(forecast = "dotdash", counterfactual = "22")
F3_ARM_SHAPE  <- c(forecast = 16, counterfactual = 15)
F3_LWD        <- list(a = 0.35, b = 0.4, c = c(0.25, 0.6))
# y caps (cases only; bands squished at the cap): a and b's Africa row, b's
# region rows (shared scale).
F3_CASES_YMAX <- c(africa = 45000, region = 30000)
# Layout (heights relative to panel a; ONI strips in mm): bottom (b / c) rows,
# the scope/arm legend row, b : c widths, c's arm padding (larger = forecast and
# counterfactual closer), the legends above c lifted by this fraction of their
# cell, the scope/arm legend's left shift (centres it under b's panels: the
# legend cell spans b's y-axis gutter but not its strips), and the "c" tag's
# height in its cell.
F3_HEIGHTS       <- c(a = 1, bottom = 2, legend = 0.12 * 1.6)
F3_WIDTHS_BC     <- c(1, 1)
F3_C_ARM_EXPAND  <- 1.2
F3_LEGEND_LIFT   <- 0.25
F3_LEGEND_SHIFT_MM <- 13
F3_TAG_C_Y       <- 0.12
F3_HEIGHT_MM     <- 200
# Compact legend theme for the cell above c.
F3_LEGEND_THEME <- theme(legend.position = "bottom", legend.justification = "left",
                         legend.box.just = "left",
                         legend.text = element_text(size = BASE_SIZE - 1),
                         legend.title = element_text(size = BASE_SIZE - 1),
                         legend.key.height = unit(2.6, "mm"), legend.key.width = unit(5.5, "mm"),
                         legend.spacing.y = unit(0.2, "mm"), legend.spacing.x = unit(2.5, "mm"),
                         legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0))

f3_scope_levels <- function(groups) c("Africa", unname(GROUP_SHORT_FIG[groups]))
f3_scope_pal    <- function(groups) c(Africa = F3_AFRICA_COL,
                                      setNames(unname(GROUP_COLOURS[groups]), GROUP_SHORT_FIG[groups]))
#' Arm labels used in legends.
f3_arm_lab <- function(aset, stitched = FALSE)
  c(forecast = if (stitched) "Factual (smoothing to handoff, forecast after)" else "Forecast",
    counterfactual = ps_cf_label(aset, capital = TRUE))

#' TRUE when the attribution artefact's factual series is stitched from the
#' smoothing (the arms launch at different months, 09's results_arm_window()).
f3_stitched <- function(f23) length(f23$stitch_months %||% NULL) > 0L

#' Quantile columns of a region / Africa ribbon for one metric -> t / q10..q90,
#' where q10 / q90 are the 2.5% / 97.5% columns: every outer band is 95%.
f3_qcols <- function(df, prefix) df %>%
  dplyr::transmute(t = month, q10 = .data[[paste0(prefix, "_q025")]], q25 = .data[[paste0(prefix, "_q25")]],
                   q50 = .data[[paste0(prefix, "_q50")]], q75 = .data[[paste0(prefix, "_q75")]],
                   q90 = .data[[paste0(prefix, "_q975")]])

#' Panel a: Africa smoothing (masked / unmasked) + observed + the free-running
#' forecast, up to the analogue handoff, on the ENSO-phase background.
#' `ymax`: y cap (NULL = none).
fig3_panel_a <- function(sm_afr, fc_df, metric, x_limits, vlines, ymax = NULL, phases = NULL,
                         legend_pos = "bottom", legend_nrow = 1L) {
  lv <- c("Smoothing (coverage-masked)", "Smoothing (unmasked)", "Observed", "Forecast")
  pal <- c(F3_MASK_COL, F3_AFRICA_COL, "black", F3_AFRICA_COL); names(pal) <- lv
  lty <- c("solid", "solid", "blank", F3_ARM_LTY[["forecast"]]); names(lty) <- lv
  sm_q <- function(r) r %>% dplyr::transmute(t = t_start, q10 = q025, q25, q50, q75, q90 = q975)
  long <- dplyr::bind_rows(
    sm_q(sm_afr$ribbon_masked)   %>% dplyr::mutate(series = lv[1]),
    sm_q(sm_afr$ribbon_unmasked) %>% dplyr::mutate(series = lv[2]),
    f3_qcols(fc_df, metric) %>% dplyr::mutate(series = lv[4])) %>%
    dplyr::filter(t <= x_limits[2]) %>%
    dplyr::mutate(series = factor(series, levels = lv))
  obs <- sm_afr$obs %>% dplyr::transmute(t = t_start, value, series = factor(lv[3], levels = lv))
  ggplot() +
    ps_phase_layer(phases) +
    geom_ribbon(data = long, aes(t, ymin = q10, ymax = q90, fill = series), alpha = 0.15) +
    geom_ribbon(data = long, aes(t, ymin = q25, ymax = q75, fill = series), alpha = 0.30) +
    geom_line(data = long, aes(t, q50, colour = series, linetype = series), linewidth = F3_LWD$a) +
    geom_point(data = obs, aes(t, value, colour = series), shape = 21, fill = "white",
               size = 0.9, stroke = 0.3) +
    geom_vline(xintercept = vlines, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    scale_fill_manual(values = pal, breaks = lv, name = NULL, drop = FALSE, guide = "none") +
    scale_linetype_manual(values = lty, breaks = lv, name = NULL, drop = FALSE, guide = "none") +
    scale_colour_manual(values = pal, breaks = lv, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(
      linetype = unname(lty), shape = c(NA, NA, 21, NA),
      fill = c(NA, NA, "white", NA), size = c(NA, NA, 1.2, NA), stroke = 0.4), nrow = legend_nrow,
      byrow = TRUE)) +
    scale_x_date(limits = x_limits, date_breaks = "2 years", date_labels = "%Y", expand = c(0, 0)) +
    scale_y_continuous(labels = PS_LABEL_K) +
    coord_cartesian(ylim = if (!is.null(ymax)) c(0, ymax)) +
    labs(x = NULL, y = METRIC_LABELS[metric]) +
    theme_ms() +
    theme(legend.position = legend_pos, legend.justification = "left",
          legend.key.width = unit(6, "mm"), legend.key.height = unit(3, "mm"))
}

#' Panel b: Africa + regions from the analogue handoff, forecast (scope colour,
#' dot-dashed) vs counterfactual (F3_CF_COL, dashed), median + 50% + 95% bands.
#' Africa on its own y scale, the regions sharing one (free_y + a blank point
#' at their common maximum); `africa_ymax` / `region_ymax` cap them (bands
#' squished; NULL = none).
#' When the factual series is stitched, [attrib_start, factual_start) is shaded
#' "smoothing-informed" and the factual handoff gets a dashed vline.
fig3_panel_b <- function(f23, metric, aset, groups, x_limits, africa_ymax = NULL, region_ymax = NULL) {
  stitched <- f3_stitched(f23)
  al <- f3_arm_lab(aset, stitched)
  reg <- f23$region_ribbon %>% dplyr::filter(enso_group %in% groups) %>%
    dplyr::mutate(scope = unname(GROUP_SHORT_FIG[enso_group]))
  afr <- f23$africa_ribbon %>% dplyr::mutate(scope = "Africa")
  lv <- f3_scope_levels(groups)
  d <- dplyr::bind_rows(reg, afr) %>%
    dplyr::filter(run == "forecast" | analog_set == aset, month >= f23$attrib_start) %>%
    dplyr::mutate(arm = factor(as.character(run), levels = names(F3_ARM_LTY)),
                  scope = factor(scope, levels = lv),
                  key = ifelse(arm == "forecast", as.character(scope), "cf"))
  d <- dplyr::bind_cols(d, f3_qcols(d, metric) %>% dplyr::select(-t))
  pal <- c(f3_scope_pal(groups), cf = F3_CF_COL)
  region_max <- max(d$q90[d$scope != "Africa"], na.rm = TRUE)
  if (!is.null(region_ymax)) {
    reg_rows <- d$scope != "Africa"
    d$q90[reg_rows] <- pmin(d$q90[reg_rows], region_ymax)
    d$q75[reg_rows] <- pmin(d$q75[reg_rows], region_ymax)
    region_max <- region_ymax
  }
  blank <- data.frame(scope = factor(setdiff(levels(d$scope), "Africa"), levels = levels(d$scope)),
                      month = min(d$month), y = region_max)
  if (!is.null(africa_ymax)) {   # cap Africa's free row like panel a: squish its band, pad to the cap
    afr_rows <- d$scope == "Africa"
    d$q90[afr_rows] <- pmin(d$q90[afr_rows], africa_ymax)
    d$q75[afr_rows] <- pmin(d$q75[afr_rows], africa_ymax)
    blank <- rbind(blank, data.frame(scope = factor("Africa", levels = levels(d$scope)),
                                     month = min(d$month), y = africa_ymax))
  }
  stitch_layers <- if (stitched) list(
    annotate("rect", xmin = as.Date(f23$attrib_start), xmax = as.Date(f23$factual_start),
             ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5),
    annotate("text", x = as.Date(f23$attrib_start), y = Inf, label = "smoothing-informed",
             hjust = -0.05, vjust = 1.4, size = (BASE_SIZE - 2) / ggplot2::.pt, colour = "grey35"),
    geom_vline(xintercept = as.Date(f23$factual_start), linetype = "dashed", colour = "grey45",
               linewidth = 0.3)) else NULL
  ggplot(d, aes(month)) +
    stitch_layers +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = key), alpha = 0.18) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = key), alpha = 0.3) +
    geom_line(aes(y = q50, colour = key, linetype = arm), linewidth = F3_LWD$b) +
    geom_blank(data = blank, aes(month, y)) +
    facet_grid(scope ~ ., switch = "y", scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_linetype_manual(values = F3_ARM_LTY, labels = al, name = NULL,
                          guide = guide_legend(override.aes = list(colour = c("grey10", F3_CF_COL)))) +
    scale_x_date(limits = x_limits, date_breaks = "6 months", date_labels = "%b\n%Y", expand = c(0.01, 0)) +
    scale_y_continuous(labels = PS_LABEL_K) +
    labs(x = NULL, y = METRIC_LABELS[metric]) +
    theme_ms() +
    # Strips outside the axis. patchwork then reserves this strip column
    # BETWEEN axis and panel in every strip-less plot of the same column, which
    # detaches their tick labels -- so a and the ONI strips get a blank strip of
    # the same width (blank_left_strip()) instead of none.
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "plain"),
          legend.position = "bottom", legend.key.width = unit(6, "mm"))
}

#' Panel c: cumulative suspected cases and outbreak onsets over the analogue
#' window (phase-09 attribution `cum_totals`: paired replicate totals, median +
#' 50% + 95%), forecast vs counterfactual. Rows = scope (aligned with b),
#' columns = metric, each on its own log10 x range with 1-3-10 ticks; metric
#' titles under the axes (on top they push b's panels down).
fig3_panel_c <- function(totals, aset, groups, attrib_start, factual_start = attrib_start) {
  stitched <- as.Date(factual_start) > as.Date(attrib_start)
  al <- f3_arm_lab(aset, stitched)
  lv <- f3_scope_levels(groups)
  d <- totals %>%
    dplyr::filter(analog_set == aset) %>%
    dplyr::mutate(scope = ifelse(scope == "Africa", "Africa", unname(GROUP_SHORT_FIG[scope]))) %>%
    dplyr::filter(scope %in% lv) %>%
    dplyr::mutate(scope = factor(scope, levels = lv),
                  metric = factor(metric, levels = c("cases", "onsets"),
                                  labels = c("Suspected cases", "Outbreak onsets")),
                  arm = factor(run, levels = names(F3_ARM_LTY)),
                  key = ifelse(arm == "forecast", as.character(scope), "cf"))
  pal <- c(f3_scope_pal(groups), cf = F3_CF_COL)
  ggplot(d, aes(x = med, y = forcats::fct_rev(arm), colour = key, shape = arm)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = F3_LWD$c[1]) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = F3_LWD$c[2]) +
    geom_point(fill = "white", size = 1.8, stroke = 0.5) +
    scale_y_discrete(expand = expansion(add = F3_C_ARM_EXPAND)) +
    ggh4x::facet_grid2(scope ~ metric, scales = "free_x", switch = "x") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_shape_manual(values = F3_ARM_SHAPE, labels = al, name = NULL,
                       guide = guide_legend(override.aes = list(colour = c("grey10", F3_CF_COL)))) +
    scale_x_log10(labels = PS_LABEL_K, breaks = c(1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5)) +
    labs(x = paste0(sprintf("Cumulative total (%s to horizon", format(as.Date(attrib_start), "%b %Y")),
                    if (stitched) sprintf("; factual from the smoothing to %s",
                                          format(as.Date(factual_start), "%b %Y")) else "", ")"),
         y = NULL) +
    theme_ms() +
    theme(legend.position = "bottom", strip.text.y = element_blank(),
          strip.text.x = element_text(size = BASE_SIZE - 0.5),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = BASE_SIZE - 1.5),
          panel.grid.major.y = element_blank(),
          strip.switch.pad.grid = unit(0, "mm"),
          axis.title.x = element_text(margin = margin(t = 0.3, unit = "mm")),
          strip.placement = "outside",
          strip.text.x.bottom = element_text(size = BASE_SIZE - 0.5, margin = margin(t = 0.5, unit = "mm")))
}

#' ONI strip sharing its panel's x axis (Figures 1a, 3a, 3b). mode =
#' "observed": the observed record over x_limits, on the phase background. mode = "forecast":
#' observed (black) + CFSv2 forecast (F3_ONI_COL, dot-dashed, 95% band) +
#' counterfactual analogue mean (F3_CF_COL, dashed, 95% band), cut at the last
#' month whose forecast still carries CFSv2 information (the neutral 0 fallback
#' beyond the CFSv2 span is not drawn).
fig3_oni_strip <- function(oni, aset, x_limits, mode = c("observed", "forecast"), phases = NULL) {
  mode <- match.arg(mode)
  obs <- oni$obs_full %>% dplyr::filter(date >= x_limits[1], date <= x_limits[2])
  p <- ggplot() +
    ps_phase_layer(phases) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25)
  if (mode == "forecast") {
    fc_end <- max(oni$env$date[oni$env$arm == "forecast" & abs(oni$env$mean) > 1e-9])
    env <- oni$env %>%
      dplyr::filter(arm == "forecast" | analog_set == aset, date >= x_limits[1], date <= min(fc_end, x_limits[2])) %>%
      dplyr::mutate(arm = factor(arm, levels = names(F3_ARM_LTY)))
    pal <- c(forecast = F3_ONI_COL, counterfactual = F3_CF_COL)
    p <- p +
      geom_ribbon(data = env, aes(date, ymin = lo95, ymax = hi95, fill = arm), alpha = 0.18) +
      geom_line(data = env, aes(date, mean, colour = arm, linetype = arm), linewidth = F3_LWD$b) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_colour_manual(values = pal, guide = "none") +
      scale_linetype_manual(values = F3_ARM_LTY, guide = "none")
  }
  # Same x breaks (hence gridlines) as the panel below: 2 years above (a),
  # 6 months above (b); no x labels on either (the panel's axis serves both).
  p <- p +
    geom_line(data = obs, aes(date, enso_index), colour = "black", linewidth = F3_LWD$b) +
    scale_y_continuous(breaks = c(-2, 0, 2)) +
    coord_cartesian(ylim = ONI_YLIM) +
    labs(x = NULL, y = "ONI") +
    theme_ms() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length.x = unit(0, "mm"),
          plot.margin = margin(1.5, 1.5, 0, 1.5, unit = "mm"))
  if (mode == "forecast") p + scale_x_date(limits = x_limits, date_breaks = "6 months", expand = c(0.01, 0))
  else                    p + scale_x_date(limits = x_limits, date_breaks = "2 years",  expand = c(0, 0))
}

#' The merged b / c legend: scope colours (filled squares) on one row, the two
#' arms as line + point on the next, centred under b's panels.
fig3_legend_bc <- function(groups, aset, stitched = FALSE) {
  al <- f3_arm_lab(aset, stitched)
  lv <- f3_scope_levels(groups)
  d_scope <- data.frame(scope = factor(lv, levels = lv), x = 1, y = 1)
  d_arm <- data.frame(arm = factor(names(F3_ARM_LTY), levels = names(F3_ARM_LTY)), x = 1, y = 1)
  src <- ggplot() +
    geom_point(data = d_scope, aes(x, y, colour = scope), shape = 15, size = 2.5) +
    geom_line(data = d_arm, aes(x, y, linetype = arm, group = arm), linewidth = 0.45) +
    geom_point(data = d_arm, aes(x, y, shape = arm), size = 1.8) +
    scale_colour_manual(values = f3_scope_pal(groups), name = NULL,
                        guide = guide_legend(order = 1, nrow = 1)) +
    scale_linetype_manual(values = F3_ARM_LTY, labels = al, name = NULL) +
    scale_shape_manual(values = F3_ARM_SHAPE, labels = al, name = NULL) +
    guides(linetype = guide_legend(order = 2, override.aes = list(colour = c("grey10", F3_CF_COL))),
           shape    = guide_legend(order = 2, override.aes = list(colour = c("grey10", F3_CF_COL)))) +
    theme_ms() +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(7, "mm"), legend.spacing.x = unit(4, "mm"),
          legend.spacing.y = unit(0.5, "mm")) +
    ps_legend_theme() +
    theme(legend.justification = "center", legend.box.just = "center",
          legend.box.margin = margin(l = F3_LEGEND_SHIFT_MM, unit = "mm"))
  cowplot::get_legend(src)
}

#' ENSO-state legend alone (Figure 1's palette / alpha), compact.
fig3_phase_legend <- function() {
  d <- data.frame(ph = factor(names(PS_ENSO_PHASE_COL), levels = names(PS_ENSO_PHASE_COL)), x = 1, y = 1)
  cowplot::get_legend(
    ggplot(d, aes(x, y, fill = ph)) + geom_tile(alpha = PS_ENSO_PHASE_ALPHA) +
      scale_fill_manual(values = PS_ENSO_PHASE_COL, labels = PS_ENSO_PHASE_LAB, name = "ENSO state",
                        guide = guide_legend(override.aes = list(alpha = 3 * PS_ENSO_PHASE_ALPHA), nrow = 1)) +
      theme_ms() + F3_LEGEND_THEME + ps_legend_theme())
}

#' Figure 3 for one metric (a / b) x one counterfactual set:
#'   a  observed-ONI strip over Africa's smoothing (masked / unmasked),
#'      observed series and the free-running forecast to the analogue handoff;
#'   b  ONI forecast strip (CFSv2 vs the counterfactual's analogues) over
#'      Africa + the four regions from the handoff, forecast vs counterfactual;
#'   c  cumulative cases and onsets over the analogue window per scope.
#' Panel a's series legend and the ENSO-state legend sit in the white cell
#' above c (which carries c's tag); scope colours + arms under b.
build_fig3 <- function(metric = "cases", aset = MAIN_SET) {
  message(sprintf("[10] Figure 3: metric=%s set=%s", metric, aset))
  f23 <- read_results("attribution", sets = SETS)
  sm  <- read_results("smoothing")
  oni <- read_results("oni", sets = SETS)
  groups <- GROUP_ORDER[GROUP_ORDER %in% f23$groups]
  # Panel a ends at the attribution start (09's results_arm_window(): the
  # earlier arm launch when the factual series is stitched from the smoothing).
  handoff_month <- as.Date(f23$attrib_start)
  x_limits_a <- c(as.Date("2011-01-01"), as.Date(f23$attrib_start))
  x_limits_b <- c(as.Date(f23$attrib_start), max(f23$region_ribbon$month))
  caps <- if (identical(metric, "cases")) F3_CASES_YMAX else c(africa = NULL, region = NULL)
  ph_a <- enso_phase_rects(oni, x_limits_a)

  sm_afr <- sm$by_scope_metric[[paste("Africa", metric, sep = "|")]]
  fc_afr <- dplyr::filter(f23$africa_ribbon, run == "forecast")
  pa <- fig3_panel_a(sm_afr, fc_afr, metric, x_limits_a, handoff_month,
                     ymax = caps[["africa"]], phases = ph_a)
  pb <- fig3_panel_b(f23, metric, aset, groups, x_limits_b,
                     africa_ymax = caps[["africa"]], region_ymax = caps[["region"]])
  pc <- fig3_panel_c(f23$cum_totals, aset, groups, f23$attrib_start,
                     f23$factual_start %||% f23$attrib_start)

  # Two ONI strips sharing their panel's x range; patchwork aligns the panel
  # areas within each grid row, so a strip sits exactly over its panel.
  oni_a <- blank_left_strip(fig3_oni_strip(oni, aset, x_limits_a, mode = "observed", phases = ph_a))
  oni_b <- blank_left_strip(fig3_oni_strip(oni, aset, x_limits_b, mode = "forecast"))
  pa <- blank_left_strip(pa) + ps_legend_theme()
  legend_bc <- fig3_legend_bc(groups, aset, f3_stitched(f23))
  pb <- pb + theme(legend.position = "none")
  pc <- pc + theme(legend.position = "none")
  # Panel a's series legend (two rows) + the ENSO-state legend, compact, in
  # the white cell beside b's strip, which also carries c's tag at its
  # bottom-left (a tag on pc itself reserves space above pc's panels, and b's
  # panels, in the same grid row, get pushed down with them). Drawn at full
  # cell height but lifted by F3_LEGEND_LIFT (clip off), towards a's axis.
  leg_a <- cowplot::get_legend(
    fig3_panel_a(sm_afr, fc_afr, metric, x_limits_a, handoff_month, legend_nrow = 2L) +
      F3_LEGEND_THEME + ps_legend_theme())
  leg_top <- cowplot::plot_grid(leg_a, fig3_phase_legend(), ncol = 1, align = "v", axis = "l",
                                rel_heights = c(1.5, 1))
  pa <- pa + theme(legend.position = "none")
  tag_c <- ggplot() + theme_void() +
    geom_blank(data = data.frame(x = 0:1, y = 0:1), aes(x, y)) +
    annotation_custom(cowplot::as_grob(leg_top), xmin = 0, xmax = 1,
                      ymin = F3_LEGEND_LIFT, ymax = 1 + F3_LEGEND_LIFT) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme(plot.tag.position = c(0, 0.02), plot.margin = margin(0, 0, 0, 4, unit = "mm"))
  fig <- patchwork::wrap_plots(A = oni_a, B = pa, C = oni_b, D = tag_c, E = pb, F = pc,
                               G = patchwork::wrap_elements(full = legend_bc, ignore_tag = TRUE),
                               design = "AA\nBB\nCD\nEF\nG#",
                               heights = grid::unit(c(PS_ONI_STRIP_MM, F3_HEIGHTS[["a"]], PS_ONI_STRIP_MM,
                                                      F3_HEIGHTS[["bottom"]], F3_HEIGHTS[["legend"]]),
                                                    c("mm", "null", "mm", "null", "null")),
                               widths = F3_WIDTHS_BC) +
    patchwork::plot_annotation(tag_levels = list(c("a", "", "b", "c", ""))) & TAG_THEME
  # An empty "" tag still reserves a tag row above its plot (a gap under both
  # ONI strips): drop the tag grob on the untagged plots; "c" low in its cell,
  # clear of c's top border.
  for (k in c(2L, 5L, 6L)) fig[[k]] <- fig[[k]] + theme(plot.tag = element_blank())
  fig[[4]] <- fig[[4]] + theme(plot.tag.position = c(-0.03, F3_TAG_C_Y))
  # ~1 mm between each ONI strip and the panel under it (the strips carry the
  # tags: a = ONI + cases, b = ONI + regions).
  for (k in c(2L, 5L, 6L)) fig[[k]] <- fig[[k]] + theme(plot.margin = margin(1, 1.5, 1.5, 1.5, unit = "mm"))
  # c's axis row is aligned with b's two-line date labels, which left a gap
  # above c's metric titles: free c's bottom spacing.
  fig[[6]] <- patchwork::free(fig[[6]], type = "space", side = "b")
  # "ONI" next to its own tick labels in a and b (free labels, not pushed out
  # to the width of the cases-axis titles).
  for (k in c(1L, 3L))
    fig[[k]] <- patchwork::free(fig[[k]] + theme(axis.title.y = element_text(margin = margin(r = 0.5, unit = "mm"))),
                                type = "label", side = "l")

  is_main <- identical(metric, "cases") && identical(aset, MAIN_SET)
  stem   <- if (is_main) "fig3_smoothing_forecast_counterfactual"
  else sprintf("fig3_smoothing_forecast_counterfactual_%s_%s", metric, aset)
  subdir <- if (is_main) NULL else "supplementary"
  save_figure(fig, stem, height_mm = F3_HEIGHT_MM, subdir = subdir)
}

# ===========================================================================
# Section 9d: Figure 4 builder — attribution (paired ratio + ADM1 maps)
# ===========================================================================
#' Blanks the x-axis on the two upper Figure 4a sub-panels; all three share one
#' x scale per column, so they stay aligned without repeating the axis.
.fig4a_strip_x <- ggplot2::theme(axis.title.x = element_blank(),
                                 axis.text.x  = element_blank(),
                                 axis.ticks.x = element_blank())

#' Figure 4a: forest stack of the paired forecast / counterfactual ratio,
#' continent -> ENSO region -> country, in two columns (suspected cases |
#' outbreak onsets) sharing the country rows (ordered by the cases ratio), one
#' set. Colour = attribution direction per metric (PROB_THRESHOLD), as Figure
#' 4b's maps; dot area = median forecast as a share of that metric's Africa
#' total (one size scale cannot hold ~200K cases and ~300 onsets). Each column
#' has its own log ratio axis (cases 0.02-50, onsets FIG4_ONSET_LIM).
#' Returns list(stack, legend, n_rows); build_fig4() sizes the figure from
#' `n_rows`.
build_fig4_panel_a <- function(ratio_stats, aset, ratio_label = "ENSO") {
  d <- ratio_stats %>% dplyr::filter(metric %in% c("cases", "onsets"), analog_set == aset) %>%
    dplyr::group_by(metric) %>%
    dplyr::mutate(med_forecast = med_forecast / med_forecast[level == "continent"][1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(metric = factor(metric, levels = c("cases", "onsets"),
                                  labels = c("Suspected cases", "Outbreak onsets")))
  # ratio_label (09's attribution_ratio_label()): "ENSO" when both arms launch
  # together; otherwise the ratio also carries the post-handoff observations.
  x_lab <- sprintf("Forecast / %s ratio (cumulative, analogue window%s)", ps_cf_label(aset),
                   if (identical(ratio_label, "ENSO")) "" else paste0("; ", ratio_label))
  region_keys <- intersect(GROUP_ORDER, unique(d$region))
  dir_name <- sprintf("Attribution\n(P > %.0f%% chance of a > %.0f%% shift\nvs. %s)",
                      100 * PROB_THRESHOLD, 100 * PROB_THRESHOLD, ps_cf_label(aset))
  # Sub-panel margins (mm, top/bottom: 0.3 + 0.3 = the 0.6 mm gap between the
  # region facets), no tick length on the blanked x axes of the upper two
  # blocks (a blank tick still reserves its length), smaller y axis text.
  sub_margin <- theme(axis.text.y = element_text(size = BASE_SIZE - 1.5, colour = "grey20"),
                      plot.margin = margin(0.3, 1.5, 0.3, 1.5, unit = "mm"),
                      axis.ticks.length.x = unit(0, "mm"),
                      panel.grid.major.y = element_blank())
  # Figure 3's left strips (plain label outside the axis, no box) ...
  strip_left <- theme(strip.placement = "outside", strip.background = element_blank(),
                      strip.text.y.left = element_text(angle = 90, face = "bold", size = BASE_SIZE))
  # ... and a transparent one on the strip-less Africa / region blocks, so
  # patchwork aligns their axes with the country block's instead of inserting
  # an empty strip column between axis and panel. (A real `blk` column, not a
  # constant: a constant row var with a data column var leaves the data-less
  # vline layers in an NA panel.)
  blank_strip <- list(ggh4x::facet_grid2(rows = vars(blk), cols = vars(metric), switch = "y", scales = "free_x"),
                      strip_left,
                      theme(strip.text.y.left = element_text(angle = 90, face = "bold",
                                                             size = BASE_SIZE, colour = "transparent")))
  # One scale for the lines' colour and the dots' fill, so they share a legend.
  scale_direction <- scale_colour_manual(values = DIRECTION_COLOURS, name = dir_name,
                                         labels = c(increase = "Increase", "no change" = "No change",
                                                    decrease = "Decrease"), drop = FALSE,
                                         aesthetics = c("colour", "fill"))
  scale_dot_size <- scale_size_continuous(
    range = FIG4_DOT_RANGE, limits = c(0, 1),
    breaks = c(0.01, 0.1, 0.5, 1), labels = scales::label_percent(accuracy = 1),
    name = "Median forecast\n(% of Africa total)")
  # 95% interval (thin), IQR (thick), median dot (area = forecast share, white
  # border); lines and dot drawn separately so size only scales the dot.
  forest_layers <- list(
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.2),
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.12),
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 0.45),
    geom_point(aes(size = med_forecast, fill = direction), shape = 21,
               colour = "white", stroke = 0.25),
    scale_direction, scale_dot_size)
  # Per-column log ratio axes (squished at the limits); added after each
  # block's facet (ggh4x needs the facet first).
  ratio_scales <- ggh4x::facetted_pos_scales(x = list(
    metric == "Suspected cases" ~ scale_x_log10(limits = c(0.02, 50), oob = scales::squish,
                                                breaks = c(0.02, 0.1, 0.5, 1, 2, 10, 50),
                                                labels = c("0.02", "0.1", "0.5", "1", "2", "10", "50")),
    metric == "Outbreak onsets" ~ scale_x_log10(limits = FIG4_ONSET_LIM, oob = scales::squish,
                                                breaks = c(0.05, 0.2, 0.5, 1, 2, 5, 20),
                                                labels = c("0.05", "0.2", "0.5", "1", "2", "5", "20"))))

  d_africa <- d %>% dplyr::filter(level == "continent") %>%
    dplyr::mutate(scope = factor(scope, levels = "Africa"), blk = "x")
  p_africa <- ggplot(d_africa, aes(x = med, y = scope, colour = direction)) +
    forest_layers + labs(x = NULL, y = NULL) + theme_ms() + .fig4a_strip_x +
    theme(legend.position = "none") + blank_strip + ratio_scales + sub_margin

  d_region <- d %>% dplyr::filter(level == "region") %>%
    dplyr::mutate(region = factor(region, levels = region_keys),
                  scope  = factor(group_short_fig(scope), levels = rev(unname(GROUP_SHORT_FIG[region_keys]))),
                  blk = "x")
  # Metric titles on the Africa block only.
  p_region <- ggplot(d_region, aes(x = med, y = scope, colour = direction)) +
    forest_layers + labs(x = NULL, y = NULL) + theme_ms() + .fig4a_strip_x +
    theme(legend.position = "none") + blank_strip + ratio_scales + sub_margin +
    theme(strip.text.x = element_blank())

  # Countries ordered ascending by the cases ratio -> highest at the top of its
  # region's facet, same order in both columns.
  d_country <- d %>% dplyr::filter(level == "country")
  ord_src <- dplyr::filter(d_country, metric == "Suspected cases")
  d_country <- d_country %>%
    dplyr::mutate(scope  = factor(scope, levels = unique(dplyr::arrange(ord_src, med)$scope)),
                  region = factor(region, levels = region_keys))
  # patchwork's `guides = "collect"` is flaky on this 3-panel stack, so the
  # legend is built once from p_country, extracted with cowplot::get_legend(),
  # and placed by build_fig4() with every panel's own legend suppressed.
  # Legends stacked vertically and left-aligned.
  p_country_legendsrc <- ggplot(d_country, aes(x = med, y = scope, colour = direction)) +
    forest_layers +
    ggh4x::facet_grid2(region ~ metric, scales = "free", space = "free_y", switch = "y",
                       strip = ggh4x::strip_vanilla(), labeller = ggplot2::as_labeller(GROUP_SHORT_FIG)) +
    ratio_scales +
    labs(x = x_lab, y = NULL) +
    # Identical colour and fill guides so ggplot merges them into one legend;
    # size keys filled dark, since the dots' own fill comes from `direction`.
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2)),
           fill   = guide_legend(order = 1, override.aes = list(size = 2)),
           size   = guide_legend(order = 2, override.aes = list(fill = "grey30"))) +
    theme_ms() +
    theme(strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title.position = "left", legend.box = "vertical",
          legend.box.just = "left", legend.justification = "left",
          legend.spacing.y = unit(1, "mm")) +
    strip_left +
    # Region facets 0.6 mm apart (theme default 5.5 pt ~ 1.9 mm); the two
    # metric columns 1.2 mm apart; legend text as the maps' legend.
    theme(panel.spacing.y = unit(0.6, "mm"),
          legend.title = element_text(size = BASE_SIZE - 1.5),
          legend.text  = element_text(size = BASE_SIZE - 2),
          legend.key.height = unit(2, "mm"), legend.key.width = unit(3.5, "mm"),
          strip.text.x = element_blank(), panel.spacing.x = unit(1.2, "mm"))
  legend_grob <- cowplot::get_legend(p_country_legendsrc)
  p_country <- p_country_legendsrc + theme(legend.position = "none") + sub_margin +
    theme(axis.ticks.length.x = unit(0.8, "mm"),   # theme_ms()'s, kept on the bottom block
          axis.text.x = element_text(size = BASE_SIZE - 1.5, colour = "grey20"),
          axis.title.x = element_text(size = BASE_SIZE - 1))

  n_rows <- dplyr::n_distinct(d_africa$scope) + dplyr::n_distinct(d_region$scope) +
    dplyr::n_distinct(d_country$scope)
  stack <- patchwork::wrap_plots(p_africa, p_region, p_country, ncol = 1,
                                 heights = c(2, length(region_keys) + 1, dplyr::n_distinct(d_country$scope)))
  list(stack = stack, legend = legend_grob, n_rows = n_rows)
}

#' Figure 4b: ADM1 map of the signed probability of a forecast-vs-
#' counterfactual shift, one metric x one set, titled with the metric.
#' Basemap, greys (PS_GREY, keyed; low-signal ADM1s, per LOW_SIGNAL, are not
#' drawn and show the "Modelled, no ADM1 estimate" grey), fill colours and
#' scale bar as Figure 2. Significant ADM1s (per PROB_THRESHOLD) get a dashed
#' black outline. `afro_sf` comes without the small islands
#' (drop_small_islands()), and the map is cropped to its bbox so the panel
#' carries no extra white. `maps`: load_africa_maps(); `basemap`:
#' load_africa_basemap(); `show_legend`: the legend (bottom, flat colourbar) or
#' none; `show_scale`: scale bar + north arrow.
build_fig4_panel_b <- function(prob_df, metric, aset, maps, basemap, ob_countries, title,
                               show_legend = TRUE, show_scale = TRUE) {
  afro_sf <- maps$afro_sf
  d <- prob_df %>% dplyr::filter(metric == !!metric, analog_set == aset, !is.na(signed_prob))
  map_sf <- maps$adm1_sf %>% dplyr::inner_join(d, by = c("GID_1" = "adm1_unit"))
  # Significant ADM1 outlines as one deduplicated line set: drawn per polygon,
  # a border shared by two significant ADM1s gets two out-of-phase dash
  # patterns that read as a solid line.
  # Planar (GEOS) union: s2's union of lines keeps coincident edges.
  s2_was <- sf::sf_use_s2(FALSE)
  sig_lines <- dplyr::filter(map_sf, significant) %>%
    sf::st_geometry() %>% sf::st_boundary() %>% sf::st_union() %>% sf::st_line_merge() %>%
    sf::st_sf(geometry = .)
  sf::sf_use_s2(s2_was)
  bb <- sf::st_bbox(afro_sf)

  sig_lab <- sprintf("(P > %.0f%% chance of\na > %.0f%% shift vs.\n%s)", 100 * PROB_THRESHOLD,
                     100 * PROB_THRESHOLD, ps_cf_label(aset))
  ggplot() +
    ps_country_grey_layer(ps_country_status(afro_sf, ob_countries), legend_order = 3) +
    geom_sf(data = map_sf, aes(fill = signed_prob), colour = "white", linewidth = 0.005) +
    map_overlay_layers(basemap, afro_sf) +
    geom_sf(data = sig_lines, aes(colour = "sig"), linewidth = 0.15, linetype = "dashed",
            key_glyph = "polygon") +
    # Figure 2's scale bar (1000 km) and arrow.
    { if (show_scale) map_scale_layers(arrow_pad_x_mm = 2) } +
    # Figure 2's red/blue (= Figure 4a's DIRECTION_COLOURS at +/-1).
    scale_fill_gradient2(
      midpoint = 0, low = "#2166ac", mid = "white", high = "#d73027",
      limits = c(-1, 1), oob = scales::squish,
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = c("1\ndecrease", "0.5", "0", "0.5", "1\nincrease"),
      name = sprintf("Probability of shift vs.\n%s", ps_cf_label(aset)),
      guide = guide_colourbar(direction = "vertical", title.position = "top", order = 1)) +
    scale_colour_manual(values = c(sig = "black"), labels = c(sig = sig_lab), name = NULL,
                        guide = guide_legend(order = 2,
                                             override.aes = list(fill = NA, linewidth = 0.3))) +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], expand = FALSE) +
    labs(title = title) +
    taxdat::map_theme() +
    ps_map_font() +
    theme(panel.border = element_blank(),
          plot.margin = margin(1, 1, 1, 2, unit = "mm"),
          text = element_text(size = BASE_SIZE - 1),
          legend.position = "inside",
          legend.position.inside = c(0.005, 0.005),
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.spacing.y = unit(0.8, "mm"),
          legend.key.height = unit(2.4, "mm"),
          legend.key.width  = unit(2, "mm"),
          legend.box.just = "left",
          legend.title = element_text(size = BASE_SIZE - 1.5),
          legend.text  = element_text(size = BASE_SIZE - 2),
          plot.title = element_text(size = BASE_SIZE - 0.5, hjust = 0.5,
                                    margin = margin(b = 0.5, unit = "mm"))) +
    { if (!show_legend) theme(legend.position = "none")
      else list(
        theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
              legend.justification = "left", legend.spacing.y = unit(0.6, "mm"),
              legend.title.position = "top", legend.key.height = unit(2.2, "mm"),
              legend.key.width = unit(3, "mm")),
        guides(fill = guide_colourbar(order = 1, direction = "horizontal", title.position = "top",
                                      theme = theme(legend.key.width = unit(30, "mm"),
                                                    legend.key.height = unit(2, "mm"))))) }
}

#' Figure 4 for one counterfactual set: forest stack (a) | the cases map above
#' the outbreak-onset map (b), with a's legend and the shared map legend in a
#' row underneath, so both maps get the full height of the forest stack.
build_fig4 <- function(aset = MAIN_SET) {
  message(sprintf("[10] Figure 4: set=%s", aset))
  f23 <- read_results("attribution", sets = SETS)
  maps <- load_africa_maps("Figure 4")
  basemap <- load_africa_basemap(maps$afro_sf)   # utils/figures_cb.R (lakes + rivers, st_intersects clip)
  ob_countries <- modelled_countries()

  pa <- build_fig4_panel_a(f23$ratio_stats, aset, f23$ratio_label %||% "ENSO")
  map_panel <- function(metric, title, legend) build_fig4_panel_b(
    f23$prob_df, metric, aset, maps, basemap, ob_countries, title,
    show_legend = legend, show_scale = !legend)
  map_cases  <- map_panel("cases", "Suspected cases", legend = FALSE)
  map_onsets <- map_panel("onsets", "Outbreak onsets", legend = TRUE)
  # The shared map legend, compact (one-line significance label, grey keys
  # stacked, shorter colourbar), placed under the maps' column.
  map_leg <- cowplot::get_legend(
    map_onsets +
      scale_colour_manual(values = c(sig = "black"),
                          labels = c(sig = sprintf("P > %.0f%% chance of a\n> %.0f%% shift",
                                                   100 * PROB_THRESHOLD, 100 * PROB_THRESHOLD)),
                          name = NULL, guide = guide_legend(order = 2, override.aes = list(fill = NA, linewidth = 0.3))) +
      theme(legend.box = "horizontal", legend.box.just = "top",
            legend.direction = "vertical",
            legend.spacing.x = unit(2, "mm"),
            legend.title = element_text(size = BASE_SIZE - 1.5),
            legend.text = element_text(size = BASE_SIZE - 2)) +
      guides(fill = guide_colourbar(order = 1, direction = "horizontal", title.position = "top",
                                    theme = theme(legend.key.width = unit(18, "mm"),
                                                  legend.key.height = unit(2, "mm")))))
  map_onsets <- map_onsets + theme(legend.position = "none")
  maps <- patchwork::wrap_plots(map_cases, map_onsets, ncol = 1)

  # Forest stack and legends laid out on their own (wrap_elements), otherwise
  # patchwork aligns the stack's blocks with the maps' rows and opens gaps.
  top <- patchwork::wrap_plots(patchwork::wrap_elements(full = pa$stack), maps, nrow = 1, widths = FIG4_WIDTHS)
  bottom <- patchwork::wrap_plots(patchwork::wrap_elements(full = pa$legend),
                                  patchwork::wrap_elements(full = map_leg), nrow = 1, widths = FIG4_LEG_WIDTHS)
  fig <- patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1, FIG4_LEG_ROW)) +
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "", "", ""))) & TAG_THEME

  is_main <- identical(aset, MAIN_SET)
  stem   <- if (is_main) "fig4_attribution" else sprintf("fig4_attribution_%s", aset)
  subdir <- if (is_main) NULL else "supplementary"
  # Height grows with the forest-plot row count so point-ranges stay legible.
  height_mm <- max(95, 20 + FIG4_ROW_MM * pa$n_rows)
  save_figure(fig, stem, height_mm = height_mm, subdir = subdir)
}

# ===========================================================================
# Section 9e: supplementary figures
# ===========================================================================
#' Supplementary — Figure 1 ERF at country granularity: population fixed effect
#' plus that country's own ENSO-slope RE, one facet per country, coloured by
#' ENSO group. Same deployed best-lag fit and x-axis as Figure 1's ERF grid.
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
    facet_wrap(~ country, ncol = 6) +
    scale_colour_manual(values = GROUP_COLOURS, name = "ENSO group", labels = GROUP_LABELS) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_y_log10() +
    coord_cartesian(xlim = ERF_XLIM) +
    labs(x = "Oceanic Ni\u00f1o Index (ONI)", y = sprintf("%s IRR", TRACK_LABEL[[track]]),
         title = sprintf("Country-level ENSO exposure-response (%s)", TRACK_LABEL[[track]]),
         subtitle = "population + region + country slope RE; same deployed best-lag fit per ENSO group as Figure 1") +
    theme_ms() +
    theme(legend.position = "bottom",
          strip.text  = element_text(size = BASE_SIZE - 1),
          axis.text.x = element_text(angle = 45, hjust = 1))
  stem <- if (identical(track, "occurrence")) "fig1_supp_country_erf" else sprintf("fig1_supp_country_erf_%s", track)
  save_figure(fig, stem,
              height_mm = ceiling(n_countries / 6L) * 26 + 25,
              width_mm = 210, subdir = "supplementary")
}

#' Supplementary — region/country parameter caterpillar for one track. Reads the
#' derive_track_params() artefact (09) and plots, for every NON-ENSO parameter
#' family in SUPP_PARAM_FAMILIES[[track]], the country-level posterior ratio
#' (exp scale, 95% CrI) with the per-group population estimate overlaid as a
#' distinct marker. Columns are the individual terms (ordered by family), rows
#' the ENSO groups; colour encodes the group, shape the tier. Fixed-only terms
#' (spillover, ttlo, seasonal, month — no country RE) show only the population
#' marker.
build_fig_supp_track_params <- function(track = "occurrence") {
  fams <- SUPP_PARAM_FAMILIES[[track]]
  if (is.null(fams)) { warning(sprintf("[10] no SUPP_PARAM_FAMILIES for track %s", track)); return(invisible(NULL)) }
  d <- read_results("track_params", track = track) %>% dplyr::filter(family %in% fams)
  if (!nrow(d)) { warning(sprintf("[10] track params (%s): nothing to plot after family filter", track)); return(invisible(NULL)) }

  term_tbl <- d %>% dplyr::distinct(family, term) %>%
    dplyr::mutate(family = factor(family, levels = fams)) %>%
    dplyr::arrange(family, term)
  country_levels <- sort(unique(d$level[d$tier == "country"]))
  pd <- d %>% dplyr::mutate(
    enso_group = factor(enso_group, levels = GROUP_ORDER),
    term       = factor(term, levels = term_tbl$term),
    ylab       = factor(ifelse(tier == "population", "Population", level),
                        levels = c(country_levels, "Population")))

  p <- ggplot(pd, aes(x = mean, y = ylab, colour = enso_group, shape = tier)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, linewidth = 0.35) +
    geom_point(size = 1.1) +
    facet_grid(enso_group ~ term, scales = "free", space = "free_y",
               labeller = ggplot2::labeller(enso_group = GROUP_SHORT_FIG)) +
    scale_x_log10() +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_shape_manual(values = c(country = 16, population = 18), name = NULL,
                       labels = c(country = "Country", population = "Population (pooled)")) +
    labs(x = "Posterior ratio (exp scale), 95% CrI", y = NULL,
         title = sprintf("Region/country parameters (%s)", TRACK_LABEL[[track]]),
         subtitle = "population + region + country REs; ENSO exposure-response is Figures 1-2") +
    theme_ms() +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = BASE_SIZE - 2),
          strip.text.y = element_text(size = BASE_SIZE - 1),
          axis.text.y  = element_text(size = BASE_SIZE - 2),
          axis.text.x  = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2))

  n_rows  <- nrow(dplyr::distinct(pd, enso_group, ylab))
  n_terms <- nlevels(pd$term)
  stem <- sprintf("figS_params_%s", track)
  save_figure(p, stem,
              height_mm = 32 + n_rows * 3.4,
              width_mm  = 45 + n_terms * 34, subdir = "supplementary")
}

#' All requested analog sets' ENSO forecast input on one panel (vs. the
#' single-set ONI strip above Figure 3b). Only meaningful when >1 set is on hand.
build_fig_oni_comparison <- function() {
  message("[10] Supplementary: ONI envelope comparison (all analog sets)")
  oni <- read_results("oni", sets = SETS)
  key_levels <- c("Forecast", set_lab(SETS))
  env <- oni$env %>%
    dplyr::mutate(key = factor(ifelse(arm == "forecast", "Forecast", set_lab(analog_set)),
                               levels = key_levels))
  pal <- c(Forecast = FC_COLOR, setNames(unname(SET_COLORS[SETS]), set_lab(SETS)))
  fig <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(data = env, aes(date, ymin = lo, ymax = hi, fill = key), alpha = 0.18) +
    geom_line(data = env, aes(date, mean, colour = key), linewidth = 0.6) +
    geom_line(data = oni$obs, aes(date, enso_index), colour = "black", linewidth = 0.45) +
    # Both arms' launches when they differ (oni v4 `handoffs`), else one line.
    geom_vline(xintercept = unique(as.Date(c(oni$handoff, unlist(oni$handoffs)), origin = "1970-01-01")),
               linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(x = NULL, y = "Oceanic Ni\u00f1o Index (ONI)",
         title = "ENSO forecast input by counterfactual analog set") +
    theme_ms()
  save_figure(fig, "fig3_supp_oni_envelope_comparison", height_mm = 80, width_mm = 200, subdir = "supplementary")
}

# ===========================================================================
# Section 9f: Figure H — hindcast forecast skill
# ===========================================================================
#' Reads the four phase-09 hindcast artefacts; every panel is a shared builder
#' from utils/figures_hindcast.R (the notebook draws the same plots).
build_fig_hindcast <- function() {
  message("[10] Figure H: hindcast forecast skill")
  runs <- read_results("hindcast_runs"); occ <- read_results("hindcast_occ")
  cas  <- read_results("hindcast_cases"); rep <- read_results("report_hindcast")
  hc <- SCOPE$hindcast
  groups <- intersect(GROUP_ORDER, unique(runs$group))
  lab <- hc_scope_labels(GROUP_SHORT_FIG, groups)
  pal <- hc_scope_palette(GROUP_COLOURS, groups)
  th  <- theme_ms() + ggplot2::theme(legend.position = "bottom")
  refs <- intersect(c(hc$headline_occ_reference, "insample_gid", "regional", "persistence"),
                    occ$meta$references)

  pa <- hc_plot_skill_bars(occ$skill, refs, hc$headline_occ_reference, "window", "crps_fair",
                           c("all", rep$headline$short_lead), lab, th) +
    ggplot2::labs(title = "Outbreak-count skill vs climatology references") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
  pb <- hc_plot_skill_by_lead(occ$skill, hc$headline_occ_reference, "window", "crps_fair",
                              occ$meta$lead_levels, lab, pal, th) +
    ggplot2::labs(title = "Skill by lead") +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(.98, .98),
                   legend.justification = c(1, 1))
  pc <- hc_plot_case_support(cas$calibration, cas$skill,
                             intersect(c(hc$headline_case_support, "obs_positive_keys", "footprint"),
                                       cas$meta$supports),
                             hc$headline_case_support, labels = lab, theme = th) +
    ggplot2::labs(title = "Case totals by support")
  pd <- hc_plot_reliability(occ$reliability, lab, th) + ggplot2::labs(title = "Reliability")
  fig <- (pa | pb) / (pc | pd) +
    patchwork::plot_layout(heights = c(1, 1.4)) +
    patchwork::plot_annotation(tag_levels = "a") & TAG_THEME
  save_figure(fig, "figH_hindcast_skill", height_mm = 230, subdir = "hindcast")

  if (SUPP) {
    supp <- function(p, stem, h) save_figure(p, stem, height_mm = h, subdir = "hindcast")
    supp(hc_plot_skill_bars(occ$skill, occ$meta$references, hc$headline_occ_reference,
                            "window", "crps_fair", occ$meta$lead_levels, lab, th),
         "figH_supp_all_references", 110)
    supp(hc_plot_skill_bars(occ$skill, c("insample_gid", hc$headline_occ_reference),
                            hc$headline_occ_reference, "calendar_year", "crps_emp", "all",
                            lab, th) + ggplot2::labs(title = "Legacy basis (calendar-year target)"),
         "figH_supp_legacy_basis", 70)
    supp((hc_plot_pit_hist(occ$pit_hist, lab, th) | hc_plot_rank_hist(occ$rank_hist, "window",
                                                                      "all", lab, th)) /
           (hc_plot_coverage_curve(occ$coverage_curve, lab, pal, th) |
              hc_plot_roc(occ$roc, dplyr::filter(occ$cell_scores, lead_class == "all",
                                                 reference == hc$headline_occ_reference),
                          lab, pal, th)) +
           patchwork::plot_annotation(tag_levels = "a") & TAG_THEME,
         "figH_supp_calibration", 200)
    supp(hc_plot_country_lollipop(cas$skill, hc$headline_case_support, theme = th),
         "figH_supp_country_cases", 120)
    supp(hc_plot_cases_interval(cas$run_totals, hc$headline_case_support, "region", "all", lab, th),
         "figH_supp_case_intervals", 120)
  }
}

# ===========================================================================
# Section 10: main dispatch
# ===========================================================================
message(sprintf("[10] figures=%s supplementary=%s main_set=%s supp_sets=%s",
                paste(FIGURES, collapse = ","), SUPP, MAIN_SET,
                paste(SUPP_SETS, collapse = ",")))

if ("1" %in% FIGURES) {
  build_fig1(STACKED_TRACKS)
  if (SUPP) for (trk in FIG1_TRACKS) build_fig1_supp_country_erf(trk)
}

if ("2" %in% FIGURES) {
  build_fig2(STACKED_TRACKS)
}

if ("3" %in% FIGURES) {
  build_fig3("cases", MAIN_SET)
  build_fig3("onsets", MAIN_SET)
  if (SUPP) {
    for (aset in SUPP_SETS) {
      for (metric in c("cases", "onsets")) {
        build_fig3(metric, aset)
      }
    }
    if (length(SETS) > 1L) build_fig_oni_comparison()
  }
}

if ("4" %in% FIGURES) {
  # Every Figure 4 shows both metrics (cases | outbreak onsets).
  build_fig4(MAIN_SET)
  if (SUPP) for (aset in SUPP_SETS) build_fig4(aset)
}

if ("H" %in% FIGURES) {
  if (isTRUE(SCOPE$hindcast$enable)) build_fig_hindcast()
  else message("[10] figure H requested but results.scope.hindcast.enable is false; skipped")
}

# Supplementary region/country parameter caterpillars (every non-ENSO term).
if (SUPP) for (trk in FIG1_TRACKS) build_fig_supp_track_params(trk)

# Fig. S5: sampler diagnostics of the deployed models (convergence artefacts).
if ("S5" %in% FIGURES) {
  rc   <- read_results("report_convergence")
  summ <- dplyr::bind_rows(lapply(FIG1_TRACKS, function(t) read_results("convergence", track = t)$summary))
  save_figure(build_figS5_convergence(rc$fits, summ, rc$thresholds, GROUP_COLOURS,
                                      GROUP_SHORT_FIG, theme_ms()),
              "figS5_convergence", height_mm = 150, subdir = "supplementary")
}

# Fig. S6: LOO lag rankings (report_lags; the ranking Figure 1c outlines).
if ("S6" %in% FIGURES) {
  rl <- read_results("report_lags")
  save_figure(build_figS6_lag_loo(rl$ranking, rl$selected, GROUP_COLOURS, GROUP_SHORT_FIG,
                                  theme_ms()),
              "figS6_lag_loo", height_mm = 120, subdir = "supplementary")
}

# Forecast vs counterfactual by country (attribution$country_ribbon; the
# 07d_simulator_figures_counterfactual.R fig6/fig7 series), main analog set.
if ("FC" %in% FIGURES) {
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
if ("SC" %in% FIGURES) {
  sm <- read_results("smoothing")
  n_ct <- length(sm$country_group)
  for (m in c("cases", "onsets"))
    save_figure(build_figS_smoothing_country(sm$by_country, sm$country_group, m, GROUP_ORDER,
                                             GROUP_SHORT_FIG, F3_MASK_COL, theme_ms(),
                                             start = as.Date(RES$derivation$smoothing_start)),
                sprintf("figS_smoothing_country_%s", m), height_mm = 22 + 32 * ceiling(n_ct / 5),
                subdir = "supplementary")
}

# Population in ADM1s with more ENSO-attributable cases (report_attrib_exposure):
# ratio-of-means map with the IPCC robustness hatching at the headline m
# (style.attrib_headline) | population shares by class, every m. Main analog set.
if ("AE" %in% FIGURES) {
  ae <- read_results("report_attrib_exposure", sets = SETS)
  hl_m <- as.numeric(RES$style$attrib_headline$magnitude)
  if (length(hl_m) != 1L || !any(abs(ae$magnitudes - hl_m) < 1e-9))
    stop("[10] figure AE: results.style.attrib_headline$magnitude must be one of thresholds.attrib_magnitudes {",
         paste(ae$magnitudes, collapse = ", "), "}.", call. = FALSE)
  maps <- load_africa_maps("Figure AE")
  basemap <- load_africa_basemap(maps$afro_sf)
  hatch <- if (pattern_fill_supported())
    list(none = mask_hatch(spacing_mm = 0.8, lwd = 0.3, col = "grey35"),
         conflicting = mask_hatch(spacing_mm = 0.8, lwd = 0.3, col = "grey35", cross = TRUE))
  if (is.null(hatch)) message("[10] figure AE: ggplot2 < 3.5, no pattern fills -- robust ADM1s outlined instead")
  grey_sf <- ps_country_status(maps$afro_sf, modelled_countries())
  for (m in c("cases", "onsets")) {
    # A fresh grey layer per figure: ggnewscale renames a layer's fill in place
    # on first use, so a reused new_scale_fill() layer breaks the second plot.
    grey <- ps_country_grey_layer(grey_sf, legend_order = 3, labels = PS_GREY_LAB_WRAP)
    # No scale bar: the map's legends take its bottom-left corner (as Figure 4).
    fig <- build_figS_attrib_exposure(
             ae$adm1, ae$exposure, m, MAIN_SET, hl_m, maps$adm1_sf, sf::st_bbox(maps$afro_sf),
             base_layers = grey, overlay_layers = map_overlay_layers(basemap, maps$afro_sf),
             limits = c(FIG2_MAP_LO, FIG2_MAP_HI), scope_order = c("Africa", unname(GROUP_SHORT[GROUP_ORDER])),
             class_cols = c(robust = "#b2182b", conflicting = "#f4a582", childs = "black"),
             hatch_fills = hatch, map_theme = taxdat::map_theme() + ps_map_font(),
             theme = theme_ms(), base_size = BASE_SIZE) +
      patchwork::plot_annotation(tag_levels = "a") & TAG_THEME
    save_figure(fig, sprintf("figS_attrib_exposure_%s", m), height_mm = 88, subdir = "supplementary")
  }
}

# Attributable fraction of the past strong ENSO events (report_event_attrib +
# report_event_af_model): simulated (both factual readings) vs analytic, per
# event and ENSO group.
if ("EV" %in% FIGURES) {
  if (!nzchar(RES$event_manifest %||% ""))
    stop("[10] figure EV: results.event_manifest is empty (no event attribution derived).", call. = FALSE)
  ea  <- read_results("report_event_attrib")
  eam <- read_results("report_event_af_model")
  for (m in c("cases", "onsets"))
    save_figure(build_figS_event_attrib(ea, eam, m, c("Africa", unname(GROUP_SHORT[GROUP_ORDER])),
                                        theme_ms(), base_size = BASE_SIZE,
                                        tag = c(cases = "A", onsets = "B")[[m]]) + TAG_THEME,
                sprintf("figS_event_attrib_%s", m), height_mm = 75, subdir = "supplementary")
}

message("[10] done.")
