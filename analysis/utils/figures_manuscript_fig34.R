# analysis/utils/figures_manuscript_fig34.R
# =============================================================================
# Figure 3 (past strong ENSO events) and Figure 4 (the 2026 forecast vs the
# neutral-ENSO counterfactual) builders of 10_manuscript_figures.R, with the
# forecast panels and ONI strip they share (fig3_panel_a/b/c(), fig3_oni_strip(),
# also Figure 1a's strip). The F3_* / F3_EV_* / FIG4_* constants live in 10.
#
# Definitions only (loads alone). These builders are 10's own: at CALL time they
# read 10's globals -- style (FONT, BASE_SIZE, theme_ms(), TAG_THEME, PS_*),
# groups (GROUP_*), config (RES, SCOPE, SETS, MAIN_SET, opt), the figure
# constants and 10's helpers (read_results(), load_africa_maps(), save_figure(),
# ...). Only 10_manuscript_figures.R sources them, in its module loop.
# =============================================================================

f3_scope_levels <- function(groups) c("Africa", unname(GROUP_SHORT_FIG[groups]))

f3_scope_pal    <- function(groups) c(Africa = F3_AFRICA_COL,
                                      stats::setNames(unname(GROUP_COLOURS[groups]), GROUP_SHORT_FIG[groups]))

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
#' `fc_df`: the forecast ribbon (NULL = no forecast series or legend entry).
#' `ymax`: y cap (NULL = none).
#' `mask_col`: the coverage-masked smoothing's colour (default F3_MASK_COL).
fig3_panel_a <- function(sm_afr, fc_df, metric, x_limits, vlines, ymax = NULL, phases = NULL,
                         legend_pos = "bottom", legend_nrow = 1L, mask_col = F3_MASK_COL) {
  lv <- c("Smoothing (coverage-masked)", "Smoothing (unmasked)", "Observed", if (!is.null(fc_df)) "Forecast")
  k <- seq_along(lv)
  pal <- c(mask_col, F3_AFRICA_COL, "black", F3_AFRICA_COL)[k]; names(pal) <- lv
  lty <- c("solid", "solid", "blank", F3_ARM_LTY[["forecast"]])[k]; names(lty) <- lv
  sm_q <- function(r) r %>% dplyr::transmute(t = t_start, q10 = q025, q25, q50, q75, q90 = q975)
  long <- dplyr::bind_rows(
    sm_q(sm_afr$ribbon_masked)   %>% dplyr::mutate(series = lv[1]),
    sm_q(sm_afr$ribbon_unmasked) %>% dplyr::mutate(series = lv[2]),
    if (!is.null(fc_df)) f3_qcols(fc_df, metric) %>% dplyr::mutate(series = lv[4])) %>%
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
      linetype = unname(lty), shape = c(NA, NA, 21, NA)[k],
      fill = c(NA, NA, "white", NA)[k], size = c(NA, NA, 1.2, NA)[k], stroke = 0.4), nrow = legend_nrow,
      byrow = TRUE)) +
    scale_x_date(limits = x_limits, date_breaks = "2 years", date_labels = "%Y", expand = c(0, 0)) +
    scale_y_continuous(labels = PS_LABEL_K) +
    coord_cartesian(ylim = if (!is.null(ymax)) c(0, ymax)) +
    labs(x = NULL, y = METRIC_LABELS[metric]) +
    theme_ms() +
    theme(legend.position = legend_pos, legend.justification = "left",
          legend.key.width = unit(6, "mm"), legend.key.height = unit(3, "mm"))
}

#' The continuation months of panel b (09's attribution$africa_window): every
#' group read past its firm horizon on the simulated analog ENSO
#' continuation gets its months after the firm horizon veiled (FIG4_CONT_VEIL),
#' a dotted line at the firm horizon and the vertical label "continuation"
#' along it (the leads are in table S28); with a common
#' Africa window the Africa row gets a dotted line at the window end. Empty
#' when no group is read past its firm horizon (legacy windows).
f3_continuation_layers <- function(f23, groups, lv) {
  aw <- f23$africa_window
  if (is.null(aw) || !length(aw$lead)) return(NULL)
  gl <- names(aw$lead)[aw$lead > 0 & names(aw$lead) %in% groups]
  if (!length(gl)) return(NULL)
  half <- 15   # days: the veil starts between the firm horizon and the next month
  cdf <- data.frame(scope = factor(unname(GROUP_SHORT_FIG[gl]), levels = lv),
                    xmin = as.Date(aw$firm[gl]) + half, xmax = as.Date(aw$group_ends[gl]) + half,
                    lab = "continuation")
  txt <- (BASE_SIZE - 2.5) / ggplot2::.pt
  out <- list(
    geom_rect(data = cdf, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
              fill = "white", alpha = FIG4_CONT_VEIL),
    geom_vline(data = cdf, aes(xintercept = xmin), inherit.aes = FALSE, linetype = "dotted", colour = "grey45",
               linewidth = 0.3),
    # Vertical, along the dotted line (on its right), from the panel top.
    geom_text(data = cdf, aes(x = xmin, y = Inf, label = lab), inherit.aes = FALSE, angle = 90, hjust = 1.08,
              vjust = -0.4, size = txt, colour = "grey30"))
  if (identical(aw$mode, "common") && "Africa" %in% lv)
    out <- c(out, list(geom_vline(data = data.frame(scope = factor("Africa", levels = lv),
                                                    x = as.Date(aw$end) + half),
                                  aes(xintercept = x), inherit.aes = FALSE, linetype = "dotted", colour = "grey45",
                                  linewidth = 0.3)))
  out
}

#' Panel b: Africa + regions from the analogue handoff, forecast (scope colour,
#' dot-dashed) vs counterfactual (F3_CF_COL, dashed), median + 50% + 95% bands.
#' Africa on its own y scale, the regions sharing one (free_y + a blank point
#' at their common maximum); `africa_ymax` / `region_ymax` cap them (bands
#' squished; NULL = none).
#' When the factual series is stitched, [attrib_start, factual_start) is shaded
#' and labelled "informed by national data" (vertical, left of the factual
#' handoff's dashed vline).
fig3_panel_b <- function(f23, metric, aset, groups, x_limits, africa_ymax = NULL, region_ymax = NULL,
                         attrib_series = FALSE) {
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
    # Vertical, left of the dashed handoff line, from the panel top.
    annotate("text", x = as.Date(f23$factual_start), y = Inf, label = "informed by\nnational data",
             angle = 90, hjust = 1.06, vjust = -0.3, lineheight = 0.9, size = (BASE_SIZE - 2.5) / ggplot2::.pt,
             colour = "grey35"),
    geom_vline(xintercept = as.Date(f23$factual_start), linetype = "dashed", colour = "grey45",
               linewidth = 0.3)) else NULL
  # The attribution factual's median (Figure 4: the factual of the paired
  # statistics when it is not the displayed forecast), dotted, scope colour.
  attrib_layer <- NULL
  if (isTRUE(attrib_series) && !is.null(f23$attrib_ribbon)) {
    da <- dplyr::bind_rows(
      f23$attrib_ribbon %>% dplyr::filter(enso_group %in% groups) %>%
        dplyr::mutate(scope = unname(GROUP_SHORT_FIG[enso_group])),
      f23$attrib_africa_ribbon %>% dplyr::mutate(scope = "Africa")) %>%
      dplyr::filter(month >= f23$attrib_start) %>%
      dplyr::mutate(scope = factor(scope, levels = lv), key = as.character(scope))
    da <- dplyr::bind_cols(da, f3_qcols(da, metric) %>% dplyr::select(-t))
    if (!is.null(region_ymax)) da$q50[da$scope != "Africa"] <- pmin(da$q50[da$scope != "Africa"], region_ymax)
    if (!is.null(africa_ymax)) da$q50[da$scope == "Africa"] <- pmin(da$q50[da$scope == "Africa"], africa_ymax)
    attrib_layer <- geom_line(data = da, aes(month, q50, colour = key), linetype = F4_ATTRIB_LTY,
                              linewidth = F3_LWD$b, inherit.aes = FALSE)
  }
  ggplot(d, aes(month)) +
    stitch_layers +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = key), alpha = 0.18) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = key), alpha = 0.3) +
    geom_line(aes(y = q50, colour = key, linetype = arm), linewidth = F3_LWD$b) +
    attrib_layer +
    f3_continuation_layers(f23, groups, lv) +
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

#' ONI strip sharing its panel's x axis (Figures 1a, 3a, 4a). mode =
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

#' The events of report_event_attrib with display labels ("El Nino 2015-16" ->
#' "El Niño 2015–16"), bracket tiers (greedy: the lowest tier whose last
#' window ends before the launch) and the ONI at the peak month.
f3ev_events <- function(ea, oni) {
  ev <- ea$events %>%
    dplyr::arrange(launch) %>%
    dplyr::mutate(label = gsub("Nina", "Niña", gsub("Nino", "Niño", label)),
                  label = sub("(\\d{4})-(\\d{2})$", "\\1–\\2", label),
                  launch = as.Date(launch), peak = as.Date(peak), window_end = as.Date(window_end))
  tier_end <- as.Date(character(0)); ev$tier <- NA_integer_
  for (i in seq_len(nrow(ev))) {
    free <- which(tier_end < ev$launch[i])
    k <- if (length(free)) free[1] else length(tier_end) + 1L
    tier_end[k] <- ev$window_end[i]; ev$tier[i] <- k
  }
  ev$peak_oni <- oni$obs_full$enso_index[match(ev$peak, oni$obs_full$date)]
  # A multi-winter event (report_event_attrib sub-windows):
  # its earlier winters' peaks, drawn as dots on its bracket (an attribute, so
  # an artefact without sub-windows leaves the strip unchanged).
  if (!is.null(ea$sub_windows)) {
    pk <- dplyr::distinct(ea$sub_windows$windows, event, peak) %>%
      dplyr::mutate(peak = as.Date(peak)) %>%
      dplyr::inner_join(dplyr::select(ev, event, tier, phase, last_peak = peak), by = "event") %>%
      dplyr::filter(peak != last_peak)
    if (nrow(pk)) attr(ev, "extra_peaks") <- pk
  }
  ev
}

#' Event strip above a's ONI: one bracket per event from launch to window end
#' (end ticks pointing down to the ONI), a dot at the ONI peak (one per winter
#' of a multi-winter event: attribute "extra_peaks" of f3ev_events()), the label over
#' the bracket's left end; tier 1 lowest. Same x scale as panel a.
f3ev_event_strip <- function(ev, x_limits) {
  n_tier <- max(ev$tier)
  tick <- 0.28
  ggplot(ev) +
    geom_segment(aes(x = launch, xend = window_end, y = tier, yend = tier, colour = phase), linewidth = 0.35) +
    geom_segment(aes(x = launch, xend = launch, y = tier, yend = tier - tick, colour = phase), linewidth = 0.35) +
    geom_segment(aes(x = window_end, xend = window_end, y = tier, yend = tier - tick, colour = phase),
                 linewidth = 0.35) +
    geom_point(aes(x = peak, y = tier, colour = phase), size = 0.9) +
    (if (!is.null(attr(ev, "extra_peaks")))
      geom_point(data = attr(ev, "extra_peaks"), aes(x = peak, y = tier, colour = phase), size = 0.9)) +
    geom_text(aes(x = launch, y = tier + 0.24, label = label, colour = phase), hjust = 0, vjust = 0,
              size = (BASE_SIZE - 1.5) / ggplot2::.pt, family = FONT) +
    scale_colour_manual(values = F3_EV_PHASE_COL, guide = "none") +
    scale_x_date(limits = x_limits, date_breaks = "2 years", expand = c(0, 0)) +
    scale_y_continuous(limits = c(1 - tick - 0.05, n_tier + 0.75), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_ms() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "mm"),
          axis.line = element_blank(), panel.border = element_blank(),
          # theme_ms() sets the major / minor grids explicitly, so blank them by name.
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(1, 1.5, 0, 1.5, unit = "mm"))
}

#' The ENSO-state legend (Section 2b's guide), compact, for a's legend row.
f3ev_phase_legend <- function() {
  d <- data.frame(ph = factor(names(PS_ENSO_PHASE_COL), levels = names(PS_ENSO_PHASE_COL)), x = 1, y = 1)
  cowplot::get_legend(
    ggplot(d, aes(x, y, fill = ph)) +
      geom_tile(alpha = PS_ENSO_PHASE_ALPHA, colour = PS_ENSO_KEY_OUTLINE$colour,
                linewidth = PS_ENSO_KEY_OUTLINE$linewidth) +
      scale_fill_manual(values = PS_ENSO_PHASE_COL, labels = PS_ENSO_PHASE_LAB, guide = ps_enso_phase_guide()) +
      theme_ms() + F3_LEGEND_THEME + ps_legend_theme())
}

#' Legend grobs side by side at their natural widths, F3_EV_LEG_GAP_MM apart,
#' centred on the plot area (WIDTH_MM minus F3_EV_GUTTER_MM): the row sits under
#' the panels instead of spanning the page from its left edge.
f3ev_legend_row <- function(legs) {
  w <- vapply(legs, function(g) grid::convertWidth(sum(g$widths), "mm", valueOnly = TRUE), numeric(1))
  free <- max(0, WIDTH_MM - F3_EV_GUTTER_MM - sum(w) - F3_EV_LEG_GAP_MM * (length(w) - 1L))
  cells <- list(NULL); rw <- F3_EV_GUTTER_MM + free / 2
  for (i in seq_along(legs)) {
    cells <- c(cells, list(legs[[i]])); rw <- c(rw, w[i])
    if (i < length(legs)) { cells <- c(cells, list(NULL)); rw <- c(rw, F3_EV_LEG_GAP_MM) }
  }
  cowplot::plot_grid(plotlist = c(cells, list(NULL)), nrow = 1, rel_widths = c(rw, free / 2))
}

#' Panel a's cells (event strip, ONI strip, Africa smoothing, legend row), tagged
#' "a": Africa's smoothing (masked / unmasked) and the observed series to the
#' forecast handoff, on the ENSO-phase background. a's axis titles and tick
#' labels take no column of their own (free "space"): they sit in the gutter of
#' b's scope labels, so the page keeps no empty title column.
f3ev_top_panels <- function(oni, sm, f23, ev, metric) {
  handoff_month <- as.Date(f23$attrib_start)
  x_limits_a <- c(as.Date("2011-01-01"), handoff_month)
  ph_a <- enso_phase_rects(oni, x_limits_a)
  sm_afr <- sm$by_scope_metric[[paste("Africa", metric, sep = "|")]]
  cap_a <- if (identical(metric, "cases")) F3_CASES_YMAX[["a"]] else NULL
  pa <- fig3_panel_a(sm_afr, NULL, metric, x_limits_a, handoff_month, ymax = cap_a, phases = ph_a,
                     mask_col = F3_EV_MASK_COL)
  leg_row_a <- f3ev_legend_row(list(
    cowplot::get_legend(pa + F3_LEGEND_THEME + ps_legend_theme() +
                          theme(legend.key.spacing.x = unit(PS_KEY_SPACING_MM, "mm"))),
    f3ev_phase_legend()))
  lm <- function(p, t, r, b) p + theme(plot.margin = margin(t, r, b, F3_EV_LEFT_MARGIN_MM, unit = "mm"))
  strip <- lm(f3ev_event_strip(ev, x_limits_a) + labs(tag = "a"), 1, 1.5, 0)
  oni_a <- patchwork::free(lm(fig3_oni_strip(oni, MAIN_SET, x_limits_a, mode = "observed", phases = ph_a),
                              1.5, 1.5, 0) +
                             theme(axis.title.y = element_text(margin = margin(r = 0.5, unit = "mm"))),
                           type = "space", side = "l")
  pa <- patchwork::free(lm(pa + theme(legend.position = "none"), 1, 1.5, 0.5) +
                          theme(axis.title.y = element_text(margin = margin(r = 1.5, unit = "mm"))),
                        type = "space", side = "l")
  list(strip, oni_a, pa, patchwork::wrap_elements(full = leg_row_a, ignore_tag = TRUE))
}

#' Panel b's table for one metric: Africa and the ENSO groups per event, the
#' bar `af_bar` = the median of the noise-adjusted posterior AF, with its 50% /
#' 95% interval (`post`: 09's event_af_posterior), in % (as fractions here,
#' scaled below); an unidentified scope takes its pooled value and raw range,
#' flagged `unidentified` (ps_af_posterior_cols()).
f3ev_scope_table <- function(ea, post, metric, events) {
  p <- post$scope %>%
    dplyr::filter(explosive_guard == F3_EV_GUARD, metric == !!metric) %>%
    ps_af_posterior_cols(scale = 1) %>%
    dplyr::select(event, level, key, af_lo, af_q25, af_bar, af_q75, af_hi, unidentified)
  ea$af_scope %>%
    dplyr::filter(factual_source == F3_EV_FACTUAL, explosive_guard == F3_EV_GUARD, metric == !!metric,
                  level %in% c("continent", "region")) %>%
    dplyr::mutate(key = ifelse(level == "continent", "Africa", enso_group)) %>%
    # the posterior interval replaces the artefact's bootstrap one
    dplyr::select(-dplyr::any_of(c("af_lo", "af_hi", "af_q025", "af_med", "af_q975"))) %>%
    dplyr::left_join(p, by = c("event", "level", "key")) %>%
    dplyr::mutate(scope = ifelse(level == "continent", "Africa", unname(GROUP_SHORT_FIG[enso_group])),
                  scope = factor(scope, levels = f3_scope_levels(GROUP_ORDER)),
                  fill_key = as.character(scope),
                  event = factor(event, levels = events),
                  dplyr::across(c(af_bar, af_lo, af_hi, af_q25, af_q75), ~ 100 * .x)) %>%
    dplyr::filter(!is.na(scope))
}

#' One column of panel b (one event): horizontal bars, Africa on top, the 50% /
#' 95% intervals (Section 2b), arrows where the 95% interval runs past `xlim`,
#' the event label as the column title in its phase colour. `show_y`: scope
#' labels (first column only); `x_title`: the AF axis title (middle column).
f3ev_panel_b <- function(d_ev, ev_row, xlim, show_y = FALSE, x_title = NULL) {
  d_ev$col <- ev_row$label
  d_ev$scope <- droplevels(d_ev$scope)
  rows <- ps_africa_rows(d_ev$scope)
  d_ev$ypos <- rows$y
  d_ev$af_bar <- pmax(pmin(d_ev$af_bar, xlim[2]), xlim[1])   # bar + point clamped, the arrow marks the cut
  ggplot(d_ev, aes(y = ypos, x = af_bar, fill = fill_key)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25) +
    rows$rule +
    geom_col(width = 0.7, orientation = "y") +
    ps_af_interval_layers("af_bar", "af_lo", "af_hi", "af_q25", "af_q75", orientation = "y",
                          pos = "ypos", clip = xlim) +
    facet_wrap(~ col) +
    scale_fill_manual(values = f3_scope_pal(GROUP_ORDER), guide = "none") +
    scale_y_continuous(breaks = rows$breaks, labels = rows$labels, expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = c(-50, 0, 50, 100), labels = true_minus(c("-50", "0", "50", "100"))) +
    coord_cartesian(xlim = xlim) +
    labs(y = NULL, x = x_title) +
    theme_ms() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = BASE_SIZE - 0.5,
                                    colour = F3_EV_PHASE_COL[[ev_row$phase]]),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = BASE_SIZE - 1.5),
          axis.text.y = element_text(size = BASE_SIZE),
          plot.margin = margin(0.5, 2.2, 0.5, 2.2, unit = "mm")) +
    { if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                         axis.ticks.length.y = unit(0, "mm")) }
}

#' Figure 3c's crop for one event: the bbox of the ADM1s with an AF estimate
#' (`af` not missing) in `d_ev`, padded by `pad` degrees on every side.
f3ev_bbox <- function(d_ev, adm1_sf, pad = F3_EV_MAP_PAD_DEG) {
  units <- unique(d_ev$adm1_unit[!is.na(d_ev$af)])
  if (!length(units)) stop("[10] Figure 3c: no ADM1 estimate to crop the map to.", call. = FALSE)
  bb <- sf::st_bbox(adm1_sf[adm1_sf$GID_1 %in% units, ])
  bb[c("xmin", "ymin")] <- bb[c("xmin", "ymin")] - pad
  bb[c("xmax", "ymax")] <- bb[c("xmax", "ymax")] + pad
  bb
}

#' One map of panel c (one event): AF (%) per ADM1 on the diverging scale
#' (squished at F3_EV_AF_LIM), the two map greys, rivers / lakes / borders;
#' cropped to `bb`. `show_legend`: the legend row (colourbar + grey keys);
#' `scale_bar`: scale bar + north arrow (the first map).
f3ev_map <- function(d_ev, maps, basemap, obc, bb, show_legend = FALSE, scale_bar = FALSE) {
  map_sf <- maps$adm1_sf %>% dplyr::inner_join(d_ev, by = c("GID_1" = "adm1_unit"))
  ggplot() +
    # A fresh grey layer per map: ggnewscale renames a layer's fill in place.
    ps_country_grey_layer(ps_country_status(maps$afro_sf, obc), legend_order = 2) +
    geom_sf(data = map_sf, aes(fill = af), colour = "white", linewidth = 0.005) +
    map_overlay_layers(basemap, maps$afro_sf) +
    { if (scale_bar) map_scale_layers() } +
    # AF cannot exceed 100%, so only the low end is marked as clipped.
    ps_diverging_fill_scale(limits = F3_EV_AF_LIM, breaks = c(-100, -50, 0, 50, 100),
                            labels = true_minus(c("≤ -100", "-50", "0", "50", "100")),
                            name = "ENSO-attributable fraction, ADM1 (%)") +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], expand = FALSE) +
    taxdat::map_theme() +
    ps_map_font() +
    theme(panel.border = element_blank(), plot.margin = margin(0.5, 1, 0.5, 1, unit = "mm"),
          text = element_text(size = BASE_SIZE - 1),
          legend.position = if (show_legend) "bottom" else "none",
          legend.box = "horizontal", legend.box.just = "top", legend.justification = "center",
          legend.direction = "horizontal") +
    ps_legend_theme() +
    theme(legend.spacing.x = unit(6, "mm"),
          legend.title = element_text(size = BASE_SIZE - 1, vjust = 0.8),
          legend.key.height = unit(2.4, "mm"), legend.key.width = unit(3, "mm"),
          legend.key.spacing.x = unit(PS_KEY_SPACING_MM, "mm"))
}

#' Inset the group map `inset` (fig2_group_map_inset()) into plot `p`, in the
#' box `box` = c(left, bottom, right, top) in npc of p's panel; NULL `box`
#' leaves p as it is.
f34_add_group_map <- function(p, inset, box) {
  if (is.null(box)) return(p)
  p + patchwork::inset_element(inset, left = box[["left"]], bottom = box[["bottom"]], right = box[["right"]],
                               top = box[["top"]], align_to = "panel", ignore_tag = TRUE)
}

#' Figure 3: attribution of the past strong ENSO events (report_event_attrib +
#' event_af_posterior; results.event_manifest), one figure per metric:
#'   a  event brackets (launch to window end, dot at the ONI peak) over the
#'      observed-ONI strip over Africa's smoothing (masked / unmasked) and the
#'      observed series, on the ENSO-phase background;
#'   b  one column per event: the ENSO-attributable fraction AF = 1 - CF / F
#'      for Africa and each ENSO group (factual = free simulation vs the
#'      neutral-ENSO counterfactual, explosive replicates dropped): median and
#'      50% / 95% of the noise-adjusted posterior;
#'   c  AF of the replicate means per ADM1 under each event's bars.
#' The group map inset in b's F3_EV_INSET_COL-th column (F3_EV_INSET_BOX).
#' Main: cases; supplementary: onsets.
build_fig3 <- function(metrics = c("cases", "onsets")) {
  if (!nzchar(RES$event_manifest %||% ""))
    stop("[10] Figure 3: results.event_manifest is empty (no event attribution derived).", call. = FALSE)
  oni  <- read_results("oni", sets = SETS)
  sm   <- read_results("smoothing")
  f23  <- read_results("attribution", sets = SETS)   # panel a ends at its attribution start
  ea   <- read_results("report_event_attrib")
  post <- read_results("event_af_posterior")
  obc  <- modelled_countries()
  maps <- load_africa_maps("Figure 3")
  basemap <- africa_basemap(maps$afro_sf)
  ev <- f3ev_events(ea, oni)
  n_ev <- nrow(ev); mid <- ceiling(n_ev / 2)
  if (!is.null(F3_EV_INSET_BOX) && !F3_EV_INSET_COL %in% seq_len(n_ev))
    stop("[10] Figure 3: F3_EV_INSET_COL ", F3_EV_INSET_COL, " is not an event column (1..", n_ev, ").", call. = FALSE)
  inset <- if (!is.null(F3_EV_INSET_BOX))
    fig2_group_map_inset(maps, read_results("modelled_units")$members)   # the fit registry's partition (09)
  for (metric in metrics) {
    message(sprintf("[10] Figure 3: metric=%s", metric))
    top <- f3ev_top_panels(oni, sm, f23, ev, metric)
    d_b <- f3ev_scope_table(ea, post, metric, ev$event)
    d_c <- ea$af_adm1 %>%
      dplyr::filter(factual_source == F3_EV_FACTUAL, explosive_guard == F3_EV_GUARD, metric == !!metric) %>%
      dplyr::mutate(af = 100 * af_of_means)
    # Each map cropped to the ADM1s carrying an estimate for its event,
    # F3_EV_MAP_PAD_DEG around (the footprints coincide, so the maps align).
    bbs <- lapply(ev$event, function(e) f3ev_bbox(dplyr::filter(d_c, event == e), maps$adm1_sf))
    pb <- lapply(seq_len(n_ev), function(i) {
      p <- f3ev_panel_b(dplyr::filter(d_b, event == ev$event[i]), ev[i, ], F3_EV_B_XLIM, show_y = i == 1L,
                        x_title = if (i == mid) sprintf("ENSO-attributable fraction of %s (%%)", metric))
      if (i == 1L) p <- p + labs(tag = "b") + theme(plot.margin = margin(0.5, 2.2, 0.5, F3_EV_LEFT_MARGIN_MM, unit = "mm"))
      if (i == F3_EV_INSET_COL) f34_add_group_map(p, inset, F3_EV_INSET_BOX) else p
    })
    pc <- lapply(seq_len(n_ev), function(i) {
      p <- f3ev_map(dplyr::filter(d_c, event == ev$event[i]), maps, basemap, obc, bbs[[i]], scale_bar = i == 1L)
      if (i == 1L) p + labs(tag = "c") + theme(plot.margin = margin(0.5, 1, 0.5, F3_EV_LEFT_MARGIN_MM, unit = "mm"))
      else p
    })
    leg_c <- cowplot::get_legend(f3ev_map(dplyr::filter(d_c, event == ev$event[1]), maps, basemap, obc, bbs[[1]],
                                          show_legend = TRUE))
    # One flat design grid (a nested patchwork with fixed-aspect maps trips
    # patchwork's simplify_gt): the four top cells span every column, then one
    # b cell and one c cell per event, then c's legend under them.
    n_top <- length(top); ids <- c(LETTERS, letters)
    rows <- c(vapply(seq_len(n_top), function(k) strrep(ids[k], n_ev), ""),
              paste(ids[n_top + seq_len(n_ev)], collapse = ""),
              paste(ids[n_top + n_ev + seq_len(n_ev)], collapse = ""),
              strrep(ids[n_top + 2 * n_ev + 1], n_ev))
    cells <- c(top, pb, pc, list(patchwork::wrap_elements(full = leg_c, ignore_tag = TRUE)))
    # Tags drawn inside each plot's own area (top-left), so the page keeps no
    # empty tag column left of the axis titles.
    fig <- patchwork::wrap_plots(cells, design = paste(rows, collapse = "\n"),
                                 heights = grid::unit(c(F3_EV_STRIP_MM,
                                                        F3_EV_HEIGHTS[c("a", "leg_a", "b", "c", "leg_bc")]),
                                                      c("mm", "mm", rep("null", 5)))) & TAG_THEME &
      theme(plot.tag.location = "plot", plot.tag.position = c(0, 1))
    # A note only when panel b has an unidentified scope.
    ni_note <- ps_af_unidentified_note(d_b)
    if (nzchar(ni_note))
      fig <- fig + patchwork::plot_annotation(caption = paste0("b: ", ni_note),
                                              theme = theme(plot.caption = element_text(size = BASE_SIZE - 1.5,
                                                                                        colour = "grey25", hjust = 0)))
    is_main <- identical(metric, "cases")
    save_figure(fig, if (is_main) "fig3_event_attribution" else sprintf("fig3_event_attribution_%s", metric),
                height_mm = F3_EV_HEIGHT_MM, subdir = if (is_main) NULL else "supplementary")
  }
}

fig4_metric_lab <- function(metric) if (metric == "cases") "cases" else "outbreak onsets"

fig4_scope_lab  <- function(level, enso_group) ifelse(level == "continent", "Africa", unname(GROUP_SHORT_FIG[enso_group]))

#' AF of Africa and the ENSO regions (af_posterior$scope): the adjusted
#' posterior median with its 50% / 95% interval, in %, scope = Figure 3's labels;
#' an unidentified scope takes its pooled value and raw range, flagged
#' `unidentified` (ps_af_posterior_cols()).
fig4_af_data <- function(post, metric, aset, groups) {
  lv <- f3_scope_levels(groups)
  post$scope %>%
    dplyr::filter(cell == FIG4_WINDOW, analog_set == aset, explosive_guard == F3_EV_GUARD, metric == !!metric,
                  level %in% c("continent", "region")) %>%
    dplyr::mutate(scope = fig4_scope_lab(level, enso_group)) %>%
    dplyr::filter(scope %in% lv) %>%
    ps_af_posterior_cols() %>%
    dplyr::transmute(scope = factor(scope, levels = lv), level, af = af_bar,
                     af_lo, af_q25, af_q75, af_hi, unidentified)
}

#' a's ratio column: ratio of means = 1 / (1 - AF) with the AF interval mapped
#' the same way (monotone for AF < 1), and the scope's IPCC class at m.
fig4_ratio_data <- function(af, ipcc, metric, aset, m) {
  inv <- function(x) ifelse(x < 100, 100 / (100 - x), Inf)
  cls <- ipcc %>%
    dplyr::filter(window == FIG4_WINDOW, explosive_guard == F3_EV_GUARD, metric == !!metric, analog_set == aset,
                  abs(magnitude - m) < 1e-9, level %in% c("continent", "region")) %>%
    dplyr::transmute(scope = fig4_scope_lab(level, enso_group), class)
  lv <- levels(af$scope)
  af %>% dplyr::mutate(scope = as.character(scope)) %>% dplyr::left_join(cls, by = "scope") %>%
    dplyr::mutate(scope = factor(scope, levels = lv),
                  ratio = inv(af), ratio_lo = inv(af_lo), ratio_hi = inv(af_hi),
                  ratio_q25 = inv(af_q25), ratio_q75 = inv(af_q75),
                  class = factor(class, levels = names(FIG4_CLASS_COLS)))
}

#' a's ratio column: one row per scope (facets aligned with the time series),
#' bar from 1 to the ratio, 50% / 95% error bars, solid line at 1, dashed at
#' 1 + m. The axis is log10(ratio) on a linear scale (labels in ratio units), so
#' the clip arrows of ps_interval_layers() sit at the log-spaced edges.
fig4_ratio_col <- function(d, m) {
  rng <- range(c(d$ratio_lo, d$ratio_hi, 1 / 2, 1 + m), na.rm = TRUE, finite = TRUE)   # 1/2: line at 1 off the border
  lims <- c(max(FIG4_RATIO_BRK[FIG4_RATIO_BRK <= rng[1]], FIG4_RATIO_BRK[1]),
            min(FIG4_RATIO_BRK[FIG4_RATIO_BRK >= rng[2]], max(FIG4_RATIO_BRK)))
  keep <- FIG4_RATIO_BRK >= lims[1] & FIG4_RATIO_BRK <= lims[2]
  lg <- function(x) log10(pmin(pmax(x, lims[1]), lims[2]))
  d <- dplyr::mutate(d, y0 = 0, x_bar = lg(ratio), x_lo = log10(ratio_lo), x_hi = log10(ratio_hi),
                     x_q25 = lg(ratio_q25), x_q75 = lg(ratio_q75))
  ggplot(d, aes(y = y0)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.25) +
    geom_rect(aes(xmin = 0, xmax = x_bar, ymin = -0.6, ymax = 0.6, fill = class), colour = "grey15",
              linewidth = 0.2, show.legend = TRUE) +
    geom_vline(xintercept = log10(1 + m), colour = "grey30", linewidth = 0.25, linetype = "22") +
    ps_af_interval_layers("x_bar", "x_lo", "x_hi", "x_q25", "x_q75", orientation = "y", pos = "y0",
                          clip = log10(lims)) +
    facet_grid(scope ~ .) +
    scale_fill_manual(values = FIG4_CLASS_COLS, labels = FIG4_CLASS_LAB, drop = FALSE,
                      name = sprintf("IPCC, m = %g%%", 100 * m),
                      guide = guide_legend(override.aes = list(linewidth = 0.2), nrow = 2, title.position = "top")) +
    scale_x_continuous(limits = log10(lims), breaks = log10(FIG4_RATIO_BRK[keep]), labels = FIG4_RATIO_LAB[keep],
                       expand = c(0.04, 0), oob = scales::oob_keep) +
    scale_y_continuous(expand = expansion(add = F3_C_ARM_EXPAND)) +
    labs(x = "Forecast ÷\ncounterfactual, mean", y = NULL) +
    theme_ms() +
    theme(strip.text.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = BASE_SIZE - 1.5), axis.title.x = element_text(size = BASE_SIZE - 1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
}

#' The legend cell over a's two right columns: the IPCC class keys and a note
#' on the interval glyphs and the dashed 1 + m line.
fig4_ratio_legend <- function(d, m) {
  leg <- cowplot::get_legend(fig4_ratio_col(d, m) + F3_LEGEND_THEME + ps_legend_theme() +
                               theme(legend.position = "bottom", legend.direction = "horizontal",
                                     legend.key.width = unit(3, "mm"), legend.key.spacing.x = unit(1.5, "mm")))
  # The cumulative column shows medians across paired replicates (09
  # attribution cum_totals); only the ratio column is posterior (af_posterior).
  ni_note <- ps_af_unidentified_note(d, sep = "\n")
  note <- grid::textGrob(paste0("Cumulative: median across replicates; ratio: posterior\n",
                                "median; error bars 50% (thick), 95% (thin)   ┆ 1 + m",
                                if (nzchar(ni_note)) paste0("\n", ni_note) else ""),
                         x = 0, hjust = 0, gp = grid::gpar(fontsize = BASE_SIZE - 1.5, col = "grey25",
                                                           fontfamily = FONT))
  cowplot::plot_grid(leg, note, ncol = 1, rel_heights = c(1.2, 1))
}

#' c: AF of Africa and the ENSO regions as horizontal bars (the posterior
#' median) with 50% / 95% error bars; Africa on top, set apart by the rule
#' (ps_africa_rows()); x from 0 rounded out to 25 within FIG4_AF_XLIM, arrows
#' where a 95% interval runs past it.
fig4_af_panel <- function(d, metric, groups) {
  xl <- range(c(0, d$af_lo, d$af_hi), na.rm = TRUE)
  xl <- pmax(pmin(c(floor(xl[1] / 25) * 25, ceiling(xl[2] / 25) * 25), FIG4_AF_XLIM[2]), FIG4_AF_XLIM[1])
  rows <- ps_africa_rows(d$scope)
  d <- dplyr::mutate(d, yv = rows$y, af_bar = pmax(pmin(af, xl[2]), xl[1]), fill_key = as.character(scope))
  ps_true_minus(
    ggplot(d, aes(y = yv, x = af_bar, fill = fill_key)) +
      geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25) +
      geom_col(width = 0.68, orientation = "y") +
      ps_af_interval_layers("af_bar", "af_lo", "af_hi", "af_q25", "af_q75", orientation = "y", pos = "yv",
                            clip = xl) +
      rows$rule +
      scale_y_continuous(breaks = rows$breaks, labels = rows$labels, expand = expansion(add = 0.6)) +
      scale_fill_manual(values = f3_scope_pal(groups), guide = "none") +
      scale_x_continuous(breaks = seq(-100, 100, 25)) +
      coord_cartesian(xlim = xl) +
      labs(y = NULL, x = sprintf("ENSO-attributable fraction of %s (%%)", fig4_metric_lab(metric))) +
      theme_ms() +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(size = BASE_SIZE), axis.text.x = element_text(size = BASE_SIZE - 1.5),
            # transparent: c's background must not cover its neighbours' labels
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = margin(4, 2, 0.5, 2, unit = "mm")))
}

#' b: the forecast / counterfactual ADM1 ratio map over the whole of Africa
#' (attrib_ratio_map_layers(), shared with fig. AE), Figure 2's colour bar, the
#' scale bar + north arrow, legends under the map.
fig4_map <- function(adm1, metric, aset, m, maps, basemap, grey_sf) {
  map_sf <- attrib_ratio_map_sf(adm1, metric, aset, m, maps$adm1_sf, window = FIG4_WINDOW, guard = F3_EV_GUARD)
  bbox <- sf::st_bbox(maps$afro_sf)
  hatch <- ipcc_hatch_fills("Figure 4b")
  fill_scale <- ps_ratio_fill_scale(sprintf("Forecast ÷ neutral\ncounterfactual (mean %s)", fig4_metric_lab(metric)))
  ggplot() +
    ps_country_grey_layer(grey_sf, legend_order = 3) +
    attrib_ratio_map_layers(map_sf, c(FIG2_MAP_LO, FIG2_MAP_HI), fill_scale, hatch, m, hatch_nrow = 1L) +
    map_overlay_layers(basemap, maps$afro_sf) +
    map_scale_layers() +
    coord_sf(xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")], expand = FALSE) +
    taxdat::map_theme() + ps_map_font() +
    theme(panel.border = element_blank(),
          plot.margin = margin(FIG4_MAP_TOP_MM, 1, FIG4_C_GAP_MM, 0, unit = "mm"),
          plot.background = element_rect(fill = NA, colour = NA),
          legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.justification = "left", legend.background = element_blank(),
          legend.title = element_text(size = BASE_SIZE - 1, vjust = 0.8, lineheight = 0.9,
                                      margin = margin(r = 4, unit = "mm")),   # clear of the bar's first label
          legend.text = element_text(size = BASE_SIZE - 1.5),
          legend.spacing.y = unit(0.6, "mm"), legend.box.spacing = unit(3, "mm"),   # clear of the scale bar
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(2.4, "mm"), legend.key.width = unit(3, "mm"))
}

#' The arm keys under a: factual / counterfactual line + point, one row, centred.
fig4_arm_legend <- function(aset, stitched = FALSE, attrib_launch = NULL) {
  al <- f3_arm_lab(aset, stitched)
  if (!is.null(attrib_launch)) {
    # Three keys: the displayed forecast (line only), the attribution factual
    # (dotted line; the filled-circle points of the cumulative panel) and the
    # counterfactual.
    lv  <- c("forecast", "attrib", "counterfactual")
    lab <- c(forecast = al[["forecast"]],
             attrib = sprintf("Factual launched %s (attribution)", format(as.Date(attrib_launch), "%b %Y")),
             counterfactual = al[["counterfactual"]])
    d_arm <- data.frame(arm = factor(lv, levels = lv), x = 1, y = 1)
    return(cowplot::get_legend(
      ggplot(d_arm) +
        geom_line(aes(x, y, colour = arm, group = arm), linewidth = 0.45) +
        geom_point(aes(x, y, colour = arm), size = 1.8) +
        scale_colour_manual(values = c(forecast = "grey10", attrib = "grey10", counterfactual = F3_CF_COL),
                            labels = lab, name = NULL) +
        # one key per row: the three keys on one row overran the legend cell
        # (clipped at the left edge of the figure), as did two rows
        guides(colour = guide_legend(ncol = 1, override.aes = list(
          linetype = unname(c(F3_ARM_LTY[["forecast"]], F4_ATTRIB_LTY, F3_ARM_LTY[["counterfactual"]])),
          shape = c(NA, F3_ARM_SHAPE[["forecast"]], F3_ARM_SHAPE[["counterfactual"]])))) +
        theme_ms() +
        theme(legend.position = "bottom", legend.key.width = unit(7, "mm"), legend.spacing.x = unit(2, "mm"),
              legend.spacing.y = unit(0.5, "mm")) +
        ps_legend_theme() +
        theme(legend.justification = "center", legend.box.margin = margin(0, 0, 0, 0))))
  }
  d_arm <- data.frame(arm = factor(names(F3_ARM_LTY), levels = names(F3_ARM_LTY)), x = 1, y = 1)
  cowplot::get_legend(
    ggplot(d_arm) +
      geom_line(aes(x, y, linetype = arm, group = arm), linewidth = 0.45) +
      geom_point(aes(x, y, shape = arm), size = 1.8) +
      scale_linetype_manual(values = F3_ARM_LTY, labels = al, name = NULL) +
      scale_shape_manual(values = F3_ARM_SHAPE, labels = al, name = NULL) +
      guides(linetype = guide_legend(nrow = 1, override.aes = list(colour = c("grey10", F3_CF_COL))),
             shape    = guide_legend(nrow = 1, override.aes = list(colour = c("grey10", F3_CF_COL)))) +
      theme_ms() +
      theme(legend.position = "bottom", legend.key.width = unit(7, "mm"), legend.spacing.x = unit(2, "mm"),
            legend.spacing.y = unit(0.5, "mm")) +
      ps_legend_theme() +
      theme(legend.justification = "center", legend.box.margin = margin(0, 0, 0, 0)))
}

#' The group keys under c (two rows of three), left edge on c's panel edge by
#' construction: an invisible copy of c's left side -- same plot margin and y
#' labels (transparent), no tag (the grid's tag column is shared) -- freed like
#' c, so its panel starts where c's does on any device; the legend sits inside
#' that panel at its left.
fig4_group_legend_plot <- function(af, groups) {
  rows <- ps_africa_rows(af$scope)
  d <- dplyr::mutate(af, yv = rows$y)
  ggplot(d, aes(x = 0, y = yv, colour = scope)) +
    geom_point(shape = 15, alpha = 0) +
    scale_colour_manual(values = f3_scope_pal(groups), name = NULL,
                        guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 2.5))) +
    scale_y_continuous(breaks = rows$breaks, labels = rows$labels) +
    theme_ms() +
    ps_legend_theme() +
    theme(axis.text.y = element_text(size = BASE_SIZE, colour = NA), axis.ticks.y = element_line(colour = NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(),
          # every panel element blank (theme_ms() sets panel.grid.major itself)
          panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = NA, colour = NA),
          plot.margin = margin(0, 2, 0, 2, unit = "mm"),
          legend.position = "inside", legend.position.inside = c(0, 0.5), legend.justification = c(0, 0.5),
          legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0),
          legend.key.spacing.x = unit(5, "mm"), legend.key.width = unit(3, "mm"))
}

build_fig4 <- function(metric = "cases", aset = MAIN_SET) {
  message(sprintf("[10] Figure 4: metric=%s set=%s", metric, aset))
  f23  <- read_results("attribution", sets = SETS)
  oni  <- read_results("oni", sets = SETS)
  ae   <- read_results("report_attrib_exposure", sets = SETS)
  am   <- read_results("report_attrib_metrics", sets = SETS)
  post <- read_results("af_posterior", sets = SETS)
  m <- as.numeric(RES$style$attrib_headline$magnitude)
  if (length(m) != 1L || !any(abs(am$magnitudes - m) < 1e-9))
    stop("[10] Figure 4: results.style.attrib_headline$magnitude must be one of thresholds.attrib_magnitudes {",
         paste(am$magnitudes, collapse = ", "), "}.", call. = FALSE)
  groups <- GROUP_ORDER[GROUP_ORDER %in% f23$groups]
  af <- fig4_af_data(post, metric, aset, groups)
  d_ratio <- fig4_ratio_data(af, am$ipcc_scope, metric, aset, m)
  if (anyNA(d_ratio$class))
    stop("[10] Figure 4: no IPCC class for scope(s) ", paste(d_ratio$scope[is.na(d_ratio$class)], collapse = ", "),
         " in report_attrib_metrics$ipcc_scope.", call. = FALSE)

  # --- a: the ONI forecast strip, Figure 3's forecast panel, cumulative + ratio
  x_limits <- c(as.Date(f23$attrib_start), max(f23$region_ribbon$month))
  caps <- if (identical(metric, "cases")) F3_CASES_YMAX else c(africa = NULL, region = NULL)
  oni_b <- blank_left_strip(ps_true_minus(fig3_oni_strip(oni, aset, x_limits, mode = "forecast")))
  oni_b <- patchwork::free(oni_b + theme(axis.title.y = element_text(margin = margin(r = 0.5, unit = "mm"))),
                           type = "label", side = "l") + labs(tag = "a")
  # With an attribution factual (09: results.inputs.attrib_manifest) the paired
  # statistics -- c's cumulative totals, the ratio column, b and c -- compare it
  # with the counterfactual; b keeps the displayed forecast and adds its median.
  attrib <- identical(f23$attrib_factual %||% "forecast", "attrib")
  attrib_launch <- if (attrib) as.Date(f23$attrib_arm_starts$factual %||% f23$attrib_start) else NULL
  pb <- fig3_panel_b(f23, metric, aset, groups, x_limits,
                     africa_ymax = caps[["africa"]], region_ymax = caps[["region"]],
                     attrib_series = attrib) +
    theme(legend.position = "none")
  # The headline cell's cumulative totals (attribution v10 cum_totals_cells;
  # cum_totals is the attribution x drop cell the legacy figures read).
  if (!FIG4_WINDOW %in% f23$cum_totals_cells$window)
    stop("[10] Figure 4: results.style.attrib_headline$window '", FIG4_WINDOW, "' is not a cell of the attribution ",
         "artefact {", paste(unique(f23$cum_totals_cells$window), collapse = ", "), "}.", call. = FALSE)
  cum_hl <- dplyr::filter(f23$cum_totals_cells, window == FIG4_WINDOW, explosive_guard == F3_EV_GUARD, metric == !!metric)
  common_win <- grepl("_common$", FIG4_WINDOW) && identical(f23$africa_window$mode, "common")
  pc <- fig3_panel_c(cum_hl, aset, groups, f23$attrib_start,
                     if (attrib) f23$attrib_start else f23$factual_start %||% f23$attrib_start) +
    scale_x_log10(labels = PS_LABEL_K, breaks = FIG4_CUM_BRK) +
    labs(x = if (common_win) sprintf("Cumulative %s,\n%s", fig4_metric_lab(metric), f23$africa_window$label)
             else sprintf("Cumulative %s,\n%s–horizon", fig4_metric_lab(metric), format(as.Date(f23$attrib_start), "%b %Y"))) +
    theme(legend.position = "none", strip.text.x.bottom = element_blank(),
          axis.title.x = element_text(size = BASE_SIZE - 1, margin = margin(t = 0.3, unit = "mm")))
  # fig3_panel_c()'s 95% / 50% lineranges become Figure 4's capped error bars in
  # the arms' colours; its point layer (shape by arm) stays, at PS_ERRORBAR's size.
  pts <- Filter(function(l) inherits(l$geom, "GeomPoint"), pc$layers)
  if (length(pts) != 1L)
    stop("[10] Figure 4: fig3_panel_c() must have exactly one point layer, found ", length(pts), call. = FALSE)
  pts <- pts[[1L]]
  pts$aes_params$size <- PS_ERRORBAR$point_size; pts$aes_params$stroke <- PS_ERRORBAR$point_stroke
  pc$layers <- c(ps_interval_layers("med", "lo", "hi", "q25", "q75", orientation = "y", colour = NULL,
                                    point = FALSE),
                 list(pts))
  pr <- fig4_ratio_col(d_ratio, m)
  leg_ratio <- patchwork::wrap_elements(full = fig4_ratio_legend(d_ratio, m), ignore_tag = TRUE)
  leg_arms  <- patchwork::wrap_elements(full = fig4_arm_legend(aset, f3_stitched(f23), attrib_launch),
                                        ignore_tag = TRUE)

  # --- b, c -----------------------------------------------------------------
  maps <- load_africa_maps("Figure 4")
  basemap <- africa_basemap(maps$afro_sf)
  grey_sf <- ps_country_status(maps$afro_sf, modelled_countries())
  # Wrapped: its fixed aspect would otherwise set the right column's width.
  map <- patchwork::wrap_elements(full = fig4_map(ae$adm1, metric, aset, m, maps, basemap, grey_sf) +
                                    labs(tag = "b") + TAG_THEME)
  af_p <- fig4_af_panel(af, metric, groups) + labs(tag = "c")
  if (!is.null(FIG4_INSET_BOX))   # the group map (Figure 3's), the fit registry's partition (09)
    af_p <- f34_add_group_map(af_p, fig2_group_map_inset(maps, read_results("modelled_units")$members),
                              FIG4_INSET_BOX)
  leg_groups <- patchwork::free(fig4_group_legend_plot(af, groups))

  a_margin <- theme(plot.margin = margin(1, 1.5, 1.5, 1.5, unit = "mm"))
  pb <- pb + a_margin; pc <- pc + a_margin; pr <- pr + a_margin
  # b's two-line date labels set the axis row: free the other columns' bottoms.
  # c is freed whole: it fills its cell, its axis title ends on a's bottom line
  # and its labels stay inside its cell (not in a gutter column a would pay for).
  pc <- patchwork::free(pc, type = "space", side = "b")
  pr <- patchwork::free(pr, type = "space", side = "b")
  af_p <- patchwork::free(af_p)
  fig <- patchwork::wrap_plots(A = oni_b, B = leg_ratio, C = map, D = pb, E = pc, F = pr, G = af_p, H = leg_arms,
                               I = leg_groups, design = "ABBC\nDEFC\nDEFG\nHHHI", widths = FIG4_WIDTHS,
                               heights = grid::unit(c(PS_ONI_STRIP_MM, FIG4_HEIGHTS),
                                                    c("mm", "null", "null", "null"))) & TAG_THEME

  is_main <- identical(metric, "cases") && identical(aset, MAIN_SET)
  stem   <- paste(c("fig4_forecast_attribution", if (metric != "cases") metric, if (aset != MAIN_SET) aset),
                  collapse = "_")
  subdir <- if (is_main) NULL else "supplementary"
  save_figure(fig, stem, height_mm = FIG4_HEIGHT_MM, subdir = subdir)
  invisible(fig)
}
