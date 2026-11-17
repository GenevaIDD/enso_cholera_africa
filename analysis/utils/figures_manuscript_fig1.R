# figures_manuscript_fig1.R — Figure 1 builders of 10_manuscript_figures.R
#   (ENSO exposure-response overview).
#
# Provides:
#   build_fig1()             a coverage timeline + ONI strip, b inset map,
#                            c seasonality grid with d (ERF) as its last row
#   build_fig1_v91()         side-by-side layout: c lag weights over the
#                            rainfall-ONI correlation, d ERF grid
#   build_fig1_inset_map(), add_fig1_inset()   Figure 1b map
#   fig1c_grid()             Figure 1c panel grid
#
# Details:
#   Definitions only (loads alone). The builders are 10's own: at CALL time
#   they read 10's globals -- style (FONT, BASE_SIZE, theme_ms(), TAG_THEME,
#   PS_*), groups (GROUP_*), config (RES, SCOPE, SETS, MAIN_SET, opt), the
#   figure constants and 10's helpers (read_results(), load_africa_maps(),
#   save_figure(), ...). Only 10_manuscript_figures.R sources them.
#   The FIG1_* layout constants stay in 10 next to its style (sections 2 / 5a).

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

  gid_sf <- read_gadm_units(rate$country) %>%
    dplyr::filter(gid %in% rate$gid) %>%
    sf::st_make_valid() %>%
    rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE) %>%
    dplyr::inner_join(rate, by = "gid")
  missing <- setdiff(rate$gid, gid_sf$gid)
  if (length(missing)) {
    message(sprintf(
      "[10] Figure 1b: %d gid(s) without a GADM polygon: %s",
      length(missing),
      paste(utils::head(missing, 10), collapse = ",")
    ))
  }

  afro_sf <- ps_country_status(
    drop_small_islands(load_africa_sf(gadm_path = FIG1_GADM_PATH)$afro_sf),
    modelled_countries(rate)
  )
  grp_sf <- afro_sf %>%
    dplyr::inner_join(
      dplyr::distinct(rate, country, enso_group),
      by = c("ISO3_CODE" = "country")
    ) %>%
    dplyr::filter(enso_group %in% groups)

  ggplot() +
    ps_country_grey_layer(
      afro_sf,
      legend_order = 2,
      labels = PS_GREY_LAB_WRAP,
      show_legend = FALSE
    ) +
    geom_sf(data = gid_sf, aes(fill = enso_group, alpha = rate), colour = NA) +
    geom_sf(data = afro_sf, fill = NA, colour = "grey55", linewidth = 0.1) +
    geom_sf(
      data = grp_sf,
      aes(colour = enso_group),
      fill = NA,
      linewidth = 0.12
    ) +
    # No outline for the uncertain-membership countries: panel a's
    # daggers mark them.
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_alpha_continuous(
      range = c(0.1, 1),
      transform = "sqrt",
      name = "Onsets per\ncovered month",
      labels = scales::label_percent(accuracy = 1),
      guide = guide_legend(order = 1, nrow = 2)
    ) +
    coord_sf(
      xlim = FIG1_INSET_BBOX[c("xmin", "xmax")],
      ylim = FIG1_INSET_BBOX[c("ymin", "ymax")],
      expand = FALSE,
      datum = NA
    ) +
    labs(tag = "b") +
    theme_void(base_size = 5, base_family = FONT) +
    # Compact legend (title above a two-row rate key) in the South Atlantic
    # corner of the inset.
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.02),
      legend.justification = c(0, 0),
      legend.key.size = unit(1.3, "mm"),
      legend.spacing.y = unit(0.4, "mm"),
      legend.text.position = "right",
      legend.box.just = "left",
      legend.title.position = "top",
      legend.title = element_text(size = 4.5, family = FONT, vjust = 0.5),
      legend.text = element_text(size = 4, family = FONT),
      plot.background = element_rect(
        fill = "white",
        colour = "grey60",
        linewidth = 0.2
      ),
      plot.tag = element_text(
        size = BASE_SIZE + 2,
        face = "bold",
        family = FONT
      ),
      plot.tag.position = c(0.06, 0.95),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
    )
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
  xr <- pp$x.range
  yr <- pp$y.range
  n <- length(country_levels)
  right <- 1 - FIG1_INSET_MARGIN
  top <- 1 - FIG1_INSET_MARGIN * FIG1_PANEL_A_ASPECT
  box_for <- function(k) {
    y_bottom <- n - k + 0.5
    bottom <- (y_bottom - yr[1]) / diff(yr)
    left <- right - (top - bottom) / (FIG1_INSET_ASPECT * FIG1_PANEL_A_ASPECT)
    covered <- country_levels[seq_len(n) > y_bottom]
    lb <- last_bubble[as.character(covered)]
    list(
      left = left,
      bottom = bottom,
      covered = covered,
      ok = all(is.na(lb) | lb < xr[1] + left * diff(xr))
    )
  }
  ks <- rev(seq_len(min(FIG1_INSET_MAX_ROWS, n)))
  k <- ks[vapply(ks, function(k) box_for(k)$ok, logical(1))][1]
  if (is.na(k)) {
    warning(
      "[10] Figure 1b: no inset size avoids panel a's case bubbles; inset omitted"
    )
    return(p_cov)
  }
  bx <- box_for(k)
  message(sprintf(
    "[10] Figure 1b inset: %d row(s) (%s), panel fractions left=%.3f bottom=%.3f",
    k,
    paste(rev(as.character(bx$covered)), collapse = ","),
    bx$left,
    bx$bottom
  ))
  p_cov +
    patchwork::inset_element(
      inset,
      left = bx$left,
      bottom = bx$bottom,
      right = right,
      top = top,
      align_to = "panel",
      clip = TRUE
    )
}

#' Figure 1a's ONI strip: Figure 3a's observed-ONI strip (fig3_oni_strip())
#' over panel a's window, carrying the "a" tag, with tighter margins and panel
#' a's tick size.
fig1_oni_strip <- function(oni, x_limits, phases) {
  fig3_oni_strip(oni, NULL, x_limits, mode = "observed", phases = phases) +
    labs(tag = "a") +
    theme(
      axis.title.y = element_text(margin = margin(r = 0.5, unit = "mm")),
      plot.margin = margin(0.5, 0.5, 0, 0.5, unit = "mm")
    ) +
    theme(axis.text.y = element_text(size = FIG1_A_TICK_PT)) +
    # The tag inside the strip's top-left corner: a tag column would also be
    # reserved left of a's country labels (white space).
    theme(
      plot.tag.position = c(0, 1),
      plot.tag.location = "plot",
      plot.tag = element_text(
        size = BASE_SIZE + 2,
        face = "bold",
        family = FONT,
        hjust = 0,
        vjust = 1
      )
    )
}

#' Group short labels for Figure 1's group strips.
fig1_group_labeller <- function(membership_info = FALSE) {
  lab <- GROUP_SHORT_FIG
  m <- if (isTRUE(membership_info)) membership() else NULL
  if (!is.null(m)) {
    g <- m$group[match(names(lab), m$group$group), , drop = FALSE]
    ok <- !is.na(g$group)
    lab[ok] <- sprintf(
      "%s\nn %d \u00b7 p\u0303 %.2f",
      lab[ok],
      g$n_hard[ok],
      g$p_max_median[ok]
    )
  }
  ggplot2::as_labeller(lab)
}

#' Horizontal colourbar for Figure 1's c / d legend row, title left of the bar.
#' `text`: an optional legend.text element (left unset when NULL -- an explicit
#' NULL in theme() would blank the element and drop the labels to the default size).
fig1_bar <- function(bar_mm = FIG1_BAR_MM, text = NULL) {
  th <- theme(
    legend.title.position = "left",
    legend.key.width = unit(bar_mm, "mm"),
    legend.key.height = unit(1.6, "mm"),
    legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.9)
  )
  if (!is.null(text)) {
    th <- th + theme(legend.text = text)
  }
  guide_colourbar(direction = "horizontal", theme = th)
}

#' The occurrence track's selected lag per group (09's lag_rank, `selected`).
#' `lag_rank`: the tracks' lag_rank tables bound with a `track` column.
fig1_best_occ_lag <- function(lag_rank, groups) {
  lag_rank %>%
    dplyr::filter(track == "occurrence", selected, enso_group %in% groups) %>%
    dplyr::transmute(
      enso_group = factor(enso_group, levels = GROUP_ORDER),
      lag = single_lag
    )
}

#' ELPD-implied lag weights of the occurrence track per group: pseudo-BMA,
#' exp(elpd_diff) normalised over the lag grid (the weight the ELPD alone
#' implies, no standard error).
fig1_elpd_weights <- function(lag_rank, groups) {
  lag_rank %>%
    dplyr::filter(track == "occurrence", enso_group %in% groups) %>%
    dplyr::group_by(enso_group) %>%
    dplyr::mutate(w = exp(elpd_diff - max(elpd_diff)), w = w / sum(w)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      enso_group = factor(enso_group, levels = GROUP_ORDER),
      lag = single_lag,
      w,
      selected,
      header = FIG1_HEADER
    )
}

#' x range shared by c's weight strip and tiles: the lag tiles, the gap, the
#' onset column and its label's room, FIG1_X_PAD beyond on both sides.
fig1_cd_xlim <- function(lag_range) {
  c(
    lag_range[1] - 0.5 - FIG1_X_PAD,
    lag_range[2] +
      FIG1_ONSET_GAP +
      FIG1_ONSET_W -
      0.5 +
      FIG1_X_PAD +
      FIG1_ONSET_LABEL_PAD
  )
}

#' The onset column's centre (lag units).
fig1_onset_x <- function(lag_range) {
  lag_range[2] + FIG1_ONSET_GAP + (FIG1_ONSET_W - 1) / 2
}

#' Frames in data coordinates: the lag range's box (FIG1_X_PAD beyond the
#' tiles, so the rainy-season outlines on the tiles' edges stay visible inside
#' it) and, with `onset = TRUE`, the onset column's; full panel height. Pair
#' with theme(panel.border = element_blank()).
fig1_frame_layers <- function(
  lag_range,
  onset = FALSE,
  lwd = FIG1_FRAME$linewidth
) {
  box <- function(xmin, xmax) {
    annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = NA,
      colour = FIG1_FRAME$colour,
      linewidth = lwd
    )
  }
  out <- list(box(
    lag_range[1] - 0.5 - FIG1_X_PAD,
    lag_range[2] + 0.5 + FIG1_X_PAD
  ))
  if (onset) {
    ox <- fig1_onset_x(lag_range)
    out <- c(out, list(box(ox - FIG1_ONSET_W / 2, ox + FIG1_ONSET_W / 2)))
  }
  out
}

#' Figure 1c, top: one tile row per group, each lag's ELPD-implied weight
#' (fig1_elpd_weights()) as the group colour's opacity on a log scale (white
#' below FIG1_W_LIM[1]), the selected lag's line, the lead rule, under the
#' nested "ENSO group" header. Its weight legend (alpha, grey keys) is built at
#' the bottom for fig1_legend_row() to lift.
fig1_weight_panel <- function(wts, best, lag_range) {
  strip_margin <- margin(t = FIG1_STRIP_GAP_MM, b = 0.6, unit = "mm")
  # Every layer carries the nested header's column, or facet_nested() opens a
  # second, header-less set of group panels for it.
  best <- dplyr::mutate(best, header = FIG1_HEADER)
  w_guide <- guide_legend(
    direction = "horizontal",
    nrow = 1,
    order = 1,
    override.aes = list(fill = "grey25"),
    theme = theme(
      legend.title.position = "left",
      legend.key.size = unit(FIG1_W_KEY_MM, "mm"),
      legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.9)
    )
  )
  ggplot(wts, aes(x = lag, y = 1)) +
    geom_tile(aes(fill = enso_group, alpha = w), width = 1, height = 1) +
    geom_vline(
      data = best,
      aes(xintercept = lag),
      colour = FIG1_BEST_COL,
      linewidth = FIG1_BEST_LWD,
      linetype = FIG1_BEST_LTY
    ) +
    geom_vline(
      xintercept = FIG1_LEAD_RULE,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    fig1_frame_layers(lag_range, lwd = FIG1_W_FRAME_LWD) +
    ggh4x::facet_nested(
      . ~ header + enso_group,
      labeller = labeller(
        enso_group = fig1_group_labeller(),
        header = label_value
      ),
      nest_line = element_line(colour = "black", linewidth = 0.25),
      strip = ggh4x::strip_nested(
        by_layer_x = TRUE,
        text_x = ggh4x::elem_list_text(
          size = c(FIG1_HEADER_PT, BASE_SIZE - 1),
          margin = list(strip_margin, strip_margin)
        )
      )
    ) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_alpha_continuous(
      range = c(0, 1),
      transform = "log10",
      limits = FIG1_W_LIM,
      oob = scales::squish,
      breaks = FIG1_W_BREAKS,
      labels = format(FIG1_W_BREAKS, scientific = FALSE, drop0trailing = TRUE),
      name = FIG1_W_TITLE,
      guide = w_guide
    ) +
    scale_x_continuous(
      breaks = FIG1_LAG_BREAKS,
      limits = fig1_cd_xlim(lag_range),
      expand = c(0, 0),
      oob = scales::oob_keep
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, tag = "c") +
    theme_ms() +
    ps_legend_theme() +
    FIG1_CD_THEME +
    FIG1_PANEL_THEME +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      strip.text.x = element_text(
        size = BASE_SIZE - 1,
        face = "plain",
        margin = strip_margin
      ),
      # No x axis (the tiles' serves both); the tag gutter on the left.
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "mm"),
      plot.margin = margin(1.5, 1.5, 0, FIG1_TAG_GUTTER_MM, unit = "mm"),
      legend.position = "bottom"
    )
}

#' Figure 1c, tiles: rainfall anomaly-ONI correlation per rainfall month (rows)
#' and lag, per ENSO teleconnection group (columns), from 09's rain_enso
#' artefact, cropped to `lag_range`; the rainy seasons outlined (keyed "Rainy
#' season"); the selected occurrence lag (`best`); right of each group's tiles
#' the occurrence seasonal baseline per calendar month (`seas`, 09's
#' occ_seasonality), labelled on the leftmost group; under the leftmost group's
#' axis the lead / lag key on the lead rule.
fig1_corr_panel <- function(
  groups,
  lag_range,
  best,
  seas,
  rain = read_results("rain_enso")
) {
  st <- rain$corr %>%
    dplyr::filter(
      series %in% groups,
      lag >= lag_range[1],
      lag <= lag_range[2]
    ) %>%
    dplyr::mutate(
      enso_group = factor(series, levels = GROUP_ORDER),
      centre = factor(month.abb[centre], levels = rev(month.abb)),
      mark = dplyr::case_when(
        p_r < .05 & agree_r >= .5 ~ "*",
        p_r < .05 ~ "·",
        TRUE ~ ""
      )
    )
  miss <- setdiff(groups, unique(st$series))
  if (length(miss)) {
    stop(
      "[10] Figure 1c: the rain_enso artefact has no series for group(s) ",
      paste(miss, collapse = ","),
      call. = FALSE
    )
  }
  # Rainy seasons: one row per season (first .. last calendar month, possibly
  # wrapping the year); one rectangle per run of consecutive months on the y
  # positions of rev(month.abb), so a Dec -> Jan season gives two.
  seasons <- rain$seasons %>%
    dplyr::filter(weighting == FIG1_RAIN_SEASON_WEIGHTING, group %in% groups)
  miss <- setdiff(groups, unique(seasons$group))
  if (length(miss)) {
    stop(
      "[10] Figure 1c: the rain_enso artefact has no ",
      FIG1_RAIN_SEASON_WEIGHTING,
      "-weighted season for group(s) ",
      paste(miss, collapse = ","),
      call. = FALSE
    )
  }
  seasons <- seasons %>%
    dplyr::rowwise() %>%
    dplyr::reframe(
      series = group,
      season = season,
      month = ((first - 1) + 0:((last - first) %% 12)) %% 12 + 1
    ) %>%
    dplyr::group_by(series, season) %>%
    dplyr::mutate(run = cumsum(c(1, diff(month) != 1))) %>%
    dplyr::mutate(y = 13 - month) %>%
    dplyr::group_by(series, season, run) %>%
    dplyr::summarise(
      ymin = min(y) - 0.5,
      ymax = max(y) + 0.5,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      enso_group = factor(series, levels = GROUP_ORDER),
      xmin = lag_range[1] - 0.5,
      xmax = lag_range[2] + 0.5,
      key = FIG1_RAIN_SEASON_LABEL
    )
  miss <- setdiff(groups, unique(seas$enso_group))
  if (length(miss)) {
    stop(
      "[10] Figure 1c: the occ_seasonality artefact has no group(s) ",
      paste(miss, collapse = ","),
      call. = FALSE
    )
  }
  ox <- fig1_onset_x(lag_range)
  onset <- seas %>%
    dplyr::filter(enso_group %in% groups) %>%
    dplyr::mutate(
      enso_group = factor(enso_group, levels = GROUP_ORDER),
      centre = factor(month.abb[month], levels = rev(month.abb)),
      x = ox
    )
  g1 <- factor(groups[1], levels = GROUP_ORDER)
  # The column's label, rotated, right of the leftmost group's column.
  lab <- data.frame(
    enso_group = g1,
    x = ox + FIG1_ONSET_W / 2 + FIG1_ONSET_LABEL_PAD / 2,
    y = 6.5
  )
  # Lead / lag key: "left | right" on one line, the bar on the lead rule;
  # arrows running out from the bar below it.
  key_text <- data.frame(
    enso_group = g1,
    x = FIG1_LEAD_RULE + c(-FIG1_KEY_BAR_GAP, 0, FIG1_KEY_BAR_GAP),
    y = FIG1_KEY_Y,
    hjust = c(1, 0.5, 0),
    label = c(FIG1_KEY_TEXT[["left"]], "|", FIG1_KEY_TEXT[["right"]])
  )
  key_arrows <- data.frame(
    enso_group = g1,
    x = FIG1_LEAD_RULE + c(-1, 1) * FIG1_KEY_BAR_GAP,
    xend = FIG1_LEAD_RULE + c(-1, 1) * FIG1_KEY_ARROW_LEN,
    y = FIG1_KEY_Y - FIG1_KEY_LINE_SKIP
  )
  onset_bar <- fig1_bar()
  onset_bar$params$order <- 3L
  corr_bar <- fig1_bar()
  corr_bar$params$order <- 1L
  ggplot(st, aes(x = lag, y = centre, fill = r)) +
    geom_tile() +
    geom_text(
      aes(label = mark),
      size = FIG1_MARK_PT / ggplot2::.pt,
      family = FONT,
      vjust = 0.75
    ) +
    geom_vline(
      xintercept = FIG1_LEAD_RULE,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    geom_vline(
      data = best,
      aes(xintercept = lag),
      colour = FIG1_BEST_COL,
      linewidth = FIG1_BEST_LWD,
      linetype = FIG1_BEST_LTY
    ) +
    geom_rect(
      data = seasons,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, linetype = key),
      inherit.aes = FALSE,
      fill = NA,
      colour = FIG1_RAIN_SEASON_COL,
      linewidth = FIG1_RAIN_SEASON_LWD
    ) +
    scale_linetype_manual(
      values = stats::setNames("solid", FIG1_RAIN_SEASON_LABEL),
      name = NULL,
      guide = guide_legend(
        order = 2,
        override.aes = list(
          fill = NA,
          colour = FIG1_RAIN_SEASON_COL,
          linewidth = 0.5
        )
      )
    ) +
    facet_grid(. ~ enso_group, labeller = fig1_group_labeller()) +
    scale_fill_gradient2(
      low = FIG1_CORR_COL[["low"]],
      mid = FIG1_CORR_COL[["mid"]],
      high = FIG1_CORR_COL[["high"]],
      midpoint = 0,
      limits = FIG1_CORR_LIM,
      oob = scales::squish,
      name = FIG1_CORR_TITLE,
      breaks = c(-0.5, 0, 0.5),
      guide = corr_bar
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = onset,
      aes(x = x, y = centre, fill = rel),
      inherit.aes = FALSE,
      colour = "grey70",
      linewidth = FIG1_ONSET_LWD,
      width = FIG1_ONSET_W
    ) +
    # viridis is sequential; with limits symmetric in log10 its middle is 1.
    scale_fill_viridis_c(
      option = FIG1_ONSET_OPTION,
      transform = "log10",
      limits = FIG1_ONSET_LIM,
      oob = scales::squish,
      breaks = FIG1_ONSET_BREAKS,
      labels = as.character(FIG1_ONSET_BREAKS),
      name = FIG1_ONSET_TITLE,
      guide = onset_bar
    ) +
    geom_text(
      data = lab,
      aes(x = x, y = y, label = FIG1_ONSET_LABEL),
      inherit.aes = FALSE,
      angle = -90,
      hjust = 0.5,
      vjust = 0.5,
      size = FIG1_ONSET_LABEL_PT / ggplot2::.pt,
      family = FONT,
      colour = "grey20"
    ) +
    fig1_frame_layers(lag_range, onset = TRUE) +
    scale_x_continuous(
      breaks = FIG1_LAG_BREAKS,
      limits = fig1_cd_xlim(lag_range),
      expand = c(0, 0),
      oob = scales::oob_keep
    ) +
    geom_text(
      data = key_text,
      aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      vjust = 0.5,
      size = FIG1_KEY_PT / ggplot2::.pt,
      family = FONT,
      colour = "grey20"
    ) +
    geom_segment(
      data = key_arrows,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "grey20",
      linewidth = FIG1_KEY_ARROW_LWD,
      arrow = grid::arrow(
        length = unit(FIG1_KEY_HEAD_MM, "mm"),
        type = "closed",
        angle = 25
      )
    ) +
    # The key sits below the 12 month rows: hold the y range there.
    coord_cartesian(ylim = c(0.5, 12.5), expand = FALSE, clip = "off") +
    scale_y_discrete(expand = expansion(add = 0)) +
    labs(x = FIG1_X_TITLE, y = "Rainfall month", tag = NULL) +
    theme_ms() +
    ps_legend_theme() +
    FIG1_CD_THEME +
    FIG1_PANEL_THEME +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_text(size = BASE_SIZE - 2.5),
      strip.text.x = element_blank(),
      plot.margin = margin(0, 1.5, 0, 0, unit = "mm"),
      axis.title.x = element_text(
        size = BASE_SIZE - 1,
        margin = margin(t = FIG1_XTITLE_GAP_MM, unit = "mm")
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = unit(3, "mm"),
      legend.key.spacing.x = unit(PS_KEY_SPACING_MM, "mm")
    )
}

#' Figure 1d: ERF grid, rows = track, columns = ENSO teleconnection group (50%
#' + 95% CrI) on the ENSO-phase shading; arrows in the first column say what
#' an effect above / below 1 means per track (e plots the cessation rate
#' itself: above 1 = ends sooner). `erf`: the tracks' erf tables bound with a
#' `track` column.
fig1_erf_panel <- function(erf) {
  erf_data <- erf %>%
    dplyr::mutate(
      enso_group = factor(enso_group, levels = GROUP_ORDER),
      track = track_factor(track, erf = TRUE, wrap = TRUE)
    )
  ann <- data.frame(
    track = factor(levels(erf_data$track), levels = levels(erf_data$track)),
    enso_group = factor(
      levels(erf_data$enso_group)[1],
      levels = levels(erf_data$enso_group)
    ),
    x = ERF_XLIM[1] + 0.12,
    up = c("more\nonsets", "shorter\noutbreaks", "more\ncases"),
    down = c("fewer\nonsets", "longer\noutbreaks", "fewer\ncases")
  )
  arrow_layer <- function(y, yend) {
    geom_segment(
      data = ann,
      aes(x = x, xend = x, y = y, yend = yend),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = unit(0.9, "mm"), type = "closed"),
      linewidth = 0.25,
      colour = "grey25"
    )
  }
  label_layer <- function(y, label, vjust) {
    geom_text(
      data = ann,
      aes(x = x + 0.28, y = y, label = .data[[label]]),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = vjust,
      lineheight = 0.85,
      size = (BASE_SIZE - 2.5) / ggplot2::.pt,
      family = FONT,
      colour = "grey25"
    )
  }
  ggplot(
    erf_data,
    aes(
      x = at,
      y = IRR,
      ymin = IRR_lo,
      ymax = IRR_hi,
      colour = enso_group,
      fill = enso_group
    )
  ) +
    # ENSO-phase bands at utils/enso.R's ±0.5 thresholds. y runs 0 .. 100 on
    # purpose: on the log10 axis 0 maps to -Inf (the panel's bottom), whereas
    # ymin = -Inf maps to NaN and drops the band; x 3 is past ERF_XLIM.
    annotate(
      "rect",
      xmin = -3,
      xmax = ENSO_PHASE_THRESHOLDS[["nina"]],
      ymin = 0,
      ymax = 100,
      fill = PS_ENSO_PHASE_COL[["nina"]],
      alpha = PS_ENSO_PHASE_ALPHA
    ) +
    annotate(
      "rect",
      xmin = ENSO_PHASE_THRESHOLDS[["nina"]],
      xmax = ENSO_PHASE_THRESHOLDS[["elnino"]],
      ymin = 0,
      ymax = 100,
      fill = PS_ENSO_PHASE_COL[["neutral"]],
      alpha = PS_ENSO_PHASE_ALPHA
    ) +
    annotate(
      "rect",
      xmin = ENSO_PHASE_THRESHOLDS[["elnino"]],
      xmax = 3,
      ymin = 0,
      ymax = 100,
      fill = PS_ENSO_PHASE_COL[["elnino"]],
      alpha = PS_ENSO_PHASE_ALPHA
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_ribbon(
      aes(ymin = IRR_lo50, ymax = IRR_hi50),
      alpha = 0.35,
      colour = NA
    ) +
    geom_line(linewidth = 0.6) +
    arrow_layer(1.15, 6) +
    arrow_layer(0.87, 0.36) +
    label_layer(6, "up", 1) +
    label_layer(0.36, "down", 0) +
    # Track strips on the left, outside the axis; no group strips (c's serve).
    facet_grid(
      track ~ enso_group,
      switch = "y",
      labeller = labeller(enso_group = fig1_group_labeller())
    ) +
    # One "ENSO group" key (colour and fill merged) under e.
    scale_colour_manual(
      values = GROUP_COLOURS,
      labels = GROUP_SHORT_FIG,
      name = GROUP_KEY_TITLE
    ) +
    scale_fill_manual(
      values = GROUP_COLOURS,
      labels = GROUP_SHORT_FIG,
      name = GROUP_KEY_TITLE
    ) +
    scale_y_log10(breaks = FIG1_ERF_BREAKS, labels = FIG1_ERF_LABELS) +
    scale_x_continuous(breaks = FIG1_ONI_BREAKS) +
    coord_cartesian(xlim = ERF_XLIM, ylim = FIG1_ERF_YLIM) +
    labs(
      x = "Oceanic Niño Index (ONI)",
      y = "Effect size on rate of...",
      tag = "d"
    ) +
    theme_ms() +
    ps_legend_theme() +
    FIG1_CD_THEME +
    FIG1_PANEL_THEME +
    theme(
      strip.text.x = element_blank(),
      # y title hugging the track strips; bottom margin aligns e's panel
      # bottom with a's.
      axis.title.y = element_text(
        size = BASE_SIZE - 1,
        margin = margin(r = 0.5, unit = "mm")
      ),
      plot.margin = margin(
        FIG1_E_TOP_MM,
        1.5,
        FIG1_E_BOTTOM_MM,
        0,
        unit = "mm"
      ),
      axis.title.x = element_text(
        margin = margin(t = FIG1_XTITLE_GAP_MM, unit = "mm")
      )
    )
}

#' Figure 1c's legends, cut from the two plots' own gtables (nothing is
#' collected: d keeps its legend in place): the weight strip's row above the
#' tiles' row (both colourbars and the rainy-season key), FIG1_LEGEND_GAP_MM
#' apart, each row centred in the block's width.
fig1_legend_row <- function(p_w, p_tiles) {
  box <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    g$grobs[[which(g$layout$name == "guide-box-bottom")]]
  }
  centre <- function(g) {
    gtable::gtable_row(
      "c",
      list(grid::nullGrob(), g, grid::nullGrob()),
      height = sum(g$heights),
      widths = grid::unit.c(unit(1, "null"), sum(g$widths), unit(1, "null"))
    )
  }
  gw <- box(p_w)
  gt <- box(p_tiles)
  gtable::gtable_col(
    "fig1_c_legends",
    list(centre(gw), grid::nullGrob(), centre(gt)),
    width = unit(1, "null"),
    heights = grid::unit.c(
      sum(gw$heights),
      unit(FIG1_LEGEND_GAP_MM, "mm"),
      sum(gt$heights)
    )
  )
}

#' Figure 1a's data-coverage panel: per-country outbreak-period timeline on the
#' ENSO-phase background, monthly suspected cases as bubbles, the inset map (b)
#' in its top-right corner. Reads the raw survival / coverage / monthly data
#' directly -- a narrow exception to the data contract above, scoped to this
#' one panel. `groups`: the ENSO groups Figure 1 shows (inset borders).
fig1_coverage_panel <- function(phases, groups) {
  enso_groups <- readr::read_csv(GROUPING_CSV, show_col_types = FALSE)
  sim_episode_check(
    preps = list(
      surv = make_surv_dat_file(opt),
      mon_data = make_monthly_dat_file(opt)
    ),
    what = "10_manuscript_figures data-coverage panel",
    episode_key = episode_key_of(opt)
  )
  surv <- readRDS(make_surv_dat_file(opt))
  cov <- readRDS(make_agg_ts_file()) %>%
    dplyr::filter(country %in% unique(surv$country)) %>%
    dplyr::inner_join(enso_groups, by = c("country" = "iso3"))

  # Monthly case bubbles (non-ADM0 rows).
  mon_pts <- readRDS(make_monthly_dat_file(opt)) %>%
    dplyr::filter(country %in% unique(cov$country), admin_level != "ADM0")
  last_bubble <- mon_pts %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(last_bubble = max(as.Date(month)), .groups = "drop")

  # Country axis order: by date of the last case bubble, latest at the bottom
  # (none -> top), no group blocks, so the top rows are the countries whose
  # data stop earliest -- room for the inset.
  cov_ctry <- cov %>%
    dplyr::distinct(country, group) %>%
    dplyr::left_join(last_bubble, by = "country") %>%
    dplyr::arrange(
      dplyr::desc(dplyr::coalesce(as.numeric(last_bubble), -Inf)),
      country
    )
  ctry_levels <- unique(cov_ctry$country)
  no_name <- setdiff(ctry_levels, names(FIG1_COUNTRY_NAMES))
  if (length(no_name)) {
    stop(
      "[10] Figure 1a: no FIG1_COUNTRY_NAMES entry for ",
      paste(no_name, collapse = ", "),
      call. = FALSE
    )
  }

  ob_pts <- mon_pts %>%
    dplyr::group_by(country, month) %>%
    dplyr::summarise(tot_cases = sum(cases_month), .groups = "drop") %>%
    dplyr::inner_join(enso_groups, by = c("country" = "iso3")) %>%
    dplyr::mutate(country = factor(country, levels = ctry_levels))
  cov <- dplyr::mutate(cov, country = factor(country, levels = ctry_levels))

  # National zero-fill months: the occurrence panel's rows added from WHO /
  # national zeros (03 --national_zero_fill; the surv panel's national_zero_fill
  # column), collapsed per country to runs of consecutive months, each drawn
  # from its first month's start to its last month's end.
  nz_runs <- if (isTRUE(any(surv[["national_zero_fill"]]))) {
    surv %>%
      dplyr::filter(national_zero_fill, country %in% ctry_levels) %>%
      dplyr::distinct(
        country,
        month = lubridate::floor_date(as.Date(t_start), "month")
      ) %>%
      dplyr::arrange(country, month) %>%
      dplyr::group_by(country) %>%
      dplyr::mutate(
        m = lubridate::year(month) * 12L + lubridate::month(month),
        run = cumsum(c(1L, diff(m) != 1L))
      ) %>%
      dplyr::group_by(country, run) %>%
      dplyr::summarise(
        TL = min(month),
        TR = lubridate::ceiling_date(max(month), "month") - 1L,
        .groups = "drop"
      ) %>%
      dplyr::mutate(country = factor(country, levels = ctry_levels))
  } else {
    NULL
  }

  # White boxes: the coverage periods plus the zero-fill runs, so every
  # zero-filled month sits inside a box (all outlines drawn before all fills,
  # so overlapping boxes merge).
  boxes <- dplyr::bind_rows(
    dplyr::select(cov, country, TL, TR),
    if (!is.null(nz_runs)) dplyr::select(nz_runs, country, TL, TR)
  )

  p_cov <- ggplot(cov) +
    geom_rect(
      data = phases,
      aes(xmin = TL, xmax = TR, ymin = -Inf, ymax = Inf, fill = phase),
      alpha = PS_ENSO_PHASE_ALPHA,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = boxes,
      aes(x = TL - 5, xend = TR + 5, y = country, yend = country),
      linewidth = 3.2,
      colour = "gray40"
    ) +
    geom_segment(
      data = boxes,
      aes(x = TL, xend = TR, y = country, yend = country),
      linewidth = 3,
      colour = "white"
    ) +
    {
      if (!is.null(nz_runs)) {
        geom_segment(
          data = nz_runs,
          aes(
            x = TL,
            xend = TR,
            y = country,
            yend = country,
            linetype = "Zero-filled (World Health Organization)"
          ),
          linewidth = 1.1,
          colour = "gray60"
        )
      }
    } +
    scale_linetype_manual(
      values = c("Zero-filled (World Health Organization)" = "solid"),
      name = NULL
    ) +
    geom_point(
      data = ob_pts,
      aes(x = month, y = country, size = tot_cases, colour = group),
      alpha = 0.45
    ) +
    scale_size_continuous(
      transform = "sqrt",
      range = c(0.1, 3.5),
      breaks = c(5000, 10000, 20000),
      labels = PS_LABEL_K
    ) +
    # Window set on the COORD, not the scale: scale limits drop (not clip)
    # every coverage segment starting before the window, which cut bars short
    # of their bubbles.
    scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y",
      expand = c(0, 0)
    ) +
    scale_y_discrete(labels = function(x) {
      paste0(
        unname(FIG1_COUNTRY_NAMES[x]),
        ifelse(x %in% membership_low_isos(), MEMBERSHIP_MARK, "")
      )
    }) +
    coord_cartesian(xlim = FIG1_A_XLIM) +
    labs(x = NULL, y = NULL, size = "Cases/month") +
    scale_colour_manual(values = GROUP_COLOURS) +
    scale_fill_manual(
      values = PS_ENSO_PHASE_COL,
      labels = PS_ENSO_PHASE_LAB,
      drop = FALSE
    ) +
    guides(
      colour = "none",
      size = guide_legend(order = 1),
      fill = ps_enso_phase_guide(order = 2),
      linetype = guide_legend(order = 3, override.aes = list(linewidth = 1.1))
    ) +
    theme_ms() +
    ps_legend_theme() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = FIG1_A_TICK_PT),
      axis.text.y = element_text(size = FIG1_A_TICK_PT),
      plot.margin = margin(0.3, 0.5, 0.5, 0.5, unit = "mm"),
      # The legends stacked under the panel, one per line (the zero-fill
      # key does not fit on one line with the other two), title left of its
      # keys, at e's legend text sizes; compact keys and spacing. The
      # block centred on the panel (not the plot with its country labels),
      # its lines left-aligned.
      legend.position = "bottom",
      legend.location = "panel",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "center",
      legend.box.spacing = unit(1, "mm"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.title.position = "left",
      legend.title = element_text(size = BASE_SIZE - 1.5, vjust = 0.5),
      legend.text = element_text(
        size = BASE_SIZE - 2,
        margin = margin(l = 0.4, unit = "mm")
      ),
      legend.spacing.x = unit(2, "mm"),
      legend.spacing.y = unit(0.3, "mm"),
      legend.key.width = unit(1.8, "mm"),
      legend.key.height = unit(2.6, "mm"),
      legend.key.spacing.x = unit(PS_KEY_SPACING_MM, "mm")
    )

  add_fig1_inset(
    p_cov,
    build_fig1_inset_map(groups),
    country_levels = ctry_levels,
    last_bubble = stats::setNames(
      as.numeric(cov_ctry$last_bubble),
      cov_ctry$country
    )
  )
}

#' Figure 1: left column = observed-ONI strip over the data-coverage timeline
#' (a) with the inset map of onsets per covered month per gid (b); right column
#' = the occurrence track's ELPD lag weights over the rainfall anomaly-ONI
#' correlation with the onset-season column (c), their legends in one block,
#' then the ERF grid (d). ENSO phases (a's background) as Figure 3: NOAA
#' episodes from the phase-09 ONI artefact. The weight/correlation layout, kept
#' side by side; build_fig1() below draws the seasonality c.
build_fig1_v91 <- function(tracks = STACKED_TRACKS) {
  message(sprintf("[10] Figure 1: tracks=%s", paste(tracks, collapse = ",")))
  rel <- read_enso_relation(tracks)
  bind_track <- function(el) {
    dplyr::bind_rows(lapply(tracks, function(trk) {
      dplyr::mutate(rel[[trk]][[el]], track = trk)
    }))
  }
  erf <- bind_track("erf")
  lag_rank <- bind_track("lag_rank")
  groups <- intersect(GROUP_ORDER, unique(as.character(erf$enso_group)))
  lag_range <- range(lag_rank$single_lag, na.rm = TRUE)
  best <- fig1_best_occ_lag(lag_rank, groups)
  seas <- read_results("occ_seasonality")$seasonality

  # True minus signs on every negative tick / colourbar label (ps_true_minus()).
  p_w <- ps_true_minus(fig1_weight_panel(
    fig1_elpd_weights(lag_rank, groups),
    best,
    lag_range
  ))
  p_t <- ps_true_minus(fig1_corr_panel(groups, lag_range, best, seas))
  legends <- fig1_legend_row(p_w, p_t)
  heights <- FIG1_RIGHT_HEIGHTS
  heights[3] <- grid::convertHeight(sum(legends$heights), "mm") + unit(1, "mm") # the measured legend block
  # One flat patchwork: a free() plot inside a nested patchwork does not render
  # (patchwork 1.3.2). The strip's and the tiles' labels are freed so
  # "Rainfall month" hugs the month names; d's likewise.
  right <- patchwork::wrap_plots(
    patchwork::free(
      p_w + theme(legend.position = "none"),
      type = "label",
      side = "l"
    ),
    patchwork::free(
      p_t + theme(legend.position = "none"),
      type = "label",
      side = "l"
    ),
    patchwork::wrap_elements(full = legends),
    patchwork::free(
      ps_true_minus(fig1_erf_panel(erf)),
      type = "label",
      side = "l"
    ),
    ncol = 1,
    heights = heights
  ) +
    # No margin on the side facing a.
    patchwork::plot_annotation(
      theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 0, unit = "pt"))
    )
  right_block <- patchwork::wrap_elements(full = right) +
    theme(plot.margin = margin(0, 0, 0, 0))

  oni <- read_results("oni", sets = SETS)
  phases <- enso_phase_rects(oni, FIG1_A_XLIM)
  # ONI strip on top of panel a (tag "a" on the strip); "ONI" next to its own
  # tick labels (free label: not pushed out to the width of a's country
  # labels). Left column laid out on its own (wrap_elements): nested, patchwork
  # would align the strip's row with c's and open a gap above a. Outer margins:
  # FIG1_A_TOP_MM on top (the "a" tag level with "c"), FIG1_A_RIGHT_MM facing
  # the right column, FIG1_A_LEFT_MM on the left, 1.94 mm (5.5 pt) below.
  left <- patchwork::wrap_plots(
    patchwork::free(
      ps_true_minus(fig1_oni_strip(oni, FIG1_A_XLIM, phases)),
      type = "label",
      side = "l"
    ),
    fig1_coverage_panel(phases, groups),
    ncol = 1,
    heights = grid::unit(c(FIG1_ONI_STRIP_MM, 1), c("mm", "null"))
  ) +
    patchwork::plot_annotation(
      theme = theme(
        plot.margin = margin(
          FIG1_A_TOP_MM,
          FIG1_A_RIGHT_MM,
          1.94,
          FIG1_A_LEFT_MM,
          unit = "mm"
        )
      )
    )
  left_block <- patchwork::wrap_elements(full = left) +
    theme(plot.margin = margin(0, 0, 0, 0))

  fig <- (left_block | right_block) +
    patchwork::plot_layout(widths = FIG1_WIDTHS)
  save_figure(fig, "fig1_stacked", height_mm = FIG1_HEIGHT_MM)
}

# ---------------------------------------------------------------------------
# Figure 1c, seasonality design (the weight / correlation variant and
# build_fig1_v91() stay above). Columns = ENSO group on a 2-year calendar axis
# (x = 1..24, FIG1C_PEAK_X = Dec Y). Rows: the composite ONI cycle over the
# FIG1C_PHASES boxes; the group's mean monthly rainfall; per track
# (FIG1C_TRACKS) the ELPD lag weights at x = FIG1C_PEAK_X + lag, deployed lag
# marked, and the seasonality relative to the annual mean. Then c's legends and
# d as the last row of the same grid, so c and d share their edges.
# ---------------------------------------------------------------------------

#' c's 2-year x scale: three-letter month names at FIG1C_X_BREAKS.
fig1c_x_scale <- function(labels = FALSE) {
  scale_x_continuous(
    limits = c(0.5, 24.5),
    expand = c(0, 0),
    breaks = FIG1C_X_BREAKS,
    labels = if (labels) month.abb[(FIG1C_X_BREAKS - 1) %% 12 + 1] else NULL
  )
}

#' The peak rule at FIG1C_PEAK_X, in d's dashed reference style.
fig1c_peak_vline <- function() {
  geom_vline(
    xintercept = FIG1C_PEAK_X,
    colour = FIG1C_REF$colour,
    linewidth = FIG1C_REF$linewidth,
    linetype = FIG1C_REF$linetype
  )
}

#' A horizontal reference line in d's dashed style.
fig1c_ref_hline <- function(y) {
  geom_hline(
    yintercept = y,
    colour = FIG1C_REF$colour,
    linewidth = FIG1C_REF$linewidth,
    linetype = FIG1C_REF$linetype
  )
}

#' The phase-box edges (dotted) and the peak rule.
fig1c_phase_rules <- function() {
  list(
    geom_vline(
      xintercept = c(
        FIG1C_PHASES$from - 0.5,
        utils::tail(FIG1C_PHASES$to, 1) + 0.5
      ),
      colour = FIG1C_PHASE_RULE$colour,
      linewidth = FIG1C_PHASE_RULE$linewidth,
      linetype = FIG1C_PHASE_RULE$linetype
    ),
    fig1c_peak_vline()
  )
}

#' The theme of every c row.
fig1c_row_theme <- function(strip = FALSE, x_text = FALSE) {
  theme_ms() +
    FIG1_CD_THEME +
    FIG1_PANEL_THEME +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(FIG1C_SPACING_MM, "mm"),
      strip.text.x = if (strip) {
        element_text(
          size = BASE_SIZE - 1,
          margin = margin(b = 0.6, unit = "mm")
        )
      } else {
        element_blank()
      },
      axis.text.x = if (x_text) {
        element_text(size = FIG1C_FONT$axis_text)
      } else {
        element_blank()
      },
      axis.ticks.x = if (x_text) {
        element_line(linewidth = FIG1C_LWD$border)
      } else {
        element_blank()
      },
      axis.ticks.length.x = unit(if (x_text) 0.6 else 0, "mm"),
      axis.title.y = element_text(size = FIG1C_FONT$row_title, angle = 90),
      plot.margin = margin(0.4, 1.5, 0.4, 1.5, unit = "mm"),
      legend.position = "none",
      panel.border = element_rect(
        fill = NA,
        colour = "grey25",
        linewidth = FIG1C_LWD$border
      ),
      axis.ticks = element_line(linewidth = FIG1C_LWD$border)
    )
}

#' A panel tag at x = tag_x (npc of the plot, negative = left of it), top-aligned.
fig1c_tag_theme <- function(tag_x = FIG1C_TAG_X) {
  theme(
    plot.tag.location = "plot",
    plot.tag.position = c(tag_x, 1),
    plot.tag = element_text(
      size = BASE_SIZE + 2,
      face = "bold",
      family = FONT,
      hjust = 0,
      vjust = 1
    )
  )
}

#' Top row: the composite ONI cycle per group column (09's oni_composite), the
#' phase boxes and names, "Y" / "Y+1", under the nested "ENSO group" header.
fig1c_oni_row <- function(cyc, groups) {
  g <- factor(groups, levels = GROUP_ORDER)
  hd <- function(d) dplyr::mutate(d, header = FIG1_HEADER)
  cyc <- hd(tidyr::crossing(
    dplyr::filter(cyc, kind %in% names(FIG1C_ONI_COL)),
    enso_group = g
  ))
  boxes <- hd(tidyr::crossing(FIG1C_PHASES, enso_group = g))
  rng <- range(cyc$mean)
  ytop <- rng[2] + 0.9
  pk <- dplyr::filter(FIG1C_PHASES, phase == "Peak")
  labs_ph <- hd(tidyr::crossing(
    tibble::tibble(
      phase = c("Rising", "Peak", "Decay"),
      x = c(pk$from - 0.7, FIG1C_PEAK_X, pk$to + 0.7),
      hj = c(1, 0.5, 0),
      y = ytop - c(0.15, 0.6, 0.15)
    ),
    enso_group = g
  ))
  years <- hd(tidyr::crossing(
    tibble::tibble(x = c(6.5, 18.5), lab = FIG1C_YEAR_LABS),
    enso_group = g
  ))
  strip_margin <- margin(t = FIG1_STRIP_GAP_MM, b = 0.6, unit = "mm")
  ggplot() +
    geom_rect(
      data = boxes,
      aes(
        xmin = from - 0.5,
        xmax = to + 0.5,
        ymin = -Inf,
        ymax = Inf,
        fill = phase
      ),
      colour = NA
    ) +
    geom_text(
      data = labs_ph,
      aes(x = x, y = y, label = phase, hjust = hj),
      size = FIG1C_FONT$annot / .pt,
      vjust = 1,
      colour = FIG1C_FONT$annot_col,
      family = FONT
    ) +
    fig1c_ref_hline(0) +
    geom_text(
      data = years,
      aes(x = x, y = -Inf, label = lab),
      size = FIG1C_FONT$annot / .pt,
      vjust = -0.4,
      colour = FIG1C_FONT$annot_col,
      family = FONT
    ) +
    geom_line(
      data = cyc,
      aes(x, mean, colour = kind),
      linewidth = FIG1C_LWD$oni
    ) +
    fig1c_peak_vline() +
    ggh4x::facet_nested(
      . ~ header + enso_group,
      labeller = labeller(
        enso_group = fig1_group_labeller(),
        header = label_value
      ),
      nest_line = element_line(colour = "black", linewidth = 0.25),
      strip = ggh4x::strip_nested(
        by_layer_x = TRUE,
        text_x = ggh4x::elem_list_text(
          size = c(FIG1_HEADER_PT, BASE_SIZE - 1),
          margin = list(strip_margin, strip_margin)
        )
      )
    ) +
    scale_fill_manual(values = FIG1C_PHASE_FILL, guide = "none") +
    scale_colour_manual(
      values = FIG1C_ONI_COL,
      labels = FIG1C_ONI_LAB,
      name = NULL
    ) +
    fig1c_x_scale() +
    scale_y_continuous(
      breaks = c(-1, 0, 1),
      limits = c(rng[1] - 0.2, ytop),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Composite\nONI", tag = "c") +
    fig1c_row_theme(strip = TRUE) +
    fig1c_tag_theme()
}

#' The rainfall tile row: the group's mean monthly rainfall (mm), one scale.
fig1c_rain_row <- function(groups, rain = read_results("rain_enso")) {
  d <- dplyr::filter(
    rain$month,
    weighting == FIG1_RAIN_SEASON_WEIGHTING,
    group %in% groups
  ) %>%
    tidyr::crossing(year = 0:1) %>%
    dplyr::transmute(
      enso_group = factor(group, levels = GROUP_ORDER),
      x = year * 12L + mo,
      mm
    )
  ggplot(d, aes(x, 1, fill = mm)) +
    geom_tile(width = 1, height = 1) +
    annotate(
      "rect",
      xmin = 0.5,
      xmax = 24.5,
      ymin = 0.5,
      ymax = 1.5,
      fill = NA,
      colour = FIG1_FRAME$colour,
      linewidth = FIG1C_LWD$frame
    ) +
    fig1c_peak_vline() +
    facet_grid(. ~ enso_group) +
    scale_fill_distiller(
      palette = FIG1C_RAIN_PALETTE,
      direction = 1,
      name = FIG1C_RAIN_TITLE,
      guide = guide_colourbar(
        title.position = "left",
        title.vjust = 0.5,
        barheight = unit(1.6, "mm"),
        barwidth = unit(22, "mm")
      )
    ) +
    fig1c_x_scale() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Rain") +
    fig1c_row_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.title.y = element_text(
        size = FIG1C_FONT$row_title,
        angle = 0,
        vjust = 0.5
      ),
      legend.position = "bottom",
      legend.title = element_text(size = FIG1C_FONT$legend_title),
      legend.text = element_text(size = FIG1C_FONT$legend_text),
      legend.margin = margin(0, 0, 0, 0)
    )
}

#' ELPD-implied lag weights of every c track (fig1_elpd_weights()'s pseudo-BMA
#' over 09's lag_rank), at x = FIG1C_PEAK_X + lag; `selected` = the deployed lag.
fig1c_lag_weights <- function(lag_rank, groups, tracks = FIG1C_TRACKS) {
  lag_rank %>%
    dplyr::filter(track %in% tracks, enso_group %in% groups) %>%
    dplyr::group_by(track, enso_group) %>%
    dplyr::mutate(w = exp(elpd_diff - max(elpd_diff)), w = w / sum(w)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      enso_group = factor(enso_group, levels = GROUP_ORDER),
      track,
      lag = single_lag,
      x = FIG1C_PEAK_X + single_lag,
      w,
      selected
    )
}

#' A track's lag-weight strip: the group colour's opacity = weight on
#' FIG1_W_*'s log scale, the deployed lag a marker above; `key`: the same /
#' following-year arrows and labels in the leftmost column (first track),
#' else the strip without its headroom.
fig1c_weight_row <- function(wts, trk, key = FALSE) {
  d <- dplyr::filter(wts, track == trk)
  sel <- dplyr::filter(d, selected)
  k <- FIG1C_KEY_ARROW
  kd <- tidyr::crossing(
    enso_group = factor(
      levels(droplevels(d$enso_group))[1],
      levels = GROUP_ORDER
    ),
    tibble::tibble(
      sgn = c(-1, 1),
      lab = unname(FIG1C_KEY_TEXT[c("left", "right")])
    )
  ) %>%
    dplyr::mutate(
      x0 = FIG1C_PEAK_X + sgn * k$start,
      xend = FIG1C_PEAK_X + sgn * (k$start + k$len),
      tx = FIG1C_PEAK_X + sgn * k$gap,
      hj = (1 - sgn) / 2
    )
  if (!key) {
    kd <- kd[0, ]
  }
  ytop <- if (key) FIG1C_W_YTOP[["key"]] else FIG1C_W_YTOP[["nokey"]]
  ggplot(d, aes(x = x, y = 1)) +
    geom_tile(aes(fill = enso_group, alpha = w), width = 1, height = 1) +
    annotate(
      "rect",
      xmin = min(d$x) - 0.5,
      xmax = max(d$x) + 0.5,
      ymin = 0.5,
      ymax = 1.5,
      fill = NA,
      colour = FIG1_FRAME$colour,
      linewidth = FIG1C_LWD$frame
    ) +
    fig1c_peak_vline() +
    geom_point(
      data = sel,
      aes(x, 1.5),
      inherit.aes = FALSE,
      shape = 25,
      size = 0.9,
      stroke = 0,
      fill = "grey15",
      position = position_nudge(y = 0.17)
    ) +
    geom_segment(
      data = kd,
      aes(x = x0, xend = xend, y = k$y, yend = k$y),
      inherit.aes = FALSE,
      linewidth = k$lwd,
      colour = k$col,
      arrow = grid::arrow(length = unit(k$head_mm, "mm"), type = "closed")
    ) +
    geom_text(
      data = kd,
      aes(x = tx, y = k$text_y, label = lab, hjust = hj),
      inherit.aes = FALSE,
      size = FIG1C_FONT$annot / .pt,
      colour = k$col,
      vjust = 0,
      family = FONT,
      lineheight = FIG1C_FONT$lineheight
    ) +
    facet_grid(. ~ enso_group) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_alpha_continuous(
      range = c(0, 1),
      transform = "log10",
      limits = FIG1_W_LIM,
      oob = scales::squish,
      guide = "none"
    ) +
    fig1c_x_scale() +
    scale_y_continuous(limits = c(0.5, ytop), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Lag\nweight") +
    fig1c_row_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.title.y = element_text(
        size = FIG1C_FONT$row_title,
        angle = 0,
        vjust = (1 - 0.5) / (ytop - 0.5)
      )
    )
}

#' d's vertical "more / fewer" key (fig1_erf_panel()) in the leftmost column of
#' a seasonality row: arrows off the reference line over FIG1C_VKEY_SPAN of the
#' log axis' half-range, labels at their tips.
fig1c_vkey_layers <- function(d, trk) {
  g1 <- factor(
    levels(droplevels(d$enso_group))[1],
    levels = levels(d$enso_group)
  )
  top <- FIG1C_SEAS_LIM[[trk]][2]
  sp <- FIG1C_VKEY_SPAN
  arr <- data.frame(
    enso_group = g1,
    x = FIG1C_VKEY_X,
    y = top^c(sp[["from"]], -sp[["from"]]),
    yend = top^c(sp[["to"]], -sp[["to"]])
  )
  lab <- data.frame(
    enso_group = g1,
    x = FIG1C_VKEY_X + FIG1C_VKEY_DX,
    y = arr$yend,
    vjust = c(1, 0),
    label = unname(FIG1C_VKEY_TEXT[[trk]][c("up", "down")])
  )
  list(
    geom_segment(
      data = arr,
      aes(x = x, xend = x, y = y, yend = yend),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = unit(0.9, "mm"), type = "closed"),
      linewidth = 0.25,
      colour = "grey25"
    ),
    geom_text(
      data = lab,
      aes(x = x, y = y, label = label, vjust = vjust),
      inherit.aes = FALSE,
      hjust = 0,
      lineheight = 0.85,
      size = (BASE_SIZE - 2.5) / ggplot2::.pt,
      family = FONT,
      colour = "grey25"
    )
  )
}

#' A track's seasonality row: 09's fitted curve relative to the annual mean and
#' its interval, on the log axis, phase rules, d's vertical key on the left.
fig1c_seas_row <- function(seas, trk, x_text = FALSE) {
  d <- dplyr::filter(seas, track == trk)
  p <- ggplot(d, aes(x, rel, colour = enso_group, group = enso_group)) +
    fig1c_ref_hline(1) +
    fig1c_phase_rules()
  if (any(is.finite(d$lo))) {
    p <- p +
      geom_ribbon(
        data = dplyr::filter(d, is.finite(lo)),
        aes(ymin = lo, ymax = hi, fill = enso_group),
        colour = NA,
        alpha = FIG1C_SEAS_BAND_ALPHA
      )
  }
  p +
    geom_line(linewidth = FIG1C_LWD$seas) +
    fig1c_vkey_layers(d, trk) +
    facet_grid(. ~ enso_group) +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_y_log10(
      limits = FIG1C_SEAS_LIM[[trk]],
      breaks = FIG1C_SEAS_BREAKS[[trk]],
      labels = FIG1C_SEAS_LABELS[[trk]],
      oob = scales::oob_squish
    ) +
    fig1c_x_scale(labels = x_text) +
    labs(x = if (x_text) FIG1C_X_TITLE else NULL, y = FIG1C_SEAS_YLAB) +
    fig1c_row_theme(x_text = x_text)
}

#' A plot carrying only the lag-weight key (fig1_weight_panel()'s legend).
fig1c_weight_key_plot <- function() {
  d <- tibble::tibble(x = seq_along(FIG1_W_BREAKS), w = FIG1_W_BREAKS)
  ggplot(d, aes(x, 1, alpha = w)) +
    geom_tile(fill = "grey25") +
    scale_alpha_continuous(
      range = c(0, 1),
      transform = "log10",
      limits = FIG1_W_LIM,
      oob = scales::squish,
      breaks = FIG1_W_BREAKS,
      labels = format(FIG1_W_BREAKS, scientific = FALSE, drop0trailing = TRUE),
      name = gsub("\n", " ", FIG1_W_TITLE),
      guide = guide_legend(
        direction = "horizontal",
        nrow = 1,
        override.aes = list(fill = "grey25"),
        theme = theme(
          legend.title.position = "left",
          legend.key.size = unit(FIG1_W_KEY_MM, "mm")
        )
      )
    ) +
    theme_ms() +
    FIG1_PANEL_THEME +
    theme(
      legend.title = element_text(size = FIG1C_FONT$legend_title, vjust = 0.9),
      legend.text = element_text(size = FIG1C_FONT$legend_text)
    )
}

#' A track's name left of its block, a vertical rule over the block's rows.
fig1c_track_label <- function(label) {
  ggplot() +
    annotate(
      "segment",
      x = 1 + FIG1C_TRACK_SHIFT,
      xend = 1 + FIG1C_TRACK_SHIFT,
      y = 0.01,
      yend = 0.99,
      colour = "grey25",
      linewidth = FIG1C_LWD$border
    ) +
    annotate(
      "text",
      x = 0.35 + FIG1C_TRACK_SHIFT,
      y = 0.5,
      label = label,
      angle = 90,
      size = FIG1C_FONT$track / .pt,
      family = FONT,
      lineheight = 0.9
    ) +
    coord_cartesian(
      xlim = c(0, 1.05),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void() +
    theme(plot.margin = margin(0.4, 0, 0.4, 0, unit = "mm"))
}

#' 09's seasonality `fitted` per group and track on the 2-year axis, the
#' tracks named as 10's (FIG1C_SEAS_TRACK); an error when a group or track
#' has no curve.
fig1c_seasonality <- function(
  groups,
  art = read_results("seasonality"),
  tracks = FIG1C_TRACKS
) {
  s <- tibble::as_tibble(art$fitted) %>%
    dplyr::filter(enso_group %in% groups) %>%
    dplyr::mutate(
      track = names(FIG1C_SEAS_TRACK)[match(track, FIG1C_SEAS_TRACK)]
    ) %>%
    dplyr::filter(track %in% tracks)
  have <- dplyr::distinct(s, enso_group, track)
  need <- tidyr::crossing(enso_group = groups, track = tracks)
  gap <- dplyr::anti_join(need, have, by = c("enso_group", "track"))
  if (nrow(gap)) {
    stop(
      "[10] Figure 1c: the seasonality artefact has no fitted curve for ",
      paste(sprintf("%s:%s", gap$enso_group, gap$track), collapse = ", "),
      " (cessation needs duration fits with duration.seasonal_term).",
      call. = FALSE
    )
  }
  s %>%
    dplyr::transmute(enso_group, track, month, rel, lo, hi) %>%
    tidyr::crossing(year = 0:1) %>%
    dplyr::mutate(
      x = year * 12L + month,
      enso_group = factor(enso_group, levels = GROUP_ORDER)
    )
}

#' Panel c (and d as the grid's last row): the ONI and rain rows, one block per
#' track (lag strip over seasonality, its name on the left), c's legends, d.
#' @return patchwork
fig1c_grid <- function(lag_rank, erf, groups) {
  wts <- fig1c_lag_weights(lag_rank, groups)
  seas <- fig1c_seasonality(groups)
  h <- FIG1C_HEIGHTS
  trks <- FIG1C_TRACKS
  rows <- list(
    ps_true_minus(fig1c_oni_row(read_results("oni_composite")$cycle, groups)),
    fig1c_rain_row(groups)
  )
  for (i in seq_along(trks)) {
    rows <- c(
      rows,
      list(
        fig1c_weight_row(wts, trks[i], key = i == 1L),
        fig1c_seas_row(seas, trks[i], x_text = i == length(trks))
      )
    )
  }
  rh <- c(
    h$oni,
    h$rain,
    unlist(lapply(seq_along(trks), function(i) {
      c(if (i == 1L) h$w else h$w_nokey, h$seas)
    }))
  )
  # c's legends cut from the ONI and rain rows' gtables, side by side, the
  # lag-weight key under them; the rows themselves show none.
  box <- function(p) {
    g <- ggplot2::ggplotGrob(p + theme(legend.position = "bottom"))
    g$grobs[[which(g$layout$name == "guide-box-bottom")]]
  }
  b1 <- box(rows[[1]])
  b2 <- box(rows[[2]])
  b3 <- box(fig1c_weight_key_plot())
  leg <- gtable::gtable_row(
    "c_legends",
    list(grid::nullGrob(), b1, grid::nullGrob(), b2, grid::nullGrob()),
    height = max(sum(b1$heights), sum(b2$heights)),
    widths = grid::unit.c(
      unit(1, "null"),
      sum(b1$widths),
      unit(4, "mm"),
      sum(b2$widths),
      unit(1, "null")
    )
  )
  row3 <- gtable::gtable_row(
    "w_legend",
    list(grid::nullGrob(), b3, grid::nullGrob()),
    height = sum(b3$heights),
    widths = grid::unit.c(unit(1, "null"), sum(b3$widths), unit(1, "null"))
  )
  leg <- gtable::gtable_col(
    "c_legends2",
    list(leg, grid::nullGrob(), row3),
    width = unit(1, "null"),
    heights = grid::unit.c(
      sum(leg$heights),
      unit(FIG1_LEGEND_GAP_MM, "mm"),
      sum(row3$heights)
    )
  )
  rows[[2]] <- rows[[2]] + theme(legend.position = "none")
  # With d in the grid, the left gutter takes d's width: c's y titles are
  # freed so each hugs its own axis, and d's labels do not widen the gutter.
  rows <- c(
    lapply(rows, function(p) patchwork::free(p, type = "label", side = "l")),
    list(patchwork::wrap_elements(full = leg))
  )
  p_d <- ps_true_minus(fig1_erf_panel(erf)) +
    theme(panel.spacing.x = unit(FIG1C_SPACING_MM, "mm")) +
    fig1c_tag_theme()
  rows <- c(
    rows,
    list(patchwork::free(
      patchwork::free(p_d, type = "label", side = "l"),
      type = "space",
      side = "l"
    ))
  )
  rh <- c(rh, h$legend, h$d)
  # Column 1: the track names, each over its block's two rows; column 2: the rows.
  design <- c(
    lapply(seq_along(rows), function(r) patchwork::area(r, 2)),
    lapply(seq_along(trks), function(i) {
      patchwork::area(2 + 2 * (i - 1) + 1, 1, 2 + 2 * i, 1)
    })
  )
  labels <- lapply(trks, function(t) fig1c_track_label(FIG1C_TRACK_LABELS[[t]]))
  patchwork::wrap_plots(
    c(rows, labels),
    design = do.call(c, design),
    widths = c(FIG1C_LABEL_W, 1),
    heights = rh,
    guides = "keep"
  )
}

#' Figure 1 (seasonality c): left column as
#' build_fig1_v91() (ONI strip over a with the inset map b); right column =
#' fig1c_grid() (c with d as its last row).
build_fig1 <- function(tracks = STACKED_TRACKS) {
  message(sprintf("[10] Figure 1: tracks=%s", paste(tracks, collapse = ",")))
  rel <- read_enso_relation(tracks)
  bind_track <- function(el) {
    dplyr::bind_rows(lapply(tracks, function(trk) {
      dplyr::mutate(rel[[trk]][[el]], track = trk)
    }))
  }
  erf <- bind_track("erf")
  lag_rank <- bind_track("lag_rank")
  groups <- intersect(GROUP_ORDER, unique(as.character(erf$enso_group)))
  right <- fig1c_grid(lag_rank, erf, groups) +
    patchwork::plot_annotation(
      theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 0, unit = "pt"))
    )
  right_block <- patchwork::wrap_elements(full = right) +
    theme(plot.margin = margin(0, 0, 0, 0))

  oni <- read_results("oni", sets = SETS)
  phases <- enso_phase_rects(oni, FIG1_A_XLIM)
  left <- patchwork::wrap_plots(
    patchwork::free(
      ps_true_minus(fig1_oni_strip(oni, FIG1_A_XLIM, phases)),
      type = "label",
      side = "l"
    ),
    fig1_coverage_panel(phases, groups),
    ncol = 1,
    heights = grid::unit(c(FIG1_ONI_STRIP_MM, 1), c("mm", "null"))
  ) +
    patchwork::plot_annotation(
      theme = theme(
        plot.margin = margin(
          FIG1_A_TOP_MM,
          FIG1_A_RIGHT_MM,
          1.94,
          FIG1_A_LEFT_MM,
          unit = "mm"
        )
      )
    )
  left_block <- patchwork::wrap_elements(full = left) +
    theme(plot.margin = margin(0, 0, 0, 0))

  fig <- (left_block | right_block) +
    patchwork::plot_layout(widths = FIG1_WIDTHS)
  save_figure(fig, "fig1_stacked", height_mm = FIG1C_HEIGHT_MM)
}
