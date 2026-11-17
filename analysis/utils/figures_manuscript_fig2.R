# figures_manuscript_fig2.R — Figure 2 builders of 10_manuscript_figures.R
#   (a ADM1 IRR maps by track and ONI condition, b ENSO-sensitive population).
#
# Provides:
#   build_fig2()             the assembled figure, group-map inset right of b
#   add_fig2_header_rules(), add_fig2_row_keys(), add_fig2_bars(),
#   add_fig2_group_map(), add_fig2_tag_a()   gtable edits
#   fig2_fit_height()        height fit of the assembled gtable
#
# Details:
#   Definitions only (loads alone). The builders are 10's own: at CALL time
#   they read 10's globals -- style (FONT, BASE_SIZE, theme_ms(), TAG_THEME,
#   PS_*), groups (GROUP_*), config (RES, SCOPE, SETS, MAIN_SET, opt), the
#   figure constants and 10's helpers (read_results(), load_africa_maps(),
#   save_figure(), ...). Only 10_manuscript_figures.R sources them.
#   The FIG2_* constants stay in 10 (section 2).

#' Per-column rule under Figure 2's ONI header row: one segment per panel
#' column at the bottom of the last top-strip row, trimmed by FIG2_RULE_TRIM at
#' each end and as wide as the phase rules above it (facet_nested()'s
#' nest_line). Returns a gtable.
add_fig2_header_rules <- function(p) {
  g <- ggplotGrob(p)
  lay <- g$layout
  st <- lay[grepl("^strip-t", lay$name), ]
  panel_cols <- sort(unique(lay$l[grepl("^panel", lay$name)]))
  rule <- grid::segmentsGrob(
    x0 = FIG2_RULE_TRIM,
    x1 = grid::unit(1, "npc") - FIG2_RULE_TRIM,
    y0 = grid::unit(0, "npc"),
    y1 = grid::unit(0, "npc"),
    gp = grid::gpar(
      col = "black",
      lwd = FIG2_RULE_LWD * ggplot2::.pt,
      lineend = "butt"
    )
  )
  for (k in seq_along(panel_cols)) {
    g <- gtable::gtable_add_grob(
      g,
      rule,
      t = max(st$b),
      b = max(st$b),
      l = panel_cols[k],
      r = panel_cols[k],
      z = Inf,
      clip = "off",
      name = paste0("oni-rule-", k)
    )
  }
  g
}

#' Narrow each phase header's fill (the strips spanning several panel columns)
#' to the rule line under it: FIG2_RULE_TRIM in from both ends. Works on the
#' gtable from add_fig2_header_rules().
trim_fig2_phase_fill <- function(g) {
  for (i in which(
    grepl("^strip-t", g$layout$name) & g$layout$l != g$layout$r
  )) {
    s <- g$grobs[[i]]
    for (k in seq_along(s$grobs)) {
      ch <- s$grobs[[k]]$children
      if (is.null(ch)) {
        next
      }
      for (j in seq_along(ch)) {
        if (inherits(ch[[j]], "rect")) {
          ch[[j]]$width <- grid::unit(1, "npc") - 2 * FIG2_RULE_TRIM
        }
      }
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
add_fig2_row_keys <- function(
  g,
  track_labels,
  width_mm = FIG2_ROW_KEY_MM,
  bar_y = c(0.22, 0.78)
) {
  rows <- sort(unique(g$layout$t[grepl("^panel", g$layout$name)]))
  stopifnot(length(rows) == length(track_labels))
  pal <- scales::div_gradient_pal("#2166ac", "white", "#d73027", "Lab")
  ramp <- rev(pal(seq(0, 1, length.out = 101))) # top (red) -> bottom (blue)
  txt <- grid::gpar(
    fontsize = BASE_SIZE - 1.5,
    fontfamily = FONT,
    col = "grey25",
    lineheight = 0.85
  )
  arw <- grid::arrow(length = grid::unit(0.9, "mm"), type = "closed")
  agp <- grid::gpar(col = "grey25", fill = "grey25", lwd = 0.25 * ggplot2::.pt)
  mid <- mean(bar_y)
  gap <- 0.015
  g <- gtable::gtable_add_cols(g, grid::unit(width_mm, "mm"), pos = -1)
  for (k in seq_along(rows)) {
    lab <- FIG2_ROW_KEY_LABELS[[track_labels[k]]]
    if (is.null(lab)) {
      stop("[10] Figure 2: no row-key labels for track '", track_labels[k], "'")
    }
    key <- grid::gTree(
      children = grid::gList(
        grid::rasterGrob(
          matrix(ramp, ncol = 1),
          x = grid::unit(1.5, "mm"),
          y = grid::unit(mid, "npc"),
          width = grid::unit(1.2, "mm"),
          height = grid::unit(diff(bar_y), "npc"),
          hjust = 0,
          interpolate = TRUE
        ),
        grid::segmentsGrob(
          x0 = grid::unit(4.4, "mm"),
          x1 = grid::unit(4.4, "mm"),
          y0 = grid::unit(mid + gap, "npc"),
          y1 = grid::unit(bar_y[2], "npc"),
          arrow = arw,
          gp = agp
        ),
        grid::segmentsGrob(
          x0 = grid::unit(4.4, "mm"),
          x1 = grid::unit(4.4, "mm"),
          y0 = grid::unit(mid - gap, "npc"),
          y1 = grid::unit(bar_y[1], "npc"),
          arrow = arw,
          gp = agp
        ),
        grid::textGrob(
          lab[["up"]],
          x = grid::unit(5.4, "mm"),
          y = grid::unit(bar_y[2], "npc"),
          hjust = 0,
          vjust = 1,
          gp = txt
        ),
        grid::textGrob(
          lab[["down"]],
          x = grid::unit(5.4, "mm"),
          y = grid::unit(bar_y[1], "npc"),
          hjust = 0,
          vjust = 0,
          gp = txt
        )
      )
    )
    g <- gtable::gtable_add_grob(
      g,
      key,
      t = rows[k],
      l = ncol(g),
      clip = "off",
      name = paste0("row-key-", k)
    )
  }
  g
}

#' ONI-condition column label as Figure 2 draws it: "La Niña (−2.0)" ->
#' "ONI −2.0".
fig2_cond_of <- function(enso_label) {
  sub("^.*\\((.*)\\)$", "ONI \\1", enso_label)
}

#' Figure 2b's data: report_enso_exposure's grid at `cell` (direction, track,
#' magnitude, cutoff; results.style.enso_headline, the manuscript's one
#' definition of ENSO-sensitive), Africa + every ENSO group, one row per
#' (scope, ONI condition). The share is the artefact's own `pct`
#' (pop / pop_universe of the scope). `scope` is a factor, Africa first, then
#' the groups in GROUP_ORDER (GROUP_SHORT_FIG labels); `key` the fill key
#' (Africa or the group id).
fig2_bar_data <- function(cell = STYLE$enso_headline) {
  ex <- read_results("report_enso_exposure")
  hit <- with(
    ex$grid,
    direction == cell$direction &
      track == cell$track &
      abs(magnitude - cell$magnitude) < 1e-9 &
      abs(cutoff - cell$cutoff) < 1e-9
  )
  if (!any(hit)) {
    stop(sprintf(
      "[10] Figure 2b: report_enso_exposure has no grid cell direction=%s track=%s magnitude=%g cutoff=%g (results.style.enso_headline)",
      cell$direction,
      cell$track,
      cell$magnitude,
      cell$cutoff
    ))
  }
  d <- ex$grid[hit, ] %>%
    dplyr::filter(level %in% c("continent", "region")) %>%
    dplyr::mutate(
      key = dplyr::if_else(level == "continent", "Africa", enso_group)
    )
  grp <- intersect(GROUP_ORDER, unique(d$enso_group))
  lab <- c(
    Africa = "Africa",
    stats::setNames(unname(GROUP_SHORT_FIG[grp]), grp)
  )
  d %>%
    dplyr::mutate(
      scope = factor(lab[key], levels = unname(lab)),
      cond = factor(
        fig2_cond_of(enso_label),
        levels = fig2_cond_of(names(ex$enso_conds))
      )
    )
}

#' One ONI condition's bars (Figure 2b): Africa on top, set apart by
#' ps_africa_rows(); 0 and 100 at the panel's ends, so the panel spans its
#' map exactly. Every condition's plot has the same scales and theme, so their
#' panels are interchangeable in the gtable.
fig2_bar_plot <- function(dk, fills) {
  rows <- ps_africa_rows(dk$scope)
  dk$y <- rows$y
  ggplot(dk, aes(x = pct, y = y, fill = key)) +
    rows$rule +
    geom_col(width = FIG2_BAR_WIDTH, orientation = "y") +
    # Share at the bar's end; inside it (white) past FIG2_LABEL_INSIDE.
    geom_text(
      aes(
        label = sprintf("%.0f%%", 100 * pct),
        hjust = dplyr::if_else(pct > FIG2_LABEL_INSIDE, 1.15, -0.15),
        colour = dplyr::if_else(pct > FIG2_LABEL_INSIDE, "white", "grey25")
      ),
      size = (BASE_SIZE - 2.5) / ggplot2::.pt,
      family = FONT
    ) +
    scale_colour_identity() +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "50", "100"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = rows$breaks,
      labels = rows$labels,
      expand = expansion(add = 0.6)
    ) +
    scale_fill_manual(values = fills, guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_ms() +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "grey25", linewidth = 0.2),
      axis.ticks.x = element_line(colour = "grey25", linewidth = 0.2),
      axis.ticks.length.x = unit(0.6, "mm"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = BASE_SIZE - 2),
      axis.text.y = element_text(
        size = BASE_SIZE,
        colour = "grey15",
        hjust = 1
      ),
      plot.margin = margin(0, 0, 0, 0)
    )
}

#' Add Figure 2b to the maps' gtable `g` (after the header-rule / row-key
#' edits): three rows under the legend row (spacer with tag "b" + title, bars,
#' bar axis); in each map column the panel and axis of that condition's bar
#' plot, spanning the column (= the map) exactly; the group labels (a bar
#' plot's y axis) right-aligned against the first column, in a column added
#' left of the track-strip column for whatever they need beyond it (so the
#' track titles stay next to the maps).
add_fig2_bars <- function(g, bars, cond_levels) {
  lay <- g$layout
  pan <- lay[grepl("^panel", lay$name), ]
  panel_cols <- sort(unique(pan$l))
  stopifnot(length(panel_cols) == length(cond_levels))
  gb <- lay[lay$name == "guide-box-bottom", ]
  if (nrow(gb) != 1L) {
    stop("[10] Figure 2b: no bottom legend in the maps' gtable")
  }
  strip_l <- lay[grepl("^strip-l", lay$name), ]

  fills <- c(Africa = PS_AFRICA_COL, GROUP_COLOURS)
  grobs <- lapply(cond_levels, function(cd) {
    ggplotGrob(fig2_bar_plot(dplyr::filter(bars, cond == cd), fills))
  })
  part <- function(gg, re) gg$grobs[[which(grepl(re, gg$layout$name))[1]]]

  g <- gtable::gtable_add_rows(
    g,
    grid::unit(c(FIG2_BAR_GAP_MM, FIG2_BAR_ROW_MM, FIG2_BAR_AXIS_MM), "mm"),
    pos = gb$b
  )
  r_gap <- gb$b + 1L
  r_bar <- gb$b + 2L
  r_ax <- gb$b + 3L
  for (k in seq_along(panel_cols)) {
    g <- gtable::gtable_add_grob(
      g,
      part(grobs[[k]], "^panel"),
      t = r_bar,
      l = panel_cols[k],
      clip = "off",
      name = paste0("bar-panel-", k)
    )
    g <- gtable::gtable_add_grob(
      g,
      part(grobs[[k]], "^axis-b"),
      t = r_ax,
      l = panel_cols[k],
      clip = "off",
      name = paste0("bar-axis-", k)
    )
  }
  ax_l <- part(grobs[[1]], "^axis-l")
  lab_w <- grid::convertWidth(grid::grobWidth(ax_l), "mm", valueOnly = TRUE)
  sc <- min(strip_l$l)
  have <- sum(grid::convertWidth(
    g$widths[sc:(panel_cols[1] - 1L)],
    "mm",
    valueOnly = TRUE
  ))
  g <- gtable::gtable_add_cols(
    g,
    grid::unit(max(0, lab_w + 0.5 - have), "mm"),
    pos = sc - 1L
  )
  panel_cols <- panel_cols + 1L
  lc <- sc
  g <- gtable::gtable_add_grob(
    g,
    grid::editGrob(
      ax_l,
      vp = grid::viewport(
        x = 1,
        just = "right",
        width = grid::unit(lab_w, "mm")
      )
    ),
    t = r_bar,
    l = lc,
    r = panel_cols[1] - 1L,
    clip = "off",
    name = "bar-labels"
  )
  # Tag "b" + title, left-aligned at the content edge, on the spacer row.
  tag <- grid::textGrob(
    "b",
    x = 0,
    y = 0.5,
    hjust = 0,
    vjust = 0.5,
    gp = PS_TAG_GP
  )
  tag_w <- grid::convertWidth(grid::grobWidth(tag), "mm", valueOnly = TRUE) +
    1.5
  tit <- grid::textGrob(
    FIG2_BAR_TITLE,
    x = grid::unit(tag_w, "mm"),
    y = 0.5,
    hjust = 0,
    vjust = 0.5,
    gp = grid::gpar(fontsize = BASE_SIZE, fontfamily = FONT)
  )
  gtable::gtable_add_grob(
    g,
    grid::gTree(children = grid::gList(tit, tag)),
    t = r_gap,
    l = lc,
    r = max(panel_cols),
    clip = "off",
    name = "bar-title"
  )
}

#' The group-map inset: fig. GM's layers (build_figS_group_map(),
#' utils/figures_supp.R: countries filled by ENSO group, PS_GREY not_modelled
#' outside the partition, white country borders) with no legend, caption,
#' short labels or membership dagger, zero margins. coord_sf() keeps the map's
#' aspect, so the grob fits its gtable cell whatever the cell's shape.
#' `members`: the fit registry's partition (09 modelled_units$members).
fig2_group_map_inset <- function(maps, members) {
  grey <- ps_country_grey_layer(
    ps_country_status(maps$afro_sf, members$iso3),
    legend_order = 2,
    labels = PS_GREY_LAB_WRAP,
    show_legend = FALSE
  )
  p <- build_figS_group_map(
    maps$afro_sf,
    members,
    GROUP_COLOURS,
    GROUP_LABELS[GROUP_ORDER],
    GROUP_SHORT_FIG[GROUP_ORDER],
    sf::st_bbox(maps$afro_sf),
    low_isos = character(0),
    base_layers = grey,
    overlay_layers = list(geom_sf(
      data = maps$afro_sf,
      fill = NA,
      colour = "white",
      linewidth = FIG2_INSET_BORDER_LWD
    )),
    map_theme = taxdat::map_theme() + ps_map_font(),
    base_size = BASE_SIZE
  )
  # Drop the short-label layer (geom_sf_label; panel b's row labels carry the names).
  p$layers <- p$layers[
    !vapply(p$layers, function(l) inherits(l$geom, "GeomLabel"), logical(1))
  ]
  p +
    labs(caption = NULL) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

#' Add the group-map inset to Figure 2's gtable `g` after add_fig2_bars(): in
#' the row-key column (no column added, so panel a keeps its width), spanning
#' the legend row down to the bar row; the map's bottom FIG2_INSET_RULE_GAP_MM
#' above panel b's Africa rule (its height from ps_africa_rows()'s geometry:
#' `n_scope` rows, PS_AFRICA_GAP_ROWS, the 0.6-row axis expansion, over
#' FIG2_BAR_ROW_MM), right-aligned FIG2_INSET_SHIFT_MM in from the column's
#' right edge, FIG2_INSET_GAP_MM narrower than the column.
add_fig2_group_map <- function(g, inset, n_scope) {
  lay <- g$layout
  bar <- lay[grepl("^bar-panel-", lay$name), ]
  if (!nrow(bar)) {
    stop(
      "[10] Figure 2 inset: panel b not in the gtable (add_fig2_bars() first)"
    )
  }
  key <- lay[grepl("^row-key-", lay$name), ]
  if (!nrow(key)) {
    stop("[10] Figure 2 inset: no row-key column in the gtable")
  }
  gb <- lay[lay$name == "guide-box-bottom", ]
  if (nrow(gb) != 1L) {
    stop("[10] Figure 2 inset: no bottom legend row in the gtable")
  }
  gap <- PS_AFRICA_GAP_ROWS
  top <- n_scope + gap + 0.6
  bot <- 1 - 0.6
  rule_y <- n_scope - 0.5 + gap / 2
  rule_mm <- FIG2_BAR_ROW_MM * (rule_y - bot) / (top - bot)
  y0 <- grid::unit(rule_mm + FIG2_INSET_RULE_GAP_MM, "mm")
  vp <- grid::viewport(
    x = grid::unit(1, "npc") - grid::unit(FIG2_INSET_SHIFT_MM, "mm"),
    y = y0,
    just = c("right", "bottom"),
    width = grid::unit(1, "npc") - grid::unit(FIG2_INSET_GAP_MM, "mm"),
    height = grid::unit(1, "npc") - y0
  )
  gtable::gtable_add_grob(
    g,
    grid::gTree(children = grid::gList(ggplotGrob(inset)), vp = vp),
    t = gb$t,
    b = max(bar$b),
    l = min(key$l),
    r = ncol(g),
    clip = "off",
    name = "group-map"
  )
}

#' Tag "a" at the maps' top-left: the leftmost content column (right of the
#' plot margin), the top header row.
add_fig2_tag_a <- function(g) {
  st <- g$layout[grepl("^strip-t", g$layout$name), ]
  gtable::gtable_add_grob(
    g,
    grid::textGrob("a", x = 0, y = 1, hjust = 0, vjust = 1, gp = PS_TAG_GP),
    t = min(st$t),
    l = 2L,
    clip = "off",
    name = "tag-a"
  )
}

#' Figure 2's height (mm) leaving no space around the maps at `width_mm`: the
#' gtable's fixed rows plus its null (map) rows at the size the maps' fixed
#' aspect gives them once the fixed columns are taken off the width. `make_g()`
#' is built on a scratch PNG device of the figure's width (PNG font metrics;
#' the cairo PDF may differ by a fraction of a mm), + 0.5 mm.
fig2_fit_height <- function(make_g, width_mm = WIDTH_MM, probe_mm = 200) {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = width_mm, height = probe_mm, units = "mm", res = DPI)
  on.exit(
    {
      grDevices::dev.off()
      unlink(f)
    },
    add = TRUE
  )
  g <- make_g()
  if (!isTRUE(g$respect)) {
    stop("[10] Figure 2: the maps' gtable has no fixed aspect")
  }
  nh <- grid::unitType(g$heights) == "null"
  nw <- grid::unitType(g$widths) == "null"
  fixed_h <- sum(grid::convertHeight(g$heights[!nh], "mm", valueOnly = TRUE))
  fixed_w <- sum(grid::convertWidth(g$widths[!nw], "mm", valueOnly = TRUE))
  per_null <- (width_mm - fixed_w) / sum(as.numeric(g$widths[nw]))
  fixed_h + per_null * sum(as.numeric(g$heights[nh])) + 0.5
}

#' Figure 2: (a) ADM1 IRR maps faceted track (rows, strips on the left) x ONI
#' condition (columns = results.derivation.enso_conds, under an ENSO-phase
#' header filled in Figure 1's phase colours, with booktabs-style rules);
#' (b) under each ONI column, the share of the population in ENSO-sensitive
#' ADM1s (fig2_bar_data(), results.style.enso_headline) for Africa and each
#' ENSO group, the panel spanning its map; the group-map
#' inset (fig2_group_map_inset(), FIG2_INSET_*) right of b and the phase header
#' row unfilled under FIG2_PHASE_HEADER_FILL FALSE.
#' Basemap layers as taxdat::output_plot_map() (country fill -> ADM1 fill ->
#' rivers -> lakes -> country borders), built here so the two country greys
#' (PS_GREY; modelled = >= 1 onset, phase-09 `onset_rate`, as Figure 1b) are
#' keyed in the legend. IRR capped to FIG2_MAP_LO/HI; colourbar and grey keys
#' in one row between the maps and panel b, a direction key right of each map
#' row; scale bar + north arrow on the top-left map. Duration IRR is inverted
#' here so "high = worse" holds on every row.
build_fig2 <- function(tracks = STACKED_TRACKS) {
  message(sprintf(
    "[10] Figure 2 (ADM1 maps, stacked): tracks=%s",
    paste(tracks, collapse = ",")
  ))
  maps <- load_africa_maps("Figure 2")
  basemap <- africa_basemap(maps$afro_sf) # utils/figures_cb.R (lakes + rivers, st_intersects clip)

  rel <- read_enso_relation(tracks)
  d <- dplyr::bind_rows(lapply(tracks, function(trk) {
    dplyr::mutate(rel[[trk]]$adm1_irr, track = trk)
  })) %>%
    dplyr::mutate(
      IRR = dplyr::if_else(track == "duration", 1 / IRR, IRR),
      log_ratio = log10(pmin(pmax(IRR, FIG2_MAP_LO), FIG2_MAP_HI))
    )
  # enso_label "La Niña (−2.0)" -> phase "La Niña" (header row) + cond "ONI −2.0".
  cond_levels <- unique(d$enso_label)
  d <- d %>%
    dplyr::mutate(
      phase = factor(
        sub(" \\(.*$", "", enso_label),
        levels = unique(sub(" \\(.*$", "", cond_levels))
      ),
      cond = factor(
        fig2_cond_of(enso_label),
        levels = fig2_cond_of(cond_levels)
      ),
      track = track_factor(track)
    )

  map_sf <- maps$adm1_sf %>% dplyr::inner_join(d, by = c("GID_1" = "adm1_gid"))
  sig_sf <- dplyr::filter(map_sf, significant)

  # Scale bar + north arrow on the top-left panel only.
  tl <- data.frame(
    track = factor(levels(droplevels(d$track))[1], levels = levels(d$track)),
    phase = factor(levels(d$phase)[1], levels = levels(d$phase)),
    cond = factor(levels(d$cond)[1], levels = levels(d$cond))
  )

  # Header strips: ggh4x orders them outer layer first -- the phase headers in
  # level order (filled), then one ONI-condition strip per column (unfilled).
  n_ph <- nlevels(d$phase)
  n_cd <- nlevels(d$cond)
  ph_fill <- unname(FIG2_PHASE_FILL[dplyr::case_when(
    grepl("Niña", levels(d$phase)) ~ "nina",
    grepl("Niño", levels(d$phase)) ~ "elnino",
    TRUE ~ "neutral"
  )])
  if (!isTRUE(FIG2_PHASE_HEADER_FILL)) {
    ph_fill <- rep(NA, n_ph)
  }

  p <- ggplot() +
    ps_country_grey_layer(
      ps_country_status(maps$afro_sf, modelled_countries()),
      legend_order = 2
    ) +
    geom_sf(
      data = map_sf,
      aes(fill = log_ratio),
      colour = "white",
      linewidth = FIG2_ADM1_LWD
    ) +
    map_overlay_layers(basemap, maps$afro_sf) +
    geom_sf(
      data = sig_sf,
      fill = NA,
      colour = "black",
      linewidth = FIG2_SIG_LWD
    ) +
    # Countries whose class is uncertain (p_max < style.membership_low).
    geom_sf(
      data = dplyr::filter(maps$afro_sf, ISO3_CODE %in% membership_low_isos()),
      fill = NA,
      colour = "black",
      linetype = "22",
      linewidth = 0.25
    ) +
    # Scale bar 1000 km; arrow centred over it.
    map_scale_layers(data = tl) +
    ps_ratio_fill_scale("Effect size") +
    ggh4x::facet_nested(
      track ~ phase + cond,
      switch = "y",
      nest_line = element_line(colour = "black", linewidth = FIG2_RULE_LWD),
      resect = FIG2_RULE_TRIM,
      strip = ggh4x::strip_nested(
        by_layer_x = FALSE,
        # colour/linewidth per strip: elem_list_rect() does not
        # inherit the theme's blank strip.background border.
        background_x = ggh4x::elem_list_rect(
          fill = c(ph_fill, rep(NA, n_cd)),
          colour = rep(NA, n_ph + n_cd),
          linewidth = rep(0, n_ph + n_cd)
        ),
        text_x = ggh4x::elem_list_text(
          size = c(rep(BASE_SIZE + 1, n_ph), rep(BASE_SIZE - 0.5, n_cd))
        )
      )
    ) +
    taxdat::map_theme() +
    ps_map_font() +
    theme(
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      text = element_text(size = BASE_SIZE - 1),
      strip.text = element_text(size = BASE_SIZE + 2),
      strip.text.y.left = element_text(
        size = BASE_SIZE,
        angle = 90,
        margin = margin(r = FIG2_STRIP_GAP_MM, unit = "mm")
      ),
      legend.justification = c(0, 0),
      legend.background = element_blank()
    ) +
    ps_legend_theme() +
    # The maps are too small to host the legends: one row under the maps
    # (colourbar + the two greys), top-aligned; panel b goes under it.
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "top",
      legend.spacing.x = unit(6, "mm"),
      legend.title = element_text(size = BASE_SIZE - 1, vjust = 0.8),
      legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(3, "mm")
    ) +
    theme(
      panel.spacing.x = unit(FIG2_SPACING_MM[1], "mm"),
      panel.spacing.y = unit(FIG2_SPACING_MM[2], "mm")
    ) +
    coord_sf(expand = FALSE) +
    theme(
      plot.margin = margin(
        FIG2_MARGIN_MM,
        FIG2_MARGIN_MM,
        FIG2_MARGIN_MM,
        FIG2_MARGIN_MM,
        "mm"
      )
    ) +
    theme(
      legend.box.spacing = unit(FIG2_LEGEND_GAP_MM, "mm"),
      legend.margin = margin(0, 0, 0, 0)
    )

  bars <- fig2_bar_data()
  miss <- setdiff(levels(d$cond), levels(bars$cond))
  if (length(miss)) {
    stop(
      "[10] Figure 2b: report_enso_exposure lacks the map condition(s) ",
      paste(miss, collapse = ", ")
    )
  }

  inset <- if (!is.null(FIG2_INSET_SHIFT_MM)) {
    fig2_group_map_inset(maps, read_results("modelled_units")$members)
  } # the fit registry's partition (09)

  make_g <- function() {
    g <- add_fig2_header_rules(p)
    g <- trim_fig2_phase_fill(g)
    g <- add_fig2_row_keys(g, levels(droplevels(d$track)))
    g <- add_fig2_bars(g, bars, levels(d$cond))
    if (!is.null(inset)) {
      g <- add_fig2_group_map(g, inset, n_scope = nlevels(bars$scope))
    }
    add_fig2_tag_a(g)
  }
  # make_g() is passed unevaluated: the gtable edits call ggplotGrob(), whose
  # layout depends on the device's font metrics, so they must run inside
  # save_figure() with the PNG/PDF device open (building the grob first gives
  # a different layout).
  save_figure(
    make_g(),
    "fig2_adm1_stacked",
    height_mm = fig2_fit_height(make_g)
  )
}
