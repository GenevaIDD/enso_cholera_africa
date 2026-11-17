# figures_supp.R — ggplot builders for the supplementary figures, drawn by
#   10_manuscript_figures.R from phase-09 artefacts only.
#
# Provides:
#   build_figS5_convergence(), build_figS6_lag_loo()   diagnostics, lag LOO
#   build_figS_<name>()     one builder per supplementary figure (forecast,
#                           smoothing, attribution, events, membership,
#                           units, groups, data coverage, cascade, hindcast,
#                           analogues, seasonality, WHO ratio)
#   attrib_ratio_map_sf(), attrib_ratio_map_layers()   ADM1 ratio map shared
#                           by fig. AE (a) and Figure 4b
#
# Details:
#   Pure builders: return ggplot / patchwork objects, never save. Colours and
#   the theme are passed in by the caller (10's GROUP_COLOURS / theme_ms()).
#   Packages via `::`: ggplot2, dplyr, tidyr, patchwork; sf, scales for maps.

# Figure -> builder input (phase-09 artefact):
#   fig. S5  sampler diagnostics of the deployed models   (convergence, report_convergence)
#   fig. S6  LOO lag rankings by track and ENSO group     (report_lags; Pareto-k flags from lag_averaging)
#   forecast vs counterfactual by country                  (attribution$country_ribbon)
#   historical reconstruction (smoothing) by country       (smoothing$by_country)
#   fig. AE  population in ADM1s with more ENSO-attributable cases
#            (forecast vs counterfactual)                  (report_attrib_exposure)
#   the forecast / counterfactual ADM1 ratio map shared by fig. AE (a) and
#            Figure 4b (attrib_ratio_map_sf(), attrib_ratio_map_layers())
#   fig. EV  attributable fraction of past strong ENSO events
#            (simulated + analytic)                        (report_event_attrib, report_event_af_model)
#   fig. MB  membership uncertainty of the deployed classes (report_membership)
#   hindcast CRPSS against climatology by group and lead   (hindcast_occ$skill, hindcast_cases$month_skill)
#   fig. MU / GM / KC  modelled units, deployed partition, kl classes by K (below)
#   fig. DC  outbreak data coverage                     (report_dataset$outbreaks)
#   fig. PC  processing cascade as a flow diagram       (report_dataset$cascade)
#   fig. SV  historical reconstruction vs observations  (smoothing$validation)
#   fig. FA  factual vs counterfactual by analogue set  (attribution ribbons)
#   fig. AA  ADM1 attribution maps by analogue set      (report_attrib_exposure$adm1)
#   fig. SE  rainfall vs onset, cessation and case seasonality by group (seasonality)

#' Fig. S5: (a) per-fit maximum R-hat and minimum bulk / tail ESS by track and
#' ENSO group, with the DIAG_THRESHOLDS lines; (b) distribution of the
#' per-parameter R-hat and bulk ESS by track.
#' @param scalars report_convergence$fits (one row per deployed fit)
#' @param summary bind_rows of the convergence artefacts' $summary
#' @param thresholds named c(rhat, ess)
#' @param group_cols named colours by enso_group; @param group_lab labels
#' @param theme a ggplot theme
build_figS5_convergence <- function(
  scalars,
  summary,
  thresholds,
  group_cols,
  group_lab,
  theme
) {
  trk_lab <- c(
    occurrence = "Occurrence",
    duration = "Duration",
    monthly = "Monthly cases"
  )
  sc <- scalars |>
    dplyr::transmute(
      track = factor(trk_lab[track], trk_lab),
      enso_group,
      `max R-hat` = max_rhat,
      `min bulk ESS` = min_ess_bulk,
      `min tail ESS` = min_ess_tail
    ) |>
    tidyr::pivot_longer(-c(track, enso_group), names_to = "stat") |>
    dplyr::mutate(
      stat = factor(stat, c("max R-hat", "min bulk ESS", "min tail ESS")),
      # dodge and legend in the display group order (group_lab's names)
      enso_group = factor(enso_group, names(group_lab))
    )
  thr <- tibble::tibble(
    stat = factor(
      c("max R-hat", "min bulk ESS", "min tail ESS"),
      levels(sc$stat)
    ),
    y = c(thresholds[["rhat"]], thresholds[["ess"]], thresholds[["ess"]])
  )
  pa <- ggplot2::ggplot(sc, ggplot2::aes(track, value, colour = enso_group)) +
    ggplot2::geom_hline(
      data = thr,
      ggplot2::aes(yintercept = y),
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 1.6
    ) +
    ggplot2::facet_wrap(~stat, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(
      values = group_cols,
      breaks = names(group_lab),
      labels = group_lab,
      name = NULL
    ) +
    ggplot2::labs(x = NULL, y = NULL, tag = "a") +
    theme +
    ggplot2::theme(legend.position = "bottom")
  ps <- summary |>
    dplyr::transmute(track = factor(trk_lab[track], trk_lab), rhat, ess_bulk) |>
    tidyr::pivot_longer(c(rhat, ess_bulk), names_to = "stat") |>
    dplyr::mutate(
      stat = factor(
        ifelse(stat == "rhat", "R-hat", "bulk ESS"),
        c("R-hat", "bulk ESS")
      )
    )
  thr2 <- tibble::tibble(
    stat = factor(c("R-hat", "bulk ESS"), c("R-hat", "bulk ESS")),
    x = c(thresholds[["rhat"]], thresholds[["ess"]])
  )
  pb <- ggplot2::ggplot(
    dplyr::filter(ps, is.finite(value)),
    ggplot2::aes(value)
  ) +
    ggplot2::geom_histogram(bins = 40, fill = "grey55", colour = NA) +
    ggplot2::geom_vline(
      data = thr2,
      ggplot2::aes(xintercept = x),
      linetype = "dashed",
      colour = "grey20",
      linewidth = 0.3
    ) +
    ggplot2::facet_grid(track ~ stat, scales = "free") +
    ggplot2::labs(x = "Per-parameter value", y = "Parameters", tag = "b") +
    theme
  patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1, 1.4))
}

#' Fig. S6: expected log predictive density relative to the best lag
#' (+-2 SE of the paired difference) at every lag -11..11, by track (rows)
#' and ENSO group (columns); the deployed lag filled, the rank-1 lag circled.
#' @param ranking report_lags$ranking; @param selected report_lags$selected
#'
#' When `selected` carries them: the lags below the shortest lag allowed (floor_lag:
#' the min_lag rule,
#' -8, and for duration / monthly also >= occurrence lag - occ_floor_offset)
#' shaded darker, with a dashed rule at it, and the rank-2 lag marked by an
#' open diamond. `k_flags` (optional; lag_averaging $lags rows: track,
#' enso_group, single_lag, n_k_gt07) adds an open triangle under every lag with
#' Pareto k > 0.7 on some rows; without it the ranking's own n_k_gt07 is used
#' when present.
build_figS6_lag_loo <- function(
  ranking,
  selected,
  group_cols,
  group_lab,
  theme,
  k_flags = NULL
) {
  trk_lab <- c(
    occurrence = "Occurrence",
    duration = "Duration",
    monthly = "Monthly cases"
  )
  keep <- intersect(
    c(
      "track",
      "enso_group",
      "deployed_lag",
      "rank1_lag",
      "rank2_lag",
      "floor_lag"
    ),
    names(selected)
  )
  d <- ranking |>
    dplyr::select(-dplyr::any_of(setdiff(keep, c("track", "enso_group")))) |>
    dplyr::left_join(
      dplyr::select(selected, dplyr::all_of(keep)),
      by = c("track", "enso_group")
    ) |>
    dplyr::mutate(
      track = factor(trk_lab[track], trk_lab),
      group = factor(group_lab[enso_group], group_lab[names(group_lab)]),
      deployed = single_lag == deployed_lag,
      rank1 = single_lag == rank1_lag
    )
  flo <- if ("floor_lag" %in% names(d)) {
    dplyr::distinct(
      dplyr::filter(d, !is.na(floor_lag)),
      track,
      group,
      floor_lag
    )
  } else {
    NULL
  }
  r2 <- if ("rank2_lag" %in% names(d)) {
    dplyr::filter(d, single_lag == rank2_lag)
  } else {
    d[0, ]
  }
  kf <- if (!is.null(k_flags) && nrow(k_flags)) {
    kk <- dplyr::filter(k_flags, n_k_gt07 > 0) |>
      dplyr::transmute(
        track = factor(trk_lab[track], trk_lab),
        enso_group,
        single_lag = as.integer(single_lag)
      )
    dplyr::semi_join(d, kk, by = c("track", "enso_group", "single_lag"))
  } else if ("n_k_gt07" %in% names(d)) {
    dplyr::filter(d, !is.na(n_k_gt07) & n_k_gt07 > 0)
  } else {
    d[0, ]
  }
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(single_lag, elpd_diff, colour = enso_group)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = -0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey93",
      colour = NA
    )
  if (!is.null(flo) && nrow(flo)) {
    p <- p +
      ggplot2::geom_rect(
        data = flo,
        ggplot2::aes(
          xmin = -Inf,
          xmax = floor_lag - 0.5,
          ymin = -Inf,
          ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "grey80",
        colour = NA
      ) +
      ggplot2::geom_vline(
        data = flo,
        ggplot2::aes(xintercept = floor_lag - 0.5),
        linetype = "dashed",
        colour = "grey30",
        linewidth = 0.3
      )
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = elpd_diff - 2 * se_diff,
        ymax = elpd_diff + 2 * se_diff
      ),
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(d, !deployed),
      shape = 1,
      size = 1.1
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(d, rank1),
      shape = 1,
      size = 3.4,
      stroke = 0.45,
      colour = "black"
    )
  if (nrow(r2)) {
    p <- p +
      ggplot2::geom_point(
        data = r2,
        shape = 5,
        size = 2.6,
        stroke = 0.4,
        colour = "black"
      )
  }
  # the deployed lag last, large and black-outlined, so it reads at print size
  p <- p +
    ggplot2::geom_point(
      data = dplyr::filter(d, deployed),
      ggplot2::aes(fill = enso_group),
      shape = 21,
      size = 2.3,
      stroke = 0.5,
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = group_cols, guide = "none")
  if (nrow(kf)) {
    p <- p +
      ggplot2::geom_point(
        data = kf,
        ggplot2::aes(y = elpd_diff - 2 * se_diff),
        shape = 2,
        size = 1.1,
        stroke = 0.35,
        colour = "firebrick"
      )
  }
  p +
    ggplot2::facet_grid(track ~ group, scales = "free_y") +
    ggplot2::scale_colour_manual(values = group_cols, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(-10, 10, 5)) +
    ggplot2::labs(
      x = "Lag (months; negative = ENSO leads cholera)",
      y = "ELPD difference to the best lag (± 2 SE)",
      caption = paste0(
        "Large filled point: deployed lag; circled: rank 1; diamond: rank 2",
        if (nrow(kf)) "; red triangle: Pareto k > 0.7 on some rows" else "",
        "\nLight shading: leads; dark shading left of the dashed rule: ",
        "lags shorter than the shortest lag allowed"
      )
    ) +
    theme
}

#' Facet order and labels of a by-country figure: countries sorted by ENSO
#' group (in `group_order`) then ISO3, labelled "ISO3 (group short label)".
#' @noRd
.country_facets <- function(country, group, group_order, group_lab) {
  ord <- order(match(group, group_order), country)
  lev <- unique(country[ord])
  lab <- stats::setNames(
    sprintf("%s (%s)", lev, group_lab[group[match(lev, country)]]),
    lev
  )
  list(levels = lev, labels = lab)
}

#' y limits of a free-y count panel: from 0, and 0-1 when the panel is all zero.
#' @noRd
.y_floor0 <- function(l) c(0, if (is.finite(l[2]) && l[2] > 0) l[2] else 1)

#' Integer-only breaks of a count axis (at most ~4 per panel).
#' @noRd
.int_breaks <- function(l) {
  b <- scales::breaks_extended(4)(l)
  b <- unique(round(b[b >= 0]))
  if (length(b) < 2L) c(0, max(1, ceiling(l[2]))) else b
}

#' Forecast vs neutral-ENSO counterfactual by country: monthly median with 50%
#' and 95% bands per arm (the 07d_simulator_figures_counterfactual.R fig6 /
#' fig7 series), each country to its firm horizon. Forecast in its ENSO group's
#' colour, counterfactual grey.
#' @param country_ribbon attribution$country_ribbon
#' @param metric "cases" | "onsets"; @param aset counterfactual analog set
#' @param vlines dates marked dashed (the attribution start / factual launch)
build_figS_forecast_country <- function(
  country_ribbon,
  metric,
  aset,
  group_order,
  group_cols,
  group_lab,
  cf_col,
  theme,
  vlines = NULL,
  ncol = 5L
) {
  q <- function(p) paste0(metric, "_", p)
  d <- dplyr::filter(country_ribbon, run == "forecast" | analog_set == aset)
  # Rows past the country's firm horizon (the Africa window's analogue
  # continuation, `continuation` = TRUE; e.g. GW and GE after Feb 2027) are not shown.
  if ("continuation" %in% names(d)) {
    d <- dplyr::filter(d, !continuation)
  }
  d <- dplyr::mutate(
    d,
    run = factor(as.character(run), c("forecast", "counterfactual")),
    lo95 = .data[[q("q025")]],
    lo50 = .data[[q("q25")]],
    med = .data[[q("q50")]],
    hi50 = .data[[q("q75")]],
    hi95 = .data[[q("q975")]]
  )
  fc <- .country_facets(d$country, d$enso_group, group_order, group_lab)
  d$country <- factor(d$country, fc$levels)
  fcst <- dplyr::filter(d, run == "forecast")
  cf <- dplyr::filter(d, run == "counterfactual")
  p <- ggplot2::ggplot(mapping = ggplot2::aes(month)) +
    ggplot2::geom_ribbon(
      data = cf,
      ggplot2::aes(ymin = lo95, ymax = hi95),
      fill = cf_col,
      alpha = 0.15
    ) +
    ggplot2::geom_ribbon(
      data = cf,
      ggplot2::aes(ymin = lo50, ymax = hi50),
      fill = cf_col,
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = cf,
      ggplot2::aes(y = med),
      colour = cf_col,
      linetype = "22",
      linewidth = 0.35
    ) +
    ggplot2::geom_ribbon(
      data = fcst,
      ggplot2::aes(ymin = lo95, ymax = hi95, fill = enso_group),
      alpha = 0.18
    ) +
    ggplot2::geom_ribbon(
      data = fcst,
      ggplot2::aes(ymin = lo50, ymax = hi50, fill = enso_group),
      alpha = 0.35
    ) +
    ggplot2::geom_line(
      data = fcst,
      ggplot2::aes(y = med, colour = enso_group),
      linewidth = 0.4
    )
  if (length(vlines)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = as.Date(vlines),
        linetype = "dashed",
        colour = "grey45",
        linewidth = 0.25
      )
  }
  p +
    ggplot2::facet_wrap(
      ~country,
      scales = "free_y",
      ncol = ncol,
      labeller = ggplot2::as_labeller(fc$labels)
    ) +
    # legend keys in the display group order (not the ids' alphabetical order)
    ggplot2::scale_fill_manual(
      values = group_cols,
      breaks = group_order,
      labels = group_lab,
      name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = group_cols,
      breaks = group_order,
      labels = group_lab,
      name = NULL
    ) +
    ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
    # y floored at 0; an all-zero panel (LBR, TGO, NAM) gets 0-1, not -0.05-0.05;
    # integer onset ticks; short-scale labels
    ggplot2::scale_y_continuous(
      limits = .y_floor0,
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      breaks = if (metric == "onsets") .int_breaks else ggplot2::waiver(),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (metric == "cases") {
        "Suspected cases / month"
      } else {
        "Outbreak onsets / month"
      }
    ) +
    theme +
    ggplot2::theme(legend.position = "bottom")
}

#' Historical reconstruction by country: the smoothing's monthly median with
#' 50% / 95% bands restricted to surveillance-covered unit-months (the
#' country's ENSO-group colour, `group_cols`),
#' its unrestricted median (grey, dashed), and the observed series (points,
#' non-zero months) --
#' the 08c_smoothing_compile.R fig3 country panels.
#' @param by_country smoothing$by_country; @param country_group smoothing$country_group
build_figS_smoothing_country <- function(
  by_country,
  country_group,
  metric,
  group_order,
  group_lab,
  group_cols,
  theme,
  start = NULL,
  ncol = 5L
) {
  keep <- by_country[vapply(
    by_country,
    function(x) identical(x$metric, metric),
    logical(1)
  )]
  rb <- dplyr::bind_rows(lapply(keep, function(x) {
    dplyr::mutate(x$ribbon_masked, country = x$scope)
  }))
  um <- dplyr::bind_rows(lapply(keep, function(x) {
    dplyr::mutate(
      dplyr::select(x$ribbon_unmasked, t_start, q50),
      country = x$scope
    )
  }))
  ob <- dplyr::bind_rows(lapply(keep, function(x) {
    dplyr::mutate(x$obs, country = x$scope)
  }))
  if (!is.null(start)) {
    rb <- dplyr::filter(rb, t_start >= start)
    um <- dplyr::filter(um, t_start >= start)
    ob <- dplyr::filter(ob, t_start >= start)
  }
  ctry <- unique(c(rb$country, um$country))
  fc <- .country_facets(
    ctry,
    unname(country_group[ctry]),
    group_order,
    group_lab
  )
  fac <- function(d) {
    dplyr::mutate(
      d,
      country = factor(country, fc$levels),
      enso_group = unname(country_group[as.character(country)])
    )
  }
  ggplot2::ggplot(mapping = ggplot2::aes(t_start)) +
    ggplot2::geom_ribbon(
      data = fac(rb),
      ggplot2::aes(ymin = q025, ymax = q975, fill = enso_group),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      data = fac(rb),
      ggplot2::aes(ymin = q25, ymax = q75, fill = enso_group),
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = fac(um),
      ggplot2::aes(y = q50),
      colour = "grey45",
      linetype = "22",
      linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = fac(rb),
      ggplot2::aes(y = q50, colour = enso_group),
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
    ggplot2::scale_colour_manual(values = group_cols, guide = "none") +
    # observed months with a non-zero value (zeros are in the artefact; drawing
    # them would hide the series)
    ggplot2::geom_point(
      data = fac(dplyr::filter(ob, value > 0)),
      ggplot2::aes(y = value),
      shape = 1,
      size = 0.6,
      stroke = 0.3
    ) +
    ggplot2::facet_wrap(
      ~country,
      scales = "free_y",
      ncol = ncol,
      labeller = ggplot2::as_labeller(fc$labels)
    ) +
    # every 5 years (3-yearly labels ran together at 5 columns); integer
    # onset ticks; y floored at 0; short-scale labels
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::scale_y_continuous(
      limits = .y_floor0,
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      breaks = if (metric == "onsets") .int_breaks else ggplot2::waiver(),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (metric == "cases") {
        "Suspected cases / month"
      } else {
        "Outbreak onsets / month"
      }
    ) +
    theme
}

#' The ADM1 rows of the forecast / counterfactual ratio map (fig. AE (a),
#' Figure 4b): report_attrib_exposure$adm1 at the headline slice (attribution
#' window, explosive replicates dropped) for one metric, analog set and IPCC
#' magnitude, ADM1s with a ratio only, joined onto their polygons.
#' @param window the headlined attribution cell (results.style.attrib_headline.window)
#' @return sf (adm1_sf rows) with ratio and class (robust | conflicting | none)
attrib_ratio_map_sf <- function(
  adm1,
  metric,
  aset,
  magnitude,
  adm1_sf,
  window = "attribution",
  guard = TRUE
) {
  d <- adm1
  win_hl <- window
  if (all(c("window", "explosive_guard") %in% names(d))) {
    d <- dplyr::filter(d, window == win_hl, explosive_guard == guard)
  }
  d <- d %>%
    dplyr::filter(
      metric == !!metric,
      analog_set == aset,
      abs(magnitude - !!magnitude) < 1e-9,
      !is.na(ratio)
    )
  dplyr::inner_join(adm1_sf, d, by = c("GID_1" = "adm1_unit"))
}

#' The ratio map's own layers: the ADM1 ratio of means as fill (log10, capped
#' at `limits`; an infinite ratio is drawn at the cap) on `fill_scale`, then the
#' IPCC classes at `magnitude` as hatching (unhatched = robust; diagonal = no
#' robust signal; crossed = conflicting) on a second fill scale, or, without
#' pattern support (`hatch_fills` NULL), the robust ADM1s outlined.
#' @param fill_scale the ratio's fill scale (on log10 values)
#' @param hatch_fills list(none =, conflicting =) pattern fills, or NULL
#' @param hatch_ncol,hatch_nrow layout of the hatch legend
#' @return a list of layers / scales (base and overlay layers are the caller's)
attrib_ratio_map_layers <- function(
  map_sf,
  limits,
  fill_scale,
  hatch_fills,
  magnitude,
  hatch_ncol = NULL,
  hatch_nrow = NULL
) {
  out <- list(
    ggplot2::geom_sf(
      data = map_sf,
      ggplot2::aes(fill = log10(pmin(pmax(ratio, limits[1]), limits[2]))),
      colour = "white",
      linewidth = 0.005
    ),
    fill_scale
  )
  cls_lab <- c(none = "No robust signal", conflicting = "Conflicting signal")
  rob_title <- sprintf("IPCC, m = %g%%\n(unhatched = robust)", 100 * magnitude)
  if (!is.null(hatch_fills)) {
    # Simplified outlines for the pattern layer only (rendering ~1000 patterned
    # full-resolution polygons is slow); the fill layer keeps full detail.
    hd <- dplyr::filter(map_sf, class %in% names(cls_lab)) %>%
      dplyr::mutate(hatch = factor(cls_lab[class], levels = cls_lab)) %>%
      sf::st_simplify(preserveTopology = TRUE, dTolerance = 2000)
    c(
      out,
      list(
        ggnewscale::new_scale_fill(),
        ggplot2::geom_sf(data = hd, ggplot2::aes(fill = hatch), colour = NA),
        ggplot2::scale_fill_manual(
          values = stats::setNames(
            list(hatch_fills$none, hatch_fills$conflicting),
            cls_lab
          ),
          name = rob_title,
          drop = FALSE,
          guide = ggplot2::guide_legend(
            order = 2,
            ncol = hatch_ncol,
            nrow = hatch_nrow
          )
        )
      )
    )
  } else {
    c(
      out,
      list(
        ggplot2::geom_sf(
          data = dplyr::filter(map_sf, class == "robust"),
          ggplot2::aes(colour = "Robust increase"),
          fill = NA,
          linewidth = 0.2
        ),
        ggplot2::scale_colour_manual(
          values = c(`Robust increase` = "black"),
          name = NULL,
          guide = ggplot2::guide_legend(order = 2)
        )
      )
    )
  }
}

#' Fig. AE: population in ADM1s with more ENSO-attributable cases (forecast vs
#' counterfactual, report_attrib_exposure), for a double-column page.
#' (a) ADM1 map of the mean forecast / mean counterfactual ratio (log scale,
#' capped at `limits`), with the IPCC robustness classes at the headline
#' magnitude as hatching; an infinite ratio is drawn at the cap, 0 vs 0 is left
#' grey. Legends sit inside the map, bottom left.
#' (b) Share of each scope's population by IPCC class (stacked bars) and the
#' Childs shares with 95% bootstrap intervals (ratio of means filled, ratio of
#' medians open), one column per magnitude m, on the headline slice of a v2+
#' artefact. Returns a patchwork (a | b); the caller adds panel tags.
#' @param adm1,exposure report_attrib_exposure$adm1 / $exposure
#' @param metric "cases" | "onsets"; @param aset counterfactual analog set
#' @param magnitude the headline m (style.attrib_headline) the map's classes use
#' @param adm1_sf ADM1 polygons keyed GID_1
#' @param bbox map extent (sf bbox, e.g. the island-free Africa outline); its
#'   north is cropped to the modelled ADM1s plus `pad_deg` (the unmodelled
#'   Sahara and Maghreb would otherwise take a fifth of the map's height) and
#'   its west widened by `west_pad_deg` of ocean for the legends
#' @param base_layers drawn first (country greys + ggnewscale::new_scale_fill())
#' @param overlay_layers drawn over the fill (rivers, lakes, borders, scale bar)
#' @param limits ratio cap c(lo, hi); @param scope_order Africa then regions
#'   (a scope missing from it is appended, never dropped)
#' @param class_cols named colours c(robust, conflicting, childs)
#' @param hatch_fills list(none =, conflicting =) pattern fills, or NULL (no
#'   pattern support: robust ADM1s outlined instead)
#' @param map_theme theme for (a); @param theme theme for (b)
#' @param base_size base font size (pt) of the legends
#' @param widths relative widths of (a) and (b)
build_figS_attrib_exposure <- function(
  adm1,
  exposure,
  metric,
  aset,
  magnitude,
  adm1_sf,
  bbox,
  base_layers,
  overlay_layers,
  limits,
  scope_order,
  class_cols,
  hatch_fills,
  map_theme,
  theme,
  base_size = 7,
  widths = c(1.3, 1),
  pad_deg = 2,
  west_pad_deg = 6,
  window = "attribution",
  guard = TRUE
) {
  win_hl <- window
  headline <- function(d) {
    if (all(c("window", "explosive_guard") %in% names(d))) {
      dplyr::filter(d, window == win_hl, explosive_guard == guard)
    } else {
      d
    }
  }
  pct_lab <- function(m) sprintf("m = %g%%", 100 * m)
  metric_lab <- if (metric == "cases") "cases" else "outbreak onsets"
  cf_lab <- if (aset == "all") {
    "counterfactual"
  } else {
    sprintf("%s counterfactual", aset)
  }
  leg_text <- ggplot2::theme(
    legend.title = ggplot2::element_text(size = base_size - 1),
    legend.text = ggplot2::element_text(size = base_size - 1.5)
  )

  # --- (a) map -------------------------------------------------------------
  map_sf <- attrib_ratio_map_sf(
    adm1,
    metric,
    aset,
    magnitude,
    adm1_sf,
    window = win_hl,
    guard = guard
  )
  bbox[["ymax"]] <- min(bbox[["ymax"]], sf::st_bbox(map_sf)[["ymax"]] + pad_deg)
  bbox[["xmin"]] <- bbox[["xmin"]] - west_pad_deg
  all_brk <- c(1 / 3, 1 / 2, 1, 2, 3)
  all_lab <- c("1/3", "1/2", "1", "2", "3")
  keep <- all_brk >= limits[1] - 1e-9 & all_brk <= limits[2] + 1e-9
  fill_scale <- ggplot2::scale_fill_gradient2(
    midpoint = 0,
    low = "#2166ac",
    mid = "white",
    high = "#d73027",
    limits = log10(limits),
    oob = scales::squish,
    breaks = log10(all_brk[keep]),
    labels = all_lab[keep],
    name = sprintf(
      "Forecast ÷ %s\n(mean %s)",
      sub(" ", "\n", cf_lab),
      metric_lab
    ),
    guide = ggplot2::guide_colourbar(
      order = 1,
      direction = "horizontal",
      title.position = "top",
      theme = ggplot2::theme(
        legend.key.width = grid::unit(16, "mm"),
        legend.key.height = grid::unit(1.8, "mm")
      )
    )
  )
  pa <- ggplot2::ggplot() +
    base_layers +
    attrib_ratio_map_layers(
      map_sf,
      limits,
      fill_scale,
      hatch_fills,
      magnitude,
      hatch_ncol = 1L
    ) +
    overlay_layers +
    ggplot2::coord_sf(
      xlim = bbox[c("xmin", "xmax")],
      ylim = bbox[c("ymin", "ymax")],
      expand = FALSE
    ) +
    map_theme +
    leg_text +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 0, 0, unit = "mm"),
      legend.position = "inside",
      legend.position.inside = c(0.005, 0.01),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.title = ggplot2::element_text(
        size = base_size - 1.5,
        lineheight = 0.9
      ),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.spacing.y = grid::unit(0.2, "mm"),
      legend.box.spacing = grid::unit(0, "mm"),
      legend.key.spacing.y = grid::unit(0, "mm"),
      legend.background = ggplot2::element_blank(),
      legend.key.height = grid::unit(2, "mm"),
      legend.key.width = grid::unit(3, "mm")
    )

  # --- (b) population shares --------------------------------------------------
  e <- headline(exposure) %>%
    dplyr::filter(metric == !!metric, analog_set == aset, level != "country")
  lev <- unique(c(scope_order, sort(unique(as.character(e$scope)))))
  lev <- lev[lev %in% e$scope]
  e <- dplyr::mutate(
    e,
    y = match(as.character(scope), rev(lev)),
    m_lab = factor(
      pct_lab(magnitude),
      levels = pct_lab(sort(unique(magnitude)))
    )
  )
  bars <- e %>%
    dplyr::filter(class %in% c("robust", "conflicting")) %>%
    dplyr::mutate(class = factor(class, levels = c("conflicting", "robust")))
  ch_lab <- c(childs = "means ≥ 1 + m", childs_median = "medians ≥ 1 + m")
  ch <- dplyr::filter(e, class %in% names(ch_lab)) %>%
    dplyr::mutate(
      est = factor(ch_lab[class], levels = ch_lab),
      yy = y + ifelse(class == "childs", 0.13, -0.13)
    )
  pb <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = bars,
      ggplot2::aes(x = pct, y = y, fill = class),
      orientation = "y",
      width = 0.62
    ) +
    ggplot2::geom_linerange(
      data = ch,
      ggplot2::aes(y = yy, xmin = pct_lo, xmax = pct_hi),
      colour = class_cols[["childs"]],
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = ch,
      ggplot2::aes(x = pct, y = yy, shape = est),
      colour = class_cols[["childs"]],
      fill = "white",
      size = 1,
      stroke = 0.4
    ) +
    ggplot2::facet_grid(~m_lab) +
    ggplot2::scale_fill_manual(
      values = class_cols[c("robust", "conflicting")],
      breaks = c("robust", "conflicting"),
      labels = c(robust = "Robust", conflicting = "Conflicting"),
      name = "IPCC",
      guide = ggplot2::guide_legend(order = 1, nrow = 1)
    ) +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(16, 21), ch_lab),
      name = "Childs et al. rule (95% bootstrap)",
      guide = ggplot2::guide_legend(order = 2, nrow = 1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(lev),
      labels = rev(lev),
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) format(100 * x),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.06))
    ) +
    ggplot2::labs(
      x = sprintf("Share of the simulated population (%%, %s)", metric_lab),
      y = NULL
    ) +
    theme +
    leg_text +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      legend.spacing.y = grid::unit(0.3, "mm"),
      legend.key.height = grid::unit(2.2, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      panel.spacing.x = grid::unit(3, "mm"),
      plot.margin = ggplot2::margin(1, 1.5, 0, 1, unit = "mm")
    )
  patchwork::wrap_plots(pa, pb, nrow = 1, widths = widths)
}


#' Fig. EV: attributable fraction of the past strong ENSO events, one panel per
#' event, Africa and each ENSO group on the y axis, AF (%) with its 95% interval
#' for three readings: simulated with the free factual run, simulated with the
#' smoothing as factual (report_event_attrib, explosive replicates dropped), and
#' the analytic crossbasis-only AF on the observed rows (report_event_af_model).
#' @param ea   report_event_attrib artefact; @param eam report_event_af_model artefact
#' @param metric "cases" | "onsets"
#' @param scope_order display order of the scopes ("Africa", group short names)
#' @param theme a ggplot theme
#' @param tag panel tag (e.g. "a") or NULL. What each interval carries is
#'   stated in the figure caption (supplement fig. "Attribution of past ENSO
#'   events"), not in the plot.
build_figS_event_attrib <- function(
  ea,
  eam,
  metric,
  scope_order,
  theme,
  base_size = 7,
  tag = NULL,
  guard = TRUE
) {
  lab_ev <- stats::setNames(
    sprintf(
      "%s\n(launch %s)",
      event_display_label(ea$events$label),
      format(ea$events$launch, "%b %Y")
    ),
    ea$events$event
  )
  READ <- c(
    simulation = "Simulated (factual run)",
    smoothing = "Simulated (reconstruction as factual)",
    analytic = "Analytic (ENSO exposure\u2013response, observed months)"
  )
  sim <- ea$af_scope %>%
    dplyr::filter(
      explosive_guard == guard,
      level != "country",
      metric == !!metric
    ) %>%
    dplyr::transmute(
      event,
      scope = ifelse(level == "continent", "Africa", region),
      reading = factual_source,
      af = af_of_means,
      lo = af_lo,
      hi = af_hi
    )
  ana <- eam$scope %>%
    dplyr::filter(level != "country", metric == !!metric) %>%
    dplyr::transmute(
      event,
      scope = ifelse(level == "continent", "Africa", region),
      reading = "analytic",
      af = af_mean,
      lo = af_lo,
      hi = af_hi
    )
  # (an empty analytic table -- sub-windows without report_event_af_model's --
  # is skipped: its zero-length ifelse() column would not bind)
  d <- dplyr::bind_rows(sim, if (nrow(ana)) ana) %>%
    dplyr::filter(!is.na(af)) %>%
    dplyr::mutate(
      event = factor(lab_ev[event], levels = lab_ev),
      scope = factor(scope, levels = rev(scope_order)),
      reading = factor(READ[reading], levels = READ),
      dplyr::across(c(af, lo, hi), ~ 100 * pmax(pmin(.x, 1), -1))
    )
  ggplot2::ggplot(d, ggplot2::aes(x = af, y = scope, colour = reading)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi),
      width = 0,
      orientation = "y",
      position = ggplot2::position_dodge(width = 0.6),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6),
      size = 1
    ) +
    ggplot2::facet_wrap(~event, nrow = 1) +
    ggplot2::scale_colour_manual(
      values = c("grey15", "grey62", "#b2182b"),
      name = NULL
    ) +
    ggplot2::labs(
      x = sprintf(
        "ENSO-attributable fraction of %s (%%, 95%% interval; clipped to \u00b1100%%)",
        metric
      ),
      y = NULL,
      tag = tag
    ) +
    theme +
    ggplot2::theme(legend.position = "bottom")
}

#' Display form of a stored event label: "Nino" / "Nina" with the tilde, an
#' en dash in a year range ("El Nino 2015-16" -> "El Niño 2015–16").
event_display_label <- function(x) {
  x <- gsub(
    "Nino",
    "Ni\u00f1o",
    gsub("Nina", "Ni\u00f1a", x, fixed = TRUE),
    fixed = TRUE
  )
  gsub("([0-9]{4})-([0-9]{2,4})", "\\1\u2013\\2", x)
}

#' The per-winter sub-windows of report_event_attrib / report_event_af_model
#' shaped as the artefacts build_figS_event_attrib()
#' reads: one pseudo-event per (event, sub-window), id "<event>:<sub-window>",
#' label "<event label>, DJF yyyy/yy". NULL when the artefacts have no
#' sub-windows.
event_sub_window_view <- function(ea, eam) {
  if (is.null(ea$sub_windows)) {
    return(NULL)
  }
  w <- dplyr::distinct(
    ea$sub_windows$windows,
    event,
    sub_window,
    k,
    peak_season_year
  ) %>%
    dplyr::left_join(
      dplyr::select(ea$events, event, label, launch),
      by = "event"
    ) %>%
    dplyr::arrange(launch, k) %>%
    dplyr::mutate(
      id = paste(event, sub_window, sep = ":"),
      label = sprintf(
        "%s, DJF %d/%02d",
        label,
        peak_season_year,
        (peak_season_year + 1L) %% 100L
      )
    )
  re_id <- function(d) {
    if (is.null(d)) {
      NULL
    } else {
      dplyr::mutate(d, event = paste(event, sub_window, sep = ":"))
    }
  }
  list(
    ea = list(
      events = dplyr::transmute(w, event = id, label, launch),
      af_scope = re_id(ea$sub_windows$af_scope)
    ),
    eam = list(
      scope = if (!is.null(eam$sub_windows)) {
        re_id(eam$sub_windows$scope)
      } else {
        eam$scope[0, , drop = FALSE]
      }
    )
  )
}

#' Display form of a class-sweep cell id: "ns4_K5" -> "df 4, K = 5" (the ENSO
#' basis df and the number of classes), "_mix" -> ", mixture prior".
cell_display_label <- function(id) {
  out <- sub("^ns([0-9]+)_K([0-9]+)", "df \\1, K = \\2", id)
  sub("_mix$", ", mixture prior", out)
}

#' Fig. MB: membership uncertainty of the deployed classes (report_membership).
#' (A) Responsibility matrix of the deployed sweep cell
#' (agreeing chains): countries (rows, grouped by deployed class) x classes,
#' fill = posterior class probability, p_max printed, countries under `low`
#' marked. (B) Alluvial-style class paths across the sweep cells: each column
#' is a cell (the deployed cell first, then by elpd difference to it), each
#' country a line through its MAP class matched to the deployed labels
#' (kl_match_classes()), coloured by its deployed class; a country leaving its
#' deployed class is drawn as a hollow point. Agreeing chains. (ggalluvial is
#' not on the cluster's R stack, so the flows are drawn as segments.)
#' @param m report_membership artefact; @param group_cols / group_lab by group id
#' @param low style.membership_low
build_figS_membership <- function(
  m,
  group_cols,
  group_lab,
  theme,
  low = 0.6,
  base_size = 7
) {
  patchwork::wrap_plots(
    membership_heatmap(m, group_lab, theme, low = low, base_size = base_size),
    membership_paths(m, group_cols, group_lab, theme, base_size = base_size),
    ncol = 2,
    widths = c(1, 2.2)
  )
}

#' Panel a of fig. MB: membership probability per country (rows, by deployed
#' group then ISO3) and class, the deployed class labelled with p_max; a
#' dagger after the ISO3 below `low`.
membership_heatmap <- function(
  m,
  group_lab,
  theme,
  low = 0.6,
  base_size = 7,
  tag = "a"
) {
  ord <- names(group_lab)
  ct <- m$country[!m$country$carried_over, , drop = FALSE]
  rc <- grep("^r_", names(ct), value = TRUE)
  long <- tidyr::pivot_longer(
    ct[, c("iso3", "group", "p_max", rc)],
    dplyr::all_of(rc),
    names_to = "class",
    values_to = "r"
  ) |>
    dplyr::mutate(class = sub("^r_", "", class))
  # rows by deployed group (display order) then ISO3, as fig. KC and the country ERF figure
  lev_c <- ct$iso3[order(match(ct$group, ord), ct$iso3)]
  long <- dplyr::mutate(
    long,
    iso3 = factor(iso3, rev(lev_c)),
    class = factor(class, intersect(ord, unique(class)))
  )
  ylab <- stats::setNames(
    ifelse(ct$p_max < low, paste0(ct$iso3, " \u2020"), ct$iso3),
    ct$iso3
  )
  ggplot2::ggplot(long, ggplot2::aes(class, iso3, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(
      data = dplyr::filter(long, as.character(class) == group),
      ggplot2::aes(label = sprintf("%.2f", p_max)),
      size = (base_size - 2) / ggplot2::.pt
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "grey25",
      limits = c(0, 1),
      name = "P(class)"
    ) +
    ggplot2::scale_x_discrete(labels = group_lab, position = "top") +
    ggplot2::scale_y_discrete(labels = ylab) +
    # the numbers (probability of the deployed class) and the dagger (below
    # `low`) are explained in the figure caption, not in the plot
    ggplot2::labs(x = NULL, y = NULL, tag = tag) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        title.position = "top",
        barwidth = grid::unit(25, "mm")
      )
    ) +
    theme +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )
}

#' Panel b of fig. MB: each country's class across the sweep cells (the
#' deployed cell first, then by elpd difference to it).
membership_paths <- function(
  m,
  group_cols,
  group_lab,
  theme,
  base_size = 7,
  tag = "b",
  ref_label = "deployed"
) {
  ord <- names(group_lab)
  cp <- dplyr::filter(m$cell_partition, chains == m$chains)
  sv <- dplyr::filter(m$sensitivity, chains == m$chains)
  by_elpd <- sv$id[order(-sv$elpd_diff)]
  cells <- c(m$cell, by_elpd[by_elpd != m$cell])
  cp <- cp |>
    dplyr::mutate(
      matched = ifelse(is.na(matched), "unmatched", matched),
      x = match(id, cells)
    ) |>
    dplyr::filter(!is.na(x)) |>
    dplyr::arrange(
      x,
      match(matched, c(ord, "unmatched")),
      match(deployed, ord),
      iso3
    ) |>
    dplyr::group_by(x) |>
    dplyr::mutate(y = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(moved = matched != deployed)
  seg <- cp |>
    dplyr::arrange(iso3, x) |>
    dplyr::group_by(iso3) |>
    dplyr::mutate(xend = dplyr::lead(x), yend = dplyr::lead(y)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(xend))
  xl <- vapply(
    cells,
    function(id) {
      s <- sv[sv$id == id, ]
      if (id == m$cell) {
        sprintf("%s\n%s", cell_display_label(id), ref_label)
      } else {
        sprintf(
          "%s\n%.0f \u00b1 %.0f",
          cell_display_label(id),
          s$elpd_diff,
          s$se_diff
        )
      }
    },
    character(1)
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = deployed),
      alpha = 0.6,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      data = cp,
      ggplot2::aes(x, y, colour = deployed, shape = moved),
      size = 1.1,
      stroke = 0.4
    ) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 1),
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = group_cols,
      labels = group_lab,
      breaks = names(group_lab),
      name = "Deployed class"
    ) +
    # ISO3 at both ends of each path (deployed cell left, last cell right)
    ggplot2::geom_text(
      data = dplyr::filter(cp, x == 1),
      ggplot2::aes(x - 0.12, y, label = iso3, colour = deployed),
      hjust = 1,
      size = (base_size - 2.5) / ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cp, x == max(x)),
      ggplot2::aes(x + 0.12, y, label = iso3, colour = deployed),
      hjust = 0,
      size = (base_size - 2.5) / ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(cells),
      labels = xl,
      expand = ggplot2::expansion(add = 0.7)
    ) +
    ggplot2::scale_y_reverse(breaks = NULL) +
    ggplot2::labs(
      x = sprintf(
        "Class model (ELPD difference to the %s one \u00b1 SE)",
        ref_label
      ),
      y = sprintf(
        "Countries, stacked by most likely class (matched to the %s classes)",
        ref_label
      ),
      tag = tag
    ) +
    theme +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = base_size - 2)
    )
}

# ---------------------------------------------------------------------------
# Supplement maps of the analysis set and the grouping
#   fig. MU  figS_modelled_units   modelled units by admin level + ADM1s of the
#                                  modelled countries without a recorded episode
#                                  (09 modelled_units, GADM polygons)
#   fig. GM  figS_group_map        the deployed partition, centroid labels
#   fig. KC  figS_kl_classes_by_K  MAP class and its probability per country,
#                                  K = 2..5 at the pinned basis (report_kl_sweep)
# ---------------------------------------------------------------------------

# Display labels / colours of the unit status codes 09 writes
# (results_derive.R::modelled_units_tables()).
MODELLED_UNIT_STATUS <- c(
  adm1 = "ADM1 unit",
  adm2 = "ADM2 unit",
  no_onset = "Modelled unit, no onset in the panel",
  no_episode = "ADM1 without a recorded episode"
)
MODELLED_UNIT_COL <- c(
  adm1 = "#3b6ea8",
  adm2 = "#9ecae1",
  no_onset = "#fdd49e",
  no_episode = "#f2f2f2"
)
MODELLED_UNIT_BORDER <- "grey30"

#' Fig. MU: modelled units by admin level. `unit_sf` (gid, status, geometry),
#' `gap_sf` (adm1, status, geometry), `afro_sf` country outlines; base /
#' overlay layers from 10 (country greys, borders).
build_figS_modelled_units <- function(
  unit_sf,
  gap_sf,
  afro_sf,
  bbox,
  base_layers = NULL,
  overlay_layers = NULL,
  map_theme = ggplot2::theme_void(),
  base_size = 7
) {
  lev <- names(MODELLED_UNIT_STATUS)
  unit_sf$status <- factor(unit_sf$status, lev)
  gap_sf$status <- factor(gap_sf$status, lev)
  # ADM1 tiles under ADM2 tiles: an ADM2 unit sits inside its ADM1 parent.
  unit_sf <- unit_sf[order(unit_sf$status != "adm1"), ]
  ggplot2::ggplot() +
    base_layers +
    ggplot2::geom_sf(
      data = gap_sf,
      ggplot2::aes(fill = status),
      colour = MODELLED_UNIT_BORDER,
      linewidth = 0.05
    ) +
    ggplot2::geom_sf(
      data = unit_sf,
      ggplot2::aes(fill = status),
      colour = MODELLED_UNIT_BORDER,
      linewidth = 0.05
    ) +
    ggplot2::scale_fill_manual(
      values = MODELLED_UNIT_COL,
      labels = MODELLED_UNIT_STATUS,
      drop = TRUE,
      name = "Analysis unit"
    ) +
    overlay_layers +
    ggplot2::coord_sf(
      xlim = bbox[c("xmin", "xmax")],
      ylim = bbox[c("ymin", "ymax")],
      expand = FALSE,
      datum = NA
    ) +
    map_theme +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.02),
      legend.justification = c(0, 0),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size)
    )
}

#' Fig. GM: the deployed partition. `members` (iso3, group id), countries of
#' `afro_sf` (ISO3_CODE) filled by group; each group's display short label at
#' a point on the union of its members; `low_isos` (p_max below
#' style.membership_low) carry `mark` (Figure 1a's dagger) at a point on their
#' surface, with no outline (as Figure 1b).
build_figS_group_map <- function(
  afro_sf,
  members,
  group_cols,
  group_lab,
  group_short,
  bbox,
  low_isos = character(0),
  low = NULL,
  mark = "\u2020",
  base_layers = NULL,
  overlay_layers = NULL,
  map_theme = ggplot2::theme_void(),
  base_size = 7
) {
  ord <- names(group_lab)
  m <- dplyr::inner_join(
    afro_sf,
    dplyr::select(members, iso3, group),
    by = c("ISO3_CODE" = "iso3")
  )
  m$group <- factor(as.character(m$group), ord)
  lab <- do.call(
    rbind,
    lapply(split(m, as.character(m$group)), function(x) {
      p <- suppressWarnings(sf::st_point_on_surface(sf::st_union(sf::st_make_valid(
        x
      ))))
      sf::st_sf(group = as.character(x$group[1]), geometry = p)
    })
  )
  lab$label <- unname(group_short[as.character(lab$group)])
  lo <- m[m$ISO3_CODE %in% low_isos, ]
  lo_pt <- if (nrow(lo)) {
    sf::st_sf(
      label = mark,
      geometry = suppressWarnings(
        sf::st_point_on_surface(sf::st_geometry(sf::st_make_valid(lo)))
      )
    )
  } else {
    NULL
  }
  ggplot2::ggplot() +
    base_layers +
    ggplot2::geom_sf(
      data = m,
      ggplot2::aes(fill = group),
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = group_cols,
      labels = group_lab,
      breaks = ord,
      name = "ENSO group"
    ) +
    overlay_layers +
    {
      if (!is.null(lo_pt)) {
        ggplot2::geom_sf_text(
          data = lo_pt,
          ggplot2::aes(label = label),
          size = (base_size + 2) / ggplot2::.pt,
          colour = "black"
        )
      }
    } +
    ggplot2::geom_sf_label(
      data = lab,
      ggplot2::aes(label = label),
      size = base_size / ggplot2::.pt,
      label.padding = grid::unit(0.6, "mm"),
      fill = "white"
    ) +
    ggplot2::coord_sf(
      xlim = bbox[c("xmin", "xmax")],
      ylim = bbox[c("ymin", "ymax")],
      expand = FALSE,
      datum = NA
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      caption = if (length(low_isos)) {
        sprintf(
          "%s probability of the deployed group below %s",
          mark,
          if (is.null(low)) "the threshold" else sprintf("%.2f", low)
        )
      } else {
        NULL
      }
    ) +
    map_theme +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.02),
      legend.justification = c(0, 0),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size)
    )
}

#' The report_kl_sweep partition rows of the K = `Ks` fits at the pinned
#' cell's basis and prior (`cell`, e.g. "ns5_K5" or "ns5_K4_mix"), with K.
kl_classes_by_K_frame <- function(partition, fits, cell, Ks = 2:5) {
  basis <- sub("_K[0-9]+.*$", "", cell)
  mix <- grepl("_mix$", cell)
  f <- fits[
    fits$basis == basis & as.logical(fits$mix) == mix & fits$K %in% Ks,
    c("id", "K"),
    drop = FALSE
  ]
  if (!nrow(f)) {
    stop(
      "kl_classes_by_K_frame: no sweep fit at basis ",
      basis,
      if (mix) " (mix)",
      " with K in ",
      paste(Ks, collapse = ","),
      ".",
      call. = FALSE
    )
  }
  dplyr::inner_join(partition, f, by = "id")
}

#' Fig. KC: MAP class (fill) and its probability (text, opacity) per country
#' (rows, ordered by deployed group then ISO3) and K (columns). Classes named
#' like a deployed group (short "east" -> id "enso_east") take its colour and
#' display label; other classes (klc<k>) are grey.
build_figS_kl_classes_by_K <- function(
  d,
  deployed,
  group_cols,
  group_lab,
  theme,
  low = NULL,
  base_size = 7
) {
  ord <- names(group_lab)
  d$class_id <- ifelse(
    paste0("enso_", d$group) %in% ord,
    paste0("enso_", d$group),
    as.character(d$group)
  )
  other <- sort(setdiff(unique(d$class_id), ord))
  cols <- c(
    group_cols[ord],
    stats::setNames(
      grDevices::grey.colors(max(1L, length(other)), 0.45, 0.8)[seq_along(
        other
      )],
      other
    )
  )
  labs <- c(group_lab[ord], stats::setNames(other, other))
  dep <- deployed[as.character(d$iso3)]
  lev <- unique(d$iso3[order(match(dep, ord), d$iso3)])
  d$iso3 <- factor(d$iso3, rev(lev))
  d$class_id <- factor(d$class_id, c(ord, other))
  d$K <- factor(paste0("K = ", d$K), paste0("K = ", sort(unique(d$K))))
  ggplot2::ggplot(d, ggplot2::aes(K, iso3)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = class_id, alpha = p),
      colour = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", p)),
      size = (base_size - 2) / ggplot2::.pt
    ) +
    ggplot2::scale_fill_manual(
      values = cols,
      labels = labs,
      drop = TRUE,
      name = "Most likely class"
    ) +
    ggplot2::scale_alpha_continuous(
      range = c(0.3, 1),
      limits = c(0.25, 1),
      oob = scales::squish,
      name = "Class probability"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid = ggplot2::element_blank()
    )
}

#' Fig. EG (figs. GM, KC and MB combined): (a) the map of the
#' deployed partition (`map_plot`, build_figS_group_map()); (b) membership
#' probability per country and class in the deployed fit (membership_heatmap());
#' (c) most probable class per country at K in `Ks` (fig. KC's tiles; K = 5 is
#' panel b); (d) with `paths = TRUE`, the class paths across the class models
#' (membership_paths()) under the map. Panels b and c share the country rows
#' of panel b.
build_figS_enso_groups <- function(
  map_plot,
  m,
  kc_d,
  deployed,
  group_cols,
  group_lab,
  theme,
  low = 0.6,
  Ks = 2:4,
  base_size = 7,
  paths = FALSE
) {
  hm <- membership_heatmap(
    m,
    group_lab,
    theme,
    low = low,
    base_size = base_size,
    tag = "b"
  )
  rows <- levels(hm$data$iso3)
  kc_d <- kc_d[kc_d$K %in% Ks & kc_d$iso3 %in% rows, , drop = FALSE]
  kc <- build_figS_kl_classes_by_K(
    kc_d,
    deployed,
    group_cols,
    group_lab,
    theme,
    base_size = base_size
  ) +
    ggplot2::scale_y_discrete(limits = rows, labels = NULL) +
    ggplot2::labs(tag = "c") +
    # classes take the group colours of panel a's legend; the numbers are the probabilities
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )
  # the dagger is explained in the figure caption
  left <- map_plot + ggplot2::labs(tag = "a", caption = NULL)
  if (paths) {
    left <- left /
      (membership_paths(
        m,
        group_cols,
        group_lab,
        theme,
        base_size = base_size,
        tag = "d",
        ref_label = "final"
      ) +
        # line colours are the group colours of panel a's legend
        ggplot2::theme(legend.position = "none")) +
      patchwork::plot_layout(heights = c(1, 1.1))
  }
  patchwork::wrap_plots(left, hm, kc, nrow = 1, widths = c(2.4, 1, 0.75))
}

#' Fig. MP: membership_paths() alone, the class of each country
#' across the class models compared (fig. MB panel b, untagged).
build_figS_membership_paths <- function(
  m,
  group_cols,
  group_lab,
  theme,
  base_size = 7
) {
  membership_paths(
    m,
    group_cols,
    group_lab,
    theme,
    base_size = base_size,
    tag = NULL
  )
}

# ---------------------------------------------------------------------------
# Supplement figures (10 --figures DC,SV,FA,AA): data coverage, smoothing validation,
# forecast and attribution
# by analogue set.
# ---------------------------------------------------------------------------

#' The two tables of fig. DC from report_dataset$outbreaks (one row per analysis
#' outbreak: country, enso_group, admin_level, year_seas, cases).
#' @return list(country = tibble(enso_group, country, year_seas, n_episodes,
#'   cases), level = tibble(admin_level, year_seas, n_episodes))
data_coverage_frames <- function(outbreaks) {
  ob <- dplyr::mutate(
    outbreaks,
    enso_group = as.character(enso_group),
    country = as.character(country),
    admin_level = as.character(admin_level),
    year_seas = as.integer(year_seas)
  ) |>
    dplyr::filter(!is.na(year_seas))
  list(
    country = ob |>
      dplyr::group_by(enso_group, country, year_seas) |>
      dplyr::summarise(
        n_episodes = dplyr::n(),
        cases = sum(cases, na.rm = TRUE),
        .groups = "drop"
      ),
    level = ob |>
      dplyr::count(admin_level, year_seas, name = "n_episodes")
  )
}

#' Fig. DC: (A) canonical episodes (or cases) per country and cholera season
#' year, log fill, countries in blocks by ENSO group (display order and labels
#' of the caller's group_display_meta()), countries by name within a block;
#' (B) episodes per season year by administrative level (stacked bars).
#' @param frames data_coverage_frames()
#' @param group_order group ids in display order; @param group_lab named
#'   display labels (block strips)
#' @param country_lab optional ISO3 -> name (rows; the ISO3 when absent)
#' @param fill "episodes" | "cases"
#' @param level_cols optional named colours by admin level
#' @return a patchwork (A / B), tagged A, B
build_figS_data_coverage <- function(
  frames,
  group_order,
  group_lab,
  theme,
  country_lab = NULL,
  fill = c("episodes", "cases"),
  level_cols = NULL,
  base_size = 7
) {
  fill <- match.arg(fill)
  d <- frames$country
  d$value <- if (fill == "episodes") d$n_episodes else d$cases
  d <- d[d$value > 0, , drop = FALSE]
  ord <- c(
    intersect(group_order, unique(d$enso_group)),
    sort(setdiff(unique(stats::na.omit(d$enso_group)), group_order))
  )
  if (anyNA(d$enso_group)) {
    d$enso_group[is.na(d$enso_group)] <- "(ungrouped)"
    ord <- c(ord, "(ungrouped)")
  }
  strip <- stats::setNames(
    ifelse(ord %in% names(group_lab), unname(group_lab[ord]), ord),
    ord
  )
  name_of <- function(x) {
    n <- if (is.null(country_lab)) x else unname(country_lab[x])
    gsub("\n", " ", ifelse(is.na(n), x, n))
  }
  d$cname <- name_of(d$country)
  # Rows top to bottom by name within a block (the facet rows keep the blocks apart).
  d$cname <- factor(d$cname, levels = rev(sort(unique(d$cname))))
  d$enso_group <- factor(d$enso_group, levels = ord)
  yrs <- range(c(d$year_seas, frames$level$year_seas), na.rm = TRUE)
  xsc <- ggplot2::scale_x_continuous(
    limits = yrs + c(-0.5, 0.5),
    breaks = seq(yrs[1], yrs[2], by = 2),
    expand = c(0, 0)
  )
  pa <- ggplot2::ggplot(d, ggplot2::aes(year_seas, cname, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::facet_grid(
      enso_group ~ .,
      scales = "free_y",
      space = "free_y",
      labeller = ggplot2::as_labeller(strip)
    ) +
    ggplot2::scale_fill_viridis_c(
      trans = "log10",
      option = "rocket",
      direction = -1,
      end = 0.92,
      name = if (fill == "episodes") "Episodes" else "Suspected cases",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      guide = ggplot2::guide_colourbar(
        barheight = grid::unit(25, "mm"),
        barwidth = grid::unit(2.5, "mm")
      )
    ) +
    xsc +
    ggplot2::labs(x = NULL, y = NULL, tag = "a") +
    theme +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(
        angle = 0,
        hjust = 0,
        size = base_size - 1
      ),
      axis.text.y = ggplot2::element_text(size = base_size - 1.5),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(1, "mm")
    )
  lv <- frames$level
  lv$admin_level <- factor(
    lv$admin_level,
    levels = sort(unique(lv$admin_level))
  )
  pb <- ggplot2::ggplot(
    lv,
    ggplot2::aes(year_seas, n_episodes, fill = admin_level)
  ) +
    ggplot2::geom_col(width = 0.8) +
    (if (is.null(level_cols)) {
      ggplot2::scale_fill_grey(
        start = 0.25,
        end = 0.75,
        name = "Administrative\nlevel"
      )
    } else {
      ggplot2::scale_fill_manual(
        values = level_cols,
        name = "Administrative\nlevel"
      )
    }) +
    xsc +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(x = "Cholera season year", y = "Episodes", tag = "b") +
    theme +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = grid::unit(3, "mm")
    )
  patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(5.5, 1))
}

#' Fig. PC frames: the processing cascade (09 report_dataset$cascade, long
#' step / metric / value) as the boxes of a flow diagram. Every number is a
#' stamp or the difference of two nested stamps (records outside the WHO
#' regions, outbreaks outside Africa); a missing stamp prints "--".
#' @param cascade data.frame(step, metric, value)
#' @return list(main = spine boxes (row, label), excl = side boxes (row,
#'   label; drawn between `row` and `row + 1`), term = dataset boxes (x, label))
processing_cascade_frames <- function(cascade) {
  v <- function(step, metric) {
    x <- cascade$value[cascade$step == step & cascade$metric == metric]
    if (length(x) == 1L && !is.na(x)) x else NA_real_
  }
  f <- function(x) {
    ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = ","))
  }
  s00 <- "00_ingest_outbreakextractr"
  s01 <- "01_preprocess_outbreaks"
  s02 <- "02_preprocess_crosswalk"
  s03 <- "03_prepare_siz"
  main <- tibble::tibble(
    row = 1:6,
    label = c(
      sprintf(
        "Cholera time-series database\n%s weekly records",
        f(v(s00, "n_rows_read"))
      ),
      sprintf(
        "WHO African and Eastern Mediterranean regions\n%s weekly records | %s location periods | %s countries",
        f(v(s00, "n_rows_afr_emr")),
        f(v(s00, "n_lp_preoutbreak")),
        f(v(s00, "n_countries_preoutbreak"))
      ),
      sprintf(
        "Location-period episodes extracted (at least 50 suspected cases)\n%s episodes | %s location periods | %s countries",
        f(v(s00, "n_outbreaks")),
        f(v(s00, "n_lp_outbreak")),
        f(v(s00, "n_countries_outbreak"))
      ),
      sprintf(
        "Location-period episodes in Africa\n%s episodes | %s location periods | %s countries",
        f(v(s01, "n_outbreaks_africa")),
        f(v(s01, "n_lp_outbreaks_africa")),
        f(v(s01, "n_countries_outbreaks_africa"))
      ),
      sprintf(
        "Location periods mapped to administrative units\n%s location periods | %s units",
        f(v(s02, "n_lp_matched")),
        f(v(s02, "n_gids_matched"))
      ),
      sprintf(
        "Outbreak episodes (per unit and cholera year)\n%s episodes from %s location-period episodes\n(%s combine several)",
        f(v(s03, "n_episodes_canonical")),
        f(v(s03, "n_members")),
        f(v(s03, "n_multi_member_episodes"))
      )
    )
  )
  # Side boxes: what leaves the spine between two boxes (differences of the
  # counts above). 02's crosswalk roster holds sub-national location periods
  # only (national ones leave at row 4); at row 5, sub-national = the episodes
  # given a season calendar (matched + January default), then 03's input after
  # the first outbreak per location period and cholera year, then the members
  # left once the composite location periods are dropped.
  sub <- v(s03, "n_seas_episodes_matched") +
    v(s03, "n_seas_episodes_defaulted_jan")
  excl <- tibble::tibble(
    row = c(1L, 3L, 4L, 5L, 5L, 5L),
    label = c(
      sprintf(
        "Outside the two WHO regions: %s records",
        f(v(s00, "n_rows_read") - v(s00, "n_rows_afr_emr"))
      ),
      sprintf(
        "Outside Africa: %s location-period episodes",
        f(v(s01, "n_outbreaks_in") - v(s01, "n_outbreaks_africa"))
      ),
      sprintf(
        "National-level location periods: %s",
        f(v(s01, "n_lp_outbreaks_africa") - v(s02, "n_lp_in"))
      ),
      sprintf(
        "National level or without population:\n%s location-period episodes",
        f(v(s01, "n_outbreaks_africa") - sub)
      ),
      sprintf(
        "Later outbreak of a location period\nin the same cholera year: %s",
        f(sub - v(s03, "n_rows_in"))
      ),
      sprintf(
        "In the %s composite location periods\n(of %s; %s cases): %s",
        f(v(s03, "n_lp_dropped")),
        f(v(s03, "n_lp_total")),
        f(v(s03, "cases_dropped")),
        f(v(s03, "n_rows_in") - v(s03, "n_members"))
      )
    )
  )
  term <- tibble::tibble(
    x = c(-1.45, 0, 1.45),
    label = c(
      sprintf(
        "Occurrence\n%s unit-months\n(%s national-zero)\n%s onsets\n(+%s in a running\nepisode, not counted)",
        f(v("03_prepare_occ", "n_rows_out")),
        f(v("03_prepare_occ", "n_nz_rows")),
        f(v("03_prepare_occ", "n_events_out")),
        f(v("03_prepare_occ", "n_onsets_masked_running_episode"))
      ),
      sprintf(
        "Duration\n%s episode-months\n%s ends\n%s censored episodes",
        f(v("03_prepare_dur", "n_rows_out")),
        f(v("03_prepare_dur", "n_events_out")),
        f(v("03_prepare_dur", "n_censored"))
      ),
      sprintf(
        "Monthly case counts\n%s episode-months\n%s suspected cases",
        f(v("03_prepare_mon", "n_rows_out")),
        f(v("03_prepare_mon", "cases_total"))
      )
    )
  )
  list(main = main, excl = excl, term = term)
}

#' Fig. PC: the processing cascade as a flow diagram (the descriptive
#' notebook's CONSORT layout): a spine of boxes, side boxes for what leaves
#' it, and a fork into the three analysis datasets.
#' @param frames processing_cascade_frames()
#' @param base_size font size in pt of the spine boxes
#' @param family font family (10's FONT; "" = the device default)
build_figS_processing_cascade <- function(frames, base_size = 7, family = "") {
  # side boxes after one spine box stack between it and the next, one step
  # apart (taller for two-line labels); the spine gap widens to hold a stack
  excl <- frames$excl |>
    dplyr::group_by(row) |>
    dplyr::mutate(
      k = dplyr::row_number(),
      n = dplyr::n(),
      step = 0.3 +
        0.17 * (max(lengths(strsplit(label, "\n", fixed = TRUE))) - 1)
    ) |>
    dplyr::ungroup()
  span <- vapply(
    seq_len(nrow(frames$main)),
    function(r) {
      e <- excl[excl$row == r, ]
      if (nrow(e)) (e$n[1] - 1) * e$step[1] + 0.6 else 0
    },
    0
  )
  gap <- 1 + pmax(0, span - 0.5)
  main <- dplyr::mutate(
    frames$main,
    y = -cumsum(c(1, utils::head(gap, -1))),
    y_next = dplyr::lead(y)
  )
  excl <- excl |>
    dplyr::left_join(
      dplyr::select(main, row, y0 = y, y1 = y_next),
      by = "row"
    ) |>
    dplyr::mutate(y = (y0 + y1) / 2 + (n - 1) * step / 2 - (k - 1) * step)
  box_h <- 0.24
  y_end <- min(main$y)
  y_frk <- y_end - 1.15
  y_bar <- y_end - 0.5
  pt <- base_size / ggplot2::.pt
  col_line <- "grey30"
  arr <- grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dplyr::filter(main, row < max(row)),
      ggplot2::aes(x = 0, xend = 0, y = y - box_h, yend = y_next + box_h),
      arrow = arr,
      linewidth = 0.35,
      colour = col_line
    ) +
    ggplot2::geom_segment(
      data = excl,
      ggplot2::aes(x = 0, xend = 1.08, y = y, yend = y),
      arrow = arr,
      linewidth = 0.3,
      colour = "grey55"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 0, y = y_end - box_h, yend = y_bar),
      linewidth = 0.35,
      colour = col_line
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min(frames$term$x),
        xend = max(frames$term$x),
        y = y_bar,
        yend = y_bar
      ),
      linewidth = 0.35,
      colour = col_line
    ) +
    ggplot2::geom_segment(
      data = frames$term,
      ggplot2::aes(x = x, xend = x, y = y_bar, yend = y_frk + 0.42),
      arrow = arr,
      linewidth = 0.35,
      colour = col_line
    ) +
    ggplot2::geom_label(
      data = main,
      ggplot2::aes(x = 0, y = y, label = label),
      size = pt,
      lineheight = 1.05,
      fill = "white",
      label.size = 0.3,
      label.padding = grid::unit(0.38, "lines"),
      family = family
    ) +
    ggplot2::geom_label(
      data = excl,
      ggplot2::aes(x = 1.12, y = y, label = label),
      size = 0.9 * pt,
      hjust = 0,
      colour = "grey15",
      fill = "grey96",
      family = family,
      label.size = 0.2,
      label.padding = grid::unit(0.3, "lines")
    ) +
    ggplot2::geom_label(
      data = frames$term,
      ggplot2::aes(x = x, y = y_frk + 0.4, label = label),
      vjust = 1,
      size = pt,
      lineheight = 1.05,
      fill = "#dff0e4",
      label.size = 0.3,
      family = family,
      label.padding = grid::unit(0.38, "lines")
    ) +
    ggplot2::scale_x_continuous(limits = c(-2.05, 3.6), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      limits = c(y_frk - 1, -0.55),
      expand = c(0, 0)
    ) +
    ggplot2::theme_void(base_size = base_size, base_family = family)
}

#' Fig. SV's PIT histogram rows: smoothing$validation$pit_hist (stored at a
#' fine resolution) re-binned to `bins` equal bins (bins must divide the
#' stored count), as density (1 = uniform).
pit_hist_rebin <- function(pit_hist, bins = 20L) {
  nb <- max(pit_hist$bin)
  if (nb %% bins != 0L) {
    stop(
      "pit_hist_rebin: ",
      bins,
      " bins do not divide the stored ",
      nb,
      ".",
      call. = FALSE
    )
  }
  pit_hist |>
    dplyr::mutate(b = (bin - 1L) %/% (nb %/% bins) + 1L) |>
    dplyr::group_by(scope, metric, b) |>
    dplyr::summarise(
      lo = min(lo),
      hi = max(hi),
      n = sum(n),
      .groups = "drop_last"
    ) |>
    dplyr::mutate(density = if (sum(n) > 0) n / sum(n) * bins else NA_real_) |>
    dplyr::ungroup() |>
    dplyr::rename(bin = b)
}

#' Fig. SV: the historical reconstruction against the observations, from
#' smoothing$validation (09, smoothing v4; scopes already display-labelled):
#' (A) randomised PIT histograms of the monthly onsets and cases with the
#' U(0, 1) reference; (B) onset frequency by calendar month; (C) episode
#' durations; (D) fraction of covered unit-months active, Africa and each
#' group; (E) onset frequency given an ongoing outbreak in the previous month
#' at three levels. A, B, C, E for `scope`; reconstruction = median with 50%
#' and 95% intervals across replicates, observed in black.
#' @param v smoothing$validation; @param scope_order scope levels (Africa first)
#' @param pit_bins display bins of A; @param max_duration C's last duration
#' @param sim_col colour of the reconstruction
build_figS_smoothing_validation <- function(
  v,
  scope_order,
  theme,
  scope = "Africa",
  pit_bins = 20L,
  max_duration = 24L,
  sim_col = "grey25",
  base_size = 7
) {
  arm_lab <- c(sim = "Reconstruction (median, 50% / 95%)", obs = "Observed")
  arm_scale <- list(
    ggplot2::scale_colour_manual(
      values = c(sim = sim_col, obs = "black"),
      labels = arm_lab,
      name = NULL,
      limits = names(arm_lab)
    ),
    ggplot2::scale_shape_manual(
      values = c(sim = 16, obs = 4),
      labels = arm_lab,
      name = NULL,
      limits = names(arm_lab)
    )
  )
  pr <- function(d, x, dodge = 0) {
    list(
      ggplot2::geom_linerange(
        data = d,
        ggplot2::aes(x = .data[[x]], ymin = q025, ymax = q975, colour = "sim"),
        linewidth = 0.3,
        position = ggplot2::position_nudge(x = -dodge)
      ),
      ggplot2::geom_linerange(
        data = d,
        ggplot2::aes(x = .data[[x]], ymin = q25, ymax = q75, colour = "sim"),
        linewidth = 0.9,
        position = ggplot2::position_nudge(x = -dodge)
      ),
      ggplot2::geom_point(
        data = d,
        ggplot2::aes(x = .data[[x]], y = q50, colour = "sim", shape = "sim"),
        size = 0.9,
        position = ggplot2::position_nudge(x = -dodge)
      ),
      ggplot2::geom_point(
        data = d,
        ggplot2::aes(x = .data[[x]], y = obs, colour = "obs", shape = "obs"),
        size = 1.3,
        stroke = 0.5,
        position = ggplot2::position_nudge(x = dodge)
      )
    )
  }
  met_lab <- c(onsets = "Monthly onsets", cases = "Monthly cases")
  ph <- pit_hist_rebin(dplyr::filter(v$pit_hist, scope == !!scope), pit_bins) |>
    dplyr::mutate(metric = factor(met_lab[metric], met_lab))
  pa <- ggplot2::ggplot(ph) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = lo, xmax = hi, ymin = 0, ymax = density),
      fill = "grey65",
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey20",
      linewidth = 0.3
    ) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      expand = ggplot2::expansion(add = 0.02)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(x = "Randomised PIT", y = "Density", tag = "a") +
    theme
  se <- dplyr::filter(v$seasonal, scope == !!scope)
  pb <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = se,
      ggplot2::aes(month, ymin = q025, ymax = q975),
      fill = sim_col,
      alpha = 0.15
    ) +
    ggplot2::geom_ribbon(
      data = se,
      ggplot2::aes(month, ymin = q25, ymax = q75),
      fill = sim_col,
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = se,
      ggplot2::aes(month, q50, colour = "sim"),
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = se,
      ggplot2::aes(month, obs, colour = "obs", shape = "obs"),
      size = 1.3,
      stroke = 0.5,
      show.legend = FALSE
    ) +
    arm_scale +
    ggplot2::scale_x_continuous(
      breaks = 1:12,
      labels = substr(month.abb, 1, 1),
      expand = c(0.02, 0)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::labs(
      x = "Calendar month",
      y = "Onsets per at-risk unit-month",
      tag = "b"
    ) +
    theme
  du <- dplyr::filter(v$duration, scope == !!scope, duration <= max_duration)
  pc <- ggplot2::ggplot() +
    pr(du, "duration", dodge = 0.15) +
    arm_scale +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      x = "Episode duration (months)",
      y = "Share of episodes",
      tag = "c"
    ) +
    theme
  ac <- dplyr::filter(v$active, scope %in% scope_order) |>
    dplyr::mutate(scope = factor(scope, levels = scope_order))
  pd <- ggplot2::ggplot() +
    pr(ac, "scope", dodge = 0.12) +
    arm_scale +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::labs(x = NULL, y = "Covered unit-months active", tag = "d") +
    theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ind_lab <- c(
    adm1 = "same first-level\nunit",
    country = "other first-level unit,\nsame country",
    tu = "other country,\nsame transmission unit"
  )
  st_lab <- c(
    ongoing = "Outbreak ongoing in the previous month",
    none = "None ongoing"
  )
  co <- dplyr::filter(v$contagion, scope == !!scope) |>
    dplyr::mutate(
      indicator = factor(ind_lab[indicator], ind_lab),
      state = factor(st_lab[state], st_lab)
    )
  pe <- ggplot2::ggplot() +
    pr(co, "indicator", dodge = 0.12) +
    arm_scale +
    ggplot2::facet_wrap(~state, nrow = 1, scales = "free_y") +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 0.01)
    ) +
    ggplot2::labs(x = NULL, y = "Onsets per at-risk unit-month", tag = "e") +
    theme
  patchwork::wrap_plots(
    pa,
    pb,
    pc,
    pd,
    pe,
    design = "AB\nCD\nEE",
    heights = c(1, 1, 1.1),
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# Hindcast skill against climatology (10 figure H, --supplementary true;
# hindcast_occ$skill and hindcast_cases$month_skill, 09)
# ---------------------------------------------------------------------------

#' Rows of the hindcast skill-vs-climatology figure: the occurrence CRPSS of
#' the onset count per unit and season window (grain "window") against
#' `occ_reference`, and the monthly case CRPSS of the country-months of
#' `case_support` against `case_reference`, on raw counts and on log(1 + y);
#' Africa and each region, by lead class, with the season-year bootstrap
#' interval. Lead classes outside `lead_levels` are dropped.
#' @param occ_skill hindcast_occ$skill; @param month_skill hindcast_cases$month_skill
#' @param occ_rule the occurrence rule (crps_fair, the fair estimator)
#' @return tibble(target (onsets | cases_raw | cases_log), scope, lead_class,
#'   crpss, crpss_lo, crpss_hi, n_scored)
hindcast_clim_skill_frame <- function(
  occ_skill,
  month_skill,
  occ_reference = "loo_gid",
  case_reference = "clim_pooled",
  case_support = "atrisk_gids",
  lead_levels = c("0-3", "4-8", "9+", "all"),
  occ_rule = "crps_fair"
) {
  keep <- c("scope", "lead_class", "crpss", "crpss_lo", "crpss_hi", "n_scored")
  occ <- occ_skill |>
    dplyr::filter(
      scope_level %in% c("africa", "region"),
      grain == "window",
      rule == occ_rule,
      reference == occ_reference,
      lead_class %in% lead_levels
    ) |>
    dplyr::select(dplyr::all_of(keep)) |>
    dplyr::mutate(target = "onsets")
  cas <- month_skill |>
    dplyr::filter(
      scope_level %in% c("africa", "region"),
      grain == "country",
      support == case_support,
      reference == case_reference,
      rule %in% c("crps_raw", "crps_log1p"),
      lead_class %in% lead_levels
    ) |>
    dplyr::mutate(
      target = ifelse(rule == "crps_raw", "cases_raw", "cases_log")
    ) |>
    dplyr::select(target, dplyr::all_of(keep))
  if (!nrow(occ) || !nrow(cas)) {
    stop(
      "[figS_hindcast_clim_skill] no rows for occurrence reference '",
      occ_reference,
      "' or case reference '",
      case_reference,
      "' / support '",
      case_support,
      "'.",
      call. = FALSE
    )
  }
  dplyr::bind_rows(occ, cas) |>
    dplyr::mutate(lead_class = factor(lead_class, levels = lead_levels)) |>
    dplyr::relocate(target)
}

#' Hindcast CRPSS against climatology, one panel per target (onsets; monthly
#' cases, raw; monthly cases, log), lead class on x, Africa and each region
#' dodged, point = CRPSS, bar = bootstrap interval. Intervals below `ylim_lo`
#' are clipped at it and marked with an open triangle.
#' @param d hindcast_clim_skill_frame(); @param scope_lab named labels by scope
#'   key ("africa" + region ids), in display order; @param scope_col colours keyed
#'   like scope_lab; @param target_lab panel titles by target
build_figS_hindcast_clim_skill <- function(
  d,
  scope_lab,
  scope_col,
  theme,
  target_lab = NULL,
  lead_lab = c(
    `0-3` = "0\u20133",
    `4-8` = "4\u20138",
    `9+` = "9+",
    all = "All leads"
  ),
  ylim_lo = -1,
  base_size = 7
) {
  if (is.null(target_lab)) {
    target_lab <- c(
      onsets = "Onsets per unit and season window",
      cases_raw = "Monthly cases per country (raw counts)",
      cases_log = "Monthly cases per country (log(1 + cases))"
    )
  }
  d <- d |>
    dplyr::filter(
      scope %in% names(scope_lab),
      target %in% names(target_lab),
      is.finite(crpss)
    ) |>
    dplyr::mutate(
      scope = factor(unname(scope_lab[scope]), levels = unname(scope_lab)),
      target = factor(unname(target_lab[target]), levels = unname(target_lab)),
      clipped = is.finite(crpss_lo) & crpss_lo < ylim_lo,
      lo = pmax(crpss_lo, ylim_lo)
    )
  cols <- stats::setNames(
    unname(scope_col[names(scope_lab)]),
    unname(scope_lab)
  )
  dodge <- ggplot2::position_dodge(width = 0.75)
  ggplot2::ggplot(
    d,
    ggplot2::aes(lead_class, crpss, colour = scope, group = scope)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = lo, ymax = crpss_hi),
      linewidth = 0.4,
      position = dodge
    ) +
    ggplot2::geom_point(size = 1.4, position = dodge) +
    # every row, so the dodge matches the intervals; NA (unclipped) rows draw nothing
    ggplot2::geom_point(
      ggplot2::aes(y = ifelse(clipped, lo, NA_real_)),
      shape = 6,
      size = 1.2,
      position = dodge,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~target, ncol = 1) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_x_discrete(labels = lead_lab) +
    ggplot2::coord_cartesian(ylim = c(ylim_lo, 1)) +
    ggplot2::labs(
      x = "Lead (months after the country's launch)",
      y = "CRPSS against climatology"
    ) +
    theme +
    ggplot2::theme(legend.position = "bottom")
}

#' Fig. FA's series from the attribution artefact: the monthly ribbons of the
#' attribution factual (attrib_ribbon / attrib_africa_ribbon when the artefact's
#' paired statistics read an attribution factual, else the displayed forecast)
#' and of the counterfactual of each analogue set in `sets`, Africa and each
#' group in `groups`, from the attribution start to `end` (NULL: the series'
#' end). The factual rows are repeated under every set.
#' @return tibble(scope ("Africa" | group id), arm (factual | counterfactual),
#'   analog_set, month, metric (cases | onsets), q025, q25, q50, q75, q975)
forecast_analogue_frame <- function(f23, groups, sets, end = NULL) {
  af <- if (is.null(f23$attrib_factual)) "forecast" else f23$attrib_factual
  attrib <- identical(af, "attrib") && !is.null(f23$attrib_ribbon)
  fac_reg <- if (attrib) {
    f23$attrib_ribbon
  } else {
    dplyr::filter(f23$region_ribbon, run == "forecast")
  }
  fac_afr <- if (attrib) {
    f23$attrib_africa_ribbon
  } else {
    dplyr::filter(f23$africa_ribbon, run == "forecast")
  }
  tag <- function(d, scope_col, arm) {
    d <- dplyr::ungroup(d)
    d$scope <- if (is.null(scope_col)) {
      "Africa"
    } else {
      as.character(d[[scope_col]])
    }
    d$arm <- arm
    d
  }
  fac <- dplyr::bind_rows(
    tag(
      dplyr::filter(fac_reg, enso_group %in% groups),
      "enso_group",
      "factual"
    ),
    tag(fac_afr, NULL, "factual")
  )
  fac <- dplyr::bind_rows(lapply(sets, function(s) {
    dplyr::mutate(fac, analog_set = s)
  }))
  cf <- dplyr::bind_rows(
    tag(
      dplyr::filter(
        f23$region_ribbon,
        run == "counterfactual",
        enso_group %in% groups
      ),
      "enso_group",
      "counterfactual"
    ),
    tag(
      dplyr::filter(f23$africa_ribbon, run == "counterfactual"),
      NULL,
      "counterfactual"
    )
  ) |>
    dplyr::filter(analog_set %in% sets)
  d <- dplyr::bind_rows(fac, cf) |>
    dplyr::mutate(month = as.Date(month)) |>
    dplyr::filter(
      month >= as.Date(f23$attrib_start),
      if (is.null(end)) TRUE else month <= as.Date(end)
    )
  q <- c("q025", "q25", "q50", "q75", "q975")
  dplyr::bind_rows(lapply(c("cases", "onsets"), function(m) {
    out <- d[, c("scope", "arm", "analog_set", "month")]
    out$metric <- m
    for (k in q) {
      out[[k]] <- unname(as.numeric(d[[paste0(m, "_", k)]]))
    }
    out
  })) |>
    dplyr::select(scope, arm, analog_set, month, metric, dplyr::all_of(q))
}

#' Fig. FA: monthly onsets (a) and cases (b) under the attribution factual
#' (the scope's colour, solid, as Figure 4) and the counterfactual (cf_col,
#' dashed) of each analogue set (columns), Africa and each ENSO group (rows):
#' median with 50% and 95% bands.
#' @param d forecast_analogue_frame()
#' @param scope_lab named display labels ("Africa" and group ids), in row order
#' @param set_lab named display labels of the sets, in column order
#' @param scope_col named colours keyed like scope_lab (the factual series;
#'   NULL: near-black for every scope)
#' @param arm_lab named c(factual =, counterfactual =) legend labels
build_figS_forecast_analogues <- function(
  d,
  scope_lab,
  set_lab,
  theme,
  cf_col = "grey35",
  scope_col = NULL,
  arm_lab = c(factual = "Factual", counterfactual = "Counterfactual"),
  metrics = c(
    onsets = "Outbreak onsets per month",
    cases = "Suspected cases per month"
  ),
  base_size = 7
) {
  if (is.null(scope_col)) {
    scope_col <- stats::setNames(
      rep("grey10", length(scope_lab)),
      names(scope_lab)
    )
  }
  cols <- c(scope_col[names(scope_lab)], counterfactual = cf_col)
  one <- function(m, tag) {
    x <- dplyr::filter(
      d,
      metric == m,
      scope %in% names(scope_lab),
      analog_set %in% names(set_lab)
    ) |>
      dplyr::mutate(
        key = ifelse(arm == "factual", as.character(scope), "counterfactual"),
        scope = factor(scope, names(scope_lab)),
        analog_set = factor(analog_set, names(set_lab)),
        arm = factor(arm, names(arm_lab))
      )
    ggplot2::ggplot(x, ggplot2::aes(month, group = arm)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = q025, ymax = q975, fill = key),
        alpha = 0.15
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = q25, ymax = q75, fill = key),
        alpha = 0.3
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = q50, colour = key, linetype = arm),
        linewidth = 0.4
      ) +
      ggplot2::facet_grid(
        scope ~ analog_set,
        scales = "free_y",
        labeller = ggplot2::labeller(
          scope = ggplot2::as_labeller(scope_lab),
          analog_set = ggplot2::as_labeller(set_lab)
        )
      ) +
      ggplot2::scale_colour_manual(values = cols, guide = "none") +
      ggplot2::scale_fill_manual(values = cols, guide = "none") +
      ggplot2::scale_linetype_manual(
        values = c(factual = "solid", counterfactual = "22"),
        labels = arm_lab,
        name = NULL,
        guide = ggplot2::guide_legend(
          override.aes = list(colour = c("grey10", cf_col))
        )
      ) +
      ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      ggplot2::labs(x = NULL, y = metrics[[m]], tag = tag) +
      theme +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(
          angle = 0,
          hjust = 0,
          size = base_size - 1
        ),
        axis.text = ggplot2::element_text(size = base_size - 1.5)
      )
  }
  patchwork::wrap_plots(
    one(names(metrics)[1], "a"),
    one(names(metrics)[2], "b"),
    ncol = 1,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Fig. FA, overlay variant (build_figS_forecast_analogues is the other): one panel
#' column per metric (a onsets, b cases), Africa and each
#' ENSO group as rows. The attribution factual (fac_col, solid, median with 50%
#' and 95% bands) once, and the counterfactual median of every analogue set
#' (set_col, dashed) overlaid, so the sets compare on one axis.
#' @param d forecast_analogue_frame()
#' @param scope_lab named display labels ("Africa" and group ids), in row order
#' @param set_lab named display labels of the sets, in legend order
#' @param set_col named colours keyed like set_lab
#' @param fac_lab legend label of the factual arm
build_figS_forecast_analogues_overlay <- function(
  d,
  scope_lab,
  set_lab,
  set_col,
  theme,
  fac_col = "grey10",
  fac_lab = "Factual",
  metrics = c(
    onsets = "Outbreak onsets per month",
    cases = "Suspected cases per month"
  ),
  base_size = 7
) {
  one <- function(m, tag) {
    x <- dplyr::filter(d, metric == m, scope %in% names(scope_lab)) |>
      dplyr::mutate(scope = factor(scope, names(scope_lab)))
    fac <- dplyr::filter(x, arm == "factual") |>
      dplyr::distinct(scope, month, .keep_all = TRUE)
    cf <- dplyr::filter(
      x,
      arm == "counterfactual",
      analog_set %in% names(set_lab)
    ) |>
      dplyr::mutate(series = factor(analog_set, names(set_lab)))
    fac$series <- factor("factual", c("factual", names(set_lab)))
    cf$series <- factor(as.character(cf$series), c("factual", names(set_lab)))
    cols <- c(factual = fac_col, set_col[names(set_lab)])
    labs_ <- c(factual = fac_lab, set_lab)
    ggplot2::ggplot(mapping = ggplot2::aes(month)) +
      ggplot2::geom_ribbon(
        data = fac,
        ggplot2::aes(ymin = q025, ymax = q975),
        fill = fac_col,
        alpha = 0.10
      ) +
      ggplot2::geom_ribbon(
        data = fac,
        ggplot2::aes(ymin = q25, ymax = q75),
        fill = fac_col,
        alpha = 0.18
      ) +
      ggplot2::geom_line(
        data = cf,
        ggplot2::aes(y = q50, colour = series, linetype = series),
        linewidth = 0.45
      ) +
      ggplot2::geom_line(
        data = fac,
        ggplot2::aes(y = q50, colour = series, linetype = series),
        linewidth = 0.5
      ) +
      ggplot2::facet_wrap(
        ~scope,
        ncol = 1,
        scales = "free_y",
        strip.position = "right",
        labeller = ggplot2::as_labeller(scope_lab)
      ) +
      ggplot2::scale_colour_manual(
        values = cols,
        labels = labs_,
        name = NULL,
        drop = FALSE
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          factual = "solid",
          stats::setNames(rep("22", length(set_lab)), names(set_lab))
        ),
        labels = labs_,
        name = NULL,
        drop = FALSE
      ) +
      ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      ggplot2::labs(x = NULL, y = metrics[[m]], tag = tag) +
      theme +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(
          angle = 0,
          hjust = 0,
          size = base_size - 1
        ),
        axis.text = ggplot2::element_text(size = base_size - 1.5)
      )
  }
  patchwork::wrap_plots(
    one(names(metrics)[1], "a"),
    one(names(metrics)[2], "b"),
    nrow = 1,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Fig. AA's ADM1 rows: report_attrib_exposure$adm1 at the headlined cell
#' (`window`, explosive `guard`) for one metric and every set in `sets`, ADM1s
#' with a ratio (0 vs 0 left out), joined onto their polygons. The ratio of
#' means and p_sign = P(factual > counterfactual) over the paired replicates do
#' not depend on the IPCC magnitude: one row per ADM1 x set.
#' @return sf (adm1_sf rows) with analog_set, ratio, p_sign
attribution_analogue_map_sf <- function(
  adm1,
  metric,
  sets,
  adm1_sf,
  window = "attribution",
  guard = TRUE
) {
  d <- adm1
  win_hl <- window
  if (all(c("window", "explosive_guard") %in% names(d))) {
    d <- dplyr::filter(d, window == win_hl, explosive_guard == guard)
  }
  d <- d |>
    dplyr::filter(metric == !!metric, analog_set %in% sets, !is.na(ratio)) |>
    dplyr::distinct(adm1_unit, analog_set, .keep_all = TRUE) |>
    dplyr::select(adm1_unit, analog_set, ratio, p_sign)
  dplyr::inner_join(adm1_sf, d, by = c("GID_1" = "adm1_unit"))
}

#' Fig. AA: per metric, two map columns -- the ADM1 ratio of mean factual to
#' mean counterfactual totals (log10, capped at `limits`; an infinite ratio at
#' the cap) and the probability that the factual exceeds the counterfactual
#' (p_sign) -- one row per analogue set.
#' @param maps named list metric -> attribution_analogue_map_sf()
#' @param set_lab named set labels, in row order
#' @param base_layers_fn function() returning fresh base layers (country greys
#'   + ggnewscale::new_scale_fill(); ggnewscale renames a used layer in place)
#' @param overlay_layers drawn over the fills (borders)
#' @param ratio_fill,prob_fill fill scales of the two columns (log10 ratio, 0-1)
#' @param metric_lab named column-title prefixes by metric
build_figS_attribution_analogues <- function(
  maps,
  set_lab,
  bbox,
  ratio_fill,
  prob_fill,
  limits = c(1 / 3, 3),
  base_layers_fn = function() NULL,
  overlay_layers = NULL,
  map_theme = ggplot2::theme_void(),
  metric_lab = c(cases = "Cases", onsets = "Outbreak onsets"),
  base_size = 7
) {
  one <- function(m, what, first) {
    x <- maps[[m]]
    x$analog_set <- factor(x$analog_set, names(set_lab))
    x$value <- if (what == "ratio") {
      log10(pmin(pmax(x$ratio, limits[1]), limits[2]))
    } else {
      x$p_sign
    }
    ggplot2::ggplot() +
      base_layers_fn() +
      ggplot2::geom_sf(
        data = x,
        ggplot2::aes(fill = value),
        colour = "white",
        linewidth = 0.005
      ) +
      (if (what == "ratio") ratio_fill else prob_fill) +
      overlay_layers +
      ggplot2::facet_wrap(
        ~analog_set,
        ncol = 1,
        strip.position = "left",
        drop = FALSE,
        labeller = ggplot2::as_labeller(set_lab)
      ) +
      ggplot2::coord_sf(
        xlim = bbox[c("xmin", "xmax")],
        ylim = bbox[c("ymin", "ymax")],
        expand = FALSE,
        datum = NA
      ) +
      ggplot2::labs(
        title = sprintf(
          "%s:\n%s",
          metric_lab[[m]],
          if (what == "ratio") {
            "ratio of means"
          } else {
            "probability of an increase"
          }
        )
      ) +
      map_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = base_size, hjust = 0.5),
        strip.text.y.left = if (first) {
          ggplot2::element_text(size = base_size, angle = 90)
        } else {
          ggplot2::element_blank()
        },
        strip.placement = "outside",
        panel.spacing = grid::unit(2, "mm"),
        legend.title = ggplot2::element_text(size = base_size - 1),
        legend.text = ggplot2::element_text(size = base_size - 1.5)
      )
  }
  ps <- list()
  for (m in names(maps)) {
    for (w in c("ratio", "prob")) {
      ps[[length(ps) + 1L]] <- one(m, w, !length(ps))
    }
  }
  # legend titles above their keys (beside them, the two-line ratio title ran
  # into the key's first label)
  patchwork::wrap_plots(ps, nrow = 1, guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.title.position = "top")
}

#' Fig. SE frames: 09's seasonality artefact in display order (pure).
#' @param art seasonality artefact (rain, rates, fitted)
#' @param group_order group ids in display order
#' @return list(rain, shade (the rainy-season months), rates, fitted), each
#'   with a factor `enso_group` in display order
seasonality_frames <- function(art, group_order) {
  g <- intersect(group_order, unique(art$rain$enso_group))
  fg <- function(d) {
    dplyr::mutate(
      dplyr::filter(d, enso_group %in% g),
      enso_group = factor(enso_group, levels = g)
    )
  }
  list(
    rain = fg(art$rain),
    shade = fg(dplyr::filter(art$rain, in_season)),
    rates = fg(art$rates),
    fitted = fg(art$fitted)
  )
}

#' Fig. SE: per ENSO group (columns), (a) the rainfall climatology with the
#' group's main rainy season shaded (Fig. 1d's months, shaded in every row),
#' then on the same calendar months (b) the observed onset rate with the
#' occurrence model's seasonal spline, (c) the observed cessation rate, with
#' the duration model's harmonic seasonality (95% band) for the groups whose
#' deployed duration fit carries it (duration.seasonal_term harmonic2; none
#' does in the deployed run, which draws no curve), (d) median cases
#' per active episode-month (IQR) with the monthly model's harmonic
#' seasonality (95% band). Rows b-d on a log axis; hollow points: months with
#' fewer than min_n events (episode-months in d). Row c: lower interval ends at
#' 0 (no event) are drawn to the row's floor, a fifth of its smallest positive
#' rate, where the axis starts. Rows b and d: the axis spans the observed
#' points and the fitted curve (and its band), widened by the factor `pad` on
#' each side; intervals reaching past it are clipped at the panel edge.
#' @param frames seasonality_frames()
#' @param group_lab named column labels; @param group_cols named colours
#' @param rain_fill fill of the rainy-season band; @param rain_col bar colour
#' @param min_n the artefact's sparse threshold (legend text)
#' @param pad multiplicative margin of the tight axes of rows b and d
#' @return a patchwork, rows tagged a-d
build_figS_seasonality <- function(
  frames,
  group_lab,
  group_cols,
  theme,
  rain_fill = "#dfe7ee",
  rain_col = "#6b7f91",
  min_n = 5L,
  base_size = 7,
  pad = 1.25
) {
  lv <- levels(frames$rain$enso_group)
  lab <- stats::setNames(unname(group_lab[lv]), lv)
  # ymin 0, not -Inf: log10(0) = -Inf is drawn at the panel edge, -Inf has no log.
  shade <- ggplot2::geom_rect(
    data = frames$shade,
    ggplot2::aes(
      xmin = month - 0.5,
      xmax = month + 0.5,
      ymin = 0,
      ymax = Inf,
      fill = "Rainy season (Fig. 1d)"
    ),
    inherit.aes = FALSE,
    colour = NA
  )
  xsc <- ggplot2::scale_x_continuous(
    breaks = 1:12,
    labels = substr(month.abb, 1, 1),
    limits = c(0.5, 12.5),
    expand = c(0, 0)
  )
  fsc <- ggplot2::scale_fill_manual(
    values = c("Rainy season (Fig. 1d)" = rain_fill),
    name = NULL
  )
  csc <- ggplot2::scale_colour_manual(values = group_cols, guide = "none")
  ssc <- ggplot2::scale_shape_manual(
    values = c(observed = 16, sparse = 1),
    name = NULL,
    drop = FALSE,
    labels = c(
      observed = "Observed",
      sparse = sprintf("Observed, fewer than %d events", min_n)
    )
  )
  lsc <- ggplot2::scale_linetype_manual(
    values = c(fitted = "solid"),
    labels = c(fitted = "Fitted (model)"),
    name = NULL
  )
  base <- function(p, ylab, tag, strip = FALSE, xlab = FALSE) {
    p +
      xsc +
      fsc +
      csc +
      ggplot2::facet_grid(
        . ~ enso_group,
        labeller = ggplot2::as_labeller(lab)
      ) +
      ggplot2::labs(
        x = if (xlab) "Calendar month" else NULL,
        y = ylab,
        tag = tag
      ) +
      theme +
      ggplot2::theme(
        strip.text.x = if (strip) {
          ggplot2::element_text(size = base_size)
        } else {
          ggplot2::element_blank()
        },
        axis.text.x = if (xlab) {
          ggplot2::element_text(size = base_size - 1.5)
        } else {
          ggplot2::element_blank()
        },
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.spacing.x = grid::unit(1.5, "mm")
      )
  }
  pa <- base(
    ggplot2::ggplot(frames$rain, ggplot2::aes(month, mm)) +
      shade +
      ggplot2::geom_col(width = 0.7, fill = rain_col) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ),
    "Rainfall\n(mm per month)",
    "a",
    strip = TRUE
  )
  rate_row <- function(
    track,
    ylab,
    tag,
    xlab = FALSE,
    legend = FALSE,
    tight = FALSE
  ) {
    r <- dplyr::filter(
      frames$rates,
      track == !!track,
      is.finite(rate),
      rate > 0
    )
    f <- dplyr::filter(frames$fitted, track == !!track)
    if (tight) {
      v <- c(r$rate, f$scaled, f$scaled_lo, f$scaled_hi)
      v <- v[is.finite(v) & v > 0]
      lim <- c(min(v) / pad, max(v) * pad)
    } else {
      lim <- c(min(r$rate) / 5, NA)
    }
    fl <- lim[1]
    r <- dplyr::mutate(
      r,
      kind = factor(
        ifelse(sparse, "sparse", "observed"),
        levels = c("observed", "sparse")
      ),
      lo = pmax(lo, fl)
    )
    p <- ggplot2::ggplot(r, ggplot2::aes(month, rate, colour = enso_group)) +
      shade
    if (any(is.finite(f$scaled_lo))) {
      p <- p +
        ggplot2::geom_ribbon(
          data = dplyr::filter(f, is.finite(scaled_lo)),
          ggplot2::aes(
            month,
            ymin = scaled_lo,
            ymax = scaled_hi,
            group = enso_group
          ),
          inherit.aes = FALSE,
          fill = "grey55",
          alpha = 0.25
        )
    }
    if (nrow(f)) {
      p <- p +
        ggplot2::geom_line(
          data = f,
          ggplot2::aes(month, scaled, group = enso_group, linetype = "fitted"),
          colour = "grey20",
          linewidth = 0.45,
          inherit.aes = FALSE
        ) +
        lsc
    }
    p <- p +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo, ymax = hi),
        linewidth = 0.3
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = kind),
        size = 1.1,
        stroke = 0.4
      ) +
      ssc +
      ggplot2::scale_y_log10(
        labels = scales::label_number(drop0trailing = TRUE),
        limits = lim,
        oob = scales::oob_squish
      )
    if (!legend) {
      p <- p + ggplot2::guides(shape = "none", linetype = "none", fill = "none")
    }
    base(p, ylab, tag, xlab = xlab)
  }
  pb <- rate_row(
    "occurrence",
    "Onsets per 1,000\nunit-months at risk",
    "b",
    legend = TRUE,
    tight = TRUE
  )
  pc <- rate_row(
    "cessation",
    "Episode ends per 100\nactive episode-months",
    "c"
  )
  pd <- rate_row(
    "cases",
    "Cases per active\nepisode-month",
    "d",
    xlab = TRUE,
    tight = TRUE
  )
  patchwork::wrap_plots(pa, pb, pc, pd, ncol = 1, guides = "collect") &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = grid::unit(3, "mm"),
      legend.text = ggplot2::element_text(size = base_size - 1)
    )
}

# --- Fig. WR: observed vs modelled ENSO ratio (09 report_who_enso_ratio) -----

#' The rows of fig. WR (pure): one group, one year definition, one analogue set
#' of report_who_enso_ratio$comparison, points with an observed ratio only.
#' @return tibble(estimator, scope, kind ("country" | "group"), countries,
#'   neutral_years, obs_ratio, obs_lo, obs_hi, model_ratio, model_lo, model_q25,
#'   model_q75, model_hi)
who_enso_ratio_frame <- function(
  comparison,
  group,
  year_definition,
  analog_set
) {
  d <- tibble::as_tibble(comparison)
  d <- d[
    d$enso_group == group &
      d$year_definition == year_definition &
      d$analog_set == analog_set &
      is.finite(d$obs_ratio),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) {
    stop(
      "who_enso_ratio_frame: no observed ratio for ",
      group,
      " / ",
      year_definition,
      " / ",
      analog_set,
      call. = FALSE
    )
  }
  d$kind <- d$level # 09's country / group flag
  d[, c(
    "estimator",
    "scope",
    "kind",
    "countries",
    "neutral_years",
    "obs_ratio",
    "obs_lo",
    "obs_hi",
    "model_ratio",
    "model_lo",
    "model_q25",
    "model_q75",
    "model_hi"
  )]
}

#' Fig. WR: per estimator panel, the observed ratio (target year / mean of the
#' neutral years, dashed bootstrap interval: not a posterior interval) against
#' the modelled ratio (smoothing / counterfactual ratio of means; thin bar 95%,
#' thick bar 50% of the per-smoothing-replicate spread), log axes shared by the
#' panels, the equality diagonal dashed. Panel b (ytd_annual) is a lower bound
#' (arrow). Countries are open circles, the group a filled diamond.
#' @param d who_enso_ratio_frame()
#' @param group_lab display label of the group; group_col its colour
#' @param titles named character, estimator -> panel title
#' @param model_lab y-axis title
build_figS_who_enso_ratio <- function(
  d,
  group_lab,
  group_col,
  theme,
  titles,
  model_lab,
  obs_lab = "Observed ratio: 2026 / mean of neutral years",
  base_size = 7,
  family = ""
) {
  ests <- intersect(c("same_window", "ytd_annual"), unique(d$estimator))
  lim_of <- function(v) {
    r <- log10(range(c(v[is.finite(v) & v > 0], 0.5, 2))) + c(-0.15, 0.15)
    10^c(floor(r[1] * 2) / 2, ceiling(r[2] * 2) / 2)
  }
  xl <- lim_of(d$obs_ratio)
  yl <- lim_of(c(d$model_lo, d$model_hi))
  brk <- function(l, steps) {
    b <- as.vector(outer(
      steps,
      10^seq(floor(log10(l[1])), ceiling(log10(l[2])))
    ))
    b[b >= l[1] & b <= l[2]]
  }
  lab <- scales::label_number(drop0trailing = TRUE, big.mark = ",")
  d$label <- ifelse(d$kind == "group", group_lab, d$scope)
  cols <- c(country = "grey15", group = group_col)
  panel <- function(e, tag, lower_bound) {
    x <- d[d$estimator == e, , drop = FALSE]
    p <- ggplot2::ggplot(x, ggplot2::aes(obs_ratio, model_ratio)) +
      ggplot2::geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.25) +
      ggplot2::geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.25) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "22",
        colour = "grey45",
        linewidth = 0.3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = obs_lo, xmax = obs_hi, colour = kind),
        orientation = "y",
        width = 0,
        linewidth = 0.3,
        linetype = "22",
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = kind, fill = kind),
        size = 1.8,
        colour = "grey15",
        stroke = 0.3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = model_lo, ymax = model_hi),
        colour = "grey10",
        width = 0.12,
        linewidth = 0.3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = model_q25, ymax = model_q75),
        colour = "grey10",
        width = 0,
        linewidth = 0.7
      ) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = label, colour = kind),
        size = (base_size - 0.5) / ggplot2::.pt,
        family = family,
        seed = 3,
        min.segment.length = 0,
        segment.size = 0.2,
        segment.colour = "grey50",
        box.padding = 0.25,
        point.padding = 0.3,
        nudge_x = -0.2,
        nudge_y = 0.17,
        hjust = 1,
        show.legend = FALSE
      ) +
      ggplot2::scale_x_log10(breaks = brk(xl, 1), labels = lab) +
      ggplot2::scale_y_log10(breaks = brk(yl, c(1, 2, 5)), labels = lab) +
      ggplot2::scale_colour_manual(values = cols, guide = "none") +
      ggplot2::scale_fill_manual(
        values = c(country = "white", group = group_col),
        name = NULL,
        labels = c(
          country = "Country",
          group = paste(group_lab, "(countries shown)")
        )
      ) +
      ggplot2::scale_shape_manual(
        values = c(country = 21, group = 23),
        name = NULL,
        labels = c(
          country = "Country",
          group = paste(group_lab, "(countries shown)")
        )
      ) +
      ggplot2::coord_cartesian(xlim = xl, ylim = yl, expand = FALSE) +
      ggplot2::labs(
        title = titles[[e]],
        x = obs_lab,
        y = model_lab,
        tag = tag
      ) +
      theme +
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = ggplot2::element_blank()
      )
    if (lower_bound) {
      p <- p +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = obs_ratio * 1.12,
            xend = obs_ratio * 1.6,
            yend = model_ratio,
            colour = kind
          ),
          linewidth = 0.3,
          arrow = grid::arrow(length = grid::unit(0.8, "mm"), type = "closed")
        )
    }
    p
  }
  ps <- lapply(seq_along(ests), function(i) {
    panel(ests[i], letters[i], identical(ests[i], "ytd_annual")) +
      (if (i > 1L) ggplot2::theme(legend.position = "none") else NULL)
  })
  patchwork::wrap_plots(ps, nrow = 1)
}
