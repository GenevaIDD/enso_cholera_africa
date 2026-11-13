# figures_hindcast.R -- ggplot builders for the hindcast-skill artefacts.
#
# Sourced-only library (no numeric prefix; see analysis/utils/README.md section 5).
# Pure: function definitions only, never source()s another module; uses ggplot2,
# dplyr and tidyr through `::` at call time.
# Sourced by: analysis/10_manuscript_figures.R (figure H),
#             notebooks/hindcast_skill.qmd.
#
# Every builder takes tidy tables from the phase-09 artefacts
# (results_hindcast_runs / _occ / _cases, built by utils/hindcast_skill.R) plus
# presentation arguments, so 10 can pass its manuscript theme and palette and
# the notebook its own:
#   labels  named character, scope key -> display label, in display order; must
#           include "africa" when Africa rows are shown (hc_scope_labels())
#   pal     named colours keyed like `labels` (NULL = ggplot default)
#   theme   a ggplot theme object
# Builders return a ggplot; they never save.

#' hc_scope_labels -- display labels for "africa" (first, as in Figures 3-4)
#' plus the region keys.
hc_scope_labels <- function(group_short, groups = names(group_short), africa = "Africa") {
  c(africa = africa, group_short[groups])
}

#' hc_scope_palette -- region colours plus a neutral Africa colour.
hc_scope_palette <- function(group_colours, groups = names(group_colours),
                             africa = "grey25") {
  c(group_colours[groups], africa = africa)
}

.hc_scope <- function(df, labels) {
  df <- dplyr::filter(df, scope %in% names(labels))
  dplyr::mutate(df, scope_lab = factor(unname(labels[scope]), levels = unname(labels)))
}

.hc_colour_scale <- function(pal, labels, aes = "colour", name = NULL) {
  if (is.null(pal)) return(NULL)
  vals <- stats::setNames(unname(pal[names(labels)]), unname(labels))
  if (aes == "fill") ggplot2::scale_fill_manual(values = vals, name = name)
  else ggplot2::scale_colour_manual(values = vals, name = name)
}

.hc_ref_label <- function(x) {
  lab <- c(insample_gid = "In-sample per unit", loo_gid = "Leave-one-window-out per unit",
           loo_gid_fullhist = "Leave-one-window-out per unit, full history",
           regional = "Leave-one-window-out regional rate",
           zero = "Always zero", persistence = "Persistence (1 y earlier)",
           loo_climatology = "LOO climatology")
  out <- unname(lab[x]); ifelse(is.na(out), x, out)
}

.hc_support_label <- function(x) {
  lab <- c(obs_positive_keys = "Unit-seasons with observed cases", exact_gid = "Exact unit",
           atrisk_gids = "At-risk units", atrisk_keys = "At-risk unit-seasons", footprint = "Whole footprint")
  out <- unname(lab[x]); ifelse(is.na(out), x, out)
}

# Lead classes for display: an en dash in a range ("0-3" -> "0–3").
.hc_lead_label <- function(x) gsub("-", "\u2013", x, fixed = TRUE)

# Fill scale of a bar chart keyed by `keys` (labelled by `lab_fun`), the
# headline key opaque and the others at `alpha_other`, the legend keys drawn
# at the same alphas; ggplot's default hues when `pal` is NULL.
.hc_keyed_fill <- function(keys, headline, lab_fun, pal, alpha_other, name, nrow = NULL) {
  al <- ifelse(keys == headline, 1, alpha_other)
  guide <- ggplot2::guide_legend(nrow = nrow, override.aes = list(alpha = al))
  if (is.null(pal)) return(ggplot2::scale_fill_discrete(name = name, guide = guide))
  ggplot2::scale_fill_manual(values = stats::setNames(unname(pal[keys]), lab_fun(keys)), name = name,
                             guide = guide)
}

#' hc_plot_skill_bars
#'
#' Occurrence CRPSS by scope for several references at one grain / rule /
#' lead class. The headline reference is drawn filled, the others outlined.
#' `pal`: colours keyed by reference (NULL: ggplot's default hues);
#' `legend_nrow`: rows of the reference legend (NULL: ggplot's choice).
hc_plot_skill_bars <- function(skill, references, headline, grain = "window",
                               rule = "crps_fair", lead_classes = "all", labels,
                               theme = ggplot2::theme_bw(), pal = NULL, legend_nrow = NULL) {
  references <- intersect(references, unique(skill$reference))
  d <- skill |>
    dplyr::filter(scope_level %in% c("region", "africa"), grain == !!grain,
                  rule == !!rule, lead_class %in% lead_classes,
                  reference %in% references) |>
    .hc_scope(labels) |>
    dplyr::mutate(ref_lab = factor(.hc_ref_label(reference), levels = .hc_ref_label(references)),
                  headline = reference == headline,
                  lead_class = factor(lead_class, levels = lead_classes))
  ggplot2::ggplot(d, ggplot2::aes(scope_lab, crpss, group = ref_lab)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::geom_col(ggplot2::aes(fill = ref_lab, alpha = headline), colour = "grey20",
                      position = ggplot2::position_dodge(width = .8), width = .75,
                      linewidth = .25) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .35), guide = "none") +
    .hc_keyed_fill(references, headline, .hc_ref_label, pal, .35, "Reference", nrow = legend_nrow) +
    ggplot2::facet_wrap(~ lead_class, nrow = 1, labeller = ggplot2::as_labeller(.hc_lead_label)) +
    ggplot2::labs(x = NULL, y = "CRPSS (outbreak count)") +
    theme
}

#' hc_plot_skill_by_lead -- CRPSS against lead class, one line per scope.
hc_plot_skill_by_lead <- function(skill, reference, grain = "window", rule = "crps_fair",
                                  lead_levels, labels, pal = NULL,
                                  theme = ggplot2::theme_bw()) {
  d <- skill |>
    dplyr::filter(scope_level %in% c("region", "africa"), grain == !!grain, rule == !!rule,
                  reference == !!reference, lead_class != "all") |>
    .hc_scope(labels) |>
    dplyr::mutate(lead_class = factor(lead_class, levels = setdiff(lead_levels, "all")))
  ggplot2::ggplot(d, ggplot2::aes(lead_class, crpss, colour = scope_lab, group = scope_lab)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 2) +
    .hc_colour_scale(pal, labels) +
    ggplot2::scale_x_discrete(labels = .hc_lead_label) +
    ggplot2::labs(x = "Lead (months after launch)", y = "CRPSS (outbreak count)",
                  colour = NULL) +
    theme
}

#' hc_plot_cell_scores -- monthly-cell BSS (vs a reference) and AUC by scope.
hc_plot_cell_scores <- function(cell_scores, reference, lead_class = "all", labels,
                                pal = NULL, theme = ggplot2::theme_bw()) {
  d <- cell_scores |>
    dplyr::filter(scope_level %in% c("region", "africa"), reference == !!reference,
                  lead_class == !!lead_class) |>
    .hc_scope(labels) |>
    dplyr::select(scope_lab, BSS = bss, AUC = auc) |>
    tidyr::pivot_longer(c(BSS, AUC), names_to = "metric") |>
    dplyr::mutate(ref0 = ifelse(metric == "AUC", .5, 0))
  ggplot2::ggplot(d, ggplot2::aes(scope_lab, value, fill = scope_lab)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = ref0), colour = "grey30") +
    ggplot2::geom_col(width = .7, show.legend = FALSE) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    .hc_colour_scale(pal, labels, "fill") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme
}

#' hc_plot_case_support
#'
#' Cases track by support: sim/obs median ratio, CRPSS vs a reference and 95 %
#' coverage, one facet per metric, bars per scope coloured by support.
hc_plot_case_support <- function(calibration, skill, supports, headline,
                                 reference = "loo_climatology", rule = "crps_raw",
                                 lead_class = "all", labels, theme = ggplot2::theme_bw(), pal = NULL) {
  supports <- intersect(supports, unique(calibration$support))
  cal <- calibration |>
    dplyr::filter(scope_level %in% c("region", "africa"), lead_class == !!lead_class,
                  support %in% supports)
  sk <- skill |>
    dplyr::filter(scope_level %in% c("region", "africa"), lead_class == !!lead_class,
                  support %in% supports, reference == !!reference, rule == !!rule) |>
    dplyr::select(scope_level, scope, support, crpss)
  d <- dplyr::left_join(cal, sk, by = c("scope_level", "scope", "support")) |>
    .hc_scope(labels) |>
    dplyr::transmute(scope_lab, support, headline = support == headline,
                     `Sim / obs (median)` = sim_obs_ratio,
                     `CRPSS vs LOO climatology` = crpss, `95% coverage` = cov95) |>
    tidyr::pivot_longer(-c(scope_lab, support, headline), names_to = "metric") |>
    dplyr::mutate(support = factor(.hc_support_label(support),
                                   levels = .hc_support_label(supports)),
                  ref = dplyr::case_when(metric == "Sim / obs (median)" ~ 1,
                                         metric == "95% coverage" ~ .95, TRUE ~ 0))
  ggplot2::ggplot(d, ggplot2::aes(scope_lab, value, fill = support)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = ref), colour = "grey30", linetype = "dashed") +
    ggplot2::geom_col(ggplot2::aes(alpha = headline), colour = "grey20", linewidth = .2,
                      position = ggplot2::position_dodge(width = .85), width = .8) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .45), guide = "none") +
    .hc_keyed_fill(supports, headline, .hc_support_label, pal, .45, "Support") +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme
}

#' hc_plot_country_lollipop -- per-country case CRPSS (worst model CRPS first).
hc_plot_country_lollipop <- function(skill, support, reference = "loo_climatology",
                                     rule = "crps_raw", min_runs = 3L,
                                     theme = ggplot2::theme_bw()) {
  d <- skill |>
    dplyr::filter(scope_level == "country", support == !!support, lead_class == "all",
                  reference == !!reference, rule == !!rule, n_runs >= min_runs,
                  is.finite(crpss)) |>
    dplyr::arrange(dplyr::desc(mean_model * n_runs)) |>
    dplyr::mutate(scope = factor(scope, levels = rev(scope)),
                  sign = ifelse(crpss >= 0, "beats reference", "worse than reference"))
  ggplot2::ggplot(d, ggplot2::aes(crpss, scope, colour = sign)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey30") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = crpss, yend = scope)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(`beats reference` = "#1b7837",
                                            `worse than reference` = "#b2182b"), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(max(-3, min(d$crpss, 0)), 1)) +
    ggplot2::labs(x = "Case CRPSS vs LOO climatology", y = NULL,
                  caption = "Ordered by total model CRPS (largest first)") +
    theme
}

#' hc_plot_cases_interval -- per-run window totals: obs vs ensemble 50/95 %.
#' The observed point is coloured by where it falls against the 95 % interval
#' (below / inside / above), read from obs vs q025 / q975 (the artefact's `pos`
#' codes the forecast's error instead: "over" = observation below). Region
#' grain has one row per group run; the artefact holds no Africa window totals.
hc_plot_cases_interval <- function(run_totals, support, grain = "region",
                                   lead_class = "all", labels,
                                   theme = ggplot2::theme_bw()) {
  pos_col <- c(below = "#2166ac", inside = "black", above = "#b2182b")
  pos_lab <- c(below = "below the 95% interval", inside = "inside the 95% interval",
               above = "above the 95% interval")
  d <- run_totals |>
    dplyr::filter(support == !!support, grain == !!grain, lead_class == !!lead_class) |>
    dplyr::mutate(obs_pos = factor(dplyr::case_when(obs < q025 ~ "below", obs > q975 ~ "above",
                                                    TRUE ~ "inside"),
                                   levels = names(pos_col)))
  d <- if (grain == "region") .hc_scope(d, labels) else dplyr::mutate(d, scope_lab = scope)
  ggplot2::ggplot(d, ggplot2::aes(factor(launch_year))) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q025, ymax = q975), colour = "grey70",
                            linewidth = 1) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q25, ymax = q75), colour = "grey35",
                            linewidth = 2.2) +
    ggplot2::geom_point(ggplot2::aes(y = med), shape = 3, size = 2) +
    # show.legend = TRUE keys every position, also one no window falls in
    ggplot2::geom_point(ggplot2::aes(y = obs, colour = obs_pos), size = 2.2, show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = pos_col, labels = pos_lab, drop = FALSE,
                                 name = "Observed total") +
    # 1-3-10 breaks (the default log1p breaks collided near the top of a
    # panel); rounded before labelling (log1p's inverse returns 999.99..., which
    # printed "1 000" instead of "1K")
    ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, outer(c(1, 3), 10^(0:6))),
                                labels = function(x)
                                  scales::label_number(scale_cut = scales::cut_short_scale())(round(x))) +
    ggplot2::facet_wrap(~ scope_lab, scales = "free_y") +
    ggplot2::labs(x = "Launch year", y = "Cases over the window (log1p)") +
    theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5))
}

#' hc_plot_monthly -- monthly cases: ensemble median + 95 % band vs observed.
hc_plot_monthly <- function(monthly, support, grain = "region", scopes = NULL,
                            labels = NULL, theme = ggplot2::theme_bw()) {
  d <- monthly |> dplyr::filter(support == !!support, grain == !!grain)
  if (!is.null(scopes)) d <- dplyr::filter(d, scope %in% scopes)
  d <- if (grain == "region" && !is.null(labels)) .hc_scope(d, labels)
       else dplyr::mutate(d, scope_lab = scope)
  ggplot2::ggplot(d, ggplot2::aes(month)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, group = launch_year),
                         fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = med, group = launch_year), colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(y = obs), colour = "#b2182b", size = .8) +
    ggplot2::facet_wrap(~ scope_lab, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Cases per month") +
    theme
}

.hc_hist <- function(d, x_lab, theme) {
  ggplot2::ggplot(d, ggplot2::aes(xmin = lo, xmax = hi, ymin = 0, ymax = density)) +
    ggplot2::geom_rect(fill = "grey65", colour = "white") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    ggplot2::facet_wrap(~ scope_lab) +
    ggplot2::labs(x = x_lab, y = "Density (1 = uniform)") +
    theme
}

#' hc_plot_pit_hist -- pooled cell-level PIT histogram by scope.
hc_plot_pit_hist <- function(pit_hist, labels, theme = ggplot2::theme_bw())
  .hc_hist(.hc_scope(pit_hist, labels), "PIT", theme)

#' hc_plot_rank_hist -- randomised verification-rank histogram by scope.
hc_plot_rank_hist <- function(rank_hist, grain = "window", lead_class = "all", labels,
                              theme = ggplot2::theme_bw())
  .hc_hist(.hc_scope(dplyr::filter(rank_hist, grain == !!grain, lead_class == !!lead_class),
                     labels), "Normalised rank of the observation", theme)

#' hc_plot_reliability -- reliability diagram (decile bins of p_hat).
# Ticks 0, s, 2s, ... up to the upper limit, s = the limit's half rounded to one
# significant digit (two or three ticks on every facet of a free-scale grid).
.hc_axis_breaks <- function(l) { s <- signif(l[2] / 2, 1); if (!is.finite(s) || s <= 0) return(0); seq(0, l[2], by = s) }
hc_plot_reliability <- function(reliability, labels, theme = ggplot2::theme_bw()) {
  d <- .hc_scope(reliability, labels)
  ggplot2::ggplot(d, ggplot2::aes(mean_p, obs_freq)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = "grey35", alpha = .8) +
    ggplot2::facet_wrap(~ scope_lab, scales = "free") +
    ggplot2::scale_size_continuous(name = "Unit-months",
      labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    # two or three round ticks per small facet from 0 (five 0.0x labels ran
    # into each other; extended breaks just outside the expanded limits left
    # one unlabelled tick on some facets)
    ggplot2::scale_x_continuous(limits = c(0, NA), breaks = .hc_axis_breaks,
                                labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE)) +
    ggplot2::scale_y_continuous(limits = c(0, NA), breaks = .hc_axis_breaks,
                                labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE)) +
    ggplot2::labs(x = "Mean forecast probability (bin)", y = "Observed onset frequency") +
    theme
}

#' hc_plot_roc -- ROC curves by scope, AUC in the legend.
hc_plot_roc <- function(roc, auc_tbl, labels, pal = NULL, theme = ggplot2::theme_bw()) {
  d <- .hc_scope(roc, labels)
  a <- .hc_scope(auc_tbl, labels)
  leg <- stats::setNames(sprintf("%s (AUC %.2f)", a$scope_lab, a$auc), as.character(a$scope_lab))
  ggplot2::ggplot(d, ggplot2::aes(fpr, tpr, colour = scope_lab)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_step(direction = "vh") +
    (if (is.null(pal)) ggplot2::scale_colour_discrete(labels = leg, name = NULL)
     else ggplot2::scale_colour_manual(values = stats::setNames(unname(pal[names(labels)]),
                                                                unname(labels)),
                                       labels = leg, name = NULL)) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False-positive rate", y = "True-positive rate") +
    theme
}

#' hc_plot_coverage_curve -- empirical vs nominal central-interval coverage.
hc_plot_coverage_curve <- function(coverage_curve, labels, pal = NULL,
                                   theme = ggplot2::theme_bw()) {
  d <- .hc_scope(coverage_curve, labels)
  ggplot2::ggplot(d, ggplot2::aes(width, coverage, colour = scope_lab)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
    .hc_colour_scale(pal, labels) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Nominal interval width", y = "Empirical coverage", colour = NULL) +
    theme
}

#' hc_plot_interval_score -- mean Winkler score across widths.
hc_plot_interval_score <- function(interval_score, grain = "window", lead_class = "all",
                                   labels, pal = NULL, theme = ggplot2::theme_bw()) {
  d <- .hc_scope(dplyr::filter(interval_score, grain == !!grain,
                               lead_class == !!lead_class), labels)
  ggplot2::ggplot(d, ggplot2::aes(width, IS, colour = scope_lab)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.6) +
    .hc_colour_scale(pal, labels) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Nominal interval width", y = "Mean interval (Winkler) score",
                  colour = NULL) +
    theme
}

#' hc_plot_season -- seasonal cycle of onset frequency, obs vs simulated band.
hc_plot_season <- function(season_cycle, labels, theme = ggplot2::theme_bw()) {
  d <- .hc_scope(season_cycle, labels)
  ggplot2::ggplot(d, ggplot2::aes(month_of_year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sim_lo, ymax = sim_hi), fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = sim_mean), colour = "grey30") +
    ggplot2::geom_line(ggplot2::aes(y = freq_obs), colour = "#b2182b") +
    ggplot2::geom_point(ggplot2::aes(y = freq_obs), colour = "#b2182b", size = 1.2) +
    ggplot2::scale_x_continuous(breaks = 1:12, labels = month.abb) +
    ggplot2::facet_wrap(~ scope_lab) +
    ggplot2::labs(x = NULL, y = "Monthly onset frequency",
                  caption = "Red = observed; grey = simulated mean and 95 % band") +
    theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5))
}
