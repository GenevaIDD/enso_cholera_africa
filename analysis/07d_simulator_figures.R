# 07d_simulator_figures.R — diagnostic figures for the cholera outbreak simulator
#
# Reads validation/metrics.rds produced by 07c_simulator_postprocess.R and produces
# six diagnostic figures:
#   1. pit_calibration.png       — PIT histogram (calibration check)
#   2. seasonal_frequency.png    — monthly outbreak frequency, obs vs sim by country
#   3. crps_annual.png           — annual CRPS distribution
#   4. ar_run_lengths.png        — simulated run/gap length histograms
#   5. contagion_conditionals.png — P(outbreak | lag status) for three spillover lags
#   6. crps_total_cases.png      — annual total cases, obs vs sim (sim = episode
#                                  sums of the mon track's cases_month)
#
# Usage:
#   Rscript analysis/07d_simulator_figures.R -C analysis/configs/sim_smoketest/sim_smoke.yml
#   Rscript analysis/07d_simulator_figures.R -C ... --fig-dir figures/sim/myrun/
#   Rscript analysis/07d_simulator_figures.R -C ... --engine spatpomp --if-validated
#     (exit 0 with a message when the run has no validation metrics: a forward
#      forecast / counterfactual run has no observations to validate against)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(cowplot)
  library(lubridate)
  library(qs2)
})

# ---------------------------------------------------------------------------
# fig_pit_calibration — PIT histogram
# ---------------------------------------------------------------------------

#' Histogram of randomised PIT values with U[0,1] reference line.
#'
#' @param pit tibble(gid, month, p_hat, y_obs, pit) from compute_pit_binom().
#' @return ggplot object.
fig_pit_calibration <- function(pit) {
  if (is.null(pit) || nrow(pit) == 0L)
    return(ggplot() + labs(title = "PIT calibration (no data)") + theme_minimal())

  n_cells <- nrow(pit)
  bin_w   <- 0.05
  ggplot(pit, aes(x = pit)) +
    geom_histogram(binwidth = bin_w, fill = "steelblue3", colour = "white", alpha = 0.85) +
    geom_hline(yintercept = n_cells * bin_w, linetype = "dashed",
               colour = "tomato2", linewidth = 0.7) +
    labs(
      title    = "PIT calibration",
      subtitle = sprintf(
        "n=%d cells   mean=%.3f  (U[0,1] ≈ 0.500)   sd=%.3f  (U[0,1] ≈ 0.289)",
        n_cells,
        mean(pit$pit, na.rm = TRUE),
        stats::sd(pit$pit, na.rm = TRUE)),
      x = "Probability integral transform",
      y = "Count",
      caption = "Dashed line: expected count under perfect calibration (U[0,1])"
    ) +
    theme_minimal(base_size = 11)
}

# ---------------------------------------------------------------------------
# fig_seasonal_frequency — monthly outbreak frequency by country
# ---------------------------------------------------------------------------

#' Observed vs simulated monthly outbreak frequency, faceted by country.
#'
#' Shows up to 16 countries ordered by peak observed frequency.
#'
#' @param monthly_freq tibble(country, month_of_year, freq_obs, freq_sim_mean,
#'   freq_sim_lo, freq_sim_hi) from compute_monthly_freq().
#' @return ggplot object.
fig_seasonal_frequency <- function(monthly_freq) {
  if (is.null(monthly_freq) || nrow(monthly_freq) == 0L)
    return(ggplot() + labs(title = "Seasonal frequency (no data)") + theme_minimal())

  top_ctry <- monthly_freq %>%
    dplyr::filter(!is.na(freq_obs)) %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(peak = max(freq_obs, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(peak)) %>%
    dplyr::pull(country)

  df <- monthly_freq %>%
    dplyr::filter(country %in% top_ctry) %>%
    dplyr::mutate(country = factor(country, levels = top_ctry))

  ggplot(df, aes(x = month_of_year)) +
    geom_ribbon(aes(ymin = freq_sim_q025, ymax = freq_sim_q975),
                fill = "steelblue3", alpha = 0.18, na.rm = TRUE) +
    geom_ribbon(aes(ymin = freq_sim_q25, ymax = freq_sim_q75),
                fill = "steelblue3", alpha = 0.35, na.rm = TRUE) +
    geom_line(aes(y = freq_sim_mean), colour = "steelblue4", linewidth = 0.6,
              na.rm = TRUE) +
    geom_line(aes(y = freq_obs), colour = "tomato2", linewidth = 0.75,
              na.rm = TRUE) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, minor_breaks = NULL) +
    facet_wrap(~country, ncol = 4, scales = "free_y") +
    labs(
      title    = "Monthly outbreak frequency — observed (red) vs simulated (blue, 50%/95% CI ribbons)",
      subtitle = "All countries, sorted by peak observed frequency",
      x = "Month", y = "Mean frequency"
    ) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}

# ---------------------------------------------------------------------------
# fig_crps_annual — annual CRPS distribution + multi-level count comparison
# ---------------------------------------------------------------------------

#' Annual outbreak-count CRPS violin (top) + sim-vs-obs count panels by
#' Africa / AFRO-region / country (bottom, when crps_count_by_geo is supplied).
#'
#' The violin shows the distribution of per-(gid,year) CRPS scores, with the
#' observed outbreak count per year annotated below the x-axis. The count
#' panels mirror the structure of fig_crps_total_cases: simulated 80% CI
#' ribbon + mean line + observed red line.
#'
#' @param crps tibble(gid, year, n_sim, y_obs, crps_count) from
#'   compute_crps_annual().
#' @param crps_count_by_geo optional named list(africa, region, country) from
#'   compute_annual_count_by_geo(). When NULL only the violin is returned.
#' @return ggplot or cowplot composite.
fig_crps_annual <- function(crps, crps_count_by_geo = NULL) {
  if (is.null(crps) || nrow(crps) == 0L)
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "Annual CRPS (no data)") +
             ggplot2::theme_minimal())

  df <- crps %>% dplyr::mutate(year = factor(year))

  # Per-year observed outbreak count (summed across gids) for annotation
  obs_per_year <- crps %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(total_obs = sum(y_obs, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(year = factor(year),
                  label = sprintf("n=%d", round(total_obs)))

  p_violin <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = crps_count)) +
    ggplot2::geom_violin(fill = "steelblue3", alpha = 0.65) +
    ggplot2::geom_boxplot(width = 0.08, fill = "white", alpha = 0.7,
                          outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_text(
      data  = obs_per_year,
      ggplot2::aes(x = year, y = -Inf, label = label),
      vjust = 1.8, size = 2.5, colour = "tomato3") +
    ggplot2::labs(
      title    = "Annual outbreak-count CRPS by year",
      subtitle = sprintf(
        "n=%d (gid, year) pairs   median=%.3f   mean=%.3f   (red=observed outbreaks/year)",
        nrow(crps),
        stats::median(crps$crps_count, na.rm = TRUE),
        mean(crps$crps_count, na.rm = TRUE)),
      x = "Year", y = "CRPS"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = ggplot2::margin(b = 20))

  if (is.null(crps_count_by_geo)) return(p_violin)

  # ── Helper: one panel per geographic level ──────────────────────────────
  make_count_panel <- function(tbl, level_label) {
    if (is.null(tbl) || nrow(tbl) == 0L)
      return(ggplot2::ggplot() +
               ggplot2::labs(title = paste(level_label, "(no data)")) +
               ggplot2::theme_minimal())

    n_geo <- dplyr::n_distinct(tbl$geography)
    ncols  <- if (n_geo <= 4) n_geo else if (n_geo <= 9) 3L else 4L

    geo_crps <- tbl %>%
      dplyr::group_by(geography) %>%
      dplyr::summarise(mean_crps = mean(crps_count, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(geo_label = sprintf("%s\n(mean CRPS=%.1f)", geography, mean_crps))
    tbl <- dplyr::left_join(tbl, geo_crps, by = "geography") %>%
      dplyr::mutate(geo_label = factor(geo_label))

    ggplot2::ggplot(tbl, ggplot2::aes(x = year)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = total_sim_q025, ymax = total_sim_q975),
        fill = "steelblue3", alpha = 0.18, na.rm = TRUE) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = total_sim_q25, ymax = total_sim_q75),
        fill = "steelblue3", alpha = 0.35, na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_sim_mean), colour = "steelblue4",
        linewidth = 0.6, na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_sim_unmasked_mean), colour = "steelblue4",
        linewidth = 0.5, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_obs), colour = "tomato2",
        linewidth = 0.75, na.rm = TRUE) +
      ggplot2::facet_wrap(~geo_label, ncol = ncols, scales = "free_y") +
      ggplot2::labs(
        title    = sprintf("Annual outbreak count — %s level", level_label),
        subtitle = "Red=observed  Blue=sim mean (50%/95% CI ribbons)  Dashed=unmasked sim mean",
        x = "Year", y = "Outbreak count"
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  p_africa  <- make_count_panel(crps_count_by_geo$africa,  "All-Africa")
  p_region  <- make_count_panel(crps_count_by_geo$region,  "AFRO region")
  p_country <- make_count_panel(crps_count_by_geo$country, "Country")

  cowplot::plot_grid(p_violin, p_africa, p_region, p_country,
                     ncol = 1, rel_heights = c(1.2, 0.8, 1.2, 3))
}

# ---------------------------------------------------------------------------
# fig_ar_run_lengths — run and gap length histograms
# ---------------------------------------------------------------------------

#' Bar charts of simulated run lengths (Y=1) and gap lengths (Y=0).
#'
#' @param ar_metrics list(run_length_table, gap_length_table, transition) from
#'   compute_ar_metrics().
#' @return ggplot object.
fig_ar_run_lengths <- function(ar_metrics) {
  if (is.null(ar_metrics))
    return(ggplot() + labs(title = "AR run/gap lengths (no data)") + theme_minimal())

  tbl_to_df <- function(tbl, label) {
    if (length(tbl) == 0L)
      return(tibble::tibble(length = integer(), count = integer(), type = character()))
    tibble::tibble(
      length = as.integer(names(tbl)),
      count  = as.integer(tbl),
      type   = label
    )
  }

  df <- dplyr::bind_rows(
    tbl_to_df(ar_metrics$run_length_table, "Run lengths (Y=1)"),
    tbl_to_df(ar_metrics$gap_length_table, "Gap lengths (Y=0)")
  )

  if (nrow(df) == 0L)
    return(ggplot() + labs(title = "AR run/gap lengths (no data)") + theme_minimal())

  ggplot(df, aes(x = length, y = count, fill = type)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("Run lengths (Y=1)" = "steelblue3",
                                  "Gap lengths (Y=0)" = "indianred3")) +
    facet_wrap(~type, scales = "free", ncol = 2) +
    labs(
      title = "Simulated run and gap lengths",
      x = "Length (months)", y = "Count"
    ) +
    theme_minimal(base_size = 11)
}

# ---------------------------------------------------------------------------
# fig_contagion_conditionals — P(outbreak | lag status) for three lags
# ---------------------------------------------------------------------------

#' Grouped bar chart of P(Y=1 | lag=1), P(Y=1 | lag=0), and unconditional P(Y)
#' for each of the three spillover lag types (ADM1, country, TU).
#'
#' @param contagion list(adm1_cond, country_cond, tu_cond, unconditional) from
#'   compute_contagion_conditionals().
#' @return ggplot object.
fig_contagion_conditionals <- function(contagion) {
  if (is.null(contagion))
    return(ggplot() + labs(title = "Contagion conditionals (no data)") + theme_minimal())

  rows <- list()
  lag_map <- c(adm1_cond    = "ADM1 lag",
               country_cond = "Country lag",
               tu_cond      = "TU lag")
  for (nm in names(lag_map)) {
    cond <- contagion[[nm]]
    if (!is.null(cond) && !is.null(cond$lag_active)) {
      rows[[nm]] <- tibble::tibble(
        lag_type  = lag_map[[nm]],
        condition = c("lag = 1", "lag = 0", "unconditional"),
        prob      = c(cond$lag_active    %||% NA_real_,
                      cond$lag_inactive  %||% NA_real_,
                      cond$unconditional %||% NA_real_)
      )
    }
  }

  df <- dplyr::bind_rows(rows)
  if (nrow(df) == 0L)
    return(ggplot() + labs(title = "Contagion conditionals (no data)") + theme_minimal())

  df <- df %>%
    dplyr::mutate(
      lag_type  = factor(lag_type,  levels = c("ADM1 lag", "Country lag", "TU lag")),
      condition = factor(condition, levels = c("lag = 1", "lag = 0", "unconditional"))
    )

  uncond_label <- sprintf("Overall P(Y) = %.4f", contagion$unconditional %||% NA_real_)

  ggplot(df, aes(x = lag_type, y = prob, fill = condition)) +
    geom_col(position = "dodge") +
    scale_fill_manual(
      values = c("lag = 1" = "steelblue4", "lag = 0" = "steelblue1",
                 "unconditional" = "grey55"),
      name = NULL
    ) +
    labs(
      title    = "Contagion conditionals — P(outbreak | lag status)",
      subtitle = uncond_label,
      x = NULL, y = "P(outbreak)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# fig_crps_total_cases — annual total-case CRPS at three geographic levels
# ---------------------------------------------------------------------------

#' Line + ribbon plots of annual total suspected cases: observed vs simulated
#' ensemble, faceted by geographic level (Africa, AFRO region, country).
#'
#' Each panel shows the simulated 80% CI ribbon (masked), the masked sim mean,
#' the unmasked sim mean (dashed), and the observed total (red line). The CRPS
#' value is annotated in the strip label.
#'
#' @param crps_total_cases named list(africa, region, country) from
#'   compute_crps_total_cases().
#' @return ggplot object (list of three faceted panels via cowplot).
fig_crps_total_cases <- function(crps_total_cases) {
  if (is.null(crps_total_cases))
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "Total-case CRPS (no data — metric 2b skipped)") +
             ggplot2::theme_minimal())

  make_panel <- function(tbl, level_label) {
    if (is.null(tbl) || nrow(tbl) == 0L)
      return(ggplot2::ggplot() +
               ggplot2::labs(title = paste(level_label, "(no data)")) +
               ggplot2::theme_minimal())

    # Determine number of geographies to pick sensible facet count
    n_geo <- dplyr::n_distinct(tbl$geography)
    ncols <- if (n_geo <= 4) n_geo else if (n_geo <= 9) 3L else 4L

    tbl <- tbl %>%
      dplyr::mutate(
        geography = as.character(geography),
        crps_label = sprintf("CRPS=%.0f", crps_cases)
      )

    # Per-geography mean CRPS for strip annotation
    geo_crps <- tbl %>%
      dplyr::group_by(geography) %>%
      dplyr::summarise(mean_crps = mean(crps_cases, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(geo_label = sprintf("%s\n(mean CRPS=%.0f)", geography, mean_crps))
    tbl <- dplyr::left_join(tbl, geo_crps, by = "geography") %>%
      dplyr::mutate(geo_label = factor(geo_label))

    p <- ggplot2::ggplot(tbl, ggplot2::aes(x = year)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = total_sim_q025, ymax = total_sim_q975),
        fill = "steelblue3", alpha = 0.18, na.rm = TRUE) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = total_sim_q25, ymax = total_sim_q75),
        fill = "steelblue3", alpha = 0.35, na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_sim_mean), colour = "steelblue4",
        linewidth = 0.6, na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_sim_unmasked_mean), colour = "steelblue4",
        linewidth = 0.5, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = total_obs), colour = "tomato2",
        linewidth = 0.75, na.rm = TRUE) +
      ggplot2::facet_wrap(~geo_label, ncol = ncols, scales = "free_y") +
      ggplot2::labs(
        title    = sprintf("Annual total cases — %s level", level_label),
        subtitle = "Red=observed  Blue=sim mean (50%/95% CI ribbons)  Dashed=unmasked sim mean",
        x = "Year", y = "Total suspected cases"
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    p
  }

  p_africa  <- make_panel(crps_total_cases$africa,  "All-Africa")
  p_region  <- make_panel(crps_total_cases$region,  "AFRO region")
  p_country <- make_panel(crps_total_cases$country, "Country")

  cowplot::plot_grid(p_africa, p_region, p_country,
                     ncol = 1, rel_heights = c(1, 1.5, 3))
}

.hindcast_band_plot <- function(band, series, y_lab, title, subtitle) {
  launch <- as.Date(series$meta$launch)
  ggplot(band, aes(x = month)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "steelblue3", alpha = 0.18) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "steelblue3", alpha = 0.30) +
    geom_line(aes(y = med), colour = "steelblue4", linewidth = 0.8) +
    geom_line(aes(y = obs), colour = "tomato2", linewidth = 0.8) +
    geom_point(aes(y = obs), colour = "tomato2", size = 1.6) +
    geom_vline(xintercept = launch, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    labs(title = title, subtitle = subtitle, x = NULL, y = y_lab,
         caption = "Blue: simulated forecast (line = median, bands = 50% / 95% across replicates). Red: observed. Dashed: launch month.") +
    theme_minimal(base_size = 11)
}

#' Monthly outbreak occurrence: simulated forecast band vs observed.
fig_hindcast_occurrence <- function(series) {
  if (is.null(series)) return(ggplot() + labs(title = "Hindcast occurrence (no data)") + theme_minimal())
  m <- series$meta
  .hindcast_band_plot(
    series$occ_band, series,
    y_lab = "Active admin units",
    title = sprintf("Hindcast vs observed — monthly outbreak occurrence (%s, launch %s)",
                    m$group, m$launch_year),
    subtitle = sprintf("%d replicates · %d admin units · %d scored (gid,month) cells",
                       m$n_rep, m$n_gid, m$n_cells))
}

#' Monthly suspected cases: simulated forecast band vs observed (at onset month).
fig_hindcast_cases <- function(series) {
  if (is.null(series) || !isTRUE(series$have_cases))
    return(ggplot() + labs(title = "Hindcast cases (monthly-cases data unavailable)") + theme_minimal())
  m <- series$meta
  .hindcast_band_plot(
    series$cases_band, series,
    y_lab = "Suspected cases / month",
    title = sprintf("Hindcast vs observed — monthly suspected cases (%s, launch %s)",
                    m$group, m$launch_year),
    subtitle = sprintf("%d replicates · %d admin units · %d scored (gid,month) cells",
                       m$n_rep, m$n_gid, m$n_cells))
}

#' Histogram of cumulative predicted cases across replicates vs observed total.
fig_hindcast_cumulative_cases <- function(series) {
  if (is.null(series) || !isTRUE(series$have_cases) || !length(series$cum_sim))
    return(ggplot() + labs(title = "Hindcast cumulative cases (monthly-cases data unavailable)") + theme_minimal())
  m   <- series$meta
  cs  <- series$cum_sim
  obs <- series$obs_cum_cases
  med <- stats::median(cs, na.rm = TRUE)
  lo  <- stats::quantile(cs, 0.025, names = FALSE, na.rm = TRUE)
  hi  <- stats::quantile(cs, 0.975, names = FALSE, na.rm = TRUE)
  pct <- 100 * mean(cs <= obs, na.rm = TRUE)   # observed's percentile in predictive
  df  <- tibble::tibble(cum = cs)
  ggplot(df, aes(x = cum)) +
    geom_histogram(bins = 40, fill = "steelblue3", colour = "white", alpha = 0.85) +
    geom_vline(xintercept = med, linetype = "dotted",
               colour = "steelblue4", linewidth = 0.7) +
    geom_vline(xintercept = obs, linetype = "dashed",
               colour = "tomato2", linewidth = 0.9) +
    labs(
      title    = sprintf("Cumulative cases: predicted distribution vs observed (%s, launch %s)",
                         m$group, m$launch_year),
      subtitle = sprintf("Predicted median %s [95%% %s, %s] · observed %s (predictive percentile %.0f%%)",
                         format(round(med), big.mark = ","),
                         format(round(lo), big.mark = ","),
                         format(round(hi), big.mark = ","),
                         format(round(obs), big.mark = ","), pct),
      x = "Cumulative suspected cases over scored window",
      y = "Replicates",
      caption = "Dashed red: observed total. Dotted blue: predicted median. Over the same coverage-masked cells.") +
    theme_minimal(base_size = 11)
}

#' Write the three hindcast figures for one run.
make_hindcast_figures <- function(series, fig_dir) {
  if (is.null(series)) return(invisible(NULL))
  figs <- list(
    hindcast_occurrence       = fig_hindcast_occurrence(series),
    hindcast_cases            = fig_hindcast_cases(series),
    hindcast_cumulative_cases = fig_hindcast_cumulative_cases(series)
  )
  widths  <- c(hindcast_occurrence = 10, hindcast_cases = 10, hindcast_cumulative_cases = 8)
  heights <- c(hindcast_occurrence = 5,  hindcast_cases = 5,  hindcast_cumulative_cases = 5)
  save <- figureSaver(fig_dir, dpi = 150)
  for (nm in names(figs)) save(figs[[nm]], nm, width = widths[[nm]], height = heights[[nm]])
  invisible(figs)
}

#' Manifest-driven hindcast figure orchestrator (mirrors run_hindcast_skill()).
#'
#' Loops the manifest's group x launch_year entries; for each run with an
#' existing merged/monthly_summary it writes the three hindcast figures under
#' fig_root/hindcast_<group>_<year>/. Respects the derive/plot split: it does
#' NOT aggregate partitions — runs without a merged summary are skipped (run
#' 07c_simulator_postprocess.R --aggregate / --hindcast first).
run_hindcast_figures <- function(manifest_path, obs_panel_path, option_list,
                                 fig_root, workers = 1L) {
  man <- read_hindcast_manifest(manifest_path, obs_panel_path)   # utils/sim_io.R

  figure_run <- function(run) {
    g <- run$group; yr <- run$year; ye <- run$entry
    monthly_data_path <- makeMonthlyDatFile(run$opt)
    series <- compute_hindcast_series(
      merged_dir = run$merged_dir, obs_panel_path = obs_panel_path,
      monthly_data_path = monthly_data_path,
      launch = ye$launch, horizon_end = ye$horizon_end,
      out_fmt = run$out_fmt, group_label = g,
      launch_year = run$launch_year)
    fig_dir <- file.path(fig_root, sprintf("hindcast_%s_%s", g, yr))
    make_hindcast_figures(series, fig_dir)
    message("[hindcast-fig]   figures -> ", fig_dir)
    NULL
  }

  it <- for_each_hindcast_run(man, option_list, "[hindcast-fig]", figure_run,
                              announce = "after", workers = workers)
  n_fig <- it$n_ok; n_skip <- it$n_skip
  message(sprintf("\n[hindcast-fig] ---------------------------------------------"))
  message(sprintf("  runs figured : %d   skipped: %d", n_fig, n_skip))
  message(sprintf("  fig root     : %s", fig_root))
  invisible(list(n_fig = n_fig, n_skip = n_skip))
}

# ---------------------------------------------------------------------------
# make_sim_figures — orchestrate all six figures
# ---------------------------------------------------------------------------

#' Generate all diagnostic figures and save as PNGs.
#'
#' @param metrics named list returned by run_validation() — keys pit,
#'   monthly_freq, crps, ar_metrics, contagion, crps_total_cases.
#' @param fig_dir character; directory to write PNG files (created if absent).
#' @return invisible named list of ggplot objects.
#' @details Depends on `figureSaver()` (`analysis/utils/figures_io.R`) being
#'   loaded by the caller -- the CLI block below sources it directly;
#'   `analysis/tests/test_simulator_validation.R`, which calls this function
#'   after `source()`-ing this file, sources it separately for the same
#'   reason (`figures_io.R` is deliberately not part of the 17-module
#'   `utils.R` facade -- see that file's header).
make_sim_figures <- function(metrics, fig_dir) {
  figs <- list(
    pit_calibration        = fig_pit_calibration(metrics$pit),
    seasonal_frequency     = fig_seasonal_frequency(metrics$monthly_freq),
    crps_annual            = fig_crps_annual(metrics$crps, metrics$crps_count_by_geo),
    ar_run_lengths         = fig_ar_run_lengths(metrics$ar_metrics),
    contagion_conditionals = fig_contagion_conditionals(metrics$contagion),
    crps_total_cases       = fig_crps_total_cases(metrics$crps_total_cases)
  )

  widths  <- c(pit_calibration = 7, seasonal_frequency = 12,
               crps_annual     = 10, ar_run_lengths     = 10,
               contagion_conditionals = 7, crps_total_cases = 10)
  heights <- c(pit_calibration = 5, seasonal_frequency = 16,
               crps_annual     = 14, ar_run_lengths     =  5,
               contagion_conditionals = 5, crps_total_cases = 14)

  # figureSaver() (utils/figures_io.R) creates fig_dir on first
  # write and replaces the bare per-figure ggsave() + message() pair below.
  save <- figureSaver(fig_dir, dpi = 150)
  for (nm in names(figs)) {
    save(figs[[nm]], nm, width = widths[[nm]], height = heights[[nm]])
  }

  invisible(figs)
}

# ---------------------------------------------------------------------------
# CLI block
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({ library(optparse); library(yaml) })
  source("analysis/utils/config.R")
  source("analysis/utils/formula.R")               # model-term assembly
  source("analysis/utils/base.R")
  source("analysis/utils/gid.R")
  source("analysis/utils/paths.R")
  source("analysis/utils/clim_compile.R")
  source("analysis/utils/clim_db.R")
  source("analysis/utils/occ_data.R")
  source("analysis/utils/siz_data.R")      # occ_data.R dependency (utils/README.md §5); no size figure is drawn here
  source("analysis/utils/provenance.R")
  source("analysis/utils/sim_io.R")
  source("analysis/utils/figures_io.R")
  source("analysis/utils/anchor.R")
  # load_observed_from_panel() / compute_hindcast_series(). Only
  # those two are called here; they need gid.R + sim_io.R, sourced above.
  source("analysis/utils/obs_sim_compare.R")   # crps_ens()/region_key_of() -- used by sim_metrics.R
  source("analysis/utils/sim_metrics.R")

  option_list <- list(
    make_option(c("--fig-dir"), default = NULL, type = "character", dest = "fig_dir",
                help = "Override output directory for PNG figures."),
    make_option("--engine", default = "v1", type = "character",
                help = "Output layout to resolve figure inputs against: 'v1' (07b_simulator_v1.R, generated_data/sim/<run_tag>/[run_id/], default) or 'spatpomp' (07b_simulator_spatpomp.R, generated_data/sim_spatpomp/<run_tag>/ — flat, NOT run_id-scoped; see makeSpatpompRunFile() in utils.R)."),
    make_option("--workers", default = 1L, type = "integer",
                help = "--hindcast: forked workers drawing runs in parallel (serial when 1) [default: %default]"),
    make_option("--hindcast", action = "store_true", default = FALSE,
                help = "Hindcast-figures mode: loop the --manifest, writing observed-vs-forecast occurrence/cases/cumulative-cases figures per run (no config needed)."),
    make_option("--manifest", default = NULL, type = "character",
                help = "hindcast_horizons.yml manifest (required with --hindcast)."),
    make_option("--obs-panel", default = NULL, type = "character", dest = "obs_panel",
                help = "Observed survival panel RDS (required with --hindcast)."),
    make_option("--if-validated", action = "store_true", default = FALSE, dest = "if_validated",
                help = paste("Skip (exit 0) when the run has no validation metrics instead of",
                             "stopping: forward forecast / counterfactual runs are never",
                             "validated, and every figure here is a validation figure.")),
    make_option("--hindcast-fig-dir", default = "figures/hindcast", type = "character",
                dest = "hindcast_fig_dir",
                help = "Root output dir for --hindcast figures [default %default].")
  )
  opt <- make_options_from_config(option_list = option_list,
                                   enforce_options = character(0),
                                   model = NULL)

  # --hindcast is manifest-driven and self-contained (no single simulator: block).
  if (isTRUE(opt$hindcast)) {
    run_hindcast_figures(manifest_path  = opt$manifest,
                         obs_panel_path = opt$obs_panel,
                         option_list    = option_list,
                         fig_root       = opt$hindcast_fig_dir,
                         workers        = opt$workers %||% 1L)
    quit(save = "no", status = 0)
  }

  if (is.null(opt$simulator))
    stop("Config must contain a `simulator:` block.")

  # Path dispatch mirrors 07c_simulator_postprocess.R's --engine handling
  # (see makeSpatpompRunFile() in utils.R) so figures resolve against the
  # same layout --aggregate/--validate wrote to.
  engine    <- opt$engine %||% "v1"
  dirs      <- makeSimEngineDirs(opt, engine)   # utils/paths.R (validates --engine)
  run_dir   <- dirs$run_dir
  valid_dir <- dirs$valid_dir

  # Mirror run_validation()'s metrics serialization: list-valued metrics use
  # rds/qs only (parquet configs fall back to qs).
  metrics_fmt  <- sim_metrics_fmt(opt$simulator$output$format)
  metrics_path <- file.path(valid_dir, paste0("metrics", sim_io_ext(metrics_fmt)))
  if (!file.exists(metrics_path) && isTRUE(opt$if_validated)) {
    message("[07d_simulator_figures] no validation metrics at ", metrics_path,
            " (a forward run is not validated): nothing to draw, skipped (--if-validated).")
    quit(save = "no", status = 0)
  }
  if (!file.exists(metrics_path))
    stop("metrics", sim_io_ext(metrics_fmt), " not found at: ", metrics_path,
         "\n  Run 07c_simulator_postprocess.R --validate first",
         if (identical(engine, "v1")) " (pass --engine spatpomp if this is a spatPomp run)." else ".")

  fig_dir <- if (!is.null(opt$fig_dir) && nzchar(opt$fig_dir)) {
    opt$fig_dir
  } else {
    file.path(run_dir, "figures")
  }

  message("[07d_simulator_figures] loading metrics from: ", metrics_path)
  metrics <- sim_io_read(metrics_path, metrics_fmt)

  make_sim_figures(metrics, fig_dir)

  message("[07d_simulator_figures] done. Figures written to: ", fig_dir)
}
