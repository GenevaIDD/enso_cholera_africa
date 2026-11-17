# 06_fit_figures_lag_compare.R — cross-lag comparison figures for single-lag
#   ENSO fits
#
# Usage:
#   Rscript analysis/06_fit_figures_lag_compare.R \
#     -d analysis/configs/enso_cb_singlelag_occ/
#   Rscript analysis/06_fit_figures_lag_compare.R -d analysis/configs/<batch>/ \
#     --format pdf
#
# Details:
#   Reads a whole single-lag sweep batch dir (track auto-detected per config).
#   Writes lag_compare_erf (ERFs of every lag overlaid), region and country
#   ERF by lag, and lag_compare_elpd (elpd by lag, best lag highlighted, from
#   the 05_fit_postprocess_loo_lags.R ranking located via utils/loo_lags.R).

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/figures_io.R")
source("analysis/utils/loo_lags.R") # latest_ranking_file()

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-d", "--batch_dir"),
    default = NULL,
    type = "character",
    help = "Directory of single-lag YAML configs (required)"
  ),
  make_option(
    c("--format"),
    default = "png",
    type = "character",
    help = "Output format: pdf | png [default: png]"
  ),
  make_option(
    c("--dpi"),
    default = 150L,
    type = "integer",
    help = "DPI for raster output [default: 150]"
  ),
  make_option(
    c("--vars"),
    default = NULL,
    type = "character",
    help = "Comma-separated CB variable names (default: all)"
  ),
  make_option(
    c("--group"),
    default = NULL,
    type = "character",
    help = paste(
      "Filter configs to this grouping group (e.g. enso_east).",
      "REQUIRED for grouped batches: without it every group's",
      "fits are overlaid on the same lag axis [default: all]"
    )
  )
)
script_opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(script_opt$batch_dir) || !dir.exists(script_opt$batch_dir)) {
  stop("Provide a valid batch directory via -d <dir/>.")
}

config_paths <- list.files(
  script_opt$batch_dir,
  pattern = "\\.yml$",
  full.names = TRUE
)
config_paths <- config_paths[basename(config_paths) != "manifest.yml"]
if (length(config_paths) == 0) {
  stop("No .yml configs found in: ", script_opt$batch_dir)
}

# Group filter. Same filename-token convention as 05_fit_postprocess_loo_lags.R:
# make_batch_configs.R strips non-alphanumerics from the group when naming
# configs, so `enso_east` appears as `_grp-ensoeast.yml`.
if (!is.null(script_opt$group)) {
  grp_fname <- gsub("[^A-Za-z0-9]", "", script_opt$group)
  keep <- grepl(paste0("_grp-", grp_fname, "\\.yml$"), basename(config_paths))
  config_paths <- config_paths[keep]
  if (length(config_paths) == 0) {
    stop(
      "No configs match group '",
      script_opt$group,
      "' (filename token: '",
      grp_fname,
      "') in: ",
      script_opt$batch_dir
    )
  }
  message(sprintf(
    "Group filter '%s': %d config(s) selected",
    script_opt$group,
    length(config_paths)
  ))
}

fmt <- match.arg(script_opt$format %||% "png", c("pdf", "png"))

# 3. Helpers -------------------------------------------------------------------

.fill_defaults <- function(opt, model) {
  defaults <- get_default_options()
  for (nm in setdiff(
    names(defaults),
    c("occurrence", "size", "duration", "monthly", "climate")
  )) {
    if (is.null(opt[[nm]])) opt[[nm]] <- defaults[[nm]]
  }
  model_defaults <- defaults[[model]] %||% list()
  for (nm in names(model_defaults)) {
    if (is.null(opt[[model]][[nm]])) opt[[model]][[nm]] <- model_defaults[[nm]]
  }
  opt
}

.resolve_single_lag <- function(model_blk) {
  dlnm_cfg <- model_blk$dlnm %||% list()
  sl <- dlnm_cfg$defaults$single_lag
  if (is.null(sl)) {
    for (vv in names(dlnm_cfg$variables %||% list())) {
      if (!is.null(dlnm_cfg$variables[[vv]]$single_lag)) {
        sl <- dlnm_cfg$variables[[vv]]$single_lag
        break
      }
    }
  }
  if (is.null(sl)) NA_integer_ else as.integer(sl)
}

# 4. Collect global ERF surfaces per lag ---------------------------------------

erf_rows <- list()
model <- NULL
for (cp in config_paths) {
  opt <- yaml::read_yaml(cp)
  mt <- {
    if (!is.null(opt$occurrence)) {
      "occurrence"
    } else if (!is.null(opt$size)) {
      "size"
    } else if (!is.null(opt$duration)) {
      "duration"
    } else if (!is.null(opt$monthly)) {
      "monthly"
    }
  }
  model <- model %||% mt
  opt <- .fill_defaults(opt, mt)
  sl <- .resolve_single_lag(opt[[mt]])
  if (is.na(sl)) {
    next
  }

  kpfx <- if (mt == "size") {
    "size_"
  } else if (mt == "duration") {
    "dur_"
  } else if (mt == "monthly") {
    "mon_"
  } else {
    ""
  }
  surf_file <- make_run_file(opt, kind = paste0(kpfx, "cb_surfaces"))
  if (!file.exists(surf_file)) {
    message("Surfaces missing, skipping lag ", sl, ": ", surf_file)
    next
  }
  surfaces <- readRDS(surf_file)
  vsel <- if (!is.null(script_opt$vars)) {
    intersect(trimws(strsplit(script_opt$vars, ",")[[1]]), names(surfaces))
  } else {
    names(surfaces)
  }
  for (v in vsel) {
    surf <- surfaces[[v]]
    if (is.null(surf) || is.null(surf$single_lag)) {
      next
    }
    erf_rows[[length(erf_rows) + 1L]] <- data.frame(
      variable = v,
      lag = sl,
      exposure = surf$at,
      logRR = surf$matfit[, 1],
      lo = surf$matlo[, 1],
      hi = surf$mathi[, 1],
      lo50 = surf$matlo50[, 1],
      hi50 = surf$mathi50[, 1],
      grouping = opt$grouping$source,
      grp = opt$grouping$group
    )
  }
}

# 5. Output directory ----------------------------------------------------------
# make_figure_dir() (utils/paths.R): passing the batch DIRECTORY as
# config_path reproduces this script's own basename(normalizePath(...))
# run_tag exactly -- derive_config_run_tag() is a no-op on an
# extension-less directory basename.
fig_dir <- make_figure_dir(
  opt = list(),
  track = "cb",
  config_path = script_opt$batch_dir
)

# Group suffix on every filename: without it a per-group loop overwrites one
# set of figures five times over. Appended to each call's stem below (not
# figure_saver()'s `prefix`, which lands BEFORE the stem) to keep the
# existing "<name><grp_sfx>.<fmt>" shape unchanged.
grp_sfx <- if (!is.null(script_opt$group)) paste0("_", script_opt$group) else ""

save_fig <- figure_saver(
  fig_dir,
  dpi = as.integer(script_opt$dpi %||% 150L),
  formats = fmt
)

batch <- derive_config_run_tag(script_opt$batch_dir) # batch dir basename

# 6. Figure 1: overlaid ERFs by lag --------------------------------------------

if (length(erf_rows) > 0) {
  erf_df <- dplyr::bind_rows(erf_rows) %>%
    dplyr::mutate(lag_f = factor(lag, levels = sort(unique(lag)))) %>%
    tibble::as_tibble()

  #   Save object
  erf_file <- here::here(
    "generated_data/",
    stringr::str_glue("erf_df_compare_lages_{batch}.rds")
  )
  saveRDS(erf_df, file = erf_file)

  p_erf <- ggplot(
    erf_df %>%
      dplyr::filter(variable == "enso"),
    aes(x = exposure, y = logRR, colour = lag_f, fill = lag_f, group = lag_f)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.4, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_grid(grp ~ lag_f) +
    scale_colour_viridis_d(name = "ENSO lag") +
    scale_fill_viridis_d(name = "ENSO lag") +
    labs(
      title = sprintf(
        "Single-lag exposure-response by lag (%s)",
        model %||% "?"
      ),
      subtitle = "Posterior median log-RR \u00b1 50% / 95% credible interval per lag",
      x = "Exposure",
      y = "log-RR"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  save_fig(p_erf, paste0("lag_compare_erf", grp_sfx), width = 12, height = 9)
} else {
  message(
    "No single-lag ERF surfaces found — skipping lag_compare_erf. ",
    "Run 05_fit_postprocess_cb.R on each config first."
  )
}

# 7. Figures 1b/1c: region and country ERF by lag ------------------------------
# Draw-based ribbons.
#
# For each lag config, loads the brms fit once to extract:
#   - population CB coefficient draws  [S × K]
#   - ADM1 RE draws                    [S × n_adm1 × K]
# Region ERF  = b_pop + r_region          (AFRO_region group, 4 regions)
# Country ERF = b_pop + r_region + r_ctry (AFRO_region:country group, 22 countries)
# ERF evaluated via Z_cen %*% t(b_draws) (centered at ENSO=0, matching extract_cb_surface cen=0).
# Uncertainty reflects posterior variability in CB coefs + mean RE per group.

if (length(erf_rows) > 0) {
  suppressPackageStartupMessages({
    library(brms)
    library(dlnm)
  })
  n_g_draws <- 200L
  # region_rows_g  <- list()
  country_rows_g <- list()

  for (cp in config_paths) {
    opt <- yaml::read_yaml(cp)
    mt <- {
      if (!is.null(opt$occurrence)) {
        "occurrence"
      } else if (!is.null(opt$size)) {
        "size"
      } else if (!is.null(opt$duration)) {
        "duration"
      } else if (!is.null(opt$monthly)) {
        "monthly"
      }
    }

    opt <- .fill_defaults(opt, mt)
    sl <- .resolve_single_lag(opt[[mt]])
    if (is.na(sl)) {
      next
    }

    kpfx <- if (mt == "size") {
      "size_"
    } else if (mt == "duration") {
      "dur_"
    } else if (mt == "monthly") {
      "mon_"
    } else {
      ""
    }
    fit_file <- make_run_file(opt, kind = paste0(kpfx, "fit"))
    meta_file <- make_run_file(opt, kind = paste0(kpfx, "cb_meta"))
    if (!file.exists(fit_file) || !file.exists(meta_file)) {
      message("Lag ", sl, ": fit/meta missing, skipping.")
      next
    }
    message(sprintf("Lag %2d: loading fit for region/country ERFs...", sl))
    brms_fit <- readRDS(fit_file)
    vm <- readRDS(meta_file)[["enso"]]

    # Country ERF = population CB draws + country RE draws, evaluated at 50
    # equally-spaced exposure points over vm$range, centered at exposure = 0.
    # extract_cb_country_surface() (utils/dlnm.R) is the shared implementation —
    # also used by 05_fit_postprocess_lag_compare_surface.R.
    country_surf <- extract_cb_country_surface(
      fit = brms_fit,
      meta = vm,
      at = NULL,
      countries = NULL,
      n_draws = n_g_draws,
      seed = 1234L,
      probs = c(0.025, 0.975),
      probs50 = c(0.25, 0.75)
    )

    for (ctry in names(country_surf)) {
      s <- country_surf[[ctry]]
      country_rows_g[[length(country_rows_g) + 1L]] <- data.frame(
        exposure = s$at,
        logRR = s$matfit,
        lo = s$matlo,
        hi = s$mathi,
        lo50 = s$matlo50,
        hi50 = s$mathi50,
        lag = sl,
        country = ctry,
        grouping = opt$grouping$source,
        group = opt$grouping$group
      )
    }

    rm(brms_fit)
    message(sprintf("Lag %2d: done.", sl))
  }

  # region_colours <- tryCatch(taxdat::colors_afro_regions(), error = function(e) NULL)

  # Figure 1b: region ERF by lag --------------------------------------------
  #
  # if (length(region_rows_g) > 0) {
  #   reg_df <- dplyr::bind_rows(region_rows_g) %>%
  #     dplyr::mutate(lag_f = factor(lag, levels = sort(unique(lag))))
  #
  #   p_erf_region <- ggplot(reg_df,
  #                          aes(x = exposure, y = logRR,
  #                              colour = AFRO_region, fill = AFRO_region,
  #                              group = AFRO_region)) +
  #     geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  #     geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  #     geom_line(linewidth = 0.7) +
  #     facet_grid(rows = vars(AFRO_region), cols = vars(lag_f)) +
  #     labs(
  #       title    = sprintf("ENSO exposure-response by region and lag (%s)", model %||% "?"),
  #       subtitle = "Pop. CB draws + region RE; ribbon = 95% CI; centered at ENSO = 0",
  #       x = "ENSO index", y = "log-RR",
  #       colour = "AFRO region", fill = "AFRO region"
  #     ) +
  #     theme_bw(base_size = 9) +
  #     theme(legend.position = "none",
  #           strip.text.x    = element_text(size = 7),
  #           strip.text.y    = element_text(size = 8, angle = 0, hjust = 0))
  #   if (!is.null(region_colours))
  #     p_erf_region <- p_erf_region +
  #     scale_colour_manual(values = region_colours) +
  #     scale_fill_manual(values  = region_colours)
  #   n_lags_reg <- dplyr::n_distinct(reg_df$lag)
  #   save_fig(p_erf_region, paste0("lag_compare_erf_region", grp_sfx),
  #            width = n_lags_reg * 1.6 + 2, height = 4 * 2.2 + 1.5)
  # } else {
  #   message("No region ERF rows — skipping lag_compare_erf_region.")
  # }

  # Figure 1c: country ERF by lag -------------------------------------------

  if (length(country_rows_g) > 0) {
    ctry_df <- dplyr::bind_rows(country_rows_g) %>%
      dplyr::mutate(lag_f = factor(lag, levels = sort(unique(lag)))) %>%
      tibble::as_tibble()

    #   Save object
    ctry_erf_file <- here::here(
      "generated_data/",
      stringr::str_glue("ctry_erf_df_compare_lages_{batch}.rds")
    )
    saveRDS(ctry_df, file = ctry_erf_file)

    n_ctry <- dplyr::n_distinct(ctry_df$country)

    p_erf_country <- ggplot(
      ctry_df,
      aes(
        x = exposure,
        y = logRR,
        colour = group,
        fill = group,
        group = country
      )
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08, colour = NA) +
      geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.16, colour = NA) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      facet_grid(rows = vars(group), cols = vars(lag_f), switch = "y") +
      labs(
        title = sprintf(
          "ENSO exposure-response by country and lag (%s)",
          model %||% "?"
        ),
        subtitle = sprintf(
          "Pop. CB draws + region + country RE (%d countries); ribbons = 50%% / 95%% CI; centered at ENSO = 0",
          n_ctry
        ),
        x = "ENSO index",
        y = "log-RR"
      ) +
      theme_bw(base_size = 8) +
      theme(
        legend.position = "none",
        strip.text.x = element_text(size = 6),
        strip.text.y.left = element_text(size = 6, angle = 0, hjust = 1),
        axis.text.x = element_text(size = 5),
        panel.spacing.y = unit(0.5, "pt")
      )
    # if (!is.null(region_colours))
    #   p_erf_country <- p_erf_country +
    #   scale_colour_manual(values = region_colours) +
    #   scale_fill_manual(values  = region_colours)
    n_lags_ctry <- dplyr::n_distinct(ctry_df$lag)
    save_fig(
      p_erf_country,
      paste0("lag_compare_erf_country", grp_sfx),
      width = n_lags_ctry * 1.6 + 2,
      height = n_ctry * 0.9 + 2,
      limitsize = FALSE
    )
  } else {
    message("No country ERF rows — skipping lag_compare_erf_country.")
  }
}

# 8. Figure 2: elpd by lag (from 05b ranking) ----------------------------------

# Under --group the ranking must be the SAME group's ranking, otherwise the
# elpd panel silently describes a different subregion than the ERF panels.
# utils/loo_lags.R: <outdir>/loo_compare/<batch>/ first, legacy in-config-dir
# fallback; latest by filename date, then mtime. group = NULL takes the
# latest ranking across any group.
rank_file <- latest_ranking_file(script_opt$batch_dir, group = script_opt$group)
if (!is.null(rank_file)) {
  message("Reading ranking: ", rank_file)
  ranking <- readRDS(rank_file)$ranking
  best_lag <- ranking$single_lag[which.max(ranking$elpd_loo)]
  plot_df <- ranking %>% dplyr::mutate(is_best = single_lag == best_lag)

  p_elpd <- ggplot(plot_df, aes(x = single_lag, y = elpd_loo)) +
    geom_errorbar(
      aes(ymin = elpd_loo - se_elpd, ymax = elpd_loo + se_elpd),
      width = 0.2,
      colour = "grey50"
    ) +
    geom_line(colour = "grey60") +
    geom_point(aes(colour = is_best), size = 3) +
    scale_colour_manual(
      values = c(`FALSE` = "black", `TRUE` = "firebrick"),
      guide = "none"
    ) +
    scale_x_continuous(breaks = sort(unique(plot_df$single_lag))) +
    labs(
      title = sprintf("ENSO single-lag elpd comparison (%s)", model %||% "?"),
      subtitle = sprintf("Best lag = %d", best_lag),
      x = "ENSO lag (months)",
      y = "elpd_loo (\u00b1 SE)"
    ) +
    theme_bw(base_size = 11)
  save_fig(p_elpd, paste0("lag_compare_elpd", grp_sfx), width = 7, height = 4.5)
} else {
  message(
    "No loo_compare_lags ranking found for ",
    script_opt$batch_dir,
    " — skipping lag_compare_elpd. Run 05_fit_postprocess_loo_lags.R first."
  )
}

message(sprintf("Done. Cross-lag figures written to %s/", fig_dir))
