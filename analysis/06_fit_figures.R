# 06_fit_figures.R — unified config-driven figures for occurrence and size tracks.
#
# Detects occurrence vs size from the YAML config block and routes to the
# appropriate plotting section. Data-overview figures (Figure 1) live in the
# standalone script 01_preprocess_figures_data_overview.R.
#
# Usage:
#   Rscript analysis/06_fit_figures.R -C analysis/configs/run1_benchmark/run1_brms_seasonal.yml
#   Rscript analysis/06_fit_figures.R -d analysis/configs/run1_benchmark/
#   Rscript analysis/06_fit_figures.R -C analysis/configs/size_run1_benchmark/size_run1_brms.yml
#   Rscript analysis/06_fit_figures.R -d analysis/configs/size_run1_benchmark/

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(optparse)
  library(sf)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/figures_io.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/provenance.R")
source("analysis/utils/plots.R")

# CLI options ----------------------------------------------------------------

option_list <- list(
  make_option(c("-C", "--config"),     default = NULL, type = "character",
              help = "Path to a single YAML config"),
  make_option(c("-d", "--config_dir"), default = NULL, type = "character",
              help = "Path to a directory of YAML configs (batch mode)"),
  make_option(c("--outdir"),           default = "figures/", type = "character",
              help = "Output directory for PNG files [default: figures/]")
)
cli_opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(xor(!is.null(cli_opt$config), !is.null(cli_opt$config_dir)))

config_paths <- if (!is.null(cli_opt$config)) {
  cli_opt$config
} else {
  list.files(cli_opt$config_dir, pattern = "\\.ya?ml$", full.names = TRUE)
}

# Figure directory creation is saveFigure()'s job (utils/figures_io.R) --
# one dir.create() per write, at the makeFigureDir()-resolved path.

# Static geodata (map figures) -----------------------------------------------
# MEMOISED ACCESSORS, not values -- built lazily so a track that draws no maps
# at all never pays for the geodata and never fails hard if it could not be
# built. Call as .afro_sf() / .adm1_sf(); each builds once and caches.
#
# giscoR >= 1.1 renamed the countries dataset gisco_countries ->
# gisco_countries_2024 (the old name errors with "data set not found"). The
# schema this code depends on is unchanged: ISO3_CODE joins to
# gisco_countrycode, which still carries `continent`.

.afro_sf <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    e <- environment()
    data(gisco_countries_2024, package = "giscoR", envir = e)
    data(gisco_countrycode,    package = "giscoR", envir = e)
    cache <<- get("gisco_countries_2024", envir = e) %>%
      inner_join(get("gisco_countrycode", envir = e) %>%
                   select(ISO3_CODE, continent) %>%
                   filter(continent == "Africa")) %>%
      nngeo::st_remove_holes() %>%
      rmapshaper::ms_filter_islands(min_area = 1e9) %>%
      rmapshaper::ms_simplify(keep = 0.5, keep_shapes = FALSE) %>%
      mutate(area = st_area(geometry)) %>%
      filter(as.numeric(area) > 1e10)
    cache
  }
})

.adm1_sf <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    cache <<- st_read("data/geodata/gadm_410-levels.gpkg",
                      layer = "ADM_1", quiet = TRUE)
    cache
  }
})

# Occurrence figures ---------------------------------------------------------

.figures_occurrence <- function(opt, tag, outdir) {
  re_stats_path <- makeRunFile(opt, kind = "re_stats")
  if (!file.exists(re_stats_path)) {
    warning("missing re_stats — run 05_fit_postprocess.R first: ", re_stats_path)
    return(invisible(NULL))
  }

  re_stats <- readRDS(re_stats_path)
  if (is.null(re_stats)) {
    message("re_stats is NULL — skipping model figures for: ", tag)
    return(invisible(NULL))
  }
  re_stats_region  <- re_stats$re_stats_region
  re_stats_country <- re_stats$re_stats_country
  re_stats_adm1    <- re_stats$re_stats_adm1

  pd <- position_dodge(.35)
  save <- figureSaver(dir = outdir, prefix = paste0(tag, "_"))

  # 1) Region-level effects --------------------------------------------------
  p_region <- re_stats_region %>%
    ggplot(aes(x = covar, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = q025, ymax = q975), alpha = .3, position = pd, width = 0) +
    ggh4x::facet_nested(. ~ var_type2 + covar_type, scales = "free_x", space = "free_x") +
    scale_color_manual(values = taxdat::colors_afro_regions()) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom", legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 9)) +
    labs(x = "covariate", y = "odds ratio")
  save(p_region, "region_effects", width = 9, height = 6)

  # 2) Country-level climatic effects ----------------------------------------
  p_country <- re_stats_country %>%
    filter(var_type2 == "climatic") %>%
    mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .35)) %>%
    ggplot(aes(x = AFRO_region, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
    geom_errorbar(aes(ymin = q025, ymax = q975, group = variable,
                      alpha = I(alpha)), position = pd, width = 0) +
    geom_point(aes(group = variable, alpha = I(alpha)), position = pd, size = 1.2) +
    ggrepel::geom_label_repel(aes(label = country, group = variable, alpha = I(alpha)),
                              position = pd,
                              min.segment.length = 0, max.iter = 1e4, max.overlaps = Inf) +
    ggh4x::facet_nested_wrap(. ~ covar_type + covar, scales = "free_x") +
    scale_color_manual(values = taxdat::colors_afro_regions()) +
    labs(x = "AFRO region", y = "hazard ratio") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 4))
  save(p_country, "country_effects_climate", width = 10, height = 8)

  # 3) Country ongoing-outbreak effects --------------------------------------
  p_ongoing <- re_stats_country %>%
    filter(var_type1 == "ongoing outbreak lag1") %>%
    mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .35)) %>%
    ggplot(aes(x = country, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
    geom_errorbar(aes(ymin = q025, ymax = q975, group = variable,
                      alpha = I(alpha)), position = pd, width = 0) +
    geom_point(aes(group = variable, alpha = I(alpha)), position = pd, size = 1.2) +
    facet_grid(~ AFRO_region, scales = "free_x", switch = "y", space = "free") +
    scale_color_manual(values = taxdat::colors_afro_regions()) +
    labs(x = "AFRO region", y = "hazard ratio") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.text = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(0, 10))
  save(p_ongoing, "country_effects_ongoing", width = 8, height = 4)

  # 4) Country time-to-last-outbreak effects ---------------------------------
  p_country_time_outbreak <- re_stats_country %>%
    filter(var_type1 == "time to last outbreak") %>%
    mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .35)) %>%
    ggplot(aes(x = var, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
    geom_errorbar(aes(ymin = q025, ymax = q975, group = variable,
                      alpha = I(alpha)), position = pd, width = 0) +
    geom_point(aes(group = variable, alpha = I(alpha)), position = pd, size = 1.2) +
    ggh4x::facet_nested(AFRO_region ~ ., scales = "free_x", switch = "y", space = "free") +
    scale_color_manual(values = taxdat::colors_afro_regions()) +
    labs(x = "AFRO region", y = "hazard ratio") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.text = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(0, 9))
  save(p_country_time_outbreak, "country_effects_time_to_last_outbreak",
       width = 6, height = 10)

  # 5) ADM1-level climatic effects -------------------------------------------
  p_adm1 <- re_stats_adm1 %>%
    filter(var_type2 == "climatic") %>%
    mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .2)) %>%
    ggplot(aes(x = country, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
    geom_errorbar(aes(ymin = q025, ymax = q975, group = adm1_unit,
                      alpha = I(alpha)), position = pd, width = 0) +
    geom_point(aes(group = adm1_unit, alpha = I(alpha)), position = pd, size = 1.2) +
    ggh4x::facet_nested_wrap(. ~ covar_type + covar, scales = "free_x") +
    scale_color_manual(values = taxdat::colors_afro_regions()) +
    labs(x = "AFRO region", y = "hazard ratio") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(0, 4))
  save(p_adm1, "adm1_effects_climate", width = 10, height = 8)

  # 6) Country-level effects map ---------------------------------------------
  res_sf <- .afro_sf() %>%
    inner_join(re_stats_country, by = c("ISO3_CODE" = "country")) %>%
    mutate(covar = var_type0) %>%
    filter(var_type2 == "climatic")

  p_map <- .afro_sf() %>%
    ggplot() +
    geom_sf(fill = "darkgray", color = "gray") +
    geom_sf(data = res_sf, aes(fill = mean)) +
    facet_grid(covar_type ~ covar) +
    scale_fill_gradient2(midpoint = 1,
                         low = scales::muted("blue"), high = scales::muted("red")) +
    taxdat::map_theme() +
    theme(panel.border = element_blank()) +
    labs(fill = "odds ratio")
  save(p_map, "country_effects_climate_map", width = 8, height = 8)

  # 7) ADM1-level effects map ------------------------------------------------
  res_adm1_sf <- .adm1_sf() %>%
    mutate(GID_1 = gidAdm1(GID_1, strict = FALSE)) %>%
    inner_join(re_stats_adm1 %>%
                 mutate(mean = if_else(mean > 2, 2, mean)),
               by = c("GID_1" = "adm1_unit")) %>%
    mutate(covar = var_type0)

  p_map_adm1 <- .afro_sf() %>%
    ggplot() +
    geom_sf(fill = "darkgray", color = "gray") +
    geom_sf(data = res_adm1_sf %>% filter(var_type2 == "climatic"),
            aes(fill = mean), lwd = .1) +
    geom_sf(color = "black", alpha = 0) +
    facet_grid(covar_type ~ covar) +
    scale_fill_gradient2(midpoint = 1,
                         low = scales::muted("blue"), high = scales::muted("red")) +
    taxdat::map_theme() +
    theme(panel.border = element_blank()) +
    labs(fill = "odds ratio")
  save(p_map_adm1, "adm1_effects_climate_map", width = 8, height = 8)

  # 8) ADM1 intercept map ----------------------------------------------------
  p_map_adm1_I <- .afro_sf() %>%
    ggplot() +
    geom_sf(fill = "darkgray", color = "gray") +
    geom_sf(data = res_adm1_sf %>% filter(var == "(Intercept)"),
            aes(fill = mean), lwd = .1) +
    geom_sf(color = "black", alpha = 0) +
    taxdat::map_theme() +
    theme(panel.border = element_blank()) +
    labs(fill = "odds ratio")
  save(p_map_adm1_I, "adm1_effects_intercept_map", width = 8, height = 8)

  # 9) Mean vs anomaly quadrant scatter --------------------------------------
  S_sqrt     <- function(x) sign(x) * sqrt(abs(x))
  IS_sqrt    <- function(x) x^2 * sign(x)
  S_sqrt_trans <- function() scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

  adm1_climatic <- re_stats_adm1 %>%
    filter(var_type2 == "climatic") %>%
    select(-variable) %>% ungroup()
  has_climatology_anomaly <- all(c("climatology", "anomaly") %in% adm1_climatic$var_type1)
  if (has_climatology_anomaly) {
    p_quadrants <- adm1_climatic %>%
      select(adm1_unit, AFRO_region, country, mean, q025, q975, var_type0, var_type1) %>%
      pivot_wider(
        id_cols     = c(adm1_unit, AFRO_region, country, var_type0),
        names_from  = var_type1,
        values_from = c(mean, q025, q975),
        names_sep   = "."
      ) %>%
      dplyr::rename(
        mean     = `mean.climatology`,
        q025     = `q025.climatology`,
        q975     = `q975.climatology`,
        mean_dev = `mean.anomaly`,
        q025_dev = `q025.anomaly`,
        q975_dev = `q975.anomaly`
      ) %>%
      mutate(covar = var_type0) %>%
      ggplot(aes(x = mean, y = mean_dev, color = AFRO_region)) +
      geom_point(alpha = .6) +
      geom_errorbar(aes(ymin = q025_dev, ymax = q975_dev), width = 0, alpha = .05) +
      geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, alpha = .05) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = .4) +
      geom_vline(aes(xintercept = 1), lty = 2, lwd = .4) +
      facet_grid(. ~ covar) +
      coord_cartesian(xlim = c(0, 3), ylim = c(0, 3)) +
      scale_color_manual(values = taxdat::colors_afro_regions()) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(x = "odds ratio of climatology", y = "odds ratio of anomaly")
    save(p_quadrants, "adm1_effects_quadrants", width = 8, height = 5)
  } else {
    message("Skipping quadrant plot: model lacks both climatology and anomaly ADM1 terms.")
  }

  # 10) Climate data tile (diagnostic) ---------------------------------------
  plot_clim_tile <- function(surv_dat_all, variable) {
    if (str_detect(variable, "temp")) {
      this_limits <- c(-5, 5)
      this_trans  <- "identity"
    } else {
      this_limits <- c(-400, 400)
      this_trans  <- S_sqrt_trans()
    }
    surv_dat_all %>%
      filter(admin_level != "ADM0") %>%
      group_by(t_start, ADM1) %>%
      slice(1) %>%
      mutate(AFRO_region = factor(AFRO_region,
                                  levels = c("Western Africa", "Central Africa",
                                             "Eastern Africa", "Southern Africa"))) %>%
      ggplot(aes(x = t_start, y = ADM1)) +
      geom_tile(aes(fill = !!rlang::sym(variable)), color = "darkgray", lwd = .01) +
      ggh4x::facet_nested(AFRO_region + country ~ ., scales = "free", space = "free_y",
                          switch = "y") +
      theme_bw() +
      scale_fill_gradient2(limits = this_limits, transform = this_trans,
                           low = scales::muted("blue"), high = scales::muted("red")) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y")
  }

  surv_dat_all <- tryCatch(
    readRDS(makeSurvDatFile(opt = list(baseline = opt$baseline))),
    error = function(e) { warning("could not load surv_dat: ", conditionMessage(e)); NULL }
  )

  if (!is.null(surv_dat_all)) {
    plot_list <- map(str_c("dev_", names(getVarTypeDic())),
                     function(x) plot_clim_tile(surv_dat_all = surv_dat_all, variable = x))
    p_clim_dat <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
    save(p_clim_dat, "clim_data_tile", width = 25, height = 15)
  }

  message("wrote occurrence figures with prefix: ", file.path(outdir, tag))
  invisible(NULL)
}

# Duration figures ------------------------------------------------------------
# Reuses the occurrence re_stats structure (computeNestedSurvReStatsBrms() is
# called unforked for duration in 05_fit_postprocess.R), so the
# same var_type1/var_type2/covar/covar_type decoration is already present.
# Only the climatic-effect blocks apply: duration has no ongoing_outbreak_*
# or time_to_last_outbreak covariates (not part of the cessation-hazard
# formula, see 04_fit_dur_brms.R), and no makeSurvDatFile()-based climate
# data tile (duration uses the active-month panel, a different data file).

.figures_duration <- function(opt, tag, outdir) {
  re_stats_path <- makeRunFile(opt, kind = "dur_re_stats")
  if (!file.exists(re_stats_path)) {
    warning("missing re_stats — run 05_fit_postprocess.R first: ", re_stats_path)
    return(invisible(NULL))
  }

  re_stats <- readRDS(re_stats_path)
  if (is.null(re_stats)) {
    message("re_stats is NULL — skipping model figures for: ", tag)
    return(invisible(NULL))
  }
  re_stats_region  <- re_stats$re_stats_region
  re_stats_country <- re_stats$re_stats_country
  re_stats_adm1    <- re_stats$re_stats_adm1

  pd <- position_dodge(.35)
  save <- figureSaver(dir = outdir, prefix = paste0(tag, "_"))

  # Existence guards for the climatic-effect blocks below. These fits can be
  # ENSO-only (duration.covariates: enso only, no precip/temp fixed effects,
  # as in prod_run_dur_grouped_enso_singlelag), in which case var_type2 is
  # "non-climatic" for every RE term and the "climatic"-filtered tibbles below
  # are legitimately 0 rows; under grouping mode (opt$grouping$enable) REStats
  # is also missing AFRO_region entirely, so re_stats_region is NULL
  # (computeNestedSurvReStatsBrms drops the region level when it isn't in the
  # formula). ggplot2/ggh4x error on a faceting variable with zero values
  # rather than drawing an empty panel, so each block needs its own guard --
  # same pattern as the quadrant-plot guard below (has_climatology_anomaly).
  has_region           <- !is.null(re_stats_region) && nrow(re_stats_region) > 0
  has_country_climatic <- !is.null(re_stats_country) &&
    any(re_stats_country$var_type2 == "climatic")
  has_adm1_climatic    <- !is.null(re_stats_adm1) &&
    any(re_stats_adm1$var_type2 == "climatic")

  # 1) Region-level effects --------------------------------------------------
  if (has_region) {
    p_region <- re_stats_region %>%
      ggplot(aes(x = covar, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
      geom_point(position = pd) +
      geom_errorbar(aes(ymin = q025, ymax = q975), alpha = .3, position = pd, width = 0) +
      ggh4x::facet_nested(. ~ var_type2 + covar_type, scales = "free_x", space = "free_x") +
      scale_color_manual(values = taxdat::colors_afro_regions()) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "bottom", legend.title = element_blank()) +
      coord_cartesian(ylim = c(0, 9)) +
      labs(x = "covariate", y = "hazard ratio")
    save(p_region, "region_effects", width = 9, height = 6)
  } else {
    message("Skipping region-level effects plot: re_stats_region is NULL/empty ",
            "(grouping mode drops AFRO_region from the RE formula).")
  }

  # 2) Country-level climatic effects ----------------------------------------
  if (has_country_climatic) {
    p_country <- re_stats_country %>%
      filter(var_type2 == "climatic") %>%
      mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .35)) %>%
      ggplot(aes(x = AFRO_region, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = variable,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = variable, alpha = I(alpha)), position = pd, size = 1.2) +
      ggrepel::geom_label_repel(aes(label = country, group = variable, alpha = I(alpha)),
                                position = pd,
                                min.segment.length = 0, max.iter = 1e4, max.overlaps = Inf) +
      ggh4x::facet_nested_wrap(. ~ covar_type + covar, scales = "free_x") +
      scale_color_manual(values = taxdat::colors_afro_regions()) +
      labs(x = "AFRO region", y = "hazard ratio") +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "bottom", legend.title = element_blank()) +
      coord_cartesian(ylim = c(0, 4))
    save(p_country, "country_effects_climate", width = 10, height = 8)
  } else {
    message("Skipping country-level climatic effects plot: model has no ",
            "climatic (var_type2 == \"climatic\") country RE terms.")
  }

  # 3) ADM1-level climatic effects -------------------------------------------
  if (has_adm1_climatic) {
    p_adm1 <- re_stats_adm1 %>%
      filter(var_type2 == "climatic") %>%
      mutate(alpha = if_else(q025 > 1 | q975 < 1, .8, .2)) %>%
      ggplot(aes(x = country, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = .35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = adm1_unit,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = adm1_unit, alpha = I(alpha)), position = pd, size = 1.2) +
      ggh4x::facet_nested_wrap(. ~ covar_type + covar, scales = "free_x") +
      scale_color_manual(values = taxdat::colors_afro_regions()) +
      labs(x = "AFRO region", y = "hazard ratio") +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      coord_cartesian(ylim = c(0, 4))
    save(p_adm1, "adm1_effects_climate", width = 10, height = 8)
  } else {
    message("Skipping ADM1-level climatic effects plot: model has no ",
            "climatic (var_type2 == \"climatic\") ADM1 RE terms.")
  }

  # 4) Country-level effects map ---------------------------------------------
  if (has_country_climatic) {
    res_sf <- .afro_sf() %>%
      inner_join(re_stats_country, by = c("ISO3_CODE" = "country")) %>%
      mutate(covar = var_type0) %>%
      filter(var_type2 == "climatic")

    p_map <- .afro_sf() %>%
      ggplot() +
      geom_sf(fill = "darkgray", color = "gray") +
      geom_sf(data = res_sf, aes(fill = mean)) +
      facet_grid(covar_type ~ covar) +
      scale_fill_gradient2(midpoint = 1,
                           low = scales::muted("blue"), high = scales::muted("red")) +
      taxdat::map_theme() +
      theme(panel.border = element_blank()) +
      labs(fill = "hazard ratio")
    save(p_map, "country_effects_climate_map", width = 8, height = 8)
  } else {
    message("Skipping country-level climatic effects map: model has no ",
            "climatic (var_type2 == \"climatic\") country RE terms.")
  }

  # 5) ADM1-level effects map -------------------------------------------------
  # res_adm1_sf is built unconditionally (block 6's intercept map needs it
  # regardless of whether any climatic terms exist); only the climatic-filtered
  # ggplot below is guarded.
  res_adm1_sf <- .adm1_sf() %>%
    mutate(GID_1 = gidAdm1(GID_1, strict = FALSE)) %>%
    inner_join(re_stats_adm1 %>%
                 mutate(mean = if_else(mean > 2, 2, mean)),
               by = c("GID_1" = "adm1_unit")) %>%
    mutate(covar = var_type0)

  if (has_adm1_climatic) {
    p_map_adm1 <- .afro_sf() %>%
      ggplot() +
      geom_sf(fill = "darkgray", color = "gray") +
      geom_sf(data = res_adm1_sf %>% filter(var_type2 == "climatic"),
              aes(fill = mean), lwd = .1) +
      geom_sf(color = "black", alpha = 0) +
      facet_grid(covar_type ~ covar) +
      scale_fill_gradient2(midpoint = 1,
                           low = scales::muted("blue"), high = scales::muted("red")) +
      taxdat::map_theme() +
      theme(panel.border = element_blank()) +
      labs(fill = "hazard ratio")
    save(p_map_adm1, "adm1_effects_climate_map", width = 8, height = 8)
  } else {
    message("Skipping ADM1-level climatic effects map: model has no ",
            "climatic (var_type2 == \"climatic\") ADM1 RE terms.")
  }

  # 6) ADM1 intercept map ----------------------------------------------------
  # var == "(Intercept)" is the stan_surv/rstanarm naming convention (its `b[]`
  # variable strings carry the parens literally, per extractSurvReDraws() in
  # occ_post.R). computeNestedSurvReStatsBrms() -- the adapter duration always
  # uses (04_fit_dur_brms.R / 04_fit_dur_brms_cb.R are both brms fits, never
  # stan_surv) -- derives var from the stripped brms parameter name `term`
  # (b_Intercept -> "Intercept", no parens), so a bare "(Intercept)" match
  # here silently filtered out every row and produced an all-gray map.
  # Matching both spellings keeps this block correct for either backend.
  p_map_adm1_I <- .afro_sf() %>%
    ggplot() +
    geom_sf(fill = "darkgray", color = "gray") +
    geom_sf(data = res_adm1_sf %>% filter(var %in% c("(Intercept)", "Intercept")),
            aes(fill = mean), lwd = .1) +
    geom_sf(color = "black", alpha = 0) +
    taxdat::map_theme() +
    theme(panel.border = element_blank()) +
    labs(fill = "hazard ratio")
  save(p_map_adm1_I, "adm1_effects_intercept_map", width = 8, height = 8)

  # 7) Mean vs anomaly quadrant scatter --------------------------------------
  S_sqrt     <- function(x) sign(x) * sqrt(abs(x))
  IS_sqrt    <- function(x) x^2 * sign(x)
  S_sqrt_trans <- function() scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

  adm1_climatic <- re_stats_adm1 %>%
    filter(var_type2 == "climatic") %>%
    select(-variable) %>% ungroup()
  has_climatology_anomaly <- all(c("climatology", "anomaly") %in% adm1_climatic$var_type1)
  if (has_climatology_anomaly) {
    p_quadrants <- adm1_climatic %>%
      select(adm1_unit, AFRO_region, country, mean, q025, q975, var_type0, var_type1) %>%
      pivot_wider(
        id_cols     = c(adm1_unit, AFRO_region, country, var_type0),
        names_from  = var_type1,
        values_from = c(mean, q025, q975),
        names_sep   = "."
      ) %>%
      dplyr::rename(
        mean     = `mean.climatology`,
        q025     = `q025.climatology`,
        q975     = `q975.climatology`,
        mean_dev = `mean.anomaly`,
        q025_dev = `q025.anomaly`,
        q975_dev = `q975.anomaly`
      ) %>%
      mutate(covar = var_type0) %>%
      ggplot(aes(x = mean, y = mean_dev, color = AFRO_region)) +
      geom_point(alpha = .6) +
      geom_errorbar(aes(ymin = q025_dev, ymax = q975_dev), width = 0, alpha = .05) +
      geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, alpha = .05) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = .4) +
      geom_vline(aes(xintercept = 1), lty = 2, lwd = .4) +
      facet_grid(. ~ covar) +
      coord_cartesian(xlim = c(0, 3), ylim = c(0, 3)) +
      scale_color_manual(values = taxdat::colors_afro_regions()) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(x = "hazard ratio of climatology", y = "hazard ratio of anomaly")
    save(p_quadrants, "adm1_effects_quadrants", width = 8, height = 5)
  } else {
    message("Skipping quadrant plot: model lacks both climatology and anomaly ADM1 terms.")
  }

  # Posterior-predictive check (duration-specific) ---------------------------
  # Two stacked panels, deliberately NOT the size track's predicted-vs-observed
  # scatter: the duration response is a per-active-month 0/1 cessation
  # indicator, so what is checkable is (a) the monthly hazard against
  # months-since-onset and (b) the IMPLIED duration distribution, which is the
  # actual endpoint. Built by makeDurationPpCheckTbl() (utils/dur_post.R).
  pp_check_path <- makeRunFile(opt, kind = "dur_pp_check")
  if (file.exists(pp_check_path)) {
    pp <- readRDS(pp_check_path)

    p_haz <- pp$hazard %>%
      ggplot(aes(x = mso)) +
      geom_ribbon(aes(ymin = rep_lo, ymax = rep_hi),
                  fill = "steelblue", alpha = .3) +
      geom_line(aes(y = rep_med), colour = "steelblue4") +
      geom_point(aes(y = obs_hazard), size = 1.6, colour = "black") +
      scale_x_continuous(breaks = scales::breaks_width(2)) +
      theme_bw() +
      labs(subtitle = "Monthly cessation hazard (points = observed, band = 95% posterior predictive)",
           x = "months since onset", y = "cessation hazard")

    p_dur <- pp$duration %>%
      ggplot(aes(x = duration_months)) +
      geom_ribbon(aes(ymin = rep_lo, ymax = rep_hi),
                  fill = "darkorange", alpha = .3) +
      geom_line(aes(y = rep_med), colour = "darkorange3") +
      geom_point(aes(y = obs_n), size = 1.6, colour = "black") +
      scale_x_continuous(breaks = scales::breaks_width(2)) +
      theme_bw() +
      labs(subtitle = paste("Implied duration distribution within the observed at-risk window",
                            "\n(points = observed; compare SHAPE, not totals — see dur_post.R header)"),
           x = "outbreak duration (months)", y = "number of episodes")

    p_pp <- cowplot::plot_grid(p_haz, p_dur, ncol = 1, align = "v")
    save(p_pp, "pp_check", width = 7, height = 8)

    message(sprintf(paste0("  pp_check: %d episodes, %d observed cessations, ",
                           "%d censored, %d draws"),
                    pp$meta$n_episodes, pp$meta$n_obs_cessations,
                    pp$meta$n_censored, pp$meta$ndraws))
  }

  message("wrote duration figures with prefix: ", file.path(outdir, tag))
  invisible(NULL)
}

# Size figures ---------------------------------------------------------------

.figures_size <- function(opt, tag, outdir) {
  re_stats_path <- makeSizeRunFile(opt, kind = "size_re_stats")
  pp_check_path <- makeSizeRunFile(opt, kind = "size_pp_check")
  if (!file.exists(re_stats_path)) {
    warning("missing re_stats — run 05_fit_postprocess.R first: ", re_stats_path)
    return(invisible(NULL))
  }

  re_stats <- readRDS(re_stats_path)
  re_region  <- decorate_covar(re_stats$region)
  re_country <- decorate_covar(re_stats$country)
  re_adm1    <- decorate_covar(re_stats$adm1)

  pd <- position_dodge(0.35)
  save <- figureSaver(dir = outdir, prefix = paste0(tag, "_"))

  # 1) Region-level effects --------------------------------------------------
  p_region <- re_region %>%
    dplyr::filter(var != "Intercept") %>%
    ggplot(aes(x = covar, y = mean, color = AFRO_region)) +
    geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = q025, ymax = q975), alpha = 0.3,
                  position = pd, width = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom", legend.title = element_blank()) +
    coord_cartesian(ylim = c(0, 9)) +
    labs(x = "covariate", y = "rate ratio")
  save(p_region, "region_effects", width = 9, height = 6)

  # 2) Country-level climatic effects ----------------------------------------
  if (nrow(re_country) > 0) {
    p_country <- re_country %>%
      dplyr::filter(var_type2 == "climatic") %>%
      dplyr::mutate(alpha = dplyr::if_else(q025 > 1 | q975 < 1, 0.8, 0.35)) %>%
      ggplot(aes(x = AFRO_region, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = country,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = country, alpha = I(alpha)),
                 position = pd, size = 1.2) +
      labs(x = "AFRO region", y = "rate ratio") +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "bottom", legend.title = element_blank()) +
      coord_cartesian(ylim = c(0, 4))
    save(p_country, "country_effects", width = 8, height = 8)
  }

  # 3) ADM1-level climatic effects -------------------------------------------
  if (nrow(re_adm1) > 0) {
    p_adm1 <- re_adm1 %>%
      dplyr::filter(var_type2 == "climatic") %>%
      dplyr::mutate(alpha = dplyr::if_else(q025 > 1 | q975 < 1, 0.8, 0.2)) %>%
      ggplot(aes(x = country, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = adm1_unit,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = adm1_unit, alpha = I(alpha)),
                 position = pd, size = 1.2) +
      labs(x = "country", y = "rate ratio") +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      coord_cartesian(ylim = c(0, 4))
    save(p_adm1, "adm1_effects", width = 10, height = 8)
  }

  # 3b) ADM1-level climatic effects map ----------------------------------------
  if (nrow(re_adm1) > 0) {
    re_adm1_climatic <- re_adm1 %>% dplyr::filter(var_type2 == "climatic")
    if (nrow(re_adm1_climatic) > 0) {
      res_adm1_sf <- .adm1_sf() %>%
        dplyr::mutate(GID_1 = gidAdm1(GID_1, strict = FALSE)) %>%
        dplyr::inner_join(re_adm1_climatic %>%
                            dplyr::mutate(mean = dplyr::if_else(mean > 2, 2, mean)),
                          by = c("GID_1" = "adm1_unit"))
      p_map_adm1 <- .afro_sf() %>%
        ggplot() +
        geom_sf(fill = "darkgray", color = "gray") +
        geom_sf(data = res_adm1_sf, aes(fill = mean), lwd = .1) +
        geom_sf(color = "black", alpha = 0) +
        facet_wrap(~ covar) +
        scale_fill_gradient2(midpoint = 1,
                         low = scales::muted("blue"), high = scales::muted("red")) +
        taxdat::map_theme() +
        theme(panel.border = element_blank()) +
        labs(fill = "rate ratio")
      save(p_map_adm1, "adm1_effects_climate_map", width = 8, height = 8)
    }
  }

  # 4) Posterior-predictive scatter ------------------------------------------
  if (file.exists(pp_check_path)) {
    pp <- readRDS(pp_check_path)
    p_pp <- pp %>%
      ggplot(aes(x = obs, y = Estimate)) +
      geom_point(alpha = 0.4) +
      geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
      geom_abline(aes(intercept = 0, slope = 1)) +
      scale_x_log10() + scale_y_log10() +
      theme_bw() +
      labs(x = "observed suspected cases", y = "predicted (posterior mean)")
    save(p_pp, "pp_check", width = 6, height = 6)
  }

  message("wrote size figures with prefix: ", file.path(outdir, tag))
  invisible(NULL)
}

# Monthly-size figures --------------------------------------------------------
# Reuses computeNestedSizeReStatsBrms() (05_fit_postprocess.R's "monthly"
# re_stats branch) and makeSizePpCheckTbl(fit, response = "cases_month")
# unchanged -- same schema as the size track's re_stats/pp_check, so the
# plotting code below is a straight copy of .figures_size() with
# makeSizeRunFile() -> makeRunFile() (the "mon_" kind prefix routes through
# makeMonthlyRunFile(), not makeSizeRunFile() -- see 05_fit_postprocess.R's
# kpfx) and the pp_check axis label updated for the episode-month response.
# Region-level effects are guarded (unlike .figures_size()'s unconditional
# block 1): the monthly track's primary usage is ENSO grouping mode
# (mon_smoke.yml), under which computeNestedSizeReStatsBrms(drop_region=TRUE)
# returns re_stats$region == NULL, mirroring duration's has_region guard.

.figures_monthly <- function(opt, tag, outdir) {
  re_stats_path <- makeRunFile(opt, kind = "mon_re_stats")
  pp_check_path <- makeRunFile(opt, kind = "mon_pp_check")
  if (!file.exists(re_stats_path)) {
    warning("missing re_stats — run 05_fit_postprocess.R first: ", re_stats_path)
    return(invisible(NULL))
  }

  re_stats <- readRDS(re_stats_path)
  # decorate_covar() calls dplyr::mutate() unconditionally and errors on NULL
  # -- guard here (re_stats$region is NULL whenever grouping mode drops
  # AFRO_region, the monthly track's primary usage per mon_smoke.yml).
  re_region  <- if (!is.null(re_stats$region)) decorate_covar(re_stats$region) else NULL
  re_country <- decorate_covar(re_stats$country)
  re_adm1    <- decorate_covar(re_stats$adm1)

  pd <- position_dodge(0.35)
  save <- figureSaver(dir = outdir, prefix = paste0(tag, "_"))

  # 1) Region-level effects --------------------------------------------------
  if (!is.null(re_region) && nrow(re_region) > 0) {
    p_region <- re_region %>%
      dplyr::filter(var != "Intercept") %>%
      ggplot(aes(x = covar, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
      geom_point(position = pd) +
      geom_errorbar(aes(ymin = q025, ymax = q975), alpha = 0.3,
                    position = pd, width = 0) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "bottom", legend.title = element_blank()) +
      coord_cartesian(ylim = c(0, 9)) +
      labs(x = "covariate", y = "rate ratio")
    save(p_region, "region_effects", width = 9, height = 6)
  } else {
    message("Skipping region-level effects plot: re_stats$region is NULL/empty ",
            "(grouping mode drops AFRO_region from the RE formula).")
  }

  # 2) Country-level climatic effects ----------------------------------------
  if (nrow(re_country) > 0) {
    p_country <- re_country %>%
      dplyr::filter(var_type2 == "climatic") %>%
      dplyr::mutate(alpha = dplyr::if_else(q025 > 1 | q975 < 1, 0.8, 0.35)) %>%
      ggplot(aes(x = AFRO_region, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = country,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = country, alpha = I(alpha)),
                 position = pd, size = 1.2) +
      labs(x = "AFRO region", y = "rate ratio") +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "bottom", legend.title = element_blank()) +
      coord_cartesian(ylim = c(0, 4))
    save(p_country, "country_effects", width = 8, height = 8)
  }

  # 3) ADM1-level climatic effects -------------------------------------------
  if (nrow(re_adm1) > 0) {
    p_adm1 <- re_adm1 %>%
      dplyr::filter(var_type2 == "climatic") %>%
      dplyr::mutate(alpha = dplyr::if_else(q025 > 1 | q975 < 1, 0.8, 0.2)) %>%
      ggplot(aes(x = country, y = mean, color = AFRO_region)) +
      geom_hline(aes(yintercept = 1), lty = 2, lwd = 0.35) +
      geom_errorbar(aes(ymin = q025, ymax = q975, group = adm1_unit,
                        alpha = I(alpha)), position = pd, width = 0) +
      geom_point(aes(group = adm1_unit, alpha = I(alpha)),
                 position = pd, size = 1.2) +
      labs(x = "country", y = "rate ratio") +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      coord_cartesian(ylim = c(0, 4))
    save(p_adm1, "adm1_effects", width = 10, height = 8)
  }

  # 3b) ADM1-level climatic effects map ----------------------------------------
  if (nrow(re_adm1) > 0) {
    re_adm1_climatic <- re_adm1 %>% dplyr::filter(var_type2 == "climatic")
    if (nrow(re_adm1_climatic) > 0) {
      res_adm1_sf <- .adm1_sf() %>%
        dplyr::mutate(GID_1 = gidAdm1(GID_1, strict = FALSE)) %>%
        dplyr::inner_join(re_adm1_climatic %>%
                            dplyr::mutate(mean = dplyr::if_else(mean > 2, 2, mean)),
                          by = c("GID_1" = "adm1_unit"))
      p_map_adm1 <- .afro_sf() %>%
        ggplot() +
        geom_sf(fill = "darkgray", color = "gray") +
        geom_sf(data = res_adm1_sf, aes(fill = mean), lwd = .1) +
        geom_sf(color = "black", alpha = 0) +
        facet_wrap(~ covar) +
        scale_fill_gradient2(midpoint = 1,
                         low = scales::muted("blue"), high = scales::muted("red")) +
        taxdat::map_theme() +
        theme(panel.border = element_blank()) +
        labs(fill = "rate ratio")
      save(p_map_adm1, "adm1_effects_climate_map", width = 8, height = 8)
    }
  }

  # 4) Posterior-predictive scatter (cases_month) ----------------------------
  if (file.exists(pp_check_path)) {
    pp <- readRDS(pp_check_path)
    p_pp <- pp %>%
      ggplot(aes(x = obs, y = Estimate)) +
      geom_point(alpha = 0.4) +
      geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
      geom_abline(aes(intercept = 0, slope = 1)) +
      scale_x_log10() + scale_y_log10() +
      theme_bw() +
      labs(x = "observed monthly cases", y = "predicted (posterior mean)")
    save(p_pp, "pp_check", width = 6, height = 6)
  }

  message("wrote monthly figures with prefix: ", file.path(outdir, tag))
  invisible(NULL)
}

# Dispatch -------------------------------------------------------------------

for (cfg in config_paths) {
  message(">>> figures: ", cfg)
  raw_cfg <- yaml::read_yaml(cfg)
  if (!is.null(raw_cfg$occurrence)) {
    model <- "occurrence"
  } else if (!is.null(raw_cfg$size)) {
    model <- "size"
  } else if (!is.null(raw_cfg$duration)) {
    model <- "duration"
  } else if (!is.null(raw_cfg$monthly)) {
    model <- "monthly"
  } else {
    warning("config has none of 'occurrence:', 'size:', 'duration:', 'monthly:' block — skipping: ", cfg)
    next
  }

  opt <- make_options_from_config(
    option_list = list(),
    cli_argv    = c("--config", cfg),
    model       = model
  )
  tag <- tools::file_path_sans_ext(basename(cfg))
  # makeFigureDir() (utils/paths.R): track = model (occurrence/size/
  # duration) groups output by track; config_path = dirname(cfg), NOT cfg
  # itself, so <batch> resolves to the config's ENCLOSING directory name
  # (derive_config_run_tag() takes basename() of whatever path it is given,
  # and is a no-op on an extension-less directory name) rather than the
  # config's own basename. That is what disambiguates two configs sharing a
  # basename across different config directories -- e.g. set_A/run_brms.yml
  # and set_B/run_brms.yml resolve to figures/occurrence/set_A/ and
  # figures/occurrence/set_B/ respectively, so they share a path with nothing
  # (the collision test_figure_path_collisions.R guards against).
  # Filenames themselves are unchanged -- figureSaver()
  # below still prefixes every stem with "<tag>_", tag being the config's
  # own basename as before.
  fig_dir <- makeFigureDir(opt = list(figdir = cli_opt$outdir), track = model,
                           config_path = dirname(cfg))

  if (model == "occurrence") {
    .figures_occurrence(opt, tag, fig_dir)
  } else if (model == "duration") {
    .figures_duration(opt, tag, fig_dir)
  } else if (model == "monthly") {
    .figures_monthly(opt, tag, fig_dir)
  } else {
    .figures_size(opt, tag, fig_dir)
  }
}

message("\nDone.")
