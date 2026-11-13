# 01_preprocess_figures_data_overview.R — Figure 1: data overview
#
# Optional figure step of the outbreak-preprocessing phase (inputs are the phase-01
# raw outbreak data, not model fits): outbreak locations by admin level over the
# Africa basemap with lakes and rivers, plus the case-threshold example series.
#   * the outbreak data directory is --indir (default utils/paths.R::
#     defaultOutbreakDataDir(), as 01_preprocess_outbreaks.R);
#   * the Africa polygons and lakes/rivers come from utils/figures_cb.R
#     (load_taxdat_afr_sf(), load_africa_basemap()), shared with 10_manuscript_figures.R;
#   * outputs go to figures/data_overview/ so outbreak_map.png does not collide with
#     the one 01_preprocess_outbreaks.R writes to figures/.
#
# Outputs (--outdir, default figures/data_overview/):
#   outbreak_map.png               — outbreak location periods by admin level (ADM1-ADM3)
#   outbreak_threshold_example.png — one location's weekly cases vs its case threshold
#                                    (skipped if the time-series RDS is absent)
# Survival-model and outbreak-size figures live in 06_fit_figures.R.
#
# Usage:
#   Rscript analysis/01_preprocess_figures_data_overview.R
#   Rscript analysis/01_preprocess_figures_data_overview.R --indir <outbreak data dir> --dpi 200

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(optparse)
  library(taxdat)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")        # defaultOutbreakDataDir()
source("analysis/utils/figures_cb.R")   # load_taxdat_afr_sf(), load_africa_basemap()
source("analysis/utils/figures_io.R")

option_list <- list(
  make_option("--indir", type = "character",
              default = defaultOutbreakDataDir(),   # utils/paths.R, the one data vintage
              help = "Folder with outbreak_shapefiles.gpkg (+ time_series_preoutbreak_extraction.rds) [default: %default]"),
  make_option("--outdir", type = "character", default = "figures/data_overview",
              help = "Figure directory [default: %default]"),
  make_option("--dpi", type = "integer", default = 300L, help = "Raster DPI [default: %default]"),
  make_option("--example_location", type = "character",
              default = "AFR::COD::Sud-Kivu::Uvira Health District",
              help = "Location for the threshold example [default: %default]"),
  make_option("--example_from", type = "character", default = "2016-01-01",
              help = "Earliest week of the threshold example [default: %default]")
)
cli <- parse_args(OptionParser(option_list = option_list))

gpkg <- file.path(cli$indir, "outbreak_shapefiles.gpkg")
if (!file.exists(gpkg)) stop("Outbreak shapefile not found: ", gpkg, " (set --indir).")

# Outbreak location periods: one polygon per location period ---------------------
ob_shp <- sf::st_read(gpkg, quiet = TRUE) %>%
  dplyr::group_by(lctn_pr) %>%
  dplyr::slice(1) %>%
  sf::st_simplify() %>%
  # the extraction occasionally yields a GEOMETRYCOLLECTION (polygon + sliver line/point);
  # keep its polygon part, as st_remove_holes() accepts only (MULTI)POLYGON
  sf::st_collection_extract("POLYGON") %>%
  nngeo::st_remove_holes()

# Basemap: Africa polygons + lakes/rivers fully inside the continent --------------
afr_sf <- load_taxdat_afr_sf()
basemap <- load_africa_basemap(sf::st_union(afr_sf), predicate = sf::st_within)

# The AFR/EMR extraction also holds non-African location periods; the figure (and every
# model dataset) is Africa-only, so keep periods whose ISO3 `country` is on the Africa basemap.
ob_shp <- dplyr::filter(ob_shp, country %in% afr_sf$country)

ob_adm <- dplyr::filter(ob_shp, admin_level != "ADM0")
p_ob <- ggplot(ob_adm) +
  geom_sf(data = afr_sf, lwd = 0.15, fill = "gray95", color = "white") +
  geom_sf(data = sf::st_union(ob_adm), fill = "gray70") +
  geom_sf(data = dplyr::filter(ob_shp, admin_level == "ADM1"), alpha = 1, aes(fill = admin_level)) +
  geom_sf(data = dplyr::filter(ob_shp, admin_level == "ADM2"), alpha = 1, aes(fill = admin_level)) +
  geom_sf(data = dplyr::filter(ob_shp, admin_level == "ADM3"), alpha = 1, aes(fill = admin_level)) +
  geom_sf(inherit.aes = FALSE, data = basemap$rivers_sf,
          color = taxdat:::color_rivers(), lwd = 0.1, alpha = 1) +
  geom_sf(inherit.aes = FALSE, data = basemap$lakes_sf,
          fill = taxdat:::color_lake_fill(), color = taxdat:::color_lake_border(),
          lwd = 0.06, alpha = 1) +
  scale_fill_manual(values = c("#F0C858", "#8F8B1A", "#423E05")) +
  geom_sf(data = afr_sf, lwd = 0.25, color = "black", alpha = 0) +
  theme_bw()
saveFigure(p_ob, "outbreak_map", cli$outdir, width = 10, height = 7, dpi = cli$dpi)

# Threshold example: one location's weekly cases vs its case threshold ----------
ts_path <- file.path(cli$indir, "time_series_preoutbreak_extraction.rds")
if (!file.exists(ts_path)) {
  message("Skipping the threshold example: ", ts_path, " not found.")
} else {
  ts <- readRDS(ts_path)
  ex <- ts %>%
    dplyr::filter(location == cli$example_location, TL > as.Date(cli$example_from))
  # threshold / pop come from the time series when it carries them (newer extractions),
  # else from the shapefile attributes (older extractions).
  missing_cols <- setdiff(c("threshold", "pop"), names(ex))
  if (length(missing_cols))
    ex <- dplyr::inner_join(ex, sf::st_drop_geometry(ob_shp) %>%
                              dplyr::select(location_period_id = lctn_pr, dplyr::all_of(missing_cols)),
                            by = "location_period_id")
  if (!nrow(ex)) {
    message("Skipping the threshold example: no rows for ", cli$example_location)
  } else {
    p_ts <- ggplot(ex, aes(x = TL, y = sCh)) +
      geom_step() +
      geom_hline(aes(yintercept = threshold * pop), lty = 2) +
      theme_bw()
    saveFigure(p_ts, "outbreak_threshold_example", cli$outdir, width = 8, height = 4, dpi = cli$dpi)
  }
}
