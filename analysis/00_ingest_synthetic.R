# 00_ingest_synthetic.R — write the raw inputs of a synthetic world.
#
# Phase 00, variant `synthetic`: a sibling of 00_ingest_outbreakextractr.R that,
# instead of converting the restricted outbreak database, writes a small FAKE
# world in the formats of every external input, so the whole pipeline runs on a
# machine without the real data (analysis/bash/tools/run_synthetic_pipeline.sh).
# Countries, admin units and their names are invented
# (analysis/groupings/synthetic_countries.csv, ISO 3166 user-assigned codes
# XA?), the polygons lie over open ocean, and every series is simulated: no
# real record, name or boundary is written. The model and its known ENSO
# effects are documented in utils/synthetic_world.R.
#
# Run it from the root of a sandbox (its own data/ and generated_data/) with
# CHOLERA_OUTBREAK_DATA_DIR pointing inside the sandbox. It overwrites real-input
# paths (data/geodata, data/who_annual_report_africa, the ENSO caches, ...), so it
# refuses to run unless the working directory holds the marker file
# .synthetic_sandbox (written by the driver) and is not a git checkout, and it
# refuses to write into ~/secure/. Then run the real 00 -> 03 steps (the driver
# does it all).
#
# Writes (relative to the working directory):
#   $CHOLERA_OUTBREAK_DATA_DIR/<vintage>/combined_outbreaks_cholera.csv + stage1_geo/*.geojson
#   data/geodata/{gadm_410-levels,gadm_valid}.gpkg, Africa_waterbody.shp, rivers_africa_37333.shp
#   data/climate_parquet/<ISO3>/Africa_<core>.parquet.gzip (+ DOWNLOAD_MANIFEST.tsv)
#   data/worldpop/*.tif, data/country_boundaries/<ISO3>_adm0.geojson
#   data/who_annual_report_africa/{who_annual_report_africa,who_dashboard_annual_africa,
#                                  who_dashboard_weekly_africa}.rds, who_sitrep_cumulative_africa.csv
#   data/best_seasonality_estimates.csv, data/hmm_country_clusters.csv
#   data/lenssen/extracted/<model>/step01_InitVerifInfo.Rda, data/precalc_plume_quantiles_ENSO_NINO34_synthetic.rds
#   generated_data/enso_{oni,nino34,nino34_oisst}_cache.rds (+ _download_meta.yml), enso_cfsv2_cache.rds
#   generated_data/wpp2024_national_totals.csv
#   generated_data/synthetic_world/truth.yml, epidemics.rds   (the simulated truth)
#   generated_data/synthetic_world/enso_events_synth.yml       (past ENSO events of the synthetic ONI)
#
# Usage (inside a sandbox; the driver sets it up):
#   bash analysis/bash/tools/run_synthetic_pipeline.sh --sandbox synthetic_world
#   cd synthetic_world && CHOLERA_OUTBREAK_DATA_DIR=$PWD/secure \
#     Rscript analysis/00_ingest_synthetic.R --seed 2 --force

suppressPackageStartupMessages({
  library(optparse)
})
for (m in c("base", "gid", "paths", "enso", "clim_compile", "synthetic_world"))
  source(sprintf("analysis/utils/%s.R", m))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--registry", default = "analysis/groupings/synthetic_countries.csv",
              help = "synthetic country registry [default %default]"),
  make_option("--seed", default = 1L, type = "integer", help = "random seed [default %default]"),
  make_option("--adm1_per_country", default = 4L, type = "integer", help = "[default %default]"),
  make_option("--adm2_per_adm1", default = 2L, type = "integer", help = "[default %default]"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "overwrite an existing synthetic world"))))

if (!file.exists(".synthetic_sandbox") || file.exists(".git"))
  stop("[00_synthetic] run from a synthetic sandbox (a directory holding .synthetic_sandbox, not a ",
       "git checkout): this script overwrites data/geodata, data/who_annual_report_africa, the ",
       "ENSO caches and other real inputs. Use analysis/bash/tools/run_synthetic_pipeline.sh.",
       call. = FALSE)
root <- sub("/+$", "", SECURE_OUTBREAK_DIR)
if (!nzchar(Sys.getenv("CHOLERA_OUTBREAK_DATA_DIR")) ||
    normalizePath(root, mustWork = FALSE) == normalizePath("~/secure/cholera_outbreak_data", mustWork = FALSE))
  stop("[00_synthetic] set CHOLERA_OUTBREAK_DATA_DIR to a folder of this sandbox: the synthetic ",
       "world never writes into ~/secure/.", call. = FALSE)
vintage_dir <- path.expand(defaultOutbreakDataDir())
if (file.exists(file.path(vintage_dir, "combined_outbreaks_cholera.csv")) && !opt$force)
  stop("[00_synthetic] ", vintage_dir, " already holds a synthetic database; pass --force.", call. = FALSE)
say <- function(...) message("[00_synthetic] ", sprintf(...))
dirs <- c(vintage_dir, "data/geodata", "data/climate_parquet", "data/worldpop", "data/country_boundaries",
          "data/who_annual_report_africa", "generated_data/synthetic_world")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

reg  <- synthRegistry(opt$registry)
gadm <- synthGadm(reg, n_adm1 = opt$adm1_per_country, n_adm2 = opt$adm2_per_adm1)
units <- synthUnits(reg, gadm, seed = opt$seed)
say("%d countries, %d ADM1, %d ADM2, %d reporting units", nrow(reg), nrow(gadm$adm1),
    nrow(gadm$adm2), sum(units$level != "ADM0"))

# --- GADM and basemap ----------------------------------------------------------
for (f in c("data/geodata/gadm_410-levels.gpkg", "data/geodata/gadm_valid.gpkg")) {
  if (file.exists(f)) file.remove(f)
  sf::st_write(gadm$adm1, f, layer = "ADM_1", quiet = TRUE)
  sf::st_write(gadm$adm2, f, layer = "ADM_2", quiet = TRUE, append = TRUE)
}
# basemap layers with the attributes taxdat::get_lakes() / get_rivers() filter on
lake <- sf::st_sf(NAME = "Synthetic Lake", Shape_area = 0.25, geometry = sf::st_sfc(
  .msquare(reg$lon0[1] + 0.25, reg$lat0[1] + 1.25, 0.5, 0.5), crs = 4326))
river <- sf::st_sf(NAME = "Synthetic River", A_Strahler = 1L, geometry = sf::st_sfc(
  sf::st_linestring(rbind(c(reg$lon0[1] + 0.1, reg$lat0[1] + 0.1), c(reg$lon0[1] + 1.9, reg$lat0[1] + 1.9))), crs = 4326))
sf::st_write(lake, "data/geodata/Africa_waterbody.shp", delete_dsn = file.exists("data/geodata/Africa_waterbody.shp"), quiet = TRUE)
sf::st_write(river, "data/geodata/rivers_africa_37333.shp", delete_dsn = file.exists("data/geodata/rivers_africa_37333.shp"), quiet = TRUE)

# --- ENSO ----------------------------------------------------------------------
n34 <- synthNino34(seed = opt$seed + 6L)
oni <- synthOni(n34)
nino34 <- tibble::tibble(year = n34$year, month = n34$month, enso_index = round(n34$nino34, 2),
                         date = n34$date, enso_phase = NA_character_)
saveRDS(oni, "generated_data/enso_oni_cache.rds")
saveRDS(nino34, "generated_data/enso_nino34_cache.rds")
saveRDS(nino34, "generated_data/enso_nino34_oisst_cache.rds")
for (s in c("oni", "nino34", "nino34_oisst"))
  yaml::write_yaml(list(source = s, url = "synthetic (00_ingest_synthetic.R)",
                        downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                   sprintf("generated_data/enso_%s_download_meta.yml", s))
# CFSv2: one issuance at the forecast handoff (the month after the last ONI
# month), 40 members persisting the last anomaly with damped noise, 10 monthly
# targets (lead 0..9; the real cache ends 2027-06, the continuation library takes over).
init <- seq(max(oni$date), by = "month", length.out = 2)[2]
set.seed(opt$seed + 21L)
cfs <- do.call(rbind, lapply(1:40, function(mb) {
  x <- utils::tail(n34$nino34, 1); out <- numeric(10)
  for (l in 1:10) { x <- 0.9 * x + stats::rnorm(1, 0, 0.2); out[l] <- x }
  data.frame(init_date = init, target_date = seq(init, by = "month", length.out = 10),
             lead = 0:9, member = mb, nino34_anomaly = round(out, 3))
}))
saveRDS(tibble::as_tibble(cfs), "generated_data/enso_cfsv2_cache.rds")
yaml::write_yaml(list(url = "synthetic (00_ingest_synthetic.R)",
                      downloaded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), variable = "anom",
                      dims = list(time = 10L, ens = 40L), init = format(init),
                      target_range = format(range(cfs$target_date)), n_members = 40L, nrow = nrow(cfs)),
                 "generated_data/enso_cfsv2_download_meta.yml")
# C3S plume stand-in: quantiles of the CFSv2 members over the first 6 leads
# (readEnsoC3SPlume() reads the companion .rds of the .png path the deploy names).
plume <- do.call(rbind, lapply(split(cfs[cfs$lead %in% 1:6, ], cfs$target_date[cfs$lead %in% 1:6]), function(d)
  data.frame(date = d$target_date[1], q_level = c(0.1, 0.25, 0.5, 0.75, 0.9),
             nino34_anomaly = round(as.numeric(stats::quantile(d$nino34_anomaly, c(0.1, 0.25, 0.5, 0.75, 0.9))), 3))))
saveRDS(tibble::as_tibble(plume), "data/precalc_plume_quantiles_ENSO_NINO34_synthetic.rds")
# Lenssen piControl stand-ins: one long monthly series per model of LENSSEN_MODELS.
for (k in seq_along(LENSSEN_MODELS)) {
  md <- file.path("data/lenssen/extracted", LENSSEN_MODELS[k])
  dir.create(md, recursive = TRUE, showWarnings = FALSE)
  pc <- synthNino34(start = "1000-01-01", end = "1199-12-01", seed = 100L + k)
  nino34 <- as.numeric(stats::filter(pc$nino34, rep(1 / 3, 3), sides = 2))
  timeMap <- cbind(pc$year, pc$month)
  save(nino34, timeMap, file = file.path(md, "step01_InitVerifInfo.Rda"))
}
# Past ENSO events of the synthetic ONI (deploy.event_attribution.registry of synth_run1)
events <- synthEnsoEvents(oni)
yaml::write_yaml(list(events = events), "generated_data/synthetic_world/enso_events_synth.yml")
say("ENSO events: %s", paste(names(events), collapse = ", "))
reg_ev <- "analysis/configs/enso_events_synth.yml"
if (file.exists(reg_ev) && !identical(yaml::read_yaml(reg_ev)$events, yaml::read_yaml(
    "generated_data/synthetic_world/enso_events_synth.yml")$events))
  warning("[00_synthetic] the synthetic ONI's events differ from ", reg_ev, " (another --seed?): ",
          "point deploy.event_attribution.registry at generated_data/synthetic_world/enso_events_synth.yml",
          call. = FALSE)
say("ENSO: ONI %s .. %s, CFSv2 init %s, %d piControl stand-ins", min(oni$date), max(oni$date),
    init, length(LENSSEN_MODELS))

# --- Epidemics and the outbreak database ----------------------------------------
epi <- synthEpidemics(units, reg, oni, gadm)
weekly <- synthWeekly(units, epi)
csv <- synthOutbreakCsv(weekly, units)
utils::write.csv(csv, file.path(vintage_dir, "combined_outbreaks_cholera.csv"), row.names = FALSE, na = "")
geo <- synthUnitGeometry(units, gadm)
gd <- file.path(vintage_dir, "stage1_geo"); dir.create(gd, showWarnings = FALSE)
gp <- file.path(gd, "stage1_geo_AFR_synthetic.geojson")
sf::st_write(geo, gp, delete_dsn = file.exists(gp), quiet = TRUE)
n_ep <- sum(epi$onset[epi$month <= as.Date("2024-12-01") & epi$month >= as.Date("2010-01-01")])
say("database: %d weekly rows, %d outbreaks 2010-2024, %s cases", nrow(csv), n_ep,
    format(sum(weekly$sCh[weekly$outbreak_number > 0]), big.mark = ","))
saveRDS(epi, "generated_data/synthetic_world/epidemics.rds")
yaml::write_yaml(list(
  registry = opt$registry, seed = opt$seed,
  groups = lapply(split(reg, reg$group), function(d)
    list(countries = d$iso3, enso_lag = d$enso_lag[1], enso_beta = d$enso_beta[1],
         season_peak = d$season_peak[1])),
  note = "negative enso_lag = ENSO leads; the onset log-hazard carries enso_beta * ONI(t + enso_lag)"),
  "generated_data/synthetic_world/truth.yml")

# --- Climate (one parquet per ADM1 and ADM2 unit) -------------------------------
gids <- c(gadm$adm1$GID_1, gadm$adm2$GID_2)
for (i in seq_along(gids)) {
  iso <- substr(gids[i], 1, 3); r <- reg[reg$iso3 == iso, ]
  dd <- file.path("data/climate_parquet", iso); dir.create(dd, showWarnings = FALSE)
  cl <- synthClimateDaily(oni, r$rain_peak, r$enso_beta, seed = opt$seed * 1000L + i)
  arrow::write_parquet(cl, file.path(dd, paste0("Africa_", makeClimParquetCore(gids[i]))),
                       compression = "gzip")
}
writeLines(c("file\tsource", "synthetic\t00_ingest_synthetic.R"),
           "data/climate_parquet/DOWNLOAD_MANIFEST.tsv")
say("climate: %d daily parquet files 1981-2024", length(gids))

# --- Population ------------------------------------------------------------------
for (iso in reg$iso3) {
  b <- sf::st_sf(shapeISO = iso, shapeName = reg$name[reg$iso3 == iso],
                 geometry = sf::st_union(sf::st_geometry(gadm$adm1[gadm$adm1$GID_0 == iso, ])))
  f <- sprintf("data/country_boundaries/%s_adm0.geojson", iso)
  sf::st_write(b, f, delete_dsn = file.exists(f), quiet = TRUE)
  for (yr in 2015:2030)
    terra::writeRaster(synthWorldPopRaster(units, gadm, iso, yr),
                       sprintf("data/worldpop/%s_pop_%d_CN_100m_R2025A_v1.tif", tolower(iso), yr),
                       overwrite = TRUE)
}
natpop <- units[units$level == "ADM0", ]
wpp <- do.call(rbind, lapply(2000:2030, function(yr)
  data.frame(iso3 = natpop$iso3, year = yr, pop_total = round(natpop$pop * 1.025^(yr - 2020)))))
utils::write.csv(wpp, "generated_data/wpp2024_national_totals.csv", row.names = FALSE)
say("population: %d rasters, WPP-like national totals 2000-2030", nrow(reg) * 16L)

# --- WHO, seasonality, clusters --------------------------------------------------
who <- synthWhoTables(units, epi, reg)
saveRDS(who$annual, "data/who_annual_report_africa/who_annual_report_africa.rds")
saveRDS(who$dashboard_annual, "data/who_annual_report_africa/who_dashboard_annual_africa.rds")
saveRDS(who$dashboard_weekly, "data/who_annual_report_africa/who_dashboard_weekly_africa.rds")
utils::write.csv(who$sitrep, "data/who_annual_report_africa/who_sitrep_cumulative_africa.csv", row.names = FALSE)
utils::write.csv(synthSeasonality(gadm, reg), "data/best_seasonality_estimates.csv", row.names = FALSE)
utils::write.csv(data.frame(cluster = reg$hmm_cluster, country = reg$name, iso3 = reg$iso3),
                 "data/hmm_country_clusters.csv", row.names = FALSE)
say("WHO: annual 2010-2022, dashboard %s .. %s", min(who$dashboard_weekly$date_wk), max(who$dashboard_weekly$date_wk))
say("done: the synthetic world's raw inputs are in place (truth in generated_data/synthetic_world/truth.yml)")
