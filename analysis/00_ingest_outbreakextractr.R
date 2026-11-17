# 00_ingest_outbreakextractr.R — convert OutbreakExtractR stage-2 outputs to
#   the legacy outbreak files.
#
# Usage:
#   Rscript analysis/00_ingest_outbreakextractr.R
#   Rscript analysis/00_ingest_outbreakextractr.R --indir <dir> --outdir <dir>
#
# Details:
#   Reshapes combined_outbreaks_cholera.csv (+ stage1_geo_*.geojson) from the
#   OutbreakExtractR HPC pipeline into the three files that
#   01_preprocess_outbreaks.R and 02_preprocess_climate.R consume:
#   time_series_preoutbreak_extraction.rds, time_series_outbreak_extraction.rds,
#   outbreak_shapefiles.rds. --indir defaults to default_outbreak_data_dir(),
#   --outdir to --indir. Run before 01_preprocess_outbreaks.R.

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(stringr)
library(lubridate)
library(purrr)
library(furrr)
library(sf)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/paths.R") # default_outbreak_data_dir()
source("analysis/utils/provenance.R") # write_step_counts() (processing cascade, table S1)

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("--indir"),
    default = default_outbreak_data_dir(), # utils/paths.R, the one data vintage
    type = "character",
    help = "Folder containing combined_outbreaks_cholera.csv and stage1_geo/ sub-folder [default: %default]"
  ),
  make_option(
    c("--outdir"),
    default = NULL,
    type = "character",
    help = "Output folder for legacy RDS files [default: same as --indir]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$outdir)) {
  opt$outdir <- opt$indir
}

stopifnot(dir.exists(opt$indir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("---- indir :", opt$indir, "\n")
cat("---- outdir:", opt$outdir, "\n")

# 3. Legacy columns expected by 01_preprocess_outbreaks.R ----------------------

# Population provenance columns: allow downstream QC/gates to distinguish
# pop_source, duplicate-geometry clusters, and the UN-adjustment inputs when the
# extraction emits them.
POP_PROVENANCE_COLS <- c(
  "pop_source",
  "pop_geom_dup_n",
  "pop_geom_dup_class",
  "pop_year_obs",
  "pop_year_raster",
  "pop_natl_ref",
  "adj_factor_flag"
)

PREOUTBREAK_COLS <- c(
  "location",
  "TL",
  "TR",
  "sCh",
  "cCh",
  "deaths",
  "epiweek",
  "start_weekday",
  "spatial_scale",
  "composite_loc",
  "date_range",
  "who_region",
  "country",
  "admin1",
  "admin2",
  "admin3",
  "admin4",
  "admin5",
  "phantom",
  "original_location_name",
  "observation_collection_id",
  "location_period_id",
  "pop",
  "adj_factors",
  POP_PROVENANCE_COLS
)

# Columns consumed by 01_preprocess_outbreaks.R from the outbreak file
OUTBREAK_COLS <- c(
  "location",
  "outbreak_number",
  "location_period_id",
  "country",
  "TL",
  "TR",
  "sCh",
  "cCh",
  "pop",
  "threshold",
  POP_PROVENANCE_COLS
)

# 4a. Read combined stage2 CSV -------------------------------------------------

csv_path <- file.path(opt$indir, "combined_outbreaks_cholera.csv")
stopifnot(file.exists(csv_path))

cat("---- Reading", csv_path, "\n")

raw <- readr::read_csv(
  csv_path,
  col_types = readr::cols(location_period_id = readr::col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    TL = as.Date(TL),
    TR = as.Date(TR)
  )

cat("     Rows read:", nrow(raw), "\n")
n_rows_read <- nrow(raw)

# 4b. Filter to AFR + EMR ------------------------------------------------------

raw <- raw %>%
  filter(who_region %in% c("AFR", "EMR"))

cat("     Rows after AFR/EMR filter:", nrow(raw), "\n")
n_rows_afr_emr <- nrow(raw)

# 4c. Country-scale location / country fix -------------------------------------
# Country-level rows carry "AFR::AGO.Angola" style names.  Apply the Step2
# fix: strip to the code prefix only, and set country to country_iso3 for all.

raw <- raw %>%
  mutate(
    location = case_when(
      spatial_scale == "country" ~ str_extract(location, "^[^.]+"),
      TRUE ~ location
    ),
    # country_iso3 can carry a sub-national suffix (e.g. "TZA::Mainland").
    # Strip everything after the first "::" to always get a plain 3-letter code.
    country = str_extract(country_iso3, "^[A-Z]{3}")
  )

# 4d. Rename adj_factor → adj_factors (legacy column name) ---------------------

if ("adj_factor" %in% names(raw) && !("adj_factors" %in% names(raw))) {
  raw <- raw %>% rename(adj_factors = adj_factor)
}

# 4e. Ensure location_period_id is character (legacy expectation) --------------

raw <- raw %>%
  mutate(location_period_id = as.character(location_period_id))

# 5. Preoutbreak file ----------------------------------------------------------

cat("---- Building time_series_preoutbreak_extraction.rds\n")

# Add composite_loc column if absent (all FALSE)
if (!"composite_loc" %in% names(raw)) {
  raw <- raw %>% mutate(composite_loc = FALSE)
}

# Select only the 24 legacy columns (in order); fill missing ones with NA
preoutbreak <- raw %>%
  select(any_of(PREOUTBREAK_COLS)) %>%
  {
    missing <- setdiff(PREOUTBREAK_COLS, names(.))
    if (length(missing) > 0) {
      cat(
        "  Note: columns absent in source (will be NA):",
        paste(missing, collapse = ", "),
        "\n"
      )
      for (col in missing) {
        .[[col]] <- NA
      }
    }
    .
  } %>%
  select(all_of(PREOUTBREAK_COLS))

cat("     Rows:", nrow(preoutbreak), " | Cols:", ncol(preoutbreak), "\n")

saveRDS(
  preoutbreak,
  file = file.path(opt$outdir, "time_series_preoutbreak_extraction.rds")
)

cat(
  "     Saved →",
  file.path(opt$outdir, "time_series_preoutbreak_extraction.rds"),
  "\n"
)

# 6. Outbreak file -------------------------------------------------------------

cat("---- Building time_series_outbreak_extraction.rds\n")

outbreaks <- raw %>%
  filter(outbreak_number > 0) %>%
  select(any_of(OUTBREAK_COLS)) %>%
  {
    missing <- setdiff(OUTBREAK_COLS, names(.))
    if (length(missing) > 0) {
      cat(
        "  Note: columns absent in source (will be NA):",
        paste(missing, collapse = ", "),
        "\n"
      )
      for (col in missing) {
        .[[col]] <- NA
      }
    }
    .
  } %>%
  select(all_of(OUTBREAK_COLS))

cat("     Rows:", nrow(outbreaks), " | Cols:", ncol(outbreaks), "\n")

saveRDS(
  outbreaks,
  file = file.path(opt$outdir, "time_series_outbreak_extraction.rds")
)

cat(
  "     Saved →",
  file.path(opt$outdir, "time_series_outbreak_extraction.rds"),
  "\n"
)

# 7. Shapefile file ------------------------------------------------------------

cat("---- Building outbreak_shapefiles.rds\n")

geo_dir <- file.path(opt$indir, "stage1_geo")

if (dir.exists(geo_dir)) {
  geojson_files <- list.files(
    geo_dir,
    pattern = "stage1_geo_.*\\.geojson$",
    full.names = TRUE
  )
  # Restrict to AFR/EMR files only — non-AFR/EMR LPs are not in the outbreak
  # data after the who_region filter, so their geometries are never used.
  geojson_files <- geojson_files[
    grepl("stage1_geo_(AFR|EMR)_", basename(geojson_files))
  ]
  cat(
    "     Found",
    length(geojson_files),
    "AFR/EMR stage1_geo_*.geojson files\n"
  )

  if (length(geojson_files) == 0L) {
    stop("No AFR/EMR stage1_geo_*.geojson files found in ", geo_dir)
  }

  # Disable S2 so that st_make_valid() uses GEOS; some geojsons have
  # degenerate edges that S2 rejects but GEOS can repair.
  sf_use_s2(FALSE)
  on.exit(sf_use_s2(TRUE), add = TRUE)

  # availableCores() honours the Slurm/cgroup allocation; detectCores() saw the
  # whole node and oversubscribed a shared compute node.
  n_cores <- min(future::availableCores(), 8L)
  cat("     Reading geojsons in parallel (", n_cores, "cores )\n")
  plan(multisession, workers = n_cores)

  ob_shp <- future_map_dfr(geojson_files, function(f) {
    # Keep only location_period_id + geometry to avoid type-mismatch across
    # geojsons with varying attribute schemas (e.g. deaths: int vs char).
    # First row per LP within the file: the global first-row dedup below keeps
    # the same row, and the duplicates never travel back to the parent.
    sf::st_read(f, quiet = TRUE) %>%
      dplyr::transmute(
        location_period_id = as.character(location_period_id)
      ) %>%
      dplyr::filter(!duplicated(location_period_id))
  }) %>%
    # Filter to unique location periods FIRST, then validate: each LP maps to
    # one shapefile, so st_make_valid()/st_area() run once per geometry instead
    # of once per duplicate row across pull windows. group_by + slice(1) dedups
    # on the character key without the slow full-row geometry comparison that
    # distinct(.keep_all = TRUE) performs on sf objects (the 4h+ bottleneck).
    dplyr::group_by(location_period_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # st_make_valid()/st_area() are per-feature, so contiguous chunks on the
  # workers give the serial result. Workers start with s2 on: switch it off
  # there too, or they would repair with s2 instead of GEOS.
  chunk_id <- cut(seq_len(nrow(ob_shp)), breaks = n_cores, labels = FALSE)
  ob_shp <- future_map_dfr(split(ob_shp, chunk_id), function(x) {
    sf::sf_use_s2(FALSE)
    x <- sf::st_make_valid(x)
    x$area_per_1km2 <- as.numeric(sf::st_area(x)) / 1e6
    x
  }) %>%
    dplyr::mutate(lctn_pr = as.character(location_period_id)) %>%
    dplyr::select(lctn_pr, area_per_1km2, geometry)
  plan(sequential)
} else {
  # Fallback: accept a single pre-merged gpkg named stage1_geo.gpkg
  gpkg_path <- file.path(opt$indir, "stage1_geo.gpkg")
  if (!file.exists(gpkg_path)) {
    stop(
      "No stage1_geo/ sub-folder and no stage1_geo.gpkg found in ",
      opt$indir,
      "\nStage geometry into one of those forms."
    )
  }
  cat("     Using pre-merged", gpkg_path, "\n")
  ob_shp <- st_read(gpkg_path, quiet = TRUE) %>%
    transmute(location_period_id = as.character(location_period_id)) %>%
    dplyr::group_by(location_period_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    st_make_valid() %>%
    mutate(
      lctn_pr = as.character(location_period_id),
      area_per_1km2 = as.numeric(st_area(geometry)) / 1e6
    ) %>%
    select(lctn_pr, area_per_1km2, geometry)
}

cat("     LPs in shapefile:", nrow(ob_shp), "\n")

saveRDS(ob_shp, file = file.path(opt$outdir, "outbreak_shapefiles.rds"))

cat("     Saved →", file.path(opt$outdir, "outbreak_shapefiles.rds"), "\n")

# 8. Processing cascade (table S1) ---------------------------------------------
# Processing cascade (table S1): counts at each output of this step. Distinct
# counts only where the column exists in the source.
n_dist <- function(d, ...) {
  cols <- c(...)
  if (!all(cols %in% names(d))) {
    return(NA_integer_)
  }
  nrow(dplyr::distinct(d, dplyr::across(dplyr::all_of(cols))))
}
write_step_counts(
  "00_ingest_outbreakextractr",
  list(
    n_rows_read = n_rows_read,
    n_rows_afr_emr = n_rows_afr_emr,
    n_rows_preoutbreak = nrow(preoutbreak),
    n_lp_preoutbreak = n_dist(preoutbreak, "location_period_id"),
    n_countries_preoutbreak = n_dist(preoutbreak, "country"),
    n_rows_outbreak = nrow(outbreaks),
    n_outbreaks = n_dist(outbreaks, "location", "outbreak_number"),
    n_lp_outbreak = n_dist(outbreaks, "location_period_id"),
    n_countries_outbreak = n_dist(outbreaks, "country"),
    n_lp_shapefile = nrow(ob_shp)
  ),
  indir = opt$indir
)

cat("---- Done. Legacy files written to", opt$outdir, "\n")
