# 02_preprocess_climate.R — per-gid monthly climate statistics and extremes
#   from CHIRPS precipitation and ERA5-Land temperature.
#
# Usage:
#   Rscript analysis/02_preprocess_climate.R -s TRUE
#   Rscript analysis/02_preprocess_climate.R --stage match
#   Rscript analysis/02_preprocess_climate.R --stage extract
#   Rscript analysis/02_preprocess_climate.R --stage db -s TRUE
#
# Details:
#   Inputs: data/climate_parquet/ (--climate_source parquet) or the legacy
#   per-country CSV corpus data/climate_data/ (csv); outbreak shapefiles under
#   --indir for the gid universe. Outputs: generated_data/climate_stats.sqlite
#   (table climate_stats) and climate_extremes.sqlite (clim_extremes_<baseline>)
#   (--stats_db / --extremes_db), over a configurable baseline (-b).
#   The defaults are the deployed build: a bare command rebuilds the deployed
#   database; writing one whose provenance records other settings needs --force.

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(glue)
library(sf)
library(furrr)
library(RSQLite)
library(DBI)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R") # transitive dependency of siz_data.R (utils/README.md)

# 2. User-supplied options -----------------------------------------------------
option_list <- list(
  make_option(
    c("-r", "--redo_gadm"),
    default = FALSE,
    action = "store",
    type = "logical",
    help = "redo data processing"
  ),
  make_option(
    c("-s", "--redo_stats"),
    default = FALSE,
    action = "store",
    type = "logical",
    help = "redo data processing"
  ),
  make_option(
    c("-b", "--baseline"),
    default = "1990_2019",
    action = "store",
    type = "character",
    help = "baseline period"
  ),
  # The data leg's one vintage (utils/paths.R::default_outbreak_data_dir()).
  make_option(
    "--indir",
    default = default_outbreak_data_dir(),
    type = "character",
    help = "Folder containing outbreak_shapefiles.gpkg [default: %default]"
  ),
  # --- gid universe source ---
  # The centroid match below reduces every location period to one point, so a
  # composite LP (a union of K admin units, e.g.
  #   EMR::SDN::Gedaref|Gezira|Kassala|Khartoum|Blue Nile|River Nile|Red Sea
  # = 7 states) contributes exactly one gid to `u_gids`. Its other 6 members are
  # never extracted, so they have no climate row -- and make_size_data()'s
  # inner_join (siz_data.R) then SILENTLY DROPS them once the
  # crosswalk starts apportioning cases to them.
  #
  # `crosswalk` unions the centroid universe with the crosswalk's member gids so
  # every apportionment target has climate. It is ADDITIVE, and the default
  # (the deployed build).
  make_option(
    "--gid_source",
    default = "crosswalk",
    action = "store",
    type = "character",
    help = paste(
      "gid universe for climate extraction: 'centroid'",
      "(legacy, one gid per LP) or 'crosswalk' (union",
      "with u_gids_crosswalk.rds) [default: %default]"
    )
  ),
  make_option(
    "--crosswalk_gids_path",
    default = "generated_data/u_gids_crosswalk.rds",
    action = "store",
    type = "character",
    help = paste(
      "gid universe written by",
      "02_preprocess_crosswalk.R [default: %default]"
    )
  ),
  # Costing gate. dbWriteTable(..., overwrite = TRUE) below means a rerun
  # rebuilds the whole climate_stats table, not a delta -- so the cost of
  # switching to `crosswalk` is the full extraction, not the new gids alone.
  # This flag reports the diff and exits before doing any of it.
  make_option(
    "--cost_only",
    default = FALSE,
    action = "store",
    type = "logical",
    help = paste(
      "report the gid-universe diff vs the existing",
      "climate_stats table and exit without extracting",
      "[default: %default]"
    )
  ),
  # --- raw climate source ---------------------------------------------------
  # `csv`     : the per-country/per-variable/per-decade CSV corpus in
  #             data/climate_data/ (legacy; covers 29 countries, and its Ghana
  #             rows are cartesian-inflated -- see
  #             clim_compile.R::compile_stats()).
  # `parquet` : the Jahn et al. (2026) re-release of the same upstream product,
  #             one file per admin unit, 1981-01-01..2024-12-31, all of Africa.
  #             Fetch with analysis/bash/tools/download_jahn2026_parquet.sh.
  # `parquet` is the default (the deployed build); the CSV
  # path stays runnable for old-vs-new comparison.
  make_option(
    "--climate_source",
    default = "parquet",
    action = "store",
    type = "character",
    help = paste(
      "raw climate source: 'csv' (legacy",
      "data/climate_data/) or 'parquet'",
      "[default: %default]"
    )
  ),
  make_option(
    "--climate_parquet_dir",
    default = "data/climate_parquet/",
    action = "store",
    type = "character",
    help = "root of the extracted parquet corpus [default: %default]"
  ),
  make_option(
    "--weighting",
    default = "dynamic",
    action = "store",
    type = "character",
    help = paste(
      "parquet spatial aggregation: 'simple'",
      "(unweighted), 'static' (pop 2000) or 'dynamic'",
      "(pop, refreshed every 10y -- closest analogue to",
      "the CSV era's *_pop_2010) [default: %default]"
    )
  ),
  # Staging targets. A rebuild is destructive (dbWriteTable(overwrite = TRUE)
  # rebuilds the whole table), so these exist to let a parquet build land in a
  # scratch database and be diffed against the live one first.
  make_option(
    "--stats_db",
    default = "generated_data/climate_stats.sqlite",
    action = "store",
    type = "character",
    help = "climate_stats sqlite path [default: %default]"
  ),
  make_option(
    "--extremes_db",
    default = "generated_data/climate_extremes.sqlite",
    action = "store",
    type = "character",
    help = "climate_extremes sqlite path [default: %default]"
  ),
  make_option(
    "--stage",
    default = "all",
    action = "store",
    type = "character",
    help = paste(
      "all | match | extract | db: match writes",
      "matched_ob_gids.rds only; extract fills the per-gid",
      "caches; db writes the tables from the caches",
      "[default: %default]"
    )
  ),
  make_option(
    "--force",
    default = FALSE,
    action = "store",
    type = "logical",
    help = paste(
      "write a database whose provenance records other",
      "source / gid / weighting / baseline settings",
      "[default: %default]"
    )
  ),
  make_option(
    "--workers",
    default = future::availableCores(),
    action = "store",
    type = "integer",
    help = paste(
      "furrr workers; the default honours the Slurm/cgroup",
      "allocation [default: %default]"
    )
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

# --stage splits the work: `match` (the outbreak-vintage-dependent centroid
# match -> matched_ob_gids.rds), `extract` (per-gid caches, independent of the
# vintage except through the gid set) and `db` (caches -> tables); `all`
# (default) runs the three in one pass.
if (!opt$stage %in% c("all", "match", "extract", "db")) {
  stop("--stage must be one of all, match, extract, db; got: ", opt$stage)
}

# The defaults are the deployed build (--climate_source parquet, --gid_source
# crosswalk, --weighting dynamic; see generated_data/climate_db_provenance.yml),
# so a bare command rebuilds the deployed database. Writing a database whose
# provenance records other source / gid / weighting settings is refused without
# --force. A staged build (--stats_db elsewhere) writes its provenance beside
# its database, not over generated_data/climate_db_provenance.yml.
# The deployed databases' provenance; a staged build keeps its own beside it.
DEFAULT_PROV_PATH <- "generated_data/climate_db_provenance.yml"
PROV_PATH <- if (
  identical(opt$stats_db, "generated_data/climate_stats.sqlite")
) {
  DEFAULT_PROV_PATH
} else {
  sub("\\.sqlite$", "_provenance.yml", opt$stats_db)
}

# 3. Validate the raw climate source -------------------------------------------
# Deliberately fatal rather than a silent fall back to the CSV corpus: falling
# back would rebuild the whole table from the 29-country CSV source under a
# `--climate_source parquet` command line, quietly reinstating both the coverage
# gap and the Ghana cartesian blow-up the switch exists to fix.
if (!opt$climate_source %in% c("csv", "parquet")) {
  stop("--climate_source must be 'csv' or 'parquet', got: ", opt$climate_source)
}

if (identical(opt$climate_source, "parquet")) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("--climate_source parquet requires the `arrow` package.")
  }

  if (!opt$weighting %in% c("simple", "static", "dynamic")) {
    stop(
      "--weighting must be 'simple', 'static' or 'dynamic', got: ",
      opt$weighting
    )
  }

  iso_dirs <- list.dirs(opt$climate_parquet_dir, recursive = FALSE)
  if (!dir.exists(opt$climate_parquet_dir) || length(iso_dirs) == 0) {
    stop(
      "--climate_source parquet was requested but no parquet corpus was ",
      "found at: ",
      opt$climate_parquet_dir,
      "\n",
      "  Run `bash analysis/bash/tools/download_jahn2026_parquet.sh` first, or pass ",
      "--climate_parquet_dir."
    )
  }

  cat(sprintf(
    "---- climate source: parquet (%s), weighting '%s', %d country dir(s)\n",
    opt$climate_parquet_dir,
    opt$weighting,
    length(iso_dirs)
  ))
} else {
  cat("---- climate source: csv (legacy data/climate_data/)\n")
}

valid_gadm_db <- "data/geodata/gadm_valid.gpkg"

# Session for furrr (was a fixed 10 workers, whatever the allocation)
plan(multisession, workers = opt$workers)

# 4. Join spatial data with GADM -----------------------------------------------
# Aim is to join the outbreak shapefiles to the GADM shapefiles. Here all admin
# level shapefiels lower than 1 are linked to the GADM adm 2 level, and admin
# 1 level is joined to the admin 1 level. We exclued country-level outbreaks for now.

if (opt$stage %in% c("all", "match")) {
  # Shapefiles of outbreaks
  ob_sf <- st_read(file.path(opt$indir, "outbreak_shapefiles.gpkg"))

  # GADM shapefiles for admin levels 1 and 2
  if (opt$redo_gadm) {
    # Raw GADM database
    gadm_db <- "data/geodata/gadm_410-levels.gpkg"

    # ADM1 level
    gadm_adm1 <- st_read(gadm_db, layer = "ADM_1") %>%
      st_make_valid()
    # Write for reuse
    st_write(gadm_adm1, valid_gadm_db, layer = "ADM_1")

    # ADM2 level
    gadm_adm2 <- st_read(gadm_db, layer = "ADM_2") %>%
      st_make_valid()
    # Write for reuse
    st_write(gadm_adm2, valid_gadm_db, layer = "ADM_2")
  } else {
    # Only the outbreak countries: the match below filters on GID_0 anyway, and
    # the whole-continent layers dominated this step's read time. ORDER BY rowid
    # keeps the order of a full-layer read.
    ob_iso3 <- unique(ob_sf$country[!is.na(ob_sf$country)])
    read_gadm <- function(layer) {
      st_read(
        valid_gadm_db,
        query = sprintf(
          "SELECT * FROM \"%s\" WHERE GID_0 IN (%s) ORDER BY rowid",
          layer,
          paste0(
            "'",
            gsub("'", "''", ob_iso3, fixed = TRUE),
            "'",
            collapse = ", "
          )
        )
      )
    }
    gadm_adm1 <- read_gadm("ADM_1")
    gadm_adm2 <- read_gadm("ADM_2")
  }

  # Intersect outbreaks on centroids
  ob_cnt <- ob_sf %>%
    filter(admin_level != "ADM0") %>%
    group_by(lctn_pr) %>%
    slice(1) %>%
    mutate(
      admin_level_simple = case_when(
        admin_level == "ADM1" ~ "ADM1",
        TRUE ~ "ADM2"
      )
    ) %>%
    st_centroid()

  # Match by centroid
  matched_gadm <- ob_cnt %>%
    group_by(admin_level_simple, country) %>%
    group_modify(function(x, y) {
      print(y)

      if (y$admin_level_simple[1] == "ADM1") {
        ga1 <- filter(gadm_adm1, GID_0 == y$country[1])
        matched_gadm <- st_join(x, ga1, join = st_intersects)
      } else {
        ga2 <- filter(gadm_adm2, GID_0 == y$country[1])
        matched_gadm <- st_join(x, ga2, join = st_intersects)
      }

      matched_gadm %>%
        mutate(lctn_pr = x$lctn_pr)
    }) %>%
    ungroup()

  # Define the gid of the matched shape
  matched_gadm <- matched_gadm %>%
    mutate(gid = case_when(admin_level_simple == "ADM2" ~ GID_2, TRUE ~ GID_1))

  rm(gadm_adm1, gadm_adm2, ob_cnt)
  gc()

  saveRDS(matched_gadm, file = "generated_data/matched_ob_gids.rds")

  if (identical(opt$stage, "match")) {
    cat("---- --stage match: wrote generated_data/matched_ob_gids.rds\n")
    quit(save = "no", status = 0)
  }
} else {
  # extract / db: the gid universe comes from the match already on disk.
  if (!file.exists("generated_data/matched_ob_gids.rds")) {
    stop(
      "--stage ",
      opt$stage,
      " needs generated_data/matched_ob_gids.rds; ",
      "run --stage match first."
    )
  }
  matched_gadm <- readRDS("generated_data/matched_ob_gids.rds")
}

# 5. Get climate data ----------------------------------------------------------

# Unique gids
u_gids <- matched_gadm %>%
  filter(!is.na(gid)) %>%
  pull(gid) %>%
  unique()

n_gids_centroid <- length(u_gids)

# 5a. Expand the gid universe with the crosswalk members -----------------------
# Additive union only. `u_gids_crosswalk.rds` is written by the sibling
# script 02_preprocess_crosswalk.R, which resolves each LP to all of its
# member gids by population overlap instead of to one centroid.
if (identical(opt$gid_source, "crosswalk")) {
  if (!file.exists(opt$crosswalk_gids_path)) {
    stop(
      "--gid_source crosswalk was requested but the crosswalk gid universe ",
      "was not found at: ",
      opt$crosswalk_gids_path,
      "\n",
      "  Run `Rscript analysis/02_preprocess_crosswalk.R` first, or pass ",
      "--crosswalk_gids_path. This is deliberately fatal rather than a ",
      "silent fall back to the centroid universe: falling back would leave ",
      "composite members without climate, and make_size_data()'s inner_join ",
      "would then drop their apportioned cases without a word."
    )
  }

  xw_gids <- readRDS(opt$crosswalk_gids_path)
  # Accept either a bare character vector or a frame with a `gid` column.
  if (is.data.frame(xw_gids)) {
    xw_gids <- xw_gids$gid
  }
  xw_gids <- unique(xw_gids[!is.na(xw_gids) & nzchar(xw_gids)])
  # 02c writes the crosswalk's gids FLATTENED ("AGO_4_1") while the centroid
  # match above carries GADM's own spelling ("AGO.4_1"): compared as strings,
  # every shared gid would enter u_gids twice (and assert_unique_clim_rows() then
  # refuses the table), and the parquet lookup needs GADM's spelling anyway.
  # Map back through GADM's gid columns (no geometry): GADM spells some gids
  # irregularly (e.g. "GHA14.11_2"), so no string rule inverts flatten_gid().
  gadm_ids <- c(
    st_read(
      valid_gadm_db,
      query = 'SELECT GID_1 AS gid FROM "ADM_1"',
      quiet = TRUE
    )$gid,
    st_read(
      valid_gadm_db,
      query = 'SELECT GID_2 AS gid FROM "ADM_2"',
      quiet = TRUE
    )$gid
  )
  gadm_ids <- gadm_ids[!is.na(gadm_ids)]
  spelled <- stats::setNames(gadm_ids, flatten_gid(gadm_ids))
  xw_spelt <- unname(spelled[flatten_gid(xw_gids)])
  if (anyNA(xw_spelt)) {
    warning(
      sum(is.na(xw_spelt)),
      " crosswalk gid(s) not in ",
      valid_gadm_db,
      ", kept as written: ",
      paste(utils::head(xw_gids[is.na(xw_spelt)], 5), collapse = ", ")
    )
  }
  xw_gids <- unique(ifelse(is.na(xw_spelt), xw_gids, xw_spelt))

  new_gids <- setdiff(xw_gids, u_gids)
  lost_gids <- setdiff(u_gids, xw_gids) # centroid-only; kept by the union

  u_gids <- union(u_gids, xw_gids)
  # One spelling per gid, whatever the inputs: the tables are keyed flattened.
  u_gids <- u_gids[!duplicated(flatten_gid(u_gids))]

  cat(sprintf(
    paste0(
      "---- gid universe: centroid %d + crosswalk %d -> union %d\n",
      "     new gids needing extraction : %d\n",
      "     centroid-only gids retained : %d\n"
    ),
    n_gids_centroid,
    length(xw_gids),
    length(u_gids),
    length(new_gids),
    length(lost_gids)
  ))
}

# 5b. Cost the expansion against what is already in the database ---------------
# A rebuild is destructive (dbWriteTable(overwrite = TRUE) rebuilds the whole
# table, not a delta); --cost_only reports the gid-universe diff and exits
# before extracting.
# Read-only. Reports how much of `u_gids` the existing climate_stats table
# already covers, and how many gids have no raw climate files at all (those
# would return NULL from compile_all_stats() and never reach the table).
{
  covered <- character(0)
  # Guard on the file: dbConnect() would otherwise CREATE an empty sqlite file
  # as a side effect of what is meant to be a read-only costing pass.
  if (file.exists(opt$stats_db)) {
    stats_db <- get_stats_db(opt$stats_db)
    if ("climate_stats" %in% DBI::dbListTables(stats_db)) {
      covered <- DBI::dbGetQuery(
        stats_db,
        "SELECT DISTINCT gid FROM climate_stats"
      )$gid
    }
    DBI::dbDisconnect(stats_db)
  }

  # climate_stats stores FLATTENED gids ("AGO_4_1"); u_gids (from the GADM
  # join / crosswalk union) carries the raw dotted GADM form ("AGO.4_1").
  # Comparing the two forms directly makes every gid look absent -- compare
  # flattened-to-flattened, then map back to the raw form for reporting.
  u_gids_flat <- flatten_gid(u_gids)
  missing_flat <- setdiff(u_gids_flat, covered)
  missing_gids <- u_gids[u_gids_flat %in% missing_flat]

  cat(sprintf(
    paste0(
      "---- climate_stats coverage: %d of %d requested gids present, ",
      "%d absent\n"
    ),
    length(u_gids) - length(missing_gids),
    length(u_gids),
    length(missing_gids)
  ))

  if (length(missing_gids)) {
    # A gid with no raw CHIRPS/ERA5 file can never enter climate_stats.
    # Name them now rather than letting make_size_data()'s inner_join eat them.
    #
    # The CSV branch can only test for a *cached* stats file, which conflates
    # "never computed" with "no raw input". The parquet branch tests the raw
    # input itself, so `has_raw` below cleanly separates absent-upstream (GAB,
    # the one African country the release does not publish) from not-yet-built.
    if (identical(opt$climate_source, "parquet")) {
      has_raw <- !vapply(
        missing_gids,
        function(g) {
          is.null(get_clim_parquet_gid(g, clim_dir = opt$climate_parquet_dir))
        },
        logical(1),
        USE.NAMES = FALSE
      )
    } else {
      has_raw <- file.exists(
        make_clim_stats_file(var = "precip", gid = flatten_gid(missing_gids))
      )
    }
    no_raw <- missing_gids[!has_raw]

    if (length(no_raw)) {
      cat(sprintf(
        paste0(
          "---- WARNING: %d of %d absent gid(s) have no raw %s input and ",
          "will be dropped\n     e.g. %s\n"
        ),
        length(no_raw),
        length(missing_gids),
        opt$climate_source,
        paste(utils::head(no_raw, 5), collapse = ", ")
      ))
    }

    utils::write.csv(
      data.frame(
        gid = missing_gids,
        in_centroid_universe = missing_gids %in%
          matched_gadm$gid[!is.na(matched_gadm$gid)],
        has_raw = has_raw
      ),
      "generated_data/qc_climate_gids_missing.csv",
      row.names = FALSE
    )
    cat("---- wrote generated_data/qc_climate_gids_missing.csv\n")
  }

  if (isTRUE(opt$cost_only)) {
    cat("---- --cost_only: stopping before any extraction\n")
    quit(save = "no", status = 0)
  }
}

# assert_unique_clim_rows() lives in utils/clim_db.R.

# 5c. Per-gid cache manifest ---------------------------------------------------
# The parquet caches (clim_stats_all_parquet_<gid>.rds and the extremes built
# from them) are keyed by gid alone, but their values depend on the weighting
# and the parquet corpus too. Record those once for the cache directory and
# refuse to read caches built under other settings.
if (identical(opt$climate_source, "parquet")) {
  cache_manifest <- file.path(
    dirname(make_clim_stats_file("x", "y")),
    "_cache_manifest.yml"
  )
  want_cache <- list(
    weighting = opt$weighting,
    parquet_dir = sub("/+$", "", opt$climate_parquet_dir),
    cum7_definition = "corrected_7day_trailing"
  )
  if (!file.exists(cache_manifest)) {
    # First run under the manifest: the caches on disk were built by the build
    # the provenance file records, if that was a parquet build.
    prov0 <- if (file.exists(DEFAULT_PROV_PATH)) {
      yaml::read_yaml(DEFAULT_PROV_PATH)
    } else {
      list()
    }
    have_cache <- if (identical(prov0$climate_source, "parquet")) {
      list(
        weighting = prov0$weighting,
        parquet_dir = sub("/+$", "", prov0$parquet_dir),
        cum7_definition = prov0$cum7_definition
      )
    } else {
      want_cache
    }
    dir.create(dirname(cache_manifest), recursive = TRUE, showWarnings = FALSE)
    yaml::write_yaml(have_cache, cache_manifest)
    cat("---- wrote", cache_manifest, "\n")
  }
  have_cache <- yaml::read_yaml(cache_manifest)
  if (!identical(have_cache[names(want_cache)], want_cache)) {
    stop(
      "The per-gid parquet caches in ",
      dirname(cache_manifest),
      " were built with ",
      paste(
        sprintf("%s=%s", names(have_cache), unlist(have_cache)),
        collapse = ", "
      ),
      ", this run asks for ",
      paste(
        sprintf("%s=%s", names(want_cache), unlist(want_cache)),
        collapse = ", "
      ),
      ".\n  Reusing them would mix the two silently. Move the cache directory aside ",
      "(and its _cache_manifest.yml) to rebuild, or rerun with the recorded settings."
    )
  }
}

compile_stats_gid <- function(g) {
  if (identical(opt$climate_source, "parquet")) {
    compile_all_stats_parquet(
      g,
      weighting = opt$weighting,
      redo = FALSE,
      clim_dir = opt$climate_parquet_dir
    )
  } else {
    compile_all_stats(g, redo = FALSE)
  }
}

compile_extremes_gid <- function(g) {
  if (identical(opt$climate_source, "parquet")) {
    compile_all_extremes_parquet(
      g,
      redo = FALSE,
      baseline = opt$baseline,
      weighting = opt$weighting,
      clim_dir = opt$climate_parquet_dir
    )
  } else {
    compile_all_extremes(g, redo = FALSE, baseline = opt$baseline)
  }
}

if (identical(opt$stage, "extract")) {
  # Fill the caches only; nothing is collected or written to the databases.
  future_walk(u_gids, function(g) {
    compile_stats_gid(g)
    compile_extremes_gid(g)
    NULL
  })
  cat(sprintf(
    "---- --stage extract: caches present for %d gid(s)\n",
    length(u_gids)
  ))
  quit(save = "no", status = 0)
}

#' assert_db_provenance
#' Refuse to overwrite the database the provenance file describes with a build
#' of other settings (source, gid universe, weighting): the deployed tables
#' would change era without anything downstream noticing. A staged database
#' (another path) is not guarded.
assert_db_provenance <- function(db_path, prov_key) {
  prov_path <- PROV_PATH
  if (isTRUE(opt$force) || !file.exists(db_path) || !file.exists(prov_path)) {
    return(invisible())
  }
  prov <- yaml::read_yaml(prov_path)
  if (
    is.null(prov[[prov_key]]) ||
      normalizePath(prov[[prov_key]], mustWork = FALSE) !=
        normalizePath(db_path, mustWork = FALSE)
  ) {
    return(invisible())
  }
  want <- list(climate_source = opt$climate_source, gid_source = opt$gid_source)
  if (identical(opt$climate_source, "parquet")) {
    want$weighting <- opt$weighting
  }
  bad <- names(want)[
    !vapply(
      names(want),
      function(k) identical(prov[[k]], want[[k]]),
      logical(1)
    )
  ]
  if (length(bad)) {
    stop(
      db_path,
      " was built with ",
      paste(
        sprintf(
          "%s=%s",
          bad,
          vapply(bad, function(k) as.character(prov[[k]] %||% "NA"), "")
        ),
        collapse = ", "
      ),
      " (",
      prov_path,
      "); this run asks for ",
      paste(sprintf("%s=%s", bad, unlist(want[bad])), collapse = ", "),
      ".\n  Pass --force TRUE to overwrite it, or --stats_db / --extremes_db to stage elsewhere."
    )
  }
  invisible()
}

if (opt$redo_stats) {
  # Get all the climate data
  clim_data <- future_map_dfr(u_gids, compile_stats_gid)

  # # Save for later use
  # saveRDS(clim_data, file = "generated_data/all_climate_stats.rds")

  assert_unique_clim_rows(clim_data, c("gid", "month"), "climate_stats")

  # Write to database for later use
  assert_db_provenance(opt$stats_db, "stats_db")
  stats_db <- get_stats_db(opt$stats_db)
  dbWriteTable(stats_db, "climate_stats", clim_data, overwrite = TRUE)
  dbExecute(stats_db, 'CREATE INDEX idx_gid ON climate_stats(gid)')
  # overwrite = TRUE drops the old table but SQLite keeps its pages: without
  # VACUUM the file only grows (8.6 GB holding ~116 MB of live rows).
  dbExecute(stats_db, "VACUUM")
  dbDisconnect(stats_db)
} else {
  # NOTE the asymmetry with the extremes block below, which always writes: an
  # invocation without `-s TRUE` computes the gid universe and the coverage
  # report, prints this line, and leaves climate_stats untouched. A parquet
  # rebuild MUST pass -s TRUE.
  cat("---- Skipping stats computation \n")
}

# 6. Compute statistics --------------------------------------------------------

# Compute all extremes
clim_extrem_data <- future_map_dfr(u_gids, compile_extremes_gid)

# # Save for later use
# saveRDS(clim_extrem_data, file = "generated_data/all_climate_extremes.rds")

assert_unique_clim_rows(
  clim_extrem_data,
  c("gid", "month", "name"),
  glue("clim_extremes_{opt$baseline}")
)

assert_db_provenance(opt$extremes_db, "extremes_db")
extrm_db <- get_extrm_db(opt$extremes_db)
dbWriteTable(
  extrm_db,
  glue("clim_extremes_{opt$baseline}"),
  clim_extrem_data,
  overwrite = TRUE
)
dbExecute(
  extrm_db,
  glue(
    'CREATE INDEX idx_gid_{opt$baseline} ON clim_extremes_{opt$baseline}(gid)'
  )
)
dbExecute(
  extrm_db,
  glue(
    'CREATE INDEX idx_name_{opt$baseline} ON clim_extremes_{opt$baseline}(name)'
  )
)
dbExecute(
  extrm_db,
  glue(
    'CREATE INDEX idx_gid_name_{opt$baseline} ON clim_extremes_{opt$baseline}(gid, name)'
  )
)
dbExecute(extrm_db, "VACUUM")
dbDisconnect(extrm_db)

# 7. Provenance ----------------------------------------------------------------
# Which raw source, which weighting and which cum7 definition produced the
# databases is not recoverable from the databases themselves, and the two eras
# are NOT interchangeable (see cum7_definition). Stamp it so a downstream script
# -- or a future reader of qc_climate_gids_missing.csv -- can tell which era it
# is looking at instead of silently mixing them.
{
  dl_manifest <- file.path(opt$climate_parquet_dir, "DOWNLOAD_MANIFEST.tsv")

  prov <- list(
    climate_source = opt$climate_source,
    baseline = opt$baseline,
    gid_source = opt$gid_source,
    redo_stats = opt$redo_stats,
    stats_db = opt$stats_db,
    extremes_db = opt$extremes_db,
    n_gids_requested = length(u_gids),
    n_gids_extremes = length(unique(clim_extrem_data$gid)),
    n_rows_extremes = nrow(clim_extrem_data),
    month_min = format(min(clim_extrem_data$month, na.rm = TRUE)),
    month_max = format(max(clim_extrem_data$month, na.rm = TRUE)),
    git_sha = tryCatch(
      system("git rev-parse --short HEAD", intern = TRUE),
      error = function(e) NA_character_
    ),
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )

  if (identical(opt$climate_source, "parquet")) {
    prov$weighting <- opt$weighting
    prov$parquet_dir <- opt$climate_parquet_dir
    # The CSV era's mean_cum7/max_cum7 summed value + lag1..lag6 + lag6 AGAIN
    # (clim_compile.R) -- eight terms, oldest day double counted. The parquet
    # path uses a true 7-day trailing window, so *_cum7 values are NOT
    # comparable across eras. This marker is the only record of that.
    prov$cum7_definition <- "corrected_7day_trailing"
    prov$download_manifest <- if (file.exists(dl_manifest)) {
      readLines(dl_manifest)
    } else {
      "MISSING"
    }
  } else {
    prov$cum7_definition <- "legacy_8term_lag6_doubled"
  }

  if (opt$redo_stats) {
    prov$n_gids_stats <- length(unique(clim_data$gid))
    prov$n_rows_stats <- nrow(clim_data)
    prov$stats_columns <- names(clim_data)
  }

  yaml::write_yaml(prov, PROV_PATH)
  cat("---- wrote", PROV_PATH, "\n")
}
