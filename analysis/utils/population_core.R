# population_core.R — kernel of the LP x gid population crosswalk: WorldPop
#   mass inside every location-period x GADM-unit intersection, per year
#
# Provides:
#   run_population_task(), run_population_country_task()   one extraction
#   load_gadm_units(), load_lp_geometries()   disk loaders of the geometries
#   derive_crosswalk_shares(), compute_gid_cover(), derive_lp_tiling_level()
#   reconcile_crosswalk(), detect_geom_dup_lps()   QC of the crosswalk
#   drop_sliver_pieces()          drop negligible LP x gid pieces
#   write_wpp_national_totals(), load_wpp_national_totals(),
#   wpp_national_total()          WPP2024 national totals (CSV cache)
#   pop_gid_cache_key(), read_gid_pop_cache(), write_gid_pop_cache()
#   assert_boundary_cached()      the UN-adjustment precondition
#
# Details:
#   Runnable entry points: 02b_preprocess_population_run.R (one (country,
#   year) task) and 02c_preprocess_population_compile.R (bind + derive + QC).
#   `extract_fn` is injectable so tests run the algebra on a stub.
#
# See: the "What this builds" comment above POP_ROLES for the algebra

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# H[l, g, t] = UN-adjusted WorldPop mass in the intersection of location
# period geometry G_l and GADM polygon A_g, in year t. Downstream is algebra
# on H:
#   P_l       = sum_g H[l,g]    LP population (denominator of an observation)
#   Q_g       = extraction over A_g, gid population (forecast unit)
#   w[l->g]   = H[l,g] / P_l    share of LP l's people in gid g (sums to 1)
#   c[g<-l]   = H[l,g] / Q_g    share of gid g covered by LP l
#   gid_cover = max_g w[l->g]   1 for a clean LP, ~1/K for a K-way composite
#   n_gid_eff = 1 / sum_g w^2   inverse-Simpson composite index
# P_direct is extracted independently over the whole LP geometry, so
# P_l == P_direct is a real reconciliation check (the main geometry QC).
# An LP whose mass spreads over more than one gid is a composite by
# construction: no string parsing, no upstream `composite_loc`.
# All geometry roles go into one extraction call (one raster load per task);
# the per-(country, year) adjustment factor cancels in w[l->g].
#
# Geometry roles concatenated into the single extraction call.
POP_ROLES <- c("gadm1", "gadm2", "lp_whole", "lp_piece")

# WorldPop constrained coverage. Outside this range the compile step rescales by
# the WPP2024 national ratio with the shares frozen (see pop_tier).
POP_RASTER_YEAR_MIN <- 2015L
POP_RASTER_YEAR_MAX <- 2030L


# ---------------------------------------------------------------------------
# WPP2024 national totals
# ---------------------------------------------------------------------------
# Cached to CSV so that nothing downstream of the extraction needs
# OutbreakExtractR installed. Two consumers:
#   * natl_ratio in the task manifest — the only per-task signal that the UN
#     adjustment actually landed (adj_factor itself is unrecoverable; see the
#     DEVIATION note in the header);
#   * the pop_tier rescale for years outside WorldPop's 2015-2030 coverage,
#     where the gid shares are frozen and only the national level moves.

#' Write generated_data/wpp2024_national_totals.csv from the bundled table.
#' Requires OutbreakExtractR; run once on a login node, never on a worker.
write_wpp_national_totals <- function(path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    return(invisible(path))
  }
  if (!requireNamespace("OutbreakExtractR", quietly = TRUE)) {
    stop(
      "write_wpp_national_totals(): OutbreakExtractR is not installed. ",
      "Copy an existing ",
      basename(path),
      " into place instead."
    )
  }
  e <- new.env(parent = emptyenv())
  utils::data("WPP2024", package = "OutbreakExtractR", envir = e)
  w <- e$WPP2024
  out <- data.frame(
    iso3 = as.character(w$ISO3_code),
    year = as.integer(w$Time),
    pop_total = as.numeric(w$PopTotal) * 1e3, # the bundled table is thousands
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$iso3) & nzchar(out$iso3) & !is.na(out$year), ]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
  invisible(path)
}

load_wpp_national_totals <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' National total for (iso3, year); NA when the table or the row is absent.
wpp_national_total <- function(wpp, iso3, year) {
  if (is.null(wpp) || !nrow(wpp)) {
    return(NA_real_)
  }
  hit <- wpp$iso3 == toupper(iso3) & wpp$year == as.integer(year)
  if (!any(hit)) {
    return(NA_real_)
  }
  as.numeric(wpp$pop_total[which(hit)[1]])
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

#' Sanitize an sfc for extraction.
#'
#' estimate_pop_for_geometries() SILENTLY DROPS non-dimension-2 geometries
#' (OutbreakExtractR add_population.R), leaving NA for that row — which for a
#' crosswalk piece is invisible mass loss, not a visible failure.
#' st_intersection() of two polygon layers routinely yields GEOMETRYCOLLECTIONs
#' containing stray points and lines along shared edges, so this must run before
#' the extraction, not after.
#'
#' @param g sfc
#' @return list(geom = sfc, repair = character vector, empty = logical vector)
sanitize_polygons <- function(g) {
  n <- length(g)
  repair <- rep("none", n)

  types <- as.character(sf::st_geometry_type(g))

  # A GEOMETRYCOLLECTION may still hold polygonal mass alongside the debris.
  gc <- which(types == "GEOMETRYCOLLECTION")
  for (i in gc) {
    ex <- try(
      sf::st_collection_extract(g[i], "POLYGON", warn = FALSE),
      silent = TRUE
    )
    g[i] <- if (inherits(ex, "try-error") || !length(ex)) {
      sf::st_polygon()
    } else {
      sf::st_union(ex)[[1]]
    }
    repair[i] <- "collection_extract"
  }

  # Anything else non-polygonal (POINT / LINESTRING / ... from an edge-only or
  # vertex-only overlap) carries no area, so it carries no people. Emptying it
  # here is what stops estimate_pop_for_geometries() dropping it silently and
  # leaving an NA that would break sum_g H == P_direct.
  drop <- which(!types %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))
  for (i in drop) {
    g[i] <- sf::st_polygon()
  }
  repair[drop] <- "non_polygonal"

  empty <- sf::st_is_empty(g)
  repair[empty & repair == "none"] <- "empty"
  repair[empty & repair == "collection_extract"] <- "collection_extract_empty"

  # Cast the non-empty ones to a single type so exact_extract sees a uniform sfc.
  ok <- !empty
  if (any(ok)) {
    g[ok] <- sf::st_cast(g[ok], "MULTIPOLYGON", warn = FALSE)
  }

  list(geom = g, repair = repair, empty = empty)
}

#' Number of polygon parts in each geometry (a composite LP has many).
count_geom_parts <- function(g) {
  vapply(
    sf::st_geometry(g),
    function(x) {
      if (inherits(x, "MULTIPOLYGON")) length(x) else 1L
    },
    integer(1)
  )
}

#' Planimetric-free area in km^2, matching the convention already used for
#' `area_per_1km2` in analysis/00_ingest_outbreakextractr.R.
geom_area_km2 <- function(g) {
  suppressWarnings(as.numeric(sf::st_area(sf::st_geometry(g))) / 1e6)
}

#' Load the GADM ADM1 + ADM2 units of one country as a single long sf.
#'
#' Reads the same validated GPKG that analysis/02_preprocess_climate.R section 4
#' uses, so the gids produced here are exactly the gids the climate panel and
#' the simulator frame are keyed by.
#'
#' Each layer is read with an SQL filter on GID_0, so a task reads one country's
#' units instead of the whole continent; `ORDER BY rowid` keeps the order of a
#' full-layer read.
#'
#' @return sf with columns gid, gadm_level ("ADM1"|"ADM2"), gid_iso3
load_gadm_units <- function(gadm_path, iso3, levels = c("ADM1", "ADM2")) {
  stopifnot(file.exists(gadm_path))
  iso3 <- toupper(iso3)
  stopifnot(grepl("^[A-Z]{3}$", iso3))
  out <- list()
  read_layer <- function(layer) {
    sf::st_read(
      gadm_path,
      quiet = TRUE,
      query = sprintf(
        "SELECT * FROM \"%s\" WHERE GID_0 = '%s' ORDER BY rowid",
        layer,
        iso3
      )
    )
  }

  if ("ADM1" %in% levels) {
    a1 <- read_layer("ADM_1")
    if (nrow(a1)) {
      out$ADM1 <- sf::st_sf(
        gid = as.character(a1$GID_1),
        gadm_level = "ADM1",
        gid_iso3 = iso3,
        geometry = sf::st_geometry(a1)
      )
    }
  }
  if ("ADM2" %in% levels) {
    a2 <- read_layer("ADM_2")
    # Not every country has an ADM2 division in GADM 4.1 — absence is normal.
    if (nrow(a2)) {
      out$ADM2 <- sf::st_sf(
        gid = as.character(a2$GID_2),
        gadm_level = "ADM2",
        gid_iso3 = iso3,
        geometry = sf::st_geometry(a2)
      )
    }
  }

  if (!length(out)) {
    stop("load_gadm_units(): no GADM units found for ", iso3, " in ", gadm_path)
  }
  do.call(rbind, unname(out))
}

#' Load LP geometries as written by analysis/00_ingest_outbreakextractr.R
#' section 7.
#' @return sf with column lctn_pr
load_lp_geometries <- function(lp_path) {
  stopifnot(file.exists(lp_path))
  shp <- if (grepl("\\.rds$", lp_path, ignore.case = TRUE)) {
    readRDS(lp_path)
  } else {
    sf::st_read(lp_path, quiet = TRUE)
  }
  if (!"lctn_pr" %in% names(shp)) {
    stop("load_lp_geometries(): expected a `lctn_pr` column in ", lp_path)
  }
  sf::st_sf(
    lctn_pr = as.character(shp$lctn_pr),
    geometry = sf::st_geometry(shp)
  )
}


# ---------------------------------------------------------------------------
# Crosswalk construction
# ---------------------------------------------------------------------------

#' Select the LPs belonging to a (country, year) task.
#'
#' The membership rule is "LPs whose geometry INTERSECTS country X", not "LPs
#' whose `country` column is X". A cross-border LP is therefore visited by both
#' countries' tasks, each piece is extracted from its own country's raster with
#' its own adjustment factor, and the compile step sums P_direct across tasks.
#' Selecting on an attribute column instead would extract the foreign part of
#' the LP from a raster that has no data there and silently under-count it.
#'
#' Membership is by POSITIVE-AREA overlap, not st_intersects. Neighbouring
#' polygons share an edge, so an LP that merely abuts the border "intersects" it
#' and would otherwise enter the neighbour's task — contributing a spurious
#' P_direct row and a spurious cross_border flag, both of which the compile step
#' would then sum in. Area is also what decides cross_border: an LP is
#' cross-border when part of its area lies outside this country, which is
#' exactly the condition under which its pieces arrive from two tasks.
#'
#' @return list(lp = sf subset, cross_border = logical aligned to that subset)
select_lps_for_country <- function(lp_sf, country_geom, area_tol = 1e-6) {
  if (!nrow(lp_sf)) {
    return(list(lp = lp_sf, cross_border = logical(0)))
  }

  cand <- lengths(sf::st_intersects(lp_sf, country_geom)) > 0L
  lp <- lp_sf[cand, , drop = FALSE]
  if (!nrow(lp)) {
    return(list(lp = lp, cross_border = logical(0)))
  }

  # st_intersection() silently drops empty results, so carry an index through it
  # rather than assuming the output is row-aligned with the input.
  lp$.lp_idx <- seq_len(nrow(lp))
  cg <- sf::st_sf(
    geometry = if (inherits(country_geom, "sfc")) {
      country_geom
    } else {
      sf::st_sfc(country_geom, crs = sf::st_crs(lp))
    }
  )
  inter <- suppressWarnings(sf::st_intersection(lp, cg))
  a_in <- numeric(nrow(lp))
  if (nrow(inter)) {
    agg <- tapply(
      suppressWarnings(as.numeric(sf::st_area(inter))),
      inter$.lp_idx,
      sum
    )
    a_in[as.integer(names(agg))] <- as.numeric(agg)
  }
  a_tot <- suppressWarnings(as.numeric(sf::st_area(lp)))

  keep <- a_in > area_tol * pmax(a_tot, .Machine$double.eps)
  lp <- lp[keep, , drop = FALSE]
  lp$.lp_idx <- NULL

  list(lp = lp, cross_border = a_in[keep] < (1 - area_tol) * a_tot[keep])
}

#' LP x gid intersection pieces, sanitized and ready for extraction.
#'
#' Vector intersection, deliberately not rasterization: exact_extract already
#' does exact partial-cell weighting, so handing it intersection polygons is
#' exact and free, whereas rasterizing GADM onto the 100 m grid would cost ~1e9
#' cells and reintroduce the very discretization bias this removes.
#'
#' @return sf with lctn_pr, gid, gadm_level, gid_iso3, n_parts, geom_repair
build_crosswalk_pieces <- function(lp_sf, gadm_sf) {
  if (!nrow(lp_sf) || !nrow(gadm_sf)) {
    return(sf::st_sf(
      lctn_pr = character(0),
      gid = character(0),
      gadm_level = character(0),
      gid_iso3 = character(0),
      n_parts = integer(0),
      geom_repair = character(0),
      geometry = sf::st_sfc(crs = sf::st_crs(gadm_sf))
    ))
  }

  # st_intersection on two sf objects returns the cross product of overlapping
  # pairs with both attribute sets; GEOS prunes candidates via an STRtree, so
  # ~60 LPs x 774 ADM2 polygons is seconds, not minutes.
  pieces <- suppressWarnings(sf::st_intersection(lp_sf, gadm_sf))
  if (!nrow(pieces)) {
    return(sf::st_sf(
      lctn_pr = character(0),
      gid = character(0),
      gadm_level = character(0),
      gid_iso3 = character(0),
      n_parts = integer(0),
      geom_repair = character(0),
      geometry = sf::st_sfc(crs = sf::st_crs(gadm_sf))
    ))
  }

  san <- sanitize_polygons(sf::st_geometry(pieces))
  out <- sf::st_sf(
    lctn_pr = as.character(pieces$lctn_pr),
    gid = as.character(pieces$gid),
    gadm_level = as.character(pieces$gadm_level),
    gid_iso3 = as.character(pieces$gid_iso3),
    geom_repair = san$repair,
    geometry = san$geom
  )
  # Empty pieces carry no mass and would only produce NA rows downstream.
  out <- out[!san$empty, , drop = FALSE]
  if (nrow(out)) {
    out$n_parts <- count_geom_parts(out)
  } else {
    out$n_parts <- integer(0)
  }
  out
}


# ---------------------------------------------------------------------------
# The task worker
# ---------------------------------------------------------------------------

# WHY THE MANIFEST HAS NO `adj_factor`
# ------------------------------------
# `adj_factor` cannot be recovered:
# estimate_pop_for_geometries() computes it internally and returns only the
# adjusted numbers, and band_adj_factor() (OutbreakExtractR
# R/get_country_boundary.R) is identity on [0.5, 2.0] and 1.0 outside,
# so the adjusted national total alone does not identify it. Two substitutes,
# both strictly better for the failure mode that matters:
#   1. a HARD precondition that the boundary cache file exists — this is the
#      exact detector for get_country_boundary() returning NULL, which SILENTLY
#      SKIPS the UN adjustment (OutbreakExtractR get_country_boundary.R) rather
#      than failing;
#   2. `natl_ratio` = sum(Q_g over ADM1) / WPP2024 national total, which is free
#      because the ADM1 polygons tile the country and are extracted anyway.

#' Assert the national-boundary cache entry exists.
#'
#' get_country_boundary() (OutbreakExtractR R/get_country_boundary.R)
#' falls back to NULL when the disk cache misses and rgeoboundaries fails, and
#' estimate_pop_for_geometries() then treats NULL as "skip the UN adjustment"
#' (OutbreakExtractR add_population.R) rather than as an error. Every number in
#' the run would look plausible and be systematically wrong by the adjustment
#' factor. Fail loudly instead: the boundary cache must be built (on a node with
#' internet) before any array is submitted.
assert_boundary_cached <- function(iso3, boundary_cache_dir) {
  f <- file.path(boundary_cache_dir, paste0(toupper(iso3), "_adm0.geojson"))
  if (!file.exists(f)) {
    stop(
      "Boundary cache miss for ",
      iso3,
      " (expected ",
      f,
      ").\n",
      "  get_country_boundary() would return NULL and the UN adjustment ",
      "would be SILENTLY SKIPPED.\n",
      "  Build the boundary cache (",
      boundary_cache_dir,
      ") first, or pass ",
      "require_boundary = FALSE if you accept unadjusted population."
    )
  }
  invisible(f)
}

# TESTABILITY: `extract_fn` is a defaulted argument, so tests substitute a
# deterministic area-proportional stub and exercise the whole crosswalk algebra
# with no raster, no network and no data/ access. run_population_task()
# likewise takes `gadm_sf` / `lp_sf` as sf objects; the disk loaders are
# separate functions.
#' Run one (country, year) extraction task.
#'
#' @param iso3        ISO3 country code; owns the WorldPop raster for this task
#' @param year        extraction year (clamped to [2015, 2030] by the extractor)
#' @param gadm_sf     sf with gid, gadm_level, gid_iso3 (see load_gadm_units)
#' @param lp_sf       sf with lctn_pr (see load_lp_geometries); subset internally
#' @param extract_fn  the extractor. Injected so tests can substitute a
#'                    deterministic stub — see test_population_crosswalk.R.
#'                    Contract: (geom_sf, country_iso3, year, raster_dir,
#'                    boundary_cache_dir) -> numeric of length nrow(geom_sf).
#' @param wpp_total   WPP2024 national total for (iso3, year), for natl_ratio.
#' @return list(crosswalk, lp_pop, gid_pop, manifest)
run_population_task <- function(
  iso3,
  year,
  gadm_sf,
  lp_sf,
  extract_fn = NULL,
  raster_dir = "data/worldpop",
  boundary_cache_dir = "data/country_boundaries",
  sliver_mass_frac = 1e-4,
  sliver_mass_min = 1.0,
  wpp_total = NA_real_,
  require_boundary = TRUE,
  quiet = FALSE,
  progress_tag = ""
) {
  iso3 <- toupper(iso3)
  year <- as.integer(year)
  t0 <- Sys.time()

  extract_fn <- .pop_extract_fn(extract_fn, "run_population_task")
  if (isTRUE(require_boundary)) {
    assert_boundary_cached(iso3, boundary_cache_dir)
  }

  say <- .pop_say(quiet, progress_tag, t0)
  geo <- .pop_country_geometry(iso3, gadm_sf, lp_sf, say, "run_population_task")
  by_role <- .pop_extract_roles(
    geo,
    iso3,
    year,
    extract_fn,
    raster_dir,
    boundary_cache_dir,
    say,
    "run_population_task"
  )
  .pop_assemble_task(
    geo,
    by_role,
    iso3,
    year,
    gid_pop = NULL,
    sliver_mass_frac = sliver_mass_frac,
    sliver_mass_min = sliver_mass_min,
    wpp_total = wpp_total,
    raster_dir = raster_dir,
    boundary_cache_dir = boundary_cache_dir,
    t0 = t0
  )
}

#' Run one country over every year: geometry once, one extraction per year.
#'
#' The per-(country, year) run_population_task() rebuilds the country footprint,
#' the LP selection and the LP x gid pieces for every year, although none of
#' them depends on the year. Here they are built once and each year reuses
#' them; per year the extraction call and the assembly are run_population_task's,
#' so a year's result is identical to the per-year task's.
#'
#' The GADM-unit populations Q_g (roles gadm1 / gadm2) do not depend on the
#' outbreak vintage. `gid_pop_cache` hands them in per year (a named list,
#' names = years, of gid_pop tibbles as returned here); a year found there
#' extracts the LP roles only, and with no LP in the country no extraction runs
#' at all. The UN adjustment factor is computed from the country boundary and
#' the year alone (OutbreakExtractR R/add_population.R), so dropping
#' the GADM roles from the call does not move the LP values.
#'
#' @param years         extraction years
#' @param wpp_totals    WPP2024 national totals, one per year (NA allowed)
#' @param gid_pop_cache named list (names = years) of cached gid_pop; NULL/empty
#'                      extracts every role for every year
#' @param workers       years extracted at once (forked processes,
#'                      future::multicore, at most one per year); 1 = one after
#'                      another. Each year is computed exactly as with 1 worker
#'                      (no RNG), so the result does not depend on it; memory
#'                      grows with it (one raster per worker).
#' @return named list (names = years) of list(crosswalk, lp_pop, gid_pop,
#'   manifest); manifest$gid_pop_cached says whether Q_g came from the cache.
run_population_country_task <- function(
  iso3,
  years,
  gadm_sf,
  lp_sf,
  extract_fn = NULL,
  raster_dir = "data/worldpop",
  boundary_cache_dir = "data/country_boundaries",
  sliver_mass_frac = 1e-4,
  sliver_mass_min = 1.0,
  wpp_totals = rep(NA_real_, length(years)),
  gid_pop_cache = list(),
  require_boundary = TRUE,
  quiet = FALSE,
  progress_tag = "",
  workers = 1L
) {
  iso3 <- toupper(iso3)
  years <- as.integer(years)
  if (length(wpp_totals) != length(years)) {
    stop(
      "run_population_country_task(): wpp_totals must have one value per year."
    )
  }
  t0 <- Sys.time()

  extract_fn <- .pop_extract_fn(extract_fn, "run_population_country_task")
  if (isTRUE(require_boundary)) {
    assert_boundary_cached(iso3, boundary_cache_dir)
  }

  say <- .pop_say(quiet, progress_tag, t0)
  geo <- .pop_country_geometry(
    iso3,
    gadm_sf,
    lp_sf,
    say,
    "run_population_country_task"
  )

  one_year <- function(i) {
    yr <- years[i]
    t_yr <- Sys.time()
    say_y <- .pop_say(quiet, paste0(progress_tag, yr, ": "), t_yr)
    cached <- gid_pop_cache[[as.character(yr)]]
    by_role <- NULL
    if (!is.null(cached)) {
      # LP roles alone can be all NA legitimately (only unpopulated slivers:
      # zero is returned as NA), which the raster-failure guard cannot tell
      # from a failed load. Then redo the year with every role, as a cold run.
      by_role <- tryCatch(
        .pop_extract_roles(
          geo,
          iso3,
          yr,
          extract_fn,
          raster_dir,
          boundary_cache_dir,
          say_y,
          "run_population_country_task",
          roles = c("lp_whole", "lp_piece")
        ),
        error = function(e) {
          if (!grepl("all-NA", conditionMessage(e), fixed = TRUE)) {
            stop(e)
          }
          say_y("LP roles all NA with Q_g cached: extracting every role")
          NULL
        }
      )
      if (is.null(by_role)) cached <- NULL
    }
    if (is.null(by_role)) {
      by_role <- .pop_extract_roles(
        geo,
        iso3,
        yr,
        extract_fn,
        raster_dir,
        boundary_cache_dir,
        say_y,
        "run_population_country_task"
      )
    }
    res <- .pop_assemble_task(
      geo,
      by_role,
      iso3,
      yr,
      gid_pop = cached,
      sliver_mass_frac = sliver_mass_frac,
      sliver_mass_min = sliver_mass_min,
      wpp_total = wpp_totals[i],
      raster_dir = raster_dir,
      boundary_cache_dir = boundary_cache_dir,
      t0 = t_yr
    )
    res$manifest$gid_pop_cached <- !is.null(cached)
    res
  }

  workers <- max(1L, min(as.integer(workers), length(years)))
  out <- if (workers > 1L) {
    old_plan <- future::plan(future::multicore, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
    say("extracting ", length(years), " year(s) on ", workers, " worker(s)")
    furrr::future_map(
      seq_along(years),
      one_year,
      .options = furrr::furrr_options(seed = NULL)
    )
  } else {
    lapply(seq_along(years), one_year)
  }
  names(out) <- as.character(years)
  out
}

# -- shared steps of the two task runners ------------------------------------

.pop_extract_fn <- function(extract_fn, caller) {
  if (!is.null(extract_fn)) {
    return(extract_fn)
  }
  if (!requireNamespace("OutbreakExtractR", quietly = TRUE)) {
    stop(
      caller,
      "(): package 'OutbreakExtractR' is not installed ",
      "and no `extract_fn` was supplied."
    )
  }
  OutbreakExtractR::estimate_pop_for_geometries
}

# message() writes to stderr and is never block-buffered, unlike cat() to a
# stdout redirected into a SLURM log file — cat() here could sit unflushed
# for the entire duration of a slow step, making a working task look stuck.
.pop_say <- function(quiet, progress_tag, t0) {
  function(...) {
    if (!quiet) {
      message(
        progress_tag,
        ...,
        " (",
        sprintf(
          "%.1fs",
          as.numeric(
            difftime(Sys.time(), t0, units = "secs")
          )
        ),
        ")"
      )
    }
  }
}

#' The year-independent geometry of one country: its GADM units, the LPs that
#' overlap it, their LP x gid pieces and the role-tagged geometries of the
#' extraction call.
.pop_country_geometry <- function(iso3, gadm_sf, lp_sf, say, caller) {
  gadm_sf <- gadm_sf[gadm_sf$gid_iso3 == iso3, , drop = FALSE]
  if (!nrow(gadm_sf)) {
    stop(caller, "(): no GADM units for ", iso3)
  }

  # The ADM1 units tile the country, so their union is the country footprint —
  # no separate ADM0 source, and pieces are guaranteed to tile by construction.
  adm1 <- gadm_sf[gadm_sf$gadm_level == "ADM1", , drop = FALSE]
  if (!nrow(adm1)) {
    stop(caller, "(): no ADM1 units for ", iso3)
  }
  country_geom <- sf::st_union(sf::st_geometry(adm1))
  say("country footprint built")

  sel <- select_lps_for_country(lp_sf, country_geom)
  lps <- sel$lp
  say(
    "LPs intersecting ",
    iso3,
    ": ",
    nrow(lps),
    " (cross-border: ",
    sum(sel$cross_border),
    ")"
  )

  pieces <- build_crosswalk_pieces(lps, gadm_sf)
  say("LP x gid pieces: ", nrow(pieces))

  lp_san <- if (nrow(lps)) sanitize_polygons(sf::st_geometry(lps)) else NULL

  bind_role <- function(role, g) {
    if (is.null(g) || !length(g)) {
      return(NULL)
    }
    sf::st_sf(role = rep(role, length(g)), geometry = g)
  }
  parts <- list(
    bind_role(
      "gadm1",
      sf::st_geometry(gadm_sf[gadm_sf$gadm_level == "ADM1", ])
    ),
    bind_role(
      "gadm2",
      sf::st_geometry(gadm_sf[gadm_sf$gadm_level == "ADM2", ])
    ),
    bind_role("lp_whole", if (!is.null(lp_san)) lp_san$geom else NULL),
    bind_role("lp_piece", if (nrow(pieces)) sf::st_geometry(pieces) else NULL)
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]

  list(
    gadm_sf = gadm_sf,
    sel = sel,
    lps = lps,
    pieces = pieces,
    geom_all = do.call(rbind, unname(parts))
  )
}

#' One extraction call over the requested roles; values split by role.
#' All four roles = run_population_task's single call. With no geometry left in
#' the requested roles no call is made.
.pop_extract_roles <- function(
  geo,
  iso3,
  year,
  extract_fn,
  raster_dir,
  boundary_cache_dir,
  say,
  caller,
  roles = POP_ROLES
) {
  geom <- geo$geom_all[geo$geom_all$role %in% roles, , drop = FALSE]
  if (!nrow(geom)) {
    say(
      "no geometry in roles ",
      paste(roles, collapse = "/"),
      ": no extraction"
    )
    return(split(numeric(0), factor(character(0), levels = POP_ROLES)))
  }

  say(
    "single extraction call: ",
    nrow(geom),
    " geometries, year ",
    year,
    " — raster load + exact_extract() next, this is usually the slow step"
  )
  extract_t0 <- Sys.time()
  extracted <- extract_fn(
    geom_sf = geom,
    country_iso3 = iso3,
    year = year,
    raster_dir = raster_dir,
    boundary_cache_dir = boundary_cache_dir
  )
  say(
    "extract_fn() returned in ",
    sprintf(
      "%.1fs",
      as.numeric(difftime(Sys.time(), extract_t0, units = "secs"))
    )
  )
  if (length(extracted) != nrow(geom)) {
    stop(
      caller,
      "(): extract_fn returned ",
      length(extracted),
      " values for ",
      nrow(geom),
      " geometries."
    )
  }
  # An all-NA return means the raster never loaded (the download-failure branch
  # in OutbreakExtractR add_population.R only `message()`s and moves on).
  # Coalescing piece-level NA to 0 below is only sound once that is ruled out
  # here.
  if (all(is.na(extracted))) {
    stop(
      caller,
      "(): extraction returned all-NA for ",
      iso3,
      " ",
      year,
      " — the WorldPop raster almost certainly failed to load."
    )
  }

  split(extracted, factor(geom$role, levels = POP_ROLES))
}

#' Assemble a task's tables from its geometry and the per-role extraction.
#' `gid_pop`, when given, replaces the Q_g derived from the gadm1/gadm2 roles.
.pop_assemble_task <- function(
  geo,
  by_role,
  iso3,
  year,
  gid_pop = NULL,
  sliver_mass_frac,
  sliver_mass_min,
  wpp_total,
  raster_dir,
  boundary_cache_dir,
  t0
) {
  gadm_sf <- geo$gadm_sf
  lps <- geo$lps
  pieces <- geo$pieces

  # -- gid populations Q_g --------------------------------------------------
  if (is.null(gid_pop)) {
    gid_pop <- dplyr::bind_rows(
      if (length(by_role$gadm1)) {
        tibble::tibble(
          gid = gadm_sf$gid[gadm_sf$gadm_level == "ADM1"],
          gadm_level = "ADM1",
          Q = as.numeric(by_role$gadm1)
        )
      },
      if (length(by_role$gadm2)) {
        tibble::tibble(
          gid = gadm_sf$gid[gadm_sf$gadm_level == "ADM2"],
          gadm_level = "ADM2",
          Q = as.numeric(by_role$gadm2)
        )
      }
    ) %>%
      dplyr::mutate(
        iso3 = iso3,
        year = year,
        area_km2 = geom_area_km2(gadm_sf)[match(
          paste(gid, gadm_level),
          paste(gadm_sf$gid, gadm_sf$gadm_level)
        )]
      ) %>%
      dplyr::select(gid, gadm_level, iso3, year, Q, area_km2)
  }

  # -- whole-LP populations P_direct (the INDEPENDENT reconciliation target) -
  lp_pop <- if (nrow(lps)) {
    tibble::tibble(
      lctn_pr = as.character(lps$lctn_pr),
      extract_iso3 = iso3,
      year = year,
      P_direct = as.numeric(by_role$lp_whole),
      lp_area_km2 = geom_area_km2(lps),
      cross_border = geo$sel$cross_border
    )
  } else {
    tibble::tibble(
      lctn_pr = character(0),
      extract_iso3 = character(0),
      year = integer(0),
      P_direct = numeric(0),
      lp_area_km2 = numeric(0),
      cross_border = logical(0)
    )
  }

  # -- crosswalk H[l,g] -----------------------------------------------------
  crosswalk <- if (nrow(pieces)) {
    tibble::tibble(
      lctn_pr = pieces$lctn_pr,
      gid = pieces$gid,
      gadm_level = pieces$gadm_level,
      gid_iso3 = pieces$gid_iso3,
      year = year,
      # Zero-to-NA is right for a whole LP ("zero is never a population",
      # OutbreakExtractR add_population.R) but wrong for a PIECE: a genuinely
      # empty sliver of an otherwise populated LP contributes 0, and leaving it
      # NA would break the additivity sum_g H == P_direct. Safe only because the
      # all-NA raster-failure case is hard-errored above.
      H = dplyr::coalesce(as.numeric(by_role$lp_piece), 0),
      piece_area_km2 = geom_area_km2(pieces),
      n_parts = pieces$n_parts,
      geom_repair = pieces$geom_repair
    )
  } else {
    tibble::tibble(
      lctn_pr = character(0),
      gid = character(0),
      gadm_level = character(0),
      gid_iso3 = character(0),
      year = integer(0),
      H = numeric(0),
      piece_area_km2 = numeric(0),
      n_parts = integer(0),
      geom_repair = character(0)
    )
  }

  crosswalk <- drop_sliver_pieces(
    crosswalk,
    sliver_mass_frac = sliver_mass_frac,
    sliver_mass_min = sliver_mass_min
  )

  # -- manifest -------------------------------------------------------------
  natl_from_adm1 <- sum(gid_pop$Q[gid_pop$gadm_level == "ADM1"], na.rm = TRUE)
  manifest <- list(
    iso3 = iso3,
    year = year,
    n_lp = nrow(lp_pop),
    n_lp_cross = sum(lp_pop$cross_border),
    n_gid = nrow(gid_pop),
    n_pieces = nrow(crosswalk),
    natl_from_adm1 = natl_from_adm1,
    wpp_total = as.numeric(wpp_total),
    # ~1 when the UN adjustment landed. Departures are the GADM-vs-geoBoundaries
    # national footprint difference (small) or a clamped/skipped adjustment.
    natl_ratio = if (is.finite(wpp_total) && wpp_total > 0) {
      natl_from_adm1 / wpp_total
    } else {
      NA_real_
    },
    raster_dir = raster_dir,
    boundary_cache = boundary_cache_dir,
    runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )

  list(
    crosswalk = crosswalk,
    lp_pop = lp_pop,
    gid_pop = gid_pop,
    manifest = manifest
  )
}

# ---------------------------------------------------------------------------
# Partition freshness and the Q_g cache
# ---------------------------------------------------------------------------

#' Is a finished partition still valid for this LP corpus?
#'
#' `_SUCCESS` alone says the writer finished, not that it read the current
#' outbreak vintage: a new vintage would silently reuse every old part. The
#' manifest records the LP corpus it was built from; a part is fresh only when
#' that matches. Manifests carry `lp_md5` since the per-country runner; older
#' ones are matched on path + mtime.
#'
#' @param part_dir      partition dir (holds _SUCCESS + manifest.yml)
#' @param lp_shapefiles the LP corpus this run reads
#' @param lp_md5        its md5 (pass it in to hash the file once per task)
#' @return TRUE when the partition can be skipped
pop_part_fresh <- function(
  part_dir,
  lp_shapefiles,
  lp_md5 = unname(tools::md5sum(path.expand(lp_shapefiles)))
) {
  man_path <- file.path(part_dir, "manifest.yml")
  if (!file.exists(file.path(part_dir, "_SUCCESS")) || !file.exists(man_path)) {
    return(FALSE)
  }
  man <- yaml::read_yaml(man_path)
  if (!is.null(man$lp_md5)) {
    return(identical(man$lp_md5, lp_md5))
  }
  identical(man$lp_shapefiles, lp_shapefiles) &&
    identical(man$lp_mtime, format(file.mtime(lp_shapefiles)))
}

#' The inputs Q_g depends on: the GADM units, the rasters, the boundary cache
#' and the extractor. Not the outbreak vintage.
pop_gid_cache_key <- function(
  iso3,
  year,
  gadm_path,
  raster_dir,
  boundary_cache
) {
  list(
    iso3 = toupper(iso3),
    year = as.integer(year),
    gadm_path = gadm_path,
    gadm_size = unname(file.size(gadm_path)),
    gadm_mtime = format(file.mtime(gadm_path)),
    raster_dir = raster_dir,
    boundary_cache = boundary_cache,
    oer_version = tryCatch(
      as.character(utils::packageVersion("OutbreakExtractR")),
      error = function(e) NA_character_
    )
  )
}

pop_gid_cache_path <- function(cache_dir, iso3, year) {
  file.path(
    cache_dir,
    sprintf("iso3=%s", toupper(iso3)),
    sprintf("year=%d.qs", as.integer(year))
  )
}

#' Cached Q_g for one (country, year), or NULL when absent or keyed on other inputs.
read_gid_pop_cache <- function(cache_dir, key) {
  f <- pop_gid_cache_path(cache_dir, key$iso3, key$year)
  if (!file.exists(f)) {
    return(NULL)
  }
  hit <- qs2::qs_read(f)
  if (identical(hit$key, key)) hit$gid_pop else NULL
}

write_gid_pop_cache <- function(cache_dir, key, gid_pop) {
  f <- pop_gid_cache_path(cache_dir, key$iso3, key$year)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp")
  qs2::qs_save(list(key = key, gid_pop = gid_pop), tmp)
  file.rename(tmp, f)
  invisible(f)
}

#' Drop sliver pieces by MASS, not area.
#'
#' GADM 4.1 and the LP geometries are different boundary vintages, so their
#' edges disagree by a few hundred metres and every LP picks up a fringe of
#' spurious neighbour pieces. Filtering on area would also delete small dense
#' urban pieces that carry real people; filtering on mass targets exactly the
#' pieces that do not matter.
#'
#' `w_dropped` (the fraction of P_l discarded) is retained per (lp, level) — if
#' it is large for some country, that is a vintage-mismatch finding to report,
#' not a threshold to tune away.
drop_sliver_pieces <- function(
  crosswalk,
  sliver_mass_frac = 1e-4,
  sliver_mass_min = 1.0
) {
  if (!nrow(crosswalk)) {
    return(dplyr::mutate(crosswalk, w_dropped = numeric(0)))
  }
  crosswalk %>%
    # ADM1 and ADM2 each tile the country independently, so an LP's pieces sum
    # to P_l SEPARATELY at each level. Grouping must include gadm_level or the
    # shares would be computed over a double-counted denominator.
    dplyr::group_by(lctn_pr, gadm_level) %>%
    dplyr::mutate(
      .P_all = sum(H, na.rm = TRUE),
      .thr = pmax(sliver_mass_min, sliver_mass_frac * .P_all),
      .keep_piece = H >= .thr
    ) %>%
    dplyr::mutate(
      w_dropped = if (all(.keep_piece)) {
        0
      } else {
        sum(H[!.keep_piece], na.rm = TRUE) /
          dplyr::first(pmax(.P_all, .Machine$double.eps))
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.keep_piece) %>%
    dplyr::select(-.P_all, -.thr, -.keep_piece)
}


# ---------------------------------------------------------------------------
# Derivation (used by the compile step and by the tests)
# ---------------------------------------------------------------------------

#' Derive P_l, Q_g, w, c, gid_cover, n_gid_eff from bound task outputs.
#'
#' Separated from run_population_task() because a cross-border LP's pieces arrive
#' from two different tasks: P_l is only correct once they are bound.
#'
#' @param crosswalk bound crosswalk rows (lctn_pr, gid, gadm_level, year, H, ...)
#' @param gid_pop   bound gid populations (gid, gadm_level, year, Q, ...)
#' @return crosswalk with P_l, Q, w, c, gid_cover, n_gid_eff, is_composite
derive_crosswalk_shares <- function(
  crosswalk,
  gid_pop,
  cover_threshold = 0.90
) {
  if (!nrow(crosswalk)) {
    return(dplyr::mutate(
      crosswalk,
      P_l = numeric(0),
      Q = numeric(0),
      w = numeric(0),
      c = numeric(0),
      gid_cover = numeric(0),
      n_gid_eff = numeric(0),
      is_composite = logical(0)
    ))
  }

  q <- gid_pop %>%
    dplyr::select(gid, gadm_level, year, Q) %>%
    dplyr::distinct(gid, gadm_level, year, .keep_all = TRUE)

  crosswalk %>%
    # Cross-border LPs contribute pieces from more than one task; summing H
    # first makes P_l whole before any share is taken.
    dplyr::group_by(lctn_pr, gadm_level, year) %>%
    dplyr::mutate(
      P_l = sum(H, na.rm = TRUE),
      w = dplyr::if_else(P_l > 0, H / P_l, NA_real_)
    ) %>%
    dplyr::mutate(
      gid_cover = suppressWarnings(max(w, na.rm = TRUE)),
      n_gid_eff = dplyr::if_else(
        sum(w^2, na.rm = TRUE) > 0,
        1 / sum(w^2, na.rm = TRUE),
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gid_cover = dplyr::if_else(is.finite(gid_cover), gid_cover, NA_real_),
      is_composite = !is.na(gid_cover) & gid_cover < cover_threshold
    ) %>%
    dplyr::left_join(q, by = c("gid", "gadm_level", "year")) %>%
    dplyr::mutate(c = dplyr::if_else(!is.na(Q) & Q > 0, H / Q, NA_real_))
}

#' Per-LP composite summary: gid_cover, n_gid_eff, n_gid, is_composite.
#' Composite status is a property of the crosswalk, never of a name string.
compute_gid_cover <- function(crosswalk, cover_threshold = 0.90) {
  stopifnot(all(c("lctn_pr", "gadm_level", "H") %in% names(crosswalk)))
  if (!nrow(crosswalk)) {
    return(tibble::tibble(
      lctn_pr = character(0),
      gadm_level = character(0),
      year = integer(0),
      P_l = numeric(0),
      n_gid = integer(0),
      gid_cover = numeric(0),
      n_gid_eff = numeric(0),
      is_composite = logical(0)
    ))
  }
  grp <- if ("year" %in% names(crosswalk)) {
    c("lctn_pr", "gadm_level", "year")
  } else {
    c("lctn_pr", "gadm_level")
  }
  crosswalk %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      P_l = sum(H, na.rm = TRUE),
      n_gid = dplyr::n_distinct(gid),
      .w = list(
        if (sum(H, na.rm = TRUE) > 0) H / sum(H, na.rm = TRUE) else NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gid_cover = vapply(
        .w,
        function(w) if (all(is.na(w))) NA_real_ else max(w),
        numeric(1)
      ),
      n_gid_eff = vapply(
        .w,
        function(w) if (all(is.na(w))) NA_real_ else 1 / sum(w^2),
        numeric(1)
      ),
      is_composite = !is.na(gid_cover) & gid_cover < cover_threshold
    ) %>%
    dplyr::select(-.w)
}

#' Derive each LP's GADM tiling depth from its own geometry.
#'
#' The corpus `admin_level` string is not GADM's admin depth (a corpus "ADM2"
#' can be a GADM ADM1), so the level is derived per LP:
#'   - population dominance decides which levels are admissible: level L is
#'     admissible iff pcov_L = max_g w[l->g] >= cover_threshold;
#'   - area decides which admissible level represents the LP (`tiebreak`);
#'   - if no level is admissible the LP is read at ADM1 and the downstream
#'     dominance test calls it composite.
#' Levels are compared explicitly (pcov_ADM1 >= pcov_ADM2 does not always
#' hold), and the decision is taken once per LP from the median cover across
#' years, so one LP is never read at two tilings.
#' Tiebreak modes:
#'   "strict" (default): the larger of acov_ADM1 / acov_ADM2; ties go to ADM2.
#'   "tolerant": ADM2 if acov_ADM2 > tiebreak_floor, else ADM2 if
#'     acov_ADM2 - acov_ADM1 > -tiebreak_tol, else ADM1. Each tier writes its
#'     own lp_level_reason ("both_admissible_acov_floor",
#'     "both_admissible_area_tiebreak_tolerant",
#'     "both_admissible_area_tiebreak_adm1").
#'
#' @param crosswalk bound crosswalk rows (lctn_pr, gid, gadm_level, year, H,
#'   piece_area_km2).
#' @param cover_threshold population dominance required for a level to be
#'   admissible. Shares the meaning of the threshold in `compute_gid_cover()`.
#' @param tiebreak "strict" (default) or "tolerant" — see above.
#' @param tiebreak_floor "tolerant" mode only: acov_ADM2 above this wins
#'   outright, no comparison to ADM1.
#' @param tiebreak_tol "tolerant" mode only: ADM2 still wins if it trails
#'   ADM1's acov by no more than this margin.
#' @return tibble, one row per lctn_pr: lp_level, lp_level_reason, and the
#'   diagnostics pcov_ADM1, pcov_ADM2, acov_ADM1, acov_ADM2, n_member_ADM1,
#'   n_member_ADM2.
derive_lp_tiling_level <- function(
  crosswalk,
  cover_threshold = 0.90,
  tiebreak = c("strict", "tolerant"),
  tiebreak_floor = 0.80,
  tiebreak_tol = 0.10
) {
  tiebreak <- match.arg(tiebreak)
  need <- c("lctn_pr", "gid", "gadm_level", "H")
  miss <- setdiff(need, names(crosswalk))
  if (length(miss)) {
    stop(
      "derive_lp_tiling_level(): crosswalk is missing column(s): ",
      paste(miss, collapse = ", ")
    )
  }
  if (
    !is.numeric(cover_threshold) ||
      length(cover_threshold) != 1L ||
      !is.finite(cover_threshold) ||
      cover_threshold < 0 ||
      cover_threshold > 1
  ) {
    stop(
      "derive_lp_tiling_level(): cover_threshold must be a single number in [0, 1]."
    )
  }

  empty <- tibble::tibble(
    lctn_pr = character(0),
    lp_level = character(0),
    lp_level_reason = character(0),
    pcov_ADM1 = numeric(0),
    pcov_ADM2 = numeric(0),
    acov_ADM1 = numeric(0),
    acov_ADM2 = numeric(0),
    n_member_ADM1 = numeric(0),
    n_member_ADM2 = numeric(0)
  )
  if (!nrow(crosswalk)) {
    return(empty)
  }

  # An area column is not strictly required: without it every level scores 0 on
  # the tie-break, which degrades to "prefer the finer level among admissible"
  # rather than erroring. Flagged rather than silent, because the tie-break is
  # what stops a population-concentrated province being read at ADM2.
  has_area <- "piece_area_km2" %in% names(crosswalk)
  if (!has_area) {
    warning(
      "derive_lp_tiling_level(): no `piece_area_km2` column; the area ",
      "tie-break between admissible levels is unavailable."
    )
  }

  cw <- crosswalk %>%
    dplyr::filter(is.finite(.data$H), .data$H > 0) %>%
    dplyr::mutate(
      lctn_pr = as.character(.data$lctn_pr),
      gadm_level = as.character(.data$gadm_level),
      .area = if (has_area) as.numeric(.data$piece_area_km2) else NA_real_,
      .year = if ("year" %in% names(crosswalk)) as.integer(.data$year) else 0L
    )
  if (!nrow(cw)) {
    return(empty)
  }

  per_level <- cw %>%
    dplyr::group_by(.data$lctn_pr, .data$gadm_level, .data$.year) %>%
    dplyr::summarise(
      n_member = dplyr::n(),
      pcov = max(.data$H) / sum(.data$H),
      acov = {
        s <- sum(.data$.area, na.rm = TRUE)
        if (is.finite(s) && s > 0) {
          max(.data$.area, na.rm = TRUE) / s
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    # Median across years: one decision per LP, immune to a single odd year.
    dplyr::group_by(.data$lctn_pr, .data$gadm_level) %>%
    dplyr::summarise(
      n_member = stats::median(.data$n_member),
      pcov = stats::median(.data$pcov),
      acov = stats::median(.data$acov, na.rm = TRUE),
      .groups = "drop"
    )

  wide <- per_level %>%
    tidyr::pivot_wider(
      names_from = "gadm_level",
      values_from = c("pcov", "acov", "n_member")
    )
  for (nm in c(
    "pcov_ADM1",
    "pcov_ADM2",
    "acov_ADM1",
    "acov_ADM2",
    "n_member_ADM1",
    "n_member_ADM2"
  )) {
    if (!nm %in% names(wide)) wide[[nm]] <- NA_real_
  }

  wide %>%
    dplyr::mutate(
      .ok1 = !is.na(.data$pcov_ADM1) & .data$pcov_ADM1 >= cover_threshold,
      .ok2 = !is.na(.data$pcov_ADM2) & .data$pcov_ADM2 >= cover_threshold,
      # NA area sorts below any real value, so a level that has no area
      # information never wins the tie-break by default.
      .a1 = dplyr::coalesce(.data$acov_ADM1, -Inf),
      .a2 = dplyr::coalesce(.data$acov_ADM2, -Inf),
      lp_level = if (tiebreak == "strict") {
        dplyr::case_when(
          .data$.ok1 & .data$.ok2 & .data$.a2 >= .data$.a1 ~ "ADM2",
          .data$.ok1 & .data$.ok2 ~ "ADM1",
          .data$.ok2 ~ "ADM2",
          .data$.ok1 ~ "ADM1",
          TRUE ~ "ADM1"
        )
      } else {
        dplyr::case_when(
          .data$.ok1 & .data$.ok2 & .data$.a2 > tiebreak_floor ~ "ADM2",
          .data$.ok1 &
            .data$.ok2 &
            (.data$.a2 - .data$.a1) > -tiebreak_tol ~ "ADM2",
          .data$.ok1 & .data$.ok2 ~ "ADM1",
          .data$.ok2 ~ "ADM2",
          .data$.ok1 ~ "ADM1",
          TRUE ~ "ADM1"
        )
      },
      lp_level_reason = if (tiebreak == "strict") {
        dplyr::case_when(
          .data$.ok1 & .data$.ok2 ~ "both_admissible_area_tiebreak",
          .data$.ok2 ~ "only_adm2_admissible",
          .data$.ok1 ~ "only_adm1_admissible",
          TRUE ~ "none_admissible_composite"
        )
      } else {
        dplyr::case_when(
          .data$.ok1 &
            .data$.ok2 &
            .data$.a2 > tiebreak_floor ~ "both_admissible_acov_floor",
          .data$.ok1 &
            .data$.ok2 &
            (.data$.a2 - .data$.a1) >
              -tiebreak_tol ~ "both_admissible_area_tiebreak_tolerant",
          .data$.ok1 & .data$.ok2 ~ "both_admissible_area_tiebreak_adm1",
          .data$.ok2 ~ "only_adm2_admissible",
          .data$.ok1 ~ "only_adm1_admissible",
          TRUE ~ "none_admissible_composite"
        )
      }
    ) %>%
    dplyr::select(
      "lctn_pr",
      "lp_level",
      "lp_level_reason",
      "pcov_ADM1",
      "pcov_ADM2",
      "acov_ADM1",
      "acov_ADM2",
      "n_member_ADM1",
      "n_member_ADM2"
    )
}

#' Reconciliation of the crosswalk against the independent whole-LP extraction.
#'
#' `P_l` is sum_g H over pieces; `P_direct` is a separate extraction over the
#' whole LP geometry. They agree only if the intersection geometry, the sliver
#' filter and the NA handling are all correct — which is why this is the primary
#' acceptance gate for Stage 1 and not a tautology.
reconcile_crosswalk <- function(crosswalk, lp_pop, tol = 0.01) {
  p_l <- crosswalk %>%
    dplyr::group_by(lctn_pr, gadm_level, year) %>%
    dplyr::summarise(P_l = sum(H, na.rm = TRUE), .groups = "drop")

  # A cross-border LP is extracted once per country; P_direct sums across tasks.
  p_d <- lp_pop %>%
    dplyr::group_by(lctn_pr, year) %>%
    dplyr::summarise(
      P_direct = sum(P_direct, na.rm = TRUE),
      cross_border = any(cross_border),
      .groups = "drop"
    )

  p_l %>%
    dplyr::left_join(p_d, by = c("lctn_pr", "year")) %>%
    dplyr::mutate(
      abs_diff = abs(P_l - P_direct),
      rel_diff = dplyr::if_else(P_direct > 0, abs_diff / P_direct, NA_real_),
      ok = !is.na(rel_diff) & rel_diff <= tol
    )
}

#' Detect location periods that share IDENTICAL crosswalk pieces with another
#' LP at the same GADM level -- the downstream-visible symptom of the upstream
#' "duplicate / parent-inherited geometry" defect (Taxonomy hands out the same
#' WKT to more than one location_period_id; e.g. 7 Conakry LPs, hash
#' ba1de5dd1a, all carrying pop = 1,632,251 and area = 413.92 km2). Upstream
#' provenance columns (`pop_geom_dup_n`, `pop_geom_dup_class`) require a
#' patched OutbreakExtractR re-extraction that has not run against this corpus
#' (see utils/population_qc.R's check_duplicate_geometry_provenance(), which no-ops
#' until those columns exist upstream). This is the same symptom made visible
#' from data already on disk: two LPs with genuinely distinct polygons do not,
#' by chance, produce the exact same (gid, H) set in every extracted year --
#' when they do, the crosswalk is repeating one polygon's mass under two (or
#' more) location_period_ids, and both would independently claim the full
#' population as their own denominator/composite split.
#'
#' Detect-and-quarantine, not a correction: there is no principled way to
#' decide, from H alone, which (if either) LP the shared polygon truly
#' belongs to (the join code is innocent; this is a Taxonomy data problem). Callers get
#' a flag to act on (report, exclude,
#' down-weight); this function never drops or reallocates anything itself.
#'
#' @param crosswalk bound crosswalk rows (lctn_pr, gid, gadm_level, year, H)
#' @param round_dp   rounding applied to H before comparing (float-noise guard)
#' @return tibble(lctn_pr, gadm_level, geom_dup_n, geom_dup_group); geom_dup_n
#'   == 1 and geom_dup_group == NA for an LP with no crosswalk twin. An LP
#'   absent from `crosswalk` at a level (nothing to compare) is absent here too
#'   -- left_join onto a fuller table yields NA, not a false "unique".
detect_geom_dup_lps <- function(crosswalk, round_dp = 3) {
  if (!nrow(crosswalk)) {
    return(tibble::tibble(
      lctn_pr = character(0),
      gadm_level = character(0),
      geom_dup_n = integer(0),
      geom_dup_group = character(0)
    ))
  }

  sig <- crosswalk %>%
    dplyr::filter(is.finite(H)) %>%
    dplyr::mutate(
      .piece = paste0(
        gid,
        "@",
        year,
        "=",
        format(round(H, round_dp), scientific = FALSE)
      )
    ) %>%
    dplyr::group_by(lctn_pr, gadm_level) %>%
    dplyr::summarise(
      .sig = paste(sort(.piece), collapse = "|"),
      .groups = "drop"
    )

  sig %>%
    dplyr::group_by(gadm_level, .sig) %>%
    dplyr::mutate(
      geom_dup_n = dplyr::n(),
      # A stable, deterministic group label (lowest lctn_pr in the cluster) --
      # same convention as the "lowest id wins" representative pick elsewhere
      # in this project (e.g. 03_prepare_siz.R's legacy id-unification step).
      geom_dup_group = if (dplyr::n() > 1L) {
        paste0(gadm_level[1], ":", sort(lctn_pr)[1])
      } else {
        NA_character_
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(lctn_pr, gadm_level, geom_dup_n, geom_dup_group)
}
