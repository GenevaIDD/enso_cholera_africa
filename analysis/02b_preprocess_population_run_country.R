# 02b_preprocess_population_run_country.R — run the LP x gid population
#   extraction one country (all years) at a time.
#
# Usage:
#   Rscript analysis/02b_preprocess_population_run_country.R -C <cfg.yml> \
#     --task <tasks/task_0001.yml>
#   Rscript analysis/02b_preprocess_population_run_country.R -C <cfg.yml> \
#     --local --countries BDI
#
# Details:
#   Reads the task YAMLs 02a writes under --task_unit country and calls
#   population_core.R::run_population_country_task(), the years in parallel on
#   --workers forked processes (default: the Slurm allocation's CPUs). Output
#   is the per-year runner's, partition for partition, so
#   02c_preprocess_population_compile.R is unchanged; fresh years are skipped.
#   Where compute nodes cannot download, the WorldPop rasters and the
#   boundary cache are prestaged first.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(qs2)
  library(sf)
  library(dplyr)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/population_core.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    "--task",
    default = NULL,
    type = "character",
    help = "path to a single country task YAML"
  ),
  make_option(
    "--local",
    action = "store_true",
    default = FALSE,
    help = "run every task in tasks/ sequentially, in-process"
  ),
  make_option(
    "--force",
    action = "store_true",
    default = FALSE,
    help = "re-run every year even if its partition is complete and fresh"
  ),
  make_option(
    "--no_gid_cache",
    action = "store_true",
    default = FALSE,
    help = "ignore the Q_g cache (extract every role; the cache is still written)"
  ),
  make_option(
    "--countries",
    default = NULL,
    type = "character",
    help = "--local only: comma-separated ISO3 subset"
  ),
  make_option(
    "--workers",
    default = 0L,
    type = "integer",
    help = paste(
      "years extracted at once (forked); 0 = future::availableCores(),",
      "i.e. the Slurm allocation's CPUs [default: %default]"
    )
  )
)
opt <- make_options_from_config(option_list = option_list, model = NULL)

if (is.null(opt$task) && !isTRUE(opt$local)) {
  stop("one of --task <task.yml> or --local is required")
}

# GEOS behaviour must match the validation the geometries were built under
# (analysis/00_ingest_outbreakextractr.R), as in 02b_preprocess_population_run.R.
sf::sf_use_s2(FALSE)

gid_cache_dir <- make_run_file(opt, kind = "pop_gid_cache_dir")
n_workers <- if (isTRUE(opt$workers > 0L)) {
  opt$workers
} else {
  as.integer(future::availableCores())
}

# 3. One country ---------------------------------------------------------------
# Per country, every year at once: GADM units, LP corpus, footprint and
# LP x gid pieces are built once. Q_g is cached per (country, year) under
# make_run_file(opt, "pop_gid_cache_dir"), keyed on GADM, raster, boundary and
# extractor version (not the outbreak vintage). Output partitions match the
# per-year runner's (generated_data/population/parts/iso3=/year=/); a fresh
# partition is skipped. Rasters and boundary cache must be prestaged.

run_country_task_file <- function(task_path, force = FALSE, use_cache = TRUE) {
  tk <- yaml::read_yaml(task_path)
  if (is.null(tk$years)) {
    stop(
      task_path,
      " is a per-(country, year) task; run it with ",
      "02b_preprocess_population_run.R, or regenerate with 02a --task_unit country"
    )
  }
  years <- as.integer(unlist(tk$years))
  part_dirs <- unlist(tk$part_dirs)
  tag <- tk$iso3
  t0 <- Sys.time()
  elapsed <- function() {
    sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }

  # Skip only parts built from THIS LP corpus.
  lp_md5 <- unname(tools::md5sum(path.expand(tk$lp_shapefiles)))
  todo <- force |
    !vapply(
      part_dirs,
      pop_part_fresh,
      logical(1),
      lp_shapefiles = tk$lp_shapefiles,
      lp_md5 = lp_md5
    )
  if (!any(todo)) {
    message("[pop] ", tag, ": every year fresh, skipping (use --force to redo)")
    return(invisible(FALSE))
  }
  years <- years[todo]
  part_dirs <- part_dirs[todo]
  for (d in part_dirs) {
    s <- file.path(d, "_SUCCESS")
    if (file.exists(s)) file.remove(s)
  }
  message(
    "[pop] ",
    tag,
    ": ",
    length(years),
    " year(s) to build: ",
    paste(years, collapse = ",")
  )

  gadm <- load_gadm_units(tk$gadm_path, tk$iso3)
  message(
    "[pop] ",
    tag,
    ": GADM units loaded, ",
    nrow(gadm),
    " rows (",
    elapsed(),
    ")"
  )
  lp_sf <- load_lp_geometries(tk$lp_shapefiles)
  message(
    "[pop] ",
    tag,
    ": LP geometries loaded, ",
    nrow(lp_sf),
    " rows (",
    elapsed(),
    ")"
  )
  # A CRS mismatch would make st_intersection() return empty for every pair.
  if (!identical(sf::st_crs(lp_sf), sf::st_crs(gadm))) {
    lp_sf <- sf::st_transform(lp_sf, sf::st_crs(gadm))
  }
  wpp <- load_wpp_national_totals(tk$wpp_table)

  keys <- lapply(years, function(y) {
    pop_gid_cache_key(
      tk$iso3,
      y,
      tk$gadm_path,
      tk$raster_dir,
      tk$boundary_cache
    )
  })
  cache <- if (use_cache) {
    lapply(keys, function(k) read_gid_pop_cache(gid_cache_dir, k))
  } else {
    list()
  }
  names(cache) <- if (use_cache) as.character(years) else NULL
  cache <- cache[!vapply(cache, is.null, logical(1))]
  message("[pop] ", tag, ": Q_g cache hits ", length(cache), "/", length(years))

  res <- run_population_country_task(
    iso3 = tk$iso3,
    years = years,
    gadm_sf = gadm,
    lp_sf = lp_sf,
    raster_dir = tk$raster_dir,
    boundary_cache_dir = tk$boundary_cache,
    sliver_mass_frac = tk$sliver_mass_frac,
    sliver_mass_min = tk$sliver_mass_min,
    wpp_totals = vapply(
      years,
      function(y) wpp_national_total(wpp, tk$iso3, y),
      numeric(1)
    ),
    gid_pop_cache = cache,
    require_boundary = isTRUE(tk$require_boundary),
    quiet = FALSE,
    progress_tag = paste0("[pop] ", tag, " "),
    workers = n_workers
  )

  for (i in seq_along(years)) {
    r <- res[[as.character(years[i])]]
    d <- part_dirs[i]
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    qs2::qs_save(r$crosswalk, file.path(d, "crosswalk.qs"))
    qs2::qs_save(r$lp_pop, file.path(d, "lp_pop.qs"))
    qs2::qs_save(r$gid_pop, file.path(d, "gid_pop.qs"))
    if (!isTRUE(r$manifest$gid_pop_cached)) {
      write_gid_pop_cache(gid_cache_dir, keys[[i]], r$gid_pop)
    }

    man <- r$manifest
    man$task_id <- tk$task_id
    man$task_yaml <- normalizePath(task_path, mustWork = FALSE)
    man$task_unit <- "country"
    man$gadm_path <- tk$gadm_path
    man$gadm_mtime <- format(file.mtime(tk$gadm_path))
    man$lp_shapefiles <- tk$lp_shapefiles
    man$lp_mtime <- format(file.mtime(tk$lp_shapefiles))
    man$lp_md5 <- lp_md5
    man$oer_version <- keys[[i]]$oer_version
    yaml::write_yaml(man, file.path(d, "manifest.yml"))
    # Sentinel LAST, per year.
    file.create(file.path(d, "_SUCCESS"))
  }

  message(sprintf(
    "[pop] %s: done in %s — %d year(s), %d Q_g from cache",
    tag,
    elapsed(),
    length(years),
    length(cache)
  ))
  invisible(TRUE)
}

# 4. Dispatch ------------------------------------------------------------------

if (!is.null(opt$task)) {
  if (!file.exists(opt$task)) {
    stop("task YAML not found: ", opt$task)
  }
  run_country_task_file(
    opt$task,
    force = isTRUE(opt$force),
    use_cache = !isTRUE(opt$no_gid_cache)
  )
  quit(status = 0)
}

tasks_dir <- make_run_file(opt, kind = "pop_tasks_dir")
task_files <- sort(list.files(
  tasks_dir,
  pattern = "^task_\\d+\\.yml$",
  full.names = TRUE
))
if (!length(task_files)) {
  stop(
    "no task YAMLs under ",
    tasks_dir,
    " — run analysis/02a_preprocess_population_tasks.R --task_unit country first."
  )
}
if (!is.null(opt$countries) && nzchar(opt$countries)) {
  keep <- toupper(trimws(strsplit(opt$countries, ",")[[1]]))
  task_files <- task_files[vapply(
    task_files,
    function(f) {
      yaml::read_yaml(f)$iso3 %in% keep
    },
    logical(1)
  )]
}

message("[pop] --local: ", length(task_files), " country task(s)")
failed <- character(0)
for (f in task_files) {
  ok <- tryCatch(
    run_country_task_file(
      f,
      force = isTRUE(opt$force),
      use_cache = !isTRUE(opt$no_gid_cache)
    ),
    error = function(e) {
      message("[pop] FAILED ", basename(f), ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ok)) failed <- c(failed, basename(f))
}
message(sprintf(
  "[pop] --local finished: %d ok, %d failed",
  length(task_files) - length(failed),
  length(failed)
))
if (length(failed)) {
  message("[pop] failed tasks: ", paste(failed, collapse = ", "))
  quit(status = 1)
}
