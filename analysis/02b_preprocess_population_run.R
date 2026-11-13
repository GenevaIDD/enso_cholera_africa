# 02b_preprocess_population_run.R — dispatcher for the LP x gid population extraction.
#
# Phase 02, orchestration verb `run`; the logic lives in the sourced kernel
# analysis/utils/population_core.R and the reduction in
# analysis/02c_preprocess_population_compile.R.
#
# Two modes:
#   --task <task.yml>   one (country, year) — what the SLURM array element runs
#   --local             every task in tasks/ in-process, sequentially
#
# A task writes into
#   generated_data/population/parts/iso3=<ISO3>/year=<YYYY>/
#     crosswalk.qs   lctn_pr | gid | gadm_level | gid_iso3 | year | H |
#                    piece_area_km2 | n_parts | geom_repair | w_dropped
#     lp_pop.qs      lctn_pr | extract_iso3 | year | P_direct | lp_area_km2 |
#                    cross_border
#     gid_pop.qs     gid | gadm_level | iso3 | year | Q | area_km2
#     manifest.yml   provenance + natl_ratio + runtime
#     _SUCCESS       written LAST, so a half-written part is never resumed over
#
# Usage:
#   Rscript analysis/02b_preprocess_population_run.R -C <cfg.yml> --task <tasks/task_0001.yml>
#   Rscript analysis/02b_preprocess_population_run.R -C <cfg.yml> --local
#   Rscript analysis/02b_preprocess_population_run.R -C <cfg.yml> --local --countries BDI --years 2020

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

option_list <- list(
  make_option("--task", default = NULL, type = "character",
              help = "path to a single task YAML"),
  make_option("--local", action = "store_true", default = FALSE,
              help = "run every task in tasks/ sequentially, in-process"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "re-run a task even if its partition is complete and fresh"),
  make_option("--countries", default = NULL, type = "character",
              help = "--local only: comma-separated ISO3 subset"),
  make_option("--years", default = NULL, type = "character",
              help = "--local only: comma-separated year subset")
)
opt <- make_options_from_config(option_list = option_list, model = NULL)

if (is.null(opt$task) && !isTRUE(opt$local))
  stop("one of --task <task.yml> or --local is required")

# GEOS behaviour must match the validation the geometries were built under
# (analysis/00_ingest_outbreakextractr.R), or st_intersection() disagrees with
# the polygons on disk about which edges are valid.
sf::sf_use_s2(FALSE)

# ---------------------------------------------------------------------------
# One task
# ---------------------------------------------------------------------------

run_task_file <- function(task_path, force = FALSE) {
  tk <- yaml::read_yaml(task_path)
  tag <- sprintf("%s %d", tk$iso3, tk$year)

  sentinel <- file.path(tk$part_dir, "_SUCCESS")
  # Skip only a part built from THIS LP corpus (popPartFresh): a part of an
  # older outbreak vintage is rebuilt, not reused.
  lp_md5 <- unname(tools::md5sum(path.expand(tk$lp_shapefiles)))
  if (!isTRUE(force) && popPartFresh(tk$part_dir, tk$lp_shapefiles, lp_md5)) {
    message("[pop] ", tag, ": fresh _SUCCESS present, skipping (use --force to redo)")
    return(invisible(FALSE))
  }
  if (file.exists(sentinel)) {
    message("[pop] ", tag, ": rebuilding a part of another LP corpus (or --force)")
    file.remove(sentinel)
  }

  message("[pop] ", tag, ": start")
  t0 <- Sys.time()
  elapsed <- function() sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

  message("[pop] ", tag, ": loading GADM units from ", tk$gadm_path, " (", elapsed(), ")")
  gadm  <- loadGadmUnits(tk$gadm_path, tk$iso3)
  message("[pop] ", tag, ": GADM units loaded, ", nrow(gadm), " rows (", elapsed(), ")")

  message("[pop] ", tag, ": loading LP geometries from ", tk$lp_shapefiles, " (", elapsed(), ")")
  lp_sf <- loadLpGeometries(tk$lp_shapefiles)
  message("[pop] ", tag, ": LP geometries loaded, ", nrow(lp_sf), " rows (", elapsed(), ")")

  # A CRS mismatch would make st_intersection() return empty for every pair —
  # i.e. an empty crosswalk and a silently zero P_l, not an error.
  if (!identical(sf::st_crs(lp_sf), sf::st_crs(gadm))) {
    message("[pop] ", tag, ": transforming LP geometries to the GADM CRS (", elapsed(), ")")
    lp_sf <- sf::st_transform(lp_sf, sf::st_crs(gadm))
    message("[pop] ", tag, ": CRS transform done (", elapsed(), ")")
  }

  message("[pop] ", tag, ": loading WPP national totals (", elapsed(), ")")
  wpp <- loadWppNationalTotals(tk$wpp_table)
  message("[pop] ", tag, ": WPP totals loaded (", elapsed(), ")")

  message("[pop] ", tag, ": entering runPopulationTask() (", elapsed(), ")")
  res <- runPopulationTask(
    iso3               = tk$iso3,
    year               = tk$year,
    gadm_sf            = gadm,
    lp_sf              = lp_sf,
    raster_dir         = tk$raster_dir,
    boundary_cache_dir = tk$boundary_cache,
    sliver_mass_frac   = tk$sliver_mass_frac,
    sliver_mass_min    = tk$sliver_mass_min,
    wpp_total          = wppTotal(wpp, tk$iso3, tk$year),
    require_boundary   = isTRUE(tk$require_boundary),
    quiet              = FALSE,
    progress_tag       = paste0("[pop] ", tag, ": ")
  )
  message("[pop] ", tag, ": runPopulationTask() returned (", elapsed(), ")")

  dir.create(tk$part_dir, recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(res$crosswalk, file.path(tk$part_dir, "crosswalk.qs"))
  qs2::qs_save(res$lp_pop,    file.path(tk$part_dir, "lp_pop.qs"))
  qs2::qs_save(res$gid_pop,   file.path(tk$part_dir, "gid_pop.qs"))

  man <- res$manifest
  man$task_id       <- tk$task_id
  man$task_yaml     <- normalizePath(task_path, mustWork = FALSE)
  man$gadm_path     <- tk$gadm_path
  man$gadm_mtime    <- format(file.mtime(tk$gadm_path))
  man$lp_shapefiles <- tk$lp_shapefiles
  man$lp_mtime      <- format(file.mtime(tk$lp_shapefiles))
  man$lp_md5        <- lp_md5
  man$oer_version   <- tryCatch(
    as.character(utils::packageVersion("OutbreakExtractR")),
    error = function(e) NA_character_)
  yaml::write_yaml(man, file.path(tk$part_dir, "manifest.yml"))

  # Sentinel LAST. --resume scans for it, so writing it before the artefacts
  # would let a killed task be treated as complete.
  file.create(sentinel)

  message(sprintf("[pop] %s: done in %.1fs — %d pieces, %d gids, %d LPs%s",
                  tag, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  nrow(res$crosswalk), nrow(res$gid_pop), nrow(res$lp_pop),
                  if (is.finite(man$natl_ratio))
                    sprintf(", natl_ratio %.3f", man$natl_ratio) else ""))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

if (!is.null(opt$task)) {
  if (!file.exists(opt$task)) stop("task YAML not found: ", opt$task)
  run_task_file(opt$task, force = isTRUE(opt$force))
  quit(status = 0)
}

tasks_dir <- makeRunFile(opt, kind = "pop_tasks_dir")
if (!dir.exists(tasks_dir))
  stop("tasks dir not found: ", tasks_dir,
       " — run analysis/02a_preprocess_population_tasks.R first.")
task_files <- sort(list.files(tasks_dir, pattern = "^task_\\d+\\.yml$",
                              full.names = TRUE))
if (!length(task_files)) stop("no task YAMLs under ", tasks_dir)

subset_keep <- function(f) {
  tk <- yaml::read_yaml(f)
  ok <- TRUE
  if (!is.null(opt$countries) && nzchar(opt$countries))
    ok <- ok && tk$iso3 %in% toupper(trimws(strsplit(opt$countries, ",")[[1]]))
  if (!is.null(opt$years) && nzchar(opt$years))
    ok <- ok && tk$year %in% as.integer(trimws(strsplit(opt$years, ",")[[1]]))
  ok
}
if ((!is.null(opt$countries) && nzchar(opt$countries)) ||
    (!is.null(opt$years) && nzchar(opt$years))) {
  task_files <- task_files[vapply(task_files, subset_keep, logical(1))]
}

message("[pop] --local: ", length(task_files), " task(s)")
n_ok <- 0L; failed <- character(0)
for (f in task_files) {
  res <- tryCatch(run_task_file(f, force = isTRUE(opt$force)),
                  error = function(e) {
                    message("[pop] FAILED ", basename(f), ": ", conditionMessage(e))
                    NULL
                  })
  if (is.null(res)) failed <- c(failed, basename(f)) else n_ok <- n_ok + 1L
}
message(sprintf("[pop] --local finished: %d ok, %d failed", n_ok, length(failed)))
if (length(failed)) {
  message("[pop] failed tasks: ", paste(failed, collapse = ", "))
  quit(status = 1)
}
