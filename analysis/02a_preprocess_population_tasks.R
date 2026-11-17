# 02a_preprocess_population_tasks.R — write the per-(country, year) task YAMLs
#   of the population-crosswalk SLURM array.
#
# Usage:
#   Rscript analysis/02a_preprocess_population_tasks.R -C <cfg.yml>
#   Rscript analysis/02a_preprocess_population_tasks.R -C <cfg.yml> \
#     --countries BDI,ETH --years 2020
#   Rscript analysis/02a_preprocess_population_tasks.R -C <cfg.yml> \
#     --task_unit country
#
# Details:
#   Called by analysis/bash/stages/02_population.sh before sbatch. Mirrors
#   07a_simulator_tasks.R, including its stdout contract: prints
#   TASKS_DIR=<path> and N_TASKS=<N> for bash capture. --task_unit country
#   writes one task per country with every year, for
#   02b_preprocess_population_run_country.R.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/population_core.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    "--countries",
    default = NULL,
    type = "character",
    help = "comma-separated ISO3 subset (default: the whole registry)"
  ),
  make_option(
    "--years",
    default = NULL,
    type = "character",
    help = "comma-separated year subset (default: population.years range)"
  ),
  make_option(
    "--task_unit",
    default = "country_year",
    type = "character",
    help = paste(
      "country_year: one task per (country, year), for",
      "02b_preprocess_population_run.R; country: one task per",
      "country with every year, for",
      "02b_preprocess_population_run_country.R [default: %default]"
    )
  )
)
opt <- make_options_from_config(option_list = option_list, model = NULL)
if (!opt$task_unit %in% c("country_year", "country")) {
  stop("--task_unit must be country_year or country, got: ", opt$task_unit)
}

pop <- opt$population
if (is.null(pop)) {
  stop(
    "Config resolved without a `population:` block. That should be impossible ",
    "— it is registered with fill_when_absent = TRUE in .CONFIG_TOPLEVEL_BLOCKS ",
    "(analysis/utils/config.R). See analysis/tests/test_config_block_registry.R."
  )
}

# 3. Country x year grid -------------------------------------------------------
# Why (country, year) and not (country): with all 16 years requested, a
# per-country task would serialise 16 WorldPop downloads and 16 extraction
# passes into one job. Per-(country, year):
#   * each WorldPop raster is owned by exactly ONE task, so the missing
#     cross-worker download lock in download_worldpop_constrained()
#     (OutbreakExtractR R/get_pop.R:1-68) becomes structurally unreachable;
#   * the UN adjustment factor is per (country, year) anyway
#     (estimate_adj_factors(), R/get_pop.R:70-133) — nothing is recomputed;
#   * resume granularity is one country-year rather than one country.
# The only shared resource left is the national-boundary cache; building it
# once before the array removes that race.
#
# --task_unit country (the per-country runner,
# 02b_preprocess_population_run_country.R) writes one task per country with
# every year instead. The download argument above does not hold on a cluster
# whose compute nodes cannot reach WorldPop: there the rasters are prestaged
# on the login node and no task downloads.
# The country's geometry is then built once instead of once per year.
iso3s <- if (!is.null(opt$countries) && nzchar(opt$countries)) {
  toupper(trimws(strsplit(opt$countries, ",")[[1]]))
} else {
  reg <- pop$country_registry
  if (is.null(reg) || !file.exists(reg)) {
    stop("population.country_registry not found: ", reg %||% "<NULL>")
  }
  sort(unique(toupper(read.csv(reg, stringsAsFactors = FALSE)$iso3)))
}
iso3s <- iso3s[nzchar(iso3s) & !is.na(iso3s)]
if (!length(iso3s)) {
  stop("no countries to plan")
}

years <- if (!is.null(opt$years) && nzchar(opt$years)) {
  as.integer(trimws(strsplit(opt$years, ",")[[1]]))
} else {
  yr <- as.integer(unlist(pop$years))
  seq.int(min(yr), max(yr))
}
years <- sort(unique(years[!is.na(years)]))
if (!length(years)) {
  stop("no years to plan")
}

tasks_dir <- make_run_file(opt, kind = "pop_tasks_dir")
parts_dir <- make_run_file(opt, kind = "pop_parts_dir")
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)

# A stale tasks/ directory silently changes the array size that the launcher
# computes from `ls | wc -l`, so clear it rather than merge into it.
old <- list.files(tasks_dir, pattern = "^task_\\d+\\.yml$", full.names = TRUE)
if (length(old)) {
  message("[pop_make_tasks] removing ", length(old), " stale task YAML(s)")
  unlink(old)
}

# WPP national totals are cached HERE, on the login node, because this is the
# last step before submission that is allowed to need OutbreakExtractR. Workers
# and every downstream script read the CSV. Missing OutbreakExtractR is a
# warning, not an error: the only per-task consequence is natl_ratio = NA in the
# manifest, and 02c_preprocess_population_compile.R recomputes that check globally.
wpp_path <- make_run_file(opt, kind = "pop_wpp")
tryCatch(
  {
    write_wpp_national_totals(wpp_path)
    message("[pop_make_tasks] WPP2024 national totals: ", wpp_path)
  },
  error = function(e) {
    warning(
      "[pop_make_tasks] could not cache WPP totals (",
      conditionMessage(e),
      ") — manifests will carry natl_ratio = NA."
    )
  }
)

grid <- expand.grid(
  year = years,
  iso3 = iso3s,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$iso3, grid$year), , drop = FALSE]
part_dir_of <- function(iso3, year) {
  file.path(parts_dir, sprintf("iso3=%s", iso3), sprintf("year=%d", year))
}

message(sprintf(
  "[pop_make_tasks] %d countries x %d years = %d %s",
  length(iso3s),
  length(years),
  if (opt$task_unit == "country") length(iso3s) else nrow(grid),
  if (opt$task_unit == "country") "country tasks" else "tasks"
))
message("[pop_make_tasks] tasks_dir = ", tasks_dir)

# 4. Write task YAMLs ----------------------------------------------------------
# Every path the worker needs is frozen into the task file, so a worker never
# re-resolves the config. This is what makes a task reproducible in isolation
# and keeps `--resume` from silently changing inputs between submissions.
if (opt$task_unit == "country") {
  paths <- character(length(iso3s))
  for (i in seq_along(iso3s)) {
    task <- list(
      task_id = i,
      iso3 = iso3s[i],
      years = as.integer(years),
      part_dirs = vapply(
        years,
        function(y) part_dir_of(iso3s[i], y),
        character(1)
      ),
      gadm_path = pop$gadm_path,
      lp_shapefiles = pop$lp_shapefiles,
      raster_dir = pop$raster_dir,
      boundary_cache = pop$boundary_cache,
      wpp_table = make_run_file(opt, kind = "pop_wpp"),
      sliver_mass_frac = pop$sliver_mass_frac,
      sliver_mass_min = pop$sliver_mass_min,
      require_boundary = TRUE
    )
    f <- file.path(tasks_dir, sprintf("task_%04d.yml", i))
    yaml::write_yaml(task, f)
    paths[i] <- f
  }
} else {
  paths <- character(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    part_dir <- part_dir_of(grid$iso3[i], grid$year[i])
    task <- list(
      task_id = i,
      iso3 = grid$iso3[i],
      year = as.integer(grid$year[i]),
      part_dir = part_dir,
      gadm_path = pop$gadm_path,
      lp_shapefiles = pop$lp_shapefiles,
      raster_dir = pop$raster_dir,
      boundary_cache = pop$boundary_cache,
      wpp_table = make_run_file(opt, kind = "pop_wpp"),
      sliver_mass_frac = pop$sliver_mass_frac,
      sliver_mass_min = pop$sliver_mass_min,
      require_boundary = TRUE
    )
    f <- file.path(tasks_dir, sprintf("task_%04d.yml", i))
    yaml::write_yaml(task, f)
    paths[i] <- f
  }
}

message(sprintf(
  "[pop_make_tasks] wrote %d task YAML(s) to %s",
  length(paths),
  tasks_dir
))

# 5. Machine-readable output for bash capture (stdout only) --------------------
cat(sprintf("TASKS_DIR=%s\n", tasks_dir))
cat(sprintf("N_TASKS=%d\n", length(paths)))
