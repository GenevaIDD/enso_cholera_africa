# 03_prepare_lenssen_analogs.R — build the Lenssen piControl analog cache
# for Scenario B of the ENSO counterfactual forecast.
#
# Reads `step01_InitVerifInfo.Rda` from each of the 11 CMIP6 model zips in
# Zenodo record 10045687 (Lenssen et al. 2024, GRL 10.1029/2023GL106988).
# Each Rda contains `nino34` (detrended, 3-month rolling-mean ONI-equivalent)
# and `timeMap` (year, month, decimal_year) for the piControl run.
# No NetCDF extraction is needed.
#
# For each model:
#   1. Selects piControl months matching target_month (May) and
#      closest in ONI to the observed value at horizon_start (2023-05).
#   2. Extracts 0..lead_months forward trajectory per analog.
#   3. Rescales variance to match observed ONI SD.
#
# Outputs:
#   generated_data/lenssen_nino34_analogs.rds   (cache read by sampleEnsoLenssenScenarios)
#   generated_data/lenssen_nino34_analogs_meta.yml
#
# Kept side-by-side with 03c / 03d.

# Preamble ----------------------------------------------------------------

library(dplyr)
library(purrr)
library(lubridate)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")

# Options -----------------------------------------------------------------

option_list <- list(
  make_option("--horizon_start",
              default = "2023-05-01", type = "character",
              help    = "Analog initialization date (must be in ONI cache) [default: %default]"),
  make_option("--n_analogs",
              default = 50L, type = "integer",
              help    = "Max analogs per model after min-sep filter [default: %default]"),
  make_option("--lead_months",
              default = 60L, type = "integer",
              help    = "Forward trajectory length (0..lead_months = lead_months+1 cols) [default: %default]"),
  make_option("--min_sep",
              default = 12L, type = "integer",
              help    = "Min months between analog init indices within a model [default: %default]"),
  make_option("--zenodo_dir",
              default = "data/lenssen/", type = "character",
              help    = "Directory containing the 11 model zip files [default: %default]"),
  make_option("--refresh",
              default = FALSE, action = "store_true",
              help    = "Re-extract Rda files even if already present [default: %default]"),
  make_option("--refresh_enso",
              default = FALSE, action = "store_true",
              help    = "Re-download ONI cache [default: %default]"),
  make_option("--oni_cache_dir",
              default = NA_character_, type = "character",
              help    = paste("Directory of the observed ONI cache (fetchEnsoIndex());",
                              "NA = --output_dir, the earlier behaviour. The event",
                              "attribution builds one cache per event in its own",
                              "--output_dir and reads the shared generated_data/ ONI",
                              "cache through this [default: %default]")),
  make_option("--output_dir",
              default = "generated_data/", type = "character",
              help    = "Directory for output RDS and metadata YAML [default: %default]"),
  # -- ENSO-state analog sets (consumed by sampleEnsoLenssenScenarios via
  #    simulator.enso.counterfactual.analog_set).  Classification is done here,
  #    once, and persisted in the cache: it is a property of the pool, it keeps the
  #    sampler a pure row-filter, and it makes the sets auditable from the meta YAML.
  make_option("--set_stat",
              default = "peak_djf", type = "character",
              help    = paste("Expected analog-set statistic name; must agree with",
                              "--set_peak_months, or 'auto' to accept whatever",
                              "--set_peak_months implies [default: %default]")),
  make_option("--set_peak_months",
              default = "12,1,2", type = "character",
              help    = "Comma-separated calendar months defining the peak season [default: %default]"),
  make_option("--set_thresholds",
              default = paste(ENSO_PHASE_THRESHOLDS, collapse = ","), type = "character",
              help    = paste("Comma-separated nina,elnino ONI cut-offs (inclusive).",
                              "Default is the operational CPC pair [default: %default]")),
  make_option("--set_peak_season_year",
              default = NA_integer_, type = "integer",
              help    = paste("Calendar year of the FIRST month of the target peak season",
                              "(e.g. 2026 for DJF 2026/27); NA = first complete season",
                              "inside the horizon [default: %default]")),
  make_option("--min_set_size",
              default = 30L, type = "integer",
              help    = "Warn (legacy) / report (library) analog sets under this size [default: %default]"),
  # -- Library matcher (deploy.enso_counterfactual). `legacy`
  #    (default) is the builder below; `library` draws the pool from
  #    the ENSO continuation library (utils/enso.R::buildLenssenCache()).
  make_option("--match",
              default = "legacy", type = "character",
              help    = "Matcher: legacy | library [default: %default]"),
  make_option("--library",
              default = "generated_data/enso_continuation_library.qs", type = "character",
              help    = "library: the continuation library (03_prepare_enso_continuation_library.R) [default: %default]"),
  make_option("--blend_months",
              default = 3, type = "double",
              help    = "library: faded-join e-folding months [default: %default]"),
  make_option("--exclude_models",
              default = "", type = "character",
              help    = "library: comma-separated models to drop (as deploy.enso_continuation.exclude_models) [default: none]"),
  make_option("--set_winters",
              default = "primary", type = "character",
              help    = "library: primary | all_read (also classify every DJF read up to --set_read_end) [default: %default]"),
  make_option("--set_read_end",
              default = "", type = "character",
              help    = "library, all_read: last ONI month the reported months read (YYYY-MM-DD) [default: none]"),
  make_option("--set_other_winters",
              default = "neutral", type = "character",
              help    = "library: neutral (other winters must be neutral, else 'mixed') | any [default: %default]"),
  make_option("--seam_at",
              default = "launch", type = "character",
              help    = paste("library: where the analog paths join the observed ONI: launch (column 1 =",
                              "the observed launch-month ONI, offset fading from there) | launch_minus_1",
                              "(joined at the month before) [default: %default]"))
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

horizon_start <- as.Date(opt$horizon_start)
target_month  <- as.integer(format(horizon_start, "%m"))   # 5 for May

message("-- Lenssen piControl analog cache")
message("   horizon_start : ", horizon_start)
message("   target_month  : ", target_month)
message("   n_analogs     : ", opt$n_analogs)
message("   lead_months   : ", opt$lead_months, "  (trajectories = 0:", opt$lead_months, ")")
message("   min_sep       : ", opt$min_sep)
message("   zenodo_dir    : ", opt$zenodo_dir)

set_peak_months <- as.integer(trimws(strsplit(opt$set_peak_months, ",")[[1L]]))
set_thresholds  <- as.numeric(trimws(strsplit(opt$set_thresholds, ",")[[1L]]))
if (length(set_thresholds) != 2L || any(!is.finite(set_thresholds)))
  stop("--set_thresholds must be two finite numbers, 'nina,elnino' (e.g. '-0.5,0.5').")
names(set_thresholds) <- c("nina", "elnino")
set_peak_season_year <- if (is.na(opt$set_peak_season_year)) NULL else
  as.integer(opt$set_peak_season_year)

message("   set peak      : months ", paste(set_peak_months, collapse = ","),
        "  thresholds nina<=", set_thresholds[["nina"]],
        " / elnino>=", set_thresholds[["elnino"]],
        "  season_year ", if (is.null(set_peak_season_year)) "auto" else set_peak_season_year)

# Model list (matches zip filenames in zenodo record 10045687; utils/enso.R)
model_names <- LENSSEN_MODELS

# Step 1: Observed ONI at horizon_start -----------------------------------

message("-- Loading observed ONI")
oni_cache_dir <- if (is.na(opt$oni_cache_dir) || !nzchar(opt$oni_cache_dir)) opt$output_dir else
  opt$oni_cache_dir
obs_oni <- fetchEnsoIndex(
  source    = "oni",
  cache_dir = oni_cache_dir,
  refresh   = isTRUE(opt$refresh_enso)
)

# Library matcher ---------------------------------------------------------
if (identical(opt$match, "library")) {
  lib <- loadEnsoContinuationLibrary(opt$library)
  oni_lkp <- stats::setNames(as.numeric(obs_oni$enso_index), format(as.Date(obs_oni$date)))
  seam <- format(seq(horizon_start, by = "-1 month", length.out = 2L)[2L])
  seam_state <- ensoSeamState(oni_lkp, seam)
  excl <- trimws(strsplit(opt$exclude_models, ",")[[1L]]); excl <- excl[nzchar(excl)]
  cf <- ensoCounterfactualCfg(list(match = "library", n_analogs = opt$n_analogs, min_sep_months = opt$min_sep,
                                   blend_months = opt$blend_months, winters = opt$set_winters,
                                   other_winters = opt$set_other_winters, min_set_size = opt$min_set_size,
                                   seam_at = opt$seam_at))
  read_end <- if (nzchar(opt$set_read_end)) as.Date(opt$set_read_end) else NULL
  if (identical(cf$winters, "all_read") && is.null(read_end))
    stop("--set_winters all_read needs --set_read_end.")
  rw <- if (!is.null(read_end)) ensoReadWinters(horizon_start, read_end, set_peak_months) else integer(0)
  prim <- set_peak_season_year %||% rw[1L]
  if (is.null(prim) || is.na(prim)) prim <- NULL
  other <- if (identical(cf$winters, "all_read")) setdiff(as.integer(rw), prim) else integer(0)
  spec <- lenssenCacheSpec(horizon_start, opt$lead_months + 1L, cf, lib, seam_state,
                           peak_months = set_peak_months, thresholds = set_thresholds,
                           primary_years = prim, other_years = other, read_end = read_end,
                           exclude_models = excl, launch_oni = ensoLaunchOni(oni_lkp, horizon_start, cf$seam_at))
  out_rds <- makeLenssenAnalogFile(opt$output_dir)
  cache <- buildLenssenCache(spec, obs_oni, lib, out_path = out_rds)
  message("-- library cache -> ", out_rds, "  (digest ", cache$meta$digest, ")")
  message(sprintf("   seam %s: ONI %.3f, tendency %+.3f, ONI 6 months earlier %+.3f",
                  seam, seam_state$v, seam_state$tend, seam_state$phase))
  message(sprintf("   seam_at %s%s", cf$seam_at,
                  if (identical(cf$seam_at, "launch")) sprintf(" (column 1 pinned to the observed %s ONI %.3f)",
                                                              format(horizon_start), spec$launch_oni) else ""))
  if (length(attr(rw, "partial")))
    message("   winters overlapping the launch (observed in part, not classified): ",
            paste(attr(rw, "partial"), collapse = ", "))
  for (s in names(cache$meta$sets)) {
    x <- cache$meta$sets[[s]]
    message(sprintf("   %-8s n=%4d  seam error %+.3f  winter means %s", s, x$n, x$seam_error %||% NA_real_,
                    paste(sprintf("%s:%+.2f", names(x$winter_means), unlist(x$winter_means)), collapse = " ")))
    if (!identical(s, "all") && x$n < opt$min_set_size)
      warning("Analog set '", s, "' has only ", x$n, " trajectories (< --min_set_size ", opt$min_set_size, ").")
  }
  message("\nDone.")
  quit(save = "no", status = 0)
}

target_oni <- obs_oni |>
  dplyr::filter(date == horizon_start) |>
  dplyr::pull(enso_index)

if (length(target_oni) != 1L || is.na(target_oni)) {
  stop("ONI value for ", horizon_start, " not found in cache. ",
       "Run fetchEnsoIndex(source='oni', refresh=TRUE) to update.")
}
message("   observed ONI at ", horizon_start, ": ", round(target_oni, 3))

obs_sd <- sd(obs_oni$enso_index, na.rm = TRUE)
message("   observed ONI SD: ", round(obs_sd, 3))

# Step 2: Per-model functions ---------------------------------------------

# In utils/enso.R: lenssenExtractRda(),
# .lenssenSelectAnalogsLegacy(), .lenssenExtractTrajectoriesLegacy(),
# .lenssenRescaleLegacy().

# Step 3: Process all models ----------------------------------------------

extract_base  <- file.path(opt$zenodo_dir, "extracted")
scale_factors <- numeric(length(model_names))
names(scale_factors) <- model_names

all_results <- vector("list", length(model_names))
names(all_results) <- model_names

# Per-analog provenance.  Keeping `init_idx` (together with the source model)
# lets downstream work audit single-model influence on a set and, later, build
# trajectory-clustered intervals.
all_init <- vector("list", length(model_names))
names(all_init) <- model_names

for (model in model_names) {
  message("\n-- Processing ", model)

  rda_path <- tryCatch(
    lenssenExtractRda(model, opt$zenodo_dir, extract_base, isTRUE(opt$refresh)),
    error = function(e) { warning(conditionMessage(e)); NULL }
  )
  if (is.null(rda_path)) {
    message("   SKIPPED (zip not found or extraction failed)")
    next
  }

  e <- new.env(parent = emptyenv())
  load(rda_path, envir = e)

  if (!all(c("nino34", "timeMap") %in% ls(e))) {
    warning(model, ": step01_InitVerifInfo.Rda missing 'nino34' or 'timeMap' — skipped")
    next
  }

  message("   piControl length: ", length(e$nino34), " months  (",
          round(length(e$nino34) / 12), " years)")

  init_idx <- .lenssenSelectAnalogsLegacy(
    nino34       = e$nino34,
    timeMap      = e$timeMap,
    target_oni   = target_oni,
    target_month = target_month,
    n_analogs    = opt$n_analogs,
    min_sep      = opt$min_sep,
    lead_months  = opt$lead_months
  )

  if (length(init_idx) == 0L) {
    message("   SKIPPED (no analogs selected)")
    next
  }

  message("   analogs selected: ", length(init_idx),
          "  (ONI at init: range [",
          round(min(e$nino34[init_idx]), 2), ", ",
          round(max(e$nino34[init_idx]), 2), "])")

  trajs <- .lenssenExtractTrajectoriesLegacy(e$nino34, init_idx, opt$lead_months)
  res   <- .lenssenRescaleLegacy(trajs, e$nino34, obs_sd)

  scale_factors[model] <- res$scale_factor
  message("   scale_factor: ", round(res$scale_factor, 4),
          "  (model ONI SD: ", round(sd(e$nino34, na.rm = TRUE), 3), ")")

  # Store as a list of numeric vectors (one per analog); assembled into matrix below
  all_results[[model]] <- res$trajectories_adj
  all_init[[model]]    <- init_idx
}

# Step 4: Assemble into matrix and save -----------------------------------

message("\n-- Assembling cache")
non_empty   <- Filter(Negate(is.null), all_results)
models_done <- names(non_empty)

if (length(non_empty) == 0L)
  stop("No models processed. Check that zip files exist in ", opt$zenodo_dir)

# Flatten all model trajectory lists into one list of vectors
all_trajectories <- unlist(non_empty, recursive = FALSE)
n_analogs_total  <- length(all_trajectories)
n_months_traj    <- opt$lead_months + 1L   # columns per trajectory (0..lead_months)

# Build [n_analogs × n_months] matrix: each row = one trajectory
oni_matrix <- do.call(rbind, all_trajectories)
stopifnot(nrow(oni_matrix) == n_analogs_total, ncol(oni_matrix) == n_months_traj)

message("   ", length(models_done), " models  |  ",
        n_analogs_total, " analogs  |  ",
        n_months_traj, " months per trajectory (leads 0 – ", opt$lead_months, ")")

# Per-analog provenance, in the same row order as `oni_matrix`.
# `times =`, NOT `each = opt$n_analogs`: select_analogs() can return fewer than
# n_analogs for a model (short piControl run, or min_sep starving the greedy pick),
# so the per-model block lengths are not guaranteed equal.
n_per_model  <- vapply(non_empty, length, integer(1))
analog_model <- rep(models_done, times = n_per_model)
analog_init  <- unlist(all_init[models_done], use.names = FALSE)
stopifnot(length(analog_model) == n_analogs_total,
          length(analog_init)  == n_analogs_total)

# Step 4b: ENSO-state classification --------------------------------------

message("\n-- Classifying analogs into ENSO-state sets")
set_cl <- classifyEnsoAnalogSets(
  traj_mat         = oni_matrix,
  horizon_start    = horizon_start,
  peak_months      = set_peak_months,
  peak_season_year = set_peak_season_year,
  thresholds       = set_thresholds
)

if (!identical(opt$set_stat, "auto") && !identical(opt$set_stat, set_cl$stat_name))
  stop(sprintf(paste0(
    "--set_stat '%s' disagrees with --set_peak_months '%s' (which implies '%s'). ",
    "Fix one of the two, or pass --set_stat auto to accept the derived name."),
    opt$set_stat, opt$set_peak_months, set_cl$stat_name))

message("   statistic  : ", set_cl$stat_name,
        "  over ", format(set_cl$window[1L]), " .. ",
        format(set_cl$window[length(set_cl$window)]),
        "  (cols ", paste(range(set_cl$columns), collapse = "-"), ")")

# Step 5: Metadata YAML ---------------------------------------------------

meta_list <- list(
  source               = "Lenssen et al. 2024 GRL piControl (Zenodo 10045687)",
  models               = as.list(models_done),
  created              = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  horizon_start        = format(horizon_start),
  target_oni           = round(as.numeric(target_oni), 4),
  target_month         = target_month,
  n_analogs_per_model  = opt$n_analogs,
  min_sep_months       = opt$min_sep,
  lead_months          = opt$lead_months,
  obs_sd_oni           = round(obs_sd, 4),
  n_models             = length(models_done),
  n_total_analogs      = n_analogs_total,
  state_matching       = "scalar_oni_same_season",
  per_model_scale_factors = as.list(
    round(scale_factors[models_done], 4)
  ),
  # Analog-set definition: makes the ENSO-state subsets auditable straight from
  # the sidecar YAML, without loading the RDS.
  set_definition = c(
    set_cl$definition,
    list(
      per_set_model_share = lapply(
        stats::setNames(nm = sort(unique(stats::na.omit(set_cl$labels)))),
        function(s) as.list(round(
          prop.table(table(analog_model[set_cl$labels == s])), 4))
      )
    )
  )
)

meta_path <- file.path(opt$output_dir, "lenssen_nino34_analogs_meta.yml")
yaml::write_yaml(meta_list, meta_path)
message("-- Metadata written to ", meta_path)

# Assemble final cache list (format expected by sampleEnsoLenssenScenarios)
cache <- list(
  oni_trajectories = oni_matrix,
  n_months         = n_months_traj,
  n_analogs        = n_analogs_total,
  # Row-aligned with `oni_trajectories`.  `set_labels` is what
  # ensoAnalogSetRows() / sampleEnsoLenssenScenarios(analog_set=) filter on;
  # `model` / `init_index` are provenance only.
  set_labels       = set_cl$labels,
  set_stat         = set_cl$stat,
  model            = analog_model,
  init_index       = analog_init,
  meta             = meta_list
)

out_rds <- makeLenssenAnalogFile(opt$output_dir)
message("-- Saving cache to ", out_rds)
saveRDS(cache, out_rds)

# Step 6: Quick verification ----------------------------------------------

message("\n-- Verification")
message("   Matrix dim       : ", n_analogs_total, " × ", n_months_traj)
message("   NA count         : ", sum(is.na(oni_matrix)))

# Lead-0 SD is small by construction (analogs selected near target_oni)
sd_lead0  <- sd(oni_matrix[, 1L], na.rm = TRUE)
sd_global <- sd(as.vector(oni_matrix), na.rm = TRUE)
message("   SD at lead 0     : ", round(sd_lead0, 3),
        "  (narrow by design: analogs initialised near target ONI)")
message("   SD across all    : ", round(sd_global, 3),
        "  (target ~", round(obs_sd, 3), " — checks variance rescaling)")

if (abs(sd_global - obs_sd) > 0.15)
  warning("Overall trajectory SD (", round(sd_global, 3), ") differs from obs_sd (",
          round(obs_sd, 3), ") by > 0.15 — check variance rescaling")

# Analog-set composition.  Two things matter for the counterfactual sets to be
# usable: each set must be big enough to sample from without pathological reuse,
# and no set may be one model's idiosyncrasy.
message("\n-- Analog sets (", set_cl$stat_name, ")")
message("   SD of ", set_cl$stat_name, " across pool: ", round(sd(set_cl$stat), 3),
        "  (cf. lead-0 SD ", round(sd_lead0, 3),
        " — why the sets cannot be defined at lead 0)")

for (s in setdiff(ENSO_ANALOG_SETS, "all")) {
  in_set <- which(set_cl$labels == s)
  n_set  <- length(in_set)
  if (n_set == 0L) {
    warning("Analog set '", s, "' is EMPTY — the counterfactual scenario for this ",
            "set cannot be run. Widen --set_thresholds or check --set_peak_months.")
    next
  }
  v          <- set_cl$stat[in_set]
  shares     <- prop.table(table(analog_model[in_set]))
  max_share  <- max(shares)
  worst      <- names(shares)[which.max(shares)]
  message(sprintf("   %-8s n=%3d (%4.1f%%)  mean=%+.3f  sd=%.3f  range [%+.2f, %+.2f]  models=%d/%d  max share %.1f%% (%s)",
                  s, n_set, 100 * n_set / n_analogs_total,
                  mean(v), sd(v), min(v), max(v),
                  length(shares), length(models_done), 100 * max_share, worst))

  if (n_set < opt$min_set_size)
    warning("Analog set '", s, "' has only ", n_set, " trajectories. The effective ",
            "Monte-Carlo sample size for ENSO-driven quantities is bounded by this ",
            "number, not by n_replicates — intervals for this set will be coarse.")
  if (max_share > 0.40)
    warning("Analog set '", s, "' is ", round(100 * max_share, 1), "% ", worst,
            " — the set is close to a single-model idiosyncrasy rather than a ",
            "multi-model ENSO state.")
}

if (anyNA(set_cl$labels))
  warning(sum(is.na(set_cl$labels)), " analog(s) could not be classified (NA ",
          set_cl$stat_name, ").")
stopifnot(sum(!is.na(set_cl$labels)) + sum(is.na(set_cl$labels)) == n_analogs_total)

message("\nDone.")
