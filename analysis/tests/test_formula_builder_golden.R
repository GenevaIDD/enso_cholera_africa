# test_formula_builder_golden.R — characterize build_formula_terms(),
# build_size_formula_parts(), build_duration_formula_parts(), whose internals
# share private helpers.
#
# HARD gate, golden-file style: this test does not assert anything about
# *correctness* of the formula strings — it asserts current output equals a
# committed snapshot (analysis/tests/golden/formula_builder_golden.rds)
# EXACTLY: a refactor of the shared internals (per-cov/per-level RE loops,
# fixed-covariate resolution, intercept RE terms) that turns it red is not
# behaviour-preserving.
#
# Two capture passes:
#   1. Every real config under analysis/configs/** that resolves to a
#      detected model (occurrence | size | duration) is run through the
#      matching builder. Occurrence is captured twice: once with no `dat`
#      (dat = NULL, the default: no data-aware pruning) and once against the
#      synthetic AMPLE frame
#      (ample distinct levels everywhere either filter checks -> a no-op pass
#      for real configs, since none of them have a single-level column when
#      given AMPLE's levels — this exercises the dat-supplied codepath across
#      every real config without erroring). Size/duration are captured
#      against the same synthetic AMPLE frame, standing in for real
#      production data, which is out of reach here (data/ is off limits and
#      no fit is run).
#   2. A one-factor-at-a-time synthetic sweep from a "rich" baseline config
#      per track (every optional covariate/RE/knob turned on), toggling one
#      knob at a time, so each knob's individual effect on term assembly is
#      captured on its own. Size/duration/occurrence are additionally
#      captured against a synthetic SPARSE frame (every group/factor column
#      collapsed to one level), to characterize
#      filter_fixed_factor_terms_for_data() and filter_re_terms_for_data()
#      explicitly — occurrence's SPARSE capture is the synthetic analogue of
#      real configs that lose RE terms (or a fixed effect) once a
#      grouping/factor column collapses to one level.
#
# Regenerate the golden file only when a formula-string change is INTENDED: delete
# analysis/tests/golden/formula_builder_golden.rds and re-run this script,
# then diff-review the new file's git diff like any other reviewed change.
#
# Pure config + string assembly — no fit, no data/ access, no network.
#
# Run: Rscript analysis/tests/test_formula_builder_golden.R

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

GOLDEN_DIR  <- "analysis/tests/golden"
GOLDEN_PATH <- file.path(GOLDEN_DIR, "formula_builder_golden.rds")

cat("# formula builder golden characterization\n")

# --- Synthetic data frames for size/duration data-aware filters --------------
# Only grouping columns (AFRO_region/country/ADM1/ADM2/gid) and factor/character
# fixed-term columns (time_to_last_outbreak) are ever inspected by
# filter_fixed_factor_terms_for_data() / filter_re_terms_for_data() — any
# column absent from `dat` is treated as "unknown, keep" by both filters, so
# climate covariates (numeric) need not be present at all. mso_band is checked
# separately (inline, build_duration_formula_parts) for the single-occupied-
# band case. msoms1-5 / seasms1-5 are the M-spline basis columns
# .dur_ms_cols() requires to exist (not to vary) when a duration basehaz is
# "ms".

.make_dat <- function(n, region, country, adm1, adm2_gid, ttlo, mso_band) {
  d <- data.frame(
    AFRO_region            = region,
    country                = country,
    ADM1                   = adm1,
    ADM2                   = adm2_gid,
    gid                    = adm2_gid,
    time_to_last_outbreak  = factor(ttlo),
    mso_band               = mso_band,
    stringsAsFactors       = FALSE
  )
  ms_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(ms_cols) <- paste0("msoms", 1:5)
  seas_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(seas_cols) <- paste0("seasms", 1:5)
  dplyr::bind_cols(d, ms_cols, seas_cols)
}

AMPLE_DAT <- .make_dat(
  n        = 16,
  region   = rep(c("East", "West"), each = 8),
  country  = rep(c("KEN", "TZA", "NGA", "GHA"), each = 4),
  adm1     = rep(paste0("A", 1:8), each = 2),
  adm2_gid = paste0("G", 1:16),
  ttlo     = rep(c("[0,3)", "[3,6)"), length.out = 16),
  mso_band = rep(c(1, 2, 3), length.out = 16)
)

SPARSE_DAT <- .make_dat(
  n        = 6,
  region   = rep("East", 6),
  country  = rep("KEN", 6),
  adm1     = rep("A1", 6),
  adm2_gid = rep("G1", 6),
  ttlo     = rep("[0,3)", 6),
  mso_band = rep(1, 6)
)

# --- 1. Real configs ----------------------------------------------------------
cat("\n[1] capturing every real config with a detected model block\n")

resolve_config_opt <- function(path) {
  suppressWarnings(suppressMessages(
    make_options_from_config(option_list = list(), cli_argv = c("--config", path))
  ))
}

# forecast_horizons.yml / hindcast_horizons.yml are generated ENSO-horizon
# bookkeeping artefacts (produced by analysis/utils/make_forecast_configs.R),
# not 04_fit_*-style per-track configs: they have no occurrence/size/duration
# block, only generated_at/grouping/analog_sets/groups: (keyed by ENSO region,
# with nested fit-file paths and counterfactual-set metadata). Their `groups:`
# shape makes make_options_from_config() throw "subscript out of bounds" rather
# than resolve to "no model detected", so they are excluded here rather than
# silently miscounted as a resolve error.
NON_MODEL_CONFIGS <- c("manifest.yml", "forecast_horizons.yml",
                       "hindcast_horizons.yml")
# Phase-09/10 results configs and deploy manifests are pipeline-orchestration
# configs, not per-track model configs. A results config with
# results.deploy_config set additionally requires analysis/utils/deploy_core.R
# to be sourced before make_options_from_config() (a deliberate guard) which
# this golden test does not load; exclude the whole results/ and deploy/ families.
NON_MODEL_DIRS <- c("results", "deploy")
config_files <- sort(Filter(function(f)
                              !(basename(f) %in% NON_MODEL_CONFIGS) &&
                              # run-specific deploy manifests, e.g. forecast_horizons_prod_run.yml,
                              # forecast_attribution_prod_run.yml (--stage gen_attrib)
                              !grepl("^((forecast|hindcast|event_attribution)_horizons|forecast_attribution).*\\.yml$", basename(f)) &&
                              !(basename(dirname(f)) %in% NON_MODEL_DIRS) &&
                              # the synthetic run's generated batches (tools/run_synthetic_pipeline.sh)
                              !grepl("^synth_run[0-9]+_", basename(dirname(f))),
                            list.files("analysis/configs", pattern = "\\.yml$",
                                       recursive = TRUE, full.names = TRUE)))

real_capture     <- list()
real_capture_dat <- list()  # occurrence configs, dat = AMPLE_DAT
n_occ <- 0L; n_size <- 0L; n_dur <- 0L; n_mon <- 0L; n_skipped <- 0L; n_err <- 0L

for (path in config_files) {
  opt <- tryCatch(resolve_config_opt(path), error = function(e) NULL)
  if (is.null(opt)) { n_err <- n_err + 1L; next }

  model <- if (!is.null(opt$occurrence)) "occurrence"
           else if (!is.null(opt$size))     "size"
           else if (!is.null(opt$duration)) "duration"
           else if (!is.null(opt$monthly))  "monthly"
           else NULL
  if (is.null(model)) { n_skipped <- n_skipped + 1L; next }

  key <- paste0("real::", path)
  real_capture[[key]] <- tryCatch({
    suppressWarnings(suppressMessages(switch(model,
      occurrence = { n_occ  <<- n_occ  + 1L; build_formula_terms(opt, n_country = 1L, n_region = 1L) },
      size       = { n_size <<- n_size + 1L; build_size_formula_parts(opt, AMPLE_DAT) },
      duration   = { n_dur  <<- n_dur  + 1L; build_duration_formula_parts(opt, AMPLE_DAT) },
      monthly    = { n_mon  <<- n_mon  + 1L; build_monthly_formula_parts(opt, AMPLE_DAT) }
    )))
  }, error = function(e) list(.error = conditionMessage(e)))

  if (identical(model, "occurrence")) {
    real_capture_dat[[key]] <- tryCatch({
      suppressWarnings(suppressMessages(
        build_formula_terms(opt, n_country = 1L, n_region = 1L, dat = AMPLE_DAT)
      ))
    }, error = function(e) list(.error = conditionMessage(e)))
  }
}

cat(sprintf("  resolved: occurrence=%d size=%d duration=%d monthly=%d | skipped (no model)=%d | resolve errors=%d\n",
            n_occ, n_size, n_dur, n_mon, n_skipped, n_err))
# >500 in a full checkout; >300 in a code-release checkout (analysis/RELEASE),
# which versions only the deployed run's three batches (prod_run, 345 configs).
n_floor <- if (file.exists("analysis/RELEASE")) 300L else 500L
expect(n_occ + n_size + n_dur + n_mon > n_floor,
       sprintf("captured formula parts for >%d real configs (got %d)", n_floor, n_occ + n_size + n_dur + n_mon))
expect(n_err == 0L, sprintf("no config failed to resolve (got %d resolve errors)", n_err))
expect(length(real_capture_dat) == n_occ,
       sprintf("captured dat-supplied formula parts for every real occurrence config (got %d, want %d)",
               length(real_capture_dat), n_occ))

# --- 2. Synthetic OFAT sweeps --------------------------------------------------
cat("\n[2] one-factor-at-a-time synthetic sweeps from a rich per-track baseline\n")

# --- Occurrence baseline: every covariate/RE/knob this builder branches on is
# turned on, so toggling one off is visible in the captured output.
occ_base <- (function() {
  d <- get_default_options()
  opt <- list(enso = d$enso, grouping = d$grouping, contagion = d$contagion,
              occurrence = d$occurrence)
  opt$enso$lag_months <- 3L
  opt$grouping$enable <- FALSE
  opt$occurrence$backend      <- "brms"
  opt$occurrence$admin_level  <- "ADM2"
  opt$occurrence$covariates <- utils::modifyList(opt$occurrence$covariates %||% list(), list(
    time_to_last_outbreak_cat      = TRUE,
    ongoing_outbreak_lag1          = TRUE,
    ongoing_outbreak_country_lag1  = TRUE,
    ongoing_outbreak_tu_lag1       = TRUE,
    mean_precip                    = TRUE,
    dev_precip                     = TRUE,
    mean_temp_daymean               = TRUE,
    dev_temp_daymean                = TRUE,
    enso                            = TRUE,
    nbr_contagion_lag1              = TRUE
  ))
  opt$occurrence$random_effects <- utils::modifyList(opt$occurrence$random_effects %||% list(), list(
    re_region_country_adm1_intercept        = TRUE,
    re_adm2_intercept                       = FALSE,
    re_adm1_mean_precip_slope               = TRUE,
    re_adm1_dev_precip_slope                = TRUE,
    re_adm1_mean_temp_daymean_slope         = TRUE,
    re_adm1_dev_temp_daymean_slope          = TRUE,
    re_adm1_time_to_last_outbreak_slope     = TRUE,
    re_region_country_ongoing_slope         = TRUE,
    re_region_country_ongoing_country_slope = TRUE,
    re_region_country_ongoing_tu_slope      = TRUE,
    re_region_country_basehaz_slope         = TRUE,
    re_region_country_enso_lag3_slope       = TRUE,
    re_adm1_enso_lag3_slope                 = TRUE,
    re_region_country_nbr_contagion_slope   = TRUE,
    re_adm1_nbr_contagion_slope             = TRUE
  ))
  opt$occurrence$brms$basehaz_poisson  <- "seasonal"
  opt$occurrence$brms$re_uncorrelated  <- FALSE
  opt
})()

occ_variants <- list(
  baseline                    = occ_base,
  enso_off                    = { o <- occ_base; o$occurrence$covariates$enso <- FALSE; o },
  grouping_on                 = { o <- occ_base; o$grouping$enable <- TRUE; o },
  re_uncorrelated_on          = { o <- occ_base; o$occurrence$brms$re_uncorrelated <- TRUE; o },
  admin_level_adm1            = { o <- occ_base; o$occurrence$admin_level <- "ADM1"; o },
  re_adm2_intercept_on        = { o <- occ_base; o$occurrence$random_effects$re_adm2_intercept <- TRUE; o },
  basehaz_factor              = { o <- occ_base; o$occurrence$brms$basehaz_poisson <- "factor"; o },
  basehaz_cyclic              = { o <- occ_base; o$occurrence$brms$basehaz_poisson <- "cyclic"; o },
  basehaz_ms                  = { o <- occ_base; o$occurrence$brms$basehaz_poisson <- "ms"; o },
  re_ongoing_slope_off        = { o <- occ_base; o$occurrence$random_effects$re_region_country_ongoing_slope <- FALSE; o },
  re_ongoing_country_slope_off= { o <- occ_base; o$occurrence$random_effects$re_region_country_ongoing_country_slope <- FALSE; o },
  re_ongoing_tu_slope_off     = { o <- occ_base; o$occurrence$random_effects$re_region_country_ongoing_tu_slope <- FALSE; o },
  re_basehaz_slope_off        = { o <- occ_base; o$occurrence$random_effects$re_region_country_basehaz_slope <- FALSE; o },
  re_enso_shallow_off         = { o <- occ_base; o$occurrence$random_effects$re_region_country_enso_lag3_slope <- FALSE; o },
  re_enso_deep_off            = { o <- occ_base; o$occurrence$random_effects$re_adm1_enso_lag3_slope <- FALSE; o },
  nbr_contagion_off           = { o <- occ_base; o$occurrence$covariates$nbr_contagion_lag1 <- FALSE; o },
  nbr_contagion_rownorm       = { o <- occ_base; o$contagion$scale <- "rownorm"; o },
  nbr_contagion_rho200        = { o <- occ_base; o$contagion$rho_km <- 200L; o },
  re_nbr_contagion_shallow_off= { o <- occ_base; o$occurrence$random_effects$re_region_country_nbr_contagion_slope <- FALSE; o },
  re_nbr_contagion_deep_off   = { o <- occ_base; o$occurrence$random_effects$re_adm1_nbr_contagion_slope <- FALSE; o },
  all_re_off                  = { o <- occ_base
    o$occurrence$random_effects <- utils::modifyList(o$occurrence$random_effects, list(
      re_region_country_adm1_intercept = FALSE, re_adm1_mean_precip_slope = FALSE,
      re_adm1_dev_precip_slope = FALSE, re_adm1_mean_temp_daymean_slope = FALSE,
      re_adm1_dev_temp_daymean_slope = FALSE, re_adm1_time_to_last_outbreak_slope = FALSE,
      re_region_country_ongoing_slope = FALSE, re_region_country_ongoing_country_slope = FALSE,
      re_region_country_ongoing_tu_slope = FALSE, re_region_country_basehaz_slope = FALSE,
      re_region_country_enso_lag3_slope = FALSE, re_adm1_enso_lag3_slope = FALSE,
      re_region_country_nbr_contagion_slope = FALSE, re_adm1_nbr_contagion_slope = FALSE))
    o }
)

occ_sweep_capture     <- list()
occ_sweep_capture_dat <- list()  # same variants, dat = ample/sparse
for (nm in names(occ_variants)) {
  occ_sweep_capture[[nm]] <- suppressWarnings(suppressMessages(
    build_formula_terms(occ_variants[[nm]], n_country = 1L, n_region = 1L)
  ))
  for (dat_nm in c("ample", "sparse")) {
    dat <- if (dat_nm == "ample") AMPLE_DAT else SPARSE_DAT
    occ_sweep_capture_dat[[paste0(nm, "::", dat_nm)]] <- suppressWarnings(suppressMessages(
      build_formula_terms(occ_variants[[nm]], n_country = 1L, n_region = 1L, dat = dat)
    ))
  }
}
expect(length(occ_sweep_capture) == length(occ_variants),
       sprintf("captured all %d occurrence OFAT variants", length(occ_variants)))
expect(length(occ_sweep_capture_dat) == length(occ_variants) * 2L,
       sprintf("captured all %d occurrence OFAT variants x 2 dat variants", length(occ_variants)))

# --- Size baseline ---
size_base <- (function() {
  d <- get_default_options()
  opt <- list(enso = d$enso, grouping = d$grouping, size = d$size)
  opt$enso$lag_months <- 3L
  opt$grouping$enable <- FALSE
  opt$size$backend     <- "brms"
  opt$size$admin_level <- "ADM2"
  opt$size$offset_pop      <- TRUE
  opt$size$offset_meanrate <- TRUE
  opt$size$seasonal_term   <- "factor_month"
  opt$size$pop_elasticity  <- "fixed"
  opt$size$covariates <- utils::modifyList(opt$size$covariates %||% list(), list(
    time_to_last_outbreak_cat = TRUE,
    mean_precip                = TRUE,
    dev_precip                 = TRUE,
    mean_temp_daymean           = TRUE,
    dev_temp_daymean            = TRUE,
    enso                        = TRUE
  ))
  opt$size$random_effects <- utils::modifyList(opt$size$random_effects %||% list(), list(
    re_region_country_adm1_intercept   = TRUE,
    re_adm2_intercept                  = FALSE,
    re_adm1_mean_precip_slope          = TRUE,
    re_adm1_dev_precip_slope           = TRUE,
    re_adm1_mean_temp_daymean_slope    = TRUE,
    re_adm1_dev_temp_daymean_slope     = TRUE,
    re_adm1_enso_lag3_slope            = TRUE
  ))
  opt$size$brms$re_uncorrelated <- FALSE
  opt
})()

size_variants <- list(
  baseline               = size_base,
  enso_off               = { o <- size_base; o$size$covariates$enso <- FALSE; o },
  grouping_on            = { o <- size_base; o$grouping$enable <- TRUE; o },
  re_uncorrelated_on     = { o <- size_base; o$size$brms$re_uncorrelated <- TRUE; o },
  admin_level_adm1       = { o <- size_base; o$size$admin_level <- "ADM1"; o },
  re_adm2_intercept_on   = { o <- size_base; o$size$random_effects$re_adm2_intercept <- TRUE; o },
  seasonal_none          = { o <- size_base; o$size$seasonal_term <- "none"; o },
  seasonal_cyclic        = { o <- size_base; o$size$seasonal_term <- "cyclic"; o },
  seasonal_harmonic2     = { o <- size_base; o$size$seasonal_term <- "harmonic2"; o },
  pop_elasticity_free    = { o <- size_base; o$size$pop_elasticity <- "free"; o },
  pop_elasticity_free_lvl= { o <- size_base; o$size$pop_elasticity <- "free_by_level"; o },
  phi_re_region_country  = { o <- size_base; o$size$phi_re <- "region_country"; o }
)

size_sweep_capture <- list()
for (nm in names(size_variants)) {
  for (dat_nm in c("ample", "sparse")) {
    dat <- if (dat_nm == "ample") AMPLE_DAT else SPARSE_DAT
    size_sweep_capture[[paste0(nm, "::", dat_nm)]] <- suppressWarnings(suppressMessages(
      build_size_formula_parts(size_variants[[nm]], dat)
    ))
  }
}
expect(length(size_sweep_capture) == length(size_variants) * 2L,
       sprintf("captured all %d size OFAT variants x 2 dat variants", length(size_variants)))

# --- Duration baseline ---
dur_base <- (function() {
  d <- get_default_options()
  opt <- list(enso = d$enso, grouping = d$grouping, duration = d$duration)
  opt$enso$lag_months <- 3L
  opt$grouping$enable <- FALSE
  opt$duration$backend          <- "brms"
  opt$duration$admin_level       <- "ADM1"
  opt$duration$basehaz_duration  <- "factor"
  opt$duration$basehaz_poisson   <- "none"
  opt$duration$covariates <- utils::modifyList(opt$duration$covariates %||% list(), list(
    time_to_last_outbreak_cat = TRUE,
    mean_precip                = TRUE,
    dev_precip                 = TRUE,
    enso                        = TRUE
  ))
  opt$duration$random_effects <- utils::modifyList(opt$duration$random_effects %||% list(), list(
    re_region_country_adm1_intercept = TRUE,
    re_adm1_mean_precip_slope        = TRUE,
    re_adm1_dev_precip_slope         = TRUE,
    re_adm1_enso_lag3_slope          = TRUE
  ))
  opt$duration$brms$re_uncorrelated <- FALSE
  opt
})()

dur_variants <- list(
  baseline           = dur_base,
  enso_off           = { o <- dur_base; o$duration$covariates$enso <- FALSE; o },
  grouping_on        = { o <- dur_base; o$grouping$enable <- TRUE; o },
  re_uncorrelated_on = { o <- dur_base; o$duration$brms$re_uncorrelated <- TRUE; o },
  basehaz_dur_none   = { o <- dur_base; o$duration$basehaz_duration <- "none"; o },
  basehaz_dur_ms     = { o <- dur_base; o$duration$basehaz_duration <- "ms"; o },
  basehaz_seas_factor= { o <- dur_base; o$duration$basehaz_poisson <- "factor"; o },
  basehaz_seas_cyclic= { o <- dur_base; o$duration$basehaz_poisson <- "cyclic"; o },
  basehaz_seas_ms    = { o <- dur_base; o$duration$basehaz_poisson <- "ms"; o }
)

dur_sweep_capture <- list()
for (nm in names(dur_variants)) {
  for (dat_nm in c("ample", "sparse")) {
    dat <- if (dat_nm == "ample") AMPLE_DAT else SPARSE_DAT
    dur_sweep_capture[[paste0(nm, "::", dat_nm)]] <- suppressWarnings(suppressMessages(
      build_duration_formula_parts(dur_variants[[nm]], dat)
    ))
  }
}
expect(length(dur_sweep_capture) == length(dur_variants) * 2L,
       sprintf("captured all %d duration OFAT variants x 2 dat variants", length(dur_variants)))

# --- Monthly baseline ---
mon_base <- (function() {
  d <- get_default_options()
  opt <- list(enso = d$enso, grouping = d$grouping, monthly = d$monthly)
  opt$enso$lag_months <- 3L
  opt$grouping$enable <- FALSE
  opt$monthly$backend       <- "brms"
  opt$monthly$admin_level   <- "ADM2"
  opt$monthly$offset_pop    <- TRUE
  opt$monthly$seasonal_term <- "harmonic2"
  opt$monthly$phi_re        <- "region_country_adm1"
  opt$monthly$covariates <- utils::modifyList(opt$monthly$covariates %||% list(), list(
    mean_precip       = TRUE,
    dev_precip        = TRUE,
    mean_temp_daymean = TRUE,
    dev_temp_daymean  = TRUE,
    enso              = TRUE
  ))
  opt$monthly$random_effects <- utils::modifyList(opt$monthly$random_effects %||% list(), list(
    re_region_country_adm1_intercept = TRUE,
    re_adm2_intercept                = FALSE,
    re_adm1_mean_precip_slope        = TRUE,
    re_adm1_dev_precip_slope         = TRUE,
    re_adm1_mean_temp_daymean_slope  = TRUE,
    re_adm1_dev_temp_daymean_slope   = TRUE,
    re_adm1_enso_lag3_slope          = TRUE
  ))
  opt$monthly$brms$re_uncorrelated <- FALSE
  opt
})()

mon_variants <- list(
  baseline              = mon_base,
  enso_off              = { o <- mon_base; o$monthly$covariates$enso <- FALSE; o },
  grouping_on           = { o <- mon_base; o$grouping$enable <- TRUE; o },
  re_uncorrelated_on    = { o <- mon_base; o$monthly$brms$re_uncorrelated <- TRUE; o },
  admin_level_adm1      = { o <- mon_base; o$monthly$admin_level <- "ADM1"; o },
  re_adm2_intercept_on  = { o <- mon_base; o$monthly$random_effects$re_adm2_intercept <- TRUE; o },
  seasonal_none         = { o <- mon_base; o$monthly$seasonal_term <- "none"; o },
  seasonal_cyclic       = { o <- mon_base; o$monthly$seasonal_term <- "cyclic"; o },
  seasonal_factor       = { o <- mon_base; o$monthly$seasonal_term <- "factor_month"; o },
  offset_pop_off        = { o <- mon_base; o$monthly$offset_pop <- FALSE; o },
  phi_re_constant       = { o <- mon_base; o$monthly$phi_re <- "constant"; o },
  phi_re_region         = { o <- mon_base; o$monthly$phi_re <- "region"; o },
  phi_re_region_country = { o <- mon_base; o$monthly$phi_re <- "region_country"; o },
  # P4: the monthly analogue of size's pop_elasticity sweep above. `free` adds
  # the level-centred `log_pop_c` deviation term beside the RETAINED
  # offset(log(pop)); `free_by_level` interacts it with admin_level. The
  # offset-off combination is a hard error, checked separately below (it cannot
  # live in this sweep, which has no tryCatch and feeds the golden file).
  pop_elasticity_free     = { o <- mon_base; o$monthly$pop_elasticity <- "free"; o },
  pop_elasticity_free_lvl = { o <- mon_base; o$monthly$pop_elasticity <- "free_by_level"; o }
)

mon_sweep_capture <- list()
for (nm in names(mon_variants)) {
  for (dat_nm in c("ample", "sparse")) {
    dat <- if (dat_nm == "ample") AMPLE_DAT else SPARSE_DAT
    mon_sweep_capture[[paste0(nm, "::", dat_nm)]] <- suppressWarnings(suppressMessages(
      build_monthly_formula_parts(mon_variants[[nm]], dat)
    ))
  }
}
expect(length(mon_sweep_capture) == length(mon_variants) * 2L,
       sprintf("captured all %d monthly OFAT variants x 2 dat variants", length(mon_variants)))

# --- 2b. Monthly P4 (pop_elasticity) behavioural assertions -------------------
# The golden file only proves "output did not change". These assertions state
# what the new monthly P4 terms are FOR, so a future edit that silently drops
# the deviation term or the offset beside it goes red with a legible message
# rather than only as an opaque golden mismatch.
cat("\n[2b] monthly pop_elasticity (P4) term assembly\n")

.mon_parts <- function(nm) mon_sweep_capture[[paste0(nm, "::ample")]]

expect(!("log_pop_c" %in% .mon_parts("baseline")$fixed_terms),
       "mon pop_elasticity='fixed' (baseline) adds no log_pop_c term")
expect("log_pop_c" %in% .mon_parts("pop_elasticity_free")$fixed_terms,
       "mon pop_elasticity='free' adds the log_pop_c deviation term")
expect("log_pop_c:admin_level" %in% .mon_parts("pop_elasticity_free_lvl")$fixed_terms,
       "mon pop_elasticity='free_by_level' adds the log_pop_c:admin_level term")

# The offset is KEPT beside the deviation term -- that is what makes the
# coefficient read as `beta - 1` and makes coefficient 0 recover `fixed`
# exactly. Dropping it would silently reparameterise the model.
for (nm in c("pop_elasticity_free", "pop_elasticity_free_lvl")) {
  expect(any(grepl("offset(log(pop))", .mon_parts(nm)$off_terms, fixed = TRUE)),
         sprintf("mon %s KEEPS offset(log(pop)) beside the deviation term", nm))
}

# Everything except the one added fixed term must be untouched by the knob.
for (nm in c("pop_elasticity_free", "pop_elasticity_free_lvl")) {
  b <- .mon_parts("baseline"); v <- .mon_parts(nm)
  expect(identical(b$re_terms, v$re_terms) && identical(b$seas_terms, v$seas_terms) &&
           identical(b$off_terms, v$off_terms) && identical(b$slope_covs, v$slope_covs) &&
           identical(b$phi_re, v$phi_re),
         sprintf("mon %s changes ONLY fixed_terms vs baseline", nm))
  expect(identical(setdiff(v$fixed_terms, b$fixed_terms),
                   if (nm == "pop_elasticity_free") "log_pop_c" else "log_pop_c:admin_level"),
         sprintf("mon %s adds exactly one fixed term", nm))
}

# Hard-error path: the deviation parameterisation is meaningless without the
# offset it deviates FROM, so the builder must refuse rather than fit a model
# whose `log_pop_c` coefficient silently means beta instead of beta - 1.
for (elas in c("free", "free_by_level")) {
  o <- mon_base; o$monthly$pop_elasticity <- elas; o$monthly$offset_pop <- FALSE
  err <- tryCatch({ build_monthly_formula_parts(o, AMPLE_DAT); NULL },
                  error = function(e) conditionMessage(e))
  expect(!is.null(err) && grepl("offset_pop", err, fixed = TRUE),
         sprintf("mon pop_elasticity='%s' with offset_pop=false is a hard error", elas))
}

# --- 3. Compare against (or create) the committed golden file -----------------
cat("\n[3] compare against committed golden file\n")

current <- list(
  real     = real_capture,
  real_dat = real_capture_dat,
  occ      = occ_sweep_capture,
  occ_dat  = occ_sweep_capture_dat,
  size     = size_sweep_capture,
  dur      = dur_sweep_capture,
  mon      = mon_sweep_capture
)

if (!file.exists(GOLDEN_PATH)) {
  dir.create(GOLDEN_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(current, GOLDEN_PATH)
  cat("  Golden file did not exist — created:", GOLDEN_PATH, "\n")
  cat("  Review it (git diff / git add) like any other change, then re-run.\n")
  expect(TRUE, "golden file created on first run")
} else {
  golden <- readRDS(GOLDEN_PATH)

  expect(identical(names(current$real), names(golden$real)),
         sprintf("real-config key set unchanged (current=%d, golden=%d)",
                 length(current$real), length(golden$real)))
  mismatched_real <- Filter(function(k) !identical(current$real[[k]], golden$real[[k]]),
                            intersect(names(current$real), names(golden$real)))
  expect(length(mismatched_real) == 0L,
         sprintf("real-config formula parts byte-identical to golden (mismatches: %s)",
                 paste(head(mismatched_real, 5), collapse = ", ")))

  expect(identical(names(current$real_dat), names(golden$real_dat)),
         sprintf("real-config (dat-supplied) key set unchanged (current=%d, golden=%d)",
                 length(current$real_dat), length(golden$real_dat)))
  mismatched_real_dat <- Filter(function(k) !identical(current$real_dat[[k]], golden$real_dat[[k]]),
                                intersect(names(current$real_dat), names(golden$real_dat)))
  expect(length(mismatched_real_dat) == 0L,
         sprintf("real-config (dat-supplied) formula parts byte-identical to golden (mismatches: %s)",
                 paste(head(mismatched_real_dat, 5), collapse = ", ")))

  for (grp in c("occ", "occ_dat", "size", "dur", "mon")) {
    expect(identical(names(current[[grp]]), names(golden[[grp]])),
           sprintf("%s sweep key set unchanged", grp))
    mism <- Filter(function(k) !identical(current[[grp]][[k]], golden[[grp]][[k]]),
                   intersect(names(current[[grp]]), names(golden[[grp]])))
    expect(length(mism) == 0L,
           sprintf("%s sweep formula parts byte-identical to golden (mismatches: %s)",
                   grp, paste(mism, collapse = ", ")))
  }
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
