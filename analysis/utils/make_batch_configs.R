# make_batch_configs.R — batch config generator: declarative BATCH_REGISTRY ->
#   expand.grid -> per-config YAMLs + manifest.yml.
#
# Usage:
#   Rscript analysis/utils/make_batch_configs.R --list   # registry summary
#   Rscript analysis/utils/make_batch_configs.R --batch <batch name>
#   Rscript analysis/utils/make_batch_configs.R --all    # all non-frozen
#   Rscript analysis/utils/make_batch_configs.R --batch <name> --force_frozen
#
# Details:
#   Also sourced (no CLI run) by deploy_prod_run.R --stage gen_batches and for
#   interactive use of generate_batch() / discover_batch() / list_batches().

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(yaml)
  library(here)
})

source("analysis/utils/base.R") # sl_token() -- the single_lag filename token
source("analysis/utils/config.R") # get_default_options(), get_default_config()
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/priors.R") # brms prior builders
source("analysis/utils/paths.R") # nbr_covariate_name() / nbr_covariate_pattern()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_config_dir <- function() here("analysis/configs/")

get_default_config <- function() {
  yaml::read_yaml(here("analysis/config_defaults.yml"))
}

#' apply_anchored_defaults
#'
#' Opt newly generated batches in to the empirical baseline-rate anchor
#' (utils/anchor.R), so they take the residual-scale Intercept prior
#' .PRIOR_DEFAULTS[[track]]$intercept_anchored. Applied once on the resolved
#' config rather than per BATCH_REGISTRY entry because:
#'   1. the backend is known only after the grid row is merged;
#'   2. a new registry entry inherits it without a flag.
#' The inherited `prior_intercept` pin is dropped only when it equals the
#' track's un-anchored default. A different pin is a deliberate choice: it is
#' kept, and .intercept_prior_string() raises its double-counting error at fit
#' time. Existing YAMLs on disk are untouched.
#'
#' @param config resolved per-row config list
#' @param batch_model "occurrence" | "size" | "duration" | "monthly" | NULL
#' @param spec BATCH_REGISTRY entry
#' @param all_spec_keys c(names(spec$overrides), names(spec$grid))
apply_anchored_defaults <- function(config, batch_model, spec, all_spec_keys) {
  if (is.null(batch_model)) {
    return(config)
  }

  # An explicit registry setting always wins — a batch may deliberately want the
  # un-anchored model (e.g. a like-for-like comparison against an archived run).
  explicit <- paste0(batch_model, ".offset_baserate") %in%
    all_spec_keys ||
    !is.null(spec$overrides[[batch_model]]$offset_baserate)
  backend <- config[[batch_model]]$backend
  if (
    explicit ||
      is.null(backend) ||
      !backend %in% c("brms", "brms_cb", "brms_kl")
  ) {
    return(config)
  }

  config[[batch_model]]$offset_baserate <- TRUE

  inherited_pin <- .PRIOR_DEFAULTS[[batch_model]]$intercept
  for (b in c("brms", "brms_cb", "brms_kl")) {
    pin <- config[[batch_model]][[b]]$prior_intercept
    if (!is.null(pin) && identical(pin, inherited_pin)) {
      config[[batch_model]][[b]]$prior_intercept <- NULL
    }
  }
  config
}

#' strip_backend_siblings
#'
#' Drop backend sub-blocks within a model block (occurrence/size/duration/
#' monthly) that the active backend does not read, plus `dlnm:` when the
#' backend is not a brms sibling (`_cb`, or occurrence `_kl` whose lag kernel
#' reads the dlnm `lags:` basis; confirmed zero references in the plain
#' 04_fit_*_brms.R scripts), and `kernel:` unless the backend is `_kl`.
#'
#' For a sibling backend, the plain-backend sibling (e.g. `brms` for `brms_cb`)
#' is merged INTO the `_cb` block (sibling values fill any gap the `_cb`
#' block itself doesn't already set) before being dropped, so nothing set
#' only on the plain sibling (e.g. re_uncorrelated) is lost. This mirrors
#' .fill_model_block()'s own precedence (utils/config.R): block[[be]]
#' (explicit) wins over block[[plain]] (sibling).
#'
#' @param model_block  config[[batch_model]] (already backend-resolved)
#' @return             model_block with irrelevant backend sub-blocks removed
strip_backend_siblings <- function(model_block) {
  backend <- model_block$backend
  if (is.null(backend)) {
    return(model_block)
  } # ambiguous/unset — leave everything

  BACKEND_SUBBLOCKS <- c("brms", "brms_cb", "brms_kl")
  sibling <- grepl("_(cb|kl)$", backend)

  if (sibling) {
    plain <- sub("_(cb|kl)$", "", backend)
    model_block[[backend]] <- .merge1(
      model_block[[backend]],
      model_block[[plain]] %||% list()
    )
  }

  strip <- setdiff(intersect(BACKEND_SUBBLOCKS, names(model_block)), backend)
  for (s in strip) {
    model_block[[s]] <- NULL
  }

  if (!sibling) {
    model_block[["dlnm"]] <- NULL
  }
  if (!endsWith(backend, "_kl")) {
    model_block[["kernel"]] <- NULL
  }

  model_block
}

# Short tokens for varying-dimension names in config filenames.
# Pattern: {batch}_{backend}[_bh-{basehaz}][_ctry-{country}][_{abbrev}{val}...].
DIM_ABBREV <- c(
  # Occurrence nested dim keys (new schema)
  `occurrence.backend` = "bknd",
  `occurrence.brms.basehaz_poisson` = "bhp",
  `occurrence.brms.init` = "init",
  # Shared / size keys
  country = "ctry",
  baseline = "base",
  sampling = "samp",
  warmup = "wup",
  chains = "ch",
  `enso.lag_months` = "elag",
  `occurrence.dlnm.defaults.single_lag` = "sl",
  `size.dlnm.defaults.single_lag` = "sl",
  `duration.dlnm.defaults.single_lag` = "sl",
  `monthly.dlnm.defaults.single_lag` = "sl",
  `size.seasonal_term` = "sst",
  `size.admin_level` = "adml",
  `size.family` = "fam",
  `grouping.group` = "grp",
  covar_combo = "cc",
  # Monthly-track structure-compare dims
  `monthly.backend` = "bknd",
  `monthly.phi_re` = "phire",
  `monthly.prev_cases` = "prevc",
  `monthly.dlnm.defaults.re_adm1_slope` = "reslp",
  # brms_kl class sweep synthetic dims (removed by .kl_sweep_row())
  `kl_sweep.basis` = "basis",
  `kl_sweep.K` = "K",
  `kl_sweep.prior` = "prior"
)

# ---------------------------------------------------------------------------
# BATCH_REGISTRY
#
# Each entry:
#   description      label shown by --list
#   grid             list of vectors, expand.grid -> one config per row; keys
#                    are dot-path config fields ("occurrence.backend")
#   overrides        values applied to every config: 2-level dot paths
#                    ("occurrence.covariates") set a sub-block, 3-level paths
#                    ("occurrence.brms.init") a scalar; size nests in
#                    size = list(...)
#   frozen           TRUE: generate_batch() skips regeneration
#   grouping_exclude groups dropped from spec$grouping (unknown names error)
# Names: occurrence.backend first, then basehaz and country; other varying
# dims via DIM_ABBREV; the backend token is omitted when the grid has one.
# ---------------------------------------------------------------------------

# Shared DLNM argvar (exposure-response basis) definitions. Defined once so the
# project-default B-spline is not duplicated across every batch, and so a basis
# switch (e.g. bs -> ns) is a single-symbol change.
#   bs: unbounded cubic tails (degree 3, df 4 -> 1 interior knot).
#   ns: natural spline, linear tails; takes NO `degree` argument.
.ARGVAR_BS_DF4 <- list(fun = "bs", degree = 3L, df = 4L)
.ARGVAR_NS_DF3 <- list(fun = "ns", df = 3L)
.ARGVAR_NS_DF5 <- list(fun = "ns", df = 5L)
# bs df 5 as the brms_kl class sweep wrote it: no `degree` leaf (dlnm's bs
# default is cubic), so the generated configs resolve leaf for leaf as the
# hand-written analysis/configs/kl_sweep_classes_prodtu_bs5/ ones.
.ARGVAR_BS_DF5 <- list(fun = "bs", df = 5L)

BATCH_REGISTRY <- list(
  # --------------------------------------------------------------------------
  # Smoke / CI tests
  # --------------------------------------------------------------------------

  smoke_test = list(
    description = "Smoke test: brms/factor, ETH, 3 warmup / 2 sampling",
    grid = list(
      `occurrence.backend` = "brms"
    ),
    overrides = list(
      `occurrence.brms.sampling` = 2L,
      `occurrence.brms.warmup` = 3L,
      `occurrence.brms.chains` = 1L,
      `occurrence.brms.cores` = 1L,
      `occurrence.brms.basehaz_poisson` = "factor",
      country = "ETH",
      seed = 42L,
      outdir = "generated_data/smoke_test/"
    )
  ),

  smoke_test_brms_variants = list(
    description = "Smoke test: brms × 3 basehaz variants, ETH, 3 warmup / 2 sampling",
    grid = list(
      `occurrence.brms.basehaz_poisson` = c("factor", "cyclic", "ms")
    ),
    overrides = list(
      `occurrence.backend` = "brms",
      `occurrence.brms.sampling` = 2L,
      `occurrence.brms.warmup` = 3L,
      `occurrence.brms.chains` = 1L,
      `occurrence.brms.cores` = 1L,
      country = "ETH",
      seed = 42L,
      outdir = "generated_data/smoke_test/"
    )
  ),

  # --------------------------------------------------------------------------
  # Run 1 — benchmark (full Africa, default covariates)
  # --------------------------------------------------------------------------

  run1_benchmark = list(
    description = "Run 1: brms/seasonal + brms/factor, full Africa, default covariates",
    grid = list(
      `occurrence.backend` = c("brms", "brms"),
      `occurrence.brms.basehaz_poisson` = c("seasonal", "factor")
    ),
    overrides = list(
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      country = "all",
      seed = 1234L
    )
  ),

  # --------------------------------------------------------------------------
  # Run 2 — basehaz sensitivity (brms only)
  # --------------------------------------------------------------------------

  run2_basehaz_sensitivity = list(
    description = "Run 2: brms × 3 basehaz variants, full Africa",
    grid = list(
      `occurrence.brms.basehaz_poisson` = c("factor", "cyclic", "seasonal")
    ),
    overrides = list(
      `occurrence.backend` = "brms",
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      country = "all",
      seed = 1234L
    )
  ),

  # --------------------------------------------------------------------------
  # Run 4 — country sensitivity (brms/seasonal × subset of countries + all)
  # --------------------------------------------------------------------------

  run4_country_sensitivity = list(
    description = "Run 4: brms/seasonal × ETH, KEN, TZA, all",
    grid = list(
      country = c("ETH", "KEN", "TZA", "all")
    ),
    overrides = list(
      `occurrence.backend` = "brms",
      `occurrence.brms.basehaz_poisson` = "seasonal",
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      seed = 1234L
    )
  ),

  # --------------------------------------------------------------------------
  # Run 5 — covariate sensitivity (brms/seasonal, full Africa, toggle
  #          additional temperature/precipitation variables)
  # --------------------------------------------------------------------------

  run5_covariate_sensitivity = list(
    description = paste(
      "Run 5: brms/seasonal, full Africa, 4 covariate sets —",
      "default | +max_precip | +temp_daymax | +both"
    ),
    grid = list(
      country = "all" # placeholder; covariates varied via overrides below
    ),
    overrides = list(
      `occurrence.backend` = "brms",
      `occurrence.brms.basehaz_poisson` = "seasonal",
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      seed = 1234L
      # NOTE: covariate overrides require manual per-config customisation after
      # generation, or use of generate_covariate_sensitivity_batch() below.
    )
  ),

  # --------------------------------------------------------------------------
  # size_run1_benchmark — outbreak-size brms NB-log, full Africa, default RE
  # --------------------------------------------------------------------------

  size_run1_benchmark = list(
    description = "Size Run 1: brms NB-log, full Africa, default RE/covariates",
    grid = list(
      `size.seasonal_term` = "none" # placeholder single-value grid
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        brms = list(
          prior_b = "normal(0, 0.5)",
          prior_intercept = "normal(-3, 1)",
          prior_sd = "normal(0, 0.5)",
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # size_seasonal_sensitivity — vary the seasonal/month term in the size NLP
  # --------------------------------------------------------------------------

  size_seasonal_sensitivity = list(
    description = "Size: seasonal_term × {none, factor_month, cyclic, harmonic2}, full Africa, brms NB-log",
    grid = list(
      `size.seasonal_term` = c("none", "factor_month", "cyclic", "harmonic2")
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        brms = list(
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # size_enso_lag_sensitivity — outbreak-size brms with ENSO ONI lags 1-11
  # --------------------------------------------------------------------------

  size_enso_lag_sensitivity = list(
    description = "Size: ENSO ONI lags 1-11, full Africa, brms NB-log",
    grid = list(
      `enso.lag_months` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      enso = list(
        lag_months = 6L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        covariates = list(enso = TRUE),
        brms = list(
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # size_re_structure — toggle slope REs to compare RE-rich vs intercept-only
  # --------------------------------------------------------------------------

  size_re_structure = list(
    description = "Size: RE structure {intercept_only, ttlo_only, full}, full Africa",
    grid = list(
      `size.admin_level` = "ADM1" # placeholder; RE variants encoded via N copies
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        brms = list(
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
          # NOTE: RE-toggle variants require manual customization after
          # generation since they vary boolean structure inside size$random_effects.
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # enso_cb_singlelag_occ — occurrence brms_cb, single-lag DLNM on ENSO
  # One model per ENSO lag (1-11): collapse the lag dimension to a single lag
  # via dlnm::onebasis (single_lag), keeping the non-linear exposure-response
  # (argvar B-spline). ADM1 random slopes on the basis columns retained.
  # Compare lags by elpd/LOO (analysis/05_fit_postprocess_loo_lags.R).
  # --------------------------------------------------------------------------
  enso_cb_singlelag_occ = list(
    description = paste(
      "Occurrence: single-lag DLNM on ENSO, ONI lags 1-11, full Africa, brms_cb;",
      "onebasis exposure-response per lag + ADM1 basis slopes"
    ),
    grid = list(
      `occurrence.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      enso = list(
        lag_months = 5L, # unused by brms_cb single-lag (single_lag drives the column)
        source = "oni",
        cache_dir = "generated_data/"
      ),
      `occurrence.backend` = "brms_cb",
      `occurrence.admin_level` = "ADM1-ADM2",
      `occurrence.covariates` = list(
        enso = FALSE, # replaced by single-lag DLNM
        time_to_last_outbreak_cat = TRUE,
        ongoing_outbreak_lag1 = TRUE,
        mean_precip = FALSE,
        dev_precip = FALSE,
        mean_temp_daymean = FALSE,
        dev_temp_daymean = FALSE
      ),
      `occurrence.random_effects` = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_region_country_ttlo_slope = TRUE,
        re_region_country_ongoing_slope = TRUE,
        re_adm1_mean_precip_slope = FALSE,
        re_adm1_dev_precip_slope = FALSE,
        re_adm1_mean_temp_daymean_slope = FALSE,
        re_adm1_dev_temp_daymean_slope = FALSE
      ),
      `occurrence.dlnm` = list(
        enable = TRUE,
        defaults = list(
          max_lag = 11L, # only used to validate the lag column exists
          argvar = .ARGVAR_BS_DF4,
          keep_linear_term = FALSE,
          re_adm1_slope = TRUE,
          single_lag = 1L # placeholder; overwritten per grid row
        ),
        variables = list(enso = list())
      ),
      `occurrence.brms_cb` = list(
        basehaz_poisson = "seasonal",
        sampling = 1000L,
        warmup = 1000L,
        chains = 4L,
        cores = 4L,
        init = 0.2,
        re_uncorrelated = TRUE,
        threads_per_chain = 4L
      )
    )
  ),

  # --------------------------------------------------------------------------
  # size_enso_cb_singlelag — size brms_cb, single-lag DLNM on ENSO
  # Size-track analogue of enso_cb_singlelag_occ (NB-log). One model per ENSO
  # lag (1-11) via dlnm::onebasis; ADM1 random slopes on the basis columns.
  # --------------------------------------------------------------------------
  size_enso_cb_singlelag = list(
    description = paste(
      "Size: single-lag DLNM on ENSO, ONI lags 1-11, full Africa, brms_cb NB-log;",
      "onebasis exposure-response per lag + ADM1 basis slopes"
    ),
    grid = list(
      `size.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      enso = list(
        lag_months = 5L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      size = list(
        backend = "brms_cb",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1-ADM2",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        phi_re = "region_country_adm1",
        covariates = list(
          enso = FALSE, # replaced by single-lag DLNM
          time_to_last_outbreak_cat = TRUE,
          mean_precip = FALSE,
          dev_precip = FALSE,
          mean_temp_daymean = FALSE,
          dev_temp_daymean = FALSE
        ),
        random_effects = list(
          re_region_country_adm1_intercept = TRUE,
          re_adm2_intercept = TRUE,
          re_region_country_ttlo_slope = TRUE,
          re_adm1_mean_precip_slope = FALSE,
          re_adm1_dev_precip_slope = FALSE,
          re_adm1_mean_temp_daymean_slope = FALSE,
          re_adm1_dev_temp_daymean_slope = FALSE
        ),
        dlnm = list(
          enable = TRUE,
          defaults = list(
            max_lag = 11L,
            argvar = .ARGVAR_BS_DF4,
            keep_linear_term = FALSE,
            re_adm1_slope = TRUE,
            single_lag = 1L # placeholder; overwritten per grid row
          ),
          variables = list(enso = list())
        ),
        brms_cb = list(
          prior_b = "normal(0, 0.5)",
          prior_intercept = "normal(-3, 1)",
          prior_sd = "normal(0, 0.5)",
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1,
          re_uncorrelated = TRUE,
          threads_per_chain = 2L
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # africa_enso_singlelag_occ — occurrence brms_cb, single-lag DLNM on ENSO,
  # all Africa, hierarchical CB slopes (region/country/ADM1).
  # Production-scale successor to enso_cb_singlelag_occ, matching the full
  # Africa occurrence spec (africa_brms_cb_enso.yml): adds ongoing_country and
  # ongoing_TU covariates + their nested slope REs. One model per ENSO lag
  # (1-11). Compare lags by elpd/LOO (analysis/05_fit_postprocess_loo_lags.R).
  # --------------------------------------------------------------------------
  africa_enso_singlelag_occ = list(
    description = paste(
      "Occurrence: single-lag DLNM on ENSO, ONI lags 1-11, all Africa, brms_cb;",
      "hierarchical CB slopes (region/country/ADM1); full ongoing-outbreak RE spec"
    ),
    grid = list(
      `occurrence.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 5L, # unused by brms_cb single-lag
        source = "oni",
        cache_dir = "generated_data/"
      ),
      `occurrence.backend` = "brms_cb",
      `occurrence.admin_level` = "ADM1-ADM2",
      `occurrence.ongoing_outbreak_re` = TRUE,
      `occurrence.covariates` = list(
        enso = FALSE,
        time_to_last_outbreak_cat = TRUE,
        ongoing_outbreak_lag1 = TRUE,
        ongoing_outbreak_country_lag1 = TRUE,
        ongoing_outbreak_tu_lag1 = TRUE,
        mean_precip = FALSE,
        dev_precip = FALSE,
        mean_temp_daymean = FALSE,
        dev_temp_daymean = FALSE
      ),
      `occurrence.random_effects` = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_region_country_ttlo_slope = TRUE,
        re_region_country_ongoing_slope = TRUE,
        re_region_country_ongoing_country_slope = TRUE,
        re_region_country_ongoing_tu_slope = TRUE,
        re_adm1_mean_precip_slope = FALSE,
        re_adm1_dev_precip_slope = FALSE,
        re_adm1_mean_temp_daymean_slope = FALSE,
        re_adm1_dev_temp_daymean_slope = FALSE,
        re_adm1_enso_lag5_slope = FALSE # replaced by dlnm re_adm1_slope
      ),
      `occurrence.dlnm` = list(
        enable = TRUE,
        defaults = list(
          max_lag = 11L,
          argvar = .ARGVAR_BS_DF4,
          keep_linear_term = FALSE,
          re_adm1_slope = TRUE,
          single_lag = 1L # placeholder; overwritten per grid row
        ),
        variables = list(enso = list())
      ),
      `occurrence.brms` = list(
        re_uncorrelated = TRUE # controls regular (non-CB) slope RE separator
      ),
      `occurrence.brms_cb` = list(
        basehaz_poisson = "seasonal",
        sampling = 100L,
        warmup = 100L,
        chains = 4L,
        cores = 4L,
        init = 0.2,
        re_uncorrelated = TRUE,
        threads_per_chain = 4L
      )
    )
  ),

  # --------------------------------------------------------------------------
  # africa_enso_singlelag_siz — size brms_cb, single-lag DLNM on ENSO,
  # all Africa, hierarchical CB slopes (region/country/ADM1).
  # Production-scale successor to size_enso_cb_singlelag, matching the full
  # Africa size spec (africa_size_brms_cb_enso.yml). One model per ENSO lag
  # (1-11).
  # --------------------------------------------------------------------------
  africa_enso_singlelag_siz = list(
    description = paste(
      "Size: single-lag DLNM on ENSO, ONI lags 1-11, all Africa, brms_cb NB-log;",
      "hierarchical CB slopes (region/country/ADM1)"
    ),
    grid = list(
      `size.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 5L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      size = list(
        backend = "brms_cb",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1-ADM2",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        phi_re = "region_country_adm1",
        covariates = list(
          enso = FALSE,
          time_to_last_outbreak_cat = TRUE,
          mean_precip = FALSE,
          dev_precip = FALSE,
          mean_temp_daymean = FALSE,
          dev_temp_daymean = FALSE
        ),
        random_effects = list(
          re_region_country_adm1_intercept = TRUE,
          re_adm2_intercept = TRUE,
          re_region_country_ttlo_slope = TRUE,
          re_adm1_mean_precip_slope = FALSE,
          re_adm1_dev_precip_slope = FALSE,
          re_adm1_mean_temp_daymean_slope = FALSE,
          re_adm1_dev_temp_daymean_slope = FALSE,
          re_adm1_enso_lag10_slope = FALSE # replaced by dlnm re_adm1_slope
        ),
        dlnm = list(
          enable = TRUE,
          defaults = list(
            max_lag = 11L,
            argvar = .ARGVAR_BS_DF4,
            keep_linear_term = FALSE,
            re_adm1_slope = TRUE,
            single_lag = 1L # placeholder; overwritten per grid row
          ),
          variables = list(enso = list())
        ),
        brms_cb = list(
          prior_b = "normal(0, 0.5)",
          prior_intercept = "normal(-3, 1)",
          prior_sd = "normal(0, 0.5)",
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1,
          re_uncorrelated = TRUE,
          threads_per_chain = 2L
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # Country-grouping example batches
  #
  # `grouping = "<registry>"` expands one config per group in
  # analysis/groupings/<registry>.csv. Each config sets grouping.enable=TRUE,
  # grouping.source, grouping.group; fit-time resolve_group_countries() derives
  # `country` from the registry. RE nesting drops AFRO_region (country/ADM1 only).
  # --------------------------------------------------------------------------

  occ_grouped_afro_subregion = list(
    description = "Occurrence: brms/seasonal per AFRO subregion (east/central/west), country/ADM1 nesting",
    grouping = "afro_subregion",
    grid = list(
      `occurrence.backend` = "brms" # placeholder single-value grid
    ),
    overrides = list(
      `occurrence.backend` = "brms",
      `occurrence.brms.basehaz_poisson` = "seasonal",
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      seed = 1234L
    )
  ),

  size_grouped_afro_subregion = list(
    description = "Size: brms NB-log per AFRO subregion (east/central/west), country/ADM1 nesting",
    grouping = "afro_subregion",
    grid = list(
      `size.seasonal_term` = "none" # placeholder single-value grid
    ),
    overrides = list(
      seed = 1234L,
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        brms = list(
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  # HMM cluster groupings (3 groups: hmmCE / hmmS / hmmW)
  occ_grouped_hmm_clusters = list(
    description = "Occurrence: brms/seasonal per HMM cluster (hmmCE/hmmS/hmmW), country/ADM1 nesting",
    grouping = "hmm_clusters",
    grid = list(
      `occurrence.backend` = "brms" # placeholder single-value grid
    ),
    overrides = list(
      `occurrence.brms.basehaz_poisson` = "seasonal",
      `occurrence.brms.sampling` = 1000L,
      `occurrence.brms.warmup` = 1000L,
      `occurrence.brms.chains` = 4L,
      `occurrence.brms.cores` = 4L,
      seed = 1234L
    )
  ),

  size_grouped_hmm_clusters = list(
    description = "Size: brms NB-log per HMM cluster (hmmCE/hmmS/hmmW), country/ADM1 nesting",
    grouping = "hmm_clusters",
    grid = list(
      `size.seasonal_term` = "none" # placeholder single-value grid
    ),
    overrides = list(
      seed = 1234L,
      size = list(
        backend = "brms",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        brms = list(
          sampling = 1000L,
          warmup = 1000L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # occ_grouped_enso_singlelag — occurrence brms_cb single-lag DLNM × 5 ENSO
  # teleconnection groups (enso_grouping.csv: enso_east/south/west/horn/central).
  # First-pass iteration budget: 300 warmup / 200 sampling.
  # 5 groups × 11 lags = 55 configs.
  # RE spec mirrors africa_enso_singlelag_occ (full ongoing-outbreak RE structure).
  # phi_re is not applicable here (occurrence model).
  # --------------------------------------------------------------------------
  occ_grouped_enso_singlelag = list(
    description = paste(
      "Occurrence: single-lag DLNM on ENSO, ONI lags 1-11, grouped by ENSO",
      "teleconnection (enso_east/south/west/horn/central), brms_cb;",
      "full ongoing-outbreak RE spec; 300 warmup / 200 sampling (first pass)"
    ),
    grouping = "enso_grouping",
    grid = list(
      `occurrence.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 5L, # unused by brms_cb single-lag
        source = "oni",
        cache_dir = "generated_data/"
      ),
      `occurrence.backend` = "brms_cb",
      `occurrence.admin_level` = "ADM1-ADM2",
      `occurrence.ongoing_outbreak_re` = TRUE,
      `occurrence.covariates` = list(
        enso = FALSE, # replaced by single-lag DLNM
        time_to_last_outbreak_cat = TRUE,
        ongoing_outbreak_lag1 = TRUE,
        ongoing_outbreak_country_lag1 = TRUE,
        ongoing_outbreak_tu_lag1 = TRUE,
        mean_precip = FALSE,
        dev_precip = FALSE,
        mean_temp_daymean = FALSE,
        dev_temp_daymean = FALSE
      ),
      `occurrence.random_effects` = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_region_country_ttlo_slope = TRUE,
        re_region_country_ongoing_slope = TRUE,
        re_region_country_ongoing_country_slope = TRUE,
        re_region_country_ongoing_tu_slope = TRUE,
        re_adm1_mean_precip_slope = FALSE,
        re_adm1_dev_precip_slope = FALSE,
        re_adm1_mean_temp_daymean_slope = FALSE,
        re_adm1_dev_temp_daymean_slope = FALSE,
        re_adm1_enso_lag5_slope = FALSE # replaced by dlnm re_adm1_slope
      ),
      `occurrence.dlnm` = list(
        enable = TRUE,
        defaults = list(
          max_lag = 11L,
          argvar = .ARGVAR_BS_DF4,
          keep_linear_term = FALSE,
          re_adm1_slope = TRUE,
          single_lag = 1L # placeholder; overwritten per grid row
        ),
        variables = list(enso = list())
      ),
      `occurrence.brms` = list(
        re_uncorrelated = TRUE
      ),
      `occurrence.brms_cb` = list(
        basehaz_poisson = "seasonal",
        sampling = 200L,
        warmup = 300L,
        chains = 4L,
        cores = 4L,
        init = 0.2,
        re_uncorrelated = TRUE,
        threads_per_chain = 4L
      )
    )
  ),

  # --------------------------------------------------------------------------
  # mon_smoketest — smallest batch-registry-generated config exercising the
  # monthly track (brms_cb backend, ENSO-grouped, one group only). Sibling in
  # spirit to the hand-written analysis/configs/mon_smoketest/mon_smoke.yml
  # (kept side-by-side; NOT overwritten -- this entry writes a differently
  # named file, mon_smoketest.yml, into the same directory since it has no
  # varying grid dimension). Uses grouping.group = enso_east directly (not
  # the auto-expanding spec$grouping registry field) since the smoke test needs
  # exactly one config / one group.
  # --------------------------------------------------------------------------
  mon_smoketest = list(
    description = "Smoke test: monthly track, brms_cb backend, ENSO-grouped (enso_east only), 200 warmup / 200 sampling",
    grid = list(
      `monthly.backend` = "brms_cb"
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 4L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      `grouping.enable` = TRUE,
      `grouping.source` = "enso_grouping",
      `grouping.group` = "enso_east",
      monthly = list(
        backend = "brms_cb",
        brms_cb = list(
          sampling = 200L,
          warmup = 200L,
          chains = 2L,
          cores = 2L,
          refresh = NULL,
          init = "random",
          re_uncorrelated = FALSE,
          adapt_delta = 0.95,
          max_treedepth = 12L
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # mon_structure_compare — structural sensitivity sweep for the monthly
  # track: {brms, brms_cb} backend x {phi_re: constant, region_country_adm1}
  # x {prev_cases: FALSE, TRUE} x {re_adm1_slope: FALSE, TRUE} (16 configs,
  # full Africa, no ENSO grouping). re_adm1_slope (dlnm.defaults.re_adm1_slope,
  # hierarchical crossbasis slopes) is only consumed by 04_fit_mon_brms_cb.R
  # -- the plain brms backend (04_fit_mon_brms.R) does not read the dlnm:
  # block at all, so the re_adm1_slope dimension is inert (duplicate configs)
  # for the four backend=brms rows; kept in the cross product anyway rather
  # than special-cased out, for a uniform 2x2x2x2 grid.
  # --------------------------------------------------------------------------
  mon_structure_compare = list(
    description = paste(
      "Monthly: structure sweep -- {brms, brms_cb} x",
      "{phi_re: constant, region_country_adm1} x {prev_cases: F, T} x",
      "{re_adm1_slope: F, T} (re_adm1_slope inert under brms); full Africa,",
      "300 warmup / 300 sampling"
    ),
    grid = list(
      `monthly.backend` = c("brms", "brms_cb"),
      `monthly.phi_re` = c("constant", "region_country_adm1"),
      `monthly.prev_cases` = c(FALSE, TRUE),
      `monthly.dlnm.defaults.re_adm1_slope` = c(FALSE, TRUE)
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 4L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      monthly = list(
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1-ADM2",
        offset_pop = TRUE,
        seasonal_term = "harmonic2",
        covariates = list(),
        random_effects = list(
          re_region_country_adm1_intercept = TRUE,
          re_adm2_intercept = TRUE
        ),
        dlnm = list(
          enable = TRUE,
          defaults = list(
            max_lag = 3L,
            argvar = .ARGVAR_BS_DF4,
            arglag = list(fun = "ns", df = 3L),
            keep_linear_term = FALSE,
            single_lag = 0L
          ),
          variables = list(
            mso = list(single_lag = 0L, argvar = .ARGVAR_BS_DF4),
            enso = list(
              single_lag = 4L,
              argvar = list(fun = "bs", degree = 3L, df = 8L)
            )
          )
        ),
        brms = list(
          sampling = 300L,
          warmup = 300L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1
        )
      )
    )
  ),

  size_grouped_enso_singlelag = list(
    description = paste(
      "Size: single-lag DLNM on ENSO, ONI lags 1-11, grouped by ENSO",
      "teleconnection (enso_east/south/west/horn/central), brms_cb NB-log;",
      "phi_re = region_country_adm1; 300 warmup / 200 sampling (first pass)"
    ),
    grouping = "enso_grouping",
    grid = list(
      `size.dlnm.defaults.single_lag` = 1:11
    ),
    overrides = list(
      country = "all",
      seed = 1234L,
      outdir = "generated_data/",
      TL = "1970-01-01",
      TR = "2024-12-31",
      enso = list(
        lag_months = 5L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      size = list(
        backend = "brms_cb",
        family = "negbinomial",
        link = "log",
        admin_level = "ADM1-ADM2",
        offset_pop = TRUE,
        offset_meanrate = TRUE,
        seasonal_term = "none",
        phi_re = "region_country_adm1",
        covariates = list(
          enso = FALSE, # replaced by single-lag DLNM
          time_to_last_outbreak_cat = TRUE,
          mean_precip = FALSE,
          dev_precip = FALSE,
          mean_temp_daymean = FALSE,
          dev_temp_daymean = FALSE
        ),
        random_effects = list(
          re_region_country_adm1_intercept = TRUE,
          re_adm2_intercept = TRUE,
          re_region_country_ttlo_slope = TRUE,
          re_adm1_mean_precip_slope = FALSE,
          re_adm1_dev_precip_slope = FALSE,
          re_adm1_mean_temp_daymean_slope = FALSE,
          re_adm1_dev_temp_daymean_slope = FALSE,
          re_adm1_enso_lag10_slope = FALSE # replaced by dlnm re_adm1_slope
        ),
        dlnm = list(
          enable = TRUE,
          defaults = list(
            max_lag = 11L,
            argvar = .ARGVAR_BS_DF4,
            keep_linear_term = FALSE,
            re_adm1_slope = TRUE,
            single_lag = 1L # placeholder; overwritten per grid row
          ),
          variables = list(enso = list())
        ),
        brms_cb = list(
          prior_b = "normal(0, 0.5)",
          prior_intercept = "normal(-3, 1)",
          prior_sd = "normal(0, 0.5)",
          sampling = 200L,
          warmup = 300L,
          chains = 4L,
          cores = 4L,
          adapt_delta = 0.95,
          max_treedepth = 12L,
          init = 0.1,
          re_uncorrelated = TRUE,
          threads_per_chain = 2L
        )
      )
    )
  )
)


# ===========================================================================
# prod_run: the deployed model (configs/deploy/prod_run.yml)
# ===========================================================================
# Three fit batches (occurrence, duration, monthly; ONI lag -11..11 x the fit
# groups of enso_grouping_prod_run) and the brms_kl class sweep that derives that
# grouping (deploy_prod_run.R --stage gen_grouping, cell pinned in the deploy).
# The specs are written out in full: a generated config is config_defaults.yml
# plus these overrides plus its grid row.

#' Drop the legacy unscoped `prior_sd` from the ACTIVE backend block of
#' a generated model block when BOTH tier keys (prior_sd_region_country,
#' prior_sd_adm1) are set there -- the tiers win and `prior_sd` is never read
#' (priors.R warns about it). Applied by generate_batch() only under
#' `spec$drop_superseded_prior_sd = TRUE` (the prod_run tracks); a block with
#' one tier or none keeps `prior_sd`.
.drop_superseded_prior_sd <- function(model_block) {
  be <- model_block[["backend"]]
  if (is.null(be) || !is.list(model_block[[be]])) {
    return(model_block)
  }
  blk <- model_block[[be]]
  if (
    !is.null(blk[["prior_sd_region_country"]]) &&
      !is.null(blk[["prior_sd_adm1"]]) &&
      "prior_sd" %in% names(blk)
  ) {
    blk[["prior_sd"]] <- NULL
    model_block[[be]] <- blk
  }
  model_block
}

#' Per-row transform of the class sweep: sets the leaves of the row's cell
#' (basis x K, from `cells`) and removes the synthetic kl_sweep.* grid keys.
.prod_run_kl_cell_row <- function(cells) {
  function(cfg) {
    key <- paste0(cfg$kl_sweep$basis, "_K", cfg$kl_sweep$K)
    cell <- cells[[key]]
    if (is.null(cell)) {
      stop("prod_run_kl_sweep_classes: unknown sweep cell ", key, call. = FALSE)
    }
    cfg$kl_sweep <- NULL
    for (d in cell) {
      if (isTRUE(d$remove)) cfg[[d$path]] <- NULL else cfg[[d$path]] <- d$value
    }
    cfg
  }
}

BATCH_REGISTRY$prod_run_occ_grouped_enso_singlelag <- list(
  description = "Occurrence (onset hazard): Poisson GLMM per unit-month at risk, seasonal cyclic spline, ENSO crossbasis (ONI at one group lag, ns df 4) with country and country:ADM1 slopes, recorded contagion and outbreak history; ONI lags -11..11 x the fit groups.",
  grouping = "enso_grouping_prod_run",
  grid = list(occurrence.dlnm.defaults.single_lag = -11:11),
  overrides = list(
    country = "all",
    seed = 1234L,
    outdir = "generated_data/",
    TL = "1970-01-01",
    TR = "2024-12-31",
    enso = list(source = "oni", cache_dir = "generated_data/"),
    occurrence.backend = "brms_cb",
    occurrence.admin_level = "ADM1-ADM2",
    occurrence.ongoing_outbreak_re = TRUE,
    occurrence.covariates = list(
      enso = FALSE,
      time_to_last_outbreak_cat = TRUE,
      ongoing_outbreak_lag1 = TRUE,
      ongoing_outbreak_country_lag1 = FALSE,
      ongoing_outbreak_tu_lag1 = FALSE,
      mean_precip = FALSE,
      dev_precip = FALSE,
      mean_temp_daymean = FALSE,
      dev_temp_daymean = FALSE,
      nbr_contagion_lag1 = TRUE
    ),
    occurrence.random_effects = list(
      re_region_country_adm1_intercept = TRUE,
      re_adm2_intercept = TRUE,
      re_region_country_ttlo_slope = TRUE,
      re_region_country_ongoing_slope = TRUE,
      re_region_country_ongoing_country_slope = FALSE,
      re_region_country_ongoing_tu_slope = FALSE,
      re_adm1_mean_precip_slope = FALSE,
      re_adm1_dev_precip_slope = FALSE,
      re_adm1_mean_temp_daymean_slope = FALSE,
      re_adm1_dev_temp_daymean_slope = FALSE,
      re_region_country_nbr_contagion_slope = TRUE,
      re_adm1_nbr_contagion_slope = TRUE
    ),
    occurrence.dlnm = list(
      enable = TRUE,
      defaults = list(
        max_lag = 11L,
        argvar = list(
          fun = "ns",
          df = 4L
        ),
        keep_linear_term = FALSE,
        re_adm1_slope = TRUE,
        single_lag = 1L
      ),
      variables = list(enso = list())
    ),
    occurrence.brms = list(re_uncorrelated = TRUE),
    occurrence.brms_cb = list(
      basehaz_poisson = "seasonal",
      sampling = 1000L,
      warmup = 350L,
      chains = 4L,
      cores = 4L,
      init = 0.2,
      re_uncorrelated = TRUE,
      threads_per_chain = 4L,
      prior_sd_region_country = "normal(0, 0.35)",
      prior_sd_adm1 = "normal(0, 0.35)",
      prior_sds = "normal(0, 0.35)"
    ),
    contagion = list(
      rho_km = 500L,
      scale = "rownorm",
      partition = "registry",
      partition_source = "hmm_clusters"
    ),
    episodes = list(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L,
      national_zero_rule = "r3",
      at_risk_start = "after_episode",
      adm1_overlap = "censor_orphan"
    ),
    occurrence.canonical_admin_level = TRUE,
    population = list(
      gid_pop_path = "generated_data/gid_population_2015-2030.qs"
    )
  ),
  drop_superseded_prior_sd = TRUE
)

BATCH_REGISTRY$prod_run_dur_grouped_enso_singlelag <- list(
  description = "Duration (cessation hazard): Poisson GLMM per active episode-month, months-since-onset bands (prior normal(0, 1.5)) and harmonic2 cessation seasonality, ENSO crossbasis with country and country:ADM1 slopes; ONI lags -11..11 x the fit groups.",
  grouping = "enso_grouping_prod_run",
  grid = list(duration.dlnm.defaults.single_lag = -11:11),
  overrides = list(
    country = "all",
    seed = 1234L,
    outdir = "generated_data/",
    TL = "1970-01-01",
    TR = "2024-12-31",
    enso = list(source = "oni", cache_dir = "generated_data/"),
    duration.backend = "brms_cb",
    duration.admin_level = "ADM1-ADM2",
    duration.basehaz_duration = "factor",
    duration.basehaz_poisson = "none",
    duration.ttlo_reset = "onset",
    duration.covariates = list(enso = FALSE),
    duration.random_effects = list(
      re_region_country_adm1_intercept = TRUE,
      re_adm1_mean_precip_slope = FALSE,
      re_adm1_dev_precip_slope = FALSE,
      re_adm1_mean_temp_daymean_slope = FALSE,
      re_adm1_dev_temp_daymean_slope = FALSE,
      re_adm2_intercept = TRUE,
      re_adm2_leaf_scope = "adm2_only"
    ),
    duration.dlnm = list(
      enable = TRUE,
      defaults = list(
        argvar = list(fun = "ns", df = 4L),
        keep_linear_term = FALSE,
        re_adm1_slope = TRUE,
        single_lag = 1L
      ),
      variables = list(enso = list())
    ),
    duration.brms_cb = list(
      sampling = 1000L,
      warmup = 350L,
      chains = 4L,
      cores = 4L,
      refresh = NULL,
      init = "random",
      re_uncorrelated = TRUE,
      adapt_delta = 0.95,
      max_treedepth = 12L,
      threads_per_chain = 4L,
      prior_sd_region_country = "normal(0, 0.35)",
      prior_sd_adm1 = "normal(0, 0.35)",
      prior_mso_band = "normal(0, 1.5)"
    ),
    episodes = list(key = "gid_season"),
    duration.canonical_admin_level = TRUE,
    population = list(
      gid_pop_path = "generated_data/gid_population_2015-2030.qs"
    ),
    duration.seasonal_term = "harmonic2"
  ),
  drop_superseded_prior_sd = TRUE
)

BATCH_REGISTRY$prod_run_mon_grouped_enso_singlelag <- list(
  description = "Monthly case counts: negative-binomial GLMM per active episode-month, population offset and elasticity, harmonic2 seasonality, months-since-onset spline, ENSO crossbasis with country and country:ADM1 slopes; ONI lags -11..11 x the fit groups.",
  grouping = "enso_grouping_prod_run",
  grid = list(monthly.dlnm.defaults.single_lag = -11:11),
  overrides = list(
    country = "all",
    seed = 1234L,
    outdir = "generated_data/",
    TL = "1970-01-01",
    TR = "2024-12-31",
    enso = list(source = "oni", cache_dir = "generated_data/"),
    monthly = list(
      backend = "brms_cb",
      family = "negbinomial",
      link = "log",
      admin_level = "ADM1-ADM2",
      offset_pop = TRUE,
      seasonal_term = "harmonic2",
      phi_re = "region_country_adm1",
      prev_cases = FALSE,
      covariates = list(),
      random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_adm2_leaf_scope = "adm2_only"
      ),
      dlnm = list(
        enable = TRUE,
        defaults = list(
          max_lag = 11L,
          argvar = list(fun = "ns", df = 4L),
          keep_linear_term = FALSE,
          re_adm1_slope = TRUE,
          single_lag = 1L
        ),
        variables = list(
          mso = list(
            single_lag = 0L,
            argvar = list(fun = "bs", degree = 3L, df = 4L),
            re_adm1_slope = TRUE,
            re_slope_levels = c("region", "country")
          ),
          enso = list()
        )
      ),
      brms_cb = list(
        sampling = 1000L,
        warmup = 350L,
        chains = 4L,
        cores = 4L,
        refresh = NULL,
        init = "random",
        re_uncorrelated = TRUE,
        adapt_delta = 0.95,
        max_treedepth = 12L,
        threads_per_chain = 4L,
        prior_sd_region_country = "normal(0, 0.35)",
        prior_sd_adm1 = "normal(0, 0.35)"
      ),
      estimand = list(
        denominator = "gid",
        rate_denominator = "gid",
        underreporting_multiplier = 1
      ),
      pop_elasticity = "free",
      canonical_admin_level = TRUE
    ),
    population = list(
      gid_pop_path = "generated_data/gid_population_2015-2030.qs"
    ),
    episodes = list(
      key = "gid_season"
    )
  ),
  drop_superseded_prior_sd = TRUE
)

.PROD_RUN_KL_CELLS <- list(
  ns4_K2 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 2L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k2",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_west = list("NER", "LBR", "TGO")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(enso_east = 0L, enso_west = -6L),
      remove = FALSE
    )
  ),
  ns4_K3 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 3L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k3",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list("MWI", "ZWE", "ZMB"),
          enso_west = list(
            "NER",
            "LBR",
            "TGO"
          )
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(enso_east = 0L, enso_south = 4.3, enso_west = -6L),
      remove = FALSE
    )
  ),
  ns4_K4 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k4",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list("MWI", "ZWE", "ZMB"),
          enso_west = list("NER", "LBR", "TGO"),
          enso_horn = list("ETH", "SDN")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(
        enso_east = 0L,
        enso_south = 4.3,
        enso_west = -6L,
        enso_horn = 9.5
      ),
      remove = FALSE
    )
  ),
  ns4_K5 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_nga",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list(
            "MWI",
            "ZWE",
            "ZMB"
          ),
          enso_west = list("NER", "LBR", "TGO"),
          enso_nga = list("NGA"),
          enso_horn = list("ETH", "SDN")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(
        enso_east = 0L,
        enso_south = 4.3,
        enso_west = -6L,
        enso_nga = 10L,
        enso_horn = 9.5
      ),
      remove = FALSE
    )
  ),
  ns5_K2 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 2L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k2",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_west = list("NER", "LBR", "TGO")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(
        enso_east = 0L,
        enso_west = -6L
      ),
      remove = FALSE
    )
  ),
  ns5_K3 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 3L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k3",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list(
            "MWI",
            "ZWE",
            "ZMB"
          ),
          enso_west = list("NER", "LBR", "TGO")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(enso_east = 0L, enso_south = 4.3, enso_west = -6L),
      remove = FALSE
    )
  ),
  ns5_K4 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 4L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_k4",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list("MWI", "ZWE", "ZMB"),
          enso_west = list("NER", "LBR", "TGO"),
          enso_horn = list("ETH", "SDN")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(
        enso_east = 0L,
        enso_south = 4.3,
        enso_west = -6L,
        enso_horn = 9.5
      ),
      remove = FALSE
    )
  ),
  ns5_K5 = list(
    list(
      path = "outdir",
      value = "generated_data/prod_run_kl_sweep_classes/",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "dlnm", "defaults", "argvar", "df"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "classes"),
      value = 5L,
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "prior_grouping", "source"),
      value = "enso_grouping_nga",
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "seeds"),
      value = list(
        iso3 = list(
          enso_east = list("TZA", "UGA", "BDI"),
          enso_south = list(
            "MWI",
            "ZWE",
            "ZMB"
          ),
          enso_west = list("NER", "LBR", "TGO"),
          enso_nga = list("NGA"),
          enso_horn = list("ETH", "SDN")
        )
      ),
      remove = FALSE
    ),
    list(
      path = c("occurrence", "kernel", "init_mu"),
      value = list(
        enso_east = 0L,
        enso_south = 4.3,
        enso_west = -6L,
        enso_nga = 10L,
        enso_horn = 9.5
      ),
      remove = FALSE
    )
  )
)

BATCH_REGISTRY$prod_run_kl_sweep_classes <- c(
  list(
    description = "brms_kl latent-class sweep that derives the fit grouping: K 2-5 x ENSO basis {ns df 4, 5}, Gaussian lag kernel seeded by core countries, ADM1 ENSO slopes, 400 warm-up + 300 sampling.",
    grid = list(kl_sweep.basis = c("ns4", "ns5"), kl_sweep.K = 2:5),
    overrides = list(
      baseline = "1990_2019",
      country = "all",
      TL = "1970-01-01",
      TR = "2024-12-31",
      seed = 1234L,
      contagion = list(
        rho_km = 500L,
        scale = "rownorm",
        partition = "registry",
        partition_source = "hmm_clusters"
      ),
      enso = list(
        lag_months = 5L,
        source = "oni",
        cache_dir = "generated_data/"
      ),
      population = list(
        gid_pop_path = "generated_data/gid_population_2015-2030.qs"
      ),
      episodes = list(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L,
        national_zero_rule = "r3",
        at_risk_start = "after_episode",
        adm1_overlap = "censor_orphan"
      ),
      occurrence.backend = "brms_kl",
      occurrence.admin_level = "ADM1-ADM2",
      occurrence.ongoing_outbreak_re = TRUE,
      occurrence.offset_baserate = TRUE,
      occurrence.canonical_admin_level = TRUE,
      occurrence.covariates = list(
        enso = FALSE,
        time_to_last_outbreak_cat = TRUE,
        ongoing_outbreak_lag1 = TRUE,
        ongoing_outbreak_country_lag1 = FALSE,
        ongoing_outbreak_tu_lag1 = FALSE,
        mean_precip = FALSE,
        dev_precip = FALSE,
        mean_temp_daymean = FALSE,
        dev_temp_daymean = FALSE,
        nbr_contagion_lag1 = TRUE
      ),
      occurrence.random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_region_country_ttlo_slope = TRUE,
        re_region_country_ongoing_slope = TRUE,
        re_region_country_ongoing_country_slope = FALSE,
        re_region_country_ongoing_tu_slope = FALSE,
        re_adm1_mean_precip_slope = FALSE,
        re_adm1_dev_precip_slope = FALSE,
        re_adm1_mean_temp_daymean_slope = FALSE,
        re_adm1_dev_temp_daymean_slope = FALSE,
        re_adm1_enso_lag5_slope = FALSE,
        re_region_country_nbr_contagion_slope = TRUE,
        re_adm1_nbr_contagion_slope = TRUE
      ),
      occurrence.dlnm = list(
        enable = TRUE,
        defaults = list(
          argvar = list(fun = "ns", df = 3L),
          keep_linear_term = FALSE,
          re_adm1_slope = FALSE
        ),
        variables = list(enso = list(lags = list(-11L, 11L)))
      ),
      occurrence.kernel = list(
        classes = 1L,
        form = "gaussian",
        identify = "seed",
        ksd_fixed = 1.5,
        prior_beta_sd = 1L,
        init = "given",
        re_slope = list(
          enable = TRUE,
          country = TRUE,
          sd_country = 0.35,
          sd_adm1 = 0.35,
          adm1 = TRUE
        ),
        prior_grouping = list(
          source = "enso_grouping_prod_run_prior",
          membership = TRUE
        ),
        spatial = list(enable = FALSE),
        alpha = 1,
        repel = list(
          enable = FALSE
        )
      ),
      occurrence.brms_kl = list(
        basehaz_poisson = "seasonal",
        sampling = 300L,
        warmup = 400L,
        chains = 8L,
        cores = 8L,
        threads_per_chain = 4L,
        init = 0.2,
        re_uncorrelated = TRUE,
        adapt_delta = 0.9,
        max_treedepth = 10L,
        prior_sd_region_country = "normal(0, 1)",
        prior_sd_adm1 = "normal(0, 0.35)",
        prior_sds = "normal(0, 0.35)",
        coef_prior_b = list(
          nbr_lag1_r500_rn_parthmmclusters = "normal(0, 2)",
          ongoing_outbreak_lag1TRUE = "normal(0, 1)"
        ),
        coef_prior_sd = list(
          "country:ADM1" = list(
            nbr_lag1_r500_rn_parthmmclusters = "normal(0, 1)"
          )
        )
      )
    )
  ),
  list(row_transform = .prod_run_kl_cell_row(.PROD_RUN_KL_CELLS))
)

# synth_run1: the prod_run specs on the synthetic world
# (00_ingest_synthetic.R; fake countries, analysis/groupings/synthetic_countries.csv):
# same model terms and priors on the synthetic registries (fit grouping
# enso_grouping_synth, contagion partition hmm_clusters_synth), a four-lag grid
# around the simulated truth (group A lag -3, group B lag -1) and a laptop
# sampler: 2 chains x (150 warm-up + 150 sampling), one thread each.
# 4 lags x 2 groups x 3 tracks = 24 fits. Deploy configs/deploy/synth_run1.yml;
# driver analysis/bash/tools/run_synthetic_pipeline.sh.
.SYNTH_RUN1_LAGS <- c(-5L, -3L, -1L, 1L)
.SYNTH_RUN1_SAMPLER <- list(
  warmup = 150L,
  sampling = 150L,
  chains = 2L,
  cores = 2L,
  threads_per_chain = 1L,
  adapt_delta = 0.90
)
.derive_synth_run1 <- function(code) {
  tr <- c(occ = "occurrence", dur = "duration", mon = "monthly")[[code]]
  run <- "prod_run"
  spec <- BATCH_REGISTRY[[sprintf("%s_%s_grouped_enso_singlelag", run, code)]]
  spec$grouping <- "enso_grouping_synth"
  if (!is.null(spec$overrides$contagion)) {
    # occurrence only
    spec$overrides$contagion$partition_source <- "hmm_clusters_synth"
  }
  key <- paste0(tr, ".dlnm.defaults.single_lag")
  if (is.null(spec$grid[[key]])) {
    stop("synth_run1: expected the prod_run lag grid ", key)
  }
  spec$grid[[key]] <- .SYNTH_RUN1_LAGS
  if (identical(tr, "monthly")) {
    b <- spec$overrides[["monthly"]][["brms_cb"]]
    b[names(.SYNTH_RUN1_SAMPLER)] <- .SYNTH_RUN1_SAMPLER
    spec$overrides[["monthly"]][["brms_cb"]] <- b
  } else {
    k <- paste0(tr, ".brms_cb")
    b <- spec$overrides[[k]]
    b[names(.SYNTH_RUN1_SAMPLER)] <- .SYNTH_RUN1_SAMPLER
    spec$overrides[[k]] <- b
  }
  spec$description <- paste0(
    "[synth run 1] ",
    tr,
    ": the ",
    run,
    " spec on the synthetic world ",
    "(fake countries), lags ",
    paste(.SYNTH_RUN1_LAGS, collapse = ", "),
    ", 2 chains x (150 + 150)"
  )
  spec
}
for (.code in c("occ", "dur", "mon")) {
  BATCH_REGISTRY[[sprintf(
    "synth_run1_%s_grouped_enso_singlelag",
    .code
  )]] <- .derive_synth_run1(.code)
}
rm(.code)


#' Why a batch cannot be generated yet, or NULL: its grouping registry is not
#' written (a sweep rollout's fit batches before gen_grouping) or its own
#' `spec$deferred()` says so. Every batch outside a sweep rollout has a
#' registry and no `deferred`, so it is never deferred.
.batch_deferred <- function(spec) {
  if (!is.null(spec$grouping) && !file.exists(registry_path(spec$grouping))) {
    return(paste0(
      "grouping registry ",
      registry_path(spec$grouping),
      " not written yet (deploy_prod_run.R --stage gen_grouping)"
    ))
  }
  if (is.function(spec$deferred)) {
    return(spec$deferred())
  }
  NULL
}


# ---------------------------------------------------------------------------
# make_config_filename
#
# Build the semantic YAML filename stem (no extension) for one grid row.
# Pattern:
#   {batch}_{backend}[_bh-{basehaz}][_ctry-{country}][_{abbrev}{val}...]
#
# Backend is read from occurrence.backend (new schema) or backend (legacy).
# Basehaz token uses the key matching the row's backend.
# ---------------------------------------------------------------------------

make_config_filename <- function(batch_name, row, varying_dims, defaults) {
  # Backend from row: try new nested key, fall back to legacy flat key
  be <- row[["occurrence.backend"]] %||% row[["backend"]] %||% "brms_cb"

  # Backend token — omit if only one backend is varying
  be_key <- if ("occurrence.backend" %in% varying_dims) {
    "occurrence.backend"
  } else if ("backend" %in% varying_dims) {
    "backend"
  } else {
    NULL
  }
  be_token <- if (!is.null(be_key)) paste0("_", be) else ""

  # Basehaz token — use the key that matches this row's backend
  bh_key_for_be <- switch(
    be,
    brms = "occurrence.brms.basehaz_poisson",
    NA_character_
  )
  if (is.na(bh_key_for_be) || !bh_key_for_be %in% varying_dims) {
    bh_key_for_be <- if ("basehaz_poisson" %in% varying_dims) {
      "basehaz_poisson"
    } else {
      NA_character_
    }
  }
  bh_token <- if (!is.na(bh_key_for_be)) {
    paste0("_bh-", row[[bh_key_for_be]])
  } else {
    ""
  }

  # Country token
  ctry_token <- if ("country" %in% varying_dims) {
    paste0("_ctry-", row[["country"]])
  } else {
    ""
  }

  # Remaining extra dims (anything not already handled above)
  handled <- c(
    "occurrence.backend",
    "backend",
    "occurrence.brms.basehaz_poisson",
    "basehaz_poisson",
    "country"
  )
  extra_dims <- setdiff(varying_dims, handled)
  extra_token <- paste(
    vapply(
      extra_dims,
      function(dim) {
        abbrev <- if (dim %in% names(DIM_ABBREV)) DIM_ABBREV[[dim]] else dim
        # single_lag goes through sl_token(): stripping the sign would make lag -3
        # and lag 3 share one filename (a lead sweep).
        val_str <- if (identical(abbrev, "sl")) {
          sl_token(row[[dim]])
        } else {
          gsub("[^A-Za-z0-9]", "", as.character(row[[dim]]))
        }
        paste0("_", abbrev, "-", val_str)
      },
      character(1L)
    ),
    collapse = ""
  )

  paste0(batch_name, be_token, bh_token, ctry_token, extra_token)
}


# ---------------------------------------------------------------------------
# generate_batch
#
# Expand one BATCH_REGISTRY entry into per-config YAMLs + manifest.yml.
# Frozen batches are skipped (returns NULL invisibly).
# ---------------------------------------------------------------------------

# Every artefact name of a generated fit config fits NAME_MAX
# (paths.R::check_run_file_names()), checked as the config is written so a too-long
# name fails here and not after its fit. A config over the limit is removed and
# the batch stops; one that does not resolve here (e.g. its registry is absent
# locally) is left to its fit worker, which runs the same check before fitting.
.check_config_run_names <- function(cfg_file, model) {
  o <- tryCatch(
    {
      out <- NULL
      invisible(capture.output(
        out <- suppressMessages(suppressWarnings(
          make_options_from_config(
            option_list = list(),
            cli_argv = c("--config", cfg_file)
          )
        ))
      ))
      out
    },
    error = function(e) NULL
  )
  if (is.null(o)) {
    return(invisible(NULL))
  }
  tryCatch(check_run_file_names(o, model), error = function(e) {
    unlink(cfg_file)
    stop(basename(cfg_file), ": ", conditionMessage(e), call. = FALSE)
  })
  invisible(NULL)
}

generate_batch <- function(batch_name, force_frozen = FALSE) {
  spec <- BATCH_REGISTRY[[batch_name]]
  if (is.null(spec)) {
    stop("Unknown batch: '", batch_name, "'")
  }
  if (isTRUE(spec$frozen) && !force_frozen) {
    cat(
      "-- Batch '",
      batch_name,
      "' is frozen; skipping (use force_frozen=TRUE to override)\n",
      sep = ""
    )
    return(invisible(NULL))
  }
  # Not generatable yet (registry not written / spec$deferred): a clear stop
  # instead of load_grouping()'s missing-file error. Never true for a batch
  # whose registry exists and that has no `deferred` (every batch outside a sweep rollout).
  why <- .batch_deferred(spec)
  if (!is.null(why)) {
    stop("Batch '", batch_name, "' is deferred: ", why, call. = FALSE)
  }

  # Country-grouping expansion: when spec$grouping names a registry, expand one
  # config per group. Read analysis/groupings/<source>.csv, derive the distinct
  # group names, inject them as a synthetic grid dimension `grouping.group`, and
  # bake grouping.enable/source into the overrides so every config carries them.
  # Per-config `country` is left to fit-time resolve_group_countries(); the
  # registry stays the single source of truth.
  if (!is.null(spec$grouping)) {
    grp_source <- spec$grouping
    grp_reg <- load_grouping(grp_source)
    grp_names <- unique(grp_reg$group)
    if (!is.null(spec$grouping_exclude)) {
      unknown <- setdiff(spec$grouping_exclude, grp_names)
      if (length(unknown) > 0L) {
        stop(
          "grouping_exclude names no group in '",
          grp_source,
          "': ",
          paste(unknown, collapse = ", "),
          call. = FALSE
        )
      }
      grp_names <- setdiff(grp_names, spec$grouping_exclude)
    }
    spec$grid[["grouping.group"]] <- grp_names
    if (is.null(spec$overrides)) {
      spec$overrides <- list()
    }
    spec$overrides[["grouping.enable"]] <- TRUE
    spec$overrides[["grouping.source"]] <- grp_source
  }

  dir_path <- file.path(get_config_dir(), batch_name)
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

  # Determine which model block this batch belongs to, so that the non-applicable
  # block(s) (inherited from config_defaults.yml) can be stripped from every
  # config. Detection heuristic: check whether spec keys start with
  # "occurrence.", "size.", "duration.", or "monthly.".
  all_spec_keys <- c(names(spec$overrides), names(spec$grid))
  batch_model <- if (any(startsWith(all_spec_keys, "occurrence."))) {
    "occurrence"
  } else if (
    any(startsWith(all_spec_keys, "size.")) || !is.null(spec$overrides$size)
  ) {
    "size"
  } else if (
    any(startsWith(all_spec_keys, "duration.")) ||
      !is.null(spec$overrides$duration)
  ) {
    "duration"
  } else if (
    any(startsWith(all_spec_keys, "monthly.")) ||
      !is.null(spec$overrides$monthly)
  ) {
    "monthly"
  } else {
    NULL # ambiguous — keep all blocks (should not happen in practice)
  }

  grid_df <- as_tibble(do.call(
    expand.grid,
    c(spec$grid, list(stringsAsFactors = FALSE))
  ))
  varying_dims <- names(spec$grid)[vapply(spec$grid, length, integer(1L)) > 1L]
  # Optional row selection (spec$rows_fn, e.g. the deployed-lag rows of the
  # deployed run's no-contagion batch), applied after the varying dimensions are
  # fixed so the filenames match the unfiltered batch's. Absent = every row.
  if (is.function(spec$rows_fn)) {
    grid_df <- spec$rows_fn(grid_df)
  }
  defaults <- get_default_config()

  # Fill occurrence.backend for filename generation if not explicitly in grid.
  # Reads from: overrides occurrence.backend key → overrides$occurrence$backend →
  #             config_defaults occurrence.backend → hard default "brms_cb".
  if (
    !"occurrence.backend" %in% names(grid_df) && !"backend" %in% names(grid_df)
  ) {
    occ_be <- spec$overrides[["occurrence.backend"]] %||%
      spec$overrides$occurrence$backend %||%
      defaults$occurrence$backend %||%
      "brms_cb"
    grid_df[["occurrence.backend"]] <- occ_be
  }

  config_names <- character(nrow(grid_df))

  for (i in seq_len(nrow(grid_df))) {
    row <- as.list(grid_df[i, ])
    config <- defaults

    # Dotted keys ("occurrence.dlnm.defaults.single_lag") set nested leaves via
    # config.R::config_set_leaf().

    # Apply batch-level overrides first
    if (!is.null(spec$overrides)) {
      for (k in names(spec$overrides)) {
        config <- config_set_leaf(config, k, spec$overrides[[k]])
      }
    }
    # Apply this row's grid values (may override batch-level overrides)
    for (k in names(row)) {
      config <- config_set_leaf(config, k, row[[k]])
    }
    # Optional per-row transform (e.g. to inject a concrete ENSO slope key)
    if (!is.null(spec$row_transform)) {
      config <- spec$row_transform(config)
    }

    config <- apply_anchored_defaults(config, batch_model, spec, all_spec_keys)

    # One config = one model run: strip the non-applicable model block(s) so
    # that 04_fit_run.R dispatcher does not misroute (e.g. size configs must
    # not carry an occurrence: block inherited from config_defaults.yml).
    if (!is.null(batch_model)) {
      strip_models <- setdiff(
        c("occurrence", "size", "duration", "monthly"),
        batch_model
      )
      for (sm in strip_models) {
        config[[sm]] <- NULL
      }

      # A generated fit-batch config is never a simulator run (make_batch_configs.R
      # never emits `simulator:`). simulator has fill_when_absent = FALSE in
      # .CONFIG_TOPLEVEL_BLOCKS (config.R), so it stays NULL after merge —
      # exactly what test_config_block_registry.R already exercises.
      config[["simulator"]] <- NULL

      # Nor a results-phase run. A dumped `results:` block states
      # results.inputs.episode_key (default lp_season), and
      # make_options_from_config() refuses a YAML whose top-level episodes.key
      # disagrees with a stated one, so a gid_season fit config would fail to
      # resolve. results has
      # fill_when_absent = TRUE, so the loader still fills its defaults.
      config[["results"]] <- NULL

      config[[batch_model]] <- strip_backend_siblings(config[[batch_model]])
      # Opt-in per spec: drop the superseded prior_sd.
      if (isTRUE(spec$drop_superseded_prior_sd)) {
        config[[batch_model]] <- .drop_superseded_prior_sd(config[[
          batch_model
        ]])
      }
    }

    cfg_stem <- make_config_filename(batch_name, row, varying_dims, defaults)
    config_names[i] <- paste0(cfg_stem, ".yml")

    # Coerce integer fields so write_yaml emits 10, not 10.0
    for (f in c("iter", "warmup", "chains", "cores", "seed")) {
      if (!is.null(config[[f]])) config[[f]] <- as.integer(config[[f]])
    }

    cfg_file <- assert_basename_length(
      file.path(dir_path, paste0(cfg_stem, ".yml")),
      "make_config_filename"
    )
    yaml::write_yaml(config, cfg_file)
    if (!is.null(batch_model)) .check_config_run_names(cfg_file, batch_model)
  }

  # Write (or update) manifest.yml — preserve original created date
  manifest_path <- file.path(dir_path, "manifest.yml")
  existing_created <- if (file.exists(manifest_path)) {
    yaml::read_yaml(manifest_path)$created
  } else {
    as.character(Sys.Date())
  }

  manifest <- list(
    batch = batch_name,
    description = if (!is.null(spec$description)) spec$description else "",
    created = existing_created,
    varies = as.list(varying_dims),
    fixed = c(
      spec$overrides,
      spec$grid[setdiff(names(spec$grid), varying_dims)]
    ),
    n_configs = nrow(grid_df),
    configs = as.list(config_names)
  )
  yaml::write_yaml(manifest, manifest_path)

  cat(sprintf(
    "-- Generated batch '%s': %d config(s) -> %s\n",
    batch_name,
    nrow(grid_df),
    dir_path
  ))
  invisible(config_names)
}


# ---------------------------------------------------------------------------
# generate_all_batches  — convenience wrapper
# ---------------------------------------------------------------------------

generate_all_batches <- function(force_frozen = FALSE) {
  # A deferred batch (.batch_deferred(): registry not written yet, or its own
  # deferred() reason) is skipped with a message instead of stopping --all.
  invisible(lapply(names(BATCH_REGISTRY), function(b) {
    why <- .batch_deferred(BATCH_REGISTRY[[b]])
    if (!is.null(why)) {
      cat("-- Batch '", b, "' deferred: ", why, "\n", sep = "")
      return(invisible(NULL))
    }
    generate_batch(b, force_frozen = force_frozen)
  }))
}


# ---------------------------------------------------------------------------
# discover_batch
#
# Return a tibble of config metadata for a batch directory.
# Reads manifest.yml if present; falls back to globbing *.yml.
# ---------------------------------------------------------------------------

discover_batch <- function(batch_name_or_dir) {
  dir_path <- if (startsWith(batch_name_or_dir, "/")) {
    batch_name_or_dir
  } else {
    file.path(get_config_dir(), batch_name_or_dir)
  }
  if (!dir.exists(dir_path)) {
    stop("Batch directory not found: ", dir_path)
  }

  manifest_path <- file.path(dir_path, "manifest.yml")
  config_files <- if (file.exists(manifest_path)) {
    manifest <- yaml::read_yaml(manifest_path)
    file.path(dir_path, unlist(manifest$configs))
  } else {
    fls <- sort(dir(dir_path, pattern = "\\.yml$", full.names = TRUE))
    fls[basename(fls) != "manifest.yml"]
  }

  purrr::map_dfr(config_files, function(f) {
    cfg <- yaml::read_yaml(f)
    be <- cfg$occurrence$backend %||% cfg$size$backend %||% "brms_cb"
    brms <- cfg$occurrence[[be]] %||% cfg$occurrence$brms
    # brms configs store post-warmup draws as `sampling` (warmup + sampling =
    # total iter); pre-migration configs may still have `iter` directly.
    brms_iter <- if (!is.null(brms$sampling)) {
      (brms$warmup %||% NA_integer_) + brms$sampling
    } else {
      brms$iter %||% NA_integer_
    }
    tibble::tibble(
      config_path = f,
      config_name = basename(f),
      backend = be,
      country = cfg$country %||% "all",
      basehaz_poisson = cfg$occurrence$brms_cb$basehaz_poisson %||%
        cfg$occurrence$brms$basehaz_poisson %||%
        NA_character_,
      iter = brms_iter %||% NA_integer_,
      warmup = cfg$occurrence$brms_cb$warmup %||%
        cfg$occurrence$brms$warmup %||%
        NA_integer_,
      chains = cfg$occurrence$brms_cb$chains %||%
        cfg$occurrence$brms$chains %||%
        NA_integer_,
      seed = cfg$seed %||% NA_integer_,
      config = list(cfg)
    )
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x


# ---------------------------------------------------------------------------
# list_batches — print registry summary
# ---------------------------------------------------------------------------

list_batches <- function() {
  cat(sprintf("%-35s  %s\n", "batch", "description"))
  cat(strrep("-", 80), "\n")
  for (nm in names(BATCH_REGISTRY)) {
    spec <- BATCH_REGISTRY[[nm]]
    frozen <- if (isTRUE(spec$frozen)) " [frozen]" else ""
    n_cfgs <- tryCatch(
      nrow(do.call(expand.grid, c(spec$grid, list(stringsAsFactors = FALSE)))),
      error = function(e) "?"
    )
    cat(sprintf(
      "%-35s  %s  (%s config(s))%s\n",
      nm,
      spec$description,
      n_cfgs,
      frozen
    ))
  }
}


# ---------------------------------------------------------------------------
# CLI entrypoint (only runs when script is executed, not when sourced)
# ---------------------------------------------------------------------------
# sys.nframe() == 0L is TRUE only at Rscript top level, never inside source() /
# sys.source() -- the same guard as make_forecast_configs.R /
# make_hindcast_configs.R, so utils/deploy_prod_run.R (--stage gen_batches) can
# load this file and call generate_batch() directly. Matching
# "make_batch_configs.R" in commandArgs() would also fire under
# `Rscript -e 'source("analysis/utils/make_batch_configs.R")'`.

if (sys.nframe() == 0L) {
  option_list <- list(
    make_option(
      c("--batch"),
      default = NULL,
      type = "character",
      help = "Name of batch to generate (see --list for options)"
    ),
    make_option(
      c("--all"),
      default = FALSE,
      action = "store_true",
      help = "Generate all non-frozen batches"
    ),
    make_option(
      c("--list"),
      default = FALSE,
      action = "store_true",
      help = "Print the batch registry and exit"
    ),
    make_option(
      c("--force_frozen"),
      default = FALSE,
      action = "store_true",
      help = "Regenerate even frozen batches (use with care)"
    )
  )

  cli_opt <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))

  if (cli_opt$list) {
    list_batches()
  } else if (!is.null(cli_opt$batch)) {
    generate_batch(cli_opt$batch, force_frozen = cli_opt$force_frozen)
  } else if (cli_opt$all) {
    generate_all_batches(force_frozen = cli_opt$force_frozen)
  } else {
    cat("Specify --batch <name>, --all, or --list.\n")
    cat("Run with --list to see available batches.\n")
  }
}
