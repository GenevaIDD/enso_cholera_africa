# formula_legacy.R — frozen pre-resolve_slopes() occurrence formula builder.
#
# Provides:
#   build_formula_terms_legacy()   reference builder, never modified
#
# Details:
#   Kept side by side: the occurrence config-merge test compares
#   build_formula_terms() against it. Callers source base.R.

# ---------------------------------------------------------------------------
# Formula builders (used by the three occurrence run scripts)
# ---------------------------------------------------------------------------

#' .occ_group_by_level
#'
#' Per-level group strings for occurrence-model formula building, accounting
#' for multi-country / multi-region data availability.
#' NULL entries indicate levels absent from the current data subset.
#'
#' Mirrors the group_by_level list in build_size_formula() but applies
#' data-availability filtering and uses occurrence-specific column names
#' (gid for ADM2).  Called by build_formula_terms().
#'
#' @param n_country integer
#' @param n_region  integer
#' @return named list(region, country, adm1, adm2), each character or NULL
.occ_group_by_level <- function(n_country, n_region) {
  multi_country <- n_country > 1
  multi_region <- n_region > 1

  if (multi_region && multi_country) {
    list(
      region = "AFRO_region",
      country = "AFRO_region:country",
      adm1 = "AFRO_region:country:ADM1",
      adm2 = "AFRO_region:country:ADM1:gid"
    )
  } else if (multi_country) {
    list(
      region = NULL,
      country = "country",
      adm1 = "country:ADM1",
      adm2 = "country:ADM1:gid"
    )
  } else {
    list(
      region = NULL,
      country = NULL,
      adm1 = "ADM1",
      adm2 = "ADM1:gid"
    )
  }
}


#' build_formula_terms_legacy
#'
#' Legacy occurrence formula builder.  Uses a hardcoded climate_slope_to_cov
#' lookup table and emits /-nested RE terms (e.g. AFRO_region/country/ADM1).
#' Kept side-by-side with build_formula_terms() per project convention.
#' New code should use build_formula_terms() instead.
#'
#' Returns a list:
#'   $fixed_terms  — character vector of covariate names (no intercept)
#'   $re_terms     — character vector of "(x - 1 | g)" strings
#'
#' @param opt         resolved options list (must have opt$occurrence)
#' @param n_country   integer, number of country levels in the data
#' @param n_region    integer, number of AFRO_region levels in the data
build_formula_terms_legacy <- function(opt, n_country, n_region) {
  occ <- opt$occurrence
  covs <- occ$covariates %||% list()
  res <- occ$random_effects %||% list()
  be <- occ$backend %||% "brms_cb"
  re_sep <- if (isTRUE(occ[[be]]$re_uncorrelated)) "||" else "|"

  # --- Fixed-effect terms ---

  # Escape-hatch non-climate covariates
  fixed_terms <- character(0)
  if (isTRUE(covs$time_to_last_outbreak_cat)) {
    fixed_terms <- c(fixed_terms, "time_to_last_outbreak")
  }
  if (isTRUE(covs$ongoing_outbreak_lag1)) {
    fixed_terms <- c(fixed_terms, "ongoing_outbreak_lag1")
  }

  # Convention 2: climate slope flags imply climate fixed effects.
  # Canonical covariate order.
  climate_slope_to_cov <- c(
    re_adm1_mean_precip_slope = "mean_precip",
    re_adm1_mean_precip_lag1_slope = "mean_precip_lag1",
    re_adm1_mean_max_precip_slope = "mean_max_precip",
    re_adm1_mean_max_precip_lag1_slope = "mean_max_precip_lag1",
    re_adm1_mean_mean_cum7_slope = "mean_mean_cum7",
    re_adm1_mean_mean_cum7_lag1_slope = "mean_mean_cum7_lag1",
    re_adm1_mean_max_cum7_slope = "mean_max_cum7",
    re_adm1_mean_max_cum7_lag1_slope = "mean_max_cum7_lag1",
    re_adm1_dev_precip_slope = "dev_precip",
    re_adm1_dev_precip_lag1_slope = "dev_precip_lag1",
    re_adm1_dev_max_precip_slope = "dev_max_precip",
    re_adm1_dev_max_precip_lag1_slope = "dev_max_precip_lag1",
    re_adm1_dev_mean_cum7_slope = "dev_mean_cum7",
    re_adm1_dev_mean_cum7_lag1_slope = "dev_mean_cum7_lag1",
    re_adm1_dev_max_cum7_slope = "dev_max_cum7",
    re_adm1_dev_max_cum7_lag1_slope = "dev_max_cum7_lag1",
    re_adm1_mean_temp_daymean_slope = "mean_temp_daymean",
    re_adm1_mean_temp_daymean_lag1_slope = "mean_temp_daymean_lag1",
    re_adm1_mean_temp_daymax_slope = "mean_temp_daymax",
    re_adm1_mean_temp_daymax_lag1_slope = "mean_temp_daymax_lag1",
    re_adm1_mean_temp_daymin_slope = "mean_temp_daymin",
    re_adm1_mean_temp_daymin_lag1_slope = "mean_temp_daymin_lag1",
    re_adm1_dev_temp_daymean_slope = "dev_temp_daymean",
    re_adm1_dev_temp_daymean_lag1_slope = "dev_temp_daymean_lag1",
    re_adm1_dev_temp_daymax_slope = "dev_temp_daymax",
    re_adm1_dev_temp_daymax_lag1_slope = "dev_temp_daymax_lag1",
    re_adm1_dev_temp_daymin_slope = "dev_temp_daymin",
    re_adm1_dev_temp_daymin_lag1_slope = "dev_temp_daymin_lag1"
  )
  for (slope_key in names(climate_slope_to_cov)) {
    if (isTRUE(res[[slope_key]])) {
      fixed_terms <- c(fixed_terms, unname(climate_slope_to_cov[slope_key]))
    }
  }

  # Calendar escape hatch
  if (isTRUE(covs$months)) {
    fixed_terms <- c(fixed_terms, "months")
  }
  if (isTRUE(covs$quarter)) {
    fixed_terms <- c(fixed_terms, "quarter")
  }

  # ENSO fixed effect — driven by covariates$enso toggle
  enso_col <- NULL
  if (isTRUE(covs$enso)) {
    lag_n <- as.integer(opt$enso$lag_months)
    enso_col <- sprintf("enso_lag%d", lag_n)
    fixed_terms <- c(fixed_terms, enso_col)
  }

  # --- RE group helpers ---
  multi_country <- n_country > 1
  multi_region <- n_region > 1

  slope_group <- function() {
    if (multi_region && multi_country) {
      "AFRO_region/country/ADM1"
    } else if (multi_country) {
      "country/ADM1"
    } else {
      "ADM1"
    }
  }
  adm2_group <- function() {
    if (multi_region && multi_country) {
      "AFRO_region/country/ADM1/gid"
    } else if (multi_country) {
      "country/ADM1/gid"
    } else {
      "ADM1/gid"
    }
  }
  country_group <- function() {
    if (multi_region && multi_country) {
      "AFRO_region/country"
    } else if (multi_country) {
      "country"
    } else {
      NA_character_
    }
  }

  # --- RE terms ---
  re_terms <- character(0)

  if (isTRUE(res$re_region_country_adm1_intercept)) {
    int_grp <- if (isTRUE(res$re_adm2_intercept)) {
      adm2_group()
    } else {
      slope_group()
    }
    re_terms <- c(re_terms, sprintf("(1 %s %s)", re_sep, int_grp))
  }
  if (isTRUE(res$re_region_country_ttlo_slope) && !is.na(country_group())) {
    re_terms <- c(
      re_terms,
      sprintf("(time_to_last_outbreak - 1 %s %s)", re_sep, country_group())
    )
  }

  # Slope REs for each climate covariate (active if slope flag = TRUE AND FE is in fixed_terms)
  for (slope_key in names(climate_slope_to_cov)) {
    cov <- unname(climate_slope_to_cov[slope_key])
    if (isTRUE(res[[slope_key]]) && cov %in% fixed_terms) {
      re_terms <- c(
        re_terms,
        sprintf("(%s - 1 %s %s)", cov, re_sep, slope_group())
      )
    }
  }

  if (
    isTRUE(res$re_region_country_ongoing_slope) &&
      isTRUE(covs$ongoing_outbreak_lag1) &&
      !is.na(country_group())
  ) {
    re_terms <- c(
      re_terms,
      sprintf("(ongoing_outbreak_lag1 - 1 %s %s)", re_sep, country_group())
    )
  }

  # ENSO slope REs — two concrete key families, both require the fixed effect active
  if (!is.null(enso_col)) {
    lag_n <- as.integer(opt$enso$lag_months)
    shallow_key <- sprintf("re_region_country_enso_lag%d_slope", lag_n)
    deep_key <- sprintf("re_adm1_enso_lag%d_slope", lag_n)
    if (isTRUE(res[[shallow_key]])) {
      grp <- if (!is.na(country_group())) country_group() else slope_group()
      re_terms <- c(re_terms, sprintf("(%s - 1 %s %s)", enso_col, re_sep, grp))
    }
    if (isTRUE(res[[deep_key]])) {
      re_terms <- c(
        re_terms,
        sprintf("(%s - 1 %s %s)", enso_col, re_sep, slope_group())
      )
    }
  }

  # Basehaz group-varying seasonality (factor basehaz only)
  occ_bh <- if (be == "brms") {
    occ$brms$basehaz_poisson %||% "seasonal"
  } else {
    NULL
  }
  if (
    isTRUE(res$re_region_country_basehaz_slope) &&
      identical(occ_bh, "factor") &&
      !is.na(country_group())
  ) {
    re_terms <- c(
      re_terms,
      sprintf("(months - 1 %s %s)", re_sep, country_group())
    )
  }

  list(fixed_terms = fixed_terms, re_terms = re_terms)
}
