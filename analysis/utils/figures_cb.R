# figures_cb.R — shared helpers for the DLNM crossbasis figures: ONI grids,
#   exposure-response draws, RE draws, IRR extractors and Africa maps.
#
# Provides:
#   build_onebasis_grid(), eval_erf(), erf_from_draws()   ONI exposure curves
#   get_re_vec(), get_re_draws_mat(), brms_group_re_scaffold()   RE draws
#   extract_adm1_logrr_draws(), extract_adm1_irr(), extract_country_irr(),
#   extract_country_erf(), country_irr_from_adm1()        extractors
#   plot_irr_map(), mask_hatch(), ps_diverging_fill_scale()   map layers
#   load_africa_sf(), load_africa_basemap(), load_taxdat_afr_sf()   basemaps
#
# Details:
#   Used by 06_fit_figures_cb_grouped*.R, 07d counterfactual figures, 01 data
#   overview, 09b / 10 and results_core.R (extractors, COMPOSITE_CONDS).
#   All functions are pure (no side effects, no global assignments); callers
#   load packages and source config.R.

# ---------------------------------------------------------------------------
# build_onebasis_grid
# ---------------------------------------------------------------------------
#' Set up the ONI evaluation grid and centred basis matrix for a single-lag fit.
#'
#' @param vm     cb_meta entry for the "enso" variable
#'               (from readRDS(make_run_file(opt, "cb_meta" | "size_cb_meta")))
#' @param n_at   number of grid points (default 50)
#' @return named list:
#'   at_grid       numeric [n_at]   — ONI evaluation grid
#'   Z_grid        matrix [n_at × K] — onebasis (uncentred)
#'   Z_cen         matrix [n_at × K] — onebasis centred at ONI = 0
#'   K             integer           — number of CB columns
#'   cb_names      character [K]     — CB column names
#'   enso_map_vals named numeric     — c("La Niña (–1.0)"=-1, "Neutral (0)"=0, "El Niño (+1.0)"=1)
#'   enso_map_idxs named integer     — nearest at_grid index per condition
#'   enso_map_ats  named numeric     — actual at_grid value used per condition
build_onebasis_grid <- function(
  vm,
  n_at = 50L,
  map_vals = c(
    "La Ni\u00f1a (\u22121.0)" = -1.0,
    "Neutral (0)" = 0.0,
    "El Ni\u00f1o (+1.0)" = 1.0
  )
) {
  if (is.null(vm$single_lag)) {
    stop(
      "build_onebasis_grid: vm$single_lag is NULL — ",
      "use this only for single-lag fits (s08c / s14)."
    )
  }

  cb_names <- vm$cb_col_names
  K <- length(cb_names)
  at_grid <- seq(vm$range[1], vm$range[2], length.out = as.integer(n_at))
  ob_args <- vm$onebasis_attrs[names(vm$onebasis_attrs) != "range"]
  # dlnm::onebasis records a spurious `degree` attribute for an `ns` basis that
  # ns() rejects on replay ("unused argument (degree)"); strip it (no-op for ns,
  # and `bs` keeps its degree). Mirrors .onebasis_replay_args() in dlnm.R.
  if (identical(ob_args$fun, "ns")) {
    ob_args$degree <- NULL
  }
  Z_grid <- do.call(dlnm::onebasis, c(list(at_grid), ob_args))
  cen_idx <- which.min(abs(at_grid - 0))
  Z_cen <- Z_grid - matrix(Z_grid[cen_idx, ], nrow(Z_grid), K, byrow = TRUE)

  enso_map_vals <- map_vals
  enso_map_idxs <- vapply(
    enso_map_vals,
    function(v) which.min(abs(at_grid - v)),
    integer(1L)
  )
  enso_map_ats <- at_grid[enso_map_idxs]

  list(
    at_grid = at_grid,
    Z_grid = Z_grid,
    Z_cen = Z_cen,
    K = K,
    cb_names = cb_names,
    enso_map_vals = enso_map_vals,
    enso_map_idxs = enso_map_idxs,
    enso_map_ats = enso_map_ats
  )
}

# ---------------------------------------------------------------------------
# eval_erf
# ---------------------------------------------------------------------------
#' Evaluate the log-IRR exposure-response curve from a coefficient vector.
#'
#' @param b     numeric [K]     — CB coefficient vector
#' @param Z_cen matrix [n × K] — centred onebasis evaluation matrix
#' @return numeric [n] — log-IRR at each grid point
eval_erf <- function(b, Z_cen) as.numeric(Z_cen %*% b)

# ---------------------------------------------------------------------------
# erf_from_draws
# ---------------------------------------------------------------------------
#' Compute an ERF tibble (median + 95% and 50% CI) from a matrix of posterior draws.
#'
#' @param b_draws matrix [n_draws × K] — posterior CB coefficient draws
#' @param Z_cen   matrix [n_at × K]   — centred onebasis evaluation matrix
#' @param at_grid numeric [n_at]       — ONI evaluation grid
#' @return tibble: at, logRR, lo, hi, IRR, IRR_lo, IRR_hi, lo50, hi50,
#'   IRR_lo50, IRR_hi50 (the 50% columns appended last, so the 95% columns keep
#'   their positions)
erf_from_draws <- function(b_draws, Z_cen, at_grid) {
  logRR_mat <- Z_cen %*% t(b_draws) # n_at × n_draws
  tibble::tibble(
    at = at_grid,
    logRR = apply(logRR_mat, 1L, median),
    lo = apply(logRR_mat, 1L, quantile, 0.025),
    hi = apply(logRR_mat, 1L, quantile, 0.975),
    IRR = exp(logRR),
    IRR_lo = exp(lo),
    IRR_hi = exp(hi),
    lo50 = apply(logRR_mat, 1L, quantile, 0.25),
    hi50 = apply(logRR_mat, 1L, quantile, 0.75),
    IRR_lo50 = exp(lo50),
    IRR_hi50 = exp(hi50)
  )
}

# ---------------------------------------------------------------------------
# get_re_vec
# ---------------------------------------------------------------------------
#' Extract the CB RE mean vector for a label from a ranef array.
#' Returns a zero vector if the array is NULL or the label is absent.
#'
#' @param arr      ranef array [n_levels × 2 × n_cols] (summary = TRUE) or NULL
#' @param lbl      character — level label to look up
#' @param K        integer   — number of CB columns
#' @param cb_names character [K] — CB column names
#' @return named numeric [K]
get_re_vec <- function(arr, lbl, K, cb_names) {
  out <- setNames(rep(0.0, K), cb_names)
  if (is.null(arr)) {
    return(out)
  }
  avail <- intersect(cb_names, dimnames(arr)[[3]])
  if (length(avail) == 0L || !lbl %in% dimnames(arr)[[1]]) {
    return(out)
  }
  out[avail] <- arr[lbl, "Estimate", avail]
  out
}

# ---------------------------------------------------------------------------
# get_re_draws_mat
# ---------------------------------------------------------------------------
#' Extract the CB RE draw matrix for a label from a ranef draws array.
#' Returns a zero matrix if the array is NULL or the label is absent.
#'
#' @param arr      ranef draws array [n_draws × n_levels × n_cols] or NULL
#' @param lbl      character — level label to look up
#' @param K        integer   — number of CB columns
#' @param cb_names character [K] — CB column names
#' @return matrix [n_draws × K]
get_re_draws_mat <- function(arr, lbl, K, cb_names) {
  n_d <- if (is.null(arr)) 1L else dim(arr)[1L]
  zero <- matrix(0.0, n_d, K, dimnames = list(NULL, cb_names))
  if (is.null(arr)) {
    return(zero)
  }
  lvls <- dimnames(arr)[[2L]]
  if (!lbl %in% lvls) {
    return(zero)
  }
  matrix(arr[, lbl, ], nrow = n_d, dimnames = list(NULL, cb_names))
}

# ---------------------------------------------------------------------------
# brms_group_re_scaffold
# ---------------------------------------------------------------------------
#' Shared brms group-RE extraction scaffold for the manuscript-results phase.
#'
#' Loads population fixed-effect means and a seeded posterior-draw subsample,
#' selects the RE group keys (AFRO_region dropped under grouping.enable) and
#' returns the per-tier ranef arrays plus a grouping-level string parser. Used
#' by derive_enso_relation()'s extractors and 09's derive_track_params().
#'
#' @param fit       a brms fit object.
#' @param variables character — the model parameters to pull. For the ENSO
#'   smooths these are the crossbasis column names (cb_names); for the parameter
#'   summaries they are the scalar term names (Intercept, climate slopes, ...).
#' @param n_draws   NULL = use all posterior draws; integer = seeded subsample
#'   (set.seed(seed); sample.int(n_total, n_draws)), matching the ENSO
#'   extractors' exact RNG sequence so routed artefacts stay byte-identical.
#' @param seed      RNG seed for the subsample (default 42L, as the ENSO path).
#' @return named list:
#'   has_region   logical — TRUE when the fit keeps the AFRO_region RE tier.
#'   grp          list(region, country, adm1) — the ranef() group names (NA on
#'                a dropped region tier).
#'   K, variables as supplied.
#'   b_pop        named numeric [K] — fixef() posterior means for `variables`.
#'   b_pop_draws  matrix [n_draws × K] — subsampled population draws.
#'   draw_idx     integer — the retained draw indices.
#'   re_mean      the full brms::ranef(summary = TRUE) list.
#'   re_mean_vec(tier, lbl)  → named numeric [K], zero-filled for absent terms.
#'   re_draws_mat(tier, lbl) → matrix [n_draws × K], zero-filled for absent
#'                             terms (tier in "region"/"country"/"adm1").
#'   parse_level(lbl, tier)  → list(region, country, country_lbl[, adm1_gid])
#'                             splitting a flattened grouping-level string.
brms_group_re_scaffold <- function(fit, variables, n_draws = NULL, seed = 42L) {
  variables <- as.character(variables)
  K <- length(variables)

  re_mean <- brms::ranef(fit, summary = TRUE)
  has_region <- "AFRO_region" %in% names(re_mean)
  grp <- list(
    region = if (has_region) "AFRO_region" else NA_character_,
    country = if (has_region) "AFRO_region:country" else "country",
    adm1 = if (has_region) "AFRO_region:country:ADM1" else "country:ADM1"
  )

  b_pop <- brms::fixef(fit)[variables, "Estimate"]
  names(b_pop) <- variables
  b_pop_draws_raw <- unclass(brms::as_draws_matrix(
    fit,
    variable = paste0("b_", variables)
  ))
  colnames(b_pop_draws_raw) <- variables
  n_total <- nrow(b_pop_draws_raw)
  if (is.null(n_draws)) {
    draw_idx <- seq_len(n_total)
  } else {
    set.seed(seed)
    draw_idx <- sample.int(n_total, n_draws)
  }
  b_pop_draws <- b_pop_draws_raw[draw_idx, , drop = FALSE]

  re_draws_all <- brms::ranef(fit, summary = FALSE)
  pull_draws <- function(tier) {
    gname <- grp[[tier]]
    if (is.na(gname) || is.null(re_draws_all[[gname]])) {
      return(NULL)
    }
    re_draws_all[[gname]][draw_idx, , , drop = FALSE]
  }
  d_region <- pull_draws("region")
  d_country <- pull_draws("country")
  d_adm1 <- pull_draws("adm1")
  rm(re_draws_all)
  gc(verbose = FALSE)

  mean_arr <- function(tier) {
    gname <- grp[[tier]]
    if (is.na(gname)) NULL else re_mean[[gname]]
  }

  list(
    has_region = has_region,
    grp = grp,
    K = K,
    variables = variables,
    b_pop = b_pop,
    b_pop_draws = b_pop_draws,
    draw_idx = draw_idx,
    re_mean = re_mean,
    re_mean_vec = function(tier, lbl) {
      get_re_vec(mean_arr(tier), lbl, K, variables)
    },
    re_draws_mat = function(tier, lbl) {
      arr <- switch(
        tier,
        region = d_region,
        country = d_country,
        adm1 = d_adm1,
        NULL
      )
      out <- matrix(0.0, length(draw_idx), K, dimnames = list(NULL, variables))
      if (is.null(arr) || !lbl %in% dimnames(arr)[[2L]]) {
        return(out)
      }
      avail <- intersect(variables, dimnames(arr)[[3L]])
      if (length(avail)) {
        out[, avail] <- matrix(
          arr[, lbl, avail],
          nrow = length(draw_idx),
          dimnames = list(NULL, avail)
        )
      }
      out
    },
    parse_level = function(lbl, tier = c("country", "adm1")) {
      tier <- match.arg(tier)
      parts <- strsplit(lbl, "_")[[1L]]
      n <- length(parts)
      if (tier == "country") {
        if (has_region) {
          list(
            region = parts[1L],
            country = parts[2L],
            country_lbl = paste(parts[1:2], collapse = "_")
          )
        } else {
          list(
            region = NA_character_,
            country = parts[1L],
            country_lbl = parts[1L]
          )
        }
      } else {
        if (has_region) {
          list(
            region = parts[1L],
            country = parts[2L],
            country_lbl = paste(parts[1:2], collapse = "_"),
            adm1_gid = paste(parts[(n - 1L):n], collapse = "_")
          )
        } else {
          list(
            region = NA_character_,
            country = parts[1L],
            country_lbl = parts[1L],
            adm1_gid = paste(parts[2L:n], collapse = "_")
          )
        }
      }
    }
  )
}

# ---------------------------------------------------------------------------
# load_africa_sf
# ---------------------------------------------------------------------------
#' Load Africa background polygon (giscoR) and ADM1 shapefile (GADM gpkg).
#'
#' GID_1 in adm1_sf is normalised to the ADM1 key: "ETH.1_1" → "ETH_1",
#' Ghana's "GHA1_2" → "GHA_1" (same rule as gid.R::gid_to_adm1(); NA if unparseable).
#'
#' @param gadm_path path to the GADM geopackage (default data/geodata/gadm_410-levels.gpkg)
#' @return list(afro_sf, adm1_sf) — adm1_sf is NULL on read failure
load_africa_sf <- function(gadm_path = "data/geodata/gadm_410-levels.gpkg") {
  # Africa country outlines
  # Try loading gisco_countries_2024 first (newer versions), fall back to gisco_countries (older versions)
  tryCatch(
    data(gisco_countries_2024, package = "giscoR", envir = environment()),
    error = function(e) {
      data(gisco_countries, package = "giscoR", envir = environment())
    }
  )
  data(gisco_countrycode, package = "giscoR", envir = environment())

  # Use whichever dataset is available
  gisco_data <- if (exists("gisco_countries_2024", inherits = FALSE)) {
    gisco_countries_2024
  } else {
    gisco_countries
  }

  afro_sf <- gisco_data %>%
    dplyr::inner_join(
      gisco_countrycode %>%
        dplyr::select(ISO3_CODE, continent) %>%
        dplyr::filter(continent == "Africa"),
      by = "ISO3_CODE"
    ) %>%
    nngeo::st_remove_holes() %>%
    rmapshaper::ms_simplify(keep = 0.5, keep_shapes = FALSE)
  # The synthetic world's fake countries (gid.R::synthetic_countries(), call-time):
  # their outlines are the union of their GADM units. Real GADM has no XA? unit,
  # so the query returns no row there and nothing is added.
  syn_iso <- if (exists("synthetic_countries", mode = "function")) {
    synthetic_countries()$iso3
  } else {
    character(0)
  }
  if (length(syn_iso)) {
    syn_sql <- sprintf(
      "SELECT * FROM ADM_1 WHERE GID_0 IN (%s)",
      paste0("'", syn_iso, "'", collapse = ", ")
    )
    syn <- tryCatch(
      sf::st_read(gadm_path, query = syn_sql, quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(syn) && nrow(syn)) {
      geo <- lapply(split(sf::st_geometry(syn), syn$GID_0), function(g) {
        sf::st_union(g)[[1]]
      })
      add <- sf::st_sf(
        ISO3_CODE = names(geo),
        NAME_ENGL = synthetic_countries()$name[match(names(geo), syn_iso)],
        continent = "Africa",
        geometry = sf::st_sfc(geo, crs = sf::st_crs(afro_sf))
      )
      afro_sf <- dplyr::bind_rows(afro_sf, add)
    }
  }

  # ADM1 polygons
  adm1_sf <- tryCatch(
    sf::st_read(gadm_path, layer = "ADM_1", quiet = TRUE) %>%
      dplyr::mutate(
        # The ADM1 key gid.R::gid_to_adm1() gives ("AGO.1_1" -> "AGO_1"), inlined
        # (see the call-time dependency note). Stripping a literal "_1$" missed
        # Ghana, whose GADM 4.1 GID_1 is "GHA1_2" (no ".", version 2).
        GID_1 = stringr::str_extract(
          sub(
            "^([A-Z]{3})([0-9])",
            "\\1_\\2",
            stringr::str_replace_all(GID_1, "\\.", "_")
          ),
          "^[A-Z]{3}_[0-9]+"
        )
      ),
    error = function(e) {
      warning(
        "load_africa_sf: could not load ADM1 layer: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  list(afro_sf = afro_sf, adm1_sf = adm1_sf)
}

#' Drop the small Indian Ocean islands that show up only as specks east of
#' Madagascar on the Africa-wide manuscript maps (Mauritius, Réunion).
#'
#' @param afro_sf country outlines with an ISO3_CODE column (load_africa_sf()$afro_sf)
#' @param iso3    ISO3 codes to drop
#' @return afro_sf without those countries
drop_small_islands <- function(afro_sf, iso3 = c("MUS", "REU")) {
  dplyr::filter(afro_sf, !ISO3_CODE %in% iso3)
}

# ===========================================================================
# ENSO crossbasis posterior extractors + map plotting (shared by
# 06_fit_figures_cb_grouped*.R, 09b_postprocess_results.R, 10_manuscript_figures.R)
#
# Lifted out of 09_postprocess_results.R (extract_adm1_irr, extract_country_erf)
# and the 06_fit_figures_cb_grouped* scripts' inline copies, so the figure
# scripts and the manuscript derivation run one implementation. All pure;
# every call is `::`-qualified, callers load nothing beyond sourcing this file.
# ===========================================================================

# ADM1 / country IRR evaluation points (Figure 2; order = left to right).
# The manuscript results phase reads its set from the config
# (results.derivation.enso_conds, bare numbers) and builds the labels with
# enso_conds_named(), so only numbers round-trip through YAML, never the unicode
# labels. COMPOSITE_CONDS is that rule applied to the default set and stays the
# default for the 06_fit_figures_cb_grouped*.R scripts.
#' Named ONI-condition vector in the composite label convention: negative values
#' "La Niña (−v)" (U+2212 minus), positive "El Niño (+v)", zero "Neutral (0.0)",
#' one decimal, in the order given. The manuscript results config states its
#' conditions as bare numbers (results.derivation.enso_conds); 09 turns them
#' into this named vector so every label is built by one rule.
#' @param vals numeric ONI values (non-empty, finite, distinct)
#' @return named numeric vector
enso_conds_named <- function(vals) {
  vals <- as.numeric(unlist(vals))
  if (
    !length(vals) || anyNA(vals) || any(!is.finite(vals)) || anyDuplicated(vals)
  ) {
    stop(
      "enso_conds_named: ONI conditions must be a non-empty set of distinct finite numbers; got ",
      paste(vals, collapse = ", "),
      call. = FALSE
    )
  }
  lab <- ifelse(
    vals < 0,
    sprintf("La Ni\u00f1a (\u2212%.1f)", abs(vals)),
    ifelse(vals > 0, sprintf("El Ni\u00f1o (+%.1f)", vals), "Neutral (0.0)")
  )
  stats::setNames(vals, enc2utf8(lab))
}

COMPOSITE_CONDS <- enso_conds_named(c(-2, -1, 1, 2))

#' Rows of the centred onebasis matrix at the ENSO conditions a grid was built
#' for -- `build_onebasis_grid(vm, n_at, map_vals = enso_conds)` already records
#' the nearest grid index per condition, so callers need not recompute it.
#'
#' @param grid_obj Result of build_onebasis_grid() called with `map_vals`.
#' @return matrix [n_conditions x K]
cb_target_matrix <- function(grid_obj) {
  grid_obj$Z_cen[grid_obj$enso_map_idxs, , drop = FALSE]
}

#' ADM1-level ENSO IRR with a 95% CI at the given ENSO conditions: b_pop +
#' region + country + ADM1 slope REs, CI from posterior draws.
#' @param n_at ONI evaluation-grid resolution (09: 50; the grouped figure
#'   script: 100).
#' @return tibble(adm1_label, adm1_gid, country, enso_label, enso_val, IRR, lo,
#'   hi, significant, enso_group), or NULL when the fit has no ADM1 RE.
#' adm1_logrr_draws
#'
#' Per-draw log rate ratios of one ADM1 level at the target ENSO conditions:
#' population draws + region + country + ADM1 random-slope draws, mapped
#' through the crossbasis target rows. One copy of the draw combination that
#' extract_adm1_irr() and 09's dir_pop_ci() (utils/results_derive.R) both use.
#'
#' @param sc brms_group_re_scaffold() result
#' @param Z_targets cb_target_matrix() of the onebasis grid (C x K)
#' @param lbl ADM1 RE level label
#' @param parsed sc$parse_level(lbl, "adm1") (passed when already parsed)
#' @return C x D matrix of log rate ratios
adm1_logrr_draws <- function(
  sc,
  Z_targets,
  lbl,
  parsed = sc$parse_level(lbl, "adm1")
) {
  r_rgn_d <- sc$re_draws_mat("region", parsed$region)
  r_ctry_d <- sc$re_draws_mat("country", parsed$country_lbl)
  r_adm1_d <- sc$re_draws_mat("adm1", lbl)
  b_comb_d <- sc$b_pop_draws + r_rgn_d + r_ctry_d + r_adm1_d
  Z_targets %*% t(b_comb_d)
}

#' Per-draw log rate ratios of every ADM1 level of one fit, as one array.
#'
#' Same draw contract as extract_adm1_irr() (grid resolution, seed, subsample,
#' b_pop + region + country + ADM1 slope REs), so the per-ADM1 C x D slices are
#' exactly the matrices extract_adm1_irr() quantiles. 09 caches the result in the
#' enso_relation artefact so the per-draw population sums (dir_pop_ci(), the
#' ENSO-sensitive population) never reload a fit.
#'
#' @return list(logrr = array [A x C x D] (dimnames: adm1_gid, names(enso_conds),
#'   NULL), adm1 = tibble(adm1_label, adm1_gid, country, enso_group)), or NULL
#'   when the fit has no ADM1 RE.
extract_adm1_logrr_draws <- function(
  brms_fit,
  meta,
  enso_conds = COMPOSITE_CONDS,
  group_label = NA_character_,
  n_draws = 400L,
  seed = 42L,
  n_at = 50L
) {
  vm <- meta[["enso"]]
  cb_names <- vm$cb_col_names
  grid_obj <- build_onebasis_grid(vm, n_at = n_at, map_vals = enso_conds)
  Z_targets <- cb_target_matrix(grid_obj)

  sc <- brms_group_re_scaffold(
    brms_fit,
    cb_names,
    n_draws = n_draws,
    seed = seed
  )
  re_adm1_arr <- sc$re_mean[[sc$grp$adm1]]
  if (is.null(re_adm1_arr)) {
    warning("ADM1 RE not found for ", group_label)
    return(NULL)
  }
  adm1_levels <- dimnames(re_adm1_arr)[[1L]]
  parsed_all <- lapply(adm1_levels, function(l) sc$parse_level(l, "adm1"))
  adm1_gids <- vapply(parsed_all, `[[`, character(1L), "adm1_gid")

  logrr <- array(
    NA_real_,
    dim = c(length(adm1_levels), length(enso_conds), n_draws),
    dimnames = list(adm1_gids, names(enso_conds), NULL)
  )
  for (i in seq_along(adm1_levels)) {
    logrr[i, , ] <- adm1_logrr_draws(
      sc,
      Z_targets,
      adm1_levels[[i]],
      parsed_all[[i]]
    )
  }
  list(
    logrr = logrr,
    adm1 = tibble::tibble(
      adm1_label = adm1_levels,
      adm1_gid = adm1_gids,
      country = vapply(parsed_all, `[[`, character(1L), "country"),
      enso_group = group_label
    )
  )
}

extract_adm1_irr <- function(
  brms_fit,
  meta,
  enso_conds = COMPOSITE_CONDS,
  group_label = NA_character_,
  n_draws = 400L,
  seed = 42L,
  n_at = 50L
) {
  vm <- meta[["enso"]]
  cb_names <- vm$cb_col_names
  K <- length(cb_names)
  grid_obj <- build_onebasis_grid(vm, n_at = n_at, map_vals = enso_conds)
  Z_targets <- cb_target_matrix(grid_obj)

  sc <- brms_group_re_scaffold(
    brms_fit,
    cb_names,
    n_draws = n_draws,
    seed = seed
  )
  re_adm1_arr <- sc$re_mean[[sc$grp$adm1]]
  if (is.null(re_adm1_arr)) {
    warning("ADM1 RE not found for ", group_label)
    return(NULL)
  }
  adm1_levels <- dimnames(re_adm1_arr)[[1L]]

  purrr::map_dfr(adm1_levels, function(lbl) {
    parsed <- sc$parse_level(lbl, "adm1")
    r_rgn <- sc$re_mean_vec("region", parsed$region)
    r_ctry <- sc$re_mean_vec("country", parsed$country_lbl)
    r_adm1 <- sc$re_mean_vec("adm1", lbl)
    b_comb <- sc$b_pop + r_rgn + r_ctry + r_adm1
    irr_v <- exp(as.numeric(Z_targets %*% b_comb))

    logRR_d <- adm1_logrr_draws(sc, Z_targets, lbl, parsed)
    lo_v <- apply(logRR_d, 1L, quantile, 0.025)
    hi_v <- apply(logRR_d, 1L, quantile, 0.975)

    tibble::tibble(
      adm1_label = lbl,
      adm1_gid = parsed$adm1_gid,
      country = parsed$country,
      enso_label = names(enso_conds),
      enso_val = unname(enso_conds),
      IRR = irr_v,
      lo = lo_v,
      hi = hi_v,
      significant = (lo_v > 0) | (hi_v < 0),
      enso_group = group_label
    )
  })
}

#' Country-level ENSO IRR with a 95% CI (b_pop + region + country slope REs; no
#' ADM1 deviation). Same contract as extract_adm1_irr() one tier up.
#' @return tibble(country, enso_label, enso_val, IRR, lo, hi, significant,
#'   enso_group), or NULL when the fit has no country RE.
extract_country_irr <- function(
  brms_fit,
  meta,
  enso_conds = COMPOSITE_CONDS,
  group_label = NA_character_,
  n_draws = 400L,
  seed = 42L,
  n_at = 50L
) {
  vm <- meta[["enso"]]
  cb_names <- vm$cb_col_names
  grid_obj <- build_onebasis_grid(vm, n_at = n_at, map_vals = enso_conds)
  Z_targets <- cb_target_matrix(grid_obj)

  sc <- brms_group_re_scaffold(
    brms_fit,
    cb_names,
    n_draws = n_draws,
    seed = seed
  )
  re_ctry_arr <- sc$re_mean[[sc$grp$country]]
  if (is.null(re_ctry_arr)) {
    warning("country RE not found for ", group_label)
    return(NULL)
  }

  purrr::map_dfr(dimnames(re_ctry_arr)[[1L]], function(lbl) {
    parsed <- sc$parse_level(lbl, "country")
    b_comb <- sc$b_pop +
      sc$re_mean_vec("region", parsed$region) +
      sc$re_mean_vec("country", lbl)
    irr_v <- exp(as.numeric(Z_targets %*% b_comb))

    b_comb_d <- sc$b_pop_draws +
      sc$re_draws_mat("region", parsed$region) +
      sc$re_draws_mat("country", lbl)
    logRR_d <- Z_targets %*% t(b_comb_d)
    lo_v <- apply(logRR_d, 1L, quantile, 0.025)
    hi_v <- apply(logRR_d, 1L, quantile, 0.975)

    tibble::tibble(
      country = parsed$country,
      enso_label = names(enso_conds),
      enso_val = unname(enso_conds),
      IRR = irr_v,
      lo = lo_v,
      hi = hi_v,
      significant = (lo_v > 0) | (hi_v < 0),
      enso_group = group_label
    )
  })
}

#' Country-level ENSO exposure-response curves over the whole ONI grid
#' (b_pop + region + country slope REs, centred at ENSO = 0). `grid` is the
#' build_onebasis_grid() result for the group's deployed basis.
#' @return tibble from erf_from_draws() plus country, enso_group, range_lo,
#'   range_hi; NULL when the fit has no country RE.
extract_country_erf <- function(
  brms_fit,
  meta,
  grid,
  group_label = NA_character_,
  n_draws = 400L,
  seed = 42L
) {
  vm <- meta[["enso"]]
  cb_names <- vm$cb_col_names
  K <- length(cb_names)
  at_grid <- grid$at_grid
  Z_cen <- grid$Z_cen

  sc <- brms_group_re_scaffold(
    brms_fit,
    cb_names,
    n_draws = n_draws,
    seed = seed
  )
  re_ctry_arr <- sc$re_mean[[sc$grp$country]]
  if (is.null(re_ctry_arr)) {
    warning("country RE not found for ", group_label)
    return(NULL)
  }
  ctry_levels <- dimnames(re_ctry_arr)[[1L]]

  purrr::map_dfr(ctry_levels, function(lbl) {
    parsed <- sc$parse_level(lbl, "country")
    r_rgn_d <- sc$re_draws_mat("region", parsed$region)
    r_ctry_d <- sc$re_draws_mat("country", lbl)
    erf_from_draws(sc$b_pop_draws + r_rgn_d + r_ctry_d, Z_cen, at_grid) %>%
      dplyr::mutate(
        country = parsed$country,
        enso_group = group_label,
        range_lo = vm$range[1],
        range_hi = vm$range[2]
      )
  })
}

#' Roll ADM1 IRR rows up to the country: median and IQR of the ADM1 point IRRs
#' per (country, ENSO condition). `n_adm1` is the number of ADM1 units pooled.
country_irr_from_adm1 <- function(adm1_irr_df) {
  adm1_irr_df %>%
    dplyr::group_by(country, enso_label) %>%
    dplyr::summarise(
      median_IRR = median(IRR, na.rm = TRUE),
      q25_IRR = quantile(IRR, 0.25, na.rm = TRUE),
      q75_IRR = quantile(IRR, 0.75, na.rm = TRUE),
      n_adm1 = dplyr::n(),
      .groups = "drop"
    )
}

#' Faceted IRR choropleth (one panel per ENSO condition), diverging blue-white-red
#' around IRR = 1 and capped at `limits`. One implementation of the map the
#' grouped and composite figure scripts each drew twice (ADM1 and country).
#'
#' @param irr_df tibble with `enso_label` plus `adm1_gid` (level "adm1") or
#'   `country` (level "country") and `value_col`; an optional logical
#'   `significant` is outlined when `outline_sig`.
#' @param geo list(afro_sf, adm1_sf) from load_africa_sf().
#' @param value_col Column mapped to the fill ("IRR" or e.g. "median_IRR").
#' @param enso_levels Facet order [default: the factor levels / order of appearance].
#' @return a ggplot, or NULL (with a warning) when nothing joins the geometry.
plot_irr_map <- function(
  irr_df,
  geo,
  level = c("adm1", "country"),
  value_col = "IRR",
  limits = c(0.4, 2.5),
  fill_name = "IRR\n(capped)",
  title = NULL,
  subtitle = NULL,
  ncol = 4L,
  outline_sig = TRUE,
  enso_levels = NULL
) {
  level <- match.arg(level)
  if (is.null(irr_df) || !nrow(irr_df)) {
    return(NULL)
  }
  afro_sf <- geo$afro_sf
  key <- if (level == "adm1") "adm1_gid" else "country"
  base <- if (level == "adm1") geo$adm1_sf else afro_sf
  if (is.null(base)) {
    stop("plot_irr_map(): geometry for level '", level, "' is not loaded.")
  }
  lv <- enso_levels %||%
    (if (is.factor(irr_df$enso_label)) {
      levels(irr_df$enso_label)
    } else {
      unique(irr_df$enso_label)
    })
  has_sig <- isTRUE(outline_sig) && "significant" %in% names(irr_df)

  d <- tibble::tibble(
    enso_label = factor(irr_df$enso_label, levels = lv),
    IRR_capped = pmin(pmax(irr_df[[value_col]], limits[1L]), limits[2L])
  )
  d[[key]] <- irr_df[[key]]
  if (has_sig) {
    d$significant <- irr_df$significant
  }

  joined <- dplyr::inner_join(
    base,
    d,
    by = stats::setNames(key, if (level == "adm1") "GID_1" else "ISO3_CODE")
  )
  if (!nrow(joined)) {
    warning("plot_irr_map(): no ", level, " units joined to the geometry.")
    return(NULL)
  }
  sig <- if (has_sig) joined[which(joined$significant %in% TRUE), ] else NULL

  # ADM1 units are filled borderless (the significance outline carries the edge);
  # country polygons keep a thin grey border.
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = afro_sf,
      fill = "gray85",
      color = "gray60",
      lwd = 0.15
    ) +
    ggplot2::geom_sf(
      data = joined,
      ggplot2::aes(fill = IRR_capped),
      color = if (level == "adm1") NA else "gray40",
      lwd = if (level == "adm1") 0.5 else 0.2
    )
  if (!is.null(sig) && nrow(sig)) {
    p <- p +
      ggplot2::geom_sf(
        data = sig,
        fill = NA,
        color = "black",
        lwd = if (level == "adm1") 0.5 else 1.0
      )
  }
  if (level == "adm1") {
    p <- p +
      ggplot2::geom_sf(data = afro_sf, fill = NA, color = "black", lwd = 0.12)
  }
  p +
    ggplot2::facet_wrap(~enso_label, ncol = ncol) +
    ggplot2::scale_fill_gradient2(
      midpoint = 1,
      low = "#2166ac",
      mid = "white",
      high = "#d73027",
      limits = limits,
      oob = scales::squish,
      name = fill_name
    ) +
    taxdat::map_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
}

#' Repeating line-hatch pattern (grid::pattern) usable as a constant ggplot2
#' fill or as a scale_fill_manual() value: diagonal lines, plus the other
#' diagonal when `cross` (the IPCC AR6 Atlas convention: diagonal = no robust
#' signal, crossed = conflicting signal). `bg`: a solid background under the
#' lines (legend keys). Needs pattern_fill_supported() and a pattern-capable
#' device (ragg / cairo / pdf).
mask_hatch <- function(
  spacing_mm = 0.6,
  lwd = 0.5,
  col = "white",
  bg = NA,
  cross = FALSE
) {
  gp <- grid::gpar(col = col, lwd = lwd, lineend = "square")
  lines <- grid::segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 1, gp = gp)
  if (cross) {
    lines <- grid::grobTree(
      lines,
      grid::segmentsGrob(x0 = 0, y0 = 1, x1 = 1, y1 = 0, gp = gp)
    )
  }
  grid::pattern(
    if (is.na(bg)) {
      lines
    } else {
      grid::grobTree(
        grid::rectGrob(gp = grid::gpar(fill = bg, col = NA)),
        lines
      )
    },
    width = grid::unit(spacing_mm, "mm"),
    height = grid::unit(spacing_mm, "mm"),
    extend = "repeat"
  )
}

#' TRUE when ggplot2 accepts grid patterns as fills (>= 3.5.0, R >= 4.1).
pattern_fill_supported <- function() {
  utils::packageVersion("ggplot2") >= "3.5.0" && getRversion() >= "4.1.0"
}

#' taxdat's `afr_sf` Africa polygons, cleaned the way the data-overview figure
#' needs them: Yemen dropped, cropped to the continent bbox, AFRO region
#' attached (factor, taxdat level order), geometry repaired.
load_taxdat_afr_sf <- function(
  drop_iso3 = "YEM",
  bbox = sf::st_bbox(sf::st_sfc(
    sf::st_point(c(-18.8, 36.6)),
    sf::st_point(c(52.2, -36.5)),
    crs = 4326
  ))
) {
  e <- new.env()
  utils::data("afr_sf", package = "taxdat", envir = e)
  e$afr_sf %>%
    janitor::clean_names() %>%
    dplyr::mutate(country_name = country) %>%
    dplyr::filter(!country_name %in% drop_iso3) %>%
    sf::st_crop(bbox) %>%
    taxdat::get_AFRO_region(ctry_col = "country") %>%
    dplyr::mutate(
      AFRO_region = factor(
        AFRO_region,
        levels = taxdat::get_AFRO_region_levels()
      )
    ) %>%
    sf::st_make_valid()
}

#' Lakes and rivers clipped to a continent geometry, for the taxdat basemap
#' (`output_plot_map(lakes_sf=, rivers_sf=)`). `clip` is any sf/sfc;
#' `predicate` is the spatial predicate for the clip (st_intersects keeps
#' anything touching Africa, st_within only what lies fully inside).
#' @return list(lakes_sf, rivers_sf)
load_africa_basemap <- function(
  clip,
  lakes_path = "data/geodata/Africa_waterbody.shp",
  rivers_path = "data/geodata/rivers_africa_37333.shp",
  predicate = sf::st_intersects
) {
  list(
    lakes_sf = taxdat::get_lakes(path = lakes_path) %>%
      sf::st_make_valid() %>%
      sf::st_filter(clip, .predicate = predicate),
    rivers_sf = taxdat::get_rivers(path = rivers_path) %>%
      sf::st_make_valid() %>%
      sf::st_filter(clip, .predicate = predicate)
  )
}

# ---------------------------------------------------------------------------
# ps_diverging_fill_scale (Figures 1-4)
# ---------------------------------------------------------------------------
#' The diverging map colourbar shared by the manuscript maps (Figure 2's): a
#' gradient2 low / white / high around `midpoint`, values beyond `limits`
#' squished onto the end colours, a horizontal key_width_mm x key_height_mm bar
#' with its title on the left. Pass literal end labels (e.g. "≤ 1/3",
#' "≥ 3", with true minus signs); `limits` / `breaks` on the scale of the
#' mapped values (e.g. log10 for ratios).
ps_diverging_fill_scale <- function(
  limits,
  breaks,
  labels,
  name,
  midpoint = 0,
  low = "#2166ac",
  mid = "white",
  high = "#d73027",
  key_width_mm = 32,
  key_height_mm = 2,
  title_position = "left",
  order = 1,
  aesthetics = "fill"
) {
  ggplot2::scale_fill_gradient2(
    midpoint = midpoint,
    low = low,
    mid = mid,
    high = high,
    limits = limits,
    oob = scales::squish,
    breaks = breaks,
    labels = labels,
    name = name,
    aesthetics = aesthetics,
    guide = ggplot2::guide_colourbar(
      direction = "horizontal",
      title.position = title_position,
      order = order,
      theme = ggplot2::theme(
        legend.key.width = grid::unit(key_width_mm, "mm"),
        legend.key.height = grid::unit(key_height_mm, "mm")
      )
    )
  )
}
