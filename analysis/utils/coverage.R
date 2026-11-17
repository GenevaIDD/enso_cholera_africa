# coverage.R — the project's surveillance-coverage definition (which gid-months
#   were under surveillance), read by occurrence prep, 07c, 08c and 09.
#
# Provides:
#   coverage_lp_bridge()    LP -> gid bridge (full member sidecar if present)
#   build_coverage_mask()   gid x month coverage mask
#
# Details:
#   national_zero_path is REQUIRED while national_zero is on (a default
#   threshold-0 file would silently mismatch a thresholded run).
#   Sources nothing. Call-time: dplyr / tibble / lubridate via ::, `%>%`,
#   paths.R (occ_panel_lp_members_path()) for coverage_lp_bridge(surv_path =).

# Why the bridge: under gid_season the survival panel carries ONE
# representative location_period_id per gid, so its bridge misses the other
# member LPs' reporting windows (under-coverage). 03_prepare_occ.R writes the
# full member bridge as the sidecar <panel>.lp_members.csv
# (paths.R::occ_panel_lp_members_path()); callers use it when it exists and
# fall back to the panel's representative LP when it does not, so artefacts of
# a panel built without the sidecar are unchanged.
#' LP -> gid bridge for build_coverage_mask(): the panel's full member bridge
#' sidecar when it exists, else the panel itself (one representative LP per gid
#' under gid_season).
#'
#' @param obs_panel the survival panel (or any frame with location_period_id,
#'   gid, country).
#' @param surv_path the panel's path; its `<panel>.lp_members.csv` sidecar is
#'   read when present. NULL = no sidecar lookup.
#' @return a frame with location_period_id, gid, country: the sidecar rows of
#'   the panel's gids, or `obs_panel` unchanged.
coverage_lp_bridge <- function(obs_panel, surv_path = NULL) {
  if (is.null(surv_path) || !nzchar(surv_path)) {
    return(obs_panel)
  }
  mp <- occ_panel_lp_members_path(surv_path)
  if (!file.exists(mp)) {
    return(obs_panel)
  }
  mem <- utils::read.csv(mp, colClasses = "character", stringsAsFactors = FALSE)
  need <- c("location_period_id", "gid", "country")
  if (!all(need %in% names(mem))) {
    stop(
      "coverage_lp_bridge: ",
      mp,
      " lacks column(s) ",
      paste(setdiff(need, names(mem)), collapse = ", "),
      call. = FALSE
    )
  }
  gids <- unique(as.character(obs_panel$gid))
  out <- mem[mem$gid %in% gids, need, drop = FALSE]
  message(sprintf(
    "[coverage] LP bridge from %s: %d LP-gid pair(s) over %d gid(s) (panel: %d representative LP(s)).",
    basename(mp),
    nrow(out),
    length(unique(out$gid)),
    length(unique(as.character(obs_panel$location_period_id)))
  ))
  tibble::as_tibble(out)
}

# ---------------------------------------------------------------------------
# build_coverage_mask — (gid, month) cells with surveillance coverage
# ---------------------------------------------------------------------------
# Canonical, single-source-of-truth definition, called by
# 07c_simulator_postprocess.R and utils/smoothing_core.R.

#' Build a coverage mask of (gid, month) cells with observed surveillance data.
#'
#' Reads the ADM2-level coverage windows from the pre-outbreak extraction time
#' series, joins to simulated gids via the survival panel (which carries both
#' `location_period_id` and `gid`), and expands each TL–TR window into monthly
#' cells. For gids absent from their own ADM2 series, falls back to their
#' country×ADM2 window from `aggregated_admin_level_ts_preoutbreak_extraction.rds`.
#'
#' @param obs_panel data frame from the survival panel RDS — must carry
#'   `location_period_id` and `gid`.
#' @param aggregated_ts_path path to `aggregated_ts_preoutbreak_extraction.rds`.
#' @param aggregated_adm_ts_path path to
#'   `aggregated_admin_level_ts_preoutbreak_extraction.rds` (fallback).
#' @param include_adm1 logical (default `FALSE`, preserving legacy behaviour).
#'   When `TRUE`, ADM1-level extraction windows are also honoured, so a gid whose
#'   only pre-outbreak series is ADM1 uses **its own ADM1 TL–TR window** instead
#'   of being pushed to the coarse country×ADM2 fallback. This eliminates the
#'   dominant coverage-mask leaker: ~112 ADM1 gids (concentrated in
#'   coarse-reporting countries — UGA/NGA/SOM/COD/CMR) that otherwise receive a
#'   whole-country union of ADM2 windows. Kept as a flag so the effect on the
#'   total-case calibration can be measured against the legacy mask.
#' @param national_zero logical (default `TRUE`): also cover every gid of a
#'   country in that country's national-zero months -- months in which national
#'   and WHO reporting put the whole country at 0 suspected cases and no
#'   sub-national report has cases (utils/national_zero.R, built by
#'   01_preprocess_national_zero_months.R). `FALSE` gives the window-only mask.
#' @param national_zero_path path to the national-zero months
#'   (`make_national_zero_months_file(<threshold>[, <rule>])` of the run's panel);
#'   REQUIRED while `national_zero` is on (no default: a threshold-0 file would
#'   silently mismatch a thresholded run); a missing file is an error.
#' @return tibble(gid, month) of covered cells; month is a Date floored to the
#'   first of the month.
build_coverage_mask <- function(
  obs_panel,
  aggregated_ts_path,
  aggregated_adm_ts_path = NULL,
  include_adm1 = FALSE,
  national_zero = TRUE,
  national_zero_path
) {
  # Bridge: location_period_id ↔ gid from the survival panel
  lp_bridge <- obs_panel %>%
    dplyr::filter(!is.na(location_period_id) & !is.na(gid)) %>%
    dplyr::distinct(
      location_period_id = as.character(location_period_id),
      gid = as.character(gid),
      country = as.character(country)
    )

  if (!file.exists(aggregated_ts_path)) {
    stop("build_coverage_mask: aggregated_ts not found: ", aggregated_ts_path)
  }

  keep_levels <- if (isTRUE(include_adm1)) {
    c("ADM1", "ADM2", "ADM3")
  } else {
    c("ADM2", "ADM3")
  }

  agg_ts <- readRDS(aggregated_ts_path) %>%
    dplyr::filter(admin_level %in% keep_levels) %>%
    dplyr::mutate(
      location_period_id = as.character(location_period_id),
      TL = as.Date(TL),
      TR = as.Date(TR)
    ) %>%
    dplyr::select(-dplyr::any_of("country")) %>% # avoid country.x/country.y on join
    dplyr::inner_join(lp_bridge, by = "location_period_id") %>%
    dplyr::select(gid, country, TL, TR)

  # Fallback: gids without their own ADM2 coverage series → country×ADM2 window
  covered_gids <- unique(agg_ts$gid)
  all_gids <- unique(lp_bridge$gid)
  missing_gids <- setdiff(all_gids, covered_gids)

  if (
    length(missing_gids) > 0 &&
      !is.null(aggregated_adm_ts_path) &&
      file.exists(aggregated_adm_ts_path)
  ) {
    adm_ts <- readRDS(aggregated_adm_ts_path) %>%
      dplyr::filter(admin_level == "ADM2") %>%
      dplyr::mutate(
        country = as.character(country),
        TL = as.Date(TL),
        TR = as.Date(TR)
      ) %>%
      dplyr::select(country, TL, TR)

    miss_meta <- lp_bridge %>%
      dplyr::filter(gid %in% missing_gids) %>%
      dplyr::distinct(gid, country)

    fallback <- dplyr::inner_join(
      miss_meta,
      adm_ts,
      by = "country",
      relationship = "many-to-many"
    )
    message(sprintf(
      "[build_coverage_mask] %d gid(s) used country×ADM2 fallback coverage window.",
      length(missing_gids)
    ))
    agg_ts <- dplyr::bind_rows(agg_ts, fallback)
  } else if (length(missing_gids) > 0) {
    message(sprintf(
      "[build_coverage_mask] %d gid(s) have no coverage window — treated as uncovered.",
      length(missing_gids)
    ))
  }

  # Expand each (gid, TL, TR) window into monthly cells
  mask_list <- vector("list", nrow(agg_ts))
  for (i in seq_len(nrow(agg_ts))) {
    tl <- lubridate::floor_date(agg_ts$TL[i], "month")
    tr <- lubridate::floor_date(agg_ts$TR[i], "month")
    if (is.na(tl) || is.na(tr) || tl > tr) {
      next
    }
    months_seq <- seq(tl, tr, by = "month")
    mask_list[[i]] <- tibble::tibble(gid = agg_ts$gid[i], month = months_seq)
  }

  mask <- dplyr::bind_rows(mask_list) %>%
    dplyr::distinct(gid, month)

  # National-zero months: the whole country observed at 0, so every panel gid
  # of the country is covered that month.
  if (isTRUE(national_zero)) {
    if (missing(national_zero_path) || is.null(national_zero_path)) {
      stop(
        "build_coverage_mask: national_zero_path is required while national_zero is on ",
        "-- pass the run's make_national_zero_months_file(<threshold>[, <rule>]), or ",
        "national_zero = FALSE.",
        call. = FALSE
      )
    }
    if (!file.exists(national_zero_path)) {
      stop(
        "build_coverage_mask: national-zero months not found: ",
        national_zero_path,
        ". Build them with\n  Rscript analysis/01_preprocess_national_zero_months.R\n",
        "or pass national_zero = FALSE.",
        call. = FALSE
      )
    }
    nz <- readRDS(national_zero_path) %>%
      dplyr::transmute(country = as.character(country), month = as.Date(month))
    nz_cells <- lp_bridge %>%
      dplyr::distinct(gid, country) %>%
      dplyr::inner_join(nz, by = "country", relationship = "many-to-many") %>%
      dplyr::distinct(gid, month)
    n_before <- nrow(mask)
    mask <- dplyr::bind_rows(mask, nz_cells) %>% dplyr::distinct(gid, month)
    message(sprintf(
      "[build_coverage_mask] national-zero months: %d cell(s), %d new to the mask.",
      nrow(nz_cells),
      nrow(mask) - n_before
    ))
  }
  mask
}
