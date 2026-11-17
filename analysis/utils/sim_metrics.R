# sim_metrics.R — simulator validation and hindcast-skill metrics
#
# Provides:
#   load_observed_from_panel()     observed side (gid, month, Y_obs)
#   compute_pit_binom(), compute_monthly_freq(), compute_crps_annual(),
#   compute_annual_count_by_geo(), compute_ar_metrics()   occurrence
#   reconstruct_lag_series(), compute_contagion_conditionals()   contagion
#   compute_episode_case_totals(), attach_episode_case_totals(),
#   compute_crps_total_cases()     case totals (mon-track cases_month)
#   compute_duration_dist(), compute_active_fraction(),
#   compute_sim_kendall_tau()      duration
#   compute_hindcast_series()      hindcast
#
# Details:
#   07c_simulator_postprocess.R and 07d_simulator_figures.R are thin CLIs over
#   this library; orchestration (run_validation(), the hindcast loop) stays in
#   07c, partition IO in sim_io.R. Callers source base.R, gid.R,
#   obs_sim_compare.R (crps_ens(), region_key_of()), sim_io.R,
#   simulator_core.R (update_lags()) and dependence.R (paired_kendall_tau()).
#   Gate: tests/test_sim_metrics_golden.R.

# ---------------------------------------------------------------------------
# load_observed_from_panel — collapse the survival panel to (gid, month, Y_obs)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Shared geography helpers (used by compute_annual_count_by_geo() and
# compute_crps_total_cases()).
# ---------------------------------------------------------------------------

#' .sim_frame_geo
#'
#' gid -> country / AFRO_region / africa geography tags from a location frame.
#' @noRd
.sim_frame_geo <- function(frame) {
  frame %>%
    dplyr::distinct(
      gid = as.character(gid),
      country = as.character(country),
      AFRO_region = as.character(AFRO_region)
    ) %>%
    dplyr::mutate(africa = "Africa")
}

#' .geo_crps_level
#'
#' Score one geography level: optional coverage pro-rating of BOTH sides by the
#' level's covered gid-month fraction, the ensemble CRPS (crps_ens(),
#' utils/obs_sim_compare.R) plus sim mean/quantiles per (geography, year), the
#' unmasked sim mean, and the covered gid-month count. Callers build the
#' per-replicate simulated totals (`sim_ann`: replicate_id, geography, year,
#' total_sim), the observed totals (`obs_ann`: geography, year, total_obs) and
#' the unmasked means (`sim_ann_u`) their own way -- annual outbreak counts for
#' compute_annual_count_by_geo(), episode case totals for
#' compute_crps_total_cases() -- and name the CRPS column (`crps_col`).
#' @noRd
.geo_crps_level <- function(
  sim_ann,
  obs_ann,
  sim_ann_u,
  covered_frac,
  frame_geo,
  geo_col,
  crps_col
) {
  # Pro-rate by coverage fraction at the aggregated level
  frac_geo <- NULL
  if (!is.null(covered_frac)) {
    geo_map <- frame_geo %>% dplyr::distinct(gid, geography = .data[[geo_col]])
    n_gids_geo <- geo_map %>% dplyr::count(geography, name = "n_gids")
    frac_geo <- covered_frac %>%
      dplyr::left_join(geo_map, by = "gid") %>%
      dplyr::filter(!is.na(geography)) %>%
      dplyr::group_by(geography, year) %>%
      dplyr::summarise(
        n_covered_gid_months = sum(n_covered_months),
        .groups = "drop"
      ) %>%
      dplyr::left_join(n_gids_geo, by = "geography") %>%
      dplyr::mutate(proratio = (n_gids * 12L) / n_covered_gid_months)

    sim_ann <- sim_ann %>%
      dplyr::left_join(
        frac_geo %>% dplyr::select(geography, year, proratio),
        by = c("geography", "year")
      ) %>%
      dplyr::mutate(total_sim = total_sim * dplyr::coalesce(proratio, 1)) %>%
      dplyr::select(-proratio)
    obs_ann <- obs_ann %>%
      dplyr::left_join(
        frac_geo %>% dplyr::select(geography, year, proratio),
        by = c("geography", "year")
      ) %>%
      dplyr::mutate(total_obs = total_obs * dplyr::coalesce(proratio, 1)) %>%
      dplyr::select(-proratio)
  }

  joined <- dplyr::inner_join(sim_ann, obs_ann, by = c("geography", "year"))
  result <- joined %>%
    dplyr::group_by(geography, year) %>%
    dplyr::summarise(
      n_sim = dplyr::n(),
      total_obs = dplyr::first(total_obs),
      !!crps_col := crps_ens(total_sim, dplyr::first(total_obs)),
      total_sim_mean = mean(total_sim, na.rm = TRUE),
      total_sim_q25 = stats::quantile(total_sim, 0.25, na.rm = TRUE),
      total_sim_q75 = stats::quantile(total_sim, 0.75, na.rm = TRUE),
      total_sim_q025 = stats::quantile(total_sim, 0.025, na.rm = TRUE),
      total_sim_q975 = stats::quantile(total_sim, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(sim_ann_u, by = c("geography", "year"))

  if (!is.null(frac_geo)) {
    result <- dplyr::left_join(
      result,
      frac_geo %>% dplyr::select(geography, year, n_covered_gid_months),
      by = c("geography", "year")
    )
  } else {
    result$n_covered_gid_months <- NA_integer_
  }
  result
}

#' Load and collapse the observed outbreak panel.
#'
#' Requires columns `gid`, `t_start`, `outbreak` in `panel`. Groups by
#' (gid, month) and takes the mean `outbreak` value per cell, collapsing
#' multi-row months (multiple location-period ids) to the fraction of active
#' lp_ids in outbreak. This is directly comparable to the per-(gid,month)
#' predicted probability from the simulator.
#'
#' @param panel data frame with at minimum gid, t_start, outbreak columns.
#' @param gids optional character vector of gids to retain (NULL = all).
#' @param horizon optional list(start, end) of Date-coercible strings to filter
#'   the month dimension (inclusive on both ends).
#' @return tibble(gid, month, Y_obs) with month as Date; Y_obs in [0, 1].
load_observed_from_panel <- function(panel, gids = NULL, horizon = NULL) {
  required <- c("gid", "t_start", "outbreak")
  missing <- setdiff(required, names(panel))
  if (length(missing) > 0) {
    stop(
      "load_observed_from_panel: panel missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  out <- panel %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = lubridate::floor_date(as.Date(t_start), "month")
    ) %>%
    dplyr::group_by(gid, month) %>%
    dplyr::summarise(Y_obs = mean(outbreak, na.rm = TRUE), .groups = "drop")

  if (!is.null(gids)) {
    out <- out %>% dplyr::filter(gid %in% as.character(gids))
  }

  if (!is.null(horizon)) {
    h_start <- as.Date(horizon$start)
    h_end <- as.Date(horizon$end)
    out <- out %>%
      dplyr::filter(month >= h_start, month <= h_end)
  }

  out
}

# ---------------------------------------------------------------------------
# build_coverage_mask — (gid, month) cells with surveillance coverage
# ---------------------------------------------------------------------------
# Canonical definition lives in utils/sim_io.R; sourced via the
# `source("analysis/utils/sim_io.R")` call in this file's CLI-entrypoint guard.

# ---------------------------------------------------------------------------
# compute_pit_binom — randomised PIT for binary outcomes
# ---------------------------------------------------------------------------

#' Compute the randomised probability integral transform for binary outcomes.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y) — simulated monthly
#'   binary outcomes.
#' @param obs_monthly tibble(gid, month, Y_obs) with month as Date.
#' @param seed integer seed for the runif calls (default 42L).
#' @return tibble(gid, month, p_hat, y_obs, pit).
compute_pit_binom <- function(sim_monthly, obs_monthly, seed = 42L) {
  # Compute p_hat = mean(Y) per (gid, month) across replicates
  p_hat_tbl <- sim_monthly %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = as.Date(month)
    ) %>%
    dplyr::group_by(gid, month) %>%
    dplyr::summarise(p_hat = mean(Y, na.rm = TRUE), .groups = "drop")

  obs_tbl <- obs_monthly %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = as.Date(month)
    )

  joined <- dplyr::inner_join(p_hat_tbl, obs_tbl, by = c("gid", "month"))

  if (nrow(joined) == 0L) {
    warning(
      "compute_pit_binom: inner join on (gid, month) produced 0 rows — check alignment."
    )
    return(tibble::tibble(
      gid = character(),
      month = as.Date(character()),
      p_hat = numeric(),
      y_obs = integer(),
      pit = numeric()
    ))
  }

  n <- nrow(joined)
  # Binomial PIT assumes binary occurrence; coerce obs to {0,1} once so the
  # branch test and the returned y_obs column agree (a fractional Y_obs would
  # otherwise hit the non-zero branch yet truncate to 0 via as.integer()).
  y_obs_int <- as.integer(joined$Y_obs > 0)
  set.seed(seed)
  pit <- ifelse(
    y_obs_int == 0,
    stats::runif(n, 0, 1 - joined$p_hat),
    stats::runif(n, 1 - joined$p_hat, 1)
  )

  tibble::tibble(
    gid = joined$gid,
    month = joined$month,
    p_hat = joined$p_hat,
    y_obs = y_obs_int,
    pit = pit
  )
}

# ---------------------------------------------------------------------------
# compute_monthly_freq — seasonal frequency comparison
# ---------------------------------------------------------------------------

#' Compare observed vs simulated monthly outbreak frequency by country.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y).
#' @param obs_monthly tibble(gid, month, Y_obs) with month as Date.
#' @param frame tibble with columns gid and country.
#' @param mask optional tibble(gid, month) from build_coverage_mask(). When
#'   supplied, both obs and sim are restricted to covered (gid, month) cells
#'   before computing monthly frequencies.
#' @return tibble(country, month_of_year, freq_obs, freq_sim_mean,
#'   freq_sim_lo, freq_sim_hi).
compute_monthly_freq <- function(sim_monthly, obs_monthly, frame, mask = NULL) {
  ctry_map <- frame %>%
    dplyr::distinct(gid = as.character(gid), country = as.character(country))

  obs_in <- obs_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))
  sim_in <- sim_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))

  if (!is.null(mask)) {
    mask_norm <- mask %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month))
    obs_in <- dplyr::semi_join(obs_in, mask_norm, by = c("gid", "month"))
    sim_in <- dplyr::semi_join(sim_in, mask_norm, by = c("gid", "month"))
  }

  obs_aug <- obs_in %>%
    dplyr::mutate(month_of_year = lubridate::month(month)) %>%
    dplyr::inner_join(ctry_map, by = "gid") %>%
    dplyr::group_by(country, month_of_year) %>%
    dplyr::summarise(freq_obs = mean(Y_obs, na.rm = TRUE), .groups = "drop")

  sim_aug <- sim_in %>%
    dplyr::mutate(month_of_year = lubridate::month(month)) %>%
    dplyr::inner_join(ctry_map, by = "gid") %>%
    dplyr::group_by(replicate_id, country, month_of_year) %>%
    dplyr::summarise(freq_sim = mean(Y, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(country, month_of_year) %>%
    dplyr::summarise(
      freq_sim_mean = mean(freq_sim, na.rm = TRUE),
      freq_sim_q25 = stats::quantile(freq_sim, 0.25, na.rm = TRUE),
      freq_sim_q75 = stats::quantile(freq_sim, 0.75, na.rm = TRUE),
      freq_sim_q025 = stats::quantile(freq_sim, 0.025, na.rm = TRUE),
      freq_sim_q975 = stats::quantile(freq_sim, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  dplyr::full_join(obs_aug, sim_aug, by = c("country", "month_of_year"))
}

# ---------------------------------------------------------------------------
# compute_crps_annual — CRPS for annual counts
# ---------------------------------------------------------------------------

#' Compute the continuous ranked probability score for annual outbreak counts.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y).
#' @param obs_monthly tibble(gid, month, Y_obs) with month as Date.
#' @param mask optional tibble(gid, month) from build_coverage_mask(). When
#'   supplied, both series are restricted to covered cells and annual counts are
#'   pro-rated to a full year by multiplying by 12 / n_covered_months.
#' @return tibble(gid, year, n_sim, y_obs, crps_count, n_covered_months).
compute_crps_annual <- function(sim_monthly, obs_monthly, mask = NULL) {
  obs_in <- obs_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))
  sim_in <- sim_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))

  if (!is.null(mask)) {
    mask_norm <- mask %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month))
    obs_in <- dplyr::semi_join(obs_in, mask_norm, by = c("gid", "month"))
    sim_in <- dplyr::semi_join(sim_in, mask_norm, by = c("gid", "month"))
  }

  # Covered months per (gid, year) for pro-rating
  covered_frac <- if (!is.null(mask)) {
    mask %>%
      dplyr::mutate(
        gid = as.character(gid),
        month = as.Date(month),
        year = lubridate::year(month)
      ) %>%
      dplyr::group_by(gid, year) %>%
      dplyr::summarise(n_covered_months = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(proratio = 12 / n_covered_months)
  } else {
    NULL
  }

  sim_ann <- sim_in %>%
    dplyr::mutate(year = lubridate::year(month)) %>%
    dplyr::group_by(replicate_id, gid, year) %>%
    dplyr::summarise(Y_ann = sum(Y, na.rm = TRUE), .groups = "drop")

  obs_ann <- obs_in %>%
    dplyr::mutate(year = lubridate::year(month)) %>%
    dplyr::group_by(gid, year) %>%
    dplyr::summarise(Y_ann_obs = sum(Y_obs, na.rm = TRUE), .groups = "drop")

  if (!is.null(covered_frac)) {
    sim_ann <- dplyr::left_join(
      sim_ann,
      covered_frac,
      by = c("gid", "year")
    ) %>%
      dplyr::mutate(Y_ann = Y_ann * dplyr::coalesce(proratio, 1)) %>%
      dplyr::select(-n_covered_months, -proratio)
    obs_ann <- dplyr::left_join(
      obs_ann,
      covered_frac,
      by = c("gid", "year")
    ) %>%
      dplyr::mutate(Y_ann_obs = Y_ann_obs * dplyr::coalesce(proratio, 1))
  }

  joined <- dplyr::inner_join(sim_ann, obs_ann, by = c("gid", "year"))

  result <- joined %>%
    dplyr::group_by(gid, year) %>%
    dplyr::summarise(
      n_sim = dplyr::n(),
      y_obs = dplyr::first(Y_ann_obs),
      crps_count = crps_ens(Y_ann, dplyr::first(Y_ann_obs)), # utils/obs_sim_compare.R
      .groups = "drop"
    )

  if (!is.null(covered_frac)) {
    result <- dplyr::left_join(
      result,
      covered_frac %>% dplyr::select(gid, year, n_covered_months),
      by = c("gid", "year")
    )
  } else {
    result$n_covered_months <- NA_integer_
  }

  result
}

# ---------------------------------------------------------------------------
# compute_annual_count_by_geo — annual outbreak-count CRPS at three levels
# ---------------------------------------------------------------------------

#' Compute annual outbreak-count summaries (sim ribbon + obs + CRPS) at the
#' Africa / AFRO-region / country level.
#'
#' Mirror of `compute_crps_total_cases` but for outbreak OCCURRENCE COUNTS
#' (binary monthly Y_obs from the survival panel) rather than case totals.
#' Input is `sim_monthly` / `obs_monthly` — no event table required.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y) — monthly sim.
#' @param obs_monthly tibble(gid, month, Y_obs) — monthly obs (binary).
#' @param frame tibble from build_location_frame(); must carry gid, country,
#'   AFRO_region.
#' @param mask optional tibble(gid, month) from build_coverage_mask().
#' @return named list(africa, region, country) each a tibble with columns
#'   geography, year, n_sim, total_obs, crps_count, total_sim_mean,
#'   total_sim_lo, total_sim_hi, total_sim_unmasked_mean, n_covered_gid_months.
compute_annual_count_by_geo <- function(
  sim_monthly,
  obs_monthly,
  frame,
  mask = NULL
) {
  frame_geo <- .sim_frame_geo(frame)

  sim_in <- sim_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))
  obs_in <- obs_monthly %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))

  # Every simulated total below is a sum of Y over (replicate, geography, year):
  # Y is summed per (replicate, gid, year) first (counts: exact in any order),
  # so the geography joins and groupings run on ~1/12 of the rows instead of
  # the replicate x gid x month grid (which would be 07c's memory peak,
  # ~10x sim_monthly). Same totals, same groups, same sort.
  per_gid_year <- function(df) {
    u <- unique(df$month)
    dt <- data.table::data.table(
      replicate_id = df$replicate_id,
      gid = df$gid,
      year = lubridate::year(u)[match(df$month, u)],
      Y = df$Y
    )
    tibble::as_tibble(as.data.frame(dt[,
      list(Y = sum(Y, na.rm = TRUE)),
      by = c("replicate_id", "gid", "year")
    ]))
  }

  # Aggregate unmasked totals at each geographic level BEFORE applying the mask.
  # Avoids holding a masked and an unmasked full copy of sim_monthly simultaneously;
  # peak usage drops from ~4x to ~2x sim_monthly size.
  {
    sim_um <- per_gid_year(sim_in) %>%
      dplyr::left_join(frame_geo, by = "gid")
    unmasked_means <- stats::setNames(
      lapply(c("africa", "AFRO_region", "country"), function(gc) {
        sim_um %>%
          dplyr::filter(!is.na(.data[[gc]])) %>%
          dplyr::group_by(replicate_id, geography = .data[[gc]], year) %>%
          dplyr::summarise(total_u = sum(Y, na.rm = TRUE), .groups = "drop") %>%
          dplyr::group_by(geography, year) %>%
          dplyr::summarise(
            total_sim_unmasked_mean = mean(total_u, na.rm = TRUE),
            .groups = "drop"
          )
      }),
      c("africa", "AFRO_region", "country")
    )
    rm(sim_um)
  }
  gc()

  covered_frac <- NULL
  if (!is.null(mask)) {
    mask_norm <- mask %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month))
    obs_in <- dplyr::semi_join(obs_in, mask_norm, by = c("gid", "month"))
    sim_in <- dplyr::semi_join(sim_in, mask_norm, by = c("gid", "month"))
    covered_frac <- mask_norm %>%
      dplyr::mutate(year = lubridate::year(month)) %>%
      dplyr::group_by(gid, year) %>%
      dplyr::summarise(n_covered_months = dplyr::n(), .groups = "drop")
  }

  # Tag with geography before aggregating (masked copy only)
  sim_geo_tagged <- per_gid_year(sim_in) %>%
    dplyr::left_join(frame_geo, by = "gid")
  rm(sim_in)

  obs_geo_tagged <- obs_in %>%
    dplyr::mutate(year = lubridate::year(month)) %>%
    dplyr::left_join(frame_geo, by = "gid")

  crps_at_level <- function(geo_col) {
    sim_ann <- sim_geo_tagged %>%
      dplyr::filter(!is.na(.data[[geo_col]])) %>%
      dplyr::group_by(replicate_id, geography = .data[[geo_col]], year) %>%
      dplyr::summarise(total_sim = sum(Y, na.rm = TRUE), .groups = "drop")

    sim_ann_u <- unmasked_means[[geo_col]]

    obs_ann <- obs_geo_tagged %>%
      dplyr::filter(!is.na(.data[[geo_col]])) %>%
      dplyr::group_by(geography = .data[[geo_col]], year) %>%
      dplyr::summarise(total_obs = sum(Y_obs, na.rm = TRUE), .groups = "drop")

    .geo_crps_level(
      sim_ann,
      obs_ann,
      sim_ann_u,
      covered_frac,
      frame_geo,
      geo_col,
      crps_col = "crps_count"
    )
  }

  list(
    africa = crps_at_level("africa"),
    region = crps_at_level("AFRO_region"),
    country = crps_at_level("country")
  )
}

# ---------------------------------------------------------------------------
# compute_ar_metrics — run/gap lengths and transition matrix
# ---------------------------------------------------------------------------

#' Compute AR (autoregressive) sequence diagnostics from simulated monthly Y.
#'
#' Uses `rle()` to characterise run and gap lengths, and counts 1-step
#' transitions across the pooled series.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y).
#' @return list(run_length_table, gap_length_table, transition).
compute_ar_metrics <- function(sim_monthly) {
  series_list <- sim_monthly %>%
    dplyr::arrange(replicate_id, gid, month) %>%
    dplyr::group_by(replicate_id, gid) %>%
    dplyr::summarise(Y_seq = list(Y), .groups = "drop") %>%
    dplyr::pull(Y_seq)

  # Collect run/gap lengths as a list → single unlist() at the end.
  # Avoids O(n²) memory from repeated c() growth in a loop over ~498K series.
  run_lens_list <- vector("list", length(series_list))
  gap_lens_list <- vector("list", length(series_list))

  # 2×2 transition count matrix [from+1, to+1] (column-major storage).
  # tabulate() accumulates counts in-place (a do.call(rbind) over per-replicate
  # data frames would be slow and memory-intensive).
  trans_counts <- matrix(0L, nrow = 2L, ncol = 2L)

  for (i in seq_along(series_list)) {
    y <- as.integer(series_list[[i]])
    if (length(y) < 2L) {
      next
    }

    r <- rle(y)
    run_lens_list[[i]] <- r$lengths[r$values == 1L]
    gap_lens_list[[i]] <- r$lengths[r$values == 0L]

    n_t <- length(y) - 1L
    fv <- y[seq_len(n_t)]
    tv <- y[seq_len(n_t) + 1L]
    ok <- !is.na(fv) & !is.na(tv) & fv %in% 0:1 & tv %in% 0:1
    if (any(ok)) {
      # column-major index: (from,to) ∈ {0,1}² → 1..4
      idx <- fv[ok] + 1L + tv[ok] * 2L
      tc <- tabulate(idx, nbins = 4L)
      trans_counts[1L] <- trans_counts[1L] + tc[1L]
      trans_counts[2L] <- trans_counts[2L] + tc[2L]
      trans_counts[3L] <- trans_counts[3L] + tc[3L]
      trans_counts[4L] <- trans_counts[4L] + tc[4L]
    }
  }

  all_run_lens <- unlist(run_lens_list)
  all_gap_lens <- unlist(gap_lens_list)

  if (sum(trans_counts) == 0L) {
    transition <- tibble::tibble(
      from = integer(),
      to = integer(),
      n = integer(),
      p = numeric()
    )
  } else {
    transition <- tibble::tibble(
      from = c(0L, 0L, 1L, 1L),
      to = c(0L, 1L, 0L, 1L),
      n = as.integer(c(
        trans_counts[1L, 1L],
        trans_counts[1L, 2L],
        trans_counts[2L, 1L],
        trans_counts[2L, 2L]
      ))
    ) %>%
      dplyr::filter(n > 0L) %>%
      dplyr::group_by(from) %>%
      dplyr::mutate(p = n / sum(n)) %>%
      dplyr::ungroup()
  }

  list(
    run_length_table = table(all_run_lens),
    gap_length_table = table(all_gap_lens),
    transition = transition
  )
}

# ---------------------------------------------------------------------------
# reconstruct_lag_series — replay the lag update over a simulated trajectory
# ---------------------------------------------------------------------------

#' Reconstruct the lag time-series (o_adm1, o_ct, o_tu) from a monthly summary.
#'
#' Applies `update_lags()` step-by-step to recover what the lag variables would
#' have been at each month, starting from zeros.
#'
#' `Y` (onset) is always the outcome column carried through unchanged, so all
#' downstream onset-based metrics (contagion conditionals' P(Y|lag), PIT, CRPS,
#' AR run-lengths, ...) keep meaning what they always meant. `active_col` only
#' controls what is FED INTO `update_lags()` to build the lag state itself:
#' legacy/default (`NULL`) replays `update_lags()` on the onset pulse `Y`.
#' When a duration-aware run
#' supplies `active_col = "active"`, the lag state is instead replayed on the
#' persistent active-location indicator — matching what `simulate_one_tu()`
#' actually feeds `update_lags()` once a duration model is attached, rather than
#' re-deriving a stale onset-only proxy.
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y[, active]).
#' @param frame tibble from build_location_frame() with gid, adm1_unit, country, tu.
#' @param max_reps integer; restrict to the first `max_reps` replicates by ID.
#' @param active_col optional column name in `sim_monthly` (e.g. "active")
#'   whose values should be fed into `update_lags()` instead of `Y`. Default
#'   `NULL` reproduces the legacy (onset-fed) behaviour unchanged. If the
#'   column is requested but absent from `sim_monthly` (e.g. a legacy
#'   non-duration run's monthly_summary, which never grows an `active`
#'   column), falls back to the legacy Y-fed behaviour with a message.
#' @return long tibble(replicate_id, gid, month, Y, o_adm1, o_ct, o_tu).
reconstruct_lag_series <- function(
  sim_monthly,
  frame,
  max_reps = NULL,
  active_col = NULL
) {
  has_active <- !is.null(active_col)
  if (has_active && !(active_col %in% names(sim_monthly))) {
    message(
      sprintf(
        "reconstruct_lag_series: active_col '%s' not found in sim_monthly -- ",
        active_col
      ),
      "falling back to legacy Y-fed lag reconstruction."
    )
    has_active <- FALSE
  }

  rep_ids <- sort(unique(sim_monthly$replicate_id))
  if (!is.null(max_reps)) {
    rep_ids <- head(rep_ids, max_reps)
  }

  months <- sort(unique(as.Date(sim_monthly$month)))
  M <- length(months)
  N <- nrow(frame)

  result_list <- vector("list", length(rep_ids))

  for (r_idx in seq_along(rep_ids)) {
    rep_id <- rep_ids[[r_idx]]

    sub <- sim_monthly %>%
      dplyr::filter(replicate_id == rep_id) %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month))

    # Build N × M integer matrices from the monthly summary: Y (onset, always
    # carried through) and, when available/requested, the persistent active
    # state that update_lags() is actually replayed on.
    Y_mat <- matrix(0L, nrow = N, ncol = M)
    State_mat <- if (has_active) matrix(0L, nrow = N, ncol = M) else NULL
    for (m_idx in seq_len(M)) {
      month_m <- months[[m_idx]]
      sub_m <- sub %>% dplyr::filter(month == month_m)
      gid_idx <- match(sub_m$gid, frame$gid)
      ok <- !is.na(gid_idx)
      Y_mat[gid_idx[ok], m_idx] <- as.integer(sub_m$Y[ok])
      if (has_active) {
        State_mat[gid_idx[ok], m_idx] <- as.integer(sub_m[[active_col]][ok])
      }
    }
    if (!has_active) {
      State_mat <- Y_mat
    }

    # Lag state matrices (record *before* update)
    o_mat <- matrix(0L, nrow = N, ncol = M)
    oct_mat <- matrix(0L, nrow = N, ncol = M)
    otu_mat <- matrix(NA_integer_, nrow = N, ncol = M)

    # Initialise
    o_prev <- rep(0L, N)
    oct_prev <- rep(0L, N)
    otu_prev <- ifelse(is.na(frame$tu), NA_integer_, 0L)

    for (m_idx in seq_len(M)) {
      o_mat[, m_idx] <- o_prev
      oct_mat[, m_idx] <- oct_prev
      otu_mat[, m_idx] <- otu_prev

      lags <- update_lags(State_mat[, m_idx], frame)
      o_prev <- lags$o
      oct_prev <- lags$o_ct
      otu_prev <- lags$o_tu
    }

    # Assemble long tibble for this replicate
    result_list[[r_idx]] <- tibble::tibble(
      replicate_id = rep_id,
      gid = rep(frame$gid, times = M),
      month = rep(months, each = N),
      Y = as.vector(Y_mat),
      o_adm1 = as.vector(o_mat),
      o_ct = as.vector(oct_mat),
      o_tu = as.vector(otu_mat)
    )
  }

  dplyr::bind_rows(result_list)
}

# ---------------------------------------------------------------------------
# compute_contagion_conditionals
# ---------------------------------------------------------------------------

#' Compute conditional outbreak probabilities given each lag type.
#'
#' Reconstructs the lag series and computes P(Y=1 | lag=1) and P(Y=1 | lag=0)
#' for each of the three lag types (ADM1, country, TU).
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y[, active]).
#' @param frame tibble from build_location_frame() with gid, adm1_unit, country, tu.
#' @param max_reps integer; passed to reconstruct_lag_series().
#' @param active_col optional column name (e.g. "active"); passed straight
#'   through to reconstruct_lag_series() -- see that function's doc comment.
#'   Default NULL (legacy Y-fed lag reconstruction).
#' @return list(adm1_cond, country_cond, tu_cond, unconditional).
compute_contagion_conditionals <- function(
  sim_monthly,
  frame,
  max_reps = 50L,
  active_col = NULL
) {
  lag_series <- reconstruct_lag_series(
    sim_monthly,
    frame,
    max_reps = max_reps,
    active_col = active_col
  )

  # Drop the first month per (replicate_id, gid) — warm-up lag is 0 by construction
  first_months <- lag_series %>%
    dplyr::group_by(replicate_id, gid) %>%
    dplyr::summarise(first_month = min(month), .groups = "drop")
  lag_series <- lag_series %>%
    dplyr::left_join(first_months, by = c("replicate_id", "gid")) %>%
    dplyr::filter(month != first_month) %>%
    dplyr::select(-first_month)

  uncond <- mean(lag_series$Y, na.rm = TRUE)

  cond_for_lag <- function(lag_col) {
    sub <- lag_series %>% dplyr::filter(!is.na(.data[[lag_col]]))
    list(
      lag_active = mean(sub$Y[sub[[lag_col]] == 1], na.rm = TRUE),
      lag_inactive = mean(sub$Y[sub[[lag_col]] == 0], na.rm = TRUE),
      unconditional = mean(sub$Y, na.rm = TRUE)
    )
  }

  list(
    adm1_cond = cond_for_lag("o_adm1"),
    country_cond = cond_for_lag("o_ct"),
    tu_cond = cond_for_lag("o_tu"),
    unconditional = uncond
  )
}

# ---------------------------------------------------------------------------
# compute_episode_case_totals — per-episode sums of the mon track's cases_month
# ---------------------------------------------------------------------------

#' Per-episode simulated case totals from the monthly summary.
#'
#' The size track (which drew one whole-episode total `size_cases` per onset) is
#' retired; the monthly (`mon`) track's per-active-month `cases_month` is the
#' sole case estimand. This sums `cases_month` over each simulated episode so
#' the episode-level metrics (total-case CRPS, simulated Kendall's tau) can be
#' keyed on the event table exactly as before.
#'
#' Episodes are delimited per (replicate_id, gid) by the onset indicator: after
#' ordering by month, `.ep = cumsum(Y)` numbers the episodes and the active
#' months (`active == 1`) carrying the same `.ep` belong to that episode. This
#' is exact because the kernel writes `active` after promotion (so the onset
#' month is active) and before cessation, `cases_month` is non-NA exactly where
#' `active == 1`, and `Y == 1` cannot occur while a unit is already active.
#' Active months before the first onset (`.ep == 0`, a warm-start episode
#' carried in from the handoff) have no event row and are excluded, as the
#' event table excludes them.
#'
#' The total is the IN-HORIZON sum: an episode still active at the horizon end
#' (censored) contributes only its simulated months, not a whole-outbreak total.
#'
#' @param sim_monthly data frame with replicate_id, gid, month, Y, active,
#'   cases_month (merged/monthly_summary of a run with duration + monthly
#'   enabled).
#' @return tibble(replicate_id, gid, start_month [character "YYYY-MM-DD"],
#'   episode_cases, n_months_cases), one row per episode; NULL when any of the
#'   required columns is missing (no-duration or no-monthly runs).
compute_episode_case_totals <- function(sim_monthly) {
  req <- c("replicate_id", "gid", "month", "Y", "active", "cases_month")
  if (is.null(sim_monthly) || !all(req %in% names(sim_monthly))) {
    return(NULL)
  }
  sim_monthly %>%
    dplyr::select(dplyr::all_of(req)) %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month)) %>%
    dplyr::arrange(replicate_id, gid, month) %>%
    dplyr::group_by(replicate_id, gid) %>%
    dplyr::mutate(.ep = cumsum(dplyr::coalesce(as.integer(Y), 0L))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(active == 1, .ep > 0L) %>%
    dplyr::group_by(replicate_id, gid, .ep) %>%
    dplyr::summarise(
      start_month = as.character(min(month)),
      episode_cases = sum(as.numeric(cases_month), na.rm = TRUE),
      n_months_cases = sum(!is.na(cases_month)),
      .groups = "drop"
    ) %>%
    dplyr::select(-.ep)
}

#' Left-join compute_episode_case_totals() output onto the event table.
#'
#' Key: one event row per onset (Y == 1) at (replicate_id, gid, start_month);
#' both sides are coerced to the same "YYYY-MM-DD" character `start_month` (the
#' kernel writes it as character; a Date round-trip normalises any parquet/qs
#' type drift). Events with no matching episode keep `episode_cases = NA`; their
#' count is returned as attribute "n_unmatched" and reported.
#'
#' @param event_table merged event table (replicate_id, gid, start_month, ...).
#' @param episode_totals tibble from compute_episode_case_totals().
#' @return event_table with `episode_cases` and `n_months_cases` appended.
attach_episode_case_totals <- function(event_table, episode_totals) {
  et <- event_table %>%
    dplyr::mutate(
      gid = as.character(gid),
      start_month = as.character(as.Date(start_month))
    ) %>%
    dplyr::left_join(
      episode_totals %>%
        dplyr::mutate(
          gid = as.character(gid),
          start_month = as.character(as.Date(start_month))
        ),
      by = c("replicate_id", "gid", "start_month")
    )
  n_unmatched <- sum(is.na(et$episode_cases))
  message(sprintf(
    "[attach_episode_case_totals] episode case totals joined onto %d event(s); %d unmatched (episode_cases = NA).",
    nrow(et),
    n_unmatched
  ))
  if (n_unmatched > 0L) {
    warning(sprintf(
      paste0(
        "[attach_episode_case_totals] %d of %d simulated event(s) have no ",
        "matching cases_month episode total (onset month not active?)."
      ),
      n_unmatched,
      nrow(et)
    ))
  }
  attr(et, "n_unmatched") <- n_unmatched
  et
}

# ---------------------------------------------------------------------------
# compute_crps_total_cases — annual total-case CRPS at three geographic levels
# ---------------------------------------------------------------------------

#' Compute CRPS for annual total suspected cases at Africa / region / country.
#'
#' Aggregates simulated per-episode case totals (`episode_cases`, from
#' compute_episode_case_totals(): the in-horizon sum of the mon track's
#' `cases_month` over each episode) from the event table and observed
#' `suspected_cases` from outbreak_metadata.rds by calendar year of ONSET at
#' three geographic levels. When a coverage mask is supplied, only (gid, month)
#' cells in the mask contribute; annual totals are pro-rated to a full year by
#' multiplying by 12 / n_covered_months. The unmasked sim totals are also
#' retained (column `total_sim_unmasked_mean`) so coverage-driven gaps can be
#' visualised.
#'
#' Estimand note (size track retired): the observed side is each outbreak's
#' whole-outbreak total, attributed to its onset year; the simulated side is
#' the sum of `cases_month` over the episode's simulated months. An episode censored at the
#' horizon end therefore contributes only its in-horizon months, so the last
#' horizon year(s) are biased low on the simulated side.
#'
#' @param event_table tibble with columns replicate_id, gid, start_month,
#'   episode_cases (merged/event_table with compute_episode_case_totals()
#'   joined on by run_validation()).
#' @param obs_metadata_path path to `~/secure/cholera_outbreak_data/outbreak_metadata.rds`.
#' @param frame tibble from build_location_frame() — must carry gid, country,
#'   AFRO_region.
#' @param obs_panel data frame from the survival panel — carries
#'   location_period_id and gid (used to bridge obs metadata to gids).
#' @param mask optional tibble(gid, month) from build_coverage_mask().
#' @param case_admin_level admin level the case (monthly) model was fitted at
#'   ("ADM1" | "ADM2" | "ADM1-ADM2"). Observed metadata is filtered to the same
#'   admin level so obs and sim compare like-for-like. Default "ADM1-ADM2" (the
#'   monthly track's default `monthly.admin_level`).
#' @return named list(africa, region, country) each a tibble with columns
#'   geography, year, n_sim, total_obs, crps_cases, total_sim_mean,
#'   total_sim_lo, total_sim_hi, n_covered_months.
compute_crps_total_cases <- function(
  event_table,
  obs_metadata_path,
  frame,
  obs_panel,
  mask = NULL,
  case_admin_level = "ADM1-ADM2",
  total_case_basis = c("exact_gid", "region_key"),
  members_path = NULL
) {
  total_case_basis <- match.arg(total_case_basis)
  if (!file.exists(obs_metadata_path)) {
    stop(
      "compute_crps_total_cases: obs_metadata_path not found: ",
      obs_metadata_path
    )
  }

  # Geo tag lookup from frame (only gid → geo tags; country/AFRO_region come
  # from outbreak_metadata directly to avoid duplicate-column joins)
  frame_geo <- .sim_frame_geo(frame)

  # ── Footprint-reconciled (rule-B) region key ──────────────────────────────
  # The simulator forecasts one admin level per region (ADM2 where nested
  # children exist, coarse ADM1 only where they don't — dropOverlappingADM1-
  # FromFrame), while the observed corpus reports the SAME geography at a mix of
  # ADM1/ADM2/ADM3. An exact-(gid) membership test silently drops observed cases
  # reported coarser than the sim forecasts (the cross-level blindspot). The
  # region key maps every gid to its ADM1 unit — the ISO3 prefix makes "NGA_20"
  # globally unique, so gid → region_key is unambiguous (no country column
  # needed) and a coarse observed report reconciles against the sim's finer
  # forecast children at the region level. Mirrors the deployed
  # notebooks/hindcast_skill.qmd §4d.
  # region_key_of(): utils/obs_sim_compare.R.
  sim_regions <- if (identical(total_case_basis, "region_key")) {
    unique(stats::na.omit(region_key_of(frame_geo$gid)))
  } else {
    NULL
  }

  # Bridge: location_period_id → gid (obs metadata lacks gid directly)
  lp_bridge <- obs_panel %>%
    dplyr::filter(!is.na(location_period_id) & !is.na(gid)) %>%
    dplyr::distinct(
      location_period_id = as.character(location_period_id),
      gid = as.character(gid)
    )
  # Under `episodes.key: gid_season` the survival panel carries ONE
  # representative location period per gid, so the bridge above would drop
  # every other member LP's observed cases. The canonical membership map
  # (utils/footprint.R::canonical_episodes()) lists every member LP; it also
  # carries the gid's GEOMETRY-derived admin level, which then replaces the
  # corpus label for the admin-level filter below -- the same level the fits
  # and the simulator frame use. NULL (every lp_season run) = unchanged.
  if (!is.null(members_path) && file.exists(members_path)) {
    lp_bridge <- readRDS(members_path) %>%
      dplyr::count(
        location_period_id = as.character(location_period_id),
        gid = as.character(gid),
        admin_level_canon = admin_level
      ) %>%
      dplyr::arrange(location_period_id, dplyr::desc(n)) %>%
      dplyr::distinct(location_period_id, .keep_all = TRUE) %>%
      dplyr::select(-n)
  }

  # ── Observed metadata ─────────────────────────────────────────────────────
  # outbreak_metadata already has country + AFRO_region. Join only gid.
  obs_meta <- readRDS(obs_metadata_path) %>%
    dplyr::mutate(
      location_period_id = as.character(location_period_id),
      suspected_cases = as.numeric(suspected_cases),
      obs_month = lubridate::floor_date(as.Date(TL), "month"),
      year = lubridate::year(as.Date(TL)),
      country = as.character(country),
      AFRO_region = as.character(AFRO_region),
      africa = "Africa"
    ) %>%
    dplyr::left_join(lp_bridge, by = "location_period_id")
  if ("admin_level_canon" %in% names(obs_meta)) {
    obs_meta <- obs_meta %>%
      dplyr::mutate(
        admin_level = dplyr::coalesce(
          as.character(admin_level_canon),
          as.character(admin_level)
        )
      ) %>%
      dplyr::select(-admin_level_canon)
  }

  # Restrict observed to the geographic footprint the simulator forecasts —
  # ensures partial-TU runs compare like-for-like (obs and sim cover the same
  # footprint).
  #   exact_gid  : keep obs only on a frame gid. Drops coarse observed reports
  #                whose gid is not itself a frame gid (the cross-level blindspot).
  #   region_key : keep obs whose ADM1 region the sim forecasts, so a coarse ADM1
  #                observed report is retained and reconciled (rule B) against the
  #                sim's finer forecast children.
  if (identical(total_case_basis, "region_key")) {
    obs_meta <- obs_meta %>%
      dplyr::mutate(region_key = region_key_of(gid)) %>%
      dplyr::filter(!is.na(region_key), region_key %in% sim_regions)
  } else {
    obs_meta <- dplyr::filter(obs_meta, gid %in% frame_geo$gid)
  }

  # Restrict observed metadata to the same admin level the case (monthly) model
  # was fitted at. Without this, ADM1 and ADM2 reports for the same area both
  # contribute, double-counting cases.
  # `ADM1-ADM2` needs an ADM1 column; derive it from gid as in 03_prepare_siz.R.
  if ("admin_level" %in% names(obs_meta)) {
    if (!"ADM1" %in% names(obs_meta)) {
      obs_meta <- dplyr::mutate(
        obs_meta,
        # flatten_gid() first: gid is the raw (dotted, "ISO3.ADM1[.ADM2]_N")
        # panel convention here, not an underscore-joined RE level string — see
        # the same idiom at 03_prepare_occ.R:561 / 03_prepare_siz.R:593.
        ADM1 = gid_to_adm1(gid, "compute_crps_total_cases")
      )
    }
    obs_meta <- apply_admin_level_filter(obs_meta, case_admin_level)
  } else {
    warning(
      "[compute_crps_total_cases] obs metadata lacks an `admin_level` ",
      "column — skipping admin-level filter."
    )
  }

  # ── Simulated event table ─────────────────────────────────────────────────
  sim_et <- event_table %>%
    dplyr::mutate(
      gid = as.character(gid),
      start_month = lubridate::floor_date(as.Date(start_month), "month"),
      year = lubridate::year(start_month),
      episode_cases = as.numeric(episode_cases),
      africa = "Africa"
    ) %>%
    # event_table already carries `country`; join only AFRO_region from frame_geo
    dplyr::left_join(frame_geo %>% dplyr::select(gid, AFRO_region), by = "gid")

  sim_et_unmasked <- sim_et # pre-mask copy for unmasked mean

  # ── Apply mask ────────────────────────────────────────────────────────────
  covered_frac <- NULL
  if (!is.null(mask)) {
    mask_norm <- mask %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month))
    if (identical(total_case_basis, "region_key")) {
      # Project coverage up to (region_key, month) and mask BOTH sides on it.
      #   OBS: a coarse observed report often sits on a gid absent from the
      #        (frame-gid-keyed) mask — its region is covered even though its gid
      #        is not, so key the coverage test on region.
      #   SIM: the ensemble forecasts a region at its FINER children (e.g. NGA_20
      #        reported coarse at ADM1, forecast at ADM2 NGA_20_1..NGA_20_k). The
      #        coverage mask carries only the coarse OBSERVED gid (NGA_20), so an
      #        exact (gid, month) semi_join drops every fine sim child and zeroes
      #        the ensemble in exactly the regions that DID report — the cross-
      #        level blindspot, leaking back in through the sim-side coverage join.
      #        Reconcile the sim side at (region_key, month) too so a fine child is
      #        covered iff its region is covered that month. Mirrors the deployed
      #        notebooks/hindcast_skill.qmd §4d coverage join (cov_region_month).
      mask_region <- mask_norm %>%
        dplyr::mutate(region_key = region_key_of(gid)) %>%
        dplyr::filter(!is.na(region_key)) %>%
        dplyr::distinct(region_key, month)
      obs_meta <- dplyr::semi_join(
        obs_meta,
        mask_region,
        by = c("region_key", "obs_month" = "month")
      )
      sim_et <- sim_et %>%
        dplyr::mutate(region_key = region_key_of(gid)) %>%
        dplyr::semi_join(
          mask_region,
          by = c("region_key", "start_month" = "month")
        ) %>%
        dplyr::select(-region_key)
    } else {
      obs_meta <- dplyr::semi_join(
        obs_meta,
        mask_norm,
        by = c("gid", "obs_month" = "month")
      )
      sim_et <- dplyr::semi_join(
        sim_et,
        mask_norm,
        by = c("gid", "start_month" = "month")
      )
    }
    covered_frac <- mask_norm %>%
      dplyr::mutate(year = lubridate::year(month)) %>%
      dplyr::group_by(gid, year) %>%
      dplyr::summarise(n_covered_months = dplyr::n(), .groups = "drop")
  }

  # ── Per-gid observed cases: sum of distinct outbreaks ────────────────────
  # outbreak_metadata.rds carries one row per outbreak `uid`
  # (= location-outbreak_number; see 01_preprocess_outbreaks.R:77-83) with its
  # suspected_cases already summed over the outbreak window. Multiple rows per
  # (gid, year) are therefore *distinct* outbreaks and must be summed — the
  # previous mean-rate-over-LPs aggregation deflated gids with several outbreaks
  # per year. Summing the absolute case counts also puts obs on the same footing
  # as the simulator's absolute `episode_cases`, removing the earlier pop vs
  # pop_gid population-basis mismatch (no population rescaling is involved).
  obs_gid_year <- obs_meta %>%
    dplyr::filter(!is.na(suspected_cases)) %>%
    dplyr::group_by(gid, year) %>%
    dplyr::summarise(
      obs_cases_gid = sum(suspected_cases, na.rm = TRUE),
      .groups = "drop"
    )

  # region_key mode: attribute observed cases to geography from the observed
  # metadata's OWN country/AFRO_region (a coarse observed gid is absent from
  # frame_geo, so the gid → geo join used by the exact_gid path would drop it).
  # Cases are still summed per distinct outbreak, exactly as obs_gid_year does.
  obs_geo_year <- NULL
  if (identical(total_case_basis, "region_key")) {
    obs_geo_year <- obs_meta %>%
      dplyr::filter(!is.na(suspected_cases)) %>%
      dplyr::group_by(africa, AFRO_region, country, year) %>%
      dplyr::summarise(
        obs_cases = sum(suspected_cases, na.rm = TRUE),
        .groups = "drop"
      )
  }

  crps_at_level <- function(geo_col) {
    sim_ann <- sim_et %>%
      dplyr::filter(!is.na(.data[[geo_col]])) %>%
      dplyr::group_by(replicate_id, geography = .data[[geo_col]], year) %>%
      dplyr::summarise(
        total_sim = sum(episode_cases, na.rm = TRUE),
        .groups = "drop"
      )

    sim_ann_u <- sim_et_unmasked %>%
      dplyr::filter(!is.na(.data[[geo_col]])) %>%
      dplyr::group_by(replicate_id, geography = .data[[geo_col]], year) %>%
      dplyr::summarise(
        total_u = sum(episode_cases, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(geography, year) %>%
      dplyr::summarise(
        total_sim_unmasked_mean = mean(total_u, na.rm = TRUE),
        .groups = "drop"
      )

    obs_ann <- if (identical(total_case_basis, "region_key")) {
      # geo columns already carried on obs_geo_year (obs metadata's own tags)
      obs_geo_year %>%
        dplyr::filter(!is.na(.data[[geo_col]])) %>%
        dplyr::group_by(geography = .data[[geo_col]], year) %>%
        dplyr::summarise(
          total_obs = sum(obs_cases, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      obs_gid_year %>%
        dplyr::left_join(
          frame_geo %>% dplyr::select(gid, dplyr::all_of(geo_col)),
          by = "gid"
        ) %>%
        dplyr::filter(!is.na(.data[[geo_col]])) %>%
        dplyr::group_by(geography = .data[[geo_col]], year) %>%
        dplyr::summarise(
          total_obs = sum(obs_cases_gid, na.rm = TRUE),
          .groups = "drop"
        )
    }

    .geo_crps_level(
      sim_ann,
      obs_ann,
      sim_ann_u,
      covered_frac,
      frame_geo,
      geo_col,
      crps_col = "crps_cases"
    )
  }

  list(
    africa = crps_at_level("africa"),
    region = crps_at_level("AFRO_region"),
    country = crps_at_level("country")
  )
}

# ---------------------------------------------------------------------------
# .km_estimate — minimal Kaplan-Meier for right-censored discrete-time data
# ---------------------------------------------------------------------------

#' Kaplan-Meier survival curve, plus median/p90 readoff.
#'
#' Deliberately minimal (no `survival` package dependency) -- this repo's
#' duration clock is already discrete (whole months), so the textbook
#' product-limit estimator collapses to a handful of vectorised steps. Ties
#' are handled the standard way (all events/censorings at time t are counted
#' together in that step's risk set).
#'
#' @param time numeric/integer vector; observed time (event time if
#'   `event == 1`, censoring time if `event == 0`).
#' @param event integer vector (0/1); 1 = observed cessation, 0 = right-censored.
#' @return list(curve = tibble(time, n_risk, n_event, surv), median, p90).
#'   `median`/`p90` are `NA_real_` if the curve never reaches that survival
#'   level (e.g. heavy censoring at the tail) -- reported as-is, not imputed.
.km_estimate <- function(time, event) {
  ok <- !is.na(time) & !is.na(event)
  time <- time[ok]
  event <- event[ok]
  if (length(time) == 0L) {
    return(list(
      curve = tibble::tibble(
        time = numeric(),
        n_risk = integer(),
        n_event = integer(),
        surv = numeric()
      ),
      median = NA_real_,
      p90 = NA_real_
    ))
  }

  event_times <- sort(unique(time[event == 1L]))
  if (length(event_times) == 0L) {
    return(list(
      curve = tibble::tibble(
        time = numeric(),
        n_risk = integer(),
        n_event = integer(),
        surv = numeric()
      ),
      median = NA_real_,
      p90 = NA_real_
    ))
  }

  n_risk <- vapply(event_times, function(t) sum(time >= t), integer(1))
  n_event <- vapply(
    event_times,
    function(t) sum(time == t & event == 1L),
    integer(1)
  )
  surv <- cumprod(1 - n_event / n_risk)

  curve <- tibble::tibble(
    time = event_times,
    n_risk = n_risk,
    n_event = n_event,
    surv = surv
  )

  readoff <- function(level) {
    below <- which(curve$surv <= level)
    if (!length(below)) {
      return(NA_real_)
    }
    curve$time[below[1]]
  }

  list(curve = curve, median = readoff(0.5), p90 = readoff(0.1))
}

# ---------------------------------------------------------------------------
# compute_duration_dist — simulated vs observed episode-duration distribution
# ---------------------------------------------------------------------------

#' Compare simulated vs observed episode-duration distributions (KM + ECDF).
#'
#' Simulated side comes from the event table produced by a duration-aware
#' kernel run (`duration_months`, `censored`; absent on legacy no-duration
#' runs, in which case this function returns NULL). Observed side, when
#' supplied, is the observed active-month panel (`analysis/03_prepare_dur.R`
#' output, one row per active month per episode) collapsed to one row per
#' episode via its `duration_months`/`censored` episode-constant columns.
#'
#' Uses Kaplan-Meier (via `.km_estimate()`) rather than a plain ECDF so
#' right-censored (horizon-truncated / still-open) episodes are not compared
#' against complete ones as if their observed duration were their true
#' duration.
#'
#' @param event_table tibble from the simulator's merged event_table (must
#'   carry duration_months, censored, country for a duration-aware run).
#' @param obs_dur_panel optional tibble read directly from the
#'   duration-panel RDS (`load_duration_data()` / `make_duration_dat_file()`
#'   output, or a raw `readRDS()` of that path). NULL skips the observed side.
#' @param enso_grouping_path path to the ENSO teleconnection registry (keyed
#'   by ISO3 `country`); used for the by-group breakdown. Missing file
#'   silently skips the by-group split (overall KM/ECDF still computed).
#' @return NULL for a legacy (no-duration) event_table. Otherwise
#'   list(sim = list(overall = km, by_group = named list of km),
#'        obs  = same structure or NULL).
compute_duration_dist <- function(
  event_table,
  obs_dur_panel = NULL,
  enso_grouping_path = "analysis/groupings/enso_grouping.csv"
) {
  if (!all(c("duration_months", "censored") %in% names(event_table))) {
    return(NULL) # legacy run: duration_months is always 1L, no censoring info
  }

  group_lut <- NULL
  if (file.exists(enso_grouping_path)) {
    grp_tbl <- utils::read.csv(enso_grouping_path, stringsAsFactors = FALSE)
    group_lut <- stats::setNames(grp_tbl$group, grp_tbl$iso3)
  }

  .side <- function(dur, cens, country) {
    overall <- .km_estimate(dur, 1L - cens)
    by_group <- NULL
    if (!is.null(group_lut)) {
      grp <- unname(group_lut[country])
      by_group <- list()
      for (g in sort(unique(stats::na.omit(grp)))) {
        sel <- which(grp == g)
        if (length(sel) < 5L) {
          next
        }
        by_group[[g]] <- .km_estimate(dur[sel], 1L - cens[sel])
      }
    }
    list(
      overall = overall,
      by_group = by_group,
      n_episodes = length(dur),
      n_censored = sum(cens == 1L, na.rm = TRUE)
    )
  }

  sim_side <- .side(
    event_table$duration_months,
    event_table$censored,
    event_table$country
  )

  obs_side <- NULL
  if (!is.null(obs_dur_panel)) {
    obs_ep <- obs_dur_panel %>%
      dplyr::distinct(episode_id, duration_months, censored, country)
    obs_side <- .side(
      obs_ep$duration_months,
      as.integer(obs_ep$censored),
      obs_ep$country
    )
  }

  list(sim = sim_side, obs = obs_side)
}

# ---------------------------------------------------------------------------
# compute_active_fraction — fraction-of-active-months by gid/country
# ---------------------------------------------------------------------------

#' Fraction of (gid, month) cells in the "active" state — the direct
#' duty-cycle check (contagion lags are fed the persistent
#' `active` state rather than the one-month onset pulse; this metric checks
#' the resulting duty cycle against the observed one directly, independent of
#' the contagion-conditional / CRPS metrics).
#'
#' @param sim_monthly tibble(replicate_id, gid, month, Y[, active]). Returns
#'   NULL if `active` is absent (legacy non-duration run).
#' @param frame tibble from build_location_frame() with gid, country.
#' @param obs_dur_panel optional tibble from the observed duration panel (one
#'   row per observed active month per episode; `gid`, `month`). NULL skips
#'   the observed side.
#' @param mask optional tibble(gid, month) from build_coverage_mask() /
#'   the panel-derived fallback in run_validation() -- supplies the
#'   observed-side denominator (total covered months per gid). Required for
#'   a meaningful observed fraction; without it the observed denominator
#'   falls back to the span of the duration panel itself, which under-covers
#'   inactive months and so overstates the observed active fraction.
#' @return NULL if `active` is absent from sim_monthly. Otherwise
#'   list(by_gid, by_country), each a tibble with sim (and, if available,
#'   obs) active fractions.
compute_active_fraction <- function(
  sim_monthly,
  frame,
  obs_dur_panel = NULL,
  mask = NULL
) {
  if (!"active" %in% names(sim_monthly)) {
    return(NULL)
  }

  ctry_map <- frame %>%
    dplyr::distinct(gid = as.character(gid), country = as.character(country))

  sim_in <- sim_monthly %>% dplyr::mutate(gid = as.character(gid))

  sim_by_gid <- sim_in %>%
    dplyr::group_by(gid) %>%
    dplyr::summarise(
      frac_active_sim = mean(active, na.rm = TRUE),
      .groups = "drop"
    )

  sim_by_country <- sim_in %>%
    dplyr::left_join(ctry_map, by = "gid") %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      frac_active_sim = mean(active, na.rm = TRUE),
      .groups = "drop"
    )

  obs_by_gid <- NULL
  obs_by_country <- NULL
  if (!is.null(obs_dur_panel)) {
    active_months_obs <- obs_dur_panel %>%
      dplyr::mutate(gid = as.character(gid), month = as.Date(month)) %>%
      dplyr::distinct(gid, month)

    denom <- if (!is.null(mask)) {
      mask %>%
        dplyr::mutate(gid = as.character(gid)) %>%
        dplyr::count(gid, name = "n_months")
    } else {
      active_months_obs %>%
        dplyr::group_by(gid) %>%
        dplyr::summarise(n_months = dplyr::n_distinct(month), .groups = "drop")
    }

    numer <- active_months_obs %>% dplyr::count(gid, name = "n_active")
    obs_by_gid <- dplyr::full_join(numer, denom, by = "gid") %>%
      dplyr::mutate(
        n_active = dplyr::coalesce(n_active, 0L),
        frac_active_obs = n_active / n_months
      )

    obs_by_country <- obs_by_gid %>%
      dplyr::left_join(ctry_map, by = "gid") %>%
      dplyr::filter(!is.na(country)) %>%
      dplyr::group_by(country) %>%
      dplyr::summarise(
        frac_active_obs = sum(n_active) / sum(n_months),
        .groups = "drop"
      )
  }

  by_gid <- if (is.null(obs_by_gid)) {
    sim_by_gid
  } else {
    dplyr::full_join(
      sim_by_gid,
      obs_by_gid %>% dplyr::select(gid, frac_active_obs),
      by = "gid"
    )
  }
  by_country <- if (is.null(obs_by_country)) {
    sim_by_country
  } else {
    dplyr::full_join(sim_by_country, obs_by_country, by = "country")
  }

  list(by_gid = by_gid, by_country = by_country)
}

# ---------------------------------------------------------------------------
# compute_sim_kendall_tau — simulated episode-cases/duration Kendall's tau
# ---------------------------------------------------------------------------

#' Simulated episode-cases <-> duration Kendall's tau, reusing `paired_kendall_tau()`
#' (`analysis/utils/dependence.R`) so the simulated and observed (Phase 5,
#' `05_fit_postprocess_dependence.R`) numbers are computed with the same
#' estimator.
#'
#' This is a RAW-value tau (pooled across all uncensored simulated episodes,
#' `stats::cor(episode_cases, duration_months, method = "kendall")` under the
#' hood), not a PIT-residual tau -- the simulated episodes are one realised
#' draw each, not posterior draws against a fixed panel, so there is no
#' PIT to compute. Kendall's tau is rank-invariant to monotonic transforms of
#' either margin on its OWN, so this is a valid tau, but it is not
#' confound-adjusted the way the fits' PIT-residual tau is: shared covariates
#' and REs can induce raw co-variation between simulated episode cases and duration
#' even when the posterior-draw-level (residual) association is zero. If the
#' fits show independence, this simulated tau is the falsifiable check that
#' the *kernel* (not just the fits) also shows no residual coupling -- read it
#' next to, not as identical to, the PIT-residual tau.
#'
#' Censored episodes are dropped before pooling (their `duration_months` is
#' the horizon-truncation point, not the true cessation time, which would
#' bias the pooled tau; their in-horizon `episode_cases` is truncated likewise).
#'
#' Estimand note (size track retired): the simulated size margin is
#' `episode_cases`, the sum of the mon track's `cases_month` over the episode
#' (compute_episode_case_totals()). Because `episode_cases` accumulates over the
#' episode's active months, it is mechanically coupled to `duration_months`, so
#' this tau is expected to be positive even when the per-month case draws are
#' independent of cessation; read it as a kernel diagnostic, not as directly
#' comparable to the Phase-5 (size-fit) PIT-residual tau.
#'
#' @param event_table tibble from the simulator's merged event_table (must
#'   carry episode_cases, duration_months, censored -- absent on no-duration
#'   or no-monthly runs, in which case this function returns NULL).
#' @param min_n passed to paired_kendall_tau().
#' @return NULL for a legacy (no-duration) event_table. Otherwise the return
#'   value of paired_kendall_tau(): list(overall, by_group = NULL, n_common).
compute_sim_kendall_tau <- function(event_table, min_n = 5L) {
  if (
    !all(
      c("episode_cases", "duration_months", "censored") %in% names(event_table)
    )
  ) {
    return(NULL)
  }
  et <- event_table %>% dplyr::filter(censored == 0L)
  if (nrow(et) < min_n) {
    message(sprintf(
      "compute_sim_kendall_tau: only %d uncensored simulated episode(s) (< min_n = %d) -- skipping.",
      nrow(et),
      min_n
    ))
    return(NULL)
  }
  ep_id <- paste(et$replicate_id, et$gid, et$start_month, sep = "|")
  pit_a <- list(
    u = matrix(as.numeric(et$episode_cases), nrow = 1),
    episode_id = ep_id
  )
  pit_b <- list(
    u = matrix(as.numeric(et$duration_months), nrow = 1),
    episode_id = ep_id
  )
  paired_kendall_tau(pit_a, pit_b, min_n = min_n)
}

# ---------------------------------------------------------------------------
# Hindcast figures (season-by-season forecast-skill test)
# ---------------------------------------------------------------------------
# These three functions are driven by the --hindcast CLI mode (mirrors the
# --hindcast mode of 07c_simulator_postprocess.R). For a single hindcast run
# (one launch_year x ENSO group) they overlay the simulated forecast band vs
# the observed record over the *scored* window, applying the SAME panel-derived
# coverage mask 07c uses (observed at-risk (gid, month) cells within the
# forecast horizon). See run_hindcast_figures() for the manifest loop.

#' Build the per-month simulated bands + observed series for one hindcast run.
#'
#' @param merged_dir directory holding merged/monthly_summary (the aggregated
#'   spatPomp forecast: replicate_id, gid, month, Y, cases_month, ...).
#' @param obs_panel_path RDS survival panel (gid, t_start, outbreak) — observed
#'   occurrence and the coverage mask are derived from it.
#' @param monthly_data_path RDS monthly-cases model frame (gid, month,
#'   cases_month — the per-month observed case count from 03_prepare_mon.R);
#'   may be NULL/absent (cases series then omitted).
#' @param launch,horizon_end YYYY-MM-01 forecast window bounds (inclusive).
#' @param out_fmt monthly_summary serialization ("qs"/"rds"/"parquet").
#' @param group_label,launch_year tags for figure titles.
#' @return list(occ_band, cases_band, cum_sim, obs_cum_cases, meta) or NULL if
#'   the merged monthly_summary is absent.
compute_hindcast_series <- function(
  merged_dir,
  obs_panel_path,
  monthly_data_path = NULL,
  launch,
  horizon_end,
  out_fmt = "qs",
  group_label = NA_character_,
  launch_year = NA_integer_
) {
  ms_path <- file.path(
    merged_dir,
    paste0("monthly_summary", sim_io_ext(out_fmt))
  )
  if (!file.exists(ms_path)) {
    return(NULL)
  }

  ms <- sim_io_read(ms_path, out_fmt) %>%
    dplyr::mutate(
      gid = as.character(gid),
      month = local({
        u <- unique(month)
        as.Date(u)[match(month, u)]
      })
    )
  h_start <- as.Date(launch)
  h_end <- as.Date(horizon_end)
  ms <- ms %>% dplyr::filter(month >= h_start, month <= h_end)
  sim_gids <- unique(ms$gid)
  n_rep <- dplyr::n_distinct(ms$replicate_id)

  # Observed occurrence + coverage mask (07c run_validation() fallback 2):
  # only observed at-risk (gid, month) cells inside the horizon are scored.
  panel <- readRDS(obs_panel_path)
  obs_occ <- load_observed_from_panel(
    panel,
    gids = sim_gids,
    horizon = list(start = launch, end = horizon_end)
  )
  mask <- obs_occ %>% dplyr::distinct(gid, month)

  qtiles <- function(x) {
    tibble::tibble(
      q025 = stats::quantile(x, 0.025, names = FALSE, na.rm = TRUE),
      q25 = stats::quantile(x, 0.25, names = FALSE, na.rm = TRUE),
      med = stats::median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q75 = stats::quantile(x, 0.75, names = FALSE, na.rm = TRUE),
      q975 = stats::quantile(x, 0.975, names = FALSE, na.rm = TRUE)
    )
  }

  # --- Occurrence band: active admin units per (replicate, month) over the mask
  sim_masked <- ms %>% dplyr::inner_join(mask, by = c("gid", "month"))
  occ_rep <- sim_masked %>%
    dplyr::group_by(replicate_id, month) %>%
    dplyr::summarise(n_active = sum(Y, na.rm = TRUE), .groups = "drop")
  occ_band <- occ_rep %>%
    dplyr::group_by(month) %>%
    dplyr::reframe(qtiles(n_active)) %>%
    dplyr::arrange(month)
  obs_occ_month <- obs_occ %>%
    dplyr::semi_join(mask, by = c("gid", "month")) %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(obs = sum(Y_obs, na.rm = TRUE), .groups = "drop")
  occ_band <- occ_band %>%
    dplyr::left_join(obs_occ_month, by = "month") %>%
    dplyr::mutate(obs = dplyr::coalesce(obs, 0))

  # --- Cases band: suspected cases per (replicate, month) over the mask.
  #     cases_month is NA on inactive cells -> treat as 0 contribution.
  cases_rep <- sim_masked %>%
    dplyr::group_by(replicate_id, month) %>%
    dplyr::summarise(cases = sum(cases_month, na.rm = TRUE), .groups = "drop")
  cases_band <- cases_rep %>%
    dplyr::group_by(month) %>%
    dplyr::reframe(qtiles(cases)) %>%
    dplyr::arrange(month)

  # Observed cases: per-month observed case count (`cases_month`) from the
  # monthly-cases model frame (03_prepare_mon.R), restricted to the same
  # coverage mask. This is a genuine monthly allocation, apples-to-apples with
  # the simulator's own per-month `cases_month`.
  obs_cases_month <- NULL
  obs_cum_cases <- NA_real_
  have_cases <- !is.null(monthly_data_path) &&
    nzchar(monthly_data_path) &&
    file.exists(monthly_data_path)
  if (have_cases) {
    mon_all <- readRDS(monthly_data_path)
    if (all(c("gid", "month", "cases_month") %in% names(mon_all))) {
      obs_cases_cell <- mon_all %>%
        dplyr::transmute(
          gid = as.character(gid),
          month = lubridate::floor_date(as.Date(month), "month"),
          cases_month = as.numeric(cases_month)
        ) %>%
        dplyr::filter(gid %in% sim_gids, month >= h_start, month <= h_end) %>%
        dplyr::semi_join(mask, by = c("gid", "month"))
      obs_cases_month <- obs_cases_cell %>%
        dplyr::group_by(month) %>%
        dplyr::summarise(obs = sum(cases_month, na.rm = TRUE), .groups = "drop")
      obs_cum_cases <- sum(obs_cases_cell$cases_month, na.rm = TRUE)
    } else {
      have_cases <- FALSE
    }
  }
  cases_band <- cases_band %>%
    dplyr::left_join(
      obs_cases_month %||%
        tibble::tibble(month = as.Date(character()), obs = numeric()),
      by = "month"
    ) %>%
    dplyr::mutate(obs = dplyr::coalesce(obs, 0))

  # --- Cumulative cases per replicate over the whole masked window.
  cum_sim <- sim_masked %>%
    dplyr::group_by(replicate_id) %>%
    dplyr::summarise(cum = sum(cases_month, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(cum)

  list(
    occ_band = occ_band,
    cases_band = cases_band,
    cum_sim = cum_sim,
    obs_cum_cases = obs_cum_cases,
    have_cases = have_cases && !is.null(obs_cases_month),
    meta = list(
      group = group_label,
      launch_year = launch_year,
      launch = launch,
      horizon_end = horizon_end,
      n_rep = n_rep,
      n_gid = length(sim_gids),
      n_cells = nrow(mask)
    )
  )
}
