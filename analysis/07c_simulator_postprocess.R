# 07c_simulator_postprocess.R — aggregate simulation partitions + validation
#   metrics
#
# Usage:
#   Rscript analysis/07c_simulator_postprocess.R \
#     -C analysis/configs/sim_smoketest/sim_smoke.yml --aggregate
#   Rscript analysis/07c_simulator_postprocess.R \
#     -C analysis/configs/sim_smoketest/sim_smoke.yml --validate \
#     --obs-panel generated_data/outbreaks_survival_data_base-1990_2019.rds
#
# Details:
#   --aggregate merges all successful partitions into merged/; --validate
#   computes acceptance metrics vs the observed panel (combinable).
#   --combine-groups (no -C; --manifest, --arm, --analog-set) combines a
#   forecast manifest's per-group runs into full Africa, truncated to
#   africa_min_horizon. --hindcast (no -C; --manifest, --obs-panel,
#   --hindcast-summary-dir) scores every hindcast run (launch-year skill).
#
# See: analysis/utils/sim_metrics.R, sim_io.R (metrics, partition IO)

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(qs2)
})

# 2. Validation report ---------------------------------------------------------
#' Write and print a plain-text validation report.
#'
#' @param metrics named list returned by run_validation().
#' @param path character; path to write validation_report.txt.
write_validation_summary <- function(metrics, path) {
  lines <- character(0)
  lines <- c(
    lines,
    "===========================================================",
    "  Cholera Simulator — Validation Report",
    sprintf("  Generated: %s", Sys.time()),
    "===========================================================",
    ""
  )

  # PIT
  if (!is.null(metrics$pit) && nrow(metrics$pit) > 0) {
    lines <- c(
      lines,
      "--- PIT (randomised probability integral transform) ---",
      sprintf("  n cells      : %d", nrow(metrics$pit)),
      sprintf(
        "  PIT mean     : %.4f  (expect ~0.5 for calibrated model)",
        mean(metrics$pit$pit, na.rm = TRUE)
      ),
      sprintf(
        "  PIT sd       : %.4f  (expect ~0.289 for U[0,1])",
        stats::sd(metrics$pit$pit, na.rm = TRUE)
      ),
      ""
    )
  }

  # CRPS (annual counts)
  if (!is.null(metrics$crps) && nrow(metrics$crps) > 0) {
    lines <- c(
      lines,
      "--- CRPS (annual outbreak counts) ---",
      sprintf("  n (gid, year): %d", nrow(metrics$crps)),
      sprintf(
        "  mean CRPS    : %.4f",
        mean(metrics$crps$crps_count, na.rm = TRUE)
      ),
      sprintf(
        "  median CRPS  : %.4f",
        stats::median(metrics$crps$crps_count, na.rm = TRUE)
      ),
      ""
    )
  }

  # CRPS (total cases)
  if (!is.null(metrics$crps_total_cases)) {
    for (lvl in c("africa", "region", "country")) {
      tbl <- metrics$crps_total_cases[[lvl]]
      if (!is.null(tbl) && nrow(tbl) > 0) {
        lines <- c(
          lines,
          sprintf(
            "--- CRPS (annual total cases, %s level; sim = in-horizon cases_month sum per episode, by onset year) ---",
            lvl
          ),
          sprintf("  n (geography, year): %d", nrow(tbl)),
          sprintf(
            "  mean CRPS          : %.1f",
            mean(tbl$crps_cases, na.rm = TRUE)
          ),
          sprintf(
            "  median CRPS        : %.1f",
            stats::median(tbl$crps_cases, na.rm = TRUE)
          ),
          ""
        )
      }
    }
  }

  # AR metrics
  if (!is.null(metrics$ar_metrics)) {
    ar <- metrics$ar_metrics
    lines <- c(
      lines,
      "--- AR sequence metrics ---",
      "  Run lengths (1-runs):",
      paste0("    ", utils::capture.output(print(ar$run_length_table))),
      "  Gap lengths (0-runs):",
      paste0("    ", utils::capture.output(print(ar$gap_length_table))),
      ""
    )
    if (!is.null(ar$transition) && nrow(ar$transition) > 0) {
      lines <- c(
        lines,
        "  Transition matrix (from → to | n | p):",
        paste0(
          "    ",
          utils::capture.output(print(as.data.frame(ar$transition)))
        ),
        ""
      )
    }
  }

  # Contagion conditionals
  if (!is.null(metrics$contagion)) {
    cg <- metrics$contagion
    fmt_cond <- function(nm, cond) {
      if (is.null(cond) || is.na(cond$lag_active)) {
        return(character(0))
      }
      sprintf(
        "  %-16s  P(Y|lag=1)=%.4f  P(Y|lag=0)=%.4f  P(Y)=%.4f",
        nm,
        cond$lag_active %||% NA_real_,
        cond$lag_inactive %||% NA_real_,
        cond$unconditional %||% NA_real_
      )
    }
    lines <- c(
      lines,
      "--- Contagion conditionals ---",
      fmt_cond("ADM1 lag", cg$adm1_cond),
      fmt_cond("Country lag", cg$country_cond),
      fmt_cond("TU lag", cg$tu_cond),
      sprintf("  Unconditional P(Y)=%.4f", cg$unconditional %||% NA_real_),
      ""
    )
  }

  # Duration distribution (sim vs obs KM)
  if (!is.null(metrics$duration_dist)) {
    dd <- metrics$duration_dist
    fmt_km <- function(nm, km) {
      if (is.null(km)) {
        return(character(0))
      }
      sprintf(
        "  %-10s  n=%-5d  n_censored=%-5d  median=%s d  p90=%s d",
        nm,
        km$n_episodes,
        km$n_censored,
        ifelse(
          is.na(km$overall$median),
          "NA",
          sprintf("%.0f", 30.44 * km$overall$median)
        ),
        ifelse(
          is.na(km$overall$p90),
          "NA",
          sprintf("%.0f", 30.44 * km$overall$p90)
        )
      )
    }
    lines <- c(
      lines,
      "--- Duration distribution (Kaplan-Meier; months converted to ~days) ---",
      fmt_km("Simulated", dd$sim),
      fmt_km("Observed", dd$obs),
      "  Target (QA panel): median 118-146 d, p90 216-335 d",
      ""
    )
  }

  # Fraction of active months (duty-cycle check)
  if (!is.null(metrics$active_fraction)) {
    af <- metrics$active_fraction
    lines <- c(
      lines,
      "--- Fraction of active months (duty-cycle check) ---",
      sprintf(
        "  mean frac_active (sim), by gid    : %.4f",
        mean(af$by_gid$frac_active_sim, na.rm = TRUE)
      ),
      sprintf(
        "  mean frac_active (sim), by country: %.4f",
        mean(af$by_country$frac_active_sim, na.rm = TRUE)
      )
    )
    if ("frac_active_obs" %in% names(af$by_gid)) {
      lines <- c(
        lines,
        sprintf(
          "  mean frac_active (obs), by gid    : %.4f",
          mean(af$by_gid$frac_active_obs, na.rm = TRUE)
        ),
        sprintf(
          "  mean frac_active (obs), by country: %.4f",
          mean(af$by_country$frac_active_obs, na.rm = TRUE)
        )
      )
    }
    lines <- c(lines, "")
  }

  # Simulated (and, if supplied, observed) episode-cases<->duration Kendall's tau
  if (!is.null(metrics$sim_kendall_tau)) {
    stk <- metrics$sim_kendall_tau$overall
    lines <- c(
      lines,
      "--- Episode cases (sum of cases_month) <-> duration Kendall's tau ---",
      sprintf(
        "  Simulated (raw, pooled, uncensored): tau = %.3f  (n_episodes = %d)",
        stk$mean,
        stk$n_episodes
      )
    )
    if (!is.null(metrics$obs_kendall_tau)) {
      otk <- metrics$obs_kendall_tau
      lines <- c(
        lines,
        sprintf(
          "  Observed (Phase 5 size fit, PIT-residual; not the same estimand): tau = %.3f  95%% CrI = [%.3f, %.3f]  (n_episodes = %d)",
          otk$mean,
          otk$q025,
          otk$q975,
          otk$n_episodes
        )
      )
    }
    lines <- c(lines, "")
  }

  lines <- c(
    lines,
    "==========================================================="
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  cat(paste(lines, collapse = "\n"), "\n")
}

# 3. Validation metrics --------------------------------------------------------
# run_validation — orchestrate all validation metrics (occurrence/AR families,
# the total-case metric from the mon track's cases_month, plus the
# duration-track metrics when the merged run carries an
# `active` column / event_table$duration_months)

#' Run all validation metrics and write outputs to `valid_dir/`.
#'
#' @param merged_dir directory with monthly_summary.rds (and event_table.rds
#'   when obs_metadata_path is supplied).
#' @param obs_panel_path path to the observed survival panel RDS.
#' @param frame tibble from build_location_frame().
#' @param valid_dir target directory for metrics.rds + validation_report.txt.
#' @param horizon optional list(start, end) for load_observed_from_panel().
#' @param max_reps integer; max replicates used for contagion reconstruction.
#' @param aggregated_ts_path path to `aggregated_ts_preoutbreak_extraction.rds`
#'   for the incidence-coverage mask; NULL skips the mask (metrics 2 and 3
#'   run unmasked).
#' @param aggregated_adm_ts_path path to
#'   `aggregated_admin_level_ts_preoutbreak_extraction.rds` (mask fallback).
#' @param obs_metadata_path path to `outbreak_metadata.rds` for the total-case
#'   CRPS (metric 2b); NULL skips it. 2b also needs a monthly-enabled run
#'   (monthly_summary carrying `active` and `cases_month`).
#' @param case_admin_level admin level of the monthly fit, forwarded to
#'   compute_crps_total_cases(). Default "ADM1-ADM2".
#' @param obs_dur_panel_path optional observed duration panel RDS
#'   (make_duration_dat_file()); enables the observed side of
#'   compute_duration_dist() / compute_active_fraction(). NULL skips it.
#' @param dependence_path optional `dependence_<date>.rds` from
#'   05_fit_postprocess_dependence.R: the observed PIT-residual Kendall's tau
#'   shown beside the simulated raw tau. NULL skips the comparison.
#' @param metric_set "full" (default) computes every metric. "hindcast"
#'   computes only what 09's hindcast scoring reads
#'   (hindcast_skill.R::hindcast_load_run()): pit, monthly_freq, crps and,
#'   when it runs, crps_total_cases; the others stay NULL. The kept metrics
#'   equal "full" (tests/test_sim_metrics_golden.R).
#' @return invisible(metrics).
run_validation <- function(
  merged_dir,
  obs_panel_path,
  frame,
  valid_dir,
  horizon = NULL,
  max_reps = 50L,
  aggregated_ts_path = NULL,
  aggregated_adm_ts_path = NULL,
  obs_metadata_path = NULL,
  out_fmt = "qs",
  case_admin_level = "ADM1-ADM2",
  total_case_basis = "exact_gid",
  obs_dur_panel_path = NULL,
  dependence_path = NULL,
  monthly_data_path = NULL,
  members_path = NULL,
  horizon_by_country = NULL,
  national_zero_path,
  metric_set = c("full", "hindcast")
) {
  metric_set <- match.arg(metric_set)
  full <- identical(metric_set, "full")
  dir.create(valid_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- sim_io_ext(out_fmt)

  message("[run_validation] loading merged monthly_summary", ext, " ...")
  # A compact (format-2) file is expanded with only the columns kept below.
  sim_monthly <- sim_io_read(
    file.path(merged_dir, paste0("monthly_summary", ext)),
    out_fmt,
    expand = FALSE
  )
  if (is_sim_monthly_compact(sim_monthly)) {
    sim_monthly <- sim_monthly_expand(
      sim_monthly,
      cols = c("Y", "active", "cases_month")
    )
  }
  sim_monthly <- sim_monthly %>%
    # drop mu_mon/detected — unused here. `active` is kept when present
    # (has_dur runs) so the duration-track metrics below can use it,
    # and `cases_month` when present (has_mon runs) for the per-episode case
    # totals behind metric 2b; any_of() is a no-op on runs without them.
    dplyr::select(dplyr::any_of(c(
      "replicate_id",
      "gid",
      "month",
      "Y",
      "active",
      "cases_month"
    ))) %>%
    # each distinct month string parsed once (as.Date() per row of a
    # character column would dominate this read's time and memory; same Dates)
    dplyr::mutate(
      gid = as.character(gid),
      month = local({
        u <- unique(month)
        as.Date(u)[match(month, u)]
      })
    )
  # A unit mixing fit groups was simulated to its max member horizon; score each
  # country only up to its own firm horizon (utils/sim_io.R).
  if (!is.null(horizon_by_country)) {
    sim_monthly <- censor_sim_by_country(sim_monthly, horizon_by_country)
  }
  has_dur <- "active" %in% names(sim_monthly)
  has_mon <- has_dur && "cases_month" %in% names(sim_monthly)
  gc()

  message("[run_validation] loading observed panel from: ", obs_panel_path)
  panel <- readRDS(obs_panel_path)
  obs_monthly <- load_observed_from_panel(
    panel,
    gids = unique(sim_monthly$gid),
    horizon = horizon
  )
  if (!is.null(horizon_by_country)) {
    obs_monthly <- censor_sim_by_country(
      obs_monthly %>% dplyr::mutate(gid = as.character(gid)),
      horizon_by_country
    )
  }

  obs_dur_panel <- NULL
  if (full && !is.null(obs_dur_panel_path) && nzchar(obs_dur_panel_path)) {
    if (file.exists(obs_dur_panel_path)) {
      message(
        "[run_validation] loading observed duration panel from: ",
        obs_dur_panel_path
      )
      obs_dur_panel <- readRDS(obs_dur_panel_path)
    } else {
      warning(
        "[run_validation] --obs-dur-panel not found: ",
        obs_dur_panel_path
      )
    }
  }

  # ── Build incidence-coverage mask ─────────────────────────────────────────
  mask <- NULL
  if (!is.null(aggregated_ts_path) && file.exists(aggregated_ts_path)) {
    message("[run_validation] building incidence-coverage mask ...")
    mask <- build_coverage_mask(
      # the panel's full LP -> gid member bridge when its sidecar exists, else
      # the panel's representative LP per gid (utils/coverage.R)
      obs_panel = coverage_lp_bridge(panel, obs_panel_path),
      aggregated_ts_path = aggregated_ts_path,
      aggregated_adm_ts_path = aggregated_adm_ts_path,
      national_zero_path = national_zero_path
    )
    if (!is.null(horizon_by_country)) {
      # the coverage mask spans the whole observed panel: censor the run's
      # countries, leave the others (never joined to this run) as they are.
      mask <- mask %>% dplyr::mutate(gid = as.character(gid))
      in_run <- gid_iso3(mask$gid) %in% names(horizon_by_country)
      mask <- dplyr::bind_rows(
        censor_sim_by_country(mask[in_run, , drop = FALSE], horizon_by_country),
        mask[!in_run, , drop = FALSE]
      )
    }
    message(sprintf(
      "[run_validation] mask: %d covered (gid, month) cells.",
      nrow(mask)
    ))
  } else {
    # Fallback panel mask: restrict simulation to the (gid, month) pairs
    # present in the observed panel. This aligns simulation and observed
    # denominators — without it the simulation spans the full horizon
    # (~174 months/gid) while the observed panel covers only ~85 months/gid,
    # inflating observed frequencies ~2× relative to simulated.
    message(
      "[run_validation] aggregated_ts not found — deriving panel mask from obs_monthly."
    )
    mask <- obs_monthly %>%
      dplyr::distinct(gid = as.character(gid), month = as.Date(month))
    message(sprintf(
      "[run_validation] panel mask: %d covered (gid, month) cells.",
      nrow(mask)
    ))
  }

  message("[run_validation] computing PIT ...")
  pit <- compute_pit_binom(sim_monthly, obs_monthly)

  message(
    "[run_validation] computing monthly frequency (mask=",
    !is.null(mask),
    ") ..."
  )
  monthly_freq <- compute_monthly_freq(
    sim_monthly,
    obs_monthly,
    frame,
    mask = mask
  )

  message(
    "[run_validation] computing annual CRPS (mask=",
    !is.null(mask),
    ") ..."
  )
  crps <- compute_crps_annual(sim_monthly, obs_monthly, mask = mask)

  # metric_set "hindcast": the event table is needed only by the total-case
  # CRPS, and only when that metric runs (obs metadata present, mon track).
  need_tc <- has_mon &&
    !is.null(obs_metadata_path) &&
    file.exists(obs_metadata_path)
  need_et <- full || need_tc
  if (!full) {
    message(
      "[run_validation] metric_set = 'hindcast': count by geography, AR, contagion, ",
      "active fraction, duration and Kendall's tau skipped",
      if (!need_tc) "; no total-case CRPS, event table not read" else "",
      "."
    )
  }

  crps_count_by_geo <- NULL
  ar_metrics <- NULL
  contagion <- NULL
  if (full) {
    message(
      "[run_validation] computing annual count by geography (mask=",
      !is.null(mask),
      ") ..."
    )
    crps_count_by_geo <- compute_annual_count_by_geo(
      sim_monthly,
      obs_monthly,
      frame,
      mask = mask
    )

    message("[run_validation] computing AR metrics ...")
    ar_metrics <- compute_ar_metrics(sim_monthly)

    message(
      "[run_validation] computing contagion conditionals (max_reps=",
      max_reps,
      ", active_col=",
      if (has_dur) "'active'" else "NULL",
      ") ..."
    )
    contagion <- compute_contagion_conditionals(
      sim_monthly,
      frame,
      max_reps = max_reps,
      active_col = if (has_dur) "active" else NULL
    )
  }

  # ── Fraction of active months (duty-cycle check) ─────────────────────────
  # Must run before sim_monthly is dropped below — needs the per-(gid,month)
  # `active` column, not the event_table.
  active_fraction <- NULL
  if (full && has_dur) {
    message("[run_validation] computing fraction of active months ...")
    active_fraction <- compute_active_fraction(
      sim_monthly,
      frame,
      obs_dur_panel = obs_dur_panel,
      mask = mask
    )
  }

  # ── Per-episode case totals (mon track) ──────────────────────────────────
  # Must also run before sim_monthly is dropped: the episode sums of
  # cases_month replace the retired size track's per-event `size_cases`.
  episode_totals <- NULL
  if (has_mon && need_et) {
    message(
      "[run_validation] computing per-episode case totals from cases_month ..."
    )
    episode_totals <- compute_episode_case_totals(sim_monthly)
  }

  rm(sim_monthly)
  gc() # free ~3 GB before loading event_table for metric 2b

  # ── Total-case CRPS (metric 2b) + duration metrics ───────────────────────
  crps_total_cases <- NULL
  duration_dist <- NULL
  sim_kendall_tau <- NULL
  obs_kendall_tau <- NULL
  et_path <- file.path(merged_dir, paste0("event_table", ext))
  event_table <- if (need_et && file.exists(et_path)) {
    sim_io_read(et_path, out_fmt)
  } else {
    NULL
  }

  if (!is.null(event_table) && !is.null(episode_totals)) {
    event_table <- attach_episode_case_totals(event_table, episode_totals)
  }
  rm(episode_totals)

  if (!has_mon && !is.null(obs_metadata_path)) {
    message(
      "[run_validation] monthly_summary carries no active/cases_month (no mon track) ",
      "— skipping total-case CRPS (metric 2b)."
    )
  } else if (
    !is.null(event_table) &&
      !is.null(obs_metadata_path) &&
      file.exists(obs_metadata_path)
  ) {
    message("[run_validation] computing total-case CRPS (metric 2b) ...")
    # region_key totals use an LP-keyed observed mask: a cell is covered iff
    # month <= its own series_end (dur_data.R::derive_series_end()). The default
    # ADM2/ADM3 mask drops coarse-reporting countries, and the include_adm1 mask
    # counts whole at-risk spans and inflates "over" outliers. Kept separate so
    # `mask` (occurrence-basis metrics) is unchanged.
    tc_mask <- mask
    if (identical(total_case_basis, "region_key")) {
      if (!is.null(monthly_data_path) && file.exists(monthly_data_path)) {
        message(
          "[run_validation] building LP-keyed total-case coverage mask (region_key) from: ",
          monthly_data_path
        )
        tc_mask <- readRDS(monthly_data_path) %>%
          dplyr::transmute(
            gid = as.character(gid),
            month = lubridate::floor_date(as.Date(month), "month"),
            series_end = lubridate::floor_date(as.Date(series_end), "month")
          ) %>%
          dplyr::filter(!is.na(series_end), month <= series_end) %>%
          dplyr::distinct(gid, month)
        message(sprintf(
          "[run_validation] region_key total-case mask: %d covered (gid, month) cells.",
          nrow(tc_mask)
        ))
      } else {
        warning(
          "[run_validation] total_case_basis='region_key' but monthly_data_path ",
          "is unavailable (",
          monthly_data_path %||% "NULL",
          ") — falling back to the ",
          "occurrence-basis mask; coarse-reporting countries (e.g. NGA) may be dropped."
        )
      }
    }
    crps_total_cases <- compute_crps_total_cases(
      event_table = event_table,
      obs_metadata_path = obs_metadata_path,
      frame = frame,
      obs_panel = panel,
      mask = tc_mask,
      case_admin_level = case_admin_level,
      total_case_basis = total_case_basis,
      members_path = members_path
    )
  } else if (!is.null(obs_metadata_path)) {
    if (!file.exists(obs_metadata_path)) {
      warning(
        "[run_validation] obs_metadata_path not found — skipping total-case CRPS."
      )
    }
    if (need_et && is.null(event_table)) {
      warning(sprintf(
        "[run_validation] event_table%s not found in merged_dir — skipping total-case CRPS.",
        ext
      ))
    }
  }

  if (
    full && !is.null(event_table) && "duration_months" %in% names(event_table)
  ) {
    message(
      "[run_validation] computing duration distribution (Kaplan-Meier) ..."
    )
    duration_dist <- compute_duration_dist(
      event_table,
      obs_dur_panel = obs_dur_panel
    )

    message(
      "[run_validation] computing simulated episode-cases<->duration Kendall's tau ..."
    )
    sim_kendall_tau <- compute_sim_kendall_tau(event_table)

    if (!is.null(dependence_path) && nzchar(dependence_path)) {
      if (file.exists(dependence_path)) {
        message(
          "[run_validation] loading observed (Phase 5) dependence result from: ",
          dependence_path
        )
        dep <- tryCatch(readRDS(dependence_path), error = function(e) NULL)
        if (!is.null(dep) && !is.null(dep$tau) && !is.null(dep$tau$overall)) {
          obs_kendall_tau <- dep$tau$overall
        } else {
          warning(
            "[run_validation] dependence_path did not contain the expected tau$overall structure — skipping observed tau."
          )
        }
      } else {
        warning(
          "[run_validation] --dependence-rds not found: ",
          dependence_path
        )
      }
    }
  }

  metrics <- list(
    pit = pit,
    monthly_freq = monthly_freq,
    crps = crps,
    crps_count_by_geo = crps_count_by_geo,
    ar_metrics = ar_metrics,
    contagion = contagion,
    crps_total_cases = crps_total_cases,
    duration_dist = duration_dist,
    active_fraction = active_fraction,
    sim_kendall_tau = sim_kendall_tau,
    obs_kendall_tau = obs_kendall_tau
  )

  # `metrics` is a non-tabular named list, so parquet cannot serialize it;
  # honor rds/qs and fall back to qs for parquet configs.
  metrics_fmt <- sim_metrics_fmt(out_fmt)
  metrics_path <- file.path(
    valid_dir,
    paste0("metrics", sim_io_ext(metrics_fmt))
  )
  sim_io_save(metrics, metrics_path, metrics_fmt)
  message("[run_validation] metrics saved to: ", metrics_path)

  write_validation_summary(
    metrics,
    file.path(valid_dir, "validation_report.txt")
  )

  invisible(metrics)
}

# 4. Validate one run ----------------------------------------------------------
# .validate_run_from_opt — build the location frame from a fully-resolved `opt`
# and run all validation metrics for one run, returning the metrics list.
# Factored out of the --validate CLI path so the --hindcast orchestrator can
# score each per-(group, year) run with the identical frame/mask/metric
# machinery instead of duplicating it. `horizon_override` (list(start?, end?))
# restricts the observed window to the forecast window [launch, horizon_end].
.validate_run_from_opt <- function(
  opt,
  merged_dir,
  valid_dir,
  obs_panel_path,
  horizon_override = NULL,
  horizon_by_country = NULL,
  metric_set = "full"
) {
  out_fmt_v <- opt$simulator$output$format %||% "qs"
  if (
    !file.exists(file.path(
      merged_dir,
      paste0("monthly_summary", sim_io_ext(out_fmt_v))
    ))
  ) {
    stop(
      "merged/monthly_summary",
      sim_io_ext(out_fmt_v),
      " not found in ",
      merged_dir,
      " — run --aggregate first."
    )
  }

  panel <- readRDS(opt$simulator$climate_panel)
  lf_opts <- opt$simulator$location_frame %||% list()
  # admin_level source for the frame — no fits are loaded on this path,
  # so this read exists solely to supply `admin_level_data`; freed immediately
  # after the call.
  .surv_overlap <- tryCatch(
    readRDS(make_surv_dat_file(list(
      baseline = opt$baseline,
      episodes = opt$episodes
    ))),
    error = function(e) NULL
  )
  # Episode definition: the observed panel scored against, the monthly frame
  # and the size frame must share one episode digest whenever any of them -- or
  # this run's episodes.key -- is gid_season.
  gid_season_v <- identical(episode_key_of(opt), "gid_season")
  sim_episode_check(
    preps = list(
      obs_panel = obs_panel_path,
      mon_data = make_monthly_dat_file(opt),
      size_data = lf_opts$size_data_path %||% make_size_dat_file(opt)
    ),
    what = "07c --validate inputs",
    episode_key = episode_key_of(opt)
  )
  .fpop <- sim_frame_pop_args(opt$simulator, opt)
  .fpop$pop_window <- sim_pop_window(opt$simulator)
  frame <- build_location_frame(
    panel,
    hmm_clusters_path = lf_opts$hmm_clusters_path,
    outbreak_metadata_path = lf_opts$outbreak_metadata_path,
    size_data_path = lf_opts$size_data_path,
    base_rate_table_paths = resolve_base_rate_paths(
      list(occ = opt, dur = opt, mon = opt),
      lf_opts
    ),
    admin_level_data = .surv_overlap,
    # Population arguments in one place (simulator_core.R::sim_frame_pop_args():
    # location_frame.gid_pop_path / wpp_path / pop_year).
    gid_pop_path = .fpop$gid_pop_path,
    wpp_path = .fpop$wpp_path,
    pop_year = .fpop$pop_year,
    pop_window = .fpop$pop_window
  )
  rm(.surv_overlap)

  horizon_v <- opt$simulator$horizon
  if (!is.null(horizon_override)) {
    if (!is.null(horizon_override$start) && nzchar(horizon_override$start)) {
      horizon_v$start <- horizon_override$start
    }
    if (!is.null(horizon_override$end) && nzchar(horizon_override$end)) {
      horizon_v$end <- horizon_override$end
    }
  }

  run_validation(
    merged_dir = merged_dir,
    obs_panel_path = obs_panel_path,
    frame = frame,
    valid_dir = valid_dir,
    horizon = horizon_v,
    aggregated_ts_path = make_agg_ts_file(),
    aggregated_adm_ts_path = make_agg_adm_ts_file(),
    obs_metadata_path = lf_opts$outbreak_metadata_path %||% NULL,
    out_fmt = out_fmt_v,
    case_admin_level = opt$monthly$admin_level %||% "ADM1-ADM2",
    total_case_basis = opt$simulator[["validation"]][["total_case_basis"]] %||%
      "exact_gid",
    obs_dur_panel_path = opt$obs_dur_panel %||% NULL,
    dependence_path = opt$dependence_path %||% NULL,
    monthly_data_path = make_monthly_dat_file(opt),
    members_path = if (gid_season_v) make_footprint_members_file(opt) else NULL,
    horizon_by_country = horizon_by_country,
    national_zero_path = national_zero_months_file_of(opt), # threshold + rule of the run's panel
    metric_set = metric_set
  )
}

# 5. Hindcast skill ------------------------------------------------------------
# run_hindcast_skill — season-by-season forecast-skill orchestrator.
#' Score every hindcast run in `manifest_path` (hindcast_horizons.yml, produced
#' by analysis/utils/make_hindcast_configs.R) against the observed panel and
#' bind the per-run skill metrics into a combined, launch-year-tagged table.
#'
#' Reuses the standard aggregation + validation machinery verbatim: each
#' per-(group, year) config is re-resolved through make_options_from_config()
#' (so its simulator block/defaults are identical to a direct -C run), its
#' partitions aggregated if needed, then scored via .validate_run_from_opt()
#' with the observed window clipped to the forecast window [launch, horizon_end].
#' run_validation() already restricts scoring to observed (gid, month) cells via
#' its coverage mask, so hindcast skill is measured only where surveillance
#' exists inside each forecast window.
#'
#' @return invisible(list(headline, crps_by_gid_year, monthly_freq)).
run_hindcast_skill <- function(
  manifest_path,
  obs_panel_path,
  summary_dir,
  option_list,
  do_aggregate = TRUE,
  workers = 1L
) {
  man <- read_hindcast_manifest(manifest_path, obs_panel_path) # utils/sim_io.R

  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  # Build a missing merged monthly_summary from the run's partitions.
  aggregate_if_missing <- if (isTRUE(do_aggregate)) {
    function(opt_run, md, out_fmt_run) {
      parts_dir <- make_run_file(opt_run, kind = "spatpomp_parts_dir")
      tus_scope <- opt_run$simulator$scope$transmission_units %||% "all"
      tus_arg <- if (identical(tus_scope, "all")) {
        NULL
      } else {
        as.character(tus_scope)
      }
      message("[hindcast] aggregating partitions -> ", md)
      aggregate_partitions(
        parts_dir,
        md,
        tus = tus_arg,
        out_fmt = out_fmt_run,
        return_long = FALSE
      )
    }
  }

  score_run <- function(run) {
    g <- run$group
    yr <- run$year
    ye <- run$entry
    opt_run <- run$opt
    md <- run$merged_dir
    valid_dir <- make_run_file(opt_run, kind = "spatpomp_valid_dir")
    # v2 (registry) manifests: the unit ran to the max member window; each
    # country is scored only inside its own fit group's firm window (always
    # censored, never kept: observed ENSO past the firm horizon would inflate
    # skill).
    metrics <- .validate_run_from_opt(
      opt_run,
      md,
      valid_dir,
      obs_panel_path,
      horizon_override = list(start = ye$launch, end = ye$horizon_end),
      horizon_by_country = run$firm_horizon_by_country,
      metric_set = "hindcast"
    ) # only what 09 reads

    ly <- run$launch_year # the season year (the key may be "<year>-<MM>")

    crps_row <- NULL
    freq_row <- NULL
    if (!is.null(metrics$crps) && nrow(metrics$crps)) {
      crps_row <- metrics$crps %>%
        dplyr::mutate(
          group = g,
          launch_year = ly,
          launch = ye$launch,
          horizon_end = ye$horizon_end
        )
    }
    if (!is.null(metrics$monthly_freq) && nrow(metrics$monthly_freq)) {
      freq_row <- metrics$monthly_freq %>%
        dplyr::mutate(group = g, launch_year = ly)
    }

    pit_ks <- NA_real_
    if (!is.null(metrics$pit) && nrow(metrics$pit) > 1L) {
      pit_ks <- suppressWarnings(
        as.numeric(stats::ks.test(metrics$pit$pit, "punif")$statistic)
      )
    }
    mean_crps <- if (!is.null(metrics$crps) && nrow(metrics$crps)) {
      mean(metrics$crps$crps_count, na.rm = TRUE)
    } else {
      NA_real_
    }
    ctc <- metrics$crps_total_cases$country
    mean_crps_cases <- if (
      !is.null(ctc) && "crps_cases" %in% names(ctc) && nrow(ctc)
    ) {
      mean(ctc$crps_cases, na.rm = TRUE)
    } else {
      NA_real_
    }

    list(
      crps = crps_row,
      freq = freq_row,
      headline = tibble::tibble(
        group = g,
        launch_year = ly,
        launch = ye$launch %||% NA_character_,
        horizon_end = ye$horizon_end %||% NA_character_,
        n_gid_year = if (!is.null(metrics$crps)) nrow(metrics$crps) else 0L,
        mean_crps_count = mean_crps,
        mean_crps_total_cases = mean_crps_cases,
        pit_ks = pit_ks,
        n_pit = if (!is.null(metrics$pit)) nrow(metrics$pit) else 0L,
        valid_dir = valid_dir
      )
    )
  }

  it <- for_each_hindcast_run(
    man,
    option_list,
    "[hindcast]",
    score_run,
    announce = "before",
    ensure_summary = aggregate_if_missing,
    workers = workers
  )
  n_run <- it$n_ok
  n_skip <- it$n_skip
  pick <- function(k) Filter(Negate(is.null), lapply(it$results, `[[`, k))
  crps_rows <- pick("crps")
  freq_rows <- pick("freq")
  headline <- pick("headline")

  if (!n_run) {
    stop(
      "[hindcast] no runs scored — check the manifest merged dirs / aggregation."
    )
  }

  headline_tbl <- dplyr::bind_rows(headline) %>%
    dplyr::arrange(group, launch_year)
  crps_tbl <- if (length(crps_rows)) {
    dplyr::bind_rows(crps_rows)
  } else {
    tibble::tibble()
  }
  freq_tbl <- if (length(freq_rows)) {
    dplyr::bind_rows(freq_rows)
  } else {
    tibble::tibble()
  }

  metrics_fmt <- sim_metrics_fmt(man$out_fmt)
  ext <- sim_io_ext(metrics_fmt)
  sim_io_save(
    headline_tbl,
    file.path(summary_dir, paste0("hindcast_skill_headline", ext)),
    metrics_fmt
  )
  if (nrow(crps_tbl)) {
    sim_io_save(
      crps_tbl,
      file.path(summary_dir, paste0("hindcast_crps_by_gid_year", ext)),
      metrics_fmt
    )
  }
  if (nrow(freq_tbl)) {
    sim_io_save(
      freq_tbl,
      file.path(summary_dir, paste0("hindcast_monthly_freq", ext)),
      metrics_fmt
    )
  }
  # A human-readable headline copy regardless of the binary format.
  utils::write.csv(
    headline_tbl,
    file.path(summary_dir, "hindcast_skill_headline.csv"),
    row.names = FALSE
  )

  message(sprintf(
    "\n[hindcast] ------------------------------------------------"
  ))
  message(sprintf("  runs scored : %d   skipped: %d", n_run, n_skip))
  message(sprintf(
    "  groups      : %s",
    paste(unique(headline_tbl$group), collapse = ", ")
  ))
  message(sprintf(
    "  launch years: %s",
    paste(range(headline_tbl$launch_year, na.rm = TRUE), collapse = "..")
  ))
  message(sprintf("  summary  -> %s", summary_dir))
  invisible(list(
    headline = headline_tbl,
    crps_by_gid_year = crps_tbl,
    monthly_freq = freq_tbl
  ))
}

# 6. CLI -----------------------------------------------------------------------

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({
    library(optparse)
    library(yaml)
  })
  source("analysis/utils/config.R")
  source("analysis/utils/formula.R") # model-term assembly
  source("analysis/utils/base.R")
  source("analysis/utils/gid.R")
  source("analysis/utils/paths.R")
  source("analysis/utils/clim_compile.R")
  source("analysis/utils/clim_db.R")
  source("analysis/utils/occ_data.R")
  source("analysis/utils/siz_data.R")
  source("analysis/utils/provenance.R")
  source("analysis/utils/sim_io.R")
  source("analysis/utils/coverage.R") # build_coverage_mask()
  source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
  source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
  source("analysis/utils/simulator_core.R")
  source("analysis/utils/dependence.R") # paired_kendall_tau() -- reused by compute_sim_kendall_tau()
  source("analysis/utils/obs_sim_compare.R") # crps_ens()/region_key_of() -- used by sim_metrics.R
  source("analysis/utils/sim_metrics.R") # the compute_* metrics run_validation() orchestrates
  source("analysis/utils/anchor.R") # resolve_base_rate_paths() -- base_rate_table_paths for build_location_frame() on the --validate path

  option_list <- list(
    make_option(
      "--aggregate",
      default = FALSE,
      action = "store_true",
      help = "Merge partitions into merged/."
    ),
    make_option(
      "--validate",
      default = FALSE,
      action = "store_true",
      help = "Compute validation metrics vs observed panel."
    ),
    make_option(
      c("--obs-panel"),
      default = NULL,
      type = "character",
      dest = "obs_panel",
      help = "Path to observed survival panel RDS (required for --validate)."
    ),
    make_option(
      "--obs-dur-panel",
      default = NULL,
      type = "character",
      dest = "obs_dur_panel",
      help = paste(
        "Path to the observed duration panel RDS (analysis/03_prepare_dur.R",
        "output / make_duration_dat_file()) -- optional; enables the",
        "simulated-vs-observed duration KM comparison and duty-cycle",
        "(fraction-of-active-months) check (Phase 7)."
      )
    ),
    make_option(
      "--dependence-rds",
      default = NULL,
      type = "character",
      dest = "dependence_path",
      help = paste(
        "Path to a dependence_<date>.rds artefact from",
        "05_fit_postprocess_dependence.R -- optional; supplies the",
        "observed (PIT-residual) size<->duration Kendall's tau to report",
        "alongside the simulated (raw, pooled) tau (Phase 7)."
      )
    ),
    make_option(
      "--combine-groups",
      default = FALSE,
      action = "store_true",
      dest = "combine_groups",
      help = "Combine per-group forecast runs into a full-Africa dataset truncated to africa_min_horizon (needs --manifest; no -C config)."
    ),
    make_option(
      "--workers",
      default = 1L,
      type = "integer",
      help = "--hindcast: forked workers scoring runs in parallel (one run each; serial when 1) [default: %default]"
    ),
    make_option(
      "--hindcast",
      default = FALSE,
      action = "store_true",
      help = paste(
        "Season-by-season forecast-skill test: aggregate + validate every",
        "per-(group, year) run in --manifest (hindcast_horizons.yml) and write",
        "a combined launch-year-tagged skill table. Needs --manifest + --obs-panel;",
        "no -C config."
      )
    ),
    make_option(
      "--hindcast-summary-dir",
      default = "generated_data/sim_spatpomp/hindcast_summary",
      type = "character",
      dest = "hindcast_summary_dir",
      help = "Output dir for the combined --hindcast skill tables."
    ),
    make_option(
      "--manifest",
      default = NULL,
      type = "character",
      help = "forecast_horizons.yml / hindcast_horizons.yml manifest (required for --combine-groups / --hindcast)."
    ),
    make_option(
      "--arm",
      default = "factual",
      type = "character",
      help = "Which arm to combine: factual | counterfactual."
    ),
    make_option(
      "--analog-set",
      default = NULL,
      type = "character",
      dest = "analog_set",
      help = paste(
        "For --arm counterfactual: which Lenssen analog set to combine",
        "(all | nina | neutral | elnino), or 'each' to combine every set",
        "present in the manifest. Default: the manifest's",
        "primary_analog_set."
      )
    ),
    make_option(
      "--out-fmt",
      default = "qs",
      type = "character",
      dest = "out_fmt",
      help = "I/O format for --combine-groups (rds | parquet | qs)."
    ),
    make_option(
      "--merged-dir",
      default = NULL,
      type = "character",
      dest = "merged_dir",
      help = "Override the merged/ directory used by --aggregate / --validate (e.g. a combined-groups dir)."
    ),
    make_option(
      "--horizon-end",
      default = NULL,
      type = "character",
      dest = "horizon_end",
      help = "Override simulator.horizon.end for --validate (e.g. africa_min_horizon)."
    ),
    make_option(
      "--engine",
      default = "v1",
      type = "character",
      help = "Output layout to resolve --aggregate/--validate paths against: 'v1' (07b_simulator_v1.R, generated_data/sim/<run_tag>/[run_id/], default) or 'spatpomp' (07b_simulator_spatpomp.R, generated_data/sim_spatpomp/<run_tag>/ — flat, NOT run_id-scoped; see make_spatpomp_run_file() in utils.R)."
    )
  )
  opt <- make_options_from_config(
    option_list = option_list,
    enforce_options = character(0),
    model = NULL
  )

  # --- Full-Africa combine mode (no simulator config required) -------------
  if (isTRUE(opt$combine_groups)) {
    if (
      is.null(opt$manifest) ||
        !nzchar(opt$manifest) ||
        !file.exists(opt$manifest)
    ) {
      stop("--combine-groups requires --manifest <forecast_horizons.yml>.")
    }
    arm <- opt$arm %||% "factual"
    if (!arm %in% c("factual", "counterfactual")) {
      stop("--arm must be 'factual' or 'counterfactual'.")
    }
    man <- yaml::read_yaml(opt$manifest)
    if (is.null(man$groups) || !length(man$groups)) {
      stop("Manifest has no `groups`: ", opt$manifest)
    }
    # The counterfactual arm is fanned out over Lenssen analog sets, so there is one
    # merged tree per (group, set). Combining only the manifest's flat
    # `counterfactual_merged_dir` alias would silently combine the primary set alone
    # and leave the other sets un-combined.
    if (identical(arm, "factual")) {
      sets_to_do <- list(NULL)
    } else {
      avail <- names(man$analog_sets_present %||% man$analog_sets %||% list())
      if (!length(avail)) {
        avail <- as.character(unlist(man$analog_sets %||% list()))
      }
      req <- opt$analog_set %||% (man$primary_analog_set %||% "all")
      if (identical(req, "each")) {
        if (!length(avail)) {
          stop(
            "--analog-set each requires `analog_sets` in the manifest; ",
            "regenerate it with analysis/utils/make_forecast_configs.R."
          )
        }
        sets_to_do <- as.list(avail)
      } else {
        if (length(avail) && !req %in% avail) {
          stop(
            "--analog-set '",
            req,
            "' is not in the manifest's analog_sets (",
            paste(avail, collapse = ", "),
            ")."
          )
        }
        sets_to_do <- list(req)
      }
    }

    # Runs are per unit: a fit group (v1) or a partition-registry unit (v2,
    # utils/sim_io.R::manifest_runs()). v2 units may mix fit groups, so each
    # country is also censored at its own firm horizon.
    hz_by_country <- if (length(man$units)) {
      manifest_country_horizons(man, "firm")
    } else {
      NULL
    }
    for (aset in sets_to_do) {
      if (is.null(aset)) {
        gmd <- lapply(manifest_runs(man, "factual"), `[[`, "merged_dir")
        combined_dir <- man$combined[["factual"]] %||%
          stop("Manifest lacks combined$factual output dir.")
        lab <- "factual"
      } else {
        gmd <- lapply(
          manifest_runs(man, "counterfactual", aset),
          `[[`,
          "merged_dir"
        )
        cc <- man$combined$counterfactual_sets[[aset]]$merged_dir %||%
          (if (identical(aset, man$primary_analog_set)) {
            man$combined$counterfactual
          } else {
            NULL
          })
        combined_dir <- cc %||%
          stop(
            "Manifest lacks a combined counterfactual dir for analog_set '",
            aset,
            "'."
          )
        lab <- paste0("counterfactual set-", aset)
      }
      gmd <- gmd[!vapply(gmd, is.null, logical(1))]
      if (!length(gmd)) {
        stop("No per-unit merged dirs found in the manifest for ", lab, ".")
      }
      combine_group_runs(
        gmd,
        combined_dir,
        man$africa_min_horizon,
        out_fmt = opt$out_fmt %||% "qs",
        horizon_by_country = hz_by_country
      )
      message("[07b] combine-groups (", lab, ") done -> ", combined_dir)
    }
    quit(save = "no", status = 0)
  }

  # --- Season-by-season forecast-skill test (no -C config; manifest-driven) --
  if (isTRUE(opt$hindcast)) {
    run_hindcast_skill(
      manifest_path = opt$manifest,
      obs_panel_path = opt$obs_panel,
      summary_dir = opt$hindcast_summary_dir %||%
        "generated_data/sim_spatpomp/hindcast_summary",
      option_list = option_list,
      do_aggregate = TRUE,
      workers = opt$workers %||% 1L
    )
    message("[07c_simulator_postprocess] hindcast skill done.")
    quit(save = "no", status = 0)
  }

  if (is.null(opt$simulator)) {
    stop("Config must contain a `simulator:` block.")
  }
  if (!isTRUE(opt$aggregate) && !isTRUE(opt$validate)) {
    stop("Specify at least one of --aggregate or --validate.")
  }

  # v1 / spatpomp run layout (shared with 07d_simulator_figures.R).
  dirs <- make_sim_engine_dirs(
    opt,
    opt$engine %||% "v1",
    merged_dir = opt$merged_dir
  ) # utils/paths.R
  run_dir <- dirs$run_dir
  parts_dir <- dirs$parts_dir
  merged_dir <- dirs$merged_dir
  valid_dir <- dirs$valid_dir

  if (isTRUE(opt$aggregate)) {
    # Pass configured TUs so aggregate_partitions only scans those directories;
    # NULL (= scan everything) only when scope is "all".
    tus_scope <- opt$simulator$scope$transmission_units %||% "all"
    tus_arg <- if (identical(tus_scope, "all")) {
      NULL
    } else {
      as.character(tus_scope)
    }
    out_fmt <- opt$simulator$output$format %||% "qs"
    aggregate_partitions(
      parts_dir,
      merged_dir,
      tus = tus_arg,
      out_fmt = out_fmt,
      return_long = FALSE
    )
  }

  if (isTRUE(opt$validate)) {
    if (is.null(opt$obs_panel) || !nzchar(opt$obs_panel)) {
      stop("--validate requires --obs-panel <path>.")
    }
    if (!file.exists(opt$obs_panel)) {
      stop("--obs-panel not found: ", opt$obs_panel)
    }

    # Frame build + all metrics are factored into .validate_run_from_opt() so
    # the --hindcast orchestrator scores each run through the identical path.
    # --manifest (a v2 forecast manifest) censors each country of a unit
    # mixing fit groups at its own firm horizon.
    hz_c <- NULL
    if (
      !is.null(opt$manifest) &&
        nzchar(opt$manifest) &&
        file.exists(opt$manifest)
    ) {
      man_v <- yaml::read_yaml(opt$manifest)
      if (length(man_v$units)) hz_c <- manifest_country_horizons(man_v, "firm")
    }
    .validate_run_from_opt(
      opt,
      merged_dir,
      valid_dir,
      opt$obs_panel,
      horizon_override = if (
        !is.null(opt$horizon_end) && nzchar(opt$horizon_end)
      ) {
        list(end = opt$horizon_end)
      } else {
        NULL
      },
      horizon_by_country = hz_c
    )
  }
  message("[07c_simulator_postprocess] done.")
}
