# 08c_smoothing_compile.R — compiler for the chunked spatPomp historical
#   smoothing: binds the 08b partitions and runs the global reductions
#
# Usage:
#   Rscript analysis/08c_smoothing_compile.R -C <cfg.yml>
#   Rscript analysis/08c_smoothing_compile.R -C <cfg.yml> --no_isweight
#
# Details:
#   Scans partitions for _SUCCESS sentinels, binds smoothing events, per-(TU,
#   rep) logLik and frame metadata, writes the per-TU forecast handoff cache,
#   then reduces globally: importance weights, coverage mask, Africa / region
#   / country ribbons, observed-data aggregation. Writes merged ribbons (.rds)
#   under <outdir>/sim_spatpomp/<run_tag>/merged/ and fig1/fig2/fig3 under
#   figures/scratchpad/<run_tag>/.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
  library(ggplot2)
  library(cowplot)
  library(qs2)
})

source("analysis/utils/smoothing_core.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("--config", "-C"),
    default = NULL,
    type = "character",
    help = "Path to YAML sim config"
  ),
  make_option(
    "--no_isweight",
    default = FALSE,
    action = "store_true",
    help = "Disable importance weighting (uniform weights 1/N)"
  )
)
args <- parse_args(OptionParser(option_list = option_list))

S <- smoothing_make_config(args)

parts_dir <- make_smooth_run_file(
  list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
  kind = "smooth_parts_dir"
)
merged_dir <- make_smooth_run_file(
  list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
  kind = "smooth_merged_dir"
)
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(S$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 3. Collect successful partitions ---------------------------------------------
# GLOBAL reductions that cannot be done per task (sections 4-7):
#   * importance weights (softmax over per-rep total logLik summed across TUs);
#   * surveillance-coverage mask (07b parity);
#   * Africa / region / country quantile-ribbon aggregation;
#   * observed-data aggregation + ribbon figures.
if (!dir.exists(parts_dir)) {
  stop("[smoothing_compile] parts dir not found: ", parts_dir)
}
sentinels <- partition_sentinels(parts_dir) # utils/paths.R
if (length(sentinels) == 0L) {
  stop(
    "[smoothing_compile] no _SUCCESS sentinels under ",
    parts_dir,
    " — run 08b_smoothing_spatpomp.R first."
  )
}
message(sprintf(
  "[smoothing_compile] %d successful partition(s) under %s",
  length(sentinels),
  parts_dir
))

# 3b. Assemble per-TU forecast handoff cache -----------------------------------
# (cache/hist_filter_tu=<TU>_*.qs)
# Bind the FULL per-replicate trajectories (states.qs, written per chunk by
# 08b_smoothing_spatpomp.R) into one flat record list per TU, in the format
# the forward forecaster's load_smoothing_distribution() globs. Done TU-by-TU so
# the large filter_traj arrays are never all resident at once. Skipped silently
# when no states.qs are present (older smoothing runs that predate this bridge).
cache_dir <- make_smooth_run_file(
  list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
  kind = "smooth_cache_dir"
)
state_files_all <- list.files(
  parts_dir,
  pattern = "^states\\.qs$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(state_files_all) == 0L) {
  message(
    "[smoothing_compile] no states.qs found under parts/ — skipping handoff cache ",
    "(re-run 08b_smoothing_spatpomp.R with the trajectory-persisting bridge to produce it)."
  )
} else {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # Map each states.qs to its TU via the parent "parts/tu=<TU>/chunk_*/" layout.
  tu_of <- function(p) sub("^tu=", "", basename(dirname(dirname(p))))
  by_tu <- split(state_files_all, vapply(state_files_all, tu_of, character(1)))
  start8 <- format(S$FILTER_START, "%Y%m%d")
  end8 <- format(S$FILTER_END, "%Y%m%d")
  for (tu in names(by_tu)) {
    recs <- unlist(
      lapply(by_tu[[tu]], qs2::qs_read),
      recursive = FALSE,
      use.names = FALSE
    )
    n_ok <- sum(vapply(recs, function(r) !is.null(r$filter_traj), logical(1)))
    # Schema token placed AFTER "_Np=<Np>_" (not between tu and Np) so the
    # existing tu-extraction regex in 07b_simulator_spatpomp.R
    # (`^hist_filter_tu=(.*?)_Np=.*$`, non-greedy up to the first "_Np=")
    # keeps capturing exactly the TU name -- inserting the token earlier
    # would fold it into the captured TU string.
    fn <- file.path(
      cache_dir,
      sprintf(
        "hist_filter_tu=%s_Np=%d_%s_nrep=%d_%s_%s.qs",
        tu,
        as.integer(S$NP),
        sm_cache_schema_token(S$COVARIATE_FEED %||% "latent"),
        length(recs),
        start8,
        end8
      )
    )
    qs2::qs_save(recs, fn)
    message(sprintf(
      "[smoothing_compile] handoff cache TU '%s': %d record(s) (%d with traj) → %s",
      tu,
      length(recs),
      n_ok,
      basename(fn)
    ))
  }
  message(sprintf(
    "[smoothing_compile] handoff caches for %d TU(s) → %s",
    length(by_tu),
    cache_dir
  ))
  # Cache-dir identity (spatpomp_core.R::smoothing_cache_identity(), over
  # smoothing_cache_meta()): the partition and fit-grouping registries (names +
  # content digests) this smoothing ran under, Np / nrep / schema, run_tag and
  # the TUs. 07b's load_smoothing_distribution() refuses a cache dir whose
  # identity differs from the forecast run's.
  # Optional identity fields (pop_year): written only when not
  # the legacy value, so a default smoothing's meta.yml is unchanged.
  # Estimated WHO scale: its digest, from 08a's meta file.
  if (
    identical(S$WHO_OBS[["scale_mode"]], "estimated") &&
      isTRUE(S$WHO_OBS$enable)
  ) {
    wsp <- smoothing_who_scale_paths(S)
    if (!file.exists(wsp$meta)) {
      stop(
        "[smoothing_compile] who_obs.scale_mode is estimated but ",
        wsp$meta,
        " is missing (written by 08a)."
      )
    }
    S$WHO_SCALE_DIGEST <- yaml::read_yaml(wsp$meta)$who_scale_digest
  }
  # Smoothing environment (table S13; provenance.R::sim_software_versions(),
  # via smoothing_core.R's source of provenance.R): informational, not identity.
  meta <- smoothing_cache_identity(
    S,
    tus = names(by_tu),
    optional = smoothing_cache_optional_fields(S),
    software = sim_software_versions()
  )
  yaml::write_yaml(meta, smoothing_cache_meta_path(cache_dir))
  message(sprintf(
    "[smoothing_compile] cache identity → %s (partition_by=%s%s)",
    smoothing_cache_meta_path(cache_dir),
    meta$partition_by,
    if (nzchar(meta$partition_source)) {
      paste0(", source=", meta$partition_source)
    } else {
      ""
    }
  ))
}

parts <- lapply(sentinels, function(s) {
  p <- file.path(dirname(s), "part.qs")
  if (!file.exists(p)) {
    warning("[smoothing_compile] missing part.qs next to ", s)
    return(NULL)
  }
  qs2::qs_read(p)
})
parts <- Filter(Negate(is.null), parts)

events_raw <- dplyr::bind_rows(lapply(parts, function(p) p$events))
# Absent (NULL element) on partitions written before this bridge existed, and
# empty on any partition where the mon track wasn't configured — both collapse
# to a 0-row bind, so the mon section below is skipped uniformly either way.
events_mon_raw <- dplyr::bind_rows(lapply(parts, function(p) p$events_mon))
rep_logliks <- dplyr::bind_rows(lapply(parts, function(p) p$rep_logliks))
gid_meta <- dplyr::bind_rows(lapply(parts, function(p) p$frame_meta)) %>%
  dplyr::mutate(gid = as.character(gid)) %>%
  dplyr::distinct(gid, .keep_all = TRUE)

if (is.null(events_raw) || nrow(events_raw) == 0L) {
  stop(
    "[smoothing_compile] no Yout=1 events across partitions — nothing to compile."
  )
}

# 4. Importance weights (GLOBAL across TUs by rep_id) --------------------------
rep_weights <- smoothing_rep_weights(rep_logliks, uniform = args$no_isweight)

# 5. Drop burn-in, build coverage mask, mask events ----------------------------
smooth_long <- events_raw %>% dplyr::filter(t_start >= S$FILTER_START)
if (nrow(smooth_long) == 0L) {
  stop("[smoothing_compile] all events fell in the burn-in window.")
}

sim_episode_check(
  preps = list(
    surv = S$SURV_PATH,
    size_data = S$SIZE_DATA_PATH,
    mon_data = S$MON_DATA_PATH
  ),
  what = "08c smoothing compile inputs",
  episode_key = S$EPISODE_KEY %||% "lp_season"
)
message("[smoothing_compile] reading survival data: ", S$SURV_PATH)
surv_dat <- readRDS(S$SURV_PATH)
surv_panel_obs <- build_obs_panel(surv_dat)

coverage_mask <- NULL
if (file.exists(S$AGG_TS_PATH)) {
  coverage_mask <- build_coverage_mask(
    # full LP -> gid member bridge when the panel's sidecar exists (utils/coverage.R)
    obs_panel = coverage_lp_bridge(surv_dat, S$SURV_PATH),
    aggregated_ts_path = S$AGG_TS_PATH,
    aggregated_adm_ts_path = if (file.exists(S$AGG_ADM_TS_PATH)) {
      S$AGG_ADM_TS_PATH
    } else {
      NULL
    },
    # the national-zero months of the run's threshold and rule (episodes.national_zero_*)
    national_zero_path = national_zero_months_file_of(S$OCC_PANEL %||% list())
  ) %>%
    dplyr::mutate(gid = as.character(gid), month = as.Date(month))
  message(sprintf(
    "[smoothing_compile] coverage mask: %d cells, %d gids, %s – %s",
    nrow(coverage_mask),
    dplyr::n_distinct(coverage_mask$gid),
    format(min(coverage_mask$month)),
    format(max(coverage_mask$month))
  ))

  # Restrict the comparison footprint to cells the filter could actually
  # constrain: inside surveillance coverage AND carrying a non-NA obs_Y (a
  # survival row). Cells in a coverage window but WITHOUT a survival row are
  # free-run by the filter (obs_Y = NA) and their fabricated cases must not be
  # counted against observations. Intersecting here propagates to the
  # aggregations (smooth_long, obs_occ, obs_cases_mon, obs_ongoing) via .mask_to_coverage(),
  # keeping the simulated and observed sides on an identical footprint.
  n_cov0 <- nrow(coverage_mask)
  coverage_mask <- coverage_mask %>%
    dplyr::inner_join(
      surv_panel_obs %>%
        dplyr::transmute(gid = as.character(gid), month = t_start),
      by = c("gid", "month")
    ) %>%
    dplyr::distinct(gid, month)
  message(sprintf(
    "[smoothing_compile] coverage ∩ obs_Y: %d -> %d cells (%.0f%% retained)",
    n_cov0,
    nrow(coverage_mask),
    100 * nrow(coverage_mask) / max(n_cov0, 1L)
  ))
} else {
  warning(
    "[smoothing_compile] aggregated_ts not found — coverage mask disabled."
  )
}

.mask_to_coverage <- function(df, label) {
  if (is.null(coverage_mask)) {
    return(df)
  }
  n0 <- nrow(df)
  out <- df %>%
    dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::semi_join(coverage_mask, by = c("gid", "t_start" = "month"))
  message(sprintf(
    "[smoothing_compile] mask %-10s: %d -> %d rows (%.0f%% retained)",
    label,
    n0,
    nrow(out),
    100 * nrow(out) / max(n0, 1L)
  ))
  out
}

smooth_long_unmasked <- smooth_long
smooth_long <- .mask_to_coverage(smooth_long, "smooth_long")
if (nrow(smooth_long) == 0L) {
  stop("[smoothing_compile] smooth_long empty after coverage mask.")
}

message(sprintf(
  "[smoothing_compile] smooth_long: %d events, %d reps, %d months",
  nrow(smooth_long),
  length(unique(smooth_long$rep_id)),
  length(unique(smooth_long$t_start))
))

# The full filtered month grid, not only the months carrying a masked simulated
# onset: a month with no simulated onset is a zero, not a missing month.
all_months_vec <- seq(
  lubridate::floor_date(S$FILTER_START, "month"),
  lubridate::floor_date(S$FILTER_END - 1, "month"),
  by = "month"
)

# Monthly case-count (mon) track — opt-in. NULL/0-row when the mon track
# wasn't configured for this smoothing run (see events_mon_raw comment above);
# `has_mon` gates every downstream mon block so the rest of the script runs
# unchanged for pre-mon runs and runs without a monthly_fit configured.
has_mon <- !is.null(events_mon_raw) && nrow(events_mon_raw) > 0L
if (has_mon) {
  smooth_long_mon <- events_mon_raw %>% dplyr::filter(t_start >= S$FILTER_START)
  smooth_long_mon_unmasked <- smooth_long_mon
  smooth_long_mon <- .mask_to_coverage(smooth_long_mon, "smooth_long_mon")
  has_mon <- nrow(smooth_long_mon) > 0L
  if (has_mon) {
    message(sprintf(
      "[smoothing_compile] smooth_long_mon: %d events, %d reps, %d months",
      nrow(smooth_long_mon),
      length(unique(smooth_long_mon$rep_id)),
      length(unique(smooth_long_mon$t_start))
    ))
  } else {
    message(
      "[smoothing_compile] smooth_long_mon empty after coverage mask — skipping mon ribbons."
    )
  }
} else {
  message(
    "[smoothing_compile] no mon-track events across partitions — skipping mon ribbons."
  )
}

# 6. Aggregate ribbons (Africa / region / country) -----------------------------
smooth_long_africa <- dplyr::mutate(smooth_long, africa = "Africa")
message("[smoothing_compile] aggregating ribbons...")
ribbon_africa <- agg_and_ribbon(
  smooth_long_africa,
  "africa",
  rep_weights,
  all_months_vec
)
ribbon_region <- agg_and_ribbon(
  smooth_long,
  "AFRO_region",
  rep_weights,
  all_months_vec
)
ribbon_country <- agg_and_ribbon(
  smooth_long,
  "country",
  rep_weights,
  all_months_vec
)

if (has_mon) {
  smooth_long_mon_africa <- dplyr::mutate(smooth_long_mon, africa = "Africa")
  ribbon_africa_mon <- agg_and_ribbon_mon(
    smooth_long_mon_africa,
    "africa",
    rep_weights,
    all_months_vec
  )
  ribbon_region_mon <- agg_and_ribbon_mon(
    smooth_long_mon,
    "AFRO_region",
    rep_weights,
    all_months_vec
  )
  ribbon_country_mon <- agg_and_ribbon_mon(
    smooth_long_mon,
    "country",
    rep_weights,
    all_months_vec
  )
} else {
  ribbon_africa_mon <- ribbon_region_mon <- ribbon_country_mon <- NULL
}

# Unmasked companions of the ribbons above -- same aggregation functions, same
# all_months_vec grid, run against the pre-.mask_to_coverage() copies captured
# in section 5 (smooth_long_unmasked / smooth_long_mon_unmasked). NULL when no
# coverage mask was ever built (masked == unmasked already, so skip the
# redundant computation) -- panel_pair()/ribbon_plot() treat a NULL
# ribbon_unmasked as "no overlay", preserving pre-mask-run figures unchanged.
if (!is.null(coverage_mask)) {
  smooth_long_unmasked_africa <- dplyr::mutate(
    smooth_long_unmasked,
    africa = "Africa"
  )
  ribbon_africa_unmasked <- agg_and_ribbon(
    smooth_long_unmasked_africa,
    "africa",
    rep_weights,
    all_months_vec
  )
  ribbon_region_unmasked <- agg_and_ribbon(
    smooth_long_unmasked,
    "AFRO_region",
    rep_weights,
    all_months_vec
  )
  ribbon_country_unmasked <- agg_and_ribbon(
    smooth_long_unmasked,
    "country",
    rep_weights,
    all_months_vec
  )

  if (has_mon) {
    smooth_long_mon_unmasked_africa <- dplyr::mutate(
      smooth_long_mon_unmasked,
      africa = "Africa"
    )
    ribbon_africa_mon_unmasked <- agg_and_ribbon_mon(
      smooth_long_mon_unmasked_africa,
      "africa",
      rep_weights,
      all_months_vec
    )
    ribbon_region_mon_unmasked <- agg_and_ribbon_mon(
      smooth_long_mon_unmasked,
      "AFRO_region",
      rep_weights,
      all_months_vec
    )
    ribbon_country_mon_unmasked <- agg_and_ribbon_mon(
      smooth_long_mon_unmasked,
      "country",
      rep_weights,
      all_months_vec
    )
  } else {
    ribbon_africa_mon_unmasked <- ribbon_region_mon_unmasked <- ribbon_country_mon_unmasked <- NULL
  }
} else {
  ribbon_africa_unmasked <- ribbon_region_unmasked <- ribbon_country_unmasked <- NULL
  ribbon_africa_mon_unmasked <- ribbon_region_mon_unmasked <- ribbon_country_mon_unmasked <- NULL
}

# 6c. Ongoing-outbreak ("currently active episode") ribbons --------------------
# Unlike the sparse Yout onset table (smooth_long, one row per episode's onset
# month only), "ongoing outbreaks" needs every month an episode is running.
# That information lives in the raw per-replicate filter trajectory (the
# duration model's `active` state), not in the compiled events table, so this
# re-reads states.qs (already scanned once in section 3b for the forecast
# handoff cache) and reconstructs it via .smoothing_active_matrix(). Degrades
# gracefully to an all-zero ribbon when the duration track wasn't configured
# for this run (the `active` slot is pinned at 0 -- see .ft_rows_or_zero()).
if (length(state_files_all) == 0L) {
  message(
    "[smoothing_compile] no states.qs found — skipping ongoing-outbreak ribbons."
  )
  ribbon_africa_active <- ribbon_region_active <- ribbon_country_active <- NULL
  ribbon_africa_active_unmasked <- ribbon_region_active_unmasked <- ribbon_country_active_unmasked <- NULL
  obs_ongoing_africa <- obs_ongoing_region <- obs_ongoing_country <- NULL
} else {
  message(
    "[smoothing_compile] reconstructing ongoing-outbreak counts from states.qs..."
  )
  tu_country_region <- gid_meta %>% dplyr::distinct(country, AFRO_region)

  # Compute both the masked and unmasked active-episode counts in a single
  # pass over states.qs (avoids re-reading the trajectory files twice): each
  # record contributes an n_active_unmasked column (raw act_m) alongside the
  # existing n_active column (act_m * mask_mat when a coverage mask exists).
  active_long_both <- purrr::map_dfr(names(by_tu), function(tu) {
    recs <- unlist(
      lapply(by_tu[[tu]], qs2::qs_read),
      recursive = FALSE,
      use.names = FALSE
    )
    purrr::map_dfr(recs, function(rec) {
      act_m <- .smoothing_active_matrix(rec)
      if (is.null(act_m)) {
        return(NULL)
      }
      # covariate_feed: detected -> the detected-activity field (NULL otherwise)
      det_m <- .smoothing_detected_matrix(rec)
      act_m_masked <- act_m
      det_m_masked <- det_m
      if (!is.null(coverage_mask)) {
        cov_sub <- coverage_mask %>% dplyr::filter(gid %in% rec$gids)
        gid_idx <- match(cov_sub$gid, rec$gids)
        month_idx <- match(cov_sub$month, rec$filt_months)
        ok <- !is.na(gid_idx) & !is.na(month_idx)
        mask_mat <- matrix(
          FALSE,
          nrow = length(rec$gids),
          ncol = length(rec$filt_months)
        )
        mask_mat[cbind(gid_idx[ok], month_idx[ok])] <- TRUE
        act_m_masked <- act_m * mask_mat
        if (!is.null(det_m)) det_m_masked <- det_m * mask_mat
      }
      # Collapse units to their country (rows of act_m are rec$gids). The TU is
      # not necessarily a country: under partition_by: enso_region it is an
      # ENSO group spanning several countries.
      gid_country <- gid_meta$country[match(
        as.character(rec$gids),
        gid_meta$gid
      )]
      keep <- !is.na(gid_country)
      if (!any(keep)) {
        return(NULL)
      }
      by_country <- function(m) {
        rowsum(m[keep, , drop = FALSE], gid_country[keep])
      }
      a_m <- by_country(act_m_masked)
      a_u <- by_country(act_m)
      n_m <- length(rec$filt_months)
      out_r <- tibble::tibble(
        rep_id = rec$rep_id,
        t_start = rep(rec$filt_months, each = nrow(a_m)),
        country = rep(rownames(a_m), times = n_m),
        n_active = as.vector(a_m),
        n_active_unmasked = as.vector(a_u)
      )
      if (!is.null(det_m)) {
        out_r$n_detected <- as.vector(by_country(det_m_masked))
        out_r$n_detected_unmasked <- as.vector(by_country(det_m))
      }
      out_r
    })
  }) %>%
    dplyr::filter(t_start >= S$FILTER_START) %>%
    dplyr::left_join(tu_country_region, by = "country") %>%
    dplyr::filter(!is.na(AFRO_region), t_start %in% all_months_vec)

  active_long <- active_long_both %>%
    dplyr::select(
      -n_active_unmasked,
      -dplyr::any_of(c("n_detected", "n_detected_unmasked"))
    )

  active_long_africa <- dplyr::mutate(active_long, africa = "Africa")
  ribbon_africa_active <- agg_and_ribbon_active(
    active_long_africa,
    "africa",
    rep_weights,
    all_months_vec
  )
  ribbon_region_active <- agg_and_ribbon_active(
    active_long,
    "AFRO_region",
    rep_weights,
    all_months_vec
  )
  ribbon_country_active <- agg_and_ribbon_active(
    active_long,
    "country",
    rep_weights,
    all_months_vec
  )

  if (!is.null(coverage_mask)) {
    active_long_unmasked <- active_long_both %>%
      dplyr::select(
        -n_active,
        -dplyr::any_of(c("n_detected", "n_detected_unmasked"))
      ) %>%
      dplyr::rename(n_active = n_active_unmasked)

    active_long_unmasked_africa <- dplyr::mutate(
      active_long_unmasked,
      africa = "Africa"
    )
    ribbon_africa_active_unmasked <- agg_and_ribbon_active(
      active_long_unmasked_africa,
      "africa",
      rep_weights,
      all_months_vec
    )
    ribbon_region_active_unmasked <- agg_and_ribbon_active(
      active_long_unmasked,
      "AFRO_region",
      rep_weights,
      all_months_vec
    )
    ribbon_country_active_unmasked <- agg_and_ribbon_active(
      active_long_unmasked,
      "country",
      rep_weights,
      all_months_vec
    )
  } else {
    ribbon_africa_active_unmasked <- ribbon_region_active_unmasked <- ribbon_country_active_unmasked <- NULL
  }

  # Observed ongoing-outbreak count: expand each episode (outbreak_size_analysis_data's
  # TL/TR) across every month it was running, then count episode-months per
  # group -- NOT distinct gids, since two concurrent episodes at the same gid
  # both count as 2 ongoing outbreaks that month.
  obs_ongoing_raw <- readRDS(S$SIZE_DATA_PATH) %>%
    dplyr::mutate(
      gid = as.character(gid),
      TL = as.Date(TL),
      TR = as.Date(TR)
    ) %>%
    dplyr::filter(!is.na(gid), !is.na(TL), !is.na(TR)) %>%
    dplyr::mutate(
      month_start = lubridate::floor_date(TL, "month"),
      month_end = lubridate::floor_date(TR, "month")
    ) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      t_start = list(seq(month_start, month_end, by = "1 month"))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(t_start) %>%
    dplyr::transmute(row_id, gid, t_start) %>%
    dplyr::filter(t_start >= S$FILTER_START, t_start < S$FILTER_END) %>%
    .mask_to_coverage("obs_ongoing") %>%
    dplyr::filter(t_start %in% all_months_vec) %>%
    dplyr::left_join(gid_meta, by = "gid") %>%
    dplyr::filter(!is.na(AFRO_region))

  obs_ongoing_africa <- obs_ongoing_raw %>%
    dplyr::group_by(t_start) %>%
    dplyr::summarise(n_active = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(africa = "Africa")
  obs_ongoing_region <- obs_ongoing_raw %>%
    dplyr::group_by(t_start, AFRO_region) %>%
    dplyr::summarise(n_active = dplyr::n(), .groups = "drop")
  obs_ongoing_country <- obs_ongoing_raw %>%
    dplyr::group_by(t_start, country) %>%
    dplyr::summarise(n_active = dplyr::n(), .groups = "drop")
}

# 6d. DETECTED-series ribbons (covariate_feed: detected) -----------------------
# The latent ribbons above are what the model believes happened; these are what
# the surveillance record would have shown under the same delta_fn/delta_fp,
# and are the like-for-like comparison with the observed series. Built with the
# same aggregation functions on the detected subsets -- onsets: events
# whose episode was detected; monthly cases: months of detected episodes;
# ongoing: the detected-activity field. NULL (and has_detected = FALSE) for a
# latent-feed run, whose event tables carry no `detected` column.
.ribbons_3 <- function(long, fn) {
  if (is.null(long) || nrow(long) == 0L) {
    return(NULL)
  }
  list(
    africa = fn(
      dplyr::mutate(long, africa = "Africa"),
      "africa",
      rep_weights,
      all_months_vec
    ),
    region = fn(long, "AFRO_region", rep_weights, all_months_vec),
    country = fn(long, "country", rep_weights, all_months_vec)
  )
}
.det_only <- function(long) {
  if (is.null(long) || !"detected" %in% names(long)) {
    NULL
  } else {
    dplyr::filter(long, detected == 1L)
  }
}
has_detected <- "detected" %in% names(smooth_long)
ribbons_detected <- if (has_detected) {
  list(
    onsets = .ribbons_3(.det_only(smooth_long), agg_and_ribbon),
    onsets_unmasked = if (!is.null(coverage_mask)) {
      .ribbons_3(.det_only(smooth_long_unmasked), agg_and_ribbon)
    },
    mon = if (has_mon) {
      .ribbons_3(.det_only(smooth_long_mon), agg_and_ribbon_mon)
    },
    mon_unmasked = if (has_mon && !is.null(coverage_mask)) {
      .ribbons_3(.det_only(smooth_long_mon_unmasked), agg_and_ribbon_mon)
    },
    active = if (
      exists("active_long_both") && "n_detected" %in% names(active_long_both)
    ) {
      .ribbons_3(
        active_long_both %>%
          dplyr::transmute(
            rep_id,
            t_start,
            country,
            AFRO_region,
            n_active = n_detected
          ),
        agg_and_ribbon_active
      )
    },
    active_unmasked = if (
      !is.null(coverage_mask) &&
        exists("active_long_both") &&
        "n_detected_unmasked" %in% names(active_long_both)
    ) {
      .ribbons_3(
        active_long_both %>%
          dplyr::transmute(
            rep_id,
            t_start,
            country,
            AFRO_region,
            n_active = n_detected_unmasked
          ),
        agg_and_ribbon_active
      )
    }
  )
} else {
  NULL
}
if (has_detected) {
  message(
    "[smoothing_compile] detected-series ribbons built (covariate_feed: detected)."
  )
}

# 7. Observed data aggregation -------------------------------------------------
obs_occ_base <- surv_panel_obs %>%
  .mask_to_coverage("obs_occ") %>%
  dplyr::filter(t_start %in% all_months_vec) %>%
  dplyr::left_join(gid_meta, by = "gid") %>%
  dplyr::filter(!is.na(AFRO_region))

obs_occ_africa <- obs_occ_base %>%
  dplyr::group_by(t_start) %>%
  dplyr::summarise(n_outbreaks = sum(obs_Y, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(africa = "Africa")
obs_occ_region <- obs_occ_base %>%
  dplyr::group_by(t_start, AFRO_region) %>%
  dplyr::summarise(n_outbreaks = sum(obs_Y, na.rm = TRUE), .groups = "drop")
obs_occ_country <- obs_occ_base %>%
  dplyr::group_by(t_start, country) %>%
  dplyr::summarise(n_outbreaks = sum(obs_Y, na.rm = TRUE), .groups = "drop")

# True monthly-resolved observed cases: outbreaks_monthly_data_base-<baseline>.rds's
# `cases_month`, spread across each episode's active months -- the estimand the
# mon track is fit to reproduce (NOT the whole-episode total lumped onto the
# onset month: the two conventions' monthly series correlate only ~0.48).
# Put on the same footprint as obs_occ: coverage mask, then restricted to the
# simulated frame (gid_meta join) using the frame's AFRO_region/country, and
# to months the simulation produced. Without the frame restriction, a sizeable
# share of observed case mass sits on gids with no simulated counterpart, so
# the observed side would be inflated relative to a band the model can never
# match. Each (gid, month) is first collapsed with the filter's footprint rule
# (sum within an admin level, max across levels; smoothing_obs_cases_mon()), so
# nested reports of one gid at two levels are not counted twice.
obs_cases_mon_africa <- obs_cases_mon_region <- obs_cases_mon_country <- NULL
if (has_mon && file.exists(S$MON_DATA_PATH)) {
  obs_cases_mon <- smoothing_obs_cases_mon(
    readRDS(S$MON_DATA_PATH),
    gid_meta,
    S$FILTER_START,
    S$FILTER_END,
    all_months_vec,
    mask = .mask_to_coverage
  )
  obs_cases_mon_africa <- obs_cases_mon$africa
  obs_cases_mon_region <- obs_cases_mon$region
  obs_cases_mon_country <- obs_cases_mon$country
} else if (has_mon) {
  warning(
    "[smoothing_compile] S$MON_DATA_PATH not found (",
    S$MON_DATA_PATH,
    ") — the mon-track cases panel is skipped in fig1/fig2/fig3."
  )
}

# 7b. WHO dashboard comparison (national scale) --------------------------------
# Unmasked smoothed monthly cases per country vs the WHO dashboard national
# counts over the who_obs window. Built whenever the monthly WHO file exists,
# whether or not the filter scored it (the comparison is also the diagnostic
# for runs without the WHO term). who_f (covered population share) travels in
# the partitions' frame_meta when the run had who_obs on; otherwise the band
# covers the simulated gids only (who_f = 1) and says so.
who_cmp <- NULL
who_path <- S$WHO_OBS$path %||% make_who_monthly_file()
if (has_mon && file.exists(who_path)) {
  who_f <- if ("who_f" %in% names(gid_meta)) {
    tapply(gid_meta$who_f, gid_meta$country, function(x) x[1])
  } else {
    NULL
  }
  # Estimated per-country scale: the bands are on the scale
  # the filter scored (who_multiplier()); NULL = fixed who_obs.scale.
  who_scale <- if ("who_scale" %in% names(gid_meta)) {
    tapply(gid_meta$who_scale, gid_meta$country, function(x) x[1])
  } else {
    NULL
  }
  if (is.null(who_f)) {
    message(
      "[smoothing_compile] WHO comparison: no who_f in the partitions (who_obs off) -- ",
      "bands cover the simulated gids only."
    )
  }
  who_cmp <- who_vs_smoothing(
    smooth_long_mon_unmasked,
    readRDS(who_path),
    rep_weights,
    start = S$WHO_OBS$start %||% as.Date("2023-01-01"),
    end = S$WHO_OBS$end %||% as.Date("2025-01-01"),
    who_f = who_f,
    # the filter's scale; NULL at fixed scale 1 keeps the
    # comparison table unchanged (1 / who_f either way)
    who_scale = if (!is.null(who_scale)) {
      who_scale
    } else if (
      isTRUE(S$WHO_OBS$enable) &&
        !isTRUE(S$WHO_OBS$scale == 1)
    ) {
      S$WHO_OBS$scale
    }
  )
  if (!is.null(who_cmp)) {
    print(as.data.frame(who_cmp$summary), row.names = FALSE, digits = 3)
    g <- plot_who_vs_smoothing(
      list(smoothing = who_cmp$ribbon),
      who_cmp$obs,
      title = sprintf("WHO dashboard vs smoothing: %s", S$run_tag)
    )
    n_ct <- length(unique(who_cmp$obs$country))
    ggplot2::ggsave(
      file.path(S$OUT_DIR, "fig4_who_country.png"),
      g,
      width = 14,
      height = 2.4 * ceiling(n_ct / 4) + 1.2,
      dpi = 150,
      bg = "white"
    )
    message("[smoothing_compile] fig4_who_country saved.")
  }
}

# 8. Persist merged ribbons ----------------------------------------------------
saveRDS(
  list(
    ribbon_africa = ribbon_africa,
    ribbon_region = ribbon_region,
    ribbon_country = ribbon_country,
    rep_weights = rep_weights,
    obs_occ_africa = obs_occ_africa,
    obs_occ_region = obs_occ_region,
    obs_occ_country = obs_occ_country,
    # mon-track cases ribbons (NULL when the mon track wasn't configured for
    # this run), compared against obs_cases_mon_* (true per-month cases_month
    # panel).
    has_mon = has_mon,
    ribbon_africa_mon = ribbon_africa_mon,
    ribbon_region_mon = ribbon_region_mon,
    ribbon_country_mon = ribbon_country_mon,
    obs_cases_mon_africa = obs_cases_mon_africa,
    obs_cases_mon_region = obs_cases_mon_region,
    obs_cases_mon_country = obs_cases_mon_country,
    # Ongoing-outbreak ("currently active episode") ribbons, reconstructed from
    # the duration model's `active` state (all-zero/NULL-degrading when the
    # duration track wasn't configured -- see section 6c above).
    ribbon_africa_active = ribbon_africa_active,
    ribbon_region_active = ribbon_region_active,
    ribbon_country_active = ribbon_country_active,
    obs_ongoing_africa = obs_ongoing_africa,
    obs_ongoing_region = obs_ongoing_region,
    obs_ongoing_country = obs_ongoing_country,
    # Unmasked companions (pre-.mask_to_coverage()) of the ribbons above, for the
    # unmasked-trajectory overlay in fig1/fig2/fig3. NULL when no coverage
    # mask was ever built for this run (see section 6).
    ribbon_africa_unmasked = ribbon_africa_unmasked,
    ribbon_region_unmasked = ribbon_region_unmasked,
    ribbon_country_unmasked = ribbon_country_unmasked,
    ribbon_africa_mon_unmasked = ribbon_africa_mon_unmasked,
    ribbon_region_mon_unmasked = ribbon_region_mon_unmasked,
    ribbon_country_mon_unmasked = ribbon_country_mon_unmasked,
    ribbon_africa_active_unmasked = ribbon_africa_active_unmasked,
    ribbon_region_active_unmasked = ribbon_region_active_unmasked,
    ribbon_country_active_unmasked = ribbon_country_active_unmasked,
    # Detected-series ribbons (section 6d): list(onsets, mon, active, *_unmasked),
    # each list(africa, region, country). NULL for latent-feed runs.
    has_detected = has_detected,
    ribbons_detected = ribbons_detected,
    # WHO dashboard comparison (section 7b): list(ribbon, obs, summary), NULL
    # when the WHO monthly file or the mon track is absent.
    who_country = who_cmp
  ),
  file.path(merged_dir, "smoothing_ribbons.rds")
)
message(
  "[smoothing_compile] merged ribbons → ",
  file.path(merged_dir, "smoothing_ribbons.rds")
)

# 9. Figures -------------------------------------------------------------------
# fig1/fig2/fig3 fold in two additional panels when available: the mon-track
# monthly-cases ribbon (has_mon, compared against the observed per-month
# `cases_month` series obs_cases_mon_*, section 7) and the ongoing-outbreak
# ("currently active episode") ribbon reconstructed from the duration model's
# `active` state (section 6c). Both degrade gracefully to a skipped panel /
# all-zero ribbon when the corresponding track wasn't configured for a given
# run. The mon track is the only case-count channel (the size track is not
# simulated).
message("[smoothing_compile] generating figures...")

# fig1/fig2/fig3 fold in the mon-track cases panel (when a monthly fit
# was configured, has_mon) and the ongoing-outbreak panel (when the duration
# track produced a non-NULL active-state ribbon) via panel_pair() -- both
# optional pairs default to NULL/skip inside panel_pair(), so a pre-mon /
# pre-duration run renders the occurrence panel alone.

# View clamp for the Africa-wide CASES panel only (fig1_africa).
# A single unit-month -- SSD_7_2_1's 2023-04 onset -- carries a q975 of ~61k and
# flattens the entire 2011-2022 series into the axis floor. coord_cartesian()
# zooms rather than filters, so the clamped ribbons run off the top edge and no
# rows are dropped. The region/country figures keep their free scales.
AFRICA_CASES_YLIM <- c(0, 35000)

fig1 <- panel_pair(
  ribbon_africa,
  obs_occ_africa,
  "africa",
  "Africa",
  ribbon_mon = if (has_mon) ribbon_africa_mon else NULL,
  obs_cases_mon = if (has_mon) obs_cases_mon_africa else NULL,
  ribbon_active = ribbon_africa_active,
  obs_active = obs_ongoing_africa,
  ylim_cases = AFRICA_CASES_YLIM,
  ribbon_unmasked = ribbon_africa_unmasked,
  ribbon_mon_unmasked = if (has_mon) ribbon_africa_mon_unmasked else NULL,
  ribbon_active_unmasked = ribbon_africa_active_unmasked
)
cowplot::save_plot(
  file.path(S$OUT_DIR, "fig1_africa.pdf"),
  fig1,
  base_width = 14,
  base_height = 12
)
cowplot::save_plot(
  file.path(S$OUT_DIR, "fig1_africa.png"),
  fig1,
  base_width = 14,
  base_height = 12,
  bg = "white"
)
message("[smoothing_compile] fig1_africa saved.")

regions <- sort(unique(ribbon_region$AFRO_region))
plots_r <- lapply(regions, function(rg) {
  panel_pair(
    ribbon_region,
    obs_occ_region,
    "AFRO_region",
    rg,
    ribbon_mon = if (has_mon) ribbon_region_mon else NULL,
    obs_cases_mon = if (has_mon) obs_cases_mon_region else NULL,
    ribbon_active = ribbon_region_active,
    obs_active = obs_ongoing_region,
    ribbon_unmasked = ribbon_region_unmasked,
    ribbon_mon_unmasked = if (has_mon) ribbon_region_mon_unmasked else NULL,
    ribbon_active_unmasked = ribbon_region_active_unmasked
  )
})
fig2 <- cowplot::plot_grid(plotlist = plots_r, ncol = length(regions))
cowplot::save_plot(
  file.path(S$OUT_DIR, "fig2_region.pdf"),
  fig2,
  base_width = 4.5 * length(regions),
  base_height = 12
)
cowplot::save_plot(
  file.path(S$OUT_DIR, "fig2_region.png"),
  fig2,
  base_width = 4.5 * length(regions),
  base_height = 12,
  bg = "white"
)
message("[smoothing_compile] fig2_region saved.")

countries <- sort(unique(ribbon_country$country))
per_page <- 6L
n_pages <- ceiling(length(countries) / per_page)
for (pg in seq_len(n_pages)) {
  ctry_pg <- countries[seq(
    (pg - 1L) * per_page + 1L,
    min(pg * per_page, length(countries))
  )]
  plots_pg <- lapply(ctry_pg, function(ct) {
    panel_pair(
      ribbon_country,
      obs_occ_country,
      "country",
      ct,
      ribbon_mon = if (has_mon) ribbon_country_mon else NULL,
      obs_cases_mon = if (has_mon) obs_cases_mon_country else NULL,
      ribbon_active = ribbon_country_active,
      obs_active = obs_ongoing_country,
      ribbon_unmasked = ribbon_country_unmasked,
      ribbon_mon_unmasked = if (has_mon) ribbon_country_mon_unmasked else NULL,
      ribbon_active_unmasked = ribbon_country_active_unmasked
    )
  })
  fig3_pg <- cowplot::plot_grid(plotlist = plots_pg, ncol = length(ctry_pg))
  fn_base <- file.path(S$OUT_DIR, sprintf("fig3_country_p%02d", pg))
  cowplot::save_plot(
    paste0(fn_base, ".pdf"),
    fig3_pg,
    base_width = 5 * length(ctry_pg),
    base_height = 12
  )
  cowplot::save_plot(
    paste0(fn_base, ".png"),
    fig3_pg,
    base_width = 5 * length(ctry_pg),
    base_height = 12,
    bg = "white"
  )
}
message(sprintf(
  "[smoothing_compile] fig3_country saved (%d page(s)).",
  n_pages
))

message("[smoothing_compile] done. Figures in: ", S$OUT_DIR)
