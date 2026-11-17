# results_core.R — shared kernel of the manuscript results phase: registries,
#   manifest/config accessors, artefact freshness and provenance contract
#
# Provides:
#   results_inputs(), results_config_path(), resolve_results_sets()
#   load_manifest(), manifest_groups_all(), results_run_configs()
#   results_write(), results_read_if_fresh(), results_require()
#                            write / read derived tables under the stamp
#   load_sim_run_monthly(), load_sim_run_stitched(), load_arm_run()
#   results_group_meta(), group_display_meta(), relabel_groups()
#   summarise_reps(), adm1_prob_stats()   replicate summaries
#   TRACK_*, GROUP_* registries
#
# Details:
#   Single-sources what both halves of the phase need:
#   09b_postprocess_results.R derives the tables from fits and simulator
#   partitions into generated_data/results/<tag>/; 10_manuscript_figures.R
#   plots those tables and nothing else (no path to a fit or a particle
#   filter). A `_core` kernel: it sources its own dependency closure.

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
  library(qs2)
})

# Why this module sources its own dependency closure, when the utils modules
# deliberately never source each other: the `_core` libraries are the
# documented exception (see analysis/utils/smoothing_core.R, which sources its
# modules with a per-line justification). A `_core` module is a kernel meant
# to be sourced standalone by several pipeline scripts, so it owns its closure
# rather than pushing it onto every caller.
source("analysis/utils/base.R")
source("analysis/utils/config.R") # .flag_val() -- required by results_config_path() to recover -C for the re-run hint
source("analysis/utils/paths.R") # make_results_run_file() -- the artefact-path builder both halves call
source("analysis/utils/provenance.R") # attach_results_provenance()/get_results_provenance() -- the stamp carried on every derived table
source("analysis/utils/enso.R") # ENSO_ANALOG_SETS + enso_group_colours() -- required by resolve_results_sets() and GROUP_COLOURS
source("analysis/utils/gid.R") # flatten_gid() -- required by load_sim_run_monthly()
source("analysis/utils/sim_io.R") # censor_sim_by_country()/manifest_*() -- per-country firm horizons of partition units (load_sim_run_monthly(), results_derive.R)
source("analysis/utils/figures_cb.R") # COMPOSITE_CONDS + the shared ENSO crossbasis extractors (extract_adm1_irr() ...)
source("analysis/utils/footprint.R") # finest_level_totals() -- observed monthly cases per gid of the event attribution (results_derive.R::event_observed_panel())

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ---------------------------------------------------------------------------
# Registries: tracks, groups, metrics
# ---------------------------------------------------------------------------

#' Per-track registry (shared across all manuscript figures). Names throughout
#' both scripts ("occurrence", "duration", "monthly") are the literal `track`
#' values dispatched by select_best_lag() (analysis/utils/loo_lags.R -- it only
#' special-cases "occurrence", so these names are usable as-is). The other
#' three tables key off the SAME per-track filesystem/archive conventions used
#' elsewhere in the repo (04_fit_<TRACK_CODE>_brms_cb.R, make_run_file()-family
#' kind prefixes, <lag_loo_batch>_<LOO_DIR_STUB>_grouped_enso_singlelag/
#' batches). Kept as one registry rather than split across the two scripts:
#' only TRACK_LABEL is plot-side, and splitting the table across files to save
# The ADVERSE direction of each track's ENSO log-IRR: +1 where a positive
# log-IRR means more cholera (occurrence onsets, monthly cases), -1 where it
# means less (duration models the CESSATION hazard, so a positive log-IRR is a
# shorter outbreak and the adverse direction is negative). Used by the ENSO
# exposure tables (results_derive.R::derive_report_enso_exposure) to put every
# track on one "adverse / protective" axis before combining them.
TRACK_ADVERSE_SIGN <- c(occurrence = 1, duration = -1, monthly = 1)
#' one label vector costs more than it saves.
#'
#' The size track is retired downstream (the monthly track's
#' `cases_month` is the sole case estimand) and has no entry here, so a
#' `results.scope.fig1_tracks` naming "size" fails on the first registry lookup.
TRACK_CODE <- c(occurrence = "occ", duration = "dur", monthly = "mon")
TRACK_KIND_PREFIX <- c(occurrence = "", duration = "dur_", monthly = "mon_")
TRACK_LOO_DIR_STUB <- c(occurrence = "occ", duration = "dur", monthly = "mon")
TRACK_LABEL <- c(
  occurrence = "Occurrence",
  duration = "Duration",
  monthly = "Monthly"
)

# Group display tables, read from the enso_grouping sidecar
# (analysis/groupings/enso_grouping.yml via utils/config.R::load_grouping_meta()).
# The names are kept: phase-09 artefacts store GROUP_SHORT as data (`scope`,
# `region`, the smoothing `by_scope_metric` keys), so the values are pinned by
# tests/test_grouping_meta.R. For another fit
# registry use results_group_meta(ctx), which reads the manifest's `grouping`.
.ENSO_GROUP_META <- load_grouping_meta("enso_grouping")
GROUP_LABELS <- .ENSO_GROUP_META$label
GROUP_SHORT <- .ENSO_GROUP_META$short
GROUP_ORDER <- .ENSO_GROUP_META$order
# Display-only short group names ("Southern" -> "South"); GROUP_SHORT itself
# stays as stored data and is only relabelled at plot time. Shared by
# 10_manuscript_figures.R and notebooks/hindcast_skill.qmd.
GROUP_SHORT_FIG <- .ENSO_GROUP_META$short_fig
# ENSO-region palette, deliberately red/blue-free so it never collides with
# the figures' red/blue ENSO-phase colours and diverging effect scales.
GROUP_COLOURS <- group_colour_map(.ENSO_GROUP_META)

#' results_group_meta
#'
#' Display metadata of the fit registry a results context reads: the
#' manifest's `grouping` (forecast_horizons.yml), else enso_grouping. Returns
#' the load_grouping_meta() list.
results_group_meta <- function(ctx = NULL) {
  src <- if (!is.null(ctx) && !is.null(ctx$MAN$grouping)) {
    ctx$MAN$grouping
  } else {
    "enso_grouping"
  }
  load_grouping_meta(src)
}

# ---------------------------------------------------------------------------
# Group display labels
#
# Group IDS (enso_east, ...) never change: they key the fit files, the
# configs, the manifests and every 09 artefact (enso_group columns; `region` /
# `scope` columns hold the sidecar's stored `short`). What a figure or a table
# PRINTS is resolved at display time in 10 / 11 from results.style.group_labels:
#   registry  the sidecar's label / short / short_fig (the default)
#   lag       a class letter ordered by the deployed OCCURRENCE lag, longest
#             lead first (ties in sidecar order): "Class A (ENSO lead 6)" /
#             short "A −6", "Class D (ENSO lag +7)" / "D +7", "Class B (ENSO
#             lag 0)" / "B 0". The lags come from the group_labels artefact
#             (09 step, read from the deployed fits), never from a constant.
#   centroid  the group's geographic centroid name from the sidecar's
#             display-only fields (centroid_label / centroid_short, display
#             order centroid_order): "West [GW]" / short "GW". No lags
#             needed.
# Colours stay keyed by id in every mode.
# ---------------------------------------------------------------------------

GROUP_LABEL_MODES <- c("registry", "lag", "centroid")

#' One group's lag as display text: "lead 6" (negative lag), "lag +7", "lag 0".
group_lag_text <- function(lag) {
  lag <- as.integer(lag)
  ifelse(
    lag < 0L,
    sprintf("lead %d", -lag),
    ifelse(lag > 0L, sprintf("lag +%d", lag), "lag 0")
  )
}

#' Short signed lag with a true minus: "−6", "+7", "0".
group_lag_signed <- function(lag) {
  lag <- as.integer(lag)
  ifelse(
    lag < 0L,
    sprintf("−%d", -lag),
    ifelse(lag > 0L, sprintf("+%d", lag), "0")
  )
}

#' group_display_meta
#'
#' Display labels of a fit registry's groups under results.style.group_labels.
#' @param meta load_grouping_meta() list of the registry (ids in meta$order)
#' @param lags named integer (group id -> deployed occurrence lag), or the
#'   group_labels artefact's `groups` table (enso_group, lag_occ); ignored
#'   under "registry" and "centroid", required (for every group of meta$order) under "lag"
#' @param mode "registry" | "lag" | "centroid"
#' @return list(mode, ids (sidecar order), order (display order: sidecar order
#'   under registry, letter order under lag, the sidecar's centroid_order under
#'   centroid), letter (named by id; NA under registry and centroid), lag (named
#'   by id or NULL), label, short, short_fig, colour (each
#'   named by id), map (stored string -> display string: ids -> label, stored
#'   short and short_fig -> display short))
group_display_meta <- function(meta, lags = NULL, mode = "registry") {
  mode <- as.character(mode %||% "registry")
  if (length(mode) != 1L || !mode %in% GROUP_LABEL_MODES) {
    stop(
      "group_display_meta: results.style.group_labels must be one of ",
      paste(GROUP_LABEL_MODES, collapse = " | "),
      "; got '",
      paste(mode, collapse = ","),
      "'.",
      call. = FALSE
    )
  }
  ids <- as.character(meta$order)
  colour <- if (!is.null(meta$colour)) {
    stats::setNames(as.character(meta$colour[ids]), ids)
  } else {
    NULL
  }
  if (identical(mode, "registry")) {
    out <- list(
      mode = mode,
      ids = ids,
      order = ids,
      letter = stats::setNames(rep(NA_character_, length(ids)), ids),
      lag = NULL,
      label = stats::setNames(as.character(meta$label[ids]), ids),
      short = stats::setNames(as.character(meta$short[ids]), ids),
      short_fig = stats::setNames(
        as.character((meta$short_fig %||% meta$short)[ids]),
        ids
      ),
      colour = colour
    )
  } else if (identical(mode, "centroid")) {
    cs <- meta[["centroid_short"]]
    cl <- meta[["centroid_label"]]
    cs <- if (is.null(cs)) {
      rep(NA_character_, length(ids))
    } else {
      as.character(cs[ids])
    }
    cl <- if (is.null(cl)) {
      rep(NA_character_, length(ids))
    } else {
      as.character(cl[ids])
    }
    bad <- ids[is.na(cs) | !nzchar(cs) | is.na(cl) | !nzchar(cl)]
    if (length(bad)) {
      stop(
        "group_display_meta: group_labels 'centroid' needs centroid_short and centroid_label ",
        "in the grouping sidecar (",
        meta$source %||% "?",
        ") for group(s) ",
        paste(bad, collapse = ", "),
        ". Add them to the sidecar, or to ",
        "kl_grouping.R::KL_CENTROID_NAMES for a derived registry (gen_grouping), ",
        "or use group_labels 'registry'.",
        call. = FALSE
      )
    }
    if (anyDuplicated(cs)) {
      stop(
        "group_display_meta: duplicated centroid_short in the grouping sidecar.",
        call. = FALSE
      )
    }
    ord <- as.character(meta[["centroid_order"]] %||% character(0))
    if (!length(ord)) {
      ord <- ids
    }
    if (!setequal(ord, ids)) {
      stop(
        "group_display_meta: centroid_order must list every group once.",
        call. = FALSE
      )
    }
    short <- stats::setNames(cs, ids)
    # Named vectors in DISPLAY order (legends / facets list GW, GG, ...).
    lab <- stats::setNames(sprintf("%s [%s]", cl, cs), ids)
    out <- list(
      mode = mode,
      ids = ids,
      order = ord,
      letter = stats::setNames(rep(NA_character_, length(ids)), ids)[ord],
      lag = NULL,
      label = lab[ord],
      short = short[ord],
      short_fig = short[ord],
      colour = colour
    )
  } else {
    if (is.data.frame(lags)) {
      lags <- stats::setNames(
        as.integer(lags$lag_occ),
        as.character(lags$enso_group)
      )
    }
    if (is.null(lags) || !length(lags)) {
      stop(
        "group_display_meta: group_labels 'lag' needs each group's deployed occurrence lag ",
        "(the group_labels artefact; 09 step group_labels).",
        call. = FALSE
      )
    }
    miss <- setdiff(ids, names(lags))
    if (length(miss)) {
      stop(
        "group_display_meta: no deployed occurrence lag for group(s) ",
        paste(miss, collapse = ", "),
        " (group_labels artefact of another registry?).",
        call. = FALSE
      )
    }
    lag <- stats::setNames(as.integer(lags[ids]), ids)
    if (anyNA(lag)) {
      stop(
        "group_display_meta: NA occurrence lag for ",
        paste(ids[is.na(lag)], collapse = ", "),
        call. = FALSE
      )
    }
    if (length(ids) > 26L) {
      stop("group_display_meta: more than 26 groups.", call. = FALSE)
    }
    # Longest lead first = ascending lag; ties keep the sidecar order.
    ord <- ids[order(lag, seq_along(ids))]
    letter <- stats::setNames(LETTERS[match(ids, ord)], ids)
    short <- stats::setNames(
      sprintf("%s %s", letter, group_lag_signed(lag)),
      ids
    )
    # Named vectors in DISPLAY order, so a legend or a facet that takes the
    # vector's names as its level order lists A, B, C, ...
    lab <- stats::setNames(
      sprintf("Class %s (ENSO %s)", letter, group_lag_text(lag)),
      ids
    )
    out <- list(
      mode = mode,
      ids = ids,
      order = ord,
      letter = letter[ord],
      lag = lag[ord],
      label = lab[ord],
      short = short[ord],
      short_fig = short[ord],
      colour = colour
    )
  }
  # Stored -> display: ids -> label; the sidecar's stored short / short_fig
  # (the strings 09 writes into region / scope columns) -> display short.
  src <- c(
    ids,
    as.character(meta$short[ids]),
    as.character((meta$short_fig %||% meta$short)[ids])
  )
  dst <- c(out$label[ids], out$short[ids], out$short_fig[ids])
  keep <- !duplicated(src)
  clash <- unique(src[duplicated(src) & !duplicated(paste(src, dst))])
  if (length(clash)) {
    stop(
      "group_display_meta: stored group string(s) ",
      paste(clash, collapse = ", "),
      " map to more than one display label (an id or short name shared by two groups).",
      call. = FALSE
    )
  }
  out$map <- stats::setNames(unname(dst[keep]), src[keep])
  out
}

#' relabel_groups
#'
#' Replace stored group strings (ids, stored short / short_fig names) by their
#' display labels (group_display_meta()$map) in the given columns -- whole-cell
#' matches only; factor levels are mapped in place (order kept). Identity
#' under mode "registry", so a default config prints the sidecar names.
#' @param df data frame (or NULL)
#' @param disp group_display_meta() result
#' @param cols columns to map; NULL = every character / factor column
#' @param rename_cols also map the column NAMES (tables pivoted wide by group)
relabel_groups <- function(df, disp, cols = NULL, rename_cols = TRUE) {
  if (is.null(df) || is.null(disp) || identical(disp$mode, "registry")) {
    return(df)
  }
  map <- disp$map
  sub1 <- function(x) {
    hit <- x %in% names(map)
    x[hit] <- unname(map[x[hit]])
    x
  }
  if (!is.data.frame(df)) {
    return(if (is.character(df)) sub1(df) else df)
  }
  cols <- cols %||%
    names(df)[vapply(
      df,
      function(x) is.character(x) || is.factor(x),
      logical(1)
    )]
  for (cc in intersect(cols, names(df))) {
    x <- df[[cc]]
    if (is.factor(x)) {
      lv <- sub1(levels(x))
      if (anyDuplicated(lv)) {
        stop(
          "relabel_groups: column ",
          cc,
          " maps two levels onto one label.",
          call. = FALSE
        )
      }
      levels(x) <- lv
    } else if (is.character(x)) {
      x <- sub1(x)
    }
    df[[cc]] <- x
  }
  if (isTRUE(rename_cols)) {
    names(df) <- sub1(names(df))
  }
  df
}

# COMPOSITE_CONDS (the ADM1 IRR map ONI conditions, Figure 2) lives in
# utils/figures_cb.R, sourced above: the grouped figure scripts need the same
# four evaluation points and are not part of the results phase.

# ---------------------------------------------------------------------------
# WHO-AFRO subregion of a country (11's with_afro()).
# A fixed copy of taxdat::get_AFRO_region() (the pipeline's AFRO_region
# factor), so 10 / 11 do not depend on taxdat: taxdat lists DJI, SOM and SDN
# under both Eastern and Western Africa and its case_when() returns the first
# match, Eastern Africa, kept here; its "GAN" typo is dropped. GMB, CPV, STP and
# COM (absent from taxdat) are added by geography. DJI, SOM and SDN are WHO
# EMRO member states, not AFRO: WHO_EMRO_ISO3 flags them for the footnote.
# Stored by report_dataset (supplementary tables S2 / S3 by subregion);
# elsewhere display only. AFRO_SUBREGION_ORDER: west to south.
# ---------------------------------------------------------------------------
AFRO_SUBREGION <- c(
  stats::setNames(
    rep("Eastern Africa", 13L),
    c(
      "BDI",
      "ETH",
      "KEN",
      "MDG",
      "RWA",
      "SDN",
      "SSD",
      "UGA",
      "TZA",
      "ERI",
      "DJI",
      "SOM",
      "COM"
    )
  ),
  stats::setNames(
    rep("Southern Africa", 9L),
    c("BWA", "MOZ", "MWI", "NAM", "SWZ", "ZMB", "ZWE", "ZAF", "LSO")
  ),
  stats::setNames(
    rep("Central Africa", 9L),
    c("AGO", "CMR", "CAF", "TCD", "COG", "COD", "GNQ", "GAB", "STP")
  ),
  stats::setNames(
    rep("Western Africa", 16L),
    c(
      "BEN",
      "BFA",
      "CIV",
      "GHA",
      "GIN",
      "GNB",
      "LBR",
      "MLI",
      "MRT",
      "NER",
      "NGA",
      "SEN",
      "SLE",
      "TGO",
      "GMB",
      "CPV"
    )
  )
)
WHO_EMRO_ISO3 <- c("DJI", "SOM", "SDN")
AFRO_SUBREGION_ORDER <- c(
  "Western Africa",
  "Central Africa",
  "Eastern Africa",
  "Southern Africa"
)

#' afro_subregion_of
#'
#' Subregion and WHO region of ISO3 codes (anything else, e.g. "Africa" or a
#' group name in a scope column, gets NA).
#' @return data.frame(afro_subregion, who_region) with one row per element
afro_subregion_of <- function(iso3) {
  x <- toupper(trimws(as.character(iso3)))
  sub <- unname(AFRO_SUBREGION[x])
  # The synthetic world's fake countries (gid.R::synthetic_countries(), call-time;
  # a no-op for real codes, which never match its XA? codes).
  if (anyNA(sub) && exists("synthetic_countries", mode = "function")) {
    syn <- synthetic_countries()
    hit <- is.na(sub) & x %in% syn$iso3
    sub[hit] <- syn$afro_region[match(x[hit], syn$iso3)]
  }
  data.frame(
    afro_subregion = sub,
    who_region = ifelse(
      is.na(sub),
      NA_character_,
      ifelse(x %in% WHO_EMRO_ISO3, "EMRO", "AFRO")
    ),
    stringsAsFactors = FALSE
  )
}

METRIC_LABELS <- c(
  cases = "Suspected cases per month",
  onsets = "Outbreak onsets per month"
)
METRIC_CUM <- c(
  cases = "cumulative suspected cases",
  onsets = "cumulative outbreak onsets"
)

# ---------------------------------------------------------------------------
# Supplementary parameter-plot registries (derive_track_params in 09;
# build_fig_supp_track_params in 10)
# ---------------------------------------------------------------------------

#' Per-track set of parameter families the supplementary caterpillar plots show.
#' ENSO is deliberately absent (it is Figures 1/2). 09 derives every family
#' comprehensively; 10 filters to the track's entry here. "other" carries the
#' terms param_family() does not name: the deployed neighbour-pressure term
#' (nbr_lag1_*), the population elasticity (log_pop_c) and the monthly
#' negative-binomial shape.
SUPP_PARAM_FAMILIES <- list(
  occurrence = c(
    "intercept",
    "climate",
    "spillover",
    "ttlo",
    "contagion",
    "other"
  ),
  duration = c(
    "intercept",
    "climate",
    "spillover",
    "seasonal",
    "month",
    "other"
  ),
  monthly = c("intercept", "climate", "spillover", "seasonal", "month", "other")
)

#' Classify a model term into a parameter family for the supplementary plots.
#' A thin wrapper over get_covar_type() (utils/plots.R, sourced by 09 alongside
#' results_core) that adds the families get_covar_type() lumps as "other": the
#' baseline intercept, the multi-scale spillover lags (ongoing_outbreak*), the
#' distance-weighted neighbour contagion (nbr_contagion*), and the seasonal
#' baseline-hazard vs. month-of-outbreak split. Returns NA for ENSO terms so
#' the caller can drop them (they are Figures 1/2).
param_family <- function(term) {
  if (length(term) != 1L) {
    return(vapply(term, param_family, character(1)))
  }
  if (is.na(term)) {
    return(NA_character_)
  }
  t <- as.character(term)
  if (grepl("enso", t, ignore.case = TRUE)) {
    return(NA_character_)
  } # ENSO -> Figures 1/2
  if (grepl("^Intercept", t)) {
    return("intercept")
  }
  if (grepl("^ongoing_outbreak", t)) {
    return("spillover")
  }
  if (grepl("nbr_contagion", t)) {
    return("contagion")
  }
  if (grepl("time_to_last_outbreak", t)) {
    return("ttlo")
  }
  if (grepl("seas", t)) {
    return("seasonal")
  } # t_start_seas / seasms
  # mso = months since onset (duration "month-of-outbreak" baseline hazard);
  # mso_band / msoms, plus calendar month-of-year (sin_month, factor(month)).
  if (grepl("mso|month", t, ignore.case = TRUE)) {
    return("month")
  }
  switch(
    get_covar_type(t),
    "climate_var" = "climate",
    "climate anomaly" = "climate",
    "climate mean" = "climate",
    "months" = "month",
    "years to last outbreak" = "ttlo",
    "enso teleconnection" = NA_character_,
    "other"
  )
}

# ---------------------------------------------------------------------------
# Artefact schema tokens
# ---------------------------------------------------------------------------

#' Per-artefact schema tokens, baked into the filename by make_results_run_file().
#'
#' Deliberately a CODE constant and not a config key: a token records a change
#' in the derived table's *columns or estimand* -- i.e. a change in the code --
#' not a user decision. In a user-editable YAML someone could silently un-bump
#' it, which is the one thing it exists to prevent.
#'
#' A token is bumped on a schema or estimand change that neither the mtime
#' freshness check (the raw inputs can be older than a stale table) nor the
#' provenance stamp can see; the stamp catches every *decision* change on its
#' own. The `cases-mon` token in the attribution filename is part of the key:
#' the cases metric is the mon track's cases_month.
# enso_relation: the `erf` table with 50% CrI columns (lo50/hi50/IRR_lo50/
# IRR_hi50, erf_from_draws()) for Figure 1d's inner ribbon, `adm1_draws`
# (per-draw ADM1 log-IRR array, extract_adm1_logrr_draws()), `enso_conds` and
# `n_draws`; the ONI conditions come from results.derivation.enso_conds.
# attribution: load_sim_run_monthly() keys countries by gid ISO3. Ribbons carry
# 2.5% / 97.5% quantiles (*_q025 / *_q975) for uniform 50% + 95% bands;
# `cum_totals` (paired replicate totals of cases and onsets over the analogue
# window, Figure 3c); per-arm launches (results.derivation.arm_alignment):
# factual_start / arm_starts / stitch_months / ratio_label; `country_ribbon`
# (forecast / counterfactual monthly ribbon per country); `cells`
# (attrib_cells(): window x explosive mode) and `africa_window` -- with a
# common window the cells attribution_common / post_fork_common /
# attribution_firm sit beside the others -- plus cum_totals_cells /
# ratio_stats_cells, the ribbons' `continuation` flag and africa_scope
# (table S28).
# smoothing: quantile ribbons (q025 / q975); `by_country` (masked / unmasked
# ribbons + observed per country and metric) and `country_group`;
# `validation` (smoothing_validation_tables(): randomised PIT histograms of
# the observed monthly onsets and cases, onset frequency by calendar month,
# episode durations, fraction of covered unit-months active, onset frequency
# given an ongoing outbreak in the previous month in the same ADM1 / another
# ADM1 of the country / another country of the TU, and season_onsets /
# max_season_onsets / n_unit_seasons, the structural check of
# simulator.onset_eligibility: season; observed vs across replicates). Inputs
# include the coverage mask's aggregated series and the HMM cluster file.
# oni: ONI envelope with lo95 / hi95, `obs_full` (the whole observed ONI
# record), `phases` (NOAA ENSO-episode runs on it, extended by the CFSv2
# ensemble mean) for the ONI strips and phase shading of Figures 1 and 3,
# and `handoffs` (per arm).
RESULTS_SCHEMA <- c(
  enso_relation = "v3",
  attribution = "v10",
  smoothing = "v5",
  oni = "v4",
  # Per-gid outbreak-onset rate over surveillance-covered
  # months, for Figure 1b's inset map (derive_onset_rate in 09).
  onset_rate = "v1",
  # Region/country parameter summaries for the supplementary
  # caterpillar plots (derive_track_params in 09; plotted by
  # build_fig_supp_track_params in 10). Per-track tidy tibble of
  # population + country parameter estimates (every term except
  # ENSO, which is Figures 1/2).
  track_params = "v1",
  # Manuscript "key numbers" report artefacts (phase 11 notebook
  # reads these; 09 derives them). Each is a bundle of tidy,
  # kable-ready tibbles keyed as noted in make_results_run_file().
  # report_dataset: coverage on OBSERVED sub-ADM1 populations
  # (build_adm1_obs_pop, recover_dropped=TRUE), episodes_by_country
  # (ENSO group x country x admin level, table S2), cascade
  # (processing counts, table S1), raw_counts / raw_table1
  # (01's outbreak set, ADM0 excluded, Zheng layout incl. ADM3),
  # outbreaks (one row per analysis outbreak with its season-year
  # DJF ENSO state, gid and afro_subregion, for 11's gtsummary
  # tables), table1_region (the analysis set by admin level incl.
  # All, for Africa and each WHO-AFRO subregion; table S2) and
  # enso_region / enso_region_p (outbreaks by season-year ENSO
  # state per block, outbreaks per cholera year of the state,
  # Kruskal-Wallis p; table S3).
  # report_effects: effects population on the sub-ADM1 basis;
  # effects_region and effects_country keep range_lo / range_hi
  # (the group crossbasis's fitted ONI range, from
  # enso_relation), so 11 daggers conditions outside it.
  report_dataset = "v5",
  report_effects = "v6",
  # report_forecast: horizon_totals with the annualised
  # incidence inc_ann_100k and enso_group/region columns,
  # scope_pop (simulated-gid population denominators), n_months
  # (simulated vs calendar months per window), factual_start /
  # arm_starts; the region peak ratio is the max of the
  # region's monthly total; Africa rows follow the Africa
  # window and pop_grid / pop_shift carry africa_def = "mixed".
  report_smoothing = "v1",
  report_forecast = "v5",
  # ENSO exposure by direction / magnitude / certainty
  # (derive_report_enso_exposure), incl. direction `either`
  # and the grid track `any_credible` (the track-wise union:
  # some single track has P >= p).
  report_enso_exposure = "v3",
  # Population in ADM1s with more ENSO-attributable cases,
  # forecast vs counterfactual (derive_report_attrib_exposure):
  # Childs ratio-of-means point + paired-bootstrap interval,
  # childs_median (ratio of medians) and IPCC robustness
  # classes, over `window` (attribution | post_fork) and
  # `explosive_guard` axes; sets-keyed. A zero counterfactual
  # statistic with a positive forecast counts (ratio = Inf);
  # only 0 vs 0 is undefined.
  report_attrib_exposure = "v4",
  # Attributable number / fraction, probability ratio (+ FAR,
  # PS, PNS) for onset_any / cases_ge, and signal-to-noise of
  # the paired runs of report_attrib_exposure, same axes, plus
  # ipcc_scope, the IPCC class of each scope's (Africa,
  # region, country) summed paired totals at every
  # thresholds.attrib_magnitudes (Figure 4a's ratio bars).
  # Sets-keyed.
  report_attrib_metrics = "v3",
  # Attribution of past strong ENSO events: the
  # report_attrib_metrics families + report_attrib_exposure
  # classes per event x factual source (simulation |
  # smoothing) x explosive guard, observed totals and the
  # observed-anchored AN. Unkeyed (the analog set is the event
  # manifest's). + explosive_summary (table S29: what the
  # explosive guard drops, per event x factual source x set x
  # scope; explosive_cap from the stamped
  # results.thresholds.explosive_case_cap) and sub_windows (one
  # cell per per-winter sub-window of a schema-3 event
  # manifest, enso_events_v2.yml). The probability-ratio tables
  # (pr_adm1 / pr_scope / pr_pop and their sub-window versions)
  # keep their outbreak event as `pr_event`.
  report_event_attrib = "v4",
  # Its analytic cross-check: crossbasis-only AF / AN of the
  # same events on the observed panel, with sub_windows.
  report_event_af_model = "v2",
  # Posterior interval of the attributable fraction / number
  # from fixed-draw runs: per-draw AF / AN with a
  # delta-method noise adjustment. af_posterior (forecast,
  # sets-keyed; cells attribution / post_fork) and
  # event_af_posterior (events, unkeyed; explosive modes / cap
  # from the results config, + sub_windows). Africa per draw
  # is the replicate-level sum of every unit (no AF_MC
  # combinations, n_mc NA), the explosive guard drops a
  # replicate from its whole unit (derivation.af_guard_scope),
  # and a scope / ADM1 whose noise share reaches
  # thresholds.af_noise_share_max is flagged unidentified
  # (adjusted interval NA; Childs population: variants
  # unidentified / unidentified_pass).
  af_posterior = "v3",
  event_af_posterior = "v4",
  # Posterior interval of the Childs population: the Childs
  # rule per draw on the noise-adjusted ADM1 AF (variant adj)
  # and on the raw per-draw means (raw), from the same
  # fixed-draw runs. Sets-keyed.
  attrib_exposure_posterior = "v3",
  # Information component of the forecast attribution:
  # displayed (stitched) factual vs attribution factual vs
  # counterfactual, per draw; AN_total = AN_information
  # + AN_enso. Sets-keyed.
  report_attrib_information = "v2",
  # Window totals and annualised incidence of the paired
  # forecast / counterfactual runs: Africa, groups and
  # countries per attribution cell. Sets-keyed.
  report_attrib_incidence = "v1",
  # Hindcast forecast skill: derived by 09 from the
  # hindcast manifest + each run's validation/metrics and
  # merged/monthly_summary through utils/hindcast_skill.R; read
  # by 10 (figure H), notebooks/hindcast_skill.qmd and 11.
  # A reference score that sums to zero up to round-off leaves
  # the CRPSS undefined (NA).
  hindcast_runs = "v1",
  hindcast_occ = "v2",
  hindcast_cases = "v2",
  report_hindcast = "v1",
  # Supplementary tables: MCMC convergence of the
  # deployed fits (convergence, track-keyed; report_convergence,
  # table S8 / fig. S5, with each fit's deploy-gate verdict --
  # gate_verdict pass | allowed | fail, gate_reasons, ... -- and
  # the deploy's convergence block `gate`), selected lags + firm
  # horizons and the LOO lag rankings (report_lags, S9 / S10 /
  # fig. S6), the baseline-rate anchoring tables (report_anchor,
  # S7), software versions (report_software, S13), and the
  # brms_kl class-sweep comparison (report_kl_sweep).
  convergence = "v1",
  report_convergence = "v2",
  # report_lags: selected floor_lag, bound_active / bound_cost
  # (+ SE: the paired elpd cost of the deploy floor when it
  # binds), rank2_lag and best_nonneg_lag with their elpd
  # difference to rank 1, the deployed lag's Pareto-k counts;
  # the ranking carries Pareto-k columns (NA when the 05
  # ranking predates them); loo_summary (table T2,
  # report_lags_loo_summary(), at rule$se_mult).
  report_lags = "v3",
  # Fit-level model specification of the deployed fits:
  # formula terms, RE groups,
  # family, coefficients and priors per fit, checked against
  # analysis/model_terms_deployed.yml (table S34). Track-keyed;
  # written by the convergence step from its fit load.
  report_model_spec = "v1",
  report_anchor = "v1",
  report_software = "v1",
  report_kl_sweep = "v1",
  # Figure 1d: rainfall anomaly-ONI correlation per
  # ENSO group x centre month x lag, with within-group agreement,
  # and the groups' main rainy seasons; from 03_prepare_rain_adm1.R's
  # ADM1 series. + group_lags (table T1b's rainfall columns per
  # weighting x group, rain_enso_lag_summary()) and mechanisms
  # (table T1a, a copy of the curated
  # results.inputs.enso_rain_mechanisms; NULL when "").
  rain_enso = "v2",
  # Group display lags + per-country season ends: each group's
  # deployed occ / dur / mon lag
  # read from the fits (style.group_labels "lag" builds the
  # class letters from it in 10 / 11), members, counts, and
  # table S30 (season_ends, season_windows). Ids only.
  group_labels = "v1",
  # Lag averaging: per track and
  # group, Pareto-k and k-robust elpd differences per lag,
  # stacking / pseudo-BMA+ weights on the `all` and `eligible`
  # lag sets, per-lag and mixture population ERFs on a common
  # grid and the IRR at enso_conds. Track-keyed.
  lag_averaging = "v1",
  # Membership uncertainty of the deployed classes: per
  # country p_max and the
  # responsibility row of the deployed sweep cell, group
  # summary, partition stability, sensitivity to the sweep
  # cell and chain set (descriptive elpd), and the
  # membership-weighted effects M1 (kl model) and M2 (deployed
  # ERFs at each group's lag, from the 05 lag sidecars).
  report_membership = "v1",
  # Per-country WHO term scale of the smoothing (table D30):
  # the estimated table 08a
  # wrote, or (legacy preview, a fixed-scale smoothing) the
  # implied WHO/DB ratios recomputed from the run's frame;
  # + national_monthly, the WHO national monthly series
  # (03_prepare_who_monthly.R, table D36).
  report_who_scale = "v2",
  # Counterfactual analog sets:
  # per arm (forecast / event) x set, size, model shares,
  # seam error, per-winter means and the min_set_size flag
  # (table D29), from the Lenssen caches the runs read;
  # + run_design (table D34:
  # replicates, draws, draw mode, handoff, analog pool and
  # clamp margin per simulation arm, from the run configs).
  report_analog_sets = "v2",
  # Occurrence-panel tables: the national-zero audit (D32),
  # each on only when the panel's sidecar exists.
  report_national_zero_audit = "v1",
  # ... and the occurrence-panel exclusions per country (D31).
  report_occ_panel_exclusions = "v1",
  # ... and the contagion term's surveillance coverage (D33).
  report_contagion_coverage = "v1",
  # Supplement maps MU / GM / KC: modelled units by level and
  # status, the ADM1s
  # of the modelled countries without a recorded episode, and
  # the fit registry's membership (derive_modelled_units).
  modelled_units = "v1",
  # Figure 1c's onset-season column: the deployed occurrence
  # fits' seasonal
  # spline per ENSO group and calendar month, relative to the
  # annual mean (derive_occ_seasonality), so 10 only reads;
  # `draws`, the per-draw calendar curves, for seasonality's band.
  occ_seasonality = "v2",
  # Supplement fig. SE: per
  # ENSO group and calendar month, the rainfall climatology
  # and rainy season (rain_enso) against the observed onset,
  # cessation and case rates, with the occurrence and monthly
  # models' fitted seasonality and each track's peak month
  # (derive_seasonality), so 10 / 11 only read; the
  # occurrence band (lo / hi at `interval`, from
  # occ_seasonality's draws).
  seasonality = "v2",
  # Figure 1c's top row:
  # the composite El Nino / La Nina ONI cycle over Y, Y+1
  # (derive_oni_composite), so 10 only reads.
  oni_composite = "v1",
  # Supplement fig. WR: observed WHO ratio of
  # the target year to the ENSO-neutral years against the
  # smoothing / ENSO-counterfactual ratio of means over their
  # common months (derive_report_who_enso_ratio).
  report_who_enso_ratio = "v1",
  # Supplementary table: ENSO contrasts of the
  # three tracks, all outbreaks vs the subset with a confirmed
  # case, from 01's outbreak table (confirmed_agreement_table).
  report_confirmed_agreement = "v2"
)

#' Cache key for a set list. An aggregate derived for one list of analog sets
#' must never be silently reused for another, and the key must not depend on the
#' order the sets were requested in or the same table lands under two names.
sets_key <- function(sets) paste(sort(sets), collapse = "-")

# ---------------------------------------------------------------------------
# Config access
# ---------------------------------------------------------------------------

#' Recover the -C/--config path the current script was invoked with.
#'
#' Needed for two things: deriving the artefact run tag (so artefacts land under
#' generated_data/results/<tag>/), and printing the exact re-run command when 10
#' finds a missing or stale artefact. Same idiom as 07a_simulator_tasks.R:47.
#' Select the reported series (results.series) from a simulator/smoothing table.
#'   "latent"   -- the table as-is (pre-existing behaviour).
#'   "detected" -- what the surveillance record would have shown: rows/values of
#'                 episodes flagged by the per-episode detection draw
#'                 (covariate_feed: detected, column `detected`). For event
#'                 tables (one row per onset or per active month) the rows are
#'                 filtered; for the v1 monthly summary, Y and cases_month are
#'                 zeroed where undetected (the grid itself must stay complete).
#' A "detected" request against a run without the column is an error naming the
#' run, never a silent fall-back to latent.
results_select_series <- function(
  tbl,
  series,
  what,
  kind = c("events", "grid")
) {
  kind <- match.arg(kind)
  if (identical(series %||% "latent", "latent") || is.null(tbl)) {
    return(tbl)
  }
  if (!identical(series, "detected")) {
    stop(
      "results.series must be 'latent' or 'detected', got '",
      series,
      "'.",
      call. = FALSE
    )
  }
  if (!"detected" %in% names(tbl)) {
    stop(
      "results.series = 'detected' but ",
      what,
      " carries no `detected` column ",
      "-- it was not simulated with simulator.filter.covariate_feed: detected.",
      call. = FALSE
    )
  }
  if (kind == "events") {
    return(tbl[tbl$detected %in% 1L, , drop = FALSE])
  }
  keep <- tbl$detected %in% 1L
  if ("Y" %in% names(tbl)) {
    tbl$Y[!keep] <- 0L
  }
  if ("cases_month" %in% names(tbl)) {
    tbl$cases_month[!keep] <- NA
  }
  tbl
}

results_config_path <- function(argv = commandArgs(TRUE)) {
  p <- .flag_val("-C", argv) %||% .flag_val("--config", argv)
  if (is.null(p) || !nzchar(p)) {
    stop(
      "[results] no -C/--config given: the manuscript phase is config-driven ",
      "and the config path also names the artefact directory. Run with ",
      "-C analysis/configs/results/manuscript.yml",
      call. = FALSE
    )
  }
  p
}

#' Resolve the analog sets in play, honouring scope$supplementary.
#' Ordered by ENSO_ANALOG_SETS so the sets key is canonical.
resolve_results_sets <- function(res) {
  main <- as.character(res$scope$main_set)
  supp <- if (isTRUE(res$scope$supplementary)) {
    setdiff(trimws(as.character(res$scope$supp_sets %||% character(0))), "")
  } else {
    character(0)
  }
  sets <- ENSO_ANALOG_SETS[ENSO_ANALOG_SETS %in% unique(c(main, supp))]
  if (!length(sets)) {
    stop(
      "[results] no valid analog sets requested (known: ",
      paste(ENSO_ANALOG_SETS, collapse = ","),
      ")",
      call. = FALSE
    )
  }
  sets
}

#' Existence guard with a phase-tagged message.
require_input <- function(path, what, phase = "results") {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(
      sprintf(
        "[%s] missing required input (%s): %s",
        phase,
        what,
        path %||% "<NULL>"
      ),
      call. = FALSE
    )
  }
  path
}

# ---------------------------------------------------------------------------
# Manifest access
# ---------------------------------------------------------------------------
# Every accessor takes its input as an explicit argument. Both scripts resolve `MAN <-
# load_manifest(...)` once.

load_manifest <- function(path) {
  require_input(path, "forecast manifest")
  man <- yaml::read_yaml(path)
  if (is.null(man$groups) || !length(man$groups)) {
    stop("[results] manifest has no groups: ", path, call. = FALSE)
  }
  man
}

#' The factual configs of a forecast manifest, sorted (one per unit / group).
.manifest_factual_configs <- function(man) {
  sort(unname(vapply(
    manifest_runs(man, "factual"),
    function(r) as.character(r$config %||% ""),
    character(1)
  )))
}

#' TRUE when the attribution manifest names other factual runs than the forecast
#' manifest (the attribution factual launched at the counterfactual handoff),
#' FALSE when both name the same runs (launches already aligned, or no
#' attribution manifest).
results_attrib_factual_differs <- function(man, man_attrib) {
  !identical(
    .manifest_factual_configs(man),
    .manifest_factual_configs(man_attrib)
  )
}

#' Refuse an attribution manifest (utils/make_forecast_attrib_configs.R) that was
#' not generated from the forecast manifest the results read, or whose units and
#' counterfactual runs differ from it: the attribution pairs its factual with the
#' SAME counterfactual runs the display reads.
results_attrib_manifest_check <- function(
  man,
  man_attrib,
  manifest_path,
  attrib_path
) {
  norm <- function(p) normalizePath(p, mustWork = FALSE)
  src <- as.character(man_attrib$source_manifest %||% "")
  if (!nzchar(src) || norm(src) != norm(manifest_path)) {
    stop(
      sprintf(
        paste0(
          "[results] the attribution manifest %s was generated from '%s', not from ",
          "results.inputs.manifest %s. Re-run deploy_prod_run.R --stage gen_attrib."
        ),
        attrib_path,
        src,
        manifest_path
      ),
      call. = FALSE
    )
  }
  ua <- manifest_units(man)
  ub <- manifest_units(man_attrib)
  if (!setequal(names(ua), names(ub))) {
    stop(
      sprintf(
        "[results] the attribution manifest %s has units %s; the forecast manifest has %s.",
        attrib_path,
        paste(sort(names(ub)), collapse = ","),
        paste(sort(names(ua)), collapse = ",")
      ),
      call. = FALSE
    )
  }
  for (aset in as.character(unlist(
    man$analog_sets %||% man$primary_analog_set
  ))) {
    a <- vapply(
      manifest_runs(man, "counterfactual", aset),
      function(r) r$config %||% "",
      character(1)
    )
    b <- vapply(
      manifest_runs(man_attrib, "counterfactual", aset),
      function(r) r$config %||% "",
      character(1)
    )
    if (!identical(a[sort(names(a))], b[sort(names(b))])) {
      stop(
        sprintf(
          "[results] the attribution manifest %s names other counterfactual runs (set %s) than %s.",
          attrib_path,
          aset,
          manifest_path
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

load_enso_grouping <- function(path) {
  require_input(path, "ENSO grouping registry")
  # Same parse/validation as every other registry read (utils/config.R::
  # load_grouping(): trimmed, blank rows dropped, duplicate iso3 refused),
  # returned as a tibble.
  tibble::as_tibble(load_grouping(
    sub("\\.csv$", "", basename(path)),
    groupings_dir = dirname(path)
  ))
}

# All groups the manifest declares (ENSO exposure-response: every group).
manifest_groups_all <- function(man) names(man$groups)

# Groups included in Figures 3-4: every manifest group but results.scope.exclude_groups
# (empty by default).
manifest_groups_inc <- function(man, exclude = character(0)) {
  setdiff(names(man$groups), exclude)
}

# build_location_frame() in the results phase (09, and the sibling
# 07d_simulator_figures_counterfactual.R) is called WITHOUT any rate table:
# the manuscript phase never evaluates a linear predictor -- it reads
# already-simulated partition outputs and uses `frame` purely as the
# country/AFRO_region and gid/ADM1 crosswalk, so the occ/dur/mon anchors
# (`offset(log(base_rate_<track>))`, utils/anchor.R) would be dead weight and
# `build_location_frame(base_rate_table_paths = NULL)` (the default) adds no
# column at all. If a figure here ever starts scoring mu, thread
# resolve_base_rate_paths() in exactly as 07b/07c do.

#' Resolve the run_tag for one group/arm/analog_set from the manifest: the
#' `manifest` primary-set
#' alias resolves through `counterfactual_config`; non-primary sets require a
#' `counterfactual_sets` block and error loudly rather than silently falling
#' back to the primary set under the wrong label.
# The manifest entry that carries a fit group's RUN fields: the group itself
# (v1), or -- v2, partition_by: registry -- the first partition unit holding
# it (a group's countries may sit in several units; loaders that need every
# run use manifest_units() / results_run_configs()).
.manifest_run_entry <- function(man, g) {
  if (!length(man$units)) {
    return(man$groups[[g]])
  }
  u <- names(man$units)[vapply(
    man$units,
    function(x) g %in% unlist(x$groups),
    logical(1)
  )]
  if (!length(u)) {
    stop(
      sprintf(
        "[results] fit group '%s' is in no partition unit of the manifest",
        g
      ),
      call. = FALSE
    )
  }
  man$units[[u[1L]]]
}

#' results_run_configs
#'
#' Every run config behind the included fit groups (forecast + each analog set):
#' per group in v1 (order unchanged), per partition unit holding any of them in
#' v2. Used for the freshness input lists.
results_run_configs <- function(man, groups, sets) {
  unlist(unname(results_unit_arm_configs(man, groups, sets)))
}

#' results_unit_arm_configs
#'
#' The run configs of every unit (registry manifests) or group (v1) holding one
#' of `groups`, per unit: c(<factual config>, <set> = <counterfactual config>
#' for each of `sets`) -- the factual one first and unnamed. Named list by unit.
results_unit_arm_configs <- function(man, groups, sets) {
  keys <- if (length(man$units)) {
    names(man$units)[vapply(
      man$units,
      function(x) any(groups %in% unlist(x$groups)),
      logical(1)
    )]
  } else {
    groups
  }
  ent <- if (length(man$units)) man$units else man$groups
  lapply(stats::setNames(keys, keys), function(k) {
    e <- ent[[k]]
    primary <- as.character(e$primary_analog_set %||% "all")
    c(
      e$factual_config,
      vapply(
        sets,
        function(a) {
          if (identical(a, primary)) {
            e$counterfactual_config %||% e$counterfactual_sets[[a]]$config
          } else {
            e$counterfactual_sets[[a]]$config
          }
        },
        character(1)
      )
    )
  })
}

group_run_tag <- function(man, g, arm, analog_set = NULL) {
  gi <- .manifest_run_entry(man, g)
  cfg_path <- if (arm == "forecast") {
    gi$factual_config
  } else {
    primary <- as.character(gi$primary_analog_set %||% "all")
    aset <- as.character(analog_set %||% primary)
    if (identical(aset, primary)) {
      gi$counterfactual_config
    } else if (is.null(gi$counterfactual_sets)) {
      stop(sprintf(
        paste(
          "[results] group '%s': analog_set '%s' requested but the manifest",
          "predates analog sets (primary is '%s')."
        ),
        g,
        aset,
        primary
      ))
    } else {
      cs <- gi$counterfactual_sets[[aset]]
      if (is.null(cs$config)) {
        stop(sprintf(
          "[results] group '%s' has no counterfactual config for analog_set '%s'",
          g,
          aset
        ))
      }
      cs$config
    }
  }
  if (is.null(cfg_path)) {
    stop(sprintf("[results] group '%s': no %s config in manifest", g, arm))
  }
  yaml::read_yaml(cfg_path)$simulator$run_tag
}

group_config_path <- function(man, g, arm, analog_set = NULL) {
  gi <- .manifest_run_entry(man, g)
  if (arm == "forecast") {
    return(gi$factual_config)
  }
  primary <- as.character(gi$primary_analog_set %||% "all")
  aset <- as.character(analog_set %||% primary)
  if (identical(aset, primary)) {
    gi$counterfactual_config
  } else {
    gi$counterfactual_sets[[aset]]$config
  }
}

#' Attribution window start (Figures 3-4): from each arm's launch
#' (utils/sim_io.R::manifest_attrib_start(): the manifest's launch keys, else
#' each config's simulator.horizon.start), uniform per arm across the included
#' groups. `alignment` is results.derivation.arm_alignment: "stitch" = the
#' earlier launch (the factual series stitched from the smoothing before its
#' own launch), "common_start" = the later one.
resolve_attrib_start <- function(
  man,
  groups,
  sets = NULL,
  alignment = "common_start"
) {
  manifest_attrib_start(man, groups, sets, alignment = alignment)
}

#' The factual arm's launch over the included groups.
resolve_factual_start <- function(man, groups, sets = NULL) {
  manifest_factual_start(man, groups, sets)
}

# ---------------------------------------------------------------------------
# Artefact freshness + provenance contract
# ---------------------------------------------------------------------------
# One predicate, two wrappers. An mtime-only rule would be blind to the
# decisions that produce the table (see attach_results_provenance()'s docstring, utils/provenance.R).
#
# The input list is declared by 09 at derivation time and CARRIED ON THE STAMP,
# rather than recomputed by each reader. That is what lets phase 10 verify
# freshness without ever touching a fit archive: enumerating the ENSO
# exposure-response inputs means resolving the best-lag manifest and the
# per-track LOO batch directory,
# and giving 10 those accessors would reopen exactly the path to a fit object
# that the split closes. 10 instead re-stats what 09 recorded.
#
# Globs are recorded alongside the resolved files for the one case re-stat'ing
# cannot see: a file ADDED to a directory that is read by pattern (the
# loo_compare_lags_*.rds batch). Re-globbing is pure filesystem work, so 10 can
# do it knowing nothing about what the files contain.

#' Resolve an input declaration to a flat character vector of paths.
#'
#' @param files character vector of explicit paths
#' @param globs list of list(dir=, pattern=, recursive=) specs, read with
#'   list.files(). `recursive` defaults to FALSE; the smoothing partition root
#'   needs TRUE (its _SUCCESS sentinels sit one level down, per partition).
results_inputs <- function(files = character(0), globs = list()) {
  files <- as.character(files)
  files <- files[!is.na(files) & nzchar(files)]
  globbed <- unlist(
    lapply(globs, function(g) {
      if (!dir.exists(g$dir)) {
        return(character(0))
      }
      list.files(
        g$dir,
        pattern = g$pattern,
        full.names = TRUE,
        recursive = isTRUE(g$recursive)
      )
    }),
    use.names = FALSE
  )
  # The glob specs ride along as an attribute so the stamp can record them
  # (results_write(globs =)) and 10 can detect inputs added since derivation.
  structure(sort(unique(c(files, globbed %||% character(0)))), globs = globs)
}

#' The mtime rule: an artefact is
#' stale if any input that exists is newer than it. A stale figure that silently
#' outlives a re-run of an upstream pipeline stage is worse than a slower
#' rebuild.
#'
#' @return NULL if fresh, else a one-line reason.
.results_stale_reason <- function(path, inputs = character(0)) {
  if (!file.exists(path)) {
    return(sprintf("artefact does not exist: %s", path))
  }
  inputs <- inputs[nzchar(inputs) & file.exists(inputs)]
  if (!length(inputs)) {
    return(NULL)
  }
  newer <- inputs[file.mtime(inputs) > file.mtime(path)]
  if (length(newer)) {
    return(sprintf("input newer than artefact: %s", newer[1]))
  }
  NULL
}

#' Reader-side freshness, using ONLY what the stamp recorded -- see the section
#' header for why 10 cannot re-enumerate the inputs itself.
#'
#' Catches: an input deleted, an input whose mtime or size changed since
#' derivation, and (via the recorded globs) an input newly ADDED to a
#' pattern-read directory.
#'
#' @return NULL if fresh, else a one-line reason.
.results_stamp_stale_reason <- function(stamp) {
  rec <- stamp$inputs
  if (is.null(rec) || !nrow(rec)) {
    return(NULL)
  }

  gone <- rec$path[!file.exists(rec$path)]
  if (length(gone)) {
    return(sprintf("input no longer exists: %s", gone[1]))
  }

  now_mtime <- file.mtime(rec$path)
  now_size <- file.size(rec$path)
  diff <- now_mtime != rec$mtime | now_size != rec$size
  # A simulator monthly summary rewritten in place in the compact format
  # (utils/compact_sim_outputs.R: identical content, mtime kept, smaller) is the
  # input the stamp recorded: its sidecar must name the recorded size and mtime
  # (utils/sim_io.R::sim_monthly_converted_from()).
  sz <- which(diff & now_mtime == rec$mtime)
  if (length(sz)) {
    diff[sz] <- !sim_monthly_converted_from(
      rec$path[sz],
      rec$mtime[sz],
      rec$size[sz]
    )
  }
  changed <- rec$path[diff]
  if (length(changed)) {
    return(sprintf("input changed since derivation: %s", changed[1]))
  }

  added <- setdiff(results_inputs(globs = stamp$globs %||% list()), rec$path)
  if (length(added)) {
    return(sprintf("input added since derivation: %s", added[1]))
  }

  NULL
}

#' The decisions stamped onto every artefact. `scope` is included because it
#' selects what is derived, not merely what is plotted -- but `style` is NOT:
#' changing a font or a DPI is presentation, and must not invalidate a table.
.RESULTS_STAMPED_BLOCKS <- c("inputs", "scope", "derivation", "thresholds")

#' Step-scoped decisions: stamped leaves (dot-path prefixes) that shape ONE
#' kind's artefact only. Every artefact still carries them in its stamp, but
#' .results_decision_diff() compares them only for that kind (the stamp's
#' `kind`), so switching an optional step on re-derives that step alone, not
#' every artefact. clim_compare is retired (the local-climate comparison is not
#' in the public pipeline): its leaves stay listed so stamps that carry them
#' are not stale for every other kind.
.RESULTS_STEP_DECISIONS <- list(
  clim_compare = c(
    "inputs.clim_compare",
    "derivation.clim_compare",
    "thresholds.clim_compare"
  )
)

.results_stamp <- function(
  res,
  config_path,
  schema,
  inputs = character(0),
  globs = list(),
  kind = NULL
) {
  inputs <- inputs[nzchar(inputs)]
  ex <- inputs[file.exists(inputs)]
  list(
    results = res[.RESULTS_STAMPED_BLOCKS],
    kind = kind,
    config_path = config_path,
    schema = schema,
    derived_at = Sys.time(),
    r_version = R.version.string,
    globs = globs,
    inputs = data.frame(
      path = ex,
      mtime = file.mtime(ex),
      size = file.size(ex),
      stringsAsFactors = FALSE
    )
  )
}

#' Every leaf path of a nested list, as dot-separated strings, paired with its
#' value -- so a stamp/config mismatch can name the exact differing key rather
#' than reporting "the blocks differ".
.results_leaves <- function(x, prefix = character(0)) {
  if (!is.list(x) || is.null(names(x)) || !length(x)) {
    out <- list(x)
    names(out) <- paste(prefix, collapse = ".")
    return(out)
  }
  do.call(
    c,
    lapply(names(x), function(nm) .results_leaves(x[[nm]], c(prefix, nm)))
  )
}

#' Which stamped decisions disagree with the currently-resolved config.
#' THIS is the check the mtime rule cannot do: it catches "you edited
#' prob_threshold but did not re-run phase 09".
# Stamped keys added after artefacts already existed, with the value that
# reproduces the earlier behaviour: an artefact stamped before the key existed
# counts as derived under that value (so adding a key does not invalidate every
# artefact on disk). Changing the key away from it still re-derives.
#
# derivation.arm_alignment is back-filled with "stitch", the default, not
# "common_start": the stamp covers the whole derivation block, so
# "common_start" would mark every artefact stale under the default although
# only attribution / oni / report_forecast depend on the key (and their
# RESULTS_SCHEMA tokens re-derive those). A config stating "common_start"
# re-derives as usual.
.RESULTS_STAMP_BACKCOMPAT <- list(
  scope.beyond_firm_horizon = "drop",
  # monthly / season case skill, bootstrap and
  # area summaries: all off.
  scope.hindcast.headline_case_skill_reference = "clim_pooled",
  derivation.hindcast.case_skill_supports = character(0),
  derivation.hindcast.case_month_grains = c("country", "adm1"),
  derivation.hindcast.case_skill_references = c(
    "clim_pooled",
    "clim_launch_state",
    "persistence_season",
    "persistence_last"
  ),
  derivation.hindcast.case_skill_rules = c("crps_raw", "crps_log1p", "twcrps"),
  derivation.hindcast.case_clim_halfwidth = 1L,
  derivation.hindcast.case_pop_scale = "none",
  derivation.hindcast.case_mask = "none",
  derivation.hindcast.case_tw_threshold_month = 50,
  derivation.hindcast.case_tw_threshold_season = 50,
  derivation.hindcast.case_season_min_ref = 5L,
  derivation.hindcast.bootstrap_n = 0L,
  derivation.hindcast.bootstrap_seed = 20261101L,
  derivation.hindcast.bootstrap_level = 0.95,
  derivation.hindcast.skill_area_summaries = FALSE,
  derivation.arm_alignment = "stitch",
  # national-zero threshold: the exact-zero rule.
  inputs.episode_national_zero_threshold = 0L,
  # national-zero rule: the legacy rule.
  inputs.episode_national_zero_rule = "",
  # at-risk start / ADM1 overlap: the default panel.
  inputs.episode_at_risk_start = "season",
  inputs.episode_adm1_overlap = "drop",
  # report_attrib_exposure thresholds: only that
  # step reads them.
  thresholds.attrib_magnitudes = c(0.25, 0.5, 1),
  thresholds.robust_exceed = 0.66,
  thresholds.robust_agree = 0.8,
  # report_attrib_metrics thresholds: only that
  # step reads them.
  thresholds.pr_case_threshold = 50,
  thresholds.pr_alpha = 0.05,
  thresholds.sn_levels = c(1, 2),
  # bootstrap resamples: B = 500.
  derivation.attrib_boot_b = 500L,
  # attribution factual: the forecast
  # manifest's (stitched) factual.
  inputs.attrib_manifest = "",
  # One Africa window: the legacy windows (no
  # common cell, firm horizons only) and no
  # continuation-skill input.
  scope.africa.window_end = "legacy",
  scope.africa.continuation_max_lead = 0L,
  scope.africa.max_counterfactual_winter = "",
  inputs.enso_continuation_skill = "",
  # T1a's curated table as a results input.
  inputs.enso_rain_mechanisms = "analysis/tables/enso_rain_mechanisms.yml",
  # Explosive guard: the constant cap 2e5 and
  # both modes.
  thresholds.explosive_case_cap = 2e5,
  derivation.explosive_modes = c("drop", "none"),
  # af_posterior guard scope and unidentified
  # threshold: back-filled with the defaults, as
  # arm_alignment: the stamp covers the whole
  # block, so ("group", 1) would mark every
  # artefact stale, although only the four
  # af_posterior-family steps read them, and
  # their schema tokens re-derive those.
  derivation.af_guard_scope = "unit",
  thresholds.af_noise_share_max = 0.8,
  # Per-gid monthly explosive cap: none.
  thresholds.explosive_gid_month_cap = Inf,
  # Figure 1d's rain_enso method: only that step
  # reads it.
  derivation.rain_enso.baseline = "1990_2019",
  derivation.rain_enso.lags = c(-11L, 11L),
  derivation.rain_enso.dry_mm_month = 10,
  derivation.rain_enso.min_n = 30L,
  derivation.rain_enso.season_share = 0.5,
  derivation.rain_enso.season_min_annual_mm = 100,
  derivation.rain_enso.season_pop_year = 2020L,
  # Supplement fig. SE's seasonality method:
  # only that step reads it.
  derivation.seasonality.admin_level = "ADM1-ADM2",
  derivation.seasonality.min_n = 5L,
  derivation.seasonality.interval = 0.95,
  # Supplement fig. WR's method:
  # only report_who_enso_ratio reads it, so every
  # other artefact stays fresh at the defaults.
  derivation.who_enso_ratio.groups = "enso_west",
  derivation.who_enso_ratio.year_definitions = c(
    "annual_mean",
    "lag_occ",
    "djf_start",
    "djf_end"
  ),
  derivation.who_enso_ratio.first_year = 2010L,
  derivation.who_enso_ratio.censor_weeks = 8L,
  derivation.who_enso_ratio.pseudo_count = 1,
  derivation.who_enso_ratio.interval = 0.95,
  derivation.who_enso_ratio.n_boot = 2000L,
  derivation.who_enso_ratio.seed = 20261107L,
  # Confirmed-subset agreement table: only
  # report_confirmed_agreement reads it.
  derivation.confirmed_agreement.min_confirmed = 1L,
  derivation.confirmed_agreement.min_neutral = 2L,
  derivation.confirmed_agreement.min_phase = 5L,
  derivation.confirmed_agreement.window = "confirmed",
  derivation.confirmed_agreement.boot_b = 1000L,
  derivation.confirmed_agreement.seed = 20261110L,
  derivation.confirmed_agreement.level = 0.95,
  # Lag averaging: only that step reads it; off.
  derivation.lag_averaging.enable = FALSE,
  derivation.lag_averaging.methods = c("stacking", "pseudobma_plus"),
  derivation.lag_averaging.lag_sets = c("all", "eligible"),
  derivation.lag_averaging.k_threshold = 0.7,
  derivation.lag_averaging.report_min_weight = 0.01,
  # Per-capita rate population: pop_ref_year.
  derivation.rate_pop_year = "ref",
  # ENSO phase rule of report_smoothing: the
  # threshold phase.
  derivation.enso_phase.smoothing = "threshold"
)

#' Fit groups of a forecast manifest (names(man$groups)), memoised on the
#' path's mtime; NULL when the manifest cannot be read.
.results_manifest_groups_memo <- new.env(parent = emptyenv())
.results_manifest_groups <- function(path) {
  path <- as.character(unlist(path))
  if (
    length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)
  ) {
    return(NULL)
  }
  key <- paste(
    normalizePath(path, mustWork = FALSE),
    format(file.mtime(path), "%Y%m%d%H%M%OS6")
  )
  if (!is.null(.results_manifest_groups_memo[[key]])) {
    return(.results_manifest_groups_memo[[key]])
  }
  g <- tryCatch(
    names(yaml::read_yaml(path)$groups) %||% character(0),
    error = function(e) NULL
  )
  if (!is.null(g)) {
    assign(key, g, envir = .results_manifest_groups_memo)
  }
  g
}

#' Stamped leaves compared after a normalisation.
#' Named list keyed by the dotted leaf path; each value is `function(x, res)`
#' returning the normalised value, where `res` is the stamped record (the
#' stamp's `results` for the artefact side, the config's for the other) so it
#' can read `res$inputs$manifest`. .results_decision_diff() applies it to both
#' sides; when it cannot run on either side (it errors or returns NULL, e.g. an
#' unreadable manifest), that leaf is compared strictly.
#'
#' scope.exclude_groups: only the names that are fit groups of the manifest
#' matter (an excluded name that is no group is a no-op), so artefacts
#' stamped ["enso_central"] (a group of the legacy registry only) stay fresh
#' under [], while a manifest that has the group still re-derives.
#'
#' inputs.episode_national_zero_threshold: compared as an integer, so a stamp
#' holding the quoted string ("60") and the unquoted YAML (60) read the same
#' decision.
.RESULTS_STAMP_NORMALISE <- list(
  scope.exclude_groups = function(x, res) {
    g <- .results_manifest_groups(res$inputs$manifest)
    if (is.null(g)) {
      return(NULL)
    }
    sort(intersect(as.character(unlist(x)), g))
  },
  inputs.episode_national_zero_threshold = function(x, res) as.integer(x),
  # YAML [] reads as list(), the default is character(0): both mean "off"
  derivation.hindcast.case_skill_supports = function(x, res) {
    as.character(unlist(x))
  }
)

#' Apply .RESULTS_STAMP_NORMALISE to the leaves `have` (stamp record
#' `res_have`) and `want` (config record `res_want`); a leaf whose normaliser
#' fails on either side is left as it is on both.
.results_normalise_leaves <- function(have, want, res_have, res_want) {
  for (k in names(.RESULTS_STAMP_NORMALISE)) {
    f <- .RESULTS_STAMP_NORMALISE[[k]]
    nh <- tryCatch(f(have[[k]], res_have), error = function(e) NULL)
    nw <- tryCatch(f(want[[k]], res_want), error = function(e) NULL)
    if (is.null(nh) || is.null(nw)) {
      next
    }
    have[k] <- list(nh)
    want[k] <- list(nw)
  }
  list(have = have, want = want)
}

.results_decision_diff <- function(stamp, res) {
  if (is.null(stamp)) {
    return("<unstamped artefact: derived before the provenance contract>")
  }
  have <- .results_leaves(stamp$results %||% list())
  want <- .results_leaves(res[.RESULTS_STAMPED_BLOCKS])
  # Back-fill only keys the current config states (it always does, through the
  # defaults): a back-filled key the config lacks would be a spurious diff.
  for (k in names(.RESULTS_STAMP_BACKCOMPAT)) {
    if (is.null(have[[k]]) && !is.null(want[[k]])) {
      have[[k]] <- .RESULTS_STAMP_BACKCOMPAT[[k]]
    }
  }
  nl <- .results_normalise_leaves(have, want, stamp$results %||% list(), res)
  have <- nl$have
  want <- nl$want
  keys <- union(names(have), names(want))
  # Step-scoped leaves of other kinds are not this artefact's decisions (an
  # artefact stamped before `kind` existed is no step-scoped kind).
  other <- unlist(.RESULTS_STEP_DECISIONS[setdiff(
    names(.RESULTS_STEP_DECISIONS),
    stamp$kind %||% ""
  )])
  if (length(other)) {
    keys <- keys[
      !vapply(
        keys,
        function(k) any(k == other | startsWith(k, paste0(other, "."))),
        logical(1)
      )
    ]
  }
  differing <- Filter(
    function(k) !isTRUE(all.equal(have[[k]], want[[k]])),
    keys
  )
  if (!length(differing)) {
    return(character(0))
  }
  vapply(
    differing,
    function(k) {
      sprintf(
        "%s (artefact: %s | config: %s)",
        k,
        paste(format(have[[k]] %||% "<absent>"), collapse = ","),
        paste(format(want[[k]] %||% "<absent>"), collapse = ",")
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' 09 side: return the cached table if it is fresh AND was derived under the
#' current decisions, else NULL so the caller re-derives.
results_read_if_fresh <- function(
  path,
  inputs = character(0),
  res = NULL,
  redo = FALSE
) {
  if (isTRUE(redo)) {
    return(NULL)
  }
  reason <- .results_stale_reason(path, inputs)
  if (!is.null(reason)) {
    if (file.exists(path)) {
      message("[09] re-deriving (", reason, "): ", path)
    }
    return(NULL)
  }
  obj <- qs2::qs_read(path)
  if (!is.null(res)) {
    diff <- .results_decision_diff(get_results_provenance(obj), res)
    if (length(diff)) {
      message(
        "[09] re-deriving (decisions changed: ",
        paste(diff, collapse = "; "),
        "): ",
        path
      )
      return(NULL)
    }
  }
  message("[09] loading results artefact -> ", path)
  obj
}

#' 09 side: stamp and write. Every artefact carries the resolved decisions, the
#' config it came from, its schema token, a path/mtime/size row per input, and
#' the glob specs those inputs were enumerated with (so phase 10 can re-check
#' all of it without knowing what any input file is).
results_write <- function(
  obj,
  path,
  res,
  config_path,
  schema,
  inputs = character(0),
  globs = list(),
  kind = NULL
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  obj <- attach_results_provenance(
    obj,
    .results_stamp(res, config_path, schema, inputs, globs, kind = kind)
  )
  qs2::qs_save(obj, path)
  message("[09] wrote results artefact -> ", path)
  obj
}

#' 10 side: read or stop. NEVER derives -- phase 10 plots what phase 09 made,
#' so a missing, stale, or differently-decided artefact is an error naming the
#' exact command that fixes it, not a silent 40-minute re-derivation inside a
#' figure script.
#'
#' Takes NO input list: everything it checks comes off the artefact's own stamp
#' (see the section header). That is what keeps 10 free of any fit-archive
#' accessor while still refusing to plot a stale table.
#' TRUE when the environment variable RESULTS_ALLOW_STALE_INPUTS is set to a
#' true value: results_require() then WARNS instead of stopping when an
#' artefact's recorded inputs changed, were added or no longer exist since it
#' was derived (e.g. a re-prep of the panels that this rerun must not pick up,
#' or artefacts pulled from the cluster without the simulator partitions).
#' The stamp, schema and decision checks still stop. Opt-in per run (10's stage
#' driver: --allow-stale-inputs); 09's freshness rule is unaffected.
results_allow_stale_inputs <- function() {
  isTRUE(as.logical(Sys.getenv("RESULTS_ALLOW_STALE_INPUTS", "false")))
}

results_require <- function(path, config_path, res = NULL) {
  rerun <- sprintf(
    "  run: Rscript analysis/09b_postprocess_results.R -C %s",
    config_path
  )

  if (!file.exists(path)) {
    stop(
      sprintf("[10] missing results artefact: %s\n%s", path, rerun),
      call. = FALSE
    )
  }

  obj <- qs2::qs_read(path)
  stamp <- get_results_provenance(obj)

  if (is.null(stamp)) {
    stop(
      sprintf(
        paste0(
          "[10] results artefact carries no provenance stamp: %s\n",
          "     (derived before the provenance contract, or written ",
          "by something other than results_write())\n%s"
        ),
        path,
        rerun
      ),
      call. = FALSE
    )
  }

  reason <- .results_stamp_stale_reason(stamp)
  if (!is.null(reason)) {
    if (!results_allow_stale_inputs()) {
      stop(
        sprintf(
          "[10] stale results artefact: %s\n     (%s)\n%s",
          path,
          reason,
          rerun
        ),
        call. = FALSE
      )
    }
    warning(
      sprintf(
        "[10] RESULTS_ALLOW_STALE_INPUTS: using results artefact despite stale inputs: %s\n     (%s)",
        path,
        reason
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  if (!is.null(res)) {
    diff <- .results_decision_diff(stamp, res)
    if (length(diff)) {
      stop(
        sprintf(
          paste0(
            "[10] results artefact was derived under different decisions: %s\n",
            "     differing: %s\n%s"
          ),
          path,
          paste(diff, collapse = "\n                "),
          rerun
        ),
        call. = FALSE
      )
    }
  }
  message("[10] loading results artefact -> ", path)
  obj
}

# ---------------------------------------------------------------------------
# Deployed-fit archive + simulator-partition aggregates
#
# Shared by the manuscript derivation (09) and the counterfactual figure script
# (07d_simulator_figures_counterfactual.R); every input is an explicit argument.
# ---------------------------------------------------------------------------

#' The archive of currently-deployed per-track/per-group best-lag fits
#' (`deploy_prod_run.R --stage deploy_fits`): `<best_lag_dir>/MANIFEST.csv`,
#' keyed by (track, group, kind). The directory rides along as an attribute so
#' track_fit_files() needs no extra argument.
load_best_lag_manifest <- function(best_lag_dir) {
  mnf <- readr::read_csv(
    file.path(best_lag_dir, "MANIFEST.csv"),
    col_types = readr::cols(.default = "c"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(lag = as.integer(lag))
  attr(mnf, "best_lag_dir") <- best_lag_dir
  mnf
}

#' Resolve one track/group's deployed fit + cb_meta path from the archive.
#' `g` is an enso_grouping-style group name ("enso_east"); the archive's
#' `group` column drops the underscore ("ensoeast").
track_fit_files <- function(
  mnf,
  track,
  g,
  best_lag_dir = attr(mnf, "best_lag_dir")
) {
  want_track <- unname(TRACK_CODE[[track]])
  want_group <- gsub("_", "", g)
  want_fit_kind <- paste0(TRACK_KIND_PREFIX[[track]], "fit")
  want_meta_kind <- paste0(TRACK_KIND_PREFIX[[track]], "cb_meta")
  row_fit <- dplyr::filter(
    mnf,
    track == want_track,
    group == want_group,
    kind == want_fit_kind
  )
  row_meta <- dplyr::filter(
    mnf,
    track == want_track,
    group == want_group,
    kind == want_meta_kind
  )
  if (!nrow(row_fit) || !nrow(row_meta)) {
    stop(sprintf(
      "[results] ENSO exposure-response (%s): no deployed fit/meta in %s for group %s",
      track,
      file.path(best_lag_dir, "MANIFEST.csv"),
      g
    ))
  }
  list(
    fit = file.path(best_lag_dir, row_fit$file[1]),
    meta = file.path(best_lag_dir, row_meta$file[1]),
    lag = row_fit$lag[1]
  )
}

#' Sums (or maxima) of `sums` (named: output name = input column) per
#' combination of `keys`, as dplyr's group_by(keys) + summarise(.groups =
#' "drop"): one row per combination, sorted by the keys as dplyr sorts them
#' (C-locale radix order, factors by level, missing keys LAST), the keys then
#' the statistics, a plain tibble. data.table, because 09's replicate-level
#' passes over a 10^4-replicate run reach ~5 x 10^7 groups: load_sim_run_monthly()'s
#' ADM1 x month pass takes ~489 s in dplyr vs ~4 s here, per arm. `op`: "sum" or "max"
#' (both data.table GForce). An empty `df`
#' takes the dplyr path (data.table returns no statistic columns for zero
#' groups). Identical output: tests/test_results_loader_fast.R.
.sum_by_keys <- function(df, keys, sums, op = c("sum", "max")) {
  op <- match.arg(op)
  if (!nrow(df)) {
    f <- if (op == "sum") quote(sum) else quote(max)
    return(
      df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
        dplyr::summarise(
          !!!lapply(sums, function(v) rlang::call2(f, rlang::sym(v))),
          .groups = "drop"
        )
    )
  }
  dt <- data.table::as.data.table(df[, c(keys, unname(sums)), drop = FALSE])
  out <- if (op == "sum") {
    dt[, lapply(.SD, sum), keyby = keys, .SDcols = unname(sums)]
  } else {
    dt[, lapply(.SD, max), keyby = keys, .SDcols = unname(sums)]
  }
  data.table::setnames(out, unname(sums), names(sums))
  data.table::setkey(out, NULL)
  # keyby puts missing keys first; dplyr puts them last.
  if (anyNA(out[, keys, with = FALSE])) {
    data.table::setorderv(out, keys, na.last = TRUE)
  }
  tibble::as_tibble(as.data.frame(out))
}

#' One gid per ADM1 of `gids` (the first), plus every gid without a parseable
#' ADM1 or ISO3 (kept whole, so the loader meets them as it did the full grid):
#' load_sim_run_monthly()'s keep_gids for utils/sim_io.R::sim_monthly_cells().
adm1_rep_gids <- function(gids) {
  gids <- as.character(gids)
  a1 <- gid_to_adm1(gids, strict = FALSE)
  iso <- stringr::str_extract(flatten_gid(gids), "^[A-Z]{3}")
  gids[is.na(a1) | is.na(iso) | !duplicated(a1)]
}

#' Load one simulator run's `monthly_summary.qs` partitions
#' (`<sim_root>/<run_tag>/parts/**/_SUCCESS`, replicate-chunk dirs included).
#' `window_start`/`window_end` bound the ADM1-level table only; the
#' country-level table is returned in full.
#'
#' Cases are the mon track's `cases_month` (one draw per ACTIVE month -- the
#' manuscript estimand; needs a run with grouping.monthly_enable: true).
#' `series` is passed to
#' results_select_series() ("latent" keeps every row).
#'
#' `prefix` (NULL, or a smoothing_prefix_spec()): months before the run's own
#' launch taken from the historical smoothing, prepended per partition as
#' gid-month rows of the same grid, paired replicate by replicate through the
#' partition's handoff_map.csv (see load_sim_run_stitched()). Only then does the
#' return carry a `stitch` element; without it the return is unchanged.
#'
#' @return list(country_month = (replicate_id, month, n_cases, n_outbreaks,
#'   country, max_gid_cases (the largest gid's cases), run), adm1_month = (replicate_id, adm1_unit, month, cases,
#'   outbreaks, run), gids = (gid, country, run): the distinct simulated gids,
#'   the only gid-level information kept, for population denominators).
load_sim_run_monthly <- function(
  sim_root,
  run_tag,
  window_start,
  window_end,
  series = "latent",
  tag = "[results]",
  horizon_by_country = NULL,
  keep_beyond = FALSE,
  prefix = NULL
) {
  parts_dir <- file.path(sim_root, run_tag, "parts")
  if (!dir.exists(parts_dir)) {
    stop(tag, " parts dir not found for run '", run_tag, "': ", parts_dir)
  }
  part_dirs <- dirname(partition_sentinels(parts_dir)) # utils/paths.R
  if (length(part_dirs) == 0L) {
    stop(tag, " no successful TU partitions under ", parts_dir)
  }
  out_country_month <- vector("list", length(part_dirs))
  out_adm1_month <- vector("list", length(part_dirs))
  out_gids <- vector("list", length(part_dirs))
  out_stitch <- vector("list", length(part_dirs))
  for (i in seq_along(part_dirs)) {
    tu_dir <- part_dirs[i]
    # Only the cells a sum can see (utils/sim_io.R::sim_monthly_cells()): a
    # cell with an onset or a case draw, plus every cell of one gid per ADM1
    # (adm1_rep_gids()), which keeps every (replicate, country / ADM1, month)
    # group of the full grid. A row subset of the long table in its order, so
    # the sums, maxima and key sets below are those of the full table. Format-2
    # files are read without expanding the grid (~94% inactive cells).
    cells <- sim_monthly_cells(
      qs2::qs_read(file.path(tu_dir, "monthly_summary.qs")),
      value_cols = c("Y", "cases_month"),
      keep_gids = adm1_rep_gids
    )
    ms_raw <- cells$rows %>%
      results_select_series(
        series,
        paste0("monthly_summary.qs under ", tu_dir),
        kind = "grid"
      )
    if (!"cases_month" %in% names(ms_raw)) {
      stop(
        tag,
        " monthly_summary.qs under ",
        tu_dir,
        " has no cases_month column",
        " -- run '",
        run_tag,
        "' was not configured with grouping.monthly_enable: true."
      )
    }
    # `month` arrives as a character column with a handful of distinct values
    # repeated across millions of rows; parse each distinct string once.
    ms <- ms_raw %>%
      dplyr::transmute(
        replicate_id = as.integer(replicate_id),
        gid = as.character(gid),
        month = local({
          u <- unique(month)
          as.Date(u, format = "%Y-%m-%d")[match(month, u)]
        }),
        Y = as.integer(Y),
        cases = ifelse(is.na(cases_month), 0, as.numeric(cases_month))
      )
    # The full grid's distinct gids / replicates / months (first-seen order):
    # `ms` holds only some of its rows.
    ms_gids <- as.character(cells$gids)
    ms_months <- as.Date(as.character(cells$months), format = "%Y-%m-%d")
    pre <- NULL
    if (!is.null(prefix)) {
      if (isTRUE(prefix$replace)) {
        # Smoothing-as-factual: the smoothing covers the partition's own months.
        pm <- sort(unique(ms_months))
        prefix$start <- min(pm)
        prefix$months <- pm
        prefix$end_excl <- seq(max(pm), by = "month", length.out = 2L)[2L]
      }
      pre <- .smoothing_prefix_rows(
        prefix,
        tu_dir,
        ms_gids,
        as.integer(cells$replicate_ids),
        tag
      )
      if (isTRUE(prefix$replace)) {
        ms <- ms[0L, , drop = FALSE]
        ms_gids <- character(0)
      } else if (any(ms_months < prefix$end_excl)) {
        stop(
          tag,
          " ",
          tu_dir,
          " simulates months before ",
          format(prefix$end_excl),
          ", the start of its own forecast: the smoothing prefix would overlap it.",
          call. = FALSE
        )
      }
      out_stitch[[i]] <- pre$info
      if (isTRUE(prefix$replace)) out_stitch[[i]]$months <- list(pm)
    }
    # Country from the GID's ISO3 prefix, exactly as adm1_unit below -- NOT
    # from the partition directory name, which names a region (tu=enso_east)
    # under partition_by = enso_region. A gid with no ISO3 prefix is an error, not a silent NA key.
    # Parsed once per distinct gid (~10^3), not per row (~10^8 at 10^4 replicates).
    # The smoothing prefix's gids (the partition's own, in its order) first.
    u_gid <- unique(c(if (!is.null(pre)) unique(pre$rows$gid), ms_gids))
    u_country <- stringr::str_extract(flatten_gid(u_gid), "^[A-Z]{3}")
    if (anyNA(u_country)) {
      stop(
        tag,
        " monthly_summary.qs under ",
        tu_dir,
        " has gid(s) with no ISO3 prefix: ",
        paste(head(u_gid[is.na(u_country)], 5L), collapse = ", ")
      )
    }
    # country is a function of gid: the distinct gids in first-seen order, not a
    # unique() over every row's (gid, country) pair (~13 min at 10^4 replicates).
    out_gids[[i]] <- tibble::tibble(gid = u_gid, country = u_country)
    # Partition units mixing fit groups (v2 manifests): each country stops at
    # its own firm horizon (utils/sim_io.R::censor_sim_by_country()); with
    # keep_beyond the later rows are kept, flagged `beyond_firm_horizon`.
    flag_cols <- if (!is.null(horizon_by_country) && isTRUE(keep_beyond)) {
      "beyond_firm_horizon"
    } else {
      character(0)
    }
    cm_keys <- c("replicate_id", "month", "country", flag_cols)
    am_keys <- c("replicate_id", "adm1_unit", "month", flag_cols)
    aggregate_rows <- function(ms) {
      country_of <- u_country[match(ms$gid, u_gid)]
      if (!is.null(horizon_by_country)) {
        ms <- censor_sim_by_country(
          ms %>% dplyr::mutate(country = country_of),
          horizon_by_country,
          keep_beyond = keep_beyond,
          country_col = "country"
        )
        country_of <- ms$country
      }
      ms$country <- country_of
      cm_i <- .sum_by_keys(ms, cm_keys, c(n_cases = "cases", n_outbreaks = "Y"))
      # The largest single gid's cases in the country-month: the per-gid monthly
      # explosive cap (results.thresholds.explosive_gid_month_cap) reads it.
      # Same keys, same sort, so the rows line up.
      cm_i$max_gid_cases <- .sum_by_keys(
        ms,
        cm_keys,
        c(max_gid_cases = "cases"),
        op = "max"
      )$max_gid_cases
      ma <- dplyr::filter(ms, month >= window_start, month <= window_end)
      ug <- unique(ma$gid)
      ma$adm1_unit <- gid_to_adm1(ug, "load_sim_run_monthly")[match(ma$gid, ug)]
      list(
        cm = dplyr::relocate(cm_i, country, .after = n_outbreaks), # same column order as before
        am = .sum_by_keys(ma, am_keys, c(cases = "cases", outbreaks = "Y"))
      )
    }
    parts <- list()
    if (nrow(ms) || is.null(pre)) {
      parts$sim <- aggregate_rows(ms)
    }
    if (!is.null(pre)) {
      # The prefix aggregated per smoothing replicate, then handed to every
      # forecast replicate paired with it: its months are disjoint from the
      # forecast's, so no (replicate, key, month) mixes the two.
      a <- aggregate_rows(pre$rows)
      parts$pre <- list(
        cm = .to_forecast_reps(a$cm, pre$map),
        am = .to_forecast_reps(a$am, pre$map)
      )
    }
    bind_sorted <- function(k, keys) {
      if (length(parts) == 1L) {
        return(parts[[1L]][[k]])
      }
      x <- dplyr::bind_rows(parts$pre[[k]], parts$sim[[k]])
      x[
        do.call(
          order,
          c(unname(as.list(x[keys])), na.last = TRUE, method = "radix")
        ),
        ,
        drop = FALSE
      ]
    }
    out_country_month[[i]] <- bind_sorted("cm", cm_keys)
    out_adm1_month[[i]] <- bind_sorted("am", am_keys)
  }
  out <- list(
    country_month = dplyr::bind_rows(out_country_month) %>%
      dplyr::mutate(run = run_tag),
    adm1_month = dplyr::bind_rows(out_adm1_month) %>%
      dplyr::mutate(run = run_tag),
    gids = dplyr::distinct(dplyr::bind_rows(out_gids)) %>%
      dplyr::mutate(run = run_tag)
  )
  if (!is.null(prefix)) {
    info <- dplyr::bind_rows(out_stitch)
    out$stitch <- if (isTRUE(prefix$replace)) {
      pm <- sort(unique(as.Date(unlist(info$months), origin = "1970-01-01")))
      message(sprintf(
        "%s %s: factual series = the smoothing over %s..%s (%d month(s), %d replicate(s) paired with %d smoothing replicate(s))",
        tag,
        run_tag,
        format(min(pm)),
        format(max(pm)),
        length(pm),
        sum(info$n_replicates),
        sum(info$n_smoothing_reps)
      ))
      list(
        months = pm,
        parts = info[, setdiff(names(info), "months"), drop = FALSE],
        replace = TRUE
      )
    } else {
      .stitch_summary(info, prefix, tag, run_tag, out)
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Factual series stitched from the smoothing
# ---------------------------------------------------------------------------
# The two forecast arms may launch at different months from the same smoothing
# cache (utils/sim_io.R::manifest_arm_starts()). Over [earlier launch, factual
# launch) the factual arm IS the smoothing output: each forecast replicate takes
# the smoothing replicate it was handed off from (07b's handoff_map.csv), so the
# factual series runs smoothing -> forecast along one replicate. Readers: 09
# (utils/results_derive.R) and 07d (07d_simulator_figures_counterfactual.R),
# both through load_arm_run().

#' read_smoothing_parts_events
#'
#' The sparse event tables of the 08b smoothing partitions (the source
#' 08c_smoothing_compile.R binds): $events, one row per Yout == 1 onset
#' (onsets metric), and $events_mon, one row per (rep, gid, ACTIVE month) with
#' `cases_month` (cases metric). Rows before `filt_start` are dropped and
#' `series` is applied with results_select_series(kind = "events").
#' $rep_ids: every smoothing replicate present before any filter.
read_smoothing_parts_events <- function(
  parts_dir,
  filt_start,
  series = "latent",
  tag = "[results]"
) {
  if (!dir.exists(parts_dir)) {
    stop(tag, " smoothing parts dir not found: ", parts_dir)
  }
  sentinels <- partition_sentinels(parts_dir) # utils/paths.R
  if (!length(sentinels)) {
    stop(
      tag,
      " no _SUCCESS sentinels under ",
      parts_dir,
      " -- run 08b_smoothing_spatpomp.R first."
    )
  }
  parts <- lapply(sentinels, function(s) {
    p <- file.path(dirname(s), "part.qs")
    if (!file.exists(p)) {
      warning(tag, " missing part.qs next to ", s)
      return(NULL)
    }
    qs2::qs_read(p)
  })
  parts <- Filter(Negate(is.null), parts)
  events_raw <- dplyr::bind_rows(lapply(parts, function(p) p$events))
  if (is.null(events_raw) || !nrow(events_raw)) {
    stop(tag, " no Yout=1 events under ", parts_dir)
  }
  events_mon_raw <- dplyr::bind_rows(lapply(parts, function(p) p$events_mon))
  if (is.null(events_mon_raw) || !nrow(events_mon_raw)) {
    stop(
      tag,
      " no events_mon rows under ",
      parts_dir,
      " -- the smoothing run was not ",
      "configured with grouping.monthly_enable: true, so the mon-track cases ",
      "metric cannot be built."
    )
  }
  list(
    events = dplyr::mutate(events_raw, gid = as.character(gid)) %>%
      dplyr::filter(t_start >= filt_start) %>%
      results_select_series(series, paste("smoothing events under", parts_dir)),
    events_mon = dplyr::mutate(events_mon_raw, gid = as.character(gid)) %>%
      dplyr::filter(t_start >= filt_start) %>%
      results_select_series(
        series,
        paste("smoothing mon events under", parts_dir)
      ),
    rep_ids = sort(unique(c(events_raw$rep_id, events_mon_raw$rep_id))),
    parts_dir = parts_dir
  )
}

#' read_handoff_map
#'
#' The forecast-replicate -> smoothing-replicate pairing 07b writes next to
#' each successful partition (handoff_map.csv): one row per (tu, replicate_id).
#' NULL when no partition of the run carries one.
read_handoff_map <- function(sim_root, run_tag) {
  dirs <- dirname(partition_sentinels(file.path(sim_root, run_tag, "parts")))
  fs <- file.path(dirs, "handoff_map.csv")
  fs <- fs[file.exists(fs)]
  if (!length(fs)) {
    return(NULL)
  }
  dplyr::bind_rows(lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
}

#' The prefix specification load_sim_run_monthly() takes: smoothing events
#' (read_smoothing_parts_events()) over [start, end_excl).
#' `replace = TRUE` (smoothing-as-factual reading): the smoothing
#' replaces EVERY month the partition simulates (start / end_excl / months are
#' then taken from the partition's own month grid), so the factual series is
#' the smoothing replicate each forecast replicate was handed off from, over
#' the run's whole window.
smoothing_prefix_spec <- function(
  smoothing_events,
  start,
  end_excl,
  replace = FALSE
) {
  start <- as.Date(start)
  end_excl <- as.Date(end_excl)
  stopifnot(start < end_excl)
  months <- seq(start, end_excl, by = "month")
  list(
    events = smoothing_events$events,
    events_mon = smoothing_events$events_mon,
    rep_ids = smoothing_events$rep_ids,
    parts_dir = smoothing_events$parts_dir,
    start = start,
    end_excl = end_excl,
    months = months[months < end_excl],
    replace = isTRUE(replace)
  )
}

#' One forecast partition's smoothing prefix as monthly_summary-grid rows
#' (replicate_id, gid, month, Y, cases) of the SMOOTHING replicates its forecast
#' replicates were handed off from (replicate_id = the smoothing replicate):
#' the grid is the partition's gids x prefix months x those smoothing
#' replicates, zero where the smoothing replicate has no event (the sparse
#' tables carry only events). `Y` counts onsets (Yout == 1 rows), `cases` sums
#' the active months' cases_month. `map` (replicate_id, smoothing_rep_id) is
#' the partition's own handoff map, one row per forecast replicate; the caller
#' aggregates the rows and hands each aggregate to the forecast replicates
#' paired with it (.to_forecast_reps()), never a forecast-replicate grid (~50x
#' larger at 5000 forecast per 100 smoothing replicates).
.smoothing_prefix_rows <- function(prefix, tu_dir, gids, replicate_ids, tag) {
  mp <- file.path(tu_dir, "handoff_map.csv")
  if (!file.exists(mp)) {
    stop(
      tag,
      " ",
      tu_dir,
      " has no handoff_map.csv: the run predates the per-arm handoff ",
      "(2026-10-20) and its replicates cannot be paired with the smoothing replicates. ",
      "Re-run it, or use results.derivation.arm_alignment: common_start.",
      call. = FALSE
    )
  }
  hm <- utils::read.csv(mp, stringsAsFactors = FALSE)
  hm <- hm[
    hm$replicate_id %in% replicate_ids,
    c("replicate_id", "smoothing_rep_id"),
    drop = FALSE
  ]
  miss_rep <- setdiff(replicate_ids, hm$replicate_id)
  if (length(miss_rep)) {
    stop(
      tag,
      " ",
      mp,
      " pairs no smoothing replicate with ",
      length(miss_rep),
      " forecast replicate(s) (e.g. ",
      paste(head(miss_rep, 5L), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  if (anyDuplicated(hm$replicate_id)) {
    stop(
      tag,
      " ",
      mp,
      " pairs forecast replicate(s) ",
      paste(
        head(unique(hm$replicate_id[duplicated(hm$replicate_id)]), 5L),
        collapse = ", "
      ),
      " with more than one smoothing replicate.",
      call. = FALSE
    )
  }
  foreign <- setdiff(unique(hm$smoothing_rep_id), prefix$rep_ids)
  if (length(foreign)) {
    stop(
      tag,
      " ",
      mp,
      " names smoothing replicate(s) ",
      paste(head(foreign, 5L), collapse = ", "),
      " absent from ",
      prefix$parts_dir,
      ": the forecast was handed off from another ",
      "smoothing run.",
      call. = FALSE
    )
  }
  in_win <- function(ev) {
    ev[
      ev$gid %in%
        gids &
        ev$t_start >= prefix$start &
        ev$t_start < prefix$end_excl,
      ,
      drop = FALSE
    ]
  }
  on <- in_win(prefix$events) %>%
    dplyr::count(rep_id, gid, month = t_start, name = "Y")
  cs <- in_win(prefix$events_mon) %>%
    dplyr::group_by(rep_id, gid, month = t_start) %>%
    dplyr::summarise(cases = sum(cases_month, na.rm = TRUE), .groups = "drop")
  grid <- tidyr::expand_grid(
    rep_id = sort(unique(hm$smoothing_rep_id)),
    gid = as.character(gids),
    month = prefix$months
  )
  rows <- grid %>%
    dplyr::left_join(on, by = c("rep_id", "gid", "month")) %>%
    dplyr::left_join(cs, by = c("rep_id", "gid", "month")) %>%
    dplyr::transmute(
      replicate_id = as.integer(rep_id),
      gid,
      month,
      Y = as.integer(dplyr::coalesce(Y, 0L)),
      cases = dplyr::coalesce(cases, 0)
    )
  list(
    rows = rows,
    map = hm,
    info = tibble::tibble(
      part = tu_dir,
      n_replicates = length(replicate_ids),
      n_smoothing_reps = length(unique(hm$smoothing_rep_id))
    )
  )
}

#' .to_forecast_reps -- rows of `df` keyed by smoothing replicate (replicate_id)
#' copied to every forecast replicate `map` pairs with it (map: replicate_id,
#' smoothing_rep_id, one row per forecast replicate), replicate_id rewritten.
.to_forecast_reps <- function(df, map) {
  rows <- split(seq_len(nrow(df)), df$replicate_id)
  take <- rows[as.character(map$smoothing_rep_id)]
  out <- df[unlist(take, use.names = FALSE), , drop = FALSE]
  out$replicate_id <- rep(as.integer(map$replicate_id), lengths(take))
  out
}

#' The `stitch` element of a stitched load: the prefix months, the replicates
#' paired, and the seam (Africa-total median cases in the last prefix month vs
#' the first forecast month; a large jump flags a handoff artefact, not an
#' error). Printed once per run.
.stitch_summary <- function(info, prefix, tag, run_tag, out) {
  cm <- out$country_month
  last_pre <- max(prefix$months)
  first_fc <- prefix$end_excl
  med <- function(m) {
    x <- cm[cm$month == m, , drop = FALSE]
    if (!nrow(x)) {
      return(NA_real_)
    }
    stats::median(tapply(x$n_cases, x$replicate_id, sum))
  }
  seam <- c(prefix = med(last_pre), forecast = med(first_fc))
  message(sprintf(
    "%s %s: factual series stitched from the smoothing over %s..%s (%d month(s), %d replicate(s) paired with %d smoothing replicate(s)); seam median cases %s -> %s",
    tag,
    run_tag,
    format(min(prefix$months)),
    format(last_pre),
    length(prefix$months),
    sum(info$n_replicates),
    sum(info$n_smoothing_reps),
    format(round(seam[["prefix"]])),
    format(round(seam[["forecast"]]))
  ))
  list(months = prefix$months, parts = info, seam = seam)
}

#' load_sim_run_stitched
#'
#' load_sim_run_monthly() for a FACTUAL run launched at `factual_start`, read
#' over a window starting at `window_start`: when the window starts earlier,
#' the months [window_start, factual_start) come from the smoothing
#' (`smoothing_events`, read_smoothing_parts_events()), paired by replicate
#' through the run's handoff_map.csv; a run without the map is refused, never
#' silently unpaired. With window_start >= factual_start the call is exactly
#' load_sim_run_monthly() (no map read, no `stitch` element).
load_sim_run_stitched <- function(
  sim_root,
  run_tag,
  window_start,
  window_end,
  factual_start,
  smoothing_events = NULL,
  ...
) {
  window_start <- as.Date(window_start)
  factual_start <- as.Date(factual_start)
  if (window_start >= factual_start) {
    return(load_sim_run_monthly(
      sim_root,
      run_tag,
      window_start,
      window_end,
      ...
    ))
  }
  if (is.null(smoothing_events)) {
    stop(
      "load_sim_run_stitched: run '",
      run_tag,
      "' launches at ",
      format(factual_start),
      ", after the window start ",
      format(window_start),
      ", and no smoothing events were ",
      "given to stitch the months between.",
      call. = FALSE
    )
  }
  load_sim_run_monthly(
    sim_root,
    run_tag,
    window_start,
    window_end,
    ...,
    prefix = smoothing_prefix_spec(
      smoothing_events,
      window_start,
      factual_start
    )
  )
}

#' load_arm_run
#'
#' The one place that decides how an arm's run is read: the factual
#' ("forecast") arm through load_sim_run_stitched() from the attribution
#' window start, a counterfactual arm plainly. Every reader of the forecast /
#' counterfactual pair (09's two loaders, 07d's two loops) calls this.
#' `arm_starts`: list(factual, counterfactual) launches
#' (utils/sim_io.R::manifest_arm_starts()), or NULL for a plain read.
load_arm_run <- function(
  sim_root,
  run_tag,
  arm,
  window_start,
  window_end,
  arm_starts = NULL,
  smoothing_events = NULL,
  factual_source = "simulation",
  ...
) {
  arm <- match.arg(arm, c("forecast", "counterfactual"))
  factual_source <- match.arg(factual_source, c("simulation", "smoothing"))
  if (identical(arm, "forecast") && identical(factual_source, "smoothing")) {
    # Smoothing-as-factual (event attribution): every month of the
    # factual run read from the smoothing replicate it was handed off from.
    if (is.null(smoothing_events)) {
      stop(
        "load_arm_run: factual_source = 'smoothing' needs the smoothing events.",
        call. = FALSE
      )
    }
    return(load_sim_run_monthly(
      sim_root,
      run_tag,
      window_start,
      window_end,
      ...,
      prefix = smoothing_prefix_spec(
        smoothing_events,
        window_start,
        seq(as.Date(window_end), by = "month", length.out = 2L)[2L],
        replace = TRUE
      )
    ))
  }
  if (is.null(arm_starts) || identical(arm, "counterfactual")) {
    return(load_sim_run_monthly(
      sim_root,
      run_tag,
      window_start,
      window_end,
      ...
    ))
  }
  load_sim_run_stitched(
    sim_root,
    run_tag,
    window_start,
    window_end,
    arm_starts$factual,
    smoothing_events = smoothing_events,
    ...
  )
}

#' Quantile ribbon (q10/q25/q50/q75/q90) across replicates, unweighted.
summarise_reps <- function(
  df,
  group_cols,
  value_cols = c("n_cases", "n_outbreaks")
) {
  fns <- lapply(
    c(q10 = 0.10, q25 = 0.25, q50 = 0.50, q75 = 0.75, q90 = 0.90),
    function(p) {
      force(p)
      function(x) stats::quantile(x, p, names = FALSE)
    }
  )
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(value_cols), fns, .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
}

#' Paired forecast-vs-counterfactual difference per ADM1 unit. `df` has
#' (run, adm1_unit, replicate_id, month, <value_col>) with run in
#' {"forecast","counterfactual"}; paired because both arms of a group share
#' posterior_strategy: paired + seed, so replicate_id indexes the same draw.
paired_adm1_stats <- function(df, value_col, ci_level = 0.50) {
  ci_alpha <- (1 - ci_level) / 2
  df %>%
    dplyr::group_by(run, adm1_unit, replicate_id) %>%
    dplyr::summarise(value = sum(.data[[value_col]]), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = run,
      values_from = value,
      values_fill = 0
    ) %>%
    dplyr::mutate(diff = forecast - counterfactual) %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::summarise(
      forecast = mean(forecast),
      counterfactual = mean(counterfactual),
      mean_diff = mean(diff),
      diff_ci_low = stats::quantile(diff, ci_alpha, names = FALSE),
      diff_ci_high = stats::quantile(diff, 1 - ci_alpha, names = FALSE),
      prob_forecast_higher = mean(diff > 0),
      n_replicates = dplyr::n(),
      .groups = "drop"
    )
}

#' Ratio / significance / crosswalk columns for a paired_adm1_stats() table.
#' Significance (paired CI excludes zero) is never trusted where both arms are
#' low-signal. `adm1_pop` (adm1_unit, pop_adm1) non-NULL adds per-1e5 incidence.
finish_adm1_ratio <- function(
  stats_df,
  low_signal_thresh,
  country_region,
  adm1_group_horizon,
  adm1_pop = NULL,
  ratio_eps = 0.25
) {
  out <- stats_df %>%
    dplyr::mutate(
      country = stringr::str_extract(adm1_unit, "^[A-Z]{3}"),
      low_signal = (forecast + counterfactual) < low_signal_thresh,
      ratio = (forecast + ratio_eps) / (counterfactual + ratio_eps),
      ratio = ifelse(low_signal, NA_real_, ratio),
      log2_ratio = log2(ratio),
      significant = !low_signal & (diff_ci_low > 0 | diff_ci_high < 0)
    )
  if (!is.null(adm1_pop)) {
    out <- out %>%
      dplyr::left_join(adm1_pop, by = "adm1_unit") %>%
      dplyr::mutate(
        incidence_forecast = 1e5 * forecast / pop_adm1,
        incidence_counterfactual = 1e5 * counterfactual / pop_adm1
      )
  }
  out %>%
    dplyr::left_join(country_region, by = "country") %>%
    dplyr::left_join(adm1_group_horizon, by = "adm1_unit")
}

#' Share of all (i, j) pairs with f[i] > g[j] (.prob_pairs_gt) or f[i] < g[j]
#' (.prob_pairs_lt): mean(outer(f, g, `>`)) / mean(outer(f, g, `<`)) by sorting
#' g and counting, O((n + m) log m) instead of an n x m matrix (10^8 cells per
#' ADM1 at 10^4 replicates per arm). Bit-identical to the outer() form
#' (tests/test_results_loader_fast.R).
.prob_pairs_gt <- function(f, g) {
  n <- as.numeric(length(f)) * length(g)
  sum(findInterval(f, sort(g), left.open = TRUE)) / n
}
.prob_pairs_lt <- function(f, g) {
  n <- as.numeric(length(f)) * length(g)
  (n - sum(findInterval(f, sort(g)))) / n
}

#' Unpaired probability-of-shift per ADM1 unit: cross-compare EVERY forecast
#' replicate against EVERY counterfactual replicate (replicate_id pairing is
#' ignored on purpose). `keep_means = FALSE` returns only prob_larger /
#' prob_smaller (the slim variant 09's sensitivity grid uses).
adm1_prob_stats <- function(df, value_col, threshold, keep_means = TRUE) {
  adm1_prob_stats_cum(adm1_prob_cum(df, value_col), threshold, keep_means)
}

#' adm1_prob_stats() in two stages, so a caller looping over thresholds sums
#' the window totals once: adm1_prob_cum() = per (run, ADM1, replicate) totals
#' of `value_col`; adm1_prob_stats_cum() = the probabilities from them.
adm1_prob_cum <- function(df, value_col) {
  .sum_by_keys(df, c("run", "adm1_unit", "replicate_id"), c(value = value_col))
}

adm1_prob_stats_cum <- function(cum, threshold, keep_means = TRUE) {
  cum %>%
    dplyr::group_by(adm1_unit) %>%
    dplyr::group_modify(function(sub, key) {
      f <- sub$value[sub$run == "forecast"]
      cc <- sub$value[sub$run == "counterfactual"]
      if (length(f) == 0L || length(cc) == 0L) {
        out <- tibble::tibble(prob_larger = NA_real_, prob_smaller = NA_real_)
        if (keep_means) {
          out <- dplyr::bind_cols(
            out,
            tibble::tibble(
              mean_forecast = NA_real_,
              mean_counterfactual = NA_real_,
              n_forecast = length(f),
              n_counterfactual = length(cc)
            )
          )
        }
        return(out)
      }
      out <- tibble::tibble(
        prob_larger = .prob_pairs_gt(f, (1 + threshold) * cc),
        prob_smaller = .prob_pairs_lt(f, (1 - threshold) * cc)
      )
      if (keep_means) {
        out <- dplyr::bind_cols(
          out,
          tibble::tibble(
            mean_forecast = mean(f),
            mean_counterfactual = mean(cc),
            n_forecast = length(f),
            n_counterfactual = length(cc)
          )
        )
      }
      out
    }) %>%
    dplyr::ungroup()
}

#' Country / low-signal / signed-probability / significance columns for an
#' adm1_prob_stats(keep_means = TRUE) table. Low-signal ADM1 units are greyed
#' out; `sig_threshold` is both the magnitude cutoff and the confidence cutoff.
finish_adm1_prob <- function(
  stats_df,
  low_signal_thresh,
  sig_threshold,
  country_region,
  adm1_group_horizon
) {
  stats_df %>%
    dplyr::mutate(
      country = stringr::str_extract(adm1_unit, "^[A-Z]{3}"),
      low_signal = (mean_forecast + mean_counterfactual) < low_signal_thresh,
      signed_prob = dplyr::case_when(
        low_signal ~ NA_real_,
        prob_larger >= prob_smaller ~ prob_larger,
        TRUE ~ -prob_smaller
      ),
      significant = !low_signal &
        (prob_larger >= sig_threshold | prob_smaller >= sig_threshold)
    ) %>%
    dplyr::left_join(country_region, by = "country") %>%
    dplyr::left_join(adm1_group_horizon, by = "adm1_unit")
}
