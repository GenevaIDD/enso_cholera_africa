# 03_prepare_siz.R — size-track model dataset: one row per outbreak episode
#   with case count, population / mean-rate offsets, ENSO and covariates.
#
# Usage:
#   Rscript analysis/03_prepare_siz.R -b 1990_2019
#   Rscript analysis/03_prepare_siz.R -b 1990_2019 --episode_key gid_season
#
# Details:
#   Runs in parallel with 03_prepare_occ.R; its episode table feeds
#   03_prepare_dur.R. Inputs: outbreak_metadata.rds (01), climate_stats /
#   climate_extremes SQLite (02), the population crosswalk and gid population
#   under --pop_source crosswalk (default). Outputs: the size data RDS
#   (SIZE_DAT_PATH), the footprint-members RDS,
#   matched_ob_gids_seasonality_siz.rds. --pop_source and --composite_policy
#   mirror 03_prepare_occ.R so both tracks reduce composite LPs the same way.

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(readr)
library(ggplot2)
library(survival)
library(furrr)
library(DBI)
library(RSQLite)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/population_core.R") # load_gadm_units() -- used by get_nearest_seas_gid() (gid.R), Point 2 spatial fallback
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/footprint.R") # unify_footprint_ids, footprint_member_rows, collapse_to_footprint
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
# User-supplied options
# The population / composite / episode surface is shared with the other three
# 03_prepare_* scripts (utils/provenance.R::prep_cli_options()), so the tracks
# cannot drift apart on a default. --pop_source defaults to `crosswalk`:
# legacy emits no P1 denominator triple, so a legacy size frame cannot serve an
# `estimand.denominator: gid` fit. `legacy` stays available.
#
# Composite policy (only meaningful under crosswalk; inert under legacy):
#   keep           apportion to K member rows. Cases and pop scale by the same
#                  w, so all K siblings carry an IDENTICAL implied attack rate
#                  Y_l/P_l and row-wise PSIS-LOO holds out a row its siblings
#                  already encode.
#   collapse       one row per LP-event at the dominant member, carrying the
#                  FULL Y_l against pop = P_l.
#   drop           remove TRUE composites, apportion the rest.
#   drop_collapse  remove TRUE composites, collapse the rest.
# --episode_key gid_season builds ONE canonical episode per (gid, year_seas)
# shared with 03_prepare_occ.R (utils/footprint.R::build_canonical_episodes()).
option_list <- c(
  list(
    make_option(
      c("-r", "--redo"),
      default = TRUE,
      action = "store",
      type = "logical",
      help = "redo data processing"
    ),
    make_option(
      c("-f", "--make_figures"),
      default = FALSE,
      action = "store",
      type = "logical",
      help = "redo data processing"
    ),
    # Overwrite guard: the size frame / footprint map of another
    # prep stamp (e.g. another data vintage) is not replaced without --force.
    prep_force_option("size frame / footprint map")
  ),
  prep_cli_options("siz")
)

opt <- parse_args(OptionParser(option_list = option_list))
opt <- validate_prep_cli(opt, "siz")

# --pop_source (legacy centroid + sum(pop) vs. crosswalk apportionment by
# w = H/P_l) and --composite_policy (keep | collapse | drop | drop_collapse)
# mirror 03_prepare_occ.R exactly so the two tracks stay reducible the same
# way; inert under --pop_source legacy (the default is crosswalk).
# --episode_key gid_season builds the canonical episode table shared with
# 03_prepare_occ.R.
use_xwalk <- identical(opt$pop_source, "crosswalk")
gid_season <- identical(opt$episode_key, "gid_season")

if (!use_xwalk && !identical(opt$composite_policy, "keep")) {
  message(
    "[siz/composite] pop_source=legacy: composite_policy=",
    opt$composite_policy,
    " is inert (no crosswalk expansion to reduce)"
  )
}

print(opt)

# Session for furrr. Reduced from 10 -> 3, mirroring 03_prepare_occ.R: the
# keep-policy occurrence run OOM-killed at anon-rss ~20.7GB on this box (30GB
# total, ~9-10GB already used by other processes). Same crosswalk-expansion
# code path, same risk under `keep`/`drop`.
plan(multisession, workers = 3)

# 3. Data ----------------------------------------------------------------------

# Data on outbreaks
ob_ts <- readRDS(make_outbreak_metadata_file()) %>%
  ungroup()

# Under --episode_key gid_season both prep scripts must build episodes from the
# SAME location-period set, so size applies occurrence's population guard here
# (03_prepare_occ.R has always dropped LPs with pop < --min_pop at read time)
# and marks TZA mainland ADM0 before the section-5 ADM0 filter -- it is ADM0
# by geometry and never enters the canonical episode table. Both are no-ops on
# the default lp_season path, which stays byte-identical.
if (gid_season) {
  ob_ts <- ob_ts %>%
    dplyr::filter(is.na(pop) | pop >= opt$min_pop) %>%
    mark_tza_mainland_adm0()
}

# Coverage data for all locations
all_ts_dat <- readRDS(make_agg_ts_file()) %>%
  filter(location_period_id %in% ob_ts$location_period_id)

# Seasonality estimates
seas <- read_csv("data/best_seasonality_estimates.csv") %>%
  mutate(gid = map_chr(gid, flatten_gid))

# 4. Match GADM to seasonality ids ---------------------------------------------

# Seasonality-match traceability counts (see the A3 fallback in section 4 /
# policy stamp below).
# Only recomputed when opt$redo rebuilds the match cascade; NA when a cached
# matched_gadm is reloaded instead (mirrors the n_rows_in / NA_integer_ pattern
# used elsewhere for xwalk-only counts).
n_seasgid_spatial_checked <- NA_integer_
n_seasgid_spatial_recovered <- NA_integer_

if (opt$redo) {
  # Matched GADM locations
  matched_gadm <- readRDS("generated_data/matched_ob_gids.rds") %>%
    mutate(gid = map_chr(gid, flatten_gid)) %>%
    distinct(lctn_pr, admin_level, gid)

  # Get missing gids
  matched_gadm <- matched_gadm %>%
    rowwise() %>%
    mutate(
      in_seas_dat = map_lgl(gid, ~ . %in% seas$gid),
      seas_gid = case_when(
        in_seas_dat ~ gid,
        TRUE ~ get_next_gid(gid, ref_gids = seas$gid)
      )
    ) %>%
    ungroup()

  # Define ADM1 as first ADM2 if model run at adm2
  seas_adm1 <- seas %>%
    distinct(gid, run_level) %>%
    mutate(gid_adm1 = map_chr(gid, ~ get_adm1_gid(.))) %>%
    group_by(gid_adm1) %>%
    slice(1) %>%
    ungroup()

  # Match to GADM1
  matched_gadm <- matched_gadm %>%
    mutate(gid_adm1 = map_chr(gid, ~ get_adm1_gid(.))) %>%
    rowwise() %>%
    mutate(
      adm1_match = is.na(seas_gid),
      seas_gid = case_when(
        !adm1_match ~ seas_gid,
        TRUE ~ get_next_gid(gid_adm1, ref_gids = seas_adm1$gid_adm1)
      )
    ) %>%
    ungroup()

  # A3. Spatial nearest-neighbour fallback for tile-level misses --------------
  # Third fallback, after the exact match and the two numeric get_next_gid() heuristics
  # above both fail. Recovers "tile-level misses": the outbreak's COUNTRY does have
  # >=1 gid in the 24-country seasonality panel, but the specific ADM1/ADM2 tile has
  # none and no in-range numeric neighbour either. Restricted to same-country
  # candidates, so an outbreak in a country the panel does not cover at all (the large
  # majority of the stage-5 exclusion -- see notebooks/descriptive_occurrence_size.Rmd's
  # cascade-seas-breakdown chunk) still falls through to the January default in
  # section 5, unchanged.
  #
  # Two-stage, mirroring the numeric cascade above:
  #   3a. own-level:  target's `gid`      vs `seas$gid`           candidates
  #   3b. ADM1 level: target's `gid_adm1` vs `seas_adm1$gid_adm1` candidates
  # so `adm1_match` keeps meaning what section 5's join keys on.
  still_na_own <- is.na(matched_gadm$seas_gid)
  if (any(still_na_own)) {
    nn_own <- get_nearest_seas_gid(
      target_gids = matched_gadm$gid[still_na_own],
      ref_gids = seas$gid,
      gadm_path = "data/geodata/gadm_valid.gpkg"
    )
    matched_gadm$seas_gid[still_na_own] <- nn_own
    matched_gadm$adm1_match[still_na_own] <- is.na(nn_own) # FALSE where recovered (own-level match)
  }

  still_na_adm1 <- is.na(matched_gadm$seas_gid)
  if (any(still_na_adm1)) {
    nn_adm1 <- get_nearest_seas_gid(
      target_gids = matched_gadm$gid_adm1[still_na_adm1],
      ref_gids = seas_adm1$gid_adm1,
      gadm_path = "data/geodata/gadm_valid.gpkg"
    )
    matched_gadm$seas_gid[still_na_adm1] <- nn_adm1
    # adm1_match is already TRUE for these rows; an ADM1-space match keeps it TRUE.
  }

  n_seasgid_spatial_checked <- sum(still_na_own)
  n_seasgid_spatial_recovered <- sum(
    still_na_own & !is.na(matched_gadm$seas_gid)
  )
  message(sprintf(
    "[siz/seas] spatial nearest-neighbour fallback: %d/%d unmatched gid(s) recovered (%d still unmatched -- country outside the seasonality panel)",
    n_seasgid_spatial_recovered,
    n_seasgid_spatial_checked,
    n_seasgid_spatial_checked - n_seasgid_spatial_recovered
  ))

  # Save to reload.
  #
  # TRACK-SPECIFIC PATH, not shared with 03_prepare_occ.R: if one track
  # expands its gid roster and the other does not, a shared file would be won
  # by whichever ran last and the `--redo FALSE` path would silently read the
  # other track's roster.
  saveRDS(
    matched_gadm,
    file = "generated_data/matched_ob_gids_seasonality_siz.rds"
  )
} else {
  seas_cache <- "generated_data/matched_ob_gids_seasonality_siz.rds"
  if (!file.exists(seas_cache)) {
    # The pre-split artefact, accepted once so an existing tree still resumes.
    legacy <- "generated_data/matched_ob_gids_seasonality.rds"
    if (!file.exists(legacy)) {
      stop("No seasonality match cache. Re-run with --redo TRUE.")
    }
    warning(
      "Reading the shared pre-split cache ",
      legacy,
      "; re-run with --redo TRUE to write the size-track copy."
    )
    seas_cache <- legacy
  }
  matched_gadm <- readRDS(seas_cache)
}

# 5. Determine cholera years ---------------------------------------------------

# Seasonal offset based on month of lowest prob of cholera
seas_offset <- seas %>%
  filter(
    gid %in%
      matched_gadm$seas_gid |
      gid %in% seas_adm1$gid[seas_adm1$gid_adm1 %in% matched_gadm$seas_gid]
  ) %>%
  group_by(gid) %>%
  summarise(min_seas = which.min(mean))

# Add the seasonal cutoff month to data
# !!Here we are loosing data for 70 outbreaks at ADM0 level
# TODO add seasonality at ADM0 level to avoid loosing country-level outbreaks
ob_ts2 <- ob_ts %>%
  filter(admin_level != "ADM0") %>%
  left_join(
    matched_gadm %>%
      mutate(location_period_id = lctn_pr) %>%
      select(location_period_id, gid, seas_gid, adm1_match)
  ) %>%
  left_join(
    bind_rows(
      seas_offset %>%
        mutate(adm1_match = FALSE),
      seas_offset %>%
        inner_join(seas_adm1) %>%
        select(gid = gid_adm1, min_seas) %>%
        mutate(adm1_match = TRUE)
    ),
    by = c("seas_gid" = "gid", "adm1_match" = "adm1_match")
  )

# Default min_seas (month of lowest cholera probability) to January rather than
# dropping the episode. Root cause (see notebooks/descriptive_occurrence_size.Rmd's
# `cascade-seas-breakdown` chunk, Table 1a-i): the outbreak corpus spans AFR/EMR (44
# countries) but best_seasonality_estimates.csv only covers 24 -- most unmatched
# episodes are a SCOPE MISMATCH, not a data gap. The spatial nearest-neighbour
# fallback in section 4 above already recovers the genuine "tile-level miss" subset;
# what remains here is almost entirely out-of-panel countries, for which January is a
# documented, neutral default rather than a silent loss. Mirrors the min_seas-default
# idiom in analysis/utils/simulator_core.R:build_location_frame() (lines 171-186).
n_seas_episodes_matched <- sum(!is.na(ob_ts2$min_seas))
n_seas_episodes_defaulted_jan <- sum(is.na(ob_ts2$min_seas))
ob_ts2 <- ob_ts2 %>%
  mutate(min_seas = dplyr::coalesce(min_seas, 1L))
message(sprintf(
  "[siz/seas] min_seas: %d episode(s) matched a seasonality estimate, %d defaulted to January (no gid match)",
  n_seas_episodes_matched,
  n_seas_episodes_defaulted_jan
))

ob_ts2 <- ob_ts2 %>%
  mutate(
    TL_seas = map2_chr(TL, -(min_seas - 1), ~ apply_offset_months(.x, .y)) %>%
      as.Date(),
    TR_seas = map2_chr(TR, -(min_seas - 1), ~ apply_offset_months(.x, .y)) %>%
      as.Date()
  )

# min_seas is never NA past this point (defaulted above), so month_seas / year_seas /
# TL_seas_min / TL_seas_max are defined for every episode -- no drop.
ob_ts2 <- ob_ts2 %>%
  mutate(month_seas = month(TL_seas), year_seas = year(TL_seas)) %>%
  group_by(location_period_id) %>%
  mutate(TL_seas_min = min(TL_seas), TL_seas_max = max(TL_seas)) %>%
  ungroup()

# Run checks
# check_outbreak_data(dat = ob_ts2, ref = ob_ts)

# 6. Filter multiple outbreaks per year ----------------------------------------

# 38 location periods have more than one outbreak per year
multi_outbreaks <- ob_ts2 %>%
  add_count(location_period_id, year_seas) %>%
  filter(n > 1)

if (opt$make_figures) {
  multi_outbreaks %>%
    ggplot() +
    geom_vline(aes(xintercept = TL)) +
    facet_wrap(~location) +
    theme_bw()
}

# !! Keep only first outbreak per year
ob_ts3 <- ob_ts2 %>%
  group_by(location_period_id, year_seas) %>%
  # Select first event
  slice_min(TL) %>%
  ungroup()


# 7. Apportion composite location periods --------------------------------------
#
# THE SEAM IS HERE, and nowhere earlier. Expanding before add_count() above
# would inflate the multi-outbreak count K-fold, and slice_min(TL) would then
# return K tied rows instead of one event. Expanding after the collapse in
# section 10 is too late — that collapse groups on `gid`, which is precisely
# the column being corrected. Between the two is the only correct place.
#
# WHAT IT DOES. A composite LP (e.g. the 7-state Sudan union) is one row whose
# `gid` is whatever single member the centroid happened to land in. It becomes K
# rows, one per member gid, with
#
#     suspected_cases_g = suspected_cases * w[l,g]      w = H[l,g] / P_l
#     pop_g             = H[l,g]
#
# Cases and population scale by the SAME w, so the attack rate is invariant per
# member and read_size_data()'s max_attack_rate guard sees exactly what it saw
# before. A single-gid LP has w = 1 and is a strict no-op.
#
# THE PRICE, stated plainly: K pseudo-observations are treated as K real ones,
# which inflates the effective sample size and hence the precision of the fit.
# That is why `apportioned` is carried through to the model frame — it is the
# hook for the two sensitivity runs (response set NA on apportioned rows; a
# per-row dispersion multiplier keyed on n_apportion_members), neither of which
# is baked into the primary fit.
#
# SEASONALITY. min_seas / year_seas were resolved above from the centroid gid
# and are NOT recomputed per member. The seasonal offset is a regional property
# and every member of a composite LP sits in the same region, so this is a
# deliberate approximation; recomputing it would push each member through the
# `filter(!is.na(min_seas))` above and silently drop the w-share of any member
# without a seasonality match.
if (use_xwalk) {
  xwalk_path <- if (nzchar(opt$crosswalk_path)) {
    opt$crosswalk_path
  } else {
    make_pop_run_file(list(), "pop_crosswalk")
  }

  xwalk_siz <- load_population_crosswalk(xwalk_path)

  # --- explicit drop of crosswalk-unmatched LPs -------------------------------
  # See the matching block in 03_prepare_occ.R for the full rationale. In brief:
  # an LP absent from the crosswalk has no population decomposition, so neither
  # reduction path can act on it -- lp_composite_status() will not call it
  # composite (no basis for the claim) and both apportion_composite_cases() and
  # collapse_composite_lps() pass it through untouched under `keep_unmatched =
  # TRUE`. On the size track that is worse than on occurrence: the row survives
  # to the model frame with its centroid gid and NO crosswalk population, i.e. a
  # multi-ADM1 union's whole case count divided by one member's denominator.
  #
  # No-op on the final row set under the default `drop_collapse` (the unmatched
  # LPs are the multi-ADM1 unions, gid_cover far below --cover_threshold), and
  # containment rather than a fix -- the upstream geometry defect is left in
  # place so the next data refresh re-surfaces it through the counts below.
  #
  # Placed BEFORE the counters so `cases_before` and the conservation gates that
  # key off it describe the matchable corpus only. Mass removed here is reported
  # separately and is deliberately NOT folded into `cases_dropped`, which must
  # keep meaning "mass removed by --composite_policy".
  st_match <- lp_composite_status(
    ob_ts3,
    xwalk_siz,
    lp_col = "location_period_id",
    year_col = "year_seas",
    level_col = "admin_level",
    tiling = opt$composite_tiling,
    cover_threshold = opt$cover_threshold
  )
  # One row per input row, arranged by .row (composite_lp.R:241) — positional
  # mask on ob_ts3.
  unmatched <- !st_match$matched
  unmatched_lp <- sort(unique(as.character(ob_ts3$location_period_id[
    unmatched
  ])))
  cases_unmatched <- sum(ob_ts3$suspected_cases[unmatched], na.rm = TRUE)
  cases_corpus <- sum(ob_ts3$suspected_cases, na.rm = TRUE)
  message(sprintf(
    "[siz/xwalk] crosswalk-unmatched: %d row(s) in %d LP(s) dropped, %.2f cases (%.1f%%)%s",
    sum(unmatched),
    length(unmatched_lp),
    cases_unmatched,
    100 * cases_unmatched / max(cases_corpus, 1),
    if (length(unmatched_lp)) {
      sprintf(" (e.g. %s)", paste(head(unmatched_lp, 5), collapse = ", "))
    } else {
      ""
    }
  ))
  ob_ts3 <- ob_ts3[!unmatched, , drop = FALSE]
  if (!nrow(ob_ts3)) {
    stop(
      "[siz/xwalk] no location period matched the crosswalk — every row was ",
      "dropped. The crosswalk at '",
      xwalk_path,
      "' is empty, keyed on a ",
      "different gid convention, or built for a different corpus."
    )
  }

  n_before <- nrow(ob_ts3)
  n_gid_before <- dplyr::n_distinct(ob_ts3$gid)
  n_lp_before <- dplyr::n_distinct(ob_ts3$location_period_id)
  cases_before <- sum(ob_ts3$suspected_cases, na.rm = TRUE)

  do_drop <- opt$composite_policy %in% c("drop", "drop_collapse")
  do_collapse <- opt$composite_policy %in% c("collapse", "drop_collapse")

  # `drop` (without collapse) filters on dominance first, then hands the
  # survivors to the unchanged apportionment path. `drop_collapse` does both in
  # one call. See lp_composite_status() for why dominance and multiplicity are not
  # interchangeable criteria.
  dropped_lp <- character(0)
  cases_kept <- cases_before
  if (do_drop && !do_collapse) {
    st <- lp_composite_status(
      ob_ts3,
      xwalk_siz,
      lp_col = "location_period_id",
      year_col = "year_seas",
      level_col = "admin_level",
      tiling = opt$composite_tiling,
      cover_threshold = opt$cover_threshold
    )
    dropped_lp <- sort(unique(as.character(ob_ts3$location_period_id[
      st$is_composite
    ])))
    ob_ts3 <- ob_ts3[!st$is_composite, , drop = FALSE]
    cases_kept <- sum(ob_ts3$suspected_cases, na.rm = TRUE)
  }

  if (do_collapse) {
    # Snapshot the PRE-collapse mass, keyed by LP. Measuring `cases_kept` after
    # the collapse instead would make the conservation gate below compare a
    # quantity to itself — vacuously true, and blind to a collapse that silently
    # scaled the response.
    pre_lp <- as.character(ob_ts3$location_period_id)
    pre_cases <- as.numeric(ob_ts3$suspected_cases)

    ob_ts3 <- collapse_composite_lps(
      ob_ts3,
      xwalk_siz,
      drop_composites = do_drop,
      cover_threshold = opt$cover_threshold,
      lp_col = "location_period_id",
      year_col = "year_seas",
      gid_col = "gid",
      level_col = "admin_level",
      tiling = opt$composite_tiling,
      keep_unmatched = TRUE
    )
    rep <- attr(ob_ts3, "collapse_report")
    dropped_lp <- attr(ob_ts3, "collapse_dropped_lp")
    n_if_keep <- rep$n_rows_if_keep
    n_lp_comp <- rep$n_lp_composite
    max_mem <- rep$max_members_orig
    # The mass the collapse was SUPPOSED to carry through: everything except the
    # LPs it dropped, measured before it ran.
    cases_kept <- sum(pre_cases[!pre_lp %in% dropped_lp], na.rm = TRUE)
  } else {
    ob_ts3 <- apportion_composite_cases(
      ob_ts3,
      xwalk_siz,
      lp_col = "location_period_id",
      cases_col = "suspected_cases",
      year_col = "year_seas",
      gid_col = "gid",
      level_col = "admin_level",
      tiling = opt$composite_tiling,
      keep_unmatched = TRUE
    )
    rep <- attr(ob_ts3, "apportion_report")
    n_if_keep <- if (do_drop) NA_integer_ else rep$n_rows_out
    n_lp_comp <- if (do_drop) length(dropped_lp) else NA_integer_
    max_mem <- rep$max_members
  }

  cases_after <- sum(ob_ts3$suspected_cases, na.rm = TRUE)
  cases_dropped <- cases_before - cases_kept

  # --- reporting --------------------------------------------------------------
  message(sprintf(
    "[siz/composite] policy=%s threshold=%.2f",
    opt$composite_policy,
    opt$cover_threshold
  ))
  message(sprintf(
    "[siz/composite] LPs: %d total, %s true composite, %d dropped",
    n_lp_before,
    if (is.na(n_lp_comp)) {
      "n/a (not classified under this policy)"
    } else {
      sprintf("%d (%.1f%%)", n_lp_comp, 100 * n_lp_comp / max(n_lp_before, 1L))
    },
    length(dropped_lp)
  ))
  message(sprintf(
    "[siz/composite] rows: %d in -> %d out%s",
    n_before,
    nrow(ob_ts3),
    if (is.na(n_if_keep)) {
      ""
    } else {
      sprintf(
        "   (keep would have given %d, %.1fx)",
        n_if_keep,
        n_if_keep / max(n_before, 1L)
      )
    }
  ))
  message(sprintf(
    "[siz/composite] gids: %d -> %d distinct; max members = %d",
    n_gid_before,
    dplyr::n_distinct(ob_ts3$gid),
    max_mem
  ))
  message(sprintf(
    "[siz/composite] cases: %.2f retained, %.2f dropped (%.1f%%)",
    cases_kept,
    cases_dropped,
    100 * cases_dropped / max(cases_before, 1)
  ))
  message(sprintf(
    "[siz/apportion] %d LP(s) matched, %d unmatched",
    rep$n_lp_matched,
    rep$n_lp_unmatched
  ))

  # Conservation gate (i), narrowed to the RETAINED LPs. Whatever the policy,
  # the reduction step must never create or destroy cases among the LPs it kept:
  # apportionment redistributes them (a failure means sum(w) != 1 for some LP,
  # exactly what the sliver filter causes if w is carried rather than re-derived)
  # and collapse passes them through untouched. Mass leaves ONLY by an explicit
  # drop, and that quantity is asserted separately below rather than absorbed
  # into the tolerance.
  if (
    !isTRUE(abs(cases_after - cases_kept) <= 1e-6 * max(1, abs(cases_kept)))
  ) {
    stop(sprintf(
      "[siz/composite] cases NOT conserved among retained LPs: %.6f before, %.6f after",
      cases_kept,
      cases_after
    ))
  }
  if (
    !isTRUE(
      abs((cases_kept + cases_dropped) - cases_before) <=
        1e-6 * max(1, abs(cases_before))
    )
  ) {
    stop(sprintf(
      paste0(
        "[siz/composite] dropped case mass does not reconcile: ",
        "%.6f retained + %.6f dropped != %.6f before"
      ),
      cases_kept,
      cases_dropped,
      cases_before
    ))
  }
  if (!do_drop && cases_dropped != 0) {
    stop(sprintf(
      "[siz/composite] policy=%s dropped %.6f cases but must drop none.",
      opt$composite_policy,
      cases_dropped
    ))
  }
  message(sprintf(
    "[siz/composite] cases conserved among retained LPs: %.2f",
    cases_after
  ))

  if (do_collapse) {
    if (nrow(ob_ts3) > n_before) {
      stop(
        "[siz/composite] collapse GREW the row count (",
        n_before,
        " -> ",
        nrow(ob_ts3),
        "): collapse must emit at most one row per input row."
      )
    }
    dup <- ob_ts3 %>% count(location_period_id, year_seas) %>% filter(n > 1L)
    if (nrow(dup)) {
      stop(
        "[siz/composite] ",
        nrow(dup),
        " (LP, season) pair(s) contribute more ",
        "than one row after collapse — the one-row-per-LP-event invariant broke."
      )
    }
  }

  if (any(!is.finite(ob_ts3$pop) | ob_ts3$pop <= 0)) {
    stop(
      "[siz/composite] ",
      sum(!is.finite(ob_ts3$pop) | ob_ts3$pop <= 0),
      " row(s) have a non-positive population after the composite reduction."
    )
  }
}


# 8. Canonical gid-season episodes (--episode_key gid_season) ------------------
# ONE episode per (gid, year_seas), built by the SAME call 03_prepare_occ.R
# makes (utils/footprint.R::build_canonical_episodes()): onset = earliest member
# TL, end = latest member TR, cases = the finest-resolution member total,
# admin_level = the gid's geometry-derived GADM level. The digest of that
# table is stamped on every prep artefact and checked by every consumer.
ep_digest <- NA_character_
ep_report <- list()
if (gid_season) {
  gid_pop_path_ep <- if (nzchar(opt$gid_pop_path)) {
    opt$gid_pop_path
  } else {
    make_pop_run_file(list(), "pop_gid")
  }
  ep_res <- build_canonical_episodes(
    ob_ts3 %>% filter(!is.na(min_seas)),
    level_of = canonical_level_lookup(gid_pop_path_ep)
  )
  ep_digest <- digest_episodes(ep_res$episodes)
  ep_report <- ep_res$report
  message(sprintf(
    paste0(
      "[siz/episodes] gid_season: %d member row(s) -> %d canonical episode(s) ",
      "(%d multi-member, %d multi-level, %d with an internal gap, %d crossing ",
      "into the next season-year); min_seas reassigned on %d row(s); digest %s"
    ),
    ep_report$n_members,
    ep_report$n_episodes_canonical,
    ep_report$n_multi_member_episodes,
    ep_report$n_multi_level_episodes,
    ep_report$n_gap_episodes,
    ep_report$n_cross_season_episodes,
    ep_report$n_min_seas_reassigned,
    ep_digest
  ))
  size_src <- ep_res$episodes
  size_grp <- c("gid", "year_seas")
} else {
  size_src <- ob_ts3
  size_grp <- c("uid", "location_period_id", "gid", "year_seas")
}

# 9. Add climate information ---------------------------------------------------
# Make object with time
size_dat <- size_src %>%
  filter(!is.na(min_seas)) %>%
  mutate(time = yday(TL_seas)) %>%
  group_by(across(all_of(size_grp))) %>%
  group_split() %>%
  future_map_dfr(~ make_size_data(., baseline = opt$baseline))

# Conservation gate (ii) — THE CLIMATE TRAP.
# make_size_data() joins climate with inner_join (analysis/utils.R:1150-1154), so a
# member gid that has no row in climate_stats is dropped WITHOUT WARNING and its
# w-share of the cases evaporates. Apportionment then stops conserving, quietly.
# This gate is deliberately loud rather than a filter: the fix is to extract
# climate for the newly-visible member gids (02_preprocess_climate.R --gid_source
# crosswalk), not to absorb the loss here.
if (use_xwalk) {
  cases_pre <- sum(size_src$suspected_cases, na.rm = TRUE)
  cases_post <- sum(size_dat$suspected_cases, na.rm = TRUE)
  lost_gids <- setdiff(unique(size_src$gid), unique(size_dat$gid))
  if (!isTRUE(abs(cases_post - cases_pre) <= 1e-6 * max(1, abs(cases_pre)))) {
    stop(sprintf(
      paste0(
        "[siz/climate] the climate inner_join lost %.2f of %.2f cases ",
        "(%.2f%%) across %d gid(s) with no climate row, e.g. %s. ",
        "Re-run 02_preprocess_climate.R with --gid_source crosswalk."
      ),
      cases_pre - cases_post,
      cases_pre,
      100 * (cases_pre - cases_post) / max(cases_pre, 1),
      length(lost_gids),
      paste(head(lost_gids, 5), collapse = ", ")
    ))
  }
  message(sprintf(
    "[siz/climate] cases survive the climate join: %.2f",
    cases_post
  ))
}


# 10. Footprint reconciliation -------------------------------------------------
# Footprint reconciliation: collapse to ONE row per (gid, year_seas, admin_level).
# ADM3 units are relabelled to their parent ADM2 GID_2 upstream (02_preprocess_climate.R),
# so an ADM2 footprint carries ~5-10 near-duplicate location-period rows per season,
# all sharing identical gid-level covariates. A representative uid/location_period_id
# is fixed per (gid, admin_level) BEFORE the collapse so the per-location_period_id
# time_to_last_outbreak recursion (grouped below) tracks each footprint's full
# cross-season series; covariates use first() (siblings identical).
#
# CASES ARE A FLOW: summing them WITHIN an admin level is correct
# (ADM3 children -> footprint total at ADM3 resolution; ADM2 report -> footprint
# total at ADM2 resolution; summing ACROSS admin levels would double-count).
#
# POPULATION IS A STOCK, and summing it is NOT correct. Two failure modes, both
# live in the corpus:
#   * two location periods that both report the whole gid have their populations
#     ADDED, so the denominator doubles (a summed NGA_19_1 reaches 209M against
#     Nigeria's ~220M);
#   * sum() over an EMPTY group returns 0, and 0 then travels downstream as a
#     population (e.g. COD_4_6_1 at pop = 0, and log(0) from there).
# The replacement is a coverage-saturated denominator,
#     pop = min( sum_l H[l,g] , Q_g ),
# so disjoint ADM3 children still sum while overlapping whole-gid reports
# saturate at the canonical gid population: no set of location periods can
# contain more people than the gid does. See reconcile_footprint_pop() in
# analysis/utils.R and analysis/tests/test_footprint_pop_reconciliation.R.
#
# THE (gid, admin_level) ID-UNIFICATION STEP IS LEGACY-ONLY. Picking one
# representative location_period_id/uid across ALL years sharing a gid is sound
# under --pop_source legacy: `gid` there comes from a stable geometric
# containment relationship (ADM3 centroid -> its ADM2 parent), so rows sharing a
# gid really are repeat reports of the SAME footprint. Under --pop_source
# crosswalk, `gid` instead comes from apportion_composite_cases()/
# collapse_composite_lps() picking the crosswalk's DOMINANT POPULATION MEMBER for
# each location period independently -- a choice that can coincide across
# genuinely unrelated location periods whenever their member candidate sets
# share the same dominant polygon (a big town swamping several distant, smaller
# districts' matches). Forcing a single id/uid onto those rows silently merges
# distinct outbreaks' identities. Example: DRC location periods
# 14287/14296/14360/14372/14398 (Lubumbashi/Vangu/Kabalo/Kasenga/Moba health
# districts, formerly Katanga province) all resolve to gid COD.3.6_1 and would
# collapse onto one id (the lowest, 14287) across different years, each then
# carrying an unrelated small district's population as its denominator. Under
# crosswalk, each row keeps its own location_period_id/uid; only the (gid, year_seas,
# admin_level) grouping below still reconciles genuine same-year, same-gid
# duplicates.
# THE MEMBERSHIP MAP MUST BE CAPTURED HERE, before either step below.
#
# The collapse takes every non-summed column with `first()`, so the surviving
# uid / location_period_id is a REPRESENTATIVE LABEL and the identity of the
# other contributing location periods is gone. That is invisible as long as
# nothing tries to attach an external per-location-period quantity to an
# episode -- but 03_prepare_mon.R does exactly that with the monthly case
# series from 01b_preprocess_outbreaks_monthly.R, and joining it on
# (uid, location_period_id) would silently return only the label LP's cases.
# On this corpus that is 314/1581 episodes and 29.7% of the case mass (see the
# header of analysis/utils/footprint.R for the full measurement).
#
# Before `unify_footprint_ids()` specifically, because under legacy that step
# overwrites uid / location_period_id with one representative across ALL
# season-years; a map built afterwards would assign every season's cases to a
# single location period. `footprint_member_rows()` does not modify `size_dat`, so
# capturing the map here leaves the size frame unchanged.
#
# Under --episode_key gid_season the rows ARE already the canonical episodes
# and the membership map is the canonical one (with `report_level` and the
# `in_total` flag of the members whose cases the episode total counts); there
# is nothing left to collapse.
if (gid_season) {
  footprint_members <- ep_res$members
} else {
  footprint_members <- footprint_member_rows(size_dat)

  if (!use_xwalk) {
    size_dat <- unify_footprint_ids(size_dat)
  }

  size_dat <- collapse_to_footprint(size_dat)
}

if (use_xwalk) {
  gid_pop_path <- if (nzchar(opt$gid_pop_path)) {
    opt$gid_pop_path
  } else {
    make_pop_run_file(list(), "pop_gid")
  }
  wpp_path_siz <- if (nzchar(opt$wpp_path)) {
    opt$wpp_path
  } else {
    make_pop_run_file(list(), "pop_wpp")
  }

  gid_pop_siz <- load_gid_population(gid_pop_path)
  wpp_siz <- load_wpp_table(wpp_path_siz)

  pop_before <- size_dat$pop
  recon <- reconcile_footprint_pop(
    size_dat$gid,
    size_dat$year_seas,
    H_sum = size_dat$pop,
    gid_pop = gid_pop_siz,
    wpp = wpp_siz,
    on_missing = "error"
  )
  size_dat$pop <- recon$pop
  size_dat$pop_tier <- recon$pop_tier

  message("[siz/footprint] population tiers after reconciliation:")
  print(table(size_dat$pop_tier, useNA = "ifany"))
  message(sprintf(
    "[siz/footprint] %d unit(s) capped at Q_g; %d rescued from a zero/absent sum",
    sum(recon$pop_tier == "gid_capped"),
    sum(recon$pop_tier == "gid_fallback")
  ))
  message(sprintf(
    "[siz/footprint] total population %.0f -> %.0f (%+.1f%%)",
    sum(pop_before, na.rm = TRUE),
    sum(size_dat$pop),
    100 * (sum(size_dat$pop) / max(sum(pop_before, na.rm = TRUE), 1) - 1)
  ))

  # Gate (iii): zero is never a population. offset(log(pop)) would be -Inf.
  if (any(!is.finite(size_dat$pop) | size_dat$pop <= 0)) {
    stop(
      "[siz/footprint] ",
      sum(!is.finite(size_dat$pop) | size_dat$pop <= 0),
      " unit(s) still have a non-positive population."
    )
  }

  # The one filter apportionment DOES perturb. Attack rate is w-invariant, so
  # read_size_data()'s max_attack_rate guard is untouched; min_pop is an ABSOLUTE
  # floor (utils.R:6049, default 10 people), so a member carrying a tiny w-share
  # of a large LP can fall under it where the un-split LP never would.
  # Reported here, at the point of creation, rather than discovered as an
  # unexplained row-count drop at fit time.
  n_below <- sum(size_dat$pop < 10, na.rm = TRUE)
  if (n_below) {
    message(sprintf(
      paste0(
        "[siz/footprint] %d unit(s) (%.2f%% of cases) fall below ",
        "the min_pop = 10 floor and will be dropped by read_size_data()"
      ),
      n_below,
      100 *
        sum(size_dat$suspected_cases[size_dat$pop < 10], na.rm = TRUE) /
        max(sum(size_dat$suspected_cases, na.rm = TRUE), 1)
    ))
  }
}

# 11. P1: carry BOTH populations, and the coverage between them ----------------
#
# A single `pop` column has its MEANING fixed here and then forgotten:
# everything downstream -- the fit's offset(log(pop)), the mean_rate
# denominator, the simulator frame -- would receive a number with no unit
# attached, so a `Q x rbar^H` mismatch could go unnoticed.
#
# So emit all three SIDE BY SIDE and let the declared estimand pick (P0/P1 in
# resolve_size_pop_column()):
#
#   pop_footprint  H_g = min( sum_l H[l,g] , Q_g )   the reporting footprint --
#                  the population of the area actually under surveillance for
#                  this outbreak. ENDOGENOUS: its extent is a function of where
#                  the outbreak was detected, so a rate built on it is
#                  conditioned on the outcome.
#   pop_gid        Q_g                               the whole admin unit.
#                  Exogenous, defined for every unit and every projection year,
#                  hence the only denominator a forward simulator can use
#                  consistently.
#   pop_cover      c = H_g / Q_g in (0, 1]           the ratio between them, and
#                  therefore the exact factor by which mixing the two rescales a
#                  prediction. Computed at 02_preprocess_crosswalk.R:255 and
#                  discarded at the collapse; recovered here.
#
# `pop` itself is LEFT UNTOUCHED (= today's value) so this block cannot change
# any existing run. The switch is opt-in downstream, never here.
{
  gid_pop_path_p1 <- if (nzchar(opt$gid_pop_path)) {
    opt$gid_pop_path
  } else {
    tryCatch(make_pop_run_file(list(), "pop_gid"), error = function(e) "")
  }

  size_dat$pop_footprint <- size_dat$pop

  if (nzchar(gid_pop_path_p1) && file.exists(gid_pop_path_p1)) {
    gp_p1 <- load_gid_population(gid_pop_path_p1)
    wpp_p1 <- load_wpp_table(
      if (nzchar(opt$wpp_path)) {
        opt$wpp_path
      } else {
        tryCatch(make_pop_run_file(list(), "pop_wpp"), error = function(e) "")
      }
    )
    # `on_missing = "na"`: a gid-year outside the canonical extraction window is
    # a KNOWN GAP, not a corruption. Erroring here would make P1 a breaking
    # change for the legacy path; instead the gap is reported and carried as NA,
    # and resolve_size_pop_column() refuses to fit on it later. Never 0.
    q_res <- resolve_gid_population(
      size_dat$gid,
      size_dat$year_seas,
      gid_pop = gp_p1,
      wpp = wpp_p1,
      on_missing = "na"
    )
    size_dat$pop_gid <- if (is.list(q_res)) q_res$pop else as.numeric(q_res)
    size_dat$pop_cover <- ifelse(
      is.finite(size_dat$pop_gid) & size_dat$pop_gid > 0,
      size_dat$pop_footprint / size_dat$pop_gid,
      NA_real_
    )

    n_na <- sum(!is.finite(size_dat$pop_gid))
    message(sprintf(
      paste0(
        "[siz/P1] pop_footprint + pop_gid + pop_cover emitted. Q_g missing for ",
        "%d/%d row(s) (%.2f%% of cases)."
      ),
      n_na,
      nrow(size_dat),
      100 *
        sum(
          size_dat$suspected_cases[!is.finite(size_dat$pop_gid)],
          na.rm = TRUE
        ) /
        max(sum(size_dat$suspected_cases, na.rm = TRUE), 1)
    ))
    cv <- size_dat$pop_cover[is.finite(size_dat$pop_cover)]
    if (length(cv)) {
      message(sprintf(
        paste0(
          "[siz/P1] coverage c = H/Q: median %.3f  IQR [%.3f, %.3f]  ",
          "case-weighted mean %.3f  |  1/c at the median = %.2fx"
        ),
        stats::median(cv),
        stats::quantile(cv, .25),
        stats::quantile(cv, .75),
        stats::weighted.mean(
          cv,
          size_dat$suspected_cases[is.finite(size_dat$pop_cover)]
        ),
        1 / stats::median(cv)
      ))
    }
    # c > 1 is impossible by construction (reconcile_footprint_pop saturates at
    # Q_g). If it appears, the two populations came from different tables and
    # every coverage number below is meaningless -- say so loudly.
    n_over <- sum(cv > 1 + 1e-6)
    if (n_over) {
      warning(sprintf(
        paste0(
          "[siz/P1] %d row(s) have coverage c > 1, which reconcile_footprint_pop() ",
          "makes impossible. pop_footprint and pop_gid are inconsistent -- do not ",
          "trust pop_cover."
        ),
        n_over
      ))
    }
  } else {
    size_dat$pop_gid <- NA_real_
    size_dat$pop_cover <- NA_real_
    message(paste0(
      "[siz/P1] no canonical gid population table available (",
      "--gid_pop_path unset and default absent); pop_gid/pop_cover ",
      "left NA. size.estimand.denominator = 'gid' will refuse to run."
    ))
  }
}


# 12. Add time to last outbreak ------------------------------------------------

size_dat <- size_dat %>%
  group_by(location_period_id) %>%
  arrange(TL) %>%
  mutate(
    rid = row_number(),
    last_outbreak_start = as.Date("2000-01-01"),
    # Every row is an episode, so every TR counts, the first one included
    # (rid > 1 here skipped it, so episode 2 read the 2000-01-01 sentinel).
    last_outbreak_start = TR,
    last_outbreak_start = cummax(as.numeric(last_outbreak_start)) %>%
      as.Date(origin = "1970-01-01") %>%
      lag(n = 1),
    time_to_last_outbreak = difftime(
      TL,
      last_outbreak_start,
      units = "days"
    ) %>%
      as.numeric()
  ) %>%
  ungroup() %>%
  mark_tza_mainland_adm0() %>%
  # flatten_gid() first: under --pop_source crosswalk, gid comes from the LP x gid
  # crosswalk in DOTTED GADM notation ("AGO.10_1"), not the underscore form this
  # regex expects -- mirrors the fix in 03_prepare_occ.R (same regex, same bug).
  # Unflattened, str_extract() returns NA for every crosswalk-sourced row.
  mutate(ADM1 = gid_to_adm1(gid, "03_prepare_siz"))

# 13. Add ENSO lag columns (enso_lag0 ... enso_lag12) --------------------------
# Joined on outbreak start month (TL). All 13 lags always written so that
# lag-sensitivity configs can select any lag without re-running prep.
enso_dat <- fetch_enso_index(
  source = opt$enso_source,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_enso)
)
size_dat <- join_enso_lags(size_dat, enso_dat, date_col = "TL")

# Add climate lag columns (lag1..lag11) for DLNM use at fit time.
# mode = "event": each row is one outbreak; lags are resolved via a batched
# SQLite lookup on (gid, month) — no season-year grouping needed.
size_dat <- make_climate_lags(
  size_dat,
  max_lag = 11L,
  mode = "event",
  baseline = opt$baseline
)

# 14. Stamp the population provenance (Step 7 ordering guard) ------------------
# build_mean_rate_table() computes every cell as suspected_cases / pop, and the
# resulting table is then persisted and read as a FROZEN sidecar by the fit
# scripts and by simulator_core.R. A population change therefore
# invalidates the whole table -- but both artefacts still look perfectly
# plausible on their own, so nothing downstream can notice the mismatch from
# the numbers alone. The stamp travels with the file and lets
# assert_pop_provenance() refuse the pairing by name.
#
# Required order, enforced by that check:
#   pop fix -> mean_rate rebuild -> size refit -> Phase-A rebuild -> forecast
size_dat <- attach_pop_provenance(
  size_dat,
  pop_provenance_stamp(
    size_dat$pop,
    pop_source = opt$pop_source,
    extra = list(
      crosswalk_path = if (use_xwalk) opt$crosswalk_path else NA_character_,
      gid_pop_path = if (use_xwalk) opt$gid_pop_path else NA_character_,
      # The composite policy is part of the population's provenance, not merely a
      # run setting: it decides which LPs contribute a denominator at all. Without
      # it, `pop_source` reads "crosswalk" under all four policies and a mismatch
      # surfaces only as an unexplained digest change.
      composite_policy = opt$composite_policy,
      cover_threshold = opt$cover_threshold,
      composite_tiling = opt$composite_tiling,
      n_apportioned = if ("apportioned" %in% names(size_dat)) {
        sum(size_dat$apportioned, na.rm = TRUE)
      } else {
        0L
      }
    )
  )
)

# 15. Write the size frame -----------------------------------------------------
SIZE_DAT_PATH <- make_size_dat_file(opt)

# Policy stamp, built BEFORE anything is written. See prep_policy_stamp_path().
# Written for the size frame AND its membership map: the two are a pair, and
# 03_prepare_mon.R checks that the members it aggregates through carry the
# same episode digest as the duration panel it joins onto.
size_stamp <- prep_policy_stamp(
  pop_source = opt$pop_source,
  composite_policy = opt$composite_policy,
  cover_threshold = opt$cover_threshold,
  composite_tiling = opt$composite_tiling,
  pop_digest = get_pop_provenance(size_dat)$pop_digest,
  episode_key = opt$episode_key,
  episode_digest = ep_digest,
  counts = c(
    list(
      track = "size",
      n_lp_total = if (use_xwalk) {
        n_lp_before
      } else {
        dplyr::n_distinct(size_dat$location_period_id)
      },
      n_lp_dropped = if (use_xwalk) length(dropped_lp) else 0L,
      n_rows_in = if (use_xwalk) n_before else NA_integer_,
      n_rows_out = nrow(size_dat),
      cases_dropped = if (use_xwalk) cases_dropped else 0,
      dropped_lp = if (use_xwalk) dropped_lp else character(0),
      # -- Point 1/2 seasonality-match traceability --------------------------
      n_seasgid_spatial_checked = n_seasgid_spatial_checked,
      n_seasgid_spatial_recovered = n_seasgid_spatial_recovered,
      n_seas_episodes_matched = n_seas_episodes_matched,
      n_seas_episodes_defaulted_jan = n_seas_episodes_defaulted_jan
    ),
    # gid_season only (empty list otherwise, so the lp_season stamp keeps
    # exactly its historical field set).
    ep_report
  )
)

# Overwrite guard (the rule of 03_prepare_occ.R): an existing
# size frame or footprint map whose stamp differs in anything but stamped_at
# (another data vintage, another policy) or that has none is not replaced
# without --force. A rebuild under an identical stamp passes.
for (.p in c(SIZE_DAT_PATH, make_footprint_members_file(opt))) {
  assert_prep_stamp_overwrite(
    .p,
    size_stamp,
    force = isTRUE(opt$force),
    rerun_hint = "Rscript analysis/03_prepare_siz.R ... --force"
  )
}

saveRDS(size_dat, SIZE_DAT_PATH)

# 16. The footprint membership map, written as a PAIR with the size frame ------
# Same invocation, same policies, so the two are only ever meaningful together
# -- hence no separate naming parameters (see make_footprint_members_file()).
# Restricted to the footprints that actually survived to the final frame: rows
# lost after the collapse (the climate join, the min_pop / attack-rate guards
# applied later by read_size_data()) would otherwise leave members pointing at
# episodes no consumer can resolve.
{
  FOOTPRINT_MEMBERS_PATH <- make_footprint_members_file(opt)
  surviving <- paste(
    as.character(size_dat$gid),
    as.character(size_dat$year_seas),
    as.character(size_dat$admin_level)
  )
  fm_out <- footprint_members %>%
    filter(paste(gid, year_seas, admin_level) %in% surviving)
  saveRDS(fm_out, FOOTPRINT_MEMBERS_PATH)
  message(sprintf(
    paste0(
      "[siz/footprint-members] wrote %s (%d member row(s) ",
      "over %d footprint(s); %d footprint(s) have >1 ",
      "contributing location period)"
    ),
    FOOTPRINT_MEMBERS_PATH,
    nrow(fm_out),
    dplyr::n_distinct(paste(fm_out$gid, fm_out$year_seas, fm_out$admin_level)),
    dplyr::n_distinct(paste(fm_out$gid, fm_out$year_seas, fm_out$admin_level)[
      fm_out$n_members > 1L
    ])
  ))
}

# Policy sidecar (size_stamp, built above), next to the artefacts it describes.
persist_prep_policy_stamp(SIZE_DAT_PATH, size_stamp)
persist_prep_policy_stamp(FOOTPRINT_MEMBERS_PATH, size_stamp)

# Standalone copy of the stamp, so any reader that only has `opt` (every
# consumer of make_size_run_file(opt, "size_mean_rate") except the size fit
# scripts themselves) can content-address the mean_rate sidecar path without
# paying for a full readRDS() of the table above. See pop_stamp_sidecar_path() /
# current_pop_digest_token() in utils.R.
persist_pop_stamp_sidecar(opt, get_pop_provenance(size_dat))

{
  .st <- get_pop_provenance(size_dat)
  message(sprintf(
    paste0(
      "[siz] population provenance: source=%s, n=%d, total=%.0f, ",
      "digest=%s, apportioned rows=%d"
    ),
    .st$pop_source,
    .st$n,
    .st$pop_total,
    .st$pop_digest,
    .st$n_apportioned
  ))
  if (use_xwalk) {
    message(
      "[siz] the mean_rate sidecar is now STALE. Delete ",
      "generated_data/<...>/size_mean_rate*.rds and re-run the size fit ",
      "BEFORE any forecast: new populations against old rates are ",
      "silently wrong."
    )
  }
}
