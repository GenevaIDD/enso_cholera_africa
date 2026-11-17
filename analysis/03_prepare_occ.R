# 03_prepare_occ.R — occurrence-track model dataset: monthly at-risk survival
#   intervals per location and season-year, with ENSO and contagion covariates.
#
# Usage:
#   Rscript analysis/03_prepare_occ.R -b 1990_2019
#   Rscript analysis/03_prepare_occ.R -b 1990_2019 --episode_key gid_season \
#     --national_zero_fill
#
# Details:
#   Inputs: outbreak_metadata.rds (01), climate_stats / climate_extremes
#   SQLite (02), the population crosswalk under --pop_source crosswalk.
#   Outputs: make_surv_dat_file() (outbreaks_survival_data_base-<baseline>
#   .rds + panel tokens), matched_ob_gids_seasonality_occ.rds, the
#   .min_seas_source.csv season-calendar sidecar and, under gid_season, the
#   adm1_overlap / discarded_onsets / adm0_dropped (/ masked_rows) sidecars.
#   Composite-LP handling and the episode key are documented above the
#   option list; the national-zero and at-risk masks above their checks.

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
source("analysis/utils/footprint.R") # build_canonical_episodes() -- --episode_key gid_season (shared with 03_prepare_siz.R)
source("analysis/utils/provenance.R")
source("analysis/utils/plots.R")
source("analysis/utils/anchor.R")
source("analysis/utils/spatial_contagion.R") # build_adm1_dist_km()/build_contagion_table() -- section 8 distance-weighted neighbour contagion
source("analysis/utils/config.R") # load_grouping() -- registry for --contagion_partition_source
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/season_calendar.R") # gid_min_seas_source() -- per-gid calendar source sidecar
source("analysis/utils/national_zero.R") # national_zero_rule_params() / national_zero_rule_digest() -- --national_zero_rule
source("analysis/utils/coverage.R") # build_coverage_mask() -- --contagion_coverage_report

# 2. User-supplied options -----------------------------------------------------
# The population / composite / episode surface (--baseline, --enso_source,
# --refresh_enso, --min_pop, --pop_source [default crosswalk], --crosswalk_path,
# --gid_pop_path, --wpp_path, --composite_policy, --cover_threshold,
# --composite_tiling, --episode_key) is shared with the other three
# 03_prepare_* scripts: utils/provenance.R::prep_cli_options(). Track-specific
# extras below.
#
# --pop_source: "legacy" = the centroid match (one arbitrary member of a
# composite LP carries the outbreak); "crosswalk" = propagate the outbreak
# indicator to EVERY member gid (see section 6a for why propagation, not
# w-scaling, is the only defensible expansion here).
#
# --composite_policy (crosswalk only): keep = expand to K member rows (K
# pseudo-independent observations); collapse = one row per LP-event at the
# dominant-population member; drop = remove TRUE composites, expand the rest;
# drop_collapse = remove TRUE composites, collapse the rest. `drop` keys on
# DOMINANCE (gid_cover < --cover_threshold, ~22% of LPs), the row inflation on
# MULTIPLICITY (~97% of rows); only collapse removes the latter. See
# lp_composite_status(). --composite_tiling: derived (the crosswalk's
# geometry-based lp_level) | corpus_label (legacy behaviour).
#
# --episode_key gid_season: ONE canonical episode per (gid, year_seas), shared
# with 03_prepare_siz.R (utils/footprint.R::build_canonical_episodes()); the
# survival history is built per gid, not per location period.
option_list <- c(
  list(
    make_option(
      c("-r", "--redo"),
      default = FALSE,
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
    # --- country/TU spillover grouping key ---------------------------------------
    # Which key the section-8 leave-one-out spillover lags (ongoing_outbreak_country_lag1,
    # ongoing_outbreak_tu_lag1) are grouped on.
    #
    #   t_end      legacy (DEFAULT). Groups on
    #              (country[, cluster], t_start, t_end). Because survival
    #              intervals are cut at the onset event, an onset row gets an
    #              outcome-determined t_end that isolates it into its own count
    #              group -- an outcome-leakage artifact.
    #   calendar   leakage-free. Collapses the ADM1 active grid to one row per
    #              (country, ADM1, t_start) first, then does the leave-one-out
    #              at calendar-month (t_start-only) resolution, so the stored
    #              o_ct/o_tu columns become honest.
    make_option(
      "--spillover_grouping",
      default = "calendar",
      type = "character",
      help = "t_end (legacy) | calendar (leakage-free) [default: %default]"
    ),
    # --- distance-weighted neighbour contagion ----------------------------------
    # All four nbr_lag1_r{200,500}[_rn] columns are always written (cf. the 13 ENSO
    # lag columns at section 9), so a rho- or scale-sensitivity config can select
    # any of them without re-running prep. Which one enters the formula is a fit-time
    # decision (occurrence.covariates.nbr_contagion_lag1 + the contagion: block).
    make_option(
      "--contagion_rhos",
      default = "200,500",
      type = "character",
      help = "comma-separated decay bandwidths in km for the neighbour-contagion kernel [default: %default]"
    ),
    make_option(
      "--contagion_gadm_path",
      default = "data/geodata/gadm_valid.gpkg",
      type = "character",
      help = "GADM geopackage supplying ADM1 centroids for the contagion kernel [default: %default]"
    ),
    make_option(
      "--refresh_contagion_dist",
      default = FALSE,
      action = "store_true",
      help = "Rebuild the cached ADM1 centroid distance matrix even if one exists [default: %default]"
    ),
    # Partition-restricted siblings nbr_lag1_r<rho>[_rn]_pg: neighbour set AND
    # row-normalising denominator limited to the unit's group in this registry
    # (analysis/groupings/<source>.csv). Empty = not written (byte-identical
    # output). Recorded in the policy stamp; 04_fit_occ_* refuses a
    # contagion.partition: group fit whose grouping.source differs.
    # A comma list: the FIRST registry still writes the _pg pair;
    # EVERY registry writes nbr_lag1_r<rho>[_rn]_part<token> siblings, the
    # columns of contagion.partition: registry fits (utils/paths.R).
    make_option(
      "--contagion_partition_source",
      default = "",
      type = "character",
      help = "grouping registry/registries (comma list) for the partition-restricted contagion columns, e.g. enso_grouping[,<tu registry>] [default: none]"
    ),
    # --- national zero-fill --------------------------------------------------
    # gid_season only: every national-zero month (01_preprocess_national_zero_months.R)
    # of a panel gid's country that is not yet a panel row becomes a no-onset row
    # (utils/occ_data.R::national_zero_surv_rows()); never at or after an onset in
    # its season-year. Changes the panel's ROWS: "_nz" filename token and a
    # national_zero_fill column.
    make_option(
      "--national_zero_fill",
      default = FALSE,
      action = "store_true",
      help = "add no-onset rows for the national-zero months (gid_season only) [default: %default]"
    ),
    make_option(
      "--national_zero_path",
      default = "",
      type = "character",
      help = "national-zero months file [default: make_national_zero_months_file(<threshold>)]"
    ),
    # Threshold variant: the national-zero months built with
    # 01_preprocess_national_zero_months.R --zero_threshold N (WHO country-years
    # under N cases also count); "_nz<N>" filename token.
    make_option(
      "--national_zero_threshold",
      default = 0L,
      type = "integer",
      help = "with --national_zero_fill: national-zero months of this WHO case threshold (0 = exact zeros) [default: %default]"
    ),
    # Named rule: the months of 01 --rule <name>
    # (analysis/national_zero_rules.yml); "_nz<N>-<name>" filename token.
    make_option(
      "--national_zero_rule",
      default = "",
      type = "character",
      help = "with --national_zero_fill: national-zero rule of analysis/national_zero_rules.yml ('' = legacy) [default: legacy]"
    ),
    # --- occurrence-panel fixes (gid_season only) ------------------------------
    make_option(
      "--at_risk_start",
      default = "season",
      type = "character",
      help = "season | after_episode: mask rows inside an earlier season's running episode (\"_ae\") [default: %default]"
    ),
    # Contagion coverage report (stamp only; the covariate
    # is unchanged): share columns nbr_cov_lag1_r<rho>[...] / own_cov_lag1, the
    # <panel>.contagion_coverage.csv sidecar and contagion_* stamp fields.
    # Overwrite guard: an existing panel whose prep stamp differs in
    # anything but stamped_at is refused unless --force.
    make_option(
      "--force",
      default = FALSE,
      action = "store_true",
      help = "overwrite an existing panel built under a different prep stamp [default: %default]"
    ),
    make_option(
      "--contagion_coverage_report",
      default = FALSE,
      action = "store_true",
      help = "report the surveillance coverage of the neighbour-contagion term (columns + sidecar + stamp; no model change) [default: %default]"
    ),
    make_option(
      "--adm1_overlap",
      default = "drop",
      type = "character",
      help = "drop | censor_orphan: censor the ADM2 children's rows inside an orphan ADM1 parent episode (\"_a1c\") [default: %default]"
    ),
    # --- ENSO leads -----------------------------------------------------------
    # Also write enso_lead1..enso_lead<k> (ONI k months after t_start;
    # utils/enso.R::join_enso_lags()) for a single-lag sweep with negative lags.
    # Adds columns only: same rows, same filename; recorded in the policy stamp.
    make_option(
      "--enso_max_lead",
      default = 0L,
      type = "integer",
      help = "ENSO lead columns enso_lead1..k to write (0 = none) [default: %default]"
    )
  ),
  prep_cli_options("occ")
)

opt <- parse_args(OptionParser(option_list = option_list))
opt <- validate_prep_cli(opt, "occ")

use_xwalk <- identical(opt$pop_source, "crosswalk")
gid_season <- identical(opt$episode_key, "gid_season")

# --national_zero_fill (gid_season only) adds no-onset rows for the
# national-zero months (01_preprocess_national_zero_months.R) of each gid's
# country; output gets the "_nz" token. --national_zero_rule <name> reads the
# months of a named rule of analysis/national_zero_rules.yml (01 --rule) and
# writes the "_nz<N>-<name>" panel; the fill then also skips months inside
# the gid's own episode spans, and every zero-fill row is asserted outside
# them.
# --at_risk_start after_episode / --adm1_overlap censor_orphan (gid_season
# only) mask rows inside an earlier season's running episode ("_ae") and
# censor the ADM2 children inside an orphan ADM1 parent episode ("_a1c").
# Always written (gid_season): the sidecars <panel>.adm1_overlap.csv,
# .discarded_onsets.csv, .adm0_dropped.csv (and .masked_rows.csv when a mask
# is on) with their n_* / cases_* stamp counts (table S1; per country D31).
if (isTRUE(opt$national_zero_fill) && !gid_season) {
  stop("--national_zero_fill requires --episode_key gid_season.")
}
# Read by the path builders (utils/paths.R::national_zero_fill_of() /
# national_zero_threshold_of()).
opt$episodes$national_zero_fill <- isTRUE(opt$national_zero_fill)
if (is.na(opt$national_zero_threshold) || opt$national_zero_threshold < 0L) {
  stop("--national_zero_threshold must be >= 0.")
}
if (opt$national_zero_threshold > 0L && !opt$episodes$national_zero_fill) {
  stop("--national_zero_threshold needs --national_zero_fill.")
}
opt$episodes$national_zero_threshold <- as.integer(opt$national_zero_threshold)
opt$episodes$national_zero_rule <- .national_zero_rule_check(
  opt$national_zero_rule
)
if (
  nzchar(opt$episodes$national_zero_rule) && !opt$episodes$national_zero_fill
) {
  stop("--national_zero_rule needs --national_zero_fill.")
}
NZ_RULE <- opt$episodes$national_zero_rule
# Item 11 panel fixes: refused without gid_season by the spec itself.
opt$episodes$at_risk_start <- opt$at_risk_start
opt$episodes$adm1_overlap <- opt$adm1_overlap
OCC_PANEL <- occ_panel_spec_of(opt) # validates the whole occurrence-panel spec (utils/paths.R)
NZ_PATH <- if (nzchar(opt$national_zero_path)) {
  opt$national_zero_path
} else {
  as.character(national_zero_months_file_of(OCC_PANEL))
}
if (opt$episodes$national_zero_fill && !file.exists(NZ_PATH)) {
  stop(
    "--national_zero_fill: national-zero months not found: ",
    NZ_PATH,
    ". Build them with Rscript analysis/01_preprocess_national_zero_months.R",
    if (opt$episodes$national_zero_threshold > 0L) {
      paste0(" --zero_threshold ", opt$episodes$national_zero_threshold)
    } else {
      ""
    },
    if (nzchar(NZ_RULE)) paste0(" --rule ", NZ_RULE) else ""
  )
}
# The months file must be the requested rule's: its 01 stamp names the rule (a
# stamp without one was built on the legacy rule) and, under a named rule,
# carries the digest of the registry parameters it was built with.
NZ_RULE_PARAMS <- NULL
if (opt$episodes$national_zero_fill) {
  nz_stamp <- read_prep_policy_stamp(NZ_PATH) %||% list()
  have_rule <- .national_zero_rule_check(nz_stamp[["rule"]] %||% "")
  if (!identical(have_rule, NZ_RULE)) {
    stop(
      "--national_zero_rule '",
      NZ_RULE,
      "' but ",
      NZ_PATH,
      " was built under rule '",
      have_rule,
      "' (01 stamp). Rebuild it with 01_preprocess_national_zero_months.R --rule ",
      if (nzchar(NZ_RULE)) NZ_RULE else "''",
      "."
    )
  }
  if (nzchar(NZ_RULE)) {
    NZ_RULE_PARAMS <- national_zero_rule_params(NZ_RULE)
    if (
      !identical(
        nz_stamp[["rule_params_digest"]],
        national_zero_rule_digest(NZ_RULE_PARAMS)
      )
    ) {
      stop(
        "national-zero rule '",
        NZ_RULE,
        "': the parameters ",
        NZ_PATH,
        " was built with (",
        nz_stamp[["rule_params_digest"]] %||% "none",
        ") differ from the registry's (",
        national_zero_rule_digest(NZ_RULE_PARAMS),
        "). Rebuild it with 01 --rule ",
        NZ_RULE,
        "."
      )
    }
  }
}
nz_report <- list()
panel_report <- list() # panel-mask counts (gid_season only; empty otherwise)

if (!opt$spillover_grouping %in% c("t_end", "calendar")) {
  stop(
    "--spillover_grouping must be 't_end' or 'calendar', got: ",
    opt$spillover_grouping
  )
}

CONTAGION_RHOS <- as.integer(trimws(strsplit(opt$contagion_rhos, ",")[[1]]))
if (
  !length(CONTAGION_RHOS) || anyNA(CONTAGION_RHOS) || any(CONTAGION_RHOS <= 0)
) {
  stop(
    "--contagion_rhos must be a comma-separated list of positive integers (km), got: ",
    opt$contagion_rhos
  )
}

# Deliberately a message, not an error: a legacy invocation carries a
# composite_policy that has nothing to act on.
if (!use_xwalk && !identical(opt$composite_policy, "keep")) {
  message(
    "[occ/composite] pop_source=legacy: composite_policy=",
    opt$composite_policy,
    " is inert (no crosswalk expansion to reduce)"
  )
}

print(opt)

plan(multisession, workers = 3)

# 3. Data ----------------------------------------------------------------------

# Data on outbreaks
ob_ts <- readRDS(make_outbreak_metadata_file()) %>%
  ungroup() %>%
  # Drop location periods whose WorldPop population is implausibly small
  # (broken raster extraction, e.g. pop ~ 5e-6). pop is not used as a model
  # covariate or offset in the occurrence track (log_offset is time-based),
  # but filtering here ensures such locations don't enter the survival dataset.
  # Locations with unknown (NA) population are kept. Mirrors the guard in
  # read_size_data() (utils.R). Threshold set by --min_pop (default 10).
  dplyr::filter(is.na(pop) | pop >= opt$min_pop)

# gid_season: TZA mainland is ADM0 by geometry and is treated as ADM0
# throughout -- marked here, before the section-5 ADM0 filter, so it never
# enters the canonical episode table (03_prepare_siz.R does the same).
if (gid_season) {
  ob_ts <- mark_tza_mainland_adm0(ob_ts)
}

# Coverage data for all locations
all_ts_dat <- readRDS(make_agg_ts_file()) %>%
  filter(location_period_id %in% ob_ts$location_period_id)

# Seasonality estimates
seas <- read_csv("data/best_seasonality_estimates.csv") %>%
  mutate(gid = map_chr(gid, flatten_gid))

# 4. Match GADM to seasonality ids ---------------------------------------------

# Define ADM1 as first ADM2 if model run at adm2.
# Built from `seas` only — always available regardless of opt$redo.
seas_adm1 <- seas %>%
  distinct(gid, run_level) %>%
  mutate(gid_adm1 = map_chr(gid, ~ get_adm1_gid(.))) %>%
  group_by(gid_adm1) %>%
  slice(1) %>%
  ungroup()

# Seasonality-match traceability counts
# NA when a cached
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

  # 4a. Spatial nearest-neighbour fallback for tile-level misses --------------
  # Third fallback, after the exact match and the two numeric get_next_gid() heuristics
  # above both fail. Recovers "tile-level misses": the outbreak's COUNTRY does have
  # >=1 gid in the 24-country seasonality panel, but the specific ADM1/ADM2 tile has
  # none and no in-range numeric neighbour either. Restricted to same-country
  # candidates, so an outbreak in a country the panel does not cover at all (the large
  # majority of the stage-5 exclusion -- see notebooks/descriptive_occurrence_size.Rmd's
  # cascade-seas-breakdown chunk) still falls through to the January default in
  # Section 5, unchanged.
  #
  # Two-stage, mirroring the numeric cascade above:
  #   3a. own-level:  target's `gid`      vs `seas$gid`           candidates
  #   3b. ADM1 level: target's `gid_adm1` vs `seas_adm1$gid_adm1` candidates
  # so `adm1_match` keeps meaning what Section 5's join keys on.
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
    "[occ/seas] spatial nearest-neighbour fallback: %d/%d unmatched gid(s) recovered (%d still unmatched -- country outside the seasonality panel)",
    n_seasgid_spatial_recovered,
    n_seasgid_spatial_checked,
    n_seasgid_spatial_checked - n_seasgid_spatial_recovered
  ))

  # Save to reload.
  saveRDS(
    matched_gadm,
    file = "generated_data/matched_ob_gids_seasonality_occ.rds"
  )
} else {
  seas_cache <- "generated_data/matched_ob_gids_seasonality_occ.rds"
  if (!file.exists(seas_cache)) {
    # The pre-split artefact, accepted once so an existing tree still resumes.
    legacy <- "generated_data/matched_ob_gids_seasonality.rds"
    if (!file.exists(legacy)) {
      stop("No seasonality match cache. Re-run with --redo TRUE.")
    }
    warning(
      "Reading the shared pre-split cache ",
      legacy,
      "; re-run with --redo TRUE to write the occurrence-track copy."
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

# National (ADM0) rows are dropped here, before any panel is built; counted per
# country for the stamp / sidecar.
adm0_tab <- adm0_episode_counts(ob_ts)

# Add the seasonal cutoff month to data
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
# dropping the episode. This covers countries for which seasonality estimates
# were not produced because base was the best-fitting model.
n_seas_episodes_matched <- sum(!is.na(ob_ts2$min_seas))
n_seas_episodes_defaulted_jan <- sum(is.na(ob_ts2$min_seas))
# Each LP's calendar SOURCE, recorded before the default fills the gaps (a
# matched estimate can itself be January): read back per gid at the save
# (season_calendar.R::gid_min_seas_source()) -> the .min_seas_source.csv
# sidecar, so per-country season ends can flag calendar-default countries. A
# separate table, so ob_ts2's columns (and the episode digest) are unchanged.
.MIN_SEAS_SOURCES <- c("estimate_gid", "estimate_adm1", "default_january")
lp_seas_source <- ob_ts2 %>%
  transmute(
    location_period_id = as.character(location_period_id),
    min_seas_source = dplyr::case_when(
      is.na(min_seas) ~ "default_january",
      adm1_match %in% TRUE ~ "estimate_adm1",
      TRUE ~ "estimate_gid"
    )
  ) %>%
  group_by(location_period_id) %>%
  summarise(
    min_seas_source = .MIN_SEAS_SOURCES[min(match(
      min_seas_source,
      .MIN_SEAS_SOURCES
    ))],
    .groups = "drop"
  )
ob_ts2 <- ob_ts2 %>%
  mutate(min_seas = dplyr::coalesce(min_seas, 1L))

message(sprintf(
  "[occ/seas] min_seas: %d episode(s) matched a seasonality estimate, %d defaulted to January (no gid match)",
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

# Some location periods have more than one outbreak per year, noting here
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
ob_ts3_first <- ob_ts3 # discarded-onset table (gid_season)

# 6a. Propagate composite location periods to their member gids ----------------
#
# THE OCCURRENCE CONTRACT DIFFERS FROM THE SIZE CONTRACT, deliberately.
#
# The size track apportions cases by w = H/P_l, because cases are an extensive
# quantity that must be conserved. The outbreak indicator is NOT. It is a 0/1
# event, and w-scaling it produces a number no Bernoulli or survival likelihood
# can consume. So it is PROPAGATED UNCHANGED: every member gid of a composite LP
# records the outbreak.
#
# This is not a convenience — it is the only expansion consistent with the
# operator already in this file. The collapse below reads
#
#     group_by(gid, t_start, t_end) %>% summarise(outbreak = any(outbreak))
#
# i.e. an observation over a footprint IS an any() over its members. A composite
# LP observation is exactly that any() at the LP scale, so unfolding it back to
# its members and letting the same any() re-collapse them is an identity, not an
# assumption.
#
# SEASONALITY. As in the size track, min_seas / year_seas are inherited from the
# centroid gid and not re-resolved per member.
#
# --composite_policy SELECTS WHICH OF THE ABOVE ACTUALLY RUNS. Everything from
# here down sees the REDUCED set: reducing ob_ts3 removes the numerator and, via
# the at-risk padding at section 7a (which inner_joins ob_ts3's distinct gids),
# the matching denominator too. That is what makes `drop` a coverage caveat
# rather than a hazard bias -- but it also means dropped units stop feeding the
# spillover lags at section 8. Documented, not fixed: it is the direct cost of
# reducing early.
if (use_xwalk) {
  xwalk_path <- if (nzchar(opt$crosswalk_path)) {
    opt$crosswalk_path
  } else {
    make_pop_run_file(list(), "pop_crosswalk")
  }

  xwalk_occ <- load_population_crosswalk(xwalk_path)

  # --- explicit drop of crosswalk-unmatched LPs -------------------------------
  # An LP the crosswalk never saw has NO population decomposition, so neither
  # reduction path can act on it: lp_composite_status() refuses to call it
  # composite (no basis for the claim), and both collapse_composite_lps() and
  # apportion_composite_cases() would pass it through untouched under
  # `keep_unmatched = TRUE`, tagged pop_tier = "unmatched_lp". Such a row then
  # reaches the model frame carrying its centroid gid and no crosswalk
  # population -- i.e. exactly the un-reduced composite that --composite_policy
  # exists to remove, smuggled past the policy by an absent measurement.
  #
  # Dropping here is a no-op on the final row set under the default
  # `drop_collapse`: the LPs that miss the crosswalk are the multi-ADM1 unions
  # (a 6-ADM1 composite has gid_cover ~ 0.17, far below --cover_threshold 0.90),
  # so the policy would drop them anyway once they were matchable. The point is
  # that the drop is stated rather than contingent on a match succeeding.
  #
  # This is CONTAINMENT, NOT A FIX. The absent crosswalk rows are an upstream
  # data defect (duplicate/parent-inherited LP geometry in the Taxonomy API);
  # the defect is deliberately left in place so the next data refresh surfaces
  # it again through the count logged below, instead of being silently absorbed.
  #
  # Placed BEFORE the counters below so every gate downstream keeps comparing
  # like with like -- n_before, sig_before and the case-conservation accounting
  # all describe the matchable corpus, and none of them has to special-case a
  # row set the crosswalk could not describe.
  st_match <- lp_composite_status(
    ob_ts3,
    xwalk_occ,
    lp_col = "location_period_id",
    year_col = "year_seas",
    level_col = "admin_level",
    tiling = opt$composite_tiling,
    cover_threshold = opt$cover_threshold
  )
  # One row per input row, arranged by .row (composite_lp.R:241), so this is a
  # positional mask on ob_ts3.
  unmatched <- !st_match$matched
  unmatched_lp <- sort(unique(as.character(ob_ts3$location_period_id[
    unmatched
  ])))
  message(sprintf(
    "[occ/xwalk] crosswalk-unmatched: %d row(s) in %d LP(s) dropped%s",
    sum(unmatched),
    length(unmatched_lp),
    if (length(unmatched_lp)) {
      sprintf(" (e.g. %s)", paste(head(unmatched_lp, 5), collapse = ", "))
    } else {
      ""
    }
  ))
  ob_ts3 <- ob_ts3[!unmatched, , drop = FALSE]
  if (!nrow(ob_ts3)) {
    stop(
      "[occ/xwalk] no location period matched the crosswalk — every row was ",
      "dropped. The crosswalk at '",
      xwalk_path,
      "' is empty, keyed on a ",
      "different gid convention, or built for a different corpus."
    )
  }

  n_before <- nrow(ob_ts3)
  n_gid_before <- dplyr::n_distinct(ob_ts3$gid)
  n_lp_before <- dplyr::n_distinct(ob_ts3$location_period_id)
  # The per-(LP, season) event signature, as it stood before propagation.
  sig_before <- ob_ts3 %>%
    distinct(location_period_id, year_seas, suspected_cases) %>%
    arrange(location_period_id, year_seas)

  do_drop <- opt$composite_policy %in% c("drop", "drop_collapse")
  do_collapse <- opt$composite_policy %in% c("collapse", "drop_collapse")

  # `drop` (without collapse) is the only branch that needs the classification
  # separately: it filters first and then hands the survivors to the unchanged
  # expansion path. `drop_collapse` gets both in one call.
  dropped_lp <- character(0)
  n_retained <- n_before
  if (do_drop && !do_collapse) {
    st <- lp_composite_status(
      ob_ts3,
      xwalk_occ,
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
    n_retained <- nrow(ob_ts3)
  }

  if (do_collapse) {
    ob_ts3 <- collapse_composite_lps(
      ob_ts3,
      xwalk_occ,
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
  } else {
    ob_ts3 <- apportion_composite_cases(
      ob_ts3,
      xwalk_occ,
      lp_col = "location_period_id",
      # NA is the switch that makes this propagation rather than apportionment:
      # no column is multiplied by w. The w-share is still recorded as
      # `apportion_w` for the sensitivity run, but nothing consumes it here.
      cases_col = NA,
      year_col = "year_seas",
      gid_col = "gid",
      level_col = "admin_level",
      tiling = opt$composite_tiling,
      keep_unmatched = TRUE
    )
    rep <- attr(ob_ts3, "apportion_report")
    # The `keep` counterfactual is only knowable for free on the collapse branch
    # (collapse_composite_lps() counts the members it declined to emit). Under
    # `keep` the output IS the counterfactual; under `drop` it is not, and we do
    # not pay for a second match to find out.
    n_if_keep <- if (do_drop) NA_integer_ else rep$n_rows_out
    n_lp_comp <- if (do_drop) length(dropped_lp) else NA_integer_
    max_mem <- rep$max_members
  }

  # --- reporting: quantify the cost without needing a second run --------------
  message(sprintf(
    "[occ/composite] policy=%s threshold=%.2f",
    opt$composite_policy,
    opt$cover_threshold
  ))
  message(sprintf(
    "[occ/composite] LPs: %d total, %s true composite, %d dropped",
    n_lp_before,
    if (is.na(n_lp_comp)) {
      "n/a (not classified under this policy)"
    } else {
      sprintf("%d (%.1f%%)", n_lp_comp, 100 * n_lp_comp / max(n_lp_before, 1L))
    },
    length(dropped_lp)
  ))
  message(sprintf(
    "[occ/composite] rows: %d in -> %d out%s",
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
    "[occ/composite] gids: %d -> %d distinct; max members = %d",
    n_gid_before,
    dplyr::n_distinct(ob_ts3$gid),
    max_mem
  ))
  message(sprintf(
    "[occ/propagate] %d LP(s) matched, %d unmatched",
    rep$n_lp_matched,
    rep$n_lp_unmatched
  ))

  # --- gates, restated per policy ---------------------------------------------
  # Under `keep` the gate is that no event may be LOST, so the row count must not
  # shrink -- a count that did not grow means the crosswalk matched nothing.
  # Under the reducing policies the invariant INVERTS: rows must never exceed the
  # input, and every retained LP-event must contribute exactly one row. Both are
  # the same underlying claim (the reduction did what it said).
  if (!do_collapse) {
    # `n_retained` is `n_before` under `keep`, and the post-drop count under
    # `drop` — so the expansion is judged against what it was actually given,
    # never against rows a policy deliberately removed upstream of it.
    if (nrow(ob_ts3) < n_retained) {
      stop(
        "[occ/propagate] row count SHRANK (",
        n_retained,
        " -> ",
        nrow(ob_ts3),
        "): propagation must never drop an at-risk unit."
      )
    }
  } else {
    if (nrow(ob_ts3) > n_before) {
      stop(
        "[occ/composite] collapse GREW the row count (",
        n_before,
        " -> ",
        nrow(ob_ts3),
        "): collapse must emit at most one row per input row."
      )
    }
    dup <- ob_ts3 %>% count(location_period_id, year_seas) %>% filter(n > 1L)
    if (nrow(dup)) {
      stop(
        "[occ/composite] ",
        nrow(dup),
        " (LP, season) pair(s) contribute more ",
        "than one row after collapse — the one-row-per-LP-event invariant broke."
      )
    }
  }

  # PROPAGATION, NOT APPORTIONMENT. Every retained row must still carry the LP's
  # own value, un-multiplied. If cases_col were ever wired up here by mistake,
  # each member would carry w * cases and this signature would not survive. Under
  # a dropping policy the comparison is restricted to the RETAINED LPs, so the
  # equality assertion itself survives intact and still catches w-scaling.
  sig_after <- ob_ts3 %>%
    distinct(location_period_id, year_seas, suspected_cases) %>%
    arrange(location_period_id, year_seas)
  sig_ref <- sig_before %>%
    filter(location_period_id %in% unique(sig_after$location_period_id)) %>%
    arrange(location_period_id, year_seas)
  if (
    !isTRUE(all.equal(
      as.data.frame(sig_ref),
      as.data.frame(sig_after),
      check.attributes = FALSE
    ))
  ) {
    stop(
      "[occ/propagate] the per-(LP, season) event signature changed. ",
      "Something scaled the response — the occurrence track must carry ",
      "the indicator unchanged, never weight it by w."
    )
  }
  message(
    "[occ/propagate] event signature unchanged on retained LPs: ",
    "indicator carried, not scaled"
  )
}

# 7. Pad no-outbreak years and build the survival panel ------------------------
ep_digest <- NA_character_
ep_report <- list()

if (!gid_season) {
  # 7a. Pad no-outbreak years -----------------------------------------------

  # Make object with time
  ob_ts_zero_padded <- ob_ts3 %>%
    filter(!is.na(min_seas)) %>%
    mutate(time = yday(TL_seas)) %>%
    group_by(location_period_id) %>%
    group_split() %>%
    future_map_dfr(~ add_no_outbreak_years(., all_ts_dat = all_ts_dat)) %>%
    mutate(uid = str_c(location_period_id, "-1")) %>%
    inner_join(
      ob_ts3 %>%
        distinct(
          country,
          location,
          location_period_id,
          gid,
          AFRO_region,
          min_seas,
          admin_level
        )
    )

  # 7b. Make survival format data -------------------------------------------

  # Survival-format data for 0-padded years
  surv_dat_zero_padded <- ob_ts_zero_padded %>%
    filter(!is.na(min_seas)) %>%
    mutate(TL_seas = TR_seas, time = yday(TL_seas)) %>%
    group_by(uid, location_period_id, gid, year_seas) %>%
    group_split() %>%
    future_map_dfr(~ make_surv_data(., baseline = opt$baseline)) %>%
    # !! set all outbreaks to null
    mutate(outbreak = FALSE)

  # Make object with time
  surv_dat <- ob_ts3 %>%
    filter(!is.na(min_seas)) %>%
    mutate(time = yday(TL_seas)) %>%
    group_by(uid, location_period_id, gid, year_seas) %>%
    group_split() %>%
    future_map_dfr(~ make_surv_data(., baseline = opt$baseline))

  # TODO check dupplicated rows with uids and admin levels
  surv_dat_all <- bind_rows(
    surv_dat,
    surv_dat_zero_padded
  ) %>%
    # Add information on country
    inner_join(
      bind_rows(
        ob_ts %>% ungroup() %>% distinct(uid, country, admin_level),
        ob_ts_zero_padded %>% ungroup() %>% distinct(uid, country, admin_level)
      )
    ) %>%
    # Add AFRO region
    taxdat::get_AFRO_region("country") %>%
    fill_synthetic_afro_region() %>% # utils/gid.R: the synthetic countries (no-op for real data)
    arrange(location_period_id, t_start)

  # Aggregate to the simulator's emission granularity: ONE union observation per
  # (gid, month). The occurrence model and the forward simulator must be
  # calibrated at the same spatial resolution. Because ADM3 reporting units are
  # relabelled to their parent ADM2 GID_2 upstream (02_preprocess_climate.R:77-105),
  # a single ADM2 gid footprint can carry many location-period rows in the same
  # month (the ADM2 report plus every nested ADM3 sub-district), all sharing
  # identical gid-level climate covariates. Left unaggregated, the fit sees
  # ~5-10 near-duplicate rows per (gid, month) and calibrates the hazard at
  # subdistrict resolution, whereas the simulator emits ONE Bernoulli per
  # (gid, month) -> systematic occurrence mismatch (see seam-gap investigation).
  # Collapse to one row per (gid, t_start, t_end): outbreak = ANY onset anywhere
  # in the footprint (union); every covariate is constant across the siblings so
  # first() is loss-free. A single representative uid/location_period_id per gid
  # is fixed BEFORE the collapse so downstream season-aware climate lags (grouped
  # by uid) and the time_to_last_outbreak recursion (grouped by
  # location_period_id) both operate per footprint.
  surv_dat_all <- surv_dat_all %>%
    group_by(gid) %>%
    mutate(
      .rep_i = order(location_period_id)[1],
      location_period_id = location_period_id[.rep_i],
      uid = uid[.rep_i]
    ) %>%
    ungroup() %>%
    select(-.rep_i) %>%
    group_by(gid, t_start, t_end) %>%
    summarise(
      outbreak = any(outbreak),
      across(-outbreak, dplyr::first),
      .groups = "drop"
    ) %>%
    arrange(location_period_id, t_start)
} else {
  # D'/E'. Gid-season survival histories (--episode_key gid_season) ---------
  #
  # ONE survival history per (gid, season-year), built from the canonical
  # episode table that 03_prepare_siz.R builds with the identical call
  # (utils/footprint.R::build_canonical_episodes()). Compared with the lp_season
  # branch above, which builds a history PER LOCATION PERIOD, censors each at
  # its own onset and only then unions the LPs per (gid, t_start, t_end):
  #   * one onset per gid-season, at the EARLIEST member onset -- a second LP's
  #     later onset inside a running outbreak is not a second event;
  #   * the gid is censored after that onset -- no "at risk, no onset" months
  #     contributed by sibling LPs while the outbreak is running;
  #   * a covered season without an episode is at risk until the UNION of the
  #     member LPs' coverage ends;
  #   * no union step, so no duplicate t_start rows per gid (the incomplete
  #     union the simulator's build_obs_panel() had to repair with a max);
  #   * admin_level is the gid's GEOMETRY-DERIVED level, one per gid; the
  #     corpus label is kept only as admin_level_corpus.
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
  ep_occ <- ep_res$episodes
  message(sprintf(
    paste0(
      "[occ/episodes] gid_season: %d member row(s) -> %d canonical episode(s) ",
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

  # At-risk window: each member LP's surveillance coverage, read on its GID's
  # canonical season calendar, unioned per (gid, season-year).
  at_risk <- gid_at_risk_years(
    ep_res$member_rows %>% distinct(location_period_id, gid, min_seas),
    all_ts_dat
  )

  surv_dat_all <- gid_season_surv_input(ep_occ, at_risk) %>%
    group_by(gid, year_seas) %>%
    group_split() %>%
    future_map_dfr(function(x) {
      make_surv_data(x, baseline = opt$baseline) %>%
        # covered-but-no-episode seasons: no event, exactly as the lp_season
        # branch zeroes its padded years
        mutate(outbreak = outbreak & x$.has_onset[1])
    })

  # National zero-fill (--national_zero_fill): no-onset rows for the
  # national-zero months of each gid's country not yet in the panel.
  if (opt$episodes$national_zero_fill) {
    nz_rows <- national_zero_surv_rows(
      surv_dat_all,
      readRDS(NZ_PATH),
      ids = ep_occ %>%
        distinct(gid, uid, location_period_id, min_seas, country),
      baseline = opt$baseline,
      map_fn = future_map_dfr,
      # a named rule also keeps the fill out of the gid's own episode spans
      episodes = if (nzchar(NZ_RULE)) ep_occ else NULL
    )
    nz_report <- list(
      national_zero_fill = TRUE,
      national_zero_threshold = opt$episodes$national_zero_threshold,
      national_zero_path = NZ_PATH,
      national_zero_mtime = format(file.mtime(NZ_PATH)),
      n_nz_candidates = attr(nz_rows, "n_candidates"),
      n_nz_excluded_after_onset = attr(nz_rows, "n_after_onset"),
      n_nz_season_years = attr(nz_rows, "n_season_years"),
      n_nz_new_season_years = attr(nz_rows, "n_new_season_years"),
      n_nz_rows = nrow(nz_rows)
    )
    # Named rule only (empty otherwise, so the legacy stamp is unchanged).
    if (nzchar(NZ_RULE)) {
      nz_report <- c(
        nz_report,
        list(
          national_zero_rule = NZ_RULE,
          national_zero_params = paste(
            sprintf(
              "%s=%s",
              setdiff(names(NZ_RULE_PARAMS), "rule"),
              unlist(NZ_RULE_PARAMS[setdiff(names(NZ_RULE_PARAMS), "rule")])
            ),
            collapse = ";"
          ),
          national_zero_params_digest = national_zero_rule_digest(
            NZ_RULE_PARAMS
          ),
          n_nz_excluded_in_episode = attr(nz_rows, "n_in_episode")
        )
      )
    }
    message(sprintf(
      paste0(
        "[occ/national-zero] %d candidate gid-month(s); %d dropped at/after an onset; ",
        "%d row(s) added over %d season-year(s) (%d new)"
      ),
      nz_report$n_nz_candidates,
      nz_report$n_nz_excluded_after_onset,
      nz_report$n_nz_rows,
      nz_report$n_nz_season_years,
      nz_report$n_nz_new_season_years
    ))
    surv_dat_all <- bind_rows(
      mutate(surv_dat_all, national_zero_fill = FALSE),
      nz_rows
    )
  }

  # One country and one (geometry-derived) admin level per gid.
  gid_meta <- ep_occ %>%
    group_by(gid) %>%
    summarise(
      n_country = n_distinct(country),
      country = first(country),
      admin_level = first(admin_level),
      admin_level_corpus = names(sort(
        table(admin_level_corpus),
        decreasing = TRUE
      ))[1],
      .groups = "drop"
    )
  if (any(gid_meta$n_country != 1L)) {
    stop(
      "[occ/episodes] ",
      sum(gid_meta$n_country != 1L),
      " gid(s) map to ",
      "more than one country."
    )
  }

  surv_dat_all <- surv_dat_all %>%
    inner_join(gid_meta %>% select(-n_country), by = "gid") %>%
    taxdat::get_AFRO_region("country") %>%
    fill_synthetic_afro_region() %>% # utils/gid.R: the synthetic countries (no-op for real data)
    arrange(location_period_id, t_start)

  # After the zero-fill: the ADM1-ADM2
  # overlap classification (always; the load-time drop itself is unchanged),
  # the discarded onsets (always), then the masks of the spec -- identical
  # panel at the defaults (season / drop).
  adm1_overlap_tab <- adm1_overlap_parents(surv_dat_all, ep_occ)
  discarded_tab <- discarded_onset_table(
    ob_ts2,
    ob_ts3_first,
    ob_ts3 %>% filter(!is.na(min_seas)),
    ep_res
  )
  ov_before <- .overlap_adm1_set(mutate(
    surv_dat_all,
    ADM1 = gid_to_adm1(gid, "03_prepare_occ")
  ))
  surv_dat_all <- mask_occ_panel_rows(
    surv_dat_all,
    ep_occ,
    OCC_PANEL,
    overlap = adm1_overlap_tab
  )
  mask_attr <- function(a) as.integer(attr(surv_dat_all, a) %||% 0L)
  masked_rows <- attr(surv_dat_all, "masked_rows")
  n_onsets_masked <- mask_attr("n_onsets_masked_running_episode")
  if (
    !identical(
      .overlap_adm1_set(mutate(
        surv_dat_all,
        ADM1 = gid_to_adm1(gid, "03_prepare_occ")
      )),
      ov_before
    )
  ) {
    stop(
      "[occ/panel] the masks changed the set of ADM1s with ADM2 rows, so the fit's load-time ",
      "ADM1-ADM2 overlap drop would change; refusing."
    )
  }
  panel_report <- list(
    n_adm0_episodes_dropped = as.integer(sum(adm0_tab$n_adm0_episodes)),
    cases_adm0_dropped = sum(adm0_tab$cases_adm0),
    n_adm1_overlap_gids = dplyr::n_distinct(adm1_overlap_tab$gid),
    n_adm1_overlap_onsets = nrow(adm1_overlap_tab),
    n_adm1_overlap_orphan = sum(adm1_overlap_tab$class == "orphan"),
    n_adm1_overlap_covered = sum(adm1_overlap_tab$class == "covered"),
    n_onsets_discarded_lp_season = sum(
      discarded_tab$reason == "lp_second_onset"
    ),
    n_onsets_discarded_after_episode_end = sum(
      discarded_tab$after_episode_end %in% TRUE
    ),
    n_onsets_discarded_lp_calendar_only = sum(
      discarded_tab$lp_calendar_only %in% TRUE
    ),
    cases_onsets_discarded_lp_season = sum(
      discarded_tab$cases[discarded_tab$reason == "lp_second_onset"],
      na.rm = TRUE
    ),
    n_members_absorbed = sum(discarded_tab$reason == "absorbed_member"),
    n_members_absorbed_after_gap = sum(
      discarded_tab$reason == "absorbed_after_gap"
    ),
    n_rows_masked_running_episode = mask_attr("n_rows_masked_running_episode"),
    n_onsets_masked_running_episode = n_onsets_masked,
    n_nz_rows_masked = mask_attr("n_nz_rows_masked"),
    n_child_rows_censored_adm1 = mask_attr("n_child_rows_censored_adm1")
  )
  # The spec fields, only when non-default (default stamps keep their field set).
  if (!identical(OCC_PANEL$at_risk_start, "season")) {
    panel_report$at_risk_start <- OCC_PANEL$at_risk_start
  }
  if (!identical(OCC_PANEL$adm1_overlap, "drop")) {
    panel_report$adm1_overlap <- OCC_PANEL$adm1_overlap
  }
  message(sprintf(
    paste0(
      "[occ/panel] ADM1-ADM2 overlap: %d parent episode(s) (%d orphan, %d covered); discarded onsets: ",
      "%d LP second onset(s) (%d after the episode ended, %d on the LP calendar only), %d absorbed member(s) ",
      "(%d after a gap); ADM0: %d row(s) dropped; masks (%s / %s): %d running-episode row(s) (%d onset(s), ",
      "%d zero-fill), %d child row(s) censored"
    ),
    panel_report$n_adm1_overlap_onsets,
    panel_report$n_adm1_overlap_orphan,
    panel_report$n_adm1_overlap_covered,
    panel_report$n_onsets_discarded_lp_season,
    panel_report$n_onsets_discarded_after_episode_end,
    panel_report$n_onsets_discarded_lp_calendar_only,
    panel_report$n_members_absorbed,
    panel_report$n_members_absorbed_after_gap,
    panel_report$n_adm0_episodes_dropped,
    OCC_PANEL$at_risk_start,
    OCC_PANEL$adm1_overlap,
    panel_report$n_rows_masked_running_episode,
    n_onsets_masked,
    panel_report$n_nz_rows_masked,
    panel_report$n_child_rows_censored_adm1
  ))

  # Post-mask gate: no zero-fill row inside
  # its gid's own episode span. Asserted under a named rule, report-only under
  # the legacy rule (its fill does not exclude a previous season's running episode).
  if (opt$episodes$national_zero_fill) {
    assert_nz_rows_outside_episodes(
      surv_dat_all,
      ep_occ,
      assert = nzchar(NZ_RULE)
    )
  }

  # Gates: one row per (gid, t_start); one onset per canonical episode (less
  # the onsets masked inside an earlier season's running episode).
  dup <- surv_dat_all %>% count(gid, t_start) %>% filter(n > 1L)
  if (nrow(dup)) {
    stop(
      "[occ/episodes] ",
      nrow(dup),
      " (gid, t_start) pair(s) carry more than ",
      "one survival row -- the gid-season builder must emit one history per gid."
    )
  }
  n_ev <- sum(as.logical(surv_dat_all$outbreak), na.rm = TRUE)
  if (n_ev != nrow(ep_occ) - n_onsets_masked) {
    stop(sprintf(
      "[occ/episodes] %d onset row(s) but %d canonical episode(s) - %d masked onset(s).",
      n_ev,
      nrow(ep_occ),
      n_onsets_masked
    ))
  }
  message(sprintf(
    "[occ/episodes] %d survival row(s), %d onset(s) = canonical episodes%s",
    nrow(surv_dat_all),
    n_ev,
    if (n_onsets_masked) {
      sprintf(" - %d masked inside a running episode", n_onsets_masked)
    } else {
      ""
    }
  ))
  if (n_onsets_masked) {
    # soft check, per country: onsets + masked onsets = canonical episodes
    pc <- tibble(country = as.character(ep_occ$country)) %>%
      count(country, name = "episodes") %>%
      left_join(
        tibble(
          country = as.character(surv_dat_all$country[
            as.logical(surv_dat_all$outbreak) %in% TRUE
          ])
        ) %>%
          count(country, name = "onsets"),
        by = "country"
      ) %>%
      left_join(
        tibble(
          country = gid_meta$country[match(
            masked_rows$gid[masked_rows$outbreak],
            gid_meta$gid
          )]
        ) %>%
          count(country, name = "masked"),
        by = "country"
      ) %>%
      mutate(across(c(onsets, masked), ~ coalesce(.x, 0L)))
    message(
      "[occ/episodes] per country, onsets + masked = episodes: ",
      paste(
        sprintf(
          "%s %d+%d=%d%s",
          pc$country,
          pc$onsets,
          pc$masked,
          pc$episodes,
          ifelse(pc$onsets + pc$masked == pc$episodes, "", " (!)")
        ),
        collapse = ", "
      )
    )
  }
  if (opt$episodes$national_zero_fill) {
    if (
      any(surv_dat_all$national_zero_fill & as.logical(surv_dat_all$outbreak))
    ) {
      stop("[occ/national-zero] a national zero-fill row carries an onset.")
    }
    nz_tab <- table(surv_dat_all$country[surv_dat_all$national_zero_fill])
    message(
      "[occ/national-zero] rows by country: ",
      paste(
        sprintf("%s=%d", names(nz_tab), as.integer(nz_tab)),
        collapse = ", "
      )
    )
  }
}

# Compute time to last outbreak
surv_dat_all <- surv_dat_all %>%
  group_by(location_period_id) %>%
  mutate(
    rid = row_number(),
    last_outbreak_start = .TTLO_ORIGIN, # base.R: as.Date("2000-01-01")
    last_outbreak_start = case_when(
      outbreak ~ t_end,
      TRUE ~ last_outbreak_start
    ),
    last_outbreak_start = cummax(as.numeric(last_outbreak_start)) %>%
      as.Date(origin = "1970-01-01") %>%
      lag(n = 1),
    time_to_last_outbreak = difftime(
      t_start,
      last_outbreak_start,
      units = "days"
    ) %>%
      as.numeric()
  ) %>%
  ungroup()

# 8. Add simultaneous outbreaks ------------------------------------------------

adm1_outbreak <- surv_dat_all %>%
  mark_tza_mainland_adm0() %>%
  ungroup() %>%
  # flatten_gid() first: under --pop_source crosswalk, gid comes from the LP x gid
  # crosswalk in DOTTED GADM notation ("AGO.10_1"), not the underscore form this
  # regex expects. Unflattened, str_extract() returns NA for every
  # crosswalk-sourced row, and NA-ADM1 rows from every country/time then match
  # each other in the left_joins below (dplyr's default NA == NA) -- a
  # many-to-many cross join across the whole panel that OOM-killed this script
  # (anon-rss ~20-28GB) before this fix. Same convention as get_stats_gid() /
  # get_extremes_gid() (utils.R), which flatten at point of use for the same
  # reason; flatten_gid() is a no-op on already-underscore gids.
  mutate(ADM1 = gid_to_adm1(gid, "03_prepare_occ")) %>%
  filter(admin_level != "ADM0") %>%
  group_by(country, ADM1, t_start, t_end) %>%
  summarise(outbreak = any(outbreak))

# Active at any time in the calendar month before t_start's month (a span
# overlap, not the one-day test at t_start - 30, which skipped February for
# March rows and missed outbreaks inside the previous month). One value per
# group of span rows; NA spans (no outbreak) count as inactive.
# active_last_month() lives in utils/occ_data.R (the contagion-coverage
# activity grid reuses it).

# Activity spans: under gid_season the CANONICAL episode spans (earliest
# member onset to latest member end, one per gid-season), so the spillover and
# contagion covariates see the same outbreaks the response does; under
# lp_season the per-LP rows, as before.
act_src <- if (gid_season) ep_occ else ob_ts3

# Pad with time series.
# Explicit by= to avoid silently picking up extra implicit keys (e.g. `outbreak`)
# if ob_ts3 ever acquires a column whose name overlaps with adm1_outbreak.
# Many-to-many is intentional: we want all ob_ts3 TL/TR values for each ADM1
# so that the summarise below can check any(TL <= t_start & TR >= t_start).
adm1_outbreak <- adm1_outbreak %>%
  ungroup() %>%
  left_join(
    act_src %>%
      mark_tza_mainland_adm0() %>%
      ungroup() %>%
      # flatten_gid() -- see comment above; keeps this side's ADM1 in the same
      # underscore form as the left side so composite/crosswalk rows aren't
      # dropped to NA here too.
      mutate(ADM1 = gid_to_adm1(gid, "03_prepare_occ")),
    by = c("country", "ADM1"),
    relationship = "many-to-many"
  ) %>%
  group_by(country, ADM1, t_start, t_end) %>%
  summarise(
    ongoing_outbreak = any(TL <= t_start & TR >= t_start),
    ongoing_outbreak_lag1 = active_last_month(TL, TR, t_start),
    outbreak = outbreak[1]
  )
#
# adm1_outbreak %>%
#   ggplot(aes(x = t_start, y = ADM1)) +
#   geom_tile(aes(fill = ongoing_outbreak)) +
#   geom_tile(data = adm1_outbreak %>% filter(outbreak), fill = "black") +
#   facet_grid(country ~ ., space = "free", scales = "free") +
#   theme_bw()

hmm_clusters <- read_csv("data/hmm_country_clusters.csv") %>%
  select(iso3, cluster)

# ADM1 x calendar-month activity grid, built outside the `calendar` branch
# because the distance-weighted neighbour contagion needs it on BOTH
# spillover_grouping paths -- it is a pure function of adm1_outbreak and
# references nothing branch-local. `ongoing_outbreak`/`ongoing_outbreak_lag1` only depend on
# t_start (TL/TR tested at t_start / t_start - 30 above), so collapsing t_end
# here loses nothing.
#
# The grid is every ADM1 of a country x every month the country's panel has,
# not only the cells with a survival row: an ADM1 leaves the risk set after its
# onset, and its running outbreak must still count for its neighbours' country /
# TU / nbr_lag1_* covariates in those months. Cells with a survival row get the
# same flags as adm1_outbreak (same spans, same tests).
act_spans <- act_src %>%
  mark_tza_mainland_adm0() %>%
  ungroup() %>%
  mutate(ADM1 = gid_to_adm1(gid, "03_prepare_occ")) %>%
  distinct(country, ADM1, TL, TR)
adm1_active_cal <- adm1_outbreak %>%
  ungroup() %>%
  distinct(country, ADM1) %>%
  inner_join(
    adm1_outbreak %>% ungroup() %>% distinct(country, t_start),
    by = "country",
    relationship = "many-to-many"
  ) %>%
  left_join(
    act_spans,
    by = c("country", "ADM1"),
    relationship = "many-to-many"
  ) %>%
  group_by(country, ADM1, t_start) %>%
  summarise(
    ongoing_outbreak = any(TL <= t_start & TR >= t_start),
    ongoing_outbreak_lag1 = active_last_month(TL, TR, t_start),
    .groups = "drop"
  )

if (identical(opt$spillover_grouping, "calendar")) {
  # --- calendar path (leakage-free) -------------------------------------------
  # Survival intervals are cut at the onset event, so grouping the leave-one-out
  # on (..., t_start, t_end) (the `t_end` path below) isolates onset rows into
  # their own count group and mechanically relabels them o_ct=0/o_tu=0 -- an
  # outcome-leakage artifact. Fix: collapse
  # the ADM1 active grid to one row per (country, ADM1, t_start) FIRST --
  # `ongoing_outbreak`/`ongoing_outbreak_lag1` only depend on t_start (TL/TR
  # tested at t_start / t_start - 30, lines 708-709) -- then do every
  # leave-one-out at that calendar-month resolution, with no t_end anywhere in
  # a group_by or join key.
  own_key <- c("ADM1", "t_start")
  country_key <- c("ADM1", "t_start")
  tu_key <- c("country", "t_start")

  # adm1_active_cal is built above the branch (the contagion kernel uses it on
  # the t_end path too).
  own_tbl <- adm1_active_cal %>%
    select(ADM1, t_start, ongoing_outbreak, ongoing_outbreak_lag1)

  # Step 1 — Country-level spillover (any *other* ADM1 in same country active at t_start - 30)
  country_lag1 <- adm1_active_cal %>%
    group_by(country, t_start) %>%
    mutate(country_count_lag1 = sum(ongoing_outbreak_lag1)) %>%
    ungroup() %>%
    mutate(
      ongoing_outbreak_country_lag1 = (country_count_lag1 -
        as.integer(ongoing_outbreak_lag1)) >
        0
    ) %>%
    select(ADM1, t_start, ongoing_outbreak_country_lag1)

  # Step 2 — TU-level spillover (any *other country* in same HMM transmission unit active at t_start - 30)
  country_active_lag1 <- adm1_active_cal %>%
    mutate(iso3 = str_extract(ADM1, "^[A-Z]{3}")) %>%
    left_join(hmm_clusters, by = "iso3") %>%
    group_by(country, cluster, t_start) %>%
    summarise(
      country_active_lag1 = any(ongoing_outbreak_lag1),
      .groups = "drop"
    )

  tu_lag1 <- country_active_lag1 %>%
    group_by(cluster, t_start) %>%
    mutate(tu_count_lag1 = sum(country_active_lag1)) %>%
    ungroup() %>%
    mutate(
      ongoing_outbreak_tu_lag1 = (tu_count_lag1 -
        as.integer(country_active_lag1)) >
        0
    ) %>%
    select(country, t_start, ongoing_outbreak_tu_lag1)
} else {
  # --- t_end path (legacy, DEFAULT) -- byte-identical to pre-flag behaviour --
  own_key <- c("ADM1", "t_start", "t_end")
  country_key <- c("ADM1", "t_start", "t_end")
  tu_key <- c("country", "t_start", "t_end")

  own_tbl <- adm1_outbreak %>%
    ungroup() %>%
    select(ADM1, t_start, t_end, ongoing_outbreak, ongoing_outbreak_lag1)

  # Step 1 — Country-level spillover (any *other* ADM1 in same country active at t_start - 30)
  country_lag1 <- adm1_outbreak %>%
    ungroup() %>%
    group_by(country, t_start, t_end) %>%
    mutate(country_count_lag1 = sum(coalesce(ongoing_outbreak_lag1, FALSE))) %>%
    ungroup() %>%
    mutate(
      ongoing_outbreak_country_lag1 = (country_count_lag1 -
        as.integer(coalesce(ongoing_outbreak_lag1, FALSE))) >
        0
    ) %>%
    select(ADM1, t_start, t_end, ongoing_outbreak_country_lag1)

  # Step 2 — TU-level spillover (any *other country* in same HMM transmission unit active at t_start - 30)
  country_active_lag1 <- adm1_outbreak %>%
    ungroup() %>%
    mutate(iso3 = str_extract(ADM1, "^[A-Z]{3}")) %>%
    left_join(hmm_clusters, by = "iso3") %>%
    group_by(country, cluster, t_start, t_end) %>%
    summarise(
      country_active_lag1 = any(coalesce(ongoing_outbreak_lag1, FALSE)),
      .groups = "drop"
    )

  tu_lag1 <- country_active_lag1 %>%
    group_by(cluster, t_start, t_end) %>%
    mutate(tu_count_lag1 = sum(country_active_lag1)) %>%
    ungroup() %>%
    mutate(
      ongoing_outbreak_tu_lag1 = (tu_count_lag1 -
        as.integer(country_active_lag1)) >
        0
    ) %>%
    select(country, t_start, t_end, ongoing_outbreak_tu_lag1)
}

# Step 3 — Distance-weighted neighbour contagion (utils/spatial_contagion.R)
#   nbr_rho(i, t) = sum_{j != i} ongoing_outbreak_lag1_j(t) * exp(-d_ij / rho)
# Built at ADM1 grain -- the only grain at which kernel weights are commensurable
# across the mixed ADM0/ADM1/ADM2 rows, and at which self-exclusion (Inf diagonal
# in build_adm1_dist_km()) cleanly separates the neighbour term from the unit's own
# `ongoing_outbreak_lag1`. ADM2 rows inherit their parent ADM1's value.
# Keyed on (ADM1, t_start) on BOTH spillover_grouping paths: these columns depend
# only on t_start, never t_end, which is exactly the leakage-free property the
# calendar path argues for above. Do not mirror own_key/country_key/tu_key.
adm1_dist_km <- build_adm1_dist_km(
  iso3_set = unique(str_extract(adm1_active_cal$ADM1, "^[A-Z]{3}")),
  gadm_path = opt$contagion_gadm_path,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_contagion_dist)
)
CONTAGION_PARTITION_SOURCES <- trimws(strsplit(
  opt$contagion_partition_source,
  ","
)[[1L]])
CONTAGION_PARTITION_SOURCES <- CONTAGION_PARTITION_SOURCES[nzchar(
  CONTAGION_PARTITION_SOURCES
)]
CONTAGION_PARTITION <- if (length(CONTAGION_PARTITION_SOURCES)) {
  load_grouping(CONTAGION_PARTITION_SOURCES[[1L]])
} else {
  NULL
}
CONTAGION_PARTITIONS <- if (length(CONTAGION_PARTITION_SOURCES)) {
  stats::setNames(
    lapply(CONTAGION_PARTITION_SOURCES, load_grouping),
    CONTAGION_PARTITION_SOURCES
  )
}
nbr_tbl <- build_contagion_table(
  adm1_active_cal,
  adm1_dist_km,
  rhos = CONTAGION_RHOS,
  partition = CONTAGION_PARTITION,
  partitions = CONTAGION_PARTITIONS
)
NBR_COLS <- setdiff(names(nbr_tbl), c("ADM1", "t_start"))

# Add information
# !! Explicit `by` and drop `outbreak`/`country` from RHS: prior version did a
# natural join that included `outbreak` as an implicit key, silently producing
# NA `ongoing_outbreak_*` for zero-padded rows whenever ADM1-level
# any(outbreak) was TRUE. The `ungroup()` is needed because `adm1_outbreak` is
# still grouped by country after the summarise — without it, dplyr re-adds
# country into the `select`, which then triggers country.x/country.y suffixes
# in the join. Keep the survival response (`outbreak`) authoritative.
surv_dat_all <- surv_dat_all %>%
  mark_tza_mainland_adm0() %>%
  # flatten_gid() -- see comment at the adm1_outbreak block above. Without it
  # this side's ADM1 is NA for every crosswalk-sourced row too, which combined
  # with an un-flattened right side is exactly the many-to-many blowup this
  # fix removes.
  mutate(ADM1 = gid_to_adm1(gid, "03_prepare_occ")) %>%
  left_join(own_tbl, by = own_key) %>%
  left_join(country_lag1, by = country_key) %>%
  left_join(tu_lag1, by = tu_key) %>%
  left_join(nbr_tbl, by = c("ADM1", "t_start")) %>%
  # Zero-padded rows and ADM1s outside the GADM/outbreak universe must read as
  # "no neighbour pressure" (0), not NA -- brms drops NA-covariate rows silently.
  mutate(across(all_of(NBR_COLS), ~ coalesce(.x, 0)))

# Contagion coverage report (--contagion_coverage_report, off by
# default = no column, no sidecar, no stamp field). Coverage "reporting": the
# panel's full LP member bridge through coverage.R::build_coverage_mask()
# (include_adm1, the panel's own national-zero months) UNION the recorded
# episode months, read at t-1 like the covariate. Columns only -- the covariate
# itself is untouched (the activity truncation is counted, not fixed).
contagion_report <- list()
COV_SIDECAR_PATH <- NA_character_
if (isTRUE(opt$contagion_coverage_report)) {
  cov_bridge <- (if (gid_season) ep_res$member_rows else ob_ts3) %>%
    ungroup() %>%
    transmute(
      location_period_id = as.character(location_period_id),
      gid = as.character(gid),
      country = as.character(country)
    ) %>%
    filter(!is.na(location_period_id), !is.na(gid)) %>%
    distinct()
  cov_gid <- contagion_coverage_gid(
    cov_bridge,
    make_agg_ts_file(),
    make_agg_adm_ts_file(),
    act_src,
    nz_on = opt$episodes$national_zero_fill,
    nz_path = NZ_PATH
  )
  cov_universe <- rownames(adm1_dist_km)
  cov_months <- sort(unique(surv_dat_all$t_start))
  cov_grid <- adm1_coverage_grid(
    cov_gid,
    tidyr::expand_grid(
      ADM1 = sort(union(cov_universe, unique(surv_dat_all$ADM1))),
      t_start = cov_months
    ),
    lag = 1L
  )
  cov_shares <- contagion_coverage_shares(
    filter(cov_grid, ADM1 %in% cov_universe),
    adm1_dist_km,
    CONTAGION_RHOS,
    CONTAGION_PARTITION,
    CONTAGION_PARTITIONS
  )
  COV_SHARE_COLS <- setdiff(names(cov_shares), c("ADM1", "t_start"))
  if (any(grepl(nbr_covariate_pattern(), c(COV_SHARE_COLS, "own_cov_lag1")))) {
    stop(
      "[occ/contagion-coverage] a share column name matches nbr_covariate_pattern()."
    )
  }
  # counts BEFORE the join (rows the covariate itself coalesced from NA to 0)
  n_na_coalesced <- sum(
    !paste(surv_dat_all$ADM1, surv_dat_all$t_start) %in%
      paste(nbr_tbl$ADM1, nbr_tbl$t_start)
  )
  # M1: neighbour activity at t-1 that the covariate's per-country grid loses
  glob_act <- adm1_activity_grid(act_spans, cov_universe, cov_months)
  lost <- glob_act %>%
    filter(ongoing_outbreak_lag1) %>%
    anti_join(
      adm1_active_cal %>% filter(ongoing_outbreak_lag1 %in% TRUE),
      by = c("ADM1", "t_start")
    )
  surv_dat_all <- surv_dat_all %>%
    left_join(cov_shares, by = c("ADM1", "t_start")) %>%
    left_join(
      cov_grid %>% transmute(ADM1, t_start, own_cov_lag1 = as.integer(covered)),
      by = c("ADM1", "t_start")
    ) %>%
    mutate(across(all_of(c(COV_SHARE_COLS, "own_cov_lag1")), ~ coalesce(.x, 0)))
  # sidecar: ADM1 x month (the t-1 month each panel month reads)
  gid_adm1 <- surv_dat_all %>% distinct(gid = as.character(gid), ADM1)
  n_cov <- cov_gid %>%
    inner_join(gid_adm1, by = "gid") %>%
    count(ADM1, month, name = "n_gid_covered")
  cov_side <- cov_grid %>%
    mutate(month = lubridate::floor_date(t_start, "month") %m-% months(1L)) %>%
    left_join(count(gid_adm1, ADM1, name = "n_gid_panel"), by = "ADM1") %>%
    left_join(n_cov, by = c("ADM1", "month")) %>%
    mutate(
      n_gid_panel = coalesce(n_gid_panel, 0L),
      n_gid_covered = coalesce(n_gid_covered, 0L),
      in_universe = ADM1 %in% cov_universe
    ) %>%
    select(
      ADM1,
      t_start,
      month,
      covered,
      n_gid_panel,
      n_gid_covered,
      in_universe
    )
  q <- function(x) {
    sprintf(
      "%.4g/%.4g/%.4g/%.4g",
      stats::quantile(x, 0.1),
      stats::quantile(x, 0.5),
      stats::quantile(x, 0.9),
      mean(x)
    )
  }
  contagion_report <- list(
    contagion_coverage_report = TRUE,
    contagion_coverage_def = "reporting: build_coverage_mask(LP member bridge, include_adm1) U episode months, at t-1",
    contagion_coverage_nz_path = if (opt$episodes$national_zero_fill) {
      NZ_PATH
    } else {
      ""
    },
    contagion_universe_n = length(cov_universe),
    contagion_universe_digest = content_digest(sort(cov_universe)),
    contagion_coverage_digest = content_digest(as.data.frame(cov_side)),
    contagion_share_cols = paste(COV_SHARE_COLS, collapse = ","),
    contagion_cov_share_q = paste(
      sprintf(
        "%s=%s",
        COV_SHARE_COLS,
        vapply(COV_SHARE_COLS, function(cc) q(surv_dat_all[[cc]]), character(1))
      ),
      collapse = ";"
    ),
    contagion_rows_zero_cover = paste(
      sprintf(
        "%s=%d",
        COV_SHARE_COLS,
        vapply(
          COV_SHARE_COLS,
          function(cc) sum(surv_dat_all[[cc]] == 0),
          integer(1)
        )
      ),
      collapse = ";"
    ),
    contagion_activity_cells_lost = nrow(lost),
    contagion_rows_na_coalesced = n_na_coalesced,
    contagion_own_uncovered_rows = sum(surv_dat_all$own_cov_lag1 == 0)
  )
  message(sprintf(
    paste0(
      "[occ/contagion-coverage] universe %d ADM1 (digest %s); %d share column(s); ",
      "own ADM1 uncovered at t-1 in %d row(s); %d neighbour-activity cell(s) lost by ",
      "the per-country grid (M1, not fixed); %d row(s) coalesced NA -> 0"
    ),
    length(cov_universe),
    contagion_report$contagion_universe_digest,
    length(COV_SHARE_COLS),
    contagion_report$contagion_own_uncovered_rows,
    nrow(lost),
    n_na_coalesced
  ))
}

# 9. Climate and ENSO lag columns ----------------------------------------------
# Add climate lag columns (lag1..lag11) for DLNM use at fit time.
# Each lag is computed within the uid × season-year group so season boundaries
# produce NAs (→ 0) rather than spurious carry-overs.
surv_dat_all <- make_climate_lags(surv_dat_all, max_lag = 11)

# Add ENSO lag columns (enso_lag0 ... enso_lag12).
# ONI (or raw Niño 3.4) is a global index joined by calendar month of t_start.
# All 13 lag columns are always written so that lag-sensitivity configs can
# select any lag without re-running prep.
enso_dat <- fetch_enso_index(
  source = opt$enso_source,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_enso)
)
surv_dat_all <- join_enso_lags(
  surv_dat_all,
  enso_dat,
  date_col = "t_start",
  max_lead = opt$enso_max_lead
)

# 10. Guard against missing climate data ---------------------------------------
# left_join in make_surv_data() (utils.R) preserves outbreak rows beyond the
# CHIRPS/ERA5 coverage window (currently ends 2024-12); flag them here so
# 04_fit_* scripts don't silently pass NA covariates into brms.
surv_dat_all <- guard_missing_climate_data(surv_dat_all)

# 11. Save ---------------------------------------------------------------------
# Outputs: generated_data/outbreaks_survival_data_base-<baseline>.rds
#   (make_surv_dat_file()); generated_data/matched_ob_gids_seasonality_occ.rds;
#   <survival panel stem>.min_seas_source.csv (each gid's season calendar and
#   its source: seasonality estimate at the gid / its ADM1, or the January
#   default; counts in the policy stamp).
# Everything this section writes is computed FIRST, the policy stamp
# is built, and an existing panel whose stamp differs in anything but
# stamped_at is refused unless --force -- before any file is written.

SURV_DAT_PATH <- make_surv_dat_file(opt = opt)

# Per-gid season-calendar source sidecar: each
# panel gid's canonical min_seas (modal over its member rows, ties -> smallest
# month, as footprint.R::canonical_season_calendar()) and whether it is a
# seasonality estimate (at the gid / its ADM1) or the January default.
MIN_SEAS_SOURCE_PATH <- sub("\\.rds$", ".min_seas_source.csv", SURV_DAT_PATH)
gid_seas_source <- gid_min_seas_source(
  ob_ts3 %>%
    filter(!is.na(min_seas)) %>%
    transmute(
      gid = as.character(gid),
      location_period_id = as.character(location_period_id),
      min_seas
    ),
  lp_seas_source
)
gid_seas_source <- gid_seas_source[
  gid_seas_source$gid %in% unique(as.character(surv_dat_all$gid)),
  ,
  drop = FALSE
]

# LP -> gid member bridge sidecar (gid_season only): the panel keeps
# ONE representative location_period_id per gid, so a coverage mask bridged
# through the panel misses the other member LPs' reporting windows. Every
# build_coverage_mask() caller reads this file when it exists
# (utils/coverage.R::coverage_lp_bridge()) and falls back to the panel otherwise.
LP_MEMBERS_PATH <- NA_character_
if (gid_season) {
  LP_MEMBERS_PATH <- occ_panel_lp_members_path(SURV_DAT_PATH)
  lp_members <- ep_res$member_rows %>%
    ungroup() %>%
    transmute(
      location_period_id = as.character(location_period_id),
      gid = as.character(gid),
      country = as.character(country)
    ) %>%
    filter(
      !is.na(location_period_id),
      !is.na(gid),
      gid %in% unique(as.character(surv_dat_all$gid))
    ) %>%
    distinct() %>%
    arrange(gid, location_period_id)
}

# Contagion coverage sidecar (--contagion_coverage_report only).
if (isTRUE(opt$contagion_coverage_report)) {
  COV_SIDECAR_PATH <- contagion_coverage_csv_path(SURV_DAT_PATH)
  contagion_report$contagion_coverage_path <- COV_SIDECAR_PATH
}

# Panel-mask sidecars (gid_season; read by 09's report_occ_panel_exclusions,
# table D31): the ADM1-overlap parent episodes, the discarded onsets, the ADM0
# rows dropped, and -- only when a mask is on -- the masked rows.
OCC_SIDECARS <- if (gid_season) {
  occ_panel_sidecar_paths(SURV_DAT_PATH)
} else {
  list()
}

# Policy sidecar. The occurrence artefact's NAME is a function of the baseline
# only, so two runs under different --composite_policy values overwrite one
# another with no trace in the filename. make_surv_dat_file() cannot be extended to
# carry the policy without touching its ~15 call sites (most pass a bare
# `list(baseline = ...)`), so the record goes in a companion file instead.
# See prep_policy_stamp_path() in utils.R.
PANEL_STAMP <- prep_policy_stamp(
  pop_source = opt$pop_source,
  composite_policy = opt$composite_policy,
  cover_threshold = opt$cover_threshold,
  composite_tiling = opt$composite_tiling,
  # The occurrence track carries NO population column into the model frame --
  # its offset is time-based (log interval length), so there is no population
  # vector to fingerprint. Recorded as NA rather than as a digest of an empty
  # vector, which would look like a real value and compare equal across
  # genuinely different runs.
  pop_digest = NA_character_,
  episode_key = opt$episode_key,
  episode_digest = ep_digest,
  counts = c(
    list(
      track = "occurrence",
      baseline = opt$baseline,
      # Which key the section-8 o_ct/o_tu leave-one-out lags were grouped on --
      # "t_end" (legacy, leaky) or "calendar" (leakage-free). See
      # --spillover_grouping above; recorded here (not in the filename) for the
      # same reason as composite_policy.
      spillover_grouping = opt$spillover_grouping,
      # Distance-weighted neighbour-contagion kernel inputs. Like spillover_grouping,
      # these do not appear in the artefact name, so they are recorded here.
      contagion_rhos = paste(CONTAGION_RHOS, collapse = ","),
      contagion_gadm_path = opt$contagion_gadm_path,
      contagion_cols = paste(NBR_COLS, collapse = ","),
      # The _pg registry (the first listed; what assert_contagion_partition()
      # compares with grouping.source), then every listed registry with its
      # content digest (the `_part<token>` columns; compared at fit time).
      contagion_partition_source = if (length(CONTAGION_PARTITION_SOURCES)) {
        CONTAGION_PARTITION_SOURCES[[1L]]
      } else {
        ""
      },
      contagion_partition_sources = paste(
        CONTAGION_PARTITION_SOURCES,
        collapse = ","
      ),
      contagion_partition_digests = paste(
        vapply(CONTAGION_PARTITION_SOURCES, registry_digest, character(1)),
        collapse = ","
      ),
      enso_max_lead = opt$enso_max_lead,
      n_lp_total = if (use_xwalk) {
        n_lp_before
      } else {
        dplyr::n_distinct(surv_dat_all$location_period_id)
      },
      n_lp_dropped = if (use_xwalk) length(dropped_lp) else 0L,
      n_rows_in = if (use_xwalk) n_before else NA_integer_,
      n_rows_out = nrow(surv_dat_all),
      n_events_out = sum(as.logical(surv_dat_all$outbreak), na.rm = TRUE),
      dropped_lp = if (use_xwalk) dropped_lp else character(0),
      # -- Point 1/2 seasonality-match traceability --------------------------
      n_seasgid_spatial_checked = n_seasgid_spatial_checked,
      n_seasgid_spatial_recovered = n_seasgid_spatial_recovered,
      n_seas_episodes_matched = n_seas_episodes_matched,
      n_seas_episodes_defaulted_jan = n_seas_episodes_defaulted_jan,
      # per-gid calendar source (sidecar above)
      min_seas_source_path = MIN_SEAS_SOURCE_PATH,
      n_gid_min_seas_estimate = sum(
        gid_seas_source$min_seas_source != "default_january"
      ),
      n_gid_min_seas_default_jan = sum(
        gid_seas_source$min_seas_source == "default_january"
      ),
      # LP -> gid member bridge sidecar (gid_season only; NA otherwise)
      lp_members_path = LP_MEMBERS_PATH,
      # panel-mask sidecars (gid_season only; NA otherwise)
      occ_panel_sidecars = if (length(OCC_SIDECARS)) {
        paste(basename(unlist(OCC_SIDECARS)), collapse = ",")
      } else {
        NA_character_
      }
    ),
    # gid_season only (empty otherwise -- lp_season keeps its field set).
    ep_report,
    # --national_zero_fill only (empty otherwise, so default stamps are
    # unchanged). Known limitation: the climate-lag season grouping
    # (occ_data.R, t_start_seas resets) can chain a partial season that follows
    # a gap onto the previous one; only the <clim>_lag* columns are affected.
    nz_report,
    # panel masks (gid_season only): n_* / cases_* counts (table S1) and the
    # spec fields when non-default.
    panel_report,
    # --contagion_coverage_report only: contagion_* fields,
    # kept OUT of table S1 (report_cascade reads n_* / cases_* / median_*).
    contagion_report
  )
)

# Overwrite guard: a panel of the same name built differently (a
# columns-only change such as --contagion_coverage_report, or another vintage)
# must not silently replace the one a deployed run was fitted on.
assert_prep_stamp_overwrite(
  SURV_DAT_PATH,
  PANEL_STAMP,
  force = isTRUE(opt$force),
  rerun_hint = "Rscript analysis/03_prepare_occ.R ... --force"
)

saveRDS(surv_dat_all, file = SURV_DAT_PATH)
utils::write.csv(gid_seas_source, MIN_SEAS_SOURCE_PATH, row.names = FALSE)
if (gid_season) {
  utils::write.csv(lp_members, LP_MEMBERS_PATH, row.names = FALSE)
  message(sprintf(
    "[occ/coverage] LP member bridge: %d LP-gid pair(s) over %d gid(s) -> %s",
    nrow(lp_members),
    dplyr::n_distinct(lp_members$gid),
    LP_MEMBERS_PATH
  ))
}
if (isTRUE(opt$contagion_coverage_report)) {
  utils::write.csv(cov_side, COV_SIDECAR_PATH, row.names = FALSE)
  message("[occ/contagion-coverage] sidecar -> ", COV_SIDECAR_PATH)
}
if (gid_season) {
  utils::write.csv(
    adm1_overlap_tab,
    OCC_SIDECARS$adm1_overlap,
    row.names = FALSE
  )
  utils::write.csv(
    discarded_tab,
    OCC_SIDECARS$discarded_onsets,
    row.names = FALSE
  )
  utils::write.csv(adm0_tab, OCC_SIDECARS$adm0_dropped, row.names = FALSE)
  if (!is.null(masked_rows)) {
    utils::write.csv(
      masked_rows %>%
        mutate(
          country = gid_meta$country[match(gid, gid_meta$gid)],
          .after = gid
        ),
      OCC_SIDECARS$masked_rows,
      row.names = FALSE
    )
  } else if (file.exists(OCC_SIDECARS$masked_rows)) {
    unlink(OCC_SIDECARS$masked_rows)
  }
  message(
    "[occ/panel] item 11 sidecars: ",
    paste(basename(unlist(OCC_SIDECARS)), collapse = ", ")
  )
}
message(sprintf(
  paste0(
    "[occ/seas] min_seas source per gid: %d estimate (gid), %d estimate (ADM1), ",
    "%d January default; %d panel gid(s) without member rows -> %s"
  ),
  sum(gid_seas_source$min_seas_source == "estimate_gid"),
  sum(gid_seas_source$min_seas_source == "estimate_adm1"),
  sum(gid_seas_source$min_seas_source == "default_january"),
  length(setdiff(unique(as.character(surv_dat_all$gid)), gid_seas_source$gid)),
  MIN_SEAS_SOURCE_PATH
))
persist_prep_policy_stamp(SURV_DAT_PATH, PANEL_STAMP)

message(sprintf(
  "[occ/composite] policy stamp written: %s",
  prep_policy_stamp_path(SURV_DAT_PATH)
))

# 12. Plots --------------------------------------------------------------------

if (opt$make_figures) {
  cat("-- Making survival data figures for all location periods \n")

  future_walk(
    unique(surv_dat_all$location_period_id),
    ~ save_plot_surv_dat(
      .,
      surv_dat_all = surv_dat_all,
      all_ts_dat = all_ts_dat
    )
  )
}
