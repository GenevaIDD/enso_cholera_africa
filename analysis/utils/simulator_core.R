# simulator_core.R — kernel of the AR forward outbreak simulator
#
# Provides:
#   build_location_frame()      one row per location (min_seas, tu, pop,
#                               admin_level rehydrated when absent)
#   build_climate_panel_for_tu(), build_warmup_state()   per-TU inputs
#   extract_occurrence_draws(), extract_duration_draws(),
#   extract_monthly_draws()     posterior draws bundles of the brms fits
#   simulate_one_tu()           one (TU, replicate-chunk) task
#   derive_min_seas_from_panel(), load_hmm_tu_mapping(),
#   load_population_per_location()   frame helpers
#   make_draws_bundle_synthetic()    deterministic test fixture
#
# Details:
#   Sourced by 07a / 07b_simulator_v1.R, 07b_simulator_spatpomp.R and
#   smoothing_core.R (the spatPomp kernels reuse its primitives),
#   03_prepare_sim_panel*.R, 07c, 07d and 09b (location frame), and
#   stand-alone by test fixtures. Notation: location i, month m, draw s.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(qs2)
})

# ---------------------------------------------------------------------------
# Location frame
# ---------------------------------------------------------------------------

#' .pop_artefact_or_null
#'
#' Resolve a population-artefact path to itself when it exists, or NULL when it
#' does not. Used at every build_location_frame() call site so a DEFAULTED path
#' (make_pop_run_file(...)) degrades to the legacy route on a tree where
#' 02c_preprocess_population_compile.R has not run, while an EXPLICITLY configured path
#' that is missing still hard-errors inside build_location_frame(). The
#' asymmetry is the point: a default is a guess, a config key is an instruction.
.pop_artefact_or_null <- function(path) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) path else NULL
}

#' .coarsen_admin_level_to_available
#'
#' Walk a per-location admin_level vector up the ADM3 -> ADM2 -> ADM1 ladder
#' until it lands on a level actually present in `available` (typically
#' `unique(brt$cell$admin_level)` from a base-rate anchor table). A location whose
#' own level is already in `available` is returned unchanged (no-op on any
#' table that covers every level, e.g. a `cross_absent_levels = TRUE` table).
#' A location whose level is NA, unrecognised, or has no coarser match in
#' `available` either is also returned unchanged -- this function only ever
#' narrows the gap `resolve_mean_rate()` would otherwise T5 on; it never
#' invents a level that doesn't correspond to a real coarser administrative
#' parent, and it never masks a genuine "no fallback exists" T5.
#'
#' Exists because P0-P4 estimand-coherent anchor tables
#' (`cross_absent_levels = FALSE`) restrict cell/group/global crossing to
#' exactly the admin levels present in the track's training corpus (e.g.
#' ADM1/ADM2 only, per `admin_level: ADM1-ADM2`) -- so a simulator location
#' frame whose `admin_level` comes from the (broader-scope) occurrence panel
#' can legitimately contain real ADM3 locations that track was never fit on.
#' See resolve_mean_rate() (T5).
#'
#' @param admin_level character vector.
#' @param available character vector of levels present in the target table.
#' @return character vector, same length as `admin_level`.
.coarsen_admin_level_to_available <- function(admin_level, available) {
  ladder <- c("ADM3", "ADM2", "ADM1")
  vapply(
    admin_level,
    function(al) {
      if (is.na(al) || al %in% available) {
        return(al)
      }
      pos <- match(al, ladder)
      if (is.na(pos) || pos >= length(ladder)) {
        return(al)
      }
      for (j in seq(pos + 1L, length(ladder))) {
        if (ladder[j] %in% available) return(ladder[j])
      }
      al
    },
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )
}

#' Resolve one empirical baseline-rate table onto the location frame.
#'
#' The three anchors (`base_rate_occ` / `base_rate_dur` / `base_rate_mon`,
#' `utils/anchor.R`) are the same object: a country x admin_level table,
#' resolved per-unit through the T1->T5 ladder, entered into the linear
#' predictor as a fixed `offset(log(.))`. They are therefore loaded,
#' provenance-checked, level-coarsened and resolved by this one function rather
#' than by near-copies of the same 60 lines -- the kind of duplication that
#' drifts silently, and where drift means the simulator resolves a different
#' rate than the fit did.
#'
#' @param frame the location frame being built (needs `country`, `admin_level`,
#'   and, when `check_pop_provenance`, `pop_tier`).
#' @param path sidecar written by `persist_base_rate_table()` /
#'   `persist_mean_rate_table()`. Guarded by the caller.
#' @param label human name for messages, e.g. "base_rate_mon".
#' @param column,tier_column columns to write onto `frame`.
#' @param check_pop_provenance TRUE only for tables whose DENOMINATOR contains
#'   population (monthly). `mu = pop x rate x exp(LP)` multiplies the two,
#'   so pairing a post-fix population with a pre-fix rate compounds rather than
#'   cancels. occ/dur divide by exposure days alone, so no population can make
#'   their table stale and the check would only raise false alarms.
#' @param gaps the running rehydration log.
#' @return list(frame, gaps).
.resolve_frame_anchor <- function(
  frame,
  path,
  label,
  column,
  tier_column,
  check_pop_provenance,
  gaps
) {
  brt <- readRDS(path)

  # --- rate x population ordering guard (Step 7) ---------------------------
  # This is the one seam where the ordering error cannot be caught by
  # arithmetic. persist_base_rate_table() reconciles at FIT time because it holds
  # both the fresh in-memory table and the sidecar. Here there is only the
  # sidecar -- the simulator has no independent estimate of the rate to compare
  # against, so a table built from the OLD populations would be used without
  # complaint against the NEW ones resolved just above.
  if (isTRUE(check_pop_provenance)) {
    prov <- pop_provenance_of(brt)
    if (is.null(prov)) {
      unverified_msg <- sprintf(
        paste0(
          label,
          " provenance UNVERIFIED (%s carries no population ",
          "stamp -- built before the guard existed; rebuild it if the ",
          "population definition has changed since)"
        ),
        basename(path)
      )
      gaps <- c(gaps, unverified_msg)
      # Layer 3 (fail-closed coverage, see pop_stamp_sidecar_path() in
      # utils.R): fit time (persist_mean_rate_table()) now HALTS by default on an
      # unstamped pairing (assert_pop_provenance()'s require_stamp = TRUE), but
      # this simulator load path has no independent "current" population
      # stamp to compare against -- only the sidecar's own stamp -- so a hard
      # stop here would fail closed on artefacts this function cannot actually
      # judge. A `warning()` (not just a `message()`-logged `gaps` entry, easy
      # to miss in a long run log) is the honest middle ground: visible by
      # default, but not fatal on a check this function cannot complete.
      warning(unverified_msg, call. = FALSE)
    } else {
      gaps <- c(
        gaps,
        sprintf(
          "%s provenance: pop_source=%s, pop_total=%.0f, built %s",
          label,
          prov$pop_source,
          prov$pop_total %||% NA_real_,
          prov$stamped_at %||% "?"
        )
      )
      # `unstamped` is the stamp build_base_rate_table() writes when the input data
      # carried no population provenance of its own: the digest is real but the
      # SOURCE is unknown, so this function cannot tell whether the table was
      # built on the legacy `pop` or the canonical Q_g. That is the same
      # "check cannot be completed" situation as the `is.null(prov)` branch
      # above, and it gets the same treatment -- a visible warning rather than a
      # stop, because failing closed here would reject artefacts this function
      # has no evidence against. Flagged separately from the hard `legacy` case
      # below, which IS evidence.
      if (
        identical(prov$pop_source, "unstamped") &&
          "pop_tier" %in% names(frame) &&
          any(frame$pop_tier %in% c("direct", "wpp_rescale", "parent_share"))
      ) {
        unstamped_msg <- sprintf(
          paste0(
            "%s provenance UNVERIFIABLE: the frame's populations came from the ",
            "canonical per-gid table (Q_g), but %s was built from data carrying ",
            "no population stamp (pop_source = 'unstamped'), so the pairing ",
            "cannot be checked. Rebuild it if the population definition has ",
            "changed since."
          ),
          label,
          basename(path)
        )
        gaps <- c(gaps, unstamped_msg)
        warning(unstamped_msg, call. = FALSE)
      }
      # A frame whose populations came from the canonical Q_g table must not be
      # paired with rates derived from the pre-crosswalk `pop`. `pop_tier` is
      # set by the population block above and is the frame's own provenance.
      if (
        identical(prov$pop_source, "legacy") &&
          "pop_tier" %in% names(frame) &&
          any(frame$pop_tier %in% c("direct", "wpp_rescale", "parent_share"))
      ) {
        stop(
          "build_location_frame: the frame's populations were resolved from ",
          "the canonical per-gid table (Q_g), but the ",
          label,
          " sidecar at ",
          path,
          " was built from the LEGACY population ",
          "(pop_source = 'legacy'). mu = pop * rate * exp(LP), so ",
          "pairing a new population with an old rate compounds rather than ",
          "cancels.\nRequired order: pop fix -> ",
          label,
          " rebuild -> ",
          "refit -> Phase-A smoothing rebuild -> forecast."
        )
      }
    }
  }

  # --- which admin_level column keys this table (P3) ------------------------
  # `build_base_rate_table()` (utils/anchor.R) canonicalises the admin_level of the
  # mon/siz training panels before it aggregates -- `canonical_admin_level`,
  # defaulting to TRUE for any non-`footprint` denominator -- while
  # `.occ_anchor_view()`/`.dur_anchor_view()` deliberately do not. So the table's own
  # estimand stamp, not the track name, decides which of the frame's two
  # admin_level columns is the right join key. A table with no stamp (occ/dur,
  # and any pre-P0 mon/siz sidecar) keeps the corpus key, as does a frame built
  # without `gid_pop_path` -- both are exactly the previous behaviour.
  est_stamp <- brt$meta$mon_estimand %||% brt$meta$size_estimand
  key_column <- "admin_level"
  if (!is.null(est_stamp) && "admin_level_canonical" %in% names(frame)) {
    canon <- est_stamp$canonical_admin_level %||%
      !identical(
        as.character(est_stamp$denominator %||% "footprint"),
        "footprint"
      )
    if (isTRUE(canon)) key_column <- "admin_level_canonical"
  }
  # occ/dur tables built with occurrence/duration.canonical_admin_level carry a
  # plain meta flag instead of an estimand stamp (utils/anchor.R).
  if (
    is.null(est_stamp) &&
      isTRUE(brt$meta$canonical_admin_level) &&
      "admin_level_canonical" %in% names(frame)
  ) {
    key_column <- "admin_level_canonical"
  }
  frame_level <- as.character(frame[[key_column]])
  if (identical(key_column, "admin_level_canonical")) {
    n_rekeyed <- sum(
      !is.na(frame_level) &
        !is.na(frame$admin_level) &
        frame_level != as.character(frame$admin_level)
    )
    gaps <- c(
      gaps,
      sprintf(
        "%s keyed on admin_level_canonical (estimand denominator '%s'); %d unit(s) differ from the corpus label",
        label,
        as.character(est_stamp$denominator %||% "?"),
        n_rekeyed
      )
    )
  }

  # Cap the LOOKUP key at the finest level the model was actually fit on -- not
  # frame$admin_level itself (used nowhere else downstream, see the comment
  # above this block). P0-P4 estimand-coherent tables
  # (cross_absent_levels=FALSE) restrict cell/group/global crossing to exactly
  # the admin levels present in the training corpus (e.g. ADM1/ADM2 only, per
  # admin_level: ADM1-ADM2), so the resolver hard T5-stops on any location
  # whose real (GADM/canonical) admin level is finer than that -- e.g. ADM3
  # locations that are in scope for the OCCURRENCE panel (admin_level_data) but
  # were never in scope for the size fit. Rather than dropping those locations
  # from the simulator entirely, borrow the coarser parent level's rate for the
  # LOOKUP only -- the same "resolve at the coarsest available match" principle
  # the resolver already applies across country/group/global for a missing
  # CELL, just applied across LEVEL for a level absent from the table
  # altogether. A no-op whenever the table already has the location's own level
  # (e.g. tables built with cross_absent_levels=TRUE).
  admin_level_for_rate <- .coarsen_admin_level_to_available(
    frame_level,
    unique(brt$cell$admin_level)
  )
  n_capped <- sum(
    !is.na(frame_level) &
      admin_level_for_rate != frame_level
  )
  if (n_capped > 0) {
    gaps <- c(
      gaps,
      sprintf(
        "%s lookup: %d location(s) coarsened to a level present in %s (own admin_level absent from the table)",
        label,
        n_capped,
        basename(path)
      )
    )
  }

  resolved <- resolve_mean_rate(frame$country, admin_level_for_rate, brt)
  frame[[column]] <- resolved$mean_rate
  frame[[tier_column]] <- resolved$mr_tier
  tier_tab <- table(resolved$mr_tier)
  gaps <- c(
    gaps,
    sprintf(
      "%s (resolved per-unit from %s; tiers: %s)",
      label,
      basename(path),
      paste(sprintf("%s=%d", names(tier_tab), tier_tab), collapse = ", ")
    )
  )

  list(frame = frame, gaps = gaps)
}

#' Derive a per-location frame from the prepared survival data.
#'
#' One row per (gid, location_period_id) with what the kernel needs: country,
#' AFRO_region, TU, `min_seas`, population, admin_level and any opt-in
#' per-track baseline-rate anchors. Columns missing from the survival panel
#' are backfilled from their upstream sources, logged in one summary line.
#'
#' @param hmm_clusters_path path to `data/hmm_country_clusters.csv` (the file
#'   `03_prepare_occ.R` reads). NULL: `tu` is NA and the TU spillover lag is
#'   inactive.
#' @param outbreak_metadata_path unused; kept for a future WorldPop × UN-WPP
#'   population scaling path.
#' @param size_data_path path to `generated_data/outbreak_size_analysis_data.rds`
#'   (as read by `read_size_data()`). Gives per-gid `pop` (mean over training
#'   rows) when `gid_pop_path` is NULL and fills gaps in `admin_level`. Still a
#'   simulator input although the size track is retired.
#' @param base_rate_table_paths optional NAMED list of anchor sidecars,
#'   `list(occ = ..., dur = ..., mon = ...)`, as written by
#'   `persist_base_rate_table(opt, track, data)`. Each goes through the
#'   country x admin-level ladder into its own `base_rate_<track>` + `_tier`
#'   columns (the tracks share one frame, their denominators differ). NULL adds
#'   nothing. An absent or empty name is skipped; a named but missing file is a
#'   hard error.
#' @param admin_level_data optional data frame with `gid` and `admin_level`
#'   (e.g. the occurrence survival panel a caller already loaded for
#'   `drop_overlapping_adm1_from_frame()`). Backfills `admin_level` when
#'   `survival_data` is a bare climate panel; falls through to
#'   `size_data_path` (mode per gid) when NULL or a gid is missing.
#' @param gid_pop_path path to the canonical Q_g table
#'   `gid_population_<Y0>-<Y1>.qs` (`make_pop_run_file(opt, "pop_gid")`). When
#'   set, `pop` is resolved per gid by `resolve_gid_population()` (no silent
#'   default). NULL keeps the legacy size-training join.
#' @param wpp_path path to `wpp2024_national_totals.csv`; enables the T2
#'   nearest-year rescale for forecast years outside WorldPop's 2015-2030
#'   coverage. NULL disables T2 (T1/T3 still apply).
#' @param pop_year integer year at which to resolve the gid population, or a
#'   mode (sim_pop_year_spec(): `season_year` / `calendar_year`). A mode
#'   resolves the static `pop` at year(pop_window[1]) and adds `pop_y<YYYY>`
#'   columns for year(pop_window[1]) - 1 to year(pop_window[2]) plus
#'   `pop_shift_m` (min_seas - 1 under season_year, 0 under calendar_year),
#'   read per month by frame_pop_at() / frame_pop_matrix(). A mode needs
#'   `gid_pop_path` and a panel without `pop`. An integer adds no column.
#' @param pop_window length-2 Date (first, last month the run reads;
#'   sim_pop_window()); required in a mode, ignored otherwise.
build_location_frame <- function(
  survival_data,
  hmm_clusters_path = NULL,
  outbreak_metadata_path = NULL,
  size_data_path = NULL,
  base_rate_table_paths = NULL,
  admin_level_data = NULL,
  gid_pop_path = NULL,
  wpp_path = NULL,
  pop_year = 2020L,
  pop_window = NULL
) {
  pop_mode <- is.character(pop_year)
  if (pop_mode) {
    if (!(length(pop_year) == 1L && pop_year %in% .SIM_POP_YEAR_MODES)) {
      stop(
        "build_location_frame: pop_year must be an integer year or one of ",
        paste(.SIM_POP_YEAR_MODES, collapse = " / "),
        "; got ",
        pop_year
      )
    }
    pop_window <- as.Date(pop_window)
    if (
      length(pop_window) != 2L ||
        anyNA(pop_window) ||
        pop_window[1] > pop_window[2]
    ) {
      stop(
        "build_location_frame: pop_year = '",
        pop_year,
        "' needs pop_window = ",
        "c(first month, last month) of the run (sim_pop_window())."
      )
    }
    if (is.null(gid_pop_path) || !nzchar(gid_pop_path)) {
      stop(
        "build_location_frame: pop_year = '",
        pop_year,
        "' needs the canonical gid ",
        "population table (gid_pop_path); the size-training route has no year."
      )
    }
    if ("pop" %in% names(survival_data)) {
      stop(
        "build_location_frame: pop_year = '",
        pop_year,
        "' cannot apply to a panel ",
        "that already carries a `pop` column (a fixed population)."
      )
    }
  }
  pop_year_static <- if (pop_mode) {
    as.integer(format(pop_window[1], "%Y"))
  } else {
    pop_year
  }
  required <- c("gid", "country", "AFRO_region")
  missing <- setdiff(required, names(survival_data))
  if (length(missing) > 0) {
    stop(
      "build_location_frame: survival data missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  has_min_seas <- "min_seas" %in% names(survival_data)
  has_tu <- "tu" %in% names(survival_data)
  has_pop <- "pop" %in% names(survival_data)
  has_admin_level <- "admin_level" %in% names(survival_data)

  # "Already a column ⇒ pass it through, don't re-resolve", for the three
  # anchors. A caller handing us a panel that a
  # loader already anchored (attach_base_rate()) keeps that exact column rather
  # than silently re-deriving it from a sidecar that may have been built under a
  # different estimator.
  anchor_pass <- intersect(
    c("base_rate_occ", "base_rate_dur", "base_rate_mon"),
    names(survival_data)
  )

  base_cols <- c("gid", "country", "AFRO_region")
  if (has_min_seas) {
    base_cols <- c(base_cols, "min_seas")
  }
  if (has_tu) {
    base_cols <- c(base_cols, "tu")
  }
  if (has_pop) {
    base_cols <- c(base_cols, "pop")
  }
  if (length(anchor_pass)) {
    base_cols <- c(base_cols, anchor_pass)
  }
  if (has_admin_level) {
    base_cols <- c(base_cols, "admin_level")
  }

  frame <- survival_data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(base_cols))) %>%
    dplyr::mutate(
      # Flatten dots to underscores before extracting, so the ADM1 key is
      # derived correctly under both gid conventions seen in this repo
      # (flat "ISO3_N" and GADM-dotted "ISO3.ADM1[.ADM2]_N") and is always
      # emitted in the canonical underscore form "ISO3_N" that the brms RE
      # level-name keys (occ_post.R) and the test fixtures expect.
      # Matching the flatten-first idiom of 03_prepare_occ.R,
      # but inlined rather than calling flatten_gid() because gid.R is not in
      # the source closure of every build_location_frame() caller
      # (e.g. 07b_simulator_spatpomp.R). Same rule as gid.R::gid_to_adm1(),
      # including the ISO3 separator GADM 4.1 omits for Ghana ("GHA5_2" ->
      # "GHA_5"); an unparseable gid stops below instead of becoming NA.
      adm1_unit = stringr::str_extract(
        sub(
          "^([A-Z]{3})([0-9])",
          "\\1_\\2",
          stringr::str_replace_all(gid, "\\.", "_")
        ),
        "^[A-Z]{3}_[0-9]+"
      ),
      gid = as.character(gid),
      country = as.character(country)
    )
  bad_gid <- unique(frame$gid[!is.na(frame$gid) & is.na(frame$adm1_unit)])
  if (length(bad_gid)) {
    stop(sprintf(
      "[build_location_frame] %d gid(s) with no parseable ADM1: %s",
      length(bad_gid),
      paste(utils::head(bad_gid, 10), collapse = ", ")
    ))
  }
  if (has_admin_level) {
    frame$admin_level <- as.character(frame$admin_level)
  }

  # TZA-mainland exclusion (`uid == "AFR::TZA::Mainland-1"`) is already
  # baked into the fit's data — the simulator inherits whatever survived
  # the fit-time filtering. Defensive guard: drop rows whose gid does not
  # match the ADM1 regex (i.e. would have produced NA spillover lags), e.g.
  # the synthetic `admin_level == "ADM0"` mainland row, which has no cell in
  # the anchor tables (valid levels are ADM1-3 only) and would otherwise
  # hard-error the resolution below. Runs ahead of the pop/anchor
  # rehydration blocks so that a row this filter would drop anyway is never
  # asked to resolve against an anchor table.
  bad_gid <- is.na(frame$adm1_unit)
  if (any(bad_gid)) {
    message(sprintf(
      "[build_location_frame] dropping %d location(s) with non-standard gid (e.g. TZA-mainland)",
      sum(bad_gid)
    ))
    frame <- frame[!bad_gid, , drop = FALSE]
  }

  gaps <- character(0)

  # --- min_seas ---------------------------------------------------------
  # Gids whose min_seas is a default (1), not their own season start: returned
  # as attr "min_seas_defaulted" (set only when non-empty), so the season onset
  # gate can refuse them (check_onset_eligibility_frame()).
  ms_defaulted <- character(0)
  if (!has_min_seas) {
    if (
      "t_start_seas" %in%
        names(survival_data) &&
        "t_start" %in% names(survival_data)
    ) {
      ms <- derive_min_seas_from_panel(survival_data)
      frame <- dplyr::left_join(frame, ms, by = "gid")
      n_nas <- sum(is.na(frame$min_seas))
      ms_defaulted <- as.character(frame$gid[is.na(frame$min_seas)])
      if (n_nas > 0) {
        frame$min_seas[is.na(frame$min_seas)] <- 1L
      }
      gaps <- c(
        gaps,
        sprintf(
          "min_seas (derived from t_start_seas resets%s)",
          if (n_nas) {
            sprintf("; %d gid(s) had no reset row → defaulted to 1", n_nas)
          } else {
            ""
          }
        )
      )
    } else {
      frame$min_seas <- 1L
      ms_defaulted <- as.character(frame$gid)
      gaps <- c(
        gaps,
        "min_seas (defaulted to 1; panel lacks t_start/t_start_seas)"
      )
    }
  }

  # --- tu ---------------------------------------------------------------
  if (!has_tu) {
    if (
      !is.null(hmm_clusters_path) &&
        nzchar(hmm_clusters_path) &&
        file.exists(hmm_clusters_path)
    ) {
      tu_map <- load_hmm_tu_mapping(hmm_clusters_path)
      frame <- frame %>%
        dplyr::mutate(.iso3 = stringr::str_extract(adm1_unit, "^[A-Z]{3}")) %>%
        dplyr::left_join(tu_map, by = c(.iso3 = "iso3")) %>%
        dplyr::select(-.iso3)
      n_nas <- sum(is.na(frame$tu))
      gaps <- c(
        gaps,
        sprintf(
          "tu (joined from %s%s)",
          basename(hmm_clusters_path),
          if (n_nas) {
            sprintf(
              "; %d location(s) NA — country absent from cluster file (TU lag inactive there, matching fit-time)",
              n_nas
            )
          } else {
            ""
          }
        )
      )
    } else {
      frame$tu <- NA_character_
      gaps <- c(
        gaps,
        "tu (NA; no HMM clusters path — TU spillover lag inactive everywhere)"
      )
    }
  }

  # --- pop --------------------------------------------------------------
  # No silent default population (a fallback pop = 1 is one person). The
  # simulator reads the canonical whole-gid population Q_g, NOT the size
  # training `pop`, which is the reporting-footprint population H (H = c * Q,
  # c < 1 under partial coverage). load_population_per_location() below is the
  # legacy route only (no gid_pop_path). gid_pop is hoisted so the
  # canonical-admin_level block reuses it.
  gid_pop <- NULL
  if (!has_pop) {
    if (!is.null(gid_pop_path) && nzchar(gid_pop_path)) {
      if (!file.exists(gid_pop_path)) {
        stop(
          "build_location_frame: gid_pop_path not found (",
          gid_pop_path,
          "). Build it with analysis/02c_preprocess_population_compile.R, or pass NULL ",
          "to fall back to the size-training join."
        )
      }
      gid_pop <- load_gid_population(gid_pop_path)
      wpp <- if (!is.null(wpp_path) && nzchar(wpp_path)) {
        load_wpp_table(wpp_path)
      } else {
        NULL
      }
      resolved <- resolve_gid_population(
        frame$gid,
        pop_year_static,
        gid_pop,
        wpp,
        on_missing = "error"
      )
      frame$pop <- resolved$pop
      frame$pop_tier <- resolved$pop_tier
      tier_tab <- table(resolved$pop_tier)
      gaps <- c(
        gaps,
        sprintf(
          "pop (resolved per-gid from %s at year %d; tiers: %s)",
          basename(gid_pop_path),
          as.integer(pop_year_static),
          paste(sprintf("%s=%d", names(tier_tab), tier_tab), collapse = ", ")
        )
      )
      if (pop_mode) {
        # Per-year columns for the per-month lookup (frame_pop_at()). One
        # vectorised resolution over every (gid, year); the year before the
        # window start covers a season year that began in the previous
        # calendar year.
        yrs <- seq(
          as.integer(format(pop_window[1], "%Y")) - 1L,
          as.integer(format(pop_window[2], "%Y"))
        )
        ng <- nrow(frame)
        res_y <- resolve_gid_population(
          rep(frame$gid, times = length(yrs)),
          rep(yrs, each = ng),
          gid_pop,
          wpp,
          on_missing = "error"
        )
        for (k in seq_along(yrs)) {
          frame[[sprintf("pop_y%04d", yrs[k])]] <- res_y$pop[
            (k - 1L) * ng + seq_len(ng)
          ]
        }
        frame$pop_shift_m <- if (identical(pop_year, "season_year")) {
          as.integer(frame$min_seas) - 1L
        } else {
          0L
        }
        tier_y <- table(res_y$pop_tier)
        gaps <- c(
          gaps,
          sprintf(
            "pop_y%d..pop_y%d (pop_year = %s: per-month population; tiers: %s)",
            min(yrs),
            max(yrs),
            pop_year,
            paste(sprintf("%s=%d", names(tier_y), tier_y), collapse = ", ")
          )
        )
      }
    } else if (
      !is.null(size_data_path) &&
        nzchar(size_data_path) &&
        file.exists(size_data_path)
    ) {
      # The legacy ROUTE survives; the legacy silent DEFAULT does not. A gid the
      # size training never saw now stops the run instead of forecasting from
      # one person.
      pop_map <- load_population_per_location(size_data_path)
      frame <- dplyr::left_join(frame, pop_map, by = "gid")
      n_nas <- sum(is.na(frame$pop))
      if (n_nas > 0) {
        stop(
          "build_location_frame: ",
          n_nas,
          " of ",
          nrow(frame),
          " frame unit(s) have no population in ",
          basename(size_data_path),
          " (e.g. ",
          paste(utils::head(frame$gid[is.na(frame$pop)], 5), collapse = ", "),
          "). The silent `pop <- 1.0` default is no longer permitted (D4). ",
          "Pass gid_pop_path = make_pop_run_file(opt, \"pop_gid\"), which ",
          "resolves every gid through the canonical Q_g table."
        )
      }
      gaps <- c(gaps, sprintf("pop (joined from %s)", basename(size_data_path)))
    } else {
      stop(
        "build_location_frame: no population source. Pass gid_pop_path ",
        "(make_pop_run_file(opt, \"pop_gid\")) or, for the legacy path, ",
        "size_data_path. The silent `pop <- 1.0` default is no longer ",
        "permitted (D4) -- it put 216 of 683 frame units into the forecast ",
        "with a population of one person."
      )
    }
  }

  # --- admin_level (needed to resolve the anchors below and by the mon
  #     `admin_levelADM2` term; not itself a rehydration target, so
  #     no `gaps` entry unless actually joined) -----------------------------
  # Priority: (1) already on survival_data: nothing to do. (2) admin_level_data
  # (typically the occurrence survival panel -- broader gid coverage, since it
  # is not restricted to the size model's min_pop/max_attack_rate/admin_level
  # training filters). (3) size_data_path (mode per gid over training rows --
  # weaker: only covers gids with an observed outbreak). A frame with no
  # admin_level source at all is legal; the anchor block below hard-errors only
  # when an anchor actually needs it.
  if (!has_admin_level) {
    src_label <- NULL
    if (
      !is.null(admin_level_data) &&
        all(c("gid", "admin_level") %in% names(admin_level_data))
    ) {
      al_map <- admin_level_data %>%
        dplyr::mutate(
          gid = as.character(gid),
          admin_level = as.character(admin_level)
        ) %>%
        dplyr::filter(!is.na(admin_level)) %>%
        dplyr::distinct(gid, admin_level) %>%
        dplyr::group_by(gid) %>%
        dplyr::slice(1) %>% # one row per gid; ties broken by first-seen
        dplyr::ungroup()
      frame <- dplyr::left_join(frame, al_map, by = "gid")
      src_label <- "admin_level_data"
    }
    still_na <- !("admin_level" %in% names(frame)) ||
      any(is.na(frame$admin_level))
    if (
      still_na &&
        !is.null(size_data_path) &&
        nzchar(size_data_path) &&
        file.exists(size_data_path)
    ) {
      al_map2 <- load_admin_level_per_location(size_data_path)
      if ("admin_level" %in% names(frame)) {
        na_gid <- frame$gid[is.na(frame$admin_level)]
        fill <- al_map2 %>%
          dplyr::filter(gid %in% na_gid) %>%
          dplyr::rename(admin_level_fill = admin_level)
        frame <- dplyr::left_join(frame, fill, by = "gid") %>%
          dplyr::mutate(
            admin_level = dplyr::coalesce(admin_level, admin_level_fill)
          ) %>%
          dplyr::select(-admin_level_fill)
      } else {
        frame <- dplyr::left_join(frame, al_map2, by = "gid")
      }
      src_label <- if (is.null(src_label)) {
        basename(size_data_path)
      } else {
        paste0(src_label, " + ", basename(size_data_path), " (gap-fill)")
      }
    }
    if (!is.null(src_label)) {
      n_nas <- if ("admin_level" %in% names(frame)) {
        sum(is.na(frame$admin_level))
      } else {
        nrow(frame)
      }
      gaps <- c(
        gaps,
        sprintf(
          "admin_level (joined from %s%s)",
          src_label,
          if (n_nas) sprintf("; %d location(s) still unresolved", n_nas) else ""
        )
      )
    }
  }

  # --- admin_level_canonical (P3, replay side) -----------------------------
  # The fit side re-keys admin_level to the gid's GADM depth
  # (apply_canonical_admin_level(), utils/population.R); this does the same on
  # the replay side for the level-keyed lookups (anchor cell, log_pop_c centre,
  # admin_levelADM2 dummy). A NEW column, not an overwrite: the occ / dur
  # anchors and drop_overlapping_adm1_from_frame() key off the corpus
  # `admin_level`, so each consumer picks its column. Inert without
  # gid_pop_path.
  if (
    !("admin_level_canonical" %in% names(frame)) &&
      !is.null(gid_pop_path) &&
      nzchar(gid_pop_path) &&
      file.exists(gid_pop_path)
  ) {
    if (is.null(gid_pop)) {
      gid_pop <- load_gid_population(gid_pop_path)
    }
    frame <- apply_canonical_admin_level(
      frame,
      gid_pop,
      overwrite = FALSE,
      tag = "build_location_frame/P3"
    )
    rep_al <- attr(frame, "admin_level_recode")
    # A gid absent from the canonical table keeps its corpus label rather than an
    # NA key -- same "unmatched gids are not silently re-labelled" rule
    # apply_canonical_admin_level() applies under `overwrite = TRUE`.
    if ("admin_level" %in% names(frame)) {
      frame$admin_level_canonical <- ifelse(
        is.na(frame$admin_level_canonical),
        as.character(frame$admin_level),
        frame$admin_level_canonical
      )
    }
    gaps <- c(
      gaps,
      sprintf(
        paste0(
          "admin_level_canonical (from %s; %d/%d unit(s) disagree with the ",
          "corpus label, %d unmatched and kept)"
        ),
        basename(gid_pop_path),
        rep_al$n_changed %||% 0L,
        nrow(frame),
        rep_al$n_unmatched %||% 0L
      )
    )
  }

  # --- the three anchors (occ / dur / mon) ----------------------------------
  # One per brms track (utils/anchor.R). Each is opt-in: a caller that passes no path for a track
  # gets no column, and the kernels' `offset_baserate` flags are FALSE for any
  # fit made without the anchor, so an un-anchored run is untouched.
  #
  # They are DISTINCT columns rather than one shared `base_rate`, because the
  # tracks share one location frame and their denominators differ
  # (exposure days for occ/dur, pop x exposure days for mon). One
  # column would silently hand the occurrence rate to the monthly kernel.
  if (!is.null(base_rate_table_paths)) {
    if (
      !is.list(base_rate_table_paths) || is.null(names(base_rate_table_paths))
    ) {
      stop(
        "build_location_frame: base_rate_table_paths must be a NAMED list, ",
        "e.g. list(occ = <path>, mon = <path>)."
      )
    }
    unknown <- setdiff(names(base_rate_table_paths), c("occ", "dur", "mon"))
    if (length(unknown)) {
      stop(
        "build_location_frame: unknown base_rate_table_paths entries: ",
        paste(unknown, collapse = ", "),
        " (expected occ / dur / mon)."
      )
    }

    if (!("admin_level" %in% names(frame)) || all(is.na(frame$admin_level))) {
      stop(
        "build_location_frame: base_rate_table_paths given but the frame ",
        "has no `admin_level` -- required to resolve a country x ",
        "admin-level anchor table. Pass `admin_level_data`."
      )
    }

    for (tr in names(base_rate_table_paths)) {
      p <- base_rate_table_paths[[tr]]
      if (is.null(p) || !nzchar(p)) {
        next
      }
      # Already carried through from `survival_data` (anchor_pass above).
      if (anchor_spec(tr)$column %in% anchor_pass) {
        next
      }
      if (!file.exists(p)) {
        stop(
          "build_location_frame: base_rate_table_paths$",
          tr,
          " not found (",
          p,
          "). No silent default (D4) -- pass the sidecar written by ",
          "persist_base_rate_table(opt, \"",
          tr,
          "\", data)."
        )
      }
      spec <- anchor_spec(tr)
      res <- .resolve_frame_anchor(
        frame,
        p,
        spec$column,
        spec$column,
        spec$tier_col,
        # Only mon's denominator contains population (pop x exposure days);
        # occ/dur divide by exposure days alone, so no population change can
        # stale their table.
        check_pop_provenance = isTRUE(spec$has_pop),
        gaps = gaps
      )
      frame <- res$frame
      gaps <- res$gaps
    }
  }

  if (length(gaps) > 0) {
    message("[build_location_frame] panel rehydration:")
    for (g in gaps) {
      message("  - ", g)
    }
  }

  frame <- tibble::as_tibble(frame)
  ms_defaulted <- intersect(ms_defaulted, as.character(frame$gid))
  if (length(ms_defaulted)) {
    attr(frame, "min_seas_defaulted") <- ms_defaulted
  }
  frame
}

# ---------------------------------------------------------------------------
# Rehydration helpers — read upstream metadata that the survival panel drops
# ---------------------------------------------------------------------------

#' Recover per-gid `min_seas` from the panel by inverting the seasonal offset.
#'
#' make_surv_data() (occ_data.R) computes `t_start_seas = yday(days_seas)-1`
#' with `days_seas = t_start %m+% months(-(min_seas-1))`. So `t_start_seas == 0`
#' occurs at the calendar month equal to `min_seas`. We pick, per gid, the row
#' with the smallest `abs(t_start_seas)` and shift its calendar month back by
#' `floor(t_start_seas / 30.44)` months to handle gids whose panel rows never
#' include the exact reset row.
derive_min_seas_from_panel <- function(panel) {
  if (!all(c("gid", "t_start", "t_start_seas") %in% names(panel))) {
    stop(
      "derive_min_seas_from_panel: need gid + t_start + t_start_seas columns."
    )
  }
  panel %>%
    dplyr::transmute(
      gid = as.character(gid),
      .month_t = lubridate::month(as.Date(t_start)),
      t_start_seas = as.numeric(t_start_seas)
    ) %>%
    dplyr::group_by(gid) %>%
    dplyr::slice(which.min(abs(t_start_seas))) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      gid = gid,
      min_seas = ((as.integer(.month_t) -
        1L -
        as.integer(floor(t_start_seas / 30.44)) +
        12L) %%
        12L) +
        1L
    )
}

#' Read the HMM transmission-unit (TU) cluster CSV.
#'
#' Mirrors `03_prepare_occ.R` — same file, same `iso3 → cluster` mapping,
#' renamed `cluster → tu` here to match the simulator's naming.
load_hmm_tu_mapping <- function(path) {
  if (!file.exists(path)) {
    stop("load_hmm_tu_mapping: HMM clusters file not found: ", path)
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("load_hmm_tu_mapping: readr required.")
  }
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::select(iso3, cluster) %>%
    dplyr::distinct() %>%
    dplyr::transmute(iso3 = as.character(iso3), tu = as.character(cluster))
}

#' Read per-gid mean population from the size-track training data.
#'
#' Reuses `generated_data/outbreak_size_analysis_data.rds` (built by
#' `03b_prepare_siz.R`) — the same file `read_size_data()` consumes at fit
#' time, so the simulator inherits the fit's `pop` baseline exactly.
load_population_per_location <- function(path) {
  if (!file.exists(path)) {
    stop("load_population_per_location: size data file not found: ", path)
  }
  size_dat <- readRDS(path)
  if (!all(c("gid", "pop") %in% names(size_dat))) {
    stop("load_population_per_location: size data missing gid/pop columns.")
  }
  size_dat %>%
    dplyr::filter(!is.na(pop), pop > 0) %>%
    dplyr::group_by(gid = as.character(gid)) %>%
    dplyr::summarise(pop = mean(pop, na.rm = TRUE), .groups = "drop")
}

#' Read per-gid admin_level from the size-track training data.
#'
#' Companion to `load_population_per_location()` — same source file, same
#' rehydration pattern — used by `build_location_frame()` to backfill
#' `admin_level` when the survival/climate panel lacks it (e.g. the
#' continuous climate panel from `03_prepare_sim_panel.R` carries no
#' `admin_level` column). `admin_level` is required to resolve the per-track
#' anchors from their country x admin-level tables. A gid should carry exactly one
#' admin_level; where duplicate/conflicting rows exist we take the mode
#' (ties broken by first-seen) rather than silently picking one at random.
load_admin_level_per_location <- function(path) {
  if (!file.exists(path)) {
    stop("load_admin_level_per_location: size data file not found: ", path)
  }
  size_dat <- readRDS(path)
  if (!all(c("gid", "admin_level") %in% names(size_dat))) {
    stop(
      "load_admin_level_per_location: size data missing gid/admin_level columns."
    )
  }
  size_dat %>%
    dplyr::mutate(
      gid = as.character(gid),
      admin_level = as.character(admin_level)
    ) %>%
    dplyr::filter(!is.na(admin_level)) %>%
    dplyr::count(gid, admin_level, name = "n") %>%
    dplyr::group_by(gid) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(gid, admin_level)
}

# ---------------------------------------------------------------------------
# Climate / covariate panel
# ---------------------------------------------------------------------------

#' Slice the climate panel down to the locations in one TU over the horizon.
#'
#' v1 only supports the
#' historical-replay `factual` scenario — the panel produced by
#' `03_prepare_occ.R` is consumed as-is. Counterfactual and CMIP6 scenario
#' branches (S9/S9b/S10) are stubs that error informatively and refer the
#' caller to the v2 scenario wiring.
#'
#' Expects `panel` to carry one row per (gid, month) with at minimum the climate
#' columns (`mean_precip`, `dev_precip`, `mean_temp_daymean`, `dev_temp_daymean`),
#' the three lag columns, and `t_start_seas` (within-season time used by the
#' basehaz). `t_start` is the first day of each month.
#'
#' @param scenario one of "factual" (v1), "counterfactual" / "ssp245" / "ssp585"
#'   (v2; errors).
#' @param warmup_months keep this many months of panel rows before
#'   `horizon$start` so the warm-up state can be read from real rows rather
#'   than recomputed. Default 12 (covers `enso_lag11`).
#' @param end_mode how `horizon$end` is read (utils/base.R window
#'   convention). The slice is `t_start <= end` and the hard
#'   coverage error (panel ends more than 31 days before `end`) is the same in
#'   every mode, so outputs are identical:
#'   "legacy" (default) = no further check (the former 31-day grace);
#'   "closed" = a forward horizon [start, end]: WARN when the end month itself
#'     is missing from the panel (a closed horizon would silently lose its last
#'     month; the generators make this a hard error at config time);
#'   "open"   = a half-open filter window [start, end): the month before `end`
#'     must be present (the end month need not be).
build_climate_panel_for_tu <- function(
  panel,
  frame_tu,
  horizon,
  scenario = "factual",
  warmup_months = 12L,
  end_mode = c("legacy", "closed", "open")
) {
  end_mode <- match.arg(end_mode)
  if (!identical(scenario, "factual")) {
    stop(
      sprintf(
        "build_climate_panel_for_tu: scenario='%s' is a v2 feature (S9/S9b/S10 in [[cholera-outbreak-simulator-suggestions]]). ",
        scenario
      ),
      "v1 only supports 'factual' historical replay; counterfactual + projection panels are not supported."
    )
  }

  required <- c(
    "gid",
    "t_start",
    "t_start_seas",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean"
  )
  missing <- setdiff(required, names(panel))
  if (length(missing) > 0) {
    stop(
      "build_climate_panel_for_tu: panel missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  gids <- unique(as.character(frame_tu$gid))
  start <- as.Date(horizon$start)
  end <- as.Date(horizon$end)
  start_warmup <- start - lubridate::days(round(warmup_months * 30.44))

  # Coverage check on the unfiltered panel — must touch the requested locations
  # AND span the requested horizon (with a small grace at the right edge).
  panel_loc <- panel %>%
    dplyr::filter(as.character(gid) %in% gids)
  if (nrow(panel_loc) == 0) {
    stop(
      "build_climate_panel_for_tu: panel has no rows for any TU location ",
      "(check the gid <-> panel join)."
    )
  }

  end_have <- max(panel_loc$t_start)
  if (end_have < end - 31) {
    stop(
      sprintf(
        "build_climate_panel_for_tu: panel ends %s, horizon end %s — gap > 1 month. ",
        end_have,
        end
      ),
      "Either narrow `simulator.horizon.end`, or rebuild the panel with extended coverage."
    )
  }
  end_month <- as.Date(format(end, "%Y-%m-01"))
  if (identical(end_mode, "closed") && end_have < end_month) {
    warning(
      sprintf(
        paste0(
          "build_climate_panel_for_tu: panel ends %s, before the closed horizon's last month %s: ",
          "that month is not simulated. Extend the panel or narrow simulator.horizon.end."
        ),
        end_have,
        format(end_month)
      ),
      call. = FALSE
    )
  }
  if (
    identical(end_mode, "open") &&
      end_have < seq(end_month, by = "-1 month", length.out = 2L)[2L]
  ) {
    stop(
      sprintf(
        paste0(
          "build_climate_panel_for_tu: panel ends %s, before the last month of the half-open ",
          "window [.., %s)."
        ),
        end_have,
        format(end_month)
      ),
      call. = FALSE
    )
  }
  start_have <- min(panel_loc$t_start)
  if (start_have > start) {
    stop(sprintf(
      "build_climate_panel_for_tu: panel starts %s, horizon start %s — panel begins after horizon.",
      start_have,
      start
    ))
  }

  out <- panel_loc %>%
    dplyr::filter(t_start >= start_warmup, t_start <= end)

  if (nrow(out) == 0) {
    stop(
      "build_climate_panel_for_tu: empty slice for TU; check horizon vs panel coverage."
    )
  }

  out
}

# ---------------------------------------------------------------------------
# Warm-up state — read lag columns from the panel, never recomputed
# ---------------------------------------------------------------------------

#' Initial state for the AR loop: u (TTLO), o (ADM1 lag), o_ct (country lag),
#' o_tu (TU lag), basehaz clock k, and active/mso for the
#' duration-aware path. Derived from the panel rows immediately before
#' `horizon$start`; falls back to zero for missing columns.
#'
#' @param has_dur logical; set TRUE by callers that will run the
#'   duration-aware branch of simulate_one_tu() (i.e. draws_dur is
#'   non-NULL). Purely diagnostic here -- it only gates the warning below --
#'   so its default FALSE leaves every existing call site's behaviour and
#'   return value unchanged (active/mso are always computed and returned,
#'   defaulting to 0 when panel_tu carries no such columns; unused, they are
#'   harmless for the legacy path).
build_warmup_state <- function(
  panel_tu,
  frame_tu,
  horizon,
  has_dur = FALSE,
  no_prior_onset = NULL
) {
  start <- as.Date(horizon$start)

  # Kernel speed-up: the same row selection on row indices (no dplyr), used
  # whenever gid and t_start carry no NA; identical warm-up state.
  ri <- .warm_row_index(panel_tu, frame_tu$gid, start)
  if (!is.null(ri)) {
    warm <- panel_tu[ri, , drop = FALSE]
  } else {
    pre <- panel_tu %>% dplyr::filter(t_start < start)
    warm <- if (nrow(pre) > 0) {
      pre %>%
        dplyr::group_by(gid) %>%
        dplyr::filter(t_start == max(t_start)) %>%
        dplyr::ungroup()
    } else {
      pre[0, , drop = FALSE]
    }

    # If we have no pre-horizon rows for some location, fall back to the first
    # post-horizon row (so the kernel can still launch in tests where no warm-up
    # tail was provided).
    if (nrow(warm) < nrow(frame_tu)) {
      missing_gids <- setdiff(frame_tu$gid, warm$gid)
      cand <- panel_tu %>% dplyr::filter(gid %in% missing_gids)
      if (nrow(cand) > 0) {
        first_post <- cand %>%
          dplyr::group_by(gid) %>%
          dplyr::filter(t_start == min(t_start)) %>%
          dplyr::ungroup()
        warm <- dplyr::bind_rows(warm, first_post)
      }
    }

    # Align rows with frame_tu order
    warm <- warm[match(frame_tu$gid, warm$gid), , drop = FALSE]
  }

  pick <- function(col, default = 0) {
    if (col %in% names(warm)) {
      v <- warm[[col]]
      v[is.na(v)] <- default
      v
    } else {
      rep(default, nrow(warm))
    }
  }

  # u: step 4 — initialise as days from the panel's factor-level TTLO column.
  ttlo_cat <- if ("time_to_last_outbreak_cat" %in% names(warm)) {
    warm[["time_to_last_outbreak_cat"]]
  } else {
    rep(NA_character_, nrow(warm))
  }
  u_days <- ttlo_init_from_panel(ttlo_cat)
  # "none" level: a unit with no recorded onset before the
  # launch (no_prior_onset, frame_tu order; ttlo_no_prior_onset_at()) whose panel
  # category is NA or [6,Inf] starts at the sentinel. NULL (default, and every
  # fit without the level) leaves u_days as above.
  if (!is.null(no_prior_onset)) {
    if (length(no_prior_onset) != nrow(frame_tu)) {
      stop(
        "build_warmup_state: no_prior_onset must have one value per frame_tu row."
      )
    }
    cat_chr <- as.character(ttlo_cat)
    seed <- no_prior_onset %in%
      TRUE &
      (is.na(cat_chr) | cat_chr == .TTLO_LEVELS[length(.TTLO_LEVELS)])
    u_days[seed] <- .TTLO_NONE_DAYS
  }

  # o_tu: locations with NA TU membership have ongoing_outbreak_tu_lag1 = NA
  # at fit time — preserve that NA so the spillover term is absent from the LP.
  o_tu_val <- as.integer(pick("ongoing_outbreak_tu_lag1", default = 0L))
  o_tu_val[is.na(frame_tu$tu)] <- NA_integer_

  # active/mso: duration-aware state. No historical source populates
  # these into panel_tu (the continuous climate panel) -- that handoff is
  # the smoothing cache's job (build_warmup_state()). Without it,
  # every location warm-starts as inactive; flag it loudly (not silently) so
  # a duration-aware run without a real warm-start is a visible choice, not
  # an accident.
  active0 <- as.integer(pick("active", default = 0L))
  mso0 <- as.integer(pick("mso", default = 0L))
  if (
    isTRUE(has_dur) &&
      !("active" %in% names(warm)) &&
      !("mso" %in% names(warm))
  ) {
    warning(
      "build_warmup_state: duration draws supplied (has_dur = TRUE) but ",
      "panel_tu carries no 'active'/'mso' warm-up columns -- every location ",
      "starts the duration-aware simulation as inactive (active = 0, ",
      "mso = 0). Expected for v1 runs with no smoothing-cache warm start ",
      "(Phase 9); if a warm start was expected, check that panel_tu / ",
      "warmup_panels_list actually carries 'active'/'mso'.",
      call. = FALSE
    )
  }

  list(
    u = u_days,
    o = as.integer(pick("ongoing_outbreak_lag1", 0)),
    o_ct = as.integer(pick("ongoing_outbreak_country_lag1", 0)),
    o_tu = o_tu_val,
    active = active0,
    mso = mso0,
    # Detection flag (covariate_feed: detected) -- only when the warm start
    # carries it (a detected-feed smoothing handoff). Consumers default it to
    # `active`: an observed-panel warm start is detected by construction.
    det = if ("det" %in% names(warm)) as.integer(pick("det", 0L)) else NULL
  )
}

#' ttlo_no_prior_onset_at
#'
#' TRUE for each gid with no recorded onset (`outbreak == 1` row of the
#' survival panel) before month `at`: the units whose TTLO is the "none" level
#' at a launch. A gid absent from the panel has no recorded
#' onset (TRUE).
# Row of panel_tu each frame gid's warm-up state is read from, exactly as the
# dplyr path of build_warmup_state() picks it: the gid's latest pre-horizon row
# (first of any tie, in panel order); when the pre-horizon rows number fewer
# than the frame (the path's own trigger), gids without one fall back to their
# earliest row. NULL (use the dplyr path) when gid or t_start carries NA.
.warm_row_index <- function(panel_tu, gids, start) {
  ts <- panel_tu$t_start
  g <- as.character(panel_tu$gid)
  if (anyNA(ts) || anyNA(g)) {
    return(NULL)
  }
  gids <- as.character(gids)
  tn <- as.numeric(ts)
  pre <- which(ts < start)
  idx <- rep(NA_integer_, length(gids))
  n_warm <- 0L
  if (length(pre)) {
    mx <- stats::ave(tn[pre], g[pre], FUN = max)
    cand <- pre[tn[pre] == mx]
    n_warm <- length(cand)
    idx <- cand[match(gids, g[cand])]
  }
  if (n_warm < length(gids)) {
    miss <- which(is.na(idx))
    post <- which(g %in% gids[miss])
    if (length(post)) {
      mn <- stats::ave(tn[post], g[post], FUN = min)
      cand2 <- post[tn[post] == mn]
      idx[miss] <- cand2[match(gids[miss], g[cand2])]
    }
  }
  idx
}

ttlo_no_prior_onset_at <- function(surv_dat, gids, at) {
  fo <- frame_first_onset(surv_dat, gids)
  is.na(fo) | fo >= as.Date(at)
}

#' frame_first_onset
#'
#' First recorded onset month (`outbreak == 1`) of each gid in the survival
#' panel, NA when none: the per-gid column `first_onset` the loaders put on the
#' location frame when the occurrence fit carries the TTLO level "none", from
#' which every kernel derives no_prior_onset at its own launch.
frame_first_onset <- function(surv_dat, gids) {
  on <- surv_dat[
    as.integer(surv_dat$outbreak) %in% 1L,
    c("gid", "t_start"),
    drop = FALSE
  ]
  fo <- tapply(as.Date(on$t_start), as.character(on$gid), min)
  # as.numeric(): tapply() returns a 1-d array (a `dim` attribute)
  as.Date(as.numeric(fo[as.character(gids)]), origin = "1970-01-01")
}

#' The no_prior_onset vector of a frame at `start` (build_warmup_state()):
#' NULL unless the run's occurrence draws carry the "none" level; then the
#' frame must carry `first_onset` (frame_first_onset(), set by the loaders).
.ttlo_no_prior_for_frame <- function(frame_tu, start, has_none) {
  if (!isTRUE(has_none)) {
    return(NULL)
  }
  if (!"first_onset" %in% names(frame_tu)) {
    stop(
      "the occurrence fit carries the TTLO level 'none' but the location frame has no ",
      "`first_onset` column (frame_first_onset(), set by the simulator / smoothing loaders).",
      call. = FALSE
    )
  }
  is.na(frame_tu$first_onset) | frame_tu$first_onset >= as.Date(start)
}

# ---------------------------------------------------------------------------
# Draws bundle — internal representation used by both extractors and the kernel
# ---------------------------------------------------------------------------
#
# Shape: list(
#   S            number of posterior draws
#   fe           S x p matrix of population-level coefs (canonical term names)
#   re_intercept list(region, country, adm1, adm2 or NULL) of S x n_level
#                matrices, columns = level ids
#   re_slope     per term, a list of S x n_level matrices
#                (e.g. re_slope$ongoing_outbreak_lag1$country)
#   basehaz      list(kind = factor | cyclic | seasonal | ms, + its draws)
#   phi          length-S monthly NegBin dispersion; NULL for occ
#   active_terms FE terms evaluated (colnames(fe) modulo TTLO unpacking))
# The kernel reads only this bundle, never a fit object.

#' Construct a deterministic synthetic draws bundle for testing.
#'
#' Builds a draws bundle with S = `S` posterior draws over the locations in
#' `frame_tu`. All RE intercepts / slopes default to zero; FE intercept and
#' slopes default to zero unless overridden by `fe_overrides`. The basehaz is
#' a flat constant `basehaz_const`.
#'
#' @param frame_tu output of build_location_frame() restricted to one TU.
#' @param S integer; number of fake posterior draws.
#' @param fe_overrides named numeric vector, e.g. c(`(Intercept)` = -3, mean_precip = .1).
#' @param re_slope_overrides named list of matrices to set non-zero RE slopes.
#' @param phi numeric scalar for the monthly (`mon`) NegBin dispersion.
#' @param basehaz_const numeric scalar; constant added to log-hazard.
make_draws_bundle_synthetic <- function(
  frame_tu,
  S = 4,
  fe_overrides = NULL,
  re_slope_overrides = NULL,
  phi = 1.0,
  basehaz_const = 0,
  active_terms = c(
    "(Intercept)",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1"
  )
) {
  regions <- unique(frame_tu$AFRO_region)
  countries <- unique(frame_tu$country)
  adm1s <- unique(frame_tu$adm1_unit)

  zero_mat <- function(cols) {
    m <- matrix(0, nrow = S, ncol = length(cols))
    colnames(m) <- cols
    m
  }

  fe <- zero_mat(active_terms)
  if (!is.null(fe_overrides)) {
    for (nm in names(fe_overrides)) {
      if (nm %in% colnames(fe)) fe[, nm] <- fe_overrides[[nm]]
    }
  }

  re_slope <- list()
  if (!is.null(re_slope_overrides)) {
    for (term in names(re_slope_overrides)) {
      re_slope[[term]] <- re_slope_overrides[[term]]
    }
  }

  list(
    S = S,
    fe = fe,
    re_intercept = list(
      region = zero_mat(regions),
      country = zero_mat(countries),
      adm1 = zero_mat(adm1s),
      adm2 = NULL
    ),
    re_slope = re_slope,
    basehaz = list(kind = "const", const = basehaz_const),
    phi = rep(phi, S),
    active_terms = active_terms
  )
}

# ---------------------------------------------------------------------------
# Backend extractors
# ---------------------------------------------------------------------------

#' Extract a draws bundle from a fitted occurrence model.
#'
#' v1 routes on `opt$occurrence$backend`. The brms path uses
#' `posterior::as_draws_df()` to gather FE (b_*) and group-level (r_*) draws.
#' Other backends raise an informative error and refer to the design note's
#' v1-scope decision.
#'
#' v1 simplification: when called on a real fit, the extractor returns at most
#' `max_draws` posterior draws (thinned uniformly) to keep memory bounded.
#'
#' A fit opt without `occurrence.backend` is refused: every
#' occurrence fit script attaches its resolved options (`attr(fit, "opt")`,
#' which carries the backend) and deploy_fits patches a fit that lacks them, so
#' a missing backend means the wrong options reached the extractor, not a
#' brms fit to guess.
extract_occurrence_draws <- function(fit, opt, max_draws = 200) {
  backend <- opt$occurrence$backend
  if (is.null(backend)) {
    stop(
      "extract_occurrence_draws: the fit's options carry no occurrence.backend ",
      "(attr(fit, \"opt\") missing, or not an occurrence fit's options). Re-attach ",
      "the fit's resolved options (deploy_prod_run.R --stage deploy_fits does so for ",
      "a deployed fit).",
      call. = FALSE
    )
  }
  switch(
    backend,
    brms = .extract_occ_draws_brms(fit, opt, max_draws),
    brms_cb = .extract_occ_draws_brms(fit, opt, max_draws),
    stop("extract_occurrence_draws: unknown occurrence backend: ", backend)
  )
}

#' The empirical baseline-rate anchor, as two fields on a draws BUNDLE.
#'
#' Same protocol as `offset_pop` / `offset_logoffset`
#' (`.extract_nb_count_draws_brms()` below): every mu block in every kernel reads
#' the flag off the bundle rather than branching on which model it is scoring,
#' so an un-anchored fit scores exactly as it did before this change.
#'
#' Deliberately keyed on the `opt` block name and NOT on `anchor_spec()`:
#' `utils/anchor.R` depends on `paths.R` / `provenance.R` / `occ_data.R`, and
#' `simulator_core.R` is sourced by callers (e.g. `07b_simulator_spatpomp.R`)
#' that do not carry that closure. Mirrors `.OFFSET_SPEC` (`utils/config.R`),
#' which is the formula-side registry of the same fact; the two are kept in
#' step by `analysis/tests/test_base_rate_table.R`, not by a source-order
#' coupling.
#'
#' SIZE IS DELIBERATELY ABSENT: the size track is not simulated, and
#' `size.offset_baserate` was a PRIOR gate
#' only (config_defaults.yml, `.offset_terms()`), so it must not be read here.
.BUNDLE_ANCHOR_COL <- list(
  occurrence = "base_rate_occ",
  duration = "base_rate_dur",
  monthly = "base_rate_mon"
)

.bundle_anchor_flags <- function(opt, block) {
  col <- .BUNDLE_ANCHOR_COL[[block]]
  on <- !is.null(col) && isTRUE((opt[[block]] %||% list())$offset_baserate)
  list(offset_baserate = on, anchor_col = if (on) col else NA_character_)
}

#' `log(anchor)` for the rows `idx` of a location frame, or the additive
#' identity `0` when the bundle's fit did not carry the anchor offset.
#'
#' Every mu block adds this unconditionally (`lp <- lp + .anchor_log_offset(...)`),
#' so an un-anchored bundle costs one comparison and changes nothing. The three
#' `stop()`s are deliberate -- a missing anchor column would otherwise silently
#' rescale every simulated prediction by the whole baseline-rate range, which is
#' exactly the failure mode the anchor exists to prevent.
#'
#' @param bundle a draws bundle carrying `offset_baserate` / `anchor_col`.
#' @param frame  the location frame (needs the `anchor_col` column).
#' @param idx    integer row index into `frame`; `NULL` means all rows.
#' @param what   label for error messages (the track name).
.anchor_log_offset <- function(bundle, frame, idx = NULL, what = "kernel") {
  if (!isTRUE(bundle$offset_baserate %||% FALSE)) {
    return(0)
  }
  col <- bundle$anchor_col
  if (is.null(col) || length(col) != 1L || is.na(col)) {
    stop(
      what,
      ": draws bundle has offset_baserate = TRUE but no `anchor_col`. ",
      "Re-extract the bundle with the current extractor."
    )
  }
  if (!(col %in% names(frame))) {
    stop(
      what,
      ": the fit carries offset(log(",
      col,
      ")) but the location ",
      "frame has no `",
      col,
      "` column. Pass it through ",
      "build_location_frame(base_rate_table_paths = list(...)) -- see ",
      "utils/anchor.R persist_base_rate_table()."
    )
  }
  v <- if (is.null(idx)) frame[[col]] else frame[[col]][idx]
  if (any(!is.finite(v)) || any(v <= 0)) {
    stop(
      what,
      ": `",
      col,
      "` must be finite and strictly positive on every ",
      "location (",
      sum(!is.finite(v) | v <= 0),
      " bad value(s)). A ",
      "non-positive anchor means the rate table did not resolve."
    )
  }
  log(v)
}

# ---------------------------------------------------------------------------
# Onset-month exposure correction
# ---------------------------------------------------------------------------
#
# The mon and dur fits use a partial-month exposure offset in an episode's
# first and last month (dur_data.R .expand_active_months()); the replay
# kernels use log(days_in_month). The fitted per-day rate is right, so only
# the replay is corrected:
#   * onset month (mso == 1): exposure x .ONSET_EXPOSURE_FRAC, the expected
#     active fraction under a uniform onset day ((D+1)/(2D) ~ 0.5; mu is
#     linear in exposure, so no Jensen term);
#   * cessation month: only under simulator.cessation_exposure = "half", where
#     every kernel draws cessation BEFORE the monthly mean (otherwise after
#     mu) and gives a ceasing unit .CESSATION_EXPOSURE_FRAC of a month (1/4
#     if it also started that month). The cessation draw still reads the
#     post-update lags, the monthly mean the pre-update ones; the C kernel's
#     RNG stream and the state convention smoothing_core.R decodes are as
#     under "full".
# mon and dur ONLY: occurrence shares the offset object but its exposure is a
# genuine at-risk month, so the correction is a per-call-site adjustment.
.ONSET_EXPOSURE_FRAC <- 0.5

#' Additive log-offset correction for an episode's onset month.
#'
#' @param mso integer vector of months-since-onset (1-based: 1 on the onset
#'   month, as set by every kernel's active/mso promotion block). Values of 0
#'   (inactive unit) and >= 2 (interior/later month) contribute nothing.
#' @param use_logoffset logical; when FALSE the caller did not add a log-offset
#'   at all, so there is nothing to correct and the additive identity is
#'   returned. Lets the mon call sites stay a one-line change under their
#'   `offset_logoffset` bundle gate.
#' @return numeric vector the same length as `mso` (or scalar 0).
.onset_log_offset_adj <- function(mso, use_logoffset = TRUE) {
  if (!isTRUE(use_logoffset)) {
    return(0)
  }
  hit <- as.integer(mso) == 1L
  out <- numeric(length(hit))
  out[which(hit)] <- log(.ONSET_EXPOSURE_FRAC)
  out[is.na(hit)] <- NA_real_
  out
}

# Cessation-month exposure fraction (simulator.cessation_exposure: half): the
# panel's mean exposure / days_in_month at cessation is 0.505;
# the same theory-derived 1/2 as the onset month. See .ONSET_EXPOSURE_FRAC.
.CESSATION_EXPOSURE_FRAC <- 0.5

#' Additive log-offset correction for an episode's cessation month.
#'
#' @param ceased logical vector: the unit was active this month and its
#'   cessation was drawn this step.
#' @param use_logoffset as in .onset_log_offset_adj(): FALSE (no log-offset
#'   added by the caller) returns the additive identity.
#' @return numeric vector the same length as `ceased` (or scalar 0).
.cessation_log_offset_adj <- function(ceased, use_logoffset = TRUE) {
  if (!isTRUE(use_logoffset)) {
    return(0)
  }
  hit <- as.logical(ceased)
  out <- numeric(length(hit))
  out[which(hit)] <- log(.CESSATION_EXPOSURE_FRAC)
  out[is.na(hit)] <- NA_real_
  out
}

#' Replay-vs-fit exposure on the monthly fit panel (the measurement in the
#' .ONSET_EXPOSURE_FRAC comment). Per panel row the replay exposure is
#' days_in_month x 1/2 at onset (mso == 1) x 1/2 on the last month of an
#' episode that ended (`ended`, not censored) under `half`; the fit exposure is
#' pmax(exposure, 1), as in log_offset. Diagnostic only, no kernel calls it.
#' @param panel monthly fit panel (load_monthly_data()'s source RDS): month,
#'   mso, n_months, ended, exposure, cases_month.
#' @return tibble: mode (none / full / half), case_weighted (sum of
#'   cases x replay / fit over sum of cases) and person_days (sum replay over
#'   sum fit).
replay_exposure_ratios <- function(panel) {
  D <- as.numeric(lubridate::days_in_month(panel$month))
  fit <- pmax(as.numeric(panel$exposure), 1)
  ons <- as.integer(panel$mso) == 1L
  cess <- as.integer(panel$mso) == as.integer(panel$n_months) &
    panel$ended %in% c(TRUE, 1L)
  w <- as.numeric(panel$cases_month)
  ok <- is.finite(w)
  one <- function(fac) {
    rep <- D * fac
    c(
      case_weighted = sum(w[ok] * rep[ok] / fit[ok]) / sum(w[ok]),
      person_days = sum(rep) / sum(fit)
    )
  }
  r <- rbind(
    none = one(rep(1, length(D))),
    full = one(ifelse(ons, .ONSET_EXPOSURE_FRAC, 1)),
    half = one(
      ifelse(ons, .ONSET_EXPOSURE_FRAC, 1) *
        ifelse(cess, .CESSATION_EXPOSURE_FRAC, 1)
    )
  )
  tibble::tibble(
    mode = rownames(r),
    case_weighted = unname(r[, 1]),
    person_days = unname(r[, 2])
  )
}

#' Normalise a brms occurrence term name to the bare covariate name the kernel
#' uses. brms encodes a logical/two-level predictor (e.g. the spillover lags
#' `ongoing_outbreak_lag1`) with a trailing factor-level suffix `TRUE`
#' (reference = FALSE), so the FE column is `ongoing_outbreak_lag1TRUE` and the
#' group-level term inside `r_<group>[level, ...]` is likewise suffixed. The
#' kernel refers to the bare `ongoing_outbreak_lag1`, so strip the suffix.
#' Continuous terms (`enso_lag5`) and factor terms (`time_to_last_outbreak12`)
#' never carry a `TRUE` suffix and are returned unchanged.
.normalize_occ_term <- function(term) {
  is_logical_lag <- grepl("^ongoing_outbreak", term) & grepl("TRUE$", term)
  ifelse(is_logical_lag, sub("TRUE$", "", term), term)
}

#' Extract brms smooth-term coefficient matrices from a posterior draws data frame.
#'
#' Detects columns matching `s_<smooth_label>[k]`, groups by smooth label, and
#' returns a named list of S × K matrices (one per smooth term). Returns an
#' empty list when no smooth columns are present.
.extract_smooth_coef_from_draws <- function(draws) {
  cols <- names(draws)
  s_pat <- "^s_([^\\[]+)\\[(\\d+)\\]$"
  s_cols <- stringr::str_subset(cols, s_pat)
  if (length(s_cols) == 0L) {
    return(list())
  }

  mat <- stringr::str_match(s_cols, s_pat)
  labels <- unique(mat[, 2])

  lapply(stats::setNames(labels, labels), function(lbl) {
    sc <- s_cols[mat[, 2] == lbl]
    idx <- as.integer(stringr::str_match(sc, "\\[(\\d+)\\]$")[, 2])
    as.matrix(draws[, sc[order(idx)], drop = FALSE])
  })
}

#' Build a smooth prediction function for `s(t_start_seas, bs="cc", k=8)` that
#' applies the same constraint absorption brms uses internally.
#'
#' @param fit a brmsfit object; used only to read `fit$data$t_start_seas`.
#' @param K   number of free spline coefficients (= ncol of the coefficient
#'            matrix extracted from the draws; 6 for k=8 cyclic cubic).
#' @return a function `function(new_t) → N×K matrix`, or NULL on failure.
.make_t_start_seas_pred_fn <- function(fit, K) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    warning(
      ".make_t_start_seas_pred_fn: mgcv required; seasonal spline disabled."
    )
    return(NULL)
  }
  t_train <- fit$data[["t_start_seas"]]
  if (is.null(t_train)) {
    warning(
      ".make_t_start_seas_pred_fn: t_start_seas absent from fit$data; seasonal spline disabled."
    )
    return(NULL)
  }
  # k = K + 1 for cyclic cubic: the raw basis has k-1 columns (cyclic constraint
  # built into the basis), absorb.cons removes 1 more for identifiability → K = k-2.
  # So k_knots = K + 2.  For K=6 (draws), k_knots=8 matches the original fit k=8.
  k_knots <- K + 2L

  # Build both raw and absorbed smooth objects: mgcv::Predict.matrix returns
  # the RAW basis, while brms passes Stan Zs_1_1 = X_raw %*% T_full, with
  #   T_abs  = qr.Q(qr(t(C)), complete = TRUE)[, -(1:nrow(C))]  (absorption)
  #   T_DR   = eigen(sm$S[[1]])$vectors %*% diag(1 / sqrt(values)), columns
  #            reversed (Demmler-Reinsch rotation)
  #   T_full = T_abs %*% T_DR_rev
  mk_sm <- function(absorb) {
    tryCatch(
      mgcv::smoothCon(
        mgcv::s(t_start_seas, bs = "cc", k = k_knots),
        data = data.frame(t_start_seas = t_train),
        absorb.cons = absorb
      )[[1]],
      error = function(e) {
        warning(sprintf(
          ".make_t_start_seas_pred_fn: smoothCon failed (%s); seasonal spline disabled.",
          conditionMessage(e)
        ))
        NULL
      }
    )
  }
  sm <- mk_sm(TRUE)
  sm_raw <- mk_sm(FALSE)
  if (is.null(sm) || is.null(sm_raw)) {
    return(NULL)
  }

  # Verify the absorbed basis has exactly K columns.
  if (ncol(sm$X) != K) {
    warning(sprintf(
      ".make_t_start_seas_pred_fn: ncol(sm$X)=%d != K=%d; seasonal spline disabled.",
      ncol(sm$X),
      K
    ))
    return(NULL)
  }

  # Step 1: constraint-absorption matrix T_abs (K_raw × K).
  T_abs <- if (!is.null(sm_raw$C) && nrow(sm_raw$C) > 0L) {
    nc <- nrow(sm_raw$C)
    Q <- qr.Q(qr(t(sm_raw$C)), complete = TRUE)
    Q[, -(1:nc), drop = FALSE]
  } else {
    diag(K)
  }

  # Step 2: Demmler-Reinsch rotation T_DR_rev (K × K) from the penalty matrix.
  # sm$S[[1]] is the K × K penalty in the absorbed parameterisation.
  T_DR_rev <- tryCatch(
    {
      eig <- eigen(sm$S[[1]], symmetric = TRUE)
      # brms reverses column order: col j of T_DR corresponds to col K+1-j of Zs_1_1
      rev_idx <- K:1L
      eig$vectors[, rev_idx, drop = FALSE] %*%
        diag(1 / sqrt(eig$values[rev_idx]), nrow = K, ncol = K)
    },
    error = function(e) {
      warning(sprintf(
        ".make_t_start_seas_pred_fn: eigen(sm$S[[1]]) failed (%s); using identity.",
        conditionMessage(e)
      ))
      diag(K)
    }
  )

  # Combined transformation: X_raw %*% T_full = Zs_1_1_new
  T_full <- T_abs %*% T_DR_rev # K_raw × K

  function(new_t) {
    X_raw <- mgcv::Predict.matrix(
      sm_raw,
      data.frame(t_start_seas = as.numeric(new_t))
    )
    X_raw %*% T_full # N × K  (= Zs_1_1_new, column means ≈ 0 by brms centering)
  }
}

.extract_occ_draws_brms <- function(fit, opt, max_draws = 200) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior package required for brms extractor")
  }

  draws <- posterior::as_draws_df(fit)
  S_full <- nrow(draws)
  if (S_full == 0) {
    stop(".extract_occ_draws_brms: empty posterior")
  }
  S <- min(S_full, max_draws)
  if (S < S_full) {
    take <- seq(1, S_full, length.out = S) |> round() |> as.integer()
    draws <- draws[take, , drop = FALSE]
  }
  bundle <- .occ_bundle_from_draws_df(draws)
  # Attached HERE rather than in .occ_bundle_from_draws_df(), which has no
  # `opt` (it is unit-tested against synthetic draws frames).
  bundle[c("offset_baserate", "anchor_col")] <-
    .bundle_anchor_flags(opt, "occurrence")
  # Attach smooth prediction functions (one per smooth label found in draws).
  # Currently the only expected smooth is s(t_start_seas, bs="cc") → label "st_start_seas_1".
  if (length(bundle$smooth_coef) > 0L) {
    bundle$smooth_pred_fn <- lapply(bundle$smooth_coef, function(coef_mat) {
      .make_t_start_seas_pred_fn(fit, ncol(coef_mat))
    })
    names(bundle$smooth_pred_fn) <- names(bundle$smooth_coef)
  }
  bundle
}

#' Parse the fixed-effect (`b_*`) and nested random-effect (`r_*`) columns of a
#' (thinned) `posterior::as_draws_df` data frame into the kernel's matrices.
#' Shared by the occurrence, duration and NegBin-count (monthly) extractors so
#' all reproduce the brms LP identically; pulled out of `.occ_bundle_from_draws_df`
#' so it can be reused by the count extractor (whose deepest nesting group is `...:ADM1:ADM2` rather than
#' the occurrence track's `...:ADM1:gid`). Returns `list(fe, re_intercept,
#' re_slope)`; callers add backend-specific fields (smooth_coef, phi, ...).
.parse_brms_re_draws <- function(draws) {
  S <- nrow(draws)
  cols <- names(draws)
  b_cols <- stringr::str_subset(cols, "^b_")
  r_cols <- stringr::str_subset(cols, "^r_")

  fe <- as.matrix(draws[, b_cols, drop = FALSE])
  colnames(fe) <- .normalize_occ_term(stringr::str_remove(b_cols, "^b_"))

  # Group r_* columns by group prefix `r_<group>[level,term]`.
  # We split each column into (group, level, term).
  parse_r <- function(nm) {
    # nm = "r_<group>[<level>,<term>]"
    m <- stringr::str_match(nm, "^r_([^\\[]+)\\[([^,]+),(.+)\\]$")
    if (any(is.na(m[1, ]))) {
      return(NULL)
    }
    list(group = m[1, 2], level = m[1, 3], term = m[1, 4])
  }

  # Lookup table: group prefix -> simulator-bucket name. The occurrence track
  # nests by `gid` at the deepest level; the size track nests by `ADM2`
  # (= gid). Both map to the "adm2" bucket — the matching
  # composite level keys are reconstructed identically by re_level_key().
  # Bare-country entries handle grouped fits (drop_region=TRUE) whose intercept
  # and spillover REs use `country/ADM1/gid` nesting without an `AFRO_region`
  # prefix. ENSO RE slopes in grouped
  # brms_cb fits still use `AFRO_region/...` nesting and are captured by the
  # AFRO entries; the dual-key lookup in evaluate_lp/.ns1 resolves both.
  group_bucket <- function(g) {
    switch(
      g,
      "AFRO_region" = "region",
      "AFRO_region:country" = "country",
      "AFRO_region:country:ADM1" = "adm1",
      "AFRO_region:country:ADM1:gid" = "adm2",
      "AFRO_region:country:ADM1:ADM2" = "adm2",
      # Grouped fits: bare country nesting (drop_region=TRUE)
      "country" = "country",
      "country:ADM1" = "adm1",
      "country:ADM1:gid" = "adm2",
      "country:ADM1:ADM2" = "adm2",
      # re_adm2_leaf_scope: adm2_only (utils/config.R) -- the leaf is its
      # own term `(0 + is_adm2 | ...:ADM1:gid_leaf)`, whose levels use the
      # same `<...>_<adm1_unit>_<gid>` form re_level_key() builds.
      "AFRO_region:country:ADM1:gid_leaf" = "adm2",
      "country:ADM1:gid_leaf" = "adm2",
      NA_character_
    )
  }

  re_intercept <- list(region = NULL, country = NULL, adm1 = NULL, adm2 = NULL)
  re_slope <- list()

  for (nm in r_cols) {
    p <- parse_r(nm)
    if (is.null(p)) {
      next
    }
    bucket <- group_bucket(p$group)
    if (is.na(bucket)) {
      next
    }

    # Encode level. brms levels for crossed nested factors use "_" joins;
    # we keep them as-is here — the kernel uses the trailing token via
    # frame_tu, see below.
    level_key <- p$level
    term <- .normalize_occ_term(p$term)
    # adm2_only leaf: `is_adm2` is 1 on every unit that HAS a leaf level, so
    # its coefficient is that unit's leaf intercept. ADM1-only units' keys
    # never match the shared "..._none" level, so they get 0 -- correct, they
    # have no leaf.
    if (identical(term, "is_adm2") && identical(bucket, "adm2")) {
      term <- "Intercept"
    }

    if (term == "Intercept") {
      if (is.null(re_intercept[[bucket]])) {
        re_intercept[[bucket]] <- matrix(numeric(0), nrow = S, ncol = 0)
      }
      re_intercept[[bucket]] <- cbind(
        re_intercept[[bucket]],
        structure(
          as.numeric(draws[[nm]]),
          dim = c(S, 1),
          dimnames = list(NULL, level_key)
        )
      )
    } else {
      if (is.null(re_slope[[term]])) {
        re_slope[[term]] <- list()
      }
      if (is.null(re_slope[[term]][[bucket]])) {
        re_slope[[term]][[bucket]] <- matrix(numeric(0), nrow = S, ncol = 0)
      }
      re_slope[[term]][[bucket]] <- cbind(
        re_slope[[term]][[bucket]],
        structure(
          as.numeric(draws[[nm]]),
          dim = c(S, 1),
          dimnames = list(NULL, level_key)
        )
      )
    }
  }

  list(fe = fe, re_intercept = re_intercept, re_slope = re_slope)
}

#' Extract per-draw phi RE matrices from a fitted brms NegBin count model
#' (the monthly `mon` track).
#'
#' Uses brms::ranef(fit, summary=FALSE, dpar="shape") to avoid relying on the
#' brms distributional-parameter column-naming conventions in the flat draws df.
#' Returns a list with four elements:
#'   fe      — numeric[S] shape FE intercept (log scale, b_shape_Intercept)
#'   region  — matrix[S × R] or NULL
#'   country — matrix[S × C] or NULL
#'   adm1    — matrix[S × A] or NULL
#' Column names of each matrix match the re_level_key() convention (underscores,
#' dots instead of spaces), ready for match() lookup in compute_phi_vec().
.extract_phi_re_draws <- function(fit, draws, take) {
  S <- nrow(draws)
  fe <- if ("b_shape_Intercept" %in% names(draws)) {
    as.numeric(draws[["b_shape_Intercept"]])
  } else {
    rep(0, S)
  }

  phi_ranef <- brms::ranef(fit, summary = FALSE, dpar = "shape")

  ranef_bucket <- c(
    "AFRO_region" = "region",
    "AFRO_region:country" = "country",
    "AFRO_region:country:ADM1" = "adm1",
    # Grouped (drop_region) fits: build_size_formula re-roots the shape RE chain
    # at country (`shape ~ 1 + (1 | country/ADM1)`, utils/formula.R), so brms
    # names the shape ranef groups `country` / `country:ADM1` with NO AFRO prefix.
    # Without these entries the phi RE matrices were silently dropped and every
    # location collapsed to exp(b_shape_Intercept) — mirrors the bare-country
    # entries in .parse_brms_re_draws()'s group_bucket() for the mean LP.
    "country" = "country",
    "country:ADM1" = "adm1"
  )
  result <- list(fe = fe, region = NULL, country = NULL, adm1 = NULL)

  for (key in names(phi_ranef)) {
    bucket <- ranef_bucket[key]
    if (is.na(bucket)) {
      next
    }
    arr <- phi_ranef[[key]][take, , 1, drop = FALSE] # thin: [S × n_groups × 1]
    grp_labels <- gsub(" ", ".", dimnames(arr)[[2]], fixed = TRUE)
    mat <- matrix(arr, nrow = S, dimnames = list(NULL, grp_labels))
    result[[bucket]] <- mat
  }
  result
}

#' Parse a (thinned) `posterior::as_draws_df` data frame into the kernel's draws
#' bundle. Pulled out of `.extract_occ_draws_brms` so it can be unit-tested with
#' a synthetic draws frame carrying real brms column names, including the
#' suffix/enso naming.
.occ_bundle_from_draws_df <- function(draws) {
  parsed <- .parse_brms_re_draws(draws)
  smooth_coef <- .extract_smooth_coef_from_draws(draws)

  .assert_lp_terms_consumed(
    colnames(parsed$fe),
    "occurrence",
    label = "extract_occurrence_draws"
  )

  list(
    S = nrow(draws),
    fe = parsed$fe,
    re_intercept = parsed$re_intercept,
    re_slope = parsed$re_slope,
    basehaz = list(kind = "brms_seasonal_in_lp"),
    phi = NULL,
    active_terms = colnames(parsed$fe),
    smooth_coef = smooth_coef, # named list of S×K matrices, one per smooth term
    smooth_pred_fn = NULL # filled by .extract_occ_draws_brms (needs fit$data)
  )
}

# ---------------------------------------------------------------------------
# Fixed-effect coverage gate
# ---------------------------------------------------------------------------

#' Regexes for every fixed-effect column the kernels actually score.
#'
#' `evaluate_lp()` consumes a HARD-CODED whitelist of term names; anything a fit
#' carries that is not on it is silently dropped, and the forecast then replays
#' a different model than the one that was fitted. That failure is invisible in
#' every diagnostic (the run completes, the numbers are plausible) -- it is how
#' `log_pop_c` shipped in four deployed size fits of an earlier run while the kernel
#' applied a hard-coded population elasticity of 1.
#'
#' This table is the machine-checkable statement of what the kernels cover, and
#' `.assert_lp_terms_consumed()` below turns "a new covariate reached a fit but
#' not the kernel" from a silent bias into a hard error at extraction time --
#' the same class of bug as the CLI-route and prior-builder gaps recorded in
#' CLAUDE.md, caught by the same kind of gate.
.LP_CONSUMED_PATTERNS <- c(
  "^\\(Intercept\\)$",
  "^Intercept$",
  # linear_terms in evaluate_lp()
  "^mean_precip$",
  "^dev_precip$",
  "^mean_temp_daymean$",
  "^dev_temp_daymean$",
  "^enso_(lag|lead)[0-9]+$",
  "_cb[0-9]+$",
  "^sin_month$",
  "^cos_month$",
  "^sin2_month$",
  "^cos2_month$",
  "^admin_levelADM2$",
  "^log_pop_c$",
  # spillover lags (brms `TRUE` suffix already stripped by .normalize_occ_term)
  "^ongoing_outbreak_lag1$",
  "^ongoing_outbreak_country_lag1$",
  "^ongoing_outbreak_tu_lag1$",
  # distance-weighted neighbour contagion (dynamic recompute; see simulate_one_tu)
  # Literal twin of paths.R::nbr_covariate_pattern(): this vector is built at
  # source time, so it cannot call the helper. test_nbr_naming.R asserts the
  # two stay identical.
  "^nbr_lag1_r[0-9]+(_rn)?(_pg|_part[a-z0-9]+)?$",
  # TTLO factor levels
  "^time_to_last_outbreak"
)

# Additional per-track consumers OUTSIDE evaluate_lp().
.LP_CONSUMED_PATTERNS_BY_TRACK <- list(
  duration = "^mso_band" # basehaz_dur_eval()
)

# Non-`mu` distributional parameters (brms `b_<dpar>_<term>`). Not part of the
# linear predictor at all -- `shape` travels as the bundle's `phi`.
.LP_DPAR_PATTERN <- paste0(
  "^(shape|sigma|phi|nu|zi|hu|zoi|coi|kappa|beta|xi|alpha|ndt|bias|disc|",
  "quantile)_"
)

#' Hard-fail when a fit carries a fixed effect no kernel will score.
#'
#' @param fe_names colnames of the extracted FE matrix.
#' @param track one of "occurrence" / "monthly" / "duration"; selects
#'   the extra consumers above.
#' @param label human name for the error message.
.assert_lp_terms_consumed <- function(fe_names, track, label = track) {
  pats <- c(
    .LP_CONSUMED_PATTERNS,
    .LP_CONSUMED_PATTERNS_BY_TRACK[[track]],
    .LP_DPAR_PATTERN
  )
  covered <- Reduce(
    `|`,
    lapply(pats, function(p) grepl(p, fe_names)),
    init = rep(FALSE, length(fe_names))
  )
  orphan <- fe_names[!covered]
  if (length(orphan) == 0L) {
    return(invisible(TRUE))
  }
  stop(
    "[",
    label,
    "] the fit carries ",
    length(orphan),
    " fixed effect(s) no ",
    "simulator/smoothing kernel scores: ",
    paste(orphan, collapse = ", "),
    ".\n  evaluate_lp() consumes a fixed whitelist of term names, so these ",
    "coefficients would be SILENTLY DROPPED at replay -- the forecast would ",
    "run a different model than the one that was fitted, with no warning ",
    "and no diagnostic that could detect it.\n  Either add the term to ",
    "evaluate_lp()'s `linear_terms` (and to .LP_CONSUMED_PATTERNS), or ",
    "re-fit without it.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Zero-filled climate gate
# ---------------------------------------------------------------------------

#' The four local-climate panel columns that 03_prepare_sim_panel_projection.R
#' zero-fills past the climate DB's last month (`CLIMATE_DB_MAX`, 2023-12): a
#' row with all four exactly 0 is a filled row, not an observed climate month.
.CLIMATE_PANEL_COLS <- c(
  "mean_precip",
  "dev_precip",
  "mean_temp_daymean",
  "dev_temp_daymean"
)

#' Fixed-effect names that read a local-climate panel column: the linear terms
#' themselves, their lag columns, and a crossbasis built on one of them.
.CLIMATE_FE_PATTERN <- "^(mean|dev)_(precip|temp_daymean)(_lag[0-9]+|_cb[0-9]+)?$"

#' Union of the fixed-effect names of any number of draws bundles.
#'
#' @param ... draws bundles (`extract_*_draws()` output, with `$fe`), NULL, or
#'   named lists of bundles (the `*_by_group` lists); NULLs are skipped.
#' @return character vector of unique FE column names.
.draws_fe_names <- function(...) {
  out <- character(0)
  add <- function(x) {
    if (is.null(x)) {
      return(invisible())
    }
    # [[ ]], not $: a group list could hold a group whose name starts "fe".
    if (!is.null(x[["fe"]])) {
      out <<- c(out, colnames(x[["fe"]]))
      return(invisible())
    }
    if (is.list(x)) {
      for (b in x) {
        add(b)
      }
    }
    invisible()
  }
  for (x in list(...)) {
    add(x)
  }
  unique(out)
}

#' The panel months a forward run reads, for the zero-filled climate gate:
#' from the earliest of the horizon start, the smoothing handoff
#' (`enso_sim_months_start`) and the hybrid-warmup start, to `horizon.end`
#' (open when unset). Shared by the simulations.R loaders and 07b v1's own.
#'
#' @param sim the `simulator:` block.
#' @param enso_sim_months_start optional Date (smoothing handoff month).
#' @return length-2 Date (start, end); NA = open.
.sim_panel_window <- function(sim, enso_sim_months_start = NULL) {
  hz <- sim[["horizon"]] %||% list()
  h_start <- if (is.null(hz[["start"]])) as.Date(NA) else as.Date(hz[["start"]])
  starts <- h_start
  if (!is.null(enso_sim_months_start)) {
    starts <- c(starts, as.Date(enso_sim_months_start))
  }
  wu <- sim[["warmup"]] %||% list()
  if (isTRUE(wu[["use_sim_trajectories"]]) && !is.na(h_start)) {
    starts <- c(
      starts,
      h_start %m-% months(as.integer(wu[["warmup_months"]] %||% 12L))
    )
  }
  start <- if (all(is.na(starts))) as.Date(NA) else min(starts, na.rm = TRUE)
  end <- if (is.null(hz[["end"]])) as.Date(NA) else as.Date(hz[["end"]])
  c(start, end)
}

#' Hard-fail when a fit with a local-climate fixed effect would be replayed on
#' zero-filled climate rows.
#'
#' The simulator / smoothing panels carry the four climate columns only up to
#' the climate DB's last month and set them to exactly 0 after it
#' (`03_prepare_sim_panel_projection.R`). The forward-fills in the kernels only
#' cover months with NO panel row, so the zero-filled rows pass through: a fit
#' with `mean_temp_daymean` would be forecast at a mean temperature of 0 with no
#' error. Inert for every fit whose fixed effects carry no climate term (every
#' deployed fit: the deployed term set is in analysis/model_terms_deployed.yml).
#'
#' @param fe_names character, fixed-effect names of every deployed track /
#'   group (`.draws_fe_names()`).
#' @param panel climate panel (columns `t_start` and the four climate columns).
#' @param window length-2 Date (start, end) of the months the run reads; an NA
#'   end is open. NULL checks the whole panel.
#' @param label human name for the error message.
#' @return invisible(TRUE) when the gate passes.
.assert_no_zero_filled_climate <- function(
  fe_names,
  panel,
  window = NULL,
  label = "load_artefacts"
) {
  clim_fe <- grep(.CLIMATE_FE_PATTERN, as.character(fe_names), value = TRUE)
  if (length(clim_fe) == 0L) {
    return(invisible(TRUE))
  }

  missing_cols <- setdiff(.CLIMATE_PANEL_COLS, names(panel))
  if (length(missing_cols) > 0L) {
    stop(
      "[",
      label,
      "] the fits carry local-climate fixed effect(s) ",
      paste(clim_fe, collapse = ", "),
      " but the climate panel lacks ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  in_win <- rep(TRUE, nrow(panel))
  if (!is.null(window)) {
    if (!"t_start" %in% names(panel)) {
      stop(
        "[",
        label,
        "] climate panel has no t_start column to window on.",
        call. = FALSE
      )
    }
    w <- as.Date(window)
    ts <- as.Date(panel$t_start)
    if (!is.na(w[1L])) {
      in_win <- in_win & ts >= w[1L]
    }
    if (length(w) > 1L && !is.na(w[2L])) in_win <- in_win & ts <= w[2L]
  }
  zero <- in_win
  for (cc in .CLIMATE_PANEL_COLS) {
    v <- panel[[cc]]
    zero <- zero & !is.na(v) & v == 0
  }
  if (!any(zero)) {
    return(invisible(TRUE))
  }

  zm <- if ("t_start" %in% names(panel)) {
    range(as.Date(panel$t_start[zero]))
  } else {
    NULL
  }
  stop(
    "[",
    label,
    "] the fits carry local-climate fixed effect(s) ",
    paste(clim_fe, collapse = ", "),
    ", but ",
    sum(zero),
    " climate-panel ",
    "row(s) in the run window have all four climate columns exactly 0",
    if (!is.null(zm)) sprintf(" (%s to %s)", zm[1L], zm[2L]) else "",
    ".\n  These are the zero-filled months past the climate DB's last month ",
    "(03_prepare_sim_panel_projection.R, CLIMATE_DB_MAX): replaying them ",
    "would forecast at a rainfall / temperature of 0 with no other error.\n  ",
    "Either fit without local-climate terms, or build the panel with real ",
    "(or climatological) climate for every month the run reads.",
    call. = FALSE
  )
}

#' Canonical admin level implied by a flattened GADM gid.
#'
#' `CMR_5_1` -> `"ADM1"`, `CMR_5_1_1` -> `"ADM2"`. GADM gid depth IS the
#' canonical admin level -- that is exactly what `apply_canonical_admin_level()`
#' (P3, utils/provenance.R) re-keys `admin_level` to before
#' `.add_log_pop_centred()` runs. Used on both sides of the P4 re-centring (the
#' training-side centre recovery in `.recover_log_pop_centres()` and the
#' frame-side lookup in `evaluate_lp()`) so the two cannot disagree: a frame
#' unit and the training rows for the same gid must resolve to the same centre
#' or the elasticity is applied around the wrong pivot.
#'
#' @param gid character vector of flattened gids.
#' @return character vector of `"ADM<k>"` labels, same length.
.admin_level_from_gid <- function(gid) {
  paste0("ADM", lengths(strsplit(as.character(gid), "_", fixed = TRUE)) - 2L)
}

#' The admin_level column a CANONICAL-KEYED consumer should read off a frame.
#'
#' `build_location_frame()` carries two level columns once `gid_pop_path` is
#' supplied: `admin_level` (the corpus label the outbreak database asserts, which
#' the occurrence/duration anchors and `drop_overlapping_adm1_from_frame()` key off)
#' and `admin_level_canonical` (GADM depth, which is what the size and monthly
#' fits were re-keyed to by `apply_canonical_admin_level()` before any level-keyed
#' quantity -- anchor cell, `log_pop_c` centre, `admin_levelADM2` dummy -- was
#' computed). Consumers that must reproduce the size/monthly FIT-side keying call
#' this; consumers that must reproduce the corpus keying read `admin_level`
#' directly. Falls back to the corpus label, then to NA, so a frame built without
#' `gid_pop_path` behaves exactly as it did before the column existed.
.frame_admin_level_canonical <- function(frame) {
  if ("admin_level_canonical" %in% names(frame)) {
    return(as.character(frame$admin_level_canonical))
  }
  if ("admin_level" %in% names(frame)) {
    return(as.character(frame$admin_level))
  }
  rep(NA_character_, NROW(frame))
}

#' Recover the P4 `log_pop_c` centring constants from a fitted model.
#'
#' `log_pop_c = log(pop) - centre[admin_level]` (`.add_log_pop_centred()`,
#' utils/siz_data.R) is the deviation term `pop_elasticity: free` fits in
#' place of the hard-coded-1 population offset. The simulator MUST re-centre its
#' own frame units on the TRAINING centres -- centring on the frame's own mean
#' would move the pivot and silently rescale every prediction by
#' `exp((beta - 1) * shift)`.
#'
#' Three sources, in order:
#'   1. `attr(fit, "<track>_estimand")$log_pop_centres` -- written by the fit
#'      scripts from `get_size_estimand(size_dat)` / `get_monthly_estimand(mon_dat)`.
#'      Authoritative. Both attribute names are tried because the two tracks
#'      stamp under DELIBERATELY distinct names (`.ESTIMAND_SPEC`,
#'      utils/provenance.R): a mon stamp and a size stamp describe different
#'      denominators on different panels, so one shared name would let a mon
#'      artefact satisfy a size assertion.
#'   2. Recovered exactly from `fit$data`: the training frame carries both
#'      `log_pop_c` and (because `offset(log(pop))` makes `pop` a model
#'      variable) `pop`, and `log(pop) - log_pop_c == centre[admin_level]` by
#'      construction. This is what makes the fix apply to fits already on disk,
#'      which predate (1).
#'   3. Hard error. A silent 0 would be exactly the bug this function exists to
#'      close: it would apply the fitted elasticity around the wrong pivot.
#'
#' @return named numeric vector, one centre per admin_level (name `"ALL"` when
#'   the fit was centred globally, matching `.add_log_pop_centred()`'s fallback).
.recover_log_pop_centres <- function(fit, tag = "monthly") {
  for (nm in c("size_estimand", "mon_estimand")) {
    stamp <- attr(fit, nm, exact = TRUE)
    if (!is.null(stamp[["log_pop_centres"]])) {
      ctr <- stamp[["log_pop_centres"]]
      return(stats::setNames(as.numeric(ctr), names(ctr)))
    }
  }

  d <- tryCatch(fit$data, error = function(e) NULL)
  if (is.null(d) || !all(c("log_pop_c", "pop") %in% names(d))) {
    stop(
      ".recover_log_pop_centres: the ",
      tag,
      " fit has a `log_pop_c` fixed ",
      "effect (pop_elasticity != 'fixed') but neither a ",
      "`size_estimand`/`mon_estimand` ",
      "attribute carrying `log_pop_centres` nor a `fit$data` with both ",
      "`log_pop_c` and `pop`. Refusing to guess the centring constants -- ",
      "applying the fitted elasticity around the wrong pivot rescales every ",
      "prediction by exp((beta - 1) * shift). Re-fit (the fit scripts now ",
      "attach the stamp) or drop the elasticity term.",
      call. = FALSE
    )
  }

  delta <- log(as.numeric(d$pop)) - as.numeric(d$log_pop_c)
  # brms keeps only the variables the FORMULA mentions, and `admin_level` is
  # not one of them under `level_factor: false` -- so the level key is rebuilt
  # from the gid (carried as `ADM2` by the size formula's leaf RE, `d$ADM2 <-
  # d$gid` in read_size_data()). GADM gid depth IS the canonical admin level --
  # that is exactly what apply_canonical_admin_level() (P3) re-keys `admin_level`
  # to before `.add_log_pop_centred()` runs -- so this reproduces the training
  # grouping rather than approximating it. Verified against all four deployed
  # size fits of an earlier run: zero within-level spread in log(pop) - log_pop_c.
  gid_col <- intersect(c("gid", "ADM2"), names(d))
  lvl <- if ("admin_level" %in% names(d)) {
    as.character(d$admin_level)
  } else if (length(gid_col) > 0L) {
    .admin_level_from_gid(d[[gid_col[1L]]])
  } else {
    rep("ALL", nrow(d))
  }
  ok <- is.finite(delta)
  if (!any(ok)) {
    stop(
      ".recover_log_pop_centres: no finite log(pop) - log_pop_c rows in the ",
      tag,
      " fit's training data.",
      call. = FALSE
    )
  }
  ctr <- tapply(delta[ok], lvl[ok], function(v) {
    if (diff(range(v)) > 1e-8) {
      stop(
        ".recover_log_pop_centres: log(pop) - log_pop_c is not constant ",
        "within an admin_level (spread ",
        format(diff(range(v))),
        "). The ",
        "fit was not centred by .add_log_pop_centred(); refusing to guess.",
        call. = FALSE
      )
    }
    v[1L]
  })
  ctr <- stats::setNames(as.numeric(ctr), names(ctr))
  message(sprintf(
    "[%s] log_pop_c centres recovered from fit$data: %s",
    tag,
    paste(sprintf("%s=%.4f", names(ctr), ctr), collapse = " ")
  ))
  ctr
}

#' Extract a draws bundle from a fitted brms negative-binomial count GLMM.
#'
#' Full posterior draws (thinned to `max_draws`), carrying the model's nested
#' random effects and per-draw NB dispersion. Sole caller:
#' `extract_monthly_draws()` (`block = "monthly"`). `block` has no default.
#'
#' @param block name of the `opt` sub-block holding the track's settings
#'   (required). The posterior parsing is track-agnostic -- only `phi_re` and
#'   the offset toggles are read per track.
.extract_nb_count_draws_brms <- function(fit, opt, max_draws = 200L, block) {
  if (missing(block) || !is.character(block) || length(block) != 1L) {
    stop(
      ".extract_nb_count_draws_brms: `block` (the opt sub-block name, ",
      "e.g. \"monthly\") is required."
    )
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms required for the NB count extractor")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior package required for the NB count extractor")
  }

  draws <- posterior::as_draws_df(fit)
  S_full <- nrow(draws)
  if (S_full == 0) {
    stop(".extract_nb_count_draws_brms: empty posterior")
  }
  S <- min(S_full, max_draws)
  if (S < S_full) {
    take <- seq(1, S_full, length.out = S) |> round() |> as.integer()
    draws <- draws[take, , drop = FALSE]
  } else {
    take <- seq_len(S_full)
  }

  parsed <- .parse_brms_re_draws(draws)

  # Per-draw NB dispersion. When phi_re == "constant" (the default), brms stores
  # a scalar shape per draw in the `shape` column — kept as a plain numeric
  # vector for backward compat. When phi_re != "constant", phi becomes a list
  # with per-level matrices; compute_phi_vec() handles both cases in the kernel.
  blk <- opt[[block]] %||% list()
  phi_re <- blk$phi_re %||% "constant"
  phi <- if (phi_re == "constant") {
    if ("shape" %in% names(draws)) as.numeric(draws[["shape"]]) else rep(1, S)
  } else {
    .extract_phi_re_draws(fit, draws, take)
  }

  anchor <- .bundle_anchor_flags(opt, block)

  # P4 (`pop_elasticity: free` / `free_by_level`). When the fit carries the
  # `log_pop_c` deviation term the kernel must rebuild that column for every
  # frame unit, around the TRAINING centres. Absent the term this is NULL and
  # every existing bundle is unchanged.
  log_pop_centres <- if ("log_pop_c" %in% colnames(parsed$fe)) {
    .recover_log_pop_centres(fit, tag = block)
  } else {
    NULL
  }

  .assert_lp_terms_consumed(
    colnames(parsed$fe),
    block,
    label = paste0("extract_", block, "_draws")
  )

  list(
    S = S,
    fe = parsed$fe,
    re_intercept = parsed$re_intercept,
    re_slope = parsed$re_slope,
    basehaz = list(kind = "none"),
    phi = phi,
    phi_re = phi_re,
    active_terms = colnames(parsed$fe),
    log_pop_centres = log_pop_centres,
    # See `.bundle_anchor_flags()`.
    offset_baserate = anchor$offset_baserate,
    anchor_col = anchor$anchor_col,
    # Offsets travel on the BUNDLE, not on the track: every mu block in the
    # kernels reads these flags rather than branching on which model it is
    # scoring.
    offset_pop = isTRUE(blk$offset_pop %||% TRUE),
    # log_offset (= log(pmax(exposure, 1))) is an UNCONDITIONAL offset of the
    # monthly formula (build_monthly_formula_parts(), utils/config.R);
    # extract_monthly_draws() pins it TRUE.
    offset_logoffset = isTRUE(blk$offset_logoffset %||% FALSE)
  )
}


# ---------------------------------------------------------------------------
# Monthly case-count (`mon`) draws extraction
# ---------------------------------------------------------------------------

#' Extract a draws bundle from a fitted monthly case-count model.
#'
#' A thin dispatcher over `.extract_nb_count_draws_brms(block = "monthly")`
#' (a negative-binomial GLMM with the shared nested-RE machinery). Shaped on
#' extract_duration_draws().
#'
#' Two deliberate `stop()`s rather than silent approximations, following the
#' duration extractor's `basehaz_duration == "ms"` precedent in
#' extract_duration_draws():
#'
#'   * `prev_cases: true` adds a `log_cases_prev` covariate whose value is the
#'     previous month's *drawn* case count. Carrying it would need a 12th
#'     spatPomp state slot, which would bump `.SM_CACHE_SCHEMA_VERSION`
#'     (spatpomp_core.R) and invalidate every existing smoothing handoff cache.
#'     It is `false` by default (config_defaults.yml) and no shipped config sets
#'     it.
#'   * a non-brms monthly backend has no extractor at all.
#'
#' @param fit       fitted brms monthly model (artefact kind "mon_fit").
#' @param opt       config list carrying the `monthly` sub-block.
#' @param max_draws posterior thinning target.
#' @return a draws bundle consumed by the kernels' mon mu blocks, with
#'   `offset_logoffset = TRUE`.
extract_monthly_draws <- function(fit, opt, max_draws = 200L) {
  mon <- opt$monthly %||% list()
  if (isTRUE(mon$prev_cases)) {
    stop(
      "extract_monthly_draws: monthly.prev_cases: true is not yet ",
      "implemented in the simulator/smoothing kernels (the log_cases_prev ",
      "covariate would require an extra latent state slot). Re-fit with ",
      "monthly.prev_cases: false, or drop `monthly_fit` from the simulator ",
      "config."
    )
  }
  # Fallback = the fitting default (get_default_options(): brms_cb); brms and
  # brms_cb share the extractor.
  backend <- mon$backend %||% "brms_cb"
  out <- switch(
    backend,
    brms = .extract_nb_count_draws_brms(
      fit,
      opt,
      max_draws = max_draws,
      block = "monthly"
    ),
    brms_cb = .extract_nb_count_draws_brms(
      fit,
      opt,
      max_draws = max_draws,
      block = "monthly"
    ),
    stop("extract_monthly_draws: unknown monthly backend: ", backend)
  )

  # The monthly offset set is a property of the FORMULA
  # (formula.R::build_monthly_formula_parts()), not a tunable, so
  # it is pinned here rather than read from the config:
  #   offset(log(pop))       -- toggled by monthly.offset_pop (already handled)
  #   offset(log_offset)     -- UNCONDITIONAL
  # test_monthly_offset_wellformed.R is the gate on this line.
  #
  # The third offset, offset(log(base_rate_mon)), IS a real monthly toggle
  # (monthly.offset_baserate) and is therefore read from the config rather than
  # pinned -- .extract_nb_count_draws_brms(block = "monthly") has already set
  # `offset_baserate` / `anchor_col` on `out` via .bundle_anchor_flags().
  out$offset_logoffset <- TRUE
  out
}

# ---------------------------------------------------------------------------
# Duration (cessation-hazard) draws extraction
# ---------------------------------------------------------------------------

#' Extract a draws bundle from a fitted duration (discrete-time cessation
#' hazard) model. Dispatch shaped on extract_occurrence_draws().
extract_duration_draws <- function(fit, opt, max_draws = 200) {
  backend <- opt$duration$backend %||% "brms"
  switch(
    backend,
    brms = .extract_dur_draws_brms(fit, opt, max_draws),
    brms_cb = .extract_dur_draws_brms(fit, opt, max_draws),
    stop("extract_duration_draws: unknown duration backend: ", backend)
  )
}

#' Extract a draws bundle from a fitted duration model (brms backend).
#'
#' Structurally identical to the occurrence extractor -- same Poisson +
#' log-offset likelihood (analysis/04_fit_dur_brms.R), same nested RE
#' machinery (deepest nesting ADM1 by default; under `duration.admin_level:
#' ADM1-ADM2` + `re_adm2_leaf_scope: adm2_only` the `country:ADM1:gid_leaf`
#' leaf lands in the adm2 bucket exactly as monthly's does, and evaluate_lp()
#' adds it) -- so this reuses .parse_brms_re_draws() / .normalize_occ_term() /
#' .extract_smooth_coef_from_draws() / .make_t_start_seas_pred_fn() verbatim.
#'
#' Unlike .occ_bundle_from_draws_df() (which always hardcodes
#' basehaz = list(kind = "brms_seasonal_in_lp"), because the occurrence
#' kernel only ever ships fits with basehaz_poisson = "seasonal"), the
#' duration bundle's basehaz is read from opt$duration -- the duration
#' track's PRIMARY clock defaults to basehaz_duration = "factor" (a real FE
#' factor term on mso_band), which basehaz_dur_eval() evaluates directly
#' rather than delegating to evaluate_lp()'s spline machinery. The
#' SECONDARY (seasonal, t_start_seas) clock defaults to "none" and, when
#' enabled as "cyclic"/"seasonal", reuses the exact same smooth_coef /
#' smooth_pred_fn / evaluate_lp() spline path the occurrence track uses --
#' that machinery is generic on xrow$t_start_seas and was never
#' occurrence-specific. "factor"/"ms" on the secondary clock are not yet
#' supported by the v1 kernel (the same gap the occurrence track already has
#' for its own basehaz_poisson = "factor"/"ms" -- never implemented because
#' every occurrence fit shipped to the kernel used "seasonal") and error
#' loudly rather than silently evaluating to a flat hazard.
.extract_dur_draws_brms <- function(fit, opt, max_draws = 200) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior package required for brms extractor")
  }

  draws <- posterior::as_draws_df(fit)
  S_full <- nrow(draws)
  if (S_full == 0) {
    stop(".extract_dur_draws_brms: empty posterior")
  }
  S <- min(S_full, max_draws)
  if (S < S_full) {
    take <- seq(1, S_full, length.out = S) |> round() |> as.integer()
    draws <- draws[take, , drop = FALSE]
  }

  parsed <- .parse_brms_re_draws(draws)

  # The kernel draws no episode-level size, so a duration fit conditioned on
  # log(size_cases) (E1, `duration.covariates: "log(size_cases)": yes`) has no
  # value to score that term with. Hard error rather than a silent zero.
  if ("logsize_cases" %in% colnames(parsed$fe)) {
    stop(
      "extract_duration_draws: the duration fit carries a `logsize_cases` ",
      "fixed effect (duration.covariates \"log(size_cases)\"), but the size ",
      "track is retired from the simulator (2026-10-06) and no episode size ",
      "is drawn any more. Re-fit the duration model without ",
      "\"log(size_cases)\".",
      call. = FALSE
    )
  }

  smooth_coef <- .extract_smooth_coef_from_draws(draws)
  dur <- opt$duration %||% list()

  .assert_lp_terms_consumed(
    colnames(parsed$fe),
    "duration",
    label = "extract_duration_draws"
  )

  bh_kind <- switch(
    dur$basehaz_duration %||% "factor",
    none = "none",
    factor = "factor",
    ms = stop(
      "extract_duration_draws: duration.basehaz_duration = 'ms' ",
      "is not yet implemented in the v1 kernel (no runtime ",
      "M-spline basis evaluator for mso -- only 'factor' and ",
      "'none' are supported). See basehaz_dur_eval()."
    ),
    stop(
      "extract_duration_draws: unknown duration.basehaz_duration: ",
      dur$basehaz_duration
    )
  )
  basehaz <- list(
    kind = bh_kind,
    breaks = as.numeric(dur$mso_breaks %||% .SIM_DUR_MSO_BREAKS)
  )

  seas_kind <- switch(
    dur$basehaz_poisson %||% "none",
    none = "none",
    cyclic = "brms_seasonal_in_lp",
    seasonal = "brms_seasonal_in_lp",
    factor = stop(
      "extract_duration_draws: duration.basehaz_poisson = ",
      "'factor' is not yet implemented in the v1 kernel (same ",
      "gap as the occurrence track's own basehaz_poisson = ",
      "'factor', which no shipped fit has ever exercised). ",
      "Only 'none' and 'cyclic'/'seasonal' are supported."
    ),
    ms = stop(
      "extract_duration_draws: duration.basehaz_poisson = 'ms' ",
      "is not yet implemented in the v1 kernel."
    ),
    stop(
      "extract_duration_draws: unknown duration.basehaz_poisson: ",
      dur$basehaz_poisson
    )
  )

  bundle <- list(
    S = nrow(draws),
    fe = parsed$fe,
    re_intercept = parsed$re_intercept,
    re_slope = parsed$re_slope,
    basehaz = basehaz,
    phi = NULL,
    active_terms = colnames(parsed$fe),
    smooth_coef = smooth_coef,
    smooth_pred_fn = NULL
  )
  bundle[c("offset_baserate", "anchor_col")] <-
    .bundle_anchor_flags(opt, "duration")
  if (
    identical(seas_kind, "brms_seasonal_in_lp") &&
      length(bundle$smooth_coef) > 0L
  ) {
    bundle$smooth_pred_fn <- lapply(bundle$smooth_coef, function(coef_mat) {
      .make_t_start_seas_pred_fn(fit, ncol(coef_mat))
    })
    names(bundle$smooth_pred_fn) <- names(bundle$smooth_coef)
  }
  bundle
}

# ---------------------------------------------------------------------------
# Hazard / size LP evaluators (pure)
# ---------------------------------------------------------------------------

#' Reconstruct the brms grouping-factor level keys that label the columns of the
#' extracted nested intercept/slope RE matrices.
#'
#' A `/`-nested grouping such as `AFRO_region/country/ADM1/gid` expands to the
#' interaction factors `AFRO_region`, `AFRO_region:country`,
#' `AFRO_region:country:ADM1` and `AFRO_region:country:ADM1:gid`. brms labels
#' each level by pasting the *ancestor* factor values with `_` and replacing
#' spaces with `.` (e.g. region "Central Africa" -> "Central.Africa", so an
#' adm2 column is "Central.Africa_AGO_AGO_4_AGO_4_1"). The bare `frame$gid`
#' (`AGO_4_1`) therefore never matches the composite colnames; reconstruct the
#' composite from the same source columns brms used.
#'
#' @param frame frame (or frame_tu) with columns AFRO_region, country,
#'   adm1_unit, gid.
#' @param bucket one of "region", "country", "adm1", "adm2".
#' @param drop_region logical — when TRUE, omit the AFRO_region prefix so the
#'   composite keys match the colnames produced by region-dropped (grouped)
#'   brms fits, whose RE nesting is `country/ADM1(/gid)` only. The "region"
#'   bucket then returns a length-N vector of `NA` (no region RE exists in a
#'   grouped fit, so any lookup against it is a guaranteed miss → 0 contribution).
#'   Default FALSE preserves the full-Africa `AFRO_region/country/...` keys.
#' @return character vector of composite level keys, one per frame row.
re_level_key <- function(frame, bucket, drop_region = FALSE) {
  if (isTRUE(drop_region)) {
    return(switch(
      bucket,
      region = rep(NA_character_, nrow(frame)),
      country = frame$country,
      adm1 = paste(frame$country, frame$adm1_unit, sep = "_"),
      adm2 = paste(frame$country, frame$adm1_unit, frame$gid, sep = "_"),
      stop("re_level_key: unknown bucket: ", bucket)
    ))
  }
  reg <- gsub(" ", ".", frame$AFRO_region, fixed = TRUE)
  switch(
    bucket,
    region = reg,
    country = paste(reg, frame$country, sep = "_"),
    adm1 = paste(reg, frame$country, frame$adm1_unit, sep = "_"),
    adm2 = paste(reg, frame$country, frame$adm1_unit, frame$gid, sep = "_"),
    stop("re_level_key: unknown bucket: ", bucket)
  )
}

# ---------------------------------------------------------------------------
# Linear-predictor lookup plan
# ---------------------------------------------------------------------------
# lp_plan() resolves evaluate_lp()'s random-effect lookups (step- and
# replicate-invariant) once per (bundle, frame rows, drop_region); .lp_coefs()
# builds per-location coefficient vectors once per draw (one-slot cache).
# Bit-identity contract: same arithmetic and term order as the per-call code
# (tests/test_simulator_kernel_golden.R, HARD). A plan is bound to its bundle
# and frame ROWS; evaluate_lp() without one builds a throwaway plan.

#' Build the replicate- and step-invariant lookup plan of evaluate_lp().
#'
#' @param draws       draws bundle (occ / dur / mon).
#' @param frame_tu    the location rows that will be scored (RE keys, admin level).
#' @param drop_region as evaluate_lp().
#' @return an environment (it carries the one-slot per-draw coefficient cache).
lp_plan <- function(draws, frame_tu, drop_region = FALSE) {
  N <- nrow(frame_tu)
  # names() of the ROW as evaluate_lp() always took it (drop = TRUE): a
  # one-column fe matrix yields an unnamed scalar, so such a bundle scores none
  # of its terms -- kept exactly (synthetic single-term test bundles rely on it).
  fen <- names(draws$fe[1L, , drop = TRUE])

  key_cache <- list()
  keys_of <- function(bucket) {
    k <- key_cache[[bucket]]
    if (is.null(k)) {
      k <- list(
        primary = re_level_key(frame_tu, bucket, drop_region),
        alt = re_level_key(frame_tu, bucket, !drop_region)
      )
      key_cache[[bucket]] <<- k
    }
    k
  }
  # match_keys() of the per-call code: primary format first, alternate format
  # for the locations that miss.
  idx_of <- function(mat, bucket) {
    if (is.null(mat) || ncol(mat) == 0) {
      return(NULL)
    }
    keys <- keys_of(bucket)
    cn <- colnames(mat)
    idx <- match(keys$primary, cn)
    miss <- is.na(idx)
    if (any(miss)) {
      idx[miss] <- match(keys$alt[miss], cn)
    }
    idx
  }
  buckets <- c("region", "country", "adm1", "adm2")

  int_idx <- lapply(stats::setNames(buckets, buckets), function(b) {
    idx_of(draws$re_intercept[[b]], b)
  })
  slope_idx <- list()
  slope_terms <- names(draws$re_slope)
  for (term in slope_terms) {
    rst <- draws$re_slope[[term]]
    if (is.null(rst)) {
      next
    }
    slope_idx[[term]] <- lapply(stats::setNames(buckets, buckets), function(b) {
      idx_of(rst[[b]], b)
    })
  }

  cb_terms <- grep("_cb[0-9]+$", fen, value = TRUE)
  linear_terms <- c(
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    grep("^enso_(lag|lead)[0-9]+$", fen, value = TRUE),
    grep(nbr_covariate_pattern(), fen, value = TRUE),
    cb_terms,
    "sin_month",
    "cos_month",
    "sin2_month",
    "cos2_month",
    "admin_levelADM2",
    "log_pop_c"
  )
  linear_terms <- linear_terms[linear_terms %in% fen]
  ttlo_cols <- grep("^time_to_last_outbreak", fen, value = TRUE)

  plan <- new.env(parent = emptyenv())
  plan$N <- N
  plan$fen <- fen
  plan$drop_region <- drop_region
  plan$int_idx <- int_idx
  plan$slope_idx <- slope_idx
  plan$linear_terms <- linear_terms
  plan$nbr_cols <- grep(nbr_covariate_pattern(), fen, value = TRUE)
  plan$spill_terms <- intersect(
    c(
      "ongoing_outbreak_lag1",
      "ongoing_outbreak_country_lag1",
      "ongoing_outbreak_tu_lag1"
    ),
    fen
  )
  plan$ttlo_cols <- ttlo_cols
  plan$ttlo_levels <- vapply(ttlo_cols, ttlo_col_to_level, character(1))
  plan$ttlo_has_none <- .TTLO_NONE_COL %in% ttlo_cols
  plan$has_ref_slope <- !is.null(draws$re_slope[["time_to_last_outbreak01"]])
  plan$smooth_names <- if (
    length(draws$smooth_pred_fn) > 0L && length(draws$smooth_coef) > 0L
  ) {
    names(draws$smooth_pred_fn)
  } else {
    character(0)
  }
  plan$log_pop_centre <- NULL # lazily: only when a call needs log_pop_c
  plan$frame_tu <- frame_tu # for the lazy log_pop_c centre
  plan$s_cached <- NA_integer_
  plan$coefs <- NULL
  plan
}

#' Per-location coefficient vectors of one posterior draw (one-slot cache).
#' Arithmetic mirrors the per-call code it replaced term for term.
.lp_coefs <- function(plan, draws, s) {
  if (identical(plan$s_cached, s)) {
    return(plan$coefs)
  }
  N <- plan$N
  fe <- draws$fe[s, , drop = TRUE]

  # nested_slope(): sum of the term's RE slope over every nested level present.
  nested_slope <- function(term) {
    rst <- draws$re_slope[[term]]
    out <- rep(0, N)
    if (is.null(rst)) {
      return(out)
    }
    si <- plan$slope_idx[[term]]
    lk <- function(mat, idx) {
      if (is.null(idx)) {
        return(rep(0, N))
      }
      ifelse(is.na(idx), 0, mat[s, idx])
    }
    out +
      lk(rst$region, si$region) +
      lk(rst$country, si$country) +
      lk(rst$adm1, si$adm1) +
      lk(rst$adm2, si$adm2)
  }

  base <- rep(0, N)
  if ("(Intercept)" %in% plan$fen) {
    base <- base + fe[["(Intercept)"]]
  }
  if ("Intercept" %in% plan$fen) {
    base <- base + fe[["Intercept"]]
  }
  for (b in c("region", "country", "adm1", "adm2")) {
    idx <- plan$int_idx[[b]]
    if (is.null(idx)) {
      next
    }
    have <- !is.na(idx)
    if (any(have)) {
      base[have] <- base[have] + draws$re_intercept[[b]][s, idx[have]]
    }
  }

  lin <- list()
  for (term in plan$linear_terms) {
    lin[[term]] <- fe[[term]] + nested_slope(term)
  }

  spill <- list()
  for (term in plan$spill_terms) {
    false_key <- paste0(term, "FALSE")
    spill[[term]] <- list(
      true = fe[[term]] + nested_slope(term),
      false = if (!is.null(draws$re_slope[[false_key]])) {
        nested_slope(false_key)
      } else {
        NULL
      }
    )
  }

  ttlo <- list()
  for (col_name in plan$ttlo_cols) {
    ttlo[[col_name]] <- list(fe = fe[[col_name]], re = nested_slope(col_name))
  }
  ttlo_ref <- if (plan$has_ref_slope) {
    nested_slope("time_to_last_outbreak01")
  } else {
    NULL
  }

  smooth <- list()
  for (sm_nm in plan$smooth_names) {
    coef_mat <- draws$smooth_coef[[sm_nm]]
    smooth[[sm_nm]] <- if (is.null(coef_mat)) NULL else coef_mat[s, ]
  }

  plan$coefs <- list(
    base = base,
    lin = lin,
    spill = spill,
    ttlo = ttlo,
    ttlo_ref = ttlo_ref,
    smooth = smooth
  )
  plan$s_cached <- s
  plan$coefs
}

# P4 log_pop_c centre per location (training centres, keyed by the frame's
# canonical admin level); errors exactly as the per-call code did.
.lp_log_pop_centre <- function(draws, frame_tu, N) {
  ctr <- draws$log_pop_centres
  if (is.null(ctr) || length(ctr) == 0L) {
    stop(
      "evaluate_lp: the fit has a `log_pop_c` fixed effect but the draws ",
      "bundle carries no `log_pop_centres`. See ",
      ".recover_log_pop_centres() -- there is no safe default (a 0 pivot ",
      "would rescale every prediction by exp((beta - 1) * mean(log pop))).",
      call. = FALSE
    )
  }
  if (!("pop" %in% names(frame_tu))) {
    stop(
      "evaluate_lp: `log_pop_c` needs `pop` on the location frame.",
      call. = FALSE
    )
  }
  # Keyed on the frame's CANONICAL admin level (P3 re-keying; see
  # .frame_admin_level_canonical()), falling back to the gid depth rule
  # .recover_log_pop_centres() uses on the training side, then "resolve at the
  # coarsest available match" (.resolve_frame_anchor()).
  lvl <- .frame_admin_level_canonical(frame_tu)
  if (length(lvl) != N) {
    lvl <- rep(lvl, length.out = N)
  }
  if (anyNA(lvl) && "gid" %in% names(frame_tu)) {
    na_lvl <- is.na(lvl)
    lvl[na_lvl] <- .admin_level_from_gid(frame_tu$gid[na_lvl])
  }
  lvl[is.na(lvl)] <- "ALL"
  lvl <- .coarsen_admin_level_to_available(lvl, names(ctr))
  c_i <- unname(ctr[lvl])
  if (length(ctr) == 1L) {
    c_i[is.na(c_i)] <- unname(ctr[[1L]])
  }
  if (anyNA(c_i)) {
    stop(
      "evaluate_lp: no training `log_pop_c` centre for admin_level(s) ",
      paste(unique(lvl[is.na(c_i)]), collapse = ", "),
      " (fit centres: ",
      paste(names(ctr), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  c_i
}

#' Evaluate the linear predictor for one posterior draw `s`, all locations
#' in TU at time-step `m`. Inputs:
#'   draws       — draws bundle
#'   s           — posterior draw index
#'   xrow        — named list of covariate vectors length N (mean_precip, ...)
#'   o, o_ct, o_tu — integer vectors length N; NA for NA-TU locations in o_tu
#'   u_days      — numeric vector length N; days since the last outbreak,
#'     binned into the fit's `time_to_last_outbreak*` levels when present.
#'     NA → reference level (0 contribution).
#'   frame_tu    — location frame restricted to TU (for level lookups).
#'   basehaz_val — numeric vector length N, added to the LP (per-location step).
#'   drop_ct_tu  — TRUE skips the o_ct/o_tu LP contribution (own `o` is never
#'     affected); threads `simulator.drop_country_tu_spillover`.
#'   nbr         — per-step neighbour-contagion pressure, or NULL.
#'   plan        — lp_plan(draws, frame_tu, drop_region) built for exactly this
#'     bundle and these rows, or NULL (a throwaway plan is built per call).
#'   smooth_X    — named list (by smooth) of the seasonal-spline design at this
#'     step's t_start_seas, or NULL (computed from draws$smooth_pred_fn).
#'
#' Terms consumed: intercepts and nested intercept REs (region -> country ->
#' adm1 -> adm2; grouped brms_cb fits mix bare-country and AFRO-prefixed keys,
#' so the alternate key is a fallback); every linear term (climate, ENSO, nbr,
#' `_cb<k>` columns, mon harmonics, admin_levelADM2, log_pop_c centred on the
#' TRAINING centres) as (FE + nested RE slopes) x value; the spillover lags
#' with their FALSE-level slopes (NA lag -> 0, as at fit time); the TTLO
#' levels incl. the reference bin's RE slope and the optional "none" level;
#' the seasonal spline.
#'
#' Returns a numeric vector length N of log-rates (log λ).
evaluate_lp <- function(
  draws,
  s,
  xrow,
  o,
  o_ct,
  o_tu,
  u_days,
  frame_tu,
  basehaz_val = 0,
  drop_region = FALSE,
  drop_ct_tu = FALSE,
  nbr = NULL,
  plan = NULL,
  smooth_X = NULL
) {
  if (is.null(plan)) {
    plan <- lp_plan(draws, frame_tu, drop_region)
  } else if (plan$N != nrow(frame_tu)) {
    stop(
      "evaluate_lp: the lookup plan was built for ",
      plan$N,
      " rows but frame_tu has ",
      nrow(frame_tu),
      ".",
      call. = FALSE
    )
  }
  N <- plan$N
  cf <- .lp_coefs(plan, draws, s)
  lp <- cf$base

  if ("log_pop_c" %in% plan$linear_terms && is.null(xrow[["log_pop_c"]])) {
    if (is.null(plan$log_pop_centre)) {
      plan$log_pop_centre <- .lp_log_pop_centre(draws, plan$frame_tu, N)
    }
    xrow[["log_pop_c"]] <- log(as.numeric(frame_tu$pop)) - plan$log_pop_centre
  }
  # Neighbour-contagion value is supplied dynamically per step: written onto the
  # detected nbr column so the linear loop scores it like any covariate.
  if (!is.null(nbr)) {
    for (ncol_ in plan$nbr_cols) {
      xrow[[ncol_]] <- nbr
    }
  }

  for (term in plan$linear_terms) {
    xv <- xrow[[term]]
    if (is.null(xv)) {
      next
    }
    if (anyNA(xv)) {
      xv[is.na(xv)] <- 0
    }
    lp <- lp + cf$lin[[term]] * xv
  }

  add_spillover <- function(lp, term, lag_vec) {
    sp <- cf$spill[[term]]
    if (is.null(sp)) {
      return(lp)
    }
    contrib_true <- sp$true * lag_vec
    contrib_true[is.na(lag_vec)] <- 0 # NA lag → term absent at fit time → 0
    # brms stores RE slopes for all factor levels, including FALSE: applied to
    # rows where lag_vec == 0 (NA treated as 0).
    contrib_false <- if (!is.null(sp$false)) {
      w <- 1L - as.integer(lag_vec)
      w[is.na(w)] <- 0L
      sp$false * w
    } else {
      rep(0, N)
    }
    lp + contrib_true + contrib_false
  }
  lp <- add_spillover(lp, "ongoing_outbreak_lag1", o)
  if (!isTRUE(drop_ct_tu)) {
    lp <- add_spillover(lp, "ongoing_outbreak_country_lag1", o_ct)
    lp <- add_spillover(lp, "ongoing_outbreak_tu_lag1", o_tu)
  }

  if (length(plan$ttlo_cols) > 0) {
    bins <- bin_ttlo(u_days, has_none = plan$ttlo_has_none) # NA for NA input
    for (j in seq_along(plan$ttlo_cols)) {
      in_bin <- !is.na(bins) & (bins == plan$ttlo_levels[[j]])
      if (!any(in_bin)) {
        next
      }
      tt <- cf$ttlo[[plan$ttlo_cols[[j]]]]
      lp[in_bin] <- lp[in_bin] + tt$fe + tt$re[in_bin]
    }
    # Reference bin [0,1): no FE term but its RE slope is stored and applied.
    if (!is.null(cf$ttlo_ref)) {
      in_ref <- !is.na(bins) & (bins == .TTLO_LEVELS[1L]) # "[0,1)"
      if (any(in_ref)) lp[in_ref] <- lp[in_ref] + cf$ttlo_ref[in_ref]
    }
  }

  # Seasonal spline (s(t_start_seas, bs = "cc")): coefficients in
  # draws$smooth_coef, design from draws$smooth_pred_fn (or the caller's
  # precomputed smooth_X for this step).
  if (length(plan$smooth_names) > 0L) {
    t_ss <- xrow[["t_start_seas"]]
    if (!is.null(t_ss)) {
      for (sm_nm in plan$smooth_names) {
        pred_fn <- draws$smooth_pred_fn[[sm_nm]]
        coef_row <- cf$smooth[[sm_nm]]
        if (is.null(pred_fn) || is.null(coef_row)) {
          next
        }
        X_s <- if (!is.null(smooth_X[[sm_nm]])) {
          smooth_X[[sm_nm]]
        } else {
          pred_fn(t_ss)
        } # N × K
        lp <- lp + as.numeric(X_s %*% coef_row)
      }
    }
  }

  lp + basehaz_val
}

#' Compute the basehaz contribution for one time step.
#' v1 supports two kinds:
#'   "const"             — flat (synthetic fixture)
#'   "brms_seasonal_in_lp" — already baked into FE via s(t_start_seas);
#'                           the kernel passes 0 here and t_start_seas via
#'                           the FE matrix (TODO: real brms needs the design
#'                           matrix; see kernel notes).
basehaz_eval <- function(draws, t_start_seas, N) {
  bh <- draws$basehaz
  switch(
    bh$kind %||% "const",
    const = rep(bh$const %||% 0, N),
    brms_seasonal_in_lp = rep(0, N),
    none = rep(0, N),
    rep(0, N)
  )
}

# ---------------------------------------------------------------------------
# Duration: months-since-onset (mso) baseline-hazard helpers
# ---------------------------------------------------------------------------

# Mirrors analysis/utils/dur_data.R's .DUR_MSO_BREAKS exactly. Duplicated
# rather than sourced -- same rationale as bin_ttlo() just above (this file
# is a pure, self-contained kernel library; data-prep modules are never
# source()-d from here). Overridden per-fit by draws_dur$basehaz$breaks,
# which carries the actual duration.mso_breaks the fit was built with.
.SIM_DUR_MSO_BREAKS <- c(1, 2, 3, 4, 6, 9, Inf)

#' Bin an integer months-since-onset state into the mso_band factor labels
#' used at fit time. Mirrors dur_data.R::bin_mso() exactly, including the
#' `right = FALSE, include.lowest = TRUE` cut() semantics and label format
#' ("1", "2", "3", "4-5", "6-8", "9+" for the default breaks).
#'
#' @param mso numeric/integer vector; months since onset (1 = onset month).
#'   NA / out-of-range values (e.g. mso = 0 for a not-yet-active location)
#'   map to NA.
#' @param breaks numeric; left edges plus a terminal Inf.
#' @return character vector of band labels, NA where mso is NA or falls
#'   outside [breaks[1], Inf).
.dur_bin_mso <- function(mso, breaks = .SIM_DUR_MSO_BREAKS) {
  labs <- character(length(breaks) - 1L)
  for (i in seq_len(length(breaks) - 1L)) {
    lo <- breaks[i]
    hi <- breaks[i + 1L]
    labs[i] <- if (is.infinite(hi)) {
      sprintf("%d+", as.integer(lo))
    } else if (hi - lo == 1) {
      sprintf("%d", as.integer(lo))
    } else {
      sprintf("%d-%d", as.integer(lo), as.integer(hi) - 1L)
    }
  }
  # .bincode() is what cut() computes its codes with; indexing the labels
  # directly skips the factor round trip (same result as.character(cut(...))).
  labs[.bincode(as.numeric(mso), breaks, right = FALSE, include.lowest = TRUE)]
}

#' Sanitize an mso_band factor level to the brms fixed-effect column suffix.
#' brms::rename() replaces "-" -> "M" and "+" -> "P" in factor-level labels
#' when building `b_*` column names; verified empirically against a fitted
#' duration model's `fixef()` rownames: levels "1","2","3","4-5","6-8","9+"
#' -> columns "mso_band2","mso_band3","mso_band4M5","mso_band6M8","mso_band9P"
#' (the reference level "1" is absorbed into the intercept, as usual).
.dur_mso_level_to_col <- function(level) {
  paste0(
    "mso_band",
    gsub("-", "M", gsub("+", "P", level, fixed = TRUE), fixed = TRUE)
  )
}

#' Compute the primary (months-since-onset) baseline-hazard contribution to
#' the duration linear predictor, for one posterior draw, one time step.
#'
#' Companion to basehaz_eval() (the occurrence track's single-clock
#' analogue), but doing real per-band FE lookup rather than delegating to a
#' pre-baked spline: the duration track's default basehaz_duration is
#' "factor" (unlike occurrence's spline-default "seasonal"), so there is no
#' existing basehaz_eval() code path that already does the right thing.
#'
#' @param draws_dur duration draws bundle (extract_duration_draws()).
#' @param s scalar posterior draw index.
#' @param mso numeric vector length N; current months-since-onset state.
#'   Callers subset to the active locations before calling (mso is
#'   meaningless for inactive locations); NA/out-of-range entries simply
#'   contribute 0 (matched to no FE column) rather than erroring, so passing
#'   the full N-length state vector is harmless.
#' @param N number of locations.
#' @return numeric vector length N.
basehaz_dur_eval <- function(draws_dur, s, mso, N) {
  bh <- draws_dur$basehaz
  switch(
    bh$kind %||% "none",
    factor = {
      out <- rep(0, N)
      lvl <- .dur_bin_mso(mso, bh$breaks %||% .SIM_DUR_MSO_BREAKS)
      ok <- !is.na(lvl)
      if (any(ok)) {
        cols <- .dur_mso_level_to_col(lvl[ok])
        idx <- match(cols, colnames(draws_dur$fe))
        got <- !is.na(idx)
        if (any(got)) out[which(ok)[got]] <- draws_dur$fe[s, idx[got]]
      }
      out
    },
    none = rep(0, N),
    rep(0, N)
  )
}

# ---------------------------------------------------------------------------
# TTLO helpers — categorical binning (step 4)
# ---------------------------------------------------------------------------

# Cutpoints mirror occ_data.R::make_time_to_last_outbreak_cat().
.TTLO_BREAKS <- c(0, 1, 2, 3, 4, 5, 6, Inf)

#' Bin continuous days-since-last-outbreak into the factor-level labels used
#' at fit time.  Mirrors `cut((pmax(0,days)+1)/365, .TTLO_BREAKS,
#' include.lowest=TRUE, right=FALSE)` from
#' occ_data.R::make_time_to_last_outbreak_cat() (occ track; size track: the
#' identical /365 binning in siz_data.R::read_size_data()).
#'
#' @param days numeric vector; days elapsed since last outbreak end.
#'   NA → NA (warm-up not available; treated as reference-level silence in LP).
#' @param has_none logical; TRUE when the fit carries the "no recorded prior
#'   onset" level (`time_to_last_outbreaknone`): days at or
#'   above base.R::.TTLO_NONE_DAYS (the sentinel seeded for a unit with no
#'   recorded onset, kept by advance_ttlo() until its first onset) bin to
#'   "none". FALSE (default): the sentinel falls in "[6,Inf]", which is how a
#'   fit without the level coded those rows.
#' @return character vector of factor-level labels: "[0,1)", "[1,2)", …, "[6,Inf]"
#'   (base.R::.TTLO_LEVELS, exactly the labels cut() derives from these breaks),
#'   plus "none" under has_none.
bin_ttlo <- function(days, has_none = FALSE) {
  years <- (pmax(0, days) + 1) / 365
  # .bincode(): cut()'s own coding step, without the factor round trip.
  out <- .TTLO_LEVELS[.bincode(
    years,
    .TTLO_BREAKS,
    right = FALSE,
    include.lowest = TRUE
  )]
  if (isTRUE(has_none)) {
    out[!is.na(days) & days >= .TTLO_NONE_DAYS] <- .TTLO_NONE_LEVEL
  }
  out
}

# The "no recorded prior onset" TTLO level: label "none", appended AFTER "[6,Inf]"
# ("[0,1)" stays the
# reference), brms column `time_to_last_outbreaknone`, sentinel
# base.R::.TTLO_NONE_DAYS days.
.TTLO_NONE_LEVEL <- "none"
.TTLO_NONE_COL <- "time_to_last_outbreaknone"

#' .draws_has_ttlo_none
#'
#' TRUE when a draws bundle carries the "none" TTLO level (FE column
#' `time_to_last_outbreaknone`). Given a list of bundles (the grouped draws),
#' they must agree: a set mixing groups with and without the level is refused
#' (one seeding rule per run). NULL / empty -> FALSE.
.draws_has_ttlo_none <- function(draws) {
  one <- function(d) {
    .TTLO_NONE_COL %in% colnames(d$fe %||% matrix(nrow = 0, ncol = 0))
  }
  if (is.null(draws)) {
    return(FALSE)
  }
  if (!is.null(draws$fe)) {
    return(one(draws))
  }
  hits <- vapply(Filter(Negate(is.null), draws), one, logical(1))
  if (!length(hits)) {
    return(FALSE)
  }
  if (length(unique(hits)) > 1L) {
    stop(
      ".draws_has_ttlo_none: the grouped draws mix fits with and without the TTLO level ",
      "'none' (",
      paste(sprintf("%s=%s", names(hits), hits), collapse = ", "),
      "); refit the groups under one occurrence.ttlo_no_prior_level.",
      call. = FALSE
    )
  }
  hits[[1L]]
}

#' Map factor-level label back to a representative days value (lower bound).
#'
#' Used to initialise `u_days` from the panel's `time_to_last_outbreak_cat`
#' column at warm-up.  Unknown or NA labels default to the top bin (6 years),
#' matching the fit's `replace_na(list(time_to_last_outbreak = "[6,Inf]"))`.
ttlo_init_from_panel <- function(ttlo_cat_vec) {
  # lower bound in years of each base.R::.TTLO_LEVELS bin
  lo_year <- stats::setNames(c(0, 1, 2, 3, 4, 5, 6), .TTLO_LEVELS)
  chr <- as.character(ttlo_cat_vec)
  days <- lo_year[chr] * 365.25
  names(days) <- NULL
  days[is.na(days)] <- 6 * 365.25 # unknown / NA → top bin
  # "none": the no-prior-onset sentinel (a handoff state that carried it).
  days[!is.na(chr) & chr == .TTLO_NONE_LEVEL] <- .TTLO_NONE_DAYS
  days
}

#' Convert a brms TTLO FE column name back to the R factor-level label.
#'
#' Brms sanitises factor-level labels: "[1,2)" → "12", "[6,Inf]" → "6Inf".
#' This reversal is needed to map `draws$fe` column names back to the bin
#' labels returned by `bin_ttlo()`.
ttlo_col_to_level <- function(col_name) {
  inner <- sub("^time_to_last_outbreak", "", col_name)
  if (identical(inner, .TTLO_NONE_LEVEL)) {
    return(.TTLO_NONE_LEVEL)
  }
  if (grepl("Inf", inner)) {
    lo <- sub("Inf$", "", inner)
    paste0("[", lo, ",Inf]")
  } else {
    n <- nchar(inner)
    lo <- substr(inner, 1, n - 1)
    hi <- substr(inner, n, n)
    paste0("[", lo, ",", hi, ")")
  }
}

# ---------------------------------------------------------------------------
# State propagation
# ---------------------------------------------------------------------------

#' Apply one step of the nested-OR contagion update.
#'
#' Inputs Y_i = simulated occurrence at month m (integer 0/1 length N).
#' Returns the next month's lag vectors:
#'   o_next     (length N) = ADM1 lag for each location
#'   o_ct_next  (length N) = country leave-one-ADM1-out lag
#'   o_tu_next  (length N) = TU leave-one-country-out lag; **NA** for locations
#'     whose country is absent from the HMM cluster file (`frame_tu$tu == NA`),
#'     matching fit-time behaviour where `ongoing_outbreak_tu_lag1 = NA`.
#'
#' frame_tu provides adm1_unit, country, tu per location.
#'
#' `maps` (optional): build_lag_maps(frame_tu), the integer codes of the same
#' aggregation, built once per frame. With maps and a 0/1 Y without NA the
#' update runs on integer codes (tabulate()); otherwise -- no maps, or any NA /
#' non-0/1 value -- the keyed tapply() path below. Both return identical lists.
update_lags <- function(Y, frame_tu, maps = NULL) {
  if (!is.null(maps) && .is01(Y)) {
    return(.update_lags_coded(Y, maps))
  }
  N <- length(Y)
  # ADM1 active state: any location in ADM1 j active at m
  A_adm1 <- tapply(Y, frame_tu$adm1_unit, max, default = 0L)
  # Country active state: any ADM1 of country c active
  adm1_to_country <- frame_tu %>%
    dplyr::distinct(adm1_unit, country) %>%
    tibble::deframe()
  # Order A_adm1 vector by names; build per-country sum + max
  adm1_names <- names(A_adm1)
  ctry_of_adm1 <- adm1_to_country[adm1_names]

  # ADM1-by-country aggregates
  A_ct_max <- tapply(as.integer(A_adm1), ctry_of_adm1, max, default = 0L)
  A_ct_sum <- tapply(as.integer(A_adm1), ctry_of_adm1, sum, default = 0L)

  # TU active state per country: any other country of TU active
  ctry_to_tu <- frame_tu %>%
    dplyr::distinct(country, tu) %>%
    tibble::deframe()
  ctry_names <- names(A_ct_max)
  tu_of_ctry <- ctry_to_tu[ctry_names]

  # tapply silently drops NA keys → NA-TU countries excluded from TU aggregates
  A_tu_sum <- tapply(as.integer(A_ct_max), tu_of_ctry, sum, default = 0L)

  # Broadcast back to per-location
  o_next <- as.integer(A_adm1[frame_tu$adm1_unit])
  # Country leave-one-ADM1-out: (sum over ADM1 of country c — self) >= 1
  ct_for_loc <- frame_tu$country
  ct_sum_loc <- as.integer(A_ct_sum[ct_for_loc])
  self_adm1 <- as.integer(A_adm1[frame_tu$adm1_unit])
  o_ct_next <- as.integer((ct_sum_loc - self_adm1) >= 1L)

  # TU leave-one-country-out: (sum over countries of TU — self_country) >= 1
  # Locations with NA TU → NA (never coerced to 0 — they match fit-time NA rows).
  tu_for_loc <- frame_tu$tu
  tu_sum_loc <- as.integer(ifelse(
    is.na(tu_for_loc),
    NA_integer_,
    A_tu_sum[tu_for_loc]
  ))
  self_country <- as.integer(A_ct_max[ct_for_loc])
  o_tu_next <- as.integer((tu_sum_loc - self_country) >= 1L)
  # NA preserved: do NOT zero-fill NA-TU rows.

  # A_adm1 (named by ADM1 unit, the lag-1 ADM1-active 0/1 grid) is returned so
  # the neighbour-contagion caller can form `K %*% Aa1` without recomputing the
  # same tapply. Existing consumers ignore the extra field.
  list(
    o = o_next,
    o_ct = o_ct_next,
    o_tu = o_tu_next,
    A_adm1 = as.integer(A_adm1) |> stats::setNames(names(A_adm1))
  )
}

# TRUE when v is all 0/1 with no NA (the integer-coded fast paths' domain).
.is01 <- function(v) !anyNA(v) && all(v == 0 | v == 1)

#' Integer codes of update_lags()' ADM1 -> country -> TU aggregation, built once
#' per frame (kernel speed-up). Level orders are the ones tapply() uses there
#' (factor() of the keys), and every map takes the FIRST row of a key, as the
#' dplyr::distinct() + deframe() lookups of the keyed path do.
build_lag_maps <- function(frame_tu) {
  adm1_f <- factor(frame_tu$adm1_unit)
  adm1_lv <- levels(adm1_f)
  first_of_adm1 <- match(adm1_lv, frame_tu$adm1_unit)
  ctry_of_adm1 <- frame_tu$country[first_of_adm1]
  ctry_lv <- levels(factor(ctry_of_adm1))
  ctry_tu <- frame_tu$tu[match(ctry_lv, frame_tu$country)]
  tu_lv <- levels(factor(ctry_tu))
  list(
    n_adm1 = length(adm1_lv),
    adm1_lv = adm1_lv,
    adm1_i = as.integer(adm1_f), # row  -> adm1 code
    adm1_ctry = match(ctry_of_adm1, ctry_lv), # adm1 -> country code
    n_ctry = length(ctry_lv),
    ctry_i = match(frame_tu$country, ctry_lv), # row  -> country code
    ctry_tu = match(ctry_tu, tu_lv), # country -> tu code (NA)
    n_tu = length(tu_lv),
    row_tu = match(frame_tu$tu, tu_lv)
  ) # row  -> tu code (NA)
}

# Per-ADM1 max of a 0/1 vector, in build_lag_maps() level order (= tapply(v, adm1, max)).
.adm1_any <- function(v, maps) {
  A <- integer(maps$n_adm1)
  A[maps$adm1_i[v == 1]] <- 1L
  A
}

# update_lags() on integer codes; caller guarantees .is01(Y).
.update_lags_coded <- function(Y, maps) {
  A_adm1 <- .adm1_any(Y, maps)
  A_ct_sum <- tabulate(maps$adm1_ctry[A_adm1 == 1L], maps$n_ctry)
  A_ct_max <- as.integer(A_ct_sum > 0L)
  on_tu <- A_ct_max == 1L & !is.na(maps$ctry_tu)
  A_tu_sum <- tabulate(maps$ctry_tu[on_tu], maps$n_tu)
  o_next <- A_adm1[maps$adm1_i]
  o_ct_next <- as.integer((A_ct_sum[maps$ctry_i] - o_next) >= 1L)
  # NA-TU rows stay NA (never coerced to 0 -- they match fit-time NA rows).
  o_tu_next <- as.integer((A_tu_sum[maps$row_tu] - A_ct_max[maps$ctry_i]) >= 1L)
  list(
    o = o_next,
    o_ct = o_ct_next,
    o_tu = o_tu_next,
    A_adm1 = stats::setNames(A_adm1, maps$adm1_lv)
  )
}

#' TTLO advancement. Step 4: tracks elapsed days; advances by `step_days`
#' (~30.44 days per monthly step) and resets to 0 on outbreak.
#'
#' @param u_days numeric vector; days elapsed since last outbreak end.
#' @param Y integer vector 0/1; simulated occurrence at this step.
#' @param step_days numeric; days per simulation step (default 365.25/12).
#' @return numeric vector; updated days counter.
advance_ttlo <- function(u_days, Y, step_days = 365.25 / 12) {
  hit <- as.integer(Y) == 1L
  out <- as.numeric(u_days + step_days) # as ifelse(): no attributes kept
  out[which(hit)] <- 0.0
  out[is.na(hit)] <- NA_real_
  out
}

#' Per-episode surveillance-detection flag -- the `covariate_feed: detected`
#' state update of the v1 kernel (simulate_one_tu()); the Csnippet emits the
#' same rule in the same unit order (spatpomp_core_c.R, detect_block), so R and
#' C stay bit-exact.
#'
#' The occurrence fit learned its spillover covariates (o, o_ct, o_tu, nbr,
#' time_to_last_outbreak) from DETECTED episodes. Feeding the kernel the latent
#' field instead (5-20x more prevalent) evaluates those coefficients far
#' outside what they were estimated on. With this flag the kernel feeds what
#' the surveillance record would have seen, using the same delta_fn/delta_fp as
#' the detect measurement:
#'   onset this month           det ~ Bernoulli(1 - delta_fn)  (kept for the episode)
#'   active, no onset           det unchanged
#'   inactive                   det ~ Bernoulli(delta_fp)      (monthly false alarm)
#' Draws happen in ascending unit order, one per unit that needs one.
#'
#' @param onset      0/1 onset indicator this month.
#' @param active_now 0/1 activity AFTER this month's onset promotion.
#' @param det_prev   0/1 flag carried from the previous month.
#' @return integer 0/1 vector: the detected-activity field for this month.
update_detection <- function(onset, active_now, det_prev, delta_fn, delta_fp) {
  det <- as.integer(det_prev)
  q <- rep(NA_real_, length(det))
  q[onset == 1L] <- 1 - delta_fn
  q[onset != 1L & active_now == 0L] <- delta_fp
  i <- which(!is.na(q))
  if (length(i)) {
    det[i] <- as.integer(stats::rbinom(length(i), 1L, q[i]))
  }
  det
}

# ---------------------------------------------------------------------------
# Phi (NegBin dispersion) helper for the monthly count kernel
# ---------------------------------------------------------------------------

#' Compute a length-N phi (overdispersion) vector for one posterior draw.
#'
#' When phi_re == "constant", returns rep(phi_scalar, N) for backward compat.
#' When phi_re is a nested level, accumulates the FE intercept + RE contributions
#' on the log scale (brms log-link for phi) and returns exp() of the sum.
#' Uses the same re_level_key() / match() pattern as evaluate_lp().
#'
#' @param draws_nb    NB count draws bundle (output of extract_monthly_draws()).
#' @param s           scalar posterior draw index.
#' @param frame_tu    location frame for the TU (N rows).
#' @return  numeric vector length N; all values are finite and > 0.
compute_phi_vec <- function(draws_nb, s, frame_tu, drop_region = FALSE) {
  phi_re <- draws_nb$phi_re %||% "constant"
  N <- nrow(frame_tu)
  phi <- draws_nb$phi

  if (phi_re == "constant") {
    phi_s <- if (is.list(phi)) phi$fe[s] else phi[s]
    return(rep(phi_s, N))
  }

  # Non-constant: accumulate log(phi) across nested levels, then exp().
  # phi$fe is b_shape_Intercept (log scale).
  log_phi <- rep(phi$fe[s], N)

  add_phi_re_level <- function(log_phi, mat, level_keys) {
    if (is.null(mat)) {
      return(log_phi)
    }
    idx <- match(level_keys, colnames(mat))
    have <- !is.na(idx)
    if (any(have)) {
      log_phi[have] <- log_phi[have] + mat[s, idx[have]]
    }
    log_phi
  }

  # Under grouping the region level is dropped; the phi shape RE chain re-roots
  # at country (see the NB count formula builders in utils/config.R), so the lookup keys must drop the
  # region prefix too.
  if (phi_re %in% c("region", "region_country", "region_country_adm1")) {
    log_phi <- add_phi_re_level(
      log_phi,
      phi$region,
      re_level_key(frame_tu, "region", drop_region)
    )
  }
  if (phi_re %in% c("region_country", "region_country_adm1")) {
    log_phi <- add_phi_re_level(
      log_phi,
      phi$country,
      re_level_key(frame_tu, "country", drop_region)
    )
  }
  if (phi_re == "region_country_adm1") {
    log_phi <- add_phi_re_level(
      log_phi,
      phi$adm1,
      re_level_key(frame_tu, "adm1", drop_region)
    )
  }

  exp(log_phi)
}

# ---------------------------------------------------------------------------
# Monthly (`mon`) case-count mean and draw -- shared by every kernel
# ---------------------------------------------------------------------------
# One recipe, two call sites: simulate_one_tu() (v1 and the spatPomp forecast
# phase) and extract_smoothing_one_tu() (smoothing post-hoc
# draw, smoothing_core.R). The Csnippet kernel (spatpomp_core_c.R mon_mu_block)
# keeps its own C expression and shares only the ceiling, monthly_mu_cap().

#' Inject the latent months-since-onset cross-basis columns into an xrow list.
#'
#' `tab` is a `build_mso_basis()` result: `$B` is a `(mso_max + 1) x K` matrix
#' whose row `j + 1` holds the basis evaluated at `mso = j`, and `$mso_max` is
#' its last row's mso. Every other cross-basis variable is a panel column and
#' arrives through the covariate table; `mso` is latent -- created inside the
#' month loop and different per particle -- so it is looked up here instead.
#'
#' Clamping to `[0, mso_max]` is MANDATORY, not defensive: the shipped basis is
#' `bs(degree = 3)`, which extrapolates as an unbounded cubic past its boundary
#' knots. The repo documents this exact failure for CFSv2 ENSO
#' (config_defaults.yml: 0.8-1.0 units past the right knot inflated mu 26x) and
#' already clamps for it in smoothing_core.R. A long-running episode's mu
#' therefore flattens at the last trained value instead of diverging.
#'
#' @param xm      named list of per-unit covariate vectors (an `xrow`).
#' @param tab     `build_mso_basis()` output, or NULL (xm returned unchanged).
#' @param mso_idx integer vector of latent mso values, one per row of `xm`.
#' @return `xm` with `$B`'s columns added under their fitted names.
.inject_mso <- function(xm, tab, mso_idx) {
  if (is.null(tab)) {
    return(xm)
  }
  rows <- pmin(pmax(mso_idx, 0L), tab$mso_max) + 1L
  cn <- colnames(tab$B)
  for (k in seq_len(ncol(tab$B))) {
    xm[[cn[k]]] <- tab$B[rows, k]
  }
  xm
}

#' Log monthly case mean for one block of units.
#'
#' evaluate_lp() on `xm` (with the latent mso basis injected) plus the three
#' bundle-driven offsets, added in this order: log(pop) (`offset_pop`),
#' log_offset + the onset-month correction (`offset_logoffset`), the base-rate
#' anchor (`offset_baserate`). For a mon bundle these are offset(log(pop)) +
#' offset(log_offset) [+ offset(log(base_rate_mon))]; see
#' build_monthly_formula_parts().
#'
#' @param draws,s      mon draws bundle and the posterior draw index.
#' @param xm           named list of per-unit covariate vectors, cross-basis
#'                     columns already in place.
#' @param o,o_ct,o_tu,u_days  spillover lags and TTLO, one per row of `frame_g`.
#' @param frame_g      location-frame rows being scored (supplies RE keys, `pop`
#'                     and the anchor column).
#' @param log_offset   log exposure days (vector or scalar).
#' @param mso          latent months-since-onset, one per row (1 = onset month).
#' @param mso_tab      build_mso_basis() table, or NULL (no mso basis).
#' @param what         label for the anchor's error messages.
#' @param pop_floor    NULL (log(pop)) or a floor, log(pmax(pop, pop_floor)).
#' @param pop          NULL (the frame's static `pop`) or this month's
#'                     population per row (frame_pop_at(); pop_year mode): it
#'                     replaces `frame_g$pop` for both the
#'                     log(pop) offset and evaluate_lp()'s `log_pop_c`.
#' @param ceased       NULL (simulator.cessation_exposure: full) or a logical
#'                     per row, TRUE where the unit ceases this month: adds
#'                     .cessation_log_offset_adj() on the log_offset path.
#' @return numeric vector, one log mean per row of `frame_g`.
monthly_log_mu <- function(
  draws,
  s,
  xm,
  o,
  o_ct,
  o_tu,
  u_days,
  frame_g,
  log_offset,
  mso,
  mso_tab = NULL,
  what = "monthly kernel",
  drop_region = FALSE,
  drop_ct_tu = FALSE,
  pop_floor = NULL,
  pop = NULL,
  plan = NULL,
  smooth_X = NULL,
  ceased = NULL
) {
  if (!is.null(pop)) {
    frame_g$pop <- pop
  }
  xm <- .inject_mso(xm, mso_tab, mso)
  lm <- evaluate_lp(
    draws,
    s,
    xm,
    o,
    o_ct,
    o_tu,
    u_days,
    frame_g,
    basehaz_val = 0,
    drop_region = drop_region,
    drop_ct_tu = drop_ct_tu,
    plan = plan,
    smooth_X = smooth_X
  )
  if (isTRUE(draws$offset_pop %||% TRUE)) {
    lm <- lm +
      log(if (is.null(pop_floor)) frame_g$pop else pmax(frame_g$pop, pop_floor))
  }
  # Onset-month exposure correction: `log_offset` is log(days_in_month(t)),
  # exact for interior months but ~2x too long for the onset month the mon
  # panel carries. See .ONSET_EXPOSURE_FRAC.
  if (isTRUE(draws$offset_logoffset %||% FALSE)) {
    lm <- lm + log_offset + .onset_log_offset_adj(mso)
  }
  # Cessation-month correction (simulator.cessation_exposure: half; NULL under
  # full, so the default adds nothing). See .CESSATION_EXPOSURE_FRAC.
  if (!is.null(ceased) && isTRUE(draws$offset_logoffset %||% FALSE)) {
    lm <- lm + .cessation_log_offset_adj(ceased)
  }
  # offset(log(base_rate_mon)): its denominator is pop x exposure days, the
  # two offsets above. 0 for an un-anchored fit.
  lm + .anchor_log_offset(draws, frame_g, NULL, what)
}

# Monthly incidence cap: monthly cases per unit <= cap x pop. The monthly NB
# has a heavy tail (shape ~0.4 in the deployed fits) and nothing else bounds the
# mean or the draw (a smoothing run drew 149,266 cases in one month for a unit
# of 182,755 people). Applied to the mean (min(mu, cap x pop)) and to the
# drawn count (min(draw, floor(cap x pop))), in every kernel.
.INCIDENCE_CAP_DEFAULT <- 0.15

#' Read and validate `simulator.incidence_cap`: the one place the config value
#' is interpreted. Missing -> .INCIDENCE_CAP_DEFAULT; `.inf` (Inf) disables.
#' @param sim the `simulator:` config block (or NULL).
#' @return a single number > 0, possibly Inf.
incidence_cap_of <- function(sim) {
  v <- (sim %||% list())$incidence_cap %||% .INCIDENCE_CAP_DEFAULT
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0) {
    stop(
      "simulator.incidence_cap must be a single number > 0 (.inf disables); got: ",
      paste(format(v), collapse = ", ")
    )
  }
  as.numeric(v)
}

# Onset eligibility: which inactive units may draw an onset.
#   inactive  every unit not in a running episode (default).
#   season    also no onset yet in the unit's current cholera season (the
#             occurrence fit models the first onset of each gid season).
# Gate: n = floor(u_days / 30.4375 + 0.5), u_days read BEFORE this month's
# update; seas_k = (month - min_seas) mod 12. Eligible iff n >= seas_k; the
# no-prior-onset sentinel (.TTLO_NONE_DAYS) and NA always are. floor(x + 0.5),
# not round(): it matches the C kernel, so both kernels gate the same units
# and consume the RNG stream identically.
.ONSET_ELIGIBILITY_DEFAULT <- "inactive"
.ONSET_ELIGIBILITY_VALUES <- c("inactive", "season")

#' Read and validate `simulator.onset_eligibility`: the one place the config
#' value is interpreted. Missing -> .ONSET_ELIGIBILITY_DEFAULT.
#' @param sim the `simulator:` config block (or NULL).
#' @return "inactive" or "season".
onset_eligibility_of <- function(sim) {
  v <- (sim %||% list())[["onset_eligibility"]] %||% .ONSET_ELIGIBILITY_DEFAULT
  check_onset_eligibility(v)
}

#' Validate an onset_eligibility value (kernel arguments and config alike).
check_onset_eligibility <- function(v, where = "simulator.onset_eligibility") {
  if (
    !is.character(v) ||
      length(v) != 1L ||
      is.na(v) ||
      !v %in% .ONSET_ELIGIBILITY_VALUES
  ) {
    stop(
      where,
      " must be one of: ",
      paste(.ONSET_ELIGIBILITY_VALUES, collapse = ", "),
      "; got: ",
      paste(format(v), collapse = ", "),
      call. = FALSE
    )
  }
  v
}

#' Refuse the season gate on a frame whose `min_seas` was defaulted for some
#' gid it still holds: the gate's season boundary would be a guess there.
#' `defaulted` = the gids build_location_frame() recorded in attr
#' "min_seas_defaulted"; read it right after the build (later frame steps need
#' not keep the attribute) and call this on the FINAL frame, after scope and
#' overlap drops, so a gid dropped anyway does not stop the run. A no-op under
#' "inactive".
check_onset_eligibility_frame <- function(
  frame,
  onset_eligibility,
  where = "build_location_frame",
  defaulted = attr(frame, "min_seas_defaulted")
) {
  if (!identical(check_onset_eligibility(onset_eligibility), "season")) {
    return(invisible(frame))
  }
  if (!"min_seas" %in% names(frame) || anyNA(frame$min_seas)) {
    stop(
      sprintf(
        "[%s] simulator.onset_eligibility: season needs min_seas for every gid.",
        where
      ),
      call. = FALSE
    )
  }
  dflt <- intersect(as.character(defaulted), as.character(frame$gid))
  if (length(dflt)) {
    stop(
      sprintf(
        paste0(
          "[%s] simulator.onset_eligibility: season needs each gid's own season ",
          "start, but min_seas was defaulted to 1 for %d gid(s) (e.g. %s). Give the ",
          "panel a min_seas (or t_start_seas) for them, or use onset_eligibility: inactive."
        ),
        where,
        length(dflt),
        paste(utils::head(dflt, 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(frame)
}

#' Month-of-season index per (unit, month): (month(m) - min_seas) mod 12,
#' 0 on the season's first month. Same rotation as t_start_seas.
#' @return integer matrix, length(min_seas) x length(months).
season_month_index <- function(min_seas, months) {
  mon <- as.integer(lubridate::month(months))
  outer(as.integer(min_seas), mon, function(ms, m) (m - ms) %% 12L)
}

#' The season gate on top of `active == 0`: TRUE where no onset that reset
#' u_days fell in the current season (see the comment on
#' .ONSET_ELIGIBILITY_DEFAULT for the derivation).
#' @param u_days numeric, pre-update days since the last onset.
#' @param seas_k integer month-of-season index of the month being simulated.
season_onset_open <- function(u_days, seas_k) {
  is.na(u_days) | floor(u_days / 30.4375 + 0.5) >= seas_k
}

# Cessation-month exposure: what the monthly case mean of a
# unit whose episode ends this month is exposed to.
#   full  a full month (the default): the cessation is drawn
#         after the monthly mean, so the kernel cannot know the month is the
#         last.
#   half  the cessation is drawn before the monthly mean and a ceasing month
#         gets .CESSATION_EXPOSURE_FRAC of a month, as the fit's panel does
#         (.expand_active_months(): the last month's exposure ends at TR). See
#         the .ONSET_EXPOSURE_FRAC comment for the order and lag reads.
.CESSATION_EXPOSURE_DEFAULT <- "full"
.CESSATION_EXPOSURE_VALUES <- c("full", "half")

#' Read and validate `simulator.cessation_exposure`: the one place the config
#' value is interpreted. Missing -> .CESSATION_EXPOSURE_DEFAULT.
#' @param sim the `simulator:` config block (or NULL).
#' @return "full" or "half".
cessation_exposure_of <- function(sim) {
  v <- (sim %||% list())[["cessation_exposure"]] %||%
    .CESSATION_EXPOSURE_DEFAULT
  check_cessation_exposure(v)
}

#' Validate a cessation_exposure value (kernel arguments and config alike).
check_cessation_exposure <- function(
  v,
  where = "simulator.cessation_exposure"
) {
  if (
    !is.character(v) ||
      length(v) != 1L ||
      is.na(v) ||
      !v %in% .CESSATION_EXPOSURE_VALUES
  ) {
    stop(
      where,
      " must be one of: ",
      paste(.CESSATION_EXPOSURE_VALUES, collapse = ", "),
      "; got: ",
      paste(format(v), collapse = ", "),
      call. = FALSE
    )
  }
  v
}

#' Per-unit ceiling on the monthly case mean: min(incidence_cap, mu_pop_frac)
#' x pop, or Inf for every unit when neither is set (or both are Inf).
#' `mu_pop_frac` is an optional extra fraction (no config key sets it).
#' @return numeric vector the length of `pop`.
monthly_mu_cap <- function(pop, incidence_cap = NULL, mu_pop_frac = NULL) {
  frac <- min(c(incidence_cap, mu_pop_frac, Inf))
  if (!is.finite(frac)) {
    return(rep(Inf, length(pop)))
  }
  pop * frac
}

# ---------------------------------------------------------------------------
# Population year of the location frame
# ---------------------------------------------------------------------------
# `simulator.location_frame.pop_year` is the year at which every gid's
# population Q_g is resolved for the simulator and the smoothing. An integer
# (default 2020: one population for every simulated
# month) or a mode: `season_year` resolves each month's population at the
# month's cholera season year (year(t - (min_seas - 1) months), the fits' own
# `year_seas` convention, 03_prepare_siz.R) and `calendar_year` at year(t).
# In a mode build_location_frame() adds the per-year columns `pop_y<YYYY>` and
# the per-gid shift `pop_shift_m`; the kernels read the month's population
# through frame_pop_at() / frame_pop_matrix(), which return the static `pop` for a
# frame without those columns, so an integer pop_year changes nothing.

#' Default `simulator.location_frame.pop_year` (= get_default_options(); pinned
#' in tests/test_sim_pop_year.R). The one literal: every reader goes through
#' sim_pop_year_spec().
.SIM_POP_YEAR_DEFAULT <- 2020L
.SIM_POP_YEAR_MODES <- c("season_year", "calendar_year")

#' sim_pop_year_spec
#'
#' The one reader of `simulator.location_frame.pop_year`. Also the guard on the
#' dead look-alike `simulator.population.base_year` (defaulted, never read):
#' any value other than 2020 is refused with a pointer to the real key.
#' @param sim the `simulator:` config block (or NULL = defaults).
#' @return an integer year in 1950-2100, or one of `.SIM_POP_YEAR_MODES`.
sim_pop_year_spec <- function(sim) {
  sim <- sim %||% list()
  by <- (sim[["population"]] %||% list())[["base_year"]]
  if (
    !is.null(by) && !(is.numeric(by) && length(by) == 1L && isTRUE(by == 2020))
  ) {
    stop(
      sprintf(
        paste0(
          "simulator.population.base_year = %s: this key is read by no code ",
          "(only its default 2020 is accepted). The simulator population ",
          "year is simulator.location_frame.pop_year (an integer year, ",
          "`season_year` or `calendar_year`)."
        ),
        paste(format(by), collapse = ",")
      ),
      call. = FALSE
    )
  }
  v <- (sim[["location_frame"]] %||% list())[["pop_year"]]
  if (is.null(v)) {
    return(.SIM_POP_YEAR_DEFAULT)
  }
  if (is.character(v) && length(v) == 1L && v %in% .SIM_POP_YEAR_MODES) {
    return(v)
  }
  if (
    is.numeric(v) &&
      length(v) == 1L &&
      !is.na(v) &&
      v == round(v) &&
      v >= 1950 &&
      v <= 2100
  ) {
    return(as.integer(v))
  }
  stop(
    sprintf(
      paste0(
        "simulator.location_frame.pop_year must be an integer year in ",
        "1950-2100 or one of %s; got: %s"
      ),
      paste(.SIM_POP_YEAR_MODES, collapse = " / "),
      paste(format(v), collapse = ",")
    ),
    call. = FALSE
  )
}

#' TRUE when a pop_year spec is a mode (a per-month population lookup).
sim_pop_year_is_mode <- function(spec) is.character(spec)

#' sim_frame_pop_args
#'
#' The population arguments of build_location_frame() for one simulator
#' config, in one place (replaces the block duplicated in 07a, 07b_v1, 07c,
#' simulations.R x2 and smoothing_core.R x2). A configured
#' `location_frame.gid_pop_path` / `wpp_path` is used as given (a missing file
#' then errors in build_location_frame(): a config key is an instruction);
#' otherwise the default artefact path, NULL when absent
#' (.pop_artefact_or_null()).
#' @param sim       the `simulator:` config block.
#' @param opt_paths an opt for make_pop_run_file() (reads `outdir` / `population`),
#'   used only when the default path is not given.
#' @param gid_pop_path,wpp_path optional default paths (the smoothing's
#'   S$GID_POP_PATH / S$WPP_PATH) used in place of make_pop_run_file(opt_paths, .).
#' @return list(gid_pop_path, wpp_path, pop_year).
sim_frame_pop_args <- function(
  sim,
  opt_paths = NULL,
  gid_pop_path = NULL,
  wpp_path = NULL
) {
  lf <- (sim %||% list())[["location_frame"]] %||% list()
  list(
    gid_pop_path = lf[["gid_pop_path"]] %||%
      .pop_artefact_or_null(
        gid_pop_path %||% make_pop_run_file(opt_paths, "pop_gid")
      ),
    wpp_path = lf[["wpp_path"]] %||%
      .pop_artefact_or_null(
        wpp_path %||% make_pop_run_file(opt_paths, "pop_wpp")
      ),
    pop_year = sim_pop_year_spec(sim)
  )
}

#' season_year_of
#'
#' Cholera season year of calendar month(s) `month` for a gid with season start
#' `min_seas`: year(month - (min_seas - 1) months). Integer arithmetic on the
#' month index, so it is vectorised over both arguments (lubridate::period()
#' with a vector shift is not; see build_covar_for_tu()). Equal to
#' occ_data.R's `seas_of()` (tests/test_sim_pop_year.R).
season_year_of <- function(month, min_seas) {
  month <- as.Date(month)
  idx <- as.integer(format(month, "%Y")) *
    12L +
    as.integer(format(month, "%m")) -
    1L
  as.integer((idx - (as.integer(min_seas) - 1L)) %/% 12L)
}

.POP_YEAR_COL_RE <- "^pop_y[0-9]{4}$"

#' TRUE when the frame carries the per-year population columns (pop_year mode).
frame_pop_has_years <- function(frame) {
  any(grepl(.POP_YEAR_COL_RE, names(frame)))
}

#' frame_pop_at
#'
#' Population of each frame row at calendar month `month` (one Date): the
#' static `pop` for a frame without `pop_y<YYYY>` columns (an integer
#' pop_year; the default), else `pop_y<y>` with y = year(month - pop_shift_m
#' months) per row. A year outside the frame's columns is an error.
frame_pop_at <- function(frame, month) {
  if (!frame_pop_has_years(frame)) {
    return(frame$pop)
  }
  if (length(month) != 1L) {
    stop("frame_pop_at: `month` must be one Date")
  }
  y <- season_year_of(
    rep(as.Date(month), nrow(frame)),
    as.integer(frame$pop_shift_m) + 1L
  )
  col <- sprintf("pop_y%04d", y)
  miss <- setdiff(unique(col), names(frame))
  if (length(miss)) {
    stop(
      sprintf(
        paste0(
          "frame_pop_at: month %s needs population year(s) %s, which the ",
          "location frame does not carry (built for another window; ",
          "pass build_location_frame(pop_window =) covering every month read)."
        ),
        format(as.Date(month)),
        paste(sort(sub("pop_y", "", miss)), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  m <- as.matrix(as.data.frame(frame)[,
    sort(grep(.POP_YEAR_COL_RE, names(frame), value = TRUE)),
    drop = FALSE
  ])
  m[cbind(seq_len(nrow(frame)), match(col, colnames(m)))]
}

#' frame_pop_matrix
#'
#' U x M population matrix of a frame over `months` (frame_pop_at() per month),
#' or NULL for a frame without `pop_y<YYYY>` columns: the kernels keep their
#' static path (and their exact arithmetic) when this is NULL.
frame_pop_matrix <- function(frame, months) {
  if (!frame_pop_has_years(frame)) {
    return(NULL)
  }
  out <- vapply(
    seq_along(months),
    function(m) as.numeric(frame_pop_at(frame, months[m])),
    numeric(nrow(frame))
  )
  matrix(out, nrow = nrow(frame), ncol = length(months))
}

#' sim_pop_window
#'
#' The months whose population a forward run reads: .sim_panel_window()
#' (horizon, ENSO-sampler start, simulated warm-up) widened to the filter's
#' warm-up / forecast windows. Passed as build_location_frame(pop_window =);
#' only read under a pop_year mode.
sim_pop_window <- function(sim, enso_sim_months_start = NULL) {
  w <- .sim_panel_window(sim, enso_sim_months_start)
  flt <- (sim %||% list())[["filter"]] %||% list()
  .dt <- function(x) if (is.null(x)) as.Date(NA) else as.Date(x)
  ww <- flt[["warmup_window"]] %||% list()
  fw <- flt[["forecast_window"]] %||% list()
  st <- c(w[1], .dt(ww[["start"]]), .dt(fw[["start"]]))
  en <- c(w[2], .dt(fw[["end"]]))
  c(
    if (all(is.na(st))) as.Date(NA) else min(st, na.rm = TRUE),
    if (all(is.na(en))) as.Date(NA) else max(en, na.rm = TRUE)
  )
}

#' Negative-binomial monthly case draw, truncated at floor(cap).
#'
#' Units with a non-finite mean or a non-positive / non-finite dispersion are
#' not drawn (NA). With cap = Inf the RNG consumption and values are those of
#' a plain rnbinom() over the drawable units.
#' @param mu,phi per-unit mean and NB size.
#' @param cap    per-unit ceiling (monthly_mu_cap()), recycled; Inf = none.
#' @return numeric vector the length of `mu`, NA where not drawn.
draw_monthly_cases <- function(mu, phi, cap = Inf) {
  out <- rep(NA_real_, length(mu))
  ok <- is.finite(phi) & phi > 0 & is.finite(mu)
  if (any(ok)) {
    out[ok] <- stats::rnbinom(sum(ok), mu = mu[ok], size = phi[ok])
    cap <- rep_len(cap, length(mu))
    if (any(is.finite(cap[ok]))) out[ok] <- pmin(out[ok], floor(cap[ok]))
  }
  out
}

# ---------------------------------------------------------------------------
# Per-TU AR kernel
# ---------------------------------------------------------------------------

# Design choices of the forward kernel:
#   * pure-R kernel — no Rcpp
#   * sequential within a TU (nested contagion rules out month-level parallel)
#   * paired posterior strategy only (one draw index per replicate, shared by
#     every track)
#   * no size track: case counts come from the optional monthly (`mon`)
#     track's `cases_month`; there is no episode-level size draw
#' Simulate one (TU, replicate-chunk) task. Writes:
#'   <output_part>/event_table.<ext>     — long table of simulated outbreaks
#'   <output_part>/monthly_summary.<ext> — per (gid, month, replicate) Y
#'                                          (+ active, mu_mon, cases_month,
#'                                          detected when those tracks are on)
#'   <output_part>/_SUCCESS              — sentinel
#'
#' @param tu transmission-unit label.
#' @param frame_tu location frame restricted to the TU.
#' @param panel_tu climate panel restricted to the TU + horizon.
#' @param draws_occ occurrence draws bundle.
#' @param draws_dur Optional duration (cessation-hazard) draws bundle
#'   (extract_duration_draws()). NULL: the onset-only path, an unconditional
#'   per-step onset draw with no `active`/`mso` state and no cessation.
#'   Non-NULL: onset only when not active, `active`/`mso` tracked,
#'   `update_lags()` fed `active`, a cessation draw each active month. The two
#'   paths are a hard branch on purpose: masking one `rbinom()` call would
#'   change the RNG stream and break bit-identical output
#'   (test_duration_optional_backcompat.R).
#' @param replicate_ids vector of replicate IDs to simulate in this task.
#' @param posterior_idx_occ length-R vector of posterior indices to use per
#'   replicate (paired strategy); will be clamped to draws_occ$S.
#' @param posterior_idx_dur same shape for the duration draws; ignored when
#'   draws_dur is NULL.
#' @param seed integer; root seed for the L'Ecuyer-CMRG stream.
#' @param output_part directory to write partition outputs into.
#' @param dt_months step size (default 1).
#' @param dry_run if TRUE, no files written.
#' @param enso_scenarios Optional named list of [M × n_rep] numeric matrices,
#'   one per lag name, from `sample_enso_forecast_trajectories()`. Column
#'   `rep_id` overrides that replicate's ENSO covariate, broadcast over all N
#'   rows. NULL (default) uses the panel values.
#' @param warmup_panels_list Optional list (length >= max(replicate_ids)) of
#'   per-replicate warmup mini-panels from `build_hybrid_warmup_panels()`,
#'   each passed to `build_warmup_state`. NULL (default) computes one shared
#'   warmup state from `panel_tu`.
#' @param force_o_tu logical; TRUE sets update_lags()'s o_tu to 1 on every
#'   valid-TU row (under partition_by = country the leave-one-country-out value
#'   is always 0), as the spatPomp rprocess variants do
#'   (`generate_rprocess_csnippet()`). NA-TU rows untouched. Default FALSE.
#' @param drop_country_tu_spillover logical; TRUE drops o_ct/o_tu state and LP
#'   terms (own `o` untouched) even if the fit has the coefficients
#'   (`simulator.drop_country_tu_spillover`). Takes precedence over
#'   `force_o_tu` (o_tu is 0). Default FALSE.
#' @param incidence_cap NULL (no cap) or a number > 0: monthly cases per
#'   location are held at incidence_cap x pop, on the mean (`mu_mon`) and on
#'   the drawn `cases_month`. Config-driven callers pass
#'   incidence_cap_of(opt$simulator) (simulator.incidence_cap, default 0.15).
#' @param onset_eligibility "inactive" (default) or "season": an inactive unit
#'   may draw an onset only if no onset that reset u_days fell in its current
#'   season (season_onset_open()). Needs draws_dur. Config-driven callers pass
#'   onset_eligibility_of(opt$simulator).
#' @param cessation_exposure "full" (default) or "half": the monthly case mean
#'   and draw follow the cessation draw, and a unit ceasing this month gets
#'   half its exposure (.ONSET_EXPOSURE_FRAC). Config-driven callers pass
#'   cessation_exposure_of(opt$simulator).
simulate_one_tu <- function(
  tu,
  frame_tu,
  panel_tu,
  draws_occ,
  replicate_ids,
  posterior_idx_occ = NULL,
  seed,
  output_part,
  horizon,
  dt_months = 1L,
  out_fmt = "qs",
  dry_run = FALSE,
  enso_scenarios = NULL,
  warmup_panels_list = NULL,
  cb_meta_occ = NULL,
  grouping = NULL,
  force_o_tu = FALSE,
  draws_dur = NULL,
  cb_meta_dur = NULL,
  posterior_idx_dur = NULL,
  draws_mon = NULL,
  cb_meta_mon = NULL,
  posterior_idx_mon = NULL,
  drop_country_tu_spillover = FALSE,
  use_nbr = FALSE,
  nbr_rho_km = 500L,
  nbr_gadm_path = "data/geodata/gadm_valid.gpkg",
  nbr_cache_dir = "generated_data/",
  nbr_spec = NULL,
  covariate_feed = "latent",
  delta_fn = 0.10,
  delta_fp = 0.05,
  clamp_enso_to_knots = FALSE,
  clamp_enso_margin = 0,
  incidence_cap = NULL,
  onset_eligibility = "inactive",
  cessation_exposure = "full"
) {
  stopifnot(nrow(frame_tu) > 0, !is.null(draws_occ))
  R <- length(replicate_ids)
  N <- nrow(frame_tu)
  # The duration-aware path is opt-in via a non-NULL draws_dur bundle.
  # See the @param doc above for the hard-branch rationale.
  has_dur <- !is.null(draws_dur)
  # Detected-activity covariate feed (update_detection()). Needs the duration
  # track: only `active` defines an episode whose detection persists.
  if (!covariate_feed %in% c("latent", "detected")) {
    stop("simulate_one_tu: covariate_feed must be 'latent' or 'detected'.")
  }
  use_det <- identical(covariate_feed, "detected")
  if (use_det && !has_dur) {
    stop("simulate_one_tu: covariate_feed = 'detected' requires draws_dur.")
  }
  # Onset eligibility (simulator.onset_eligibility). The season gate
  # sits on top of the duration-aware `active` gate, so it needs draws_dur.
  gate_season <- identical(
    check_onset_eligibility(
      onset_eligibility,
      "simulate_one_tu: onset_eligibility"
    ),
    "season"
  )
  if (gate_season && !has_dur) {
    stop("simulate_one_tu: onset_eligibility = 'season' requires draws_dur.")
  }
  # Cessation-month exposure (simulator.cessation_exposure): under
  # half the monthly case step runs after the cessation draw. No effect without
  # the mon track (nothing else reads the exposure of a ceasing month).
  cess_half <- identical(
    check_cessation_exposure(
      cessation_exposure,
      "simulate_one_tu: cessation_exposure"
    ),
    "half"
  )
  # Monthly case-count (`mon`) path: opt-in via a non-NULL draws_mon bundle.
  # The sole case estimand.
  has_mon <- !is.null(draws_mon)
  if (has_mon && !has_dur) {
    stop(
      "simulate_one_tu: the monthly case-count track requires duration draws. ",
      "`mso` and `active` are latent states produced by the duration model, ",
      "with no defensible fallback -- set `simulator.duration_fit` (or ",
      "`simulator.grouping.duration_enable`) alongside `simulator.monthly_fit`."
    )
  }

  # Country-grouping routing. When `grouping` is supplied, each
  # country in this TU is fit by its group's own occ (+ dur/mon) model (region
  # dropped from the RE nesting). The LP evaluation below is partitioned by
  # group: every group's rows are scored with that group's draws bundle,
  # cross-basis metadata and region-dropped RE key space, then scattered back
  # into the length-N vectors. `grouping` carries `draws_occ_by_group`,
  # `cb_meta_occ_by_group` and, per enabled track, `draws_{dur,mon}_by_group` /
  # `cb_meta_{dur,mon}_by_group` (named by group); the per-row
  # group label comes from `frame_tu$group`. Lag bookkeeping (update_lags) stays
  # TU-level so cross-group spillover within a TU is preserved. When NULL the
  # single-bundle path below is byte-identical to the pre-grouping kernel.
  is_grouped <- !is.null(grouping)
  if (is_grouped) {
    if (!("group" %in% names(frame_tu))) {
      stop(
        "simulate_one_tu: grouped routing requires a `group` column on frame_tu."
      )
    }
    grp_vec <- as.character(frame_tu$group)
    ctx_groups <- unique(grp_vec)
    if (any(is.na(grp_vec))) {
      stop(
        "simulate_one_tu: grouped routing found NA group labels in frame_tu."
      )
    }
  } else {
    grp_vec <- rep("__ALL__", N) # non-empty sentinel (valid list name)
    ctx_groups <- "__ALL__"
  }

  # Raw (pre-clamp) per-replicate draw index. Under grouping each group bundle
  # has its own S, so the clamp is applied per group inside the loop; the raw
  # index keeps the same replicate→draw pairing across all groups.
  raw_idx_occ <- if (is.null(posterior_idx_occ)) {
    replicate_ids
  } else {
    as.integer(posterior_idx_occ)
  }
  if (is.null(posterior_idx_occ)) {
    posterior_idx_occ <- ((replicate_ids - 1L) %% draws_occ$S) + 1L
  }
  # Clamp supplied indices to the available number of draws — the caller may
  # not know how many draws each backend exposes.
  posterior_idx_occ <- ((as.integer(posterior_idx_occ) - 1L) %% draws_occ$S) +
    1L

  # Duration draw index — same paired-strategy clamp, only meaningful when
  # has_dur. raw_idx_dur/posterior_idx_dur stay NULL on the legacy path so
  # nothing downstream can accidentally reference them.
  raw_idx_dur <- NULL
  if (has_dur) {
    raw_idx_dur <- if (is.null(posterior_idx_dur)) {
      replicate_ids
    } else {
      as.integer(posterior_idx_dur)
    }
    if (is.null(posterior_idx_dur)) {
      posterior_idx_dur <- ((replicate_ids - 1L) %% draws_dur$S) + 1L
    }
    posterior_idx_dur <- ((as.integer(posterior_idx_dur) - 1L) %% draws_dur$S) +
      1L
  }

  # Monthly draw index — same paired-strategy clamp as occ/dur.
  raw_idx_mon <- NULL
  if (has_mon) {
    raw_idx_mon <- if (is.null(posterior_idx_mon)) {
      replicate_ids
    } else {
      as.integer(posterior_idx_mon)
    }
    if (is.null(posterior_idx_mon)) {
      posterior_idx_mon <- ((replicate_ids - 1L) %% draws_mon$S) + 1L
    }
    posterior_idx_mon <- ((as.integer(posterior_idx_mon) - 1L) %% draws_mon$S) +
      1L
  }

  # Months sequence (one row per gid per month). The panel may be irregular;
  # we drive the simulation by the unique sorted t_start values per gid.
  # The warm-up tail (t_start < horizon start) is read by build_warmup_state
  # to seed the AR state but must NOT be simulated/emitted: restrict the
  # simulated grid to the horizon proper so events fall only inside [start, end].
  months <- sort(unique(panel_tu$t_start))
  months <- months[months >= as.Date(horizon$start)]
  M <- length(months)
  if (M == 0) {
    warning(
      sprintf(
        "[simulate_one_tu] TU '%s': no panel months at or after horizon start %s — ",
        tu,
        horizon$start
      ),
      "all panel rows fall within the warmup window. Skipping this TU."
    )
    return(invisible(list(
      event_table = tibble::tibble(),
      monthly_summary = tibble::tibble()
    )))
  }

  # Wide panel: rows = months, columns = climate vars stacked per gid.
  # Cheaper: index per (gid, month) directly via a key. Warm-up rows get
  # `.m_idx = NA` (their t_start is below the trimmed grid) and are dropped.
  panel_key <- panel_tu %>%
    dplyr::mutate(
      .gid_idx = match(as.character(gid), frame_tu$gid),
      .m_idx = match(t_start, .env$months)
    ) %>%
    dplyr::filter(!is.na(.gid_idx), !is.na(.m_idx))

  build_clim <- function(col, default = 0) {
    M_mat <- matrix(default, nrow = N, ncol = M)
    M_mat[cbind(panel_key$.gid_idx, panel_key$.m_idx)] <- panel_key[[col]]
    M_mat
  }
  climate <- list(
    mean_precip = build_clim("mean_precip"),
    dev_precip = build_clim("dev_precip"),
    mean_temp_daymean = build_clim("mean_temp_daymean"),
    dev_temp_daymean = build_clim("dev_temp_daymean")
  )
  # ENSO (`enso_lag<n>`) covariate(s), if present in the panel — fed to the
  # hazard LP under the same name as the brms FE column.
  # enso_lag<n> and enso_lead<k> (deployed ENSO-lead fits; base.R::enso_lag_col_name()).
  enso_cols <- grep("^enso_(lag|lead)[0-9]+$", names(panel_tu), value = TRUE)
  for (ec in enso_cols) {
    climate[[ec]] <- build_clim(ec)
  }
  # ENSO is a global (location-invariant) index: all gids share the same value
  # in a given month. build_clim() defaults absent gid-months to 0, so gids
  # whose panel data ends before the simulation horizon (e.g. countries with no
  # recent observed outbreaks) silently receive ENSO = 0. With a multi-year
  # contpanel warmup this accumulates into systematically wrong TTLO trajectories
  # for those countries. Fix: for each month, propagate the ENSO value from any
  # present gid to all N rows (take the element with largest |value| so the
  # correct non-zero value wins over the 0 default from absent gids).
  for (ec in enso_cols) {
    mat <- climate[[ec]]
    month_vals <- apply(mat, 2L, function(col) col[which.max(abs(col))])
    climate[[ec]] <- matrix(month_vals, nrow = N, ncol = M, byrow = TRUE)
  }
  # Location-specific climate variables (mean_precip, dev_precip,
  # mean_temp_daymean, dev_temp_daymean): gids whose panel coverage ends before
  # the simulation horizon receive 0 from build_clim(). Forward-fill each gid's
  # last observed value across trailing missing months so absent gids carry a
  # plausible climatological state rather than a spurious zero.
  # Note: current elag occurrence configs use only enso_lag<n>; this fix is
  # relevant for the duration/monthly models and future occurrence configs including these terms.
  {
    has_data <- matrix(FALSE, nrow = N, ncol = M)
    if (nrow(panel_key) > 0L) {
      has_data[cbind(panel_key$.gid_idx, panel_key$.m_idx)] <- TRUE
    }
    loc_clim_cols <- c(
      "mean_precip",
      "dev_precip",
      "mean_temp_daymean",
      "dev_temp_daymean"
    )
    for (lcc in loc_clim_cols) {
      if (!(lcc %in% names(climate))) {
        next
      }
      mat <- climate[[lcc]]
      for (i in seq_len(N)) {
        real_ms <- which(has_data[i, ])
        if (length(real_ms) == 0L || max(real_ms) == M) {
          next
        }
        last_m <- max(real_ms)
        mat[i, (last_m + 1L):M] <- mat[i, last_m]
      }
      climate[[lcc]] <- mat
    }
  }

  # DLNM cross-basis design --------------------------------------------------
  # Reconstruct `<var>_cb<k>` covariate matrices for any cross-basis fit, driven
  # entirely by the cb_meta sidecars (self-describing: argvar/arglag/lag/
  # cb_col_names). Each track (occ / dur / mon) is namespaced separately because
  # the fits may carry different bases for the same variable name. Two regimes:
  #   * climate vars  — replicate-invariant; built once here (below).
  #   * "enso"        — replicate-variant (the enso_lag* matrices may be swapped
  #                     per replicate by enso_scenarios); rebuilt in the loop.
  # `cb_meta` entry "enso" routes make_lag_matrix() to the enso_lag* branch; any
  # other name reads the climate `<var>_lag*` columns (emitted by 03c/03d).
  split_cb <- function(cb_meta) {
    if (is.null(cb_meta) || length(cb_meta) == 0) {
      return(list(enso = NULL, clim = NULL))
    }
    enso <- cb_meta[intersect("enso", names(cb_meta))]
    clim <- cb_meta[setdiff(names(cb_meta), "enso")]
    list(
      enso = if (length(enso)) enso else NULL,
      clim = if (length(clim)) clim else NULL
    )
  }
  # Build cross-basis N×M matrices from a set of lag matrices (one per lag,
  # lag0..max_lag) flattened column-major, then reshaped back. Returns a named
  # list keyed by cb_col_names. `lag_getter(col)` returns an N×M matrix for the
  # named lag column (build_clim for climate, climate[[..]] for enso).
  build_cb_mats <- function(
    meta_list,
    lag_getter,
    loc_invariant = FALSE,
    clamp = FALSE,
    clamp_margin = 0
  ) {
    mats <- list()
    if (is.null(meta_list)) {
      return(mats)
    }
    for (v in names(meta_list)) {
      # Latent cross-basis variables (utils/dlnm.R) are not panel columns and so
      # have no N×M matrix: `mso` differs per particle and is only known inside
      # the month loop. The kernel builds their basis with a direct
      # build_cb_columns_from_Q() lookup table instead (see B_mso below).
      if (v %in% .LATENT_CB_VARS) {
        next
      }
      meta <- meta_list[[v]]
      if (!is.null(meta$single_lag)) {
        sl <- as.integer(meta$single_lag)
        lag_cols <- if (identical(v, "enso")) {
          # negative = enso_lead<k>
          enso_lag_col_name(sl)
        } else if (sl == 0L) {
          v
        } else {
          sprintf("%s_lag%d", v, sl)
        }
        n_lag_cols <- 1L
      } else if (!is.null(meta$lags)) {
        stop(
          "build_cb_mats: multi-lag cross-basis meta (`lags`, brms_kl) is not supported ",
          "by the simulator kernel (variable '",
          v,
          "').",
          call. = FALSE
        )
      } else {
        max_lag <- as.integer(meta$lag[2])
        lag_cols <- if (identical(v, "enso")) {
          sprintf("enso_lag%d", seq(0L, max_lag))
        } else {
          c(v, sprintf("%s_lag%d", v, seq_len(max_lag)))
        }
        n_lag_cols <- max_lag + 1L
      }
      lag_mats <- lapply(lag_cols, lag_getter)
      if (loc_invariant) {
        # ENSO is location-invariant: transform M rows then broadcast.
        Q <- matrix(0, nrow = M, ncol = n_lag_cols)
        for (k in seq_along(lag_mats)) {
          Q[, k] <- lag_mats[[k]][1L, ]
        }
        Q[is.na(Q)] <- 0
        # simulator.enso.clamp_to_knots: hold ENSO beyond the fitted boundary
        # knots widened by simulator.enso.clamp_margin (.clamp_cb_range(), dlnm.R).
        if (clamp) {
          Q <- .clamp_cb_range(Q, meta$range, clamp_margin)
        }
        cb <- build_cb_columns_from_Q(Q, meta) # M × ncb
        for (j in seq_len(ncol(cb))) {
          mats[[meta$cb_col_names[j]]] <-
            matrix(cb[, j], nrow = N, ncol = M, byrow = TRUE)
        }
      } else {
        Q <- matrix(0, nrow = N * M, ncol = n_lag_cols)
        for (k in seq_along(lag_mats)) {
          Q[, k] <- as.numeric(lag_mats[[k]])
        }
        Q[is.na(Q)] <- 0
        cb <- build_cb_columns_from_Q(Q, meta) # (N*M) × ncb
        for (j in seq_len(ncol(cb))) {
          mats[[meta$cb_col_names[j]]] <- matrix(cb[, j], nrow = N, ncol = M)
        }
      }
    }
    mats
  }

  # LP contexts. One per group under grouping; a single context
  # spanning all rows otherwise. Each context owns its draws bundles, cross-basis
  # split metadata, and replicate-invariant climate cross-basis matrices (built
  # once here, indexed by the context's rows at LP time). ENSO cross-basis is
  # replicate-variant and rebuilt inside the loop per context. The single-context
  # case reproduces the pre-grouping `clim_cb_occ` + `occ_split`
  # exactly, so the math is unchanged when `grouping` is NULL.
  make_context <- function(gname) {
    idx <- if (is_grouped) which(grp_vec == gname) else seq_len(N)
    if (is_grouped) {
      d_occ <- grouping$draws_occ_by_group[[gname]]
      cbm_occ <- grouping$cb_meta_occ_by_group[[gname]]
      if (is.null(d_occ)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing occ draws for group '%s'.",
          gname
        ))
      }
    } else {
      d_occ <- draws_occ
      cbm_occ <- cb_meta_occ
    }
    o_split <- split_cb(cbm_occ)
    # Duration draws AND cross-basis, routed per group exactly like occ.
    # `04_fit_dur_brms_cb.R` fits carry `enso_cb<k>` fixed effects, and the
    # duration track's best lag is chosen independently of the other tracks
    # (Horn: occ=11, siz=9, mon=11, dur=8), so the duration basis must be built
    # from its own cb_meta -- reusing the occurrence columns would evaluate the
    # duration coefficients on the wrong lag and the wrong spline basis.
    # `d_split` is NULL when the duration fit carries no cross-basis, which
    # makes build_cb_mats() a no-op and reproduces the pre-fix behaviour.
    d_dur <- NULL
    d_split <- NULL
    if (has_dur) {
      d_dur <- if (is_grouped) {
        grouping$draws_dur_by_group[[gname]]
      } else {
        draws_dur
      }
      if (is_grouped && is.null(d_dur)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing duration draws for group '%s'.",
          gname
        ))
      }
      cbm_dur <- if (is_grouped) {
        grouping$cb_meta_dur_by_group[[gname]]
      } else {
        cb_meta_dur
      }
      d_split <- split_cb(cbm_dur)
    }
    # Monthly bundle + its cross-basis, routed per group exactly like
    # occ. `mso` is stripped out of `m_split$clim` by build_cb_mats()'s
    # .LATENT_CB_VARS skip and handled by the B_mso lookup table instead.
    d_mon <- NULL
    m_split <- NULL
    mso_tab <- NULL
    if (has_mon) {
      d_mon <- if (is_grouped) {
        grouping$draws_mon_by_group[[gname]]
      } else {
        draws_mon
      }
      cbm_mon <- if (is_grouped) {
        grouping$cb_meta_mon_by_group[[gname]]
      } else {
        cb_meta_mon
      }
      if (is.null(d_mon)) {
        stop(sprintf(
          "simulate_one_tu: grouped routing missing monthly draws for group '%s'.",
          gname
        ))
      }
      m_split <- split_cb(cbm_mon)
      # NULL when the fit has no mso basis; else $B is (mso_max + 1) x K, rows
      # indexed by mso = 0..mso_max (read by .inject_mso()).
      mso_tab <- build_mso_basis(cbm_mon)
    }
    list(
      gname = gname,
      idx = idx,
      frame_g = frame_tu[idx, , drop = FALSE],
      draws_occ = d_occ,
      draws_dur = d_dur,
      draws_mon = d_mon,
      drop_region = is_grouped,
      occ_split = o_split,
      dur_split = d_split,
      mon_split = m_split,
      mso_tab = mso_tab,
      clim_cb_occ = build_cb_mats(
        o_split$clim,
        build_clim,
        loc_invariant = FALSE
      ),
      clim_cb_dur = if (has_dur) {
        build_cb_mats(d_split$clim, build_clim, loc_invariant = FALSE)
      } else {
        NULL
      },
      clim_cb_mon = if (has_mon) {
        build_cb_mats(m_split$clim, build_clim, loc_invariant = FALSE)
      } else {
        NULL
      }
    )
  }
  contexts <- lapply(ctx_groups, make_context)
  names(contexts) <- ctx_groups

  # t_start_seas is purely deterministic: yday(t_start shifted by min_seas) - 1.
  # Reading it from the panel (build_clim) defaults to 0 for gid-months that
  # lack a panel row (e.g. gids whose data starts mid-horizon), which forces the
  # seasonal spline to evaluate at the season-start value for all those months
  # and washes out the seasonal signal in the aggregate. Instead, compute it
  # from frame_tu$min_seas and the months grid — always correct.
  # Using lapply over gids (rows) to avoid outer() coercing the Date vector to
  # numeric, which would require an origin argument for as.Date().
  t_start_seas_mat <- do.call(
    rbind,
    lapply(as.integer(frame_tu$min_seas), function(ms) {
      as.integer(
        lubridate::yday(months %m+% lubridate::period(-(ms - 1L), "months")) -
          1L
      )
    })
  ) # N × M
  # Month-of-season index for the season onset gate (NULL under "inactive").
  seas_k_mat <- if (gate_season) {
    season_month_index(frame_tu$min_seas, months)
  } else {
    NULL
  }

  # Kernel speed-up: per context and track, the evaluate_lp() lookup plan
  # (lp_plan(): RE indices resolved once, coefficients cached per draw) and the
  # seasonal-spline design of every month (smooth_X; t_start_seas is fixed per
  # row and month). Both are pure functions of the bundle and the rows, so the
  # LP is bit-identical to building them inside every call.
  .smooth_designs <- function(d, idx) {
    nms <- if (length(d$smooth_pred_fn) > 0L && length(d$smooth_coef) > 0L) {
      names(d$smooth_pred_fn)
    } else {
      character(0)
    }
    if (!length(nms)) {
      return(NULL)
    }
    lapply(seq_len(M), function(m) {
      t_ss <- t_start_seas_mat[idx, m]
      out <- lapply(nms, function(nm) {
        f <- d$smooth_pred_fn[[nm]]
        if (is.null(f)) NULL else f(t_ss)
      })
      stats::setNames(out, nms)
    })
  }
  for (g in names(contexts)) {
    ctx <- contexts[[g]]
    ctx$plan_occ <- lp_plan(ctx$draws_occ, ctx$frame_g, ctx$drop_region)
    ctx$smooth_occ <- .smooth_designs(ctx$draws_occ, ctx$idx)
    if (has_dur) {
      ctx$plan_dur <- lp_plan(ctx$draws_dur, ctx$frame_g, ctx$drop_region)
      ctx$smooth_dur <- .smooth_designs(ctx$draws_dur, ctx$idx)
    }
    if (has_mon) {
      ctx$plan_mon <- lp_plan(ctx$draws_mon, ctx$frame_g, ctx$drop_region)
      ctx$smooth_mon <- .smooth_designs(ctx$draws_mon, ctx$idx)
    }
    # The covariate columns each LP actually reads (evaluate_lp() takes only its
    # linear terms and t_start_seas from xrow); the per-step lists below carry
    # just these instead of every climate / enso_lag column of the panel.
    ctx$need_occ <- c(ctx$plan_occ$linear_terms, "t_start_seas")
    ctx$need_dur <- if (has_dur) {
      c(ctx$plan_dur$linear_terms, "t_start_seas")
    } else {
      NULL
    }
    ctx$need_mon <- if (has_mon) {
      c(ctx$plan_mon$linear_terms, "t_start_seas")
    } else {
      NULL
    }
    contexts[[g]] <- ctx
  }
  clim_need <- intersect(
    names(climate),
    unique(unlist(lapply(contexts, function(ctx) {
      c(ctx$need_occ, ctx$need_dur, ctx$need_mon)
    })))
  )

  # Exposure offset — mirrors `offset(log_offset)` in the occurrence fit:
  #   exposure = pmax(t_end_seas - t_start_seas, 1L)
  #   log_offset = log(exposure)
  # brms parameters encode log-rate per day; the probability kernel needs
  # lambda = exp(eta) * exposure_days, not exp(eta) * dt_months.
  # Exposure = days_in_month(t_start) — the same for all gids in a given month.
  log_offset_mat <- matrix(
    log(as.integer(lubridate::days_in_month(months))),
    nrow = N,
    ncol = M,
    byrow = TRUE
  )

  # Shared warmup state: computed once when warmup_panels_list is NULL
  # (backward-compatible default); otherwise deferred to per-replicate calls inside
  # the loop so each replicate can start from a distinct importance-sampled
  # simulation trajectory.
  # TTLO level "none": seeded only when the occurrence fit
  # carries it; NULL otherwise.
  npo <- .ttlo_no_prior_for_frame(
    frame_tu,
    horizon$start,
    .draws_has_ttlo_none(
      if (!is.null(grouping)) grouping$draws_occ_by_group else draws_occ
    )
  )
  shared_warm <- if (is.null(warmup_panels_list)) {
    build_warmup_state(
      panel_tu,
      frame_tu,
      horizon,
      has_dur = has_dur,
      no_prior_onset = npo
    )
  } else {
    NULL
  }

  # Output accumulators
  event_rows <- vector("list", R)
  monthly_Y <- array(0L, dim = c(R, N, M))
  # Per-(replicate, location, month) "is this location in an active
  # episode this month" indicator. NULL on the legacy path (has_dur = FALSE)
  # so it costs nothing and never enters monthly_summary's column set.
  monthly_A <- if (has_dur) array(0L, dim = c(R, N, M)) else NULL
  # covariate_feed = "detected": the detected-activity field per month (the
  # `det` state after update_detection()); NULL otherwise, so latent runs keep
  # their monthly_summary column set.
  monthly_D <- if (use_det) array(0L, dim = c(R, N, M)) else NULL
  # Monthly case-count (mon) track: mu and the NegBin case draw at every ACTIVE
  # month. NULL on the legacy path for the
  # same reason as monthly_A. Left NA (not 0) where a location is inactive, so
  # "no episode" is distinguishable from "active episode that drew zero cases".
  monthly_MU <- if (has_mon) array(NA_real_, dim = c(R, N, M)) else NULL
  # Per-location monthly incidence ceiling (cap x pop), frame_tu row order;
  # Inf everywhere when incidence_cap is NULL (direct callers, tests).
  mu_cap <- if (has_mon) monthly_mu_cap(frame_tu$pop, incidence_cap) else NULL
  # pop_year mode: the U x M per-month population and the
  # per-month ceiling. NULL for a frame without pop_y<YYYY> columns, which
  # keeps the static `pop` / `mu_cap` above (the default, unchanged).
  pop_mat <- if (has_mon) frame_pop_matrix(frame_tu, months) else NULL
  mu_cap_mat <- if (!is.null(pop_mat)) {
    matrix(
      monthly_mu_cap(as.vector(pop_mat), incidence_cap),
      nrow = nrow(pop_mat)
    )
  } else {
    NULL
  }
  monthly_C <- if (has_mon) array(NA_integer_, dim = c(R, N, M)) else NULL

  # Use L'Ecuyer-CMRG streams seeded from (seed, replicate_id) — independent of
  # which array slot a replicate lands in (S7).
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1]), add = TRUE)

  # Snapshot the baseline ENSO covariate matrices BEFORE the replicate loop.
  # The per-replicate override overwrites `climate[[ec]]` in place;
  # without this snapshot the mutation would persist into the next replicate,
  # so a replicate that does not receive an override (e.g. rep_id beyond the
  # available scenario columns) would silently reuse the previous replicate's
  # ENSO trajectory instead of the baseline panel.
  enso_base <- if (!is.null(enso_scenarios)) {
    climate[intersect(names(enso_scenarios), names(climate))]
  } else {
    NULL
  }

  # Neighbour-contagion kernel (within-partition / within-country). Built once
  # per TU, replicate-invariant, only when `use_nbr` (flag on AND the occurrence
  # fit carries an nbr_* FE, resolved by the caller). `K_adm1_nbr` is a dense
  # n_adm1 x n_adm1 matrix in `unique(frame_tu$adm1_unit)` order with zero
  # diagonal; each step reconstructs the lag-1 ADM1-active vector from the `o`
  # state and forms `K %*% Aa1`. use_nbr=FALSE ⇒ NULL ⇒ every nbr block below is
  # a no-op (nbr stays NULL into evaluate_lp).
  K_adm1_nbr <- if (isTRUE(use_nbr)) {
    # Spec resolved by the driver (resolve_nbr_spec()); legacy callers that pass
    # only the scalar nbr_* args get the same resolution here, so the kernel's
    # scale/partition always follow the fit's own nbr FE name.
    nbr_spec <- nbr_spec %||%
      resolve_nbr_spec(
        list(
          enable = TRUE,
          rho_km = nbr_rho_km,
          gadm_path = nbr_gadm_path,
          cache_dir = nbr_cache_dir
        ),
        draws_occ,
        grouping
      )
    nbr_kernel_for_frame(frame_tu, nbr_spec)
  } else {
    NULL
  }

  # Integer codes of the lag aggregation (update_lags(maps =)) and, with the
  # neighbour kernel on, the positions linking its ADM1 order to the frame.
  lag_maps <- build_lag_maps(frame_tu)
  if (!is.null(K_adm1_nbr)) {
    nbr_k_pos <- match(rownames(K_adm1_nbr), lag_maps$adm1_lv)
    nbr_loc_pos <- match(as.character(frame_tu$adm1_unit), rownames(K_adm1_nbr))
  }

  # Monthly case-count step (`mon` track) for month m of the current replicate:
  # the log mean of every unit, the capped mean and the NegBin draw on the
  # units active this month. Called before update_lags() and the cessation draw
  # under simulator.cessation_exposure: full (the historical position), after
  # them under half, in both cases with the month's pre-update lags and its
  # pre-cessation active / mso. `ceased` (half only): TRUE where the unit's
  # cessation was drawn this month. Returns the active indices, their mean and
  # their drawn cases.
  mon_case_step <- function(
    m,
    active,
    mso,
    o,
    o_ct,
    o_tu,
    u_days,
    ceased = NULL
  ) {
    log_mu_mon <- numeric(N)
    phi_mon_vec <- numeric(N)

    # Latent mso cross-basis: one row of the per-TU lookup table per
    # location, clamped into the fitted boundary knots (build_mso_basis()).
    # Inactive locations (mso = 0) are scored too -- the result is simply
    # discarded by the `active == 1` gate below -- so the vectorised index
    # never sees an out-of-range value.
    xrow_mon <- xrow
    xrow_mon[["month_num"]] <- rep(as.integer(lubridate::month(months[m])), N)
    # monthly.seasonal_term: "harmonic2" (the default) names these four
    # columns; load_monthly_data() (utils/mon_data.R) builds them the
    # same way at fit time. Harmless when the fit used seasonal_term: none
    # -- evaluate_lp() drops any term absent from `fe`.
    mnum <- as.numeric(lubridate::month(months[m]))
    xrow_mon[["sin_month"]] <- rep(sin(2 * pi * mnum / 12), N)
    xrow_mon[["cos_month"]] <- rep(cos(2 * pi * mnum / 12), N)
    xrow_mon[["sin2_month"]] <- rep(sin(4 * pi * mnum / 12), N)
    xrow_mon[["cos2_month"]] <- rep(cos(4 * pi * mnum / 12), N)
    # monthly.admin_level: "ADM1-ADM2" makes admin_level a 2-level factor at
    # fit time; brms names the dummy "admin_levelADM2". Absent from the
    # frame (or from `fe`) => contributes nothing. Keyed on the CANONICAL
    # level, because the factor the fit levelled was the P3-re-keyed one
    # (.frame_admin_level_canonical()).
    if ("admin_level" %in% names(frame_tu)) {
      xrow_mon[["admin_levelADM2"]] <-
        as.numeric(.frame_admin_level_canonical(frame_tu) == "ADM2")
    }

    for (ctx in contexts) {
      idx <- ctx$idx
      xm <- lapply(
        xrow_mon[intersect(names(xrow_mon), ctx$need_mon)],
        function(v) v[idx]
      )
      cbm <- cb_mats_mon_byctx[[ctx$gname]]
      for (cn in names(cbm)) {
        xm[[cn]] <- cbm[[cn]][idx, m]
      }
      # Offsets are read off the bundle (monthly_log_mu()). In forecast mode
      # there is no observed `exposure` column, so the monthly exposure is
      # synthesised as log(days_in_month(t)) -- the same log_offset_mat the
      # occurrence and duration kernels already use, and the same quantity
      # the fit's panel carries for an INTERIOR month. The ONSET month is
      # partial in reality and is corrected to half a month (mso was
      # promoted to 1 above, at the new_onset block, so this reads the right
      # value). The CESSATION month keeps full-month exposure under
      # cessation_exposure: full (called before the cessation draw, `ceased`
      # NULL); under half this runs after the draw and `ceased` halves it.
      # See .ONSET_EXPOSURE_FRAC.
      lm_g <- monthly_log_mu(
        ctx$draws_mon,
        s_mon_byctx[[ctx$gname]],
        xm,
        o[idx],
        o_ct[idx],
        o_tu[idx],
        u_days[idx],
        ctx$frame_g,
        log_offset = log_offset_mat[idx, m],
        mso = mso[idx],
        mso_tab = ctx$mso_tab,
        what = "monthly kernel",
        drop_region = ctx$drop_region,
        drop_ct_tu = drop_country_tu_spillover,
        pop = if (!is.null(pop_mat)) pop_mat[idx, m] else NULL,
        plan = ctx$plan_mon,
        smooth_X = ctx$smooth_mon[[m]],
        ceased = if (is.null(ceased)) NULL else ceased[idx]
      )
      log_mu_mon[idx] <- lm_g
      phi_mon_vec[idx] <- phi_mon_byctx[[ctx$gname]]
    }

    # Monthly incidence cap (simulator.incidence_cap): the mean, and below
    # the drawn count, are held at cap x pop. Inf when incidence_cap is NULL.
    if (!is.null(mu_cap_mat)) {
      mu_cap <- mu_cap_mat[, m]
    }
    mu_mon <- exp(log_mu_mon)
    if (any(is.finite(mu_cap))) {
      mu_mon <- pmin(mu_mon, mu_cap)
    }
    act_i <- which(active == 1L)
    list(
      act_i = act_i,
      mu = mu_mon[act_i],
      cm = draw_monthly_cases(mu_mon[act_i], phi_mon_vec[act_i], mu_cap[act_i])
    )
  }

  for (r_idx in seq_len(R)) {
    rep_id <- replicate_ids[r_idx]
    set.seed(seed + rep_id)

    # Per-context posterior draw indices (each group bundle has its own S; clamp
    # the shared raw index to that group's draw count; one context gives
    # posterior_idx_occ[r_idx] exactly). Also scatter to
    # length-N per-row vectors for the event-table `posterior_draw_*` columns.
    s_occ_byctx <- vapply(
      contexts,
      function(ctx) {
        as.integer((raw_idx_occ[r_idx] - 1L) %% ctx$draws_occ$S) + 1L
      },
      integer(1)
    )
    s_occ_byrow <- integer(N)
    for (ctx in contexts) {
      s_occ_byrow[ctx$idx] <- s_occ_byctx[[ctx$gname]]
    }
    # Duration draw index per context/row — only meaningful when has_dur.
    s_dur_byctx <- NULL
    s_dur_byrow <- integer(N)
    if (has_dur) {
      s_dur_byctx <- vapply(
        contexts,
        function(ctx) {
          as.integer((raw_idx_dur[r_idx] - 1L) %% ctx$draws_dur$S) + 1L
        },
        integer(1)
      )
      for (ctx in contexts) {
        s_dur_byrow[ctx$idx] <- s_dur_byctx[[ctx$gname]]
      }
    }
    # Monthly draw index per context — only meaningful when has_mon.
    s_mon_byctx <- NULL
    if (has_mon) {
      s_mon_byctx <- vapply(
        contexts,
        function(ctx) {
          as.integer((raw_idx_mon[r_idx] - 1L) %% ctx$draws_mon$S) + 1L
        },
        integer(1)
      )
      # NegBin shape per location: a function of the draw and the rows only, so
      # computed once per replicate rather than every active month.
      phi_mon_byctx <- lapply(contexts, function(ctx) {
        compute_phi_vec(
          ctx$draws_mon,
          s_mon_byctx[[ctx$gname]],
          ctx$frame_g,
          drop_region = ctx$drop_region
        )
      })
    }

    # Per-replicate warmup state.
    # When warmup_panels_list is provided, each replicate draws its initial AR
    # state from the importance-sampled mini-panel for that global replicate ID.
    # Indexed by rep_id (global), not r_idx (within-chunk position), so chunks
    # always draw from the correct element of the full n_rep list.
    warm <- if (!is.null(warmup_panels_list)) {
      build_warmup_state(
        warmup_panels_list[[rep_id]],
        frame_tu,
        horizon,
        has_dur = has_dur,
        no_prior_onset = npo
      )
    } else {
      shared_warm
    }

    o <- warm$o
    o_ct <- warm$o_ct
    o_tu <- warm$o_tu
    if (isTRUE(drop_country_tu_spillover)) {
      o_ct[!is.na(o_ct)] <- 0L
      o_tu[!is.na(o_tu)] <- 0L
    }
    u_days <- warm$u # numeric days vector (step 4)

    # Duration-aware state (NULL/unused on the legacy path). `active`/
    # `mso` seed from build_warmup_state() (currently always 0 — see that
    # function's doc comment on the smoothing-cache handoff).
    # `open_row` carries the row index
    # (within this replicate's event vectors, ev_*) of each active location's
    # not-yet-closed event row, so cessation/censoring can patch it in place.
    active <- if (has_dur) as.integer(warm$active) else NULL
    mso <- if (has_dur) as.integer(warm$mso) else NULL
    det <- if (use_det) as.integer(warm$det %||% warm$active) else NULL
    open_row <- if (has_dur) rep(NA_integer_, N) else NULL
    # This replicate's event rows, kept as column vectors (row = position) and
    # turned into one data.frame after the month loop (kernel speed-up; the
    # per-month rbind() it replaces built the same rows in the same order).
    ev_loc <- integer(0)
    ev_m <- integer(0)
    ev_dur <- integer(0)
    ev_end <- integer(0)
    ev_cens <- integer(0)

    # Per-replicate ENSO override.
    # Replace the shared panel's ENSO covariate matrices with the replicate-
    # specific forecast trajectory from `enso_scenarios`. ENSO is location-
    # invariant, so the M-length vector is broadcast to all N rows. Indexing
    # uses rep_id (global replicate ID) so each chunk accesses the correct column
    # of the [M × n_rep] matrix regardless of its within-chunk position.
    if (!is.null(enso_scenarios)) {
      for (ec in names(enso_scenarios)) {
        if (!(ec %in% names(climate))) {
          next
        }
        # Restore the baseline ENSO matrix first so a previous replicate's
        # override never leaks into this one.
        climate[[ec]] <- enso_base[[ec]]
        if (rep_id > ncol(enso_scenarios[[ec]])) {
          stop(sprintf(
            "simulate_one_tu: replicate id %d exceeds available ENSO scenario columns (%d) for '%s' — n_replicates and the ENSO scenario matrix are out of sync.",
            rep_id,
            ncol(enso_scenarios[[ec]]),
            ec
          ))
        }
        scenario_vec <- enso_scenarios[[ec]][, rep_id]
        # Row 1 is the horizon start; the scenario runs to the end of the shared
        # panel, which can be past this TU's horizon. Keep the first M months:
        # recycling a longer vector row-wise would give each unit a different
        # (time-shifted) slice.
        if (length(scenario_vec) < M) {
          stop(sprintf(
            "simulate_one_tu: ENSO scenario '%s' has %d months but TU '%s' simulates %d.",
            ec,
            length(scenario_vec),
            tu,
            M
          ))
        }
        scenario_vec <- scenario_vec[seq_len(M)]
        climate[[ec]] <- matrix(scenario_vec, nrow = N, ncol = M, byrow = TRUE)
      }
    }

    # ENSO cross-basis: rebuilt per replicate from the (possibly scenario-swapped)
    # enso_lag* matrices. Location-invariant, so transform once per month then
    # broadcast. Combined with each context's replicate-invariant climate
    # cross-basis. Under grouping, every group may carry its own basis/max_lag,
    # so the cb columns are namespaced per context.
    enso_getter <- function(col) climate[[col]]
    cb_mats_occ_byctx <- lapply(contexts, function(ctx) {
      c(
        ctx$clim_cb_occ,
        build_cb_mats(
          ctx$occ_split$enso,
          enso_getter,
          loc_invariant = TRUE,
          clamp = clamp_enso_to_knots,
          clamp_margin = clamp_enso_margin
        )
      )
    })
    cb_mats_mon_byctx <- if (!has_mon) {
      NULL
    } else {
      lapply(contexts, function(ctx) {
        c(
          ctx$clim_cb_mon,
          build_cb_mats(
            ctx$mon_split$enso,
            enso_getter,
            loc_invariant = TRUE,
            clamp = clamp_enso_to_knots,
            clamp_margin = clamp_enso_margin
          )
        )
      })
    }
    cb_mats_dur_byctx <- if (!has_dur) {
      NULL
    } else {
      lapply(contexts, function(ctx) {
        c(
          ctx$clim_cb_dur,
          build_cb_mats(
            ctx$dur_split$enso,
            enso_getter,
            loc_invariant = TRUE,
            clamp = clamp_enso_to_knots,
            clamp_margin = clamp_enso_margin
          )
        )
      })
    }

    for (m in seq_len(M)) {
      # One named covariate vector per climate/enso column for this month
      # (full-N; shared across contexts since the climate panel is TU-wide).
      xrow <- lapply(climate[clim_need], function(M_mat) M_mat[, m])
      xrow[["t_start_seas"]] <- t_start_seas_mat[, m]

      # Neighbour-contagion pressure for this step. `o` entering the step IS the
      # lag-1 ADM1-active indicator broadcast to each location (seeded from the
      # warmup on step 1, carried from the previous step's update_lags otherwise,
      # and already keyed on `active` not `Y` on the duration path). Collapse it
      # back to one value per ADM1 (in kernel order), apply the kernel, broadcast
      # to locations by adm1_unit (ADM2 inherits its parent). NULL when !use_nbr.
      nbr_loc <- NULL
      if (isTRUE(use_nbr)) {
        o_i <- as.integer(o)
        Aa1 <- if (.is01(o_i)) {
          .adm1_any(o_i, lag_maps)[nbr_k_pos]
        } else {
          tapply(o_i, frame_tu$adm1_unit, max, default = 0L)[rownames(
            K_adm1_nbr
          )]
        }
        Aa1 <- as.numeric(Aa1)
        Aa1[is.na(Aa1)] <- 0
        nbr_adm1 <- nbr_from_adm1_active(Aa1, K_adm1_nbr)
        nbr_loc <- as.numeric(nbr_adm1)[nbr_loc_pos]
        nbr_loc[is.na(nbr_loc)] <- 0
      }

      # Occurrence LP, partitioned by context: each group's rows are scored with
      # that group's draws bundle, cross-basis columns, posterior draw and
      # region-dropped RE key space, then scattered into the length-N log_lambda.
      log_lambda <- numeric(N)
      for (ctx in contexts) {
        idx <- ctx$idx
        xo <- lapply(xrow[intersect(names(xrow), ctx$need_occ)], function(v) {
          v[idx]
        })
        cbo <- cb_mats_occ_byctx[[ctx$gname]]
        for (cn in names(cbo)) {
          xo[[cn]] <- cbo[[cn]][idx, m]
        }
        bh_g <- basehaz_eval(
          ctx$draws_occ,
          t_start_seas_mat[idx, m],
          length(idx)
        )
        log_lambda[idx] <- evaluate_lp(
          ctx$draws_occ,
          s_occ_byctx[[ctx$gname]],
          xo,
          o[idx],
          o_ct[idx],
          o_tu[idx],
          u_days[idx],
          ctx$frame_g,
          basehaz_val = bh_g,
          drop_region = ctx$drop_region,
          drop_ct_tu = drop_country_tu_spillover,
          nbr = if (is.null(nbr_loc)) NULL else nbr_loc[idx],
          plan = ctx$plan_occ,
          smooth_X = ctx$smooth_occ[[m]]
        ) +
          # Empirical baseline-hazard anchor, offset(log(base_rate_occ)). Its
          # denominator is exposure days -- exactly the log_offset_mat term
          # added below -- so with the anchor on, exp(Intercept) is a residual
          # near 1 rather than a whole-corpus hazard scale. Identically 0 for
          # any fit made without it.
          .anchor_log_offset(ctx$draws_occ, frame_tu, idx, "occurrence kernel")
      }
      # Per-step probability: 1 - exp(-λ) via -expm1(-λ).
      # λ = rate_per_day × exposure_days; the offset (log exposure) is added to
      # the LP rather than multiplying exp(eta) by dt_months, because the brms
      # parameters are on the per-day rate scale (offset(log_offset) in the fit).
      lambda <- exp(log_lambda + log_offset_mat[, m])
      p <- -expm1(-lambda)

      # Onset draw. Legacy path: unconditional per-step Bernoulli, exactly as
      # before (byte-identical RNG stream). Duration-aware path: gated on
      # not-already-active (`at_risk`) so a location that is mid-episode
      # cannot "re-onset" every month — see the @param draws_dur doc above.
      # Hard branch,
      # not a masked single call, so the legacy RNG stream is untouched.
      active_before <- if (has_dur) active else NULL
      if (has_dur) {
        at_risk <- active_before == 0L
        # Season gate (onset_eligibility: season): no onset that reset u_days
        # in this unit's current season. Same expression as the C kernel.
        if (gate_season) {
          at_risk <- at_risk & season_onset_open(u_days, seas_k_mat[, m])
        }
        Y <- integer(N)
        if (any(at_risk)) {
          Y[at_risk] <- as.integer(stats::rbinom(sum(at_risk), 1L, p[at_risk]))
        }
      } else {
        Y <- as.integer(stats::rbinom(N, 1L, p))
      }

      monthly_Y[r_idx, , m] <- Y

      if (has_dur) {
        # active/mso bookkeeping. new_onset locations start mso = 1 this
        # month; already-active locations (active_before == 1, necessarily
        # Y == 0 there since at_risk excluded them) advance mso by one.
        new_onset <- which(Y == 1L)
        still_active <- which(active_before == 1L)
        mso[still_active] <- mso[still_active] + 1L
        mso[new_onset] <- 1L
        active[new_onset] <- 1L
        monthly_A[r_idx, , m] <- active
      }

      if (any(Y == 1L)) {
        # Event rows
        idx_trigger <- which(Y == 1L)
        if (has_dur) {
          # duration_months/end_month/censored are placeholders here, patched
          # in place (by row index, via `open_row`) at cessation or at the
          # end-of-horizon censoring pass below — the row is created at onset,
          # before the episode's length is known.
          nrow_before <- length(ev_loc)
          k <- length(idx_trigger)
          ev_loc <- c(ev_loc, idx_trigger)
          ev_m <- c(ev_m, rep(m, k))
          ev_dur <- c(ev_dur, rep(NA_integer_, k)) # filled at cessation/censoring
          ev_end <- c(ev_end, rep(NA_integer_, k))
          ev_cens <- c(ev_cens, rep(NA_integer_, k))
          open_row[idx_trigger] <- nrow_before + seq_along(idx_trigger)
        } else {
          ev_loc <- c(ev_loc, idx_trigger)
          ev_m <- c(ev_m, rep(m, length(idx_trigger)))
        }
      }

      # Monthly case-count draw (`mon` track). Placement is load-bearing:
      # AFTER promotion (an onsetting unit emits its first month, mso = 1) and
      # BEFORE update_lags() and the cessation draw (spillover lags match the
      # occurrence LP; a ceasing unit emits its final month). Under
      # cessation_exposure: half the call moves after the cessation draw, fed
      # saved pre-update values. Offsets come from the bundle; mu is gated on
      # `active`, not on `Y`.
      if (has_mon && !cess_half && any(active == 1L)) {
        mc <- mon_case_step(m, active, mso, o, o_ct, o_tu, u_days)
        monthly_MU[r_idx, mc$act_i, m] <- mc$mu
        ok <- !is.na(mc$cm)
        if (any(ok)) monthly_C[r_idx, mc$act_i[ok], m] <- mc$cm[ok]
      }
      # cessation_exposure: half -- this month's activity, mso and pre-update
      # lags, kept for mon_case_step() after the cessation draw below.
      if (has_mon && cess_half) {
        act_m <- active
        mso_m <- mso
        o_m <- o
        o_ct_m <- o_ct
        o_tu_m <- o_tu
        u_days_m <- u_days
      }

      # Update state for next step. Duration-aware path feeds `active`
      # (ongoing-episode indicator) instead of `Y` (onset-only indicator) to
      # update_lags(), so the contagion covariates reflect the true duty
      # cycle of an outbreak rather than firing for exactly the onset month
      # — see the @param draws_dur doc above ("the contagion-covariate
      # duty-cycle bug"). advance_ttlo() stays keyed on Y on both paths: TTLO
      # resets at onset regardless of how long the episode then runs.
      # covariate_feed = "detected": spillover covariates and TTLO follow the
      # detected-activity field (update_detection()), not the latent one.
      if (use_det) {
        det <- update_detection(Y, active, det, delta_fn, delta_fp)
        monthly_D[r_idx, , m] <- det
      }
      lags <- update_lags(
        if (use_det) {
          det
        } else if (has_dur) {
          active
        } else {
          Y
        },
        frame_tu,
        maps = lag_maps
      )
      if (isTRUE(force_o_tu)) {
        lags$o_tu <- ifelse(is.na(lags$o_tu), NA_integer_, 1L)
      }
      # drop_country_tu_spillover takes PRECEDENCE over force_o_tu: o_ct/o_tu
      # are forced to 0 (skipping the leave-one-out aggregation entirely),
      # not 1, when both flags are set.
      if (isTRUE(drop_country_tu_spillover)) {
        lags$o_ct <- ifelse(is.na(lags$o_ct), NA_integer_, 0L)
        lags$o_tu <- ifelse(is.na(lags$o_tu), NA_integer_, 0L)
      }
      o <- lags$o
      o_ct <- lags$o_ct
      o_tu <- lags$o_tu
      u_days <- advance_ttlo(u_days, if (use_det) Y * det else Y)

      if (has_dur) {
        # Cessation hazard, evaluated for every context-partitioned row (like
        # the occurrence LP above) but only drawn/applied for currently
        # active locations. xrow_dur reuses this month's climate covariates;
        # mso feeds the primary months-since-onset baseline hazard via
        # basehaz_dur_eval().
        log_lambda_dur <- numeric(N)
        xrow_dur <- xrow
        # duration.seasonal_term: "harmonic2" names these four columns;
        # load_duration_data() builds them the same way at fit time. Inert for
        # a fit without them: ctx$need_dur (the fit's own terms) selects the
        # columns below, so they never reach evaluate_lp().
        mnum_dur <- as.numeric(lubridate::month(months[m]))
        xrow_dur[["sin_month"]] <- rep(sin(2 * pi * mnum_dur / 12), N)
        xrow_dur[["cos_month"]] <- rep(cos(2 * pi * mnum_dur / 12), N)
        xrow_dur[["sin2_month"]] <- rep(sin(4 * pi * mnum_dur / 12), N)
        xrow_dur[["cos2_month"]] <- rep(cos(4 * pi * mnum_dur / 12), N)
        for (ctx in contexts) {
          idx <- ctx$idx
          xd <- lapply(
            xrow_dur[intersect(names(xrow_dur), ctx$need_dur)],
            function(v) v[idx]
          )
          # Duration cross-basis columns, exactly as occ/mon above. The
          # deployed per-group best-lag duration fits are `_cb` fits carrying
          # `enso_cb<k>` fixed effects; without this injection evaluate_lp()
          # silently skips them (absent covariate) and the track is ENSO-blind.
          # NULL/empty for a duration fit with no cross-basis.
          cbd <- cb_mats_dur_byctx[[ctx$gname]]
          for (cn in names(cbd)) {
            xd[[cn]] <- cbd[[cn]][idx, m]
          }
          bh_g <- basehaz_dur_eval(
            ctx$draws_dur,
            s_dur_byctx[[ctx$gname]],
            mso[idx],
            length(idx)
          )
          log_lambda_dur[idx] <- evaluate_lp(
            ctx$draws_dur,
            s_dur_byctx[[ctx$gname]],
            xd,
            o[idx],
            o_ct[idx],
            o_tu[idx],
            u_days[idx],
            ctx$frame_g,
            basehaz_val = bh_g,
            drop_region = ctx$drop_region,
            drop_ct_tu = drop_country_tu_spillover,
            plan = ctx$plan_dur,
            smooth_X = ctx$smooth_dur[[m]]
          ) +
            # offset(log(base_rate_dur)) -- denominator is exposure days, the
            # same log_offset_mat term added below. 0 for un-anchored fits.
            .anchor_log_offset(ctx$draws_dur, frame_tu, idx, "duration kernel")
        }
        # Same Poisson-with-log-offset link as occurrence (04_fit_dur_brms.R /
        # build_duration_formula_parts(): offset(log_offset) on days-in-month
        # exposure, never optional for this track -- hence no bundle gate on the
        # onset correction, unlike the mon block above). Duration shares the mon
        # track's episode-overlap exposure (both panels come from
        # .expand_active_months()), so it carries the same onset defect: without
        # the correction the cessation hazard for an episode's FIRST month is
        # applied to ~31 days instead of ~15, making episodes die too early.
        lambda_dur <- exp(
          log_lambda_dur + log_offset_mat[, m] + .onset_log_offset_adj(mso)
        )
        p_dur <- -expm1(-lambda_dur)

        act_idx <- which(active == 1L)
        ended <- integer(N)
        if (length(act_idx) > 0L) {
          ended[act_idx] <- as.integer(stats::rbinom(
            length(act_idx),
            1L,
            p_dur[act_idx]
          ))
        }

        end_idx <- which(ended == 1L)
        if (length(end_idx) > 0L) {
          rn_ok <- !is.na(open_row[end_idx])
          rn <- open_row[end_idx][rn_ok]
          ev_dur[rn] <- mso[end_idx][rn_ok]
          ev_end[rn] <- m
          ev_cens[rn] <- 0L
          active[end_idx] <- 0L
          mso[end_idx] <- 0L
          open_row[end_idx] <- NA_integer_
        }
      }

      # Monthly case-count draw under cessation_exposure: half -- after the
      # cessation draw, on the units active this month (act_m), with the
      # month's pre-update lags; a unit that ceased this month gets half its
      # exposure. Same quantities as the `full` call above, plus `ceased`.
      if (has_mon && cess_half && any(act_m == 1L)) {
        mc <- mon_case_step(
          m,
          act_m,
          mso_m,
          o_m,
          o_ct_m,
          o_tu_m,
          u_days_m,
          ceased = act_m == 1L & active == 0L
        )
        monthly_MU[r_idx, mc$act_i, m] <- mc$mu
        ok <- !is.na(mc$cm)
        if (any(ok)) monthly_C[r_idx, mc$act_i[ok], m] <- mc$cm[ok]
      }
    }

    if (has_dur) {
      # End-of-horizon administrative censoring: any location still active
      # when the simulated grid runs out gets its still-open event row
      # patched as censored at the last simulated month, mirroring the
      # fit-time panel's own administrative-censoring convention
      # (analysis/utils/dur_data.R).
      still_open <- which(active == 1L & !is.na(open_row))
      if (length(still_open) > 0L) {
        rn <- open_row[still_open]
        ev_dur[rn] <- mso[still_open]
        ev_end[rn] <- M
        ev_cens[rn] <- 1L
      }
    }

    if (length(ev_loc)) {
      event_rows[[r_idx]] <- if (has_dur) {
        data.frame(
          run_tag = NA_character_, # filled at write-time by 09_simulator
          replicate_id = rep_id,
          gid = frame_tu$gid[ev_loc],
          country = frame_tu$country[ev_loc],
          tu = frame_tu$tu[ev_loc],
          start_month = as.character(months[ev_m]),
          duration_months = ev_dur,
          end_month = as.character(months[ev_end]),
          censored = ev_cens,
          posterior_draw_o = s_occ_byrow[ev_loc],
          posterior_draw_d = s_dur_byrow[ev_loc],
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          run_tag = NA_character_, # filled at write-time by 09_simulator
          replicate_id = rep_id,
          gid = frame_tu$gid[ev_loc],
          country = frame_tu$country[ev_loc],
          tu = frame_tu$tu[ev_loc],
          start_month = as.character(months[ev_m]),
          duration_months = 1L, # legacy path: no duration model
          posterior_draw_o = s_occ_byrow[ev_loc],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  event_table <- dplyr::bind_rows(event_rows)

  # Monthly summary: collapse the arrays to long form
  rid <- rep(replicate_ids, each = N * M)
  gid_v <- rep(rep(frame_tu$gid, each = M), times = R)
  m_v <- rep(rep(as.character(months), times = N), times = R)
  # aperm reorders dim c(R, N, M) → c(M, N, R) so as.vector varies month fastest,
  # matching the rid/gid_v/m_v index construction (replicate slowest, month fastest).
  # Without this, as.vector varies R fastest (column-major), misaligning every
  # (replicate_id, gid, month) label with its Y value and flattening all seasonal curves.
  monthly_summary <- tibble::tibble(
    replicate_id = rid,
    gid = gid_v,
    month = m_v,
    Y = as.vector(aperm(monthly_Y, c(3L, 2L, 1L)))
  )
  # Optional channels appended in flag order, so the legacy (has_dur = FALSE)
  # and duration-only column sets are unchanged from the two-branch form this
  # replaced. Each uses the same aperm() reordering as the columns above.
  if (has_dur) {
    monthly_summary$active <- as.vector(aperm(monthly_A, c(3L, 2L, 1L)))
  }
  if (has_mon) {
    monthly_summary$mu_mon <- as.vector(aperm(monthly_MU, c(3L, 2L, 1L)))
    monthly_summary$cases_month <- as.vector(aperm(monthly_C, c(3L, 2L, 1L)))
  }
  if (use_det) {
    monthly_summary$detected <- as.vector(aperm(monthly_D, c(3L, 2L, 1L)))
  }

  if (!dry_run) {
    dir.create(output_part, recursive = TRUE, showWarnings = FALSE)
    ext <- sim_io_ext(out_fmt)
    sim_io_save(
      event_table,
      file.path(output_part, paste0("event_table", ext)),
      out_fmt
    )
    # Format 2 (utils/sim_io.R::sim_monthly_save()): the same table, stored
    # sparsely; every reader gets it back identical through sim_io_read().
    sim_monthly_save(
      monthly_summary,
      file.path(output_part, paste0("monthly_summary", ext)),
      out_fmt
    )
    assignment_obj <- list(
      tu = tu,
      replicate_ids = replicate_ids,
      posterior_idx_occ = posterior_idx_occ,
      seed = seed,
      N_locations = N,
      M_steps = M,
      generated_at = Sys.time()
    )
    if (has_dur) {
      assignment_obj$posterior_idx_dur <- posterior_idx_dur
    }
    if (has_mon) {
      assignment_obj$posterior_idx_mon <- posterior_idx_mon
    }
    qs_save(assignment_obj, file.path(output_part, "assignment.qs"))
    writeLines(format(Sys.time()), file.path(output_part, "_SUCCESS"))
  }

  invisible(list(event_table = event_table, monthly_summary = monthly_summary))
}
