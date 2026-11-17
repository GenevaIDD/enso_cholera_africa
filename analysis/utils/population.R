# population.R — read side of the LP x gid population crosswalk: loaders,
#   population resolution, composite apportionment, canonical admin levels
#
# Provides:
#   load_population_crosswalk(), load_gid_population(), load_wpp_table()
#   resolve_gid_population(), resolve_lp_population()   tiered population
#   rescale_population_wpp()       move a population between years (WPP)
#   build_lp_observation_operator(), expected_lp_from_gid_rates()
#   apportion_composite_cases()    split a composite LP's cases over its gids
#   reconcile_footprint_pop()      footprint population of a gid-season
#   canonical_admin_level_table(), canonical_level_lookup(),
#   apply_canonical_admin_level()  one source of truth for a gid's level
#
# Details:
#   Reads the crosswalk artefacts that population_core.R builds and
#   02c_preprocess_population_compile.R compiles.
#
# See: the "accessors" section below for the crosswalk algebra

# ===========================================================================
# LP x gid population crosswalk -- accessors
# ===========================================================================
#
# Globally-keyed sidecars that fit time and simulate time must read through
# the SAME function so the two sides cannot drift.
#
#   P_l = sum_g H[l,g]   population of a LOCATION PERIOD -- the denominator of
#                        an OBSERVATION (size fits)
#   Q_g                  population of a GADM UNIT -- the denominator of a
#                        FORECAST UNIT (simulator, smoothing)
#   c[g<-l] = H[l,g]/Q_g share of gid g covered by LP l -- the observation
#                        operator that maps simulated gid cases into LP space
#
# Composite detection is NOT re-implemented here: `gid_cover`, `n_gid_eff` and
# `is_composite` are materialised onto the crosswalk by
# derive_crosswalk_shares() at compile time, so consumers join rather than
# recompute. compute_gid_cover() lives in analysis/utils/population_core.R for the
# rare consumer that needs it from raw H.

if (!exists(".POP_XWALK_CACHE", envir = globalenv(), inherits = FALSE)) {
  .POP_XWALK_CACHE <- new.env(parent = emptyenv())
}

#' load_population_crosswalk
#'
#' Memoised read of `population_crosswalk_<Y0>-<Y1>.qs`. Memoised because
#' 03_prepare_siz.R and the simulator both hit it inside loops over years, and
#' the artefact is tens of MB.
#'
#' @param path .qs written by 02c_preprocess_population_compile.R.
#' @param year optional integer vector; filters after the cached read.
#' @return tibble with at least lctn_pr, gid, gadm_level, year, H, P_l, Q, w, c.
load_population_crosswalk <- function(path, year = NULL) {
  if (is.null(path) || !nzchar(path)) {
    stop(
      "load_population_crosswalk: no path given. Set population.crosswalk_path ",
      "or run analysis/02c_preprocess_population_compile.R."
    )
  }
  if (!file.exists(path)) {
    stop(
      "load_population_crosswalk: crosswalk not found: ",
      path,
      "\n  Build it with: Rscript analysis/02c_preprocess_population_compile.R -C <cfg>"
    )
  }
  key <- paste0("xw:", normalizePath(path, mustWork = FALSE))
  if (is.null(.POP_XWALK_CACHE[[key]])) {
    .POP_XWALK_CACHE[[key]] <- qs2::qs_read(path)
  }
  out <- .POP_XWALK_CACHE[[key]]
  if (!is.null(year)) {
    out <- out[out$year %in% as.integer(year), , drop = FALSE]
  }
  out
}

#' load_gid_population
#'
#' Memoised read of `gid_population_<Y0>-<Y1>.qs` (the canonical Q_g table).
load_gid_population <- function(path, year = NULL) {
  if (is.null(path) || !nzchar(path)) {
    stop(
      "load_gid_population: no path given. Set population.gid_pop_path ",
      "or run analysis/02c_preprocess_population_compile.R."
    )
  }
  if (!file.exists(path)) {
    stop(
      "load_gid_population: gid population table not found: ",
      path,
      "\n  Build it with: Rscript analysis/02c_preprocess_population_compile.R -C <cfg>"
    )
  }
  key <- paste0("gp:", normalizePath(path, mustWork = FALSE))
  if (is.null(.POP_XWALK_CACHE[[key]])) {
    .POP_XWALK_CACHE[[key]] <- qs2::qs_read(path)
  }
  out <- .POP_XWALK_CACHE[[key]]
  if (!is.null(year)) {
    out <- out[out$year %in% as.integer(year), , drop = FALSE]
  }
  out
}

#' load_wpp_table
#'
#' Read `wpp2024_national_totals.csv` (iso3, year, pop_total). Returns NULL when
#' absent, which downgrades T2 rescaling to unavailable rather than erroring --
#' the caller decides whether that is fatal.
load_wpp_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  key <- paste0("wpp:", normalizePath(path, mustWork = FALSE))
  if (is.null(.POP_XWALK_CACHE[[key]])) {
    .POP_XWALK_CACHE[[key]] <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  .POP_XWALK_CACHE[[key]]
}

#' rescale_population_wpp
#'
#' Move a population from `from_year` to `to_year` by the WPP2024 national
#' ratio, holding the sub-national shares frozen. This is the ONLY extrapolation
#' used outside WorldPop's 2015-2030 constrained coverage: it changes the level
#' and never the geography, so a gid's share of its country is preserved exactly.
#'
#' A MALFORMED table is not the same thing as an ABSENT one. `wpp = NULL` means
#' "T2 is unavailable" and returns NA, which the caller may legitimately absorb.
#' A table that is present but lacks `pop_total` is a broken cache file, and is
#' fatal here: `wpp$pop_total` on such a frame yields NULL, `ifelse()` then
#' returns `logical(0)`, and the length-0 result surfaces four frames up as
#' `argument is of length zero` inside resolve_gid_population()'s T2 loop. Name it
#' at the source instead.
#'
#' @return numeric, same length as `pop`; NA where the WPP lookup fails.
rescale_population_wpp <- function(pop, iso3, from_year, to_year, wpp) {
  if (is.null(wpp)) {
    return(rep(NA_real_, length(pop)))
  }
  miss <- setdiff(c("iso3", "year", "pop_total"), names(wpp))
  if (length(miss)) {
    stop(
      "rescale_population_wpp: the WPP table is missing required column(s) ",
      paste(miss, collapse = ", "),
      " (have: ",
      paste(names(wpp), collapse = ", "),
      "). Expected the ",
      "wpp2024_national_totals.csv schema (iso3, year, pop_total). Pass ",
      "wpp = NULL to disable T2 rescaling deliberately.",
      call. = FALSE
    )
  }
  n <- max(length(pop), length(iso3), length(from_year), length(to_year))
  pop <- rep_len(pop, n)
  iso3 <- rep_len(toupper(as.character(iso3)), n)
  from_year <- rep_len(as.integer(from_year), n)
  to_year <- rep_len(as.integer(to_year), n)

  k <- paste(wpp$iso3, wpp$year)
  a <- wpp$pop_total[match(paste(iso3, from_year), k)]
  b <- wpp$pop_total[match(paste(iso3, to_year), k)]
  ifelse(is.finite(a) & a > 0 & is.finite(b), pop * b / a, NA_real_)
}

#' resolve_gid_population
#'
#' The ONLY lookup path a consumer should use to obtain a canonical gid
#' population. Tiered exactly like resolve_mean_rate(), and for the same reason:
#' a silent default is not permitted (an unmatched gid must not, say, forecast
#' from a population of 1 with no warning).
#'
#'   T1 "direct"        exact (gid, year) in the table
#'   T2 "wpp_rescale"   same gid, nearest extracted year, rescaled by the
#'                      WPP2024 national ratio (the pre-2015 / post-2030 tail)
#'   T3 "parent_share"  gid absent but its parent ADM1 present: assign the
#'                      residual Q(parent) - sum(Q(present siblings)), split
#'                      evenly across the missing siblings. Only fires for a
#'                      partially-extracted ADM2 division.
#'   T4                 hard stop()
#'
#' @param gid    character vector of GADM ids.
#' @param year   integer vector, recycled to length(gid).
#' @param gid_pop table from load_gid_population().
#' @param wpp    table from load_wpp_table(), or NULL to disable T2.
#' @param on_missing "error" (default) or "na" -- "na" exists only for QC
#'   reporting, never for a model input.
#' @return tibble(gid, year, pop, pop_tier), in input order.
resolve_gid_population <- function(
  gid,
  year,
  gid_pop,
  wpp = NULL,
  on_missing = c("error", "na")
) {
  on_missing <- match.arg(on_missing)
  n <- length(gid)
  year <- rep_len(as.integer(year), n)
  gid <- as.character(gid)

  out <- tibble::tibble(
    gid = gid,
    year = year,
    pop = NA_real_,
    pop_tier = NA_character_
  )

  # --- T1: exact (gid, year) -------------------------------------------------
  k1 <- match(paste(gid, year), paste(gid_pop$gid, gid_pop$year))
  hit <- !is.na(k1)
  out$pop[hit] <- gid_pop$Q[k1[hit]]
  out$pop_tier[hit] <- "direct"

  # --- T2: same gid, nearest extracted year, WPP-rescaled --------------------
  need <- is.na(out$pop)
  if (any(need) && !is.null(wpp)) {
    # nearest extracted year per needed pair (first row on a tie), then one
    # vectorised rescale: a per-pair call re-pasted the whole WPP key each time
    # (2,860 pairs: 25 s -> 0.06 s, identical; tests/test_pop_censor_golden.R)
    ii <- which(need)
    rows_by <- split(seq_len(nrow(gid_pop)), gid_pop$gid)
    j <- vapply(
      ii,
      function(i) {
        rows <- rows_by[[match(gid[i], names(rows_by))]]
        if (is.null(rows)) {
          return(NA_integer_)
        }
        jj <- rows[which.min(abs(gid_pop$year[rows] - year[i]))]
        if (length(jj)) jj else NA_integer_
      },
      integer(1)
    )
    ok <- !is.na(j)
    resc <- rep(NA_real_, length(ii))
    resc[ok] <- rescale_population_wpp(
      gid_pop$Q[j[ok]],
      gid_pop$iso3[j[ok]],
      gid_pop$year[j[ok]],
      year[ii[ok]],
      wpp
    )
    fin <- is.finite(resc)
    out$pop[ii[fin]] <- resc[fin]
    out$pop_tier[ii[fin]] <- "wpp_rescale"
  }

  # --- T3: parent ADM1 residual ---------------------------------------------
  need <- is.na(out$pop)
  if (any(need)) {
    # Keys are the FLAT underscore form on both sides: 02c normalises the
    # canonical artefacts with flatten_gid() at the write seam, and the
    # outbreak-side gids arrive flattened from 03_prepare_*. So an ADM2 id is
    # "AGO_1_2_1" and its parent ADM1 is "AGO_1_1" (ISO3_ADM1n_suffix).
    # The "_" after the ISO3 is optional: GADM 4.1 omits it for Ghana, whose
    # ADM2 "GHA14_13_2" then has parent "GHA14_2" (the form gid_pop keys on).
    parent_of <- function(g) {
      sub("^([A-Z]{3}_?[0-9]+)_[0-9]+(_.*)$", "\\1\\2", g)
    }
    for (i in which(need)) {
      p <- parent_of(gid[i])
      if (identical(p, gid[i])) {
        next
      } # not an ADM2 id
      pj <- which(gid_pop$gid == p & gid_pop$year == year[i])
      if (!length(pj)) {
        next
      }
      # Sibling prefix = parent minus its trailing suffix segment, e.g.
      # "AGO_1_1" -> "AGO_1_", which matches the ADM2 children AGO_1_<k>_1
      # without bleeding into the neighbouring ADM1 "AGO_10_1" (whose next
      # character is "0", not "_"). Must strip only the LAST segment: the
      # dotted-era `sub("_.*$", "", p)` worked because a dotted ADM1 id had
      # exactly one underscore, but a flattened one has several, so cutting at
      # the first would widen the prefix to "AGO_" and sum every ADM2 in the
      # country into the residual.
      sib <- gid_pop[
        gid_pop$gadm_level == "ADM2" &
          gid_pop$year == year[i] &
          startsWith(gid_pop$gid, paste0(sub("_[^_]*$", "", p), "_")),
      ]
      residual <- gid_pop$Q[pj[1]] - sum(sib$Q, na.rm = TRUE)
      if (is.finite(residual) && residual > 0) {
        out$pop[i] <- residual
        out$pop_tier[i] <- "parent_share"
      }
    }
  }

  # --- T4: refuse to invent a denominator ------------------------------------
  bad <- is.na(out$pop)
  if (any(bad) && identical(on_missing, "error")) {
    stop(
      "resolve_gid_population: no population resolves for ",
      sum(bad),
      " (gid, year) pair(s): ",
      paste(
        utils::head(unique(paste(out$gid[bad], out$year[bad], sep = "@")), 8),
        collapse = ", "
      ),
      "\n  Refusing to default (cf. the pop = 1.0 fallback this replaces). ",
      "Either extend the extraction years (population.years, 02a) or ",
      "give the WPP table (simulator.location_frame.wpp_path) for the T2 ",
      "rescale of years outside them."
    )
  }
  if (any(bad)) {
    out$pop_tier[bad] <- "missing"
  }
  out
}

#' resolve_lp_population
#'
#' P_l for a location period. Same tier ladder as resolve_gid_population() with
#' one deliberate difference: it NEVER returns NA. `pop` enters the size model
#' as `offset(log(pop))`, so a single NA is not a missing covariate -- it drops
#' the row, or produces a non-finite offset that test_offset_wellformed.R
#' rejects for the whole corpus. The upstream `pop` column is therefore a
#' guaranteed last resort, tagged so the fallback count is a reported QC number
#' rather than an invisible substitution.
#'
#' @param lctn_pr      character vector of location-period ids.
#' @param year         integer vector, recycled.
#' @param xwalk        table from load_population_crosswalk().
#' @param pop_upstream numeric vector of the legacy upstream `pop`, recycled.
#' @param gadm_level   which tiling to read P_l from ("ADM1" or "ADM2"); the two
#'                     agree to within the sliver tolerance by construction.
#' @param wpp          table from load_wpp_table(), or NULL to disable T2.
#' @return tibble(lctn_pr, year, pop, pop_tier), in input order.
resolve_lp_population <- function(
  lctn_pr,
  year,
  xwalk,
  pop_upstream = NA_real_,
  gadm_level = "ADM1",
  wpp = NULL
) {
  n <- length(lctn_pr)
  lctn_pr <- as.character(lctn_pr)
  year <- rep_len(as.integer(year), n)
  pop_upstream <- rep_len(as.numeric(pop_upstream), n)

  lvl <- xwalk[xwalk$gadm_level == gadm_level, , drop = FALSE]
  if (!nrow(lvl)) {
    lvl <- xwalk
  } # a corpus extracted at one level only

  p <- lvl %>%
    dplyr::group_by(lctn_pr, year) %>%
    dplyr::summarise(
      P_l = sum(H, na.rm = TRUE),
      iso3 = dplyr::first(gid_iso3),
      .groups = "drop"
    )

  out <- tibble::tibble(
    lctn_pr = lctn_pr,
    year = year,
    pop = NA_real_,
    pop_tier = NA_character_
  )

  k1 <- match(paste(lctn_pr, year), paste(p$lctn_pr, p$year))
  hit <- !is.na(k1) & is.finite(p$P_l[k1]) & (p$P_l[k1] > 0)
  out$pop[hit] <- p$P_l[k1[hit]]
  out$pop_tier[hit] <- "direct"

  need <- is.na(out$pop)
  if (any(need) && !is.null(wpp)) {
    for (i in which(need)) {
      rows <- which(p$lctn_pr == lctn_pr[i] & is.finite(p$P_l) & p$P_l > 0)
      if (!length(rows)) {
        next
      }
      j <- rows[which.min(abs(p$year[rows] - year[i]))]
      resc <- rescale_population_wpp(
        p$P_l[j],
        p$iso3[j],
        p$year[j],
        year[i],
        wpp
      )
      if (is.finite(resc) && resc > 0) {
        out$pop[i] <- resc
        out$pop_tier[i] <- "wpp_rescale"
      }
    }
  }

  fb <- is.na(out$pop)
  out$pop[fb] <- pop_upstream[fb]
  out$pop_tier[fb] <- "upstream_fallback"

  # The contract above: no NA leaves this function.
  still <- !is.finite(out$pop) | out$pop <= 0
  if (any(still)) {
    stop(
      "resolve_lp_population: ",
      sum(still),
      " location period(s) have neither ",
      "a crosswalk P_l nor a usable upstream pop, e.g. ",
      paste(utils::head(unique(out$lctn_pr[still]), 5), collapse = ", "),
      ". offset(log(pop)) would be non-finite."
    )
  }
  out
}

#' build_lp_observation_operator
#'
#' The sparse matrix c[g<-l] = H[l,g] / Q_g, as a long tibble. Mapping simulated
#' gid cases into LP space is a coverage-weighted sum, NOT a share-weighted one:
#' the simulator emits absolute cases Q_g * r_g per gid, so
#'
#'   E[y_l] = sum_g H[l,g] * r_g = sum_g c[g<-l] * sim_cases_g
#'
#' which stays correct when an LP covers only PART of a gid (c < 1).
#'
#' @param xwalk table from load_population_crosswalk().
#' @param year  integer; the population year to take coverage from.
#' @param gadm_level tiling to build the operator on.
#' @return tibble(lctn_pr, gid, c)
build_lp_observation_operator <- function(xwalk, year, gadm_level = "ADM1") {
  op <- xwalk %>%
    dplyr::filter(
      .data$year == as.integer(!!year),
      .data$gadm_level == !!gadm_level,
      !is.na(.data$c),
      .data$c > 0
    ) %>%
    dplyr::select(lctn_pr, gid, c)
  if (!nrow(op)) {
    stop(
      "build_lp_observation_operator: no coverage rows for year ",
      year,
      " at ",
      gadm_level,
      ". Was the crosswalk built for that year?"
    )
  }
  op
}

#' expected_lp_from_gid_rates
#'
#' Apply the operator: aggregate simulated per-gid cases into LP space.
#'
#' @param xwalk        table from load_population_crosswalk().
#' @param cases_by_gid tibble(gid, cases[, ...grouping cols]).
#' @param year         population year for the coverage weights.
#' @param by           optional extra grouping columns present in cases_by_gid
#'                     (e.g. c("rep_id", "month")) carried through the sum.
#' @return tibble(lctn_pr, <by>, expected_cases)
expected_lp_from_gid_rates <- function(
  xwalk,
  cases_by_gid,
  year,
  gadm_level = "ADM1",
  by = character(0)
) {
  op <- build_lp_observation_operator(
    xwalk,
    year = year,
    gadm_level = gadm_level
  )
  cases_by_gid %>%
    dplyr::mutate(gid = as.character(gid)) %>%
    dplyr::inner_join(op, by = "gid", relationship = "many-to-many") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("lctn_pr", by)))) %>%
    dplyr::summarise(
      expected_cases = sum(.data$c * .data$cases, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---------------------------------------------------------------------------
# Composite apportionment and footprint population reconciliation
#
# These two functions exist to delete two silent defects that live between the
# LP-keyed epi corpus and the gid-keyed model frame:
#
#   * the CENTROID COLLAPSE (analysis/02_preprocess_climate.R:72-106) reduces every
#     LP to one st_centroid(), so a composite LP -- a union of K admin units --
#     dumps all K units' cases AND all K units' population onto one arbitrary
#     gid, while the other K-1 are recorded as never having had an outbreak;
#   * `pop = sum(pop)` (analysis/03_prepare_siz.R:190) aggregates a STOCK with a
#     flow rule, so two location periods reaching the same (gid, season) bucket
#     have their populations ADDED, and an empty group makes sum() return 0.
#
# Both are pure functions of the crosswalk built by 02c_preprocess_population_compile.R, so
# both are unit-testable without a raster, a database, or data/ access -- see
# analysis/tests/test_case_apportionment.R and
# analysis/tests/test_footprint_pop_reconciliation.R.
# ---------------------------------------------------------------------------

#' .nearest_crosswalk_year
#'
#' The crosswalk is materialised only inside WorldPop's constrained coverage
#' (POP_RASTER_YEAR_MIN..MAX, currently 2015-2030) while the epi corpus starts
#' in 2010, so a row's season year frequently has no crosswalk year of its own.
#' `w` is a SHARE and is near year-invariant (the geometry does not move; only
#' the national total does), so snapping to the nearest extracted year is
#' correct for apportionment in a way it would NOT be for an absolute count --
#' that is what resolve_gid_population()'s T2 WPP rescale is for.
#'
#' @param want   integer vector of desired years.
#' @param have   integer vector of available years (crosswalk years).
#' @return integer vector, same length as `want`; NA where `have` is empty.
.nearest_crosswalk_year <- function(want, have) {
  have <- sort(unique(as.integer(have[is.finite(have)])))
  if (!length(have)) {
    return(rep(NA_integer_, length(want)))
  }
  vapply(
    as.integer(want),
    function(y) {
      if (!is.finite(y)) {
        return(have[1L])
      }
      have[which.min(abs(have - y))]
    },
    integer(1)
  )
}

#' apportion_composite_cases
#'
#' Expand each location-period row into one row per member GADM gid, splitting
#' the case count by the crosswalk population share w[l->g] = H[l,g] / P_l and
#' giving each member row `pop = H[l,g]`.
#'
#' WHY BOTH SCALE BY THE SAME w. Cases and population are divided by the same
#' factor, so the implied attack rate `cases / pop` is INVARIANT per member and
#' identical to the LP's. read_size_data()'s max_attack_rate guard therefore sees
#' exactly what it saw before, and the apportionment cannot manufacture an
#' implausible incidence. `min_pop` is different: it is an absolute floor, so
#' small-w members legitimately fall below it -- callers should count those.
#'
#' WHAT THIS DOES NOT CLAIM. Splitting one aggregate observation into K rows
#' fabricates K pseudo-observations that a likelihood will treat as K
#' independent data points, overstating precision. That is the known price of
#' the S-ii route (the alternative, S-iv, relocates the same assumption to the
#' smoothing handoff). The `apportioned` / `n_apportion_members` columns are the
#' hook for the mitigation: a sensitivity run can NA the response on those rows,
#' or widen the measurement dispersion in proportion to the member count.
#'
#' @param dat        tibble; one row per (location period, event).
#' @param xwalk      crosswalk from load_population_crosswalk().
#' @param lp_col     column of `dat` holding the location-period id. Note the
#'                   crosswalk calls it `lctn_pr` while the size/occurrence
#'                   prep scripts call it `location_period_id`.
#' @param cases_col  column to apportion, or NA_character_ to apportion nothing
#'                   (the occurrence track: a 0/1 indicator must be PROPAGATED,
#'                   not scaled -- see 03_prepare_occ.R).
#' @param year_col   column holding the season year, used to pick the crosswalk
#'                   year. Ignored if absent.
#' @param gid_col    column to overwrite with the member gid. The incoming value
#'                   (the centroid match) is preserved as `gid_centroid`.
#' @param level_col  optional column of `dat` giving the row's admin level; any
#'                   value other than "ADM1" selects the ADM2 tiling, mirroring
#'                   `admin_level_simple` at 02_preprocess_climate.R:77-78. When
#'                   NULL or absent, `gadm_level` is used for every row.
#' @param gadm_level scalar tiling used when `level_col` is unavailable.
#' @param keep_unmatched if TRUE (default) LPs absent from the crosswalk pass
#'                   through unchanged, tagged pop_tier = "unmatched_lp"; if
#'                   FALSE they are dropped. Never silently either way -- the
#'                   count is returned in attr(, "apportion_report").
#' @return `dat` with one row per (original row x member gid), carrying
#'   `gid` (member), `gid_centroid`, `pop`, `pop_lp`, `apportion_w`,
#'   `n_apportion_members`, `apportioned`, `pop_tier`. A summary tibble is
#'   attached as attr(, "apportion_report").
apportion_composite_cases <- function(
  dat,
  xwalk,
  lp_col = "location_period_id",
  cases_col = "suspected_cases",
  year_col = "year_seas",
  gid_col = "gid",
  level_col = "admin_level",
  gadm_level = "ADM1",
  tiling = c("derived", "corpus_label"),
  keep_unmatched = TRUE
) {
  stopifnot(is.data.frame(dat), is.data.frame(xwalk))
  # Validated here too, not only inside .match_lp_members(), because the empty-`dat`
  # short-circuit below must still reject a malformed crosswalk (as it always has).
  need <- c("lctn_pr", "gid", "gadm_level", "year", "H")
  miss <- setdiff(need, names(xwalk))
  if (length(miss)) {
    stop(
      "apportion_composite_cases: crosswalk is missing column(s): ",
      paste(miss, collapse = ", "),
      ". Was it built by 02c_preprocess_population_compile.R?"
    )
  }
  if (!lp_col %in% names(dat)) {
    stop("apportion_composite_cases: `dat` has no column '", lp_col, "'.")
  }
  if (!nrow(dat)) {
    return(dat)
  }

  mm <- .match_lp_members(
    dat,
    xwalk,
    lp_col = lp_col,
    year_col = year_col,
    level_col = level_col,
    gadm_level = gadm_level,
    tiling = tiling,
    keep_unmatched = keep_unmatched,
    .fn = "apportion_composite_cases"
  )
  matched <- mm$matched
  miss_rows <- mm$miss_rows

  base_cols <- function(d) {
    d$gid_centroid <- if (gid_col %in% names(d)) {
      as.character(d[[gid_col]])
    } else {
      NA_character_
    }
    d$pop_lp <- if ("pop" %in% names(d)) as.numeric(d$pop) else NA_real_
    d
  }

  out_hit <- dat[matched$.row, , drop = FALSE] %>% base_cols()
  out_hit[[gid_col]] <- matched$member_gid
  out_hit$pop <- matched$H
  # base_cols() captured pop_lp from `dat$pop` BEFORE `pop` above was
  # overwritten with the member share H -- i.e. from whatever pop column
  # (if any) already existed on the input, such as a legacy population
  # estimate baked into outbreak_metadata.rds. That is not what `pop_lp` is
  # documented to mean (the whole-LP crosswalk population, P_l). Overwrite it
  # here with the real P_l from the match, now that it is available.
  out_hit$pop_lp <- matched$P_l
  out_hit$apportion_w <- matched$apportion_w
  out_hit$n_apportion_members <- as.integer(matched$n_apportion_members)
  out_hit$apportioned <- matched$n_apportion_members > 1L
  out_hit$pop_tier <- "crosswalk_H"
  out_hit$pop_year <- matched$year
  if (
    !is.null(cases_col) && !is.na(cases_col) && cases_col %in% names(out_hit)
  ) {
    out_hit[[cases_col]] <- as.numeric(out_hit[[cases_col]]) *
      matched$apportion_w
  }

  out_miss <- dat[miss_rows, , drop = FALSE] %>% base_cols()
  if (nrow(out_miss)) {
    out_miss$apportion_w <- 1
    out_miss$n_apportion_members <- 1L
    out_miss$apportioned <- FALSE
    out_miss$pop_tier <- "unmatched_lp"
    out_miss$pop_year <- NA_integer_
    if (!"pop" %in% names(out_miss)) out_miss$pop <- NA_real_
  }

  out <- if (keep_unmatched) dplyr::bind_rows(out_hit, out_miss) else out_hit

  attr(out, "apportion_report") <- tibble::tibble(
    n_rows_in = nrow(dat),
    n_rows_out = nrow(out),
    n_lp_matched = dplyr::n_distinct(matched$lctn_pr),
    n_lp_unmatched = dplyr::n_distinct(as.character(dat[[lp_col]])[miss_rows]),
    n_rows_composite = sum(out$apportioned, na.rm = TRUE),
    max_members = suppressWarnings(max(
      c(0L, out$n_apportion_members),
      na.rm = TRUE
    ))
  )
  out
}

#' reconcile_footprint_pop
#'
#' The denominator for a (gid, season, admin_level) bucket after the footprint
#' collapse. Replaces `pop = sum(pop, na.rm = TRUE)`
#' (analysis/03_prepare_siz.R:190), which is a flow rule applied to a stock.
#'
#'   pop = min( sum_l H[l,g] , Q_g )
#'
#' Two location periods that tile the gid (ADM3 children relabelled to their
#' ADM2 parent) sum correctly, because their H pieces are disjoint. Two location
#' periods that BOTH report the whole gid saturate at Q_g instead of doubling
#' it. The cap is the entire point: no set of location periods can contain more
#' people than the gid does.
#'
#' The `H_sum <= 0` branch is the COD_4_6_1 defect: `sum()` of an empty group is
#' 0, and 0 then travelled downstream as a population. Zero is never a
#' population -- fall back to the canonical Q_g.
#'
#' @param gid      character vector of GADM ids (one per bucket).
#' @param year     integer vector, recycled.
#' @param H_sum    numeric vector; sum of the crosswalk mass in the bucket.
#' @param gid_pop  table from load_gid_population().
#' @param wpp      table from load_wpp_table(), or NULL to disable the T2 rescale.
#' @param on_missing "error" (default) or "na", for buckets where neither a
#'                 positive H_sum nor a resolvable Q_g exists.
#' @return tibble(gid, year, pop, pop_tier) in input order. pop_tier is one of
#'   "footprint" (H_sum below Q, the normal case), "gid_capped" (H_sum exceeded
#'   Q -- overlapping reports), "gid_fallback" (no LP mass; Q used), or
#'   "lp_mass_only" (no canonical Q; H_sum used).
reconcile_footprint_pop <- function(
  gid,
  year,
  H_sum,
  gid_pop,
  wpp = NULL,
  on_missing = c("error", "na")
) {
  on_missing <- match.arg(on_missing)
  n <- length(gid)
  gid <- as.character(gid)
  year <- rep_len(as.integer(year), n)
  H_sum <- rep_len(as.numeric(H_sum), n)

  Q <- resolve_gid_population(
    gid,
    year,
    gid_pop,
    wpp = wpp,
    on_missing = "na"
  )$pop

  has_H <- is.finite(H_sum) & H_sum > 0
  has_Q <- is.finite(Q) & Q > 0

  pop <- rep(NA_real_, n)
  tier <- rep(NA_character_, n)

  both <- has_H & has_Q
  pop[both] <- pmin(H_sum[both], Q[both])
  # A tiny relative slack so floating-point equality reports "footprint" rather
  # than a spurious cap.
  tier[both] <- ifelse(
    H_sum[both] > Q[both] * (1 + 1e-9),
    "gid_capped",
    "footprint"
  )

  only_h <- has_H & !has_Q
  pop[only_h] <- H_sum[only_h]
  tier[only_h] <- "lp_mass_only"

  only_q <- !has_H & has_Q
  pop[only_q] <- Q[only_q]
  tier[only_q] <- "gid_fallback"

  bad <- !is.finite(pop) | pop <= 0
  if (any(bad) && identical(on_missing, "error")) {
    stop(
      "reconcile_footprint_pop: ",
      sum(bad),
      " (gid, year) bucket(s) have ",
      "neither crosswalk mass nor a canonical Q_g, e.g. ",
      paste(
        utils::head(unique(paste(gid[bad], year[bad], sep = "@")), 5),
        collapse = ", "
      ),
      ".\n  Refusing to return 0 or NA as a population (cf. the ",
      "`pop = sum(pop)` empty-group zero this replaces). Either extend the ",
      "population extraction years or drop these units explicitly."
    )
  }
  if (any(bad)) {
    pop[bad] <- NA_real_
    tier[bad] <- "missing"
  }

  tibble::tibble(gid = gid, year = year, pop = pop, pop_tier = tier)
}

# ===========================================================================
# One source of truth for a gid's admin level
# ===========================================================================
#
# `admin_level` arrives from the outbreak corpus as a LABEL, and it is wrong
# often enough to matter: measured against the canonical `gadm_level` in
# `gid_population_<Y0>-<Y1>.qs`, 157 of 616 simulator-frame gids (25.5%) and
# 468 of 1582 size rows disagree. That matters because `admin_level` is not
# cosmetic -- it is the POOLING KEY for `build_mean_rate_table()` and the
# stratifying key for the population elasticity, so a mis-keyed row is pooled
# with the wrong cell and gets the wrong elasticity.
#
# The sharpest symptom: the canonical table contains NO ADM3 UNITS AT ALL,
# while the frame labels 77 gids "ADM3" (canonically 7 ADM1 + 70 ADM2). Every
# ADM3 quantity downstream -- the ADM3 `mean_rate` cells, the ADM3
# `level_factor`, the "at ADM3 size is independent of population" elasticity
# reading -- is therefore computed on a PHANTOM level. Re-keying moves the
# aggregate simulated mass by a factor of 0.708.
#
# WHY NOT parse the gid string. The obvious fix -- count separators, since
# "AGO_10_1" looks like it encodes depth -- DOES NOT WORK. Gids carry a
# trailing GADM version component, so `ADM(depth - 1)` assigns ADM1 to nothing
# and mislabels everything else.
#
# WHY NOT fix it in place at the source. The corpus label is what the outbreak
# database asserts, and it is also what `filter_overlapping_adm1()` and
# `drop_overlapping_adm1_from_frame()` currently key off. Overwriting it silently
# would change which rows those filters drop, with no record. So this function
# writes a NEW column, keeps the original as `admin_level_corpus`, and reports
# the disagreement count -- the correction is opt-in and auditable.

#' canonical_admin_level_table
#'
#' The gid -> canonical `gadm_level` map, one row per gid. `gid_population` is
#' a (gid x year) panel, so it is de-duplicated here; a gid whose level varies
#' across years would be a corrupt table and is reported rather than silently
#' resolved by taking the first row.
#'
#' @param gid_pop the canonical table (from `load_gid_population()`), or a path.
#' @return tibble(gid, admin_level_canonical)
canonical_admin_level_table <- function(gid_pop) {
  if (is.character(gid_pop)) {
    gid_pop <- load_gid_population(gid_pop)
  }
  need <- c("gid", "gadm_level")
  miss <- setdiff(need, names(gid_pop))
  if (length(miss)) {
    stop(
      "canonical_admin_level_table: the canonical population table lacks ",
      paste(miss, collapse = ", "),
      " (have: ",
      paste(names(gid_pop), collapse = ", "),
      "). Expected the ",
      "gid_population_<Y0>-<Y1>.qs schema.",
      call. = FALSE
    )
  }

  tab <- tibble::tibble(
    gid = flatten_gid(as.character(gid_pop$gid)),
    gadm_level = as.character(gid_pop$gadm_level)
  )
  tab <- unique(tab[!is.na(tab$gid) & !is.na(tab$gadm_level), , drop = FALSE])

  dup <- tab$gid[duplicated(tab$gid)]
  if (length(dup)) {
    stop(
      "canonical_admin_level_table: ",
      length(unique(dup)),
      " gid(s) carry ",
      "MORE THAN ONE gadm_level in the canonical table, e.g. ",
      paste(utils::head(unique(dup), 5), collapse = ", "),
      ".\n  A gid's admin level cannot depend on the year; this table is ",
      "corrupt. Refusing to pick one arbitrarily.",
      call. = FALSE
    )
  }

  names(tab) <- c("gid", "admin_level_canonical")
  tab
}

#' canonical_level_lookup
#'
#' A `level_of(gid)` closure over the canonical table, for
#' footprint.R::build_canonical_episodes() (which stays dependency-free and so
#' takes the lookup as an argument). Unmatched gids are an ERROR: under
#' `--episode_key gid_season` the canonical level IS the episode's
#' admin_level, and there is no defensible default.
canonical_level_lookup <- function(gid_pop) {
  tab <- canonical_admin_level_table(gid_pop)
  function(gid) {
    key <- flatten_gid(as.character(gid))
    lv <- tab$admin_level_canonical[match(key, tab$gid)]
    if (anyNA(lv)) {
      stop(
        "canonical_level_lookup: ",
        length(unique(key[is.na(lv)])),
        " gid(s) ",
        "are absent from the canonical population table, e.g. ",
        paste(utils::head(unique(key[is.na(lv)]), 5), collapse = ", "),
        ". Refusing to assign an admin level by default.",
        call. = FALSE
      )
    }
    lv
  }
}

#' apply_canonical_admin_level
#'
#' Attach the canonical admin level to `dat` and (when `overwrite = TRUE`)
#' promote it to `admin_level`, preserving the corpus label as
#' `admin_level_corpus`.
#'
#' UNMATCHED GIDS ARE NOT SILENTLY RE-LABELLED. A gid absent from the canonical
#' table keeps its corpus label and is counted in the report. Defaulting an
#' unmatched row to any level would be exactly the "silent default" the
#' provenance gates exist to forbid; `on_unmatched = "error"` is available for callers
#' that require full coverage.
#'
#' @param dat data frame with a `gid` column (dotted or underscore form).
#' @param gid_pop canonical population table, or a path to it.
#' @param overwrite promote the canonical level to `admin_level`. FALSE only
#'   attaches `admin_level_canonical`, for measurement without behaviour change.
#' @param on_unmatched "keep" (default) | "error".
#' @param tag log prefix.
#' @return `dat` plus `admin_level_canonical`, `admin_level_corpus`, and (when
#'   `overwrite`) a rewritten `admin_level`. Carries an
#'   `attr(., "admin_level_recode")` report.
apply_canonical_admin_level <- function(
  dat,
  gid_pop,
  overwrite = TRUE,
  on_unmatched = c("keep", "error"),
  tag = "P3"
) {
  on_unmatched <- match.arg(on_unmatched)
  if (!"gid" %in% names(dat)) {
    stop(
      "apply_canonical_admin_level: `dat` has no `gid` column.",
      call. = FALSE
    )
  }

  tab <- canonical_admin_level_table(gid_pop)
  key <- flatten_gid(as.character(dat$gid))
  can <- tab$admin_level_canonical[match(key, tab$gid)]

  corpus <- if ("admin_level" %in% names(dat)) {
    as.character(dat$admin_level)
  } else {
    rep(NA_character_, nrow(dat))
  }

  unmatched <- is.na(can)
  if (any(unmatched) && identical(on_unmatched, "error")) {
    stop(
      "apply_canonical_admin_level: ",
      sum(unmatched),
      " row(s) / ",
      length(unique(key[unmatched])),
      " gid(s) are absent from the ",
      "canonical population table, e.g. ",
      paste(utils::head(unique(key[unmatched]), 5), collapse = ", "),
      ".\n  Refusing to assign an admin level by default.",
      call. = FALSE
    )
  }

  changed <- !unmatched & !is.na(corpus) & can != corpus
  rep <- list(
    n_rows = nrow(dat),
    n_unmatched = sum(unmatched),
    n_changed = sum(changed),
    gids_unmatched = length(unique(key[unmatched])),
    gids_changed = length(unique(key[changed])),
    overwrite = isTRUE(overwrite),
    crosstab = if (any(changed)) {
      ct <- as.data.frame(
        table(corpus = corpus[changed], canonical = can[changed]),
        stringsAsFactors = FALSE
      )
      ct[ct$Freq > 0L, , drop = FALSE] # table() emits empty cells
    } else {
      NULL
    }
  )

  dat$admin_level_corpus <- corpus
  dat$admin_level_canonical <- can
  if (isTRUE(overwrite)) {
    dat$admin_level <- ifelse(unmatched, corpus, can)
  }

  message(sprintf(
    paste0(
      "[%s] canonical admin_level: %d/%d row(s) re-keyed (%d gid(s)); ",
      "%d row(s) unmatched (label kept); overwrite=%s"
    ),
    tag,
    rep$n_changed,
    rep$n_rows,
    rep$gids_changed,
    rep$n_unmatched,
    isTRUE(overwrite)
  ))
  if (!is.null(rep$crosstab)) {
    message(paste0(
      "[",
      tag,
      "]   ",
      paste(
        sprintf(
          "%s->%s: %d",
          rep$crosstab$corpus,
          rep$crosstab$canonical,
          rep$crosstab$Freq
        ),
        collapse = "  "
      )
    ))
  }

  attr(dat, "admin_level_recode") <- rep
  dat
}
