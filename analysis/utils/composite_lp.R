
#' .matchLpMembers
#'
#' The member-matching core shared by `apportionCompositeCases()` (which expands
#' one LP row into K member rows) and `collapseCompositeLPs()` (which reduces one
#' LP row to its dominant member). Extracted verbatim from
#' `apportionCompositeCases()` so the two policies cannot drift apart in how they
#' pick a tiling, snap a year, or renormalise the population shares.
#'
#' Does three things and nothing else:
#'   1. per-row tiling (ADM1 vs ADM2) — from the crosswalk's geometry-derived
#'      `lp_level` under `tiling = "derived"`, or from the corpus's own
#'      `admin_level` string under `"corpus_label"` — with a fallback to
#'      whatever tiling the crosswalk actually carries;
#'   2. nearest-year snapping per (LP, tiling) via `.nearestCrosswalkYear()`;
#'   3. the join to the member pieces, renormalising `apportion_w = H / sum(H)`
#'      over the SURVIVING pieces (the sliver filter in population_core.R can
#'      drop pieces after `w` was derived, and conservation depends on
#'      `sum(w) == 1` exactly).
#'
#' It makes no decision about what to do with the match -- no scaling, no
#' dropping, no column rewriting. That is the caller's job.
#'
#' @param dat        tibble; one row per (location period, event).
#' @param xwalk      crosswalk from loadPopulationCrosswalk().
#' @param lp_col     column of `dat` holding the location-period id.
#' @param year_col   column holding the season year, used to pick the crosswalk
#'                   year. Ignored if absent.
#' @param level_col  optional column of `dat` giving the row's admin level; any
#'                   value other than "ADM1" selects the ADM2 tiling. Only
#'                   consulted under `tiling = "corpus_label"`, or for LPs the
#'                   crosswalk carries no `lp_level` for.
#' @param gadm_level scalar tiling used when `level_col` is unavailable.
#' @param tiling     "derived" (default) reads each LP at the crosswalk's frozen
#'                   `lp_level`; "corpus_label" trusts the corpus's
#'                   `admin_level` string.
#'                   Ignored (falls back to `level_col`) when the crosswalk has
#'                   no `lp_level` column.
#' @param keep_unmatched recorded in the report only; the caller decides.
#' @param .fn        label used in error messages, so the caller owns them.
#' @return list(matched, miss_rows, report):
#'   `matched` — one row per (input row x member gid), with `.row` (index into
#'     `dat`), `lctn_pr`, `gadm_level`, `year`, `member_gid`, `H`, `apportion_w`,
#'     `n_apportion_members`, `P_l`;
#'   `miss_rows` — integer indices of input rows with no crosswalk match;
#'   `report` — a one-row tibble of counts.
.matchLpMembers <- function(dat, xwalk,
                            lp_col     = "location_period_id",
                            year_col   = "year_seas",
                            level_col  = "admin_level",
                            gadm_level = "ADM1",
                            tiling     = c("derived", "corpus_label"),
                            keep_unmatched = TRUE,
                            .fn = ".matchLpMembers") {
  stopifnot(is.data.frame(dat), is.data.frame(xwalk))
  need <- c("lctn_pr", "gid", "gadm_level", "year", "H")
  miss <- setdiff(need, names(xwalk))
  if (length(miss))
    stop(.fn, ": crosswalk is missing column(s): ",
         paste(miss, collapse = ", "),
         ". Was it built by 02c_preprocess_population_compile.R?")
  if (!lp_col %in% names(dat))
    stop(.fn, ": `dat` has no column '", lp_col, "'.")

  # --- per-row tiling -------------------------------------------------------
  # `corpus_label` is kept for comparison and is WRONG in general: the
  # outbreak corpus's admin depth is not GADM's. Malawi's corpus "ADM2" is a
  # district, which GADM calls ADM1 (GADM's ADM2 is the Traditional Authority),
  # so reading it at ADM2 shattered one district into 8-25 units and `gid_cover`
  # measured within-district population concentration instead of compositeness.
  # `derived` uses the geometry-based `lp_level` frozen onto the crosswalk by
  # 02c (deriveLpTilingLevel()). Kept side-by-side, not deleted, so the earlier
  # results remain reproducible.
  tiling <- match.arg(tiling, c("derived", "corpus_label"))
  use_level <- !is.null(level_col) && level_col %in% names(dat)
  row_level <- if (use_level) {
    ifelse(as.character(dat[[level_col]]) == "ADM1", "ADM1", "ADM2")
  } else {
    rep(as.character(gadm_level), nrow(dat))
  }
  if (identical(tiling, "derived") && "lp_level" %in% names(xwalk)) {
    lp_lvl <- xwalk %>%
      dplyr::distinct(lctn_pr = as.character(.data$lctn_pr),
                      lp_level = as.character(.data$lp_level)) %>%
      dplyr::filter(!is.na(.data$lp_level))
    if (anyDuplicated(lp_lvl$lctn_pr))
      stop(.fn, ": the crosswalk carries more than one `lp_level` for some ",
           "location period(s) — it must be one decision per LP.")
    hit <- lp_lvl$lp_level[match(as.character(dat[[lp_col]]), lp_lvl$lctn_pr)]
    # LPs the crosswalk has no level for keep the corpus-label guess: there is
    # nothing better to fall back to, and it is what they got before.
    row_level <- ifelse(is.na(hit), row_level, hit)
  }
  # Fall back to the tiling that actually exists in the crosswalk, so a corpus
  # extracted at one level only does not silently apportion nothing.
  have_lvl <- unique(as.character(xwalk$gadm_level))
  row_level[!row_level %in% have_lvl] <- have_lvl[1L]

  row_year <- if (!is.null(year_col) && year_col %in% names(dat)) {
    suppressWarnings(as.integer(dat[[year_col]]))
  } else {
    rep(NA_integer_, nrow(dat))
  }

  # --- snap each row to the nearest crosswalk year for its (LP, level) -------
  xw <- xwalk %>%
    dplyr::mutate(lctn_pr    = as.character(.data$lctn_pr),
                  gid        = as.character(.data$gid),
                  gadm_level = as.character(.data$gadm_level),
                  year       = as.integer(.data$year)) %>%
    dplyr::filter(is.finite(.data$H), .data$H > 0)

  key <- tibble::tibble(.row = seq_len(nrow(dat)),
                        lctn_pr = as.character(dat[[lp_col]]),
                        gadm_level = row_level,
                        .want_year = row_year)
  yrs <- xw %>%
    dplyr::group_by(.data$lctn_pr, .data$gadm_level) %>%
    dplyr::summarise(.have = list(unique(.data$year)), .groups = "drop")
  key <- key %>%
    dplyr::left_join(yrs, by = c("lctn_pr", "gadm_level")) %>%
    dplyr::mutate(year = vapply(seq_len(dplyr::n()), function(i) {
      h <- .data$.have[[i]]
      if (is.null(h)) NA_integer_ else .nearestCrosswalkYear(.data$.want_year[i], h)
    }, integer(1))) %>%
    dplyr::select(-".have", -".want_year")

  # --- attach the member pieces ---------------------------------------------
  pieces <- xw %>%
    dplyr::select("lctn_pr", "gadm_level", "year", member_gid = "gid", "H") %>%
    dplyr::group_by(.data$lctn_pr, .data$gadm_level, .data$year) %>%
    # Renormalise over the surviving pieces rather than trusting a carried `w`:
    # the sliver filter in population_core.R can drop pieces after w was
    # derived, and conservation of cases depends on sum(w) == 1 EXACTLY.
    dplyr::mutate(apportion_w         = .data$H / sum(.data$H),
                  n_apportion_members = dplyr::n(),
                  P_l                 = sum(.data$H)) %>%
    dplyr::ungroup()

  matched <- key %>%
    dplyr::inner_join(pieces, by = c("lctn_pr", "gadm_level", "year"),
                      relationship = "many-to-many")

  hit_rows  <- sort(unique(matched$.row))
  miss_rows <- setdiff(seq_len(nrow(dat)), hit_rows)

  list(
    matched   = matched,
    miss_rows = miss_rows,
    report    = tibble::tibble(
      n_rows_in       = nrow(dat),
      n_rows_matched  = nrow(matched),
      n_lp_matched    = dplyr::n_distinct(matched$lctn_pr),
      n_lp_unmatched  = dplyr::n_distinct(as.character(dat[[lp_col]])[miss_rows]),
      max_members     = suppressWarnings(max(c(0L, matched$n_apportion_members),
                                             na.rm = TRUE)),
      keep_unmatched  = isTRUE(keep_unmatched)
    )
  )
}

#' lpCompositeStatus
#'
#' Per-INPUT-ROW composite classification, using the dominance criterion
#'   gid_cover = max_g w[l->g] ,   is_composite := gid_cover < cover_threshold
#' i.e. exactly `computeGidCover()` (population_core.R:662) and
#' `deriveCrosswalkShares()` (:646-654), but
#'
#'   (a) evaluated at RUN time against a caller-supplied `cover_threshold`,
#'       rather than read back from the `is_composite` column frozen into the
#'       crosswalk at 02b compile time; and
#'   (b) computed over the pieces that SURVIVE `.matchLpMembers()` (finite,
#'       positive H), so it describes the same member set the caller is about to
#'       act on. Zero-mass pieces cannot change `max(w)` anyway, so this agrees
#'       with `computeGidCover()` on a clean crosswalk -- it only diverges where
#'       the crosswalk itself is degenerate.
#'
#' TWO NON-EQUIVALENT NOTIONS OF "COMPOSITE" COEXIST IN THIS CODEBASE and must
#' not be conflated:
#'   is_composite (DOMINANCE, here)   -- gid_cover < 0.90; no member owns the LP
#'   apportioned  (MULTIPLICITY)      -- n_apportion_members > 1; any fringe at all
#' The second is far more common (a dominant gid plus a sliver is multi-member
#' but NOT composite). `--composite_policy drop` keys on the first.
#'
#' @inheritParams .matchLpMembers
#' @param cover_threshold gid_cover strictly below this => true composite.
#' @return tibble with one row per row of `dat`: `.row`, `lctn_pr`, `pop_year`,
#'   `P_l`, `n_member`, `gid_cover`, `n_gid_eff`, `is_composite`, `matched`.
#'   Unmatched rows get `matched = FALSE`, `is_composite = FALSE` (never
#'   dropped on a classification we could not make) and NA statistics.
lpCompositeStatus <- function(dat, xwalk,
                              lp_col     = "location_period_id",
                              year_col   = "year_seas",
                              level_col  = "admin_level",
                              gadm_level = "ADM1",
                              tiling     = c("derived", "corpus_label"),
                              cover_threshold = 0.90) {
  empty <- tibble::tibble(.row = integer(0), lctn_pr = character(0),
                          pop_year = integer(0), P_l = numeric(0),
                          n_member = integer(0), gid_cover = numeric(0),
                          n_gid_eff = numeric(0), is_composite = logical(0),
                          matched = logical(0))
  if (!nrow(dat)) return(empty)

  mm <- .matchLpMembers(dat, xwalk, lp_col = lp_col, year_col = year_col,
                        level_col = level_col, gadm_level = gadm_level,
                        tiling = tiling, .fn = "lpCompositeStatus")

  hit <- mm$matched %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::summarise(
      lctn_pr   = dplyr::first(.data$lctn_pr),
      pop_year  = dplyr::first(.data$year),
      P_l       = dplyr::first(.data$P_l),
      n_member  = as.integer(dplyr::first(.data$n_apportion_members)),
      gid_cover = suppressWarnings(max(.data$apportion_w, na.rm = TRUE)),
      n_gid_eff = { s <- sum(.data$apportion_w^2, na.rm = TRUE)
                    if (s > 0) 1 / s else NA_real_ },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gid_cover    = dplyr::if_else(is.finite(.data$gid_cover), .data$gid_cover, NA_real_),
      is_composite = !is.na(.data$gid_cover) & .data$gid_cover < cover_threshold,
      matched      = TRUE
    )

  miss <- tibble::tibble(
    .row      = as.integer(mm$miss_rows),
    lctn_pr   = as.character(dat[[lp_col]])[mm$miss_rows],
    pop_year  = NA_integer_,
    P_l       = NA_real_,
    n_member  = 1L,
    gid_cover = NA_real_,
    n_gid_eff = NA_real_,
    # An LP the crosswalk never saw is NOT declared composite: we have no basis
    # for the claim, and `drop` must not remove rows on an absent measurement.
    is_composite = FALSE,
    matched      = FALSE
  )

  dplyr::bind_rows(hit, miss) %>% dplyr::arrange(.data$.row)
}

#' collapseCompositeLPs
#'
#' The REDUCING counterpart to `apportionCompositeCases()`: instead of expanding
#' one location-period row into K member rows, keep exactly ONE row per input
#' row, relocated to the LP's dominant member gid and carrying the WHOLE LP.
#'
#' WHY THIS EXISTS. Expansion turns one real observation into K pseudo-
#' observations that the likelihood treats as K independent data points. On the
#' occurrence track the indicator is propagated, so "at least one member had an
#' outbreak" is fitted as "all K members had one"; on the size track cases and
#' population are scaled by the same w, so all K siblings carry an IDENTICAL
#' implied attack rate and row-wise PSIS-LOO is optimistic by a factor that
#' scales with K -- and K is geographically clustered. This function removes the
#' pseudo-replication at source.
#'
#' WHAT IT IS, IN ONE LINE. Legacy centroid matching, but pick the
#' dominant-POPULATION member instead of the polygon centroid, and use the
#' crosswalk population. A strict improvement on the legacy path, which put all
#' 7 Sudanese states' cases AND population on one arbitrary gid
#' (02_preprocess_crosswalk.R:6-24).
#'
#' NOTHING IS SCALED. The occurrence indicator and the size count both pass
#' through untouched, which is why there is no `cases_col` argument. The
#' denominator is `P_l` (the whole LP), not the dominant member's `H`, so the
#' implied attack rate stays exactly `Y_l / P_l` -- identical to what the LP
#' reported and to what `apportionCompositeCases()` yields per member.
#'
#' TRUE COMPOSITES ARE DROPPED, NOT COLLAPSED (`drop_composites = TRUE`). When
#' no member owns the LP (`gid_cover < cover_threshold`) there is no defensible
#' representative, so the honest treatment is to remove the observation. Callers
#' must remove the matching AT-RISK rows too, or dropping the numerator while
#' keeping the denominator biases the hazard downward.
#'
#' @inheritParams .matchLpMembers
#' @param drop_composites drop rows whose LP is a true composite by the
#'   dominance criterion (see `lpCompositeStatus()`).
#' @param cover_threshold gid_cover strictly below this => true composite.
#' @param gid_col column to overwrite with the dominant member gid. The incoming
#'   value is preserved as `gid_centroid`, same convention as
#'   `apportionCompositeCases()`.
#' @return `dat` with AT MOST one row per input row -- never more -- carrying
#'   `gid` (dominant member), `gid_centroid`, `pop` (= P_l), `pop_lp`,
#'   `pop_year`, `pop_tier`, plus `apportion_w = 1`, `n_apportion_members = 1L`,
#'   `apportioned = FALSE` for schema parity with the expansion path, plus
#'   `gid_cover`, `n_gid_eff`, `n_member_orig`, `collapsed`. A summary tibble is
#'   attached as attr(, "collapse_report") and the dropped LP ids as
#'   attr(, "collapse_dropped_lp").
collapseCompositeLPs <- function(dat, xwalk,
                                 drop_composites = TRUE,
                                 cover_threshold = 0.90,
                                 lp_col     = "location_period_id",
                                 year_col   = "year_seas",
                                 gid_col    = "gid",
                                 level_col  = "admin_level",
                                 gadm_level = "ADM1",
                                 tiling     = c("derived", "corpus_label"),
                                 keep_unmatched = TRUE) {
  stopifnot(is.data.frame(dat), is.data.frame(xwalk))
  need <- c("lctn_pr", "gid", "gadm_level", "year", "H")
  miss <- setdiff(need, names(xwalk))
  if (length(miss))
    stop("collapseCompositeLPs: crosswalk is missing column(s): ",
         paste(miss, collapse = ", "),
         ". Was it built by 02c_preprocess_population_compile.R?")
  if (!lp_col %in% names(dat))
    stop("collapseCompositeLPs: `dat` has no column '", lp_col, "'.")
  if (!is.numeric(cover_threshold) || length(cover_threshold) != 1L ||
      !is.finite(cover_threshold) || cover_threshold < 0 || cover_threshold > 1)
    stop("collapseCompositeLPs: cover_threshold must be a single number in [0, 1].")
  if (!nrow(dat)) return(dat)

  mm        <- .matchLpMembers(dat, xwalk, lp_col = lp_col, year_col = year_col,
                               level_col = level_col, gadm_level = gadm_level,
                               tiling = tiling, keep_unmatched = keep_unmatched,
                               .fn = "collapseCompositeLPs")
  matched   <- mm$matched
  miss_rows <- mm$miss_rows

  # --- one dominant member per input row ------------------------------------
  # Ties are broken by gid string order, so the choice is deterministic and
  # platform-independent -- never by whatever order the join happened to emit.
  dom <- matched %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::mutate(
      gid_cover = suppressWarnings(max(.data$apportion_w, na.rm = TRUE)),
      n_gid_eff = { s <- sum(.data$apportion_w^2, na.rm = TRUE)
                    if (s > 0) 1 / s else NA_real_ }
    ) %>%
    dplyr::arrange(dplyr::desc(.data$H), .data$member_gid, .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gid_cover    = dplyr::if_else(is.finite(.data$gid_cover), .data$gid_cover, NA_real_),
      is_composite = !is.na(.data$gid_cover) & .data$gid_cover < cover_threshold
    ) %>%
    dplyr::arrange(.data$.row)

  n_lp_composite <- dplyr::n_distinct(dom$lctn_pr[dom$is_composite])
  dropped_lp     <- sort(unique(as.character(dom$lctn_pr[dom$is_composite])))
  n_rows_dropped <- if (drop_composites) sum(dom$is_composite) else 0L

  if (drop_composites) dom <- dom[!dom$is_composite, , drop = FALSE]

  base_cols <- function(d) {
    d$gid_centroid <- if (gid_col %in% names(d)) as.character(d[[gid_col]]) else NA_character_
    d$pop_lp       <- if ("pop" %in% names(d)) as.numeric(d$pop) else NA_real_
    d
  }

  out_hit <- dat[dom$.row, , drop = FALSE] %>% base_cols()
  out_hit[[gid_col]]          <- dom$member_gid
  # P_l, not H: the collapsed row stands for the WHOLE location period, so the
  # attack rate must remain Y_l / P_l exactly.
  out_hit$pop                 <- dom$P_l
  # base_cols() captured pop_lp from `dat$pop` BEFORE `pop` above was
  # overwritten with P_l -- see the identical note in apportionCompositeCases().
  # Overwrite with the real P_l now that it is available; pop_lp == pop here by
  # construction (documented as such), kept as a separate column for schema
  # parity with the expansion path.
  out_hit$pop_lp               <- dom$P_l
  out_hit$apportion_w         <- 1
  out_hit$n_apportion_members <- 1L
  out_hit$apportioned         <- FALSE
  out_hit$pop_tier            <- "lp_collapsed"
  out_hit$pop_year            <- dom$year
  out_hit$gid_cover           <- dom$gid_cover
  out_hit$n_gid_eff           <- dom$n_gid_eff
  out_hit$n_member_orig       <- as.integer(dom$n_apportion_members)
  out_hit$collapsed           <- TRUE

  out_miss <- dat[miss_rows, , drop = FALSE] %>% base_cols()
  if (nrow(out_miss)) {
    out_miss$apportion_w         <- 1
    out_miss$n_apportion_members <- 1L
    out_miss$apportioned         <- FALSE
    out_miss$pop_tier            <- "unmatched_lp"
    out_miss$pop_year            <- NA_integer_
    out_miss$gid_cover           <- NA_real_
    out_miss$n_gid_eff           <- NA_real_
    out_miss$n_member_orig       <- 1L
    out_miss$collapsed           <- FALSE
    if (!"pop" %in% names(out_miss)) out_miss$pop <- NA_real_
  }

  out <- if (keep_unmatched) dplyr::bind_rows(out_hit, out_miss) else out_hit

  attr(out, "collapse_report") <- tibble::tibble(
    n_rows_in        = nrow(dat),
    n_rows_out       = nrow(out),
    # What --composite_policy keep would have produced, free of a second run.
    n_rows_if_keep   = nrow(matched) + (if (keep_unmatched) length(miss_rows) else 0L),
    n_lp_total       = dplyr::n_distinct(as.character(dat[[lp_col]])),
    n_lp_matched     = mm$report$n_lp_matched,
    n_lp_unmatched   = mm$report$n_lp_unmatched,
    n_lp_composite   = n_lp_composite,
    n_lp_dropped     = if (drop_composites) length(dropped_lp) else 0L,
    n_rows_dropped   = as.integer(n_rows_dropped),
    n_gid_out        = dplyr::n_distinct(out[[gid_col]]),
    max_members_orig = mm$report$max_members,
    cover_threshold  = cover_threshold,
    drop_composites  = isTRUE(drop_composites)
  )
  attr(out, "collapse_dropped_lp") <- if (drop_composites) dropped_lp else character(0)
  out
}
