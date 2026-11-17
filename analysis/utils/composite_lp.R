# composite_lp.R — classify composite location periods (LPs spanning several
#   GADM units) and collapse each LP to its dominant member gid.
#
# Provides:
#   lp_composite_status()     per-row gid_cover, n_gid_eff, is_composite
#   collapse_composite_lps()  one row per LP: dominant gid, true composites
#                             dropped (drop_composites = TRUE)
#
# Details:
#   Shares its member-matching core (.match_lp_members()) with
#   population.R::apportion_composite_cases(), so the two cannot drift apart.
#   Dropping a composite's cases requires the caller to drop its at-risk rows.

#' .match_lp_members
#'
#' The member-matching core shared by `apportion_composite_cases()` (one LP row
#' -> K member rows) and `collapse_composite_lps()` (one LP row -> its dominant
#' member), so the two policies pick tilings, snap years and renormalise shares
#' the same way. It does three things:
#'   1. per-row tiling (ADM1 vs ADM2): the crosswalk's `lp_level` under
#'      `tiling = "derived"`, the corpus `admin_level` under "corpus_label",
#'      falling back to whatever tiling the crosswalk carries;
#'   2. nearest-year snapping per (LP, tiling), `.nearest_crosswalk_year()`;
#'   3. the join to member pieces, renormalising `apportion_w = H / sum(H)`
#'      over the surviving pieces (the sliver filter can drop pieces, and
#'      conservation needs `sum(w) == 1` exactly).
#' It makes no decision about the match (no scaling, dropping or rewriting);
#' that is the caller's job.
#'
#' @param dat        tibble; one row per (location period, event).
#' @param xwalk      crosswalk from load_population_crosswalk().
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
.match_lp_members <- function(
  dat,
  xwalk,
  lp_col = "location_period_id",
  year_col = "year_seas",
  level_col = "admin_level",
  gadm_level = "ADM1",
  tiling = c("derived", "corpus_label"),
  keep_unmatched = TRUE,
  .fn = ".match_lp_members"
) {
  stopifnot(is.data.frame(dat), is.data.frame(xwalk))
  need <- c("lctn_pr", "gid", "gadm_level", "year", "H")
  miss <- setdiff(need, names(xwalk))
  if (length(miss)) {
    stop(
      .fn,
      ": crosswalk is missing column(s): ",
      paste(miss, collapse = ", "),
      ". Was it built by 02c_preprocess_population_compile.R?"
    )
  }
  if (!lp_col %in% names(dat)) {
    stop(.fn, ": `dat` has no column '", lp_col, "'.")
  }

  # --- per-row tiling -------------------------------------------------------
  # `corpus_label` is kept for comparison and is WRONG in general: the
  # outbreak corpus's admin depth is not GADM's. Malawi's corpus "ADM2" is a
  # district, which GADM calls ADM1 (GADM's ADM2 is the Traditional Authority),
  # so reading it at ADM2 shattered one district into 8-25 units and `gid_cover`
  # measured within-district population concentration instead of compositeness.
  # `derived` uses the geometry-based `lp_level` frozen onto the crosswalk by
  # 02c (derive_lp_tiling_level()). Kept side-by-side, not deleted, so the earlier
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
      dplyr::distinct(
        lctn_pr = as.character(.data$lctn_pr),
        lp_level = as.character(.data$lp_level)
      ) %>%
      dplyr::filter(!is.na(.data$lp_level))
    if (anyDuplicated(lp_lvl$lctn_pr)) {
      stop(
        .fn,
        ": the crosswalk carries more than one `lp_level` for some ",
        "location period(s) — it must be one decision per LP."
      )
    }
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
    dplyr::mutate(
      lctn_pr = as.character(.data$lctn_pr),
      gid = as.character(.data$gid),
      gadm_level = as.character(.data$gadm_level),
      year = as.integer(.data$year)
    ) %>%
    dplyr::filter(is.finite(.data$H), .data$H > 0)

  key <- tibble::tibble(
    .row = seq_len(nrow(dat)),
    lctn_pr = as.character(dat[[lp_col]]),
    gadm_level = row_level,
    .want_year = row_year
  )
  yrs <- xw %>%
    dplyr::group_by(.data$lctn_pr, .data$gadm_level) %>%
    dplyr::summarise(.have = list(unique(.data$year)), .groups = "drop")
  key <- key %>%
    dplyr::left_join(yrs, by = c("lctn_pr", "gadm_level")) %>%
    dplyr::mutate(
      year = vapply(
        seq_len(dplyr::n()),
        function(i) {
          h <- .data$.have[[i]]
          if (is.null(h)) {
            NA_integer_
          } else {
            .nearest_crosswalk_year(.data$.want_year[i], h)
          }
        },
        integer(1)
      )
    ) %>%
    dplyr::select(-".have", -".want_year")

  # --- attach the member pieces ---------------------------------------------
  pieces <- xw %>%
    dplyr::select("lctn_pr", "gadm_level", "year", member_gid = "gid", "H") %>%
    dplyr::group_by(.data$lctn_pr, .data$gadm_level, .data$year) %>%
    # Renormalise over the surviving pieces rather than trusting a carried `w`:
    # the sliver filter in population_core.R can drop pieces after w was
    # derived, and conservation of cases depends on sum(w) == 1 exactly.
    dplyr::mutate(
      apportion_w = .data$H / sum(.data$H),
      n_apportion_members = dplyr::n(),
      P_l = sum(.data$H)
    ) %>%
    dplyr::ungroup()

  matched <- key %>%
    dplyr::inner_join(
      pieces,
      by = c("lctn_pr", "gadm_level", "year"),
      relationship = "many-to-many"
    )

  hit_rows <- sort(unique(matched$.row))
  miss_rows <- setdiff(seq_len(nrow(dat)), hit_rows)

  list(
    matched = matched,
    miss_rows = miss_rows,
    report = tibble::tibble(
      n_rows_in = nrow(dat),
      n_rows_matched = nrow(matched),
      n_lp_matched = dplyr::n_distinct(matched$lctn_pr),
      n_lp_unmatched = dplyr::n_distinct(as.character(dat[[lp_col]])[
        miss_rows
      ]),
      max_members = suppressWarnings(max(
        c(0L, matched$n_apportion_members),
        na.rm = TRUE
      )),
      keep_unmatched = isTRUE(keep_unmatched)
    )
  )
}

#' lp_composite_status
#'
#' Per-INPUT-ROW composite classification, using the dominance criterion
#'   gid_cover = max_g w[l->g] ,   is_composite := gid_cover < cover_threshold
#' i.e. exactly `compute_gid_cover()` and
#' `derive_crosswalk_shares()` (population_core.R), but
#'
#'   (a) evaluated at RUN time against a caller-supplied `cover_threshold`,
#'       rather than read back from the `is_composite` column frozen into the
#'       crosswalk at 02b compile time; and
#'   (b) computed over the pieces that SURVIVE `.match_lp_members()` (finite,
#'       positive H), so it describes the same member set the caller is about to
#'       act on. Zero-mass pieces cannot change `max(w)` anyway, so this agrees
#'       with `compute_gid_cover()` on a clean crosswalk -- it only diverges where
#'       the crosswalk itself is degenerate.
#'
#' TWO NON-EQUIVALENT NOTIONS OF "COMPOSITE" COEXIST IN THIS CODEBASE and must
#' not be conflated:
#'   is_composite (DOMINANCE, here)   -- gid_cover < 0.90; no member owns the LP
#'   apportioned  (MULTIPLICITY)      -- n_apportion_members > 1; any fringe at all
#' The second is far more common (a dominant gid plus a sliver is multi-member
#' but NOT composite). `--composite_policy drop` keys on the first.
#'
#' @inheritParams .match_lp_members
#' @param cover_threshold gid_cover strictly below this => true composite.
#' @return tibble with one row per row of `dat`: `.row`, `lctn_pr`, `pop_year`,
#'   `P_l`, `n_member`, `gid_cover`, `n_gid_eff`, `is_composite`, `matched`.
#'   Unmatched rows get `matched = FALSE`, `is_composite = FALSE` (never
#'   dropped on a classification we could not make) and NA statistics.
lp_composite_status <- function(
  dat,
  xwalk,
  lp_col = "location_period_id",
  year_col = "year_seas",
  level_col = "admin_level",
  gadm_level = "ADM1",
  tiling = c("derived", "corpus_label"),
  cover_threshold = 0.90
) {
  empty <- tibble::tibble(
    .row = integer(0),
    lctn_pr = character(0),
    pop_year = integer(0),
    P_l = numeric(0),
    n_member = integer(0),
    gid_cover = numeric(0),
    n_gid_eff = numeric(0),
    is_composite = logical(0),
    matched = logical(0)
  )
  if (!nrow(dat)) {
    return(empty)
  }

  mm <- .match_lp_members(
    dat,
    xwalk,
    lp_col = lp_col,
    year_col = year_col,
    level_col = level_col,
    gadm_level = gadm_level,
    tiling = tiling,
    .fn = "lp_composite_status"
  )

  hit <- mm$matched %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::summarise(
      lctn_pr = dplyr::first(.data$lctn_pr),
      pop_year = dplyr::first(.data$year),
      P_l = dplyr::first(.data$P_l),
      n_member = as.integer(dplyr::first(.data$n_apportion_members)),
      gid_cover = suppressWarnings(max(.data$apportion_w, na.rm = TRUE)),
      n_gid_eff = {
        s <- sum(.data$apportion_w^2, na.rm = TRUE)
        if (s > 0) 1 / s else NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gid_cover = dplyr::if_else(
        is.finite(.data$gid_cover),
        .data$gid_cover,
        NA_real_
      ),
      is_composite = !is.na(.data$gid_cover) &
        .data$gid_cover < cover_threshold,
      matched = TRUE
    )

  miss <- tibble::tibble(
    .row = as.integer(mm$miss_rows),
    lctn_pr = as.character(dat[[lp_col]])[mm$miss_rows],
    pop_year = NA_integer_,
    P_l = NA_real_,
    n_member = 1L,
    gid_cover = NA_real_,
    n_gid_eff = NA_real_,
    # An LP the crosswalk never saw is NOT declared composite: we have no basis
    # for the claim, and `drop` must not remove rows on an absent measurement.
    is_composite = FALSE,
    matched = FALSE
  )

  dplyr::bind_rows(hit, miss) %>% dplyr::arrange(.data$.row)
}

#' collapse_composite_lps
#'
#' The reducing counterpart to `apportion_composite_cases()`: one row per input
#' row, relocated to the LP's dominant-population member gid and carrying the
#' whole LP. This avoids the K pseudo-observations that expansion creates.
#' Nothing is scaled: the occurrence indicator and the case count pass through,
#' and the denominator is `P_l`, so the attack rate stays `Y_l / P_l`.
#' True composites (`gid_cover < cover_threshold`) are dropped, not collapsed
#' (`drop_composites = TRUE`); callers must drop the matching at-risk rows too,
#' or the hazard is biased downward.
#'
#' @inheritParams .match_lp_members
#' @param drop_composites drop rows whose LP is a true composite by the
#'   dominance criterion (see `lp_composite_status()`).
#' @param cover_threshold gid_cover strictly below this => true composite.
#' @param gid_col column to overwrite with the dominant member gid. The incoming
#'   value is preserved as `gid_centroid`, same convention as
#'   `apportion_composite_cases()`.
#' @return `dat` with AT MOST one row per input row -- never more -- carrying
#'   `gid` (dominant member), `gid_centroid`, `pop` (= P_l), `pop_lp`,
#'   `pop_year`, `pop_tier`, plus `apportion_w = 1`, `n_apportion_members = 1L`,
#'   `apportioned = FALSE` for schema parity with the expansion path, plus
#'   `gid_cover`, `n_gid_eff`, `n_member_orig`, `collapsed`. A summary tibble is
#'   attached as attr(, "collapse_report") and the dropped LP ids as
#'   attr(, "collapse_dropped_lp").
collapse_composite_lps <- function(
  dat,
  xwalk,
  drop_composites = TRUE,
  cover_threshold = 0.90,
  lp_col = "location_period_id",
  year_col = "year_seas",
  gid_col = "gid",
  level_col = "admin_level",
  gadm_level = "ADM1",
  tiling = c("derived", "corpus_label"),
  keep_unmatched = TRUE
) {
  stopifnot(is.data.frame(dat), is.data.frame(xwalk))
  need <- c("lctn_pr", "gid", "gadm_level", "year", "H")
  miss <- setdiff(need, names(xwalk))
  if (length(miss)) {
    stop(
      "collapse_composite_lps: crosswalk is missing column(s): ",
      paste(miss, collapse = ", "),
      ". Was it built by 02c_preprocess_population_compile.R?"
    )
  }
  if (!lp_col %in% names(dat)) {
    stop("collapse_composite_lps: `dat` has no column '", lp_col, "'.")
  }
  if (
    !is.numeric(cover_threshold) ||
      length(cover_threshold) != 1L ||
      !is.finite(cover_threshold) ||
      cover_threshold < 0 ||
      cover_threshold > 1
  ) {
    stop(
      "collapse_composite_lps: cover_threshold must be a single number in [0, 1]."
    )
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
    .fn = "collapse_composite_lps"
  )
  matched <- mm$matched
  miss_rows <- mm$miss_rows

  # --- one dominant member per input row ------------------------------------
  # Ties are broken by gid string order, so the choice is deterministic and
  # platform-independent -- never by whatever order the join happened to emit.
  dom <- matched %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::mutate(
      gid_cover = suppressWarnings(max(.data$apportion_w, na.rm = TRUE)),
      n_gid_eff = {
        s <- sum(.data$apportion_w^2, na.rm = TRUE)
        if (s > 0) 1 / s else NA_real_
      }
    ) %>%
    dplyr::arrange(dplyr::desc(.data$H), .data$member_gid, .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gid_cover = dplyr::if_else(
        is.finite(.data$gid_cover),
        .data$gid_cover,
        NA_real_
      ),
      is_composite = !is.na(.data$gid_cover) & .data$gid_cover < cover_threshold
    ) %>%
    dplyr::arrange(.data$.row)

  n_lp_composite <- dplyr::n_distinct(dom$lctn_pr[dom$is_composite])
  dropped_lp <- sort(unique(as.character(dom$lctn_pr[dom$is_composite])))
  n_rows_dropped <- if (drop_composites) sum(dom$is_composite) else 0L

  if (drop_composites) {
    dom <- dom[!dom$is_composite, , drop = FALSE]
  }

  base_cols <- function(d) {
    d$gid_centroid <- if (gid_col %in% names(d)) {
      as.character(d[[gid_col]])
    } else {
      NA_character_
    }
    d$pop_lp <- if ("pop" %in% names(d)) as.numeric(d$pop) else NA_real_
    d
  }

  out_hit <- dat[dom$.row, , drop = FALSE] %>% base_cols()
  out_hit[[gid_col]] <- dom$member_gid
  # P_l, not H: the collapsed row stands for the whole location period, so the
  # attack rate must remain Y_l / P_l exactly.
  out_hit$pop <- dom$P_l
  # base_cols() captured pop_lp from `dat$pop` BEFORE `pop` above was
  # overwritten with P_l -- see the identical note in apportion_composite_cases().
  # Overwrite with the real P_l now that it is available; pop_lp == pop here by
  # construction (documented as such), kept as a separate column for schema
  # parity with the expansion path.
  out_hit$pop_lp <- dom$P_l
  out_hit$apportion_w <- 1
  out_hit$n_apportion_members <- 1L
  out_hit$apportioned <- FALSE
  out_hit$pop_tier <- "lp_collapsed"
  out_hit$pop_year <- dom$year
  out_hit$gid_cover <- dom$gid_cover
  out_hit$n_gid_eff <- dom$n_gid_eff
  out_hit$n_member_orig <- as.integer(dom$n_apportion_members)
  out_hit$collapsed <- TRUE

  out_miss <- dat[miss_rows, , drop = FALSE] %>% base_cols()
  if (nrow(out_miss)) {
    out_miss$apportion_w <- 1
    out_miss$n_apportion_members <- 1L
    out_miss$apportioned <- FALSE
    out_miss$pop_tier <- "unmatched_lp"
    out_miss$pop_year <- NA_integer_
    out_miss$gid_cover <- NA_real_
    out_miss$n_gid_eff <- NA_real_
    out_miss$n_member_orig <- 1L
    out_miss$collapsed <- FALSE
    if (!"pop" %in% names(out_miss)) out_miss$pop <- NA_real_
  }

  out <- if (keep_unmatched) dplyr::bind_rows(out_hit, out_miss) else out_hit

  attr(out, "collapse_report") <- tibble::tibble(
    n_rows_in = nrow(dat),
    n_rows_out = nrow(out),
    # What --composite_policy keep would have produced, free of a second run.
    n_rows_if_keep = nrow(matched) +
      (if (keep_unmatched) length(miss_rows) else 0L),
    n_lp_total = dplyr::n_distinct(as.character(dat[[lp_col]])),
    n_lp_matched = mm$report$n_lp_matched,
    n_lp_unmatched = mm$report$n_lp_unmatched,
    n_lp_composite = n_lp_composite,
    n_lp_dropped = if (drop_composites) length(dropped_lp) else 0L,
    n_rows_dropped = as.integer(n_rows_dropped),
    n_gid_out = dplyr::n_distinct(out[[gid_col]]),
    max_members_orig = mm$report$max_members,
    cover_threshold = cover_threshold,
    drop_composites = isTRUE(drop_composites)
  )
  attr(out, "collapse_dropped_lp") <- if (drop_composites) {
    dropped_lp
  } else {
    character(0)
  }
  out
}
