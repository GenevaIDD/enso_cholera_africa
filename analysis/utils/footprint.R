# footprint.R — the (gid, year_seas, admin_level) footprint collapse and the
#   canonical episode table shared by the size, duration and monthly tracks.
#
# Provides:
#   FOOTPRINT_KEY, footprint_key()   the collapse key, named once
#   unify_footprint_ids(), collapse_to_footprint()   one row per footprint
#   footprint_member_rows()          which LPs fed each footprint
#   aggregate_monthly_to_footprint() monthly series on the footprint key
#   build_canonical_episodes(), canonical_episodes()   cross-track episodes
#   canonical_season_calendar(), representative_ids_by_gid()
#   adm0_episode_counts(), discarded_onset_table()
#
# Details:
#   The size frame must stay bit-identical so archived fits remain comparable:
#   unify_footprint_ids() / collapse_to_footprint() are not to be re-specified;
#   footprint_member_rows() is purely additive. Callers source base.R; dplyr /
#   tibble via `::`.

# Why footprints: 03_prepare_siz.R collapses episodes to one row per
# (gid, year_seas, admin_level), summing cases over every location period (LP)
# that resolves to that gid and season-year. The surviving uid /
# location_period_id is a representative label, not the source of the cases.
# So anything that attaches a per-LP quantity to a size/duration episode (e.g.
# the 01b monthly case series) must join on the footprint key through the
# membership map this module returns, never on (uid, location_period_id).

# The collapse key. Named once so a consumer cannot spell it differently from
# the producer -- the exact failure this module exists to prevent.
FOOTPRINT_KEY <- c("gid", "year_seas", "admin_level")


#' .footprint_key_chr
#'
#' Coerce the three key columns to character so a factor `admin_level` on one
#' side and a character `admin_level` on the other cannot silently fail to
#' match. 01_preprocess_outbreaks.R makes `admin_level` a FACTOR
#' (`factor(admin_level, str_c("ADM", 0:3))`), 01b_preprocess_outbreaks_monthly.R
#' makes it a character, and dplyr's join would either error or coerce
#' depending on version. Normalise at every boundary instead of trusting it.
#'
#' @param d data frame carrying FOOTPRINT_KEY
#' @return `d` with the three key columns as character
.footprint_key_chr <- function(
  d,
  .fn = ".footprint_key_chr",
  key = FOOTPRINT_KEY
) {
  miss <- setdiff(key, names(d))
  if (length(miss)) {
    stop(
      .fn,
      ": missing footprint key column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  for (k in key) {
    d[[k]] <- as.character(d[[k]])
  }
  d
}


#' unify_footprint_ids
#'
#' LEGACY-ONLY id unification, moved verbatim from 03_prepare_siz.R.
#' Fixes one representative location_period_id / uid per (gid, admin_level)
#' across all season-years, so the per-location_period_id
#' `time_to_last_outbreak` recursion downstream tracks each footprint's full
#' cross-season series instead of restarting per location period.
#'
#' WHY IT IS LEGACY-ONLY (the original comment, condensed -- the full version
#' stays at the call site because it documents a decision about that pipeline,
#' not about this function). Under `--pop_source legacy`, `gid` comes from a
#' stable geometric containment relationship (ADM3 centroid -> its ADM2
#' parent), so rows sharing a gid really are repeat reports of the same
#' footprint. Under `--pop_source crosswalk`, `gid` is instead the crosswalk's
#' dominant-population member chosen per location period independently, which
#' can coincide across genuinely unrelated location periods; forcing one id onto
#' those would merge distinct outbreaks' identities (DRC 14287/14296/14360/
#' 14372/14398 all resolve to COD.3.6_1). Callers must therefore gate this on
#' `!use_xwalk` -- the guard is left with the caller rather than hidden in here,
#' so the condition stays visible in the pipeline script.
#'
#' @param dat episode table carrying `gid`, `admin_level`, `location_period_id`,
#'   `uid`
#' @return `dat` with `location_period_id` / `uid` replaced by the per-(gid,
#'   admin_level) representative
unify_footprint_ids <- function(dat) {
  need <- c("gid", "admin_level", "location_period_id", "uid")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "unify_footprint_ids: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  dat %>%
    dplyr::group_by(gid, admin_level) %>%
    dplyr::mutate(
      .rep_i = order(location_period_id)[1],
      location_period_id = location_period_id[.rep_i],
      uid = uid[.rep_i]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.rep_i)
}


#' footprint_member_rows
#'
#' Capture which (uid, location_period_id) rows feed each footprint, before the
#' collapse loses that information. `dat` itself is not modified.
#'
#' Must be called BEFORE `unify_footprint_ids()`: under `--pop_source legacy`
#' that overwrites `uid` / `location_period_id` with one representative, so a
#' later map would assign every season's cases to one location period.
#'
#' `.member_i == 1L` marks the row `dplyr::first()` picks in
#' collapse_to_footprint(). `apportion_w` (w = H/P_l; 1 under legacy) is
#' carried because an external per-LP series is not pre-scaled: summed over
#' the K member rows of a `--composite_policy keep` LP it would count K times.
#'
#' @param dat episode table, post-climate-join, pre-collapse
#' @return tibble, one row per input row: FOOTPRINT_KEY (character) +
#'   `uid`, `location_period_id`, `TL`, `TR`, `suspected_cases`, `apportion_w`,
#'   `.member_i`, `n_members`
footprint_member_rows <- function(dat) {
  need <- c(FOOTPRINT_KEY, "uid", "location_period_id")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "footprint_member_rows: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  out <- tibble::tibble(
    gid = as.character(dat$gid),
    year_seas = as.character(dat$year_seas),
    admin_level = as.character(dat$admin_level),
    uid = as.character(dat$uid),
    location_period_id = as.character(dat$location_period_id),
    TL = if ("TL" %in% names(dat)) as.Date(dat$TL) else as.Date(NA),
    TR = if ("TR" %in% names(dat)) as.Date(dat$TR) else as.Date(NA),
    suspected_cases = if ("suspected_cases" %in% names(dat)) {
      as.numeric(dat$suspected_cases)
    } else {
      NA_real_
    },
    # 1 under legacy / collapse / drop_collapse; H/P_l under `keep`.
    apportion_w = if ("apportion_w" %in% names(dat)) {
      as.numeric(dat$apportion_w)
    } else {
      1
    }
  )

  out %>%
    dplyr::group_by(gid, year_seas, admin_level) %>%
    dplyr::mutate(.member_i = dplyr::row_number(), n_members = dplyr::n()) %>%
    dplyr::ungroup()
}


#' collapse_to_footprint
#'
#' The footprint collapse, moved verbatim from 03_prepare_siz.R.
#'
#' CASES ARE A FLOW: summing them WITHIN an admin level is correct (ADM3
#' children -> footprint total at ADM3 resolution; ADM2 report -> footprint
#' total at ADM2 resolution). Summing ACROSS admin levels would double-count,
#' which is why `admin_level` is part of the key.
#'
#' POPULATION IS A STOCK and `sum(pop)` here is NOT the final answer -- the
#' caller replaces it with the coverage-saturated min(sum_l H, Q_g) via
#' `reconcile_footprint_pop()`. That step stays at the call site: it needs the
#' population tables, which this module deliberately does not depend on.
#'
#' @param dat episode table carrying FOOTPRINT_KEY, `suspected_cases`, `pop`
#' @return one row per footprint
collapse_to_footprint <- function(dat) {
  need <- c(FOOTPRINT_KEY, "suspected_cases", "pop")
  miss <- setdiff(need, names(dat))
  if (length(miss)) {
    stop(
      "collapse_to_footprint: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  dat %>%
    dplyr::group_by(gid, year_seas, admin_level) %>%
    dplyr::summarise(
      suspected_cases = sum(suspected_cases, na.rm = TRUE),
      # sum_l H[l,g] under --pop_source crosswalk (pop is H per member after
      # apportionment); the legacy sum(pop) otherwise.
      pop = sum(pop, na.rm = TRUE),
      dplyr::across(-c(suspected_cases, pop), dplyr::first),
      .groups = "drop"
    )
}


#' aggregate_monthly_to_footprint
#'
#' Push a per-(uid, location_period_id, month) series through the same
#' membership the footprint collapse used, so the result is on the same
#' footprint as the episode's `suspected_cases`, `pop` and `exposure`.
#'
#' THE CONSERVATION GATE IS THE POINT. `sum_m cases_month` over a footprint must
#' equal that footprint's collapsed `suspected_cases`, because both are sums of
#' the same member location periods' cases -- one grouped by month, one not.
#' Any mismatch means a member's monthly series is missing or came from a
#' different corpus vintage, and the honest response is to say which footprint
#' and stop, not to coalesce the gap to zero. This is returned rather than
#' enforced here so the caller decides fatality (and can report before dying).
#'
#' @param monthly   per-(uid, location_period_id, month) tibble from
#'   `make_monthly_cases_file()`; needs `uid`, `location_period_id`, `month`,
#'   `cases_month`, and optionally `obs_days`, `n_intervals`
#' @param members   output of `footprint_member_rows()`
#' @param footprints collapsed frame (output of `collapse_to_footprint()`, or any
#'   frame with FOOTPRINT_KEY + `suspected_cases`) used for the conservation
#'   check. NULL to skip the check.
#' @param tolerance relative tolerance for conservation
#' @param key the episode key (`footprint_key()`); FOOTPRINT_KEY by default, so
#'   every existing call is unchanged. `c("gid", "year_seas")` under
#'   `episodes.key: gid_season`.
#' @return list(monthly = <footprint x month tibble>, recon = <per-footprint
#'   reconciliation tibble>, report = <one-row summary tibble>)
aggregate_monthly_to_footprint <- function(
  monthly,
  members,
  footprints = NULL,
  tolerance = 1e-6,
  key = FOOTPRINT_KEY
) {
  need_m <- c("uid", "location_period_id", "month", "cases_month")
  miss <- setdiff(need_m, names(monthly))
  if (length(miss)) {
    stop(
      "aggregate_monthly_to_footprint: `monthly` is missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(c(key, "uid", "location_period_id") %in% names(members))) {
    stop(
      "aggregate_monthly_to_footprint: `members` is not a footprint_member_rows() ",
      "output.",
      call. = FALSE
    )
  }

  # Select down to exactly the columns this function uses BEFORE the join.
  # `monthly` (01b's output) carries its own `admin_level` (and `location`,
  # `country`) that would otherwise collide with the footprint's `admin_level`
  # from `members` and get suffixed to `admin_level.x`/`admin_level.y` --
  # silently breaking the group_by() below instead of erroring loudly, except
  # that dplyr does error on group_by() not finding the bare name, which is
  # how this was caught.
  mon <- monthly %>%
    dplyr::transmute(
      uid = as.character(uid),
      location_period_id = as.character(location_period_id),
      month = as.Date(month),
      cases_month = cases_month,
      obs_days = if ("obs_days" %in% names(monthly)) {
        obs_days
      } else {
        NA_real_
      },
      n_intervals = if ("n_intervals" %in% names(monthly)) {
        n_intervals
      } else {
        NA_integer_
      }
    )

  # Members whose location period has NO monthly series at all. Counted and
  # returned, never silently dropped: this is the corpus-vintage tripwire.
  mem_key <- paste(members$uid, members$location_period_id)
  mon_key <- unique(paste(mon$uid, mon$location_period_id))
  n_member_unmatched <- sum(!mem_key %in% mon_key)

  # `relationship = "many-to-many"` is DECLARED, not accidental. One member row
  # matches that LP's several months (one-to-many), and under
  # --composite_policy keep one (uid, lp) can also appear as K member rows of K
  # different footprints (many-to-one). Both are intended; declaring it stops
  # dplyr warning about a fan-out that is the whole mechanism.
  joined <- members %>%
    dplyr::select(dplyr::all_of(key), uid, location_period_id, apportion_w) %>%
    dplyr::inner_join(
      mon,
      by = c("uid", "location_period_id"),
      relationship = "many-to-many"
    )

  out <- joined %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(key, "month")))) %>%
    dplyr::summarise(
      # apportion_w is 1 unless --composite_policy keep split this LP across
      # member gids; then the monthly series must be split by the same w the
      # episode total was, or the footprint totals stop matching.
      cases_month = sum(cases_month * apportion_w, na.rm = TRUE),
      obs_days = sum(obs_days, na.rm = TRUE),
      n_lp = dplyr::n_distinct(location_period_id),
      n_intervals = sum(n_intervals, na.rm = TRUE),
      .groups = "drop"
    )

  recon <- NULL
  if (!is.null(footprints)) {
    fp <- .footprint_key_chr(
      footprints,
      "aggregate_monthly_to_footprint",
      key
    ) %>%
      dplyr::select(dplyr::all_of(key), episode_cases = suspected_cases)
    recon <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
      dplyr::summarise(monthly_cases = sum(cases_month), .groups = "drop") %>%
      dplyr::full_join(fp, by = key) %>%
      dplyr::mutate(
        monthly_cases = dplyr::coalesce(monthly_cases, 0),
        episode_cases = dplyr::coalesce(episode_cases, 0),
        abs_err = abs(monthly_cases - episode_cases),
        rel_err = abs_err / pmax(1, abs(episode_cases)),
        conserved = rel_err <= tolerance
      )
  }

  report <- tibble::tibble(
    n_members = nrow(members),
    n_member_unmatched = as.integer(n_member_unmatched),
    n_footprints_in = nrow(dplyr::distinct(members[, key, drop = FALSE])),
    n_footprints_out = nrow(dplyr::distinct(out[, key, drop = FALSE])),
    n_rows_out = nrow(out),
    cases_out = sum(out$cases_month),
    n_not_conserved = if (is.null(recon)) {
      NA_integer_
    } else {
      as.integer(sum(!recon$conserved))
    },
    max_rel_err = if (is.null(recon)) NA_real_ else max(recon$rel_err)
  )

  list(monthly = out, recon = recon, report = report)
}


# ===========================================================================
# Canonical gid-season episodes (`episodes.key: gid_season`)
# ===========================================================================
#
# One table keyed (gid, year_seas), built by one function that
# 03_prepare_occ.R and 03_prepare_siz.R call with the same arguments, so every
# track trains on the same episodes:
#   onset = earliest member TL; end = latest member TR;
#   cases = the finest-resolution footprint total: members grouped by
#           reporting level (`report_level`), summed within a level, largest
#           total wins (ties -> finer); the per-episode form of
#           spatpomp_core.R::.aggregate_footprint_counts(), so the monthly
#           conservation gate stays exact;
#   admin_level = the geometry-derived GADM level of the gid (caller-supplied
#           `level_of`); the corpus label only decides nesting;
#   uid / location_period_id = one representative per gid across seasons
#           (uid <-> gid 1:1, asserted).
# Base-only: no population/provenance calls. The caller supplies `level_of`
# and computes the digest (provenance.R).

#' footprint_key
#'
#' The episode key for a given `episodes.key`. `lp_season` keeps
#' FOOTPRINT_KEY exactly, so every existing caller is unchanged.
footprint_key <- function(episode_key = "lp_season") {
  if (identical(episode_key, "gid_season")) {
    c("gid", "year_seas")
  } else {
    FOOTPRINT_KEY
  }
}

.ADMIN_LEVEL_RANK <- c(ADM0 = 0L, ADM1 = 1L, ADM2 = 2L, ADM3 = 3L)

.admin_rank <- function(level) {
  r <- unname(.ADMIN_LEVEL_RANK[as.character(level)])
  r[is.na(r)] <- -1L
  r
}

# `%m+%` (lubridate::add_with_rollback) by a per-row number of months -- the
# same arithmetic as occ_data.R::apply_offset_months(), vectorised over k.
.shift_months <- function(d, k) {
  d <- as.Date(d)
  k <- as.integer(k)
  out <- d
  for (u in unique(k[!is.na(k)])) {
    i <- which(k == u)
    out[i] <- lubridate::add_with_rollback(d[i], lubridate::period(u, "month"))
  }
  out
}

#' canonical_season_calendar
#'
#' One cholera-season calendar per gid. `min_seas` is resolved per LOCATION
#' PERIOD upstream (from the LP's centroid seasonality match), and under the
#' crosswalk collapse `gid` is overwritten with the dominant member while
#' `min_seas` is not -- so two LPs on one gid can sit on different `year_seas`
#' calendars, and "the same gid-season" would not even be well defined.
#' Canonical `min_seas` = the gid's most frequent value over its member rows
#' (ties -> smallest month). TL_seas / TR_seas / year_seas / month_seas are
#' re-derived for every row (a no-op where min_seas did not change), and the
#' first-outbreak-per-(LP, gid, season) rule of section C is re-applied,
#' because re-dating can move two outbreaks of one LP into the same season.
#'
#' @param members post-composite-reduction episode rows (`ob_ts3`).
#' @return `members`, re-dated, with attributes `n_min_seas_reassigned` and
#'   `n_redate_dropped`.
canonical_season_calendar <- function(members) {
  need <- c("gid", "location_period_id", "uid", "TL", "TR", "min_seas")
  miss <- setdiff(need, names(members))
  if (length(miss)) {
    stop(
      "canonical_season_calendar: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  canon <- members %>%
    dplyr::count(gid, min_seas, name = ".n") %>%
    dplyr::arrange(gid, dplyr::desc(.n), min_seas) %>%
    dplyr::distinct(gid, .keep_all = TRUE) %>%
    dplyr::select(gid, .min_seas_canon = min_seas)

  out <- members %>%
    dplyr::left_join(canon, by = "gid")
  n_reassigned <- sum(out$min_seas != out$.min_seas_canon, na.rm = TRUE)

  out <- out %>%
    dplyr::mutate(
      min_seas = .min_seas_canon,
      TL_seas = .shift_months(TL, -(min_seas - 1L)),
      TR_seas = .shift_months(TR, -(min_seas - 1L)),
      month_seas = lubridate::month(TL_seas),
      year_seas = lubridate::year(TL_seas)
    ) %>%
    dplyr::select(-.min_seas_canon)
  if ("TL_seas_min" %in% names(out)) {
    out <- out %>%
      dplyr::group_by(location_period_id) %>%
      dplyr::mutate(TL_seas_min = min(TL_seas), TL_seas_max = max(TL_seas)) %>%
      dplyr::ungroup()
  }

  n_in <- nrow(out)
  out <- out %>%
    dplyr::arrange(TL, uid) %>%
    dplyr::distinct(location_period_id, gid, year_seas, .keep_all = TRUE)

  attr(out, "n_min_seas_reassigned") <- as.integer(n_reassigned)
  attr(out, "n_redate_dropped") <- as.integer(n_in - nrow(out))
  out
}

#' finest_level_totals
#'
#' Per key: sum `value_col` within each level, then keep the level with the
#' largest total (ties -> finer level). Returns the total, the level used and
#' how many distinct levels the key had. Same value as
#' spatpomp_core.R::.aggregate_footprint_counts(..., keys = keys) -- pinned by
#' analysis/tests/test_episode_key.R -- but also says WHICH level won, which
#' the episode builder needs to flag the members it counted.
finest_level_totals <- function(
  df,
  value_col,
  keys,
  level_col = "admin_level"
) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(keys, level_col)))) %>%
    dplyr::summarise(
      .v = sum(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(.rank = .admin_rank(.data[[level_col]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::arrange(dplyr::desc(.v), dplyr::desc(.rank), .by_group = TRUE) %>%
    dplyr::summarise(
      total = dplyr::first(.v),
      level_used = as.character(dplyr::first(.data[[level_col]])),
      n_levels = dplyr::n(),
      .groups = "drop"
    )
}

#' finest_level_rows
#'
#' The rows of `df` that finest_level_totals() counts: per key, only the rows of
#' the level it picks (largest total, ties -> finer). Row order and columns are
#' kept, so a sum over the result runs over the same rows in the same order as a
#' sum over `df` wherever each key has one level (bitwise equal there). Used for
#' observed monthly cases summed above the gid (08c, 09 smoothing), where
#' summing every row double-counts nested reports of one gid at two levels.
#' Without `level_col` in `df`, `df` is returned unchanged.
finest_level_rows <- function(df, value_col, keys, level_col = "admin_level") {
  if (!level_col %in% names(df)) {
    return(df)
  }
  used <- finest_level_totals(df, value_col, keys, level_col) %>%
    dplyr::select(dplyr::all_of(keys), .level_used = level_used)
  df %>%
    dplyr::mutate(.level_used = as.character(.data[[level_col]])) %>%
    dplyr::semi_join(used, by = c(keys, ".level_used")) %>%
    dplyr::select(-".level_used")
}

#' representative_ids_by_gid
#'
#' One (uid, location_period_id) per gid across all seasons: the row with the
#' smallest location_period_id -- exactly the rule 03_prepare_occ.R has always
#' used before its gid union (`order(location_period_id)[1]`), now shared by
#' both prep scripts.
representative_ids_by_gid <- function(dat) {
  dat %>%
    dplyr::group_by(gid) %>%
    dplyr::summarise(
      .rep_i = order(location_period_id)[1],
      rep_location_period_id = location_period_id[.rep_i],
      rep_uid = uid[.rep_i],
      .groups = "drop"
    ) %>%
    dplyr::select(-.rep_i)
}

#' canonical_episodes
#'
#' @param members re-dated member rows (canonical_season_calendar() output):
#'   gid, year_seas, uid, location_period_id, admin_level (the member's
#'   REPORTING level, corpus label), TL, TR, TL_seas, TR_seas, `cases_col`,
#'   optionally `pop_col` and `apportion_w`.
#' @param level_of function(gid) -> canonical GADM level (no NA allowed).
#' @return list(episodes, members, report). `episodes`: one row per
#'   (gid, year_seas) with every member column (taken from the ONSET member),
#'   overwritten by the canonical TL/TR/TL_seas/TR_seas/cases/pop,
#'   `admin_level` (geometry), `admin_level_corpus` (the counted reporting
#'   level), representative uid/location_period_id, `onset_lp`,
#'   `censor_lp_id`, `n_members`, `n_levels`. `members`: the footprint_member_rows()
#'   schema keyed on the episode, plus `report_level` and `in_total`.
canonical_episodes <- function(
  members,
  level_of,
  cases_col = "suspected_cases",
  pop_col = "pop"
) {
  need <- c(
    "gid",
    "year_seas",
    "uid",
    "location_period_id",
    "admin_level",
    "TL",
    "TR",
    "TL_seas",
    "TR_seas",
    cases_col
  )
  miss <- setdiff(need, names(members))
  if (length(miss)) {
    stop(
      "canonical_episodes: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  key <- c("gid", "year_seas")
  has_pop <- pop_col %in% names(members)

  m <- members %>%
    dplyr::mutate(
      report_level = as.character(admin_level),
      TL = as.Date(TL),
      TR = as.Date(TR),
      TL_seas = as.Date(TL_seas),
      TR_seas = as.Date(TR_seas)
    )
  lvl <- finest_level_totals(m, cases_col, key, "report_level")
  m <- m %>%
    dplyr::left_join(
      lvl %>% dplyr::select(dplyr::all_of(key), level_used, n_levels),
      by = key
    ) %>%
    dplyr::mutate(in_total = report_level == level_used)

  # Onset member: earliest TL; among ties prefer a counted member, then the
  # finer reporting level, then uid -- deterministic.
  m <- m %>%
    dplyr::mutate(.rank = .admin_rank(report_level)) %>%
    dplyr::arrange(
      gid,
      year_seas,
      TL,
      dplyr::desc(in_total),
      dplyr::desc(.rank),
      uid
    )

  agg <- m %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
    dplyr::summarise(
      .TL = min(TL),
      .TR = max(TR),
      .TL_seas = min(TL_seas),
      .TR_seas = max(TR_seas),
      .cases = sum(.data[[cases_col]][in_total], na.rm = TRUE),
      .pop = if (has_pop) {
        sum(.data[[pop_col]][in_total], na.rm = TRUE)
      } else {
        NA_real_
      },
      onset_lp = dplyr::first(location_period_id),
      # censoring is a property of the member whose report ENDS the episode;
      # ties broken on the smallest location_period_id.
      censor_lp_id = location_period_id[order(
        -as.numeric(TR),
        as.character(location_period_id)
      )[1]],
      n_members = dplyr::n(),
      admin_level_corpus = dplyr::first(level_used),
      # reporting gap inside the episode: a member starting more than a month
      # after every earlier member has ended.
      .gap = {
        o <- order(TL)
        tl <- TL[o]
        tr <- TR[o]
        length(tl) > 1L &&
          any(tl[-1] > cummax(as.numeric(tr))[-length(tr)] + 31)
      },
      .groups = "drop"
    )

  onset <- m %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(key)), .keep_all = TRUE) %>%
    dplyr::select(-in_total, -.rank, -report_level)

  ep <- onset %>%
    dplyr::inner_join(agg, by = key) %>%
    dplyr::mutate(
      TL = .TL,
      TR = .TR,
      TL_seas = .TL_seas,
      TR_seas = .TR_seas,
      month_seas = lubridate::month(TL_seas),
      admin_level = as.character(level_of(gid))
    )
  ep[[cases_col]] <- ep$.cases
  if (has_pop) {
    ep[[pop_col]] <- ep$.pop
  }
  ep <- ep %>% dplyr::select(-.TL, -.TR, -.TL_seas, -.TR_seas, -.cases, -.pop)

  rep_ids <- representative_ids_by_gid(m)
  ep <- ep %>%
    dplyr::left_join(rep_ids, by = "gid") %>%
    dplyr::mutate(
      uid = rep_uid,
      location_period_id = rep_location_period_id
    ) %>%
    dplyr::select(-rep_uid, -rep_location_period_id)
  if ("TL_seas_min" %in% names(ep)) {
    ep <- ep %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(TL_seas_min = min(TL_seas), TL_seas_max = max(TL_seas)) %>%
      dplyr::ungroup()
  }

  # --- invariants -----------------------------------------------------------
  if (anyDuplicated(ep[, key])) {
    stop("canonical_episodes: (gid, year_seas) is not unique.", call. = FALSE)
  }
  if (any(ep$TL > ep$TR, na.rm = TRUE)) {
    stop(
      "canonical_episodes: ",
      sum(ep$TL > ep$TR, na.rm = TRUE),
      " episode(s) end before they start.",
      call. = FALSE
    )
  }
  uid_gids <- tapply(ep$gid, ep$uid, function(g) length(unique(g)))
  if (any(uid_gids > 1L)) {
    stop(
      "canonical_episodes: representative uid <-> gid is not 1:1 (",
      sum(uid_gids > 1L),
      " uid(s) represent several gids, e.g. ",
      paste(utils::head(names(uid_gids)[uid_gids > 1L], 3), collapse = ", "),
      "). This happens when one location period spans several gids ",
      "(--composite_policy keep/drop); use collapse or drop_collapse.",
      call. = FALSE
    )
  }

  # --- members, in the footprint_member_rows() schema ------------------------------
  mem <- m %>%
    dplyr::left_join(
      ep %>%
        dplyr::select(dplyr::all_of(key), episode_admin_level = admin_level),
      by = key
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
    dplyr::mutate(.member_i = dplyr::row_number(), n_members = dplyr::n()) %>%
    dplyr::ungroup()
  mem <- tibble::tibble(
    gid = as.character(mem$gid),
    year_seas = as.character(mem$year_seas),
    admin_level = as.character(mem$episode_admin_level),
    uid = as.character(mem$uid),
    location_period_id = as.character(mem$location_period_id),
    TL = as.Date(mem$TL),
    TR = as.Date(mem$TR),
    suspected_cases = as.numeric(mem[[cases_col]]),
    apportion_w = if ("apportion_w" %in% names(mem)) {
      dplyr::coalesce(as.numeric(mem$apportion_w), 1)
    } else {
      1
    },
    report_level = mem$report_level,
    in_total = mem$in_total,
    .member_i = mem$.member_i,
    n_members = mem$n_members
  )

  yr_tr <- lubridate::year(ep$TR_seas)
  nxt <- ep %>%
    dplyr::arrange(gid, TL) %>%
    dplyr::group_by(gid) %>%
    dplyr::mutate(.overlap = dplyr::lead(TL) <= TR) %>%
    dplyr::ungroup()
  report <- list(
    n_episodes_canonical = nrow(ep),
    n_members = nrow(mem),
    n_multi_member_episodes = sum(ep$n_members > 1L),
    n_multi_level_episodes = sum(ep$n_levels > 1L),
    n_gap_episodes = sum(ep$.gap),
    n_cross_season_episodes = sum(
      yr_tr > as.integer(ep$year_seas),
      na.rm = TRUE
    ),
    n_overlapping_next = sum(nxt$.overlap, na.rm = TRUE),
    cases_members = sum(as.numeric(members[[cases_col]]), na.rm = TRUE),
    cases_episodes = sum(ep[[cases_col]], na.rm = TRUE)
  )
  ep <- ep %>% dplyr::select(-.gap)

  list(episodes = ep, members = mem, report = report)
}

#' build_canonical_episodes
#'
#' THE shared entry point: 03_prepare_occ.R and 03_prepare_siz.R call this
#' with identical arguments on their post-composite-reduction `ob_ts3`, so the
#' two tracks (and duration/monthly, which inherit size's episodes) build the
#' same episode table -- the episode digest then proves it.
#'
#' @param ob_ts3 post-C2 episode rows (one per LP-event per member gid).
#' @param level_of function(gid) -> canonical level (population.R::canonical_level_lookup()).
#' @return canonical_episodes() output, with `report` extended by the calendar
#'   counts.
build_canonical_episodes <- function(
  ob_ts3,
  level_of,
  cases_col = "suspected_cases",
  pop_col = "pop"
) {
  cal <- canonical_season_calendar(ob_ts3)
  res <- canonical_episodes(
    cal,
    level_of,
    cases_col = cases_col,
    pop_col = pop_col
  )
  res$report <- c(
    res$report,
    list(
      n_min_seas_reassigned = attr(cal, "n_min_seas_reassigned"),
      n_redate_dropped = attr(cal, "n_redate_dropped")
    )
  )
  res$member_rows <- cal
  res
}

# ---------------------------------------------------------------------------
# Discarded onsets and the ADM0 drop
# ---------------------------------------------------------------------------

#' adm0_episode_counts
#'
#' National (ADM0) database outbreak rows, which 03_prepare_occ.R drops before
#' any panel is built: per country, the number of rows and their cases.
#' @param ob_ts outbreak rows (country, admin_level, suspected_cases).
#' @return tibble(country, n_adm0_episodes, cases_adm0).
adm0_episode_counts <- function(ob_ts) {
  ob_ts %>%
    dplyr::ungroup() %>%
    dplyr::filter(as.character(admin_level) == "ADM0") %>%
    dplyr::group_by(country = as.character(country)) %>%
    dplyr::summarise(
      n_adm0_episodes = dplyr::n(),
      cases_adm0 = sum(as.numeric(suspected_cases), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(country)
}

#' discarded_onset_table
#'
#' Every database onset that does not become a canonical episode onset of the
#' occurrence panel, one row each, by reason:
#'   lp_second_onset    a later onset of the same location period in the same
#'                      LP-calendar season-year, removed by the first-onset
#'                      filter (03 section C); `after_episode_end` TRUE when it
#'                      starts after the kept episode ended (a second wave);
#'                      `lp_calendar_only` TRUE when it falls in another season-
#'                      year of the gid's CANONICAL calendar than the kept one
#'                      (discarded only because of the LP calendar);
#'   redate_duplicate   dropped as a duplicate after re-dating onto the gid's
#'                      canonical calendar (canonical_season_calendar());
#'   absorbed_member    a later member of a canonical episode (its onset is the
#'                      episode's earliest member onset);
#'   absorbed_after_gap the same, starting more than 31 days after every
#'                      earlier member ended (a re-onset after a gap).
#'
#' @param ob_ts2 rows before the first-onset filter (location_period_id, gid,
#'   country, TL, TR, suspected_cases, year_seas).
#' @param ob_ts3_first rows right after it (same columns).
#' @param ob_ts3_final the rows build_canonical_episodes() was given (after the
#'   composite reduction; rows with a NA min_seas are ignored).
#' @param ep_res build_canonical_episodes() output (episodes, members, member_rows).
#' @return tibble(reason, location_period_id, gid, country, TL, TR, cases,
#'   year_seas_lp, year_seas_gid, after_episode_end, lp_calendar_only).
discarded_onset_table <- function(ob_ts2, ob_ts3_first, ob_ts3_final, ep_res) {
  k <- function(d) {
    paste(d$location_period_id, d$gid, as.Date(d$TL), as.Date(d$TR))
  }
  norm <- function(d) {
    d %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        location_period_id = as.character(location_period_id),
        gid = as.character(gid),
        country = as.character(country),
        TL = as.Date(TL),
        TR = as.Date(TR)
      )
  }
  a <- norm(ob_ts2)
  b <- norm(ob_ts3_first)
  f <- norm(ob_ts3_final)
  ms <- ep_res$episodes %>%
    dplyr::distinct(gid = as.character(gid), min_seas = as.integer(min_seas))
  gyear <- function(gid, TL) {
    m <- ms$min_seas[match(gid, ms$gid)]
    ifelse(
      is.na(m),
      NA_integer_,
      as.integer(lubridate::year(as.Date(TL) %m+% months(-(m - 1L))))
    )
  }
  # 1. LP-calendar first-onset filter
  lp2 <- a[!k(a) %in% k(b), , drop = FALSE]
  kept <- b %>%
    dplyr::group_by(location_period_id, year_seas) %>%
    dplyr::summarise(k_TL = min(TL), k_TR = max(TR), .groups = "drop")
  lp2 <- lp2 %>%
    dplyr::left_join(kept, by = c("location_period_id", "year_seas"))
  t1 <- tibble::tibble(
    reason = "lp_second_onset",
    location_period_id = lp2$location_period_id,
    gid = lp2$gid,
    country = lp2$country,
    TL = lp2$TL,
    TR = lp2$TR,
    cases = as.numeric(lp2$suspected_cases),
    year_seas_lp = as.integer(lp2$year_seas),
    year_seas_gid = gyear(lp2$gid, lp2$TL),
    after_episode_end = !is.na(lp2$k_TR) & lp2$TL > lp2$k_TR,
    lp_calendar_only = !is.na(lp2$k_TL) &
      gyear(lp2$gid, lp2$TL) != gyear(lp2$gid, lp2$k_TL)
  )
  # 2. re-dating duplicates (rows of ob_ts3 the canonical calendar dropped)
  cal <- norm(ep_res$member_rows)
  b_ok <- f[!is.na(f$min_seas), , drop = FALSE]
  rd <- b_ok[!k(b_ok) %in% k(cal), , drop = FALSE]
  t2 <- tibble::tibble(
    reason = "redate_duplicate",
    location_period_id = rd$location_period_id,
    gid = rd$gid,
    country = rd$country,
    TL = rd$TL,
    TR = rd$TR,
    cases = as.numeric(rd$suspected_cases),
    year_seas_lp = as.integer(rd$year_seas),
    year_seas_gid = gyear(rd$gid, rd$TL),
    after_episode_end = NA,
    lp_calendar_only = NA
  )
  # 3. absorbed members of a canonical episode
  mem <- ep_res$members %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gid = as.character(gid),
      TL = as.Date(TL),
      TR = as.Date(TR)
    ) %>%
    dplyr::arrange(gid, year_seas, TL) %>%
    dplyr::group_by(gid, year_seas) %>%
    dplyr::mutate(
      .on_TL = min(TL),
      .prev_TR = dplyr::lag(cummax(as.numeric(TR))),
      .gap = !is.na(.prev_TR) & as.numeric(TL) > .prev_TR + 31
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(TL > .on_TL)
  ctry <- ep_res$episodes %>%
    dplyr::distinct(gid = as.character(gid), country = as.character(country))
  t3 <- tibble::tibble(
    reason = as.character(ifelse(
      mem$.gap,
      "absorbed_after_gap",
      "absorbed_member"
    )),
    location_period_id = as.character(mem$location_period_id),
    gid = mem$gid,
    country = ctry$country[match(mem$gid, ctry$gid)],
    TL = mem$TL,
    TR = mem$TR,
    cases = as.numeric(mem$suspected_cases),
    year_seas_lp = NA_integer_,
    year_seas_gid = as.integer(mem$year_seas),
    after_episode_end = NA,
    lp_calendar_only = NA
  )
  dplyr::bind_rows(t1, t2, t3) %>% dplyr::arrange(reason, country, gid, TL)
}
