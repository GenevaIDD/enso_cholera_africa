# season_calendar.R — per-country cholera-season calendar from each gid's
#   `min_seas`, and the hindcast / event window rules that use it
#
# Provides:
#   season_end_month(), season_end_by_country(), first_season_end()
#   event_window_end(), hindcast_window_end()   window-end rules
#   onset_rate_troughs(), season_trough_by_country(), season_window_of()
#                                       the --season_window rule
#   read_panel_season_calendar(), panel_season_calendar(), season_end_source()
#   season_end_spread_warnings(), write_season_ends_csv()   diagnostics, CSV
#   gid_min_seas_source()               per-gid calendar source (03 sidecar)
#
# Details:
#   `min_seas` is the first month of a gid's cholera year (03_prepare_occ.R).
#   Used by make_hindcast_configs.R and make_event_attribution_configs.R
#   (--season_end by_country_from_panel, --season_window). Pure definitions,
#   base R only (digest optional); sources nothing.
#   Gate: tests/test_season_calendar.R.
#
# See: the table above season_end_month() for every rule in one place

# Month index (year * 12 + month - 1) of a date, and back to a first-of-month Date.
.sc_month_index <- function(d) {
  d <- as.Date(d)
  as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m")) - 1L
}
.sc_from_index <- function(i) {
  i <- as.integer(i)
  as.Date(sprintf("%04d-%02d-01", i %/% 12L, i %% 12L + 1L))
}
.sc_check_month <- function(m, what) {
  m <- as.integer(m)
  if (!length(m) || anyNA(m) || any(m < 1L | m > 12L)) {
    stop(
      sprintf("season_calendar: %s must be months 1-12.", what),
      call. = FALSE
    )
  }
  m
}

# The rules in one place:
#
#   season_end_month(min_seas)          the end of the cholera year: min_seas - 1
#                                       (mod 12; min_seas 1 -> December)
#   season_end_by_country(panel)        gid -> country: the modal min_seas (ties
#                                       -> the later season end)
#   season_end_span(months)             circular spread of season-end months
#   first_season_end(from, end, min)    first `end` month at least `min` months
#                                       after `from`
#   event_window_end(peak, end, lags)   event rule: first end on or after the
#                                       peak + the largest positive lag
#   hindcast_window_end(launch, end, min_lag, min_window_months)
#                                       hindcast rule: first end at least
#                                       `min_window_months` after the launch
#                                       + max(min_lag, 0)
#   onset_rate_troughs(surv, ...)       per-country onset-rate trough + its
#                                       bootstrap stability (occurrence panel)
#   season_trough_by_country(cal, ...)  trough per country: onset rate where
#                                       supported, else the modal min_seas
#   season_window_of(trough, year, ...) --season_window rule: month after the
#                                       trough (floored at floor_month) to the
#                                       next trough + max(min_lag, 0)
#   season_cut_share(onsets, start, fl) onset share of the months the floor cuts
#   read_panel_season_calendar(path)    (gid, country, min_seas) of a sim panel,
#                                       with its path and content digest
#                                       (panel_season_calendar() for a panel in
#                                       memory; season_end_source() its record)
#   season_end_spread_warnings()        by_group diagnostic: a group whose
#                                       members' ends span > 4 months
#   season_ends_csv_path() / write_season_ends_csv()
#                                       the season_ends_<manifest stem>.csv
#                                       written next to each manifest
#   gid_min_seas_source(members, lp_source)
#                                       per-gid calendar SOURCE (seasonality
#                                       estimate vs the January default), for
#                                       03_prepare_occ.R's sidecar
#
#' season_end_month
#'
#' The last month of the cholera year that starts in `min_seas`:
#' min_seas - 1, wrapping 1 -> 12 (a year starting in January ends in December).
season_end_month <- function(min_seas) {
  m <- .sc_check_month(min_seas, "min_seas")
  ((m - 2L) %% 12L) + 1L
}

#' season_end_by_country
#'
#' One season end per country from a gid-level calendar. Each gid carries one
#' `min_seas` (a gid with several values, or in several countries, is refused);
#' a country takes its modal `min_seas` over its gids, ties going to the LATER
#' season end, i.e. the larger end month 1-12 (the longer window when a window
#' is cut at the first end after a date in the same year).
#'
#' @param panel data frame with `gid`, `min_seas` and a country column (ISO3);
#'   any number of rows per gid (a sim panel has one per month).
#' @param aggregate "mode" (the only rule implemented).
#' @param country_col name of the ISO3 column (default "country").
#' @return data.frame(iso3, season_end_month, min_seas_mode, modal_share, n_gid,
#'   n_distinct), one row per country, sorted by iso3. modal_share = share of the
#'   country's gids on the modal value; n_distinct = distinct min_seas values.
season_end_by_country <- function(
  panel,
  aggregate = "mode",
  country_col = "country"
) {
  if (!identical(aggregate, "mode")) {
    stop(
      "season_end_by_country: aggregate must be \"mode\" (the only rule implemented).",
      call. = FALSE
    )
  }
  need <- c("gid", "min_seas", country_col)
  miss <- setdiff(need, names(panel))
  if (length(miss)) {
    stop(
      "season_end_by_country: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  g <- unique(data.frame(
    gid = as.character(panel$gid),
    iso3 = as.character(panel[[country_col]]),
    min_seas = as.integer(panel$min_seas),
    stringsAsFactors = FALSE
  ))
  if (anyNA(g$min_seas) || anyNA(g$iso3)) {
    stop("season_end_by_country: NA min_seas or country.", call. = FALSE)
  }
  .sc_check_month(g$min_seas, "min_seas")
  dup <- unique(g$gid[duplicated(g$gid)])
  if (length(dup)) {
    stop(
      sprintf(
        "season_end_by_country: %d gid(s) carry several min_seas or countries, e.g. %s.",
        length(dup),
        paste(utils::head(dup, 3L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  rows <- lapply(split(g, g$iso3), function(x) {
    tab <- table(x$min_seas)
    ms <- as.integer(names(tab))
    n <- as.integer(tab)
    cand <- ms[n == max(n)]
    pick <- cand[which.max(season_end_month(cand))]
    data.frame(
      iso3 = x$iso3[1L],
      season_end_month = season_end_month(pick),
      min_seas_mode = pick,
      modal_share = max(n) / nrow(x),
      n_gid = nrow(x),
      n_distinct = length(ms),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$iso3), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' season_end_span
#'
#' Circular spread of a set of season-end months (1-12): the length, in months,
#' of the shortest arc of the year covering them all (0 for one month; {4, 11}
#' gives 5, the arc Nov-Apr).
season_end_span <- function(months) {
  m <- sort(unique(.sc_check_month(months, "season-end months")))
  if (length(m) < 2L) {
    return(0L)
  }
  gaps <- c(diff(m), m[1L] + 12L - m[length(m)])
  as.integer(12L - max(gaps))
}

#' first_season_end
#'
#' The first first-of-month date whose month is `end_month` and that falls at
#' least `min_months` months after `from` (min_months = 0: on or after `from`;
#' 1: strictly after a first-of-month `from`).
first_season_end <- function(from, end_month, min_months = 0L) {
  e <- .sc_check_month(end_month, "end_month")
  i0 <- .sc_month_index(from) + as.integer(min_months)
  .sc_from_index(i0 + (e - 1L - i0 %% 12L) %% 12L)
}

#' event_window_end
#'
#' Past-event rule: the season end is the first `end_month` on or after the
#' event peak; the window runs to it + the largest positive lag in `lags`.
#' @return list(season_end, horizon_end (Dates), max_lag)
event_window_end <- function(peak, end_month, lags) {
  se <- first_season_end(peak, end_month, 0L)
  max_lag <- max(c(as.integer(lags), 0L))
  list(
    season_end = se,
    horizon_end = .sc_from_index(.sc_month_index(se) + max_lag),
    max_lag = max_lag
  )
}

#' hindcast_window_end
#'
#' Hindcast rule under per-country season ends: the season end is the first
#' `end_month` at least `min_window_months` after the launch (cholera-year ends
#' sit close to a boreal-spring launch for many countries, and the "strictly
#' after" rule of the group season ends would give 1-3-month stubs); the window
#' runs to it + max(min_lag, 0) when `add_best_lag` (a lead never shortens the
#' window: the hindcast reads observed ENSO).
#' @return list(season_end, horizon_end) (Dates)
hindcast_window_end <- function(
  launch,
  end_month,
  min_lag,
  min_window_months = 6L,
  add_best_lag = TRUE
) {
  se <- first_season_end(launch, end_month, as.integer(min_window_months))
  ext <- if (isTRUE(add_best_lag)) max(as.integer(min_lag), 0L) else 0L
  list(season_end = se, horizon_end = .sc_from_index(.sc_month_index(se) + ext))
}

# ---------------------------------------------------------------------------
# Season windows (hindcast --season_window): one cholera season per
# country and season year, from the month after its seasonal TROUGH to the next
# trough, the start floored at `floor_month` (after the ONI spring
# predictability barrier). The trough is the country's onset-rate minimum where
# its occurrence panel supports one (enough onsets, a trough stable under a
# bootstrap over years), else its modal min_seas (the month of lowest cholera
# probability of the seasonality estimates, 03_prepare_occ.R).
# ---------------------------------------------------------------------------

# Centred circular moving average over k (odd) months; k = 1 is the identity.
.sc_circ_smooth <- function(v, k = 3L) {
  k <- as.integer(k)
  if (k <= 1L) {
    return(as.numeric(v))
  }
  h <- (k - 1L) %/% 2L
  idx <- outer(seq_len(12L), -h:h, function(i, j) ((i - 1L + j) %% 12L) + 1L)
  rowMeans(matrix(as.numeric(v)[idx], nrow = 12L))
}
# Circular distance between months (0-6).
.sc_circ_dist <- function(a, b) {
  d <- abs(as.integer(a) - as.integer(b)) %% 12L
  pmin(d, 12L - d)
}

#' onset_rate_troughs
#'
#' Per-country seasonal trough of the onset rate (onsets / at-risk unit-months
#' by calendar month, circularly smoothed over `smooth_months`), and its
#' stability: the share of `bootstrap` resamples of the country's years (with
#' replacement; seeded per country, so a country's value does not depend on the
#' others) whose trough lies within one month of the point estimate.
#' @param surv at-risk rows of an occurrence (survival) panel: t_start,
#'   outbreak (1 = onset) and a country column.
#' @param years integer vector of calendar years kept (NULL = all).
#' @return data.frame(iso3, n_onsets, n_atrisk, onset_trough, stability), one
#'   row per country (onset_trough / stability NA without onsets), with
#'   attribute `onsets_by_month` (ISO3 x 12 matrix of raw onset counts).
onset_rate_troughs <- function(
  surv,
  years = NULL,
  smooth_months = 3L,
  bootstrap = 500L,
  seed = 1L,
  country_col = "country"
) {
  miss <- setdiff(c("t_start", "outbreak", country_col), names(surv))
  if (length(miss)) {
    stop(
      "onset_rate_troughs: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  d <- as.Date(surv$t_start)
  x <- data.frame(
    iso3 = as.character(surv[[country_col]]),
    yr = as.integer(format(d, "%Y")),
    m = as.integer(format(d, "%m")),
    on = as.integer(!is.na(surv$outbreak) & surv$outbreak == 1),
    stringsAsFactors = FALSE
  )
  x <- x[!is.na(x$iso3) & !is.na(x$m), , drop = FALSE]
  if (!is.null(years)) {
    x <- x[x$yr %in% as.integer(years), , drop = FALSE]
  }
  trough_of <- function(on, ar) {
    which.min(.sc_circ_smooth(on / pmax(ar, 1), smooth_months))
  }
  isos <- sort(unique(x$iso3))
  onm <- matrix(0, length(isos), 12L, dimnames = list(isos, month.abb))
  rows <- lapply(isos, function(iso) {
    xi <- x[x$iso3 == iso, , drop = FALSE]
    yrs <- sort(unique(xi$yr))
    # year x month onset and at-risk counts
    on <- matrix(0, length(yrs), 12L)
    ar <- on
    k <- cbind(match(xi$yr, yrs), xi$m)
    for (r in seq_len(nrow(xi))) {
      on[k[r, 1L], k[r, 2L]] <- on[k[r, 1L], k[r, 2L]] + xi$on[r]
      ar[k[r, 1L], k[r, 2L]] <- ar[k[r, 1L], k[r, 2L]] + 1
    }
    onm[iso, ] <<- colSums(on)
    n_on <- sum(on)
    tr <- NA_integer_
    st <- NA_real_
    if (n_on > 0) {
      tr <- trough_of(colSums(on), colSums(ar))
      set.seed(as.integer(seed))
      bt <- vapply(
        seq_len(as.integer(bootstrap)),
        function(b) {
          s <- sample.int(length(yrs), length(yrs), replace = TRUE)
          trough_of(
            colSums(on[s, , drop = FALSE]),
            colSums(ar[s, , drop = FALSE])
          )
        },
        integer(1)
      )
      st <- if (length(bt)) mean(.sc_circ_dist(bt, tr) <= 1L) else NA_real_
    }
    data.frame(
      iso3 = iso,
      n_onsets = as.integer(n_on),
      n_atrisk = as.integer(sum(ar)),
      onset_trough = as.integer(tr),
      stability = st,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      iso3 = character(0),
      n_onsets = integer(0),
      n_atrisk = integer(0),
      onset_trough = integer(0),
      stability = numeric(0)
    )
  }
  rownames(out) <- NULL
  attr(out, "onsets_by_month") <- onm
  out
}

#' season_trough_by_country
#'
#' One trough month per country. rule "onset_rate": the onset-rate trough where
#' the country has >= min_onsets onsets and stability >= min_stability, else its
#' modal min_seas; rule "min_seas": always the modal min_seas.
#' @param cal season_end_by_country() output (iso3, min_seas_mode, modal_share)
#' @param onsets onset_rate_troughs() output (NULL under rule "min_seas")
#' @return data.frame(iso3, trough_month, trough_source, min_seas_mode,
#'   modal_share, n_onsets, onset_trough, stability), sorted by iso3, with the
#'   `onsets_by_month` attribute carried over.
season_trough_by_country <- function(
  cal,
  onsets = NULL,
  rule = c("onset_rate", "min_seas"),
  min_onsets = 30L,
  min_stability = 0.8
) {
  rule <- match.arg(rule)
  if (identical(rule, "onset_rate") && is.null(onsets)) {
    stop(
      "season_trough_by_country: rule \"onset_rate\" needs `onsets`.",
      call. = FALSE
    )
  }
  o <- if (is.null(onsets)) {
    NULL
  } else {
    onsets[match(cal$iso3, onsets$iso3), , drop = FALSE]
  }
  n_on <- if (is.null(o)) {
    rep(NA_integer_, nrow(cal))
  } else {
    as.integer(o$n_onsets)
  }
  o_tr <- if (is.null(o)) {
    rep(NA_integer_, nrow(cal))
  } else {
    as.integer(o$onset_trough)
  }
  st <- if (is.null(o)) rep(NA_real_, nrow(cal)) else as.numeric(o$stability)
  use <- identical(rule, "onset_rate") &
    !is.na(o_tr) &
    !is.na(n_on) &
    n_on >= min_onsets &
    !is.na(st) &
    st >= min_stability
  out <- data.frame(
    iso3 = cal$iso3,
    trough_month = ifelse(use, o_tr, as.integer(cal$min_seas_mode)),
    trough_source = ifelse(use, "onset_rate", "min_seas"),
    min_seas_mode = as.integer(cal$min_seas_mode),
    modal_share = cal$modal_share,
    n_onsets = n_on,
    onset_trough = o_tr,
    stability = st,
    stringsAsFactors = FALSE
  )
  out$trough_month <- as.integer(out$trough_month)
  out <- out[order(out$iso3), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "onsets_by_month") <- attr(onsets, "onsets_by_month")
  out
}

#' season_window_of
#'
#' The hindcast window of the season of `season_year` for a trough month: the
#' season starts the month after the trough; the launch is that month, or
#' `floor_month` when the start falls in January .. floor_month - 1 (the start
#' months before the floor are cut: `months_cut`); the season ends at the first
#' trough month after the launch; the firm horizon adds max(min_lag, 0) when
#' `add_best_lag` (a lead never shortens it: the hindcast reads observed ENSO).
#' @return list(start_month, launch_month, floored, months_cut, launch,
#'   season_end, horizon_end) (dates first-of-month)
season_window_of <- function(
  trough_month,
  season_year,
  floor_month = 5L,
  min_lag = 0L,
  add_best_lag = TRUE
) {
  tr <- .sc_check_month(trough_month, "trough_month")
  fl <- .sc_check_month(floor_month, "floor_month")
  start <- tr %% 12L + 1L
  floored <- start < fl
  lm <- if (floored) fl else start
  launch <- as.Date(sprintf("%04d-%02d-01", as.integer(season_year), lm))
  se <- first_season_end(launch, tr, 1L)
  ext <- if (isTRUE(add_best_lag)) max(as.integer(min_lag), 0L) else 0L
  list(
    start_month = start,
    launch_month = lm,
    floored = floored,
    months_cut = if (floored) fl - start else 0L,
    launch = launch,
    season_end = se,
    horizon_end = .sc_from_index(.sc_month_index(se) + ext)
  )
}

#' season_cut_share
#'
#' Share of a country's onsets (12 monthly counts) in the season months the
#' launch floor cuts (start_month .. floor_month - 1); 0 when nothing is cut,
#' NA without onsets.
season_cut_share <- function(onsets_by_month, start_month, floor_month = 5L) {
  v <- as.numeric(onsets_by_month)
  if (length(v) != 12L) {
    stop("season_cut_share: 12 monthly counts expected.", call. = FALSE)
  }
  if (sum(v) <= 0) {
    return(NA_real_)
  }
  if (start_month >= floor_month) {
    return(0)
  }
  sum(v[start_month:(floor_month - 1L)]) / sum(v)
}

#' panel_season_calendar / read_panel_season_calendar
#'
#' The gid calendar of a sim climate panel (03_prepare_sim_panel.R carries each
#' gid's min_seas and country): one row per gid, with attributes `path` and
#' `digest` (xxhash64 of the sorted gid / country / min_seas table, NA when the
#' digest package is absent) for the manifests' `season_end_source`.
#' panel_season_calendar() takes a panel already in memory.
panel_season_calendar <- function(p, path = NA_character_) {
  miss <- setdiff(c("gid", "country", "min_seas"), names(p))
  if (length(miss)) {
    stop(
      "panel_season_calendar: the climate panel ",
      path,
      " lacks column(s) ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  g <- unique(data.frame(
    gid = as.character(p$gid),
    country = as.character(p$country),
    min_seas = as.integer(p$min_seas),
    stringsAsFactors = FALSE
  ))
  g <- g[order(g$gid, g$country, g$min_seas), , drop = FALSE]
  rownames(g) <- NULL
  attr(g, "path") <- path
  attr(g, "digest") <- if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(g[, c("gid", "country", "min_seas")], algo = "xxhash64")
  } else {
    NA_character_
  }
  g
}
read_panel_season_calendar <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(
      "read_panel_season_calendar: climate panel not found: '",
      path %||% "",
      "'.",
      call. = FALSE
    )
  }
  panel_season_calendar(readRDS(path), path)
}

#' season_end_source
#'
#' The manifest record of where per-country season ends came from.
season_end_source <- function(cal) {
  list(
    panel = as.character(attr(cal, "path")),
    min_seas_digest = as.character(attr(cal, "digest")),
    rule = "season_end_month = min_seas - 1; gid -> country: modal min_seas, ties -> the later end"
  )
}

#' season_end_spread_warnings
#'
#' Under the group season ends (by_group): one message per fit group whose
#' member countries' cholera-year ends (from the panel calendar) span more than
#' `max_span` months of the year -- a group mixing hemispheres, whose single
#' season end misplaces some members' seasons. character(0) when none.
#' @param members named list group -> ISO3 vector
#' @param se_tbl season_end_by_country() output
season_end_spread_warnings <- function(members, se_tbl, max_span = 4L) {
  out <- character(0)
  for (g in names(members)) {
    x <- se_tbl[se_tbl$iso3 %in% members[[g]], , drop = FALSE]
    if (nrow(x) < 2L) {
      next
    }
    sp <- season_end_span(x$season_end_month)
    if (sp > max_span) {
      out <- c(
        out,
        sprintf(
          paste0(
            "group '%s': its countries' cholera-year ends span %d months ",
            "(%s); one group season end misplaces some of them ",
            "(see --season_end by_country_from_panel)."
          ),
          g,
          sp,
          paste(sprintf("%s=%d", x$iso3, x$season_end_month), collapse = ", ")
        )
      )
    }
  }
  out
}

#' season_ends_csv_path / write_season_ends_csv
#'
#' `season_ends_<manifest stem>.csv` next to a manifest: one row per country x
#' window (launch year or event) with the country's season end and its firm and
#' unit horizons, the per-country record behind the manifest's `countries` block.
season_ends_csv_path <- function(manifest) {
  file.path(
    dirname(manifest),
    sprintf("season_ends_%s.csv", sub("\\.[^.]*$", "", basename(manifest)))
  )
}
write_season_ends_csv <- function(rows, manifest) {
  p <- season_ends_csv_path(manifest)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, p, row.names = FALSE)
  invisible(p)
}

#' gid_min_seas_source
#'
#' Where each gid's season calendar comes from. `min_seas` is resolved per
#' location period in 03_prepare_occ.R (a seasonality estimate matched at the
#' gid or its ADM1, else the January default) and made one per gid by the
#' canonical rule (modal value over the gid's member rows, ties -> smallest
#' month; footprint.R::canonical_season_calendar()). This reproduces that value
#' and labels it.
#' @param members rows (gid, location_period_id, min_seas): the post-composite
#'   episode rows the calendar is resolved from.
#' @param lp_source (location_period_id, min_seas_source) with min_seas_source in
#'   estimate_gid | estimate_adm1 | default_january.
#' @return data.frame(gid, min_seas, min_seas_source, n_rows, n_rows_default):
#'   min_seas_source of the canonical value = the best source among the rows
#'   carrying it (estimate_gid > estimate_adm1 > default_january).
gid_min_seas_source <- function(members, lp_source) {
  lvl <- c("estimate_gid", "estimate_adm1", "default_january")
  m <- data.frame(
    gid = as.character(members$gid),
    location_period_id = as.character(members$location_period_id),
    min_seas = as.integer(members$min_seas),
    stringsAsFactors = FALSE
  )
  m <- m[!is.na(m$min_seas), , drop = FALSE]
  src <- stats::setNames(
    as.character(lp_source$min_seas_source),
    as.character(lp_source$location_period_id)
  )
  m$src <- unname(src[m$location_period_id])
  if (anyNA(m$src) || !all(m$src %in% lvl)) {
    stop(
      "gid_min_seas_source: location period(s) without a known min_seas source.",
      call. = FALSE
    )
  }
  rows <- lapply(split(m, m$gid), function(x) {
    tab <- table(x$min_seas)
    ms <- as.integer(names(tab))
    n <- as.integer(tab)
    canon <- min(ms[n == max(n)])
    s <- x$src[x$min_seas == canon]
    data.frame(
      gid = x$gid[1L],
      min_seas = canon,
      min_seas_source = lvl[min(match(s, lvl))],
      n_rows = nrow(x),
      n_rows_default = sum(x$src == "default_january"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      gid = character(0),
      min_seas = integer(0),
      min_seas_source = character(0),
      n_rows = integer(0),
      n_rows_default = integer(0)
    )
  }
  out <- out[order(out$gid), , drop = FALSE]
  rownames(out) <- NULL
  out
}
