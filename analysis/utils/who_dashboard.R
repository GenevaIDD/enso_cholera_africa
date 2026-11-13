# WHO Global Cholera & AWD dashboard: national weekly counts ------------------
#
# Library for the WHO dashboard weekly series
# (`data/who_annual_report_africa/who_dashboard_weekly_africa.rds`, one row per
# country x epi week, national level only, 2023 onwards) and for the
# year-to-date cumulative snapshots transcribed from the WHO multi-country
# cholera epidemiological updates
# (`data/who_annual_report_africa/who_sitrep_cumulative_africa.csv`, one row
# per country x report cut), which fill the months the weekly feed lacks. Used by:
#   - 03_prepare_who_monthly.R   (weekly -> monthly, cumulative -> monthly,
#                                 gap-fill merge; written once)
#   - smoothing_core.R           (observation table + coverage fraction + per-unit
#                                 share for the `_who` smoothing measurement)
#   - 08c_smoothing_compile.R    (WHO vs smoothing comparison)
# Sources nothing (utils module purity); needs dplyr at call time.

#' whoWeeklyToMonthly
#'
#' Weekly national counts -> calendar-month counts by day overlap: each week is
#' the 7 days starting at `date_wk`, each day carries `cases / 7`, and days are
#' summed by calendar month (the same day-overlap apportionment as
#' 01b_preprocess_outbreaks_monthly.R).
#'
#' `date_wk` is a Monday for almost every row; SDN dates its weeks on Sundays
#' (a few Saturdays) and SSD has a handful of odd weekdays. They are all taken
#' as day 1 of their week: a one-day shift is immaterial at monthly grain.
#'
#' A week that is absent, or present with NA cases, is NOT reported: its days
#' do not count towards `n_days`. The dashboard reports zero weeks explicitly,
#' so an unreported week must never be read as zero. A month is complete when
#' every one of its days is covered by a reported week (`partial = FALSE`).
#'
#' @param wk data frame with `iso_3_code`, `date_wk` (Date), `cases`.
#' @param duplicates "error" (default): two rows for one (country, week) with
#'   different cases stop; exact duplicates are collapsed. "sum": add them.
#' @return tibble(iso3, t_start, cases, n_days, partial), one row per country
#'   and calendar month with at least one reported day. Months with no reported
#'   day are absent (callers treat absence as missing, never as zero).
whoWeeklyToMonthly <- function(wk, duplicates = c("error", "sum")) {
  duplicates <- match.arg(duplicates)
  miss <- setdiff(c("iso_3_code", "date_wk", "cases"), names(wk))
  if (length(miss))
    stop("whoWeeklyToMonthly: missing column(s): ", paste(miss, collapse = ", "))
  wk <- dplyr::transmute(wk, iso3 = as.character(iso_3_code),
                         date_wk = as.Date(date_wk), cases = as.numeric(cases))
  wk <- wk[!is.na(wk$cases), , drop = FALSE]
  wk <- dplyr::distinct(wk)
  dup <- duplicated(wk[, c("iso3", "date_wk")])
  if (any(dup)) {
    if (duplicates == "error")
      stop(sprintf("whoWeeklyToMonthly: %d (country, week) pair(s) reported twice with different cases (e.g. %s %s); pass duplicates = \"sum\" to add them.",
                   sum(dup), wk$iso3[dup][1], format(wk$date_wk[dup][1])))
    wk <- dplyr::ungroup(dplyr::summarise(dplyr::group_by(wk, iso3, date_wk),
                                          cases = sum(cases), .groups = "drop"))
  }
  if (nrow(wk) == 0L) return(.who_empty_monthly())
  .who_intervals_to_monthly(dplyr::tibble(iso3 = wk$iso3, start = wk$date_wk,
                                          end = wk$date_wk + 6L, cases = wk$cases))
}

.who_empty_monthly <- function()
  dplyr::tibble(iso3 = character(0), t_start = as.Date(character(0)),
                cases = numeric(0), n_days = integer(0), partial = logical(0))

#' .who_intervals_to_monthly
#'
#' Day-overlap apportionment shared by whoWeeklyToMonthly() and
#' whoCumulativeToMonthly(): each reported interval `[start, end]` (both
#' inclusive) spreads `cases` evenly over its days, and days are summed by
#' calendar month. `n_days` counts the distinct covered days, `partial` flags a
#' month with an uncovered day.
#'
#' @param iv tibble(iso3, start, end, cases); Dates, end >= start.
#' @return tibble(iso3, t_start, cases, n_days, partial), sorted.
.who_intervals_to_monthly <- function(iv) {
  if (nrow(iv) == 0L) return(.who_empty_monthly())
  len <- as.integer(iv$end - iv$start) + 1L
  if (any(len < 1L)) stop(".who_intervals_to_monthly: end before start.")
  days <- data.frame(
    iso3  = rep(as.character(iv$iso3), times = len),
    day   = as.Date(unlist(Map(function(s, n) s + seq_len(n) - 1L, iv$start, len)),
                    origin = "1970-01-01"),
    cases = rep(as.numeric(iv$cases) / len, times = len),
    stringsAsFactors = FALSE)
  days$t_start <- as.Date(format(days$day, "%Y-%m-01"))
  # Overlapping intervals (e.g. a Sunday-dated week next to a Monday-dated one)
  # put two reports on one day: count the day once for coverage, keep both cases.
  # group_by(), not `.by =`: the cluster's R 4.2.1 library may predate dplyr 1.1.
  out <- dplyr::summarise(dplyr::group_by(days, iso3, t_start),
                          cases  = sum(cases),
                          n_days = dplyr::n_distinct(day),
                          .groups = "drop")
  out$partial <- out$n_days < .who_days_in_month(out$t_start)
  dplyr::arrange(dplyr::as_tibble(out), iso3, t_start)
}

.who_days_in_month <- function(t_start) {
  nxt <- as.Date(vapply(t_start, function(d)
    format(seq(d, by = "month", length.out = 2L)[2L]), character(1)))
  as.integer(nxt - t_start)
}

#' whoCumulativeToMonthly
#'
#' Calendar-month national counts from the year-to-date cumulative snapshots of
#' the WHO multi-country cholera epidemiological updates (Table 1: "1 January
#' to <cut>" cases, plus the "last 28 days" column). Per country, consecutive
#' cuts bound an interval `(prev_cut, cut]` whose cases are the difference of
#' the two cumulatives; the first cut of a year (or a country's first cut) runs
#' from 1 January, since the cumulative restarts there. Interval cases are then
#' spread over days and summed by calendar month exactly as the weekly series
#' (.who_intervals_to_monthly()).
#'
#' An unchanged cumulative is ambiguous: WHO's footnote says missing data do
#' not imply absence, and a country that stopped reporting keeps its last
#' total. It is read as a reported zero only when the country is in
#' `zero_iso3` (a stated per-country decision, e.g. TZA after the Ministry
#' declared the outbreak controlled on 25 Feb 2026) or when the table printed an
#' explicit 0 in the last-28-days column; otherwise the interval is dropped
#' (missing, never zero). A negative difference (WHO revises totals
#' retrospectively) is dropped with a warning.
#'
#' An interval longer than `max_interval_days` is dropped too: a country's
#' first appearance in a mid-year table, or a lone end-of-year cut, carries a
#' year-to-date total with no monthly information, and spreading it evenly
#' would invent months. The updates are issued every 4-5 weeks, so 45 days
#' keeps every consecutive-cut interval and the January cut.
#'
#' @param snap data frame with `iso3`, `table_as_of` (Date), `cum_cases`,
#'   `last28_cases` (NA for "-").
#' @param zero_iso3 ISO3 codes whose unchanged cumulatives are reported zeros.
#' @param max_interval_days longest interval (days) spread over months.
#' @return tibble(iso3, t_start, cases, n_days, partial), as whoWeeklyToMonthly().
whoCumulativeToMonthly <- function(snap, zero_iso3 = character(0), max_interval_days = 45L) {
  miss <- setdiff(c("iso3", "table_as_of", "cum_cases", "last28_cases"), names(snap))
  if (length(miss))
    stop("whoCumulativeToMonthly: missing column(s): ", paste(miss, collapse = ", "))
  s <- dplyr::transmute(snap, iso3 = as.character(iso3), cut = as.Date(table_as_of),
                        cum = as.numeric(cum_cases), last28 = as.numeric(last28_cases))
  s <- s[!is.na(s$iso3) & !is.na(s$cut) & !is.na(s$cum), , drop = FALSE]
  s <- dplyr::distinct(s)
  dup <- duplicated(s[, c("iso3", "cut")])
  if (any(dup))
    stop(sprintf("whoCumulativeToMonthly: %d (country, cut) pair(s) with different totals (e.g. %s %s).",
                 sum(dup), s$iso3[dup][1], format(s$cut[dup][1])))
  if (nrow(s) == 0L) return(.who_empty_monthly())
  s <- dplyr::arrange(s, iso3, cut)
  s$year <- as.integer(format(s$cut, "%Y"))
  prev_same <- c(FALSE, s$iso3[-1] == s$iso3[-nrow(s)] & s$year[-1] == s$year[-nrow(s)])
  s$start <- as.Date(ifelse(prev_same, c(NA, s$cut[-nrow(s)]) + 1L,
                            as.Date(sprintf("%d-01-01", s$year))), origin = "1970-01-01")
  s$cases <- ifelse(prev_same, s$cum - c(NA, s$cum[-nrow(s)]), s$cum)
  neg <- s$cases < 0
  if (any(neg)) {
    warning(sprintf("whoCumulativeToMonthly: %d interval(s) dropped, cumulative revised down: %s",
                    sum(neg), paste(unique(sprintf("%s %s", s$iso3[neg], format(s$cut[neg]))),
                                    collapse = ", ")))
    s <- s[!neg, , drop = FALSE]
  }
  s <- s[as.integer(s$cut - s$start) + 1L <= as.integer(max_interval_days), , drop = FALSE]
  zero_ok <- s$iso3 %in% zero_iso3 | (!is.na(s$last28) & s$last28 == 0)
  s <- s[s$cases > 0 | zero_ok, , drop = FALSE]
  .who_intervals_to_monthly(dplyr::tibble(iso3 = s$iso3, start = s$start,
                                          end = s$cut, cases = s$cases))
}

#' whoMergeMonthly
#'
#' Fill the weekly-dashboard monthly table from the cumulative-derived one:
#' every weekly row is kept as is (a partial weekly month still wins), and a
#' cumulative row enters only for a (country, month) with no weekly row. Adds
#' `source` ("dashboard_weekly" / "sitrep_cumulative").
#'
#' @param weekly output of whoWeeklyToMonthly().
#' @param cumulative output of whoCumulativeToMonthly(), or NULL (weekly only).
#' @return tibble(iso3, t_start, cases, n_days, partial, source), sorted.
whoMergeMonthly <- function(weekly, cumulative = NULL) {
  weekly$source <- rep("dashboard_weekly", nrow(weekly))
  if (!is.null(cumulative) && nrow(cumulative)) {
    key_w <- paste(weekly$iso3, format(as.Date(weekly$t_start)))
    key_c <- paste(cumulative$iso3, format(as.Date(cumulative$t_start)))
    add <- cumulative[!key_c %in% key_w, , drop = FALSE]
    add$source <- rep("sitrep_cumulative", nrow(add))
    weekly <- dplyr::bind_rows(weekly, add)
  }
  dplyr::arrange(dplyr::as_tibble(weekly), iso3, t_start)
}

#' whoObsTable
#'
#' The scored WHO observations: complete months in `[start, end)`.
#'
#' @param who_monthly output of whoWeeklyToMonthly() (or its saved file).
#' @return tibble(country, t_start, cases_who).
whoObsTable <- function(who_monthly, start, end) {
  start <- as.Date(start); end <- as.Date(end)
  x <- who_monthly[!who_monthly$partial &
                     who_monthly$t_start >= start & who_monthly$t_start < end, , drop = FALSE]
  dplyr::tibble(country = as.character(x$iso3), t_start = as.Date(x$t_start),
                cases_who = as.numeric(x$cases))
}

#' whoCoverageFraction
#'
#' Share of each country's population covered by the simulated gids:
#' sum of the frame's `pop` over the country's gids / the WPP national total in
#' `pop_year`. The WHO count is national, the simulated count covers only the
#' gids in the frame, so the simulated count is scaled by 1 / who_f.
#'
#' Clamped to (0, 1]: gid populations and WPP totals come from different
#' products, so a fully covered country can exceed 1 by rounding. A country
#' missing from WPP is an error (no silent default).
#'
#' @param frame location frame with `country` (ISO3) and `pop`.
#' @param wpp   loadWppTable() table (`iso3`, `year`, `pop_total`).
#' @return tibble(country, pop_covered, pop_national, who_f).
whoCoverageFraction <- function(frame, wpp, pop_year) {
  if (is.null(wpp))
    stop("whoCoverageFraction: no WPP table (simulator WPP path missing); ",
         "the WHO term needs national populations.")
  if (!"pop" %in% names(frame) || anyNA(frame$pop))
    stop("whoCoverageFraction: frame needs a complete `pop` column.")
  cov <- dplyr::summarise(dplyr::group_by(dplyr::tibble(country = as.character(frame$country),
                                                         pop = as.numeric(frame$pop)), country),
                          pop_covered = sum(pop), .groups = "drop")
  w <- wpp[as.integer(wpp$year) == as.integer(pop_year), c("iso3", "pop_total")]
  cov$pop_national <- as.numeric(w$pop_total[match(cov$country, w$iso3)])
  bad <- cov$country[is.na(cov$pop_national) | cov$pop_national <= 0]
  if (length(bad))
    stop(sprintf("whoCoverageFraction: no WPP %d total for: %s",
                 as.integer(pop_year), paste(bad, collapse = ", ")))
  cov$who_f <- pmin(1, cov$pop_covered / cov$pop_national)
  if (any(cov$who_f <= 0))
    stop("whoCoverageFraction: zero covered population for: ",
         paste(cov$country[cov$who_f <= 0], collapse = ", "))
  cov
}

#' whoBlockShare
#'
#' Each unit's share of its country's WHO term, summing to 1 within the
#' country. With `block_granularity: country` the split does not matter (the
#' whole country sits in one bpfilter block, so the term enters it exactly
#' once). With ADM1 blocks (who_obs.allow_adm1_blocks) it is the tempering
#' weight of each block.
#'
#' Contribution weights (a unit's share of the simulated count) are refused:
#' a state-dependent multiplier on a log-pmf <= 0 favours particles with no
#' simulated cases over one that matches the WHO count.
#'
#' @param frame location frame (`country`, `pop`).
#' @param block_share "pop" (population share) or "equal".
#' @return numeric vector, frame row order.
whoBlockShare <- function(frame, block_share = c("pop", "equal")) {
  block_share <- match.arg(block_share)
  ct <- as.character(frame$country)
  w  <- if (block_share == "pop") as.numeric(frame$pop) else rep(1, nrow(frame))
  if (anyNA(w) || any(w < 0)) stop("whoBlockShare: invalid `pop` in frame.")
  tot <- stats::ave(w, ct, FUN = sum)
  # A country whose units all have pop 0 falls back to equal shares.
  zero <- tot <= 0
  w[zero] <- 1; tot[zero] <- stats::ave(w, ct, FUN = sum)[zero]
  w / tot
}

# ---------------------------------------------------------------------------
# Per-country WHO term scale
# ---------------------------------------------------------------------------
# The WHO term scores the simulated national count times a multiplier against
# the WHO dashboard count. Fixed mode: multiplier =
# who_obs.scale / who_f, i.e. the uncovered population is assumed to have the
# covered gids' incidence and WHO to report like the database. Estimated mode
# (who_obs.scale_mode: estimated): scale_c is calibrated per country on the
# months both sources observe, so multiplier = scale_c / who_f = W_c / D_c
# shrunk towards the pooled ratio (who_f cancels for an identified country).

#' whoDbObservedMonths
#'
#' The database side of the overlap, on the filter's own observation: per
#' (country, month) of the window, the observed monthly cases summed over the
#' frame's gids (`obs_panel` = smoothing_core's `surv_panel_obs`, whose
#' `cases_month` is already footprint-aligned per gid: sum within admin level,
#' max across levels), and whether the database observes the country that
#' month (at least one frame gid with a non-NA `cases_month`).
#' @param obs_panel tibble(gid, t_start, cases_month, ...).
#' @param frame location frame (gid, country) -- after the ADM1/ADM2 overlap drop.
#' @return tibble(country, t_start, db_cases, db_observed), one row per
#'   (country, month) with any frame-gid panel row in [start, end).
whoDbObservedMonths <- function(obs_panel, frame, start, end) {
  start <- as.Date(start); end <- as.Date(end)
  if (!"cases_month" %in% names(obs_panel))
    stop("whoDbObservedMonths: the observation panel has no `cases_month` column ",
         "(simulator.filter.monthly_case_obs is 'none'); the WHO scale needs the ",
         "database's monthly cases.", call. = FALSE)
  ctry <- stats::setNames(as.character(frame$country), as.character(frame$gid))
  x <- obs_panel[as.character(obs_panel$gid) %in% names(ctry) &
                   as.Date(obs_panel$t_start) >= start & as.Date(obs_panel$t_start) < end, , drop = FALSE]
  x <- dplyr::tibble(country = unname(ctry[as.character(x$gid)]), t_start = as.Date(x$t_start),
                     cases = as.numeric(x$cases_month))
  out <- dplyr::summarise(dplyr::group_by(x, country, t_start),
                          db_cases = sum(cases, na.rm = TRUE),
                          db_observed = any(!is.na(cases)), .groups = "drop")
  dplyr::arrange(out, country, t_start)
}

#' whoScaleEstimate
#'
#' Per-country WHO term scale. Over the months O_c in [start, end) that the
#' WHO table scores (complete months, whoObsTable()) AND the database observes
#' (whoDbObservedMonths()): W_c = WHO cases, D_c = database cases, E_c = D_c /
#' who_f_c (the database count extrapolated to the national population).
#' Pooled ratio over the identified countries s0 = sum W / sum E; per country
#' scale_c = (W_c + kappa s0) / (E_c + kappa) (kappa = `shrink_cases`), for a
#' country with D_c >= `min_db_cases` and |O_c| >= `min_months`; any other
#' country takes the fallback (`pooled`: s0; `fixed`: `fixed_scale`). Clamped
#' to `clamp`. multiplier = scale / who_f is what the filter applies.
#'
#' @param db_months whoDbObservedMonths() table.
#' @param who_tbl   whoObsTable() table (country, t_start, cases_who).
#' @param cov       whoCoverageFraction() table (country, pop_covered,
#'   pop_national, who_f): one row per frame country.
#' @return tibble, one row per frame country: country, who_f, pop_covered,
#'   pop_national, est_months, who_cases, db_cases, db_cases_ext, ratio_raw
#'   (W / D, the implied WHO/DB ratio), scale_raw (W / E), scale_pooled (s0),
#'   scale, multiplier, source (estimated | fallback_pooled | fallback_fixed),
#'   clamped, months_who_pos_db_zero, months_db_pos_who_zero (diagnostics
#'   over the overlap months).
whoScaleEstimate <- function(db_months, who_tbl, cov, start, end,
                             min_db_cases = 100, min_months = 6L, shrink_cases = 0,
                             fallback = c("pooled", "fixed"), fixed_scale = 1,
                             clamp = c(0.1, 10)) {
  fallback <- match.arg(fallback)
  start <- as.Date(start); end <- as.Date(end)
  clamp <- as.numeric(unlist(clamp))
  if (length(clamp) != 2L || !all(is.finite(clamp)) || clamp[1] <= 0 || clamp[1] >= clamp[2])
    stop("whoScaleEstimate: clamp must be c(lo, hi) with 0 < lo < hi.", call. = FALSE)
  w <- who_tbl[as.Date(who_tbl$t_start) >= start & as.Date(who_tbl$t_start) < end, , drop = FALSE]
  d <- db_months[db_months$db_observed %in% TRUE, , drop = FALSE]
  key_w <- paste(w$country, format(as.Date(w$t_start)))
  key_d <- paste(d$country, format(as.Date(d$t_start)))
  ov <- dplyr::tibble(country = as.character(w$country), t_start = as.Date(w$t_start),
                      W = as.numeric(w$cases_who),
                      D = as.numeric(d$db_cases)[match(key_w, key_d)])
  ov <- ov[!is.na(ov$D), , drop = FALSE]
  agg <- dplyr::summarise(dplyr::group_by(ov, country),
                          est_months = dplyr::n(), who_cases = sum(W), db_cases = sum(D),
                          months_who_pos_db_zero = sum(W > 0 & D == 0),
                          months_db_pos_who_zero = sum(D > 0 & W == 0), .groups = "drop")
  out <- dplyr::tibble(country = as.character(cov$country), who_f = as.numeric(cov$who_f),
                       pop_covered = as.numeric(cov$pop_covered),
                       pop_national = as.numeric(cov$pop_national))
  out <- dplyr::left_join(out, agg, by = "country")
  for (k in c("est_months", "months_who_pos_db_zero", "months_db_pos_who_zero"))
    out[[k]] <- as.integer(dplyr::coalesce(as.integer(out[[k]]), 0L))
  out$who_cases <- dplyr::coalesce(as.numeric(out$who_cases), 0)
  out$db_cases  <- dplyr::coalesce(as.numeric(out$db_cases), 0)
  out$db_cases_ext <- out$db_cases / out$who_f
  out$ratio_raw <- ifelse(out$db_cases > 0, out$who_cases / out$db_cases, NA_real_)
  out$scale_raw <- ifelse(out$db_cases_ext > 0, out$who_cases / out$db_cases_ext, NA_real_)
  ident <- out$db_cases >= min_db_cases & out$est_months >= as.integer(min_months)
  s0 <- if (any(ident)) sum(out$who_cases[ident]) / sum(out$db_cases_ext[ident]) else NA_real_
  if (identical(fallback, "pooled") && !is.finite(s0))
    stop(sprintf(paste0("whoScaleEstimate: no country identified (db_cases >= %s and >= %d ",
                        "overlap months in [%s, %s)); a pooled scale needs at least one. Use ",
                        "fallback: fixed or widen scale_est."),
                 format(min_db_cases), as.integer(min_months), format(start), format(end)),
         call. = FALSE)
  out$scale_pooled <- s0
  kappa <- as.numeric(shrink_cases)
  est <- if (is.finite(s0)) (out$who_cases + kappa * s0) / (out$db_cases_ext + kappa)
         else out$scale_raw
  out$scale <- ifelse(ident, est, if (identical(fallback, "pooled")) s0 else as.numeric(fixed_scale))
  out$source <- ifelse(ident, "estimated", paste0("fallback_", fallback))
  out$clamped <- out$scale < clamp[1] | out$scale > clamp[2]
  out$scale <- pmin(pmax(out$scale, clamp[1]), clamp[2])
  out$multiplier <- whoMultiplier(out$who_f, out$scale)
  out[, c("country", "who_f", "pop_covered", "pop_national", "est_months", "who_cases",
          "db_cases", "db_cases_ext", "ratio_raw", "scale_raw", "scale_pooled", "scale",
          "multiplier", "source", "clamped", "months_who_pos_db_zero", "months_db_pos_who_zero")]
}

#' whoMultiplier
#'
#' The one definition of the WHO term's national multiplier: the simulated
#' covered-gid count times scale / who_f is compared with the WHO national
#' count. `scale` is the fixed who_obs.scale (a scalar, the fixed-mode term)
#' or the per-country estimated scale (vector, aligned with `who_f`). Used by
#' the C filter builder, who_vs_smoothing() (08c) and 09.
whoMultiplier <- function(who_f, scale = 1) scale / who_f

#' whoDropDbOverlap
#'
#' db_overlap: skip_who -- drop the WHO rows of every (country, month) the
#' database observes, so no country-month is scored by both sources (the
#' estimated scale is calibrated on exactly those months).
#' @return `who_tbl` without those rows.
whoDropDbOverlap <- function(who_tbl, db_months) {
  obs <- db_months[db_months$db_observed %in% TRUE, , drop = FALSE]
  drop <- paste(who_tbl$country, format(as.Date(who_tbl$t_start))) %in%
    paste(obs$country, format(as.Date(obs$t_start)))
  who_tbl[!drop, , drop = FALSE]
}

#' whoScaleDigest
#'
#' Content digest (base.R::content_digest(), call-time) of a
#' whoScaleEstimate() table, rows sorted by country: the identity 08a writes
#' into each task, 08b re-checks and 08c records in the smoothing cache
#' meta.yml.
whoScaleDigest <- function(tbl) {
  tbl <- as.data.frame(tbl)[order(as.character(tbl$country), method = "radix"), , drop = FALSE]
  rownames(tbl) <- NULL
  content_digest(tbl)
}
