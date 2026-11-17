# loo_lags.R — the shared ENSO single-lag selection rule, ranking-file
#   lookup, Pareto-k summaries and per-lag LOO sidecars.
#
# Provides:
#   select_best_lag()            dispatch over the rules and tracks
#   select_best_lag_occurrence(), select_best_lag_downstream(),
#   select_best_lag_best_elpd()  the rules (lag-rule note below)
#   lag_selection_cfg(), check_occ_floor_offset()   deploy.lag_selection
#   best_lags_by_group()         selected lag per (track, group)
#   latest_ranking_file(), read_latest_ranking()    loo_compare_lags_*.rds
#   pareto_k_summary(), ranking_pareto()            Pareto-k of each lag
#   build_lag_loo_sidecar(), read_lag_loo_sidecar() per-lag LOO sidecar
#
# Details:
#   One selection rule for 05_fit_postprocess_loo_lags.R, deploy_core.R
#   (deploy_select_lags()) and notebooks/loo_lag_results_by_group.Rmd, and one
#   ranking-file "latest" rule for every reader of the ranking caches.

# Lag-rule note (the roxygen below refers to it).
#
# Selection rule (`ranking`: one row per single_lag, from loo_compare_lags_*.rds
# `$ranking`; elpd_diff / se_diff relative to the rank-1 fit):
#   occurrence: the MAXIMUM single_lag within `se_mult` SEs (se_diff) of the
#     rank-1 fit.
#   size / duration / monthly: hybrid with the group's occurrence lag occ_lag:
#     1. the track's own occurrence-style selection;
#     2. if that lag >= occ_lag, keep it;
#     3. else the MAX single_lag >= occ_lag within `se_mult` SEs of the best fit
#        in that subset (select_best_lag_within_eligible(), approximated since
#        the cached diffs are relative to the possibly excluded rank-1).
#     `enforce_occ_floor = FALSE` returns step 1 only.
#   leads: a negative single_lag is an ENSO LEAD. Every rule takes `min_lag`
#     (default 0). The SE window stays anchored on the rank-1 fit (cached paired
#     SEs), then is filtered to eligible lags; if none is in the window, the
#     best-ranked eligible lag is taken. Do NOT use
#     select_best_lag_within_eligible() here: its total-se quadrature admits
#     nearly every lag.
#   best_elpd (`rule = "best_elpd"`): occurrence takes its best-ranked lag
#     (min_lag defaults to -Inf); the other tracks the best-ranked lag
#     >= occ_lag - occ_floor_offset (no floor when NULL). No SE window. See
#     select_best_lag_best_elpd() and lag_selection_cfg().

#' Select the maximum `single_lag` among fits within `se_mult` SEs of
#' `elpd_diff` from the rank-1 fit — the self-contained "occurrence" rule,
#' also used as step 1 of the downstream-track hybrid rule.
#'
#' @param ranking A ranking data.frame/tibble (see the lag-rule note for the
#'   expected shape). Must contain exactly one `rank == 1` row.
#' @param se_mult Number of SEs of `se_diff` defining "statistically
#'   indistinguishable from the best fit" [default 2].
#' @param min_lag Smallest eligible single_lag [default 0: leads (< 0) are
#'   ranked but never selected; see the lag-rule note].
#' @return A list: `lag` (selected single_lag, integer), `rank` (its rank in
#'   `ranking`), `within_se` (the eligible subset of `ranking` this was
#'   chosen from, for provenance/reporting).
select_best_lag_occurrence <- function(ranking, se_mult = 2, min_lag = 0L) {
  stopifnot(
    is.data.frame(ranking),
    all(c("rank", "single_lag", "elpd_diff", "se_diff") %in% names(ranking))
  )
  if (sum(ranking$rank == 1) != 1L) {
    stop(
      "select_best_lag_occurrence(): `ranking` must have exactly one rank == 1 row."
    )
  }
  elig <- ranking[ranking$single_lag >= min_lag, , drop = FALSE]
  if (nrow(elig) == 0) {
    stop(sprintf(
      "select_best_lag_occurrence(): no single_lag >= min_lag (%s) in `ranking`.",
      format(min_lag)
    ))
  }
  within <- elig[abs(elig$elpd_diff) <= se_mult * elig$se_diff, , drop = FALSE]
  # No eligible lag within the window of the (lead) rank-1 fit: best-ranked eligible lag.
  if (nrow(within) == 0) {
    within <- elig[which.min(elig$rank), , drop = FALSE]
  }
  sel_lag <- max(within$single_lag)
  sel_row <- within[within$single_lag == sel_lag, , drop = FALSE][1, ]
  list(lag = as.integer(sel_lag), rank = sel_row$rank, within_se = within)
}

#' Select the maximum `single_lag` among fits with `single_lag >= occ_lag`
#' that are within `se_mult` SEs of the best-ranked fit *within that eligible
#' subset*. Used as the fallback step of the downstream-track hybrid rule,
#' when the track's unconstrained pick (select_best_lag_occurrence()) falls
#' below `occ_lag`.
#'
#' The cached `elpd_diff`/`se_diff` columns are exact only relative to the
#' ranking's actual rank-1 fit (loo::loo_compare()'s convention) — which the
#' `occ_lag` floor may exclude, as it does whenever this fallback triggers.
#' There is no cached pairwise SE between two arbitrary (non-rank-1) fits, so
#' the "within SE of the eligible-best fit" test here is approximated from
#' each fit's own marginal `se_elpd` combined in quadrature
#' (`sqrt(se_elpd_i^2 + se_elpd_ref^2)`, an independence assumption). This is
#' conservative (wider) relative to a true paired-LOO SE for positively
#' correlated fits, which — for this "prefer the longest indistinguishable
#' lag" rule — errs toward including more candidate lags rather than fewer.
#'
#' @param ranking A ranking data.frame/tibble (see the lag-rule note for shape).
#' @param occ_lag Integer; the occurrence lag selected for the same group.
#' @param se_mult Number of (approximated) SEs defining "statistically
#'   indistinguishable from the best eligible fit" [default 2].
#' @return A list: `lag` (selected single_lag, integer), `rank` (its rank),
#'   `within_se` (the eligible subset this was chosen from), `ref` (the
#'   best-ranked eligible row used as the SE-window anchor).
select_best_lag_within_eligible <- function(ranking, occ_lag, se_mult = 2) {
  stopifnot(
    is.data.frame(ranking),
    all(c("rank", "single_lag", "elpd_loo", "se_elpd") %in% names(ranking))
  )
  if (is.null(occ_lag) || is.na(occ_lag)) {
    stop("select_best_lag_within_eligible(): `occ_lag` is required.")
  }
  elig <- ranking[ranking$single_lag >= occ_lag, , drop = FALSE]
  if (nrow(elig) == 0) {
    stop(sprintf(
      "select_best_lag_within_eligible(): no single_lag >= occ_lag (%d) in `ranking`.",
      occ_lag
    ))
  }
  ref <- elig[which.min(elig$rank), , drop = FALSE][1, ]
  d_elpd <- elig$elpd_loo - ref$elpd_loo
  d_se <- sqrt(elig$se_elpd^2 + ref$se_elpd^2)
  within <- elig[abs(d_elpd) <= se_mult * d_se, , drop = FALSE]
  sel_lag <- max(within$single_lag)
  sel_row <- within[within$single_lag == sel_lag, , drop = FALSE][1, ]
  list(
    lag = as.integer(sel_lag),
    rank = sel_row$rank,
    within_se = within,
    ref = ref
  )
}

#' Select a downstream-track (size / duration / monthly) lag via the hybrid
#' rule described in the lag-rule note: the track's own unconstrained
#' SE-window pick if it already clears `occ_lag`, otherwise the SE-window
#' pick within the `single_lag >= occ_lag` eligible subset.
#'
#' @param ranking A ranking data.frame/tibble (see the lag-rule note for shape).
#' @param occ_lag Integer; the occurrence lag selected for the same group
#'   (the floor below which downstream lags are not eligible, when the
#'   unconstrained pick doesn't already clear it). Required when
#'   `enforce_occ_floor = TRUE`; ignored (may be `NULL`) otherwise.
#' @param enforce_occ_floor Logical; when `TRUE` (default), apply the hybrid
#'   rule above. When `FALSE`, ignore occurrence entirely and return the
#'   track's unconstrained `select_best_lag_occurrence()` pick.
#' @param se_mult SE multiplier used at every step (unconstrained pick and,
#'   if triggered, the eligible-subset fallback) [default 2].
#' @param min_lag Smallest eligible single_lag [default 0; see the lag-rule note].
#' @return A list: `lag` (selected single_lag, integer), `rank` (its rank).
select_best_lag_downstream <- function(
  ranking,
  occ_lag = NULL,
  enforce_occ_floor = TRUE,
  se_mult = 2,
  min_lag = 0L
) {
  stopifnot(
    is.data.frame(ranking),
    all(c("rank", "single_lag") %in% names(ranking))
  )
  if (!isTRUE(enforce_occ_floor)) {
    return(select_best_lag_occurrence(
      ranking,
      se_mult = se_mult,
      min_lag = min_lag
    ))
  }
  if (is.null(occ_lag) || is.na(occ_lag)) {
    stop(
      "select_best_lag_downstream(): `occ_lag` is required when enforce_occ_floor = TRUE."
    )
  }
  sel_unconstrained <- select_best_lag_occurrence(
    ranking,
    se_mult = se_mult,
    min_lag = min_lag
  )
  if (sel_unconstrained$lag >= occ_lag) {
    return(sel_unconstrained)
  }
  select_best_lag_within_eligible(
    ranking,
    occ_lag = max(occ_lag, min_lag),
    se_mult = se_mult
  )
}

#' Best-ranked eligible lag (rule "best_elpd").
#'
#' Eligible: `single_lag >= min_lag` and, when both `occ_lag` and
#' `occ_floor_offset` are given, `single_lag >= occ_lag - occ_floor_offset`
#' (downstream tracks: their ENSO signal may precede the occurrence signal by
#' at most one outbreak duration). Among the eligible lags, the lowest
#' `rank` (highest elpd_loo) wins. Leads (negative single_lag) are eligible
#' unless `min_lag` excludes them.
#'
#' @param ranking ranking data.frame (see the lag-rule note).
#' @param min_lag smallest eligible single_lag [default -Inf: whole grid].
#' @param occ_lag occurrence lag for the same group (downstream tracks), or NULL.
#' @param occ_floor_offset months the downstream lag may fall below occ_lag,
#'   or NULL for no occurrence floor.
#' @return list(lag, rank, within_se = the eligible subset).
select_best_lag_best_elpd <- function(
  ranking,
  min_lag = -Inf,
  occ_lag = NULL,
  occ_floor_offset = NULL
) {
  stopifnot(
    is.data.frame(ranking),
    all(c("rank", "single_lag") %in% names(ranking))
  )
  floor <- min_lag
  if (!is.null(occ_lag) && !is.na(occ_lag) && !is.null(occ_floor_offset)) {
    floor <- max(floor, occ_lag - as.numeric(occ_floor_offset))
  }
  elig <- ranking[ranking$single_lag >= floor, , drop = FALSE]
  if (nrow(elig) == 0) {
    stop(sprintf(
      "select_best_lag_best_elpd(): no single_lag >= %s in `ranking`.",
      format(floor)
    ))
  }
  best <- elig[which.min(elig$rank), , drop = FALSE]
  list(lag = as.integer(best$single_lag), rank = best$rank, within_se = elig)
}

#' Selection rules (deploy.lag_selection.rule).
LAG_SELECTION_RULES <- c("longest_within_se", "best_elpd")

#' lag_selection_cfg -- normalise a deploy.lag_selection block.
#'
#' Missing leaves, min_lag "rule_default" and occ_floor_offset "none" fall back
#' to the current rule's defaults, so a deploy
#' without the block selects exactly as before: rule longest_within_se,
#' se_mult 2, min_lag 0 (leads ranked, never selected), no occurrence floor.
#' For rule best_elpd, min_lag defaults to -Inf (leads eligible).
#' @return list(rule, se_mult, min_lag, occ_floor_offset).
lag_selection_cfg <- function(x = NULL) {
  x <- x %||% list()
  rule <- x[["rule"]] %||% "longest_within_se"
  if (!rule %in% LAG_SELECTION_RULES) {
    stop(
      "lag_selection.rule must be one of ",
      paste(LAG_SELECTION_RULES, collapse = " | "),
      " (got '",
      rule,
      "').",
      call. = FALSE
    )
  }
  ml <- x[["min_lag"]]
  if (identical(ml, "rule_default")) {
    ml <- NULL
  }
  off <- x[["occ_floor_offset"]]
  if (identical(off, "none")) {
    off <- NULL
  }
  if (!is.null(ml) && (!is.numeric(ml) || length(ml) != 1L || is.na(ml))) {
    stop(
      "lag_selection.min_lag must be a number or 'rule_default'.",
      call. = FALSE
    )
  }
  if (
    !is.null(off) &&
      (!is.numeric(off) || length(off) != 1L || is.na(off) || off < 0)
  ) {
    stop(
      "lag_selection.occ_floor_offset must be a non-negative number or 'none'.",
      call. = FALSE
    )
  }
  list(
    rule = rule,
    se_mult = as.numeric(x[["se_mult"]] %||% 2),
    min_lag = if (is.null(ml)) {
      (if (identical(rule, "best_elpd")) -Inf else 0L)
    } else {
      as.numeric(ml)
    },
    occ_floor_offset = off
  )
}

#' check_occ_floor_offset -- the deploy lag floor against the episode span.
#'
#' `occ_floor_offset` (months a downstream lag may fall below the occurrence
#' lag) is meant to be the median whole-month episode span: onset month to
#' last active month, `n_months - 1` (duration_span_summary() in
#' dur_data.R, stamped by 03_prepare_dur.R as `median_span_months`). Warns,
#' never stops: the floor is a deploy decision and may differ on purpose.
#'
#' @param cfg normalised lag selection, lag_selection_cfg()'s return value.
#' @param span the duration panel's policy stamp (read_prep_policy_stamp() on the
#'   duration panel), or NULL.
#' @return invisible status: "no_floor" (occ_floor_offset unset; nothing to
#'   check), "match" (silent), "mismatch" (warning()) or "missing" (message():
#'   no stamp, or a stamp without `median_span_months`, e.g. one written
#'   before the field existed: re-run 03_prepare_dur.R --restamp).
check_occ_floor_offset <- function(cfg, span) {
  off <- cfg[["occ_floor_offset"]]
  if (is.null(off)) {
    return(invisible("no_floor"))
  }
  med <- (span %||% list())[["median_span_months"]]
  if (is.null(med) || length(med) != 1L || is.na(med)) {
    message(
      "[lag selection] occ_floor_offset = ",
      format(off),
      ": the duration ",
      "policy stamp has no median_span_months to check it against ",
      "(re-run analysis/03_prepare_dur.R --restamp)."
    )
    return(invisible("missing"))
  }
  if (isTRUE(all.equal(as.numeric(off), as.numeric(med)))) {
    return(invisible("match"))
  }
  warning(
    sprintf(
      paste0(
        "[lag selection] occ_floor_offset = %s differs from the median whole-month ",
        "episode span %s (duration policy stamp, median_span_months; IQR %s-%s). ",
        "The floor is meant to be that span; keep it only if the difference is deliberate."
      ),
      format(off),
      format(med),
      format(span[["span_months_q25"]] %||% NA),
      format(span[["span_months_q75"]] %||% NA)
    ),
    call. = FALSE
  )
  invisible("mismatch")
}

#' Dispatch to `select_best_lag_occurrence()` / `select_best_lag_downstream()`
#' by track name — the single entry point pipeline callers should use.
#'
#' @param ranking A ranking data.frame/tibble (see the lag-rule note for shape).
#' @param track One of "occurrence", "size", "duration", "monthly".
#' @param occ_lag Integer; required (and used) for every track except
#'   "occurrence" when `enforce_occ_floor = TRUE`; ignored otherwise.
#' @param se_mult SE multiplier used throughout (both the occurrence rule and
#'   every step of the downstream hybrid rule) [default 2].
#' @param enforce_occ_floor Logical [default TRUE]; passed through to
#'   `select_best_lag_downstream()` for non-occurrence tracks (ignored for
#'   "occurrence"). Set `FALSE` to select a downstream track's unconstrained
#'   SE-window lag, independent of the occurrence lag.
#' @param min_lag Smallest eligible single_lag [default 0; see the lag-rule note].
#' @return A list with at least `lag` and `rank` (see the functions above).
select_best_lag <- function(
  ranking,
  track,
  occ_lag = NULL,
  se_mult = 2,
  enforce_occ_floor = TRUE,
  min_lag = 0L,
  rule = "longest_within_se",
  occ_floor_offset = NULL
) {
  # rule best_elpd (see select_best_lag_best_elpd()): occurrence on the whole
  # eligible grid, downstream floored at occ_lag - occ_floor_offset. The
  # longest_within_se arguments (se_mult, enforce_occ_floor) do not apply.
  if (identical(rule, "best_elpd")) {
    occ <- identical(track, "occurrence") || track %in% c("occ")
    return(select_best_lag_best_elpd(
      ranking,
      min_lag = min_lag,
      occ_lag = if (occ) NULL else occ_lag,
      occ_floor_offset = if (occ) NULL else occ_floor_offset
    ))
  }
  if (!identical(rule, "longest_within_se")) {
    stop("select_best_lag(): unknown rule '", rule, "'.", call. = FALSE)
  }
  if (identical(track, "occurrence")) {
    select_best_lag_occurrence(ranking, se_mult = se_mult, min_lag = min_lag)
  } else {
    select_best_lag_downstream(
      ranking,
      occ_lag = occ_lag,
      enforce_occ_floor = enforce_occ_floor,
      se_mult = se_mult,
      min_lag = min_lag
    )
  }
}

# ---------------------------------------------------------------------------
# Ranking-file location
#
# Where 05_fit_postprocess_loo_lags.R writes loo_compare_lags_*.rds: one
# location (make_loo_compare_run_file(): <outdir>/loo_compare/<batch>/, legacy
# analysis/configs/<batch>/ read-only fallback) and one "latest" rule (filename
# YYYYMMDD descending, then mtime; mtime alone resets on checkout or cp).
# Callers must source paths.R and base.R first (modules never source each
# other).
# ---------------------------------------------------------------------------

#' Filename shape of a ranking cache: `loo_compare_lags[_<group>]_<YYYYMMDD>.rds`.
#' Capture 3 = group token (empty for the ungrouped file), capture 4 = date.
LOO_COMPARE_RANKING_RE <- "^loo_compare_lags(_(.+))?_([0-9]{8})\\.rds$"

# Legacy (in-config-dir) location: `batch` itself when it is a directory,
# else analysis/configs/<batch> for a bare batch name.
.loo_compare_legacy_dir <- function(batch, legacy_dir = NULL) {
  if (!is.null(legacy_dir)) {
    return(legacy_dir)
  }
  if (dir.exists(batch)) batch else file.path("analysis/configs", batch)
}

.loo_compare_empty <- function() {
  data.frame(
    file = character(0),
    group = character(0),
    date = character(0),
    mtime = as.POSIXct(character(0)),
    legacy = logical(0),
    stringsAsFactors = FALSE
  )
}

# All ranking caches in one directory, parsed; NULL when none.
.loo_compare_list <- function(dir, legacy) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) {
    return(NULL)
  }
  files <- list.files(dir, pattern = LOO_COMPARE_RANKING_RE, full.names = TRUE)
  if (!length(files)) {
    return(NULL)
  }
  m <- regmatches(
    basename(files),
    regexec(LOO_COMPARE_RANKING_RE, basename(files))
  )
  data.frame(
    file = files,
    group = vapply(
      m,
      function(x) if (nzchar(x[3])) x[3] else NA_character_,
      character(1)
    ),
    date = vapply(m, `[`, character(1), 4L),
    mtime = file.mtime(files),
    legacy = rep(legacy, length(files)),
    stringsAsFactors = FALSE
  )
}

#' Ranking caches for a batch: new location first, legacy location as fallback.
#'
#' @param batch Batch config dir (`analysis/configs/<batch>/`) or bare batch
#'   name (`"<batch>"`); both resolve through `make_loo_compare_run_file()`.
#' @param group NULL = every group; "" = the ungrouped file only; a group
#'   token (e.g. "enso_east") = that group only. Group tokens are matched
#'   verbatim (with underscores) -- ranking filenames carry the group string
#'   as written by `05_fit_postprocess_loo_lags.R --group`, unlike the YAML
#'   config filenames which strip non-alphanumerics (`_grp-ensoeast`).
#' @param opt Options list; only `opt$outdir` is used (via make_loo_compare_run_file).
#' @param legacy_dir Override the fallback directory [default: `batch` if it
#'   is a directory, else `analysis/configs/<batch>`].
#' @param latest_only TRUE = one row per group (the latest by date, then
#'   mtime); FALSE = every file, e.g. for a freshness-input list.
#' @return data.frame(file, group [NA = ungrouped], date [chr YYYYMMDD],
#'   mtime, legacy [lgl: came from the fallback dir]) ordered date desc then
#'   mtime desc; zero rows if nothing matched in either location.
loo_compare_ranking_files <- function(
  batch,
  group = NULL,
  opt = list(),
  legacy_dir = NULL,
  latest_only = TRUE
) {
  new_dir <- make_loo_compare_run_file(opt, "loo_compare_dir", batch = batch)
  tbl <- .loo_compare_list(new_dir, legacy = FALSE)
  if (is.null(tbl)) {
    tbl <- .loo_compare_list(
      .loo_compare_legacy_dir(batch, legacy_dir),
      legacy = TRUE
    )
  }
  if (is.null(tbl)) {
    return(.loo_compare_empty())
  }
  if (!is.null(group)) {
    tbl <- if (nzchar(group)) {
      tbl[!is.na(tbl$group) & tbl$group == group, , drop = FALSE]
    } else {
      tbl[is.na(tbl$group), , drop = FALSE]
    }
  }
  tbl <- tbl[order(-xtfrm(tbl$date), -as.numeric(tbl$mtime)), , drop = FALSE]
  if (latest_only && nrow(tbl)) {
    key <- ifelse(is.na(tbl$group), "", tbl$group)
    tbl <- tbl[!duplicated(key), , drop = FALSE]
  }
  rownames(tbl) <- NULL
  tbl
}

#' The single latest ranking cache for `group` (or, with `group = NULL`, the
#' latest across every group), or NULL when none exists in either location.
latest_ranking_file <- function(
  batch,
  group = NULL,
  opt = list(),
  legacy_dir = NULL
) {
  tbl <- loo_compare_ranking_files(
    batch,
    group = group,
    opt = opt,
    legacy_dir = legacy_dir,
    latest_only = TRUE
  )
  if (!nrow(tbl)) {
    return(NULL)
  }
  tbl$file[1L]
}

#' `readRDS(latest_ranking_file(...))$ranking` -- the ranking tibble whose
#' shape is documented in this file's header. `required = TRUE` stop()s when
#' no cache exists; `required = FALSE` returns NULL so callers can warn/skip.
read_latest_ranking <- function(
  batch,
  group = NULL,
  opt = list(),
  legacy_dir = NULL,
  required = TRUE
) {
  f <- latest_ranking_file(
    batch,
    group = group,
    opt = opt,
    legacy_dir = legacy_dir
  )
  if (is.null(f)) {
    if (required) {
      stop(
        sprintf(
          "[loo_lags] no loo_compare_lags ranking for group '%s' in %s (nor legacy %s)",
          group %||% "<any>",
          make_loo_compare_run_file(opt, "loo_compare_dir", batch = batch),
          .loo_compare_legacy_dir(batch, legacy_dir)
        ),
        call. = FALSE
      )
    }
    return(NULL)
  }
  readRDS(f)$ranking
}

#' Best single lag per ENSO group for one track: the latest ranking cache per
#' group (read_latest_ranking()) through the canonical select_best_lag() rule.
#'
#' Single owner of "which lag does a figure/driver pick for a group", so the
#' figure scripts agree with the deployed fits (deploy_core.R) and the
#' manuscript (09b_postprocess_results.R) instead of re-reading the raw rank-1
#' row of the ranking (`ranking$single_lag[1]`), which can differ from the
#' SE-window pick.
#'
#' @param batch Batch config dir or bare batch name (see loo_compare_ranking_files()).
#' @param groups Character vector of group tokens (e.g. "enso_east").
#' @param track "occurrence" | "size" | "duration" | "monthly" (select_best_lag()'s dispatch).
#' @param occ_lags Optional named integer vector (names = groups) of occurrence
#'   lags; only used by a downstream track with `enforce_occ_floor = TRUE`.
#' @param se_mult SE multiplier of the selection rule [default 2; 0 reduces the
#'   rule to the raw rank-1 lag].
#' @param enforce_occ_floor Passed to select_best_lag() [default FALSE: each
#'   track's own unconstrained pick, as deploy_select_lags()].
#' @param overrides Named integer vector (names = groups); wins over the ranking
#'   and needs no ranking file.
#' @param opt Options list; only `opt$outdir` is used.
#' @param min_lag Smallest eligible single_lag [default 0; see the lag-rule note].
#' @return Named integer vector, one entry per group that resolved. Groups with
#'   neither an override nor a ranking cache are dropped with a warning.
best_lags_by_group <- function(
  batch,
  groups,
  track,
  occ_lags = NULL,
  se_mult = 2,
  enforce_occ_floor = FALSE,
  overrides = integer(0),
  opt = list(),
  min_lag = 0L,
  rule = "longest_within_se",
  occ_floor_offset = NULL
) {
  out <- integer(0)
  for (g in groups) {
    if (g %in% names(overrides)) {
      out[g] <- as.integer(overrides[[g]])
      next
    }
    rk <- read_latest_ranking(batch, group = g, opt = opt, required = FALSE)
    if (is.null(rk)) {
      warning(
        sprintf(
          "[loo_lags] no ranking cache for group '%s' in batch '%s'; skipped",
          g,
          batch
        ),
        call. = FALSE
      )
      next
    }
    occ_lag <- if (!is.null(occ_lags) && g %in% names(occ_lags)) {
      occ_lags[[g]]
    } else {
      NULL
    }
    out[g] <- select_best_lag(
      rk,
      track = track,
      occ_lag = occ_lag,
      se_mult = se_mult,
      enforce_occ_floor = enforce_occ_floor,
      min_lag = min_lag,
      rule = rule,
      occ_floor_offset = occ_floor_offset
    )$lag
  }
  out
}

# ---------------------------------------------------------------------------
# Pareto-k diagnostics of the per-lag PSIS-LOO
#
# Every single-lag fit's psis_loo object (cached by 05_fit_postprocess_loo_lags.R
# via make_run_file(opt, "<kpfx>loo")) holds its Pareto-k values, but the ranking
# caches kept none of them. pareto_k_summary() reduces one loo object to a
# one-row table; ranking_pareto() appends those columns to a ranking (passing
# them through when the ranking already has them, else computing them from the
# per-fit loo objects, else NA), recording where they came from in
# attr(, "pareto_source"). The columns are additive: every selection rule
# above reads only its named columns (rank, single_lag, elpd_loo, se_elpd,
# elpd_diff, se_diff), so the picks are identical with or without them.
# ---------------------------------------------------------------------------

#' Pareto-k columns appended to a ranking (in this order).
PARETO_K_COLS <- c(
  "n_obs",
  "k_threshold",
  "n_k_bad",
  "n_k_gt07",
  "n_k_gt1",
  "k_max",
  "share_k_bad"
)

#' Pareto-k threshold used when the installed loo exports no threshold function.
PARETO_K_DEFAULT_THRESHOLD <- 0.7

# loo::pareto_k_threshold() where the installed loo exports it, else 0.7. (loo
# 2.10 keeps its sample-size-dependent threshold internal; the fixed 0.7 keeps
# the summary identical across the local and the cluster loo versions.)
.pareto_k_threshold <- function(loo_obj) {
  if ("pareto_k_threshold" %in% getNamespaceExports("loo")) {
    thr <- tryCatch(
      getExportedValue("loo", "pareto_k_threshold")(loo_obj),
      error = function(e) NULL
    )
    if (is.numeric(thr) && length(thr) >= 1L && is.finite(thr[1L])) {
      return(as.numeric(thr[1L]))
    }
  }
  PARETO_K_DEFAULT_THRESHOLD
}

# The Pareto-k vector of a psis_loo object, or of a synthetic list carrying
# `$diagnostics$pareto_k` (or a bare `$pareto_k`).
.pareto_k_values <- function(loo_obj) {
  k <- loo_obj$diagnostics$pareto_k %||% loo_obj$pareto_k
  if (is.null(k)) {
    stop(
      "pareto_k_summary(): no diagnostics$pareto_k in the loo object.",
      call. = FALSE
    )
  }
  as.numeric(k)
}

#' One-row Pareto-k summary of one PSIS-LOO object.
#'
#' @param loo A `psis_loo` object (or a list with `$diagnostics$pareto_k`).
#' @return data.frame with PARETO_K_COLS: n_obs, k_threshold (loo's exported
#'   threshold where available, else 0.7), n_k_bad (k > k_threshold), n_k_gt07
#'   (k > 0.7), n_k_gt1 (k > 1), k_max, share_k_bad (n_k_bad / n_obs). An NA k
#'   counts in n_obs but in none of the k counts.
pareto_k_summary <- function(loo) {
  k <- .pareto_k_values(loo)
  thr <- .pareto_k_threshold(loo)
  n <- length(k)
  n_bad <- sum(k > thr, na.rm = TRUE)
  data.frame(
    n_obs = as.integer(n),
    k_threshold = thr,
    n_k_bad = as.integer(n_bad),
    n_k_gt07 = as.integer(sum(k > 0.7, na.rm = TRUE)),
    n_k_gt1 = as.integer(sum(k > 1, na.rm = TRUE)),
    k_max = if (any(!is.na(k))) max(k, na.rm = TRUE) else NA_real_,
    share_k_bad = if (n > 0L) n_bad / n else NA_real_,
    stringsAsFactors = FALSE
  )
}

.pareto_na_row <- function() {
  data.frame(
    n_obs = NA_integer_,
    k_threshold = NA_real_,
    n_k_bad = NA_integer_,
    n_k_gt07 = NA_integer_,
    n_k_gt1 = NA_integer_,
    k_max = NA_real_,
    share_k_bad = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Append the Pareto-k columns to a lag ranking.
#'
#' Three sources, in order:
#'   1. the ranking already carries every PARETO_K_COLS column (a ranking
#'      written with them): returned unchanged, pareto_source "ranking";
#'   2. `loo_files` is given: each lag's pareto_k_summary() from its per-fit loo
#'      object, pareto_source "loo_cache" (or "loo_cache_partial" when some
#'      lags' files are missing or unreadable; those rows get NA);
#'   3. otherwise NA columns, pareto_source "none".
#'
#' @param ranking Ranking data.frame (lag-rule note shape; needs `single_lag`).
#' @param loo_files NULL, or a vector / list named by single_lag (as character,
#'   e.g. "-3", "8") of per-fit loo cache paths, or of loo objects already in
#'   memory (what 05 passes).
#' @return `ranking` with PARETO_K_COLS appended (a partial set is replaced),
#'   rows and every other column unchanged, and attr(, "pareto_source").
ranking_pareto <- function(ranking, loo_files = NULL) {
  stopifnot(is.data.frame(ranking), "single_lag" %in% names(ranking))
  if (all(PARETO_K_COLS %in% names(ranking))) {
    attr(ranking, "pareto_source") <- "ranking"
    return(ranking)
  }
  out <- ranking[, setdiff(names(ranking), PARETO_K_COLS), drop = FALSE]
  keys <- as.character(as.integer(out$single_lag))
  rows <- vector("list", length(keys))
  n_ok <- 0L
  for (i in seq_along(keys)) {
    src <- if (!is.null(loo_files) && keys[i] %in% names(loo_files)) {
      loo_files[[keys[i]]]
    } else {
      NULL
    }
    obj <- if (is.character(src) && length(src) == 1L && !is.na(src)) {
      if (file.exists(src)) {
        tryCatch(readRDS(src), error = function(e) NULL)
      } else {
        NULL
      }
    } else {
      src
    }
    row <- if (is.null(obj)) {
      NULL
    } else {
      tryCatch(pareto_k_summary(obj), error = function(e) NULL)
    }
    if (!is.null(row)) {
      n_ok <- n_ok + 1L
    }
    rows[[i]] <- row %||% .pareto_na_row()
  }
  pk <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    .pareto_na_row()[0, , drop = FALSE]
  }
  for (cc in PARETO_K_COLS) {
    out[[cc]] <- pk[[cc]]
  }
  attr(out, "pareto_source") <- if (is.null(loo_files) || n_ok == 0L) {
    "none"
  } else if (n_ok == length(keys)) {
    "loo_cache"
  } else {
    "loo_cache_partial"
  }
  out
}

# ---------------------------------------------------------------------------
# Per-(batch, group) lag sidecar (05_fit_postprocess_loo_lags.R --lag_sidecar)
#
# lag_loo_sidecar[_<group>]_<YYYYMMDD>.rds under <outdir>/loo_compare/<batch>/;
# its name never matches LOO_COMPARE_RANKING_RE, so no ranking reader sees it.
# Holds what lag averaging / stacking and per-lag ERFs need without a fit:
#   elpd_pointwise, pareto_k  N x L matrices (rows = data rows, cols = lags)
#   pareto, estimates         per-lag pareto_k_summary() / elpd_loo, se, p_loo
#   cb                        per lag: b_draws (b_enso_cb*), cb_meta, sources
#   row_digest                digest of the row key, checked on write and read
#   track, group, batch       explicit, so readers never parse file names
# ---------------------------------------------------------------------------

#' Population-coefficient pattern of an exposure's crossbasis columns:
#' `b_<var>_cb<k>` (`^b_enso_cb[0-9]+$` for the ENSO basis).
cb_coef_pattern <- function(var) sprintf("^b_%s_cb[0-9]+$", var)

#' Filename shape of a lag sidecar: `lag_loo_sidecar[_<group>]_<YYYYMMDD>.rds`.
LAG_LOO_SIDECAR_RE <- "^lag_loo_sidecar(_(.+))?_([0-9]{8})\\.rds$"

#' Sidecar schema tag.
LAG_LOO_SIDECAR_SCHEMA <- "lag_loo_sidecar_v1"

#' Digest of a fit's row key: every model-data column except `exclude_cols`
#' (the ENSO crossbasis columns, the only ones that change with the lag),
#' factors as character, row names and attributes dropped (xxhash64).
#'
#' @param data The fit's model data (`fit$data`).
#' @param exclude_cols Column names left out of the key.
#' @return Character digest, with attributes n_rows and key_cols.
fit_row_digest <- function(data, exclude_cols = character(0)) {
  stopifnot(is.data.frame(data))
  keep <- setdiff(names(data), exclude_cols)
  key <- lapply(as.list(data)[keep], function(x) {
    if (is.factor(x)) as.character(x) else as.vector(unclass(x))
  })
  d <- digest::digest(
    list(n = nrow(data), cols = keep, key = unname(key)),
    algo = "xxhash64"
  )
  structure(d, n_rows = nrow(data), key_cols = keep)
}

#' Population draws of the ENSO crossbasis coefficients (`b_<cb_names>`) as a
#' plain S x K matrix with columns `cb_names`. Reads a brmsfit as figures_cb.R's
#' extractors do; anything else goes through posterior::as_draws_matrix()
#' (the tests pass a draws_matrix).
fit_cb_pop_draws <- function(fit, cb_names) {
  vars <- paste0("b_", cb_names)
  dm <- if (inherits(fit, "brmsfit")) {
    brms::as_draws_matrix(fit, variable = vars)
  } else {
    posterior::subset_draws(posterior::as_draws_matrix(fit), variable = vars)
  }
  m <- unclass(dm)[, vars, drop = FALSE]
  matrix(as.numeric(m), nrow = nrow(m), dimnames = list(NULL, cb_names))
}

#' One lag's sidecar entry (pure: the fit has already been read).
#'
#' @param single_lag Integer lag.
#' @param label Column label (05's `lag<NN>` / `lagm<NN>`).
#' @param loo The fit's psis_loo object.
#' @param fit_data The fit's model data (`fit$data`), for the row digest.
#' @param b_draws S x K population crossbasis draws (fit_cb_pop_draws()).
#' @param cb_meta_enso The fit's cb_meta "enso" entry, or NULL.
#' @param exclude_cols Columns left out of the row key [default: the
#'   crossbasis columns, from b_draws and cb_meta_enso].
#' @param files Named character of provenance paths (fit, loo, cb_meta).
lag_sidecar_entry <- function(
  single_lag,
  label,
  loo,
  fit_data,
  b_draws,
  cb_meta_enso = NULL,
  exclude_cols = NULL,
  files = character(0)
) {
  pw <- loo$pointwise
  if (is.null(pw) || !"elpd_loo" %in% colnames(pw)) {
    stop(
      "lag_sidecar_entry(): the loo object has no pointwise elpd_loo (lag ",
      single_lag,
      ").",
      call. = FALSE
    )
  }
  exclude_cols <- exclude_cols %||%
    unique(c(colnames(b_draws), cb_meta_enso$cb_col_names))
  list(
    single_lag = as.integer(single_lag),
    label = label,
    elpd = as.numeric(pw[, "elpd_loo"]),
    pareto_k = .pareto_k_values(loo),
    estimates = data.frame(
      single_lag = as.integer(single_lag),
      elpd_loo = loo$estimates["elpd_loo", "Estimate"],
      se_elpd = loo$estimates["elpd_loo", "SE"],
      p_loo = loo$estimates["p_loo", "Estimate"]
    ),
    pareto = cbind(
      data.frame(single_lag = as.integer(single_lag)),
      pareto_k_summary(loo)
    ),
    row_digest = as.character(fit_row_digest(fit_data, exclude_cols)),
    n_rows = nrow(fit_data),
    b_draws = b_draws,
    cb_meta = cb_meta_enso,
    files = files
  )
}

# Row alignment of a list of lag_sidecar_entry() results.
.check_lag_entries_rows <- function(entries) {
  labs <- vapply(entries, `[[`, "", "label")
  ns <- vapply(
    entries,
    function(e) c(length(e$elpd), length(e$pareto_k), e$n_rows),
    numeric(3)
  )
  if (length(unique(as.vector(ns))) != 1L) {
    stop(
      sprintf(
        "lag sidecar: mismatched N across lags (elpd/k/data rows: %s).",
        paste(
          sprintf("%s %s", labs, apply(ns, 2, paste, collapse = "/")),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }
  .check_lag_digests(stats::setNames(
    vapply(entries, `[[`, "", "row_digest"),
    labs
  ))
}

.check_lag_digests <- function(dg) {
  if (length(unique(dg)) != 1L) {
    stop(
      sprintf(
        "lag sidecar: the row-key digest differs across lags (%s): the fits do not share their rows.",
        paste(sprintf("%s=%s", names(dg), dg), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that an assembled sidecar's rows line up across lags: one N for the
#' pointwise elpd, Pareto-k and n_obs, and one row digest. stop()s otherwise.
check_lag_sidecar_rows <- function(sidecar) {
  n <- c(nrow(sidecar$elpd_pointwise), nrow(sidecar$pareto_k), sidecar$n_obs)
  if (
    length(unique(n)) != 1L ||
      !identical(colnames(sidecar$elpd_pointwise), colnames(sidecar$pareto_k))
  ) {
    stop(
      "lag sidecar: pointwise elpd, Pareto-k and n_obs disagree on N or on the lag columns.",
      call. = FALSE
    )
  }
  .check_lag_digests(sidecar$row_digest)
}

#' Assemble a lag sidecar from per-lag entries (lag_sidecar_entry()), sorted by
#' lag. stop()s when the rows do not line up across lags (N or row digest).
#'
#' @param entries List of lag_sidecar_entry() results (>= 1).
#' @param track "occurrence" | "duration" | "monthly" | "size".
#' @param group Group token (e.g. "enso_east"), or NA / NULL for an ungrouped batch.
#' @param batch Batch config dir path or bare batch name.
build_lag_loo_sidecar <- function(
  entries,
  track,
  group = NA_character_,
  batch
) {
  stopifnot(length(entries) >= 1L)
  entries <- entries[order(vapply(entries, `[[`, integer(1), "single_lag"))]
  .check_lag_entries_rows(entries)
  labs <- vapply(entries, `[[`, "", "label")
  mat <- function(f) {
    m <- do.call(cbind, lapply(entries, `[[`, f))
    colnames(m) <- labs
    m
  }
  cb <- lapply(entries, function(e) {
    list(
      single_lag = e$single_lag,
      b_draws = e$b_draws,
      cb_meta = e$cb_meta,
      files = e$files
    )
  })
  names(cb) <- labs
  out <- list(
    schema = LAG_LOO_SIDECAR_SCHEMA,
    track = track,
    group = if (is.null(group)) NA_character_ else group,
    batch = derive_config_run_tag(batch),
    batch_dir = batch,
    lags = stats::setNames(
      vapply(entries, `[[`, integer(1), "single_lag"),
      labs
    ),
    n_obs = length(entries[[1L]]$elpd),
    elpd_pointwise = mat("elpd"),
    pareto_k = mat("pareto_k"),
    pareto = do.call(rbind, lapply(entries, `[[`, "pareto")),
    estimates = do.call(rbind, lapply(entries, `[[`, "estimates")),
    row_digest = stats::setNames(vapply(entries, `[[`, "", "row_digest"), labs),
    cb = cb,
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  rownames(out$pareto) <- NULL
  rownames(out$estimates) <- NULL
  out
}

#' Latest lag sidecar file of (batch, group) by filename date, then mtime, or
#' NULL. `group` NULL / NA / "" = the ungrouped file. New location only
#' (sidecars never lived in the config dir).
latest_lag_sidecar_file <- function(batch, group = NULL, opt = list()) {
  dir <- make_loo_compare_run_file(opt, "loo_compare_dir", batch = batch)
  if (!dir.exists(dir)) {
    return(NULL)
  }
  f <- list.files(dir, pattern = LAG_LOO_SIDECAR_RE, full.names = TRUE)
  if (!length(f)) {
    return(NULL)
  }
  m <- regmatches(basename(f), regexec(LAG_LOO_SIDECAR_RE, basename(f)))
  want <- if (is.null(group) || is.na(group)) "" else group
  keep <- vapply(m, `[`, "", 3L) == want
  if (!any(keep)) {
    return(NULL)
  }
  f <- f[keep]
  d <- vapply(m[keep], `[`, "", 4L)
  f[order(-xtfrm(d), -as.numeric(file.mtime(f)))][1L]
}

#' Read the latest lag sidecar of (batch, group), checking its schema and row
#' alignment; NULL (or a stop() when `required`) when there is none. The file
#' path is attached as attr(, "file").
read_lag_loo_sidecar <- function(
  batch,
  group = NULL,
  opt = list(),
  required = TRUE
) {
  f <- latest_lag_sidecar_file(batch, group = group, opt = opt)
  if (is.null(f)) {
    if (required) {
      stop(
        sprintf(
          "[loo_lags] no lag_loo_sidecar for group '%s' in %s (run 05_fit_postprocess_loo_lags.R --lag_sidecar).",
          group %||% "",
          make_loo_compare_run_file(opt, "loo_compare_dir", batch = batch)
        ),
        call. = FALSE
      )
    }
    return(NULL)
  }
  x <- readRDS(f)
  if (!identical(x$schema, LAG_LOO_SIDECAR_SCHEMA)) {
    stop(
      "[loo_lags] ",
      f,
      ": unknown sidecar schema '",
      x$schema %||% "<none>",
      "'.",
      call. = FALSE
    )
  }
  check_lag_sidecar_rows(x)
  attr(x, "file") <- f
  x
}
