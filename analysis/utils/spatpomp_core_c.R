# utils/spatpomp_core_c.R — the Csnippet spatPomp builder.
#
# Public surface:
#   build_spatpomp_tu_c(frame_tu, panel_tu_warmup, draws_occ, window,
#                       n_draws_in_filter, obs_data, fixed_draw_idx)
#     -> spatPomp object with Csnippet rprocess / rinit.
#
# Design:
#   * Sources spatpomp_core.R so all shared helpers (build_obs_for_tu, build_covar_for_tu,
#     make_block_list, run_bpfilter_tu, filtered_states_to_warmup, etc.) are
#     inherited without modification.
#   * fixed_draw_idx is REQUIRED (errors if NULL). The N-independent-runs design
#     always provides it; the C variant is not supported without it.
#   * LP is split into:
#       lp_static  — precomputed once per (unit, month); captures intercept FE,
#                    intercept REs, climate FE + RE slopes × covariate values,
#                    seasonal spline, and basehaz.
#       ttlo_coef_b{0..6} — per-unit TTLO RE contributions (7 bins).
#       spill_{o,oct,otu}_{t,f} — per-unit spillover TRUE/FALSE RE contributions.
#   * Topology arrays (adm1/country/tu membership) are embedded as static const
#     arrays in the generated Csnippet string.
#   * Warmup state is passed via additional per-unit covariate columns read by
#     the rinit Csnippet.

suppressPackageStartupMessages(library(glue))

# Source the sibling core module only when not already loaded (avoid
# double-load when the dispatcher has already sourced it before sourcing
# this file).
if (!exists("build_obs_for_tu", mode = "function")) {
  source(file.path(dirname(sys.frame(1)$ofile), "spatpomp_core.R"))
}

# ---------------------------------------------------------------------------
# Nested-slope helper (draws_sub with S=1; always uses row 1)
# ---------------------------------------------------------------------------

# Evaluate the sum of RE slopes across all nested levels for one term.
# Replicates the `nested_slope()` local closure in evaluate_lp() but operates
# on draws_sub directly (S=1 fixed-draw subset) and is callable outside LP.
#
# @param draws_sub  1-draw draws bundle (output of .subset_draws with S=1).
# @param term       character; brms term name (e.g. "mean_precip").
# @param frame_tu   location frame for the TU (N rows).
# @return numeric vector of length N.
.ns1 <- function(draws_sub, term, frame_tu, drop_region = FALSE) {
  N   <- nrow(frame_tu)
  rst <- draws_sub$re_slope[[term]]
  if (is.null(rst)) return(rep(0.0, N))
  # Dual-key lookup: try primary (drop_region) format first; fall back to the
  # alternate format for mixed-nesting fits (e.g. grouped brms_cb) where
  # intercept REs use bare keys but ENSO RE slopes use AFRO-prefixed keys.
  re_keys2 <- function(bucket) list(
    primary = re_level_key(frame_tu, bucket, drop_region),
    alt     = re_level_key(frame_tu, bucket, !drop_region))
  lk <- function(mat, bucket) {
    if (is.null(mat) || ncol(mat) == 0L) return(rep(0.0, N))
    keys <- re_keys2(bucket)
    idx  <- match(keys$primary, colnames(mat))
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(keys$alt[miss], colnames(mat))
    ifelse(is.na(idx), 0.0, mat[1L, idx])
  }
  lk(rst$region,  "region")  +
  lk(rst$country, "country") +
  lk(rst$adm1,    "adm1")    +
  lk(rst$adm2,    "adm2")
}

# ---------------------------------------------------------------------------
# Precompute lp_static — static LP components (all except TTLO and spillover)
# ---------------------------------------------------------------------------

#' Precompute the static part of the linear predictor for every (unit, time)
#' in the covariate table.
#'
#' Calls evaluate_lp() with NA spillover and NA u_days so that TTLO and
#' spillover contributions are suppressed (both return 0 for NA inputs).
#' The result captures: intercept FE, nested intercept REs, climate FE + RE
#' slope × covariate values, seasonal spline, and basehaz contribution.
#'
#' @param draws_sub  1-draw draws bundle.
#' @param frame_tu   location frame for the TU.
#' @param covar_tbl  covariate table from build_covar_for_tu (long, sorted by
#'                   time then unit).
#' @param enso_names character vector of ENSO covariate column names to pass.
#' @param cb_names   DLNM cross-basis columns to inject. Cross-basis terms are
#'                   neither TTLO- nor spillover-dependent, so they are fully
#'                   absorbed into the static LP here (no Csnippet change needed).
#'                   Empty when the fit carries no cross-basis.
#'
#'                   Accepts either form: an UNNAMED vector of covariate column
#'                   names (each maps to the identically-named model term -- the
#'                   occurrence track's case, where the bare `enso_cb<k>` columns
#'                   ARE that fit's basis), or a NAMED vector mapping MODEL TERM
#'                   -> COVARIATE COLUMN. The named form is what the monthly
#'                   and duration tracks must pass: every track selects its ENSO
#'                   single-lag independently (Horn: occ 11, mon 11, dur 8), so
#'                   their columns live under the `.SP_MON_CB_PREFIX` /
#'                   `.SP_DUR_CB_PREFIX` namespaces and have to be renamed back to
#'                   the bare model term before `evaluate_lp()`'s `_cb[0-9]+$`
#'                   whitelist will match them against the fitted coefficients.
#'                   Passing the unnamed form for those tracks silently evaluates
#'                   their coefficients against the OCCURRENCE basis.
#' @return numeric vector of length nrow(covar_tbl).
precompute_lp_static <- function(draws_sub, frame_tu, covar_tbl, enso_names,
                                 cb_names = character(0), drop_region = FALSE) {
  U         <- nrow(frame_tu)
  gids      <- as.character(frame_tu$gid)
  clim_cols <- intersect(c("mean_precip", "dev_precip",
                            "mean_temp_daymean", "dev_temp_daymean"),
                          names(covar_tbl))

  # Empirical baseline-rate anchor, offset(log(base_rate_occ)) /
  # offset(log(base_rate_dur)) -- precompute_dur_lp_static() delegates here, so
  # both tracks are covered by this one line, each with its own bundle. Static
  # per-unit, exactly like the mon block's pop offset below, so it
  # folds into lp_static at R time and the Csnippet needs no change (and
  # .SM_CACHE_SCHEMA_VERSION does not move).
  static_offset <- .anchor_log_offset(draws_sub, frame_tu, NULL,
                                      "precompute_lp_static")

  # Resolve cb_names into a term -> column map (see @param). Unnamed entries map
  # to themselves, which is byte-identical to the previous plain-vector handling.
  cb_col  <- as.character(cb_names)
  cb_term <- names(cb_names) %||% rep("", length(cb_names))
  cb_term <- ifelse(nzchar(cb_term), cb_term, cb_col)
  cb_keep <- cb_col %in% names(covar_tbl)
  cb_col  <- cb_col[cb_keep]; cb_term <- cb_term[cb_keep]

  lp_out <- rep(NA_real_, nrow(covar_tbl))
  times  <- sort(unique(covar_tbl$time))

  for (t in times) {
    mask_t  <- covar_tbl$time == t
    rows_t  <- covar_tbl[mask_t, , drop = FALSE]

    # Reorder to match frame_tu unit order (evaluate_lp expects frame_tu row
    # order for xrow vectors).
    perm <- match(gids, rows_t$unit)
    rows_ord <- rows_t[perm, , drop = FALSE]

    xrow <- setNames(
      lapply(clim_cols, function(col) rows_ord[[col]]),
      clim_cols
    )
    xrow$t_start_seas <- rows_ord$t_start_seas
    for (ec in intersect(enso_names, names(rows_t)))
      xrow[[ec]] <- rows_ord[[ec]]
    for (i in seq_along(cb_col))
      xrow[[cb_term[i]]] <- rows_ord[[cb_col[i]]]

    bh   <- basehaz_eval(draws_sub, xrow$t_start_seas, N = U)
    lp_t <- evaluate_lp(
      draws_sub, s = 1L, xrow,
      o      = rep(NA_integer_, U),
      o_ct   = rep(NA_integer_, U),
      o_tu   = rep(NA_integer_, U),
      u_days = rep(NA_real_,    U),
      frame_tu    = frame_tu,
      basehaz_val = bh,
      drop_region = drop_region
    )
    lp_t <- lp_t + static_offset

    # Map lp_t (frame_tu order) back to covar_tbl row order for time t.
    inv_perm <- match(rows_t$unit, gids)
    lp_out[mask_t] <- lp_t[inv_perm]
  }

  lp_out
}

# ---------------------------------------------------------------------------
# Precompute TTLO coefficients — per-unit, per-bin (7 bins)
# ---------------------------------------------------------------------------

#' Precompute per-unit TTLO RE + FE contributions for each of the 7 bins.
#'
#' Bin labels and column name mapping (matching ttlo_col_to_level in 09a):
#'   bin 0 [0,1):   ref (no FE); RE key = "time_to_last_outbreak01"
#'   bin 1 [1,2):   FE col "time_to_last_outbreak12"
#'   bin 2 [2,3):   FE col "time_to_last_outbreak23"
#'   ...
#'   bin 6 [6,Inf]: FE col "time_to_last_outbreak6Inf"
#'
#'   bin 7 none:    FE col "time_to_last_outbreaknone" -- only when the fit
#'                  carries the "no recorded prior onset" level
#'
#' @return data.frame with columns ttlo_coef_b0..ttlo_coef_b6 (U rows), plus
#'   ttlo_coef_b7 when the fit has the level "none" (absent otherwise, so the
#'   covariate table of a fit without it is unchanged).
precompute_ttlo_coefs <- function(draws_sub, frame_tu, drop_region = FALSE) {
  N  <- nrow(frame_tu)
  fe <- draws_sub$fe[1L, , drop = TRUE]

  bin_col <- c(
    b0 = "time_to_last_outbreak01",
    b1 = "time_to_last_outbreak12",
    b2 = "time_to_last_outbreak23",
    b3 = "time_to_last_outbreak34",
    b4 = "time_to_last_outbreak45",
    b5 = "time_to_last_outbreak56",
    b6 = "time_to_last_outbreak6Inf"
  )

  out <- vector("list", 7L)
  names(out) <- paste0("ttlo_coef_b", 0:6)

  for (k in 0:6) {
    col_nm <- bin_col[[paste0("b", k)]]
    # Bin 0 is the reference level: FE contribution = 0 (absorbed into intercept)
    fe_val <- if (k == 0L) 0.0
              else if (col_nm %in% names(fe)) fe[[col_nm]]
              else 0.0
    re_val <- .ns1(draws_sub, col_nm, frame_tu, drop_region)
    out[[paste0("ttlo_coef_b", k)]] <- fe_val + re_val
  }
  if (.TTLO_NONE_COL %in% names(fe))
    out[["ttlo_coef_b7"]] <- fe[[.TTLO_NONE_COL]] +
      .ns1(draws_sub, .TTLO_NONE_COL, frame_tu, drop_region)

  as.data.frame(out)
}

#' .ttlo_bin_c
#'
#' The C chain mapping a days-since-last-onset value to its TTLO bin index
#' (0..6 = .TTLO_LEVELS; 7 = "none" when `has_none`), emitted in the three
#' kernel blocks (occurrence, monthly, duration). `var` is the C expression of
#' the days value, `indent` the leading spaces, `spaced` whether `+ 1.0` is
#' written with spaces (the blocks differ). Without the level the two lines are
#' pinned by a golden (tests/test_ttlo_none_level.R); with it, a value at
#' or above base.R::.TTLO_NONE_DAYS is bin 7.
.ttlo_bin_c <- function(var, indent = "    ", spaced = TRUE, has_none = FALSE) {
  plus <- if (isTRUE(spaced)) sprintf("(%s + 1.0)", var) else sprintf("(%s+1.0)", var)
  chain <- "(yrs<1.0)?0:(yrs<2.0)?1:(yrs<3.0)?2:(yrs<4.0)?3:(yrs<5.0)?4:(yrs<6.0)?5:6"
  if (isTRUE(has_none))
    chain <- sprintf("(%s >= %s) ? 7 : %s", var, sprintf("%.1f", .TTLO_NONE_DAYS), chain)
  paste0(indent, "double yrs = (", var, " < 0.0 ? 1.0/365.0 : ", plus, "/365.0);\n",
         indent, "int bin = ", chain, ";\n")
}

# ---------------------------------------------------------------------------
# Precompute spillover coefficients — per-unit, TRUE/FALSE levels
# ---------------------------------------------------------------------------

#' Precompute per-unit spillover (ongoing outbreak) RE + FE contributions.
#'
#' For each of the three spillover terms (o, o_ct, o_tu) and both factor
#' levels (TRUE/FALSE):
#'   TRUE  contribution = fe[[term]] + nested_slope(term)[u]
#'   FALSE contribution = nested_slope(paste0(term, "FALSE"))[u]
#'
#' NA-TU units have o_tu contributions forced to 0 (they match fit-time NA
#' behaviour where ongoing_outbreak_tu_lag1 = NA → term absent → no weight).
#'
#' @return data.frame with 6 columns, U rows:
#'         spill_o_t, spill_o_f, spill_oct_t, spill_oct_f,
#'         spill_otu_t, spill_otu_f.
precompute_spillover_coefs <- function(draws_sub, frame_tu, drop_region = FALSE) {
  N  <- nrow(frame_tu)
  fe <- draws_sub$fe[1L, , drop = TRUE]

  make_pair <- function(term) {
    if (!(term %in% names(fe)))
      return(list(t = rep(0.0, N), f = rep(0.0, N)))
    fe_val <- fe[[term]]
    re_t   <- .ns1(draws_sub, term,                    frame_tu, drop_region)
    re_f   <- .ns1(draws_sub, paste0(term, "FALSE"),   frame_tu, drop_region)
    list(t = fe_val + re_t, f = re_f)
  }

  sp_o   <- make_pair("ongoing_outbreak_lag1")
  sp_oct <- make_pair("ongoing_outbreak_country_lag1")
  sp_otu <- make_pair("ongoing_outbreak_tu_lag1")

  # Force o_tu coefficients to 0 for NA-TU locations.
  na_tu <- is.na(frame_tu$tu)
  if (any(na_tu)) {
    sp_otu$t[na_tu] <- 0.0
    sp_otu$f[na_tu] <- 0.0
  }

  data.frame(
    spill_o_t   = sp_o$t,
    spill_o_f   = sp_o$f,
    spill_oct_t = sp_oct$t,
    spill_oct_f = sp_oct$f,
    spill_otu_t = sp_otu$t,
    spill_otu_f = sp_otu$f
  )
}

#' Precompute per-unit neighbour-contagion coefficients beta_nbr(u) for the C
#' rprocess. Mirrors precompute_spillover_coefs(): the fixed effect plus the
#' summed nested RE slope (region/country/adm1/adm2 via .ns1), so grouped
#' Design-B fits get the correct per-unit slope. Units whose (group) fit carries
#' no nbr_* FE get beta_nbr = 0 (the grep yields nothing), exactly as
#' make_pair()'s `!(term %in% names(fe))` guard returns 0 for an absent
#' spillover term. Unlike spillover this is a PLAIN continuous linear term, so
#' there is no FALSE-factor branch; the per-step neighbour PRESSURE is computed
#' in the Csnippet from the `o` state and the embedded WNBR kernel.
#'
#' @return numeric vector length nrow(frame_tu): beta_nbr per unit.
precompute_nbr_coefs <- function(draws_sub, frame_tu, drop_region = FALSE) {
  N  <- nrow(frame_tu)
  fe <- draws_sub$fe[1L, , drop = TRUE]
  nbr_col <- grep(nbrCovariatePattern(), names(fe), value = TRUE)
  if (length(nbr_col) == 0L) return(rep(0.0, N))
  nbr_col <- nbr_col[1L]
  as.numeric(fe[[nbr_col]]) + .ns1(draws_sub, nbr_col, frame_tu, drop_region)
}

# ---------------------------------------------------------------------------
# Exogenous observed o_tu fold (country-partition mode only)
# ---------------------------------------------------------------------------

#' Fold OBSERVED cross-country contagion into lp_static for single-country
#' partitions.
#'
#' Under partition_by=country the internal `tu` indicator is collapsed to a
#' constant, so update_lags() yields o_tu = 0 for every unit at every step and
#' the rprocess always takes the FALSE spillover branch (adds spill_otu_f). The
#' positive cross-country contagion term is therefore lost, clipping the
#' continent-wide occurrence peaks.
#'
#' This restores it WITHOUT touching the (frozen) 07a kernel or the Csnippet:
#' for every (unit, time) whose lag-1 calendar month had an observed outbreak in
#' ANOTHER country of the same TRUE transmission unit, we add
#'   delta[unit] = spill_otu_t[unit] - spill_otu_f[unit]
#' to lp_static. Combined with the rprocess's unconditional spill_otu_f, the
#' effective contribution becomes spill_otu_t exactly when the observed term is
#' active — algebraically identical to the TU-partitioned run, but driven by the
#' observed panel instead of the (single-country) simulated state.
#'
#' No-op unless `otu_obs` is supplied AND frame_tu carries a `tu_true` column
#' (set only by smoothing_apply_partition in country mode); the TU-partitioned
#' and ungrouped paths are completely unchanged.
#'
#' @param lp_static   numeric, length nrow(covar_tbl) (frame_tu/covar_tbl order).
#' @param covar_tbl   long covariate table (columns unit, time, ...).
#' @param filt_months sorted Dates; filt_months[t] is the calendar month for
#'                    covar time index t.
#' @param frame_tu    partition location frame (must carry `gid`; `tu_true` when
#'                    folding applies).
#' @param spill_coefs data.frame from precompute (spill_otu_t / spill_otu_f).
#' @param otu_obs     data.frame(t_start, o_tu_obs): observed cross-country
#'                    activity (0/1) by calendar month, or NULL.
#' @return lp_static with the observed o_tu term folded in (or unchanged).
.fold_observed_otu_into_lp <- function(lp_static, covar_tbl, filt_months,
                                       frame_tu, spill_coefs, otu_obs) {
  if (is.null(otu_obs) || is.null(frame_tu$tu_true)) return(lp_static)
  if (nrow(otu_obs) == 0L) return(lp_static)

  # Per-unit contribution gained when o_tu flips 0 -> 1.
  delta <- spill_coefs$spill_otu_t - spill_coefs$spill_otu_f
  names(delta) <- as.character(frame_tu$gid)

  # Observed cross-country activity keyed by floored calendar month.
  obs_lookup <- stats::setNames(
    as.numeric(otu_obs$o_tu_obs),
    as.character(lubridate::floor_date(as.Date(otu_obs$t_start), "month")))

  # Lag-1: the o_tu read at filter month m reflects activity at month (m - 1).
  # Previous month first-day = floor( first-of-month - 1 day ).
  fm      <- lubridate::floor_date(as.Date(filt_months), "month")
  lag_key <- as.character(lubridate::floor_date(fm - 1, "month"))
  o_tu_by_time <- unname(obs_lookup[lag_key])
  o_tu_by_time[is.na(o_tu_by_time)] <- 0

  active_times <- which(o_tu_by_time > 0.5)
  if (length(active_times) == 0L) return(lp_static)
  for (t in active_times) {
    mask_t <- covar_tbl$time == t
    if (!any(mask_t)) next
    units_t <- as.character(covar_tbl$unit[mask_t])
    lp_static[mask_t] <- lp_static[mask_t] + delta[units_t]
  }
  lp_static
}

# ---------------------------------------------------------------------------
# Group-aware precompute dispatch
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Grouped precompute stitching: one copy of the per-group loop shared by
# precompute_lp_components() / precompute_dur_lp_components() /
# precompute_mon_lp_components().
# ---------------------------------------------------------------------------

#' .stitch_grouped_lp
#'
#' For each ENSO group (frame_tu$group, in order of first appearance): take the
#' group's draws subset (stop with `missing_msg` + group when absent), compute
#' its static LP on the group's covariate rows and its time-since-last-outbreak
#' (TTLO) and spillover coefficient rows, and stitch them into full-length
#' outputs (TTLO/spillover via an NA template shaped like the first group's
#' table). `extra_fn(draws_g, frame_g[, g])` -- a track's additional per-group
#' quantity -- is collected by group name.
#'
#' @param lp_fn function(draws_g, frame_g, covar_g) -> per-row static LP
#' @param ttlo_fn,spill_fn function(draws_g, frame_g, drop_region = TRUE)
#' @return list(lp, ttlo, spill, extra (named by group), groups, grp_vec)
#' @noRd
.stitch_grouped_lp <- function(frame_tu, covar_tbl, draws_by_group, lp_fn,
                               ttlo_fn, spill_fn, missing_msg, extra_fn = NULL) {
  U          <- nrow(frame_tu)
  grp_vec    <- as.character(frame_tu$group)
  groups     <- unique(grp_vec)
  lp_full    <- rep(NA_real_, nrow(covar_tbl))
  ttlo_full  <- NULL
  spill_full <- NULL
  extra      <- setNames(vector("list", length(groups)), groups)

  for (g in groups) {
    row_mask_g <- grp_vec == g
    frame_g    <- frame_tu[row_mask_g, , drop = FALSE]
    gids_g     <- as.character(frame_g$gid)
    draws_g    <- draws_by_group[[g]]
    if (is.null(draws_g)) stop(missing_msg, g)

    cov_mask_g <- covar_tbl$unit %in% gids_g
    lp_full[cov_mask_g] <- lp_fn(draws_g, frame_g, covar_tbl[cov_mask_g, , drop = FALSE])

    ttlo_g  <- ttlo_fn(draws_g, frame_g, drop_region = TRUE)
    spill_g <- spill_fn(draws_g, frame_g, drop_region = TRUE)

    if (is.null(ttlo_full)) {
      ttlo_full            <- ttlo_g[rep(1L, U), , drop = FALSE]
      ttlo_full[]          <- NA_real_
      spill_full           <- spill_g[rep(1L, U), , drop = FALSE]
      spill_full[]         <- NA_real_
      rownames(ttlo_full)  <- NULL
      rownames(spill_full) <- NULL
    }
    ttlo_full[row_mask_g, ]  <- ttlo_g
    spill_full[row_mask_g, ] <- spill_g

    if (!is.null(extra_fn)) {
      x <- if (length(formals(extra_fn)) >= 3L) extra_fn(draws_g, frame_g, g)
           else extra_fn(draws_g, frame_g)
      if (!is.null(x)) extra[[g]] <- x
    }
  }
  list(lp = lp_full, ttlo = ttlo_full, spill = spill_full, extra = extra,
       groups = groups, grp_vec = grp_vec)
}

#' Precompute lp_static / TTLO / spillover, with optional multi-group routing.
#'
#' Ungrouped (grouping_sub = NULL): identical to the three standalone
#' precompute_* calls — zero behaviour change.
#'
#' Grouped (Design B, multi-group per TU): segments frame_tu rows by
#' frame_tu$group and runs the existing precompute_* functions per group against
#' that group's single-draw bundle, with drop_region = TRUE (grouped fits drop
#' the AFRO_region RE level). Results are stitched back into full-frame order:
#'   - lp_static: assigned into a length-nrow(covar_tbl) vector via the per-group
#'     covar-row mask (covar_tbl$unit %in% group gids); precompute_lp_static
#'     returns rows in the masked-subset order, so the assignment is order-safe.
#'   - ttlo / spillover: full U-row data.frames filled where frame_tu$group == g.
#'     frame_g preserves frame_tu row order for group g, so the per-group result
#'     aligns with the logical mask. build_covar_for_tu_c later joins these by
#'     unit, so per-group fill order is safe downstream.
#'
#' The union of enso_names / cb_names is passed to every group: evaluate_lp
#' ignores terms absent from a group's draws, so no per-group name filtering is
#' required.
#'
#' @param frame_tu     location frame for the TU (must carry a `group` column
#'                     under grouping).
#' @param covar_tbl    covariate table (long, from build_covar_for_tu).
#' @param enso_names   union of ENSO covariate column names.
#' @param cb_names     union of cross-basis covariate column names.
#' @param draws_sub    single-draw bundle (representative; used only ungrouped).
#' @param grouping_sub list(draws_sub_by_group=...) from .build_grouping_sub, or
#'                     NULL for the ungrouped path.
#' @return list(lp = <numeric>, ttlo = <data.frame>, spill = <data.frame>).
precompute_lp_components <- function(frame_tu, covar_tbl, enso_names, cb_names,
                                     draws_sub, grouping_sub = NULL,
                                     use_nbr = FALSE) {
  if (is.null(grouping_sub)) {
    out <- list(
      lp    = precompute_lp_static(draws_sub, frame_tu, covar_tbl,
                                   enso_names, cb_names),
      ttlo  = precompute_ttlo_coefs(draws_sub, frame_tu),
      spill = precompute_spillover_coefs(draws_sub, frame_tu)
    )
    if (isTRUE(use_nbr))
      out$nbr <- precompute_nbr_coefs(draws_sub, frame_tu)
    return(out)
  }

  if (!("group" %in% names(frame_tu)))
    stop("precompute_lp_components: grouped mode requires a `group` column on ",
         "frame_tu.")

  st <- .stitch_grouped_lp(
    frame_tu, covar_tbl, grouping_sub$draws_sub_by_group,
    lp_fn    = function(d, f, cv) precompute_lp_static(d, f, cv, enso_names, cb_names,
                                                       drop_region = TRUE),
    ttlo_fn  = precompute_ttlo_coefs, spill_fn = precompute_spillover_coefs,
    missing_msg = "precompute_lp_components: no draws subset for group ",
    extra_fn = if (isTRUE(use_nbr)) function(d, f)
      precompute_nbr_coefs(d, f, drop_region = TRUE))

  out <- list(lp = st$lp, ttlo = st$ttlo, spill = st$spill)
  # beta_nbr is a plain per-unit vector (not a multi-column data.frame like
  # ttlo/spill); 0 is the correct value for a group whose fit lacks the FE, so a
  # zero-initialised full vector doubles as the "absent FE" default.
  if (isTRUE(use_nbr)) {
    nbr_full <- rep(0.0, nrow(frame_tu))
    for (g in st$groups) nbr_full[st$grp_vec == g] <- st$extra[[g]]
    out$nbr <- nbr_full
  }
  out
}

# ---------------------------------------------------------------------------
# Sibling precompute functions for the DURATION (cessation-hazard) model (C)
# ---------------------------------------------------------------------------
# Mirrors the occurrence precompute trio exactly. Differences
# from occurrence:
#   * basehaz_val is always 0 in precompute_dur_lp_static() -- NOT because it
#     is folded in statically (unlike occurrence's t_start_seas basehaz), but
#     because basehaz_eval() (the generic dispatcher precompute_lp_static()
#     calls) does not recognise the duration model's "factor" mso-band kind
#     and silently returns 0 for it -- exactly the value the R
#     reference wants here, since the mso-band hazard is STATE-dependent
#     (mso changes every step) and can never be folded into a per-(unit,time)
#     build-time constant the way occurrence's t_start_seas-keyed basehaz can.
#     The real mso-band contribution is computed at Csnippet RUNTIME from the
#     DBH lookup table below (mirrors the TTLO yrs-bin ternary chain).
#   * TTLO/spillover coefficients are drawn from the DURATION draws bundle's
#     own fe/re_slope (a real duration model can, and typically will, carry
#     different TTLO/spillover coefficients than the occurrence model whose
#     outbreak episodes it shares) -- reuses precompute_ttlo_coefs /
#     precompute_spillover_coefs verbatim (already generic over any draws_sub).

#' Precompute the static part of the duration LP (climate + intercept + RE +
#' ENSO/cross-basis; basehaz always 0 here -- see block comment above).
precompute_dur_lp_static <- function(draws_dur_sub, frame_tu, covar_tbl,
                                     enso_names, cb_names = character(0),
                                     drop_region = FALSE) {
  precompute_lp_static(draws_dur_sub, frame_tu, covar_tbl, enso_names,
                       cb_names, drop_region)
}

#' Sibling of precompute_ttlo_coefs for the duration model.
precompute_dur_ttlo_coefs <- function(draws_dur_sub, frame_tu, drop_region = FALSE) {
  precompute_ttlo_coefs(draws_dur_sub, frame_tu, drop_region)
}

#' Sibling of precompute_spillover_coefs for the duration model.
precompute_dur_spillover_coefs <- function(draws_dur_sub, frame_tu, drop_region = FALSE) {
  precompute_spillover_coefs(draws_dur_sub, frame_tu, drop_region)
}

#' Precompute the months-since-onset (mso) baseline-hazard lookup table.
#'
#' basehaz_dur_eval() is an FE-only lookup (no per-unit RE -- confirmed by its
#' source: draws_dur$fe[s, ] indexed by mso_band column, with no
#' re_intercept/re_slope term at all), so the mso-band contribution is a
#' GLOBAL scalar per bin, not per-unit -- unlike TTLO's per-unit RE-adjusted
#' bins. This lets the whole table collapse to a small `static const double`
#' array in the generated Csnippet, indexed by a runtime bin computed from
#' the breakpoints (mirrors the TTLO yrs-bin ternary chain, generated
#' dynamically here instead of hardcoded since duration.mso_breaks is a
#' config value, unlike the fixed 7-bin .TTLO_BREAKS).
#'
#' @param draws_dur_sub 1-draw duration draws bundle.
#' @return list(kind = "factor"|"none", breaks = <numeric or NULL>,
#'              table = <numeric, length K, or empty>).
#'         kind = "none": no mso-band term in this fit; Csnippet generates a
#'         constant-0 contribution (mirrors basehaz_dur_eval()'s "none"/
#'         unrecognised-kind fallback).
precompute_dur_basehaz_table <- function(draws_dur_sub) {
  bh   <- draws_dur_sub$basehaz
  kind <- bh$kind %||% "none"
  if (!identical(kind, "factor"))
    return(list(kind = "none", breaks = NULL, table = numeric(0)))

  breaks     <- bh$breaks %||% .SIM_DUR_MSO_BREAKS
  left_edges <- breaks[-length(breaks)]   # one representative value per bin
  tbl        <- basehaz_dur_eval(draws_dur_sub, s = 1L, mso = left_edges,
                                 N = length(left_edges))
  list(kind = "factor", breaks = breaks, table = as.numeric(tbl))
}

#' Group-aware dispatch for duration precomputes (sibling of
#' precompute_lp_components).
#'
#' Ungrouped: standalone calls for lp, ttlo, spill, plus a single mso-hazard
#' table.
#' Grouped: per-group dispatch with drop_region = TRUE for lp/ttlo/spill,
#' stitched to full-frame order exactly like the occurrence sibling.
#' The mso-hazard table is genuinely per-group (each group's duration draws
#' bundle carries its own FE coefficients), so this additionally returns
#' `grp_id` (0-based group index per unit, matching `groups` order) and
#' `bh_by_group` (list of per-group precompute_dur_basehaz_table() results) --
#' the Csnippet generator embeds one DBH[] row per group and looks up
#' DBH[grp_id[u]][bin]. All groups must share identical `breaks` (a single
#' `duration.mso_breaks` config value, never group-varying in practice); this
#' is checked explicitly and errors loudly rather than silently mismatching
#' bin boundaries against per-group coefficient tables.
#'
#' @return list(lp, ttlo, spill, bh, grp_id = NULL, groups = NULL) ungrouped;
#'         list(lp, ttlo, spill, bh_by_group, grp_id, groups) grouped.
precompute_dur_lp_components <- function(frame_tu, covar_tbl, enso_names,
                                         cb_names, draws_dur_sub,
                                         grouping_dur_sub = NULL) {
  if (is.null(grouping_dur_sub)) {
    return(list(
      lp    = precompute_dur_lp_static(draws_dur_sub, frame_tu, covar_tbl,
                                       enso_names, cb_names),
      ttlo  = precompute_dur_ttlo_coefs(draws_dur_sub, frame_tu),
      spill = precompute_dur_spillover_coefs(draws_dur_sub, frame_tu),
      bh    = precompute_dur_basehaz_table(draws_dur_sub),
      grp_id = NULL,
      groups = NULL
    ))
  }

  if (!("group" %in% names(frame_tu)))
    stop("precompute_dur_lp_components: grouped mode requires `group` column on frame_tu.")

  st <- .stitch_grouped_lp(
    frame_tu, covar_tbl, grouping_dur_sub$draws_sub_by_group,
    lp_fn    = function(d, f, cv) precompute_dur_lp_static(d, f, cv, enso_names, cb_names,
                                                           drop_region = TRUE),
    ttlo_fn  = precompute_dur_ttlo_coefs, spill_fn = precompute_dur_spillover_coefs,
    missing_msg = "precompute_dur_lp_components: no duration draws subset for group ",
    extra_fn = function(d, f) precompute_dur_basehaz_table(d))
  groups <- st$groups; grp_vec <- st$grp_vec
  lp_full <- st$lp; ttlo_full <- st$ttlo; spill_full <- st$spill
  grp_id      <- match(grp_vec, groups) - 1L
  bh_by_group <- st$extra

  # All groups must share identical mso-band breakpoints (single global
  # duration.mso_breaks config value; see docstring). Loud error, no silent
  # mismatch -- same convention as extract_duration_draws()'s unknown-backend
  # guard (test_duration_kernel.R #5).
  kinds <- vapply(bh_by_group, `[[`, character(1), "kind")
  if (length(unique(kinds)) > 1L)
    stop("precompute_dur_lp_components: groups disagree on basehaz$kind ",
         "(", paste(unique(kinds), collapse = ", "), ") -- ",
         "duration.mso_breaks/basehaz kind must be identical across groups.")
  if (identical(kinds[[1]], "factor")) {
    brk1 <- bh_by_group[[1]]$breaks
    same_breaks <- vapply(bh_by_group, function(b) isTRUE(all.equal(b$breaks, brk1)),
                          logical(1))
    if (!all(same_breaks))
      stop("precompute_dur_lp_components: groups disagree on basehaz$breaks -- ",
           "duration.mso_breaks must be identical across groups.")
  }

  list(lp = lp_full, ttlo = ttlo_full, spill = spill_full,
       bh_by_group = bh_by_group, grp_id = grp_id, groups = groups)
}

# ---------------------------------------------------------------------------
# Sibling precompute functions for the MONTHLY-CASE (`mon`) model (C variant)
# ---------------------------------------------------------------------------
# Mirrors the occurrence trio (precompute_lp_static, precompute_ttlo_coefs,
# precompute_spillover_coefs) with basehaz_val = 0 (a count model, not a
# survival model) and three differences that are all forced by
# build_monthly_formula_parts() (utils/config.R) rather than chosen here:
#
#   * OFFSETS. mon is offset(log(pop)) [toggled] + offset(log_offset)
#     [toggled] + its own empirical anchor offset(log(base_rate_mon)).
#     log(pop) and the anchor are static per-unit constants and fold in once
#     -- the C rprocess never reads pop -- but log_offset is a per-(unit, time)
#     covariate column, so it is added inside the per-time loop. Every toggle
#     is read off the DRAWS BUNDLE (offset_pop / offset_logoffset), never
#     hardcoded -- the same bundle-driven rule the v1 kernel follows, so the
#     offset set is asserted in exactly one place (the extractor) and tested
#     once. `pop_elasticity: free` is NOT folded here: the `log_pop_c`
#     deviation term is a fitted COEFFICIENT, not an offset, so it is derived
#     inside evaluate_lp() from `frame_tu$pop` and the bundle's
#     `log_pop_centres` (simulator_core.R) -- one derivation serves every
#     kernel. `off_pop` below stays the elasticity-1 part.
#
#   * COVARIATE NAMING. mon's panel cross-bases live in the covariate table
#     under the .SP_MON_CB_PREFIX namespace ("monCB_enso_cb1", ...) because
#     occurrence and monthly cases pick their ENSO single-lag independently and
#     would otherwise collide on one `enso_cb<k>` column. evaluate_lp() must
#     see them under the BARE fitted term, so callers pass a term -> column
#     map (`mon_covar_map`) instead of a plain `cb_names` vector.
#
#   * THE mso CROSS-BASIS IS LATENT. It has no panel column at all, so it
#     cannot be folded into lp_static; it is precomputed as a per-unit x
#     per-mso-value coefficient table (precompute_mon_mso_coefs) and looked up
#     at Csnippet runtime from the latent mv_new[u]. Structurally this is
#     precompute_ttlo_coefs()'s per-unit-RE bin table, NOT the duration
#     basehaz's global `static const double` array (which works only because
#     basehaz_dur_eval() is FE-only).

#' Precompute the static (non-mso) part of the monthly-case LP for every
#' (unit, time) in covar_tbl.
#'
#' @param draws_mon_sub 1-draw monthly draws bundle.
#' @param frame_tu      location frame for the TU.
#' @param covar_tbl     covariate table from build_covar_for_tu (long).
#' @param mon_covar_map named character vector mapping MODEL TERM -> COVARIATE
#'   COLUMN for mon's panel covariates (harmonics map to themselves; prefixed
#'   cross-bases map monCB_enso_cb1 -> enso_cb1). Entries whose column is
#'   absent from covar_tbl are dropped silently, matching evaluate_lp()'s own
#'   "term not present => contributes nothing" convention.
#' @return numeric vector of length nrow(covar_tbl).
precompute_mon_lp_static <- function(draws_mon_sub, frame_tu, covar_tbl,
                                     mon_covar_map = character(0),
                                     drop_region = FALSE) {
  U         <- nrow(frame_tu)
  gids      <- as.character(frame_tu$gid)
  clim_cols <- intersect(c("mean_precip", "dev_precip",
                            "mean_temp_daymean", "dev_temp_daymean"),
                          names(covar_tbl))

  # Offsets, read off the bundle (see block comment above).
  off_pop  <- if (isTRUE(draws_mon_sub$offset_pop      %||% TRUE))
                log(frame_tu$pop)       else rep(0.0, U)
  # offset(log(base_rate_mon)) -- mon's OWN anchor column.
  off_anch <- .anchor_log_offset(draws_mon_sub, frame_tu, NULL,
                                 "precompute_mon_lp_static")
  static_offset <- off_pop + off_anch
  use_logoff    <- isTRUE(draws_mon_sub$offset_logoffset %||% FALSE)

  # Normalise the term -> column map: unnamed entries map to themselves.
  map_nm <- names(mon_covar_map) %||% rep("", length(mon_covar_map))
  map_nm <- ifelse(nzchar(map_nm), map_nm, as.character(mon_covar_map))
  keep   <- as.character(mon_covar_map) %in% names(covar_tbl)
  mon_covar_map <- as.character(mon_covar_map)[keep]
  map_nm        <- map_nm[keep]

  # admin_level dummy: static per unit (a frame attribute, not a panel column).
  # Mirrors the v1 kernel's admin_adm2_vec resolution exactly, including the
  # CANONICAL key -- the mon fit's admin_level factor is re-keyed to GADM
  # depth before it is levelled (.frameAdminLevelCanonical(), simulator_core.R).
  adm2_vec <- if ("admin_level" %in% names(frame_tu))
    as.numeric(.frameAdminLevelCanonical(frame_tu) == "ADM2") else NULL

  lp_out <- rep(NA_real_, nrow(covar_tbl))
  times  <- sort(unique(covar_tbl$time))
  # pop_year mode: a `.pop_t` column (this row's month's
  # population, set by build_spatpomp_tu_c()) replaces the static `pop` in the
  # log(pop) offset and in evaluate_lp()'s `log_pop_c`. Absent by default.
  pop_tv <- ".pop_t" %in% names(covar_tbl)

  for (t in times) {
    mask_t   <- covar_tbl$time == t
    rows_t   <- covar_tbl[mask_t, , drop = FALSE]
    perm     <- match(gids, rows_t$unit)
    rows_ord <- rows_t[perm, , drop = FALSE]
    frame_t  <- frame_tu
    off_t    <- static_offset
    if (pop_tv) {
      frame_t$pop <- rows_ord$.pop_t
      off_t <- (if (isTRUE(draws_mon_sub$offset_pop %||% TRUE)) log(frame_t$pop)
                else rep(0.0, U)) + off_anch
    }

    xrow <- setNames(lapply(clim_cols, function(col) rows_ord[[col]]), clim_cols)
    xrow$t_start_seas <- rows_ord$t_start_seas
    for (i in seq_along(mon_covar_map))
      xrow[[map_nm[i]]] <- rows_ord[[mon_covar_map[i]]]
    if (!is.null(adm2_vec)) xrow[["admin_levelADM2"]] <- adm2_vec

    lp_t <- evaluate_lp(
      draws_mon_sub, s = 1L, xrow,
      o      = rep(NA_integer_, U),
      o_ct   = rep(NA_integer_, U),
      o_tu   = rep(NA_integer_, U),
      u_days = rep(NA_real_,    U),
      frame_tu    = frame_t,
      basehaz_val = 0,
      drop_region = drop_region
    )
    lp_t <- lp_t + off_t
    if (use_logoff) lp_t <- lp_t + rows_ord$log_offset

    inv_perm <- match(rows_t$unit, gids)
    lp_out[mask_t] <- lp_t[inv_perm]
  }

  lp_out
}

#' Did precompute_mon_lp_static() actually fold `log_offset` into lp_static?
#'
#' The generated Csnippet needs to know, because its onset-month exposure
#' correction is only valid if there is an offset there to correct (see
#' .ONSET_EXPOSURE_FRAC, simulator_core.R). This reads the SAME bundle key
#' precompute_mon_lp_static() reads at :789, so the two cannot disagree -- the
#' failure mode being guarded against is the recurring "config computed but not
#' threaded to its consumer" pattern, where the fold and the correction get
#' their answer from two different places.
#'
#' `offset_logoffset` is a property of the mon FIT, so it is identical across
#' groups; the grouped branch reads the first group's bundle and asserts the
#' rest agree rather than silently taking one group's word for all of them.
#'
#' @param draws_mon_sub    ungrouped 1-draw mon bundle, or NULL.
#' @param grouping_mon_sub grouped mon subset (`draws_sub_by_group`), or NULL.
#' @return TRUE/FALSE; FALSE when there is no mon bundle at all.
.mon_logoffset_flag <- function(draws_mon_sub, grouping_mon_sub = NULL) {
  bundles <- if (!is.null(grouping_mon_sub))
    grouping_mon_sub$draws_sub_by_group %||% list()
  else if (!is.null(draws_mon_sub)) list(draws_mon_sub)
  else list()

  if (!length(bundles)) return(FALSE)
  flags <- vapply(bundles,
                  function(b) isTRUE((b %||% list())$offset_logoffset %||% FALSE),
                  logical(1))
  if (length(unique(flags)) > 1L)
    stop(".mon_logoffset_flag: the per-group mon bundles disagree on ",
         "`offset_logoffset` (", sum(flags), " of ", length(flags), " TRUE). ",
         "That is a property of the fitted formula and must be uniform; ",
         "re-extract the bundles with the current extractor.")
  flags[[1L]]
}

#' Sibling of precompute_ttlo_coefs for the monthly-case model.
precompute_mon_ttlo_coefs <- function(draws_mon_sub, frame_tu,
                                      drop_region = FALSE) {
  precompute_ttlo_coefs(draws_mon_sub, frame_tu, drop_region)
}

#' Sibling of precompute_spillover_coefs for the monthly-case model.
precompute_mon_spillover_coefs <- function(draws_mon_sub, frame_tu,
                                           drop_region = FALSE) {
  precompute_spillover_coefs(draws_mon_sub, frame_tu, drop_region)
}

#' Precompute static per-unit NegBin dispersion (phi) for the monthly model.
precompute_mon_phi <- function(draws_mon_sub, frame_tu, drop_region = FALSE) {
  compute_phi_vec(draws_mon_sub, 1L, frame_tu, drop_region = drop_region)
}

#' Precompute the latent months-since-onset (mso) cross-basis coefficient table.
#'
#' The `mso` cross-basis carries per-unit hierarchical RE slopes
#' (monthly.re_adm1_slope defaults to TRUE), so unlike the duration basehaz it
#' cannot collapse to a global constant array. It CAN collapse over the basis
#' index k, because the basis itself is a build-time lookup table
#' (build_mso_basis()) evaluated on the integers 0..mso_max:
#'
#'   entry[u, j] = sum_k (fe[["mso_cb<k>"]] + ns1("mso_cb<k>")[u]) * B[j+1, k]
#'
#' so the Csnippet does ONE covariate read indexed by the clamped latent
#' mv_new[u] -- zero basis arithmetic in C. Structurally a clone of
#' precompute_ttlo_coefs(), with (mso_max + 1) bins instead of 7.
#'
#' Clamping to `mso_max` happens twice, by construction: build_mso_basis()
#' evaluates only 0..mso_max, and the Csnippet clamps its index into that
#' range. This is mandatory, not defensive -- bs(degree = 3) extrapolates as an
#' unbounded cubic past its boundary knots (0.8-1.0 units past the right knot
#' can inflate mu 26x; see config_defaults.yml).
#'
#' @param draws_mon_sub 1-draw monthly draws bundle.
#' @param frame_tu      location frame for the TU.
#' @param mso_basis     list(B, mso_max) from build_mso_basis(), or NULL.
#' @return data.frame with columns mon_mso_coef_j0..j<mso_max> (U rows), or
#'   NULL when the fit carries no mso cross-basis.
precompute_mon_mso_coefs <- function(draws_mon_sub, frame_tu, mso_basis,
                                     drop_region = FALSE) {
  if (is.null(mso_basis)) return(NULL)
  U       <- nrow(frame_tu)
  fe      <- draws_mon_sub$fe[1L, , drop = TRUE]
  B       <- mso_basis$B
  mso_max <- as.integer(mso_basis$mso_max)
  cn      <- colnames(B)

  # Per-unit total coefficient (FE + summed nested RE slopes) for each basis
  # column k, then collapse against the basis rows.
  coef_k <- lapply(seq_along(cn), function(k) {
    fe_val <- if (cn[k] %in% names(fe)) fe[[cn[k]]] else 0.0
    fe_val + .ns1(draws_mon_sub, cn[k], frame_tu, drop_region)
  })

  out <- vector("list", mso_max + 1L)
  names(out) <- paste0("mon_mso_coef_j", 0:mso_max)
  for (j in 0:mso_max) {
    acc <- rep(0.0, U)
    for (k in seq_along(cn)) acc <- acc + coef_k[[k]] * B[j + 1L, k]
    out[[j + 1L]] <- acc
  }

  as.data.frame(out)
}

#' Group-aware dispatch for monthly-case precomputes.
#'
#' Ungrouped: standalone calls for lp, ttlo, spill, phi, mso.
#' Grouped: per-group dispatch with drop_region = TRUE, stitched to full-frame
#' order exactly like the occurrence/duration siblings.
#'
#' Groups may carry DIFFERENT mso ranges (each group's `mon` fit has its own
#' cb_meta), but the Csnippet needs one rectangular table. The stitched table
#' is therefore widened to the global max and each group's rows are padded by
#' repeating that group's own clamped last column -- i.e. every group keeps its
#' own clamp point, and the shared width is purely a layout detail.
#'
#' @param mso_basis          list(B, mso_max) -- ungrouped mode only.
#' @param mso_basis_by_group named list of per-group mso bases -- grouped only.
#' @return list(lp, ttlo, spill, phi, mso = <df or NULL>).
precompute_mon_lp_components <- function(frame_tu, covar_tbl, mon_covar_map,
                                          draws_mon_sub,
                                          grouping_mon_sub   = NULL,
                                          mso_basis          = NULL,
                                          mso_basis_by_group = NULL) {
  if (is.null(grouping_mon_sub)) {
    return(list(
      lp    = precompute_mon_lp_static(draws_mon_sub, frame_tu, covar_tbl,
                                       mon_covar_map),
      ttlo  = precompute_mon_ttlo_coefs(draws_mon_sub, frame_tu),
      spill = precompute_mon_spillover_coefs(draws_mon_sub, frame_tu),
      phi   = precompute_mon_phi(draws_mon_sub, frame_tu),
      mso   = precompute_mon_mso_coefs(draws_mon_sub, frame_tu, mso_basis)
    ))
  }

  if (!("group" %in% names(frame_tu)))
    stop("precompute_mon_lp_components: grouped mode requires `group` column on frame_tu.")

  st <- .stitch_grouped_lp(
    frame_tu, covar_tbl, grouping_mon_sub$draws_sub_by_group,
    lp_fn    = function(d, f, cv) precompute_mon_lp_static(d, f, cv, mon_covar_map,
                                                           drop_region = TRUE),
    ttlo_fn  = precompute_mon_ttlo_coefs, spill_fn = precompute_mon_spillover_coefs,
    missing_msg = "precompute_mon_lp_components: no monthly draws subset for group ",
    extra_fn = function(d, f, g) list(
      phi = precompute_mon_phi(d, f, drop_region = TRUE),
      mso = precompute_mon_mso_coefs(d, f, (mso_basis_by_group %||% list())[[g]],
                                     drop_region = TRUE)))
  U <- nrow(frame_tu); groups <- st$groups; grp_vec <- st$grp_vec
  lp_full <- st$lp; ttlo_full <- st$ttlo; spill_full <- st$spill
  phi_full <- rep(NA_real_, U)
  for (g in groups) phi_full[grp_vec == g] <- st$extra[[g]]$phi
  mso_by_g <- lapply(st$extra, `[[`, "mso")

  # Stitch the per-group mso tables into one rectangular U x (Jmax+1) frame,
  # padding each group's short rows with its own last (clamped) column.
  mso_full <- NULL
  if (any(!vapply(mso_by_g, is.null, logical(1)))) {
    ncols    <- vapply(mso_by_g, function(d) if (is.null(d)) 0L else ncol(d), integer(1))
    jmax     <- max(ncols) - 1L
    mso_full <- as.data.frame(setNames(
      rep(list(rep(0.0, U)), jmax + 1L), paste0("mon_mso_coef_j", 0:jmax)))
    for (g in groups) {
      d <- mso_by_g[[g]]
      if (is.null(d)) next     # group with no mso basis: contributes 0
      row_mask_g <- grp_vec == g
      for (j in 0:jmax)
        mso_full[row_mask_g, j + 1L] <- d[[min(j + 1L, ncol(d))]]
    }
  }

  list(lp = lp_full, ttlo = ttlo_full, spill = spill_full,
       phi = phi_full, mso = mso_full)
}

# ---------------------------------------------------------------------------
# Augment covariate table with precomputed LP components + init state
# ---------------------------------------------------------------------------

#' Augment the base covariate table with all Csnippet-needed columns.
#'
#' Adds:
#'   lp_static        — time-varying precomputed LP (one value per covar row)
#'   ttlo_coef_b0..6  — per-unit TTLO contributions (constant in time)
#'   spill_{o,oct,otu}_{t,f} — per-unit spillover contributions
#'   u_days_init, o_init, o_ct_init, o_tu_init — warmup-state init for rinit
#'
#' @param covar_tbl      base covariate table (from build_covar_for_tu).
#' @param lp_static_vec  numeric vector length nrow(covar_tbl).
#' @param ttlo_coefs_df  data.frame from precompute_ttlo_coefs (U rows).
#' @param spill_coefs_df data.frame from precompute_spillover_coefs (U rows).
#' @param warmup_state   list(u, o, o_ct, o_tu) from build_warmup_state.
#' @param frame_tu       location frame for the TU.
#' @param size_phi_vec   optional per-unit mon NegBin dispersion; emitted as the
#'                       `size_phi_init` column the rinit Csnippet reads to seed
#'                       the `size_phi` state slot ST(8,u). The name is the
#'                       historical state-slot name (see .SP_UNIT_STATENAMES).
#' @param mon_mu_cap_vec optional per-unit ceiling on the mon mean (cap x pop),
#'                       emitted as the `mon_mu_cap` column.
#' @return augmented covariate tibble.
build_covar_for_tu_c <- function(covar_tbl, lp_static_vec,
                                  ttlo_coefs_df, spill_coefs_df,
                                  warmup_state, frame_tu,
                                  size_phi_vec        = NULL,
                                  dur_lp_static_vec   = NULL,
                                  dur_ttlo_coefs_df   = NULL,
                                  dur_spill_coefs_df  = NULL,
                                  mon_lp_static_vec   = NULL,
                                  mon_ttlo_coefs_df   = NULL,
                                  mon_spill_coefs_df  = NULL,
                                  mon_mso_coefs_df    = NULL,
                                  nbr_coefs_vec       = NULL,
                                  mon_mu_cap_vec      = NULL,
                                  mon_mu_cap_by_row   = FALSE) {
  gids     <- as.character(frame_tu$gid)
  sentinel <- .SP_OTU_SENTINEL

  # 1. lp_static: time-varying, already in covar_tbl row order.
  covar_tbl$lp_static <- lp_static_vec

  # 2. Per-unit constant columns: TTLO + spillover.
  #    Join by unit so replication across time steps is handled automatically.
  ttlo_aug  <- cbind(data.frame(unit = gids, stringsAsFactors = FALSE),
                      ttlo_coefs_df)
  spill_aug <- cbind(data.frame(unit = gids, stringsAsFactors = FALSE),
                      spill_coefs_df)

  covar_tbl <- dplyr::left_join(covar_tbl, ttlo_aug,  by = "unit")
  covar_tbl <- dplyr::left_join(covar_tbl, spill_aug, by = "unit")

  # 2b. Neighbour-contagion coefficient beta_nbr(u), one extra per-unit covar
  #     column. Only present when use_nbr (nbr_coefs_vec supplied); the Csnippet
  #     unpacks `nbr_coef` and multiplies it by the per-step kernel matvec. No
  #     new STATE slot — neighbour pressure is derived from the `o` state.
  if (!is.null(nbr_coefs_vec)) {
    nbr_aug <- data.frame(unit = gids, nbr_coef = as.numeric(nbr_coefs_vec),
                          stringsAsFactors = FALSE)
    covar_tbl <- dplyr::left_join(covar_tbl, nbr_aug, by = "unit")
  }

  # 3. Warmup-state init columns (read by rinit Csnippet at t0).
  init_aug <- data.frame(
    unit        = gids,
    u_days_init = as.numeric(warmup_state$u),
    o_init      = as.numeric(warmup_state$o),
    o_ct_init   = as.numeric(warmup_state$o_ct),
    o_tu_init   = ifelse(is.na(warmup_state$o_tu),
                          sentinel,
                          as.numeric(warmup_state$o_tu)),
    active_init = as.numeric(warmup_state$active %||% rep(0L, length(gids))),
    mso_init    = as.numeric(warmup_state$mso    %||% rep(0L, length(gids))),
    stringsAsFactors = FALSE
  )
  # Detection-flag seed: present only under covariate_feed = "detected"
  # (build_spatpomp_tu_c() sets warmup_state$det), so a latent run's covariate
  # table -- and therefore every cov_bases offset -- is unchanged.
  if (!is.null(warmup_state$det)) init_aug$det_init <- as.numeric(warmup_state$det)

  covar_tbl <- dplyr::left_join(covar_tbl, init_aug, by = "unit")

  # 4. Duration (cessation-hazard model): lp_static, TTLO, spillover.
  #    Column names carry "dur_" prefix. No dur_phi (duration has no NegBin
  #    dispersion) and no basehaz table here -- the mso-band hazard table has no
  #    per-unit RE (FE-only) and is embedded as a static C constant array by
  #    the rprocess Csnippet generator, not joined as a per-unit covariate.
  if (!is.null(dur_lp_static_vec)) {
    covar_tbl$dur_lp_static <- dur_lp_static_vec
  }
  if (!is.null(dur_ttlo_coefs_df)) {
    dur_ttlo_aug <- cbind(
      data.frame(unit = gids, stringsAsFactors = FALSE),
      setNames(dur_ttlo_coefs_df, paste0("dur_", names(dur_ttlo_coefs_df)))
    )
    covar_tbl <- dplyr::left_join(covar_tbl, dur_ttlo_aug, by = "unit")
  }
  if (!is.null(dur_spill_coefs_df)) {
    dur_spill_aug <- cbind(
      data.frame(unit = gids, stringsAsFactors = FALSE),
      setNames(dur_spill_coefs_df, paste0("dur_", names(dur_spill_coefs_df)))
    )
    covar_tbl <- dplyr::left_join(covar_tbl, dur_spill_aug, by = "unit")
  }

  # 5. Monthly cases (`mon`): lp_static, TTLO, spillover, the latent-mso
  #    coefficient table, and phi. Column names carry a "mon_" prefix, except
  #    phi, which is the `size_phi_init` column the rinit Csnippet reads to
  #    seed ST(8,u). mon_mso_coef_j* is per-unit constant (the latent mso
  #    index that selects among its columns is a STATE, resolved at Csnippet
  #    runtime).
  if (!is.null(mon_lp_static_vec)) {
    covar_tbl$mon_lp_static <- mon_lp_static_vec
  }
  if (!is.null(mon_ttlo_coefs_df)) {
    mon_ttlo_aug <- cbind(
      data.frame(unit = gids, stringsAsFactors = FALSE),
      setNames(mon_ttlo_coefs_df, paste0("mon_", names(mon_ttlo_coefs_df)))
    )
    covar_tbl <- dplyr::left_join(covar_tbl, mon_ttlo_aug, by = "unit")
  }
  if (!is.null(mon_spill_coefs_df)) {
    mon_spill_aug <- cbind(
      data.frame(unit = gids, stringsAsFactors = FALSE),
      setNames(mon_spill_coefs_df, paste0("mon_", names(mon_spill_coefs_df)))
    )
    covar_tbl <- dplyr::left_join(covar_tbl, mon_spill_aug, by = "unit")
  }
  if (!is.null(mon_mso_coefs_df)) {
    mon_mso_aug <- cbind(data.frame(unit = gids, stringsAsFactors = FALSE),
                          mon_mso_coefs_df)
    covar_tbl <- dplyr::left_join(covar_tbl, mon_mso_aug, by = "unit")
  }
  if (!is.null(size_phi_vec)) {
    phi_aug <- data.frame(unit = gids, size_phi_init = as.numeric(size_phi_vec),
                           stringsAsFactors = FALSE)
    covar_tbl <- dplyr::left_join(covar_tbl, phi_aug, by = "unit")
  }
  # Monthly incidence ceiling, cap x pop per unit (monthlyMuCap()); read by
  # mon_mu_block as fmin(exp(mn_lp), mon_mu_cap). Absent unless a finite cap
  # is set, so uncapped builds emit the same covariate table as before.
  # mon_mu_cap_by_row (pop_year mode): the ceiling is already
  # one value per covar_tbl row (cap x the row's month's population), in the
  # row order every join above preserves; the column lands in the same place.
  if (!is.null(mon_mu_cap_vec) && isTRUE(mon_mu_cap_by_row)) {
    if (length(mon_mu_cap_vec) != nrow(covar_tbl))
      stop("build_covar_for_tu_c: mon_mu_cap_by_row needs one ceiling per covariate row")
    covar_tbl$mon_mu_cap <- as.numeric(mon_mu_cap_vec)
  } else if (!is.null(mon_mu_cap_vec)) {
    cap_aug <- data.frame(unit = gids, mon_mu_cap = as.numeric(mon_mu_cap_vec),
                          stringsAsFactors = FALSE)
    covar_tbl <- dplyr::left_join(covar_tbl, cap_aug, by = "unit")
  }

  covar_tbl
}

#' De-digit covariate column names to avoid spatPomp's internal flatten collision.
#'
#' spatPomp v1.1.0 flattens its covariate table with paste0(covariate_name,
#' unit_index), where unit_index is the bare positional integer 1..U (see
#' construct_spatPomp()'s .local() body; docstring at spatpomp_core.R:449-452).
#' That concatenation is ambiguous whenever one covariate name is another name
#' plus digits -- concretely `mon_mso_coef_j1` vs `mon_mso_coef_j11`
#' (precompute_mon_mso_coefs() emits j0..j<mso_max>, and the deployed mon fits
#' have mso_max 12-26), which collides for every TU with U >= 11 and makes
#' tidyr::spread() throw "Each row of output must be identified by a unique
#' combination of keys". Appending a non-digit sentinel to any digit-ending
#' name makes the concatenation injective for any U: the maximal trailing
#' digit-run of paste0(name, ui) then always recovers ui exactly, and the
#' remainder always recovers name exactly.
#'
#' Offset-neutral: spatPomp re-imposes OUR column order right after the
#' spread (`pomp_covar[, c(timename, cov_col_order)]` with cov_col_order built
#' by looping unit_covarnames x seq_len(U)), so renaming changes no column's
#' position and downstream `cov_bases = (k-1) * U` offsets stay valid.
#' Numerically bit-identical for runs that already work.
#'
#' Must be applied to `covar_tbl` only AFTER every name-based consumer of its
#' columns (cov_bases, init_bases, the generate_*_csnippet() calls) has run --
#' see call site in build_spatpomp_tu_c().
#'
#' @param covar_tbl long covariate table (columns unit, time, ...).
#' @param exclude   column names never renamed (the join keys).
#' @return covar_tbl with digit-ending covariate columns renamed.
.sp_dedigit_covar_names <- function(covar_tbl, exclude = c("unit", "time")) {
  nm  <- names(covar_tbl)
  hit <- grepl("[0-9]$", nm) & !(nm %in% exclude)
  names(covar_tbl)[hit] <- paste0(nm[hit], "_")
  # spatPomp validates uniqueness of `data` names only, never of covariate
  # names, so a duplicate here would silently shift every __covars offset
  # instead of erroring. Cannot happen with the current column set (no
  # "<digit-ending>_" twin exists), but the check is one line and turns a
  # silent-wrong-answer class into a loud failure.
  if (anyDuplicated(names(covar_tbl)))
    stop(".sp_dedigit_covar_names: de-digiting produced duplicate covariate ",
         "column names (", paste(unique(names(covar_tbl)[duplicated(names(covar_tbl))]),
                                  collapse = ", "), "); ",
         "spatPomp would silently mis-index __covars. Rename the offending ",
         "source column instead.")
  covar_tbl
}

# ---------------------------------------------------------------------------
# Membership arrays for topology (embedded as static C constants)
# ---------------------------------------------------------------------------

#' Build 0-based membership arrays for the rprocess Csnippet.
#'
#' @return named list:
#'   adm1_id   — 0-based ADM1 index per unit
#'   ct_id     — 0-based country index per unit
#'   tu_id     — 0-based TU index per unit (-1 for NA-TU units)
#'   adm1_to_ct — country index per ADM1
#'   ct_to_tu  — TU index per country (-1 for NA-TU countries)
#'   n_adm1, n_ct, n_tu  — level counts
make_membership_arrays_c <- function(frame_tu) {
  adm1_lev <- unique(frame_tu$adm1_unit)
  ct_lev   <- unique(frame_tu$country)
  tu_lev   <- unique(stats::na.omit(frame_tu$tu))

  adm1_id <- match(frame_tu$adm1_unit, adm1_lev) - 1L
  ct_id   <- match(frame_tu$country,   ct_lev)   - 1L
  tu_id_r <- match(frame_tu$tu,        tu_lev)
  tu_id   <- ifelse(is.na(tu_id_r), -1L, tu_id_r - 1L)

  # ADM1 → country (one distinct mapping)
  adm1_ct_map <- frame_tu %>%
    dplyr::distinct(adm1_unit, country)
  adm1_to_ct <- match(
    adm1_ct_map$country[match(adm1_lev, adm1_ct_map$adm1_unit)],
    ct_lev
  ) - 1L

  # Country → TU (-1 for NA-TU countries)
  ct_tu_map  <- frame_tu %>% dplyr::distinct(country, tu)
  ct_tu_vals <- ct_tu_map$tu[match(ct_lev, ct_tu_map$country)]
  ct_to_tu   <- ifelse(is.na(ct_tu_vals), -1L, match(ct_tu_vals, tu_lev) - 1L)

  list(
    adm1_id    = adm1_id,
    ct_id      = ct_id,
    tu_id      = tu_id,
    adm1_to_ct = adm1_to_ct,
    ct_to_tu   = ct_to_tu,
    n_adm1     = length(adm1_lev),
    n_ct       = length(ct_lev),
    n_tu       = length(tu_lev)
  )
}

# ---------------------------------------------------------------------------
# Csnippet generators — rinit
# ---------------------------------------------------------------------------

#' Generate the rinit Csnippet for N units.
#'
#' Reads warmup state from per-unit covariate columns (u_days_init{u}, etc.)
#' and initialises the state vector.  draw_idx is set to fixed_draw (always 1
#' in the N-independent-runs design).
#'
#' @param U integer; number of units.
#' @param fixed_draw integer; draw index to embed (default 1L).
#' @return pomp::Csnippet object.
generate_rinit_csnippet <- function(U, init_bases, fixed_draw = 1L,
                                     size_phi_base = NULL, zero_slots = integer(0)) {
  # zero_slots: 0-based state indices seeded 0.0 (the WHO scratch slots cm/ccm,
  # appended after every other slot; empty for every other measurement, which
  # keeps the emitted C unchanged).
  # det(11) is seeded only when init_bases carries `det_init`, which
  # build_covar_for_tu_c() adds only under covariate_feed = "detected".
  det_init_base <- if ("det_init" %in% names(init_bases))
    as.integer(init_bases[["det_init"]]) else NULL
  # init_bases: named integer vector giving the absolute (0-based) position of
  #   each init covariate's FIRST-UNIT copy in the covars array passed to rinit.
  #   Accessed as covars[base + u] (0-based u), bypassing covindex.
  #
  # State index mapping (fixed by .SP_UNIT_STATENAMES order, 11 states total):
  #   k=0 u_days  k=1 o  k=2 o_ct  k=3 o_tu  k=4 Yout
  #   k=5 draw_idx  k=6 p  k=7 mu  k=8 size_phi  k=9 active  k=10 mso
  # spatPomp uses __stateindex[k*U+u] for the flat position of state k, unit u.
  #
  # size_phi_base: integer base offset in __covars for the size_phi_init column
  #   (the mon dispersion), or NULL (Bernoulli path → seed 0.0).
  #
  # active_init/mso_init are always present in init_bases (build_covar_for_tu_c
  # unconditionally adds them, defaulting to 0 via warmup_state$active %||% 0)
  # so they are written unconditionally here, same convention as u_days/o/o_ct/o_tu.
  # Legacy (!use_dur) runs: rprocess never writes ST(9,u)/ST(10,u) again, so
  # these seeded-0 values are simply held in place for the whole run (pomp's
  # __x state buffer persists in-place across steps -- same precedent as
  # size_phi/ST(8,u), which rprocess also never unconditionally rewrites).
  fixed_d_c      <- sprintf("%.1f", as.double(fixed_draw))
  ud_init_base   <- as.integer(init_bases[["u_days_init"]])
  o_init_base    <- as.integer(init_bases[["o_init"]])
  oct_init_base  <- as.integer(init_bases[["o_ct_init"]])
  otu_init_base  <- as.integer(init_bases[["o_tu_init"]])
  active_init_base <- as.integer(init_bases[["active_init"]])
  mso_init_base    <- as.integer(init_bases[["mso_init"]])
  # k * U offsets (R-time literals embedded in the generated C string).
  sv_offsets <- as.integer(0:10) * as.integer(U)   # 0*U, 1*U, ..., 10*U
  size_phi_line <- if (!is.null(size_phi_base)) {
    paste0("    __x[__stateindex[", sv_offsets[9], " + u]] = __covars[",
           as.integer(size_phi_base), " + u];\n")
  } else {
    paste0("    __x[__stateindex[", sv_offsets[9], " + u]] = 0.0;\n")
  }
  code <- paste0(
    "\n/* --- rinit: __x[__stateindex[k*U+u]] for state var k, unit u (O(1) code) --- */\n",
    "/* .SP_UNIT_STATENAMES order: u_days(0),o(1),o_ct(2),o_tu(3),Yout(4),draw_idx(5),p(6),mu(7),size_phi(8),active(9),mso(10) */\n",
    "{ int u;\n",
    "  for (u = 0; u < ", U, "; u++) {\n",
    "    __x[__stateindex[", sv_offsets[1], " + u]] = __covars[", ud_init_base,  " + u];\n",
    "    __x[__stateindex[", sv_offsets[2], " + u]] = __covars[", o_init_base,   " + u];\n",
    "    __x[__stateindex[", sv_offsets[3], " + u]] = __covars[", oct_init_base, " + u];\n",
    "    __x[__stateindex[", sv_offsets[4], " + u]] = __covars[", otu_init_base, " + u];\n",
    "    __x[__stateindex[", sv_offsets[5], " + u]] = 0.0;\n",       # Yout = 0
    "    __x[__stateindex[", sv_offsets[6], " + u]] = ", fixed_d_c, ";\n",  # draw_idx
    "    __x[__stateindex[", sv_offsets[7], " + u]] = 0.5;\n",       # p = 0.5
    "    __x[__stateindex[", sv_offsets[8], " + u]] = 0.0;\n",       # mu = 0 (overwritten by rprocess)
    size_phi_line,                                                     # size_phi: from covars or 0
    "    __x[__stateindex[", sv_offsets[10], " + u]] = __covars[", active_init_base, " + u];\n",
    "    __x[__stateindex[", sv_offsets[11], " + u]] = __covars[", mso_init_base,    " + u];\n",
    if (!is.null(det_init_base))
      paste0("    __x[__stateindex[", 11L * as.integer(U), " + u]] = __covars[", det_init_base, " + u];\n")
    else "",
    # if-guarded: paste0() recycles a zero-length vector to "" and would emit
    # a stray `__x[__stateindex[ + u]] = 0.0;` (valid C that zeroes u_days).
    if (length(zero_slots))
      paste0("    __x[__stateindex[", as.integer(zero_slots) * as.integer(U), " + u]] = 0.0;\n",
             collapse = "")
    else "",
    "  }\n",
    "}\n"
  )
  pomp::Csnippet(code)
}

# ---------------------------------------------------------------------------
# Csnippet generators — rprocess
# ---------------------------------------------------------------------------

#' Generate the rprocess Csnippet for N units.
#'
#' The generated C code:
#'   1. Embeds topology as static const int arrays (initialised once).
#'   2. Unpacks current state via __stateindex loop (O(1) code).
#'   3. Unpacks covariates via __covarindex loop (O(1) code).
#'   4. For each unit: computes LP = lp_static + TTLO bin coef +
#'      spillover coef, then hazard p and draws Yout.
#'   5. Computes ADM1 → country → TU aggregates for lag update.
#'   6. Writes updated state back via __stateindex loop (O(1) code).
#'
#' @param U          integer; number of units.
#' @param memb       list from make_membership_arrays_c().
#' @param fixed_draw integer; draw index to write into draw_idx state (default 1L).
#' @param force_o_tu logical; when TRUE, ST(3,u) (next-step o_tu) is written as
#'   an unconditional 1.0 for every valid-TU unit instead of the computed
#'   leave-one-country-out expression. Under partition_by=country a TU is a
#'   single country, so the computed expression is always 0 (see
#'   .fold_observed_otu_into_lp docstring) — this flag is the simulated-period
#'   counterpart, used when there is no observed panel to fold in (forecast /
#'   counterfactual), or when the caller prefers the simpler always-active
#'   assumption over the observed fold during smoothing. NA-TU units (TUID<0)
#'   still get the -1.0 sentinel regardless. Default FALSE preserves prior
#'   behaviour exactly.
#' @param use_dur   logical; opt-in for the duration-aware (cessation
#'   hazard) path, mirroring simulate_one_tu()'s `use_dur` gate. FALSE (default)
#'   emits no active/mso code at all — states 9/10 are never written, so they
#'   hold their rinit-seeded 0 for the whole run, as size_phi/ST(8,u) do.
#' @param dur_bh    list(kind, breaks, table) from precompute_dur_basehaz_table();
#'   ungrouped duration basehaz lookup. NULL when duration is grouped or absent.
#' @param dur_bh_by_group named list of per-group dur_bh lists (grouped duration
#'   mode only); all groups must share identical kind/breaks (enforced by
#'   precompute_dur_lp_components()). NULL when ungrouped or duration absent.
#' @param dur_grp_id integer vector length U (0-based, in the SAME unit order as
#'   `memb`/frame_sorted — NOT frame_tu row order) mapping each unit to its
#'   row in dur_bh_by_group's tables. Required iff dur_bh_by_group is supplied.
#' @param dur_groups character vector of group labels giving the row order of
#'   dur_bh_by_group's tables (i.e. dur_grp_id indexes into this). Required iff
#'   dur_bh_by_group is supplied.
#' @param use_nbr logical; emit the distance-weighted neighbour-contagion term.
#'   When TRUE, the preamble embeds the WNBR[N_ADM1*N_ADM1] kernel, phase 1
#'   reconstructs the lag-1 ADM1-active vector from the `o` state and adds
#'   `bnbr[u] * (WNBR %*% Aa1)[ADID[u]]` to the occurrence LP. Gated ONLY by this
#'   flag -- independent of force_o_tu / drop_country_tu_spillover. FALSE
#'   (default) emits byte-identical C to the pre-nbr code.
#' @param nbr_kernel dense ADM1 x ADM1 kernel (nrow = ncol = memb$n_adm1), in the
#'   SAME 0-based ADM1 order as memb$adm1_id. Required iff use_nbr. From
#'   build_partition_nbr_kernel(frame_sorted, ...).
#' @param mon_mso_ncol integer; number of mon_mso_coef_j* covariate columns
#'   (= mso_max + 1). 0 when the mon fit carries no mso cross-basis, in which
#'   case the generated block simply omits the mso contribution. Ignored unless
#'   the measurement selects the mon path.
#' @param mon_mu_cap logical; TRUE (a finite simulator.incidence_cap on the mon
#'   path) caps mu at the per-unit `mon_mu_cap` covariate, fmin(exp(mn_lp), cap).
#'   FALSE emits the uncapped expression unchanged.
#' @param onset_season logical; TRUE (simulator.onset_eligibility: season)
#'   gates the onset draw on no u_days-resetting onset in the unit's current
#'   season, reading the `seas_k` covariate as sk[u] (seasonOnsetOpen() in
#'   simulator_core.R is the R side). Requires use_dur. FALSE emits the
#'   previous C unchanged.
#' @param cess_half logical; TRUE (simulator.cessation_exposure: half) emits
#' mon_mu_block AFTER dur_cessation_block, reading a snapshot
#'   of this month's pre-cessation active / mso (av_act / mv_act) and halving
#'   the exposure of a unit that ceased this step. No-op without the mon path.
#'   FALSE emits the previous C unchanged.
#' @return pomp::Csnippet object.
generate_rprocess_csnippet <- function(U, memb, cov_bases, fixed_draw = 1L,
                                        measurement = "occurrence_bernoulli",
                                        force_o_tu = FALSE,
                                        use_dur = FALSE,
                                        dur_bh = NULL,
                                        dur_bh_by_group = NULL,
                                        dur_grp_id = NULL,
                                        dur_groups = NULL,
                                        mon_mso_ncol = 0L,
                                        mon_use_logoffset = FALSE,
                                        mon_mu_cap = FALSE,
                                        drop_country_tu_spillover = FALSE,
                                        use_nbr = FALSE,
                                        nbr_kernel = NULL,
                                        covariate_feed = "latent",
                                        delta_fn = .SP_STATE_DELTA_FN,
                                        delta_fp = .SP_STATE_DELTA_FP,
                                        who_slots = NULL,
                                        ttlo_none = FALSE,
                                        onset_season = FALSE,
                                        cess_half = FALSE) {
  # ttlo_none: the occurrence fit carries the TTLO level
  # "none" -- an 8th occurrence bin (ttlo_coef_b7) and the sentinel branch in
  # its bin chain. FALSE emits exactly the previous C.
  # who_slots: NULL, or c(cm = k, ccm = k + 1), the 0-based indices of the WHO
  # scratch slots (`_who` measurement). Non-NULL emits the simulated-cases block
  # after mon_mu_block; NULL emits exactly the previous C.
  # Detected-activity covariate feed: C twin of simulator_core.R::
  # update_detection(); needs the duration track's `active`.
  use_det <- identical(covariate_feed, "detected")
  if (use_det && !isTRUE(use_dur))
    stop("generate_rprocess_csnippet: covariate_feed = 'detected' requires use_dur.")

  # Neighbour-contagion build-time contract: the WNBR matvec needs the dense
  # ADM1 x ADM1 kernel, indexed in the SAME 0-based ADM1 order as memb$adm1_id
  # (i.e. unique(frame_sorted$adm1_unit)). Caller builds it via
  # build_partition_nbr_kernel(frame_sorted, ...). Fail loud here rather than
  # emit C that reads a mis-sized/absent WNBR.
  if (isTRUE(use_nbr)) {
    if (is.null(nbr_kernel))
      stop("generate_rprocess_csnippet: use_nbr=TRUE requires nbr_kernel.")
    if (!is.matrix(nbr_kernel) || nrow(nbr_kernel) != memb$n_adm1 ||
        ncol(nbr_kernel) != memb$n_adm1)
      stop("generate_rprocess_csnippet: nbr_kernel must be an ", memb$n_adm1,
           "x", memb$n_adm1, " matrix (= memb$n_adm1); got ",
           paste(dim(nbr_kernel), collapse = "x"), ".")
  }

  use_mon    <- .is_mon_meas(measurement)
  is_dur_grouped <- isTRUE(use_dur) && !is.null(dur_bh_by_group)

  # The mon mu block reads av_new/mv_new, which only the duration promotion
  # block creates. Callers hard-require a duration fit on the mon path; this is
  # the generator-level restatement of that contract, so a mis-wired caller
  # fails here rather than emitting C that will not compile.
  if (use_mon && !isTRUE(use_dur))
    stop("generate_rprocess_csnippet: measurement='", measurement, "' requires ",
         "use_dur=TRUE (mso/active are duration-model latent states).")

  # --- helper: collapse integer vector to C initializer string ---
  csv_int <- function(v) paste(as.integer(v), collapse = ", ")

  # --- helper: collapse double vector to C initializer string, round-trip
  #     exact at %.17g so the embedded WNBR kernel cannot drift from the R
  #     reference (simulator_core.R). ---
  csv_dbl <- function(v) paste(sprintf("%.17g", as.double(v)), collapse = ", ")

  # Onset-month exposure correction, emitted as the SAME double the v1 kernel
  # gets from log(.ONSET_EXPOSURE_FRAC) -- printed at %.17g (round-trip exact)
  # rather than written as M_LN2, so the C and R kernels cannot drift.
  # See .ONSET_EXPOSURE_FRAC (simulator_core.R) for what this corrects.
  onset_adj_c <- sprintf("%.17g", log(.ONSET_EXPOSURE_FRAC))
  # Cessation-month twin (cessation_exposure: half), printed the same way.
  cess_adj_c  <- sprintf("%.17g", log(.CESSATION_EXPOSURE_FRAC))

  # --- helper: generate an if/else-if chain binning a double C expression
  #     into a 0-based band index, mirroring .dur_bin_mso()'s cut(breaks,
  #     right=FALSE, include.lowest=TRUE) semantics exactly. breaks is a
  #     numeric vector (last element Inf); band i (0-based) covers
  #     [breaks[i+1], breaks[i+2]). Values below breaks[1] (e.g. mso=0, not
  #     yet active) fall through to bin_var = -1 (NA), matched downstream to
  #     a 0.0 contribution -- same "no column matched" semantics as
  #     basehaz_dur_eval()'s R reference.
  dur_bin_lookup_c <- function(breaks, val_expr, bin_var) {
    k <- length(breaks) - 1L
    lines <- character(k)
    for (i in seq_len(k)) {
      lo <- breaks[i]; hi <- breaks[i + 1L]
      cond <- if (is.infinite(hi)) {
        sprintf("(%s) >= %s", val_expr, format(lo, scientific = FALSE))
      } else {
        sprintf("(%s) >= %s && (%s) < %s", val_expr, format(lo, scientific = FALSE),
                 val_expr, format(hi, scientific = FALSE))
      }
      lines[i] <- paste0(if (i == 1L) "if" else "else if", " (", cond, ") { ",
                          bin_var, " = ", i - 1L, "; }\n")
    }
    paste0("int ", bin_var, " = -1;\n", paste(lines, collapse = ""))
  }

  n_adm1  <- memb$n_adm1
  n_ct    <- memb$n_ct
  n_tu    <- memb$n_tu
  # Guard against zero-length A_tu array (all units have NA TU).
  n_tu_nz <- max(n_tu, 1L)

  fixed_d_c <- sprintf("%.1f", as.double(fixed_draw))  # "1.0"

  # Extract covariate base offsets (0-based positions in __covars).
  # cov_bases["name"] = group_index(name) * U; group_index is 0-based rank
  # of "name" among all non-unit/non-time columns of the augmented covariate
  # table (long format). spatPomp lays out __covars as U consecutive values
  # per covariate group, in column order.
  lps_base  <- as.integer(cov_bases[["lp_static"]])
  loff_base <- as.integer(cov_bases[["log_offset"]])
  n_tb      <- if (isTRUE(ttlo_none)) 8L else 7L   # occurrence TTLO bins
  tb_bases  <- as.integer(cov_bases[paste0("ttlo_coef_b", seq_len(n_tb) - 1L)])
  if (anyNA(tb_bases))
    stop("generate_rprocess_csnippet: missing ttlo_coef_b* covariate column(s) for ", n_tb, " bins.")
  sot_base  <- as.integer(cov_bases[["spill_o_t"]])
  sof_base  <- as.integer(cov_bases[["spill_o_f"]])
  soct_base <- as.integer(cov_bases[["spill_oct_t"]])
  socf_base <- as.integer(cov_bases[["spill_oct_f"]])
  sout_base <- as.integer(cov_bases[["spill_otu_t"]])
  souf_base <- as.integer(cov_bases[["spill_otu_f"]])

  # Neighbour-contagion per-unit coefficient column (occurrence track only).
  # Absent unless use_nbr, in which case build_covar_for_tu_c() joined a
  # `nbr_coef` column (beta_nbr(u) = fe[nbr] + nested-slope sum).
  nbr_base  <- if (isTRUE(use_nbr)) as.integer(cov_bases[["nbr_coef"]]) else NA_integer_


  # Duration-model covariate bases (duration-aware path only).
  if (isTRUE(onset_season) && !isTRUE(use_dur))
    stop("generate_rprocess_csnippet: onset_season requires use_dur.")
  sk_base <- if (isTRUE(onset_season)) as.integer(cov_bases[["seas_k"]]) else NA_integer_
  dur_lps_base  <- if (use_dur) as.integer(cov_bases[["dur_lp_static"]]) else NA_integer_
  dur_tb_bases  <- if (use_dur) as.integer(cov_bases[paste0("dur_ttlo_coef_b", 0:6)]) else rep(NA_integer_, 7L)
  dur_sot_base  <- if (use_dur) as.integer(cov_bases[["dur_spill_o_t"]])   else NA_integer_
  dur_sof_base  <- if (use_dur) as.integer(cov_bases[["dur_spill_o_f"]])   else NA_integer_
  dur_soct_base <- if (use_dur) as.integer(cov_bases[["dur_spill_oct_t"]]) else NA_integer_
  dur_socf_base <- if (use_dur) as.integer(cov_bases[["dur_spill_oct_f"]]) else NA_integer_
  dur_sout_base <- if (use_dur) as.integer(cov_bases[["dur_spill_otu_t"]]) else NA_integer_
  dur_souf_base <- if (use_dur) as.integer(cov_bases[["dur_spill_otu_f"]]) else NA_integer_

  # Monthly-case covariate bases (mon path only). mon_lp_static already folds
  # log(pop) and the per-(unit,time) log_offset -- see precompute_mon_lp_static.
  mon_lps_base  <- if (use_mon) as.integer(cov_bases[["mon_lp_static"]]) else NA_integer_
  mon_tb_bases  <- if (use_mon) as.integer(cov_bases[paste0("mon_ttlo_coef_b", 0:6)]) else rep(NA_integer_, 7L)
  mon_sot_base  <- if (use_mon) as.integer(cov_bases[["mon_spill_o_t"]])   else NA_integer_
  mon_sof_base  <- if (use_mon) as.integer(cov_bases[["mon_spill_o_f"]])   else NA_integer_
  mon_soct_base <- if (use_mon) as.integer(cov_bases[["mon_spill_oct_t"]]) else NA_integer_
  mon_socf_base <- if (use_mon) as.integer(cov_bases[["mon_spill_oct_f"]]) else NA_integer_
  mon_sout_base <- if (use_mon) as.integer(cov_bases[["mon_spill_otu_t"]]) else NA_integer_
  mon_souf_base <- if (use_mon) as.integer(cov_bases[["mon_spill_otu_f"]]) else NA_integer_
  mon_mso_ncol  <- if (use_mon) as.integer(mon_mso_ncol) else 0L
  mon_mso_bases <- if (use_mon && mon_mso_ncol > 0L)
    as.integer(cov_bases[paste0("mon_mso_coef_j", seq_len(mon_mso_ncol) - 1L)]) else integer(0)
  # Monthly incidence ceiling (mon_mu_cap covariate); only when a cap is set.
  use_mu_cap    <- use_mon && isTRUE(mon_mu_cap)
  mon_cap_base  <- if (use_mu_cap) as.integer(cov_bases["mon_mu_cap"]) else NA_integer_
  if (use_mu_cap && is.na(mon_cap_base))
    stop("generate_rprocess_csnippet: mon_mu_cap = TRUE but no `mon_mu_cap` covariate column.")

  # ---- Section 1: topology preamble ----
  preamble <- glue::glue(
    "/* --- topology (precomputed at build time) --- */\n",
    "#define N_UNITS  {U}\n",
    "#define N_ADM1   {n_adm1}\n",
    "#define N_CT     {n_ct}\n",
    "#define N_TU     {n_tu}\n",
    "#define N_TU_NZ  {n_tu_nz}\n",
    "static const int ADID[N_UNITS]  = {{ {csv_int(memb$adm1_id)} }};\n",
    "static const int CTID[N_UNITS]  = {{ {csv_int(memb$ct_id)} }};\n",
    "static const int TUID[N_UNITS]  = {{ {csv_int(memb$tu_id)} }};\n",
    "static const int A2C[N_ADM1]    = {{ {csv_int(memb$adm1_to_ct)} }};\n",
    "static const int C2T[N_CT]      = {{ {csv_int(memb$ct_to_tu)} }};\n",
    "/* ST(k,u): access state variable k (0-based in unit_statenames order) for unit u. */\n",
    "/* The state vector is alphabetically sorted by variable name, so stateindex[k] is */\n",
    "/* NOT the group base; stateindex[k*N_UNITS+u] gives the correct flat position.    */\n",
    "/* unit_statenames: u_days(0),o(1),o_ct(2),o_tu(3),Yout(4),draw_idx(5),p(6),mu(7),size_phi(8),active(9),mso(10) */\n",
    "#define ST(k,u) __x[__stateindex[(k)*N_UNITS+(u)]]\n",
    .trim = FALSE
  )

  # Neighbour-contagion kernel, embedded row-major in the SAME 0-based ADM1
  # order as ADID[] (both derive from unique(frame_sorted$adm1_unit)). t() makes
  # the column-major R matrix row-major for WNBR[a*N_ADM1+b]. Symmetric apart
  # from the zeroed diagonal, so transpose is cosmetic, but kept explicit so the
  # indexing contract is unambiguous. Emitted only under use_nbr.
  nbr_preamble <- if (isTRUE(use_nbr)) paste0(
    "/* --- neighbour-contagion kernel (ADM1 x ADM1, row-major, partition order) --- */\n",
    "static const double WNBR[N_ADM1*N_ADM1] = { ",
    csv_dbl(as.numeric(t(nbr_kernel))), " };\n"
  ) else ""

  # ---- Section 2: unpack state via ST(k,u) macro (O(1) code) ----
  # ST(k,u) = __x[__stateindex[k*N_UNITS+u]] defined in preamble.
  # k is 0-based in .SP_UNIT_STATENAMES order:
  #   0=u_days  1=o  2=o_ct  3=o_tu  4=Yout  5=draw_idx  6=p
  # The state vector is sorted alphabetically by spatPomp/pomp, so
  # __stateindex[k] ≠ k*U; we must use __stateindex[k*N_UNITS+u].
  unpack_state <- paste0(
    "/* --- unpack state via ST(k,u) = __x[__stateindex[k*N_UNITS+u]] --- */\n",
    "double ud[N_UNITS], ov[N_UNITS], oc[N_UNITS], ot[N_UNITS];\n",
    "{ int u;\n",
    "  for (u = 0; u < N_UNITS; u++) {\n",
    "    ud[u] = ST(0, u);\n",
    "    ov[u] = ST(1, u);\n",
    "    oc[u] = ST(2, u);\n",
    "    ot[u] = ST(3, u);\n",
    "  }\n",
    "}\n",
    if (isTRUE(use_dur)) paste0(
      "/* --- Phase 9: unpack active(9)/mso(10) -- only meaningful when use_dur --- */\n",
      "double av[N_UNITS], mv[N_UNITS];\n",
      "{ int u;\n",
      "  for (u = 0; u < N_UNITS; u++) {\n",
      "    av[u] = ST(9, u);\n",
      "    mv[u] = ST(10, u);\n",
      "  }\n",
      "}\n"
    ) else "",
    if (use_det) paste0(
      "/* --- detected feed: unpack det(11) --- */\n",
      "double dv[N_UNITS];\n",
      "{ int u; for (u = 0; u < N_UNITS; u++) dv[u] = ST(11, u); }\n"
    ) else ""
  )

  # ---- Section 3: unpack covariates via direct __covars access (O(1) code) ----
  # __covars is the full covariate vector at time t (all U units, all groups).
  # Bases are embedded at R build time: base = group_index * U.
  # This bypasses __covindex entirely (which is per-snippet and string-unindexable).
  unpack_cov <- paste0(
    "/* --- unpack covariates via direct __covars access (O(1) code) --- */\n",
    "/* Bases embedded at R build time from augmented covariate column order. */\n",
    "double lps[N_UNITS], loff[N_UNITS];\n",
    sprintf("double tb[%d][N_UNITS];\n", n_tb),
    "double sot[N_UNITS], sof[N_UNITS];\n",
    "double soct[N_UNITS], socf[N_UNITS];\n",
    "double sout[N_UNITS], souf[N_UNITS];\n",
    "{ int u;\n",
    "  for (u = 0; u < N_UNITS; u++) {\n",
    "    lps[u]   = __covars[", lps_base,  " + u];\n",
    "    loff[u]  = __covars[", loff_base, " + u];\n",
    paste0("    tb[", seq_len(n_tb) - 1L, "][u] = __covars[", tb_bases, " + u];\n", collapse = ""),
    "    sot[u]   = __covars[", sot_base,  " + u];\n",
    "    sof[u]   = __covars[", sof_base,  " + u];\n",
    "    soct[u]  = __covars[", soct_base, " + u];\n",
    "    socf[u]  = __covars[", socf_base, " + u];\n",
    "    sout[u]  = __covars[", sout_base, " + u];\n",
    "    souf[u]  = __covars[", souf_base, " + u];\n",
    "  }\n",
    "}\n",
    sprintf("const double *ttb[%d] = { %s };\n", n_tb,
            paste0("tb[", seq_len(n_tb) - 1L, "]", collapse = ",")),
    # Month-of-season index (season onset gate only).
    if (isTRUE(onset_season)) paste0(
      "/* --- month-of-season index (onset_eligibility: season) --- */\n",
      "double sk[N_UNITS];\n",
      "{ int u; for (u = 0; u < N_UNITS; u++) sk[u] = __covars[", sk_base, " + u]; }\n"
    ) else ""
  )

  # ---- Section 3b (nbr only): unpack per-unit neighbour-contagion coef ----
  # bnbr[u] = beta_nbr(u), the fixed + nested-slope sum precomputed at build
  # time and joined as the `nbr_coef` covar column. Emitted only under use_nbr.
  nbr_unpack <- if (isTRUE(use_nbr)) paste0(
    "/* --- unpack neighbour-contagion per-unit coef (occurrence track) --- */\n",
    "double bnbr[N_UNITS];\n",
    "{ int u; for (u = 0; u < N_UNITS; u++) bnbr[u] = __covars[", nbr_base, " + u]; }\n"
  ) else ""

  # ---- Section 4: phase 1 — LP + hazard + Yout ----
  # Duration-aware onset gate: the RNG call itself is skipped (not masked)
  # for already-active units, in ascending unit order -- reproducibility-
  # critical, mirrors simulate_one_tu()'s identical hard
  # branch (a masked always-draw would change the RNG stream the moment any
  # unit's `active` first flips to 1). Legacy path (use_dur=FALSE): identical
  # unconditional draw, byte-for-byte the onset-only code.
  #
  # Bernoulli draw via rbinom(1.0, p) (Rmath.h, available in Csnippets), NOT
  # `unif_rand() < p`: the two are both valid Bernoulli(p) samplers but do NOT
  # consume/interpret the underlying uniform draw the same way -- confirmed
  # empirically (matched-seed comparison of stats::rbinom() vs a compiled
  # unif_rand()<p reference: they silently disagree once p is away from the
  # 0/1 extremes) and confirmed that C's Rmath rbinom() IS byte-identical to
  # R's stats::rbinom() under a matched seed. The v1 kernel draws Yn via
  # stats::rbinom() (simulate_one_tu()), so the C side must
  # call rbinom() too, or R/C trajectories silently desync the first time a
  # non-degenerate probability is hit.
  # Season gate (onset_eligibility: season): also skip the
  # draw when an onset that reset u_days fell in the unit's current season.
  # ud[u] is the pre-update value, so floor(ud/30.4375 + 0.5) is the months
  # since that onset minus one; eligible iff it is >= sk[u]. NaN is eligible.
  # Same expression as simulator_core.R::seasonOnsetOpen().
  yn_draw_line <- if (isTRUE(use_dur) && isTRUE(onset_season)) {
    paste0("    if (av[u] < 0.5 && (ISNAN(ud[u]) || floor(ud[u] / 30.4375 + 0.5) >= sk[u])) ",
           "{ Yn[u] = rbinom(1.0, p); } else { Yn[u] = 0.0; }\n")
  } else if (isTRUE(use_dur)) {
    "    if (av[u] < 0.5) { Yn[u] = rbinom(1.0, p); } else { Yn[u] = 0.0; }\n"
  } else {
    "    Yn[u] = rbinom(1.0, p);\n"
  }
  # Neighbour-contagion reconstruction, ONCE per step before the per-unit loop.
  # Rebuild the lag-1 ADM1-active 0/1 vector Aa1 from the `o` state ov[u] (o is
  # written in phase 2 as Aa[ADID[u]], so every unit of an ADM1 already carries
  # that ADM1's OR-aggregate -- the assignment is exact, no max needed), then
  # nbr1[a] = sum_b WNBR[a*N_ADM1+b] * Aa1[b]. Mirrors nbr_from_adm1_active()
  # (K %*% Aa1) exactly.
  nbr_recon <- if (isTRUE(use_nbr)) paste0(
    "  int Aa1[N_ADM1];\n",
    "  { int a; for (a = 0; a < N_ADM1; a++) Aa1[a] = 0; }\n",
    "  { int uu; for (uu = 0; uu < N_UNITS; uu++) Aa1[ADID[uu]] = (int)(ov[uu] + 0.5); }\n",
    "  double nbr1[N_ADM1];\n",
    "  { int a, b;\n",
    "    for (a = 0; a < N_ADM1; a++) {\n",
    "      double s = 0.0;\n",
    "      for (b = 0; b < N_ADM1; b++) s += WNBR[a*N_ADM1 + b] * (double)Aa1[b];\n",
    "      nbr1[a] = s;\n",
    "    }\n",
    "  }\n"
  ) else ""

  # Neighbour-contagion LP contribution, gated ONLY by use_nbr. Deliberately NOT
  # under drop_country_tu_spillover / force_o_tu: nbr is an independent within-
  # ADM1-neighbourhood term, not a country/TU spillover.
  nbr_lp_line <- if (isTRUE(use_nbr))
    "    lp += bnbr[u] * nbr1[ADID[u]];\n" else ""

  phase1 <- paste0(
    "/* --- phase 1: LP, hazard, Yout --- */\n",
    "double Yn[N_UNITS], pn[N_UNITS];\n",
    "{\n",
    "  int u;\n",
    nbr_recon,
    "  for (u = 0; u < N_UNITS; u++) {\n",
    "    double lp = lps[u];\n",
    "    /* TTLO bin: yrs = (max(0, u_days) + 1) / 365 */\n",
    .ttlo_bin_c("ud[u]", "    ", spaced = TRUE, has_none = ttlo_none),
    "    lp += ttb[bin][u];\n",
    "    /* spillover */\n",
    "    int oi   = (int)(ov[u] + 0.5);\n",
    "    int oict = (int)(oc[u] + 0.5);\n",
    "    double otr = ot[u];\n",
    "    int na_tu_u = (otr < -0.5);\n",
    "    int oitu    = na_tu_u ? 0 : (int)(otr + 0.5);\n",
    "    lp += (oi   == 1) ? sot[u]  : sof[u];\n",
    if (isTRUE(drop_country_tu_spillover)) "" else
      "    lp += (oict == 1) ? soct[u] : socf[u];\n",
    if (isTRUE(drop_country_tu_spillover)) "" else
      "    if (!na_tu_u) lp += (oitu == 1) ? sout[u] : souf[u];\n",
    nbr_lp_line,
    "    /* hazard: p = 1 - exp(-exp(lp + log_offset)) */\n",
    "    double p = -expm1(-exp(lp + loff[u]));\n",
    "    if (p < 1e-9) p = 1e-9; else if (p > 1.0-1e-9) p = 1.0-1e-9;\n",
    "    pn[u] = p;\n",
    yn_draw_line,
    "  }\n",
    "}\n"
  )

  # ---- Section 4c (duration only): active/mso promotion, BEFORE lag update ----
  # Mirrors simulate_one_tu()'s ordering exactly: new-onset units (Yn==1, only
  # possible for at-risk units per the phase-1 gate) start mso=1 this month;
  # already-active units (av==1) advance mso by one. Mutually exclusive by
  # construction (phase 1 never draws Yn=1 for an already-active unit), so
  # the if/else-if below matches R's disjoint new_onset/still_active index sets.
  promotion_block <- if (isTRUE(use_dur)) paste0(
    "/* --- phase 1c (duration): active/mso promotion, BEFORE lag update --- */\n",
    "double av_new[N_UNITS], mv_new[N_UNITS];\n",
    "{ int u;\n",
    "  for (u = 0; u < N_UNITS; u++) {\n",
    "    av_new[u] = av[u]; mv_new[u] = mv[u];\n",
    "    if (Yn[u] > 0.5) { av_new[u] = 1.0; mv_new[u] = 1.0; }\n",
    "    else if (av[u] > 0.5) { mv_new[u] = mv[u] + 1.0; }\n",
    "  }\n",
    "}\n"
  ) else ""

  # ---- Section 4d (mon only): monthly-case LP + mu, AFTER promotion ----
  # Position is load-bearing: this block reads the LATENT mv_new/av_new, which
  # only promotion_block creates. Placed AFTER promotion so a unit onsetting
  # this month is seen at mso = 1 rather than 0, and BEFORE dur_cessation_block
  # so a unit whose episode ends this month still emits its final month's mean.
  # Matches simulate_one_tu()'s ordering exactly.
  # simulator.cessation_exposure: half (cess_half) moves it AFTER
  # dur_cessation_block instead, reading the av_act / mv_act snapshot so the
  # final month is still emitted, at half its exposure.
  #
  # ud/ov/oc/ot are deliberately the PRE-update values (same as the occurrence
  # LP and the v1 kernel), NOT the post-lag-update ST(k,u) reads that
  # dur_cessation_block uses.
  #
  # Writes the ST(7,u) = mu slot (rinit seeds it 0; nothing else writes it),
  # so no state vector / cache schema change. The generated C comment strings
  # below are left verbatim so the emitted code stays byte-identical.
  # Inactive units get mu = 0: the measurement gates its NegBin term on
  # `mu > 0.0`, NOT the recorded `active` state (ST(9,u)), because
  # dur_cessation_block runs AFTER this block and can reset active/mso to 0
  # in the SAME step for a unit whose episode ends this month -- gating on
  # `active` would silently drop that final month's case observation. `mu`
  # computed here is deterministically 0.0 iff av_new[u] < 0.5 (this month,
  # pre-cessation), so it is an exact, order-independent "active this month"
  # signal. mu is also consumed downstream (smoothing, forecast summaries)
  # and must not hold a counterfactual "cases if this unit were active" value.
  # Under cess_half the block runs after dur_cessation_block, so it reads the
  # snapshot av_act / mv_act (this month's activity before cessation) instead
  # of av_new / mv_new, which the cessation draw may already have reset; a
  # unit with av_act = 1 and av_new = 0 ceased this step. ud/ov/oc/ot are the
  # unpacked pre-update copies either way.
  cess_half <- isTRUE(cess_half) && isTRUE(use_mon) && isTRUE(use_dur)
  mn_av <- if (cess_half) "av_act" else "av_new"
  mn_mv <- if (cess_half) "mv_act" else "mv_new"
  mon_mu_block <- if (use_mon) {
    mso_decl <- if (mon_mso_ncol > 0L) paste0(
      "  double mn_mso[", mon_mso_ncol, "][N_UNITS];\n"
    ) else ""
    mso_read <- if (mon_mso_ncol > 0L) paste0(
      "      mn_mso[", seq_len(mon_mso_ncol) - 1L, "][u] = __covars[",
      mon_mso_bases, " + u];\n", collapse = ""
    ) else ""
    # Clamp into the fitted mso range. Mandatory, not defensive: bs(degree=3)
    # extrapolates as an unbounded cubic past its boundary knots.
    mso_lookup <- if (mon_mso_ncol > 0L) paste0(
      "      int jm = (int)(", mn_mv, "[u] + 0.5);\n",
      "      if (jm < 0) jm = 0; else if (jm > ", mon_mso_ncol - 1L,
      ") jm = ", mon_mso_ncol - 1L, ";\n",
      "      mn_lp += mn_mso[jm][u];\n"
    ) else ""
    # Onset-month exposure correction. precompute_mon_lp_static() folds the
    # per-(unit, time) log_offset into mn_lps at BUILD time, where the latent
    # mso is unknowable, so the correction cannot live there -- it has to be
    # applied here, against mv_new[u]. Gated on the same bundle flag that
    # decided whether log_offset was folded in at all (mon_use_logoffset,
    # threaded from the caller): correcting an offset that was never added
    # would be a spurious halving. Cessation is deliberately NOT corrected --
    # dur_cessation_block runs after this one. See .ONSET_EXPOSURE_FRAC.
    mon_onset_adj <- if (isTRUE(mon_use_logoffset)) paste0(
      "      if ((int)(", mn_mv, "[u] + 0.5) == 1) mn_lp += ", onset_adj_c, ";\n",
      # Cessation-month correction (cessation_exposure: half): the unit was
      # active this month and its cessation was drawn this step.
      if (cess_half) paste0("      if (av_new[u] < 0.5) mn_lp += ", cess_adj_c,
                            ";   /* ceased this step: half a month */\n") else ""
    ) else ""
    paste0(
      "/* --- phase 1d: monthly-size LP + mu (mon path, after promotion) --- */\n",
      "/* mon_lp_static already folds log(pop) + log_offset; NO mean_rate. */\n",
      "{\n",
      "  double mn_lps[N_UNITS], mn_tb[7][N_UNITS];\n",
      "  double mn_sot[N_UNITS], mn_sof[N_UNITS];\n",
      "  double mn_soct[N_UNITS], mn_socf[N_UNITS];\n",
      "  double mn_sout[N_UNITS], mn_souf[N_UNITS];\n",
      mso_decl,
      "  { int u;\n",
      "    for (u = 0; u < N_UNITS; u++) {\n",
      "      mn_lps[u]  = __covars[", mon_lps_base,  " + u];\n",
      paste0("      mn_tb[", 0:6, "][u] = __covars[", mon_tb_bases, " + u];\n", collapse = ""),
      "      mn_sot[u]  = __covars[", mon_sot_base,  " + u];\n",
      "      mn_sof[u]  = __covars[", mon_sof_base,  " + u];\n",
      "      mn_soct[u] = __covars[", mon_soct_base, " + u];\n",
      "      mn_socf[u] = __covars[", mon_socf_base, " + u];\n",
      "      mn_sout[u] = __covars[", mon_sout_base, " + u];\n",
      "      mn_souf[u] = __covars[", mon_souf_base, " + u];\n",
      mso_read,
      "    }\n",
      "  }\n",
      "  const double *mn_ttb[7] = { mn_tb[0],mn_tb[1],mn_tb[2],mn_tb[3],mn_tb[4],mn_tb[5],mn_tb[6] };\n",
      "  { int u;\n",
      "    for (u = 0; u < N_UNITS; u++) {\n",
      "      if (", mn_av, "[u] < 0.5) { ST(7, u) = 0.0; continue; }\n",
      "      double mn_lp = mn_lps[u];\n",
      .ttlo_bin_c("ud[u]", "      ", spaced = FALSE),
      "      mn_lp += mn_ttb[bin][u];\n",
      "      int oi   = (int)(ov[u]+0.5);\n",
      "      int oict = (int)(oc[u]+0.5);\n",
      "      double otr = ot[u]; int na_tu_u = (otr < -0.5);\n",
      "      int oitu = na_tu_u ? 0 : (int)(otr+0.5);\n",
      "      mn_lp += (oi   == 1) ? mn_sot[u]  : mn_sof[u];\n",
      if (isTRUE(drop_country_tu_spillover)) "" else
        "      mn_lp += (oict == 1) ? mn_soct[u] : mn_socf[u];\n",
      if (isTRUE(drop_country_tu_spillover)) "" else
        "      if (!na_tu_u) mn_lp += (oitu == 1) ? mn_sout[u] : mn_souf[u];\n",
      mon_onset_adj,
      mso_lookup,
      if (use_mu_cap) paste0(
        "      ST(7, u) = fmin(exp(mn_lp), __covars[", mon_cap_base,
        " + u]);   /* mu, held at incidence_cap x pop */\n")
      else "      ST(7, u) = exp(mn_lp);   /* mu: this month's NegBin mean */\n",
      "    }\n",
      "  }\n",
      "}\n"
    )
  } else ""

  # update_lags()'s source array: duration-aware path feeds the (post-
  # promotion) ongoing-episode indicator `active` instead of onset-only
  # `Yout`, fixing the duty-cycle bug -- see simulate_one_tu()'s identical
  # fix. TTLO (ST(0,u) below) stays keyed on Yn regardless of use_dur.
  active_src_var <- if (use_det) "dn" else if (isTRUE(use_dur)) "av_new" else "Yn"
  # TTLO resets on a DETECTED onset under the detected feed (fit-time ttlo is
  # built from recorded episodes); on any onset otherwise (unchanged).
  ttlo_reset_c <- if (use_det) "(Yn[u] > 0.5 && dn[u] > 0.5)" else "Yn[u]"

  # ---- Section 5: phase 2 — update_lags ----
  tu_block <- if (n_tu > 0L) {
    paste0(
      "  int At[N_TU_NZ]; memset(At, 0, sizeof(At));\n",
      "  { int c; for (c=0; c<N_CT; c++) { int t=C2T[c]; if(t>=0) At[t]+=Cm[c]; } }\n"
    )
  } else {
    paste0(
      "  int At[N_TU_NZ]; memset(At, 0, sizeof(At));\n"
    )
  }

  phase2 <- paste0(
    "/* --- phase 2: lag aggregation (ADM1 -> country -> TU) --- */\n",
    "{\n",
    "  int Aa[N_ADM1]; memset(Aa, 0, sizeof(Aa));\n",
    "  { int u; for (u=0; u<N_UNITS; u++) if ((int)", active_src_var, "[u]) Aa[ADID[u]]=1; }\n",
    "  int Cs[N_CT], Cm[N_CT];\n",
    "  memset(Cs, 0, sizeof(Cs)); memset(Cm, 0, sizeof(Cm));\n",
    "  {\n",
    "    int j; for (j=0; j<N_ADM1; j++) {\n",
    "      int c = A2C[j]; Cs[c] += Aa[j];\n",
    "      if (Aa[j] > Cm[c]) Cm[c] = Aa[j];\n",
    "    }\n",
    "  }\n",
    tu_block,
    "  /* --- writeback via ST(k,u) macro (k=0..6; mu=ST(7) in phase1b; size_phi=ST(8) static) --- */\n",
    "  { int u;\n",
    "    for (u = 0; u < N_UNITS; u++) {\n",
    "      ST(0, u) = ", ttlo_reset_c, " ? 0.0 : ud[u] + 365.25/12.0;\n",
    "      ST(1, u) = (double)Aa[ADID[u]];\n",
    if (isTRUE(drop_country_tu_spillover)) {
      "      ST(2, u) = 0.0;\n"
    } else {
      "      ST(2, u) = (double)((Cs[CTID[u]] - Aa[ADID[u]]) >= 1);\n"
    },
    "      ST(3, u) = (TUID[u] < 0) ? -1.0 :\n",
    # drop_country_tu_spillover takes PRECEDENCE over force_o_tu: o_tu forced
    # to 0.0 (skipping the At/Cm aggregation entirely), not 1.0, when both
    # flags are set. Mirrors simulator_core.R::simulate_one_tu()'s gate.
    if (isTRUE(drop_country_tu_spillover)) {
      "                    0.0;\n"
    } else if (isTRUE(force_o_tu)) {
      "                    1.0;\n"
    } else {
      "                    (double)((At[TUID[u]] - Cm[CTID[u]]) >= 1);\n"
    },
    "      ST(4, u) = Yn[u];\n",
    "      ST(5, u) = ", fixed_d_c, ";\n",
    "      ST(6, u) = pn[u];\n",
    "    }\n",
    "  }\n",
    "}\n"
  )

  # ---- Section 6 (duration only): cessation hazard, AFTER lag update ----
  # Runs after phase2's ST(0..3) writeback, so ST(0,u)/ST(1,u)/ST(2,u)/ST(3,u)
  # already hold the POST-update u_days/o/o_ct/o_tu values -- read directly
  # instead of recomputing, exactly matching the v1 kernel's use of the
  # post-update_lags() u_days/o/o_ct/o_tu.
  # RNG draw gated on av_new[u]>0.5 (currently-active units only), skipping
  # the rbinom() call entirely for inactive units -- same reproducibility
  # rationale as the phase-1 onset gate. No clamp on p_dur (deliberate
  # asymmetry vs. occurrence's p, matching the v1 kernel
  # exactly). Uses rbinom(1.0, p_dur) rather than `unif_rand() < p_dur` for
  # the same reason as the phase-1 onset draw above -- see that comment.
  dur_cessation_block <- if (isTRUE(use_dur)) {
    # DBH (mso-band basehaz lookup table) constants: declared ONCE, outside
    # the per-unit loop below (static const arrays with compile-time-constant
    # initializers -- no benefit and worse readability if placed inside the
    # loop). dur_bh_lookup_stmt (per-unit bin index + table read) is placed
    # inside the loop, where `u` and `mso_u` are in scope.
    dur_bh_const_decl <- ""
    dur_bh_lookup_stmt <- "double dbh_val = 0.0;\n"
    if (is_dur_grouped) {
      breaks <- dur_bh_by_group[[1]]$breaks
      kind1  <- dur_bh_by_group[[1]]$kind %||% "none"
      if (identical(kind1, "factor")) {
        k_mso  <- length(breaks) - 1L
        n_grp  <- length(dur_bh_by_group)
        tbl_rows <- vapply(dur_groups, function(g) {
          paste0("{ ", paste(sprintf("%.17g", dur_bh_by_group[[g]]$table), collapse = ", "), " }")
        }, character(1))
        dur_bh_const_decl <- paste0(
          "#define N_GRP  ", n_grp, "\n",
          "#define K_MSO  ", k_mso, "\n",
          "static const int GRPID[N_UNITS] = { ", csv_int(dur_grp_id), " };\n",
          "static const double DBH[N_GRP][K_MSO] = { ", paste(tbl_rows, collapse = ", "), " };\n"
        )
        dur_bh_lookup_stmt <- paste0(
          dur_bin_lookup_c(breaks, "mso_u", "bin_m"),
          "double dbh_val = (bin_m >= 0) ? DBH[GRPID[u]][bin_m] : 0.0;\n"
        )
      }
    } else if (!is.null(dur_bh) && identical(dur_bh$kind %||% "none", "factor")) {
      breaks <- dur_bh$breaks
      k_mso  <- length(breaks) - 1L
      dur_bh_const_decl <- paste0(
        "#define K_MSO  ", k_mso, "\n",
        "static const double DBH[K_MSO] = { ",
        paste(sprintf("%.17g", dur_bh$table), collapse = ", "), " };\n"
      )
      dur_bh_lookup_stmt <- paste0(
        dur_bin_lookup_c(breaks, "mso_u", "bin_m"),
        "double dbh_val = (bin_m >= 0) ? DBH[bin_m] : 0.0;\n"
      )
    }

    paste0(
      "/* --- phase 3 (duration): cessation hazard, AFTER lag update --- */\n",
      "{\n",
      dur_bh_const_decl,
      "  double dur_lps[N_UNITS];\n",
      "  double dur_tb[7][N_UNITS];\n",
      "  double dur_sot[N_UNITS], dur_sof[N_UNITS];\n",
      "  double dur_soct[N_UNITS], dur_socf[N_UNITS];\n",
      "  double dur_sout[N_UNITS], dur_souf[N_UNITS];\n",
      "  { int u;\n",
      "    for (u = 0; u < N_UNITS; u++) {\n",
      "      dur_lps[u]  = __covars[", dur_lps_base,  " + u];\n",
      paste0("      dur_tb[", 0:6, "][u] = __covars[", dur_tb_bases, " + u];\n", collapse = ""),
      "      dur_sot[u]  = __covars[", dur_sot_base,  " + u];\n",
      "      dur_sof[u]  = __covars[", dur_sof_base,  " + u];\n",
      "      dur_soct[u] = __covars[", dur_soct_base, " + u];\n",
      "      dur_socf[u] = __covars[", dur_socf_base, " + u];\n",
      "      dur_sout[u] = __covars[", dur_sout_base, " + u];\n",
      "      dur_souf[u] = __covars[", dur_souf_base, " + u];\n",
      "    }\n",
      "  }\n",
      "  const double *dur_ttb[7] = { dur_tb[0],dur_tb[1],dur_tb[2],dur_tb[3],dur_tb[4],dur_tb[5],dur_tb[6] };\n",
      "  { int u;\n",
      "    for (u = 0; u < N_UNITS; u++) {\n",
      "      if (av_new[u] < 0.5) continue;\n",
      "      double dlp = dur_lps[u];\n",
      "      double ud_n = ST(0, u);\n",
      .ttlo_bin_c("ud_n", "      ", spaced = TRUE),
      "      dlp += dur_ttb[bin][u];\n",
      "      int oi   = (int)(ST(1, u) + 0.5);\n",
      "      int oict = (int)(ST(2, u) + 0.5);\n",
      "      double otr = ST(3, u);\n",
      "      int na_tu_u = (otr < -0.5);\n",
      "      int oitu    = na_tu_u ? 0 : (int)(otr + 0.5);\n",
      "      dlp += (oi   == 1) ? dur_sot[u]  : dur_sof[u];\n",
      if (isTRUE(drop_country_tu_spillover)) "" else
        "      dlp += (oict == 1) ? dur_soct[u] : dur_socf[u];\n",
      if (isTRUE(drop_country_tu_spillover)) "" else
        "      if (!na_tu_u) dlp += (oitu == 1) ? dur_sout[u] : dur_souf[u];\n",
      "      double mso_u = mv_new[u];\n",
      dur_bh_lookup_stmt,
      "      dlp += dbh_val;\n",
      # Onset-month exposure correction. Unconditional (duration's
      # offset(log_offset) is never optional), and keyed off mso_u = mv_new[u],
      # the PROMOTED value -- so an episode onsetting this month sees ~half a
      # month of cessation exposure rather than a full one. Matches
      # the v1 kernel's `.onset_log_offset_adj(mso)`. See .ONSET_EXPOSURE_FRAC.
      "      double lambda_dur = exp(dlp + loff[u] +\n",
      "                              ((int)(mso_u + 0.5) == 1 ? ", onset_adj_c, " : 0.0));\n",
      "      double p_dur = -expm1(-lambda_dur);\n",  # deliberately no clamp
      "      if (rbinom(1.0, p_dur) > 0.5) { av_new[u] = 0.0; mv_new[u] = 0.0; }\n",
      "    }\n",
      "  }\n",
      "}\n",
      "/* --- writeback active(9)/mso(10) (duration only) --- */\n",
      "{ int u;\n",
      "  for (u = 0; u < N_UNITS; u++) {\n",
      "    ST(9, u)  = av_new[u];\n",
      "    ST(10, u) = mv_new[u];\n",
      "  }\n",
      "}\n"
    )
  } else ""

  # ---- Section 4c' (detected feed only): detection flag, right after promotion ----
  # Same rule, same unit order, same position in the step (after promotion,
  # before any other RNG use) as update_detection() in the v1 kernel, so R and
  # C consume the RNG stream identically. dn[] feeds phase 2's lag aggregation
  # and the TTLO reset; it is written back to det(11) here.
  detect_block <- if (use_det) sprintf(paste0(
    "/* --- detected feed: per-episode detection flag --- */\n",
    "double dn[N_UNITS];\n",
    "{ int u;\n",
    "  for (u = 0; u < N_UNITS; u++) {\n",
    "    if (Yn[u] > 0.5)          dn[u] = rbinom(1.0, %.17g);\n",
    "    else if (av_new[u] > 0.5) dn[u] = dv[u];\n",
    "    else                      dn[u] = rbinom(1.0, %.17g);\n",
    "    ST(11, u) = dn[u];\n",
    "  }\n",
    "}\n"), 1 - delta_fn, delta_fp) else ""

  # ---- Section 4d (WHO measurement only): simulated monthly cases ----
  # cm(u) is the C twin of simulator_core.R::drawMonthlyCases(): a NegBin draw
  # with mean mu = ST(7,u) and size size_phi = ST(8,u) for every unit active
  # this month (mu > 0 iff active, see mon_mu_block), held at floor(cap) under
  # the incidence cap, 0 otherwise. R's rnbinom(mu =) calls the same
  # rnbinom_mu(). ccm(u) is the country sum of cm, replicated on the country's
  # units for the per-unit dunit. Runs after mon_mu_block and before phase 2 /
  # cessation, so it reads this month's mu; emitted only for `_who`, so every
  # other run draws exactly the same random numbers as before.
  who_block <- if (!is.null(who_slots)) {
    if (!use_mon) stop("generate_rprocess_csnippet: who_slots needs the mon path.")
    paste0(
      "/* --- WHO term: simulated monthly cases cm(", who_slots[["cm"]],
      ") + country sum ccm(", who_slots[["ccm"]], ") --- */\n",
      "{ double cs[N_CT]; int u;\n",
      "  for (u = 0; u < N_CT; u++) cs[u] = 0.0;\n",
      "  for (u = 0; u < N_UNITS; u++) {\n",
      "    double mu_u = ST(7, u), phi_u = ST(8, u), c_u = 0.0;\n",
      "    if (mu_u > 0.0 && phi_u > 0.0) {\n",
      "      c_u = rnbinom_mu(phi_u, mu_u);\n",
      if (use_mu_cap) paste0(
        "      c_u = fmin(c_u, floor(__covars[", mon_cap_base, " + u]));\n") else "",
      "    }\n",
      "    ST(", who_slots[["cm"]], ", u) = c_u;\n",
      "    cs[CTID[u]] += c_u;\n",
      "  }\n",
      "  for (u = 0; u < N_UNITS; u++) ST(", who_slots[["ccm"]], ", u) = cs[CTID[u]];\n",
      "}\n")
  } else ""

  # ---- cessation_exposure: half -- snapshot, then cessation before the mean ----
  # The snapshot is this month's activity before the cessation draw (the
  # av_new predicate the full-order mon_mu_block gates on). Order under half:
  # onset, promotion, detection, snapshot, phase 2, cessation (unchanged: it
  # still reads the post-update lags ST(0..3)), mon mean, WHO draw. The mon
  # block draws nothing, so without `_who` the RNG stream is the full order's.
  if (cess_half) {
    snapshot_block <- paste0(
      "/* --- cessation_exposure half: this month's active/mso before cessation --- */\n",
      "double av_act[N_UNITS], mv_act[N_UNITS];\n",
      "{ int u; for (u = 0; u < N_UNITS; u++) { av_act[u] = av_new[u]; mv_act[u] = mv_new[u]; } }\n")
    full_code <- paste0(preamble,    "\n",
                         nbr_preamble,
                         unpack_state, "\n",
                         unpack_cov,   "\n",
                         nbr_unpack,
                         phase1,       "\n",
                         promotion_block, "\n",
                         detect_block,
                         snapshot_block,
                         phase2,       "\n",
                         dur_cessation_block,
                         mon_mu_block,
                         who_block)
    return(pomp::Csnippet(full_code))
  }

  full_code <- paste0(preamble,    "\n",
                       nbr_preamble,
                       unpack_state, "\n",
                       unpack_cov,   "\n",
                       nbr_unpack,
                       phase1,       "\n",
                       promotion_block, "\n",
                       detect_block,
                       mon_mu_block,
                       who_block,
                       phase2,       "\n",
                       dur_cessation_block)

  pomp::Csnippet(full_code)
}

# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

#' Build a spatPomp object using Csnippet rprocess/rinit (C variant).
#'
#' Wraps the occurrence (and optional duration / monthly-case) kernel in
#' compiled Csnippet rprocess / rinit, with Csnippet measurement functions that
#' read the latent hazard `p` stored in the state vector.
#'
#' REQUIRES fixed_draw_idx (errors if NULL). The N-independent-runs design
#' always provides it; this variant is not supported for the random-draw path.
#'
#' @param frame_tu          location frame for this TU (from build_location_frame).
#' @param panel_tu_warmup   climate + outbreak panel for this TU, covering at
#'                          least the warmup window (t_start >= window$start - tail
#'                          needed for build_warmup_state + the window itself).
#' @param draws_occ         full occurrence draws bundle.
#' @param window            list(start="YYYY-MM-DD", end="YYYY-MM-DD") — the
#'                          filter (warmup) window. Observation times = months
#'                          in [window$start, window$end).
#' @param n_draws_in_filter integer; ignored when fixed_draw_idx is provided (the
#'                          draw subset is always size 1 in that case).
#' @param fixed_draw_idx    integer (required): the posterior draw this object
#'                          commits to (N-independent-runs design: one
#'                          bpfilter run per replicate, each on one draw).
#' @param draws_dur         Optional duration (cessation-hazard) draws bundle
#'                          (`extract_duration_draws()`). Ignored when
#'                          `grouping$draws_dur_by_group` is supplied.
#' @param draws_mon         Optional monthly-case (`mon`) draws bundle
#'                          (`extract_monthly_draws()`), used when
#'                          `measurement` is one of `.SP_MON_MEAS`.
#' @param cb_meta_mon       Optional cb_meta sidecar for the mon fit (ungrouped).
#' @param cb_meta_dur       Optional cb_meta sidecar for the `dur` fit
#'                          (ungrouped); its panel cross-bases go under the
#'                          `.SP_DUR_CB_PREFIX` namespace.
#' @param mon_case_col      Optional column name in `obs_data`/the panel holding
#'                          OBSERVED MONTHLY cases (emitted as `obs_CM`).
#' @param force_o_tu logical; forwarded to generate_rprocess_csnippet() and
#'   used to gate .fold_observed_otu_into_lp(). See that function's docstrings
#'   for the partition_by=country rationale. Default FALSE (no change).
#' @param drop_country_tu_spillover logical; forwarded to
#'   generate_rprocess_csnippet() (suppresses o_ct/o_tu LP contribution and
#'   forces both to 0 post-update) and used to zero the o_ct/o_tu warmup-init
#'   state and gate .fold_observed_otu_into_lp() off. Takes PRECEDENCE over
#'   force_o_tu. Default FALSE (no change).
#' @param use_nbr logical; enable the distance-weighted neighbour-contagion term
#'   (occurrence track). Forwarded to generate_rprocess_csnippet() and used to
#'   build the per-partition kernel + beta_nbr coef column. Default FALSE (no
#'   change; byte-identical to pre-nbr).
#' @param nbr_rho_km numeric; kernel decay bandwidth in km. Default 500.
#' @param nbr_gadm_path path to the validated GADM geopackage (ADM1 centroids).
#' @param nbr_cache_dir directory for the cached ADM1 distance matrix.
#' @return spatPomp object with Csnippet rprocess and rinit.
build_spatpomp_tu_c <- function(frame_tu,
                                 panel_tu_warmup,
                                 draws_occ,
                                 window,
                                 n_draws_in_filter = 1L,
                                 obs_data          = NULL,
                                 fixed_draw_idx    = NULL,
                                 na_fill_obs       = NA_integer_,
                                 cb_meta           = NULL,
                                 grouping          = NULL,
                                 measurement       = "occurrence_bernoulli",
                                 otu_obs           = NULL,
                                 force_o_tu        = FALSE,
                                 draws_dur         = NULL,
                                 draws_mon         = NULL,
                                 cb_meta_mon       = NULL,
                                 cb_meta_dur       = NULL,
                                 mon_case_col      = NULL,
                                 delta_fn          = .SP_STATE_DELTA_FN,
                                 delta_fp          = .SP_STATE_DELTA_FP,
                                 # Inactive-state case mixture of the `_actn`
                                 # measurements (eps, mu_0). Inert
                                 # for every other measurement.
                                 case_fp_rate      = .SP_CASE_FP_RATE,
                                 case_fp_mean      = .SP_MON_MU_FLOOR,
                                 drop_country_tu_spillover = FALSE,
                                 use_nbr           = FALSE,
                                 nbr_rho_km        = 500L,
                                 nbr_gadm_path     = "data/geodata/gadm_valid.gpkg",
                                 nbr_cache_dir     = "generated_data/",
                                 # NULL or a number > 0: the mon mean is held at
                                 # incidence_cap x pop (a per-unit `mon_mu_cap`
                                 # covariate read by mon_mu_block). NULL or Inf
                                 # emits nothing, byte-identical to an uncapped
                                 # build. See simulator.incidence_cap.
                                 incidence_cap     = NULL,
                                 nbr_spec          = NULL,
                                 covariate_feed    = "latent",
                                 # NULL, or the WHO dashboard term inputs
                                 # (`_who` measurement only): list(tbl =
                                 # tibble(country, t_start, cases_who), w =
                                 # per-unit share, f = per-unit covered
                                 # population share of its country (both in
                                 # frame_tu row order), phi, mu_floor, scale).
                                 who_obs           = NULL,
                                 # "inactive" (default, unchanged) or "season":
                                 # simulator.onset_eligibility.
                                 onset_eligibility = "inactive",
                                 # "full" (default, unchanged) or "half":
                                 # simulator.cessation_exposure.
                                 cessation_exposure = "full") {

  if (is.null(fixed_draw_idx))
    stop("build_spatpomp_tu_c: fixed_draw_idx is required (C variant only ",
         "supports the N-independent-runs design.")

  stopifnot(nrow(frame_tu) > 0, !is.null(draws_occ))
  .check_measurement(measurement, "build_spatpomp_tu_c")

  U      <- nrow(frame_tu)

  # WHO dashboard term: the measurement and its inputs go together.
  use_who <- .is_who_meas(measurement)
  if (use_who && is.null(who_obs))
    stop("build_spatpomp_tu_c: measurement '", measurement, "' needs who_obs ",
         "(simulator.filter.who_obs).")
  if (!use_who && !is.null(who_obs))
    stop("build_spatpomp_tu_c: who_obs given but measurement '", measurement,
         "' does not score it.")
  if (use_who) {
    if (length(who_obs$w) != U || length(who_obs$f) != U)
      stop("build_spatpomp_tu_c: who_obs$w / who_obs$f must have one value per frame_tu row.")
    # Per-unit estimated scale (who_obs.scale_mode: estimated).
    # Read with [[ ]]: `who_obs$scale_u` is safe, but a `$s`-style prefix would
    # partial-match `scale`.
    su <- who_obs[["scale_u"]]
    if (!is.null(su) && (length(su) != U || any(!is.finite(su)) || any(su <= 0)))
      stop("build_spatpomp_tu_c: who_obs$scale_u must be one finite value > 0 per frame_tu row.")
    if (is.null(who_obs$tbl))
      who_obs$tbl <- data.frame(country = character(0), t_start = as.Date(character(0)),
                                cases_who = numeric(0))
  }

  # Filter-window months
  filt_months <- .sp_filter_window_months(panel_tu_warmup, window, frame_tu,
                                          "build_spatpomp_tu_c")
  M <- length(filt_months)

  # Single-draw subset (fixed_draw_idx path always used here).
  # Under grouping the representative draws_occ may have fewer draws than some
  # group bundle; per-group clamping in .build_grouping_sub guarantees a valid
  # index per group, so the upper-bound check is relaxed to the grouped path.
  d_idx <- as.integer(fixed_draw_idx)
  if (is.null(grouping)) stopifnot(d_idx >= 1L, d_idx <= draws_occ$S)
  else                   stopifnot(d_idx >= 1L)
  draws_sub <- .subset_draws(draws_occ, min(d_idx, draws_occ$S))
  n_sub     <- 1L

  # Country-grouping mode: per-group single-draw subsets (shared helper in
  # spatpomp_core.R). NULL when ungrouped (zero behaviour change).
  grouping_sub <- .build_grouping_sub(grouping, fixed_draw_idx,
                                       n_draws_in_filter, frame_tu)

  # Duration (cessation-hazard) draws: single-draw subset on the occurrence
  # draw index. Grouped takes priority when the caller supplies a grouped
  # duration bundle; ungrouped falls back to a single-draw subset clamped to
  # draws_dur$S. NULL
  # draws_dur/grouping$draws_dur_by_group => use_dur stays FALSE, zero
  # behaviour change vs. the onset-only path.
  draws_dur_sub    <- NULL
  grouping_dur_sub <- NULL
  if (!is.null(grouping) && !is.null(grouping$draws_dur_by_group)) {
    grouping_dur_sub <- .build_grouping_dur_sub(grouping, fixed_draw_idx,
                                                n_draws_in_filter, frame_tu)
  } else if (!is.null(draws_dur)) {
    dur_idx       <- min(d_idx, draws_dur$S)
    draws_dur_sub <- .subset_draws(draws_dur, dur_idx)
  }
  use_dur <- !is.null(draws_dur_sub) || !is.null(grouping_dur_sub)
  onset_season <- identical(checkOnsetEligibility(onset_eligibility,
                                                  "build_spatpomp_tu_c: onset_eligibility"), "season")
  if (onset_season && !use_dur)
    stop("build_spatpomp_tu_c: onset_eligibility = 'season' requires a duration fit.")
  cess_half <- identical(checkCessationExposure(cessation_exposure,
                                                "build_spatpomp_tu_c: cessation_exposure"), "half")

  # Monthly-case (`mon`) draws + the latent-mso basis table. Mirrors the
  # duration block above; additionally hard-requires a duration fit, because `active`
  # and `mso` are real latent states owned by the duration model and there is
  # no invented fallback for them.
  use_mon            <- .is_mon_meas(measurement)
  draws_mon_sub      <- NULL
  grouping_mon_sub   <- NULL
  size_phi_vec       <- NULL    # mon dispersion -> the `size_phi` state slot
  mso_basis          <- NULL
  mso_basis_by_group <- NULL
  if (use_mon) {
    if (!use_dur)
      stop("build_spatpomp_tu_c: measurement='", measurement, "' requires a ",
           "duration fit (`simulator.duration_fit`) -- `active`/`mso` are ",
           "duration-model latent states.")
    if (is.null(draws_mon) && is.null(grouping))
      stop("build_spatpomp_tu_c: measurement='", measurement, "' requires ",
           "draws_mon (ungrouped) or grouping$draws_mon_by_group (grouped).")
    if (!is.null(grouping) && !is.null(grouping$draws_mon_by_group)) {
      grouping_mon_sub   <- .build_grouping_mon_sub(grouping, fixed_draw_idx,
                                                    n_draws_in_filter, frame_tu)
      mso_basis_by_group <- lapply(grouping$cb_meta_mon_by_group %||% list(),
                                    build_mso_basis)
    } else {
      mon_idx       <- min(d_idx, draws_mon$S)
      draws_mon_sub <- .subset_draws(draws_mon, mon_idx)
      mso_basis     <- build_mso_basis(cb_meta_mon)
    }
  }

  # Warmup state. TTLO level "none": units with no recorded
  # onset before the window start are seeded at the sentinel when the
  # occurrence fit carries the level; NULL otherwise (unchanged).
  ttlo_none <- .draws_has_ttlo_none(if (!is.null(grouping)) grouping$draws_occ_by_group else draws_occ)
  warmup_state <- build_warmup_state(panel_tu_warmup, frame_tu,
                                      horizon = list(start = window$start,
                                                     end   = window$end),
                                      no_prior_onset = .ttlo_no_prior_for_frame(
                                        frame_tu, window$start, ttlo_none))
  # Detected feed: seed det(11) from the warm start's own flag when it carries
  # one, else from `active` (an observed-panel warm start is detected).
  if (identical(covariate_feed, "detected"))
    warmup_state$det <- warmup_state$det %||% warmup_state$active
  if (isTRUE(drop_country_tu_spillover)) {
    warmup_state$o_ct[!is.na(warmup_state$o_ct)] <- 0L
    warmup_state$o_tu[!is.na(warmup_state$o_tu)] <- 0L
  }

  # Observation table — include obs_CM column only for the mon path
  obs_tbl <- build_obs_for_tu(panel_tu_warmup, frame_tu, window, filt_months,
                               obs_data     = obs_data, na_fill = na_fill_obs,
                               mon_case_col = if (use_mon) mon_case_col else NULL,
                               who_tbl      = if (use_who) who_obs$tbl else NULL)
  # The mon dunit Csnippet references obs_CM; spatPomp only declares it when the
  # column is present in obs_tbl. Ensure it exists (all-NA when absent) so the
  # case term contributes nothing -- backward-compat with Bernoulli.
  if (use_mon && !("obs_CM" %in% names(obs_tbl)))
    obs_tbl$obs_CM <- NA_real_

  # ENSO column names. Under grouping
  # take the union across groups so every group's lag column is present.
  enso_names <- .sp_enso_covar_names(panel_tu_warmup, draws_occ, grouping)

  # DLNM cross-basis covariate column names (occurrence model). Cross-basis is
  # baked into lp_static below, so no Csnippet change is needed. Under grouping
  # take the union across groups. Empty when no cross-basis (zero change).
  cb_covar_names <- .sp_cb_covar_names(cb_meta, grouping)

  # cb_meta for the covariate table. Under grouping, pass the per-group map so
  # build_covar_for_tu builds each unit's cross-basis from ITS group's cb_meta
  # (groups may carry different single_lags / bases — a single shared cb_meta
  # would mis-parse all but one group). Ungrouped: pass the single cb_meta.
  cb_meta_by_group <- if (!is.null(grouping)) grouping$cb_meta_occ_by_group else NULL

  # duration.seasonal_term: harmonic2 -- the calendar-month harmonics a duration
  # fit carries (any group's bundle). character(0) for every fit without them,
  # which leaves the covariate table and dur_covar_map exactly as before.
  dur_harm_terms <- if (use_dur)
    intersect(c("sin_month", "cos_month", "sin2_month", "cos2_month"),
              .draws_fe_names(draws_dur_sub, grouping_dur_sub)) else character(0)

  # Base covariate table (includes cross-basis columns when cb_meta is non-NULL)
  covar_tbl <- build_covar_for_tu(panel_tu_warmup, frame_tu, filt_months,
                                   enso_cols        = enso_names,
                                   cb_meta          = if (is.null(grouping)) cb_meta else NULL,
                                   cb_meta_by_group = cb_meta_by_group,
                                   # calendar harmonics: mon's, and duration's when its fit has them
                                   mon_harmonics    = use_mon || length(dur_harm_terms) > 0L,
                                   cb_meta_mon      = if (use_mon && is.null(grouping_mon_sub)) cb_meta_mon else NULL,
                                   cb_meta_mon_by_group = if (use_mon) grouping$cb_meta_mon_by_group else NULL,
                                   cb_meta_dur      = if (use_dur && is.null(grouping_dur_sub)) cb_meta_dur else NULL,
                                   cb_meta_dur_by_group  = if (use_dur) grouping$cb_meta_dur_by_group else NULL,
                                   season_k         = onset_season)

  # Model-term -> covariate-column maps for the non-occurrence tracks' panel
  # covariates. mon's harmonics map to themselves; every cross-basis is
  # namespaced in the covariate table (each track picks its ENSO single-lag
  # independently and they would otherwise collide on one `enso_cb<k>` column)
  # and must reach evaluate_lp() under the bare fitted term. character(0)
  # whenever the track is off, which restores the pre-namespace behaviour.
  mon_covar_map <- if (use_mon) {
    harm <- intersect(c("sin_month", "cos_month", "sin2_month", "cos2_month"),
                       names(covar_tbl))
    c(stats::setNames(harm, harm), attr(covar_tbl, "mon_cb_map") %||% character(0))
  } else character(0)
  dur_covar_map <- attr(covar_tbl, "dur_cb_map") %||% character(0)
  # Duration harmonics map to themselves, like mon's, and fold into the static
  # duration LP (precompute_dur_lp_static()): calendar-month terms depend on
  # neither the latent state nor the lags, so the Csnippet needs no change.
  dur_harm <- intersect(dur_harm_terms, names(covar_tbl))
  if (length(dur_harm))
    dur_covar_map <- c(stats::setNames(dur_harm, dur_harm), dur_covar_map)

  # Precompute LP components (one-time O(M×U) R calls). Group-aware dispatch:
  # ungrouped → unchanged; grouped → per-group precompute with drop_region=TRUE,
  # stitched into full-frame order.
  lp_comp     <- precompute_lp_components(frame_tu, covar_tbl, enso_names,
                                          cb_covar_names, draws_sub, grouping_sub,
                                          use_nbr = use_nbr)
  lp_static   <- lp_comp$lp
  ttlo_coefs  <- lp_comp$ttlo
  spill_coefs <- lp_comp$spill
  # beta_nbr(u) per-unit coef (occurrence track; NULL unless use_nbr). Joined as
  # the `nbr_coef` covar column; multiplied in-snippet by the per-step WNBR matvec.
  nbr_coefs   <- lp_comp$nbr

  # Country-partition mode: fold OBSERVED cross-country o_tu contagion into the
  # static occurrence LP (no-op unless otu_obs supplied + frame_tu$tu_true set).
  # Restores the spillover term that update_lags() zeroes out within a
  # single-country partition; see .fold_observed_otu_into_lp.
  # Skipped entirely when force_o_tu is active: forcing ST(3,u) to 1.0
  # unconditionally already selects the TRUE branch every step, so folding the
  # observed-months delta on top would double-count the spillover contribution.
  # Also skipped when drop_country_tu_spillover is active: that flag
  # deliberately suppresses the o_tu contribution altogether, so re-adding an
  # observed-months delta here would defeat its purpose.
  if (!isTRUE(force_o_tu) && !isTRUE(drop_country_tu_spillover)) {
    lp_static <- .fold_observed_otu_into_lp(lp_static, covar_tbl, filt_months,
                                            frame_tu, spill_coefs, otu_obs)
  }

  # Duration precompute. Mirrors the occurrence precompute above:
  # group-aware dispatch via precompute_dur_lp_components(); the mso-band
  # cessation-hazard table (dur_bh / dur_bh_by_group) is a small static
  # lookup embedded directly into the generated Csnippet, not a covariate.
  dur_lp_static_vec  <- NULL
  dur_ttlo_coefs_df  <- NULL
  dur_spill_coefs_df <- NULL
  dur_bh             <- NULL
  dur_bh_by_group    <- NULL
  dur_groups         <- NULL
  if (use_dur) {
    dur_comp <- precompute_dur_lp_components(
      # `dur_covar_map` (term -> namespaced column), NOT `cb_covar_names` (the
      # OCCURRENCE fit's bare columns): duration picks its ENSO lag
      # independently. Falls back to the occurrence names only when duration
      # carries no cross-basis of its own, in which case the map is empty and
      # both are a no-op. The deployed per-group best-lag duration fits ARE
      # `_cb` fits (Horn: dur lag 8 vs occ 11), so this is the live case.
      # The duration harmonics (dur_harm, when the fit has them) ride in front
      # of either, so the fallback is decided on the cross-basis part alone.
      frame_tu, covar_tbl, enso_names,
      if (!length(dur_harm)) {
        if (length(dur_covar_map)) dur_covar_map else cb_covar_names
      } else {
        dur_cb <- dur_covar_map[setdiff(names(dur_covar_map), dur_harm)]
        c(dur_covar_map[dur_harm], if (length(dur_cb)) dur_cb else cb_covar_names)
      },
      draws_dur_sub     = if (is.null(grouping_dur_sub)) draws_dur_sub else NULL,
      grouping_dur_sub  = grouping_dur_sub
    )
    dur_lp_static_vec  <- dur_comp$lp
    dur_ttlo_coefs_df  <- dur_comp$ttlo
    dur_spill_coefs_df <- dur_comp$spill
    dur_groups         <- dur_comp$groups
    if (is.null(grouping_dur_sub)) {
      dur_bh <- dur_comp$bh
    } else {
      dur_bh_by_group <- dur_comp$bh_by_group
    }
  }

  # Monthly-case precompute (mon path only). Mirrors the duration precompute; the
  # latent-mso cross-basis additionally collapses to a per-unit x per-mso-value
  # coefficient table, looked up at Csnippet runtime from mv_new[u].
  mon_lp_static_vec  <- NULL
  mon_ttlo_coefs_df  <- NULL
  mon_spill_coefs_df <- NULL
  mon_mso_coefs_df   <- NULL
  # pop_year mode: each covariate row's month's population,
  # as a temporary `.pop_t` column read by precompute_mon_lp_static() and the
  # per-row ceiling below, dropped before the covariate table is finalised.
  # Covariate row `time` k holds filt_months[k + 1] (time M repeats month M;
  # build_covar_for_tu()). NULL for a frame without pop_y<YYYY> columns.
  pop_row <- NULL
  if (use_mon && framePopHasYears(frame_tu)) {
    pop_mat <- framePopMatrix(frame_tu, filt_months)
    pop_row <- pop_mat[cbind(match(covar_tbl$unit, as.character(frame_tu$gid)),
                             pmin(covar_tbl$time + 1L, M))]
    covar_tbl$.pop_t <- pop_row
  }
  if (use_mon) {
    mon_comp <- precompute_mon_lp_components(
      frame_tu, covar_tbl, mon_covar_map,
      draws_mon_sub      = if (is.null(grouping_mon_sub)) draws_mon_sub else NULL,
      grouping_mon_sub   = grouping_mon_sub,
      mso_basis          = mso_basis,
      mso_basis_by_group = mso_basis_by_group
    )
    mon_lp_static_vec  <- mon_comp$lp
    mon_ttlo_coefs_df  <- mon_comp$ttlo
    mon_spill_coefs_df <- mon_comp$spill
    mon_mso_coefs_df   <- mon_comp$mso
    size_phi_vec       <- mon_comp$phi
  }
  # Monthly incidence ceiling: only with the mon track and a finite cap.
  mon_mu_cap_vec <- if (use_mon && !is.null(incidence_cap) && is.finite(incidence_cap))
    monthlyMuCap(if (is.null(pop_row)) frame_tu$pop else pop_row, incidence_cap) else NULL
  if (!is.null(pop_row)) covar_tbl$.pop_t <- NULL

  # Augment covar table with precomputed columns + init state + (duration) dur
  # cols + (mon) monthly cols and the mon dispersion seed
  covar_tbl <- build_covar_for_tu_c(covar_tbl, lp_static, ttlo_coefs, spill_coefs,
                                      warmup_state, frame_tu,
                                      size_phi_vec        = if (use_mon) size_phi_vec else NULL,
                                      dur_lp_static_vec   = dur_lp_static_vec,
                                      dur_ttlo_coefs_df   = dur_ttlo_coefs_df,
                                      dur_spill_coefs_df  = dur_spill_coefs_df,
                                      mon_lp_static_vec   = mon_lp_static_vec,
                                      mon_ttlo_coefs_df   = mon_ttlo_coefs_df,
                                      mon_spill_coefs_df  = mon_spill_coefs_df,
                                      mon_mso_coefs_df    = mon_mso_coefs_df,
                                      nbr_coefs_vec       = nbr_coefs,
                                      mon_mu_cap_vec      = mon_mu_cap_vec,
                                      mon_mu_cap_by_row   = !is.null(pop_row))

  # Membership arrays for topology.
  # spatPomp sorts unit_names alphabetically from the data/covariate table.
  # ADID/CTID/TUID must be in unit_names (alphabetical) order, not frame_tu row order.
  # method = "radix" pins this to byte-order (C-locale-equivalent) collation,
  # matching the dplyr::arrange(..., .locale = "C") that built the pomp's
  # unit axis in spatpomp_core.R — a plain sort() under a system ICU locale
  # (e.g. en_US.UTF-8) disagrees with it (e.g. "BDI_1_1" vs "BDI_11_1") and
  # silently swaps two units' state trajectories under their labels.
  sorted_gids   <- sort(as.character(frame_tu$gid), method = "radix")
  frame_sorted  <- frame_tu[match(sorted_gids, as.character(frame_tu$gid)), ]
  memb <- make_membership_arrays_c(frame_sorted)

  # Neighbour-contagion kernel (occurrence track). Built from frame_sorted so the
  # returned matrix is indexed in unique(frame_sorted$adm1_unit) order -- exactly
  # the 0-based ADM1 order make_membership_arrays_c() assigned to memb$adm1_id /
  # ADID[]. Under partition_by=country this is the within-country kernel for free
  # (cross-country neighbours are simply not in the partition). NULL unless use_nbr.
  nbr_kernel_c <- if (isTRUE(use_nbr)) {
  # Spec resolved by the driver (resolve_nbr_spec()); legacy callers that pass
  # only the scalar nbr_* args get the same resolution here, so the kernel's
  # scale/partition always follow the fit's own nbr FE name.
  nbr_spec <- nbr_spec %||% resolve_nbr_spec(
    list(enable = TRUE, rho_km = nbr_rho_km, gadm_path = nbr_gadm_path,
         cache_dir = nbr_cache_dir), draws_occ, grouping)
  nbr_kernel_for_frame(frame_sorted, nbr_spec)
  } else NULL

  # Grouped-duration topology: dur_grp_id must be in frame_sorted (alphabetical
  # unit) order to match memb/ADID/CTID/TUID -- NOT frame_tu row order (unlike
  # precompute_dur_lp_components()'s internal grp_id, which is frame_tu-ordered
  # and used only for stitching lp_full against covar_tbl$unit row masks).
  dur_grp_id_sorted <- if (!is.null(dur_bh_by_group)) {
    as.integer(match(as.character(frame_sorted$group), dur_groups) - 1L)
  } else {
    NULL
  }

  # Compute covariate bases for Csnippet access.
  # covar_tbl is in long format (one row per unit-time pair). spatPomp reshapes
  # to wide: for each non-unit/non-time column k (0-based), U consecutive slots
  # in __covars hold the U unit values. Base for column k = k * U.
  cov_cols  <- setdiff(names(covar_tbl), c("unit", "time"))
  cov_bases <- setNames(
    as.integer((seq_along(cov_cols) - 1L) * as.integer(U)),
    cov_cols
  )
  init_bases <- cov_bases[intersect(c("u_days_init", "o_init", "o_ct_init", "o_tu_init",
                                       "active_init", "mso_init", "det_init"),
                                     names(cov_bases))]

  # Generate Csnippets. size_phi_base seeds state k=8 from the size_phi_init
  # column (the mon dispersion; absent off the mon path -> seeded 0).
  size_phi_base_val <- if (use_mon && !is.null(size_phi_vec))
    cov_bases[["size_phi_init"]] else NULL
  # WHO scratch slots: appended after every other slot (0-based indices).
  unit_statenames <- .sp_unit_statenames(covariate_feed, who = use_who)
  who_slots <- if (use_who)
    c(cm = match("cm", unit_statenames) - 1L, ccm = match("ccm", unit_statenames) - 1L)
  else NULL
  rinit_cs <- generate_rinit_csnippet(U, init_bases, fixed_draw = d_idx,
                                       size_phi_base = size_phi_base_val,
                                       zero_slots = if (use_who) unname(who_slots) else integer(0))
  rproc_cs <- generate_rprocess_csnippet(U, memb, cov_bases, fixed_draw = d_idx,
                                          ttlo_none   = ttlo_none,
                                          measurement = measurement,
                                          force_o_tu  = force_o_tu,
                                          drop_country_tu_spillover = drop_country_tu_spillover,
                                          use_dur          = use_dur,
                                          dur_bh           = dur_bh,
                                          dur_bh_by_group  = dur_bh_by_group,
                                          dur_grp_id       = dur_grp_id_sorted,
                                          dur_groups       = dur_groups,
                                          mon_mso_ncol     = if (!is.null(mon_mso_coefs_df))
                                            ncol(mon_mso_coefs_df) else 0L,
                                          mon_use_logoffset = .mon_logoffset_flag(
                                            draws_mon_sub, grouping_mon_sub),
                                          mon_mu_cap  = !is.null(mon_mu_cap_vec),
                                          use_nbr     = use_nbr,
                                          nbr_kernel  = nbr_kernel_c,
                                          covariate_feed = covariate_feed,
                                          delta_fn    = delta_fn,
                                          delta_fp    = delta_fp,
                                          who_slots   = who_slots,
                                          onset_season = onset_season,
                                          cess_half   = cess_half)

  # Measurement Csnippet: mon hurdle vs Bernoulli, hazard- vs state-conditioning
  # (.select_dunit()). runit_measure: always .SP_RUNIT — bpfilter never
  # calls rmeasure.
  # The WHO term's per-unit constants in spatPomp unit order (frame_sorted).
  who_c <- if (use_who) {
    ord <- match(sorted_gids, as.character(frame_tu$gid))
    # fixed mode: the original expression (Csnippet text unchanged);
    # estimated: whoMultiplier() per unit.
    who_s <- if (is.null(who_obs[["scale_u"]])) who_obs$scale / who_obs$f[ord]
             else whoMultiplier(who_obs$f[ord], who_obs[["scale_u"]][ord])
    # `_actn_who`: the country-inactive branch as a mixture.
    if (.is_actn_meas(measurement))
      .sp_case_who_norm_c(who_w = who_obs$w[ord], who_s = who_s,
                          phi = who_obs$phi, mu_floor = who_obs$mu_floor,
                          fp_rate = who_obs[["fp_rate"]] %||% .SP_CASE_FP_RATE)
    else
      .sp_case_who_c(who_w = who_obs$w[ord], who_s = who_s,
                     phi = who_obs$phi, mu_floor = who_obs$mu_floor)
  } else NULL
  dunit_use <- .select_dunit(measurement, delta_fn, delta_fp, who_c = who_c,
                             case_fp_rate = case_fp_rate,
                             case_fp_mean = case_fp_mean)

  # De-digit any covariate column ending in a digit (e.g. mon_mso_coef_j1 vs
  # j11) so spatPomp's internal paste0(name, unit_index) flatten cannot
  # collide at any unit count. Safe here: every name-based consumer of
  # covar_tbl's headers (cov_bases, init_bases, both Csnippet generators) has
  # already run; see .sp_dedigit_covar_names()'s docstring for the full proof.
  covar_tbl <- .sp_dedigit_covar_names(covar_tbl)

  # Build spatPomp object
  obj <- spatPomp::spatPomp(
    data             = obs_tbl,
    units            = "unit",
    times            = "time",
    t0               = 0L,
    unit_statenames  = unit_statenames,
    covar            = covar_tbl,
    rinit            = rinit_cs,
    rprocess         = pomp::discrete_time(rproc_cs, delta.t = 1),
    dunit_measure    = dunit_use,
    eunit_measure    = .SP_EUNIT,
    vunit_measure    = .SP_VUNIT,
    runit_measure    = .SP_RUNIT,
    paramnames       = character(0)
  )

  # Attach the object's metadata
  attr(obj, "frame_tu")          <- frame_tu
  attr(obj, "filt_months")       <- filt_months
  attr(obj, "warmup_state")      <- warmup_state
  attr(obj, "draws_sub")         <- draws_sub
  attr(obj, "draws_sub_by_group") <- if (!is.null(grouping_sub))
    grouping_sub$draws_sub_by_group else NULL
  attr(obj, "n_draws_in_filter") <- n_sub
  attr(obj, "window")            <- window
  attr(obj, "enso_names")        <- enso_names
  attr(obj, "fixed_draw_idx")    <- fixed_draw_idx
  attr(obj, "core_variant")      <- "c"
  attr(obj, "measurement")       <- measurement
  attr(obj, "size_phi")          <- size_phi_vec     # mon dispersion; NULL off the mon path
  attr(obj, "draws_dur_sub")     <- draws_dur_sub
  attr(obj, "grouping_dur_sub")  <- grouping_dur_sub
  attr(obj, "use_dur")           <- use_dur
  # Set only under the season gate, so a default build's attributes are unchanged.
  if (onset_season) attr(obj, "onset_eligibility") <- "season"
  # Set only under half, so a default build's attributes are unchanged.
  if (cess_half) attr(obj, "cessation_exposure") <- "half"
  attr(obj, "draws_mon_sub")     <- draws_mon_sub
  attr(obj, "grouping_mon_sub")  <- grouping_mon_sub
  attr(obj, "mso_basis")         <- mso_basis
  attr(obj, "mso_basis_by_group") <- mso_basis_by_group
  attr(obj, "mon_covar_map")     <- mon_covar_map
  attr(obj, "dur_covar_map")     <- dur_covar_map
  attr(obj, "who_obs")           <- if (use_who) c(who_obs[c("phi", "mu_floor", "scale")],
                                                    # per-unit scale in unit order, estimated mode only
                                                    if (!is.null(who_obs[["scale_u"]]))
                                                      list(scale_u = who_obs[["scale_u"]][match(sorted_gids, as.character(frame_tu$gid))])) else NULL

  obj
}
