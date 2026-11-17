# anchor.R — empirical baseline-rate anchors (offset(log(base_rate_*))) for
#   the brms tracks, so that exp(intercept) ~ 1 is the right prior scale.
#
# Provides:
#   anchor_spec()               one track's registry entry (.ANCHOR_SPEC)
#   anchor_enabled()            is the track's base-rate offset on?
#   build_base_rate_table()     per-cell pooled / geometric anchor table
#   resolve_base_rate()         cell -> group -> global fallback lookup
#   attach_base_rate()          add the anchor column to a panel
#   persist_base_rate_table()   write the table sidecar
#   resolve_base_rate_paths()   sidecar paths of a fit
#
# Details:
#   Contract per track: denominator(anchor) == product of the track's OTHER
#   offsets, so mu = <other offsets> * anchor * exp(lp).
#   Table schema = provenance.R's mean-rate table (cell / group / global /
#   meta, tiers T1-T5), so resolve_base_rate() delegates to resolve_mean_rate().
#   Callers source base.R, paths.R, provenance.R and occ_data.R.
#
# See: analysis/utils/README.md

# Why anchors exist. Every track's linear predictor is
# `mu = <offsets> * exp(intercept + ...)`, and the `Intercept` prior is a
# statement about `exp(intercept)`. Without an anchor, offsets that multiply by
# exposure-days and/or raw population make the prior-predictive total
# overshoot the observed total by 10^2-10^3 before any random effect
# contributes. The size track was well behaved because it carries an empirical
# per-row rate anchor, `offset(log(mean_rate))`, so its intercept only corrects
# a small residual. This module generalises that mechanism to every track, so
# that for each one
#
#     denominator(anchor) == product of that track's OTHER offsets
#
# and `exp(intercept) ~ 1` is the right prior scale everywhere: one intercept
# default (`normal(0, 1)`).
#
# Relation to provenance.R: `build_mean_rate_table()` / `resolve_mean_rate()` /
# `persist_mean_rate_table()` are the size-track originals and stay untouched.
# The table schema here is deliberately identical (`cell` / `group` / `global`
# / `meta`, same T1-T5 tier labels, same `mean_rate` column name), which lets
# `resolve_base_rate()` delegate to `resolve_mean_rate()` instead of
# duplicating the fallback ladder, and lets the size track be built through
# either function and compared cell-for-cell.

# ---------------------------------------------------------------------------
# Track registry
# ---------------------------------------------------------------------------

#' Per-track anchor specification.
#'
#' One row per brms track. The two load-bearing fields are `view_fn` and
#' `column`:
#'
#'   `view_fn`  reduces that track's prepared panel to the four columns the
#'              shared builder needs -- `country`, `admin_level`, `num`, `den`.
#'              `den` MUST be exactly the product of the track's other offsets,
#'              since that is the identity the whole mechanism rests on.
#'   `column`   the name the resolved anchor takes as a data column. Size keeps
#'              `mean_rate` (its offset term and its 88-file kernel footprint
#'              already say `mean_rate`); the new tracks get distinct names so a
#'              frame carrying several of them cannot alias.
#'
#' Helper functions are stored as STRINGS and resolved with `match.fun()` at call
#' time, so this list can be defined before the functions it names exist -- it is
#' a top-level assignment in a sourced-only module, and `analysis/utils.R`'s
#' source order is not something a registry should depend on.
.ANCHOR_SPEC <- list(
  occ = list(
    label = "occurrence",
    opt_block = "occurrence",
    cfg_key = "base_rate",
    column = "base_rate_occ",
    tier_col = "base_rate_occ_tier",
    dat_fn = "make_surv_dat_file",
    kind = "occ_base_rate",
    attr_name = "base_rate_table",
    view_fn = ".occ_anchor_view",
    estimator = "pooled",
    has_pop = FALSE,
    den_desc = "exposure days (exp(log_offset))"
  ),
  dur = list(
    label = "duration",
    opt_block = "duration",
    cfg_key = "base_rate",
    column = "base_rate_dur",
    tier_col = "base_rate_dur_tier",
    dat_fn = "make_duration_dat_file",
    kind = "dur_base_rate",
    attr_name = "base_rate_table",
    view_fn = ".dur_anchor_view",
    estimator = "pooled",
    has_pop = FALSE,
    den_desc = "exposure days (exp(log_offset))"
  ),
  mon = list(
    label = "monthly",
    opt_block = "monthly",
    cfg_key = "base_rate",
    column = "base_rate_mon",
    tier_col = "base_rate_mon_tier",
    dat_fn = "make_monthly_dat_file",
    kind = "mon_base_rate",
    attr_name = "base_rate_table",
    view_fn = ".mon_anchor_view",
    estimator = "pooled",
    has_pop = TRUE,
    den_desc = "pop * exposure days (exp(log(pop) + log_offset))"
  ),
  siz = list(
    label = "size",
    opt_block = "size",
    cfg_key = "mean_rate",
    column = "mean_rate",
    tier_col = "mr_tier",
    dat_fn = ".siz_anchor_dat_file",
    kind = "size_mean_rate",
    # Size keeps the ORIGINAL attribute name for the same reason it keeps the
    # `mean_rate` column name: provenance.R's persist_mean_rate_table() and
    # tests/gates_estimand_coherence.R read `attr(d, "mean_rate_table")`, and
    # renaming it would break them for no gain. The attribute name is per-track
    # metadata, not a constant.
    attr_name = "mean_rate_table",
    view_fn = ".siz_anchor_view",
    # POOLED: all four tracks share one estimator. Against `geometric` it moves the
    # frame aggregate by +47.8% (median cell ratio 1.245, 24.4% of cells >2x) but
    # changes the tier of ZERO cells and zero panel rows -- values move, the
    # information source does not. The aggregate is absorbed by a -0.391
    # intercept shift, well inside the anchored normal(0, 1).
    # `geometric` remains fully supported and is proven identical to the size
    # original build_mean_rate_table() cell-for-cell (relative difference exactly 0),
    # so a like-for-like refit of the archive needs only
    # `size.mean_rate.estimator: geometric` in the config.
    estimator = "pooled",
    has_pop = TRUE,
    den_desc = "pop (exp(log(pop)))"
  )
)

#' Look up (and validate) one track's anchor spec.
#'
#' @param track one of "occ", "dur", "mon", "siz".
#' @return the `.ANCHOR_SPEC` entry, with `track` added.
anchor_spec <- function(track) {
  if (
    !is.character(track) ||
      length(track) != 1L ||
      !track %in% names(.ANCHOR_SPEC)
  ) {
    stop(
      "anchor_spec: unknown track '",
      paste(track, collapse = ","),
      "' (expected one of ",
      paste(names(.ANCHOR_SPEC), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  c(.ANCHOR_SPEC[[track]], list(track = track))
}

#' Per-track anchor config block, e.g. `opt$duration$base_rate`.
.anchor_cfg <- function(opt, spec) {
  (opt[[spec$opt_block]] %||% list())[[spec$cfg_key]] %||% list()
}

#' Is the baseline-rate offset switched on for this track?
#'
#' Reads `opt$<block>$offset_baserate`. Defaults to FALSE so every config that
#' exists today stays bit-inert -- no anchor column, no formula change, no path
#' change. Size is the exception: its anchor is not gated, it has always been on.
anchor_enabled <- function(opt, track) {
  spec <- anchor_spec(track)
  if (identical(track, "siz")) {
    return(TRUE)
  }
  isTRUE((opt[[spec$opt_block]] %||% list())$offset_baserate)
}


# ---------------------------------------------------------------------------
# Per-track views: (country, admin_level, num, den)
# ---------------------------------------------------------------------------
#
# Each view applies that track's own quality guards but deliberately NOT the
# country or admin_level subsetting its fit-time loader applies. Same deviation,
# and same reason, as `build_mean_rate_table()` (provenance.R:769-774): the anchor
# table is GLOBAL -- built once, shared by every country/group/backend config --
# so a per-run subset in here would make one run's table disagree with another's
# while both resolve to the same filename.

#' Occurrence view: events per exposure-day.
#'
#' `num`/`den` come from `.occ_poisson_view()` (occ_data.R), the same helper the
#' fit scripts use to build `log_offset` -- so `den == exp(log_offset)` holds by
#' construction rather than by a comment.
.occ_anchor_view <- function(opt, d) {
  d <- .occ_poisson_view(d)
  # Same pooling-key rule as the mon/siz views: when occurrence.
  # canonical_admin_level is on, the fit's rows are canonically keyed, so the
  # anchor must be too. Default OFF -- unchanged for every existing config.
  # Inline key test (not canonical_admin_level_enabled()) so the default path
  # keeps needing nothing from utils/population.R.
  if (isTRUE((opt$occurrence %||% list())$canonical_admin_level)) {
    d <- .apply_p3(d, opt, "occurrence", tag = "build_base_rate_table/occ/P3")
  }
  d %>%
    dplyr::filter(!is.na(t_start), !is.na(t_end)) %>%
    dplyr::transmute(
      country = as.character(country),
      admin_level = as.character(admin_level),
      num = as.numeric(y),
      den = as.numeric(exposure)
    )
}

#' Duration view: cessations per exposure-day.
#'
#' `den` reproduces `dur_data.R:350`'s `log_offset = log(pmax(exposure, 1))`.
.dur_anchor_view <- function(opt, d) {
  if (isTRUE((opt$duration %||% list())$canonical_admin_level)) {
    d <- .apply_p3(d, opt, "duration", tag = "build_base_rate_table/dur/P3")
  }
  d %>%
    dplyr::transmute(
      country = as.character(country),
      admin_level = as.character(admin_level),
      num = as.numeric(ended),
      den = as.numeric(pmax(exposure, 1))
    )
}

#' Monthly view: cases per person per exposure-day, on the DECLARED population.
#'
#' `num` is ROUNDED, matching `load_monthly_data()` (mon_data.R) -- the
#' negbinomial response is the rounded count, so the anchor must be a rate for
#' the same quantity. `den` is the product of both monthly offsets.
#'
#' Estimand-aware for the same reason `.siz_anchor_view()` is: `base_rate_mon`
#' enters the fit as `offset(log(base_rate_mon))` alongside `offset(log(pop))`,
#' so a rate built on the footprint H paired with an offset built on the whole
#' unit Q is the Q x rbar^H dimensional error, and both offsets being fixed
#' means the likelihood cannot detect it. The P3 re-key is here for the same
#' reason it is in `load_monthly_data()`: `admin_level` is this table's POOLING
#' key, so a mis-keyed row is averaged into the wrong cell and handed back to
#' every consumer of that cell.
#'
#' `pop` is emitted ALONGSIDE `den` (not folded into it) purely so
#' `build_base_rate_table()` can digest the population itself for the provenance
#' stamp. Stamping `den` would digest `pop * exposure_days`, which changes when
#' the episode windows change and is therefore not a population fingerprint at
#' all -- the Step-7 ordering guard would then compare two things that are not
#' the same quantity.
#'
#' Unlike size there are no `min_pop` / `max_attack_rate` guards: those are size
#' keys, and inventing monthly equivalents here would silently change every
#' existing mon anchor cell. Left as a deliberate open question.
.mon_anchor_view <- function(opt, d) {
  mon <- opt$monthly %||% list()
  estimand <- resolve_monthly_estimand(opt)
  d <- resolve_size_pop_column(d, estimand, tag = "build_base_rate_table/mon")

  d <- .apply_p3(
    d,
    opt,
    "monthly",
    estimand,
    tag = "build_base_rate_table/mon/P3"
  )

  d %>%
    dplyr::filter(is.finite(pop), pop > 0) %>%
    dplyr::transmute(
      country = as.character(country),
      admin_level = as.character(admin_level),
      num = as.numeric(round(cases_month)),
      pop = as.numeric(pop),
      den = as.numeric(pop) * as.numeric(pmax(exposure, 1))
    )
}

#' Size panel path (the one literal `build_mean_rate_table()` defaults to).
.siz_anchor_dat_file <- function(opt) {
  opt$size$analysis_data_path %||% make_size_dat_file(opt)
}

#' Size view: cases per person, on the DECLARED population denominator.
#'
#' Reproduces `build_mean_rate_table()`'s P0/P2/P3 preamble (provenance.R:741-784)
#' verbatim in effect: estimand-resolved `pop` column, canonical admin level, and
#' the `min_pop` / `max_attack_rate` quality guards. Those are genuinely
#' size-specific (they police a population denominator the other three tracks
#' either do not have or do not divide by alone), which is why they live in the
#' view rather than in the shared builder.
.siz_anchor_view <- function(opt, d) {
  min_pop <- opt$size$min_pop %||% 10
  max_ar <- opt$size$max_attack_rate %||% 0.10

  estimand <- resolve_size_estimand(opt)
  d <- resolve_size_pop_column(d, estimand, tag = "build_base_rate_table/siz")

  d <- .apply_p3(d, opt, "size", estimand, tag = "build_base_rate_table/siz/P3")

  d %>%
    dplyr::filter(!is.na(TL), !is.na(TR)) %>%
    dplyr::filter(is.na(pop) | pop >= min_pop) %>%
    dplyr::filter(is.na(pop) | !(pop > 0) | suspected_cases <= max_ar * pop) %>%
    # `pop` is emitted redundantly (it IS `den` here) so that every view whose
    # denominator contains a population carries that population under one name,
    # which is what build_base_rate_table()'s provenance digest keys on. Provably a
    # no-op for this view: identical vector, and `.anchor_engine()` reads only
    # country/admin_level/num/den.
    dplyr::transmute(
      country = as.character(country),
      admin_level = as.character(admin_level),
      num = as.numeric(suspected_cases),
      pop = as.numeric(pop),
      den = as.numeric(pop)
    )
}


# ---------------------------------------------------------------------------
# Estimator engines
# ---------------------------------------------------------------------------

# The estimator. The size original estimates each cell by a GEOMETRIC mean of
# per-row `log(cases/pop)`. That is impossible for occ / dur / mon: occ and dur
# responses are 0/1 and mon has structural zero months, so per-row `log(rate)`
# is `-Inf` on most rows. They need a POOLED estimator,
# `log(sum(num) / sum(den))` per cell, with the EB shrinkage weight rederived
# for it (below). All four tracks default to `pooled` (see the `siz` spec), so
# there is one estimator to reason about. `geometric` remains supported
# per-config and reproduces the size original cell-for-cell, which a
# like-for-like refit of an archived size fit needs.
#' Build the estimator engine for a view.
#'
#' Both estimators reduce to the same three primitives, which is what lets the
#' table assembly below be written once:
#'
#'   `units`    the rows the parents aggregate over -- observations for
#'              `geometric`, (country, admin_level) cells for `pooled`.
#'   `agg_adj`  a subset of `units` -> its ADM1-EQUIVALENT log rate (used for the
#'              country / ENSO-group / global parents, which pool across levels
#'              and so must be level-adjusted first).
#'   `agg_raw`  a subset of `units` -> its log rate AS OBSERVED (used for own-cell
#'              and own-(group, level) estimates, which are already single-level).
#'
#' GEOMETRIC (size original, provenance.R:822-953)
#'   Cell estimate is `mean(log(num/den))`; the EB weight is `n/(n + k)` with
#'   `k = sigma2_within / tau2_between` -- a variance-ratio shrinkage in which
#'   `n` counts OBSERVATIONS.
#'
#' POOLED (occ / dur / mon)
#'   Cell estimate is `log(Y_c / E_c)` with `Y_c = sum(num)`, `E_c = sum(den)`.
#'   There is no `sigma2_within` analogue -- a cell's sampling variance is not
#'   estimated from a spread of per-row log rates (there are none: most rows are
#'   0), it is implied by the Poisson likelihood. The delta-method variance of
#'   `log(Y/E)` under `Y ~ Poisson(E r)` is `1/Y`, so the EB weight is
#'
#'       w = tau2 / (tau2 + 1/Y)  =  Y / (Y + 1/tau2)
#'
#'   which is EXACTLY the existing `n/(n + k)` form under two substitutions:
#'   `n := Y_c` (event count, not row count) and `k := 1/tau2_between`. The
#'   shrinkage code below is therefore shared verbatim; only these two inputs
#'   change. `tau2_between` is a method-of-moments estimate over cells with at
#'   least one event: `var(log r_c) - mean(1/Y_c)`, i.e. observed spread minus
#'   expected sampling noise, floored at 1e-6.
#'
#'   Zero-event cells (`Y_c = 0`) get `w = 0` and fall back to the parent
#'   EXACTLY, tier "T1-zero" -- a distinct label from T2 ("cell absent entirely")
#'   so diagnostics can tell "observed, saw nothing" from "never observed". No
#'   continuity correction: a `+0.5` would bias each cell by an amount varying
#'   with its own `E_c`, whereas the parents already pool the whole corpus.
#'
#' @param d view tibble (country, admin_level, num, den), post-guards.
#' @param estimator "geometric" | "arithmetic" | "pooled".
#' @return list(units, agg_adj, agg_raw, n_of, raw_of, k, tau2, label).
.anchor_engine <- function(d, estimator) {
  if (identical(estimator, "pooled")) {
    units <- d %>%
      dplyr::group_by(country, admin_level) %>%
      dplyr::summarise(
        a = sum(num),
        b = sum(den),
        n_rows = dplyr::n(),
        .groups = "drop"
      )

    seen <- units %>% dplyr::filter(a > 0)
    if (nrow(seen) >= 2L) {
      log_r <- log(seen$a / seen$b)
      tau2 <- max(1e-6, stats::var(log_r) - mean(1 / seen$a))
    } else {
      warning(
        "build_base_rate_table: fewer than 2 cells with >=1 event -- ",
        "tau2_between not identifiable, falling back to tau2 = 1 ",
        "(equal-weight shrinkage)."
      )
      tau2 <- 1
    }
    k <- 1 / tau2

    return(list(
      units = units,
      agg_adj = function(u, lf) {
        log(sum(u$a) / sum(u$b * exp(lf[u$admin_level])))
      },
      agg_raw = function(u) {
        if (sum(u$a) > 0) log(sum(u$a) / sum(u$b)) else NA_real_
      },
      n_of = function(u) sum(u$a), # EB weight counts EVENTS
      raw_of = function(u) if (sum(u$a) > 0) sum(u$a) / sum(u$b) else 0,
      k = k,
      tau2 = tau2,
      label = "pooled"
    ))
  }

  # --- geometric / arithmetic: observation-level, log scale ------------------
  units <- d %>%
    dplyr::mutate(rate = num / den, a = log(num / den), b = 1)

  cell_var <- units %>%
    dplyr::group_by(country, admin_level) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_log = mean(a),
      var_log = if (dplyr::n() > 1) stats::var(a) else NA_real_,
      .groups = "drop"
    )
  within_ok <- !is.na(cell_var$var_log)
  sigma2_within <- if (any(within_ok)) {
    stats::weighted.mean(
      cell_var$var_log[within_ok],
      w = cell_var$n[within_ok] - 1
    )
  } else {
    NA_real_
  }
  tau2 <- stats::var(cell_var$mean_log)
  if (!is.finite(sigma2_within) || sigma2_within <= 0) {
    sigma2_within <- tau2
  }
  k <- sigma2_within / tau2
  if (!is.finite(k) || k <= 0) {
    k <- 1
  } # degenerate corpus: equal-weight shrink

  list(
    units = units,
    agg_adj = function(u, lf) mean(u$a - lf[u$admin_level]),
    agg_raw = function(u) mean(u$a),
    n_of = function(u) nrow(u),
    raw_of = function(u) {
      if (identical(estimator, "arithmetic")) {
        mean(u$rate)
      } else {
        exp(mean(u$a))
      }
    },
    k = k,
    tau2 = tau2,
    label = estimator
  )
}

#' Cross-admin-level contrast `L(level)`, CI-gated.
#'
#' Same device, same gate and same default as `build_mean_rate_table()`
#' (provenance.R:837-897): a fallback/borrowing correction that says "finer admin
#' units have higher events per unit exposure, so scale by this much when
#' borrowing a rate across levels". It is estimated country-demeaned so that
#' between-country variation cannot masquerade as a level effect, and any
#' coefficient whose CI crosses 0 is gated to 0.
#'
#' The pooled branch cannot demean per-row log rates (there are none), so it
#' estimates the contrast from a Poisson GLM on the cell table with a country
#' fixed effect -- the country FE does exactly the job the demeaning does -- and
#' `offset(log E_c)` making the coefficients log-rate ratios.
#'
#' The `lm` branch keeps `confint()` (t-quantile intervals) so the geometric
#' estimator reproduces `build_mean_rate_table()` exactly. The pooled `glm` branch
#' uses `confint.default()` explicitly: plain `confint()` on a `glm` dispatches
#' to MASS's profile-likelihood method only if MASS happens to be attached, which
#' would make a coefficient gate depend on package load order.
#'
#' @return list(level_factor named over `allowed_levels`, level_factor_ci).
.anchor_level_factor <- function(
  eng,
  allowed_levels,
  valid_levels,
  use_level_factor
) {
  level_factor <- stats::setNames(
    rep(0, length(allowed_levels)),
    allowed_levels
  )
  if (!isTRUE(use_level_factor)) {
    message(
      "[build_base_rate_table] level_factor DISABLED: all levels 0; parents ",
      "are unadjusted and borrowed rates are not re-inflated across levels."
    )
    return(list(level_factor = level_factor, level_factor_ci = NULL))
  }

  multi <- eng$units %>%
    dplyr::group_by(country) %>%
    dplyr::filter(dplyr::n_distinct(admin_level) > 1) %>%
    dplyr::ungroup()
  if (!(nrow(multi) > 0 && dplyr::n_distinct(multi$country) >= 2)) {
    warning(
      "build_base_rate_table: fewer than 2 countries span multiple admin ",
      "levels -- level_factor left at 0 for all levels."
    )
    return(list(level_factor = level_factor, level_factor_ci = NULL))
  }

  multi <- multi %>%
    dplyr::mutate(admin_level = factor(admin_level, levels = valid_levels))

  ci <- if (identical(eng$label, "pooled")) {
    fit <- stats::glm(
      a ~ admin_level + country + offset(log(b)),
      family = stats::poisson(),
      data = multi
    )
    stats::confint.default(fit)
  } else {
    dm <- multi %>%
      dplyr::group_by(country) %>%
      dplyr::mutate(a = a - mean(a)) %>%
      dplyr::ungroup()
    fit <- stats::lm(a ~ admin_level, data = dm)
    stats::confint(fit)
  }
  for (lv in setdiff(allowed_levels, "ADM1")) {
    coef_name <- paste0("admin_level", lv)
    if (coef_name %in% names(stats::coef(fit))) {
      ci_lv <- ci[coef_name, ]
      if (ci_lv[1] <= 0 && ci_lv[2] >= 0) {
        warning(
          "build_base_rate_table: level factor for ",
          lv,
          " CI crosses 0 (",
          round(ci_lv[1], 3),
          ", ",
          round(ci_lv[2], 3),
          ") -- gated to 0."
        )
      } else {
        level_factor[[lv]] <- stats::coef(fit)[[coef_name]]
      }
    }
  }
  list(level_factor = level_factor, level_factor_ci = ci)
}


# ---------------------------------------------------------------------------
# build_base_rate_table
# ---------------------------------------------------------------------------

#' Build the shared empirical baseline-rate table for one track.
#'
#' The generalisation of `build_mean_rate_table()` (provenance.R:718) to all four
#' tracks. Same output schema, same five-tier fallback ladder, same
#' `cross_absent_levels` switch, same reconciliation contract -- the differences
#' are exactly two: the per-track view (which numerator over which denominator)
#' and the estimator (see `.anchor_engine()`).
#'
#' Tiers, cheapest first:
#'   T1       own (country, admin_level) cell, EB-shrunk toward its country parent
#'   T1-zero  cell observed but with zero events -- parent exactly (pooled only)
#'   T2       country seen, this level not -- country parent, level-adjusted
#'   T3-own   country unseen -- its ENSO group's own (group, level) estimate
#'   T3-parent  ... or the group parent, level-adjusted
#'   T4       group unseen -- global parent, level-adjusted
#'   T5       nothing resolves -- hard error in `resolve_base_rate()`, never a
#'            silent default
#'
#' Config (under `opt$<block>$base_rate`, or `opt$size$mean_rate` for size):
#'   `estimator`            "pooled" | "geometric" | "arithmetic"
#'   `level_factor`         TRUE (default) -- see `.anchor_level_factor()`
#'   `cross_absent_levels`  TRUE (default) -- cross levels no row occupies, so a
#'                          query at an unobserved level degrades to T2/T3/T4
#'                          rather than erroring at T5
#'
#' @param opt options list. NOTE `opt$country` is deliberately IGNORED: this
#'   table is global, see the views' comment above.
#' @param track "occ" | "dur" | "mon" | "siz".
#' @param path override for the prepared panel (default: the track's own
#'   `make_*_dat_file(opt)`).
#' @param enso_grouping_path iso3 -> ENSO-group registry, the T3 tier's key.
#' @return list(cell, group, global, meta) -- schema identical to
#'   `build_mean_rate_table()`, so `resolve_mean_rate()` accepts it unchanged.
build_base_rate_table <- function(
  opt,
  track,
  path = NULL,
  enso_grouping_path = "analysis/groupings/enso_grouping.csv"
) {
  spec <- anchor_spec(track)
  cfg <- .anchor_cfg(opt, spec)

  if (is.null(path)) {
    path <- match.fun(spec$dat_fn)(opt)
  }
  path <- as.character(path)
  if (!file.exists(path)) {
    stop(
      "build_base_rate_table(",
      track,
      "): missing input ",
      path,
      call. = FALSE
    )
  }
  d_in <- readRDS(path)

  # Capture the population stamp BEFORE the guards, for the same reason the size
  # original does (provenance.R:724-731): `summarise()` builds a fresh frame and
  # drops non-standard attributes, so reading it here keeps this function
  # independent of which side of that line any dplyr verb falls on today.
  pop_prov_in <- if (isTRUE(spec$has_pop)) get_pop_provenance(d_in) else NULL

  # `siz` goes through the shared resolver (size_mean_rate_estimator(), paths.R) so
  # this function and make_size_run_file()/build_mean_rate_table() cannot disagree
  # about which estimator an under-specified config asked for -- see that
  # helper's doc comment for the silent wrong-table failure it closes. occ/dur/
  # mon have no such split: their path builders (paths.R:193/615/695) already
  # hardcode the same "pooled" that `spec$estimator` carries.
  estimator <- if (identical(track, "siz")) {
    size_mean_rate_estimator(opt)
  } else {
    cfg$estimator %||% spec$estimator
  }
  if (!estimator %in% c("geometric", "arithmetic", "pooled")) {
    stop(
      "build_base_rate_table(",
      track,
      "): unknown estimator '",
      estimator,
      "' (must be geometric | arithmetic | pooled)",
      call. = FALSE
    )
  }

  d <- match.fun(spec$view_fn)(opt, d_in)

  # --- shared guards --------------------------------------------------------
  # The denominator must be strictly positive under every estimator -- it is the
  # offset product, so a zero would make `log(den)` itself undefined. The
  # NUMERATOR guard is estimator-dependent, and that difference is the whole
  # reason the pooled branch exists: `geometric` takes `log(num/den)` per row and
  # so cannot see a zero-event row at all, while `pooled` MUST keep those rows --
  # they carry exposure, and dropping them would inflate every rate by exactly
  # the fraction of the corpus that saw nothing.
  n_before <- nrow(d)
  d <- dplyr::filter(d, is.finite(den), den > 0, is.finite(num), num >= 0)
  if (!identical(estimator, "pooled")) {
    d <- dplyr::filter(d, num > 0)
  }
  n_dropped <- n_before - nrow(d)
  if (n_dropped > 0) {
    warning(
      "build_base_rate_table(",
      track,
      "): dropped ",
      n_dropped,
      " row(s) with ",
      "non-positive / non-finite num or den."
    )
  }
  if (!nrow(d)) {
    stop(
      "build_base_rate_table(",
      track,
      "): no rows survive the guards",
      call. = FALSE
    )
  }

  allowed_levels <- c("ADM1", "ADM2", "ADM3")
  if (!all(d$admin_level %in% allowed_levels)) {
    bad <- setdiff(unique(d$admin_level), allowed_levels)
    stop(
      "build_base_rate_table(",
      track,
      "): unexpected admin_level value(s): ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  cross_absent <- cfg$cross_absent_levels %||% TRUE
  valid_levels <- if (isTRUE(cross_absent)) {
    allowed_levels
  } else {
    intersect(allowed_levels, unique(d$admin_level))
  }
  if (!identical(valid_levels, allowed_levels)) {
    message(sprintf(
      paste0(
        "[build_base_rate_table/%s] levels present: %s (of %s) -- absent level(s) ",
        "are NOT crossed into the cell/group/global tables ",
        "(%s.%s.cross_absent_levels = FALSE)."
      ),
      track,
      paste(valid_levels, collapse = "/"),
      paste(allowed_levels, collapse = "/"),
      spec$opt_block,
      spec$cfg_key
    ))
  }

  eng <- .anchor_engine(d, estimator)
  k <- eng$k
  lfr <- .anchor_level_factor(
    eng,
    allowed_levels,
    valid_levels,
    cfg$level_factor %||% TRUE
  )
  level_factor <- lfr$level_factor

  units <- eng$units

  # --- own cells (T1 / T1-zero) ---------------------------------------------
  own_cell <- units %>%
    dplyr::group_by(country, admin_level) %>%
    dplyr::group_modify(
      ~ tibble::tibble(
        n_cell = eng$n_of(.x),
        log_r_raw = eng$agg_raw(.x),
        mr_raw = eng$raw_of(.x)
      )
    ) %>%
    dplyr::ungroup()

  # --- country parents (T2) -------------------------------------------------
  country_tbl <- units %>%
    dplyr::group_by(country) %>%
    dplyr::group_modify(
      ~ tibble::tibble(
        log_r_country = eng$agg_adj(.x, level_factor),
        n_country = eng$n_of(.x)
      )
    ) %>%
    dplyr::ungroup()

  # --- ENSO group registry + parents (T3) -----------------------------------
  if (!file.exists(enso_grouping_path)) {
    stop(
      "build_base_rate_table(",
      track,
      "): ENSO grouping registry not found: ",
      enso_grouping_path,
      call. = FALSE
    )
  }
  grouping <- utils::read.csv(enso_grouping_path, stringsAsFactors = FALSE) %>%
    dplyr::transmute(
      country = as.character(iso3),
      enso_group = as.character(group)
    )
  u_grp <- units %>% dplyr::left_join(grouping, by = "country")
  if (any(is.na(u_grp$enso_group))) {
    missing_ctry <- unique(u_grp$country[is.na(u_grp$enso_group)])
    warning(
      "build_base_rate_table(",
      track,
      "): country/countries not in ENSO ",
      "grouping registry (T3 unavailable for them if they later drop out ",
      "of the cell table): ",
      paste(missing_ctry, collapse = ", ")
    )
  }
  u_grp <- dplyr::filter(u_grp, !is.na(enso_group))

  group_tbl <- u_grp %>%
    dplyr::group_by(enso_group) %>%
    dplyr::group_modify(
      ~ tibble::tibble(
        log_r_group = eng$agg_adj(.x, level_factor),
        n_group = eng$n_of(.x)
      )
    ) %>%
    dplyr::ungroup()
  group_level_tbl <- u_grp %>%
    dplyr::group_by(enso_group, admin_level) %>%
    dplyr::group_modify(
      ~ tibble::tibble(
        log_r = eng$agg_raw(.x),
        n = eng$n_of(.x)
      )
    ) %>%
    dplyr::ungroup()

  # --- global parent (T4) ---------------------------------------------------
  log_r_global <- eng$agg_adj(units, level_factor)

  # --- zero-event parent fall-through ---------------------------------------
  # Pooled-only by construction: `geometric` filters zero-event rows out before
  # the estimator ever sees them (the `num > 0` guard above), so a `geometric`
  # parent is always a mean over finite log-rates and provenance.R's builder has
  # no analogue of this block. Under `pooled` a country -- or an ENSO group, or a
  # (group, level) cell -- can carry exposure and zero events, giving
  # log(0 / E) = -Inf. That is NOT a crash: it flows into `mean_rate = 0` and
  # then into `offset(log(0)) = -Inf`, i.e. a silent mu = 0 for every row of that
  # country. So walk the ladder the tiers already describe -- country -> its ENSO
  # group -> global -- taking the first FINITE value at each step.
  #
  # Own cells need no repair: `n_cell == 0` iff `log_r_raw == -Inf` under pooled
  # (both are functions of `sum(num)`), and the `n_cell <= 0` branch of the
  # assembly below already routes those to `log_r_parent` as tier "T1-zero".
  if (!is.finite(log_r_global)) {
    stop(
      "build_base_rate_table(",
      track,
      "): the whole corpus has zero events ",
      "(global parent log-rate = ",
      log_r_global,
      ") -- no baseline rate is ",
      "estimable at any tier.",
      call. = FALSE
    )
  }

  group_tbl <- dplyr::mutate(
    group_tbl,
    log_r_group = ifelse(is.finite(log_r_group), log_r_group, log_r_global)
  )
  # A zero-event (group, level) cell has no own estimate, so it must fall to the
  # group parent + level factor -- which is exactly what an NA `log_r` means to
  # the `group` assembly below (tier "T3-parent"), as an absent cell does.
  group_level_tbl <- dplyr::mutate(
    group_level_tbl,
    log_r = ifelse(is.finite(log_r), log_r, NA_real_)
  )
  # Re-select the original three columns so the `cell` join below is unchanged.
  # A country absent from the ENSO registry has `enso_group = NA`, hence
  # `log_r_group = NA`, and correctly falls straight through to the global parent.
  country_tbl <- country_tbl %>%
    dplyr::left_join(grouping, by = "country") %>%
    dplyr::left_join(
      dplyr::select(group_tbl, enso_group, log_r_group),
      by = "enso_group"
    ) %>%
    dplyr::mutate(
      log_r_country = dplyr::case_when(
        is.finite(log_r_country) ~ log_r_country,
        is.finite(log_r_group) ~ log_r_group,
        TRUE ~ log_r_global
      )
    ) %>%
    dplyr::select(country, log_r_country, n_country)

  # --- assembly -------------------------------------------------------------
  cell <- tidyr::crossing(
    country = unique(units$country),
    admin_level = valid_levels
  ) %>%
    dplyr::left_join(own_cell, by = c("country", "admin_level")) %>%
    dplyr::left_join(country_tbl, by = "country") %>%
    dplyr::mutate(
      # unname(): `level_factor` is a named vector, so `[admin_level]` carries
      # the names through into log_r_parent and on into the `mean_rate` COLUMN.
      # The size original (provenance.R) produces an unnamed column, and a named
      # column is a latent nuisance regardless -- it survives as.data.frame() and
      # shows up in CSV/print output for no reason.
      lvl_fac = unname(level_factor[admin_level]),
      log_r_parent = log_r_country + lvl_fac,
      # `n_cell` counts observations under `geometric` and EVENTS under
      # `pooled`; `k` is `sigma2_within/tau2` and `1/tau2` respectively. The
      # expression is shared because the two parameterisations were chosen to
      # make it so -- see `.anchor_engine()`.
      w_shrink = ifelse(is.na(n_cell), NA_real_, n_cell / (n_cell + k)),
      mr_tier = dplyr::case_when(
        is.na(n_cell) ~ "T2",
        n_cell <= 0 ~ "T1-zero",
        TRUE ~ "T1"
      ),
      log_mean_rate = dplyr::case_when(
        is.na(n_cell) ~ log_r_parent,
        n_cell <= 0 ~ log_r_parent,
        TRUE ~ w_shrink * log_r_raw + (1 - w_shrink) * log_r_parent
      ),
      mean_rate = exp(log_mean_rate),
      mr_parent = exp(log_r_parent),
      estimator = estimator,
      k = k,
      level_factor = lvl_fac
    ) %>%
    dplyr::select(
      country,
      admin_level,
      mean_rate,
      mr_tier,
      n_cell,
      w_shrink,
      mr_raw,
      mr_parent,
      estimator,
      k,
      level_factor
    )

  group <- tidyr::crossing(
    enso_group = unique(group_tbl$enso_group),
    admin_level = valid_levels
  ) %>%
    dplyr::left_join(group_level_tbl, by = c("enso_group", "admin_level")) %>%
    dplyr::left_join(group_tbl, by = "enso_group") %>%
    dplyr::mutate(
      lvl_fac = unname(level_factor[admin_level]), # see cell table above
      log_mean_rate = ifelse(!is.na(log_r), log_r, log_r_group + lvl_fac),
      mean_rate = exp(log_mean_rate),
      mr_tier = ifelse(!is.na(log_r), "T3-own", "T3-parent")
    ) %>%
    dplyr::select(enso_group, admin_level, mean_rate, mr_tier, n_group)

  global <- tibble::tibble(
    admin_level = valid_levels,
    mean_rate = exp(log_r_global + level_factor[valid_levels]),
    mr_tier = "T4"
  )

  meta <- list(
    track = track,
    column = spec$column,
    estimator = estimator,
    k = k,
    tau2 = eng$tau2,
    level_factor = level_factor,
    level_factor_ci = lfr$level_factor_ci,
    n_rows = nrow(d),
    n_cells = nrow(own_cell),
    n_countries = length(unique(units$country)),
    country_group = grouping,
    build_time = Sys.time()
  )
  # occ/dur carry no estimand stamp, so when their views re-keyed admin_level
  # (occurrence/duration.canonical_admin_level) the table says so here, for the
  # replay-side join (.resolve_frame_anchor(), utils/simulator_core.R). Only set
  # when TRUE, so every existing occ/dur table is unchanged.
  if (
    track %in%
      c("occ", "dur") &&
      isTRUE(
        (opt[[c(occ = "occurrence", dur = "duration")[[track]]]] %||%
          list())$canonical_admin_level
      )
  ) {
    meta$canonical_admin_level <- TRUE
  }

  if (isTRUE(spec$has_pop)) {
    # Digest the POPULATION, not the denominator. For `siz` the two are the same
    # vector, but `mon`'s `den` is `pop * pmax(exposure, 1)` -- a compound that
    # moves when the episode windows move. Stamping that would make the Step-7
    # ordering guard (.resolve_frame_anchor(), simulator_core.R) compare a
    # population fingerprint on the frame side against a population-times-time
    # fingerprint on the table side, i.e. a check that can never agree and
    # therefore never catch anything. `.mon_anchor_view()`/`.siz_anchor_view()` both
    # emit `pop`; the `%||% d$den` keeps any view that does not.
    pop_prov <- pop_prov_in %||%
      pop_provenance_stamp(d[["pop"]] %||% d$den, pop_source = "unstamped")
    meta$pop_provenance <- pop_prov
  }
  if (identical(track, "mon")) {
    # Structural, not an attribute -- same rationale as the `siz` branch below.
    # This is the side of the seam load_monthly_data()'s own stamp is checked
    # against, so a mon anchor built on the footprint and a mon fit declaring
    # `gid` cannot pass each other silently.
    meta$mon_estimand <- monthly_estimand_stamp(
      opt,
      extra = list(
        use_level_factor = isTRUE(cfg$level_factor %||% TRUE),
        levels_present = paste(valid_levels, collapse = "/"),
        # Which admin_level column this table's cells are keyed on. Recorded
        # because the replay side (.resolve_frame_anchor(), utils/simulator_core.R)
        # must join on the SAME convention -- the frame carries both the corpus
        # label and the canonical one, and picking the wrong column silently
        # resolves a different cell. Same expression `.mon_anchor_view()` gates on.
        canonical_admin_level = isTRUE(
          (opt$monthly %||% list())$canonical_admin_level %||%
            !identical(resolve_monthly_estimand(opt)$denominator, "footprint")
        ),
        n_rows = nrow(d)
      )
    )
  }
  if (identical(track, "siz")) {
    # Structural, not an attribute -- `summarise()` drops attributes, and this
    # table is the side of the seam most likely to be reached via a rebuilt
    # object. Mirrors provenance.R:997-1003.
    meta$size_estimand <- build_size_estimand_stamp(
      opt,
      extra = list(
        use_level_factor = isTRUE(cfg$level_factor %||% TRUE),
        levels_present = paste(valid_levels, collapse = "/"),
        # See the mon branch above -- the replay-side join key.
        canonical_admin_level = isTRUE(
          (opt$size %||% list())$canonical_admin_level %||%
            !identical(resolve_size_estimand(opt)$denominator, "footprint")
        ),
        n_rows = nrow(d)
      )
    )
  }

  brt <- list(cell = cell, group = group, global = global, meta = meta)
  if (isTRUE(spec$has_pop)) {
    attach_pop_provenance(brt, meta$pop_provenance)
  } else {
    brt
  }
}


# ---------------------------------------------------------------------------
# resolve_base_rate / attach_base_rate / persist_base_rate_table
# ---------------------------------------------------------------------------

#' The ONLY lookup path a consumer should use to obtain a baseline-rate anchor.
#'
#' A thin delegate to `resolve_mean_rate()` (provenance.R:1025). This is not
#' laziness: `build_base_rate_table()` emits the SAME `cell`/`group`/`global`/`meta`
#' schema as `build_mean_rate_table()` precisely so the T1->T5 fallback ladder does
#' not have to exist twice. Two copies of a five-tier fallback chain is exactly
#' the kind of duplication that drifts silently -- fit time and simulate time
#' already depend on resolving identically, on all four tracks.
#'
#' @param country character vector of ISO3 codes.
#' @param admin_level character vector, same length.
#' @param brt list from `build_base_rate_table()` (or its persisted sidecar).
#' @return tibble(base_rate, br_tier), one row per input pair, in input order.
resolve_base_rate <- function(country, admin_level, brt) {
  out <- resolve_mean_rate(country, admin_level, brt)
  tibble::tibble(base_rate = out$mean_rate, br_tier = out$mr_tier)
}

#' Join a track's resolved anchor onto a model frame.
#'
#' Adds `spec$column` (and `spec$tier_col`) and stashes the table itself on the
#' frame as `attr(d, spec$attr_name)`, mirroring `siz_data.R:225-230`, so
#' `persist_base_rate_table()` can reconcile in-memory against on-disk without the
#' caller having to carry the table separately. The attribute is
#' `"base_rate_table"` on occ/dur/mon and `"mean_rate_table"` on siz, which is
#' what keeps provenance.R's `persist_mean_rate_table()` and
#' `tests/gates_estimand_coherence.R:156` working against a frame the shared
#' builder produced.
#'
#' @param d model frame (must carry `country` and `admin_level`).
#' @param opt options list.
#' @param track "occ" | "dur" | "mon" | "siz".
#' @param brt prebuilt table; built here when NULL.
attach_base_rate <- function(d, opt, track, brt = NULL) {
  spec <- anchor_spec(track)
  if (is.null(brt)) {
    brt <- build_base_rate_table(opt, track)
  }
  res <- resolve_base_rate(
    as.character(d$country),
    as.character(d$admin_level),
    brt
  )
  d[[spec$column]] <- res$base_rate
  d[[spec$tier_col]] <- res$br_tier
  attr(d, spec$attr_name) <- brt
  d
}

#' Persist a track's anchor table and assert it reconciles against the frame.
#'
#' Generalises `persist_mean_rate_table()` (provenance.R:1076). The sidecar path is
#' run-INVARIANT by construction (`make_run_file(opt, spec$kind)` -> the
#' `*_base_rate` kinds, which carry no country/group/family/lag tokens), so every
#' fit config that shares the same estimator -- and, where the denominator
#' contains population, the same population digest -- resolves to the identical
#' file and none of them recomputes a second copy.
#'
#' The `file.exists()` skip is safe for the same reason it is safe for size: a
#' given path's content is immutable, and any input that would change the content
#' is in the path, so a change MISSES the cache and rebuilds rather than silently
#' reusing a stale table.
#'
#' @param opt options list.
#' @param track "occ" | "dur" | "mon" | "siz".
#' @param data frame returned by the track's loader, carrying the
#'   `base_rate_table` attribute and the resolved anchor column.
#' @return invisibly, the sidecar path.
persist_base_rate_table <- function(opt, track, data) {
  spec <- anchor_spec(track)
  brt <- attr(data, spec$attr_name)
  if (is.null(brt)) {
    stop(
      "persist_base_rate_table(",
      track,
      "): `data` carries no ",
      spec$attr_name,
      " attribute -- was the anchor attached by attach_base_rate()?",
      call. = FALSE
    )
  }

  path <- make_run_file(opt, kind = spec$kind)
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(brt, path)
    message(spec$kind, " (global table) -> ", path)
  }
  on_disk <- readRDS(path)

  # Provenance BEFORE arithmetic: the numeric reconciliation below DOES detect a
  # stale sidecar, but it reports the symptom and leaves the operator to infer
  # the cause. Checking the stamp first names the cause, and the rebuild order.
  if (isTRUE(spec$has_pop)) {
    assert_pop_provenance(
      data_stamp = pop_provenance_of(brt),
      sidecar_stamp = pop_provenance_of(on_disk),
      what = paste0("the ", spec$kind, " sidecar at ", path)
    )
  }

  check <- resolve_base_rate(
    as.character(data$country),
    as.character(data$admin_level),
    on_disk
  )
  if (
    !isTRUE(all.equal(check$base_rate, data[[spec$column]], tolerance = 1e-8))
  ) {
    stop(
      "persist_base_rate_table(",
      track,
      "): reconciliation FAILED -- data$",
      spec$column,
      " does not match a fresh lookup against the persisted ",
      "sidecar at ",
      path,
      ". The in-memory table used to build this fit's ",
      "data has drifted from what is on disk; delete the stale sidecar and ",
      "re-run.",
      call. = FALSE
    )
  }
  invisible(path)
}

#' Resolve the `base_rate_table_paths` argument for `build_location_frame()`.
#'
#' The simulator / smoothing / figure entry points all face the same question:
#' "which anchor sidecars, if any, does this run need?". Answering it per call
#' site is exactly how the `mean_rate_table_path` wiring ended up duplicated
#' across six files, so it is answered once here.
#'
#' A track contributes a path only when BOTH hold:
#'   * an `opt` for it was supplied (i.e. that track is part of the run), and
#'   * its fit was made with the anchor on (`<track>.offset_baserate: true`).
#' Everything else yields NULL for that track, and a fully-NULL result is
#' returned as NULL so `build_location_frame()` takes its inert path.
#'
#' An explicit `simulator.location_frame.base_rate_table_paths.<track>` entry
#' always wins, mirroring how `lf_opts$mean_rate_table_path` overrides the
#' derived size path. A track that is anchored but whose sidecar is missing is
#' a hard error here rather than a silent NULL -- `build_location_frame()`
#' would otherwise hand the kernel an un-anchored frame for an anchored fit,
#' which is a silent multiplicative bias, not a crash.
#'
#' @param track_opts named list of per-track option lists, e.g.
#'   `list(occ = opt_occ, dur = opt_dur, mon = opt_mon)`. A NULL entry means
#'   that track is not part of this run.
#' @param lf_opts the `simulator$location_frame` config block (may be NULL).
#' @return a named list of paths for the anchored tracks, or NULL when none.
resolve_base_rate_paths <- function(track_opts, lf_opts = NULL) {
  overrides <- (lf_opts %||% list())$base_rate_table_paths %||% list()
  out <- list()
  for (tr in c("occ", "dur", "mon")) {
    o <- track_opts[[tr]]
    if (is.null(o)) {
      next
    }
    if (!anchor_enabled(o, tr)) {
      next
    }
    p <- overrides[[tr]]
    if (is.null(p) || !nzchar(p)) {
      p <- make_run_file(o, kind = anchor_spec(tr)$kind)
    }
    if (!file.exists(p)) {
      stop(
        "resolve_base_rate_paths: the ",
        tr,
        " fit was made with ",
        tr,
        ".offset_baserate: true but its anchor sidecar is missing (",
        p,
        "). Re-run that track's 03_prepare_* step, or point ",
        "simulator.location_frame.base_rate_table_paths.",
        tr,
        " at the right file.",
        call. = FALSE
      )
    }
    out[[tr]] <- p
  }
  if (length(out) == 0L) NULL else out
}
