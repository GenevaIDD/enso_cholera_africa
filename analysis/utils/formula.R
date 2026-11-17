# formula.R — model-term assembly (fixed, random, crossbasis and offset terms)
#   for every track.
#
# Provides:
#   build_formula_terms()                          occurrence (brms)
#   build_duration_formula_parts(), build_duration_formula()
#   build_monthly_formula_parts(), build_monthly_formula()
#   build_size_formula_parts(), build_size_formula()
#   integrate_cb_terms()                           crossbasis terms
#   re_group_by_level(), nested_slope_re_terms(), filter_re_terms_for_data(),
#   filter_fixed_factor_terms_for_data()           RE / fixed-term helpers
#   validate_size_options(), resolve_slopes()      slope / convention resolution
#
# Details:
#   Internal but shared: .OFFSET_SPEC / .offset_terms() / .anchor_gate_on()
#   (offsets + base-rate anchor gate), .occ_basehaz() (also read by priors.R).
#   Callers source base.R; at call time paths.R (nbr_covariate_name(),
#   .leaf_scope_tag() when the leaf scope is not "all") and config.R
#   (get_default_options() when a size / monthly block has no backend).
#   Golden gate: analysis/tests/test_formula_builder_golden.R.

# ---------------------------------------------------------------------------
# Occurrence basehaz label (shared with utils/priors.R)
# ---------------------------------------------------------------------------

#' .occ_basehaz
#'
#' Resolve an occurrence brms/brms_cb backend block's `basehaz_poisson` to its
#' INTERNAL label, mapping the preferred config spelling "seasonal" onto the
#' identical-in-every-respect "cyclic" (the same `s(t_start_seas, bs = 'cc',
#' k = 8)` spline; see 04_fit_occ_brms.R's header). Called by
#' 04_fit_occ_brms.R and 04_fit_occ_brms_cb.R and by .has_smooth_terms() (the
#' class-"sds" prior gate). That last consumer is why this is a function and
#' not a copy: the prior gate has to agree with
#' the fit script about whether a smooth is in the formula, and a `set_prior(
#' class = "sds")` emitted against a formula that turned out to have no smooth
#' is a hard brms error, not a silently-ignored prior.
#'
#' Deliberately NOT folded in here:
#'   - `.validate_basehaz()` (line ~95) reads the RAW value per backend and
#'     needs the un-normalised string for its error messages.
#'   - `make_run_file()`'s basehaz token (utils/paths.R) tokenises the RAW
#'     value: "seasonal" and "cyclic" are the same model but must stay distinct
#'     artefact names, since configs in the repo use both spellings.
#'   - build_formula_terms()'s / build_formula_terms_legacy()'s
#'     `re_region_country_basehaz_slope` gate (formula.R / formula_legacy.R)
#'     dispatches on the BACKEND from `opt` rather than taking a
#'     resolved block, and only ever compares against "factor", for which
#'     normalisation is a no-op.
#'
#' @param brms_cfg the occurrence backend block the caller already resolved:
#'   `opt$occurrence$brms`, or `opt$occurrence$brms_cb %||% opt$occurrence$brms`
#'   at a `_cb` call site. NULL is treated as an empty block (-> the default).
#' @return character(1): "factor" | "cyclic" | "ms"
.occ_basehaz <- function(brms_cfg) {
  raw <- (brms_cfg %||% list())$basehaz_poisson %||% "seasonal"
  if (identical(raw, "seasonal")) "cyclic" else raw
}

#' .occ_baseline_block
#'
#' The occurrence brms baseline-hazard term for a resolved `basehaz`
#' (`.occ_basehaz()`): month factor, cyclic seasonal spline, or the M-spline
#' columns built by prepare_occ_brms_data(). Single definition for the
#' 04_fit_occ_brms*.R workers.
#' @param basehaz  "factor" | "cyclic" | "ms"
#' @param ms_basis the M-spline basis matrix (only read when basehaz == "ms")
#' @return character(1) formula fragment
.occ_baseline_block <- function(basehaz, ms_basis = NULL) {
  switch(
    basehaz,
    factor = "factor(months)",
    cyclic = "s(t_start_seas, bs = 'cc', k = 8)",
    ms = paste(colnames(ms_basis), collapse = " + ")
  )
}


#' integrate_cb_terms
#'
#' Fold DLNM basis columns (build_cb_block()) into occurrence formula terms:
#' drop the linear FE + RE of every dlnm variable with keep_linear_term: false,
#' append the basis columns of `linear_vars` as fixed effects, and add nested
#' region/country/ADM1 slope REs for those with re_adm1_slope: true. Lifted
#' verbatim from 04_fit_occ_brms_cb.R; also used by 04_fit_occ_brms_kl.R, whose
#' kernel variable is excluded from `linear_vars` because its columns enter the
#' non-linear kernel term instead.
#'
#' @param terms       build_formula_terms() output
#' @param dlnm_cfg    occurrence `dlnm:` block
#' @param cb_cols     all basis columns, in build order (build_cb_block()$cb_cols)
#' @param cb_meta     build_cb_block()$cb_meta
#' @param re_sep      "|" or "||"
#' @param drop_region TRUE in grouping mode (slopes collapse to country/ADM1)
#' @param linear_vars dlnm variables whose columns enter the linear predictor
#' @return `terms` with fixed_terms / re_terms updated
integrate_cb_terms <- function(
  terms,
  dlnm_cfg,
  cb_cols,
  cb_meta,
  re_sep,
  drop_region,
  linear_vars = names(cb_meta)
) {
  if (length(cb_cols) == 0) {
    return(terms)
  }
  dlnm_defaults <- dlnm_cfg$defaults %||% list()

  # Drop linear FE + RE for CB variables with keep_linear_term: false
  for (v in names(dlnm_cfg$variables)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    klt <- var_spec$keep_linear_term %||%
      dlnm_defaults$keep_linear_term %||%
      FALSE
    if (!isTRUE(klt)) {
      terms$fixed_terms <- terms$fixed_terms[terms$fixed_terms != v]
      terms$re_terms <- terms$re_terms[
        !grepl(sprintf("^\\(%s[[:space:]]", v), terms$re_terms)
      ]
    }
  }

  # Append CB fixed-effect terms (build order kept)
  lin_cols <- unlist(lapply(linear_vars, function(v) cb_meta[[v]]$cb_col_names))
  terms$fixed_terms <- c(terms$fixed_terms, cb_cols[cb_cols %in% lin_cols])

  # Append CB RE terms -- hierarchical (region/country/ADM1) slopes, opt-in per
  # variable via re_adm1_slope: true. Grouping mode (drop_region) collapses
  # them to country/ADM1.
  cb_groups <- re_group_by_level("gid", drop_region = isTRUE(drop_region))
  for (v in intersect(names(dlnm_cfg$variables), linear_vars)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    re_slope <- var_spec$re_adm1_slope %||%
      dlnm_defaults$re_adm1_slope %||%
      FALSE
    if (!isTRUE(re_slope)) {
      next
    }
    v_cb_names <- cb_meta[[v]]$cb_col_names
    if (is.null(v_cb_names)) {
      next
    }
    terms$re_terms <- c(
      terms$re_terms,
      nested_slope_re_terms(v_cb_names, cb_groups, re_sep)
    )
  }
  terms
}

# ---------------------------------------------------------------------------
# Formula builders: random-effect / fixed-term helpers + occurrence
# ---------------------------------------------------------------------------

#' re_group_by_level
#'
#' Per-level brms RE group strings using the ALWAYS-full-nesting convention
#' `AFRO_region / country / ADM1 [/ <adm2_col>]`.  Unlike .occ_group_by_level()
#' this does NOT collapse for single-country / single-region data: every level
#' is emitted so occurrence and size tracks share one resolver and slopes are
#' fully hierarchical.  The ADM2 leaf column differs by track (occurrence uses
#' `gid`, size uses `ADM2`).
#'
#' When `drop_region = TRUE` (country-grouping mode) the AFRO_region level is
#' removed entirely and the hierarchy roots at `country` (no group-level RE is
#' added — the fitted group *is* the model scope). The per-level loops in the
#' formula builders already `next` on a NULL group, so a NULL `region` cleanly
#' drops all region-level intercept and slope terms.
#'
#' @param adm2_col character — leaf column for the ADM2 group ("gid" | "ADM2")
#' @param drop_region logical — drop the AFRO_region top level (grouping mode)
#' @return named list(region, country, adm1, adm2); region is NULL when dropped
re_group_by_level <- function(adm2_col = "gid", drop_region = FALSE) {
  if (isTRUE(drop_region)) {
    list(
      region = NULL,
      country = "country",
      adm1 = "country:ADM1",
      adm2 = sprintf("country:ADM1:%s", adm2_col)
    )
  } else {
    list(
      region = "AFRO_region",
      country = "AFRO_region:country",
      adm1 = "AFRO_region:country:ADM1",
      adm2 = sprintf("AFRO_region:country:ADM1:%s", adm2_col)
    )
  }
}


#' nested_slope_re_terms
#'
#' Build hierarchical (multilevel) random-slope RE term strings for a set of
#' columns, emitting one `(col - 1 <re_sep> <group>)` term per column × level.
#' Together these give a separate variance component at each nesting level
#' (partial pooling region -> country -> ADM1).  Used for DLNM crossbasis
#' columns so both tracks build CB slopes identically.
#'
#' @param cols    character vector of column names (e.g. crossbasis columns)
#' @param groups  named list from re_group_by_level()
#' @param re_sep  "||" (uncorrelated) or "|"
#' @param levels  character vector of group levels to emit
#' @return character vector of RE term strings
nested_slope_re_terms <- function(
  cols,
  groups,
  re_sep,
  levels = c("region", "country", "adm1")
) {
  if (length(cols) == 0) {
    return(character(0))
  }
  terms <- character(0)
  for (cn in cols) {
    for (lvl in levels) {
      grp <- groups[[lvl]]
      if (is.null(grp)) {
        next
      }
      terms <- c(terms, sprintf("(%s - 1 %s %s)", cn, re_sep, grp))
    }
  }
  terms
}


#' filter_re_terms_for_data
#'
#' Drops or reshapes RE terms whose grouping columns have fewer than 2 distinct
#' values in `dat`.
#'   - Slash (`country/ADM1`, intercept nesting): each column needs >= 2 levels;
#'     a single-level column is removed from the nesting (`(1 | country/ADM1)`
#'     -> `(1 | ADM1)` for one country); with none left the term is dropped.
#'   - Colon (`country:ADM1`, slope groups): counts distinct combinations; a
#'     single-level compound key is dropped, otherwise kept.
#' `keep` columns are never trimmed (.single_country_keep() keeps `country` in a
#' one-country grouped fit, so group names stay what the simulator and figure
#' extractors map).
#'
#' @param re_terms character vector of "(... | group)" RE strings
#' @param dat      data.frame on which the model will be fit
#' @param keep     grouping columns exempt from the <2-levels rule
#' @return filtered (and possibly trimmed) character vector
filter_re_terms_for_data <- function(re_terms, dat, keep = character(0)) {
  if (length(re_terms) == 0L) {
    return(re_terms)
  }

  result <- character(0)

  for (term in re_terms) {
    # Extract group string: after the last "|" or "||", before the closing ")"
    grp_str <- sub(".*\\|\\|?\\s+", "", term)
    grp_str <- trimws(sub("\\s*\\)\\s*$", "", grp_str))

    if (grepl("/", grp_str, fixed = TRUE)) {
      # Slash notation: intercept nesting hierarchy.
      # Each component is its own column; check individually.
      components <- strsplit(grp_str, "/", fixed = TRUE)[[1]]
      valid <- vapply(
        components,
        function(comp) {
          if (!comp %in% names(dat)) {
            return(TRUE)
          } # unknown col — keep
          if (comp %in% keep) {
            return(TRUE)
          } # exempt (single-country fit)
          dplyr::n_distinct(dat[[comp]], na.rm = TRUE) >= 2L
        },
        logical(1)
      )

      if (!any(valid)) {
        message(sprintf(
          "[filter_re] Dropping intercept RE (all nesting levels have <2 distinct values): %s",
          term
        ))
        next
      }
      kept <- components[valid]
      if (length(kept) < length(components)) {
        dropped <- components[!valid]
        new_grp <- paste(kept, collapse = "/")
        new_term <- sub(grp_str, new_grp, term, fixed = TRUE)
        message(sprintf(
          "[filter_re] Trimmed intercept RE (removed sparse level(s): %s): %s  ->  %s",
          paste(dropped, collapse = ", "),
          term,
          new_term
        ))
        result <- c(result, new_term)
      } else {
        result <- c(result, term)
      }
    } else {
      # Colon notation: slope RE compound group key (or single variable).
      # Drop the term only when the *compound* key has <2 distinct values.
      cols <- strsplit(grp_str, ":", fixed = TRUE)[[1]]
      known_cols <- cols[cols %in% names(dat)]
      if (length(known_cols) == 0L || all(known_cols %in% keep)) {
        result <- c(result, term) # cannot check, or exempt — keep
        next
      }
      combo <- do.call(
        paste,
        c(lapply(known_cols, function(cc) dat[[cc]]), list(sep = ":"))
      )
      if (dplyr::n_distinct(combo, na.rm = TRUE) < 2L) {
        message(sprintf(
          "[filter_re] Dropping slope RE (group has <2 distinct values): %s",
          term
        ))
        next
      }
      result <- c(result, term)
    }
  }

  result
}


#' filter_fixed_factor_terms_for_data
#'
#' Drop fixed-effect terms whose data column is a factor (or character) with
#' fewer than 2 distinct levels in `dat`, and cascade that removal to any RE
#' slope terms that reference the same column.
#'
#' Context: `time_to_last_outbreak` is stored as a factor (bins from
#' `make_time_to_last_outbreak_cat()`).  In a country-grouped fit whose subset
#' contains outbreaks from only one time-to-last-outbreak bin (e.g. all
#' Southern-Africa outbreaks in the "[6,Inf]" category), the factor has a
#' single level and brms raises "contrasts can be applied only to factors with
#' 2 or more levels" when constructing the model matrix.
#'
#' Cascade rule: any RE slope term whose LHS predictor string contains the
#' dropped column name is also removed.  E.g. dropping "time_to_last_outbreak"
#' removes `(time_to_last_outbreak - 1 || country)` and
#' `(time_to_last_outbreak - 1 || country:ADM1)`.
#'
#' @param fixed_terms character vector of fixed-effect column names
#' @param re_terms    character vector of "(... | group)" RE strings
#' @param dat         data.frame on which the model will be fit
#' @return list(fixed_terms, re_terms) after filtering
filter_fixed_factor_terms_for_data <- function(fixed_terms, re_terms, dat) {
  if (length(fixed_terms) == 0L) {
    return(list(fixed_terms = fixed_terms, re_terms = re_terms))
  }

  keep <- vapply(
    fixed_terms,
    function(col) {
      if (!col %in% names(dat)) {
        return(TRUE)
      } # unknown col — keep
      v <- dat[[col]]
      if (!is.factor(v) && !is.character(v)) {
        return(TRUE)
      } # numeric — keep
      dplyr::n_distinct(v, na.rm = TRUE) >= 2L
    },
    logical(1)
  )

  dropped <- fixed_terms[!keep]
  if (length(dropped) > 0L) {
    for (col in dropped) {
      message(sprintf(
        "[filter_fixed] Dropping fixed term '%s': factor/character with <2 distinct levels.",
        col
      ))
    }
    # Cascade: drop RE slope terms whose LHS references the dropped column.
    if (length(re_terms) > 0L) {
      drop_re <- vapply(
        re_terms,
        function(t) {
          any(vapply(
            dropped,
            function(col) grepl(col, t, fixed = TRUE),
            logical(1)
          ))
        },
        logical(1)
      )
      if (any(drop_re)) {
        for (t in re_terms[drop_re]) {
          message(sprintf(
            "[filter_fixed] Dropping RE term (references dropped fixed term): %s",
            t
          ))
        }
        re_terms <- re_terms[!drop_re]
      }
    }
  }

  list(fixed_terms = fixed_terms[keep], re_terms = re_terms)
}


# ============================================================================
# Shared formula-builder helpers
# ============================================================================
#
# Extracted from build_formula_terms() (occurrence), build_size_formula_parts()
# and build_duration_formula_parts(): the internals below were duplicated
# near-verbatim across two or all three builders. Only the internals are
# merged here -- each builder keeps its own signature, return shape, and
# track-specific gates (calendar escape hatches, ENSO/spillover special
# cases, seasonal terms, baseline hazards, pop-elasticity, phi_re). The golden
# characterization test (analysis/tests/test_formula_builder_golden.R) is the
# byte-identical guarantee that motivated extracting these instead of leaving
# the duplication in place -- it must stay green after every one of these
# helpers is wired in.

#' .resolve_fixed_covariates
#'
#' Canonicalize a track's random-effects flags into per-covariate slope
#' levels (via resolve_slopes()), then derive the fixed-effect covariate set
#' as the union of (a) explicitly-toggled covariates and (b) slope-implied
#' covariates (Convention 2). Shared by the size and duration builders, whose
#' fixed-covariate resolution is line-for-line identical; occurrence's is not
#' (explicit skip_covs, calendar escape hatches, no "1" prepend) and is left
#' as its own block in build_formula_terms().
#'
#' @param block       opt$size or opt$duration (must have $random_effects,
#'                    $covariates)
#' @param admin_level passed through to resolve_slopes()
#' @param enso_lag    opt$enso$lag_months; only read when block$covariates$enso
#'                    is TRUE
#' @return list(slopes, slope_covs, fixed_covs) -- fixed_covs excludes the
#'   leading "1" intercept term; callers prepend it.
.resolve_fixed_covariates <- function(block, admin_level, enso_lag) {
  slopes <- resolve_slopes(block$random_effects, admin_level)
  slope_covs <- names(slopes)

  enso_on <- isTRUE(block$covariates$enso)
  if (!enso_on) {
    slope_covs <- slope_covs[!grepl("^enso_lag\\d+$", slope_covs)]
  }

  covs_map <- block$covariates %||% list()
  covs_on <- names(covs_map)[vapply(
    names(covs_map),
    function(k) isTRUE(covs_map[[k]]),
    logical(1)
  )]
  covs_on <- covs_on[covs_on != "enso"]
  covs_on <- sub(
    "^time_to_last_outbreak_cat$",
    "time_to_last_outbreak",
    covs_on
  )
  fixed_covs <- unique(c(covs_on, slope_covs))

  if (enso_on) {
    fixed_covs <- unique(c(
      fixed_covs,
      sprintf("enso_lag%d", as.integer(enso_lag))
    ))
  }

  list(slopes = slopes, slope_covs = slope_covs, fixed_covs = fixed_covs)
}

#' .slope_re_terms
#'
#' The per-(covariate, level) RE-term loop: for each covariate in `covs`,
#' emit "(cov - 1 <sep> group)" for each level in `levels` whose slope flag
#' (slopes[[cov]][[lvl]]) is TRUE and whose group string (gbl[[lvl]]) is not
#' NULL (NULL means the level was dropped, e.g. AFRO_region under grouping
#' mode). Shared by all three builders' per-covariate slope-RE loops,
#' including occurrence's single-covariate time_to_last_outbreak slope (pass
#' a length-1 `covs`).
#'
#' @param covs   character vector of covariate names to loop over
#' @param slopes named list from resolve_slopes() (or equivalent), keyed by
#'               covariate then level
#' @param gbl    named list of group strings from re_group_by_level()
#' @param re_sep "|" or "||"
#' @param levels character vector, the level names to check (subset of
#'               names(gbl))
#' @return character vector of RE term strings (possibly empty)
.slope_re_terms <- function(covs, slopes, gbl, re_sep, levels) {
  out <- character(0)
  for (cov in covs) {
    for (lvl in levels) {
      if (!isTRUE(slopes[[cov]][[lvl]])) {
        next
      }
      grp <- gbl[[lvl]]
      if (is.null(grp)) {
        next
      }
      out <- c(out, sprintf("(%s - 1 %s %s)", cov, re_sep, grp))
    }
  }
  out
}

#' .re_terms_for_column
#'
#' Emit "(col - 1 <sep> group)" for every level in `levels` whose group
#' string is not NULL -- no per-level flag lookup (the caller has already
#' gated on a single flag before calling this). Shared by occurrence's
#' special-cased RE loops (ongoing outbreak, the two spillover lags, ENSO
#' shallow/deep, basehaz group-varying seasonality), which are all this
#' pattern behind a different upstream gate.
#'
#' @param col    data column name (single covariate)
#' @param levels character vector, the level names to check
#' @param gbl    named list of group strings from re_group_by_level()
#' @param re_sep "|" or "||"
#' @return character vector of RE term strings (possibly empty)
.re_terms_for_column <- function(col, levels, gbl, re_sep) {
  out <- character(0)
  for (lvl in levels) {
    grp <- gbl[[lvl]]
    if (is.null(grp)) {
      next
    }
    out <- c(out, sprintf("(%s - 1 %s %s)", col, re_sep, grp))
  }
  out
}

#' .intercept_re_term
#'
#' The re_region_country_adm1_intercept RE block, including the
#' re_adm2_intercept-forced-FALSE-under-ADM1 warning. Shared by all four
#' tracks (duration via .scoped_intercept_re_terms(); its re_adm2_intercept
#' defaults to FALSE, which keeps every pre-ADM1-ADM2 duration config at
#' the adm1 intercept).
#'
#' Preserves the separator asymmetry verbatim: occurrence passes its
#' resolved re_sep ("|" or "||"); size and duration hardcode "|" (harmless
#' for an intercept-only term, but a literal difference in the emitted
#' formula string that must not be "fixed" here -- see the golden test).
#'
#' @param re_flags    opt$occurrence$random_effects / opt$size$random_effects
#'                    / opt$duration$random_effects
#' @param admin_level the track's admin_level ("ADM1" | "ADM2")
#' @param gbl         named list of group strings from re_group_by_level()
#' @param re_sep      "|" or "||" (occurrence) / "|" (size, duration)
#' @param label       appended to the warning message (e.g. " (size model)");
#'                    "" for occurrence
#' @return character vector of length 0 or 1
.intercept_re_term <- function(
  re_flags,
  admin_level,
  gbl,
  re_sep,
  label = "",
  leaf_scope = "all"
) {
  if (!isTRUE(re_flags$re_region_country_adm1_intercept)) {
    return(character(0))
  }
  re_adm2i <- isTRUE(re_flags$re_adm2_intercept)
  if (re_adm2i && identical(admin_level, "ADM1")) {
    warning("re_adm2_intercept forced FALSE: admin_level == ADM1", label)
    re_adm2i <- FALSE
  }
  # `re_adm2_leaf_scope: adm2_only` (size / monthly / duration): the leaf exists only on
  # canonical-ADM2 units. Under the historical `all`, an ADM1-only unit's gid
  # also gets a leaf -- a singleton under its own ADM1 node, so two REs describe
  # one geography and only their sum is identified, while the singletons still
  # inform sd(leaf) as if they were ADM2 units. Here the nested intercept stops
  # at ADM1 and the leaf is a separate term whose design column `is_adm2` is 0
  # for ADM1-only rows (their single shared "none" level is prior-only). `gbl`
  # must then come from re_group_by_level("gid_leaf") -- see .leaf_columns().
  if (re_adm2i && identical(leaf_scope, "adm2_only")) {
    return(c(
      sprintf("(1 %s %s)", re_sep, gsub(":", "/", gbl$adm1)),
      sprintf("(0 + is_adm2 %s %s)", re_sep, gbl$adm2)
    ))
  }
  int_grp <- if (re_adm2i) {
    gsub(":", "/", gbl$adm2)
  } else {
    gsub(":", "/", gbl$adm1)
  }
  sprintf("(1 %s %s)", re_sep, int_grp)
}

#' .scoped_intercept_re_terms
#'
#' The intercept RE block of a track under its `re_adm2_leaf_scope`, shared by
#' the size / monthly builders (.build_count_track_parts()) and the duration
#' builder: under `adm2_only` the group strings come from
#' re_group_by_level("gid_leaf"), otherwise from the track's own `gbl`.
#'
#' @param cfg the track's config block (its random_effects are read).
#' @param admin_level the track's admin_level.
#' @param gbl group strings from re_group_by_level(<track leaf column>).
#' @param drop_region grouping mode (AFRO_region dropped).
#' @param label appended to .intercept_re_term()'s warning.
#' @return character vector of intercept RE terms.
.scoped_intercept_re_terms <- function(
  cfg,
  admin_level,
  gbl,
  drop_region,
  label = ""
) {
  leaf_scope <- .leaf_scope_of(cfg)
  int_gbl <- if (identical(leaf_scope, "adm2_only")) {
    re_group_by_level("gid_leaf", drop_region = drop_region)
  } else {
    gbl
  }
  .intercept_re_term(
    cfg$random_effects,
    admin_level,
    int_gbl,
    "|",
    label,
    leaf_scope = leaf_scope
  )
}

#' .leaf_scope_of
#'
#' `random_effects.re_adm2_leaf_scope` for a size / monthly / duration block. The default
#' "all" is returned without touching utils/paths.R (where the validating
#' `.leaf_scope_tag()` lives), so the default formula path keeps its existing
#' source closure.
.leaf_scope_of <- function(blk) {
  s <- ((blk %||% list())$random_effects %||% list())$re_adm2_leaf_scope %||%
    "all"
  if (identical(s, "all")) "all" else .leaf_scope_tag(blk)
}

#' .leaf_columns
#'
#' Data side of `re_adm2_leaf_scope: adm2_only`: `is_adm2` (1 on canonical
#' ADM2 rows) and `gid_leaf` (the gid on ADM2 rows, "none" elsewhere). Needs
#' the GEOMETRY-derived admin level -- the corpus label is exactly what
#' mislabels whole-state LPs "ADM2" -- so it refuses unless P3 has run
#' (`canonical_admin_level`) or the panel was prepped under
#' `episodes.key: gid_season`, whose admin_level is canonical by construction.
.leaf_columns <- function(d, canonical_ok, tag = "leaf") {
  if (!isTRUE(canonical_ok)) {
    stop(
      tag,
      ": random_effects.re_adm2_leaf_scope = 'adm2_only' needs the ",
      "geometry-derived admin level. Set canonical_admin_level: true (or ",
      "prep under --episode_key gid_season).",
      call. = FALSE
    )
  }
  d$is_adm2 <- as.integer(as.character(d$admin_level) == "ADM2")
  d$gid_leaf <- ifelse(d$is_adm2 == 1L, as.character(d$gid), "none")
  d
}


# ---------------------------------------------------------------------------
# Offsets (shared across all four tracks)
# ---------------------------------------------------------------------------
#
# All four tracks' offsets are built here, so each Intercept prior stays in
# step with its offsets (see anchor.R's header). Emission order is pop, then
# person-time, then the anchor: brms ignores it but the formula strings and
# their goldens do not. The anchor is a joined column in its own
# offset(log(...)) term, never folded into log_offset, which would hide it from
# pp_check / postprocessing and break the simulator's per-offset flags.
.OFFSET_SPEC <- list(
  # `time`       — emit the unconditional person-time offset `offset(log_offset)`
  # `pop_gate`   — config key under the track's block toggling `offset(log(pop))`
  # `anchor_gate`— config key toggling the empirical baseline-rate anchor
  # `anchor_col` — the joined column name (utils/anchor.R `.ANCHOR_SPEC$column`)
  occ = list(
    block = "occurrence",
    pop_gate = NULL,
    time = TRUE,
    anchor_gate = "offset_baserate",
    anchor_col = "base_rate_occ"
  ),
  siz = list(
    block = "size",
    pop_gate = "offset_pop",
    time = FALSE,
    anchor_gate = "offset_meanrate",
    anchor_col = "mean_rate"
  ),
  dur = list(
    block = "duration",
    pop_gate = NULL,
    time = TRUE,
    anchor_gate = "offset_baserate",
    anchor_col = "base_rate_dur"
  ),
  mon = list(
    block = "monthly",
    pop_gate = "offset_pop",
    time = TRUE,
    anchor_gate = "offset_baserate",
    anchor_col = "base_rate_mon"
  )
)

#' .offset_terms
#'
#' The offset block for one track, as brms formula strings.
#'
#' Size reaches its anchor through `size.offset_meanrate` (on by default, and on
#' in every config that exists); occurrence / duration / monthly reach theirs
#' through `<track>.offset_baserate` (off by default, so every pre-existing
#' config emits exactly the offsets it emitted before).
#'
#' Deliberately reads `opt` directly rather than calling `anchor_enabled()`:
#' `utils/config.R` is sourced by everything and must not depend on
#' `utils/anchor.R`, which itself depends on `paths.R` / `provenance.R` /
#' `occ_data.R`. The two gates are kept in step by
#' `analysis/tests/test_base_rate_table.R`, not by a source-order coupling.
#'
#' @param opt   parsed options list.
#' @param track "occ" | "siz" | "dur" | "mon".
#' @return character vector of `offset(...)` strings, possibly length 0.
.offset_terms <- function(opt, track) {
  spec <- .offset_spec(track)
  block <- opt[[spec$block]] %||% list()
  off_terms <- character(0)

  if (!is.null(spec$pop_gate) && isTRUE(block[[spec$pop_gate]])) {
    off_terms <- c(off_terms, "offset(log(pop))")
  }

  if (isTRUE(spec$time)) {
    off_terms <- c(off_terms, "offset(log_offset)")
  }

  if (isTRUE(block[[spec$anchor_gate]])) {
    off_terms <- c(off_terms, sprintf("offset(log(%s))", spec$anchor_col))
  }

  off_terms
}

#' Look up a track's offset spec, with an actionable error on a typo.
.offset_spec <- function(track) {
  spec <- .OFFSET_SPEC[[track]]
  if (is.null(spec)) {
    stop(
      ".offset_spec: unknown track '",
      paste(track, collapse = ","),
      "' (expected one of ",
      paste(names(.OFFSET_SPEC), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  spec
}

#' .anchor_gate_on
#'
#' Does this track use the unified residual-scale `Intercept` prior?
#'
#' Read by the prior builders, off by default on all four tracks, and keyed on
#' `<track>.offset_baserate` for all four -- including size, whose FORMULA gate
#' is the older, separate `size.offset_meanrate` (see `.OFFSET_SPEC`). That
#' split is deliberate and is the one place the formula gate and the prior gate
#' are allowed to differ:
#'
#'   - occurrence / duration / monthly had no anchor at all, so one flag has to
#'     turn on both the `offset(log(<anchor>))` term and the residual-scale
#'     intercept. Turning on only one of the two would reintroduce exactly the
#'     mis-calibration the anchor exists to remove, so they share `offset_baserate`.
#'   - size carries `offset(log(mean_rate))` independently of this
#'     mechanism, and its prior-predictive checks are already well-calibrated
#'     under `normal(-3, 1)`. Its offset is therefore NOT the thing being
#'     switched on; only the move to the unified `normal(0, 1)` is, and that is
#'     a numeric change to a deployed model. It gets its own opt-in so it is a
#'     deliberate choice, never a side-effect.
#'
#' @return logical(1)
.anchor_gate_on <- function(opt, track) {
  spec <- .offset_spec(track)
  block <- opt[[spec$block]] %||% list()
  on <- isTRUE(block$offset_baserate)

  # An anchored intercept with no anchor offset is the same double-counting bug
  # in the other direction: the prior would sit at residual scale while the
  # linear predictor still has to span the full offset scale.
  if (on && !isTRUE(block[[spec$anchor_gate]])) {
    stop(
      spec$block,
      ".offset_baserate is TRUE but ",
      spec$block,
      ".",
      spec$anchor_gate,
      " is not, so the Intercept prior would move to a ",
      "residual scale with no offset(log(",
      spec$anchor_col,
      ")) term in ",
      "the formula to justify it. Set both, or neither.",
      call. = FALSE
    )
  }
  on
}


#' build_formula_terms
#'
#' Occurrence formula terms via the size-model resolve_slopes() convention,
#' with explicit per-level RE terms (AFRO_region, AFRO_region:country, ...)
#' instead of `/`-nesting. Active slope RE flags imply their fixed effect;
#' opt$occurrence$covariates switches time_to_last_outbreak_cat,
#' ongoing_outbreak_lag1, enso, months and quarter. Special-cased: the intercept
#' REs, re_region_country_ongoing_slope, the ENSO slope families (need
#' covariates$enso) and re_region_country_basehaz_slope (factor basehaz only).
#' Returns list($fixed_terms, $re_terms ("(x - 1 | g)" strings)).
#'
#' @param opt         resolved options list (must have opt$occurrence)
#' @param n_country   integer, number of country levels in the data. Unread in
#'                    the body (kept for interface parity with
#'                    build_formula_terms_legacy(), which uses it).
#' @param n_region    integer, number of AFRO_region levels (see n_country)
#' @param dat         optional data frame. When supplied, applies
#'                    filter_fixed_factor_terms_for_data() (drops single-level
#'                    factor fixed effects and their RE slopes), then
#'                    filter_re_terms_for_data(). NULL (default) skips both, so
#'                    full nesting survives (test_occurrence_config_merge.R's
#'                    single-country case relies on it).
build_formula_terms <- function(opt, n_country, n_region, dat = NULL) {
  occ <- opt$occurrence
  covs <- occ$covariates %||% list()
  re_raw <- occ$random_effects %||% list()
  be <- occ$backend %||% "brms_cb"
  re_sep <- if (isTRUE(occ[[be]]$re_uncorrelated)) "||" else "|"

  # --- Resolve slopes via the canonical resolver (size convention) ---
  # Slopes cap at ADM1: occurrence adm2-level intercept is handled separately.
  slopes <- resolve_slopes(re_raw, "ADM1")

  # Covariates to skip from the generic climate loop:
  #   time_to_last_outbreak — FE comes from the escape-hatch, not slope implication
  #   ongoing               — column name mismatch (data col = ongoing_outbreak_lag1)
  #   enso_lag*             — dynamic column handled on its own ENSO branch
  #   nbr_contagion         — dynamic column (rho/scale-dependent), own branch.
  #                           This is the resolve_slopes()-decoded token from
  #                           re_{region_country,adm1}_nbr_contagion_slope, not
  #                           the covariates$ toggle name — same shape as
  #                           "ongoing" above.
  skip_covs <- c(
    "time_to_last_outbreak",
    "ongoing",
    "ongoing_country",
    "ongoing_tu",
    "nbr_contagion",
    grep("^enso_lag", names(slopes), value = TRUE)
  )
  climate_covs <- setdiff(names(slopes), skip_covs)

  # --- Fixed-effect terms ---
  fixed_terms <- character(0)
  if (isTRUE(covs$time_to_last_outbreak_cat)) {
    fixed_terms <- c(fixed_terms, "time_to_last_outbreak")
  }
  if (isTRUE(covs$ongoing_outbreak_lag1)) {
    fixed_terms <- c(fixed_terms, "ongoing_outbreak_lag1")
  }
  if (isTRUE(covs$ongoing_outbreak_country_lag1)) {
    fixed_terms <- c(fixed_terms, "ongoing_outbreak_country_lag1")
  }
  if (isTRUE(covs$ongoing_outbreak_tu_lag1)) {
    fixed_terms <- c(fixed_terms, "ongoing_outbreak_tu_lag1")
  }

  # Convention 2: active climate slope flags imply climate fixed effects
  fixed_terms <- c(fixed_terms, climate_covs)

  # Calendar escape hatch
  if (isTRUE(covs$months)) {
    fixed_terms <- c(fixed_terms, "months")
  }
  if (isTRUE(covs$quarter)) {
    fixed_terms <- c(fixed_terms, "quarter")
  }

  # ENSO fixed effect — driven by covariates$enso toggle
  enso_col <- NULL
  if (isTRUE(covs$enso)) {
    lag_n <- as.integer(opt$enso$lag_months)
    enso_col <- sprintf("enso_lag%d", lag_n)
    fixed_terms <- c(fixed_terms, enso_col)
  }

  # Distance-weighted neighbour contagion — driven by the nbr_contagion_lag1
  # toggle. Dynamic column name, exactly like ENSO: the rho and scale live in
  # the top-level contagion: block and must match the columns 03_prepare_occ.R
  # materialised (nbr_lag1_r<rho>[_rn]).
  nbr_col <- NULL
  if (isTRUE(covs$nbr_contagion_lag1)) {
    nbr_col <- nbr_covariate_name(
      opt$contagion$rho_km,
      opt$contagion$scale,
      opt$contagion$partition,
      opt$contagion$partition_source
    )
    fixed_terms <- c(fixed_terms, nbr_col)
  }

  # --- Per-level group strings (always full nesting; shared with size track) ---
  # Country-grouping mode drops the AFRO_region level (region-rooted terms vanish
  # because gbl$region becomes NULL and the per-level loops skip NULL groups).
  drop_region <- isTRUE(opt$grouping$enable)
  gbl <- re_group_by_level("gid", drop_region = drop_region)

  # --- RE terms ---
  re_terms <- character(0)

  # Intercept RE — keep nested-slash form (same as size model)
  re_terms <- c(
    re_terms,
    .intercept_re_term(re_raw, occ$admin_level, gbl, re_sep)
  )

  # Climate slope REs — explicit per-level terms (size convention)
  re_terms <- c(
    re_terms,
    .slope_re_terms(
      climate_covs,
      slopes,
      gbl,
      re_sep,
      c("region", "country", "adm1")
    )
  )

  # ttlo slope RE — explicit per-level (region + country), via resolve_slopes()
  re_terms <- c(
    re_terms,
    .slope_re_terms(
      "time_to_last_outbreak",
      slopes,
      gbl,
      re_sep,
      c("region", "country")
    )
  )

  # Ongoing slope RE — special-cased (data column = ongoing_outbreak_lag1)
  if (
    isTRUE(re_raw$re_region_country_ongoing_slope) &&
      isTRUE(covs$ongoing_outbreak_lag1)
  ) {
    re_terms <- c(
      re_terms,
      .re_terms_for_column(
        "ongoing_outbreak_lag1",
        c("region", "country"),
        gbl,
        re_sep
      )
    )
  }

  # Spillover slope REs — special-cased (column name mismatch with resolve_slopes() decoding)
  for (.new_col in c(
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1"
  )) {
    .flag <- if (.new_col == "ongoing_outbreak_country_lag1") {
      "re_region_country_ongoing_country_slope"
    } else {
      "re_region_country_ongoing_tu_slope"
    }
    if (isTRUE(re_raw[[.flag]]) && isTRUE(covs[[.new_col]])) {
      re_terms <- c(
        re_terms,
        .re_terms_for_column(.new_col, c("region", "country"), gbl, re_sep)
      )
    }
  }

  # ENSO slope REs — special-cased (dynamic column; requires covariates$enso = TRUE)
  if (!is.null(enso_col)) {
    lag_n <- as.integer(opt$enso$lag_months)
    shallow_key <- sprintf("re_region_country_enso_lag%d_slope", lag_n)
    deep_key <- sprintf("re_adm1_enso_lag%d_slope", lag_n)
    if (isTRUE(re_raw[[shallow_key]])) {
      re_terms <- c(
        re_terms,
        .re_terms_for_column(enso_col, c("region", "country"), gbl, re_sep)
      )
    }
    if (isTRUE(re_raw[[deep_key]])) {
      re_terms <- c(
        re_terms,
        .re_terms_for_column(
          enso_col,
          c("region", "country", "adm1"),
          gbl,
          re_sep
        )
      )
    }
  }

  # Neighbour-contagion slope REs — special-cased (dynamic column; requires
  # covariates$nbr_contagion_lag1 = TRUE). Unlike the ENSO pair above, the deep
  # flag adds ONLY the adm1 level, so setting both flags yields the shallow and
  # deep terms exactly once each rather than duplicating region/country.
  if (!is.null(nbr_col)) {
    if (isTRUE(re_raw$re_region_country_nbr_contagion_slope)) {
      re_terms <- c(
        re_terms,
        .re_terms_for_column(nbr_col, c("region", "country"), gbl, re_sep)
      )
    }
    if (isTRUE(re_raw$re_adm1_nbr_contagion_slope)) {
      re_terms <- c(
        re_terms,
        .re_terms_for_column(nbr_col, "adm1", gbl, re_sep)
      )
    }
  }

  # Basehaz group-varying seasonality — special-cased (factor basehaz only)
  occ_bh <- if (be == "brms") {
    occ$brms$basehaz_poisson %||% "seasonal"
  } else {
    NULL
  }
  if (
    isTRUE(re_raw$re_region_country_basehaz_slope) &&
      identical(occ_bh, "factor")
  ) {
    re_terms <- c(
      re_terms,
      .re_terms_for_column("months", c("region", "country"), gbl, re_sep)
    )
  }

  # Data-aware pruning — opt-in via dat, same order size/duration
  # use: fixed-factor filter first (its RE cascade must run before the RE
  # filter below sees the surviving terms), then the RE-group filter. NULL
  # (the default) skips both.
  if (!is.null(dat)) {
    fft <- filter_fixed_factor_terms_for_data(fixed_terms, re_terms, dat)
    fixed_terms <- fft$fixed_terms
    re_terms <- filter_re_terms_for_data(
      fft$re_terms,
      dat,
      keep = .single_country_keep(opt, dat)
    )
  }

  # Offsets — shared with the other three tracks (.offset_terms()). With the
  # anchor gate off the emitted string is "+ offset(log_offset)".
  off_terms <- .offset_terms(opt, "occ")

  .flag_single_country(
    list(fixed_terms = fixed_terms, re_terms = re_terms, off_terms = off_terms),
    opt,
    dat
  )
}

#' .single_country_keep
#'
#' Grouping columns filter_re_terms_for_data() must keep for this fit:
#' "country" when a grouped fit (grouping.enable) has exactly one country in
#' `dat`, else character(0). With one country the country-level REs are kept
#' (one level each) under a tight sd prior (priors.R::.SINGLE_COUNTRY_SD,
#' applied when the parts carry `single_country`), instead of being trimmed:
#' trimming renamed the nested intercept groups to `ADM1` / `ADM1:gid`, which
#' simulator_core.R::.parse_brms_re_draws() and .extract_phi_re_draws() do not
#' map, so those REs were dropped at replay without a warning.
#' @noRd
.single_country_keep <- function(opt, dat) {
  if (
    !isTRUE(opt$grouping$enable) || is.null(dat) || !"country" %in% names(dat)
  ) {
    return(character(0))
  }
  if (dplyr::n_distinct(dat$country, na.rm = TRUE) == 1L) {
    "country"
  } else {
    character(0)
  }
}

#' .flag_single_country -- `single_country = TRUE` on a parts/terms list when
#' .single_country_keep() applies; the list is unchanged otherwise (no new
#' element, so multi-country outputs are byte-identical).
#' @noRd
.flag_single_country <- function(parts, opt, dat) {
  if (length(.single_country_keep(opt, dat))) {
    parts$single_country <- TRUE
  }
  parts
}


# ============================================================================
# Outbreak-size model helpers
# ============================================================================

validate_size_options <- function(opt) {
  # Advisory only: Convention 2 means build_size_formula() already pulls in
  # the fixed effect for any active slope flag, so a mismatch here does not
  # break the model.  The warning is kept as a discoverability aid; the
  # occurrence side has no equivalent (build_formula_terms() handles it
  # silently via resolve_slopes()).
  stopifnot(!is.null(opt$size))

  valid_phi_re <- c(
    "constant",
    "region",
    "region_country",
    "region_country_adm1"
  )
  phi_re <- opt$size$phi_re %||% "constant"
  if (!phi_re %in% valid_phi_re) {
    stop(sprintf(
      "size.phi_re must be one of: %s (got '%s')",
      paste(valid_phi_re, collapse = ", "),
      phi_re
    ))
  }

  slopes <- resolve_slopes(opt$size$random_effects, opt$size$admin_level)
  size_covs <- opt$size$covariates %||% list()
  enso_on <- isTRUE(opt$size$covariates$enso)

  for (cov in names(slopes)) {
    if (grepl("^enso_lag\\d+$", cov)) {
      if (!enso_on) {
        warning(
          sprintf(
            "size.random_effects slope on %s requested but size.covariates$enso is FALSE. ",
            cov
          ),
          "Set size.covariates.enso: true to enable."
        )
      }
      next
    }
    is_cov_on <- isTRUE(size_covs[[cov]])
    if (!is_cov_on) {
      warning(sprintf(
        "size.random_effects slope on %s requested but size.covariates$%s is FALSE.",
        cov,
        cov
      ))
    }
  }
  invisible(opt)
}

#' resolve_slopes
#'
#' Canonicalise a random_effects block into a per-covariate, per-level
#' specification of slope random effects. Works for both occurrence and size.
#' Returns a list keyed by covariate name, each value a list of
#' list(region, country, adm1, adm2) booleans.
#'
#' @param re_block    the random_effects list (e.g. opt$size$random_effects)
#' @param admin_level "ADM1" | "ADM2"
resolve_slopes <- function(re_block, admin_level = "ADM1") {
  re <- re_block %||% list()
  out <- list()

  # Guard against typo'd flags. Slope keys are decoded by regex below, so a
  # near-miss such as `re_adm1_mean_precip_sloep` matches nothing and is
  # dropped without a trace -- taking with it the implied fixed effect under
  # Convention 2. The model then fits without the covariate at all.
  nms <- names(re) %||% character(0)
  known_re_keys <- c(
    "re_region_country_adm1_intercept",
    "re_adm2_intercept",
    "re_adm2_leaf_scope",
    "slopes"
  )
  unknown <- setdiff(
    nms[
      !grepl("^re_adm1_.*_slope$", nms) &
        !grepl("^re_region_country_.*_slope$", nms)
    ],
    known_re_keys
  )
  if (length(unknown) > 0L) {
    warning(
      "resolve_slopes(): ignoring unrecognised random_effects key(s): ",
      paste(unknown, collapse = ", "),
      ". Expected re_adm1_<cov>_slope, re_region_country_<cov>_slope, or one of: ",
      paste(known_re_keys, collapse = ", "),
      ". A typo here silently drops both the random slope and the ",
      "fixed effect it implies.",
      call. = FALSE
    )
  }

  legacy_adm1_levels <- c("region", "country", "adm1")
  if (identical(admin_level, "ADM2")) {
    legacy_adm1_levels <- c(legacy_adm1_levels, "adm2")
  }

  norm_cov <- function(c) {
    if (identical(c, "ttlo")) "time_to_last_outbreak" else c
  }

  for (nm in names(re)) {
    if (!isTRUE(re[[nm]])) {
      next
    }
    if (grepl("^re_adm1_.*_slope$", nm)) {
      cov <- norm_cov(sub("^re_adm1_(.*)_slope$", "\\1", nm))
      out[[cov]] <- .lvl_set(legacy_adm1_levels)
    } else if (grepl("^re_region_country_.*_slope$", nm)) {
      cov <- norm_cov(sub("^re_region_country_(.*)_slope$", "\\1", nm))
      # Skip non-covariate slopes (intercept, basehaz, ongoing, ttlo, adm1_intercept)
      if (cov %in% c("adm1_intercept", "basehaz")) {
        next
      }
      out[[cov]] <- .lvl_set(c("region", "country"))
    }
  }

  # Explicit per-level map — overrides at level granularity
  if (!is.null(re$slopes) && is.list(re$slopes)) {
    for (cov in names(re$slopes)) {
      spec <- re$slopes[[cov]] %||% list()
      cur <- out[[cov]] %||% .lvl_set(character(0))
      for (lvl in c("region", "country", "adm1", "adm2")) {
        if (!is.null(spec[[lvl]])) cur[[lvl]] <- isTRUE(spec[[lvl]])
      }
      out[[cov]] <- cur
    }
  }

  # Drop entries with no active level
  keep <- vapply(
    out,
    function(s) any(vapply(s, isTRUE, logical(1))),
    logical(1)
  )
  out[keep]
}

#' .lvl_set — helper: list(region, country, adm1, adm2) all FALSE except the
#' levels named in `active`.
.lvl_set <- function(active) {
  out <- list(region = FALSE, country = FALSE, adm1 = FALSE, adm2 = FALSE)
  for (lvl in active) {
    out[[lvl]] <- TRUE
  }
  out
}

#' .seasonal_terms
#'
#' Resolve the calendar-month seasonal fixed-effect term(s), shared by the
#' size and monthly formula builders.
#'
#' @param term  "none" | "factor_month" | "cyclic" | "harmonic2" (NULL treated as "none")
#' @param label appended to the "Unknown seasonal_term" warning for caller context
#' @return character vector of term strings (character(0) for "none")
.seasonal_terms <- function(term, label = "") {
  switch(
    term %||% "none",
    "none" = character(0),
    "factor_month" = "factor(month)",
    "cyclic" = "s(month, bs = 'cc', k = 6)",
    "harmonic2" = c("sin_month", "cos_month", "sin2_month", "cos2_month"),
    {
      warning("Unknown seasonal_term", label, ": ", term)
      character(0)
    }
  )
}

#' .phi_re_term
#'
#' Resolve the negbinomial dispersion (phi/shape) random-effect term shared
#' by the size and monthly formula builders, under either the
#' AFRO_region-rooted hierarchy or the country-rooted (grouping-mode)
#' hierarchy, dropping levels sparse in `dat`.
#'
#' @param phi_re      "constant" | "region" | "region_country" | "region_country_adm1"
#' @param drop_region logical — TRUE under grouping mode (re-roots at country)
#' @param dat         data frame to check RE-group sparsity against
#' @return character(1) RE term string (e.g. "(1 | AFRO_region/country)"), or
#'         character(0) if phi_re == "constant" or all candidate levels are
#'         sparse (<2 distinct values) in `dat`.
.phi_re_term <- function(phi_re, drop_region, dat, keep = character(0)) {
  if (identical(phi_re, "constant") || is.null(phi_re)) {
    return(character(0))
  }

  phi_rhs_str <- if (isTRUE(drop_region)) {
    switch(
      phi_re,
      region = "(1 | country)",
      region_country = "(1 | country)",
      region_country_adm1 = "(1 | country/ADM1)"
    )
  } else {
    switch(
      phi_re,
      region = "(1 | AFRO_region)",
      region_country = "(1 | AFRO_region/country)",
      region_country_adm1 = "(1 | AFRO_region/country/ADM1)"
    )
  }

  phi_rhs_vec <- filter_re_terms_for_data(c(phi_rhs_str), dat, keep = keep)
  if (length(phi_rhs_vec) == 0L) {
    message(
      "[phi_re] All phi grouping levels sparse (<2 distinct values); using constant phi."
    )
    character(0)
  } else {
    phi_rhs_vec[1L]
  }
}

#' .pop_elasticity_terms
#'
#' Frees the population elasticity (size and monthly parts builders).
#' `offset(log(pop))` imposes an elasticity of exactly 1. The offset is kept and
#' a deviation term `log_pop_c` added beside it, so the coefficient reads as
#' `beta - 1` and `beta = 1` reproduces the offset-only model. `log_pop_c` must
#' be centred within admin_level (.add_log_pop_centred(); monthly:
#' load_monthly_data()), or that identity breaks. `free_by_level` is
#' implemented but not the reported model.
#'
#' @param cfg      the track's config block (opt$size | opt$monthly)
#' @param blk_name "size" | "monthly" -- used only in the error message
#' @return character fixed-effect terms to append (empty when "fixed")
#' @noRd
.pop_elasticity_terms <- function(cfg, blk_name) {
  pop_elasticity <- cfg$pop_elasticity %||% "fixed"
  if (identical(pop_elasticity, "fixed")) {
    return(character())
  }
  if (!isTRUE(cfg$offset_pop)) {
    stop(
      blk_name,
      ".pop_elasticity = '",
      pop_elasticity,
      "' requires ",
      blk_name,
      ".offset_pop = true: the elasticity term is parameterised as a ",
      "DEVIATION from the offset, so without the offset it estimates ",
      "beta - 1 while the model carries no beta at all.",
      call. = FALSE
    )
  }
  switch(
    pop_elasticity,
    free = "log_pop_c",
    free_by_level = "log_pop_c:admin_level"
  )
}

#' .build_count_track_parts
#'
#' Shared body of build_size_formula_parts() / build_monthly_formula_parts(). The two
#' negative-binomial tracks differ only in:
#'   - the config block and the anchor/offset tag ("size"/"siz", "monthly"/"mon");
#'   - the RE leaf column: size's ADM2 vs monthly's gid (the GAM's
#'     s(gid, bs = "re"));
#'   - the label used in error messages (" (size model)" / " (monthly model)");
#'   - monthly's `extra_check`: cyclic s(month, bs = 'cc') cannot carry RE
#'     slopes, checked after the seasonal term and before the offsets.
#' Proven byte-identical by tests/test_formula_builder_golden.R (every real
#' config + OFAT sweeps + sparse frames, both tracks).
#'
#' @noRd
.build_count_track_parts <- function(
  opt,
  dat,
  block,
  tag,
  leaf_col,
  label,
  extra_check = NULL
) {
  cfg <- opt[[block, exact = FALSE]]
  # The track's own runtime default (size brms, monthly brms_cb), as the fit
  # scripts resolve it; resolved configs always carry the key.
  be <- cfg[["backend"]] %||% get_default_options()[[block]][["backend"]]
  re_sep <- if (isTRUE(cfg[[be]]$re_uncorrelated)) "||" else "|"

  rfc <- .resolve_fixed_covariates(cfg, cfg$admin_level, opt$enso$lag_months)
  slopes <- rfc$slopes
  slope_covs <- rfc$slope_covs
  fixed_terms <- c("1", rfc$fixed_covs)

  # Per-level RE group strings (shared resolver). Country-grouping mode drops
  # the AFRO_region level (region entry -> NULL).
  drop_region <- isTRUE(opt$grouping$enable)
  group_by_level <- re_group_by_level(leaf_col, drop_region = drop_region)
  re_terms <- character()

  re_terms <- c(
    re_terms,
    .scoped_intercept_re_terms(
      cfg,
      cfg$admin_level,
      group_by_level,
      drop_region,
      label
    )
  )

  re_terms <- c(
    re_terms,
    .slope_re_terms(
      slope_covs,
      slopes,
      group_by_level,
      re_sep,
      c("region", "country", "adm1", "adm2")
    )
  )

  # Seasonal term
  seas_terms <- .seasonal_terms(cfg$seasonal_term, label)

  if (!is.null(extra_check)) {
    extra_check(cfg, slope_covs)
  }

  # Offsets -- shared with the other three tracks (.offset_terms()).
  off_terms <- .offset_terms(opt, tag)

  # P4: free the population elasticity (see .pop_elasticity_terms()).
  fixed_terms <- c(fixed_terms, .pop_elasticity_terms(cfg, block))

  # Drop fixed-effect factor/character columns with <2 distinct levels in the
  # data (e.g. single-TTLO-bin enso_south subset), and cascade removal to any
  # RE slope terms referencing the dropped column. Then drop or trim RE terms
  # whose grouping columns have < 2 distinct values in the data (e.g. a
  # country-grouped fit with outbreaks from only one country): slash-nested
  # intercept REs are trimmed to the deepest valid level; single-level slope
  # RE groups are dropped entirely.
  fft <- filter_fixed_factor_terms_for_data(fixed_terms, re_terms, dat)
  fixed_terms <- fft$fixed_terms
  # Single-country keep: monthly only (size is retired downstream and keeps
  # the plain trimming).
  keep <- if (identical(block, "monthly")) {
    .single_country_keep(opt, dat)
  } else {
    character(0)
  }
  re_terms <- filter_re_terms_for_data(fft$re_terms, dat, keep = keep)

  out <- list(
    fixed_terms = fixed_terms,
    seas_terms = seas_terms,
    re_terms = re_terms,
    off_terms = off_terms,
    slope_covs = slope_covs,
    phi_re = cfg$phi_re %||% "constant"
  )
  if (length(keep)) {
    out$single_country <- TRUE
  }
  out
}

#' .monthly_cyclic_slope_check
#'
#' Cyclic s(month, bs='cc') cannot carry RE slopes -- plain harmonic columns
#' can. A RE slope on a seasonal column with cyclic selected would imply a
#' fixed effect on a column .seasonal_terms() never builds for "cyclic",
#' which would otherwise fail downstream as an opaque missing-column error.
#' @noRd
.monthly_cyclic_slope_check <- function(cfg, slope_covs) {
  if (identical(cfg$seasonal_term, "cyclic")) {
    seasonal_cols <- c("sin_month", "cos_month", "sin2_month", "cos2_month")
    bad <- intersect(slope_covs, seasonal_cols)
    if (length(bad) > 0L) {
      stop(
        "monthly.seasonal_term = 'cyclic' cannot carry a random slope on ",
        "seasonal column(s): ",
        paste(bad, collapse = ", "),
        ". Use seasonal_term = 'harmonic2' for a seasonal RE slope, or ",
        "drop the re_adm1_<col>_slope flag."
      )
    }
  }
  invisible(NULL)
}

#' .build_nb_formula
#'
#' Shared body of build_size_formula() / build_monthly_formula():
#' assemble the brms negative-binomial formula from a parts list, plus the
#' optional overdispersion (`shape ~ 1 + ...`) RE chain. Country-grouping mode
#' drops the (meaningless) AFRO_region root from that chain, re-rooting at
#' country.
#' @noRd
.build_nb_formula <- function(opt, dat, parts, response) {
  phi_re <- parts$phi_re

  rhs <- paste(
    c(parts$fixed_terms, parts$seas_terms, parts$re_terms, parts$off_terms),
    collapse = " + "
  )
  bf_main <- brms::bf(stats::as.formula(paste(response, "~", rhs)))

  drop_region <- isTRUE(opt$grouping$enable)
  phi_term <- .phi_re_term(phi_re, drop_region, dat)
  if (length(phi_term) == 0L) {
    bf_main
  } else {
    bf_main + brms::lf(stats::as.formula(paste("shape ~ 1 +", phi_term)))
  }
}


#' build_size_formula_parts
#'
#' Extract the formula components for the size model without assembling the
#' brms::bf() object.  Intended as a building block for build_size_formula()
#' and for 04d_run_siz_brms_cb.R, which needs to modify the terms before
#' assembling the final formula.
#'
#' @param opt  parsed options list (from make_options_from_config)
#' @param dat  size dataset (currently unused; kept for interface parity)
#' @return named list:
#'   $fixed_terms  character — "1" + FE covariate names
#'   $seas_terms   character — seasonal term(s) or character(0)
#'   $re_terms     character — "(x - 1 | g)" RE strings
#'   $off_terms    character — "offset(log(...))" strings
#'   $slope_covs   character — covariate names with active slope REs
#'   $phi_re       character(1) — "constant" | "region" | ...
build_size_formula_parts <- function(opt, dat) {
  stopifnot(!is.null(opt$size))
  # Size leaf col = ADM2. Shared body: .build_count_track_parts().
  .build_count_track_parts(
    opt,
    dat,
    "size",
    "siz",
    leaf_col = "ADM2",
    label = " (size model)"
  )
}


#' build_size_formula
#'
#' Construct the brms::bf() object for the size model from `opt`.
#' Convention 2: fixed-effect covariates = union of (a) size.covariates
#' toggles set TRUE and (b) slope RE implied covariates.
#' ENSO fixed effect when size.covariates.enso = TRUE.
#' Priors from opt$size$brms$prior_*.
#'
#' Delegates term construction to build_size_formula_parts(); callers that
#' need to modify terms before assembly (e.g. 04d_run_siz_brms_cb.R) should
#' call build_size_formula_parts() directly.
build_size_formula <- function(opt, dat) {
  .build_nb_formula(
    opt,
    dat,
    build_size_formula_parts(opt, dat),
    "suspected_cases"
  )
}


# ============================================================================
# Outbreak-duration model helpers (cessation hazard)
# ============================================================================

# Column-name prefixes for the two optional M-spline bases. The bases themselves
# are built in 04_fit_dur_brms.R (splines2::mSpline) and bound onto the data
# before the formula is assembled, exactly as 04_fit_occ_brms.R does for
# the occurrence baseline. The builder below only NAMES them, so it stays a pure
# function of (opt, dat) and is testable without splines2.
.DUR_MS_PREFIX <- "msoms" # basehaz_duration = "ms", basis on mso
.DUR_SEAS_MS_PREFIX <- "seasms" # basehaz_poisson  = "ms", basis on t_start_seas

.dur_ms_cols <- function(dat, prefix, what, cfg_key) {
  cols <- grep(sprintf("^%s[0-9]+$", prefix), names(dat), value = TRUE)
  if (length(cols) == 0L) {
    stop(
      sprintf(
        paste(
          "duration.%s = 'ms' but no %s* columns are present in the",
          "data. The %s M-spline basis must be bound onto the frame",
          "before build_duration_formula*() is called (see",
          "04_fit_dur_brms.R)."
        ),
        cfg_key,
        prefix,
        what
      ),
      call. = FALSE
    )
  }
  sort(cols)
}

#' build_duration_formula_parts
#'
#' Term construction for the discrete-time cessation hazard. Templated on
#' build_formula_terms() (occurrence) for the covariate/slope conventions, but
#' following build_size_formula_parts() in ending with the two data-aware
#' filters — filter_fixed_factor_terms_for_data() and filter_re_terms_for_data()
#' — which the occurrence builder omits. Those matter more here than on either
#' existing track: the duration panel is one row per active month of 1582
#' episodes, so a per-country or per-group subset can easily leave a factor
#' covariate, an RE grouping column, or an mso band with a single level.
#'
#' Two baseline hazards on two different clocks, independently switchable:
#'   $bh_terms    duration.basehaz_duration on mso        (primary)
#'   $seas_terms  duration.basehaz_poisson  on t_start_seas (default off),
#'                then duration.seasonal_term's calendar-month harmonics
#'                (sin_month .. cos2_month; default off)
#'
#' The offset is TIME, not population: log_offset = log(pmax(exposure, 1)) with
#' `exposure` the days of the calendar month intersected with [TL, TR]
#' (analysis/utils/dur_data.R). It is what makes the Poisson fit a per-day
#' hazard rather than a per-row count, so it is never optional.
#'
#' @param opt resolved options list (must have opt$duration)
#' @param dat duration panel from load_duration_data()
#' @return named list(fixed_terms, bh_terms, seas_terms, re_terms, off_terms,
#'   slope_covs)
build_duration_formula_parts <- function(opt, dat) {
  stopifnot(!is.null(opt$duration))
  dur <- opt$duration
  be <- dur$backend %||% "brms"
  re_sep <- if (isTRUE(dur[[be]]$re_uncorrelated)) "||" else "|"

  # --- Covariates: same union rule as the size track (convention 2) ---------
  rfc <- .resolve_fixed_covariates(
    dur,
    dur$admin_level %||% "ADM1",
    opt$enso$lag_months
  )
  slopes <- rfc$slopes
  slope_covs <- rfc$slope_covs
  fixed_terms <- c("1", rfc$fixed_covs)

  # --- Primary baseline hazard: months since onset --------------------------
  # `mso_band` is rebuilt from duration.mso_breaks at load time
  # (load_duration_data), so a band sensitivity run needs no data re-prep.
  bh_terms <- switch(
    dur$basehaz_duration %||% "factor",
    none = character(0),
    factor = "mso_band",
    ms = .dur_ms_cols(dat, .DUR_MS_PREFIX, "mso", "basehaz_duration"),
    {
      warning("Unknown duration.basehaz_duration: ", dur$basehaz_duration)
      character(0)
    }
  )
  # A subset with a single occupied band carries no baseline contrast; keeping
  # the term would make the design matrix rank-deficient against the intercept.
  if (
    identical(bh_terms, "mso_band") &&
      length(unique(stats::na.omit(dat[["mso_band"]]))) < 2L
  ) {
    message(
      "[duration] mso_band has <2 occupied levels in this subset; ",
      "dropping the mso baseline (constant hazard for this fit)."
    )
    bh_terms <- character(0)
  }

  # --- Secondary baseline hazard: seasonal time -----------------------------
  seas_terms <- switch(
    dur$basehaz_poisson %||% "none",
    none = character(0),
    factor = "factor(month_seas)",
    cyclic = ,
    seasonal = "s(t_start_seas, bs = 'cc', k = 8)",
    ms = .dur_ms_cols(dat, .DUR_SEAS_MS_PREFIX, "seasonal", "basehaz_poisson"),
    {
      warning("Unknown duration.basehaz_poisson: ", dur$basehaz_poisson)
      character(0)
    }
  )
  # Calendar-month seasonality (duration.seasonal_term), the monthly track's
  # harmonic2 columns; load_duration_data() builds them. Only none | harmonic2:
  # the simulator kernels and 09's seasonality step score the four harmonics
  # and nothing else.
  dur_seas <- dur[["seasonal_term"]] %||% "none"
  if (!dur_seas %in% c("none", "harmonic2")) {
    stop(
      "duration.seasonal_term must be 'none' or 'harmonic2'; got '",
      dur_seas,
      "'.",
      call. = FALSE
    )
  }
  seas_terms <- c(seas_terms, .seasonal_terms(dur_seas, " (duration)"))

  # --- Random effects -------------------------------------------------------
  drop_region <- isTRUE(opt$grouping$enable)
  group_by_level <- re_group_by_level("gid", drop_region = drop_region)
  re_terms <- character()

  # `admin_level: ADM1-ADM2` + `re_adm2_leaf_scope: adm2_only` adds the
  # monthly track's `(0 + is_adm2 | country:ADM1:gid_leaf)` leaf.
  re_terms <- c(
    re_terms,
    .scoped_intercept_re_terms(
      dur,
      dur$admin_level %||% "ADM1",
      group_by_level,
      drop_region
    )
  )

  re_terms <- c(
    re_terms,
    .slope_re_terms(
      slope_covs,
      slopes,
      group_by_level,
      re_sep,
      c("region", "country", "adm1")
    )
  )

  # --- Offset ---------------------------------------------------------------
  # Shared with the other three tracks (.offset_terms()). With the anchor gate
  # off this is exactly the "offset(log_offset)" string it hardcoded before.
  off_terms <- .offset_terms(opt, "dur")

  # Data-aware pruning (see the roxygen note above).
  fft <- filter_fixed_factor_terms_for_data(fixed_terms, re_terms, dat)
  fixed_terms <- fft$fixed_terms
  re_terms <- filter_re_terms_for_data(
    fft$re_terms,
    dat,
    keep = .single_country_keep(opt, dat)
  )

  .flag_single_country(
    list(
      fixed_terms = fixed_terms,
      bh_terms = bh_terms,
      seas_terms = seas_terms,
      re_terms = re_terms,
      off_terms = off_terms,
      slope_covs = slope_covs
    ),
    opt,
    dat
  )
}

#' build_duration_formula
#'
#' Assemble the brms::bf() object for the cessation hazard. Response is
#' `y = as.integer(ended)` (set by load_duration_data), family Poisson with the
#' time offset — the Holford / Laird-Olivier equivalent of a discrete-time
#' hazard, which is why the occurrence backend's offset apparatus transfers
#' unchanged (04_fit_occ_brms.R).
#'
#' Callers that need to modify terms before assembly should call
#' build_duration_formula_parts() directly, as 04_fit_siz_brms_cb.R does on the
#' size side.
build_duration_formula <- function(opt, dat) {
  parts <- build_duration_formula_parts(opt, dat)
  rhs <- paste(
    c(
      parts$fixed_terms,
      parts$bh_terms,
      parts$seas_terms,
      parts$re_terms,
      parts$off_terms
    ),
    collapse = " + "
  )
  brms::bf(stats::as.formula(paste("y ~", rhs)))
}


# ============================================================================
# Outbreak monthly-size model helpers (per-calendar-month size, response
# cases_month) — fourth model track alongside occurrence/size/duration.
# Bayesian hierarchical version of a validated exploratory GAM.
# Reuses the same private helpers the size and duration builders already
# share (see the docstring above .resolve_fixed_covariates()); this is the
# third caller.
# ============================================================================

#' build_monthly_formula_parts
#'
#' Term construction for the monthly outbreak-size model. Templated on
#' build_size_formula_parts(), with two genuinely new pieces: the offset is
#' pop (offset_pop toggle, size's key) x log_offset (unconditional, duration's
#' person-time offset -- already on the panel as
#' log_offset = log(pmax(exposure, 1))), and the ADM2-leaf RE column is `gid`
#' (re_group_by_level("gid", ...)) rather than size's `"ADM2"`, matching the
#' GAM's s(gid, bs = "re").
#'
#' @param opt resolved options list (must have opt$monthly)
#' @param dat monthly panel from load_monthly_data()
#' @return named list(fixed_terms, seas_terms, re_terms, off_terms,
#'   slope_covs, phi_re)
build_monthly_formula_parts <- function(opt, dat) {
  stopifnot(!is.null(opt$monthly))
  # Leaf col = gid (the GAM's s(gid, bs="re")), unlike size's ADM2; offset =
  # pop (toggle) x person-time (unconditional, duration's behaviour). Shared
  # body: .build_count_track_parts(), plus the cyclic-slope guard.
  .build_count_track_parts(
    opt,
    dat,
    "monthly",
    "mon",
    leaf_col = "gid",
    label = " (monthly model)",
    extra_check = .monthly_cyclic_slope_check
  )
}

#' build_monthly_formula
#'
#' Construct the brms::bf() object for the monthly outbreak-size model from
#' `opt`. Response is cases_month (per-calendar-month outbreak size).
#' Delegates term construction to build_monthly_formula_parts(); callers that
#' need to modify terms before assembly (e.g. 04_fit_mon_brms_cb.R, which
#' appends DLNM crossbasis columns) should call build_monthly_formula_parts()
#' directly, as 04_fit_siz_brms_cb.R does on the size side.
build_monthly_formula <- function(opt, dat) {
  .build_nb_formula(
    opt,
    dat,
    build_monthly_formula_parts(opt, dat),
    "cases_month"
  )
}
