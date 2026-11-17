# priors.R — brms prior construction for every track
#
# Provides:
#   build_occ_priors()        occurrence
#   build_duration_priors()   duration
#   build_monthly_priors()    monthly case counts
#   build_size_priors()       size (retired downstream, kept)
#   .PRIOR_DEFAULTS           default prior strings per track (also read by
#                             make_batch_configs.R)
#
# Details:
#   Callers source base.R (`%||%`) and formula.R (build_*_formula_parts() when
#   `parts` is NULL, resolve_slopes(), .anchor_gate_on(), .occ_basehaz()).
#   tests/test_prior_builder_golden.R pins the output byte for byte.

# ============================================================================
# Shared prior-builder helpers (brms prior consolidation)
# ============================================================================
#
# Internals shared by every track's prior builder. Same pattern as the
# formula builders: only internals are shared, each track keeps its own public
# builder, its own gates, and its own `c()` assembly order.
# analysis/tests/test_prior_builder_golden.R pins the output byte for byte.

#' .PRIOR_DEFAULTS
#'
#' One table of default brms prior strings per track. `sd` is split into two
#' group-scoped tiers, `sd_region_country` and `sd_adm1` (see
#' .re_group_tiers()/.tiered_sd_priors() below) -- ADM2 counts as the ADM1
#' tier. Both tiers default to the same value (2.5 for occurrence/duration,
#' 0.5 for size/monthly). A config can differentiate the two tiers via
#' `prior_sd_region_country` / `prior_sd_adm1`; the legacy `prior_sd` key
#' still works and seeds both tiers uniformly (see .prior_strings()). The 5x
#' gap between occurrence/duration and size/monthly, and the `lkj(2)` vs
#' no-`cor`-prior split, remain deliberate, unexplained asymmetries on
#' structurally identical crossbasis random slopes (nested_slope_re_terms()).
#' Do not change these values as part of a refactor commit.
#'
#' `intercept_anchored` is the intercept prior used when the track's empirical
#' baseline-rate anchor is switched on (`.anchor_gate_on()`). It is deliberately
#' the SAME for all four tracks, because with the anchor on the intercept no
#' longer has to span that track's offset scale -- the anchor's denominator is
#' exactly the product of the track's other offsets, so
#' `mu = <other offsets> x anchor x exp(lp)` and `exp(Intercept) ~ 1` is the
#' correct prior centre everywhere. The four DIFFERENT un-anchored values above
#' are four different corrections for four different offset scales; the anchored
#' value is one residual-scale prior. See utils/anchor.R's header for the
#' prior-predictive blowups that motivated it.
#'
#' `sds` is the wiggliness prior on a penalised-spline basis --
#' the occurrence/duration `basehaz_poisson: cyclic|seasonal` term
#' `s(t_start_seas, bs = 'cc', k = 8)` and the size/monthly `seasonal_term:
#' cyclic` term `s(month, bs = 'cc', k = 6)`. It is entered here at EXACTLY
#' brms's own package default, `student_t(3, 0, 2.5)` (verified via
#' `brms::get_prior()` for both the Poisson and negbinomial families, and at
#' the Stan level: adding this prior explicitly produces byte-identical
#' `make_stancode()` output -- `lprior += student_t_lpdf(sds_1 | 3, 0, 2.5)` --
#' to leaving it implicit). So this key is a no-op at its default; it makes
#' the value visible and overridable (`prior_sds:`) instead of inherited from
#' brms.
#'
#' That visibility is the point. In prior predictive checks with the
#' baseline-rate anchor on, the occurrence track overshoots the observed total
#' far more than the other three tracks, and this term is the cause: a
#' `factor` basehaz instead of the spline removes most of the overshoot, while
#' tightening `prior_b` or both sd tiers does not. Unlike `sd` below, `sds` is not a
#' nuisance variance: it
#' is the smoothing penalty on the seasonal baseline hazard, i.e. directly on
#' the project's headline scientific object. Whether 2.5 is defensible is an
#' open question; changing it would alter every deployed occurrence config, so
#' -- exactly like the `sd` values above -- do not change it as part of a
#' refactor commit.
.PRIOR_DEFAULTS <- list(
  occurrence = list(
    intercept = "normal(0, 1.5)",
    intercept_anchored = "normal(0, 1)",
    b = "normal(0, 0.5)",
    sd_region_country = "normal(0, 2.5)",
    sd_adm1 = "normal(0, 2.5)",
    sds = "student_t(3, 0, 2.5)",
    cor = "lkj(2)"
  ),
  duration = list(
    intercept = "normal(0, 1.5)",
    intercept_anchored = "normal(0, 1)",
    b = "normal(0, 0.5)",
    sd_region_country = "normal(0, 2.5)",
    sd_adm1 = "normal(0, 2.5)",
    sds = "student_t(3, 0, 2.5)",
    cor = "lkj(2)"
  ),
  size = list(
    intercept = "normal(-3, 1)",
    intercept_anchored = "normal(0, 1)",
    b = "normal(0, 0.5)",
    sd_region_country = "normal(0, 0.5)",
    sd_adm1 = "normal(0, 0.5)",
    sds = "student_t(3, 0, 2.5)",
    cor = NULL,
    phi_intercept = "normal(1, 1)",
    phi_sd = "normal(0, 0.5)"
  ),
  monthly = list(
    intercept = "normal(-3, 1)",
    intercept_anchored = "normal(0, 1)",
    b = "normal(0, 0.5)",
    sd_region_country = "normal(0, 0.5)",
    sd_adm1 = "normal(0, 0.5)",
    sds = "student_t(3, 0, 2.5)",
    cor = NULL,
    phi_intercept = "normal(1, 1)",
    phi_sd = "normal(0, 0.5)"
  )
)

#' .prior_strings
#'
#' Resolve each of a track's prior keys as
#' `brms_block[[paste0("prior_", k)]] %||% .PRIOR_DEFAULTS[[track]][[k]]` --
#' config can override, but the hardcoded default always exists.
#' `opt$occurrence$brms$prior_sd` would be honoured if someone added it, but
#' config_defaults.yml/get_default_options() do not.
#'
#' The `sd_region_country` / `sd_adm1` tier keys additionally fall back to the
#' legacy blanket `prior_sd` key before the hardcoded default, so a config that sets
#' `prior_sd:` keeps working unchanged and seeds both tiers identically.
#'
#' Superseded legacy key: when BOTH tier keys are set,
#' `prior_sd` is never read. The local .prior_sd_superseded() then warns
#' (class "prior_sd_superseded"; a warning, not an error: the tiers are the
#' explicit choice) if `prior_sd` is not the track's default. Same
#' value-as-signal rule as .intercept_prior_string() below: a merged
#' `brms_block` cannot tell "written in the YAML" from "came down from the
#' defaults", and the legacy key does come down from them (config_defaults.yml
#' / get_default_options() set `prior_sd` for size, duration and monthly, and
#' the batch generator writes it into the configs beside the tier keys), so
#' a `prior_sd` equal to the default carries no information and is superseded
#' silently. The default is the blanket value, which both tier defaults in
#' .PRIOR_DEFAULTS restate; compared with whitespace removed. Local rather than
#' top-level so the module's symbol table (tests/test_config_split_symbols.R)
#' is unchanged.
#'
#' @param brms_block opt$<track>$brms (or the brms_cb %||% brms merge already
#'   computed by a `_cb` fit script's call site)
#' @param track one of "occurrence" | "duration" | "size" | "monthly"
#' @return named list, same keys as .PRIOR_DEFAULTS[[track]]
.prior_strings <- function(brms_block, track) {
  .prior_sd_superseded <- function(b, track) {
    legacy <- b[["prior_sd"]]
    if (
      is.null(legacy) ||
        is.null(b[["prior_sd_region_country"]]) ||
        is.null(b[["prior_sd_adm1"]])
    ) {
      return(invisible(FALSE))
    }
    default <- .PRIOR_DEFAULTS[[track]][["sd_region_country"]]
    if (identical(gsub("\\s+", "", legacy), gsub("\\s+", "", default))) {
      return(invisible(FALSE))
    }
    msg <- paste0(
      track,
      ".brms.prior_sd = '",
      legacy,
      "' is ignored: both tier ",
      "keys are set (prior_sd_region_country = '",
      b[["prior_sd_region_country"]],
      "', prior_sd_adm1 = '",
      b[["prior_sd_adm1"]],
      "') and win. Drop `prior_sd:`, or set it ",
      "to the track default '",
      default,
      "'."
    )
    warning(structure(
      class = c("prior_sd_superseded", "warning", "condition"),
      list(message = msg, call = NULL)
    ))
    invisible(TRUE)
  }

  defaults <- .PRIOR_DEFAULTS[[track]]
  brms_block <- brms_block %||% list()
  legacy_sd <- brms_block[["prior_sd"]] # back-compat: pre-split blanket key seeds
  # both tiers when a config hasn't split them
  .prior_sd_superseded(brms_block, track)
  stats::setNames(
    lapply(names(defaults), function(k) {
      val <- brms_block[[paste0("prior_", k)]]
      if (is.null(val) && k %in% c("sd_region_country", "sd_adm1")) {
        val <- legacy_sd
      }
      val %||% defaults[[k]]
    }),
    names(defaults)
  )
}

#' .intercept_prior_string
#'
#' The `Intercept` prior for a track, honouring the anchor gate.
#'
#' The config-pin trap this exists to close: 322 configs in the repo carry an
#' explicit `prior_intercept:` key, and `config_defaults.yml` /
#' `get_default_options()` set one for size / duration / monthly too. Flipping
#' `offset_baserate: true` on such a config while the pin still resolved would
#' put an un-anchored intercept prior on an anchored formula -- the offset
#' correction applied twice, silently, with a perfectly plausible-looking fit
#' coming out the other end.
#'
#' A merged `brms_block` cannot distinguish "the user wrote this in their YAML"
#' from "this came down from the defaults", so the value itself is the signal:
#'
#'   - pin == the track's own UN-ANCHORED default -> inherited (or a redundant
#'     restatement of it). It carries no information the default does not, so it
#'     is superseded by `intercept_anchored` silently. This is the case for all
#'     322 configs (233 x `normal(-3, 1)`, 89 x `normal(0, 1.5)`), which is
#'     what lets an existing config be anchored by flipping one flag.
#'   - pin != that default -> a deliberate, hand-tuned choice whose intent under
#'     an anchored offset is genuinely unknowable. HARD ERROR rather than guess.
#'     `prior_intercept_anchored:` is the escape hatch and resolves through the
#'     ordinary `.prior_strings()` mechanism.
#'
#' @param s resolved prior strings from .prior_strings()
#' @param brms_block the same brms block that was passed to .prior_strings()
#' @param track "occurrence" | "duration" | "size" | "monthly"
#' @param anchored logical(1) from .anchor_gate_on()
#' @return character(1) prior string for class = "Intercept"
.intercept_prior_string <- function(s, brms_block, track, anchored) {
  if (!isTRUE(anchored)) {
    return(s$intercept)
  }

  pinned <- (brms_block %||% list())[["prior_intercept"]]
  unanchored_default <- .PRIOR_DEFAULTS[[track]]$intercept
  if (!is.null(pinned) && !identical(pinned, unanchored_default)) {
    stop(
      track,
      ".brms.prior_intercept is pinned to '",
      pinned,
      "', but this ",
      "run has the baseline-rate anchor switched on, which puts ",
      "offset(log(<anchor>)) in the formula and moves the intercept to a ",
      "RESIDUAL scale. Applying the pinned prior there would double-count ",
      "the offset correction. Either drop the `prior_intercept:` key (the ",
      "anchored default '",
      .PRIOR_DEFAULTS[[track]]$intercept_anchored,
      "' then applies), or, if the pin was chosen FOR the anchored scale, ",
      "rename it to `prior_intercept_anchored:`.",
      call. = FALSE
    )
  }
  s$intercept_anchored
}

#' .prior_formula_gates
#'
#' The has_b / has_sd / has_cor gate logic copy-pasted across the occurrence
#' and duration inline prior blocks (04_fit_occ_brms.R:158,
#' 04_fit_dur_brms.R:201-204). Occurrence and duration differ in which gates
#' they actually use (see build_occ_priors()/build_duration_priors()), but the
#' computation itself is identical.
#'
#' @param parts list with $fixed_terms, $re_terms, and optionally $bh_terms /
#'   $seas_terms (absent treated as character(0))
#' @param re_uncorrelated brms_cfg$re_uncorrelated (logical or NULL)
#' @return list(has_b, has_sd, has_cor)
.prior_formula_gates <- function(parts, re_uncorrelated) {
  bh_terms <- parts$bh_terms %||% character(0)
  seas_terms <- parts$seas_terms %||% character(0)
  list(
    has_b = length(c(setdiff(parts$fixed_terms, "1"), bh_terms, seas_terms)) >
      0L,
    has_sd = length(parts$re_terms) > 0L,
    has_cor = !isTRUE(re_uncorrelated) && any(grepl("\\+.*\\|", parts$re_terms))
  )
}

#' .re_term_group_names
#'
#' Expand each RE term's grouping string into the brms-internal group
#' name(s) it actually produces. Slash-nested groups ("A/B/C", used only by
#' .intercept_re_term()) expand into brms's own "/" nesting semantics:
#' separate variance components named "A", "A:B", "A:B:C". Colon-compound
#' groups ("A:B", every slope RE built by .slope_re_terms() /
#' nested_slope_re_terms()) and bare single-column groups are already the
#' literal brms group name and pass through unchanged. Operates on
#' already data-filtered re_terms (see filter_re_terms_for_data()), so a
#' slash-nested term that was trimmed to a single remaining level (e.g.
#' "(1 | ADM1)" after a single-country/-region ENSO-grouped subset) is
#' handled identically to a term that was never nested.
#'
#' @param re_terms character vector of "(... | group)" RE strings
#' @return character vector of brms group names (deduplicated)
.re_term_group_names <- function(re_terms) {
  if (length(re_terms) == 0L) {
    return(character(0))
  }
  names_list <- lapply(re_terms, function(term) {
    grp_str <- sub(".*\\|\\|?\\s+", "", term)
    grp_str <- trimws(sub("\\s*\\)\\s*$", "", grp_str))
    if (grepl("/", grp_str, fixed = TRUE)) {
      comps <- strsplit(grp_str, "/", fixed = TRUE)[[1]]
      vapply(
        seq_along(comps),
        function(i) paste(comps[seq_len(i)], collapse = ":"),
        character(1)
      )
    } else {
      grp_str
    }
  })
  unique(unlist(names_list))
}

#' .re_group_tiers
#'
#' Classify the group names implied by `re_terms` into the region/country
#' tier and the ADM1/ADM2 tier. re_group_by_level() (line ~1321) always
#' includes the literal "ADM1" path segment in the adm1 and adm2 group
#' strings, and never in the region/country ones -- true at every
#' admin_level, under drop_region, and after filter_re_terms_for_data()
#' sparsity trimming -- so a substring test on the already-built,
#' already-filtered group name is sufficient and avoids re-threading
#' re_group_by_level()'s output (which may not match what actually
#' survived data-aware filtering) through every builder.
#' set_prior(group = X) for an X that is not an actual model parameter is a
#' hard brms error (verified empirically: ".validate_prior(): priors do not
#' correspond to any model parameter"), so this MUST operate on the final,
#' post-filtering re_terms -- never a static prediction from
#' re_group_by_level() alone.
#'
#' @param re_terms character vector of "(... | group)" RE strings (final,
#'   post data-aware filtering)
#' @return list(region_country = character, adm1 = character) -- brms group
#'   names for each tier; either may be character(0)
.re_group_tiers <- function(re_terms) {
  groups <- .re_term_group_names(re_terms)
  is_adm1 <- grepl("ADM1", groups, fixed = TRUE)
  list(
    region_country = groups[!is_adm1],
    adm1 = groups[is_adm1]
  )
}

#' .tiered_sd_priors
#'
#' Emit group-scoped class="sd" priors from .re_group_tiers() output instead
#' of one unscoped blanket prior, so the region/country tier and the
#' ADM1/ADM2 tier (intercepts, climate/ENSO slopes, and DLNM crossbasis
#' slopes alike -- every RE in this codebase resolves to one of these two
#' tiers) can be given different variance priors. Never emits a prior for
#' an absent tier: an unscoped/absent-group class="sd" prior on a model
#' with zero matching group-level effects is a hard brms error (verified
#' empirically), so a trimmed-away or never-present tier must stay silent.
#'
#' @param s prior strings from .prior_strings(); must include
#'   $sd_region_country, $sd_adm1
#' @param tiers result of .re_group_tiers(re_terms)
#' @param nlpar optional brms non-linear parameter name (e.g. "lin" for the
#'   brms_kl backend's linear sub-predictor); NULL = brms default ("").
#' @return brmsprior, or NULL when both tiers are empty (safe to c() either
#'   way -- c(x, NULL) is a no-op). Built via do.call(c, <list>) rather than
#'   an accumulating `priors <- c(priors, ...)` loop: c(NULL, <brmsprior>)
#'   dispatches on NULL's (missing) class and silently degrades to a plain
#'   list instead of brms:::c.brmsprior() -- confirmed empirically -- which
#'   would then hard-error one level up when the caller tries to c() it
#'   alongside a real brmsprior ("Cannot add 'list' objects to the prior").
.tiered_sd_priors <- function(s, tiers, nlpar = NULL) {
  nl <- nlpar %||% "" # brms::set_prior()'s own default: NULL is a no-op
  entries <- c(
    lapply(tiers$region_country, function(g) {
      brms::set_prior(s$sd_region_country, class = "sd", group = g, nlpar = nl)
    }),
    lapply(tiers$adm1, function(g) {
      brms::set_prior(s$sd_adm1, class = "sd", group = g, nlpar = nl)
    })
  )
  if (length(entries) == 0L) {
    return(NULL)
  }
  do.call(c, entries)
}

#' .has_smooth_terms
#'
#' TRUE iff the track's model formula contains a penalised mgcv smooth, i.e.
#' iff brms will create a class="sds" parameter. Gating is MANDATORY, not
#' defensive: `set_prior(class = "sds")` against a formula with no smooth is
#' the same hard `.validate_prior()` failure as an absent sd tier ("The
#' following priors do not correspond to any model parameter"), verified
#' empirically. 99 occurrence/duration configs in the repo have no smooth.
#'
#' Two resolution paths, because the smooth enters the formula at two
#' different layers:
#'   - occurrence: the baseline-hazard block is assembled by the FIT SCRIPT
#'     (04_fit_occ_brms.R:149-153 / 04_fit_occ_brms_cb.R:201-205), not by
#'     build_formula_terms(), so it is absent from `parts` and there is nothing
#'     to scan. Resolved from the same backend block the script used, through
#'     the same .occ_basehaz() call -- that shared helper is what guarantees
#'     the prior and the formula cannot disagree about whether the spline is
#'     there.
#'   - duration / size / monthly: the smooth IS in `parts` (as $seas_terms,
#'     from .seasonal_terms() or the duration seasonal-basehaz switch), so it
#'     is detected by the term string itself rather than re-resolved from
#'     config -- the sibling of .re_group_tiers()'s "classify the terms you
#'     actually built, never a static prediction" rule, and for the same
#'     reason: any drift between prediction and formula is a hard brms error.
#'
#' @param parts    formula-parts list; $bh_terms / $seas_terms / $fixed_terms
#'   are scanned when present (absent treated as character(0)). Unused for
#'   occurrence.
#' @param brms_cfg the resolved backend block. Unused except for occurrence.
#' @param track    "occurrence" | "duration" | "size" | "monthly"
#' @return logical(1)
.has_smooth_terms <- function(parts, brms_cfg, track) {
  if (identical(track, "occurrence")) {
    return(identical(.occ_basehaz(brms_cfg), "cyclic"))
  }
  terms <- c(
    parts$bh_terms %||% character(0),
    parts$seas_terms %||% character(0),
    parts$fixed_terms %||% character(0)
  )
  any(grepl("^s\\(", terms))
}

#' .prior_config_slopes
#'
#' has_slopes: TRUE iff the model has any population-level (class="b") fixed
#' effect, from any of four sources -- an active slope RE (resolve_slopes()
#' non-empty), an explicit covariate toggle, covariates$enso, or a
#' seasonal_term that itself resolves to plain class="b" columns
#' ("factor_month" / "harmonic2" -- see .seasonal_terms()). The seasonal
#' check is deliberately narrower than "any non-'none' seasonal_term":
#' build_size_formula_parts()/build_monthly_formula_parts() splice
#' `parts$seas_terms` into the model RHS unconditionally (regardless of
#' `covariates`), so a config with covariates: [] but seasonal_term:
#' harmonic2/factor_month has real class="b" parameters in the fit; without
#' has_slopes brms would silently use its own default flat prior on them in a
#' normal (data) fit, and error outright under sample_prior = "only" ("no
#' proper priors ... parameter 'b'"). "cyclic" is excluded:
#' it resolves to an mgcv smooth term (s(month, bs="cc", ...)), a different
#' parameter-class structure (smooth basis + class="sds"): a `bs="cc"` basis has an
#' empty null space, so brms creates no class="b" coefficient for it at all
#' (`get_prior()` on `s(month, bs="cc", k=6)` yields class="b" rows only for
#' the ordinary covariates -- contrast a thin-plate `s(x, k=5)`, which does add
#' a class="b" `sx_1`). A cyclic seasonal_term therefore genuinely does not
#' make has_slopes TRUE. Its class="sds" parameter is handled separately by
#' .has_smooth_terms() below.
#' Shared by build_size_priors() and build_monthly_priors().
#'
#' @param block opt$size or opt$monthly
#' @param admin_level block$admin_level
#' @return logical(1)
.prior_config_slopes <- function(block, admin_level) {
  slopes <- resolve_slopes(block$random_effects, admin_level)
  covs <- block$covariates %||% list()
  length(slopes) > 0 ||
    length(names(Filter(isTRUE, covs))) > 0 ||
    isTRUE(block$covariates$enso) ||
    (block$seasonal_term %||% "none") %in% c("factor_month", "harmonic2")
}

#' .shape_dpar_priors
#'
#' The `dpar = "shape"` Intercept+sd pair for the negbinomial overdispersion
#' sub-model, or NULL when phi_re == "constant" (no shape sub-model at all).
#' Shared by size and monthly, whose phi-prior blocks are identical.
#'
#' @param strings result of .prior_strings(brms_block, track) -- must include
#'   $phi_intercept, $phi_sd
#' @param phi_re "constant" | "region" | "region_country" | "region_country_adm1"
#' @return brmsprior of length 2, or NULL
# Prior on the country-level RE sds of a grouped fit of ONE country
# (formula.R::.single_country_keep(): the country terms are kept with a single
# level instead of being trimmed). One level carries no information on its sd,
# so a tight prior keeps the country deviation near 0 (the intercept and fixed
# effects absorb the country) while the group names downstream code maps stay.
.SINGLE_COUNTRY_SD <- "normal(0, 0.05)"

#' .single_country_strings -- `s` with the region/country tier set to
#' .SINGLE_COUNTRY_SD when `parts` is flagged single_country (grouped fits
#' have no region level, so that tier is exactly the `country` group).
#' @noRd
.single_country_strings <- function(s, parts) {
  if (isTRUE(parts$single_country)) {
    s$sd_region_country <- .SINGLE_COUNTRY_SD
  }
  s
}

.shape_dpar_priors <- function(strings, phi_re) {
  if (identical(phi_re %||% "constant", "constant")) {
    return(NULL)
  }
  c(
    brms::set_prior(strings$phi_intercept, class = "Intercept", dpar = "shape"),
    brms::set_prior(strings$phi_sd, class = "sd", dpar = "shape")
  )
}

#' .build_hazard_track_priors
#'
#' Shared body of build_occ_priors() / build_duration_priors(): the two Poisson hazard
#' tracks build Intercept, [b],
#' [tiered sd], [sds], [cor] identically except that duration gates the `b`
#' prior on `parts` having fixed effects (`gate_b = TRUE`; brms errors on a
#' prior for an empty class, and the constant-hazard null variant has none)
#' while occurrence always emits it (`gate_b = FALSE` -- do not "harmonise":
#' an occurrence config with no fixed effects would change output).
#' Output pinned byte for byte by tests/test_prior_builder_golden.R.
#'
#' @param block config block name ("occurrence" | "duration")
#' @param tag   anchor-gate track tag ("occ" | "dur")
#' @param nlpar optional brms non-linear parameter carrying the linear
#'   predictor (brms_kl: "lin"). NULL leaves the output unchanged.
#'   When set, the intercept becomes `class = "b", coef =
#'   "Intercept"` (nl models have no Intercept class) and b/sd/sds carry the
#'   nlpar; cor priors never take one (brms scopes them by group only).
#' @noRd
.build_hazard_track_priors <- function(
  opt,
  block,
  tag,
  parts,
  brms_cfg,
  gate_b,
  nlpar = NULL
) {
  brms_cfg <- brms_cfg %||% opt[[block, exact = FALSE]]$brms
  s <- .single_country_strings(.prior_strings(brms_cfg, block), parts)
  gates <- .prior_formula_gates(parts, brms_cfg$re_uncorrelated)
  intercept <- .intercept_prior_string(
    s,
    brms_cfg,
    block,
    .anchor_gate_on(opt, tag)
  )
  nl <- nlpar %||% ""

  priors <- if (is.null(nlpar)) {
    brms::set_prior(intercept, class = "Intercept")
  } else {
    brms::set_prior(intercept, class = "b", coef = "Intercept", nlpar = nl)
  }
  if (!gate_b || gates$has_b) {
    priors <- c(priors, brms::set_prior(s$b, class = "b", nlpar = nl))
  }
  priors <- c(
    priors,
    .tiered_sd_priors(s, .re_group_tiers(parts$re_terms), nlpar)
  )
  if (.has_smooth_terms(parts, brms_cfg, block)) {
    priors <- c(priors, brms::set_prior(s$sds, class = "sds", nlpar = nl))
  }
  if (gates$has_cor) {
    priors <- c(priors, brms::set_prior(s$cor, class = "cor"))
  }
  priors
}

#' .build_count_track_priors
#'
#' Shared body of build_size_priors() / build_monthly_priors().
#' The two negative-binomial tracks differ only in the config block, the
#' anchor-gate tag and the parts builder used when `parts` is NULL. `b` is
#' gated on config slopes and PREPENDED, and the
#' shape/phi priors come last. `opt[[block, exact = FALSE]]` keeps the
#' original `opt$size` / `opt$monthly` partial-matching semantics.
#' Proven byte-identical by tests/test_prior_builder_golden.R.
#'
#' @param block    config block name ("size" | "monthly")
#' @param tag      anchor-gate track tag ("siz" | "mon")
#' @param parts_fn build_size_formula_parts | build_monthly_formula_parts
#' @noRd
.build_count_track_priors <- function(
  opt,
  block,
  tag,
  parts_fn,
  brms_cfg,
  parts
) {
  cfg <- opt[[block, exact = FALSE]]
  parts <- parts %||% parts_fn(opt, NULL)
  has_slopes <- .prior_config_slopes(cfg, cfg$admin_level)
  brms_cfg <- brms_cfg %||% cfg$brms
  s <- .single_country_strings(.prior_strings(brms_cfg, block), parts)
  intercept <- .intercept_prior_string(
    s,
    brms_cfg,
    block,
    .anchor_gate_on(opt, tag)
  )

  priors <- c(
    brms::set_prior(intercept, class = "Intercept"),
    .tiered_sd_priors(s, .re_group_tiers(parts$re_terms))
  )
  if (has_slopes) {
    priors <- c(brms::set_prior(s$b, class = "b"), priors)
  }
  if (.has_smooth_terms(parts, brms_cfg, block)) {
    priors <- c(priors, brms::set_prior(s$sds, class = "sds"))
  }
  phi_priors <- .shape_dpar_priors(s, cfg$phi_re)
  if (!is.null(phi_priors)) {
    priors <- c(priors, phi_priors)
  }
  # single country: the shape RE keeps its `country` level (grouped fits root it
  # at country, formula.R::.phi_re_term()); tight prior on that sd as well.
  if (isTRUE(parts$single_country) && !is.null(phi_priors)) {
    priors <- c(
      priors,
      brms::set_prior(
        .SINGLE_COUNTRY_SD,
        class = "sd",
        group = "country",
        dpar = "shape"
      )
    )
  }
  priors
}

#' build_occ_priors
#'
#' Occurrence-track brms priors. Two carve-outs vs. build_duration_priors():
#'   - has_b is forced TRUE, not derived from `parts`. The baseline-hazard
#'     block is always non-empty but is not
#'     represented in `parts`, so the class-"b" gate cannot be derived from
#'     it. `parts$re_terms` is consumed for has_cor AND for the group-scoped
#'     sd priors via .re_group_tiers() -- occurrence always has at least the
#'     region/country/ADM1 intercept RE in every real config; a synthetic
#'     zero-RE `parts` emits zero sd priors.
#'   - Occurrence has no config surface for its own priors: no
#'     `opt$occurrence$brms$prior_*` key exists in config_defaults.yml or
#'     get_default_options(). `.prior_strings()` would honour one if added.
#'
#' @param opt resolved options list (must have opt$occurrence)
#' @param parts formula-parts list from build_formula_terms() -- $re_terms is
#'   read for has_cor and for the tiered sd priors
#' @param brms_cfg optional override; NULL resolves to opt$occurrence$brms,
#'   preserving every existing call site's behaviour
#' @param nlpar optional non-linear parameter name for the linear predictor
#'   (brms_kl backend: "lin"); NULL = unchanged output
#' @return brmsprior: Intercept, b, [sd tiers], [sds], [cor]
build_occ_priors <- function(opt, parts, brms_cfg = NULL, nlpar = NULL) {
  .build_hazard_track_priors(
    opt,
    "occurrence",
    "occ",
    parts,
    brms_cfg,
    gate_b = FALSE,
    nlpar = nlpar
  )
}

#' build_duration_priors
#'
#' Duration-track brms priors, config-driven via opt$duration$brms$prior_*.
#' Uses all three gates from `parts` (unlike build_occ_priors()): has_b exists
#' because brms errors on a prior for an empty class, and the constant-hazard
#' null variant (basehaz_duration: none, no covariates, no slopes) genuinely
#' has no population-level effects beyond the intercept
#' (04_fit_dur_brms.R:194-199).
#'
#' `prior_mso_band` (brms block, default NULL): its own prior on the
#' episode-age band effects, one `class = "b", coef = "mso_band<level>"` row
#' per non-reference band occupied in `dat`. Those rows override `prior_b` for
#' the band coefficients only. Unset, the output is unchanged and `dat` is not
#' read. Motivation: under prior_b normal(0, 0.5) the band contrasts against
#' the sparse first month are shrunk well below the empirical per-day ratios.
#'
#' @param opt resolved options list (must have opt$duration)
#' @param parts formula-parts list from build_duration_formula_parts()
#' @param brms_cfg optional override; NULL resolves to opt$duration$brms
#' @param dat the duration panel the fit uses (its `mso_band` factor); read only
#'   when `prior_mso_band` is set
#' @return brmsprior: Intercept, [b], [mso_band b], [sd tiers], [sds], [cor]
build_duration_priors <- function(opt, parts, brms_cfg = NULL, dat = NULL) {
  priors <- .build_hazard_track_priors(
    opt,
    "duration",
    "dur",
    parts,
    brms_cfg,
    gate_b = TRUE
  )
  brms_cfg <- brms_cfg %||% opt[["duration"]][["brms"]]
  band_prior <- brms_cfg[["prior_mso_band"]]
  if (is.null(band_prior) || !identical(parts$bh_terms, "mso_band")) {
    return(priors)
  }
  if (is.null(dat) || is.null(dat[["mso_band"]])) {
    stop(
      "build_duration_priors: duration prior_mso_band is set but `dat` with an ",
      "`mso_band` column was not passed.",
      call. = FALSE
    )
  }
  lv <- levels(droplevels(as.factor(dat[["mso_band"]])))[-1L]
  # brms coefficient names for factor levels: "4-5" -> "4M5", "9+" -> "9P".
  coefs <- paste0(
    "mso_band",
    gsub("+", "P", gsub("-", "M", lv, fixed = TRUE), fixed = TRUE)
  )
  do.call(
    c,
    c(
      list(priors),
      lapply(coefs, function(k) {
        brms::set_prior(band_prior, class = "b", coef = k)
      })
    )
  )
}

#' build_size_priors
#'
#' Build brms priors from opt$size$brms$prior_* string specs (or `brms_cfg`,
#' when supplied explicitly).
#'
#' @param opt resolved options list (must have opt$size)
#' @param brms_cfg optional override; NULL resolves to opt$size$brms
#'   (the single-arg call shape used by 04_fit_siz_brms.R).
#'   04_fit_siz_brms_cb.R passes `opt$size$brms_cb %||% opt$size$brms`
#'   explicitly so a `brms_cb.prior_*` override is honoured.
#' @param parts optional formula-parts list from build_size_formula_parts()
#'   -- when supplied by a `_cb` call site, `parts`
#'   already has DLNM/ENSO crossbasis slope RE terms appended
#'   (04_fit_siz_brms_cb.R), so those slopes get correctly tiered sd priors
#'   too. NULL (every other call site) recomputes `parts` from `opt` alone;
#'   build_size_formula_parts() tolerates dat = NULL (every internal
#'   `%in% names(dat)` sparsity check treats an unknown column as "keep,
#'   can't check").
build_size_priors <- function(opt, brms_cfg = NULL, parts = NULL) {
  .build_count_track_priors(
    opt,
    "size",
    "siz",
    build_size_formula_parts,
    brms_cfg,
    parts
  )
}

#' build_monthly_priors
#'
#' Structural copy of build_size_priors(): brms priors from
#' opt$monthly$brms$prior_* string specs (or opt$monthly$brms_cb$prior_*
#' when `brms_cfg` is supplied by a brms_cb call site, see
#' 04_fit_mon_brms_cb.R). Internals go through the shared
#' .prior_config_slopes()/.prior_strings()/.shape_dpar_priors() helpers;
#' pinned by test_prior_builder_golden.R.
#'
#' @param parts optional formula-parts list from build_monthly_formula_parts()
#'   -- see build_size_priors()'s `parts` doc for the
#'   identical rationale (CB/ENSO slope RE terms already appended by
#'   04_fit_mon_brms_cb.R when supplied; NULL recomputes from `opt` alone).
build_monthly_priors <- function(opt, brms_cfg = NULL, parts = NULL) {
  .build_count_track_priors(
    opt,
    "monthly",
    "mon",
    build_monthly_formula_parts,
    brms_cfg,
    parts
  )
}
