# test_hardcoded_default_parity.R — HARD: hardcoded `%||%` literals in
#   paths.R, 04_fit_run.R and every module agree with get_default_options().
#
# Usage:
#   Rscript analysis/tests/test_hardcoded_default_parity.R
#
# Details:
#   Guards the "config computed but not threaded to its consumer" bug class:
#   a literal that drifts from get_default_options() makes a writer and a
#   reader name different artefacts. [1]-[5]: the hand-kept SITE_MAP over
#   paths.R / 04_fit_run.R (new-unmapped, stale-mapped and mapped-but-unequal
#   fail). [6]: simulator and smoothing modules; [7]: every other analysis/
#   and utils/ file, both read off the parse tree by expression root.

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat(
  "# hardcoded default parity (paths.R / 04_fit_run.R / simulator modules vs get_default_options())\n"
)

PATHS_R <- "analysis/utils/paths.R"

# Why this exists: test_config_defaults_parity.R ratchets
# `config_defaults.yml` vs `get_default_options()`, the two *declared*
# default sources. It does not cover a third, undeclared source:
# consumer-side literals hardcoded directly in `paths.R` filename builders
# (`DEFAULT_BACKEND` and friends, plus every bare `%||% "literal"` fallback)
# and in `04_fit_run.R`'s dispatcher. When one of these literals silently
# drifts from `get_default_options()`, the artefact a script WRITES (named
# from the canonical default) and the artefact a second script READS (named
# from the stale literal) disagree, and `file.exists()` guards turn that into
# a silent wrong-file read rather than an error. HARD on the precedent of
# `cli_route_completeness`, hard because its allowlist (`KNOWN_UNROUTED`) is
# `character(0)`.
#
# Design: hybrid, registry-led (same shape as
# `test_cli_route_completeness.R`): a regex mechanically *finds* every
# `%||% "literal"` site and `DEFAULT_*` constant in `paths.R`; it cannot *map*
# them (`mon$backend` -> `monthly.backend`, `opt$baseline` -> top-level
# `baseline`, and some sites have no `get_default_options()` counterpart at
# all). So SITE_MAP hand-maintains the mapping — to a `get_default_options()`
# key path for sites that assert equality, or to an explicit "no counterpart"
# reason for sites that structurally cannot have one. New-unmapped fails;
# stale-mapped fails; mapped-but-unequal fails.
#
# --- Extraction --------------------------------------------------------------

#' Every `<expr> %||% "<literal>"` occurrence in a file, skipping comment
#' lines (trimmed line starts with `#`) so roxygen prose that merely
#' *mentions* the pattern (e.g. paths.R's `make_figure_dir()` docstring) is not
#' picked up as a code site. Returns a data.frame(line, expr, literal),
#' deduplicated on (expr, literal) — the same expression/fallback pair
#' appearing at several call sites (e.g. `opt$outdir %||% "generated_data/"`,
#' 8 sites) is one fact to verify, not eight.
.extract_pct_or_sites <- function(path) {
  txt <- readLines(path, warn = FALSE)
  pat <- "([A-Za-z0-9_$.]+(?:\\[\\[[^]]*\\]\\])?)\\s*%\\|\\|%\\s*\"([^\"]*)\""
  out <- list()
  for (i in seq_along(txt)) {
    if (grepl("^\\s*#", txt[i])) {
      next
    }
    m <- gregexpr(pat, txt[i], perl = TRUE)
    parts <- regmatches(txt[i], m)[[1]]
    if (!length(parts)) {
      next
    }
    grp <- regmatches(txt[i], m)
    mm <- attr(m[[1]], "capture.start")
    for (k in seq_along(parts)) {
      cap <- regmatches(parts[k], regexec(pat, parts[k], perl = TRUE))[[1]]
      out[[length(out) + 1L]] <- data.frame(
        line = i,
        expr = cap[2],
        literal = cap[3],
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, out)
  df[!duplicated(df[c("expr", "literal")]), ]
}

#' Every `DEFAULT_<NAME> <- "<literal>"` constant assignment.
.extract_default_constants <- function(path) {
  txt <- readLines(path, warn = FALSE)
  pat <- "^\\s*(DEFAULT_[A-Z_]+)\\s*<-\\s*\"([^\"]*)\""
  out <- list()
  for (i in seq_along(txt)) {
    cap <- regmatches(txt[i], regexec(pat, txt[i], perl = TRUE))[[1]]
    if (length(cap) == 0L) {
      next
    }
    out[[length(out) + 1L]] <- data.frame(
      line = i,
      expr = cap[2],
      literal = cap[3],
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) {
    return(data.frame(
      line = integer(0),
      expr = character(0),
      literal = character(0)
    ))
  }
  do.call(rbind, out)
}

pct_sites <- .extract_pct_or_sites(PATHS_R)
dflt_sites <- .extract_default_constants(PATHS_R)

cat(sprintf(
  "\n  %%||%% \"literal\" sites (deduped)     : %d\n",
  nrow(pct_sites)
))
cat(sprintf("  DEFAULT_* constant assignments    : %d\n", nrow(dflt_sites)))

# --- SITE_MAP ------------------------------------------------------------
#
# One entry per (expr, literal) tuple found above. `key_path` (or
# `key_paths`, plural, when the same site is exercised for more than one
# model block) names the get_default_options() location that must hold the
# identical value; `normalize` (optional) is applied to both sides before
# comparing (`tolower`, for the one site whose caller already normalizes
# case itself). `reason` (no `key_path`) marks a site with no defaults-side
# counterpart at all — allowlisted the way `.base_rate_flag_token`'s `%||% TRUE`
# already is in the design doc, not drift to be repaired here.

get_by_path <- function(d, path) {
  for (p in path) {
    if (is.null(d)) {
      return(NULL)
    }
    d <- d[[p]]
  }
  d
}

SITE_MAP <- list(
  # -- occurrence -------------------------------------------------------
  list(
    expr = "opt$occurrence$base_rate$estimator",
    literal = "pooled",
    reason = "no `occurrence.base_rate` block in get_default_options() -- the occ base-rate anchor table has no config surface of its own (unlike size$mean_rate), same class as `.base_rate_flag_token`'s `%||% TRUE`."
  ),
  list(
    expr = "occ$admin_level",
    literal = "adm2",
    key_path = c("occurrence", "admin_level"),
    normalize = tolower,
    note = "call site wraps in tolower(); default is \"ADM2\", case-insensitive match is the correct comparison."
  ),
  list(
    expr = "brms_sub$init",
    literal = "random",
    key_path = c("occurrence", "brms", "init")
  ),
  list(
    expr = 'lr[["scale"]]',
    literal = "month",
    key_path = c("occurrence", "kernel", "lag_re", "scale"),
    note = ".kl_token(): kernel.lag_re.scale picks the -lre (month) / -lrl (logit) token."
  ),

  # -- shared top-level ---------------------------------------------------
  list(expr = "opt$baseline", literal = "1990_2019", key_path = c("baseline")),
  list(expr = "opt$country", literal = "all", key_path = c("country")),
  list(
    expr = "opt$outdir",
    literal = "generated_data/",
    key_path = c("outdir")
  ),
  list(
    expr = "opt$enso$source",
    literal = "oni",
    key_path = c("enso", "source")
  ),
  # -- cross-track episode key (episode_key_of()) and the size/monthly RE leaf
  #    scope (.leaf_scope_tag()) --------------------------------------------
  list(expr = "ep$key", literal = "lp_season", key_path = c("episodes", "key")),
  list(
    expr = "re$re_adm2_leaf_scope",
    literal = "all",
    key_paths = list(
      c("size", "random_effects", "re_adm2_leaf_scope"),
      c("monthly", "random_effects", "re_adm2_leaf_scope")
    )
  ),

  # -- estimand token (shared helper, called for block="size" and "monthly") --
  list(
    expr = "blk$estimand$denominator",
    literal = "footprint",
    key_paths = list(
      c("size", "estimand", "denominator"),
      c("monthly", "estimand", "denominator")
    )
  ),
  list(
    expr = "blk$pop_elasticity",
    literal = "fixed",
    key_paths = list(
      c("size", "pop_elasticity"),
      c("monthly", "pop_elasticity")
    )
  ),

  # -- dlnm basis token (shared helper .dlnm_basis_token, called for occ/dur/mon) --
  list(
    expr = "av$fun",
    literal = "bs",
    key_paths = list(
      c("occurrence", "dlnm", "defaults", "argvar", "fun"),
      c("duration", "dlnm", "defaults", "argvar", "fun"),
      c("monthly", "dlnm", "defaults", "argvar", "fun")
    )
  ),
  list(
    expr = "av$df",
    literal = "",
    reason = "formatting sentinel, not a modeling default: the argvar df default is 4L (get_default_options()$<track>$dlnm$defaults$argvar$df), but this %||% \"\" only fires for an argvar that omits df entirely, in which case the basis token drops the df suffix (\"ns\" not \"nsNULL\") rather than inventing a number. Same sentinel class as mon$seasonal_term above -- a fallback for a hand-built config, never reached via make_options_from_config()."
  ),

  # -- size -----------------------------------------------------------
  list(
    expr = "opt$size$backend",
    literal = "brms",
    key_path = c("size", "backend")
  ),
  list(
    expr = "opt$size$family",
    literal = "negbinomial",
    key_path = c("size", "family")
  ),
  list(
    expr = "opt$size$admin_level",
    literal = "ADM1",
    key_path = c("size", "admin_level")
  ),
  list(
    expr = "opt$size$phi_re",
    literal = "constant",
    key_path = c("size", "phi_re")
  ),

  # -- duration ---------------------------------------------------------
  list(
    expr = "dur$base_rate$estimator",
    literal = "pooled",
    reason = "no `duration.base_rate` block in get_default_options() -- same as the occurrence base_rate site above."
  ),
  list(
    expr = "dur$backend",
    literal = "brms",
    key_path = c("duration", "backend")
  ),
  list(
    expr = "dur$basehaz_duration",
    literal = "factor",
    key_path = c("duration", "basehaz_duration")
  ),
  list(
    expr = "dur$basehaz_poisson",
    literal = "none",
    key_path = c("duration", "basehaz_poisson")
  ),
  list(
    expr = 'dur[["seasonal_term"]]',
    literal = "none",
    key_path = c("duration", "seasonal_term")
  ),
  list(
    expr = "dur$admin_level",
    literal = "ADM1",
    key_path = c("duration", "admin_level")
  ),

  # -- monthly ------------------------------------------------------------
  list(
    expr = "mon$base_rate$estimator",
    literal = "pooled",
    reason = "no `monthly.base_rate` block in get_default_options() -- same as the occurrence base_rate site above."
  ),
  list(
    expr = "mon$backend",
    literal = "brms_cb",
    key_path = c("monthly", "backend")
  ),
  list(
    expr = "mon$family",
    literal = "negbinomial",
    key_path = c("monthly", "family")
  ),
  list(
    expr = "mon$admin_level",
    literal = "ADM1-ADM2",
    key_path = c("monthly", "admin_level")
  ),
  list(
    expr = "mon$seasonal_term",
    literal = "none",
    reason = "sentinel, not the modeling default: get_default_options()$monthly$seasonal_term is \"harmonic2\" (confirmed on every deployed monthly fit filename, which carries \"seas-harmonic2\"). This %||% only fires for a hand-built opt missing seasonal_term entirely (never true for a config routed through make_options_from_config(), which fills it before make_monthly_run_file() is reached) -- \"none\" there means \"assume no seasonal term\", not \"assume the canonical default\". Pre-existing behavior, not one of the six findings this gate was built to close; flagged for the record, not fixed here (out of this unit's scope)."
  ),
  list(
    expr = "mon$phi_re",
    literal = "constant",
    reason = "sentinel, not the modeling default: get_default_options()$monthly$phi_re is \"region_country_adm1\" (confirmed on every deployed monthly fit filename, which carries \"phire-region-country-adm1\"). Same reasoning as mon$seasonal_term above -- pre-existing behavior, not one of the six findings, flagged for the record and not fixed in this unit."
  ),

  # -- figures / simulator run-root builders -------------------------------
  list(
    expr = "opt$figdir",
    literal = "figures",
    reason = "get_default_options()$figdir is NULL by design -- figdir is not in the shared top-level key list (config.R's `for (nm in c(\"baseline\", ...))` loop) that make_options_from_config() always fills, so a config that omits it stays NULL and this fallback is the only default that ever applies."
  ),
  list(
    expr = "opt$simulator$run_tag",
    literal = "default",
    reason = "get_default_options()$simulator$run_tag is NULL by design (a directory-naming fallback, not a modeling default); make_sim_run_file() and make_spatpomp_run_file() independently choose the same literal \"default\" for it."
  ),
  list(
    expr = "opt$simulator$run_tag",
    literal = "s10_historical_smoothing",
    reason = "get_default_options()$simulator$run_tag is NULL by design (see the \"default\" entry above); make_smooth_run_file() deliberately chooses a different fallback literal than the v1/spatpomp builders because it roots a different pipeline (08b) with its own conventional directory name."
  ),
  list(
    expr = "opt$simulator$output$format",
    literal = "qs",
    key_path = c("simulator", "output", "format")
  )
)

# DEFAULT_* constants (paths.R:230-234) map onto get_default_options() the
# same way make_occ_run_file()'s own comment already claims ("must mirror
# get_default_options()") -- this gate is what actually checks that claim.
DEFAULT_CONST_MAP <- list(
  DEFAULT_BACKEND = c("occurrence", "backend"),
  DEFAULT_BASEHAZ_POISSON = c("occurrence", "brms", "basehaz_poisson"),
  DEFAULT_COUNTRY = c("country")
)

defaults <- get_default_options()

# --- 1. no NEW %||% site is unmapped -----------------------------------------
cat("\n[1] every extracted %||% \"literal\" site is in SITE_MAP\n")

.site_key <- function(expr, literal) paste0(expr, "\u0001", literal)
mapped_keys <- vapply(
  SITE_MAP,
  function(s) .site_key(s$expr, s$literal),
  character(1)
)
extracted_keys <- .site_key(pct_sites$expr, pct_sites$literal)

new_sites <- setdiff(extracted_keys, mapped_keys)
if (length(new_sites)) {
  cat("        new, unmapped sites (add to SITE_MAP):\n")
  cat(
    paste0("          - ", gsub("\u0001", " %||% \"", new_sites), "\"\n"),
    sep = ""
  )
}
expect(
  length(new_sites) == 0L,
  "no %||% \"literal\" site in paths.R is missing from SITE_MAP"
)

# --- 2. no SITE_MAP entry is stale -------------------------------------------
cat("\n[2] no SITE_MAP entry has gone stale\n")

stale_sites <- setdiff(mapped_keys, extracted_keys)
if (length(stale_sites)) {
  cat(
    "        stale entries (site no longer in paths.R -- delete from SITE_MAP):\n"
  )
  cat(
    paste0("          - ", gsub("\u0001", " %||% \"", stale_sites), "\"\n"),
    sep = ""
  )
}
expect(
  length(stale_sites) == 0L,
  "every SITE_MAP entry still has a matching site in paths.R"
)

# --- 3. mapped sites agree with get_default_options() ------------------------
cat("\n[3] mapped sites agree with get_default_options()\n")

for (s in SITE_MAP) {
  label <- sprintf("%s %%||%% \"%s\"", s$expr, s$literal)
  if (!is.null(s$reason)) {
    expect(TRUE, sprintf("%s: no counterpart (%s)", label, s$reason))
    next
  }
  norm <- if (!is.null(s$normalize)) s$normalize else identity
  if (!is.null(s$key_paths)) {
    ok <- all(vapply(
      s$key_paths,
      function(kp) {
        identical(norm(get_by_path(defaults, kp)), norm(s$literal))
      },
      logical(1)
    ))
    expect(
      ok,
      sprintf(
        "%s: matches get_default_options() at all of %s",
        label,
        paste(
          vapply(s$key_paths, paste, character(1), collapse = "."),
          collapse = ", "
        )
      )
    )
  } else {
    val <- get_by_path(defaults, s$key_path)
    expect(
      identical(norm(val), norm(s$literal)),
      sprintf(
        "%s: matches get_default_options()$%s (got \"%s\")",
        label,
        paste(s$key_path, collapse = "$"),
        format(val)
      )
    )
  }
}

# --- 4. DEFAULT_* constants agree with get_default_options() -----------------
cat(
  "\n[4] DEFAULT_* constants (paths.R:230-234) agree with get_default_options()\n"
)

expect(
  setequal(dflt_sites$expr, names(DEFAULT_CONST_MAP)),
  sprintf(
    "DEFAULT_* constants on disk match DEFAULT_CONST_MAP (on disk: %s | map: %s)",
    paste(sort(dflt_sites$expr), collapse = ","),
    paste(sort(names(DEFAULT_CONST_MAP)), collapse = ",")
  )
)

for (i in seq_len(nrow(dflt_sites))) {
  nm <- dflt_sites$expr[i]
  lit <- dflt_sites$literal[i]
  kp <- DEFAULT_CONST_MAP[[nm]]
  if (is.null(kp)) {
    next
  } # already failed by the setequal check above
  val <- get_by_path(defaults, kp)
  expect(
    identical(val, lit),
    sprintf(
      "%s <- \"%s\": matches get_default_options()$%s (got \"%s\")",
      nm,
      lit,
      paste(kp, collapse = "$"),
      format(val)
    )
  )
}

# --- 5. 04_fit_run.R dispatcher reads the canonical per-track default --------
cat("\n[5] 04_fit_run.R backend fallback reads get_default_options()\n")

fit_run_txt <- paste(
  readLines("analysis/04_fit_run.R", warn = FALSE),
  collapse = "\n"
)
expect(
  grepl(
    "get_default_options\\(\\)\\[\\[model\\]\\]\\$backend",
    fit_run_txt,
    perl = TRUE
  ),
  "04_fit_run.R resolves an omitted --backend via get_default_options()[[model]]$backend"
)
expect(
  !grepl(
    'if\\s*\\(model == "occurrence"\\)\\s*"stansurv"\\s*else\\s*"brms"',
    fit_run_txt,
    perl = TRUE
  ),
  "04_fit_run.R no longer hardcodes the pre-fix two-way (occurrence/else) backend fallback"
)

# =============================================================================
# 6. Simulator modules
# =============================================================================
# Every `<expr> %||% <literal>` in the simulator / smoothing modules whose expr
# resolves to a key that get_default_options() defines must equal that default.
# Literals are strings, numbers (incl. `L`), TRUE/FALSE, NA*, Inf. Unlike
# SITE_MAP above, sites are not listed one by one (there are ~200): the
# expression's ROOT is mapped to a key-path prefix (`flt` -> simulator.filter,
# `sim` -> simulator, `opt` -> top level, ...), and the rest of the `$` chain
# is the key. So the registries are:
#   SIM_EXPR_MAP     exact expressions mapped directly (the smoothing state S$*)
#   SIM_ROOT_MAP     root -> key-path prefix(es) (file-specific or "*")
#   SIM_IGNORED      roots that are not config lists (draw bundles, stamps,
#                    already-resolved function arguments), with the reason
#   SIM_ALLOWLIST    (file, expr, literal) sites allowed to differ, with the reason
# A root in none of them fails [6a]; an allowlist entry whose site is gone, or
# whose key now has an equal default, fails [6c]. A key whose default is NULL
# ("unset") is not compared (its literal is the consumer's own fallback for the
# unset case), except simulator.seed, whose effective default is the top-level
# `seed` (make_options_from_config() copies it into a NULL simulator.seed).
# Sentinel literals (NA*, "", "(none)", ...: helpers_parse.R::NULLOR_SENTINELS)
# are never compared. Sites come from the parse tree (which also catches a
# literal on the line after its `%||%`).
cat("\n[6] simulator-module %||% literals agree with get_default_options()\n")
source("analysis/tests/helpers_parse.R")

# simulator_core_legacy.R: the frozen pre-speed-up kernel copy, same literals.
SIM_FILES <- c(
  "analysis/utils/simulator_core.R",
  "analysis/utils/simulator_core_legacy.R",
  "analysis/utils/spatpomp_core.R",
  "analysis/utils/spatpomp_core_c.R",
  "analysis/utils/smoothing_core.R",
  "analysis/utils/simulations.R",
  "analysis/07b_simulator_v1.R",
  "analysis/07b_simulator_spatpomp.R"
)

# (file, expr, literal) per operand of a literal-terminated `%||%` chain,
# deduplicated per file; `(x %||% list())$k` and `x[["k"]]` read as `x$k`.
sim_sites <- do.call(rbind, lapply(SIM_FILES, nullor_literal_sites))
cat(sprintf(
  "\n  %d (file, expr, literal) sites in %d simulator modules\n",
  nrow(sim_sites),
  length(SIM_FILES)
))

# Exact expressions -> key path. The smoothing state `S` (smoothing_core.R) is
# not a config list, but these fields are set straight from config keys and
# their `%||%` re-states the key's fallback.
SIM_EXPR_MAP <- list(
  `S$BLOCK_GRANULARITY` = c("simulator", "filter", "block_granularity"),
  `S$PARTITION_BY` = c("simulator", "filter", "partition_by"),
  `S$PARTITION_SOURCE` = c("simulator", "filter", "partition_source"),
  `S$COVARIATE_FEED` = c("simulator", "filter", "covariate_feed"),
  `S$MONTHLY_CASE_OBS` = c("simulator", "filter", "monthly_case_obs"),
  `S$NBR_RHO_KM` = c("simulator", "nbr_contagion", "rho_km"),
  `S$NBR_GADM_PATH` = c("simulator", "nbr_contagion", "gadm_path"),
  `S$NBR_CACHE_DIR` = c("simulator", "nbr_contagion", "cache_dir"),
  `S$NBR_PARTITION_SOURCE` = c(
    "simulator",
    "nbr_contagion",
    "partition_source"
  ),
  `S$EPISODE_KEY` = c("episodes", "key"),
  `S$NATIONAL_ZERO_THRESHOLD` = c("episodes", "national_zero_threshold"),
  `S$POP_YEAR` = c("simulator", "location_frame", "pop_year"),
  `S$ONSET_ELIGIBILITY` = c("simulator", "onset_eligibility"),
  `S$CESSATION_EXPOSURE` = c("simulator", "cessation_exposure"),
  `S$DRAW_DESIGN` = c("simulator", "filter", "draw_design")
)

# Root -> key-path prefix(es). file "*" = every simulator module. A root with
# several prefixes is a helper called for several blocks; the literal must
# match all of them (or be allowlisted).
SIM_ROOT_MAP <- list(
  list(file = "*", root = "opt$simulator", prefix = list("simulator")),
  list(file = "*", root = "opt", prefix = list(character(0))),
  # the occurrence fit's own resolved options (loaders: opt_occ$baseline %||% opt$baseline %||% ...)
  list(file = "*", root = "opt_occ", prefix = list(character(0))),
  list(file = "*", root = "sim", prefix = list("simulator")),
  list(file = "*", root = "sim_scope", prefix = list(c("simulator", "scope"))),
  list(file = "*", root = "flt", prefix = list(c("simulator", "filter"))),
  list(
    file = "*",
    root = "lf_opts",
    prefix = list(c("simulator", "location_frame"))
  ),
  list(
    file = "*",
    root = "nbr",
    prefix = list(c("simulator", "nbr_contagion"))
  ),
  list(
    file = "*",
    root = "nbr_cfg",
    prefix = list(c("simulator", "nbr_contagion"))
  ),
  list(
    file = "*",
    root = "nbr_cfg_08",
    prefix = list(c("simulator", "nbr_contagion"))
  ),
  list(file = "*", root = "wu", prefix = list(c("simulator", "warmup"))),
  list(file = "*", root = "dur", prefix = list("duration")),
  list(file = "*", root = "mon", prefix = list("monthly")),
  # .resolve_who_obs_cfg(w = flt$who_obs, ...)
  list(
    file = "spatpomp_core.R",
    root = "w",
    prefix = list(c("simulator", "filter", "who_obs"))
  ),
  # .resolve_who_obs_cfg(): se <- w[["scale_est"]]
  list(
    file = "spatpomp_core.R",
    root = "se",
    prefix = list(c("simulator", "filter", "who_obs", "scale_est"))
  ),
  # smoothing_load_config(): cfg is the raw (unmerged) YAML, top level
  list(file = "smoothing_core.R", root = "cfg", prefix = list(character(0))),
  # extract_monthly_draws()'s shared NB helper: blk <- opt[[block]], block = size | monthly
  list(
    file = "simulator_core.R",
    root = "blk",
    prefix = list("size", "monthly")
  )
)

SIM_IGNORED <- list(
  `S` = "smoothing state (smoothing_default_config()); its config-keyed fields are in SIM_EXPR_MAP",
  `prov` = "population provenance stamp, not config",
  `est_stamp` = "estimand stamp read from a prep artefact, not config",
  `rep_al` = "replicate-alignment result, not config",
  `bundle` = "posterior draws bundle, not config",
  `draws` = "posterior draws bundle, not config",
  `draws_nb` = "posterior draws bundle, not config",
  `draws_mon_sub` = "posterior draws bundle, not config",
  `bh` = "baseline-hazard spec carried by a draws bundle, not config",
  `dur_bh` = "baseline-hazard spec carried by a draws bundle, not config",
  `kind` = "baseline-hazard kind of a draws bundle, not config",
  `b` = "draws-bundle element, not config",
  `d` = "draws-bundle element, not config",
  `r` = "replicate / trajectory record, not config",
  `x` = "trajectory record, not config",
  `rep_id` = "replicate id variable, not config",
  `draw_d` = "draw index variable, not config",
  `replicate_offset` = "function argument (07b chunk offset), resolved by the caller",
  `offset_logoffset` = "function argument (C builder), resolved from the draws bundle",
  `partition_source` = "function argument of smoothing_cache_meta(), resolved by the caller",
  `grouping_source` = "function argument of smoothing_cache_meta(), resolved by the caller",
  `fits_dir` = "function argument of smoothing_cache_meta(), resolved by the caller",
  `got` = "cache meta.yml read back for a comparison, not config",
  `expect` = "expected cache meta for a comparison, not config",
  `surv_path` = "file path variable (\"\" / \"NULL\" are file.exists() / message guards)",
  `args` = "optparse CLI options of 08a/08b (args$chunk_size), not config"
)
SIM_IGNORED_L <- lapply(names(SIM_IGNORED), function(r) {
  list(file = "*", root = r, reason = SIM_IGNORED[[r]])
})

.SIM_ALLOW_SEED <- paste(
  "simulator.seed defaults to NULL and make_options_from_config() copies the top-level `seed` (1234) into it,",
  "so this literal is reached only by a hand-built opt. The two backends' literals disagree (07b spatPomp",
  "20260101, v1 / loaders 20260603); every deployed simulator config sets simulator.seed. Left as is (E4: changing",
  "a fallback is out of scope)."
)
SIM_ALLOWLIST <- list(
  # -- simulator.seed (effective default: top-level seed) ------------------------
  list(
    file = "simulations.R",
    expr = "sim$seed",
    literal = "20260603L",
    reason = .SIM_ALLOW_SEED
  ),
  list(
    file = "07b_simulator_v1.R",
    expr = "sim$seed",
    literal = "20260603L",
    reason = .SIM_ALLOW_SEED
  ),
  list(
    file = "07b_simulator_spatpomp.R",
    expr = "opt$simulator$seed",
    literal = "20260101L",
    reason = .SIM_ALLOW_SEED
  ),
  # -- filter.enable ----------------------------------------------------------------
  list(
    file = "07b_simulator_spatpomp.R",
    expr = "flt$enable",
    literal = "TRUE",
    reason = "get_default_options() follows the only live reader (simulations.R: isTRUE(), i.e. FALSE when unset); this literal feeds `filter_enabled`, which 07b never reads, and is unreachable for a resolved opt. Commented at the site."
  ),
  # -- smoothing raw-config outdir ----------------------------------------------------
  list(
    file = "smoothing_core.R",
    expr = "cfg$outdir",
    literal = "\"generated_data\"",
    reason = "cfg is the raw YAML (live for a config omitting outdir); same directory as the default \"generated_data/\" (no trailing slash: file.path() joins either). Left as is."
  ),
  # -- fit-opt dispatch fallbacks (opt = the fit's own resolved options) --------------
  # (the occurrence backend fallback is a stop(), the monthly one the fitting
  # default brms_cb, 07b's horizon.end fallback get_default_options().)
  list(
    file = "simulator_core.R",
    expr = "blk$phi_re",
    literal = "\"constant\"",
    reason = "shared NB helper for size (default \"constant\") and monthly (default \"region_country_adm1\"); the fit opt always carries phi_re. Same sentinel class as paths.R's mon$phi_re. Left as is."
  )
)

defaults <- get_default_options()

# Effective default at a key path: simulator.seed reads the top-level seed.
.sim_effective <- function(kp) {
  if (identical(kp, c("simulator", "seed"))) {
    return(defaults$seed)
  }
  config_get_leaf(defaults, paste(kp, collapse = "."))
}

#' Classify `sites` and check them against an allowlist: the labels of the
#' unclassified sites, of the mismatches outside the allowlist, and of the
#' allowlist entries whose site is gone or no longer differs (stale).
.nullor_gate <- function(sites, root_map, expr_map, ignored, allowlist) {
  cls <- nullor_classify(
    sites,
    defaults,
    root_map,
    expr_map,
    ignored,
    effective = .sim_effective
  )
  key <- function(file, expr, literal) {
    paste(file, expr, literal, sep = "\u0001")
  }
  ak <- vapply(
    allowlist,
    function(a) key(a$file, a$expr, a$literal),
    character(1)
  )
  sk <- key(cls$file, cls$expr, cls$literal)
  lab <- sprintf(
    "%s:%d %s %%||%% %s%s",
    cls$file,
    cls$line,
    cls$expr,
    cls$literal,
    ifelse(nzchar(cls$detail), sprintf("  [default %s]", cls$detail), "")
  )
  list(
    cls = cls,
    unclassified = lab[cls$status == "UNCLASSIFIED"],
    mismatched = lab[cls$status == "MISMATCH" & !(sk %in% ak)],
    gone = gsub("\u0001", " | ", setdiff(ak, sk)),
    allowed_equal = lab[sk %in% ak & cls$status != "MISMATCH"]
  )
}

g6 <- .nullor_gate(
  sim_sites,
  SIM_ROOT_MAP,
  SIM_EXPR_MAP,
  SIM_IGNORED_L,
  SIM_ALLOWLIST
)
unclassified <- g6$unclassified
mismatched <- g6$mismatched
gone <- g6$gone
allowed_equal <- g6$allowed_equal
st6 <- table(factor(
  g6$cls$status,
  levels = c(
    "MATCH",
    "MISMATCH",
    "UNSET",
    "NOKEY",
    "SENTINEL",
    "IGNORED",
    "BARE",
    "NONCHAIN",
    "UNCLASSIFIED"
  )
))
cat(sprintf(
  "  %s\n",
  paste(sprintf("%s %d", names(st6), as.integer(st6)), collapse = ", ")
))

cat("\n[6a] every site's root is classified\n")
if (length(unclassified)) {
  cat("        unclassified (add the root to SIM_ROOT_MAP or SIM_IGNORED):\n")
  cat(paste0("          - ", unclassified, "\n"), sep = "")
}
expect(
  length(unclassified) == 0L,
  "every simulator-module %||% literal site has a mapped or ignored root"
)

cat("\n[6b] mapped literals equal get_default_options() (or are allowlisted)\n")
if (length(mismatched)) {
  cat("        literal differs from the default:\n")
  cat(paste0("          - ", mismatched, "\n"), sep = "")
}
expect(
  length(mismatched) == 0L,
  "no simulator-module %||% literal disagrees with get_default_options() outside SIM_ALLOWLIST"
)

cat("\n[6c] SIM_ALLOWLIST is not stale\n")
if (length(gone)) {
  cat("        allowlisted site no longer in the code (delete the entry):\n")
  cat(paste0("          - ", gone, "\n"), sep = "")
}
if (length(allowed_equal)) {
  cat("        allowlisted site now equals its default (delete the entry):\n")
  cat(paste0("          - ", allowed_equal, "\n"), sep = "")
}
expect(length(gone) == 0L, "every SIM_ALLOWLIST entry still has its site")
expect(
  length(allowed_equal) == 0L,
  "no SIM_ALLOWLIST entry now equals its default"
)

cat(
  "\n[6d] no pop_year %||% literal remains (tier-2 item 9: one reader, sim_pop_year_spec())\n"
)
pop_sites <- sim_sites[grepl("pop_year$|POP_YEAR$", sim_sites$expr), ]
if (nrow(pop_sites)) {
  cat(
    paste0(
      "          - ",
      pop_sites$file,
      ":",
      pop_sites$line,
      " ",
      pop_sites$expr,
      " %||% ",
      pop_sites$literal,
      "\n"
    ),
    sep = ""
  )
}
expect(
  nrow(pop_sites) == 0L,
  "no `pop_year %||% <literal>` site in the simulator modules (read via sim_pop_year_spec())"
)

cat(
  "\n[6e] extract_occurrence_draws() refuses a fit opt without occurrence.backend (item 17 C3)\n"
)
sc_env <- new.env()
sys.source("analysis/utils/simulator_core.R", sc_env)
e6 <- tryCatch(
  {
    sc_env$extract_occurrence_draws(NULL, list())
    ""
  },
  error = conditionMessage
)
expect(
  grepl("occurrence.backend", e6, fixed = TRUE),
  "an opt with no occurrence.backend is refused (no silent brms fallback)"
)
e6s <- tryCatch(
  {
    sc_env$extract_occurrence_draws(
      NULL,
      list(occurrence = list(backend = "no_such_backend"))
    )
    ""
  },
  error = conditionMessage
)
expect(
  grepl("unknown occurrence backend", e6s, fixed = TRUE),
  "a stated backend still dispatches"
)

# =============================================================================
# 7. Every other module and script
# =============================================================================
# The [6] rule for analysis/*.R and analysis/utils/*.R minus the SIM_FILES:
# every operand of a literal-terminated `%||%` chain whose root maps to a key
# that get_default_options() defines must equal that default. paths.R sites
# that SITE_MAP above already covers ([1]-[3], string literals) are skipped.
# Registries:
#   NULLOR_ROOT_MAP   root -> key-path prefix(es), file-specific or "*"
#   NULLOR_EXPR_MAP   exact expressions -> key path (SIM_EXPR_MAP plus the
#                     locals that hold one config value)
#   NULLOR_IGNORED    roots that are not config lists, each with the reason
#   NULLOR_ALLOWLIST  (file, expr, literal) sites allowed to differ, with the reason
# [7a] every chain operand's root is mapped or ignored; [7b] no mismatch
# outside the allowlist; [7c] the allowlist is not stale; [7d] the pop_year
# sites (routed through sim_frame_pop_args()) leave no
# pop_year literal behind. A key whose default is NULL is not compared, nor a
# sentinel literal. Operands that are no key chain (a call, an index) are not
# config reads and are counted, not checked.
cat(
  "\n[7] %||% literals in every other module / script agree with get_default_options()\n"
)

NULLOR_FILES <- setdiff(
  c(
    list.files("analysis", pattern = "\\.R$", full.names = TRUE),
    list.files("analysis/utils", pattern = "\\.R$", full.names = TRUE)
  ),
  SIM_FILES
)
nullor_sites <- do.call(rbind, lapply(NULLOR_FILES, nullor_literal_sites))
# paths.R string-literal sites are SITE_MAP's (asserted in [1]-[3]).
site_map_keys <- vapply(
  SITE_MAP,
  function(s) paste(s$expr, s$literal, sep = "\u0001"),
  character(1)
)
str_lit <- function(l) {
  v <- nullor_literal_value(str2lang(l))$value
  if (is.character(v)) v else NA_character_
}
in_site_map <- nullor_sites$file == "paths.R" &
  paste(
    nullor_sites$expr,
    vapply(nullor_sites$literal, str_lit, character(1)),
    sep = "\u0001"
  ) %in%
    site_map_keys
nullor_sites <- nullor_sites[!in_site_map, , drop = FALSE]
cat(sprintf(
  "\n  %d (file, expr, literal) sites in %d files (%d paths.R sites left to SITE_MAP)\n",
  nrow(nullor_sites),
  length(unique(nullor_sites$file)),
  sum(in_site_map)
))

.tr <- function(t, ...) {
  lapply(c(t, ...), function(x) strsplit(x, ".", fixed = TRUE)[[1]])
}
.DLNM_DEFAULTS <- .tr(
  "occurrence.dlnm.defaults",
  "duration.dlnm.defaults",
  "monthly.dlnm.defaults",
  "size.dlnm.defaults"
)
.BRMS_ALL <- .tr(
  "occurrence.brms",
  "duration.brms",
  "monthly.brms",
  "size.brms"
)

NULLOR_ROOT_MAP <- list(
  # -- shared names (opt = resolved options, or a raw config / CLI list read at top level)
  list(file = "*", root = "opt", prefix = list(character(0))),
  list(file = "*", root = "occ", prefix = .tr("occurrence")),
  list(file = "*", root = "dur", prefix = .tr("duration")),
  list(file = "*", root = "mon", prefix = .tr("monthly")),
  list(file = "*", root = "sim", prefix = .tr("simulator")),
  list(file = "*", root = "RES", prefix = .tr("results")),
  list(file = "*", root = "STYLE", prefix = .tr("results.style")),
  list(file = "*", root = "ctx$RES", prefix = .tr("results")),
  list(file = "*", root = "ctx$STYLE", prefix = .tr("results.style")),
  # -- population / prep
  list(
    file = "02c_preprocess_population_compile.R",
    root = "pop",
    prefix = .tr("population")
  ),
  list(file = "population_qc.R", root = "opt", prefix = .tr("population")),
  # -- fit scripts: the track's backend blocks
  list(
    file = "04_fit_occ_brms.R",
    root = "brms_cfg",
    prefix = .tr("occurrence.brms")
  ),
  list(
    file = "04_fit_occ_brms_cb.R",
    root = "brms_cfg",
    prefix = .tr("occurrence.brms")
  ),
  list(
    file = "04_fit_occ_brms_kl.R",
    root = "brms_cfg",
    prefix = .tr("occurrence.brms")
  ),
  list(
    file = "04_fit_occ_brms_kl.R",
    root = "kcfg",
    prefix = .tr("occurrence.kernel")
  ),
  list(file = "04_fit_dur_brms.R", root = "dur_cfg", prefix = .tr("duration")),
  list(
    file = "04_fit_dur_brms.R",
    root = "brms_cfg",
    prefix = .tr("duration.brms")
  ),
  list(
    file = "04_fit_dur_brms_cb.R",
    root = "dur_cfg",
    prefix = .tr("duration")
  ),
  list(
    file = "04_fit_dur_brms_cb.R",
    root = "brms_cfg",
    prefix = .tr("duration.brms")
  ),
  list(
    file = "04_fit_dur_brms_cb.R",
    root = "dlnm_defaults",
    prefix = .tr("duration.dlnm.defaults")
  ),
  list(
    file = "04_fit_dur_brms_cb.R",
    root = "var_spec",
    prefix = .tr("duration.dlnm.defaults")
  ),
  list(
    file = "04_fit_mon_brms.R",
    root = "brms_cfg",
    prefix = .tr("monthly.brms")
  ),
  list(
    file = "04_fit_mon_brms_cb.R",
    root = "brms_cfg",
    prefix = .tr("monthly.brms")
  ),
  list(
    file = "04_fit_mon_brms_cb.R",
    root = "dlnm_defaults",
    prefix = .tr("monthly.dlnm.defaults")
  ),
  list(
    file = "04_fit_mon_brms_cb.R",
    root = "var_spec",
    prefix = .tr("monthly.dlnm.defaults")
  ),
  list(
    file = "04_fit_siz_brms.R",
    root = "brms_cfg",
    prefix = .tr("size.brms")
  ),
  list(
    file = "04_fit_siz_brms_cb.R",
    root = "brms_cfg",
    prefix = .tr("size.brms")
  ),
  list(
    file = "04_fit_siz_brms_cb.R",
    root = "dlnm_defaults",
    prefix = .tr("size.dlnm.defaults")
  ),
  list(
    file = "04_fit_siz_brms_cb.R",
    root = "var_spec",
    prefix = .tr("size.dlnm.defaults")
  ),
  # -- postprocess / figures / simulator post (raw or resolved configs)
  list(
    file = "05_fit_postprocess_dependence.R",
    root = "size_opt",
    prefix = list(character(0))
  ),
  list(file = "*", root = "opt_run", prefix = list(character(0))),
  # -- shared formula / prior / path builders (called for several blocks)
  list(
    file = "config.R",
    root = "lf",
    prefix = .tr("simulator.location_frame")
  ),
  # resolve_brms_init(): called by the occurrence brms_cb / brms_kl workers only
  list(file = "config.R", root = "brms_cfg", prefix = .tr("occurrence.brms")),
  list(
    file = "config.R",
    root = "block",
    prefix = .tr("occurrence", "size", "duration", "monthly")
  ),
  list(file = "formula.R", root = "brms_cfg", prefix = .tr("occurrence.brms")),
  list(file = "formula.R", root = "dlnm_defaults", prefix = .DLNM_DEFAULTS),
  list(file = "formula.R", root = "var_spec", prefix = .DLNM_DEFAULTS),
  list(file = "formula.R", root = "blk", prefix = .tr("size", "monthly")),
  list(file = "formula.R", root = "cfg", prefix = .tr("size", "monthly")),
  list(file = "dlnm.R", root = "spec", prefix = .DLNM_DEFAULTS),
  list(file = "priors.R", root = "block", prefix = .tr("size", "monthly")),
  list(file = "paths.R", root = "ep", prefix = .tr("episodes")),
  list(file = "paths.R", root = "brms_sub", prefix = .tr("occurrence.brms")),
  list(file = "paths.R", root = "lr", prefix = .tr("occurrence.kernel.lag_re")),
  list(file = "paths.R", root = "blk", prefix = .tr("size", "monthly")),
  list(
    file = "paths.R",
    root = "re",
    prefix = .tr("size.random_effects", "monthly.random_effects")
  ),
  list(
    file = "paths.R",
    root = "br",
    prefix = .tr(
      "occurrence.base_rate",
      "duration.base_rate",
      "monthly.base_rate"
    )
  ),
  list(
    file = "paths.R",
    root = "av",
    prefix = .tr(
      "occurrence.dlnm.defaults.argvar",
      "duration.dlnm.defaults.argvar",
      "monthly.dlnm.defaults.argvar"
    )
  ),
  list(file = "anchor.R", root = "cfg", prefix = .tr("size.mean_rate")),
  list(file = "provenance.R", root = "cfg", prefix = list(character(0))),
  list(file = "provenance.R", root = "o", prefix = list(character(0))),
  list(file = "provenance.R", root = "blk", prefix = .tr("size", "monthly")),
  list(
    file = "provenance.R",
    root = "es",
    prefix = .tr("size.estimand", "monthly.estimand")
  ),
  # -- brms_kl kernel blocks
  list(file = "*", root = "kcfg", prefix = .tr("occurrence.kernel")),
  list(file = "dlnm_kernel.R", root = "k", prefix = .tr("occurrence.kernel")),
  list(file = "kl_sweep.R", root = "kc", prefix = .tr("occurrence.kernel")),
  list(
    file = "*",
    root = "cg",
    prefix = .tr("occurrence.kernel.prior_grouping.contiguity")
  ),
  list(file = "*", root = "dc", prefix = .tr("occurrence.kernel.discrete")),
  list(file = "*", root = "pn", prefix = .tr("occurrence.kernel.pin")),
  list(
    file = "dlnm_kernel.R",
    root = "pc",
    prefix = .tr("occurrence.kernel.pin.pilot")
  ),
  list(
    file = "dlnm_kernel.R",
    root = "hs",
    prefix = .tr("occurrence.kernel.hs")
  ),
  list(
    file = "dlnm_kernel.R",
    root = "lr",
    prefix = .tr("occurrence.kernel.lag_re")
  ),
  list(
    file = "dlnm_kernel.R",
    root = "brms_cfg",
    prefix = .tr("occurrence.brms")
  ),
  list(file = "dlnm_kernel.R", root = "dep", prefix = .tr("deploy")),
  list(
    file = "kl_sweep.R",
    root = "argvar",
    prefix = .tr("occurrence.dlnm.defaults.argvar")
  ),
  # -- deploy
  list(file = "deploy_core.R", root = "dep", prefix = .tr("deploy")),
  list(file = "deploy_core.R", root = "cv", prefix = .tr("deploy.convergence")),
  list(
    file = "deploy_core.R",
    root = "ls_cfg",
    prefix = .tr("deploy.lag_selection")
  ),
  list(
    file = "deploy_core.R",
    root = "ea",
    prefix = .tr("deploy.event_attribution")
  ),
  list(file = "deploy_core.R", root = "o1", prefix = list(character(0))),
  list(file = "deploy_core.R", root = "t", prefix = list(character(0))),
  list(file = "deploy_prod_run.R", root = "DEP", prefix = .tr("deploy")),
  list(file = "deploy_prod_run.R", root = "dep", prefix = .tr("deploy")),
  list(
    file = "deploy_prod_run.R",
    root = "ea",
    prefix = .tr("deploy.event_attribution")
  ),
  list(
    file = "deploy_prod_run.R",
    root = "ks",
    prefix = .tr("deploy.kl_sweep")
  ),
  list(
    file = "deploy_prod_run.R",
    root = ".hc_sw",
    prefix = .tr("deploy.hindcast.season_window")
  ),
  list(file = "deploy_prod_run.R", root = "thist", prefix = list(character(0))),
  list(
    file = "kl_grouping.R",
    root = "cfg",
    prefix = .tr("deploy.convergence")
  ),
  list(file = "loo_lags.R", root = "x", prefix = .tr("deploy.lag_selection")),
  # -- config generators: templates / generated configs are full sim configs;
  #    opt$enso_continuation is the deploy block passed through
  list(
    file = "make_forecast_configs.R",
    root = "factual_tpl",
    prefix = list(character(0))
  ),
  list(
    file = "make_forecast_configs.R",
    root = "opt$enso_continuation",
    prefix = .tr("deploy.enso_continuation")
  ),
  list(
    file = "make_forecast_configs.R",
    root = "cont",
    prefix = .tr("deploy.enso_continuation")
  ),
  list(
    file = "make_hindcast_configs.R",
    root = "tpl",
    prefix = list(character(0))
  ),
  list(
    file = "make_hindcast_configs.R",
    root = "cfg0",
    prefix = list(character(0))
  ), # tpl + episodes key
  list(
    file = "make_hindcast_configs.R",
    root = "cfg",
    prefix = list(character(0))
  ), # hindcast_param_draw_mode(): a patched sim config
  list(
    file = "make_event_attribution_configs.R",
    root = "tpl",
    prefix = list(character(0))
  ),
  list(
    file = "make_event_attribution_configs.R",
    root = "cfg",
    prefix = list(character(0))
  ),
  list(
    file = "make_event_attribution_configs.R",
    root = "cf_block",
    prefix = .tr("simulator.enso.counterfactual")
  ),
  list(
    file = "make_forecast_attrib_configs.R",
    root = "cfg",
    prefix = list(character(0))
  ),
  list(
    file = "make_batch_configs.R",
    root = "cfg",
    prefix = list(character(0))
  ),
  # -- ENSO
  list(
    file = "enso.R",
    root = "cont",
    prefix = .tr("simulator.enso.continuation")
  ),
  list(
    file = "enso.R",
    root = "cf",
    prefix = .tr("simulator.enso.counterfactual")
  ),
  list(
    file = "enso_continuation_skill.R",
    root = "d",
    prefix = .tr("simulator.enso.continuation")
  ),
  list(
    file = "spatial_contagion.R",
    root = "cfg",
    prefix = .tr("simulator.nbr_contagion")
  ),
  # -- results
  list(
    file = "results_derive.R",
    root = "DERIVATION",
    prefix = .tr("results.derivation")
  ),
  list(
    file = "results_derive.R",
    root = "THRESHOLDS",
    prefix = .tr("results.thresholds")
  ),
  list(
    file = "results_derive.R",
    root = "INPUTS",
    prefix = .tr("results.inputs")
  ),
  list(
    file = "results_derive.R",
    root = "africa",
    prefix = .tr("results.scope.africa")
  ),
  list(file = "results_derive.R", root = "fac_cfg", prefix = .tr("simulator")),
  list(
    file = "results_derive.R",
    root = "cont",
    prefix = .tr("simulator.enso.continuation")
  ),
  list(
    file = "results_derive.R",
    root = "sim1",
    prefix = .tr("simulator.enso")
  ),
  # run_design_row() (D34): a run config's simulator.enso block.
  list(
    file = "results_derive.R",
    root = "enso",
    prefix = .tr("simulator.enso")
  ),
  # .sim_param_draw_mode(): a run config's simulator.forecast block.
  list(
    file = "results_derive.R",
    root = "fc",
    prefix = .tr("simulator.forecast")
  ),
  list(file = "results_report.R", root = "res", prefix = .tr("results")),
  list(
    file = "hindcast_skill.R",
    root = "cfg",
    prefix = .tr("results.derivation.hindcast")
  ),
  list(file = "hindcast_skill.R", root = "op", prefix = .tr("episodes"))
)

# Exact expressions: the smoothing state S$* (as in [6]) and locals holding
# one config value.
NULLOR_EXPR_MAP <- c(
  SIM_EXPR_MAP,
  list(
    `S$WHO_OBS$db_overlap` = c("simulator", "filter", "who_obs", "db_overlap"),
    cfg_chunk = c("simulator", "parallel", "chunk_size"), # 08a: cfg$simulator$parallel$chunk_size
    kfix = c("occurrence", "kernel", "prior_grouping", "kappa"), # dlnm_kernel.R
    tpl_nzt = c("episodes", "national_zero_threshold")
  )
) # provenance.R: cfg$episodes$...

# File-keyed, so a new file reusing one of these names is classified afresh.
.ig <- function(file, roots, reason) {
  lapply(roots, function(r) list(file = file, root = r, reason = reason))
}
.R_CONTEXT <- "results / smoothing context environment; its config-keyed fields are mapped (ctx$RES, ctx$STYLE; S$* in NULLOR_EXPR_MAP)"
.R_CLI <- "optparse CLI options of the script, not config"
.R_MANIFEST <- "forecast / hindcast / event manifest entry or cache metadata, not config"
.R_STAMP <- "prep / provenance stamp, registry sidecar, recorded fit metadata or task YAML, not config"
.R_RECORD <- "artefact, fit, metrics or table record, or a generic local, not config"
.R_ARG <- "function argument or local resolved by the caller (the config value is threaded in and its default applied there)"
.R_BATCH <- "batch-generator spec / grid row or its defaults list (make_batch_configs.R's own registry), not get_default_options()"
NULLOR_IGNORED <- c(
  # Classified sites:
  .ig("deploy_core.R", "v", .R_RECORD),
  .ig("make_event_attribution_configs.R", "ev", .R_MANIFEST),
  .ig("results_derive.R", "w", .R_MANIFEST),
  .ig("08a_smoothing_tasks.R", c("S", "args"), .R_CONTEXT),
  .ig("08b_smoothing_spatpomp.R", "S", .R_CONTEXT),
  .ig("08c_smoothing_compile.R", "S", .R_CONTEXT),
  .ig("results_derive.R", "ctx", .R_CONTEXT),
  .ig("05_fit_postprocess.R", "script_opt", .R_CLI),
  .ig("05_fit_postprocess_lag_compare_surface.R", "script_opt", .R_CLI),
  .ig("05_fit_postprocess_loo_lags.R", "script_opt", .R_CLI),
  .ig("06_fit_figures_lag_compare.R", "script_opt", .R_CLI),
  .ig("07c_simulator_postprocess.R", c("man", "ye"), .R_MANIFEST),
  .ig("07d_simulator_figures_counterfactual.R", "man", .R_MANIFEST),
  .ig(
    "results_derive.R",
    c("man", "u", "e", "em", "hc", "si", "win", "m", "def"),
    .R_MANIFEST
  ),
  .ig("results_core.R", c("man_attrib", "gi", "e", "r"), .R_MANIFEST),
  .ig("sim_io.R", c("man", "u", "expected", "pm_d", "r"), .R_MANIFEST),
  .ig("make_event_attribution_configs.R", c("e", "m"), .R_MANIFEST),
  .ig("make_forecast_attrib_configs.R", "u", .R_MANIFEST),
  .ig("hindcast_skill.R", c("idx", "run"), .R_MANIFEST),
  .ig("03_prepare_dur.R", "old", .R_STAMP),
  .ig("03_prepare_occ.R", "nz_stamp", .R_STAMP),
  .ig("03_prepare_lenssen_analogs.R", "x", .R_STAMP),
  .ig("deploy_core.R", "st", .R_STAMP),
  .ig("occ_data.R", "st", .R_STAMP),
  .ig("kl_grouping.R", c("st", "settings", "meta", "r"), .R_STAMP),
  .ig("provenance.R", c("st", "stamp", "a", "b", "r", "estimand"), .R_STAMP),
  .ig("results_derive.R", "st", .R_STAMP),
  .ig("loo_lags.R", "span", .R_STAMP),
  .ig("figures_kl.R", "meta", .R_STAMP),
  .ig("config.R", "side", .R_STAMP),
  .ig("08b_smoothing_spatpomp.R", c("task", "got_dg", "r"), .R_STAMP),
  .ig("07c_simulator_postprocess.R", c("cond", "cg"), .R_RECORD),
  .ig("07d_simulator_figures.R", c("cond", "contagion"), .R_RECORD),
  .ig(
    "results_derive.R",
    c("fam", "spec_track", "x", "cf", "n_smoothing", "clamp_margin"),
    .R_RECORD
  ),
  .ig("diagnostics.R", "fit", .R_RECORD),
  .ig("kl_sweep.R", "fit_sub", .R_RECORD),
  .ig("10_manuscript_figures.R", c("f23", "oni"), .R_RECORD),
  .ig("figures_manuscript_fig34.R", "f23", .R_RECORD), # 10's Figure 3/4 builders
  .ig("10_manuscript_figures_legacy.R", "f23", .R_RECORD),
  .ig("make_pipeline_docs.R", "s", .R_RECORD),
  .ig("deploy_prod_run.R", "res", .R_RECORD),
  .ig(
    "enso.R",
    c(
      "lib",
      "analog_library",
      "library",
      "path",
      "max_lead",
      "max_per_model",
      "analog_max_per_model"
    ),
    .R_ARG
  ),
  .ig("dlnm_kernel.R", c("x", "threads"), .R_ARG),
  .ig("national_zero.R", c("x", "rule"), .R_ARG),
  .ig(
    "make_event_attribution_configs.R",
    c("x", "psy", "mode", "n_smoothing"),
    .R_ARG
  ),
  .ig("make_hindcast_configs.R", "n", .R_ARG), # hindcast_n_replicates(): its argument, 0 = keep
  .ig(
    "deploy_prod_run.R",
    "o",
    "generator option list built by the deploy tool (.forecast_opt() / .event_opt(); CLI-shaped), not a config block"
  ),
  .ig(
    "make_forecast_attrib_configs.R",
    "o",
    "generator option list (deploy_prod_run.R's .forecast_opt()), not a config block"
  ),
  .ig("paths.R", "threshold", .R_ARG),
  .ig("base.R", "workers", .R_ARG),
  .ig("gid.R", "admin_level_opt", .R_ARG),
  .ig("formula.R", "term", .R_ARG),
  .ig("priors.R", c("nlpar", "phi_re"), .R_ARG),
  .ig("loo_lags.R", "group", .R_ARG),
  .ig(
    "sim_io.R",
    c("path", "fmt", "out_fmt", "smoothing_input", "cfg"),
    .R_ARG
  ),
  .ig("season_calendar.R", "path", .R_ARG),
  .ig("results_core.R", c("mode", "series", "path"), .R_ARG),
  .ig("provenance.R", "national_zero_threshold", .R_ARG),
  .ig(
    "kl_grouping.R",
    c("fill_from", "meta_from", "cell", "chains", "n_agree", "mode_tol"),
    .R_ARG
  ),
  .ig("make_forecast_configs.R", "cf_panel", .R_ARG),
  .ig("07c_simulator_postprocess.R", "monthly_data_path", .R_ARG),
  .ig("06_fit_figures_cb_country.R", ".cfg_path", .R_ARG),
  .ig("06_fit_figures_lag_compare.R", "model", .R_ARG),
  .ig("02a_preprocess_population_tasks.R", "reg", .R_ARG),
  .ig(
    "config.R",
    "defaults",
    "the track's defaults list inside .fill_model_block() (itself get_default_options()[[model]])"
  ),
  .ig(
    "make_batch_configs.R",
    c("defaults", "row", "spec", "brms", "brms_iter"),
    .R_BATCH
  )
)

.R_LEGACY <- "legacy meaning of a missing key in old fits: a block written before the key existed had no such term, so the fallback names that, not the current default (same class as SITE_MAP's mon$seasonal_term / mon$phi_re in paths.R). Left as is."
.R_PENDING <- "found by T3-B, decision pending"
NULLOR_ALLOWLIST <- list(
  # -- legacy-fit group ---------------------------------
  list(
    file = "formula.R",
    expr = "cfg$phi_re",
    literal = "\"constant\"",
    reason = .R_LEGACY
  ),
  list(
    file = "priors.R",
    expr = "block$seasonal_term",
    literal = "\"none\"",
    reason = .R_LEGACY
  ),
  list(
    file = "mon_data.R",
    expr = "mon$seasonal_term",
    literal = "\"none\"",
    reason = .R_LEGACY
  ),
  # -- unreachable / out-of-scope fallbacks ----------------------------------------------
  list(
    file = "make_batch_configs.R",
    expr = "cfg$size$backend",
    literal = "\"brms_cb\"",
    reason = "summary-table operand (`cfg$occurrence$backend %||% cfg$size$backend %||% \"brms_cb\"`): the literal is the occurrence default, reached only by a config with neither block. Left as is."
  ),
  list(
    file = "config.R",
    expr = "block$backend",
    literal = "\"brms_cb\"",
    reason = ".fill_model_block(): last fallback after the track's own defaults$backend, which get_default_options() always sets; \"brms_cb\" is the occurrence default, so the literal is unreachable for every track. Left as is."
  ),
  list(
    file = "results_derive.R",
    expr = "fac_cfg$seed",
    literal = "20260603L",
    reason = .SIM_ALLOW_SEED
  ),
  # -- live drift, reported, not fixed ------------------------
  list(
    file = "08a_smoothing_tasks.R",
    expr = "cfg_chunk",
    literal = "10L",
    reason = paste(
      .R_PENDING,
      "-- 08a reads the RAW smoothing YAML, so a config without simulator.parallel.chunk_size gets 10 here while get_default_options() says 25 (live for such a config)."
    )
  )
)


g7 <- .nullor_gate(
  nullor_sites,
  NULLOR_ROOT_MAP,
  NULLOR_EXPR_MAP,
  NULLOR_IGNORED,
  NULLOR_ALLOWLIST
)
st7 <- table(factor(g7$cls$status, levels = names(st6)))
cat(sprintf(
  "  %s\n",
  paste(sprintf("%s %d", names(st7), as.integer(st7)), collapse = ", ")
))

cat("\n[7a] every site's root is classified\n")
if (length(g7$unclassified)) {
  cat(
    "        unclassified (add the root to NULLOR_ROOT_MAP or NULLOR_IGNORED):\n"
  )
  cat(paste0("          - ", g7$unclassified, "\n"), sep = "")
}
expect(
  length(g7$unclassified) == 0L,
  "every %||% literal site outside the simulator modules has a mapped or ignored root"
)

cat("\n[7b] mapped literals equal get_default_options() (or are allowlisted)\n")
if (length(g7$mismatched)) {
  cat("        literal differs from the default:\n")
  cat(paste0("          - ", g7$mismatched, "\n"), sep = "")
}
expect(
  length(g7$mismatched) == 0L,
  "no %||% literal outside the simulator modules disagrees with get_default_options() outside NULLOR_ALLOWLIST"
)

cat("\n[7c] NULLOR_ALLOWLIST is not stale\n")
if (length(g7$gone)) {
  cat("        allowlisted site no longer in the code (delete the entry):\n")
  cat(paste0("          - ", g7$gone, "\n"), sep = "")
}
if (length(g7$allowed_equal)) {
  cat(
    "        allowlisted site no longer differs from its default (delete the entry):\n"
  )
  cat(paste0("          - ", g7$allowed_equal, "\n"), sep = "")
}
expect(length(g7$gone) == 0L, "every NULLOR_ALLOWLIST entry still has its site")
expect(
  length(g7$allowed_equal) == 0L,
  "no NULLOR_ALLOWLIST entry now equals its default"
)

cat(
  "\n[7d] the former pop_year sites read the population year through sim_frame_pop_args()\n"
)
pop7 <- nullor_sites[grepl("(^|\\$)(pop_year|POP_YEAR)$", nullor_sites$expr), ] # not rate_pop_year
if (nrow(pop7)) {
  cat(
    paste0(
      "          - ",
      pop7$file,
      ":",
      pop7$line,
      " ",
      pop7$expr,
      " %||% ",
      pop7$literal,
      "\n"
    ),
    sep = ""
  )
}
expect(
  nrow(pop7) == 0L,
  "no `pop_year %||% <literal>` site outside the simulator modules"
)
fpa <- callers_of(
  "sim_frame_pop_args",
  c("analysis/07a_simulator_tasks.R", "analysis/07c_simulator_postprocess.R")
)
expect(
  setequal(
    unique(basename(fpa$file[fpa$kind == "call"])),
    c("07a_simulator_tasks.R", "07c_simulator_postprocess.R")
  ),
  "07a and 07c (the two former pop_year literal sites outside [6]) call sim_frame_pop_args()"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
