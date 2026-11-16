# utils/config.R — Config-driven option parsing for the cholera pipeline
#
# Three entry points used by all run + postprocess scripts:
#
#   make_options_from_config(option_list, enforce_options = NULL,
#                             cli_argv = NULL, model = NULL)
#     Parse CLI, merge with YAML if --config is given, fill remaining
#     keys from defaults. Returns a complete named list. `model` hints
#     which model block ("occurrence" | "size") to construct when the
#     YAML has no model block (pure-CLI invocation).
#
#   get_default_options()
#     Returns the canonical defaults list (mirrors config_defaults.yml).
#
#   print_options(opt)
#     Pretty-print the resolved options at the top of each run log.
#
# Schema: top-level shared keys + exactly one of occurrence: | size:.
# Convention 2: active slope RE flags imply climate fixed effects.
# Reference: uvira-modeling/analysis/config.R (same pattern).
#
# This module is the OPTIONS CORE (defaults, YAML merge, CLI routing, grouping,
# dotted-leaf helpers, simulator CLI helpers). Model-term assembly lives in
# sibling modules:
#   utils/formula.R         build_formula_terms(), build_{size,duration,monthly}_formula_parts() & co.
#   utils/formula_legacy.R  build_formula_terms_legacy() (frozen, side-by-side)
#   utils/priors.R          build_{occ,size,duration,monthly}_priors() & co.
# Callers that build formulas/priors source those modules next to this one.
# Dependencies (caller-sourced, never source()d here): base.R (`%||%`);
# deploy_core.R at call time only (make_options_from_config() ->
# resolve_results_inputs_from_deploy(), exists()-guarded).


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


#' Return the value following a command-line flag, or NULL if absent.
#'
#' Shared helper for the simulator run scripts (07b_simulator_v1.R,
#' 07a_simulator_tasks.R, 07b_simulator_spatpomp.R), which each need to
#' re-read a handful of CLI flags that make_options_from_config() routed to a
#' top-level key. Returns NULL when the flag is absent or has no following
#' token, so callers can `%||%`-default it.
#'
#' @param flag  flag string to look up (e.g. "--n_replicates")
#' @param args  character vector of CLI args (default: commandArgs(trailingOnly))
#' @return the token after the flag, or NULL
.flag_val <- function(flag, args = commandArgs(trailingOnly = TRUE)) {
  idx <- which(args == flag)
  if (length(idx) == 0L || idx[1] == length(args)) return(NULL)
  args[idx[1] + 1L]
}

#' Common optparse options shared by the simulator offspring runners
#' (07b_simulator_v1.R, 07b_simulator_spatpomp.R). Each offspring extends this with
#' its own backend-specific flags.
#'
#' @return list of optparse make_option() objects
simulator_common_options <- function() {
  list(
    optparse::make_option(c("--dry-run"),
      default = FALSE, action = "store_true", dest = "dry_run",
      help = "Validate config + artefacts + paths, then exit without running."),
    optparse::make_option(c("--overwrite"),
      default = NULL, action = "store_true",
      help = "Overwrite existing simulator outputs."),
    optparse::make_option(c("--n_replicates"),
      default = NULL, type = "integer",
      help = "Override simulator.n_replicates."),
    optparse::make_option(c("--workers"),
      default = NULL, type = "integer",
      help = "Override simulator.parallel.workers.")
  )
}

#' Apply the CLI overrides common to both simulator offspring for the nested
#' `simulator` keys that make_options_from_config() routed to top-level flags.
#'
#' @param opt  resolved options list (must contain `simulator`)
#' @return opt with simulator.n_replicates / simulator.parallel.workers overridden
apply_simulator_common_cli <- function(opt) {
  if (!is.null(.flag_val("--n_replicates")))
    opt$simulator$n_replicates <- as.integer(.flag_val("--n_replicates"))
  if (!is.null(.flag_val("--workers")))
    opt$simulator$parallel$workers <- as.integer(.flag_val("--workers"))
  opt
}

#' Guard: the config must carry a `simulator:` block.
#'
#' @param opt  resolved options list
require_simulator_block <- function(opt) {
  if (is.null(opt$simulator))
    stop("Config must contain a `simulator:` block. Got: ",
         paste(names(opt), collapse = ", "))
  invisible(opt)
}

#' Validate basehaz value against the allowed set for the active backend.
#'
#' @param opt  resolved options list
#' @return opt (invisibly); stops on invalid value
.validate_basehaz <- function(opt) {
  poisson_ok  <- c("factor", "cyclic", "seasonal", "ms")

  occ <- opt$occurrence
  if (!is.null(occ)) {
    be <- occ$backend %||% "brms_cb"
    if (be == "brms") {
      bh <- occ$brms$basehaz_poisson %||% "seasonal"
      if (!bh %in% poisson_ok)
        stop(sprintf("Invalid occurrence.brms.basehaz_poisson '%s'. Must be one of: %s",
                     bh, paste(poisson_ok, collapse = ", ")))
    } else if (be %in% c("brms_cb", "brms_kl")) {
      # brms siblings: their own block, else the plain brms block
      bh <- occ[[be]]$basehaz_poisson %||% occ$brms$basehaz_poisson %||% "seasonal"
      if (!bh %in% poisson_ok)
        stop(sprintf("Invalid occurrence.%s.basehaz_poisson '%s'. Must be one of: %s",
                     be, bh, paste(poisson_ok, collapse = ", ")))
    }
  }

  # Duration track: two independent baseline hazards on two different clocks.
  # Both live at model level (not in the backend sub-block) because the track
  # has one backend, so there is nothing to disambiguate.  "none" is admissible
  # for the seasonal one and, as a sensitivity variant only, for mso.
  dur <- opt$duration
  if (!is.null(dur)) {
    dur_ok <- c("factor", "ms", "none")
    bhd <- dur$basehaz_duration %||% "factor"
    if (!bhd %in% dur_ok)
      stop(sprintf("Invalid duration.basehaz_duration '%s'. Must be one of: %s",
                   bhd, paste(dur_ok, collapse = ", ")))
    bhp <- dur$basehaz_poisson %||% "none"
    if (!bhp %in% c("none", poisson_ok))
      stop(sprintf("Invalid duration.basehaz_poisson '%s'. Must be one of: %s",
                   bhp, paste(c("none", poisson_ok), collapse = ", ")))
    brk <- as.numeric(unlist(dur$mso_breaks %||% c(1, Inf)))
    if (length(brk) < 2L || !is.infinite(brk[length(brk)]) ||
        is.unsorted(brk, strictly = TRUE))
      stop("duration.mso_breaks must be strictly increasing and end in .inf; got: ",
           paste(brk, collapse = ", "))
    tr <- dur$ttlo_reset %||% "onset"
    if (!tr %in% c("onset", "end"))
      stop(sprintf("Invalid duration.ttlo_reset '%s'. Must be one of: onset, end", tr))
    st <- dur[["seasonal_term"]] %||% "none"
    if (!st %in% c("none", "harmonic2"))
      stop(sprintf("Invalid duration.seasonal_term '%s'. Must be one of: none, harmonic2", st))
  }

  # Monthly track. `estimand` + `pop_elasticity` are supported (the size
  # track's denominator-coherence apparatus; see
  # .resolveTrackEstimand()/.ESTIMAND_SPEC in utils/provenance.R). Both
  # default to footprint/fixed.
  #
  # `offset_meanrate` remains refused: mon's anchor is `base_rate_mon`
  # (offset_baserate, a rate per person-day) and it has no `mean_rate` table at
  # all, so the key would name a column that is never built. Refusing loudly
  # beats a silent no-op in the size builders.
  mon <- opt$monthly
  if (!is.null(mon)) {
    unsupported <- intersect(c("offset_meanrate"), names(mon))
    if (length(unsupported) > 0L)
      stop(sprintf(
        "monthly.%s is not supported: mon's baseline-rate anchor is `base_rate_mon` (monthly.offset_baserate), not size's `mean_rate` table.",
        paste(unsupported, collapse = ", monthly.")))
  }

  # Simulator: the size track is retired downstream -- the monthly
  # (`mon`) track's cases_month is the sole case estimand. A simulator config
  # still naming a size fit, its cross-basis sidecar, a `simulator.size:` block
  # or the size `mean_rate` table would otherwise be a silent no-op, so refuse
  # it. Only the `simulator:` block is checked: the fit-side top-level `size:`
  # model block (04_fit_siz_*) is untouched. A YAML `~` reads as NULL and
  # passes. `posterior_strategy: independent` permuted size draws; only
  # "paired" remains.
  sim <- opt$simulator
  if (!is.null(sim)) {
    # `[[` (exact matching), not `$`: sim$size would partial-match size_fit.
    retired <- c(
      "simulator.size_fit"     = !is.null(sim[["size_fit"]]),
      "simulator.size_cb_meta" = !is.null(sim[["size_cb_meta"]]),
      "simulator.size"         = !is.null(sim[["size"]]),
      "simulator.location_frame.mean_rate_table_path" =
        !is.null((sim[["location_frame"]] %||% list())[["mean_rate_table_path"]]))
    if (any(retired))
      stop(sprintf(paste0(
        "%s: the size track is retired downstream of fitting (2026-10-06); the ",
        "simulator models case counts with the monthly track only ",
        "(simulator.monthly_fit / grouping.monthly_enable). Delete the key(s) -- ",
        "for a generated config, regenerate it (deploy_prod_run.R --stage ",
        "gen_templates / gen_forecast / gen_hindcast)."),
        paste(names(retired)[retired], collapse = ", ")), call. = FALSE)
    ps <- sim$posterior_strategy
    if (!is.null(ps) && !identical(ps, "paired"))
      stop(sprintf(paste0(
        "simulator.posterior_strategy = '%s' is not supported: only \"paired\" ",
        "remains (the `independent` branch permuted the retired size track's ",
        "draws)."), paste(ps, collapse = ",")), call. = FALSE)
  }
  invisible(opt)
}

#' Deep-merge two lists: fill any missing keys in `x` from `y` (one level deep).
#'
#' @param x  target list (partial)
#' @param y  defaults list (complete)
#' @return   merged list
.merge1 <- function(x, y) {
  if (is.null(x)) return(y)
  for (nm in names(y)) {
    if (is.null(x[[nm]])) x[[nm]] <- y[[nm]]
  }
  x
}

#' Recursive variant of .merge1: fills missing keys at any nesting depth.
#' Used for the simulator block which has 2-3 level nesting (horizon, scope,
#' output, parallel) and no backend variants.
#'
#' @param x  target list (partial; may be NULL)
#' @param y  defaults list (complete)
#' @return   merged list
.merge_recursive <- function(x, y) {
  if (is.null(x)) return(y)
  if (!is.list(x) || !is.list(y)) return(x)
  for (nm in names(y)) {
    if (is.null(x[[nm]])) {
      x[[nm]] <- y[[nm]]
    } else if (is.list(y[[nm]]) && is.list(x[[nm]])) {
      x[[nm]] <- .merge_recursive(x[[nm]], y[[nm]])
    }
  }
  x
}

# ---------------------------------------------------------------------------
# Dotted-leaf helpers over nested config lists. Shared by make_batch_configs.R
# (grid overrides), deploy_core.R (deploy-key lookup, sim-template overlay and
# leaf diffs) and the config-coverage tests. A dotted key addresses nesting:
# "simulator.filter.measurement" -> x[["simulator"]][["filter"]][["measurement"]].
# ---------------------------------------------------------------------------

#' Get the value at a dotted path, or NULL if any level is missing.
config_get_leaf <- function(x, dotted) {
  Reduce(function(x, k) if (is.list(x)) x[[k]] else NULL,
         strsplit(dotted, ".", fixed = TRUE)[[1]], x)
}

#' Set the value at a dotted path, creating intermediate lists as needed.
#' `value = NULL` deletes the leaf (base-R `[[<-` semantics).
config_set_leaf <- function(x, dotted, value) {
  parts <- strsplit(dotted, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1L) {
    x[[dotted]] <- value
    return(x)
  }
  head_key <- parts[[1]]
  if (is.null(x[[head_key]])) x[[head_key]] <- list()
  x[[head_key]] <- config_set_leaf(x[[head_key]],
                                   paste(parts[-1L], collapse = "."), value)
  x
}

#' Every leaf path of a nested list, as dot-separated strings. A "leaf" is
#' anything that is not a non-empty named list, so a character vector, an
#' unnamed list (YAML sequence) and an empty mapping (`{}`) each stay one leaf.
config_leaf_paths <- function(x, prefix = character(0)) {
  if (!is.list(x) || is.null(names(x)) || !length(x)) {
    return(paste(prefix, collapse = "."))
  }
  unlist(lapply(names(x), function(nm) config_leaf_paths(x[[nm]], c(prefix, nm))),
         use.names = FALSE)
}

#' Flatten a nested list to a named list of dotted leaf path -> value.
config_flatten_leaves <- function(x) {
  paths <- config_leaf_paths(x)
  paths <- paths[nzchar(paths)]
  stats::setNames(lapply(paths, function(p) config_get_leaf(x, p)), paths)
}

#' Registry of the shared top-level config blocks (everything that is not a
#' model block).  `make_options_from_config()` merges these in a loop, so
#' adding a new top-level block is a one-line edit here plus the matching entry
#' in `get_default_options()` — no change to the merge code itself.
#'
#' Fields:
#'   merge            "flat"      -> .merge1        (one level deep)
#'                    "recursive" -> .merge_recursive (any nesting depth)
#'   fill_when_absent TRUE  -> block is materialised from defaults even when the
#'                             YAML omits it entirely
#'                    FALSE -> a missing block stays NULL (used by `simulator`,
#'                             where downstream code tests `is.null(opt$simulator)`
#'                             to decide whether this is a simulator run at all)
#'
#' Completeness against `get_default_options()` is asserted by
#' analysis/tests/test_config_block_registry.R — a block added to the defaults
#' but not here would silently never be filled.
.CONFIG_TOPLEVEL_BLOCKS <- list(
  enso       = list(merge = "flat",      fill_when_absent = TRUE),
  contagion  = list(merge = "flat",      fill_when_absent = TRUE),
  grouping   = list(merge = "flat",      fill_when_absent = TRUE),
  population = list(merge = "recursive", fill_when_absent = TRUE),
  episodes   = list(merge = "flat",      fill_when_absent = TRUE),
  simulator  = list(merge = "recursive", fill_when_absent = FALSE),
  results    = list(merge = "recursive", fill_when_absent = TRUE),
  deploy     = list(merge = "recursive", fill_when_absent = TRUE)
)

#' Fill a model block (occurrence | size) from its defaults.
#' Three levels: model → covariates/random_effects (deep) + backend sub-block.
#'
#' @param block    list (from YAML or empty list)
#' @param defaults list (from get_default_options()[[model]])
#' @return         merged block
.fill_model_block <- function(block, defaults) {
  if (is.null(block)) block <- list()

  # Top-level model keys (non-nested)
  nested_keys <- c("covariates", "random_effects", "slopes",
                   "brms", "dlnm", "kernel")
  for (nm in setdiff(names(defaults), nested_keys)) {
    if (is.null(block[[nm]])) block[[nm]] <- defaults[[nm]]
  }

  # covariates: one-level merge
  block$covariates <- .merge1(block$covariates, defaults$covariates %||% list())

  # random_effects: one-level merge (slopes sub-key stays as-is)
  if (is.null(block$random_effects)) {
    block$random_effects <- defaults$random_effects %||% list()
  } else {
    re_def <- defaults$random_effects %||% list()
    for (nm in names(re_def)) {
      if (nm == "slopes") next   # granular escape hatch — don't overwrite with defaults
      if (is.null(block$random_effects[[nm]])) block$random_effects[[nm]] <- re_def[[nm]]
    }
    if (is.null(block$random_effects$slopes)) block$random_effects$slopes <- list()
  }

  # dlnm sub-block: fill 'enable' + 'defaults' from canonical defaults.
  # 'variables' is config-specific and never overwritten.
  # Backwards-compat: if the old single-variable scalar key 'variable' is present
  # and 'variables' is absent, wrap it into the new list form.
  if (!is.null(defaults$dlnm)) {
    if (is.null(block$dlnm)) {
      block$dlnm <- defaults$dlnm
    } else {
      if (is.null(block$dlnm$enable))
        block$dlnm$enable <- defaults$dlnm$enable
      block$dlnm$defaults <- .merge1(block$dlnm$defaults,
                                     defaults$dlnm$defaults %||% list())
      if (is.null(block$dlnm$variables)) block$dlnm$variables <- list()
    }
    # Backwards-compat shim: if the old scalar key 'variable' (a character string)
    # is present and 'variables' is absent or empty, wrap it as a one-entry named
    # list.  Only fires when block$dlnm$variable is a non-empty character string.
    .old_dlnm_var <- block$dlnm[["variable"]]
    if (is.character(.old_dlnm_var) && length(.old_dlnm_var) == 1L && nzchar(.old_dlnm_var)) {
      if (!is.list(block$dlnm$variables) || length(block$dlnm$variables) == 0L)
        block$dlnm$variables <- stats::setNames(list(list()), .old_dlnm_var)
      block$dlnm[["variable"]] <- NULL
    }
    rm(.old_dlnm_var)
  }

  # kernel sub-block (occurrence brms_kl backend): nested (prior_grouping,
  # spatial), so filled recursively; absent defaults -> left untouched.
  if (!is.null(defaults$kernel))
    block$kernel <- .merge_recursive(block$kernel, defaults$kernel)

  # Backend sub-block: fill the active backend's keys from its defaults.
  # brms siblings (_cb crossbasis, _kl lag kernel) fall back to brms defaults
  # (same MCMC settings).
  be <- block$backend %||% defaults$backend %||% "brms_cb"
  be_defaults <- defaults[[be]] %||%
    if (be %in% c("brms_cb", "brms_kl")) defaults[["brms"]] else NULL
  if (grepl("_(cb|kl)$", be)) {
    # Inherit from the config file's own plain-backend sibling (if present)
    # before falling to static canonical defaults, so a value set only under
    # the plain sibling (e.g. re_uncorrelated) still reaches the active _cb
    # block. Precedence: block[[be]] (explicit) > block[[plain]] (sibling,
    # if present) > defaults[[be]] (canonical).
    plain <- sub("_(cb|kl)$", "", be)
    be_defaults <- .merge1(block[[plain]], be_defaults %||% list())
  }
  if (!is.null(be_defaults)) {
    block[[be]] <- .merge1(block[[be]], be_defaults)
  }

  block
}


# Outbreak-data root of the two defaults below that live there: the same rule
# as paths.R::SECURE_OUTBREAK_ROOT (CHOLERA_OUTBREAK_DATA_DIR, else
# ~/secure/cholera_outbreak_data), restated because modules never source each
# other. Unset, the defaults equal config_defaults.yml's literals.
.outbreak_data_root <- function() {
  d <- Sys.getenv("CHOLERA_OUTBREAK_DATA_DIR")
  if (nzchar(d)) sub("/+$", "", d) else "~/secure/cholera_outbreak_data"
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' get_default_options
#'
#' Returns the full set of default option values, mirroring config_defaults.yml.
#' Pure-CLI invocations (no --config) still get all keys filled.
#'
#' @return named list of defaults
get_default_options <- function() {
  list(
    # Top-level shared
    baseline          = "1990_2019",
    country           = "all",
    TL                = "1970-01-01",
    TR                = "2024-12-31",
    seed              = 1234L,
    outdir            = "generated_data/",
    overwrite         = FALSE,

    # ENSO data-prep (shared)
    enso = list(
      lag_months = 11L,
      source     = "oni",
      cache_dir  = "generated_data/"
    ),

    # Distance-weighted neighbour contagion (shared). Produced by
    # 03_prepare_occ.R via utils/spatial_contagion.R as the four columns
    # nbr_lag1_r{200,500}[_rn]; selected into the occurrence formula by
    # occurrence$covariates$nbr_contagion_lag1. Split top-level for exactly the
    # reason `enso` is: the prep script has no occurrence$ opt, the formula
    # builder does. `scale`: "raw" (default) | "rownorm" -> the _rn sibling;
    # `partition`: "none" (default) | "group" -> the _pg sibling.
    contagion = list(
      rho_km    = 500L,
      scale     = "raw",
      # "none" | "group": the partition-restricted _pg sibling (neighbours and
      # row-normalising denominator limited to the unit's grouping.source group)
      # | "registry": the _part<token> sibling, limited to the unit's group in
      # the registry `partition_source` (independent of grouping.source).
      partition = "none",
      partition_source = "",
      gadm_path = "data/geodata/gadm_valid.gpkg",
      refresh   = FALSE
    ),

    # Country-grouping (shared). When enable=TRUE, the run is restricted to the
    # countries of `group` in registry analysis/groupings/<source>.csv (cols
    # iso3,group) and random effects nest only country/ADM1 (region dropped).
    grouping = list(
      enable = FALSE,
      source = NULL,
      group  = NULL
    ),

    # Population crosswalk (shared). Drives the 02_*_population extraction
    # pipeline (LP x gid WorldPop mass H[l,g,t]) and the downstream accessors.
    # Every default here reproduces current behaviour: the artefacts are only
    # consumed where a caller opts in (size.pop_source, simulator.population).
    # Cross-track episode definition (shared). lp_season = historical
    # per-location-period episodes; gid_season = one canonical episode per
    # (gid, year_seas) shared by all four tracks. Selects the "_ep-gs" artefact
    # names (utils/paths.R) and the loaders' episode guard (assertEpisodeKey()).
    episodes = list(
      key           = "lp_season",
      # gid_season only: occurrence panel with national-zero month rows
      # (03_prepare_occ.R --national_zero_fill; "_nz" artefact token).
      national_zero_fill = FALSE,
      # national_zero_fill only: WHO country-years under this many cases also
      # count as national zeros (01 --zero_threshold); "_nz<N>" token. 0 = exact.
      national_zero_threshold = 0L,
      # national_zero_fill only: named rule of analysis/national_zero_rules.yml
      # (01 --rule; 03 --national_zero_rule); "" = legacy, else "_nz<N>-<rule>".
      national_zero_rule = "",
      # gid_season only: "season" = at risk from the season
      # start; "after_episode" = rows inside an earlier season's running episode
      # are masked (03 --at_risk_start; "_ae" token).
      at_risk_start = "season",
      # gid_season only: "drop" = the load-time ADM1-ADM2
      # overlap drop alone; "censor_orphan" = also censor the ADM2 children's
      # rows inside an orphan ADM1 parent episode (03 --adm1_overlap; "_a1c").
      adm1_overlap = "drop"
    ),

    population = list(
      crosswalk_path   = NULL,   # NULL -> makePopRunFile(opt, "pop_crosswalk")
      gid_pop_path     = NULL,   # NULL -> makePopRunFile(opt, "pop_gid")
      years            = c(2015L, 2030L),   # inclusive extraction range
      scaling          = "un_wpp2024",
      wpp_table        = "generated_data/wpp2024_national_totals.csv",
      raster_dir       = "data/worldpop",
      boundary_cache   = "data/country_boundaries",
      gadm_path        = "data/geodata/gadm_valid.gpkg",
      # The vintage directory must equal paths.R::defaultOutbreakDataDir()
      # (config.R may not source paths.R; tests/test_lp_shapefiles_vintage.R).
      lp_shapefiles    = file.path(.outbreak_data_root(), "outbreak_data_Nov_2026", "outbreak_shapefiles.rds"),
      country_registry = "analysis/groupings/enso_grouping.csv",
      sliver_mass_frac = 1.0e-4,
      sliver_mass_min  = 1.0,
      reconcile_tol    = 0.01,
      composite = list(
        cover_threshold = 0.90,
        c_min           = 0.05,
        # Which GADM tiling an LP is read at: the geometry-derived `lp_level`
        # frozen by 02c, or the outbreak corpus's own admin_level string (whose
        # admin depth is not GADM's).
        tiling          = "derived",  # derived | corpus_label
        # Area tie-break used by deriveLpTilingLevel() when both GADM levels
        # are population-admissible. "tolerant" is the default; "strict" is
        # kept side-by-side for reproducibility of earlier results.
        tiebreak        = "tolerant", # strict | tolerant
        tiebreak_floor  = 0.80,
        tiebreak_tol    = 0.10
      ),
      slurm = list(
        partition         = "public-cpu",
        time              = "04:00:00",
        mem_per_task      = "16G",
        cpus_per_task     = 1L,
        array_concurrency = 1000L,
        log_dir           = "logs"
      )
    ),

    # Occurrence model defaults
    occurrence = list(
      backend             = "brms_cb",
      admin_level         = "ADM2",
      ongoing_outbreak_re = TRUE,
      # Baseline-rate anchor (utils/anchor.R); see config_defaults.yml for the
      # full rationale. FALSE keeps every pre-existing config bit-inert.
      offset_baserate     = FALSE,
      # P3 for occurrence: re-key admin_level from the geometry-derived GADM
      # level before the admin-level filter (utils/population.R::.applyP3()).
      # NULL = corpus labels (current behaviour).
      canonical_admin_level = NULL,
      # time_to_last_outbreak + level "none" (no recorded prior onset of the
      # location), a load-time recode appended after "[6,Inf]"
      # (occ_data.R::makeTimeToLastOutbreakCat(no_prior_level =)); fit token "_ttn".
      ttlo_no_prior_level = FALSE,

      covariates = list(
        time_to_last_outbreak_cat = TRUE,
        ongoing_outbreak_lag1          = TRUE,
        ongoing_outbreak_country_lag1  = TRUE,
        ongoing_outbreak_tu_lag1       = TRUE,
        # Distance-weighted neighbour contagion (utils/spatial_contagion.R).
        # FALSE is load-bearing: makeOccRunFile()'s nbr_token is gated on this
        # flag, so every pre-existing config keeps its exact artefact names.
        nbr_contagion_lag1             = FALSE,
        enso                           = FALSE,
        months                    = FALSE,
        quarter                   = FALSE
      ),

      # Active-TRUE flags only. Under the resolve_slopes() regex decoder any
      # re_adm1_<col>_slope key set TRUE by a YAML config is automatically
      # wired — the full false-flag enumeration is kept as a comment in
      # config_defaults.yml for discoverability but is non-load-bearing here.
      random_effects = list(
        re_region_country_adm1_intercept   = TRUE,
        re_region_country_ttlo_slope       = TRUE,
        re_adm1_mean_precip_slope          = TRUE,
        re_adm1_dev_precip_slope           = TRUE,
        re_adm1_mean_temp_daymean_slope    = TRUE,
        re_adm1_dev_temp_daymean_slope     = TRUE,
        re_region_country_ongoing_slope         = TRUE,
        re_region_country_ongoing_country_slope = TRUE,
        re_region_country_ongoing_tu_slope      = TRUE,
        re_adm2_intercept                       = FALSE,
        re_region_country_basehaz_slope    = FALSE,
        # Neighbour-contagion slope REs. Shallow/deep pair, same split as the
        # ENSO slope flags; both gated on covariates$nbr_contagion_lag1.
        re_region_country_nbr_contagion_slope = FALSE,
        re_adm1_nbr_contagion_slope           = FALSE,
        slopes                             = list()
      ),

      # DLNM crossbasis defaults (consumed by 04b_run_occ_brms_cb.R)
      dlnm = list(
        enable   = FALSE,
        defaults = list(
          max_lag          = 3L,
          argvar           = list(fun = "bs", degree = 3L, df = 4L),
          arglag           = list(fun = "ns", df = 3L),
          keep_linear_term = FALSE,
          re_adm1_slope    = FALSE,
          single_lag       = NULL,
          lags             = NULL
        ),
        variables = list()
      ),

      # Gaussian lag-kernel ENSO term + latent classes (brms_kl backend;
      # utils/dlnm_kernel.R). Read only when occurrence.backend == "brms_kl".
      kernel = list(
        cen            = 0,
        sigma_lb       = 0.3,
        prior_mu       = NULL,
        prior_mu_sd    = 3,
        prior_sigma_sd = 2,
        prior_beta_sd  = 0.5,
        classes        = 1L,
        engine         = "stan",
        form           = "gaussian",
        lag_df         = 10L,
        lag_prior      = "horseshoe",
        hs             = list(scale_global = 0.1, slab_scale = 1, slab_df = 4),
        discrete       = list(widths = 0, centre_step = 1L, combo_prior = "exchangeable"),
        pin            = list(enable = FALSE, seeds_n = 1L, seeds = NULL, arrangements = "pilot",
                              cover = 0.95, max_runs = 6L, parallel_runs = 2L,
                              order_source = "pilot", sweep_batch = NULL,
                              sweep_outdir = "generated_data/",
                              pilot = list(chains = 2L, warmup = 200L, sampling = 200L)),
        re_slope       = list(enable = FALSE, country = TRUE, adm1 = TRUE, sd_country = 0.35, sd_adm1 = 0.35,
                              sd_country_fixed = NULL),
        init_mu        = NULL,
        ksd_fixed      = NULL,
        seeds          = list(per_group = 1L, iso3 = NULL),
        sweep          = list(batch = NULL, outdir = "generated_data/"),
        grainsize      = 1L,
        grainsize_rows = 1L,
        split_big      = TRUE,
        pair_collapse  = TRUE,
        identify       = "auto",
        repel          = list(enable = FALSE, lag_scale = 3, erf_scale = 0.5, erf_grid = c(-1.5, -0.75, 0.75, 1.5)),
        lag_re         = list(enable = FALSE, sd = 2, scale = "month"),
        erf_at         = 1.5,
        alpha          = 1,
        init           = "auto",
        init_starts    = 8L,
        drop_region    = NULL,
        gq_row_loglik  = FALSE,
        keep_re_draws  = FALSE,
        partition_name = NULL,
        prior_grouping = list(
          source         = NULL,
          kappa          = NULL,
          kappa_prior_sd = 2,
          deploy_config  = NULL,
          membership     = NULL,
          contiguity     = list(enable = FALSE, weight = 0.5)
        ),
        spatial = list(
          enable       = FALSE,
          prior        = "icar",
          adjacency    = "country_adjacency",
          tau_prior_sd = 1,
          rho_prior    = c(2, 2)
        )
      ),

      brms = list(
        basehaz_poisson = "seasonal",
        sampling        = 1000L,
        warmup          = 1000L,
        chains          = 4L,
        cores           = 4L,
        refresh         = NULL,
        # Matches size/duration's brms sampler settings (brms default is
        # adapt_delta 0.8 / max_treedepth 10).
        adapt_delta     = 0.95,
        max_treedepth   = 12L,
        init            = "random",
        re_init_r       = NULL,
        re_uncorrelated = FALSE
      )
    ),

    # Outbreak-size model defaults
    size = list(
      backend         = "brms",
      family          = "negbinomial",
      link            = "log",
      admin_level     = "ADM1",
      offset_pop      = TRUE,
      offset_meanrate = TRUE,
      # Size's anchor offset predates this key (offset_meanrate, on by default),
      # so here offset_baserate is purely the *prior* gate. FALSE = the 322
      # configs pinning prior_intercept keep fitting the model they fit before;
      # flipping it is a modelling decision, not a refactor side-effect.
      offset_baserate = FALSE,
      seasonal_term   = "none",

      min_pop         = 10L,
      max_attack_rate = 0.10,
      # `pooled`: one estimator across all four
      # tracks (utils/anchor.R `.ANCHOR_SPEC`). Set "geometric" here to refit an
      # archived size fit like-for-like; the sidecar path carries est-<estimator>
      # so the two tables can never collide on disk.
      mean_rate       = list(estimator          = "pooled",
                             level_factor       = TRUE,
                             cross_absent_levels = TRUE),

      # Which population the size estimand is denominated in, on BOTH sides of
      # mu = pop x mean_rate x exp(eta). `denominator` and `rate_denominator`
      # must agree: they are the same quantity, and resolveSizeEstimand() hard-
      # errors on a mismatch precisely because the likelihood cannot detect one
      # (both offsets are fixed, so a units error is pure multiplicative bias).
      estimand = list(
        denominator               = "footprint",   # gid | footprint
        rate_denominator          = "footprint",   # must equal `denominator`
        underreporting_multiplier = 1.0
      ),
      # `fixed` imposes beta == 1 via offset(log(pop)); `free`/`free_by_level`
      # add the centred deviation term log_pop_c so the coefficient reads as
      # beta - 1 and beta == 1 recovers the offset-only fit exactly.
      pop_elasticity  = "fixed",
      # NULL = keep the corpus `admin_level` labels (current behaviour); TRUE =
      # re-key from the canonical `gadm_level` in the gid population table.
      canonical_admin_level = NULL,

      # Overdispersion RE: random intercepts on phi (log link). constant = scalar
      # phi (current default); other levels add (1 | group).
      phi_re = "constant",   # constant | region | region_country | region_country_adm1

      covariates = list(),

      random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm1_ttlo_slope                = TRUE,
        re_adm1_mean_precip_slope         = TRUE,
        re_adm1_dev_precip_slope          = TRUE,
        re_adm1_mean_temp_daymean_slope   = TRUE,
        re_adm1_dev_temp_daymean_slope    = TRUE,
        re_region_country_ongoing_slope   = FALSE,
        re_adm2_intercept                 = FALSE,
        # all = historical leaf over every unit; adm2_only = leaf on canonical
        # ADM2 units only (see .intercept_re_term()).
        re_adm2_leaf_scope                = "all",
        slopes                            = list()
      ),

      brms = list(
        sampling            = 1000L,
        warmup              = 1000L,
        chains              = 4L,
        cores               = 4L,
        adapt_delta         = 0.95,
        max_treedepth       = 12L,
        init                = 0.1,
        prior_b             = "normal(0, 0.5)",
        prior_intercept     = "normal(-3, 1)",
        prior_sd            = "normal(0, 0.5)",
        prior_phi_intercept = "normal(1, 1)",
        prior_phi_sd        = "normal(0, 0.5)"
      ),

      # DLNM crossbasis defaults for the size model (consumed by
      # 04d_run_siz_brms_cb.R).  Schema mirrors occurrence.dlnm.
      dlnm = list(
        enable   = FALSE,
        defaults = list(
          max_lag          = 3L,
          argvar           = list(fun = "bs", degree = 3L, df = 4L),
          arglag           = list(fun = "ns", df = 3L),
          keep_linear_term = FALSE,
          re_adm1_slope    = FALSE,
          single_lag       = NULL
        ),
        variables = list()
      )
    ),

    # Outbreak-duration model defaults -- the discrete-time cessation hazard,
    # fitted as Poisson-with-log-offset on one row per ACTIVE MONTH per episode
    # (analysis/utils/dur_data.R). Third model block alongside occurrence/size;
    # merged by .fill_model_block(), NOT by .CONFIG_TOPLEVEL_BLOCKS.
    duration = list(
      backend     = "brms",
      family      = "poisson",
      link        = "log",
      # ADM1 | ADM2 = intercept-RE depth only, no row filter (historical);
      # ADM1-ADM2 = the monthly episode set (dur_data.R::durationEpisodeRows()).
      admin_level = "ADM1",
      # Baseline-rate anchor; denominator is pmax(exposure, 1), exactly
      # exp(log_offset), the track's only other offset. See config_defaults.yml.
      offset_baserate = FALSE,
      # P3 for duration (see occurrence.canonical_admin_level). NULL = off.
      canonical_admin_level = NULL,

      # PRIMARY baseline hazard, on months-since-onset (mso). This is the
      # cessation hazard's natural clock: a 1-month flare and a 9-month
      # entrenched epidemic do not stop at the same rate. Measured on the
      # prepared panel the monthly hazard runs 0.0006 -> 0.31 across bands
      # (see .DUR_MSO_BREAKS in dur_data.R), so "none" is a sensitivity
      # variant, never a default.
      basehaz_duration = "factor",   # factor | ms | none
      # SECONDARY baseline hazard, on seasonal time (t_start_seas) -- the
      # "outbreaks end when the season ends" hypothesis. Reuses the occurrence
      # branch verbatim. Off by default: mso already carries the time axis the
      # likelihood needs, and stacking both is a modelling choice, not a default.
      basehaz_poisson  = "none",     # none | factor | cyclic | seasonal | ms
      # Calendar-month seasonality of the cessation hazard, built exactly as
      # monthly.seasonal_term "harmonic2": sin/cos of the active month's
      # calendar month at the annual and semi-annual period, plain class-b
      # fixed effects (load_duration_data() builds the columns). Off by
      # default; the deployed duration fits carry no seasonal term.
      seasonal_term    = "none",     # none | harmonic2
      # mso band edges for bin_mso(). MUST end in Inf. Rebuilt at load time by
      # load_duration_data(), so a band sensitivity run does not require
      # re-running 03_prepare_dur.R.
      mso_breaks       = c(1, 2, 3, 4, 6, 9, Inf),
      # Which event resets time_to_last_outbreak in the v1 kernel. `onset`
      # reproduces 03_prepare_occ.R:650 (the occurrence track's fit-time rule);
      # `end` matches the size track's TR-based reset (03_prepare_siz.R:768-780).
      # The two tracks genuinely disagree; default `onset` changes nothing.
      ttlo_reset       = "onset",    # onset | end

      # Standalone toggles (does not inherit from occurrence or size). Empty
      # {} => only the slope-implied covariates reach the formula (convention 2).
      covariates = list(),

      random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept                = FALSE, # needs admin_level ADM1-ADM2 (or ADM2)
        re_adm2_leaf_scope               = "all", # all | adm2_only
        re_adm1_mean_precip_slope        = TRUE,
        re_adm1_dev_precip_slope         = TRUE,
        re_adm1_mean_temp_daymean_slope  = TRUE,
        re_adm1_dev_temp_daymean_slope   = TRUE,
        slopes                           = list()
      ),

      # DLNM crossbasis defaults for the duration model (consumed by
      # 04_fit_dur_brms_cb.R). Schema mirrors occurrence.dlnm / size.dlnm.
      dlnm = list(
        enable   = FALSE,
        defaults = list(
          max_lag          = 3L,
          argvar           = list(fun = "bs", degree = 3L, df = 4L),
          arglag           = list(fun = "ns", df = 3L),
          keep_linear_term = FALSE,
          re_adm1_slope    = FALSE,
          single_lag       = NULL
        ),
        variables = list()
      ),

      brms = list(
        sampling        = 1000L,
        warmup          = 1000L,
        chains          = 4L,
        cores           = 4L,
        adapt_delta     = 0.95,
        max_treedepth   = 12L,
        init            = "random",
        re_uncorrelated = FALSE,
        prior_b         = "normal(0, 0.5)",
        prior_intercept = "normal(0, 1.5)",
        prior_sd        = "normal(0, 2.5)",
        # Own prior on the episode-age band effects (mso_band2 .. mso_band<last>),
        # overriding prior_b for those coefficients only. NULL = they share
        # prior_b. Fit-name token "_mbp<sd>".
        prior_mso_band  = NULL
      )
    ),

    # Outbreak monthly-size model defaults -- fourth model track alongside
    # occurrence/size/duration, response cases_month (per-calendar-month
    # outbreak size). Bayesian hierarchical version of a validated negative-
    # binomial GAM. merged by .fill_model_block(), NOT by
    # .CONFIG_TOPLEVEL_BLOCKS.
    monthly = list(
      backend     = "brms_cb",       # brms | brms_cb (brms_cb is the primary model)
      family      = "negbinomial",
      link        = "log",
      admin_level = "ADM1-ADM2",     # ADM1 | ADM2 | ADM1-ADM2 (slope cap + row filter)
      offset_pop  = TRUE,            # offset(log(pop))
      # offset(log_offset) is UNCONDITIONAL (duration's behaviour): log_offset =
      # log(pmax(exposure, 1)) is already on the panel (dur_data.R).
      # Baseline-rate anchor; mon is the one anchored track with TWO other
      # offsets, so its denominator is pop * pmax(exposure, 1). See
      # config_defaults.yml.
      offset_baserate = FALSE,

      # Which population the monthly estimand is denominated in. Same apparatus
      # as size.estimand -- resolved by the shared .resolveTrackEstimand()
      # (utils/provenance.R) so the two tracks cannot drift. `rate_denominator`
      # is the population base_rate_mon was built on (mon's analogue of size's
      # mean_rate); it must agree with `denominator` for the same reason it must
      # on size -- both offsets are fixed, so a units error is pure bias.
      estimand = list(
        denominator               = "footprint",   # gid | footprint
        rate_denominator          = "footprint",   # must equal `denominator`
        underreporting_multiplier = 1.0
      ),
      # `fixed` imposes beta == 1 via offset(log(pop)); `free`/`free_by_level`
      # add the centred deviation term log_pop_c.
      pop_elasticity  = "fixed",
      # NULL = keep the corpus `admin_level` labels (current behaviour); TRUE =
      # re-key from the canonical `gadm_level` in the gid population table.
      canonical_admin_level = NULL,

      seasonal_term = "harmonic2",   # none | factor_month | cyclic | harmonic2
      phi_re        = "region_country_adm1", # constant | region | region_country | region_country_adm1

      # Flat data-shaping key, same style as duration.mso_breaks/ttlo_reset.
      # TRUE => mon_data.R builds cases_prev/log_cases_prev.
      prev_cases = FALSE,

      # Standalone, does not inherit from size/occurrence. Empty {} => only the
      # slope-implied covariates reach the formula (convention 2).
      covariates = list(),

      random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept                = TRUE, # leaf col = gid => (1 | .../ADM1/gid)
        re_adm2_leaf_scope               = "all", # all | adm2_only
        slopes                           = list()
      ),

      # DLNM crossbasis defaults for the monthly model. Schema mirrors
      # occurrence.dlnm / size.dlnm / duration.dlnm, plus the monthly-only
      # `re_slope_levels` key below.
      dlnm = list(
        enable   = TRUE,
        defaults = list(
          max_lag          = 3L,
          argvar           = list(fun = "bs", degree = 3L, df = 4L),
          arglag           = list(fun = "ns", df = 3L),
          keep_linear_term = FALSE,
          re_adm1_slope    = FALSE,
          # Which re_group_by_level() levels nested_slope_re_terms() emits for a
          # CB variable's random slopes. Monthly-only, and read by
          # 04_fit_mon_brms_cb.R alone -- 04_fit_mon_brms.R never touches the
          # dlnm: block, and the occ/siz/dur _cb scripts keep the three-level
          # default hardcoded in nested_slope_re_terms().
          re_slope_levels  = c("region", "country", "adm1"),
          single_lag       = 0L  # 0 = onebasis on the bare column (mso, month_num, ...)
        ),
        variables = list(
          # mso drops the country:ADM1 slope level: 4 cb columns x 2 nesting
          # levels are not identified on monthly per-ADM1 data, and tail draws of
          # the ADM1 slopes blow up the simulator's per-unit mu.
          mso  = list(single_lag = 0L, argvar = list(fun = "bs", degree = 3L, df = 4L), re_adm1_slope = TRUE,
                      re_slope_levels = c("region", "country")),
          enso = list(single_lag = 4L, argvar = list(fun = "bs", degree = 3L, df = 8L), re_adm1_slope = TRUE)
        )
      ),

      # --- brms backend sub-block --- (byte-identical key set to size.brms).
      # No brms_cb: sub-block -- size doesn't have one either; the fit script
      # falls back with opt$monthly$brms_cb %||% opt$monthly$brms %||% list().
      brms = list(
        sampling            = 1000L,
        warmup              = 1000L,
        chains              = 4L,
        cores               = 4L,
        adapt_delta         = 0.95,
        max_treedepth       = 12L,
        init                = 0.1,
        prior_b             = "normal(0, 0.5)",
        prior_intercept     = "normal(-3, 1)",
        prior_sd            = "normal(0, 0.5)",
        prior_phi_intercept = "normal(1, 1)",
        prior_phi_sd        = "normal(0, 0.5)"
      )
    ),

    # Simulator defaults (consumed by 07_simulator_run.R dispatcher + offspring;
    # mirrors config_defaults.yml)
    simulator = list(
      backend             = "v1",
      run_tag             = NULL,
      occurrence_fit      = NULL,
      # Monthly case-count (`mon`) track: a case count at every active month.
      # Requires simulator.duration_fit. The size track is retired downstream:
      # size_fit / size_cb_meta / size.* / location_frame.mean_rate_table_path
      # are refused by .validate_basehaz().
      monthly_fit         = NULL,
      monthly_cb_meta     = NULL,
      climate_panel       = NULL,
      # Grouped mode: load one occ (+ dur/mon) fit per group from `fits_dir`, route
      # each country to its group's draws via registry <source>.csv.
      grouping            = list(
        enable         = FALSE,
        source         = NULL,
        fits_dir       = NULL,
        fits_baseline  = NULL,
        monthly_enable = FALSE
      ),
      # Filter (spatPomp bpfilter) defaults. Each key's value is the EFFECTIVE
      # runtime default, i.e. the consumer's `%||%` fallback, so filling it here
      # changes no run (the runtime R defaults
      # win; config_defaults.yml restates these values and
      # tests/test_config_defaults_parity.R + test_hardcoded_default_parity.R
      # gate both). Where two consumers disagreed, the value is the one the
      # deployed path uses and the other literal stays, commented, at its site:
      #   enable        FALSE: the only live reader is simulations.R
      #                 (isTRUE(sim$filter$enable) -> survival panel for the
      #                 warmup obs); 07b_simulator_spatpomp.R's
      #                 `flt$enable %||% TRUE` feeds a variable nothing reads.
      #   warmup_window 07b_simulator_spatpomp.R's fallback; every forecast /
      #                 hindcast / event config omits it and runs on it, while
      #                 every smoothing (historical) config sets it, so
      #                 smoothing_core.R's `%||% list()` (= keep its own
      #                 2008-01 burn-in) is reached only by a hand-built opt.
      #   forecast_window / smoothing_input: NULL = unset. 07b derives the
      #                 forecast window from warmup_window$end and
      #                 horizon$end, and smoothing keeps its own FILTER_START /
      #                 FILTER_END, so no static value reproduces both.
      # Np is read by 07b only: 08 smoothing reads it from the raw, unmerged
      # config (default 25000, smoothing_core.R::smoothing_default_config()).
      # save_states / forecast_engine are not read; configs that carry them
      # are unaffected (no key check).
      filter              = list(
        enable            = FALSE,
        Np                = 2000L,
        block_granularity = "country",
        partition_by      = "transmission_unit",
        partition_source  = "",
        warmup_window     = list(start = "2021-07-01", end = "2023-05-01"),
        forecast_window   = NULL,
        measurement       = "occurrence_duration_mon_NB",
        n_draws_in_filter = 1L,
        # Smoothing draw design: "per_tu_random" = per-unit
        # random draws; "shared_cycle" = replicate j
        # runs draw j in every unit (see config_defaults.yml).
        draw_design       = "per_tu_random",
        smoothing_input   = NULL,
        # delta_fn/delta_fp are inert for every non-`*_state_detect`
        # measurement. covariate_feed / monthly_case_obs: "latent" / "none" ARE
        # the baseline (resolve_filter_cfg(), spatpomp_core.R).
        # `who_obs` disabled is the baseline, and enabling it is
        # refused unless the `_who` measurement is selected too.
        delta_fn         = 0.10,
        delta_fp         = 0.05,
        # Inactive-state case mixture (eps, mu_0) of the `_actn` measurements
        # (resolve_filter_cfg()); inert for every other measurement.
        case_fp_rate     = 1.0,
        case_fp_mean     = 0.1,
        covariate_feed   = "latent",
        monthly_case_obs = "none",
        # filter_traj storage of run_one_bpfilter_replicate() (spatpomp_core.R):
        # "pruned" (the default) / "pruned_float32" =
        # .bpfilter_traj_compact(), pruned every traj_prune_every steps;
        # "full" = spatPomp's own bpfilter. See config_defaults.yml.
        traj_storage     = "pruned",
        traj_prune_every = 1L,
        who_obs          = list(
          enable            = FALSE,
          path              = "generated_data/who_dashboard_monthly_africa.rds",
          start             = "2023-01-01",
          end               = "2025-01-01",
          phi               = 5,
          mu_floor          = 1,
          fp_rate           = 1.0,
          scale             = 1,
          block_share       = "pop",
          allow_adm1_blocks = FALSE,
          # Per-country scale: "fixed" (scale / who_f)
          # or "estimated" (who_dashboard.R::
          # whoScaleEstimate() over scale_est); db_overlap "score_both"
          # or "skip_who". See config_defaults.yml.
          scale_mode        = "fixed",
          scale_est         = list(start = "2023-01-01", end = "2025-01-01",
                                   min_db_cases = 100L, min_months = 6L,
                                   shrink_cases = 0, fallback = "pooled",
                                   clamp = c(0.1, 10)),
          db_overlap        = "score_both"
        )
      ),
      location_frame      = list(
        hmm_clusters_path      = "data/hmm_country_clusters.csv",
        outbreak_metadata_path = file.path(.outbreak_data_root(), "outbreak_metadata.rds"),
        size_data_path         = "generated_data/outbreak_size_analysis_data.rds",
        # Year at which the simulator / smoothing population Q_g is resolved:
        # an integer year, or `season_year` / `calendar_year` (a per-month
        # population). Read only through simulator_core.R::simPopYearSpec()
        # 2020 = one fixed population year.
        pop_year               = 2020L
      ),
      horizon             = list(start = "2024-01-01", end = "2034-12-31"),
      n_replicates        = 100L,
      draws_per_replicate = 1L,
      posterior_strategy  = "paired",
      seed                = NULL,
      scope = list(
        transmission_units = "all",
        countries          = "all",
        admin_level        = "ADM2",
        unseen_adm1        = "restrict"
      ),
      # Read by no code. base_year != 2020 is refused by simPopYearSpec(),
      # which points to location_frame.pop_year.
      population = list(
        source    = "worldpop",
        base_year = 2020L,
        scaling   = "un_wpp2024"
      ),
      mean_rate = "frozen",
      # Suppress o_ct/o_tu state tracking AND LP contribution in every kernel,
      # independent of the fit. Default FALSE = today's behaviour unchanged.
      # See config_defaults.yml for the full rationale.
      drop_country_tu_spillover = FALSE,
      # Monthly cases per unit <= incidence_cap x pop, on the mean and the drawn
      # count, in every kernel and the smoothing; .inf disables. Read and
      # validated by simulator_core.R::incidenceCapOf(). See config_defaults.yml.
      incidence_cap = 0.15,
      # Which inactive units may draw an onset: "inactive" (every unit not in a
      # running episode) or "season" (also no recorded onset yet this season,
      # matching the occurrence fit's risk set). Read and
      # validated by simulator_core.R::onsetEligibilityOf(). See config_defaults.yml.
      onset_eligibility = "inactive",
      # Monthly case exposure of a unit's cessation month: "full" (a full month)
      # or "half" (the cessation is drawn before
      # the monthly mean and a ceasing month gets half a month, as in the fit
      # panel). Read and validated by
      # simulator_core.R::cessationExposureOf(). See config_defaults.yml.
      cessation_exposure = "full",
      # Distance-weighted neighbour contagion (nbr_lag1_r<rho>) consumed in the
      # v1 and spatPomp kernels. enable=FALSE (default) OR a fit that carries no
      # nbr_* fixed effect ⇒ zero contagion machinery built/executed, byte-
      # identical to pre-change. See config_defaults.yml for the full rationale.
      nbr_contagion = list(
        enable    = FALSE,
        rho_km    = 500L,
        gadm_path = "data/geodata/gadm_valid.gpkg",
        cache_dir = "generated_data/"
      ),
      # Forecast-phase posterior-draw handling (07b forecast handoff only).
      # param_draw_mode: "random" (default) = fresh chunk-safe random draw per
      # replicate, shared across a region's members; "cycle" = deterministic
      # 1..S spread; "inherited" = reuse handoff particle's draw (collapses to
      # one draw under n_draws_in_filter=1); "joint" = the
      # draw the handoff smoothing replicate ran (shared_cycle cache).
      # spread_param_draws is a retained
      # legacy alias, honoured only when param_draw_mode is blank/unset. See
      # config_defaults.yml for the full rationale.
      forecast = list(
        param_draw_mode = "random",
        spread_param_draws = FALSE,
        # Handoff state selection among the smoothing replicates: "uniform"
        # (default; balanced, by global replicate id),
        # "shuffled" (balanced per block of m ids, seeded permutation per
        # block) or "loglik" (SIR by exp(logLik), ESS ~1).
        # spatpomp_core.R::select_handoff_warmup(); see config_defaults.yml.
        handoff_weighting = "uniform",
        # N4 diagnostic (hindcast-only): anchor the forecast handoff's active/mso
        # to the OBSERVED active footprint at the launch month, bypassing the PF
        # cold-start under-recruitment. FALSE = today's behaviour (smoothed IC).
        # See config_defaults.yml for the full rationale.
        observed_active_handoff = FALSE
      ),
      # Validation-phase (07c) scoring options.
      # total_case_basis: "exact_gid" (today's behaviour) scores total cases on
      # an exact-(gid) observed↔frame join; "region_key" reconciles the observed
      # corpus to the sim's forecast footprint at the ADM1 region level (rule B),
      # recovering coarse-reported observed cases the exact-gid join drops. See
      # config_defaults.yml for the full rationale.
      validation = list(
        total_case_basis = "exact_gid"
      ),
      scenario = list(
        climate = "factual",
        enso    = "historical"
      ),
      output = list(
        format               = "qs",
        aggregate            = "adm2_month",
        keep_event_table     = TRUE,
        keep_monthly_summary = TRUE
      ),
      overwrite = FALSE,
      parallel = list(
        backend    = "furrr",
        over       = "transmission_unit",
        chunk_size = 25L,
        workers    = 1L
      ),
      # Per-replicate ENSO forecast source for the factual branch (analysis/utils/
      # simulations.R). `cfsv2_target_peak` is NULL = no target-peak correction,
      # the loader's fallback (`%||% NULL`); 2.5 is a recommended per-config
      # override (config_defaults.yml's simulator.enso.cfsv2_target_peak
      # comment), never a default.
      enso = list(
        forecast_source     = "cfsv2",
        cfsv2_target_peak   = NULL,
        cache_dir           = "generated_data/",
        cfsv2_url           = "https://www.cpc.ncep.noaa.gov/products/CFSv2/dataInd3/nino34Mon.nc",
        cfsv2_init          = "latest",
        cfsv2_rebase        = "seam",
        cfsv2_sd_correction = FALSE,
        cfsv2_calibration   = "",
        cfsv2_calibration_target = "ersst",
        clamp_to_knots      = FALSE,
        clamp_margin        = 0.25,
        # ENSO past the CFSv2 span (enso.R::sampleEnsoCFSv2Trajectories()).
        continuation        = list(method = "neutral", phi = NULL,
                                   fit_window = "1991_2020", legacy_seam = FALSE,
                                   library = "generated_data/enso_continuation_library.qs",
                                   k = 30L, blend_months = 3, min_sep = 24L,
                                   max_per_model = Inf, exclude_models = c("CM2.1-NMME"),
                                   # analog seam: last_cfsv2 | last_full_oni
                                   seam = "last_cfsv2")
      ),
      # SLURM array-job defaults, consumed by the v1 simulator's bash
      # stage and payload (not read
      # from `opt` in R — the bash scripts hardcode the same values directly).
      # Documented here so config_defaults.yml and the actual submission
      # scripts cannot drift without this parity gate catching it.
      slurm = list(
        partition         = "public-cpu",
        time              = "12:00:00",
        mem_per_task      = "8G",
        cpus_per_task     = 1L,
        modules           = c("GCC/11.3.0", "OpenMPI/4.1.4", "R/4.2.1"),
        array_concurrency = 1000L,
        log_dir           = "logs/"
      )
    ),

    # Manuscript results defaults (phases 09 postprocess + 10 figures).
    #
    # This block is the authoritative record of the *output decisions* behind
    # the manuscript figures: which ENSO analog set is the main scenario, which
    # groups are dropped, the attribution thresholds, the IRR subsample size,
    # and the Figure-1 lag-selection rule, so a figure can be traced back to
    # the decisions that produced it.
    #
    # `redo_cache` is deliberately NOT here: it is an operational override, not
    # an output decision, and stays a CLI-only flag on 09 (landing at
    # opt$redo_cache via the CLI-only merge below).
    results = list(
      # Optional pointer to a deploy config (analysis/configs/deploy/<stem>.yml).
      # "" (default) = the results block is fully self-contained, today's
      # behaviour unchanged. (Empty-string, not NULL, is the "unset" sentinel:
      # .merge_recursive cannot round-trip a NULL leaf -- `x[[nm]] <- NULL`
      # deletes the key -- so a NULL default could never be restated in a config
      # and still satisfy the self-containment coverage test. This mirrors the
      # deploy block's own launch_years="" inert sentinel.) When set to a path,
      # the overlapping `inputs.*` leaves are DERIVED from the pointed-to deploy
      # block by resolve_results_inputs_from_deploy() with explicit-over-
      # inherited precedence (any key physically present in the results YAML
      # still wins), and the deploy config path is added to the freshness
      # `inputs` vector so a change to it invalidates the 09 artefacts. See
      # utils/deploy_core.R.
      deploy_config = "",
      # Names of results.inputs leaves (manifest, grouping, sim_root,
      # smoothing_parts_dir, lag_loo_batch, best_lag_dir, hindcast_manifest, episode_key,
      # episode_national_zero_fill) that MAY be set
      # explicitly to a value differing from what `deploy_config` derives. With
      # deploy_config set, any OTHER explicit leaf that disagrees is refused by
      # resolve_results_inputs_from_deploy() (utils/deploy_core.R): a mixed-run
      # results block must be a stated decision, never a leftover. Empty = none.
      # Not a stamped block (results_core.R::.RESULTS_STAMPED_BLOCKS): it is a
      # validation policy, not an output decision.
      allow_overrides = character(0),
      # Forked worker processes for 09's per-group fit loads and per-run
      # hindcast scoring; artefacts are identical for any value. Not a stamped
      # block (speed and memory only). 09's --workers overrides it.
      workers = 1L,
      # brms_kl class-sweep comparison dir written by 05_fit_postprocess_kl_sweep.R
      # (generated_data/kl_sweep_compare/<batch>); "" = the report_kl_sweep step is
      # off. Not a stamped block: only report_kl_sweep reads it, and its files are
      # that artefact's freshness inputs. Not deploy-derived.
      kl_sweep_dir = "",
      # The sweep cell (fit id, e.g. "ns5_K5") the deployed grouping was derived
      # from; with kl_sweep_dir it turns on report_membership (membership
      # uncertainty of the deployed classes). Top level and not
      # stamped, like kl_sweep_dir; aligned from deploy.kl_sweep.cell by
      # --stage align_results.
      kl_sweep_cell = "",
      # Past-ENSO-event attribution manifest written by deploy --stage gen_events
      # (deploy.event_attribution.manifest); "" = the report_event_attrib step is
      # off. Not a stamped block: only that step reads it, and the manifest and
      # its run configs are the artefact's freshness inputs. Aligned from the
      # deploy by --stage align_results when the deploy sets events.
      event_manifest = "",
      # Posterior interval of the attributable fraction / number (09 steps
      # af_posterior, event_af_posterior): needs fixed-draw runs (param_draw_mode
      # cycle or joint); off for runs with random draws. Not a stamped block: it only
      # selects the two steps, so turning it on leaves every other artefact fresh.
      af_posterior = FALSE,
      # Replicates each (unit, posterior draw) needs after the explosive drop in
      # those steps (results_derive.R AF_POST_MIN_REPS). Unstamped.
      af_posterior_min_reps = 20L,
      # Which simulated series the manuscript tables report: "latent" (default,
      # the pre-existing tables) or "detected" (what surveillance would have
      # recorded; needs runs simulated with simulator.filter.covariate_feed:
      # detected). See results_core.R::results_select_series().
      series = "latent",
      # Raw inputs read by 09 only. 10 never touches a fit or a filter.
      inputs = list(
        manifest            = "analysis/configs/forecast_horizons.yml",
        grouping            = "analysis/groupings/enso_grouping.csv",
        sim_root            = "generated_data/sim_spatpomp",
        # MUST stay in step with simulator.filter.smoothing_input in
        # analysis/configs/forecast_country_by_group/*.yml -- Figure 3 splices
        # this smoothing arm onto the forecast arm built from those configs, so
        # a stale value silently glues two different filter runs together.
        # Run-specific: "" (unset sentinel) by default; the results YAML sets it
        # or derives it from `deploy_config`, and 09 refuses an unset value.
        smoothing_parts_dir = "",
        # Per-track batch prefix; the track stub is interpolated ->
        # analysis/configs/<batch>_<stub>_grouped_enso_singlelag/. `default`
        # applies to any track without its own entry (09b_postprocess_results.R's
        # LAG_LOO_DIR()). Each track's batch MUST be the one its best_lag_dir
        # archive entry was selected from. Run-specific: "" by default, as above.
        lag_loo_batch       = list(default = "", occurrence = "", duration = "",
                                   monthly = ""),
        best_lag_dir        = "generated_data/enso_best_lag_fits",
        # Attribution manifest (deploy-derived from
        # deploy.outputs.forecast_attrib_manifest; make_forecast_attrib_configs.R):
        # the factual of every attribution step. "" = the forecast manifest's
        # (stitched) factual. See config_defaults.yml.
        attrib_manifest     = "",
        # Skill of the analog ENSO continuation by lead (CRPSS vs climatology;
        # 03_prepare_enso_continuation_library.R --validate), read for table S28
        # (the CRPSS at each group's longest continuation lead). "" = no CRPSS
        # column. Stamped (an input path); not deploy-derived.
        enso_continuation_skill = "",
        # Curated literature table T1a (ENSO-rainfall link per group): the 09
        # rain_enso step copies it into its artefact (`mechanisms`); 11 prints
        # it. "" = no T1a. Stamped (an input path); not deploy-derived.
        enso_rain_mechanisms = "analysis/tables/enso_rain_mechanisms.yml",
        # Hindcast manifest (make_hindcast_configs.R output; deploy-derived from
        # deploy.outputs.hindcast_manifest). 09 scores every run it lists from
        # the run's validation/metrics + merged/monthly_summary
        # (utils/hindcast_skill.R); 07c's hindcast_summary/ is not read.
        hindcast_manifest   = "analysis/configs/hindcast_horizons.yml",
        # Cross-track episode definition every 09/10/11 panel is read under
        # (lp_season | gid_season). Deploy-derived from deploy.episodes.key and
        # copied onto the top-level episodes.key, so the survival / monthly /
        # duration / size panels match the episodes the deployed fits were
        # estimated on. A results input, hence stamped: changing it re-derives.
        episode_key         = "lp_season",
        # national zero-fill occurrence panel ("_nz", gid_season only). Deploy-derived
        # from deploy.episodes.national_zero_fill and copied onto the top-level
        # episodes.national_zero_fill alongside episode_key. Stamped, like episode_key.
        episode_national_zero_fill = FALSE,
        # its threshold ("_nz<N>"), deploy-derived from
        # deploy.episodes.national_zero_threshold. Stamped, like episode_key.
        episode_national_zero_threshold = 0L,
        # its rule ("_nz<N>-<rule>"; "" = legacy), deploy-derived from
        # deploy.episodes.national_zero_rule. Stamped (back-compat "").
        episode_national_zero_rule = "",
        # the panel's at-risk start and ADM1-overlap treatment
        # ("_ae" / "_a1c"), deploy-derived from deploy.episodes.*. Stamped
        # (back-compat "season" / "drop").
        episode_at_risk_start = "season",
        episode_adm1_overlap = "drop"
      ),
      # Read by BOTH scripts: selects what is derived AND what is plotted.
      scope = list(
        main_set       = "neutral",
        supp_sets      = c("all", "nina", "elnino"),
        supplementary  = FALSE,
        # Fit groups dropped from Figures 3-4; Figures 1-2 always show every
        # group. Every name must be a group of inputs.grouping (09 refuses an
        # unknown name). Empty: enso_central
        # (numeric instability) exists only in the legacy enso_grouping registry.
        exclude_groups = character(0),
        fig1_tracks    = c("occurrence", "duration", "monthly"),
        build_figures  = c("1", "2", "3", "4", "H"),   # 10 only; "H" = hindcast skill (S5 / S6 / FC / SC / AE: 10 --figures)
        # Hindcast skill (09 derives, 10 figure H, notebooks/hindcast_skill.qmd).
        # The headline variants are the fair ones (leave-one-window-out per-gid
        # climatology; surveillance-defined case support); the legacy in-sample
        # / obs-positive variants are always derived alongside.
        hindcast = list(
          enable                 = TRUE,
          headline_occ_reference = "loo_gid",
          headline_case_support  = "atrisk_gids",
          # reference of the monthly / season case-skill headline columns
          headline_case_skill_reference = "clim_pooled"
        ),
        # v2 (partition_by: registry) manifests only: a unit is simulated to its
        # longest member horizon. "drop" scores each country only up to its own
        # fit group's firm horizon; "flag" keeps the later months, flagged
        # beyond_firm_horizon. No effect on v1 manifests.
        beyond_firm_horizon = "drop",
        # One Africa window for every attribution step: with
        # window_end a first-of-month date, Africa is every included group over
        # [attribution start, window_end], a group whose firm horizon ends earlier
        # being read on into the analog ENSO continuation up to
        # continuation_max_lead months past it (0 = firm only; "to_end"), and,
        # with max_counterfactual_winter ("YYYY/YY"), no further than its last
        # month whose ONI reads stay before the next DJF. "legacy" = the former
        # definitions (no common cell). See config_defaults.yml.
        africa = list(window_end = "legacy", continuation_max_lead = 0L,
                      max_counterfactual_winter = "")
      ),
      # 09 only.
      derivation = list(
        smoothing_start   = "2011-01-01",
        n_irr_draws       = 400L,
        # ENSO exposure-response lag selection: for duration/monthly,
        # require the select_best_lag() pick to be >= the occurrence track's pick
        # for the same group. FALSE uses each track's own unconstrained pick.
        enforce_occ_floor = FALSE,
        # Forecast arms launching at different months: "stitch" (window from the
        # earlier launch, factual series from the smoothing up to its own launch)
        # or "common_start" (window from the later launch). See config_defaults.yml.
        arm_alignment = "stitch",
        # Paired-bootstrap resamples of report_attrib_exposure /
        # report_attrib_metrics / report_event_attrib (Monte-Carlo precision of
        # the replicate means). 0 = no bootstrap: interval columns NA, point
        # estimates unchanged.
        attrib_boot_b = 500L,
        # Explosive-replicate modes every attribution cell is derived under:
        # "drop" (a pair is dropped when either arm's group window cases reach
        # thresholds.explosive_case_cap) and/or "none" (every pair kept).
        explosive_modes = c("drop", "none"),
        # Explosive guard scope of the af_posterior steps:
        # "unit" drops a replicate explosive in any group from its whole
        # simulation unit; "group" from its group only.
        af_guard_scope = "unit",
        # Phase-11 "key numbers" report. pop_ref_year is the year the canonical
        # gid_population table is joined at for the population-covered and
        # significant-effect-population tables; duration_ratio_horizon is the
        # K-month truncation for the duration restricted-mean-duration (RMD)
        # ratio.
        pop_ref_year           = 2025L,
        # Population of the per-capita rates (inc_ann_100k, an_per_100k):
        # "ref" = pop_ref_year; "window" = the
        # month-weighted mean of the canonical gid population over the months
        # each table reads. "People living in" counts stay on
        # pop_ref_year either way.
        rate_pop_year          = "ref",
        duration_ratio_horizon = 24L,
        # ONI values (vs neutral 0) the composite ADM1 IRRs, the significant-
        # and ENSO-sensitive-population tables and the duration RMD ratios are
        # evaluated at. Labels come from figures_cb.R::enso_conds_named(); the
        # default equals the code constant COMPOSITE_CONDS.
        enso_conds             = c(-2, -1, 1, 2),
        # Figure 1d's rainfall tables (rain_enso step). See config_defaults.yml.
        rain_enso = list(
          baseline             = "1990_2019",
          lags                 = c(-11L, 11L),
          dry_mm_month         = 10,
          min_n                = 30L,
          season_share         = 0.5,
          season_min_annual_mm = 100,
          season_pop_year      = 2020L
        ),
        # Supplement fig. SE (seasonality step): the analysis rows of the three
        # tracks' observed monthly rates (applyAdminLevelFilter() rule of the
        # fits), the events (cases: episode-months) below which a month is
        # sparse and takes no part in the peak, and the interval level. See
        # config_defaults.yml.
        seasonality = list(
          admin_level = "ADM1-ADM2",
          min_n       = 5L,
          interval    = 0.95
        ),
        # Supplement fig. WR (report_who_enso_ratio step): the
        # observed ratio of 2026 WHO cases to the mean of the ENSO-neutral years
        # against the modelled ratio of the smoothing to the ENSO counterfactual
        # over their common months, for the countries of `groups`. Years from
        # `first_year`; a year's ENSO state under each of `year_definitions`
        # (enso.R::ENSO_YEAR_DEFINITIONS); a target-year WHO series stopping more
        # than `censor_weeks` before the cut-off is censored; `pseudo_count`
        # cases on both sides of every ratio of means; `interval` level of the
        # bootstraps (`n_boot`, `seed`) and of the per-replicate spread. See
        # config_defaults.yml.
        who_enso_ratio = list(
          groups           = "enso_west",
          year_definitions = c("annual_mean", "lag_occ", "djf_start", "djf_end"),
          first_year       = 2010L,
          censor_weeks     = 8L,
          pseudo_count     = 1,
          interval         = 0.95,
          n_boot           = 2000L,
          seed             = 20261107L
        ),
        # Supplementary table of ENSO contrasts, all episodes of the analysis set vs those with
        # at least `min_confirmed` confirmed cases (report_confirmed_agreement):
        # a confirmed contrast needs `min_neutral` neutral and
        # `min_phase` phase onsets; both sets restricted to the confirmed
        # onsets' window ("confirmed") or the full record ("full"); duration and
        # monthly-case intervals from `boot_b` location-period bootstraps
        # (`seed`) at `level`. See config_defaults.yml.
        confirmed_agreement = list(
          min_confirmed = 1L,
          min_neutral   = 2L,
          min_phase     = 5L,
          window        = "confirmed",
          boot_b        = 1000L,
          seed          = 20261110L,
          level         = 0.95
        ),
        # ENSO phase rule of report_smoothing's by-phase table:
        # "threshold" (monthly ONI at +/-0.5, no persistence; legacy) or
        # "episode" (NOAA CPC episodes, the runs Figures 1 and 3 shade).
        enso_phase = list(smoothing = "threshold"),
        # Lag averaging (09 step lag_averaging): from
        # the 05 lag sidecars (05_fit_postprocess_loo_lags.R --lag_sidecar),
        # Pareto-k per lag, the k-robust elpd difference, stacking / pseudo-BMA+
        # weights over the lag sets ("all" = every fitted lag, "eligible" = the
        # deploy's lag floor), the stacked population ERF and the IRR at
        # enso_conds. enable false = the step is off. report_min_weight is the
        # display cut of the weight tables only (the mixture is never truncated).
        lag_averaging = list(
          enable            = FALSE,
          methods           = c("stacking", "pseudobma_plus"),
          lag_sets          = c("all", "eligible"),
          k_threshold       = 0.7,
          report_min_weight = 0.01
        ),
        # Hindcast scoring method (utils/hindcast_skill.R::hindcast_skill_bundle()).
        # lead_breaks: lead classes in months after launch (7, 12 -> 0-7, 8-12,
        # 13+). occ_*: occurrence target grain ("window" headline,
        # "calendar_year" legacy), climatology references and CRPS estimators.
        # case_*: case-total supports, references and rules. zero_floor floors
        # log-score probabilities at 1/(2 n_rep).
        hindcast = list(
          lead_breaks      = c(7L, 12L),
          occ_grains       = c("window", "calendar_year"),
          occ_references   = c("insample_gid", "loo_gid", "loo_gid_fullhist",
                               "regional", "zero", "persistence"),
          occ_rules        = c("crps_emp", "crps_fair"),
          case_supports    = c("obs_positive_keys", "exact_gid", "atrisk_gids",
                               "atrisk_keys", "footprint"),
          case_references  = c("loo_climatology", "zero", "persistence"),
          case_rules       = c("crps_raw", "crps_log1p", "pinball"),
          pinball_taus     = c(0.05, 0.25, 0.5, 0.75, 0.95),
          interval_widths  = c(0.5, 0.8, 0.9, 0.95, 0.99),
          reliability_bins = 10L,
          zero_floor       = TRUE,
          rank_seed        = 42L,
          # Monthly and season case skill (hindcast_skill.R section 5b); off
          # while case_skill_supports is empty. Units per (area, month) in each
          # country's window (case_month_grains) and per (country, lead class)
          # season total; references on the run's fixed support
          # (case_skill_references), fair CRPS / threshold-weighted CRPS rules
          # (thresholds in cases per month / per season total), climatology
          # window +/- case_clim_halfwidth months, population scaling ("wpp" |
          # "none"), surveillance mask ("coverage" | "none"), minimum reference
          # values of a season total.
          case_skill_supports      = character(0),
          case_month_grains        = c("country", "adm1"),
          case_skill_references    = c("clim_pooled", "clim_launch_state",
                                       "persistence_season", "persistence_last"),
          case_skill_rules         = c("crps_raw", "crps_log1p", "twcrps"),
          case_clim_halfwidth      = 1L,
          case_pop_scale           = "none",
          case_mask                = "none",
          case_tw_threshold_month  = 50,
          case_tw_threshold_season = 50,
          case_season_min_ref      = 5L,
          # Season-year bootstrap CIs on every CRPSS (0 = none) and per-area
          # median / positive-share summaries.
          bootstrap_n          = 0L,
          bootstrap_seed       = 20261101L,
          bootstrap_level      = 0.95,
          skill_area_summaries = FALSE
        )
      ),
      # BOTH: 09 computes against these, 10 states them in captions.
      thresholds = list(
        prob_threshold = 0.50,
        low_signal     = 2.0,
        # Phase-11 forecast threshold-population grid: for each (cutoff p ×
        # magnitude m) pair, population in admin units with >= p posterior
        # probability of a >= m fractional increase (or decrease) in the horizon.
        pop_prob_cutoffs = c(0.5, 0.8, 0.95),
        pop_magnitudes   = c(0.25, 0.5),
        # Phase-11 ENSO-sensitive population: an ADM1 is sensitive at a
        # condition when ANY track's IRR is beyond 1 +/- this fraction, in
        # either direction (|log IRR| > log(1 + sensitive_irr)).
        sensitive_irr    = 0.10,
        # Phase-11 ENSO exposure tables: magnitude axis (IRR beyond 1 + m in
        # the adverse or protective direction; 0 = any effect). The certainty
        # axis reuses pop_prob_cutoffs.
        enso_magnitudes  = c(0, 0.10, 0.25, 0.50),
        # report_attrib_exposure: Childs magnitude axis (mean forecast >= (1 + m)
        # x mean counterfactual) and the IPCC robustness shares.
        attrib_magnitudes = c(0.25, 0.5, 1),
        robust_exceed     = 0.66,
        robust_agree      = 0.80,
        # report_attrib_metrics: the cases_ge event threshold (the database's
        # outbreak-size convention: OutbreakExtractR episodes have >= 50
        # suspected cases; no code enforces it), the one-sided exact level of
        # the PR lower bound when the counterfactual never has the event, and
        # the S/N emergence levels.
        pr_case_threshold = 50,
        pr_alpha          = 0.05,
        sn_levels         = c(1, 2),
        # The explosive-replicate cap of the "drop" mode (group window cases).
        explosive_case_cap = 2e5,
        # The per-gid monthly cap of the "drop" mode: a replicate with any gid
        # of the group at or above it in any window month is explosive too.
        # Inf = none.
        explosive_gid_month_cap = Inf,
        # af_posterior steps: a scope or ADM1 whose Monte-Carlo
        # noise share is not below it has an unidentified posterior spread
        # (adjusted interval NA, reported apart).
        af_noise_share_max = 0.8
      ),
      # 10 only.
      style = list(
        out_dir  = "figures/manuscript",
        formats  = c("pdf", "png"),
        font     = "Arial",
        dpi      = 600L,
        map_lo   = 0.3,
        map_hi   = 3,
        erf_xlim = c(-2.5, 2.5),
        # 11 only: one CSV per printed table (Science's machine-readable tables).
        tables_dir = "figures/manuscript/tables",
        # 11 only: the report_enso_exposure grid cell the manuscript headlines
        # (population where P(effect beyond 1 + magnitude in direction) >= cutoff).
        enso_headline = list(direction = "adverse", track = "any_credible",
                             magnitude = 0.25, cutoff = 0.95),
        # 10 + 11: the report_attrib_exposure magnitude the manuscript headlines.
        attrib_headline = list(magnitude = 0.5, window = "attribution", explosive = "drop"),
        # 11 only: the report_attrib_metrics cells the manuscript headlines.
        attrib_metrics_headline = list(sn_level = 1.0, pr_event = "onset_any"),
        # 10 + 11: group display labels (results_core.R::group_display_meta()):
        # "registry" = the grouping sidecar's names; "lag" = class letters
        # ordered by the deployed occurrence lag, longest lead first ("Class A
        # (ENSO lead 6)", short "A -6"), from the group_labels artefact;
        # "centroid" = the sidecar's display-only centroid names ("West [GW]",
        # short "GW", order centroid_order). Ids and colours never change. Not
        # stamped (presentation).
        group_labels = "registry",
        # 10 + 11: p_max (largest posterior class probability of a country in
        # the deployed sweep cell, report_membership) below which a country is
        # marked as uncertainly assigned. Display only.
        membership_low = 0.6,
        # 10 only, supplement fig. WR: the year definition and the
        # counterfactual analogue set drawn (every one is derived). Not stamped.
        who_enso_ratio = list(year_definition = "annual_mean", analog_set = "neutral"),
        # 11 only: data release of each D / T table:
        # "supp" = Supplementary Data with the manuscript (policy-relevant
        # outputs), "github" = the repository's data folder. 11 copies each
        # CSV into <tables_dir>/supp_data/ or <tables_dir>/github/ and writes
        # data_release_manifest.csv (results_report.R::report_release_split()).
        # Not stamped (presentation).
        data_release = list(
          D1 = "github", D2 = "supp", D3 = "github", D4 = "github", D5 = "github", D6 = "github",
          D7 = "github", D8 = "github", D9 = "github", D10 = "supp", D11 = "supp", D12 = "supp",
          D13 = "supp", D14 = "supp", D15 = "supp", D16 = "supp", D17 = "supp", D18 = "supp",
          D19 = "supp", D20 = "supp", D21 = "supp", D22 = "github", D23 = "github", D24 = "github",
          D25 = "github", D26 = "github", D27 = "github", D28 = "github", D29 = "github", D30 = "supp",
          D31 = "github", D32 = "github", D33 = "github", D34 = "github", D35 = "github", D36 = "github", D37 = "github", D38 = "github", D39 = "github", D40 = "github", T1 = "github",
          T2 = "github", T3 = "supp", T4 = "supp", T5 = "supp"),
        # 11 only: table slugs withheld from the release (listed in the
        # manifest as "withheld", not copied): the ADM1 observed counts of
        # the past-event attribution. Not stamped (presentation).
        data_release_exclude = c("D16_event_attrib_observed_adm1")
      )
    ),

    # -------------------------------------------------------------------------
    # deploy: master contract for deploying a prod_run<N> occurrence batch into
    # the shared enso_best_lag_fits archive and threading its run-tag + the
    # nbr-contagion enable flag through historical smoothing (08), forecast (07)
    # and hindcast (07). Consumed by analysis/utils/deploy_prod_run.R; schema +
    # logic live in utils/deploy_core.R. Mirrors the self-contained
    # manuscript-results contract (the `results` block above): the deploy YAML
    # restates every leaf here, so the file alone is the whole deployment-
    # decision record. Defaults are deliberately INERT (nbr_contagion.enable =
    # FALSE, hindcast.launch_years = "") and name NO run: every run-specific
    # leaf is "" and deploy_core.R::deploy_require_keys() refuses a deploy YAML
    # that leaves one unset. Each rollout lives entirely in
    # analysis/configs/deploy/prod_run<N>.yml.
    # -------------------------------------------------------------------------
    deploy = list(
      version = 1L,

      # RUN IDENTITY. Every leaf that names a prod_run<N> (fits.*_dir,
      # templates.*, run_tags.*, outputs.*_config_outdir, smoothing_input)
      # defaults to "" -- the unset sentinel (a NULL leaf cannot round-trip
      # .merge_recursive; see results$deploy_config above). The run is named
      # ONLY in analysis/configs/deploy/prod_run<N>.yml, and
      # deploy_core.R::deploy_require_keys() refuses a YAML that leaves one
      # unset, so a stale run can never be inherited from these defaults.
      #
      # Source batch dirs (analysis/configs/<dir>/) the best-lag fits are
      # collected FROM, plus the shared archive they are copied TO and the
      # groups to collect. Tracks may come from different runs (a batch dir is
      # reused until its track is refit).
      fits = list(
        occ_dir = "",
        dur_dir = "",
        mon_dir = "",
        dest    = "generated_data/enso_best_lag_fits",
        groups  = c("east", "horn", "south", "west")
      ),

      # ENSO best-lag selection rule (utils/loo_lags.R::lag_selection_cfg()):
      # longest_within_se (the rule every run up to 14 deployed: the longest
      # lag >= min_lag within se_mult SEs of the rank-1 fit, leads never
      # selected) | best_elpd (best-ranked lag >= min_lag, leads included when
      # min_lag is null; downstream tracks floored at occ_lag - occ_floor_offset).
      # min_lag "rule_default" = the rule's default (0 / -Inf); occ_floor_offset
      # "none" = no occurrence floor (strings, not null: the config merge drops
      # null leaves).
      lag_selection = list(
        rule             = "longest_within_se",
        se_mult          = 2,
        min_lag          = "rule_default",
        occ_floor_offset = "none"
      ),

      # Convergence gate on deploy_fits / gen_grouping
      # (deploy_core.R::deploy_convergence_cfg()). action: off (nothing is read,
      # printed or written) | warn | refuse (stop
      # before the prune, dest untouched). scope: population (every parameter
      # but ^(r_|z_|zs_)) | all. rhat_max / ess_min "default" =
      # diagnostics.R::DIAG_THRESHOLDS (1.01 / 400; never restated here).
      # allow: <track>:<short group> ids that pass although they fail (stated
      # in Methods / S8).
      convergence = list(
        action        = "off",
        scope         = "population",
        rhat_max      = "default",
        ess_min       = "default",
        max_divergent = 0L,
        allow         = character(0)
      ),

      # Distance-weighted neighbour-contagion covariate parameters. These drive
      # the occurrence artefact filename token via paths.R (scale="raw" ->
      # "_nbr500", scale="rownorm" -> "_nbr500n"), so deploy_collect_fits()
      # reconstructs the SAME canonical names the fits were written with.
      contagion = list(
        rho_km    = 500L,
        scale     = "raw",
        partition = "none",
        # partition = "registry": the registry restricting the neighbour set
        # (covariate nbr_lag1_r<rho>[_rn]_part<token>, utils/paths.R). Also
        # written into the generated sim configs as
        # simulator.nbr_contagion.partition_source. "" otherwise; under
        # partition_by = registry it must refine partition_source (default:
        # the same registry).
        partition_source = ""
      ),

      # Gate for the spatPomp/v1 neighbour-contagion kernel. Threaded into the
      # generated smoothing/forecast/hindcast templates as
      # simulator.nbr_contagion.enable; resolve_use_nbr() (spatial_contagion.R)
      # activates the kernel iff this is TRUE AND a deployed fit carries the nbr
      # FE. FALSE by default => no neighbour contagion.
      nbr_contagion = list(
        enable = FALSE,
        rho_km = 500L
      ),

      # Cross-track episode definition for the WHOLE rollout (see the top-level
      # `episodes:` block). The single authority: deploy_collect_fits() refuses
      # per-lag configs / fits / prep artefacts of another key, the smoothing
      # stage refuses a historical template declaring another key, and the
      # forecast/hindcast generators write it into every generated config
      # (utils/provenance.R::applyEpisodeKeyToSimConfig()). "lp_season" keeps
      # every existing rollout byte-identical. national_zero_fill: the rollout's
      # occurrence fits are on the national zero-fill panel ("_nz"; gid_season only),
      # written into every generated sim config and checked against the
      # per-lag occurrence configs.
      # national_zero_rule / at_risk_start / adm1_overlap:
      # the occurrence-panel spec fields (paths.R::occPanelSpecOf());
      # ""  / "season" / "drop" = today's panel, written into nothing.
      episodes = list(
        key           = "lp_season",
        national_zero_fill = FALSE,
        national_zero_threshold = 0L,
        national_zero_rule = "",
        at_risk_start = "season",
        adm1_overlap  = "drop"
      ),

      # The three spatPomp templates the HPC stages consume. GENERATED by
      # `deploy_prod_run.R --stage gen_templates` (skeleton in
      # analysis/configs/sim_templates/ + the `sim:` overlay below) at
      # analysis/configs/<run_tags.<arm>>/<run_tags.<arm>>.yml; these paths
      # must equal that location (deploy_core.R::deploy_sim_template_path()).
      templates = list(
        historical     = "",
        forecast       = "",
        counterfactual = ""
      ),

      # Base run_tags for the four spatPomp phases (historical smoothing,
      # forecast, counterfactual, hindcast). The generators suffix these with
      # _<group>[_<year>] per run.
      run_tags = list(
        historical     = "",
        forecast       = "",
        counterfactual = "",
        hindcast       = ""
      ),

      # Where the generators WRITE and the manifests they overwrite. The
      # manifests are SHARED with manuscript.yml's results.inputs.manifest, so
      # the 09/10 phase reads this deployment immediately after gen_forecast /
      # gen_hindcast.
      outputs = list(
        forecast_config_outdir = "",
        forecast_manifest      = "analysis/configs/forecast_horizons.yml",
        # Attribution manifest written by --stage gen_attrib
        # (utils/make_forecast_attrib_configs.R): the forecast manifest with the
        # factual arm launched at the counterfactual handoff, read by 09's
        # attribution steps through results.inputs.attrib_manifest. "" = no
        # attribution factual (the stitched forecast factual is used).
        forecast_attrib_manifest = "",
        hindcast_config_outdir = "",
        hindcast_manifest      = "analysis/configs/hindcast_horizons.yml"
      ),

      # The historical-smoothing cache the forecast/counterfactual/hindcast arms
      # hand off from. MUST be the cache produced by the historical template
      # above (== generated_data/sim_spatpomp/<run_tags.historical>/cache).
      smoothing_input = "",

      grouping           = "enso_grouping",
      # brms_kl latent-class sweep that DERIVES the fit registry `grouping`
      # (deploy_prod_run.R --stage gen_grouping, utils/kl_grouping.R; relay
      # phase `classes` in bash/tools/run_pipeline.sh). batch "" = no sweep:
      # `grouping` is an existing registry, used as is. cell = the pinned fit id
      # (kl_sweep.R::kl_sweep_inventory() `id`, e.g. "ns5_K4"); never chosen
      # automatically. chains: "agree" (the agreeing chains recorded by
      # 05_fit_postprocess_kl_sweep.R) | "all" | ids "2,6". fill_from: registry
      # whose countries absent from the fit are carried over ("" = none).
      # meta_from: registry whose sidecar labels / colours / order are reused
      # for the derived groups it shares. n_agree / mode_tol: how the
      # comparison picks each fit's agreeing chains (kl_sweep.R::
      # kl_agreeing_chains(); 05 --n_agree / --mode_tol): n_agree 0 = the
      # dominant mode (chains whose mean class centres lie within mode_tol
      # lags), N > 0 = the N chains whose centres agree best. Passed to 05 by
      # stages/04_fit_kl.sh --deploy; gen_grouping refuses a comparison run
      # with other values when chains = "agree".
      kl_sweep = list(
        batch     = "",
        cell      = "",
        chains    = "agree",
        n_agree   = 0L,
        mode_tol  = 1.5,
        fill_from = "",
        meta_from = "enso_grouping"
      ),
      analog_sets        = "all",
      primary_analog_set = "all",

      # Forward-partition for the forecast/counterfactual/hindcast arms.
      # "country" (default) = today's per-country SLURM-array dispatch, byte-for-byte
      # identical generated configs. "enso_region" = one forward partition per ENSO
      # group (ADM1 contagion kernel couples member countries across borders); the
      # generators write filter.partition_by and the deploy tool prints
      # `--partition_by enso_region` on the submit command. Smoothing stays per-country.
      partition_by = "country",
      # partition_by = "registry": the country-level registry
      # (analysis/groupings/<name>.csv) whose groups are the simulation units
      # (TUs), independent of `grouping` (the fit registry). One forward run per
      # unit, simulated to the max member firm horizon; v2 manifests carry
      # units + per-country firm horizons. "" unless partition_by = "registry".
      partition_source = "",

      # ENSO past the CFSv2 span (utils/enso.R): "neutral" (0; the historical
      # behaviour), "ar1" (each replicate's member decays from its last value
      # with the AR(1) coefficient of the observed series) or "analog" (k-NN
      # analog continuation from the piControl `library`; k / blend_months /
      # min_sep as simulator.enso.continuation). Written into the forward
      # configs as simulator.enso.continuation when not "neutral".
      # min_sep: minimum months between two analog starts within one model
      # (the continuation 24; the counterfactual pool's
      # deploy.enso_counterfactual.min_sep_months 12, which allows adjacent
      # same-season years; utils/enso.R::ENSO_COUNTERFACTUAL_DEFAULTS).
      # seam: last_cfsv2 (the analog continuation joins at
      # the last CFSv2 month; written into nothing) | last_full_oni (at the last
      # full 3-month CFSv2 mean; needs method analog).
      enso_continuation = list(method = "neutral",
                               library = "generated_data/enso_continuation_library.qs",
                               k = 30L, blend_months = 3, min_sep = 24L,
                               max_per_model = Inf, exclude_models = c("CM2.1-NMME"),
                               seam = "last_cfsv2"),

      # Counterfactual ENSO analog pool
      # (enso.R::ensoCounterfactualCfg(), ENSO_COUNTERFACTUAL_DEFAULTS). match:
      # legacy (the shared hand-built Lenssen cache; nothing is built or
      # recorded) | library (per-run cache drawn from the continuation library
      # by gen_forecast / gen_events: n_analogs per model, min_sep_months
      # between analog starts within one model, faded join over blend_months,
      # sets labelled on the `winters` read (primary | all_read) over read_span
      # (direct | unit), other_winters neutral | any; a primary set under
      # min_set_size is refused; seam_at: launch | launch_minus_1, library only).
      enso_counterfactual = list(
        match          = "legacy",
        n_analogs      = 50L,
        min_sep_months = 12L,
        blend_months   = 3,
        winters        = "primary",
        read_span      = "direct",
        other_winters  = "neutral",
        min_set_size   = 30L,
        seam_at        = "launch"   # library: join the analog paths to the observed ONI of the launch month (launch_minus_1 = the month before)
      ),

      # Hindcast-only knobs (make_hindcast_configs.R). launch_years EXPLICIT in
      # the deploy YAML (default "" here = auto-derive from the smoothing cache,
      # which only works once the cache exists) so all config-wiring runs
      # locally before the HPC smoothing job finishes.
      # season_end: "by_group" = the
      # season_end_by_group month of each country's fit group (unchanged
      # output); "by_country_from_panel" = each country's own cholera-year end
      # (min_seas - 1, modal over its gids in climate_panel; ties -> the later
      # end), the window ending at the first such month at least
      # season_min_window_months after the launch (utils/season_calendar.R).
      # season_window: enable = one run per (unit, launch month,
      # season year), each country's season from the month after its trough
      # (start floored at launch_floor_month) to the next trough + max(0, min
      # lag); replaces launch_month / season_end for the hindcast only (the
      # event attribution still reads season_end). trough_rule onset_rate =
      # the country's onset-rate trough in its occurrence panel where it has
      # >= min_onsets onsets and stability >= min_stability (bootstrap over
      # years), else its modal min_seas; trough_years "" = every panel year;
      # survival_panel "" = the template episodes' panel.
      hindcast = list(
        launch_month        = "05-01",
        launch_years        = "",
        season_end          = "by_group",
        season_min_window_months = 6L,
        season_end_by_group = "",
        horizon_months      = 12L,
        add_best_lag        = TRUE,
        # posterior-draw mode and replicates of the hindcasts: "" / 0 =
        # patch_hindcast_config()'s random at the forecast template's count;
        # "joint" = the joint draw design.
        param_draw_mode     = "",
        n_replicates        = 0L,
        climate_panel       = "generated_data/sim_climate_panel_base-1990_2019.rds",
        season_window       = list(
          enable             = FALSE,
          launch_floor_month = 5L,
          trough_rule        = "onset_rate",
          min_onsets         = 30L,
          min_stability      = 0.8,
          trough_years       = "",
          smooth_months      = 3L,
          bootstrap          = 500L,
          seed               = 1L,
          survival_panel     = ""
        )
      ),

      # Past-ENSO-event attribution (--stage gen_events,
      # utils/make_event_attribution_configs.R; 09 step report_event_attrib).
      # events: comma-separated ids of `registry` (analysis/configs/enso_events_v2.yml);
      # "" = none, and gen_events refuses to run. Per event and unit, a factual
      # (observed ENSO) and a counterfactual (Lenssen `analog_set` from a per-event
      # cache under lenssen_dir/<event>/) config handing off from the smoothing at
      # the event's launch. climate_panel must cover every window (the smoothing's
      # extended panel); season_end "" = hindcast.season_end, season_end_by_group
      # "" = hindcast.season_end_by_group.
      # run_tag_base / config_outdir / manifest name a run: required when events
      # is set (no default).
      event_attribution = list(
        events              = "",
        registry            = "analysis/configs/enso_events_v2.yml",
        analog_set          = "neutral",
        lenssen_dir         = "generated_data/lenssen_events",
        climate_panel       = "generated_data/sim_climate_panel_base-1990_2019.rds",
        season_end          = "",
        season_end_by_group = "",
        run_tag_base        = "",
        config_outdir       = "",
        manifest            = "",
        # Posterior-draw mode of both arms: "" = the hindcast default (random);
        # "cycle" = a fixed draws x replicates design (max_draws from the
        # templates, n_replicates below; max_draws must be coprime to the
        # smoothing's replicate count).
        param_draw_mode     = "",
        # Replicates of every event run (both arms): 0 = the forecast
        # template's simulator.n_replicates.
        n_replicates        = 0L
      ),

      # Per-run overlay for the generated spatPomp templates (--stage
      # gen_templates). Each arm is built as
      #   analysis/configs/sim_templates/<arm>.yml + sim.shared + sim.<arm>
      # (later layers win) + the leaves derived from this deploy block. Keys are
      # DOTTED leaf paths into the template (e.g. simulator.filter.measurement),
      # so every per-run decision is one auditable line. Contents are free-form:
      # test_deploy_config_coverage.R treats sim.<arm> as opaque. Empty by
      # default; gen_templates refuses a run whose overlay omits a required leaf
      # (deploy_core.R::.DEPLOY_SIM_REQUIRED_LEAVES).
      sim = list(
        shared         = list(),
        historical     = list(),
        forecast       = list(),
        counterfactual = list()
      )
    )
  )
}


#' load_grouping
#'
#' Read a country-grouping registry CSV (columns `iso3`, `group`) from
#' `analysis/groupings/<source>.csv`. Defined in utils/config.R (not utils.R) because
#' utils/config.R is the lowest-level shared module, sourced first by all entry points
#' (fit scripts, batch maker, simulator, tests), so both config-resolution and
#' the batch maker can use it without depending on utils.R load order.
#'
#' @param source character — registry name (no extension)
#' @param groupings_dir character — directory holding registries
#' @return data.frame(iso3, group), both character
load_grouping <- function(source, groupings_dir = "analysis/groupings") {
  reg  <- .read_registry_csv(source, groupings_dir, c("iso3", "group"), "load_grouping")
  path <- attr(reg, "path")
  dup <- reg$iso3[duplicated(reg$iso3)]
  if (length(dup))
    stop("load_grouping: duplicate iso3 in registry ", path, ": ",
         paste(unique(dup), collapse = ", "))
  rownames(reg) <- NULL
  attr(reg, "path") <- NULL
  reg
}

#' .read_registry_csv
#'
#' Shared reader of the `analysis/groupings/<source>.csv` registries
#' (load_grouping(), load_adjacency()): existence check, all-character read,
#' required-column check, whitespace trim, rows with an empty required cell
#' dropped. Returns the required columns only, with the file path as
#' attribute "path" for the caller's messages.
#' @noRd
.read_registry_csv <- function(source, groupings_dir, cols, who) {
  if (is.null(source) || !nzchar(source))
    stop(who, ": `source` is required")
  path <- file.path(groupings_dir, paste0(source, ".csv"))
  if (!file.exists(path))
    stop(who, ": registry not found: ", path)
  reg <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character",
                         comment.char = "#")
  if (!all(cols %in% names(reg)))
    stop(who, ": registry must have columns ", paste(cols, collapse = ","), ": ", path)
  for (cl in cols) reg[[cl]] <- trimws(reg[[cl]])
  keep <- Reduce(`&`, lapply(cols, function(cl) nzchar(reg[[cl]])))
  reg  <- reg[keep, cols, drop = FALSE]
  attr(reg, "path") <- path
  reg
}

#' load_adjacency
#'
#' Read a country-adjacency registry (`analysis/groupings/<source>.csv`,
#' columns `iso3_a,iso3_b`, one undirected edge per row; written by
#' utils/make_country_adjacency.R). Used by the brms_kl latent-class backend's
#' spatial (ICAR) membership term. Self-loops are dropped and each edge is
#' returned once with iso3_a < iso3_b.
#'
#' @param source character — registry name (no extension)
#' @param groupings_dir character — directory holding registries
#' @return data.frame(iso3_a, iso3_b), both character, sorted
load_adjacency <- function(source, groupings_dir = "analysis/groupings") {
  adj <- .read_registry_csv(source, groupings_dir, c("iso3_a", "iso3_b"), "load_adjacency")
  a <- pmin(adj$iso3_a, adj$iso3_b)
  b <- pmax(adj$iso3_a, adj$iso3_b)
  out <- unique(data.frame(iso3_a = a, iso3_b = b, stringsAsFactors = FALSE)[a != b, , drop = FALSE])
  out <- out[order(out$iso3_a, out$iso3_b), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------------------
# Registry metadata, identity and label helpers
# ---------------------------------------------------------------------------
# A grouping registry `analysis/groupings/<source>.csv` (iso3,group) may carry
# an optional sidecar `analysis/groupings/<source>.yml`:
#
#   label_prefix: "enso_"      # resolves short forms ("east" -> "enso_east")
#   order: [enso_east, ...]    # display / iteration order
#   groups:
#     enso_east: {label: "East Africa", short: "East", short_fig: "East",
#                 colour: "#1b9e77"}
#
# Every field falls back to the registry's own group labels, so a registry
# without a sidecar works everywhere (label = short = short_fig = the group
# string, order = sorted labels, colours from a fixed qualitative palette).

#' registry_path
#' @return path of `<groupings_dir>/<source>.<ext>`
registry_path <- function(source, groupings_dir = "analysis/groupings", ext = "csv") {
  if (is.null(source) || !nzchar(source))
    stop("registry_path: `source` is required", call. = FALSE)
  file.path(groupings_dir, paste0(source, ".", ext))
}

#' registry_digest
#'
#' Content hash (xxhash64 of the CSV bytes) of a grouping registry, the
#' identity stamped wherever a registry shapes an artefact (03 contagion
#' provenance, smoothing cache `meta.yml`, forecast/hindcast manifests, sim
#' run_meta). Delegates to base.R::content_digest() on the file's raw bytes
#' (call-time dependency; base.R is the universal floor every caller sources),
#' which is xxhash64 of the bytes as they are: the value equals the former
#' `digest::digest(file = path, algo = "xxhash64")` for every registry
#' (pinned in tests/test_t2_shared_infra.R). NA when the optional `digest`
#' package is absent (a missing package weakens the check, never breaks the
#' pipeline).
registry_digest <- function(source, groupings_dir = "analysis/groupings") {
  path <- registry_path(source, groupings_dir)
  if (!file.exists(path))
    stop("registry_digest: registry not found: ", path, call. = FALSE)
  content_digest(readBin(path, what = "raw", n = file.size(path)))
}

#' grouping_label_token
#'
#' Filename-safe token of a group label: non-alphanumerics stripped, the rule
#' make_batch_configs.R::make_config_filename() applies to every extra grid
#' dimension ("enso_east" -> "ensoeast"; batch configs `_grp-ensoeast.yml`).
grouping_label_token <- function(x) gsub("[^A-Za-z0-9]", "", as.character(x))

.GROUP_FALLBACK_COLOURS <- c("#1b9e77", "#d95f02", "#66a61e", "#7570b3",
                             "#e7298a", "#e6ab02", "#a6761d", "#666666")

#' load_grouping_meta
#'
#' Registry + sidecar metadata. Returns list(source, labels (character, in
#' `order`), label_prefix, order, label, short, short_fig, colour,
#' centroid_short, centroid_label, centroid_order) where label .. centroid_label
#' are named character vectors keyed by group label (the centroid fields NA
#' when the sidecar has none) and centroid_order is the sidecar's display order
#' of group_labels "centroid" (character(0) when absent).
load_grouping_meta <- function(source, groupings_dir = "analysis/groupings") {
  reg    <- load_grouping(source, groupings_dir)
  labels <- sort(unique(reg$group))
  side_p <- registry_path(source, groupings_dir, "yml")
  side   <- if (file.exists(side_p)) yaml::read_yaml(side_p) else list()
  unknown <- setdiff(c(unlist(side$order), names(side$groups)), labels)
  if (length(unknown))
    stop(sprintf("load_grouping_meta: sidecar %s names group(s) absent from the registry: %s",
                 side_p, paste(unknown, collapse = ", ")), call. = FALSE)
  ord <- as.character(unlist(side$order))
  ord <- c(ord, setdiff(labels, ord))
  field <- function(f, fallback) {
    v <- vapply(ord, function(g) {
      x <- side$groups[[g]][[f]]
      if (is.null(x) || !nzchar(as.character(x))) NA_character_ else as.character(x)
    }, character(1))
    miss <- is.na(v)
    v[miss] <- fallback[miss]
    stats::setNames(v, ord)
  }
  label <- field("label", ord)
  short <- field("short", label)
  colour_fb <- rep_len(.GROUP_FALLBACK_COLOURS, length(ord))
  # Display-only centroid names (results.style.group_labels "centroid";
  # results_core.R::group_display_meta()): NA / empty when the sidecar has none.
  na_fb <- rep(NA_character_, length(ord))
  c_ord <- as.character(unlist(side[["centroid_order"]]))
  if (length(c_ord) && !setequal(c_ord, ord))
    stop(sprintf("load_grouping_meta: sidecar %s centroid_order must list every group once (groups: %s).",
                 side_p, paste(ord, collapse = ", ")), call. = FALSE)
  list(source         = source,
       labels         = ord,
       label_prefix   = as.character(side$label_prefix %||% ""),
       order          = ord,
       label          = label,
       short          = short,
       short_fig      = field("short_fig", short),
       colour         = field("colour", colour_fb),
       centroid_short = field("centroid_short", na_fb),
       centroid_label = field("centroid_label", na_fb),
       centroid_order = c_ord)
}

#' resolve_group_label
#'
#' Map group names that may be short forms ("east") or full labels
#' ("enso_east") onto registry labels. A full label wins; otherwise
#' `label_prefix` is prepended. Errors on a name that resolves to no label.
resolve_group_label <- function(meta, x) {
  x <- as.character(x)
  out <- ifelse(x %in% meta$labels, x, paste0(meta$label_prefix, x))
  bad <- x[!out %in% meta$labels]
  if (length(bad))
    stop(sprintf("resolve_group_label: group(s) %s not in registry '%s' (labels: %s)",
                 paste(bad, collapse = ", "), meta$source,
                 paste(meta$labels, collapse = ", ")), call. = FALSE)
  out
}

#' grouping_source_from_token
#'
#' The unique registry whose source name has alphanumeric stem `token`
#' (grouping_label_token(source) == token). Inverse of the `_part<token>`
#' neighbour-contagion covariate suffix (utils/paths.R::nbrCovariateName()).
#' Errors on zero or several matches, so a token always names one registry.
grouping_source_from_token <- function(token, groupings_dir = "analysis/groupings") {
  srcs <- sub("\\.csv$", "", list.files(groupings_dir, pattern = "\\.csv$"))
  hit  <- srcs[tolower(grouping_label_token(srcs)) == tolower(token)]
  if (length(hit) != 1L)
    stop(sprintf("grouping_source_from_token: token '%s' matches %d registr%s in %s%s",
                 token, length(hit), if (length(hit) == 1L) "y" else "ies",
                 groupings_dir,
                 if (length(hit)) paste0(": ", paste(hit, collapse = ", ")) else ""),
         call. = FALSE)
  hit
}

#' assert_partition_labels
#'
#' Partition (TU) labels become path components (`parts/tu=<label>/`,
#' `hist_filter_tu=<label>_Np=...qs`, `<prefix>_<label>[_set-<aset>].yml`),
#' so they must be path- and glob-safe and unambiguous inside those names.
assert_partition_labels <- function(labels, source = "") {
  labels <- unique(as.character(labels))
  bad <- labels[!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", labels) |
                  labels == "_NA_TU_" |
                  grepl("_Np=", labels, fixed = TRUE) |
                  grepl("_set-", labels, fixed = TRUE)]
  if (length(bad))
    stop(sprintf(paste0("assert_partition_labels: partition registry '%s' has label(s) ",
                        "unusable as TU keys: %s (allowed: ^[A-Za-z0-9][A-Za-z0-9_-]*$, ",
                        "not '_NA_TU_', no '_Np=' or '_set-')."),
                 source, paste(bad, collapse = ", ")), call. = FALSE)
  invisible(labels)
}


#' resolve_group_countries
#'
#' When `opt$grouping$enable` is TRUE, derive the country subset for the
#' configured `group` from the registry and write it into `opt$country` so all
#' existing data-subset logic (load_surv_data / readSizeData) is reused
#' unchanged. Grouping wins over an explicit `country` (with a warning).
#'
#' @param opt resolved options list
#' @param groupings_dir directory holding registries
#' @return opt, possibly with `country` overwritten
resolve_group_countries <- function(opt, groupings_dir = "analysis/groupings") {
  if (!isTRUE(opt$grouping$enable)) return(opt)
  src <- opt$grouping$source
  grp <- opt$grouping$group
  if (is.null(src) || is.null(grp) || !nzchar(src) || !nzchar(grp))
    stop("grouping.enable=TRUE requires grouping.source and grouping.group")
  reg  <- load_grouping(src, groupings_dir)
  isos <- reg$iso3[reg$group == grp]
  if (!length(isos))
    stop("resolve_group_countries: group '", grp,
         "' not found in registry '", src, "'")
  if (!is.null(opt$country) && !identical(opt$country, "all"))
    warning("grouping.enable=TRUE: overriding country='", opt$country,
            "' with group '", grp, "' membership")
  opt$country <- paste(isos, collapse = ",")
  cat("-- Grouping '", src, "' / group '", grp, "': country = ",
      opt$country, "\n", sep = "")
  opt
}


#' make_options_from_config
#'
#' Parse CLI options, optionally load a YAML config, then fill any remaining
#' keys from get_default_options(). Priority (high → low):
#'   CLI explicit flag  >  YAML config value  >  default
#'
#' The model blocks (occurrence / size) are merged recursively so that a
#' partial YAML inherits all unspecified keys from the defaults.
#'
#' @param option_list    optparse option list for the calling script
#' @param enforce_options  character vector of flat CLI keys that override YAML
#' @param cli_argv       character vector of CLI args (for testability; default NULL)
#' @param model          "occurrence" | "size" | NULL (auto-detect from YAML)
#' @return fully-resolved named list of options
make_options_from_config <- function(option_list,
                                      enforce_options = NULL,
                                      cli_argv        = NULL,
                                      model           = NULL) {

  # Inject standard flags shared by all scripts
  caller_flags <- vapply(option_list, function(o) o@long_flag, character(1))
  shared_flags <- list()
  if (!"--config" %in% caller_flags) {
    shared_flags <- c(shared_flags, list(optparse::make_option(c("-C", "--config"),
      default = NULL, type = "character",
      help    = "Path to YAML config file. CLI flags override config values. [default: none]")))
  }
  if (!"--update_config" %in% caller_flags) {
    shared_flags <- c(shared_flags, list(optparse::make_option(c("--update_config"),
      default = FALSE, action = "store_true",
      help    = "Re-write the YAML config with the resolved options after parsing. [default: FALSE]")))
  }
  full_option_list <- c(option_list, shared_flags)

  # Parse CLI
  parser  <- optparse::OptionParser(option_list = full_option_list)
  opt_cli <- if (!is.null(cli_argv)) {
    optparse::parse_args(parser, args = cli_argv)
  } else {
    optparse::parse_args(parser)
  }

  # Load YAML if --config provided, else start from CLI
  if (!is.null(opt_cli$config)) {
    if (!file.exists(opt_cli$config)) stop("Config file not found: ", opt_cli$config)
    opt <- yaml::read_yaml(opt_cli$config)
  } else {
    opt <- opt_cli
  }

  # Apply enforce_options: CLI values win over YAML for these flat keys
  # (top-level shared keys like country, seed, etc.)
  if (!is.null(enforce_options)) {
    cat("-- Enforcing CLI options:", paste(enforce_options, collapse = ", "), "\n")
    for (nm in enforce_options) {
      if (!is.null(opt_cli[[nm]])) opt[[nm]] <- opt_cli[[nm]]
    }
  }

  # Merge any CLI-only flat keys not present in the YAML
  if (!is.null(opt_cli$config)) {
    for (nm in names(opt_cli)) {
      if (nm %in% c("config", "update_config")) next
      if (is.null(opt[[nm]])) opt[[nm]] <- opt_cli[[nm]]
    }
  }

  # Auto-detect or apply model hint
  detected_model <- if (!is.null(opt$occurrence)) "occurrence"
                     else if (!is.null(opt$size))     "size"
                     else if (!is.null(opt$duration)) "duration"
                     else if (!is.null(opt$monthly))  "monthly"
                     else model   # fall back to hint

  # Fill defaults
  defaults <- get_default_options()

  # Top-level shared keys
  for (nm in c("baseline", "country", "TL", "TR", "seed", "outdir", "overwrite")) {
    if (is.null(opt[[nm]])) {
      if (!nm %in% c("config", "update_config"))
        cat("-- Setting default:", nm, "=", format(defaults[[nm]]), "\n")
      opt[[nm]] <- defaults[[nm]]
    }
  }

  # Shared top-level blocks, driven by .CONFIG_TOPLEVEL_BLOCKS so that adding a
  # block never requires editing this loop.  Behaviour is identical to the
  # hand-written lines this replaced: enso/grouping are flat-merged and always
  # materialised, simulator is recursively merged only when present.
  for (.blk in names(.CONFIG_TOPLEVEL_BLOCKS)) {
    .spec <- .CONFIG_TOPLEVEL_BLOCKS[[.blk]]
    if (is.null(opt[[.blk]]) && !isTRUE(.spec$fill_when_absent)) next
    opt[[.blk]] <- if (identical(.spec$merge, "recursive")) {
      .merge_recursive(opt[[.blk]], defaults[[.blk]])
    } else {
      .merge1(opt[[.blk]], defaults[[.blk]])
    }
  }

  # Inherit top-level seed if simulator.seed is null
  if (!is.null(opt$simulator) && is.null(opt$simulator$seed)) {
    opt$simulator$seed <- opt$seed
  }

  # Model block
  if (!is.null(detected_model)) {
    opt[[detected_model]] <- .fill_model_block(opt[[detected_model]],
                                                defaults[[detected_model]])
  }

  # Route legacy CLI flat keys into the nested model block.
  # Only when a config file was loaded AND the model block was detected
  # (so pure-CLI invocation still routes correctly).
  if (!is.null(detected_model)) {
    be <- opt[[detected_model]]$backend
    opt <- .route_cli_to_nested(opt, opt_cli, detected_model, be, enforce_options)
  }

  # Grouping: derive country subset from the registry (overrides `country`).
  # Done after model-block fill so the subset is in place for downstream use.
  opt <- resolve_group_countries(opt)

  # Results<-deploy coupling: when results.deploy_config points at a deploy YAML,
  # derive the overlapping results.inputs.* leaves from it (explicit-over-
  # inherited precedence). The logic lives in utils/deploy_core.R, which sources
  # THIS module (config.R) and so cannot be sourced back here (modules never
  # source each other). A caller that sets results.deploy_config must therefore
  # also source utils/deploy_core.R; if it hasn't, fail loudly rather than
  # silently ignoring the pointer (the "dormant branch / missing source" trap).
  if (!is.null(opt$results) && nzchar(opt$results$deploy_config %||% "")) {
    if (!exists("resolve_results_inputs_from_deploy", mode = "function"))
      stop("results.deploy_config is set (", opt$results$deploy_config,
           ") but resolve_results_inputs_from_deploy() is not available -- ",
           "source analysis/utils/deploy_core.R before make_options_from_config().")
    # Thread the known config path (opt_cli$config) so explicit-key detection
    # reads the config actually being resolved, not whatever commandArgs() the
    # resolver would otherwise fall back to. Matters when make_options_from_config
    # is called programmatically (tests) where commandArgs() is empty -- without
    # this the "explicit-wins" guard silently disengages and deploy-derived leaves
    # overwrite the explicit ones.
    opt <- resolve_results_inputs_from_deploy(opt, results_config_path = opt_cli$config)
  }

  # A results config reads every observed panel under ONE episode definition:
  # results.inputs.episode_key (stamped, deploy-derived) is copied onto the
  # top-level episodes.key that makeSurvDatFile() & co. key off. A top-level
  # episodes.key written in the results YAML must agree with it. Applies only
  # to a YAML that STATES results.inputs.episode_key (the manuscript results
  # config must restate every results leaf, test_results_config_coverage.R):
  # opt$results is materialised with defaults for every config, and batch-
  # generated fit configs carry a dumped results: block, so neither may have
  # its own episodes.key overridden by a default.
  raw_cfg <- if (!is.null(opt_cli$config) && file.exists(opt_cli$config))
    yaml::read_yaml(opt_cli$config) else NULL
  if (!is.null(raw_cfg$results$inputs$episode_key)) {
    ek <- opt$results$inputs$episode_key
    raw_top <- (raw_cfg$episodes %||% list())$key
    if (!is.null(raw_top) && !identical(raw_top, ek))
      stop("results config ", opt_cli$config, ": top-level episodes.key = '", raw_top,
           "' disagrees with results.inputs.episode_key = '", ek,
           "'. Drop the top-level key; the results phase takes it from results.inputs.",
           call. = FALSE)
    # Same value set as paths.R::EPISODE_KEYS (not sourced by every caller).
    if (!ek %in% c("lp_season", "gid_season"))
      stop("results.inputs.episode_key = '", ek, "' is not one of {lp_season, gid_season}.",
           call. = FALSE)
    # national zero-fill: same rule, results.inputs.episode_national_zero_fill ->
    # top-level episodes.national_zero_fill.
    nzf <- isTRUE(opt$results$inputs$episode_national_zero_fill)
    raw_nzf <- (raw_cfg$episodes %||% list())$national_zero_fill
    if (!is.null(raw_nzf) && !identical(isTRUE(raw_nzf), nzf))
      stop("results config ", opt_cli$config, ": top-level episodes.national_zero_fill = ",
           tolower(isTRUE(raw_nzf)), " disagrees with results.inputs.episode_national_zero_fill = ",
           tolower(nzf), ". Drop the top-level key; the results phase takes it from ",
           "results.inputs.", call. = FALSE)
    if (nzf && !identical(ek, "gid_season"))
      stop("results.inputs.episode_national_zero_fill = true requires ",
           "results.inputs.episode_key = 'gid_season'.", call. = FALSE)
    nzt <- as.integer(opt$results$inputs$episode_national_zero_threshold %||% 0L)
    # Typed back, so a new stamp records the integer (a quoted "60" in an
    # older YAML would otherwise be stamped as a string).
    opt$results$inputs$episode_national_zero_threshold <- nzt
    raw_nzt <- (raw_cfg$episodes %||% list())$national_zero_threshold
    if (!is.null(raw_nzt) && !identical(as.integer(raw_nzt), nzt))
      stop("results config ", opt_cli$config, ": top-level episodes.national_zero_threshold = ",
           raw_nzt, " disagrees with results.inputs.episode_national_zero_threshold = ",
           nzt, ". Drop the top-level key; the results phase takes it from ",
           "results.inputs.", call. = FALSE)
    if (nzt > 0L && !nzf)
      stop("results.inputs.episode_national_zero_threshold = ", nzt, " requires ",
           "results.inputs.episode_national_zero_fill = true.", call. = FALSE)
    # national-zero rule: same rule again.
    nzr <- as.character(opt$results$inputs[["episode_national_zero_rule"]] %||% "")
    if (identical(nzr, "legacy")) nzr <- ""
    raw_nzr <- (raw_cfg$episodes %||% list())[["national_zero_rule"]]
    if (!is.null(raw_nzr) && !identical(as.character(raw_nzr), nzr))
      stop("results config ", opt_cli$config, ": top-level episodes.national_zero_rule = '",
           raw_nzr, "' disagrees with results.inputs.episode_national_zero_rule = '", nzr,
           "'. Drop the top-level key; the results phase takes it from results.inputs.",
           call. = FALSE)
    if (nzchar(nzr) && !nzf)
      stop("results.inputs.episode_national_zero_rule = '", nzr, "' requires ",
           "results.inputs.episode_national_zero_fill = true.", call. = FALSE)
    # at-risk start / ADM1 overlap: same rule again.
    ars <- as.character(opt$results$inputs[["episode_at_risk_start"]] %||% "season")
    a1o <- as.character(opt$results$inputs[["episode_adm1_overlap"]] %||% "drop")
    for (kv in list(c("at_risk_start", ars, "season"), c("adm1_overlap", a1o, "drop"))) {
      raw_v <- (raw_cfg$episodes %||% list())[[kv[[1]]]]
      if (!is.null(raw_v) && !identical(as.character(raw_v), kv[[2]]))
        stop("results config ", opt_cli$config, ": top-level episodes.", kv[[1]], " = '", raw_v,
             "' disagrees with results.inputs.episode_", kv[[1]], " = '", kv[[2]],
             "'. Drop the top-level key; the results phase takes it from results.inputs.",
             call. = FALSE)
      if (!identical(kv[[2]], kv[[3]]) && !identical(ek, "gid_season"))
        stop("results.inputs.episode_", kv[[1]], " = '", kv[[2]], "' requires ",
             "results.inputs.episode_key = 'gid_season'.", call. = FALSE)
    }
    opt$episodes <- utils::modifyList(opt$episodes %||% list(),
                                      list(key = ek, national_zero_fill = nzf,
                                           national_zero_threshold = nzt,
                                           national_zero_rule = nzr,
                                           at_risk_start = ars, adm1_overlap = a1o))
  }

  # Optionally write back the resolved config
  if (isTRUE(opt_cli$update_config) && !is.null(opt_cli$config)) {
    bak <- paste0(opt_cli$config, ".bak")
    file.copy(opt_cli$config, bak, overwrite = TRUE)
    yaml::write_yaml(opt, opt_cli$config)
    cat("-- Config updated:", opt_cli$config, "(backup:", bak, ")\n")
  }

  # Validate
  .validate_basehaz(opt)

  opt
}


# CLI routing registry: flat CLI key → c(backend_sub_block, key_in_sub_block),
# or NULL for a model-level key that is not inside a backend sub-block.
#
# The MCMC keys are identical for every brms-family (model, backend) pair, so
# they live in one fragment stamped with the sub-block name by .route_fragment().
# Hand-maintaining three parallel tables let the occurrence/brms `sampling` slot
# go missing, which silently dropped `--sampling` and ran the config's draw
# count instead; test_cli_route_completeness.R is the gate against a recurrence.
.MCMC_ROUTE_KEYS <- c(
  sampling        = "sampling",
  warmup          = "warmup",
  chains          = "chains",
  cores           = "cores",
  refresh         = "refresh",
  brms_init       = "init",
  re_uncorrelated = "re_uncorrelated",
  re_init_r       = "re_init_r"
)

#' Build one (model, backend) route entry.
#'
#' @param keys       named character vector: flat CLI key -> key in sub-block
#' @param sub_block  backend sub-block name the keys live under
#' @param extra      extra route entries appended verbatim (incl. NULL paths)
#' @return named list of c(sub_block, key) paths, plus `extra`
.route_fragment <- function(keys, sub_block, extra = list()) {
  paths <- lapply(unname(keys), function(k) c(sub_block, k))
  names(paths) <- names(keys)
  c(paths, extra)
}

.CLI_ROUTE <- list(
  occurrence = list(
    brms = .route_fragment(
      .MCMC_ROUTE_KEYS, "brms",
      list(basehaz_poisson = c("brms", "basehaz_poisson"),
           admin_level     = NULL)),
    brms_cb = .route_fragment(
      .MCMC_ROUTE_KEYS, "brms_cb",
      list(basehaz_poisson = c("brms_cb", "basehaz_poisson"),
           admin_level     = NULL)),
    brms_kl = .route_fragment(
      .MCMC_ROUTE_KEYS, "brms_kl",
      list(basehaz_poisson = c("brms_kl", "basehaz_poisson"),
           admin_level     = NULL))
  ),
  size = list(
    brms    = .route_fragment(.MCMC_ROUTE_KEYS, "brms",
                              list(admin_level = NULL)),
    brms_cb = .route_fragment(.MCMC_ROUTE_KEYS, "brms_cb",
                              list(admin_level = NULL))
  ),
  # `basehaz_duration` / `basehaz_poisson` are model-level on the duration track
  # (see .validate_basehaz), so they route with a NULL path like admin_level.
  duration = list(
    brms = .route_fragment(
      .MCMC_ROUTE_KEYS, "brms",
      list(basehaz_duration = NULL,
           basehaz_poisson  = NULL,
           admin_level      = NULL)),
    brms_cb = .route_fragment(
      .MCMC_ROUTE_KEYS, "brms_cb",
      list(basehaz_duration = NULL,
           basehaz_poisson  = NULL,
           admin_level      = NULL))
  ),
  # Structural copy of `size` — monthly has no basehaz/model-specific route keys.
  monthly = list(
    brms    = .route_fragment(.MCMC_ROUTE_KEYS, "brms",
                              list(admin_level = NULL)),
    brms_cb = .route_fragment(.MCMC_ROUTE_KEYS, "brms_cb",
                              list(admin_level = NULL))
  )
)

# Retained names — callers and tests may reference the per-track tables.
.OCC_CLI_ROUTE  <- .CLI_ROUTE$occurrence
.SIZE_CLI_ROUTE <- .CLI_ROUTE$size
.DUR_CLI_ROUTE  <- .CLI_ROUTE$duration
.MON_CLI_ROUTE  <- .CLI_ROUTE$monthly

#' Route flat CLI key-value pairs into the correct nested model block path.
.route_cli_to_nested <- function(opt, opt_cli, model, backend, enforce_options) {
  route <- .CLI_ROUTE[[model]][[backend]]
  if (is.null(route)) return(opt)

  keys_to_route <- if (!is.null(enforce_options)) {
    intersect(enforce_options, names(route))
  } else {
    names(route)
  }

  for (flat_key in keys_to_route) {
    cli_val <- opt_cli[[flat_key]]
    if (is.null(cli_val)) next
    path <- route[[flat_key]]
    if (is.null(path)) {
      # model-level key (not backend sub-block), e.g. admin_level
      opt[[model]][[flat_key]] <- cli_val
    } else {
      opt[[model]][[path[1]]][[path[2]]] <- cli_val
    }
  }
  opt
}


#' resolve_brms_mcmc
#'
#' Resolve brms iteration counts from a backend sub-block.
#'
#' @details
#' brms takes `warmup` plus a *total* `iter`, while the configs specify the
#' post-warmup draw count as `sampling`, so `iter = warmup + sampling`. Configs
#' predating that convention set `iter` directly; there `sampling` is recovered
#' as `iter - warmup`. That subtraction is why the key cannot be renamed
#' textually: reading a legacy `iter: 2000, warmup: 1000` as `sampling: 2000`
#' would double the retained draws. The legacy branch announces itself so such
#' a config is visible in the run log.
#'
#' @param brms_cfg         the `brms` / `brms_cb` sub-block of the model block
#' @param label            backend name used as the log-message prefix
#' @param default_warmup   warmup when the block sets none
#' @param default_sampling post-warmup draws when the block sets neither
#'                         `sampling` nor `iter`
#' @return list(warmup, sampling, iter, refresh)
resolve_brms_mcmc <- function(brms_cfg, label = "brms",
                              default_warmup   = 1000L,
                              default_sampling = 1000L) {
  warmup_val <- as.integer(brms_cfg$warmup %||% default_warmup)
  if (!is.null(brms_cfg$sampling)) {
    sampling_val <- as.integer(brms_cfg$sampling)
  } else if (!is.null(brms_cfg$iter)) {
    sampling_val <- as.integer(brms_cfg$iter) - warmup_val
    message(sprintf("%s: 'iter' key found (legacy config); deriving sampling = %d - %d = %d",
                    label, as.integer(brms_cfg$iter), warmup_val, sampling_val))
  } else {
    sampling_val <- as.integer(default_sampling)
  }
  iter_val <- warmup_val + sampling_val
  message(sprintf("%s MCMC: warmup = %d  sampling = %d  total iter = %d",
                  label, warmup_val, sampling_val, iter_val))
  list(warmup   = warmup_val,
       sampling = sampling_val,
       iter     = iter_val,
       refresh  = brms_cfg$refresh %||% max(sampling_val %/% 10L, 1L))
}


#' resolve_brms_init
#'
#' brms `init` argument from a backend block's `init` key: a number is passed
#' through, "random" / "0" map to brms's values, anything else warns and falls
#' back to "random". Lifted from the occurrence brms workers.
#' @return "random" or a number
resolve_brms_init <- function(brms_cfg) {
  brms_init_str <- brms_cfg$init %||% "random"
  if (is.numeric(brms_init_str)) {
    brms_init_str
  } else switch(brms_init_str,
    random = "random",
    `0`    = 0,
    {
      warning("Unknown brms_init value '", brms_init_str, "' — using random.")
      "random"
    }
  )
}

#' resolve_brms_threads
#'
#' brms `threads` argument: brms::threading(n) when the backend block sets
#' threads_per_chain > 1, else NULL (no within-chain threading).
resolve_brms_threads <- function(brms_cfg) {
  threads_per_chain <- brms_cfg$threads_per_chain %||% NULL
  if (!is.null(threads_per_chain) && as.integer(threads_per_chain) > 1L) {
    message(sprintf("Multithreading enabled: %d threads per chain", as.integer(threads_per_chain)))
    brms::threading(as.integer(threads_per_chain))
  } else {
    NULL
  }
}


#' print_options
#'
#' Pretty-print the resolved options list at the top of a run log.
#'
#' @param opt  resolved options list from make_options_from_config()
print_options <- function(opt) {
  cat("-------------- Run options ---------------\n")
  # Identify which model block is present
  model <- if (!is.null(opt$occurrence))   "occurrence"
            else if (!is.null(opt$size))     "size"
            else if (!is.null(opt$duration)) "duration"
            else NULL

  simple_keys <- setdiff(names(opt),
                         c("occurrence", "size", "duration", "enso", "simulator"))
  for (nm in simple_keys) {
    if (nm %in% c("config", "update_config")) next
    cat(sprintf("  %-25s %s\n", nm, format(opt[[nm]])))
  }
  if (!is.null(opt$enso))
    cat(sprintf("  %-25s lag=%s source=%s\n", "enso",
                opt$enso$lag_months, opt$enso$source))

  if (!is.null(opt$simulator)) {
    sim <- opt$simulator
    cat(sprintf("  %-25s run_tag=%s\n", "simulator", sim$run_tag %||% "(auto)"))
    cat(sprintf("  %-25s occ=%s\n",      "  fit (occurrence)",   sim$occurrence_fit %||% "(none)"))
    cat(sprintf("  %-25s dur=%s\n",      "  fit (duration)",     sim$duration_fit   %||% "(none)"))
    cat(sprintf("  %-25s mon=%s\n",      "  fit (monthly)",      sim$monthly_fit    %||% "(none)"))
    cat(sprintf("  %-25s %s\n",          "  climate_panel",      sim$climate_panel  %||% "(none)"))
    lf <- sim$location_frame %||% list()
    cat(sprintf("  %-25s hmm=%s pop=%s\n",
                "  location_frame",
                lf$hmm_clusters_path %||% "(none)",
                lf$size_data_path    %||% "(none)"))
    cat(sprintf("  %-25s %s -> %s\n",    "  horizon",            sim$horizon$start, sim$horizon$end))
    cat(sprintf("  %-25s R=%s draws/r=%s strategy=%s seed=%s\n",
                "  replicates",
                sim$n_replicates, sim$draws_per_replicate,
                sim$posterior_strategy, sim$seed %||% "(top-level)"))
    cat(sprintf("  %-25s TUs=%s ctry=%s adm=%s unseen=%s\n",
                "  scope",
                paste(sim$scope$transmission_units, collapse = ","),
                paste(sim$scope$countries,          collapse = ","),
                sim$scope$admin_level, sim$scope$unseen_adm1))
    cat(sprintf("  %-25s climate=%s enso=%s\n",
                "  scenario",
                sim$scenario$climate, sim$scenario$enso))
    cat(sprintf("  %-25s aggregate=%s format=%s\n",
                "  output",
                sim$output$aggregate, sim$output$format))
    cat(sprintf("  %-25s backend=%s workers=%s chunk_size=%s\n",
                "  parallel",
                sim$parallel$backend, sim$parallel$workers, sim$parallel$chunk_size))
  }

  if (!is.null(model)) {
    block <- opt[[model]]
    be    <- block$backend
    cat(sprintf("  %-25s %s (backend: %s)\n", model, "", be))
    # Covariates
    active_covs <- names(Filter(isTRUE, block$covariates %||% list()))
    cat(sprintf("  %-25s %s\n", "  covariates (active)",
                paste(active_covs, collapse = ", ")))
    # REs
    re_flat <- block$random_effects %||% list()
    re_flat$slopes <- NULL
    active_res <- names(Filter(isTRUE, re_flat))
    cat(sprintf("  %-25s %s\n", "  random_effects (active)",
                paste(active_res, collapse = ", ")))
    # Backend sub-block
    be_block <- block[[be]] %||% list()
    for (k in names(be_block)) {
      cat(sprintf("  %-25s %s\n", paste0("  ", model, ".", be, ".", k),
                  format(be_block[[k]])))
    }
  }
  cat("------------------------------------------\n")
}

#' Which model block a raw config YAML carries: "occurrence" | "size" |
#' "duration" | "monthly" (first match in that order), else stop().
#'
#' Single owner of the raw-YAML block test the figure scripts each
#' re-implemented (06_fit_figures_cb*.R `.detect_model()`, the promoted
#' 06_fit_figures_cb_grouped*.R). Reads the file as written, before
#' make_options_from_config() merges defaults in, so default blocks never
#' count as "present".
#'
#' @param cfg_path Path to a config YAML.
detect_model_from_yaml <- function(cfg_path) {
  if (is.null(cfg_path) || !file.exists(cfg_path))
    stop("--config / -C is required and must point to an existing YAML.", call. = FALSE)
  raw <- yaml::read_yaml(cfg_path)
  if (!is.null(raw$occurrence))     "occurrence"
  else if (!is.null(raw$size))      "size"
  else if (!is.null(raw$duration))  "duration"
  else if (!is.null(raw$monthly))   "monthly"
  else stop("Config must have an occurrence:, size:, duration: or monthly: block: ", cfg_path,
            call. = FALSE)
}
