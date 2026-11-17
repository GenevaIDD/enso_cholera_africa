#!/usr/bin/env Rscript
# make_forecast_configs.R — per-unit factual + counterfactual forecast configs,
#   each run to the longest horizon the ENSO forecast and deployed lags allow
#
# Usage:
#   Rscript analysis/utils/deploy_prod_run.R -C <D.yml> --stage gen_forecast
#   Rscript analysis/utils/make_forecast_configs.R \
#     --factual_template <tpl.yml> --counterfactual_template <cf tpl.yml> \
#     [--fits_dir <dir>]
#
# Details:
#   Inputs (read-only): the deployed best-fit directory (--fits_dir), the
#   enso_grouping registry, the cached CFSv2 metadata (--cfsv2_meta) and the
#   factual + counterfactual forecast templates. Outputs: one factual + one
#   counterfactual config per unit (--config_outdir) and a manifest (group ->
#   horizon + merged dirs) read by 07c_simulator_postprocess.R
#   --combine-groups. Normally run by the deploy stage, which supplies every
#   run-specific argument; generate_forecast_configs() is the in-process entry.
#
# See: the "Pass 1" section for the horizon arithmetic

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(lubridate)
})

source("analysis/utils/config.R") # load_grouping
source("analysis/utils/base.R") # %||%
source("analysis/utils/paths.R") # make_run_file (sim_* routing), read_deployed_lag()
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R") # resolve_base_rate_paths()/anchor_spec() -- base_rate_table_paths emission
source("analysis/utils/sim_io.R") # assert_arm_launch() -- per-arm launch (handoff) month
source("analysis/utils/simulator_core.R") # sim_pop_year_spec() -- the arms' population year (manifest provenance)

# The factual/counterfactual TEMPLATES have NO default. smoothing_input is
# inherited VERBATIM from them (the generator never rewrites it), so a stale
# default silently pairs the wrong historical-smoothing cache with a defaults
# regeneration (the factual + set-all arms then fail
# validate_spatpomp_options()). The run is named
# only in analysis/configs/deploy/prod_run<N>.yml (deploy.templates.*); the
# standalone CLI below refuses to run without both --*_template flags.

option_list <- list(
  make_option(
    "--fits_dir",
    type = "character",
    default = "generated_data/enso_best_lag_fits",
    help = "Deployed best-fit artifact dir (== grouping.fits_dir); lag read from the sl<N> filename tag"
  ),
  make_option(
    "--grouping",
    type = "character",
    default = "enso_grouping",
    help = "Grouping registry name (analysis/groupings/<name>.csv)"
  ),
  make_option(
    "--cfsv2_meta",
    type = "character",
    default = "generated_data/enso_cfsv2_download_meta.yml",
    help = "Cached CFSv2 metadata YAML (uses max(target_range))"
  ),
  make_option(
    "--factual_template",
    type = "character",
    default = NULL,
    help = "Factual forecast template (REQUIRED; deploy.templates.forecast)"
  ),
  make_option(
    "--counterfactual_template",
    type = "character",
    default = NULL,
    help = "Counterfactual template (REQUIRED; deploy.templates.counterfactual)"
  ),
  make_option(
    "--config_outdir",
    type = "character",
    default = "analysis/configs/forecast_country_by_group",
    help = "Directory to write the generated per-group configs"
  ),
  make_option(
    "--manifest",
    type = "character",
    default = "analysis/configs/forecast_horizons.yml",
    help = "Path to write the horizons/merged-dir manifest"
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "generated_data/",
    help = "Simulator output root (for computing per-group merged dirs)"
  ),
  make_option(
    "--episode_key",
    type = "character",
    default = "",
    help = paste(
      "Episode definition written into every generated config",
      "(lp_season | gid_season). '' (default) inherits the",
      "template untouched; deploy_prod_run.R passes",
      "deploy.episodes.key. See apply_episode_key_to_sim_config()."
    )
  ),
  make_option(
    "--partition_by",
    type = "character",
    default = "",
    help = paste(
      "Forward-run top-level partition: '' (default, inherit",
      "template; country-partitioned) or 'enso_region' to run",
      "each ENSO group as one region TU (ADM1 contagion couples",
      "across borders; handoff IC assembled from the per-country",
      "smoothing caches). scope.countries/smoothing_input are",
      "unchanged — only filter.partition_by is written."
    )
  ),
  make_option(
    "--partition_source",
    type = "character",
    default = "",
    help = paste(
      "Partition registry (analysis/groupings/<name>.csv) for",
      "--partition_by registry: one config per registry group",
      "(unit), possibly mixing fit groups; unit horizon = max",
      "member firm horizon; v2 manifest with units/countries."
    )
  ),
  make_option(
    "--nbr_partition_source",
    type = "character",
    default = NULL,
    help = "Neighbour-contagion partition registry [default: --partition_source]"
  ),
  make_option(
    "--analog_sets",
    type = "character",
    default = "all",
    help = paste(
      "Comma-separated Lenssen analog sets to emit a counterfactual",
      "config for (all|nina|neutral|elnino). One config per",
      "group x set; the FACTUAL arm is emitted once per group and",
      "is never set-tokenised. [default: %default]"
    )
  ),
  make_option(
    "--primary_analog_set",
    type = "character",
    default = "all",
    help = paste(
      "Which set the flat manifest keys counterfactual_config /",
      "counterfactual_merged_dir alias, so existing consumers",
      "(s22/s18) keep working unchanged [default: %default]"
    )
  ),
  make_option(
    "--prune_stale",
    action = "store_true",
    default = FALSE,
    help = paste(
      "Delete pre-existing counterfactual_*.yml in --config_outdir",
      "that this run does not regenerate. Untokenised leftovers",
      "otherwise still match the SLURM launcher's glob and would",
      "submit phantom jobs."
    )
  )
)
# ---------------------------------------------------------------------------
# generate_forecast_configs(opt) — library entry point
# ---------------------------------------------------------------------------
# The whole generation body lives in this function so the deploy orchestrator
# (analysis/utils/deploy_prod_run.R) can call it IN-PROCESS with an `opt` built from
# the deploy block, instead of shelling out to Rscript. The thin CLI wrapper at
# the bottom (guarded by sys.nframe() == 0L) preserves the standalone
# `Rscript analysis/utils/make_forecast_configs.R ...` entrypoint byte-for-byte.
# The per-run helpers below (merged_dir_for, write_config, patch_common) are
# defined inside so merged_dir_for closes over opt$outdir exactly as before.
# The deployed lag of each track comes from utils/paths.R::read_deployed_lag().
generate_forecast_configs <- function(opt) {
  analog_sets <- trimws(strsplit(opt$analog_sets, ",")[[1L]])
  analog_sets <- analog_sets[nzchar(analog_sets)]
  if (!length(analog_sets)) {
    stop("--analog_sets resolved to nothing.")
  }
  bad_sets <- setdiff(analog_sets, ENSO_ANALOG_SETS)
  if (length(bad_sets)) {
    stop(
      "--analog_sets contains unknown set(s): ",
      paste(bad_sets, collapse = ", "),
      ". Valid: ",
      paste(ENSO_ANALOG_SETS, collapse = ", ")
    )
  }
  if (!opt$primary_analog_set %in% analog_sets) {
    stop(
      "--primary_analog_set '",
      opt$primary_analog_set,
      "' is not among --analog_sets (",
      paste(analog_sets, collapse = ", "),
      "). The flat manifest keys must alias a set ",
      "that actually gets generated."
    )
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  # Merged dir where `07b --aggregate` will write for a given run_tag/run_id —
  # computed exactly as 07b does so the combine step reads from the same place.
  # These are spatPomp-backend runs, so 07b uses the `spatpomp_merged_dir` layout
  # (generated_data/sim_spatpomp/<run_tag>/merged, flat — run_id is ignored), NOT
  # the v1 `sim_merged_dir` layout. Using the wrong kind made the manifest point at
  # empty generated_data/sim/<run_tag>/merged/<run_id> dirs and broke --combine-groups.
  merged_dir_for <- function(run_tag, run_id, fmt = "qs") {
    o <- list(
      outdir = opt$outdir,
      simulator = list(
        run_tag = run_tag,
        run_id = run_id,
        output = list(format = fmt)
      )
    )
    make_run_file(o, kind = "spatpomp_merged_dir")
  }

  # Configs are written with provenance.R::write_yaml_with_header().
  write_config <- write_yaml_with_header

  # ---------------------------------------------------------------------------
  # Resolve inputs
  # ---------------------------------------------------------------------------

  reg <- load_grouping(opt$grouping)
  groups <- sort(unique(reg$group))

  if (!file.exists(opt$cfsv2_meta)) {
    stop(
      "CFSv2 metadata not found: ",
      opt$cfsv2_meta,
      " (run the simulator's ENSO fetch first, or pass --cfsv2_meta)."
    )
  }
  meta <- yaml::read_yaml(opt$cfsv2_meta)
  if (is.null(meta$target_range)) {
    stop("`target_range` absent from ", opt$cfsv2_meta)
  }
  last_cfsv2_month <- max(as.Date(unlist(meta$target_range)))
  message(sprintf(
    "[make_forecast_configs] last CFSv2 forecast month: %s",
    format(last_cfsv2_month)
  ))
  # Firm horizons are tied to the ENSO seam: the last month
  # whose ONI owes nothing to the continuation (enso.R::enso_seam_month()). The
  # default seam "last_cfsv2" is the legacy rule (firm ONI end = last CFSv2 month);
  # "last_full_oni" moves every firm horizon back one month.
  firm_seam <- as.character(
    (opt$enso_continuation %||% list())[["seam"]] %||% "last_cfsv2"
  )
  firm_oni_end <- enso_seam_month(last_cfsv2_month, firm_seam)
  firm_legacy <- identical(firm_seam, "last_cfsv2")
  # Only the analog continuation has a seam state (enso.R::sample_enso_cfsv2_trajectories()
  # refuses another seam otherwise): refuse it here, at config time.
  if (
    !firm_legacy &&
      !identical(
        (opt$enso_continuation %||% list())[["method"]] %||% "neutral",
        "analog"
      )
  ) {
    stop(
      "[make_forecast_configs] enso_continuation.seam '",
      firm_seam,
      "' needs enso_continuation.method ",
      "'analog' (only the analog continuation has a seam state).",
      call. = FALSE
    )
  }
  if (!firm_legacy) {
    message(sprintf(
      "[make_forecast_configs] seam %s: firm ONI end %s (last CFSv2 - 1 month).",
      firm_seam,
      format(firm_oni_end)
    ))
  }

  factual_tpl <- yaml::read_yaml(opt$factual_template)
  counterfactual_tpl <- yaml::read_yaml(opt$counterfactual_template)
  # Episode definition the deployed fits are looked up under (the "_ep-gs"
  # filename filter of read_deployed_lag()): the deploy's key when given, else
  # the factual template's own episodes.key.
  ep_key <- if (nzchar(opt$episode_key %||% "")) {
    opt$episode_key
  } else {
    (factual_tpl$episodes %||% list())$key %||% "lp_season"
  }
  # national zero-fill occurrence fits ("_nz"): the template's own declaration
  # (gen_templates writes deploy.episodes.national_zero_fill into it; the deploy
  # re-checks it in deploy_check_sim_template()).
  ep_nz <- isTRUE((factual_tpl$episodes %||% list())$national_zero_fill)
  ep_nzt <- .episodes_nz_threshold(factual_tpl$episodes)
  # The whole occurrence-panel spec (utils/paths.R, with its variant tokens).
  ep_spec <- .occ_panel_spec_from_episodes(utils::modifyList(
    factual_tpl$episodes %||% list(),
    list(key = ep_key)
  ))
  base_tag_fac <- factual_tpl$simulator$run_tag
  base_tag_cf <- counterfactual_tpl$simulator$run_tag
  # Population year (simulator.location_frame.pop_year): both
  # arms must resolve the same one (gen_templates refuses otherwise); recorded
  # in the manifest only when not the default.
  pop_year_arms <- list(
    factual = sim_pop_year_spec(factual_tpl$simulator),
    counterfactual = sim_pop_year_spec(counterfactual_tpl$simulator)
  )
  if (!identical(pop_year_arms$factual, pop_year_arms$counterfactual)) {
    stop(
      sprintf(
        paste0(
          "[make_forecast_configs] the factual and counterfactual templates resolve ",
          "different simulator.location_frame.pop_year (%s vs %s); the arms must ",
          "share one population year."
        ),
        format(pop_year_arms$factual),
        format(pop_year_arms$counterfactual)
      ),
      call. = FALSE
    )
  }

  # Each arm's launch (= its smoothing handoff month; patch_common() never
  # rewrites horizon.start, so every generated config of an arm shares it).
  launch <- c(
    factual = assert_arm_launch(factual_tpl$simulator, opt$factual_template),
    counterfactual = assert_arm_launch(
      counterfactual_tpl$simulator,
      opt$counterfactual_template
    )
  )
  launch_chr <- vapply(launch, format, character(1))
  if (launch[["counterfactual"]] > launch[["factual"]]) {
    stop(
      sprintf(
        paste0(
          "[make_forecast_configs] the counterfactual template launches at %s, after ",
          "the factual template (%s). Only an earlier (or equal) counterfactual launch ",
          "is supported: the results readers stitch the factual series from the ",
          "smoothing back to the counterfactual launch, never the reverse."
        ),
        launch_chr[["counterfactual"]],
        launch_chr[["factual"]]
      ),
      call. = FALSE
    )
  }
  if (launch[["counterfactual"]] < launch[["factual"]]) {
    message(sprintf(
      paste0(
        "[make_forecast_configs] arms launch at different months: factual %s, ",
        "counterfactual %s; attribution window starts at %s (factual series ",
        "stitched from the smoothing output before its launch)."
      ),
      launch_chr[["factual"]],
      launch_chr[["counterfactual"]],
      launch_chr[["counterfactual"]]
    ))
  }

  date_str <- format(Sys.Date(), "%Y%m%d")
  out_fmt <- factual_tpl$simulator$output$format %||% "qs"

  dir.create(opt$config_outdir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Per-group config generation
  # ---------------------------------------------------------------------------
  # Handoff: each arm launches at its template's simulator.horizon.start, which
  # must equal filter.forecast_window.start and (when set) filter.handoff_month
  # (sim_io.R::assert_arm_launch(), the 07b invariant). The two arms may launch at
  # different months from the same smoothing cache (the counterfactual at the
  # Lenssen analog anchor, the factual later); the manifest records each arm's
  # `launch` and the `attribution_start` (the earlier launch), which the results
  # readers use to stitch the factual series from the smoothing output up to its
  # own launch. A counterfactual launching after the factual arm is refused.

  manifest_groups <- list()
  manifest_units <- list()
  manifest_countries <- list()
  horizons <- c()
  skipped <- c()
  generated_cf_paths <- character(0)

  patch_common <- function(
    cfg,
    run_tag,
    run_id,
    isos,
    horizon_end_chr,
    base_rate_table_paths = NULL
  ) {
    cfg$simulator$run_tag <- run_tag
    cfg$simulator$run_id <- run_id
    cfg$simulator$scope$countries <- paste(isos, collapse = ",")
    cfg$simulator$horizon$end <- horizon_end_chr
    cfg$simulator$filter$forecast_window$end <- horizon_end_chr
    # Optional region partition. Only written when --partition_by is given, so the
    # default (country-partitioned) output stays byte-for-byte identical. Under
    # enso_region the forward run couples member-country ADM1 units across borders
    # via the contagion kernel; scope.countries (the region member list) and
    # smoothing_input (per-country caches) are left untouched — the region handoff
    # IC is assembled from those same per-country caches by 07b.
    # %||% "" keeps the library entry point NULL-safe for any caller whose opt
    # carries no partition_by key (deploy_prod_run.R's .forecast_opt() always
    # sets one): nzchar(NULL) is logical(0) and `if (logical(0))` errors.
    if (nzchar(opt$partition_by %||% "")) {
      cfg$simulator$filter$partition_by <- opt$partition_by
    }
    # Registry partition: the partition registry, and the neighbour-contagion
    # partition registry (default: the same registry, so contagion blocks refine
    # the TUs -- spatpomp_core.R::check_partition_vs_grouping()).
    if (identical(opt$partition_by %||% "", "registry")) {
      cfg$simulator$filter$partition_source <- opt$partition_source
      cfg$simulator$nbr_contagion$partition_source <-
        opt$nbr_partition_source %||% opt$partition_source
    }
    # ENSO continuation past the CFSv2 span (utils/enso.R); written only when set.
    if (
      !is.null(opt$enso_continuation) &&
        !identical(opt$enso_continuation$method %||% "neutral", "neutral")
    ) {
      cont <- opt$enso_continuation
      # The analog seam is written only when not the default, so
      # a deploy that states seam: last_cfsv2 regenerates its configs identically.
      if (
        identical(as.character(cont[["seam"]] %||% "last_cfsv2"), "last_cfsv2")
      ) {
        cont[["seam"]] <- NULL
      }
      cfg$simulator$enso$continuation <- cont
    }
    # The occ/dur/mon anchors (utils/anchor.R), resolved below at generation
    # time from each deployed fit's own training opt. NULL whenever
    # no deployed fit carries `<track>.offset_baserate: true`; the key is then
    # absent, so generated configs do not change until a track is fitted with
    # the anchor.
    if (!is.null(base_rate_table_paths) && length(base_rate_table_paths)) {
      cfg$simulator$location_frame$base_rate_table_paths <- base_rate_table_paths
    }
    # Episode definition: the deploy's key, written into the config and checked
    # against the prep artefacts it will read (climate panel, size frame -- still
    # the location frame's population / admin_level source -- survival panel)
    # -- now, not at SLURM run time.
    cfg <- apply_episode_key_to_sim_config(
      cfg,
      opt$episode_key %||% "",
      label = run_tag
    )
    preflight_sim_episode_inputs(cfg, label = run_tag)
    cfg
  }

  # ---- Pass 1: per FIT GROUP -- deployed lags, firm horizon, anchors ---------
  # Horizon: the covariate at month M uses ONI at M - L, so with ENSO known to
  # firm_oni_end a group forecasts to firm_oni_end + its min lag (occ/dur/mon):
  #   firm_oni_end       = enso.R::enso_seam_month(last CFSv2 month, seam)
  #   horizon_end_g      = firm_oni_end %m+% months(min_lag_g)  (closed)
  #   africa_min_horizon = min_g horizon_end_g  (full-Africa cutoff)
  #   projection_cal_end = max_g horizon_end_g  (projection panel must cover it)
  # Lags come from the DEPLOYED fits in --fits_dir
  # (paths.R::read_deployed_lag()), not the raw LOO ranking: the deploy's
  # deploy_select_lags() may pick another lag.
  group_info <- list()
  for (g in groups) {
    # Deployed occ/dur/mon fits of this group (utils/paths.R::read_deployed_lag();
    # NULL when the track has no deployed fit for the group).
    bl <- lapply(stats::setNames(nm = c("occ", "dur", "mon")), function(tr) {
      read_deployed_lag(
        opt$fits_dir,
        g,
        opt$grouping,
        tr,
        episode_key = ep_key,
        national_zero_fill = ep_nz,
        national_zero_threshold = ep_nzt,
        occ_panel = ep_spec
      )
    })
    missing_tr <- names(bl)[vapply(bl, is.null, logical(1))]
    if (length(missing_tr)) {
      message(sprintf(
        "[make_forecast_configs] SKIP group '%s' (no deployed %s fit)",
        g,
        paste(missing_tr, collapse = "/")
      ))
      skipped <- c(skipped, g)
      next
    }

    lags <- vapply(bl, `[[`, integer(1), "lag")
    min_lag <- min(lags)
    horizon_end <- firm_oni_end %m+% months(min_lag)
    horizon_chr <- format(horizon_end)
    # A deployed ENSO lead (negative lag; deploy.lag_selection) moves the firm
    # horizon BEFORE the firm ONI end. The horizon is CLOSED and the launch month
    # is simulated (utils/base.R), so a firm horizon equal to an
    # arm's launch still gives it one firm month; only a firm horizon before the
    # launch leaves no month with forecast ENSO: refuse it (the later launch
    # binds; it is the factual arm whenever the launches differ).
    late_arm <- names(launch)[which.max(launch)]
    if (!length(horizon_month_seq(max(launch), horizon_end))) {
      stop(
        sprintf(
          paste0(
            "[make_forecast_configs] group '%s': firm horizon %s (%s %s + min lag %d; ",
            "lags occ/dur/mon = %s) is before the %s arm's launch %s. Lower the leads ",
            "this deploy may select (deploy.lag_selection) or drop the group."
          ),
          g,
          horizon_chr,
          if (firm_legacy) "last CFSv2" else "firm ONI end",
          format(firm_oni_end),
          min_lag,
          paste(lags, collapse = "/"),
          late_arm,
          launch_chr[[late_arm]]
        ),
        call. = FALSE
      )
    }
    isos <- reg$iso3[reg$group == g]

    # base_rate_table_paths: resolve HERE, at config-generation time, from each
    # deployed fit's own training opt, so a missing sidecar stops this script
    # instead of a SLURM array. resolve_base_rate_paths() does the anchor_enabled()
    # gate and the existence check, and returns NULL when no deployed fit was
    # trained with the anchor.
    base_rate_table_paths <- resolve_base_rate_paths(
      lapply(bl, function(b) attr(readRDS(b$file), "opt"))
    )

    group_info[[g]] <- list(
      bl = bl,
      lags = lags,
      min_lag = min_lag,
      horizon_chr = horizon_chr,
      isos = isos,
      base_rate_table_paths = base_rate_table_paths
    )
  }

  # ---- Units: the simulation partition --------------------------------------
  # Without a partition registry a unit IS a fit group (byte-identical configs,
  # tags and manifest to before). Under partition_by: registry a unit is a group
  # of opt$partition_source: its members are the registry's countries whose fit
  # group has deployed fits, it may mix fit groups (each country keeps its own
  # group's draws and lag) and it is simulated to the MAX member firm horizon;
  # readers censor each country at its own firm horizon (manifest `countries`).
  registry_mode <- identical(opt$partition_by %||% "", "registry")
  units <- if (!registry_mode) {
    lapply(stats::setNames(nm = names(group_info)), function(g) {
      list(members = group_info[[g]]$isos, groups = g)
    })
  } else {
    if (!nzchar(opt$partition_source %||% "")) {
      stop(
        "[make_forecast_configs] --partition_by registry requires --partition_source.",
        call. = FALSE
      )
    }
    preg <- load_grouping(opt$partition_source)
    assert_partition_labels(preg$group, opt$partition_source)
    fit_of <- stats::setNames(reg$group, reg$iso3)
    uncov <- setdiff(reg$iso3[reg$group %in% names(group_info)], preg$iso3)
    if (length(uncov)) {
      stop(
        sprintf(
          paste0(
            "[make_forecast_configs] partition registry '%s' does not cover ",
            "%d country(ies) with deployed fits: %s"
          ),
          opt$partition_source,
          length(uncov),
          paste(uncov, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    u <- list()
    for (lab in sort(unique(preg$group))) {
      mem <- preg$iso3[preg$group == lab]
      nofit <- mem[!mem %in% names(fit_of)]
      if (length(nofit)) {
        message(sprintf(
          "[make_forecast_configs] unit '%s': %s not in fit registry '%s' -- dropped",
          lab,
          paste(nofit, collapse = ","),
          opt$grouping
        ))
      }
      mem <- mem[mem %in% names(fit_of)]
      skip <- mem[!fit_of[mem] %in% names(group_info)]
      if (length(skip)) {
        message(sprintf(
          "[make_forecast_configs] unit '%s': %s dropped (fit group without deployed fits)",
          lab,
          paste(skip, collapse = ",")
        ))
      }
      mem <- setdiff(mem, skip)
      if (!length(mem)) {
        next
      }
      u[[lab]] <- list(
        members = mem,
        groups = sort(unique(unname(fit_of[mem])))
      )
    }
    u
  }

  for (un in names(units)) {
    U <- units[[un]]
    isos <- U$members
    gi <- group_info[U$groups]
    firm_chr <- vapply(gi, `[[`, character(1), "horizon_chr")
    horizon_chr <- max(firm_chr) # ISO dates sort lexically
    # Anchors are per-track files shared across groups only when identical: one
    # config carries one base_rate_table_paths.
    brt <- unique(lapply(gi, `[[`, "base_rate_table_paths"))
    if (length(brt) > 1L) {
      stop(
        sprintf(
          paste0(
            "[make_forecast_configs] unit '%s' mixes fit groups (%s) whose ",
            "deployed fits use different anchor sidecars; one config carries ",
            "one simulator.location_frame.base_rate_table_paths."
          ),
          un,
          paste(U$groups, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    base_rate_table_paths <- brt[[1L]]

    # run tags / ids (unit-suffixed → isolate parts/merged dirs)
    tag_fac <- paste0(base_tag_fac, "_", un)
    rid_fac <- paste0(date_str, "_", un, "_forecast")

    lag_lines <- vapply(
      names(gi),
      function(g) {
        x <- gi[[g]]
        sprintf(
          "# best_lag_occ=%d, best_lag_dur=%d, best_lag_mon=%d -> min_lag=%d; horizon.end=%s (%s %s + %d mo).",
          x$lags[["occ"]],
          x$lags[["dur"]],
          x$lags[["mon"]],
          x$min_lag,
          x$horizon_chr,
          if (firm_legacy) "last CFSv2" else "firm ONI end, last CFSv2 - 1",
          format(firm_oni_end),
          x$min_lag
        )
      },
      character(1)
    )
    hdr <- function(kind, extra = character(0)) {
      generated_yaml_header(
        "analysis/utils/make_forecast_configs.R",
        c(
          if (!registry_mode) {
            sprintf(
              "# %s forecast for ENSO group '%s' (%d countries).",
              kind,
              un,
              length(isos)
            )
          } else {
            c(
              sprintf(
                "# %s forecast for partition unit '%s' of registry '%s' (%d countries; fit groups %s).",
                kind,
                un,
                opt$partition_source,
                length(isos),
                paste(U$groups, collapse = ", ")
              ),
              sprintf(
                "# Unit horizon.end = max member firm horizon = %s; each country is censored at its",
                horizon_chr
              ),
              "# own fit group's firm horizon by the readers (manifest `countries`)."
            )
          },
          if (!registry_mode) {
            lag_lines
          } else {
            paste0(sprintf("# [%s] ", names(gi)), sub("^# ", "", lag_lines))
          },
          extra,
          sprintf(
            "# Launch (smoothing handoff) %s: horizon.start = forecast_window.start = handoff month.",
            launch_chr[[
              if (identical(kind, "Factual")) "factual" else "counterfactual"
            ]]
          )
        )
      )
    }

    # -- Factual arm: ONE config per unit, deliberately NOT set-tokenised.
    # The factual arm draws from the CFSv2 forecast, not from the Lenssen analog pool,
    # so it has no analog set. Tokenising it would force either N identical redundant
    # runs or N aliases of a single output tree; instead every counterfactual set is
    # contrasted against the same factual run.
    cfg_fac <- patch_common(
      factual_tpl,
      tag_fac,
      rid_fac,
      isos,
      horizon_chr,
      base_rate_table_paths
    )
    path_fac <- file.path(opt$config_outdir, paste0("forecast_", un, ".yml"))
    assert_basename_length(c(path_fac, tag_fac), "make_forecast_configs")
    write_config(cfg_fac, path_fac, hdr("Factual"))

    # -- Counterfactual arm: one config per analog set.
    cf_sets <- list()
    for (aset in analog_sets) {
      tag_cf <- paste0(base_tag_cf, "_", un, "_set-", aset)
      rid_cf <- paste0(date_str, "_", un, "_counterfactual_set-", aset)
      path_cf <- file.path(
        opt$config_outdir,
        paste0("counterfactual_", un, "_set-", aset, ".yml")
      )

      cfg_cf <- patch_common(
        counterfactual_tpl,
        tag_cf,
        rid_cf,
        isos,
        horizon_chr,
        base_rate_table_paths
      )
      cfg_cf$simulator$enso$counterfactual$analog_set <- aset

      write_config(
        cfg_cf,
        path_cf,
        hdr(
          sprintf("Counterfactual (analog_set=%s)", aset),
          sprintf(
            "# Lenssen analog set '%s'; MC precision for ENSO-driven quantities is bounded by the SET size.",
            aset
          )
        )
      )

      cf_sets[[aset]] <- list(
        analog_set = aset,
        run_tag = tag_cf,
        run_id = rid_cf,
        config = path_cf,
        merged_dir = merged_dir_for(tag_cf, rid_cf, out_fmt),
        launch = launch_chr[["counterfactual"]]
      )
      generated_cf_paths <- c(generated_cf_paths, path_cf)
    }

    primary <- cf_sets[[opt$primary_analog_set]]
    run_fields <- list(
      factual_config = path_fac,
      factual_merged_dir = merged_dir_for(tag_fac, rid_fac, out_fmt),
      # Per-arm launch (= smoothing handoff) months; read by
      # utils/sim_io.R::manifest_arm_launches().
      factual_launch = launch_chr[["factual"]],
      analog_sets = as.list(analog_sets),
      primary_analog_set = opt$primary_analog_set,
      counterfactual_sets = cf_sets,
      # Flat aliases for the primary set. KEEP THESE: the results-phase run-tag
      # resolver reads `counterfactual_config` / `counterfactual_merged_dir`
      # directly, and the default (`all`) invocation stays byte-identical.
      counterfactual_config = primary$config,
      counterfactual_merged_dir = primary$merged_dir,
      counterfactual_launch = primary$launch
    )

    if (!registry_mode) {
      g <- un
      x <- gi[[g]]
      manifest_groups[[g]] <- c(
        list(
          best_lag_occ = x$lags[["occ"]],
          best_lag_dur = x$lags[["dur"]],
          best_lag_mon = x$lags[["mon"]],
          min_lag = x$min_lag,
          horizon_end = x$horizon_chr,
          n_countries = length(isos),
          iso3 = as.list(isos),
          occ_fit_file = x$bl$occ$file,
          dur_fit_file = x$bl$dur$file,
          mon_fit_file = x$bl$mon$file
        ),
        run_fields
      )
      horizons[g] <- x$horizon_chr
      message(sprintf(
        "[make_forecast_configs] %-12s occ=%d dur=%d mon=%d min=%d horizon_end=%s  sets: %s",
        g,
        x$lags[["occ"]],
        x$lags[["dur"]],
        x$lags[["mon"]],
        x$min_lag,
        x$horizon_chr,
        paste(analog_sets, collapse = ",")
      ))
    } else {
      manifest_units[[un]] <- c(
        list(
          members = as.list(isos),
          groups = as.list(U$groups),
          horizon_end = horizon_chr
        ),
        run_fields
      )
      message(sprintf(
        "[make_forecast_configs] unit %-12s groups=%s members=%s horizon_end=%s (firm: %s)  sets: %s",
        un,
        paste(U$groups, collapse = ","),
        paste(isos, collapse = ","),
        horizon_chr,
        paste(sprintf("%s=%s", names(firm_chr), firm_chr), collapse = ","),
        paste(analog_sets, collapse = ",")
      ))
    }
  }

  # Registry mode: per fit group (no run fields -- runs are per unit), per country.
  if (registry_mode) {
    for (g in names(group_info)) {
      x <- group_info[[g]]
      manifest_groups[[g]] <- list(
        best_lag_occ = x$lags[["occ"]],
        best_lag_dur = x$lags[["dur"]],
        best_lag_mon = x$lags[["mon"]],
        min_lag = x$min_lag,
        horizon_end = x$horizon_chr,
        n_countries = length(x$isos),
        iso3 = as.list(x$isos),
        occ_fit_file = x$bl$occ$file,
        dur_fit_file = x$bl$dur$file,
        mon_fit_file = x$bl$mon$file,
        units = as.list(names(units)[vapply(
          units,
          function(U) g %in% U$groups,
          logical(1)
        )])
      )
      horizons[g] <- x$horizon_chr
    }
    for (un in names(units)) {
      for (iso in units[[un]]$members) {
        g <- reg$group[reg$iso3 == iso]
        manifest_countries[[iso]] <- list(
          group = g,
          unit = un,
          firm_horizon_end = group_info[[g]]$horizon_chr,
          unit_horizon_end = manifest_units[[un]]$horizon_end
        )
      }
    }
  }

  if (!length(manifest_groups)) {
    stop(
      "No group has a deployed occurrence, duration AND monthly fit in ",
      opt$fits_dir,
      " (episodes.key = ",
      ep_key,
      ") -- nothing generated."
    )
  }

  horizon_dates <- as.Date(unname(horizons))
  africa_min_horizon <- min(horizon_dates)
  projection_cal_end <- max(horizon_dates)

  combined_tag_fac <- paste0(base_tag_fac, "_combined")

  # One combined counterfactual tree per analog set (`_set-all` for all sets).
  combined_cf <- lapply(stats::setNames(nm = analog_sets), function(aset) {
    tag <- paste0(base_tag_cf, "_combined_set-", aset)
    list(
      analog_set = aset,
      run_tag = tag,
      merged_dir = merged_dir_for(tag, NULL, out_fmt)
    )
  })

  manifest <- list(
    generated_at = format(Sys.time()),
    grouping = opt$grouping,
    last_cfsv2_month = format(last_cfsv2_month),
    africa_min_horizon = format(africa_min_horizon),
    projection_cal_end = format(projection_cal_end),
    skipped_groups = as.list(skipped),
    analog_sets = as.list(analog_sets),
    primary_analog_set = opt$primary_analog_set,
    # Attribution window start: the earlier arm launch (sim_io.R::manifest_attrib_start()).
    attribution_start = format(min(launch)),
    combined = list(
      factual = merged_dir_for(combined_tag_fac, NULL, out_fmt),
      # Flat alias for the primary set (see the per-group note above).
      counterfactual = combined_cf[[opt$primary_analog_set]]$merged_dir,
      counterfactual_sets = combined_cf
    ),
    groups = manifest_groups
  )
  # Schema v2 (registry mode only; additive -- a v1 reader still finds `groups`
  # with each fit group's firm horizon). utils/sim_io.R::manifest_*() read both.
  if (registry_mode) {
    unit_h <- as.Date(vapply(manifest_units, `[[`, character(1), "horizon_end"))
    manifest <- c(list(schema = 2L), manifest)
    manifest$grouping_digest <- registry_digest(opt$grouping)
    manifest$partition_by <- "registry"
    manifest$partition_source <- opt$partition_source
    manifest$partition_digest <- registry_digest(opt$partition_source)
    manifest$contagion_partition_source <- opt$nbr_partition_source %||%
      opt$partition_source
    # method only for neutral / ar1 (as before); analog adds the library path,
    # its content digest and the k-NN knobs (enso.R::enso_continuation_provenance()).
    manifest$enso_continuation <- if (
      identical((opt$enso_continuation %||% list())[["method"]], "analog")
    ) {
      enso_continuation_provenance(opt$enso_continuation)
    } else {
      opt$enso_continuation[intersect(
        names(opt$enso_continuation),
        c("method", "phi")
      )] %||%
        list(method = "neutral")
    }
    manifest$africa_max_horizon <- format(max(unit_h))
    manifest$projection_cal_end <- format(max(
      as.Date(manifest$projection_cal_end),
      max(unit_h)
    ))
    manifest$units <- manifest_units
    manifest$countries <- manifest_countries
  }
  # Firm-horizon rule: recorded only when not the legacy rule, so
  # a legacy run regenerates identically (provenance.R::manifest_add_provenance();
  # readers: sim_io.R::manifest_last_firm_oni()).
  manifest <- manifest_add_provenance(
    manifest,
    "firm_oni_end",
    format(firm_oni_end),
    format(last_cfsv2_month)
  )
  manifest <- manifest_add_provenance(
    manifest,
    "firm_horizon_rule",
    paste0(firm_seam, " + min lag"),
    "last_cfsv2 + min lag"
  )
  # Closed horizons: the arms' climate panels must hold the last simulated month
  # (sim_io.R::assert_panel_covers_horizon(); a hard error here, a warning in 07b).
  for (pp in unique(c(
    factual_tpl$simulator$climate_panel,
    counterfactual_tpl$simulator$climate_panel
  ))) {
    assert_panel_covers_horizon(
      pp,
      manifest$projection_cal_end,
      "[make_forecast_configs]"
    )
  }

  # Counterfactual analog cache, library mode (deploy
  # enso_counterfactual.match: library): built or verified HERE, in-process, at the
  # counterfactual template's cache_path, with every DJF the reported months read
  # classified (winters: all_read). The reported months of a group can reach its
  # unit's horizon (the Africa window reads groups on past their firm horizon), so
  # each group's window end is its unit's horizon; read_span "unit" also applies
  # the unit's minimum lag (censored members still feed the unit's contagion).
  # Legacy mode (default): nothing is built and the manifest is unchanged.
  cf_match <- enso_counterfactual_cfg(opt$enso_counterfactual)
  if (!identical(cf_match$match, "legacy")) {
    lcf <- lenssen_counterfactual_args(
      counterfactual_tpl$simulator$enso$counterfactual,
      launch[["counterfactual"]]
    )
    win_rows <- do.call(
      rbind,
      lapply(names(units), function(un) {
        data.frame(
          group = units[[un]]$groups,
          unit = un,
          end = as.Date(
            manifest_units[[un]]$horizon_end %||%
              max(vapply(
                group_info[units[[un]]$groups],
                `[[`,
                character(1),
                "horizon_chr"
              ))
          ),
          stringsAsFactors = FALSE
        )
      })
    )
    cf_read_end <- enso_read_end(
      win_rows,
      lapply(group_info, `[[`, "lags"),
      cf_match$read_span
    )
    cf_panel <- counterfactual_tpl$simulator$climate_panel
    cf_n <- if (!is.null(cf_panel) && file.exists(cf_panel)) {
      pnl <- readRDS(cf_panel)
      lenssen_n_lead_months(
        lcf$horizon_start,
        max(as.Date(pnl$t_start)),
        panel_max_enso_lead(pnl)
      ) +
        1L
    } else {
      message(
        "[make_forecast_configs] counterfactual panel '",
        cf_panel %||% "",
        "' not found here: ",
        "cache length from projection_cal_end + 12 lead months."
      )
      lenssen_n_lead_months(
        lcf$horizon_start,
        manifest$projection_cal_end,
        12L
      ) +
        1L
    }
    cf_cache <- ensure_lenssen_counterfactual_cache(
      lcf$cache_path,
      lcf$horizon_start,
      cf_n,
      cf_match,
      opt$enso_continuation,
      cf_read_end,
      sets = analog_sets,
      primary_set = opt$primary_analog_set,
      oni_dir = opt$outdir %||% "generated_data/",
      label = "[make_forecast_configs]"
    )
    message(sprintf(
      "[make_forecast_configs] counterfactual winters: primary DJF %d, constrained %s (read end %s, span %s)",
      cf_cache$winters$primary,
      if (length(cf_cache$winters$other)) {
        paste(cf_cache$winters$other, collapse = ",")
      } else {
        "none"
      },
      format(cf_read_end),
      cf_match$read_span
    ))
    manifest <- manifest_add_provenance(
      manifest,
      "enso_counterfactual",
      c(
        cf_cache$provenance,
        list(read_end = format(cf_read_end), read_span = cf_match$read_span)
      ),
      NULL
    )
  }

  manifest <- manifest_add_provenance(
    manifest,
    "pop_year",
    pop_year_arms$factual,
    .SIM_POP_YEAR_DEFAULT
  )
  dir.create(dirname(opt$manifest), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(manifest, opt$manifest)

  # Stale counterfactual configs. This script only ever WRITES, so an earlier run's
  # untokenised `counterfactual_<group>.yml` would survive here and still match the
  # SLURM launcher's `counterfactual_*.yml` glob -> phantom job submissions against
  # an output tree nothing else references.
  stale_cf <- setdiff(
    list.files(
      opt$config_outdir,
      pattern = "^counterfactual_.*\\.yml$",
      full.names = TRUE
    ),
    generated_cf_paths
  )
  if (length(stale_cf)) {
    if (isTRUE(opt$prune_stale)) {
      file.remove(stale_cf)
      message(sprintf(
        "[make_forecast_configs] pruned %d stale counterfactual config(s):\n    %s",
        length(stale_cf),
        paste(basename(stale_cf), collapse = "\n    ")
      ))
    } else {
      warning(sprintf(
        paste0(
          "%d stale counterfactual config(s) left in %s and NOT regenerated by this run:\n    %s\n",
          "  These still match the launcher glob in analysis/bash/stages/07_forecast_groups.sh\n",
          "  and would submit jobs for output trees the manifest no longer references.\n",
          "  Re-run with --prune_stale to delete them."
        ),
        length(stale_cf),
        opt$config_outdir,
        paste(basename(stale_cf), collapse = "\n    ")
      ))
    }
  }

  message(
    "\n[make_forecast_configs] ----------------------------------------------"
  )
  message(sprintf(
    "  groups generated : %d  (%s)",
    length(manifest_groups),
    paste(names(manifest_groups), collapse = ", ")
  ))
  message(sprintf(
    "  analog sets      : %s  (primary: %s)",
    paste(analog_sets, collapse = ", "),
    opt$primary_analog_set
  ))
  message(sprintf(
    "  configs written  : %d factual + %d counterfactual (%d %s x %d sets)",
    length(units),
    length(generated_cf_paths),
    length(units),
    if (registry_mode) "units" else "groups",
    length(analog_sets)
  ))
  if (registry_mode) {
    message(sprintf(
      "  units (registry %s): %s",
      opt$partition_source,
      paste(
        sprintf(
          "%s{%s}",
          names(units),
          vapply(
            units,
            function(U) paste(U$members, collapse = ","),
            character(1)
          )
        ),
        collapse = " "
      )
    ))
  }
  if (length(skipped)) {
    message(sprintf("  groups skipped   : %s", paste(skipped, collapse = ", ")))
  }
  message(sprintf("  africa_min_horizon: %s", format(africa_min_horizon)))
  message(sprintf(
    "  projection_cal_end: %s  (rebuild the projection panel to this)",
    manifest$projection_cal_end
  ))
  if (registry_mode) {
    message(sprintf("  africa_max_horizon: %s", manifest$africa_max_horizon))
  }
  message(sprintf("  configs -> %s", opt$config_outdir))
  message(sprintf("  manifest -> %s", opt$manifest))
  message("[make_forecast_configs] done.")

  invisible(manifest)
}

# ---------------------------------------------------------------------------
# Thin CLI wrapper — standalone `Rscript analysis/utils/make_forecast_configs.R`
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  opt <- parse_args(OptionParser(option_list = option_list))
  missing_tpl <- c("--factual_template", "--counterfactual_template")[
    c(is.null(opt$factual_template), is.null(opt$counterfactual_template))
  ]
  if (length(missing_tpl)) {
    stop(
      "[make_forecast_configs] ",
      paste(missing_tpl, collapse = " and "),
      " required (no run-specific default). Normally run via: Rscript ",
      "analysis/utils/deploy_prod_run.R -C analysis/configs/deploy/",
      "prod_run<N>.yml --stage gen_forecast",
      call. = FALSE
    )
  }
  generate_forecast_configs(opt)
}
