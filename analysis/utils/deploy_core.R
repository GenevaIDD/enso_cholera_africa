# deploy_core.R — domain logic and artefact naming of the deploy orchestrator
#   (deploy_prod_run.R); bash stages keep only the SLURM boundary.
#
# Provides:
#   deploy_config_path(), load_deploy_config()   locate + merge a deploy YAML
#   deploy_collect_fits()     select best lags, copy fits, write MANIFEST.csv
#   deploy_check_convergence(), deploy_convergence_cfg()   convergence gate
#   resolve_results_inputs_from_deploy()   results.inputs.* from a deploy
#   compare_deploy_results(), align_results_to_deploy()    results vs deploy
#   deploy_plan_batches(), deploy_check_batch_groupings()  fit batches
#   deploy_plan_sim_templates(), deploy_write_sim_template()   sim templates
#   deploy_resolved_keys(), deploy_emit_keys()   KEY=value paths for bash
#   deploy_require_keys(), deploy_write_provenance()
#
# Details:
#   Sourced by deploy_prod_run.R and (resolve / compare helpers) by
#   09b_postprocess_results.R, 10_manuscript_figures.R and tests. Unlike other
#   modules it sources config.R, base.R, paths.R, provenance.R, loo_lags.R and
#   diagnostics.R itself; none of them depends back on it.
#
# See: analysis/utils/deploy_prod_run.R

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/config.R") # make_options_from_config, .merge_recursive, .flag_val
source("analysis/utils/base.R") # %||%
source("analysis/utils/paths.R") # make_*_run_file family, make_deploy_run_file
source("analysis/utils/provenance.R") # attach_deploy_provenance
source("analysis/utils/loo_lags.R") # select_best_lag_occurrence / _downstream
source("analysis/utils/diagnostics.R") # fit_diagnostics, fit_diagnostics_gate (convergence gate)
source("analysis/utils/simulator_core.R") # sim_pop_year_spec() -- gen_templates' cross-arm population-year check
source("analysis/utils/national_zero.R") # national_zero_rule_params() -- deploy_check_national_zero_rule()

# ---------------------------------------------------------------------------
# Config access
# ---------------------------------------------------------------------------

#' Recover the -C/--config path the deploy orchestrator was invoked with.
#' The path also names the artefact run tag (generated_data/deploy/<stem>/),
#' exactly like results_config_path() in results_core.R.
deploy_config_path <- function(argv = commandArgs(TRUE)) {
  p <- .flag_val("-C", argv) %||% .flag_val("--config", argv)
  if (is.null(p) || !nzchar(p)) {
    stop(
      "[deploy] no -C/--config given: the deploy phase is config-driven and ",
      "the config path also names the artefact directory. Run with ",
      "-C analysis/configs/deploy/prod_run<N>.yml",
      call. = FALSE
    )
  }
  p
}

#' Load and resolve a deploy config's `deploy:` block.
#'
#' Merges the YAML's `deploy:` block recursively over get_default_options()$deploy
#' (so a partial deploy YAML inherits unspecified leaves) WITHOUT the full
#' make_options_from_config() machinery (no model detection / basehaz validation
#' a deploy file would trip over).
#'
#' @param path Path to the deploy YAML.
#' @return The resolved `deploy` block (named list).
load_deploy_config <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("[deploy] deploy config not found: ", path, call. = FALSE)
  }
  raw <- yaml::read_yaml(path)
  defaults <- get_default_options()
  .merge_recursive(raw$deploy, defaults$deploy)
}

# ---------------------------------------------------------------------------
# Best-lag fit collection
# ---------------------------------------------------------------------------

# Per-track registry. `make_fn`/`kinds` single-source the canonical filenames via
# paths.R; `model` drives make_options_from_config()'s block detection;
# `anchor_of` pulls the track's offset_baserate flag for the opt-attr patch gate.
# The deployed tracks are occ/dur/mon: size is retired downstream,
# the monthly track carries the case counts. NO hardcoded filename strings.
.deploy_track_specs <- function() {
  list(
    occ = list(
      track = "occ",
      model = "occurrence",
      prefix = "outbreaks_survival_",
      batch_key = "occ_dir",
      make_fn = make_run_file,
      kinds = c("fit", "re_stats", "loo", "cb_meta"),
      anchor_of = function(o) isTRUE(o$occurrence$offset_baserate)
    ),
    dur = list(
      track = "dur",
      model = "duration",
      prefix = "outbreak_duration_",
      batch_key = "dur_dir",
      make_fn = make_duration_run_file,
      kinds = c("dur_fit", "dur_re_stats", "dur_loo", "dur_cb_meta"),
      anchor_of = function(o) isTRUE(o$duration$offset_baserate)
    ),
    mon = list(
      track = "mon",
      model = "monthly",
      prefix = "outbreak_monthly_",
      batch_key = "mon_dir",
      make_fn = make_monthly_run_file,
      kinds = c("mon_fit", "mon_re_stats", "mon_loo", "mon_cb_meta"),
      anchor_of = function(o) isTRUE(o$monthly$offset_baserate)
    )
  )
}

# Latest ranking cache for a (batch, group): utils/loo_lags.R::
# read_latest_ranking() -- <outdir>/loo_compare/<batch>/ first, legacy
# analysis/configs/<batch>/ fallback; latest by filename date, then mtime.
# `dep$fits$*_dir` holds bare batch names; deploy groups are the short form
# ("east"), ranking filenames carry the full label ("enso_east").
.deploy_ranking <- function(dep, batch_key, group) {
  read_latest_ranking(
    dep$fits[[batch_key]],
    group = .deploy_group_full(dep, group)
  )
}

# Group-label resolution against the deploy's fit registry (dep$grouping):
# deploy.fits.groups may hold short forms ("east") or full labels
# ("enso_east"); the registry sidecar's label_prefix maps one onto the other
# (utils/config.R::load_grouping_meta() / resolve_group_label()).
.deploy_group_meta <- function(dep) {
  load_grouping_meta(dep$grouping %||% "enso_grouping")
}
.deploy_group_full <- function(dep, group) {
  resolve_group_label(.deploy_group_meta(dep), group)
}

# Number of countries of deploy group `group` in the deploy's grouping registry.
.deploy_group_n_countries <- function(dep, group) {
  reg <- load_grouping(dep$grouping)
  sum(reg$group == .deploy_group_full(dep, group))
}

# The brms_kl class sweep that derives `grouping` (deploy.kl_sweep.batch;
# "" = none, the registry is used as is) and its 05_fit_postprocess_kl_sweep.R
# comparison dir, generated_data/kl_sweep_compare/<batch> (paths.R::
# make_kl_sweep_run_file(), the path 05 writes; outdir without a trailing slash so
# the path reads as manuscript.yml's results.kl_sweep_dir).
.deploy_kl_sweep_batch <- function(dep) as.character(dep$kl_sweep$batch %||% "")
deploy_kl_sweep_compare_dir <- function(dep) {
  b <- .deploy_kl_sweep_batch(dep)
  if (!nzchar(b)) {
    return("")
  }
  make_kl_sweep_run_file(
    list(outdir = "generated_data"),
    "kl_sweep_dir",
    batch = b
  )
}

# Batch stem from a batch dir name: "prod_run_occ_grouped_enso_singlelag" ->
# "prod_run". Shared by the lag_loo_batch derivation below.
.deploy_batch_stem <- function(dir_name) {
  sub("_(occ|dur|duration|mon|monthly)_grouped.*$", "", dir_name)
}

# Per-lag config path for (track_spec, group, lag), matching the batch layout:
# analysis/configs/<batch>/<batch>_sl-<lag>_grp-<label token>.yml, `group`
# being the FULL registry label ("enso_east" -> "_grp-ensoeast").
.deploy_cfg_path <- function(batch_dir, group, lag) {
  stem <- basename(batch_dir)
  assert_basename_length(
    file.path(
      batch_dir,
      sprintf(
        "%s_sl-%s_grp-%s.yml",
        stem,
        sl_token(lag),
        grouping_label_token(group)
      )
    ),
    ".deploy_cfg_path"
  )
}

# Artefacts in `files` that belong to the OTHER episode definition than
# `ep_key` -- by the "_ep-gs" filename token (utils/paths.R::
# .episode_key_path_token()) -- or, among occurrence artefacts (`occ_prefix`), to
# the other national zero-fill variant ("_nz", .national_zero_fill_file_match()).
# deploy_collect_fits() refuses to prune these.
.deploy_foreign_episode_artefacts <- function(
  files,
  ep_key,
  national_zero_fill = FALSE,
  occ_prefix = "outbreaks_survival_",
  national_zero_threshold = 0L,
  occ_panel = NULL
) {
  b <- basename(files)
  other_key <- grepl("_ep-gs", b, fixed = TRUE) !=
    identical(ep_key, "gid_season")
  # occ_panel (an occurrence-panel spec, utils/paths.R) replaces the pair and
  # carries every panel-variant token: exact whole-suffix match.
  other_nz <- startsWith(b, occ_prefix) &
    !(if (is.null(occ_panel)) {
      .national_zero_fill_file_match(
        files,
        national_zero_fill,
        national_zero_threshold
      )
    } else {
      .occ_panel_variant_file_match(files, occ_panel)
    })
  files[other_key | other_nz]
}

#' Check a deploy's national-zero rule: the rule exists in
#' analysis/national_zero_rules.yml and the local national-zero months file,
#' when present, was built with exactly the registry's parameters (its 01
#' stamp's rule_params_digest). Silent for the legacy rule ("").
#' @return invisibly, the rule's parameters (NULL for the legacy rule).
deploy_check_national_zero_rule <- function(
  dep,
  say = function(...) invisible(NULL)
) {
  spec <- .occ_panel_spec_from_episodes(dep$episodes)
  rule <- spec$national_zero_rule
  if (!nzchar(rule)) {
    return(invisible(NULL))
  }
  params <- national_zero_rule_params(rule)
  path <- as.character(national_zero_months_file_of(spec))
  if (file.exists(path)) {
    st <- read_prep_policy_stamp(path) %||% list()
    if (
      !identical(st[["rule_params_digest"]], national_zero_rule_digest(params))
    ) {
      stop(
        "[deploy] deploy.episodes.national_zero_rule = '",
        rule,
        "': ",
        path,
        " was built with other parameters (",
        st[["rule_params_digest"]] %||% "no rule stamp",
        ") than ",
        NATIONAL_ZERO_RULES_PATH,
        " (",
        national_zero_rule_digest(params),
        "). Rebuild it with 01_preprocess_national_zero_months.R --rule ",
        rule,
        ".",
        call. = FALSE
      )
    }
  } else {
    say(
      "[deploy] national-zero months %s not present here: rule parameters not checked",
      path
    )
  }
  say(
    "[deploy] national-zero rule %s: %s",
    rule,
    paste(
      sprintf(
        "%s=%s",
        setdiff(names(params), "rule"),
        unlist(params[setdiff(names(params), "rule")])
      ),
      collapse = " "
    )
  )
  invisible(params)
}

#' Select the best (occ, duration, monthly) lag triple per group.
#'
#' Occurrence = max single_lag within 2 SEs of the rank-1 fit
#' (select_best_lag_occurrence); duration/monthly = each track's own
#' unconstrained SE-window pick (select_best_lag_downstream with
#' enforce_occ_floor = FALSE).
#'
#' @param dep Resolved deploy block.
#' @return Named list keyed by group; each element c(occ=, dur=, mon=) int
#'   (names and order from .deploy_track_specs()).
#' @details With deploy.lag_selection.occ_floor_offset set, the floor is
#'   checked against the median whole-month episode span stamped on the
#'   deploy's duration panel (loo_lags.R::check_occ_floor_offset()): a
#'   mismatch warns, a missing stamp only messages; the selection is
#'   unchanged either way.
deploy_select_lags <- function(dep) {
  specs <- .deploy_track_specs()
  out <- list()
  # deploy.lag_selection (loo_lags.R::lag_selection_cfg()); absent / default =
  # the longest_within_se rule with se_mult 2 and min_lag 0, exactly as before.
  ls_cfg <- lag_selection_cfg(dep$lag_selection)
  best <- identical(ls_cfg$rule, "best_elpd")
  if (!is.null(ls_cfg$occ_floor_offset)) {
    check_occ_floor_offset(
      ls_cfg,
      read_prep_policy_stamp(make_duration_dat_file(list(
        baseline = get_default_options()$baseline,
        episodes = list(key = dep$episodes$key %||% "lp_season")
      )))
    )
  }
  for (g in dep$fits$groups) {
    occ_rank <- .deploy_ranking(dep, specs$occ$batch_key, g)
    occ_lag <- if (best) {
      select_best_lag_best_elpd(occ_rank, min_lag = ls_cfg$min_lag)$lag
    } else {
      select_best_lag_occurrence(
        occ_rank,
        se_mult = ls_cfg$se_mult,
        min_lag = ls_cfg$min_lag
      )$lag
    }

    lag_of <- function(sp) {
      if (sp$track == "occ") {
        return(occ_lag)
      }
      rank <- .deploy_ranking(dep, sp$batch_key, g)
      if (best) {
        return(
          select_best_lag_best_elpd(
            rank,
            min_lag = ls_cfg$min_lag,
            occ_lag = occ_lag,
            occ_floor_offset = ls_cfg$occ_floor_offset
          )$lag
        )
      }
      select_best_lag_downstream(
        rank,
        occ_lag,
        enforce_occ_floor = FALSE,
        se_mult = ls_cfg$se_mult,
        min_lag = ls_cfg$min_lag
      )$lag
    }
    out[[g]] <- vapply(specs, function(sp) as.integer(lag_of(sp)), integer(1))
  }
  out
}

# ---------------------------------------------------------------------------
# Convergence gate on the source fits
# ---------------------------------------------------------------------------

.DEPLOY_CONVERGENCE_ACTIONS <- c("off", "warn", "refuse")

# Short form of deploy group `g` ("east" from "east" or "enso_east"): the full
# registry label without the sidecar's label_prefix.
.deploy_group_short <- function(dep, g) {
  full <- .deploy_group_full(dep, g)
  pre <- .deploy_group_meta(dep)$label_prefix %||% ""
  if (nzchar(pre)) sub(paste0("^", pre), "", full) else full
}

# A "default" sentinel (or absent leaf) resolves to `default`; anything else
# must be one finite number.
.deploy_conv_number <- function(x, default, what) {
  if (is.null(x) || identical(as.character(x), "default")) {
    return(unname(default))
  }
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || !is.finite(v)) {
    stop(
      "[deploy] deploy.convergence.",
      what,
      " must be \"default\" or a number (got '",
      paste(x, collapse = ","),
      "').",
      call. = FALSE
    )
  }
  v
}

#' deploy_convergence_cfg
#'
#' Resolve and validate deploy.convergence. A missing block (or leaf) takes
#' the in-function default, which is today's behaviour: action `off` (no fit
#' is read, nothing is printed or written), scope `population`, thresholds
#' "default" (= diagnostics.R::DIAG_THRESHOLDS, never restated here),
#' max_divergent 0, allow empty. `allow` ids are `<track>:<short group>` over
#' .deploy_track_specs() x deploy.fits.groups (e.g. "occ:west"); an unknown id
#' is an error.
#'
#' @param validate_allow FALSE skips the `allow` id check (gen_grouping: the
#'   fit registry the ids refer to may not be written yet).
#' @return list(action, scope, rhat_max, ess_min, max_divergent, allow,
#'   rhat_any, ebfmi_min, ids (every valid id; empty when allow is empty),
#'   raw (the YAML threshold leaves))
deploy_convergence_cfg <- function(dep, validate_allow = TRUE) {
  cv <- dep[["convergence"]] %||% list()
  action <- as.character(cv[["action"]] %||% "off")
  if (length(action) != 1L || !action %in% .DEPLOY_CONVERGENCE_ACTIONS) {
    stop(
      "[deploy] deploy.convergence.action must be one of {",
      paste(.DEPLOY_CONVERGENCE_ACTIONS, collapse = ", "),
      "} (got '",
      paste(action, collapse = ","),
      "').",
      call. = FALSE
    )
  }
  scope <- as.character(cv[["scope"]] %||% "population")
  if (length(scope) != 1L || !scope %in% names(DIAG_SCOPES)) {
    stop(
      "[deploy] deploy.convergence.scope must be one of {",
      paste(names(DIAG_SCOPES), collapse = ", "),
      "} (got '",
      paste(scope, collapse = ","),
      "').",
      call. = FALSE
    )
  }
  rhat_max <- .deploy_conv_number(
    cv[["rhat_max"]],
    DIAG_THRESHOLDS[["rhat"]],
    "rhat_max"
  )
  ess_min <- .deploy_conv_number(
    cv[["ess_min"]],
    DIAG_THRESHOLDS[["ess"]],
    "ess_min"
  )
  md <- suppressWarnings(as.numeric(cv[["max_divergent"]] %||% 0L))
  if (length(md) != 1L || !is.finite(md) || md < 0 || md != round(md)) {
    stop(
      "[deploy] deploy.convergence.max_divergent must be a non-negative integer (got '",
      paste(cv[["max_divergent"]], collapse = ","),
      "').",
      call. = FALSE
    )
  }
  allow <- as.character(unlist(cv[["allow"]] %||% character(0)))
  ids <- character(0)
  if (length(allow) && validate_allow) {
    ids <- as.vector(outer(
      names(.deploy_track_specs()),
      vapply(
        dep$fits$groups,
        function(g) .deploy_group_short(dep, g),
        character(1)
      ),
      paste,
      sep = ":"
    ))
    bad <- setdiff(allow, ids)
    if (length(bad)) {
      stop(
        "[deploy] deploy.convergence.allow names unknown fit(s): ",
        paste(bad, collapse = ", "),
        ". Ids are <track>:<short group> over this deploy's fits: ",
        paste(ids, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  list(
    action = action,
    scope = scope,
    rhat_max = rhat_max,
    ess_min = ess_min,
    max_divergent = as.integer(md),
    allow = unique(allow),
    rhat_any = DIAG_GATE[["rhat_any"]],
    ebfmi_min = DIAG_GATE[["ebfmi"]],
    ids = ids,
    raw = list(
      rhat_max = cv[["rhat_max"]] %||% "default",
      ess_min = cv[["ess_min"]] %||% "default"
    )
  )
}

# One read per existing source fit of `plan` (is_fit rows, in plan order),
# shared by the episode check and the convergence gate: the fit's episode
# record (fit_episode_stamp_of()) when `episode`, its fit_diagnostics() when
# `diagnostics`. Neither -> nothing is read. Both lists are keyed by plan$file.
.deploy_read_source_fits <- function(
  plan,
  episode = FALSE,
  diagnostics = FALSE,
  read_fn = readRDS,
  diag_fn = fit_diagnostics
) {
  out <- list(episode = list(), diagnostics = list())
  if (!episode && !diagnostics) {
    return(out)
  }
  for (i in which(plan$is_fit)) {
    if (!file.exists(plan$src[[i]])) {
      next
    }
    f <- plan$file[[i]]
    fit <- read_fn(plan$src[[i]])
    if (episode) {
      out$episode[[f]] <- fit_episode_stamp_of(fit, f)
    }
    if (diagnostics) {
      out$diagnostics[[f]] <- diag_fn(fit)
    }
    rm(fit)
  }
  out
}

#' deploy_check_convergence
#'
#' The convergence gate of `--stage deploy_fits`: every source fit of `plan`
#' through diagnostics.R::fit_diagnostics_gate() under deploy.convergence
#' (deploy_convergence_cfg()). Run by deploy_collect_fits() BEFORE anything in
#' dest is pruned or copied, so a refusal leaves the current deployment as it
#' is. A failing fit listed in `allow` passes as `allowed`; a failing best lag
#' never falls back to another lag. `off`: returns NULL without reading a fit.
#' `warn`: warning() naming the failing fits. `refuse`: stop() naming them,
#' except under `dry_run`, which prints "WOULD refuse" and never stops (the
#' preview gate). A source fit not produced yet is reported `missing` (the
#' copy stops on it).
#'
#' @param plan  deploy_collect_fits()'s plan (track, group, lag, kind, src, file, is_fit)
#' @param diags precomputed fit_diagnostics() results keyed by plan$file (the
#'   shared read); NULL = read here through read_fn / diag_fn
#' @return NULL (off) or list(cfg, table (one row per fit: track, group, id,
#'   lag, file, the gate columns, verdict, action, allowed), failed, allowed)
deploy_check_convergence <- function(
  plan,
  dep,
  dry_run = FALSE,
  say = function(...) invisible(NULL),
  read_fn = readRDS,
  diag_fn = fit_diagnostics,
  diags = NULL,
  cfg = deploy_convergence_cfg(dep)
) {
  if (identical(cfg$action, "off")) {
    return(NULL)
  }
  fits <- plan[plan$is_fit, , drop = FALSE]
  if (is.null(diags)) {
    diags <- .deploy_read_source_fits(
      fits,
      diagnostics = TRUE,
      read_fn = read_fn,
      diag_fn = diag_fn
    )$diagnostics
  }
  short_of <- stats::setNames(
    vapply(
      dep$fits$groups,
      function(g) .deploy_group_short(dep, g),
      character(1)
    ),
    vapply(
      dep$fits$groups,
      function(g) grouping_label_token(.deploy_group_full(dep, g)),
      character(1)
    )
  )
  rows <- lapply(seq_len(nrow(fits)), function(i) {
    id <- paste(fits$track[[i]], short_of[[fits$group[[i]]]], sep = ":")
    d <- diags[[fits$file[[i]]]]
    g <- if (is.null(d)) {
      NULL
    } else {
      fit_diagnostics_gate(
        d,
        scope = cfg$scope,
        rhat_max = cfg$rhat_max,
        ess_min = cfg$ess_min,
        max_divergent = cfg$max_divergent,
        rhat_any = cfg$rhat_any,
        ebfmi_min = cfg$ebfmi_min
      )
    }
    if (is.null(g)) {
      g <- tibble::tibble(
        scope = cfg$scope,
        n_params_scope = NA_integer_,
        max_rhat_scope = NA_real_,
        max_rhat_variable = NA_character_,
        min_ess_bulk_scope = NA_real_,
        min_ess_tail_scope = NA_real_,
        max_rhat_all = NA_real_,
        divergences = NA_integer_,
        treedepth_hits = NA_integer_,
        min_ebfmi = NA_real_,
        pass = NA,
        reasons = "missing"
      )
    }
    allowed <- id %in% cfg$allow
    verdict <- if (is.na(g$pass)) {
      "missing"
    } else if (g$pass) {
      "pass"
    } else if (allowed) {
      "allowed"
    } else {
      "fail"
    }
    cbind(
      data.frame(
        track = fits$track[[i]],
        group = fits$group[[i]],
        id = id,
        lag = fits$lag[[i]],
        file = fits$file[[i]],
        stringsAsFactors = FALSE
      ),
      as.data.frame(g, stringsAsFactors = FALSE),
      data.frame(
        verdict = verdict,
        action = cfg$action,
        allowed = allowed,
        stringsAsFactors = FALSE
      )
    )
  })
  tab <- do.call(rbind, rows)
  say(
    "[deploy] convergence gate: action=%s scope=%s rhat_max=%s ess_min=%s max_divergent=%d (+ R-hat <= %s on every parameter, E-BFMI >= %s)",
    cfg$action,
    cfg$scope,
    format(cfg$rhat_max),
    format(cfg$ess_min),
    cfg$max_divergent,
    format(cfg$rhat_any),
    format(cfg$ebfmi_min)
  )
  for (i in seq_len(nrow(tab))) {
    if (identical(tab$verdict[[i]], "missing")) {
      say(
        "[deploy]   %-8s %-10s lag=%3s  (source fit not produced yet)",
        "MISSING",
        tab$id[[i]],
        tab$lag[[i]]
      )
    } else {
      say(
        "[deploy]   %-8s %-10s lag=%3s  rhat=%s (%s) ess_bulk=%s ess_tail=%s rhat_all=%s div=%s td=%s ebfmi=%s  %s%s",
        toupper(tab$verdict[[i]]),
        tab$id[[i]],
        tab$lag[[i]],
        formatC(tab$max_rhat_scope[[i]], format = "f", digits = 4),
        tab$max_rhat_variable[[i]] %||% NA,
        formatC(tab$min_ess_bulk_scope[[i]], format = "f", digits = 0),
        formatC(tab$min_ess_tail_scope[[i]], format = "f", digits = 0),
        formatC(tab$max_rhat_all[[i]], format = "f", digits = 4),
        tab$divergences[[i]],
        tab$treedepth_hits[[i]],
        formatC(tab$min_ebfmi[[i]], format = "f", digits = 2),
        tab$reasons[[i]],
        if (tab$allowed[[i]]) " [allow-listed]" else ""
      )
    }
  }
  failed <- tab$id[tab$verdict == "fail"]
  allowed <- tab$id[tab$verdict == "allowed"]
  unused <- setdiff(cfg$allow, allowed)
  if (length(unused)) {
    say(
      "[deploy]   allow-listed but not failing: %s",
      paste(unused, collapse = ", ")
    )
  }
  say(
    "[deploy] convergence: %d pass, %d allowed, %d fail, %d missing (of %d fits)",
    sum(tab$verdict == "pass"),
    length(allowed),
    length(failed),
    sum(tab$verdict == "missing"),
    nrow(tab)
  )
  if (length(failed)) {
    msg <- sprintf(
      paste0(
        "%d deployed fit(s) fail the convergence gate (scope %s): %s. Refit them (longer chains) ",
        "or list them in deploy.convergence.allow, stated in Methods / Table S8."
      ),
      length(failed),
      cfg$scope,
      paste(failed, collapse = ", ")
    )
    if (identical(cfg$action, "refuse")) {
      if (dry_run) {
        say("[deploy] WOULD refuse: %s", msg)
      } else {
        stop(
          "[deploy] refusing to deploy: ",
          msg,
          " Nothing in ",
          dep$fits$dest,
          " was pruned or copied.",
          call. = FALSE
        )
      }
    } else {
      warning("[deploy] ", msg, call. = FALSE)
    }
  }
  list(cfg = cfg, table = tab, failed = failed, allowed = allowed)
}

# Canonical filenames are built by the same paths.R code (make_*_run_file) the
# fits were written with, so tokens such as the occurrence `_nbr500` are never
# re-hardcoded in bash printf templates.
#' Collect the canonical per-group best-lag fit artefacts into dep$fits$dest.
#'
#' @param dep         Resolved deploy block (opt$deploy).
#' @param dry_run     If TRUE, report the planned prune/copy and write nothing.
#' @param verbose     If TRUE, print progress.
#' @param read_fn,diag_fn how a source fit is read and diagnosed (injectable
#'   for tests; the convergence gate, deploy_check_convergence()).
#' @return A list: `selected` (per-group occ/dur/mon lag triples), `manifest`
#'   (data.frame written to MANIFEST.csv), `copied`/`pruned`/`patched` counts,
#'   `dest`, `convergence` (deploy_check_convergence(); NULL under action
#'   off). The orchestrator stamps this into the provenance manifest.
deploy_collect_fits <- function(
  dep,
  dry_run = FALSE,
  verbose = TRUE,
  read_fn = readRDS,
  diag_fn = fit_diagnostics
) {
  specs <- .deploy_track_specs()
  dest <- dep$fits$dest
  say <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n", sep = "")
  # The rollout's episode definition (deploy.episodes.key) -- checked against
  # every per-lag config (below, also under --dry-run), every copied fit's own
  # record and the prep artefacts they were trained on (after the copy).
  ep_key <- dep$episodes$key %||% "lp_season"
  if (!ep_key %in% EPISODE_KEYS) {
    stop(
      "[deploy] deploy.episodes.key = '",
      ep_key,
      "' is not one of {",
      paste(EPISODE_KEYS, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
  say("[deploy] episodes.key = %s", ep_key)
  ep_nz <- isTRUE(dep$episodes$national_zero_fill)
  if (ep_nz && !identical(ep_key, "gid_season")) {
    stop(
      "[deploy] deploy.episodes.national_zero_fill = true requires deploy.episodes.key = ",
      "'gid_season'.",
      call. = FALSE
    )
  }
  ep_nzt <- .episodes_nz_threshold(dep$episodes)
  if (ep_nzt > 0L && !ep_nz) {
    stop(
      "[deploy] deploy.episodes.national_zero_threshold = ",
      ep_nzt,
      " requires ",
      "deploy.episodes.national_zero_fill = true.",
      call. = FALSE
    )
  }
  if (ep_nz) {
    say(
      "[deploy] episodes.national_zero_fill = true (occurrence fits on the _nz%s panel)",
      if (ep_nzt > 0L) ep_nzt else ""
    )
  }
  # The whole occurrence-panel spec: its keys are read in-function with their
  # defaults when the deploy block does not state them.
  ep_spec <- .occ_panel_spec_from_episodes(dep$episodes)
  if (
    nzchar(ep_spec$national_zero_rule) ||
      !identical(ep_spec$at_risk_start, "season") ||
      !identical(ep_spec$adm1_overlap, "drop")
  ) {
    say(
      "[deploy] occurrence panel variant %s",
      .occ_panel_variant_token(ep_spec)
    )
  }
  deploy_check_national_zero_rule(dep, say = say)

  ls_cfg <- lag_selection_cfg(dep$lag_selection)
  say(
    "[deploy] lag selection: rule=%s min_lag=%s se_mult=%s occ_floor_offset=%s",
    ls_cfg$rule,
    format(ls_cfg$min_lag),
    format(ls_cfg$se_mult),
    format(ls_cfg$occ_floor_offset %||% "none")
  )
  selected <- deploy_select_lags(dep)
  for (g in names(selected)) {
    say(
      "[deploy] group=%-6s occ=%d dur=%d mon=%d",
      g,
      selected[[g]][["occ"]],
      selected[[g]][["dur"]],
      selected[[g]][["mon"]]
    )
  }

  # Build the expected artefact set: for each (track, group), load the per-lag
  # config -> opt -> canonical source paths via make_*_run_file. `_nbr500` emerges
  # automatically for occurrence (paths.R reads the config's contagion block).
  plan_rows <- list() # one row per (track, group, kind)
  fit_opts <- list() # opt + anchor flag keyed by "<track>:<group>" for the fit patch
  for (sp in specs) {
    for (g in dep$fits$groups) {
      lag <- selected[[g]][[sp$track]]
      cfg <- .deploy_cfg_path(
        file.path("analysis/configs", dep$fits[[sp$batch_key]]),
        .deploy_group_full(dep, g),
        lag
      )
      if (!file.exists(cfg)) {
        stop(
          "[deploy] per-lag config not found for ",
          sp$track,
          " group '",
          g,
          "' lag ",
          lag,
          ": ",
          cfg,
          call. = FALSE
        )
      }
      opt <- make_options_from_config(
        list(),
        cli_argv = c("-C", cfg),
        model = sp$model
      )
      if (!identical(episode_key_of(opt), ep_key)) {
        stop(
          "[deploy] per-lag config ",
          cfg,
          " declares episodes.key = '",
          episode_key_of(opt),
          "' but deploy.episodes.key = '",
          ep_key,
          "'. Point deploy.fits.",
          sp$batch_key,
          " at the batch fitted under ",
          "the deploy's episode definition.",
          call. = FALSE
        )
      }
      # national zero-fill: only the occurrence fits are estimated on that panel.
      if (
        identical(sp$track, "occ") &&
          !identical(national_zero_fill_of(opt), ep_nz)
      ) {
        stop(
          "[deploy] per-lag config ",
          cfg,
          " declares episodes.national_zero_fill = ",
          tolower(national_zero_fill_of(opt)),
          " but deploy.episodes.national_zero_fill = ",
          tolower(ep_nz),
          ". Point deploy.fits.",
          sp$batch_key,
          " at the batch ",
          "fitted on the matching occurrence panel.",
          call. = FALSE
        )
      }
      if (
        identical(sp$track, "occ") &&
          !identical(national_zero_threshold_of(opt), ep_nzt)
      ) {
        stop(
          "[deploy] per-lag config ",
          cfg,
          " declares episodes.national_zero_threshold = ",
          national_zero_threshold_of(opt),
          " but deploy.episodes.national_zero_threshold = ",
          ep_nzt,
          ". Point deploy.fits.",
          sp$batch_key,
          " at the batch ",
          "fitted on the matching occurrence panel.",
          call. = FALSE
        )
      }
      if (identical(sp$track, "occ")) {
        o_spec <- occ_panel_spec_of(opt)
        for (f in c("national_zero_rule", "at_risk_start", "adm1_overlap")) {
          if (!identical(o_spec[[f]], ep_spec[[f]])) {
            stop(
              "[deploy] per-lag config ",
              cfg,
              " declares episodes.",
              f,
              " = '",
              o_spec[[f]],
              "' but deploy.episodes.",
              f,
              " = '",
              ep_spec[[f]],
              "'. Point deploy.fits.",
              sp$batch_key,
              " at the batch ",
              "fitted on the matching occurrence panel.",
              call. = FALSE
            )
          }
        }
      }
      # Keyed on the same group token as plan$group (full label, tokenised), so
      # the opt-attr patch below finds it whether the deploy lists "east" or
      # "enso_east", under any registry.
      fit_opts[[paste(
        sp$track,
        grouping_label_token(.deploy_group_full(dep, g)),
        sep = ":"
      )]] <-
        list(opt = opt, anchor = sp$anchor_of(opt))
      for (kind in sp$kinds) {
        # A single-country group has no nested RE stats: the 04_fit_*_brms_cb
        # scripts skip compute_nested_surv_re_stats_brms() when n_countries < 2.
        if (
          grepl("re_stats$", kind) && .deploy_group_n_countries(dep, g) < 2L
        ) {
          say("[deploy] group=%s: single country -- no %s to collect", g, kind)
          next
        }
        src <- sp$make_fn(opt, kind)
        plan_rows[[length(plan_rows) + 1L]] <- data.frame(
          track = sp$track,
          group = grouping_label_token(.deploy_group_full(dep, g)),
          lag = lag,
          kind = kind,
          src = src,
          file = basename(src),
          is_fit = grepl("_fit_", basename(src)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  plan <- do.call(rbind, plan_rows)

  # One read per source fit, shared by the episode check (gid_season) and the
  # convergence gate (deploy.convergence.action != off); neither -> no read,
  # exactly as before the gate existed.
  conv_cfg <- deploy_convergence_cfg(dep)
  src_reads <- .deploy_read_source_fits(
    plan,
    episode = identical(ep_key, "gid_season"),
    diagnostics = !identical(conv_cfg$action, "off"),
    read_fn = read_fn,
    diag_fn = diag_fn
  )

  # --- episode definition: fits + the prep artefacts they trained on --------
  # BEFORE anything in dest is pruned or overwritten, so a mismatch leaves the
  # current deployment untouched. Every source fit must carry the deploy's
  # episode record and -- under gid_season -- the same episode digest as the
  # four prep artefacts currently on disk for that key (the size frame stays a
  # simulator input -- population / admin_level fill and 03_prepare_dur.R's
  # source -- although size is not fit), i.e. the fits agree with
  # each other AND with the panels the smoothing / forecast / hindcast runs will
  # read. Runs under --dry-run too (sources not produced yet are skipped there;
  # the copy below still stops on them). Silent for an lp_season rollout whose
  # fits and preps carry no gid_season stamp -- so the fit reads are skipped
  # there (the per-lag config check above already pins its key) and a default
  # deploy costs exactly what it did before.
  if (identical(ep_key, "gid_season")) {
    ep_recs <- src_reads$episode
    o1 <- fit_opts[[1L]]$opt
    ep_opt <- list(
      baseline = o1$baseline %||% "1990_2019",
      episodes = .episodes_block(ep_spec, always = "national_zero_fill")
    )
    preps <- list(
      surv = make_surv_dat_file(ep_opt),
      size_data = make_size_dat_file(ep_opt),
      dur_data = make_duration_dat_file(ep_opt),
      mon_data = make_monthly_dat_file(ep_opt)
    )
    if (file.exists(preps$surv)) {
      assert_occ_panel_variant(
        ep_opt,
        preps$surv,
        track = "[deploy] occurrence prep"
      )
    }
    for (n in names(preps)) {
      if (file.exists(preps[[n]])) ep_recs[[n]] <- episode_stamp_of(preps[[n]])
    }
    assert_episode_coherence(
      ep_recs,
      what = "deploy fits and prep artefacts",
      require_gid_season = identical(ep_key, "gid_season")
    )
    say(
      "[deploy] OK: %d fit(s) + %d prep artefact(s) share episodes.key = %s",
      sum(file.exists(plan$src[plan$is_fit])),
      sum(vapply(preps, file.exists, logical(1))),
      ep_key
    )
  }

  # --- convergence gate (deploy.convergence) --------------------------------
  # Also BEFORE the prune / copy: a refusal leaves dest untouched. Off (the
  # default) reads, prints and writes nothing.
  conv <- deploy_check_convergence(
    plan,
    dep,
    dry_run = dry_run,
    say = say,
    diags = src_reads$diagnostics,
    cfg = conv_cfg
  )

  # --- prune stale: any fit-prefixed rds in dest not in the expected set ----
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  prefixes <- vapply(specs, function(s) s$prefix, character(1))
  existing <- unlist(lapply(prefixes, function(p) {
    list.files(dest, pattern = sprintf("^%s.*\\.rds$", p), full.names = TRUE)
  }))
  stale <- existing[!(basename(existing) %in% plan$file)]
  # Never prune a fit of the OTHER episode definition: that is another
  # rollout's deployment, not a stale leftover of this one (the "_ep-gs" token
  # is utils/paths.R::.episode_key_path_token()). Give each episode definition its
  # own deploy.fits.dest.
  foreign <- .deploy_foreign_episode_artefacts(
    stale,
    ep_key,
    ep_nz,
    occ_prefix = specs$occ$prefix,
    national_zero_threshold = ep_nzt,
    occ_panel = ep_spec
  )
  if (length(foreign)) {
    stop(
      sprintf(
        paste0(
          "[deploy] %s holds %d artefact(s) of the other episode definition (this ",
          "deploy is episodes.key = '%s'), e.g. %s. Refusing to prune another ",
          "rollout's fits -- set a separate deploy.fits.dest."
        ),
        dest,
        length(foreign),
        ep_key,
        paste(utils::head(basename(foreign), 3), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  n_pruned <- 0L
  for (f in stale) {
    say(
      "[deploy] %s stale artefact: %s",
      if (dry_run) "WOULD prune" else "pruning",
      basename(f)
    )
    if (!dry_run) {
      file.remove(f)
    }
    n_pruned <- n_pruned + 1L
  }

  # --- copy ----------------------------------------------------------------
  n_copied <- 0L
  for (i in seq_len(nrow(plan))) {
    src <- plan$src[[i]]
    if (!file.exists(src)) {
      stop(
        "[deploy] source artefact missing (fit not produced yet?): ",
        src,
        call. = FALSE
      )
    }
    dst <- file.path(dest, plan$file[[i]])
    say("[deploy] %s %s", if (dry_run) "WOULD copy" else "copy", plan$file[[i]])
    if (!dry_run) {
      file.copy(src, dst, overwrite = TRUE)
    }
    n_copied <- n_copied + 1L
  }

  # --- stamp the `opt` attribute onto freshly-copied fits missing it --------
  # In-process equivalent of patch_opt_attr.R: occ/monthly fits already
  # carry opt; older duration fits may not. Gate on the anchor, using the opt
  # we already resolved for naming.
  n_patched <- 0L
  if (!dry_run) {
    for (i in which(plan$is_fit)) {
      dst <- file.path(dest, plan$file[[i]])
      fit <- readRDS(dst)
      if (!is.null(attr(fit, "opt"))) {
        next
      }
      key <- paste(plan$track[[i]], plan$group[[i]], sep = ":")
      fo <- fit_opts[[key]]
      if (is.null(fo) || !isTRUE(fo$anchor)) {
        next
      } # anchor-gated, like patch_opt_attr.R
      attr(fit, "opt") <- fo$opt
      saveRDS(fit, dst)
      say("[deploy] stamped opt attr onto %s", plan$file[[i]])
      n_patched <- n_patched + 1L
    }
  }

  # --- MANIFEST.csv --------------------------------------------------------
  manifest <- data.frame(
    track = plan$track,
    group = plan$group,
    lag = plan$lag,
    kind = plan$kind,
    file = plan$file,
    stringsAsFactors = FALSE
  )
  if (!dry_run) {
    utils::write.csv(
      manifest,
      file.path(dest, "MANIFEST.csv"),
      row.names = FALSE
    )
    # The rule the lags in MANIFEST.csv were selected by (a sidecar, so the
    # MANIFEST readers -- results_core.R::load_best_lag_manifest() -- are unchanged).
    yaml::write_yaml(
      list(
        rule = ls_cfg$rule,
        se_mult = ls_cfg$se_mult,
        min_lag = if (is.finite(ls_cfg$min_lag)) ls_cfg$min_lag else "-Inf",
        occ_floor_offset = ls_cfg$occ_floor_offset
      ),
      file.path(dest, "LAG_SELECTION.yml")
    )
    # The convergence verdicts the deployed fits passed (or were allowed)
    # under; under action off a CONVERGENCE.csv of an earlier deploy would
    # misdescribe this one, so it is removed.
    conv_csv <- file.path(dest, "CONVERGENCE.csv")
    if (!is.null(conv)) {
      utils::write.csv(conv$table, conv_csv, row.names = FALSE)
    } else if (file.exists(conv_csv)) {
      file.remove(conv_csv)
      say(
        "[deploy] removed stale %s (deploy.convergence.action = off)",
        conv_csv
      )
    }
  }

  # --- post-deploy assertion: exactly one fit per (track, group) -----------
  if (!dry_run) {
    for (g in dep$fits$groups) {
      for (sp in specs) {
        pat <- sprintf(
          "^%s.*_fit_.*_grp-%s-%s_.*\\.rds$",
          sp$prefix,
          dep$grouping %||% "enso_grouping",
          .deploy_group_full(dep, g)
        )
        hits <- list.files(dest, pattern = pat)
        if (length(hits) != 1L) {
          stop(
            sprintf(
              "[deploy] %d %s fit(s) match group '%s' in %s -- expected 1:\n  %s",
              length(hits),
              sp$prefix,
              g,
              dest,
              paste(hits, collapse = "\n  ")
            ),
            call. = FALSE
          )
        }
      }
    }
    say(
      "[deploy] OK: every (track, group) resolves to exactly one fit in %s",
      dest
    )
  }

  list(
    selected = selected,
    manifest = manifest,
    dest = dest,
    copied = n_copied,
    pruned = n_pruned,
    patched = n_patched,
    dry_run = dry_run,
    convergence = conv
  )
}

# ---------------------------------------------------------------------------
# Deploy <-> results coupling
# ---------------------------------------------------------------------------

# The ONLY auto-derived results.inputs leaves (deploy leaf -> results leaf).
# scope.* is deliberately NOT derived -- it encodes reporting choices that may
# legitimately differ (compare-only).
# Root of every spatPomp run tree, as the deploy YAML and results.inputs spell
# it (no trailing slash -- 08c writes under file.path(cfg$outdir %||%
# "generated_data", "sim_spatpomp", run_tag), so this literal matches the cache
# path byte for byte; paths.R::make_smooth_run_file() with a resolved opt would
# give "generated_data//sim_spatpomp").
.DEPLOY_SIM_ROOT <- "generated_data/sim_spatpomp"

#' The historical smoothing run's `parts` or `cache` dir for this deploy.
#' `cache` is what deploy.smoothing_input must equal.
deploy_smoothing_dir <- function(dep, sub = c("parts", "cache")) {
  file.path(.DEPLOY_SIM_ROOT, dep$run_tags$historical, match.arg(sub))
}

.deploy_derived_results_inputs <- function(dep, config_path) {
  stem <- .deploy_batch_stem
  list(
    best_lag_dir = dep$fits$dest,
    manifest = dep$outputs$forecast_manifest,
    # The attribution factual's manifest (--stage gen_attrib): 09's attribution
    # steps read the factual arm from it; "" = the forecast manifest's
    # (stitched) factual.
    attrib_manifest = as.character(
      dep$outputs$forecast_attrib_manifest %||% ""
    ),
    hindcast_manifest = dep$outputs$hindcast_manifest,
    episode_key = dep$episodes$key %||% "lp_season",
    episode_national_zero_fill = isTRUE(dep$episodes$national_zero_fill),
    episode_national_zero_threshold = .episodes_nz_threshold(dep$episodes),
    # read in-function (default "" = legacy) when the deploy block omits it.
    episode_national_zero_rule = .occ_panel_spec_from_episodes(
      dep$episodes
    )$national_zero_rule,
    # read in-function (defaults season / drop) likewise.
    episode_at_risk_start = .occ_panel_spec_from_episodes(
      dep$episodes
    )$at_risk_start,
    episode_adm1_overlap = .occ_panel_spec_from_episodes(
      dep$episodes
    )$adm1_overlap,
    grouping = file.path("analysis/groupings", paste0(dep$grouping, ".csv")),
    sim_root = .DEPLOY_SIM_ROOT,
    smoothing_parts_dir = deploy_smoothing_dir(dep, "parts"),
    lag_loo_batch = list(
      default = stem(dep$fits$occ_dir),
      occurrence = stem(dep$fits$occ_dir),
      duration = stem(dep$fits$dur_dir),
      monthly = stem(dep$fits$mon_dir)
    )
  )
}

#' Derive the overlapping results.inputs.* leaves from the deploy config the
#' results block points at (opt$results$deploy_config), with EXPLICIT-over-
#' inherited precedence: any inputs leaf physically present in the results YAML
#' wins; everything else is taken from the deploy block.
#'
#' GUARD: an explicit leaf whose value DISAGREES with what the deploy derives
#' is refused unless its name is listed under results.allow_overrides, so a
#' results block cannot silently mix runs. A mixed-run results block is
#' allowed, but only as a stated decision.
#'
#' Called from make_options_from_config() post-processing when
#' results.deploy_config is non-null (config.R guards on the function existing).
#'
#' @param opt Resolved options list (has opt$results and opt$results$deploy_config).
#' @param results_config_path Path to the results YAML (-C), for explicit-key
#'   detection; defaults to the active -C/--config flag.
#' @return opt with opt$results$inputs updated.
resolve_results_inputs_from_deploy <- function(
  opt,
  results_config_path = NULL
) {
  dep <- load_deploy_config(opt$results$deploy_config)
  derived <- .deploy_derived_results_inputs(dep, opt$results$deploy_config)

  results_config_path <- results_config_path %||%
    .flag_val("-C") %||%
    .flag_val("--config")
  explicit <- character(0)
  if (!is.null(results_config_path) && file.exists(results_config_path)) {
    raw <- yaml::read_yaml(results_config_path)
    explicit <- names(raw$results$inputs %||% list())
  }

  allowed <- as.character(unlist(opt$results$allow_overrides %||% character(0)))
  unknown <- setdiff(allowed, names(derived))
  if (length(unknown)) {
    stop(
      "[results<-deploy] results.allow_overrides names leaves that are not ",
      "deploy-derived inputs: ",
      paste(unknown, collapse = ", "),
      ". Valid names: ",
      paste(names(derived), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  conflicts <- character(0)
  for (leaf in names(derived)) {
    if (leaf %in% explicit) {
      # explicit-in-YAML wins ...
      have <- opt$results$inputs[[leaf]]
      if (
        !.deploy_leaf_same(have, derived[[leaf]]) &&
          !(leaf %in% allowed)
      ) {
        # ... but only when it agrees, or is allowed
        conflicts <- c(
          conflicts,
          sprintf(
            "  inputs.%s\n    results: %s%s\n    deploy : %s",
            leaf,
            .deploy_leaf_str(have),
            .deploy_leaf_type_hint(have, derived[[leaf]]),
            .deploy_leaf_str(derived[[leaf]])
          )
        )
      }
      next
    }
    opt$results$inputs[[leaf]] <- derived[[leaf]]
  }
  if (length(conflicts)) {
    stop(
      sprintf(
        paste0(
          "[results<-deploy] %d explicit results.inputs leaf/leaves disagree with the ",
          "deploy config they point at (%s):\n%s\n",
          "Either align the results YAML to the deploy\n",
          "  Rscript analysis/utils/deploy_prod_run.R -C %s --stage align_results\n",
          "or, if the mixed-run choice is deliberate, name the leaf under ",
          "results.allow_overrides."
        ),
        length(conflicts),
        opt$results$deploy_config,
        paste(conflicts, collapse = "\n"),
        opt$results$deploy_config
      ),
      call. = FALSE
    )
  }
  opt
}

# Fully-resolved results.inputs for a results YAML, absent any deploy coupling:
# explicit-over-default (recursive merge over get_default_options()$results$inputs).
.resolved_results_inputs <- function(results_path) {
  raw <- yaml::read_yaml(results_path)
  defaults <- get_default_options()$results$inputs
  .merge_recursive(raw$results$inputs, defaults)
}

# Does an explicit results.inputs leaf agree with its deploy-derived value?
# Compared as .deploy_leaf_str() strings, except that a deploy-derived LOGICAL
# must also be logical in the results YAML: a quoted "TRUE" prints the same as
# TRUE but reads as FALSE through isTRUE() (episode_national_zero_fill was
# accepted here and then silently read as off). A deploy-derived NUMBER must
# be a number in the results YAML and equal as a number (a quoted "60" is
# refused); int == double still agrees.
.deploy_leaf_same <- function(have, want) {
  if (is.numeric(want)) {
    return(is.numeric(have) && identical(as.numeric(have), as.numeric(want)))
  }
  identical(.deploy_leaf_str(have), .deploy_leaf_str(want)) &&
    (!is.logical(want) || is.logical(have))
}

# The refusal's fix for a leaf that prints like the deploy's but has the wrong
# YAML type ("" when the type is not the problem).
.deploy_leaf_type_hint <- function(have, want) {
  if (is.logical(want) && !is.logical(have)) {
    return(" (a string; write the YAML boolean true/false, unquoted)")
  }
  if (is.numeric(want) && !is.numeric(have)) {
    return(" (a string; write the YAML number, unquoted)")
  }
  ""
}

# Flatten a leaf value to a single comparable string. A named list (e.g.
# lag_loo_batch) prints as k=v pairs, an unnamed list (YAML sequence) as its
# elements, an empty mapping as "{}", NULL as "<absent>"; atomics as-is.
.deploy_leaf_str <- function(v) {
  if (is.null(v)) {
    return("<absent>")
  }
  if (is.list(v)) {
    if (!length(v)) {
      return("{}")
    }
    if (is.null(names(v))) {
      return(paste(vapply(v, .deploy_leaf_str, ""), collapse = ","))
    }
    return(paste(sprintf("%s=%s", names(v), unlist(v)), collapse = ","))
  }
  paste(as.character(v), collapse = ",")
}

#' Leaf-by-leaf comparison of two (nested) lists.
#'
#' Values are compared as `.deploy_leaf_str()` strings, so the comparison is
#' order-insensitive and blind to integer-vs-double (both serialise the same in
#' YAML) but not to TRUE-vs-"yes".
#'
#' @param a,b     Lists to compare.
#' @param ignore  Leaf paths to drop (a path also drops everything below it).
#' @param flatten TRUE: compare every dotted leaf (config_flatten_leaves()).
#'                FALSE: compare the top-level entries as given.
#' @return data.frame(leaf, a, b, status) with status one of match | differs |
#'   only_a | only_b, rows in first-appearance order (a's leaves, then b's).
deploy_leaf_diff <- function(a, b, ignore = character(0), flatten = TRUE) {
  fa <- if (flatten) config_flatten_leaves(a %||% list()) else (a %||% list())
  fb <- if (flatten) config_flatten_leaves(b %||% list()) else (b %||% list())
  leaves <- union(names(fa), names(fb))
  dropped <- vapply(
    leaves,
    function(p) {
      any(p == ignore | startsWith(p, paste0(ignore, ".")))
    },
    logical(1)
  )
  leaves <- leaves[!dropped]
  if (!length(leaves)) {
    return(data.frame(
      leaf = character(0),
      a = character(0),
      b = character(0),
      status = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(leaves, function(leaf) {
    in_a <- leaf %in% names(fa)
    in_b <- leaf %in% names(fb)
    sa <- if (in_a) .deploy_leaf_str(fa[[leaf]]) else "<absent>"
    sb <- if (in_b) .deploy_leaf_str(fb[[leaf]]) else "<absent>"
    st <- if (!in_a) {
      "only_b"
    } else if (!in_b) {
      "only_a"
    } else if (identical(sa, sb)) {
      "match"
    } else {
      "differs"
    }
    data.frame(
      leaf = leaf,
      a = sa,
      b = sb,
      status = st,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Print the non-matching rows of a deploy_leaf_diff() table.
deploy_print_leaf_diff <- function(
  tab,
  a_label = "a",
  b_label = "b",
  indent = "    "
) {
  tab <- tab[tab$status != "match", , drop = FALSE]
  for (i in seq_len(nrow(tab))) {
    cat(sprintf(
      "%s%-7s %s\n%s        %s: %s\n%s        %s: %s\n",
      indent,
      tab$status[[i]],
      tab$leaf[[i]],
      indent,
      a_label,
      tab$a[[i]],
      indent,
      b_label,
      tab$b[[i]]
    ))
  }
  invisible(tab)
}

#' Print a leaf-diff table over the auto-derived keys: deploy-derived value,
#' results effective value, and match/mismatch/override status. Runnable both
#' standalone and from the orchestrator after a deploy stage.
#'
#' @param deploy_path  Path to the deploy YAML.
#' @param results_path Path to the results YAML.
#' @return (invisibly) a data.frame of the comparison.
compare_deploy_results <- function(deploy_path, results_path) {
  dep <- load_deploy_config(deploy_path)
  derived <- .deploy_derived_results_inputs(dep, deploy_path)
  res_in <- .resolved_results_inputs(results_path)
  raw <- yaml::read_yaml(results_path)
  explicit <- names(raw$results$inputs %||% list())

  allowed <- as.character(unlist(raw$results$allow_overrides %||% character(0)))

  d <- deploy_leaf_diff(
    derived,
    res_in[intersect(names(derived), names(res_in))],
    flatten = FALSE
  )
  d <- d[d$leaf %in% names(derived), , drop = FALSE]
  # Status on the typed compare the resolve guard uses (.deploy_leaf_same()),
  # so a quoted "60" or "TRUE" is a MISMATCH here too.
  same <- d$status != "only_a" &
    vapply(
      d$leaf,
      function(leaf) .deploy_leaf_same(res_in[[leaf]], derived[[leaf]]),
      logical(1),
      USE.NAMES = FALSE
    )
  # An explicit leaf that disagrees is an "override" only when
  # results.allow_overrides names it; otherwise 09/10 refuse it
  # (resolve_results_inputs_from_deploy()).
  tab <- data.frame(
    leaf = d$leaf,
    deploy = d$a,
    results = ifelse(d$status == "only_a", "", d$b),
    status = ifelse(
      same,
      "match",
      ifelse(
        d$leaf %in% explicit & d$leaf %in% allowed,
        "override (allowed)",
        "MISMATCH"
      )
    ),
    stringsAsFactors = FALSE
  )

  cat("\n[deploy<->results] leaf comparison\n")
  cat(sprintf("  deploy : %s\n  results: %s\n\n", deploy_path, results_path))
  for (i in seq_len(nrow(tab))) {
    cat(sprintf(
      "  %-20s %-19s deploy=%s | results=%s\n",
      tab$leaf[[i]],
      tab$status[[i]],
      tab$deploy[[i]],
      tab$results[[i]]
    ))
  }
  if (any(tab$status == "MISMATCH")) {
    cat(
      "\n  [!] MISMATCH leaves are refused by 09/10 when results.deploy_config points at ",
      "this deploy.\n      Align: Rscript analysis/utils/deploy_prod_run.R -C ",
      deploy_path,
      " --stage align_results\n      or name the leaf under results.allow_overrides.\n",
      sep = ""
    )
  }
  invisible(tab)
}

# --- Results YAML alignment: a TEXT edit, never a yaml round trip ------------
# manuscript.yml is the manuscript's decision record and is mostly comments.
# A yaml::read_yaml() -> yaml::write_yaml() round trip would drop every one of
# them and re-quote every string. The aligner edits the file line by line: one line
# per leaf it owns, every other line byte-identical, then re-reads the result
# and asserts each leaf parsed to the intended value before it says "aligned".

# Line number of a dotted leaf under the top-level `results:` mapping, found by
# an indentation walk (blank and comment lines skipped; only lines at the
# block's own child indentation can match, so a same-named grandchild is never
# picked). NA_integer_ when absent.
.yaml_leaf_line <- function(lines, dotted) {
  indent_of <- function(ln) nchar(sub("^( *).*$", "\\1", ln))
  skip <- function(ln) !nzchar(trimws(ln)) || grepl("^\\s*#", ln)
  pos <- 0L
  parent_indent <- -1L
  for (k in c("results", strsplit(dotted, ".", fixed = TRUE)[[1]])) {
    child_indent <- NA_integer_
    hit <- NA_integer_
    i <- pos + 1L
    while (i <= length(lines)) {
      ln <- lines[i]
      if (!skip(ln)) {
        ind <- indent_of(ln)
        if (ind <= parent_indent) {
          break
        } # left the parent block
        if (is.na(child_indent)) {
          child_indent <- ind
        }
        if (ind == child_indent && grepl(sprintf("^ *%s:( |$)", k), ln)) {
          hit <- i
          break
        }
      }
      i <- i + 1L
    }
    if (is.na(hit)) {
      return(NA_integer_)
    }
    pos <- hit
    parent_indent <- indent_of(lines[hit])
  }
  pos
}

# `key:<pad>old_value[<pad># comment]` -> same key, same padding, same trailing
# comment, value replaced by a double-quoted scalar (the file's house style);
# logicals and numbers are written bare (typed).
.yaml_set_scalar_line <- function(line, value) {
  m <- regmatches(
    line,
    regexec("^( *[A-Za-z0-9_]+:)( *)([^#]*?)( *#.*)?$", line, perl = TRUE)
  )[[1]]
  if (!length(m) || !nzchar(trimws(m[4]))) {
    stop(
      "[results<-deploy] not a scalar `key: value` line: ",
      sQuote(line),
      call. = FALSE
    )
  }
  # Logicals are written as YAML booleans (true / false), never quoted: a
  # quoted "TRUE" reads back as a string, which isTRUE() treats as FALSE (the
  # national zero-fill flag was silently read as off that way).
  # Numbers likewise bare: a quoted "60" is a string to every reader.
  v <- if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    tolower(as.character(value))
  } else if (is.numeric(value) && length(value) == 1L && !is.na(value)) {
    format(value, scientific = FALSE, digits = 15L)
  } else {
    sprintf('"%s"', as.character(value))
  }
  paste0(m[2], m[3], v, m[5])
}

# Comment lines that still name a run other than the deploy's (e.g. "prod_run<M>"
# or "_run<M>/" while aligning to prod_run<N>). Prose cannot be rewritten
# safely, so these are REPORTED for a human edit, never changed. A run id may
# carry one letter suffix ("run<N>b", a re-forecast of run<N>), read only when no
# further letter follows, so "run<N>b" is its own run but "run<N>occ" names
# run<N>.
.RUN_ID_RE <- "(?<![A-Za-z])run([0-9]+(?:[a-z](?![A-Za-z]))?)"
.results_stale_run_comments <- function(lines, deploy_path) {
  tag <- derive_config_run_tag(deploy_path)
  n <- sub("^.*run", "", tag)
  is_comment <- grepl("#", lines, fixed = TRUE)
  hits <- regmatches(lines, gregexpr(.RUN_ID_RE, lines, perl = TRUE))
  stale <- vapply(
    seq_along(lines),
    function(i) {
      if (!is_comment[i] || !length(hits[[i]])) {
        return(FALSE)
      }
      # only the part after the first `#` is comment text
      cm <- regmatches(
        sub("^[^#]*#", "", lines[i]),
        gregexpr(.RUN_ID_RE, sub("^[^#]*#", "", lines[i]), perl = TRUE)
      )[[1]]
      any(sub("^run", "", cm) != n)
    },
    logical(1)
  )
  data.frame(
    line = which(stale),
    text = trimws(lines[stale]),
    stringsAsFactors = FALSE
  )
}

# Where the pre-edit copy goes: under the deploy's provenance dir in
# generated_data/ (gitignored), never beside the config.
.deploy_results_backup_path <- function(results_path, deploy_path) {
  run_dir <- make_deploy_run_file(
    list(outdir = "generated_data/"),
    "deploy",
    config_path = deploy_path
  )
  file.path(
    run_dir,
    sprintf(
      "%s.%s.bak",
      basename(results_path),
      format(Sys.time(), "%Y%m%d-%H%M%S")
    )
  )
}

#' Plan the line-level edits that align a results YAML to a deploy YAML:
#' `results.deploy_config` <- deploy_path, and every auto-derived inputs.*
#' leaf (.deploy_derived_results_inputs()) <- its deploy-derived value.
#' Pure: reads both files, writes nothing.
#'
#' @return list(lines, new_lines, targets = named list(dotted leaf -> value),
#'   changes = data.frame(line, leaf, old, new) for the lines that differ,
#'   stale_comments = data.frame(line, text) of comments naming another run).
deploy_plan_results_alignment <- function(results_path, deploy_path) {
  dep <- load_deploy_config(deploy_path)
  derived <- .deploy_derived_results_inputs(dep, deploy_path)
  targets <- c(
    list(deploy_config = deploy_path),
    config_flatten_leaves(list(inputs = derived))
  )
  # The sweep comparison read by report_kl_sweep: a top-level, unstamped
  # results key, aligned only when this deploy derives its grouping from one.
  if (nzchar(.deploy_kl_sweep_batch(dep))) {
    targets$kl_sweep_dir <- deploy_kl_sweep_compare_dir(dep)
    # ... and its deployed cell (report_membership).
    cell <- as.character(dep$kl_sweep$cell %||% "")
    if (nzchar(cell)) targets$kl_sweep_cell <- cell
  }
  # The past-ENSO-event manifest read by report_event_attrib: top level,
  # unstamped, aligned only when this deploy sets events.
  ea <- dep$event_attribution %||% list()
  if (
    nzchar(as.character(ea$events %||% "")) &&
      nzchar(as.character(ea$manifest %||% ""))
  ) {
    targets$event_manifest <- as.character(ea$manifest)
  }

  lines <- readLines(results_path, warn = FALSE)
  new <- lines
  rows <- list()
  for (leaf in names(targets)) {
    ln <- .yaml_leaf_line(lines, leaf)
    if (is.na(ln)) {
      stop(
        "[results<-deploy] ",
        results_path,
        " has no `results.",
        leaf,
        "` line to rewrite -- add the key (the aligner never inserts one).",
        call. = FALSE
      )
    }
    cand <- .yaml_set_scalar_line(lines[ln], targets[[leaf]])
    if (!identical(cand, lines[ln])) {
      new[ln] <- cand
      rows[[length(rows) + 1L]] <- data.frame(
        line = ln,
        leaf = leaf,
        old = lines[ln],
        new = cand,
        stringsAsFactors = FALSE
      )
    }
  }
  changes <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      line = integer(0),
      leaf = character(0),
      old = character(0),
      new = character(0),
      stringsAsFactors = FALSE
    )
  }
  list(
    lines = lines,
    new_lines = new,
    targets = targets,
    changes = changes,
    stale_comments = .results_stale_run_comments(new, deploy_path)
  )
}

.print_results_alignment_plan <- function(plan, results_path, deploy_path) {
  cat(sprintf(
    "\n[deploy<->results] align %s -> %s\n",
    results_path,
    deploy_path
  ))
  if (!nrow(plan$changes)) {
    cat("  every owned leaf already matches; nothing to change.\n")
  } else {
    cat(sprintf("  %d line(s) to change:\n", nrow(plan$changes)))
    for (i in seq_len(nrow(plan$changes))) {
      cat(sprintf(
        "  L%-4d results.%s\n      - %s\n      + %s\n",
        plan$changes$line[i],
        plan$changes$leaf[i],
        trimws(plan$changes$old[i]),
        trimws(plan$changes$new[i])
      ))
    }
  }
  if (nrow(plan$stale_comments)) {
    cat(sprintf(
      "  [!] %d comment line(s) still name another run -- edit by hand:\n",
      nrow(plan$stale_comments)
    ))
    for (i in seq_len(nrow(plan$stale_comments))) {
      cat(sprintf(
        "      L%-4d %s\n",
        plan$stale_comments$line[i],
        plan$stale_comments$text[i]
      ))
    }
  }
  invisible(plan)
}

#' Align a results YAML to a deploy YAML (deploy_plan_results_alignment()) as a
#' TEXT edit. Default is a dry-run report; write = TRUE persists, after copying
#' the file to generated_data/deploy/<deploy stem>/, then re-reads the result
#' and REFUSES (restoring the copy) unless every owned leaf parsed to its
#' intended value and compare_deploy_results() is all-match.
#'
#' @param results_path Path to the results YAML (rewritten when write = TRUE).
#' @param deploy_path  Path to the deploy YAML.
#' @param write        If TRUE, overwrite results_path; else report only.
#' @return (invisibly) the plan (see deploy_plan_results_alignment()).
align_results_to_deploy <- function(results_path, deploy_path, write = FALSE) {
  plan <- deploy_plan_results_alignment(results_path, deploy_path)
  .print_results_alignment_plan(plan, results_path, deploy_path)

  if (!isTRUE(write)) {
    cat("\n[deploy<->results] dry-run (write = FALSE): no file changed.\n")
    return(invisible(plan))
  }
  if (nrow(plan$changes)) {
    bak <- .deploy_results_backup_path(results_path, deploy_path)
    dir.create(dirname(bak), recursive = TRUE, showWarnings = FALSE)
    file.copy(results_path, bak, overwrite = TRUE)
    writeLines(plan$new_lines, results_path)
    # Verify by re-reading: every owned leaf must parse to the intended value.
    raw <- yaml::read_yaml(results_path)
    bad <- Filter(
      function(leaf) {
        !.deploy_leaf_same(
          config_get_leaf(raw$results, leaf),
          plan$targets[[leaf]]
        )
      },
      names(plan$targets)
    )
    if (length(bad)) {
      file.copy(bak, results_path, overwrite = TRUE)
      stop(
        "[results<-deploy] re-read of ",
        results_path,
        " does not carry the intended ",
        "value for results.",
        paste(bad, collapse = ", results."),
        " -- edit refused and the file restored from ",
        bak,
        call. = FALSE
      )
    }
    cat(sprintf(
      "\n[deploy<->results] aligned %s to %s (pre-edit copy: %s)\n",
      results_path,
      deploy_path,
      bak
    ))
  }
  tab <- compare_deploy_results(deploy_path, results_path)
  if (any(tab$status != "match")) {
    stop(
      "[results<-deploy] post-align comparison is not all-match: ",
      paste(tab$leaf[tab$status != "match"], collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(plan$stale_comments)) {
    cat(sprintf(
      "\n  [!] %d comment line(s) above still name another run; update them by hand.\n",
      nrow(plan$stale_comments)
    ))
  }
  invisible(plan)
}

# ---------------------------------------------------------------------------
# Sim-template generation (--stage gen_templates)
# ---------------------------------------------------------------------------
# The three spatPomp templates of a rollout (deploy.templates.*) are GENERATED,
# not hand-written: each arm is
#     analysis/configs/sim_templates/<arm>.yml        (run-agnostic skeleton)
#   + deploy sim.shared, then sim.<arm>               (dotted-leaf overlay; later wins)
#   + .deploy_sim_derived_leaves()                    (pure functions of the deploy block)
#   + apply_episode_key_to_sim_config(deploy.episodes.key)
# Derived leaves cannot be overlaid (refused), required leaves must be
# overlaid (refused when absent), and simulator.run_id is provenance only: it
# never makes an existing template "differ".

.DEPLOY_SIM_ARMS <- c("historical", "forecast", "counterfactual")
DEPLOY_SIM_SKELETON_DIR <- "analysis/configs/sim_templates"
.DEPLOY_SIM_GENERATOR <- "analysis/utils/deploy_prod_run.R --stage gen_templates"

# Leaves every run must choose in its `sim:` overlay (the skeletons omit them).
.DEPLOY_SIM_REQUIRED_LEAVES <- list(
  historical = c(
    "simulator.filter.measurement",
    "simulator.filter.partition_by"
  ),
  forecast = "simulator.filter.measurement",
  counterfactual = "simulator.filter.measurement"
)

# Leaves never compared when deciding unchanged/differs.
.DEPLOY_SIM_IGNORED_LEAVES <- "simulator.run_id"

deploy_sim_skeleton_path <- function(arm) {
  file.path(
    DEPLOY_SIM_SKELETON_DIR,
    paste0(match.arg(arm, .DEPLOY_SIM_ARMS), ".yml")
  )
}

#' Where an arm's template is written: analysis/configs/<tag>/<tag>.yml with
#' tag = run_tags.<arm>. deploy.templates.<arm> (what bash and the other stages
#' read) must name exactly that file.
deploy_sim_template_path <- function(dep, arm) {
  arm <- match.arg(arm, .DEPLOY_SIM_ARMS)
  tag <- dep$run_tags[[arm]]
  path <- assert_basename_length(
    file.path("analysis/configs", tag, paste0(tag, ".yml")),
    "deploy_sim_template_path"
  )
  if (!identical(path, dep$templates[[arm]])) {
    stop(
      "[deploy] deploy.templates.",
      arm,
      " is '",
      dep$templates[[arm]],
      "' but gen_templates writes the ",
      arm,
      " template to '",
      path,
      "' (analysis/configs/<run_tags.",
      arm,
      ">/<run_tags.",
      arm,
      ">.yml). ",
      "Make the two agree in the deploy YAML.",
      call. = FALSE
    )
  }
  path
}

#' simulator.run_id for a generated template: <yyyymmdd>_<arm>_<fit lineage>,
#' the lineage grouping the three tracks by the run that fit them in occ, dur,
#' mon order (e.g. an occurrence refit on an earlier run's duration and monthly
#' fits: prodrunBocc_prodrunAdurmon).
#' Provenance only -- spatPomp output paths depend on run_tag alone.
deploy_sim_run_id <- function(dep, arm, date = Sys.Date()) {
  stems <- gsub(
    "_",
    "",
    vapply(
      c(occ = dep$fits$occ_dir, dur = dep$fits$dur_dir, mon = dep$fits$mon_dir),
      .deploy_batch_stem,
      ""
    )
  )
  lineage <- vapply(
    unique(stems),
    function(st) {
      paste0(st, paste(names(stems)[stems == st], collapse = ""))
    },
    ""
  )
  paste(
    format(as.Date(date), "%Y%m%d"),
    arm,
    paste(lineage, collapse = "_"),
    sep = "_"
  )
}

#' The leaves of an arm's template that are pure functions of the deploy block.
#' Named list: dotted leaf path -> value. The historical arm derives nothing
#' beyond the shared leaves; the forward arms add their smoothing hand-off.
.deploy_sim_derived_leaves <- function(dep, arm, date = Sys.Date()) {
  d <- list(
    simulator.run_tag = dep$run_tags[[arm]],
    simulator.run_id = deploy_sim_run_id(dep, arm, date),
    simulator.grouping.source = dep$grouping,
    simulator.grouping.fits_dir = dep$fits$dest,
    simulator.nbr_contagion.enable = isTRUE(dep$nbr_contagion$enable),
    simulator.nbr_contagion.rho_km = dep$nbr_contagion$rho_km
  )
  # Registry partition (partition_by: registry): the partition registry and the
  # neighbour-contagion registry (default: the same, so contagion blocks refine
  # the TUs) are deploy decisions, derived into every arm. Absent otherwise, so
  # existing country / enso_region templates are unchanged leaf for leaf.
  if (identical(dep$partition_by %||% "country", "registry")) {
    ps <- deploy_partition_source(dep)
    d$simulator.filter.partition_source <- ps
    d$simulator.nbr_contagion.partition_source <- deploy_nbr_partition_source(
      dep
    )
  }
  if (!identical(arm, "historical")) {
    cache <- deploy_smoothing_dir(dep, "cache")
    if (!identical(cache, dep$smoothing_input)) {
      stop(
        "[deploy] deploy.smoothing_input is '",
        dep$smoothing_input,
        "' but the ",
        "historical smoothing run (run_tags.historical) writes its cache to '",
        cache,
        "'. The forecast/counterfactual arms would hand off from the wrong cache.",
        call. = FALSE
      )
    }
    d$simulator.filter.smoothing_input <- cache
  }
  d
}

#' deploy.partition_source, checked: required (and must name an existing,
#' path-safe registry) when deploy.partition_by is "registry".
deploy_partition_source <- function(dep) {
  ps <- as.character(dep$partition_source %||% "")
  if (!identical(dep$partition_by %||% "country", "registry")) {
    return(ps)
  }
  if (!nzchar(ps)) {
    stop(
      "[deploy] deploy.partition_by is 'registry' but deploy.partition_source is empty: ",
      "name the partition registry (analysis/groupings/<name>.csv).",
      call. = FALSE
    )
  }
  reg <- load_grouping(ps) # errors on a missing registry
  assert_partition_labels(reg$group, ps)
  ps
}

#' The neighbour-contagion partition registry of a registry-partitioned deploy:
#' deploy.contagion.partition_source, else the partition registry.
deploy_nbr_partition_source <- function(dep) {
  x <- as.character(dep$contagion$partition_source %||% "")
  if (!nzchar(x)) {
    return(deploy_partition_source(dep))
  }
  # A distinct contagion registry must REFINE the partition registry: every
  # contagion block's countries sit in one partition unit (the sim-time kernel
  # is cut to the unit; spatpomp_core.R::check_partition_vs_grouping()).
  if (identical(dep$partition_by %||% "country", "registry")) {
    preg <- load_grouping(deploy_partition_source(dep))
    creg <- load_grouping(x)
    unit <- stats::setNames(preg$group, preg$iso3)
    both <- creg[creg$iso3 %in% preg$iso3, , drop = FALSE]
    span <- vapply(
      split(unname(unit[both$iso3]), both$group),
      function(u) length(unique(u)),
      integer(1)
    )
    bad <- names(span)[span > 1L]
    if (length(bad)) {
      stop(
        sprintf(
          paste0(
            "[deploy] contagion registry '%s' does not refine partition registry ",
            "'%s': block(s) %s span several units."
          ),
          x,
          deploy_partition_source(dep),
          paste(bad, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  x
}

#' Build one arm's template (a config list; nothing is written).
deploy_build_sim_template <- function(
  dep,
  arm,
  date = Sys.Date(),
  skeleton = deploy_sim_skeleton_path(arm)
) {
  arm <- match.arg(arm, .DEPLOY_SIM_ARMS)
  cfg <- yaml::read_yaml(skeleton)
  derived <- .deploy_sim_derived_leaves(dep, arm, date)
  protected <- c(names(derived), "episodes")
  hits <- function(keys) {
    keys[vapply(
      keys,
      function(k) {
        any(
          k == protected |
            startsWith(k, paste0(protected, ".")) |
            startsWith(protected, paste0(k, "."))
        )
      },
      logical(1)
    )]
  }

  in_skel <- hits(config_leaf_paths(cfg))
  if (length(in_skel)) {
    stop(
      "[deploy] skeleton ",
      skeleton,
      " sets derived leaves: ",
      paste(in_skel, collapse = ", "),
      call. = FALSE
    )
  }

  overlay <- c(dep$sim$shared %||% list(), dep$sim[[arm]] %||% list())
  bad <- hits(names(overlay))
  if (length(bad)) {
    stop(
      "[deploy] deploy sim.shared / sim.",
      arm,
      " sets leaves derived from the ",
      "deploy block: ",
      paste(unique(bad), collapse = ", "),
      ". Set them through ",
      "their deploy keys (run_tags.*, fits.dest, grouping, nbr_contagion.*, ",
      "episodes.key, smoothing_input) instead.",
      call. = FALSE
    )
  }
  for (i in seq_along(overlay)) {
    # in order: sim.<arm> wins over sim.shared
    cfg <- config_set_leaf(cfg, names(overlay)[[i]], overlay[[i]])
  }
  for (k in names(derived)) {
    cfg <- config_set_leaf(cfg, k, derived[[k]])
  }
  cfg <- apply_episode_key_to_sim_config(
    cfg,
    dep$episodes$key,
    label = paste(arm, "template"),
    national_zero_fill = isTRUE(dep$episodes$national_zero_fill),
    national_zero_threshold = .episodes_nz_threshold(dep$episodes),
    occ_panel = .occ_panel_spec_from_episodes(dep$episodes)
  )

  req <- .DEPLOY_SIM_REQUIRED_LEAVES[[arm]]
  missing <- req[vapply(
    req,
    function(k) is.null(config_get_leaf(cfg, k)),
    logical(1)
  )]
  if (length(missing)) {
    stop(
      "[deploy] the ",
      arm,
      " template needs ",
      paste(missing, collapse = ", "),
      " -- set it in the deploy YAML's sim.shared or sim.",
      arm,
      " block.",
      call. = FALSE
    )
  }

  front <- intersect(c("backend", "run_tag", "run_id"), names(cfg$simulator))
  cfg$simulator <- cfg$simulator[c(front, setdiff(names(cfg$simulator), front))]
  cfg
}

#' TRUE when two templates agree on every leaf except simulator.run_id.
deploy_sim_templates_equal <- function(a, b) {
  all(
    deploy_leaf_diff(a, b, ignore = .DEPLOY_SIM_IGNORED_LEAVES)$status ==
      "match"
  )
}

#' Build all three arms and classify each against what is on disk.
#'
#' @return named list (one per arm) of list(arm, path, status, cfg, diff):
#'   status "absent" (no file), "unchanged" (equal but for run_id) or "differs";
#'   diff is the deploy_leaf_diff(on-disk, generated) table (NULL when absent).
deploy_plan_sim_templates <- function(dep, date = Sys.Date()) {
  plan <- stats::setNames(
    lapply(.DEPLOY_SIM_ARMS, function(arm) {
      path <- deploy_sim_template_path(dep, arm)
      cfg <- deploy_build_sim_template(dep, arm, date)
      if (!file.exists(path)) {
        return(list(
          arm = arm,
          path = path,
          status = "absent",
          cfg = cfg,
          diff = NULL
        ))
      }
      diff <- deploy_leaf_diff(
        yaml::read_yaml(path),
        cfg,
        ignore = .DEPLOY_SIM_IGNORED_LEAVES
      )
      list(
        arm = arm,
        path = path,
        cfg = cfg,
        diff = diff,
        status = if (all(diff$status == "match")) "unchanged" else "differs"
      )
    }),
    .DEPLOY_SIM_ARMS
  )
  deploy_assert_sim_pop_year(
    lapply(plan, function(p) p$cfg),
    "generated templates"
  )
  plan
}

#' The arms of a deploy must resolve one simulator.location_frame.pop_year:
#' the forecast arms hand off from the smoothing and the attribution
#' compares them, so a population-year difference between arms
#' would read as an effect. Set it once in the deploy YAML's `sim.shared`.
#' @param cfgs named list (arm -> parsed template).
#' @return the common spec (invisibly).
deploy_assert_sim_pop_year <- function(cfgs, what = "templates") {
  specs <- lapply(cfgs, function(cfg) {
    sim_pop_year_spec((cfg %||% list())$simulator)
  })
  vals <- vapply(specs, function(v) as.character(v), character(1))
  if (length(unique(vals)) > 1L) {
    stop(
      sprintf(
        paste0(
          "[deploy] the %s resolve different simulator.location_frame.pop_year ",
          "(%s). Set it once in the deploy YAML's sim.shared block."
        ),
        what,
        paste(sprintf("%s=%s", names(vals), vals), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(specs[[1L]])
}

deploy_sim_template_header <- function(
  config_path,
  arm,
  skeleton = deploy_sim_skeleton_path(arm)
) {
  generated_yaml_header(
    .DEPLOY_SIM_GENERATOR,
    c(
      sprintf("# spatPomp %s template of deploy %s", arm, config_path),
      sprintf(
        "#   = %s + deploy sim.shared + sim.%s + leaves derived from the deploy block.",
        skeleton,
        arm
      ),
      "# To change it, edit the deploy YAML's `sim:` block (or the skeleton), not this",
      "# file, and re-run (a differing existing template is refused unless --force):",
      sprintf(
        "#   Rscript analysis/utils/deploy_prod_run.R -C %s --stage gen_templates",
        config_path
      )
    )
  )
}

#' Write one planned arm (an element of deploy_plan_sim_templates()).
deploy_write_sim_template <- function(p, config_path) {
  dir.create(dirname(p$path), recursive = TRUE, showWarnings = FALSE)
  write_yaml_with_header(
    p$cfg,
    p$path,
    deploy_sim_template_header(config_path, p$arm)
  )
}

#' Assert an arm's template (on disk) agrees with the deploy block: run_tag,
#' neighbour-contagion gate, episode key, fits dir and -- forecast /
#' counterfactual -- the smoothing cache it hands off from. Used by the
#' smoothing / gen_forecast / gen_hindcast stages, so a template edited (or
#' hand-written) behind gen_templates' back is caught before it is consumed.
#'
#' @return the parsed template (invisibly).
deploy_check_sim_template <- function(dep, arm, path = dep$templates[[arm]]) {
  arm <- match.arg(arm, .DEPLOY_SIM_ARMS)
  if (!file.exists(path)) {
    stop(
      "[deploy] ",
      arm,
      " template not found: ",
      path,
      ". Generate it with ",
      "--stage gen_templates.",
      call. = FALSE
    )
  }
  t <- yaml::read_yaml(path)
  sim <- t$simulator %||% list()
  got_tag <- sim$run_tag %||% "<none>"
  if (!identical(got_tag, dep$run_tags[[arm]])) {
    stop(
      "[deploy] ",
      arm,
      " template run_tag mismatch: template has '",
      got_tag,
      "', deploy expects '",
      dep$run_tags[[arm]],
      "'",
      call. = FALSE
    )
  }
  if (isTRUE(dep$nbr_contagion$enable) && !isTRUE(sim$nbr_contagion$enable)) {
    stop(
      "[deploy] deploy.nbr_contagion.enable is TRUE but the ",
      arm,
      " template '",
      path,
      "' has simulator.nbr_contagion.enable != true -- it would run WITHOUT ",
      "the neighbour-contagion term the deployed fits carry.",
      call. = FALSE
    )
  }
  ep_key <- dep$episodes$key %||% "lp_season"
  tpl_key <- (t$episodes %||% list())$key %||% "lp_season"
  if (!identical(tpl_key, ep_key)) {
    stop(
      "[deploy] deploy.episodes.key is '",
      ep_key,
      "' but the ",
      arm,
      " template '",
      path,
      "' declares episodes.key = '",
      tpl_key,
      "'. Regenerate it with --stage ",
      "gen_templates (it must be built on the same episode definition as the ",
      "deployed fits).",
      call. = FALSE
    )
  }
  dep_nz <- isTRUE(dep$episodes$national_zero_fill)
  tpl_nz <- isTRUE((t$episodes %||% list())$national_zero_fill)
  if (!identical(tpl_nz, dep_nz)) {
    stop(
      "[deploy] deploy.episodes.national_zero_fill is ",
      tolower(dep_nz),
      " but the ",
      arm,
      " template '",
      path,
      "' declares ",
      tolower(tpl_nz),
      ". Regenerate it with ",
      "--stage gen_templates.",
      call. = FALSE
    )
  }
  dep_nzt <- .episodes_nz_threshold(dep$episodes)
  tpl_nzt <- .episodes_nz_threshold(t$episodes)
  if (!identical(tpl_nzt, dep_nzt)) {
    stop(
      "[deploy] deploy.episodes.national_zero_threshold is ",
      dep_nzt,
      " but the ",
      arm,
      " template '",
      path,
      "' declares ",
      tpl_nzt,
      ". Regenerate it with ",
      "--stage gen_templates.",
      call. = FALSE
    )
  }
  dep_spec <- .occ_panel_spec_from_episodes(dep$episodes)
  tpl_spec <- .occ_panel_spec_from_episodes(t$episodes)
  for (f in c("national_zero_rule", "at_risk_start", "adm1_overlap")) {
    if (!identical(tpl_spec[[f]], dep_spec[[f]])) {
      stop(
        "[deploy] deploy.episodes.",
        f,
        " is '",
        dep_spec[[f]],
        "' but the ",
        arm,
        " template '",
        path,
        "' declares '",
        tpl_spec[[f]],
        "'. Regenerate it with ",
        "--stage gen_templates.",
        call. = FALSE
      )
    }
  }
  fd <- (sim$grouping %||% list())$fits_dir
  if (!is.null(fd) && !identical(fd, dep$fits$dest)) {
    stop(
      "[deploy] ",
      arm,
      " template grouping.fits_dir is '",
      fd,
      "' but deploy.fits.dest ",
      "is '",
      dep$fits$dest,
      "'.",
      call. = FALSE
    )
  }
  # Registry partition: the template must carry this deploy's partition and
  # contagion registries. The historical (smoothing) partition may be the
  # registry itself (native per-unit caches) or per-country (stitch source).
  if (identical(dep$partition_by %||% "country", "registry")) {
    ps <- deploy_partition_source(dep)
    got <- (sim$filter %||% list())$partition_source %||% ""
    if (!identical(got, ps)) {
      stop(
        "[deploy] deploy.partition_source is '",
        ps,
        "' but the ",
        arm,
        " template '",
        path,
        "' has filter.partition_source = '",
        got,
        "'. Regenerate it with ",
        "--stage gen_templates.",
        call. = FALSE
      )
    }
    nps <- (sim$nbr_contagion %||% list())$partition_source %||% ""
    if (!identical(nps, deploy_nbr_partition_source(dep))) {
      stop(
        "[deploy] the ",
        arm,
        " template's nbr_contagion.partition_source is '",
        nps,
        "' but the deploy's contagion registry is '",
        deploy_nbr_partition_source(dep),
        "'.",
        call. = FALSE
      )
    }
    if (identical(arm, "historical")) {
      hp <- (sim$filter %||% list())$partition_by %||% "transmission_unit"
      if (!hp %in% c("registry", "country")) {
        stop(
          "[deploy] a registry-partitioned deploy needs the historical smoothing ",
          "partitioned by 'registry' (native per-unit caches) or 'country' (stitched), ",
          "not '",
          hp,
          "'.",
          call. = FALSE
        )
      }
    }
  }
  # One population year across the arms: a forward arm against
  # the historical template on disk.
  hist_path <- (dep$templates %||% list())[["historical"]]
  if (
    !identical(arm, "historical") &&
      !is.null(hist_path) &&
      file.exists(hist_path)
  ) {
    deploy_assert_sim_pop_year(
      stats::setNames(
        list(yaml::read_yaml(hist_path), t),
        c("historical", arm)
      ),
      sprintf("on-disk templates (%s)", path)
    )
  }
  if (!identical(arm, "historical")) {
    si <- (sim$filter %||% list())$smoothing_input %||% "<none>"
    if (!identical(si, dep$smoothing_input)) {
      stop(
        "[deploy] ",
        arm,
        " template filter.smoothing_input is '",
        si,
        "' but ",
        "deploy.smoothing_input is '",
        dep$smoothing_input,
        "' -- the forward arm ",
        "would hand off from another run's smoothing cache.",
        call. = FALSE
      )
    }
    # The arm's launch: horizon.start == forecast_window.start (== handoff_month
    # when set), checked here so a bad overlay fails at deploy time, not in 07b
    # (utils/sim_io.R::assert_arm_launch(), call-time; the same check as
    # make_forecast_configs.R's).
    launch <- assert_arm_launch(
      sim,
      sprintf("[deploy] %s template '%s'", arm, path)
    )
    # The analog anchor may differ from the handoff (a counterfactual handed off
    # from a later smoothed state with analog ENSO from the anchor): the
    # simulator only reports it (utils/simulations.R), so does the deploy.
    if (identical(arm, "counterfactual")) {
      lhs <- ((sim$enso %||% list())$counterfactual %||%
        list())$lenssen_horizon_start
      if (
        !is.null(lhs) &&
          nzchar(as.character(lhs)) &&
          as.Date(as.character(lhs)) != launch
      ) {
        message(sprintf(
          "[deploy] counterfactual template launches at %s, analog anchor (lenssen_horizon_start) %s.",
          format(launch),
          lhs
        ))
      }
      # Library-mode counterfactual analogs (deploy.enso_counterfactual.match;
      # read in-function, absent = legacy): the per-run cache is built by
      # gen_forecast, so it may never be the shared legacy cache.
      cf_match <- as.character(
        (dep[["enso_counterfactual"]] %||% list())[["match"]] %||% "legacy"
      )
      cp <- as.character(
        ((sim$enso %||% list())$counterfactual %||% list())[["cache_path"]] %||%
          "generated_data/lenssen_nino34_analogs.rds"
      )
      if (
        !identical(cf_match, "legacy") &&
          normalizePath(cp, mustWork = FALSE) ==
            normalizePath(
              "generated_data/lenssen_nino34_analogs.rds",
              mustWork = FALSE
            )
      ) {
        stop(
          "[deploy] counterfactual template '",
          path,
          "' reads the shared Lenssen cache '",
          cp,
          "', but deploy.enso_counterfactual.match is '",
          cf_match,
          "': overlay a per-run ",
          "simulator.enso.counterfactual.cache_path (sim.counterfactual), e.g. ",
          "generated_data/lenssen_<run>/lenssen_nino34_analogs.rds.",
          call. = FALSE
        )
      }
    }
  }
  invisible(t)
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

#' Build + persist the deploy provenance manifest (generated_data/deploy/<stem>/
#' deploy_manifest.qs). Mirrors results_write()'s stamp contract.
#'
#' @param dep         Resolved deploy block.
#' @param config_path The deploy YAML path (names the run tag).
#' @param stage       The orchestrator stage just completed.
#' @param inputs      data.frame(path, mtime, size) of deployed inputs (or NULL).
#' @param generated   Named list describing generated outputs (or NULL).
#' @return (invisibly) the stamped manifest object.
deploy_write_provenance <- function(
  dep,
  config_path,
  stage,
  inputs = NULL,
  generated = NULL
) {
  stamp <- list(
    deploy = dep,
    config_path = config_path,
    schema = dep$version %||% 1L,
    stage = stage,
    derived_at = format(Sys.time()),
    r_version = R.version.string,
    inputs = inputs,
    generated = generated
  )
  obj <- attach_deploy_provenance(stamp, stamp)
  path <- make_deploy_run_file(
    list(outdir = "generated_data/"),
    kind = "deploy_manifest",
    config_path = config_path
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(obj, path)
  } else {
    saveRDS(obj, sub("\\.qs$", ".rds", path))
  }
  invisible(obj)
}

# ---------------------------------------------------------------------------
# Run identity: required keys + the KEY=value contract read by bash drivers
# ---------------------------------------------------------------------------

# Deploy leaves that name a specific prod_run<N>. get_default_options()$deploy
# leaves them NULL on purpose: the run is identified ONLY by the deploy YAML, so
# a YAML that omits one must fail loudly here rather than silently inherit a
# stale run's directory from the defaults.
.DEPLOY_REQUIRED_KEYS <- c(
  "fits.occ_dir",
  "fits.dur_dir",
  "fits.mon_dir",
  "templates.historical",
  "templates.forecast",
  "templates.counterfactual",
  "run_tags.historical",
  "run_tags.forecast",
  "run_tags.counterfactual",
  "run_tags.hindcast",
  "outputs.forecast_config_outdir",
  "outputs.hindcast_config_outdir",
  "smoothing_input"
)

.deploy_leaf <- function(dep, dotted) config_get_leaf(dep, dotted) # config.R

#' Stop, naming every run-specific deploy key that is NULL or empty.
#'
#' @param dep         Resolved deploy block.
#' @param config_path The deploy YAML path (for the error message only).
#' @return (invisibly) TRUE.
deploy_require_keys <- function(dep, config_path = "<deploy config>") {
  empty <- vapply(
    .DEPLOY_REQUIRED_KEYS,
    function(k) {
      v <- .deploy_leaf(dep, k)
      is.null(v) || !length(v) || all(is.na(v)) || !any(nzchar(v))
    },
    logical(1)
  )
  if (any(empty)) {
    stop(
      "[deploy] ",
      config_path,
      " leaves run-specific key(s) unset: ",
      paste0("deploy.", .DEPLOY_REQUIRED_KEYS[empty], collapse = ", "),
      ". These have no default (they name a prod_run<N>); set them in the ",
      "deploy YAML.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The resolved run paths a bash driver needs, as an ordered named list of
#' scalar strings (KEY -> value). Pure: values come straight from `dep`; the
#' only derivations are the batch-dir prefix, the YAML stem and OBS_PANEL (the
#' observed survival panel matching `episodes.key`, via make_surv_dat_file() at
#' the default baseline), so bash never builds a path itself. Emitted by
#' `deploy_prod_run.R --stage resolve`.
#'
#' @param dep         Resolved deploy block.
#' @param config_path The deploy YAML path (DEPLOY_CONFIG / DEPLOY_TAG).
#' @return Named list of length-1 character values.
deploy_resolved_keys <- function(dep, config_path) {
  batch_dir <- function(b) file.path("analysis/configs", b)
  c(
    list(
      DEPLOY_CONFIG = config_path,
      DEPLOY_TAG = derive_config_run_tag(config_path),
      OCC_BATCH = dep$fits$occ_dir,
      DUR_BATCH = dep$fits$dur_dir,
      MON_BATCH = dep$fits$mon_dir,
      OCC_DIR = batch_dir(dep$fits$occ_dir),
      DUR_DIR = batch_dir(dep$fits$dur_dir),
      MON_DIR = batch_dir(dep$fits$mon_dir),
      FITS_DEST = dep$fits$dest,
      # Full registry labels; "" while a kl sweep has yet to derive the registry
      # (--stage gen_grouping), so resolve works before the `classes` phase.
      GROUPS = if (
        file.exists(registry_path(dep$grouping %||% "enso_grouping"))
      ) {
        paste(.deploy_group_full(dep, dep$fits$groups), collapse = " ")
      } else if (nzchar(.deploy_kl_sweep_batch(dep))) {
        ""
      } else {
        stop(
          "[deploy] registry not found: ",
          registry_path(dep$grouping %||% "enso_grouping"),
          " (deploy.grouping; no deploy.kl_sweep.batch derives it).",
          call. = FALSE
        )
      },
      HISTORICAL_TEMPLATE = dep$templates$historical,
      SMOOTHING_INPUT = dep$smoothing_input,
      FORECAST_CONFIG_DIR = dep$outputs$forecast_config_outdir,
      FORECAST_MANIFEST = dep$outputs$forecast_manifest,
      # Attribution manifest (--stage gen_attrib); "" = no attribution factual.
      FORECAST_ATTRIB_MANIFEST = as.character(
        dep$outputs$forecast_attrib_manifest %||% ""
      ),
      HINDCAST_CONFIG_DIR = dep$outputs$hindcast_config_outdir,
      HINDCAST_MANIFEST = dep$outputs$hindcast_manifest,
      # The counterfactual arm's Lenssen cache (per-run in library mode).
      COUNTERFACTUAL_CACHE = .deploy_cf_cache_path(dep),
      # Past-ENSO-event attribution (--stage gen_events); "" when unset. EVENTS
      # (the event ids) is the on/off switch run_pipeline.sh's sim phase reads.
      EVENTS = trimws(as.character(
        (dep$event_attribution %||% list())$events %||% ""
      )),
      EVENT_CONFIG_DIR = as.character(
        (dep$event_attribution %||% list())$config_outdir %||% ""
      ),
      EVENT_MANIFEST = as.character(
        (dep$event_attribution %||% list())$manifest %||% ""
      ),
      PARTITION_BY = dep$partition_by %||% "country",
      PARTITION_SOURCE = dep$partition_source %||% "",
      GROUPING = dep$grouping %||% "enso_grouping",
      KL_SWEEP_BATCH = .deploy_kl_sweep_batch(dep),
      KL_SWEEP_DIR = if (nzchar(.deploy_kl_sweep_batch(dep))) {
        batch_dir(.deploy_kl_sweep_batch(dep))
      } else {
        ""
      },
      KL_SWEEP_COMPARE_DIR = deploy_kl_sweep_compare_dir(dep),
      KL_SWEEP_CELL = as.character(dep$kl_sweep$cell %||% ""),
      KL_SWEEP_N_AGREE = as.character(as.integer(dep$kl_sweep$n_agree %||% 0L)),
      KL_SWEEP_MODE_TOL = as.character(as.numeric(
        dep$kl_sweep$mode_tol %||% 1.5
      ))
    ),
    list(
      # 07c/07d --obs-panel: must share the run's episode definition, or
      # assert_episode_coherence() refuses the hindcast validation.
      OBS_PANEL = as.character(make_surv_dat_file(list(
        baseline = get_default_options()$baseline,
        episodes = .episodes_block(
          .occ_panel_spec_from_episodes(dep$episodes),
          always = c("national_zero_fill", "national_zero_threshold")
        )
      )))
    )
  )
}

#' The counterfactual Lenssen cache the run's counterfactual template reads
#' (T-int 7(c)): the deploy overlay's simulator.enso.counterfactual.cache_path
#' (sim.counterfactual, else sim.shared), else the skeleton's
#' (analysis/configs/sim_templates/counterfactual.yml; the shared legacy cache).
.deploy_cf_cache_path <- function(dep) {
  leaf <- "simulator.enso.counterfactual.cache_path"
  ov <- dep$sim %||% list()
  v <- (ov$counterfactual %||% list())[[leaf]] %||%
    (ov$shared %||% list())[[leaf]]
  if (is.null(v)) {
    sk <- file.path(DEPLOY_SIM_SKELETON_DIR, "counterfactual.yml")
    if (file.exists(sk)) v <- config_get_leaf(yaml::read_yaml(sk), leaf)
  }
  as.character(v %||% "generated_data/lenssen_nino34_analogs.rds")
}

#' Print `KEY=value` lines on stdout -- the contract bash drivers parse with
#' `grep '^KEY=' | cut -d= -f2-` (same as 07a's TASKS_DIR= / N_TASKS=).
deploy_emit_keys <- function(keys) {
  for (k in names(keys)) {
    cat(k, "=", keys[[k]] %||% "", "\n", sep = "")
  }
  invisible(keys)
}

# ---------------------------------------------------------------------------
# Fit-batch generation (deploy_prod_run.R --stage gen_batches)
# ---------------------------------------------------------------------------

#' Load make_batch_configs.R into its own environment. The file redefines
#' `%||%` / get_default_config() and re-sources config.R at top level; an
#' isolated env keeps those out of the caller. Its sys.nframe() CLI guard makes
#' loading a no-op apart from the definitions.
deploy_batch_env <- function(path = "analysis/utils/make_batch_configs.R") {
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

#' Plan which batch config dirs to generate: the fit batches named by
#' `dep$fits$*_dir` and, when set, the kl class sweep `dep$kl_sweep$batch`
#' (track "kl").
#'
#' A batch dir that already holds a manifest.yml is left alone ("exists"): the
#' same dir is shared across rollouts (one run's dur / mon batches can back a
#' later occurrence refit), and a
#' regeneration would rewrite configs whose fits already exist. A registered
#' batch whose grouping registry (BATCH_REGISTRY spec$grouping) is not on disk
#' yet is "deferred": in a sweep rollout the fit dirs are generated after
#' --stage gen_grouping has written it (generate_batch() reads the registry to
#' expand one config per group).
#'
#' @param dep Resolved deploy block.
#' @param env Environment from deploy_batch_env().
#' @return data.frame(track, batch, dir, status, grouping) with status in
#'   "exists" | "generate" | "deferred" | "unregistered"; grouping = the spec's
#'   registry (NA when it has none or is unregistered).
deploy_plan_batches <- function(dep, env) {
  specs <- .deploy_track_specs()
  batches <- vapply(specs, function(s) dep$fits[[s$batch_key]], character(1))
  if (nzchar(.deploy_kl_sweep_batch(dep))) {
    batches <- c(batches, kl = .deploy_kl_sweep_batch(dep))
  }
  keep <- !duplicated(batches)
  batches <- batches[keep]
  dirs <- file.path(env$get_config_dir(), batches)
  grouping <- vapply(
    batches,
    function(b) {
      g <- env$BATCH_REGISTRY[[b]]$grouping
      if (is.null(g)) NA_character_ else as.character(g)
    },
    character(1)
  )
  reg_absent <- vapply(
    grouping,
    function(g) !is.na(g) && !file.exists(registry_path(g)),
    logical(1)
  )
  status <- ifelse(
    file.exists(file.path(dirs, "manifest.yml")),
    "exists",
    ifelse(
      !batches %in% names(env$BATCH_REGISTRY),
      "unregistered",
      ifelse(reg_absent, "deferred", "generate")
    )
  )
  data.frame(
    track = names(batches),
    batch = unname(batches),
    dir = dirs,
    status = status,
    grouping = unname(grouping),
    stringsAsFactors = FALSE
  )
}

#' Refuse a fit batch (occ/dur/mon row of deploy_plan_batches()) whose
#' BATCH_REGISTRY grouping registry is not `deploy.grouping`: its configs
#' would be expanded over another registry's groups than the run deploys
#' (fit collection, templates and results all read deploy.grouping). Rows
#' without a spec grouping (the kl sweep, ungrouped batches) are exempt.
#'
#' @return (invisibly) TRUE.
deploy_check_batch_groupings <- function(dep, plan) {
  fit_rows <- plan[
    plan$track %in% c("occ", "dur", "mon") & !is.na(plan$grouping),
    ,
    drop = FALSE
  ]
  bad <- fit_rows[
    fit_rows$grouping != (dep$grouping %||% "enso_grouping"),
    ,
    drop = FALSE
  ]
  if (nrow(bad)) {
    stop(
      "[deploy] fit batch(es) grouped by another registry than deploy.grouping '",
      dep$grouping %||% "enso_grouping",
      "': ",
      paste(
        sprintf("%s (BATCH_REGISTRY grouping '%s')", bad$batch, bad$grouping),
        collapse = ", "
      ),
      ". Set the batch spec's grouping (make_batch_configs.R) or deploy.grouping so they agree.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
