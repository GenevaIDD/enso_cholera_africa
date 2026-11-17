#!/usr/bin/env Rscript
# make_event_attribution_configs.R — paired factual / counterfactual spatPomp
#   configs of the past-ENSO-event attribution, one pair per event and unit
#
# Usage:
#   Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_events
#   Rscript analysis/utils/make_event_attribution_configs.R \
#     --template <tpl> --counterfactual_template <cf tpl> --events <ids> \
#     --run_tag_base <base> --config_outdir <dir> --manifest <manifest.yml>
#
# Details:
#   Events: analysis/configs/enso_events_v2.yml, selected by
#   deploy.event_attribution.events. Factual = a hindcast with OBSERVED ENSO;
#   counterfactual = the same config + simulator.enso.counterfactual (Lenssen
#   analogs from a per-event cache, 03_prepare_lenssen_analogs.R). Both arms
#   hand off from the same smoothing state at the event's launch month.
#   Writes the configs and a manifest that 09's report_event_attrib reads.
#   Normally run by the deploy stage, which supplies every argument.
#
# See: generate_event_attribution_configs() for the pairing and the windows

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(lubridate)
})

source("analysis/utils/config.R") # load_grouping, registry_digest, assert_partition_labels
source("analysis/utils/base.R") # %||%
source("analysis/utils/paths.R") # make_run_file, read_deployed_lag, smoothing_cache_coverage, make_event_lenssen_dir
source("analysis/utils/enso.R") # ENSO_ANALOG_SETS, panel_max_enso_lead
source("analysis/utils/provenance.R") # generated_yaml_header, write_yaml_with_header, episode helpers
source("analysis/utils/anchor.R") # resolve_base_rate_paths()
source("analysis/utils/season_calendar.R") # per-country season ends, window rules
source("analysis/utils/make_hindcast_configs.R") # patch_hindcast_config() (CLI guarded)

option_list <- list(
  make_option(
    "--template",
    type = "character",
    default = NULL,
    help = "Forecast factual template (REQUIRED; deploy.templates.forecast)"
  ),
  make_option(
    "--counterfactual_template",
    type = "character",
    default = NULL,
    help = "Counterfactual template (REQUIRED; deploy.templates.counterfactual): its enso.counterfactual block"
  ),
  make_option(
    "--events",
    type = "character",
    default = "",
    help = "Comma-separated event ids of --events_registry (REQUIRED)"
  ),
  make_option(
    "--events_registry",
    type = "character",
    default = "analysis/configs/enso_events_v2.yml",
    help = "Event registry [default: %default]"
  ),
  make_option(
    "--analog_set",
    type = "character",
    default = "neutral",
    help = "Lenssen analog set of the counterfactual arm [default: %default]"
  ),
  make_option(
    "--lenssen_dir",
    type = "character",
    default = "generated_data/lenssen_events",
    help = "Root of the per-event Lenssen caches (<dir>/<event>/) [default: %default]"
  ),
  make_option(
    "--fits_dir",
    type = "character",
    default = "generated_data/enso_best_lag_fits",
    help = "Deployed best-fit dir (== grouping.fits_dir)"
  ),
  make_option(
    "--grouping",
    type = "character",
    default = "enso_grouping",
    help = "Fit registry"
  ),
  make_option(
    "--climate_panel",
    type = "character",
    default = "generated_data/sim_climate_panel_base-1990_2019.rds",
    help = "HISTORICAL continuous panel (observed climate + ENSO) covering every window"
  ),
  make_option(
    "--season_end",
    type = "character",
    default = "by_group",
    help = paste(
      "Season ends: by_group (--season_end_by_group; default, unchanged output)",
      "| by_country_from_panel (each country's min_seas - 1 from --climate_panel)"
    )
  ),
  make_option(
    "--season_end_by_group",
    type = "character",
    default = "",
    help = "Per-group season-end month 'grp=MM,...' (REQUIRED for every fit group under --season_end by_group)"
  ),
  make_option(
    "--season_min_window_months",
    type = "integer",
    default = 6L,
    help = "Accepted for symmetry with make_hindcast_configs.R; unused (the event window ends at the season containing the peak) [default: %default]"
  ),
  make_option(
    "--smoothing_input",
    type = "character",
    default = "",
    help = "Smoothing cache dir (default: inherit from the template)"
  ),
  make_option(
    "--episode_key",
    type = "character",
    default = "",
    help = "Episode definition written into every config ('' = template's)"
  ),
  make_option(
    "--partition_by",
    type = "character",
    default = "",
    help = "'' (template) | enso_region | registry"
  ),
  make_option(
    "--partition_source",
    type = "character",
    default = "",
    help = "Partition registry for --partition_by registry"
  ),
  make_option(
    "--nbr_partition_source",
    type = "character",
    default = NULL,
    help = "Neighbour-contagion partition registry [default: --partition_source]"
  ),
  make_option(
    "--run_tag_base",
    type = "character",
    default = NULL,
    help = "Base run_tag (REQUIRED); <base>_<event>_<unit>[_cf-<set>]"
  ),
  make_option(
    "--config_outdir",
    type = "character",
    default = NULL,
    help = "Directory for the generated configs (REQUIRED)"
  ),
  make_option(
    "--manifest",
    type = "character",
    default = NULL,
    help = "Event manifest path (REQUIRED)"
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "generated_data/",
    help = "Simulator output root (per-run merged dirs)"
  ),
  make_option(
    "--param_draw_mode",
    type = "character",
    default = "",
    help = "Posterior-draw mode of both arms: '' (patch_hindcast_config's random) | random | cycle | inherited | joint"
  ),
  make_option(
    "--n_replicates",
    type = "integer",
    default = 0L,
    help = "Replicates of both arms; 0 = the template's simulator.n_replicates [default: %default]"
  )
)

# ---------------------------------------------------------------------------
# Pure helpers (tests/test_event_attribution_configs.R)
# ---------------------------------------------------------------------------

#' load_enso_events
#'
#' The selected events of the registry, validated: phase, months, launch before
#' the peak, the annual rule's `relaunch_of` + 12 months. Named list event id ->
#' list(id, phase, label, launch, peak (Dates), peak_season_year, launch_rule,
#' relaunch_of, peaks, peak_season_years, sub_windows).
#'
#' A multi-winter event (enso_events_v2.yml) gives `peaks` and
#' `peak_season_years` (same length, strictly increasing) instead of `peak` /
#' `peak_season_year`, and `sub_windows` (none | per_winter; per_winter needs
#' two peaks or more). Its `peak` is then the LAST peak (it sets the
#' whole-event window end) and `peak_season_year` the FIRST winter; a
#' single-peak event has peaks = peak, peak_season_years = peak_season_year,
#' sub_windows "none".
load_enso_events <- function(path, ids) {
  reg <- yaml::read_yaml(path)$events
  if (!length(reg)) {
    stop("[events] no `events:` in ", path, call. = FALSE)
  }
  ids <- trimws(ids)
  ids <- ids[nzchar(ids)]
  if (!length(ids)) {
    stop(
      "[events] no event selected (deploy.event_attribution.events).",
      call. = FALSE
    )
  }
  bad <- setdiff(ids, names(reg))
  if (length(bad)) {
    stop(
      "[events] unknown event id(s) in ",
      path,
      ": ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  first_of_month <- function(x, what, id) {
    d <- as.Date(as.character(x %||% ""))
    if (is.na(d) || format(d, "%d") != "01") {
      stop(
        sprintf(
          "[events] %s: %s '%s' is not a first-of-month date.",
          id,
          what,
          x %||% ""
        ),
        call. = FALSE
      )
    }
    d
  }
  out <- lapply(stats::setNames(nm = ids), function(id) {
    e <- reg[[id]]
    if (!e$phase %in% c("elnino", "lanina")) {
      stop(
        sprintf("[events] %s: phase must be elnino or lanina.", id),
        call. = FALSE
      )
    }
    launch <- first_of_month(e$launch, "launch", id)
    multi <- !is.null(e[["peaks"]])
    if (multi && !is.null(e[["peak"]])) {
      stop(
        sprintf("[events] %s: give either `peak` or `peaks`, not both.", id),
        call. = FALSE
      )
    }
    peaks <- if (multi) {
      do.call(c, lapply(e[["peaks"]], first_of_month, what = "peaks", id = id))
    } else {
      first_of_month(e[["peak"]], "peak", id)
    }
    if (length(peaks) > 1L && any(diff(peaks) <= 0)) {
      stop(
        sprintf("[events] %s: peaks must be strictly increasing.", id),
        call. = FALSE
      )
    }
    peak <- peaks[length(peaks)]
    if (launch >= peaks[1L]) {
      stop(
        sprintf(
          "[events] %s: launch %s is not before the peak %s.",
          id,
          launch,
          peaks[1L]
        ),
        call. = FALSE
      )
    }
    rule <- e$launch_rule %||% "crossing"
    if (!rule %in% c("crossing", "annual")) {
      stop(
        sprintf("[events] %s: launch_rule must be crossing or annual.", id),
        call. = FALSE
      )
    }
    if (identical(rule, "annual")) {
      ref <- reg[[e$relaunch_of %||% ""]]
      if (is.null(ref)) {
        stop(
          sprintf("[events] %s: launch_rule annual needs relaunch_of.", id),
          call. = FALSE
        )
      }
      if (as.Date(ref$launch) %m+% months(12L) != launch) {
        stop(
          sprintf(
            "[events] %s: launch %s is not %s's launch + 12 months.",
            id,
            launch,
            e$relaunch_of
          ),
          call. = FALSE
        )
      }
    }
    psys <- if (multi) {
      if (!is.null(e[["peak_season_year"]])) {
        stop(
          sprintf(
            "[events] %s: a `peaks` event gives `peak_season_years`, not `peak_season_year`.",
            id
          ),
          call. = FALSE
        )
      }
      as.integer(unlist(e[["peak_season_years"]] %||% NA))
    } else {
      as.integer(e[["peak_season_year"]] %||% NA)
    }
    if (!length(psys) || anyNA(psys)) {
      stop(
        sprintf(
          "[events] %s: peak_season_year%s is required.",
          id,
          if (multi) "s" else ""
        ),
        call. = FALSE
      )
    }
    if (length(psys) != length(peaks)) {
      stop(
        sprintf(
          "[events] %s: %d peaks but %d peak_season_years.",
          id,
          length(peaks),
          length(psys)
        ),
        call. = FALSE
      )
    }
    if (length(psys) > 1L && any(diff(psys) <= 0)) {
      stop(
        sprintf(
          "[events] %s: peak_season_years must be strictly increasing.",
          id
        ),
        call. = FALSE
      )
    }
    sw <- as.character(e[["sub_windows"]] %||% "none")
    if (length(sw) != 1L || !sw %in% c("none", "per_winter")) {
      stop(
        sprintf("[events] %s: sub_windows must be none or per_winter.", id),
        call. = FALSE
      )
    }
    if (identical(sw, "per_winter") && length(peaks) < 2L) {
      stop(
        sprintf(
          "[events] %s: sub_windows per_winter needs two peaks or more.",
          id
        ),
        call. = FALSE
      )
    }
    ev <- list(
      id = id,
      phase = e$phase,
      label = e$label %||% id,
      launch = launch,
      peak = peak,
      peak_season_year = psys[1L],
      launch_rule = rule,
      relaunch_of = e$relaunch_of,
      peaks = peaks,
      peak_season_years = psys,
      sub_windows = sw
    )
    # Event-level constraint scope: `other_winters`
    # overrides deploy.enso_counterfactual.other_winters for this event only.
    # Absent = the deploy's rule; the element is then not added, so an event
    # without the key loads exactly as before.
    ow <- e[["other_winters"]]
    if (!is.null(ow)) {
      ow <- as.character(ow)
      if (length(ow) != 1L || !ow %in% c("neutral", "any")) {
        stop(
          sprintf("[events] %s: other_winters must be neutral or any.", id),
          call. = FALSE
        )
      }
      ev$other_winters <- ow
    }
    ev
  })
  out
}

#' event_counterfactual_cfg
#'
#' The counterfactual matcher settings of one event: deploy.enso_counterfactual
#' (enso_counterfactual_cfg()) with the event's registry `other_winters`, when it
#' sets one, in place of the deploy's. "any" constrains only the event's own
#' (primary) winters: the later winters its windows read are still classified
#' and recorded (cache stat_matrix; role read_unconstrained) but label no
#' trajectory "mixed". An event without the key gets the deploy's settings.
event_counterfactual_cfg <- function(cf, ev) {
  cf <- enso_counterfactual_cfg(cf)
  ow <- ev[["other_winters"]]
  if (is.null(ow)) {
    return(cf)
  }
  cf$other_winters <- ow
  enso_counterfactual_cfg(cf)
}

#' event_group_window
#'
#' One fit group's (or, under by_country_from_panel, one country's) window for
#' one event: the season end is the first `season_end_month` on or after the
#' peak month, extended by the group's largest positive deployed lag
#' (utils/season_calendar.R::event_window_end()). Returns list(season_end,
#' horizon_end (Dates), max_lag).
event_group_window <- function(launch, peak, season_end_month, lags) {
  launch <- as.Date(launch)
  peak <- as.Date(peak)
  sem <- as.integer(season_end_month)
  if (length(sem) != 1L || is.na(sem) || sem < 1L || sem > 12L) {
    stop("event_group_window: season_end_month must be 1-12.", call. = FALSE)
  }
  w <- event_window_end(peak, sem, lags)
  if (w$horizon_end <= launch) {
    stop(
      sprintf(
        "event_group_window: window end %s is not after the launch %s.",
        w$horizon_end,
        launch
      ),
      call. = FALSE
    )
  }
  w
}

#' event_sub_windows
#'
#' The consecutive per-winter sub-windows of a multi-winter event for one fit
#' group (or, under by_country_from_panel, one country): sub-window k ends at
#' event_group_window(launch, peaks[k], season_end_month, lags)$horizon_end
#' (the cholera season containing peak k + the largest positive deployed lag);
#' k = 1 starts at the launch, k + 1 the month after k's end. Ends must strictly
#' increase (two peaks in one season would give an empty sub-window): stops
#' otherwise. The last end is the whole-event window end of the last peak.
#' @return data.frame(k, start, end) (Dates)
event_sub_windows <- function(launch, peaks, season_end_month, lags) {
  launch <- as.Date(launch)
  ends <- do.call(
    c,
    lapply(as.Date(peaks), function(p) {
      event_group_window(launch, p, season_end_month, lags)$horizon_end
    })
  )
  .event_sub_window_frame(
    launch,
    ends,
    sprintf("season end month %d", as.integer(season_end_month))
  )
}

#' Consecutive sub-windows from their ends: k = 1 starts at the launch, k + 1
#' the month after k's end; ends must strictly increase.
.event_sub_window_frame <- function(launch, ends, what = "") {
  if (length(ends) > 1L && any(diff(ends) <= 0)) {
    stop(
      sprintf(
        "event_sub_windows: sub-window ends %s do not increase (%s): two peaks fall in one cholera season.",
        paste(format(ends), collapse = ", "),
        what
      ),
      call. = FALSE
    )
  }
  starts <- c(
    as.Date(launch),
    if (length(ends) > 1L) ends[-length(ends)] %m+% months(1L)
  )
  data.frame(k = seq_along(ends), start = starts, end = ends)
}

#' The sub-window ids of an event: "winter_<peak season year>".
event_sub_window_ids <- function(ev) {
  sprintf("winter_%d", as.integer(ev$peak_season_years))
}

#' Months the counterfactual sampler draws ENSO for (wrapper of
#' utils/enso.R::lenssen_n_lead_months()).
event_lenssen_lead_months <- function(launch, panel_end, max_lead) {
  lenssen_n_lead_months(launch, panel_end, max_lead)
}

#' The 03_prepare_lenssen_analogs.R command building one event's cache
#' (wrapper of utils/enso.R::lenssen_build_command(); the legacy spec gives the
#' former string exactly).
event_lenssen_command <- function(
  ev,
  lenssen_dir,
  lead_months,
  cf = NULL,
  ...
) {
  lenssen_build_command(
    ev$launch,
    lead_months,
    ev$peak_season_year,
    make_event_lenssen_dir(lenssen_dir, ev$id),
    cf = cf,
    ...
  )
}

#' event_lenssen_status
#'
#' Whether an event's Lenssen cache exists and fits the event: meta
#' horizon_start = launch, trajectories long enough, analogs classified at the
#' event's peak DJF, and the requested set non-empty. list(ok, reason, path, n_set).
event_lenssen_status <- function(
  ev,
  lenssen_dir,
  lead_months_needed,
  analog_set
) {
  dir <- make_event_lenssen_dir(lenssen_dir, ev$id)
  path <- make_lenssen_analog_file(dir)
  meta_path <- file.path(dir, "lenssen_nino34_analogs_meta.yml")
  res <- function(ok, reason, n_set = NA_integer_) {
    list(ok = ok, reason = reason, path = path, n_set = n_set)
  }
  if (!file.exists(path) || !file.exists(meta_path)) {
    return(res(FALSE, "cache absent"))
  }
  m <- yaml::read_yaml(meta_path)
  if (!identical(as.character(m$horizon_start), format(ev$launch))) {
    return(res(
      FALSE,
      sprintf(
        "cache horizon_start %s != launch %s",
        m$horizon_start,
        format(ev$launch)
      )
    ))
  }
  if (as.integer(m$lead_months %||% 0L) < lead_months_needed) {
    return(res(
      FALSE,
      sprintf(
        "cache lead_months %s < %d needed",
        m$lead_months,
        lead_months_needed
      )
    ))
  }
  psy <- m$set_definition$peak_season_year
  if (is.null(psy) || as.integer(psy) != ev$peak_season_year) {
    return(res(
      FALSE,
      sprintf(
        "cache sets classified at season year %s, not %d",
        psy %||% "auto",
        ev$peak_season_year
      )
    ))
  }
  n_set <- if (identical(analog_set, "all")) {
    as.integer(m$n_total_analogs %||% NA)
  } else {
    as.integer((m$set_definition$counts %||% list())[[analog_set]] %||% NA)
  }
  if (!is.na(n_set) && n_set < 1L) {
    return(res(FALSE, sprintf("analog set '%s' is empty", analog_set), n_set))
  }
  res(TRUE, "ok", n_set)
}

#' event_param_draw_mode
#'
#' Set the forward arms' posterior-draw mode after patch_hindcast_config(),
#' which writes "random". `mode` "" leaves the config as it is. Under "cycle"
#' replicate r runs draw ((r - 1) mod max_draws) + 1 and, with uniform handoff,
#' starts from smoothing replicate ((r - 1) mod n_smoothing) + 1. Under the
#' config's uniform handoff (simulator.forecast.handoff_weighting) a max_draws
#' sharing a factor with n_smoothing would tie each draw to a few smoothing
#' states (200 vs 100: one), so it is refused. "shuffled" pairs the states by a
#' seeded permutation per block of n_smoothing ids, so any max_draws passes.
#' Under "joint" replicate r starts from smoothing replicate
#' ((r - 1) mod n_smoothing) + 1 and runs its draw (07b checks the cache is
#' shared_cycle): the handoff must be uniform and max_draws equal n_smoothing.
#' @param cfg a patched event config
#' @param mode "" | "random" | "cycle" | "inherited" | "joint"
#' @param n_smoothing replicates of the smoothing cache the arms hand off from
#'   (its meta.yml `nrep`); read under "cycle" with uniform handoff and "joint"
#' @return cfg
event_param_draw_mode <- function(cfg, mode, n_smoothing = NULL) {
  mode <- as.character(mode %||% "")
  if (!nzchar(mode)) {
    return(cfg)
  }
  if (!mode %in% c("random", "cycle", "inherited", "joint")) {
    stop(
      "[make_event_attribution_configs] --param_draw_mode must be random, cycle, inherited or joint; got '",
      mode,
      "'.",
      call. = FALSE
    )
  }
  handoff <- as.character(
    cfg$simulator$forecast[["handoff_weighting"]] %||% "uniform"
  )
  if (identical(mode, "joint")) {
    if (!identical(handoff, "uniform")) {
      stop(
        "[make_event_attribution_configs] param_draw_mode joint pairs each replicate with its smoothing ",
        "replicate itself: the template's simulator.forecast.handoff_weighting must be uniform, not '",
        handoff,
        "'.",
        call. = FALSE
      )
    }
    S <- as.integer(cfg$simulator$max_draws %||% 200L)
    n <- as.integer(n_smoothing %||% NA)
    if (is.na(n) || n < 1L) {
      stop(
        "[make_event_attribution_configs] param_draw_mode joint needs the smoothing cache's replicate ",
        "count (meta.yml nrep) to check max_draws against it.",
        call. = FALSE
      )
    }
    if (S != n) {
      stop(
        sprintf(
          paste0(
            "[make_event_attribution_configs] param_draw_mode joint: max_draws %d but %d ",
            "smoothing replicates; the joint design runs one draw per smoothing replicate."
          ),
          S,
          n
        ),
        call. = FALSE
      )
    }
  }
  if (identical(mode, "cycle") && identical(handoff, "uniform")) {
    S <- as.integer(cfg$simulator$max_draws %||% 200L)
    n <- as.integer(n_smoothing %||% NA)
    if (is.na(n) || n < 1L) {
      stop(
        "[make_event_attribution_configs] param_draw_mode cycle needs the smoothing cache's replicate ",
        "count (meta.yml nrep) to check max_draws against it.",
        call. = FALSE
      )
    }
    gcd <- function(a, b) if (b == 0L) a else gcd(b, a %% b)
    if (gcd(S, n) != 1L) {
      ok <- S
      while (ok > 1L && gcd(ok, n) != 1L) {
        ok <- ok - 1L
      }
      stop(
        sprintf(
          paste0(
            "[make_event_attribution_configs] param_draw_mode cycle: max_draws %d shares a ",
            "factor with the %d smoothing replicates (gcd %d), so each draw would start ",
            "from only %d smoothing states. Use a max_draws coprime to %d (e.g. %d), or ",
            "simulator.forecast.handoff_weighting: shuffled."
          ),
          S,
          n,
          gcd(S, n),
          n %/% gcd(S, n),
          n,
          ok
        ),
        call. = FALSE
      )
    }
  }
  cfg$simulator$forecast$param_draw_mode <- mode
  cfg
}

#' event_n_replicates
#'
#' Set the forward arms' replicate count (deploy.event_attribution.n_replicates).
#' 0 or NULL leaves the template's simulator.n_replicates.
#' @param cfg a patched event config
#' @param n integer >= 0
#' @return cfg
event_n_replicates <- function(cfg, n) {
  if (is.null(n)) {
    return(cfg)
  }
  n <- suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop(
      "[make_event_attribution_configs] --n_replicates must be an integer >= 0; got '",
      paste(n, collapse = ","),
      "'.",
      call. = FALSE
    )
  }
  if (n == 0L) {
    return(cfg)
  }
  cfg$simulator$n_replicates <- n
  cfg
}

# ---------------------------------------------------------------------------
# generate_event_attribution_configs(opt) -- library entry point
# ---------------------------------------------------------------------------
# The two arms (in-sample):
#   factual        = a hindcast (patch_hindcast_config(): interior handoff,
#                    observed ENSO, historical climate panel)
#   counterfactual = the same config + simulator.enso.counterfactual (Lenssen
#                    piControl analogs matched to the launch ONI, one analog
#                    set, default `neutral`), cached per event by
#                    03_prepare_lenssen_analogs.R --horizon_start <launch>.
# Both arms share seed, handoff weighting and param_draw_mode
# (event_param_draw_mode()), so replicate r pairs across arms and F_r - C_r
# isolates the ENSO trajectory.
# Window per fit group: launch .. the end of the season holding the event
# peak + the group's largest positive deployed lag. A unit runs to its
# members' latest window; 09 censors each country at its own group's window.
generate_event_attribution_configs <- function(opt) {
  for (k in c(
    "template",
    "counterfactual_template",
    "run_tag_base",
    "config_outdir",
    "manifest"
  )) {
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) {
      stop(
        "[make_event_attribution_configs] --",
        k,
        " is required.",
        call. = FALSE
      )
    }
  }
  if (!opt$analog_set %in% ENSO_ANALOG_SETS) {
    stop(
      "[make_event_attribution_configs] unknown --analog_set '",
      opt$analog_set,
      "'.",
      call. = FALSE
    )
  }
  events <- load_enso_events(
    opt$events_registry,
    strsplit(opt$events, ",")[[1L]]
  )

  reg <- load_grouping(opt$grouping)
  groups <- sort(unique(reg$group))
  tpl <- yaml::read_yaml(opt$template)
  cf_tpl <- yaml::read_yaml(opt$counterfactual_template)
  cf_block <- cf_tpl$simulator$enso$counterfactual
  if (!identical(cf_block$scenario, "lenssen")) {
    stop(
      "[make_event_attribution_configs] the counterfactual template's enso.counterfactual.scenario is '",
      cf_block$scenario %||% "",
      "', not lenssen.",
      call. = FALSE
    )
  }

  ep_key <- if (nzchar(opt$episode_key %||% "")) {
    opt$episode_key
  } else {
    (tpl$episodes %||% list())$key %||% "lp_season"
  }
  ep_nz <- isTRUE((tpl$episodes %||% list())$national_zero_fill)
  ep_nzt <- .episodes_nz_threshold(tpl$episodes)
  ep_spec <- .occ_panel_spec_from_episodes(utils::modifyList(
    tpl$episodes %||% list(),
    list(key = ep_key)
  ))
  smoothing_input <- if (nzchar(opt$smoothing_input %||% "")) {
    opt$smoothing_input
  } else {
    tpl$simulator$filter$smoothing_input
  }
  if (is.null(smoothing_input) || !nzchar(smoothing_input)) {
    stop("[make_event_attribution_configs] no smoothing_input.", call. = FALSE)
  }
  # The smoothing cache's replicate count (meta.yml nrep, written by 08c): only
  # needed to check a cycle draw design under uniform handoff, or a joint one
  # (event_param_draw_mode()).
  n_smoothing <- if ((opt$param_draw_mode %||% "") %in% c("cycle", "joint")) {
    mp <- file.path(smoothing_input, "meta.yml")
    if (file.exists(mp)) {
      as.integer(yaml::read_yaml(mp)$nrep %||% NA)
    } else {
      NA_integer_
    }
  }

  # Season ends: every fit group needs one (no fallback length here).
  sem_map <- local({
    s <- trimws(opt$season_end_by_group)
    out <- integer(0)
    for (p in trimws(strsplit(s, ",")[[1L]])) {
      if (nzchar(p)) {
        kv <- trimws(strsplit(p, "=")[[1L]])
        if (length(kv) != 2L) {
          stop("Bad key=value in --season_end_by_group: ", p, call. = FALSE)
        }
        out[kv[1L]] <- as.integer(kv[2L])
      }
    }
    out
  })
  # --season_end by_country_from_panel: the season end is each COUNTRY's own
  # cholera-year end (min_seas - 1 of its modal gid calendar in
  # --climate_panel, utils/season_calendar.R), so the window is computed per
  # country with its group's lags; a group's horizon_end is the max over its
  # members, a unit runs to its latest member country, and countries.<iso>
  # carries the country's season end, min_seas mode, modal share and firm
  # horizon (the readers' censor). The manifest records season_end_source
  # (panel + digest) and season_ends_<manifest stem>.csv is written next to
  # it. The default, by_group, leaves the output unchanged and warns when a
  # group's members' ends span more than 4 months.
  season_mode <- as.character(opt[["season_end"]] %||% "by_group")
  if (!season_mode %in% c("by_group", "by_country_from_panel")) {
    stop(
      "[make_event_attribution_configs] --season_end must be by_group or by_country_from_panel; got '",
      season_mode,
      "'.",
      call. = FALSE
    )
  }
  by_country <- identical(season_mode, "by_country_from_panel")
  miss_sem <- setdiff(groups, names(sem_map))
  if (!by_country && length(miss_sem)) {
    stop(
      "[make_event_attribution_configs] --season_end_by_group lacks fit group(s): ",
      paste(miss_sem, collapse = ", "),
      call. = FALSE
    )
  }

  # Coverage: smoothing cache (handoff IC; smoothing-as-factual reading) and
  # climate panel (both arms read observed climate; the counterfactual sampler
  # draws ENSO for every panel month after the launch).
  cov <- smoothing_cache_coverage(smoothing_input)
  if (!file.exists(opt$climate_panel)) {
    stop(
      "[make_event_attribution_configs] climate panel not found: ",
      opt$climate_panel,
      call. = FALSE
    )
  }
  panel <- readRDS(opt$climate_panel)
  panel_end <- max(as.Date(panel$t_start))
  max_lead <- panel_max_enso_lead(panel)
  # The panel's gid calendar: each country's season end (by_country_from_panel),
  # or only the hemisphere-mix diagnostic (by_group).
  se_cal <- if (by_country) {
    panel_season_calendar(panel, opt$climate_panel)
  } else {
    tryCatch(
      panel_season_calendar(panel, opt$climate_panel),
      error = function(e) {
        message(
          "[make_event_attribution_configs] season-end diagnostic skipped: ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }
  se_tbl <- if (!is.null(se_cal)) season_end_by_country(se_cal) else NULL
  rm(panel)

  # Deployed lags + anchors per fit group (as make_hindcast_configs.R).
  group_info <- list()
  for (g in groups) {
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
    if (any(vapply(bl, is.null, logical(1)))) {
      message(sprintf(
        "[make_event_attribution_configs] SKIP group '%s' (missing deployed fit)",
        g
      ))
      next
    }
    group_info[[g]] <- list(
      bl = bl,
      lags = vapply(bl, `[[`, integer(1), "lag"),
      isos = reg$iso3[reg$group == g],
      base_rate_table_paths = resolve_base_rate_paths(
        lapply(bl, function(b) attr(readRDS(b$file), "opt"))
      )
    )
  }
  if (!length(group_info)) {
    stop(
      "[make_event_attribution_configs] no group has deployed occ/dur/mon fits in ",
      opt$fits_dir,
      call. = FALSE
    )
  }
  sem_c <- NULL # by_country_from_panel: ISO3 -> season-end month
  if (by_country) {
    all_isos <- unlist(lapply(group_info, `[[`, "isos"), use.names = FALSE)
    miss_c <- setdiff(all_isos, se_tbl$iso3)
    if (length(miss_c)) {
      stop(
        sprintf(
          paste0(
            "[make_event_attribution_configs] no gid of %s in the climate panel %s, so no ",
            "season end can be derived (--season_end by_country_from_panel)."
          ),
          paste(miss_c, collapse = ", "),
          opt$climate_panel
        ),
        call. = FALSE
      )
    }
    sem_c <- stats::setNames(
      se_tbl$season_end_month[match(all_isos, se_tbl$iso3)],
      all_isos
    )
  } else if (!is.null(se_tbl)) {
    for (w in season_end_spread_warnings(
      lapply(group_info, `[[`, "isos"),
      se_tbl
    )) {
      warning("[make_event_attribution_configs] ", w, call. = FALSE)
    }
  }
  fit_of <- stats::setNames(reg$group, reg$iso3)

  registry_mode <- identical(opt$partition_by %||% "", "registry")
  units <- if (!registry_mode) {
    lapply(stats::setNames(nm = names(group_info)), function(g) {
      list(members = group_info[[g]]$isos, groups = g)
    })
  } else {
    if (!nzchar(opt$partition_source %||% "")) {
      stop(
        "[make_event_attribution_configs] --partition_by registry requires --partition_source.",
        call. = FALSE
      )
    }
    preg <- load_grouping(opt$partition_source)
    assert_partition_labels(preg$group, opt$partition_source)
    fit_of <- stats::setNames(reg$group, reg$iso3)
    u <- list()
    for (lab in sort(unique(preg$group))) {
      mem <- preg$iso3[preg$group == lab]
      mem <- mem[mem %in% names(fit_of) & fit_of[mem] %in% names(group_info)]
      if (length(mem)) {
        u[[lab]] <- list(
          members = mem,
          groups = sort(unique(unname(fit_of[mem])))
        )
      }
    }
    u
  }

  date_str <- format(Sys.Date(), "%Y%m%d")
  out_fmt <- tpl$simulator$output$format %||% "qs"
  merged_dir_for <- function(run_tag, run_id) {
    make_run_file(
      list(
        outdir = opt$outdir,
        simulator = list(
          run_tag = run_tag,
          run_id = run_id,
          output = list(format = out_fmt)
        )
      ),
      kind = "spatpomp_merged_dir"
    )
  }
  dir.create(opt$config_outdir, recursive = TRUE, showWarnings = FALSE)
  aset <- opt$analog_set
  n_written <- 0L
  man_events <- list()
  build_cmds <- character(0)
  # Counterfactual analog matcher (deploy enso_counterfactual).
  cf_ev <- enso_counterfactual_cfg(opt$enso_counterfactual)
  season_rows <- NULL # by_country_from_panel: season_ends_<manifest stem>.csv

  for (ev in events) {
    launch_chr <- format(ev$launch)
    if (
      !is.null(cov) &&
        (ev$launch %m-% months(1L) < cov$start || ev$launch > cov$end)
    ) {
      stop(
        sprintf(
          paste0(
            "[make_event_attribution_configs] %s: the handoff IC month %s is outside the ",
            "smoothing coverage [%s, %s)."
          ),
          ev$id,
          format(ev$launch %m-% months(1L)),
          format(cov$start),
          format(cov$end)
        ),
        call. = FALSE
      )
    }
    # by_country_from_panel: one window per country (its own season end, its
    # group's lags); the group's window is the latest of its members'.
    cwin <- if (by_country) {
      lapply(stats::setNames(nm = names(sem_c)), function(iso) {
        event_group_window(
          ev$launch,
          ev$peak,
          sem_c[[iso]],
          group_info[[fit_of[[iso]]]]$lags
        )
      })
    } else {
      NULL
    }
    gwin <- lapply(stats::setNames(nm = names(group_info)), function(g) {
      if (!by_country) {
        return(event_group_window(
          ev$launch,
          ev$peak,
          sem_map[[g]],
          group_info[[g]]$lags
        ))
      }
      cw <- cwin[group_info[[g]]$isos]
      list(
        season_end = max(do.call(c, lapply(cw, `[[`, "season_end"))),
        horizon_end = max(do.call(c, lapply(cw, `[[`, "horizon_end"))),
        max_lag = cw[[1L]]$max_lag
      )
    })
    for (g in names(gwin)) {
      h <- gwin[[g]]$horizon_end
      # Closed window (utils/base.R): its last month must be in the panel.
      if (h > panel_end) {
        stop(
          sprintf(
            paste0(
              "[make_event_attribution_configs] %s / %s: window end %s is past the climate ",
              "panel (%s); the closed window's last month would not be simulated."
            ),
            ev$id,
            g,
            format(h),
            format(panel_end)
          ),
          call. = FALSE
        )
      }
      if (!is.null(cov) && h > cov$end) {
        warning(
          sprintf(
            paste0(
              "[make_event_attribution_configs] %s / %s: window end %s is past the smoothing ",
              "coverage (%s); the smoothing-as-factual reading is truncated there."
            ),
            ev$id,
            g,
            format(h),
            format(cov$end)
          ),
          call. = FALSE
        )
      }
    }
    # Per-winter sub-windows (enso_events_v2.yml): per group
    # (consecutive; under per-country season ends the group's k-th end is its
    # latest member's) and, by_country_from_panel, per country.
    sw_tab <- NULL
    if (identical(ev$sub_windows %||% "none", "per_winter")) {
      sw_c <- if (by_country) {
        lapply(stats::setNames(nm = names(sem_c)), function(iso) {
          event_sub_windows(
            ev$launch,
            ev$peaks,
            sem_c[[iso]],
            group_info[[fit_of[[iso]]]]$lags
          )
        })
      } else {
        NULL
      }
      sw_g <- lapply(stats::setNames(nm = names(group_info)), function(g) {
        if (!by_country) {
          return(event_sub_windows(
            ev$launch,
            ev$peaks,
            sem_map[[g]],
            group_info[[g]]$lags
          ))
        }
        ends <- do.call(
          c,
          lapply(seq_along(ev$peaks), function(k) {
            max(do.call(
              c,
              lapply(sw_c[group_info[[g]]$isos], function(x) x$end[k])
            ))
          })
        )
        .event_sub_window_frame(ev$launch, ends, g)
      })
      for (g in names(sw_g)) {
        if (sw_g[[g]]$end[nrow(sw_g[[g]])] != gwin[[g]]$horizon_end) {
          stop(
            sprintf(
              "[make_event_attribution_configs] %s / %s: the last sub-window ends %s, not at the window end %s.",
              ev$id,
              g,
              format(sw_g[[g]]$end[nrow(sw_g[[g]])]),
              format(gwin[[g]]$horizon_end)
            ),
            call. = FALSE
          )
        }
      }
      ids_sw <- event_sub_window_ids(ev)
      bounds <- function(x, k) {
        list(start = format(x$start[k]), end = format(x$end[k]))
      }
      sw_tab <- lapply(seq_along(ev$peaks), function(k) {
        out <- list(
          id = ids_sw[k],
          peak = format(ev$peaks[k]),
          peak_season_year = ev$peak_season_years[k],
          groups = lapply(sw_g, bounds, k = k)
        )
        if (by_country) {
          out$countries <- lapply(sw_c, bounds, k = k)
        }
        out
      })
    }
    lead_needed <- event_lenssen_lead_months(ev$launch, panel_end, max_lead)
    ev_cache <- NULL
    if (identical(cf_ev$match, "legacy") && length(ev$peak_season_years) > 1L) {
      stop(
        sprintf(
          paste0(
            "[make_event_attribution_configs] %s: a multi-winter event (peak seasons %s) needs ",
            "enso_counterfactual.match 'library' (the legacy builder labels one DJF)."
          ),
          ev$id,
          paste(ev$peak_season_years, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (identical(cf_ev$match, "legacy")) {
      lst <- event_lenssen_status(ev, opt$lenssen_dir, lead_needed, aset)
      if (!isTRUE(lst$ok)) {
        build_cmds <- c(
          build_cmds,
          event_lenssen_command(ev, opt$lenssen_dir, lead_needed)
        )
      }
    } else {
      # Library mode: built or verified in-process under the
      # event's lenssen dir; every DJF the event windows read is classified
      # (primary = the event's peak winter; the others neutral under winters:
      # all_read; an event's registry `other_winters` may relax that for the
      # event alone, event_counterfactual_cfg()).
      ev_win <- do.call(
        rbind,
        lapply(names(units), function(un) {
          data.frame(
            group = units[[un]]$groups,
            unit = un,
            end = do.call(
              c,
              lapply(units[[un]]$groups, function(g) gwin[[g]]$horizon_end)
            ),
            stringsAsFactors = FALSE
          )
        })
      )
      ev_read_end <- enso_read_end(
        ev_win,
        lapply(group_info, `[[`, "lags"),
        cf_ev$read_span
      )
      ev_path <- make_lenssen_analog_file(make_event_lenssen_dir(
        opt$lenssen_dir,
        ev$id
      ))
      ev_cache <- ensure_lenssen_counterfactual_cache(
        ev_path,
        ev$launch,
        lead_needed + 1L,
        event_counterfactual_cfg(cf_ev, ev),
        opt$enso_continuation,
        ev_read_end,
        primary_year = ev$peak_season_years %||% ev$peak_season_year,
        sets = aset,
        primary_set = aset,
        oni_dir = opt$outdir %||% "generated_data/",
        label = sprintf("[make_event_attribution_configs] %s", ev$id)
      )
      ev_cache$read_end <- ev_read_end
      lst <- list(
        ok = TRUE,
        reason = ev_cache$status$reason,
        path = ev_path,
        n_set = as.integer(ev_cache$status$counts[[aset]])
      )
    }

    man_units <- list()
    for (un in names(units)) {
      U <- units[[un]]
      gi <- group_info[U$groups]
      brt <- unique(lapply(gi, `[[`, "base_rate_table_paths"))
      if (length(brt) > 1L) {
        stop(
          sprintf(
            "[make_event_attribution_configs] unit '%s' mixes groups with different anchors.",
            un
          ),
          call. = FALSE
        )
      }
      horizon_chr <- if (!by_country) {
        max(vapply(
          U$groups,
          function(g) format(gwin[[g]]$horizon_end),
          character(1)
        ))
      } else {
        # the unit's latest MEMBER country (a group's latest may sit in another
        # unit); the readers' censor only truncates, so it must cover them all
        ch_u <- vapply(
          U$members,
          function(iso) format(cwin[[iso]]$horizon_end),
          character(1)
        )
        h <- max(ch_u)
        if (any(as.Date(ch_u) > as.Date(h))) {
          stop(
            sprintf(
              "[make_event_attribution_configs] %s / %s: unit horizon %s does not cover every member.",
              ev$id,
              un,
              h
            ),
            call. = FALSE
          )
        }
        h
      }
      tag_f <- sprintf("%s_%s_%s", opt$run_tag_base, ev$id, un)
      rid_f <- sprintf("%s_%s_%s_factual", date_str, ev$id, un)
      tag_c <- sprintf("%s_cf-%s", tag_f, aset)
      rid_c <- sprintf("%s_%s_%s_cf-%s", date_str, ev$id, un, aset)
      path_f <- file.path(
        opt$config_outdir,
        sprintf("event_%s_%s.yml", ev$id, un)
      )
      path_c <- file.path(
        opt$config_outdir,
        sprintf("event_%s_%s_cf-%s.yml", ev$id, un, aset)
      )
      assert_basename_length(
        c(path_f, path_c, tag_f, tag_c, rid_f, rid_c),
        "make_event_attribution_configs"
      )

      cfg_f <- patch_hindcast_config(
        tpl,
        opt,
        tag_f,
        rid_f,
        U$members,
        launch_chr,
        horizon_chr,
        brt[[1L]]
      )
      cfg_c <- patch_hindcast_config(
        tpl,
        opt,
        tag_c,
        rid_c,
        U$members,
        launch_chr,
        horizon_chr,
        brt[[1L]]
      )
      cfg_f <- event_param_draw_mode(cfg_f, opt$param_draw_mode, n_smoothing)
      cfg_c <- event_param_draw_mode(cfg_c, opt$param_draw_mode, n_smoothing)
      cfg_f <- event_n_replicates(cfg_f, opt$n_replicates)
      cfg_c <- event_n_replicates(cfg_c, opt$n_replicates)
      cfg_c$simulator$enso$counterfactual <- utils::modifyList(
        cf_block,
        list(
          scenario = "lenssen",
          cache_path = lst$path,
          lenssen_horizon_start = launch_chr,
          analog_set = aset
        )
      )

      win_lines <- sprintf(
        "# [%s window end %s: season end %s + max lag %d (occ %d, dur %d, mon %d)]",
        U$groups,
        vapply(
          U$groups,
          function(g) format(gwin[[g]]$horizon_end),
          character(1)
        ),
        vapply(
          U$groups,
          function(g) format(gwin[[g]]$season_end),
          character(1)
        ),
        vapply(U$groups, function(g) gwin[[g]]$max_lag, integer(1)),
        vapply(gi, function(x) x$lags[["occ"]], integer(1)),
        vapply(gi, function(x) x$lags[["dur"]], integer(1)),
        vapply(gi, function(x) x$lags[["mon"]], integer(1))
      )
      if (by_country) {
        win_lines <- c(
          sub("season end", "latest member season end", win_lines),
          sprintf(
            "# [%s (%s) season end %s (cholera-year end month %d) -> firm end %s]",
            U$members,
            vapply(U$members, function(i) fit_of[[i]], character(1)),
            vapply(
              U$members,
              function(i) format(cwin[[i]]$season_end),
              character(1)
            ),
            vapply(U$members, function(i) sem_c[[i]], integer(1)),
            vapply(
              U$members,
              function(i) format(cwin[[i]]$horizon_end),
              character(1)
            )
          )
        )
      }
      if (!is.null(sw_tab)) {
        win_lines <- c(
          win_lines,
          vapply(
            sw_tab,
            function(w) {
              sprintf(
                "# [sub-window %s (peak %s): %s]",
                w$id,
                w$peak,
                paste(
                  sprintf(
                    "%s %s..%s",
                    U$groups,
                    vapply(w$groups[U$groups], `[[`, "", "start"),
                    vapply(w$groups[U$groups], `[[`, "", "end")
                  ),
                  collapse = "; "
                )
              )
            },
            character(1)
          )
        )
      }
      common <- c(
        sprintf(
          "# %s (%s): launch %s, peak %s; unit '%s' (%d countries; fit groups %s).",
          ev$label,
          ev$id,
          launch_chr,
          format(ev$peak),
          un,
          length(U$members),
          paste(U$groups, collapse = ", ")
        ),
        win_lines,
        sprintf("# smoothing_input: %s", smoothing_input)
      )
      write_yaml_with_header(
        cfg_f,
        path_f,
        generated_yaml_header(
          "analysis/utils/make_event_attribution_configs.R",
          c(
            "# Event attribution, FACTUAL arm: observed ENSO from the launch (hindcast seams).",
            common
          )
        )
      )
      write_yaml_with_header(
        cfg_c,
        path_c,
        generated_yaml_header(
          "analysis/utils/make_event_attribution_configs.R",
          c(
            sprintf(
              "# Event attribution, COUNTERFACTUAL arm: Lenssen analogs (set %s) from the launch.",
              aset
            ),
            common,
            sprintf("# Lenssen cache: %s", lst$path)
          )
        )
      )
      n_written <- n_written + 2L

      man_units[[un]] <- list(
        members = as.list(U$members),
        groups = as.list(U$groups),
        horizon_end = horizon_chr,
        factual_config = path_f,
        factual_run_tag = tag_f,
        factual_run_id = rid_f,
        factual_merged_dir = merged_dir_for(tag_f, rid_f),
        factual_launch = launch_chr,
        analog_sets = list(aset),
        primary_analog_set = aset,
        counterfactual_sets = stats::setNames(
          list(list(
            analog_set = aset,
            run_tag = tag_c,
            run_id = rid_c,
            config = path_c,
            merged_dir = merged_dir_for(tag_c, rid_c),
            launch = launch_chr
          )),
          aset
        )
      )
    }

    man_groups <- lapply(stats::setNames(nm = names(group_info)), function(g) {
      x <- group_info[[g]]
      out <- list(
        best_lag_occ = x$lags[["occ"]],
        best_lag_dur = x$lags[["dur"]],
        best_lag_mon = x$lags[["mon"]],
        max_lag = gwin[[g]]$max_lag,
        season_end_month = if (!by_country) sem_map[[g]],
        season_end = format(gwin[[g]]$season_end),
        horizon_end = format(gwin[[g]]$horizon_end),
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
      # by_country_from_panel: no group season end (each member's is under
      # countries.<iso>); season_end / horizon_end are the latest member's.
      if (by_country) {
        out$season_end_month <- NULL
      }
      out
    })
    countries <- list()
    for (un in names(units)) {
      for (iso in units[[un]]$members) {
        g <- reg$group[reg$iso3 == iso]
        countries[[iso]] <- if (!by_country) {
          list(
            group = g,
            unit = un,
            firm_horizon_end = format(gwin[[g]]$horizon_end),
            unit_horizon_end = man_units[[un]]$horizon_end
          )
        } else {
          r <- se_tbl[se_tbl$iso3 == iso, , drop = FALSE]
          list(
            group = g,
            unit = un,
            season_end_month = r$season_end_month,
            min_seas_mode = r$min_seas_mode,
            modal_share = r$modal_share,
            season_end = format(cwin[[iso]]$season_end),
            firm_horizon_end = format(cwin[[iso]]$horizon_end),
            unit_horizon_end = man_units[[un]]$horizon_end
          )
        }
        if (by_country) {
          r <- se_tbl[se_tbl$iso3 == iso, , drop = FALSE]
          season_rows <- rbind(
            season_rows,
            data.frame(
              iso3 = iso,
              group = g,
              unit = un,
              event = ev$id,
              launch = launch_chr,
              peak = format(ev$peak),
              season_end_month = r$season_end_month,
              min_seas_mode = r$min_seas_mode,
              modal_share = r$modal_share,
              n_gid = r$n_gid,
              n_distinct = r$n_distinct,
              group_season_end_month = if (g %in% names(sem_map)) {
                sem_map[[g]]
              } else {
                NA_integer_
              },
              season_end = format(cwin[[iso]]$season_end),
              firm_horizon_end = format(cwin[[iso]]$horizon_end),
              unit_horizon_end = man_units[[un]]$horizon_end,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
    man_events[[ev$id]] <- list(
      phase = ev$phase,
      label = ev$label,
      launch = launch_chr,
      peak = format(ev$peak),
      peak_season_year = ev$peak_season_year,
      launch_rule = ev$launch_rule,
      lenssen_cache = lst$path,
      lenssen_ready = isTRUE(lst$ok),
      lenssen_status = lst$reason,
      lenssen_lead_months = lead_needed,
      lenssen_n_set = lst$n_set,
      africa_max_horizon = max(vapply(
        man_units,
        `[[`,
        character(1),
        "horizon_end"
      )),
      groups = man_groups,
      units = man_units,
      countries = countries
    )
    # Library-mode cache identity: absent in legacy mode, so a
    # legacy event manifest regenerates identically.
    if (!is.null(ev_cache)) {
      man_events[[ev$id]] <- c(
        man_events[[ev$id]],
        list(
          lenssen_digest = ev_cache$provenance$digest,
          lenssen_spec_digest = ev_cache$provenance$spec_digest,
          lenssen_winters = c(
            list(
              primary = ev_cache$winters$primary,
              other = as.list(ev_cache$winters$other),
              partial = as.list(ev_cache$winters$partial %||% integer(0)),
              read_end = format(ev_cache$read_end)
            ),
            # the event-level constraint scope, only when the registry sets one
            if (!is.null(ev[["other_winters"]])) {
              list(other_winters = ev[["other_winters"]])
            }
          ),
          lenssen_set_counts = as.list(ev_cache$status$counts)
        )
      )
    }
    # Per-winter sub-windows (schema 3): absent for a single-window event, so a
    # v1-registry manifest regenerates identically.
    if (!is.null(sw_tab)) {
      man_events[[ev$id]] <- c(
        man_events[[ev$id]],
        list(
          peaks = as.list(format(ev$peaks)),
          peak_season_years = as.list(ev$peak_season_years),
          sub_windows = sw_tab
        )
      )
    }
    message(sprintf(
      "[make_event_attribution_configs] %-12s launch %s peak %s -> %s; Lenssen %s",
      ev$id,
      launch_chr,
      format(ev$peak),
      paste(
        sprintf(
          "%s=%s",
          names(man_units),
          vapply(man_units, `[[`, character(1), "horizon_end")
        ),
        collapse = " "
      ),
      lst$reason
    ))
  }

  # Schema 3 only when an event has per-winter sub-windows.
  has_sw <- any(vapply(
    man_events,
    function(e) !is.null(e$sub_windows),
    logical(1)
  ))
  manifest <- list(
    schema = if (has_sw) 3L else 2L,
    kind = "event_attribution",
    generated_at = format(Sys.time()),
    events_registry = opt$events_registry,
    grouping = opt$grouping,
    grouping_digest = registry_digest(opt$grouping),
    analog_set = aset,
    lenssen_dir = opt$lenssen_dir,
    smoothing_input = smoothing_input,
    climate_panel = opt$climate_panel,
    panel_end = format(panel_end),
    episode_key = ep_key,
    partition_by = if (registry_mode) {
      "registry"
    } else {
      (opt$partition_by %||% "")
    },
    partition_source = opt$partition_source %||% "",
    partition_digest = if (registry_mode) {
      registry_digest(opt$partition_source)
    } else {
      NULL
    },
    contagion_partition_source = opt$nbr_partition_source %||%
      opt$partition_source %||%
      "",
    # Only when the deploy sets one, so a default manifest is unchanged.
    param_draw_mode = if (nzchar(opt$param_draw_mode %||% "")) {
      opt$param_draw_mode
    } else {
      NULL
    },
    n_replicates = if (isTRUE(opt$n_replicates > 0L)) {
      as.integer(opt$n_replicates)
    } else {
      NULL
    },
    events = man_events
  )
  # Only under by_country_from_panel, so a by_group manifest is unchanged.
  if (by_country) {
    manifest$season_end <- season_mode
    manifest$season_end_source <- season_end_source(se_cal)
  }
  # Population year: both templates must agree; recorded only
  # when not the default.
  py_f <- sim_pop_year_spec(tpl$simulator)
  py_cf <- sim_pop_year_spec(cf_tpl$simulator)
  if (!identical(py_f, py_cf)) {
    stop(
      sprintf(
        paste0(
          "[make_event_attribution_configs] factual and counterfactual templates ",
          "resolve different simulator.location_frame.pop_year (%s vs %s)."
        ),
        format(py_f),
        format(py_cf)
      ),
      call. = FALSE
    )
  }
  manifest <- manifest_add_provenance(
    manifest,
    "pop_year",
    py_f,
    .SIM_POP_YEAR_DEFAULT
  )
  dir.create(dirname(opt$manifest), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(manifest, opt$manifest)
  if (!is.null(season_rows)) {
    message(sprintf(
      "[make_event_attribution_configs] season ends -> %s",
      write_season_ends_csv(season_rows, opt$manifest)
    ))
  }

  message(sprintf(
    "[make_event_attribution_configs] %d configs (%d events x %d units x 2 arms) -> %s",
    n_written,
    length(events),
    length(units),
    opt$config_outdir
  ))
  message(sprintf(
    "[make_event_attribution_configs] manifest -> %s",
    opt$manifest
  ))
  if (length(build_cmds)) {
    message(
      "[make_event_attribution_configs] Lenssen caches to build before the counterfactual runs:"
    )
    for (cmd in build_cmds) {
      message("  ", cmd)
    }
  }
  invisible(list(manifest = manifest, lenssen_build = build_cmds))
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  opt <- parse_args(OptionParser(option_list = option_list))
  generate_event_attribution_configs(opt)
}
