#!/usr/bin/env Rscript
# make_forecast_attrib_configs.R — attribution-factual configs of a forecast
#   rollout, launched at the counterfactual's handoff month, + their manifest
#
# Usage:
#   Rscript analysis/utils/deploy_prod_run.R -C <deploy.yml> --stage gen_attrib
#   Rscript analysis/utils/make_forecast_attrib_configs.R \
#     --manifest <forecast manifest.yml> --attrib_manifest <attrib.yml> \
#     [--dry-run] [--force]
#
# Details:
#   One factual config per simulation unit, written next to the factual
#   configs as forecast_<unit>_attrib.yml (the 07 submitter's `--arm attrib`),
#   plus an attribution manifest pointing the factual arm at those runs. 09
#   reads every attribution statistic from it (results.inputs.attrib_manifest);
#   the forecast manifest's factual stays the DISPLAYED forecast.
#
# See: generate_forecast_attrib_configs() for the pairing and overwrite policy

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/base.R") # %||%
source("analysis/utils/config.R") # config_flatten_leaves, config_set_leaf
source("analysis/utils/paths.R") # make_run_file (spatpomp_merged_dir)
source("analysis/utils/sim_io.R") # manifest_units, manifest_arm_launches, assert_arm_launch
source("analysis/utils/provenance.R") # generated_yaml_header, write_yaml_with_header

.ATTRIB_GENERATOR <- "analysis/utils/make_forecast_attrib_configs.R"
.ATTRIB_SUFFIX <- "_attrib"

# ---------------------------------------------------------------------------
# Pure helpers (tests/test_forecast_attrib_configs.R)
# ---------------------------------------------------------------------------

#' The launch months of a forecast manifest's two arms: list(factual,
#' counterfactual), each one Date. Refuses arms that do not launch at one month
#' across units and analog sets (the attribution factual must share the handoff
#' of every counterfactual set).
attrib_arm_launches <- function(man) {
  al <- manifest_arm_launches(man)
  one <- function(arm) {
    x <- unique(al$launch[al$arm == arm])
    if (length(x) != 1L) {
      stop(
        sprintf(
          "[attrib] the %s arm of the forecast manifest launches at %s; one month is required.",
          arm,
          if (length(x)) paste(format(x), collapse = ", ") else "no month"
        ),
        call. = FALSE
      )
    }
    x
  }
  list(factual = one("factual"), counterfactual = one("counterfactual"))
}

# Each attribution config is the generated factual config with three dates
# patched (simulator.horizon.start, filter.forecast_window.start,
# filter.handoff_month) and `_attrib` appended to run_tag and run_id. Nothing
# else changes, so the arm cannot drift from the factual template.
#' The attribution config: the factual config with its launch moved to
#' `launch` and `_attrib` appended to run_tag and run_id.
attrib_config_from_factual <- function(cfg, launch) {
  launch <- format(as.Date(launch))
  sim <- cfg$simulator
  if (is.null(sim)) {
    stop("[attrib] the factual config has no simulator block.", call. = FALSE)
  }
  if (!nzchar(as.character(sim$run_tag %||% ""))) {
    stop("[attrib] the factual config has no simulator.run_tag.", call. = FALSE)
  }
  sim$run_tag <- paste0(sim$run_tag, .ATTRIB_SUFFIX)
  sim$run_id <- paste0(sim$run_id %||% "", .ATTRIB_SUFFIX)
  sim$horizon$start <- launch
  sim$filter$forecast_window$start <- launch
  sim$filter$handoff_month <- launch
  cfg$simulator <- sim
  assert_arm_launch(
    cfg$simulator,
    sprintf("attribution config '%s'", sim$run_tag)
  )
  cfg
}

#' Leaves of two configs that differ, simulator.run_id excepted (it carries the
#' generation date). character(0) = equal.
attrib_config_diff <- function(a, b) {
  fa <- config_flatten_leaves(a)
  fb <- config_flatten_leaves(b)
  leaves <- setdiff(union(names(fa), names(fb)), "simulator.run_id")
  leaves[!vapply(leaves, function(l) identical(fa[[l]], fb[[l]]), logical(1))]
}

#' The merged dir 07b's aggregate writes for a run tag (spatPomp layout; the run
#' id is ignored, see make_forecast_configs.R::merged_dir_for()).
attrib_merged_dir <- function(run_tag, outdir = "generated_data/", fmt = "qs") {
  make_run_file(
    list(
      outdir = outdir,
      simulator = list(
        run_tag = run_tag,
        run_id = "",
        output = list(format = fmt)
      )
    ),
    kind = "spatpomp_merged_dir"
  )
}

#' The attribution manifest: the forecast manifest with each unit's (v2) or
#' group's (v1) factual run replaced by `runs[[entry]]` (list(config,
#' merged_dir)) launching at `launch`, `combined$factual` by `combined_factual`,
#' plus provenance keys. `runs = NULL` keeps the forecast manifest's factual runs
#' (aligned launches).
attrib_manifest_from <- function(
  man,
  runs,
  launch,
  display_launch,
  source_manifest,
  combined_factual = NULL
) {
  launch <- format(as.Date(launch))
  swap <- function(x, en) {
    if (is.null(runs)) {
      return(x)
    }
    r <- runs[[en]]
    if (is.null(r)) {
      stop(
        sprintf("[attrib] no attribution run for manifest entry '%s'.", en),
        call. = FALSE
      )
    }
    x$factual_config <- r$config
    x$factual_merged_dir <- r$merged_dir
    x$factual_launch <- launch
    x
  }
  out <- man
  if (length(man$units)) {
    for (en in names(man$units)) {
      out$units[[en]] <- swap(man$units[[en]], en)
    }
  } else {
    for (en in names(man$groups)) {
      out$groups[[en]] <- swap(man$groups[[en]], en)
    }
  }
  if (!is.null(runs) && !is.null(out$combined) && !is.null(combined_factual)) {
    out$combined$factual <- combined_factual
  }
  out$attribution_start <- launch
  out$kind <- "forecast_attribution"
  out$source_manifest <- source_manifest
  out$display_factual_launch <- format(as.Date(display_launch))
  out$attribution_factual <- if (is.null(runs)) "forecast" else "attrib"
  out
}

#' The combined-run merged dir of the attribution arm: the forecast manifest's
#' `combined$factual` with `_attrib` appended to its run-tag directory.
attrib_combined_dir <- function(combined_factual) {
  if (is.null(combined_factual) || !nzchar(combined_factual)) {
    return(NULL)
  }
  run_dir <- dirname(combined_factual)
  file.path(paste0(run_dir, .ATTRIB_SUFFIX), basename(combined_factual))
}

# ---------------------------------------------------------------------------
# Entry point (deploy_prod_run.R --stage gen_attrib, or the CLI below)
# ---------------------------------------------------------------------------

# The attribution factual launches at the counterfactual handoff, from the
# same smoothing replicate, seed, replicates and posterior-draw design, so
# factual - counterfactual isolates the ENSO trajectory replicate by replicate.
# The gap between the two factuals is the information component (09's
# report_attrib_information).
# Overwrite policy (as gen_templates): absent -> written; equal but for
# simulator.run_id -> kept; different -> refused unless --force. Configs of
# units no longer in the manifest are pruned. If the factual already launches
# at the counterfactual handoff, no config is written.
#' generate_forecast_attrib_configs
#'
#' @param o list(manifest, attrib_manifest, dry_run, force).
#' @return invisible list(launch, aligned, written, unchanged, pruned, configs).
generate_forecast_attrib_configs <- function(o) {
  tag <- "[make_forecast_attrib_configs]"
  if (!nzchar(o$manifest %||% "") || !file.exists(o$manifest)) {
    stop(
      tag,
      " forecast manifest not found: '",
      o$manifest %||% "",
      "' (run --stage gen_forecast first).",
      call. = FALSE
    )
  }
  if (!nzchar(o$attrib_manifest %||% "")) {
    stop(tag, " --attrib_manifest is required.", call. = FALSE)
  }
  man <- yaml::read_yaml(o$manifest)
  lau <- attrib_arm_launches(man)
  units <- manifest_units(man)
  dry <- isTRUE(o$dry_run)
  if (lau$counterfactual > lau$factual) {
    stop(
      tag,
      " the counterfactual launches after the factual; unsupported.",
      call. = FALSE
    )
  }
  aligned <- lau$counterfactual == lau$factual

  fac_paths <- vapply(
    units,
    function(u) as.character(u$factual_config %||% ""),
    character(1)
  )
  if (any(!nzchar(fac_paths) | !file.exists(fac_paths))) {
    stop(
      tag,
      " factual config(s) named by the manifest are missing: ",
      paste(
        fac_paths[!nzchar(fac_paths) | !file.exists(fac_paths)],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  cfg_dirs <- unique(dirname(fac_paths))
  if (length(cfg_dirs) != 1L) {
    stop(
      tag,
      " the factual configs sit in several directories: ",
      paste(cfg_dirs, collapse = ", "),
      call. = FALSE
    )
  }

  plan <- list()
  written <- character(0)
  unchanged <- character(0)
  refused <- list()
  if (!aligned) {
    for (en in names(units)) {
      fac <- yaml::read_yaml(fac_paths[[en]])
      cfg <- attrib_config_from_factual(fac, lau$counterfactual)
      path <- file.path(
        cfg_dirs,
        sub(
          "\\.yml$",
          paste0(.ATTRIB_SUFFIX, ".yml"),
          basename(fac_paths[[en]])
        )
      )
      assert_basename_length(
        c(path, cfg$simulator$run_tag),
        "make_forecast_attrib_configs"
      )
      status <- "absent"
      if (file.exists(path)) {
        d <- attrib_config_diff(yaml::read_yaml(path), cfg)
        status <- if (length(d)) "differs" else "unchanged"
        if (length(d)) refused[[path]] <- d
      }
      plan[[en]] <- list(
        path = path,
        cfg = cfg,
        status = status,
        merged_dir = attrib_merged_dir(
          cfg$simulator$run_tag,
          cfg$outdir %||% "generated_data/",
          cfg$simulator$output$format %||% "qs"
        )
      )
    }
    if (length(refused) && !isTRUE(o$force)) {
      for (p in names(refused)) {
        message(sprintf(
          "%s %s differs from the generated config on: %s",
          tag,
          p,
          paste(refused[[p]], collapse = ", ")
        ))
      }
      stop(
        tag,
        " ",
        length(refused),
        " existing attribution config(s) differ; re-run with --force to overwrite.",
        call. = FALSE
      )
    }
  }
  stale <- setdiff(
    list.files(
      cfg_dirs,
      pattern = paste0("^forecast_.*", .ATTRIB_SUFFIX, "\\.yml$"),
      full.names = TRUE
    ),
    vapply(plan, `[[`, "", "path")
  )

  message(sprintf(
    "%s factual launch %s, counterfactual launch %s%s",
    tag,
    format(lau$factual),
    format(lau$counterfactual),
    if (aligned) " (aligned: no attribution configs needed)" else ""
  ))
  for (en in names(plan)) {
    message(sprintf(
      "%s   %-8s %-9s %s",
      tag,
      en,
      plan[[en]]$status,
      plan[[en]]$path
    ))
  }
  if (length(stale)) {
    message(sprintf(
      "%s   stale    %s",
      tag,
      paste(basename(stale), collapse = ", ")
    ))
  }

  am <- attrib_manifest_from(
    man,
    runs = if (aligned) {
      NULL
    } else {
      lapply(plan, function(p) list(config = p$path, merged_dir = p$merged_dir))
    },
    launch = lau$counterfactual,
    display_launch = lau$factual,
    source_manifest = o$manifest,
    combined_factual = attrib_combined_dir(man$combined$factual)
  )
  am$generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  if (dry) {
    message(sprintf(
      "%s dry-run: would write %d config(s) and %s",
      tag,
      sum(vapply(plan, function(p) p$status != "unchanged", logical(1))),
      o$attrib_manifest
    ))
  } else {
    for (en in names(plan)) {
      p <- plan[[en]]
      if (identical(p$status, "unchanged")) {
        unchanged <- c(unchanged, p$path)
        next
      }
      write_yaml_with_header(
        p$cfg,
        p$path,
        generated_yaml_header(
          .ATTRIB_GENERATOR,
          c(
            sprintf(
              "# Attribution factual for unit '%s': the factual config %s",
              en,
              fac_paths[[en]]
            ),
            sprintf(
              "# launched at the counterfactual handoff %s (factual forecast launch %s).",
              format(lau$counterfactual),
              format(lau$factual)
            ),
            "# Only horizon.start, filter.forecast_window.start, filter.handoff_month, run_tag and run_id differ."
          )
        )
      )
      written <- c(written, p$path)
    }
    if (length(stale)) {
      file.remove(stale)
    }
    dir.create(
      dirname(o$attrib_manifest),
      recursive = TRUE,
      showWarnings = FALSE
    )
    write_yaml_with_header(
      am,
      o$attrib_manifest,
      generated_yaml_header(
        .ATTRIB_GENERATOR,
        c(
          sprintf(
            "# Attribution manifest of %s: the factual arm %s.",
            o$manifest,
            if (aligned) {
              "is the forecast manifest's own (launches already aligned)"
            } else {
              sprintf(
                "launches at %s (the counterfactual handoff)",
                format(lau$counterfactual)
              )
            }
          ),
          "# Read by 09 through results.inputs.attrib_manifest."
        )
      )
    )
    message(sprintf(
      "%s wrote %d config(s), %d unchanged, %d pruned; manifest %s",
      tag,
      length(written),
      length(unchanged),
      length(stale),
      o$attrib_manifest
    ))
  }
  invisible(list(
    launch = lau$counterfactual,
    aligned = aligned,
    written = written,
    unchanged = unchanged,
    pruned = if (dry) character(0) else stale,
    configs = vapply(plan, `[[`, "", "path"),
    config_dir = cfg_dirs
  ))
}

if (sys.nframe() == 0L) {
  opt <- optparse::parse_args(optparse::OptionParser(
    option_list = list(
      make_option(
        "--manifest",
        type = "character",
        default = "",
        help = "Forecast manifest written by make_forecast_configs.R (REQUIRED)"
      ),
      make_option(
        "--attrib_manifest",
        type = "character",
        default = "",
        help = "Attribution manifest to write (REQUIRED)"
      ),
      make_option(
        "--dry-run",
        action = "store_true",
        default = FALSE,
        dest = "dry_run",
        help = "Report the plan; write nothing"
      ),
      make_option(
        "--force",
        action = "store_true",
        default = FALSE,
        help = "Overwrite existing attribution configs that differ"
      )
    )
  ))
  generate_forecast_attrib_configs(opt)
}
