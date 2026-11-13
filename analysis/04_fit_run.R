# 04_fit_run.R — backend dispatcher
#
# Reads the model block from the YAML config and routes to the appropriate
# run script, passing all arguments through unchanged.
#
# Config must contain exactly one of `occurrence:`, `size:`, `duration:`,
# or `monthly:`. Backend is read from `cfg$<block>$backend`.
#
# For CLI-only invocations (no -C config), `--backend` implies the occurrence
# model and defaults to "brms_cb".
#
# Dispatch is programmatic: the target script path is built with sprintf()
# from the track ("occ"/"siz"/"dur"/"mon") and the literal backend string, e.g.
# backend = "brms_cb" -> analysis/04_fit_occ_brms_cb.R. Adding a new backend
# is "add a file" — no dispatch table to edit here.
#
# Usage:
#   Rscript analysis/04_fit_run.R -C analysis/configs/prod_run_occ_grouped_enso_singlelag/<config>.yml
#   Rscript analysis/04_fit_run.R -C analysis/configs/dur_smoketest/dur_smoke.yml
#   Rscript analysis/04_fit_run.R -C analysis/configs/mon_smoketest/mon_smoke.yml
#   Rscript analysis/04_fit_run.R --backend brms --country ETH --basehaz_poisson seasonal

library(yaml)
source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")

TRACK <- c(occurrence = "occ", size = "siz", duration = "dur", monthly = "mon")

# Pull raw CLI args before optparse touches them
args <- commandArgs(trailingOnly = TRUE)

# Determine model + backend from config, then CLI overrides
model   <- NULL
backend <- NULL
cfg     <- NULL

config_idx <- which(args %in% c("--config", "-C"))
if (length(config_idx) > 0 && config_idx[1] < length(args)) {
  cfg_path <- args[config_idx[1] + 1]
  if (file.exists(cfg_path)) {
    cfg <- yaml::read_yaml(cfg_path)
    if (!is.null(cfg$occurrence)) {
      model   <- "occurrence"
      backend <- cfg$occurrence$backend
    } else if (!is.null(cfg$size)) {
      model   <- "size"
      backend <- cfg$size$backend
    } else if (!is.null(cfg$duration)) {
      model   <- "duration"
      backend <- cfg$duration$backend
    } else if (!is.null(cfg$monthly)) {
      model   <- "monthly"
      backend <- cfg$monthly$backend
    } else {
      stop("Config must contain exactly one of `occurrence:`, `size:`, ",
           "`duration:` or `monthly:`: ", cfg_path)
    }
  } else {
    stop("Config file not found: ", cfg_path)
  }
}

# CLI --backend override (occurrence assumed for CLI-only invocations)
backend_idx <- which(args == "--backend")
if (length(backend_idx) > 0 && backend_idx[1] < length(args)) {
  backend <- args[backend_idx[1] + 1]
  if (is.null(model)) model <- "occurrence"
}

# Final defaults
if (is.null(model))   model   <- "occurrence"
if (is.null(backend)) backend <- get_default_options()[[model]]$backend
backend <- tolower(trimws(backend))

if (!model %in% names(TRACK)) {
  stop(sprintf("Unknown model '%s'. Must be one of: %s",
               model, paste(names(TRACK), collapse = ", ")))
}

script <- sprintf("analysis/04_fit_%s_%s.R", TRACK[[model]], backend)

if (!file.exists(script)) {
  stop(sprintf("Unknown backend '%s' for model '%s' (no such script: %s)",
               backend, model, script))
}

message(sprintf("[04_fit_run.R] model = %s  backend = %s  ->  %s", model, backend, script))

# --backend is the dispatcher's own flag; the workers do not define it.
if (length(backend_idx) > 0) {
  drop <- backend_idx
  drop <- c(drop, drop[drop < length(args)] + 1L)
  args <- args[-drop]
}
status <- system2("Rscript", c(script, args))
quit(save = "no", status = status)
