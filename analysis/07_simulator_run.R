# 07_simulator_run.R — simulator backend dispatcher
#
# Reads the `simulator:` block from the YAML config and routes to the
# appropriate forward-simulator offspring, passing all arguments through
# unchanged (minus `--backend`, which is consumed here).
#
# Config must contain a `simulator:` block. Backend is read from
# `cfg$simulator$backend` (default "v1"); a `--backend` CLI flag overrides it.
#
# Dispatch is programmatic: the target script path is built with sprintf()
# from the literal backend string, e.g. backend = "spatpomp" ->
# analysis/07b_simulator_spatpomp.R. Adding a new backend is "add a file" —
# no dispatch table to edit here.
#
# Backends:
#   v1        -> analysis/07b_simulator_v1.R        (pure-R historical-replay kernel)
#   spatpomp  -> analysis/07b_simulator_spatpomp.R  (spatPomp bpfilter forecast)
#
# Run order: phase numbers are identifiers, not run order (analysis/PIPELINE.md).
# The spatpomp backend needs the phase-08 smoothing cache
# (simulator.filter.smoothing_input), so 08a/08b/08c run BEFORE it; v1 does not.
#
# Usage:
#   Rscript analysis/07_simulator_run.R -C analysis/configs/sim_smoketest/sim_smoke.yml
#   Rscript analysis/07_simulator_run.R -C <spatpomp cfg>            # backend: spatpomp in YAML
#   Rscript analysis/07_simulator_run.R -C <cfg> --backend spatpomp  # CLI override
#   Rscript analysis/07_simulator_run.R -C <cfg> --dry-run --task ...  # flags forwarded to offspring

library(yaml)

# Pull raw CLI args before optparse touches them
args <- commandArgs(trailingOnly = TRUE)

# Determine backend from config, then CLI override
backend <- NULL

config_idx <- which(args %in% c("--config", "-C"))
if (length(config_idx) > 0 && config_idx[1] < length(args)) {
  cfg_path <- args[config_idx[1] + 1]
  if (file.exists(cfg_path)) {
    cfg <- yaml::read_yaml(cfg_path)
    if (is.null(cfg$simulator)) {
      stop("Config must contain a `simulator:` block: ", cfg_path)
    }
    backend <- cfg$simulator$backend
  } else {
    stop("Config file not found: ", cfg_path)
  }
}

# CLI --backend override (and strip it from the forwarded args, since the
# offspring build their own optparse lists that error on unknown flags)
backend_idx <- which(args == "--backend")
if (length(backend_idx) > 0 && backend_idx[1] < length(args)) {
  backend <- args[backend_idx[1] + 1]
  args    <- args[-c(backend_idx[1], backend_idx[1] + 1)]
}

# Default
if (is.null(backend)) backend <- "v1"
backend <- tolower(trimws(backend))

script <- sprintf("analysis/07b_simulator_%s.R", backend)

if (!file.exists(script)) {
  stop(sprintf("Unknown simulator backend '%s' (no such script: %s)", backend, script))
}

message(sprintf("[07_simulator_run.R] backend = %s  ->  %s", backend, script))

status <- system2("Rscript", c(script, args))
quit(save = "no", status = status)
