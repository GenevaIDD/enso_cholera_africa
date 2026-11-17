# 07_simulator_run.R — simulator backend dispatcher
#
# Usage:
#   Rscript analysis/07_simulator_run.R \
#     -C analysis/configs/sim_smoketest/sim_smoke.yml
#   Rscript analysis/07_simulator_run.R -C <cfg> --backend spatpomp
#   Rscript analysis/07_simulator_run.R -C <cfg> --dry-run
#
# Details:
#   Reads cfg$simulator$backend (default "v1"; --backend overrides and is
#   consumed here) and runs analysis/07b_simulator_<backend>.R with all other
#   arguments passed through: v1 = pure-R historical-replay kernel, spatpomp =
#   spatPomp bpfilter forecast. The spatpomp backend needs the phase-08
#   smoothing cache (simulator.filter.smoothing_input), so 08a/08b/08c run
#   BEFORE it; v1 does not.
#
# See: analysis/PIPELINE.md (run order)

# 1. Preamble ------------------------------------------------------------------
library(yaml)

# 2. Resolve the backend -------------------------------------------------------
# Pull raw CLI args before optparse touches them
# Dispatch is programmatic: the target script path is built with sprintf()
# from the literal backend string, e.g. backend = "spatpomp" ->
# analysis/07b_simulator_spatpomp.R. Adding a new backend is "add a file" —
# no dispatch table to edit here.
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
  args <- args[-c(backend_idx[1], backend_idx[1] + 1)]
}

# Default
if (is.null(backend)) {
  backend <- "v1"
}
backend <- tolower(trimws(backend))

# 3. Dispatch ------------------------------------------------------------------
script <- sprintf("analysis/07b_simulator_%s.R", backend)

if (!file.exists(script)) {
  stop(sprintf(
    "Unknown simulator backend '%s' (no such script: %s)",
    backend,
    script
  ))
}

message(sprintf("[07_simulator_run.R] backend = %s  ->  %s", backend, script))

status <- system2("Rscript", c(script, args))
quit(save = "no", status = status)
