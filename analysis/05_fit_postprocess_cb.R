# 05_fit_postprocess_cb.R — extract crossbasis exposure-response surfaces
#
# Usage:
#   Rscript analysis/05_fit_postprocess_cb.R -C <config.yml>
#   Rscript analysis/05_fit_postprocess_cb.R -C <config.yml> --cumul --n_at 100
#
# Details:
#   Reads the brms_cb fit + cb_meta sidecar and applies dlnm::crosspred() per
#   draw to the posterior coefficients. Writes make_run_file kind
#   "cb_surfaces": a list with one surface object per dlnm variable.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/dlnm.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-c", "--country"),
    default = NULL,
    type = "character",
    help = "ISO3 code(s) or 'all' — overrides config country [default: from config]"
  ),
  make_option(
    c("--cumul"),
    default = FALSE,
    action = "store_true",
    help = "Also compute cumulative exposure-response curves [default: FALSE]"
  ),
  make_option(
    c("--n_at"),
    default = 50L,
    type = "integer",
    help = "Number of exposure grid points for prediction [default: 50]"
  ),
  make_option(
    c("--vars"),
    default = NULL,
    type = "character",
    help = "Comma-separated variable names to process (default: all in meta)"
  )
)

# Auto-detect occurrence vs size so size single-lag fits get surfaces too.
.detect_model <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  cfg <- NULL
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) {
      cfg <- argv[i + 1L]
      break
    }
  }
  if (is.null(cfg) || !file.exists(cfg)) {
    stop("--config / -C is required and must point to an existing YAML.")
  }
  raw <- yaml::read_yaml(cfg)
  if (!is.null(raw$occurrence)) {
    "occurrence"
  } else if (!is.null(raw$size)) {
    "size"
  } else if (!is.null(raw$duration)) {
    "duration"
  } else if (!is.null(raw$monthly)) {
    "monthly"
  } else {
    stop(
      "Config must have an occurrence:, size:, duration:, or monthly: block."
    )
  }
}
model_type <- .detect_model()

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = "country",
  model = model_type
)
stopifnot(opt[[model_type]]$backend == "brms_cb")
print_options(opt)

# File-kind prefix: size track uses size_fit / size_cb_meta / size_cb_surfaces;
# duration uses dur_fit / dur_cb_meta / dur_cb_surfaces; monthly uses
# mon_fit / mon_cb_meta / mon_cb_surfaces.
kpfx <- if (model_type == "size") {
  "size_"
} else if (model_type == "duration") {
  "dur_"
} else if (model_type == "monthly") {
  "mon_"
} else {
  ""
}

# 3. Locate files --------------------------------------------------------------

fit_file <- make_run_file(opt, kind = paste0(kpfx, "fit"))
meta_file <- make_run_file(opt, kind = paste0(kpfx, "cb_meta"))
out_file <- make_run_file(opt, kind = paste0(kpfx, "cb_surfaces"))

if (!file.exists(fit_file)) {
  stop(
    "Fit file not found: ",
    fit_file,
    "\n  Run 04_fit_run.R with backend: brms_cb first."
  )
}
if (!file.exists(meta_file)) {
  stop(
    "CB meta file not found: ",
    meta_file,
    "\n  Expected to be written by 04_fit_occ_brms_cb.R alongside the fit."
  )
}

message("Loading fit: ", fit_file)
brms_fit <- readRDS(fit_file)

message("Loading CB meta: ", meta_file)
cb_meta <- readRDS(meta_file)

# Which variables to process
vars <- if (!is.null(opt$vars)) {
  trimws(strsplit(opt$vars, ",")[[1]])
} else {
  NULL # all
}

# 4. Extract surfaces ----------------------------------------------------------
# Output (make_run_file kind = "cb_surfaces"): a named list, one entry per
# variable in the dlnm.variables block. Each entry contains:
#   $at, $lag, $matfit, $matlo, $mathi  — median + 95% CI surfaces
#   $allfit                             — full posterior surface array
#   $cumfit, $cumlo, $cumhi             — cumulative curves (if --cumul)
#   $meta                               — CB metadata for the variable

n_at <- as.integer(opt$n_at %||% 50L)
cumul <- isTRUE(opt$cumul)

message(sprintf(
  "Extracting CB surfaces: vars = %s  n_at = %d  cumul = %s",
  paste(vars %||% names(cb_meta), collapse = ", "),
  n_at,
  cumul
))

surfaces <- extract_cb_surface(
  fit = brms_fit,
  meta = cb_meta,
  vars = vars,
  at = NULL, # use stored range; override via at= if needed
  cumul = cumul
)

if (length(surfaces) == 0) {
  warning(
    "No surfaces extracted. Check that the fit contains the expected CB coefficient names."
  )
  quit(save = "no", status = 1)
}

message(sprintf(
  "Extracted surfaces for: %s",
  paste(names(surfaces), collapse = ", ")
))

saveRDS(surfaces, out_file)
message("Saved surfaces -> ", out_file)
