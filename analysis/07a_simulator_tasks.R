# 07a_simulator_tasks.R — materialise per-run task YAMLs for SLURM dispatch.
#
# Run before the v1 simulator array is submitted.
# Loads the config and the climate panel (to derive the TU list via the location
# frame), then writes one task_NNNN.yml per (TU, replicate-chunk) under
# <outdir>/sim/<run_tag>/tasks/.
#
# The fitted models are intentionally NOT loaded here — task planning
# only needs the TU membership, which comes from the panel + HMM cluster file.
#
# Emits on stdout (for bash capture):
#   TASKS_DIR=<path>
#   N_TASKS=<N>
#
# Usage:
#   Rscript analysis/07a_simulator_tasks.R -C analysis/configs/sim_smoketest/sim_smoke.yml

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/population.R")   # loadGidPopulation()/resolveGidPopulation() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
source("analysis/utils/composite_lp.R")   # mutual dependency with population.R (README: "source both together")
source("analysis/utils/simulator_core.R")
source("analysis/utils/anchor.R")

opt <- make_options_from_config(option_list = list(), model = NULL)

if (is.null(opt$simulator))
  stop("Config must contain a `simulator:` block. Got: ",
       paste(names(opt), collapse = ", "))

# Resolve run_tag (same logic as 09_simulator.R).
if (is.null(opt$simulator$run_tag)) {
  # .flag_val() is defined in utils/config.R (shared across the simulator run scripts).
  cfg_path <- .flag_val("-C") %||% .flag_val("--config")
  opt$simulator$run_tag <- if (!is.null(cfg_path))
    tools::file_path_sans_ext(basename(cfg_path))
  else
    format(Sys.time(), "sim_%Y%m%d_%H%M%S")
  message("[09_make_sim_tasks] derived run_tag = ", opt$simulator$run_tag)
}

tasks_dir <- makeRunFile(opt, kind = "sim_tasks_dir")
parts_dir <- makeRunFile(opt, kind = "sim_parts_dir")
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
message("[09_make_sim_tasks] tasks_dir = ", tasks_dir)

# ---------------------------------------------------------------------------
# Load climate panel + build location frame (for accurate TU list).
# Fits are not loaded — task planning only needs TU membership.
# ---------------------------------------------------------------------------
panel_path <- opt$simulator$climate_panel
if (is.null(panel_path) || !nzchar(panel_path))
  stop("simulator.climate_panel is required (got NULL/empty).")
if (!file.exists(panel_path))
  stop("simulator.climate_panel not found: ", panel_path)

message("[09_make_sim_tasks] reading climate panel: ", panel_path)
panel <- readRDS(panel_path)

# admin_level source for anchor resolution — no fits are loaded in
# this script (task planning only needs TU membership), so this read exists
# solely to supply `admin_level_data`; freed immediately after the call.
.surv_overlap <- tryCatch(
  readRDS(makeSurvDatFile(list(baseline = opt$baseline))),
  error = function(e) NULL)

lf_opts <- opt$simulator$location_frame %||% list()
.fpop <- simFramePopArgs(opt$simulator, opt)
.fpop$pop_window <- simPopWindow(opt$simulator)
frame <- build_location_frame(
  panel,
  hmm_clusters_path      = lf_opts$hmm_clusters_path,
  outbreak_metadata_path = lf_opts$outbreak_metadata_path,
  size_data_path         = lf_opts$size_data_path,
  # No fits are loaded here, so the three anchor tracks resolve off this
  # script's own `opt`.
  # NULL (fully inert) unless a track sets `offset_baserate: true`.
  base_rate_table_paths  = resolveBaseRatePaths(
                             list(occ = opt, dur = opt, mon = opt), lf_opts),
  admin_level_data       = .surv_overlap,
  # Population arguments in one place (simulator_core.R::simFramePopArgs():
  # location_frame.gid_pop_path / wpp_path / pop_year).
  gid_pop_path           = .fpop$gid_pop_path,
  wpp_path               = .fpop$wpp_path,
  pop_year               = .fpop$pop_year,
  pop_window             = .fpop$pop_window
)
rm(.surv_overlap)
message(sprintf("[09_make_sim_tasks] location frame: %d locations, %d TU(s) (incl. NA)",
                nrow(frame), length(unique(frame$tu))))

# Fallback: collapse all-NA TU to a synthetic bucket (same as 09_simulator.R).
if (all(is.na(frame$tu))) {
  message("[09_make_sim_tasks] all TUs are NA — collapsing to _ALL_")
  frame$tu <- "_ALL_"
}

# ---------------------------------------------------------------------------
# Build task plan.
# build_task_plan() / write_task_yamls() are the shared canonical planner
# defined in analysis/utils.R (sourced above) — the same functions the local
# furrr path (07b_simulator_v1.R) uses, so the SLURM task YAMLs are guaranteed
# identical (correct seed stride + posterior_strategy handling). The TU `frame`
# built here is fit-free (build_location_frame only, no fitted models loaded);
# build_task_plan() consumes it without loading fits. It resolves the parts
# output path internally via makeRunFile(opt, "sim_parts_dir") — identical to
# the `parts_dir` computed above.
# ---------------------------------------------------------------------------
plan  <- build_task_plan(opt, frame = frame)
paths <- write_task_yamls(plan, tasks_dir)
message(sprintf("[09_make_sim_tasks] wrote %d task YAML(s) to %s",
                length(paths), tasks_dir))

# ---------------------------------------------------------------------------
# Machine-readable output for bash capture (stdout only).
# ---------------------------------------------------------------------------
cat(sprintf("TASKS_DIR=%s\n", tasks_dir))
cat(sprintf("N_TASKS=%d\n", length(paths)))
