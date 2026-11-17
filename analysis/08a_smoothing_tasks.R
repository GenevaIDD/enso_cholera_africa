# 08a_smoothing_tasks.R — materialise per-(TU, replicate-chunk) task YAMLs for
#   the chunked spatPomp historical smoothing (SLURM array dispatch)
#
# Usage:
#   Rscript analysis/08a_smoothing_tasks.R -C <cfg.yml>
#   Rscript analysis/08a_smoothing_tasks.R -C <cfg.yml> --n_reps 100 \
#     --chunk_size 10
#
# Details:
#   Runs BEFORE 07b_simulator_spatpomp.R: 08c's handoff cache is the
#   forecast's smoothing_input. Writes one task_NNNN.yml per (TU, chunk) under
#   <outdir>/sim_spatpomp/<run_tag>/tasks/ and prints TASKS_DIR=<path> and
#   N_TASKS=<N> on stdout for bash capture.
#
# See: analysis/PIPELINE.md (run order)

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/smoothing_core.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("--config", "-C"),
    default = NULL,
    type = "character",
    help = "Path to YAML sim config"
  ),
  make_option(
    "--n_reps",
    default = NULL,
    type = "integer",
    help = "Number of independent bpfilter runs [config or default]"
  ),
  make_option(
    "--Np",
    default = NULL,
    type = "integer",
    help = "Particles per bpfilter run [config or default]"
  ),
  make_option(
    "--workers",
    default = NULL,
    type = "integer",
    help = "furrr workers per array task [config or default]"
  ),
  make_option(
    "--chunk_size",
    default = NULL,
    type = "integer",
    help = "Replicates per array task [config simulator.parallel.chunk_size or 10]"
  ),
  make_option(
    "--tus",
    default = NULL,
    type = "character",
    help = "Comma-separated TU subset (e.g. hmmW,hmmCE)"
  )
)
args <- parse_args(OptionParser(option_list = option_list))

# 3. Load artefacts ------------------------------------------------------------
# Unlike the 07a planner, this one DOES load the occurrence fits: it needs the
# true global `max_draw_s` (under grouping, the min S across all groups) to
# reproduce the single-job per-(TU, rep) occurrence draw indices exactly, so a
# chunked run recompiles to the same result as a single-job run.
S <- smoothing_make_config(args)
art <- smoothing_load_artefacts(S)
art <- smoothing_apply_partition(art, S) # stamp `part` / collapse o_tu in country mode
# Estimated WHO scale (who_obs.scale_mode: estimated): written
# once here, next to tasks/ (who_scale.csv + who_scale_meta.yml), shown before
# the array is submitted; its digest goes into every task and 08b re-checks it.
# NULL (nothing written, tasks unchanged) in fixed mode.
who_scale_digest <- smoothing_write_who_scale(S, art)

tu_subset <- if (!is.null(args$tus)) {
  trimws(strsplit(args$tus, ",")[[1]])
} else {
  NULL
}
tus <- smoothing_tu_list(art$frame, tu_subset)
# Seeds and draw indices are keyed on the FULL partition universe, so a
# --tus re-run gives each TU the same streams as the full run.
tus_all <- smoothing_tu_list(art$frame)
message(sprintf(
  "[smoothing_make_tasks] %s partitions: %s",
  S$PARTITION_BY,
  paste(tus, collapse = ", ")
))

# 4. Chunk size and output paths -----------------------------------------------
# Resolve chunk_size: CLI > config simulator.parallel.chunk_size > 10.
cfg_chunk <- NULL
if (!is.null(args$config)) {
  cfg <- yaml::read_yaml(args$config)
  cfg_chunk <- cfg$simulator$parallel$chunk_size
}
chunk_size <- args$chunk_size %||% cfg_chunk %||% 10L
chunk_size <- as.integer(chunk_size)
if (chunk_size < 1L) {
  stop("[smoothing_make_tasks] chunk_size must be >= 1")
}

# Path layout (rooted at <outdir>/sim_spatpomp/<run_tag>/).
tasks_dir <- make_smooth_run_file(
  list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
  kind = "smooth_tasks_dir"
)
parts_dir <- make_smooth_run_file(
  list(outdir = S$outdir, simulator = list(run_tag = S$run_tag)),
  kind = "smooth_parts_dir"
)
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
message("[smoothing_make_tasks] tasks_dir = ", tasks_dir)
message("[smoothing_make_tasks] parts_dir = ", parts_dir)

# 5. Build task plan -----------------------------------------------------------
# Per-(TU, rep) occurrence draw indices + task seeds (identical to the
# single-job assignment; see smoothing_core.R).
seed_rows <- smoothing_task_seeds(S, tus_all)

plan <- list()
i <- 0L
for (tu in tus) {
  draw_idx_tu <- smoothing_draw_idx_for_tu(S, tu, tus_all, art$max_draw_s) # length N_REPS
  chunks <- split(seq_len(S$N_REPS), ceiling(seq_len(S$N_REPS) / chunk_size))
  for (k in seq_along(chunks)) {
    i <- i + 1L
    ri <- chunks[[k]]
    seeds_tu <- seed_rows$task_seed[seed_rows$tu == tu] # length N_REPS, rep-ordered
    plan[[i]] <- list(
      task_id = sprintf("task_%04d", i),
      tu = tu,
      chunk_id = sprintf("chunk_%04d", k),
      replicate_ids = as.integer(ri),
      draw_idx = as.integer(draw_idx_tu[ri]),
      seeds = as.integer(seeds_tu[ri]),
      Np = as.integer(S$NP),
      measurement = S$MEASUREMENT,
      covariate_feed = S$COVARIATE_FEED %||% "latent",
      monthly_case_obs = S$MONTHLY_CASE_OBS %||% "none",
      who_obs_enable = isTRUE(S$WHO_OBS$enable),
      block_granularity = S$BLOCK_GRANULARITY %||% "country",
      partition_by = S$PARTITION_BY %||% "transmission_unit",
      partition_source = S$PARTITION_SOURCE %||% "",
      grouped = isTRUE(S$GROUPED),
      window = list(start = S$hist_window$start, end = S$hist_window$end),
      output_part = file.path(
        parts_dir,
        paste0("tu=", tu),
        sprintf("chunk_%04d", k)
      )
    )
    # Only when estimated / skip_who, so a default run's task YAMLs are unchanged.
    if (!is.null(who_scale_digest)) {
      plan[[i]]$who_scale_digest <- who_scale_digest
    }
    if (identical(S$WHO_OBS[["db_overlap"]], "skip_who")) {
      plan[[i]]$who_db_overlap <- "skip_who"
    }
  }
}

# 6. Write task YAMLs ----------------------------------------------------------
paths <- character(length(plan))
for (j in seq_along(plan)) {
  paths[j] <- file.path(tasks_dir, paste0(plan[[j]]$task_id, ".yml"))
  yaml::write_yaml(plan[[j]], paths[j])
}
message(sprintf(
  "[smoothing_make_tasks] wrote %d task YAML(s) (chunk_size=%d) to %s",
  length(paths),
  chunk_size,
  tasks_dir
))

cat(sprintf("TASKS_DIR=%s\n", tasks_dir))
cat(sprintf("N_TASKS=%d\n", length(paths)))
