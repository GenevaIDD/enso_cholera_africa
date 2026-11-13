# 08b_smoothing_spatpomp.R — runner for the chunked spatPomp historical-
# smoothing pipeline. Companion to 07b_simulator_spatpomp.R (forward
# forecaster); this script does the historical fit (Phase A bpfilter over the
# full observed history + Phase B mon-track monthly case regeneration) for ONE
# (TU, replicate-chunk), and writes a partition for 08c_smoothing_compile.R.
#
# Two modes (mirrors 07b_simulator_v1.R):
#   --task task.yml   single-task mode (SLURM array slot): run only the chunk's
#                     replicates for the task's TU, write <output_part>/part.qs
#                     + a _SUCCESS sentinel.
#   (no --task)       local-loop mode: run every TU × replicate in-process (small
#                     Np) and write one partition per TU — for smoke tests.
#
# Usage:
#   Rscript analysis/08b_smoothing_spatpomp.R -C <cfg.yml> --task <task.yml>
#   Rscript analysis/08b_smoothing_spatpomp.R -C <cfg.yml> --n_reps 4 --Np 200 --tus hmmCE

suppressPackageStartupMessages({ library(optparse); library(qs2) })

option_list <- list(
  make_option(c("--config", "-C"), default = NULL, type = "character",
              help = "Path to YAML sim config"),
  make_option("--task",    default = NULL, type = "character",
              help = "Path to a task_NNNN.yml (single-task array mode)"),
  make_option("--n_reps",  default = NULL, type = "integer", help = "Local mode: n reps"),
  make_option("--Np",      default = NULL, type = "integer", help = "Particles per run"),
  make_option("--workers", default = NULL, type = "integer", help = "furrr workers"),
  make_option("--tus",     default = NULL, type = "character", help = "Local mode: TU subset"),
  make_option("--block_granularity", default = NULL, type = "character",
              help = "bpfilter block grouping: country (default) | adm1"),
  make_option("--partition_by", default = NULL, type = "character",
              help = "Top-level partition unit: transmission_unit (default) | country | enso_region | registry"),
  make_option("--partition_source", default = NULL, type = "character",
              help = "Partition registry under analysis/groupings/ (partition_by = registry)")
)
args <- parse_args(OptionParser(option_list = option_list))

source("analysis/utils/smoothing_core.R")

S   <- smoothing_make_config(args)
art <- smoothing_load_artefacts(S)

message(sprintf("[smoothing_run] N_REPS=%d Np=%d workers=%d measurement=%s grouped=%s partition=%s",
                S$N_REPS, S$NP, S$WORKERS, S$MEASUREMENT, S$GROUPED, S$PARTITION_BY))

# Partitions written without _SUCCESS because a bpfilter run failed (see
# .write_partition); a non-zero count makes this script exit 1.
.failed_parts <- character(0)

# Write one partition (.qs) + _SUCCESS sentinel for a single (TU, chunk).
.write_partition <- function(tu, chunk_id, results, sm, output_part) {
  events     <- smoothing_events_long(sm)
  # NULL when the mon track isn't configured for this run (opt-in, per
  # smoothing_events_mon_long()'s own docstring) — 08c_smoothing_compile.R
  # treats a missing/empty events_mon the same way.
  events_mon <- smoothing_events_mon_long(sm)
  rep_logliks <- data.frame(
    rep_id = vapply(results, function(r) as.integer(r$rep_id %||% NA_integer_), integer(1)),
    tu     = tu,
    logLik = vapply(results, function(r) as.numeric(r$logLik %||% NA_real_), numeric(1)),
    stringsAsFactors = FALSE
  )
  frame_meta <- if (!is.null(sm)) {
    # who_f: covered population share of the country, present when who_obs
    # is on (08c puts the WHO comparison on the national scale with it).
    # who_scale: the estimated per-country scale, only then.
    dplyr::select(sm$frame_sorted, gid, AFRO_region, country, dplyr::any_of(c("who_f", "who_scale"))) %>%
      dplyr::mutate(gid = as.character(gid)) %>%
      dplyr::distinct(gid, .keep_all = TRUE)
  } else NULL

  dir.create(output_part, recursive = TRUE, showWarnings = FALSE)
  part <- list(tu = tu, chunk_id = chunk_id, events = events, events_mon = events_mon,
               rep_logliks = rep_logliks, frame_meta = frame_meta)
  qs2::qs_save(part, file.path(output_part, "part.qs"))

  # Persist the FULL per-replicate smoothing trajectories (whole filter_traj
  # array, not just the terminal handoff slice) so the forward forecaster can
  # seed its initial conditions from the historical smoothing distribution.
  # Kept in a separate states.qs (not part.qs) so 08b's ribbon compile — which
  # binds every partition into memory at once — never has to load these large
  # arrays; the per-TU handoff cache is assembled TU-by-TU in 08b instead. Each
  # record already carries {tu, rep_id, draw_d, filt_months, filter_traj,
  # logLik} — exactly the fields load_smoothing_distribution() consumes.
  # The `_who` measurement's scratch slots (cm/ccm) were consumed by Phase B
  # above and are not part of the smoothing state: drop them so states.qs and
  # the handoff caches keep the sv11/sv12 layouts (spatpomp_core.R).
  results <- lapply(results, function(r) {
    r$filter_traj <- .sp_strip_scratch_slots(r$filter_traj); r })
  qs2::qs_save(results, file.path(output_part, "states.qs"))

  n_ok     <- sum(vapply(results, function(r) !is.null(r$filter_traj), logical(1)))
  # _SUCCESS only when every filter run succeeded: a failed replicate has no
  # trajectory for this TU, so the partition is incomplete. Without the
  # sentinel, 08_smoothing.sh --resume resubmits it and the afterok compile
  # does not run on a partial set. part.qs / states.qs stay for inspection.
  sentinel <- file.path(output_part, "_SUCCESS")
  if (n_ok == length(results)) {
    file.create(sentinel)
  } else {
    if (file.exists(sentinel)) file.remove(sentinel)
    .failed_parts <<- c(.failed_parts, output_part)
    warning(sprintf("[smoothing_run] TU '%s' %s: %d of %d filter runs failed; no _SUCCESS written.",
                    tu, chunk_id, length(results) - n_ok, length(results)), call. = FALSE)
  }
  n_ev     <- if (is.null(events))     0L else nrow(events)
  n_ev_mon <- if (is.null(events_mon)) 0L else nrow(events_mon)
  message(sprintf("[smoothing_run] TU '%s' %s: %d/%d filter runs ok, %d events (%d mon) → %s (+states.qs)",
                  tu, chunk_id, n_ok, length(results), n_ev, n_ev_mon,
                  file.path(output_part, "part.qs")))
}

# Run Phase A (chunk) + Phase B for one TU and a set of replicates.
.run_one_chunk <- function(tu, rep_ids, draw_idx, seeds, chunk_id, output_part) {
  pc <- smoothing_build_precomp(S, art, tu)[[tu]]
  if (is.null(pc)) { message(sprintf("[smoothing_run] TU '%s' empty — skipping.", tu)); return(invisible()) }

  # art$draws_dur is NULL unless simulator.duration_fit (ungrouped) /
  # simulator.grouping.duration_enable (grouped) is configured; NULL keeps
  # use_dur FALSE inside build_spatpomp_tu_c() and reproduces the previous,
  # duration-free smoothing behaviour exactly.
  # art$draws_mon follows the same opt-in rule (simulator.monthly_fit /
  # simulator.grouping.monthly_enable); NULL leaves the onset table only (no
  # case counts).
  results <- smoothing_filter_chunk(
    S, pc, art$draws_occ, tu,
    rep_ids = rep_ids, draw_idx = draw_idx, seeds = seeds, workers = S$WORKERS,
    draws_dur = art$draws_dur, draws_mon = art$draws_mon)

  sm <- extract_smoothing_one_tu(
    results, pc$frame_tu, pc$panel_tu,
    FILTER_START = S$FILTER_START, FILTER_END = S$FILTER_END, grouping = pc$grouping,
    draws_mon = art$draws_mon, cb_meta_mon = pc$cb_meta_mon,
    incidence_cap = S$INCIDENCE_CAP,
    cessation_exposure = S$CESSATION_EXPOSURE %||% "full")

  .write_partition(tu, chunk_id, results, sm, output_part)
}

if (!is.null(args$task)) {
  # ─── Single-task (array) mode ──────────────────────────────────────────────
  if (!file.exists(args$task)) stop("[smoothing_run] task not found: ", args$task)
  task <- yaml::read_yaml(args$task)
  message(sprintf("[smoothing_run] task %s: TU='%s' %s (%d reps)",
                  task$task_id, task$tu, task$chunk_id, length(task$replicate_ids)))

  # Task-carried scalars are authoritative for this run.
  S$NP          <- as.integer(task$Np %||% S$NP)
  S$MEASUREMENT <- task$measurement %||% S$MEASUREMENT
  S$BLOCK_GRANULARITY <- task$block_granularity %||% S$BLOCK_GRANULARITY
  S$PARTITION_BY <- task$partition_by %||% S$PARTITION_BY
  S$PARTITION_SOURCE <- task$partition_source %||% S$PARTITION_SOURCE
  # who_obs comes from the config (like covariate_feed / monthly_case_obs); the
  # task records whether it was on, so a task/config mismatch fails loudly.
  if (!is.null(task$who_obs_enable) &&
      !identical(isTRUE(task$who_obs_enable), isTRUE(S$WHO_OBS$enable)))
    stop("[smoothing_run] task ", task$task_id, " was generated with who_obs.enable = ",
         task$who_obs_enable, " but the config now says ", isTRUE(S$WHO_OBS$enable),
         "; regenerate the tasks (08a).")

  # Estimated WHO scale: this task must score the table 08a
  # wrote (and showed); a different estimate (inputs changed since 08a) is
  # refused. Inert in fixed mode (no digest on either side).
  got_dg <- if (is.null(art$who_scale)) NULL else whoScaleDigest(art$who_scale)
  if (!identical(got_dg, task$who_scale_digest))
    stop("[smoothing_run] task ", task$task_id, " was generated with WHO-scale digest ",
         task$who_scale_digest %||% "<none>", " but this run estimates ", got_dg %||% "<none>",
         "; regenerate the tasks (08a) and review who_scale.csv.")
  if (!identical(task$who_db_overlap %||% "score_both",
                 if (isTRUE(S$WHO_OBS$enable)) S$WHO_OBS[["db_overlap"]] %||% "score_both" else "score_both"))
    stop("[smoothing_run] task ", task$task_id, " was generated with who_obs.db_overlap = ",
         task$who_db_overlap %||% "score_both", " but the config now differs; regenerate the tasks (08a).")

  # Stamp partition column / collapse o_tu now that S is fully resolved.
  art <- smoothing_apply_partition(art, S)

  .run_one_chunk(
    tu          = task$tu,
    rep_ids     = as.integer(task$replicate_ids),
    draw_idx    = as.integer(task$draw_idx),
    seeds       = as.integer(task$seeds),
    chunk_id    = task$chunk_id,
    output_part = task$output_part)

} else {
  # ─── Local-loop mode (smoke parity) ────────────────────────────────────────
  smoothing_write_who_scale(S, art)   # estimated WHO scale; no-op in fixed mode
  art       <- smoothing_apply_partition(art, S)
  tu_subset <- if (!is.null(args$tus)) trimws(strsplit(args$tus, ",")[[1]]) else NULL
  tus       <- smoothing_tu_list(art$frame, tu_subset)
  # Seeds and draw indices are keyed on the FULL partition universe, so a
  # --tus re-run gives each TU the same streams as the full run.
  tus_all   <- smoothing_tu_list(art$frame)
  message(sprintf("[smoothing_run] local mode over %s partitions: %s",
                  S$PARTITION_BY, paste(tus, collapse = ", ")))

  parts_dir <- makeSmoothRunFile(list(outdir = S$outdir,
                                      simulator = list(run_tag = S$run_tag)),
                                 kind = "smooth_parts_dir")
  seed_rows <- smoothing_task_seeds(S, tus_all)

  for (tu in tus) {
    draw_idx_tu <- smoothing_draw_idx_for_tu(S, tu, tus_all, art$max_draw_s)
    seeds_tu    <- seed_rows$task_seed[seed_rows$tu == tu]
    output_part <- file.path(parts_dir, paste0("tu=", tu), "chunk_0001")
    .run_one_chunk(
      tu          = tu,
      rep_ids     = seq_len(S$N_REPS),
      draw_idx    = as.integer(draw_idx_tu),
      seeds       = as.integer(seeds_tu),
      chunk_id    = "chunk_0001",
      output_part = output_part)
  }
}

if (length(.failed_parts)) {
  message(sprintf("[smoothing_run] %d partition(s) incomplete (failed bpfilter runs), no _SUCCESS: %s",
                  length(.failed_parts), paste(.failed_parts, collapse = ", ")))
  quit(save = "no", status = 1L)
}
message("[smoothing_run] done.")
