#!/usr/bin/env Rscript
# compact_sim_outputs.R — rewrite existing forward-simulation monthly
#   summaries in the compact format 2 (sim_io.R::sim_monthly_compact()).
#
# Usage:
#   Rscript analysis/utils/compact_sim_outputs.R \
#     --root generated_data/sim_spatpomp \
#     --pattern '<run glob>' --dry-run   # list candidates + sizes, read nothing
#   Rscript analysis/utils/compact_sim_outputs.R --root ... --pattern ... \
#     [--log out.csv]
#
# Details:
#   Non-step tool: converts files 07b / 07c already wrote, derives nothing.
#   A file is replaced only when the compact object reads back identical() to
#   the long table; a `<file>.format2` sidecar keeps the 09 stamps fresh.
#   Memory: one file at a time, about 3x its long table (read + compact +
#   check): a few GB for a forecast run, ~10 GB for a large event run.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(qs2)
})
source("analysis/utils/base.R")
source("analysis/utils/paths.R") # partition_sentinels()
source("analysis/utils/sim_io.R")

#' A merged / combined long table as compact blocks: itself when it is one
#' dense grid, else the blocks of the run's (already compact) partitions,
#' bound in aggregate_partitions()'s order -- kept only when they expand to
#' exactly this table. NULL otherwise.
compact_bound_table <- function(ms, run_dir) {
  cx <- sim_monthly_compact(ms)
  if (!is.null(cx) && identical(sim_monthly_expand(cx), ms)) {
    return(cx)
  }
  parts_dir <- file.path(run_dir, "parts")
  if (!dir.exists(parts_dir)) {
    return(NULL)
  }
  fs <- file.path(dirname(partition_sentinels(parts_dir)), "monthly_summary.qs")
  if (!length(fs) || !all(file.exists(fs))) {
    return(NULL)
  }
  parts <- lapply(fs, sim_io_read, fmt = "qs", expand = FALSE)
  if (!all(vapply(parts, is_sim_monthly_compact, logical(1)))) {
    return(NULL)
  }
  out <- sim_monthly_bind(parts)
  if (identical(sim_monthly_expand(out), ms)) out else NULL
}

# For every monthly_summary.qs / .rds under <root>/<run>/parts/** and
# <root>/<run>/merged/ (and combined/) of the runs matching --pattern:
#   * already compact                    -> left alone;
#   * long, a dense grid (one partition) -> compacted;
#   * long, a bind of partitions (merged/) -> the blocks of the run's compact
#     partitions, as aggregate_partitions() writes it, when they bind back
#     to exactly this table (partitions are converted first);
#   * anything else                      -> left long (reported).
# A file is replaced only when the compact object reads back identical() to
# the long table (sim_io_read()), through a temp file renamed over it; the
# original modification time and permissions are kept, and a sidecar
# `<file>.format2` records the original size and mtime: the 09 stamps that
# list these files (the hindcast step) record size too, and read the sidecar
# to stay fresh (results_core.R::.results_stamp_stale_reason()).
#' Convert one file; returns a one-row tibble (path, status, bytes_before, bytes_after).
compact_one_file <- function(path, dry_run = FALSE) {
  fmt <- if (grepl("\\.rds$", path)) "rds" else "qs"
  before <- file.size(path)
  row <- function(status, after = NA_real_) {
    tibble::tibble(
      path = path,
      status = status,
      bytes_before = before,
      bytes_after = after
    )
  }
  if (dry_run) {
    return(row("candidate"))
  }
  x <- sim_io_read(path, fmt, expand = FALSE)
  if (is_sim_monthly_compact(x)) {
    return(row("already_compact", before))
  }
  bound <- basename(dirname(path)) %in% c("merged", "combined")
  cx <- if (bound) {
    compact_bound_table(x, dirname(dirname(path)))
  } else {
    c1 <- sim_monthly_compact(x)
    if (!is.null(c1) && identical(sim_monthly_expand(c1), x)) c1
  }
  if (is.null(cx)) {
    return(row("left_long", before))
  }
  tmp <- paste0(path, ".compact_tmp")
  sim_io_save(cx, tmp, fmt)
  # the file as written must read back as the long table it replaces
  if (!identical(sim_io_read(tmp, fmt), x)) {
    unlink(tmp)
    return(row("verify_failed", before))
  }
  mt <- file.mtime(path)
  md <- file.mode(path)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    return(row("rename_failed", before))
  }
  Sys.chmod(path, md)
  Sys.setFileTime(path, mt)
  # the 09 stamps record size as well as mtime: the sidecar lets them recognise
  # this file as the one they recorded (results_core.R::.results_stamp_stale_reason())
  sim_monthly_write_sidecar(path, before, mt)
  row("compacted", file.size(path))
}

#' The monthly summaries of the runs under `root` matching `pattern`, each
#' run's partitions first (a merged file is rebuilt from its compact parts).
compact_candidates <- function(root, pattern) {
  runs <- Sys.glob(file.path(root, pattern))
  runs <- runs[dir.exists(runs)]
  unlist(
    lapply(runs, function(r) {
      f <- list.files(
        r,
        pattern = "^monthly_summary\\.(qs|rds)$",
        recursive = TRUE,
        full.names = TRUE
      )
      in_parts <- startsWith(f, file.path(r, "parts", ""))
      c(f[in_parts], f[!in_parts])
    }),
    use.names = FALSE
  )
}

if (sys.nframe() == 0L) {
  opt <- parse_args(OptionParser(
    option_list = list(
      make_option(
        "--root",
        default = "generated_data/sim_spatpomp",
        help = "simulator output root [default: %default]"
      ),
      make_option(
        "--pattern",
        default = NULL,
        help = "glob of run directories under --root (required), e.g. '*_prod_run*'"
      ),
      make_option(
        "--dry-run",
        action = "store_true",
        default = FALSE,
        dest = "dry_run",
        help = "list the candidate files and their sizes; read and write nothing"
      ),
      make_option(
        "--log",
        default = NULL,
        help = "write the per-file table to this CSV"
      )
    )
  ))
  if (is.null(opt$pattern) || !nzchar(opt$pattern)) {
    stop("--pattern is required (a glob of run directories).", call. = FALSE)
  }
  files <- compact_candidates(opt$root, opt$pattern)
  if (!length(files)) {
    stop(
      "no monthly_summary files under ",
      file.path(opt$root, opt$pattern),
      call. = FALSE
    )
  }
  message(sprintf(
    "[compact_sim_outputs] %d file(s), %.1f GB%s",
    length(files),
    sum(file.size(files)) / 1e9,
    if (opt$dry_run) " (dry run)" else ""
  ))
  res <- list()
  for (f in files) {
    r <- tryCatch(compact_one_file(f, opt$dry_run), error = function(e) {
      tibble::tibble(
        path = f,
        status = paste("error:", conditionMessage(e)),
        bytes_before = file.size(f),
        bytes_after = NA_real_
      )
    })
    message(sprintf(
      "  %-16s %8.1f -> %8.1f MB  %s",
      r$status,
      r$bytes_before / 1e6,
      r$bytes_after / 1e6,
      f
    ))
    res[[length(res) + 1L]] <- r
    invisible(gc())
  }
  res <- dplyr::bind_rows(res)
  s <- dplyr::count(res, status)
  message(
    "[compact_sim_outputs] ",
    paste(sprintf("%s %d", s$status, s$n), collapse = ", "),
    sprintf(
      " | %.2f -> %.2f GB",
      sum(res$bytes_before) / 1e9,
      sum(res$bytes_after, na.rm = TRUE) / 1e9
    )
  )
  if (!is.null(opt$log)) {
    utils::write.csv(res, opt$log, row.names = FALSE)
  }
  if (any(grepl("^(error|verify_failed|rename_failed)", res$status))) {
    quit(status = 1)
  }
}
