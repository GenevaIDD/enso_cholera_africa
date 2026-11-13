# kl_compile_batch.R — compile every Stan program of a brms_kl config batch once
#
# Non-step tool (no phase number; see analysis/utils/README.md §3). Runs
# `04_fit_run.R -C <config> --compile_only` for each config of a batch dir, in
# parallel. Configs that share a program (same sha, e.g. one ENSO basis at
# every K) meet at its compile lock (utils/dlnm_kernel.R::kl_with_compile_lock()):
# the first compiles, the others wait and get the finished executable. Configs
# whose fit already exists are skipped by the fit script. Submitted by
# bash/stages/04_fit_kl.sh as one job ahead of the fit array (afterok), so no
# array task compiles or races on an executable another task is writing.
#
# Run (repo root):
#   Rscript analysis/utils/kl_compile_batch.R -d analysis/configs/<batch> [--workers 12]

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "brms_kl batch config directory (manifest.yml is skipped)"),
  make_option("--workers", type = "integer", default = NA_integer_,
              help = "Parallel compile processes [default: one per config, at most SLURM_CPUS_PER_TASK]")
)))
if (is.null(opt$dir) || !dir.exists(opt$dir)) stop("--dir must be an existing batch config directory.")

cfgs <- sort(list.files(opt$dir, pattern = "\\.yml$", full.names = TRUE))
cfgs <- cfgs[basename(cfgs) != "manifest.yml"]
if (!length(cfgs)) stop("no configs in ", opt$dir)
cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
workers <- if (!is.na(opt$workers)) opt$workers else min(length(cfgs), if (is.na(cpus)) 4L else cpus)
workers <- max(1L, as.integer(workers))
message(sprintf("[kl_compile_batch] %d config(s) in %s, %d worker(s)", length(cfgs), opt$dir, workers))

one <- function(cfg) {
  log <- tempfile("kl_compile_", fileext = ".log")
  st <- system2("Rscript", c("analysis/04_fit_run.R", "-C", cfg, "--compile_only"), stdout = log, stderr = log)
  out <- readLines(log, warn = FALSE); unlink(log)
  list(cfg = basename(cfg), status = st, tail = utils::tail(out, 15L),
       compiled = grep("--compile_only", out, value = TRUE, fixed = TRUE))
}
res <- parallel::mclapply(cfgs, one, mc.cores = workers, mc.preschedule = FALSE)

failed <- FALSE
for (r in res) {
  if (inherits(r, "try-error") || is.null(r$status)) {
    message("[kl_compile_batch] worker error: ", as.character(r)); failed <- TRUE; next
  }
  if (r$status != 0L) {
    failed <- TRUE
    message(sprintf("[kl_compile_batch] FAILED %s (status %d):\n  %s", r$cfg, r$status,
                    paste(r$tail, collapse = "\n  ")))
  } else {
    message(sprintf("[kl_compile_batch] ok %s: %s", r$cfg,
                    if (length(r$compiled)) r$compiled[[1]] else "fit exists or nothing to compile"))
  }
}
if (failed) quit(save = "no", status = 1)
message("[kl_compile_batch] done")
