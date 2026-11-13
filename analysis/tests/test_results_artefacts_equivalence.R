# test_results_artefacts_equivalence.R -- 09's manuscript artefacts are unchanged
# by a refactor.
#
# HARD gate, full suite only (see analysis/tests/run_all.R); SKIPs (exit 0) when
# no baseline exists. Guards refactors of the derivers in
# utils/results_derive.R: every artefact 09 writes must stay identical.
#
# Workflow:
#   1. On the pre-refactor code, run 09 once and copy its output:
#        Rscript analysis/09b_postprocess_results.R -C analysis/configs/results/manuscript.yml --redo_cache TRUE
#        cp -r generated_data/results/manuscript generated_data/results/manuscript_pre_b7c
#   2. After each refactor step:
#        Rscript analysis/tests/test_results_artefacts_equivalence.R --run
#      (--run re-derives with --redo_cache TRUE first; without it, the current
#      artefacts are compared as they are.)
#
# Per *.qs in the baseline dir: the payload (object minus its
# "results_provenance" attribute) must be identical(); the stamp must agree on
# the stamped results blocks, config path, schema token, globs and the input
# paths + sizes (derived_at, r_version and input mtimes are expected to move).
# On failure all.equal() is printed so a names-attribute-only difference is
# distinguishable from a value difference.
#
# A deliberate change to a kind's version, file name, payload or stamped keys
# invalidates the baseline: retake it after such a change.
#
# Reads generated_data/ (fits, simulator partitions, smoothing parts) -- the
# pipeline outputs 09 consumes. No model fit.

suppressPackageStartupMessages({ library(optparse); library(qs2) })
cli <- parse_args(OptionParser(option_list = list(
  make_option("--baseline-dir", dest = "baseline_dir",
              default = "generated_data/results/manuscript_pre_b7c"),
  make_option("--current-dir", dest = "current_dir",
              default = "generated_data/results/manuscript"),
  make_option("--config", default = "analysis/configs/results/manuscript.yml"),
  make_option("--run", action = "store_true", default = FALSE,
              help = "re-derive the current artefacts with 09 --redo_cache TRUE first"))))

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/provenance.R")   # getResultsProvenance()

cat("# 09 manuscript artefacts: current vs pre-refactor baseline\n")
base_files <- list.files(cli$baseline_dir, pattern = "\\.qs$")
if (!dir.exists(cli$baseline_dir) || !length(base_files)) {
  cat("  SKIP  no baseline artefacts in", cli$baseline_dir, "(see the header for how to make one)\n")
  quit(status = 0)
}

if (cli$run) {
  r <- run_rscript("analysis/09b_postprocess_results.R",
                   c("-C", cli$config, "--redo_cache", "TRUE"))
  expect(identical(r$status, 0L), "09b_postprocess_results.R re-derived without error")
}

strip <- function(x) { attr(x, "results_provenance") <- NULL; x }
stamp_core <- function(x) {
  s <- getResultsProvenance(x)
  if (is.null(s)) return(NULL)
  list(results = s$results, config_path = s$config_path, schema = s$schema,
       globs = s$globs,
       inputs = if (is.null(s$inputs)) NULL else s$inputs[, c("path", "size"), drop = FALSE])
}

cur_files <- list.files(cli$current_dir, pattern = "\\.qs$")
for (f in base_files) {
  cur <- file.path(cli$current_dir, f)
  if (!file.exists(cur)) { expect(FALSE, sprintf("%s: present in current output", f)); next }
  a <- qs_read(file.path(cli$baseline_dir, f)); b <- qs_read(cur)
  ok <- identical(strip(a), strip(b))
  expect(ok, sprintf("%s: payload identical", f))
  if (!ok) print(all.equal(strip(a), strip(b)))
  ok_s <- identical(stamp_core(a), stamp_core(b))
  expect(ok_s, sprintf("%s: provenance stamp agrees (modulo time/R version/mtimes)", f))
  if (!ok_s) print(all.equal(stamp_core(a), stamp_core(b)))
}
extra <- setdiff(cur_files, base_files)
expect(length(extra) == 0L, sprintf("no artefact written that the baseline lacks (extra: %s)",
                                    if (length(extra)) paste(extra, collapse = ", ") else "none"))

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
