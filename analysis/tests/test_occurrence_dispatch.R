# test_occurrence_dispatch.R — static check that 04_fit_run.R routes correctly.
#
# Hard-fail gate. 04_fit_run.R detects whether the config carries an
# `occurrence:` or `size:` block, reads the corresponding `backend:`, and
# dispatches *programmatically* — the target script path is built with
# sprintf() from the track ("occ"/"siz") and the literal backend string
# (e.g. model=occurrence, backend=brms_cb -> analysis/04_fit_occ_brms_cb.R)
# and validated with file.exists() rather than looked up in a static table.
# This means "add a backend" is "add a file", but it also means we can no
# longer assert on literal target-script strings in the source (there are
# none) — instead we assert on the dispatch *mechanism*.
#
# Checks:
#   1. dispatcher reads cfg$occurrence$backend and cfg$size$backend
#   2. dispatch is built via sprintf() from a track token ("occ"/"siz") and
#      the backend string, gated by file.exists() (no hardcoded DISPATCH
#      table to go stale)
#   3. the script stop()s when neither block is present

suppressPackageStartupMessages({})

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

path <- "analysis/04_fit_run.R"
if (!file.exists(path)) fail(paste("missing dispatcher:", path))

parsed <- tryCatch(parse(file = path, keep.source = TRUE),
                   error = function(e) fail(paste("parse error:", e$message)))
src <- paste(as.character(parsed), collapse = "\n")

# --- 1. reads the correct nested backend keys --------------------------------
if (!grepl("cfg\\$occurrence\\$backend", src, perl = TRUE))
  fail("04_fit_run.R does not read cfg$occurrence$backend")
if (!grepl("cfg\\$size\\$backend", src, perl = TRUE))
  fail("04_fit_run.R does not read cfg$size$backend")

# --- 2. dispatch is programmatic (sprintf + file.exists), not a static table -
if (!grepl("sprintf\\(\\s*\"analysis/04_fit_%s_%s\\.R\"", src, perl = TRUE))
  fail("04_fit_run.R: dispatch is not built via sprintf(\"analysis/04_fit_%s_%s.R\", ...)")
if (!grepl("file.exists(script)", src, fixed = TRUE))
  fail("04_fit_run.R: dispatch target is not validated via file.exists(script)")
# The track map (occurrence -> occ, size -> siz) that feeds the sprintf() call.
if (!grepl("occurrence\\s*=\\s*\"occ\"", src, perl = TRUE) ||
    !grepl("size\\s*=\\s*\"siz\"", src, perl = TRUE))
  fail("04_fit_run.R: missing occurrence/size -> occ/siz track map")

# --- 3. stop() when neither block present ------------------------------------
if (!grepl("stop\\(", src))
  fail("04_fit_run.R: no stop() call — missing guard for configs with neither occurrence: nor size:")

message("PASS  test_occurrence_dispatch")
