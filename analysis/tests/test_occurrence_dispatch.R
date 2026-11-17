# test_occurrence_dispatch.R — HARD: static check that 04_fit_run.R routes
#   correctly.
#
# Usage:
#   Rscript analysis/tests/test_occurrence_dispatch.R
#
# Details:
#   04_fit_run.R reads the track block's `backend:` and dispatches
#   *programmatically*: the target path is built with sprintf() from the
#   track token and the backend string (e.g. occurrence + brms_cb ->
#   analysis/04_fit_occ_brms_cb.R) and validated with file.exists(), not
#   looked up in a static table. "Add a backend" is "add a file", but there
#   are no literal target strings to assert on, so the test asserts the
#   dispatch *mechanism*: the backend keys read, sprintf + file.exists, and
#   a stop() when neither block is present.

suppressPackageStartupMessages({})

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

path <- "analysis/04_fit_run.R"
if (!file.exists(path)) {
  fail(paste("missing dispatcher:", path))
}

parsed <- tryCatch(parse(file = path, keep.source = TRUE), error = function(e) {
  fail(paste("parse error:", e$message))
})
src <- paste(as.character(parsed), collapse = "\n")

# --- 1. reads the correct nested backend keys --------------------------------
if (!grepl("cfg\\$occurrence\\$backend", src, perl = TRUE)) {
  fail("04_fit_run.R does not read cfg$occurrence$backend")
}
if (!grepl("cfg\\$size\\$backend", src, perl = TRUE)) {
  fail("04_fit_run.R does not read cfg$size$backend")
}

# --- 2. dispatch is programmatic (sprintf + file.exists), not a static table -
if (!grepl("sprintf\\(\\s*\"analysis/04_fit_%s_%s\\.R\"", src, perl = TRUE)) {
  fail(
    "04_fit_run.R: dispatch is not built via sprintf(\"analysis/04_fit_%s_%s.R\", ...)"
  )
}
if (!grepl("file.exists(script)", src, fixed = TRUE)) {
  fail("04_fit_run.R: dispatch target is not validated via file.exists(script)")
}
# The track map (occurrence -> occ, size -> siz) that feeds the sprintf() call.
if (
  !grepl("occurrence\\s*=\\s*\"occ\"", src, perl = TRUE) ||
    !grepl("size\\s*=\\s*\"siz\"", src, perl = TRUE)
) {
  fail("04_fit_run.R: missing occurrence/size -> occ/siz track map")
}

# --- 3. stop() when neither block present ------------------------------------
if (!grepl("stop\\(", src)) {
  fail(
    "04_fit_run.R: no stop() call — missing guard for configs with neither occurrence: nor size:"
  )
}

message("PASS  test_occurrence_dispatch")
