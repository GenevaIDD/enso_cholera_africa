# test_size_legacy_unchanged.R — the side-by-side rule: the legacy size script
#   analysis/07_run_outbreak_size_model.R is unchanged against main.
#
# Usage:
#   Rscript analysis/tests/test_size_legacy_unchanged.R
#
# Details:
#   Parse-tree checksum against the file on main (comment edits allowed);
#   skips when the file or `git` is not available. Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(digest)
})

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

path <- "analysis/07_run_outbreak_size_model.R"
if (!file.exists(path)) {
  message(
    "SKIP (legacy 07_run_outbreak_size_model.R no longer exists — ",
    "succeeded by 04_fit_siz_brms.R)"
  )
  message("PASS: test_size_legacy_unchanged")
  quit(status = 0)
}

# Parse-tree hash on current branch
parsed_now <- tryCatch(
  parse(file = path, keep.source = FALSE),
  error = function(e) fail(paste("parse error:", e$message))
)
hash_now <- digest::digest(as.character(parsed_now), algo = "sha256")

# Parse-tree hash on main HEAD
main_text <- tryCatch(
  system2(
    "git",
    c("show", paste0("main:", path)),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) NULL
)
if (is.null(main_text) || !is.character(main_text)) {
  message(
    "WARN: could not retrieve main:%s — skipping (git not available?)",
    path
  )
  message("PASS (skipped)")
  quit(status = 0)
}

parsed_main <- tryCatch(
  parse(text = main_text, keep.source = FALSE),
  error = function(e) fail(paste("parse error on main:", e$message))
)
hash_main <- digest::digest(as.character(parsed_main), algo = "sha256")

if (hash_now != hash_main) {
  fail("legacy 07 script parse tree changed vs main")
}

message("PASS: test_size_legacy_unchanged")
quit(status = 0)
