# test_no_scratchpad_callers.R — HARD: nothing on the production path runs a
#   scratchpad.
#
# Usage:
#   Rscript analysis/tests/test_no_scratchpad_callers.R
#
# Details:
#   "A number means it runs; no number means it's sourced" breaks when a
#   numbered step, a sourced library, or an HPC driver in analysis/bash/
#   executes analysis/scratchpads/*: a reader can no longer tell which
#   scratchpads are load-bearing. Load-bearing code lives in numbered
#   scripts; this test keeps it that way. Pure text/parse -- no fit, no
#   data, no network.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_parse.R") # parse_data_of()

cat("# no production-path file invokes analysis/scratchpads/\n")

# What counts as a caller:
#   * R files (analysis/[0-9][0-9]*.R and analysis/utils/*.R): any STRING
#     CONSTANT that contains "scratchpads/" -- a source() path, an
#     Rscript/system() command, a file path. Comments are COMMENT tokens in
#     R's parse data, so a comment naming scratchpads/ stays legal.
#   * Shell files (analysis/bash/*.sh): any non-comment line containing
#     "scratchpads/".
# Ratchet: KNOWN_SCRATCHPAD_CALLERS is empty; any new caller fails, and
# an entry that stops calling must be deleted (a stale entry also fails).
KNOWN_SCRATCHPAD_CALLERS <- character(0) # "relative/path:kind" entries; empty by design

PATTERN <- "scratchpads/"

r_files <- c(
  list.files("analysis", pattern = "^[0-9][0-9].*\\.R$", full.names = TRUE),
  list.files("analysis/utils", pattern = "\\.R$", full.names = TRUE)
)
# Recursive: bash/ is tiered (stages/ campaigns/ slurm/ tools/); the flat
# top-level *.sh are back-compat symlinks, skipped so each script is scanned once.
sh_files <- list.files(
  "analysis/bash",
  pattern = "\\.sh$",
  full.names = TRUE,
  recursive = TRUE
)
sh_files <- sh_files[!nzchar(Sys.readlink(sh_files))]
expect(
  length(r_files) > 20L && length(sh_files) > 0L,
  "scanned the numbered scripts, utils and bash drivers"
)

callers <- character(0)

for (f in r_files) {
  pd <- parse_data_of(f)
  if (is.null(pd)) {
    callers <- c(callers, sprintf("%s:unparseable", f))
    next
  }
  strs <- pd$text[pd$token == "STR_CONST"]
  if (any(grepl(PATTERN, strs, fixed = TRUE))) {
    callers <- c(callers, sprintf("%s:R", f))
  }
}

for (f in sh_files) {
  lines <- readLines(f, warn = FALSE)
  code <- sub("(^|[[:space:]])#.*$", "", lines) # drop comments (also trailing ones)
  if (any(grepl(PATTERN, code, fixed = TRUE))) {
    callers <- c(callers, sprintf("%s:sh", f))
  }
}

new_callers <- setdiff(callers, KNOWN_SCRATCHPAD_CALLERS)
expect(
  length(new_callers) == 0L,
  sprintf(
    "no file invokes a scratchpad (offenders: %s)",
    if (length(new_callers)) paste(new_callers, collapse = ", ") else "none"
  )
)
stale <- setdiff(KNOWN_SCRATCHPAD_CALLERS, callers)
expect(
  length(stale) == 0L,
  sprintf(
    "no stale KNOWN_SCRATCHPAD_CALLERS entry (stale: %s)",
    if (length(stale)) paste(stale, collapse = ", ") else "none"
  )
)

# The promoted scripts exist and parse ------------------------------------------------
promoted <- c(
  "analysis/06_fit_figures_cb_grouped.R",
  "analysis/06_fit_figures_cb_grouped_composite.R",
  "analysis/07d_simulator_figures_counterfactual.R",
  "analysis/01_preprocess_figures_data_overview.R",
  "analysis/05_fit_postprocess_kl_sweep.R",
  "analysis/06_fit_figures_kl_sweep.R"
)
for (f in promoted) {
  ok <- file.exists(f) && !inherits(try(parse(f), silent = TRUE), "try-error")
  expect(ok, sprintf("promoted script exists and parses: %s", f))
}
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
