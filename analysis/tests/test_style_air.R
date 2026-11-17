# test_style_air.R — HARD: the R code is formatted as air formats it (air.toml).
#
# Usage:
#   Rscript analysis/tests/test_style_air.R
#
# Details:
#   Runs `air format --check` over the repo (air.toml: line width 80,
#   scratchpads excluded) and fails on any file air would change; fix with
#   `bash analysis/bash/tools/style.sh --fix`. SKIPs (passes) when air is not
#   installed, so a machine without the formatter can still run the suite.
#
# See: analysis/bash/tools/style.sh

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# R code formatted by air\n")

air <- Sys.which("air")
if (!nzchar(air)) {
  cat("  SKIP  air is not on PATH (https://posit-dev.github.io/air/)\n")
  quit(status = 0L)
}
ver <- suppressWarnings(system2(air, "--version", stdout = TRUE, stderr = TRUE))
cat(sprintf("  %s\n", ver[1]))

expect(file.exists("air.toml"), "air.toml at the repo root")
out <- suppressWarnings(system2(
  air,
  c("format", "--check", "."),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(out, "status") %||% 0L
offenders <- trimws(sub("^\\s*-\\s*", "", grep("^\\s*- ", out, value = TRUE)))
expect(
  identical(as.integer(status), 0L),
  sprintf(
    "air format --check: every file formatted (%d to reformat%s)",
    length(offenders),
    if (length(offenders)) {
      paste0(": ", paste(head(offenders, 10), collapse = ", "))
    } else {
      ""
    }
  )
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
quit(status = as.integer(FAIL > 0L))
