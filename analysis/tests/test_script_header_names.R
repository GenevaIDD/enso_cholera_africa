# test_script_header_names.R — HARD: numbered scripts describe themselves under
#   their own name.
#
# Usage:
#   Rscript analysis/tests/test_script_header_names.R
#
# Details:
#   Pure text, no fit, no data, no network. A rename can leave a header behind,
#   naming a file that no longer exists in the header or the usage lines; two
#   cheap regex checks catch that (the first comment line names the file;
#   every `Rscript analysis/NN*.R` command resolves).

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# numbered-script headers name their own file; Rscript commands resolve\n")

numbered <- list.files(
  "analysis",
  pattern = "^[0-9][0-9].*\\.R$",
  full.names = TRUE
)
expect(
  length(numbered) > 40L,
  sprintf("scanned %d numbered scripts", length(numbered))
)

# 1. the first comment line of every analysis/[0-9][0-9]*.R -- skipping a
#    shebang and `# ====` / `# ----` rule lines -- starts with
#    `# <own filename>` (or `# analysis/<own filename>`).
# 1. first comment line = own filename ------------------------------------------
first_comment <- function(f) {
  lines <- readLines(f, warn = FALSE, n = 50L)
  lines <- lines[
    grepl("^#", lines) &
      !grepl("^#!", lines) &
      !grepl("^#\\s*[=-]{5,}\\s*$", lines)
  ]
  if (length(lines)) lines[[1]] else ""
}
bad <- character(0)
for (f in numbered) {
  own <- basename(f)
  fc <- first_comment(f)
  ok <- startsWith(fc, paste0("# ", own)) ||
    startsWith(fc, paste0("# analysis/", own))
  if (!ok) bad <- c(bad, sprintf("%s (%s)", own, fc))
}
expect(
  length(bad) == 0L,
  sprintf(
    "first header line names the script itself (offenders: %s)",
    if (length(bad)) paste(bad, collapse = "; ") else "none"
  )
)

# 2. every `Rscript analysis/NN*.R` command written in a numbered script, a
#    utils/ module or a maintained bash driver (stages/ campaigns/ tools/
#    slurm/; archive/ is exempt) names a file that exists.
# 2. Rscript commands resolve ---------------------------------------------------
sh_files <- list.files(
  file.path("analysis/bash", c("stages", "campaigns", "tools", "slurm")),
  pattern = "\\.sh$",
  full.names = TRUE
)
scan_files <- c(
  numbered,
  list.files("analysis/utils", pattern = "\\.R$", full.names = TRUE),
  sh_files
)
dangling <- character(0)
for (f in scan_files) {
  lines <- readLines(f, warn = FALSE)
  hits <- regmatches(
    lines,
    gregexpr("Rscript +analysis/[0-9]{2}[a-z]?_[A-Za-z0-9_]+\\.R", lines)
  )
  for (k in which(lengths(hits) > 0L)) {
    for (h in hits[[k]]) {
      target <- sub("^Rscript +", "", h)
      if (!file.exists(target)) {
        dangling <- c(dangling, sprintf("%s:%d -> %s", f, k, target))
      }
    }
  }
}
expect(
  length(dangling) == 0L,
  sprintf(
    "every `Rscript analysis/NN*.R` names an existing script (dangling: %s)",
    if (length(dangling)) paste(dangling, collapse = "; ") else "none"
  )
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
