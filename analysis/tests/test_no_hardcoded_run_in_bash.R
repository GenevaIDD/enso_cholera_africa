# test_no_hardcoded_run_in_bash.R — HARD: no maintained bash driver names a
#   prod run.
#
# Usage:
#   Rscript analysis/tests/test_no_hardcoded_run_in_bash.R
#
# Details:
#   A production run is named ONLY in analysis/configs/deploy/prod_run<N>.yml;
#   drivers that need one take `--deploy <yaml>` and read its paths from
#   `deploy_prod_run.R --stage resolve`. Hardcoded run defaults must be
#   bumped by hand per rollout and drift from the run the deploy YAML names.
#   Also asserts every deploy-aware driver still resolves through the deploy
#   tool. Pure text -- no fit, no data, no network.
#
# See: analysis/bash/README.md ("Run identity")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# no maintained bash driver hardcodes a prod_run<N>\n")

# What counts: any non-comment line of a real (non-symlink) *.sh under
# analysis/bash/{stages,campaigns,tools,slurm}/ containing `prod_run<digit>` or
# `_run<digit>`. Comments may cite runs (history, examples); `prod_run<N>`
# placeholders do not match. archive/ is exempt (run-pinned by definition).
# Ratchet: KNOWN_RUN_PINNED is empty; any new hit fails, and an entry that
# stops matching must be deleted (a stale entry also fails).
KNOWN_RUN_PINNED <- character(0) # "relative/path:line" entries; empty by design

PATTERN <- "prod_run[0-9]|_run[0-9]"
TIERS <- file.path("analysis/bash", c("stages", "campaigns", "tools", "slurm"))

sh_files <- unlist(lapply(
  TIERS,
  list.files,
  pattern = "\\.sh$",
  full.names = TRUE
))
sh_files <- sh_files[!nzchar(Sys.readlink(sh_files))]
# Tools that pin one run by design, not drivers of a prod run: the synthetic
# pipeline runs its own deploy (configs/deploy/synth_run1.yml).
RUN_PINNED_TOOLS <- c("analysis/bash/tools/run_synthetic_pipeline.sh")
sh_files <- setdiff(sh_files, RUN_PINNED_TOOLS)
expect(
  length(sh_files) > 20L,
  sprintf("scanned the tiered bash drivers (%d files)", length(sh_files))
)

hits <- character(0)
for (f in sh_files) {
  lines <- readLines(f, warn = FALSE)
  code <- sub("(^|[[:space:]])#.*$", "", lines) # drop comments (also trailing ones)
  idx <- grep(PATTERN, code)
  if (length(idx)) hits <- c(hits, sprintf("%s:%d", f, idx))
}

new_hits <- setdiff(hits, KNOWN_RUN_PINNED)
expect(
  length(new_hits) == 0L,
  sprintf(
    "no non-comment run-numbered token (offenders: %s)",
    if (length(new_hits)) paste(new_hits, collapse = ", ") else "none"
  )
)
stale <- setdiff(KNOWN_RUN_PINNED, hits)
expect(
  length(stale) == 0L,
  sprintf(
    "no stale KNOWN_RUN_PINNED entry (stale: %s)",
    if (length(stale)) paste(stale, collapse = ", ") else "none"
  )
)

# Deploy-aware drivers resolve their run through the deploy tool ---------------
mentions <- function(f, s) {
  any(grepl(s, readLines(f, warn = FALSE), fixed = TRUE))
}
DEPLOY_AWARE <- c(
  "analysis/bash/stages/04_fit_kl.sh",
  "analysis/bash/stages/05_loo_compare.sh",
  "analysis/bash/stages/07_forecast_groups.sh",
  "analysis/bash/stages/07_hindcast.sh",
  "analysis/bash/campaigns/enso_grouped_singlelag.sh",
  "analysis/bash/campaigns/dur_enso_grouped_singlelag.sh",
  "analysis/bash/campaigns/mon_enso_grouped_singlelag.sh",
  "analysis/bash/campaigns/postfit.sh",
  "analysis/bash/campaigns/forecast_groups_pipeline.sh",
  "analysis/bash/tools/check_forecast_groups_files.sh"
)
# A code release (analysis/RELEASE) ships neither of the last two.
if (file.exists("analysis/RELEASE")) {
  DEPLOY_AWARE <- DEPLOY_AWARE[file.exists(DEPLOY_AWARE)]
}
for (f in DEPLOY_AWARE) {
  expect(
    file.exists(f) &&
      mentions(f, "--deploy)") &&
      mentions(f, "deploy_prod_run.R -C \"$DEPLOY\" --stage resolve"),
    sprintf(
      "%s takes --deploy and resolves via deploy_prod_run.R --stage resolve",
      f
    )
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
