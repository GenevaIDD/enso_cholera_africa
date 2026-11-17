# test_offset_wellformed.R — HARD: the Poisson offset of the saved survival
#   dataset is well-defined on every row.
#
# Usage:
#   Rscript analysis/tests/test_offset_wellformed.R
#   Rscript analysis/tests/test_offset_wellformed.R --baseline 1990_2019
#
# Details:
#   Asserts min(t_end_seas - t_start_seas) > 0 and that the computed
#   `exposure` is integer-valued and >= 1. Run before any Poisson fit; exits
#   non-zero on failure. Writes generated_data/test_reports/
#   offset_wellformed.md.

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

source("analysis/utils/base.R") # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/paths.R") # make_surv_dat_file() -- used directly
source("analysis/utils/provenance.R") # transitive dependency of paths.R (utils/README.md)

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("-b", "--baseline"),
      default = "1990_2019",
      type = "character",
      help = "baseline period (matches 03_prepare_occ.R)"
    )
  )
))

report_dir <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir, "offset_wellformed.md")

con <- file(report, "w")
on.exit(close(con))
cat("# test_offset_wellformed\n\n", file = con)
cat(sprintf("Baseline: `%s`\n\n", opt$baseline), file = con)

surv_file <- make_surv_dat_file(opt = opt)
if (!file.exists(surv_file)) {
  cat(
    "FAIL: survival data file not found at ",
    surv_file,
    "\n",
    "Run analysis/03_prepare_occ.R first.\n",
    sep = "",
    file = con
  )
  # A code-release checkout (analysis/RELEASE) has no outbreak data: nothing to
  # check there; the synthetic world's panel is checked by test_synthetic_world.R.
  if (file.exists("analysis/RELEASE")) {
    message(
      "test_offset_wellformed: SKIP (release checkout, no survival data: ",
      surv_file,
      ")"
    )
    quit(status = 0)
  }
  message("test_offset_wellformed: FAIL (missing survival data)")
  quit(status = 1)
}

surv_dat_all <- readRDS(surv_file)

# Match the filtering applied in 04_run_survival_models.R before fitting:
# zero-length intervals at ADM0 rows are dropped there too.
dat <- surv_dat_all %>%
  filter(admin_level != "ADM0", !is.na(t_start), !is.na(t_end))

interval_days <- as.integer(dat$t_end_seas) - as.integer(dat$t_start_seas)
exposure <- pmax(interval_days, 1L)

cat("## Diagnostics\n\n", file = con)
cat(sprintf("- nrow (post-filter): %d\n", nrow(dat)), file = con)
cat(
  sprintf("- min(t_end_seas - t_start_seas) = %d\n", min(interval_days)),
  file = con
)
cat(
  sprintf("- max(t_end_seas - t_start_seas) = %d\n", max(interval_days)),
  file = con
)
cat(
  sprintf("- # rows with interval <= 0:       %d\n", sum(interval_days <= 0)),
  file = con
)
cat(
  sprintf(
    "- # rows with non-integer interval: %d\n",
    sum(interval_days != round(interval_days))
  ),
  file = con
)
cat("\n", file = con)

failures <- c()
if (min(interval_days) <= 0) {
  failures <- c(
    failures,
    sprintf(
      "min interval is %d (need > 0); use exposure = pmax(..., 1L) or drop zero-length rows",
      min(interval_days)
    )
  )
}
if (any(exposure < 1L)) {
  failures <- c(failures, "exposure has values < 1 after pmax")
}
if (!is.integer(exposure)) {
  failures <- c(failures, "exposure is not integer-valued")
}

cat("## Verdict\n\n", file = con)
if (length(failures) == 0) {
  cat("PASS\n", file = con)
  message("test_offset_wellformed: PASS")
  quit(status = 0)
} else {
  cat("FAIL:\n", file = con)
  for (f in failures) {
    cat("- ", f, "\n", sep = "", file = con)
  }
  message("test_offset_wellformed: FAIL")
  for (f in failures) {
    message("  - ", f)
  }
  quit(status = 1)
}
