# test_diagnostics.R — generic MCMC diagnostics of a saved fit (brmsfit or
#   brms_kl kl_fit).
#
# Usage:
#   Rscript analysis/tests/test_diagnostics.R --fit path/to/fit.rds
#
# Details:
#   Goes through the shared utils/diagnostics.R::fit_diagnostics(). Soft:
#   divergences in sparse-data fits are expected, so it reports but does not
#   exit non-zero. Thresholds (DIAG_THRESHOLDS): Rhat < 1.01, bulk-ESS > 400,
#   tail-ESS > 400, zero divergences post-warmup, no iteration at the maximum
#   tree depth.

suppressPackageStartupMessages({
  library(optparse)
})
source("analysis/utils/base.R")
source("analysis/utils/diagnostics.R")

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("--fit"),
      default = NA_character_,
      type = "character",
      help = "path to saved fit (.rds)"
    )
  )
))
if (is.na(opt$fit)) {
  stop("--fit is required")
}

report_dir <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
tag <- tools::file_path_sans_ext(basename(opt$fit))
report <- file.path(report_dir, sprintf("diagnostics_%s.md", tag))
con <- file(report, "w")
on.exit(close(con), add = TRUE)
cat(sprintf("# test_diagnostics — %s\n\n", tag), file = con)

d <- fit_diagnostics(readRDS(opt$fit))
th <- d$thresholds

cat(sprintf("- n_parameters: %d\n", d$n_params), file = con)
cat(sprintf("- max Rhat:   %.3f\n", d$max_rhat), file = con)
cat(sprintf("- min bulk-ESS: %.0f\n", d$min_ess_bulk), file = con)
cat(sprintf("- min tail-ESS: %.0f\n", d$min_ess_tail), file = con)
cat(
  sprintf("- params with Rhat > %.2f: %d\n", th[["rhat"]], d$n_rhat_gt),
  file = con
)
cat(
  sprintf("- params with bulk-ESS < %.0f: %d\n", th[["ess"]], d$n_ess_bulk_lt),
  file = con
)
cat(
  sprintf("- params with tail-ESS < %.0f: %d\n", th[["ess"]], d$n_ess_tail_lt),
  file = con
)
cat(sprintf("- divergent transitions: %s\n", d$divergences), file = con)
cat(sprintf("- max-treedepth hits:    %s\n", d$treedepth_hits), file = con)
cat(sprintf("- min E-BFMI:            %.3f\n", d$min_ebfmi), file = con)

cat("\n## Verdict (soft)\n\n", file = con)
if (
  d$n_rhat_gt == 0 &&
    d$n_ess_bulk_lt == 0 &&
    d$n_ess_tail_lt == 0 &&
    isTRUE(d$divergences %in% c(0L, NA_integer_))
) {
  cat("PASS\n", file = con)
  message("test_diagnostics: PASS (", tag, ")")
} else {
  cat("SOFT-FAIL (see counts above)\n", file = con)
  message("test_diagnostics: SOFT-FAIL (", tag, ")")
}
quit(status = 0)
