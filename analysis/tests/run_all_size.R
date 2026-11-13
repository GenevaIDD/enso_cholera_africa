# run_all_size.R — outbreak-size pipeline validation harness.
#
# Mirror of analysis/tests/run_all.R for the size pipeline.
#
# CLI:
#   Rscript analysis/tests/run_all_size.R [--quick] [--baseline 1990_2019]

suppressPackageStartupMessages({
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-b", "--baseline"), default = "1990_2019", type = "character"),
  make_option(c("--quick"), default = FALSE, type = "logical",
              help = "skip slow tests (diagnostics)")
)))

report_dir   <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(report_dir, "SIZE_SUMMARY.md")

run_one <- function(script, args = character(), hard_fail = FALSE) {
  cmd <- paste(c("Rscript", script, args), collapse = " ")
  message("\n>>> ", cmd)
  status <- system(cmd)
  list(script = script, args = paste(args, collapse = " "),
       status = status, hard_fail = hard_fail)
}

results <- list()

# Hard-fail gates --------------------------------------------------------
results$response_invariant <- run_one(
  "analysis/tests/test_size_response_invariant.R",
  hard_fail = TRUE
)

results$offset_wellformed <- run_one(
  "analysis/tests/test_size_offset_wellformed.R",
  args = sprintf("--baseline %s", opt$baseline),
  hard_fail = TRUE
)

results$filename_roundtrip <- run_one(
  "analysis/tests/test_size_filename_roundtrip.R",
  hard_fail = TRUE
)

results$legacy_unchanged <- run_one(
  "analysis/tests/test_size_legacy_unchanged.R",
  hard_fail = TRUE
)

results$postprocessing_adapter <- run_one(
  "analysis/tests/test_size_postprocessing_adapter.R",
  args = sprintf("--baseline %s --country ETH", opt$baseline),
  hard_fail = TRUE
)

results$nested_slopes <- run_one(
  "analysis/tests/test_size_nested_slopes.R",
  hard_fail = TRUE
)

results$backend_dispatch <- run_one(
  "analysis/tests/test_size_backend_dispatch.R",
  hard_fail = TRUE
)

results$family_dispatch <- run_one(
  "analysis/tests/test_size_family_dispatch.R",
  hard_fail = TRUE
)

results$priors <- run_one(
  "analysis/tests/test_size_priors.R",
  hard_fail = TRUE
)

results$config_deep_merge <- run_one(
  "analysis/tests/test_size_config_deep_merge.R",
  hard_fail = TRUE
)

results$batch_config_determinism <- run_one(
  "analysis/tests/test_size_batch_config_determinism.R",
  hard_fail = TRUE
)

results$seasonal_term <- run_one(
  "analysis/tests/test_size_seasonal_term.R",
  hard_fail = TRUE
)

results$simulator_dlnm <- run_one(
  "analysis/tests/test_simulator_dlnm.R",
  hard_fail = TRUE
)

# Soft-fail tests --------------------------------------------------------
if (!opt$quick) {
  results$diagnostics <- run_one(
    "analysis/tests/test_size_diagnostics.R",
    hard_fail = FALSE
  )
}

# Summary ----------------------------------------------------------------

con <- file(summary_path, "w"); on.exit(close(con), add = TRUE)
cat("# Size test harness summary\n\n", file = con)
cat(sprintf("Generated: %s\n\n", Sys.time()), file = con)
cat("| test | status | hard_fail | result |\n", file = con)
cat("|------|--------|-----------|--------|\n", file = con)
hard_failed <- FALSE
for (nm in names(results)) {
  r <- results[[nm]]
  res <- if (r$status == 0) "PASS" else "FAIL"
  if (r$status != 0 && r$hard_fail) hard_failed <- TRUE
  cat(sprintf("| %s | %d | %s | %s |\n",
              nm, r$status, r$hard_fail, res), file = con)
}

if (hard_failed) {
  message("\nrun_all_size: HARD FAILURE — see ", summary_path)
  quit(status = 1)
} else {
  message("\nrun_all_size: all hard gates PASSED — see ", summary_path)
  quit(status = 0)
}
