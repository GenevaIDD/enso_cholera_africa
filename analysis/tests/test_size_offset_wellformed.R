# test_size_offset_wellformed.R — the size fit's log offsets are finite.
#
# Usage:
#   Rscript analysis/tests/test_size_offset_wellformed.R --baseline 1990_2019
#
# Details:
#   After the size runner's filter (admin_level != "ADM0", non-NA TL/TR),
#   log(pop) and log(mean_rate) are finite for every row that would be fed to
#   brm(). Skips when the prep output is missing. Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("-b", "--baseline"),
      default = "1990_2019",
      type = "character"
    )
  )
))

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

path <- "generated_data/outbreak_size_analysis_data.rds"
if (!file.exists(path)) {
  message("WARN: prep output missing — skipping")
  message("PASS (skipped)")
  quit(status = 0)
}

d <- readRDS(path) %>%
  filter(admin_level != "ADM0", !is.na(TL), !is.na(TR)) %>%
  ungroup() %>%
  mutate(mean_rate = mean(suspected_cases / pop, na.rm = TRUE))

if (any(!is.finite(log(d$pop)))) {
  fail("log(pop) not finite for some rows fed to fit")
}

if (!is.finite(unique(d$mean_rate))) {
  fail("mean_rate is not finite")
}

if (any(!is.finite(log(d$mean_rate)))) {
  fail("log(mean_rate) not finite")
}

message("PASS: test_size_offset_wellformed")
quit(status = 0)
