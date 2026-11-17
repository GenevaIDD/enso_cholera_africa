# test_size_filename_roundtrip.R — deterministic size-track filenames.
#
# Usage:
#   Rscript analysis/tests/test_size_filename_roundtrip.R
#
# Details:
#   make_run_file(opt, "size_fit") is deterministic; configs differing only in
#   defaulted fields share a basename, any non-defaulted change alters it;
#   each `size_*` kind has a unique basename for the same opt. Not registered
#   in analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R") # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/paths.R") # make_run_file() -- used directly
source("analysis/utils/provenance.R") # transitive dependency of paths.R (utils/README.md)

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

# Base opt: load defaults, surface size namespace.
opt_base <- make_options_from_config(
  option_list = list(),
  cli_argv = character(),
  model = "size"
)

# 1) determinism
p1 <- make_run_file(opt_base, kind = "size_fit")
p2 <- make_run_file(opt_base, kind = "size_fit")
if (!identical(p1, p2)) {
  fail("size_fit filename not deterministic")
}

# 2) defaulted-field invariance: changing IO outdir alone changes dirname, not basename
opt_alt <- opt_base
opt_alt$outdir <- "generated_data/alt/"
if (basename(make_run_file(opt_alt, "size_fit")) != basename(p1)) {
  fail("changing only outdir altered the basename")
}

# 3) any non-defaulted change must alter the basename
opt_eth <- opt_base
opt_eth$country <- "ETH"
if (basename(make_run_file(opt_eth, "size_fit")) == basename(p1)) {
  fail("country change did not alter basename")
}

# 4) each kind unique
kinds <- c(
  "size_fit",
  "size_re_stats",
  "size_coef_summary",
  "size_loo",
  "size_pp_check"
)
bn <- vapply(
  kinds,
  function(k) basename(make_run_file(opt_base, k)),
  character(1)
)
if (anyDuplicated(bn)) {
  fail("two size_* kinds collide on basename")
}

message("PASS: test_size_filename_roundtrip")
quit(status = 0)
