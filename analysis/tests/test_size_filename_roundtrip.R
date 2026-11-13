# test_size_filename_roundtrip.R — deterministic-filename guarantees.
#
# Hard-fail gate. Asserts:
#   1. makeRunFile(opt, "size_fit") is deterministic given identical opts.
#   2. Two configs differing only in defaulted fields produce identical basenames.
#   3. Two configs differing in any non-defaulted field produce distinct basenames.
#   4. Each `size_*` kind produces a unique basename for the same opt.
#

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")   # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/paths.R")   # makeRunFile() -- used directly
source("analysis/utils/provenance.R")   # transitive dependency of paths.R (utils/README.md)

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

# Base opt: load defaults, surface size namespace.
opt_base <- make_options_from_config(option_list = list(),
                                     cli_argv       = character(),
                                     model          = "size")

# 1) determinism
p1 <- makeRunFile(opt_base, kind = "size_fit")
p2 <- makeRunFile(opt_base, kind = "size_fit")
if (!identical(p1, p2)) fail("size_fit filename not deterministic")

# 2) defaulted-field invariance: changing IO outdir alone changes dirname, not basename
opt_alt <- opt_base
opt_alt$outdir <- "generated_data/alt/"
if (basename(makeRunFile(opt_alt, "size_fit")) != basename(p1)) {
  fail("changing only outdir altered the basename")
}

# 3) any non-defaulted change must alter the basename
opt_eth <- opt_base
opt_eth$country <- "ETH"
if (basename(makeRunFile(opt_eth, "size_fit")) == basename(p1)) {
  fail("country change did not alter basename")
}

# 4) each kind unique
kinds <- c("size_fit", "size_re_stats", "size_coef_summary", "size_loo", "size_pp_check")
bn    <- vapply(kinds, function(k) basename(makeRunFile(opt_base, k)), character(1))
if (anyDuplicated(bn)) fail("two size_* kinds collide on basename")

message("PASS: test_size_filename_roundtrip")
quit(status = 0)
