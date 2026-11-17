# test_monthly_config_merge.R — .fill_model_block() sibling-inheritance
#   precedence on the monthly track.
#
# Usage:
#   Rscript analysis/tests/test_monthly_config_merge.R
#
# Details:
#   The monthly counterpart of the occurrence / duration config-merge tests.
#   Guarded rule: a field set only under the plain-backend sibling (e.g.
#   monthly.brms.re_uncorrelated) must still resolve on the active `_cb`
#   block (monthly$brms_cb$re_uncorrelated), so make_batch_configs.R's
#   strip_backend_siblings() can safely drop the plain `brms:` sibling once
#   its value is merged forward. Precedence: block[[be]] (explicit) >
#   block[[plain]] (sibling) > defaults[[be]] (canonical).

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

defaults <- get_default_options()
stopifnot(
  !is.null(defaults$monthly),
  identical(defaults$monthly$backend, "brms_cb")
)

resolve_yaml <- function(yaml_body) {
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(yaml_body, tmp)
  make_options_from_config(
    option_list = list(),
    cli_argv = c("--config", tmp),
    model = "monthly"
  )
}

# --- 1. Sanity: monthly block resolves from an empty (backend-only) YAML ---
cat("\n[1] empty YAML (backend only) -> monthly defaults populated\n")

opt <- resolve_yaml("monthly:\n  backend: brms_cb\n")
expect(!is.null(opt$monthly), "opt$monthly is present")
expect(identical(opt$monthly$backend, "brms_cb"), "backend matches YAML value")
expect(
  !is.null(opt$monthly$brms_cb$sampling),
  "active backend (brms_cb) sub-block is filled from defaults"
)

# --- 2. .fill_model_block(): brms_cb inherits re_uncorrelated from the -----
#        plain brms sibling ------------------------------------------------
cat(
  "\n[2] .fill_model_block(): monthly/brms_cb inherits re_uncorrelated from the plain brms sibling\n"
)

opt2a <- resolve_yaml(paste(
  "monthly:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  sep = "\n"
))
expect(
  isTRUE(opt2a$monthly$brms_cb$re_uncorrelated),
  "value set only on the plain brms sibling reaches monthly$brms_cb$re_uncorrelated"
)

opt2b <- resolve_yaml(paste(
  "monthly:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  "  brms_cb:",
  "    re_uncorrelated: false",
  sep = "\n"
))
expect(
  identical(opt2b$monthly$brms_cb$re_uncorrelated, FALSE),
  "explicit brms_cb$re_uncorrelated wins over a conflicting plain-sibling value"
)

# --- 3. Plain-backend monthly config is unaffected by the _cb-only branch --
cat("\n[3] backend: brms -> no sibling-inheritance branch triggered\n")

opt3 <- resolve_yaml(paste(
  "monthly:",
  "  backend: brms",
  "  brms:",
  "    re_uncorrelated: true",
  sep = "\n"
))
expect(
  isTRUE(opt3$monthly$brms$re_uncorrelated),
  "plain brms backend: re_uncorrelated set directly on brms is preserved"
)
expect(
  is.null(opt3$monthly$brms_cb),
  "plain brms backend: no brms_cb block materialised"
)

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  quit(save = "no", status = 1)
}
message("PASS  test_monthly_config_merge")
