# test_size_config_deep_merge.R — verify deep-merge semantics of the
# size-namespace config layer.
#
# Hard-fail gate. make_options_from_config() promises that nested blocks —
# specifically opt$size, opt$size$random_effects,
# opt$size$covariates, opt$covariates and opt$enso — are
# *recursively* merged with config_defaults.yml, so that a partial YAML
# override only changes the keys it sets and leaves siblings alone.
#
# Failure modes this guards against:
#   - a refactor of the merge loop that replaces nested blocks wholesale
#     (loses every default key the user didn't restate);
#   - a key-name typo or path mismatch in nested_blocks that drops one of
#     the blocks from the recursive branch;
#   - regression where a single TRUE flag inside size$random_effects
#     silently zeroes the rest of the random_effects map.
#
#
# All checks are pure config-resolution — no fit, no data needed.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

defaults <- get_default_options()
# Sanity: defaults carry a non-empty size block, otherwise the rest of the
# test is meaningless.
stopifnot(!is.null(defaults$size),
          is.list(defaults$size$random_effects),
          is.list(defaults$size$covariates) ||
            is.null(defaults$size$covariates))

# Helper: write a temp YAML and resolve it through the public entry point.
resolve_yaml <- function(yaml_body, model = "size") {
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(yaml_body, tmp)
  make_options_from_config(
    option_list = list(),
    cli_argv    = c("--config", tmp),
    model       = model
  )
}

# --- 1. Empty config: every default surfaces unchanged --------------------
cat("\n[1] empty YAML config -> all defaults populated\n")

opt <- resolve_yaml("size:\n  backend: brms\n")   # minimum legal size config
# .fill_model_block fills only the active backend sub-block, not all backends.
# Check individual keys rather than an exact whole-block comparison.
expect(identical(opt$size$backend,        defaults$size$backend),        "opt$size$backend matches default")
expect(identical(opt$size$family,         defaults$size$family),         "opt$size$family matches default")
expect(identical(opt$size$link,           defaults$size$link),           "opt$size$link matches default")
expect(identical(opt$size$random_effects, defaults$size$random_effects), "opt$size$random_effects matches default (deep merge)")
expect(!is.null(opt$size$brms$sampling),
       "active backend (brms) sub-block is filled from defaults")
expect(identical(opt$enso, defaults$enso),
       "opt$enso matches defaults$enso")

# --- 2. Partial size block: only the touched key changes ------------------
cat("\n[2] partial size: override only family -> siblings preserved\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "size:",
  "  family: poisson",
  sep = "\n"))
expect(identical(opt$size$family, "poisson"),
       "user-set family wins")
expect(identical(opt$size$link, defaults$size$link),
       "size$link default is preserved (not wiped)")
expect(identical(opt$size$admin_level, defaults$size$admin_level),
       "size$admin_level default is preserved")
expect(identical(opt$size$offset_pop, defaults$size$offset_pop),
       "size$offset_pop default is preserved")
expect(identical(opt$size$random_effects, defaults$size$random_effects),
       "size$random_effects map is preserved (recursive merge respects siblings)")

# --- 3. Partial random_effects map: 2-deep merge fills sibling defaults ---
# make_options_from_config() now merges size$random_effects against
# defaults$size$random_effects (two levels deep), so a partial user-provided
# block only changes the keys it sets and leaves sibling defaults in place.
cat("\n[3] partial size$random_effects: 2-deep merge fills sibling defaults\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "size:",
  "  random_effects:",
  "    re_adm1_mean_precip_slope: true",
  sep = "\n"))
expect(isTRUE(opt$size$random_effects$re_adm1_mean_precip_slope),
       "user-set re_adm1_mean_precip_slope is TRUE")
# 2-deep merge: sibling defaults ARE auto-filled when partial random_effects given.
expect(isTRUE(opt$size$random_effects$re_region_country_adm1_intercept),
       "2-deep merge: sibling re_region_country_adm1_intercept default preserved")

# --- 4. Partial size$covariates: only listed keys change ------------------
cat("\n[4] partial size$covariates\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "size:",
  "  covariates:",
  "    mean_precip: true",
  sep = "\n"))
expect(isTRUE(opt$size$covariates$mean_precip),
       "user-set size$covariates$mean_precip is TRUE")

# --- 5. Top-level covariates partial: untouched flags inherit defaults ---
# Pick a key that defaults to FALSE so flipping it to TRUE doesn't trip the
# RE/covariate consistency validator (which fires on mean_precip pairs).
cat("\n[5] partial top-level covariates\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "covariates:",
  "  mean_precip_lag1: true",
  sep = "\n"))
expect(identical(opt$covariates$mean_precip_lag1, TRUE),
       "user-set covariates$mean_precip_lag1 is TRUE")
# Any other default-covariate flag should survive.
other_cov <- setdiff(names(defaults$covariates), "mean_precip_lag1")
if (length(other_cov) > 0) {
  k <- other_cov[1]
  expect(identical(opt$covariates[[k]], defaults$covariates[[k]]),
         sprintf("covariates$%s default preserved alongside the partial override", k))
}

# --- 6. enso partial block and covariates$enso toggle --------------------
cat("\n[6] partial enso block + covariates.enso\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "covariates:",
  "  enso: true",
  sep = "\n"))
expect(isTRUE(opt$covariates$enso),
       "user-set covariates$enso is TRUE")
if (!is.null(defaults$enso$lag_months)) {
  expect(identical(opt$enso$lag_months, defaults$enso$lag_months),
         "enso$lag_months default preserved (enso block is data-prep params only)")
}
expect(is.null(defaults$enso$use),
       "defaults$enso no longer carries a 'use' key")

# --- 7. CLI -> YAML -> defaults priority chain ----------------------------
# Documented contract per make_options_from_config() header:
#   CLI explicit flag > YAML > default
# CURRENT BEHAVIOR (see comment in utils/config.R around the merge loop): the
# CLI-vs-YAML promotion uses a "key present in CLI but absent from YAML"
# heuristic, so optparse-default and optparse-explicit CLI values are
# indistinguishable. As a result, an EXPLICIT CLI flag does NOT currently
# beat a YAML-set value for the same key. This test pins that contract;
# flip the assertion once the heuristic is replaced by an explicit-vs-
# default detector.
cat("\n[7] CLI vs YAML priority (current heuristic)\n")

tmp <- tempfile(fileext = ".yml")
writeLines(c("backend: brms", "country: ETH"), tmp)
opt <- make_options_from_config(
  option_list = list(
    optparse::make_option("--country", default = NULL, type = "character")
  ),
  cli_argv = c("--config", tmp, "--country", "BEN"),
  model = "size"
)
unlink(tmp)
# Current behaviour: YAML wins because the merge heuristic only promotes
# CLI keys that are absent from the YAML (cf. utils/config.R, the "Simple
# heuristic" comment block). Pinning this ensures the limitation is visible.
expect(identical(opt$country, "ETH"),
       "CURRENT: YAML country wins over CLI --country (1-deep heuristic)")
# YAML clearly wins over default for backend (top-level key from YAML survives).
expect(identical(opt$backend, "brms"),
       "YAML backend overrides default (non-default value)")

# Sanity: CLI-only key absent from YAML DOES get promoted.
tmp2 <- tempfile(fileext = ".yml")
writeLines("backend: brms", tmp2)
opt <- make_options_from_config(
  option_list = list(
    optparse::make_option("--country", default = NULL, type = "character")
  ),
  cli_argv = c("--config", tmp2, "--country", "BEN"),
  model = "size"
)
unlink(tmp2)
expect(identical(opt$country, "BEN"),
       "CLI --country promoted when YAML doesn't set country")

# --- 8. Granular slopes: map merges and survives ---------------------------
# Sanity that the recursive merge does not flatten the slopes: list-of-lists.
cat("\n[8] granular size$random_effects$slopes preserved\n")

opt <- resolve_yaml(paste(
  "backend: brms",
  "size:",
  "  random_effects:",
  "    slopes:",
  "      mean_precip_lag1:",
  "        adm1: true",
  sep = "\n"))
sl <- opt$size$random_effects$slopes
expect(is.list(sl) && !is.null(sl$mean_precip_lag1),
       "slopes map preserved")
expect(isTRUE(sl$mean_precip_lag1$adm1),
       "slopes$mean_precip_lag1$adm1 == TRUE survives the merge")

# --- 9. phi_re key survives deep merge ------------------------------------
cat("\n[9] phi_re key survives config merge\n")
opt <- resolve_yaml(paste(
  "backend: brms",
  "size:",
  "  phi_re: region_country_adm1",
  sep = "\n"))
expect(identical(opt$size$phi_re, "region_country_adm1"),
       "size$phi_re user value preserved through deep merge")

# --- 10. .fill_model_block() sibling-inheritance precedence for brms_cb ----
# re_uncorrelated must resolve from block[[be]] (explicit) >
# block[[plain]] (sibling) > defaults[[be]] (canonical), so
# make_batch_configs.R's strip_backend_siblings() can safely drop the plain
# `brms:` sibling from a brms_cb config once its value is merged forward.
cat("\n[10] .fill_model_block(): size/brms_cb inherits re_uncorrelated from the plain brms sibling\n")

opt10a <- resolve_yaml(paste(
  "size:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  sep = "\n"))
expect(isTRUE(opt10a$size$brms_cb$re_uncorrelated),
       "value set only on the plain brms sibling reaches size$brms_cb$re_uncorrelated")

opt10b <- resolve_yaml(paste(
  "size:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  "  brms_cb:",
  "    re_uncorrelated: false",
  sep = "\n"))
expect(identical(opt10b$size$brms_cb$re_uncorrelated, FALSE),
       "explicit brms_cb$re_uncorrelated wins over a conflicting plain-sibling value")

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
message("PASS  test_size_config_deep_merge")
