# test_occurrence_config_merge.R — HARD: occurrence-model config-layer merge
#   and convention-2 slope implication.
#
# Usage:
#   Rscript analysis/tests/test_occurrence_config_merge.R
#
# Details:
#   Pure config resolution -- no fit, no data. Pins the occurrence defaults
#   of make_options_from_config(model = "occurrence"), the deep merge of
#   occurrence.random_effects, convention-2 (slope flags imply climate fixed
#   effects), the covariates escape hatch, model-block stripping, explicit
#   per-level RE terms, cross-model flag equivalence and _cb backend-sibling
#   inheritance.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/formula_legacy.R") # build_formula_terms_legacy()
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

defaults <- get_default_options()
stopifnot(
  !is.null(defaults$occurrence),
  is.list(defaults$occurrence$random_effects),
  !is.null(defaults$occurrence$backend)
)

resolve_yaml <- function(yaml_body) {
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(yaml_body, tmp)
  make_options_from_config(
    option_list = list(),
    cli_argv = c("--config", tmp),
    model = "occurrence"
  )
}

# --- 1. Empty config: occurrence defaults populated --------------------------
cat("\n[1] empty YAML -> occurrence defaults populated\n")

opt <- resolve_yaml("occurrence:\n  backend: brms_cb\n")
expect(!is.null(opt$occurrence), "opt$occurrence is present")
expect(
  identical(opt$occurrence$backend, "brms_cb"),
  "backend matches YAML value"
)
expect(
  !is.null(opt$occurrence$random_effects),
  "occurrence$random_effects filled from defaults"
)
expect(
  !is.null(opt$occurrence$covariates),
  "occurrence$covariates filled from defaults"
)
expect(!is.null(opt$enso), "opt$enso (top-level shared) filled from defaults")

# --- 2. Partial random_effects: 2-deep merge preserves sibling defaults ------
cat("\n[2] partial occurrence$random_effects -> siblings preserved\n")

opt <- resolve_yaml(paste(
  "occurrence:",
  "  backend: brms",
  "  random_effects:",
  "    re_adm1_mean_precip_slope: true",
  sep = "\n"
))

expect(
  isTRUE(opt$occurrence$random_effects$re_adm1_mean_precip_slope),
  "user-set slope flag is TRUE"
)
expect(
  !is.null(opt$occurrence$random_effects$re_region_country_adm1_intercept),
  "sibling intercept RE default preserved (2-deep merge)"
)
expect(
  !is.null(opt$occurrence$random_effects$re_adm1_dev_precip_slope),
  "sibling dev_precip slope default preserved"
)

# --- 3. Convention-2: slope flags imply climate fixed effects ----------------
# build_formula_terms() reads occurrence$random_effects for climate covariates,
# not occurrence$covariates. A slope flag set to TRUE must appear in the
# returned fixed_terms WITHOUT the cov being listed under occurrence$covariates.
# Convention-2: build_formula_terms() emits climate fixed effects for active
# slope RE flags — no need to list covariates separately.
cat("\n[3] convention-2: slope flags imply climate fixed effects\n")

opt3 <- list(
  enso = defaults$enso,
  occurrence = list(
    backend = "brms",
    covariates = list(
      time_to_last_outbreak_cat = FALSE,
      ongoing_outbreak_lag1 = FALSE,
      enso = FALSE
    ),
    random_effects = list(
      re_adm1_mean_precip_slope = TRUE,
      re_adm1_dev_precip_slope = TRUE,
      re_adm1_mean_temp_daymean_slope = FALSE,
      re_adm1_dev_temp_daymean_slope = FALSE
    ),
    brms = list(re_uncorrelated = FALSE)
  )
)
# build_formula_terms is the occurrence-model formula assembler
terms <- build_formula_terms(opt3, n_country = 3L, n_region = 2L)
fixed <- terms$fixed_terms
re <- terms$re_terms

expect(
  "mean_precip" %in% fixed,
  "re_adm1_mean_precip_slope implies mean_precip as fixed effect"
)
expect(
  "dev_precip" %in% fixed,
  "re_adm1_dev_precip_slope implies dev_precip as fixed effect"
)
expect(
  !"mean_temp_daymean" %in% fixed,
  "inactive slope flag -> mean_temp_daymean NOT in fixed effects"
)
expect(
  !"dev_temp_daymean" %in% fixed,
  "inactive slope flag -> dev_temp_daymean NOT in fixed effects"
)

# --- 4. occurrence$covariates escape hatch for non-climate keys ---------------
# occurrence.covariates is the escape hatch for non-climate keys only
# (time_to_last_outbreak_cat, ongoing_outbreak_lag1, enso).
cat("\n[4] occurrence$covariates escape hatch for non-climate keys\n")

opt4 <- list(
  enso = defaults$enso,
  occurrence = list(
    backend = "brms",
    covariates = list(
      time_to_last_outbreak_cat = TRUE,
      ongoing_outbreak_lag1 = TRUE,
      enso = FALSE
    ),
    random_effects = list(
      re_adm1_mean_precip_slope = FALSE,
      re_adm1_dev_precip_slope = FALSE
    ),
    brms = list(re_uncorrelated = FALSE)
  )
)
terms4 <- build_formula_terms(opt4, n_country = 3L, n_region = 2L)
fixed4 <- terms4$fixed_terms
expect(
  "time_to_last_outbreak" %in% fixed4,
  "time_to_last_outbreak_cat in covariates -> fixed effect added"
)
expect(
  "ongoing_outbreak_lag1" %in% fixed4,
  "ongoing_outbreak_lag1 in covariates -> fixed effect added"
)
expect(
  !"mean_precip" %in% fixed4,
  "climate cov NOT in fixed effects when slope flag is FALSE"
)

# ENSO via covariates$enso toggle
opt4e <- opt4
opt4e$occurrence$covariates$enso <- TRUE
terms4e <- build_formula_terms(opt4e, n_country = 3L, n_region = 2L)
expect(
  any(grepl("^enso_lag", terms4e$fixed_terms)),
  "covariates$enso = TRUE -> enso_lag<N> in fixed effects"
)

# --- 5. No size: block when only occurrence: is present ---------------------
cat("\n[5] occurrence-only YAML -> no size: block in opt\n")

opt5 <- resolve_yaml(paste(
  "occurrence:",
  "  backend: brms",
  sep = "\n"
))
expect(
  is.null(opt5$size),
  "opt$size is NULL when YAML has only occurrence: block"
)

# --- 6. build_formula_terms() emits explicit per-level RE terms ---------------
# The new builder (resolve_slopes convention) must emit three separate terms for
# a re_adm1_* flag in a multi-region multi-country context, rather than one
# /-nested term.  The legacy builder emits the /-nested form.
cat("\n[6] build_formula_terms() emits explicit per-level RE terms\n")

opt6 <- list(
  enso = defaults$enso,
  occurrence = list(
    backend = "brms",
    covariates = list(
      time_to_last_outbreak_cat = FALSE,
      ongoing_outbreak_lag1 = FALSE,
      enso = FALSE
    ),
    random_effects = list(
      re_adm1_mean_precip_slope = TRUE,
      re_adm1_dev_precip_slope = FALSE,
      re_adm1_mean_temp_daymean_slope = FALSE,
      re_adm1_dev_temp_daymean_slope = FALSE
    ),
    brms = list(re_uncorrelated = FALSE)
  )
)
terms6 <- build_formula_terms(opt6, n_country = 3L, n_region = 2L)
terms6l <- build_formula_terms_legacy(opt6, n_country = 3L, n_region = 2L)
re6 <- terms6$re_terms
re6l <- terms6l$re_terms

# New builder: three explicit per-level terms
expect(
  any(grepl("mean_precip - 1 | AFRO_region$", re6, fixed = TRUE)) ||
    any(grepl("mean_precip - 1 | AFRO_region", re6, fixed = TRUE)),
  "new builder: region-level slope on mean_precip emitted"
)
expect(
  any(grepl("AFRO_region:country", re6, fixed = TRUE)),
  "new builder: country-level slope on mean_precip emitted"
)
expect(
  any(grepl("AFRO_region:country:ADM1", re6, fixed = TRUE)),
  "new builder: adm1-level slope on mean_precip emitted"
)
expect(
  !any(grepl("AFRO_region/country/ADM1", re6, fixed = TRUE)),
  "new builder: no /-nested form in RE terms"
)

# Legacy builder: single /-nested term
expect(
  any(grepl("AFRO_region/country/ADM1", re6l, fixed = TRUE)),
  "legacy builder: /-nested slope term present"
)

# Fixed effects are identical between the two builders
expect(
  identical(sort(terms6$fixed_terms), sort(terms6l$fixed_terms)),
  "new and legacy builders: identical fixed_terms"
)

# Single-country context: always full nesting (re_group_by_level does not
# collapse to bare ADM1) — region/country/adm1 slope terms are all emitted.
terms6_sc <- build_formula_terms(opt6, n_country = 1L, n_region = 1L)
re6_sc <- terms6_sc$re_terms
expect(
  any(grepl("mean_precip", re6_sc, fixed = TRUE)),
  "single-country: at least one mean_precip slope term emitted"
)
expect(
  any(grepl(
    "mean_precip - 1 | AFRO_region:country:ADM1",
    re6_sc,
    fixed = TRUE
  )),
  "single-country: full-nesting adm1 group emitted"
)
expect(
  any(grepl("mean_precip - 1 | AFRO_region)", re6_sc, fixed = TRUE)),
  "single-country: region group emitted (always full nesting)"
)

# --- 7. Cross-model flag equivalence: occurrence and size same covariate ------
# re_adm1_mean_precip_slope = TRUE on occurrence and size must produce
# the same fixed effect (mean_precip) and the same ADM1-level RE group string.
cat("\n[7] cross-model flag equivalence for shared covariate (mean_precip)\n")

opt7_occ <- list(
  enso = defaults$enso,
  occurrence = list(
    backend = "brms",
    covariates = list(
      time_to_last_outbreak_cat = FALSE,
      ongoing_outbreak_lag1 = FALSE,
      enso = FALSE
    ),
    random_effects = list(re_adm1_mean_precip_slope = TRUE),
    brms = list(re_uncorrelated = FALSE)
  )
)
opt7_size <- list(
  enso = defaults$enso,
  size = list(
    backend = "brms",
    family = "negbinomial",
    link = "log",
    admin_level = "ADM1",
    offset_pop = FALSE,
    offset_meanrate = FALSE,
    seasonal_term = "none",
    covariates = list(),
    random_effects = list(
      re_adm1_mean_precip_slope = TRUE,
      re_region_country_adm1_intercept = FALSE
    ),
    brms = defaults$size$brms
  )
)
source("analysis/utils/config.R") # ensure build_size_formula is available

occ_terms7 <- build_formula_terms(opt7_occ, n_country = 3L, n_region = 2L)
size_slopes7 <- resolve_slopes(opt7_size$size$random_effects, "ADM1")

expect(
  "mean_precip" %in% occ_terms7$fixed_terms,
  "cross-model: occurrence emits mean_precip as fixed effect"
)
expect(
  "mean_precip" %in% names(size_slopes7),
  "cross-model: resolve_slopes() returns mean_precip for size"
)

# Both models should produce an ADM1-level RE term for the same group string
occ_adm1_term <- grep(
  "mean_precip.*AFRO_region:country:ADM1",
  occ_terms7$re_terms,
  value = TRUE
)
expect(
  length(occ_adm1_term) > 0,
  "cross-model: occurrence emits AFRO_region:country:ADM1 slope term for mean_precip"
)
expect(
  isTRUE(size_slopes7$mean_precip$adm1),
  "cross-model: resolve_slopes() sets adm1=TRUE for size re_adm1_mean_precip_slope"
)

# --- 8. .fill_model_block() sibling-inheritance precedence for _cb backends --
# re_uncorrelated (and any field with the same shape) must resolve
# from block[[be]] (explicit) > block[[plain]] (sibling) > defaults[[be]]
# (canonical), so make_batch_configs.R's strip_backend_siblings() can safely
# drop the plain `brms:` sibling from a brms_cb config once its value has been
# merged forward -- see test_batch_config_block_stripping.R §[3] for the
# generator-side regression guard this precedence rule protects.
cat(
  "\n[8] .fill_model_block(): brms_cb inherits re_uncorrelated from the plain brms sibling\n"
)

opt8a <- resolve_yaml(paste(
  "occurrence:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  sep = "\n"
))
expect(
  isTRUE(opt8a$occurrence$brms_cb$re_uncorrelated),
  "value set only on the plain brms sibling reaches occurrence$brms_cb$re_uncorrelated"
)

opt8b <- resolve_yaml(paste(
  "occurrence:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  "  brms_cb:",
  "    re_uncorrelated: false",
  sep = "\n"
))
expect(
  identical(opt8b$occurrence$brms_cb$re_uncorrelated, FALSE),
  "explicit brms_cb$re_uncorrelated wins over a conflicting plain-sibling value"
)

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  quit(save = "no", status = 1)
}
message("PASS  test_occurrence_config_merge")
