# test_size_priors.R — HARD: build_size_priors() prior-class assembly.
#
# Usage:
#   Rscript analysis/tests/test_size_priors.R
#
# Details:
#   Always present: Intercept (class "Intercept"), sd (class "sd"). The "b"
#   class is present ONLY when there is at least one fixed-effect slope: brms
#   rejects a "b" prior otherwise ("prior class b not in model"), so a bug in
#   the has_slopes branch silently breaks every fit. A fixed-effect slope
#   comes from (a) any prefix-flag slope RE (the slope cov is also a fixed
#   effect), (b) any size.covariates toggle set to TRUE, or (c)
#   covariates$enso = TRUE.

suppressPackageStartupMessages({
  library(brms)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/priors.R") # brms prior builders
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

mk_opt <- function(re = list(), covariates = list(), covariates_enso = FALSE) {
  list(
    enso = list(
      lag_months = 11L,
      source = "oni",
      cache_dir = "generated_data/"
    ),
    size = list(
      backend = "brms",
      family = "negbinomial",
      link = "log",
      admin_level = "ADM1",
      offset_pop = TRUE,
      offset_meanrate = TRUE,
      brms = list(
        prior_b = "normal(0, 0.5)",
        prior_intercept = "normal(-3, 1)",
        prior_sd = "normal(0, 0.5)"
      ),
      random_effects = re,
      covariates = modifyList(list(enso = covariates_enso), covariates)
    )
  )
}
classes <- function(p) sort(unique(as.character(p$class)))

# --- 1. No slopes, no covariates_override: NO 'b' class ------------------
cat("\n[1] intercept-only model\n")
opt <- mk_opt(re = list(re_region_country_adm1_intercept = TRUE))
p <- build_size_priors(opt)
cls <- classes(p)
expect(
  identical(cls, c("Intercept", "sd")),
  "intercept-only: only Intercept + sd classes (no 'b')"
)

# --- 2. Slope RE (prefix flag) alone: 'b' class IS present ---------------
cat("\n[2] slope from prefix flag\n")
opt <- mk_opt(
  re = list(
    re_region_country_adm1_intercept = TRUE,
    re_adm1_mean_precip_slope = TRUE
  )
)
p <- build_size_priors(opt)
cls <- classes(p)
expect("b" %in% cls, "re_adm1_mean_precip_slope -> 'b' class present")
expect(
  all(c("Intercept", "sd") %in% cls),
  "Intercept + sd still present alongside 'b'"
)

# --- 3. Slope from the granular slopes: map: 'b' class IS present --------
cat("\n[3] slope from granular slopes: map\n")
opt <- mk_opt(
  re = list(
    slopes = list(
      mean_precip_lag1 = list(adm1 = TRUE)
    )
  )
)
p <- build_size_priors(opt)
expect(
  "b" %in% classes(p),
  "granular slopes.<cov>.adm1=TRUE -> 'b' class present"
)

# --- 4. size.covariates alone (no slope RE): 'b' class IS present ---------
cat("\n[4] fixed effect via size.covariates only\n")
opt <- mk_opt(
  re = list(re_region_country_adm1_intercept = TRUE),
  covariates = list(mean_precip = TRUE)
)
p <- build_size_priors(opt)
expect("b" %in% classes(p), "size.covariates=TRUE alone -> 'b' class present")

# --- 5. All-FALSE in size.covariates: NO 'b' class -----------------------
cat("\n[5] size.covariates with no TRUE entry\n")
opt <- mk_opt(
  re = list(re_region_country_adm1_intercept = TRUE),
  covariates = list(mean_precip = FALSE, dev_precip = FALSE)
)
p <- build_size_priors(opt)
expect(!"b" %in% classes(p), "all-FALSE size.covariates -> no 'b' class")

# --- 6. ENSO fixed effect via covariates$enso emits 'b' class ------------
cat("\n[6] ENSO fixed effect flag\n")
opt <- mk_opt(
  re = list(
    re_region_country_adm1_intercept = TRUE,
    re_region_country_enso_lag11_slope = TRUE
  ),
  covariates_enso = TRUE
)
p <- build_size_priors(opt)
expect(
  "b" %in% classes(p),
  "re_region_country_enso_lag11_slope with covariates$enso=TRUE -> 'b' class present"
)

# --- 7. Custom prior strings round-trip ----------------------------------
cat("\n[7] custom prior strings\n")
opt <- mk_opt(
  re = list(
    re_region_country_adm1_intercept = TRUE,
    re_adm1_mean_precip_slope = TRUE
  )
)
opt$size$brms$prior_b <- "normal(0, 0.25)"
opt$size$brms$prior_intercept <- "normal(-5, 2)"
opt$size$brms$prior_sd <- "exponential(1)"
p <- build_size_priors(opt)
# brms stores the prior string under $prior; check ours are present.
priors_str <- as.character(p$prior)
expect(
  any(grepl("normal(0, 0.25)", priors_str, fixed = TRUE)),
  "custom prior_b string round-trips"
)
expect(
  any(grepl("normal(-5, 2)", priors_str, fixed = TRUE)),
  "custom prior_intercept string round-trips"
)
expect(
  any(grepl("exponential(1)", priors_str, fixed = TRUE)),
  "custom prior_sd string round-trips"
)

# --- 8. phi_re = "region_country_adm1" adds shape-class priors ----------------
cat("\n[8] phi_re = region_country_adm1 adds Intercept+sd with dpar=shape\n")
opt <- mk_opt(re = list(re_region_country_adm1_intercept = TRUE))
opt$size$phi_re <- "region_country_adm1"
p <- build_size_priors(opt)
phi_rows <- p[as.character(p$dpar) == "shape", ]
phi_cls <- sort(unique(as.character(phi_rows$class)))
expect(
  identical(phi_cls, c("Intercept", "sd")),
  "phi_re=region_country_adm1: Intercept(dpar=shape) + sd(dpar=shape) present"
)
main_rows <- p[as.character(p$dpar) == "", ]
expect(
  !"shape" %in% as.character(main_rows$class),
  "main priors unaffected — no spurious shape contamination"
)

# --- 9. phi_re = "constant" (default) emits no shape-class priors -----------
cat("\n[9] phi_re = constant: no shape-class priors\n")
opt <- mk_opt(re = list(re_region_country_adm1_intercept = TRUE))
opt$size$phi_re <- "constant"
p <- build_size_priors(opt)
phi_rows <- p[as.character(p$dpar) == "shape", ]
expect(nrow(phi_rows) == 0L, "phi_re=constant: no shape-dpar priors emitted")

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  quit(save = "no", status = 1)
}
message("PASS  test_size_priors")
