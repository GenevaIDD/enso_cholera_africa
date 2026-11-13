# test_size_nested_slopes.R — verify per-cov × per-level slope wiring
#
# Pure-formula test (no MCMC): exercises
#   - resolve_slopes()            (size convention)
#   - build_size_formula()        (brms RHS)
# across the flag families that mirror the occurrence-side convention:
#
#   re_adm1_<cov>_slope            -> region + country + adm1
#   re_region_country_<cov>_slope  -> region + country
#   random_effects$slopes$<cov>    -> arbitrary per-level subset (escape hatch)
#
# Also covers:
#   - ENSO wiring (top-level opt$enso$use) → enso_lag<N> column emission
#   - admin_level == "ADM2" promotes re_adm1_* to also include adm2 in brms
#   - ttlo alias maps to time_to_last_outbreak
#   - Slope flag with no covariate set produces a validate_size_options() warning
#
# Run:  Rscript analysis/tests/test_size_nested_slopes.R
#
# Conventional CI exit: 0 = all asserts passed; 1 = any failure.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

# Local micro-harness — avoid external test framework dep.
PASS <- 0L; FAIL <- 0L; FAILS <- character(0)
expect <- function(cond, label) {
  if (isTRUE(cond)) {
    PASS <<- PASS + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    FAILS <<- c(FAILS, label)
    cat(sprintf("  FAIL  %s\n", label))
  }
}

# Helper — extract the RHS of the brms formula as a single string.
# deparse() inserts line-continuation padding; squish multiple spaces to one
# so simple `fixed = TRUE` substring matches work.
.squish <- function(s) gsub("\\s+", " ", paste(s, collapse = " "))
brms_rhs <- function(bf_obj) .squish(deparse(bf_obj$formula[[3L]]))

# Size-convention slope resolver, applied to an opt list.
size_slopes <- function(opt) resolve_slopes(opt$size$random_effects, opt$size$admin_level)

# Minimal opt builder. The size pipeline only ever reads slot names, never
# touches dat beyond passing it through, so the tiny tibble below is enough
# for brms formula assembly.
make_opt <- function(re = list(),
                     covariates_enso = FALSE,
                     enso_lag = 11L,
                     admin_level = "ADM1",
                     covariates = list()) {
  list(
    enso = list(lag_months = enso_lag, source = "oni",
                cache_dir  = "generated_data/"),
    size = list(
      backend = "brms", family = "negbinomial", link = "log",
      admin_level = admin_level,
      offset_pop = TRUE, offset_meanrate = TRUE,
      random_effects = re,
      # ENSO toggle lives in size.covariates (new schema); merge with caller overrides.
      covariates     = modifyList(list(enso = covariates_enso), covariates)
    )
  )
}

# Tiny synthetic size_dat with the columns the size formula builder reads.
mk_dat <- function() {
  n <- 60L
  tibble::tibble(
    suspected_cases       = sample(0:50, n, TRUE),
    AFRO_region           = rep(c("Western Africa", "Eastern Africa"), each = 30),
    country               = rep(c("BEN", "TGO", "ETH"), each = 20),
    ADM1                  = paste0(rep(c("BEN", "TGO", "ETH"), each = 20), "_",
                                   rep(1:4, length.out = n)),
    time_to_last_outbreak = factor(sample(c("[0,1]","(1,2]","(2,3]","(6,Inf]"),
                                          n, TRUE),
                                   levels = c("[0,1]","(1,2]","(2,3]","(6,Inf]")),
    ongoing_outbreak_lag1 = rbinom(n, 1, 0.3),
    mean_precip           = rnorm(n), dev_precip       = rnorm(n),
    mean_temp_daymean     = rnorm(n), dev_temp_daymean = rnorm(n),
    mean_precip_lag1      = rnorm(n),
    enso_lag6  = rnorm(n), enso_lag11 = rnorm(n),
    pop = 1e5, mean_rate = 1e-3
  )
}

# ============================================================================
# 1. Canonicaliser: prefix flag families resolve to expected level sets.
# ============================================================================
cat("\n[1] resolve_slopes() — prefix family expansion\n")

opt <- make_opt(list(
  re_adm1_mean_precip_slope        = TRUE,
  re_region_country_ongoing_slope  = TRUE,
  re_adm1_ttlo_slope               = TRUE
))
slopes <- size_slopes(opt)
expect(identical(sort(names(slopes)),
                 sort(c("mean_precip", "ongoing", "time_to_last_outbreak"))),
       "covariate names: re_adm1_* and re_region_country_* both detected, ttlo aliased")
expect(isTRUE(slopes$mean_precip$region) &&
       isTRUE(slopes$mean_precip$country) &&
       isTRUE(slopes$mean_precip$adm1) &&
       !isTRUE(slopes$mean_precip$adm2),
       "re_adm1_mean_precip_slope -> region+country+adm1 (no adm2 at admin_level=ADM1)")
expect(isTRUE(slopes$ongoing$region) &&
       isTRUE(slopes$ongoing$country) &&
       !isTRUE(slopes$ongoing$adm1) &&
       !isTRUE(slopes$ongoing$adm2),
       "re_region_country_ongoing_slope -> region+country only (shallow)")
expect(isTRUE(slopes$time_to_last_outbreak$adm1),
       "ttlo alias under re_adm1_ttlo_slope resolves to time_to_last_outbreak")

# admin_level == "ADM2" promotes re_adm1_* to add adm2
opt2 <- make_opt(list(re_adm1_mean_precip_slope = TRUE), admin_level = "ADM2")
sl2  <- size_slopes(opt2)
expect(isTRUE(sl2$mean_precip$adm2),
       "admin_level=ADM2: re_adm1_mean_precip_slope adds adm2 level")

# ============================================================================
# 2. Granular slopes: sub-map overrides prefix flags at the level granularity.
# ============================================================================
cat("\n[2] resolve_slopes() — granular slopes: map\n")

opt <- make_opt(list(
  re_adm1_mean_precip_slope = TRUE,
  slopes = list(
    mean_precip      = list(region = FALSE, country = FALSE, adm1 = TRUE),
    mean_precip_lag1 = list(country = TRUE, adm1 = TRUE)
  )
))
slopes <- size_slopes(opt)
expect(!isTRUE(slopes$mean_precip$region) &&
       !isTRUE(slopes$mean_precip$country) &&
       isTRUE(slopes$mean_precip$adm1),
       "slopes.mean_precip overrides legacy flag at level granularity")
expect(!isTRUE(slopes$mean_precip_lag1$region) &&
       isTRUE(slopes$mean_precip_lag1$country) &&
       isTRUE(slopes$mean_precip_lag1$adm1),
       "slopes.mean_precip_lag1 enables a covariate the prefix flags don't cover")

# A slopes entry with all FALSE drops out.
opt <- make_opt(list(slopes = list(
  mean_precip = list(region = FALSE, country = FALSE, adm1 = FALSE)
)))
expect(length(size_slopes(opt)) == 0L,
       "all-FALSE slopes entry is dropped")

# ============================================================================
# 3. brms cascade: each active (cov, level) pair appears as its own RE term.
# ============================================================================
cat("\n[3] build_size_formula() — brms RE term emission\n")

# 3a. Deep family produces three terms with the per-level groups.
opt <- make_opt(list(
  re_region_country_adm1_intercept = TRUE,
  re_adm1_mean_precip_slope        = TRUE
))
bf  <- build_size_formula(opt, mk_dat())
rhs <- brms_rhs(bf)
expect(grepl("(1 | AFRO_region/country/ADM1)", rhs, fixed = TRUE),
       "intercept RE preserves nested-slash form")
expect(grepl("(mean_precip - 1 | AFRO_region)", rhs, fixed = TRUE),
       "brms: region-level slope on mean_precip")
expect(grepl("(mean_precip - 1 | AFRO_region:country)", rhs, fixed = TRUE),
       "brms: country-level slope on mean_precip")
expect(grepl("(mean_precip - 1 | AFRO_region:country:ADM1)", rhs, fixed = TRUE),
       "brms: adm1-level slope on mean_precip")

# 3b. Shallow family stops at country.
opt <- make_opt(list(
  re_region_country_ongoing_slope = TRUE
))
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(grepl("(ongoing_outbreak_lag1 - 1 | AFRO_region)", rhs, fixed = TRUE) ||
       # In our test opt, the cov name from `re_region_country_ongoing_slope` is "ongoing".
       grepl("(ongoing - 1 | AFRO_region)", rhs, fixed = TRUE),
       "brms: shallow-family region-level slope emitted")
expect(!grepl(":ADM1", rhs, fixed = TRUE),
       "brms: shallow-family does NOT emit an ADM1-level slope")

# 3c. Granular slopes: country-only emits one term, no region, no adm1.
opt <- make_opt(list(slopes = list(
  mean_precip_lag1 = list(region = FALSE, country = TRUE, adm1 = FALSE)
)))
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(grepl("(mean_precip_lag1 - 1 | AFRO_region:country)", rhs, fixed = TRUE),
       "brms: granular country-only slope emitted")
expect(!grepl("(mean_precip_lag1 - 1 | AFRO_region)", rhs, fixed = TRUE) &&
       !grepl("ADM1", paste(grep("mean_precip_lag1", strsplit(rhs, " \\+ ")[[1]], value = TRUE), collapse = " "), fixed = TRUE),
       "brms: granular country-only slope does NOT emit region or adm1 terms")

# 3d. ADM2 admin level: re_adm1_* additionally produces the adm2 group.
opt <- make_opt(list(re_adm1_dev_precip_slope = TRUE), admin_level = "ADM2")
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(grepl("(dev_precip - 1 | AFRO_region:country:ADM1:ADM2)", rhs, fixed = TRUE),
       "brms: admin_level=ADM2 emits the adm2-level slope")

# ============================================================================
# 4. ENSO wiring on the brms side — mirrors occurrence-side.
# ============================================================================
cat("\n[4] build_size_formula() — ENSO\n")

# 4a. covariates$enso FALSE: enso slope flag is ignored.
opt <- make_opt(list(re_region_country_enso_lag11_slope = TRUE),
                covariates_enso = FALSE, enso_lag = 11L)
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(!grepl("enso", rhs, fixed = TRUE),
       "brms: enso slope ignored when covariates$enso is FALSE (with console message)")

# 4b. covariates$enso TRUE, lag_months = 11: emits enso_lag11 as fixed effect AND
#     the requested shallow-scope slope at region+country.
opt <- make_opt(list(re_region_country_enso_lag11_slope = TRUE),
                covariates_enso = TRUE, enso_lag = 11L)
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(grepl("enso_lag11", rhs, fixed = TRUE),
       "brms: ENSO active -> enso_lag11 column in fixed effects")
expect(grepl("(enso_lag11 - 1 | AFRO_region)", rhs, fixed = TRUE) &&
       grepl("(enso_lag11 - 1 | AFRO_region:country)", rhs, fixed = TRUE),
       "brms: ENSO slope cascades to region + country (shallow scope)")
expect(!grepl("(enso_lag11 - 1 | AFRO_region:country:ADM1)", rhs, fixed = TRUE),
       "brms: ENSO shallow-scope slope does NOT cascade to ADM1")

# 4c. covariates$enso TRUE alone (no slope flag): only fixed effect, no slope RE.
opt <- make_opt(list(), covariates_enso = TRUE, enso_lag = 6L)
rhs <- brms_rhs(build_size_formula(opt, mk_dat()))
expect(grepl("enso_lag6", rhs, fixed = TRUE) &&
       !grepl("(enso_lag6 - 1 |", rhs, fixed = TRUE),
       "brms: ENSO as fixed effect only when no slope flag set")

# ============================================================================
# 6. validate_size_options() — warns when slope cov is FALSE.
# ============================================================================
cat("\n[6] validate_size_options() — warnings\n")

# 6a. Slope on a covariate not listed in size$covariates should warn.
opt <- make_opt(list(re_adm1_mean_precip_slope = TRUE))
# size$covariates is empty (make_opt default), so mean_precip is not on
w <- tryCatch(withCallingHandlers(
  { validate_size_options(opt); character(0) },
  warning = function(w) { invokeRestart("muffleWarning"); conditionMessage(w) }),
  warning = function(w) conditionMessage(w))
got_warn <- length(w) > 0 && any(grepl("mean_precip", w))
# tryCatch flattens to single value; do a fresh capture:
ws <- character(0)
withCallingHandlers(validate_size_options(opt),
  warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
expect(any(grepl("mean_precip", ws)),
       "validate_size_options warns when slope cov is FALSE")

# 6b. ENSO slope requested but covariates$enso FALSE -> warn.
opt <- make_opt(list(re_region_country_enso_lag11_slope = TRUE),
                covariates_enso = FALSE, enso_lag = 11L)
ws <- character(0)
withCallingHandlers(validate_size_options(opt),
  warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
expect(any(grepl("enso", ws)),
       "validate_size_options warns when enso slope requested but enso.use FALSE")

# ============================================================================
# Summary
# ============================================================================
cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  cat("Failures:\n"); for (m in FAILS) cat("  - ", m, "\n", sep = "")
  quit(save = "no", status = 1)
}
