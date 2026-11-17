# test_size_phi_re_formula.R — build_size_formula() phi RE wiring and the
#   validate_size_options() phi_re guard.
#
# Usage:
#   Rscript analysis/tests/test_size_phi_re_formula.R
#
# Details:
#   phi_re = "constant" gives a plain brmsformula with no phi sub-formula;
#   "region" / "region_country" / "region_country_adm1" add a phi sub-formula
#   with the expected grouping term; a bad value throws. Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(brms)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

# Minimal opt builder — no data needed for formula construction.
mk_opt <- function(phi_re = "constant", backend = "brms") {
  list(
    enso = list(
      lag_months = 11L,
      source = "oni",
      cache_dir = "generated_data/"
    ),
    size = list(
      backend = backend,
      family = "negbinomial",
      link = "log",
      admin_level = "ADM1",
      offset_pop = TRUE,
      offset_meanrate = TRUE,
      seasonal_term = "none",
      phi_re = phi_re,
      brms = list(re_uncorrelated = FALSE),
      random_effects = list(re_region_country_adm1_intercept = TRUE),
      covariates = list()
    )
  )
}

# Dummy data — 2 regions x 2 countries x 2 ADM1. Column values matter: the
# formula builders drop RE grouping levels with <2 distinct values
# (filter_re_terms_for_data), so a 0-row frame would strip every phi term.
dummy_dat <- expand.grid(
  ADM1 = c("A1", "A2"),
  country = c("C1", "C2"),
  AFRO_region = c("R1", "R2"),
  stringsAsFactors = FALSE
)
dummy_dat$country <- paste(dummy_dat$AFRO_region, dummy_dat$country, sep = "_")
dummy_dat$ADM1 <- paste(dummy_dat$country, dummy_dat$ADM1, sep = "_")
dummy_dat$suspected_cases <- 100L
dummy_dat$pop <- 1e5
dummy_dat$mean_rate <- 1e-3

# Helper: extract the text of the phi sub-formula from a brmsformula, or NULL.
# brmsformula stores distributional formulas in bf_obj$pforms as plain formula
# objects (class "formula"), not as sub-lists. brms names the negbinomial
# dispersion dpar `shape`, so build_size_formula() attaches `shape ~ ...`.
phi_formula_text <- function(bf_obj) {
  pf <- bf_obj$pforms
  if (is.null(pf) || is.null(pf[["shape"]])) {
    return(NULL)
  }
  paste(deparse(pf[["shape"]]), collapse = " ")
}

# --- 1. constant: plain brmsformula, no phi sub-formula -------------------
cat("\n[1] phi_re = constant: no phi sub-formula\n")
opt <- mk_opt("constant")
bf <- build_size_formula(opt, dummy_dat)
expect(inherits(bf, "brmsformula"), "constant: result is brmsformula")
expect(is.null(phi_formula_text(bf)), "constant: no phi sub-formula attached")

# --- 2. region: phi ~ 1 + (1 | AFRO_region) --------------------------------
cat("\n[2] phi_re = region\n")
opt <- mk_opt("region")
bf <- build_size_formula(opt, dummy_dat)
expect(inherits(bf, "brmsformula"), "region: result is brmsformula")
phi_txt <- phi_formula_text(bf)
expect(!is.null(phi_txt), "region: phi sub-formula present")
expect(
  grepl("AFRO_region", phi_txt, fixed = TRUE),
  "region: grouping term contains AFRO_region"
)
expect(
  !grepl("country", phi_txt, fixed = TRUE),
  "region: grouping term does NOT contain country"
)

# --- 3. region_country: phi ~ 1 + (1 | AFRO_region/country) ---------------
cat("\n[3] phi_re = region_country\n")
opt <- mk_opt("region_country")
bf <- build_size_formula(opt, dummy_dat)
phi_txt <- phi_formula_text(bf)
expect(!is.null(phi_txt), "region_country: phi sub-formula present")
expect(
  grepl("AFRO_region/country", phi_txt, fixed = TRUE),
  "region_country: grouping is AFRO_region/country"
)
expect(
  !grepl("ADM1", phi_txt, fixed = TRUE),
  "region_country: grouping does NOT contain ADM1"
)

# --- 4. region_country_adm1: phi ~ 1 + (1 | AFRO_region/country/ADM1) ----
cat("\n[4] phi_re = region_country_adm1\n")
opt <- mk_opt("region_country_adm1")
bf <- build_size_formula(opt, dummy_dat)
phi_txt <- phi_formula_text(bf)
expect(!is.null(phi_txt), "region_country_adm1: phi sub-formula present")
expect(
  grepl("AFRO_region/country/ADM1", phi_txt, fixed = TRUE),
  "region_country_adm1: grouping is AFRO_region/country/ADM1"
)

# --- 5. validate_size_options: invalid phi_re throws ----------------------
cat("\n[5] validate_size_options: bad phi_re value throws\n")
opt <- mk_opt("bad_value")
threw <- tryCatch(
  {
    validate_size_options(opt)
    FALSE
  },
  error = function(e) TRUE
)
expect(isTRUE(threw), "phi_re='bad_value' raises an error")

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  quit(save = "no", status = 1)
}
message("PASS  test_size_phi_re_formula")
