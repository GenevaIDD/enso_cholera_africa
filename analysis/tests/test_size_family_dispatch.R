# test_size_family_dispatch.R — verify size$family routes correctly.
#
# Hard-fail gate. The size pipeline supports two response families:
#
#   negbinomial -> brms::negbinomial("log")
#   poisson     -> brms::poisson("log")
#
# 04d (brms) calls match.fun(opt$size$family %||% ...)(opt$size$link %||% ...).
#
# Tests:
#   1. brms (07a) static: match.fun(<family>)(<link>) is the call form used
#      at the brm() site, and both supported family names resolve to a
#      brmsfamily object with the requested link.
#   2. filename roundtrip: makeRunFile() includes the family in the basename
#      so swapping family produces a distinct output path.
#
# All checks are static (no MCMC) so the test is fast.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/base.R")   # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/paths.R")   # makeRunFile() -- used directly
source("analysis/utils/provenance.R")   # transitive dependency of paths.R (utils/README.md)

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

# Local micro-harness for per-assertion reporting.
PASS <- 0L; FAIL <- 0L

# --- 1. brms runner (04d) -------------------------------------------------
cat("\n[1] 04d (brms) family routing\n")

path_07a <- "analysis/04_fit_siz_brms.R"
src_07a  <- paste(readLines(path_07a), collapse = "\n")

expect(grepl("match\\.fun\\(opt\\$size\\$family",
             src_07a),
       "04d uses match.fun(opt$size$family...) at brm() site")

# match.fun resolution should succeed for both supported families and return
# a brmsfamily object with the requested link.
fam_nb  <- match.fun("negbinomial")("log")
expect(inherits(fam_nb, "brmsfamily") && identical(fam_nb$family, "negbinomial"),
       "match.fun('negbinomial')('log') -> brmsfamily(negbinomial)")
expect(identical(fam_nb$link, "log"),
       "negbinomial family carries link = 'log'")

fam_po  <- match.fun("poisson")("log")
expect(inherits(fam_po, "family") &&
       identical(fam_po$family, "poisson") &&
       identical(fam_po$link, "log"),
       "match.fun('poisson')('log') -> family(poisson, link=log)")

# opt builder for the makeRunFile check below.
mk_opt <- function(fam) list(
  enso = list(lag_months = 11L, source = "oni", cache_dir = "generated_data/"),
  covariates = list(enso = FALSE, mean_precip = TRUE),
  size = list(backend="brms", family=fam, link="log", admin_level="ADM1",
              offset_pop=TRUE, offset_meanrate=TRUE,
              random_effects = list(re_region_country_adm1_intercept = TRUE,
                                    re_adm1_mean_precip_slope = TRUE),
              covariates = list())
)

# --- 2. makeRunFile basename carries the family --------------------------
cat("\n[2] makeRunFile basename uniqueness by family\n")

opt_nb <- mk_opt("negbinomial"); opt_nb$baseline <- "1990_2019"
opt_nb$country <- "ETH"; opt_nb$outdir <- "generated_data/"
opt_po <- opt_nb; opt_po$size$family <- "poisson"

fp_nb <- makeRunFile(opt_nb, kind = "size_fit")
fp_po <- makeRunFile(opt_po, kind = "size_fit")
expect(!identical(fp_nb, fp_po),
       "makeRunFile: family difference produces distinct basenames")
expect(grepl("negbinomial", fp_nb) && grepl("poisson", fp_po),
       "makeRunFile: family name appears in the basename")

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
message("PASS  test_size_family_dispatch")
