# test_size_seasonal_term.R — verify the opt$size$seasonal_term contract
# across all four documented values.
#
# Hard-fail gate. opt$size$seasonal_term is documented in 07a's CLI help and
# in config_defaults.yml as one of:
#
#   "none" | "factor_month" | "cyclic" | "harmonic2"
#
# The size_seasonal_sensitivity BATCH_REGISTRY entry generates one config
# per value, so any code path that branches on seasonal_term must accept
# all four without errors.
#
# Failure modes this guards against:
#   - a refactor of makeRunFile that drops the "seas-<term>" filename
#     suffix and silently collapses 4 distinct sensitivity runs to 1 path;
#   - a future build_size_formula() change
#     that errors on a non-"none" value (the current implementation treats
#     non-"none" as a no-op in the formula RHS — this test pins that, so
#     when the real implementation lands, the regression of "term silently
#     ignored" surfaces here);
#   - a regression in the size_seasonal_sensitivity batch that drops one
#     of the four documented values.
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/base.R")   # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/paths.R")   # makeSizeRunFile() -- used directly
source("analysis/utils/provenance.R")   # transitive dependency of paths.R (utils/README.md)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

TERMS <- c("none", "factor_month", "cyclic", "harmonic2")

mk_opt <- function(term, backend = "brms") {
  list(
    baseline = "1990_2019",
    country  = "all",
    outdir   = "generated_data/",
    enso     = list(use = FALSE),
    covariates = list(mean_precip = TRUE),
    size = list(
      backend         = backend,
      family          = "negbinomial",
      link            = "log",
      admin_level     = "ADM1",
      seasonal_term   = term,
      offset_pop      = TRUE,
      offset_meanrate = TRUE,
      prior_b         = "normal(0, 0.5)",
      prior_intercept = "normal(-3, 1)",
      prior_sd        = "normal(0, 0.5)",
      random_effects  = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm1_mean_precip_slope        = TRUE
      ),
      covariates_override = list()
    )
  )
}

mk_dat <- function() {
  set.seed(7)
  n <- 60L
  tibble::tibble(
    suspected_cases = sample(0:20, n, TRUE),
    AFRO_region     = rep(c("Western Africa", "Eastern Africa"), each = 30),
    country         = rep(c("BEN", "TGO", "ETH", "KEN"), each = 15),
    ADM1            = paste0(rep(c("BEN","TGO","ETH","KEN"), each = 15),
                             "_", rep(1:5, length.out = n)),
    time_to_last_outbreak = factor(
      sample(c("[0,1]","(6,Inf]"), n, TRUE),
      levels = c("[0,1]","(6,Inf]")),
    mean_precip = rnorm(n),
    pop = 1e5, mean_rate = 1e-3,
    # Month column — populate so any future seasonal_term implementation
    # has a real predictor to consume.
    month       = sample.int(12L, n, replace = TRUE)
  )
}

# --- 1. Filename roundtrip per term ---------------------------------------
# makeSizeRunFile encodes seasonal_term as a "seas-<term>" suffix when
# different from the default ("none"). Each of the four documented values
# must produce a distinct path.
cat("\n[1] makeSizeRunFile encodes seasonal_term uniquely\n")

paths <- vapply(TERMS, function(t) {
  makeSizeRunFile(mk_opt(t), kind = "size_fit")
}, character(1))

expect(length(unique(paths)) == length(TERMS),
       "all four seasonal_term values produce distinct paths")

# "none" is the default — its path must NOT carry the seas- token.
expect(!grepl("seas-none", paths[["none"]], fixed = TRUE),
       "seasonal_term='none' omits the seas- suffix (default)")

for (t in setdiff(TERMS, "none")) {
  expect(grepl(paste0("seas-", t), paths[[t]], fixed = TRUE),
         sprintf("seasonal_term='%s' adds 'seas-%s' to path", t, t))
}

# --- 2. build_size_formula accepts every documented term without error ---
# Current contract: non-"none" values do not yet contribute terms to the
# formula RHS (the implementation is a known TODO). The
# function MUST NOT error on any of the four values — otherwise the
# size_seasonal_sensitivity batch can't run at all.
cat("\n[2] build_size_formula accepts every documented term\n")

dat <- mk_dat()
for (t in TERMS) {
  bf_obj <- tryCatch(build_size_formula(mk_opt(t, "brms"), dat),
                     error = function(e) e)
  expect(!inherits(bf_obj, "error"),
         sprintf("build_size_formula(seasonal_term='%s') does not error", t))
  if (!inherits(bf_obj, "error")) {
    rhs_str <- paste(deparse(bf_obj$formula[[3L]]), collapse = " ")
    rhs_str <- gsub("\\s+", " ", rhs_str)
    if (t == "none") {
      expect(!grepl("month", rhs_str, fixed = TRUE),
             "none: no 'month' term in formula RHS")
    } else if (t == "factor_month") {
      expect(grepl("factor(month)", rhs_str, fixed = TRUE),
             "factor_month: factor(month) in formula RHS")
    } else if (t == "cyclic") {
      expect(grepl("s(month", rhs_str, fixed = TRUE),
             "cyclic: s(month, ...) cyclic spline in formula RHS")
    } else if (t == "harmonic2") {
      expect(grepl("sin_month", rhs_str, fixed = TRUE),
             "harmonic2: sin_month harmonic term in formula RHS")
    }
  }
}

# --- 4. BATCH_REGISTRY size_seasonal_sensitivity carries all four terms --
# Regression guard against the batch entry losing one of the documented
# values (which would silently shrink the sensitivity sweep).
cat("\n[4] BATCH_REGISTRY$size_seasonal_sensitivity covers all four terms\n")

# Load the batch registry without triggering CLI side-effects.
.local_env <- new.env()
sys.source("analysis/utils/make_batch_configs.R", envir = .local_env)
spec <- .local_env$BATCH_REGISTRY[["size_seasonal_sensitivity"]]
expect(!is.null(spec),
       "BATCH_REGISTRY contains 'size_seasonal_sensitivity'")
if (!is.null(spec)) {
  grid_terms <- spec$grid[["size.seasonal_term"]]
  expect(setequal(grid_terms, TERMS),
         sprintf("grid$size.seasonal_term == {%s}", paste(TERMS, collapse = ", ")))
}

# --- 5. Committed batch dir has one config per term ----------------------
# Each generated YAML's $size$seasonal_term field must match the filename
# token, otherwise downstream postprocessing keyed on filename will read
# the wrong term.
cat("\n[5] committed YAMLs in size_seasonal_sensitivity round-trip\n")

batch_dir <- "analysis/configs/size_seasonal_sensitivity"
if (dir.exists(batch_dir)) {
  yml_files <- list.files(batch_dir, pattern = "\\.yml$", full.names = TRUE)
  yml_files <- setdiff(yml_files, file.path(batch_dir, "manifest.yml"))
  expect(length(yml_files) == length(TERMS),
         sprintf("%d YAMLs in size_seasonal_sensitivity (one per term)",
                 length(TERMS)))
  seen <- character()
  for (f in yml_files) {
    cfg <- yaml::read_yaml(f)
    term_in_cfg  <- cfg$size$seasonal_term
    # Filename token (after gsub of non-alnum, the names use lowercase).
    fn_token <- sub(".*_sst-([^.]+)\\.yml$", "\\1", basename(f))
    # The DIM_ABBREV makes "factor_month" -> "factormonth" (underscore
    # stripped); normalise back for comparison.
    fn_norm  <- gsub("[^A-Za-z0-9]", "", as.character(term_in_cfg))
    expect(identical(fn_token, fn_norm),
           sprintf("%s: filename token '%s' matches cfg$size$seasonal_term '%s'",
                   basename(f), fn_token, term_in_cfg))
    seen <- c(seen, term_in_cfg)
  }
  expect(setequal(seen, TERMS),
         sprintf("committed batch covers all four terms: %s",
                 paste(sort(seen), collapse = ", ")))
} else {
  cat("  SKIP  analysis/configs/size_seasonal_sensitivity not present\n")
}

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
message("PASS  test_size_seasonal_term")
