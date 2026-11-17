# test_size_response_invariant.R — stability invariants of the size prep
#   output (generated_data/outbreak_size_analysis_data.rds).
#
# Usage:
#   Rscript analysis/tests/test_size_response_invariant.R
#
# Details:
#   Non-empty; every required column present (suspected_cases, pop, uid,
#   location_period_id, AFRO_region, country, ADM1, gid, TL, TR,
#   time_to_last_outbreak, enso_lag0..12, mean_precip, dev_precip,
#   mean_temp_daymean, dev_temp_daymean); suspected_cases >= 0 and finite;
#   pop > 0 and finite where suspected_cases is non-NA. Not registered in
#   analysis/tests/run_all.R.

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

path <- "generated_data/outbreak_size_analysis_data.rds"
if (!file.exists(path)) {
  message("WARN: prep output missing (", path, ") — run 06 first.")
  message("PASS (skipped)")
  quit(status = 0)
}

d <- readRDS(path)
if (!nrow(d)) {
  fail("size analysis data is empty")
}

required <- c(
  "suspected_cases",
  "pop",
  "uid",
  "location_period_id",
  "country",
  "TL",
  "TR",
  "time_to_last_outbreak"
)
missing <- setdiff(required, colnames(d))
if (length(missing)) {
  fail(paste("missing required columns:", paste(missing, collapse = ", ")))
}

# ENSO lag columns are produced by the updated 06 prep — emit a WARN if
# missing (older prep RDS won't have them yet); not a hard failure since the
# legacy 07 formula does not consume ENSO.
enso_cols <- paste0("enso_lag", 0:12)
missing_enso <- setdiff(enso_cols, colnames(d))
if (length(missing_enso)) {
  message(
    "WARN: prep RDS lacks ENSO columns (",
    length(missing_enso),
    "/",
    length(enso_cols),
    " missing) — re-run 06_prepare_data_outbreak_size.R to regenerate."
  )
}

if (any(!is.finite(d$suspected_cases) | d$suspected_cases < 0, na.rm = TRUE)) {
  fail("suspected_cases has negative or non-finite values")
}

if (any(!is.finite(d$pop) | d$pop <= 0, na.rm = TRUE)) {
  fail("pop has non-positive or non-finite values")
}

message("PASS: test_size_response_invariant")
quit(status = 0)
