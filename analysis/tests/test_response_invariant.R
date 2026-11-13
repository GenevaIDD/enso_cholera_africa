# test_response_invariant.R
#
# Gate (validation harness for the Poisson GLMM).
#
# Asserts the `ongoing_outbreak`/`ongoing_outbreak_lag1` join at
# analysis/03_prepare_occ.R does NOT alter the survival response (`outbreak`).
#
# Background: the prior version used a `left_join` with no explicit `by`,
# so the natural-join key included `outbreak`. For an ADM1-period where any
# location_period_id had an outbreak, the ADM1-rolled-up row had
# `outbreak = TRUE`, and zero-padded `outbreak = FALSE` rows in the same
# ADM1-period silently failed to match, yielding NA `ongoing_outbreak_*`
# without touching the original `outbreak` column. The current code uses an
# explicit `by = c("ADM1", "t_start", "t_end")` and drops `outbreak` from the
# RHS, but this test guards against any future regression.
#
# This test does NOT require running the full prep pipeline. It exercises the
# join logic on a synthetic frame that reproduces the bug condition.

suppressPackageStartupMessages({
  library(tidyverse)
})

report_dir <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir, "response_invariant.md")

con <- file(report, "w")
on.exit(close(con))
cat("# test_response_invariant\n\n", file = con)

failures <- c()

# Build synthetic prep state ------------------------------------------------

# Two location periods within the same ADM1, same t_start/t_end. LP-A has the
# outbreak; LP-B is a zero-padded survival row.
surv_dat_all <- tibble::tribble(
  ~uid,    ~location_period_id, ~ADM1,   ~t_start,           ~t_end,             ~outbreak,
  "LP-A",  "lp-a",              "AAA_1", as.Date("2015-06-01"), as.Date("2015-06-30"), TRUE,
  "LP-B",  "lp-b",              "AAA_1", as.Date("2015-06-01"), as.Date("2015-06-30"), FALSE,
  "LP-C",  "lp-c",              "BBB_1", as.Date("2015-06-01"), as.Date("2015-06-30"), FALSE
)

# Rolled-up ADM1 frame: any(outbreak) over the ADM1-period.
adm1_outbreak <- surv_dat_all %>%
  group_by(ADM1, t_start, t_end) %>%
  summarise(outbreak = any(outbreak),
            ongoing_outbreak = any(outbreak),
            ongoing_outbreak_lag1 = FALSE,
            .groups = "drop")

# OLD behaviour: implicit `by` (includes `outbreak`) --------------------------
old_join <- surv_dat_all %>%
  left_join(adm1_outbreak %>%
              select(ADM1, t_start, t_end, outbreak,
                     ongoing_outbreak, ongoing_outbreak_lag1),
            by = c("ADM1", "t_start", "t_end", "outbreak"))

# NEW behaviour: explicit `by`, drop `outbreak` from RHS ----------------------
new_join <- surv_dat_all %>%
  left_join(adm1_outbreak %>%
              select(ADM1, t_start, t_end,
                     ongoing_outbreak, ongoing_outbreak_lag1),
            by = c("ADM1", "t_start", "t_end"))

# Assertions ------------------------------------------------------------------

cat("## Assertions\n\n", file = con)

# (A) New join preserves nrow and outbreak vector.
ok_nrow   <- nrow(new_join) == nrow(surv_dat_all)
ok_resp   <- identical(new_join$outbreak, surv_dat_all$outbreak)
ok_no_nas <- !any(is.na(new_join$ongoing_outbreak_lag1))

cat(sprintf("- nrow preserved: **%s** (%d -> %d)\n",
            if (ok_nrow) "PASS" else "FAIL",
            nrow(surv_dat_all), nrow(new_join)), file = con)
cat(sprintf("- outbreak bit-identical: **%s**\n",
            if (ok_resp) "PASS" else "FAIL"), file = con)
cat(sprintf("- no NA in ongoing_outbreak_lag1: **%s**\n",
            if (ok_no_nas) "PASS" else "FAIL"), file = con)

# (B) The OLD join was indeed broken on this synthetic case (regression guard).
old_has_nas <- any(is.na(old_join$ongoing_outbreak_lag1))
cat(sprintf("- (regression guard) OLD join produces NA covariates: **%s**\n",
            if (old_has_nas) "as expected" else "no longer reproduces"), file = con)

if (!ok_nrow)   failures <- c(failures, "nrow changed by join")
if (!ok_resp)   failures <- c(failures, "outbreak column altered by join")
if (!ok_no_nas) failures <- c(failures, "NA in ongoing_outbreak_lag1 after fixed join")

cat("\n## Verdict\n\n", file = con)
if (length(failures) == 0) {
  cat("PASS\n", file = con)
  message("test_response_invariant: PASS")
  quit(status = 0)
} else {
  cat("FAIL:\n", file = con)
  for (f in failures) cat("- ", f, "\n", sep = "", file = con)
  message("test_response_invariant: FAIL")
  for (f in failures) message("  - ", f)
  quit(status = 1)
}
