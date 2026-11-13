# test_results_occ_panel_tables.R -- the 09 occurrence-panel tables.
#
# HARD gate (see analysis/tests/run_all.R), quick: synthetic sidecars only.
#
# WHAT THIS GUARDS
# ----------------
#   * report_national_zero_audit (table D32): the country-year audit that 01
#     writes beside the panel's months file, tagged with rule + threshold; its
#     D32 CSV columns are pinned (the supplement reads them by name).
#   * every occurrence-panel table step is OFF when its sidecar does not exist,
#     so a panel without the sidecars (an older run) derives unchanged.
#
# Run: Rscript analysis/tests/test_results_occ_panel_tables.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "config", "paths", "provenance", "results_core", "national_zero"))
  source(sprintf("analysis/utils/%s.R", m))
source("analysis/utils/results_derive.R")

td <- tempfile("occtab"); dir.create(td)
opt <- list(baseline = "1990_2019", outdir = td,
            episodes = list(key = "gid_season", national_zero_fill = TRUE,
                            national_zero_threshold = 60L, national_zero_rule = "r2"))

cat("\n[D32] report_national_zero_audit\n")
au <- tibble(country = c("AAA", "AAA"), year = c(2019L, 2020L), who_cases_max = c(900, 30),
             who_exact_zero = FALSE, db_cases_deep = c(800, NA), db_who_ratio = c(0.9, NA),
             db_exceeds_who = FALSE, months_national_and_who = 0L, months_national_who_missing = 0L,
             months_who_no_national = 0L, months_who_below_threshold = c(0L, 8L),
             months_claimed = c(0L, 8L), head_unclaimed = c(0L, 4L), tail_unclaimed = 0L,
             auto_excluded_exact = 0L, auto_excluded_threshold = 0L, static_excluded = FALSE)
ap <- file.path(td, "national_zero_months_t60_r2.audit.csv")
utils::write.csv(au, ap, row.names = FALSE)
ctx <- list2env(list(opt = opt, NZ_AUDIT_PATH = ap), parent = emptyenv())
r <- derive_report_national_zero_audit(ctx)
D32_COLS <- c("rule", "zero_threshold", "country", "year", "who_cases_max", "who_exact_zero",
              "db_cases_deep", "db_who_ratio", "db_exceeds_who", "months_national_and_who",
              "months_national_who_missing", "months_who_no_national", "months_who_below_threshold",
              "months_claimed", "head_unclaimed", "tail_unclaimed", "auto_excluded_exact",
              "auto_excluded_threshold", "static_excluded")
expect(identical(names(r$audit), D32_COLS), "D32 columns pinned (rule, zero_threshold + the 01 audit)")
expect(all(r$audit$rule == "r2") && all(r$audit$zero_threshold == 60L), "rows tagged with rule and threshold")
expect(r$totals$months_claimed == 8L && r$totals$head_unclaimed == 4L && r$totals$country_years == 2L,
       "totals over the audit")
expect(identical(names(r$params)[1:2], c("rule", "zero_threshold")), "params carry the rule")
expect(identical(names(nationalZeroAudit(tibble(country = character(), month = as.Date(character()),
                                                source = character()),
                                         tibble(country = "AAA", year = 2020L, suspected_cases = 1,
                                                is_zero = FALSE, suspected_cases_max = 1),
                                         NULL, 2020L)),
                 D32_COLS[-(1:2)]), "nationalZeroAudit() writes exactly the D32 audit columns")
expect(identical(.RESULTS_STAMP_BACKCOMPAT[["inputs.episode_national_zero_rule"]], "") &&
         identical(get_default_options()$results$inputs$episode_national_zero_rule, ""),
       "results.inputs.episode_national_zero_rule: default and back-compat \"\"")

cat("\n[off] steps absent without sidecars\n")
ids_on  <- names(results_steps(list2env(list(FIG1_TRACKS = "occurrence", HC_SCOPE = list(enable = FALSE),
                                             LAG_AVG = list(enable = FALSE), NZ_AUDIT_PATH = ap),
                                        parent = emptyenv())))
ids_off <- names(results_steps(list2env(list(FIG1_TRACKS = "occurrence", HC_SCOPE = list(enable = FALSE),
                                             LAG_AVG = list(enable = FALSE), NZ_AUDIT_PATH = ""),
                                        parent = emptyenv())))
expect("report_national_zero_audit" %in% ids_on && !"report_national_zero_audit" %in% ids_off &&
         identical(setdiff(ids_on, ids_off), "report_national_zero_audit"),
       "D32 step on only with the audit; nothing else changes")

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
