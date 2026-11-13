# analysis/tests/test_duration_span_stamp.R
#
# HARD gate (quick, synthetic data only; no fit, no data/ access, no
# generated_data/): the duration policy stamp's episode-span and guard-count
# fields and the deploy lag-floor check.
#
#   *    The deploy floor `deploy.lag_selection.occ_floor_offset` is meant to
#        be the median whole-month episode span (onset month to last active
#        month, n_months - 1). 03_prepare_dur.R stamps it
#        (durationSpanSummary(), utils/dur_data.R) and
#        loo_lags.R::check_occ_floor_offset() compares the two.
#   *    occ_data.R::guardOutbreakEpisodes() has no caller; the duration track
#        stamps what it WOULD drop (outbreakEpisodeGuardFlags(), dur_data.R).
#
#   [1] durationSpanSummary(): span quantiles and elapsed months; refusal.
#   [2] check_occ_floor_offset(): silent on a match or with no floor, warning
#       on a mismatch, message when the stamp lacks the field.
#   [3] outbreakEpisodeGuardFlags(): `any` selects exactly the rows
#       guardOutbreakEpisodes() drops; the classes are the guard's.
#   [4] no caller of guardOutbreakEpisodes() in analysis/ outside
#       analysis/scratchpads/ (helpers_parse.R::callers_of()), its own
#       definition and this test.
#   [5] 03_prepare_dur.R --restamp on a synthetic panel in a temp dir: the
#       stamp gains the span and guard fields and keeps every other one, the
#       panel is untouched, and a stamp that does not describe the panel is
#       refused.
#   [6] deploy_core.R::deploy_select_lags() calls the check:
#       on mocked rankings and a duration panel stamp in a temp dir, a
#       mismatching span warns, a matching one is silent, a missing stamp
#       only messages, no occ_floor_offset reads nothing; the lags are the
#       same in every case.
#   [7] results_derive.R::report_cascade() keeps the span fields of the
#       duration stamp (median_span_months, span_months_q25 / _q75) for S1.
#
# Run: Rscript analysis/tests/test_duration_span_stamp.R

suppressPackageStartupMessages({
  library(dplyr)
})

source("analysis/utils/base.R")
source("analysis/utils/occ_data.R")   # guardOutbreakEpisodes() (the reference)
source("analysis/utils/dur_data.R")   # durationSpanSummary(), outbreakEpisodeGuardFlags()
source("analysis/utils/loo_lags.R")   # lag_selection_cfg(), check_occ_floor_offset()
source("analysis/tests/helpers_parse.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# duration stamp: episode span, guard counts, lag-floor check\n")

# Capture the conditions an expression signals, muffled.
conds_of <- function(expr) {
  w <- character(0); m <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(value = val, warnings = w, messages = m)
}

# --- [1] durationSpanSummary ------------------------------------------------
cat("\n[1] durationSpanSummary()\n")
d0 <- as.Date("2020-01-15")
ep <- tibble::tibble(
  n_months = c(1L, 2L, 5L, 5L, 8L),            # span 0, 1, 4, 4, 7
  TL = d0,
  TR = d0 + c(10, 40, 120, 130, 220))
s <- durationSpanSummary(ep)
expect(identical(names(s), c("median_span_months", "span_months_q25",
                             "span_months_q75", "median_elapsed_months")),
       "returns the four stamp fields, in order")
expect(s$median_span_months == 4 && s$span_months_q25 == 1 && s$span_months_q75 == 4,
       sprintf("span = n_months - 1: median 4, IQR 1-4 (got %s, %s-%s)",
               s$median_span_months, s$span_months_q25, s$span_months_q75))
expect(isTRUE(all.equal(s$median_elapsed_months, 120 / 30.44)),
       sprintf("median elapsed months = median(TR - TL) / 30.44 (got %.4f)", s$median_elapsed_months))
expect(all(vapply(s, function(x) is.numeric(x) && length(x) == 1L, logical(1))),
       "every field is a numeric scalar (report_cascade() keeps those)")
err <- tryCatch({ durationSpanSummary(ep[, c("n_months", "TL")]); NULL },
                error = function(e) conditionMessage(e))
expect(!is.null(err) && grepl("TR", err), "a missing column is refused, naming it")

# --- [2] check_occ_floor_offset ---------------------------------------------
cat("\n[2] check_occ_floor_offset()\n")
cfg4 <- lag_selection_cfg(list(rule = "best_elpd", occ_floor_offset = 4))
stamp <- list(median_span_months = 4, span_months_q25 = 3, span_months_q75 = 7)

r <- conds_of(check_occ_floor_offset(cfg4, stamp))
expect(identical(r$value, "match") && !length(r$warnings) && !length(r$messages),
       "offset 4 vs stamped span 4: silent, 'match'")
r <- conds_of(check_occ_floor_offset(cfg4, utils::modifyList(stamp, list(median_span_months = 5))))
expect(identical(r$value, "mismatch") && length(r$warnings) == 1L && !length(r$messages),
       "offset 4 vs stamped span 5: one warning, 'mismatch'")
expect(length(r$warnings) == 1L && grepl("median_span_months", r$warnings) && grepl("3-7", r$warnings),
       "the warning names the stamp field and the IQR")
r <- conds_of(check_occ_floor_offset(cfg4, list(n_episodes = 10L, median_duration_d = 5)))
expect(identical(r$value, "missing") && !length(r$warnings) && length(r$messages) == 1L,
       "a stamp without median_span_months: one message, no warning, 'missing'")
r <- conds_of(check_occ_floor_offset(cfg4, NULL))
expect(identical(r$value, "missing") && !length(r$warnings) && length(r$messages) == 1L,
       "no stamp at all: one message, 'missing'")
r <- conds_of(check_occ_floor_offset(lag_selection_cfg(NULL), list(median_span_months = 5)))
expect(identical(r$value, "no_floor") && !length(r$warnings) && !length(r$messages),
       "no occ_floor_offset (the default 'none'): nothing to check, silent")
r <- conds_of(check_occ_floor_offset(lag_selection_cfg(list(occ_floor_offset = 4L)), stamp))
expect(identical(r$value, "match") && !length(r$warnings),
       "an integer offset matches a double span")

# --- [3] outbreakEpisodeGuardFlags vs guardOutbreakEpisodes -----------------
cat("\n[3] outbreakEpisodeGuardFlags() vs guardOutbreakEpisodes()\n")
eps <- tibble::tibble(
  location_period_id = sprintf("lp%02d", 1:12),
  suspected_cases    = c(-3, 0, NA, 1, 49, 50, 1000, 49, 0, 200, 300, 60),
  TL                 = d0,
  TR                 = d0 + c(30, 30, 30, 30, 30, 30, 730, 731, 900, 731, NA, 0))
fl <- outbreakEpisodeGuardFlags(eps)
expect(nrow(fl) == nrow(eps) && identical(names(fl), c("cases_le0", "cases_lt50", "gt730d", "any")),
       "one row per episode; columns cases_le0, cases_lt50, gt730d, any")
kept <- suppressWarnings(suppressMessages(guardOutbreakEpisodes(eps)))
dropped <- setdiff(eps$location_period_id, kept$location_period_id)
expect(identical(sort(eps$location_period_id[fl$any]), sort(dropped)),
       sprintf("`any` selects exactly the rows the guard drops (%d: %s)",
               length(dropped), paste(dropped, collapse = ", ")))
expect(identical(which(fl$cases_le0), c(1L, 2L, 9L)), "cases_le0: suspected_cases <= 0 (NA never flags)")
expect(identical(which(fl$cases_lt50), c(4L, 5L, 8L)), "cases_lt50: 0 < suspected_cases < 50, disjoint from cases_le0")
expect(identical(which(fl$gt730d), c(8L, 9L, 10L)), "gt730d: TR - TL > 730 days (730 itself and NA TR do not flag)")
# The guard's own warning carries the per-class counts: they must agree.
w <- tryCatch({ guardOutbreakEpisodes(eps); NA_character_ }, warning = function(c) conditionMessage(c))
expect(grepl(sprintf("%d suspected_cases<=0, %d suspected_cases<50, %d duration>730d",
                     sum(fl$cases_le0), sum(fl$cases_lt50), sum(fl$gt730d)), w, fixed = TRUE),
       "per-class counts equal the ones guardOutbreakEpisodes() reports")
fl2 <- outbreakEpisodeGuardFlags(eps, min_cases = 10, max_duration_days = 100)
kept2 <- suppressWarnings(suppressMessages(guardOutbreakEpisodes(eps, min_cases = 10, max_duration_days = 100)))
expect(identical(sort(eps$location_period_id[fl2$any]),
                 sort(setdiff(eps$location_period_id, kept2$location_period_id))),
       "the same rows under non-default thresholds")

# --- [4] no caller of guardOutbreakEpisodes() --------------------------------
cat("\n[4] guardOutbreakEpisodes() has no caller in analysis/ outside scratchpads\n")
r_files <- list.files("analysis", pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!startsWith(r_files, "analysis/scratchpads/")]
expect(length(r_files) > 100L, sprintf("scanned %d R files outside scratchpads", length(r_files)))
uses <- callers_of("guardOutbreakEpisodes", r_files)
self <- "analysis/tests/test_duration_span_stamp.R"
defn <- uses$kind == "definition" & uses$file == "analysis/utils/occ_data.R"
expect(sum(defn) == 1L, "found its one definition in utils/occ_data.R")
expect(!any(uses$kind == "unparseable"),
       sprintf("every file outside scratchpads parses (unparseable: %s)",
               paste(uses$file[uses$kind == "unparseable"], collapse = ", ")))
offenders <- uses[!defn & uses$file != self & uses$kind != "unparseable", , drop = FALSE]
expect(nrow(offenders) == 0L,
       sprintf("no call, reference or quoted name elsewhere (offenders: %s)",
               if (nrow(offenders)) paste(sprintf("%s:%d", offenders$file, offenders$line), collapse = ", ")
               else "none"))
# callers_of() itself sees a call when there is one (positive control).
expect(any(uses$file == self & uses$kind == "call"), "callers_of() finds this test's own calls")

# --- [5] 03_prepare_dur.R --restamp on a synthetic panel ---------------------
cat("\n[5] 03_prepare_dur.R --restamp (synthetic panel, temp dir)\n")
source("analysis/utils/provenance.R")   # prepPolicyStamp() & co. (the script's own writers)
work <- file.path(tempdir(), "dur_restamp")
unlink(work, recursive = TRUE)
dir.create(file.path(work, "generated_data"), recursive = TRUE)
invisible(file.symlink(normalizePath("analysis"), file.path(work, "analysis")))
size_path  <- "generated_data/outbreak_size_analysis_data.rds"          # makeSizeDatFile(lp_season)
panel_path <- "generated_data/outbreaks_duration_data_base-1990_2019.rds" # makeDurationDatFile()
saveRDS(eps, file.path(work, size_path))
panel <- tibble::tibble(
  episode_id      = c("a", "a", "b", "b", "b", "c"),
  country         = "KEN",
  duration_months = c(2L, 2L, 3L, 3L, 3L, 1L),
  censored        = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  ended           = c(0L, 1L, 0L, 0L, 1L, 0L),
  TL              = d0 + c(0, 0, 10, 10, 10, 5),
  TR              = d0 + c(40, 40, 70, 70, 70, 9))
saveRDS(panel, file.path(work, panel_path))
old <- prepPolicyStamp(pop_source = "crosswalk", composite_policy = "drop_collapse",
                       cover_threshold = 0.9, counts = list(
                         track = "duration", baseline = "1990_2019", source_episodes = size_path,
                         n_episodes = 3L, n_censored = 1L, n_rows_out = 6L, n_events_out = 2L,
                         median_duration_d = 2))
persistPrepPolicyStamp(file.path(work, panel_path), old)
panel_md5 <- unname(tools::md5sum(file.path(work, panel_path)))

run_restamp <- function() {
  log <- file.path(work, "restamp.log")
  owd <- setwd(work); on.exit(setwd(owd))
  st <- system2("Rscript", c("analysis/03_prepare_dur.R", "-b", "1990_2019", "--restamp"),
                stdout = log, stderr = log)
  list(status = st, log = readLines(log, warn = FALSE))
}
r <- run_restamp()
expect(r$status == 0L, sprintf("--restamp exits 0 (log tail: %s)", paste(tail(r$log, 3), collapse = " | ")))
new <- readRDS(file.path(work, sub("\\.rds$", ".policy.rds", panel_path)))
want_span <- durationSpanSummary(tibble::tibble(n_months = c(2L, 3L, 1L), TL = d0 + c(0, 10, 5),
                                                TR = d0 + c(40, 70, 9)))
expect(identical(new[names(want_span)], want_span),
       sprintf("span fields from the panel (median span %s)", want_span$median_span_months))
want_guard <- list(n_guard_cases_le0 = sum(fl$cases_le0), n_guard_cases_lt50 = sum(fl$cases_lt50),
                   n_guard_gt730d = sum(fl$gt730d), n_guard_any = sum(fl$any))
expect(identical(new[names(want_guard)], want_guard),
       sprintf("guard counts from the episode table (any = %d)", want_guard$n_guard_any))
expect(identical(new[names(old)], old), "every existing field kept, stamped_at included")
expect(is.character(new$restamped_at) && nzchar(new$restamped_at), "restamped_at recorded")
expect(identical(unname(tools::md5sum(file.path(work, panel_path))), panel_md5),
       "the panel itself is untouched")

# A stamp that does not describe the panel is refused, and left as it was.
bad <- utils::modifyList(old, list(n_rows_out = 7L))
persistPrepPolicyStamp(file.path(work, panel_path), bad)
r <- run_restamp()
expect(r$status != 0L && any(grepl("n_rows_out", r$log)),
       "a stamp naming another row count is refused, naming the field")
expect(identical(readRDS(file.path(work, sub("\\.rds$", ".policy.rds", panel_path))), bad),
       "a refused restamp writes nothing")
unlink(work, recursive = TRUE)

# --- [6] deploy_select_lags() calls check_occ_floor_offset() ---------------
cat("\n[6] deploy_select_lags() checks occ_floor_offset against the stamped span\n")
source("analysis/utils/deploy_core.R")
# Mocks: one group, rankings over lags -2..6 (occ best at 5, dur best at 0).
read_latest_ranking <- function(dir, group) {
  best <- if (grepl("occ", dir)) 5L else 0L
  lags <- -2:6
  data.frame(single_lag = lags, rank = rank(abs(lags - best), ties.method = "first"))
}
.deploy_group_full <- function(dep, group) group
dep6 <- list(fits = list(groups = list("east"), occ_dir = "b_occ", dur_dir = "b_dur", mon_dir = "b_mon"),
             episodes = list(key = "gid_season"),
             lag_selection = list(rule = "best_elpd", occ_floor_offset = 4))
work6 <- file.path(tempdir(), "dur_floor_check"); unlink(work6, recursive = TRUE)
dir.create(file.path(work6, "generated_data"), recursive = TRUE)
owd6 <- setwd(work6)
panel6 <- as.character(makeDurationDatFile(list(baseline = get_default_options()$baseline,
                                                episodes = list(key = "gid_season"))))
put_stamp <- function(med) saveRDS(list(n_episodes = 9L, median_span_months = med,
                                        span_months_q25 = 3, span_months_q75 = 7),
                                   prepPolicyStampPath(panel6))
r_miss <- conds_of(deploy_select_lags(dep6))
put_stamp(5); r_bad <- conds_of(deploy_select_lags(dep6))
put_stamp(4); r_ok  <- conds_of(deploy_select_lags(dep6))
dep6n <- dep6; dep6n$lag_selection$occ_floor_offset <- "none"
put_stamp(5); r_none <- conds_of(deploy_select_lags(dep6n))
setwd(owd6); unlink(work6, recursive = TRUE)
expect(grepl("_ep-gs", panel6), sprintf("the gid_season duration panel is read (%s)", panel6))
expect(!length(r_miss$warnings) && sum(grepl("median_span_months", r_miss$messages)) == 1L,
       "no duration stamp: one message, no warning")
expect(length(r_bad$warnings) == 1L && grepl("occ_floor_offset = 4", r_bad$warnings) &&
         grepl("span 5", r_bad$warnings, fixed = TRUE),
       "stamped span 5 vs offset 4: one warning naming both")
expect(!length(r_ok$warnings) && !length(r_ok$messages), "stamped span 4 vs offset 4: silent")
expect(!length(r_none$warnings) && !length(r_none$messages), "occ_floor_offset none: nothing checked, silent")
expect(identical(r_miss$value, r_bad$value) && identical(r_bad$value, r_ok$value) &&
         r_ok$value$east[["occ"]] == 5L && r_ok$value$east[["dur"]] == 1L,
       sprintf("the selection does not depend on the check (occ %d, dur %d floored at 5 - 4)",
               r_ok$value$east[["occ"]], r_ok$value$east[["dur"]]))

# --- [7] report_cascade() carries the span IQR --------------------------------
cat("\n[7] report_cascade() keeps the span fields\n")
renv7 <- new.env(); sys.source("analysis/utils/results_derive.R", renv7)
work7 <- tempfile("cascade_"); dir.create(work7)
dur7 <- file.path(work7, "dur.rds")
saveRDS(list(pop_source = "legacy", n_episodes = 9L, median_duration_d = 5,
             median_span_months = 4, span_months_q25 = 3, span_months_q75 = 7,
             median_elapsed_months = 4.57, n_guard_any = 5L), prepPolicyStampPath(dur7))
ctx7 <- list(PREP_COUNTS_DIR = work7, SIZE_DATA_PATH = file.path(work7, "none_siz.rds"),
             SURV_PATH = file.path(work7, "none_occ.rds"), DUR_DATA_PATH = dur7,
             MON_DATA_PATH = file.path(work7, "none_mon.rds"))
c7 <- suppressWarnings(renv7$report_cascade(ctx7))
m7 <- c7$metric[c7$step == "03_prepare_dur"]
expect(all(c("median_span_months", "span_months_q25", "span_months_q75", "median_elapsed_months",
             "n_guard_any") %in% m7) && !"pop_source" %in% m7,
       sprintf("03_prepare_dur rows: %s", paste(m7, collapse = ", ")))
expect(c7$value[c7$step == "03_prepare_dur" & c7$metric == "span_months_q75"] == 7,
       "span_months_q75 reaches the cascade with its value")
unlink(work7, recursive = TRUE)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
