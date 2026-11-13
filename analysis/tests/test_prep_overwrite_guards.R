# test_prep_overwrite_guards.R — overwrite guards of the 03 scripts that write
# fixed-name artefacts.
#
#   [1] assertPrepStampOverwrite(): the P6 defaults unchanged (no artefact; no
#       stamp; same stamp but stamped_at; a changed field; --force); the new
#       `ignore` (restamped_at of 03_prepare_dur.R --restamp) and `old_stamp`
#       arguments, on synthetic files with real prepPolicyStamp() shapes
#   [2] whoMonthlyGuardStamp(): the meta written to and read back from YAML
#       compares equal to the in-memory one; `written` / input mtimes are
#       ignored; a content or allowlist change is caught; no meta = refused
#   [3] every guarded script (siz, dur, mon, sim panel, projection, extended,
#       WHO monthly): a --force option, the stamp built and the guard run
#       BEFORE the artefact, its stamp or its meta is written
#
# Synthetic files only, no data. HARD, quick.
# Run: Rscript analysis/tests/test_prep_overwrite_guards.R

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

td <- tempfile("guards"); dir.create(td)

cat("\n[1] assertPrepStampOverwrite(): defaults and the new arguments\n")
pp <- file.path(td, "outbreak_duration_data_base-1990_2019_ep-gs.rds")
st <- prepPolicyStamp("crosswalk", "drop_collapse", 0.9, "derived", pop_digest = NA_character_,
                      episode_key = "gid_season", episode_digest = "abc",
                      counts = list(track = "duration", n_rows_out = 10L, n_episodes = 3L))
expect(identical(errmsg(assertPrepStampOverwrite(pp, st)), ""), "no artefact: nothing to guard")
saveRDS(1, pp)
expect(grepl("no stamp", errmsg(assertPrepStampOverwrite(pp, st))), "an artefact without a stamp is refused")
expect(identical(suppressMessages(assertPrepStampOverwrite(pp, st, force = TRUE)), "<no stamp>"),
       "--force overrides a missing stamp")
persistPrepPolicyStamp(pp, st)
st2 <- st; st2$stamped_at <- "2030-01-01T00:00:00"
expect(identical(errmsg(assertPrepStampOverwrite(pp, st2)), ""), "identical stamp but stamped_at: passes (default rerun)")
st_v <- st; st_v$episode_digest <- "def"; st_v$n_rows_out <- 12L
e <- errmsg(assertPrepStampOverwrite(pp, st_v))
expect(grepl("Refusing", e) && grepl("episode_digest", e) && grepl("n_rows_out", e),
       "another vintage (episode digest, row count): refused, naming the fields")
expect(identical(sort(suppressMessages(assertPrepStampOverwrite(pp, st_v, force = TRUE))),
                 c("episode_digest", "n_rows_out")), "--force overrides and reports the fields")
# 03_prepare_dur.R --restamp adds restamped_at to the stored stamp.
persistPrepPolicyStamp(pp, c(st, list(restamped_at = "2026-10-27T10:00:00")))
expect(grepl("restamped_at", errmsg(assertPrepStampOverwrite(pp, st))),
       "default ignore: restamped_at is a difference (03_prepare_occ.R's P6 rule unchanged)")
expect(identical(errmsg(assertPrepStampOverwrite(pp, st, ignore = c("stamped_at", "restamped_at"))), ""),
       "dur's ignore = stamped_at + restamped_at: a restamped panel rebuilt identically passes")
expect(identical(errmsg(assertPrepStampOverwrite(pp, st, old_stamp = st2)), ""),
       "old_stamp: compared instead of the .policy.rds")
expect(grepl("no stamp", errmsg(assertPrepStampOverwrite(pp, st, old_stamp = NULL))), "old_stamp NULL: refused")

cat("\n[2] whoMonthlyGuardStamp()\n")
mon <- data.frame(iso3 = c("KEN", "KEN", "TZA"), t_start = as.Date(c("2024-01-01", "2024-02-01", "2024-01-01")),
                  cases = c(10, 0, 5.5), n_days = c(31L, 29L, 31L), partial = c(FALSE, FALSE, TRUE),
                  source = c("dashboard_weekly", "dashboard_weekly", "sitrep_cumulative"), stringsAsFactors = FALSE)
meta_of <- function(mon, zero = character(0), written = "2026-10-01T00:00:00", mt = "2026-09-01T00:00:00",
                    cum = TRUE) list(
  generator = "analysis/03_prepare_who_monthly.R", written = written,
  inputs = list(weekly = list(path = "w.rds", mtime = mt, n_rows = 40L),
                cumulative = if (cum) list(path = "c.csv", mtime = mt, n_rows = 7L) else NULL),
  zero_iso3 = as.list(zero), max_interval_days = 45L,
  counts = list(country_months = nrow(mon), complete = sum(!mon$partial), by_source = as.list(table(mon$source))))
for (cum in c(TRUE, FALSE)) {
  out <- file.path(td, sprintf("who_dashboard_monthly_africa_%s.rds", cum)); mp <- sub("\\.rds$", "_meta.yml", out)
  m0 <- meta_of(mon, cum = cum)
  saveRDS(mon, out); yaml::write_yaml(m0, mp)
  old <- whoMonthlyGuardStamp(yaml::read_yaml(mp), readRDS(out))
  expect(identical(old, whoMonthlyGuardStamp(meta_of(mon, written = "2026-10-29T12:00:00", mt = "2026-10-29T11:00:00",
                                                     cum = cum), mon)),
         sprintf("cumulative %s: on-disk meta == in-memory meta once written / mtimes are dropped", cum))
  expect(identical(errmsg(assertPrepStampOverwrite(out, whoMonthlyGuardStamp(meta_of(mon, cum = cum), mon),
                                                   old_stamp = old)), ""),
         sprintf("cumulative %s: identical rebuild passes", cum))
}
mon2 <- mon; mon2$cases[1] <- 11
e <- errmsg(assertPrepStampOverwrite(out, whoMonthlyGuardStamp(meta_of(mon2, cum = FALSE), mon2), old_stamp = old))
expect(grepl("content_digest", e), "a refreshed feed (content change): refused, naming content_digest")
e <- errmsg(assertPrepStampOverwrite(out, whoMonthlyGuardStamp(meta_of(mon, zero = "TZA", cum = FALSE), mon),
                                     old_stamp = old))
expect(grepl("zero_iso3", e), "another zero allowlist: refused")
expect(is.null(whoMonthlyGuardStamp(NULL, mon)), "no meta: NULL (an existing table without meta is refused)")

cat("\n[3] every guarded 03 script: --force, stamp and guard before the first write\n")
first_line <- function(src, pat) { i <- grep(pat, src); if (length(i)) min(i) else NA_integer_ }
scripts <- list(
  `03_prepare_siz.R` = list(stamp = "^size_stamp <- prepPolicyStamp\\(", guard = "assertPrepStampOverwrite\\(\\.p, size_stamp",
                            writes = c("^saveRDS\\(size_dat, SIZE_DAT_PATH\\)", "saveRDS\\(fm_out", "^persistPrepPolicyStamp\\(SIZE_DAT_PATH",
                                       "^persistPopStampSidecar\\(")),
  `03_prepare_dur.R` = list(stamp = "^dur_stamp <- prepPolicyStamp\\(", guard = "^assertPrepStampOverwrite\\(DUR_DAT_PATH, dur_stamp",
                            writes = c("^saveRDS\\(dur_dat", "^persistPrepPolicyStamp\\(DUR_DAT_PATH, dur_stamp")),
  `03_prepare_mon.R` = list(stamp = "^mon_stamp <- prepPolicyStamp\\(", guard = "^assertPrepStampOverwrite\\(MON_DAT_PATH, mon_stamp",
                            writes = c("^saveRDS\\(mon_dat", "^persistPrepPolicyStamp\\(MON_DAT_PATH")),
  `03_prepare_sim_panel.R` = list(stamp = "^\\.stamp <- prepPolicyStamp\\(", guard = "^assertPrepStampOverwrite\\(out_path, \\.stamp",
                                  writes = c("^saveRDS\\(panel", "^persistPrepPolicyStamp\\(out_path")),
  `03_prepare_sim_panel_projection.R` = list(stamp = "^\\.stamp <- prepPolicyStamp\\(", guard = "^assertPrepStampOverwrite\\(out_path, \\.stamp",
                                             writes = c("^saveRDS\\(panel", "^persistPrepPolicyStamp\\(out_path")),
  `03_prepare_sim_panel_extended.R` = list(stamp = "^\\.stamp <- prepPolicyStamp\\(", guard = "^assertPrepStampOverwrite\\(out, \\.stamp",
                                           writes = c("^saveRDS\\(panel", "^persistPrepPolicyStamp\\(out")),
  `03_prepare_who_monthly.R` = list(stamp = "^meta <- list\\(", guard = "^assertPrepStampOverwrite\\(out, whoMonthlyGuardStamp",
                                    writes = c("^saveRDS\\(mon, out\\)", "^yaml::write_yaml\\(meta, meta_path\\)")))
for (s in names(scripts)) {
  src <- readLines(file.path("analysis", s)); x <- scripts[[s]]
  i_s <- first_line(src, x$stamp); i_g <- first_line(src, x$guard)
  i_w <- vapply(x$writes, function(p) first_line(src, p), integer(1))
  expect(any(grepl("prepForceOption\\(", src)) && any(grepl("force = isTRUE\\(opt\\$force\\)", src)),
         sprintf("%s: --force option threaded to the guard", s))
  expect(!is.na(i_s) && !is.na(i_g) && !anyNA(i_w) && i_s < i_g && all(i_g < i_w),
         sprintf("%s: stamp (l.%s) -> guard (l.%s) -> every write (l.%s)", s, i_s, i_g, paste(i_w, collapse = ",")))
}
# The P6 guard of 03_prepare_occ.R keeps its own option and default arguments.
src_occ <- readLines("analysis/03_prepare_occ.R")
expect(length(grep("^assertPrepStampOverwrite\\(SURV_DAT_PATH, PANEL_STAMP, force = isTRUE\\(opt\\$force\\),$", src_occ)) == 1L,
       "03_prepare_occ.R: P6 call unchanged (default ignore / old_stamp)")

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
