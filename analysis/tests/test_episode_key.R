# test_episode_key.R -- the cross-track canonical episode definition.
#
# HARD gate (see analysis/tests/run_all.R).
#
# WHAT THIS GUARDS
# ----------------
# `episodes.key: gid_season` makes ONE outbreak episode per (gid, year_seas)
# the shared object of all four tracks (utils/footprint.R::
# buildCanonicalEpisodes()): onset = earliest member TL, end = latest member TR,
# cases = the finest-resolution member total, admin_level = the gid's
# geometry-derived level. Every prep artefact is stamped with a digest of that
# table and every loader / consumer checks it (utils/provenance.R). The default
# `lp_season` must keep every historical behaviour and artefact name.
#
# Pure synthetic tables. No sf, no raster, no SQLite, no network, no data/.
# No MCMC.
#
# Run: Rscript analysis/tests/test_episode_key.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lubridate)
  library(purrr)
})

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/footprint.R")
source("analysis/utils/composite_lp.R")   # transitive dependency of population.R
source("analysis/utils/population.R")
source("analysis/utils/occ_data.R")

# .aggregate_footprint_counts() is extracted by PARSING spatpomp_core.R rather
# than sourcing it, so this test does not need spatPomp / pomp installed.
local({
  ex <- parse("analysis/utils/spatpomp_core.R", keep.source = FALSE)
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(as.character(e[[2]]), ".aggregate_footprint_counts"),
                as.list(ex))
  stopifnot(length(hit) == 1L)
  eval(hit[[1]], envir = globalenv())
})

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

level_of <- function(gid) ifelse(grepl("^AAA_1_", gid), "ADM2", "ADM1")

# ---------------------------------------------------------------------------
# Fixture: gid AAA_1_2_1 (canonical ADM2), season 2010 (min_seas = 1).
#   LP 101  reported at ADM2   2010-03-05 -> 2010-05-20   100 cases
#   LP 102  reported at ADM3   2010-04-10 -> 2010-07-01    60 cases
#   LP 103  reported at ADM3   2010-03-20 -> 2010-04-30    50 cases
# Plus a second gid BBB_2_1 (canonical ADM1) with one LP in 2011.
# ---------------------------------------------------------------------------
mk <- function(lp, gid, lvl, tl, tr, cases, ms = 1L, pop = 1000) {
  tibble(uid = paste0(lp, "-1"), location_period_id = lp, gid = gid,
         country = substr(gid, 1, 3), location = paste0("loc", lp),
         admin_level = lvl, TL = as.Date(tl), TR = as.Date(tr),
         suspected_cases = cases, pop = pop, min_seas = ms)
}
rows <- bind_rows(
  mk(101, "AAA_1_2_1", "ADM2", "2010-03-05", "2010-05-20", 100),
  mk(102, "AAA_1_2_1", "ADM3", "2010-04-10", "2010-07-01",  60, pop = 300),
  mk(103, "AAA_1_2_1", "ADM3", "2010-03-20", "2010-04-30",  50, pop = 200),
  mk(201, "BBB_2_1",   "ADM1", "2011-02-01", "2011-03-15",  40))
cal_in <- rows %>%
  mutate(TL_seas = .shiftMonths(TL, -(min_seas - 1L)),
         TR_seas = .shiftMonths(TR, -(min_seas - 1L)),
         year_seas = year(TL_seas), month_seas = month(TL_seas))

cat("# canonical gid-season episodes\n")

# --- 1. the canonical episode -----------------------------------------------
cat("\n[1] one episode per (gid, year_seas): min TL, max TR, finest-level total\n")
res <- buildCanonicalEpisodes(cal_in, level_of)
ep  <- res$episodes %>% arrange(gid)
a   <- ep %>% filter(gid == "AAA_1_2_1")
expect(nrow(ep) == 2L, "two gids -> two episodes")
expect(a$TL == as.Date("2010-03-05") && a$TR == as.Date("2010-07-01"),
       "onset = earliest member TL, end = latest member TR")
expect(a$suspected_cases == 110, "cases = finest-level total (ADM3: 60 + 50), not 210")
expect(a$pop == 500, "pop = the same counted members' population (300 + 200)")
expect(a$admin_level == "ADM2", "admin_level is the GEOMETRY level of the gid")
expect(a$admin_level_corpus == "ADM3", "admin_level_corpus records the counted reporting level")
expect(as.character(a$censor_lp_id) == "102", "censoring member = the one whose report ends last")
expect(a$n_members == 3L && a$n_levels == 2L, "membership counts")
m <- res$members %>% filter(gid == "AAA_1_2_1") %>% arrange(location_period_id)
expect(identical(m$in_total, c(FALSE, TRUE, TRUE)), "members flag exactly the counted ones")
expect(all(m$admin_level == "ADM2") && identical(m$report_level, c("ADM2", "ADM3", "ADM3")),
       "members carry the episode level AND their own reporting level")
expect(a$uid == "101-1" && a$location_period_id == 101,
       "representative uid/LP per gid = smallest location_period_id")
expect(res$report$n_episodes_canonical == 2L && res$report$n_multi_level_episodes == 1L,
       "report counts")

cat("\n[1b] a tie between levels goes to the FINER level\n")
tie <- cal_in %>% mutate(suspected_cases = ifelse(location_period_id == 101, 110, suspected_cases))
rt  <- buildCanonicalEpisodes(tie, level_of)$episodes %>% filter(gid == "AAA_1_2_1")
expect(rt$admin_level_corpus == "ADM3" && rt$suspected_cases == 110, "tie -> ADM3")

cat("\n[1c] a gid absent from the canonical table is an error, never a default\n")
lk <- function(g) { v <- level_of(g); v[g == "BBB_2_1"] <- NA; if (anyNA(v)) stop("absent") ; v }
expect(throws(buildCanonicalEpisodes(cal_in, lk)), "missing canonical level stops")

# --- 2. parity with the simulator's observation rule ------------------------
cat("\n[2] finestLevelTotals() == spatpomp_core.R::.aggregate_footprint_counts()\n")
fl <- finestLevelTotals(cal_in, "suspected_cases", c("gid", "year_seas"))
ag <- .aggregate_footprint_counts(cal_in, "suspected_cases", "v", keys = c("gid", "year_seas"))
j  <- inner_join(fl, ag, by = c("gid", "year_seas"))
expect(nrow(j) == 2L && all(j$total == j$v), "same totals on every key")

# --- 3. one season calendar per gid ------------------------------------------
cat("\n[3] canonicalSeasonCalendar(): mode of min_seas per gid, re-dated\n")
cc <- bind_rows(
  mk(301, "CCC_1_1", "ADM1", "2012-02-10", "2012-03-01", 5, ms = 3L),
  mk(302, "CCC_1_1", "ADM1", "2012-05-10", "2012-06-01", 5, ms = 3L),
  mk(303, "CCC_1_1", "ADM1", "2012-02-20", "2012-03-20", 5, ms = 1L)) %>%
  mutate(TL_seas = .shiftMonths(TL, -(min_seas - 1L)),
         TR_seas = .shiftMonths(TR, -(min_seas - 1L)),
         year_seas = year(TL_seas))
ccal <- canonicalSeasonCalendar(cc)
expect(all(ccal$min_seas == 3L), "mode wins (3 x2 vs 1 x1)")
expect(attr(ccal, "n_min_seas_reassigned") == 1L, "one row reassigned")
expect(ccal$year_seas[ccal$location_period_id == 303] == 2011L,
       "a February onset on a March calendar falls in the previous season-year")
tie2 <- cc %>% filter(location_period_id != 302)
expect(all(canonicalSeasonCalendar(tie2)$min_seas == 1L), "ties -> smallest month")
dup <- bind_rows(mk(401, "DDD_1_1", "ADM1", "2013-01-05", "2013-01-20", 1, ms = 1L),
                 mk(401, "DDD_1_1", "ADM1", "2013-02-05", "2013-02-20", 1, ms = 1L)) %>%
  mutate(uid = c("401-1", "401-2"), TL_seas = TL, TR_seas = TR, year_seas = year(TL))
dcal <- canonicalSeasonCalendar(dup)
expect(nrow(dcal) == 1L && attr(dcal, "n_redate_dropped") == 1L &&
         dcal$TL == as.Date("2013-01-05"),
       "first outbreak per (LP, gid, season) re-applied after re-dating")

# --- 4. the gid-season occurrence builder ----------------------------------
cat("\n[4] gidAtRiskYears() / gidSeasonSurvInput() / makeSurvData()\n")
all_ts_dat <- tibble(
  location_period_id = c(101, 102, 103, 201),
  TL = as.Date(c("2009-01-01", "2009-06-01", "2010-01-01", "2011-01-01")),
  TR = as.Date(c("2010-12-31", "2012-06-30", "2010-06-30", "2011-12-31")))
lp_gid <- res$member_rows %>% distinct(location_period_id, gid, min_seas)
ar <- gidAtRiskYears(lp_gid, all_ts_dat) %>% arrange(gid, year_seas)
arA <- ar %>% filter(gid == "AAA_1_2_1")
expect(identical(as.integer(arA$year_seas), 2009:2012), "union of member coverage years")
expect(arA$TR_seas[arA$year_seas == 2012] == as.Date("2012-06-30"),
       "a no-episode year is at risk until the LATEST member coverage end")
si <- gidSeasonSurvInput(res$episodes, ar)
expect(sum(si$.has_onset) == 2L && nrow(si) == nrow(distinct(si, gid, year_seas)),
       "one input row per (gid, season); onsets = episodes")
# makeSurvData() with stubbed climate lookups.
getPrecipVars <- function() "precip"
getTempVars   <- function() "temp"
getStatsGid   <- function(gid) {
  mo <- seq(as.Date("2008-01-01"), as.Date("2013-12-01"), by = "month")
  tibble(gid = gid, month = mo, precip = 1, temp = 20)
}
getExtremesGid <- function(gid, baseline) {
  mo <- seq(as.Date("2008-01-01"), as.Date("2013-12-01"), by = "month")
  bind_rows(tibble(gid = gid, month = mo, name = "precip", mean = 1, sd = 1, hilo_q = 0L, hilo_sd = 0L),
            tibble(gid = gid, month = mo, name = "temp",   mean = 20, sd = 1, hilo_q = 0L, hilo_sd = 0L))
}
sv <- si %>% group_by(gid, year_seas) %>% group_split() %>%
  map_dfr(function(x) makeSurvData(x, baseline = "1990_2019") %>%
            mutate(outbreak = outbreak & x$.has_onset[1]))
svA <- sv %>% filter(gid == "AAA_1_2_1")
expect(sum(svA$outbreak) == 1L, "exactly one onset for the 3-LP gid-season")
expect(max(svA$t_end[year(svA$t_start) == 2010]) == as.Date("2010-03-05"),
       "the gid is censored at the canonical onset -- no post-onset rows that season")
expect(nrow(sv) == nrow(distinct(sv, gid, t_start)), "one row per (gid, t_start)")

cat("\n[4b] addNoOutbreakYears() still behaves as before (lpCoverageSeasonYears extraction)\n")
x101 <- cal_in %>% filter(location_period_id == 101)
pad  <- addNoOutbreakYears(x101, all_ts_dat) %>% arrange(year_seas)
expect(identical(as.integer(pad$year_seas), 2009:2010), "padded years = covered years")
expect(pad$TR_seas[pad$year_seas == 2010] == x101$TL_seas,
       "the outbreak year's padded window ends at the outbreak's own onset")
expect(pad$TR_seas[pad$year_seas == 2009] == as.Date("2009-12-31"), "a full covered year")

# --- 5. monthly aggregation through the canonical membership ----------------
cat("\n[5] aggregateMonthlyToFootprint(key = gid_season) conserves the episode total\n")
monthly <- tibble(
  uid = c("101-1", "101-1", "102-1", "102-1", "103-1", "201-1"),
  location_period_id = c(101, 101, 102, 102, 103, 201),
  month = as.Date(c("2010-03-01", "2010-04-01", "2010-04-01", "2010-06-01",
                    "2010-03-01", "2011-02-01")),
  cases_month = c(70, 30, 40, 20, 50, 40))
mem_in <- res$members %>% filter(in_total)
ag5 <- aggregateMonthlyToFootprint(monthly, mem_in,
                                   footprints = res$episodes %>%
                                     select(gid, year_seas, suspected_cases),
                                   key = footprintKey("gid_season"))
expect(all(ag5$recon$conserved), "every canonical episode reconciles exactly")
expect(sum(ag5$monthly$cases_month[ag5$monthly$gid == "AAA_1_2_1"]) == 110,
       "the ADM2 report is not added on top of its ADM3 children")
cat("\n[5b] the default key is unchanged (gid, year_seas, admin_level)\n")
expect(identical(footprintKey(), FOOTPRINT_KEY) &&
         identical(footprintKey("lp_season"), c("gid", "year_seas", "admin_level")),
       "lp_season key = FOOTPRINT_KEY")
mem_lp <- footprintMembers(cal_in)
ag5b <- aggregateMonthlyToFootprint(monthly, mem_lp)
expect(nrow(ag5b$monthly %>% distinct(gid, year_seas, admin_level)) == 3L,
       "default key still splits the AAA gid-season by reporting level")

# --- 6. the digest -----------------------------------------------------------
cat("\n[6] episodeDigest()\n")
d1 <- episodeDigest(res$episodes)
d2 <- episodeDigest(res$episodes[nrow(res$episodes):1, ] %>%
                      mutate(gid = sub("^AAA_1_2_1$", "AAA.1.2_1", gid)))
expect(!is.na(d1) && identical(d1, d2), "invariant to row order and dotted/flat gid spelling")
d3 <- episodeDigest(res$episodes %>% mutate(TL = TL + 1))
expect(!identical(d1, d3), "changes when an onset changes")

# --- 7. guards -----------------------------------------------------------------
cat("\n[7] assertEpisodeKey() / assertEpisodeCoherence()\n")
td <- tempfile("epk"); dir.create(td)
f_lp <- file.path(td, "lp.rds");  saveRDS(1, f_lp)
f_gs <- file.path(td, "gs.rds");  saveRDS(1, f_gs)
f_gs2 <- file.path(td, "gs2.rds"); saveRDS(1, f_gs2)
f_none <- file.path(td, "none.rds"); saveRDS(1, f_none)
persistPrepPolicyStamp(f_lp, prepPolicyStamp("crosswalk", "drop_collapse", .9, "derived"))
persistPrepPolicyStamp(f_gs, prepPolicyStamp("crosswalk", "drop_collapse", .9, "derived",
                                             episode_key = "gid_season", episode_digest = d1))
persistPrepPolicyStamp(f_gs2, prepPolicyStamp("legacy", "drop_collapse", .9, "derived",
                                              episode_key = "gid_season", episode_digest = d3))
opt_lp <- list(baseline = "1990_2019")
opt_gs <- list(baseline = "1990_2019", episodes = list(key = "gid_season"))
expect(!throws(assertEpisodeKey(opt_lp, f_none, "occurrence")),
       "default config + unstamped (historical) artefact passes")
expect(!throws(assertEpisodeKey(opt_lp, f_lp, "size")), "lp config + lp artefact passes")
expect(throws(assertEpisodeKey(opt_lp, f_gs, "size"), "Re-run"),
       "lp config + gid_season artefact stops with the re-run command")
expect(throws(assertEpisodeKey(opt_gs, f_lp, "occurrence"), "03_prepare_occ.R"),
       "gid_season config + lp artefact stops naming the prep script")
expect(identical(assertEpisodeKey(opt_gs, f_gs, "monthly")$digest, d1),
       "gid_season config + gid_season artefact returns the record")
expect(!throws(assertEpisodeCoherence(list(a = episodeStampOf(f_lp), b = episodeStampOf(f_none)))),
       "all lp_season -> silent")
expect(throws(assertEpisodeCoherence(list(a = episodeStampOf(f_lp), b = episodeStampOf(f_gs)))),
       "a lp/gid mix stops")
expect(throws(assertEpisodeCoherence(list(occ = episodeStampOf(f_gs), siz = episodeStampOf(f_gs2))),
              "pop_source=legacy"),
       "digest mismatch stops and names each side's pop_source")
expect(throws(assertEpisodeCoherence(list(surv = episodeStampOf(f_gs),
                                          fit = fitEpisodeStampOf(structure(list(), class = "x"))))),
       "a gid_season artefact paired with an unstamped (lp_season) fit stops")
fit_ok <- structure(list(), episode_stamp = list(key = "gid_season", digest = d1))
expect(!throws(assertEpisodeCoherence(list(surv = episodeStampOf(f_gs), fit = fitEpisodeStampOf(fit_ok)))),
       "matching fit record passes")

# --- 8. paths ---------------------------------------------------------------
cat("\n[8] artefact names: unchanged by default, _ep-gs under gid_season\n")
expect(makeSurvDatFile(opt_lp) == "generated_data/outbreaks_survival_data_base-1990_2019.rds" &&
         makeSizeDatFile() == "generated_data/outbreak_size_analysis_data.rds" &&
         makeFootprintMembersFile() == "generated_data/outbreak_footprint_members.rds" &&
         makeDurationDatFile(opt_lp) == "generated_data/outbreaks_duration_data_base-1990_2019.rds" &&
         makeMonthlyDatFile(opt_lp) == "generated_data/outbreaks_monthly_data_base-1990_2019.rds",
       "default prep paths are the historical literals")
expect(makeSurvDatFile(opt_gs) == "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs.rds" &&
         makeSizeDatFile(opt_gs) == "generated_data/outbreak_size_analysis_data_ep-gs.rds" &&
         makeFootprintMembersFile(opt_gs) == "generated_data/outbreak_footprint_members_ep-gs.rds",
       "gid_season prep paths carry _ep-gs")
expect(throws(episodeKeyOf(list(episodes = list(key = "gid"))), "Refusing"),
       "an unknown key is an error, not a fallback")

# --- 9. TZA mainland ------------------------------------------------------------
cat("\n[9] markTzaMainlandAdm0()\n")
tz <- tibble(uid = c("AFR::TZA::Mainland-1", "x-1"), admin_level = c("ADM1", "ADM2"))
expect(identical(markTzaMainlandAdm0(tz)$admin_level, c("ADM0", "ADM2")),
       "TZA mainland -> ADM0, everything else untouched")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
