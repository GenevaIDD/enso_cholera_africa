# test_occ_panel_masks.R -- occurrence-panel masks, ADM1-overlap
# orphans, discarded onsets, the TTLO level "none".
#
# HARD gate (see analysis/tests/run_all.R), quick: synthetic tables only.
#
# WHAT THIS GUARDS
#   [1] maskOccPanelRows() at the default spec returns the panel IDENTICAL.
#   [2] at_risk_start = after_episode: rows of a later season inside (TL, TR] of
#       an earlier-season episode are masked, the episode's own onset row is
#       kept, an onset there is masked and counted, zero-fill rows counted.
#   [3] adm1_overlap = censor_orphan: adm1OverlapParents() classifies parent
#       episodes orphan / covered at month level; only the ADM2 children's rows
#       inside an ORPHAN span leave; no onset is ever censored.
#   [4] .overlapAdm1Set() keeps filterOverlappingADM1() /
#       dropOverlappingADM1FromFrame() identical to their pre-P4 bodies.
#   [5] makeTimeToLastOutbreakCat(no_prior_level = TRUE): "none" appended after
#       "[6,Inf]", "[0,1)" stays the reference; ttloNoPriorOnset();
#       canonicalize_term() maps time_to_last_outbreaknone; fit token "_ttn".
#   [6] discardedOnsetTable() reasons; adm0EpisodeCounts().
#   [7] report_occ_panel_exclusions table D31 columns (pinned).
#
# Run: Rscript analysis/tests/test_occ_panel_masks.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(lubridate) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "gid", "paths", "provenance", "occ_data", "footprint", "occ_post"))
  source(sprintf("analysis/utils/%s.R", m))
errmsg <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
d <- as.Date
mrow <- function(gid, m, onset = NA, al = "ADM2", nz = FALSE) {
  t0 <- d(m)
  tibble(gid = gid, t_start = t0,
         t_end = if (is.na(onset)) ceiling_date(t0, "month") - 1 else d(onset),
         outbreak = !is.na(onset), admin_level = al, country = "AAA", national_zero_fill = nz)
}

cat("\n[1] default spec: identical\n")
ep <- tibble(gid = "AAA.1.1_1", year_seas = c(2015L, 2016L), TL = d(c("2015-11-10", "2016-02-05")),
             TR = d(c("2016-02-20", "2016-03-31")), min_seas = 1L, country = "AAA", admin_level = "ADM2")
surv <- bind_rows(
  mrow("AAA.1.1_1", sprintf("2015-%02d-01", 1:10)),
  mrow("AAA.1.1_1", "2015-11-01", onset = "2015-11-10"),
  mrow("AAA.1.1_1", "2016-01-01"),
  mrow("AAA.1.1_1", "2016-02-01", onset = "2016-02-05"),
  mrow("AAA.1.1_1", "2016-12-01", nz = TRUE))
sp0 <- .occPanelSpec("gid_season")
expect(identical(maskOccPanelRows(surv, ep, sp0), surv), "season / drop: panel returned identical, no attributes")

cat("\n[2] at_risk_start = after_episode\n")
spA <- .occPanelSpec("gid_season", at_risk_start = "after_episode")
mA <- maskOccPanelRows(surv, ep, spA)
expect(!any(mA$t_start %in% d(c("2016-01-01", "2016-02-01"))) &&
         any(mA$t_start == d("2015-11-01") & mA$outbreak),
       "Jan / Feb 2016 (inside the 2015 episode, a later season) masked; the 2015 onset row kept")
expect(attr(mA, "n_rows_masked_running_episode") == 2L && attr(mA, "n_onsets_masked_running_episode") == 1L,
       "2 rows masked, 1 of them an onset (the 2016 episode starts inside the running 2015 one)")
expect(any(mA$t_start == d("2016-12-01")) && attr(mA, "n_nz_rows_masked") == 0L,
       "a zero-fill row after the episode ended is kept")
surv_nz <- bind_rows(surv, mrow("AAA.1.1_1", "2016-03-01", nz = TRUE))
expect(attr(maskOccPanelRows(surv_nz, ep, spA), "n_nz_rows_masked") == 0L,
       "a same-season zero-fill row is not the running-episode mask's (the fill excludes it)")
ep_run <- tibble(gid = "AAA.1.1_1", year_seas = 2015L, TL = d("2015-11-10"), TR = d("2016-03-31"), min_seas = 1L)
surv_nz2 <- bind_rows(surv[1:11, ], mrow("AAA.1.1_1", "2016-03-01", nz = TRUE))
mz <- maskOccPanelRows(surv_nz2, ep_run, spA)
expect(attr(mz, "n_nz_rows_masked") == 1L && !any(mz$t_start == d("2016-03-01")),
       "a zero-fill row inside the previous season's running episode is masked and counted")
mr <- attr(mA, "masked_rows")
expect(identical(names(mr), c("gid", "t_start", "outbreak", "national_zero_fill", "reason")) &&
         all(mr$reason == "running_episode"), "masked_rows sidecar rows")

cat("\n[3] adm1_overlap = censor_orphan\n")
par <- "AAA.2_1"; c1 <- "AAA.2.1_1"; c2 <- "AAA.2.2_1"; solo <- "AAA.3_1"
ep3 <- tibble(gid = c(par, par, c1, solo), year_seas = c(2016L, 2017L, 2017L, 2016L),
              TL = d(c("2016-05-01", "2017-03-01", "2017-04-10", "2016-05-01")),
              TR = d(c("2016-06-30", "2017-04-30", "2017-05-31", "2016-05-31")), min_seas = 1L)
surv3 <- bind_rows(
  mrow(par, "2016-04-01", al = "ADM1"), mrow(par, "2016-05-01", onset = "2016-05-01", al = "ADM1"),
  mrow(c1, sprintf("2016-%02d-01", 4:7)), mrow(c2, sprintf("2016-%02d-01", 4:7)),
  mrow(c1, "2017-03-01"), mrow(c1, "2017-04-01", onset = "2017-04-10"),
  mrow(solo, "2016-05-01", onset = "2016-05-01", al = "ADM1"))
ov <- adm1OverlapParents(surv3, ep3)
expect(nrow(ov) == 2L && identical(ov$class, c("orphan", "covered")) && all(ov$gid == par),
       "parent episodes: 2016 orphan, 2017 covered (child c1's episode overlaps); the ADM1-only unit is no parent")
spC <- .occPanelSpec("gid_season", adm1_overlap = "censor_orphan")
mC <- maskOccPanelRows(surv3, ep3, spC, overlap = ov)
gone <- anti_join(surv3, mC, by = c("gid", "t_start"))
expect(setequal(paste(gone$gid, gone$t_start),
                paste(rep(c(c1, c2), each = 2), d(c("2016-05-01", "2016-06-01")))) &&
         attr(mC, "n_child_rows_censored_adm1") == 4L,
       "the children's May-Jun 2016 rows (the orphan span) are censored; the parent and 2017 rows stay")
expect(!any(gone$outbreak), "no onset censored")
bad_ov <- mutate(ov, class = "orphan")
expect(grepl("onset row", errmsg(maskOccPanelRows(surv3, ep3, spC, overlap = bad_ov))),
       "an onset inside an 'orphan' span is a hard error")
expect(grepl("needs adm1OverlapParents", errmsg(maskOccPanelRows(surv3, ep3, spC))),
       "censor_orphan without the classification refused")

cat("\n[4] .overlapAdm1Set() refactor\n")
old_filter <- function(dat) {
  adm1_with_adm2 <- dat %>% dplyr::filter(admin_level == "ADM2", !is.na(ADM1)) %>% dplyr::pull(ADM1) %>% unique()
  dplyr::filter(dat, !(admin_level == "ADM1" & ADM1 %in% adm1_with_adm2))
}
dat4 <- tibble(admin_level = factor(c("ADM1", "ADM2", "ADM1", "ADM2", "ADM1")),
               ADM1 = factor(c("A_1", "A_1", "A_2", NA, "A_3")), gid = letters[1:5])
expect(identical(filterOverlappingADM1(dat4), old_filter(dat4)) &&
         identical(filterOverlappingADM1(mutate(dat4, across(everything(), as.character))),
                   old_filter(mutate(dat4, across(everything(), as.character)))),
       "filterOverlappingADM1() identical to its pre-P4 body (factor and character columns)")
expect(identical(nrow(dropOverlappingADM1FromFrame(tibble(gid = letters[1:5]), dat4)), 4L),
       "dropOverlappingADM1FromFrame() drops the one overlapping ADM1 gid")

cat("\n[5] TTLO level none\n")
tt <- tibble(t_start = d(c("2010-01-01", "2010-02-01", "2012-03-01", "2012-04-01")),
             time_to_last_outbreak = c(NA, as.numeric(d("2010-02-01") - .TTLO_ORIGIN), 400, 40))
c0 <- makeTimeToLastOutbreakCat(tt)
c1n <- makeTimeToLastOutbreakCat(tt, no_prior_level = TRUE)
expect(identical(c0, makeTimeToLastOutbreakCat(tt, no_prior_level = FALSE)), "flag off: unchanged")
expect(identical(levels(c1n$time_to_last_outbreak), c(.TTLO_LEVELS, "none")) &&
         levels(c1n$time_to_last_outbreak)[1] == "[0,1)",
       "levels: the seven bins, then \"none\" (reference [0,1) unchanged)")
expect(identical(as.character(c1n$time_to_last_outbreak), c("none", "none", "[1,2)", "[0,1)")),
       "NA and origin rows -> none; other rows keep their bin")
expect(identical(ttloNoPriorOnset(tt$t_start, tt$time_to_last_outbreak), c(TRUE, TRUE, FALSE, FALSE)),
       "ttloNoPriorOnset(): is.na(ttlo) | t_start - ttlo == .TTLO_ORIGIN")
expect(identical(canonicalize_term(c("time_to_last_outbreaknone", "time_to_last_outbreak12")),
                 c("ttlo_none", "ttlo_1_2")), "canonicalize_term(): time_to_last_outbreaknone -> ttlo_none")
o <- list(baseline = "1990_2019", occurrence = list(backend = "brms"),
          episodes = list(key = "gid_season", national_zero_fill = TRUE, national_zero_threshold = 60L))
f0 <- basename(makeOccRunFile(o, "fit"))
o$occurrence$ttlo_no_prior_level <- TRUE
f1 <- basename(makeOccRunFile(o, "fit"))
expect(!grepl("_ttn", f0) && identical(sub("_ttn_ep-gs", "_ep-gs", f1), f0) && grepl("_ttn_ep-gs_nz60\\.rds$", f1),
       "fit token _ttn only with the flag, just before the panel token")

cat("\n[6] discarded onsets, ADM0\n")
ob2 <- tibble(location_period_id = c("L1", "L1", "L1", "L2"), gid = c("g1", "g1", "g1", "g2"),
              country = "AAA", TL = d(c("2015-02-01", "2015-08-01", "2015-03-01", "2015-02-01")),
              TR = d(c("2015-04-30", "2015-09-30", "2015-03-20", "2015-03-31")),
              suspected_cases = c(100, 80, 60, 70), year_seas = 2015L, min_seas = 1L)
ob3 <- ob2[c(1, 4), ]
mem <- tibble(gid = c("g1", "g2"), year_seas = 2015L, uid = c("u1", "u2"), location_period_id = c("L1", "L2"),
              TL = ob3$TL, TR = ob3$TR, suspected_cases = c(100, 70))
er <- list(episodes = tibble(gid = c("g1", "g2"), min_seas = 1L, country = "AAA", year_seas = 2015L),
           members = mem, member_rows = ob3)
dt <- discardedOnsetTable(ob2, ob3, ob3, er)
expect(nrow(dt) == 2L && all(dt$reason == "lp_second_onset") &&
         identical(dt$after_episode_end[order(dt$TL)], c(FALSE, TRUE)),
       "two LP second onsets, one after the kept episode ended (a second wave)")
expect(sum(dt$cases) == 140, "their cases")
mem2 <- bind_rows(mem, tibble(gid = "g1", year_seas = 2015L, uid = "u3", location_period_id = "L3",
                              TL = d(c("2015-04-10")), TR = d("2015-05-30"), suspected_cases = 5),
                  tibble(gid = "g2", year_seas = 2015L, uid = "u4", location_period_id = "L4",
                         TL = d(c("2015-06-15")), TR = d("2015-07-30"), suspected_cases = 9))
er2 <- er; er2$members <- mem2
dt2 <- discardedOnsetTable(ob2, ob3, ob3, er2)
expect(setequal(dt2$reason[dt2$reason != "lp_second_onset"], c("absorbed_member", "absorbed_after_gap")),
       "absorbed members: within the episode, and after a > 31-day gap")
a0 <- adm0EpisodeCounts(tibble(country = c("AAA", "AAA", "BBB"), admin_level = c("ADM0", "ADM0", "ADM1"),
                               suspected_cases = c(10, 20, 5)))
expect(identical(a0$country, "AAA") && a0$n_adm0_episodes == 2L && a0$cases_adm0 == 30,
       "adm0EpisodeCounts(): ADM0 rows per country")

cat("\n[7] D31 table\n")
if (exists("occ_panel_exclusions_table", mode = "function") || TRUE) {
  for (m in c("config", "results_core", "national_zero")) source(sprintf("analysis/utils/%s.R", m))
  source("analysis/utils/results_derive.R")
  s37 <- occ_panel_exclusions_table(ov %>% mutate(country = "AAA"), dt2, a0,
                                    bind_rows(mr, attr(mC, "masked_rows")) %>% mutate(country = "AAA"))
  expect(identical(names(s37), OCC_PANEL_EXCLUSION_COLS) && identical(OCC_PANEL_EXCLUSION_COLS[1:4],
         c("country", "adm1_overlap_parent_episodes", "adm1_overlap_orphan", "adm1_overlap_covered")),
         "D31 columns pinned")
  a <- s37[s37$country == "AAA", ]
  expect(a$adm1_overlap_orphan == 1 && a$child_rows_censored == 4 && a$onsets_masked_running_episode == 1 &&
           a$onsets_discarded_lp_second_onset == 2 && a$members_absorbed_after_gap == 1 &&
           a$adm0_rows_dropped == 2 && a$cases_adm0_dropped == 30,
         "D31 counts per country from the sidecars")
  s37e <- occ_panel_exclusions_table(NULL, NULL, NULL, NULL)
  expect(nrow(s37e) == 0L && identical(names(s37e), OCC_PANEL_EXCLUSION_COLS), "no sidecar rows -> empty table")
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
