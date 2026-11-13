# ===========================================================================
# 01b_preprocess_outbreaks_monthly.R -- per-MONTH outbreak case counts
# ===========================================================================
#
# WHY THIS STEP EXISTS. 01_preprocess_outbreaks.R reads the sub-weekly outbreak
# time series and immediately collapses it to one row per episode
# (01_preprocess_outbreaks.R:121-126):
#
#     group_by(uid, location, location_period_id, country) %>%
#     summarise(TL = min(TL), TR = max(TR), suspected_cases = sum(sCh))
#
# That is the right frame for the size track, but it destroys the within-outbreak
# temporal resolution, and NOTHING downstream preserves it -- the duration panel
# (utils/dur_data.R:356-368) carries the episode TOTAL repeated on every active
# month, which is useless as a monthly response. So a model of "cases in month m
# of outbreak i" has no response variable to fit.
#
# This script re-reads the same raw file and aggregates to MONTH grain instead,
# leaving 01_preprocess_outbreaks.R untouched (side-by-side rule: the two are
# alternative aggregations of one source, not a replacement).
#
# THE KEY IS (uid, location_period_id, month), NOT (uid, month). On the current
# corpus six `uid`s carry two `location_period_id`s each -- 01's grouping is by
# (uid, location, location_period_id, country), so one uid observed under two
# location periods yields two episodes. ALL SIX PAIRS OVERLAP IN CALENDAR MONTHS
# (e.g. AFR::UGA::Kotido-1 runs 2010-05-03..07-18 under lp 16743 and
# 2010-07-19..08-29 under lp 6454, both active in July 2010). Keying on
# (uid, month) would not fan the join out -- a row-count gate would pass -- but
# it would hand the same month's cases to two different episodes, and even a
# uid-level conservation check would still pass. Silent double counting.
# Measured on the duration panel: 3 duplicated (uid, month), 0 duplicated
# (uid, location_period_id, month).
#
# APPORTIONMENT. Reporting intervals are sub-weekly and straddle month
# boundaries, so an interval's `sCh` is split across the months it overlaps in
# proportion to day-overlap. The weights are normalised WITHIN each interval
# (w = overlap / sum(overlap)), not against the interval's nominal length, which
# makes conservation of cases algebraic rather than approximate and makes
# .expandActiveMonths()'s inclusive "+1 day" convention irrelevant to the split.
#
# `sCh` IS ALREADY FRACTIONAL UPSTREAM. 19.1% of episode totals in
# outbreak_metadata.rds are non-integer, with denominators of 2, 3, 4 and 6 --
# the signature of multi-period aggregate reports already divided evenly across
# their constituent intervals by OutbreakExtractR. Monthly aggregation is
# therefore legitimate, but for those episodes the monthly profile is partly an
# artefact of that upstream split rather than observed epidemiology. The
# `cases_straddling_frac` column below flags where the boundary split of THIS
# script contributed; the upstream split is not separable here and is a documented
# caveat, not something this script can undo.
#
# READS ~/secure/cholera_outbreak_data/ (the only monthly-track step that
# does). Writes ~/secure/cholera_outbreak_data/outbreak_monthly_cases.rds + a
# .policy.rds sidecar.
#
# Usage:
#   Rscript analysis/01b_preprocess_outbreaks_monthly.R
#   Rscript analysis/01b_preprocess_outbreaks_monthly.R --indir <outbreak data dir>
#   (default utils/paths.R::defaultOutbreakDataDir(), the data leg's one vintage)
# ===========================================================================

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(tibble)
  library(optparse)
})

source("analysis/utils/base.R")
source("analysis/utils/gid.R")          # filterAfrica
source("analysis/utils/paths.R")        # makeMonthlyCasesFile
source("analysis/utils/provenance.R")   # prepPolicyStamp, persistPrepPolicyStamp
source("analysis/utils/dur_data.R")     # .expandActiveMonths
source("analysis/utils/occ_data.R")
source("analysis/utils/anchor.R")

option_list <- list(
  make_option("--indir",
              default = defaultOutbreakDataDir(), type = "character",   # utils/paths.R, the one data vintage
              help    = "folder with time_series_outbreak_extraction.rds [default: %default]"),
  make_option("--metadata",
              default = makeOutbreakMetadataFile(), type = "character",
              help    = "episode table written by 01_preprocess_outbreaks.R, for gate B [default: %default]"),
  make_option("--tolerance",
              default = 1e-6, type = "double",
              help    = "relative tolerance for the case-conservation gates [default: %default]"),
  make_option("--drop_cases_le0",
              default = TRUE, action = "store", type = "logical",
              help    = "drop outbreaks whose total cases are <= 0, as 01_preprocess_outbreaks.R does [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Africa ISO3 filter. The same africaIso3() as 01_preprocess_outbreaks.R --
# gate B below compares totals against that script's output, so any divergence
# in the country filter would surface as a spurious conservation failure.
africa_iso3 <- africaIso3()   # utils/gid.R: giscoR Africa + the synthetic countries


# 1. Read the raw interval series -----------------------------------------

ts_path <- file.path(opt$indir, "time_series_outbreak_extraction.rds")
if (!file.exists(ts_path))
  stop("01b: raw outbreak time series not found: ", ts_path)

cat("---- Reading", ts_path, "\n")

raw <- readRDS(ts_path) %>%
  ungroup() %>%
  mutate(
    location_period_id = as.character(location_period_id),
    TL  = as.Date(TL),
    TR  = as.Date(TR),
    sCh = as.numeric(sCh),
    # Same uid construction as 01_preprocess_outbreaks.R:119.
    uid = str_c(location, "-", outbreak_number)
  ) %>%
  filterAfrica(africa_iso3, "01b/outbreak-intervals")

cat("     intervals read (post-Africa filter):", nrow(raw), "\n")

# Row-level sanity, reported rather than silently dropped.
n_bad_dates <- sum(is.na(raw$TL) | is.na(raw$TR))
n_rev_dates <- sum(!is.na(raw$TL) & !is.na(raw$TR) & raw$TR < raw$TL)
n_na_cases  <- sum(is.na(raw$sCh))
n_neg_cases <- sum(!is.na(raw$sCh) & raw$sCh < 0)
cat(sprintf("     dropped: %d NA-date, %d reversed-date, %d NA-cases, %d negative-cases row(s)\n",
            n_bad_dates, n_rev_dates, n_na_cases, n_neg_cases))
if (n_neg_cases > 0)
  warning("01b: negative sCh values present and dropped -- inspect the extraction.",
          call. = FALSE)

intervals <- raw %>%
  filter(!is.na(TL), !is.na(TR), TR >= TL, !is.na(sCh), sCh >= 0) %>%
  mutate(interval_id = sprintf("IV%08d", row_number()))

# Non-positive outbreaks: the same drop as 01_preprocess_outbreaks.R
# (--drop_cases_le0), on the (uid, location_period_id) totals gate B compares,
# so the two episode universes stay equal.
if (isTRUE(opt$drop_cases_le0)) {
  le0 <- intervals %>%
    group_by(uid, location_period_id) %>%
    summarise(total = sum(sCh), .groups = "drop") %>%
    filter(total <= 0)
  if (nrow(le0) > 0) {
    warning(sprintf("01b: dropping %d outbreak(s) with total cases <= 0: %s",
                    nrow(le0), paste(le0$uid, collapse = ", ")), call. = FALSE)
    intervals <- anti_join(intervals, le0, by = c("uid", "location_period_id"))
  }
}

cat("     intervals kept:", nrow(intervals), "\n")


# 2. Apportion each interval's cases across the months it overlaps --------
#
# .expandActiveMonths() (utils/dur_data.R:236-265) is reused verbatim: it is
# fully vectorised, needs only (episode_id, TL, TR), and already implements the
# day-overlap arithmetic with the project's inclusive "+1 day" convention. Handing
# it interval_id AS episode_id keeps that arithmetic single-sourced with the
# duration track instead of forking a second copy here. Its `mso` / `n_months`
# outputs are meaningless at interval grain and are dropped.

month_parts <- .expandActiveMonths(
  tibble(episode_id = intervals$interval_id,
         TL         = intervals$TL,
         TR         = intervals$TR)
) %>%
  select(interval_id = episode_id, month, overlap_days = exposure)

# Weights normalised WITHIN interval => sum(w) == 1 exactly, so
# sum(sCh * w) == sCh by construction, whatever the day convention.
month_parts <- month_parts %>%
  group_by(interval_id) %>%
  mutate(w = overlap_days / sum(overlap_days)) %>%
  ungroup()

stopifnot(all(is.finite(month_parts$w)), all(month_parts$w > 0))

apportioned <- month_parts %>%
  inner_join(
    intervals %>% select(interval_id, uid, location, location_period_id,
                         country, sCh, TL, TR),
    by = "interval_id"
  ) %>%
  mutate(
    cases_part  = sCh * w,
    # An interval that touches more than one month had its cases split by THIS
    # script; one wholly inside a month did not.
    straddling  = w < 1
  )


# 3. Collapse to (uid, location_period_id, month) -------------------------
#
# `obs_days` is the number of days in the month actually covered by a reporting
# interval. It is NOT the duration panel's `exposure`, which counts calendar days
# between TL = min and TR = max and therefore silently includes months in which
# nothing was reported at all (the raw intervals can have holes). Carrying both
# lets the fit choose an honest denominator and lets us measure the gap.

monthly <- apportioned %>%
  group_by(uid, location, location_period_id, country, month) %>%
  summarise(
    cases_month           = sum(cases_part),
    obs_days              = sum(overlap_days),
    n_intervals           = n_distinct(interval_id),
    cases_straddling_frac = if (sum(cases_part) > 0)
                              sum(cases_part[straddling]) / sum(cases_part)
                            else 0,
    .groups = "drop"
  )

# Admin level, derived exactly as 01_preprocess_outbreaks.R:127-132 (including
# the TZA::Mainland shift-down). Computed on distinct locations and joined back:
# taxdat::get_admin_level() is a per-string call and there are ~2.6k locations
# against ~50k intervals.
admin_lut <- monthly %>%
  distinct(location) %>%
  mutate(admin_level = map_chr(location, ~ str_c("ADM", taxdat::get_admin_level(.)))) %>%
  mutate(admin_level = if_else(
    str_detect(location, "TZA::Mainland"),
    str_c("ADM", as.integer(str_extract(admin_level, "[0-9]+")) - 1L),
    admin_level
  ))

monthly <- monthly %>%
  left_join(admin_lut, by = "location") %>%
  select(uid, location, location_period_id, country, admin_level, month,
         cases_month, obs_days, n_intervals, cases_straddling_frac) %>%
  arrange(uid, location_period_id, month)


# 4. Gates ----------------------------------------------------------------

cat("\n================ monthly-cases gates ================\n")

# C. Structural integrity of the output. Checked first: the conservation gates
# below are only meaningful on a well-formed key.
n_dup <- sum(duplicated(monthly[, c("uid", "location_period_id", "month")]))
gC <- n_dup == 0 &&
      !any(is.na(monthly$cases_month)) &&
      all(monthly$cases_month >= 0) &&
      !any(is.na(monthly$month)) &&
      !any(is.na(monthly$admin_level))
cat(sprintf("C key + values    : %d duplicate (uid,lp,month), %d NA cases, %d negative ... %s\n",
            n_dup, sum(is.na(monthly$cases_month)),
            sum(monthly$cases_month < 0, na.rm = TRUE), if (gC) "OK" else "FAIL"))
if (!gC) stop("01b: gate C failed -- output key is not unique or values are invalid.")

# A. The apportionment conserves cases, per (uid, lp), against the SAME frame it
# was computed from. This isolates the month-split from any filter difference.
recon_a <- monthly %>%
  group_by(uid, location_period_id) %>%
  summarise(monthly_total = sum(cases_month), .groups = "drop") %>%
  full_join(
    intervals %>%
      group_by(uid, location_period_id) %>%
      summarise(interval_total = sum(sCh), .groups = "drop"),
    by = c("uid", "location_period_id")
  ) %>%
  mutate(
    monthly_total  = coalesce(monthly_total, 0),
    interval_total = coalesce(interval_total, 0),
    rel_err = abs(monthly_total - interval_total) / pmax(1, abs(interval_total))
  )

max_rel_a <- max(recon_a$rel_err)
gA <- max_rel_a <= opt$tolerance
cat(sprintf("A apportionment   : %d (uid,lp) group(s), max rel err %.3e (tol %.1e) ... %s\n",
            nrow(recon_a), max_rel_a, opt$tolerance, if (gA) "OK" else "FAIL"))
if (!gA) {
  print(recon_a %>% arrange(desc(rel_err)) %>% head(10) %>% as.data.frame())
  stop("01b: gate A failed -- the month apportionment does not conserve cases.")
}

# B. This script's episode universe matches 01_preprocess_outbreaks.R's. Tests
# the FILTERS (Africa, date/case validity), not the arithmetic. A failure here
# most often means `--indir` disagrees with the vintage that produced the
# metadata file, which is exactly the drift worth catching.
if (file.exists(opt$metadata)) {
  md <- readRDS(opt$metadata) %>%
    ungroup() %>%
    mutate(location_period_id = as.character(location_period_id)) %>%
    select(uid, location_period_id, episode_total = suspected_cases)

  recon_b <- recon_a %>%
    select(uid, location_period_id, monthly_total) %>%
    full_join(md, by = c("uid", "location_period_id")) %>%
    mutate(
      in_monthly = !is.na(monthly_total),
      in_meta    = !is.na(episode_total),
      rel_err    = abs(coalesce(monthly_total, 0) - coalesce(episode_total, 0)) /
                     pmax(1, abs(coalesce(episode_total, 0)))
    )

  n_only_monthly <- sum(recon_b$in_monthly & !recon_b$in_meta)
  n_only_meta    <- sum(!recon_b$in_monthly & recon_b$in_meta)
  both           <- recon_b %>% filter(in_monthly, in_meta)
  max_rel_b      <- if (nrow(both)) max(both$rel_err) else NA_real_

  gB <- n_only_monthly == 0 && n_only_meta == 0 &&
        !is.na(max_rel_b) && max_rel_b <= opt$tolerance
  cat(sprintf("B vs 01 metadata  : %d shared, %d monthly-only, %d meta-only, max rel err %.3e ... %s\n",
              nrow(both), n_only_monthly, n_only_meta, max_rel_b,
              if (gB) "OK" else "FAIL"))
  if (!gB) {
    cat("   worst mismatches:\n")
    print(both %>% arrange(desc(rel_err)) %>% head(10) %>% as.data.frame())
    if (n_only_meta > 0) {
      cat("   present in metadata but not here (first 10):\n")
      print(recon_b %>% filter(in_meta, !in_monthly) %>% head(10) %>% as.data.frame())
    }
    stop("01b: gate B failed -- this script's episode universe or totals diverge ",
         "from ", opt$metadata, ". Check that --indir matches the vintage that ",
         "produced it (01_preprocess_outbreaks.R defaults to ",
         defaultOutbreakDataDir(), ").")
  }
} else {
  max_rel_b <- NA_real_
  cat(sprintf("B vs 01 metadata  : SKIPPED -- %s not found\n", opt$metadata))
  warning("01b: gate B skipped, ", opt$metadata, " is missing. Run ",
          "01_preprocess_outbreaks.R to enable the cross-check.", call. = FALSE)
}

# D. Report (not enforce) how much of the monthly signal is an artefact of the
# boundary split, and how far obs_days falls short of the calendar month. Both
# are read at fit time, not here.
overlap_check <- monthly %>%
  mutate(days_in_month = as.numeric(days_in_month(month)),
         over = obs_days > days_in_month + 1e-9)
cat(sprintf("D provenance      : %.1f%% of cases from straddling intervals (case-weighted); median obs_days %.0f\n",
            100 * sum(monthly$cases_month * monthly$cases_straddling_frac) /
              max(1, sum(monthly$cases_month)),
            median(monthly$obs_days)))
cat(sprintf("   months where obs_days exceeds the calendar month: %d %s\n",
            sum(overlap_check$over),
            if (sum(overlap_check$over) > 0) "(overlapping reporting intervals!)" else ""))
if (any(overlap_check$over))
  warning("01b: some (uid, lp, month) cells have obs_days > days in month, which ",
          "means the raw reporting intervals OVERLAP. obs_days is then a ",
          "double count and must not be used as an offset. Investigate before ",
          "relying on it.", call. = FALSE)

cat("\n-- rows per month-count bucket\n")
print(monthly %>%
        count(n_intervals) %>%
        arrange(n_intervals) %>%
        head(10) %>%
        as.data.frame())

cat("=====================================================\n\n")


# 5. Save -----------------------------------------------------------------

OUT_PATH <- makeMonthlyCasesFile()
saveRDS(monthly, file = OUT_PATH)
message(sprintf("[01b] wrote %s (%d rows x %d cols, %d (uid,lp) episode(s))",
                OUT_PATH, nrow(monthly), ncol(monthly),
                nrow(distinct(monthly, uid, location_period_id))))

# Policy sidecar. The artefact name is a constant (makeMonthlyCasesFile() takes
# no opt), so --indir and the tolerance would otherwise leave no trace at all --
# and --indir is exactly the parameter that silently changes the contents. Same
# pattern and reasoning as 03_prepare_dur.R:183.
persistPrepPolicyStamp(OUT_PATH, prepPolicyStamp(
  pop_source       = NA_character_,
  composite_policy = NA_character_,
  cover_threshold  = NA_real_,
  composite_tiling = NA_character_,
  pop_digest       = NA_character_,
  counts = list(
    track                 = "monthly_cases",
    indir                 = opt$indir,
    source_file           = ts_path,
    metadata_file         = opt$metadata,
    tolerance             = opt$tolerance,
    n_intervals_in        = nrow(raw),
    n_intervals_kept      = nrow(intervals),
    n_rows_out            = nrow(monthly),
    n_uid_lp              = nrow(distinct(monthly, uid, location_period_id)),
    cases_total           = sum(monthly$cases_month),
    max_recon_rel_err_a   = max_rel_a,
    max_recon_rel_err_b   = max_rel_b,
    frac_cases_straddling = sum(monthly$cases_month * monthly$cases_straddling_frac) /
                              max(1, sum(monthly$cases_month)))))

message(sprintf("[01b] policy stamp written: %s", prepPolicyStampPath(OUT_PATH)))
