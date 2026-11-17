# 03_prepare_mon.R — build the monthly case-count model frame: a per-month
#   case count on the duration track's episode x active-month skeleton.
#
# Usage:
#   Rscript analysis/03_prepare_mon.R -b 1990_2019
#
# Details:
#   Runs after 03_prepare_dur.R. Adds one column of substance (`cases_month`)
#   and its denominators; every covariate (mso, exposure, pop, ENSO and
#   climate lags) is inherited from the duration panel unchanged, so the
#   match cascade and lag joins are not copied. Cases from
#   01b_preprocess_outbreaks_monthly.R are joined through the footprint
#   membership map on (gid, year_seas, admin_level), not on
#   (uid, location_period_id, month); gates M1-M4 check the join.
#
# See: analysis/utils/footprint.R

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(tibble)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/provenance.R") # prep_policy_stamp, persist_prep_policy_stamp
source("analysis/utils/paths.R") # make_duration_dat_file, make_monthly_cases_file, ...
source("analysis/utils/footprint.R") # aggregate_monthly_to_footprint, FOOTPRINT_KEY

# 2. Options -------------------------------------------------------------------
option_list <- c(
  prep_cli_options("mon"),
  list(
    make_option(
      "--dur_data",
      default = "",
      type = "character",
      help = "duration panel to attach cases to [default: make_duration_dat_file()]"
    ),
    make_option(
      "--monthly_cases",
      default = "",
      type = "character",
      help = "per-(uid, lp, month) case counts [default: make_monthly_cases_file()]"
    ),
    make_option(
      "--members",
      default = "",
      type = "character",
      help = "footprint membership map [default: make_footprint_members_file()]"
    ),
    make_option(
      "--tolerance",
      default = 1e-6,
      type = "double",
      help = "relative tolerance for the conservation gates [default: %default]"
    ),
    # Out-of-window case mass (gate M4 below) is a REAL property of the footprint
    # collapse, not a bug to be tuned away: an episode's [TL, TR] is its
    # representative location period's window, while the cases summed into it come
    # from siblings whose windows differ. The threshold exists so the fraction has
    # to be looked at once and consciously accepted, rather than discovered later.
    make_option(
      "--max_out_of_window_frac",
      default = 0.05,
      type = "double",
      help = "fail if more than this share of footprint case mass falls outside the episode window [default: %default]"
    ),
    # Overwrite guard: a panel of another prep stamp is not replaced
    # without --force.
    prep_force_option("monthly panel")
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
opt <- validate_prep_cli(opt, "mon")

# The episode key is inherited from the duration panel's stamp; an explicit
# --episode_key selects the matching "_ep-gs" inputs and asserts the stamp.
opt_in <- c(
  opt,
  list(
    episodes = list(
      key = if (identical(opt$episode_key, "inherit")) {
        "lp_season"
      } else {
        opt$episode_key
      }
    )
  )
)
DUR_PATH <- if (nzchar(opt$dur_data)) {
  opt$dur_data
} else {
  make_duration_dat_file(opt_in)
}
opt <- resolve_inherited_episode_key(opt, DUR_PATH, "monthly")
gid_season <- identical(opt$episodes$key, "gid_season")
MON_PATH <- if (nzchar(opt$monthly_cases)) {
  opt$monthly_cases
} else {
  make_monthly_cases_file()
}
MEM_PATH <- if (nzchar(opt$members)) {
  opt$members
} else {
  make_footprint_members_file(opt)
}
# The episode key: (gid, year_seas, admin_level), or (gid, year_seas) under
# gid_season, where an episode does not split across admin levels.
KEY <- footprint_key(opt$episodes$key)

print(opt)

# 3. Data ----------------------------------------------------------------------

for (p in c(DUR_PATH, MON_PATH, MEM_PATH)) {
  if (!file.exists(p)) {
    stop(
      "03_prepare_mon.R: input not found: ",
      p,
      "\n  needs, in order: analysis/01b_preprocess_outbreaks_monthly.R,",
      "\n                   analysis/03_prepare_siz.R (writes the membership map),",
      "\n                   analysis/03_prepare_dur.R"
    )
  }
}

dur_dat <- readRDS(DUR_PATH)
monthly <- readRDS(MON_PATH)
members <- readRDS(MEM_PATH)

message(sprintf(
  "[mon] duration panel : %s (%d rows, %d episodes)",
  DUR_PATH,
  nrow(dur_dat),
  n_distinct(dur_dat$episode_id)
))
message(sprintf(
  "[mon] monthly cases  : %s (%d rows, %d (uid, lp))",
  MON_PATH,
  nrow(monthly),
  n_distinct(paste(monthly$uid, monthly$location_period_id))
))
message(sprintf(
  "[mon] members map    : %s (%d member rows, %d footprints)",
  MEM_PATH,
  nrow(members),
  n_distinct(paste(members$gid, members$year_seas, members$admin_level))
))

# 3a. Freshness ----------------------------------------------------------------
# The membership map is written as a pair with the size frame, and the duration
# panel is built FROM that size frame. If the panel is older than the map, it
# describes a different episode universe and every gate below is measuring the
# wrong thing (e.g. differing `pop` on every row, or episodes the refreshed
# corpus no longer has). Reported, not fatal -- the numeric gates are the real defence,
# and this line is what makes
# their failure interpretable.
{
  mt <- file.mtime(c(DUR_PATH, MON_PATH, MEM_PATH))
  names(mt) <- c("dur_panel", "monthly_cases", "members")
  cat("\n-- input mtimes\n")
  print(mt)
  if (mt[["dur_panel"]] < mt[["members"]]) {
    warning(
      "03_prepare_mon.R: the duration panel (",
      format(mt[["dur_panel"]]),
      ") is OLDER than the footprint membership map (",
      format(mt[["members"]]),
      "). Re-run 03_prepare_dur.R before ",
      "trusting this frame.",
      call. = FALSE
    )
  }
}
# Under gid_season the pairing is checked by CONTENT, not mtime: the members
# map and the duration panel must carry the same episode digest.
if (gid_season) {
  assert_episode_coherence(
    list(
      dur = episode_stamp_of(DUR_PATH),
      members = episode_stamp_of(MEM_PATH)
    ),
    what = "03_prepare_mon.R inputs",
    require_gid_season = TRUE
  )
}

# Only the members whose cases the canonical episode total counts (the finest
# reporting level, `in_total`) feed the monthly series -- the same members the
# size/duration totals were built from, so gate M3 stays an equality.
if (gid_season) {
  members <- members %>% filter(in_total)
}

# 4. Build ---------------------------------------------------------------------
# Monthly cases within an episode, for a model of the form
#   cases_{i,m} ~ NegBin(mu_{i,m}, theta),
#   mu_{i,m} = exp(eta_{i,m}) * pop_i * exposure_days_{i,m}.
# Covariates come from the duration panel (one copy of the lag joins).
# Monthly cases (01b, keyed on uid, location_period_id, month) are joined
# through 03_prepare_siz.R's footprint membership map on (gid, year_seas,
# admin_level), NOT on the episode's uid / location_period_id: those are a
# representative label, and that join would silently keep only the label LP's
# cases (footprint.R).

# Episode-level footprint totals, taken from the panel itself so the
# conservation gate compares the monthly series against the very numbers the
# model's siblings (size, duration) were fitted on.
episodes <- dur_dat %>%
  distinct(
    across(all_of(KEY)),
    episode_id,
    uid,
    location_period_id,
    TL,
    TR,
    suspected_cases,
    pop,
    country
  ) %>%
  mutate(across(all_of(KEY), as.character))

if (anyDuplicated(episodes[, KEY])) {
  stop(
    "03_prepare_mon.R: (",
    paste(KEY, collapse = ", "),
    ") is not unique on ",
    "the duration panel -- it is not the episode key it is assumed to be ",
    "here. Do not proceed; the join below would fan out."
  )
}

agg <- aggregate_monthly_to_footprint(
  monthly = monthly,
  members = members,
  footprints = episodes %>% select(all_of(KEY), suspected_cases),
  tolerance = opt$tolerance,
  key = KEY
)

fp_monthly <- agg$monthly
recon <- agg$recon
cat("\n-- footprint aggregation report\n")
print(as.data.frame(agg$report))

# 4a. Out-of-window months -----------------------------------------------------
# The footprint's monthly series can extend beyond the episode's [TL, TR],
# because that window is `first()`-ed from the representative location period
# while the cases come from siblings with different windows. Those months have
# no `mso` and no `exposure`, so there is nowhere on the panel to put them.
# Quantified here and dropped; gate M4 refuses to proceed if the mass is large.
in_window <- fp_monthly %>%
  inner_join(
    episodes %>% select(all_of(KEY), episode_id, TL, TR),
    by = KEY,
    relationship = "many-to-one"
  ) %>%
  mutate(
    month_end = ceiling_date(month, "month") - 1L,
    overlaps = month_end >= as.Date(TL) & month <= as.Date(TR)
  )

n_oow <- sum(!in_window$overlaps)
mass_oow <- sum(in_window$cases_month[!in_window$overlaps])
frac_oow <- mass_oow / max(1, sum(in_window$cases_month))

# 4b. Join onto the panel ------------------------------------------------------

n_before <- nrow(dur_dat)
mon_dat <- dur_dat %>%
  mutate(across(all_of(KEY), as.character)) %>%
  left_join(
    fp_monthly %>%
      select(all_of(KEY), month, cases_month, obs_days, n_lp, n_intervals),
    by = c(KEY, "month"),
    relationship = "many-to-one"
  )

# Panel months the footprint reports nothing for. These are REPORTING GAPS
# inside [TL, TR] -- the raw interval series can have holes while the episode's
# TL = min and TR = max span them -- not missing data to be imputed. Counted
# and named before being set to 0, because the failure mode to avoid is an
# implicit filter(cases_month > 0) silently deleting a fifth of the panel.
gap <- mon_dat %>% filter(is.na(cases_month))
cat(sprintf(
  "\n-- reporting-gap months: %d panel row(s) over %d episode(s) have no reported cases\n",
  nrow(gap),
  n_distinct(gap$episode_id)
))
if (nrow(gap)) {
  print(gap %>% count(country, sort = TRUE) %>% head(10) %>% as.data.frame())
}

mon_dat <- mon_dat %>%
  mutate(
    cases_gap = is.na(cases_month),
    cases_month = coalesce(cases_month, 0),
    obs_days = coalesce(obs_days, 0),
    n_lp = coalesce(n_lp, 0L),
    n_intervals = coalesce(n_intervals, 0L)
  )

# 5. Verification gates --------------------------------------------------------
# Same rationale as 03_prepare_dur.R's G1-G4: these are properties of THIS
# corpus, so they run on every invocation rather than living in a test file. A
# corpus refresh is exactly when they start failing.

cat("\n================ monthly panel gates ================\n")

# M1. The join is a lookup, not a fan-out.
m1 <- nrow(mon_dat) == n_before
cat(sprintf(
  "M1 row count      : %d rows in, %d out ... %s\n",
  n_before,
  nrow(mon_dat),
  if (m1) "OK" else "FAIL"
))
if (!m1) {
  stop(
    "03_prepare_mon.R: M1 failed -- the footprint join changed the row count."
  )
}

# M2. Response validity. NegBin cannot take a negative or missing count.
m2 <- !any(is.na(mon_dat$cases_month)) && all(mon_dat$cases_month >= 0)
cat(sprintf(
  "M2 response       : %d NA, %d negative ... %s\n",
  sum(is.na(mon_dat$cases_month)),
  sum(mon_dat$cases_month < 0),
  if (m2) "OK" else "FAIL"
))
if (!m2) {
  stop("03_prepare_mon.R: M2 failed -- cases_month is not a valid count.")
}

# M3. THE LOAD-BEARING GATE. Every footprint's monthly cases must sum to the
# episode total the size/duration tracks were fitted on. Both are sums of the
# same member location periods' cases -- one grouped by month, one not -- so
# anything other than equality means a member's monthly series is missing or
# came from a different corpus vintage. This is what makes the (uid, lp) ->
# footprint remapping verifiable rather than merely plausible.
n_bad <- sum(!recon$conserved)
m3 <- n_bad == 0L
cat(sprintf(
  "M3 conservation   : %d/%d footprint(s) off, max rel err %.3e (tol %.1e) ... %s\n",
  n_bad,
  nrow(recon),
  max(recon$rel_err),
  opt$tolerance,
  if (m3) "OK" else "FAIL"
))
if (!m3) {
  print(
    recon %>%
      arrange(desc(rel_err)) %>%
      select(all_of(KEY), monthly_cases, episode_cases, rel_err) %>%
      head(15) %>%
      as.data.frame()
  )
  stop(
    "03_prepare_mon.R: M3 failed -- the monthly series does not reconcile to ",
    "the episode totals. Most likely the duration panel and ",
    basename(MON_PATH),
    " come from different corpus vintages: rebuild ",
    "01b -> 03_prepare_siz.R -> 03_prepare_dur.R in order."
  )
}

# M4. Out-of-window mass. Not an error condition per se -- see the flag's help
# text -- but it is case mass the model will never see, so it has to be small
# and it has to be looked at.
m4 <- frac_oow <= opt$max_out_of_window_frac
cat(sprintf(
  "M4 out-of-window  : %d month-row(s), %.0f cases (%.2f%% of mass; max %.2f%%) ... %s\n",
  n_oow,
  mass_oow,
  100 * frac_oow,
  100 * opt$max_out_of_window_frac,
  if (m4) "OK" else "FAIL"
))
if (!m4) {
  stop(
    "03_prepare_mon.R: M4 failed -- too much footprint case mass falls ",
    "outside the episode windows. The (gid, year_seas, admin_level) collapse ",
    "is merging location periods whose active periods barely overlap; ",
    "inspect before raising --max_out_of_window_frac."
  )
}

# M5. The response must actually vary within an outbreak, otherwise s(mso) has
# nothing to fit and the whole point of moving to month grain is lost.
within <- mon_dat %>%
  group_by(episode_id) %>%
  filter(n() > 1) %>%
  summarise(
    cv = sd(cases_month) / pmax(1e-9, mean(cases_month)),
    .groups = "drop"
  )
cat(sprintf(
  "M5 within-episode : %d multi-month episode(s), median CV of cases_month %.2f\n",
  nrow(within),
  median(within$cv, na.rm = TRUE)
))
if (nrow(within) && median(within$cv, na.rm = TRUE) < 0.05) {
  warning(
    "03_prepare_mon.R: M5 -- cases are nearly flat within outbreaks ",
    "(median CV < 0.05). s(mso) will have little to fit."
  )
}

# 6. Descriptives --------------------------------------------------------------

cat("\n-- cases_month by months-since-onset band\n")
print(
  mon_dat %>%
    group_by(mso_band) %>%
    summarise(
      n_rows = n(),
      n_gap = sum(cases_gap),
      cases = sum(cases_month),
      median = median(cases_month),
      p90 = quantile(cases_month, .9),
      per_1e5_pd = 1e5 * sum(cases_month) / sum(as.numeric(pop) * exposure),
      .groups = "drop"
    ) %>%
    as.data.frame()
)

cat("\n-- footprint multiplicity actually exercised\n")
print(
  mon_dat %>%
    distinct(episode_id, n_lp) %>%
    count(n_lp) %>%
    as.data.frame()
)

cat("=====================================================\n\n")

# 7. Save ----------------------------------------------------------------------

MON_DAT_PATH <- make_monthly_dat_file(opt = opt)

# Policy sidecar. Mirrors the 03_prepare_dur.R one -- the artefact name is a
# function of the baseline alone, so the inputs and the two judgement calls
# (gap months set to 0, out-of-window mass dropped) would otherwise leave no
# trace. No population denominator is chosen here, hence pop_digest = NA; the
# population / composite / episode fields are copied from the duration stamp
# (which copied them from size), so all four tracks carry one record shape.
# Built before the panel is written (same fields and values), so
# the overwrite guard can compare it.
up <- prep_stamp_inherit(DUR_PATH)
mon_stamp <- prep_policy_stamp(
  pop_source = up$pop_source,
  composite_policy = up$composite_policy,
  cover_threshold = up$cover_threshold,
  composite_tiling = up$composite_tiling,
  pop_digest = NA_character_,
  episode_key = up$episode_key,
  episode_digest = up$episode_digest,
  counts = list(
    track = "monthly",
    baseline = opt$baseline,
    source_dur_panel = DUR_PATH,
    source_monthly_cases = MON_PATH,
    source_members = MEM_PATH,
    n_rows_out = nrow(mon_dat),
    n_episodes = n_distinct(mon_dat$episode_id),
    cases_total = sum(mon_dat$cases_month),
    n_gap_rows = sum(mon_dat$cases_gap),
    n_gap_episodes = n_distinct(mon_dat$episode_id[mon_dat$cases_gap]),
    n_out_of_window_rows = as.integer(n_oow),
    cases_out_of_window = mass_oow,
    frac_out_of_window = frac_oow,
    n_member_unmatched = agg$report$n_member_unmatched,
    max_recon_rel_err = agg$report$max_rel_err,
    n_multi_lp_episodes = n_distinct(mon_dat$episode_id[mon_dat$n_lp > 1L])
  )
)

# Overwrite guard (the rule of 03_prepare_occ.R): an existing
# panel whose stamp differs in anything but stamped_at (another data vintage)
# or that has none is not replaced without --force.
assert_prep_stamp_overwrite(
  MON_DAT_PATH,
  mon_stamp,
  force = isTRUE(opt$force),
  rerun_hint = "Rscript analysis/03_prepare_mon.R ... --force"
)

saveRDS(mon_dat, file = MON_DAT_PATH)
message(sprintf(
  "[mon] wrote %s (%d rows x %d cols)",
  MON_DAT_PATH,
  nrow(mon_dat),
  ncol(mon_dat)
))
persist_prep_policy_stamp(MON_DAT_PATH, mon_stamp)

message(sprintf(
  "[mon] policy stamp written: %s",
  prep_policy_stamp_path(MON_DAT_PATH)
))
