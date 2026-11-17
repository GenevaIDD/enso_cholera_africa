# 03_prepare_dur.R — build the outbreak-duration model frame: one row per
#   active episode-month, flagging the month the outbreak ceased.
#
# Usage:
#   Rscript analysis/03_prepare_dur.R -b 1990_2019
#   Rscript analysis/03_prepare_dur.R -b 1990_2019 --episode_key gid_season
#   Rscript analysis/03_prepare_dur.R -b 1990_2019 --episode_key gid_season \
#     --restamp
#
# Details:
#   Feeds the discrete-time cessation hazard (Poisson with log offset). Runs
#   after 03_prepare_siz.R: it expands that script's episode table, which
#   already carries gid / ADM1 / min_seas from the seasonality-match cascade.
#   The episode key is inherited from the size artefact's policy stamp unless
#   --episode_key names one. --restamp rewrites only the .policy.rds stamp.
#
# See: analysis/utils/dur_data.R (likelihood and design rationale)

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(readr)
library(tibble)
library(DBI)
library(RSQLite)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R") # make_climate_lags, apply_offset_months, guard_missing_climate_data
source("analysis/utils/siz_data.R") # .make_size_clim_stats
source("analysis/utils/enso.R")
source("analysis/utils/dur_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
option_list <- c(
  prep_cli_options("dur"),
  list(
    make_option(
      "--size_data",
      default = "",
      type = "character",
      help = "episode table to expand [default: make_size_dat_file() for the episode key]"
    ),
    make_option(
      "--censor_tol_days",
      default = 73L,
      type = "integer",
      help = paste(
        "TR within this many days of the series end => censored",
        "(31 d + the 6-week closing run, which is not part of the outbreak)",
        "[default: %default]"
      )
    ),
    make_option(
      "--max_lag",
      default = 11L,
      type = "integer",
      help = "climate lag columns to build, for DLNM use at fit time [default: %default]"
    ),
    make_option(
      "--enso_source",
      default = "oni",
      type = "character",
      help = "ENSO index source: 'oni' (default) or 'nino34'"
    ),
    make_option(
      "--refresh_enso",
      default = FALSE,
      action = "store_true",
      help = "re-download the ENSO index even if a cache exists"
    ),
    # Also write enso_lead1..k (ONI k months after the active month);
    # the monthly frame inherits them from this panel. Columns only.
    make_option(
      "--enso_max_lead",
      default = 0L,
      type = "integer",
      help = "ENSO lead columns enso_lead1..k to write (0 = none) [default: %default]"
    ),
    make_option(
      "--long_episode_days",
      default = 730L,
      type = "integer",
      help = "episodes longer than this are WARNED about, never dropped [default: %default]"
    ),
    # Overwrite guard: a panel of another prep stamp is not replaced
    # without --force (restamped_at counts as volatile).
    prep_force_option("duration panel"),
    make_option(
      "--restamp",
      default = FALSE,
      action = "store_true",
      help = paste0(
        "rewrite only the existing panel's policy stamp (span + guard ",
        "fields); no panel rebuild [default: %default]"
      )
    )
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
opt <- validate_prep_cli(opt, "dur")
if (!nzchar(opt$size_data)) {
  opt$size_data <- as.character(make_size_dat_file(list(
    episodes = list(
      key = if (identical(opt$episode_key, "inherit")) {
        "lp_season"
      } else {
        opt$episode_key
      }
    )
  )))
}
# The episode key is INHERITED from the size artefact's policy stamp
# (--episode_key inherit, the default, reads the default size frame); passing
# --episode_key gid_season selects the "_ep-gs" size frame and asserts it.
opt <- resolve_inherited_episode_key(opt, opt$size_data, "duration")
gid_season <- identical(opt$episodes$key, "gid_season")
print(opt)

# 3. Data ----------------------------------------------------------------------

if (!file.exists(opt$size_data)) {
  stop(
    "03_prepare_dur.R: episode table not found: ",
    opt$size_data,
    " -- run analysis/03_prepare_siz.R first."
  )
}

episodes <- readRDS(opt$size_data)
message(sprintf(
  "[dur] %d episode(s) read from %s",
  nrow(episodes),
  opt$size_data
))

# What occ_data.R::guard_outbreak_episodes() WOULD drop, counted, never applied
# (dur_data.R header: the 730-day drop would cut the cessation tail). Stamped
# as n_guard_* so the counts reach table S1.
guard <- outbreak_episode_guard_flags(episodes)
guard_counts <- list(
  n_guard_cases_le0 = sum(guard$cases_le0),
  n_guard_cases_lt50 = sum(guard$cases_lt50),
  n_guard_gt730d = sum(guard$gt730d),
  n_guard_any = sum(guard$any)
)
message(sprintf(
  paste0(
    "[dur] guard_outbreak_episodes() would drop %d / %d episode(s) ",
    "(cases<=0 %d, cases<50 %d, >730 d %d); all kept"
  ),
  guard_counts$n_guard_any,
  nrow(episodes),
  guard_counts$n_guard_cases_le0,
  guard_counts$n_guard_cases_lt50,
  guard_counts$n_guard_gt730d
))

# One row per episode of the panel: the gates below, the span summary and the
# stamp's counts all read it. TL / TR feed duration_span_summary().
episode_level <- function(dur_dat) {
  dur_dat %>%
    group_by(episode_id) %>%
    summarise(
      country = first(country),
      n_rows = n(),
      n_months = first(duration_months),
      censored = first(censored),
      n_ended = sum(ended),
      TL = first(TL),
      TR = first(TR),
      .groups = "drop"
    )
}

# 4. --restamp: rewrite only the existing panel's policy stamp -----------------
# read_prep_policy_stamp() / persist_prep_policy_stamp() on the panel the normal run
# wrote; every existing field is kept, the span and guard fields are added or
# refreshed, and `restamped_at` records when. Refused when the panel, its stamp
# or the episode table it was built from does not match. No panel rebuild, no
# ENSO fetch, no climate read.
if (isTRUE(opt$restamp)) {
  DUR_DAT_PATH <- make_duration_dat_file(opt = opt)
  if (!file.exists(DUR_DAT_PATH)) {
    stop(
      "03_prepare_dur.R --restamp: no panel at ",
      DUR_DAT_PATH,
      " -- run without --restamp first."
    )
  }
  old <- read_prep_policy_stamp(DUR_DAT_PATH)
  if (is.null(old)) {
    stop(
      "03_prepare_dur.R --restamp: no readable policy stamp at ",
      prep_policy_stamp_path(DUR_DAT_PATH),
      " -- run without --restamp."
    )
  }
  dur_dat <- readRDS(DUR_DAT_PATH)
  ep_lvl <- episode_level(dur_dat)
  up <- prep_stamp_inherit(opt$size_data)
  chk <- c(
    track = identical(old[["track"]], "duration"),
    source_episodes = identical(old[["source_episodes"]], opt$size_data),
    n_episodes = identical(as.integer(old[["n_episodes"]]), nrow(ep_lvl)),
    n_rows_out = identical(as.integer(old[["n_rows_out"]]), nrow(dur_dat)),
    episode_digest = is.na(old[["episode_digest"]] %||% NA) ||
      is.na(up$episode_digest) ||
      identical(old[["episode_digest"]], up$episode_digest)
  )
  if (!all(chk)) {
    stop(
      "03_prepare_dur.R --restamp: the stamp at ",
      prep_policy_stamp_path(DUR_DAT_PATH),
      " does not describe this panel / episode table (mismatch: ",
      paste(names(chk)[!chk], collapse = ", "),
      "). Rebuild without --restamp."
    )
  }
  new <- utils::modifyList(
    old,
    c(
      duration_span_summary(ep_lvl),
      guard_counts,
      list(restamped_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
    )
  )
  persist_prep_policy_stamp(DUR_DAT_PATH, new)
  changed <- names(new)[
    !vapply(names(new), function(k) identical(old[[k]], new[[k]]), logical(1))
  ]
  message(sprintf(
    "[dur] --restamp: policy stamp rewritten (%s); fields added / changed: %s",
    prep_policy_stamp_path(DUR_DAT_PATH),
    paste(changed, collapse = ", ")
  ))
  quit(save = "no", status = 0)
}

# 5. Build ---------------------------------------------------------------------
enso_dat <- fetch_enso_index(
  source = opt$enso_source,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_enso)
)

dur_dat <- make_duration_panel(
  episodes = episodes,
  baseline = opt$baseline,
  enso_dat = enso_dat,
  enso_max_lead = opt$enso_max_lead,
  max_lag = opt$max_lag,
  censor_tol_days = opt$censor_tol_days,
  long_episode_days = opt$long_episode_days,
  # gid_season: censoring is decided by the member whose report ENDS the
  # canonical episode, not by the per-gid representative LP.
  censor_col = if (gid_season) "censor_lp_id" else "location_period_id"
)

# 6. Verification gates --------------------------------------------------------
# These are the panel's acceptance checks. They run on every invocation rather
# than living in a test file because they are properties of THIS corpus, and a
# corpus refresh is exactly when they would start failing.

ep_lvl <- episode_level(dur_dat)

cat("\n================ duration panel gates ================\n")

# G1. Row count is exactly the sum of active months.
g1 <- nrow(dur_dat) == sum(ep_lvl$n_months)
cat(sprintf(
  "G1 row count      : %d rows vs sum(duration_months) = %d ... %s\n",
  nrow(dur_dat),
  sum(ep_lvl$n_months),
  if (g1) "OK" else "FAIL"
))
if (!g1) {
  stop(
    "03_prepare_dur.R: G1 failed -- the month expansion lost or duplicated rows."
  )
}

# G2. `ended` fires exactly once per uncensored episode and never for a
# censored one. This is the discrete-time-hazard likelihood contract.
n_unc <- sum(!ep_lvl$censored)
g2 <- sum(dur_dat$ended) == n_unc &&
  all(ep_lvl$n_ended == as.integer(!ep_lvl$censored))
cat(sprintf(
  "G2 event count    : sum(ended) = %d vs %d uncensored episode(s) ... %s\n",
  sum(dur_dat$ended),
  n_unc,
  if (g2) "OK" else "FAIL"
))
if (!g2) {
  stop(
    "03_prepare_dur.R: G2 failed -- `ended` does not mark exactly one ",
    "month per uncensored episode."
  )
}

# G3. TZA must come back flagged. Its database is documented as ending in 2019
# (see CLAUDE.md, Data sources), so if the censoring logic cannot find TZA it is
# not finding anything, and every downstream duration estimate is biased short.
# Hard stop, not a warning. Skipped only for the synthetic world (every country
# a synthetic code): real data without TZA fail here.
tza <- ep_lvl %>% filter(country == "TZA")
synth_only <- nrow(ep_lvl) > 0 &&
  all(ep_lvl$country %in% synthetic_countries()$iso3)
g3 <- (synth_only && !nrow(tza)) || sum(tza$censored) > 0
cat(sprintf(
  "G3 TZA censoring  : %d / %d TZA episode(s) flagged ... %s\n",
  sum(tza$censored),
  nrow(tza),
  if (synth_only && !nrow(tza)) {
    "SKIP (synthetic world)"
  } else if (g3) {
    "OK"
  } else {
    "FAIL"
  }
))
if (!g3) {
  stop(
    "03_prepare_dur.R: G3 failed -- TZA has no right-censored episodes, ",
    "but its time series is documented to end in 2019. The series-end ",
    "derivation or the tolerance is wrong. Fix before fitting."
  )
}

# G4. The empirical cessation hazard must vary with months-since-onset. A flat
# hazard would mean a geometric duration distribution, the mso baseline would be
# pointless, and a constant-hazard model would be the honest choice.
haz <- dur_dat %>%
  group_by(mso_band) %>%
  summarise(
    at_risk = n(),
    events = sum(ended),
    exposure_d = sum(exposure),
    haz_month = events / at_risk,
    haz_day = events / exposure_d,
    .groups = "drop"
  )
cat("G4 hazard by mso band:\n")
print(as.data.frame(haz))
rng <- range(haz$haz_month[haz$at_risk >= 30])
cat(sprintf(
  "   monthly hazard range over bands with >=30 at-risk rows: %.3f - %.3f (ratio %.2fx)\n",
  rng[1],
  rng[2],
  rng[2] / rng[1]
))
if (rng[2] / rng[1] < 1.2) {
  warning(
    "03_prepare_dur.R: G4 -- the cessation hazard is close to flat across ",
    "mso bands (ratio < 1.2x). Reconsider duration.basehaz_duration = 'none'."
  )
}

# Observed descriptives the simulator validation compares against.
cat("\n-- observed duration (days) by AFRO_region\n")
print(
  dur_dat %>%
    distinct(episode_id, AFRO_region, duration_days, censored) %>%
    group_by(AFRO_region) %>%
    summarise(
      n = n(),
      n_cens = sum(censored),
      median = median(duration_days),
      p90 = quantile(duration_days, .9),
      .groups = "drop"
    ) %>%
    as.data.frame()
)

cat("=====================================================\n\n")

# 7. Save ----------------------------------------------------------------------

DUR_DAT_PATH <- make_duration_dat_file(opt = opt)

# Policy sidecar. Mirrors 03_prepare_occ.R:804 -- the artefact name is a
# function of the baseline alone, so the parameters that shaped it (censoring
# tolerance above all) would otherwise leave no trace. The duration track has no
# population denominator, hence pop_digest = NA. The population / composite /
# episode fields are COPIED from the size stamp (prep_stamp_inherit()), so all
# four tracks carry the same record and assert_episode_coherence() can name the
# policy behind each side. Built before the panel is written, so the overwrite
# guard can compare it.
up <- prep_stamp_inherit(opt$size_data)
dur_stamp <- prep_policy_stamp(
  pop_source = up$pop_source,
  composite_policy = up$composite_policy,
  cover_threshold = up$cover_threshold,
  composite_tiling = up$composite_tiling,
  pop_digest = NA_character_,
  episode_key = up$episode_key,
  episode_digest = up$episode_digest,
  counts = c(
    list(
      track = "duration",
      baseline = opt$baseline,
      source_episodes = opt$size_data,
      censor_tol_days = opt$censor_tol_days,
      max_lag = opt$max_lag,
      enso_source = opt$enso_source,
      enso_max_lead = opt$enso_max_lead,
      n_episodes = nrow(ep_lvl),
      n_censored = sum(ep_lvl$censored),
      n_rows_out = nrow(dur_dat),
      n_events_out = sum(dur_dat$ended),
      # NB: despite the name, the median number of calendar months TOUCHED
      # (n_months = span + 1), not days. Kept as is (S1 shows it); the span the
      # deploy lag floor is checked against is median_span_months below.
      median_duration_d = median(ep_lvl$n_months)
    ),
    duration_span_summary(ep_lvl), # median_span_months, span_months_q25 / _q75,
    # median_elapsed_months
    guard_counts
  )
) # n_guard_*

# Overwrite guard (the rule of 03_prepare_occ.R): an existing
# panel whose stamp differs in anything but stamped_at / restamped_at (another
# data vintage, another censoring tolerance) or that has none is not replaced
# without --force. A rebuild under an identical stamp passes.
assert_prep_stamp_overwrite(
  DUR_DAT_PATH,
  dur_stamp,
  force = isTRUE(opt$force),
  rerun_hint = "Rscript analysis/03_prepare_dur.R ... --force",
  ignore = c("stamped_at", "restamped_at")
)

saveRDS(dur_dat, file = DUR_DAT_PATH)
message(sprintf(
  "[dur] wrote %s (%d rows x %d cols)",
  DUR_DAT_PATH,
  nrow(dur_dat),
  ncol(dur_dat)
))
persist_prep_policy_stamp(DUR_DAT_PATH, dur_stamp)

message(sprintf(
  "[dur] policy stamp written: %s",
  prep_policy_stamp_path(DUR_DAT_PATH)
))
