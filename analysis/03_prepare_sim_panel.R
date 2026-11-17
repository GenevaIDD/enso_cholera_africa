# 03_prepare_sim_panel.R — build a continuous per-(gid, month) climate and
#   ENSO panel for the cholera simulator over the observational range.
#
# Usage:
#   Rscript analysis/03_prepare_sim_panel.R -b 1990_2019
#   Rscript analysis/03_prepare_sim_panel.R -b 1990_2019 \
#     --episode_key gid_season --enso_max_lead 11
#
# Details:
#   One row per (gid, calendar month) over [cal_start - warmup_months,
#   cal_end], with every column the simulator consumer contract requires.
#   Writes generated_data/sim_climate_panel_base-<baseline>.rds
#   (make_sim_panel_file()) + its .policy.rds stamp. 03_prepare_occ.R still
#   feeds model fitting and the observation panel.

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(DBI)
library(RSQLite)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/simulator_core.R") # derive_min_seas_from_panel
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------

option_list <- list(
  make_option(
    c("-b", "--baseline"),
    default = "1990_2019",
    type = "character",
    help = "Baseline period; must match the occurrence fit [default: %default]"
  ),
  make_option(
    "--cal_start",
    default = "2010-01-01",
    type = "character",
    help = "First calendar month of the output window [default: %default]"
  ),
  make_option(
    "--cal_end",
    default = "2024-12-01",
    type = "character",
    help = "Last calendar month of the output window [default: %default]"
  ),
  make_option(
    "--warmup_months",
    default = 12L,
    type = "integer",
    help = "Months prepended before cal_start for warm-up state [default: %default]"
  ),
  make_option(
    "--enso_source",
    default = "oni",
    type = "character",
    help = "ENSO index source: 'oni' (default) or 'nino34' [default: %default]"
  ),
  make_option(
    "--enso_max_lead",
    default = 0L,
    type = "integer",
    help = paste(
      "Also write enso_lead1..K (ONI K months after the month; base.R::enso_lag_col_name()),",
      "for deployed ENSO-lead fits. 0 = none [default: %default]"
    )
  ),
  make_option(
    "--refresh_enso",
    default = TRUE,
    action = "store_true",
    help = "Re-download the ENSO index even if a cache exists [default: %default]"
  ),
  make_option(
    "--no_refresh_enso",
    action = "store_false",
    dest = "refresh_enso",
    help = "Use the cached ENSO index (offline runs, the synthetic world)"
  ),
  make_option(
    "--skip_climate_check",
    default = TRUE,
    action = "store_true",
    help = paste(
      "Warn instead of error when climate columns have NAs in the",
      "calendar window. Use when only ENSO covariates are needed",
      "and climate DB coverage is shorter than cal_end [default: %default]"
    )
  )
)

option_list <- c(
  option_list,
  list(
    make_option(
      "--episode_key",
      default = "lp_season",
      type = "character",
      help = paste0(
        "lp_season | gid_season: which survival panel defines ",
        "the gid scope and per-gid season calendar; the ",
        "output carries the matching _ep-gs token [default: %default]"
      )
    ),
    # Overwrite guard: a panel of another prep stamp is not replaced
    # without --force.
    prep_force_option("sim panel")
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$episode_key %in% EPISODE_KEYS) {
  stop("--episode_key must be one of ", paste(EPISODE_KEYS, collapse = " | "))
}
opt$episodes <- list(key = opt$episode_key)
print(opt)

cal_start <- as.Date(opt$cal_start)
cal_end <- as.Date(opt$cal_end)
warmup_start <- cal_start %m-% months(as.integer(opt$warmup_months))

message("-- Building continuous sim climate panel")
message("   baseline:      ", opt$baseline)
message(
  "   window:        ",
  warmup_start,
  " – ",
  cal_end,
  "  (warmup: ",
  warmup_start,
  " – ",
  cal_start - 1L,
  ")"
)

# 3. Scope + min_seas from survival panel --------------------------------------
# Why a separate panel: the survival-fitting dataset (03_prepare_occ.R
# output) only carries rows for each location's at-risk months, so per-gid
# coverage ends at the last at-risk month (~2019-07 median). The simulator's
# build_climate_panel_for_tu() hard-errors when the panel ends > 1 month
# before horizon.end, hence one row per (gid, calendar month) here.

message("-- Loading survival panel for gid scope and min_seas")
surv_panel <- readRDS(make_surv_dat_file(opt))

gids_meta <- surv_panel %>% distinct(gid, country, AFRO_region)
min_seas_tbl <- derive_min_seas_from_panel(surv_panel) # tibble(gid, min_seas)

message("   ", nrow(gids_meta), " gids in scope")

# 4. Continuous month grid -----------------------------------------------------

message("-- Building month grid")
months_seq <- seq.Date(warmup_start, cal_end, by = "1 month")

panel <- expand_grid(gid = gids_meta$gid, t_start = months_seq) %>%
  left_join(gids_meta, by = "gid") %>%
  left_join(min_seas_tbl, by = "gid") %>%
  mutate(
    t_start_seas = yday(t_start %m+% months(-(min_seas - 1L))) - 1L
  )

message(
  "   ",
  nrow(panel),
  " rows  (",
  nrow(gids_meta),
  " gids × ",
  length(months_seq),
  " months)"
)

# 5. Batch-query climate -------------------------------------------------------

message("-- Querying climate_stats and climate_extremes for all gids")
all_gids <- unique(panel$gid)
clim_stats <- get_stats_gid(gids = all_gids)
clim_extrms <- get_extremes_gid(gids = all_gids, baseline = opt$baseline)

clim_dev <- clim_stats %>%
  pivot_longer(cols = c(get_precip_vars(), get_temp_vars())) %>%
  inner_join(
    clim_extrms %>% select(gid, month, name, mean, sd),
    by = c("gid", "month", "name")
  ) %>%
  mutate(dev = value - mean, dev_sd = (value - mean) / sd) %>%
  select(gid, month, name, mean, dev, dev_sd) %>%
  pivot_wider(names_from = "name", values_from = c("mean", "dev", "dev_sd"))

# 6. Left-join climate onto panel ----------------------------------------------
# NOT inner_join: silent drops would re-create the coverage gap that this
# script is designed to fix.

message("-- Joining climate onto panel")
panel <- panel %>%
  left_join(clim_stats, by = c("gid", "t_start" = "month")) %>%
  left_join(clim_dev, by = c("gid", "t_start" = "month"))

# 7. Assert no NAs in the required kernel columns ------------------------------
# The 4 required kernel columns, within the calendar window (warmup rows may
# legitimately have NAs before 1980-01-01).

req_cols <- c(
  "mean_precip",
  "dev_precip",
  "mean_temp_daymean",
  "dev_temp_daymean"
)
in_window <- panel %>% filter(t_start >= cal_start, t_start <= cal_end)
na_counts <- in_window %>%
  summarise(across(all_of(req_cols), ~ sum(is.na(.)), .names = "na_{.col}"))
if (any(na_counts > 0)) {
  print(na_counts)
  if (isTRUE(opt$skip_climate_check)) {
    warning(
      "NA climate values present in calendar window — proceeding because ",
      "--skip_climate_check is set. Climate columns will be NA for uncovered months."
    )
  } else {
    stop(
      "NA climate values in required columns within the calendar window — ",
      "check DB coverage (climate_stats.sqlite / climate_extremes.sqlite)"
    )
  }
}
message("   NA check passed for ", nrow(in_window), " in-window rows")

# 8. ENSO lags (enso_lag0 … enso_lag12) ----------------------------------------

message("-- Fetching ENSO index and joining lags")
enso_dat <- fetch_enso_index(
  source = opt$enso_source,
  cache_dir = "generated_data/",
  refresh = isTRUE(opt$refresh_enso)
)
panel <- join_enso_lags(
  panel,
  enso_dat,
  date_col = "t_start",
  max_lead = opt$enso_max_lead %||% 0L
)

# 9. Warm-up lag columns from survival panel -----------------------------------
# Left-join by (gid, t_start): rows outside the at-risk window get NA,
# which is safe — build_warmup_state() falls back to 0 via pick() and
# ttlo_init_from_panel() maps NA → top TTLO bin.

message("-- Joining warm-up lag columns from survival panel")
lag_cols <- c(
  "ongoing_outbreak_lag1",
  "ongoing_outbreak_country_lag1",
  "ongoing_outbreak_tu_lag1",
  "time_to_last_outbreak"
)
# Deduplicate to one row per (gid, t_start): the binary lag indicators are
# keyed at the ADM1/country/TU level so they are identical across uids for
# the same (gid, t_start). For time_to_last_outbreak we take the minimum
# (most recent outbreak across location periods at the same ADM2).
# ongoing_outbreak_tu_lag1 may be NA for gids with no TU membership — preserve
# that NA so the spillover term is absent from the LP in build_warmup_state.
surv_lags <- surv_panel %>%
  select(gid, t_start, any_of(lag_cols)) %>%
  group_by(gid, t_start) %>%
  summarise(
    across(
      any_of(c("ongoing_outbreak_lag1", "ongoing_outbreak_country_lag1")),
      ~ as.integer(max(coalesce(., 0L))),
      .names = "{.col}"
    ),
    across(
      any_of("ongoing_outbreak_tu_lag1"),
      ~ {
        v <- .[!is.na(.)]
        if (length(v) == 0) NA_integer_ else as.integer(max(v))
      },
      .names = "{.col}"
    ),
    across(
      any_of("time_to_last_outbreak"),
      ~ suppressWarnings(min(., na.rm = TRUE)),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

panel <- panel %>%
  left_join(surv_lags, by = c("gid", "t_start")) %>%
  make_time_to_last_outbreak_cat()

# 10. Climate lag columns for DLNM cross-basis reconstruction at sim time ------
# Add <var>_lag1 … _lag11 for every (mean|dev)_ climate column (excluding
# dev_sd_* and any existing _lag\d+ columns). Consumed by build_cb_mats()
# in simulator_core.R when cb_meta names a climate variable (not "enso").
# Continuous calendar lags per gid — no season-year reset — matching the
# ENSO/join_enso_lags calendar convention. Lag k = NA for the first k rows of
# each gid; build_cb_mats() replaces NA with 0 before crossbasis().

message("-- Computing climate lag columns (lag1..lag11)")
clim_lag_vars <- grep("^(mean|dev)_", names(panel), value = TRUE)
clim_lag_vars <- clim_lag_vars[!grepl("^dev_sd_", clim_lag_vars)]
clim_lag_vars <- clim_lag_vars[!grepl("_lag\\d+$", clim_lag_vars)]

max_lag_clim <- 11L
lag_fns <- setNames(
  lapply(seq_len(max_lag_clim), function(k) {
    force(k)
    function(x) dplyr::lag(x, n = k)
  }),
  paste0("lag", seq_len(max_lag_clim))
)

panel <- panel %>%
  dplyr::arrange(gid, t_start) %>%
  dplyr::group_by(gid) %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(clim_lag_vars),
    lag_fns,
    .names = "{.col}_{.fn}"
  )) %>%
  dplyr::ungroup()

message(
  "   Added ",
  max_lag_clim * length(clim_lag_vars),
  " lag columns (",
  length(clim_lag_vars),
  " vars × lag1..",
  max_lag_clim,
  ")"
)

# 11. Save ---------------------------------------------------------------------

out_path <- make_sim_panel_file(opt)
message("-- Saving to ", out_path)
# Episode record of the survival panel this panel was scoped from, so the
# simulator can refuse to pair it with fits/panels of another episode
# definition (utils/provenance.R::sim_episode_check()). Built before the panel is
# written.
.up <- prep_stamp_inherit(make_surv_dat_file(opt))
.stamp <- prep_policy_stamp(
  pop_source = .up$pop_source,
  composite_policy = .up$composite_policy,
  cover_threshold = .up$cover_threshold,
  composite_tiling = .up$composite_tiling,
  episode_key = .up$episode_key,
  episode_digest = .up$episode_digest,
  counts = list(
    track = "sim_panel",
    baseline = opt$baseline,
    enso_max_lead = as.integer(opt$enso_max_lead %||% 0L),
    source_surv_panel = as.character(make_surv_dat_file(opt))
  )
)
# Overwrite guard (the rule of 03_prepare_occ.R): a panel of
# another stamp (another vintage's survival panel, another lead count) or with
# none is not replaced without --force. The stamp carries no ENSO content: an
# ENSO-only change under the same stamp passes (back up and diff the panels).
assert_prep_stamp_overwrite(
  out_path,
  .stamp,
  force = isTRUE(opt$force),
  rerun_hint = "Rscript analysis/03_prepare_sim_panel.R ... --force"
)
saveRDS(panel, file = out_path)
persist_prep_policy_stamp(out_path, .stamp)
message("   Done. ", nrow(panel), " rows written.")
