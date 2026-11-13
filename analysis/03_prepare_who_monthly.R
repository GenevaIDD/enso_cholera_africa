# 03_prepare_who_monthly.R — WHO national cholera counts → monthly counts
#
# Purpose: monthly suspected cholera cases per country, the observation of the
#   WHO term in the historical smoothing (simulator.filter.who_obs, measurement
#   occurrence_duration_mon_NB_state_detect_act_who). Two sources:
#   (1) the WHO Global Cholera & AWD dashboard weekly series, each week's cases
#       split across calendar months by day overlap
#       (utils/who_dashboard.R::whoWeeklyToMonthly());
#   (2) the year-to-date cumulative snapshots transcribed from the WHO
#       multi-country cholera epidemiological updates (Table 1), differenced
#       between report cuts and split the same way
#       (whoCumulativeToMonthly()). They fill only the (country, month) pairs
#       the weekly feed lacks (whoMergeMonthly()); a weekly row always wins.
#   Unreported weeks / cuts are gaps, never zeros; a month is complete when
#   every day is covered. An unchanged cumulative is read as a reported zero
#   only for the countries in --zero_iso3 (or when the table printed an
#   explicit 0 in its last-28-days column); a year-to-date span longer than
#   --max_interval_days (a country's first mid-year table, a lone December
#   cut) has no monthly information and is dropped.
# Inputs:  data/who_annual_report_africa/who_dashboard_weekly_africa.rds
#          data/who_annual_report_africa/who_sitrep_cumulative_africa.csv
#            (iso3, table_as_of, cum_cases, last28_cases, report_id, source_url;
#             seeded from parsed WHO situation-report tables)
# Outputs: generated_data/who_dashboard_monthly_africa.rds (makeWhoMonthlyFile())
#          columns iso3, t_start, cases, n_days, partial, source
#          generated_data/who_dashboard_monthly_africa_meta.yml (inputs, zero
#          allowlist, row counts by source)
# Usage:   Rscript analysis/03_prepare_who_monthly.R --zero_iso3 TZA
#          Rscript analysis/03_prepare_who_monthly.R --input <weekly.rds> --cumulative <snapshots.csv> --outdir generated_data
#          Rscript analysis/03_prepare_who_monthly.R --cumulative ""   # weekly only

suppressMessages({
  library(dplyr)
  library(optparse)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/who_dashboard.R")
source("analysis/utils/provenance.R")   # assertPrepStampOverwrite, whoMonthlyGuardStamp, prepForceOption

option_list <- list(
  make_option("--input",
              default = "data/who_annual_report_africa/who_dashboard_weekly_africa.rds",
              type = "character", help = "WHO dashboard weekly rds [default: %default]"),
  make_option("--cumulative",
              default = "data/who_annual_report_africa/who_sitrep_cumulative_africa.csv",
              type = "character",
              help = "WHO multi-country update cumulative snapshots csv; \"\" = weekly only [default: %default]"),
  make_option("--zero_iso3", default = "", type = "character",
              help = "Comma-separated ISO3 whose unchanged cumulative is a reported zero [default: none]"),
  make_option("--max_interval_days", default = 45L, type = "integer",
              help = "Longest cut-to-cut interval spread over months; longer year-to-date spans (a country's first mid-year table, a lone December cut) are dropped [default: %default]"),
  make_option("--outdir", default = "generated_data", type = "character",
              help = "Output directory [default: %default]"),
  # Overwrite guard: a table whose meta (less its timestamps) or
  # content differs is not replaced without --force.
  prepForceOption("WHO monthly table")
)
opt <- parse_args(OptionParser(option_list = option_list))
zero_iso3 <- toupper(trimws(strsplit(opt$zero_iso3, ",")[[1]]))
zero_iso3 <- zero_iso3[nzchar(zero_iso3)]

if (!file.exists(opt$input)) stop("[who_monthly] input not found: ", opt$input)
wk  <- readRDS(opt$input)
mon_wk <- whoWeeklyToMonthly(wk)

mon_cum <- NULL
use_cum <- nzchar(opt$cumulative)
if (use_cum) {
  if (!file.exists(opt$cumulative))
    stop("[who_monthly] cumulative snapshots not found: ", opt$cumulative,
         " (pass --cumulative \"\" for the weekly series alone)")
  snap <- read.csv(opt$cumulative, stringsAsFactors = FALSE)
  mon_cum <- whoCumulativeToMonthly(snap, zero_iso3 = zero_iso3,
                                    max_interval_days = opt$max_interval_days)
  message(sprintf("[who_monthly] cumulative: %d snapshots, %d countries, cuts %s to %s -> %d country-months (zero allowlist: %s)",
                  nrow(snap), n_distinct(snap$iso3), min(as.Date(snap$table_as_of)),
                  max(as.Date(snap$table_as_of)), nrow(mon_cum),
                  if (length(zero_iso3)) paste(zero_iso3, collapse = ",") else "none"))
}
mon <- whoMergeMonthly(mon_wk, mon_cum)

summ <- mon %>%
  group_by(iso3) %>%
  summarise(first_complete = suppressWarnings(min(t_start[!partial])),
            last_complete  = suppressWarnings(max(t_start[!partial])),
            n_complete = sum(!partial), n_partial = sum(partial),
            n_filled = sum(source == "sitrep_cumulative"),
            n_filled_complete = sum(source == "sitrep_cumulative" & !partial),
            cases_complete = round(sum(cases[!partial])), .groups = "drop")
print(as.data.frame(summ), row.names = FALSE)
message(sprintf("[who_monthly] %d weeks, %d countries -> %d country-months (%d complete; %d filled from the cumulative updates, %d of them complete)",
                nrow(wk), n_distinct(mon$iso3), nrow(mon), sum(!mon$partial),
                sum(mon$source == "sitrep_cumulative"),
                sum(mon$source == "sitrep_cumulative" & !mon$partial)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- makeWhoMonthlyFile(opt$outdir)

# Provenance sidecar: the two inputs and the one judgement call (which
# countries' unchanged cumulatives are zeros) would otherwise leave no trace.
# Built before the table is written (same fields and values).
meta <- list(
  generator = "analysis/03_prepare_who_monthly.R",
  written   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  inputs    = list(
    weekly     = list(path = opt$input,
                      mtime = format(file.mtime(opt$input), "%Y-%m-%dT%H:%M:%S"),
                      n_rows = nrow(wk)),
    cumulative = if (use_cum) list(path = opt$cumulative,
                                   mtime = format(file.mtime(opt$cumulative), "%Y-%m-%dT%H:%M:%S"),
                                   n_rows = nrow(snap)) else NULL),
  zero_iso3 = as.list(zero_iso3),
  max_interval_days = opt$max_interval_days,
  counts    = list(country_months = nrow(mon), complete = sum(!mon$partial),
                   by_source = as.list(table(mon$source)))
)
meta_path <- sub("\\.rds$", "_meta.yml", out)

# Overwrite guard (the rule of 03_prepare_occ.R, on the table's
# own record): an existing table whose meta -- less `written` and the inputs'
# mtimes -- or content differs (a refreshed WHO feed, another zero allowlist),
# or that has no meta, is not replaced without --force. An identical rebuild
# passes. utils/provenance.R::whoMonthlyGuardStamp().
old_who <- if (file.exists(out) && file.exists(meta_path))
  whoMonthlyGuardStamp(yaml::read_yaml(meta_path), readRDS(out)) else NULL
assertPrepStampOverwrite(out, whoMonthlyGuardStamp(meta, mon), force = isTRUE(opt$force),
                         rerun_hint = "Rscript analysis/03_prepare_who_monthly.R ... --force",
                         old_stamp = old_who)

saveRDS(mon, out)
message("[who_monthly] wrote ", out)
yaml::write_yaml(meta, meta_path)
message("[who_monthly] wrote ", meta_path)
