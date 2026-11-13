# test_sim_metrics_golden.R -- golden characterization of the simulator
# validation metrics (07c_simulator_postprocess.R / 07d's compute_hindcast_series).
#
# HARD gate (see analysis/tests/run_all.R). The pure metric functions
# (utils/sim_metrics.R), the partition IO (utils/sim_io.R) and the shared CRPS /
# geography / hindcast helpers must leave every metric bit-identical, so
# this test runs ONE deterministic synthetic world through:
#   * run_validation() end to end (every compute_* metric, the panel-mask
#     fallback, attach_episode_case_totals(), write_validation_summary()) under
#     both total-case bases (exact_gid; region_key with an LP-keyed mask);
#   * aggregate_partitions() / combine_group_runs() on a temp partition tree;
#   * compute_hindcast_series() on the merged run;
#   * a few direct calls (load_observed_from_panel() with a horizon,
#     compute_pit_binom() seed, reconstruct_lag_series() with active_col);
# and compares the captured list with identical() against
# analysis/tests/golden/sim_metrics_golden.rds.
#
#   Rscript analysis/tests/test_sim_metrics_golden.R --capture   # once, from a reference version
#   Rscript analysis/tests/test_sim_metrics_golden.R             # compare
#
# Synthetic data only (temp dirs) -- no fit, no generated_data, no network.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tibble); library(tidyr)
  library(lubridate); library(purrr); library(qs2)
})
cli <- parse_args(OptionParser(option_list = list(
  make_option("--capture", action = "store_true", default = FALSE,
              help = "write the golden file instead of comparing"))))

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

# The 07c closure (07c's own source() calls sit behind its sys.nframe() guard).
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/formula.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/dependence.R")
source("analysis/utils/figures_io.R")
source("analysis/utils/anchor.R")
for (.f in c("analysis/utils/obs_sim_compare.R", "analysis/utils/sim_metrics.R"))
  if (file.exists(.f)) source(.f)
source("analysis/07c_simulator_postprocess.R")
source("analysis/07d_simulator_figures.R")  # compute_hindcast_series() (CLI guarded)

GOLDEN_PATH <- "analysis/tests/golden/sim_metrics_golden.rds"
cat("# simulator validation metrics golden (mode:", if (cli$capture) "capture" else "compare", ")\n")

# ---------------------------------------------------------------------------
# Synthetic world
# ---------------------------------------------------------------------------
set.seed(20261006)
frame <- expand.grid(tu = c("hmmW", "hmmE"), country_i = 1:2, adm1n = 1:2, adm2n = 1:2,
                     stringsAsFactors = FALSE) %>%
  mutate(country     = c(hmmW = "NGA", hmmE = "KEN")[tu],
         country     = ifelse(country_i == 2, c(NGA = "GHA", KEN = "MWI")[country], country),
         gid         = sprintf("%s_%d_%d", country, adm1n, adm2n),
         adm1_unit   = sprintf("%s_%d", country, adm1n),
         AFRO_region = ifelse(tu == "hmmW", "WEST", "EAST"),
         min_seas    = 1L, pop = 1e5) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  as_tibble()

months <- seq(as.Date("2017-01-01"), as.Date("2018-12-01"), by = "month")
n_rep  <- 12L

# Kernel-like episodes: start w.p. 0.12 when inactive, persist w.p. 0.6.
sim_one <- function(r, g) {
  act <- integer(length(months)); y <- integer(length(months))
  for (m in seq_along(months)) {
    prev <- if (m == 1L) 0L else act[m - 1L]
    if (prev == 1L) act[m] <- as.integer(runif(1) < 0.6)
    else if (runif(1) < 0.12) { act[m] <- 1L; y[m] <- 1L }
  }
  tibble(replicate_id = r, gid = g, month = months, Y = y, active = act,
         cases_month = ifelse(act == 1L, as.numeric(rpois(length(months), 40)), NA_real_))
}
ms <- bind_rows(lapply(seq_len(n_rep), function(r) bind_rows(lapply(frame$gid, function(g) sim_one(r, g)))))

# Event table: one row per (replicate, gid, episode); censored if active in the last month.
episodes <- function(d) {
  rl <- rle(d$active); ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1L
  keep <- rl$values == 1L
  if (!any(keep)) return(NULL)
  tibble(replicate_id = d$replicate_id[1], gid = d$gid[1],
         start_month = format(d$month[starts[keep]]),
         duration_months = rl$lengths[keep],
         censored = as.integer(ends[keep] == nrow(d)))
}
et <- ms %>% group_by(replicate_id, gid) %>% group_split() %>% lapply(episodes) %>% bind_rows() %>%
  left_join(frame %>% select(gid, country), by = "gid")

# Observed survival panel: at-risk months Mar-2017..Oct-2018 for all but one gid,
# plus a coarse ADM1 report (gid NGA_1) feeding the region_key total-case basis.
obs_gids <- setdiff(frame$gid, "MWI_2_2")
obs_months <- months[months >= as.Date("2017-03-01") & months <= as.Date("2018-10-01")]
panel <- expand.grid(gid = obs_gids, t_start = obs_months, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  mutate(outbreak = as.integer(runif(n()) < 0.15),
         location_period_id = paste0("lp_", gid))
panel <- bind_rows(panel, tibble(gid = "NGA_1", t_start = as.Date("2018-02-01"),
                                 outbreak = 1L, location_period_id = "lp_NGA_1"))

obs_meta <- panel %>% filter(outbreak == 1L) %>%
  group_by(location_period_id, gid) %>%
  summarise(TL = min(t_start), .groups = "drop") %>%
  mutate(suspected_cases = as.numeric(rpois(n(), 300)),
         country = sub("_.*$", "", gid),
         AFRO_region = ifelse(country %in% c("NGA", "GHA"), "WEST", "EAST"),
         admin_level = ifelse(gid == "NGA_1", "ADM1", "ADM2")) %>%
  select(location_period_id, suspected_cases, TL, country, AFRO_region, admin_level)

obs_dur <- panel %>% filter(outbreak == 1L, gid != "NGA_1") %>%
  mutate(episode_id = paste0(gid, "_", format(t_start)),
         country = sub("_.*$", "", gid), month = format(t_start),
         duration_months = 1L + (seq_len(n()) %% 4L), censored = as.integer(seq_len(n()) %% 7L == 0L)) %>%
  select(episode_id, gid, country, month, duration_months, censored)

mon_data <- expand.grid(gid = obs_gids, month = months, stringsAsFactors = FALSE) %>% as_tibble() %>%
  mutate(series_end = as.Date("2018-08-01"),
         cases_month = as.numeric(rpois(n(), 20)))

tmp <- tempfile("simgolden_"); dir.create(tmp)
merged_dir <- file.path(tmp, "merged"); dir.create(merged_dir)
sim_io_save(ms, file.path(merged_dir, "monthly_summary.qs"), "qs")
sim_io_save(et, file.path(merged_dir, "event_table.qs"), "qs")
paths <- list(panel = file.path(tmp, "panel.rds"), meta = file.path(tmp, "meta.rds"),
              dur = file.path(tmp, "dur.rds"), mon = file.path(tmp, "mon.rds"),
              dep = file.path(tmp, "dep.rds"))
saveRDS(panel, paths$panel); saveRDS(obs_meta, paths$meta); saveRDS(obs_dur, paths$dur)
saveRDS(mon_data, paths$mon); saveRDS(list(tau = list(overall = list(tau = 0.31))), paths$dep)

# ---------------------------------------------------------------------------
# Capture
# ---------------------------------------------------------------------------
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
out <- list()

for (basis in c("exact_gid", "region_key")) {
  vdir <- file.path(tmp, paste0("valid_", basis))
  set.seed(1L)
  out[[paste0("validation_", basis)]] <- quiet(run_validation(
    merged_dir = merged_dir, obs_panel_path = paths$panel, frame = frame,
    valid_dir = vdir, max_reps = 8L, obs_metadata_path = paths$meta,
    total_case_basis = basis, obs_dur_panel_path = paths$dur,
    dependence_path = paths$dep,
    monthly_data_path = if (basis == "region_key") paths$mon else NULL))
  out[[paste0("report_", basis)]] <-   # the report's "Generated: <Sys.time()>" line is masked
    sub("^(\\s*Generated: ).*$", "\\1<TIME>", readLines(file.path(vdir, "validation_report.txt")))
}

# metric_set = "hindcast" (07c --hindcast): every metric 09 reads is identical to
# the full set's; the rest stay NULL. Checked against the full-set run above
# rather than stored, so the golden file is unchanged.
hc_read  <- c("pit", "monthly_freq", "crps", "crps_total_cases")
hc_skip  <- c("crps_count_by_geo", "ar_metrics", "contagion", "active_fraction",
              "duration_dist", "sim_kendall_tau", "obs_kendall_tau")
for (basis in c("exact_gid", "region_key")) {
  set.seed(1L)
  hc <- quiet(run_validation(
    merged_dir = merged_dir, obs_panel_path = paths$panel, frame = frame,
    valid_dir = file.path(tmp, paste0("valid_hc_", basis)), max_reps = 8L,
    obs_metadata_path = paths$meta, total_case_basis = basis,
    obs_dur_panel_path = paths$dur, dependence_path = paths$dep,
    monthly_data_path = if (basis == "region_key") paths$mon else NULL,
    metric_set = "hindcast"))
  fullm <- out[[paste0("validation_", basis)]]
  expect(identical(names(hc), names(fullm)), sprintf("hindcast set keeps the metrics list shape (%s)", basis))
  for (k in hc_read)
    expect(identical(hc[[k]], fullm[[k]]), sprintf("hindcast set: %s identical to the full set (%s)", k, basis))
  expect(all(vapply(hc[hc_skip], is.null, logical(1))), sprintf("hindcast set: unread metrics NULL (%s)", basis))
}
# no observed metadata (the clusters' case): no total-case CRPS, event table not read
hc_nometa <- quiet(run_validation(
  merged_dir = merged_dir, obs_panel_path = paths$panel, frame = frame,
  valid_dir = file.path(tmp, "valid_hc_nometa"), max_reps = 8L,
  obs_metadata_path = file.path(tmp, "absent_meta.rds"), metric_set = "hindcast"))
expect(is.null(hc_nometa$crps_total_cases) &&
         all(vapply(c("pit", "monthly_freq", "crps"), function(k)
           identical(hc_nometa[[k]], out$validation_exact_gid[[k]]), logical(1))),
       "hindcast set without metadata: pit / monthly_freq / crps unchanged, no total-case CRPS")

# Partition tree: 3 TUs x 2 chunks, one partition without a sentinel (skipped).
parts <- file.path(tmp, "parts")
k <- 0L
for (tu in unique(frame$tu)) for (chunk in 1:2) {
  k <- k + 1L
  pd <- file.path(parts, paste0("tu=", tu), paste0("chunk=", chunk)); dir.create(pd, recursive = TRUE)
  reps <- if (chunk == 1L) 1:6 else 7:12
  gids <- frame$gid[frame$tu == tu]
  qs_save(et %>% filter(replicate_id %in% reps, gid %in% gids), file.path(pd, "event_table.qs"))
  qs_save(ms %>% filter(replicate_id %in% reps, gid %in% gids), file.path(pd, "monthly_summary.qs"))
  if (k != 3L) writeLines("done", file.path(pd, "_SUCCESS"))
}
agg_dir <- file.path(tmp, "agg")
out$aggregate_all <- quiet(aggregate_partitions(parts, agg_dir))
out$aggregate_tu  <- quiet(aggregate_partitions(parts, file.path(tmp, "agg_tu"), tus = "hmmW"))

grp_dirs <- c(enso_west = agg_dir, enso_east = file.path(tmp, "agg_tu"))
combined <- file.path(tmp, "combined")
out$combine <- quiet(combine_group_runs(grp_dirs, combined, cutoff_month = "2018-06-01"))

out$hindcast <- quiet(compute_hindcast_series(
  merged_dir, paths$panel, monthly_data_path = paths$mon,
  launch = "2018-01-01", horizon_end = "2018-09-01",
  group_label = "enso_west", launch_year = 2018L))
out$hindcast_no_cases <- quiet(compute_hindcast_series(
  merged_dir, paths$panel, launch = "2017-06-01", horizon_end = "2017-12-01"))

out$obs_horizon <- load_observed_from_panel(panel, gids = frame$gid,
                                            horizon = list(start = "2017-06-01", end = "2018-03-01"))
out$pit_seed7   <- compute_pit_binom(ms, load_observed_from_panel(panel), seed = 7L)
out$lags_active <- quiet(reconstruct_lag_series(ms, frame, max_reps = 3L, active_col = "active"))
out$km          <- .km_estimate(c(1, 2, 2, 3, 5, 8), c(1, 1, 0, 1, 0, 1))

# Paths written into results (temp dirs) differ between runs -> strip them.
strip_paths <- function(x) {
  if (is.list(x) && !is.data.frame(x)) return(lapply(x, strip_paths))
  if (is.character(x)) return(gsub(tmp, "<TMP>", x, fixed = TRUE))
  x
}
out <- strip_paths(out)

expect(length(out) >= 12L && !is.null(out$validation_exact_gid$crps_total_cases) &&
         !is.null(out$validation_region_key$crps_total_cases) &&
         !is.null(out$validation_exact_gid$duration_dist) &&
         !is.null(out$validation_exact_gid$active_fraction),
       "fixture exercises every metric family (total-case CRPS, duration, active fraction)")

if (cli$capture) {
  dir.create(dirname(GOLDEN_PATH), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, GOLDEN_PATH)
  cat("  captured", length(out), "entries ->", GOLDEN_PATH, "\n")
} else {
  expect(file.exists(GOLDEN_PATH), paste("golden file exists:", GOLDEN_PATH,
                                         "(run once with --capture on the pre-refactor code)"))
  if (file.exists(GOLDEN_PATH)) {
    gold <- readRDS(GOLDEN_PATH)
    expect(identical(sort(names(gold)), sort(names(out))), "same captured entries as the golden")
    for (nm in intersect(names(gold), names(out))) {
      ok <- identical(gold[[nm]], out[[nm]])
      expect(ok, sprintf("identical to golden: %s", nm))
      if (!ok) print(all.equal(gold[[nm]], out[[nm]]))
    }
  }
}

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
