# test_bpfilter_traj_compact.R -- simulator.filter.traj_storage: the compact
# filter_traj store of run_one_bpfilter_replicate() (spatpomp_core.R::
# .bpfilter_traj_compact(), a copy of spatPomp 1.1.0's bpfilter.internal with
# ancestry pruning and optional float32 storage).
#
# HARD gate (see analysis/tests/run_all.R). Synthetic, no MCMC, ~1 min.
#
#   C  config: resolve_filter_cfg() defaults (full, 1) and refusals; the
#      smoothing config threads both keys; both callers pass them on.
#   G  golden: "full" reproduces the trajectory and logLik that
#      spatPomp::bpfilter(filter_traj = TRUE) gave before the change
#      (analysis/tests/golden/bpfilter_traj_full_golden.rds).
#   P  "pruned" is identical() to "full" (trajectory and logLik, same seed)
#      for prune_every 1, 3 and > ntimes, on country and ADM1 blocks and on
#      the duration + monthly measurement; pruning really drops columns.
#   R  regression: with few particles over a long window a level is
#      compacted twice while the level above still holds dead columns whose
#      parent is already 0; those entries must keep their place in the remap,
#      or "pruned" walks wrong lineages silently. "pruned" identical() to
#      "full" on 90 such runs.
#   F  "pruned_float32": logLik identical, integer-valued states exact,
#      every other state within 1e-6 relative, same dims / dimnames.
#
# The fixture's onset hazard reads each particle's own lag (ongoing_outbreak
# _lag1), so the particle weights differ and the resampling coalesces the
# ancestry; with a state-independent hazard every weight is equal, the
# resampling is the identity, and a pruning bug would go unseen.
#
# Run:     Rscript analysis/tests/test_bpfilter_traj_compact.R
# Capture: Rscript analysis/tests/test_bpfilter_traj_compact.R --capture
#          (only from code whose "full" path is spatPomp::bpfilter unchanged)

suppressPackageStartupMessages({
  library(tibble); library(dplyr); library(lubridate); library(tidyr)
})
source("analysis/utils/smoothing_core.R")   # sources simulator_core.R, spatpomp_core{,_c}.R
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

GOLDEN_PATH <- "analysis/tests/golden/bpfilter_traj_full_golden.rds"
CAPTURE     <- "--capture" %in% commandArgs(trailingOnly = TRUE)

cat("# test_bpfilter_traj_compact\n\n")
is_err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# ===========================================================================
# Fixture: n_ctry countries x n_adm1 ADM1, Bernoulli obs at 20 %
# ===========================================================================
make_fixture <- function(n_ctry = 3L, n_adm1 = 3L, n_months = 30L,
                         measurement = "occurrence_bernoulli", seed = 1L) {
  ctry <- sprintf("C%02d", seq_len(n_ctry))
  frame_tu <- tidyr::expand_grid(country = ctry, adm1n = seq_len(n_adm1)) %>%
    dplyr::mutate(gid = sprintf("%s_%d_1", country, adm1n),
                  adm1_unit = sprintf("%s_%d", country, adm1n),
                  AFRO_region = "REG_W", tu = "hmmW", min_seas = 1L, pop = 1e5) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
  months <- as.Date("2015-01-01") %m+% months(0:(n_months - 1L))
  set.seed(seed)
  panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::left_join(frame_tu %>% dplyr::select(gid, min_seas), by = "gid") %>%
    dplyr::mutate(country = sub("_[0-9]+_[0-9]+$", "", gid),
                  t_start_seas = as.integer(lubridate::yday(t_start) - 1L),
                  mean_precip = 0, dev_precip = 0, mean_temp_daymean = 0, dev_temp_daymean = 0,
                  outbreak = stats::rbinom(dplyr::n(), 1L, 0.2),
                  time_to_last_outbreak_cat = "[6,Inf]",
                  ongoing_outbreak_lag1 = 0L, ongoing_outbreak_country_lag1 = 0L,
                  ongoing_outbreak_tu_lag1 = 0L) %>%
    tibble::as_tibble()
  window <- list(start = format(min(months)), end = format(max(months) %m+% months(1L)))
  # Per-day hazard: monthly onset p ~ 0.18 cold, ~ 0.8 after an active month.
  d_occ <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c(
    "(Intercept)" = -5, "ongoing_outbreak_lag1" = 2, "ongoing_outbreak_country_lag1" = 0.7))
  extra <- list()
  if (!identical(measurement, "occurrence_bernoulli")) {
    d_dur <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -3),
                                         active_terms = "(Intercept)")
    d_dur$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
    d_mon <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -9),
                                         active_terms = "(Intercept)", phi = 0.45)
    d_mon$offset_pop <- TRUE; d_mon$offset_logoffset <- TRUE
    extra <- list(draws_dur = d_dur, draws_mon = d_mon)
  }
  obj <- do.call(build_spatpomp_tu_c, c(list(frame_tu, panel_tu, d_occ, window,
    n_draws_in_filter = 1L, fixed_draw_idx = 1L, measurement = measurement), extra))
  list(obj = obj, frame_tu = frame_tu)
}

fx_occ <- make_fixture()
fx_mon <- make_fixture(n_ctry = 2L, n_adm1 = 2L, n_months = 24L,
                       measurement = "occurrence_duration_mon_NB", seed = 2L)
CASES <- list(
  occ_country = list(fx = fx_occ, gran = "country", Np = 300L, seed = 11L),
  occ_adm1    = list(fx = fx_occ, gran = "adm1",    Np = 200L, seed = 12L),
  mon_country = list(fx = fx_mon, gran = "country", Np = 200L, seed = 13L))
blocks_of <- function(cs) make_block_list(cs$fx$obj, cs$fx$frame_tu, granularity = cs$gran)

# ===========================================================================
# G -- golden of the pre-change "full" path
# ===========================================================================
if (CAPTURE) {
  # The pre-change run_one_bpfilter_replicate() body, verbatim.
  gold <- lapply(CASES, function(cs) {
    set.seed(cs$seed)
    bpf <- spatPomp::bpfilter(cs$fx$obj, Np = cs$Np, block_list = blocks_of(cs),
                              filter_traj = TRUE)
    list(filter_traj = bpf@filter.traj, logLik = pomp::logLik(bpf))
  })
  dir.create(dirname(GOLDEN_PATH), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gold, GOLDEN_PATH)
  cat("captured", length(gold), "cases ->", GOLDEN_PATH, "\n")
  quit(status = 0L)
}
if (!file.exists(GOLDEN_PATH)) stop("golden missing: run with --capture first")
gold <- readRDS(GOLDEN_PATH)

run_mode <- function(cs, mode, k = 1L)
  run_one_bpfilter_replicate(cs$fx$obj, Np = cs$Np, block_list = blocks_of(cs),
                             seed = cs$seed, traj_storage = mode, prune_every = k)

# ===========================================================================
# C -- config
# ===========================================================================
cat("## C: config\n")
fc0 <- resolve_filter_cfg(list())
expect(identical(fc0$traj_storage, "pruned") && identical(fc0$traj_prune_every, 1L),
       "C1 resolve_filter_cfg() defaults: traj_storage pruned, traj_prune_every 1")
d_flt <- get_default_options()$simulator$filter
expect(identical(d_flt$traj_storage, "pruned") && identical(d_flt$traj_prune_every, 1L),
       "C2 get_default_options() carries the same defaults")
fc1 <- resolve_filter_cfg(list(traj_storage = "pruned_float32", traj_prune_every = 4))
expect(identical(fc1$traj_storage, "pruned_float32") && identical(fc1$traj_prune_every, 4L),
       "C3 resolve_filter_cfg() passes both keys through")
expect(is_err(resolve_filter_cfg(list(traj_storage = "compact"))) &&
       is_err(resolve_filter_cfg(list(traj_prune_every = 0))) &&
       is_err(resolve_filter_cfg(list(traj_prune_every = c(1, 2)))),
       "C4 an unknown traj_storage and a non-positive / vector traj_prune_every are refused")
expect(identical(as.character(utils::packageVersion("spatPomp")),
                 .SP_TRAJ_COMPACT_SPATPOMP_VERSION),
       "C5 installed spatPomp is the version .bpfilter_traj_compact() was copied from")
cfg_path <- tempfile(fileext = ".yml")
yaml::write_yaml(list(simulator = list(backend = "spatpomp", filter = list(
  measurement = "occurrence_bernoulli", traj_storage = "full", traj_prune_every = 3L))), cfg_path)
S1 <- smoothing_make_config(list(config = cfg_path))
yaml::write_yaml(list(simulator = list(backend = "spatpomp", filter = list(
  measurement = "occurrence_bernoulli"))), cfg_path)
S0 <- smoothing_make_config(list(config = cfg_path))
expect(identical(S1$TRAJ_STORAGE, "full") && identical(S1$TRAJ_PRUNE_EVERY, 3L) &&
       identical(S0$TRAJ_STORAGE, "pruned") && identical(S0$TRAJ_PRUNE_EVERY, 1L) &&
       identical(smoothing_default_config()$TRAJ_STORAGE, "pruned"),
       "C6 smoothing_make_config() threads the keys; a config without them gets pruned / 1")
call_passes <- function(path, fn = NULL) {
  txt <- if (is.null(fn)) readLines(path) else deparse(get(fn))
  txt <- paste(txt, collapse = " ")
  grepl("run_one_bpfilter_replicate\\([^)]*traj_storage\\s*=[^)]*prune_every\\s*=", txt)
}
expect(call_passes(NULL, "smoothing_filter_one") &&
       call_passes("analysis/07b_simulator_spatpomp.R"),
       "C7 smoothing_filter_one() and 07b pass traj_storage / prune_every to the filter")

# ===========================================================================
# G / P / F
# ===========================================================================
for (nm in names(CASES)) {
  cs <- CASES[[nm]]
  cat(sprintf("\n## %s (Np = %d, %s blocks)\n", nm, cs$Np, cs$gran))
  full <- run_mode(cs, "full")
  expect(identical(full, gold[[nm]]),
         sprintf("G1 %s: full == pre-change golden (trajectory and logLik)", nm))
  ntimes <- dim(full$filter_traj)[3] - 1L
  for (k in c(1L, 3L, ntimes + 5L)) {
    expect(identical(run_mode(cs, "pruned", k), full),
           sprintf("P1 %s: pruned (prune_every = %d) identical() to full", nm, k))
  }
  # Pruning really drops columns (the identity above is not vacuous).
  e <- new.env()
  set.seed(cs$seed)
  invisible(.bpfilter_traj_compact(cs$fx$obj, Np = cs$Np,
                                   block_list = lapply(blocks_of(cs), as.integer),
                                   diag_env = e, verbose = FALSE))
  kc <- e$kept_cols
  expect(all(kc[nrow(kc), ] == cs$Np) && sum(kc) < 0.5 * length(kc) * cs$Np,
         sprintf("P2 %s: top level full, store < half of full (kept %.1f %%)",
                 nm, 100 * sum(kc) / (length(kc) * cs$Np)))

  f32 <- run_mode(cs, "pruned_float32")
  a <- full$filter_traj; b <- f32$filter_traj
  is_int <- a == round(a)
  rel <- abs(a - b) / pmax(abs(a), .Machine$double.xmin)
  expect(identical(f32$logLik, full$logLik), sprintf("F1 %s: float32 logLik identical", nm))
  expect(identical(dim(a), dim(b)) && identical(dimnames(a), dimnames(b)),
         sprintf("F2 %s: float32 trajectory dims / dimnames identical", nm))
  expect(identical(a[is_int], b[is_int]),
         sprintf("F3 %s: float32 integer-valued states exact (%d entries)", nm, sum(is_int)))
  expect(any(!is_int) && max(rel[!is_int]) < 1e-6,
         sprintf("F4 %s: float32 non-integer states within 1e-6 relative (max %.2e, %d entries)",
                 nm, max(rel[!is_int]), sum(!is_int)))
  for (k in c(3L)) {
    expect(identical(run_mode(cs, "pruned_float32", k), f32),
           sprintf("F5 %s: float32 result does not depend on prune_every", nm))
  }
}

# ===========================================================================
# R -- regression: repeated compaction under heavy coalescence
# ===========================================================================
cat("\n## R: few particles, long window\n")
fx_long <- make_fixture(n_ctry = 3L, n_adm1 = 3L, n_months = 120L, seed = 3L)
for (gran in c("country", "adm1")) {
  bl <- make_block_list(fx_long$obj, fx_long$frame_tu, granularity = gran)
  n_diff <- 0L; n_run <- 0L
  for (Np in c(5L, 10L, 20L)) for (sd in 1:15) {
    full_r <- run_one_bpfilter_replicate(fx_long$obj, Np = Np, block_list = bl, seed = sd,
                                         traj_storage = "full")
    pr_r <- tryCatch(run_one_bpfilter_replicate(fx_long$obj, Np = Np, block_list = bl, seed = sd,
                                                traj_storage = "pruned"),
                     error = function(e) NULL)
    n_run <- n_run + 1L
    if (!identical(pr_r, full_r)) n_diff <- n_diff + 1L
  }
  expect(n_diff == 0L,
         sprintf("R1 %s blocks: pruned identical() to full on %d runs (Np 5/10/20, 120 months; %d differ)",
                 gran, n_run, n_diff))
}

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
