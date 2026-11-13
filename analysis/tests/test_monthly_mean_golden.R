# test_monthly_mean_golden.R -- golden characterization of every site that
# assembles and draws the monthly (`mon`) case-count mean.
#
# HARD gate (see analysis/tests/run_all.R). Built for the monthly-incidence-cap
# refactor, which routes the four inline copies of the monthly mean recipe
# (evaluate_lp() + mso basis lookup + log(pop) + log_offset/onset correction +
# base-rate anchor) through simulator_core.R::monthlyLogMu(), and the two NB
# draws through drawMonthlyCases(). With no cap configured every site must
# produce the same numbers as before. (The rk_* sites of the removed R filter
# kernel, make_step_fn(), were dropped; the sm_* sites run on C-kernel
# bpfilter trajectories.)
#
#   v1_ungrouped / v1_grouped   simulate_one_tu(): mu_mon + cases_month
#                               (v1 kernel AND the spatPomp forecast phase)
#   sm_ungrouped / sm_grouped   extract_smoothing_one_tu() on a C-kernel
#                               bpfilter trajectory: cases_mon_mat (the
#                               smoothing's post-hoc draw)
#
# Every fixture turns on all three offsets (pop, log_offset with the onset
# correction, a base-rate anchor), a latent mso cross-basis, nonzero RE
# intercepts and a spillover slope, so each term of the recipe moves mu.
# No MCMC.
#
#   Rscript analysis/tests/test_monthly_mean_golden.R --capture   # once, pre-refactor
#   Rscript analysis/tests/test_monthly_mean_golden.R             # compare

suppressPackageStartupMessages({
  library(optparse); library(tibble); library(dplyr); library(lubridate); library(tidyr)
})
cli <- parse_args(OptionParser(option_list = list(
  make_option("--capture", action = "store_true", default = FALSE,
              help = "write the golden file instead of comparing"))))

source("analysis/utils/smoothing_core.R")   # sources simulator_core.R, spatpomp_core{,_c}.R, dlnm.R
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

GOLDEN_PATH <- "analysis/tests/golden/monthly_mean_golden.rds"
cat("# monthly mean golden (mode:", if (cli$capture) "capture" else "compare", ")\n")

# ---------------------------------------------------------------------------
# Fixture: 2 countries x 2 ADM1 x 2 ADM2, unequal populations, 2 groups.
# ---------------------------------------------------------------------------
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2, adm2n = 1:2) %>%
  dplyr::mutate(
    gid           = sprintf("%s_%d_%d", country, adm1n, adm2n),
    adm1_unit     = sprintf("%s_%d", country, adm1n),
    AFRO_region   = "REG_W",
    tu            = "hmmW",
    min_seas      = 1L,
    pop           = c(2e4, 5e4, 1e5, 3e5, 8e4, 1.5e5, 4e4, 6e5),
    base_rate_mon = c(1.2, 0.8, 1.0, 1.5, 0.9, 1.1, 0.7, 1.3),
    group         = dplyr::if_else(country == "AAA", "ga", "gb")
  ) %>%
  dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit,
                base_rate_mon, group) %>%
  tibble::as_tibble()

months  <- as.Date("2024-01-01") %m+% months(0:11)
horizon <- list(start = "2024-01-01", end = "2024-12-31")
window  <- list(start = "2024-01-01", end = "2025-01-01")
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
  dplyr::left_join(frame_tu %>% dplyr::select(gid, min_seas), by = "gid") %>%
  dplyr::mutate(
    country           = sub("_[0-9]+_[0-9]+$", "", gid),
    t_start_seas      = as.integer(lubridate::yday(t_start) - 1L),
    mean_precip       = 0.0, dev_precip = 0.0,
    mean_temp_daymean = 0.0, dev_temp_daymean = 0.0,
    outbreak          = 0L,
    time_to_last_outbreak_cat     = "[6,Inf]",
    ongoing_outbreak_lag1         = 0L,
    ongoing_outbreak_country_lag1 = 0L,
    ongoing_outbreak_tu_lag1      = 0L
  ) %>%
  tibble::as_tibble()

cb_meta_mon <- list(mso = list(single_lag = 0L, range = c(0, 8),
                               onebasis_attrs = list(fun = "lin"),
                               cb_col_names = "mso_cb1"))

zm <- function(cols, S, vals = 0) {
  m <- matrix(vals, nrow = S, ncol = length(cols), byrow = TRUE); colnames(m) <- cols; m
}
# Bundle with explicit FE values and (optionally) seeded RE intercepts.
# `drop_region` = TRUE keys country/ADM1 as the grouped (region-dropped) fits do.
make_bundle <- function(frame_g, fe_vals, S = 2L, re_sd = 0, seed = 1L,
                        drop_region = FALSE, basehaz = list(kind = "const", const = 0),
                        phi = 1.0) {
  set.seed(seed)
  fe <- zm(names(fe_vals), S, unname(fe_vals))
  ctry <- unique(as.character(frame_g$country))
  adm1 <- if (drop_region) unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))
          else unique(frame_g$adm1_unit)
  list(S = S, fe = fe,
       re_intercept = list(
         region  = if (drop_region) NULL else zm(unique(frame_g$AFRO_region), S),
         country = zm(ctry, S, stats::rnorm(length(ctry), 0, re_sd)),
         adm1    = zm(adm1, S, stats::rnorm(length(adm1), 0, re_sd)),
         adm2    = NULL),
       re_slope = list(), basehaz = basehaz, phi = rep(phi, S),
       active_terms = names(fe_vals))
}
occ_fe <- c("(Intercept)" = log(0.006), "ongoing_outbreak_lag1" = 0.3)
dur_fe <- c("(Intercept)" = log(0.008))
mon_fe <- c("(Intercept)" = log(4e-7), "mso_cb1" = -0.12, "ongoing_outbreak_lag1" = 0.2)
mon_flags <- function(b) {
  b$offset_pop <- TRUE; b$offset_logoffset <- TRUE
  b$offset_baserate <- TRUE; b$anchor_col <- "base_rate_mon"
  b
}
dur_bh <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)

# Ungrouped bundles.
d_occ <- make_bundle(frame_tu, occ_fe, re_sd = 0.4, seed = 11L)
d_dur <- make_bundle(frame_tu, dur_fe, re_sd = 0.2, seed = 12L, basehaz = dur_bh)
d_mon <- mon_flags(make_bundle(frame_tu, mon_fe, re_sd = 0.6, seed = 13L, phi = 0.45))

# Grouped bundles (region dropped), distinct per group.
fr <- split(frame_tu, frame_tu$group)
g_occ <- list(ga = make_bundle(fr$ga, occ_fe, re_sd = 0.4, seed = 21L, drop_region = TRUE),
              gb = make_bundle(fr$gb, occ_fe, re_sd = 0.4, seed = 22L, drop_region = TRUE))
g_dur <- list(ga = make_bundle(fr$ga, dur_fe, re_sd = 0.2, seed = 23L, drop_region = TRUE, basehaz = dur_bh),
              gb = make_bundle(fr$gb, dur_fe, re_sd = 0.2, seed = 24L, drop_region = TRUE, basehaz = dur_bh))
g_mon <- list(ga = mon_flags(make_bundle(fr$ga, mon_fe, re_sd = 0.6, seed = 25L,
                                         drop_region = TRUE, phi = 0.45)),
              gb = mon_flags(make_bundle(fr$gb, c(mon_fe[1] + 0.7, mon_fe[-1]), re_sd = 0.6,
                                         seed = 26L, drop_region = TRUE, phi = 0.6)))
grouping <- list(
  draws_occ_by_group   = g_occ,
  draws_dur_by_group   = g_dur,
  draws_mon_by_group   = g_mon,
  cb_meta_occ_by_group = list(ga = NULL, gb = NULL),
  cb_meta_dur_by_group = list(ga = NULL, gb = NULL),
  cb_meta_mon_by_group = list(ga = cb_meta_mon, gb = cb_meta_mon)
)

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------
run_v1 <- function(grouped, ...) {
  res <- simulate_one_tu(
    tu = "hmmW", frame_tu = frame_tu, panel_tu = panel_tu,
    draws_occ = if (grouped) g_occ$ga else d_occ,
    draws_dur = if (grouped) g_dur$ga else d_dur,
    draws_mon = if (grouped) g_mon$ga else d_mon,
    cb_meta_mon = if (grouped) NULL else cb_meta_mon,
    grouping = if (grouped) grouping else NULL,
    replicate_ids = 1:4, posterior_idx_occ = c(1L, 2L, 1L, 2L),
    seed = 4242L, output_part = tempfile("mmg_"), horizon = horizon,
    dry_run = TRUE, ...)
  ms <- res$monthly_summary
  ms[, c("replicate_id", "gid", "month", "active", "mu_mon", "cases_month")]
}

# Smoothing: one C-kernel bpfilter replicate per draw (as smoothing_filter_one()
# does), then the post-hoc monthly draw.
run_sm <- function(grouped, ...) {
  results <- lapply(1:2, function(r) {
    obj <- build_spatpomp_tu_c(
      frame_tu, panel_tu, if (grouped) g_occ$ga else d_occ, window,
      n_draws_in_filter = 1L, fixed_draw_idx = r,
      draws_dur = if (grouped) NULL else d_dur,
      draws_mon = if (grouped) NULL else d_mon,
      cb_meta_mon = if (grouped) NULL else cb_meta_mon,
      grouping = if (grouped) grouping else NULL,
      measurement = "occurrence_duration_mon_NB", ...)
    res <- run_one_bpfilter_replicate(obj, Np = 40L,
                                      block_list = make_block_list(obj, frame_tu),
                                      seed = 900L + r)
    list(tu = "hmmW", rep_id = r, draw_d = r, filt_months = attr(obj, "filt_months"),
         filter_traj = res$filter_traj, logLik = res$logLik,
         gids = sort(as.character(frame_tu$gid), method = "radix"))
  })
  set.seed(31337L)
  sm <- extract_smoothing_one_tu(results, frame_tu, panel_tu,
                                 FILTER_START = as.Date("2024-01-01"),
                                 FILTER_END   = as.Date("2025-01-01"),
                                 grouping     = if (grouped) grouping else NULL,
                                 draws_mon    = if (grouped) NULL else d_mon,
                                 cb_meta_mon  = if (grouped) NULL else cb_meta_mon)
  list(cases = lapply(sm$reps, `[[`, "cases_mon_mat"),
       yout  = lapply(sm$reps, `[[`, "Yout_mat"))
}

out <- list(
  v1_ungrouped = run_v1(FALSE),
  v1_grouped   = run_v1(TRUE),
  sm_ungrouped = run_sm(FALSE),
  sm_grouped   = run_sm(TRUE)
)

# Fixture sanity: every site must actually produce live monthly means/draws,
# or the golden would pin nothing.
expect(sum(!is.na(out$v1_ungrouped$cases_month)) > 20, "fixture: v1 ungrouped draws cases")
expect(sum(!is.na(out$v1_grouped$cases_month))   > 20, "fixture: v1 grouped draws cases")
expect(sum(!is.na(unlist(out$sm_ungrouped$cases))) > 10, "fixture: smoothing ungrouped draws cases")
expect(sum(!is.na(unlist(out$sm_grouped$cases)))   > 10, "fixture: smoothing grouped draws cases")

if (cli$capture) {
  dir.create(dirname(GOLDEN_PATH), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, GOLDEN_PATH)
  cat("captured", length(out), "sites ->", GOLDEN_PATH, "\n")
} else {
  if (!file.exists(GOLDEN_PATH)) stop("golden missing: run with --capture first")
  gold <- readRDS(GOLDEN_PATH)
  expect(identical(sort(names(gold)), sort(names(out))), "same sites as the golden")
  for (nm in names(out))
    expect(identical(out[[nm]], gold[[nm]]), sprintf("%s identical to golden", nm))
}

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
