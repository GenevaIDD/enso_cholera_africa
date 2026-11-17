# test_hindcast_skill_primitives.R — HARD: utils/hindcast_skill.R scoring
#   primitives.
#
# Usage:
#   Rscript analysis/tests/test_hindcast_skill_primitives.R
#
# Details:
#   Synthetic temp files only -- no fit, no generated_data. Pins the pieces
#   the fair hindcast scoring rests on: the CRPS estimators, the lead
#   classes, the leave-one-window-out reference (must not see the scored
#   window), the case supports' observed attribution, the episode-key refusal
#   of the run index, the PIT-vs-panel support guard of the per-run loader,
#   and the horizon cohorts of per-country season ends.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(qs2)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c(
  "base",
  "config",
  "paths",
  "provenance",
  "gid",
  "occ_data",
  "sim_io",
  "sim_metrics",
  "obs_sim_compare",
  "population",
  "hindcast_skill"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
near <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))

cat("# scoring rules\n")
set.seed(1)
x <- rpois(25, 2)
y <- 3
brute_fair <- mean(abs(x - y)) - sum(abs(outer(x, x, "-"))) / (2 * 25 * 24)
expect(
  near(crps_fair(x, y), brute_fair),
  "crps_fair equals the brute-force unbiased form"
)
expect(
  near(crps_ens(x, y), mean(abs(x - y)) - .5 * mean(abs(outer(x, x, "-")))),
  "crps_ens equals the notebook's crps_emp"
)
expect(crps_fair(5, 2) == 3, "crps_fair of a 1-member ensemble is |x - y|")
expect(
  near(crps_discrete(c(0, 0, 1), 0), 2) &&
    near(crps_discrete(c(0, 0, 1), 2), 0),
  "crps_discrete of a point mass at k is |y - k|"
)
big <- rbinom(2e5, 6, .3)
expect(
  abs(crps_binom(6, .3, 2) - crps_fair(big, 2)) < .01,
  "crps_binom (exact) agrees with the sampled fair CRPS"
)
expect(
  near(pinball_loss(c(1, 5), 3, c(.9, .1)), c(.9 * 2, .9 * 2)),
  "pinball_loss is tau*(y-q) above and (1-tau)*(q-y) below"
)
expect(
  near(skill_ratio(c(1, 1, 5), c(2, 2, NA)), .5),
  "skill_ratio ignores rows with no reference score"
)
expect(
  near(interval_score(1, 3, 4, .8), 2 + 10 * 1),
  "interval_score adds the 2/alpha penalty"
)

cat("# lead classes\n")
expect(
  identical(
    hindcast_lead_class(c(0, 7, 8, 12, 13, 20)),
    c("0-7", "0-7", "8-12", "8-12", "13+", "13+")
  ),
  "lead-class boundaries"
)
expect(
  identical(hindcast_lead_levels(), c("all", "0-7", "8-12", "13+")),
  "lead-class levels"
)

cat("# occurrence references\n")
mk <- function(m) as.Date(sprintf("2015-%02d-01", m))
pool <- tibble(
  group = "g",
  gid = "A",
  month = mk(1:12),
  y = c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0)
)
windows <- tibble(
  group = "g",
  launch_year = 2015L,
  launch = mk(5),
  horizon_end = mk(7)
)
units <- tibble(group = "g", launch_year = 2015L, gid = "A")
r <- occ_reference_rate("loo_gid", pool, obs_panel = NULL, windows, units)
expect(
  near(r$p_ref, 0) && r$m_ref == 9,
  "loo_gid leaves out every month of the scored window"
)
r2 <- occ_reference_rate(
  "loo_gid_fullhist",
  NULL,
  obs_panel = transmute(pool, gid, month, y = c(1, rep(0, 11))),
  windows,
  units
)
expect(
  near(r2$p_ref, 1 / 9),
  "loo_gid_fullhist uses the panel outside the window"
)
u <- tibble(
  group = "g",
  launch_year = 2015:2017,
  grain = "window",
  gid = "A",
  unit = 2015:2017,
  lead_class = "all",
  n_atrisk = 3L,
  y_obs = c(0L, 1L, 3L),
  months = list(mk(5:7), mk(5:7), mk(5:7))
)
sr <- score_occurrence_references(
  u,
  c("insample_gid", "zero"),
  "crps_emp",
  pool,
  NULL,
  windows
)
ins <- filter(sr, reference == "insample_gid")
expect(
  near(
    ins$score_ref,
    vapply(u$y_obs, function(v) crps_ens(u$y_obs, v), numeric(1))
  ),
  "insample_gid scores the gid's own values (legacy fig-9 baseline)"
)
expect(
  near(filter(sr, reference == "zero")$score_ref, c(0, 1, 3)),
  "zero reference CRPS = y"
)

cat("# case supports\n")
ms <- expand.grid(
  replicate_id = 1:2,
  gid = c("AAA_1_1", "AAA_1_2", "AAA_2_1"),
  month = c("2015-05-01", "2015-06-01"),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  mutate(Y = 0L, cases_month = 1)
mon_all <- build_mon_obs_frame(tibble(
  gid = "AAA_1",
  admin_level = "ADM1",
  month = as.Date("2015-05-01"),
  cases_month = 10
))
gid_pop <- tibble(
  gid = c("AAA_1", "AAA_1_1", "AAA_1_2", "AAA_2", "AAA_2_1"),
  Q = c(100, 30, 20, 50, 50),
  gadm_level = c("ADM1", "ADM2", "ADM2", "ADM1", "ADM2"),
  year = 2015L
)
run <- list(
  bases = obs_sim_bases(ms, mon_all, gid_pop),
  cells = tibble(gid = "AAA_1_1", month = as.Date("2015-05-01"))
)
obs_of <- function(sp) {
  sum(case_support(run, sp, mon_all, gid_pop)$obs_month$cases_month)
}
expect(
  near(obs_of("obs_positive_keys"), run$bases$obs_B) &&
    near(obs_of("obs_positive_keys"), 5),
  "obs_positive_keys reproduces basis B (10 x 50/100)"
)
expect(
  near(obs_of("atrisk_gids"), 3),
  "atrisk_gids scales by the at-risk gids' population share"
)
expect(
  near(obs_of("footprint"), 5) && near(obs_of("atrisk_keys"), 5),
  "footprint / atrisk_keys keep the whole key's population share"
)
expect(
  identical(
    sort(case_support(run, "footprint", mon_all, gid_pop)$gids),
    c("AAA_1_1", "AAA_1_2", "AAA_2_1")
  ),
  "footprint keeps every sim gid"
)
expect(
  !is.null(run$bases$obs_rk) && attr(mon_all, "n_na_lvl") == 0,
  "obs_sim_bases exposes obs_rk; build_mon_obs_frame counts NA levels"
)

cat("# run index + loader guards\n")
tmp <- tempfile("hcsk_")
dir.create(tmp)
cfg <- file.path(tmp, "r.yml")
yaml::write_yaml(
  list(
    outdir = file.path(tmp, "gd/"),
    simulator = list(run_tag = "r", output = list(format = "qs"))
  ),
  cfg
)
md <- file.path(tmp, "merged")
dir.create(md)
qs_save(ms, file.path(md, "monthly_summary.qs"))
man <- list(
  groups = list(
    g = list(
      years = list(
        `2015` = list(
          config = cfg,
          merged_dir = md,
          launch = "2015-05-01",
          horizon_end = "2015-06-01"
        )
      )
    )
  )
)
err <- tryCatch(
  {
    suppressMessages(hindcast_run_index(man, "gid_season"))
    ""
  },
  error = function(e) conditionMessage(e)
)
expect(
  grepl("episodes.key", err),
  "hindcast_run_index refuses a run under another episodes.key"
)
idx <- suppressMessages(tryCatch(
  hindcast_run_index(man, "lp_season"),
  error = function(e) e
))
expect(
  inherits(idx, "error") && grepl("metrics", conditionMessage(idx)),
  "hindcast_run_index stops when validation/metrics is missing"
)
vd <- file.path(tmp, "valid")
dir.create(vd)
pit <- tibble(
  gid = "AAA_1_1",
  month = as.Date(c("2015-05-01", "2015-06-01")),
  p_hat = .5,
  y_obs = 0L,
  pit = .3
)
qs_save(list(pit = pit), file.path(vd, "metrics.qs"))
row <- tibble(
  group = "g",
  launch_year = 2015L,
  launch = as.Date("2015-05-01"),
  horizon_end = as.Date("2015-06-01"),
  run_tag = "r",
  episode_key = "lp_season",
  monthly_summary_path = file.path(md, "monthly_summary.qs"),
  metrics_path = file.path(vd, "metrics.qs"),
  out_fmt = "qs"
)
panel_ok <- tibble(gid = "AAA_1_1", month = pit$month, y = 0L)
panel_bad <- panel_ok[1, ]
lr <- suppressMessages(hindcast_load_run(row, panel_ok, mon_all, gid_pop))
expect(
  nrow(lr$cells) == 2 && identical(lr$cells$lead, c(0L, 1L)),
  "hindcast_load_run accepts PIT cells equal to the panel's at-risk cells"
)
err <- tryCatch(
  {
    suppressMessages(hindcast_load_run(row, panel_bad, mon_all, gid_pop))
    ""
  },
  error = function(e) conditionMessage(e)
)
expect(
  grepl("at-risk", err),
  "hindcast_load_run stops when PIT cells differ from the panel"
)

# Each country is read up to its own firm horizon inside one (group, launch
# year) row; identical to a single-window read when the members share one
# horizon.
cat("# horizon cohorts (per-country season ends)\n")
# One (group, launch year) row whose countries end at different months: each is
# read up to its own firm horizon; without per-country horizons, the row's.
mo <- seq(as.Date("2015-05-01"), as.Date("2015-09-01"), by = "month")
ms2 <- tidyr::expand_grid(
  replicate_id = 1:2,
  gid = c("AAA_1_1", "BBB_1_1"),
  month = mo
) |>
  mutate(Y = 0L, cases_month = 0)
md2 <- file.path(tmp, "merged2")
dir.create(md2)
qs_save(ms2, file.path(md2, "monthly_summary.qs"))
opt2 <- make_options_from_config(
  option_list = list(),
  enforce_options = character(0),
  cli_argv = c("--config", cfg),
  model = NULL
)
vd2 <- make_run_file(opt2, kind = "spatpomp_valid_dir")
dir.create(vd2, recursive = TRUE)
pit2 <- tidyr::expand_grid(gid = c("AAA_1_1", "BBB_1_1"), month = mo) |>
  mutate(p_hat = .5, y_obs = 0L, pit = .3)
qs_save(list(pit = pit2), file.path(vd2, "metrics.qs"))
yrs <- function(h) {
  list(
    `2015` = list(
      config = cfg,
      merged_dir = md2,
      launch = "2015-05-01",
      horizon_end = h
    )
  )
}
man_c <- list(
  schema = 2L,
  groups = list(
    g1 = list(
      iso3 = list("AAA", "BBB"),
      units = list("u"),
      years = list(
        `2015` = list(launch = "2015-05-01", horizon_end = "2015-08-01")
      )
    )
  ),
  units = list(
    u = list(
      members = list("AAA", "BBB"),
      groups = list("g1"),
      years = yrs("2015-08-01")
    )
  ),
  countries = list(
    AAA = list(
      group = "g1",
      unit = "u",
      years = list(`2015` = list(firm_horizon_end = "2015-06-01"))
    ),
    BBB = list(
      group = "g1",
      unit = "u",
      years = list(`2015` = list(firm_horizon_end = "2015-08-01"))
    )
  )
)
ix <- suppressMessages(hindcast_run_index(man_c, "lp_season"))
expect(
  nrow(ix) == 1L &&
    ix$horizon_end == as.Date("2015-08-01") &&
    identical(ix$countries[[1]], c("AAA", "BBB")) &&
    identical(unname(ix$country_horizons[[1]]), c("2015-06-01", "2015-08-01")),
  "v2 index: one row per (group, year), horizon_end = the group's latest, per-country horizons"
)
rp <- .hindcast_read_part(ix[1, ])
expect(
  max(rp$ms$month[startsWith(rp$ms$gid, "AAA")]) == as.Date("2015-06-01") &&
    max(rp$ms$month[startsWith(rp$ms$gid, "BBB")]) == as.Date("2015-08-01") &&
    max(rp$pit$month[startsWith(rp$pit$gid, "AAA")]) == as.Date("2015-06-01"),
  "each country read up to its own firm horizon (sim and PIT)"
)
ix0 <- ix
ix0$country_horizons <- list(NULL)
rp0 <- .hindcast_read_part(ix0[1, ])
expect(
  max(rp0$ms$month[startsWith(rp0$ms$gid, "AAA")]) == as.Date("2015-08-01"),
  "without per-country horizons: the row's window (former behaviour)"
)
man_s <- man_c
man_s$countries$AAA$years$`2015`$firm_horizon_end <- "2015-08-01"
ix_s <- suppressMessages(hindcast_run_index(man_s, "lp_season"))
rp_s <- .hindcast_read_part(ix_s[1, ])
rp_s0 <- .hindcast_read_part(ix0[1, ])
expect(
  identical(rp_s$ms, rp_s0$ms) && identical(rp_s$pit, rp_s0$pit),
  "one shared horizon: identical to the single-window read"
)
man_v1 <- list(
  groups = list(
    g1 = list(iso3 = list("AAA", "BBB"), years = yrs("2015-08-01"))
  ),
  countries = man_c$countries
)
ix1 <- suppressMessages(hindcast_run_index(man_v1, "lp_season"))
expect(
  nrow(ix1) == 1L &&
    identical(unname(ix1$country_horizons[[1]]), c("2015-06-01", "2015-08-01")),
  "v1 layout (unit = group) with per-country horizons: members from the group's iso3"
)

cat("# season windows (make_hindcast_configs.R --season_window)\n")
# Unit u (AAA, BBB; fit group g1), season 2015: run "2015-05" scores AAA
# (May-Jun), run "2015-07" scores BBB (Jul-Aug). Each run simulates both
# countries from its launch; 07c scored only its own country (the PIT).
mk_run <- function(tag, months, pit_gid, pit_months) {
  cf <- file.path(tmp, paste0(tag, ".yml"))
  yaml::write_yaml(
    list(
      outdir = file.path(tmp, "gd/"),
      simulator = list(run_tag = tag, output = list(format = "qs"))
    ),
    cf
  )
  mdr <- file.path(tmp, paste0("merged_", tag))
  dir.create(mdr)
  qs_save(
    tidyr::expand_grid(
      replicate_id = 1:2,
      gid = c("AAA_1_1", "BBB_1_1"),
      month = months
    ) |>
      mutate(Y = 0L, cases_month = 1),
    file.path(mdr, "monthly_summary.qs")
  )
  o <- make_options_from_config(
    option_list = list(),
    enforce_options = character(0),
    cli_argv = c("--config", cf),
    model = NULL
  )
  v <- make_run_file(o, kind = "spatpomp_valid_dir")
  dir.create(v, recursive = TRUE)
  qs_save(
    list(
      pit = tibble(
        gid = pit_gid,
        month = pit_months,
        p_hat = .5,
        y_obs = 0L,
        pit = .3
      )
    ),
    file.path(v, "metrics.qs")
  )
  list(config = cf, merged_dir = mdr)
}
ra <- mk_run(
  "ra",
  seq(as.Date("2015-05-01"), as.Date("2015-09-01"), by = "month"),
  "AAA_1_1",
  as.Date(c("2015-05-01", "2015-06-01"))
)
rb <- mk_run(
  "rb",
  seq(as.Date("2015-07-01"), as.Date("2015-09-01"), by = "month"),
  "BBB_1_1",
  as.Date(c("2015-07-01", "2015-08-01"))
)
man_w <- list(
  schema = 2L,
  season_end = "season_window",
  groups = list(
    g1 = list(
      iso3 = list("AAA", "BBB"),
      units = list("u"),
      years = list(
        `2015` = list(launch = "2015-05-01", horizon_end = "2015-08-01")
      )
    )
  ),
  units = list(
    u = list(
      members = list("AAA", "BBB"),
      groups = list("g1"),
      years = list(
        `2015-05` = c(
          ra,
          list(
            launch = "2015-05-01",
            horizon_end = "2015-06-01",
            launch_year = 2015L,
            scored = list("AAA")
          )
        ),
        `2015-07` = c(
          rb,
          list(
            launch = "2015-07-01",
            horizon_end = "2015-08-01",
            launch_year = 2015L,
            scored = list("BBB")
          )
        )
      )
    )
  ),
  countries = list(
    AAA = list(
      group = "g1",
      unit = "u",
      floored = TRUE,
      years = list(
        `2015-05` = list(
          season_end = "2015-06-01",
          firm_horizon_end = "2015-06-01"
        )
      )
    ),
    BBB = list(
      group = "g1",
      unit = "u",
      floored = FALSE,
      years = list(
        `2015-07` = list(
          season_end = "2015-07-01",
          firm_horizon_end = "2015-08-01"
        )
      )
    )
  )
)
ixw <- suppressMessages(hindcast_run_index(man_w, "lp_season"))
expect(
  nrow(ixw) == 2L &&
    all(ixw$group == "g1") &&
    all(ixw$launch_year == 2015L) &&
    all(ixw$horizon_end == as.Date("2015-08-01")) &&
    identical(ixw$countries, list("AAA", "BBB")) &&
    identical(
      ixw$country_launch,
      list(c(AAA = "2015-05-01"), c(BBB = "2015-07-01"))
    ),
  "index: one row per run, each with only its scored country, its launch and the group horizon"
)
panel_w <- tibble(
  gid = c("AAA_1_1", "AAA_1_1", "BBB_1_1", "BBB_1_1"),
  month = as.Date(c("2015-05-01", "2015-06-01", "2015-07-01", "2015-08-01")),
  y = 0L
)
mon_w <- build_mon_obs_frame(tibble(
  gid = c("AAA_1", "AAA_1", "BBB_1"),
  admin_level = "ADM1",
  month = as.Date(c("2015-05-01", "2015-08-01", "2015-07-01")),
  cases_month = c(10, 7, 4)
))
lw <- suppressMessages(hindcast_load_run(ixw, panel_w, mon_w, gid_pop))
lc <- arrange(lw$cells, gid, month)
expect(
  identical(lc$lead, c(0L, 1L, 0L, 1L)) &&
    lw$meta$launch == as.Date("2015-05-01") &&
    identical(
      lw$country_launch,
      c(AAA = as.Date("2015-05-01"), BBB = as.Date("2015-07-01"))
    ),
  "pooled season: each country's leads count from its own launch; meta launch = the earliest"
)
expect(
  identical(
    lw$country_horizon,
    c(AAA = as.Date("2015-06-01"), BBB = as.Date("2015-08-01"))
  ),
  "run$country_horizon: each country's firm horizon (season mode)"
)
expect(
  identical(
    sort(paste(lw$mon_obs$gid, lw$mon_obs$month)),
    c("AAA_1 2015-05-01", "BBB_1 2015-07-01")
  ) &&
    max(lw$sim$month[startsWith(lw$sim$gid, "AAA")]) == as.Date("2015-06-01") &&
    min(lw$sim$month[startsWith(lw$sim$gid, "BBB")]) == as.Date("2015-07-01"),
  "observed cases and sim cut to each country's own window (AAA's August case dropped)"
)
expect(
  identical(
    hindcast_lead_months(
      as.Date(c("2015-07-01", "2016-01-01")),
      as.Date("2015-05-01")
    ),
    c(2L, 8L)
  ) &&
    identical(
      hindcast_lead_months(
        as.Date(c("2015-07-01", "2015-08-01")),
        as.Date("2015-05-01"),
        c(AAA = as.Date("2015-05-01"), BBB = as.Date("2015-07-01")),
        c("AAA", "BBB")
      ),
      c(2L, 1L)
    ),
  "hindcast_lead_months(): from the run launch, or from each row's country launch"
)
expect(
  is.null(lr$country_launch) && is.null(lr$mon_obs),
  "outside season windows the loaded run carries no country launches (former output)"
)
cw <- arrange(lw$country_windows, iso3)
expect(
  identical(cw$iso3, c("AAA", "BBB")) &&
    identical(cw$launch, as.Date(c("2015-05-01", "2015-07-01"))) &&
    identical(cw$season_end, as.Date(c("2015-06-01", "2015-07-01"))) &&
    identical(cw$firm_horizon_end, as.Date(c("2015-06-01", "2015-08-01"))) &&
    identical(cw$floored, c(TRUE, FALSE)),
  "run$country_windows: each country's launch, season end, firm horizon and floor flag"
)
expect(
  identical(lr$country_windows$launch, as.Date("2015-05-01")) &&
    identical(lr$country_windows$firm_horizon_end, as.Date("2015-06-01")) &&
    is.na(lr$country_windows$season_end),
  "run$country_windows outside season windows: the run's launch and horizon"
)
bo <- lw$bases$obs_month
bs <- lw$bases$sim_in
expect(
  identical(bo$lead[bo$country == "BBB"], 0L) &&
    identical(bo$lead[bo$country == "AAA"], 0L) &&
    all(c("lead", "lead_class") %in% names(bs)) &&
    identical(range(bs$lead[bs$country == "BBB"]), c(0L, 1L)),
  "bases$obs_month / bases$sim_in carry lead and lead_class from each country's launch"
)

sc_w <- score_cases_run(lw, "footprint", mon_w, gid_pop, c(7L, 12L), 0.5)
expect(
  nrow(sc_w) > 0 &&
    all(sc_w$lead_class %in% c("all", "0-7")) &&
    near(
      sc_w$obs[
        sc_w$grain == "country" & sc_w$scope == "BBB" & sc_w$lead_class == "all"
      ],
      4
    ),
  "score_cases_run on a pooled season: leads from each country launch, obs in its window"
)

cat("# vectorised kernel\n")
set.seed(3)
ens <- lapply(c(1, 2, 7, 40, 0), function(n) rnbinom(n, mu = 20, size = .5))
yv <- c(3, 0, 25, 11, 4)
S <- sort_rows(pad_rows(ens))
nn <- rowSums(!is.na(S))
expect(
  identical(nn, c(1, 2, 7, 40, 0)) &&
    all(diff(S[4, ]) >= 0) &&
    all(is.na(S[5, ])),
  "pad_rows / sort_rows: NA-padded, each row sorted"
)
ref_fair <- vapply(
  1:5,
  function(i) if (length(ens[[i]])) crps_fair(ens[[i]], yv[i]) else NA_real_,
  0
)
ref_ens <- vapply(1:4, function(i) crps_ens(ens[[i]], yv[i]), 0)
expect(
  near(crps_sorted(S, yv, nn), ref_fair, 1e-12),
  "crps_sorted (fair) == crps_fair, ragged rows, empty row NA"
)
expect(
  near(crps_sorted(S[1:4, ], yv[1:4], nn[1:4], fair = FALSE), ref_ens, 1e-12),
  "crps_sorted (energy) == crps_ens"
)
pr <- c(.025, .25, .5, .9, .975)
expect(
  near(
    quantile7_sorted(S[1:4, ], pr, nn[1:4]),
    t(vapply(
      1:4,
      function(i) stats::quantile(ens[[i]], pr, names = FALSE),
      numeric(5)
    )),
    1e-12
  ),
  "quantile7_sorted == stats::quantile(type = 7)"
)
tw <- 10
twi <- function(x, y, t) {
  # brute force: integral over z > t of (F(z) - 1{y <= z})^2
  z <- seq(t, max(c(x, y, t)) + 1, by = 1e-3)
  sum((ecdf(x)(z) - (y <= z))^2) * 1e-3
}
sc <- score_sorted(
  S[3:4, ],
  yv[3:4],
  c("twcrps", "crps_energy"),
  nn[3:4],
  tw = tw
)
expect(
  near(
    crps_sorted(pmax(S[3:4, ], tw), pmax(yv[3:4], tw), nn[3:4], fair = FALSE),
    c(twi(ens[[3]], yv[3], tw), twi(ens[[4]], yv[4], tw)),
    1e-2
  ),
  "threshold-weighted CRPS == the integral above the threshold (chaining identity)"
)
expect(
  near(
    score_sorted(
      S[2:4, ],
      yv[2:4],
      "pinball",
      nn[2:4],
      taus = c(.1, .5, .9)
    )$pinball,
    vapply(
      2:4,
      function(i) {
        mean(pinball_loss(
          stats::quantile(ens[[i]], c(.1, .5, .9), names = FALSE),
          yv[i],
          c(.1, .5, .9)
        ))
      },
      0
    ),
    1e-12
  ),
  "score_sorted pinball == mean pinball loss of the type-7 quantiles"
)
expect(
  identical(
    rows_by_unit(c(2L, 1L, 2L), c(5, 7, 9), 3L)[, 1:2],
    matrix(c(7, 5, NA, NA, 9, NA), 3, 2)
  ),
  "rows_by_unit: one NA-padded row per unit"
)

cat("# bootstrap and area summaries\n")
df <- tibble(
  scope = "s",
  launch_year = rep(2012:2015, each = 3),
  area = rep(c("a", "b", "c"), 4),
  group = "g",
  score_model = c(1, 2, 3, 2, 2, 2, 5, 1, 1, 3, 3, 3),
  score_ref = c(2, 2, 2, 3, 1, 4, 2, 2, 2, NA, 3, 6)
)
boot <- hindcast_boot_spec(
  df$launch_year,
  list(bootstrap_n = 200L, bootstrap_seed = 7L)
)
set.seed(99)
r1 <- runif(1)
set.seed(99)
b0 <- hindcast_boot_spec(
  df$launch_year,
  list(bootstrap_n = 200L, bootstrap_seed = 7L)
)
expect(
  identical(runif(1), r1) && identical(boot, b0),
  "hindcast_boot_spec: own seed, global RNG state restored"
)
sk <- .skill_rollup(df, "scope", .skill_stats_units(), boot, "area")
brute <- vapply(
  seq_len(nrow(boot$counts)),
  function(b) {
    w <- boot$counts[b, match(df$launch_year, boot$years)]
    ok <- is.finite(df$score_model) & is.finite(df$score_ref)
    1 - sum(w[ok] * df$score_model[ok]) / sum(w[ok] * df$score_ref[ok])
  },
  0
)
expect(
  near(sk$crpss, skill_ratio(df$score_model, df$score_ref)) &&
    near(
      c(sk$crpss_lo, sk$crpss_hi),
      stats::quantile(brute, c(.025, .975), names = FALSE),
      1e-12
    ),
  "season-year bootstrap CI == brute-force resampling of whole years"
)
area_cr <- vapply(
  c("a", "b", "c"),
  function(a) with(df[df$area == a, ], skill_ratio(score_model, score_ref)),
  0
)
expect(
  near(sk$median_area_crpss, median(area_cr)) &&
    near(sk$share_area_skill, mean(area_cr > 0)),
  "per-area median CRPSS and positive share"
)
expect(
  is.null(hindcast_boot_spec(2012:2015, list())),
  "no bootstrap_n: no CIs (former output)"
)

cat("# monthly / season case references (fixed support)\n")
# One country AAA, one ADM1 key AAA_1 (sim gid AAA_1_1, half its population),
# season May-Jun 2015; observed ADM1 cases in every month 2010-2017 = 100 x
# (year - 2009) + month, so each reference value is traceable.
hm <- seq(as.Date("2010-01-01"), as.Date("2017-12-01"), by = "month")
yr <- as.integer(format(hm, "%Y"))
mo <- as.integer(format(hm, "%m"))
mon_h <- build_mon_obs_frame(tibble(
  gid = "AAA_1",
  admin_level = "ADM1",
  month = hm,
  cases_month = 100 * (yr - 2009) + mo
))
gp_h <- tibble(
  gid = c("AAA_1", "AAA_1_1", "AAA_1_2"),
  gadm_level = c("ADM1", "ADM2", "ADM2"),
  iso3 = "AAA",
  Q0 = c(100, 50, 50)
) |>
  tidyr::expand_grid(year = 2015:2030) |>
  mutate(Q = Q0 * 1.1^(year - 2015)) |>
  select(-Q0)
wpp_h <- tibble(
  iso3 = "AAA",
  year = 2000:2030,
  pop_total = 1.1^(2000:2030 - 2015)
)
win <- as.Date(c("2015-05-01", "2015-06-01"))
ms_h <- tidyr::expand_grid(replicate_id = 1:4, gid = "AAA_1_1", month = win) |>
  mutate(Y = 0L, cases_month = c(300, 310, 320, 330, 340, 350, 360, 370))
run_h <- list(
  meta = tibble(group = "g", launch_year = 2015L, launch = win[1]),
  sim = ms_h,
  bases = obs_sim_bases(ms_h, filter(mon_h, month %in% win), gp_h),
  cells = tibble(gid = "AAA_1_1", month = win),
  country_windows = tibble(
    iso3 = "AAA",
    launch = win[1],
    season_end = win[2],
    firm_horizon_end = win[2],
    floored = NA
  )
)
run_h$bases$sim_in <- mutate(run_h$bases$sim_in, cases_month = cases_month)
cfg_h <- list(
  case_skill_supports = "atrisk_gids",
  case_month_grains = "adm1",
  case_skill_references = c(
    "clim_pooled",
    "clim_launch_state",
    "persistence_season",
    "persistence_last"
  ),
  case_skill_rules = c("crps_raw", "crps_log1p", "twcrps"),
  case_clim_halfwidth = 1L,
  case_pop_scale = "none",
  case_mask = "none",
  case_tw_threshold_month = 400,
  case_tw_threshold_season = 400,
  case_season_min_ref = 3L
)
A_h <- list(
  atrisk_gids = case_unit_arrays(
    run_h,
    "atrisk_gids",
    filter(mon_h, month %in% win),
    gp_h,
    c(3L, 8L)
  )
)
H0 <- case_obs_history(mon_h, tibble(month = hm), gp_h, cfg_h)
u0 <- case_skill_units(run_h, A_h, H0, cfg_h, c(3L, 8L))
mu <- u0$month[u0$month$month == win[1], ]
v <- function(y, m) 0.5 * (100 * (y - 2009) + m) # observed on the support (share 1/2)
clim <- unlist(lapply(setdiff(2010:2017, 2015), function(y) v(y, 4:6)))
expect(
  nrow(u0$month) == 2 && near(mu$obs, v(2015, 5)) && mu$clim_pooled__n == 21,
  "month unit: observed on the fixed support; clim_pooled = Apr-Jun of the 7 other years"
)
expect(
  near(mu$clim_pooled__crps_raw, crps_fair(clim, mu$obs)) &&
    near(mu$clim_pooled__crps_log1p, crps_fair(log1p(clim), log1p(mu$obs))) &&
    near(mu$clim_pooled__twcrps, crps_fair(pmax(clim, 400), pmax(mu$obs, 400))),
  "clim_pooled scores: fair CRPS raw / log1p / threshold-weighted"
)
expect(
  near(mu$model__crps_raw, crps_fair(c(300, 320, 340, 360), mu$obs)),
  "model score: the support's replicate values"
)
expect(
  near(mu$persistence_season__crps_raw, abs(v(2014, 5) - mu$obs)) &&
    near(mu$persistence_last__crps_raw, abs(v(2015, 4) - mu$obs)),
  "persistence_season = May 2014; persistence_last = April 2015"
)
expect(
  mu$clim_launch_state__n == 21,
  "launch-state climatology: every year active in April (all kept)"
)
s_all <- u0$season[
  u0$season$grain == "country" & u0$season$lead_class == "all",
]
sclim <- vapply(setdiff(2010:2017, 2015), function(y) v(y, 5) + v(y, 6), 0)
expect(
  near(s_all$obs, v(2015, 5) + v(2015, 6)) &&
    s_all$clim_pooled__n == 7 &&
    near(s_all$clim_pooled__crps_raw, crps_fair(sclim, s_all$obs)),
  "season total: May+Jun of the 7 other years"
)
expect(
  s_all$persistence_season__n == 1L &&
    near(
      s_all$persistence_season__crps_raw,
      abs(v(2014, 5) + v(2014, 6) - s_all$obs)
    ),
  "season persistence: May+Jun 2014, not cut by the season minimum (point forecast)"
)
# a 14-month window (May 2015 - Jun 2016) overlaps the neighbouring seasons:
# May 2014 - Jun 2017 is left out, so Apr-Jun 2010-2013 and Apr 2014 remain
run_l <- run_h
run_l$country_windows$firm_horizon_end <- as.Date("2016-06-01")
u_l <- case_skill_units(run_l, A_h, H0, cfg_h, c(3L, 8L))
expect(
  u_l$month$clim_pooled__n[u_l$month$month == win[1]] == 13,
  "leave-season-out: the seasons overlapping the scored window are dropped"
)
# coverage: April-June 2012 unsurveilled (no case rows either) -> 2012 drops out
mon_c <- filter(
  mon_h,
  !(month %in% as.Date(c("2012-04-01", "2012-05-01", "2012-06-01")))
)
cov <- tibble(gid = "AAA_1_1", month = hm) |>
  filter(!month %in% as.Date(c("2012-04-01", "2012-05-01", "2012-06-01")))
Hc <- case_obs_history(
  mon_c,
  tibble(month = hm),
  gp_h,
  modifyList(cfg_h, list(case_mask = "coverage")),
  coverage = cov
)
uc <- case_skill_units(
  run_h,
  A_h,
  Hc,
  modifyList(cfg_h, list(case_mask = "coverage")),
  c(3L, 8L)
)
expect(
  uc$month$clim_pooled__n[uc$month$month == win[1]] == 18 &&
    uc$season$clim_pooled__n[
      uc$season$grain == "country" & uc$season$lead_class == "all"
    ] ==
      6,
  "surveillance mask: unsurveilled months are no reference values (months and season totals)"
)
expect(
  inherits(
    tryCatch(
      case_obs_history(
        mon_h,
        tibble(month = hm),
        gp_h,
        modifyList(cfg_h, list(case_pop_scale = "wpp"))
      ),
      error = identity
    ),
    "error"
  ),
  "case_pop_scale 'wpp' without the WPP table stops (no silent nearest-year fallback)"
)
cfg_w <- modifyList(cfg_h, list(case_pop_scale = "wpp"))
Hw <- case_obs_history(mon_h, tibble(month = hm), gp_h, cfg_w, wpp = wpp_h)
uw <- case_skill_units(run_h, A_h, Hw, cfg_w, c(3L, 8L))
climw <- unlist(lapply(setdiff(2010:2017, 2015), function(y) {
  v(y, 4:6) * 1.1^(2015 - y)
}))
expect(
  near(
    uw$month$clim_pooled__crps_raw[uw$month$month == win[1]],
    crps_fair(climw, mu$obs)
  ),
  "population scaling: values moved to the scored year's population (WPP growth before 2015)"
)

# season windows: units count leads from each country's own launch
Hs <- case_obs_history(mon_w, panel_w, gid_pop, cfg_h)
As <- list(
  footprint = case_unit_arrays(lw, "footprint", mon_w, gid_pop, c(3L, 8L))
)
us <- case_skill_units(
  lw,
  As,
  Hs,
  modifyList(cfg_h, list(case_month_grains = "country")),
  c(3L, 8L)
)
mb <- arrange(us$month, country, month)
expect(
  identical(mb$country, c("AAA", "AAA", "BBB", "BBB")) &&
    identical(mb$lead, c(0L, 1L, 0L, 1L)) &&
    identical(
      mb$month,
      as.Date(c("2015-05-01", "2015-06-01", "2015-07-01", "2015-08-01"))
    ),
  "pooled season: month units inside each country's window, leads from its launch"
)

cat("# sim_cells: ragged replicate frames (missing, duplicated, NA rows)\n")
# The bundle's golden worlds are complete gid x month x replicate grids; these
# are the sim_cells() branches they do not reach, against plain dplyr sums.
set.seed(11)
mr <- as.Date(c("2015-05-01", "2015-06-01", "2015-07-01"))
rg <- tidyr::expand_grid(
  replicate_id = 1:5,
  gid = c("AAA_1_1", "AAA_1_2"),
  month = mr
) |>
  mutate(Y = rbinom(n(), 1, .5), cases_month = rpois(n(), 20) + 0)
rg <- rg[-c(2, 9, 17), ] # missing (cell, replicate) rows
rg <- rg[!(rg$gid == "AAA_1_2" & rg$replicate_id == 4L), ] # a gid without replicate 4
rg <- bind_rows(rg, rg[c(1, 5, 5), ]) # duplicated rows
rg$cases_month[c(3, 8)] <- NA
rg <- rg[sample(nrow(rg)), ]
sc <- sim_cells(rg, "cases_month", na.rm = TRUE)
ref <- rg |>
  group_by(gid, month, replicate_id) |>
  summarise(v = sum(cases_month, na.rm = TRUE), .groups = "drop")
ci <- match(paste(ref$gid, as.integer(ref$month)), paste(sc$gid, sc$day))
Xr <- matrix(0, length(sc$gid), length(sc$reps))
Xr[cbind(ci, match(ref$replicate_id, sc$reps))] <- ref$v
Pr <- matrix(FALSE, length(sc$gid), length(sc$reps))
Pr[cbind(ci, match(ref$replicate_id, sc$reps))] <- TRUE
expect(
  identical(sc$X, Xr) && identical(sc$present, Pr) && identical(sc$reps, 1:5),
  "sim_cells: duplicate rows summed (NA dropped), 0 and present = FALSE where a row is missing"
)
expect(
  is.null(
    sim_cells(
      distinct(rg, gid, month, replicate_id, .keep_all = TRUE) |>
        tidyr::complete(gid, month, replicate_id, fill = list(Y = 0L)),
      "Y"
    )$present
  ),
  "sim_cells: present is NULL on a complete grid"
)
run_r <- list(
  meta = tibble(launch_year = 2015L),
  sim = rg,
  cells = tibble(
    gid = rep(c("AAA_1_1", "AAA_1_2"), each = 3),
    month = rep(mr, 2),
    lead_class = rep(c("0-3", "0-3", "4-8"), 2),
    y_obs = 0L,
    year = 2015L
  )
)
oc <- onset_window_counts(run_r, "window")
oref <- bind_rows(mutate(run_r$cells, lead_class = "all"), run_r$cells) |>
  inner_join(rg, by = c("gid", "month"), relationship = "many-to-many") |>
  group_by(gid, lead_class, replicate_id) |>
  summarise(sc = sum(Y), .groups = "drop")
er <- matrix(NA_real_, nrow(oc$obs), length(oc$reps))
er[cbind(
  match(paste(oref$gid, oref$lead_class), paste(oc$obs$gid, oc$obs$lead_class)),
  match(oref$replicate_id, oc$reps)
)] <- oref$sc
expect(
  identical(oc$ens, er),
  "onset_window_counts: per-unit replicate sums over present rows, NA where a unit has none"
)

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
