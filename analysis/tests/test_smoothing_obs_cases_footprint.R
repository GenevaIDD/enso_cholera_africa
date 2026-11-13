# test_smoothing_obs_cases_footprint.R -- 08c's observed monthly-cases series
# (utils/smoothing_core.R::smoothing_obs_cases_mon()) uses the filter's footprint
# rule per (gid, month).
#
# HARD gate (see analysis/tests/run_all.R). 08c summed every monthly-frame row
# per month x group, while the filter scores one value per (gid, month): the sum
# within an admin level, the maximum across levels
# (spatpomp_core.R::.aggregate_footprint_counts()). Under lp_season a gid-month
# reported at two levels was counted twice on the observed side of the cases
# ribbons (+8.5% on the 1990_2019 lp_season frame).
#   [1] one admin level per gid (gid_season): identical to the pre-fix 08c code,
#       kept verbatim below as the reference;
#   [2] two levels at one gid-month: the level maximum, not the sum;
#   [3] per (gid, month) equal to build_obs_panel()'s cases_month (the filter);
#   [4] 08c calls smoothing_obs_cases_mon();
#   [5] footprint.R::finestLevelRows() (shared by 08c and 09's derive_smoothing()):
#       one level per key returns the input unchanged; its per-key sums equal
#       finestLevelTotals(); 09 applies it to the observed monthly cases.
#
# Synthetic data only -- no filter run, no data files.

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
# smoothing_core.R / spatpomp_core.R source() their kernel closures; evaluate
# only the definitions needed here.
.defs <- function(path, want) {
  for (e in parse(path, keep.source = FALSE))
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% want) eval(e, globalenv())
}
.defs("analysis/utils/smoothing_core.R", "smoothing_obs_cases_mon")
.defs("analysis/utils/spatpomp_core.R", ".aggregate_footprint_counts")
.defs("analysis/utils/sim_io.R", "build_obs_panel")
.defs("analysis/utils/footprint.R", c("finestLevelRows", "finestLevelTotals", ".adminRank",
                                      ".ADMIN_LEVEL_RANK"))

cat("# 08c observed monthly cases: footprint rule\n")
.near <- function(a, b, label)
  expect(isTRUE(all.equal(a, b, tolerance = 1e-9)),
         paste0(label, if (!isTRUE(all.equal(a, b, tolerance = 1e-9)))
                  sprintf("  [got=%s want=%s]", deparse(a), deparse(b)) else ""))

# Pre-fix 08c section 5 (obs_cases_mon_*), verbatim apart from the arguments.
old_obs_cases_mon <- function(mon_dat, gid_meta, start, end, months, mask) {
  raw <- mon_dat %>%
    dplyr::mutate(gid = as.character(gid),
                  t_start = as.Date(month),
                  cases   = as.numeric(cases_month)) %>%
    dplyr::filter(!is.na(t_start), !is.na(cases),
                  t_start >= start, t_start < end) %>%
    mask("obs_cases_mon") %>%
    dplyr::select(-dplyr::any_of(c("AFRO_region", "country"))) %>%
    dplyr::filter(t_start %in% months) %>%
    dplyr::left_join(gid_meta, by = "gid") %>%
    dplyr::filter(!is.na(AFRO_region))
  list(
    africa = raw %>%
      dplyr::group_by(t_start) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(africa = "Africa"),
    region = raw %>%
      dplyr::group_by(t_start, AFRO_region) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop"),
    country = raw %>%
      dplyr::group_by(t_start, country) %>%
      dplyr::summarise(n_cases = sum(cases, na.rm = TRUE), .groups = "drop"))
}

set.seed(20261030)
months   <- seq(as.Date("2020-01-01"), by = "month", length.out = 10)
start    <- months[2]; end <- months[10]          # window [start, end)
gids     <- sprintf("AAA.%d_1", 1:8)
gid_meta <- tibble(gid = gids[1:7],               # gids[8] is outside the simulated frame
                   AFRO_region = rep(c("WEST", "EAST"), length.out = 7),
                   country     = rep(c("AAA", "BBB", "CCC"), length.out = 7))
lvl_of   <- setNames(rep(c("ADM1", "ADM2"), length.out = 8), gids)

# One level per gid (gid_season); duplicate rows in a gid-month (overlapping
# episodes), NA and fractional cases, months outside the window, a data
# region/country that disagrees with gid_meta.
n <- 140
mon_gs <- tibble(gid = sample(gids, n, TRUE),
                 month = sample(c(months, months[1] - 31), n, TRUE),
                 cases_month = ifelse(runif(n) < 0.1, NA, round(rexp(n, 1 / 40), 2)),
                 AFRO_region = "XX", country = "ZZZ") %>%
  mutate(admin_level = lvl_of[gid])
mask_ok <- tidyr::expand_grid(gid = gids, month = months) %>%      # 08c's coverage cells
  dplyr::slice_sample(prop = 0.6)
mask_fn <- function(df, label)
  dplyr::semi_join(dplyr::mutate(df, gid = as.character(gid)), mask_ok,
                   by = c("gid", "t_start" = "month"))
ident <- function(df, label) df

cat("\n[1] one admin level per gid: identical to the pre-fix 08c\n")
for (mk in list(list("no mask", ident), list("coverage mask", mask_fn))) {
  new <- smoothing_obs_cases_mon(mon_gs, gid_meta, start, end, months, mask = mk[[2]])
  old <- old_obs_cases_mon(mon_gs, gid_meta, start, end, months, mask = mk[[2]])
  for (nm in c("africa", "region", "country")) {
    ok <- identical(new[[nm]], old[[nm]])
    expect(ok, sprintf("%s, %s identical", mk[[1]], nm))
    if (!ok) print(all.equal(old[[nm]], new[[nm]]))
  }
}
expect(nrow(smoothing_obs_cases_mon(mon_gs[0, ], gid_meta, start, end, months)$country) == 0L,
       "empty frame -> empty tables")

cat("\n[2] nested reports at two levels: level maximum, not the sum\n")
m3 <- months[3]
mon_ml <- tibble(gid = gids[1], month = m3,
                 cases_month = c(100, 30, 40, NA),
                 admin_level = c("ADM1", "ADM2", "ADM2", "ADM2"))
mon_ml <- bind_rows(mon_ml, tibble(gid = gids[2], month = m3, cases_month = 25, admin_level = "ADM2"))
new <- smoothing_obs_cases_mon(mon_ml, gid_meta, start, end, months)
old <- old_obs_cases_mon(mon_ml, gid_meta, start, end, months, ident)
.near(new$africa$n_cases, 125, "Africa = max(100, 30 + 40) + 25")
.near(old$africa$n_cases, 195, "pre-fix 08c summed both levels (195)")
.near(new$country$n_cases[new$country$country == "AAA"], 100, "country AAA = 100")
mon_ml2 <- mutate(mon_ml, cases_month = replace(cases_month, 1, 50))
.near(smoothing_obs_cases_mon(mon_ml2, gid_meta, start, end, months)$africa$n_cases, 95,
             "a finer level with the larger total wins: 70 + 25")

cat("\n[3] per (gid, month) equal to the filter's build_obs_panel() cases_month\n")
mon_mix <- bind_rows(mon_gs, mon_ml)
surv <- tidyr::expand_grid(gid = gids, t_start = months) %>% mutate(outbreak = 0L)
obs  <- build_obs_panel(surv, mon_mix) %>% filter(!is.na(cases_month))
per_gid <- function(g, m)
  smoothing_obs_cases_mon(filter(mon_mix, gid == g, month == m),
                          tibble(gid = g, AFRO_region = "R", country = "C"),
                          months[1], months[10] + 1, months)$africa$n_cases
chk <- obs %>% filter(t_start %in% months) %>% slice_head(n = 40)
vals <- mapply(per_gid, chk$gid, chk$t_start)
vals <- vapply(vals, function(v) if (length(v)) v else NA_real_, numeric(1))
.near(unname(vals), chk$cases_month, "40 gid-months match the filter")
.near(per_gid(gids[1], m3), obs$cases_month[obs$gid == gids[1] & obs$t_start == m3],
             "the two-level gid-month matches the filter")

cat("\n[4] 08c uses the shared helper\n")
src <- readLines("analysis/08c_smoothing_compile.R")
expect(any(grepl("smoothing_obs_cases_mon(", src, fixed = TRUE)), "08c calls smoothing_obs_cases_mon()")
expect(!any(grepl("obs_cases_mon_raw", src, fixed = TRUE)), "no row-level obs_cases_mon_raw sum left in 08c")

cat("\n[5] finestLevelRows(): shared by 08c and 09\n")
gs <- mon_gs %>% mutate(t_start = as.Date(month)) %>% filter(!is.na(cases_month))
expect(identical(finestLevelRows(gs, "cases_month", c("gid", "t_start")), gs),
       "one level per key: input returned unchanged (rows, order, columns)")
expect(identical(finestLevelRows(select(gs, -admin_level), "cases_month", c("gid", "t_start")),
                 select(gs, -admin_level)), "no admin_level column: unchanged")
mix <- bind_rows(gs, mon_ml %>% mutate(t_start = as.Date(month)) %>% filter(!is.na(cases_month)))
kept <- finestLevelRows(mix, "cases_month", c("gid", "t_start"))
tot  <- finestLevelTotals(mix, "cases_month", c("gid", "t_start"))
chk5 <- kept %>% group_by(gid, t_start) %>% summarise(s = sum(cases_month), nl = n_distinct(admin_level),
                                                      .groups = "drop") %>%
  inner_join(tot, by = c("gid", "t_start"))
expect(nrow(chk5) == nrow(tot) && all(chk5$nl == 1L), "one level kept per key")
.near(chk5$s, chk5$total, "kept rows sum to finestLevelTotals()")
rd <- readLines("analysis/utils/results_derive.R")
i0 <- grep("^derive_smoothing <- function", rd); i1 <- grep("^ensure_smoothing <- function", rd)
expect(any(grepl('finestLevelRows("obs_cases"', rd[i0:i1], fixed = TRUE)),
       "09 derive_smoothing() keeps one level of the observed monthly cases")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
