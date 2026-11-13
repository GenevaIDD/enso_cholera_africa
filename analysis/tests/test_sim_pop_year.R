# test_sim_pop_year.R -- simulator.location_frame.pop_year:
# the simulator / smoothing population year, an integer or a per-month mode.
#
# HARD gate (see analysis/tests/run_all.R). Synthetic, no MCMC, no data.
#
#   A  simPopYearSpec(): default = get_default_options(), validation, the
#      simulator.population.base_year guard; simFramePopArgs().
#   B  seasonYearOf() = occ_data.R's seas_of(); framePopAt() / framePopMatrix().
#   C  build_location_frame(): an integer pop_year adds no column; a mode adds
#      pop_y<YYYY> = resolveGidPopulation() per year (T2 years through WPP
#      included) and pop_shift_m.
#   D  constant-population equivalence: a mode frame whose per-year columns all
#      equal `pop` gives output identical() to the scalar frame in every kernel
#      (simulate_one_tu, the C filter build / precompute_mon_lp_static, the
#      smoothing post-hoc draw), under a binding incidence cap too.
#   E  a population doubling from one season year on moves the mon log mean by
#      exactly log 2 (elasticity-1 offset) on those months only.
#   F  smoothing threading: smoothing_make_config() reads location_frame.pop_year
#      / gid_pop_path / wpp_path (the "not threaded" bug).
#   G  smoothing cache identity: pop_year legacy / match / mismatch; the deploy
#      cross-arm check.
#
# Run: Rscript analysis/tests/test_sim_pop_year.R   (from the project root)

suppressPackageStartupMessages({
  library(tibble); library(dplyr); library(lubridate); library(tidyr)
})
source("analysis/utils/smoothing_core.R")   # sources config, simulator_core, spatpomp_core{,_c}, population
source("analysis/utils/occ_data.R")         # seas_of() lives inside nationalZeroSurvRows(); re-stated below
source("analysis/utils/deploy_core.R")      # deploy_assert_sim_pop_year()
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

cat("# test_sim_pop_year\n\n")

# ===========================================================================
# A -- simPopYearSpec / simFramePopArgs
# ===========================================================================
cat("## A: simPopYearSpec()\n")
dflt <- get_default_options()$simulator
expect(identical(dflt$location_frame$pop_year, .SIM_POP_YEAR_DEFAULT) &&
       identical(.SIM_POP_YEAR_DEFAULT, 2020L),
       "A1 get_default_options() pop_year = .SIM_POP_YEAR_DEFAULT = 2020L")
expect(identical(simPopYearSpec(NULL), 2020L) && identical(simPopYearSpec(dflt), 2020L),
       "A2 unset / default -> 2020L")
expect(identical(simPopYearSpec(list(location_frame = list(pop_year = 2027))), 2027L) &&
       identical(simPopYearSpec(list(location_frame = list(pop_year = "season_year"))), "season_year") &&
       identical(simPopYearSpec(list(location_frame = list(pop_year = "calendar_year"))), "calendar_year"),
       "A3 integer year (YAML double) -> integer; the two modes pass")
bad <- list("2020", "season", 2020.5, 1900, NA, c(2020, 2021), TRUE)
expect(all(vapply(bad, function(v) !is.null(err_of(simPopYearSpec(list(location_frame = list(pop_year = v))))),
                  logical(1))), "A4 refuses strings, fractions, out-of-range, NA, vectors, logicals")
e <- err_of(simPopYearSpec(list(population = list(base_year = 2025))))
expect(!is.null(e) && grepl("location_frame.pop_year", e, fixed = TRUE),
       "A5 population.base_year != 2020 refused, pointing to location_frame.pop_year")
expect(is.null(err_of(simPopYearSpec(list(population = list(base_year = 2020L))))),
       "A6 population.base_year = 2020 accepted")
fpa <- simFramePopArgs(list(location_frame = list(gid_pop_path = "x.qs", wpp_path = "w.csv",
                                                  pop_year = "season_year")))
expect(identical(fpa, list(gid_pop_path = "x.qs", wpp_path = "w.csv", pop_year = "season_year")),
       "A7 simFramePopArgs: configured paths used as given")
fpa0 <- simFramePopArgs(list(), gid_pop_path = tempfile(), wpp_path = tempfile())
expect(is.null(fpa0$gid_pop_path) && is.null(fpa0$wpp_path) && identical(fpa0$pop_year, 2020L),
       "A8 simFramePopArgs: absent default paths -> NULL (legacy route), default year")

# ===========================================================================
# B -- seasonYearOf / framePopAt
# ===========================================================================
cat("\n## B: seasonYearOf(), framePopAt()\n")
seas_of_ref <- function(month, ms) as.integer(lubridate::year(month %m+% months(-(ms - 1L))))
grid <- tidyr::expand_grid(month = as.Date("2023-01-01") %m+% months(0:23), ms = 1:12)
expect(identical(seasonYearOf(grid$month, grid$ms),
                 mapply(function(m, s) seas_of_ref(m, s), grid$month, grid$ms, USE.NAMES = FALSE)),
       "B1 seasonYearOf() = occ_data.R seas_of() on a 24-month x 12-season grid")
fr_b <- tibble::tibble(gid = c("A", "B"), pop = c(1, 2), pop_shift_m = c(0L, 6L),
                       pop_y2023 = c(10, 20), pop_y2024 = c(11, 21))
expect(identical(framePopAt(fr_b, as.Date("2024-03-01")), c(11, 20)),
       "B2 framePopAt: calendar row reads 2024, season row (shift 6) reads 2023")
expect(identical(framePopAt(fr_b[, 1:2], as.Date("2024-03-01")), c(1, 2)) &&
       is.null(framePopMatrix(fr_b[, 1:2], as.Date("2024-03-01"))),
       "B3 no pop_y columns: framePopAt = static pop, framePopMatrix = NULL")
e <- err_of(framePopAt(fr_b, as.Date("2026-01-01")))
expect(!is.null(e) && grepl("2026", e), "B4 a year outside the frame's columns is an error")
pm <- framePopMatrix(fr_b, as.Date(c("2024-03-01", "2024-08-01")))
expect(identical(pm, matrix(c(11, 20, 11, 21), 2)), "B5 framePopMatrix is U x M of framePopAt")

# ===========================================================================
# C -- build_location_frame(pop_year =)
# ===========================================================================
cat("\n## C: build_location_frame()\n")
gids <- c("AAA_1_1", "AAA_2_1", "BBB_1_1")
gp <- tidyr::expand_grid(gid = gids, year = 2015:2030) %>%
  dplyr::mutate(iso3 = substr(gid, 1, 3), gadm_level = "ADM2",
                Q = 1e4 * (match(gid, gids)) * 1.02^(year - 2015))
gp_path <- tempfile(fileext = ".qs"); qs2::qs_save(gp, gp_path)
wpp <- tidyr::expand_grid(iso3 = c("AAA", "BBB"), year = 2000:2040) %>%
  dplyr::mutate(pop_total = 1e6 * 1.03^(year - 2000))
wpp_path <- tempfile(fileext = ".csv"); utils::write.csv(wpp, wpp_path, row.names = FALSE)
panel_c <- tibble::tibble(gid = gids, country = substr(gids, 1, 3), AFRO_region = "R",
                          min_seas = c(1L, 7L, 3L), tu = "T")
f_def <- build_location_frame(panel_c, gid_pop_path = gp_path, wpp_path = wpp_path)
f_int <- build_location_frame(panel_c, gid_pop_path = gp_path, wpp_path = wpp_path, pop_year = 2020L,
                              pop_window = as.Date(c("2010-01-01", "2012-01-01")))
expect(identical(f_def, f_int) && !framePopHasYears(f_int) && !"pop_shift_m" %in% names(f_int),
       "C1 integer pop_year: no new column, pop_window ignored, = the default call")
win <- as.Date(c("2013-05-01", "2016-02-01"))
f_sy <- build_location_frame(panel_c, gid_pop_path = gp_path, wpp_path = wpp_path,
                             pop_year = "season_year", pop_window = win)
yrs <- 2012:2016
expect(identical(grep("^pop_y", names(f_sy), value = TRUE), sprintf("pop_y%04d", yrs)),
       "C2 season_year: pop_y<year> for year(start) - 1 .. year(end)")
ok <- all(vapply(yrs, function(y) identical(f_sy[[sprintf("pop_y%04d", y)]],
  resolveGidPopulation(f_sy$gid, y, gp, loadWppTable(wpp_path), on_missing = "error")$pop), logical(1)))
expect(ok, "C3 each pop_y<year> = resolveGidPopulation() at that year (T2 WPP years included)")
expect(identical(f_sy$pop, resolveGidPopulation(f_sy$gid, 2013L, gp, loadWppTable(wpp_path))$pop),
       "C4 static pop at year(pop_window[1])")
expect(identical(f_sy$pop_shift_m, f_sy$min_seas - 1L), "C5 season_year: pop_shift_m = min_seas - 1")
f_cy <- build_location_frame(panel_c, gid_pop_path = gp_path, wpp_path = wpp_path,
                             pop_year = "calendar_year", pop_window = win)
expect(all(f_cy$pop_shift_m == 0L), "C6 calendar_year: pop_shift_m = 0")
expect(!is.null(err_of(build_location_frame(panel_c, gid_pop_path = gp_path, pop_year = "season_year"))),
       "C7 a mode without pop_window is refused")
expect(!is.null(err_of(build_location_frame(dplyr::mutate(panel_c, pop = 1), gid_pop_path = gp_path,
                                            pop_year = "season_year", pop_window = win))),
       "C8 a mode on a panel with a fixed `pop` is refused")

# ===========================================================================
# Fixture for D / E: test_incidence_cap.R's TU (2 countries x 2 ADM1)
# ===========================================================================
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2) %>%
  dplyr::mutate(gid = sprintf("%s_%d_1", country, adm1n),
                adm1_unit = sprintf("%s_%d", country, adm1n),
                AFRO_region = "REG_W", tu = "hmmW", min_seas = c(1L, 4L, 7L, 10L),
                pop = c(2e4, 1.8e5, 6e4, 3e5)) %>%
  dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  tibble::as_tibble()
# Mode frame with a constant population: every pop_y column = pop.
frame_mode <- frame_tu
for (y in 2023:2024) frame_mode[[sprintf("pop_y%04d", y)]] <- frame_tu$pop
frame_mode$pop_shift_m <- frame_tu$min_seas - 1L
months <- as.Date("2024-01-01") %m+% months(0:5)
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
  dplyr::left_join(frame_tu %>% dplyr::select(gid, min_seas), by = "gid") %>%
  dplyr::mutate(country = sub("_[0-9]+_[0-9]+$", "", gid),
                t_start_seas = as.integer(lubridate::yday(t_start) - 1L),
                mean_precip = 0, dev_precip = 0, mean_temp_daymean = 0, dev_temp_daymean = 0,
                outbreak = 0L, time_to_last_outbreak_cat = "[6,Inf]",
                ongoing_outbreak_lag1 = 0L, ongoing_outbreak_country_lag1 = 0L,
                ongoing_outbreak_tu_lag1 = 0L) %>%
  tibble::as_tibble()
horizon <- list(start = "2024-01-01", end = "2024-06-30")
window  <- list(start = "2024-01-01", end = "2024-07-01")
d_occ <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = 5))
d_dur <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -30),
                                     active_terms = "(Intercept)")
d_dur$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
d_mon <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -5),
                                     active_terms = "(Intercept)", phi = 0.45)
d_mon$offset_pop <- TRUE; d_mon$offset_logoffset <- TRUE

# ===========================================================================
# D -- constant-population equivalence, mode vs scalar
# ===========================================================================
cat("\n## D: constant-population equivalence\n")
run_v1 <- function(fr, cap = NULL) suppressWarnings(
  simulate_one_tu(tu = "hmmW", frame_tu = fr, panel_tu = panel_tu,
                  draws_occ = d_occ, draws_dur = d_dur, draws_mon = d_mon,
                  replicate_ids = 1:4, seed = 77L, output_part = tempfile("py_"),
                  horizon = horizon, dry_run = TRUE, incidence_cap = cap)$monthly_summary)
expect(identical(run_v1(frame_tu), run_v1(frame_mode)), "D1 simulate_one_tu: mode == scalar (no cap)")
expect(identical(run_v1(frame_tu, 1e-4), run_v1(frame_mode, 1e-4)),
       "D2 simulate_one_tu: mode == scalar under a binding cap (per-month cap path)")

build_c <- function(fr, cap = NULL)
  build_spatpomp_tu_c(fr, panel_tu, d_occ, window, n_draws_in_filter = 1L, fixed_draw_idx = 1L,
                      draws_dur = d_dur, draws_mon = d_mon,
                      measurement = "occurrence_duration_mon_NB", incidence_cap = cap)
o_s <- build_c(frame_tu, 0.15); o_m <- build_c(frame_mode, 0.15)
expect(identical(o_s@covar@table, o_m@covar@table) &&
       identical(o_s@covar@times, o_m@covar@times),
       "D3 C filter build: covariate table (mon_lp_static, mon_mu_cap) identical, mode == scalar")
sim_df <- function(o) { set.seed(2024L); spatPomp::simulate(o, format = "data.frame") }
expect(identical(sim_df(o_s), sim_df(o_m)), "D4 C filter build: seeded simulate() identical")
o_s0 <- build_c(frame_tu); o_m0 <- build_c(frame_mode)
expect(identical(o_s0@covar@table, o_m0@covar@table), "D5 C filter build without a cap: identical")

results <- lapply(1:2, function(r) {
  obj <- build_c(frame_tu, 0.15)
  res <- run_one_bpfilter_replicate(obj, Np = 20L, block_list = make_block_list(obj, frame_tu),
                                    seed = 700L + r)
  list(tu = "hmmW", rep_id = r, draw_d = 1L, filt_months = attr(obj, "filt_months"),
       filter_traj = res$filter_traj, logLik = res$logLik,
       gids = sort(as.character(frame_tu$gid), method = "radix"))
})
run_sm <- function(fr, cap = NULL) {
  set.seed(31L)
  sm <- extract_smoothing_one_tu(results, fr, panel_tu,
                                 FILTER_START = as.Date("2024-01-01"),
                                 FILTER_END = as.Date("2024-07-01"),
                                 draws_mon = d_mon, incidence_cap = cap)
  lapply(sm$reps, `[[`, "cases_mon_mat")
}
sm_s <- run_sm(frame_tu, 1e-4); sm_m <- run_sm(frame_mode, 1e-4)
expect(sum(!is.na(unlist(sm_s))) > 5 && identical(sm_s, sm_m),
       "D6 extract_smoothing_one_tu: mode == scalar under a binding cap")
expect(identical(run_sm(frame_tu), run_sm(frame_mode)), "D7 extract_smoothing_one_tu: mode == scalar, no cap")

# ===========================================================================
# E -- a population step moves the log mean by log 2 exactly where it applies
# ===========================================================================
cat("\n## E: population step\n")
covar_tbl <- build_covar_for_tu(panel_tu, frame_tu, months, mon_harmonics = TRUE)
lp0 <- precompute_mon_lp_static(d_mon, frame_tu, covar_tbl)
ct  <- covar_tbl
mi  <- pmin(ct$time + 1L, length(months))
gi  <- match(ct$unit, frame_tu$gid)
step <- months[mi] >= as.Date("2024-04-01")
ct$.pop_t <- frame_tu$pop[gi] * ifelse(step, 2, 1)
lp1 <- precompute_mon_lp_static(d_mon, frame_tu, ct)
expect(isTRUE(all.equal(lp1 - lp0, ifelse(step, log(2), 0))),
       "E1 precompute_mon_lp_static: pop x2 from 2024-04 moves the log mean by log 2 there only")
ct$.pop_t <- frame_tu$pop[gi]
expect(identical(precompute_mon_lp_static(d_mon, frame_tu, ct), lp0),
       "E2 precompute_mon_lp_static: `.pop_t` = pop reproduces the static lp exactly")

# ===========================================================================
# F -- smoothing threading
# ===========================================================================
cat("\n## F: smoothing_make_config() threading\n")
cfg_path <- tempfile(fileext = ".yml")
yaml::write_yaml(list(simulator = list(
  backend = "spatpomp",
  location_frame = list(pop_year = "season_year", gid_pop_path = "gp_x.qs", wpp_path = "wpp_x.csv"),
  filter = list(measurement = "occurrence_bernoulli"))), cfg_path)
S <- smoothing_make_config(list(config = cfg_path))
expect(identical(S$POP_YEAR, "season_year") && identical(S$GID_POP_PATH, "gp_x.qs") &&
       identical(S$WPP_PATH, "wpp_x.csv"),
       "F1 location_frame.pop_year / gid_pop_path / wpp_path reach the smoothing state")
pa <- smoothing_frame_pop_args(S)
expect(identical(pa$pop_year, "season_year") && identical(pa$gid_pop_path, "gp_x.qs") &&
       identical(pa$pop_window, c(S$BURNIN_START, S$FILTER_END)),
       "F2 smoothing_frame_pop_args(): configured path, threaded year, burn-in..filter-end window")
yaml::write_yaml(list(simulator = list(backend = "spatpomp",
  filter = list(measurement = "occurrence_bernoulli"))), cfg_path)
S0 <- smoothing_make_config(list(config = cfg_path))
expect(identical(S0$POP_YEAR, 2020L) &&
       identical(S0$GID_POP_PATH, smoothing_default_config()$GID_POP_PATH),
       "F3 a config without the keys keeps the defaults")

# ===========================================================================
# G -- cache identity, deploy cross-arm check
# ===========================================================================
cat("\n## G: identity and refusals\n")
dm <- tempfile("py_cache_"); dir.create(dm)
exp_c <- smoothing_cache_meta("country")
yaml::write_yaml(exp_c, smoothing_cache_meta_path(dm))
expect(is.null(err_of(check_smoothing_cache_meta(dm, c(exp_c, list(pop_year = 2020L))))),
       "G1 legacy cache (no pop_year) accepted by a default run")
e <- err_of(check_smoothing_cache_meta(dm, c(exp_c, list(pop_year = "season_year"))))
expect(!is.null(e) && grepl("pop_year", e), "G2 legacy cache refused by a season_year run")
S_m <- list(N_REPS = 2L, NP = 10L, PARTITION_BY = "country", run_tag = "x", POP_YEAR = "season_year")
meta_m <- smoothing_cache_identity(S_m, tus = "a", optional = smoothing_cache_optional_fields(S_m))
expect(identical(meta_m$pop_year, "season_year"), "G3 08c writes pop_year when not 2020")
S_d <- modifyList(S_m, list(POP_YEAR = 2020L))
expect(is.null(smoothing_cache_identity(S_d, tus = "a",
                                        optional = smoothing_cache_optional_fields(S_d))$pop_year),
       "G4 08c writes no pop_year at the default")
yaml::write_yaml(meta_m, smoothing_cache_meta_path(dm))
expect(is.null(err_of(check_smoothing_cache_meta(dm, c(exp_c, list(pop_year = "season_year"))))),
       "G5 season_year cache accepted by a season_year run")
expect(!is.null(err_of(check_smoothing_cache_meta(dm, c(exp_c, list(pop_year = 2020L))))),
       "G6 season_year cache refused by a default run")
unlink(dm, recursive = TRUE)
tp <- function(v) list(simulator = list(location_frame = if (is.null(v)) list() else list(pop_year = v)))
expect(is.null(err_of(deploy_assert_sim_pop_year(list(historical = tp(NULL), forecast = tp(2020L),
                                                      counterfactual = tp(NULL))))),
       "G7 deploy: arms at the default agree")
e <- err_of(deploy_assert_sim_pop_year(list(historical = tp("season_year"), forecast = tp(NULL),
                                            counterfactual = tp("season_year"))))
expect(!is.null(e) && grepl("sim.shared", e), "G8 deploy: arms with different pop_year refused")

unlink(c(gp_path, wpp_path, cfg_path))
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
