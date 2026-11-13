# test_results_af_posterior.R -- the posterior interval of the attributable
# fraction / number (utils/results_derive.R: af_posterior_tables() and the
# af_posterior / event_af_posterior steps), and the bootstrap switch
# (results.derivation.attrib_boot_b = 0).
#
# HARD gate, quick-safe: no fit, no data, no simulation output.
#   [1] recovery: a synthetic fixed-draw design (99 draws x 50 replicates, two
#       units of one group, negative-binomial epidemic noise) with a KNOWN
#       per-draw AF: the noise-adjusted 95% interval recovers the true spread,
#       the raw one is wider, the noise share is in (0, 1);
#   [2] definitions: region AF_pool / AN_pool, Africa per draw from the
#       replicate-level unit totals (the "independent" combination kept side by
#       side), country rows, per-draw table, explosive guard and n_min;
#   [2b] two groups in one unit with an exact joint AF, and the
#       unit-level explosive guard;
#   [2c] a pure-noise scope is flagged unidentified (adjusted
#       interval NA), an identified one is unchanged;
#   [3] design refusals: random draws, not the cycle index, arms with different
#       maps, missing draws and too few replicates per draw; a warning when
#       max_draws shares a factor with the smoothing's replicate count;
#   [4] bootstrap off: attrib_metrics_tables() / attrib_exposure_tables() with
#       n_boot = 0 give NA intervals and the same point estimates;
#   [5] af_posterior_window_reps(): window sums per country and unit, the
#       explosive keys, and a v1 manifest (no unit column) keyed by group.
#
# Run: Rscript analysis/tests/test_results_af_posterior.R

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
cat("# posterior attributable fraction\n")
near <- function(a, b, tol) isTRUE(abs(a - b) <= tol)

#' A fixed-draw design: D draws x R replicates per unit, group g over units
#' u1 (country AAA) and u2 (BBB), group h on unit u3 (CCC). Per draw d, the
#' expected factual total mu_f[d] and the true AF af_true[d] are shared by the
#' units of a group (one parameter vector); epidemic noise is negative binomial.
synth <- function(D = 99L, R = 50L, nb_size = 2, seed = 1L) {
  set.seed(seed)
  af_true <- list(g = rnorm(D, 0.45, 0.12), h = rnorm(D, 0.10, 0.05))
  mu_f    <- list(g = exp(rnorm(D, log(3000), 0.3)), h = exp(rnorm(D, log(800), 0.3)))
  n <- D * R; rid <- seq_len(n); draw <- ((rid - 1L) %% D) + 1L
  units <- list(u1 = c("g", "AAA", 0.6), u2 = c("g", "BBB", 0.4), u3 = c("h", "CCC", 1))
  rows <- lapply(names(units), function(u) {
    grp <- units[[u]][1]; s <- as.numeric(units[[u]][3])
    f  <- rnbinom(n, mu = s * mu_f[[grp]][draw], size = nb_size)
    cf <- rnbinom(n, mu = s * mu_f[[grp]][draw] * (1 - af_true[[grp]][draw]), size = nb_size)
    tibble::tibble(cell = "w", enso_group = grp, unit = u, country = units[[u]][2], replicate_id = rid,
                   draw = draw, f = f, cf = cf)
  })
  p <- dplyr::bind_rows(rows)
  arm <- function(run, aset, v) dplyr::transmute(p, cell, enso_group, unit, country, run = run,
                                                 analog_set = aset, replicate_id, cases = .data[[v]],
                                                 onsets = as.numeric(.data[[v]] > 0))
  reps <- dplyr::bind_rows(arm("forecast", NA_character_, "f"), arm("counterfactual", "neutral", "cf"))
  list(reps = reps, explosive = reps[0, c("cell", "enso_group", "run", "analog_set", "replicate_id")],
       draws = dplyr::distinct(p, cell, unit, replicate_id, draw), af_true = af_true, mu_f = mu_f, p = p)
}
gs <- c(g = "G", h = "H")

# --- [1] recovery -------------------------------------------------------------
cat("# [1] recovery of a known posterior\n")
sy <- synth()
out <- af_posterior_tables(sy$reps, sy$explosive, sy$draws, "neutral", gs, guards = TRUE)
g <- dplyr::filter(out$scope, metric == "cases", level == "region", key == "g")
tq <- unname(quantile(sy$af_true$g, c(.025, .975)))
expect(nrow(g) == 1L && g$n_draws == 99L && g$n_min == 50L, "region g: 99 draws, 50 replicates per unit and draw")
expect(g$af_noise_share > 0 && g$af_noise_share < 1, sprintf("noise share in (0, 1): %.2f", g$af_noise_share))
expect(near(g$af_adj_lo, tq[1], 0.06) && near(g$af_adj_hi, tq[2], 0.06),
       sprintf("adjusted 95%% [%.2f, %.2f] recovers the true [%.2f, %.2f] (+/- 0.06)",
               g$af_adj_lo, g$af_adj_hi, tq[1], tq[2]))
expect(g$af_raw_lo < g$af_adj_lo && g$af_raw_hi > g$af_adj_hi, "the raw interval is wider (it carries the noise)")
expect(g$af_adj_lo <= g$af_adj_med && g$af_adj_med <= g$af_adj_hi && g$af_adj_q25 <= g$af_adj_q75,
       "adjusted quantiles ordered")
sy_q <- synth(nb_size = 50)   # little epidemic noise: raw ~ adjusted ~ truth
gq <- dplyr::filter(af_posterior_tables(sy_q$reps, sy_q$explosive, sy_q$draws, "neutral", gs, guards = TRUE)$scope,
                    metric == "cases", level == "region", key == "g")
expect(gq$af_noise_share < g$af_noise_share && near(gq$af_raw_hi, gq$af_adj_hi, 0.03),
       sprintf("less epidemic noise -> smaller share (%.2f < %.2f), raw ~ adjusted", gq$af_noise_share,
               g$af_noise_share))

# --- [2] definitions ----------------------------------------------------------
cat("# [2] definitions\n")
mean_by_draw <- sy$p %>% dplyr::group_by(enso_group, unit, draw) %>%
  dplyr::summarise(f = mean(f), cf = mean(cf), .groups = "drop") %>%
  dplyr::group_by(enso_group, draw) %>% dplyr::summarise(f = sum(f), cf = sum(cf), .groups = "drop")
mg <- dplyr::filter(mean_by_draw, enso_group == "g")
expect(near(g$af_pool, 1 - sum(mg$cf) / sum(mg$f), 1e-12), "region AF_pool = 1 - sum CF_d / sum F_d")
expect(near(g$an_pool, mean(mg$f - mg$cf), 1e-9), "region AN_pool = mean over draws of F_d - CF_d")
dg <- dplyr::filter(out$draws, metric == "cases", level == "region", key == "g") %>% dplyr::arrange(draw)
expect(isTRUE(all.equal(dg$af, 1 - mg$cf / mg$f)) && isTRUE(all.equal(dg$an, mg$f - mg$cf)),
       "per-draw table: AF_d, AN_d of the unit means summed at the same draw")
afr <- dplyr::filter(out$scope, metric == "cases", level == "continent")
gm <- mean_by_draw %>% dplyr::group_by(enso_group) %>% dplyr::summarise(f = mean(f), cf = mean(cf))
expect(near(afr$af_pool, 1 - sum(gm$cf) / sum(gm$f), 1e-12) && afr$n_draws == 99L && is.na(afr$n_mc),
       "Africa (joint): pool of the group means; n_draws = the design's draws, no n_mc")
expect(is.na(g$n_mc) && is.finite(g$af_noise_share_iqr) && g$af_noise_share_iqr > 0,
       "region rows: no n_mc; robust (IQR) noise share as a diagnostic column")
da <- dplyr::filter(out$draws, metric == "cases", level == "continent") %>% dplyr::arrange(draw)
ma <- mean_by_draw %>% dplyr::group_by(draw) %>% dplyr::summarise(f = sum(f), cf = sum(cf)) %>% dplyr::arrange(draw)
expect(nrow(da) == 99L && isTRUE(all.equal(da$af, 1 - ma$cf / ma$f)),
       "Africa (joint) per-draw rows: AF_d of the summed unit means at the same draw d")
oi <- af_posterior_tables(sy$reps, sy$explosive, sy$draws, "neutral", gs, guards = TRUE, africa = "independent")
ai <- dplyr::filter(oi$scope, metric == "cases", level == "continent")
expect(ai$n_mc == AF_POST_MC && nrow(dplyr::filter(oi$draws, metric == "cases", level == "continent")) == AF_POST_MC,
       "Africa (independent, kept side by side): AF_MC combinations")
expect(near(ai$af_adj_lo, afr$af_adj_lo, 0.03) && near(ai$af_adj_hi, afr$af_adj_hi, 0.03),
       "units of different groups are independent here: joint ~ independent Africa interval")
ctry <- dplyr::filter(out$scope, metric == "cases", level == "country")
expect(setequal(ctry$key, c("AAA", "BBB", "CCC")) && all(ctry$region[ctry$key == "AAA"] == "G"),
       "country rows keyed by ISO3, region label from group_short")
o2 <- af_posterior_tables(sy$reps, sy$explosive, sy$draws, "neutral", gs, guards = TRUE)
expect(identical(out, o2), "deterministic (fixed Africa seed)")
ex <- dplyr::distinct(dplyr::filter(sy$reps, enso_group == "g", replicate_id %in% 1:3, run == "forecast"),
                      cell, enso_group, run, analog_set, replicate_id)
og <- af_posterior_tables(sy$reps, ex, sy$draws, "neutral", gs, guards = c(TRUE, FALSE))
gg <- dplyr::filter(og$scope, metric == "cases", level == "region", key == "g")
expect(gg$n_min[gg$explosive_guard] == 49L && gg$n_min[!gg$explosive_guard] == 50L,
       "explosive guard drops replicates 1-3 of group g (draws 1-3 keep 49), not of h")

# --- [2b] dependent groups in one unit ---------------------------------------
cat("# [2b] joint Africa interval\n")
#' Groups g and h simulated together in unit u1: per replicate the factual is
#' 1000 in each group and the counterfactuals are 1000 (1 - a) and 1000 a, so
#' each group's AF spreads over (0.1, 0.9) while Africa's is exactly 0.5 in
#' every replicate (the issue's example). No epidemic noise.
dep <- function(D = 40L, R = 25L) {
  set.seed(3L)
  n <- D * R; rid <- seq_len(n); draw <- ((rid - 1L) %% D) + 1L
  a <- stats::runif(D, 0.1, 0.9)[draw]
  p <- dplyr::bind_rows(
    tibble::tibble(enso_group = "g", country = "AAA", f = 1000, cf = 1000 * (1 - a)),
    tibble::tibble(enso_group = "h", country = "CCC", f = 1000, cf = 1000 * a)) %>%
    dplyr::mutate(cell = "w", unit = "u1", replicate_id = rep(rid, 2L), draw = rep(draw, 2L))
  arm <- function(run, aset, v) dplyr::transmute(p, cell, enso_group, unit, country, run = run,
                                                 analog_set = aset, replicate_id, cases = .data[[v]],
                                                 onsets = as.numeric(.data[[v]] > 0))
  reps <- dplyr::bind_rows(arm("forecast", NA_character_, "f"), arm("counterfactual", "neutral", "cf"))
  list(reps = reps, explosive = reps[0, c("cell", "enso_group", "run", "analog_set", "replicate_id")],
       draws = dplyr::distinct(p, cell, unit, replicate_id, draw))
}
dy <- dep()
aj <- dplyr::filter(af_posterior_tables(dy$reps, dy$explosive, dy$draws, "neutral", gs, guards = TRUE)$scope,
                    metric == "cases", level == "continent")
an <- dplyr::filter(af_posterior_tables(dy$reps, dy$explosive, dy$draws, "neutral", gs, guards = TRUE,
                                        africa = "independent")$scope, metric == "cases", level == "continent")
expect(near(aj$af_raw_lo, 0.5, 1e-12) && near(aj$af_raw_hi, 0.5, 1e-12),
       sprintf("joint: Africa interval [%.3f, %.3f] = the exact 0.5 of every replicate", aj$af_raw_lo, aj$af_raw_hi))
expect(an$af_raw_hi - an$af_raw_lo > 0.3,
       sprintf("independent recombination invents a spread [%.2f, %.2f] (the defect)", an$af_raw_lo, an$af_raw_hi))
# Unit-level guard: replicate 1 explosive in group g also leaves group h of the same unit.
ex1 <- dplyr::distinct(dplyr::filter(dy$reps, enso_group == "g", replicate_id == 1L, run == "forecast"),
                       cell, enso_group, run, analog_set, replicate_id)
pu <- af_posterior_pairs(dy$reps, ex1, dy$draws, "neutral", "cases", guard = TRUE)
pg <- af_posterior_pairs(dy$reps, ex1, dy$draws, "neutral", "cases", guard = TRUE, guard_scope = "group")
expect(!any(pu$replicate_id == 1L) && any(pg$replicate_id == 1L & pg$enso_group == "h") &&
         !any(pg$replicate_id == 1L & pg$enso_group == "g"),
       "guard_scope unit drops the replicate from every group of its unit; group drops it from g only")

# --- [2c] unidentified posterior spread ---------------------------------------
cat("# [2c] unidentified spread\n")
#' Every draw has the same expected totals per arm (its unit-arm mean), so
#' the same true AF: the spread across draws is pure noise.
sy_n <- synth()
sy_n$reps <- sy_n$reps %>% dplyr::group_by(enso_group, unit, run) %>%
  dplyr::mutate(cases = stats::rnbinom(dplyr::n(), mu = mean(cases), size = 2)) %>% dplyr::ungroup()
on <- af_posterior_tables(sy_n$reps, sy_n$explosive, sy_n$draws, "neutral", gs, guards = TRUE,
                          noise_share_max = 0.8)
gn <- dplyr::filter(on$scope, metric == "cases", level == "region", key == "g")
expect(gn$af_noise_share >= 0.8 && identical(gn$af_identified, FALSE) && is.na(gn$af_adj_lo) && is.na(gn$af_adj_hi),
       sprintf("pure-noise scope (share %.2f): unidentified, adjusted interval NA", gn$af_noise_share))
expect(is.finite(gn$af_raw_lo) && is.finite(gn$af_raw_hi) && is.finite(gn$af_pool) && gn$af_var_between >= 0,
       "unidentified scope keeps its pooled value, raw interval and between-draw variance")
g8 <- dplyr::filter(af_posterior_tables(sy$reps, sy$explosive, sy$draws, "neutral", gs, guards = TRUE,
                                        noise_share_max = 0.8)$scope, metric == "cases", level == "region", key == "g")
expect(isTRUE(g8$af_identified) && near(g8$af_adj_lo, g$af_adj_lo, 1e-12) && near(g8$af_adj_hi, g$af_adj_hi, 1e-12),
       "an identified scope: same adjusted interval as without the threshold")
expect(!"af_identified" %in% names(out$scope), "noise_share_max = Inf (default of the core): no flag columns")

# --- [3] design refusals ------------------------------------------------------
cat("# [3] design refusals\n")
cfg <- function(mode = "cycle", S = 99L) list(simulator = list(max_draws = S, forecast = list(param_draw_mode = mode)))
d_ok <- dplyr::distinct(dplyr::filter(sy$draws, unit == "u1"), replicate_id, draw)
chk <- af_posterior_check_design("u1", list(factual = cfg(), neutral = cfg()), list(factual = d_ok, neutral = d_ok))
expect(attr(chk, "max_draws") == 99L && nrow(chk) == nrow(d_ok), "a cycle design passes, S = max_draws")
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
expect(err(af_posterior_check_design("u1", list(factual = cfg("random"), neutral = cfg()),
                                     list(factual = d_ok, neutral = d_ok))),
       "refused: an arm with param_draw_mode random")
d_rand <- dplyr::mutate(d_ok, draw = sample(draw))
expect(err(af_posterior_check_design("u1", list(factual = cfg(), neutral = cfg()),
                                     list(factual = d_rand, neutral = d_rand))),
       "refused: draws that are not the cycle index")
d_shift <- dplyr::mutate(d_ok, replicate_id = replicate_id + 99L, draw = ((replicate_id - 1L) %% 99L) + 1L)
expect(err(af_posterior_check_design("u1", list(factual = cfg(), neutral = cfg()),
                                     list(factual = d_ok, neutral = d_shift))),
       "refused: arms with different replicate -> draw maps")
expect(err(af_posterior_check_design("u1", list(factual = cfg(S = 99L), neutral = cfg(S = 199L)),
                                     list(factual = d_ok, neutral = d_ok))),
       "refused: arms with different max_draws")
warns <- function(expr) { w <- NULL; withCallingHandlers(expr, warning = function(x) {
  w <<- conditionMessage(x); invokeRestart("muffleWarning") }); w }
d200 <- tibble::tibble(replicate_id = 1:400, draw = ((1:400 - 1L) %% 200L) + 1L)
expect(!is.null(warns(af_posterior_check_design("u1", list(factual = cfg(S = 200L)), list(factual = d200),
                                                n_smoothing = 100L))),
       "warned: max_draws 200 shares a factor with 100 smoothing replicates")
expect(is.null(warns(af_posterior_check_design("u1", list(factual = cfg()), list(factual = d_ok),
                                               n_smoothing = 100L))),
       "no warning: max_draws 99 coprime to 100")
cfg_sh <- function(h) { cf <- cfg(S = 200L); cf$simulator$forecast$handoff_weighting <- h; cf }
expect(is.null(warns(af_posterior_check_design("u1", list(factual = cfg_sh("shuffled"), neutral = cfg_sh("shuffled")),
                                               list(factual = d200, neutral = d200), n_smoothing = 100L))),
       "no warning: max_draws 200 vs 100 under shuffled handoff")
expect(!is.null(warns(af_posterior_check_design("u1", list(factual = cfg_sh("uniform")), list(factual = d200),
                                                n_smoothing = 100L))),
       "warned: max_draws 200 vs 100 under an explicit uniform handoff")
pairs <- af_posterior_pairs(sy$reps, sy$explosive, sy$draws, "neutral", "cases")
expect(af_posterior_check_reps(pairs, sy$draws) == 50L, "check_reps: smallest count 50")
expect(err(af_posterior_check_reps(dplyr::filter(pairs, !(unit == "u2" & draw == 7L)), sy$draws)),
       "refused: a draw of the design missing from a unit")
expect(err(af_posterior_check_reps(dplyr::filter(pairs, !(unit == "u3" & draw == 5L & replicate_id > 5 * 99)),
                                   sy$draws, min_reps = 20L)),
       "refused: fewer than min_reps replicates for a draw")

# --- [4] bootstrap off --------------------------------------------------------
cat("# [4] attrib_boot_b = 0\n")
a0 <- as.Date("2026-05-01")
am <- expand.grid(replicate_id = 1:4, adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
                  run = c("forecast", "counterfactual"), stringsAsFactors = FALSE) %>%
  dplyr::mutate(month = a0, analog_set = ifelse(run == "forecast", NA_character_, "all"),
                enso_group = ifelse(substr(adm1_unit, 1, 3) == "KEN", "enso_east", "enso_horn"),
                group_horizon_end = as.Date("2026-06-01"),
                cases = ifelse(run == "counterfactual", 10, c(30, 12, 30, 0)[replicate_id]),
                outbreaks = as.integer(cases > 10)) %>% tibble::as_tibble()
cm <- am %>% dplyr::mutate(country = substr(adm1_unit, 1, 3)) %>%
  dplyr::group_by(replicate_id, month, country, run, analog_set, enso_group, group_horizon_end) %>%
  dplyr::summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
pop <- tibble::tibble(adm1_unit = c("KEN_1", "KEN_2", "UGA_1"), pop = c(100, 40, 60))
ge <- c(enso_east = "East", enso_horn = "Horn")
m50 <- attrib_metrics_tables(am, cm, pop, ge, "all", a0, 20, 0.05, c(1, 2), n_boot = 50L, cap = 2e5)
m0  <- attrib_metrics_tables(am, cm, pop, ge, "all", a0, 20, 0.05, c(1, 2), n_boot = 0L, cap = 2e5)
expect(all(is.na(m0$af_scope$af_lo)) && all(is.na(m0$af_scope$an_hi)) && all(is.na(m0$pr_scope$pr_hi)),
       "n_boot 0: AF / AN / PR intervals NA")
expect(identical(m0$af_scope$af_of_means, m50$af_scope$af_of_means) && identical(m0$pr_scope$pr, m50$pr_scope$pr) &&
         identical(m0$sn_scope$sn, m50$sn_scope$sn),
       "n_boot 0: AF / PR / S/N point estimates unchanged")
x50 <- suppressMessages(attrib_exposure_tables(am, cm, pop, ge, "all", a0, 0.25, 0.66, 0.8, n_boot = 50L, cap = 2e5))
x0  <- suppressMessages(attrib_exposure_tables(am, cm, pop, ge, "all", a0, 0.25, 0.66, 0.8, n_boot = 0L, cap = 2e5))
expect(identical(x0$exposure$pop, x50$exposure$pop) && identical(x0$adm1$class, x50$adm1$class),
       "n_boot 0: exposure point populations and IPCC classes unchanged")

# --- [5] window totals --------------------------------------------------------
cat("# [5] af_posterior_window_reps()\n")
am2 <- dplyr::bind_rows(dplyr::mutate(am, unit = "U1"),
                        dplyr::mutate(am, unit = "U1", month = as.Date("2026-04-01"), cases = 1000))
cm2 <- dplyr::bind_rows(cm, dplyr::mutate(cm, month = as.Date("2026-04-01"), n_cases = 1000))
wr <- af_posterior_window_reps(am2, cm2, "attribution", a0, cap = 2e5)
ken1 <- dplyr::filter(wr$reps, country == "KEN", run == "forecast", replicate_id == 1L)
expect(nrow(ken1) == 1L && ken1$cases == 30 + 30 && ken1$unit == "U1" && ken1$cell == "attribution",
       "country window total sums its ADM1s from the window start (the earlier month excluded)")
expect(nrow(wr$explosive) == 0L, "no explosive replicate in the window")
wr1 <- af_posterior_window_reps(am, cm, "w", a0, cap = 2e5)
expect(identical(sort(unique(wr1$reps$unit)), c("enso_east", "enso_horn")), "no unit column (v1): unit = group")
cells2 <- dplyr::bind_rows(attrib_cell_table(a0, cap = 2e5),
                           dplyr::mutate(attrib_cell_table(as.Date("2026-04-01"), cap = 2e5), window = "earlier"))
st <- af_posterior_stash_build(list(groups = c("enso_east", "enso_horn"), adm1_month = am2, country_month = cm2),
                               cells2)
expect(setequal(unique(st$reps$cell), c("attribution", "earlier")) &&
         max(st$reps$cases[st$reps$cell == "earlier"]) > max(st$reps$cases[st$reps$cell == "attribution"]),
       "stash: one window-total set per cell")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
