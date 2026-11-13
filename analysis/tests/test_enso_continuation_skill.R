# test_enso_continuation_skill.R -- skill by lead of the ENSO continuation
# (utils/enso_continuation_skill.R; written by
# `03_prepare_enso_continuation_library.R --validate`). HARD gate, quick,
# synthetic inputs only (no library file, no ONI cache, no network).
#
#   [1] fair CRPS: a degenerate (all-equal) ensemble and a one-member forecast
#       score |error|; the pairwise term lowers a spread ensemble's score.
#   [2] seams: history (6 months) and horizon respected; seam state features.
#   [3] lead indexing: an oracle method that returns obs[s + k] at lead k
#       scores 0 at every lead; a one-month-late oracle does not; ar1 rows are
#       phi^k * ONI(s).
#   [4] references: CRPSS of a reference against itself is 0 (per seed and
#       over seeds); climatology draws only same-calendar-month values outside
#       the +- exclude window.
#   [5] analog method: [H x n_ens] from continueEnsoAnalog(), seed-reproducible;
#       ensoContinuationLibraryRehorizon() recomputes the candidate rule.
#   [6] output schema of ensoContinuationSkill() (columns, one row per method x
#       lead, n_seeds, seed range ordered).
#   [7] settings: defaults' simulator.enso.continuation, overlaid by a deploy
#       YAML's deploy.enso_continuation (.inf -> Inf, exclude_models).
#
# Run: Rscript analysis/tests/test_enso_continuation_skill.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
source("analysis/utils/obs_sim_compare.R")
source("analysis/utils/enso_continuation_skill.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

# --- synthetic observed ONI ----------------------------------------------------
obs_d <- seq(as.Date("1960-01-01"), as.Date("2000-12-01"), by = "month")
set.seed(3)
obs   <- as.numeric(stats::arima.sim(list(ar = 0.9), n = length(obs_d), sd = 0.3))
obs_m <- as.integer(format(obs_d, "%m"))
H <- 6L

cat("\n[1] fair CRPS\n")
expect(isTRUE(all.equal(crps_fair(rep(0.7, 25), -0.4), 1.1)), "degenerate ensemble: CRPS = |error|")
expect(isTRUE(all.equal(crps_fair(0.7, -0.4), 1.1)), "one-member forecast: CRPS = |error|")
expect(crps_fair(c(-0.4, 0.6, 0.8, 1.0), -0.4) < mean(abs(c(-0.4, 0.6, 0.8, 1.0) + 0.4)),
       "spread ensemble: pairwise term lowers the score")

cat("\n[2] seams\n")
seams <- ensoSkillSeams(obs_d, "1960-01-01", "2000-12-01", horizon = H)
expect(min(seams) == 7L && max(seams) == length(obs_d) - H,
       "seams need 6 months of history and H months of future")
s <- ensoSkillSeams(obs_d, "1970-01-01", "1970-12-01", horizon = H)
expect(length(s) == 12L && all(format(obs_d[s], "%Y") == "1970"), "window restricts the seams")
st <- ensoSkillSeamState(obs, obs_m, 100L)
expect(identical(st, list(mon = obs_m[100], v = obs[100], tend = obs[100] - obs[97], phase = obs[94])),
       "seam state: v, tend = v - ONI(s-3), phase = ONI(s-6), calendar month")
expect(!is.null(err_of(ensoContinuationHindcast(obs, obs_d, obs_m, list(n = ensoSkillMethodNeutral()),
                                                length(obs) - 2L, horizon = H))),
       "a seam without H months of future is refused")

cat("\n[3] lead indexing\n")
oracle <- function(shift) function(ctx) matrix(ctx$obs[ctx$s + seq_len(ctx$H) + shift], ctx$H, 3L)
s_test <- seams[seq(1, length(seams), by = 37)]
sc <- ensoContinuationHindcast(obs, obs_d, obs_m,
                               list(oracle = oracle(0L), late = oracle(1L), ar1 = ensoSkillMethodAr1(0.9)),
                               s_test[s_test + H + 1L <= length(obs)], horizon = H, progress = 0L)
expect(all(sc$crps[sc$method == "oracle"] == 0), "oracle (obs[s + k] at lead k) scores 0 at every lead")
expect(all(sc$crps[sc$method == "late"] > 0), "one-month-late oracle scores > 0")
expect(identical(sort(unique(sc$lead)), seq_len(H)), "leads 1..H")
a1 <- sc[sc$method == "ar1", ]
expect(isTRUE(all.equal(a1$mean, 0.9^a1$lead * a1$x0)) && all(a1$sd == 0), "ar1 lead k = phi^k * ONI(s), sd 0")
o1 <- sc[sc$method == "oracle", ]
expect(isTRUE(all.equal(o1$y, obs[match(o1$seam, obs_d) + o1$lead])), "y at lead k = observed ONI(s + k)")

cat("\n[4] references\n")
clim <- ensoSkillMethodClim(n_ens = 50L, pool_from = "1960-01-01", exclude = 24L)
set.seed(9)
ctx <- list(s = 200L, x0 = obs[200], state = ensoSkillSeamState(obs, obs_m, 200L), H = H,
            obs = seq_along(obs) + 0, obs_month = obs_m, obs_date = obs_d)   # values = indices
e <- clim(ctx)
expect(identical(dim(e), c(H, 50L)), "climatology [H x n_ens]")
expect(all(abs(e - 200) > 24), "climatology leaves out +- exclude months around the seam")
expect(all(vapply(seq_len(H), function(k) all(obs_m[e[k, ]] == ((obs_m[200] - 1L + k) %% 12L) + 1L), TRUE)),
       "climatology at lead k draws the target calendar month only")
meths <- list(neutral = ensoSkillMethodNeutral(), clim = clim, ar1 = ensoSkillMethodAr1(0.9))
sc2 <- ensoContinuationHindcast(obs, obs_d, obs_m, meths, s_test, horizon = H, seed = 1L, progress = 0L)
by <- ensoSkillByLead(sc2, references = c("clim", "neutral"))
expect(all(by$crpss_vs_clim[by$method == "clim"] == 0) && all(by$crpss_vs_neutral[by$method == "neutral"] == 0),
       "a reference's CRPSS against itself is 0")
expect(all(is.na(by$r_seam[by$method == "neutral"])), "r_seam NA for a constant forecast")
set.seed(1); sc3 <- ensoContinuationHindcast(obs, obs_d, obs_m, meths, s_test, horizon = H, progress = 0L)
expect(identical(sc2, sc3), "seed = s equals set.seed(s) before the call (one RNG stream)")

cat("\n[5] analog method\n")
mk <- function(n) as.numeric(stats::arima.sim(list(ar = 0.9), n = n, sd = 0.3))
set.seed(11); xa <- mk(400); xb <- mk(300)
lib <- list(series = list(A = xa, B = xb),
            features = rbind(ensoContinuationFeatures(xa, rep_len(1:12, 400), "A", 24L),
                             ensoContinuationFeatures(xb, rep_len(1:12, 300), "B", 24L)))
okk <- lib$features$ok
lib$feature_sd <- vapply(c("v", "tend", "phase"), function(f) sd(lib$features[[f]][okk]), numeric(1))
lib$meta <- list(horizon_months = 24L)
cand <- lib$features[okk, c("model", "i", "mon", "v", "tend", "phase")]
lib$by_month <- split(cand, factor(cand$mon, levels = 1:12))
re <- ensoContinuationLibraryRehorizon(lib, 24L)
expect(identical(re$features, lib$features) && isTRUE(all.equal(re$feature_sd, lib$feature_sd)) &&
         identical(re$by_month, lib$by_month), "rehorizon to the same horizon is the identity")
re12 <- ensoContinuationLibraryRehorizon(lib, 12L)
expect(sum(re12$features$ok) == sum(okk) + 2L * 12L && re12$meta$horizon_months == 12L,
       "rehorizon to 12 months adds 12 candidate starts per series")
an <- ensoSkillMethodAnalog(lib, n_ens = 20L, k = 5L, blend_months = 3, min_sep = 24L,
                            exclude_models = "B")
ctx$obs <- obs
set.seed(4); m1 <- an(ctx); set.seed(4); m2 <- an(ctx)
set.seed(4)
m3 <- continueEnsoAnalog(lib, rep(ctx$state$v, 20L), rep(ctx$state$tend, 20L), rep(ctx$state$phase, 20L),
                         ctx$state$mon, H, k = 5L, blend_months = 3, min_sep = 24L, exclude_models = "B")
expect(identical(dim(m1), c(H, 20L)) && identical(m1, m2) && identical(m1, m3),
       "analog = continueEnsoAnalog() on the seam state, [H x n_ens], seed-reproducible")

cat("\n[6] output schema\n")
meths$analog <- ensoSkillMethodAnalog(lib, n_ens = 20L, k = 5L)
sk <- ensoContinuationSkill(obs, obs_d, obs_m, meths, s_test, horizon = H, seeds = c(1L, 2L))
want_cols <- c("method", "lead", "n_seams", "n_seeds", "crps", "rmse", "spread", "r_seam",
               "crps_clim", "crpss_vs_clim", "crpss_vs_clim_min", "crpss_vs_clim_max",
               "crps_neutral", "crpss_vs_neutral", "crpss_vs_neutral_min", "crpss_vs_neutral_max")
expect(identical(names(sk), want_cols), "skill table columns")
expect(nrow(sk) == 4L * H && !anyDuplicated(paste(sk$method, sk$lead)), "one row per method x lead")
expect(all(sk$n_seeds == 2L) && all(sk$n_seams == length(s_test)), "n_seeds, n_seams")
expect(all(sk$crpss_vs_clim[sk$method == "clim"] == 0), "clim vs itself 0 over seeds")
expect(all(sk$crpss_vs_clim_min <= sk$crpss_vs_clim_max), "seed range ordered")
one <- ensoSkillByLead(ensoContinuationHindcast(obs, obs_d, obs_m, meths, s_test, horizon = H, seed = 1L,
                                                progress = 0L))
sk1 <- ensoContinuationSkill(obs, obs_d, obs_m, meths, s_test, horizon = H, seeds = 1L)
expect(isTRUE(all.equal(sk1$crps, one$crps)) && isTRUE(all.equal(sk1$crpss_vs_clim, one$crpss_vs_clim)),
       "one seed: the skill table equals ensoSkillByLead() of that seed")

cat("\n[7] settings\n")
d <- ensoContinuationSkillSettings("analysis/config_defaults.yml")
expect(identical(d$exclude_models, "CM2.1-NMME") && is.infinite(d$max_per_model) && d$k == 30L &&
         d$min_sep == 24L && d$blend_months == 3, "defaults: k 30, blend 3, min_sep 24, no cap, CM2.1-NMME out")
dep <- tempfile(fileext = ".yml")
yaml::write_yaml(list(deploy = list(enso_continuation = list(method = "analog", k = 10L,
                                                             exclude_models = list()))), dep)
o <- ensoContinuationSkillSettings("analysis/config_defaults.yml", dep)
expect(o$method == "analog" && o$k == 10L && length(o$exclude_models) == 0L && o$min_sep == 24L &&
         o$source == dep, "deploy block overlays the defaults leaf by leaf")
yaml::write_yaml(list(deploy = list(other = 1)), dep)
expect(!is.null(err_of(ensoContinuationSkillSettings("analysis/config_defaults.yml", dep))),
       "a deploy YAML without deploy.enso_continuation is refused")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
