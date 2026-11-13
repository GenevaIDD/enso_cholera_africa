# test_enso_counterfactual_match.R — counterfactual ENSO analogs.
#
# HARD gate, quick, synthetic data only (no data/, no generated_data/).
#   [1] the legacy builder's helpers in utils/enso.R equal frozen
#       copies of the in-script originals of 03_prepare_lenssen_analogs.R
#       (trajectories, rescale, init indices) on synthetic piControl series
#   [2]-[5] neighbours, paths, library matcher, cache, seam, generators
#   [6] launch-month seam pin (seam_at): "launch_minus_1" is the
#       first library builder bit for bit (spec, trajectories, digest);
#       "launch" pins column 1 to the observed launch ONI (seam error 0), fades
#       the offset from there, keeps the selection, and enters the identity
#
# Run: Rscript analysis/tests/test_enso_counterfactual_match.R

suppressPackageStartupMessages({ library(dplyr); library(lubridate) })
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

# --- frozen copies of 03_prepare_lenssen_analogs.R (at c032a12) -------------
old_select_analogs <- function(nino34, timeMap, target_oni, target_month,
                            n_analogs, min_sep, lead_months) {
  same_season_idx <- which(timeMap[, 2L] == target_month)
  last_valid_idx  <- length(nino34) - lead_months
  same_season_idx <- same_season_idx[same_season_idx <= last_valid_idx]
  if (length(same_season_idx) == 0L) {
    warning("No valid analog months for model (target_month=", target_month, ")")
    return(integer(0))
  }
  dist_vec <- abs(nino34[same_season_idx] - target_oni)
  ranked   <- same_season_idx[order(dist_vec)]
  selected <- integer(0)
  for (idx in ranked) {
    if (length(selected) >= n_analogs) break
    if (length(selected) == 0L || all(abs(idx - selected) >= min_sep)) {
      selected <- c(selected, idx)
    }
  }
  selected
}
old_extract_trajectories <- function(nino34, init_indices, lead_months) {
  purrr::map(init_indices, function(idx) {
    nino34[idx + seq(0L, lead_months)]
  })
}
old_rescale_trajectories <- function(trajectories, model_nino34, obs_sd) {
  model_sd     <- sd(model_nino34, na.rm = TRUE)
  scale_factor <- obs_sd / model_sd
  list(
    trajectories_adj = lapply(trajectories, function(traj) traj * scale_factor),
    scale_factor     = scale_factor
  )
}

# --- synthetic piControl series ---------------------------------------------
set.seed(8)
mk_series <- function(n_years, sd_scale) {
  n <- 12L * n_years
  x <- as.numeric(stats::filter(rnorm(n, 0, 0.25), 0.93, method = "recursive")) * sd_scale
  list(nino34 = x, timeMap = cbind(year = rep(seq_len(n_years), each = 12L), month = rep(1:12, n_years)))
}
models <- list(A = mk_series(300, 1.0), B = mk_series(200, 1.6), C = mk_series(150, 0.7))

cat("# [1] 8.1 legacy helpers moved verbatim\n")
expect(identical(LENSSEN_MODELS, c("CanESM5", "CCSM4-NMME", "CESM1-1-CAM5-CMIP5", "CESM1-NMME",
                                   "CESM2", "CM2.1-NMME", "CM2.5-NMME", "CM4", "ESM4", "GISS-E21G", "MIROC6")),
       "LENSSEN_MODELS = the 11 zip names, in the builders' order")
for (m in names(models)) {
  e <- models[[m]]
  for (cfg in list(list(t = 0.24, mon = 9L, n = 50L, sep = 12L, lead = 60L),
                   list(t = -1.2, mon = 5L, n = 20L, sep = 24L, lead = 120L))) {
    a <- old_select_analogs(e$nino34, e$timeMap, cfg$t, cfg$mon, cfg$n, cfg$sep, cfg$lead)
    b <- .lenssenSelectAnalogsLegacy(e$nino34, e$timeMap, cfg$t, cfg$mon, cfg$n, cfg$sep, cfg$lead)
    expect(identical(a, b), sprintf("model %s target %.2f: init indices identical", m, cfg$t))
    ta <- old_extract_trajectories(e$nino34, a, cfg$lead)
    tb <- .lenssenExtractTrajectoriesLegacy(e$nino34, b, cfg$lead)
    expect(identical(ta, tb), sprintf("model %s: trajectories identical", m))
    expect(identical(old_rescale_trajectories(ta, e$nino34, 0.7858),
                     .lenssenRescaleLegacy(tb, e$nino34, 0.7858)),
           sprintf("model %s: rescale identical (trajectories + factor)", m))
  }
}

# --- [2] neighbours(horizon =), faded-join paths, seam state -------------
cat("# [2] 8.2 ensoAnalogNeighbours(horizon =), ensoAnalogPaths, continueEnsoAnalog RNG-identical\n")
old_neighbours <- function(lib, state, k = 30L, min_sep = 24L,
                           feats = c("v", "tend", "phase"),
                           max_per_model = Inf, exclude_models = character(0)) {
  cand <- lib$by_month[[as.character(state$mon)]]
  if (length(exclude_models)) {
    cand <- cand[!cand$model %in% exclude_models, , drop = FALSE]
  }
  max_per_model <- as.numeric(max_per_model %||% Inf)
  d <- 0
  for (f in feats) d <- d + ((cand[[f]] - state[[f]]) / lib$feature_sd[[f]])^2
  cand <- cand[order(d), c("model", "i")]
  sel <- integer(0)
  for (r in seq_len(nrow(cand))) {
    if (length(sel) >= k) break
    same <- cand$model[sel] == cand$model[r]
    if (sum(same) < max_per_model &&
        (!any(same) || all(abs(cand$i[sel][same] - cand$i[r]) >= min_sep))) sel <- c(sel, r)
  }
  cand[sel, , drop = FALSE]
}
old_continue <- function(lib, seam_v, seam_tend, seam_phase, seam_month, n_months,
                         k = 30L, blend_months = 3, min_sep = 24L,
                         max_per_model = Inf, exclude_models = character(0)) {
  n_rep <- length(seam_v)
  fade <- if (is.infinite(blend_months)) rep(1, n_months) else exp(-seq_len(n_months) / blend_months)
  key  <- sprintf("%.10g|%.10g|%.10g", seam_v, seam_tend, seam_phase)
  nn_cache <- list()
  out <- matrix(NA_real_, nrow = n_months, ncol = n_rep)
  for (r in seq_len(n_rep)) {
    nn <- nn_cache[[key[r]]]
    if (is.null(nn)) {
      nn <- old_neighbours(lib, list(mon = seam_month, v = seam_v[r], tend = seam_tend[r],
                                     phase = seam_phase[r]), k = k, min_sep = min_sep,
                           max_per_model = max_per_model, exclude_models = exclude_models)
      nn_cache[[key[r]]] <- nn
    }
    p <- sample.int(nrow(nn), 1L)
    s <- lib$series[[nn$model[p]]]; i <- nn$i[p]
    out[, r] <- s[i + seq_len(n_months)] + (seam_v[r] - s[i]) * fade
  }
  out
}

# A synthetic continuation library: three models with unequal raw scales,
# rescaled to one sd, features with a 48-month candidate rule (as the deployed one).
mk_lib <- function(models, H = 48L, sd_obs = 0.8175) {
  series <- list(); feats <- list()
  for (m in names(models)) {
    x <- models[[m]]$nino34; mon <- models[[m]]$timeMap[, 2L]
    x <- x * sd_obs / stats::sd(x)
    series[[m]] <- x
    feats[[m]] <- ensoContinuationFeatures(x, mon, model = m, horizon = H)
  }
  features <- do.call(rbind, feats); rownames(features) <- NULL
  ok <- features$ok
  lib <- list(series = series, features = features,
              feature_sd = vapply(c("v", "tend", "phase"), function(f) stats::sd(features[[f]][ok]), numeric(1)),
              meta = list(horizon_months = H, n_months = lapply(series, length),
                          rescale = list(window = "1991_2020", obs_sd = sd_obs)))
  lib$meta$digest <- ensoContinuationLibraryDigest(lib)
  cand <- lib$features[lib$features$ok, c("model", "i", "mon", "v", "tend", "phase")]
  lib$by_month <- split(cand, factor(cand$mon, levels = 1:12))
  lib
}
lib <- mk_lib(models)
st <- list(mon = 8L, v = 0.4, tend = 0.3, phase = -0.2)
for (args in list(list(k = 30L, min_sep = 24L), list(k = 10L, min_sep = 12L, max_per_model = 3),
                  list(k = 20L, min_sep = 24L, exclude_models = "B"))) {
  a <- do.call(old_neighbours, c(list(lib, st), args))
  b <- do.call(ensoAnalogNeighbours, c(list(lib, st), args))
  expect(identical(a, b), sprintf("ensoAnalogNeighbours(horizon = NULL) == frozen (k %d)", args$k))
}
nn_h <- ensoAnalogNeighbours(lib, st, k = Inf, min_sep = 12L, horizon = 400L)
expect(all(nn_h$i + 400L <= vapply(lib$series, length, integer(1))[nn_h$model]),
       "horizon = 400: every start has 400 months of future")
expect(nrow(nn_h) < nrow(ensoAnalogNeighbours(lib, st, k = Inf, min_sep = 12L)),
       "horizon = 400 drops late starts")
seam_v <- c(0.4, 0.4, 1.1); seam_t <- c(0.3, 0.3, 0.5); seam_p <- c(-0.2, -0.2, 0.1)
for (bm in c(3, Inf)) {
  set.seed(42); a <- old_continue(lib, seam_v, seam_t, seam_p, 8L, 24L, k = 10L, blend_months = bm)
  set.seed(42); b <- continueEnsoAnalog(lib, seam_v, seam_t, seam_p, 8L, 24L, k = 10L, blend_months = bm)
  expect(identical(a, b), sprintf("continueEnsoAnalog(blend %s) identical to the frozen copy, same seed", bm))
}
nn <- ensoAnalogNeighbours(lib, st, k = 5L)
P <- ensoAnalogPaths(lib, nn, 0.4, 12L, 3)
s1 <- lib$series[[nn$model[1]]]; i1 <- nn$i[1]
expect(isTRUE(all.equal(P[1, ], s1[i1 + 1:12] + (0.4 - s1[i1]) * exp(-(1:12) / 3))),
       "ensoAnalogPaths: path(j) = analog(i + j) + (seam - analog(i)) rho^j")
expect(isTRUE(all.equal(ensoAnalogPaths(lib, nn, 0.4, 1L, 1e-9)[, 1], vapply(seq_len(nrow(nn)), function(r)
  lib$series[[nn$model[r]]][nn$i[r] + 1L], numeric(1)))), "blend -> 0: the analog itself")
oni_lkp <- stats::setNames(c(0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1.0),
                           format(seq(as.Date("2026-01-01"), by = "1 month", length.out = 8)))
ss <- ensoSeamState(oni_lkp, "2026-08-01")
expect(isTRUE(all.equal(ss, list(mon = 8L, v = 1.0, tend = 0.4, phase = 0.2))),
       "ensoSeamState: v = ONI(m), tend = v - ONI(m-3), phase = ONI(m-6)")
expect(!is.null(errmsg(ensoSeamState(oni_lkp, "2026-05-01"))), "ensoSeamState: missing ONI stops")

# --- [3] library matcher, multi-winter sets, read winters, cache identity --
cat("# [3] 8.3 library matcher and cache\n")
# Unequal raw model scales (B 1.6x, C 0.7x) are exactly what the legacy matcher
# gets wrong: it matches the raw series, then rescales.
oni_dates <- seq(as.Date("2024-01-01"), as.Date("2026-08-01"), by = "1 month")
oni <- data.frame(date = oni_dates, enso_index = round(0.9 * sin(seq_along(oni_dates) / 5), 3))
launch <- as.Date("2026-05-01")
oni_lkp <- stats::setNames(oni$enso_index, format(oni$date))
st <- ensoSeamState(oni_lkp, "2026-04-01")
pool <- buildLenssenCounterfactualPool(lib, st, n_months = 30L, n_per_model = 15L, min_sep = 12L,
                                       exclude_models = "C", blend_months = 3)
tab <- table(pool$model)
expect(all(tab <= 15L) && !"C" %in% names(tab), "per-model quota 15 and exclusion of C respected")
sep_ok <- all(vapply(split(pool$init_index, pool$model), function(i) all(diff(sort(i)) >= 12L), logical(1)))
expect(sep_ok, "within-model starts at least min_sep = 12 months apart")
expect(all(lib$features$mon[match(paste(pool$model, pool$init_index),
                                  paste(lib$features$model, lib$features$i))] == 4L),
       "every start is in the seam's calendar month (April)")
nn <- data.frame(model = pool$model, i = pool$init_index)
expect(identical(pool$traj, ensoAnalogPaths(lib, nn, st$v, 30L, 3)),
       "pool trajectories = faded-join paths from the seam (column 1 = the launch month = path(1))")
expect(identical(buildLenssenCounterfactualPool(lib, st, 30L, 15L, 12L, "C", 3), pool), "deterministic")
# Level matching: library pool's seam-month level is near the target; the legacy
# matcher (raw distance, rescale after) is biased for scaled models.
bias <- function(m, legacy) {
  if (legacy) {
    e <- models[[m]]
    idx <- .lenssenSelectAnalogsLegacy(e$nino34, e$timeMap, st$v, 4L, 15L, 12L, 30L)
    return(abs(mean(e$nino34[idx] * 0.8175 / sd(e$nino34)) - st$v))
  }
  r <- which(pool$model == m)
  abs(mean(vapply(r, function(k) lib$series[[m]][pool$init_index[k]], numeric(1))) - st$v)
}
b_leg <- vapply(c("A", "B"), bias, numeric(1), legacy = TRUE)
b_lib <- vapply(c("A", "B"), bias, numeric(1), legacy = FALSE)
expect(max(b_lib) < max(b_leg),
       sprintf("per-model seam-level bias: library %s < legacy %s (raw matching, rescale after)",
               paste(round(b_lib, 3), collapse = "/"), paste(round(b_leg, 3), collapse = "/")))

# Multi-winter classification.
tm <- matrix(0, nrow = 4, ncol = 30)   # launch 2026-05: DJF 2026/27 = cols 8:10, DJF 2027/28 = cols 20:22
tm[1, 8:10] <- 1.0; tm[2, 8:10] <- -1.0; tm[3, 8:10] <- 0.0; tm[4, 8:10] <- 0.1
tm[4, 20:22] <- 1.2                    # neutral first winter, El Nino second
one <- classifyEnsoAnalogSets(tm, launch, peak_season_year = 2026L)
expect(identical(classifyEnsoAnalogSetsMulti(tm, launch, primary_years = 2026L), one),
       "no other winter: classifyEnsoAnalogSetsMulti == classifyEnsoAnalogSets")
mw <- classifyEnsoAnalogSetsMulti(tm, launch, primary_years = 2026L, other_years = 2027L)
expect(identical(mw$labels, c("elnino", "nina", "neutral", "mixed")),
       "other_winters neutral: the row whose DJF 2027/28 is El Nino is 'mixed', the others keep their label")
mw_any <- classifyEnsoAnalogSetsMulti(tm, launch, primary_years = 2026L, other_years = 2027L, other_winters = "any")
expect(identical(mw_any$labels, one$labels), "other_winters any: primary labels only")
expect(identical(dim(mw$stat_matrix), c(4L, 2L)) && identical(colnames(mw$stat_matrix), c("2026", "2027")),
       "stat_matrix per winter")

cat("# [3b] read winters / read end\n")
expect(identical(as.integer(ensoReadWinters("2026-05-01", "2027-06-01")), 2026L),
       "forecast launch 2026-05, reads to 2027-06: DJF 2026/27 only")
expect(identical(as.integer(ensoReadWinters("2026-05-01", "2028-06-01")), c(2026L, 2027L)),
       "reads to 2028-06 (West, Africa window to 2027-10): DJF 2026/27 and 2027/28")
rw15 <- ensoReadWinters("2014-09-01", "2017-07-01")
expect(identical(as.integer(rw15), c(2014L, 2015L, 2016L)),
       "elnino_2015 (launch 2014-09, reads to 2017-07): DJF 2014/15, 2015/16, 2016/17")
expect(identical(attr(ensoReadWinters("2015-01-01", "2016-03-01"), "partial"), 2014L),
       "a winter begun before the launch is partial (not classifiable)")
win <- data.frame(group = c("enso_west", "enso_east"), unit = c("hmmW", "hmmW"),
                  end = as.Date(c("2027-10-01", "2028-01-01")))
lags <- list(enso_west = c(-6L, -4L, -8L), enso_east = c(0L, -1L, 2L))
expect(identical(ensoReadEnd(win, lags, "direct"), as.Date("2028-06-01")),
       "direct: West 2027-10 - (-8) = 2028-06 (East 2028-01 + 1 = 2028-02)")
expect(identical(ensoReadEnd(win, lags, "unit"), as.Date("2028-09-01")),
       "unit: unit end 2028-01 with the unit's min lag -8 = 2028-09")

cat("# [3c] cache spec / digest / status / provenance\n")
# seam_at launch_minus_1: the first library builder (the identities below were
# written for it; the launch pin is section [6]).
cfg <- list(match = "library", n_analogs = 15L, min_sep_months = 12L, winters = "all_read",
            seam_at = "launch_minus_1")
spec <- lenssenCacheSpec(launch, 30L, cfg, lib, st, primary_years = 2026L, other_years = 2027L,
                         read_end = "2028-06-01", exclude_models = "C")
tmp <- tempfile("lenssen_"); dir.create(tmp); path <- file.path(tmp, "lenssen_nino34_analogs.rds")
cache <- buildLenssenCache(spec, oni, lib, out_path = path)
expect(file.exists(path) && file.exists(file.path(tmp, "lenssen_nino34_analogs_meta.yml")), "cache + meta written")
expect(identical(cache$oni_trajectories, pool$traj), "cache trajectories = the pool")
expect(all(c("oni_trajectories", "n_months", "n_analogs", "set_labels", "set_stat", "model", "init_index", "meta")
           %in% names(cache)), "legacy cache contract kept")
cache2 <- buildLenssenCache(spec, oni, lib)
expect(identical(cache2$meta$digest, cache$meta$digest) && identical(cache2$set_labels, cache$set_labels),
       "deterministic: same digest and labels")
stt <- lenssenCacheStatus(path, spec, sets = c("all", "neutral"), min_set_size = 5L)
expect(isTRUE(stt$ok) && identical(stt$digest, cache$meta$digest), "status ok for its own spec")
spec_b <- spec; spec_b$n_analogs <- 20L
expect(!isTRUE(lenssenCacheStatus(path, spec_b)$ok) &&
         grepl("n_analogs", lenssenCacheStatus(path, spec_b)$reason), "spec mismatch -> not ok, naming the field")
stt_big <- lenssenCacheStatus(path, spec, sets = "neutral", min_set_size = 10000L)
expect(isTRUE(stt_big$ok) && identical(stt_big$small, "neutral"), "a set under min_set_size is reported")
legacy_path <- file.path(tmp, "legacy.rds")
saveRDS(list(oni_trajectories = matrix(0, 2, 2), meta = list(horizon_start = "2026-05-01")), legacy_path)
expect(!isTRUE(lenssenCacheStatus(legacy_path, spec)$ok), "a legacy cache never satisfies a library spec")
expect(is.null(lenssenCacheProvenance(legacy_path)), "provenance of a legacy cache: NULL (nothing recorded)")
pv <- lenssenCacheProvenance(path)
expect(identical(pv$digest, cache$meta$digest) && identical(pv$match, "library"), "provenance carries the digest")
ev <- list(id = "elnino_2015", launch = as.Date("2014-09-01"), peak_season_year = 2015L)
old_cmd <- sprintf(paste("Rscript analysis/03_prepare_lenssen_analogs.R --horizon_start %s --lead_months %d",
                         "--set_peak_season_year %d --output_dir %s --oni_cache_dir generated_data/"),
                   "2014-09-01", 150L, 2015L, "generated_data/lenssen_events/elnino_2015")
expect(identical(lenssenBuildCommand(ev$launch, 150L, 2015L, "generated_data/lenssen_events/elnino_2015"), old_cmd),
       "lenssenBuildCommand(legacy) == the former event command")
expect(grepl("--match library", lenssenBuildCommand(ev$launch, 150L, 2015L, "x", cf = list(match = "library"))),
       "library mode adds the matcher flags")
expect(!is.null(errmsg(ensoCounterfactualCfg(list(match = "bogus")))), "unknown match refused")
unlink(tmp, recursive = TRUE)

# --- [4] analog seam in the CFSv2 sampler, counterfactual arguments -------
cat("# [4] 8.4 analog_seam, lenssenCounterfactualArgs\n")
T_oni <- as.Date("2026-08-01")
mm <- function(a, b) seq.Date(a, b, by = "1 month")
odates <- mm(as.Date("2020-01-01"), T_oni)
ovals <- 1.5 * sin(seq_along(odates) / 6)
hist_oni <- data.frame(year = as.integer(format(odates, "%Y")), month = as.integer(format(odates, "%m")),
                       enso_index = ovals)
hist_n34 <- hist_oni; hist_n34$enso_index <- hist_n34$enso_index + 0.05
cfs_t <- mm(T_oni %m-% months(2L), T_oni %m+% months(8L))
cfs <- expand.grid(target_date = cfs_t, member = 1:5)
cfs$init_date <- T_oni %m-% months(2L)
cfs$lead <- as.integer(lubridate::interval(cfs$init_date, cfs$target_date) %/% months(1))
cfs$nino34_anomaly <- 1.8 + 0.1 * (cfs$member - 3) - 0.05 * cfs$lead
simm <- mm(T_oni %m-% months(12L), T_oni %m+% months(20L))
samp <- function(...) sampleEnsoCFSv2Trajectories(cfsv2_data = cfs, n_rep = 12L, historical_oni = hist_oni,
  historical_nino34 = hist_n34, sim_months = simm, max_lag = 3L, rebase = "seam", seed = 3L,
  beyond_cfsv2 = "analog", analog_library = lib, analog_k = 10L, ...)
a_def <- suppressMessages(samp())
a_exp <- suppressMessages(samp(analog_seam = "last_cfsv2"))
expect(identical(a_def, a_exp), "analog_seam default == explicit last_cfsv2 (bit-identical)")
a_full <- suppressMessages(samp(analog_seam = "last_full_oni"))
last_cfs <- max(cfs_t)
pre <- simm < last_cfs
expect(identical(a_full$enso_lag0[pre, ], a_def$enso_lag0[pre, ]),
       "last_full_oni: every month before the last CFSv2 month unchanged")
at_last <- simm == last_cfs
expect(!identical(a_full$enso_lag0[at_last, ], a_def$enso_lag0[at_last, ]) &&
         all(is.finite(a_full$enso_lag0[at_last, ])),
       "last_full_oni: the last CFSv2 month (a 2-month ONI mean before) comes from the continuation")
e <- errmsg(suppressMessages(sampleEnsoCFSv2Trajectories(cfsv2_data = cfs, n_rep = 4L, historical_oni = hist_oni,
  historical_nino34 = hist_n34, sim_months = simm, max_lag = 3L, seed = 3L, beyond_cfsv2 = "neutral",
  analog_seam = "last_full_oni")))
expect(!is.null(e) && grepl("needs beyond_cfsv2 = 'analog'", e, fixed = TRUE), "a non-analog continuation refuses last_full_oni")
expect(is.null(ensoContinuationProvenance(list(method = "neutral"))$seam) &&
         identical(ensoContinuationProvenance(list(method = "neutral", seam = "last_full_oni"))$seam, "last_full_oni"),
       "ensoContinuationProvenance records the seam only off the default")

la <- lenssenCounterfactualArgs(NULL, as.Date("2026-05-01"))
expect(identical(la, list(cache_path = "generated_data/lenssen_nino34_analogs.rds",
                          horizon_start = as.Date("2026-05-01"), mean_offset = 0, log_var_scaling = 0,
                          analog_set = "all")), "lenssenCounterfactualArgs defaults")
lb <- lenssenCounterfactualArgs(list(cache_path = "x.rds", lenssen_horizon_start = "2026-04-01",
                                     analog_set = "neutral", mean_offset = 1L), as.Date("2026-05-01"))
expect(identical(lb$horizon_start, as.Date("2026-04-01")) && identical(lb$cache_path, "x.rds") &&
         identical(lb$analog_set, "neutral") && identical(lb$mean_offset, 1),
       "lenssenCounterfactualArgs reads the block")
expect(identical(lenssenCounterfactualArgs(list(lenssen_horizon_start = ""), as.Date("2026-05-01"))$horizon_start,
                 as.Date("2026-05-01")), "empty anchor = the arm's horizon.start")

# --- [5] in-process cache build / verify (generators) -------------------
cat("# [5] 8.4b ensureLenssenCounterfactualCache\n")
gd <- tempfile("gd_"); dir.create(gd)
saveRDS(data.frame(year = as.integer(format(oni$date, "%Y")), month = as.integer(format(oni$date, "%m")),
                   enso_index = oni$enso_index, date = oni$date, enso_phase = ensoPhase(oni$enso_index, "episode")),
        file.path(gd, "enso_oni_cache.rds"))
lib_path <- file.path(gd, "lib.qs")
lib_save <- lib; lib_save$by_month <- NULL; lib_save$path <- NULL
qs2::qs_save(lib_save, lib_path)
cont <- list(library = lib_path, exclude_models = "C")
cf_lib <- list(match = "library", n_analogs = 15L, min_sep_months = 12L, winters = "all_read", min_set_size = 3L)
cpath <- file.path(gd, "lenssen_prod_run", "lenssen_nino34_analogs.rds")
expect(is.null(ensureLenssenCounterfactualCache(cpath, launch, 30L, NULL, cont, "2028-06-01", oni_dir = gd)),
       "legacy (no enso_counterfactual block): nothing built, NULL")
expect(!file.exists(cpath), "legacy: no file written")
r1 <- suppressMessages(ensureLenssenCounterfactualCache(cpath, launch, 30L, cf_lib, cont, "2028-06-01",
                                                        sets = c("all", "neutral"), primary_set = "all", oni_dir = gd))
expect(file.exists(cpath) && isTRUE(r1$status$ok), "library: built in-process and verified")
expect(identical(r1$winters$primary, 2026L) && identical(as.integer(r1$winters$other), 2027L),
       "winters: primary DJF 2026/27, constrained DJF 2027/28 (read end 2028-06)")
expect(!"C" %in% readRDS(cpath)$model, "continuation exclusion applied to the counterfactual pool (D7: both uses)")
mt <- file.mtime(cpath)
Sys.sleep(1.1)
r2 <- suppressMessages(ensureLenssenCounterfactualCache(cpath, launch, 30L, cf_lib, cont, "2028-06-01",
                                                        sets = "all", oni_dir = gd))
expect(identical(file.mtime(cpath), mt) && identical(r2$provenance$digest, r1$provenance$digest),
       "second call verifies the cache (no rebuild, same digest)")
e <- errmsg(suppressMessages(ensureLenssenCounterfactualCache("generated_data/lenssen_nino34_analogs.rds", launch, 30L,
                                                              cf_lib, cont, "2028-06-01", oni_dir = gd)))
expect(!is.null(e) && grepl("shared legacy cache", e), "library mode refuses the shared legacy cache path")
e <- errmsg(suppressMessages(ensureLenssenCounterfactualCache(cpath, launch, 30L,
  modifyList(cf_lib, list(min_set_size = 10000L)), cont, "2028-06-01", sets = c("all", "neutral"),
  primary_set = "neutral", oni_dir = gd)))
expect(!is.null(e) && grepl("under enso_counterfactual.min_set_size", e), "primary set under min_set_size: refused")
expect(identical(r1$spec$seam_at, "launch") &&
         isTRUE(all.equal(readRDS(cpath)$oni_trajectories[, 1],
                          rep(oni_lkp[["2026-05-01"]], nrow(readRDS(cpath)$oni_trajectories)))),
       "generators: the in-process build pins column 1 by default (seam_at launch)")
unlink(gd, recursive = TRUE)

# T5: the deploy defaults restate the in-function ones.
cat("\n[T5] get_default_options()$deploy = ENSO_COUNTERFACTUAL_DEFAULTS / seam last_cfsv2\n")
source("analysis/utils/config.R")
dflt <- get_default_options()$deploy
expect(identical(dflt$enso_counterfactual, ENSO_COUNTERFACTUAL_DEFAULTS),
       "deploy.enso_counterfactual default = ENSO_COUNTERFACTUAL_DEFAULTS")
expect(identical(ensoCounterfactualCfg(dflt$enso_counterfactual), ensoCounterfactualCfg(NULL)),
       "the default block resolves as a missing block (legacy)")
expect(identical(dflt$enso_continuation$seam, "last_cfsv2"), "deploy.enso_continuation.seam default last_cfsv2")

# --- [6] launch-month seam pin (seam_at) ---------------------------------------
cat("# [6] seam_at launch / launch_minus_1\n")
old_spec <- function(launch, n_months, cf, lib, seam_state,
                     peak_months = c(12L, 1L, 2L), thresholds = ENSO_PHASE_THRESHOLDS,
                     primary_years = NULL, other_years = integer(0), read_end = NULL,
                     exclude_models = character(0)) {
  cf <- ensoCounterfactualCfg(cf)
  launch <- as.Date(launch)
  list(match = "library", horizon_start = format(launch),
       seam_month = format(seq(launch, by = "-1 month", length.out = 2L)[2L]),
       seam_state = lapply(seam_state[c("v", "tend", "phase")], function(v) round(as.numeric(v), 6)),
       n_months = as.integer(n_months), n_analogs = cf$n_analogs, min_sep_months = cf$min_sep_months,
       blend_months = cf$blend_months,
       exclude_models = as.list(sort(as.character(unlist(exclude_models)))),
       library = lib$path %||% "", library_digest = lib$meta$digest %||% ensoContinuationLibraryDigest(lib),
       sd_window = lib$meta$rescale$window %||% "", obs_sd = lib$meta$rescale$obs_sd %||% NA_real_,
       peak_months = as.integer(peak_months), thresholds = as.list(thresholds),
       winters = cf$winters, read_span = cf$read_span, other_winters = cf$other_winters,
       read_end = if (is.null(read_end)) "" else format(as.Date(read_end)),
       primary_years = as.list(as.integer(primary_years %||% integer(0))),
       other_years = as.list(as.integer(other_years)))
}
expect(identical(ensoCounterfactualCfg(NULL)$seam_at, "launch"), "default seam_at is launch (library mode)")
expect(!is.null(errmsg(ensoCounterfactualCfg(list(seam_at = "peak")))), "unknown seam_at refused")
cf_lm1 <- list(match = "library", n_analogs = 15L, min_sep_months = 12L, winters = "all_read",
               seam_at = "launch_minus_1")
cf_l <- modifyList(cf_lm1, list(seam_at = "launch"))
a_old <- old_spec(launch, 30L, cf_lm1, lib, st, primary_years = 2026L, other_years = 2027L,
                  read_end = "2028-06-01", exclude_models = "C")
a_lm1 <- lenssenCacheSpec(launch, 30L, cf_lm1, lib, st, primary_years = 2026L, other_years = 2027L,
                          read_end = "2028-06-01", exclude_models = "C")
expect(identical(a_old, a_lm1), "launch_minus_1: the spec is the pre-pin spec, field for field (same digests)")
obs_l <- ensoLaunchOni(oni_lkp, launch)
expect(is.null(ensoLaunchOni(oni_lkp, launch, "launch_minus_1")) && identical(obs_l, oni_lkp[["2026-05-01"]]),
       "ensoLaunchOni: the published launch ONI (NULL under launch_minus_1)")
expect(!is.null(errmsg(ensoLaunchOni(oni_lkp, "2027-01-01"))), "ensoLaunchOni: an unpublished launch month stops")
expect(!is.null(errmsg(lenssenCacheSpec(launch, 30L, cf_l, lib, st, primary_years = 2026L))),
       "seam_at launch without launch_oni: refused")
a_l <- lenssenCacheSpec(launch, 30L, cf_l, lib, st, primary_years = 2026L, other_years = 2027L,
                        read_end = "2028-06-01", exclude_models = "C", launch_oni = obs_l)
expect(identical(a_l[names(a_lm1)], a_lm1) && identical(a_l$seam_at, "launch") &&
         identical(a_l$launch_oni, round(obs_l, 6)), "launch: the pre-pin spec + seam_at + launch_oni")
c_lm1 <- buildLenssenCache(a_lm1, oni, lib)
c_l <- buildLenssenCache(a_l, oni, lib)
pool_lm1 <- buildLenssenCounterfactualPool(lib, st, 30L, 15L, 12L, "C", 3)
expect(identical(c_lm1$oni_trajectories, pool_lm1$traj) && is.null(c_lm1$meta$seam_at),
       "launch_minus_1 cache = the unpinned pool, no seam_at in the meta")
expect(identical(c_l$model, c_lm1$model) && identical(c_l$init_index, c_lm1$init_index) &&
         identical(c_l$distance, c_lm1$distance), "launch: same analogs (selection on the launch - 1 state)")
expect(all(c_l$oni_trajectories[, 1] == obs_l), "launch: column 1 = the observed launch-month ONI, every row")
r <- 3L; s_r <- lib$series[[c_l$model[r]]]; i_r <- c_l$init_index[r]
expect(isTRUE(all.equal(c_l$oni_trajectories[r, 2:30],
                        s_r[i_r + 2:30] + (obs_l - s_r[i_r + 1L]) * exp(-(1:29) / 3))),
       "launch: column j = analog(i + j) + (obs - analog(i + 1)) rho^(j - 1)")
expect(all(abs(vapply(c_l$meta$sets, `[[`, numeric(1), "seam_error")) < 1e-9),
       "launch: seam error 0 in every set (S35 / D29)")
expect(any(abs(vapply(c_lm1$meta$sets, `[[`, numeric(1), "seam_error")) > 1e-3),
       "launch_minus_1: a non-zero seam error (what the pin removes)")
expect(!identical(c_l$meta$digest, c_lm1$meta$digest) && identical(c_l$meta$seam_at, "launch"),
       "the pin enters the content digest and the meta")
tp <- tempfile("pin_"); dir.create(tp); pp <- file.path(tp, "lenssen_nino34_analogs.rds")
invisible(buildLenssenCache(a_l, oni, lib, out_path = pp))
expect(isTRUE(lenssenCacheStatus(pp, a_l)$ok) && grepl("seam_at", lenssenCacheStatus(pp, a_lm1)$reason),
       "a pinned cache does not satisfy a launch_minus_1 spec (and vice versa by field)")
expect(identical(lenssenCacheProvenance(pp)$seam_at, "launch"), "provenance records seam_at launch")
pp2 <- file.path(tp, "lm1", "lenssen_nino34_analogs.rds")
invisible(buildLenssenCache(a_lm1, oni, lib, out_path = pp2))
expect(!"seam_at" %in% names(lenssenCacheProvenance(pp2)), "launch_minus_1 provenance: no seam_at element")
pin1 <- buildLenssenCounterfactualPool(lib, st, 1L, 15L, 12L, "C", 3, pin = 0.7)
expect(identical(dim(pin1$traj), c(nrow(pool_lm1$traj), 1L)) && all(pin1$traj == 0.7), "one-month pinned pool")
expect(grepl("--seam_at launch", lenssenBuildCommand(launch, 30L, 2026L, "x", cf = list(match = "library"))) &&
         grepl("--seam_at launch_minus_1", lenssenBuildCommand(launch, 30L, 2026L, "x", cf = cf_lm1)),
       "the library build command carries --seam_at")
unlink(tp, recursive = TRUE)

cat(sprintf("\n[test_enso_counterfactual_match] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
