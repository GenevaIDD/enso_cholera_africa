# test_lag_averaging.R -- Pareto-k of the per-lag PSIS-LOO
# and the per-(batch, group) lag sidecar of 05_fit_postprocess_loo_lags.R.
#
# HARD gate (see analysis/tests/run_all.R). No fit file, no data, no MCMC:
# loo objects come from loo::loo() on simulated log-likelihood matrices, the
# "fits" are a data.frame + a posterior::draws_matrix.
#   [1] pareto_k_summary(): counts / threshold / max / share on synthetic loo
#       lists and on a real psis_loo object.
#   [2] ranking_pareto(): k columns present (pass-through) / computed from
#       per-fit loo caches (paths or in-memory objects; partial) / absent (NA).
#   [3] every selection rule picks the same lag with or without the k columns.
#   [4] sidecar: entries -> build_lag_loo_sidecar() (sorted by lag, N x L
#       matrices, draws, cb_meta, digests) -> write under the
#       makeLooCompareRunFile() sidecar kind -> read_lag_loo_sidecar().
#   [5] row alignment: a row-key digest or N mismatch across lags stops (on
#       build and on read); the crossbasis columns are not part of the key.
# Run: Rscript analysis/tests/test_lag_averaging.R

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/loo_lags.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
stops <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

cat("# lag averaging, part (a): Pareto-k + lag sidecar\n")

# Simulated pointwise log-likelihood (S draws x N rows) -> a real psis_loo.
sim_loo <- function(S = 200L, N = 30L, seed = 1L, shift = 0) {
  set.seed(seed)
  mu <- rnorm(N, shift, 0.3)
  ll <- sapply(mu, function(m) dnorm(rnorm(S, m, 0.2), 0, 1, log = TRUE))
  suppressWarnings(loo::loo(ll, r_eff = rep(1, N)))
}

# [1] pareto_k_summary ---------------------------------------------------------
cat("\n[1] pareto_k_summary\n")
syn <- list(diagnostics = list(pareto_k = c(0.1, 0.5, 0.71, 0.9, 1.2, NA)))
ps  <- pareto_k_summary(syn)
expect(identical(names(ps), PARETO_K_COLS), "columns = PARETO_K_COLS, in order")
expect(ps$n_obs == 6L && ps$n_k_gt07 == 3L && ps$n_k_gt1 == 1L && ps$n_k_bad == 3L,
       "counts: n_obs 6 (NA counted), k > 0.7: 3, k > 1: 1")
expect(isTRUE(all.equal(ps$k_max, 1.2)) && isTRUE(all.equal(ps$share_k_bad, 0.5)),
       "k_max 1.2, share_k_bad 0.5")
expect(identical(ps$k_threshold, if ("pareto_k_threshold" %in% getNamespaceExports("loo")) ps$k_threshold else 0.7),
       "threshold: loo's exported pareto_k_threshold() if any, else 0.7")
expect(pareto_k_summary(list(pareto_k = c(0.2, 0.3)))$n_k_bad == 0L, "bare $pareto_k accepted; all good -> 0 bad")
expect(is.na(pareto_k_summary(list(pareto_k = c(NA_real_, NA_real_)))$k_max), "all-NA k -> k_max NA")
expect(stops(pareto_k_summary(list(estimates = 1))), "no Pareto-k in the object -> stop")
lo <- sim_loo()
pr <- pareto_k_summary(lo)
expect(pr$n_obs == 30L && isTRUE(all.equal(pr$k_max, max(loo::pareto_k_values(lo)))),
       "real psis_loo: n_obs and k_max agree with loo::pareto_k_values()")

# [2] ranking_pareto -------------------------------------------------------------
cat("\n[2] ranking_pareto\n")
mk_ranking <- function(lags, elpd, se_diff = 1) {
  best <- max(elpd)
  out <- data.frame(single_lag = as.integer(lags), label = paste0("lag", lags),
                    elpd_loo = elpd, se_elpd = 5, elpd_diff = elpd - best,
                    se_diff = ifelse(elpd == best, 0, se_diff), p_loo = 10)
  out$rank <- rank(-out$elpd_loo, ties.method = "first")
  out <- out[order(out$rank), c("rank", setdiff(names(out), "rank"))]
  rownames(out) <- NULL
  out
}
rk <- mk_ranking(c(-2, 0, 3, 7), c(-100, -101, -102, -130))
tmp <- tempfile("lag_avg_"); dir.create(tmp)
loos  <- list(`-2` = sim_loo(seed = 1), `0` = sim_loo(seed = 2), `3` = sim_loo(seed = 3), `7` = sim_loo(seed = 4))
files <- vapply(names(loos), function(l) { f <- file.path(tmp, paste0("loo_", l, ".rds")); saveRDS(loos[[l]], f); f },
                character(1))

r_none <- ranking_pareto(rk)
expect(identical(attr(r_none, "pareto_source"), "none") && all(PARETO_K_COLS %in% names(r_none)) &&
         all(is.na(r_none$k_max)), "no loo_files -> NA columns, source 'none'")
expect(identical(r_none[, names(rk)], rk), "original columns and rows unchanged")

r_file <- ranking_pareto(rk, loo_files = files)
expect(identical(attr(r_file, "pareto_source"), "loo_cache"), "per-fit loo cache paths -> source 'loo_cache'")
expect(isTRUE(all.equal(r_file$k_max[r_file$single_lag == 3],
                        max(loo::pareto_k_values(loos[["3"]])))), "k_max matched to its lag (not to the row order)")
r_mem <- ranking_pareto(rk, loo_files = loos)
expect(isTRUE(all.equal(r_mem, r_file, check.attributes = FALSE)), "in-memory loo objects == the same loo read from files")

r_part <- ranking_pareto(rk, loo_files = files[c("-2", "0")])
expect(identical(attr(r_part, "pareto_source"), "loo_cache_partial") &&
         sum(is.na(r_part$n_obs)) == 2L, "missing lags -> NA rows, source 'loo_cache_partial'")
r_bad <- ranking_pareto(rk, loo_files = c(`-2` = file.path(tmp, "nope.rds")))
expect(identical(attr(r_bad, "pareto_source"), "none"), "no readable file -> source 'none'")

r_pass <- ranking_pareto(r_file)
expect(identical(attr(r_pass, "pareto_source"), "ranking") &&
         isTRUE(all.equal(r_pass, r_file, check.attributes = FALSE)),
       "k columns already present -> passed through unchanged, source 'ranking'")
r_half <- r_file; r_half$k_max <- NULL
expect(identical(attr(ranking_pareto(r_half, loo_files = loos), "pareto_source"), "loo_cache") &&
         identical(names(ranking_pareto(r_half, loo_files = loos)), names(r_file)),
       "a partial k-column set is recomputed (column order restored)")

# [3] selection is unchanged by the k columns ------------------------------------
cat("\n[3] selection identical with / without the Pareto-k columns\n")
rks <- list(rk,
            mk_ranking(1:4, c(-100, -101, -102, -130)),
            mk_ranking(-2:3, c(-104, -100, -101, -101.5, -102, -130), se_diff = 1.5),
            mk_ranking(c(-6, 7, 9, -10, 2, 11), c(0, -3, -5, -7, -9, -12), se_diff = 5))
same <- TRUE
for (r0 in rks) {
  r1 <- ranking_pareto(r0, loo_files = setNames(rep(list(sim_loo()), nrow(r0)), as.character(r0$single_lag)))
  picks <- function(r) c(
    select_best_lag(r, "occurrence")$lag,
    select_best_lag(r, "occurrence", min_lag = -Inf)$lag,
    select_best_lag(r, "duration", occ_lag = 2L)$lag,
    select_best_lag(r, "monthly", enforce_occ_floor = FALSE)$lag,
    select_best_lag(r, "occurrence", rule = "best_elpd", min_lag = -Inf)$lag,
    select_best_lag(r, "duration", occ_lag = 3L, rule = "best_elpd", min_lag = -Inf, occ_floor_offset = 4)$lag,
    select_best_lag_within_eligible(r, occ_lag = 0L)$lag)
  same <- same && identical(picks(r0), picks(r1))
}
expect(same, "every rule (longest_within_se / best_elpd / within_eligible, all tracks) picks the same lag")

# [4] sidecar: build, write, read ------------------------------------------------
cat("\n[4] lag sidecar\n")
N <- 30L; S <- 200L; cbn <- paste0("enso_cb", 1:3)
fake_fit <- function(lag, seed) {
  set.seed(seed)
  d <- data.frame(y = rep(0:1, length.out = N), country = factor(rep(c("A", "B", "C"), length.out = N)),
                  log_offset = log(seq_len(N)))
  for (cc in cbn) d[[cc]] <- rnorm(N)                       # the crossbasis differs by lag
  dm <- posterior::as_draws_matrix(matrix(rnorm(S * 4), S, 4,
                                          dimnames = list(NULL, c("b_Intercept", paste0("b_", cbn)))))
  list(data = d, draws = dm)
}
cb_meta_for <- function(lag) list(single_lag = as.integer(lag), cb_col_names = cbn, range = c(-2, 2.5),
                                  onebasis_attrs = list(fun = "ns", knots = c(-0.5, 0.5), Boundary.knots = c(-2, 2.5),
                                                        intercept = FALSE, range = c(-2, 2.5)))
lags <- c(2L, -1L, 0L)                                       # deliberately unsorted
entries <- lapply(seq_along(lags), function(i) {
  ff <- fake_fit(lags[i], 10 + i)
  lag_sidecar_entry(lags[i], if (lags[i] < 0) paste0("lagm", -lags[i]) else sprintf("lag%02d", lags[i]),
                    loo = sim_loo(S = S, N = N, seed = 20 + i), fit_data = ff$data,
                    b_draws = fit_cb_pop_draws(ff$draws, cbn), cb_meta_enso = cb_meta_for(lags[i]),
                    files = c(fit = "f.rds", loo = "l.rds", cb_meta = "m.rds"))
})
expect(identical(dim(entries[[1]]$b_draws), c(S, 3L)) && identical(colnames(entries[[1]]$b_draws), cbn),
       "fit_cb_pop_draws(): S x K matrix with the cb column names (b_Intercept dropped)")
sc <- build_lag_loo_sidecar(entries, track = "duration", group = "enso_west", batch = "analysis/configs/fake_dur_batch/")
expect(identical(unname(sc$lags), c(-1L, 0L, 2L)) && identical(colnames(sc$elpd_pointwise), c("lagm1", "lag00", "lag02")),
       "lags sorted, matrix columns labelled as 05 labels them")
expect(identical(dim(sc$elpd_pointwise), c(N, 3L)) && identical(dim(sc$pareto_k), c(N, 3L)) && sc$n_obs == N,
       "N x L pointwise elpd and Pareto-k matrices")
expect(isTRUE(all.equal(unname(colSums(sc$elpd_pointwise)), sc$estimates$elpd_loo)),
       "column sums of the pointwise elpd = each lag's elpd_loo")
expect(identical(sc$track, "duration") && identical(sc$group, "enso_west") && identical(sc$batch, "fake_dur_batch"),
       "track, group and batch tag stored explicitly")
expect(identical(names(sc$cb), colnames(sc$elpd_pointwise)) && identical(sc$cb$lag02$cb_meta$single_lag, 2L) &&
         identical(dim(sc$cb$lagm1$b_draws), c(S, 3L)), "per-lag draws + cb_meta keyed by lag label")
expect(length(unique(sc$row_digest)) == 1L, "one row digest although the crossbasis columns differ by lag")
expect(identical(names(sc$pareto), c("single_lag", PARETO_K_COLS)) && nrow(sc$pareto) == 3L,
       "per-lag Pareto-k summary table")

opt <- list(outdir = tmp)
f_sc <- makeLooCompareRunFile(opt, "loo_compare_lag_sidecar", batch = "fake_dur_batch", group = "enso_west",
                              date = as.Date("2026-10-27"))
dir.create(dirname(f_sc), recursive = TRUE, showWarnings = FALSE)
saveRDS(sc, f_sc)
expect(grepl(LAG_LOO_SIDECAR_RE, basename(f_sc)) && !grepl(LOO_COMPARE_RANKING_RE, basename(f_sc)),
       "sidecar filename matches LAG_LOO_SIDECAR_RE, not the ranking regex")
rd <- read_lag_loo_sidecar("fake_dur_batch", "enso_west", opt = opt)
expect(identical(attr(rd, "file"), f_sc) && isTRUE(all.equal(rd$elpd_pointwise, sc$elpd_pointwise)),
       "read_lag_loo_sidecar() finds and returns it")
expect(nrow(loo_compare_ranking_files("fake_dur_batch", opt = opt, latest_only = FALSE)) == 0L,
       "the ranking locator does not see the sidecar")
expect(is.null(read_lag_loo_sidecar("fake_dur_batch", "enso_east", opt = opt, required = FALSE)) &&
         stops(read_lag_loo_sidecar("fake_dur_batch", "enso_east", opt = opt)),
       "no sidecar for the group -> NULL / stop when required")
f_old <- makeLooCompareRunFile(opt, "loo_compare_lag_sidecar", batch = "fake_dur_batch", group = "enso_west",
                               date = as.Date("2026-01-01"))
saveRDS(list(schema = "old"), f_old); Sys.setFileTime(f_old, Sys.time() + 3600)
expect(identical(latest_lag_sidecar_file("fake_dur_batch", "enso_west", opt = opt), f_sc),
       "latest sidecar by filename date, not mtime")

# [5] row alignment ---------------------------------------------------------------
cat("\n[5] row alignment across lags\n")
expect(identical(fit_row_digest(data.frame(a = factor(c("x", "y")), b = 1:2)),
                 fit_row_digest(data.frame(a = c("x", "y"), b = 1:2))),
       "row digest: factor and character keys agree")
d0 <- fake_fit(0, 1)$data
expect(!identical(as.character(fit_row_digest(d0, cbn)), as.character(fit_row_digest(d0))),
       "the crossbasis columns are left out of the key only when excluded")
bad <- entries
ff  <- fake_fit(0L, 99); ff$data$country[3] <- "A"; ff$data$country[1] <- "B"
bad[[3]] <- lag_sidecar_entry(0L, "lag00", loo = sim_loo(S = S, N = N, seed = 23), fit_data = ff$data,
                              b_draws = fit_cb_pop_draws(ff$draws, cbn), cb_meta_enso = cb_meta_for(0L))
expect(stops(build_lag_loo_sidecar(bad, track = "duration", group = "enso_west", batch = "x")),
       "a row-key digest mismatch across lags stops the build")
short <- entries
ff2 <- fake_fit(2L, 5); ff2$data <- ff2$data[-1, ]
short[[1]] <- lag_sidecar_entry(2L, "lag02", loo = sim_loo(S = S, N = N - 1L, seed = 21), fit_data = ff2$data,
                                b_draws = fit_cb_pop_draws(ff2$draws, cbn), cb_meta_enso = cb_meta_for(2L))
expect(stops(build_lag_loo_sidecar(short, track = "duration", group = "enso_west", batch = "x")),
       "a different N in one lag stops the build")
tampered <- sc; tampered$row_digest[2] <- "0000000000000000"
saveRDS(tampered, f_sc)
expect(stops(read_lag_loo_sidecar("fake_dur_batch", "enso_west", opt = opt)),
       "a stored digest mismatch stops the reader")
expect(stops(lag_sidecar_entry(0L, "lag00", loo = list(estimates = NULL), fit_data = d0,
                               b_draws = entries[[1]]$b_draws)),
       "a loo object without pointwise elpd stops the entry")

# ===========================================================================
# Parts (b) and (c): report_lags v2 rows, the lag floor, and the lag_averaging
# step (results_derive.R) on synthetic sidecars.
# ===========================================================================
source("analysis/utils/results_core.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/results_derive.R")
cat("\n# lag averaging, parts (b) / (c)\n")

# [6] the deploy floor ------------------------------------------------------------
cat("\n[6] lag floor\n")
r_occ4 <- lag_selection_cfg(list(rule = "best_elpd", se_mult = 2, min_lag = -8, occ_floor_offset = 4))
expect(lag_floor_of(r_occ4, "occurrence", NULL) == -8 && lag_floor_of(r_occ4, "duration", -6L) == -8 &&
         lag_floor_of(r_occ4, "monthly", 0L) == -4,
       "min_lag -8, occ_floor_offset 4: occurrence -8; downstream max(min_lag, occ - 4) (West -8, East -4)")
expect(identical(lag_floor_of(lag_selection_cfg(NULL), "occurrence"), 0) &&
         lag_floor_of(lag_selection_cfg(NULL), "duration", 3L, enforce_occ_floor = TRUE) == 3 &&
         lag_floor_of(lag_selection_cfg(NULL), "duration", 3L) == 0,
       "longest_within_se: min_lag 0; occ lag floor only under enforce_occ_floor")
expect(lag_floor_of(lag_selection_cfg(list(rule = "best_elpd")), "occurrence") == -Inf,
       "best_elpd without min_lag: unbounded (-Inf)")
rkd <- mk_ranking(-11:11, -abs((-11:11) + 9) * 3 - 100)
for (ol in c(-6L, 0L, 4L)) {
  fl <- lag_floor_of(r_occ4, "duration", ol)
  el <- select_best_lag_best_elpd(rkd, min_lag = -8, occ_lag = ol, occ_floor_offset = 4)$within_se$single_lag
  expect(setequal(rkd$single_lag[rkd$single_lag >= fl], el),
         sprintf("floor reproduces select_best_lag_best_elpd's eligible set (occ lag %d)", ol))
}

# [7] report_lags v2 selected row ---------------------------------------------------
cat("\n[7] report_lags v2 rows\n")
rkw <- mk_ranking(-11:11, c(-120, -110, -100, -101, -103, -104, -105, -106, -107, -108, -109, -111,
                            -112, -113, -114, -115, -116, -117, -118, -119, -121, -122, -123), se_diff = 2)
rkw$selected <- rkw$single_lag == -8L
row <- report_lags_selected_row(rkw, "duration", "enso_west", -8L, -8)
expect(row$rank1_lag == -9L && row$floor_lag == -8L && isTRUE(row$bound_active) &&
         isTRUE(all.equal(row$bound_cost, -1)) && isTRUE(all.equal(row$bound_cost_se, 2)),
       "rank-1 -9 below the floor -8: bound active, cost = the deployed lag's paired elpd diff")
expect(row$rank2_lag == -8L && isTRUE(all.equal(row$elpd_diff_rank2, -1)) &&
         row$best_nonneg_lag == 0L && isTRUE(all.equal(row$elpd_diff_nonneg, -11)),
       "rank-2 lag and the best lag >= 0 with their elpd difference to rank 1")
expect(is.na(row$n_k_gt07_deployed) && all(c(names(row)[1:8]) == c("track", "enso_group", "deployed_lag",
         "selected_lag", "rank1_lag", "elpd_diff_deployed", "se_diff_deployed", "n_lags")),
       "v1 columns first, unchanged; no k columns in the ranking -> NA")
rkk <- ranking_pareto(rkw, loo_files = setNames(rep(list(sim_loo()), nrow(rkw)), as.character(rkw$single_lag)))
rkk$selected <- rkw$selected[match(rkk$single_lag, rkw$single_lag)]
row_k <- report_lags_selected_row(rkk, "duration", "enso_west", -8L, -8)
expect(identical(row_k$n_k_gt07_deployed, as.integer(rkk$n_k_gt07[rkk$single_lag == -8])),
       "the deployed lag's Pareto-k counts when the ranking has them")
row_u <- report_lags_selected_row(rkw, "occurrence", "enso_west", -8L, -Inf)
expect(is.na(row_u$floor_lag) && !row_u$bound_active && is.na(row_u$bound_cost),
       "unbounded rule: no floor, never active")

# [8] weights, differences, k-robust --------------------------------------------------
cat("\n[8] weights and paired differences\n")
set.seed(5)
Nw <- 80L
base <- rnorm(Nw, -2, 0.5)
Ew <- cbind(a = base, b = base - 0.5 + rnorm(Nw, 0, 0.05), c = base - 1)
w_st <- lag_weights(Ew, 1:3, "stacking")
expect(w_st[1] > 0.95 && isTRUE(all.equal(sum(w_st), 1)), "a dominant lag gets stacking weight ~1")
Eeq <- cbind(base, base, base)
w_pb <- lag_weights(Eeq, 1:3, "pseudobma_plus")
expect(isTRUE(all.equal(w_pb, rep(1 / 3, 3))), "identical pointwise elpd -> equal pseudo-BMA+ weights")
w_el <- lag_weights(Ew, 2:3, "stacking")
expect(w_el[1] == 0 && isTRUE(all.equal(sum(w_el), 1)), "a lag outside the set gets weight 0")
expect(identical(lag_weights(Ew, 2L, "pseudobma_plus"), c(0, 1, 0)), "a one-lag set gets weight 1")
expect(identical(lag_weights(Ew, 1:3, "pseudobma_plus", seed = 3), lag_weights(Ew, 1:3, "pseudobma_plus", seed = 3)),
       "pseudo-BMA+ (Bayesian bootstrap) is seeded")
lo3 <- list(sim_loo(seed = 31, N = 40), sim_loo(seed = 32, N = 40, shift = 0.2), sim_loo(seed = 33, N = 40))
E3 <- do.call(cbind, lapply(lo3, function(l) l$pointwise[, "elpd_loo"]))
d3 <- lag_elpd_diff(E3, 1L)
cmp <- loo::loo_compare(lo3[[1]], lo3[[2]])
expect(isTRUE(all.equal(abs(d3$diff[2]), abs(cmp[2, "elpd_diff"]))) &&
         isTRUE(all.equal(d3$se[2], cmp[2, "se_diff"])),
       "paired elpd difference and SE = loo::loo_compare()")
K3 <- matrix(0.1, 40, 3); K3[1:3, 2] <- 0.9; K3[5, 3] <- 1.2
kr <- lag_krobust(E3, K3, 0.7)
expect(kr$n_dropped == 4L && isTRUE(all.equal(kr$diff[kr$rank1], 0)),
       "k-robust: the union of rows with k > 0.7 at any lag is dropped; rank 1 re-found")

# [9] ERFs: grid, centring, mixture ----------------------------------------------------
cat("\n[9] ERF grid and mixture\n")
gr <- lag_common_grid(list(c(-2, 2.5), c(-1.5, 3), c(-2.2, 2)))
expect(min(gr) == -1.5 && max(gr) == 2 && length(gr) == LAG_AVG_GRID_N && 0 %in% gr,
       "common grid = intersection of the ranges, 200 points incl. 0")
vm <- cb_meta_for(0L); b <- sc$cb$lag00$b_draws
lr <- lag_erf_logrr(vm, b, c(-1, 0, 1))
expect(identical(dim(lr), c(3L, nrow(b))) && all(lr[2, ] == 0), "log-IRR draws: x by S, exactly 0 at ONI 0")
m1 <- lag_erf_logrr(vm, b[1:50, ], gr); m2 <- lag_erf_logrr(cb_meta_for(2L), sc$cb$lag02$b_draws[1:50, ], gr)
mx <- lag_mixture_draws(list(m1, m2), c(1, 0), 50L)
expect(all(attr(mx, "lag_index") == 1L) &&
         identical(.lag_erf_summary(mx, gr), .lag_erf_summary(m1, gr)),
       "mixture with weights (1, 0) = the single lag's ERF")
mh <- lag_mixture_draws(list(m1, m2), c(0.3, 0.7), 50L)
expect(sum(attr(mh, "lag_index") == 1L) == 15L, "systematic sampling: 30% of 50 draws from lag 1")

# [10] lag_averaging_group / derive_lag_averaging on the synthetic sidecar -------------
cat("\n[10] lag_averaging step\n")
conds <- enso_conds_named(c(-1, 1.5))
cfg <- lag_avg_cfg(list(enable = TRUE))
la <- lag_averaging_group(sc, 0L, 0, cfg, conds, 60L)
expect(identical(la$lags$single_lag, c(-1L, 0L, 2L)) && identical(la$lags$eligible, c(FALSE, TRUE, TRUE)) &&
         la$lags$w_stacking_eligible[1] == 0 &&
         isTRUE(all.equal(sum(la$lags$w_stacking_all), 1)) && isTRUE(all.equal(sum(la$lags$w_pseudobma_plus_eligible), 1)),
       "per-lag table: eligible = lags >= the floor; weights per method x set, 0 outside the set")
expect(isTRUE(all.equal(la$lags$elpd_diff_deployed[la$lags$deployed], 0)) &&
         la$roles$single_lag[la$roles$role == "deployed"] == 0L,
       "differences to the deployed lag; roles table")
expect(nrow(dplyr::filter(la$erf, variant == "mixture")) == 4L * length(la$grid) &&
         all(abs(dplyr::filter(la$erf, at == 0)$logRR) < 1e-12),
       "ERF rows: 3 lags + 4 mixtures on the grid, logRR(0) = 0")
expect(nrow(la$irr) == (3L + 4L) * 2L && all(la$irr$p_gt1 >= 0 & la$irr$p_gt1 <= 1),
       "IRR at every ENSO condition for every lag and mixture")
expect(stops(lag_averaging_group(sc, 5L, 0, cfg, conds, 60L)), "a deployed lag absent from the sidecar stops")
sc_bad <- sc; sc_bad$row_digest[2] <- "zz"
expect(stops(lag_averaging_group(sc_bad, 0L, 0, cfg, conds, 60L)), "mismatched row digests stop")

tmp2 <- tempfile("la_ctx_"); dir.create(tmp2)
opt2 <- list(outdir = tmp2)
f2 <- makeLooCompareRunFile(opt2, "loo_compare_lag_sidecar", batch = "fake_dur_grouped_enso_singlelag",
                            group = "enso_west", date = as.Date("2026-10-27"))
dir.create(dirname(f2), recursive = TRUE, showWarnings = FALSE); saveRDS(sc, f2)
ctx <- list2env(list(opt = opt2, LAG_LOO_BATCH = list(default = "fake"), LAG_AVG = cfg,
                     MAN = list(groups = list(enso_west = list(best_lag_occ = 0L, best_lag_dur = 0L,
                                                               best_lag_mon = 0L),
                                              enso_east = list(best_lag_occ = 1L, best_lag_dur = 1L))),
                     ENSO_CONDS = conds, N_IRR_DRAWS = 60L, RES = list(deploy_config = ""),
                     ENFORCE_OCC_FLOOR = FALSE, ENSO_BEST_LAG_DIR = file.path(tmp2, "none")),
                parent = emptyenv())
dla <- suppressMessages(derive_lag_averaging(ctx, "duration"))
expect(identical(unique(dla$lags$enso_group), "enso_west") && identical(dla$missing_groups, "enso_east") &&
         all(dla$lags$track == "duration") && dla$meta$floor_lag == 0L,
       "derive_lag_averaging(): groups with a sidecar, missing ones listed, the deploy floor")
msg <- character(0)
res <- withCallingHandlers(ensure_lag_averaging(ctx, "monthly"),
                           message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
expect(is.null(res) && any(grepl("no lag sidecar.*skipped", msg)),
       "no sidecar for a track: the step is skipped with a message, not an error")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
