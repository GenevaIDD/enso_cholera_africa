# test_results_membership.R -- the 09 report_membership step: membership uncertainty of
# the deployed classes from the extended 05
# kl-sweep comparison.
#
# HARD gate (see analysis/tests/run_all.R). Synthetic compare dir + registry,
# no fit, no MCMC. Pins:
#   [1] membership_map_partition(): argmax mean responsibility, p_max
#   [2] the MAP == registry assertion stops, naming the countries; a missing
#       cell stops
#   [3] membership_elpd_diff(): sum and SE = loo::loo_compare() on the same
#       pointwise values; Pareto k > 0.7 counts
#   [4] country / group tables: responsibility row, p_max, carried-over
#       countries NA, soft size = sum of r; sensitivity: ARI 1 and full
#       agreement for the deployed cell, label-free matching for another cell
#   [5] M2: degenerate responsibilities -> the hard ERF; linear in r
#   [6] derive_report_membership() end to end (M2 skipped without sidecars, and
#       computed with them)
#   [7] (soft, local) the deployed ns5_K5 cell's p_max equals gen_grouping's
#       membership CSV when both are on disk
#
# Run: Rscript analysis/tests/test_results_membership.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(tidyr) })
source("analysis/utils/results_core.R")
for (.m in c("loo_lags", "kl_sweep", "dlnm")) source(file.path("analysis/utils", paste0(.m, ".R")))
source("analysis/utils/results_derive.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

# Synthetic sweep: cells A (deployed, 2 classes) and B (2 classes, labels swapped
# names), 5 fitted countries; registry adds a carried-over country (ZZZ).
isos <- c("AAA", "BBB", "CCC", "DDD", "EEE")
rmat <- rbind(AAA = c(0.9, 0.1), BBB = c(0.55, 0.45), CCC = c(0.2, 0.8), DDD = c(0.05, 0.95), EEE = c(1, 0))
mk_resp <- function(id, R, labels, chains) {
  tibble(id = id, chains = chains, iso3 = rep(rownames(R), each = ncol(R)), class = rep(seq_len(ncol(R)), nrow(R)),
         label = rep(labels, nrow(R)), mean = as.vector(t(R)), sd = 0.1, q2.5 = 0, q97.5 = 1,
         prior_group = "g_x", identify = "seed")
}
RB <- rbind(AAA = c(0.2, 0.8), BBB = c(0.3, 0.7), CCC = c(0.9, 0.1), DDD = c(0.6, 0.4), EEE = c(0.1, 0.9))
resp <- bind_rows(mk_resp("A", rmat, c("g_one", "g_two"), "agree"), mk_resp("A", rmat, c("g_one", "g_two"), "all"),
                  mk_resp("B", RB, c("h_p", "h_q"), "agree"), mk_resp("B", RB, c("h_p", "h_q"), "all"))
reg <- tibble(iso3 = c(isos, "ZZZ"), group = c("g_one", "g_one", "g_two", "g_two", "g_one", "g_two"))

cat("# [1] MAP partition\n")
maps <- membership_map_partition(resp)
ma <- maps[maps$id == "A" & maps$chains == "agree", ]
expect(identical(ma$label, c("g_one", "g_one", "g_two", "g_two", "g_one")) &&
         isTRUE(all.equal(ma$p_max, c(0.9, 0.55, 0.8, 0.95, 1))), "argmax label and p_max per country")

cat("# [2] MAP == registry\n")
rg <- setNames(reg$group, reg$iso3)
expect(!throws(membership_check_map(ma, rg)), "the deployed cell's MAP is the registry: passes")
rg2 <- rg; rg2[["BBB"]] <- "g_two"
expect(throws(membership_check_map(ma, rg2), "BBB \\(cell g_one, registry g_two\\)"),
       "a MAP differing from the registry stops, naming the country")

cat("# [3] elpd difference vs loo_compare\n")
sim_loo <- function(seed, N = 5L, shift = 0) {
  set.seed(seed); mu <- rnorm(N, shift, 0.3)
  ll <- sapply(mu, function(m) dnorm(rnorm(300, m, 0.2), 0, 1, log = TRUE))
  suppressWarnings(loo::loo(ll, r_eff = rep(1, N)))
}
la <- sim_loo(1); lb <- sim_loo(2, shift = 0.3)
lc <- bind_rows(tibble(chains = "agree", id = "A", iso3 = isos, elpd = la$pointwise[, "elpd_loo"],
                       pareto_k = c(0.1, 0.8, 0.2, 0.3, 1.1)),
                tibble(chains = "agree", id = "B", iso3 = isos, elpd = lb$pointwise[, "elpd_loo"],
                       pareto_k = c(0.1, 0.2, 0.2, 0.3, 0.4)))
lc <- bind_rows(lc, mutate(lc, chains = "all"))
e <- membership_elpd_diff(lc, "B", "A", "agree")
cmp <- loo::loo_compare(la, lb)
expect(isTRUE(all.equal(abs(e$diff), abs(cmp[2, "elpd_diff"]))) && isTRUE(all.equal(e$se, cmp[2, "se_diff"])),
       "pairwise elpd difference and SE = loo::loo_compare()")
expect(e$n_k07 == 0L && e$n_k07_ref == 2L, "Pareto k > 0.7 counts of both cells")

cat("# [4] country, group, sensitivity tables\n")
m1 <- bind_rows(tibble(id = "A", iso3 = isos, quantity = "kcon", mix_mean = 1:5 / 10, mix_q2.5 = 0, mix_q97.5 = 1,
                       map_mean = 0.5, diff_mean = 1:5 / 10 - 0.5),
                tibble(id = "A", iso3 = isos, quantity = "kmu", mix_mean = 2, mix_q2.5 = 1, mix_q97.5 = 3,
                       map_mean = 2, diff_mean = 0))
ct <- membership_country_table(resp[resp$id == "A" & resp$chains == "agree", ], m1, reg)
expect(identical(sort(ct$iso3), sort(reg$iso3)) && isTRUE(ct$carried_over[ct$iso3 == "ZZZ"]) &&
         is.na(ct$p_max[ct$iso3 == "ZZZ"]) && is.na(ct$r_g_one[ct$iso3 == "ZZZ"]) &&
         is.na(ct$m1_kcon_mix[ct$iso3 == "ZZZ"]), "a carried-over country is listed with NA values")
expect(isTRUE(all.equal(ct$r_g_one[ct$iso3 == "BBB"], 0.55)) && isTRUE(all.equal(ct$m1_kcon_mix[ct$iso3 == "CCC"], 0.3)),
       "responsibility row and M1 mixture per country")
gt <- membership_group_table(ct)
expect(identical(gt$n_hard, c(3L, 2L)) && identical(gt$n_carried_over, c(0L, 1L)) &&
         isTRUE(all.equal(gt$soft_size, c(sum(rmat[, 1]), sum(rmat[, 2])))) &&
         isTRUE(all.equal(gt$p_max_min, c(0.55, 0.8))),
       "group summary: hard / carried-over members, soft size, min p_max")
lt <- tibble(id = c("A", "B", "A", "B"), chains = c("agree", "agree", "all", "all"), stacking_weight = c(0.7, 0.3, 0.6, 0.4))
sens <- membership_sensitivity(maps, rg[isos], lc, lt, "A")
sa <- sens[sens$id == "A" & sens$chains == "agree", ]
sb <- sens[sens$id == "B" & sens$chains == "agree", ]
expect(isTRUE(all.equal(sa$ari, 1)) && sa$n_agree_matched == 5L && sa$moved == "" && sa$elpd_diff == 0,
       "the deployed cell agrees with itself (ARI 1, nothing moves)")
# B's MAP: AAA h_q, BBB h_q, CCC h_p, DDD h_p, EEE h_q -> matched h_q->g_one, h_p->g_two
expect(sb$n_agree_label == 0L && sb$n_agree_matched == 5L && isTRUE(all.equal(sb$ari, 1)) &&
         isTRUE(all.equal(sb$stacking_weight, 0.3)),
       "another cell with other labels but the same partition: 0 by label, 5 after matching, ARI 1")

cat("# [5] M2\n")
x <- c(-1, 0, 1, 1.5); conds <- enso_conds_named(c(1, 1.5))
set.seed(3)
lr <- list(g_one = matrix(rnorm(4 * 50, 0.2), 4), g_two = matrix(rnorm(4 * 50, -0.3), 4))
R1 <- rbind(AAA = c(g_one = 1, g_two = 0), CCC = c(g_one = 0, g_two = 1))
m2 <- membership_m2(R1, lr, c(AAA = "g_one", CCC = "g_two"), x, conds)
mi <- tidyr::pivot_wider(dplyr::select(m2$irr, iso3, erf, enso_label, irr), names_from = erf, values_from = irr)
expect(isTRUE(all.equal(mi$mixture, mi$hard)), "degenerate responsibilities: mixture = the hard ERF")
R2 <- rbind(BBB = c(g_one = 0.25, g_two = 0.75))
m2b <- membership_m2(R2, lr, c(BBB = "g_one"), x, conds)
expect(isTRUE(all.equal(m2b$erf$logRR, apply(0.25 * lr$g_one[1:2, ] + 0.75 * lr$g_two[1:2, ], 1, median))),
       "mixture = sum_k r_ck logRR_k per draw (ERF on the non-condition points)")
expect(throws(membership_m2(R2, lr["g_one"], c(BBB = "g_one"), x, conds), "no ERF draws"),
       "a group without draws stops")

cat("# [6] derive_report_membership()\n")
td <- tempfile("mb"); dir.create(td); kd <- file.path(td, "cmp"); dir.create(kd)
wr <- function(d, f) write.csv(d, file.path(kd, f), row.names = FALSE)
wr(resp, "responsibilities.csv"); wr(m1, "country_mixture.csv"); wr(lc, "loo_country.csv"); wr(lt, "loo.csv")
wr(tibble(id = c("A", "B"), chains = 6, n_modes = 1, mode_sizes = "6", agreeing = "1,2,3"), "fits.csv")
wr(tibble(id = c("A", "B"), chains = "agree", n_draws = 100, ari_mean = c(0.9, 0.8)), "stability.csv")
wr(tibble(id = "A", iso3 = isos, group = "one", p = 1), "partition.csv")
saveRDS(list(settings = list(n_agree = 0L, mode_tol = 1.5)), file.path(kd, "loo.rds"))
gd <- file.path(td, "groupings"); dir.create(gd); write.csv(reg, file.path(gd, "reg_x.csv"), row.names = FALSE)
ctx <- list2env(list(KL_SWEEP_DIR = kd, KL_SWEEP_CELL = "A", GROUPING_PATH = file.path(gd, "reg_x.csv"),
                     MAN = list(groups = list(g_one = list(best_lag_occ = 0L), g_two = list(best_lag_occ = 2L))),
                     FIG1_TRACKS = "occurrence", opt = list(outdir = td), LAG_LOO_BATCH = list(default = "fake"),
                     ENSO_CONDS = enso_conds_named(c(-1, 1)), N_IRR_DRAWS = 40L), parent = emptyenv())
mb <- suppressMessages(derive_report_membership(ctx))
expect(identical(mb$cell, "A") && nrow(mb$country) == 6L && is.null(mb$mixture_m2) &&
         identical(mb$m2_skipped_tracks, "occurrence") && nrow(mb$sensitivity) == 4L,
       "end to end; M2 skipped (no sidecar) with the track recorded")
ctx$KL_SWEEP_CELL <- "C"
expect(throws(derive_report_membership(ctx), "not a fit"), "a cell absent from the sweep stops")
ctx$KL_SWEEP_CELL <- "A"
# With sidecars: one per group at the deployed lag (+ a second lag).
cbn <- paste0("enso_cb", 1:3)
vm <- function(l) list(single_lag = l, cb_col_names = cbn, range = c(-2, 2.5),
                       onebasis_attrs = list(fun = "ns", knots = c(-0.5, 0.5), Boundary.knots = c(-2, 2.5),
                                             intercept = FALSE, range = c(-2, 2.5)))
for (g in c("g_one", "g_two")) {
  lags <- c(0L, 2L); labs <- c("lag00", "lag02")
  sc <- list(schema = LAG_LOO_SIDECAR_SCHEMA, track = "occurrence", group = g, lags = setNames(lags, labs), n_obs = 3L,
             elpd_pointwise = matrix(-1, 3, 2, dimnames = list(NULL, labs)),
             pareto_k = matrix(0.1, 3, 2, dimnames = list(NULL, labs)),
             row_digest = setNames(c("d", "d"), labs),
             cb = setNames(lapply(lags, function(l) list(single_lag = l, cb_meta = vm(l),
                   b_draws = matrix(rnorm(60 * 3, if (g == "g_one") 0.2 else -0.2, 0.05), 60, 3,
                                    dimnames = list(NULL, cbn)))), labs))
  f <- makeLooCompareRunFile(ctx$opt, "loo_compare_lag_sidecar", batch = "fake_occ_grouped_enso_singlelag",
                             group = g, date = as.Date("2026-10-27"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); saveRDS(sc, f)
}
mb2 <- suppressMessages(derive_report_membership(ctx))
h <- mb2$mixture_m2
expect(!is.null(h) && setequal(unique(h$iso3), isos) && all(h$track == "occurrence") &&
         !length(mb2$m2_skipped_tracks) && nrow(h) == 5L * 2L * 2L,
       "M2 with sidecars: per fitted country x {mixture, hard} x ENSO condition")
eee <- tidyr::pivot_wider(dplyr::select(dplyr::filter(h, iso3 == "EEE"), erf, enso_label, irr),
                          names_from = erf, values_from = irr)
expect(isTRUE(all.equal(eee$mixture, eee$hard)), "a country with r = (1, 0) has mixture = hard")

cat("# [7] (soft) the deployed cell on disk\n")
kd_real <- "generated_data/kl_sweep_compare/kl_sweep_classes_prodtu_nz60"
mc <- "generated_data/kl_groupings/enso_grouping_kl_nz60_ns5k5_membership.csv"
if (file.exists(file.path(kd_real, "responsibilities.csv")) && file.exists(mc)) {
  r <- read.csv(file.path(kd_real, "responsibilities.csv"))
  mp <- membership_map_partition(r[r$id == "ns5_K5" & r$chains == "agree", ])
  m <- read.csv(mc)
  expect(max(abs(m$p_max - mp$p_max[match(m$iso3, mp$iso3)])) < 1e-12,
         "ns5_K5 agreeing-chain p_max = gen_grouping's membership CSV (to 1e-12)")
} else cat("  SKIP  local sweep outputs absent\n")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
