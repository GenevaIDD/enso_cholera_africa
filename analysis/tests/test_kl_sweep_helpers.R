# test_kl_sweep_helpers.R — HARD: utils/kl_sweep.R on Stan-free fake kl_fit
#   objects.
#
# Usage:
#   Rscript analysis/tests/test_kl_sweep_helpers.R
#
# Details:
#   Fake fits from tests/helpers_synth.R::make_fake_kl_fit(); no sampling.
#   Covers chain modes, agreeing chains, per-fit rows, the comparability
#   guard, the country-block LOO comparison with stacking weights,
#   membership uncertainty (responsibility rows, the M1 country mixture,
#   adjusted Rand index, class matching, partition stability) and the
#   make_kl_sweep_run_file() paths.

suppressPackageStartupMessages(library(dplyr))
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c(
  "base",
  "config",
  "paths",
  "dlnm",
  "dlnm_kernel",
  "diagnostics",
  "kl_sweep"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
source("analysis/tests/helpers_synth.R")

# 6 chains: chains 1, 2, 4, 6 at centres (0, 9), chains 3 and 5 at (-6, 4).
cen <- rbind(c(0, 9), c(0, 9), c(-6, 4), c(0.2, 9), c(-6, 4.2), c(0.1, 9))
fa <- make_fake_kl_fit(
  cen,
  lp = c(-1000, -1000, -990, -1000, -990, -1000),
  ll_shift = 2,
  seed = 1
)
fb <- make_fake_kl_fit(cen, ll_shift = 0, seed = 2)

cat("# [1] chain modes / agreeing chains\n")
md <- kl_chain_modes(fa, mode_tol = 1.5)
expect(identical(md$sizes, c(4L, 2L)), "two modes of sizes 4 and 2")
expect(
  identical(unname(md$mode[c("1", "2", "4", "6")]), rep(1L, 4)),
  "dominant mode = the 4 chains at (0, 9)"
)
ag <- kl_agreeing_chains(fa, 1.5, 0L, md)
expect(
  identical(ag$chains, c(1L, 2L, 4L, 6L)),
  "agreeing = dominant mode (size beats lp)"
)
ag2 <- kl_agreeing_chains(fa, 1.5, 2L, md)
expect(
  length(ag2$chains) == 2L && ag2$kmu_spread <= 0.25,
  "n_agree = 2 picks the closest pair"
)
expect(
  inherits(
    tryCatch(kl_agreeing_chains(fa, 1.5, 7L, md), error = identity),
    "error"
  ),
  "n_agree > chains is an error"
)

cat("# [2] per-fit row and diagnostics tables\n")
r <- tibble(id = "x_K2", basis = "ns5", mix = FALSE, model = "cfg")
row <- kl_fit_row(fa, r, md, ag)
expect(
  nrow(row) == 1L && row$n_modes == 2L && row$mode_sizes == "4+2",
  "row: modes"
)
expect(
  isTRUE(row$max_rhat_agree < row$max_rhat_all),
  "row: agreeing chains lower R-hat"
)
dt <- kl_fit_diagnostics_tables(fa, "x_K2", ag)
expect(
  setequal(unique(dt$variables$chains), c("all", "agree")),
  "diagnostics: all + agree"
)
expect(
  sum(dt$chains$agreeing) == 4L && nrow(dt$chains) == 6L,
  "diagnostics_chains: 6 rows, 4 agreeing"
)
expect(
  !any(grepl("^kl_log", dt$variables$variable)),
  "KL_SWEEP_EXCL drops kl_log_*"
)
sub <- kl_subset_chains(fa, ag$chains)
pr <- kl_partition_rows(sub, r)
expect(
  nrow(pr) == length(fa$country_levels) && all(pr$p > 0.8),
  "partition: one row per country, p = 0.9"
)

cat("# [3] LOO comparison\n")
lo <- kl_sweep_loo(list(a = fa, b = fb), list(a = ag$chains, b = ag$chains))
ta <- filter(lo$table, chains == "agree")
expect(abs(sum(ta$stacking_weight) - 1) < 1e-6, "stacking weights sum to 1")
expect(ta$elpd_diff_loo[ta$id == "a"] == 0, "the better fit has elpd_diff 0")
expect(
  nrow(lo$country) == 2L * 2L * length(fa$country_levels),
  "per-country rows: chains x fits x countries"
)
fc <- make_fake_kl_fit(cen, countries = sprintf("D%02d", 1:6))
expect(
  inherits(
    tryCatch(kl_sweep_check_comparable(list(fa, fc)), error = identity),
    "error"
  ),
  "different countries are refused"
)

cat(
  "# [4] membership: responsibilities, M1 mixture, ARI, matching, stability\n"
)
prior <- tibble(
  iso3 = fa$country_levels,
  group = rep(c("g1", "g2"), length.out = 6)
)
rr <- kl_responsibility_rows(sub, r, "agree", prior)
expect(
  nrow(rr) == fa$K * length(fa$country_levels),
  "responsibilities: K x C rows"
)
sums <- tapply(rr$mean, rr$iso3, sum)
expect(all(abs(sums - 1) < 1e-12), "responsibilities: each country sums to 1")
expect(
  identical(rr$label, fa$class_labels[rr$class]),
  "responsibilities: label = class_labels[class]"
)
expect(
  all(rr$prior_group == prior$group[match(rr$iso3, prior$iso3)]),
  "responsibilities: prior_group"
)
top <- rr %>% group_by(iso3) %>% slice_max(mean, n = 1) %>% ungroup()
expect(
  all(top$class == ((match(top$iso3, fa$country_levels) - 1L) %% 2L + 1L)),
  "responsibilities: the high-responsibility class is the generating one"
)
fl <- make_fake_kl_fit(cen, labels = c("enso_east", "enso_west"))
expect(
  all(
    kl_responsibility_rows(fl, r, "all")$label %in% c("enso_east", "enso_west")
  ) &&
    all(is.na(kl_responsibility_rows(fl, r, "all")$prior_group)),
  "responsibilities: fit labels, NA prior"
)

mx <- kl_country_mixture_rows(sub, r, prior)
expect(
  setequal(unique(mx$quantity), c("kcon", "kmu")) && nrow(mx) == 2L * 6L,
  "mixture: kcon + kmu per country"
)
R <- exp(.kl_draw_array(sub, "kl_log_resp", c(sub$K, 6)))
KC <- .kl_draw_array(sub, "kcon", sub$K)
want <- vapply(1:6, function(c) mean(rowSums(R[,, c] * KC)), 1)
got <- mx$mix_mean[mx$quantity == "kcon"]
expect(
  max(abs(got - want)) < 1e-12,
  "mixture: mean of sum_k r_ck kcon_k per draw"
)
mapc <- (0:5 %% 2) + 1
expect(
  max(abs(mx$map_mean[mx$quantity == "kcon"] - colMeans(KC)[mapc])) < 1e-12,
  "mixture: MAP = kcon of the argmax class"
)
expect(
  max(abs(mx$diff_mean - (mx$mix_mean - mx$map_mean))) < 1e-12,
  "mixture: diff = mixture - MAP"
)
# 0.9 / 0.1 responsibilities, kcon = 0.5 k: mixture - MAP = +-0.1 * 0.5
dk <- mx$diff_mean[mx$quantity == "kcon"]
expect(
  all(abs(abs(dk) - 0.05) < 0.02),
  "mixture: |mixture - MAP| = 0.1 x the kcon gap"
)
nokmu <- sub
nokmu$form <- "dlnm"
nokmu$draws <- posterior::subset_draws(
  sub$draws,
  variable = "^kmu",
  regex = TRUE,
  exclude = TRUE
)
mn <- kl_country_mixture_rows(nokmu, r)
expect(
  all(is.na(mn$mix_mean[mn$quantity == "kmu"])) &&
    all(grepl("no kmu", mn$note[mn$quantity == "kmu"])) &&
    !anyNA(mn$mix_mean[mn$quantity == "kcon"]),
  "mixture: absent kmu -> NA with a reason, kcon kept"
)

a1 <- c(1, 1, 2, 2, 3, 3)
expect(
  adjusted_rand_index(a1, c("b", "b", "c", "c", "a", "a")) == 1,
  "ARI = 1 under relabelling"
)
expect(
  abs(adjusted_rand_index(c(0, 0, 1, 1), c(0, 0, 1, 2)) - 4 / 7) < 1e-12,
  "ARI textbook value 4/7"
)
expect(
  abs(adjusted_rand_index(c(0, 0, 1, 2), c(0, 0, 1, 1)) - 4 / 7) < 1e-12,
  "ARI symmetric"
)
expect(
  adjusted_rand_index(rep(1, 5), rep(2, 5)) == 1,
  "ARI: identical trivial partitions = 1"
)
expect(
  adjusted_rand_index(rep(1, 4), 1:4) == 0,
  "ARI: one cluster vs singletons = 0"
)

mc <- kl_match_classes(
  c("A", "A", "B", "B", "C", "C"),
  c("z", "z", "x", "x", "y", "x")
)
expect(
  identical(unname(mc$map[c("A", "B", "C")]), c("z", "x", "y")) &&
    mc$n_agree == 5L,
  "matching: optimal one-to-one map, 5 of 6 agree"
)
mu <- kl_match_classes(c("A", "A", "B", "B", "C"), c("u", "u", "v", "v", "v"))
expect(
  is.na(mu$map[["C"]]) && mu$n_agree == 4L,
  "matching: unequal class counts leave a label unmatched"
)
expect(
  nrow(.kl_perms(4)) == 24L &&
    !anyDuplicated(apply(.kl_perms(4), 1, paste, collapse = "")),
  "perms: 4! distinct"
)

st <- kl_partition_stability(sub, r, "agree", seed = 3)
expect(
  st$ari_q97.5 <= 1 && st$ari_q2.5 < 1 && st$moved_mean > 0,
  "stability: soft responsibilities move countries"
)
st2 <- kl_partition_stability(sub, r, "agree", seed = 3)
expect(identical(st, st2), "stability: reproducible under its seed")
set.seed(99)
u0 <- runif(1)
set.seed(99)
invisible(kl_partition_stability(sub, r, seed = 3))
u1 <- runif(1)
expect(u0 == u1, "stability: caller RNG state restored")
dg <- sub
dd <- as.data.frame(dg$draws)
for (k in 1:2) {
  for (c in 1:6) {
    dd[[sprintf("kl_log_resp[%d,%d]", k, c)]] <- if (mapc[c] == k) 0 else -Inf
  }
}
dg$draws <- posterior::as_draws_df(
  dd[, setdiff(names(dd), c(".chain", ".iteration", ".draw"))] |>
    cbind(.chain = dd$.chain, .iteration = dd$.iteration)
)
sd1 <- kl_partition_stability(dg, r, "agree")
expect(
  sd1$ari_q2.5 == 1 && sd1$p_ari1 == 1 && sd1$moved_mean == 0,
  "stability: degenerate responsibilities -> ARI 1"
)

cat("# [5] paths\n")
p <- vapply(
  c(
    "kl_sweep_fits",
    "kl_sweep_loo",
    "kl_sweep_diagnostics",
    "kl_sweep_rds",
    "kl_sweep_responsibilities",
    "kl_sweep_country_mixture",
    "kl_sweep_stability"
  ),
  function(k) {
    make_kl_sweep_run_file(list(), k, batch = "analysis/configs/kl_sweep_x")
  },
  ""
)
expect(
  all(
    startsWith(p, "generated_data//kl_sweep_compare/kl_sweep_x/") |
      startsWith(p, "generated_data/kl_sweep_compare/kl_sweep_x/")
  ),
  "paths under kl_sweep_compare/<tag>/"
)
expect(
  !anyDuplicated(p) && grepl("loo\\.rds$", p[["kl_sweep_rds"]]),
  "distinct paths; rds = loo.rds"
)

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
quit(status = if (FAIL > 0L) 1L else 0L)
