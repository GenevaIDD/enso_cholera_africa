# test_kl_lp_equivalence.R — the patched brms_kl program equals brms at K = 1
#
# For one kernel (kernel.classes = 1, no prior partition, no ICAR), the
# program produced by utils/dlnm_kernel.R::patch_kl_latent_class_stancode()
# must have the same log density as brms's own program for the non-linear
# kernel formula: same likelihood (one class: log_sum_exp of one term), same
# priors, same transforms (brms's bounded kmu and the patched ordered[1]
# mapped through inv_logit carry the same Jacobian). This pins the whole
# patch -- data assembly, kernel basis matrix, weights, linear predictor --
# to brms, without sampling.
#
# Compiles both programs (fresh, in tempdir -- see below), then
# evaluates log_prob (with Jacobian) at several random parameter points.
# HARD gate (deterministic); registered under !--quick in run_all.R because
# of the two compilations. Skips (exit 0) without CmdStan.
#
# Also: kernel.form discrete -- the program's sum over ordered (ties allowed)
# centre combinations equals a brute-force R sum over every labelled tuple
# (plus seeds and random slopes); kernel.form gaussian on the reduce_sum
# country path (seeds, fixed width, random slopes) equals a brute force;
# identify repel (label-free registry pull, pairwise repulsion, the Jacobian
# of its bounded centres and kappa) equals a brute force on both paths.
#
# Run: Rscript analysis/tests/test_kl_lp_equivalence.R

suppressPackageStartupMessages({ library(dplyr); library(lubridate) })
if (inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")) {
  message("SKIP test_kl_lp_equivalence: CmdStan not available"); quit(status = 0)
}
for (m in c("base", "config", "formula", "priors", "paths", "enso", "dlnm", "dlnm_kernel"))
  source(sprintf("analysis/utils/%s.R", m))
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_synth.R")
PASS <- 0L; FAIL <- 0L
cat("# brms_kl: patched program == brms program at K = 1 (log density)\n")

# Fresh directory every run: cmdstanr can only attach log_prob() model methods
# to an executable compiled in the current session.
stan_dir <- file.path(tempdir(), "stan_kl_lp")
dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)

panel <- make_synthetic_kl_panel(n_country = 4L, adm1_per_country = 2L, n_year = 3L,
                                 lag_min = -3L, lag_max = 3L, K = 1L, mu_true = 1,
                                 sigma_true = 1, beta_true = 0.5, log_h0_intercept = -4)
opt  <- kl_synthetic_opt(lag_min = -3L, lag_max = 3L)
kcfg <- resolve_kernel_cfg(opt)
m    <- kl_sort_by_country(assemble_kl_model(
  opt, panel, list(re_uncorrelated = TRUE), .occ_baseline_block("cyclic"), kcfg))
J <- m$meta$df

compile <- function(code, stem) {
  f <- file.path(stan_dir, sprintf("%s_%s.stan", stem,
                                   substr(digest::digest(code, "sha1", serialize = FALSE), 1, 12)))
  if (!file.exists(f)) writeLines(code, f)
  cmdstanr::cmdstan_model(f, compile_model_methods = TRUE, quiet = TRUE)
}
fit1 <- function(mod, data) {
  f <- mod$sample(data = data, chains = 1, iter_warmup = 1, iter_sampling = 1,
                  refresh = 0, seed = 1, show_messages = FALSE, show_exceptions = FALSE,
                  diagnostics = NULL)
  f$init_model_methods(verbose = FALSE)
  f
}

sc_brms <- brms::stancode(m$formula, data = m$data, family = m$family, prior = m$priors)
sd_brms <- unclass(brms::standata(m$formula, data = m$data, family = m$family, prior = m$priors))
code_kl <- patch_kl_latent_class_stancode(sc_brms, m$meta)
pg      <- kl_prior_partition(kcfg, levels(droplevels(m$data$country)))
sd_kl   <- build_kl_latent_standata(m, pg, NULL)
expect(pg$K == 1L && !pg$use, "K = 1, no prior partition")

f_brms <- fit1(compile(as.character(sc_brms), "brms_k1"), sd_brms)
f_kl   <- fit1(compile(as.character(code_kl), "occ_kl_k1"), sd_kl)

set.seed(42)
lin_point <- function() {
  p <- list(b_lin = rnorm(sd_brms$K_lin, 0, 0.5))
  if (!is.null(sd_brms$nb_lin_1)) {
    p$zs_lin_1_1 <- rnorm(sd_brms$knots_lin_1[1])
    p$sds_lin_1  <- array(runif(sd_brms$nb_lin_1, 0.2, 1.5))
  }
  for (g in 1:3) {
    M <- sd_brms[[sprintf("M_%d", g)]]; N <- sd_brms[[sprintf("N_%d", g)]]
    p[[sprintf("sd_%d", g)]] <- array(runif(M, 0.1, 1))
    p[[sprintf("z_%d", g)]]  <- matrix(rnorm(M * N), M, N)
  }
  p
}
lo <- kcfg$lag_min; hi <- kcfg$lag_max
for (rep in 1:4) {
  lp   <- lin_point()
  kmu  <- runif(1, lo + 0.2, hi - 0.2); ksd <- runif(1, 0.5, 3); kb <- rnorm(J, 0, 0.5)
  pb <- c(lp, list(b_kmu = array(kmu), b_ksd = array(ksd)),
          stats::setNames(lapply(kb, array), sprintf("b_kb%d", seq_len(J))))
  pk <- c(lp, list(kmu_ord = array(qlogis((kmu - lo) / (hi - lo))), kmu_free = numeric(0),
                   ksd = array(ksd), kbeta_free = matrix(kb, J, 1), kpi = array(1),
                   kc = numeric(0), keta = matrix(numeric(0), J - 1, 0),
                   kz = matrix(numeric(0), 0, 1), khs_tau = numeric(0),
                   khs_lambda = matrix(numeric(0), 0, 1), khs_c2 = numeric(0),
                   kkappa_free = numeric(0), kphi = matrix(numeric(0), 0, 0), ktau = numeric(0)))
  lpb <- f_brms$log_prob(f_brms$unconstrain_variables(pb), jacobian = TRUE)
  lpk <- f_kl$log_prob(f_kl$unconstrain_variables(pk), jacobian = TRUE)
  cat(sprintf("  point %d: brms %.8f  patched %.8f  diff %.2e\n", rep, lpb, lpk, lpb - lpk))
  expect(abs(lpb - lpk) < 1e-6 * max(1, abs(lpb)),
         sprintf("point %d: log density equal (|diff| < 1e-6 relative)", rep))
}

# --- [2] form discrete: the ordered sum == brute force over labelled tuples --------
# For K = 2 classes, widths {0, 1.5} and 7 candidate centres, the program sums
# over the 28 ordered (ties allowed) centre pairs with prior log(K!/prod n_i!).
# That must equal the brute-force sum over all 7^2 labelled pairs l, class k
# taking the k-th smallest centre of l (each multiset is hit K!/prod n_i! times).
# R computes it from scratch (lin part X_lin b_lin + offset: random effects and
# spline coefficients at 0), and the difference between two parameter points
# that share the lin part is compared with the program's (jacobian = FALSE,
# Dirichlet(1) constant, so only the kbeta prior and the likelihood move).
cat("\n# form discrete: ordered-combination sum == labelled brute force\n")
optd  <- kl_synthetic_opt(lag_min = -3L, lag_max = 3L,
                          kernel = list(classes = 2L, form = "discrete",
                                        discrete = list(widths = c(0, 1.5))))
kcfgd <- resolve_kernel_cfg(optd)
md    <- kl_sort_by_country(assemble_kl_model(
  optd, panel, list(re_uncorrelated = TRUE), .occ_baseline_block("cyclic"), kcfgd))
pgd   <- kl_prior_partition(kcfgd, levels(droplevels(md$data$country)))
sdd   <- build_kl_latent_standata(md, pgd, NULL)
mod_d <- compile(as.character(patch_kl_latent_class_stancode(
  brms::stancode(md$formula, data = md$data, family = md$family, prior = md$priors), md$meta)),
  "occ_kl_disc")
f_d   <- fit1(mod_d, sdd)
expect(sdd$kl_M == choose(7 + 1, 2) && sdd$kl_S == 2L, "discrete: 28 ordered pairs, 2 widths")

KB <- as.matrix(md$data[, md$meta$cb_col_names])            # N x (NL * J), lag-major
Jd <- md$meta$df; NLd <- length(md$meta$lags)
ctry <- as.integer(droplevels(md$data$country)); Cd <- max(ctry)
brute_ll <- function(b_lin, kbeta, kpi, fix = rep(0L, Cd), D = NULL) {
  lin <- as.numeric(sdd$X_lin %*% as.numeric(b_lin)) + as.numeric(sdd$offsets_lin)
  # membership log-probabilities per country: log kpi, or a point mass on the
  # seed's class (kernel.pin)
  LSM <- t(vapply(fix, function(f) if (f == 0L) log(kpi) else replace(rep(-Inf, 2), f, 0), numeric(2)))
  # country log-lik of class k under centre index i, width s
  LL <- array(NA_real_, c(Cd, 2, sdd$kl_NCt, sdd$kl_S))
  for (s in seq_len(sdd$kl_S)) for (i in seq_len(sdd$kl_NCt)) {
    w <- sdd$kl_Wq[, (s - 1) * sdd$kl_NCt + i]
    Bq <- Reduce(`+`, lapply(seq_len(NLd), function(l) w[l] * KB[, (l - 1) * Jd + seq_len(Jd)]))
    for (k in 1:2) {
      eta <- lin + as.numeric(Bq %*% kbeta[, k]) + (if (is.null(D)) 0 else rowSums(Bq * D))
      LL[, k, i, s] <- rowsum(dpois(sdd$Y, exp(eta), log = TRUE), ctry)[, 1]
    }
  }
  terms <- c()
  for (s in seq_len(sdd$kl_S)) for (a in seq_len(sdd$kl_NCt)) for (b in seq_len(sdd$kl_NCt)) {
    srt <- sort(c(a, b))
    lc  <- vapply(seq_len(Cd), function(c) {
      v <- LSM[c, ] + c(LL[c, 1, srt[1], s], LL[c, 2, srt[2], s]); max(v) + log(sum(exp(v - max(v))))
    }, numeric(1))
    terms <- c(terms, -log(sdd$kl_S) + sum(lc))
  }
  max(terms) + log(sum(exp(terms - max(terms))))
}
pd <- function(b_lin, kbeta, kpi, f = f_d, re = NULL) {
  p <- list(b_lin = b_lin, zs_lin_1_1 = rep(0, sdd$knots_lin_1[1]),
            sds_lin_1 = array(rep(0.5, sdd$nb_lin_1)),
            kmu_ord = numeric(0), kmu_free = numeric(0), ksd = numeric(0),
            kbeta_free = kbeta, kpi = kpi, kc = numeric(0), keta = matrix(numeric(0), Jd - 1, 0),
            kz = matrix(numeric(0), 0, 2), khs_tau = numeric(0), khs_lambda = matrix(numeric(0), 0, 2),
            khs_c2 = numeric(0), kkappa_free = numeric(0), kphi = matrix(numeric(0), 0, 2),
            ktau = numeric(0), krho = numeric(0),
            kre_zc = matrix(numeric(0), 0, Cd), kre_za = matrix(numeric(0), 0, sdd$kl_NA),
            kre_sdc = numeric(0), kre_sda = numeric(0))
  if (!is.null(re)) p[names(re)] <- re
  for (g in 1:2) {
    M <- sdd[[sprintf("M_%d", g)]]; N <- sdd[[sprintf("N_%d", g)]]
    p[[sprintf("sd_%d", g)]] <- array(rep(0.5, M)); p[[sprintf("z_%d", g)]] <- matrix(0, M, N)
  }
  f$log_prob(f$unconstrain_variables(p), jacobian = FALSE)
}
set.seed(7)
b_lin <- array(c(-4))
for (rep in 1:3) {
  kb1 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2); kb2 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2)
  pi1 <- { x <- runif(2); x / sum(x) }; pi2 <- { x <- runif(2); x / sum(x) }
  d_stan  <- pd(b_lin, kb2, pi2) - pd(b_lin, kb1, pi1)
  d_brute <- (brute_ll(b_lin, kb2, pi2) + sum(dnorm(kb2, 0, kcfgd$prior_beta_sd, log = TRUE))) -
             (brute_ll(b_lin, kb1, pi1) + sum(dnorm(kb1, 0, kcfgd$prior_beta_sd, log = TRUE)))
  cat(sprintf("  point %d: program %.8f  brute force %.8f  diff %.2e\n", rep, d_stan, d_brute,
              d_stan - d_brute))
  expect(abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
         sprintf("discrete point %d: log-density difference equals the labelled brute force", rep))
}

# Seed anchoring (kernel.pin): countries 1 and 4 fixed to classes 2 and 1.
fixd  <- c(2L, 0L, 0L, 1L)
f_dfx <- fit1(mod_d, build_kl_latent_standata(md, pgd, NULL, fix = fixd))   # same program, seeded data
for (rep in 1:2) {
  kb1 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2); kb2 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2)
  pi1 <- { x <- runif(2); x / sum(x) }; pi2 <- { x <- runif(2); x / sum(x) }
  d_stan  <- pd(b_lin, kb2, pi2, f_dfx) - pd(b_lin, kb1, pi1, f_dfx)
  d_brute <- (brute_ll(b_lin, kb2, pi2, fixd) + sum(dnorm(kb2, 0, kcfgd$prior_beta_sd, log = TRUE))) -
             (brute_ll(b_lin, kb1, pi1, fixd) + sum(dnorm(kb1, 0, kcfgd$prior_beta_sd, log = TRUE)))
  cat(sprintf("  seeded point %d: program %.8f  brute force %.8f  diff %.2e\n", rep, d_stan, d_brute,
              d_stan - d_brute))
  expect(abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
         sprintf("discrete + seeds point %d: log-density difference equals the brute force", rep))
}

# Random slopes (kernel.re_slope): country + ADM1 deviations of the exposure
# coefficients, applied at each class's lag. Fixed deviations at both points, so
# their priors cancel and the brute force adds rowSums(Bq * D) to each row.
kcfg_re <- kcfgd; kcfg_re$re_slope <- list(enable = TRUE, sd_country = 0.35, sd_adm1 = 0.35)
sd_re <- build_kl_latent_standata(modifyList(md, list(kcfg = kcfg_re)), pgd, NULL)
f_dre <- fit1(mod_d, sd_re)                                     # same program, re switched on
re <- list(kre_zc = matrix(rnorm(Jd * Cd), Jd, Cd), kre_za = matrix(rnorm(Jd * sd_re$kl_NA), Jd, sd_re$kl_NA),
           kre_sdc = array(runif(Jd, 0.1, 0.5)), kre_sda = array(runif(Jd, 0.1, 0.5)))
Dre <- t(as.numeric(re$kre_sdc) * re$kre_zc)[ctry, , drop = FALSE] +
       t(as.numeric(re$kre_sda) * re$kre_za)[as.integer(sd_re$kl_adm1), , drop = FALSE]
for (rep in 1:2) {
  kb1 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2); kb2 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2)
  pi1 <- { x <- runif(2); x / sum(x) }; pi2 <- { x <- runif(2); x / sum(x) }
  d_stan  <- pd(b_lin, kb2, pi2, f_dre, re) - pd(b_lin, kb1, pi1, f_dre, re)
  d_brute <- (brute_ll(b_lin, kb2, pi2, D = Dre) + sum(dnorm(kb2, 0, kcfgd$prior_beta_sd, log = TRUE))) -
             (brute_ll(b_lin, kb1, pi1, D = Dre) + sum(dnorm(kb1, 0, kcfgd$prior_beta_sd, log = TRUE)))
  cat(sprintf("  random-slope point %d: program %.8f  brute force %.8f  diff %.2e\n", rep, d_stan, d_brute,
              d_stan - d_brute))
  expect(abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
         sprintf("discrete + random slopes point %d: log-density difference equals the brute force", rep))
}
# deviations on vs off at the same point: program difference = likelihood
# difference (brute force) + the standard-normal prior of the deviations
re0 <- modifyList(re, list(kre_zc = 0 * re$kre_zc, kre_za = 0 * re$kre_za))
d0  <- pd(b_lin, kb2, pi2, f_dre, re) - pd(b_lin, kb2, pi2, f_dre, re0)
d0b <- brute_ll(b_lin, kb2, pi2, D = Dre) - brute_ll(b_lin, kb2, pi2) +
       sum(dnorm(re$kre_zc, log = TRUE)) + sum(dnorm(re$kre_za, log = TRUE)) -
       sum(dnorm(0 * re$kre_zc, log = TRUE)) - sum(dnorm(0 * re$kre_za, log = TRUE))
cat(sprintf("  deviations on vs off: program %.8f  brute force %.8f  diff %.2e\n", d0, d0b, d0 - d0b))
expect(abs(d0) > 1e-3 && abs(d0 - d0b) < 1e-6 * max(1, abs(d0b)),
       "random slopes: on-vs-off log-density difference equals the brute force (deviations used)")

# --- [3] form gaussian, country path (reduce_sum) + seeds + fixed width + slopes --
# K = 2, kernel width fixed at 1.5, seed countries fixed to their class, country
# and ADM1 deviations on the exposure coefficients: the reduce_sum country path
# (kl_rs) against an R brute force at two parameter points sharing the lin
# part, kmu and the deviations (only kbeta and kpi move).
cat("\n# form gaussian, country path: reduce_sum log density == brute force\n")
optg  <- kl_synthetic_opt(lag_min = -3L, lag_max = 3L,
                          kernel = list(classes = 2L, form = "gaussian", ksd_fixed = 1.5,
                                        re_slope = list(enable = TRUE)))
kcfgg <- resolve_kernel_cfg(optg)
mg    <- kl_sort_by_country(assemble_kl_model(
  optg, panel, list(re_uncorrelated = TRUE), .occ_baseline_block("cyclic"), kcfgg))
pgg   <- kl_prior_partition(kcfgg, levels(droplevels(mg$data$country)))
fixg  <- c(1L, 0L, 0L, 2L)
sdg   <- build_kl_latent_standata(mg, pgg, NULL, fix = fixg)
expect(sdg$kl_rs == 1L && sdg$kl_re == 1L && sdg$kl_ksd_fix == 1L && sdg$kl_mu_ordered == 1L,
       "gaussian: country path, slopes, fixed width (no prior partition -> ordered centres)")
mod_g <- compile(as.character(patch_kl_latent_class_stancode(
  brms::stancode(mg$formula, data = mg$data, family = mg$family, prior = mg$priors), mg$meta)),
  "occ_kl_gauss_rs")
f_g <- fit1(mod_g, sdg)
# the same model with every country split by rows (threads = 8 > N / (2 * 8) rule):
# the per-class row reduce_sums must give the same log density
sdg_split <- build_kl_latent_standata(mg, pgg, NULL, fix = fixg, threads = 8L)
expect(sdg_split$kl_nbig == 4L && sdg_split$kl_nsmall == 0L && sdg$kl_nbig == 0L,
       "split variant: all 4 countries split by rows; default none")
f_gs <- fit1(mod_g, sdg_split)
KBg <- as.matrix(mg$data[, mg$meta$cb_col_names]); NLg <- length(mg$meta$lags)
ctg <- as.integer(droplevels(mg$data$country)); Cg <- max(ctg)
reg_ <- list(kre_zc = matrix(rnorm(Jd * Cg), Jd, Cg), kre_za = matrix(rnorm(Jd * sdg$kl_NA), Jd, sdg$kl_NA),
             kre_sdc = array(runif(Jd, 0.1, 0.5)), kre_sda = array(runif(Jd, 0.1, 0.5)))
Dg <- t(as.numeric(reg_$kre_sdc) * reg_$kre_zc)[ctg, , drop = FALSE] +
      t(as.numeric(reg_$kre_sda) * reg_$kre_za)[as.integer(sdg$kl_adm1), , drop = FALSE]
kmug <- c(-1.3, 1.7)
brute_g <- function(kbeta, kpi) {
  lin <- as.numeric(sdg$X_lin %*% as.numeric(b_lin)) + as.numeric(sdg$offsets_lin)
  W <- kl_weights(kmug, c(1.5, 1.5), mg$meta$lags)
  LSM <- t(vapply(fixg, function(f) if (f == 0L) log(kpi) else replace(rep(-Inf, 2), f, 0), numeric(2)))
  LL <- sapply(1:2, function(k) {
    Bw <- Reduce(`+`, lapply(seq_len(NLg), function(l) W[k, l] * KBg[, (l - 1) * Jd + seq_len(Jd)]))
    eta <- lin + as.numeric(Bw %*% kbeta[, k]) + rowSums(Bw * Dg)
    rowsum(dpois(sdg$Y, exp(eta), log = TRUE), ctg)[, 1]
  })
  sum(vapply(seq_len(Cg), function(c) { v <- LSM[c, ] + LL[c, ]; max(v) + log(sum(exp(v - max(v)))) }, 0))
}
pg_ <- function(kbeta, kpi, f = f_g) {
  p <- list(b_lin = b_lin, zs_lin_1_1 = rep(0, sdg$knots_lin_1[1]), sds_lin_1 = array(rep(0.5, sdg$nb_lin_1)),
            kmu_ord = array(qlogis((kmug + 3) / 6)), kmu_free = numeric(0), ksd = numeric(0), kbeta_free = kbeta, kpi = kpi,
            kc = numeric(0), keta = matrix(numeric(0), Jd - 1, 0), kz = matrix(numeric(0), 0, 2),
            khs_tau = numeric(0), khs_lambda = matrix(numeric(0), 0, 2), khs_c2 = numeric(0),
            kkappa_free = numeric(0), kphi = matrix(numeric(0), 0, 1), ktau = numeric(0), krho = numeric(0))
  for (g in 1:2) {
    M <- sdg[[sprintf("M_%d", g)]]; N <- sdg[[sprintf("N_%d", g)]]
    p[[sprintf("sd_%d", g)]] <- array(rep(0.5, M)); p[[sprintf("z_%d", g)]] <- matrix(0, M, N)
  }
  p[names(reg_)] <- reg_
  f$log_prob(f$unconstrain_variables(p), jacobian = FALSE)
}
for (rep in 1:2) {
  kb1 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2); kb2 <- matrix(rnorm(2 * Jd, 0, 0.6), Jd, 2)
  pi1 <- { x <- runif(2); x / sum(x) }; pi2 <- { x <- runif(2); x / sum(x) }
  d_stan  <- pg_(kb2, pi2) - pg_(kb1, pi1)
  d_brute <- (brute_g(kb2, pi2) + sum(dnorm(kb2, 0, kcfgg$prior_beta_sd, log = TRUE))) -
             (brute_g(kb1, pi1) + sum(dnorm(kb1, 0, kcfgg$prior_beta_sd, log = TRUE)))
  cat(sprintf("  gaussian point %d: program %.8f  brute force %.8f  diff %.2e\n", rep, d_stan, d_brute,
              d_stan - d_brute))
  expect(abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
         sprintf("gaussian country path point %d: log-density difference equals the brute force", rep))
  d_split <- pg_(kb2, pi2, f_gs) - pg_(kb1, pi1, f_gs)
  cat(sprintf("  gaussian split point %d: program %.8f  brute force %.8f  diff %.2e\n", rep, d_split, d_brute,
              d_split - d_brute))
  expect(abs(d_split - d_brute) < 1e-6 * max(1, abs(d_brute)),
         sprintf("gaussian row-split path point %d: log-density difference equals the brute force", rep))
  expect(abs(pg_(kb2, pi2, f_gs) - pg_(kb2, pi2)) < 1e-6 * max(1, abs(pg_(kb2, pi2))),
         sprintf("gaussian point %d: row-split and country log density identical", rep))
}

# --- [4] identify repel: label-free registry pull + repulsive prior ---------------
# K = 3 exchangeable classes, kernel width fixed 1.5, registry pull over 2 groups
# (C01, C02 -> g_east; C03 -> g_south; C04 outside). Brute force: the mixture
# likelihood with klogit = log kpi + kappa * A Z, plus every prior term that
# moves between two points sharing the lin part (kbeta normal, kappa half-normal,
# Z standard normal, the pairwise repulsion); Dirichlet(1), the uniform kmu and
# the half-normal normalisation are constant. Serial and threaded (pair) paths.
# Then the Jacobian: (jacobian TRUE - FALSE) differences between the two points
# = the bounded kmu's log((hi - lo) s (1 - s)) + kappa's log kappa (same kpi).
cat("\n# identify repel: program == brute force (likelihood, pull, repulsion, Jacobian)\n")
gdir_lp <- file.path(tempdir(), "kl_grp_lp"); dir.create(gdir_lp, showWarnings = FALSE)
writeLines(c("iso3,group", "C01,g_east", "C02,g_east", "C03,g_south"), file.path(gdir_lp, "kl_lp_toy.csv"))
optr  <- kl_synthetic_opt(lag_min = -3L, lag_max = 3L,
                          kernel = list(classes = 3L, form = "gaussian", identify = "repel", ksd_fixed = 1.5,
                                        prior_grouping = list(source = "kl_lp_toy", membership = TRUE)))
kcfgr <- resolve_kernel_cfg(optr)
mr    <- kl_sort_by_country(assemble_kl_model(
  optr, panel, list(re_uncorrelated = TRUE), .occ_baseline_block("cyclic"), kcfgr))
pgr   <- kl_prior_partition(kcfgr, levels(droplevels(mr$data$country)), groupings_dir = gdir_lp)
sdr   <- build_kl_latent_standata(mr, pgr, NULL)
sdr_t <- build_kl_latent_standata(mr, pgr, NULL, threads = 2L)
expect(pgr$K == 3L && sdr$kl_rep_G == 2L && sdr$kl_use_pg == 1L && sdr$kl_mu_ordered == 0L &&
       sdr$kl_rs == 0L && sdr_t$kl_rs == 1L && sdr_t$kl_pairs == 1L,
       "repel: K 3, 2 registry groups pulled, free centres; serial and threaded pair paths")
mod_r <- compile(as.character(patch_kl_latent_class_stancode(
  brms::stancode(mr$formula, data = mr$data, family = mr$family, prior = mr$priors), mr$meta, repel = TRUE)),
  "occ_kl_repel")
f_r <- fit1(mod_r, sdr); f_rt <- fit1(mod_r, sdr_t)
KBr <- as.matrix(mr$data[, mr$meta$cb_col_names]); NLr <- length(mr$meta$lags); Jr <- mr$meta$df
ctr <- as.integer(droplevels(mr$data$country)); Cr <- max(ctr)
lse <- function(v) max(v) + log(sum(exp(v - max(v))))
brute_r <- function(kmu, kbeta, kpi, kap, Z) {
  lin <- as.numeric(sdr$X_lin %*% as.numeric(b_lin)) + as.numeric(sdr$offsets_lin)
  W <- kl_weights(kmu, rep(1.5, 3), mr$meta$lags)
  LL <- sapply(1:3, function(k) {
    Bw <- Reduce(`+`, lapply(seq_len(NLr), function(l) W[k, l] * KBr[, (l - 1) * Jr + seq_len(Jr)]))
    rowsum(dpois(sdr$Y, exp(lin + as.numeric(Bw %*% kbeta[, k])), log = TRUE), ctr)[, 1]
  })
  lg  <- matrix(log(kpi), Cr, 3, byrow = TRUE) + kap * (sdr$kl_rep_A %*% Z)
  lsm <- lg - apply(lg, 1, lse)
  ll  <- sum(vapply(seq_len(Cr), function(c) lse(lsm[c, ] + LL[c, ]), 0))
  E   <- sdr$kl_rep_V %*% kbeta
  rep <- 0
  for (j in 1:2) for (k in (j + 1):3)
    rep <- rep + log(1 - exp(-0.5 * (((kmu[j] - kmu[k]) / 3)^2 + sum((E[, j] - E[, k])^2) / (4 * 0.5^2))))
  ll + sum(dnorm(kbeta, 0, kcfgr$prior_beta_sd, log = TRUE)) + dnorm(kap, 0, 2, log = TRUE) +
    sum(dnorm(Z, log = TRUE)) + rep
}
pr_ <- function(kmu, kbeta, kpi, kap, Z, f = f_r, jac = FALSE) {
  p <- list(b_lin = b_lin, zs_lin_1_1 = rep(0, sdr$knots_lin_1[1]), sds_lin_1 = array(rep(0.5, sdr$nb_lin_1)),
            kmu_ord = numeric(0), kmu_free = array(kmu), ksd = numeric(0), kbeta_free = kbeta, kpi = kpi,
            kc = numeric(0), keta = matrix(numeric(0), Jr - 1, 0), kz = matrix(numeric(0), 0, 3),
            khs_tau = numeric(0), khs_lambda = matrix(numeric(0), 0, 3), khs_c2 = numeric(0),
            kre_zc = matrix(numeric(0), 0, Cr), kre_za = matrix(numeric(0), 0, sdr$kl_NA),
            kre_sdc = numeric(0), kre_sda = numeric(0),
            kkappa_free = array(kap), kgrp = Z, kphi = matrix(numeric(0), 0, 3), ktau = numeric(0), krho = numeric(0))
  for (g in 1:2) {
    M <- sdr[[sprintf("M_%d", g)]]; N <- sdr[[sprintf("N_%d", g)]]
    p[[sprintf("sd_%d", g)]] <- array(rep(0.5, M)); p[[sprintf("z_%d", g)]] <- matrix(0, M, N)
  }
  f$log_prob(f$unconstrain_variables(p), jacobian = jac)
}
pt <- function() list(kmu = runif(3, -2.8, 2.8), kbeta = matrix(rnorm(3 * Jr, 0, 0.6), Jr, 3),
                      kpi = { x <- runif(3); x / sum(x) }, kap = runif(1, 0.2, 2.5), Z = matrix(rnorm(6), 2, 3))
set.seed(11)
for (rep in 1:3) {
  a <- pt(); b <- pt()
  d_brute <- do.call(brute_r, b) - do.call(brute_r, a)
  for (path in c("serial", "threaded")) {
    f <- if (path == "serial") f_r else f_rt
    d_stan <- do.call(pr_, c(b, list(f = f))) - do.call(pr_, c(a, list(f = f)))
    cat(sprintf("  repel %s point %d: program %.8f  brute force %.8f  diff %.2e\n", path, rep, d_stan, d_brute,
                d_stan - d_brute))
    expect(abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
           sprintf("repel %s point %d: log-density difference equals the brute force", path, rep))
  }
}
# repulsion alone: two classes brought within 0.01 lag and the same response
# -> the program's log density falls by the brute-force repulsion difference
a <- pt(); b <- a; b$kmu[2] <- b$kmu[1] + 0.01; b$kbeta[, 2] <- b$kbeta[, 1] + 1e-3
d_stan <- pr_(b$kmu, b$kbeta, b$kpi, b$kap, b$Z) - pr_(a$kmu, a$kbeta, a$kpi, a$kap, a$Z)
d_brute <- do.call(brute_r, b) - do.call(brute_r, a)
expect(d_stan < -5 && abs(d_stan - d_brute) < 1e-6 * max(1, abs(d_brute)),
       sprintf("repel: near-duplicate classes penalised (%.1f log units), as the brute force", d_stan))
lo_r <- kcfgr$lag_min; hi_r <- kcfgr$lag_max
jac_r <- function(q) sum(log(hi_r - lo_r) + log((q$kmu - lo_r) / (hi_r - lo_r)) + log(1 - (q$kmu - lo_r) / (hi_r - lo_r))) +
  log(q$kap)
for (rep in 1:2) {
  a <- pt(); b <- pt(); b$kpi <- a$kpi
  dj_stan <- (do.call(pr_, c(b, list(jac = TRUE))) - do.call(pr_, c(b, list(jac = FALSE)))) -
             (do.call(pr_, c(a, list(jac = TRUE))) - do.call(pr_, c(a, list(jac = FALSE))))
  dj <- jac_r(b) - jac_r(a)
  cat(sprintf("  repel Jacobian point %d: program %.8f  analytic %.8f  diff %.2e\n", rep, dj_stan, dj, dj_stan - dj))
  expect(abs(dj_stan - dj) < 1e-6 * max(1, abs(dj)),
         sprintf("repel point %d: Jacobian of the bounded centres and kappa as analytic", rep))
}

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1)
cat("PASS  test_kl_lp_equivalence\n")
