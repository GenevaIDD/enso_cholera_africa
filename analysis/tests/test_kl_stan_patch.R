# test_kl_stan_patch.R — HARD: brms_kl building blocks, no sampling.
#
# Usage:
#   Rscript analysis/tests/test_kl_stan_patch.R
#   Rscript analysis/tests/test_kl_stan_patch.R --capture
#
# Details:
#   The lag basis, the stage-1 program, the patched latent-class program
#   (byte identity with tests/golden/kl_latent_program.stan; regenerate with
#   --capture and review the diff), the helpers and every kernel option,
#   artefact names and config refusals. The syntax check needs CmdStan's
#   stanc (skipped with a note if absent); everything else is pure R.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})
for (m in c(
  "base",
  "config",
  "formula",
  "priors",
  "paths",
  "enso",
  "dlnm",
  "dlnm_kernel",
  "loo_lags"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_synth.R")
PASS <- 0L
FAIL <- 0L
CAPTURE <- "--capture" %in% commandArgs(TRUE)
GOLDEN <- "analysis/tests/golden/kl_latent_program.stan"

# [1] multi-lag basis ----------------------------------------------------------
# [1] dlnm `lags:` basis: shared pooled basis replayed per lag; design replay
cat("\n[1] dlnm lags: basis\n")
set.seed(3)
n <- 120
d <- data.frame(
  enso_lead2 = rnorm(n),
  enso_lead1 = rnorm(n),
  enso_lag0 = rnorm(n),
  enso_lag1 = c(NA, rnorm(n - 1)),
  enso_lag2 = rnorm(n)
)
cfg <- list(
  enable = TRUE,
  defaults = list(argvar = list(fun = "ns", df = 3L)),
  variables = list(enso = list(lags = list(-2, 2)))
)
cbb <- suppressMessages(build_cb_block(d, cfg))
mt <- cbb$cb_meta$enso
expect(
  identical(mt$lags, -2:2) && mt$df == 3L && length(mt$cb_col_names) == 15L,
  "lags [-2, 2] -> 5 lags x 3 basis columns"
)
expect(
  identical(
    mt$cb_col_names[1:4],
    c("enso_cb1_lm2", "enso_cb2_lm2", "enso_cb3_lm2", "enso_cb1_lm1")
  ),
  "columns lag-major, lead token m<k>"
)
X <- sapply(mt$lags, function(l) {
  x <- d[[enso_lag_col_name(l)]]
  x[is.na(x)] <- 0
  x
})
pooled <- unclass(dlnm::onebasis(as.vector(X), fun = "ns", df = 3L))
B <- as.matrix(cbb$dat[, mt$cb_col_names])
expect(
  all(vapply(
    seq_along(mt$lags),
    function(i) {
      isTRUE(all.equal(
        unname(B[, (i - 1) * 3 + 1:3]),
        unname(pooled[(i - 1) * n + 1:n, ]),
        check.attributes = FALSE
      ))
    },
    logical(1)
  )),
  "each lag's columns = the pooled basis rows (one shared basis)"
)
expect(
  isTRUE(all.equal(
    as.matrix(build_cb_design(d, cbb$cb_meta)),
    B,
    check.attributes = FALSE
  )),
  "build_cb_design() replays the fit-time basis"
)
expect(
  inherits(
    try(
      suppressMessages(build_cb_block(
        d,
        modifyList(
          cfg,
          list(
            variables = list(enso = list(lags = list(-2, 2), single_lag = 1))
          )
        )
      )),
      silent = TRUE
    ),
    "try-error"
  ),
  "lags and single_lag together are refused"
)

# [2] stage-1 program ------------------------------------------------------------
# [2] kernel expression and the brms-generated stage-1 program: every patch
#     landmark exactly once, C_1.. = the basis columns in stored order
cat("\n[2] kernel formula and brms program\n")
panel <- make_synthetic_kl_panel(
  n_country = 4L,
  adm1_per_country = 2L,
  n_year = 3L,
  lag_min = -3L,
  lag_max = 3L,
  K = 2L,
  mu_true = c(-1, 2),
  sigma_true = c(1, 1),
  beta_true = c(0.5, -0.5)
)
opt <- kl_synthetic_opt(
  lag_min = -3L,
  lag_max = 3L,
  kernel = list(classes = 2L)
)
kcfg <- resolve_kernel_cfg(opt)
m <- kl_sort_by_country(suppressMessages(assemble_kl_model(
  opt,
  panel,
  list(re_uncorrelated = TRUE),
  .occ_baseline_block("cyclic"),
  kcfg
)))
J <- m$meta$df
NL <- length(m$meta$lags)
ex <- build_kl_expr(m$meta)
expect(
  all(vapply(
    m$meta$cb_col_names,
    function(cn) {
      lengths(regmatches(
        ex,
        gregexpr(paste0("\\b", cn, "\\b"), ex, perl = TRUE)
      )) ==
        1L
    },
    logical(1)
  )),
  "kernel expression mentions every basis column exactly once"
)
expect(
  all(c("drop_region") %in% names(kcfg)) && isTRUE(kcfg$drop_region),
  "K >= 2 drops the region level by default"
)
expect(
  !any(grepl("AFRO_region", m$terms$re_terms)),
  "region level absent from lin REs"
)
sc <- suppressWarnings(brms::stancode(
  m$formula,
  data = m$data,
  family = m$family,
  prior = m$priors
))
L <- strsplit(sc, "\n", fixed = TRUE)[[1]]
lm <- kl_stan_landmarks(J)
for (nm in names(lm)) {
  hits <- tryCatch(.kl_find1(L, lm[[nm]], nm), error = function(e) NA_integer_)
  expect(!is.na(hits), sprintf("landmark '%s' found exactly once", nm))
}
sd1 <- suppressWarnings(unclass(brms::standata(
  m$formula,
  data = m$data,
  family = m$family,
  prior = m$priors
)))
expect(
  sum(grepl("^C_[0-9]+$", names(sd1))) == J * NL,
  "brms passes J x NL kernel covariates"
)
expect(
  all(vapply(
    seq_len(J * NL),
    function(i) {
      isTRUE(all.equal(
        as.numeric(sd1[[paste0("C_", i)]]),
        m$data[[m$meta$cb_col_names[i]]]
      ))
    },
    logical(1)
  )),
  "C_i == the i-th stored basis column (lag-major order)"
)

# [3] patched program --------------------------------------------------------------
# [3] patched program: stanc syntax, data identifiers == the R data list,
#     byte identity with the pinned golden (tests/golden/kl_latent_program.stan;
#     regenerate with --capture and review the diff)
cat("\n[3] patched latent-class program\n")
code <- patch_kl_latent_class_stancode(sc, m$meta)
expect(
  identical(
    attr(code, "brms_version"),
    as.character(utils::packageVersion("brms"))
  ),
  "header brms version == installed brms"
)
lv <- levels(droplevels(m$data$country))
pg <- kl_prior_partition(kcfg, lv)
adj <- kl_adjacency_data(attr(panel, "adjacency"), lv)
kcfg_sp <- kcfg
kcfg_sp$spatial$enable <- TRUE
sdk <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kcfg_sp)),
  pg,
  adj
))
Lp <- strsplit(code, "\n", fixed = TRUE)[[1]]
data_block <- Lp[
  (which(Lp == "data {") + 1):(which(Lp == "transformed data {") - 2)
]
decl <- data_block[!grepl("^\\s*//", data_block) & grepl(";", data_block)]
ids <- sub(".*?([A-Za-z_][A-Za-z0-9_]*);.*$", "\\1", decl, perl = TRUE)
expect(
  all(ids %in% names(sdk)),
  sprintf(
    "every declared data identifier is supplied (missing: %s)",
    paste(setdiff(ids, names(sdk)), collapse = ",")
  )
)
expect(
  all(setdiff(names(sdk), ids) %in% names(sd1)),
  sprintf(
    "undeclared entries are brms's own unused data only (%s)",
    paste(setdiff(names(sdk), ids), collapse = ",")
  )
)
expect(
  sdk$kl_use_sp == 1L && sdk$kl_n_edges == 3L && sdk$kl_n_comp == 1L,
  "ICAR data: 4 countries on a line -> 3 edges, 1 component"
)
if (!inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")) {
  f <- tempfile(fileext = ".stan")
  writeLines(code, f)
  ok <- cmdstanr::cmdstan_model(f, compile = FALSE)$check_syntax(quiet = TRUE)
  expect(isTRUE(ok), "patched program passes stanc")
} else {
  cat("  NOTE  CmdStan absent: stanc syntax check skipped\n")
}
if (CAPTURE) {
  writeLines(code, GOLDEN)
  cat("  captured ", GOLDEN, "\n")
} else {
  expect(
    file.exists(GOLDEN) &&
      identical(readLines(GOLDEN), strsplit(code, "\n", fixed = TRUE)[[1]]),
    "patched program == pinned golden (regenerate with --capture, review the diff)"
  )
}

# [4] helpers ----------------------------------------------------------------
# [4] adjacency components, prior partition, lag weights; then [4b]-[4k]:
#     identification; class lag-basis form; discrete centres (ordered
#     combinations, prior, tie relabelling); pinned arrangements (seeds,
#     fixed classes, orders); identify repel (data, program = default + repel
#     fragments, inits, artefact token, per-draw relabelling and registry
#     naming); kernel.lag_re (country lag offsets: config, refusals, program
#     variant = default with stacked kernel terms, data on every path, inits,
#     token)
cat("\n[4] adjacency, prior partition, weights\n")
toy <- data.frame(
  iso3_a = c("A", "B", "D", "Z"),
  iso3_b = c("B", "C", "E", "A")
)
ad <- kl_adjacency_data(toy, c("A", "B", "C", "D", "E", "F"))
expect(
  ad$n_edges == 3L &&
    ad$n_comp == 3L &&
    identical(as.integer(ad$comp), c(1L, 1L, 1L, 2L, 2L, 3L)),
  "edges to unfitted countries dropped; components A-B-C | D-E | F"
)
expect(
  identical(ad$isolated, "F") && identical(as.numeric(ad$comp_n), c(3, 2, 1)),
  "isolated country and component sizes"
)
expect(
  identical(as.integer(ad$n_par), c(0L, 1L, 1L, 0L, 1L, 0L)) &&
    identical(as.integer(ad$par_nodes), c(1L, 2L, 4L)),
  "DAGAR parents = lower-indexed neighbours (A<-, B<-A, C<-B, D<-, E<-D, F<-)"
)
gdir <- file.path(tempdir(), "kl_grp")
dir.create(gdir, showWarnings = FALSE)
writeLines(
  c("iso3,group", "C01,g_east", "C02,g_east", "C03,g_south", "XXX,g_west"),
  file.path(gdir, "kl_toy.csv")
)
k2 <- kcfg
k2$prior_grouping <- list(source = "kl_toy")
k2$classes <- 3L
pp <- kl_prior_partition(
  k2,
  c("C01", "C02", "C03", "C04"),
  groupings_dir = gdir
)
expect(
  pp$use &&
    pp$K == 3L &&
    identical(pp$class_labels, c("g_east", "g_south", "klc3")),
  "prior partition: registry groups present in data, then unanchored extra class"
)
expect(
  identical(unname(rowSums(pp$A)), c(1, 1, 1, 0)) &&
    pp$A["C03", "g_south"] == 1,
  "one-hot rows; country absent from the registry unanchored"
)
expect(!pp$mu_normal, "no deploy_config / prior_mu -> uniform kmu prior")
# Contiguity mix: C01 borders C02 (east), C03 (south) and XXX (g_west, not a
# class: no fitted member) -> shares east 1/2, south 1/2; C02 borders only C01
# (interior); C04 is outside the registry.
writeLines(
  c("iso3_a,iso3_b", "C01,C02", "C01,C03", "C01,XXX", "C03,C04"),
  file.path(gdir, "kl_toy_adj.csv")
)
k3 <- k2
k3$spatial$adjacency <- "kl_toy_adj"
k3$prior_grouping$contiguity <- list(enable = TRUE, weight = 0)
expect(
  identical(
    suppressMessages(kl_prior_partition(
      k3,
      c("C01", "C02", "C03", "C04"),
      groupings_dir = gdir
    ))$A,
    pp$A
  ),
  "contiguity mix, weight 0: A identical to the one-hot partition"
)
k3$prior_grouping$contiguity$weight <- 0.5
pc <- suppressMessages(kl_prior_partition(
  k3,
  c("C01", "C02", "C03", "C04"),
  groupings_dir = gdir
))
expect(
  isTRUE(all.equal(unname(pc$A["C01", ]), c(0.75, 0.25, 0))) &&
    identical(pc$A["C02", ], pp$A["C02", ]) &&
    identical(pc$A["C04", ], pp$A["C04", ]),
  "contiguity mix, weight 0.5: border row (0.75, 0.25, 0); interior and unregistered rows unchanged"
)
expect(
  isTRUE(all.equal(unname(pc$A["C03", ]), c(0.5, 0.5, 0))),
  "contiguity mix: an unregistered neighbour (C04) is ignored (C03's only counted neighbour is C01, east)"
)
expect(
  isTRUE(all.equal(unname(rowSums(pc$A)), c(1, 1, 1, 0))),
  "contiguity mix rows stay on the simplex"
)
k4 <- k3
k4$identify <- "seed"
k4$prior_grouping$membership <- NULL
expect(
  inherits(
    tryCatch(
      kl_prior_partition(k4, c("C01", "C02"), groupings_dir = gdir),
      error = identity
    ),
    "error"
  ),
  "contiguity under identify seed without membership: refused (A would be unused)"
)
expect(
  grepl(
    "-mem-cmix0p5",
    .kl_token(list(
      classes = 1L,
      prior_grouping = list(
        source = "enso_grouping",
        membership = TRUE,
        contiguity = list(enable = TRUE, weight = 0.5)
      )
    ))
  ) &&
    !grepl(
      "cmix",
      .kl_token(list(
        classes = 1L,
        prior_grouping = list(
          source = "enso_grouping",
          contiguity = list(enable = FALSE, weight = 0.5)
        )
      ))
    ),
  "artefact token -cmix<w> only when contiguity is enabled"
)
W <- kl_weights(c(-1, 4), c(0.5, 3), -3:3)
expect(
  isTRUE(all.equal(rowSums(W), c(1, 1))) && which.max(W[1, ]) == 3L,
  "kernel weights normalised, peak at the centre"
)
writeLines(
  c("# test adjacency", "iso3_a,iso3_b", "B,A", "A,B", "C,C", "C,D"),
  file.path(gdir, "kl_adj.csv")
)
la <- load_adjacency("kl_adj", gdir)
expect(
  identical(la$iso3_a, c("A", "C")) && identical(la$iso3_b, c("B", "D")),
  "load_adjacency: comment skipped, edges undirected + unique, self-loop dropped"
)

cat("\n[4b] identification (kernel.identify)\n")
ev <- kl_erf_direction(m$meta, kcfg)
expect(
  isTRUE(all.equal(crossprod(ev$Q), diag(ncol(ev$Q)))) &&
    max(abs(crossprod(ev$Q, ev$v))) < 1e-10,
  "erf: Q orthonormal and orthogonal to the contrast v"
)
kc0 <- c(-0.4, 1.1)
eta0 <- matrix(c(0.3, -0.2, 0.5, 0.1), J - 1, 2)
kb0 <- sapply(1:2, function(k) kc0[k] * ev$v / sum(ev$v^2) + ev$Q %*% eta0[, k])
expect(
  isTRUE(all.equal(as.numeric(crossprod(ev$v, kb0)), kc0)),
  "erf: v' kbeta reproduces the ordered contrast kc"
)
idf <- function(id, pgsrc = NULL) {
  k <- kcfg
  k$identify <- id
  k$prior_grouping$source <- pgsrc
  kl_identify(k)
}
expect(
  idf("auto") == "kmu" &&
    idf("auto", "x") == "anchor" &&
    idf("erf", "x") == "erf",
  "identify: auto -> kmu / anchor; erf kept"
)
expect(
  inherits(try(idf("anchor"), silent = TRUE), "try-error"),
  "identify anchor without a prior partition refused"
)
kerf <- kcfg
kerf$identify <- "erf"
sde <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kerf)),
  pg,
  NULL
))
expect(
  sde$kl_erf_ordered == 1L && sde$kl_mu_ordered == 0L && sde$kl_use_pg == 0L,
  "identify erf: ERF ordered, kmu free, no anchoring term"
)
# naming of rank-identified classes from a registry
A3 <- matrix(
  c(1, 1, 0, 0, 0, 1),
  3,
  2,
  dimnames = list(c("X", "Y", "Z"), c("g1", "g2"))
)
lr <- log(rbind(c(0.1, 0.9), c(0.2, 0.8), c(0.95, 0.05))) # rows = countries X, Y, Z; cols = classes
fake <- posterior::as_draws_df(matrix(
  rep(as.vector(t(lr)), each = 4),
  4,
  dimnames = list(
    NULL,
    sprintf("kl_log_resp[%d,%d]", rep(1:2, 3), rep(1:3, each = 2))
  )
))
expect(
  identical(kl_name_classes(fake, A3, c("g1", "g2"), 2L), c("g2", "g1")),
  "naming: class whose members are g2's countries is called g2"
)

cat("\n[4c] class lag-basis form (kernel.form dlnm) and post-hoc relabelling\n")
Lb <- kl_lag_basis(-11:11, 7L)
expect(
  identical(dim(Lb), c(23L, 7L)) &&
    isTRUE(all.equal(unname(rowSums(Lb)), rep(1, 23))),
  "lag basis: 23 lags x 7 B-splines, rows sum to 1"
)
kd <- kcfg
kd$form <- "dlnm"
kd$identify <- "auto"
expect(kl_identify(kd) == "relabel", "form dlnm: identify auto -> relabel")
sdd2 <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kd)),
  pg,
  NULL
))
expect(
  sdd2$kl_dlnm == 1L &&
    sdd2$kl_hs == 1L &&
    sdd2$kl_field_cols == pg$K &&
    sdd2$kl_mu_ordered == 0L &&
    sdd2$kl_erf_ordered == 0L,
  "form dlnm data: horseshoe on, symmetric field (K columns), no ordering"
)
expect(
  all(ids %in% names(sdd2)),
  "form dlnm: every declared data identifier supplied"
)
kb2 <- kd
kb2$engine <- "brms"
expect(
  inherits(try(kl_engine(kb2), silent = TRUE), "try-error"),
  "form dlnm with engine brms refused"
)
kx <- kd
kx$identify <- "erf"
expect(
  inherits(try(kl_identify(kx), silent = TRUE), "try-error"),
  "form dlnm with identify erf refused"
)
# relabelling: two classes stored in swapped order in half of the draws
mk <- function(k1, k2) {
  c(
    `kcon[1]` = k1,
    `kcon[2]` = k2,
    `kpi[1]` = 0.3,
    `kpi[2]` = 0.7,
    `kbeta[1,1]` = k1,
    `kbeta[1,2]` = k2,
    `kl_log_resp[1,1]` = log(0.9),
    `kl_log_resp[2,1]` = log(0.1),
    `klogit[1,1]` = 1,
    `klogit[1,2]` = 2
  )
}
toy <- posterior::as_draws_df(as.data.frame(rbind(
  mk(-1, 2),
  mk(2, -1),
  mk(-1.2, 1.8),
  mk(1.9, -0.8)
)))
rl <- kl_relabel_draws(toy, 2L)
rd <- as.data.frame(rl$draws)
expect(
  all(rd$`kcon[1]` < rd$`kcon[2]`),
  "relabel: kcon increasing in every draw"
)
expect(
  identical(rd$`kbeta[1,1]`, rd$`kcon[1]`) &&
    identical(rd$`kpi[1]`, c(0.3, 0.7, 0.3, 0.7)),
  "relabel: kbeta and kpi permuted with their class"
)
expect(
  isTRUE(all.equal(exp(rd$`kl_log_resp[1,1]`), c(0.9, 0.1, 0.9, 0.1))) &&
    identical(rd$`klogit[1,1]`, c(1, 2, 1, 2)),
  "relabel: responsibilities (class first) and logits (class second) permuted"
)
expect(
  nrow(rl$table) == 2L && all(rl$table$share == 0.5),
  "relabel: permutation table"
)

cat(
  "\n[4d] discrete centres (kernel.form discrete): ordered combinations, prior, relabel\n"
)
cb <- kl_disc_combos(7L, 2L)
expect(
  nrow(cb) == 28L &&
    all(cb[, 1] <= cb[, 2]) &&
    !anyDuplicated(apply(cb, 1, paste, collapse = ",")),
  "7 centres, K = 2: 28 distinct non-decreasing pairs"
)
expect(
  nrow(kl_disc_combos(23L, 4L)) == choose(26, 4) &&
    nrow(kl_disc_combos(23L, 1L)) == 23L,
  "23 centres: K = 4 -> choose(26, 4) = 14950; K = 1 -> 23"
)
lab <- as.matrix(expand.grid(1:7, 1:7, 1:7))
mult <- table(apply(t(apply(lab, 1, sort)), 1, paste, collapse = ","))
cb3 <- kl_disc_combos(7L, 3L)
kd3 <- kcfg
kd3$form <- "discrete"
kd3$discrete <- list(widths = c(0, 1.5))
ds3 <- kl_disc_spec(kd3, -3:3, 3L)
expect(
  isTRUE(all.equal(
    exp(ds3$combo_lp),
    as.numeric(mult[apply(cb3, 1, paste, collapse = ",")]),
    check.attributes = FALSE
  )),
  "exchangeable prior exp(lp) = number of labelled tuples sorting to each combination"
)
expect(
  isTRUE(all.equal(unname(colSums(ds3$Wq)), rep(1, 14))) &&
    all(ds3$Wq[, 1:7] %in% c(0, 1)) &&
    identical(unname(apply(ds3$Wq[, 1:7], 2, which.max)), 1:7),
  "width 0 = single lag at the centre; every kernel column sums to 1"
)
kd3s <- kd3
kd3s$discrete$centre_step <- 2L
kd3s$discrete$combo_prior <- "uniform"
ds3s <- kl_disc_spec(kd3s, -3:3, 2L)
expect(
  identical(ds3s$centres, c(-3, -1, 1, 3)) && all(ds3s$combo_lp == 0),
  "centre_step 2 -> every other lag; combo_prior uniform -> 0"
)
dsn <- kl_disc_spec(
  kd3,
  -3:3,
  2L,
  mu_normal = TRUE,
  mu_mean = c(-2, 2),
  mu_sd = 0.5
)
expect(
  identical(
    as.numeric(ds3$centres[dsn$combos[which.max(dsn$combo_lp), ]]),
    c(-2, 2)
  ),
  "position priors: the most probable combination sits at the (sorted) prior means"
)
expect(
  inherits(
    try(
      kl_disc_spec(
        modifyList(kd3, list(discrete = list(widths = 1:6))),
        -11:11,
        5L
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  "oversized grid refused (S * M cap)"
)
expect(
  isTRUE(all.equal(kl_weights(1, 0, -3:3)[1, ], c(0, 0, 0, 0, 1, 0, 0))),
  "kl_weights: width 0 -> single lag"
)
expect(
  kl_identify(kd3) == "lagorder",
  "form discrete: identify auto -> lagorder"
)
kx3 <- kd3
kx3$identify <- "kmu"
expect(
  inherits(try(kl_identify(kx3), silent = TRUE), "try-error"),
  "form discrete with identify kmu refused"
)
kg <- kcfg
kg$identify <- "lagorder"
expect(
  inherits(try(kl_identify(kg), silent = TRUE), "try-error"),
  "identify lagorder without form discrete refused"
)
kb3 <- kd3
kb3$engine <- "brms"
expect(
  inherits(try(kl_engine(kb3), silent = TRUE), "try-error"),
  "form discrete with engine brms refused"
)
sdx <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kd3)),
  pg,
  NULL
))
expect(
  sdx$kl_disc == 1L &&
    sdx$kl_mu_ordered == 0L &&
    sdx$kl_erf_ordered == 0L &&
    sdx$kl_use_pg == 0L &&
    sdx$kl_field_cols == pg$K &&
    sdx$kl_M == 28L &&
    sdx$kl_S == 2L &&
    identical(dim(sdx$kl_Wq), c(NL, 14L)),
  "form discrete data: no kmu/ERF ordering, no anchoring, symmetric field, 28 x 2 terms"
)
expect(
  all(ids %in% names(sdx)),
  "form discrete: every declared data identifier supplied"
)
expect(
  sdk$kl_disc == 0L && sdk$kl_M == 1L,
  "other forms: discrete switched off (dummy grid)"
)
# relabelling within ties only
mkd <- function(l1, l2, c1, c2) {
  c(
    `klag[1]` = l1,
    `klag[2]` = l2,
    `kcon[1]` = c1,
    `kcon[2]` = c2,
    `kpi[1]` = 0.3,
    `kpi[2]` = 0.7,
    `kbeta[1,1]` = c1,
    `kbeta[1,2]` = c2,
    `kl_lag_prob[1,1]` = 0.9,
    `kl_lag_prob[2,1]` = 0.1,
    `kwidth[1]` = 2
  )
}
toyd <- posterior::as_draws_df(as.data.frame(rbind(
  mkd(2, 2, 1, -1),
  mkd(-3, 5, 1, -1),
  mkd(2, 2, -1, 1)
)))
rd2 <- as.data.frame(kl_relabel_draws(toyd, 2L, lag = "klag")$draws)
expect(
  identical(rd2$`kcon[1]`, c(-1, 1, -1)) &&
    identical(rd2$`kpi[1]`, c(0.7, 0.3, 0.3)) &&
    identical(rd2$`kl_lag_prob[1,1]`, c(0.1, 0.9, 0.9)),
  "lagorder relabel: tied classes ordered by kcon, untied classes left in lag order"
)
rn <- as.data.frame(kl_disc_rename_draws(toyd, 2L))
expect(
  identical(rn$`kmu[2]`, c(2, 5, 2)) &&
    identical(rn$`ksd[1]`, c(2, 2, 2)) &&
    !"klag[1]" %in% names(rn),
  "rename: klag -> kmu, kwidth -> ksd per class"
)

cat(
  "\n[4e] pinned arrangements (kernel.pin): seeds, fixed classes, candidate orders\n"
)
dseed <- data.frame(
  country = factor(c("A", "A", "B", "C", "C", "C", "D")),
  y = c(1, 0, 3, 1, 1, 0, 5)
)
regs <- data.frame(
  iso3 = c("A", "B", "C", "D"),
  group = c("g1", "g1", "g2", "g2")
)
sd0 <- kl_pin_seeds(dseed, regs, c("g1", "g2"), n = 1L)
expect(
  identical(sd0, list(g1 = "B", g2 = "D")),
  "seeds: each group's country with the most events"
)
expect(
  identical(
    kl_pin_seeds(dseed, regs, c("g1", "g2"), n = 2L, seeds = list(g2 = "C"))$g2,
    "C"
  ) &&
    identical(kl_pin_seeds(dseed, regs, c("g1", "g2"), n = 2L)$g1, c("B", "A")),
  "seeds: explicit list overrides; n = 2 by events"
)
expect(
  identical(
    kl_pin_fix(c("A", "B", "C", "D"), sd0, c("g2", "g1")),
    c(0L, 2L, 0L, 1L)
  ),
  "fix: seeds of sigma[k] -> class k, others free"
)
fake_pilot <- function(g, probs) list(group = g, lags = -3:3, lag_prob = probs)
pa <- list(
  fake_pilot("ga", c(0, 1, 0, 0, 0, 0, 0)),
  fake_pilot("gb", c(0, 0, 0, 0, 0, 0.5, 0.5)),
  fake_pilot("gc", c(0, 0, 0, 0, 0, 0.5, 0.5))
)
ar <- kl_pin_arrangements(pa, "pilot", cover = 0.95)
expect(
  nrow(ar) == 2L &&
    setequal(ar$arrangement, c("ga>gb>gc", "ga>gc>gb")) &&
    isTRUE(all.equal(sum(ar$prob), 1, tolerance = 0.02)),
  "orders: ga surely shortest; gb/gc overlapping -> both orders kept"
)
ara <- kl_pin_arrangements(pa, "all")
expect(
  nrow(ara) == 6L && all(lengths(ara$sigma) == 3L),
  "orders: mode all -> 3! permutations"
)
kp <- kd3
kp$pin <- list(enable = TRUE)
expect(
  inherits(try(kl_pin_cfg(kp), silent = TRUE), "try-error"),
  "pin without prior_grouping refused"
)
kp$prior_grouping$source <- "x"
expect(
  isTRUE(kl_pin_cfg(kp)$enable),
  "pin with form discrete + prior_grouping accepted"
)
expect(
  is.null(
    kl_pin_cfg(modifyList(
      kp,
      list(pin = list(enable = TRUE, seeds_n = 2L))
    ))$seeds
  ),
  "pin: absent seeds stay NULL (not partially matched to seeds_n)"
)
kpg <- kp
kpg$form <- "gaussian"
expect(
  inherits(try(kl_pin_cfg(kpg), silent = TRUE), "try-error"),
  "pin with form gaussian refused"
)
# order_source sweep: lag weights from a sweep's LOO ranking caches
sw_out <- file.path(tempdir(), "kl_sweep_out")
mk_rank <- function(g, el) {
  f <- make_loo_compare_run_file(
    list(outdir = sw_out),
    "loo_compare_ranking",
    batch = "kl_fake_sweep",
    group = g,
    date = "2026-10-12"
  )
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(ranking = data.frame(single_lag = -3:3, elpd_loo = el)), f)
}
mk_rank("ga", c(-50, -50, -50, -10, -50, -50, -50)) # ga: lag 0
mk_rank("gb", c(-50, -50, -50, -50, -50, -10, -10)) # gb: lag 2 or 3
so <- kl_pin_sweep_orders("kl_fake_sweep", c("ga", "gb"), outdir = sw_out)
expect(
  isTRUE(all.equal(
    unname(so$ga$lag_prob[so$ga$lags == 0]),
    1,
    tolerance = 1e-12
  )) &&
    isTRUE(all.equal(
      unname(so$gb$lag_prob[so$gb$lags %in% 2:3]),
      c(0.5, 0.5),
      tolerance = 1e-12
    )),
  "sweep: lag weights exp(elpd - max) per group"
)
expect(
  identical(kl_pin_arrangements(so, "pilot")$arrangement, "ga>gb"),
  "sweep: weights feed the same order selection"
)
expect(
  inherits(
    try(
      kl_pin_sweep_orders("kl_fake_sweep", "gz", outdir = sw_out),
      silent = TRUE
    ),
    "try-error"
  ),
  "sweep: missing group ranking refused"
)
expect(
  identical(kl_pin_run_file("a/x_kl_pin_b.rds", 3), "a/x_kl_pin_b_arr03.rds") &&
    identical(
      kl_pin_run_file("a/x.rds", 1:2),
      c("a/x_arr01.rds", "a/x_arr02.rds")
    ),
  "SLURM split: per-order run file = kl_pin path + _arr<i>"
)
kps <- kp
kps$pin$order_source <- "sweep"
expect(
  inherits(try(kl_pin_cfg(kps), silent = TRUE), "try-error"),
  "sweep without sweep_batch or deploy_config refused"
)
kps$pin$sweep_batch <- "kl_fake_sweep"
expect(
  identical(kl_pin_cfg(kps)$order_source, "sweep"),
  "sweep with sweep_batch accepted"
)
sdf <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kd3)),
  pg,
  NULL,
  fix = c(2L, 0L, 0L, 1L)
))
expect(
  identical(as.integer(sdf$kl_fix), c(2L, 0L, 0L, 1L)) &&
    identical(as.integer(sdx$kl_fix), rep(0L, 4L)),
  "data: kl_fix passed through; all free by default"
)

cat("\n[4f] random slopes on the exposure coefficients (kernel.re_slope)\n")
kre <- kd3
kre$re_slope <- list(enable = TRUE, sd_country = 0.2, sd_adm1 = 0.3)
sdr <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kre)),
  pg,
  NULL
))
a1 <- interaction(
  as.character(m$data$country),
  as.character(m$data$ADM1),
  drop = TRUE,
  lex.order = TRUE
)
expect(
  sdr$kl_re == 1L &&
    sdr$kl_NA == nlevels(a1) &&
    identical(as.integer(sdr$kl_adm1), as.integer(a1)) &&
    identical(as.numeric(sdr$kl_re_sd), c(0.2, 0.3)),
  "re_slope data: country:ADM1 units, sd priors passed"
)
expect(sdx$kl_re == 0L, "re_slope off by default")
krg <- kcfg
krg$re_slope <- list(enable = TRUE)
krg$form <- "dlnm"
expect(
  inherits(try(kl_engine(krg), silent = TRUE), "try-error"),
  "re_slope with form dlnm refused"
)

cat(
  "\n[4g] Gaussian + seeds: fixed width, threaded country path, sweep modes\n"
)
kgs <- kcfg
kgs$identify <- "seed"
expect(
  inherits(try(kl_identify(kgs), silent = TRUE), "try-error"),
  "identify seed without prior_grouping refused"
)
kgs$prior_grouping$source <- "x"
expect(
  kl_identify(kgs) == "seed",
  "identify seed accepted with a prior partition"
)
kgs$ksd_fixed <- 1.5
kgs$re_slope <- list(enable = TRUE)
sdg <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kgs)),
  pg,
  NULL,
  threads = 4L
))
expect(
  sdg$kl_rs == 1L &&
    sdg$kl_re == 1L &&
    sdg$kl_ksd_fix == 1L &&
    sdg$kl_ksd_value == 1.5 &&
    sdg$kl_mu_ordered == 0L &&
    sdg$kl_use_pg == 0L &&
    sdg$kl_maxn == max(c(1L, table(m$data$country)[sdg$kl_small])),
  "gaussian + seed data: country path, random slopes, fixed width, unordered centres, no anchoring"
)
expect(
  all(ids %in% names(sdg)),
  "gaussian + seed: every declared data identifier supplied"
)
kga <- kgs
kga$re_slope <- list(enable = TRUE, country = FALSE)
sdga <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kga)),
  pg,
  NULL,
  threads = 4L
))
ina <- make_kl_class_inits(sdga, pg, kga, 1L, strategy = "spread")[[1]]
expect(
  sdga$kl_re == 1L &&
    sdga$kl_re_ctry == 0L &&
    sdg$kl_re_ctry == 1L &&
    is.null(ina$kre_zc) &&
    is.null(ina$kre_sdc) &&
    length(ina$kre_sda) == sdga$kl_J,
  "re_slope country = no: ADM1 deviations only (no country block, none in the inits)"
)
kgc <- kgs
kgc$re_slope <- list(enable = TRUE, adm1 = FALSE)
sdgc <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kgc)),
  pg,
  NULL,
  threads = 4L
))
inc <- make_kl_class_inits(sdgc, pg, kgc, 1L, strategy = "spread")[[1]]
expect(
  sdgc$kl_re == 1L &&
    sdgc$kl_re_adm1 == 0L &&
    sdgc$kl_re_ctry == 1L &&
    sdg$kl_re_adm1 == 1L &&
    is.null(inc$kre_za) &&
    is.null(inc$kre_sda) &&
    length(inc$kre_sdc) == sdgc$kl_J,
  "re_slope adm1 = no: country deviations only (no ADM1 block, none in the inits)"
)
kg0 <- kgs
kg0$re_slope <- list(enable = TRUE, country = FALSE, adm1 = FALSE)
expect(
  inherits(try(kl_engine(kg0), silent = TRUE), "try-error"),
  "re_slope with country and adm1 both off refused"
)
expect(
  grepl(
    "kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA)",
    code,
    fixed = TRUE
  ) &&
    !grepl("diag_pre_multiply(kre_sda, kre_za)", code, fixed = TRUE),
  "program: every ADM1 deviation goes through kl_adm1_dev() (zero when re_slope.adm1 is off)"
)
sdg1 <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kcfg)),
  pg,
  NULL
))
expect(
  sdg1$kl_rs == 0L && sdg1$kl_ksd_fix == 0L,
  "default gaussian: serial path, sampled width"
)
sp8 <- kl_country_split(c(100, 10, 10, 30, 10), 8L)
expect(
  identical(sp8$big, c(1L, 4L)) &&
    identical(sp8$small, c(2L, 3L, 5L)) &&
    identical(kl_country_split(c(100, 10), 1L)$big, integer(0)) &&
    identical(kl_country_split(c(100, 10), 8L, enable = FALSE)$big, integer(0)),
  "split: countries above N / (2 threads) rows split by rows; none unthreaded or disabled"
)
expect(
  sdg$kl_nsmall + sdg$kl_nbig == sdg$kl_C &&
    setequal(c(sdg$kl_small, sdg$kl_big), seq_len(sdg$kl_C)) &&
    sdg1$kl_nbig == 0L &&
    sdg1$kl_nsmall == sdg1$kl_C,
  "split data: every country summed once; unthreaded -> none split"
)
kgf <- kgs
kgf$re_slope <- list(enable = TRUE, country = TRUE, sd_country_fixed = 0.2)
sdgf <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kgf)),
  pg,
  NULL,
  threads = 4L
))
inf <- make_kl_class_inits(sdgf, pg, kgf, 1L, strategy = "spread")[[1]]
expect(
  sdgf$kl_re_ctry == 1L &&
    sdgf$kl_re_sdc_fix == 1L &&
    sdgf$kl_re_sdc_value == 0.2 &&
    !is.null(inf$kre_zc) &&
    is.null(inf$kre_sdc) &&
    sdg$kl_re_sdc_fix == 0L,
  "re_slope sd_country_fixed: country deviations kept, their sd fixed (no kre_sdc parameter or init)"
)
ing <- make_kl_class_inits(
  sdg,
  pg,
  kgs,
  2L,
  strategy = "given",
  mu_modes = list(-2, 1.5)
)
expect(
  abs(ing[[1]]$kmu_free[1] + 2) < 0.3 && abs(ing[[2]]$kmu_free[2] - 1.5) < 0.3,
  "init given: every chain starts each class at its init_mu (+- jitter)"
)
kgr <- kgs
kgr$gq_row_loglik <- TRUE
expect(
  inherits(
    try(
      suppressWarnings(build_kl_latent_standata(
        modifyList(m, list(kcfg = kgr)),
        pg,
        NULL,
        threads = 2L
      )),
      silent = TRUE
    ),
    "try-error"
  ),
  "gq_row_loglik on the country path refused"
)
kgd <- kd3
kgd$ksd_fixed <- 1
expect(
  inherits(try(kl_engine(kgd), silent = TRUE), "try-error"),
  "ksd_fixed with form discrete refused"
)
mk_rank("gc", c(-10, -11, -30, -30, -30, -10.5, -30)) # gc: modes -3 and 2
sm <- kl_sweep_modes("kl_fake_sweep", c("ga", "gc", "gz"), outdir = sw_out)
expect(
  identical(sm$ga, 0) && identical(sm$gc, c(-3, 2)) && is.null(sm$gz),
  "sweep modes: top lag, then separated modes with weight >= 0.05; missing group NULL"
)
in4 <- make_kl_class_inits(
  sdg,
  pg,
  kgs,
  4L,
  strategy = "sweep",
  mu_modes = list(c(-3, 2), 5)
)
expect(
  abs(in4[[1]]$kmu_free[1] + 3) < 0.3 &&
    abs(in4[[2]]$kmu_free[1] - 2) < 0.3 &&
    abs(in4[[3]]$kmu_free[1] + 3) < 0.3 &&
    is.null(in4[[1]]$ksd) &&
    !is.null(in4[[1]]$kre_zc),
  "sweep inits: chains alternate across a class's modes; no ksd init when fixed; slopes at 0"
)

# [5] artefact names ------------------------------------------------------------------
cat(
  "\n[4h] identify repel: repulsive prior, label-free registry pull, post-hoc naming\n"
)
krp <- kcfg
krp$identify <- "repel"
krp$classes <- 3L
expect(kl_identify(krp) == "repel", "identify repel accepted for form gaussian")
for (bad in list(
  list(form = "dlnm"),
  list(form = "discrete"),
  list(seeds = list(iso3 = list(g_east = "C01"))),
  list(classes = 1L),
  list(repel = list(lag_scale = 0)),
  list(repel = list(erf_grid = c(0, 1)))
)) {
  expect(
    inherits(
      try(kl_identify(modifyList(krp, bad)), silent = TRUE),
      "try-error"
    ),
    sprintf(
      "identify repel refused with %s",
      paste(names(unlist(bad)), unlist(bad), collapse = ", ")
    )
  )
}
rc <- kl_repel_cfg(krp)
expect(
  identical(
    rc,
    list(lag_scale = 3, erf_scale = 0.5, erf_grid = c(-1.5, -0.75, 0.75, 1.5))
  ),
  "repel defaults: lag_scale 3, erf_scale 0.5, grid -1.5 / -0.75 / 0.75 / 1.5"
)
expect(
  identical(
    kl_repel_cfg(
      kcfg_repel_def <- modifyList(
        krp,
        list(repel = get_default_options()$occurrence$kernel$repel)
      )
    ),
    rc
  ),
  "kl_repel_cfg() fallbacks == get_default_options() occurrence.kernel.repel"
)
krp$prior_grouping <- list(source = "kl_toy", membership = TRUE)
lv4 <- levels(droplevels(m$data$country))
pr <- kl_prior_partition(krp, lv4, groupings_dir = gdir)
expect(
  pr$use &&
    pr$K == 3L &&
    identical(pr$class_labels, c("klc1", "klc2", "klc3")) &&
    identical(colnames(pr$A), c("g_east", "g_south")) &&
    identical(pr$groups, c("g_east", "g_south")),
  "repel partition: K = kernel.classes (3 classes, 2 groups), A = countries x groups, labels klc<k>"
)
pr2 <- kl_prior_partition(
  modifyList(krp, list(classes = 2L)),
  lv4,
  groupings_dir = gdir
)
expect(pr2$K == 2L && ncol(pr2$A) == 2L, "repel partition: K = 2 with 2 groups")
krc <- krp
krc$spatial$adjacency <- "kl_toy_adj"
krc$prior_grouping$contiguity <- list(enable = TRUE, weight = 0.5)
pcr <- suppressMessages(kl_prior_partition(krc, lv4, groupings_dir = gdir))
expect(
  isTRUE(all.equal(unname(pcr$A["C01", ]), c(0.75, 0.25))),
  "repel + membership: contiguity mix accepted"
)
krn <- krc
krn$prior_grouping$membership <- NULL
expect(
  inherits(
    try(
      suppressMessages(kl_prior_partition(krn, lv4, groupings_dir = gdir)),
      silent = TRUE
    ),
    "try-error"
  ),
  "repel without membership: contiguity refused (A unused)"
)
mr <- modifyList(m, list(kcfg = krp))
sdr <- suppressWarnings(build_kl_latent_standata(mr, pr, NULL))
expect(
  sdr$kl_use_pg == 1L &&
    all(sdr$kl_A == 0) &&
    identical(dim(sdr$kl_A), c(4L, 3L)) &&
    sdr$kl_rep_G == 2L &&
    identical(dim(sdr$kl_rep_A), c(4L, 2L)) &&
    sdr$kl_mu_ordered == 0L &&
    sdr$kl_erf_ordered == 0L &&
    sdr$kl_field_cols == 3L &&
    sdr$kl_mu_normal == 0L,
  "repel data: kappa on, kappa * A zero, C x G pull matrix, free centres, symmetric field, uniform kmu"
)
expect(
  isTRUE(all.equal(
    sdr$kl_rep_V,
    .kl_basis_centred(m$meta, rc$erf_grid, 0),
    check.attributes = FALSE
  )) &&
    sdr$kl_rep_NX == 4L &&
    sdr$kl_rep_lag_scale == 3 &&
    sdr$kl_rep_erf_scale == 0.5,
  "repel data: centred exposure basis at the grid, scales"
)
krp0 <- krp
krp0$prior_grouping$membership <- NULL
sdr0 <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = krp0)),
  pr,
  NULL
))
expect(
  sdr0$kl_use_pg == 0L &&
    sdr0$kl_rep_G == 0L &&
    identical(dim(sdr0$kl_rep_A), c(4L, 0L)),
  "repel, membership off: no pull (registry only names the classes)"
)
expect(
  is.null(
    sdx_def <- suppressWarnings(build_kl_latent_standata(m, pg, NULL))$kl_rep_NX
  ),
  "default data carries no repel entries (the default program declares none)"
)
code_r <- patch_kl_latent_class_stancode(sc, m$meta, repel = TRUE)
Lr <- strsplit(code_r, "\n", fixed = TRUE)[[1]]
db_r <- Lr[(which(Lr == "data {") + 1):(which(Lr == "transformed data {") - 2)]
ids_r <- sub(
  ".*?([A-Za-z_][A-Za-z0-9_]*);.*$",
  "\\1",
  db_r[!grepl("^\\s*//", db_r) & grepl(";", db_r)],
  perl = TRUE
)
expect(
  all(ids_r %in% names(sdr)) && all(c("kl_rep_V", "kl_rep_A") %in% ids_r),
  sprintf(
    "repel program: every declared data identifier supplied (missing: %s)",
    paste(setdiff(ids_r, names(sdr)), collapse = ",")
  )
)
strip <- as.character(code_r)
for (fr in c(.KL_STAN_REPEL_DATA, .KL_STAN_REPEL_PARAMS, .KL_STAN_REPEL_TP)) {
  strip <- sub(paste0("\n", fr), "", strip, fixed = TRUE)
}
expect(
  identical(
    strsplit(strip, "\n", fixed = TRUE)[[1]][-1],
    strsplit(code, "\n", fixed = TRUE)[[1]][-1]
  ) &&
    grepl("identify repel", Lr[1]),
  "repel program = the default program plus the repel fragments (header notes identify repel)"
)
if (!inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")) {
  f <- tempfile(fileext = ".stan")
  writeLines(code_r, f)
  expect(
    isTRUE(cmdstanr::cmdstan_model(f, compile = FALSE)$check_syntax(
      quiet = TRUE
    )),
    "repel program passes stanc"
  )
}
inr <- make_kl_class_inits(sdr, pr, krp, 2L, strategy = "spread")
expect(
  identical(dim(inr[[1]]$kgrp), c(2L, 3L)) &&
    all(inr[[1]]$kgrp == 0) &&
    length(inr[[1]]$kmu_free) == 3L &&
    length(unique(round(inr[[1]]$kmu_free, 6))) == 3L,
  "repel inits: zero group logits, three distinct free centres"
)
expect(kl_init_strategy(krp, sdr) == "spread", "repel: init auto = spread")
expect(
  grepl("-K3-pgkltoy-mem.*-idrepel$", .kl_token(krp)) &&
    grepl(
      "-idrepel-rl2-rs0p25-rgm1_1-a0p05$",
      .kl_token(modifyList(
        krp,
        list(
          alpha = 0.05,
          repel = list(lag_scale = 2, erf_scale = 0.25, erf_grid = c(-1, 1))
        )
      ))
    ) &&
    !grepl(
      "-a",
      sub(".*-idrepel", "", .kl_token(modifyList(krp, list(alpha = 1))))
    ) &&
    grepl(
      "-idseed$",
      .kl_token(modifyList(krp, list(identify = "seed", alpha = 0.05)))
    ),
  "artefact token: -idrepel, then repel scales / grid / alpha only when not default (alpha only under repel)"
)

cat("\n[4i] identify seed + kernel.repel.enable: seeded labels, repulsion on\n")
ksr <- kgs
ksr$repel <- list(enable = TRUE)
expect(
  kl_identify(ksr) == "seed" &&
    kl_repel_on(ksr) &&
    !kl_repel_on(kgs) &&
    kl_repel_on(krp) &&
    !kl_repel_on(modifyList(krp, list(identify = "anchor"))),
  "kl_repel_on: identify repel always, seed only with repel.enable"
)
for (bad in list(
  list(identify = "anchor"),
  list(identify = "kmu"),
  list(form = "dlnm"),
  list(repel = list(lag_scale = -1))
)) {
  expect(
    inherits(
      try(kl_identify(modifyList(ksr, bad)), silent = TRUE),
      "try-error"
    ),
    sprintf(
      "repel.enable refused with %s",
      paste(names(unlist(bad)), unlist(bad), collapse = ", ")
    )
  )
}
sds <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = ksr)),
  pg,
  NULL,
  threads = 4L
))
expect(
  sds$kl_rep_NX == 4L &&
    sds$kl_rep_G == 0L &&
    identical(dim(sds$kl_rep_A), c(nlevels(droplevels(m$data$country)), 0L)) &&
    identical(sds$kl_A, sdg$kl_A) &&
    sds$kl_use_pg == sdg$kl_use_pg &&
    sds$kl_field_cols == sdg$kl_field_cols &&
    identical(as.integer(sds$kl_fix), as.integer(sdg$kl_fix)),
  "seed + repel data: repulsion grid / scales, no label-free pull (G = 0); seeds, labelled pull, field as seed"
)
expect(
  all(ids_r %in% names(sds)),
  sprintf(
    "seed + repel: every data identifier of the repel program supplied (missing: %s)",
    paste(setdiff(ids_r, names(sds)), collapse = ",")
  )
)
ins <- make_kl_class_inits(sds, pg, ksr, 1L, strategy = "spread")[[1]]
expect(
  is.null(ins$kgrp),
  "seed + repel inits: no group logits (no label-free pull)"
)
expect(
  grepl(
    "-idseed-rp-a0p5$",
    .kl_token(modifyList(
      krp,
      list(identify = "seed", alpha = 0.5, repel = list(enable = TRUE))
    ))
  ) &&
    grepl(
      "-idseed-rp-rl2$",
      .kl_token(modifyList(
        krp,
        list(identify = "seed", repel = list(enable = TRUE, lag_scale = 2))
      ))
    ) &&
    grepl(
      "-idseed$",
      .kl_token(modifyList(
        krp,
        list(identify = "seed", alpha = 0.5, repel = list(enable = FALSE))
      ))
    ),
  "artefact token: -idseed-rp (+ repel scales / alpha when not default) only with repel.enable"
)
ksi <- ksr
ksi$spatial$enable <- TRUE
ksi$spatial$prior <- "icar"
sdsi <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = ksi)),
  pg,
  adj,
  threads = 4L
))
expect(
  sdsi$kl_use_sp == 1L &&
    sdsi$kl_sp_dagar == 0L &&
    sdsi$kl_field_cols == pg$K - 1L &&
    sdsi$kl_n_edges == 3L &&
    sdsi$kl_rep_NX == 4L &&
    sdsi$kl_rep_G == 0L &&
    identical(as.integer(sdsi$kl_fix), as.integer(sds$kl_fix)),
  "seed + repel + ICAR data: field with a reference class (K - 1 columns), repulsion and seeds as seed + repel"
)
expect(
  all(ids_r %in% names(sdsi)),
  sprintf(
    "seed + repel + ICAR: every data identifier of the repel program supplied (missing: %s)",
    paste(setdiff(ids_r, names(sdsi)), collapse = ",")
  )
)
expect(
  grepl("-idseed-rp-a0p5$", .kl_token(modifyList(ksi, list(alpha = 0.5)))),
  "seed + repel + ICAR: kernel token as seed + repel (the -sp token comes from the spatial block)"
)

cat("\n[4j] kernel.lag_re: country lag offsets (program variant)\n")
kll <- kgs
kll$lag_re <- list(enable = TRUE, sd = 2)
expect(
  kl_lag_re_on(kll) &&
    !kl_lag_re_on(kgs) &&
    identical(kl_lag_re_cfg(kll)$sd, 2) &&
    !inherits(try(kl_engine(kll), silent = TRUE), "try-error"),
  "kl_lag_re_on / kl_lag_re_cfg: on with lag_re.enable, sd 2, accepted by kl_engine"
)
for (bad in list(
  list(form = "dlnm"),
  list(form = "discrete"),
  list(lag_re = list(enable = TRUE, sd = 0)),
  list(lag_re = list(enable = TRUE, sd = -1))
)) {
  expect(
    inherits(
      try(kl_lag_re_cfg(modifyList(kll, bad)), silent = TRUE),
      "try-error"
    ),
    sprintf(
      "lag_re refused with %s",
      paste(names(unlist(bad)), unlist(bad), collapse = ", ")
    )
  )
}
k1 <- resolve_kernel_cfg(kl_synthetic_opt(
  kernel = list(lag_re = list(enable = TRUE))
))
expect(
  inherits(try(kl_engine(k1), silent = TRUE), "try-error"),
  "lag_re without latent classes refused"
)
k2d <- resolve_kernel_cfg(kl_synthetic_opt(kernel = list(classes = 2L)))
expect(
  !kl_lag_re_on(k2d) &&
    isTRUE(all.equal(as.numeric(k2d$lag_re[["sd"]] %||% 2), 2)),
  "resolve_kernel_cfg: the default lag_re (enable no, sd 2) is carried and off"
)
code_l <- patch_kl_latent_class_stancode(sc, m$meta, lag_re = TRUE)
Ll <- strsplit(code_l, "\n", fixed = TRUE)[[1]]
expect(
  grepl("; lag_re\\)$", Ll[1]) &&
    !grepl("lag_re", strsplit(code, "\n", fixed = TRUE)[[1]][1]),
  "lag_re program header notes lag_re; the default header does not"
)
expect(
  identical(
    as.character(patch_kl_latent_class_stancode(sc, m$meta, lag_re = FALSE)),
    as.character(code)
  ),
  "lag_re = FALSE: the default program, byte-identical"
)
lik_l <- .kl_lag_re_likelihood(.KL_STAN_LIKELIHOOD, "model")
expect(
  !grepl("kl_KBu \\* ", lik_l) &&
    lengths(regmatches(
      lik_l,
      gregexpr("kl_lag_re_bw(", lik_l, fixed = TRUE)
    )) ==
      5L &&
    grepl("kl_ptc, kl_pa", lik_l, fixed = TRUE) &&
    grepl("kl_small, kl_tixc", lik_l, fixed = TRUE) &&
    grepl("kl_tixc, kl_Yv, kl_U, kl_BW", lik_l, fixed = TRUE),
  "lag_re likelihood: every per-class kernel term stacked by country (pair, country small / big, plain paths)"
)
expect(
  inherits(
    try(
      .kl_lag_re_likelihood(
        sub(
          "kl_U = kl_KBu * ktheta;",
          "kl_U = ktheta;",
          .KL_STAN_LIKELIHOOD,
          fixed = TRUE
        ),
        "model"
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  "lag_re substitution: a changed shared likelihood text fails loudly (pattern counts)"
)
for (fr in c(
  .KL_STAN_LAG_RE_FUNCTIONS,
  .KL_STAN_LAG_RE_DATA,
  .KL_STAN_LAG_RE_TDATA_DECL,
  .KL_STAN_LAG_RE_TDATA_STMT,
  .KL_STAN_LAG_RE_PARAMS,
  .KL_STAN_LAG_RE_TPARAMS_DECL,
  .KL_STAN_LAG_RE_TP,
  .KL_STAN_LAG_RE_MODEL
)) {
  expect(
    grepl(sub("^\n", "", fr), code_l, fixed = TRUE),
    sprintf(
      "lag_re program carries fragment %s",
      sub(
        "^\\s*//\\s*",
        "",
        strsplit(sub("^\n", "", fr), "\n", fixed = TRUE)[[1]][1]
      )
    )
  )
}
db_l <- Ll[(which(Ll == "data {") + 1):(which(Ll == "transformed data {") - 2)]
ids_l <- sub(
  ".*?([A-Za-z_][A-Za-z0-9_]*);.*$",
  "\\1",
  db_l[!grepl("^\\s*//", db_l) & grepl(";", db_l)],
  perl = TRUE
)
kpl <- kll
kpl$pair_collapse <- FALSE # country path (no pairs)
kpn <- kll
kpn$re_slope <- list(enable = FALSE) # pair / country path without slopes
sdl <- list(
  pairs_re = suppressWarnings(build_kl_latent_standata(
    modifyList(m, list(kcfg = kll)),
    pg,
    NULL,
    threads = 4L
  )),
  country_re = suppressWarnings(build_kl_latent_standata(
    modifyList(m, list(kcfg = kpl)),
    pg,
    NULL,
    threads = 4L
  )),
  pairs = suppressWarnings(build_kl_latent_standata(
    modifyList(m, list(kcfg = kpn)),
    pg,
    NULL,
    threads = 4L
  )),
  plain = suppressWarnings(build_kl_latent_standata(
    modifyList(m, list(kcfg = kpn)),
    pg,
    NULL,
    threads = 1L
  ))
)
expect(
  sdl$pairs_re$kl_pairs == 1L &&
    sdl$pairs_re$kl_re == 1L &&
    sdl$country_re$kl_pairs == 0L &&
    sdl$country_re$kl_rs == 1L &&
    sdl$pairs$kl_re == 0L &&
    sdl$pairs$kl_pairs == 1L &&
    sdl$plain$kl_rs == 0L,
  "lag_re data on every path: pairs + slopes, country + slopes, pairs, plain"
)
for (nm in names(sdl)) {
  expect(
    all(ids_l %in% names(sdl[[nm]])) && identical(sdl[[nm]]$kl_lag_re_sd, 2),
    sprintf(
      "lag_re program (%s): every declared data identifier supplied (missing: %s), kl_lag_re_sd 2",
      nm,
      paste(setdiff(ids_l, names(sdl[[nm]])), collapse = ",")
    )
  )
}
expect(
  is.null(sdg$kl_lag_re_sd),
  "default data carries no lag_re entries (the default program declares none)"
)
if (!inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")) {
  f <- tempfile(fileext = ".stan")
  writeLines(code_l, f)
  expect(
    isTRUE(cmdstanr::cmdstan_model(f, compile = FALSE)$check_syntax(
      quiet = TRUE
    )),
    "lag_re program passes stanc"
  )
  code_lr <- patch_kl_latent_class_stancode(
    sc,
    m$meta,
    repel = TRUE,
    lag_re = TRUE
  )
  f2 <- tempfile(fileext = ".stan")
  writeLines(code_lr, f2)
  expect(
    isTRUE(cmdstanr::cmdstan_model(f2, compile = FALSE)$check_syntax(
      quiet = TRUE
    )) &&
      grepl(
        "identify repel; lag_re",
        strsplit(code_lr, "\n", fixed = TRUE)[[1]][1],
        fixed = TRUE
      ),
    "repel + lag_re program passes stanc"
  )
}
inl <- make_kl_class_inits(sdl$pairs_re, pg, kll, 2L, strategy = "given")[[1]]
expect(
  identical(inl$klag_z, as.array(rep(0, sdl$pairs_re$kl_C))) &&
    identical(inl$klag_sd, as.array(0.5)) &&
    is.null(
      make_kl_class_inits(sdg, pg, kgs, 1L, strategy = "spread")[[1]]$klag_z
    ),
  "lag_re inits: zero offsets, sd 0.5; none without lag_re"
)
expect(
  grepl(
    "-idseed-lre2$",
    .kl_token(modifyList(kgs, list(lag_re = list(enable = TRUE))))
  ) &&
    grepl(
      "-idseed-lre1p5$",
      .kl_token(modifyList(kgs, list(lag_re = list(enable = TRUE, sd = 1.5))))
    ) &&
    !grepl(
      "-lre",
      .kl_token(modifyList(kgs, list(lag_re = list(enable = FALSE, sd = 1.5))))
    ),
  "artefact token: -lre<sd> only with lag_re.enable"
)

cat(
  "\n[4k] kernel.lag_re.scale: month (pinned) vs logit (centres inside the lag window)\n"
)
GOLDEN_LR <- "analysis/tests/golden/kl_latent_program_lag_re.stan"
expect(
  identical(
    paste(readLines(GOLDEN_LR), collapse = "\n"),
    sub("\n$", "", as.character(code_l))
  ) &&
    identical(
      as.character(patch_kl_latent_class_stancode(
        sc,
        m$meta,
        lag_re = TRUE,
        lag_re_scale = "month"
      )),
      as.character(code_l)
    ),
  "scale month: the lag_re program byte-identical to the pinned golden (kl_latent_program_lag_re.stan)"
)
expect(
  identical(kl_lag_re_cfg(kll)$scale, "month") &&
    identical(
      kl_lag_re_cfg(modifyList(
        kll,
        list(lag_re = list(scale = "logit"))
      ))$scale,
      "logit"
    ) &&
    identical(get_default_options()$occurrence$kernel$lag_re$scale, "month"),
  "kl_lag_re_cfg: scale month by default (and in get_default_options), logit accepted"
)
for (bad in list("months", c("month", "logit"), 1)) {
  expect(
    inherits(
      try(
        kl_lag_re_cfg(modifyList(kll, list(lag_re = list(scale = bad)))),
        silent = TRUE
      ),
      "try-error"
    ),
    sprintf("lag_re.scale '%s' refused", paste(bad, collapse = ","))
  )
}
kllg <- modifyList(kll, list(lag_re = list(scale = "logit", sd = 0.2)))
code_g <- patch_kl_latent_class_stancode(
  sc,
  m$meta,
  lag_re = TRUE,
  lag_re_scale = "logit"
)
Lg <- strsplit(code_g, "\n", fixed = TRUE)[[1]]
expect(
  grepl("; lag_re logit\\)$", Lg[1]),
  "logit program header notes 'lag_re logit'"
)
for (fr in c(
  .KL_STAN_LAG_RE_FUNCTIONS_LOGIT,
  .KL_STAN_LAG_RE_DATA_LOGIT,
  .KL_STAN_LAG_RE_TDATA_DECL,
  .KL_STAN_LAG_RE_TDATA_STMT,
  .KL_STAN_LAG_RE_PARAMS_LOGIT,
  .KL_STAN_LAG_RE_TPARAMS_DECL_LOGIT,
  .KL_STAN_LAG_RE_TP_LOGIT,
  .KL_STAN_LAG_RE_MODEL
)) {
  expect(
    grepl(sub("^\n", "", fr), code_g, fixed = TRUE),
    sprintf(
      "logit program carries fragment %s",
      sub(
        "^\\s*//\\s*",
        "",
        strsplit(sub("^\n", "", fr), "\n", fixed = TRUE)[[1]][1]
      )
    )
  )
}
nbw <- function(x, f) lengths(regmatches(x, gregexpr(f, x, fixed = TRUE)))
expect(
  nbw(code_g, "kl_lag_re_bw_cent(kl_KBl, kl_lags, kccent, kksd)") == 8L &&
    nbw(code_g, "kl_lag_re_bw(kl_KBl") == 0L &&
    !grepl("kmu[k] - d[c]", code_g, fixed = TRUE),
  "logit program: every kernel term reads the centres kccent (model + gq), no month kernel left"
)
lik_g <- .kl_lag_re_likelihood(.KL_STAN_LIKELIHOOD, "model", "logit")
expect(
  nbw(lik_g, "kl_lag_re_bw_cent(") == nbw(lik_l, "kl_lag_re_bw(") &&
    !grepl("kl_KBu \\* ", lik_g),
  "logit likelihood: the same substitutions as month, on kl_lag_re_bw_cent"
)
# R mirror of the Stan centre map (the formula the summaries rely on)
cent <- function(kmu, d, lo = -11, hi = 11) {
  p <- pmin(pmax((kmu - lo) / (hi - lo), 1e-9), 1 - 1e-9)
  lo + (hi - lo) * stats::plogis(stats::qlogis(p) + d)
}
expect(
  isTRUE(all.equal(cent(3, 0), 3)) &&
    all(cent(c(-11, -5, 0, 9, 11), c(-50, 3, 8, 20, 50)) >= -11) &&
    all(cent(c(-11, -5, 0, 9, 11), c(-50, 3, 8, 20, 50)) <= 11) &&
    is.finite(cent(11, -1)) &&
    cent(11, -1) < 11 &&
    isTRUE(all.equal(
      (cent(0, 1e-4) - cent(0, -1e-4)) / 2e-4,
      22 / 4,
      tolerance = 1e-6
    )),
  "centre map: zero offset = kmu, always inside [-11, 11], finite at the window edge, slope (hi - lo) / 4 at the centre"
)
sdg_l <- suppressWarnings(build_kl_latent_standata(
  modifyList(m, list(kcfg = kllg)),
  pg,
  NULL,
  threads = 4L
))
dbg <- Lg[(which(Lg == "data {") + 1):(which(Lg == "transformed data {") - 2)]
ids_g <- sub(
  ".*?([A-Za-z_][A-Za-z0-9_]*);.*$",
  "\\1",
  dbg[!grepl("^\\s*//", dbg) & grepl(";", dbg)],
  perl = TRUE
)
expect(
  all(ids_g %in% names(sdg_l)) && identical(sdg_l$kl_lag_re_sd, 0.2),
  "logit program: every declared data identifier supplied, kl_lag_re_sd 0.2 (logit units)"
)
if (!inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")) {
  f3 <- tempfile(fileext = ".stan")
  writeLines(code_g, f3)
  expect(
    isTRUE(cmdstanr::cmdstan_model(f3, compile = FALSE)$check_syntax(
      quiet = TRUE
    )),
    "logit lag_re program passes stanc"
  )
  code_gr <- patch_kl_latent_class_stancode(
    sc,
    m$meta,
    repel = TRUE,
    lag_re = TRUE,
    lag_re_scale = "logit"
  )
  f4 <- tempfile(fileext = ".stan")
  writeLines(code_gr, f4)
  expect(
    isTRUE(cmdstanr::cmdstan_model(f4, compile = FALSE)$check_syntax(
      quiet = TRUE
    )),
    "repel + logit lag_re program passes stanc"
  )
}
expect(
  grepl("-idseed-lrl0p2$", .kl_token(kllg)) &&
    grepl("-idseed-lre2$", .kl_token(kll)) &&
    !grepl(
      "-lr",
      .kl_token(modifyList(kllg, list(lag_re = list(enable = FALSE))))
    ),
  "artefact token: -lrl<sd> under scale logit, -lre<sd> under month, none when off"
)
# summaries / figure: the logit draws carry kccent, read in place of kmu + kclag
fk <- list(
  K = 2L,
  country_levels = c("AAA", "BBB"),
  class_labels = c("enso_a", "enso_b"),
  draws = posterior::as_draws_df(data.frame(
    `kmu[1]` = c(-2, -2),
    `kmu[2]` = c(5, 5),
    `kclag[1]` = c(30, 30),
    `kclag[2]` = c(-1, -1),
    `kccent[1,1]` = cent(-2, 30),
    `kccent[2,1]` = cent(5, 30),
    `kccent[1,2]` = cent(-2, -1),
    `kccent[2,2]` = cent(5, -1),
    `kl_log_resp[1,1]` = log(c(0.999, 0.999)),
    `kl_log_resp[2,1]` = log(c(0.001, 0.001)),
    `kl_log_resp[1,2]` = log(c(0.001, 0.001)),
    `kl_log_resp[2,2]` = log(c(0.999, 0.999)),
    `klag_sd[1]` = c(1, 1),
    check.names = FALSE
  ))
)
cs <- kl_country_lag_summary(fk)
expect(
  all(cs$scale == "logit") &&
    all(cs$lag_mean <= 11) &&
    cs$lag_mean[1] > 10.9 &&
    isTRUE(all.equal(
      cs$lag_mean[2],
      0.001 * cent(-2, -1) + 0.999 * cent(5, -1)
    )),
  "kl_country_lag_summary (logit): country lag from kccent, inside the window"
)
fm <- fk
fm$draws <- posterior::subset_draws(
  fk$draws,
  variable = "kccent",
  exclude = TRUE
)
expect(
  all(kl_country_lag_summary(fm)$scale == "month") &&
    kl_country_lag_summary(fm)$lag_mean[1] > 27,
  "kl_country_lag_summary (month): kmu + kclag, unbounded"
)

# Naming: exchangeable classes relabelled per draw, then named after the registry.
# 3 classes over 5 countries; true classes: A = {X1, X2}, B = {X3}, C = {X4, X5}.
# Half the draws hold the classes in another order (label switching).
mk_draws <- function(perms, lp, kmu_true) {
  truth <- rbind(
    c(.9, .9, .05, .05, .05),
    c(.05, .05, .9, .05, .05),
    c(.05, .05, .05, .9, .9)
  )
  rows <- lapply(seq_len(nrow(perms)), function(i) {
    p <- perms[i, ] # sampled class s holds true class p[s]
    v <- c(lp__ = lp[i])
    for (s in 1:3) {
      for (c in 1:5) {
        v[sprintf("kl_log_resp[%d,%d]", s, c)] <- log(truth[p[s], c])
      }
    }
    for (s in 1:3) {
      v[sprintf("kmu[%d]", s)] <- kmu_true[p[s]]
      v[sprintf("kcon[%d]", s)] <- p[s] / 10
      v[sprintf("kbeta[1,%d]", s)] <- p[s]
      v[sprintf("kpi[%d]", s)] <- c(.4, .2, .4)[p[s]]
    }
    v
  })
  posterior::as_draws_df(as.data.frame(do.call(rbind, rows)))
}
pm <- rbind(
  matrix(1:3, 6, 3, byrow = TRUE),
  matrix(c(3, 1, 2), 6, 3, byrow = TRUE)
)
dr <- mk_draws(pm, lp = c(-10, rep(-12, 11)), kmu_true = c(4, -5, 9))
Areg <- matrix(
  c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1),
  5,
  2,
  dimnames = list(sprintf("X%d", 1:5), c("g_a", "g_c"))
)
rr <- kl_repel_relabel(dr, 3L, Areg)
rd <- as.data.frame(rr$draws)
expect(
  identical(rr$labels, c("g_a", "g_c", "klc3")),
  "naming, K > G: classes named after their registry group in registry order, the extra class klc3"
)
expect(
  all(rd[["kmu[1]"]] == 4) &&
    all(rd[["kmu[2]"]] == 9) &&
    all(rd[["kmu[3]"]] == -5) &&
    all(rd[["kbeta[1,2]"]] == 3) &&
    all(abs(rd[["kl_log_resp[3,3]"]] - log(.9)) < 1e-12),
  "relabelling: every draw aligned (kmu, kbeta, responsibilities follow their class)"
)
expect(
  identical(sort(rr$table$n), c(6L, 6L)) &&
    rr$naming$passes <= 3L &&
    identical(rr$naming$assignment, c("g_a", "g_c", NA)),
  "relabel table: two sampled labellings; naming record stored"
)
rr_b <- kl_repel_relabel(dr, 3L, Areg)
expect(
  identical(rr_b$labels, rr$labels) && identical(as.data.frame(rr_b$draws), rd),
  "relabelling and naming deterministic"
)
A3g <- cbind(Areg, g_z = 0)
A3g[, "g_c"] <- c(0, 0, 0, 1, 1)
A3g[, "g_z"] <- c(0, 0, 1, 0, 0)
A4g <- cbind(A3g, g_w = 0)
rr4 <- kl_repel_relabel(dr, 3L, A4g)
expect(
  identical(rr4$labels, c("g_a", "g_c", "g_z")),
  "naming, K < G: three classes take the three overlapping groups, the empty group names nothing"
)
rr0 <- kl_repel_relabel(dr, 3L, NULL)
expect(
  identical(rr0$labels, c("klc1", "klc2", "klc3")) &&
    identical(
      vapply(
        1:3,
        function(k) unique(as.data.frame(rr0$draws)[[sprintf("kmu[%d]", k)]]),
        0
      ),
      c(-5, 4, 9)
    ),
  "no registry: klc<k> by posterior median lag"
)
expect(
  identical(.kl_assign_max(matrix(c(1, 1, 0, 2), 2)), c(1L, 2L)) &&
    identical(.kl_assign_max(matrix(c(5, 4), 2, 1)), c(1L, NA)) &&
    nrow(.kl_perm_matrix(4L)) == 24L &&
    !anyDuplicated(apply(.kl_perm_matrix(4L), 1, paste, collapse = "")),
  "assignment helper: one-to-one maximum, left-over class NA; permutations complete"
)

# [5] artefact names distinct across lag range / engine / K / prior / spatial
cat("\n[5] artefact names\n")
o <- get_default_options()
o$baseline <- "1990_2019"
o$country <- "A-B"
o$occurrence$backend <- "brms_kl"
o$occurrence$dlnm <- list(
  enable = TRUE,
  defaults = list(argvar = list(fun = "ns", df = 3L)),
  variables = list(enso = list(lags = list(-11, 11)))
)
variant <- function(f) {
  x <- o
  f(x)
}
paths <- c(
  base = make_run_file(o, "fit"),
  lag6 = make_run_file(
    variant(function(x) {
      x$occurrence$dlnm$variables$enso$lags <- list(-6, 6)
      x
    }),
    "fit"
  ),
  brms = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$engine <- "brms"
      x
    }),
    "fit"
  ),
  k2 = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x
    }),
    "fit"
  ),
  pg = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$prior_grouping$source <- "enso_grouping"
      x
    }),
    "fit"
  ),
  sp = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$spatial$enable <- TRUE
      x
    }),
    "fit"
  ),
  iderf = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "erf"
      x
    }),
    "fit"
  ),
  dlnm = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$form <- "dlnm"
      x
    }),
    "fit"
  ),
  disc = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$form <- "discrete"
      x
    }),
    "fit"
  ),
  pin = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$form <- "discrete"
      x$occurrence$kernel$pin$enable <- TRUE
      x
    }),
    "fit"
  ),
  re = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$form <- "discrete"
      x$occurrence$kernel$re_slope$enable <- TRUE
      x
    }),
    "fit"
  ),
  gseed = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x
    }),
    "fit"
  ),
  gseeda = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x$occurrence$kernel$re_slope$country <- FALSE
      x
    }),
    "fit"
  ),
  gsdcf = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x$occurrence$kernel$re_slope$sd_country_fixed <- 0.2
      x
    }),
    "fit"
  ),
  gsdcp = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x$occurrence$kernel$re_slope$sd_country <- 0.15
      x
    }),
    "fit"
  ),
  gsda1 = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x$occurrence$kernel$re_slope$country <- FALSE
      x$occurrence$kernel$re_slope$sd_adm1 <- 1
      x
    }),
    "fit"
  ),
  gsdac = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$identify <- "seed"
      x$occurrence$kernel$ksd_fixed <- 1.5
      x$occurrence$kernel$re_slope$enable <- TRUE
      x$occurrence$kernel$re_slope$adm1 <- FALSE
      x$occurrence$kernel$re_slope$sd_adm1 <- 1
      x
    }),
    "fit"
  ),
  discw = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$form <- "discrete"
      x$occurrence$kernel$discrete$widths <- c(1, 2.5)
      x
    }),
    "fit"
  ),
  dagar = make_run_file(
    variant(function(x) {
      x$occurrence$kernel$classes <- 2L
      x$occurrence$kernel$spatial$enable <- TRUE
      x$occurrence$kernel$spatial$prior <- "dagar"
      x
    }),
    "fit"
  ),
  cb = make_run_file(
    variant(function(x) {
      x$occurrence$backend <- "brms_cb"
      x
    }),
    "fit"
  )
)
expect(
  !anyDuplicated(paths),
  "fit paths distinct across lag range / engine / K / prior / ICAR / DAGAR / backend"
)
expect(
  grepl("_cb-enso_lags-m11to11_ns3_kl\\.rds$", basename(paths[["base"]])),
  "kernel tokens in the name"
)
expect(
  grepl("-rec[-_]", basename(paths[["gsdac"]])) &&
    !grepl("-sda", basename(paths[["gsdac"]])) &&
    grepl("-re[-_]", basename(paths[["gseed"]])),
  "re_slope.adm1 = no: -rec token (country deviations only), no -sda token"
)
expect(
  grepl("-re-sdcp0p15", basename(paths[["gsdcp"]])) &&
    !grepl("-sdcp", basename(paths[["gseed"]])),
  "country sd prior scale: -sdcp<s> when not 0.35, no token at the default"
)
expect(
  make_run_file(o, "kl_resp") != make_run_file(o, "kl_partition"),
  "kl_* kinds distinct"
)

# [6] refusals -----------------------------------------------------------------------
cat(
  "\n[5b] defaults are not shadowed by longer keys (R `$` partial matching)\n"
)
od <- get_default_options()
od$occurrence$kernel$prior_mu <- NULL
od$occurrence$kernel$prior_grouping$kappa <- NULL # as after a merge
od$occurrence$dlnm <- list(
  enable = TRUE,
  defaults = list(argvar = list(fun = "ns", df = 3L)),
  variables = list(enso = list(lags = list(-11, 11)))
)
kd <- resolve_kernel_cfg(od)
expect(
  is.null(kd$prior_mu) && .kl_prior_strings(kd)$mu == "uniform(-11, 11)",
  "absent prior_mu -> uniform kmu prior (not prior_mu_sd)"
)
od$occurrence$backend <- "brms_kl"
od$occurrence$kernel$prior_grouping$source <- "enso_grouping"
expect(
  !grepl("-kap", basename(make_run_file(od, "fit"))),
  "absent kappa -> no fixed-kappa token (not kappa_prior_sd)"
)
sdd <- suppressWarnings(build_kl_latent_standata(
  m,
  kl_prior_partition(kcfg, lv),
  NULL
))
expect(
  sdd$kl_kappa_fixed == 0L && sdd$kl_mu_normal == 0L,
  "data: kappa estimated, kmu prior uniform by default"
)

cat("\n[5c] rows brms would drop are dropped up front\n")
pna <- panel
pna$t_start_seas[c(5, 40)] <- NA
mna <- kl_sort_by_country(suppressMessages(assemble_kl_model(
  opt,
  pna,
  list(re_uncorrelated = TRUE),
  .occ_baseline_block("cyclic"),
  kcfg
)))
sdna <- suppressWarnings(build_kl_latent_standata(
  mna,
  kl_prior_partition(kcfg, lv),
  NULL
))
expect(
  nrow(mna$data) == nrow(panel) - 2L &&
    sdna$N == nrow(mna$data) &&
    length(sdna$kl_ctry) == sdna$N,
  "NA rows dropped; brms N == country index length"
)

# [6] config refusals
cat("\n[6] config refusals\n")
bad <- function(f) inherits(try(f(), silent = TRUE), "try-error")
o2 <- opt
o2$occurrence$covariates$enso <- TRUE
expect(
  bad(function() resolve_kernel_cfg(o2)),
  "linear covariates.enso refused (double count)"
)
o3 <- opt
o3$occurrence$dlnm$variables <- list(mean_precip = list(lags = list(0, 3)))
expect(
  bad(function() resolve_kernel_cfg(o3)),
  "kernel variable other than enso refused"
)
k4 <- kcfg
k4$engine <- "brms"
expect(bad(function() kl_engine(k4)), "engine brms with latent classes refused")
k5 <- resolve_kernel_cfg(kl_synthetic_opt(
  kernel = list(spatial = list(enable = TRUE))
))
expect(
  bad(function() kl_engine(k5)),
  "spatial term without latent classes refused"
)
k6 <- kcfg
k6$spatial$prior <- "car"
expect(bad(function() kl_spatial_prior(k6)), "unknown spatial prior refused")

# [7] TTLO level none (occurrence.ttlo_no_prior_level): the 8-level factor
#     reaches the brms_kl model and its standata
cat("\n[7] TTLO level none (occurrence.ttlo_no_prior_level; kl sweep option)\n")
# The load-time recode of occ_data.R::load_surv_data() on a synthetic panel
# (some rows with no recorded prior onset: NA, or still at the TTLO origin),
# then the brms_kl model + latent-class standata, no sampling.
source("analysis/utils/occ_data.R")
p8 <- panel
set.seed(7)
d8 <- sample(c(10, 200, 500, 900, 1400, 1900, 2300, 3000), nrow(p8), TRUE)
i8 <- sample.int(nrow(p8), 30L)
d8[i8[1:15]] <- NA
d8[i8[16:30]] <- as.numeric(as.Date(p8$t_start[i8[16:30]])) -
  as.numeric(.TTLO_ORIGIN)
p8$time_to_last_outbreak <- d8
p8 <- make_time_to_last_outbreak_cat(p8, no_prior_level = TRUE)
expect(
  identical(levels(p8$time_to_last_outbreak), c(.TTLO_LEVELS, "none")) &&
    sum(p8$time_to_last_outbreak == "none") == 30L,
  "8-level factor, 30 rows at level none"
)
o8 <- opt
o8$occurrence$ttlo_no_prior_level <- TRUE
o8$occurrence$covariates$time_to_last_outbreak_cat <- TRUE
o8$occurrence$random_effects$re_region_country_ttlo_slope <- TRUE
k8 <- resolve_kernel_cfg(o8)
m8 <- kl_sort_by_country(suppressMessages(assemble_kl_model(
  o8,
  p8,
  list(re_uncorrelated = TRUE),
  .occ_baseline_block("cyclic"),
  k8
)))
expect(
  "time_to_last_outbreak" %in%
    m8$terms$fixed_terms &&
    any(grepl("^\\(time_to_last_outbreak - 1 \\|", m8$terms$re_terms)),
  "TTLO fixed term and country slope in the linear part"
)
sd8 <- suppressWarnings(build_kl_latent_standata(
  m8,
  kl_prior_partition(k8, lv),
  NULL
))
expect(
  sd8$N == nrow(m8$data) &&
    "time_to_last_outbreaknone" %in% colnames(sd8$X_lin) &&
    !any(grepl("^time_to_last_outbreak01", colnames(sd8$X_lin))),
  "standata: every row kept, a 'none' column in X_lin, [0,1) the reference"
)
expect(
  grepl(
    "_ttn$",
    tools::file_path_sans_ext(basename(make_occ_run_file(o8, "fit")))
  ) &&
    !grepl("_ttn", basename(make_occ_run_file(opt, "fit"))),
  "fit name token _ttn only with the level on"
)

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1)
}
cat("PASS  test_kl_stan_patch\n")
