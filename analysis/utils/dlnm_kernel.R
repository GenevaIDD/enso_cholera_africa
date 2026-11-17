# dlnm_kernel.R — Gaussian lag-kernel ENSO term and latent country classes for
#   the occurrence brms_kl backend (04_fit_occ_brms_kl.R, sandbox).
#
# Provides:
#   resolve_kernel_cfg(), kl_engine()   validated `occurrence.kernel` block
#   assemble_kl_model(), build_kl_formula(), build_kl_priors()   brms model
#   patch_kl_latent_class_stancode()    K-class country mixture program
#   build_kl_latent_standata(), fit_kl_stan()   data + fit (class kl_fit)
#   kl_relabel_draws(), kl_name_classes()       label classes after sampling
#   kl_pin_run(), kl_pin_plan(), kl_pin_fit_one(), kl_pin_finish()   pinning
#   summarise_kl_fit(), summarise_kl_latent(), kl_responsibilities()
#   kl_chain_agreement(), kl_map_partition(), write_kl_partition()
#   kl_loo_country()                    country-level LOO
#
# Details:
#   Callers source base.R, dlnm.R, formula.R, priors.R and config.R; also
#   deploy_core.R / paths.R only with kernel.prior_grouping.deploy_config,
#   and loo_lags.R only with kernel.pin.order_source sweep.
#   Packages via `::`: brms, dlnm, posterior, cmdstanr, loo.
#
# See: analysis/utils/README.md

# Stage 1 (single kernel): the single-lag ENSO onebasis of brms_cb is
# replaced by a normalised Gaussian kernel over the dlnm `lags:` range,
#
#   kl_i = sum_l w_l(kmu, ksd) * sum_j kb_j * B_j(ENSO at lag l, row i),
#   w_l  = exp(-(l - kmu)^2 / (2 ksd^2)) / sum_l' exp(-(l' - kmu)^2 / (2 ksd^2)),
#
# written as a brms non-linear formula (`lin` = the usual occurrence linear
# predictor; nlpars kmu, ksd, kb1..kbJ -- `mu`/`sigma` are reserved brms dpar
# names). B is ONE onebasis shared by every lag (build_cb_block() `lags:` mode,
# utils/dlnm.R), so kb_j means the same exposure function at every lag and
# ksd -> 0 recovers the single-lag model at lag kmu.
#
# Stage 2 (kernel.classes >= 2): K latent country classes, each with its own
# (kmu_k, ksd_k, kb_k), marginalised at country level. The Stan program is the
# stage-1 brms program (brms::stancode()) patched at fixed landmarks
# (patch_kl_latent_class_stancode()); membership logits can carry a prior
# partition (kernel.prior_grouping) and an ICAR spatial term
# (kernel.spatial). See the stage-2 section below.
#
# Forms (kernel.form): "gaussian" (above; kmu, ksd sampled), "dlnm" (class
# B-spline lag basis, group horseshoe) and "discrete": Gaussian kernels whose
# centres sit on the lag grid and whose widths come from a short list
# (kernel.discrete), the class centres ordered (ties allowed) and summed out
# exactly -- see "Form discrete" in the stage-2 section.
.KL_NLPAR <- list(mu = "kmu", sigma = "ksd", beta = "kb")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

#' Merged dlnm spec of one variable (shared defaults, per-variable wins), the
#' same merge build_cb_block() applies.
#' @noRd
.kl_var_spec <- function(dlnm_cfg, v) {
  spec <- dlnm_cfg$defaults %||% list()
  vs <- dlnm_cfg$variables[[v]] %||% list()
  for (k in names(vs)) {
    spec[[k]] <- vs[[k]]
  }
  spec
}

#' resolve_kernel_cfg
#'
#' Validated lag-kernel configuration from `opt$occurrence` (`kernel:` block +
#' the dlnm `lags:` of the kernel variable). Stops on configurations the
#' backend cannot fit: dlnm disabled, no or several `lags:` variables, a
#' kernel variable other than enso (lag columns come from join_enso_lags()), or
#' the linear `covariates.enso` term switched on (it would double-count ENSO).
#'
#' @return list(variable, lags, lag_min, lag_max, cen, sigma_lb, prior_mu,
#'   prior_mu_sd, prior_sigma_sd, prior_beta_sd, classes, alpha, init,
#'   drop_region, gq_row_loglik, keep_re_draws, partition_name,
#'   prior_grouping, spatial)
resolve_kernel_cfg <- function(opt) {
  occ <- opt$occurrence %||% list()
  dl <- occ$dlnm %||% list()
  k <- occ$kernel %||% list()
  if (!isTRUE(dl$enable)) {
    stop(
      "brms_kl: occurrence.dlnm.enable must be true -- the lag kernel reads ",
      "the dlnm `lags:` basis."
    )
  }
  lag_vars <- Filter(
    function(v) !is.null(.kl_var_spec(dl, v)$lags),
    names(dl$variables)
  )
  if (!identical(lag_vars, "enso")) {
    stop(
      "brms_kl: exactly one dlnm variable must set `lags:` and it must be ",
      "`enso` (found: ",
      if (length(lag_vars)) paste(lag_vars, collapse = ", ") else "none",
      ")."
    )
  }
  if (isTRUE(occ$covariates$enso)) {
    stop(
      "brms_kl: occurrence.covariates.enso must be false -- the linear ",
      "enso_lag<N> term would double-count the kernel term."
    )
  }
  spec <- .kl_var_spec(dl, "enso")
  if (isTRUE(spec$re_adm1_slope)) {
    warning(
      "brms_kl: re_adm1_slope is ignored for the kernel variable (out of scope)."
    )
  }
  lags <- .cb_lags_spec(spec$lags)
  pg <- k$prior_grouping %||% list()
  sp <- k$spatial %||% list()
  K <- as.integer(k$classes %||% 1L)
  list(
    variable = "enso",
    lags = lags,
    lag_min = min(lags),
    lag_max = max(lags),
    cen = as.numeric(k$cen %||% 0),
    sigma_lb = as.numeric(k$sigma_lb %||% 0.3),
    prior_mu = k[["prior_mu"]], # [[ ]]: `$` would partially match prior_mu_sd
    prior_mu_sd = as.numeric(k$prior_mu_sd %||% 3),
    prior_sigma_sd = as.numeric(k$prior_sigma_sd %||% 2),
    prior_beta_sd = as.numeric(k$prior_beta_sd %||% 0.5),
    classes = K,
    engine = k$engine %||% "stan",
    form = k$form %||% "gaussian",
    lag_df = as.integer(k$lag_df %||% 10L),
    lag_prior = k$lag_prior %||% "horseshoe",
    hs = k$hs %||% list(),
    discrete = k$discrete %||% list(),
    pin = k$pin %||% list(),
    re_slope = k$re_slope %||% list(),
    ksd_fixed = k[["ksd_fixed"]], # [[ ]]: exact key
    seeds = k$seeds %||% list(),
    sweep = k$sweep %||% list(),
    grainsize = as.integer(k$grainsize %||% 1L),
    grainsize_rows = as.integer(k$grainsize_rows %||% 1L),
    init_mu = k[["init_mu"]],
    split_big = !isFALSE(k$split_big),
    pair_collapse = !isFALSE(k[["pair_collapse"]]),
    identify = k$identify %||% "auto",
    repel = k$repel %||% list(),
    lag_re = k[["lag_re"]] %||% list(), # [[ ]]: exact key
    erf_at = as.numeric(k$erf_at %||% 1.5),
    alpha = as.numeric(k$alpha %||% 1),
    init = k$init %||% "auto",
    init_starts = as.integer(k$init_starts %||% 8L),
    drop_region = k$drop_region %||% (K >= 2L || !is.null(pg$source)),
    gq_row_loglik = isTRUE(k$gq_row_loglik),
    keep_re_draws = isTRUE(k$keep_re_draws),
    partition_name = k$partition_name,
    prior_grouping = pg,
    spatial = sp
  )
}

#' TRUE when the config asks for latent classes: two or more classes, or a
#' prior partition (whose groups set K).
kl_is_latent <- function(kcfg) {
  kcfg$classes >= 2L || !is.null(kcfg$prior_grouping$source)
}

#' Spatial membership prior: "icar" (intrinsic CAR, soft sum-to-zero per
#' component) or "dagar" (proper directed-acyclic-graph autoregressive field,
#' non-centred, beta prior on rho; isolated countries get an iid N(0, tau)).
kl_spatial_prior <- function(kcfg) {
  pr <- kcfg$spatial$prior %||% "icar"
  if (!pr %in% c("icar", "dagar")) {
    stop(
      "brms_kl: kernel.spatial.prior must be 'icar' or 'dagar' (got '",
      pr,
      "')."
    )
  }
  pr
}

#' kl_identify
#'
#' How the latent classes are identified (label switching broken):
#'   "kmu"    ordered lag centres (classes = lag rank);
#'   "anchor" prior-partition pull kappa * A + class kmu priors at the groups'
#'            deployed lags (class k = registry group k); needs prior_grouping;
#'   "erf"    ordered exposure-response contrast: log RR at kernel.erf_at vs
#'            cen (classes = rank of the El Nino response). With a prior
#'            partition, the registry sets K and names classes after sampling
#'            (majority of member countries), but adds no anchoring term.
#'   "relabel" (form dlnm) exchangeable classes, each draw relabelled after
#'            sampling by sorting on kcon (log RR at erf_at vs cen); the
#'            registry, if any, names classes (see kl_relabel_draws()).
#'   "lagorder" (form discrete) class centres ordered by lag, ties allowed;
#'            classes tied on a lag are relabelled after sampling by kcon
#'            (kl_relabel_draws(lag = "klag")); the registry, if any, names
#'            classes.
#'   "seed"   (forms gaussian / dlnm) class k holds the seed country (or
#'            countries) of prior group k (kernel.seeds; fixed membership,
#'            kl_fix) -- labels fixed by the seeds, centres unordered, one run;
#'            every other country is free. Needs prior_grouping. With
#'            prior_grouping.membership: yes the other countries also get the
#'            kappa * A pull towards their registry group (kmu priors unchanged).
#'   "repel"  (form gaussian) exchangeable classes, no seeds, kept apart by a
#'            pairwise repulsive prior on their distance in (lag centre,
#'            exposure-response) space (kernel.repel; kl_repel_cfg()), so two
#'            classes may share a lag when their ENSO responses differ. With a
#'            prior partition and prior_grouping.membership: yes, the registry
#'            pull is label-free (kappa * group-level class logits, any K);
#'            each draw is relabelled after sampling and the classes are named
#'            after the registry (kl_repel_relabel()).
#'   "auto"   = "lagorder" for form discrete, "relabel" for form dlnm;
#'            otherwise "anchor" with a prior partition, else "kmu" (the
#'            original rule).
kl_identify <- function(kcfg) {
  id <- kcfg$identify %||% "auto"
  if (
    !id %in%
      c("auto", "kmu", "anchor", "erf", "relabel", "lagorder", "seed", "repel")
  ) {
    stop(
      "brms_kl: kernel.identify must be auto | kmu | anchor | erf | relabel | lagorder | seed | repel (got '",
      id,
      "')."
    )
  }
  pg <- !is.null(kcfg$prior_grouping$source)
  if (identical(kcfg$form %||% "gaussian", "discrete")) {
    # Discrete centres: the ordering of the centres is the identification.
    if (identical(id, "auto")) {
      id <- "lagorder"
    }
    if (!identical(id, "lagorder")) {
      stop(
        "brms_kl: kernel.form discrete is identified by its ordered lag centres ",
        "(identify lagorder | auto; got '",
        id,
        "')."
      )
    }
    return(id)
  }
  if (identical(id, "lagorder")) {
    stop("brms_kl: kernel.identify lagorder is for kernel.form discrete.")
  }
  if (isTRUE(kcfg$repel[["enable"]]) && !id %in% c("seed", "repel")) {
    stop(
      "brms_kl: kernel.repel.enable is for identify seed (identify repel always repels; got '",
      id,
      "')."
    )
  }
  if (identical(id, "seed")) {
    if (!pg) {
      stop(
        "brms_kl: kernel.identify seed needs kernel.prior_grouping.source (the groups the seeds come from)."
      )
    }
    if (isTRUE(kcfg$repel[["enable"]])) {
      if (!identical(kcfg$form %||% "gaussian", "gaussian")) {
        stop(
          "brms_kl: kernel.repel.enable under identify seed is for kernel.form gaussian (got '",
          kcfg$form,
          "')."
        )
      }
      kl_repel_cfg(kcfg) # validates kernel.repel
    }
    return(id)
  }
  if (identical(id, "repel")) {
    if (!identical(kcfg$form %||% "gaussian", "gaussian")) {
      stop(
        "brms_kl: kernel.identify repel is for kernel.form gaussian (got '",
        kcfg$form,
        "')."
      )
    }
    if (!is.null(kcfg$seeds[["iso3"]])) {
      stop(
        "brms_kl: kernel.identify repel has no seed countries -- unset kernel.seeds.iso3."
      )
    }
    if ((kcfg$classes %||% 1L) < 2L) {
      stop(
        "brms_kl: kernel.identify repel needs kernel.classes >= 2 (K is not set by the registry)."
      )
    }
    kl_repel_cfg(kcfg) # validates kernel.repel
    return(id)
  }
  if (identical(kcfg$form %||% "gaussian", "dlnm")) {
    # Class lag-basis form: no kernel centre to order, and a constrained
    # contrast parametrisation would fight the non-centred horseshoe. Classes
    # are made exchangeable (symmetric spatial field, no class-specific priors)
    # and each draw is relabelled by sorting on kcon after sampling.
    if (identical(id, "auto")) {
      id <- "relabel"
    }
    if (!id %in% c("relabel", "anchor")) {
      stop(
        "brms_kl: kernel.form dlnm supports identify relabel | anchor (got '",
        id,
        "')."
      )
    }
  }
  if (
    identical(id, "relabel") && !identical(kcfg$form %||% "gaussian", "dlnm")
  ) {
    stop("brms_kl: kernel.identify relabel is for kernel.form dlnm.")
  }
  if (identical(id, "auto")) {
    id <- if (pg) "anchor" else "kmu"
  }
  if (identical(id, "anchor") && !pg) {
    stop(
      "brms_kl: kernel.identify 'anchor' needs kernel.prior_grouping.source."
    )
  }
  id
}

#' kl_repel_cfg
#'
#' The repulsive prior of identify repel (kernel.repel), validated. Every pair
#' of classes (j, k) adds log(1 - exp(-d2 / 2)) to the log prior, with
#'   d2 = ((kmu_j - kmu_k) / lag_scale)^2
#'        + mean_x (erf_j(x) - erf_k(x))^2 / erf_scale^2,
#' erf_k(x) = log RR of a sustained ENSO value x vs cen in class k, x over
#' erf_grid. Two identical classes have density 0; classes lag_scale months
#' apart, or whose responses differ by erf_scale (root mean square over the
#' grid), pay log(1 - exp(-1/2)) = -0.93; well-separated classes pay ~0.
#' @return list(lag_scale, erf_scale, erf_grid)
kl_repel_cfg <- function(kcfg) {
  out <- list(
    lag_scale = as.numeric(kcfg$repel$lag_scale %||% 3),
    erf_scale = as.numeric(kcfg$repel$erf_scale %||% 0.5),
    erf_grid = as.numeric(unlist(
      kcfg$repel$erf_grid %||% c(-1.5, -0.75, 0.75, 1.5)
    ))
  )
  for (nm in c("lag_scale", "erf_scale")) {
    if (length(out[[nm]]) != 1L || is.na(out[[nm]]) || out[[nm]] <= 0) {
      stop("brms_kl: kernel.repel.", nm, " must be one number > 0.")
    }
  }
  if (
    !length(out$erf_grid) ||
      anyNA(out$erf_grid) ||
      any(out$erf_grid == (kcfg$cen %||% 0))
  ) {
    stop("brms_kl: kernel.repel.erf_grid must be ENSO values away from cen.")
  }
  out
}

#' kl_repel_on
#'
#' Is the pairwise repulsion (kl_repel_cfg()) in the program? Always under
#' identify repel; under identify seed when kernel.repel.enable is yes (labels
#' stay fixed by the seeds, no label-free pull: kl_rep_G = 0).
kl_repel_on <- function(kcfg) {
  id <- kcfg$identify %||% "auto"
  identical(id, "repel") ||
    (identical(id, "seed") && isTRUE(kcfg$repel[["enable"]]))
}

#' kl_lag_re_on -- is the country lag offset (kernel.lag_re) in the program?
kl_lag_re_on <- function(kcfg) isTRUE(kcfg[["lag_re"]][["enable"]])

#' kl_lag_re_scale -- the scale of the country lag offset (kernel.lag_re.scale):
#' "month" (centre kmu_k + d_c, lags) or "logit" (centre inside the lag window).
kl_lag_re_scale <- function(kcfg) kcfg[["lag_re"]][["scale"]] %||% "month"

#' kl_lag_re_cfg
#'
#' The country lag offset of the Gaussian kernel (kernel.lag_re), validated:
#' one offset per country, d_c = sigma * z_c, z_c ~ normal(0, 1), sigma ~
#' half-normal(0, sd), shared by the classes (membership is summed out per
#' country). Country c in class k has lag centre
#'   scale month: kmu_k + d_c (d_c in lags; the centre can leave the lag
#'                window, where the kernel weights pile on the edge lag and the
#'                likelihood is flat in d_c);
#'   scale logit: lo + (hi - lo) * inv_logit(logit((kmu_k - lo) / (hi - lo)) + d_c)
#'                (d_c in logit units of the window [lo, hi] = the first / last
#'                lag; always inside it; 1 unit ~ (hi - lo) / 4 lags at the
#'                window centre, fewer towards the edges).
#' sd is in the units of d_c (lags for month, logit units for logit). Latent
#' classes and form gaussian only (form dlnm / discrete have no kernel centre
#' to shift).
#' @return list(sd, scale)
kl_lag_re_cfg <- function(kcfg) {
  lr <- kcfg[["lag_re"]] %||% list()
  if (!identical(kcfg$form %||% "gaussian", "gaussian")) {
    stop(
      "brms_kl: kernel.lag_re is for kernel.form gaussian (got '",
      kcfg$form,
      "')."
    )
  }
  if (!kl_is_latent(kcfg)) {
    stop(
      "brms_kl: kernel.lag_re needs latent classes (classes >= 2 or prior_grouping)."
    )
  }
  sd <- as.numeric(lr[["sd"]] %||% 2)
  if (length(sd) != 1L || is.na(sd) || sd <= 0) {
    stop(
      "brms_kl: kernel.lag_re.sd must be one number > 0 (lags, or logit units under scale logit)."
    )
  }
  scale <- kl_lag_re_scale(kcfg)
  if (
    !(is.character(scale) &&
      length(scale) == 1L &&
      scale %in% c("month", "logit"))
  ) {
    stop(
      "brms_kl: kernel.lag_re.scale must be 'month' or 'logit' (got '",
      paste(scale, collapse = ","),
      "')."
    )
  }
  list(sd = sd, scale = scale)
}

# Engines (kernel.engine). "stan" (default): the brms program for the kernel
# formula is patched (stage-2 section) into a vectorised program -- the kernel
# weights are computed once per gradient and the term is one N x (NL*J)
# matrix-vector product -- and sampled with cmdstanr; K = 1 is the single
# kernel. "brms": brms::brm() on the non-linear formula itself, which
# recomputes every lag weight row by row (~10 ms per gradient on 1.7k rows):
# a slow cross-check for small fits, single kernel only. At K = 1 the two
# programs have the same log density (tests/test_kl_lp_equivalence.R).
#' Resolved engine: "stan" (patched vectorised program, any K) or "brms"
#' (non-linear brms formula; single kernel only).
kl_engine <- function(kcfg) {
  eng <- kcfg$engine %||% "stan"
  if (!eng %in% c("stan", "brms")) {
    stop("brms_kl: kernel.engine must be 'stan' or 'brms' (got '", eng, "').")
  }
  if (identical(eng, "brms") && kl_is_latent(kcfg)) {
    stop(
      "brms_kl: kernel.engine 'brms' fits a single kernel only; latent classes ",
      "need engine 'stan'."
    )
  }
  if (
    identical(eng, "brms") &&
      (kcfg$form %||% "gaussian") %in% c("dlnm", "discrete")
  ) {
    stop("brms_kl: kernel.form ", kcfg$form, " needs engine 'stan'.")
  }
  if (!(kcfg$form %||% "gaussian") %in% c("gaussian", "dlnm", "discrete")) {
    stop("brms_kl: kernel.form must be gaussian | dlnm | discrete.")
  }
  if (!(kcfg$lag_prior %||% "horseshoe") %in% c("horseshoe", "normal")) {
    stop("brms_kl: kernel.lag_prior must be horseshoe | normal.")
  }
  if (
    isTRUE(kcfg$re_slope$enable) && identical(kcfg$form %||% "gaussian", "dlnm")
  ) {
    stop("brms_kl: kernel.re_slope needs kernel.form discrete or gaussian.")
  }
  if (
    isTRUE(kcfg$re_slope$enable) &&
      isFALSE(kcfg$re_slope[["country"]]) &&
      isFALSE(kcfg$re_slope[["adm1"]])
  ) {
    stop(
      "brms_kl: kernel.re_slope.country and .adm1 are both off; set kernel.re_slope.enable: no instead."
    )
  }
  if (
    !is.null(kcfg[["ksd_fixed"]]) &&
      !identical(kcfg$form %||% "gaussian", "gaussian")
  ) {
    stop(
      "brms_kl: kernel.ksd_fixed is for kernel.form gaussian (discrete: kernel.discrete.widths)."
    )
  }
  if (
    !is.null(kcfg[["ksd_fixed"]]) &&
      !(is.numeric(kcfg[["ksd_fixed"]]) && kcfg[["ksd_fixed"]] > 0)
  ) {
    stop("brms_kl: kernel.ksd_fixed must be a positive number of lags.")
  }
  if (isTRUE(kcfg$spatial$enable) && !kl_is_latent(kcfg)) {
    stop(
      "brms_kl: kernel.spatial needs latent classes (classes >= 2 or prior_grouping)."
    )
  }
  if (kl_lag_re_on(kcfg)) {
    kl_lag_re_cfg(kcfg)
  } # validates kernel.lag_re
  eng
}

# ---------------------------------------------------------------------------
# Stage 1: brms non-linear formula and priors
# ---------------------------------------------------------------------------

#' build_kl_expr
#'
#' Non-linear kernel term as a brms formula string over the J x n_lag basis
#' columns (`meta` = cb_meta entry of the kernel variable, `lags:` mode):
#' `(sum_l w_l * (kb1*c_1l + ... + kbJ*c_Jl)) / (sum_l w_l)`. Each column is
#' mentioned exactly once, lag-major, so brms numbers them C_1.. in
#' meta$cb_col_names order.
build_kl_expr <- function(meta) {
  J <- meta$df
  lags <- meta$lags
  w <- sprintf("exp(-((%d - kmu) / ksd)^2 / 2)", as.integer(lags))
  inner <- vapply(
    seq_along(lags),
    function(i) {
      cols <- meta$col_matrix[, i]
      sprintf(
        "%s * (%s)",
        w[i],
        paste(sprintf("kb%d * %s", seq_len(J), cols), collapse = " + ")
      )
    },
    character(1)
  )
  sprintf(
    "(%s) / (%s)",
    paste(inner, collapse = " + "),
    paste(w, collapse = " + ")
  )
}

#' build_kl_formula
#'
#' brms non-linear formula: y ~ lin + <kernel>, lin ~ 1 + baseline + fixed + RE
#' + offset, kmu ~ 1, ksd ~ 1, kb1 ~ 1 .. kbJ ~ 1.
build_kl_formula <- function(meta, terms, baseline_block) {
  lin_rhs <- paste(
    c("1", baseline_block, terms$fixed_terms, terms$re_terms, terms$off_terms),
    collapse = " + "
  )
  nlp <- kl_nlpars(meta$df)
  subs <- c(
    list(stats::as.formula(paste("lin ~", lin_rhs))),
    lapply(nlp, function(p) stats::as.formula(paste(p, "~ 1")))
  )
  do.call(
    brms::bf,
    c(
      list(stats::as.formula(paste("y ~ lin +", build_kl_expr(meta)))),
      subs,
      list(nl = TRUE)
    )
  )
}

#' Kernel nlpar names: kmu, ksd, kb1..kbJ.
kl_nlpars <- function(J) {
  c(.KL_NLPAR$mu, .KL_NLPAR$sigma, sprintf("%s%d", .KL_NLPAR$beta, seq_len(J)))
}

#' Kernel-parameter prior strings from the resolved kernel config. `kmu`
#' defaults to uniform over the lag range (bounds are set separately).
#' @noRd
.kl_prior_strings <- function(kcfg) {
  list(
    mu = if (is.null(kcfg[["prior_mu"]])) {
      sprintf("uniform(%d, %d)", kcfg$lag_min, kcfg$lag_max)
    } else {
      sprintf(
        "normal(%g, %g)",
        as.numeric(kcfg[["prior_mu"]]),
        kcfg$prior_mu_sd
      )
    },
    sigma = sprintf("normal(0, %g)", kcfg$prior_sigma_sd),
    beta = sprintf("normal(0, %g)", kcfg$prior_beta_sd)
  )
}

#' .kl_coef_priors
#'
#' Per-coefficient overrides of the linear-part priors, from the brms_kl block:
#'   coef_prior_b:  {<coef>: "<prior>"}               class b, nlpar lin
#'   coef_prior_sd: {<group>: {<coef>: "<prior>"}}    class sd, nlpar lin
#' Coef names as brms writes them (factor level appended, e.g.
#' ongoing_outbreak_lag1TRUE). Named so as not to share a prefix with the
#' prior_* keys (`$` partial matching in .prior_strings()). NULL when unset.
.kl_coef_priors <- function(brms_cfg) {
  bc <- brms_cfg[["coef_prior_b"]]
  sc <- brms_cfg[["coef_prior_sd"]]
  ent <- c(
    lapply(names(bc), function(cf) {
      brms::set_prior(bc[[cf]], class = "b", coef = cf, nlpar = "lin")
    }),
    unlist(
      lapply(names(sc), function(g) {
        lapply(names(sc[[g]]), function(cf) {
          brms::set_prior(
            sc[[g]][[cf]],
            class = "sd",
            group = g,
            coef = cf,
            nlpar = "lin"
          )
        })
      }),
      recursive = FALSE
    )
  )
  if (!length(ent)) {
    return(NULL)
  }
  do.call(c, ent)
}

#' build_kl_priors
#'
#' The occurrence priors of the linear part (build_occ_priors(nlpar = "lin"):
#' intercept, b, tiered sd, sds, cor -- identical to brms_cb's) plus the
#' kernel priors: kmu bounded to the lag range, ksd bounded below at
#' kernel.sigma_lb (half-normal), kb_j normal.
build_kl_priors <- function(opt, terms, brms_cfg, kcfg, J) {
  s <- .kl_prior_strings(kcfg)
  c(
    build_occ_priors(opt, terms, brms_cfg, nlpar = "lin"),
    .kl_coef_priors(brms_cfg),
    brms::set_prior(
      s$mu,
      class = "b",
      nlpar = .KL_NLPAR$mu,
      lb = kcfg$lag_min,
      ub = kcfg$lag_max
    ),
    brms::set_prior(
      s$sigma,
      class = "b",
      nlpar = .KL_NLPAR$sigma,
      lb = kcfg$sigma_lb
    ),
    do.call(
      c,
      lapply(seq_len(J), function(j) {
        brms::set_prior(
          s$beta,
          class = "b",
          nlpar = sprintf("%s%d", .KL_NLPAR$beta, j)
        )
      })
    )
  )
}

#' assemble_kl_model
#'
#' Everything brms needs for the kernel model, and the seam stage 2 patches:
#' the dlnm basis (build_cb_block(), `lags:` mode for the kernel variable;
#' any other dlnm variable enters `lin` exactly as in brms_cb via
#' integrate_cb_terms()), the formula terms (build_formula_terms(); region
#' level dropped when `drop_region`), the nl formula and the priors.
#'
#' @param opt            resolved options (after load_surv_data())
#' @param pois_dat       Poisson view (prepare_occ_brms_data())
#' @param brms_cfg       backend block (opt$occurrence$brms_kl %||% $brms)
#' @param baseline_block .occ_baseline_block() output
#' @param kcfg           resolve_kernel_cfg(opt)
#' @return list(formula, priors, data, family, terms, cb_meta, meta, kcfg, opt_f)
assemble_kl_model <- function(
  opt,
  pois_dat,
  brms_cfg,
  baseline_block,
  kcfg = resolve_kernel_cfg(opt)
) {
  dlnm_cfg <- opt$occurrence$dlnm
  cbb <- build_cb_block(pois_dat, dlnm_cfg, " (occurrence lag-kernel model)")
  dat <- cbb$dat

  drop_region <- isTRUE(opt$grouping$enable) || isTRUE(kcfg$drop_region)
  opt_f <- opt
  if (drop_region) {
    opt_f$grouping$enable <- TRUE
  } # formula/prior builders only

  terms <- build_formula_terms(
    opt_f,
    n_country = nlevels(dat$country),
    n_region = length(unique(dat$AFRO_region)),
    dat = dat
  )
  re_sep <- if (isTRUE(brms_cfg$re_uncorrelated)) "||" else "|"
  terms <- integrate_cb_terms(
    terms,
    dlnm_cfg,
    cbb$cb_cols,
    cbb$cb_meta,
    re_sep,
    drop_region = drop_region,
    linear_vars = setdiff(names(cbb$cb_meta), kcfg$variable)
  )
  meta <- cbb$cb_meta[[kcfg$variable]]
  form <- build_kl_formula(meta, terms, baseline_block)

  # brms silently drops rows with NA in any model variable; drop them here so
  # brms's data, the country index of the latent-class program and every
  # summary see the same rows.
  vars <- intersect(
    unique(unlist(lapply(c(list(form$formula), form$pforms), all.vars))),
    names(dat)
  )
  cc <- stats::complete.cases(as.data.frame(dat)[, vars, drop = FALSE])
  if (!all(cc)) {
    message(sprintf(
      "[brms_kl] dropping %d row(s) with NA in model variables (%s), as brms would",
      sum(!cc),
      paste(
        vars[vapply(vars, function(v) anyNA(dat[[v]]), logical(1))],
        collapse = ", "
      )
    ))
    dat <- dat[cc, , drop = FALSE]
  }
  list(
    formula = form,
    priors = build_kl_priors(opt_f, terms, brms_cfg, kcfg, meta$df),
    data = dat,
    family = stats::poisson(),
    terms = terms,
    cb_meta = cbb$cb_meta,
    meta = meta,
    kcfg = kcfg,
    opt_f = opt_f
  )
}

# ---------------------------------------------------------------------------
# Summaries (both stages)
# ---------------------------------------------------------------------------

#' kl_weights
#'
#' Normalised Gaussian lag weights, one row per (kmu, ksd) pair. A width of
#' 0 puts all the weight on the lag nearest the centre (single-lag model).
#' @return matrix length(mu) x length(lags); rows sum to 1
kl_weights <- function(mu, sigma, lags) {
  mu <- as.numeric(mu)
  sigma <- rep_len(as.numeric(sigma), length(mu))
  D <- outer(-mu, as.numeric(lags), `+`)
  W <- exp(-0.5 * D^2 / sigma^2)
  z <- sigma == 0
  if (any(z)) {
    W[z, ] <- t(vapply(
      which(z),
      function(i) {
        w <- numeric(length(lags))
        w[which.min(abs(D[i, ]))] <- 1
        w
      },
      numeric(length(lags))
    ))
  }
  W / rowSums(W)
}

#' Basis of the kernel variable at exposure values `x`, centred at `cen`.
#' @noRd
.kl_basis_centred <- function(meta, x, cen) {
  args <- .onebasis_replay_args(meta$onebasis_attrs)
  B <- unclass(do.call(dlnm::onebasis, c(list(x), args)))
  B0 <- unclass(do.call(dlnm::onebasis, c(list(cen), args)))
  sweep(matrix(B, nrow = length(x)), 2, as.numeric(B0))
}

#' summarise_kernel_draws
#'
#' Posterior summaries of one kernel from its draws: parameters, lag weights
#' and the exposure-response (log-RR of a SUSTAINED exposure x across all
#' lags vs `cen`; since the weights sum to 1 this is sum_j kb_j (B_j(x) -
#' B_j(cen))).
#'
#' @param kmu,ksd numeric vectors (one value per draw)
#' @param beta    matrix draws x J
#' @param meta    cb_meta entry of the kernel variable
#' @param cen     reference exposure
#' @param n_grid  exposure grid size
#' @return list(params, weights, erf) of tibbles (mean, q2.5, q50, q97.5)
summarise_kernel_draws <- function(
  kmu,
  ksd,
  beta,
  meta,
  cen = 0,
  probs = c(0.025, 0.5, 0.975),
  n_grid = 50L
) {
  qs <- function(m) {
    m <- as.matrix(m)
    out <- cbind(
      mean = colMeans(m),
      t(apply(m, 2, stats::quantile, probs = probs, names = FALSE))
    )
    colnames(out) <- c("mean", sprintf("q%s", probs * 100))
    tibble::as_tibble(out)
  }
  beta <- as.matrix(beta)
  params <- dplyr::bind_cols(
    tibble::tibble(par = c("kmu", "ksd", sprintf("kb%d", seq_len(ncol(beta))))),
    qs(cbind(kmu, ksd, beta))
  )
  W <- kl_weights(kmu, ksd, meta$lags)
  weights <- dplyr::bind_cols(tibble::tibble(lag = meta$lags), qs(W))
  at <- seq(meta$range[1], meta$range[2], length.out = n_grid)
  Bc <- .kl_basis_centred(meta, at, cen)
  erf <- dplyr::bind_cols(tibble::tibble(exposure = at), qs(beta %*% t(Bc)))
  list(params = params, weights = weights, erf = erf)
}

#' summarise_kl_fit
#'
#' Stage-1 summary of a brms_kl brmsfit: summarise_kernel_draws() on the
#' b_kmu_Intercept / b_ksd_Intercept / b_kb<j>_Intercept draws.
summarise_kl_fit <- function(fit, meta, cen = 0) {
  J <- meta$df
  vn <- sprintf("b_%s_Intercept", kl_nlpars(J))
  dm <- posterior::as_draws_matrix(fit, variable = vn)
  summarise_kernel_draws(
    kmu = dm[, vn[1]],
    ksd = dm[, vn[2]],
    beta = dm[, vn[-(1:2)], drop = FALSE],
    meta = meta,
    cen = cen
  )
}

# ---------------------------------------------------------------------------
# Stage 2: latent country classes (patched brms Stan program, cmdstanr)
# ---------------------------------------------------------------------------
#
# Model. Country c has a latent class z_c in 1..K; class k owns its kernel
# (kmu[k], ksd[k], kbeta[, k]). The linear part `lin` (baseline, covariates,
# random effects, offset) is shared and kept verbatim from brms. Per country,
#
#   log p(y_c) = log_sum_exp_k( log_softmax(klogit[c])[k]
#                               + sum_{i in c} poisson_log_lpmf(y_i | lin_i + kl_i(k)) )
#   klogit[c, k] = log kpi[k] + kkappa * A[c, k] + kfield[c, k]
#
# A = one-hot prior partition, or its contiguity-softened version (kernel.prior_grouping; kkappa estimated or
# fixed), kfield = spatial field on the country adjacency (kernel.spatial;
# one column per class, the last class is the reference with field 0):
# ICAR (soft sum-to-zero per connected component, island singletons near 0)
# or DAGAR (proper, non-centred; ported from archive/cholera_regime_sim).
# Without a prior partition kmu is `ordered` (labels fixed by lag order);
# with one, class k is anchored to prior group k and kmu is free in the lag
# range. The label-invariant co-membership matrix is reported either way.
#
# Form discrete. Each class's kernel centre is one of the candidate centres
# (the lag grid, every kernel.discrete.centre_step-th lag) and the kernel
# width one of kernel.discrete.widths, shared by the classes (0 = single
# lag). Centres are summed out exactly: the class centres are ORDERED with
# ties allowed (l_1 <= ... <= l_K), so the sum runs over the
# choose(NC + K - 1, K) multisets of centres instead of the NC^K labelled
# tuples, and class k is the k-th shortest lag (labels fixed by
# construction). The countries couple the classes (each country is a mixture
# over all of them), so the sum does not factorise further:
#
#   log p(y | beta, pi) = log_sum_exp_{s, m} ( log p(s) + log p(m)
#        + sum_c log_sum_exp_k( log_softmax(klogit[c])[k] + LL[k][c, q(s, m_k)] ) )
#
# with LL[k][c, q] the country log-likelihood of class k under kernel q
# (width s, centre m_k), precomputed once per gradient (kernel-weighted
# exposure bases are data). log p(m) = log(K! / prod_i n_i!) under
# combo_prior "exchangeable" (n_i = classes at centre i: the prior the
# exchangeable model with iid uniform centres induces on the ordered
# representation, so ties are not favoured) or 0 under "uniform"; plus, when
# kmu priors are set (prior_mu or the prior partition's deployed lags, sorted
# ascending), a discretised normal per position. Generated quantities draw
# (s, m) given the parameters (klag, kwidth) and report the conditional
# P(centre of class k = i) (kl_lag_prob); responsibilities and the country
# log-likelihood are conditional on the drawn (s, m). Two classes tied on a
# centre are exchangeable, so draws are relabelled within ties by kcon.
#
# Program. brms::stancode() of the stage-1 model is patched at landmarks
# that are each asserted to occur exactly once (a brms upgrade that moves
# one fails loudly, naming it). brms's data block and linear predictor are
# untouched; the kernel nlpar parameters, their priors, the scalar non-linear
# predictor and the single likelihood statement are replaced. Threading
# (reduce_sum) is not supported on this path.

# Stan fragments ------------------------------------------------------------

.KL_STAN_FUNCTIONS <- '
  // ---- brms_kl latent classes (utils/dlnm_kernel.R) ----
  // Row log-likelihood of every row under every class: eta = lin + U[tix, k],
  // U = KBu * Theta the kernel term of every distinct exposure row (theta_k =
  // the class per-lag exposure coefficients, lag-major, whatever the form --
  // Gaussian kernel or class lag basis -- that produced them), tix = the
  // distinct exposure row of each row.
  matrix kl_class_rowll(vector nlp_lin, matrix U, array[] int tix,
                        vector Yv, vector lgY1) {
    int N = rows(nlp_lin);
    int K = cols(U);
    matrix[N, K] ll;
    for (k in 1:K) {
      vector[N] eta = nlp_lin + col(U, k)[tix];
      ll[:, k] = Yv .* eta - exp(eta) - lgY1;
    }
    return ll;
  }
  // Non-centred DAGAR field (Datta et al. 2019), one field per column of z:
  // a deterministic map from standard normals, so no Jacobian / log-det term.
  // Ported verbatim from archive/cholera_regime_sim (stan/regime_mixture_nc.stan,
  // dagar_field()). A node with no parent (first node of a component, island)
  // gets sd tau.
  matrix kl_dagar_field(matrix z, real rho, real tau,
                        array[] int n_par, array[] int par_nodes,
                        array[] int par_start) {
    int Ii = rows(z);
    int Uc = cols(z);
    matrix[Ii, Uc] out;
    real rho2 = square(rho);
    for (i in 1:Ii) {
      real denom = 1 + (n_par[i] - 1) * rho2;  // == 1 - rho2 when n_par[i]==0
      real b     = rho / denom;
      real sd_i  = tau * sqrt((1 - rho2) / denom);
      row_vector[Uc] parent_sum = rep_row_vector(0.0, Uc);
      for (p in 1:n_par[i])
        parent_sum += out[par_nodes[par_start[i] + p - 1], ];
      out[i, ] = b * parent_sum + sd_i * z[i, ];
    }
    return out;
  }
  // Membership log-probabilities of one country: log_softmax of its logits, or
  // a point mass on class `fix` when the country is a seed (kernel.pin).
  row_vector kl_member_logp(row_vector logit, int fix) {
    row_vector[cols(logit)] out;
    if (fix == 0) return log_softmax(logit\')\';
    out = rep_row_vector(negative_infinity(), cols(logit));
    out[fix] = 0;
    return out;
  }
  // Sum row log-likelihoods within each country (rows sorted by country).
  matrix kl_country_ll(matrix ll, array[] int cs, array[] int ce) {
    int C = size(cs);
    int K = cols(ll);
    matrix[C, K] out;
    for (c in 1:C) for (k in 1:K) out[c, k] = sum(ll[cs[c]:ce[c], k]);
    return out;
  }
  // Form discrete: country log-likelihood of every class under every kernel
  // q (column q of KBq[j] = exposure basis j weighted over the lags by kernel
  // q). out[k][c, q] = sum_{i in c} log Poisson(y_i | lin_i + KBq[, q] kbeta[, k]
  // + E[i, q]); lgY1c = country sums of lgamma(y + 1); E (used when use_E) =
  // the country/ADM1 deviations of the exposure coefficients times the kernel
  // basis (kl_disc_re_term()), the same for every class.
  array[] matrix kl_disc_country_ll(vector nlp_lin, array[] matrix KBq, matrix kbeta,
                                    vector Yv, vector lgY1c, array[] int cs, array[] int ce,
                                    matrix E, int use_E) {
    int N = rows(nlp_lin);
    int Q = cols(KBq[1]);
    int J = size(KBq);
    int K = cols(kbeta);
    int C = size(cs);
    array[K] matrix[C, Q] out;
    matrix[N, Q] eta0 = rep_matrix(nlp_lin, Q);
    if (use_E) eta0 += E;   // random-slope deviations (class-independent)
    for (k in 1:K) {
      matrix[N, Q] eta = eta0;
      matrix[N, Q] ll;
      for (j in 1:J) eta += kbeta[j, k] * KBq[j];
      ll = diag_pre_multiply(Yv, eta) - exp(eta);
      for (c in 1:C)
        out[k][c] = rep_row_vector(1.0, ce[c] - cs[c] + 1) * ll[cs[c]:ce[c], ] - lgY1c[c];
    }
    return out;
  }
  // Forms gaussian / dlnm, country path (kernel threading or random slopes):
  // class log-probabilities of one country, lsm[k] + sum_i log Poisson(y_i |
  // lin_i + B_i theta_k), or with random slopes lin_i + (B_i T_k)(kbeta_k + D_i),
  // where B_i T_k = the exposure basis of row i weighted over the lags by the
  // class k kernel (T_k maps lag-major basis columns to J exposure columns).
  // The kernel term depends on the row only through its exposure basis B_i,
  // which takes few distinct values (ENSO is one series: one per calendar
  // month), so it is computed once per distinct row t by the caller and
  // indexed here by tix: U[t, k] = B_t theta_k, BW[k][t] = B_t T_k.
  vector kl_ctry_logp(vector lin, array[] int tix, vector y, real lgy, row_vector lsm,
                    matrix U, array[] matrix BW, matrix kbeta, int use_re, matrix D) {
    int K = cols(lsm);
    int n = rows(lin);
    vector[K] lp;
    for (k in 1:K) {
      vector[n] eta;
      if (use_re) {
        matrix[n, rows(kbeta)] Bw = BW[k][tix, ];
        eta = lin + Bw * col(kbeta, k) + rows_dot_product(Bw, D);
      } else {
        eta = lin + col(U, k)[tix];
      }
      lp[k] = lsm[k] + dot_product(y, eta) - sum(exp(eta)) - lgy;
    }
    return lp;
  }
  // reduce_sum partial sum over countries cidx[start..end]: sum of the country
  // mixture log-likelihoods. lin_pad = the linear predictor of each of those
  // countries, padded to the largest one (the sliced argument, so it is not
  // copied per task).
  real kl_partial_ctry(array[] vector lin_pad, int start, int end, array[] int cidx,
                       array[] int tix, vector Yv,
                       vector lgY1c, array[] int cs, array[] int ce, matrix lsm,
                       matrix U, array[] matrix BW, matrix kbeta, int use_re,
                       matrix kre_c, matrix kre_a, array[] int adm1) {
    real tot = 0;
    for (i in start:end) {
      int c = cidx[i];
      int n = ce[c] - cs[c] + 1;
      matrix[use_re ? n : 0, use_re ? rows(kbeta) : 0] D;
      if (use_re)
        D = (rep_matrix(col(kre_c, c), n) + kre_a[:, adm1[cs[c]:ce[c]]])\';
      tot += log_sum_exp(kl_ctry_logp(lin_pad[i - start + 1][1:n], tix[cs[c]:ce[c]],
                                    Yv[cs[c]:ce[c]], lgY1c[c], lsm[c], U, BW, kbeta,
                                    use_re, D));
    }
    return tot;
  }
  // reduce_sum partial sum over rows r0 + start - 1 .. r0 + end - 1 of one large
  // country c, for class k: sum_i [y_i eta_ik - exp(eta_ik)] (lgamma term and
  // membership added by the caller). Large countries are split by rows so one
  // country does not bound the threading (kernel.split_big).
  real kl_partial_rows(array[] real lin_s, int start, int end, int r0, int k, int c,
                       array[] int tix, vector Yv, matrix U, array[] matrix BW, matrix kbeta,
                       int use_re, matrix kre_c, matrix kre_a, array[] int adm1) {
    int a = r0 + start - 1;
    int b = r0 + end - 1;
    int n = end - start + 1;
    vector[n] eta;
    if (use_re) {
      matrix[n, rows(kbeta)] Bw = BW[k][tix[a:b], ];
      matrix[n, rows(kbeta)] D = (rep_matrix(col(kre_c, c), n) + kre_a[:, adm1[a:b]])\';
      eta = to_vector(lin_s) + Bw * col(kbeta, k) + rows_dot_product(Bw, D);
    } else {
      eta = to_vector(lin_s) + col(U, k)[tix[a:b]];
    }
    return dot_product(Yv[a:b], eta) - sum(exp(eta));
  }
  // Pair path (kernel.pair_collapse): the kernel term of a row depends on the
  // row only through its pair p = (distinct exposure row t, ADM1 unit a) --
  // (t, country) without random slopes -- so with g_pk the class k kernel term
  // of pair p, sum_{i in p} [y_i (lin_i + g_pk) - exp(lin_i + g_pk)]
  //   = sum_{i in p} y_i lin_i + g_pk Yp_p - exp(g_pk) S_p,
  // with Yp_p = sum y_i and S_p = sum exp(lin_i) over the pair rows. This
  // returns the class-dependent part over the pairs of country c:
  // sum_p [g_pk Yp_p - exp(g_pk) S_p] for every k (the row part
  // sum_i y_i lin_i - lgamma(y_i + 1) is the same for every class and is added
  // by the caller). g_pk = U[t_p, k], or with random slopes
  // BW[k][t_p] (kbeta_k + kre_c[, c] + kre_a[, a_p]).
  vector kl_pair_logp(vector S, vector Yp, array[] int pt, array[] int pa, int c,
                      matrix U, array[] matrix BW, matrix kbeta, int use_re,
                      matrix kre_c, matrix kre_a) {
    int K = cols(kbeta);
    int np = rows(S);
    vector[K] out;
    matrix[use_re ? np : 0, use_re ? rows(kbeta) : 0] D;
    if (use_re) D = (rep_matrix(col(kre_c, c), np) + kre_a[:, pa])\';
    for (k in 1:K) {
      vector[np] g;
      if (use_re) g = rows_dot_product(BW[k][pt, ], D + rep_matrix(col(kbeta, k)\', np));
      else g = col(U, k)[pt];
      out[k] = dot_product(Yp, g) - dot_product(S, exp(g));
    }
    return out;
  }
  // reduce_sum partial sum over countries cidx[start..end] on the pair path:
  // sum of log_sum_exp_k(lsm[c, k] + kl_pair_logp(...)[k]). S_pad = the pair
  // sums S_p of each of those countries, padded to the largest (sliced).
  real kl_partial_pairs(array[] vector S_pad, int start, int end, array[] int cidx,
                        vector Ypair, array[] int pcs, array[] int pce,
                        array[] int pt, array[] int pa, matrix lsm,
                        matrix U, array[] matrix BW, matrix kbeta, int use_re,
                        matrix kre_c, matrix kre_a) {
    real tot = 0;
    for (i in start:end) {
      int c = cidx[i];
      int np = pce[c] - pcs[c] + 1;
      tot += log_sum_exp(lsm[c]\' + kl_pair_logp(S_pad[i - start + 1][1:np], Ypair[pcs[c]:pce[c]],
                                                pt[pcs[c]:pce[c]], pa[pcs[c]:pce[c]], c,
                                                U, BW, kbeta, use_re, kre_c, kre_a));
    }
    return tot;
  }
  // Form discrete, kernel.re_slope: E[i, q] = sum_j KBq[j][i, q] D[i, j], with D
  // the country + ADM1 deviation of the row for exposure coefficient j.
  matrix kl_disc_re_term(array[] matrix KBq, matrix D) {
    int N = rows(D);
    int Q = cols(KBq[1]);
    matrix[N, Q] E = rep_matrix(0, N, Q);
    for (j in 1:size(KBq)) E += diag_pre_multiply(col(D, j), KBq[j]);
    return E;
  }
  // kernel.re_slope: ADM1 deviations of the ENSO exposure coefficients (J x NA),
  // zero when that tier is off (kernel.re_slope.adm1: false; sda, za empty).
  matrix kl_adm1_dev(vector sda, matrix za, int use, int J, int NA) {
    if (use) return diag_pre_multiply(sda, za);
    return rep_matrix(0, J, NA);
  }
  // Form discrete: log-likelihood of every (width s, centre combination m),
  // classes marginalised per country: element (s - 1) * M + m is
  //   sum_c log_sum_exp_k( lsm[c, k] + LLq[k][c, qidx[s, k, m]] ).
  // Vectorised over combinations: the log_sum_exp over k is shifted by its
  // elementwise max (the shift cancels in the gradient).
  vector kl_disc_combo_ll(array[] matrix LLq, matrix lsm, array[,,] int qidx) {
    int S = dims(qidx)[1];
    int K = dims(qidx)[2];
    int M = dims(qidx)[3];
    int C = rows(lsm);
    vector[S * M] tot;
    for (s in 1:S) {
      array[K] matrix[C, M] G;
      matrix[C, M] mx;
      matrix[C, M] acc = rep_matrix(0, C, M);
      for (k in 1:K) G[k] = LLq[k][:, qidx[s, k]] + rep_matrix(col(lsm, k), M);
      mx = G[1];
      for (k in 2:K) mx = fmax(mx, G[k]);
      for (k in 1:K) acc += exp(G[k] - mx);
      tot[((s - 1) * M + 1):(s * M)] = (rep_row_vector(1.0, C) * (mx + log(acc)))\';
    }
    return tot;
  }'

.KL_STAN_DATA <- '
  // ---- brms_kl latent classes (utils/dlnm_kernel.R) ----
  int<lower=1> kl_K;                                // classes (1 = single kernel)
  int<lower=1> kl_C;                                // countries
  array[N] int<lower=1, upper=kl_C> kl_ctry;        // country of each row (rows sorted)
  array[kl_C] int<lower=1> kl_cs;                   // first row of each country
  array[kl_C] int<lower=1> kl_ce;                   // last row of each country
  int<lower=1> kl_J;                                // exposure basis functions
  int<lower=2> kl_NL;                               // lags
  int<lower=1> kl_T;                                // distinct rows of the kernel exposure basis (ENSO: calendar months)
  array[N] int<lower=1, upper=kl_T> kl_tix;         // distinct exposure-basis row of each row
  int<lower=0, upper=1> kl_pairs;                   // country path summed per pair (kernel.pair_collapse)
  int<lower=1> kl_P;                                // pairs: (exposure row, ADM1) with random slopes, else (exposure row, country)
  array[N] int<lower=1, upper=kl_P> kl_pix;         // pair of each row
  array[kl_P] int<lower=1, upper=kl_T> kl_pt;       // exposure row of each pair
  array[kl_P] int<lower=1> kl_pa;                   // ADM1 unit of each pair (pairs sorted by country)
  array[kl_C] int<lower=1> kl_pcs;                  // first pair of each country
  array[kl_C] int<lower=1> kl_pce;                  // last pair of each country
  int<lower=1> kl_maxp;                             // pairs of the country with the most
  vector[kl_NL] kl_lags;                            // lag grid (ascending)
  int<lower=0, upper=1> kl_dlnm;                    // form: 0 = Gaussian kernel, 1 = class lag basis
  int<lower=0, upper=1> kl_disc;                    // form discrete: grid centres, ordered (ties allowed), summed out
  int<lower=1> kl_S;                                // discrete: kernel widths (shared by the classes, summed out)
  int<lower=1> kl_NCt;                              // discrete: candidate centres
  vector[kl_NCt] kl_cent;                           // discrete: centre values (lags)
  vector[kl_S] kl_width;                            // discrete: width values (lags; 0 = single lag)
  matrix[kl_NL, kl_S * kl_NCt] kl_Wq;               // discrete: lag weights of kernel q = (s - 1) * NCt + centre
  int<lower=1> kl_M;                                // discrete: ordered centre combinations
  array[kl_M, kl_K] int<lower=1, upper=kl_NCt> kl_combo;   // discrete: centre index per class, non-decreasing
  vector[kl_M] kl_combo_lp;                         // discrete: log prior of each combination
  int<lower=0, upper=1> kl_re;                      // discrete: country + ADM1 deviations of the exposure coefficients
  int<lower=1> kl_NA;                               // ADM1 units (country:ADM1)
  array[N] int<lower=1, upper=kl_NA> kl_adm1;       // ADM1 unit of each row
  vector<lower=0>[2] kl_re_sd;                      // half-normal sds of the country / ADM1 deviation sds
  int<lower=0, upper=1> kl_re_ctry;                 // country-level deviations on (0 = ADM1 only)
  int<lower=0, upper=1> kl_re_adm1;                 // ADM1-level deviations on (0 = country only)
  int<lower=0, upper=1> kl_re_sdc_fix;              // country deviation sd fixed (kl_re_sdc_value), not estimated
  real<lower=0> kl_re_sdc_value;
  int<lower=0, upper=1> kl_rs;                      // gaussian/dlnm: country path kl_partial_ctry (reduce_sum)
  int<lower=1> kl_grain;                            // reduce_sum grainsize (countries)
  int<lower=1> kl_maxn;                             // rows of the largest country summed whole
  int<lower=0> kl_nsmall;                           // country path: countries summed whole (country reduce_sum)
  array[kl_nsmall] int<lower=1, upper=kl_C> kl_small;
  int<lower=0> kl_nbig;                             // country path: countries split by rows (one reduce_sum per class)
  array[kl_nbig] int<lower=1, upper=kl_C> kl_big;
  int<lower=1> kl_grain_rows;                       // reduce_sum grainsize (rows) for the split countries
  int<lower=0, upper=1> kl_ksd_fix;                 // gaussian: kernel width fixed (kl_ksd_value)
  real<lower=0> kl_ksd_value;
  int<lower=1> kl_BL;                               // lag-basis functions (form dlnm)
  matrix[kl_NL, kl_BL] kl_Lb;                       // B-spline lag basis, rows = lags (form dlnm)
  int<lower=0, upper=1> kl_hs;                      // dlnm prior: 1 = regularised group horseshoe, 0 = normal
  real<lower=0> kl_hs_scale_global;                 // half-Cauchy scale of the per-class global shrinkage
  real<lower=0> kl_hs_slab_scale;                   // slab scale (caps active coefficients)
  real<lower=0> kl_hs_slab_df;
  real<lower=0> kl_alpha;                           // Dirichlet concentration
  int<lower=0, upper=1> kl_mu_normal;               // normal (1) or uniform (0) kmu prior
  vector[kl_K] kl_mu_mean;
  real<lower=0> kl_mu_sd;
  real<lower=0> kl_sd_lb;                           // lower bound of ksd
  real<lower=0> kl_sd_sd;                           // half-normal sd of ksd
  real<lower=0> kl_beta_sd;
  int<lower=0, upper=1> kl_use_pg;                  // prior-partition anchoring (kappa * A) on?
  int<lower=0, upper=1> kl_mu_ordered;              // identify by ordered kmu (lag centre)
  int<lower=0, upper=1> kl_erf_ordered;             // identify by ordered ERF contrast
  vector[kl_J] kl_erf_v;                            // B(erf_at) - B(cen): ERF contrast direction
  matrix[kl_J, kl_J - 1] kl_erf_Q;                  // orthonormal complement of kl_erf_v
  matrix[kl_C, kl_K] kl_A;                          // one-hot prior membership
  int<lower=0, upper=1> kl_kappa_fixed;
  real<lower=0> kl_kappa_value;
  real<lower=0> kl_kappa_sd;
  int<lower=0, upper=1> kl_use_sp;                  // spatial membership term on?
  int<lower=0> kl_field_cols;                       // field columns: K - 1 (reference class) or K (symmetric)
  int<lower=0> kl_n_edges;
  array[kl_n_edges] int<lower=1, upper=kl_C> kl_node1;
  array[kl_n_edges] int<lower=1, upper=kl_C> kl_node2;
  int<lower=0> kl_n_comp;                           // connected components
  array[kl_C] int<lower=0, upper=kl_n_comp> kl_comp;
  vector[kl_n_comp] kl_comp_n;
  real<lower=0> kl_tau_sd;
  int<lower=0, upper=1> kl_sp_dagar;                // spatial prior: 0 = ICAR, 1 = DAGAR
  int<lower=0> kl_NP;                               // DAGAR parent links
  array[kl_C] int<lower=0> kl_n_par;
  array[kl_NP] int<lower=1, upper=kl_C> kl_par_nodes;
  array[kl_C] int<lower=1> kl_par_start;
  real<lower=0> kl_rho_a;                           // beta prior on the DAGAR rho
  real<lower=0> kl_rho_b;
  int<lower=0, upper=1> kl_gq_row_loglik;
  array[kl_C] int<lower=0, upper=kl_K> kl_fix;      // seed anchoring (kernel.pin): 0 = free, k = country fixed to class k'

.KL_STAN_TDATA_DECL <- '
  // ---- brms_kl latent classes ----
  real kl_lo = kl_lags[1];
  real kl_hi = kl_lags[kl_NL];
  vector[N] kl_Yv = to_vector(Y);
  vector[N] kl_lgY1;
  matrix[N, kl_NL * kl_J] kl_KB;   // exposure basis at every lag, lag-major (brms C_1..C_<NL*J>)
  matrix[kl_T, kl_NL * kl_J] kl_KBu;   // its distinct rows (kl_KB[n] = kl_KBu[kl_tix[n]])
  vector[kl_P] kl_Ypair = rep_vector(0, kl_P);   // pair sums of Y
  array[kl_C] int kl_call;         // 1..C (pair path: every country summed whole)
  vector[kl_C] kl_lgY1c;           // country sums of lgamma(Y + 1)
  array[kl_disc ? kl_J : 0] matrix[N, kl_S * kl_NCt] kl_KBq;   // discrete: kernel-weighted exposure basis
  array[kl_disc ? kl_S : 0, kl_K, kl_M] int kl_qidx;          // discrete: kernel of class k in (s, m)
  vector[kl_disc ? kl_S * kl_M : 0] kl_qlp;                   // discrete: log prior of (s, m)'

.KL_STAN_TDATA_STMT <- '
  for (c in 1:kl_C) kl_lgY1c[c] = sum(kl_lgY1[kl_cs[c]:kl_ce[c]]);
  for (n in 1:N) kl_KBu[kl_tix[n]] = kl_KB[n];
  for (n in 1:N) kl_Ypair[kl_pix[n]] += kl_Yv[n];
  for (c in 1:kl_C) kl_call[c] = c;
  if (kl_disc) {
    for (j in 1:kl_J) {
      matrix[N, kl_NL] bj;
      for (l in 1:kl_NL) bj[:, l] = col(kl_KB, (l - 1) * kl_J + j);
      kl_KBq[j] = bj * kl_Wq;
    }
    for (s in 1:kl_S) for (k in 1:kl_K) for (m in 1:kl_M)
      kl_qidx[s, k, m] = (s - 1) * kl_NCt + kl_combo[m, k];
    for (s in 1:kl_S) kl_qlp[((s - 1) * kl_M + 1):(s * kl_M)] = kl_combo_lp - log(kl_S);
  }'

.KL_STAN_PARAMS <- '
  // ---- brms_kl latent classes ----
  // Gaussian-kernel form (kl_dlnm = 0; form discrete keeps only kbeta_free)
  ordered[(kl_mu_ordered && !kl_dlnm) ? kl_K : 0] kmu_ord;   // identify = kmu: lag order fixes labels
  vector<lower=kl_lo, upper=kl_hi>[(kl_mu_ordered || kl_dlnm || kl_disc) ? 0 : kl_K] kmu_free;
  vector<lower=kl_sd_lb>[(kl_dlnm || kl_disc || kl_ksd_fix) ? 0 : kl_K] ksd;
  matrix[kl_J, (kl_erf_ordered || kl_dlnm) ? 0 : kl_K] kbeta_free;
  ordered[kl_erf_ordered ? kl_K : 0] kc;            // identify = erf: log RR at erf_at vs cen, increasing
  matrix[kl_J - 1, kl_erf_ordered ? kl_K : 0] keta; // free basis directions orthogonal to the contrast
  // Class lag-basis form (kl_dlnm = 1): non-centred coefficients, column-major
  // J x BL per class (exposure basis j within lag-basis function b)
  matrix[kl_dlnm ? kl_J * kl_BL : 0, kl_K] kz;
  vector<lower=0>[(kl_dlnm && kl_hs) ? kl_K : 0] khs_tau;              // per-class global shrinkage
  matrix<lower=0>[(kl_dlnm && kl_hs) ? kl_BL : 0, kl_K] khs_lambda;    // per lag window (group)
  vector<lower=0>[(kl_dlnm && kl_hs) ? kl_K : 0] khs_c2;               // slab
  // Random slopes on the exposure coefficients (form discrete, kernel.re_slope)
  matrix[(kl_re && kl_re_ctry) ? kl_J : 0, kl_C] kre_zc;   // country deviations, standardised
  matrix[(kl_re && kl_re_adm1) ? kl_J : 0, kl_NA] kre_za;   // ADM1 deviations, standardised
  vector<lower=0>[(kl_re && kl_re_ctry && !kl_re_sdc_fix) ? kl_J : 0] kre_sdc;   // country sd per exposure coefficient
  vector<lower=0>[(kl_re && kl_re_adm1) ? kl_J : 0] kre_sda;   // ADM1 sd per exposure coefficient
  // Membership
  simplex[kl_K] kpi;
  vector<lower=0>[(kl_use_pg && !kl_kappa_fixed) ? 1 : 0] kkappa_free;
  matrix[kl_use_sp ? kl_C : 0, kl_field_cols] kphi; // ICAR field / DAGAR innovations
  vector<lower=0>[kl_use_sp ? 1 : 0] ktau;          // spatial scale
  vector<lower=0, upper=1>[(kl_use_sp && kl_sp_dagar) ? 1 : 0] krho;   // DAGAR dependence'

.KL_STAN_TPARAMS_DECL <- '
  // ---- brms_kl latent classes ----
  vector[(kl_dlnm || kl_disc) ? 0 : kl_K] kmu;      // Gaussian form: class lag centres
  matrix[(kl_dlnm || kl_disc) ? 0 : kl_NL, kl_K] kw;   // Gaussian form: class lag weights
  vector[(kl_dlnm || kl_disc) ? 0 : kl_K] kksd;     // Gaussian form: kernel widths (fixed or ksd)
  array[(kl_rs && kl_re && !kl_dlnm && !kl_disc) ? kl_K : 0] matrix[kl_NL * kl_J, kl_J] kTw;   // lag-weighting maps
  matrix[kl_J, kl_K] kbeta;                         // sustained-exposure coefficients (sum over lags)
  matrix[kl_disc ? 0 : kl_NL * kl_J, kl_K] ktheta;  // per-lag exposure coefficients, lag-major
  vector[kl_K] kcon;                                // log RR of sustained ENSO erf_at vs cen
  real kkappa;                                      // prior-partition strength
  matrix[kl_C, kl_K] klogit;                        // membership logits
  matrix[kl_re ? kl_J : 0, kl_C] kre_c;             // country deviations of the exposure coefficients
  matrix[kl_use_sp ? kl_C : 0, kl_field_cols] kfield;   // spatial membership field (scaled)'

.KL_STAN_TPARAMS_STMT <- '
  // ---- brms_kl latent classes ----
  if (kl_disc) {
    // Discrete centres: the lag weights are data (kl_Wq); only the exposure
    // coefficients are parameters (centres and width are summed out).
    kbeta = kbeta_free;
  } else if (!kl_dlnm) {
    // Gaussian kernel: theta at lag l = w_l(kmu, ksd) * kbeta (weights sum to 1)
    if (kl_mu_ordered) {
      kmu = kl_lo + (kl_hi - kl_lo) * inv_logit(kmu_ord);
      lprior += sum(log(kl_hi - kl_lo) + log_inv_logit(kmu_ord) + log1m_inv_logit(kmu_ord));  // Jacobian
    } else {
      kmu = kmu_free;
    }
    // kbeta: free, or the linear (constant-Jacobian) map from the ordered ERF
    // contrast kc and the orthogonal directions keta, so dot(v, kbeta[, k]) = kc[k].
    if (kl_erf_ordered) {
      for (k in 1:kl_K)
        kbeta[:, k] = kc[k] * kl_erf_v / dot_self(kl_erf_v) + kl_erf_Q * keta[:, k];
    } else {
      kbeta = kbeta_free;
    }
    if (kl_ksd_fix) kksd = rep_vector(kl_ksd_value, kl_K);
    else kksd = ksd;
    for (k in 1:kl_K) {
      kw[:, k] = softmax(-0.5 * square((kl_lags - kmu[k]) / kksd[k]));
      for (l in 1:kl_NL) ktheta[((l - 1) * kl_J + 1):(l * kl_J), k] = kw[l, k] * kbeta[:, k];
    }
    if (kl_rs && kl_re)
      for (k in 1:kl_K) {
        kTw[k] = rep_matrix(0, kl_NL * kl_J, kl_J);
        for (l in 1:kl_NL) for (j in 1:kl_J) kTw[k][(l - 1) * kl_J + j, j] = kw[l, k];
      }
  } else {
    // Class lag basis: J x NL coefficient surface = (J x BL coefficients) * Lb\'.
    // Regularised group horseshoe (Piironen & Vehtari 2017): one local scale per
    // lag-basis function b, shared by its J exposure coefficients, so a lag
    // window is either active or shrunk to zero as a whole.
    for (k in 1:kl_K) {
      matrix[kl_J, kl_BL] th;
      matrix[kl_J, kl_NL] bl;
      if (kl_hs) {
        vector[kl_BL] lt = sqrt(khs_c2[k] * square(khs_lambda[:, k])
                                ./ (khs_c2[k] + square(khs_tau[k]) * square(khs_lambda[:, k])));
        th = to_matrix(kz[:, k], kl_J, kl_BL) .* rep_matrix((khs_tau[k] * lt)\', kl_J);
      } else {
        th = kl_beta_sd * to_matrix(kz[:, k], kl_J, kl_BL);
      }
      bl = th * kl_Lb\';
      ktheta[:, k] = to_vector(bl);
      kbeta[:, k]  = bl * rep_vector(1.0, kl_NL);
    }
  }
  for (k in 1:kl_K) kcon[k] = dot_product(kl_erf_v, kbeta[:, k]);
  if (kl_re) {
    if (!kl_re_ctry) kre_c = rep_matrix(0, kl_J, kl_C);
    else if (kl_re_sdc_fix) kre_c = kl_re_sdc_value * kre_zc;
    else kre_c = diag_pre_multiply(kre_sdc, kre_zc);
  }
  kkappa = 0;
  if (kl_use_pg) kkappa = kl_kappa_fixed ? kl_kappa_value : kkappa_free[1];
  if (kl_use_sp) {
    if (kl_sp_dagar)
      kfield = kl_dagar_field(kphi, krho[1], ktau[1], kl_n_par, kl_par_nodes, kl_par_start);
    else
      kfield = ktau[1] * kphi;
  }
  for (c in 1:kl_C) {
    klogit[c] = log(kpi)\';
    if (kl_use_pg) klogit[c] = klogit[c] + kkappa * kl_A[c];
    if (kl_use_sp && kl_field_cols > 0)
      klogit[c, 1:kl_field_cols] = klogit[c, 1:kl_field_cols] + kfield[c];
  }
  if (kl_re) {
    if (kl_re_ctry && !kl_re_sdc_fix)
      lprior += normal_lpdf(kre_sdc | 0, kl_re_sd[1]) - kl_J * normal_lccdf(0 | 0, kl_re_sd[1]);
    if (kl_re_adm1)
      lprior += normal_lpdf(kre_sda | 0, kl_re_sd[2]) - kl_J * normal_lccdf(0 | 0, kl_re_sd[2]);
  }
  if (kl_disc) {
    lprior += normal_lpdf(to_vector(kbeta) | 0, kl_beta_sd);   // centre priors: kl_combo_lp
  } else if (!kl_dlnm) {
    if (kl_mu_normal) lprior += normal_lpdf(kmu | kl_mu_mean, kl_mu_sd);
    else lprior += -kl_K * log(kl_hi - kl_lo);      // uniform on the lag range (as brms)
    if (!kl_ksd_fix)
      lprior += normal_lpdf(ksd | 0, kl_sd_sd) - kl_K * normal_lccdf(kl_sd_lb | 0, kl_sd_sd);
    lprior += normal_lpdf(to_vector(kbeta) | 0, kl_beta_sd);
  } else {
    lprior += std_normal_lpdf(to_vector(kz));
    if (kl_hs) {
      lprior += student_t_lpdf(khs_tau | 1, 0, kl_hs_scale_global)
                - kl_K * student_t_lccdf(0 | 1, 0, kl_hs_scale_global);           // half-Cauchy
      lprior += student_t_lpdf(to_vector(khs_lambda) | 1, 0, 1)
                - kl_K * kl_BL * student_t_lccdf(0 | 1, 0, 1);                    // half-Cauchy
      lprior += inv_gamma_lpdf(khs_c2 | 0.5 * kl_hs_slab_df,
                               0.5 * kl_hs_slab_df * square(kl_hs_slab_scale));
    }
  }
  lprior += dirichlet_lpdf(kpi | rep_vector(kl_alpha, kl_K));
  if (kl_use_pg && !kl_kappa_fixed)
    lprior += normal_lpdf(kkappa_free | 0, kl_kappa_sd) - normal_lccdf(0 | 0, kl_kappa_sd);
  if (kl_use_sp)
    lprior += normal_lpdf(ktau | 0, kl_tau_sd) - normal_lccdf(0 | 0, kl_tau_sd);
  if (kl_use_sp && kl_sp_dagar)
    lprior += beta_lpdf(krho | kl_rho_a, kl_rho_b);'

# identify = repel only (patch_kl_latent_class_stancode(repel = TRUE)): the
# default program does not carry these, so it stays byte-identical.
.KL_STAN_REPEL_DATA <- '
  // ---- identify = repel (utils/dlnm_kernel.R::kl_repel_cfg()) ----
  int<lower=1> kl_rep_NX;                           // exposure-response grid points
  matrix[kl_rep_NX, kl_J] kl_rep_V;                 // exposure basis at the grid minus at cen
  real<lower=0> kl_rep_lag_scale;                   // lag-centre distance scale (lags)
  real<lower=0> kl_rep_erf_scale;                   // response distance scale (log RR, rms over the grid)
  int<lower=0> kl_rep_G;                            // registry groups of the membership pull (0 = none)
  matrix[kl_C, kl_rep_G] kl_rep_A;                  // registry membership of each country'

.KL_STAN_REPEL_PARAMS <- '
  matrix[kl_rep_G, kl_K] kgrp;                      // identify = repel: group-level class logits, standardised'

.KL_STAN_REPEL_TP <- '
  // ---- identify = repel: label-free registry pull, pairwise repulsion ----
  // kl_use_pg = 1 with kl_A = 0, so the kappa * A term above is zero and kkappa
  // (estimated or fixed) is the sd of the group-level class logits.
  if (kl_rep_G > 0) {
    matrix[kl_C, kl_K] kpull = kl_rep_A * kgrp;
    for (c in 1:kl_C) klogit[c] = klogit[c] + kkappa * kpull[c];
    lprior += std_normal_lpdf(to_vector(kgrp));
  }
  {
    matrix[kl_rep_NX, kl_K] kerf = kl_rep_V * kbeta;   // log RR at the grid, per class
    for (j in 1:(kl_K - 1)) for (k in (j + 1):kl_K)
      lprior += log1m_exp(-0.5 * (square((kmu[j] - kmu[k]) / kl_rep_lag_scale)
                                  + dot_self(col(kerf, j) - col(kerf, k))
                                    / (kl_rep_NX * square(kl_rep_erf_scale))));
  }'

# kernel.lag_re only (patch_kl_latent_class_stancode(lag_re = TRUE)): a lag
# offset per country, shared by the classes. The kernel term of country c under
# class k is computed on the country-stacked distinct exposure rows
# (c - 1) * kl_T + t, and the row / pair indices kl_tixc / kl_ptc point into
# that stack, so the likelihood functions are the shared ones, called with the
# stacked U / BW instead of the per-class ones (.kl_lag_re_likelihood()). The
# default program does not carry these, so it stays byte-identical.
.KL_STAN_LAG_RE_FUNCTIONS <- '
  // ---- kernel.lag_re: country lag offsets (utils/dlnm_kernel.R) ----
  // Lag-weighted exposure basis of every class on the country-stacked distinct
  // exposure rows: BW[k][(c - 1) * T + t] = sum_l w_ckl KBl[l][t], with
  // w_ck = softmax(-0.5 ((lags - kmu_k - d_c) / ksd_k)^2) (KBl[l] = the J
  // exposure-basis columns of lag l).
  array[] matrix kl_lag_re_bw(array[] matrix KBl, vector lags, vector kmu, vector kksd, vector d) {
    int NL = size(KBl);
    int T = rows(KBl[1]);
    int J = cols(KBl[1]);
    int K = rows(kmu);
    int C = rows(d);
    array[K] matrix[C * T, J] BW;
    for (k in 1:K) {
      for (c in 1:C) {
        vector[NL] w = softmax(-0.5 * square((lags - kmu[k] - d[c]) / kksd[k]));
        matrix[T, J] acc = w[1] * KBl[1];
        for (l in 2:NL) acc += w[l] * KBl[l];
        BW[k][((c - 1) * T + 1):(c * T)] = acc;
      }
    }
    return BW;
  }
  // Kernel term without random slopes: U[(c - 1) * T + t, k] = BW[k][.] kbeta[, k].
  matrix kl_lag_re_u(array[] matrix BW, matrix kbeta) {
    int K = cols(kbeta);
    matrix[rows(BW[1]), K] U;
    for (k in 1:K) U[:, k] = BW[k] * col(kbeta, k);
    return U;
  }'

.KL_STAN_LAG_RE_DATA <- '
  // ---- kernel.lag_re (utils/dlnm_kernel.R::kl_lag_re_cfg()) ----
  real<lower=0> kl_lag_re_sd;                       // half-normal sd of the country lag-offset sd (lags)'

.KL_STAN_LAG_RE_TDATA_DECL <- '
  // ---- kernel.lag_re ----
  array[kl_NL] matrix[kl_T, kl_J] kl_KBl;           // exposure basis of each lag (J columns of kl_KBu)
  array[N] int kl_tixc;                             // row -> country-stacked distinct exposure row
  array[kl_P] int kl_ptc;                           // pair -> country-stacked distinct exposure row'

.KL_STAN_LAG_RE_TDATA_STMT <- '
  // ---- kernel.lag_re ----
  for (l in 1:kl_NL) kl_KBl[l] = kl_KBu[:, ((l - 1) * kl_J + 1):(l * kl_J)];
  for (n in 1:N) kl_tixc[n] = kl_tix[n] + (kl_ctry[n] - 1) * kl_T;
  for (c in 1:kl_C) for (q in kl_pcs[c]:kl_pce[c]) kl_ptc[q] = kl_pt[q] + (c - 1) * kl_T;'

.KL_STAN_LAG_RE_PARAMS <- '
  vector[kl_C] klag_z;                              // kernel.lag_re: country lag offsets, standardised
  vector<lower=0>[1] klag_sd;                       // kernel.lag_re: their sd (lags)'

.KL_STAN_LAG_RE_TPARAMS_DECL <- '
  vector[kl_C] kclag;                               // kernel.lag_re: country lag offsets (lags)'

.KL_STAN_LAG_RE_TP <- '
  // ---- kernel.lag_re: country lag offsets (non-centred) ----
  kclag = klag_sd[1] * klag_z;
  lprior += normal_lpdf(klag_sd | 0, kl_lag_re_sd) - normal_lccdf(0 | 0, kl_lag_re_sd);'

.KL_STAN_LAG_RE_MODEL <- '
  // ---- kernel.lag_re: standardised country lag offsets ----
  target += std_normal_lpdf(klag_z);'

# kernel.lag_re.scale logit only, in place of / next to the month fragments
# above (the month program is unchanged): the offset moves the class centre on
# the logit scale of the lag window, so every country centre stays inside it.
# The centres are a transformed parameter (kccent[k, c], lags) read by the
# kernel term and by the R summaries.
.KL_STAN_LAG_RE_FUNCTIONS_LOGIT <- '
  // ---- kernel.lag_re (scale logit): country lag centres (utils/dlnm_kernel.R) ----
  // As kl_lag_re_bw(), with the centre of country c under class k given:
  // w_ck = softmax(-0.5 ((lags - cent[k, c]) / ksd_k)^2).
  array[] matrix kl_lag_re_bw_cent(array[] matrix KBl, vector lags, matrix cent, vector kksd) {
    int NL = size(KBl);
    int T = rows(KBl[1]);
    int J = cols(KBl[1]);
    int K = rows(cent);
    int C = cols(cent);
    array[K] matrix[C * T, J] BW;
    for (k in 1:K) {
      for (c in 1:C) {
        vector[NL] w = softmax(-0.5 * square((lags - cent[k, c]) / kksd[k]));
        matrix[T, J] acc = w[1] * KBl[1];
        for (l in 2:NL) acc += w[l] * KBl[l];
        BW[k][((c - 1) * T + 1):(c * T)] = acc;
      }
    }
    return BW;
  }
  // Kernel term without random slopes: U[(c - 1) * T + t, k] = BW[k][.] kbeta[, k].
  matrix kl_lag_re_u(array[] matrix BW, matrix kbeta) {
    int K = cols(kbeta);
    matrix[rows(BW[1]), K] U;
    for (k in 1:K) U[:, k] = BW[k] * col(kbeta, k);
    return U;
  }'

.KL_STAN_LAG_RE_DATA_LOGIT <- '
  // ---- kernel.lag_re (utils/dlnm_kernel.R::kl_lag_re_cfg()) ----
  real<lower=0> kl_lag_re_sd;                       // half-normal sd of the country lag-offset sd (logit units of the lag window)'

.KL_STAN_LAG_RE_PARAMS_LOGIT <- '
  vector[kl_C] klag_z;                              // kernel.lag_re: country lag offsets, standardised
  vector<lower=0>[1] klag_sd;                       // kernel.lag_re: their sd (logit units)'

.KL_STAN_LAG_RE_TPARAMS_DECL_LOGIT <- '
  vector[kl_C] kclag;                               // kernel.lag_re: country lag offsets (logit units)
  matrix[kl_K, kl_C] kccent;                        // kernel.lag_re: lag centre of country c under class k (lags)'

.KL_STAN_LAG_RE_TP_LOGIT <- '
  // ---- kernel.lag_re (scale logit): country lag offsets (non-centred) ----
  kclag = klag_sd[1] * klag_z;
  lprior += normal_lpdf(klag_sd | 0, kl_lag_re_sd) - normal_lccdf(0 | 0, kl_lag_re_sd);
  for (k in 1:kl_K) {
    // class centre as a share of the lag window, kept off 0 / 1 (logit finite)
    real pk = fmin(fmax((kmu[k] - kl_lo) / (kl_hi - kl_lo), 1e-9), 1 - 1e-9);
    kccent[k] = to_row_vector(kl_lo + (kl_hi - kl_lo) * inv_logit(logit(pk) + kclag));
  }'

#' The lag_re Stan fragments of a scale (month: the original ones).
#' @noRd
.kl_lag_re_fragments <- function(scale = "month") {
  if (identical(scale, "logit")) {
    list(
      functions = .KL_STAN_LAG_RE_FUNCTIONS_LOGIT,
      data = .KL_STAN_LAG_RE_DATA_LOGIT,
      params = .KL_STAN_LAG_RE_PARAMS_LOGIT,
      tparams_decl = .KL_STAN_LAG_RE_TPARAMS_DECL_LOGIT,
      tp = .KL_STAN_LAG_RE_TP_LOGIT
    )
  } else {
    list(
      functions = .KL_STAN_LAG_RE_FUNCTIONS,
      data = .KL_STAN_LAG_RE_DATA,
      params = .KL_STAN_LAG_RE_PARAMS,
      tparams_decl = .KL_STAN_LAG_RE_TPARAMS_DECL,
      tp = .KL_STAN_LAG_RE_TP
    )
  }
}

#' The likelihood / generated-quantities text of the lag_re variant: the shared
#' text with the per-class kernel terms (kl_U, kl_BW over kl_T distinct rows)
#' replaced by their country-stacked versions and the row / pair indices
#' kl_tix / kl_pt by kl_tixc / kl_ptc. Every pattern must match exactly the
#' expected number of times, so a change of the shared text fails loudly.
#' @param text .KL_STAN_LIKELIHOOD or the formatted .KL_STAN_GQ_STMT
#' @param where "model" or "gq" (the expected counts differ)
#' @param scale kernel.lag_re.scale: "month" (kl_lag_re_bw on kmu + kclag) or
#'   "logit" (kl_lag_re_bw_cent on the centres kccent)
#' @noRd
.kl_lag_re_likelihood <- function(
  text,
  where = c("model", "gq"),
  scale = "month"
) {
  where <- match.arg(where)
  bw <- if (identical(scale, "logit")) {
    "kl_lag_re_bw_cent(kl_KBl, kl_lags, kccent, kksd)"
  } else {
    "kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag)"
  }
  subs <- list(
    list(
      "matrix[kl_re ? 0 : kl_T, kl_re ? 0 : kl_K] kl_U;",
      "matrix[kl_re ? 0 : kl_C * kl_T, kl_re ? 0 : kl_K] kl_U;",
      c(model = 2L, gq = 1L)
    ),
    list(
      "array[kl_re ? kl_K : 0] matrix[kl_T, kl_J] kl_BW;",
      "array[kl_re ? kl_K : 0] matrix[kl_C * kl_T, kl_J] kl_BW;",
      c(model = 2L, gq = 1L)
    ),
    list(
      "for (k in 1:kl_K) kl_BW[k] = kl_KBu * kTw[k];",
      paste0("kl_BW = ", bw, ";"),
      c(model = 2L, gq = 1L)
    ),
    list(
      "kl_U = kl_KBu * ktheta;",
      paste0("kl_U = kl_lag_re_u(", bw, ", kbeta);"),
      c(model = 2L, gq = 1L)
    ),
    list(
      "kl_class_rowll(nlp_lin, kl_KBu * ktheta, kl_tix, kl_Yv, kl_lgY1)",
      paste0(
        "kl_class_rowll(nlp_lin, kl_lag_re_u(",
        bw,
        ", kbeta), kl_tixc, kl_Yv, kl_lgY1)"
      ),
      c(model = 1L, gq = 1L)
    ),
    # pair path: reduce_sum (model), per country (gq)
    list(
      "kl_pt, kl_pa, kl_lsm, kl_U, kl_BW",
      "kl_ptc, kl_pa, kl_lsm, kl_U, kl_BW",
      c(model = 1L, gq = 0L)
    ),
    list(
      "kl_pt[kl_pcs[c]:kl_pce[c]]",
      "kl_ptc[kl_pcs[c]:kl_pce[c]]",
      c(model = 0L, gq = 1L)
    ),
    # country path: small countries, large countries split by rows (model); per country (gq)
    list(
      "kl_small, kl_tix, kl_Yv",
      "kl_small, kl_tixc, kl_Yv",
      c(model = 1L, gq = 0L)
    ),
    list(
      "kl_tix, kl_Yv, kl_U, kl_BW",
      "kl_tixc, kl_Yv, kl_U, kl_BW",
      c(model = 1L, gq = 0L)
    ),
    list(
      "kl_tix[kl_cs[c]:kl_ce[c]]",
      "kl_tixc[kl_cs[c]:kl_ce[c]]",
      c(model = 0L, gq = 1L)
    )
  )
  for (x in subs) {
    n <- lengths(regmatches(text, gregexpr(x[[1]], text, fixed = TRUE)))
    if (n != x[[3]][[where]]) {
      stop(
        ".kl_lag_re_likelihood (",
        where,
        "): '",
        x[[1]],
        "' found ",
        n,
        " time(s), expected ",
        x[[3]][[where]],
        " -- the shared likelihood text changed; update the lag_re patterns.",
        call. = FALSE
      )
    }
    if (n) text <- gsub(x[[1]], x[[2]], text, fixed = TRUE)
  }
  left <- lengths(regmatches(
    text,
    gregexpr("kl_tix[^c]|kl_pt[^c]|kl_KBu \\* ", text)
  ))
  if (left) {
    stop(
      ".kl_lag_re_likelihood (",
      where,
      "): ",
      left,
      " per-class kernel term(s) / index use(s) left unreplaced.",
      call. = FALSE
    )
  }
  text
}

.KL_STAN_LIKELIHOOD <- '
    // ---- brms_kl latent classes: country-level mixture ----
    if (kl_disc) {
      // discrete centres: classes per country, then (width, combination) summed out
      matrix[kl_C, kl_K] kl_lsm;
      matrix[kl_re ? N : 0, kl_re ? kl_S * kl_NCt : 0] kl_E;
      if (kl_re)
        kl_E = kl_disc_re_term(kl_KBq, (kre_c[:, kl_ctry] + kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA)[:, kl_adm1])\');
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      target += log_sum_exp(kl_qlp + kl_disc_combo_ll(
        kl_disc_country_ll(nlp_lin, kl_KBq, kbeta, kl_Yv, kl_lgY1c, kl_cs, kl_ce, kl_E, kl_re),
        kl_lsm, kl_qidx));
    } else if (kl_rs && kl_pairs) {
      // pair path (kernel.pair_collapse): the row part once, the class part per
      // (exposure row, ADM1) pair, reduce_sum over countries (no row split)
      matrix[kl_C, kl_K] kl_lsm;
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      matrix[kl_re ? 0 : kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_T, kl_J] kl_BW;
      vector[kl_P] kl_esum = rep_vector(0, kl_P);
      array[kl_C] vector[kl_maxp] kl_esum_pad;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        for (k in 1:kl_K) kl_BW[k] = kl_KBu * kTw[k];
      } else {
        kl_U = kl_KBu * ktheta;
      }
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      {
        vector[N] kl_el = exp(nlp_lin);
        for (n in 1:N) kl_esum[kl_pix[n]] += kl_el[n];
      }
      for (c in 1:kl_C) {
        kl_esum_pad[c] = rep_vector(0, kl_maxp);
        kl_esum_pad[c][1:(kl_pce[c] - kl_pcs[c] + 1)] = kl_esum[kl_pcs[c]:kl_pce[c]];
      }
      target += dot_product(kl_Yv, nlp_lin) - sum(kl_lgY1c);
      target += reduce_sum(kl_partial_pairs, kl_esum_pad, kl_grain, kl_call, kl_Ypair, kl_pcs, kl_pce,
                           kl_pt, kl_pa, kl_lsm, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a);
    } else if (kl_rs) {
      // country path: reduce_sum over countries (threaded when compiled with
      // STAN_THREADS), random slopes on the exposure coefficients when kl_re
      matrix[kl_C, kl_K] kl_lsm;
      array[kl_nsmall] vector[kl_maxn] kl_lin_pad;
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      // kernel term once per distinct exposure row, not per row
      matrix[kl_re ? 0 : kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_T, kl_J] kl_BW;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        for (k in 1:kl_K) kl_BW[k] = kl_KBu * kTw[k];
      } else {
        kl_U = kl_KBu * ktheta;
      }
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      for (i in 1:kl_nsmall) {
        int c = kl_small[i];
        kl_lin_pad[i] = rep_vector(0, kl_maxn);
        kl_lin_pad[i][1:(kl_ce[c] - kl_cs[c] + 1)] = nlp_lin[kl_cs[c]:kl_ce[c]];
      }
      if (kl_nsmall > 0)
        target += reduce_sum(kl_partial_ctry, kl_lin_pad, kl_grain, kl_small, kl_tix, kl_Yv, kl_lgY1c,
                             kl_cs, kl_ce, kl_lsm, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a, kl_adm1);
      // large countries: each class log-likelihood summed over row blocks in
      // parallel, then the country mixture
      for (b in 1:kl_nbig) {
        int c = kl_big[b];
        array[kl_ce[c] - kl_cs[c] + 1] real kl_lin_s = to_array_1d(nlp_lin[kl_cs[c]:kl_ce[c]]);
        vector[kl_K] kl_lpk;
        for (k in 1:kl_K)
          kl_lpk[k] = kl_lsm[c, k] - kl_lgY1c[c]
                      + reduce_sum(kl_partial_rows, kl_lin_s, kl_grain_rows, kl_cs[c], k, c,
                                   kl_tix, kl_Yv, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a, kl_adm1);
        target += log_sum_exp(kl_lpk);
      }
    } else {
      matrix[kl_C, kl_K] kl_llc = kl_country_ll(
        kl_class_rowll(nlp_lin, kl_KBu * ktheta, kl_tix, kl_Yv, kl_lgY1),
        kl_cs, kl_ce);
      for (c in 1:kl_C) target += log_sum_exp(kl_member_logp(klogit[c], kl_fix[c])\' + kl_llc[c]\');
    }'

.KL_STAN_ICAR <- '
  // ---- brms_kl latent classes: random-slope deviations (non-centred) ----
  if (kl_re) target += std_normal_lpdf(to_vector(kre_zc)) + std_normal_lpdf(to_vector(kre_za));
  // ---- brms_kl latent classes: spatial membership prior ----
  if (kl_use_sp && kl_sp_dagar) target += std_normal_lpdf(to_vector(kphi));   // DAGAR innovations
  if (kl_use_sp && !kl_sp_dagar) {                                             // ICAR
    for (k in 1:kl_field_cols) {
      vector[kl_n_comp] csum = rep_vector(0, kl_n_comp);
      target += -0.5 * dot_self(col(kphi, k)[kl_node1] - col(kphi, k)[kl_node2]);
      for (c in 1:kl_C) csum[kl_comp[c]] += kphi[c, k];
      // Soft sum-to-zero per connected component (island singletons -> ~0).
      // sd 0.01 * size, not the usual 0.001: with a mixture likelihood the field
      // is weakly informed and 0.001 made every transition hit max treedepth
      // (synthetic case D); 0.01 removed that and moves the field mean by ~0.01.
      target += normal_lpdf(csum | 0, 0.01 * kl_comp_n);
    }
  }'

.KL_STAN_GQ_DECL <- '
  // ---- brms_kl latent classes ----
  matrix[kl_K, kl_C] kl_log_resp;               // log class responsibilities
  vector[kl_C] kl_log_lik_ctry;                 // country log-likelihood (LOCO unit)
  vector[kl_gq_row_loglik ? N : 0] log_lik;     // optional row-level log-lik
  vector[kl_disc ? kl_K : 0] klag;              // discrete: drawn class centres (non-decreasing)
  vector[kl_disc ? 1 : 0] kwidth;               // discrete: drawn kernel width
  matrix[kl_disc ? kl_K : 0, kl_NCt] kl_lag_prob;   // discrete: P(centre of class k = i | parameters)'

# GQ statements: %s = brms's linear-predictor statements (copied from model).
.KL_STAN_GQ_STMT <- '
  // ---- brms_kl latent classes ----
  {
    vector[N] nlp_lin = rep_vector(0.0, N);
    matrix[N, kl_K] kl_ll;
    matrix[kl_C, kl_K] kl_llc;
%s
    if (kl_disc) {
      // draw (width, combination) given the parameters; class quantities below
      // are conditional on the draw
      matrix[kl_re ? N : 0, kl_re ? kl_J : 0] kl_D;
      matrix[kl_re ? N : 0, kl_re ? kl_S * kl_NCt : 0] kl_E;
      array[kl_K] matrix[kl_C, kl_S * kl_NCt] kl_LLq;
      if (kl_re) {
        kl_D = (kre_c[:, kl_ctry] + kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA)[:, kl_adm1])\';
        kl_E = kl_disc_re_term(kl_KBq, kl_D);
      }
      kl_LLq = kl_disc_country_ll(nlp_lin, kl_KBq, kbeta, kl_Yv, kl_lgY1c, kl_cs, kl_ce, kl_E, kl_re);
      matrix[kl_C, kl_K] kl_lsm;
      vector[kl_S * kl_M] kl_w;
      int kl_d;
      int kl_s = 1;
      int kl_m;
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      kl_w = softmax(kl_qlp + kl_disc_combo_ll(kl_LLq, kl_lsm, kl_qidx));
      kl_d = categorical_rng(kl_w);
      while (kl_d > kl_s * kl_M) kl_s += 1;
      kl_m = kl_d - (kl_s - 1) * kl_M;
      kwidth[1] = kl_width[kl_s];
      kl_lag_prob = rep_matrix(0, kl_K, kl_NCt);
      for (s in 1:kl_S) for (m in 1:kl_M) for (k in 1:kl_K)
        kl_lag_prob[k, kl_combo[m, k]] += kl_w[(s - 1) * kl_M + m];
      for (k in 1:kl_K) {
        klag[k] = kl_cent[kl_combo[kl_m, k]];
        kl_llc[:, k] = col(kl_LLq[k], kl_qidx[kl_s, k, kl_m]);
        if (kl_gq_row_loglik) {
          vector[N] eta = nlp_lin;
          for (j in 1:kl_J) eta += kbeta[j, k] * col(kl_KBq[j], kl_qidx[kl_s, k, kl_m]);
          if (kl_re) eta += col(kl_E, kl_qidx[kl_s, k, kl_m]);
          kl_ll[:, k] = kl_Yv .* eta - exp(eta) - kl_lgY1;
        }
      }
    } else if (kl_rs) {
      // country path (row log-lik not available: kernel.gq_row_loglik refused)
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      matrix[kl_re ? 0 : kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_T, kl_J] kl_BW;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        for (k in 1:kl_K) kl_BW[k] = kl_KBu * kTw[k];
      } else {
        kl_U = kl_KBu * ktheta;
      }
      if (kl_pairs) {
        vector[kl_P] kl_esum = rep_vector(0, kl_P);
        vector[N] kl_el = exp(nlp_lin);
        for (n in 1:N) kl_esum[kl_pix[n]] += kl_el[n];
        for (c in 1:kl_C)
          kl_llc[c] = (dot_product(kl_Yv[kl_cs[c]:kl_ce[c]], nlp_lin[kl_cs[c]:kl_ce[c]]) - kl_lgY1c[c]
                       + kl_pair_logp(kl_esum[kl_pcs[c]:kl_pce[c]], kl_Ypair[kl_pcs[c]:kl_pce[c]],
                                      kl_pt[kl_pcs[c]:kl_pce[c]], kl_pa[kl_pcs[c]:kl_pce[c]], c,
                                      kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a))\';
      } else {
        for (c in 1:kl_C) {
          int n = kl_ce[c] - kl_cs[c] + 1;
          matrix[kl_re ? n : 0, kl_re ? kl_J : 0] D;
          if (kl_re) D = (rep_matrix(col(kre_c, c), n) + kl_kre_a[:, kl_adm1[kl_cs[c]:kl_ce[c]]])\';
          kl_llc[c] = kl_ctry_logp(nlp_lin[kl_cs[c]:kl_ce[c]], kl_tix[kl_cs[c]:kl_ce[c]],
                                 kl_Yv[kl_cs[c]:kl_ce[c]], kl_lgY1c[c], rep_row_vector(0, kl_K),
                                 kl_U, kl_BW, kbeta, kl_re, D)\';
        }
      }
    } else {
      kl_ll  = kl_class_rowll(nlp_lin, kl_KBu * ktheta, kl_tix, kl_Yv, kl_lgY1);
      kl_llc = kl_country_ll(kl_ll, kl_cs, kl_ce);
    }
    for (c in 1:kl_C) {
      vector[kl_K] lp = kl_member_logp(klogit[c], kl_fix[c])\' + kl_llc[c]\';
      kl_log_lik_ctry[c] = log_sum_exp(lp);
      kl_log_resp[:, c] = lp - kl_log_lik_ctry[c];
    }
    if (kl_gq_row_loglik)
      for (n in 1:N) log_lik[n] = log_sum_exp(kl_log_resp[:, kl_ctry[n]] + kl_ll[n]\');
  }'

# Patcher ---------------------------------------------------------------------

#' The fixed strings the patcher keys on, for a kernel with J basis functions.
#' Each must occur exactly once in brms::stancode() of the stage-1 model.
#' (Statement landmarks for the nlpar priors are matched by regex in
#' patch_kl_latent_class_stancode().)
kl_stan_landmarks <- function(J) {
  nlp <- kl_nlpars(J)
  c(
    header = "// generated with brms ",
    functions = "functions {",
    prior_only = "  int prior_only;  // should the likelihood be ignored?",
    tdata = "transformed data {",
    params = "parameters {",
    tparams = "transformed parameters {",
    model = "model {",
    gq = "generated quantities {",
    if_prior_only = "  if (!prior_only) {",
    lin_decl = "    vector[N] nlp_lin = rep_vector(0.0, N);",
    mu_decl = "    vector[N] mu;",
    nl_loop = "      // compute non-linear predictor values",
    likelihood = "    target += poisson_log_lpmf(Y | mu);",
    lprior_total = "  target += lprior;",
    stats::setNames(
      sprintf("    vector[N] nlp_%s = rep_vector(0.0, N);", nlp),
      paste0("nlp_decl_", nlp)
    ),
    stats::setNames(
      sprintf("    nlp_%s += X_%s * b_%s;", nlp, nlp, nlp),
      paste0("nlp_add_", nlp)
    )
  )
}

# Index of the single line matching landmark `s`: a regex, or else an exact
# line, or else (for prefixes such as the version header) a line prefix.
.kl_find1 <- function(L, s, what = s, regex = FALSE) {
  idx <- if (regex) grep(s, L, perl = TRUE) else which(L == s)
  if (!regex && !length(idx)) {
    idx <- which(startsWith(L, s))
  }
  if (length(idx) != 1L) {
    stop(
      "patch_kl_latent_class_stancode: landmark '",
      what,
      "' found ",
      length(idx),
      " time(s) (expected 1) -- the brms-generated program ",
      "changed layout (brms ",
      as.character(utils::packageVersion("brms")),
      "); update kl_stan_landmarks().",
      call. = FALSE
    )
  }
  idx
}

#' patch_kl_latent_class_stancode
#'
#' Turn the stage-1 brms program (one kernel) into the K-class country
#' mixture. See the section header for the model; `stancode` is
#' brms::stancode() of assemble_kl_model()'s formula/data/prior, generated
#' WITHOUT threading.
#'
#' @param stancode character(1)
#' @param meta     cb_meta entry of the kernel variable (J = meta$df, lags)
#' @param repel    TRUE = also the identify repel data, parameters and prior
#'   (.KL_STAN_REPEL_*); FALSE leaves the shared program unchanged
#' @param lag_re   TRUE = also the country lag offsets (kernel.lag_re,
#'   .KL_STAN_LAG_RE_*, .kl_lag_re_likelihood()); FALSE leaves the shared
#'   program unchanged
#' @param lag_re_scale kernel.lag_re.scale, "month" (the original variant) or
#'   "logit" (.KL_STAN_LAG_RE_*_LOGIT: centres inside the lag window)
#' @return character(1) patched program; attribute "brms_version"
patch_kl_latent_class_stancode <- function(
  stancode,
  meta,
  repel = FALSE,
  lag_re = FALSE,
  lag_re_scale = "month"
) {
  J <- meta$df
  NL <- length(meta$lags)
  L <- strsplit(stancode, "\n", fixed = TRUE)[[1]]
  lm <- kl_stan_landmarks(J)
  nlp <- kl_nlpars(J)

  # Everything that must exist exactly once, checked up front.
  for (nm in names(lm)) {
    .kl_find1(L, lm[[nm]], nm)
  }
  if (any(grepl("reduce_sum|partial_log_lik", L))) {
    stop(
      "patch_kl_latent_class_stancode: threaded (reduce_sum) brms code is ",
      "not supported -- generate the program without threads.",
      call. = FALSE
    )
  }
  brms_version <- sub(
    "^// generated with brms ",
    "",
    L[.kl_find1(L, lm[["header"]])]
  )
  n_C <- sum(grepl("^  vector\\[N\\] C_[0-9]+;$", L))
  if (n_C != NL * J) {
    stop(
      "patch_kl_latent_class_stancode: expected ",
      NL * J,
      " non-linear ",
      "covariates C_1..C_",
      NL * J,
      " (J = ",
      J,
      " x ",
      NL,
      " lags), found ",
      n_C,
      " -- the kernel term must be the only non-linear covariate user.",
      call. = FALSE
    )
  }

  drop <- integer(0)
  # (1) kernel nlpar parameter declarations
  for (p in nlp) {
    drop <- c(
      drop,
      .kl_find1(
        L,
        sprintf("^  vector(<[^>]*>)?\\[K_%s\\] b_%s;", p, p),
        paste0("param b_", p),
        regex = TRUE
      )
    )
  }
  # (2) their prior statements (possibly multi-line, up to the first ';')
  for (p in nlp) {
    i <- .kl_find1(
      L,
      sprintf("^  lprior \\+= [a-z_]+\\(b_%s \\|", p),
      paste0("lprior b_", p),
      regex = TRUE
    )
    j <- i
    while (!grepl(";\\s*$", L[j])) {
      j <- j + 1L
    }
    drop <- c(drop, i:j)
  }
  # (3) model block: nlpar predictors, `mu`, the non-linear loop, likelihood
  for (p in nlp) {
    drop <- c(
      drop,
      .kl_find1(L, lm[[paste0("nlp_decl_", p)]]),
      .kl_find1(L, lm[[paste0("nlp_add_", p)]])
    )
  }
  drop <- c(drop, .kl_find1(L, lm[["mu_decl"]]))
  i <- .kl_find1(L, lm[["nl_loop"]])
  if (
    !identical(trimws(L[i - 1L]), "for (n in 1:N) {") ||
      !startsWith(trimws(L[i + 1L]), "mu[n] = ") ||
      !identical(trimws(L[i + 2L]), "}")
  ) {
    stop(
      "patch_kl_latent_class_stancode: unexpected layout around the ",
      "non-linear predictor loop.",
      call. = FALSE
    )
  }
  drop <- c(drop, (i - 1L):(i + 2L))

  # brms linear-predictor statements (copied into generated quantities):
  # model-block lines after the nlp_lin declaration up to the likelihood,
  # excluding dropped lines, declarations and comment-only lines.
  i_lin <- .kl_find1(L, lm[["lin_decl"]])
  i_lik <- .kl_find1(L, lm[["likelihood"]])
  lin_idx <- setdiff(seq.int(i_lin + 1L, i_lik - 1L), drop)
  lin_stmt <- L[lin_idx]
  lin_stmt <- lin_stmt[
    !grepl("^\\s*//", lin_stmt) &
      nzchar(trimws(lin_stmt)) &
      !grepl("^\\s*vector\\[N\\] ", lin_stmt)
  ]

  kb_assign <- sprintf(
    "  kl_KB[:, %d] = C_%d;",
    seq_len(NL * J),
    seq_len(NL * J)
  )
  tdata_stmt <- c(
    "  // ---- brms_kl latent classes ----",
    "  for (n in 1:N) kl_lgY1[n] = lgamma(Y[n] + 1);",
    kb_assign,
    sub("^\n", "", .KL_STAN_TDATA_STMT)
  )

  # Insertions: (after-line index -> text). Block ends are the line before the
  # next block header.
  i_fun <- .kl_find1(L, lm[["functions"]])
  i_po <- .kl_find1(L, lm[["prior_only"]])
  i_td <- .kl_find1(L, lm[["tdata"]])
  i_pa <- .kl_find1(L, lm[["params"]])
  i_tp <- .kl_find1(L, lm[["tparams"]])
  i_mo <- .kl_find1(L, lm[["model"]])
  i_gq <- .kl_find1(L, lm[["gq"]])
  i_lp <- .kl_find1(L, lm[["lprior_total"]])
  lik <- .KL_STAN_LIKELIHOOD
  gq_stm <- sprintf(.KL_STAN_GQ_STMT, paste0("  ", lin_stmt, collapse = "\n"))
  if (isTRUE(lag_re)) {
    lik <- .kl_lag_re_likelihood(lik, "model", lag_re_scale)
    gq_stm <- .kl_lag_re_likelihood(gq_stm, "gq", lag_re_scale)
  }
  ins <- list(
    list(after = i_fun, text = .KL_STAN_FUNCTIONS),
    list(after = i_po, text = .KL_STAN_DATA),
    list(after = i_td, text = .KL_STAN_TDATA_DECL),
    list(after = i_pa - 2L, text = paste(tdata_stmt, collapse = "\n")), # before TD's "}"
    list(after = i_pa, text = .KL_STAN_PARAMS),
    list(after = i_tp, text = .KL_STAN_TPARAMS_DECL),
    list(after = i_mo - 2L, text = .KL_STAN_TPARAMS_STMT), # before TP's "}"
    list(after = i_lik, text = lik), # replaces likelihood
    list(after = i_lp, text = .KL_STAN_ICAR),
    list(after = i_gq, text = .KL_STAN_GQ_DECL),
    list(after = length(L) - 1L, text = gq_stm)
  )
  # identify = repel: its data, parameters and prior after the shared ones
  if (isTRUE(repel)) {
    ins <- c(
      ins,
      list(
        list(after = i_po, text = .KL_STAN_REPEL_DATA),
        list(after = i_pa, text = .KL_STAN_REPEL_PARAMS),
        list(after = i_mo - 2L, text = .KL_STAN_REPEL_TP)
      )
    )
  }
  # kernel.lag_re: its functions, data, transformed data, parameters and prior
  # after the shared ones (transformed data statements after kl_KBu is filled)
  if (isTRUE(lag_re)) {
    fr <- .kl_lag_re_fragments(lag_re_scale)
    ins <- c(
      ins,
      list(
        list(after = i_fun, text = fr$functions),
        list(after = i_po, text = fr$data),
        list(after = i_td, text = .KL_STAN_LAG_RE_TDATA_DECL),
        list(after = i_pa - 2L, text = .KL_STAN_LAG_RE_TDATA_STMT),
        list(after = i_pa, text = fr$params),
        list(after = i_tp, text = fr$tparams_decl),
        list(after = i_mo - 2L, text = fr$tp),
        list(after = i_lp, text = .KL_STAN_LAG_RE_MODEL)
      )
    )
  }
  # The block-closing "}" lines the "- 2L"/"- 1L" offsets rely on:
  for (k in c(i_pa - 1L, i_mo - 1L, length(L))) {
    if (!identical(trimws(L[k]), "}")) {
      stop(
        "patch_kl_latent_class_stancode: expected a closing brace at line ",
        k,
        " of the brms program.",
        call. = FALSE
      )
    }
  }
  drop <- c(drop, i_lik)

  extra <- vector("list", length(L))
  for (x in ins) {
    extra[[x$after]] <- c(extra[[x$after]], x$text)
  }
  out <- character(0)
  for (k in seq_along(L)) {
    if (!k %in% drop) {
      out <- c(out, L[k])
    }
    if (length(extra[[k]])) out <- c(out, extra[[k]])
  }
  out <- c(
    sprintf(
      "// patched by utils/dlnm_kernel.R::patch_kl_latent_class_stancode (brms_kl latent classes; J = %d, %d lags%s)",
      J,
      NL,
      paste0(
        if (isTRUE(repel)) "; identify repel" else "",
        if (isTRUE(lag_re)) "; lag_re" else "",
        if (isTRUE(lag_re) && identical(lag_re_scale, "logit")) " logit" else ""
      )
    ),
    out
  )
  res <- paste(out, collapse = "\n")
  attr(res, "brms_version") <- brms_version
  res
}

# Data ------------------------------------------------------------------------

#' kl_prior_partition
#'
#' Prior class membership from a groupings registry (kernel.prior_grouping):
#' class k <- the k-th registry group present among the fitted countries
#' (registry order); K = max(kernel.classes, number of such groups), extra
#' classes unanchored. Each anchored class's kmu prior is centred on the
#' group's deployed occurrence lag when `deploy_config` is set (deployed fits
#' found with paths.R::read_deployed_lag()), else on kernel.prior_mu, else
#' the kmu prior stays uniform. With kernel.prior_grouping.contiguity.enable
#' the one-hot rows of A are softened towards the groups of each country's
#' neighbours (kl_contiguity_mix()).
#'
#' @param kcfg resolve_kernel_cfg() output
#' @param country_levels fitted countries (factor-level order)
#' Under identify repel K = kernel.classes and A is C x G (columns = the
#' registry groups present, also returned as `groups`); class labels klc<k>.
#'
#' @return list(use, K, A (C x K; C x G under repel), class_labels, mu_normal,
#'   mu_mean (K), lags_found[, groups])
kl_prior_partition <- function(
  kcfg,
  country_levels,
  groupings_dir = "analysis/groupings"
) {
  K0 <- max(1L, kcfg$classes)
  mid <- (kcfg$lag_min + kcfg$lag_max) / 2
  src <- kcfg$prior_grouping$source
  mu_normal <- !is.null(kcfg[["prior_mu"]])
  base_mean <- if (mu_normal) as.numeric(kcfg[["prior_mu"]]) else mid
  if (is.null(src) || !nzchar(src)) {
    return(list(
      use = FALSE,
      K = K0,
      A = matrix(0, length(country_levels), K0),
      class_labels = sprintf("klc%d", seq_len(K0)),
      mu_normal = mu_normal,
      mu_mean = rep(base_mean, K0),
      lags_found = NULL
    ))
  }
  reg_all <- load_grouping(src, groupings_dir)
  reg <- reg_all[reg_all$iso3 %in% country_levels, , drop = FALSE]
  groups <- unique(reg$group)
  if (!length(groups)) {
    stop("kl_prior_partition: no fitted country is in registry '", src, "'.")
  }
  # identify repel: the registry does not label the classes. K = kernel.classes
  # (fewer, as many or more classes than groups), A = C x G registry membership
  # (columns = groups), classes klc<k> until named after sampling.
  repel <- identical(kcfg$identify, "repel")
  if (repel && kcfg$classes < 2L) {
    stop(
      "brms_kl: kernel.identify repel needs kernel.classes >= 2 (the registry does not set K)."
    )
  }
  K <- if (repel) kcfg$classes else max(kcfg$classes, length(groups))
  labels <- if (repel) {
    sprintf("klc%d", seq_len(K))
  } else {
    c(
      groups,
      if (K > length(groups)) sprintf("klc%d", seq.int(length(groups) + 1L, K))
    )
  }
  cols <- if (repel) groups else labels
  A <- matrix(
    0,
    length(country_levels),
    length(cols),
    dimnames = list(country_levels, cols)
  )
  A[cbind(match(reg$iso3, country_levels), match(reg$group, cols))] <- 1
  cg <- kcfg$prior_grouping[["contiguity"]] %||% list()
  if (isTRUE(cg[["enable"]])) {
    id <- kl_identify(kcfg)
    if (
      !(identical(id, "anchor") ||
        (id %in%
          c("seed", "repel") &&
          isTRUE(kcfg$prior_grouping[["membership"]])))
    ) {
      stop(
        "brms_kl: kernel.prior_grouping.contiguity acts on the kappa * A membership pull, which is ",
        "only on under identify anchor, or identify seed / repel with prior_grouping.membership: yes."
      )
    }
    w <- as.numeric(cg[["weight"]] %||% 0.5)
    A0 <- A
    A <- kl_contiguity_mix(
      A,
      reg_all,
      load_adjacency(
        kcfg$spatial$adjacency %||% "country_adjacency",
        groupings_dir
      ),
      w
    )
    message(sprintf(
      "[brms_kl] contiguity prior (mix, weight %s): %d of %d anchored rows softened",
      format(w),
      sum(rowSums(abs(A - A0)) > 1e-12),
      sum(rowSums(A0) > 0)
    ))
  }

  mu_mean <- rep(base_mean, K)
  lags_found <- NULL
  if (repel) {
    return(list(
      use = TRUE,
      K = K,
      A = A,
      groups = groups,
      class_labels = labels,
      mu_normal = mu_normal,
      mu_mean = pmin(pmax(mu_mean, kcfg$lag_min), kcfg$lag_max),
      lags_found = NULL
    ))
  }
  dc <- kcfg$prior_grouping$deploy_config
  if (!is.null(dc) && nzchar(dc)) {
    dep <- yaml::read_yaml(dc)$deploy
    lags_found <- vapply(
      groups,
      function(g) {
        r <- read_deployed_lag(
          dep$fits$dest,
          g,
          src,
          "occ",
          episode_key = dep$episodes$key %||% "lp_season",
          national_zero_fill = isTRUE(dep$episodes$national_zero_fill),
          national_zero_threshold = .episodes_nz_threshold(dep$episodes),
          occ_panel = .occ_panel_spec_from_episodes(dep$episodes)
        )
        if (is.null(r)) NA_real_ else as.numeric(r$lag)
      },
      numeric(1)
    )
    if (any(is.na(lags_found))) {
      message(
        "[brms_kl] no deployed occurrence fit for group(s) ",
        paste(groups[is.na(lags_found)], collapse = ", "),
        " in ",
        dep$fits$dest,
        " -- their kmu prior stays centred at ",
        base_mean
      )
    }
    hit <- !is.na(lags_found)
    mu_mean[seq_along(groups)][hit] <- lags_found[hit]
    mu_normal <- mu_normal || any(hit)
  }
  list(
    use = TRUE,
    K = K,
    A = A,
    class_labels = labels,
    mu_normal = mu_normal,
    mu_mean = pmin(pmax(mu_mean, kcfg$lag_min), kcfg$lag_max),
    lags_found = lags_found
  )
}

#' kl_contiguity_mix
#'
#' Contiguity-softened prior membership (kernel.prior_grouping.contiguity):
#' each anchored row of the one-hot prior partition becomes
#'   A[c, ] <- (1 - weight) * A[c, ] + weight * s_c,
#' s_c = the share of c's neighbours (adjacency registry; every registry
#' country counts, fitted or not) in each anchored group. Neighbours outside
#' the registry, or in a registry group that is not a class (e.g. a group no
#' fitted country belongs to), are ignored; a country with no such neighbour,
#' or outside the registry (zero row), keeps its row. Interior countries (all
#' neighbours in their own group) are unchanged, and weight 0 returns A
#' exactly. Rows stay on the simplex, so kappa keeps its meaning: the logit
#' gap between a country's own group and a non-contiguous one. Seed rows are
#' irrelevant (their class is fixed by kl_fix).
#'
#' @param A one-hot C x K prior partition (rownames = countries, colnames = class labels)
#' @param reg full registry, load_grouping() (not restricted to fitted countries)
#' @param adj load_adjacency() output (iso3_a, iso3_b)
#' @param weight neighbour share in [0, 1]
#' @return C x K matrix, same dimnames
kl_contiguity_mix <- function(A, reg, adj, weight = 0.5) {
  if (
    !is.numeric(weight) ||
      length(weight) != 1L ||
      is.na(weight) ||
      weight < 0 ||
      weight > 1
  ) {
    stop(
      "kl_contiguity_mix: kernel.prior_grouping.contiguity.weight must be one number in [0, 1] (got ",
      paste(weight, collapse = ", "),
      ")."
    )
  }
  grp <- stats::setNames(reg$group, reg$iso3)
  for (c in rownames(A)[rowSums(A) > 0]) {
    nb <- c(adj$iso3_b[adj$iso3_a == c], adj$iso3_a[adj$iso3_b == c])
    g <- grp[intersect(nb, names(grp))]
    g <- g[g %in% colnames(A)]
    if (!length(g)) {
      next
    }
    s <- as.numeric(table(factor(g, levels = colnames(A)))) / length(g)
    A[c, ] <- (1 - weight) * A[c, ] + weight * s
  }
  A
}

#' kl_adjacency_data
#'
#' Stan inputs of the spatial membership term: the registry edges restricted
#' to the fitted countries, as index pairs, connected components (BFS) and
#' their sizes (ICAR), and the DAGAR parent lists under the country order.
#' `adj` = load_adjacency() output (or NULL: no edges).
#' @return list(n_edges, node1, node2, n_comp, comp, comp_n, isolated,
#'   n_par, par_nodes, par_start)
kl_adjacency_data <- function(adj, country_levels) {
  C <- length(country_levels)
  if (is.null(adj)) {
    adj <- data.frame(iso3_a = character(0), iso3_b = character(0))
  }
  keep <- adj$iso3_a %in% country_levels & adj$iso3_b %in% country_levels
  n1 <- match(adj$iso3_a[keep], country_levels)
  n2 <- match(adj$iso3_b[keep], country_levels)
  nb <- lapply(seq_len(C), function(i) c(n2[n1 == i], n1[n2 == i]))
  comp <- integer(C)
  g <- 0L
  for (s in seq_len(C)) {
    if (comp[s] > 0L) {
      next
    }
    g <- g + 1L
    queue <- s
    comp[s] <- g
    while (length(queue)) {
      v <- queue[1L]
      queue <- queue[-1L]
      for (w in nb[[v]]) {
        if (comp[w] == 0L) {
          comp[w] <- g
          queue <- c(queue, w)
        }
      }
    }
  }
  # DAGAR parent structure under the natural (country_levels) ordering:
  # parents of i = neighbours j < i (port of archive/cholera_regime_sim
  # R/00_utils.R::dagar_parents()).
  parents <- lapply(seq_len(C), function(i) sort(unique(nb[[i]][nb[[i]] < i])))
  n_par <- vapply(parents, length, integer(1))
  par_start <- c(1L, cumsum(n_par) + 1L)[seq_len(C)]
  list(
    n_edges = length(n1),
    node1 = as.array(as.integer(n1)),
    node2 = as.array(as.integer(n2)),
    n_comp = g,
    comp = as.array(comp),
    comp_n = as.array(as.numeric(tabulate(comp, g))),
    isolated = country_levels[lengths(nb) == 0L],
    n_par = as.array(as.integer(n_par)),
    par_nodes = as.array(as.integer(unlist(parents))),
    par_start = as.array(as.integer(par_start))
  )
}

#' kl_lag_basis
#'
#' Cubic B-spline basis over the lag grid for kernel.form dlnm (intercept
#' included, so rows sum to 1 and a constant lag-response is representable).
#' Local support: each function spans a few lags, so a zero coefficient group
#' means no response over that window -- what the group horseshoe needs.
#' @return NL x lag_df matrix
kl_lag_basis <- function(lags, lag_df = 7L) {
  lag_df <- max(4L, as.integer(lag_df)) # knots ~ (NL - 1) / (lag_df - 3) lags apart
  B <- splines::bs(as.numeric(lags), df = lag_df, degree = 3, intercept = TRUE)
  matrix(
    as.numeric(B),
    nrow = length(lags),
    dimnames = list(as.character(lags), NULL)
  )
}

#' kl_erf_direction
#'
#' The exposure-response contrast used by kernel.identify = "erf": v =
#' B(erf_at) - B(cen) (so v' kbeta[, k] is class k's log RR of a sustained
#' ENSO value erf_at vs cen) and an orthonormal basis Q of its complement, so
#' kbeta[, k] = kc[k] v / |v|^2 + Q eta[, k] is a bijection.
#' @return list(v = numeric(J), Q = J x (J - 1) matrix)
kl_erf_direction <- function(meta, kcfg) {
  v <- as.numeric(.kl_basis_centred(
    meta,
    kcfg$erf_at %||% 1.5,
    kcfg$cen %||% 0
  ))
  if (sum(v^2) < 1e-12) {
    stop(
      "brms_kl: kernel.erf_at gives a zero ERF contrast; choose a value away from cen."
    )
  }
  Q <- qr.Q(qr(matrix(v, ncol = 1)), complete = TRUE)[, -1, drop = FALSE]
  list(v = v, Q = Q)
}

#' kl_disc_combos
#'
#' Every non-decreasing K-tuple of centre indices 1..NC (ties allowed), in
#' lexicographic order: the ordered representation of K exchangeable class
#' centres. M = choose(NC + K - 1, K) rows (strictly increasing K-subsets of
#' 1..NC+K-1 shifted by 0..K-1).
#' @return integer matrix M x K
kl_disc_combos <- function(NC, K) {
  NC <- as.integer(NC)
  K <- as.integer(K)
  if (K == 1L) {
    return(matrix(seq_len(NC), ncol = 1L))
  }
  cb <- t(utils::combn(NC + K - 1L, K)) -
    matrix(seq_len(K) - 1L, choose(NC + K - 1L, K), K, byrow = TRUE)
  storage.mode(cb) <- "integer"
  cb
}

#' kl_disc_spec
#'
#' Kernel grid of kernel.form discrete (kernel.discrete): candidate centres
#' (every centre_step-th lag of the grid), widths, the NL x (S * NC) lag
#' weights (column (s - 1) * NC + i = kernel of width s centred on centre i;
#' width 0 = all weight on that lag), the ordered centre combinations and
#' their log prior: log(K! / prod_i n_i!) under combo_prior "exchangeable"
#' (the prior iid uniform centres of exchangeable classes induce on the
#' ordered representation), 0 under "uniform"; plus, when `mu_normal`, a
#' discretised normal(mu_mean[k], mu_sd) on the centre of position k
#' (`mu_mean` must be ascending, like the positions).
#' Refuses grids above `max_terms` (S * M) -- the per-gradient cost is
#' ~ S * M * C * K.
#' @return list(centres, widths, Wq, combos, combo_lp, M, centre_step, combo_prior)
kl_disc_spec <- function(
  kcfg,
  lags,
  K,
  mu_normal = FALSE,
  mu_mean = NULL,
  mu_sd = 3,
  max_terms = 2e5
) {
  dc <- kcfg$discrete %||% list()
  widths <- as.numeric(unlist(dc$widths %||% 0))
  step <- as.integer(dc$centre_step %||% 1L)
  cprior <- dc$combo_prior %||% "exchangeable"
  if (!length(widths) || any(!is.finite(widths)) || any(widths < 0)) {
    stop("brms_kl: kernel.discrete.widths must be non-negative numbers.")
  }
  if (anyDuplicated(widths)) {
    stop("brms_kl: kernel.discrete.widths has duplicates.")
  }
  if (is.na(step) || step < 1L) {
    stop("brms_kl: kernel.discrete.centre_step must be >= 1.")
  }
  if (!cprior %in% c("exchangeable", "uniform")) {
    stop("brms_kl: kernel.discrete.combo_prior must be exchangeable | uniform.")
  }
  lags <- as.numeric(lags)
  cent <- lags[seq(1L, length(lags), by = step)]
  NC <- length(cent)
  M <- choose(NC + K - 1, K)
  if (length(widths) * M > max_terms) {
    stop(sprintf(
      paste0(
        "brms_kl: kernel.form discrete with K = %d, %d centres and %d width(s) ",
        "sums %s terms per gradient (> %s): raise kernel.discrete.centre_step, ",
        "use fewer widths or fewer classes."
      ),
      K,
      NC,
      length(widths),
      format(length(widths) * M, big.mark = ","),
      format(max_terms, big.mark = ",")
    ))
  }
  Wq <- do.call(
    cbind,
    lapply(widths, function(w) t(kl_weights(cent, rep(w, NC), lags)))
  )
  dimnames(Wq) <- list(
    as.character(lags),
    sprintf("w%g_c%g", rep(widths, each = NC), rep(cent, length(widths)))
  )
  combos <- kl_disc_combos(NC, K)
  lp <- if (identical(cprior, "exchangeable")) {
    lfactorial(K) -
      apply(combos, 1, function(r) sum(lfactorial(tabulate(r, NC))))
  } else {
    numeric(nrow(combos))
  }
  if (isTRUE(mu_normal)) {
    mu_mean <- rep_len(as.numeric(mu_mean), K)
    for (k in seq_len(K)) {
      pk <- stats::dnorm(cent, mu_mean[k], mu_sd, log = TRUE)
      pk <- pk - (max(pk) + log(sum(exp(pk - max(pk)))))
      lp <- lp + pk[combos[, k]]
    }
  }
  list(
    centres = cent,
    widths = widths,
    Wq = Wq,
    combos = combos,
    combo_lp = as.numeric(lp),
    M = nrow(combos),
    centre_step = step,
    combo_prior = cprior
  )
}

#' build_kl_latent_standata
#'
#' brms::standata() of the stage-1 model plus the kl_* entries of the patched
#' program. `m$data` must already be sorted by country (kl_sort_by_country()).
#' @param m   assemble_kl_model() output (sorted)
#' @param pg  kl_prior_partition() output
#' @param adj kl_adjacency_data() output, or NULL (spatial term off)
#' @param fix integer per country (country order): 0 = free, k = fixed to class
#'   k (seed anchoring, kernel.pin); NULL = all free
#' @param threads within-chain threads (brms_kl.threads_per_chain): > 1 routes
#'   forms gaussian / dlnm through the reduce_sum country path (kl_rs)
build_kl_latent_standata <- function(
  m,
  pg,
  adj = NULL,
  fix = NULL,
  threads = 1L
) {
  kcfg <- m$kcfg
  dat <- m$data
  ctry <- droplevels(as.factor(dat$country))
  ci <- as.integer(ctry)
  if (is.unsorted(ci)) {
    stop("build_kl_latent_standata: rows must be sorted by country.")
  }
  r <- rle(ci)
  ce <- cumsum(r$lengths)
  cs <- ce - r$lengths + 1L
  sd <- brms::standata(
    m$formula,
    data = dat,
    family = m$family,
    prior = m$priors
  )
  sd <- unclass(sd)
  if (sd$N != nrow(dat)) {
    stop(
      "build_kl_latent_standata: brms kept ",
      sd$N,
      " of ",
      nrow(dat),
      " rows (NA ",
      "in a model variable?) -- the country index would be misaligned. ",
      "assemble_kl_model() should have dropped them."
    )
  }
  use_sp <- isTRUE(kcfg$spatial$enable) && !is.null(adj)
  if (!use_sp) {
    adj <- list(
      n_edges = 0L,
      node1 = integer(0),
      node2 = integer(0),
      n_comp = 0L,
      comp = rep(0L, nlevels(ctry)),
      comp_n = numeric(0),
      n_par = rep(0L, nlevels(ctry)),
      par_nodes = integer(0),
      par_start = rep(1L, nlevels(ctry))
    )
  }
  sp_prior <- kl_spatial_prior(kcfg)
  ident <- kl_identify(kcfg)
  anchor <- isTRUE(pg$use) && identical(ident, "anchor")
  # kernel.prior_grouping.membership: the kappa * A membership pull also under
  # identify seed (class k = registry group k there too); kmu priors unchanged.
  # Under identify repel the pull is label-free (.KL_STAN_REPEL_TP: kappa *
  # group-level class logits); kl_A is then zero and kl_use_pg only switches
  # kkappa on.
  mem <- isTRUE(kcfg$prior_grouping[["membership"]])
  if (mem && !ident %in% c("seed", "anchor", "repel")) {
    stop(
      "brms_kl: kernel.prior_grouping.membership needs identify seed, anchor or repel (got '",
      ident,
      "')."
    )
  }
  repel <- identical(ident, "repel")
  use_pg <- anchor || (mem && isTRUE(pg$use))
  # Without anchoring, class k is not registry group k: drop the group-centred
  # kmu priors (kernel.prior_mu or uniform instead).
  mu_normal <- if (anchor) pg$mu_normal else !is.null(kcfg[["prior_mu"]])
  mu_mean <- if (anchor) {
    pg$mu_mean
  } else {
    rep(
      as.numeric(kcfg[["prior_mu"]] %||% ((kcfg$lag_min + kcfg$lag_max) / 2)),
      pg$K
    )
  }
  ev <- kl_erf_direction(m$meta, kcfg)
  dlnm <- identical(kcfg$form %||% "gaussian", "dlnm")
  disc <- identical(kcfg$form %||% "gaussian", "discrete")
  # kernel.re_slope: ADM1 units are country:ADM1 (ADM1 labels can repeat across countries)
  gauss <- identical(kcfg$form %||% "gaussian", "gaussian")
  use_re <- (disc || gauss) && isTRUE(kcfg$re_slope$enable)
  # forms gaussian / dlnm: country path when threaded or with random slopes
  use_rs <- !disc && (as.integer(threads %||% 1L) > 1L || use_re)
  if (use_rs && isTRUE(kcfg$gq_row_loglik)) {
    stop(
      "brms_kl: kernel.gq_row_loglik is not available on the threaded / random-slope path."
    )
  }
  ksd_fix <- gauss && !is.null(kcfg[["ksd_fixed"]])
  # Threaded country path: countries above N / (2 * threads) rows are split by
  # rows (one reduce_sum per class) so that the largest country does not bound
  # the speed-up; the others are summed whole by the country reduce_sum.
  split <- kl_country_split(
    r$lengths,
    if (use_rs) as.integer(threads %||% 1L) else 1L,
    enable = !isFALSE(kcfg$split_big)
  )
  adm1 <- interaction(
    as.character(dat$country),
    as.character(dat$ADM1),
    drop = TRUE,
    lex.order = TRUE
  )
  # Discrete centres: positions are lag ranks, so any centre priors (prior_mu,
  # or the prior partition's deployed lags) are sorted onto the positions.
  ds <- if (disc) {
    kl_disc_spec(
      kcfg,
      m$meta$lags,
      pg$K,
      mu_normal = pg$mu_normal,
      mu_mean = sort(pg$mu_mean),
      mu_sd = kcfg$prior_mu_sd
    )
  } else {
    list(
      centres = m$meta$lags[1],
      widths = 1,
      M = 1L,
      Wq = matrix(0, length(m$meta$lags), 1L),
      combos = matrix(1L, 1L, pg$K),
      combo_lp = 0
    )
  }
  Lb <- kl_lag_basis(m$meta$lags, kcfg$lag_df %||% 10L)
  hs <- kcfg$hs %||% list()
  rho_ab <- as.numeric(unlist(kcfg$spatial$rho_prior %||% c(2, 2)))
  kfix <- kcfg$prior_grouping[["kappa"]] # [[ ]]: `$` would partially match kappa_prior_sd
  tix <- .kl_exposure_rows(sd, m$meta$df * length(m$meta$lags))
  prs <- .kl_pairs(ci, as.integer(adm1), tix, by_adm1 = use_re)
  # repulsion data: identify repel, or identify seed with kernel.repel.enable
  # (whose registry pull stays the labelled kappa * A: no label-free pull, G = 0)
  rep_data <- if (kl_repel_on(kcfg)) {
    rc <- kl_repel_cfg(kcfg)
    rg <- if (repel && use_pg) pg$A else matrix(0, nlevels(ctry), 0L)
    list(
      kl_rep_NX = length(rc$erf_grid),
      kl_rep_V = matrix(
        .kl_basis_centred(m$meta, rc$erf_grid, kcfg$cen %||% 0),
        length(rc$erf_grid)
      ),
      kl_rep_lag_scale = rc$lag_scale,
      kl_rep_erf_scale = rc$erf_scale,
      kl_rep_G = ncol(rg),
      kl_rep_A = unname(rg)
    )
  }
  # country lag offsets (kernel.lag_re): only their prior sd is extra data
  lag_re_data <- if (kl_lag_re_on(kcfg)) {
    list(kl_lag_re_sd = kl_lag_re_cfg(kcfg)$sd)
  }
  c(
    sd,
    rep_data,
    lag_re_data,
    list(
      kl_K = pg$K,
      kl_C = nlevels(ctry),
      kl_ctry = ci,
      kl_cs = as.array(cs),
      kl_ce = as.array(ce),
      kl_J = m$meta$df,
      kl_NL = length(m$meta$lags),
      kl_lags = as.array(as.numeric(m$meta$lags)),
      kl_T = max(tix),
      kl_tix = as.array(tix),
      kl_pairs = as.integer(use_rs && !isFALSE(kcfg[["pair_collapse"]])),
      kl_P = length(prs$pt),
      kl_pix = as.array(prs$pix),
      kl_pt = as.array(prs$pt),
      kl_pa = as.array(prs$pa),
      kl_pcs = as.array(prs$pcs),
      kl_pce = as.array(prs$pce),
      kl_maxp = as.integer(max(prs$pce - prs$pcs + 1L)),
      kl_dlnm = as.integer(dlnm),
      kl_BL = ncol(Lb),
      kl_Lb = Lb,
      kl_disc = as.integer(disc),
      kl_S = length(ds$widths),
      kl_NCt = length(ds$centres),
      kl_cent = as.array(as.numeric(ds$centres)),
      kl_width = as.array(as.numeric(ds$widths)),
      kl_Wq = unname(ds$Wq),
      kl_M = as.integer(ds$M),
      kl_combo = unname(ds$combos),
      kl_combo_lp = as.array(as.numeric(ds$combo_lp)),
      kl_re = as.integer(use_re),
      kl_NA = nlevels(adm1),
      kl_adm1 = as.array(as.integer(adm1)),
      kl_re_ctry = as.integer(!isFALSE(kcfg$re_slope[["country"]])),
      kl_re_adm1 = as.integer(!isFALSE(kcfg$re_slope[["adm1"]])),
      kl_re_sdc_fix = as.integer(!is.null(kcfg$re_slope[["sd_country_fixed"]])),
      kl_re_sdc_value = as.numeric(kcfg$re_slope[["sd_country_fixed"]] %||% 0),
      kl_re_sd = as.array(c(
        as.numeric(kcfg$re_slope[["sd_country"]] %||% 0.35),
        as.numeric(kcfg$re_slope[["sd_adm1"]] %||% 0.35)
      )),
      kl_rs = as.integer(use_rs),
      kl_grain = max(1L, as.integer(kcfg$grainsize %||% 1L)),
      kl_maxn = as.integer(max(c(1L, r$lengths[split$small]))),
      kl_nsmall = length(split$small),
      kl_small = as.array(as.integer(split$small)),
      kl_nbig = length(split$big),
      kl_big = as.array(as.integer(split$big)),
      kl_grain_rows = max(1L, as.integer(kcfg$grainsize_rows %||% 1L)),
      kl_ksd_fix = as.integer(ksd_fix),
      kl_ksd_value = as.numeric(kcfg[["ksd_fixed"]] %||% 1),
      kl_hs = as.integer(
        dlnm && identical(kcfg$lag_prior %||% "horseshoe", "horseshoe")
      ),
      kl_hs_scale_global = as.numeric(hs$scale_global %||% 0.1),
      kl_hs_slab_scale = as.numeric(hs$slab_scale %||% 1),
      kl_hs_slab_df = as.numeric(hs$slab_df %||% 4),
      kl_alpha = kcfg$alpha,
      kl_mu_normal = as.integer(mu_normal),
      kl_mu_mean = as.array(as.numeric(mu_mean)),
      kl_mu_sd = kcfg$prior_mu_sd,
      kl_sd_lb = kcfg$sigma_lb,
      kl_sd_sd = kcfg$prior_sigma_sd,
      kl_beta_sd = kcfg$prior_beta_sd,
      kl_use_pg = as.integer(use_pg),
      kl_A = if (repel) matrix(0, nlevels(ctry), pg$K) else pg$A,
      kl_mu_ordered = as.integer(identical(ident, "kmu") && !dlnm),
      kl_erf_ordered = as.integer(identical(ident, "erf")),
      kl_erf_v = as.array(ev$v),
      kl_erf_Q = ev$Q,
      kl_kappa_fixed = as.integer(!is.null(kfix)),
      kl_kappa_value = as.numeric(kfix %||% 0),
      kl_kappa_sd = as.numeric(kcfg$prior_grouping$kappa_prior_sd %||% 2),
      kl_use_sp = as.integer(use_sp),
      # symmetric field (K columns) when classes are exchangeable (relabel, repel; and
      # lagorder, whose tied classes are); otherwise class K is the reference
      kl_field_cols = as.integer(
        if (ident %in% c("relabel", "lagorder", "repel")) pg$K else pg$K - 1L
      ),
      kl_n_edges = adj$n_edges,
      kl_node1 = as.array(adj$node1),
      kl_node2 = as.array(adj$node2),
      kl_n_comp = adj$n_comp,
      kl_comp = as.array(as.integer(adj$comp)),
      kl_comp_n = as.array(adj$comp_n),
      kl_tau_sd = as.numeric(kcfg$spatial$tau_prior_sd %||% 1),
      kl_sp_dagar = as.integer(identical(sp_prior, "dagar")),
      kl_NP = length(adj$par_nodes),
      kl_n_par = as.array(as.integer(adj$n_par)),
      kl_par_nodes = as.array(as.integer(adj$par_nodes)),
      kl_par_start = as.array(as.integer(adj$par_start)),
      kl_rho_a = rho_ab[1],
      kl_rho_b = rho_ab[2],
      kl_gq_row_loglik = as.integer(kcfg$gq_row_loglik),
      kl_fix = as.array(as.integer(fix %||% rep(0L, nlevels(ctry))))
    )
  )
}

#' Distinct rows of the kernel exposure basis (brms C_1..C_<n_col>): the index
#' of each row's distinct row, in order of first appearance. Rows are matched
#' on the exact bits of every value ("%a"), so kl_KBu[kl_tix[n]] == kl_KB[n]
#' exactly whatever the exposure; for ENSO (one series) a distinct row is a
#' calendar month.
#' @noRd
.kl_exposure_rows <- function(sd, n_col) {
  cols <- sprintf("C_%d", seq_len(n_col))
  if (!all(cols %in% names(sd))) {
    stop(
      "build_kl_latent_standata: brms data lacks the kernel covariates C_1..C_",
      n_col,
      "."
    )
  }
  key <- do.call(
    paste,
    c(lapply(sd[cols], function(v) sprintf("%a", as.numeric(v))), sep = "|")
  )
  match(key, unique(key))
}

#' Pairs of the pair path (kernel.pair_collapse): (exposure row, ADM1 unit)
#' with random slopes (the kernel term varies by ADM1), else (exposure row,
#' country). Sorted by country, unit, exposure row, so each country's pairs
#' are contiguous (pcs..pce).
#' @param ci country index of each row (1..C, rows sorted by country)
#' @param a1 ADM1 unit index of each row (country:ADM1)
#' @param tix distinct exposure row of each row (.kl_exposure_rows())
#' @return list(pix (pair of each row), pt, pa (exposure row / ADM1 unit of each pair), pcs, pce)
#' @noRd
.kl_pairs <- function(ci, a1, tix, by_adm1) {
  g <- if (by_adm1) a1 else rep(0L, length(ci))
  key <- paste(ci, g, tix)
  f <- which(!duplicated(key))
  f <- f[order(ci[f], g[f], tix[f])]
  pix <- match(key, key[f])
  r <- rle(ci[f])
  if (!identical(r$values, seq_len(max(ci)))) {
    stop(
      ".kl_pairs: every country must have rows (country index 1..C, sorted)."
    )
  }
  pce <- cumsum(r$lengths)
  list(
    pix = pix,
    pt = tix[f],
    pa = a1[f],
    pcs = pce - r$lengths + 1L,
    pce = pce
  )
}

#' kl_country_split -- which countries the threaded country path sums whole
#' (small) and which it splits by rows (big): big = more than N / (2 threads)
#' rows, only when threads > 1 and `enable`.
#' @param n rows per country (country order)
#' @return list(small, big) of country indices
kl_country_split <- function(n, threads = 1L, enable = TRUE) {
  big <- if (enable && threads > 1L) {
    which(n > sum(n) / (2 * threads))
  } else {
    integer(0)
  }
  list(small = setdiff(seq_along(n), big), big = big)
}

#' Sort the model data by country (stable) -- the mixture sums contiguous rows.
kl_sort_by_country <- function(m) {
  o <- order(as.integer(as.factor(m$data$country)), seq_len(nrow(m$data)))
  m$data <- m$data[o, , drop = FALSE]
  m
}

#' Neutral starts for brms's own `lin` parameters, sized from its data list:
#' standardised group effects z_g = 0 with sd_g = 0.1, penalised-spline
#' coefficients zs = 0 with sds = 0.5. Left to cmdstanr's default (-2, 2)
#' unconstrained draws they made most multi-start optimisations fail on
#' non-finite gradients. Correlated-RE Cholesky factors are left to the default.
#' @noRd
.kl_brms_lin_inits <- function(sd) {
  out <- list()
  for (g in as.integer(sub(
    "^M_",
    "",
    grep("^M_[0-9]+$", names(sd), value = TRUE)
  ))) {
    M <- sd[[sprintf("M_%d", g)]]
    N <- sd[[sprintf("N_%d", g)]]
    out[[sprintf("sd_%d", g)]] <- as.array(rep(0.1, M))
    out[[sprintf("z_%d", g)]] <- matrix(0, M, N)
  }
  for (i in as.integer(sub(
    "^nb_lin_",
    "",
    grep("^nb_lin_[0-9]+$", names(sd), value = TRUE)
  ))) {
    kn <- sd[[sprintf("knots_lin_%d", i)]]
    out[[sprintf("sds_lin_%d", i)]] <- as.array(rep(
      0.5,
      sd[[sprintf("nb_lin_%d", i)]]
    ))
    for (j in seq_along(kn)) {
      out[[sprintf("zs_lin_%d_%d", i, j)]] <- as.array(rep(0, kn[j]))
    }
  }
  out
}

#' make_kl_class_inits
#'
#' Per-chain partial init lists for the class parameters (everything else is
#' cmdstanr's default random init). strategy "spread": class lag centres at
#' evenly spaced positions of the lag range, jittered by chain; "wide": every
#' class centred mid-range with a kernel a quarter of the range wide, so each
#' class sees all lags at the start (a narrow kernel far from the signal gets
#' almost no gradient and can settle on a null effect -- synthetic case E);
#' "deployed": centred on pg$mu_mean (the prior partition's deployed lags)
#' where known; "random": no class inits; "auto" (default): "wide" under
#' identify = erf, else "spread". The `lin` coefficients start at 0 except the
#' intercept, set to the log of the crude rate.
#' "sweep": each class centre starts at one of its group's sweep lag modes
#' (`mu_modes`, a list per class from kl_sweep_modes(); chain i takes mode
#' ((i - 1) mod n) + 1, so chains split across a multimodal group); classes
#' without modes start mid-range.
make_kl_class_inits <- function(
  sd,
  pg,
  kcfg,
  n_chains,
  strategy = kcfg$init,
  seed = 1L,
  mu_modes = NULL
) {
  # "auto": wide kernels when classes are identified by the ERF contrast (their
  # centres are free and must find the lag from anywhere), else "spread".
  if (is.null(strategy) || strategy %in% c("auto", "optimize")) {
    strategy <- if (sd$kl_erf_ordered == 1L) "wide" else "spread"
  }
  if (identical(strategy, "random")) {
    return(NULL)
  }
  K <- pg$K
  J <- sd$kl_J
  lo <- kcfg$lag_min
  hi <- kcfg$lag_max
  rate0 <- log(max(sum(sd$Y), 1) / sum(exp(sd$offsets_lin)))
  lapply(seq_len(n_chains), function(i) {
    set.seed(seed + i)
    base <- lo + (hi - lo) * (seq_len(K) - 0.5) / K
    if (identical(strategy, "deployed") && pg$use && sd$kl_use_pg == 1L) {
      base <- pg$mu_mean
    }
    if (identical(strategy, "wide")) {
      base <- rep((lo + hi) / 2, K)
    }
    jit <- stats::runif(K, -0.75, 0.75)
    if (strategy %in% c("sweep", "given")) {
      base <- vapply(
        seq_len(K),
        function(k) {
          md <- mu_modes[[k]]
          if (is.null(md) || !length(md)) {
            (lo + hi) / 2
          } else {
            md[(i - 1L) %% length(md) + 1L]
          }
        },
        numeric(1)
      )
      jit <- stats::runif(K, -0.25, 0.25)
    }
    mu <- pmin(pmax(base + jit, lo + 0.05), hi - 0.05)
    ksd0 <- if (identical(strategy, "wide")) {
      (hi - lo) / 4
    } else {
      max(1.5, kcfg$sigma_lb + 0.5)
    }
    b_lin <- rep(0, sd$K_lin)
    b_lin[1] <- rate0
    init <- c(
      list(kpi = rep(1 / K, K), b_lin = as.array(b_lin)),
      .kl_brms_lin_inits(sd)
    )
    if (sd$kl_re == 1L) {
      # random slopes start at 0, small sds
      # country block only when on: a 0-row init matrix is written to JSON as a
      # flat [] and rejected by Stan (dims (0, C) vs (0))
      if (sd$kl_re_ctry == 1L) {
        init$kre_zc <- matrix(0, J, sd$kl_C)
        if (sd$kl_re_sdc_fix == 0L) init$kre_sdc <- as.array(rep(0.1, J))
      }
      if (sd$kl_re_adm1 == 1L) {
        init$kre_za <- matrix(0, J, sd$kl_NA)
        init$kre_sda <- as.array(rep(0.1, J))
      }
    }
    if (isTRUE(sd$kl_rep_G > 0L)) {
      init$kgrp <- matrix(0, sd$kl_rep_G, K)
    } # identify repel: no pull at the start
    if (!is.null(sd$kl_lag_re_sd)) {
      # kernel.lag_re: no country offset at the start
      init$klag_z <- as.array(rep(0, sd$kl_C))
      init$klag_sd <- as.array(0.5)
    }
    if (sd$kl_dlnm == 1L) {
      BL <- sd$kl_BL
      init$kz <- matrix(stats::rnorm(J * BL * K, 0, 0.3), J * BL, K)
      if (sd$kl_hs == 1L) {
        init$khs_tau <- as.array(rep(sd$kl_hs_scale_global, K))
        init$khs_lambda <- matrix(1, BL, K)
        init$khs_c2 <- as.array(rep(sd$kl_hs_slab_scale^2, K))
      }
      return(init)
    }
    if (sd$kl_disc == 1L) {
      # centres and width are summed out
      init$kbeta_free <- matrix(stats::rnorm(J * K, 0, 0.05), J, K)
      return(init)
    }
    if (sd$kl_ksd_fix == 0L) {
      init$ksd <- rep(ksd0, K)
    }
    if (sd$kl_mu_ordered == 1L) {
      init$kmu_ord <- as.array(stats::qlogis((sort(mu) - lo) / (hi - lo)))
    } else {
      init$kmu_free <- as.array(mu)
    }
    if (sd$kl_erf_ordered == 1L) {
      init$kc <- as.array(sort(
        seq(-0.5, 0.5, length.out = K) + stats::runif(K, -0.05, 0.05)
      ))
      init$keta <- matrix(stats::rnorm((J - 1) * K, 0, 0.05), J - 1, K)
    } else {
      init$kbeta_free <- matrix(stats::rnorm(J * K, 0, 0.05), J, K)
    }
    init
  })
}

# Fit -------------------------------------------------------------------------

#' kl_sweep_modes
#'
#' Lag modes of each group in an existing single-lag sweep (kernel.init
#' sweep): the lags of its LOO weights exp(elpd - max) (kl_pin_sweep_orders()),
#' taken greedily by weight, keeping those with weight >= min_w and at least
#' `sep` lags from a mode already kept (the top lag always). Groups without a
#' ranking get NULL (mid-range start).
#' @return named list group -> numeric lags (by decreasing weight)
kl_sweep_modes <- function(
  batch,
  groups,
  outdir = "generated_data/",
  min_w = 0.05,
  sep = 3
) {
  out <- lapply(groups, function(g) {
    p <- tryCatch(
      kl_pin_sweep_orders(batch, g, outdir = outdir)[[1]],
      error = function(e) NULL
    )
    if (is.null(p)) {
      return(NULL)
    }
    o <- order(-p$lag_prob)
    keep <- numeric(0)
    for (i in o) {
      if (length(keep) && p$lag_prob[i] < min_w) {
        break
      }
      if (!length(keep) || all(abs(p$lags[i] - keep) >= sep)) {
        keep <- c(keep, p$lags[i])
      }
    }
    keep
  })
  names(out) <- groups
  out
}

#' kl_init_strategy -- resolve kernel.init: "auto" is "optimize" under
#' identify = erf (free class centres: multimodal, see kl_multistart_optimize())
#' and "spread" otherwise.
kl_init_strategy <- function(kcfg, sd) {
  st <- kcfg$init %||% "auto"
  if (
    !st %in%
      c(
        "auto",
        "optimize",
        "spread",
        "wide",
        "deployed",
        "random",
        "sweep",
        "given"
      )
  ) {
    stop(
      "brms_kl: kernel.init must be auto | optimize | spread | wide | deployed | random | sweep | given."
    )
  }
  if (identical(st, "auto")) {
    st <- if (sd$kl_erf_ordered == 1L) "optimize" else "spread"
  }
  st
}

#' kl_multistart_optimize
#'
#' Multi-start posterior-mode search used to initialise the chains
#' (kernel.init = optimize). With free class centres a class whose effect is
#' small can settle with its kernel at the edge of the lag range and a null
#' response, a separate local mode the sampler does not leave (synthetic case
#' E: those modes were 127-142 log units below the correct one, i.e. negligible
#' posterior mass). L-BFGS (jacobian = TRUE, the mode of the sampled density)
#' runs from n starts alternating "wide" and "spread" class inits, in parallel;
#' failed starts are tolerated. Every chain then starts at the best mode.
#'
#' @return list(init = list of per-chain init lists (identical), table =
#'   tibble(start, strategy, lp, kmu[k], kc[k]) over all starts, best = index)
kl_multistart_optimize <- function(
  model,
  sd,
  pg,
  kcfg,
  n_starts,
  n_chains,
  seed = 1L,
  cores = 1L
) {
  pinfo <- model$variables()$parameters
  pars <- names(pinfo)
  one <- function(i) {
    strat <- c("wide", "spread")[(i - 1L) %% 2L + 1L]
    ini <- make_kl_class_inits(
      sd,
      pg,
      kcfg,
      1L,
      strategy = strat,
      seed = seed + 1000L + i
    )
    tryCatch(
      {
        o <- suppressMessages(model$optimize(
          data = sd,
          init = ini,
          seed = seed + i,
          jacobian = TRUE,
          iter = 10000L,
          refresh = 0L
        ))
        d <- posterior::as_draws_df(o$draws())
        list(strategy = strat, lp = o$lp(), draws = d)
      },
      error = function(e) list(strategy = strat, lp = NA_real_, draws = NULL)
    )
  }
  old <- options(cmdstanr_warn_inits = FALSE)
  on.exit(options(old), add = TRUE)
  res <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(n_starts), one, mc.cores = min(cores, n_starts))
  } else {
    lapply(seq_len(n_starts), one)
  }
  K <- sd$kl_K
  tab <- dplyr::bind_rows(lapply(seq_along(res), function(i) {
    r <- res[[i]]
    row <- tibble::tibble(start = i, strategy = r$strategy, lp = r$lp)
    if (!is.null(r$draws)) {
      for (v in c(
        sprintf("kmu[%d]", seq_len(K)),
        sprintf("kc[%d]", seq_len(K))
      )) {
        if (v %in% names(r$draws)) row[[v]] <- r$draws[[v]][1]
      }
    }
    row
  }))
  if (all(is.na(tab$lp))) {
    stop("kl_multistart_optimize: every optimisation start failed.")
  }
  best <- which.max(tab$lp)
  rv <- posterior::as_draws_rvars(res[[best]]$draws)
  init1 <- lapply(intersect(pars, names(rv)), function(v) {
    a <- posterior::draws_of(rv[[v]])
    if (pinfo[[v]]$dimensions == 0L) {
      # declared scalar
      a[1]
    } else {
      array(as.numeric(a), dim = dim(a)[-1])
    } # vector / matrix / array, shape kept
  })
  names(init1) <- intersect(pars, names(rv))
  message(sprintf(
    "[brms_kl] multi-start optimisation: %d/%d starts converged; best lp %.1f (start %d); others %s",
    sum(!is.na(tab$lp)),
    n_starts,
    tab$lp[best],
    best,
    paste(
      sprintf("%.1f", sort(tab$lp[-best], decreasing = TRUE)),
      collapse = ", "
    )
  ))
  list(init = rep(list(init1), n_chains), table = tab, best = best)
}

#' kl_program_file -- patched program of model `m`, written once to
#' `<stan_dir>/occ_kl_latent_<sha1>.stan` (cmdstanr caches the executable next
#' to it). @return list(code, sha, file)
kl_program_file <- function(
  m,
  stan_dir = "generated_data/stan_kl",
  suffix = ""
) {
  code <- patch_kl_latent_class_stancode(
    brms::stancode(
      m$formula,
      data = m$data,
      family = m$family,
      prior = m$priors
    ),
    m$meta,
    repel = kl_repel_on(m$kcfg),
    lag_re = kl_lag_re_on(m$kcfg),
    lag_re_scale = if (kl_lag_re_on(m$kcfg)) {
      kl_lag_re_cfg(m$kcfg)$scale
    } else {
      "month"
    }
  )
  sha <- digest::digest(as.character(code), algo = "sha1", serialize = FALSE)
  dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(
    stan_dir,
    sprintf("occ_kl_latent_%s%s.stan", substr(sha, 1, 12), suffix)
  )
  if (!file.exists(file)) {
    # Written to a temp file and renamed (atomic on one file system): array
    # tasks building the same program never read a half-written .stan.
    tmp <- tempfile(paste0(basename(file), "."), tmpdir = stan_dir)
    writeLines(as.character(code), tmp)
    if (!file.rename(tmp, file)) unlink(tmp)
  }
  list(code = code, sha = sha, file = file)
}

#' kl_with_compile_lock
#'
#' Run `expr` (the cmdstanr compile of `file`) holding a lock on that program.
#' Array tasks of one batch share stan_dir: without the lock several compile
#' the same executable at once and one can exec it while another is still
#' writing it ("Permission denied"). The
#' lock is a directory (mkdir is atomic, no extra package); a waiter gets the
#' finished executable, which cmdstanr then reuses without recompiling. A lock
#' older than `stale_min` minutes (a killed task) is broken.
kl_with_compile_lock <- function(
  file,
  expr,
  wait_min = 90,
  stale_min = 60,
  poll_sec = 5
) {
  lock <- paste0(file, ".lock")
  t0 <- Sys.time()
  repeat {
    if (dir.create(lock, showWarnings = FALSE)) {
      break
    }
    age <- as.numeric(difftime(
      Sys.time(),
      file.info(lock)$mtime,
      units = "mins"
    ))
    if (!is.na(age) && age > stale_min) {
      unlink(lock, recursive = TRUE)
      next
    }
    if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > wait_min) {
      stop(
        "kl_with_compile_lock: waited ",
        wait_min,
        " min for ",
        lock,
        call. = FALSE
      )
    }
    Sys.sleep(poll_sec)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  force(expr)
}

#' kl_compile_model
#'
#' Write and compile the patched program of model `m` under the compile lock:
#' the one compile path, shared by fit_kl_stan() and the batch precompile
#' (04_fit_occ_brms_kl.R --compile_only, utils/kl_compile_batch.R). Within-chain
#' threads (brms_kl.threads_per_chain): reduce_sum over countries on the
#' gaussian / dlnm path; the executable is built with STAN_THREADS and kept
#' apart from the unthreaded one (file suffix _thr).
#' @return kl_program_file() output plus `model` (CmdStanModel) and `threads`.
kl_compile_model <- function(m, brms_cfg, stan_dir = "generated_data/stan_kl") {
  kcfg <- m$kcfg
  threads <- max(1L, as.integer(brms_cfg$threads_per_chain %||% 1L))
  if (identical(kcfg$form %||% "gaussian", "discrete")) {
    threads <- 1L
  }
  prog <- kl_program_file(
    m,
    stan_dir,
    suffix = if (threads > 1L) "_thr" else ""
  )
  model <- kl_with_compile_lock(
    prog$file,
    if (threads > 1L) {
      cmdstanr::cmdstan_model(
        prog$file,
        cpp_options = list(stan_threads = TRUE)
      )
    } else {
      cmdstanr::cmdstan_model(prog$file)
    }
  )
  c(prog, list(model = model, threads = threads))
}

#' fit_kl_stan
#'
#' Generate (brms::stancode), patch and compile the kernel program (any K;
#' K = 1 is the single kernel), then sample with cmdstanr. The patched program is written to
#' `<stan_dir>/occ_kl_latent_<sha1>.stan` (cmdstanr caches the executable
#' next to it, so an identical program is compiled once).
#'
#' @param m        assemble_kl_model() output (sorted: kl_sort_by_country())
#' @param brms_cfg backend block (MCMC settings: warmup/sampling/chains/cores,
#'                 adapt_delta, max_treedepth, refresh)
#' @param pg,adj   kl_prior_partition() / kl_adjacency_data() (adj NULL = no ICAR)
#' @return object of class "kl_fit": list(draws (posterior draws_df), program,
#'   program_sha1, program_file, country_levels, class_labels, K, chains,
#'   diagnostics, inits, data_dims, prior_partition, adjacency)
#' @param fix,labels seed anchoring (kernel.pin): per-country fixed class (see
#'   build_kl_latent_standata()) and the class labels it implies; with seeds the
#'   classes are not exchangeable, so no post-hoc relabelling is applied
fit_kl_stan <- function(
  m,
  brms_cfg,
  pg,
  adj = NULL,
  seed = 1L,
  stan_dir = "generated_data/stan_kl",
  fix = NULL,
  labels = NULL
) {
  kcfg <- m$kcfg
  prog <- kl_compile_model(m, brms_cfg, stan_dir)
  code <- prog$code
  sha <- prog$sha
  file <- prog$file
  model <- prog$model
  threads <- prog$threads

  sd <- build_kl_latent_standata(m, pg, adj, fix = fix, threads = threads)
  pinned <- any(sd$kl_fix > 0L)
  mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms_kl (engine stan)")
  chains <- as.integer(brms_cfg$chains %||% 4L)
  strategy <- kl_init_strategy(kcfg, sd)
  init_table <- NULL
  if (identical(strategy, "optimize")) {
    ms <- kl_multistart_optimize(
      model,
      sd,
      pg,
      kcfg,
      n_starts = kcfg$init_starts %||% 8L,
      n_chains = chains,
      seed = seed,
      cores = as.integer(brms_cfg$cores %||% chains)
    )
    inits <- ms$init
    init_table <- ms$table
  } else {
    modes <- NULL
    if (identical(strategy, "given")) {
      # kernel.init_mu: starting centre per class label (lags); classes not listed start mid-range
      labs <- labels %||% pg$class_labels
      modes <- lapply(labs, function(l) {
        v <- kcfg$init_mu[[l]]
        if (is.null(v)) NULL else as.numeric(unlist(v))
      })
      message(
        "[brms_kl] given init centres: ",
        paste(
          sprintf(
            "%s=%s",
            labs,
            vapply(modes, function(x) paste(x %||% "mid", collapse = "/"), "")
          ),
          collapse = "; "
        )
      )
    }
    if (identical(strategy, "sweep")) {
      batch <- kcfg$sweep$batch %||% kcfg$pin$sweep_batch
      if (is.null(batch)) {
        stop("brms_kl: kernel.init sweep needs kernel.sweep.batch.")
      }
      labs <- labels %||% pg$class_labels
      modes <- kl_sweep_modes(
        batch,
        labs,
        outdir = kcfg$sweep$outdir %||% "generated_data/"
      )
      message(
        "[brms_kl] sweep init modes: ",
        paste(
          sprintf(
            "%s=%s",
            labs,
            vapply(modes, function(x) paste(x %||% "mid", collapse = "/"), "")
          ),
          collapse = "; "
        )
      )
      modes <- unname(modes)
    }
    inits <- make_kl_class_inits(
      sd,
      pg,
      kcfg,
      chains,
      strategy = strategy,
      seed = seed,
      mu_modes = modes
    )
  }
  # Partial inits by design (class parameters + lin coefficients only).
  old_opt <- options(cmdstanr_warn_inits = FALSE)
  on.exit(options(old_opt), add = TRUE)
  fit <- model$sample(
    data = sd,
    chains = chains,
    parallel_chains = as.integer(brms_cfg$cores %||% chains),
    iter_warmup = mcmc$warmup,
    iter_sampling = mcmc$sampling,
    seed = seed,
    init = inits %||% NULL,
    refresh = mcmc$refresh,
    adapt_delta = brms_cfg$adapt_delta %||% 0.95,
    max_treedepth = brms_cfg$max_treedepth %||% 12L,
    threads_per_chain = if (threads > 1L) threads else NULL,
    show_messages = TRUE
  )

  vars <- fit$metadata()$stan_variables
  keep <- vars[grepl(
    "^(kmu|ksd|kbeta|kc|kcon|ktheta|khs_tau|khs_lambda|khs_c2|kpi|kkappa|kgrp|ktau|krho|kfield|klogit|kw|klag|kwidth|kl_lag_prob|kre_sdc|kre_sda|kre_c|kl_log_resp|kl_log_lik_ctry|kclag|kccent|klag_sd|b_lin|sd_[0-9]+|sds_lin_[0-9]+|lp__|lprior)$",
    vars
  )]
  if (isTRUE(kcfg$keep_re_draws)) {
    keep <- c(keep, vars[grepl("^r_[0-9]+_lin_", vars)])
  }
  if (isTRUE(kcfg$gq_row_loglik)) {
    keep <- c(keep, "log_lik")
  }
  draws <- posterior::as_draws_df(fit$draws(variables = keep))
  if (sd$kl_ksd_fix == 1L) {
    # fixed width: expose it as ksd[k] for the summaries
    dd <- as.data.frame(draws)
    for (k in seq_len(pg$K)) {
      dd[[sprintf("ksd[%d]", k)]] <- sd$kl_ksd_value
    }
    draws <- posterior::as_draws_df(dd)
  }
  diag <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) {
    NULL
  })
  ctry_levels <- levels(droplevels(as.factor(m$data$country)))
  ident <- kl_identify(kcfg)
  relabel_log <- NULL
  if (identical(ident, "relabel") && pg$K >= 2L) {
    rl <- kl_relabel_draws(draws, pg$K)
    draws <- rl$draws
    relabel_log <- rl$table
  }
  naming <- NULL
  if (identical(ident, "repel") && pg$K >= 2L) {
    # exchangeable classes: every draw aligned to one labelling, then the
    # classes named after the registry (kl_repel_relabel())
    rr <- kl_repel_relabel(draws, pg$K, if (isTRUE(pg$use)) pg$A else NULL)
    draws <- rr$draws
    relabel_log <- rr$table
    naming <- rr$naming
    if (is.null(labels)) labels <- rr$labels
  }
  disc_spec <- NULL
  if (identical(ident, "lagorder")) {
    # tied centres -> order by kcon; then expose the drawn centres / width as
    # kmu / ksd so the Gaussian-form summaries apply unchanged
    if (pg$K >= 2L && !pinned) {
      rl <- kl_relabel_draws(draws, pg$K, lag = "klag")
      draws <- rl$draws
      relabel_log <- rl$table
    }
    draws <- kl_disc_rename_draws(draws, pg$K)
    disc_spec <- list(
      centres = as.numeric(sd$kl_cent),
      widths = as.numeric(sd$kl_width),
      M = sd$kl_M,
      centre_step = kcfg$discrete$centre_step %||% 1L,
      combo_prior = kcfg$discrete$combo_prior %||% "exchangeable"
    )
  }
  if (is.null(labels)) {
    labels <- pg$class_labels
    if (!identical(ident, "anchor")) {
      labels <- if (isTRUE(pg$use)) {
        kl_name_classes(draws, pg$A, pg$class_labels, pg$K)
      } else {
        sprintf("klc%d", seq_len(pg$K))
      }
    }
  }
  structure(
    list(
      draws = draws,
      program = as.character(code),
      program_sha1 = sha,
      program_file = file,
      country_levels = ctry_levels,
      class_labels = labels,
      K = pg$K,
      identify = ident,
      chains = chains,
      diagnostics = diag,
      inits = if (identical(strategy, "optimize")) NULL else inits,
      init_strategy = strategy,
      init_table = init_table,
      form = kcfg$form %||% "gaussian",
      relabel = relabel_log,
      naming = naming,
      discrete = disc_spec,
      fix = as.integer(sd$kl_fix),
      threads = threads,
      data_dims = list(N = sd$N, C = sd$kl_C, J = sd$kl_J, NL = sd$kl_NL),
      prior_partition = pg[intersect(
        c(
          "use",
          "K",
          "class_labels",
          "mu_normal",
          "mu_mean",
          "lags_found",
          "groups"
        ),
        names(pg)
      )],
      adjacency = if (!is.null(adj)) {
        adj[c("n_edges", "n_comp", "isolated")]
      } else {
        NULL
      }
    ),
    class = "kl_fit"
  )
}

#' kl_relabel_draws
#'
#' Post-hoc identification for exchangeable classes (kernel.identify =
#' relabel): in every draw, classes are re-indexed in increasing order of
#' kcon (log RR of a sustained ENSO value erf_at vs cen), and every
#' class-indexed quantity is permuted with them. Valid because the model is
#' symmetric in the class labels under relabel (symmetric spatial field, no
#' class-specific priors, no anchoring), so each permutation of a mode is an
#' equivalent mode; chains that settled on different labellings are mapped to
#' the same one. Classes whose kcon posteriors overlap stay hard to separate --
#' the returned table reports how often each permutation occurred.
#'
#' With `lag` (form discrete, identify = lagorder: the name of the drawn
#' centre variable, "klag"), classes are sorted by centre first and kcon
#' second, so only classes tied on a centre are reordered -- they are the
#' exchangeable ones (lag order already fixes the others).
#'
#' With `perm` (a draws x K matrix, perm[i, new k] = old class; identify
#' repel, kl_repel_relabel()) that permutation is applied as given.
#'
#' @param draws posterior draws_df
#' @return list(draws = relabelled draws_df, table = tibble(permutation, n, share))
kl_relabel_draws <- function(draws, K, lag = NULL, perm = NULL) {
  d <- as.data.frame(posterior::as_draws_df(draws))
  kc <- as.matrix(d[, sprintf("kcon[%d]", seq_len(K)), drop = FALSE])
  if (is.null(perm)) {
    perm <- if (is.null(lag)) {
      # perm[i, new k] = old class
      t(apply(kc, 1, order))
    } else {
      lg <- as.matrix(d[, sprintf("%s[%d]", lag, seq_len(K)), drop = FALSE])
      matrix(
        t(vapply(
          seq_len(nrow(kc)),
          function(i) order(lg[i, ], kc[i, ]),
          integer(K)
        )),
        ncol = K
      )
    }
  }
  key <- apply(perm, 1, paste, collapse = ",")
  # class position of every class-indexed family (kmu / ksd / kw / kgrp exist
  # only on the Gaussian form, i.e. under identify repel here)
  pos <- c(
    kpi = 1L,
    kcon = 1L,
    khs_tau = 1L,
    khs_c2 = 1L,
    kl_log_resp = 1L,
    klag = 1L,
    kl_lag_prob = 1L,
    kmu = 1L,
    ksd = 1L,
    kccent = 1L,
    kbeta = 2L,
    ktheta = 2L,
    khs_lambda = 2L,
    klogit = 2L,
    kfield = 2L,
    kw = 2L,
    kgrp = 2L
  )
  vars <- setdiff(names(d), c(".chain", ".iteration", ".draw"))
  fam <- sub("\\[.*$", "", vars)
  keep <- fam %in% names(pos)
  cls <- vars[keep]
  cfam <- fam[keep]
  idx <- lapply(strsplit(sub("^[^[]+\\[(.*)\\]$", "\\1", cls), ","), as.integer)
  out <- d
  for (kk in unique(key)) {
    rows <- which(key == kk)
    p <- perm[rows[1], ]
    if (identical(as.integer(p), seq_len(K))) {
      next
    }
    src <- vapply(
      seq_along(cls),
      function(i) {
        ix <- idx[[i]]
        ix[pos[[cfam[i]]]] <- p[ix[pos[[cfam[i]]]]]
        sprintf("%s[%s]", cfam[i], paste(ix, collapse = ","))
      },
      character(1)
    )
    out[rows, cls] <- d[rows, src]
  }
  tab <- tibble::as_tibble(as.data.frame(
    table(permutation = key),
    stringsAsFactors = FALSE
  ))
  tab <- dplyr::arrange(
    dplyr::mutate(dplyr::rename(tab, n = Freq), share = n / sum(n)),
    dplyr::desc(n)
  )
  list(draws = posterior::as_draws_df(out), table = tab)
}

#' All permutations of 1..n in lexicographic order (n! x n matrix).
#' @noRd
.kl_perm_matrix <- function(n) {
  if (n <= 1L) {
    return(matrix(seq_len(n), 1L))
  }
  sub <- .kl_perm_matrix(n - 1L)
  unname(do.call(
    rbind,
    lapply(seq_len(n), function(i) {
      cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow(sub)))
    })
  ))
}

#' One-to-one assignment of the rows of S (classes) to its columns (groups)
#' maximising the summed S, by exhaustive search (max(K, G) <= 8); the first
#' maximum in lexicographic order wins. @return integer(K): group of each
#' class, NA for classes left over when K > G.
#' @noRd
.kl_assign_max <- function(S) {
  K <- nrow(S)
  G <- ncol(S)
  n <- max(K, G)
  if (n > 8L) {
    stop(
      "kl_repel_relabel: more than 8 classes or registry groups (exhaustive assignment)."
    )
  }
  Sp <- matrix(0, n, n)
  Sp[seq_len(K), seq_len(G)] <- S
  P <- .kl_perm_matrix(n)
  sc <- vapply(
    seq_len(nrow(P)),
    function(i) sum(Sp[cbind(seq_len(n), P[i, ])]),
    numeric(1)
  )
  a <- P[which.max(sc), seq_len(K)]
  a[a > G] <- NA_integer_
  a
}

#' kl_repel_relabel
#'
#' Post-hoc labels of identify repel, whose classes are exchangeable (no
#' seeds, no class-specific prior, symmetric field and registry pull):
#'   1. Alignment: each draw's classes are permuted to agree with one
#'      labelling, the permutation maximising sum_k <R_i[perm[k], ], T[k, ]>
#'      (R_i the draw's K x C responsibilities, exp(kl_log_resp)); T starts
#'      at R of the draw of highest lp__ (the first on ties) and becomes the
#'      mean aligned responsibilities, until no draw changes permutation
#'      (at most 50 passes). Within- and across-chain label switching both
#'      map to the same labelling.
#'   2. Naming, with a registry (A, C x G): the one-to-one class -> group
#'      assignment maximising the posterior-mean number of a class's members
#'      in its group (S = Rbar A, K x G; .kl_assign_max()). With K > G the
#'      classes left over are klc<k>; with K < G the groups left over name no
#'      class. Without a registry every class is klc<k>.
#'   3. Order: named classes in registry group order, then the others by
#'      posterior median lag centre (kmu); klc<k> = the class in position k.
#' Deterministic given the draws.
#'
#' @param draws posterior draws_df (kl_log_resp, lp__, kmu, every class family)
#' @param A     registry membership (C x G, columns = groups) or NULL
#' @return list(draws (relabelled), labels (K), table (permutation counts,
#'   kl_relabel_draws()), naming = list(groups, assignment (group of each
#'   final class, NA = none), overlap (S, final order x groups), passes))
kl_repel_relabel <- function(draws, K, A = NULL) {
  d <- as.data.frame(posterior::as_draws_df(draws))
  n <- nrow(d)
  C <- length(grep("^kl_log_resp\\[1,", names(d)))
  R <- array(NA_real_, c(n, K, C))
  for (k in seq_len(K)) {
    R[, k, ] <- exp(as.matrix(d[,
      sprintf("kl_log_resp[%d,%d]", k, seq_len(C)),
      drop = FALSE
    ]))
  }
  Rk <- lapply(seq_len(K), function(k) matrix(R[, k, ], n, C)) # draws x C, per sampled class
  P <- .kl_perm_matrix(K)
  Tm <- matrix(R[which.max(d$lp__), , ], K, C)
  cur <- rep(0L, n)
  passes <- 0L
  repeat {
    passes <- passes + 1L
    # M[[b]][, a] = <R_i[a, ], T[b, ]>; score of permutation p = sum_k M[[k]][, p[k]]
    M <- lapply(seq_len(K), function(b) {
      vapply(Rk, function(x) as.numeric(x %*% Tm[b, ]), numeric(n))
    })
    M <- lapply(M, function(x) matrix(x, n, K))
    sc <- vapply(
      seq_len(nrow(P)),
      function(j) {
        Reduce(`+`, lapply(seq_len(K), function(k) M[[k]][, P[j, k]]))
      },
      numeric(n)
    )
    new <- max.col(matrix(sc, n), ties.method = "first")
    Tm <- matrix(0, K, C)
    for (j in unique(new)) {
      rows <- which(new == j)
      for (k in seq_len(K)) {
        Tm[k, ] <- Tm[k, ] + colSums(Rk[[P[j, k]]][rows, , drop = FALSE])
      }
    }
    Tm <- Tm / n
    if (identical(new, cur) || passes >= 50L) {
      break
    }
    cur <- new
  }
  perm <- P[cur, , drop = FALSE] # perm[i, aligned k] = sampled class
  groups <- if (is.null(A)) character(0) else colnames(A)
  assign <- if (length(groups)) {
    .kl_assign_max(Tm %*% A)
  } else {
    rep(NA_integer_, K)
  }
  kmu_d <- if (all(sprintf("kmu[%d]", seq_len(K)) %in% names(d))) {
    as.matrix(d[, sprintf("kmu[%d]", seq_len(K)), drop = FALSE])
  } else {
    NULL
  }
  kmu <- vapply(
    seq_len(K),
    function(k) {
      if (is.null(kmu_d)) {
        NA_real_
      } else {
        stats::median(kmu_d[cbind(seq_len(n), perm[, k])])
      }
    },
    numeric(1)
  )
  ord <- order(is.na(assign), assign, kmu)
  fin <- perm[, ord, drop = FALSE]
  rl <- kl_relabel_draws(draws, K, perm = fin)
  a_f <- assign[ord]
  labels <- ifelse(is.na(a_f), sprintf("klc%d", seq_len(K)), groups[a_f])
  S <- NULL
  if (length(groups)) {
    S <- (Tm %*% A)[ord, , drop = FALSE]
    rownames(S) <- labels
  }
  list(
    draws = rl$draws,
    labels = labels,
    table = rl$table,
    naming = list(
      groups = groups,
      assignment = if (length(groups)) groups[a_f] else rep(NA_character_, K),
      overlap = S,
      passes = passes
    )
  )
}

#' kl_disc_rename_draws -- form discrete: the drawn class centres klag[k]
#' become kmu[k] and the drawn shared width kwidth[1] becomes ksd[k] (every
#' class), so summarise_kernel_draws() and the reports read them as the
#' Gaussian form's kernel parameters.
kl_disc_rename_draws <- function(draws, K) {
  d <- as.data.frame(posterior::as_draws_df(draws))
  for (k in seq_len(K)) {
    d[[sprintf("kmu[%d]", k)]] <- d[[sprintf("klag[%d]", k)]]
    d[[sprintf("ksd[%d]", k)]] <- d[["kwidth[1]"]]
  }
  d <- d[, setdiff(names(d), sprintf("klag[%d]", seq_len(K))), drop = FALSE]
  posterior::as_draws_df(d)
}

#' kl_name_classes
#'
#' Names for classes identified by rank (kernel.identify kmu / erf) from a
#' prior partition: the one-to-one class -> registry-group assignment that
#' maximises the posterior-mean number of member countries in their prior
#' group (sum_c r_ck A_cg), by exhaustive search over assignments (K <= 7).
#' Classes left without a group keep klc<k>.
#' @param draws draws_df with kl_log_resp; A C x K one-hot prior matrix whose
#'   columns are the registry groups (`groups`, first ones) then extra classes
#' @return character(K)
kl_name_classes <- function(draws, A, groups, K) {
  C <- nrow(A)
  dm <- posterior::as_draws_matrix(posterior::subset_draws(
    draws,
    variable = "kl_log_resp"
  ))
  R <- matrix(
    vapply(
      seq_len(K * C),
      function(i) {
        idx <- arrayInd(i, c(K, C))
        mean(exp(dm[, sprintf("kl_log_resp[%d,%d]", idx[1], idx[2])]))
      },
      numeric(1)
    ),
    K,
    C
  )
  ng <- sum(colSums(A) > 0)
  S <- R %*% A[, seq_len(ng), drop = FALSE] # K x ng
  if (K > 7L) {
    stop("kl_name_classes: K > 7 not supported by exhaustive assignment.")
  }
  perms <- function(v) {
    if (length(v) <= 1) {
      list(v)
    } else {
      do.call(
        c,
        lapply(seq_along(v), function(i) {
          lapply(perms(v[-i]), function(p) c(v[i], p))
        })
      )
    }
  }
  best <- NULL
  best_s <- -Inf
  for (p in perms(seq_len(K))) {
    # p[g] = class given group g
    sc <- sum(S[cbind(p[seq_len(ng)], seq_len(ng))])
    if (sc > best_s) {
      best_s <- sc
      best <- p
    }
  }
  lab <- sprintf("klc%d", seq_len(K))
  lab[best[seq_len(ng)]] <- groups[seq_len(ng)]
  lab
}

# ---------------------------------------------------------------------------
# Pinned arrangements (kernel.pin; form discrete)
# ---------------------------------------------------------------------------
#
# Under form discrete a class is a position in lag order, but which class (its
# curve and member countries) fills which position is left to the sampler: a
# class can settle in a position whose lag range excludes its true lag and stay
# there (synthetic case H: the lag-8 class in position 1, fitted at lag -8 with
# the opposite sign). kernel.pin fixes the class -> position correspondence by
# seed anchoring and tests the possible correspondences against each other:
#   1. pilots        each group's lag distribution, from one of
#                    order_source "pilot": a single-class discrete fit per
#                      prior group (registry kernel.prior_grouping.source, that
#                      group's countries only) -> its lag posterior;
#                    order_source "sweep": an existing single-lag sweep's LOO
#                      ranking per group (05_fit_postprocess_loo_lags.R caches,
#                      read with loo_lags.R::read_latest_ranking()) -> weights
#                      exp(elpd_loo - max) over its lags (no extra fits; the
#                      sweep's own model, not this one);
#   2. arrangements  orders of the groups along the lag positions, with their
#                    probability when the pilot lags are drawn independently
#                    (ties broken at random); "pilot" keeps the most probable
#                    (up to `cover` probability, at most `max_runs`), "all" every
#                    permutation;
#   3. runs          one latent-class fit per arrangement sigma: the seed
#                    countries of group sigma[k] are fixed to position k
#                    (kl_fix), every other country's class is inferred; runs in
#                    parallel;
#   4. combine       stacking weights on each run's country-level PSIS-LOO
#                    (kl_log_lik_ctry). Runs differ only in their seeds, so a
#                    wrong order shows up as seeds forced into a lag position
#                    that fits them badly.
# Seeds default to each group's `seeds_n` countries with the most events
# (kernel.pin.seeds: explicit iso3 per group overrides).

#' kl_pin_cfg -- validated kernel.pin block (enable FALSE -> other fields unused).
kl_pin_cfg <- function(kcfg) {
  pn <- kcfg$pin %||% list()
  pc <- pn$pilot %||% list()
  out <- list(
    enable = isTRUE(pn$enable),
    seeds_n = as.integer(pn$seeds_n %||% 1L),
    seeds = pn[["seeds"]],
    arrangements = pn$arrangements %||% "pilot",
    cover = as.numeric(pn$cover %||% 0.95),
    max_runs = as.integer(pn$max_runs %||% 6L),
    parallel_runs = as.integer(pn$parallel_runs %||% 2L),
    order_source = pn$order_source %||% "pilot",
    sweep_batch = pn[["sweep_batch"]],
    sweep_outdir = pn$sweep_outdir %||% "generated_data/",
    pilot = list(
      chains = as.integer(pc$chains %||% 2L),
      warmup = as.integer(pc$warmup %||% 200L),
      sampling = as.integer(pc$sampling %||% 200L)
    )
  )
  if (!out$enable) {
    return(out)
  }
  if (!identical(kcfg$form %||% "gaussian", "discrete")) {
    stop("brms_kl: kernel.pin needs kernel.form discrete.")
  }
  if (is.null(kcfg$prior_grouping$source)) {
    stop(
      "brms_kl: kernel.pin needs kernel.prior_grouping.source (the groups to pin)."
    )
  }
  if (!out$arrangements %in% c("pilot", "all")) {
    stop("brms_kl: kernel.pin.arrangements must be pilot | all.")
  }
  if (!out$order_source %in% c("pilot", "sweep")) {
    stop("brms_kl: kernel.pin.order_source must be pilot | sweep.")
  }
  if (
    identical(out$order_source, "sweep") &&
      is.null(out$sweep_batch) &&
      is.null(kcfg$prior_grouping$deploy_config)
  ) {
    stop(
      "brms_kl: kernel.pin.order_source sweep needs kernel.pin.sweep_batch ",
      "(or kernel.prior_grouping.deploy_config, whose fits.occ_dir is used)."
    )
  }
  if (is.na(out$seeds_n) || out$seeds_n < 1L) {
    stop("brms_kl: kernel.pin.seeds_n must be >= 1.")
  }
  out
}

#' All permutations of a vector (list of vectors). @noRd
.kl_perms <- function(v) {
  if (length(v) <= 1L) {
    list(v)
  } else {
    do.call(
      c,
      lapply(seq_along(v), function(i) {
        lapply(.kl_perms(v[-i]), function(p) c(v[i], p))
      })
    )
  }
}

#' kl_pin_seeds -- seed countries per group: `seeds[[g]]` when given, else the
#' group's n fitted countries with the most events (sum of y; ties by iso3).
#' @return named list group -> iso3 vector
kl_pin_seeds <- function(dat, reg, groups, n = 1L, seeds = NULL) {
  ev <- tapply(dat$y, as.character(dat$country), sum)
  out <- lapply(groups, function(g) {
    if (!is.null(seeds[[g]])) {
      s <- intersect(as.character(unlist(seeds[[g]])), names(ev))
      if (!length(s)) {
        stop(
          "kl_pin_seeds: no fitted country among the seeds given for group ",
          g,
          "."
        )
      }
      return(s)
    }
    cs <- intersect(reg$iso3[reg$group == g], names(ev))
    utils::head(cs[order(-ev[cs], cs)], n)
  })
  names(out) <- groups
  out
}

#' kl_pin_fix -- per-country fixed class for arrangement `sigma` (sigma[k] =
#' group at position k): seeds of group sigma[k] -> k, other countries 0.
kl_pin_fix <- function(country_levels, seeds, sigma) {
  fix <- integer(length(country_levels))
  for (k in seq_along(sigma)) {
    fix[match(seeds[[sigma[k]]], country_levels)] <- k
  }
  fix
}

#' kl_pin_pilots
#'
#' One single-class discrete fit per group on that group's countries (no prior
#' partition, no spatial term, region level dropped, MCMC from kernel.pin.pilot),
#' programs compiled up front, fits in parallel.
#' @param dat Poisson view before assembly (prepare_occ_brms_data()$pois_dat)
#' @return list per group: group, countries, lags, lag_prob, kcon (mean, q2.5, q97.5), fit
kl_pin_pilots <- function(
  opt,
  dat,
  brms_cfg,
  baseline_block,
  reg,
  groups,
  seed = 1L,
  stan_dir = "generated_data/stan_kl",
  cores = 2L
) {
  pc <- kl_pin_cfg(resolve_kernel_cfg(opt))$pilot
  bcfg <- utils::modifyList(
    brms_cfg,
    list(
      chains = pc$chains,
      cores = pc$chains,
      warmup = pc$warmup,
      sampling = pc$sampling
    )
  )
  opt_p <- opt
  k <- opt_p$occurrence$kernel
  k$classes <- 1L
  k$prior_grouping <- list(source = NULL)
  k$spatial <- list(enable = FALSE)
  k$pin <- list(enable = FALSE)
  k$drop_region <- TRUE
  opt_p$occurrence$kernel <- k
  kc <- resolve_kernel_cfg(opt_p)
  models <- lapply(groups, function(g) {
    sub <- dat[
      as.character(dat$country) %in% reg$iso3[reg$group == g],
      ,
      drop = FALSE
    ]
    if (is.factor(sub$country)) {
      sub$country <- droplevels(sub$country)
    }
    kl_sort_by_country(suppressMessages(assemble_kl_model(
      opt_p,
      sub,
      bcfg,
      baseline_block,
      kc
    )))
  })
  for (mm in models) {
    # compile once, serially
    f <- kl_program_file(mm, stan_dir)$file
    kl_with_compile_lock(f, cmdstanr::cmdstan_model(f))
  }
  one <- function(i) {
    mm <- models[[i]]
    lv <- levels(droplevels(as.factor(mm$data$country)))
    f <- fit_kl_stan(
      mm,
      bcfg,
      kl_prior_partition(kc, lv),
      NULL,
      seed = seed + i,
      stan_dir = stan_dir
    )
    dm <- posterior::as_draws_matrix(f$draws)
    cn <- f$discrete$centres
    kco <- dm[, "kcon[1]"]
    list(
      group = groups[i],
      countries = lv,
      lags = cn,
      lag_prob = colMeans(dm[,
        sprintf("kl_lag_prob[1,%d]", seq_along(cn)),
        drop = FALSE
      ]),
      kcon = c(
        mean = mean(kco),
        stats::quantile(kco, c(0.025, 0.975), names = FALSE)
      ),
      fit = f
    )
  }
  res <- parallel::mclapply(
    seq_along(groups),
    one,
    mc.preschedule = FALSE,
    mc.cores = max(1L, min(length(groups), cores %/% pc$chains))
  )
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) {
    stop(
      "kl_pin_pilots: pilot fit failed for group(s) ",
      paste(groups[bad], collapse = ", "),
      ": ",
      as.character(res[bad][[1]])
    )
  }
  names(res) <- groups
  res
}

#' kl_pin_sweep_orders
#'
#' Group lag distributions from an existing single-lag sweep (order_source
#' "sweep"): the latest LOO ranking cache of each group in `batch`
#' (loo_lags.R::read_latest_ranking(), under `<outdir>/loo_compare/<batch>/`),
#' turned into weights exp(elpd_loo - max) over the sweep's lags. These weights
#' ignore the SE of the differences (the caches hold totals only), so they are
#' somewhat overconfident -- raise kernel.pin.cover to compensate. Same shape as
#' kl_pin_pilots() (no fit; kcon NA).
#' @return list per group: group, lags, lag_prob, kcon, elpd_best, source
kl_pin_sweep_orders <- function(
  batch,
  groups,
  outdir = "generated_data/",
  lags = NULL
) {
  res <- lapply(groups, function(g) {
    r <- read_latest_ranking(batch, group = g, opt = list(outdir = outdir))
    if (!all(c("single_lag", "elpd_loo") %in% names(r))) {
      stop(
        "kl_pin_sweep_orders: ranking for ",
        g,
        " lacks single_lag / elpd_loo."
      )
    }
    o <- order(r$single_lag)
    lg <- as.numeric(r$single_lag[o])
    el <- as.numeric(r$elpd_loo[o])
    if (!is.null(lags) && !all(lags %in% lg)) {
      message(sprintf(
        "[kl_pin] sweep %s / %s covers lags %g..%g, the model %g..%g",
        batch,
        g,
        min(lg),
        max(lg),
        min(lags),
        max(lags)
      ))
    }
    w <- exp(el - max(el))
    list(
      group = g,
      lags = lg,
      lag_prob = w / sum(w),
      kcon = c(NA_real_, NA_real_, NA_real_),
      elpd_best = max(el),
      source = sprintf("sweep %s", batch)
    )
  })
  names(res) <- groups
  res
}

#' kl_pin_arrangements
#'
#' Orders of the groups along the lag positions and their probability under
#' independent pilot lag posteriors (Monte Carlo, ties broken at random).
#' mode "pilot": the most probable orders up to cumulative probability `cover`
#' (at most `max_runs`); "all": every permutation.
#' @return tibble(arrangement "g1>g2>...", prob, sigma (list))
kl_pin_arrangements <- function(
  pilots,
  mode = "pilot",
  cover = 0.95,
  max_runs = 6L,
  n_mc = 20000L,
  seed = 1L
) {
  groups <- vapply(pilots, `[[`, character(1), "group")
  if (length(groups) > 6L) {
    stop("kl_pin_arrangements: more than 6 groups (720 orders) not supported.")
  }
  set.seed(seed)
  L <- vapply(
    pilots,
    function(p) {
      sample(p$lags, n_mc, replace = TRUE, prob = p$lag_prob) +
        stats::runif(n_mc, 0, 1e-3)
    },
    numeric(n_mc)
  )
  L <- matrix(L, nrow = n_mc)
  key <- apply(L, 1, function(r) paste(groups[order(r)], collapse = ">"))
  pr <- table(key) / n_mc
  all <- vapply(.kl_perms(groups), paste, character(1), collapse = ">")
  arr <- tibble::tibble(
    arrangement = all,
    prob = as.numeric(ifelse(all %in% names(pr), pr[all], 0))
  )
  arr <- arr[order(-arr$prob, arr$arrangement), , drop = FALSE]
  if (identical(mode, "pilot")) {
    n <- min(
      max_runs,
      which(cumsum(arr$prob) >= cover - 1e-12)[1],
      sum(arr$prob > 0)
    )
    arr <- arr[seq_len(max(1L, n)), , drop = FALSE]
  }
  arr$sigma <- strsplit(arr$arrangement, ">", fixed = TRUE)
  arr
}

#' kl_pin_fit_one -- the latent-class fit of arrangement i (seeds of group
#' arr$sigma[[i]][k] fixed to class k, labels = sigma, seed + 100 i): one run of
#' kl_pin_fit_arrangements(), or one SLURM array task (--pin_stage run).
kl_pin_fit_one <- function(
  m,
  brms_cfg,
  pg,
  adj,
  seeds,
  arr,
  i,
  seed = 1L,
  stan_dir = "generated_data/stan_kl"
) {
  lv <- levels(droplevels(as.factor(m$data$country)))
  sigma <- arr$sigma[[i]]
  f <- fit_kl_stan(
    m,
    brms_cfg,
    pg,
    adj,
    seed = seed + 100L * i,
    stan_dir = stan_dir,
    fix = kl_pin_fix(lv, seeds, sigma),
    labels = sigma
  )
  f$arrangement <- arr$arrangement[i]
  f$pin_task <- i
  f
}

#' kl_pin_run_file -- per-arrangement fit file of the split (SLURM array)
#' path: the kl_pin artefact path with _arr<i> before the extension.
kl_pin_run_file <- function(pin_file, i) {
  vapply(
    as.integer(i),
    function(k) sub("\\.rds$", sprintf("_arr%02d.rds", k), pin_file),
    character(1)
  )
}

#' kl_pin_fit_arrangements -- one latent-class fit per arrangement (seeds
#' pinned, labels = sigma), the shared program compiled once, runs in parallel.
#' @return list of kl_fit (each with $arrangement)
kl_pin_fit_arrangements <- function(
  m,
  brms_cfg,
  pg,
  adj,
  seeds,
  arr,
  seed = 1L,
  stan_dir = "generated_data/stan_kl",
  parallel_runs = 2L
) {
  f <- kl_program_file(m, stan_dir)$file
  kl_with_compile_lock(f, cmdstanr::cmdstan_model(f))
  one <- function(i) {
    kl_pin_fit_one(
      m,
      brms_cfg,
      pg,
      adj,
      seeds,
      arr,
      i,
      seed = seed,
      stan_dir = stan_dir
    )
  }
  res <- parallel::mclapply(
    seq_len(nrow(arr)),
    one,
    mc.preschedule = FALSE,
    mc.cores = max(1L, min(nrow(arr), parallel_runs))
  )
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) {
    stop(
      "kl_pin_fit_arrangements: run(s) failed: ",
      paste(arr$arrangement[bad], collapse = ", "),
      ": ",
      as.character(res[bad][[1]])
    )
  }
  res
}

#' kl_pin_combine -- country-level PSIS-LOO per run and stacking weights.
#' @return tibble(arrangement, elpd_loo, se, elpd_diff (vs best), weight,
#'   n_high_k (Pareto k > 0.7), max_rhat_kcon)
kl_pin_combine <- function(fits) {
  loos <- lapply(fits, function(f) suppressWarnings(kl_loo_country(f)))
  lpd <- vapply(
    loos,
    function(l) l$pointwise[, "elpd_loo"],
    numeric(nrow(loos[[1]]$pointwise))
  )
  lpd <- matrix(lpd, ncol = length(fits))
  w <- if (length(fits) == 1L) 1 else as.numeric(loo::stacking_weights(lpd))
  el <- vapply(
    loos,
    function(l) l$estimates["elpd_loo", "Estimate"],
    numeric(1)
  )
  rh <- vapply(
    fits,
    function(f) {
      max(vapply(
        sprintf("kcon[%d]", seq_len(f$K)),
        function(v) {
          posterior::rhat(posterior::extract_variable_matrix(f$draws, v))
        },
        numeric(1)
      ))
    },
    numeric(1)
  )
  tibble::tibble(
    arrangement = vapply(fits, `[[`, character(1), "arrangement"),
    elpd_loo = el,
    se = vapply(loos, function(l) l$estimates["elpd_loo", "SE"], numeric(1)),
    elpd_diff = el - max(el),
    weight = w,
    n_high_k = vapply(
      loos,
      function(l) sum(l$diagnostics$pareto_k > 0.7),
      integer(1)
    ),
    max_rhat_kcon = rh
  )
}

#' kl_pin_plan
#'
#' kernel.pin up to the runs: seeds, group lag distributions (pilot fits or
#' sweep rankings) and the arrangements to fit. Pilot fit objects are dropped
#' (their summaries are kept), so the plan is small enough to hand to SLURM
#' array tasks (--pin_stage plan).
#' @return list(groups, seeds, pilots, arrangements, order_source)
kl_pin_plan <- function(
  opt,
  dat,
  brms_cfg,
  baseline_block,
  m,
  pg,
  seed = 1L,
  stan_dir = "generated_data/stan_kl",
  groupings_dir = "analysis/groupings",
  cores = 4L
) {
  kcfg <- m$kcfg
  pc <- kl_pin_cfg(kcfg)
  reg <- load_grouping(kcfg$prior_grouping$source, groupings_dir)
  groups <- pg$class_labels[pg$class_labels %in% reg$group]
  if (length(groups) != pg$K) {
    stop(
      "kl_pin_plan: kernel.pin needs one class per prior group (K = ",
      pg$K,
      ", groups: ",
      paste(groups, collapse = ", "),
      ")."
    )
  }
  seeds <- kl_pin_seeds(m$data, reg, groups, pc$seeds_n, pc[["seeds"]])
  message(
    "[kl_pin] seeds: ",
    paste(
      sprintf("%s=%s", groups, vapply(seeds, paste, "", collapse = "+")),
      collapse = "; "
    )
  )
  if (identical(pc$order_source, "sweep")) {
    batch <- pc$sweep_batch %||%
      yaml::read_yaml(kcfg$prior_grouping$deploy_config)$deploy$fits$occ_dir
    pilots <- kl_pin_sweep_orders(
      batch,
      groups,
      outdir = pc$sweep_outdir,
      lags = m$meta$lags
    )
  } else {
    pilots <- kl_pin_pilots(
      opt,
      dat,
      brms_cfg,
      baseline_block,
      reg,
      groups,
      seed = seed,
      stan_dir = stan_dir,
      cores = cores
    )
  }
  for (p in pilots) {
    message(sprintf(
      "[kl_pin] %s %s: lag mode %g, P(mode) %.2f%s",
      pc$order_source,
      p$group,
      p$lags[which.max(p$lag_prob)],
      max(p$lag_prob),
      if (is.na(p$kcon[1])) {
        ""
      } else {
        sprintf(", kcon %.2f [%.2f, %.2f]", p$kcon[1], p$kcon[2], p$kcon[3])
      }
    ))
  }
  arr <- kl_pin_arrangements(
    pilots,
    pc$arrangements,
    pc$cover,
    pc$max_runs,
    seed = seed
  )
  message(
    "[kl_pin] arrangements: ",
    paste(sprintf("%s (%.3f)", arr$arrangement, arr$prob), collapse = "; ")
  )
  list(
    groups = groups,
    seeds = seeds,
    pilots = lapply(pilots, function(p) p[setdiff(names(p), "fit")]),
    arrangements = arr,
    order_source = pc$order_source
  )
}

#' kl_pin_finish -- stacking over the arrangement fits of a plan.
#' @return list(seeds, pilots, arrangements, fits, table, best)
kl_pin_finish <- function(plan, fits) {
  tab <- kl_pin_combine(fits)
  arr <- plan$arrangements
  tab$pilot_prob <- arr$prob[match(tab$arrangement, arr$arrangement)]
  list(
    seeds = plan$seeds,
    pilots = plan$pilots,
    arrangements = arr,
    fits = fits,
    table = tab,
    best = which.max(tab$weight)
  )
}

#' kl_pin_run
#'
#' kernel.pin end to end in one process: kl_pin_plan(), the runs in parallel
#' (kl_pin_fit_arrangements()), kl_pin_finish(). The SLURM split runs the same
#' three pieces as plan / array / combine jobs (04_fit_occ_brms_kl.R --pin_stage).
#' @param dat Poisson view before assembly (pilots re-assemble per group)
#' @param m,pg,adj the full latent-class model (sorted), prior partition, adjacency
#' @return list(seeds, pilots, arrangements, fits, table, best (index of the
#'   highest stacking weight))
kl_pin_run <- function(
  opt,
  dat,
  brms_cfg,
  baseline_block,
  m,
  pg,
  adj = NULL,
  seed = 1L,
  stan_dir = "generated_data/stan_kl",
  groupings_dir = "analysis/groupings",
  cores = 4L
) {
  plan <- kl_pin_plan(
    opt,
    dat,
    brms_cfg,
    baseline_block,
    m,
    pg,
    seed = seed,
    stan_dir = stan_dir,
    groupings_dir = groupings_dir,
    cores = cores
  )
  fits <- kl_pin_fit_arrangements(
    m,
    brms_cfg,
    pg,
    adj,
    plan$seeds,
    plan$arrangements,
    seed = seed,
    stan_dir = stan_dir,
    parallel_runs = kl_pin_cfg(m$kcfg)$parallel_runs
  )
  kl_pin_finish(plan, fits)
}

# Summaries -------------------------------------------------------------------

#' Draws matrix of a K x C (or C) indexed Stan variable, reshaped per draw.
#' @noRd
.kl_draw_array <- function(fit, var, dims) {
  dm <- posterior::as_draws_matrix(posterior::subset_draws(
    fit$draws,
    variable = var
  ))
  idx <- as.matrix(expand.grid(lapply(dims, seq_len)))
  cols <- sprintf("%s[%s]", var, apply(idx, 1, paste, collapse = ","))
  arr <- dm[, cols, drop = FALSE]
  array(as.numeric(arr), dim = c(nrow(arr), dims))
}

#' kl_responsibilities
#'
#' Posterior class responsibilities per country (exp of kl_log_resp).
#' @return tibble(country, class, label, mean, sd, q2.5, q97.5)
kl_responsibilities <- function(fit) {
  R <- exp(.kl_draw_array(
    fit,
    "kl_log_resp",
    c(fit$K, length(fit$country_levels))
  ))
  grid <- expand.grid(class = seq_len(fit$K), c = seq_along(fit$country_levels))
  vals <- t(vapply(
    seq_len(nrow(grid)),
    function(i) {
      x <- R[, grid$class[i], grid$c[i]]
      c(
        mean(x),
        stats::sd(x),
        stats::quantile(x, c(0.025, 0.975), names = FALSE)
      )
    },
    numeric(4)
  ))
  tibble::tibble(
    country = fit$country_levels[grid$c],
    class = grid$class,
    label = fit$class_labels[grid$class],
    mean = vals[, 1],
    sd = vals[, 2],
    q2.5 = vals[, 3],
    q97.5 = vals[, 4]
  )
}

#' kl_comembership
#'
#' Posterior probability that two different countries share a class: per
#' draw sum_k r_ck r_dk (r = responsibilities; the classes of distinct
#' countries are conditionally independent given the parameters), averaged
#' over draws (or within each chain when by_chain = TRUE); diagonal = 1.
#' Label-invariant.
#' @return C x C matrix (dimnames = countries), or a list of them per chain
kl_comembership <- function(fit, by_chain = FALSE) {
  C <- length(fit$country_levels)
  R <- exp(.kl_draw_array(fit, "kl_log_resp", c(fit$K, C)))
  comem <- function(rows) {
    M <- matrix(
      0,
      C,
      C,
      dimnames = list(fit$country_levels, fit$country_levels)
    )
    for (d in rows) {
      M <- M + crossprod(matrix(R[d, , ], nrow = fit$K))
    } # K x C, kept 2-D at K = 1
    M <- M / length(rows)
    diag(M) <- 1 # a country always shares its own class (sum_k r_ck^2 < 1 is not co-membership)
    M
  }
  if (!by_chain) {
    return(comem(seq_len(dim(R)[1])))
  }
  ch <- fit$draws$.chain
  lapply(split(seq_along(ch), ch), comem)
}

#' kl_chain_agreement
#'
#' Pairwise agreement of the per-chain co-membership matrices.
#' @return tibble(chain_a, chain_b, max_abs_diff, mean_abs_diff)
kl_chain_agreement <- function(fit) {
  M <- kl_comembership(fit, by_chain = TRUE)
  if (length(M) < 2L) {
    return(tibble::tibble(
      chain_a = integer(0),
      chain_b = integer(0),
      max_abs_diff = numeric(0),
      mean_abs_diff = numeric(0)
    ))
  }
  pr <- utils::combn(names(M), 2)
  tibble::tibble(
    chain_a = as.integer(pr[1, ]),
    chain_b = as.integer(pr[2, ]),
    max_abs_diff = apply(pr, 2, function(p) max(abs(M[[p[1]]] - M[[p[2]]]))),
    mean_abs_diff = apply(pr, 2, function(p) mean(abs(M[[p[1]]] - M[[p[2]]])))
  )
}

#' kl_map_partition
#'
#' Point partition of the countries. "argmax": each country's class with the
#' highest mean responsibility (labels = class labels: the anchoring registry
#' groups under a prior partition, klc<k> otherwise). "hclust": average-link
#' clustering of 1 - co-membership cut at K groups (label-invariant), each
#' cluster labelled by the majority argmax class of its members.
#' @return tibble(iso3, group)
kl_map_partition <- function(
  fit,
  method = c("argmax", "hclust"),
  resp = kl_responsibilities(fit)
) {
  method <- match.arg(method)
  am <- resp %>%
    dplyr::group_by(country) %>%
    dplyr::slice_max(mean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  am <- am[match(fit$country_levels, am$country), ]
  if (method == "argmax") {
    return(tibble::tibble(iso3 = am$country, group = am$label))
  }
  M <- kl_comembership(fit)
  cl <- stats::cutree(
    stats::hclust(stats::as.dist(1 - M), method = "average"),
    k = min(fit$K, length(fit$country_levels))
  )
  lab <- vapply(
    split(am$label, cl),
    function(x) names(which.max(table(x))),
    character(1)
  )
  tibble::tibble(
    iso3 = fit$country_levels,
    group = unname(lab[as.character(cl)])
  )
}

#' kl_subset_chains
#'
#' A kl_fit restricted to some of its chains (e.g. the ones that converged to
#' the same mode), for downstream use such as utils/make_kl_grouping.R. Draws
#' are kept as they are; chains are renumbered 1..n in increasing source order
#' (posterior::subset_draws()) and `source_chains` records the original ids.
#' Per-chain fields (diagnostics, inits) are subset to match.
#'
#' @param fit    kl_fit (fit_kl_stan() output, as saved by 04_fit_occ_brms_kl.R)
#' @param chains integer ids of the chains to keep
kl_subset_chains <- function(fit, chains) {
  chains <- sort(unique(as.integer(chains)))
  have <- sort(unique(fit$draws$.chain))
  if (!length(chains) || !all(chains %in% have)) {
    stop(
      "kl_subset_chains: chains ",
      paste(setdiff(chains, have), collapse = ", "),
      " not in the fit (it has ",
      paste(have, collapse = ", "),
      ")."
    )
  }
  n0 <- length(have)
  fit$draws <- posterior::subset_draws(fit$draws, chain = chains)
  fit$diagnostics <- lapply(fit$diagnostics, function(x) {
    if (length(x) == n0) x[chains] else x
  })
  if (length(fit$inits) == n0) {
    fit$inits <- fit$inits[chains]
  }
  fit$source_chains <- c(fit$source_chains %||% have)[chains]
  fit$chains <- length(chains)
  fit
}

#' write_kl_partition -- iso3,group CSV in the groupings-registry format,
#' validated by reading it back with load_grouping().
write_kl_partition <- function(partition, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    partition[, c("iso3", "group")],
    path,
    row.names = FALSE,
    quote = FALSE
  )
  invisible(load_grouping(sub("\\.csv$", "", basename(path)), dirname(path)))
}

#' kl_loo_country
#'
#' PSIS-LOO over countries on kl_log_lik_ctry (leave-one-country-out unit).
#' Each country is a large block, so Pareto-k is often high; treat as
#' descriptive and prefer exact refits for decisions.
kl_loo_country <- function(fit) {
  ll <- .kl_draw_array(fit, "kl_log_lik_ctry", length(fit$country_levels))
  ll <- matrix(ll, nrow = dim(ll)[1])
  r_eff <- loo::relative_eff(exp(ll), chain_id = fit$draws$.chain)
  out <- loo::loo(ll, r_eff = r_eff)
  attr(
    out,
    "note"
  ) <- "country-block PSIS-LOO: high Pareto-k expected; descriptive only"
  out
}

#' summarise_kl_dlnm_class
#'
#' Summaries of one class of the lag-basis form (kernel.form dlnm):
#'   lagresp  lag-response at the ERF contrast: g_l = sum_j (B_j(erf_at) -
#'            B_j(cen)) theta_{j,l}, the log RR at lag l of ENSO at erf_at vs
#'            cen (mean, 95% interval, P(g_l > 0));
#'   erf      sustained-exposure log RR over the exposure grid (from kbeta);
#'   params   kcon (= sum_l g_l), peak lag (largest |posterior-mean g_l|), effective
#'            width (inverse participation ratio (sum |g|)^2 / sum g^2: 1 = one
#'            lag, NL = flat), and the horseshoe global scale when present.
#' @return list(params, lagresp, erf) of tibbles
summarise_kl_dlnm_class <- function(
  dm,
  k,
  meta,
  cen = 0,
  erf_at = 1.5,
  probs = c(0.025, 0.5, 0.975),
  n_grid = 50L
) {
  J <- meta$df
  lags <- meta$lags
  NL <- length(lags)
  v <- as.numeric(.kl_basis_centred(meta, erf_at, cen))
  TH <- dm[, sprintf("ktheta[%d,%d]", seq_len(NL * J), k), drop = FALSE]
  G <- vapply(
    seq_len(NL),
    function(l) TH[, (l - 1) * J + seq_len(J), drop = FALSE] %*% v,
    numeric(nrow(TH))
  )
  G <- matrix(G, nrow = nrow(TH))
  qs <- function(m) {
    m <- as.matrix(m)
    out <- cbind(
      mean = colMeans(m),
      t(apply(m, 2, stats::quantile, probs = probs, names = FALSE))
    )
    colnames(out) <- c("mean", sprintf("q%s", probs * 100))
    tibble::as_tibble(out)
  }
  lagresp <- dplyr::bind_cols(
    tibble::tibble(lag = lags),
    qs(G),
    tibble::tibble(p_pos = colMeans(G > 0))
  )
  ipr <- rowSums(abs(G))^2 / pmax(rowSums(G^2), 1e-300)
  pars <- cbind(kcon = dm[, sprintf("kcon[%d]", k)], width_ipr = ipr)
  hs <- sprintf("khs_tau[%d]", k)
  if (hs %in% colnames(dm)) {
    pars <- cbind(pars, hs_tau = dm[, hs])
  }
  params <- dplyr::bind_cols(tibble::tibble(par = colnames(pars)), qs(pars))
  # Peak of the posterior-mean lag-response (largest |mean g_l|). The per-draw
  # argmax is not used: a weakly constrained boundary lag with a wide interval
  # wins it in many draws and drags the summary to the edge of the range.
  gm <- colMeans(G)
  params <- dplyr::bind_rows(
    params,
    tibble::tibble(par = "peak_lag", mean = lags[which.max(abs(gm))])
  )
  beta <- dm[, sprintf("kbeta[%d,%d]", seq_len(J), k), drop = FALSE]
  at <- seq(meta$range[1], meta$range[2], length.out = n_grid)
  Bc <- .kl_basis_centred(meta, at, cen)
  erf <- dplyr::bind_cols(tibble::tibble(exposure = at), qs(beta %*% t(Bc)))
  list(params = params, lagresp = lagresp, erf = erf)
}

#' summarise_kl_latent
#'
#' Per-class kernel summaries (summarise_kernel_draws()) plus class weights,
#' kappa and tau.
#' @return list(classes = list per class of params/weights/erf, global = tibble)
summarise_kl_latent <- function(fit, meta, cen = 0, erf_at = 1.5) {
  K <- fit$K
  J <- meta$df
  dm <- posterior::as_draws_matrix(fit$draws)
  if (identical(fit$form, "dlnm")) {
    classes <- lapply(seq_len(K), function(k) {
      s <- summarise_kl_dlnm_class(dm, k, meta, cen = cen, erf_at = erf_at)
      lapply(s, function(t) {
        dplyr::mutate(t, class = k, label = fit$class_labels[k], .before = 1)
      })
    })
    gv <- intersect(
      colnames(dm),
      c(sprintf("kpi[%d]", seq_len(K)), "kkappa", "ktau[1]", "krho[1]")
    )
    global <- tibble::tibble(
      par = gv,
      mean = colMeans(dm[, gv, drop = FALSE]),
      q2.5 = apply(dm[, gv, drop = FALSE], 2, stats::quantile, 0.025),
      q97.5 = apply(dm[, gv, drop = FALSE], 2, stats::quantile, 0.975)
    )
    return(list(classes = classes, global = global))
  }
  classes <- lapply(seq_len(K), function(k) {
    beta <- dm[, sprintf("kbeta[%d,%d]", seq_len(J), k), drop = FALSE]
    s <- summarise_kernel_draws(
      dm[, sprintf("kmu[%d]", k)],
      dm[, sprintf("ksd[%d]", k)],
      beta,
      meta,
      cen
    )
    if (!is.null(fit$discrete)) {
      # form discrete: posterior probability of each candidate centre
      # (Rao-Blackwellised: mean over draws of P(centre | parameters))
      cn <- fit$discrete$centres
      s$lagprob <- tibble::tibble(
        lag = cn,
        prob = colMeans(dm[,
          sprintf("kl_lag_prob[%d,%d]", k, seq_along(cn)),
          drop = FALSE
        ])
      )
    }
    lapply(s, function(t) {
      dplyr::mutate(t, class = k, label = fit$class_labels[k], .before = 1)
    })
  })
  gv <- intersect(
    colnames(dm),
    c(
      sprintf("kpi[%d]", seq_len(K)),
      "kkappa",
      "ktau[1]",
      "krho[1]",
      "kwidth[1]",
      sprintf("kre_sdc[%d]", seq_len(J)),
      sprintf("kre_sda[%d]", seq_len(J)),
      "klag_sd[1]"
    )
  )
  global <- tibble::tibble(
    par = gv,
    mean = colMeans(dm[, gv, drop = FALSE]),
    q2.5 = apply(dm[, gv, drop = FALSE], 2, stats::quantile, 0.025),
    q97.5 = apply(dm[, gv, drop = FALSE], 2, stats::quantile, 0.975)
  )
  out <- list(classes = classes, global = global)
  cv <- grep("^kclag\\[", colnames(dm), value = TRUE)
  if (length(cv)) {
    out$country_lag <- kl_country_lag_summary(fit)
  } # kernel.lag_re
  out
}

#' kl_country_lag_summary
#'
#' Country lag offsets of a kernel.lag_re fit: per country the offset d_c
#' (kclag) and its own lag centre, membership-weighted per draw,
#' sum_k r_ck cent_kc, with the modal class. cent_kc = kmu_k + d_c under
#' kernel.lag_re.scale month (offset in lags), the drawn kccent[k, c] under
#' scale logit (offset in logit units of the lag window; the scale is told by
#' the draws: only the logit program has kccent).
#' @return tibble(iso3, class, scale, offset mean / q2.5 / q97.5, lag mean / q2.5 / q97.5)
kl_country_lag_summary <- function(fit) {
  d <- as.data.frame(fit$draws)
  K <- fit$K
  C <- fit$country_levels
  if (!sprintf("kclag[%d]", 1L) %in% names(d)) {
    return(NULL)
  }
  logit <- "kccent[1,1]" %in% names(d)
  q <- function(v) {
    c(mean(v), stats::quantile(v, c(0.025, 0.975), names = FALSE))
  }
  dplyr::bind_rows(lapply(seq_along(C), function(c) {
    off <- d[[sprintf("kclag[%d]", c)]]
    P <- sapply(seq_len(K), function(k) {
      exp(d[[sprintf("kl_log_resp[%d,%d]", k, c)]])
    })
    lag <- if (logit) {
      rowSums(
        P * sapply(seq_len(K), function(k) d[[sprintf("kccent[%d,%d]", k, c)]])
      )
    } else {
      rowSums(P * sapply(seq_len(K), function(k) d[[sprintf("kmu[%d]", k)]])) +
        off
    }
    a <- q(off)
    b <- q(lag)
    tibble::tibble(
      iso3 = C[c],
      class = fit$class_labels[which.max(colMeans(P))],
      scale = if (logit) "logit" else "month",
      offset_mean = a[1],
      offset_q2.5 = a[2],
      offset_q97.5 = a[3],
      lag_mean = b[1],
      lag_q2.5 = b[2],
      lag_q97.5 = b[3]
    )
  }))
}
