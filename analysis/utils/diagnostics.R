# diagnostics.R — MCMC convergence diagnostics of one saved fit (brmsfit or
#   brms_kl kl_fit).
#
# Provides:
#   fit_diagnostics()        R-hat, bulk / tail ESS, per-chain sampler table,
#                            settings, package versions, scalar summaries
#   fit_diagnostics_row()    one-row summary of a fit_diagnostics() result
#   diag_scope_summary()     summary restricted to a DIAG_SCOPES scope
#   fit_diagnostics_gate()   pure verdict on a scope summary (DIAG_GATE)
#
# Details:
#   Consumers: results_derive.R (09 convergence step), kl_sweep.R, and
#   deploy_core.R's convergence gate (deploy_prod_run.R --stage deploy_fits).
#   Callers source base.R; packages via `::`: posterior, brms, tibble.

#' Thresholds reported with every summary (the tests/test_diagnostics.R
#' convention): R-hat above `rhat` and bulk / tail ESS below `ess` are counted.
DIAG_THRESHOLDS <- c(rhat = 1.01, ess = 400)

#' Variables left out of the per-parameter summary: the log density and brms's
#' log prior are not model parameters.
DIAG_EXCLUDE <- "^(lp__|lprior)$"

#' Parameter scopes of the deploy convergence gate (fit_diagnostics_gate()):
#' an EXCLUSION regex over the per-parameter summary, NULL = every parameter.
#' `population` drops the group-level levels (`r_`), their standardised
#' deviates (`z_`) and the smooth deviates (`zs_`): thousands of levels sit at
#' R-hat 1.010-1.021 with bulk ESS above 2,000 (multiple-testing noise), while
#' the RE standard deviations (`sd_`, `sds_`), fixed effects and shapes stay.
DIAG_SCOPES <- list(all = NULL, population = "^(r_|z_|zs_)")

#' Bounds of the deploy gate that are not reported thresholds: R-hat on every
#' parameter (the simulator draws the `r_` levels too) and the minimum E-BFMI
#' per chain (Stan's convention).
DIAG_GATE <- c(rhat_any = 1.05, ebfmi = 0.3)

#' E-BFMI of one chain from its energy draws (Stan's definition).
#' @noRd
.ebfmi <- function(E) {
  E <- as.numeric(E)
  if (length(E) < 3L || !all(is.finite(E))) {
    return(NA_real_)
  }
  den <- sum((E - mean(E))^2)
  if (den <= 0) {
    return(NA_real_)
  }
  sum(diff(E)^2) / den
}

#' Per-parameter R-hat / ESS of a draws object, `exclude` variables dropped.
#' @noRd
.diag_summary <- function(draws, exclude = DIAG_EXCLUDE) {
  vars <- posterior::variables(draws)
  vars <- vars[!grepl(exclude, vars)]
  s <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = vars),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  tibble::tibble(
    variable = s$variable,
    rhat = s$rhat,
    ess_bulk = s$ess_bulk,
    ess_tail = s$ess_tail
  )
}

#' Sampler table of a brmsfit from brms::nuts_params() (post-warmup).
#' @noRd
.diag_sampler_brms <- function(fit, max_treedepth) {
  np <- brms::nuts_params(fit)
  chains <- sort(unique(np$Chain))
  rows <- lapply(chains, function(ch) {
    v <- function(p) np$Value[np$Chain == ch & np$Parameter == p]
    td <- v("treedepth__")
    tibble::tibble(
      chain = as.integer(ch),
      divergent = as.integer(sum(v("divergent__"))),
      treedepth_hits = if (is.na(max_treedepth)) {
        NA_integer_
      } else {
        as.integer(sum(td >= max_treedepth))
      },
      ebfmi = .ebfmi(v("energy__"))
    )
  })
  do.call(rbind, rows)
}

#' Sampler settings stored in a brmsfit (stan_args of the first chain).
#' @noRd
.diag_settings_brms <- function(fit) {
  sa <- tryCatch(fit$fit@stan_args[[1]], error = function(e) list())
  ctl <- sa$control %||% list()
  num <- function(x) {
    if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x)[1]
  }
  tibble::tibble(
    backend = as.character(fit$backend %||% NA_character_),
    n_chains = as.integer(brms::nchains(fit)),
    iter = num(sa$iter),
    warmup = num(sa$warmup),
    n_draws = as.integer(brms::ndraws(fit)),
    adapt_delta = num(ctl$adapt_delta),
    max_treedepth = num(ctl$max_treedepth)
  )
}

#' Package versions a brmsfit records (fit$version), one column per entry.
#' @noRd
.diag_versions <- function(v = NULL) {
  keys <- c("brms", "rstan", "stanHeaders", "cmdstanr", "cmdstan")
  vals <- vapply(
    keys,
    function(k) {
      x <- v[[k]]
      if (is.null(x)) NA_character_ else as.character(x)
    },
    character(1)
  )
  tibble::as_tibble(as.list(vals))
}

#' Scalar summaries of a summary + sampler table.
#' @noRd
.diag_scalars <- function(summary, sampler, thresholds) {
  mx <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  mn <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  sm <- function(x) {
    if (all(is.na(x))) NA_integer_ else as.integer(sum(x, na.rm = TRUE))
  }
  list(
    n_params = nrow(summary),
    max_rhat = mx(summary$rhat),
    min_ess_bulk = mn(summary$ess_bulk),
    min_ess_tail = mn(summary$ess_tail),
    n_rhat_gt = sum(summary$rhat > thresholds[["rhat"]], na.rm = TRUE),
    n_ess_bulk_lt = sum(summary$ess_bulk < thresholds[["ess"]], na.rm = TRUE),
    n_ess_tail_lt = sum(summary$ess_tail < thresholds[["ess"]], na.rm = TRUE),
    divergences = sm(sampler$divergent),
    treedepth_hits = sm(sampler$treedepth_hits),
    min_ebfmi = mn(sampler$ebfmi)
  )
}

#' fit_diagnostics
#'
#' @param fit a brmsfit or a kl_fit.
#' @param chains kl_fit only: chain ids to keep (e.g. the chains of the
#'   dominant mode); NULL = all chains. The sampler table keeps the same chains.
#' @param exclude regex of variables left out of the per-parameter summary.
#' @param thresholds named c(rhat, ess), see DIAG_THRESHOLDS.
#' @return list(summary = tibble(variable, rhat, ess_bulk, ess_tail),
#'   sampler = tibble(chain, divergent, treedepth_hits, ebfmi),
#'   settings = tibble(backend, n_chains, iter, warmup, n_draws, adapt_delta,
#'   max_treedepth), versions = tibble(brms, rstan, stanHeaders, cmdstanr,
#'   cmdstan) (NA for a kl_fit), thresholds, and the scalars n_params,
#'   max_rhat, min_ess_bulk, min_ess_tail, n_rhat_gt, n_ess_bulk_lt,
#'   n_ess_tail_lt, divergences, treedepth_hits, min_ebfmi.
fit_diagnostics <- function(
  fit,
  chains = NULL,
  exclude = DIAG_EXCLUDE,
  thresholds = DIAG_THRESHOLDS
) {
  if (inherits(fit, "brmsfit")) {
    if (!is.null(chains)) {
      stop("fit_diagnostics: `chains` is supported for kl_fit only.")
    }
    settings <- .diag_settings_brms(fit)
    summary <- .diag_summary(posterior::as_draws_array(fit), exclude)
    sampler <- .diag_sampler_brms(fit, settings$max_treedepth)
    versions <- .diag_versions(fit$version)
  } else if (inherits(fit, "kl_fit")) {
    draws <- fit$draws
    all_ch <- sort(unique(draws$.chain))
    keep <- if (is.null(chains)) all_ch else as.integer(chains)
    if (!all(keep %in% all_ch)) {
      stop(
        "fit_diagnostics: chains ",
        paste(setdiff(keep, all_ch), collapse = ","),
        " not in the fit."
      )
    }
    if (!is.null(chains)) {
      draws <- posterior::subset_draws(draws, chain = keep)
    }
    summary <- .diag_summary(draws, exclude)
    dg <- fit$diagnostics
    pick <- function(x) {
      if (is.null(x) || length(x) < max(keep)) {
        rep(NA_real_, length(keep))
      } else {
        as.numeric(x[keep])
      }
    }
    sampler <- tibble::tibble(
      chain = as.integer(keep),
      divergent = as.integer(pick(dg$num_divergent)),
      treedepth_hits = as.integer(pick(dg$num_max_treedepth)),
      ebfmi = pick(dg$ebfmi)
    )
    settings <- tibble::tibble(
      backend = "brms_kl",
      n_chains = length(keep),
      iter = NA_real_,
      warmup = NA_real_,
      n_draws = as.integer(posterior::ndraws(draws)),
      adapt_delta = NA_real_,
      max_treedepth = NA_real_
    )
    versions <- .diag_versions(NULL)
  } else {
    stop(
      "fit_diagnostics: unsupported fit class ",
      paste(class(fit), collapse = "/"),
      " (brmsfit or kl_fit)."
    )
  }
  c(
    list(
      summary = summary,
      sampler = sampler,
      settings = settings,
      versions = versions,
      thresholds = thresholds
    ),
    .diag_scalars(summary, sampler, thresholds)
  )
}

#' One-row tibble of the scalar fields of a fit_diagnostics() result.
fit_diagnostics_row <- function(d) {
  tibble::tibble(
    n_chains = d$settings$n_chains,
    n_draws = d$settings$n_draws,
    n_params = d$n_params,
    max_rhat = d$max_rhat,
    min_ess_bulk = d$min_ess_bulk,
    min_ess_tail = d$min_ess_tail,
    n_rhat_gt = d$n_rhat_gt,
    n_ess_bulk_lt = d$n_ess_bulk_lt,
    n_ess_tail_lt = d$n_ess_tail_lt,
    divergences = d$divergences,
    treedepth_hits = d$treedepth_hits,
    min_ebfmi = d$min_ebfmi
  )
}

#' The per-parameter summary restricted to a gate scope (DIAG_SCOPES).
#' @param summary tibble(variable, rhat, ess_bulk, ess_tail) (fit_diagnostics()$summary)
#' @param scope   a name of DIAG_SCOPES
diag_scope_summary <- function(summary, scope = "population") {
  if (
    !(is.character(scope) &&
      length(scope) == 1L &&
      scope %in% names(DIAG_SCOPES))
  ) {
    stop(
      "diag_scope_summary: scope must be one of {",
      paste(names(DIAG_SCOPES), collapse = ", "),
      "} (got '",
      paste(scope, collapse = ","),
      "').",
      call. = FALSE
    )
  }
  ex <- DIAG_SCOPES[[scope]]
  if (is.null(ex)) {
    return(summary)
  }
  summary[!grepl(ex, summary$variable), , drop = FALSE]
}

#' fit_diagnostics_gate
#'
#' Pure pass / fail verdict of one fit_diagnostics() result (the deploy
#' convergence gate). The scoped parameters (DIAG_SCOPES) are held to R-hat <=
#' `rhat_max` and bulk and tail ESS >= `ess_min`; every parameter to R-hat <=
#' `rhat_any`; the summed divergences to <= `max_divergent`; every chain's
#' E-BFMI to >= `ebfmi_min`. Max-treedepth hits are reported, never gated.
#' Fails closed: an empty or all-NA scope, or an NA divergence / E-BFMI
#' summary, fails with reason `na`.
#'
#' @param d a fit_diagnostics() result (needs `summary` and `sampler`).
#' @return one-row tibble(scope, n_params_scope, max_rhat_scope,
#'   max_rhat_variable, min_ess_bulk_scope, min_ess_tail_scope, max_rhat_all,
#'   divergences, treedepth_hits, min_ebfmi, pass, reasons); `reasons` is ""
#'   on a pass, else ";"-joined from rhat, ess_bulk, ess_tail, rhat_any,
#'   divergences, ebfmi, na.
fit_diagnostics_gate <- function(
  d,
  scope = "population",
  rhat_max = DIAG_THRESHOLDS[["rhat"]],
  ess_min = DIAG_THRESHOLDS[["ess"]],
  max_divergent = 0L,
  rhat_any = DIAG_GATE[["rhat_any"]],
  ebfmi_min = DIAG_GATE[["ebfmi"]]
) {
  ss <- diag_scope_summary(d$summary, scope)
  th <- c(rhat = rhat_max, ess = ess_min)
  sc <- .diag_scalars(ss, d$sampler, th)
  sa <- .diag_scalars(d$summary, d$sampler, th)
  worst <- if (nrow(ss) && any(!is.na(ss$rhat))) {
    ss$variable[which.max(ss$rhat)]
  } else {
    NA_character_
  }
  is_na <- sc$n_params == 0L ||
    is.na(sc$max_rhat) ||
    is.na(sc$min_ess_bulk) ||
    is.na(sc$min_ess_tail) ||
    is.na(sa$max_rhat) ||
    is.na(sc$divergences) ||
    is.na(sc$min_ebfmi)
  gt <- function(x, b) isTRUE(x > b)
  lt <- function(x, b) isTRUE(x < b)
  reasons <- c(
    if (gt(sc$max_rhat, rhat_max)) "rhat",
    if (lt(sc$min_ess_bulk, ess_min)) "ess_bulk",
    if (lt(sc$min_ess_tail, ess_min)) "ess_tail",
    if (gt(sa$max_rhat, rhat_any)) "rhat_any",
    if (gt(sc$divergences, max_divergent)) "divergences",
    if (lt(sc$min_ebfmi, ebfmi_min)) "ebfmi",
    if (is_na) "na"
  )
  tibble::tibble(
    scope = scope,
    n_params_scope = sc$n_params,
    max_rhat_scope = sc$max_rhat,
    max_rhat_variable = worst,
    min_ess_bulk_scope = sc$min_ess_bulk,
    min_ess_tail_scope = sc$min_ess_tail,
    max_rhat_all = sa$max_rhat,
    divergences = sc$divergences,
    treedepth_hits = sc$treedepth_hits,
    min_ebfmi = sc$min_ebfmi,
    pass = !length(reasons),
    reasons = paste(reasons, collapse = ";")
  )
}
