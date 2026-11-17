# kl_sweep.R — comparison of a brms_kl latent-class sweep, behind
#   05_fit_postprocess_kl_sweep.R and 06_fit_figures_kl_sweep.R.
#
# Provides:
#   kl_sweep_inventory()            fit inventory of the sweep config dirs
#   kl_chain_modes(), kl_agreeing_chains()   chain modes, agreeing chains
#   kl_fit_row(), kl_fit_diagnostics_tables()   R-hat / ESS, all + agreeing
#   kl_class_rows(), kl_partition_rows()     class and MAP-partition rows
#   kl_responsibility_rows(), kl_country_mixture_rows(),
#   kl_partition_stability()        membership uncertainty
#   adjusted_rand_index(), kl_match_classes()   label-free helpers
#   kl_sweep_loo(), kl_sweep_compare()   country-block LOO / WAIC, stacking
#
# Details:
#   The country mixture is the M1 sum_k r_ck * kcon_k / kmu_k vs the hard-MAP
#   class. Callers source base.R, config.R, paths.R, dlnm_kernel.R and
#   diagnostics.R; packages via `::`: dplyr, tibble, posterior, loo.

#' Variables left out of the per-fit R-hat / ESS summaries: per-row and
#' per-country log likelihoods and logits, the membership field, raw RE
#' deviations, kernel weights and the log prior (lp__ is kept).
KL_SWEEP_EXCL <- "^kl_log|^klogit|^kfield|^kre_c|^kre_z|^lprior|^kw|^ktheta"

#' Strip the "enso_" prefix of a class label.
#' @noRd
.kl_short <- function(x) sub("^enso_", "", x)

#' Basis token of a dlnm argvar spec: fun + df ("ns5", "bs5"). A missing fun
#' falls back to the runtime default (get_default_options()'s dlnm argvar fun,
#' "bs"); resolved configs always carry it.
#' @noRd
.kl_basis_token <- function(argvar) {
  paste0(argvar[["fun"]] %||% "bs", argvar[["df"]] %||% "")
}

#' kl_sweep_inventory
#'
#' One row per config of the given sweep dirs (manifest.yml excluded), with
#' the expected fit path. K = kernel.classes under identify repel, else the
#' number of seeded classes (the prior registry's group count without
#' seeds), basis from the ENSO argvar, mix from the
#' contiguity prior -- read from the resolved config, never from the filename.
#' @param batch_dirs character vector of config dirs.
#' @return tibble(config, model, basis, K, mix, id, fit_file, present)
kl_sweep_inventory <- function(batch_dirs) {
  rows <- lapply(batch_dirs, function(d) {
    cfgs <- sort(list.files(d, "\\.yml$", full.names = TRUE))
    cfgs <- cfgs[basename(cfgs) != "manifest.yml"]
    lapply(cfgs, function(cfg) {
      opt <- NULL
      invisible(utils::capture.output(
        opt <- suppressMessages(make_options_from_config(
          option_list = list(),
          cli_argv = c("--config", cfg),
          model = "occurrence"
        ))
      ))
      if (!identical(opt$occurrence$backend, "brms_kl")) {
        stop(
          "kl_sweep_inventory: ",
          cfg,
          " is not a brms_kl config (backend '",
          opt$occurrence$backend,
          "')."
        )
      }
      kc <- opt$occurrence$kernel %||% list()
      dl <- opt$occurrence$dlnm
      av <- (dl$variables$enso %||% list())$argvar %||%
        dl$defaults$argvar %||%
        list()
      # K = classes the fit will have: kernel.classes under identify repel
      # (exchangeable classes; the registry only pulls), one per seeded group
      # under seed identification (registry groups without a fitted country,
      # e.g. enso_central, get none), else the registry's group count.
      src <- kc$prior_grouping$source %||% ""
      seeds <- (kc$seeds %||% list())$iso3 %||% list()
      K <- if (identical(kc[["identify"]], "repel")) {
        as.integer(kc[["classes"]])
      } else if (length(seeds)) {
        length(seeds)
      } else if (nzchar(src)) {
        length(unique(load_grouping(src)$group))
      } else {
        as.integer(kc$classes %||% NA)
      }
      tibble::tibble(
        config = cfg,
        model = sub("\\.yml$", "", basename(cfg)),
        basis = .kl_basis_token(av),
        K = as.integer(K),
        mix = isTRUE(kc$prior_grouping$contiguity$enable),
        fit_file = make_run_file(opt, "fit")
      )
    })
  })
  inv <- dplyr::bind_rows(unlist(rows, recursive = FALSE))
  if (!nrow(inv)) {
    stop(
      "kl_sweep_inventory: no configs in ",
      paste(batch_dirs, collapse = ", ")
    )
  }
  inv$present <- file.exists(inv$fit_file)
  inv$id <- sprintf("%s_K%d%s", inv$basis, inv$K, ifelse(inv$mix, "_mix", ""))
  if (anyDuplicated(inv$id)) {
    stop(
      "kl_sweep_inventory: duplicate fit ids (",
      paste(unique(inv$id[duplicated(inv$id)]), collapse = ", "),
      "): two configs share basis, K and prior."
    )
  }
  inv
}

#' kl_chain_modes
#'
#' Chain modes of a kl_fit: chains whose per-chain mean class centres kmu all
#' lie within `mode_tol` lags of each other (complete linkage on the maximum
#' distance). Modes are ordered by size, then by mean lp__.
#' @return list(chains, centres (chains x K), lp (per chain), mode (per chain,
#'   1 = dominant), sizes)
kl_chain_modes <- function(fit, mode_tol = 1.5) {
  dr <- fit$draws
  ch <- dr$.chain
  chains <- sort(unique(ch))
  K <- fit$K
  cm <- vapply(
    seq_len(K),
    function(k) {
      as.numeric(tapply(
        posterior::extract_variable(dr, sprintf("kmu[%d]", k)),
        ch,
        mean
      ))
    },
    numeric(length(chains))
  )
  cm <- matrix(cm, nrow = length(chains), dimnames = list(chains, NULL))
  lpc <- tapply(posterior::extract_variable(dr, "lp__"), ch, mean)
  grp <- if (length(chains) > 1L) {
    stats::cutree(
      stats::hclust(stats::dist(cm, "maximum"), "complete"),
      h = mode_tol
    )
  } else {
    stats::setNames(1L, chains)
  }
  sz <- table(grp)
  lpm <- tapply(lpc[names(grp)], grp, mean)
  ord <- names(sz)[order(-as.integer(sz), -lpm[names(sz)])]
  mode <- match(as.character(grp), ord)
  list(
    chains = as.integer(chains),
    centres = cm,
    lp = stats::setNames(as.numeric(lpc), chains),
    mode = stats::setNames(mode, chains),
    sizes = as.integer(sz[ord])
  )
}

#' kl_agreeing_chains
#'
#' The chains a sweep fit is summarised on: the dominant mode
#' (kl_chain_modes()), or with n_agree > 0 the n chains whose per-chain mean
#' centres agree best (smallest maximum pairwise |kmu| difference; ties:
#' smaller maximum co-membership difference, then higher mean lp__), whether
#' or not they form one mode.
#' @return list(chains = integer ids, kmu_spread = max |kmu| difference among them)
kl_agreeing_chains <- function(
  fit,
  mode_tol = 1.5,
  n_agree = 0L,
  modes = kl_chain_modes(fit, mode_tol)
) {
  chains <- modes$chains
  cm <- modes$centres
  D2 <- as.matrix(stats::dist(cm, "maximum"))
  if (n_agree > 0L) {
    if (n_agree > length(chains)) {
      stop(
        "kl_agreeing_chains: n_agree ",
        n_agree,
        " > ",
        length(chains),
        " chains."
      )
    }
    M <- kl_comembership(fit, by_chain = TRUE)
    D1 <- outer(
      seq_along(chains),
      seq_along(chains),
      Vectorize(function(a, b) mean(abs(M[[a]] - M[[b]])))
    )
    cs <- utils::combn(length(chains), n_agree)
    sc <- apply(cs, 2, function(s) {
      c(max(D2[s, s]), max(D1[s, s]), -mean(modes$lp[s]))
    })
    sel <- cs[, order(sc[1, ], sc[2, ], sc[3, ])[1]]
  } else {
    sel <- which(modes$mode == 1L)
  }
  list(
    chains = chains[sel],
    kmu_spread = if (length(sel) > 1L) max(D2[sel, sel]) else 0
  )
}

#' kl_fit_row
#'
#' One summary row of a sweep fit: modes, agreeing chains, sampler
#' diagnostics, R-hat / ESS on all chains and on the agreeing chains.
#' @param r one inventory row; @param modes kl_chain_modes(); @param ag
#'   kl_agreeing_chains().
kl_fit_row <- function(fit, r, modes, ag, exclude = KL_SWEEP_EXCL) {
  da <- fit_diagnostics(fit, exclude = exclude)
  dg <- if (length(ag$chains) >= 2L) {
    fit_diagnostics(fit, chains = ag$chains, exclude = exclude)
  } else {
    NULL
  }
  g <- function(d, f) if (is.null(d)) NA_real_ else d[[f]]
  tibble::tibble(
    id = r$id,
    basis = r$basis,
    K = fit$K,
    mix = r$mix,
    config = r$model,
    chains = length(modes$chains),
    n_modes = length(modes$sizes),
    mode_sizes = paste(modes$sizes, collapse = "+"),
    agreeing = paste(ag$chains, collapse = ","),
    kmu_spread_agree = ag$kmu_spread,
    n_params = da$n_params,
    divergences = da$divergences,
    treedepth_hits = da$treedepth_hits,
    min_ebfmi = da$min_ebfmi,
    max_rhat_all = da$max_rhat,
    n_rhat105_all = sum(da$summary$rhat > 1.05, na.rm = TRUE),
    min_ess_bulk_all = da$min_ess_bulk,
    min_ess_tail_all = da$min_ess_tail,
    max_rhat_agree = g(dg, "max_rhat"),
    n_rhat105_agree = if (is.null(dg)) {
      NA_integer_
    } else {
      sum(dg$summary$rhat > 1.05, na.rm = TRUE)
    },
    min_ess_bulk_agree = g(dg, "min_ess_bulk"),
    min_ess_tail_agree = g(dg, "min_ess_tail"),
    lp_agree = mean(modes$lp[as.character(ag$chains)]),
    lp_spread = diff(range(modes$lp))
  )
}

#' Per-variable and per-chain diagnostics tables of one fit (all chains and
#' the agreeing chains), for diagnostics.csv / diagnostics_chains.csv.
#' @return list(variables = tibble(id, chains, variable, rhat, ess_bulk,
#'   ess_tail), chains = tibble(id, chain, agreeing, divergent,
#'   treedepth_hits, ebfmi))
kl_fit_diagnostics_tables <- function(fit, id, ag, exclude = KL_SWEEP_EXCL) {
  da <- fit_diagnostics(fit, exclude = exclude)
  v <- dplyr::mutate(da$summary, id = id, chains = "all", .before = 1)
  if (length(ag$chains) >= 2L) {
    v <- dplyr::bind_rows(
      v,
      dplyr::mutate(
        fit_diagnostics(fit, chains = ag$chains, exclude = exclude)$summary,
        id = id,
        chains = "agree",
        .before = 1
      )
    )
  }
  list(
    variables = v,
    chains = dplyr::mutate(
      da$sampler,
      id = id,
      agreeing = chain %in% ag$chains,
      .before = 1
    )
  )
}

#' Class rows (agreeing chains): kmu / kcon / kpi per class and kkappa, with
#' mean, 95% interval and R-hat.
kl_class_rows <- function(fit_sub, r) {
  ds <- fit_sub$draws
  K <- fit_sub$K
  key <- intersect(
    c(
      sprintf("kmu[%d]", 1:K),
      sprintf("kcon[%d]", 1:K),
      sprintf("kpi[%d]", 1:K),
      "kkappa"
    ),
    posterior::variables(ds)
  )
  s <- posterior::summarise_draws(
    posterior::subset_draws(ds, variable = key),
    "mean",
    ~ stats::quantile(.x, c(.025, .975)),
    "rhat"
  )
  idx <- suppressWarnings(as.integer(sub(".*\\[(\\d+)\\]", "\\1", s$variable)))
  tibble::tibble(
    id = r$id,
    basis = r$basis,
    K = K,
    mix = r$mix,
    param = sub("\\[.*", "", s$variable),
    class = ifelse(
      grepl("\\[", s$variable),
      .kl_short(fit_sub$class_labels)[idx],
      ""
    ),
    mean = s$mean,
    q2.5 = s$`2.5%`,
    q97.5 = s$`97.5%`,
    rhat = s$rhat
  )
}

#' MAP partition rows (agreeing chains): each country's argmax class and its
#' posterior probability.
kl_partition_rows <- function(fit_sub, r) {
  resp <- kl_responsibilities(fit_sub)
  pmax <- dplyr::summarise(
    dplyr::group_by(resp, country),
    p = max(mean),
    .groups = "drop"
  )
  kl_map_partition(fit_sub, "argmax", resp) |>
    dplyr::left_join(pmax, by = c(iso3 = "country")) |>
    dplyr::mutate(id = r$id, group = .kl_short(group), .before = 1)
}

# Membership uncertainty ------------------------------------------------------
#
# Class index -> label: every row below names class k by fit$class_labels[k],
# the fit's own names. Under identify: seed these are the prior registry
# groups the seed countries fix (dlnm_kernel.R::fit_kl_stan(), labels =
# kl_name_classes(); no per-draw relabelling), so averaging over chains that
# share a mode is consistent and the labels are the ones
# kl_grouping.R::kl_grouping_from_fit() writes into a derived registry. Under
# identify repel the draws were relabelled per draw when the fit was made
# (dlnm_kernel.R::kl_repel_relabel()) and the labels are registry groups by
# overlap (klc<k> for classes left over). The fit's `identify` is recorded in
# every row; nothing is relabelled here.

#' Draws x K x C responsibility array of a kl_fit (exp of kl_log_resp).
#' @noRd
.kl_resp_array <- function(fit) {
  exp(.kl_draw_array(fit, "kl_log_resp", c(fit$K, length(fit$country_levels))))
}

#' MAP class index per country: the class of largest mean responsibility
#' (first on ties, as kl_map_partition(method = "argmax")).
#' @noRd
.kl_map_class <- function(R) apply(apply(R, c(2, 3), mean), 2, which.max)

#' The prior registry of a sweep config (kernel.prior_grouping.source), or
#' NULL when it has none.
#' @return tibble(iso3, group) or NULL
kl_sweep_prior_grouping <- function(config) {
  opt <- NULL
  invisible(utils::capture.output(
    opt <- suppressMessages(make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", config),
      model = "occurrence"
    ))
  ))
  src <- (opt$occurrence$kernel %||% list())$prior_grouping$source %||% ""
  if (!nzchar(src)) {
    return(NULL)
  }
  g <- load_grouping(src)
  tibble::tibble(
    iso3 = as.character(g$iso3),
    group = as.character(g$group),
    source = src
  )
}

#' Each country's prior registry group (NA when the country or the registry
#' is absent).
#' @noRd
.kl_prior_of <- function(iso3, prior) {
  if (is.null(prior)) {
    return(rep(NA_character_, length(iso3)))
  }
  as.character(prior$group[match(iso3, prior$iso3)])
}

#' kl_responsibility_rows
#'
#' The full country x class responsibility matrix of one fit (posterior mean,
#' sd and 95% interval of r_ck = exp(kl_log_resp[k, c])), on the chains of
#' `fit_sub` (tagged `chains`, "all" or "agree"), with each country's prior
#' registry group and the fit's identification rule. Rows of one country sum
#' to 1 in `mean`.
#' @param fit_sub kl_fit (or a kl_subset_chains() of it); @param r one
#'   inventory row (id); @param prior tibble(iso3, group) or NULL.
#' @return tibble(id, chains, iso3, class, label, mean, sd, q2.5, q97.5,
#'   prior_group, identify)
kl_responsibility_rows <- function(
  fit_sub,
  r,
  chains = c("agree", "all"),
  prior = NULL
) {
  chains <- match.arg(chains)
  resp <- kl_responsibilities(fit_sub)
  tibble::tibble(
    id = r$id,
    chains = chains,
    iso3 = resp$country,
    class = resp$class,
    label = resp$label,
    mean = resp$mean,
    sd = resp$sd,
    q2.5 = resp$q2.5,
    q97.5 = resp$q97.5,
    prior_group = .kl_prior_of(resp$country, prior),
    identify = fit_sub$identify %||% NA_character_
  )
}

#' kl_country_mixture_rows
#'
#' Membership-propagated class effects per country (M1), on the chains of
#' `fit_sub`: per draw the mixture sum_k r_ck * x_k and the hard-MAP value
#' x_MAP(c), for x = kcon (El Nino contrast: log RR of a sustained ENSO value
#' erf_at vs cen, dlnm_kernel.R tparams) and x = kmu (class lag centre; the
#' Gaussian form, and the discrete form's drawn centre exposed as kmu). Both
#' summarised with their difference (mixture - MAP). A quantity the fit's
#' form does not have (kmu under form dlnm) gets NA with the reason in `note`.
#' The class-independent country deviation of kernel.re_slope adds the same
#' term to both and is not included (it cancels in the difference).
#' @return tibble(id, iso3, map_class, map_label, p_max, prior_group,
#'   quantity, mix_mean, mix_sd, mix_q2.5, mix_q97.5, map_mean, map_q2.5,
#'   map_q97.5, diff_mean, diff_q2.5, diff_q97.5, identify, note)
kl_country_mixture_rows <- function(fit_sub, r, prior = NULL) {
  R <- .kl_resp_array(fit_sub)
  K <- fit_sub$K
  ctry <- fit_sub$country_levels
  map <- .kl_map_class(R)
  pmax <- apply(apply(R, c(2, 3), mean), 2, max)
  vars <- posterior::variables(fit_sub$draws)
  q <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE)
  one <- function(par) {
    have <- all(sprintf("%s[%d]", par, seq_len(K)) %in% vars)
    base <- tibble::tibble(
      id = r$id,
      iso3 = ctry,
      map_class = map,
      map_label = fit_sub$class_labels[map],
      p_max = pmax,
      prior_group = .kl_prior_of(ctry, prior),
      quantity = par
    )
    if (!have) {
      na <- rep(NA_real_, length(ctry))
      return(dplyr::mutate(
        base,
        mix_mean = na,
        mix_sd = na,
        mix_q2.5 = na,
        mix_q97.5 = na,
        map_mean = na,
        map_q2.5 = na,
        map_q97.5 = na,
        diff_mean = na,
        diff_q2.5 = na,
        diff_q97.5 = na,
        identify = fit_sub$identify %||% NA_character_,
        note = sprintf(
          "no %s in this fit (form %s)",
          par,
          fit_sub$form %||% "unknown"
        )
      ))
    }
    X <- .kl_draw_array(fit_sub, par, K) # draws x K
    nd <- nrow(X)
    mix <- matrix(
      vapply(
        seq_along(ctry),
        function(c) rowSums(matrix(R[,, c], nrow = nd) * X),
        numeric(nd)
      ),
      nrow = nd
    )
    hard <- X[, map, drop = FALSE]
    dif <- mix - hard
    s <- function(M) t(apply(M, 2, function(x) c(mean(x), stats::sd(x), q(x))))
    sm <- s(mix)
    sh <- s(hard)
    sd_ <- s(dif)
    dplyr::mutate(
      base,
      mix_mean = sm[, 1],
      mix_sd = sm[, 2],
      mix_q2.5 = sm[, 3],
      mix_q97.5 = sm[, 4],
      map_mean = sh[, 1],
      map_q2.5 = sh[, 3],
      map_q97.5 = sh[, 4],
      diff_mean = sd_[, 1],
      diff_q2.5 = sd_[, 3],
      diff_q97.5 = sd_[, 4],
      identify = fit_sub$identify %||% NA_character_,
      note = ""
    )
  }
  dplyr::bind_rows(one("kcon"), one("kmu"))
}

#' adjusted_rand_index
#'
#' Hubert-Arabie adjusted Rand index of two partitions of the same items
#' (label-free). 1 when both are the same trivial partition (all in one
#' cluster, or all singletons), where the index is 0 / 0.
#' @param a,b vectors of cluster labels, same length.
#' @return numeric(1)
adjusted_rand_index <- function(a, b) {
  if (length(a) != length(b)) {
    stop("adjusted_rand_index: partitions differ in length.")
  }
  if (anyNA(a) || anyNA(b)) {
    stop("adjusted_rand_index: NA labels.")
  }
  n <- length(a)
  if (n < 2L) {
    return(1)
  }
  tab <- table(as.character(a), as.character(b))
  c2 <- function(x) sum(x * (x - 1) / 2)
  sij <- c2(tab)
  sa <- c2(rowSums(tab))
  sb <- c2(colSums(tab))
  ex <- sa * sb / (n * (n - 1) / 2)
  mx <- (sa + sb) / 2
  if (mx - ex == 0) {
    return(1)
  } # only when both are the same trivial partition
  (sij - ex) / (mx - ex)
}

#' All permutations of 1..n (rows), n <= 7.
#' @noRd
.kl_perms <- function(n) {
  if (n <= 1L) {
    return(matrix(seq_len(n), nrow = 1L))
  }
  p <- .kl_perms(n - 1L)
  do.call(
    rbind,
    lapply(seq_len(n), function(i) cbind(i, ifelse(p >= i, p + 1L, p)))
  )
}

#' kl_match_classes
#'
#' Optimal one-to-one matching of the class labels of partition `x` to those
#' of partition `y` (same items), maximising the number of items whose matched
#' labels agree; exhaustive over permutations (at most 7 labels on the larger
#' side). Labels left without a partner (unequal class counts) map to NA.
#' @return list(map = named character (x label -> y label), n_agree,
#'   agree = logical per item, n = items)
kl_match_classes <- function(x, y) {
  if (length(x) != length(y)) {
    stop("kl_match_classes: partitions differ in length.")
  }
  x <- as.character(x)
  y <- as.character(y)
  lx <- sort(unique(x))
  ly <- sort(unique(y))
  m <- max(length(lx), length(ly))
  if (m > 7L) {
    stop(
      "kl_match_classes: more than 7 classes, exhaustive matching not supported."
    )
  }
  tab <- matrix(0, m, m)
  tab[seq_along(lx), seq_along(ly)] <- table(factor(x, lx), factor(y, ly))
  P <- .kl_perms(m)
  sc <- apply(P, 1, function(p) sum(tab[cbind(seq_len(m), p)]))
  best <- P[which.max(sc), ]
  to <- ly[best[seq_along(lx)]] # NA where the partner is padding
  map <- stats::setNames(to, lx)
  agree <- !is.na(map[x]) & map[x] == y
  list(map = map, n_agree = sum(agree), agree = unname(agree), n = length(x))
}

#' kl_partition_stability
#'
#' Partition stability under membership uncertainty, on the chains of
#' `fit_sub`: per draw, one partition sampled country by country from
#' Cat(r_c) (that draw's responsibilities), compared with the MAP partition
#' (largest mean responsibility) by the adjusted Rand index and by the number
#' of countries in another class. Fixed `seed`; the caller's RNG state is
#' restored.
#' @return tibble(id, chains, n_draws, n_countries, ari_mean, ari_q2.5,
#'   ari_q25, ari_q50, ari_q75, ari_q97.5, p_ari1, moved_mean, moved_q97.5,
#'   seed)
kl_partition_stability <- function(
  fit_sub,
  r,
  chains = c("agree", "all"),
  seed = 1L
) {
  chains <- match.arg(chains)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }
  set.seed(seed)
  R <- .kl_resp_array(fit_sub)
  nd <- dim(R)[1]
  K <- dim(R)[2]
  C <- dim(R)[3]
  map <- .kl_map_class(R)
  U <- matrix(stats::runif(nd * C), nd, C)
  ari <- numeric(nd)
  moved <- integer(nd)
  for (d in seq_len(nd)) {
    cp <- apply(matrix(R[d, , ], nrow = K), 2, cumsum) # K x C
    cp <- matrix(cp, nrow = K)
    cp[K, ] <- 1 # guard rounding
    z <- colSums(sweep(cp, 2, U[d, ], "<")) + 1L
    ari[d] <- adjusted_rand_index(z, map)
    moved[d] <- sum(z != map)
  }
  qa <- stats::quantile(ari, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  tibble::tibble(
    id = r$id,
    chains = chains,
    n_draws = nd,
    n_countries = C,
    ari_mean = mean(ari),
    ari_q2.5 = qa[1],
    ari_q25 = qa[2],
    ari_q50 = qa[3],
    ari_q75 = qa[4],
    ari_q97.5 = qa[5],
    p_ari1 = mean(ari >= 1 - 1e-12),
    moved_mean = mean(moved),
    moved_q97.5 = stats::quantile(moved, 0.975, names = FALSE),
    seed = as.integer(seed)
  )
}

#' Stop unless all fits share countries and rows (LOO comparable).
kl_sweep_check_comparable <- function(fits) {
  lv <- lapply(fits, `[[`, "country_levels")
  nn <- vapply(fits, function(f) as.numeric(f$data_dims$N), numeric(1))
  if (length(unique(lv)) != 1L || length(unique(nn)) != 1L) {
    stop(
      "kl_sweep: fits differ in countries or rows, LOO not comparable (N = ",
      paste(unique(nn), collapse = ", "),
      ")."
    )
  }
  invisible(lv[[1]])
}

#' Draws x countries log-likelihood matrix of a kl_fit.
#' @noRd
.kl_ll_matrix <- function(f) {
  a <- .kl_draw_array(f, "kl_log_lik_ctry", length(f$country_levels))
  matrix(a, nrow = dim(a)[1])
}

#' kl_sweep_loo
#'
#' Country-block PSIS-LOO (kl_loo_country(), descriptive: high Pareto k is
#' expected) and WAIC of every fit, on all chains and on the agreeing chains;
#' loo_compare differences and stacking weights (loo::stacking_weights() on
#' the pointwise elpd, as kl_pin_combine()).
#' @param fits named list of kl_fit; @param agree named list of chain ids.
#' @return list(res = list(all, agree) of list(L, W, ll), table, country)
kl_sweep_loo <- function(fits, agree) {
  ctry <- kl_sweep_check_comparable(fits)
  one <- function(sel) {
    fs <- lapply(names(fits), function(n) {
      if (sel) kl_subset_chains(fits[[n]], agree[[n]]) else fits[[n]]
    })
    names(fs) <- names(fits)
    list(
      L = lapply(fs, function(f) {
        suppressWarnings(suppressMessages(kl_loo_country(f)))
      }),
      W = lapply(fs, function(f) suppressWarnings(loo::waic(.kl_ll_matrix(f)))),
      ll = vapply(fs, function(f) mean(rowSums(.kl_ll_matrix(f))), numeric(1))
    )
  }
  res <- list(all = one(FALSE), agree = one(TRUE))
  # loo_compare() rows: model names in rownames (older loo) or a `model` column (newer)
  cmp <- function(x, nm) {
    d <- as.data.frame(x)
    key <- if ("model" %in% names(d)) as.character(d$model) else rownames(d)
    d[match(nm, key), , drop = FALSE]
  }
  tab <- dplyr::bind_rows(lapply(names(res), function(s) {
    L <- res[[s]]$L
    W <- res[[s]]$W
    nm <- names(L)
    if (length(L) >= 2L) {
      cl <- cmp(loo::loo_compare(L), nm)
      cw <- cmp(loo::loo_compare(W), nm)
      pw <- do.call(cbind, lapply(L, function(l) l$pointwise[, "elpd_loo"]))
      wt <- as.numeric(loo::stacking_weights(pw))
    } else {
      cl <- cw <- data.frame(elpd_diff = 0, se_diff = 0)
      wt <- 1
    }
    tibble::tibble(
      id = nm,
      chains = s,
      elpd_loo = vapply(L, function(l) l$estimates["elpd_loo", 1], 1),
      se_elpd_loo = vapply(L, function(l) l$estimates["elpd_loo", 2], 1),
      p_loo = vapply(L, function(l) l$estimates["p_loo", 1], 1),
      n_k07 = vapply(L, function(l) sum(l$diagnostics$pareto_k > 0.7), 1),
      elpd_diff_loo = cl$elpd_diff,
      se_diff_loo = cl$se_diff,
      elpd_diff_waic = cw$elpd_diff,
      se_diff_waic = cw$se_diff,
      stacking_weight = wt,
      loglik_insample = res[[s]]$ll[nm]
    )
  }))
  country <- dplyr::bind_rows(lapply(names(res), function(s) {
    dplyr::bind_rows(lapply(names(res[[s]]$L), function(n) {
      tibble::tibble(
        chains = s,
        id = n,
        iso3 = ctry,
        elpd = res[[s]]$L[[n]]$pointwise[, "elpd_loo"],
        pareto_k = res[[s]]$L[[n]]$diagnostics$pareto_k
      )
    }))
  }))
  list(res = res, table = tab, country = country)
}

#' kl_sweep_compare
#'
#' The whole comparison of one sweep: inventory -> per-fit modes, rows,
#' diagnostics, classes, partition and membership uncertainty
#' (responsibilities on all and agreeing chains, the M1 country mixture and
#' partition stability on the agreeing chains) -> cross-fit LOO. Fits absent
#' on disk are listed and skipped.
#' @param stability_seed seed of kl_partition_stability() (per fit).
#' @return list(inv, agree, fits (table), classes, partition, responsibilities,
#'   country_mixture, stability, loo, loo_country, diagnostics,
#'   diagnostics_chains, res, settings)
kl_sweep_compare <- function(
  batch_dirs,
  mode_tol = 1.5,
  n_agree = 0L,
  stability_seed = 1L
) {
  inv <- kl_sweep_inventory(batch_dirs)
  message(sprintf(
    "[kl_sweep] fits present: %d of %d",
    sum(inv$present),
    nrow(inv)
  ))
  if (any(!inv$present)) {
    message(
      "[kl_sweep]   missing: ",
      paste(inv$id[!inv$present], collapse = ", ")
    )
  }
  use <- inv[inv$present, ]
  if (!nrow(use)) {
    stop("kl_sweep: no fits present.")
  }
  fits <- list()
  agree <- list()
  rows <- list()
  cls <- list()
  part <- list()
  dvar <- list()
  dch <- list()
  resp <- list()
  mixr <- list()
  stab <- list()
  for (i in seq_len(nrow(use))) {
    r <- use[i, ]
    f <- readRDS(r$fit_file)
    if (!inherits(f, "kl_fit")) {
      stop("kl_sweep: ", r$fit_file, " is not a kl_fit.")
    }
    modes <- kl_chain_modes(f, mode_tol)
    ag <- kl_agreeing_chains(f, mode_tol, n_agree, modes)
    fits[[r$id]] <- f
    agree[[r$id]] <- ag$chains
    rows[[r$id]] <- kl_fit_row(f, r, modes, ag)
    dt <- kl_fit_diagnostics_tables(f, r$id, ag)
    dvar[[r$id]] <- dt$variables
    dch[[r$id]] <- dt$chains
    sub <- kl_subset_chains(f, ag$chains)
    cls[[r$id]] <- kl_class_rows(sub, r)
    part[[r$id]] <- kl_partition_rows(sub, r)
    prior <- kl_sweep_prior_grouping(r$config)
    resp[[r$id]] <- dplyr::bind_rows(
      kl_responsibility_rows(f, r, "all", prior),
      kl_responsibility_rows(sub, r, "agree", prior)
    )
    mixr[[r$id]] <- kl_country_mixture_rows(sub, r, prior)
    stab[[r$id]] <- dplyr::bind_rows(
      kl_partition_stability(f, r, "all", stability_seed),
      kl_partition_stability(sub, r, "agree", stability_seed)
    )
  }
  lo <- kl_sweep_loo(fits, agree)
  tab_fit <- dplyr::bind_rows(rows)
  list(
    inv = inv,
    agree = agree,
    fits = tab_fit,
    classes = dplyr::bind_rows(cls),
    partition = dplyr::bind_rows(part),
    responsibilities = dplyr::bind_rows(resp),
    country_mixture = dplyr::bind_rows(mixr),
    stability = dplyr::bind_rows(stab),
    loo = dplyr::left_join(
      lo$table,
      dplyr::select(tab_fit, id, basis, K, mix),
      by = "id"
    ),
    loo_country = lo$country,
    diagnostics = dplyr::bind_rows(dvar),
    diagnostics_chains = dplyr::bind_rows(dch),
    res = lo$res,
    settings = list(
      batch_dirs = batch_dirs,
      mode_tol = mode_tol,
      n_agree = n_agree,
      exclude = KL_SWEEP_EXCL
    )
  )
}
