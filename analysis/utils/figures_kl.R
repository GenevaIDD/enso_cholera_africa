# utils/figures_kl.R -- ggplot builders for the brms_kl class-sweep figures
#
# Pure builders (they return ggplot objects and never save) behind
# 06_fit_figures_kl_sweep.R. Batch-level figures read only the tables written
# by 05_fit_postprocess_kl_sweep.R (utils/paths.R::makeKlSweepRunFile());
# per-fit figures read one kl_fit plus its kl_meta sidecar.
#
# Dependencies (call-time; modules never source() each other): base.R (%||%),
# config.R (load_grouping()), dlnm_kernel.R (.kl_draw_array(),
# .kl_basis_centred()), kl_sweep.R (kl_chain_modes()), figures_cb.R
# (load_taxdat_afr_sf(), membership maps only). Packages via `::`: ggplot2,
# dplyr, tidyr, tibble, posterior, sf.

#' Class colours, keyed by the class label without "enso_": the enso_grouping
#' (_kl) sidecar colours, so a class has one colour in every sweep figure.
KL_CLASS_COL <- c(east = "#1b9e77", horn = "#d95f02", west = "#66a61e",
                  nga = "#e6ab02", south = "#7570b3", central = "#a6761d")
KL_PRIOR_COL <- c(reg = "#2a78d6", mix = "#eb6834")
KL_CHAIN_COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#eda100", "#7a5cc7",
                  "#52514e", "#9a9892")

#' Shared theme of the kl figures.
kl_theme <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#e6e5e0", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   legend.position = "bottom")
}

#' Model labels for the sweep tables: K, prior (registry / + contiguity mix)
#' and basis as factors in a stable order.
#' @noRd
.kl_sweep_labels <- function(d) {
  dplyr::mutate(d,
    prior = factor(ifelse(as.logical(mix), "mix", "reg"), c("reg", "mix")),
    basis = factor(basis, intersect(c("ns3", "ns4", "ns5", "bs3", "bs4", "bs5"), unique(basis))),
    K = as.integer(K))
}

# --- batch-level -------------------------------------------------------------------

#' Country-block elpd difference to the best fit (+-2 SE of the difference) by
#' K, prior and basis; `chains` "agree" (default) or "all".
plot_kl_sweep_elpd <- function(loo, chains = "agree") {
  d <- .kl_sweep_labels(dplyr::filter(loo, .data$chains == !!chains))
  ggplot2::ggplot(d, ggplot2::aes(K, elpd_diff_loo, colour = prior)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = elpd_diff_loo - 2 * se_diff_loo,
                                          ymax = elpd_diff_loo + 2 * se_diff_loo),
                             position = ggplot2::position_dodge(width = 0.4), size = 0.3) +
    ggplot2::facet_wrap(~basis, nrow = 1) +
    ggplot2::scale_colour_manual(values = KL_PRIOR_COL, name = "membership prior",
                                 labels = c(reg = "registry", mix = "+ contiguity mix")) +
    ggplot2::labs(x = "number of classes K", y = "elpd difference to the best fit (± 2 SE)",
                  title = "Country-block PSIS-LOO of the class sweep",
                  subtitle = sprintf("%s chains; descriptive (one block per country, high Pareto k)",
                                     if (chains == "agree") "agreeing" else "all")) +
    kl_theme()
}

#' Stacking weights of the sweep fits.
plot_kl_sweep_weights <- function(loo, chains = "agree") {
  d <- .kl_sweep_labels(dplyr::filter(loo, .data$chains == !!chains))
  ggplot2::ggplot(d, ggplot2::aes(factor(K), stacking_weight, fill = prior)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
    ggplot2::facet_wrap(~basis, nrow = 1) +
    ggplot2::scale_fill_manual(values = KL_PRIOR_COL, name = "membership prior",
                               labels = c(reg = "registry", mix = "+ contiguity mix")) +
    ggplot2::labs(x = "number of classes K", y = "stacking weight",
                  title = "Stacking weights over the sweep fits") +
    kl_theme()
}

#' Class centres (kmu, mean and 95% interval, agreeing chains) per fit.
plot_kl_sweep_centres <- function(classes) {
  d <- .kl_sweep_labels(dplyr::filter(classes, param == "kmu"))
  ggplot2::ggplot(d, ggplot2::aes(interaction(prior, K, sep = " K"), mean, colour = class)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = q2.5, ymax = q97.5),
                             position = ggplot2::position_dodge(width = 0.6), size = 0.25) +
    ggplot2::facet_wrap(~basis, nrow = 1) +
    ggplot2::scale_colour_manual(values = KL_CLASS_COL, name = "class") +
    ggplot2::scale_y_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 5)) +
    ggplot2::labs(x = "membership prior / K", y = "class centre (lag, months; negative = lead)",
                  title = "Class lag centres by fit") +
    kl_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Convergence by fit: max R-hat and min bulk ESS, all chains vs agreeing chains.
plot_kl_sweep_convergence <- function(fits, rhat_line = 1.01) {
  d <- .kl_sweep_labels(fits) |>
    dplyr::transmute(id, basis, K, prior,
                     `max R-hat, all chains` = max_rhat_all, `max R-hat, agreeing chains` = max_rhat_agree,
                     `min bulk ESS, all chains` = min_ess_bulk_all,
                     `min bulk ESS, agreeing chains` = min_ess_bulk_agree) |>
    tidyr::pivot_longer(-c(id, basis, K, prior)) |>
    dplyr::mutate(stat = sub(",.*", "", name), chains = sub(".*, ", "", name))
  ggplot2::ggplot(d, ggplot2::aes(interaction(prior, K, sep = " K"), value, colour = chains)) +
    ggplot2::geom_point(size = 1.6, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_hline(data = tibble::tibble(stat = "max R-hat", y = rhat_line),
                        ggplot2::aes(yintercept = y), linetype = "dashed", colour = "#52514e") +
    ggplot2::facet_grid(stat ~ basis, scales = "free_y") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c(`all chains` = "#9a9892", `agreeing chains` = "#2a78d6"),
                                 name = NULL) +
    ggplot2::labs(x = "membership prior / K", y = NULL, title = "Convergence of the sweep fits") +
    kl_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' MAP partition per fit (country x fit; fill = class, opacity = probability).
plot_kl_sweep_partition <- function(partition, fits) {
  ord <- .kl_sweep_labels(fits) |> dplyr::arrange(basis, K, prior)
  d <- dplyr::mutate(partition, id = factor(id, ord$id),
                     group = factor(group, intersect(names(KL_CLASS_COL), unique(group))))
  ggplot2::ggplot(d, ggplot2::aes(id, iso3, fill = group, alpha = p)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = KL_CLASS_COL, name = "class") +
    ggplot2::scale_alpha_continuous(range = c(0.25, 1), limits = c(0.25, 1), oob = scales::squish,
                                    name = "probability") +
    ggplot2::labs(x = NULL, y = NULL, title = "Most probable class per country (agreeing chains)") +
    kl_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                                panel.grid.major = ggplot2::element_blank())
}

# --- per fit -------------------------------------------------------------------------

#' Short class labels of a fit.
#' @noRd
.kl_fit_labels <- function(fit) sub("^enso_", "", fit$class_labels)

#' Traces of the class centres and El Nino contrasts, by chain.
plot_kl_traces <- function(fit) {
  K <- fit$K; lab <- .kl_fit_labels(fit)
  tr <- posterior::as_draws_df(posterior::subset_draws(
          fit$draws, variable = c(sprintf("kmu[%d]", 1:K), sprintf("kcon[%d]", 1:K)))) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(dplyr::starts_with("k"), names_to = "variable") |>
    dplyr::mutate(k = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable)),
                  what = ifelse(startsWith(variable, "kmu"), "centre (lag)", "El Niño contrast (log RR)"),
                  class = factor(lab[k], lab), chain = factor(.chain))
  ggplot2::ggplot(tr, ggplot2::aes(.iteration, value, colour = chain)) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_grid(what ~ class, scales = "free_y") +
    ggplot2::scale_colour_manual(values = KL_CHAIN_COL, name = "chain") +
    ggplot2::labs(x = "draw (post-warmup)", y = NULL,
                  title = "Traces of the class centres and El Niño contrasts") +
    kl_theme()
}

#' Posterior of each class centre, by chain.
plot_kl_centres_by_chain <- function(fit) {
  K <- fit$K; lab <- .kl_fit_labels(fit)
  d <- posterior::as_draws_df(posterior::subset_draws(fit$draws, variable = sprintf("kmu[%d]", 1:K))) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(dplyr::starts_with("kmu"), names_to = "variable") |>
    dplyr::mutate(k = as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable)),
                  class = factor(lab[k], lab), chain = factor(.chain))
  ggplot2::ggplot(d, ggplot2::aes(value, fill = chain)) +
    ggplot2::geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.55, colour = NA) +
    ggplot2::facet_wrap(~class, nrow = 1) +
    ggplot2::scale_fill_manual(values = KL_CHAIN_COL, name = "chain") +
    ggplot2::scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 5)) +
    ggplot2::labs(x = "class centre (lag, months; negative = ENSO leads)", y = "draws",
                  title = "Posterior of each class centre, by chain") +
    kl_theme()
}

#' Class membership maps: prior partition, all chains, dominant mode, and each
#' chain. Fill = most probable class, opacity = its probability, dark outline
#' = seed countries.
#' @param meta the fit's kl_meta sidecar; @param afr load_taxdat_afr_sf().
plot_kl_membership_maps <- function(fit, meta, afr, mode_tol = 1.5) {
  K <- fit$K; lab <- fit$class_labels; ctry <- fit$country_levels; C <- length(ctry)
  ch <- fit$draws$.chain
  md <- kl_chain_modes(fit, mode_tol)
  R  <- exp(.kl_draw_array(fit, "kl_log_resp", c(K, C)))
  resp_by <- function(rows, panel) {
    Rm <- apply(R[rows, , , drop = FALSE], c(2, 3), mean)
    tibble::tibble(iso3 = ctry, group = lab[apply(Rm, 2, which.max)], p = apply(Rm, 2, max), panel = panel)
  }
  src <- meta$kcfg$prior_grouping$source %||% "enso_grouping"
  reg <- load_grouping(src)
  dom <- md$chains[md$mode == 1L]
  maps <- dplyr::bind_rows(
    tibble::tibble(iso3 = ctry, group = reg$group[match(ctry, reg$iso3)], p = 1,
                   panel = sprintf("prior partition (%s)", src)),
    resp_by(seq_along(ch), "posterior, all chains"),
    resp_by(which(ch %in% dom), sprintf("dominant mode (chains %s)", paste(dom, collapse = ","))),
    dplyr::bind_rows(lapply(md$chains, function(i)
      resp_by(which(ch == i), sprintf("chain %d (mode %d)", i, md$mode[[as.character(i)]])))))
  panels <- unique(maps$panel)
  geo <- tidyr::crossing(sf::st_drop_geometry(afr)[, "country", drop = FALSE], panel = panels) |>
    dplyr::rename(iso3 = country) |>
    dplyr::left_join(maps, by = c("iso3", "panel")) |>
    dplyr::left_join(dplyr::select(afr, iso3 = country), by = "iso3") |>
    sf::st_as_sf() |>
    dplyr::mutate(panel = factor(panel, panels), class = sub("^enso_", "", group))
  seeds <- unlist(meta$seeds)
  ggplot2::ggplot(geo) +
    ggplot2::geom_sf(ggplot2::aes(fill = class, alpha = p), colour = "white", linewidth = 0.15) +
    ggplot2::geom_sf(data = dplyr::filter(geo, iso3 %in% seeds), fill = NA, colour = "#0b0b0b",
                     linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = KL_CLASS_COL, na.value = "#e6e5e0", name = "class",
                               na.translate = FALSE) +
    ggplot2::scale_alpha_continuous(range = c(0.2, 1), limits = c(0.25, 1), oob = scales::squish,
                                    name = "probability") +
    ggplot2::facet_wrap(~panel, ncol = 3) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' Class exposure-response (log RR of a sustained ENSO value vs 0 under each
#' class's own kernel), dominant-mode chains.
plot_kl_class_erf <- function(fit, meta, mode_tol = 1.5, n = 60L) {
  K <- fit$K; lab <- .kl_fit_labels(fit); J <- fit$data_dims$J
  md  <- kl_chain_modes(fit, mode_tol)
  dom <- md$chains[md$mode == 1L]
  rng <- meta$meta$range
  at  <- seq(max(-2, rng[1]), min(2.5, rng[2]), length.out = n)
  Bc  <- .kl_basis_centred(meta$meta, at, 0)
  dm  <- posterior::as_draws_matrix(posterior::subset_draws(fit$draws, chain = dom))
  erf <- dplyr::bind_rows(lapply(seq_len(K), function(k) {
    e <- dm[, sprintf("kbeta[%d,%d]", seq_len(J), k), drop = FALSE] %*% t(Bc)
    tibble::tibble(class = lab[k], x = at, mean = colMeans(e),
                   lo = apply(e, 2, stats::quantile, .025), hi = apply(e, 2, stats::quantile, .975))
  })) |> dplyr::mutate(class = factor(class, lab))
  ggplot2::ggplot(erf, ggplot2::aes(x, mean)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = class), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(colour = class), linewidth = 0.8) +
    ggplot2::facet_wrap(~class, nrow = 1) +
    ggplot2::scale_colour_manual(values = KL_CLASS_COL, guide = "none") +
    ggplot2::scale_fill_manual(values = KL_CLASS_COL, guide = "none") +
    ggplot2::labs(x = "sustained ENSO index (ONI)", y = "log RR vs ONI 0",
                  title = sprintf("Class exposure-response, dominant mode (chains %s)",
                                  paste(dom, collapse = ","))) +
    kl_theme()
}

#' Posterior of the country / ADM1 deviation sds of the class ENSO
#' coefficients (dashed: the ADM1 half-normal prior, or the country one when
#' kernel.re_slope.adm1 is off).
plot_kl_re_sds <- function(fit, meta) {
  vv <- grep("^kre_sd[ca]\\[", posterior::variables(fit$draws), value = TRUE)
  if (!length(vv)) return(NULL)
  adm1    <- any(startsWith(vv, "kre_sda"))
  sd_pr   <- as.numeric(meta$kcfg$re_slope[[if (adm1) "sd_adm1" else "sd_country"]] %||% 0.35)
  re <- posterior::as_draws_df(posterior::subset_draws(fit$draws, variable = vv)) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(dplyr::starts_with("kre"), names_to = "variable") |>
    dplyr::mutate(level = ifelse(startsWith(variable, "kre_sdc"), "country", "ADM1"),
                  coef = paste("coefficient", sub(".*\\[(\\d+)\\]", "\\1", variable)))
  ggplot2::ggplot(re, ggplot2::aes(value)) +
    ggplot2::geom_density(ggplot2::aes(fill = level, colour = level), alpha = 0.3) +
    ggplot2::stat_function(fun = function(x) 2 * stats::dnorm(x, 0, sd_pr),
                           linetype = "dashed", colour = "#52514e", n = 200) +
    ggplot2::facet_wrap(~coef, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(country = "#2a78d6", ADM1 = "#eb6834"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(country = "#2a78d6", ADM1 = "#eb6834"), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 2.2)) +
    ggplot2::labs(x = "sd of the deviations", y = "density",
                  title = sprintf("ENSO-coefficient deviation sds (dashed: %s half-normal(0, %s) prior)",
                                  if (adm1) "ADM1" else "country", sd_pr)) +
    kl_theme()
}

#' Country lag centres of a kernel.lag_re fit: per country its own lag,
#' membership-weighted per draw, sum_k r_ck cent_kc (mean and 95% interval),
#' coloured by the country's modal class; cent_kc = kmu_k + d_c (d_c = kclag,
#' lags) under kernel.lag_re.scale month, the drawn kccent[k, c] under scale
#' logit (only its program has kccent; d_c then in logit units). NULL when the
#' fit has no country lag offsets.
plot_kl_country_lag <- function(fit) {
  d <- as.data.frame(fit$draws); K <- fit$K; C <- fit$country_levels
  if (!"kclag[1]" %in% names(d)) return(NULL)
  logit <- "kccent[1,1]" %in% names(d)
  lab <- sub("^enso_", "", fit$class_labels)
  rows <- lapply(seq_along(C), function(c) {
    P   <- sapply(seq_len(K), function(k) exp(d[[sprintf("kl_log_resp[%d,%d]", k, c)]]))
    lag <- if (logit) rowSums(P * sapply(seq_len(K), function(k) d[[sprintf("kccent[%d,%d]", k, c)]]))
           else rowSums(P * sapply(seq_len(K), function(k) d[[sprintf("kmu[%d]", k)]])) + d[[sprintf("kclag[%d]", c)]]
    tibble::tibble(iso3 = C[c], class = lab[which.max(colMeans(P))], lag = mean(lag),
                   lo = stats::quantile(lag, 0.025, names = FALSE), hi = stats::quantile(lag, 0.975, names = FALSE),
                   offset = mean(d[[sprintf("kclag[%d]", c)]]))
  })
  df <- dplyr::bind_rows(rows)
  df$iso3 <- factor(df$iso3, levels = df$iso3[order(df$class, df$lag)])
  sdl <- if ("klag_sd[1]" %in% names(d)) stats::median(d[["klag_sd[1]"]]) else NA_real_
  ggplot2::ggplot(df, ggplot2::aes(lag, iso3, colour = class)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#c3c2bd") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "country lag centre (class lag + country offset; months)", y = NULL, colour = "modal class",
                  title = sprintf("Country lag centres (kernel.lag_re %s; offset sd median %.2f%s)",
                                  if (logit) "logit" else "month", sdl, if (logit) " logit units" else " lags")) +
    kl_theme()
}

