# 06_fit_figures_cb_adm1.R — per-ADM1 crossbasis exposure-response figures
#
# Sibling to 06_fit_figures_cb.R (unchanged by side-by-side rule).
# Requires a brms_cb fit with re_adm1_slope: true for at least one CB variable.
#
# For each CB variable that has ADM1 random slopes, combines the global fixed-
# effect draws with per-ADM1 RE draws to produce location-specific surfaces:
#   1. Faceted lag-response at median exposure (all ADM1 units on one plot)
#   2. Faceted cumulative exposure-response (all ADM1 units on one plot)
#
# Output: figures/cb/<run_tag>/adm1/<var>_lag_response_adm1.{png,pdf}
#                               adm1/<var>_cumulative_adm1.{png,pdf}
#
# Usage:
#   Rscript analysis/06_fit_figures_cb_adm1.R \
#     -C analysis/configs/test_cb_eth_reslope/test_cb_eth_reslope.yml
#   Rscript analysis/06_fit_figures_cb_adm1.R -C ... --format pdf --n_at 100

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(optparse)
  library(brms)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/plots.R")
source("analysis/utils/figures_io.R")

option_list <- list(
  make_option(c("-c", "--country"),
              default = NULL, type = "character",
              help = "ISO3 code(s) or 'all' — overrides config country [default: from config]"),
  make_option(c("--format"),
              default = "png", type = "character",
              help = "Output format: pdf | png [default: png]"),
  make_option(c("--dpi"),
              default = 150L, type = "integer",
              help = "DPI for raster output [default: 150]"),
  make_option(c("--n_at"),
              default = 50L, type = "integer",
              help = "Number of exposure grid points [default: 50]"),
  make_option(c("--vars"),
              default = NULL, type = "character",
              help = "Comma-separated CB variable names to plot (default: all)")
)

# Auto-detect occurrence vs size vs duration from the config (size and
# duration single-lag fits also carry ADM1 basis slopes and deserve per-ADM1
# figures — mirrors the dispatch in 06_fit_figures_cb.R).
.detect_model <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  cfg  <- NULL
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) { cfg <- argv[i + 1L]; break }
  }
  if (is.null(cfg) || !file.exists(cfg))
    stop("--config / -C is required and must point to an existing YAML.")
  detect_model_from_yaml(cfg)   # utils/config.R
}
model_type <- .detect_model()

opt <- make_options_from_config(option_list     = option_list,
                                enforce_options = "country",
                                model           = model_type)
stopifnot(opt[[model_type]]$backend == "brms_cb")
print_options(opt)

# File-kind prefix: size track uses size_fit / size_cb_meta; duration uses
# dur_fit / dur_cb_meta; monthly uses mon_fit / mon_cb_meta (routed to
# makeDurationRunFile() / makeMonthlyRunFile() by makeRunFile()).
kpfx <- if (model_type == "size") "size_" else if (model_type == "duration") "dur_" else if (model_type == "monthly") "mon_" else ""

# Locate files ----------------------------------------------------------------

fit_file  <- makeRunFile(opt, kind = paste0(kpfx, "fit"))
meta_file <- makeRunFile(opt, kind = paste0(kpfx, "cb_meta"))

if (!file.exists(fit_file))
  stop("Fit file not found: ", fit_file,
       "\n  Run 04_fit_run.R with backend: brms_cb first.")
if (!file.exists(meta_file))
  stop("CB meta file not found: ", meta_file)

message("Loading fit: ", fit_file)
brms_fit <- readRDS(fit_file)

message("Loading CB meta: ", meta_file)
cb_meta <- readRDS(meta_file)

# Output directory ------------------------------------------------------------

fmt     <- match.arg(opt$format %||% "png", c("pdf", "png"))
# makeFigureDir() (utils/paths.R) recovers -C/--config from argv
# itself when config_path is omitted; subdir = "adm1" reproduces this
# script's figures/cb/<run_tag>/adm1/ nesting exactly.
fig_dir <- makeFigureDir(opt = list(), track = "cb", subdir = "adm1")

save_fig <- figureSaver(fig_dir, width = 11, height = 8,
                         dpi = as.integer(opt$dpi %||% 150L), formats = fmt)

# Which variables to process --------------------------------------------------

vars_req <- if (!is.null(opt$vars)) {
  intersect(trimws(strsplit(opt$vars, ",")[[1]]), names(cb_meta))
} else {
  names(cb_meta)
}
if (length(vars_req) == 0) stop("No matching CB variables.")

# ADM1 RE structure -----------------------------------------------------------

re_groups <- tryCatch(names(brms::ranef(brms_fit)), error = function(e) character(0))
# Accept bare "ADM1" (occurrence uncorrelated slopes) or any nested group that
# ends with ":ADM1" (size: "AFRO_region:country:ADM1").
adm1_group <- if ("ADM1" %in% re_groups) {
  "ADM1"
} else {
  candidates <- re_groups[endsWith(re_groups, ":ADM1") | endsWith(re_groups, "/ADM1")]
  if (length(candidates) == 0L)
    stop("No ADM1 random effects found in fit.\n",
         "Re-run with re_adm1_slope: true for at least one CB variable.")
  candidates[which.max(nchar(candidates))]   # deepest nesting
}
message("Using RE group for ADM1 slopes: ", adm1_group)

# Parent groups for the hierarchical-slope sum. With full nesting the per-ADM1
# slope is FE + r_region + r_country + r_adm1, so we also need the region and
# country group names (derived from the ADM1 group string — everything but
# the last colon/slash-separated part is the country group, everything but
# the last two is the region group). If a parent group is absent (e.g.
# grouping mode drops region, so a 2-part "country:ADM1" group has no region
# parent; or a flat single-level fit has no parent at all), its contribution
# is treated as 0. Grouping-mode fits (adm1_group = "country:ADM1", 2 parts)
# resolve country_group to bare "country" here (>=2 parts, not >=3) -- a
# >=3 threshold falls back to NA for a 2-part grouping-mode group, silently
# dropping the country-level RE contribution from every grouping-mode
# per-ADM1 estimate.
grp_parts     <- strsplit(adm1_group, "[:/]")[[1]]
n_parts       <- length(grp_parts)
country_group <- if (n_parts >= 2) paste(grp_parts[seq_len(n_parts - 1)], collapse = ":") else NA_character_
region_group  <- if (n_parts >= 3) paste(grp_parts[seq_len(n_parts - 2)], collapse = ":") else NA_character_

adm1_levels <- rownames(brms::ranef(brms_fit)[[adm1_group]])
message(sprintf("ADM1 levels (%d): %s",
                length(adm1_levels), paste(adm1_levels, collapse = ", ")))

n_at    <- as.integer(opt$n_at %||% 50L)
probs   <- c(0.025, 0.975)   # outer 95% CI
probs50 <- c(0.25, 0.75)     # inner 50% CrI

# Per-variable loop -----------------------------------------------------------

for (v in vars_req) {
  vm           <- cb_meta[[v]]
  cb_col_names <- vm$cb_col_names
  if (is.null(cb_col_names) || length(cb_col_names) == 0) next
  n_cols <- length(cb_col_names)

  # Check that this variable actually has ADM1 RE slopes in the fit.
  # brms replaces spaces with dots in posterior parameter names, so we must
  # apply the same substitution when building lookup strings.
  adm1_pname <- function(lvl) gsub(" ", ".", lvl, fixed = TRUE)
  test_re_name <- sprintf("r_%s[%s,%s]",
                          adm1_group, adm1_pname(adm1_levels[1]), cb_col_names[1])
  re_check <- tryCatch(
    brms::as_draws_matrix(brms_fit, variable = test_re_name),
    error = function(e) NULL
  )
  if (is.null(re_check) || nrow(re_check) == 0) {
    message(sprintf("[%s] No ADM1 RE slopes found — skipping per-ADM1 figures.", v))
    next
  }

  message(sprintf("[%s] extracting per-ADM1 surfaces ...", v))

  # Global FE draws: n_draws × n_cols
  b_names  <- sprintf("b_%s", cb_col_names)
  fe_draws <- unclass(brms::as_draws_matrix(brms_fit, variable = b_names))
  n_draws  <- nrow(fe_draws)

  # Hierarchical RE draws for one ADM1 level: sum the region + country + ADM1
  # slope deviations (each n_draws × n_cols). Missing parent groups -> 0.
  grab_re <- function(group, level) {
    vnames <- sprintf("r_%s[%s,%s]", group, adm1_pname(level), cb_col_names)
    tryCatch(unclass(brms::as_draws_matrix(brms_fit, variable = vnames)),
             error = function(e) matrix(0, n_draws, n_cols))
  }
  nested_re_draws <- function(adm1_level) {
    lin <- adm1_lineage_levels(adm1_level)
    acc <- grab_re(adm1_group, lin$adm1)
    if (!is.na(country_group) && country_group %in% re_groups)
      acc <- acc + grab_re(country_group, lin$country)
    if (!is.na(region_group) && region_group %in% re_groups)
      acc <- acc + grab_re(region_group, lin$region)
    acc
  }

  # Prediction grid
  rng     <- vm$range %||% c(NA_real_, NA_real_)
  at_vals <- seq(rng[1], rng[2], length.out = n_at)
  mid_idx <- which.min(abs(at_vals - median(at_vals)))

  # ---- Single-lag fit: no lag axis → faceted per-ADM1 exposure-response ----
  if (!is.null(vm$single_lag)) {
    L <- as.integer(vm$single_lag)
    # onebasis stub: reattach stored attributes so crosspred can reconstruct it
    ob_stub <- structure(matrix(0, nrow = 1L, ncol = n_cols),
                          class = c("onebasis", "matrix"))
    for (a in names(vm$onebasis_attrs)) attr(ob_stub, a) <- vm$onebasis_attrs[[a]]

    erf_rows <- list()
    for (adm1 in adm1_levels) {
      re_draws <- nested_re_draws(adm1)
      combined <- fe_draws + re_draws
      draw_preds <- vector("list", n_draws)
      for (i in seq_len(n_draws)) {
        draw_preds[[i]] <- tryCatch(
          dlnm::crosspred(ob_stub,
                          coef = as.numeric(combined[i, ]),
                          vcov = diag(0, n_cols),
                          at   = at_vals,
                          cen  = 0),
          error = function(e) NULL
        )
      }
      draw_preds <- Filter(Negate(is.null), draw_preds)
      if (length(draw_preds) == 0) {
        warning(sprintf("All draws failed for ADM1='%s' — skipping.", adm1)); next
      }
      erf_mat <- do.call(cbind, lapply(draw_preds, function(x) as.numeric(x$matfit)))
      erf_rows[[adm1]] <- data.frame(
        ADM1     = adm1,
        exposure = at_vals,
        logRR    = apply(erf_mat, 1, median),
        lo       = apply(erf_mat, 1, quantile, probs = probs[1]),
        hi       = apply(erf_mat, 1, quantile, probs = probs[2]),
        lo50     = apply(erf_mat, 1, quantile, probs = probs50[1]),
        hi50     = apply(erf_mat, 1, quantile, probs = probs50[2])
      )
      message(sprintf("  %s: %d/%d successful draws", adm1, length(draw_preds), n_draws))
    }
    if (length(erf_rows) == 0) {
      warning(sprintf("[%s] No ADM1 ERFs extracted.", v)); next
    }
    erf_df <- dplyr::bind_rows(erf_rows)
    p_erf <- ggplot(erf_df, aes(x = exposure)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
      geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#2166ac") +
      geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      facet_wrap(~ ADM1) +
      labs(
        title    = sprintf("Exposure-response per ADM1 (single lag %d): %s", L, v),
        subtitle = "Posterior median + 50% / 95% CrI",
        x = v, y = "log-RR"
      ) +
      theme_bw(base_size = 10)
    save_fig(p_erf, sprintf("%s_erf_adm1_%s", v, slTag(L)))
    next
  }

  # CB stub
  cb_stub <- structure(
    matrix(0, nrow = 1L, ncol = n_cols),
    class  = c("crossbasis", "matrix"),
    argvar = vm$argvar,
    arglag = vm$arglag,
    lag    = vm$lag,
    df     = vm$df,
    range  = vm$range
  )

  lag_rows  <- list()
  cum_rows  <- list()
  heat_rows <- list()

  for (adm1 in adm1_levels) {
    # Hierarchical ADM1 RE draws: FE + region + country + ADM1 (n_draws × n_cols)
    re_draws <- nested_re_draws(adm1)

    # Combine FE + RE
    combined <- fe_draws + re_draws

    # Per-draw crosspred
    draw_preds <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      draw_preds[[i]] <- tryCatch(
        dlnm::crosspred(cb_stub,
                        coef  = as.numeric(combined[i, ]),
                        vcov  = diag(0, n_cols),
                        at    = at_vals,
                        cen   = 0,
                        cumul = TRUE),
        error = function(e) NULL
      )
    }
    draw_preds <- Filter(Negate(is.null), draw_preds)
    if (length(draw_preds) == 0) {
      warning(sprintf("All draws failed for ADM1='%s' — skipping.", adm1))
      next
    }

    # Shared quantities across the three figure types
    mat_list <- lapply(draw_preds, function(x) x$matfit)
    mat_arr  <- simplify2array(mat_list)
    lag_seq  <- as.numeric(sub("^lag", "", colnames(draw_preds[[1]]$matfit)))
    mat_med  <- apply(mat_arr, c(1, 2), median)   # at × lag posterior median surface

    # Lag-response at median exposure
    lag_rows[[adm1]] <- data.frame(
      ADM1  = adm1,
      lag   = lag_seq,
      logRR = apply(mat_arr, c(1, 2), median)[mid_idx, ],
      lo    = apply(mat_arr, c(1, 2), quantile, probs = probs[1])[mid_idx, ],
      hi    = apply(mat_arr, c(1, 2), quantile, probs = probs[2])[mid_idx, ],
      lo50  = apply(mat_arr, c(1, 2), quantile, probs = probs50[1])[mid_idx, ],
      hi50  = apply(mat_arr, c(1, 2), quantile, probs = probs50[2])[mid_idx, ]
    )

    # Cumulative exposure-response
    cum_list <- lapply(draw_preds, function(x) x$allfit)
    cum_arr  <- do.call(cbind, lapply(cum_list, as.numeric))
    cum_rows[[adm1]] <- data.frame(
      ADM1     = adm1,
      exposure = at_vals,
      cumRR    = apply(cum_arr, 1, median),
      lo       = apply(cum_arr, 1, quantile, probs = probs[1]),
      hi       = apply(cum_arr, 1, quantile, probs = probs[2]),
      lo50     = apply(cum_arr, 1, quantile, probs = probs50[1]),
      hi50     = apply(cum_arr, 1, quantile, probs = probs50[2])
    )

    # Heatmap surface (full at × lag grid)
    colnames(mat_med) <- lag_seq
    rownames(mat_med) <- round(at_vals, 3)
    heat_rows[[adm1]] <- as.data.frame(mat_med) %>%
      tibble::rownames_to_column("exposure") %>%
      tidyr::pivot_longer(-exposure, names_to = "lag", values_to = "logRR") %>%
      dplyr::mutate(
        ADM1     = adm1,
        exposure = as.numeric(exposure),
        lag      = as.numeric(lag)
      )

    message(sprintf("  %s: %d/%d successful draws", adm1, length(draw_preds), n_draws))
  }

  if (length(lag_rows) == 0) {
    warning(sprintf("[%s] No ADM1 surfaces extracted.", v))
    next
  }

  lag_df  <- dplyr::bind_rows(lag_rows)
  cum_df  <- dplyr::bind_rows(cum_rows)
  heat_df <- dplyr::bind_rows(heat_rows)
  med_exp <- round(at_vals[mid_idx], 2)

  # Shared colour scale limits across all ADM1 panels (symmetric around 0)
  heat_lim <- max(abs(heat_df$logRR), na.rm = TRUE)

  # ---- Faceted lag-response ----
  p_lag <- ggplot(lag_df, aes(x = lag)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#2166ac") +
    geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ ADM1) +
    labs(
      title    = sprintf("Lag-response per ADM1: %s", v),
      subtitle = sprintf("At median exposure (%s = %.2f) — posterior median + 50%% / 95%% CrI",
                         v, med_exp),
      x = "Lag (months)",
      y = "log-RR"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_lag, sprintf("%s_lag_response_adm1", v))

  # ---- Faceted cumulative exposure-response ----
  p_cum <- ggplot(cum_df, aes(x = exposure)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#d73027") +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#d73027") +
    geom_line(aes(y = cumRR), colour = "#d73027", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ ADM1) +
    labs(
      title    = sprintf("Cumulative exposure-response per ADM1: %s", v),
      subtitle = "Summed over all lags — posterior median + 50% / 95% CrI",
      x        = v,
      y        = "Cumulative log-RR"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_cum, sprintf("%s_cumulative_adm1", v))

  # ---- Faceted heatmap (exposure × lag surface per ADM1) ----
  p_heat <- ggplot(heat_df, aes(x = exposure, y = lag, fill = logRR)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d73027",
      midpoint = 0,
      limits   = c(-heat_lim, heat_lim),
      name     = "log-RR"
    ) +
    facet_wrap(~ ADM1) +
    labs(
      title    = sprintf("Exposure-lag surface per ADM1: %s", v),
      subtitle = "Posterior median log-RR",
      x        = v,
      y        = "Lag (months)"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_heat, sprintf("%s_heatmap_adm1", v), width = 13, height = 9)
}

message(sprintf("Done. Per-ADM1 figures written to %s/", fig_dir))
