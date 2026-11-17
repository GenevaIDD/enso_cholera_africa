# 06_fit_figures_cb_country.R — per-country crossbasis exposure-response figures
#
# Usage:
#   Rscript analysis/06_fit_figures_cb_country.R \
#     -C analysis/configs/test_cb_eth_reslope/test_cb_eth_reslope.yml
#   Rscript analysis/06_fit_figures_cb_country.R -C <config.yml> --format pdf \
#     --n_at 100
#
# Details:
#   Requires a brms_cb fit with re_adm1_slope: true for at least one CB
#   variable (nested_slope_re_terms() then also emits a country-level slope).
#   Per variable with country slopes: faceted lag-response, cumulative ERF,
#   heatmap; single-lag fits: faceted ERF only. Writes
#   figures/cb/<run_tag>/country/<var>_*_country*.{png,pdf}.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(optparse)
  library(brms)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/plots.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-c", "--country"),
    default = NULL,
    type = "character",
    help = "ISO3 code(s) or 'all' — overrides config country [default: from config]"
  ),
  make_option(
    c("--format"),
    default = "png",
    type = "character",
    help = "Output format: pdf | png [default: png]"
  ),
  make_option(
    c("--dpi"),
    default = 150L,
    type = "integer",
    help = "DPI for raster output [default: 150]"
  ),
  make_option(
    c("--n_at"),
    default = 50L,
    type = "integer",
    help = "Number of exposure grid points [default: 50]"
  ),
  make_option(
    c("--vars"),
    default = NULL,
    type = "character",
    help = "Comma-separated CB variable names to plot (default: all)"
  )
)

# Auto-detect occurrence vs size vs duration from the config — mirrors the
# dispatch in 06_fit_figures_cb.R / 06_fit_figures_cb_adm1.R.
.detect_model <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  cfg <- NULL
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) {
      cfg <- argv[i + 1L]
      break
    }
  }
  if (is.null(cfg) || !file.exists(cfg)) {
    stop("--config / -C is required and must point to an existing YAML.")
  }
  detect_model_from_yaml(cfg) # utils/config.R
}
model_type <- .detect_model()

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = "country",
  model = model_type
)
stopifnot(opt[[model_type]]$backend == "brms_cb")
print_options(opt)

# File-kind prefix: size track uses size_fit / size_cb_meta; duration uses
# dur_fit / dur_cb_meta; monthly uses mon_fit / mon_cb_meta (routed to
# make_duration_run_file() / make_monthly_run_file() by make_run_file()).
kpfx <- if (model_type == "size") {
  "size_"
} else if (model_type == "duration") {
  "dur_"
} else if (model_type == "monthly") {
  "mon_"
} else {
  ""
}

# 3. Locate files --------------------------------------------------------------

fit_file <- make_run_file(opt, kind = paste0(kpfx, "fit"))
meta_file <- make_run_file(opt, kind = paste0(kpfx, "cb_meta"))

if (!file.exists(fit_file)) {
  stop(
    "Fit file not found: ",
    fit_file,
    "\n  Run 04_fit_run.R with backend: brms_cb first."
  )
}
if (!file.exists(meta_file)) {
  stop("CB meta file not found: ", meta_file)
}

message("Loading fit: ", fit_file)
brms_fit <- readRDS(fit_file)

message("Loading CB meta: ", meta_file)
cb_meta <- readRDS(meta_file)

# 4. Output directory ----------------------------------------------------------

fmt <- match.arg(opt$format %||% "png", c("pdf", "png"))
.cfg_path <- local({
  argv <- commandArgs(trailingOnly = TRUE)
  cfg <- NULL
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) {
      cfg <- argv[i + 1L]
      break
    }
  }
  cfg
})
run_tag <- basename(tools::file_path_sans_ext(.cfg_path %||% "run_cb"))
fig_dir <- file.path("figures", "cb", run_tag, "country")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

save_fig <- function(p, name, width = 11, height = 8) {
  path <- file.path(fig_dir, paste0(name, ".", fmt))
  ggplot2::ggsave(
    path,
    p,
    width = width,
    height = height,
    dpi = as.integer(opt$dpi %||% 150L)
  )
  message("Saved: ", path)
  invisible(path)
}

# 5. Which variables to process ------------------------------------------------

vars_req <- if (!is.null(opt$vars)) {
  intersect(trimws(strsplit(opt$vars, ",")[[1]]), names(cb_meta))
} else {
  names(cb_meta)
}
if (length(vars_req) == 0) {
  stop("No matching CB variables.")
}

# 6. Country RE structure ------------------------------------------------------
# Country surfaces combine the global fixed-effect draws with per-country RE
# draws (region + country, i.e. one level shallower than the ADM1 script — no
# ADM1-level deviation is added).

re_groups <- tryCatch(names(brms::ranef(brms_fit)), error = function(e) {
  character(0)
})
# Accept bare "country" (grouping mode — AFRO_region dropped from the RE
# formula, so "country" is the top-level slope group; ENSO best-lag fits are
# all like this) or any nested group that ends with ":country" (ungrouped
# runs: "AFRO_region:country").
country_group <- if ("country" %in% re_groups) {
  "country"
} else {
  candidates <- re_groups[
    endsWith(re_groups, ":country") | endsWith(re_groups, "/country")
  ]
  if (length(candidates) == 0L) {
    stop(
      "No country-level random effects found in fit.\n",
      "Re-run with re_adm1_slope: true for at least one CB variable — ",
      "nested_slope_re_terms() emits a country-level slope term alongside ",
      "the ADM1-level one."
    )
  }
  candidates[which.max(nchar(candidates))] # deepest nesting
}
message("Using RE group for country slopes: ", country_group)

# Parent group for the hierarchical-slope sum. With region present the
# per-country slope is FE + r_region + r_country (derived from the country
# group string). If the parent group is absent (grouping mode, or a flat
# single-level fit), its contribution is treated as 0.
grp_parts <- strsplit(country_group, "[:/]")[[1]]
region_group <- if (length(grp_parts) >= 2) grp_parts[1] else NA_character_

country_levels <- rownames(brms::ranef(brms_fit)[[country_group]])
message(sprintf(
  "Country levels (%d): %s",
  length(country_levels),
  paste(country_levels, collapse = ", ")
))

n_at <- as.integer(opt$n_at %||% 50L)
probs <- c(0.025, 0.975) # outer 95% CI
probs50 <- c(0.25, 0.75) # inner 50% CrI

# 7. Per-variable loop ---------------------------------------------------------
# Figures per variable (figures/cb/<run_tag>/country/):
#   <var>_lag_response_country  — faceted lag-response at median exposure
#   <var>_cumulative_country    — faceted cumulative exposure-response
#   <var>_heatmap_country       — faceted exposure x lag surface
#   <var>_erf_country_sl<N>     — single-lag fits: faceted ERF (no lag axis)

for (v in vars_req) {
  vm <- cb_meta[[v]]
  cb_col_names <- vm$cb_col_names
  if (is.null(cb_col_names) || length(cb_col_names) == 0) {
    next
  }
  n_cols <- length(cb_col_names)

  # Check that this variable actually has country-level RE slopes in the fit.
  # brms replaces spaces with dots in posterior parameter names, so we must
  # apply the same substitution when building lookup strings.
  country_pname <- function(lvl) gsub(" ", ".", lvl, fixed = TRUE)
  test_re_name <- sprintf(
    "r_%s[%s,%s]",
    country_group,
    country_pname(country_levels[1]),
    cb_col_names[1]
  )
  re_check <- tryCatch(
    brms::as_draws_matrix(brms_fit, variable = test_re_name),
    error = function(e) NULL
  )
  if (is.null(re_check) || nrow(re_check) == 0) {
    message(sprintf(
      "[%s] No country-level RE slopes found — skipping per-country figures.",
      v
    ))
    next
  }

  message(sprintf("[%s] extracting per-country surfaces ...", v))

  # Global FE draws: n_draws × n_cols
  b_names <- sprintf("b_%s", cb_col_names)
  fe_draws <- unclass(brms::as_draws_matrix(brms_fit, variable = b_names))
  n_draws <- nrow(fe_draws)

  # Hierarchical RE draws for one country level: sum the region + country
  # slope deviations (each n_draws × n_cols). Missing parent group -> 0.
  # Region prefix stripped from a nested "<region>_<country>" level string
  # (ungrouped runs only — grouping-mode "country" levels are bare ISO3 and
  # need no stripping, str_remove is then a no-op).
  grab_re <- function(group, level) {
    vnames <- sprintf("r_%s[%s,%s]", group, country_pname(level), cb_col_names)
    tryCatch(
      unclass(brms::as_draws_matrix(brms_fit, variable = vnames)),
      error = function(e) matrix(0, n_draws, n_cols)
    )
  }
  nested_re_draws <- function(country_level) {
    acc <- grab_re(country_group, country_level)
    if (!is.na(region_group) && region_group %in% re_groups) {
      region_level <- stringr::str_remove(country_level, "_[A-Z]{3}$")
      acc <- acc + grab_re(region_group, region_level)
    }
    acc
  }

  # Prediction grid
  rng <- vm$range %||% c(NA_real_, NA_real_)
  at_vals <- seq(rng[1], rng[2], length.out = n_at)
  mid_idx <- which.min(abs(at_vals - median(at_vals)))

  # ---- Single-lag fit: no lag axis → faceted per-country exposure-response ----
  if (!is.null(vm$single_lag)) {
    L <- as.integer(vm$single_lag)
    # onebasis stub: reattach stored attributes so crosspred can reconstruct it
    ob_stub <- structure(
      matrix(0, nrow = 1L, ncol = n_cols),
      class = c("onebasis", "matrix")
    )
    for (a in names(vm$onebasis_attrs)) {
      attr(ob_stub, a) <- vm$onebasis_attrs[[a]]
    }

    erf_rows <- list()
    for (ctry in country_levels) {
      re_draws <- nested_re_draws(ctry)
      combined <- fe_draws + re_draws
      draw_preds <- vector("list", n_draws)
      for (i in seq_len(n_draws)) {
        draw_preds[[i]] <- tryCatch(
          dlnm::crosspred(
            ob_stub,
            coef = as.numeric(combined[i, ]),
            vcov = diag(0, n_cols),
            at = at_vals,
            cen = 0
          ),
          error = function(e) NULL
        )
      }
      draw_preds <- Filter(Negate(is.null), draw_preds)
      if (length(draw_preds) == 0) {
        warning(sprintf("All draws failed for country='%s' — skipping.", ctry))
        next
      }
      erf_mat <- do.call(
        cbind,
        lapply(draw_preds, function(x) as.numeric(x$matfit))
      )
      erf_rows[[ctry]] <- data.frame(
        country = ctry,
        exposure = at_vals,
        logRR = apply(erf_mat, 1, median),
        lo = apply(erf_mat, 1, quantile, probs = probs[1]),
        hi = apply(erf_mat, 1, quantile, probs = probs[2]),
        lo50 = apply(erf_mat, 1, quantile, probs = probs50[1]),
        hi50 = apply(erf_mat, 1, quantile, probs = probs50[2])
      )
      message(sprintf(
        "  %s: %d/%d successful draws",
        ctry,
        length(draw_preds),
        n_draws
      ))
    }
    if (length(erf_rows) == 0) {
      warning(sprintf("[%s] No country ERFs extracted.", v))
      next
    }
    erf_df <- dplyr::bind_rows(erf_rows)
    p_erf <- ggplot(erf_df, aes(x = exposure)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
      geom_ribbon(
        aes(ymin = lo50, ymax = hi50),
        alpha = 0.35,
        fill = "#2166ac"
      ) +
      geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      facet_wrap(~country) +
      labs(
        title = sprintf(
          "Exposure-response per country (single lag %d): %s",
          L,
          v
        ),
        subtitle = "Posterior median + 50% / 95% CrI",
        x = v,
        y = "log-RR"
      ) +
      theme_bw(base_size = 10)
    save_fig(p_erf, sprintf("%s_erf_country_%s", v, sl_tag(L)))
    next
  }

  # CB stub
  cb_stub <- structure(
    matrix(0, nrow = 1L, ncol = n_cols),
    class = c("crossbasis", "matrix"),
    argvar = vm$argvar,
    arglag = vm$arglag,
    lag = vm$lag,
    df = vm$df,
    range = vm$range
  )

  lag_rows <- list()
  cum_rows <- list()
  heat_rows <- list()

  for (ctry in country_levels) {
    # Hierarchical country RE draws: FE + region + country (n_draws × n_cols)
    re_draws <- nested_re_draws(ctry)

    # Combine FE + RE
    combined <- fe_draws + re_draws

    # Per-draw crosspred
    draw_preds <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      draw_preds[[i]] <- tryCatch(
        dlnm::crosspred(
          cb_stub,
          coef = as.numeric(combined[i, ]),
          vcov = diag(0, n_cols),
          at = at_vals,
          cen = 0,
          cumul = TRUE
        ),
        error = function(e) NULL
      )
    }
    draw_preds <- Filter(Negate(is.null), draw_preds)
    if (length(draw_preds) == 0) {
      warning(sprintf("All draws failed for country='%s' — skipping.", ctry))
      next
    }

    # Shared quantities across the three figure types
    mat_list <- lapply(draw_preds, function(x) x$matfit)
    mat_arr <- simplify2array(mat_list)
    lag_seq <- as.numeric(sub("^lag", "", colnames(draw_preds[[1]]$matfit)))
    mat_med <- apply(mat_arr, c(1, 2), median) # at × lag posterior median surface

    # Lag-response at median exposure
    lag_rows[[ctry]] <- data.frame(
      country = ctry,
      lag = lag_seq,
      logRR = apply(mat_arr, c(1, 2), median)[mid_idx, ],
      lo = apply(mat_arr, c(1, 2), quantile, probs = probs[1])[mid_idx, ],
      hi = apply(mat_arr, c(1, 2), quantile, probs = probs[2])[mid_idx, ],
      lo50 = apply(mat_arr, c(1, 2), quantile, probs = probs50[1])[mid_idx, ],
      hi50 = apply(mat_arr, c(1, 2), quantile, probs = probs50[2])[mid_idx, ]
    )

    # Cumulative exposure-response
    cum_list <- lapply(draw_preds, function(x) x$allfit)
    cum_arr <- do.call(cbind, lapply(cum_list, as.numeric))
    cum_rows[[ctry]] <- data.frame(
      country = ctry,
      exposure = at_vals,
      cumRR = apply(cum_arr, 1, median),
      lo = apply(cum_arr, 1, quantile, probs = probs[1]),
      hi = apply(cum_arr, 1, quantile, probs = probs[2]),
      lo50 = apply(cum_arr, 1, quantile, probs = probs50[1]),
      hi50 = apply(cum_arr, 1, quantile, probs = probs50[2])
    )

    # Heatmap surface (full at × lag grid)
    colnames(mat_med) <- lag_seq
    rownames(mat_med) <- round(at_vals, 3)
    heat_rows[[ctry]] <- as.data.frame(mat_med) %>%
      tibble::rownames_to_column("exposure") %>%
      tidyr::pivot_longer(-exposure, names_to = "lag", values_to = "logRR") %>%
      dplyr::mutate(
        country = ctry,
        exposure = as.numeric(exposure),
        lag = as.numeric(lag)
      )

    message(sprintf(
      "  %s: %d/%d successful draws",
      ctry,
      length(draw_preds),
      n_draws
    ))
  }

  if (length(lag_rows) == 0) {
    warning(sprintf("[%s] No country surfaces extracted.", v))
    next
  }

  lag_df <- dplyr::bind_rows(lag_rows)
  cum_df <- dplyr::bind_rows(cum_rows)
  heat_df <- dplyr::bind_rows(heat_rows)
  med_exp <- round(at_vals[mid_idx], 2)

  # Shared colour scale limits across all country panels (symmetric around 0)
  heat_lim <- max(abs(heat_df$logRR), na.rm = TRUE)

  # ---- Faceted lag-response ----
  p_lag <- ggplot(lag_df, aes(x = lag)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#2166ac") +
    geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~country) +
    labs(
      title = sprintf("Lag-response per country: %s", v),
      subtitle = sprintf(
        "At median exposure (%s = %.2f) — posterior median + 50%% / 95%% CrI",
        v,
        med_exp
      ),
      x = "Lag (months)",
      y = "log-RR"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_lag, sprintf("%s_lag_response_country", v))

  # ---- Faceted cumulative exposure-response ----
  p_cum <- ggplot(cum_df, aes(x = exposure)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#d73027") +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#d73027") +
    geom_line(aes(y = cumRR), colour = "#d73027", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~country) +
    labs(
      title = sprintf("Cumulative exposure-response per country: %s", v),
      subtitle = "Summed over all lags — posterior median + 50% / 95% CrI",
      x = v,
      y = "Cumulative log-RR"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_cum, sprintf("%s_cumulative_country", v))

  # ---- Faceted heatmap (exposure × lag surface per country) ----
  p_heat <- ggplot(heat_df, aes(x = exposure, y = lag, fill = logRR)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      limits = c(-heat_lim, heat_lim),
      name = "log-RR"
    ) +
    facet_wrap(~country) +
    labs(
      title = sprintf("Exposure-lag surface per country: %s", v),
      subtitle = "Posterior median log-RR",
      x = v,
      y = "Lag (months)"
    ) +
    theme_bw(base_size = 10)
  save_fig(p_heat, sprintf("%s_heatmap_country", v), width = 13, height = 9)
}

message(sprintf("Done. Per-country figures written to %s/", fig_dir))
