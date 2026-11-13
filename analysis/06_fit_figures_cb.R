# 06_fit_figures_cb.R — crossbasis exposure-response surface figures
#
# Sibling to 06_fit_figures.R (unchanged by side-by-side rule).
# Reads the cb_surfaces output from 05_fit_postprocess_cb.R and generates:
#
#   Per variable:
#     1. Heatmap of the log-RR surface over (exposure × lag)
#     2. Cumulative exposure-response curve (summed over all lags, with CI)
#     3. Lag-response curve at the median exposure value (with CI)
#
# Output: figures/cb/<run_tag>/<var>_surface.pdf  (or .png, --format)
#
# Usage:
#   Rscript analysis/06_fit_figures_cb.R -C analysis/configs/run_cb_precip/run_cb_precip.yml
#   Rscript analysis/06_fit_figures_cb.R -C ... --format png --dpi 200

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(optparse)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/figures_io.R")

option_list <- list(
  make_option(c("-C", "--config"),
              default = NULL, type = "character",
              help = "Single-lag crossbasis fit config YAML (required)"),
  make_option(c("-c", "--country"),
              default = NULL, type = "character",
              help = "ISO3 code(s) or 'all' — overrides config country [default: from config]"),
  make_option(c("--format"),
              default = "pdf", type = "character",
              help = "Output format: pdf | png [default: pdf]"),
  make_option(c("--dpi"),
              default = 150L, type = "integer",
              help = "DPI for raster output [default: 150]"),
  make_option(c("--vars"),
              default = NULL, type = "character",
              help = "Comma-separated variable names to plot (default: all)")
)

# Auto-detect occurrence vs size from the config so size single-lag fits get
# figures too.
.detect_model <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  cfg  <- NULL
  op_cfg <- parse_args(OptionParser(option_list = option_list))$config
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) { cfg <- argv[i + 1L]; break }
  } 
  
  if (is.null(cfg)) {
    if (!is.null(op_cfg)) {
      cfg <- op_cfg
    }
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

# File-kind prefix: size track uses size_cb_surfaces etc; duration uses
# dur_cb_surfaces; monthly uses mon_cb_surfaces.
kpfx <- if (model_type == "size") "size_" else if (model_type == "duration") "dur_" else if (model_type == "monthly") "mon_" else ""

# Locate files ------------------------------------------------------------

surf_file <- makeRunFile(opt, kind = paste0(kpfx, "cb_surfaces"))
if (!file.exists(surf_file))
  stop("Surfaces file not found: ", surf_file,
       "\n  Run 05_fit_postprocess_cb.R first.")

message("Loading surfaces: ", surf_file)
surfaces <- readRDS(surf_file)

if (length(surfaces) == 0)
  stop("No surfaces in file: ", surf_file)

vars <- if (!is.null(opt$vars)) {
  intersect(trimws(strsplit(opt$vars, ",")[[1]]), names(surfaces))
} else {
  names(surfaces)
}
if (length(vars) == 0) stop("No matching variables to plot.")

# Output directory --------------------------------------------------------

fmt <- match.arg(opt$format %||% "pdf", c("pdf", "png"))
# makeFigureDir() (utils/paths.R) recovers -C/--config from argv
# itself (derive_config_run_tag() -> .flag_val()) when config_path is
# omitted, reproducing the .cfg_path scan this replaces byte-for-byte.
fig_dir <- makeFigureDir(opt = list(), track = "cb")

# Helper: save ggplot (figureSaver(), utils/figures_io.R)
save_fig <- figureSaver(fig_dir, dpi = as.integer(opt$dpi %||% 150L),
                        formats = fmt)

# Figures -----------------------------------------------------------------

for (v in vars) {
  surf <- surfaces[[v]]
  if (is.null(surf)) next
  
  at_vals  <- surf$at
  matfit   <- surf$matfit   # length(at) × n_lags, median log-RR
  matlo    <- surf$matlo     # 95% CI
  mathi    <- surf$mathi
  matlo50  <- surf$matlo50   # 50% CrI
  mathi50  <- surf$mathi50
  
  # ---- Single-lag fit: no lag axis → plot the exposure-response curve only ----
  if (!is.null(surf$single_lag)) {
    L <- as.integer(surf$single_lag)
    erf_df <- data.frame(
      exposure = at_vals,
      logRR    = matfit[, 1],
      lo       = matlo[, 1],
      hi       = mathi[, 1],
      lo50     = matlo50[, 1],
      hi50     = mathi50[, 1]
    )
    p_erf <- ggplot(erf_df, aes(x = exposure)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
      geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#2166ac") +
      geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(
        title    = sprintf("Exposure-response (single lag %d): %s", L, v),
        subtitle = "Posterior median + 50% / 95% CrI",
        x        = v,
        y        = "log-RR"
      ) +
      theme_bw(base_size = 11)
    save_fig(p_erf, sprintf("%s_erf_%s", v, slTag(L)), width = 6, height = 4)
    next
  }
  
  lag_axis <- as.numeric(surf$lag)
  
  # ---- 1. Surface heatmap (log-RR over exposure × lag) ----
  # Reshape matfit to long format
  colnames(matfit) <- lag_axis
  rownames(matfit) <- round(at_vals, 3)
  
  heat_df <- as.data.frame(matfit) %>%
    tibble::rownames_to_column("exposure") %>%
    tidyr::pivot_longer(-exposure, names_to = "lag", values_to = "logRR") %>%
    dplyr::mutate(
      exposure = as.numeric(exposure),
      lag      = as.numeric(lag)
    )
  
  p_heat <- ggplot(heat_df, aes(x = exposure, y = lag, fill = logRR)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d73027",
      midpoint = 0,
      name     = "log-RR"
    ) +
    labs(
      title    = sprintf("Exposure-lag surface: %s", v),
      subtitle = "Posterior median log-RR",
      x        = v,
      y        = "Lag (months)"
    ) +
    theme_bw(base_size = 11)
  save_fig(p_heat, sprintf("%s_heatmap", v), width = 7, height = 4.5)
  
  # ---- 2. Lag-response at median exposure ----
  mid_idx  <- which.min(abs(at_vals - median(at_vals)))
  lag_df   <- data.frame(
    lag    = lag_axis,
    logRR  = matfit[mid_idx, ],
    lo     = matlo[mid_idx, ],
    hi     = mathi[mid_idx, ],
    lo50   = matlo50[mid_idx, ],
    hi50   = mathi50[mid_idx, ]
  )

  p_lag <- ggplot(lag_df, aes(x = lag)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#2166ac") +
    geom_line(aes(y = logRR), colour = "#2166ac", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(
      title    = sprintf("Lag-response at median exposure: %s = %.2f",
                         v, at_vals[mid_idx]),
      subtitle = "Posterior median + 50% / 95% CrI",
      x        = "Lag (months)",
      y        = "log-RR"
    ) +
    theme_bw(base_size = 11)
  save_fig(p_lag, sprintf("%s_lag_response", v), width = 6, height = 4)
  
  # ---- 3. Cumulative exposure-response (if computed) ----
  if (!is.null(surf$cumfit)) {
    cum_df <- data.frame(
      exposure = at_vals,
      cumRR    = surf$cumfit,
      lo       = surf$cumlo,
      hi       = surf$cumhi,
      lo50     = surf$cumlo50,
      hi50     = surf$cumhi50
    )

    p_cum <- ggplot(cum_df, aes(x = exposure)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#d73027") +
      geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.35, fill = "#d73027") +
      geom_line(aes(y = cumRR), colour = "#d73027", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(
        title    = sprintf("Cumulative exposure-response: %s", v),
        subtitle = "Summed over all lags — posterior median + 50% / 95% CrI",
        x        = v,
        y        = "Cumulative log-RR"
      ) +
      theme_bw(base_size = 11)
    save_fig(p_cum, sprintf("%s_cumulative", v), width = 6, height = 4)
  }
}

message(sprintf("Done. Figures written to %s/", fig_dir))
