# 06_fit_figures_cb_grouped.R — ERF / IRR figures for one grouped ENSO single-lag fit
#
# Sibling to 06_fit_figures_cb.R / _cb_adm1.R / _cb_country.R. The extraction,
# map and lag-selection logic lives in shared utils modules
# (utils/figures_cb.R, utils/loo_lags.R) that 09b_postprocess_results.R also uses.
#
# For a grouped (grouping.enable) single-lag DLNM fit of ANY track (occurrence |
# size | duration | monthly), the ENSO-slope RE nesting is auto-detected from the
# ranef group names:
#   region kept    : AFRO_region / AFRO_region:country / AFRO_region:country:ADM1
#   region dropped : country / country:ADM1   (grouped prod fits; region term is 0)
#
# Figures (figures/cb/<config-stem>/grouped/):
#   03_country_erf.png      — country ERF (b_pop + region + country slope RE), 95% CI
#   04_adm1_irr_map.png     — ADM1 IRR map at the four composite ENSO conditions
#   05_country_irr_map.png  — country median IRR (across ADM1 units), same conditions
#
# Two ways to name the fit:
#   -C <config.yml>                       exactly that config
#   -d <batch_dir> --group <enso_group>   the group's BEST lag in the batch, resolved
#       like the deployed fits and the manuscript: latest loo_compare_lags ranking
#       -> select_best_lag() (utils/loo_lags.R::best_lags_by_group()), then the
#       <batch>_sl-<lag>_grp-<group>.yml config (utils/paths.R::groupedLagConfigPath()).
#       The track is detected from the batch's configs.
#
# Usage:
#   Rscript analysis/06_fit_figures_cb_grouped.R \
#     -C analysis/configs/prod_run_occ_grouped_enso_singlelag/<config>.yml
#   Rscript analysis/06_fit_figures_cb_grouped.R \
#     -d analysis/configs/prod_run_occ_grouped_enso_singlelag --group enso_east
#   ... --format pdf --dpi 300 --n_draws 200 --fig_dir figures/my_dir

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(brms)
})

source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/loo_lags.R")
source("analysis/utils/figures_cb.R")
source("analysis/utils/figures_io.R")

option_list <- list(
  make_option(c("-d", "--batch_dir"), default = NULL, type = "character",
              help = "Grouped single-lag batch config dir (with --group): use the group's best lag"),
  make_option(c("--group"), default = NULL, type = "character",
              help = "ENSO group token for -d, e.g. enso_east"),
  make_option(c("--lag_rule"), default = "longest_within_se", type = "character",
              help = "With -d: best-lag rule, longest_within_se | best_elpd [default: longest_within_se]"),
  make_option(c("--format"), default = "png", type = "character",
              help = "Output format: pdf | png [default: png]"),
  make_option(c("--dpi"), default = 200L, type = "integer",
              help = "DPI for raster output [default: 200]"),
  make_option(c("--n_draws"), default = 400L, type = "integer",
              help = "Posterior draws for the CIs [default: 400]"),
  # NOT --outdir: that is the shared artefact-root key makeRunFile() reads.
  make_option(c("--fig_dir"), default = NULL, type = "character",
              help = "Figure directory [default: figures/cb/<config-stem>/grouped]")
)

# Resolve the config: -C wins; else -d/--group -> best lag -> config path ---------
argv     <- commandArgs(trailingOnly = TRUE)
if (any(c("-h", "--help") %in% argv)) {
  print_help(OptionParser(option_list = option_list,
                          usage = "Rscript analysis/06_fit_figures_cb_grouped.R (-C <config.yml> | -d <batch_dir> --group <enso_group>) [options]"))
  quit(status = 0)
}
cfg_path <- .flag_val("-C", argv) %||% .flag_val("--config", argv)
if (is.null(cfg_path)) {
  batch_dir <- .flag_val("-d", argv) %||% .flag_val("--batch_dir", argv)
  group     <- .flag_val("--group", argv)
  if (is.null(batch_dir) || is.null(group))
    stop("Supply -C <config.yml>, or -d <batch_dir> together with --group <enso_group>.")
  cfgs <- list.files(batch_dir, pattern = "\\.yml$", full.names = TRUE)
  cfgs <- cfgs[basename(cfgs) != "manifest.yml"]
  if (!length(cfgs)) stop("No configs in ", batch_dir)
  track <- detect_model_from_yaml(cfgs[1L])
  ls_rule <- lag_selection_cfg(list(rule = .flag_val("--lag_rule", argv) %||% "longest_within_se"))
  lag   <- best_lags_by_group(batch_dir, group, track, min_lag = ls_rule$min_lag, rule = ls_rule$rule)
  if (!length(lag)) stop("No LOO ranking for group '", group, "' in ", batch_dir,
                         " -- run 05_fit_postprocess_loo_lags.R first.")
  cfg_path <- groupedLagConfigPath(batch_dir, group, lag[[group]])
  if (is.null(cfg_path)) stop("No config for group ", group, " at lag ", lag[[group]])
  message(sprintf("Best lag (%s, %s): %d -> %s", track, group, lag[[group]], cfg_path))
}

# Drop the flags handled above so make_options_from_config() sees -C only.
.drop_flags <- function(argv, flags) {
  keep <- rep(TRUE, length(argv))
  for (i in seq_along(argv)) if (argv[i] %in% flags) keep[i:min(i + 1L, length(argv))] <- FALSE
  argv[keep]
}
cli_argv <- c("-C", cfg_path, .drop_flags(argv, c("-C", "--config", "-d", "--batch_dir", "--group", "--lag_rule")))

model_type <- detect_model_from_yaml(cfg_path)
opt <- make_options_from_config(option_list = option_list, cli_argv = cli_argv,
                                model = model_type)
stopifnot(opt[[model_type]]$backend == "brms_cb",
          isTRUE(opt$grouping$enable), !is.null(opt$grouping$group))
grp_label   <- opt$grouping$group
model_label <- c(occurrence = "outbreak occurrence", size = "outbreak size",
                 duration = "outbreak cessation", monthly = "monthly incidence")[[model_type]]
message("Config: ", cfg_path, " | model: ", model_type, " | group: ", grp_label)

# Fit + CB metadata --------------------------------------------------------------
kpfx      <- c(occurrence = "", size = "size_", duration = "dur_", monthly = "mon_")[[model_type]]
fit_file  <- makeRunFile(opt, kind = paste0(kpfx, "fit"))
meta_file <- makeRunFile(opt, kind = paste0(kpfx, "cb_meta"))
if (!file.exists(fit_file))  stop("Fit not found: ", fit_file, "\n  Run 04_fit_run.R with backend: brms_cb first.")
if (!file.exists(meta_file)) stop("CB meta not found: ", meta_file)
message("Loading fit: ", fit_file)
brms_fit <- readRDS(fit_file)
meta     <- readRDS(meta_file)

vm <- meta[["enso"]]
if (is.null(vm$single_lag)) stop("vm$single_lag is NULL -- this script is for single-lag fits only.")
lag_val <- vm$single_lag
message(sprintf("Single-lag fit: lag = %d", lag_val))

# Output directory ----------------------------------------------------------------
fmt     <- match.arg(opt$format %||% "png", c("pdf", "png"))
fig_dir <- opt$fig_dir %||% makeFigureDir(opt, track = "cb", config_path = cfg_path, subdir = "grouped")
save_fig <- figureSaver(fig_dir, width = 18, height = 8, dpi = as.integer(opt$dpi), formats = fmt)
message("Figures -> ", fig_dir)

n_draws <- as.integer(opt$n_draws %||% 400L)
seed    <- opt$seed %||% 1234L
grid    <- build_onebasis_grid(vm, n_at = 100L, map_vals = COMPOSITE_CONDS)
geo     <- load_africa_sf()

# Section A: country-level ERF (03_country_erf) ------------------------------------
message("\n=== Section A: country-level ERF ===")
erf_country <- extract_country_erf(brms_fit, meta, grid, group_label = grp_label,
                                   n_draws = n_draws, seed = seed)
if (is.null(erf_country) || !nrow(erf_country)) {
  message("Skipping country ERF (country RE not found).")
} else {
  region_lookup <- taxdat::get_AFRO_region(dplyr::distinct(erf_country, country), ctry_col = "country")
  erf_country <- erf_country %>%
    dplyr::left_join(region_lookup, by = "country") %>%
    dplyr::mutate(
      AFRO_region = factor(AFRO_region, levels = sort(unique(AFRO_region))),
      country     = factor(country, levels = {
        tmp <- dplyr::distinct(., AFRO_region, country)
        as.character(tmp[order(tmp$AFRO_region, tmp$country), ]$country)
      }))
  n_countries <- nlevels(erf_country$country)
  p_erf <- ggplot(erf_country, aes(x = at)) +
    geom_hline(yintercept = 1, lty = 2, lwd = 0.3, color = "gray40") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.3, color = "gray40") +
    geom_ribbon(aes(ymin = IRR_lo, ymax = IRR_hi, fill = AFRO_region), alpha = 0.20, colour = NA) +
    geom_line(aes(y = IRR, color = AFRO_region), lwd = 0.7) +
    facet_wrap(~ country, ncol = 4) +
    scale_y_log10() +
    scale_color_manual(values = taxdat::colors_afro_regions(), name = "AFRO region") +
    scale_fill_manual(values = taxdat::colors_afro_regions(), guide = "none") +
    labs(x = "ENSO ONI", y = "IRR (log scale)",
         title = sprintf("Country-level ENSO ERF — %s (%s, lag %d)", model_label, grp_label, lag_val),
         subtitle = paste("b_pop + region + country slope RE; shaded = 95% CI;",
                          "colored by AFRO region; centered at ENSO = 0")) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_fig(p_erf, "03_country_erf", width = 12,
           height = ceiling(n_countries / 4L) * 3.0 + 2.0, limitsize = FALSE)
}

# Section B: ADM1 IRR map (04_adm1_irr_map) -----------------------------------------
message("\n=== Section B: ADM1-level ENSO IRR map ===")
adm1_irr <- extract_adm1_irr(brms_fit, meta, enso_conds = COMPOSITE_CONDS, group_label = grp_label,
                             n_draws = n_draws, seed = seed, n_at = 100L)
if (is.null(adm1_irr) || !nrow(adm1_irr)) {
  message("Skipping ADM1 map (ADM1 RE not found).")
} else {
  message(sprintf("ADM1 units: %d | IRR range [%.2f, %.2f] | significant per condition: %s",
                  dplyr::n_distinct(adm1_irr$adm1_label), min(adm1_irr$IRR), max(adm1_irr$IRR),
                  paste(sprintf("%s=%d", names(COMPOSITE_CONDS),
                                tapply(adm1_irr$significant, factor(adm1_irr$enso_label,
                                                                    levels = names(COMPOSITE_CONDS)), sum)),
                        collapse = "; ")))
  p_adm1 <- plot_irr_map(
    adm1_irr, geo, level = "adm1", enso_levels = names(COMPOSITE_CONDS),
    fill_name = sprintf("IRR\nat lag %d\n(capped)", lag_val),
    title = sprintf("ADM1-level ENSO IRR at lag %d — %s (%s)", lag_val, model_label, grp_label),
    subtitle = paste("IRR = exp(pop + region + country + ADM1 RE);",
                     "outlined = 95% CI excludes 1; blue = protective, red = risk"))
  if (!is.null(p_adm1)) save_fig(p_adm1, "04_adm1_irr_map")

  # Section C: country median IRR map (05_country_irr_map) ----------------------------
  message("\n=== Section C: country-level IRR map ===")
  p_ctry <- plot_irr_map(
    country_irr_from_adm1(adm1_irr), geo, level = "country", value_col = "median_IRR",
    enso_levels = names(COMPOSITE_CONDS), outline_sig = FALSE,
    fill_name = sprintf("Median\nIRR\nat lag %d\n(capped)", lag_val),
    title = sprintf("Country-level ENSO IRR at lag %d — %s (%s)", lag_val, model_label, grp_label),
    subtitle = "Median IRR across ADM1 units; gray = not in model group")
  if (!is.null(p_ctry)) save_fig(p_ctry, "05_country_irr_map")
}

message("\nAll figures written to: ", fig_dir)
