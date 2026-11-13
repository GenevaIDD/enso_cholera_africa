# 06_fit_figures_cb_grouped_composite.R — composite whole-Africa ENSO IRR maps
#
# Stitches the best-lag grouped single-lag ENSO fits of every group into one
# Africa-wide map per track. Sibling to 06_fit_figures_cb_grouped.R. The IRR
# extraction, map drawing and lag selection live in the shared utils modules
# (utils/figures_cb.R, utils/loo_lags.R) that 09b_postprocess_results.R also uses.
#
# Best lag per group is the CANONICAL rule (utils/loo_lags.R::best_lags_by_group():
# latest loo_compare_lags ranking -> select_best_lag(), the longest lag within
# --se_mult SEs of the best), i.e. the same pick as the deployed fits and the
# manuscript. --se_mult 0 takes the raw rank-1 lag.
# Groups with no ranking (and no override) are dropped; enso_central is not a
# candidate (no countries with data).
#
# Figures (figures/cb/<batch-dir-name>/composite/):
#   01_enso_grouping_map.png            — country fill = ENSO teleconnection group
#   02_occ_composite_irr.png            — ADM1 IRR at the four composite conditions (occurrence)
#   03_siz_composite_irr.png            — ADM1 IRR (size)
#   04_occ_composite_country_irr.png    — country IRR (occurrence)
#   05_siz_composite_country_irr.png    — country IRR (size)
#
# Usage:
#   Rscript analysis/06_fit_figures_cb_grouped_composite.R \
#     --occ_dir analysis/configs/prod_run_occ_grouped_enso_singlelag \
#     --siz_dir analysis/configs/<size batch> \
#     --grouping analysis/groupings/enso_grouping_prod_run.csv
#   # override specific groups (the rest are resolved from the rankings):
#   #   --occ_lags "enso_east=2"  --siz_lags "enso_horn=4,enso_west=5"
#   # raw rank-1 lag:  --se_mult 0

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(brms)
  library(sf)
})

# results_core.R sources base/config/paths/provenance/enso/gid/figures_cb itself.
source("analysis/utils/results_core.R")
source("analysis/utils/loo_lags.R")
source("analysis/utils/figures_io.R")

CANDIDATE_GROUPS <- c("enso_east", "enso_south", "enso_west", "enso_horn")

option_list <- list(
  make_option("--occ_dir", type = "character", default = NULL,
              help = "Occurrence grouped batch config directory"),
  make_option("--siz_dir", type = "character", default = NULL,
              help = "Size grouped batch config directory"),
  make_option("--grouping", type = "character", default = "analysis/groupings/enso_grouping.csv",
              help = "enso_grouping.csv [default: %default]"),
  make_option("--occ_lags", type = "character", default = NULL,
              help = "Occurrence lag overrides, 'group=lag,...' (e.g. 'enso_east=2')"),
  make_option("--siz_lags", type = "character", default = NULL,
              help = "Size lag overrides, 'group=lag,...' (e.g. 'enso_horn=4')"),
  make_option("--se_mult", type = "double", default = 2,
              help = "SE multiplier of the best-lag rule; 0 = raw rank-1 lag [default: %default]"),
  make_option("--lag_rule", type = "character", default = "longest_within_se",
              help = "Best-lag rule: longest_within_se | best_elpd (deploy.lag_selection.rule) [default: %default]"),
  make_option("--occ_floor_offset", type = "double", default = NA,
              help = "best_elpd only: size pick >= occurrence lag minus this (NA = no floor) [default: %default]"),
  make_option("--enforce_occ_floor", type = "logical", default = FALSE,
              help = "Floor the size lag at the group's occurrence lag [default: %default]"),
  make_option("--n_draws", type = "integer", default = 400L,
              help = "Posterior draws for the CIs [default: %default]"),
  make_option("--dpi", type = "integer", default = 200L, help = "Raster DPI [default: %default]"),
  make_option("--fig_dir", type = "character", default = NULL,
              help = "Figure directory [default: figures/cb/<batch-dir-name>/composite]")
)
cli <- parse_args(OptionParser(option_list = option_list))
if (is.null(cli$occ_dir) && is.null(cli$siz_dir))
  stop("Supply at least one of --occ_dir / --siz_dir.")

fig_dir  <- cli$fig_dir %||% makeFigureDir(list(), track = "cb",
                                           config_path = cli$occ_dir %||% cli$siz_dir,
                                           subdir = "composite")
save_fig <- figureSaver(fig_dir, width = 14, height = 8, dpi = cli$dpi)
message("Figures -> ", fig_dir)

# 'group=lag,group=lag' -> named integer vector
parse_lag_overrides <- function(s) {
  if (is.null(s) || !nzchar(s)) return(stats::setNames(integer(0), character(0)))
  kv  <- strsplit(trimws(strsplit(s, ",", fixed = TRUE)[[1L]]), "=", fixed = TRUE)
  if (any(vapply(kv, length, integer(1)) != 2L)) stop("Malformed lag override (expected 'group=lag'): ", s)
  stats::setNames(vapply(kv, function(x) as.integer(trimws(x[2L])), integer(1)),
                  vapply(kv, function(x) trimws(x[1L]), character(1)))
}

# Best lags (occurrence first: the size pick may be floored at it) ------------------
resolve_lags <- function(dir, track, overrides, occ_lags = NULL) {
  if (is.null(dir)) return(stats::setNames(integer(0), character(0)))
  ls_rule <- lag_selection_cfg(list(rule = cli$lag_rule, se_mult = cli$se_mult,
                                    occ_floor_offset = if (is.na(cli$occ_floor_offset)) "none" else cli$occ_floor_offset))
  best_lags_by_group(dir, CANDIDATE_GROUPS, track, occ_lags = occ_lags, se_mult = cli$se_mult,
                     enforce_occ_floor = isTRUE(cli$enforce_occ_floor) && !is.null(occ_lags),
                     overrides = overrides, min_lag = ls_rule$min_lag, rule = ls_rule$rule,
                     occ_floor_offset = ls_rule$occ_floor_offset)
}
OCC_LAGS <- resolve_lags(cli$occ_dir, "occurrence", parse_lag_overrides(cli$occ_lags))
SIZ_LAGS <- resolve_lags(cli$siz_dir, "size", parse_lag_overrides(cli$siz_lags), occ_lags = OCC_LAGS)
for (nm in c("OCC", "SIZ")) {
  lg <- get(paste0(nm, "_LAGS"))
  if (length(lg)) message(nm, " best lags: ", paste(sprintf("%s=%d", names(lg), lg), collapse = ", "))
}

# One fit load per (track, group): ADM1 and country IRR tables together -------------
load_group_tables <- function(batch_dir, group, lag, track) {
  cfg <- groupedLagConfigPath(batch_dir, group, lag)
  if (is.null(cfg)) return(NULL)
  opt  <- suppressMessages(make_options_from_config(option_list = list(),
                                                    cli_argv = c("-C", cfg), model = track))
  kpfx <- c(occurrence = "", size = "size_")[[track]]
  fit_file  <- makeRunFile(opt, kind = paste0(kpfx, "fit"))
  meta_file <- makeRunFile(opt, kind = paste0(kpfx, "cb_meta"))
  if (!file.exists(fit_file) || !file.exists(meta_file)) {
    warning(sprintf("Fit/meta not found (%s, %s, lag %d): %s", track, group, lag, fit_file))
    return(NULL)
  }
  message(sprintf("  Loading %s | %s | lag %d ...", track, group, lag))
  fit <- readRDS(fit_file); meta <- readRDS(meta_file)
  out <- list(adm1    = extract_adm1_irr(fit, meta, group_label = group, n_draws = cli$n_draws),
              country = extract_country_irr(fit, meta, group_label = group, n_draws = cli$n_draws))
  rm(fit); gc(verbose = FALSE)
  out
}

# Geometries ------------------------------------------------------------------------
message("\n=== Loading Africa geometries ===")
geo <- load_africa_sf()
if (is.null(geo$adm1_sf)) stop("ADM1 shapefile required for the composite map -- check gadm_410-levels.gpkg.")

# Section 0: ENSO grouping membership choropleth ------------------------------------
message("\n=== Section 0: ENSO grouping map ===")
enso_grp   <- load_enso_grouping(cli$grouping)
# Labels, order and palette from the registry's own sidecar
# (analysis/groupings/<source>.yml; utils/config.R::load_grouping_meta()).
grp_meta   <- load_grouping_meta(sub("\\.csv$", "", basename(cli$grouping)),
                                 groupings_dir = dirname(cli$grouping))
grp_labels <- stats::setNames(paste0(grp_meta$short[grp_meta$order], " (", grp_meta$order, ")"),
                              grp_meta$order)
grp_sf <- geo$afro_sf %>%
  dplyr::left_join(enso_grp %>% dplyr::mutate(group_label = grp_labels[group]),
                   by = c("ISO3_CODE" = "iso3")) %>%
  dplyr::mutate(group_label = factor(tidyr::replace_na(group_label, "Not assigned"),
                                     levels = c(unname(grp_labels), "Not assigned")))
# Palette indexed BY NAME (a positional setNames() once swapped South and Horn).
group_pal <- c(stats::setNames(grp_meta$colour[grp_meta$order], unname(grp_labels)),
               "Not assigned" = "gray80")
p_grp <- ggplot() +
  geom_sf(data = grp_sf, aes(fill = group_label), color = "white", lwd = 0.2) +
  scale_fill_manual(values = group_pal, name = "ENSO group", na.value = "gray80") +
  taxdat::map_theme() +
  theme(panel.border = element_blank(), legend.position = "right") +
  labs(title = "ENSO teleconnection grouping — Africa",
       subtitle = sprintf("%d countries in %d groups (enso_central excluded from models)",
                          nrow(enso_grp), dplyr::n_distinct(enso_grp$group)))
save_fig(p_grp, "01_enso_grouping_map", width = 9, height = 7)

# Sections 1-2: composite maps per track --------------------------------------------
composite_for_track <- function(batch_dir, lags, track, track_label, stems) {
  if (is.null(batch_dir) || !dir.exists(batch_dir) || !length(lags)) {
    message("--", if (track == "occurrence") "occ" else "siz", "_dir not supplied/found or no lags; skipping ", track, ".")
    return(invisible(NULL))
  }
  message("\n=== ", track_label, " composite IRR ===")
  tabs <- Filter(Negate(is.null), purrr::imap(lags, function(lag, grp) load_group_tables(batch_dir, grp, lag, track)))
  adm1_all <- dplyr::bind_rows(lapply(tabs, `[[`, "adm1"))
  ctry_all <- dplyr::bind_rows(lapply(tabs, `[[`, "country"))
  lag_txt  <- paste(sprintf("%s=lag%d", names(lags), lags), collapse = "; ")

  one_map <- function(df, level, stem, unit_col, unit_name, title) {
    if (!nrow(df)) { message("No ", level, " IRR rows for ", track_label, "; skipping."); return(invisible(NULL)) }
    n_tot <- dplyr::n_distinct(df[[unit_col]])
    n_sig <- dplyr::n_distinct(df[[unit_col]][df$significant])
    message(sprintf("  Significant %s (any condition): %d / %d", unit_name, n_sig, n_tot))
    p <- plot_irr_map(df, geo, level = level, enso_levels = names(COMPOSITE_CONDS),
                      title = sprintf("%s — %s", title, track_label),
                      subtitle = sprintf("Best lags: %s\nBlack borders = 95%% CI excludes 1 (%d / %d %s)",
                                         lag_txt, n_sig, n_tot, unit_name))
    if (!is.null(p)) save_fig(p, stem)
  }
  one_map(adm1_all, "adm1",    stems[[1L]], "adm1_gid", "ADM1 units", "Composite ENSO IRR map")
  one_map(ctry_all, "country", stems[[2L]], "country",  "countries",  "Composite ENSO IRR map (country)")
}
composite_for_track(cli$occ_dir, OCC_LAGS, "occurrence", "Outbreak occurrence",
                    c("02_occ_composite_irr", "04_occ_composite_country_irr"))
composite_for_track(cli$siz_dir, SIZ_LAGS, "size", "Outbreak size",
                    c("03_siz_composite_irr", "05_siz_composite_country_irr"))

message("\nAll figures written to: ", fig_dir)
