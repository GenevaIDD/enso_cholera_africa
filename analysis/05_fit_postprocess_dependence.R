# 05_fit_postprocess_dependence.R — size <-> duration posterior-dependence check
#
# Usage:
#   Rscript analysis/05_fit_postprocess_dependence.R \
#     --size_config analysis/configs/dur_dependence_smoketest/size_smoke.yml \
#     --dur_config  analysis/configs/dur_smoketest/dur_smoke.yml
#
# Details:
#   Pairs a size and a duration fit trained on the same episodes and asks
#   whether the two outcomes still co-vary after each model's covariates (the
#   simulator draws them independently). Evaluates and reports a
#   pre-registered decision rule; does not implement any coupling.
#   Writes under the size config's outdir dependence_<YYYYMMDD>.rds (full
#   result) and dependence_summary.txt (also printed).
#
# See: analysis/utils/dependence.R

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/dur_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/siz_post.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/dependence.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
script_options <- list(
  make_option(
    c("--size_config"),
    default = NULL,
    type = "character",
    help = "Path to the size-track YAML config (required)"
  ),
  make_option(
    c("--dur_config"),
    default = NULL,
    type = "character",
    help = "Path to the duration-track YAML config (required)"
  ),
  make_option(
    c("--min_n"),
    default = 5L,
    type = "integer",
    help = "Minimum paired episodes required per group [default: %default]"
  )
)
opt_parser <- OptionParser(option_list = script_options)
script_opt <- parse_args(opt_parser)

if (is.null(script_opt$size_config) || !file.exists(script_opt$size_config)) {
  stop("Provide a valid --size_config <path.yml>.")
}
if (is.null(script_opt$dur_config) || !file.exists(script_opt$dur_config)) {
  stop("Provide a valid --dur_config <path.yml>.")
}

# 3. Helpers -------------------------------------------------------------------
# Method (analysis/utils/dependence.R):
#   1. Randomized PIT per episode for each fit (pit_size_nb / pit_duration_episode).
#   2. Kendall's tau between the paired PIT sequences, computed PER POSTERIOR
#      DRAW (not on posterior means) -> a posterior distribution over tau,
#      overall and (optionally) split by ENSO teleconnection group
#      (analysis/groupings/enso_grouping.csv).
#   3. Decomposition: correlate the two fits' ADM1-level intercept estimates
#      (re_stats artefacts) to see whether any dependence is a location
#      effect (some ADM1s are just bigger AND longer) rather than an
#      episode-level residual dependence.

# Backfill an opt list read from YAML with config_defaults (mirrors
# 05_fit_postprocess_loo_lags.R's .fill_defaults()).
.fill_defaults <- function(opt, model) {
  defaults <- get_default_options()
  for (nm in setdiff(names(defaults), c("occurrence", "size", "duration"))) {
    if (is.null(opt[[nm]])) opt[[nm]] <- defaults[[nm]]
  }
  model_defaults <- defaults[[model]] %||% list()
  for (nm in names(model_defaults)) {
    if (is.null(opt[[model]][[nm]])) opt[[model]][[nm]] <- model_defaults[[nm]]
  }
  opt
}

# 4. Load the fits -------------------------------------------------------------
message("=== Size fit ===")
size_opt <- yaml::read_yaml(script_opt$size_config)
size_opt <- .fill_defaults(size_opt, "size")
size_fit_path <- make_run_file(size_opt, kind = "size_fit")
if (!file.exists(size_fit_path)) {
  stop("Size fit not found: ", size_fit_path)
}
message("Loading: ", size_fit_path)
size_fit <- readRDS(size_fit_path)
size_dat <- read_size_data(size_opt)

message("\n=== Duration fit ===")
dur_opt <- yaml::read_yaml(script_opt$dur_config)
dur_opt <- .fill_defaults(dur_opt, "duration")
dur_fit_path <- make_run_file(dur_opt, kind = "dur_fit")
if (!file.exists(dur_fit_path)) {
  stop("Duration fit not found: ", dur_fit_path)
}
message("Loading: ", dur_fit_path)
dur_fit <- readRDS(dur_fit_path)
dur_dat <- load_duration_data(dur_opt)

message(sprintf(
  "\nSize episodes: %d   Duration episodes: %d   Common: %d",
  dplyr::n_distinct(size_dat$episode_id),
  dplyr::n_distinct(dur_dat$episode_id),
  length(intersect(size_dat$episode_id, dur_dat$episode_id))
))

# 5. PIT + Kendall's tau (method steps 1-2) ------------------------------------

message("\nComputing size-side PIT (randomized NB quantile)...")
pit_a <- pit_size_nb(size_fit, size_dat)

message("Computing duration-side PIT (chained censored-hazard quantile)...")
pit_b <- pit_duration_episode(dur_fit, dur_dat)

# ENSO teleconnection group, keyed by country (analysis/groupings/enso_grouping.csv
# -- not under data/, a repo-tracked registry; see CLAUDE.md "ENSO teleconnection
# groupings"). Mapped from episode_id -> group via size_dat's own country column
# (episode_id -> country is 1:1 within a track; size_dat carries `country`).
enso_grp_csv <- "analysis/groupings/enso_grouping.csv"
group_vec <- NULL
if (file.exists(enso_grp_csv) && "country" %in% names(size_dat)) {
  grp_tbl <- utils::read.csv(enso_grp_csv, stringsAsFactors = FALSE)
  ep_country <- stats::setNames(size_dat$country, size_dat$episode_id)
  common_ep <- intersect(pit_a$episode_id, pit_b$episode_id)
  ctry_lut <- stats::setNames(grp_tbl$group, grp_tbl$iso3)
  group_vec <- stats::setNames(
    unname(ctry_lut[ep_country[common_ep]]),
    common_ep
  )
} else {
  message(
    "No enso_grouping.csv / no `country` column found -- skipping by-group tau."
  )
}

message("\nComputing paired Kendall's tau (propagated through draws)...")
tau_res <- paired_kendall_tau(
  pit_a,
  pit_b,
  group = group_vec,
  min_n = script_opt$min_n
)

# 6. RE-intercept decomposition (method step 3) --------------------------------
#
# Track-specific extraction: size's adm1 table keys the intercept row on
# `var == "Intercept"`; duration's re_stats_adm1 table (reused occurrence
# adapter) keys it on `term == "Intercept"`. Both share `country` +
# `adm1_unit` as the join key (see analysis/utils/dependence.R's
# re_intercept_correlation() doc comment for why this extraction is NOT
# pushed into that shared helper).

message("\nComputing ADM1-level intercept correlation (decomposition)...")
re_corr <- tryCatch(
  {
    size_re_path <- make_run_file(size_opt, kind = "size_re_stats")
    dur_re_path <- make_run_file(dur_opt, kind = "dur_re_stats")
    if (!file.exists(size_re_path) || !file.exists(dur_re_path)) {
      message(
        "  Skipping: re_stats artefact(s) not found (single-country / smoke fit?)."
      )
      message(sprintf(
        "    size: %s (%s)",
        size_re_path,
        file.exists(size_re_path)
      ))
      message(sprintf(
        "    dur:  %s (%s)",
        dur_re_path,
        file.exists(dur_re_path)
      ))
      NULL
    } else {
      size_re <- readRDS(size_re_path)
      dur_re <- readRDS(dur_re_path)

      size_adm1 <- size_re$adm1 %>%
        dplyr::filter(var == "Intercept") %>%
        dplyr::transmute(country, adm1_unit, value = mean)

      dur_adm1 <- dur_re$re_stats_adm1 %>%
        dplyr::filter(term == "Intercept") %>%
        dplyr::transmute(country, adm1_unit, value = mean)

      re_intercept_correlation(size_adm1, dur_adm1)
    }
  },
  error = function(e) {
    message("  RE-intercept correlation failed: ", conditionMessage(e))
    NULL
  }
)

# 7. Decision rule -------------------------------------------------------------
# Pre-registered decision rule:
#   |tau| < 0.1 AND the 95% CrI on tau contains 0  => independent draws are
#     adequate for the simulator (no change to the joint sampling).
#   otherwise => couple the draws, first by adding log(size) as a covariate to
#     the duration model. This script only evaluates and reports the rule; the
#     size-model-with-log(duration_months) sensitivity configs
#     (analysis/configs/dur_dependence_smoketest/) give the mechanical-confound
#     check (fit unchanged, config-only covariate add).

ov <- tau_res$overall
independent <- (abs(ov$mean) < 0.1) && (ov$q025 <= 0) && (ov$q975 >= 0)

summary_lines <- c(
  "=== Phase 5: size <-> duration dependence ===",
  sprintf("Size fit:     %s", size_fit_path),
  sprintf("Duration fit: %s", dur_fit_path),
  sprintf("Paired episodes: %d", tau_res$n_common),
  "",
  "Overall Kendall's tau (posterior distribution over draws):",
  sprintf(
    "  mean = %.3f  sd = %.3f  95%% CrI = [%.3f, %.3f]  (n_episodes = %d)",
    ov$mean,
    ov$sd,
    ov$q025,
    ov$q975,
    ov$n_episodes
  )
)

if (!is.null(tau_res$by_group)) {
  summary_lines <- c(summary_lines, "", "By ENSO group:")
  for (g in names(tau_res$by_group)) {
    gg <- tau_res$by_group[[g]]
    summary_lines <- c(
      summary_lines,
      sprintf(
        "  %-14s mean = %.3f  95%% CrI = [%.3f, %.3f]  (n = %d)",
        g,
        gg$mean,
        gg$q025,
        gg$q975,
        gg$n_episodes
      )
    )
  }
}

if (!is.null(re_corr)) {
  summary_lines <- c(
    summary_lines,
    "",
    sprintf(
      "ADM1 intercept correlation (size vs duration, response scale): r = %.3f (n = %d ADM1 units)",
      re_corr$r,
      re_corr$n
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  strrep("-", 60),
  sprintf(
    "Decision rule: |tau| < 0.1 AND 95%% CrI contains 0  =>  %s",
    if (independent) {
      "INDEPENDENT (no coupling needed)"
    } else {
      "COUPLE draws (plan's E1: add log(size) to the duration model)"
    }
  ),
  strrep("-", 60)
)

message("\n", paste(summary_lines, collapse = "\n"))

# 8. Write outputs -------------------------------------------------------------

outdir <- size_opt$outdir %||% "generated_data/"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

out_rds <- file.path(
  outdir,
  sprintf("dependence_%s.rds", format(Sys.Date(), "%Y%m%d"))
)
saveRDS(
  list(
    tau = tau_res,
    re_corr = re_corr,
    independent = independent,
    size_fit_path = size_fit_path,
    dur_fit_path = dur_fit_path
  ),
  out_rds
)
message("\nSaved -> ", out_rds)

out_txt <- file.path(outdir, "dependence_summary.txt")
writeLines(summary_lines, out_txt)
message("Saved -> ", out_txt)
