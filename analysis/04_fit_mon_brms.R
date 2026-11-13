# 04_fit_mon_brms.R — brms monthly-size (within-episode, per-calendar-month)
# negative-binomial model, parametric only (no dlnm crossbasis)
#
# Fourth model track, alongside occurrence/size/duration. Sibling to
# 04_fit_mon_brms_cb.R (kept intact per side-by-side rule). Activated via
# monthly.backend: brms in the YAML config. This is the parametric-only
# baseline used by the structure-comparison LOO arm (mon_structure_compare
# batch) — no distributed-/single-lag crossbasis terms, just the linear
# fixed-effect covariates build_monthly_formula_parts() resolves from
# monthly.covariates.
#
#   cases_month_{i,m} ~ NegBin(mu_{i,m}, shape)
#   mu_{i,m} = exp( admin_level + enso_lag<N> + [seasonal]
#                   + (RE | AFRO_region/country/ADM1/gid) )
#              * pop_i * exposure_days_{i,m}
#
# Output files (via makeRunFile / makeMonthlyRunFile):
#   kind = "mon_fit" — brms fit object
#
# No persistFitPopStamp() — see 04_fit_mon_brms_cb.R's header for why.
# RE stats / coef summary are centralised in 05_fit_postprocess.R.
#
# MCMC settings: opt$monthly$brms.
#
# Usage:
#   Rscript analysis/04_fit_run.R -C analysis/configs/<batch>/<cfg>.yml

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(lubridate)  # %m-% -- monthly.prev_cases: yes (mon_data.R:67)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/priors.R")                # brms prior builders
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/dur_data.R")
source("analysis/utils/mon_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R")

option_list <- list(
  make_option(c("-C", "--config"),
              default = NULL, type = "character",
              help = "Path to YAML config file"),
  make_option(c("-a", "--admin_level"),
              default = NULL, type = "character",
              help = "Admin level (ADM1 | ADM2 | ADM1-ADM2) [monthly default: ADM1-ADM2]"),
  make_option(c("-b", "--baseline"),
              default = NULL, type = "character",
              help = "Baseline period [default: 1990_2019]"),
  make_option(c("-c", "--country"),
              default = NULL, type = "character",
              help = "ISO3 code(s), comma-separated, or 'all'"),
  make_option(c("--sampling"),
              default = NULL, type = "integer",
              help = "Post-warmup sampling iterations"),
  make_option(c("--warmup"),
              default = NULL, type = "integer",
              help = "Warmup iterations"),
  make_option(c("--chains"),
              default = NULL, type = "integer",
              help = "MCMC chains"),
  make_option(c("--cores"),
              default = NULL, type = "integer",
              help = "Parallel cores"),
  make_option(c("--seed"),
              default = NULL, type = "integer",
              help = "RNG seed"),
  make_option(c("--overwrite"),
              default = NULL, action = "store_true",
              help = "Re-fit even if output file already exists")
)

opt <- make_options_from_config(
  option_list     = option_list,
  enforce_options = c("country", "baseline", "admin_level",
                      "sampling", "warmup", "chains", "cores",
                      "seed", "overwrite"),
  model           = "monthly"
)
stopifnot(!is.null(opt$monthly), opt$monthly$backend == "brms")
print_options(opt)

brms_cfg <- opt$monthly$brms %||% list()

# Early exit if fit already exists ----------------------------------------

fit_file <- makeRunFile(opt, kind = "mon_fit")
checkRunFileNames(opt, "monthly")   # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("mon_fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# Load data -----------------------------------------------------------------

mon_dat <- load_monthly_data(opt)
# Episode record of the prep artefact this fit trains on (assertEpisodeKey()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(mon_dat, "episode_stamp")

# Anchor sidecar (utils/anchor.R) -- the monthly analogue of
# `persistMeanRateTable(opt, size_dat)` (04_fit_siz_brms.R:89). Gated on
# monthly.offset_baserate, so it is inert for every config that predates it.
if (anchorEnabled(opt, "mon")) persistBaseRateTable(opt, "mon", mon_dat)

message(sprintf("Monthly panel: %d rows, %d episodes, %d total cases",
                nrow(mon_dat),
                dplyr::n_distinct(mon_dat$episode_id),
                round(sum(mon_dat$cases_month, na.rm = TRUE))))
if (nrow(mon_dat) == 0L)
  stop("Monthly panel is empty after subsetting — check `country` / `grouping`.")

# Build formula -------------------------------------------------------------

drop_region <- isTRUE(opt$grouping$enable)
parts       <- build_monthly_formula_parts(opt, mon_dat)
bf_obj      <- build_monthly_formula(opt, mon_dat)

message("Formula:\n  ", deparse1(bf_obj$formula))

# Priors ------------------------------------------------------------------

prior_obj <- build_monthly_priors(opt, parts = parts)

# Fit ---------------------------------------------------------------------

mcmc         <- resolve_brms_mcmc(brms_cfg, label = "brms")
warmup_val   <- mcmc$warmup
sampling_val <- mcmc$sampling
iter_val     <- mcmc$iter
refresh_val  <- mcmc$refresh
chains_val   <- brms_cfg$chains %||% 4L
cores_val    <- brms_cfg$cores  %||% 4L

threads_per_chain <- brms_cfg$threads_per_chain %||% NULL
threads_arg <- if (!is.null(threads_per_chain) &&
                   as.integer(threads_per_chain) > 1L) {
  message(sprintf("Multithreading enabled: %d threads per chain",
                  as.integer(threads_per_chain)))
  brms::threading(as.integer(threads_per_chain))
} else {
  NULL
}

t0 <- Sys.time()
brms_fit <- brms::brm(
  formula  = bf_obj,
  data     = mon_dat,
  family   = match.fun(opt$monthly$family %||% "negbinomial")(opt$monthly$link %||% "log"),
  prior    = prior_obj,
  chains   = chains_val,
  iter     = iter_val,
  warmup   = warmup_val,
  cores    = cores_val,
  seed     = opt$seed,
  init     = brms_cfg$init %||% 0.1,
  refresh  = refresh_val,
  backend  = "cmdstanr",
  threads  = threads_arg,
  control  = list(adapt_delta   = brms_cfg$adapt_delta   %||% 0.95,
                  max_treedepth = brms_cfg$max_treedepth %||% 12L)
)
t1 <- Sys.time()
message(sprintf("Fit time: %.1f min", as.numeric(difftime(t1, t0, units = "mins"))))

# Attach training config so downstream scripts can recover opt without YAML
attr(brms_fit, "opt") <- opt

# Attach the estimand stamp beside it. `opt` alone is NOT sufficient: the stamp
# also carries `log_pop_centres`, the P4 pivots the simulator must re-centre
# frame units on. Recomputing them from a frame's own mean would apply the
# elasticity around the wrong population and silently rescale every prediction.
attr(brms_fit, "mon_estimand") <- getMonthlyEstimand(mon_dat)

# Save fit ------------------------------------------------------------------

attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved mon_fit -> ", fit_file)
