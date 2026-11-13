# 04_fit_siz_brms.R — config-driven brms outbreak-size fit
#
# Status: SCAFFOLD. Implements the config /
# CLI parsing and the high-level flow. Helper functions
# (build_size_formula(), build_size_priors(), readSizeData(),
# computeNestedSizeReStatsBrms()) are stubbed in analysis/utils/config.R and
# analysis/utils.R.
#
# Legacy: analysis/07_run_outbreak_size_model.R (kept intact).
#
# Examples:
#   Rscript analysis/04_fit_siz_brms.R \
#     -C analysis/configs/smoke_test/smoke_size_brms.yml
#   Rscript analysis/04_fit_siz_brms.R \
#     -C analysis/configs/size_run1_benchmark/size_run1_brms_default.yml --country ETH

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/priors.R")                # brms prior builders
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/siz_post.R")
source("analysis/utils/anchor.R")   # transitive dependency of siz_data.R (utils/README.md)

option_list <- list(
  make_option(c("-C", "--config"),
              default = NULL, type = "character",
              help = "Path to YAML config file"),
  make_option(c("-a", "--admin_level"),
              default = NULL, type = "character",
              help = "Admin level (ADM1 | ADM2) [size default: ADM1]"),
  make_option(c("-b", "--baseline"),
              default = NULL, type = "character",
              help = "Baseline period [default: 1990_2019]"),
  make_option(c("-c", "--country"),
              default = NULL, type = "character",
              help = "ISO3 code(s), comma-separated, or 'all'"),
  make_option(c("--family"),
              default = NULL, type = "character",
              help = "brms family: negbinomial | poisson [size default: negbinomial]"),
  make_option(c("--seasonal_term"),
              default = NULL, type = "character",
              help = "none | factor_month | cyclic | harmonic2 [size default: none]"),
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
  enforce_options = c("country", "baseline", "admin_level"),
  model           = "size"
)
print(opt)

# Validate covariate ↔ RE consistency (size-specific rule)
validate_size_options(opt)

# Load data ---------------------------------------------------------------

size_dat <- readSizeData(opt)
# Episode record of the prep artefact this fit trains on (assertEpisodeKey()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(size_dat, "episode_stamp")
persistMeanRateTable(opt, size_dat)

# Build model spec --------------------------------------------------------

parts     <- build_size_formula_parts(opt, size_dat)
bf_obj    <- build_size_formula(opt, size_dat)
prior_obj <- build_size_priors(opt, parts = parts)

fit_path  <- makeRunFile(opt, kind = "size_fit")
checkRunFileNames(opt, "size")   # every artefact name fits NAME_MAX before the fit (paths.R)
message("size_fit output -> ", fit_path)

# Fit ---------------------------------------------------------------------

brms_cfg <- opt$size$brms %||% list()

if (file.exists(fit_path) && !isTRUE(opt$overwrite)) {
  message("size_fit already exists; skipping (use --overwrite to refit).")
  stan_fit <- readRDS(fit_path)
} else {
  mcmc         <- resolve_brms_mcmc(brms_cfg, label = "brms")
  warmup_val   <- mcmc$warmup
  sampling_val <- mcmc$sampling
  iter_val     <- mcmc$iter
  refresh_val  <- mcmc$refresh
  chains_val  <- brms_cfg$chains  %||% 4L
  cores_val   <- brms_cfg$cores   %||% 4L

  threads_per_chain <- brms_cfg$threads_per_chain %||% NULL
  threads_arg <- if (!is.null(threads_per_chain) && as.integer(threads_per_chain) > 1L) {
    message(sprintf("Multithreading enabled: %d threads per chain", as.integer(threads_per_chain)))
    brms::threading(as.integer(threads_per_chain))
  } else {
    NULL
  }

  t0 <- Sys.time()
  stan_fit <- brms::brm(
    formula  = bf_obj,
    data     = size_dat,
    family   = match.fun(opt$size$family %||% "negbinomial")(opt$size$link %||% "log"),
    iter     = iter_val,
    warmup   = warmup_val,
    chains   = chains_val,
    cores    = cores_val,
    refresh  = refresh_val,
    control  = list(adapt_delta   = brms_cfg$adapt_delta   %||% 0.95,
                    max_treedepth = brms_cfg$max_treedepth %||% 12L),
    init     = brms_cfg$init %||% 0.1,
    prior    = prior_obj,
    seed     = opt$seed,
    backend  = "cmdstanr",
    threads  = threads_arg,
    file     = tools::file_path_sans_ext(fit_path),
    save_model = file.path(dirname(fit_path),
                           paste0(tools::file_path_sans_ext(basename(fit_path)), ".stan"))
  )
  t1 <- Sys.time()
  message(sprintf("Fit time: %.1f min", as.numeric(difftime(t1, t0, units = "mins"))))

  # Layer 4 (artefact chain, see popStampSidecarPath() in utils.R): record the
  # population THIS fit was actually trained on, in a small companion sidecar
  # next to fit_path. Deliberately inside this branch only (a fresh fit) --
  # the skip-reuse branch above loads an OLD fit trained on a population that
  # may not be size_dat's, so stamping there would misattribute the CURRENT
  # run's population to a fit that never saw it.
  persistFitPopStamp(fit_path, getPopProvenance(size_dat))
}

# Attach training config as attribute so 09_simulator.R can recover offset flags,
# admin_level, and backend without re-parsing the YAML.
attr(stan_fit, "opt") <- opt
# ... and the estimand stamp, which carries `log_pop_centres`. Under
# `pop_elasticity: free` the simulator MUST re-centre its own frame units on
# these training centres (`.recover_log_pop_centres()`, utils/simulator_core.R);
# brms keeps only formula variables in `fit$data`, so without this the centres
# have to be reverse-engineered from log(pop) - log_pop_c.
attr(stan_fit, "size_estimand") <- getSizeEstimand(size_dat)
attr(stan_fit, "episode_stamp") <- ep_rec
saveRDS(stan_fit, fit_path)

message("done: ", fit_path)
