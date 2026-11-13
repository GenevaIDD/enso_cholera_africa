# 04_fit_occ_brms.R — config-driven brms Poisson-GLMM fits
#
# Equivalence argument: a survival model on monthly counting-process rows is
# algebraically equivalent to a Poisson GLM with offset(log(interval_length))
# under the Holford-form baseline (Holford 1980; Laird & Olivier 1981).
# Three baseline-hazard variants:
#   basehaz_poisson = "factor"   Holford-exact: + factor(months)
#   basehaz_poisson = "cyclic"   Cyclic spline:  + s(t_start_seas, bs="cc", k=8)
#   basehaz_poisson = "ms"       M-spline basis as fixed cols
#   basehaz_poisson = "seasonal" Cyclic spline (same as cyclic; preferred label)
#
# Reads all parameters from a YAML config (--config / -C) or from CLI flags.
# Priority: explicit CLI flag > YAML config value > default (config_defaults.yml).
#
# Init dispatch (brms_init):
#   "random"  standard brms random initialisation (default)
#   "0"       initialise all parameters at 0
#
# Examples:
#   Rscript analysis/04_fit_occ_brms.R \
#     -C analysis/configs/<batch>/<config>.yml
#   Rscript analysis/04_fit_occ_brms.R \
#     -C analysis/configs/<batch>/<config>.yml --country ETH
#   Rscript analysis/04_fit_occ_brms.R \
#     --country all --basehaz_poisson seasonal --iter 4000
#
# Legacy: analysis/04b_run_poisson_models_legacy.R is the original script.

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(splines2)
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
source("analysis/utils/composite_lp.R")   # transitive dependency of population.R
source("analysis/utils/population.R")     # .applyP3() -- occurrence/duration.canonical_admin_level (load_*_data)
source("analysis/utils/siz_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/plots.R")
source("analysis/utils/anchor.R")

option_list <- list(
  make_option(c("-a", "--admin_level"),
              default = NULL, type = "character",
              help = "Admin level (ADM1 | ADM2) [config default: ADM2]"),
  make_option(c("-b", "--baseline"),
              default = NULL, type = "character",
              help = "Baseline period [config default: 1990_2019]"),
  make_option(c("-c", "--country"),
              default = NULL, type = "character",
              help = "ISO3 code(s), comma-separated (e.g. 'ETH' or 'ETH,KEN'), or 'all' [config default: ETH]"),
  make_option(c("--basehaz_poisson"),
              default = NULL, type = "character",
              help = "Baseline hazard: factor|cyclic|seasonal|ms [config default: seasonal]"),
  make_option(c("--sampling"),
              default = NULL, type = "integer",
              help = "Post-warmup sampling iterations [config default: 1000]"),
  make_option(c("--warmup"),
              default = NULL, type = "integer",
              help = "Warmup iterations [config default: 1000]"),
  make_option(c("--chains"),
              default = NULL, type = "integer",
              help = "MCMC chains [config default: 4]"),
  make_option(c("--cores"),
              default = NULL, type = "integer",
              help = "Parallel cores [config default: 4]"),
  make_option(c("--seed"),
              default = NULL, type = "integer",
              help = "RNG seed [config default: 1234]"),
  make_option(c("--brms_init"),
              default = NULL, type = "character",
              help = "Init strategy: random | 0 [config default: random]"),
  make_option(c("--overwrite"),
              default = NULL, action = "store_true",
              help = "Re-fit even if output file already exists [default: FALSE]")
)

opt <- make_options_from_config(option_list     = option_list,
                                 enforce_options = c("country", "admin_level", "basehaz_poisson",
                                                     "sampling", "warmup", "chains",
                                                     "cores", "seed", "brms_init",
                                                     "overwrite"),
                                 model           = "occurrence")
stopifnot(opt$occurrence$backend == "brms")
print_options(opt)

brms_cfg <- opt$occurrence$brms %||% list()

# Map "seasonal" -> "cyclic" internally (same spline, preferred label in config).
# Shared with 04_fit_occ_brms_cb.R and with the class-"sds" prior gate, which
# must agree with this line about whether the formula has a smooth in it.
basehaz <- .occ_basehaz(brms_cfg)
stopifnot(basehaz %in% c("factor", "cyclic", "ms"))

# Early exit if fit already exists ----------------------------------------

fit_file <- makeRunFile(opt, kind = "fit")
checkRunFileNames(opt, "occurrence")   # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("Fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# Load data + Poisson view -------------------------------------------------
# prepare_occ_brms_data() (utils/occ_data.R): load_surv_data() -> episode
# stamp -> base-rate anchor sidecar (when occurrence.offset_baserate) ->
# .occPoissonView() (the single definition of y / exposure / log_offset, which
# the anchor denominator shares) -> optional M-spline baseline columns.

prep     <- prepare_occ_brms_data(opt, basehaz)
opt      <- prep$opt
pois_dat <- prep$pois_dat
ms_basis <- prep$ms_basis
ep_rec   <- prep$ep_rec

# Build formula from config -----------------------------------------------

n_country <- nlevels(pois_dat$country)
n_region  <- length(unique(pois_dat$AFRO_region))

terms <- build_formula_terms(opt, n_country = n_country, n_region = n_region, dat = pois_dat)

baseline_block <- .occ_baseline_block(basehaz, ms_basis)

rhs <- paste(
  c("1", baseline_block, terms$fixed_terms, terms$re_terms, terms$off_terms),
  collapse = " + "
)
form <- as.formula(paste("y ~", rhs))
message("Formula:\n  ", deparse1(form))

# Priors ------------------------------------------------------------------

priors <- build_occ_priors(opt, terms, brms_cfg)

# Init dispatch -----------------------------------------------------------

brms_init_str <- brms_cfg$init %||% "random"
brms_init_val <- if (is.numeric(brms_init_str)) {
  brms_init_str
} else switch(brms_init_str,
  random = "random",
  `0`    = 0,
  {
    warning("Unknown brms_init value '", brms_init_str, "' — using random.")
    "random"
  }
)

# Fit ---------------------------------------------------------------------

mcmc         <- resolve_brms_mcmc(brms_cfg, label = "brms")
warmup_val   <- mcmc$warmup
sampling_val <- mcmc$sampling
iter_val     <- mcmc$iter
refresh_val  <- mcmc$refresh
chains_val  <- brms_cfg$chains  %||% 4L
cores_val   <- brms_cfg$cores   %||% 4L

threads_arg <- resolve_brms_threads(brms_cfg)

t0 <- Sys.time()
brms_fit <- brms::brm(
  formula = form,
  data    = pois_dat,
  family  = poisson(),
  prior   = priors,
  chains  = chains_val,
  iter    = iter_val,
  warmup  = warmup_val,
  cores   = cores_val,
  seed    = opt$seed,
  init    = brms_init_val,
  refresh = refresh_val,
  backend = "cmdstanr",
  threads = threads_arg,
  control = list(adapt_delta   = brms_cfg$adapt_delta   %||% 0.95,
                 max_treedepth = brms_cfg$max_treedepth %||% 12L)
)
t1 <- Sys.time()
message(sprintf("Fit time: %.1f min", as.numeric(difftime(t1, t0, units = "mins"))))

# Attach training config so downstream scripts can recover opt without YAML
attr(brms_fit, "opt") <- opt

# Save fit ----------------------------------------------------------------
attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved fit -> ", fit_file)

# Postprocess -------------------------------------------------------------
# Only attempt if the nested RE groups the adapter expects are present.

n_countries <- dplyr::n_distinct(pois_dat$country)

re_groups     <- if (n_countries >= 2L) names(brms::ranef(brms_fit)) else character(0)
drop_region   <- isTRUE(opt$grouping$enable)
needed_groups <- if (drop_region) {
  c("country", "country:ADM1")
} else {
  c("AFRO_region", "AFRO_region:country", "AFRO_region:country:ADM1")
}
# "AFRO_region:country:ADM1:gid" / "country:ADM1:gid" is optional (re_adm2_intercept).

if (n_countries < 2L) {
  message("Skipping computeNestedSurvReStatsBrms: single-country / partial nesting.")
} else if (all(needed_groups %in% re_groups)) {
  re_stats <- computeNestedSurvReStatsBrms(brms_fit, pois_dat, drop_region = drop_region)
  re_file  <- makeRunFile(opt, kind = "re_stats")
  saveRDS(re_stats, re_file)
  message("Saved RE stats -> ", re_file)
} else {
  message("Skipping computeNestedSurvReStatsBrms: needed groups not present (",
          paste(setdiff(needed_groups, re_groups), collapse = ", "),
          "). Likely a single-country / Phase-1 fit.")
}
