# 04_fit_occ_brms_cb.R — brms Poisson-GLMM occurrence fit with dlnm
#   crossbasis predictors (the deployed occurrence backend).
#
# Usage:
#   Rscript analysis/04_fit_run.R -C analysis/configs/<batch>/<cfg>.yml
#   Rscript analysis/04_fit_occ_brms_cb.R -C <cfg.yml> \
#     --chains 2 --sampling 100 --warmup 100
#
# Details:
#   Activated via occurrence.backend: brms_cb. One dlnm::crossbasis() per
#   variable in occurrence.dlnm.variables (each may override
#   occurrence.dlnm.defaults). Writes make_run_file() kinds fit, cb_meta and
#   re_stats (if full nesting is present).

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(splines2)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/priors.R") # brms prior builders
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/composite_lp.R") # transitive dependency of population.R
source("analysis/utils/population.R") # .apply_p3() -- occurrence/duration.canonical_admin_level (load_*_data)
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/plots.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-a", "--admin_level"),
    default = NULL,
    type = "character",
    help = "Admin level (ADM1 | ADM2) [config default: ADM2]"
  ),
  make_option(
    c("-b", "--baseline"),
    default = NULL,
    type = "character",
    help = "Baseline period [config default: 1990_2019]"
  ),
  make_option(
    c("-c", "--country"),
    default = NULL,
    type = "character",
    help = "ISO3 code(s), comma-separated (e.g. 'ETH' or 'ETH,KEN'), or 'all'"
  ),
  make_option(
    c("--basehaz_poisson"),
    default = NULL,
    type = "character",
    help = "Baseline hazard: factor|cyclic|seasonal|ms [config default: seasonal]"
  ),
  make_option(
    c("--sampling"),
    default = NULL,
    type = "integer",
    help = "Post-warmup sampling iterations [config default: 1000]"
  ),
  make_option(
    c("--warmup"),
    default = NULL,
    type = "integer",
    help = "Warmup iterations [config default: 1000]"
  ),
  make_option(
    c("--chains"),
    default = NULL,
    type = "integer",
    help = "MCMC chains [config default: 4]"
  ),
  make_option(
    c("--cores"),
    default = NULL,
    type = "integer",
    help = "Parallel cores [config default: 4]"
  ),
  make_option(
    c("--seed"),
    default = NULL,
    type = "integer",
    help = "RNG seed [config default: 1234]"
  ),
  make_option(
    c("--brms_init"),
    default = NULL,
    type = "character",
    help = "Init strategy: random | 0 [config default: random]"
  ),
  make_option(
    c("--overwrite"),
    default = NULL,
    action = "store_true",
    help = "Re-fit even if output file already exists [default: FALSE]"
  )
)

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = c(
    "country",
    "admin_level",
    "basehaz_poisson",
    "sampling",
    "warmup",
    "chains",
    "cores",
    "seed",
    "brms_init",
    "overwrite"
  ),
  model = "occurrence"
)
stopifnot(opt$occurrence$backend == "brms_cb")
print_options(opt)

# brms_cb falls back to brms sub-block for MCMC settings
brms_cfg <- opt$occurrence$brms_cb %||% opt$occurrence$brms %||% list()

# Map "seasonal" -> "cyclic" internally (same spline, preferred label in config).
# Shared with 04_fit_occ_brms.R and with the class-"sds" prior gate, which must
# agree with this line about whether the formula has a smooth in it.
basehaz <- .occ_basehaz(brms_cfg)
stopifnot(basehaz %in% c("factor", "cyclic", "ms"))

# 3. Early exit if fit already exists ------------------------------------------

fit_file <- make_run_file(opt, kind = "fit")
check_run_file_names(opt, "occurrence") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("Fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data + Poisson view --------------------------------------------------
# prepare_occ_brms_data() (utils/occ_data.R): load_surv_data() -> episode
# stamp -> base-rate anchor sidecar (when occurrence.offset_baserate) ->
# .occ_poisson_view() (the single definition of y / exposure / log_offset, which
# the anchor denominator shares) -> optional M-spline baseline columns.

prep <- prepare_occ_brms_data(opt, basehaz)
opt <- prep$opt
pois_dat <- prep$pois_dat
ms_basis <- prep$ms_basis
ep_rec <- prep$ep_rec

# 5. Guard: no data after grouping/country filter ------------------------------

if (nrow(pois_dat) == 0L) {
  message("No rows in pois_dat after country/grouping filter — skipping fit.")
  quit(save = "no", status = 0)
}

# 6. Crossbasis predictors -----------------------------------------------------
# One dlnm::crossbasis() per variable in occurrence.dlnm.variables.
# CB matrices are injected as flat columns, same pattern as M-spline basis.
# Each variable may override the shared occurrence.dlnm.defaults (max_lag,
# argvar, arglag, keep_linear_term, re_adm1_slope).
# Special variable stem "enso": reuses the pre-built
# enso_lag0..enso_lag{max_lag} columns from join_enso_lags()
# (calendar-aligned, no season reset). All other stems: season-aware lagging
# via make_lag_matrix() (boundary → 0).

dlnm_cfg <- opt$occurrence$dlnm %||% list()
cbb <- build_cb_block(pois_dat, dlnm_cfg, " (occurrence model)")
pois_dat <- cbb$dat
cb_cols <- cbb$cb_cols
cb_meta <- cbb$cb_meta

# 7. Build formula from config -------------------------------------------------

n_country <- nlevels(pois_dat$country)
n_region <- length(unique(pois_dat$AFRO_region))
re_sep <- if (isTRUE(brms_cfg$re_uncorrelated)) "||" else "|"

terms <- build_formula_terms(
  opt,
  n_country = n_country,
  n_region = n_region,
  dat = pois_dat
)

# Integrate crossbasis terms into the formula (utils/formula.R): drop linear
# terms with keep_linear_term: false, append CB fixed effects, and add nested
# slope REs where re_adm1_slope: true (grouping mode collapses them).
terms <- integrate_cb_terms(
  terms,
  dlnm_cfg,
  cb_cols,
  cb_meta,
  re_sep,
  drop_region = isTRUE(opt$grouping$enable)
)

baseline_block <- .occ_baseline_block(basehaz, ms_basis)

rhs <- paste(
  c("1", baseline_block, terms$fixed_terms, terms$re_terms, terms$off_terms),
  collapse = " + "
)
form <- as.formula(paste("y ~", rhs))
message("Formula:\n  ", deparse1(form))

# 8. Priors --------------------------------------------------------------------

priors <- build_occ_priors(opt, terms, brms_cfg)

# 9. Init dispatch -------------------------------------------------------------

brms_init_val <- resolve_brms_init(brms_cfg)

# 10. Fit ----------------------------------------------------------------------

mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms_cb")
warmup_val <- mcmc$warmup
sampling_val <- mcmc$sampling
iter_val <- mcmc$iter
refresh_val <- mcmc$refresh
chains_val <- brms_cfg$chains %||% 4L
cores_val <- brms_cfg$cores %||% 4L

threads_arg <- resolve_brms_threads(brms_cfg)

t0 <- Sys.time()
brms_fit <- brms::brm(
  formula = form,
  data = pois_dat,
  family = poisson(),
  prior = priors,
  chains = chains_val,
  iter = iter_val,
  warmup = warmup_val,
  cores = cores_val,
  seed = opt$seed,
  init = brms_init_val,
  refresh = refresh_val,
  backend = "cmdstanr",
  threads = threads_arg,
  control = list(
    adapt_delta = brms_cfg$adapt_delta %||% 0.95,
    max_treedepth = brms_cfg$max_treedepth %||% 12L
  )
)
t1 <- Sys.time()
message(sprintf(
  "Fit time: %.1f min",
  as.numeric(difftime(t1, t0, units = "mins"))
))

# Attach training config so downstream scripts can recover opt without YAML
attr(brms_fit, "opt") <- opt

# 11. Save fit + CB metadata ---------------------------------------------------

attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved fit -> ", fit_file)

if (length(cb_meta) > 0) {
  meta_file <- make_run_file(opt, kind = "cb_meta")
  saveRDS(cb_meta, meta_file)
  message("Saved CB meta -> ", meta_file)
}

# 12. Postprocess RE stats -----------------------------------------------------
# Only attempt if the nested RE groups the adapter expects are present.

n_countries <- dplyr::n_distinct(pois_dat$country)
re_groups <- if (n_countries >= 2L) {
  names(brms::ranef(brms_fit))
} else {
  character(0)
}
drop_region <- isTRUE(opt$grouping$enable)
needed_groups <- if (drop_region) {
  c("country", "country:ADM1")
} else {
  c("AFRO_region", "AFRO_region:country", "AFRO_region:country:ADM1")
}

if (n_countries < 2L) {
  message(
    "Skipping compute_nested_surv_re_stats_brms: single-country / partial nesting."
  )
} else if (all(needed_groups %in% re_groups)) {
  re_stats <- compute_nested_surv_re_stats_brms(
    brms_fit,
    pois_dat,
    drop_region = drop_region
  )
  re_file <- make_run_file(opt, kind = "re_stats")
  saveRDS(re_stats, re_file)
  message("Saved RE stats -> ", re_file)
} else {
  message(
    "Skipping compute_nested_surv_re_stats_brms: needed groups not present (",
    paste(setdiff(needed_groups, re_groups), collapse = ", "),
    ")."
  )
}
