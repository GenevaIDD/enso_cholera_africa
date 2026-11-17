# 04_fit_siz_brms_cb.R — brms negative-binomial size model with dlnm crossbasis
#
# Usage:
#   Rscript analysis/04_fit_run.R \
#     -C analysis/configs/smoke_cb_precip_siz/smoke_cb_precip_siz.yml
#
# Details:
#   Selected by size.backend: brms_cb. One dlnm::crossbasis() per variable in
#   size.dlnm.variables (per-variable spec overrides size.dlnm.defaults).
#   MCMC settings: opt$size$brms_cb, falling back to opt$size$brms.
#   Writes (make_run_file kinds) size_fit (with attr "opt"), size_cb_meta (CB
#   metadata per variable), size_re_stats (if full nesting is present).

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
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
source("analysis/utils/siz_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/siz_post.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-C", "--config"),
    default = NULL,
    type = "character",
    help = "Path to YAML config file"
  ),
  make_option(
    c("-a", "--admin_level"),
    default = NULL,
    type = "character",
    help = "Admin level (ADM1 | ADM2) [size default: ADM1]"
  ),
  make_option(
    c("-b", "--baseline"),
    default = NULL,
    type = "character",
    help = "Baseline period [default: 1990_2019]"
  ),
  make_option(
    c("-c", "--country"),
    default = NULL,
    type = "character",
    help = "ISO3 code(s), comma-separated, or 'all'"
  ),
  make_option(
    c("--family"),
    default = NULL,
    type = "character",
    help = "brms family: negbinomial | poisson [size default: negbinomial]"
  ),
  make_option(
    c("--sampling"),
    default = NULL,
    type = "integer",
    help = "Post-warmup sampling iterations"
  ),
  make_option(
    c("--warmup"),
    default = NULL,
    type = "integer",
    help = "Warmup iterations"
  ),
  make_option(
    c("--chains"),
    default = NULL,
    type = "integer",
    help = "MCMC chains"
  ),
  make_option(
    c("--cores"),
    default = NULL,
    type = "integer",
    help = "Parallel cores"
  ),
  make_option(c("--seed"), default = NULL, type = "integer", help = "RNG seed"),
  make_option(
    c("--overwrite"),
    default = NULL,
    action = "store_true",
    help = "Re-fit even if output file already exists"
  )
)

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = c(
    "country",
    "baseline",
    "admin_level",
    "sampling",
    "warmup",
    "chains",
    "cores",
    "seed",
    "overwrite"
  ),
  model = "size"
)
stopifnot(opt$size$backend == "brms_cb")
print_options(opt)

# brms_cb falls back to brms sub-block for MCMC settings (same pattern as
# 04b_run_occ_brms_cb.R)
brms_cfg <- opt$size$brms_cb %||% opt$size$brms %||% list()

# 3. Early exit if fit already exists ------------------------------------------

fit_file <- make_run_file(opt, kind = "size_fit")
check_run_file_names(opt, "size") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("size_fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data -----------------------------------------------------------------

size_dat <- read_size_data(opt)
# Episode record of the prep artefact this fit trains on (assert_episode_key()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(size_dat, "episode_stamp")
persist_mean_rate_table(opt, size_dat)
# Captured HERE, not at save time: build_cb_block() below rebuilds `size_dat`
# and any non-standard attribute (this stamp, which carries the P4
# `log_pop_centres` the simulator re-centres on) would not survive the rebuild.
size_estimand_stamp <- get_size_estimand(size_dat)

# 5. Crossbasis predictors -----------------------------------------------------
# Lag matrices are read from pre-computed <var>_lag1 … <var>_lag{max_lag}
# columns written by 03b_prepare_siz.R (make_climate_lags mode="event").
# If a required column is absent, make_lag_matrix() raises a clear error with
# the regeneration command.
# One dlnm::crossbasis() per variable in size.dlnm.variables.
# CB matrices are flattened to <var>_cb1 … <var>_cbN columns.

dlnm_cfg <- opt$size$dlnm %||% list()
cbb <- build_cb_block(size_dat, dlnm_cfg, " (size model)")
size_dat <- cbb$dat
cb_cols <- cbb$cb_cols
cb_meta <- cbb$cb_meta

# 6. Build formula from parts --------------------------------------------------

drop_region <- isTRUE(opt$grouping$enable)
parts <- build_size_formula_parts(opt, size_dat)

# Integrate crossbasis terms into the formula
if (length(cb_cols) > 0) {
  dlnm_defaults <- dlnm_cfg$defaults %||% list()

  # Drop plain FE + slope RE for CB variables with keep_linear_term: false
  for (v in names(dlnm_cfg$variables)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    klt <- var_spec$keep_linear_term %||%
      dlnm_defaults$keep_linear_term %||%
      FALSE
    if (!isTRUE(klt)) {
      parts$fixed_terms <- parts$fixed_terms[parts$fixed_terms != v]
      parts$re_terms <- parts$re_terms[
        !grepl(sprintf("^\\(%s[[:space:]]", v), parts$re_terms)
      ]
    }
  }

  # Append CB fixed-effect terms
  parts$fixed_terms <- c(parts$fixed_terms, cb_cols)

  # Append CB slope RE terms — hierarchical (region/country/ADM1) slopes, opt-in
  # per variable via re_adm1_slope: true. Shared plumbing with the occ CB script.
  # drop_region = TRUE (grouping mode) collapses to country/ADM1 so groups with a
  # single AFRO_region (e.g. enso_south) don't hit the 1-level contrast error.
  # filter_re_terms_for_data() then drops any remaining terms whose grouping
  # column has <2 distinct values (e.g. single-country enso_south size data).
  cb_groups <- re_group_by_level("ADM2", drop_region = drop_region)
  re_sep <- if (isTRUE(brms_cfg$re_uncorrelated)) "||" else "|"
  for (v in names(dlnm_cfg$variables)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    re_slope <- var_spec$re_adm1_slope %||%
      (dlnm_defaults$re_adm1_slope %||% FALSE)
    if (!isTRUE(re_slope)) {
      next
    }
    v_cb_names <- cb_meta[[v]]$cb_col_names
    if (is.null(v_cb_names)) {
      next
    }
    new_cb_re <- filter_re_terms_for_data(
      nested_slope_re_terms(v_cb_names, cb_groups, re_sep),
      size_dat
    )
    parts$re_terms <- c(parts$re_terms, new_cb_re)
  }
}

# Assemble brms::bf()
phi_re <- parts$phi_re
rhs <- paste(
  c(parts$fixed_terms, parts$seas_terms, parts$re_terms, parts$off_terms),
  collapse = " + "
)
form_str <- paste("suspected_cases ~", rhs)
bf_main <- brms::bf(stats::as.formula(form_str))

# Mirrors build_size_formula() in utils/config.R: country-grouping mode (drop_region)
# re-roots the phi hierarchy at country instead of AFRO_region; sparse levels
# fall back to constant phi via .phi_re_term()'s filter_re_terms_for_data() guard.
phi_term <- .phi_re_term(phi_re, drop_region, size_dat)
bf_obj <- if (length(phi_term) == 0L) {
  bf_main
} else {
  bf_main + brms::lf(stats::as.formula(paste("shape ~ 1 +", phi_term)))
}

message("Formula:\n  ", deparse1(bf_obj$formula))

# 7. Priors --------------------------------------------------------------------

prior_obj <- build_size_priors(opt, brms_cfg, parts)

# 8. Fit -----------------------------------------------------------------------

mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms_cb")
warmup_val <- mcmc$warmup
sampling_val <- mcmc$sampling
iter_val <- mcmc$iter
refresh_val <- mcmc$refresh
chains_val <- brms_cfg$chains %||% 4L
cores_val <- brms_cfg$cores %||% 4L

threads_per_chain <- brms_cfg$threads_per_chain %||% NULL
threads_arg <- if (
  !is.null(threads_per_chain) &&
    as.integer(threads_per_chain) > 1L
) {
  message(sprintf(
    "Multithreading enabled: %d threads per chain",
    as.integer(threads_per_chain)
  ))
  brms::threading(as.integer(threads_per_chain))
} else {
  NULL
}

t0 <- Sys.time()
brms_fit <- brms::brm(
  formula = bf_obj,
  data = size_dat,
  family = match.fun(opt$size$family %||% "negbinomial")(
    opt$size$link %||% "log"
  ),
  prior = prior_obj,
  chains = chains_val,
  iter = iter_val,
  warmup = warmup_val,
  cores = cores_val,
  seed = opt$seed,
  init = brms_cfg$init %||% 0.1,
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
# ... and the estimand stamp, which carries `log_pop_centres`. Under
# `pop_elasticity: free` the simulator MUST re-centre its own frame units on
# these training centres (`.recover_log_pop_centres()`, utils/simulator_core.R);
# brms keeps only formula variables in `fit$data`, so without this the centres
# have to be reverse-engineered from log(pop) - log_pop_c.
attr(brms_fit, "size_estimand") <- size_estimand_stamp

# 9. Save fit + CB metadata ----------------------------------------------------

attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved size_fit -> ", fit_file)

# Layer 4 (artefact chain, see pop_stamp_sidecar_path() in utils.R): record the
# population THIS fit was actually trained on. Safe unconditionally here --
# the skip-reuse case above exits early via quit(), so every path that
# reaches this line just performed a fresh fit against size_dat.
persist_fit_pop_stamp(fit_file, get_pop_provenance(size_dat))

if (length(cb_meta) > 0) {
  meta_file <- make_run_file(opt, kind = "size_cb_meta")
  saveRDS(cb_meta, meta_file)
  message("Saved CB meta -> ", meta_file)
}

# 10. Postprocess RE stats -----------------------------------------------------
# Only attempt if the nested RE groups the adapter expects are present.

n_countries <- dplyr::n_distinct(size_dat$country)
re_groups <- if (n_countries >= 2L) {
  names(brms::ranef(brms_fit))
} else {
  character(0)
}
needed_groups <- if (drop_region) {
  c("country", "country:ADM1")
} else {
  c("AFRO_region", "AFRO_region:country", "AFRO_region:country:ADM1")
}

if (n_countries < 2L) {
  message(
    "Skipping compute_nested_size_re_stats_brms: single-country / partial nesting."
  )
} else if (all(needed_groups %in% re_groups)) {
  re_stats <- compute_nested_size_re_stats_brms(
    brms_fit,
    opt,
    drop_region = drop_region
  )
  re_file <- make_run_file(opt, kind = "size_re_stats")
  saveRDS(re_stats, re_file)
  message("Saved RE stats -> ", re_file)
} else {
  message(
    "Skipping compute_nested_size_re_stats_brms: needed groups not present (",
    paste(setdiff(needed_groups, re_groups), collapse = ", "),
    ")."
  )
}
