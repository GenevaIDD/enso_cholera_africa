# 04_fit_dur_brms_cb.R — brms cessation-hazard duration model with dlnm
#   crossbasis predictors (the deployed duration backend).
#
# Usage:
#   Rscript analysis/04_fit_run.R -C analysis/configs/<batch>/<cfg>.yml
#   Rscript analysis/04_fit_dur_brms_cb.R -C <cfg.yml> --country ETH \
#     --chains 2 --sampling 100 --warmup 100
#
# Details:
#   Activated via duration.backend: brms_cb. Same discrete-time Poisson
#   cessation hazard as 04_fit_dur_brms.R, with covariates entered as a
#   dlnm::crossbasis() or single-lag dlnm::onebasis() (duration.dlnm.variables).
#   Writes make_duration_run_file() kinds dur_fit, dur_cb_meta and
#   dur_re_stats (if full nesting is present). MCMC settings:
#   duration.brms_cb, falling back to duration.brms.
#
# See: analysis/04_fit_dur_brms.R ("Load data": the Holford equivalence)

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
source("analysis/utils/dur_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/plots.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-a", "--admin_level"),
    default = NULL,
    type = "character",
    help = "Admin level (ADM1 | ADM2 | ADM1-ADM2) [config default: ADM1]"
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
    help = "ISO3 code(s), comma-separated, or 'all' [config default: all]"
  ),
  make_option(
    c("--basehaz_duration"),
    default = NULL,
    type = "character",
    help = "Baseline hazard on months-since-onset: factor|ms|none [config default: factor]"
  ),
  make_option(
    c("--basehaz_poisson"),
    default = NULL,
    type = "character",
    help = "Baseline hazard on seasonal time: none|factor|cyclic|seasonal|ms [config default: none]"
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
    "basehaz_duration",
    "basehaz_poisson",
    "sampling",
    "warmup",
    "chains",
    "cores",
    "seed",
    "brms_init",
    "overwrite"
  ),
  model = "duration"
)
stopifnot(!is.null(opt$duration), opt$duration$backend == "brms_cb")
print_options(opt)

dur_cfg <- opt$duration

# brms_cb falls back to brms sub-block for MCMC settings (same pattern as
# 04_fit_occ_brms_cb.R / 04_fit_siz_brms_cb.R).
brms_cfg <- dur_cfg$brms_cb %||% dur_cfg$brms %||% list()

# 3. Early exit if fit already exists ------------------------------------------

fit_file <- make_duration_run_file(opt, kind = "dur_fit")
check_run_file_names(opt, "duration") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("dur_fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data -----------------------------------------------------------------

dur_dat <- load_duration_data(opt)
# Episode record of the prep artefact this fit trains on (assert_episode_key()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(dur_dat, "episode_stamp")

# Anchor sidecar (utils/anchor.R) -- the duration analogue of
# `persist_mean_rate_table(opt, size_dat)` (04_fit_siz_brms_cb.R:110). Gated on
# duration.offset_baserate, so it is inert for every config that predates it.
if (anchor_enabled(opt, "dur")) {
  persist_base_rate_table(opt, "dur", dur_dat)
}

message(sprintf(
  "Duration panel: %d active-month rows, %d episodes, %d cessations",
  nrow(dur_dat),
  dplyr::n_distinct(dur_dat$episode_id),
  sum(dur_dat$y, na.rm = TRUE)
))
if (nrow(dur_dat) == 0L) {
  stop(
    "Duration panel is empty after subsetting — check `country` / `grouping`."
  )
}

if ("climate_missing" %in% names(dur_dat)) {
  n_miss <- sum(dur_dat$climate_missing, na.rm = TRUE)
  if (n_miss > 0L) {
    message(sprintf(
      paste0(
        "  %d row(s) (%.2f%%) lie outside climate coverage and ",
        "will be dropped by brms, carrying %d cessation event(s); ",
        "t_start %s to %s."
      ),
      n_miss,
      100 * n_miss / nrow(dur_dat),
      sum(dur_dat$y[dur_dat$climate_missing], na.rm = TRUE),
      format(min(dur_dat$t_start[dur_dat$climate_missing])),
      format(max(dur_dat$t_start[dur_dat$climate_missing]))
    ))
  }
}

# 5. M-spline bases as fixed columns -------------------------------------------
# Same construction as 04_fit_dur_brms.R — independent of the CB machinery
# below, since basehaz_duration / basehaz_poisson are never crossbasis inputs.

if (identical(dur_cfg$basehaz_duration %||% "factor", "ms")) {
  mso_basis <- splines2::mSpline(
    dur_dat$mso,
    degree = 3,
    df = 5,
    intercept = FALSE
  )
  colnames(mso_basis) <- paste0(.DUR_MS_PREFIX, seq_len(ncol(mso_basis)))
  dur_dat <- dplyr::bind_cols(dur_dat, tibble::as_tibble(mso_basis))
}

if (identical(dur_cfg$basehaz_poisson %||% "none", "ms")) {
  seas_basis <- splines2::mSpline(
    dur_dat$t_start_seas,
    degree = 3,
    df = 5,
    intercept = FALSE
  )
  colnames(seas_basis) <- paste0(.DUR_SEAS_MS_PREFIX, seq_len(ncol(seas_basis)))
  dur_dat <- dplyr::bind_cols(dur_dat, tibble::as_tibble(seas_basis))
}

# 6. Crossbasis predictors -----------------------------------------------------
# The only difference from 04_fit_dur_brms.R: one or more covariates enter as
# a dlnm::crossbasis() (distributed-lag, non-linear in both exposure and lag)
# or dlnm::onebasis() (single-lag, non-linear in exposure only) instead of a
# plain linear term. Typical use: duration.dlnm.variables.enso with
# single_lag: <N> to test a NON-LINEAR effect of ENSO at a fixed lag, as a
# counterpart to the linear enso_lag<N> fixed effect fit by 04_fit_dur_brms.R.
#
# Lag matrices (distributed-lag mode) are read from pre-computed
# <var>_lag1 … <var>_lag{max_lag} columns; single-lag mode reads the single
# pre-computed <var>_lag<N> column (for ENSO: enso_lag<N>, written by
# join_enso_lags() — see analysis/utils/enso.R). If a required column is
# absent, make_lag_matrix() / the single-lag lookup below raises a clear error.
#
# One dlnm::crossbasis() (distributed-lag) or dlnm::onebasis() (single-lag) per
# variable in duration.dlnm.variables. CB matrices are flattened to
# <var>_cb1 … <var>_cbN columns. Verbatim port of the construction loop shared
# by 04_fit_occ_brms_cb.R / 04_fit_siz_brms_cb.R, applied to dur_dat.

dlnm_cfg <- dur_cfg$dlnm %||% list()
cbb <- build_cb_block(dur_dat, dlnm_cfg, " (duration model)")
dur_dat <- cbb$dat
cb_cols <- cbb$cb_cols
cb_meta <- cbb$cb_meta

# 7. Build formula from parts --------------------------------------------------
# build_duration_formula_parts() is called directly (per its own docstring's
# explicit recommendation for callers that must modify terms before assembly,
# as 04_fit_siz_brms_cb.R does on the size side) instead of
# build_duration_formula(), which would re-derive parts without the CB terms.

drop_region <- isTRUE(opt$grouping$enable)
parts <- build_duration_formula_parts(opt, dur_dat)

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
  # per variable via re_adm1_slope: true. Shared plumbing with the occ/size CB
  # scripts. drop_region = TRUE (grouping mode) collapses to country/ADM1 so
  # groups with a single AFRO_region (e.g. enso_south) don't hit the 1-level
  # contrast error. filter_re_terms_for_data() then drops any remaining terms
  # whose grouping column has <2 distinct values (e.g. single-country subsets).
  # "gid" is the duration/occurrence leaf-column convention (size uses "ADM2").
  cb_groups <- re_group_by_level("gid", drop_region = drop_region)
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
      dur_dat,
      keep = .single_country_keep(opt, dur_dat) # one country: keep the country level
    )
    parts$re_terms <- c(parts$re_terms, new_cb_re)
  }
}

form <- brms::bf(stats::as.formula(paste(
  "y ~",
  paste(
    c(
      parts$fixed_terms,
      parts$bh_terms,
      parts$seas_terms,
      parts$re_terms,
      parts$off_terms
    ),
    collapse = " + "
  )
)))
message("Formula:\n  ", deparse1(brms::brmsterms(form)$allvars))
message(
  "  ",
  paste(
    c(
      parts$fixed_terms,
      parts$bh_terms,
      parts$seas_terms,
      parts$re_terms,
      parts$off_terms
    ),
    collapse = " + "
  )
)

# 8. Priors --------------------------------------------------------------------
#
# Identical conditional-prior logic to 04_fit_dur_brms.R: the `b` prior is
# skipped when there are no population-level effects beyond the intercept
# (constant-hazard null variant), and `cor` only applies when correlated RE
# terms are present.

priors <- build_duration_priors(opt, parts, brms_cfg, dat = dur_dat)

# 9. Init ----------------------------------------------------------------------

brms_init_str <- brms_cfg$init %||% "random"
brms_init_val <- if (is.numeric(brms_init_str)) {
  brms_init_str
} else {
  switch(as.character(brms_init_str), random = "random", `0` = 0, {
    warning("Unknown brms_init value '", brms_init_str, "' — using random.")
    "random"
  })
}

# 10. Fit ----------------------------------------------------------------------

mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms_cb")
warmup_val <- mcmc$warmup
sampling_val <- mcmc$sampling
iter_val <- mcmc$iter
refresh_val <- mcmc$refresh

chains_val <- brms_cfg$chains %||% 4L
cores_val <- brms_cfg$cores %||% 4L

threads_per_chain <- brms_cfg$threads_per_chain %||% NULL
threads_arg <- if (
  !is.null(threads_per_chain) && as.integer(threads_per_chain) > 1L
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
  formula = form,
  data = dur_dat,
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
message("Saved dur_fit -> ", fit_file)

if (length(cb_meta) > 0) {
  meta_file <- make_duration_run_file(opt, kind = "dur_cb_meta")
  saveRDS(cb_meta, meta_file)
  message("Saved CB meta -> ", meta_file)
}

# 12. Postprocess --------------------------------------------------------------
#
# compute_nested_surv_re_stats_brms() is REUSED, not forked — same rationale as
# 04_fit_dur_brms.R: the response is Poisson with the identical nested RE
# structure, and the adapter enumerates RE levels from ranef() rather than
# from anything occurrence-specific or CB-specific.

n_countries <- dplyr::n_distinct(dur_dat$country)

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
    "Skipping compute_nested_surv_re_stats_brms: single-country / partial nesting."
  )
} else if (all(needed_groups %in% re_groups)) {
  re_stats <- compute_nested_surv_re_stats_brms(
    brms_fit,
    dur_dat,
    drop_region = drop_region
  )
  re_file <- make_duration_run_file(opt, kind = "dur_re_stats")
  saveRDS(re_stats, re_file)
  message("Saved RE stats -> ", re_file)
} else {
  message(
    "Skipping compute_nested_surv_re_stats_brms: needed groups not present (",
    paste(setdiff(needed_groups, re_groups), collapse = ", "),
    "). Likely a single-country / smoke fit."
  )
}
