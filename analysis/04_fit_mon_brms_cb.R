# 04_fit_mon_brms_cb.R — brms monthly case-count (within-episode, per calendar
#   month) negative-binomial model with dlnm crossbasis predictors.
#
# Usage:
#   Rscript analysis/04_fit_run.R -C analysis/configs/<batch>/<cfg>.yml
#   Rscript analysis/04_fit_mon_brms_cb.R -C <cfg.yml> --country ETH \
#     --chains 2 --sampling 100 --warmup 100
#
# Details:
#   Activated via monthly.backend: brms_cb (the deployed monthly backend).
#   One crossbasis / onebasis per variable in monthly.dlnm.variables, built by
#   build_cb_block() (utils/dlnm.R) as on the other tracks. Writes
#   make_monthly_run_file() kinds mon_fit and mon_cb_meta; RE stats are
#   computed in 05_fit_postprocess.R. MCMC settings: monthly.brms_cb, falling
#   back to monthly.brms.

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(dlnm)
  library(lubridate) # %m-% -- monthly.prev_cases: yes (mon_data.R:67)
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
source("analysis/utils/dur_data.R")
source("analysis/utils/mon_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
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
    help = "Admin level (ADM1 | ADM2 | ADM1-ADM2) [monthly default: ADM1-ADM2]"
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
  model = "monthly"
)
stopifnot(!is.null(opt$monthly), opt$monthly$backend == "brms_cb")
print_options(opt)

# brms_cb falls back to brms sub-block for MCMC settings (same pattern as
# 04_fit_occ_brms_cb.R / 04_fit_siz_brms_cb.R / 04_fit_dur_brms_cb.R).
brms_cfg <- opt$monthly$brms_cb %||% opt$monthly$brms %||% list()

# 3. Early exit if fit already exists ------------------------------------------

fit_file <- make_run_file(opt, kind = "mon_fit")
check_run_file_names(opt, "monthly") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("mon_fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data -----------------------------------------------------------------

mon_dat <- load_monthly_data(opt)
# Episode record of the prep artefact this fit trains on (assert_episode_key()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(mon_dat, "episode_stamp")

# Capture the estimand stamp NOW, before build_cb_block() below: that helper
# rebuilds `mon_dat` through dplyr verbs, which drop non-standard attributes,
# so reading the stamp at save time would silently yield NULL -- taking the P4
# `log_pop_centres` with it. Exactly the trap 04_fit_siz_brms_cb.R hit.
mon_estimand_stamp <- get_monthly_estimand(mon_dat)

# Anchor sidecar (utils/anchor.R) -- the monthly analogue of
# `persist_mean_rate_table(opt, size_dat)` (04_fit_siz_brms_cb.R). Gated on
# monthly.offset_baserate, so it is inert for every config that predates it.
if (anchor_enabled(opt, "mon")) {
  persist_base_rate_table(opt, "mon", mon_dat)
}

message(sprintf(
  "Monthly panel: %d rows, %d episodes, %d total cases",
  nrow(mon_dat),
  dplyr::n_distinct(mon_dat$episode_id),
  round(sum(mon_dat$cases_month, na.rm = TRUE))
))
if (nrow(mon_dat) == 0L) {
  stop(
    "Monthly panel is empty after subsetting — check `country` / `grouping`."
  )
}

# 5. Crossbasis predictors -----------------------------------------------------
#
#   cases_month_{i,m} ~ NegBin(mu_{i,m}, shape)
#   mu_{i,m} = exp( admin_level + s(enso_lag<N>_m) + s(mso_m) + [seasonal]
#                   + (RE | AFRO_region/country/ADM1/gid) )
#              * pop_i * exposure_days_{i,m}
#
# The reference config carries mso and enso as single-lag onebasis terms
# (single_lag: 0 / 4).
#
# One dlnm::crossbasis() (distributed-lag) or dlnm::onebasis() (single-lag) per
# variable in monthly.dlnm.variables. The loop shared (build_cb_block())
# with the occ/size/duration _cb scripts, applied to mon_dat.

dlnm_cfg <- opt$monthly$dlnm %||% list()
cbb <- build_cb_block(mon_dat, dlnm_cfg, " (monthly model)")
mon_dat <- cbb$dat
cb_cols <- cbb$cb_cols
cb_meta <- cbb$cb_meta

# 6. Build formula from parts --------------------------------------------------
# build_monthly_formula_parts() is called directly (as the size/duration _cb
# scripts do) instead of build_monthly_formula(), which would re-derive parts
# without the CB terms.

drop_region <- isTRUE(opt$grouping$enable)
parts <- build_monthly_formula_parts(opt, mon_dat)

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
  # per variable via re_adm1_slope: true. Shared plumbing with the occ/size/dur
  # CB scripts. drop_region = TRUE (grouping mode) collapses to country/ADM1 so
  # groups with a single AFRO_region (e.g. enso_south) don't hit the 1-level
  # contrast error. filter_re_terms_for_data() then drops any remaining terms
  # whose grouping column has <2 distinct values. "gid" is the duration/
  # occurrence leaf-column convention (size uses "ADM2").
  #
  # Monthly-only extra knob: re_slope_levels picks WHICH of those nesting levels
  # a variable's slopes are emitted at (the occ/size/dur scripts always take
  # nested_slope_re_terms()'s three-level default). mso ships as
  # [region, country] — dropping its country:ADM1 slopes, which are unidentified
  # on monthly per-ADM1 data.
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
    # nested_slope_re_terms() silently skips a level it can't find, so an
    # unrecognised name here would quietly delete every slope for this variable.
    re_levels <- as.character(unlist(
      var_spec$re_slope_levels %||%
        dlnm_defaults$re_slope_levels %||%
        c("region", "country", "adm1")
    ))
    bad_levels <- setdiff(re_levels, names(cb_groups))
    if (length(bad_levels) > 0L) {
      stop(
        "monthly.dlnm variable '",
        v,
        "': unknown re_slope_levels entry/entries ",
        paste(bad_levels, collapse = ", "),
        "; allowed: ",
        paste(names(cb_groups), collapse = ", ")
      )
    }
    new_cb_re <- filter_re_terms_for_data(
      nested_slope_re_terms(v_cb_names, cb_groups, re_sep, levels = re_levels),
      mon_dat,
      keep = .single_country_keep(opt, mon_dat) # one country: keep the country level
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
form_str <- paste("cases_month ~", rhs)
bf_main <- brms::bf(stats::as.formula(form_str))

phi_term <- .phi_re_term(
  phi_re,
  drop_region,
  mon_dat,
  keep = .single_country_keep(opt, mon_dat)
)
bf_obj <- if (length(phi_term) == 0L) {
  bf_main
} else {
  bf_main + brms::lf(stats::as.formula(paste("shape ~ 1 +", phi_term)))
}

message("Formula:\n  ", deparse1(bf_obj$formula))

# 7. Priors --------------------------------------------------------------------

prior_obj <- build_monthly_priors(opt, brms_cfg, parts)

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
  data = mon_dat,
  family = match.fun(opt$monthly$family %||% "negbinomial")(
    opt$monthly$link %||% "log"
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
  save_pars = brms::save_pars(all = TRUE),
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

# ...and the estimand stamp captured before build_cb_block(). `opt` alone is NOT
# sufficient: the stamp carries `log_pop_centres`, the P4 pivots the simulator
# must re-centre frame units on.
attr(brms_fit, "mon_estimand") <- mon_estimand_stamp

# 9. Save fit + CB metadata ----------------------------------------------------
# No persist_fit_pop_stamp() — the panel carries a prep-policy stamp
# (03_prepare_mon.R), not a pop_provenance attribute; 04_fit_dur_brms_cb.R
# stamps nothing either (same rationale: no pop-crosswalk apparatus here).
#
# RE stats / coef summary are NOT computed inline — they are centralised in
# 05_fit_postprocess.R (unlike the size/duration _cb scripts), because the
# monthly intercept RE nests one level deeper (…/ADM1/gid) than the adapter
# functions' 3-level (AFRO_region/country/ADM1) expectation.

attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved mon_fit -> ", fit_file)

if (length(cb_meta) > 0) {
  meta_file <- make_run_file(opt, kind = "mon_cb_meta")
  saveRDS(cb_meta, meta_file)
  message("Saved CB meta -> ", meta_file)
}
