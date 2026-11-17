# 04_fit_dur_brms.R — config-driven brms fit of the outbreak-duration track
#   (discrete-time cessation hazard, Poisson with log offset).
#
# Usage:
#   Rscript analysis/04_fit_run.R \
#     -C analysis/configs/dur_smoketest/dur_smoke.yml
#   Rscript analysis/04_fit_dur_brms.R -C <cfg> --country ETH --sampling 200
#   Rscript analysis/04_fit_dur_brms.R -C <cfg> --basehaz_duration none
#
# Details:
#   One row per active episode-month (03_prepare_dur.R); y = 1 on an
#   uncensored episode's last active month. The offset is TIME (days active in
#   the month): dropping it silently turns a per-day hazard into a per-row
#   count. Two configurable baseline hazards: duration.basehaz_duration on
#   months since onset (primary) and duration.basehaz_poisson on seasonal time.
#   Parameters: CLI flag > YAML config (-C) > config_defaults.yml.

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(splines2)
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
stopifnot(!is.null(opt$duration), opt$duration$backend == "brms")
print_options(opt)

dur_cfg <- opt$duration
brms_cfg <- dur_cfg$brms %||% list()

# 3. Early exit if fit already exists ------------------------------------------
#
# NOTE this is exactly the early exit that makes an artefact-name collision
# silent: a config that resolves to an existing name is reported as a success
# and leaves the OTHER config's fit on disk. make_duration_run_file() is gated by
# analysis/tests/test_duration_filename_roundtrip.R for that reason.

fit_file <- make_duration_run_file(opt, kind = "dur_fit")
check_run_file_names(opt, "duration") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("Fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data -----------------------------------------------------------------
#
# Duration = discrete-time cessation hazard, one row per active month per
# episode (03_prepare_dur.R -> utils/dur_data.R); y = 1 on an uncensored
# episode's last active month; log_offset = log(days of the month inside
# [TL, TR]). Fitted as Poisson with a log offset (the Holford form), like
# occurrence, so the occurrence RE postprocessing is reused. Scaling, the
# `mso_band` rebuild from duration.mso_breaks and country subsetting happen in
# load_duration_data().

dur_dat <- load_duration_data(opt)
# Episode record of the prep artefact this fit trains on (assert_episode_key()
# in the loader); stamped onto the fit below so the simulator can refuse a fit
# whose episode definition differs from the artefacts it is paired with.
ep_rec <- attr(dur_dat, "episode_stamp")

# Anchor sidecar (utils/anchor.R) -- the duration analogue of
# `persist_mean_rate_table(opt, size_dat)` (04_fit_siz_brms.R). Gated on
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

# guard_missing_climate_data() KEEPS rows outside CHIRPS/ERA5 coverage and flags
# them (occ_data.R), leaving the filtering to the fit script. brms will
# drop them with a generic "Rows containing NAs were excluded" warning that
# says nothing about how many CESSATIONS went with them — and a dropped `y = 1`
# row is a lost event, not just lost exposure. Report it explicitly.
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
#
# Same construction as the occurrence track's `ms` basehaz
# (04_fit_occ_brms.R), but on each baseline's own clock. The formula
# builder only NAMES these columns (^msoms[0-9]+$ / ^seasms[0-9]+$) and errors
# if they are absent, so the basis must be bound on before it is called.

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

# 6. Build formula from config -------------------------------------------------
#
# The offset is TIME (days active in the month). Dropping it does not error:
# it silently turns a per-day hazard into a per-row count. It is emitted by
# build_duration_formula_parts() and checked by test_duration_config_merge.R.
# Two baseline hazards, configured independently:
#   duration.basehaz_duration = factor | ms | none on `mso` (primary);
#   duration.basehaz_poisson = none | factor | cyclic | seasonal | ms on
#     `t_start_seas` (default none: season-end as a testable alternative).
# build_duration_formula_parts() prunes levels and RE terms the subset cannot
# support, which is what makes a single-country smoke run possible.

parts <- build_duration_formula_parts(opt, dur_dat)
form <- build_duration_formula(opt, dur_dat)
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

# 7. Priors --------------------------------------------------------------------
#
# Strings live in the config (duration.brms.prior_*) so a prior-sensitivity
# batch is a config change, not a code change. The `b` prior is conditional:
# brms errors on a prior for a class with no parameters, and the constant-
# hazard null variant (basehaz_duration: none, no covariates, no slopes) has
# no population-level effects beyond the intercept.

priors <- build_duration_priors(opt, parts, brms_cfg, dat = dur_dat)

# 8. Init ----------------------------------------------------------------------

brms_init_str <- brms_cfg$init %||% "random"
brms_init_val <- if (is.numeric(brms_init_str)) {
  brms_init_str
} else {
  switch(as.character(brms_init_str), random = "random", `0` = 0, {
    warning("Unknown brms_init value '", brms_init_str, "' — using random.")
    "random"
  })
}

# 9. Fit -----------------------------------------------------------------------

mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms")
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

# 10. Save fit -----------------------------------------------------------------
attr(brms_fit, "episode_stamp") <- ep_rec
saveRDS(brms_fit, fit_file)
message("Saved fit -> ", fit_file)

# 11. Postprocess --------------------------------------------------------------
#
# compute_nested_surv_re_stats_brms() is REUSED, not forked: the response is Poisson
# with the identical nested RE structure, and the adapter enumerates RE levels
# from ranef() rather than from anything occurrence-specific. Forking it would
# be the most expensive duplication in this track — do it only if this fails.

n_countries <- dplyr::n_distinct(dur_dat$country)

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
