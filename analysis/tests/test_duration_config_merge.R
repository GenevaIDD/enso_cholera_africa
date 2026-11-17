# test_duration_config_merge.R — HARD: the `duration:` config block and its
#   formula builders.
#
# Usage:
#   Rscript analysis/tests/test_duration_config_merge.R
#
# Details:
#   Guards two things that fail in opposite ways: the block MERGE (section 1)
#   and the FORMULA (section 3); see the notes above those sections. Fully
#   synthetic: no prepared panel, no climate DBs, no fit. The data frame below
#   carries only the columns the builder inspects.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/dur_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/anchor.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

resolve_yaml <- function(...) {
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(...), tmp)
  suppressMessages(utils::capture.output(
    opt <- make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", tmp)
    )
  ))
  opt
}

# A minimal duration panel: two regions x two countries x two ADM1s, six mso
# bands, everything the builder looks at and nothing else.
mk_dat <- function(n_country = 2L, n_region = 2L) {
  n <- 240L
  set.seed(11)
  data.frame(
    y = rbinom(n, 1L, 0.2),
    mso = rep(c(1, 2, 3, 5, 7, 10), length.out = n),
    mso_band = bin_mso(rep(c(1, 2, 3, 5, 7, 10), length.out = n)),
    log_offset = log(rep(c(17, 31, 30, 28), length.out = n)),
    t_start_seas = rep(seq_len(12), length.out = n),
    month_seas = factor(rep(seq_len(12), length.out = n)),
    AFRO_region = rep(head(c("East", "West"), n_region), length.out = n),
    country = rep(head(c("AAA", "BBB"), n_country), length.out = n),
    ADM1 = rep(c("a1", "a2", "a3", "a4"), length.out = n),
    gid = rep(paste0("g", 1:8), length.out = n),
    mean_precip = rnorm(n),
    dev_precip = rnorm(n),
    mean_temp_daymean = rnorm(n),
    dev_temp_daymean = rnorm(n),
    enso_lag8 = rnorm(n),
    stringsAsFactors = FALSE
  )
}
dat <- mk_dat()

cat("== test_duration_config_merge ==\n")

# ---------------------------------------------------------------------------
# 1. Block merge
# ---------------------------------------------------------------------------
# MERGE. `duration` is a MODEL block: it is filled by .fill_model_block(),
# deliberately NOT registered in .CONFIG_TOPLEVEL_BLOCKS. Getting that wrong is
# silent in both directions: a duration config whose model-detection arm is
# missing resolves every unspecified key to NULL (no error, and the fit then
# runs with basehaz_duration = NULL, i.e. no baseline hazard at all), while
# registering it in the top-level registry would double-merge it and break
# test_config_block_registry.R's completeness assertion.
cat("\n-- merge\n")

defaults <- get_default_options()
opt <- resolve_yaml("duration:", "  backend: brms")

expect(
  !is.null(opt$duration),
  "a `duration:` config resolves a duration block at all"
)
expect(
  identical(opt$duration$basehaz_duration, "factor") &&
    identical(opt$duration$basehaz_poisson, "none"),
  "model-level defaults are filled (basehaz_duration/basehaz_poisson)"
)
expect(
  identical(as.numeric(unlist(opt$duration$mso_breaks)), .DUR_MSO_BREAKS),
  "mso_breaks defaults to .DUR_MSO_BREAKS, terminal Inf included"
)
expect(
  identical(opt$duration$ttlo_reset, "onset"),
  "ttlo_reset defaults to `onset` (the occurrence track's fit-time rule)"
)
expect(
  isTRUE(opt$duration$random_effects$re_adm1_mean_precip_slope) &&
    is.list(opt$duration$random_effects$slopes),
  "random_effects filled, slopes escape hatch materialised as a list"
)
expect(
  !is.null(opt$duration$brms$sampling) && !is.null(opt$duration$brms$prior_sd),
  "the brms sub-block is filled from the active backend's defaults"
)

expect(
  !"duration" %in% names(.CONFIG_TOPLEVEL_BLOCKS),
  "duration is NOT in .CONFIG_TOPLEVEL_BLOCKS (it is a model block)"
)

# YAML wins, siblings still default. The failure this catches is a
# .fill_model_block() nested_keys omission: `covariates` merged as a flat key
# would replace the user's map wholesale rather than filling into it.
opt_p <- resolve_yaml(
  "duration:",
  "  basehaz_duration: ms",
  "  covariates:",
  "    enso: true",
  "  brms:",
  "    chains: 2"
)
expect(
  identical(opt_p$duration$basehaz_duration, "ms"),
  "explicit basehaz_duration overrides the default"
)
expect(
  isTRUE(opt_p$duration$covariates$enso),
  "explicit covariate toggle survives the merge"
)
expect(
  identical(as.integer(opt_p$duration$brms$chains), 2L),
  "explicit brms.chains overrides the default"
)
expect(
  !is.null(opt_p$duration$brms$prior_sd) &&
    identical(opt_p$duration$admin_level, "ADM1"),
  "unspecified siblings (brms.prior_sd, admin_level) still default"
)

# ---------------------------------------------------------------------------
# 2. Validation
# ---------------------------------------------------------------------------
cat("\n-- validation\n")

bad <- function(...) {
  inherits(try(resolve_yaml(...), silent = TRUE), "try-error")
}

expect(
  bad("duration:", "  basehaz_duration: cyclic"),
  "an out-of-set basehaz_duration is rejected (mso takes factor|ms|none)"
)
expect(
  bad("duration:", "  basehaz_poisson: mspline"),
  "an out-of-set basehaz_poisson is rejected"
)
expect(
  bad("duration:", "  mso_breaks: [1, 2, 3]"),
  "mso_breaks not ending in .inf is rejected (bin_mso would error later)"
)
expect(
  bad("duration:", "  mso_breaks: [1, 4, 3, .inf]"),
  "non-monotone mso_breaks is rejected"
)
expect(
  bad("duration:", "  ttlo_reset: TR"),
  "an out-of-set ttlo_reset is rejected"
)
expect(
  !bad(
    "duration:",
    "  basehaz_poisson: seasonal",
    "  mso_breaks: [1, 2, 3, 4, .inf]"
  ),
  "the parsimony break set and a seasonal basehaz are both admissible"
)

# ---------------------------------------------------------------------------
# 3. Formula construction
# ---------------------------------------------------------------------------
# FORMULA. build_duration_formula_parts() decides what the cessation hazard is
# actually conditioned on. Every failure mode here produces a formula that fits
# perfectly well: a dropped `mso_band` gives a constant hazard (the model the
# measured band range rejects), a dropped offset turns a per-day hazard into a
# per-row count, and a mis-set mso_breaks silently re-bins the baseline. None
# of them error.
cat("\n-- formula\n")

rhs_of <- function(opt, d = dat) {
  p <- build_duration_formula_parts(opt, d)
  paste(
    c(p$fixed_terms, p$bh_terms, p$seas_terms, p$re_terms, p$off_terms),
    collapse = " + "
  )
}

r <- rhs_of(opt)
expect(
  grepl("\\bmso_band\\b", r),
  "the default formula carries the mso baseline (`factor` => mso_band)"
)
expect(
  grepl("offset\\(log_offset\\)", r),
  "the TIME offset is present -- without it the Poisson is a count, not a hazard"
)
expect(
  grepl("\\(1 \\| AFRO_region/country/ADM1\\)", r),
  "the nested intercept RE is emitted in slash form"
)
expect(
  all(vapply(
    c("mean_precip", "dev_precip", "mean_temp_daymean", "dev_temp_daymean"),
    function(cv) grepl(sprintf("\\b%s\\b", cv), r),
    logical(1)
  )),
  "convention 2 holds: active slope flags imply the climate fixed effects"
)
expect(
  !grepl("enso_lag", r),
  "no ENSO term when duration.covariates.enso is unset"
)
expect(
  !grepl("t_start_seas|month_seas", r),
  "no seasonal baseline by default (basehaz_poisson = none)"
)

r_none <- rhs_of(resolve_yaml("duration:", "  basehaz_duration: none"))
expect(
  !grepl("\\bmso_band\\b", r_none),
  "basehaz_duration = none removes the mso baseline (constant-hazard variant)"
)
expect(
  grepl("offset\\(log_offset\\)", r_none),
  "... and leaves the offset in place"
)

r_seas <- rhs_of(resolve_yaml("duration:", "  basehaz_poisson: seasonal"))
expect(
  grepl("s\\(t_start_seas, bs = 'cc', k = 8\\)", r_seas) &&
    grepl("\\bmso_band\\b", r_seas),
  "the two baselines stack: seasonal cyclic spline ALONGSIDE mso_band"
)
expect(
  grepl(
    "factor\\(month_seas\\)",
    rhs_of(resolve_yaml("duration:", "  basehaz_poisson: factor"))
  ),
  "basehaz_poisson = factor uses the seasonal month, not the calendar one"
)

r_enso <- rhs_of(resolve_yaml(
  "duration:",
  "  covariates:",
  "    enso: true",
  "enso:",
  "  lag_months: 8"
))
expect(
  grepl("\\benso_lag8\\b", r_enso),
  "covariates.enso = true pulls in enso_lag<N> at the configured lag"
)

# `ms` needs its basis bound onto the frame first. A silent character(0) here
# would drop the baseline entirely, which is why this is an error, not a warning.
opt_ms <- resolve_yaml("duration:", "  basehaz_duration: ms")
expect(
  inherits(
    try(build_duration_formula_parts(opt_ms, dat), silent = TRUE),
    "try-error"
  ),
  "basehaz_duration = 'ms' errors when the mso spline columns are absent"
)
dat_ms <- cbind(dat, msoms1 = rnorm(nrow(dat)), msoms2 = rnorm(nrow(dat)))
p_ms <- build_duration_formula_parts(opt_ms, dat_ms)
expect(
  identical(p_ms$bh_terms, c("msoms1", "msoms2")),
  "... and picks up exactly the msoms* columns once they are bound on"
)

# Grouping mode re-roots the hierarchy at country (region has no meaning inside
# a single fitted group).
r_grp <- rhs_of(resolve_yaml(
  "duration:",
  "  backend: brms",
  "grouping:",
  "  enable: true",
  "  source: enso_grouping",
  "  group: enso_east"
))
expect(
  !grepl("AFRO_region", r_grp) && grepl("\\(1 \\| country/ADM1\\)", r_grp),
  "grouping mode drops AFRO_region and re-roots the RE chain at country"
)

# Data-aware pruning: a single-country subset must not emit country-level REs.
dat_1c <- mk_dat(n_country = 1L, n_region = 1L)
r_1c <- rhs_of(opt, dat_1c)
expect(
  !grepl("AFRO_region:country\\)", r_1c),
  "single-country data: country-level slope REs are pruned, not left to fail in brms"
)

# A subset occupying one mso band carries no baseline contrast.
dat_1b <- dat
dat_1b$mso <- 3
dat_1b$mso_band <- bin_mso(rep(3, nrow(dat_1b)))
expect(
  !grepl("\\bmso_band\\b", suppressMessages(rhs_of(opt, dat_1b))),
  "single-band data: mso_band is dropped rather than left rank-deficient"
)

# ---------------------------------------------------------------------------
# 4. mso_breaks reaches the binning, not just the config
# ---------------------------------------------------------------------------
# The band set is the baseline hazard's whole parameterisation and is rebuilt
# from config at load time. If load_duration_data() ignored it, a parsimony run
# would silently fit the default six bands under the parsimony filename.
cat("\n-- mso_breaks round-trip\n")

opt_brk <- resolve_yaml("duration:", "  mso_breaks: [1, 2, 3, 4, .inf]")
b <- bin_mso(
  c(1, 2, 3, 5, 9),
  breaks = as.numeric(unlist(opt_brk$duration$mso_breaks))
)
expect(
  identical(levels(b), c("1", "2", "3", "4+")),
  "config mso_breaks drive bin_mso's level set"
)
expect(
  identical(as.character(b), c("1", "2", "3", "4+", "4+")),
  "... and its assignment"
)

# ---------------------------------------------------------------------------
# 5. .fill_model_block() sibling-inheritance precedence for brms_cb
# ---------------------------------------------------------------------------
# re_uncorrelated must resolve from block[[be]] (explicit) >
# block[[plain]] (sibling) > defaults[[be]] (canonical), so
# make_batch_configs.R's strip_backend_siblings() can safely drop the plain
# `brms:` sibling from a brms_cb config once its value is merged forward.
cat(
  "\n-- .fill_model_block(): duration/brms_cb inherits re_uncorrelated from the plain brms sibling\n"
)

opt5a <- resolve_yaml(
  "duration:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true"
)
expect(
  isTRUE(opt5a$duration$brms_cb$re_uncorrelated),
  "value set only on the plain brms sibling reaches duration$brms_cb$re_uncorrelated"
)

opt5b <- resolve_yaml(
  "duration:",
  "  backend: brms_cb",
  "  brms:",
  "    re_uncorrelated: true",
  "  brms_cb:",
  "    re_uncorrelated: false"
)
expect(
  identical(opt5b$duration$brms_cb$re_uncorrelated, FALSE),
  "explicit brms_cb$re_uncorrelated wins over a conflicting plain-sibling value"
)

cat(sprintf(
  "\n== test_duration_config_merge: %d passed, %d failed ==\n",
  PASS,
  FAIL
))
if (FAIL > 0) {
  quit(status = 1)
}
