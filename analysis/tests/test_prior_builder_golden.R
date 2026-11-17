# test_prior_builder_golden.R — HARD: golden characterization of each
#   track's brms prior construction (analysis/utils/priors.R).
#
# Usage:
#   Rscript analysis/tests/test_prior_builder_golden.R
#   Rscript analysis/tests/test_prior_builder_golden.R --capture
#
# Details:
#   build_occ_priors() / build_duration_priors() / build_size_priors() /
#   build_monthly_priors() and the shared .PRIOR_DEFAULTS table. Asserts the
#   current output equals tests/golden/prior_builder_golden.rds EXACTLY
#   (identical()), not that a prior is *correct*; --capture writes the
#   golden (review the diff and commit the RDS). Sections [4]-[6] add
#   positive assertions. Pure config + string assembly -- no fit, no data/
#   access, no network.

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(brms)
  library(optparse)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/priors.R") # brms prior builders
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cli <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--capture",
      action = "store_true",
      default = FALSE,
      help = "Run the pre-refactor reference implementations and (re)write the golden file"
    )
  )
))
USE_NEW <- !isTRUE(cli$capture)

GOLDEN_DIR <- "analysis/tests/golden"
GOLDEN_PATH <- file.path(GOLDEN_DIR, "prior_builder_golden.rds")

cat(
  "# prior builder golden characterization (mode:",
  if (USE_NEW) "compare" else "capture",
  ")\n"
)

# The reference implementations .ref_occ_priors() / .ref_dur_priors() (one
# unscoped `sd` prior, no tiers) stay in this file as inert documentation; no
# track routes through them. The golden pins the two group-scoped sd tiers
# (region/country vs. ADM1, via .re_group_tiers()/.tiered_sd_priors()) and the
# class="sds" spline wiggliness prior with its gate (.has_smooth_terms()).
#
# --- Reference implementations (frozen) ---------------------------------------

#' .ref_occ_priors — copied from 04_fit_occ_brms.R:158-164 (== 04_fit_occ_brms_cb.R:210-216)
.ref_occ_priors <- function(terms, brms_cfg) {
  has_cor_params <- !isTRUE(brms_cfg$re_uncorrelated) &&
    any(grepl("\\+.*\\|", terms$re_terms))
  priors <- c(
    prior(normal(0, 1.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 2.5), class = "sd")
  )
  if (has_cor_params) {
    priors <- c(priors, prior(lkj(2), class = "cor"))
  }
  priors
}

#' .ref_dur_priors — copied from 04_fit_dur_brms.R:201-215 (== 04_fit_dur_brms_cb.R:238-252)
.ref_dur_priors <- function(parts, brms_cfg) {
  has_b <- length(c(
    setdiff(parts$fixed_terms, "1"),
    parts$bh_terms,
    parts$seas_terms
  )) >
    0L
  has_cor_params <- !isTRUE(brms_cfg$re_uncorrelated) &&
    any(grepl("\\+.*\\|", parts$re_terms))

  priors <- brms::set_prior(
    brms_cfg$prior_intercept %||% "normal(0, 1.5)",
    class = "Intercept"
  )
  if (has_b) {
    priors <- c(
      priors,
      brms::set_prior(brms_cfg$prior_b %||% "normal(0, 0.5)", class = "b")
    )
  }
  if (length(parts$re_terms)) {
    priors <- c(
      priors,
      brms::set_prior(brms_cfg$prior_sd %||% "normal(0, 2.5)", class = "sd")
    )
  }
  if (has_cor_params) {
    priors <- c(priors, brms::set_prior("lkj(2)", class = "cor"))
  }
  priors
}

# occ/dur, like size/monthly below, are called directly in BOTH modes (see
# the note above the frozen reference implementations).
occ_priors_fn <- function(opt, terms, brms_cfg) {
  build_occ_priors(opt, terms, brms_cfg)
}
dur_priors_fn <- function(opt, parts, brms_cfg) {
  build_duration_priors(opt, parts, brms_cfg, dat = AMPLE_DAT)
}
# build_size_priors()/build_monthly_priors() are called directly in both
# modes.

# --- Synthetic data frame for size/duration/monthly data-aware filters -------
# Mirrors test_formula_builder_golden.R's AMPLE_DAT (ample distinct levels
# everywhere, so filter_fixed_factor_terms_for_data()/filter_re_terms_for_data()
# are a no-op) -- priors depend on term SHAPE, not term-filtering behaviour,
# which is already gated by test_formula_builder_golden.R itself.
.make_dat <- function(n, region, country, adm1, adm2_gid, ttlo, mso_band) {
  d <- data.frame(
    AFRO_region = region,
    country = country,
    ADM1 = adm1,
    ADM2 = adm2_gid,
    gid = adm2_gid,
    time_to_last_outbreak = factor(ttlo),
    mso_band = mso_band,
    stringsAsFactors = FALSE
  )
  ms_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(ms_cols) <- paste0("msoms", 1:5)
  seas_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(seas_cols) <- paste0("seasms", 1:5)
  dplyr::bind_cols(d, ms_cols, seas_cols)
}

AMPLE_DAT <- .make_dat(
  n = 16,
  region = rep(c("East", "West"), each = 8),
  country = rep(c("KEN", "TZA", "NGA", "GHA"), each = 4),
  adm1 = rep(paste0("A", 1:8), each = 2),
  adm2_gid = paste0("G", 1:16),
  ttlo = rep(c("[0,3)", "[3,6)"), length.out = 16),
  mso_band = rep(c(1, 2, 3), length.out = 16)
)

# The brms_cb prior lookup (build_size_priors()/build_monthly_priors() read
# the brms_cb sub-block from a `brms_cb:` call site) is exercised, not
# assumed: for every real size/monthly config the brms_cfg passed in mirrors
# exactly what the deployed 04_fit_<track>_brms[_cb].R script resolves
# (opt$<track>$brms for backend brms; opt$<track>$brms_cb %||% opt$<track>$brms
# for backend brms_cb — see the `brms_cfg <-` lines in those four scripts).
# Occurrence/duration resolve brms_cb in the fit script's own scope, so there
# is no asymmetry to exercise there.
#
# Scoping note: formula parts/terms below are built via the plain formula
# builders (build_formula_terms() / build_size_formula_parts() / etc.),
# BEFORE any DLNM crossbasis columns/RE terms are spliced in -- identical
# scoping to test_formula_builder_golden.R, which characterizes the same
# builders the same way. CB slope RE terms are always the single-covariate
# "(col - 1 | group)" form (nested_slope_re_terms(), utils/formula.R) and can
# never trip the has_cor "+.*|" regex, so this scoping does not under-test the
# cor gate for brms_cb configs.
#
# --- 1. Real configs ----------------------------------------------------------
cat(
  "\n[1] capturing priors for every real config with a brms/brms_cb backend\n"
)

resolve_config_opt <- function(path) {
  suppressWarnings(suppressMessages(
    make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", path)
    )
  ))
}

# forecast_horizons.yml / hindcast_horizons.yml are generated ENSO-horizon
# bookkeeping artefacts, not 04_fit_*-style per-track configs; they throw
# "subscript out of bounds" in make_options_from_config() rather than resolving
# to "no model detected". Phase-09/10 results configs and deploy manifests are
# pipeline-orchestration configs, not model configs, and a results config with
# results.deploy_config set additionally requires analysis/utils/deploy_core.R
# to be sourced first (a deliberate guard) which this test does not load.
# Exclude all of them rather than miscounting them as resolve errors.
NON_MODEL_CONFIGS <- c(
  "manifest.yml",
  "forecast_horizons.yml",
  "hindcast_horizons.yml"
)
NON_MODEL_DIRS <- c("results", "deploy")
config_files <- sort(Filter(
  function(f) {
    !(basename(f) %in% NON_MODEL_CONFIGS) &&
      # run-specific deploy manifests, e.g. forecast_horizons_prod_run.yml,
      # forecast_attribution_prod_run.yml (--stage gen_attrib)
      !grepl(
        "^((forecast|hindcast|event_attribution)_horizons|forecast_attribution).*\\.yml$",
        basename(f)
      ) &&
      !(basename(dirname(f)) %in% NON_MODEL_DIRS) &&
      # the synthetic run's generated batches (tools/run_synthetic_pipeline.sh)
      !grepl("^synth_run[0-9]+_", basename(dirname(f)))
  },
  list.files(
    "analysis/configs",
    pattern = "\\.yml$",
    recursive = TRUE,
    full.names = TRUE
  )
))

real_occ <- list()
real_dur <- list()
real_size <- list()
real_mon <- list()
# Expected class="sds" gate per real config, resolved from the config itself
# rather than read back off the emitted priors, so Section [5] compares two
# independent derivations. Deliberately NOT part of `current` below: the golden
# file's payload stays exactly the four prior lists plus the fixtures.
smooth_occ <- list()
smooth_dur <- list()
smooth_size <- list()
smooth_mon <- list()
n_occ <- 0L
n_dur <- 0L
n_size <- 0L
n_mon <- 0L
n_skipped <- 0L
n_err <- 0L

# Section [6]: every "prior_sd_superseded" warning raised while a real
# config's priors are built, keyed by track::config (muffled here, asserted
# absent there), and the number of real brms blocks that set both tier keys
# AND the legacy prior_sd -- the case the warning inspects. Pure observation:
# the builders' return values (the golden payload) are untouched.
SD_SUPERSEDED <- character(0)
N_TIERS_AND_LEGACY <- 0L
.sd_watch <- function(key, brms_block, expr) {
  b <- brms_block %||% list()
  if (
    !is.null(b[["prior_sd"]]) &&
      !is.null(b[["prior_sd_region_country"]]) &&
      !is.null(b[["prior_sd_adm1"]])
  ) {
    N_TIERS_AND_LEGACY <<- N_TIERS_AND_LEGACY + 1L
  }
  withCallingHandlers(expr, prior_sd_superseded = function(w) {
    SD_SUPERSEDED[[key]] <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
}

for (path in config_files) {
  opt <- tryCatch(resolve_config_opt(path), error = function(e) NULL)
  if (is.null(opt)) {
    n_err <- n_err + 1L
    next
  }

  key <- paste0("real::", path)

  if (
    !is.null(opt$occurrence) && opt$occurrence$backend %in% c("brms", "brms_cb")
  ) {
    n_occ <- n_occ + 1L
    real_occ[[key]] <- tryCatch(
      {
        terms <- suppressWarnings(suppressMessages(
          build_formula_terms(
            opt,
            n_country = 1L,
            n_region = 1L,
            dat = AMPLE_DAT
          )
        ))
        brms_cfg <- if (identical(opt$occurrence$backend, "brms_cb")) {
          opt$occurrence$brms_cb %||% opt$occurrence$brms
        } else {
          opt$occurrence$brms
        }
        # Independent of the builder: mirrors 04_fit_occ_brms.R:148-153's
        # baseline_block switch, which is where the smooth actually enters.
        smooth_occ[[key]] <- identical(.occ_basehaz(brms_cfg), "cyclic")
        .sd_watch(
          paste0("occ::", path),
          brms_cfg,
          occ_priors_fn(opt, terms, brms_cfg)
        )
      },
      error = function(e) list(.error = conditionMessage(e))
    )
  }

  if (
    !is.null(opt$duration) && opt$duration$backend %in% c("brms", "brms_cb")
  ) {
    n_dur <- n_dur + 1L
    real_dur[[key]] <- tryCatch(
      {
        parts <- suppressWarnings(suppressMessages(
          build_duration_formula_parts(opt, AMPLE_DAT)
        ))
        dur_cfg <- opt$duration
        brms_cfg <- if (identical(dur_cfg$backend, "brms_cb")) {
          dur_cfg$brms_cb %||% dur_cfg$brms
        } else {
          dur_cfg$brms
        }
        # config.R's duration seasonal-basehaz switch: cyclic|seasonal are the
        # two spellings that resolve to s(t_start_seas, bs = 'cc', k = 8).
        smooth_dur[[key]] <- (dur_cfg$basehaz_poisson %||% "none") %in%
          c("cyclic", "seasonal")
        .sd_watch(
          paste0("dur::", path),
          brms_cfg,
          dur_priors_fn(opt, parts, brms_cfg)
        )
      },
      error = function(e) list(.error = conditionMessage(e))
    )
  }

  if (!is.null(opt$size) && opt$size$backend %in% c("brms", "brms_cb")) {
    n_size <- n_size + 1L
    real_size[[key]] <- tryCatch(
      {
        # .seasonal_terms(): "cyclic" is the only value yielding s(month, bs='cc').
        smooth_size[[key]] <- identical(
          opt$size$seasonal_term %||% "none",
          "cyclic"
        )
        if (identical(opt$size$backend, "brms_cb")) {
          .sd_watch(
            paste0("size::", path),
            opt$size$brms_cb %||% opt$size$brms,
            build_size_priors(opt, opt$size$brms_cb %||% opt$size$brms)
          )
        } else {
          .sd_watch(
            paste0("size::", path),
            opt$size$brms,
            build_size_priors(opt)
          )
        }
      },
      error = function(e) list(.error = conditionMessage(e))
    )
  }

  if (!is.null(opt$monthly) && opt$monthly$backend %in% c("brms", "brms_cb")) {
    n_mon <- n_mon + 1L
    real_mon[[key]] <- tryCatch(
      {
        smooth_mon[[key]] <- identical(
          opt$monthly$seasonal_term %||% "none",
          "cyclic"
        )
        if (identical(opt$monthly$backend, "brms_cb")) {
          .sd_watch(
            paste0("mon::", path),
            opt$monthly$brms_cb %||% opt$monthly$brms,
            build_monthly_priors(opt, opt$monthly$brms_cb %||% opt$monthly$brms)
          )
        } else {
          .sd_watch(
            paste0("mon::", path),
            opt$monthly$brms,
            build_monthly_priors(opt)
          )
        }
      },
      error = function(e) list(.error = conditionMessage(e))
    )
  }
}

cat(sprintf(
  "  resolved: occurrence=%d duration=%d size=%d monthly=%d | resolve errors=%d\n",
  n_occ,
  n_dur,
  n_size,
  n_mon,
  n_err
))
expect(
  n_err == 0L,
  sprintf("no config failed to resolve (got %d resolve errors)", n_err)
)
expect(
  n_occ + n_dur + n_size + n_mon > 300L,
  sprintf(
    "captured priors for >300 real brms/brms_cb configs (got %d)",
    n_occ + n_dur + n_size + n_mon
  )
)

# --- 2. Synthetic gate-combination fixtures -----------------------------------
cat(
  "\n[2] synthetic gate-combination fixtures (occ/dur) and phi_re x has_slopes (siz/mon)\n"
)

# Four gate combinations, shared by occ and dur:
#   no_re           has_b=T(dur)/forced(occ) has_sd=F has_cor=F
#   re_uncorrelated has_b=T                  has_sd=T has_cor=F  ("||" terms)
#   re_correlated   has_b=T                  has_sd=T has_cor=T  ("|" + "+" terms)
#   constant_null   has_b=F(dur)/forced(occ) has_sd=F has_cor=F  (dur's null variant)
GATE_FIXTURES <- list(
  no_re = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = character(0),
    seas_terms = character(0),
    re_terms = character(0)
  ),
  re_uncorrelated = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = "mso_band",
    seas_terms = character(0),
    re_terms = c("(mean_precip - 1 || AFRO_region)")
  ),
  re_correlated = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = "mso_band",
    seas_terms = character(0),
    re_terms = c("(1 + mean_precip | AFRO_region)")
  ),
  constant_null = list(
    fixed_terms = "1",
    bh_terms = character(0),
    seas_terms = character(0),
    re_terms = character(0)
  ),
  # Isolate each sd-prior tier so the golden file (and
  # Section [4] below) pin .re_group_tiers()'s classification, not just the
  # has_sd on/off gate the four fixtures above already covered.
  region_country_only = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = character(0),
    seas_terms = character(0),
    re_terms = c("(mean_precip - 1 | AFRO_region:country)")
  ),
  adm1_only = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = character(0),
    seas_terms = character(0),
    re_terms = c("(mean_precip - 1 | AFRO_region:country:ADM1)")
  ),
  # sds gate: a seasonal basehaz that resolves to an mgcv smooth.
  # This is the DURATION side of the gate (.has_smooth_terms() scans
  # parts$seas_terms for the other three tracks); the occurrence side is driven
  # by basehaz_poisson in the brms_cfg variants below, since the occurrence
  # baseline block never reaches `parts` at all.
  smooth_seas = list(
    fixed_terms = c("1", "mean_precip"),
    bh_terms = character(0),
    seas_terms = "s(t_start_seas, bs = 'cc', k = 8)",
    re_terms = c("(1 | AFRO_region/country/ADM1)")
  )
)

BRMS_CFG_VARIANTS <- list(
  default = list(),
  re_uncorrelated_flag = list(re_uncorrelated = TRUE),
  custom_priors = list(
    prior_intercept = "normal(-1, 2)",
    prior_b = "normal(0, 1)",
    prior_sd = "cauchy(0, 1)",
    prior_sds = "normal(0, 1)"
  ),
  # Occurrence's sds gate reads basehaz_poisson off the backend block. `default`
  # above leaves it absent, which resolves to the config default "seasonal" ->
  # smooth present, so this variant is the only occurrence fixture WITHOUT a
  # smooth -- the negative case Section [5] needs. (It is inert for duration,
  # whose own seasonal basehaz lives at model level, not in the brms block.)
  basehaz_factor = list(basehaz_poisson = "factor")
)

gate_occ <- list()
gate_dur <- list()
for (fx_nm in names(GATE_FIXTURES)) {
  for (cfg_nm in names(BRMS_CFG_VARIANTS)) {
    key <- paste0(fx_nm, "::", cfg_nm)
    gate_occ[[key]] <- occ_priors_fn(
      list(),
      GATE_FIXTURES[[fx_nm]],
      BRMS_CFG_VARIANTS[[cfg_nm]]
    )
    gate_dur[[key]] <- dur_priors_fn(
      list(),
      GATE_FIXTURES[[fx_nm]],
      BRMS_CFG_VARIANTS[[cfg_nm]]
    )
  }
}
expect(
  length(gate_occ) == length(GATE_FIXTURES) * length(BRMS_CFG_VARIANTS),
  "captured all occ gate x brms_cfg combinations"
)
expect(
  length(gate_dur) == length(GATE_FIXTURES) * length(BRMS_CFG_VARIANTS),
  "captured all dur gate x brms_cfg combinations"
)

# phi_re on/off x has_slopes on/off, for size and monthly.
mk_slope_opt <- function(track, has_slopes, phi_re, seasonal_term = "none") {
  block <- list(
    admin_level = "ADM1",
    brms = list(),
    random_effects = if (has_slopes) {
      list(re_adm1_mean_precip_slope = TRUE)
    } else {
      list()
    },
    covariates = list(enso = FALSE),
    phi_re = phi_re,
    seasonal_term = seasonal_term
  )
  setNames(list(block), track)
}

slope_size <- list()
slope_mon <- list()
for (has_slopes in c(FALSE, TRUE)) {
  for (phi_re in c("constant", "region_country_adm1")) {
    key <- sprintf("slopes-%s::phi-%s", has_slopes, phi_re)
    slope_size[[key]] <- build_size_priors(mk_slope_opt(
      "size",
      has_slopes,
      phi_re
    ))
    slope_mon[[key]] <- build_monthly_priors(mk_slope_opt(
      "monthly",
      has_slopes,
      phi_re
    ))
  }
}
# sds gate for size/monthly: seasonal_term "cyclic" is the only value
# .seasonal_terms() turns into s(month, bs='cc', k=6). Appended as extra keys
# rather than folded in as a third loop axis, so the four keys above (and their
# golden entries) keep their existing names. has_slopes stays FALSE here on
# purpose: a cyclic smooth adds NO class="b" parameter (empty null space), so
# these fixtures must show an `sds` row and no `b` row -- which is also the
# positive control for .prior_config_slopes()'s deliberate "cyclic" exclusion.
for (st in c("cyclic", "harmonic2")) {
  key <- sprintf("slopes-FALSE::phi-constant::seas-%s", st)
  slope_size[[key]] <- build_size_priors(mk_slope_opt(
    "size",
    FALSE,
    "constant",
    st
  ))
  slope_mon[[key]] <- build_monthly_priors(mk_slope_opt(
    "monthly",
    FALSE,
    "constant",
    st
  ))
}
expect(
  length(slope_size) == 6L,
  "captured all size phi_re x has_slopes x seasonal_term combinations"
)
expect(
  length(slope_mon) == 6L,
  "captured all monthly phi_re x has_slopes x seasonal_term combinations"
)

# --- 3. Compare against (or create) the committed golden file -----------------
# One deliberate exception to "current == golden": the gate_occ "custom_priors"
# fixture. The frozen .ref_occ_priors() ignores brms_cfg$prior_*, while
# build_occ_priors() honours it, so this one synthetic combination must differ
# from the reference. No real occurrence config sets prior_intercept, prior_b
# or prior_sd (they set the tier keys prior_sd_region_country / prior_sd_adm1 /
# prior_sds), so real_occ is unaffected. This section carves this one key out
# of the strict identity gate and asserts the override took effect instead.
cat("\n[3] compare against committed golden file\n")

current <- list(
  real_occ = real_occ,
  real_dur = real_dur,
  real_size = real_size,
  real_mon = real_mon,
  gate_occ = gate_occ,
  gate_dur = gate_dur,
  slope_size = slope_size,
  slope_mon = slope_mon
)

if (isTRUE(cli$capture)) {
  dir.create(GOLDEN_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(current, GOLDEN_PATH)
  cat("  Golden file (re)written:", GOLDEN_PATH, "\n")
  cat(
    "  Review it (git diff / git add) like any other change, then re-run without --capture.\n"
  )
  expect(TRUE, "golden file captured")
} else {
  if (!file.exists(GOLDEN_PATH)) {
    fail(sprintf(
      "golden file missing: %s -- run with --capture first",
      GOLDEN_PATH
    ))
  }
  golden <- readRDS(GOLDEN_PATH)

  # Occurrence's "custom_priors" gate fixture (see the note above [3]) is the one
  # combination that MUST diverge from the frozen reference, which ignores
  # brms_cfg$prior_*. Carve it out of the
  # strict identity gate below and check it on its own terms instead.
  KNOWN_DIVERGENT <- list(
    gate_occ = grep("::custom_priors$", names(gate_occ), value = TRUE)
  )

  for (grp in names(current)) {
    expect(
      identical(names(current[[grp]]), names(golden[[grp]])),
      sprintf(
        "%s key set unchanged (current=%d, golden=%d)",
        grp,
        length(current[[grp]]),
        length(golden[[grp]])
      )
    )
    known <- KNOWN_DIVERGENT[[grp]] %||% character(0)
    check_keys <- setdiff(
      intersect(names(current[[grp]]), names(golden[[grp]])),
      known
    )
    mism <- Filter(
      function(k) !identical(current[[grp]][[k]], golden[[grp]][[k]]),
      check_keys
    )
    expect(
      length(mism) == 0L,
      sprintf(
        "%s priors byte-identical to golden (mismatches: %s)",
        grp,
        paste(head(mism, 5), collapse = ", ")
      )
    )
  }

  # Positive check for the carved-out occurrence override case: the NEW
  # builder must actually honour the custom prior_* strings (proving the
  # divergence above is the intended fix, not silent breakage).
  custom <- BRMS_CFG_VARIANTS$custom_priors
  for (key in KNOWN_DIVERGENT$gate_occ %||% character(0)) {
    p <- gate_occ[[key]]
    got_intercept <- p$prior[p$class == "Intercept"]
    got_b <- p$prior[p$class == "b"]
    expect(
      identical(got_intercept, custom$prior_intercept) &&
        (length(got_b) == 0L || identical(got_b, custom$prior_b)),
      sprintf(
        "gate_occ::%s honours brms_cfg$prior_* overrides (consequence 3 fix)",
        key
      )
    )
  }
}

# --- 4. Positive assertions for the region/country vs. ADM1 tier split -------
# A byte-identical golden diff documents THAT something changed but not what
# the new shape actually is. These assertions pin the tiering behaviour
# itself (.re_group_tiers()'s classification), independent of the golden file.
cat("\n[4] region/country vs. ADM1 sd-prior tier split\n")

.sd_groups <- function(p) {
  if (is.data.frame(p)) as.character(p$group[p$class == "sd"]) else character(0)
}

for (trk_nm in c("gate_occ", "gate_dur")) {
  trk <- get(trk_nm)

  rc_groups <- .sd_groups(trk[["region_country_only::default"]])
  expect(
    identical(rc_groups, "AFRO_region:country"),
    sprintf(
      "%s region_country_only fixture: exactly one sd row scoped to AFRO_region:country (got %s)",
      trk_nm,
      paste(rc_groups, collapse = ", ")
    )
  )

  adm1_groups <- .sd_groups(trk[["adm1_only::default"]])
  expect(
    identical(adm1_groups, "AFRO_region:country:ADM1"),
    sprintf(
      "%s adm1_only fixture: exactly one sd row scoped to AFRO_region:country:ADM1 (got %s)",
      trk_nm,
      paste(adm1_groups, collapse = ", ")
    )
  )
}

# A fixture spanning both tiers (real configs with nested intercept + slope REs
# routinely do) must show 2+ distinct sd groups, split across the region/country
# and ADM1 tiers -- proving the split isn't just "one group renamed" but an
# actual partition of the variance components.
multi_tier_hit <- FALSE
for (trk_nm in c("real_occ", "real_dur", "real_size", "real_mon")) {
  trk <- get(trk_nm)
  for (p in trk) {
    groups <- .sd_groups(p)
    if (length(unique(groups)) < 2L) {
      next
    }
    is_adm1 <- grepl("ADM1", groups, fixed = TRUE)
    if (any(is_adm1) && any(!is_adm1)) {
      multi_tier_hit <- TRUE
      break
    }
  }
  if (multi_tier_hit) break
}
expect(
  multi_tier_hit,
  "at least one real config shows both region/country- and ADM1-tier sd rows"
)

# Every occurrence / duration config whose basehaz_poisson resolves to
# cyclic/seasonal, and every size/monthly config with seasonal_term: cyclic,
# carries exactly one `sds` row at brms's own package default
# (student_t(3, 0, 2.5)), so the emitted model is numerically the same as
# leaving it implicit (make_stancode() is byte-identical either way). This
# section asserts the gate both ways, since a byte-identical-to-golden check
# can only say that output is consistent, not correct -- and an over-eager
# gate is not cosmetic: set_prior(class="sds") against a formula with no
# smooth is a hard brms .validate_prior() failure.
#
# --- 5. class="sds" spline wiggliness prior gate -------------------------------
# Same rationale as Section [4]: the golden comparison proves the emitted priors
# changed *consistently*, these assertions prove they changed *correctly*. The
# stakes are higher here than for a wrong-but-harmless prior string -- an sds
# prior emitted against a formula with no smooth is a hard brms
# .validate_prior() error, so a false positive on this gate breaks the fit
# outright. Both directions are therefore asserted, on fixtures and on all ~400
# real configs, against an expectation derived from the config rather than read
# back off the priors.
cat("\n[5] class=\"sds\" smooth-prior gate\n")

.sds_rows <- function(p) {
  if (is.data.frame(p)) {
    as.character(p$prior[p$class == "sds"])
  } else {
    character(0)
  }
}

# -- fixtures: occurrence (gate driven by basehaz_poisson in the brms block) ---
occ_smooth_keys <- grep(
  "::(default|re_uncorrelated_flag)$",
  names(gate_occ),
  value = TRUE
)
expect(
  all(vapply(
    gate_occ[occ_smooth_keys],
    function(p) length(.sds_rows(p)) == 1L,
    logical(1)
  )),
  "gate_occ: every fixture with an absent basehaz_poisson (-> default 'seasonal') emits exactly one sds row"
)
occ_factor_keys <- grep("::basehaz_factor$", names(gate_occ), value = TRUE)
expect(
  all(vapply(
    gate_occ[occ_factor_keys],
    function(p) length(.sds_rows(p)) == 0L,
    logical(1)
  )),
  "gate_occ: basehaz_poisson = 'factor' emits no sds row"
)
expect(
  identical(
    unique(unlist(lapply(gate_occ[occ_smooth_keys], .sds_rows))),
    .PRIOR_DEFAULTS$occurrence$sds
  ),
  sprintf(
    "gate_occ: sds default is brms's own package default (%s)",
    .PRIOR_DEFAULTS$occurrence$sds
  )
)

# -- fixtures: duration (gate driven by the seas_terms string in `parts`) ------
dur_smooth_keys <- grep("^smooth_seas::", names(gate_dur), value = TRUE)
dur_flat_keys <- setdiff(names(gate_dur), dur_smooth_keys)
expect(
  all(vapply(
    gate_dur[dur_smooth_keys],
    function(p) length(.sds_rows(p)) == 1L,
    logical(1)
  )),
  "gate_dur: the smooth_seas fixture emits exactly one sds row"
)
expect(
  all(vapply(
    gate_dur[dur_flat_keys],
    function(p) length(.sds_rows(p)) == 0L,
    logical(1)
  )),
  "gate_dur: every fixture without a smooth seas_term emits no sds row"
)

# -- fixtures: size / monthly (gate driven by seasonal_term) -------------------
for (trk_nm in c("slope_size", "slope_mon")) {
  trk <- get(trk_nm)
  expect(
    length(.sds_rows(trk[["slopes-FALSE::phi-constant::seas-cyclic"]])) == 1L,
    sprintf("%s: seasonal_term = 'cyclic' emits exactly one sds row", trk_nm)
  )
  expect(
    length(.sds_rows(trk[["slopes-FALSE::phi-constant::seas-harmonic2"]])) ==
      0L,
    sprintf("%s: seasonal_term = 'harmonic2' emits no sds row", trk_nm)
  )
  # Pinned: a bs='cc' basis has an empty null space, so brms
  # creates no class="b" coefficient for it -- .prior_config_slopes()'s
  # exclusion of "cyclic" from has_slopes is correct, not an oversight.
  p_cyc <- trk[["slopes-FALSE::phi-constant::seas-cyclic"]]
  expect(
    sum(p_cyc$class == "b") == 0L,
    sprintf(
      "%s: a cyclic seasonal_term adds an sds row but no class='b' row",
      trk_nm
    )
  )
}

# -- config override -----------------------------------------------------------
expect(
  identical(
    .sds_rows(gate_occ[["no_re::custom_priors"]]),
    BRMS_CFG_VARIANTS$custom_priors$prior_sds
  ),
  "gate_occ: brms_cfg$prior_sds override is honoured"
)
expect(
  identical(
    .sds_rows(gate_dur[["smooth_seas::custom_priors"]]),
    BRMS_CFG_VARIANTS$custom_priors$prior_sds
  ),
  "gate_dur: brms_cfg$prior_sds override is honoured"
)

# -- all real configs, both directions ----------------------------------------
for (pair in list(
  c("real_occ", "smooth_occ"),
  c("real_dur", "smooth_dur"),
  c("real_size", "smooth_size"),
  c("real_mon", "smooth_mon")
)) {
  priors_l <- get(pair[1])
  expect_l <- get(pair[2])
  keys <- intersect(names(priors_l), names(expect_l))
  got <- vapply(
    priors_l[keys],
    function(p) length(.sds_rows(p)) > 0L,
    logical(1)
  )
  want <- vapply(expect_l[keys], isTRUE, logical(1))
  mism <- keys[got != want]
  expect(
    length(mism) == 0L,
    sprintf(
      "%s: sds row present iff the config's formula has a smooth (%d/%d with smooth; mismatches: %s)",
      pair[1],
      sum(want),
      length(keys),
      if (length(mism)) {
        paste(head(basename(mism), 3), collapse = ", ")
      } else {
        "none"
      }
    )
  )
}

# --- 6. Legacy prior_sd superseded by both tier keys --------------------------
# .prior_sd_superseded() (utils/priors.R, called inside .prior_strings()) warns
# when both tier keys are set -- so the blanket prior_sd is never read -- AND
# prior_sd is not the track's default. The defaults write prior_sd into the
# configs (e.g. the duration configs' normal(0, 2.5) beside 0.35 tiers), so the
# warning must stay silent on every real config, and the resolved strings must
# not move (the golden comparison above already pins that).
cat("\n[6] legacy prior_sd superseded by both tier keys (E9)\n")

expect(
  length(SD_SUPERSEDED) == 0L,
  sprintf(
    "no real config raises prior_sd_superseded (%d did: %s)",
    length(SD_SUPERSEDED),
    paste(head(names(SD_SUPERSEDED), 3), collapse = ", ")
  )
)
expect(
  N_TIERS_AND_LEGACY > 0L,
  sprintf(
    "the silent case is exercised: %d real brms blocks set both tiers and prior_sd",
    N_TIERS_AND_LEGACY
  )
)

.n_sd_warn <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, prior_sd_superseded = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  n
}
TIERS <- list(
  prior_sd_region_country = "normal(0, 0.35)",
  prior_sd_adm1 = "normal(0, 0.2)"
)

for (trk in names(.PRIOR_DEFAULTS)) {
  stale <- c(TIERS, list(prior_sd = "cauchy(0, 1)"))
  n <- .n_sd_warn(s <- .prior_strings(stale, trk))
  expect(
    n == 1L,
    sprintf(
      "%s: both tiers + prior_sd cauchy(0, 1) -> one warning (got %d)",
      trk,
      n
    )
  )
  expect(
    identical(s$sd_region_country, TIERS$prior_sd_region_country) &&
      identical(s$sd_adm1, TIERS$prior_sd_adm1),
    sprintf("%s: the tiers still win over the warned prior_sd", trk)
  )
  # The default, also with other spacing, is silent.
  dflt <- .PRIOR_DEFAULTS[[trk]]$sd_region_country
  for (v in c(dflt, gsub(" ", "", dflt, fixed = TRUE))) {
    n <- .n_sd_warn(.prior_strings(c(TIERS, list(prior_sd = v)), trk))
    expect(
      n == 0L,
      sprintf(
        "%s: both tiers + prior_sd '%s' (the default) -> silent (got %d)",
        trk,
        v,
        n
      )
    )
  }
  # One tier only: prior_sd still seeds the other tier, so it is not dead.
  n <- .n_sd_warn(.prior_strings(
    list(prior_sd_adm1 = "normal(0, 0.2)", prior_sd = "cauchy(0, 1)"),
    trk
  ))
  expect(
    n == 0L,
    sprintf("%s: one tier + prior_sd cauchy(0, 1) -> silent (got %d)", trk, n)
  )
}

# Through a public builder: exactly one warning per build.
n <- .n_sd_warn(dur_priors_fn(
  list(),
  GATE_FIXTURES$re_uncorrelated,
  c(TIERS, list(prior_sd = "cauchy(0, 1)"))
))
expect(
  n == 1L,
  sprintf(
    "build_duration_priors(): both tiers + prior_sd cauchy(0, 1) -> one warning (got %d)",
    n
  )
)

# The legacy-only custom_priors fixture (prior_sd set, no tier key) is the
# pre-split shape: prior_sd seeds both tiers and is live, so no warning.
for (fx_nm in names(GATE_FIXTURES)) {
  n <- .n_sd_warn(occ_priors_fn(
    list(),
    GATE_FIXTURES[[fx_nm]],
    BRMS_CFG_VARIANTS$custom_priors
  )) +
    .n_sd_warn(dur_priors_fn(
      list(),
      GATE_FIXTURES[[fx_nm]],
      BRMS_CFG_VARIANTS$custom_priors
    ))
  expect(
    n == 0L,
    sprintf(
      "custom_priors (legacy prior_sd only) x %s -> silent (got %d)",
      fx_nm,
      n
    )
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
message("PASS  test_prior_builder_golden")
