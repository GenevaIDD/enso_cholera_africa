# test_mon_filename_roundtrip.R — HARD: make_monthly_run_file() names never
#   collide across load-bearing monthly-track fields.
#
# Usage:
#   Rscript analysis/tests/test_mon_filename_roundtrip.R
#
# Details:
#   The failure mode is COLLISION, silent in the worst way: every fit script
#   early-exits when its output path exists (`if (file.exists(fit_file) &&
#   !isTRUE(opt$overwrite)) quit(status = 0)`), so two configs that differ
#   in a load-bearing field but share a filename do not overwrite or error --
#   the second reports success and leaves the FIRST one's fit on disk under
#   the second one's name. Pure path arithmetic: no data, no fit, no writes.

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

# Minimal resolved-option shape; make_monthly_run_file() reads nothing else.
mk_opt <- function(...) {
  o <- list(
    baseline = "1990_2019",
    country = "all",
    outdir = "generated_data/",
    enso = list(lag_months = 5L),
    grouping = list(enable = FALSE),
    monthly = list(
      backend = "brms_cb",
      family = "negbinomial",
      admin_level = "ADM1-ADM2",
      seasonal_term = "harmonic2",
      phi_re = "constant",
      prev_cases = FALSE,
      covariates = list(),
      random_effects = list(re_adm2_intercept = FALSE),
      dlnm = list(enable = FALSE, defaults = list(), variables = list())
    )
  )
  mods <- list(...)
  for (nm in names(mods)) {
    if (nm %in% names(o$monthly)) {
      o$monthly[[nm]] <- mods[[nm]]
    } else {
      o[[nm]] <- mods[[nm]]
    }
  }
  o
}
base_opt <- mk_opt()
nm <- function(opt, kind = "mon_fit") basename(make_monthly_run_file(opt, kind))

cat("== test_mon_filename_roundtrip ==\n")

# ---------------------------------------------------------------------------
# 1. Shape
# ---------------------------------------------------------------------------
cat("\n-- shape\n")

p <- make_monthly_run_file(base_opt, "mon_fit")
expect(grepl("^generated_data/", p), "path is rooted at opt$outdir")
expect(grepl("\\.rds$", p), "path carries the .rds extension")
expect(
  identical(
    nm(base_opt),
    "outbreak_monthly_brms_cb_fit_base-1990_2019_ctry-all_fam-negbinomial_seas-harmonic2_adm-ADM1-ADM2.rds"
  ),
  "default config resolves to the documented canonical name"
)

expect(
  identical(
    make_monthly_run_file(mk_opt(), "mon_fit"),
    file.path("generated_data/", nm(base_opt))
  ),
  "opt$outdir is honoured verbatim"
)

opt_outdir <- base_opt
opt_outdir$outdir <- "/tmp/x"
expect(
  identical(
    make_monthly_run_file(opt_outdir, "mon_fit"),
    file.path("/tmp/x", nm(base_opt))
  ),
  "opt$outdir override is honoured verbatim"
)

expect(
  inherits(
    try(make_monthly_run_file(base_opt, "size_fit"), silent = TRUE),
    "try-error"
  ),
  "an unknown kind is rejected by match.arg rather than silently accepted"
)

# 05_fit_postprocess.R reaches every track through the one make_run_file() facade
# with a per-track kind prefix. Without the `^mon_` arm, a monthly kind falls
# through to make_occ_run_file() and silently yields an OCCURRENCE path -- which
# exists, so the postprocessor would load the wrong fit rather than fail.
expect(
  identical(
    make_run_file(base_opt, "mon_fit"),
    make_monthly_run_file(base_opt, "mon_fit")
  ),
  "make_run_file() routes mon_* kinds to make_monthly_run_file()"
)
expect(
  all(vapply(
    c(
      "mon_re_stats",
      "mon_coef_summary",
      "mon_loo",
      "mon_kfold",
      "mon_pp_check",
      "mon_cb_meta",
      "mon_cb_surfaces"
    ),
    function(k) {
      identical(make_run_file(base_opt, k), make_monthly_run_file(base_opt, k))
    },
    logical(1)
  )),
  "... for every monthly kind, not just mon_fit"
)

# ---------------------------------------------------------------------------
# 2. Every kind is distinct
# ---------------------------------------------------------------------------
cat("\n-- kinds\n")

kinds <- c(
  "mon_fit",
  "mon_re_stats",
  "mon_coef_summary",
  "mon_loo",
  "mon_kfold",
  "mon_pp_check",
  "mon_cb_meta",
  "mon_cb_surfaces"
)
paths <- vapply(kinds, function(k) nm(base_opt, k), character(1))
expect(
  length(unique(paths)) == length(kinds),
  "the eight artefact kinds resolve to eight distinct filenames"
)
expect(
  all(mapply(
    function(k, f) grepl(sub("^mon_", "", k), f, fixed = TRUE),
    kinds,
    paths
  )),
  "each filename contains its own kind label"
)

# ---------------------------------------------------------------------------
# 3. Load-bearing fields separate
# ---------------------------------------------------------------------------
# Each entry is a config that differs from the default in exactly ONE field.
# None may collide with the default, or with each other.
# The fields that MUST separate: baseline, country/group scope, backend,
# kind, family, seasonal_term, admin_level, the re_adm2_intercept flag,
# prev_cases, the ENSO lag, the DLNM single_lag token, and phi_re.
cat("\n-- separation\n")

variants <- list(
  baseline = mk_opt(baseline = "1980_2009"),
  country = mk_opt(country = "ETH"),
  grouped = mk_opt(
    grouping = list(
      enable = TRUE,
      source = "enso_grouping",
      group = "enso_east"
    )
  ),
  family = mk_opt(family = "poisson"),
  seasonal_term = mk_opt(seasonal_term = "none"),
  admin_level = mk_opt(admin_level = "ADM2"),
  re_adm2_intercept = mk_opt(random_effects = list(re_adm2_intercept = TRUE)),
  prev_cases = mk_opt(prev_cases = TRUE),
  enso = mk_opt(covariates = list(enso = TRUE)),
  single_lag = mk_opt(
    backend = "brms_cb",
    dlnm = list(
      enable = TRUE,
      defaults = list(single_lag = 4L),
      variables = list()
    )
  ),
  phi_re = mk_opt(phi_re = "region_country_adm1")
)
for (v in names(variants)) {
  expect(
    !identical(nm(variants[[v]]), nm(base_opt)),
    sprintf("changing `%s` alone changes the filename", v)
  )
}
all_names <- c(nm(base_opt), vapply(variants, nm, character(1)))
expect(
  length(unique(all_names)) == length(all_names),
  "no two single-field variants collide with each other"
)

# The two fields specific to this track's collision history (mirrors size's
# prev_cases / phi_re tokens) are called out separately so a regression names
# itself rather than hiding in the loop above.
expect(
  grepl("_prevcases(_|\\.rds$)", nm(variants$prev_cases)),
  "prev_cases = TRUE is encoded as _prevcases"
)
expect(
  !grepl("_prevcases", nm(base_opt)),
  "the DEFAULT prev_cases = FALSE earns no token, so existing artefact names are stable"
)

expect(
  grepl("_phire-region-country-adm1", nm(variants$phi_re)),
  "a non-constant phi_re is encoded as _phire-<level>"
)
expect(
  !grepl("_phire-", nm(base_opt)),
  "the DEFAULT phi_re = constant earns no token, so existing artefact names are stable"
)

expect(
  grepl("_sl4\\.rds$", nm(variants$single_lag)),
  "a DLNM single_lag is encoded as _sl<N> (mirrors size/duration)"
)
expect(
  !grepl("_sl[0-9]", nm(base_opt)),
  "no single_lag configured (dlnm disabled) earns no _sl token"
)

# Grouping mode replaces the country token rather than appending to it -- the
# grouped simulator discovers fits by the `grp-<source>-<group>` key.
expect(
  grepl("_grp-enso_grouping-enso_east_", nm(variants$grouped)) &&
    !grepl("ctry-", nm(variants$grouped)),
  "grouping mode swaps ctry-<...> for grp-<source>-<group>"
)

# ---------------------------------------------------------------------------
# 4. Determinism
# ---------------------------------------------------------------------------
cat("\n-- determinism\n")

expect(
  identical(nm(mk_opt()), nm(mk_opt())),
  "the same config resolves to the same path on repeated calls"
)

cat(sprintf(
  "\n== test_mon_filename_roundtrip: %d passed, %d failed ==\n",
  PASS,
  FAIL
))
if (FAIL > 0) {
  quit(status = 1)
}
