# analysis/tests/test_duration_filename_roundtrip.R
#
# Gate on makeDurationRunFile() (analysis/utils/paths.R).
#
# The failure mode is COLLISION, and it is silent in the worst possible way:
# every fit script in this repo early-exits when its output path already exists
# (`if (file.exists(fit_file) && !isTRUE(opt$overwrite)) quit(status = 0)`), so
# two configs that differ in a load-bearing field but resolve to the same
# filename do not overwrite each other and do not error -- the second one
# reports success and leaves the FIRST one's fit on disk under the second one's
# name. Every downstream artefact then describes a model that was never fitted.
#
# On the size track a legacy- and a crosswalk-built mean_rate table could
# collide on one filename (see .sizeEstimandPathToken() in paths.R); the
# duration track gets the gate up front.
#
# The fields that MUST separate: baseline, country/group scope, backend, kind,
# both baseline hazards, admin_level, the ENSO lag, and -- specific to this
# track -- `mso_breaks`, because load_duration_data() rebuilds `mso_band` from
# config at every load, so the band set is a property of the FIT and not of the
# prepared panel.
#
# Pure path arithmetic: no data, no fit, no filesystem writes.
#
# Run: Rscript analysis/tests/test_duration_filename_roundtrip.R

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

# Minimal resolved-option shape; makeDurationRunFile() reads nothing else.
mk_opt <- function(...) {
  o <- list(
    baseline = "1990_2019",
    country  = "all",
    outdir   = "generated_data/",
    enso     = list(lag_months = 11L),
    grouping = list(enable = FALSE),
    duration = list(
      backend          = "brms",
      admin_level      = "ADM1",
      basehaz_duration = "factor",
      basehaz_poisson  = "none",
      mso_breaks       = c(1, 2, 3, 4, 6, 9, Inf),
      covariates       = list(),
      random_effects   = list()
    ))
  mods <- list(...)
  for (nm in names(mods)) {
    if (nm %in% names(o$duration)) o$duration[[nm]] <- mods[[nm]]
    else                            o[[nm]]          <- mods[[nm]]
  }
  o
}
base_opt <- mk_opt()
nm <- function(opt, kind = "dur_fit") basename(makeDurationRunFile(opt, kind))

cat("== test_duration_filename_roundtrip ==\n")

# ---------------------------------------------------------------------------
# 1. Shape
# ---------------------------------------------------------------------------
cat("\n-- shape\n")

p <- makeDurationRunFile(base_opt, "dur_fit")
expect(grepl("^generated_data/", p), "path is rooted at opt$outdir")
expect(grepl("\\.rds$", p),          "path carries the .rds extension")
expect(identical(nm(base_opt),
  "outbreak_duration_brms_fit_base-1990_2019_ctry-all_bhd-factor_adm-ADM1.rds"),
  "default config resolves to the documented canonical name")

expect(identical(makeDurationRunFile(mk_opt(outdir = "/tmp/x"), "dur_fit"),
                 file.path("/tmp/x", nm(base_opt))),
       "opt$outdir is honoured verbatim")

expect(inherits(try(makeDurationRunFile(base_opt, "size_fit"), silent = TRUE),
                "try-error"),
       "an unknown kind is rejected by match.arg rather than silently accepted")

# 05_fit_postprocess.R reaches every track through the one makeRunFile() facade
# with a per-track kind prefix. Without the `^dur_` arm, a duration kind falls
# through to makeOccRunFile() and silently yields an OCCURRENCE path -- which
# exists, so the postprocessor would load the wrong fit rather than fail.
expect(identical(makeRunFile(base_opt, "dur_fit"),
                 makeDurationRunFile(base_opt, "dur_fit")),
       "makeRunFile() routes dur_* kinds to makeDurationRunFile()")
expect(all(vapply(c("dur_re_stats", "dur_coef_summary", "dur_loo", "dur_pp_check"),
                  function(k) identical(makeRunFile(base_opt, k),
                                        makeDurationRunFile(base_opt, k)),
                  logical(1))),
       "... for every duration kind, not just dur_fit")

# ---------------------------------------------------------------------------
# 2. Every kind is distinct
# ---------------------------------------------------------------------------
cat("\n-- kinds\n")

kinds <- c("dur_fit", "dur_re_stats", "dur_coef_summary", "dur_loo", "dur_pp_check")
paths <- vapply(kinds, function(k) nm(base_opt, k), character(1))
expect(length(unique(paths)) == length(kinds),
       "the five artefact kinds resolve to five distinct filenames")
expect(all(mapply(function(k, f) grepl(sub("^dur_", "", k), f, fixed = TRUE),
                  kinds, paths)),
       "each filename contains its own kind label")

# ---------------------------------------------------------------------------
# 3. Load-bearing fields separate
# ---------------------------------------------------------------------------
# Each entry is a config that differs from the default in exactly ONE field.
# None may collide with the default, or with each other.
cat("\n-- separation\n")

variants <- list(
  baseline         = mk_opt(baseline = "1980_2009"),
  country          = mk_opt(country = "ETH"),
  grouped          = mk_opt(grouping = list(enable = TRUE,
                                            source = "enso_grouping",
                                            group  = "enso_east")),
  basehaz_duration = mk_opt(basehaz_duration = "ms"),
  basehaz_none     = mk_opt(basehaz_duration = "none"),
  basehaz_poisson  = mk_opt(basehaz_poisson = "seasonal"),
  mso_breaks       = mk_opt(mso_breaks = c(1, 2, 3, 4, Inf)),
  admin_level      = mk_opt(admin_level = "ADM2"),
  admin_pooled     = mk_opt(admin_level = "ADM1-ADM2"),
  adm2_intercept   = mk_opt(random_effects = list(re_adm2_intercept = TRUE)),
  leaf_adm2        = mk_opt(random_effects = list(re_adm2_leaf_scope = "adm2_only")),
  enso             = mk_opt(covariates = list(enso = TRUE))
)
for (v in names(variants)) {
  expect(!identical(nm(variants[[v]]), nm(base_opt)),
         sprintf("changing `%s` alone changes the filename", v))
}
all_names <- c(nm(base_opt), vapply(variants, nm, character(1)))
expect(length(unique(all_names)) == length(all_names),
       "no two single-field variants collide with each other")

# The one that is specific to this track and easiest to forget. Called out
# separately so a regression names itself rather than hiding in the loop above.
expect(grepl("_mso-1-2-3-4_", nm(variants$mso_breaks)),
       "a non-default mso break set is encoded as _mso-<edges> (finite edges only)")
expect(!grepl("_mso-", nm(base_opt)),
       "the DEFAULT break set earns no token, so existing artefact names are stable")

# Same for the second baseline hazard: `none` is the default and must stay inert.
expect(!grepl("_bhp-", nm(base_opt)) &&
       grepl("_bhp-seasonal_", nm(variants$basehaz_poisson)),
       "basehaz_poisson tokenises only when it is not `none`")

# ADM1-ADM2 rows + ADM2 leaf (monthly's tokens, monthly's order); the defaults
# (re_adm2_intercept FALSE, leaf scope all) earn no token.
expect(grepl("_adm-ADM1-ADM2_adm2i_leaf-adm2", nm(mk_opt(admin_level = "ADM1-ADM2",
         random_effects = list(re_adm2_intercept = TRUE, re_adm2_leaf_scope = "adm2_only")))) &&
       !grepl("_adm2i|_leaf-", nm(base_opt)) &&
       identical(nm(mk_opt(random_effects = list(re_adm2_intercept = FALSE,
                                                 re_adm2_leaf_scope = "all"))), nm(base_opt)),
       "ADM2 intercept / leaf tokenise only when switched on (_adm-ADM1-ADM2_adm2i_leaf-adm2)")

# Grouping mode replaces the country token rather than appending to it -- the
# grouped simulator discovers fits by the `grp-<source>-<group>` key.
expect(grepl("_grp-enso_grouping-enso_east_", nm(variants$grouped)) &&
       !grepl("ctry-", nm(variants$grouped)),
       "grouping mode swaps ctry-<...> for grp-<source>-<group>")

# ---------------------------------------------------------------------------
# 4. Determinism
# ---------------------------------------------------------------------------
cat("\n-- determinism\n")

expect(identical(nm(mk_opt()), nm(mk_opt())),
       "the same config resolves to the same path on repeated calls")
expect(identical(nm(mk_opt(mso_breaks = list(1, 2, 3, 4, 6, 9, Inf))),
                 nm(base_opt)),
       "YAML list-of-scalars mso_breaks resolves identically to the numeric vector")

cat(sprintf("\n== test_duration_filename_roundtrip: %d passed, %d failed ==\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
