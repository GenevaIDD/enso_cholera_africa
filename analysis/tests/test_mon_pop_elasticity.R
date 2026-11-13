# analysis/tests/test_mon_pop_elasticity.R
#
# Gate on the monthly track's denominator/elasticity apparatus -- the size
# track's machinery, which guards against a ~76x monthly overestimate in CMR.
#
# The failure this file guards is NOT "the formula string changed" -- that is
# test_formula_builder_golden.R's job. It is the two SILENT ones:
#
#   1. A mon fit is run with `pop_elasticity: free`, the kernel cannot find the
#      training centres, and re-centres frame units on their own mean instead.
#      Every prediction is then rescaled by exp((beta - 1) * shift) with no
#      diagnostic that could detect it. The mon stamp is written under a
#      DIFFERENT attribute name (`mon_estimand`) than size's (`size_estimand`),
#      so a `.recover_log_pop_centres()` that only knew the size name would fall
#      through to the fit$data path -- or, on a fit whose data was rebuilt by a
#      dplyr verb, hard-error at forecast time.
#   2. A mon artefact built under one denominator silently overwrites one built
#      under another, because `file.exists()` early-exit means colliding names
#      neither overwrite nor error (see test_mon_filename_roundtrip.R's header).
#
# Plus the invariant that makes `fixed` and `free` NESTED rather than merely
# similar: log_pop_c is centred WITHIN admin_level, so a coefficient of 0
# reproduces the offset-only model EXACTLY.
#
# Companions: test_formula_builder_golden.R (term assembly + hard-error path),
# test_mon_filename_roundtrip.R (the other filename-separation fields),
# test_monthly_offset_wellformed.R (offset wiring), gates_estimand_coherence.R
# (the size-side frame gates).
#
# Pure config/attribute/vector arithmetic: no data/ access, no fit, no
# filesystem writes.
#
# Run: Rscript analysis/tests/test_mon_pop_elasticity.R

suppressPackageStartupMessages({
  library(dplyr)
})

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/siz_data.R")        # .addLogPopCentred()
source("analysis/utils/dlnm.R")
source("analysis/utils/simulator_core.R")  # .recover_log_pop_centres(), .adminLevelFromGid()

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# monthly pop-elasticity / estimand apparatus\n")

# ---------------------------------------------------------------------------
# [1] P0 -- resolveMonthlyEstimand() defaults and validation
# ---------------------------------------------------------------------------
cat("\n[1] P0: monthly estimand resolution\n")

# Absent block must resolve to TODAY's deployed mon behaviour, bit-for-bit.
es0 <- resolveMonthlyEstimand(list(monthly = list()))
expect(identical(es0$denominator, "footprint") &&
         identical(es0$rate_denominator, "footprint") &&
         identical(es0$underreporting_multiplier, 1.0) &&
         identical(es0$pop_elasticity, "fixed"),
       "absent monthly.estimand resolves to footprint/footprint/1/fixed")
expect(identical(es0$track, "monthly"),
       "monthly estimand records track = 'monthly' (not 'size')")

# The two tracks must not share a resolver result: a mon config's keys must not
# be read out of opt$size, nor vice versa.
opt_split <- list(size    = list(estimand = list(denominator = "gid")),
                  monthly = list())
expect(identical(resolveSizeEstimand(opt_split)$denominator, "gid") &&
         identical(resolveMonthlyEstimand(opt_split)$denominator, "footprint"),
       "size.estimand does not leak into the monthly resolution")

# denominator != rate_denominator is the whole point of the block: hard error.
err <- tryCatch({
  resolveMonthlyEstimand(list(monthly = list(estimand = list(
    denominator = "gid", rate_denominator = "footprint")))); NULL
}, error = function(e) conditionMessage(e))
expect(!is.null(err) && grepl("rate_denominator", err, fixed = TRUE),
       "monthly denominator != rate_denominator is a hard error")
expect(!is.null(err) && grepl("monthly", err, fixed = TRUE),
       "the mismatch error names the `monthly` block (not `size`)")

# A typo must NOT silently fall back to the default.
for (bad in list(list(denominator = "footprnt"), list(underreporting_multiplier = 0))) {
  e <- tryCatch({ resolveMonthlyEstimand(list(monthly = list(estimand = bad))); NULL },
                error = function(e) conditionMessage(e))
  expect(!is.null(e),
         sprintf("monthly.estimand %s is refused, not defaulted", names(bad)[1]))
}
e <- tryCatch({ resolveMonthlyEstimand(list(monthly = list(pop_elasticity = "loose"))); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e), "monthly.pop_elasticity typo is refused, not defaulted")

# ---------------------------------------------------------------------------
# [2] P1 -- resolveSizePopColumn() on a monthly estimand
# ---------------------------------------------------------------------------
cat("\n[2] P1: pop column resolution under a monthly estimand\n")

d_p1 <- data.frame(pop           = c(10, 20),
                   pop_footprint = c(10, 20),
                   pop_gid       = c(100, 200))

d_fp <- suppressMessages(resolveSizePopColumn(d_p1, es0, tag = "test"))
expect(identical(d_fp$pop, c(10, 20)), "denominator 'footprint' points pop at pop_footprint")

es_gid <- resolveMonthlyEstimand(list(monthly = list(estimand = list(denominator = "gid"))))
d_gid  <- suppressMessages(resolveSizePopColumn(d_p1, es_gid, tag = "test"))
expect(identical(d_gid$pop, c(100, 200)), "denominator 'gid' points pop at pop_gid")

# Asymmetric fallback: footprint tolerates a pre-P1 panel (back-compat), gid
# must NOT -- silently handing back a footprint population under a gid stamp is
# the exact substitution this machinery exists to prevent.
d_pre <- data.frame(pop = c(10, 20))
expect(identical(suppressMessages(resolveSizePopColumn(d_pre, es0, tag = "test"))$pop,
                 c(10, 20)),
       "pre-P1 panel + footprint falls back to `pop` (existing mon configs unchanged)")
e <- tryCatch({ suppressMessages(resolveSizePopColumn(d_pre, es_gid, tag = "test")); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("pop_gid", e, fixed = TRUE),
       "pre-P1 panel + gid is a hard error, never a silent footprint fallback")
expect(!is.null(e) && grepl("03_prepare_mon", e, fixed = TRUE),
       "the error names 03_prepare_mon.R as the prep script to re-run")

# ---------------------------------------------------------------------------
# [3] P4 -- .addLogPopCentred() nesting invariant
# ---------------------------------------------------------------------------
cat("\n[3] P4: log_pop_c centring\n")

d_p4 <- data.frame(
  admin_level = c("ADM1", "ADM1", "ADM1", "ADM2", "ADM2"),
  pop         = c(1e5, 2e5, 4e5, 1e4, 9e4)
)
d_c <- suppressMessages(.addLogPopCentred(d_p4, tag = "test/P4"))
ctr <- attr(d_c, "log_pop_centres")

expect(setequal(names(ctr), c("ADM1", "ADM2")),
       "centres are computed per admin_level, one per level present")
expect(isTRUE(all.equal(unname(ctr[["ADM1"]]), mean(log(c(1e5, 2e5, 4e5))))) &&
         isTRUE(all.equal(unname(ctr[["ADM2"]]), mean(log(c(1e4, 9e4))))),
       "each centre is the mean log(pop) of ITS level (not the global mean)")

# The nesting invariant: within each level the deviation sums to zero, so the
# term cannot move the intercept and coefficient 0 recovers `fixed` exactly.
by_lvl <- tapply(d_c$log_pop_c, d_c$admin_level, sum)
expect(all(abs(by_lvl) < 1e-10),
       "log_pop_c sums to 0 WITHIN each level (fixed and free stay nested)")
expect(isTRUE(all.equal(log(d_c$pop) - d_c$log_pop_c,
                        unname(ctr[as.character(d_c$admin_level)]))),
       "log(pop) - log_pop_c == centre[admin_level] by construction")

# The fallback path (no admin_level column) must still be usable and must use
# the "ALL" key .recover_log_pop_centres() expects.
d_all <- suppressMessages(.addLogPopCentred(data.frame(pop = c(1e5, 4e5)), tag = "test/P4"))
expect(identical(names(attr(d_all, "log_pop_centres")), "ALL"),
       "no admin_level column -> a single centre keyed 'ALL'")

# ---------------------------------------------------------------------------
# [4] P4 carriage -- the mon stamp, and .recover_log_pop_centres() reading it
# ---------------------------------------------------------------------------
cat("\n[4] P4 carriage: mon_estimand stamp -> kernel centres\n")

opt_free <- list(monthly = list(pop_elasticity = "free"))
stamp <- monthlyEstimandStamp(opt_free, extra = list(log_pop_centres = ctr))

expect(identical(stamp$track, "monthly") && identical(stamp$pop_elasticity, "free"),
       "monthlyEstimandStamp() records the monthly track and its elasticity mode")
expect(identical(stamp$pop_column, "pop_footprint"),
       "the stamp records WHICH population column the run resolved to")

# Attribute name: mon and size stamp under deliberately DISTINCT names, so a
# mon artefact cannot satisfy a size assertion.
obj  <- attachSizeEstimand(list(), stamp, track = "monthly")
expect(!is.null(getMonthlyEstimand(obj)) && is.null(getSizeEstimand(obj)),
       "a monthly stamp is readable as monthly and NOT as size")

# (1) stamp path -- the one the fit scripts now write.
fit_stamped <- structure(list(), mon_estimand = stamp)
rec <- .recover_log_pop_centres(fit_stamped, tag = "monthly")
expect(isTRUE(all.equal(rec, stats::setNames(as.numeric(ctr), names(ctr)))),
       ".recover_log_pop_centres() reads the `mon_estimand` stamp (not just size)")

fit_siz <- structure(list(), size_estimand = list(log_pop_centres = ctr))
expect(isTRUE(all.equal(.recover_log_pop_centres(fit_siz, tag = "size"),
                        stats::setNames(as.numeric(ctr), names(ctr)))),
       "the size stamp path still works (the mon lookup is additive)")

# (2) fit$data path -- what makes the fix apply to mon fits already on disk,
# which predate the stamp. brms keeps only formula variables, so admin_level
# may be absent and the level key is rebuilt from the gid.
fit_data <- structure(list(data = data.frame(
  pop       = d_p4$pop,
  log_pop_c = d_c$log_pop_c,
  gid       = c("CMR_5_1", "CMR_5_2", "CMR_5_3", "CMR_5_1_1", "CMR_5_1_2"))),
  class = "list")
expect(identical(.adminLevelFromGid(c("CMR_5_1", "CMR_5_1_1")), c("ADM1", "ADM2")),
       ".adminLevelFromGid(): gid depth is the canonical admin level")
rec2 <- .recover_log_pop_centres(fit_data, tag = "monthly")
expect(isTRUE(all.equal(rec2[c("ADM1", "ADM2")],
                        stats::setNames(as.numeric(ctr[c("ADM1", "ADM2")]),
                                        c("ADM1", "ADM2")))),
       "centres recovered exactly from fit$data via the gid-derived level key")

# (3) neither source -- hard error, never a silent 0.
e <- tryCatch({ .recover_log_pop_centres(structure(list(), class = "list"),
                                         tag = "monthly"); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("mon_estimand", e, fixed = TRUE),
       "no stamp and no usable fit$data is a hard error naming mon_estimand")

# ---------------------------------------------------------------------------
# [5] The kernel whitelist covers log_pop_c on the monthly track
# ---------------------------------------------------------------------------
cat("\n[5] kernel consumption of the log_pop_c fixed effect\n")

expect(isTRUE(.assert_lp_terms_consumed(
        c("Intercept", "mean_precip", "log_pop_c"), "monthly", "mon test")),
       "a monthly fit carrying log_pop_c passes .assert_lp_terms_consumed()")
e <- tryCatch({ .assert_lp_terms_consumed(c("Intercept", "log_pop_UNSCORED"),
                                          "monthly", "mon test"); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("log_pop_UNSCORED", e, fixed = TRUE),
       "an unscored monthly fixed effect is still a hard error (gate is live)")

# ---------------------------------------------------------------------------
# [6] Filename separation -- estimand/elasticity must reach the artefact name
# ---------------------------------------------------------------------------
cat("\n[6] _estnd- path token on monthly artefacts\n")

mk_opt <- function(...) {
  o <- list(
    baseline = "1990_2019", country = "all", outdir = "generated_data/",
    enso = list(lag_months = 0L), grouping = list(enable = FALSE),
    monthly = list(backend = "brms", family = "negbinomial",
                   seasonal_term = "harmonic2", admin_level = "ADM1-ADM2",
                   offset_pop = TRUE, phi_re = "constant"))
  o$monthly <- utils::modifyList(o$monthly, list(...))
  o
}

base_file <- makeMonthlyRunFile(mk_opt(), "mon_fit")
expect(!grepl("_estnd-", base_file, fixed = TRUE),
       "default (footprint/fixed) monthly run file carries NO _estnd- token")

free_file <- makeMonthlyRunFile(mk_opt(pop_elasticity = "free"), "mon_fit")
gid_file  <- makeMonthlyRunFile(
  mk_opt(estimand = list(denominator = "gid", rate_denominator = "gid")), "mon_fit")
both_file <- makeMonthlyRunFile(
  mk_opt(estimand = list(denominator = "gid", rate_denominator = "gid"),
         pop_elasticity = "free", canonical_admin_level = TRUE), "mon_fit")

expect(grepl("_estnd-free", free_file, fixed = TRUE),
       "pop_elasticity 'free' is encoded in the monthly run file name")
expect(grepl("_estnd-gid", gid_file, fixed = TRUE),
       "denominator 'gid' is encoded in the monthly run file name")
expect(grepl("_estnd-gid-free-canon", both_file, fixed = TRUE),
       "denominator + elasticity + canonical level compose into one token")
expect(length(unique(c(base_file, free_file, gid_file, both_file))) == 4L,
       "all four monthly estimand variants resolve to DISTINCT filenames")

# The anchor table is keyed on the same population, so it must separate too --
# a footprint-based base_rate_mon must not be reused by a gid run.
base_anchor <- makeMonthlyRunFile(mk_opt(offset_baserate = TRUE), "mon_base_rate")
gid_anchor  <- makeMonthlyRunFile(
  mk_opt(offset_baserate = TRUE,
         estimand = list(denominator = "gid", rate_denominator = "gid")),
  "mon_base_rate")
expect(!identical(base_anchor, gid_anchor),
       "the mon base-rate anchor file separates footprint from gid runs")

# The size token must be untouched by the generalisation (it is the documented
# token in every archived size artefact name).
expect(identical(.sizeEstimandPathToken(list(size = list())), "") &&
         identical(.sizeEstimandPathToken(list(size = list(pop_elasticity = "free"))),
                   "_estnd-free"),
       "the size token is unchanged by the track generalisation")
expect(identical(.estimandPathToken(list(monthly = list(
         mean_rate = list(level_factor = FALSE))), "monthly"), ""),
       "size-only mean_rate knobs are NOT encoded into the monthly token")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
