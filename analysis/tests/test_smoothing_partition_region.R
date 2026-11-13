# test_smoothing_partition_region.R — region-level (enso_region) bpfilter
# smoothing partition.
#
# Guards the "couple member countries during smoothing" feature: run the
# spatPomp historical smoothing (phase 08) as ONE enso_region partition per ENSO
# group instead of per country, so the distance-weighted neighbour-contagion
# kernel couples ADM1 units ACROSS country borders inside the particle filter
# (not only later in the forward assemble-from-members handoff). Four moving
# parts, all exercised here with NO fit / NO data / NO network:
#
#   1. The shared partition -> frame-column helpers promoted into
#      spatpomp_core.R (partition_frame_column / assert_partition_column), the
#      single source of truth consumed by both the forward simulator (07b) and
#      smoothing_apply_partition() below. Locks the three valid modes + the
#      unknown-value stop() moved out of 07b's old local partition_tu_column(),
#      and the enso_region grouping.enable hint.
#   2. smoothing_apply_partition()'s enso_region branch (smoothing_core.R): sets
#      frame$part = frame$group, leaves frame$tu INTACT (no collapse, no
#      tu_true), so the leave-one-country-out o_tu term resolves within-region
#      and force_o_tu / otu_obs both stay inert. Country mode (the regression
#      control) still collapses tu and stamps tu_true.
#   3. smoothing_tu_list() enumerates ENSO groups, and a two-country one-group
#      frame yields ONE partition whose per-partition frame holds BOTH countries'
#      ADM1 units (so build_partition_nbr_kernel will span them).
#   4. The forward direct-vs-assemble SELECTOR in load_smoothing_distribution()
#      (07b): the exact Sys.glob() probe that picks a NATIVE region cache
#      (direct load) over the per-country assemble fallback. 07b is a pipeline
#      script that runs on source(), so — exactly as test_smoothing_cache_schema.R
#      does for the same file's globs — we prove the real probe against real
#      tmpdir files rather than sourcing 07b.
#
# Run: Rscript analysis/tests/test_smoothing_partition_region.R

suppressPackageStartupMessages({
  library(qs2)
})

# spatpomp_core.R defines build_obs_for_tu(); sourcing it BEFORE smoothing_core.R
# satisfies spatpomp_core_c.R's `if (!exists("build_obs_for_tu"))` guard, so the
# transitive `source(dirname(sys.frame(1)$ofile), ...)` there is skipped (that
# ofile lookup is NULL under Rscript and would otherwise error). simulator_core.R
# is spatpomp_core.R's own prerequisite (same order as test_smoothing_cache_schema.R).
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core.R")   # partition_frame_column / assert_partition_column
source("analysis/utils/smoothing_core.R")  # smoothing_apply_partition / smoothing_tu_list

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

# A minimal two-country, one-ENSO-group frame. enso_east spans NGA (2 ADM1) and
# COD (1 ADM1); the two NGA rows share HMM TU 1, COD is TU 2.
frame2 <- function() data.frame(
  tu        = c(1L, 1L, 2L),
  country   = c("NGA", "NGA", "COD"),
  group     = c("enso_east", "enso_east", "enso_east"),
  adm1_unit = c("NGA.1", "NGA.2", "COD.1"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
cat("\n[1] shared partition -> frame-column helpers (spatpomp_core.R)\n")

expect(identical(partition_frame_column("country"), "country"),
       "partition_frame_column(country) == 'country'")
expect(identical(partition_frame_column("enso_region"), "group"),
       "partition_frame_column(enso_region) == 'group'")
expect(identical(partition_frame_column("transmission_unit"), "tu"),
       "partition_frame_column(transmission_unit) == 'tu'")
expect(inherits(try(partition_frame_column("bogus"), silent = TRUE), "try-error"),
       "partition_frame_column(<unknown>) stop()s")

expect(identical(assert_partition_column(frame2(), "enso_region"), "group"),
       "assert_partition_column returns 'group' when the column is present")

e_absent <- tryCatch(assert_partition_column(data.frame(tu = 1L), "enso_region"),
                     error = function(e) conditionMessage(e))
expect(is.character(e_absent) && grepl("grouping.enable", e_absent, fixed = TRUE),
       "assert_partition_column errors with the grouping.enable hint when `group` is absent")

e_allna <- tryCatch(
  assert_partition_column(data.frame(group = NA_character_), "enso_region"),
  error = function(e) conditionMessage(e))
expect(is.character(e_allna) && grepl("grouping.enable", e_allna, fixed = TRUE),
       "assert_partition_column errors when `group` is all-NA")

# ---------------------------------------------------------------------------
cat("\n[2] smoothing_apply_partition() enso_region branch (smoothing_core.R)\n")

fr  <- frame2()
out <- smoothing_apply_partition(list(frame = fr), list(PARTITION_BY = "enso_region"))

expect(identical(out$frame$part, as.character(fr$group)),
       "enso_region: frame$part == frame$group")
expect(identical(out$frame$tu, fr$tu),
       "enso_region: frame$tu left INTACT (not collapsed)")
expect(is.null(out$frame$tu_true),
       "enso_region: no tu_true stamped (o_tu resolves within-region; fold stays inert)")

e_noreg <- tryCatch(
  smoothing_apply_partition(list(frame = data.frame(tu = 1L, country = "NGA")),
                            list(PARTITION_BY = "enso_region")),
  error = function(e) conditionMessage(e))
expect(is.character(e_noreg) && grepl("grouping.enable", e_noreg, fixed = TRUE),
       "enso_region: errors (via assert_partition_column) when `group` is absent")

# Regression control: country mode still collapses tu and stamps tu_true.
outc <- smoothing_apply_partition(list(frame = frame2()),
                                  list(PARTITION_BY = "country"))
expect(identical(outc$frame$part, c("NGA", "NGA", "COD")),
       "country (control): frame$part == frame$country")
expect(length(unique(outc$frame$tu)) == 1L,
       "country (control): tu collapsed to a single constant")
expect(!is.null(outc$frame$tu_true) &&
         identical(as.character(outc$frame$tu_true), as.character(frame2()$tu)),
       "country (control): tu_true preserves the pre-collapse tu label")

# ---------------------------------------------------------------------------
cat("\n[3] partition enumeration + per-partition frame spans both countries\n")

expect(identical(smoothing_tu_list(out$frame), "enso_east"),
       "smoothing_tu_list enumerates ENSO groups (one 'enso_east'), not countries/TUs")

# The per-partition frame is frame[part == tu] (the exact subset
# smoothing_build_precomp uses). For the single enso_east partition it must hold
# BOTH NGA and COD ADM1 units — that cross-country frame is what makes
# build_partition_nbr_kernel couple ADM1 units across the border.
part     <- out$frame$part
frame_tu <- out$frame[!is.na(part) & part == "enso_east", , drop = FALSE]
expect(nrow(frame_tu) == 3L,
       "one enso_east partition holds all 3 ADM1 rows")
expect(setequal(unique(frame_tu$country), c("NGA", "COD")),
       "the enso_east partition frame spans BOTH member countries")
expect(setequal(frame_tu$adm1_unit, c("NGA.1", "NGA.2", "COD.1")),
       "the enso_east partition frame holds every member country's ADM1 unit")

# ---------------------------------------------------------------------------
cat("\n[4] forward direct-vs-assemble selector (07b load_smoothing_distribution probe)\n")

# The selector in load_smoothing_distribution() is
#   length(Sys.glob(file.path(dir, sprintf("hist_filter_tu=%s_*.qs", tu)))) > 0L
# => direct load of a native region cache; else assemble from per-country caches.
# Prove the real probe (not a paraphrase) against real files in a tmpdir.
probe <- function(dir, tu)
  length(Sys.glob(file.path(dir, sprintf("hist_filter_tu=%s_*.qs", tu)))) > 0L

svtok <- sm_cache_schema_token()   # current-schema token so the names are realistic

# (a) native region cache present -> direct branch fires for the region key.
dir_native <- tempfile("sm_native_"); dir.create(dir_native)
on.exit(unlink(dir_native, recursive = TRUE), add = TRUE)
invisible(file.create(file.path(dir_native,
  sprintf("hist_filter_tu=enso_east_Np=10000_%s_nrep=2_20090101_20110101.qs", svtok))))
expect(probe(dir_native, "enso_east"),
       "direct: a native hist_filter_tu=enso_east_*.qs makes the region probe hit")

# (b) only per-country caches present -> region probe misses, assemble fallback.
dir_members <- tempfile("sm_members_"); dir.create(dir_members)
on.exit(unlink(dir_members, recursive = TRUE), add = TRUE)
for (iso in c("NGA", "COD"))
  file.create(file.path(dir_members,
    sprintf("hist_filter_tu=%s_Np=10000_%s_nrep=2_20090101_20110101.qs", iso, svtok)))
expect(!probe(dir_members, "enso_east"),
       "assemble: per-country caches do NOT match the region key (region probe misses)")
expect(probe(dir_members, "NGA") && probe(dir_members, "COD"),
       "assemble: each per-country key still resolves (members are loadable to assemble)")

# The underscore in 'enso_east' is delimiter-safe: the region cache's own key
# resolves and does not get confused with a per-country one.
expect(!probe(dir_native, "enso_southcentral"),
       "a different region key does not spuriously match the enso_east cache")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
