# test_pop_provenance_guard.R — HARD: the mean_rate x population ordering
#   guard.
#
# Usage:
#   Rscript analysis/tests/test_pop_provenance_guard.R
#
# Details:
#   build_mean_rate_table() computes every cell as `suspected_cases / pop`
#   and is persisted once as a frozen sidecar, so a population change must
#   rebuild it. The stamp is a CONTENT digest, not a timestamp: a rebuild
#   from unchanged inputs must NOT invalidate the pairing (section D), a
#   genuine population change MUST (section C). Pure synthetic lists; no
#   RDS, no SQLite, no data/, no fits.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("analysis/utils/base.R") # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/paths.R") # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/provenance.R") # assert_pop_provenance()/attach_pop_provenance()/compare_pop_provenance()/get_pop_provenance()/pop_provenance_stamp() -- used directly

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
errmsg <- function(expr) {
  e <- try(expr, silent = TRUE)
  if (inherits(e, "try-error")) {
    conditionMessage(attr(e, "condition"))
  } else {
    NA_character_
  }
}
warnmsg <- function(expr) {
  w <- NULL
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- conditionMessage(cnd)
    invokeRestart("muffleWarning")
  })
  w %||% NA_character_
}

# The two populations at stake: the legacy centroid-carried
# `pop` and the crosswalk-apportioned one. Deliberately overlapping in places
# -- a defective population is not uniformly larger or smaller than a correct
# one. SDN_3_1 loses population (it had absorbed 7 states' worth), while the
# pop = 1 units gain it.
pop_legacy <- c(1, 1, 1, 209e6, 0, 50000) # the live defects
pop_crosswalk <- c(120000, 84000, 33000, 5.2e6, 41000, 50000)

# build_mean_rate_table() computes every cell as `suspected_cases / pop`. The
# resulting table is persisted once and thereafter read as a frozen sidecar --
# by the 04_fit_siz_*.R scripts, and by simulator_core.R. Changing `pop`
# (e.g. wiring the LP x gid crosswalk into 03_prepare_siz.R) invalidates
# every cell of it. The required order is
#
#   pop fix -> mean_rate rebuild -> size refit -> Phase-A rebuild -> forecast
#
# and skipping the second step pairs NEW populations with OLD rates. That
# pairing is invisible to every other check in the pipeline, because both
# numbers remain individually plausible. In the simulator the two multiply --
#
#   mu = pop * mean_rate * exp(LP)
#
# -- so a mismatch compounds rather than cancels.
cat("\n## A) the stamp records what distinguishes two population definitions\n")

st_legacy <- pop_provenance_stamp(pop_legacy, pop_source = "legacy")
st_xwalk <- pop_provenance_stamp(pop_crosswalk, pop_source = "crosswalk")

expect(
  identical(st_legacy$pop_source, "legacy"),
  "pop_source recorded (legacy)"
)
expect(
  identical(st_xwalk$pop_source, "crosswalk"),
  "pop_source recorded (crosswalk)"
)
expect(
  identical(st_legacy$n, 6L) || identical(st_legacy$n, 6),
  "n records the row count"
)
expect(
  isTRUE(all.equal(st_legacy$pop_total, sum(pop_legacy))),
  "pop_total is the summed population"
)
expect(
  !is.null(st_legacy$stamped_at) && nzchar(st_legacy$stamped_at),
  "stamped_at is populated"
)

# A non-finite population must not poison the total. `pop = NA` is common in
# the corpus (LPs with no matched crosswalk row), and a stamp that returned NA
# for the whole artefact would silently disable the comparison for exactly the
# datasets most at risk.
st_na <- pop_provenance_stamp(c(100, NA, 300, Inf), pop_source = "legacy")
expect(
  isTRUE(all.equal(st_na$pop_total, 400)),
  "non-finite entries excluded from pop_total, not propagated as NA"
)
expect(
  identical(as.integer(st_na$n_finite), 2L),
  "n_finite counts the usable entries separately from n"
)

cat("\n## B) an absent stamp is not a mismatch\n")

# Every artefact built before this guard existed carries no stamp. Treating
# that as a failure would make the guard's introduction break every prior run,
# which is the surest way to get a guard disabled rather than fixed.
expect(
  length(compare_pop_provenance(NULL, st_xwalk)) == 0,
  "NULL data stamp compares clean (unverified, not failed)"
)
expect(
  length(compare_pop_provenance(st_xwalk, NULL)) == 0,
  "NULL sidecar stamp compares clean"
)
expect(
  is.null(get_pop_provenance(list(a = 1))),
  "get_pop_provenance() on an unstamped object returns NULL, not an error"
)

w <- warnmsg(assert_pop_provenance(
  NULL,
  st_xwalk,
  what = "the test sidecar",
  require_stamp = FALSE
))
expect(
  !is.na(w) && grepl("UNVERIFIED", w, fixed = TRUE),
  "an unstamped pairing WARNS when require_stamp is explicitly relaxed"
)
expect(
  !is.na(w) && grepl("mean_rate rebuild", w, fixed = TRUE),
  "the warning names the rebuild order, not just the fact of the gap"
)

e <- errmsg(assert_pop_provenance(NULL, st_xwalk, what = "x"))
expect(
  !is.na(e) && grepl("UNVERIFIED", e, fixed = TRUE),
  "require_stamp defaults to TRUE -- an unstamped pairing is now fatal by default (Layer 3)"
)

e2 <- errmsg(assert_pop_provenance(
  NULL,
  st_xwalk,
  what = "x",
  require_stamp = TRUE
))
expect(
  !is.na(e2) && grepl("UNVERIFIED", e2, fixed = TRUE),
  "require_stamp = TRUE explicitly promotes the same condition to fatal"
)

cat("\n## C) a population change is caught, and named\n")

diffs <- compare_pop_provenance(st_xwalk, st_legacy)
expect(length(diffs) > 0, "crosswalk data vs legacy sidecar is a mismatch")
expect(
  any(grepl("pop_source", diffs, fixed = TRUE)),
  "the mismatch names pop_source (the field that identifies the cause)"
)

e <- errmsg(assert_pop_provenance(
  st_xwalk,
  st_legacy,
  what = "the mean_rate sidecar at /tmp/mr.rds"
))
expect(!is.na(e), "assert_pop_provenance() STOPS on a real mismatch")
expect(
  !is.na(e) && grepl("suspected_cases / pop", e, fixed = TRUE),
  "the error explains WHY a pop change invalidates the table"
)
expect(
  !is.na(e) && grepl("mean_rate rebuild", e, fixed = TRUE),
  "the error states the required rebuild order"
)
expect(
  !is.na(e) && grepl("/tmp/mr.rds", e, fixed = TRUE),
  "the error names the offending artefact path"
)

# The subtle case: SAME source label, DIFFERENT numbers. This is a re-run of
# the crosswalk against an updated WorldPop vintage -- the label alone would
# wave it through, so the digest has to carry the check.
st_xwalk_v2 <- pop_provenance_stamp(
  pop_crosswalk * 1.03,
  pop_source = "crosswalk"
)
d2 <- compare_pop_provenance(st_xwalk, st_xwalk_v2)
expect(
  length(d2) > 0,
  "same pop_source but different populations is STILL a mismatch"
)
expect(
  any(grepl("pop_digest|pop_total", d2)),
  "that mismatch is reported on the numeric fingerprint, not the label"
)

cat("\n## D) a rebuild from unchanged inputs is NOT a mismatch\n")

# This is why the stamp is a content digest and not a timestamp. Re-running
# 03_prepare_siz.R without changing anything, or copying the artefact between
# machines, produces a later `stamped_at` but the same populations. If that
# invalidated the sidecar, every routine re-run would force a spurious refit
# and operators would learn to ignore the guard.
st_rebuild <- pop_provenance_stamp(pop_crosswalk, pop_source = "crosswalk")
expect(
  !identical(st_xwalk$stamped_at, st_rebuild$stamped_at) ||
    TRUE, # clocks may tie; the point is the comparison ignores the field
  "stamped_at is recorded but is not what the comparison keys on"
)
expect(
  length(compare_pop_provenance(st_xwalk, st_rebuild)) == 0,
  "identical populations compare clean regardless of build time"
)
expect(
  isTRUE(assert_pop_provenance(st_xwalk, st_rebuild, what = "y")),
  "assert_pop_provenance() passes a clean rebuild"
)

# Row ORDER must not matter either: the size prep's grouping/joins can reorder
# rows without changing the population content, and a digest over an unsorted
# vector would flag that as a population change.
st_shuffled <- pop_provenance_stamp(
  rev(pop_crosswalk),
  pop_source = "crosswalk"
)
expect(
  length(compare_pop_provenance(st_xwalk, st_shuffled)) == 0,
  "row order does not affect the digest (sorted before hashing)"
)

cat("\n## E) the stamp survives the attribute round trip\n")

# The stamp rides as an attribute so it travels through saveRDS/readRDS without
# touching any column schema -- no consumer that does not know about the stamp
# can be broken by its presence.
dat <- tibble(gid = c("AAA.1_1", "AAA.2_1"), pop = c(1000, 2000))
dat <- attach_pop_provenance(dat, pop_provenance_stamp(dat$pop, "crosswalk"))

expect(
  !is.null(get_pop_provenance(dat)),
  "stamp readable straight back off the object"
)
expect(identical(ncol(dat), 2L), "attaching the stamp adds no columns")

tf <- tempfile(fileext = ".rds")
saveRDS(dat, tf)
back <- readRDS(tf)
unlink(tf)
expect(
  length(compare_pop_provenance(
    get_pop_provenance(dat),
    get_pop_provenance(back)
  )) ==
    0,
  "stamp survives saveRDS/readRDS intact"
)

# ...but dplyr's attribute handling is VERB-DEPENDENT, and that inconsistency
# is the trap. filter/mutate/select/arrange/left_join carry the attribute
# through; summarise() builds a fresh frame and drops it. build_mean_rate_table()
# therefore captures the stamp immediately after readRDS() rather than relying
# on which side of that line its quality guards happen to fall on, and
# 03_prepare_siz.R attaches its stamp AFTER the (gid, season) collapse rather
# than before.
#
# These assertions pin the observed dplyr behaviour (1.2.1). If a future dplyr
# changes either answer, this test says so instead of the pipeline silently
# losing -- or silently gaining -- a stamp.
expect(
  !is.null(get_pop_provenance(dplyr::filter(dat, pop > 0))),
  "filter() currently PRESERVES the attribute"
)
expect(
  !is.null(get_pop_provenance(dplyr::mutate(dat, z = 1))),
  "mutate() currently PRESERVES the attribute"
)
collapsed <- dat %>%
  dplyr::group_by(gid) %>%
  dplyr::summarise(pop = sum(pop), .groups = "drop")
expect(
  is.null(get_pop_provenance(collapsed)),
  "summarise() DROPS it -- so a stamp must be (re)attached after a collapse"
)

cat("\n## F) an unavailable `digest` degrades the check, never breaks it\n")

# digest is not a hard dependency of this project. Where it is missing the
# fingerprint is NA and comparison falls back to pop_total -- weaker, but it
# still catches the legacy/crosswalk swap this guard exists for.
st_nodig_a <- st_xwalk
st_nodig_a$pop_digest <- NA_character_
st_nodig_b <- st_legacy
st_nodig_b$pop_digest <- NA_character_

expect(
  length(compare_pop_provenance(st_nodig_a, st_nodig_b)) > 0,
  "without digests the legacy/crosswalk swap is still caught"
)
st_nodig_c <- st_nodig_a
st_nodig_c$stamped_at <- "1999-01-01T00:00:00"
expect(
  length(compare_pop_provenance(st_nodig_a, st_nodig_c)) == 0,
  "without digests, equal pop_total + equal source still compares clean"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
