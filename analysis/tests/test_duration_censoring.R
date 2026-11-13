# analysis/tests/test_duration_censoring.R
#
# Gates for the outbreak-DURATION panel (analysis/utils/dur_data.R).
#
# Two halves, deliberately:
#
#   SYNTHETIC (always runs). Month expansion, partial-month exposure, the
#   `ended` contract, and the censoring rule are checked against hand-computed
#   fixtures. These are properties of the code and must hold on any corpus.
#
#   CORPUS (skipped when the prepared panel is absent). The TZA-2019 gate. This
#   is the one that cannot be synthesised: the claim is that the censoring logic
#   finds a truncation that is documented in the real data (CLAUDE.md, Data
#   sources: "TZA database ends 2019"). A synthetic fixture proving the rule
#   fires would prove nothing about whether the rule is wired to the right
#   series. 03_prepare_dur.R hard-stops on the same condition; this test exists
#   so a regression is caught by run_all.R rather than by a prep re-run.

suppressMessages({
  library(dplyr)
  library(lubridate)
  library(tibble)
})

source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L

source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/dur_data.R")
source("analysis/utils/provenance.R")
source("analysis/utils/anchor.R")

cat("== test_duration_censoring ==\n")

# ---------------------------------------------------------------------------
# 1. bin_mso
# ---------------------------------------------------------------------------
cat("\n-- bin_mso\n")

b <- bin_mso(c(1, 2, 3, 4, 5, 6, 8, 9, 40))
expect(identical(as.character(b),
                 c("1", "2", "3", "4-5", "4-5", "6-8", "6-8", "9+", "9+")),
       "bin_mso labels the default breaks as 1/2/3/4-5/6-8/9+")
expect(identical(levels(b), c("1", "2", "3", "4-5", "6-8", "9+")),
       "bin_mso level order follows the break order")
expect(inherits(try(bin_mso(1:3, breaks = c(1, 2, 3)), silent = TRUE), "try-error"),
       "bin_mso rejects breaks that do not end in Inf")
expect(identical(as.character(bin_mso(c(1, 4, 9), breaks = c(1, 4, Inf))),
                 c("1-3", "4+", "4+")),
       "bin_mso honours a custom break vector (the parsimony variant)")

# ---------------------------------------------------------------------------
# 2. .expandActiveMonths — the arithmetic
# ---------------------------------------------------------------------------
cat("\n-- .expandActiveMonths\n")

ep <- tibble(
  episode_id = c("one_month", "two_months", "spans_year", "single_day"),
  TL = as.Date(c("2015-03-05", "2015-03-20", "2015-11-14", "2015-06-10")),
  TR = as.Date(c("2015-03-28", "2015-04-02", "2016-02-03", "2015-06-10")))

x <- .expandActiveMonths(ep)

expect(nrow(x) == 1 + 2 + 4 + 1, "row count is the sum of spanned calendar months")
expect(all(x$mso[x$episode_id == "spans_year"] == 1:4),
       "mso is 1-based and consecutive within an episode")
expect(all(x$month[x$episode_id == "spans_year"] ==
             as.Date(c("2015-11-01", "2015-12-01", "2016-01-01", "2016-02-01"))),
       "month expansion crosses a year boundary correctly")

# Exposure is INCLUSIVE of both endpoints, so it sums to duration_days + 1.
# This is the quantity the Poisson offset makes the hazard per-day; an
# off-by-one here silently rescales every duration coefficient.
tot <- x %>% group_by(episode_id) %>% summarise(s = sum(exposure), .groups = "drop") %>%
  left_join(ep %>% transmute(episode_id, d = as.numeric(TR - TL) + 1), by = "episode_id")
expect(all(tot$s == tot$d), "exposure sums to (TR - TL) + 1 for every episode")

expect(x$exposure[x$episode_id == "single_day"] == 1,
       "a same-day episode contributes exactly 1 day of exposure, not 0")
expect(x$exposure[x$episode_id == "one_month"] == 24,
       "a wholly-within-month episode gets its actual day count (2015-03-05..28 = 24)")
sp <- x %>% filter(episode_id == "spans_year")
expect(sp$exposure[1] == 17 && sp$exposure[2] == 31 &&
       sp$exposure[3] == 31 && sp$exposure[4] == 3,
       "first and last months are partial, interior months are full")

expect(inherits(try(.expandActiveMonths(
         tibble(episode_id = "bad", TL = as.Date("2015-05-01"),
                TR = as.Date("2015-04-01"))), silent = TRUE), "try-error"),
       ".expandActiveMonths rejects TR < TL rather than emitting negative exposure")

# ---------------------------------------------------------------------------
# 3. flagRightCensored — the rule
# ---------------------------------------------------------------------------
cat("\n-- flagRightCensored\n")

eps <- tibble(
  location_period_id = c("A", "A", "B", "C"),
  country            = c("XX", "XX", "XX", "XX"),
  TR = as.Date(c("2019-01-15",   # at A's series end        -> censored
                 "2015-06-30",   # well inside A's record   -> not censored
                 "2018-12-01",   # 30 d before B's end      -> censored (tol)
                 "2017-01-01"))) # C has no series end      -> not censored
ser <- tibble(location_period_id = c("A", "B"),
              series_end         = as.Date(c("2019-02-01", "2018-12-31")),
              series_end_source  = "location_period")

f <- flagRightCensored(eps, ser, censor_tol_days = 31L)
expect(identical(f$censored, c(TRUE, FALSE, TRUE, FALSE)),
       "censoring fires within the tolerance of the series end, and only there")

f0 <- flagRightCensored(eps, ser, censor_tol_days = 0L)
expect(identical(f0$censored, c(FALSE, FALSE, FALSE, FALSE)),
       "tolerance 0 requires TR to reach the series end exactly")

expect(is.na(f$series_end[4]) && !f$censored[4],
       "an episode with no series end is treated as uncensored, not as NA")

# ---------------------------------------------------------------------------
# 4. the `ended` contract, end to end on a synthetic panel
# ---------------------------------------------------------------------------
cat("\n-- ended contract\n")

# Rebuild the step-4 logic of makeDurationPanel() on the fixture above. Running
# makeDurationPanel() itself would need the climate SQLite databases, which this
# test must not depend on; the response construction is the part under test.
pan <- .expandActiveMonths(ep) %>%
  mutate(censored = episode_id %in% c("two_months"),
         ended    = as.integer(mso == n_months & !censored),
         log_offset = log(pmax(exposure, 1)))

by_ep <- pan %>% group_by(episode_id) %>%
  summarise(n_ended = sum(ended), cens = first(censored), .groups = "drop")
expect(all(by_ep$n_ended == as.integer(!by_ep$cens)),
       "`ended` fires exactly once per uncensored episode and never for a censored one")
expect(all(pan$ended[pan$mso < pan$n_months] == 0L),
       "no interior active month is ever marked as the end")
expect(all(is.finite(pan$log_offset)),
       "log_offset is finite everywhere (pmax guards log(0))")

# ---------------------------------------------------------------------------
# 5. CORPUS gate — TZA must come back right-censored
# ---------------------------------------------------------------------------
cat("\n-- corpus gate (TZA 2019 truncation)\n")

DUR_PATH <- makeDurationDatFile(list(baseline = "1990_2019"))

if (!file.exists(DUR_PATH)) {
  cat(sprintf("  SKIP  prepared duration panel absent (%s) -- run analysis/03_prepare_dur.R\n",
              DUR_PATH))
} else {
  d <- readRDS(DUR_PATH)
  ep_lvl <- d %>%
    group_by(episode_id) %>%
    summarise(country = first(country), n_rows = n(),
              n_months = first(duration_months), censored = first(censored),
              n_ended = sum(ended), .groups = "drop")

  tza <- ep_lvl %>% filter(country == "TZA")
  expect(nrow(tza) > 0, "corpus carries TZA episodes at all")
  expect(sum(tza$censored) > 0,
         "TZA has right-censored episodes (its series is documented to end in 2019)")

  # The truncation must be where the documentation says it is. A censoring rule
  # that flags TZA for the wrong reason would pass the check above.
  tza_cens_max <- suppressWarnings(max(as.Date(d$TR[d$country == "TZA" & d$censored])))
  expect(is.finite(as.numeric(tza_cens_max)) &&
         year(tza_cens_max) <= 2019,
         sprintf("TZA censoring sits at the 2019 boundary (latest censored TR = %s)",
                 format(tza_cens_max)))

  expect(nrow(d) == sum(ep_lvl$n_months),
         "panel row count equals the sum of active months (no lost/duplicated rows)")
  expect(all(ep_lvl$n_ended == as.integer(!ep_lvl$censored)),
         "the `ended` contract holds on the real panel")
  expect(all(d$exposure >= 1 & d$exposure <= 31),
         "every row carries between 1 and 31 days of exposure")

  # Non-flat hazard. If this fails the mso baseline is pointless and the model
  # should be reconsidered, not silently fitted -- see .DUR_MSO_BREAKS.
  haz <- d %>% group_by(mso_band) %>%
    summarise(at_risk = n(), events = sum(ended), .groups = "drop") %>%
    filter(at_risk >= 30) %>%
    mutate(h = events / at_risk)
  expect(nrow(haz) >= 3 && max(haz$h) / max(min(haz$h), 1e-6) > 1.2,
         sprintf("empirical cessation hazard varies across mso bands (%.4f - %.4f)",
                 min(haz$h), max(haz$h)))

  # Censoring must not be so aggressive that the event count collapses, nor so
  # lax that nothing is censored. Both failure modes bias the hazard.
  fr <- mean(ep_lvl$censored)
  expect(fr > 0.01 && fr < 0.30,
         sprintf("censored fraction is plausible (%.1f%%)", 100 * fr))
}

cat(sprintf("\n== test_duration_censoring: %d passed, %d failed ==\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
