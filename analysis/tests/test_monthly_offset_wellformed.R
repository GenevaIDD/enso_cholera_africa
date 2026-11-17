# test_monthly_offset_wellformed.R — HARD: the monthly case-count (`mon`)
#   offset contract.
#
# Usage:
#   Rscript analysis/tests/test_monthly_offset_wellformed.R
#
# Details:
#   Per build_monthly_formula_parts() and extract_monthly_draws():
#     offset(log(pop))       -- toggled by monthly.offset_pop
#     offset(log_offset)     -- UNCONDITIONAL person-time exposure
#     offset(log(mean_rate)) -- ABSENT (no bundle field, no frame column)
#   Part A pins the draw extractor's flags; Part B pins that
#   simulate_one_tu()'s mon mu block honours them. No fitted model needed.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/dlnm.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")

.passes <- 0L
.fails <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# test_monthly_offset_wellformed — the mon-track offset contract\n\n")

# ===========================================================================
# Part A -- extract_monthly_draws() hardcodes the override, not the inner
# extractor's own reading of monthly.offset_logoffset (which does not even
# exist as a monthly.* config key).
# ===========================================================================
# Part A: extract_monthly_draws() unconditionally sets offset_logoffset <- TRUE
#         on whatever .extract_nb_count_draws_brms(block = "monthly") returns
#         -- exercised by monkey-patching the (heavy, brms-fit-dependent) inner
#         extractor with a deliberately-wrong stub, so this test runs in
#         milliseconds and needs no fitted model. The inner extractor's `block`
#         argument is required.
cat(
  "## Part A: extract_monthly_draws() overrides offset flags unconditionally\n"
)
{
  orig_fn <- .extract_nb_count_draws_brms
  seen_block <- NULL

  # Deliberately WRONG bundle: if the override line in extract_monthly_draws()
  # were ever deleted, this stub's value would leak straight through.
  .extract_nb_count_draws_brms <<- function(fit, opt, max_draws = 200L, block) {
    seen_block <<- block
    list(
      S = 1L,
      fe = matrix(
        0,
        nrow = 1L,
        ncol = 1L,
        dimnames = list(NULL, "(Intercept)")
      ),
      re_intercept = list(
        region = NULL,
        country = NULL,
        adm1 = NULL,
        adm2 = NULL
      ),
      re_slope = list(),
      basehaz = list(kind = "none"),
      phi = 1,
      phi_re = "constant",
      active_terms = "(Intercept)",
      offset_pop = TRUE,
      offset_logoffset = FALSE # WRONG for mon -- extract_monthly_draws() must clobber this
    )
  }

  out <- tryCatch(
    extract_monthly_draws(fit = NULL, opt = list(monthly = list())),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )

  .extract_nb_count_draws_brms <<- orig_fn # restore before any assertion can short-circuit

  expect_true(
    !is.null(out),
    "extract_monthly_draws() runs against a stubbed inner extractor"
  )
  expect_true(
    identical(seen_block, "monthly"),
    "extract_monthly_draws() calls the NB count extractor with block = \"monthly\""
  )
  if (!is.null(out)) {
    expect_true(
      !("offset_meanrate" %in% names(out)),
      "the bundle carries no offset_meanrate field (size track's mean_rate offset retired)"
    )
    expect_true(
      identical(out$offset_logoffset, TRUE),
      "extract_monthly_draws() forces offset_logoffset = TRUE regardless of the inner extractor"
    )
    expect_true(
      identical(out$offset_pop, TRUE),
      "offset_pop passes through from the inner extractor unmodified (monthly.offset_pop is a real toggle)"
    )
  }
}

# --- the NB count extractor has no default block (formerly "size") -----------
cat("\n## Part A1: .extract_nb_count_draws_brms() requires `block`\n")
{
  err <- tryCatch(
    .extract_nb_count_draws_brms(fit = NULL, opt = list()),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.character(err) && grepl("block", err),
    ".extract_nb_count_draws_brms() without `block` raises a clear error"
  )
  expect_true(
    !exists(".extract_size_draws_brms") && !exists("extract_size_draws"),
    "the size extractors are deleted, not kept side-by-side"
  )
}

# --- prev_cases: true still stops loudly (adjacent one-liner, cheap to pin here) --
cat(
  "\n## Part A2: monthly.prev_cases = TRUE stops rather than silently approximating\n"
)
{
  err <- tryCatch(
    extract_monthly_draws(
      fit = NULL,
      opt = list(monthly = list(prev_cases = TRUE))
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    is.character(err) && grepl("prev_cases", err),
    "monthly.prev_cases = TRUE raises a clear error before touching `fit`"
  )
}

# ===========================================================================
# Part B -- kernel wiring: simulate_one_tu()'s mon mu block never reads
# mean_rate, and does read pop / log_offset when their flags are on.
# ===========================================================================
# Part B: with a bundle carrying the real-world flags (offset_pop=TRUE,
#         offset_logoffset=TRUE), simulate_one_tu()'s mon mu block actually
#         behaves accordingly: a stray `mean_rate` frame column and a stray
#         `offset_meanrate = TRUE` bundle field are both ignored (mu_mon and the
#         drawn cases_month are byte-identical), while pop and log_offset each
#         do move mu when their flag is on -- ruling out the trivial "the block
#         is reading nothing at all" false pass.
cat(
  "\n## Part B: simulate_one_tu() mon mu block honours the bundle's offset flags\n"
)

# `mean_rate` is a STRAY column here (the location frame no longer carries
# one); NULL omits it entirely.
fixture_frame <- function(pop = 1e5, mean_rate = NULL) {
  fr <- tibble::tibble(
    gid = "A_1_1",
    country = "AAA",
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = 1L,
    pop = pop,
    adm1_unit = "A_1"
  )
  if (!is.null(mean_rate)) {
    fr$mean_rate <- mean_rate
  }
  fr
}
make_panel <- function(frame_tu, months) {
  tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = lubridate::yday(t_start),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L,
      time_to_last_outbreak_cat = NA
    )
}

months <- as.Date("2024-01-01") # single month -> forced onset is active in month 1
horizon <- list(start = "2024-01-01", end = "2024-01-31")

# Forced onset (occ intercept huge), cessation hazard ~0 (dur intercept very
# negative) -> active = 1 in the one simulated month for every run below, so
# mu_mon/cases_month are always populated and directly comparable.
#
# `run_months` / `run_horizon` default to the single-month fixture above, so
# B1-B3 are unchanged; B4 overrides them to get interior months as well as the
# onset month (see the .ONSET_EXPOSURE_FRAC correction it now checks).
run_mon <- function(
  pop,
  mean_rate = NULL,
  mon_overrides = list(),
  seed = 42L,
  run_months = months,
  run_horizon = horizon
) {
  frame_tu <- fixture_frame(pop = pop, mean_rate = mean_rate)
  panel <- make_panel(frame_tu, run_months)

  draws_occ <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = 30)
  )
  draws_dur <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = -30),
    active_terms = c("(Intercept)")
  )
  draws_mon <- make_draws_bundle_synthetic(
    frame_tu,
    S = 1,
    fe_overrides = c(`(Intercept)` = log(20)),
    active_terms = c("(Intercept)")
  )
  draws_mon$offset_pop <- TRUE
  draws_mon$offset_logoffset <- TRUE
  for (nm in names(mon_overrides)) {
    draws_mon[[nm]] <- mon_overrides[[nm]]
  }

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    replicate_ids = 1L,
    seed = seed,
    output_part = tempfile("mon_offset_test_"),
    horizon = run_horizon,
    dry_run = TRUE
  )
  res$monthly_summary
}

# --- B1: a stray mean_rate frame column is never read ------------------------
ms_none <- run_mon(pop = 1e5)
ms_lo <- run_mon(pop = 1e5, mean_rate = 1e-3)
ms_hi <- run_mon(pop = 1e5, mean_rate = 50) # wildly different mean_rate

expect_true(
  all(!is.na(ms_lo$mu_mon)) && all(!is.na(ms_hi$mu_mon)),
  "mu_mon populated (active=1 forced) in both mean_rate variants"
)
expect_equal(
  ms_lo$mu_mon,
  ms_hi$mu_mon,
  "mu_mon is byte-identical across mean_rate=1e-3 vs 50: mean_rate never read"
)
expect_equal(
  ms_none$mu_mon,
  ms_lo$mu_mon,
  "mu_mon identical with and without a mean_rate column on the frame"
)
expect_equal(
  ms_lo$cases_month,
  ms_hi$cases_month,
  "cases_month draw is identical too (same RNG stream, same mu -> mean_rate genuinely inert, not just averaged out)"
)

# --- B2: a stray offset_meanrate = TRUE bundle field is ignored too -----------
# (the mon block's dead offset_meanrate branch was deleted with the size track)
ms_lo_on <- run_mon(
  pop = 1e5,
  mean_rate = 1e-3,
  mon_overrides = list(offset_meanrate = TRUE)
)
ms_hi_on <- run_mon(
  pop = 1e5,
  mean_rate = 50,
  mon_overrides = list(offset_meanrate = TRUE)
)

expect_equal(
  ms_lo_on$mu_mon,
  ms_hi_on$mu_mon,
  "a stray offset_meanrate = TRUE on the bundle does not make the kernel read mean_rate"
)
expect_equal(
  ms_lo_on$mu_mon,
  ms_lo$mu_mon,
  "offset_meanrate = TRUE vs absent: mu_mon unchanged"
)

# --- B3: offset_pop (real-world default TRUE) moves mu proportionally to pop ---
# (rules out "the block reads nothing")
ms_pop_lo <- run_mon(pop = 1e5)
ms_pop_hi <- run_mon(pop = 2e5)

expect_true(
  isTRUE(all.equal(
    ms_pop_hi$mu_mon / ms_pop_lo$mu_mon,
    rep(2, nrow(ms_pop_lo))
  )),
  "offset_pop=TRUE: doubling pop doubles mu_mon (log-link offset(log(pop)))"
)

# --- B4: offset_logoffset adds log(days_in_month), MINUS the onset-month
# ---     half-month exposure correction (.ONSET_EXPOSURE_FRAC).
#
# A month in the replay is not always a full month of person-time: the episode
# starts partway through its onset month. The fit's log_offset column is the
# real overlap (dur_data.R:257), so the per-day rate is estimated having already
# divided that out; only the replay, which multiplies by days_in_month(), is
# wrong. .onset_log_offset_adj() (simulator_core.R) adds log(0.5) on mso == 1.
#
# Three months here, not one: month 1 is the onset month (mso == 1, forced by
# the huge occ intercept), months 2-3 are interior (mso 2, 3 -- cessation hazard
# is ~0 and onset cannot re-fire on an already-active unit, since the at-risk
# set is `active_before == 0`). So one assertion pins all three facts at once:
# the offset is applied, it is halved exactly on the onset month, and it is
# exactly days_in_month on interior months.
months3 <- seq(as.Date("2024-01-01"), by = "month", length.out = 3L)
horizon3 <- list(start = "2024-01-01", end = "2024-03-31")

ms_logoff_on <- run_mon(
  pop = 1e5,
  mon_overrides = list(offset_logoffset = TRUE),
  run_months = months3,
  run_horizon = horizon3
)
ms_logoff_off <- run_mon(
  pop = 1e5,
  mon_overrides = list(offset_logoffset = FALSE),
  run_months = months3,
  run_horizon = horizon3
)

ord <- order(ms_logoff_on$month)
days3 <- as.integer(lubridate::days_in_month(months3)) # 31, 29 (leap), 31
expected <- days3 * c(.ONSET_EXPOSURE_FRAC, 1, 1) # 15.5, 29, 31

expect_true(
  nrow(ms_logoff_on) == 3L && nrow(ms_logoff_off) == 3L,
  "B4 fixture: all 3 months are active (forced onset, ~zero cessation hazard)"
)
expect_true(
  isTRUE(all.equal(
    ms_logoff_on$mu_mon[ord] / ms_logoff_off$mu_mon[ord],
    expected
  )),
  sprintf(
    paste0(
      "offset_logoffset=TRUE multiplies mu_mon by days_in_month, ",
      "halved on the onset month: %s vs. FALSE"
    ),
    paste(expected, collapse = ", ")
  )
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
