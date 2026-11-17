# test_results_explosive_guard.R — HARD: the explosive-replicate guard as a
#   results decision (utils/results_derive.R).
#
# Usage:
#   Rscript analysis/tests/test_results_explosive_guard.R
#
# Details:
#   Quick, synthetic inputs only. Checks explosive_pair_summary() (table S29),
#   the per-window summary cells, that the cap and modes come from the context
#   (no silent constant), Figure 4a's group-level drops, stamp back-compat and
#   the per-gid monthly cap (results.thresholds.explosive_gid_month_cap).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(a, b, tolerance = tol))
err <- function(e) {
  inherits(try(suppressMessages(e), silent = TRUE), "try-error")
}
D <- as.Date

# [1] explosive_pair_summary() (table S29) by hand: per-arm counts, F-only /
#     CF-only / both drops, kept pairs, the exact sign test, AN over the
#     dropped pairs, AF none - AF drop; Africa compounds (a replicate
#     explosive in one group drops the Africa pair).
cat("# [1] explosive_pair_summary()\n")
# Two groups, 6 replicates. Group A: F = 10 except r1, r2 (300: F only), r3
# (300 in both); CF = 5 except r3 (300). Group B: F = 4, CF = 2, r6 CF 300.
cum <- bind_rows(
  tibble(
    enso_group = "A",
    run = "forecast",
    analog_set = NA_character_,
    replicate_id = 1:6,
    cum_cases = c(300, 300, 300, 10, 10, 10)
  ),
  tibble(
    enso_group = "A",
    run = "counterfactual",
    analog_set = "neutral",
    replicate_id = 1:6,
    cum_cases = c(5, 5, 300, 5, 5, 5)
  ),
  tibble(
    enso_group = "B",
    run = "forecast",
    analog_set = NA_character_,
    replicate_id = 1:6,
    cum_cases = 4
  ),
  tibble(
    enso_group = "B",
    run = "counterfactual",
    analog_set = "neutral",
    replicate_id = 1:6,
    cum_cases = c(2, 2, 2, 2, 2, 300)
  ),
  tibble(
    enso_group = "A",
    run = "counterfactual",
    analog_set = "other",
    replicate_id = 1:6,
    cum_cases = 1
  )
)
s <- explosive_pair_summary(cum, 100, "neutral")
a <- filter(s, scope == "A")
b <- filter(s, scope == "B")
af <- filter(s, scope == "Africa")
expect(
  a$n_rep_f == 6L &&
    a$n_expl_f == 3L &&
    a$n_expl_cf == 1L &&
    a$n_drop_f_only == 2L &&
    a$n_drop_cf_only == 0L &&
    a$n_drop_both == 1L &&
    a$n_pairs_kept == 3L &&
    near(a$share_dropped, 0.5),
  "group A: counts per arm, F-only / CF-only / both, kept pairs, share"
)
expect(
  near(a$sign_p, stats::binom.test(2, 2)$p.value),
  "group A: exact sign test of F-only vs CF-only drops"
)
expect(
  near(a$an_dropped, (300 - 5) * 2 + 0),
  "group A: AN over the dropped pairs"
)
expect(
  near(
    a$af_none,
    1 - mean(c(5, 5, 300, 5, 5, 5)) / mean(c(300, 300, 300, 10, 10, 10))
  ) &&
    near(a$af_drop, 1 - 5 / 10) &&
    near(a$af_none_minus_drop, a$af_none - a$af_drop),
  "group A: AF with every pair, with the kept pairs, difference"
)
expect(
  b$n_drop_cf_only == 1L &&
    b$n_drop_f_only == 0L &&
    near(b$an_dropped, 4 - 300),
  "group B: one CF-only drop"
)
expect(
  af$n_pairs == 6L &&
    af$n_drop_f_only == 2L &&
    af$n_drop_both == 1L &&
    af$n_drop_cf_only == 1L &&
    af$n_pairs_kept == 2L,
  "Africa compounds: a pair explosive in any group or arm is dropped (4 of 6)"
)
expect(
  near(af$af_drop, 1 - (5 + 2) / (10 + 4)),
  "Africa AF over the kept pairs sums the groups"
)
expect(
  setequal(unique(s$analog_set), "neutral") && nrow(s) == 3L,
  "one row per scope of the requested set"
)
s0 <- explosive_pair_summary(cum, Inf, "neutral")
expect(
  all(s0$n_pairs_kept == s0$n_pairs) &&
    all(is.na(s0$sign_p)) &&
    all(abs(s0$af_none_minus_drop) < 1e-12),
  "cap Inf: nothing dropped, no sign test, AF unchanged"
)

# [2] explosive_summary_cells() per window (cells cut by .cell_frame()); no
#     "drop" mode -> no summary; attrib_information_explosive() for the
#     stitched vs counterfactual / stitched vs May comparisons.
cat("# [2] explosive_summary_cells(), attrib_information_explosive()\n")
months <- seq(D("2026-05-01"), by = "month", length.out = 3)
cm <- tidyr::expand_grid(
  tibble(run = c("forecast", "counterfactual"), analog_set = c(NA, "neutral")),
  country = c("AAA", "BBB"),
  replicate_id = 1:5,
  month = months
) %>%
  mutate(
    enso_group = ifelse(country == "AAA", "A", "B"),
    group_horizon_end = D("2026-07-01"),
    n_cases = ifelse(run == "forecast", 20, 10),
    n_outbreaks = 1L
  )
cm$n_cases[
  cm$run == "forecast" &
    cm$replicate_id == 2 &
    cm$country == "AAA" &
    cm$month == months[3]
] <- 500
cells <- attrib_cell_table(D("2026-05-01"), D("2026-07-01"), cap = 200)
es <- explosive_summary_cells(cm, cells, "neutral", "forecast")
expect(
  setequal(es$window, c("attribution", "post_fork")) &&
    all(es$factual == "forecast"),
  "one summary per window, tagged with the factual"
)
expect(
  filter(es, window == "attribution", scope == "A")$n_drop_f_only == 1L &&
    filter(es, window == "post_fork", scope == "A")$n_drop_f_only == 1L &&
    filter(es, window == "attribution", scope == "Africa")$n_pairs_kept == 4L,
  "the explosive month counts in both windows; Africa drops the pair"
)
es2 <- explosive_summary_cells(
  cm,
  attrib_cell_table(D("2026-05-01"), D("2026-07-01"), cap = 600),
  "neutral",
  "attrib"
)
expect(all(es2$n_pairs_kept == es2$n_pairs), "a higher cap drops nothing")
expect(
  !nrow(explosive_summary_cells(
    cm,
    attrib_cell_table(D("2026-05-01"), cap = 200, modes = "none"),
    "neutral",
    "x"
  )),
  "no drop mode: no summary"
)
rp <- function(run, aset, cases) {
  tibble(
    cell = "attribution",
    enso_group = "A",
    unit = "U",
    country = "AAA",
    run = run,
    analog_set = aset,
    replicate_id = 1:4,
    cases = cases,
    onsets = 1
  )
}
ie <- attrib_information_explosive(
  rp("forecast", NA, c(300, 10, 10, 10)),
  bind_rows(
    rp("forecast", NA, c(10, 300, 10, 10)),
    rp("counterfactual", "neutral", c(5, 5, 5, 300))
  ),
  100,
  "neutral"
)
expect(
  setequal(ie$comparison, c("stitched_vs_cf", "stitched_vs_may")),
  "both comparisons"
)
expect(
  filter(ie, comparison == "stitched_vs_may", scope == "A")$n_drop_f_only ==
    1L &&
    filter(ie, comparison == "stitched_vs_may", scope == "A")$n_drop_cf_only ==
      1L &&
    filter(ie, comparison == "stitched_vs_cf", scope == "A")$n_drop_cf_only ==
      1L,
  "stitched vs May: the May factual is the reference arm"
)

# [3] the cap and modes reach every cell from the context: modes {drop, none}
#     / {none} / an unknown mode refused; a context without a cap is refused
#     (no silent constant); the "none" cell keeps every replicate.
cat("# [3] cap and modes from the context\n")
ctx <- list2env(list(EXPLOSIVE_CAP = 123, EXPLOSIVE_MODES = c("none", "drop")))
cl <- attrib_cells(ctx, D("2026-05-01"), D("2026-06-01"))
expect(
  identical(cl$explosive_mode, rep(c("none", "drop"), 2)) &&
    identical(cl$cap, rep(c(Inf, 123), 2)) &&
    identical(cl$explosive_guard, rep(c(FALSE, TRUE), 2)),
  "attrib_cells(): the context's cap and modes, in order"
)
expect(
  identical(unique(event_cells(ctx, D("2015-03-01"), "ev")$cap), c(Inf, 123)),
  "event cells too"
)
ctx$EXPLOSIVE_MODES <- "none"
expect(
  identical(attrib_cells(ctx, D("2026-05-01"))$explosive_guard, FALSE),
  "modes [none]: one cell per window"
)
expect(
  err(attrib_cell_table(D("2026-05-01"), cap = 1, modes = c("drop", "winsor"))),
  "refused: unknown mode (no winsor)"
)
expect(
  err(attrib_cells(list2env(list()), D("2026-05-01"))),
  "refused: no cap in the context"
)
expect(
  err(drop_explosive_replicates(tibble(cum_cases = 1))),
  "drop_explosive_replicates(): the cap is required"
)
am <- cm %>%
  transmute(
    run,
    analog_set,
    enso_group,
    adm1_unit = paste0(country, "_1"),
    replicate_id,
    month,
    group_horizon_end,
    cases = n_cases,
    outbreaks = n_outbreaks
  )
pop <- tibble(adm1_unit = c("AAA_1", "BBB_1"), pop = c(10, 20))
mw <- suppressMessages(attrib_metrics_windows(
  am,
  cm,
  pop,
  c(A = "Ay", B = "Be"),
  "neutral",
  attrib_cell_table(D("2026-05-01"), cap = 200),
  50,
  .05,
  1,
  n_boot = 0L
))
na <- filter(mw$af_adm1, adm1_unit == "AAA_1", metric == "cases")
expect(
  na$n_pairs[na$explosive_guard] == 4L && na$n_pairs[!na$explosive_guard] == 5L,
  "drop cell: the explosive pair dropped; none cell: every pair kept"
)

# [4] Figure 4a's country ratios drop the GROUP-level explosive replicates
#     (country_cum_kept()), not the country's own cap.
cat("# [4] Figure 4a's country ratios: group-level keys\n")
cm4 <- cm %>%
  filter(month == months[1]) %>%
  bind_rows(
    tibble(
      run = "forecast",
      analog_set = NA,
      country = "CCC",
      replicate_id = 1:5,
      month = months[1],
      enso_group = "A",
      group_horizon_end = D("2026-07-01"),
      n_cases = c(150, 1, 1, 1, 1),
      n_outbreaks = 1L
    ),
    tibble(
      run = "counterfactual",
      analog_set = "neutral",
      country = "CCC",
      replicate_id = 1:5,
      month = months[1],
      enso_group = "A",
      group_horizon_end = D("2026-07-01"),
      n_cases = 1,
      n_outbreaks = 1L
    )
  )
kc <- country_cum_kept(cm4, D("2026-05-01"), 160)
# Group A replicate 1: AAA 20 + CCC 150 = 170 >= 160 -> dropped from AAA and
# CCC alike, though neither country reaches the cap on its own.
expect(
  !any(kc$replicate_id == 1 & kc$enso_group == "A" & kc$run == "forecast") &&
    any(kc$replicate_id == 1 & kc$enso_group == "B"),
  "a group-level explosive replicate leaves every country of its group only"
)
expect(
  nrow(country_cum_kept(cm4, D("2026-05-01"), Inf)) ==
    nrow(window_cum_per_rep(cm4, D("2026-05-01"), c("country", "enso_group"))),
  "cap Inf: every country row kept"
)

# [5] stamp back-compat: an artefact stamped before the keys is fresh at cap
#     2e5 and modes [drop, none], stale at another cap.
cat("# [5] stamp back-compat\n")
defs <- get_default_options()$results
expect(
  identical(defs$thresholds$explosive_case_cap, 2e5) &&
    identical(defs$derivation$explosive_modes, c("drop", "none")) &&
    identical(defs$style$attrib_headline$explosive, "drop"),
  "defaults: cap 2e5, modes [drop, none], headline drop"
)
st <- list(results = defs[.RESULTS_STAMPED_BLOCKS])
st$results$thresholds$explosive_case_cap <- NULL
st$results$derivation$explosive_modes <- NULL
expect(
  !length(.results_decision_diff(st, defs)),
  "an artefact stamped before the keys is fresh at 2e5 / [drop, none]"
)
r1 <- defs
r1$thresholds$explosive_case_cap <- 1e5
expect(length(.results_decision_diff(st, r1)) == 1L, "another cap re-derives")
r2 <- defs
r2$style$attrib_headline$explosive <- "none"
expect(
  !length(.results_decision_diff(st, r2)),
  "the headline mode is style (unstamped)"
)

# [6] the per-gid monthly cap (results.thresholds.explosive_gid_month_cap):
#     window_cum_per_rep() carries the largest gid-month; a replicate is
#     explosive at the group window cap OR the per-gid cap; Inf disables the
#     per-gid cap; a finite per-gid cap without max_gid_cases is refused; S29
#     counts it; cells carry it; stamp back-compat at Inf.
cat("# [6] per-gid monthly cap\n")
# Group A, 3 replicates x 2 months, 2 countries; max_gid_cases is the largest
# gid of the country-month. r1: a 90-case gid-month (total 100); r2: total 150
# spread thin (gids <= 20); r3: small.
cm6 <- tibble(
  run = "forecast",
  analog_set = NA_character_,
  enso_group = "A",
  group_horizon_end = D("2026-06-01"),
  replicate_id = rep(1:3, each = 4),
  month = rep(D(c("2026-05-01", "2026-05-01", "2026-06-01", "2026-06-01")), 3),
  country = rep(c("AAA", "BBB"), 6),
  n_cases = c(90, 5, 3, 2, 40, 40, 35, 35, 2, 2, 2, 2),
  n_outbreaks = 1L,
  max_gid_cases = c(90, 5, 3, 2, 20, 20, 18, 18, 1, 1, 1, 1)
)
w6 <- window_cum_per_rep(cm6, D("2026-05-01"), "enso_group")
expect(
  identical(w6$cum_cases, c(100, 150, 8)) &&
    identical(w6$max_gid_cases, c(90, 20, 1)),
  "window_cum_per_rep: group totals and the window's largest gid-month"
)
expect(
  identical(drop_explosive_replicates(w6, 120)$replicate_id, c(1L, 3L)) &&
    identical(
      drop_explosive_replicates(
        w6,
        c(case = 120, gid_month = Inf)
      )$replicate_id,
      c(1L, 3L)
    ),
  "a bare cap / gid cap Inf: the group window rule alone"
)
expect(
  identical(
    drop_explosive_replicates(w6, c(case = 1e6, gid_month = 50))$replicate_id,
    c(2L, 3L)
  ),
  "per-gid cap: the 90-case gid-month drops r1, the thin 150 total stays"
)
expect(
  identical(
    drop_explosive_replicates(w6, c(case = 120, gid_month = 50))$replicate_id,
    3L
  ),
  "either cap drops"
)
expect(
  identical(
    explosive_replicate_keys(
      cm6,
      D("2026-05-01"),
      c(case = 1e6, gid_month = 50)
    )$replicate_id,
    1L
  ),
  "explosive_replicate_keys: the per-gid explosive replicate"
)
expect(
  err(drop_explosive_replicates(
    select(w6, -max_gid_cases),
    c(case = 1e6, gid_month = 50)
  )),
  "a finite per-gid cap without max_gid_cases is refused"
)
expect(
  nrow(drop_explosive_replicates(
    select(w6, -max_gid_cases),
    c(case = 1e6, gid_month = Inf)
  )) ==
    3L,
  "no max_gid_cases is fine without a per-gid cap"
)
cum6 <- bind_rows(
  mutate(w6, run = "forecast", analog_set = NA_character_),
  mutate(w6, run = "counterfactual", analog_set = "neutral", max_gid_cases = 1)
)
s6 <- explosive_pair_summary(cum6, c(case = 1e6, gid_month = 50), "neutral")
s6a <- filter(s6, scope == "A")
expect(
  s6a$n_expl_f == 1L &&
    s6a$n_expl_cf == 0L &&
    s6a$n_drop_f_only == 1L &&
    s6a$n_pairs_kept == 2L &&
    s6a$cap == 1e6 &&
    s6a$gid_month_cap == 50,
  "S29 counts the per-gid explosive replicate and records both caps"
)
ct6 <- attrib_cell_table(D("2026-05-01"), cap = c(case = 1e6, gid_month = 2e5))
expect(
  identical(ct6$cap, c(1e6, Inf)) &&
    identical(ct6$gid_cap, c(2e5, Inf)) &&
    identical(.cell_cap(ct6[1, ]), c(case = 1e6, gid_month = 2e5)),
  "cells: drop carries both caps, none neither"
)
ctx6 <- list(EXPLOSIVE_CAP = 1e6, EXPLOSIVE_GID_CAP = 2e5)
expect(
  identical(.ctx_explosive_cap(ctx6), c(case = 1e6, gid_month = 2e5)),
  ".ctx_explosive_cap: both caps"
)
expect(
  is.infinite(defs$thresholds$explosive_gid_month_cap),
  "default per-gid cap: Inf (none)"
)
st6 <- list(results = defs[.RESULTS_STAMPED_BLOCKS])
st6$results$thresholds$explosive_gid_month_cap <- NULL
expect(
  !length(.results_decision_diff(st6, defs)),
  "an artefact stamped before the per-gid cap is fresh at Inf"
)
r6 <- defs
r6$thresholds$explosive_gid_month_cap <- 2e5
expect(
  length(.results_decision_diff(st6, r6)) == 1L,
  "a finite per-gid cap re-derives"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
