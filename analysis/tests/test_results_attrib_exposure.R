# test_results_attrib_exposure.R — HARD: report_attrib_exposure and the
#   helpers extracted for it, each against a literal copy of the old code.
#
# Usage:
#   Rscript analysis/tests/test_results_attrib_exposure.R
#
# Details:
#   Quick-safe: no fit, no data, no simulation output. report_attrib_exposure
#   (utils/results_derive.R) is the population in ADM1s with more
#   ENSO-attributable cases, forecast vs counterfactual. Pins the extracted
#   helpers (identical() to the pre-extraction code), the classification
#   rules, attrib_exposure_tables() end to end and the window cells.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/results_core.R") # base, figures_cb (enso_conds_named), registries
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# report_attrib_exposure + extracted helpers\n")

# ---- literal pre-extraction code -------------------------------------------
legacy_pop_summary_rows <- function(mat, enso_conds) {
  qs <- apply(mat, 1L, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE) # 3 x C
  tibble::tibble(
    enso_label = names(enso_conds),
    enso_val = unname(enso_conds),
    pop_mean = rowMeans(mat),
    pop_lo = qs[1L, ],
    pop_med = qs[2L, ],
    pop_hi = qs[3L, ]
  )
}
legacy_mask_pop_rows <- function(masks, dr, point = NULL) {
  sc <- enso_scopes(dr$universe)
  rows <- list()
  for (t in names(masks)) {
    for (k in seq_len(nrow(sc))) {
      idx <- sc$idx[[k]]
      mk <- masks[[t]][idx, , , drop = FALSE]
      r <- dplyr::bind_cols(
        .scope_cols(dr, sc, k),
        legacy_pop_summary_rows(pop_mass_by_draw(mk, dr$pop0[idx]), dr$conds)
      ) %>%
        dplyr::mutate(track = t, n_adm1_mean = rowMeans(colSums(mk, dims = 1L)))
      if (!is.null(point)) {
        pt <- point[[t]][idx, , drop = FALSE]
        r <- dplyr::mutate(
          r,
          n_adm1_point = colSums(pt),
          pop_point = colSums(pt * dr$pop0[idx])
        )
      }
      rows[[length(rows) + 1L]] <- r
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::relocate(
      track,
      level,
      scope,
      enso_group,
      region,
      country,
      enso_label,
      enso_val
    ) %>%
    add_pop_shares()
}
legacy_prob_grid_rows <- function(masks, dr, cutoffs, credible_union = NULL) {
  sc <- enso_scopes(dr$universe)
  rows <- list()
  prob <- lapply(masks, function(m) apply(m, c(1L, 2L), mean)) # U x C each
  hits <- lapply(stats::setNames(names(masks), names(masks)), function(t) {
    function(p) prob[[t]] >= p
  })
  if (length(credible_union)) {
    stopifnot(all(credible_union %in% names(masks)))
    hits$any_credible <- function(p) {
      Reduce(`|`, lapply(credible_union, function(t) prob[[t]] >= p))
    }
  }
  for (t in names(hits)) {
    for (k in seq_len(nrow(sc))) {
      for (p in cutoffs) {
        idx <- sc$idx[[k]]
        hit <- hits[[t]](p)[idx, , drop = FALSE]
        rows[[length(rows) + 1L]] <- dplyr::bind_cols(
          .scope_cols(dr, sc, k),
          tibble::tibble(
            enso_label = names(dr$conds),
            enso_val = unname(dr$conds),
            cutoff = p,
            n_adm1 = unname(colSums(hit)),
            pop = unname(colSums(hit * dr$pop0[idx]))
          )
        ) %>%
          dplyr::mutate(track = t)
      }
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::relocate(
      track,
      level,
      scope,
      enso_group,
      region,
      country,
      enso_label,
      enso_val,
      cutoff
    ) %>%
    add_pop_shares()
}
legacy_cum_group <- function(cm, attrib_start) {
  cm %>%
    dplyr::filter(month >= attrib_start, month <= group_horizon_end) %>%
    dplyr::group_by(enso_group, run, analog_set, replicate_id) %>%
    dplyr::summarise(
      cum_cases = sum(n_cases),
      cum_outbreaks = sum(n_outbreaks),
      .groups = "drop"
    ) %>%
    dplyr::filter(cum_cases < 2e5)
}
legacy_cum_country <- function(cm, attrib_start) {
  cm %>%
    dplyr::filter(month >= attrib_start, month <= group_horizon_end) %>%
    dplyr::group_by(country, enso_group, run, analog_set, replicate_id) %>%
    dplyr::summarise(
      cum_cases = sum(n_cases),
      cum_outbreaks = sum(n_outbreaks),
      .groups = "drop"
    ) %>%
    dplyr::filter(cum_cases < 2e5)
}
legacy_pair <- function(adm1_month_both, vc, aset) {
  cum_adm1 <- adm1_month_both %>%
    dplyr::filter(run == "forecast" | analog_set == aset) %>%
    dplyr::group_by(run, adm1_unit, replicate_id) %>%
    dplyr::summarise(value = sum(.data[[vc]]), .groups = "drop")
  dplyr::inner_join(
    dplyr::filter(cum_adm1, run == "forecast") %>%
      dplyr::transmute(adm1_unit, replicate_id, f = value),
    dplyr::filter(cum_adm1, run == "counterfactual") %>%
      dplyr::transmute(adm1_unit, replicate_id, cf = value),
    by = c("adm1_unit", "replicate_id")
  )
}

# ---- synthetic per-draw universe (as test_results_sensitive_pop.R) -----------
conds <- enso_conds_named(c(-1, 1))
set.seed(1)
arr <- array(
  rnorm(3 * 2 * 4, sd = 0.2),
  dim = c(3, 2, 4),
  dimnames = list(c("KEN_1", "KEN_2", "UGA_1"), names(conds), NULL)
)
uni <- tibble::tibble(
  adm1_gid = c("KEN_1", "KEN_2", "UGA_1"),
  country = c("KEN", "KEN", "UGA"),
  enso_group = c("enso_east", "enso_east", "enso_east"),
  region = "East",
  pop = c(100, NA, 50),
  pop_missing = c(FALSE, TRUE, FALSE)
)
dr <- list(
  tracks = c("a", "b"),
  conds = conds,
  n_draws = 4L,
  universe = uni,
  U = uni$adm1_gid,
  pop0 = c(100, 0, 50),
  logrr = list(a = arr, b = -arr)
)
masks <- list(a = arr > 0, b = arr > 0.1)
point <- list(a = arr[,, 1] > 0, b = arr[,, 1] > 0.1)

# [1] refactor gates (identical() to the pre-extraction code):
#     window_cum_per_rep() + drop_explosive_replicates() == attribution's two
#     inline cumulative blocks (incl. a replicate at exactly the cap);
#     explosive_replicate_keys() is their exact complement;
#     pair_adm1_cum() == report_forecast's inline pop_shift pairing;
#     mask_pop_rows() (with / without `point`), prob_grid_rows() and
#     .pop_summary_rows() == their pre-extraction bodies;
cat("# [1] refactor gates\n")
expect(
  identical(mask_pop_rows(masks, dr), legacy_mask_pop_rows(masks, dr)),
  "mask_pop_rows() == pre-extraction body"
)
expect(
  identical(
    mask_pop_rows(masks, dr, point),
    legacy_mask_pop_rows(masks, dr, point)
  ),
  "mask_pop_rows(point =) == pre-extraction body"
)
expect(
  identical(
    prob_grid_rows(masks, dr, c(0.5, 1), credible_union = c("a", "b")),
    legacy_prob_grid_rows(masks, dr, c(0.5, 1), credible_union = c("a", "b"))
  ),
  "prob_grid_rows() (credible_union) == pre-extraction body"
)
expect(
  identical(
    prob_grid_rows(masks, dr, 0.5),
    legacy_prob_grid_rows(masks, dr, 0.5)
  ),
  "prob_grid_rows() == pre-extraction body"
)
mat <- matrix(c(1, 5, 2, 7, 3, 9, 4, 11), 2)
expect(
  identical(.pop_summary_rows(mat, conds), legacy_pop_summary_rows(mat, conds)),
  ".pop_summary_rows() == pre-extraction body"
)

# country_month: 2 groups (KEN east, UGA horn), 4 replicates, 3 months, forecast + 1 set
a0 <- as.Date("2026-05-01")
cm <- expand.grid(
  replicate_id = 1:4,
  month = as.Date(c("2026-04-01", "2026-05-01", "2026-06-01")),
  country = c("KEN", "UGA"),
  run = c("forecast", "counterfactual"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    analog_set = ifelse(run == "forecast", NA_character_, "all"),
    enso_group = ifelse(country == "KEN", "enso_east", "enso_horn"),
    group_horizon_end = as.Date("2026-06-01"),
    n_cases = 1000 * replicate_id,
    n_outbreaks = replicate_id %% 3L
  ) %>%
  tibble::as_tibble()
# replicate 4 of KEN forecast at exactly the cap in-window; replicate 3 above it
cm$n_cases[
  cm$country == "KEN" &
    cm$run == "forecast" &
    cm$replicate_id == 4 &
    cm$month >= a0
] <- 1e5
cm$n_cases[
  cm$country == "KEN" &
    cm$run == "forecast" &
    cm$replicate_id == 3 &
    cm$month >= a0
] <- 2e5
cm_f <- dplyr::mutate(
  cm,
  run = factor(run, levels = c("forecast", "counterfactual"))
)
expect(
  identical(
    drop_explosive_replicates(window_cum_per_rep(cm_f, a0, "enso_group"), 2e5),
    legacy_cum_group(cm_f, a0)
  ),
  "group cumulative + explosive drop == attribution's inline block (factor run)"
)
expect(
  identical(
    drop_explosive_replicates(
      window_cum_per_rep(cm_f, a0, c("country", "enso_group")),
      2e5
    ),
    legacy_cum_country(cm_f, a0)
  ),
  "country cumulative + explosive drop == attribution's inline block"
)
keys <- explosive_replicate_keys(cm, a0, 2e5)
expect(
  nrow(keys) + nrow(legacy_cum_group(cm, a0)) ==
    nrow(window_cum_per_rep(cm, a0, "enso_group")),
  "explosive_replicate_keys() is the exact complement"
)
expect(
  identical(sort(keys$replicate_id), c(3L, 4L)) &&
    all(keys$enso_group == "enso_east") &&
    all(keys$run == "forecast"),
  "a replicate AT the cap is explosive (cum_cases < cap kept), only its group/arm"
)

am <- cm %>%
  dplyr::transmute(
    replicate_id,
    adm1_unit = paste0(country, "_1"),
    month,
    cases = n_cases,
    outbreaks = n_outbreaks,
    run,
    analog_set,
    enso_group
  )
expect(
  identical(pair_adm1_cum(am, "cases", "all"), legacy_pair(am, "cases", "all")),
  "pair_adm1_cum() == report_forecast's inline pop_shift pairing"
)
expect(
  identical(
    pair_adm1_cum(am, "outbreaks", "all"),
    legacy_pair(am, "outbreaks", "all")
  ),
  "pair_adm1_cum() (onsets) == inline pairing"
)

# [2] rules: childs_hit(), attrib_exposure_classes() (ties disagree,
#     mean_cf == 0, exact 0.66 / 0.80 boundaries, m = 1), attrib_class_hits(),
#     boot_weights(), boot_adm1_stats() with identity weights (means and
#     medians equal the point statistics);
cat("# [2] rules\n")
expect(
  identical(
    childs_hit(c(3, 2, 1, 5, 0, 1), c(2, 1, 0, 0, 0, 1), 0.5),
    c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  ),
  "childs_hit: >= (1+m) x cf; positive vs zero counts (infinite ratio), 0 vs 0 does not"
)
pr <- tibble::tibble(
  adm1_unit = rep(c("A", "B", "C", "D"), each = 5L),
  replicate_id = rep(1:5, 4),
  f = c(
    30,
    30,
    30,
    30,
    10,
    10,
    10,
    10,
    10,
    10,
    5,
    0,
    0,
    0,
    0,
    30,
    30,
    30,
    1,
    1
  ),
  cf = c(
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    0,
    0,
    0,
    0,
    0,
    10,
    10,
    10,
    5,
    5
  )
)
cl <- attrib_exposure_classes(
  pr,
  m = 1,
  robust_exceed = 0.6,
  robust_agree = 0.8
)
expect(
  identical(cl$p_sign, c(0.8, 0, 0.2, 0.6)),
  "p_sign strict: ties disagree"
)
expect(
  identical(cl$p_exceed, c(0.8, 0, 0.2, 0.6)),
  "p_exceed at m = 1: f > 2 cf"
)
expect(
  identical(cl$class, c("robust", "none", "none", "conflicting")),
  "classes: robust at the exact 0.8 boundary; >= exceed but < agree = conflicting"
)
expect(
  identical(cl$ratio[3], Inf) && cl$childs[3],
  "mean_cf == 0 < mean_f: ratio Inf, counts as Childs"
)
expect(
  isTRUE(all.equal(cl$ratio[c(1, 2, 4)], c(26 / 10, 1, 92 / 40))),
  "ratio of means"
)
hh <- attrib_class_hits(
  dplyr::mutate(cl, analog_set = "all"),
  c("A", "B", "C", "D", "E"),
  c("all", "nina")
)
expect(
  identical(dim(hh$robust), c(5L, 2L)) &&
    !any(hh$robust[, "nina"]) &&
    !any(hh$paired["E", ]),
  "class hits: U x S; absent ADM1 / set is FALSE"
)
expect(
  !any(hh$undefined[, "all"]) &&
    identical(
      unname(hh$undefined_median[, "all"]),
      c(FALSE, FALSE, TRUE, FALSE, FALSE)
    ),
  "undefined = both statistics 0 (C: means 1 vs 0 defined, medians 0 vs 0 undefined)"
)
w1 <- .with_seed(7, boot_weights(5L, 10L))
w2 <- .with_seed(7, boot_weights(5L, 10L))
expect(
  identical(w1, w2) && all(colSums(w1) == 5L) && identical(dim(w1), c(5L, 10L)),
  "boot_weights: n x B counts summing to n, reproducible under the seed"
)
expect(
  identical(dim(boot_weights(0L, 4L)), c(0L, 4L)) &&
    identical(dim(boot_weights(1L, 3L)), c(1L, 3L)),
  "boot_weights: n = 0 and n = 1 keep the matrix shape"
)
set.seed(99)
r0 <- runif(1)
set.seed(99)
invisible(.with_seed(1, runif(5)))
expect(identical(runif(1), r0), ".with_seed restores the caller's RNG state")
mats <- list(
  F = matrix(pr$f, 4, byrow = TRUE),
  CF = matrix(pr$cf, 4, byrow = TRUE)
)
bs <- boot_adm1_stats(mats, matrix(1L, 5, 1))
expect(
  identical(childs_hit(bs$mean_f, bs$mean_cf, 1)[, 1], cl$childs),
  "boot_adm1_stats with identity weights: Childs (means) == the point hits"
)
expect(
  identical(bs$med_f[, 1], cl$med_f) &&
    identical(bs$med_cf[, 1], cl$med_cf) &&
    identical(childs_hit(bs$med_f, bs$med_cf, 1)[, 1], cl$childs_median),
  "boot_adm1_stats with identity weights: medians == the point medians"
)
w2 <- matrix(c(2L, 0L, 1L, 1L, 1L), 5, 1) # replicate 1 twice, replicate 2 dropped
bs2 <- boot_adm1_stats(mats, w2)
expect(
  isTRUE(all.equal(
    bs2$med_f[, 1],
    apply(mats$F[, c(1, 1, 3, 4, 5)], 1, median)
  )),
  "boot_adm1_stats: resample medians repeat replicates by their weights"
)
expect(
  identical(cl$ratio_median, c(3, 1, NA, 3)) &&
    identical(cl$childs_median, c(TRUE, FALSE, FALSE, TRUE)),
  "ratio of medians: robust to the tail (A: 3, D: 3 vs mean 2.3), median cf 0 -> NA"
)

# [3] attrib_exposure_tables() end to end on 2 groups x 3 ADM1 x 4 replicates;
cat("# [3] attrib_exposure_tables() end to end\n")
# 3 ADM1 (KEN_1, KEN_2 east; UGA_1 horn), 4 replicates; KEN forecast replicate 4
# explosive at group level (2e5 over KEN's two ADM1s)
am3 <- expand.grid(
  replicate_id = 1:4,
  adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
  run = c("forecast", "counterfactual"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    month = a0,
    analog_set = ifelse(run == "forecast", NA_character_, "all"),
    enso_group = ifelse(
      substr(adm1_unit, 1, 3) == "KEN",
      "enso_east",
      "enso_horn"
    ),
    group_horizon_end = as.Date("2026-06-01"),
    cases = dplyr::case_when(
      run == "counterfactual" ~ 10,
      adm1_unit == "KEN_1" ~ 30,
      adm1_unit == "KEN_2" ~ 12,
      TRUE ~ c(30, 30, 30, 0)[replicate_id]
    ),
    outbreaks = as.integer(cases > 10)
  ) %>%
  tibble::as_tibble()
am3$cases[
  am3$run == "forecast" & am3$replicate_id == 4 & am3$adm1_unit == "KEN_1"
] <- 2e5
cm3 <- am3 %>%
  dplyr::mutate(country = substr(adm1_unit, 1, 3)) %>%
  dplyr::group_by(
    replicate_id,
    month,
    country,
    run,
    analog_set,
    enso_group,
    group_horizon_end
  ) %>%
  dplyr::summarise(
    n_cases = sum(cases),
    n_outbreaks = sum(outbreaks),
    .groups = "drop"
  )
pop3 <- tibble::tibble(
  adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
  pop = c(100, 40, NA)
)
gs <- c(enso_east = "East", enso_horn = "Horn")
out <- attrib_exposure_tables(
  am3,
  cm3,
  pop3,
  gs,
  "all",
  a0,
  cap = 2e5,
  mags = c(0.5, 1),
  robust_exceed = 0.66,
  robust_agree = 0.8,
  n_boot = 50L
)
expect(
  identical(out$explosive$replicate_id, 4L) &&
    out$explosive$enso_group == "enso_east",
  "explosive replicate keyed to its group and arm"
)
a <- dplyr::filter(out$adm1, metric == "cases", magnitude == 0.5)
expect(
  identical(a$n_pairs, c(3L, 3L, 4L)),
  "explosive replicate dropped from every ADM1 of its group only"
)
expect(
  identical(a$class, c("robust", "none", "conflicting")) &&
    identical(a$childs, c(TRUE, FALSE, TRUE)),
  "per-ADM1 classes and Childs hits (UGA: 3 of 4 above 1.5x, 3 of 4 agree)"
)
ex <- dplyr::filter(out$exposure, metric == "cases", magnitude == 0.5)
afr <- dplyr::filter(ex, level == "continent")
reg <- dplyr::filter(ex, level == "region")
expect(
  isTRUE(all.equal(
    afr$pop[order(afr$class)],
    (reg %>%
      dplyr::group_by(class) %>%
      dplyr::summarise(p = sum(pop)) %>%
      dplyr::arrange(class))$p
  )),
  "Africa population = sum of regions, every class"
)
ch <- dplyr::filter(afr, class == "childs")
expect(
  ch$n_adm1 == 2L &&
    ch$pop == 100 &&
    ch$n_adm1 == ch$n_adm1_point &&
    ch$n_pop_missing == 1L,
  "Childs Africa: KEN_1 + UGA_1 (pop missing -> 0, counted)"
)
expect(isTRUE(all.equal(ch$pct, 100 / 140)), "pct = pop / pop_universe")
expect(
  ch$pop_lo <= ch$pop_point && ch$pop_point <= ch$pop_hi,
  "bootstrap interval brackets the point"
)
expect(
  all(is.na(
    dplyr::filter(afr, !class %in% c("childs", "childs_median"))$pop_lo
  )),
  "no interval on the IPCC class rows"
)
out2 <- attrib_exposure_tables(
  am3,
  cm3,
  pop3,
  gs,
  "all",
  a0,
  cap = 2e5,
  mags = c(0.5, 1),
  robust_exceed = 0.66,
  robust_agree = 0.8,
  n_boot = 50L
)
expect(identical(out, out2), "deterministic (fixed bootstrap seed)")
expect(
  identical(
    sort(unique(out$exposure$class)),
    c("childs", "childs_median", "conflicting", "none", "robust")
  ) &&
    setequal(unique(out$exposure$level), c("continent", "region", "country")),
  "five classes at continent / region / country level"
)
chm <- dplyr::filter(afr, class == "childs_median")
expect(
  chm$n_adm1 == chm$n_adm1_point &&
    !is.na(chm$pop_lo) &&
    chm$n_adm1_undefined == 0L,
  "childs_median rows: point counts + bootstrap interval"
)

# [4] attrib_exposure_windows(): window x explosive-guard cells.
cat("# [4] attrib_exposure_windows()\n")
# a later factual launch: post_fork keeps only months >= 2026-06-01
am4 <- dplyr::bind_rows(
  am3,
  dplyr::mutate(
    am3,
    month = as.Date("2026-06-01"),
    cases = ifelse(run == "forecast", 10, 10)
  )
)
cm4 <- am4 %>%
  dplyr::mutate(country = substr(adm1_unit, 1, 3)) %>%
  dplyr::group_by(
    replicate_id,
    month,
    country,
    run,
    analog_set,
    enso_group,
    group_horizon_end
  ) %>%
  dplyr::summarise(
    n_cases = sum(cases),
    n_outbreaks = sum(outbreaks),
    .groups = "drop"
  )
ow <- attrib_exposure_windows(
  am4,
  cm4,
  pop3,
  gs,
  "all",
  attrib_cell_table(a0, as.Date("2026-06-01"), cap = 2e5),
  mags = 0.5,
  robust_exceed = 0.66,
  robust_agree = 0.8,
  n_boot = 20L
)
expect(
  setequal(
    unique(paste(ow$exposure$window, ow$exposure$explosive_guard)),
    c(
      "attribution TRUE",
      "attribution FALSE",
      "post_fork TRUE",
      "post_fork FALSE"
    )
  ),
  "four window x guard cells"
)
pf <- dplyr::filter(
  ow$exposure,
  window == "post_fork",
  explosive_guard,
  level == "continent",
  class == "childs",
  metric == "cases"
)
expect(
  pf$n_adm1 == 0L,
  "post_fork: equal arms after the fork -> no Childs ADM1"
)
ng <- dplyr::filter(
  ow$adm1,
  window == "attribution",
  !explosive_guard,
  metric == "cases"
)
expect(all(ng$n_pairs == 4L), "explosive_guard FALSE keeps every replicate")
expect(
  all(ow$explosive$explosive_guard) &&
    all(ow$explosive$window == "attribution"),
  "explosive keys reported for the guarded cells (none post-fork here)"
)
same <- attrib_exposure_windows(
  am3,
  cm3,
  pop3,
  gs,
  "all",
  attrib_cell_table(a0, a0, cap = 2e5),
  mags = c(0.5, 1),
  robust_exceed = 0.66,
  robust_agree = 0.8,
  n_boot = 50L
)
strip <- function(d) dplyr::select(d, -window, -explosive_guard)
expect(
  identical(
    strip(dplyr::filter(
      same$exposure,
      window == "attribution",
      explosive_guard
    )),
    out$exposure
  ),
  "attribution x guarded slice == attrib_exposure_tables() (the v1 rows)"
)
expect(
  identical(
    strip(dplyr::filter(same$exposure, window == "post_fork", explosive_guard)),
    out$exposure
  ),
  "arms launched together: post_fork == attribution"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
