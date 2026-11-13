# test_loo_lags_helpers.R -- synthetic checks for the per-group lag resolver and
# the grouped-config helpers shared by the figure scripts.
#
# HARD gate (see analysis/tests/run_all.R). No fit, no data, no MCMC.
#   * loo_compare_ranking_files(): the filename DATE beats mtime (a git
#     checkout / cp without -p resets mtime, so mtime alone can flip a pick).
#   * best_lags_by_group(): equals select_best_lag() on the latest ranking,
#     differs from the raw rank-1 lag on a ranking built so they disagree,
#     honours overrides, and drops (with a warning) a group with no cache.
#   * groupedLagConfigPath() / detect_model_from_yaml().
#   * rankings carrying the Pareto-k columns select the same lag.
#
# Run: Rscript analysis/tests/test_loo_lags_helpers.R

source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/loo_lags.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# loo_lags helpers (latest ranking / per-group best lag / grouped config path)\n")

# Ranking whose raw rank-1 lag (1) differs from the SE-window pick (3):
# lag 1 is best; lags 2 and 3 are within 2 SE of it; lag 4 is not.
mk_ranking <- function(elpd) {
  best <- max(elpd)
  tibble::tibble(single_lag = seq_along(elpd), elpd_loo = elpd, se_elpd = 5,
                 elpd_diff = elpd - best,
                 se_diff = ifelse(elpd == best, 0, 1), p_loo = 10) |>
    dplyr::mutate(rank = rank(-elpd_loo, ties.method = "first")) |>
    dplyr::arrange(rank)
}
rk_a <- mk_ranking(c(-100, -101, -102, -130))
rk_b <- mk_ranking(c(-100, -101, -131, -130))   # newer file: SE-window pick = 2

tmp   <- tempfile("loo_helpers_"); dir.create(tmp)
opt   <- list(outdir = tmp)
batch <- "fake_occ_grouped_enso_singlelag"
dir   <- file.path(tmp, "loo_compare", batch); dir.create(dir, recursive = TRUE)

save_rk <- function(rk, group, date) {
  f <- file.path(dir, sprintf("loo_compare_lags_%s_%s.rds", group, date))
  saveRDS(list(ranking = rk), f); f
}
invisible(f_old <- save_rk(rk_a, "enso_east", "20260101"))
f_new <- save_rk(rk_b, "enso_east", "20260301")
# mtime deliberately inverted: the OLD-dated file is the most recently modified.
Sys.setFileTime(f_new, Sys.time() - 86400 * 30)
Sys.setFileTime(f_old, Sys.time())

lf <- loo_compare_ranking_files(batch, group = "enso_east", opt = opt)
expect(nrow(lf) == 1L && basename(lf$file) == basename(f_new),
       "latest ranking is picked by filename date, not mtime")
expect(identical(latest_ranking_file(batch, group = "enso_east", opt = opt), f_new),
       "latest_ranking_file() agrees")

# best_lags_by_group == select_best_lag() on the latest ranking
bl <- best_lags_by_group(batch, "enso_east", "occurrence", opt = opt)
expect(identical(unname(bl), select_best_lag(rk_b, "occurrence")$lag),
       "best_lags_by_group() == select_best_lag() on the latest ranking")

# the canonical rule differs from the raw rank-1 row (ranking$single_lag[1]) on rk_a
invisible(save_rk(rk_a, "enso_west", "20260301"))
bl_w <- best_lags_by_group(batch, "enso_west", "occurrence", opt = opt)
expect(unname(bl_w) == 3L && rk_a$single_lag[1L] == 1L,
       "SE-window lag (3) differs from raw rank-1 lag (1)")
bl_w0 <- best_lags_by_group(batch, "enso_west", "occurrence", se_mult = 0, opt = opt)
expect(unname(bl_w0) == 1L, "se_mult = 0 reproduces the raw rank-1 lag")

# overrides win and need no file; a group with no cache is dropped with a warning
bl_o <- best_lags_by_group(batch, c("enso_west", "enso_horn"), "occurrence",
                           overrides = c(enso_horn = 7L), opt = opt)
expect(identical(bl_o, c(enso_west = 3L, enso_horn = 7L)), "overrides win, no file needed")
warned <- FALSE
bl_m <- withCallingHandlers(
  best_lags_by_group(batch, c("enso_west", "enso_nowhere"), "occurrence", opt = opt),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
expect(warned && identical(names(bl_m), "enso_west"), "missing group dropped with a warning")

# downstream track honours the occurrence floor when asked
bl_s <- best_lags_by_group(batch, "enso_west", "size", occ_lags = c(enso_west = 4L),
                           enforce_occ_floor = TRUE, opt = opt)
expect(unname(bl_s) >= 4L, "enforce_occ_floor lifts a downstream lag to the occurrence floor")

# ENSO leads (-11..11 sweep): ranked, never selected (min_lag = 0) ---------------
# lags -2..3; the lead -1 is rank 1, lags 0..2 within 2 SE of it, lag 3 not.
mk_lead_ranking <- function() {
  lag  <- -2:3
  elpd <- c(-104, -100, -101, -101.5, -102, -130)
  se_e <- rep(5, length(lag))
  best <- max(elpd)
  tibble::tibble(single_lag = lag, elpd_loo = elpd, se_elpd = se_e,
                 elpd_diff = elpd - best, se_diff = ifelse(elpd == best, 0, 1.5),
                 p_loo = 10) |>
    dplyr::mutate(rank = rank(-elpd_loo, ties.method = "first")) |>
    dplyr::arrange(rank)
}
rk_lead <- mk_lead_ranking()
sel_lead <- select_best_lag(rk_lead, "occurrence")
expect(sel_lead$lag >= 0L, "a lead at rank 1 is never selected (min_lag = 0)")
expect(identical(sel_lead$lag, 2L), "with a lead at rank 1, the longest eligible lag within SE is picked")
expect(identical(select_best_lag(rk_lead, "occurrence", min_lag = -Inf)$lag, 2L),
       "min_lag = -Inf: the SE-window rule over all lags (max lag in window) is unchanged")
expect(select_best_lag(rk_lead, "duration", enforce_occ_floor = FALSE)$lag >= 0L,
       "downstream tracks honour min_lag too")
# the window stays anchored on the lead rank-1 fit (paired SEs), filtered to
# lags >= 0 -- not the wide quadrature-SE fallback, which picked the max lag
rk_wide <- rk_lead; rk_wide$se_elpd <- 100   # total se_elpd >> paired se_diff
expect(identical(select_best_lag(rk_wide, "occurrence")$lag, 2L),
       "a wide total se_elpd does not widen the window (lag 3, outside it, is not picked)")
# no lag >= 0 within the window of a lead: the best-ranked eligible lag
rk_far <- rk_lead
rk_far$elpd_diff[rk_far$single_lag >= 0] <- -50; rk_far$se_diff[rk_far$single_lag >= 0] <- 2
rk_far$elpd_diff[rk_far$single_lag == 1] <- -40
rk_far$elpd_loo <- -100 + rk_far$elpd_diff
rk_far$rank <- rank(-rk_far$elpd_loo, ties.method = "first")
expect(identical(select_best_lag(rk_far, "occurrence")$lag, 1L),
       "no eligible lag in the lead's window: the best-ranked eligible lag")
# min_lag = 0 is a no-op on every lead-free ranking
for (rk in list(rk_a, rk_b))
  expect(identical(select_best_lag(rk, "occurrence")$lag,
                   select_best_lag(rk, "occurrence", min_lag = -Inf)$lag),
         "min_lag = 0 is a no-op on a ranking without leads")
writeLines("occurrence:\n  backend: brms_cb",
           file.path(tmp, "fake_occ_grouped_enso_singlelag_sl-m3_grp-ensowest.yml"))
expect(grepl("_sl-m3_grp-ensowest\\.yml$", groupedLagConfigPath(tmp, "enso_west", -3L) %||% ""),
       "groupedLagConfigPath resolves a lead config (_sl-m3_)")

# groupedLagConfigPath / detect_model_from_yaml
cdir <- file.path(tmp, "configs", "fake_occ_grouped_enso_singlelag"); dir.create(cdir, recursive = TRUE)
writeLines("occurrence:\n  backend: brms_cb", file.path(cdir, "fake_occ_grouped_enso_singlelag_sl-3_grp-ensowest.yml"))
writeLines("size:\n  backend: brms_cb",       file.path(cdir, "fake_size_x_sl-3_grp-ensoeast.yml"))
writeLines("batch: x",                        file.path(cdir, "manifest.yml"))
p <- groupedLagConfigPath(cdir, "enso_west", 3L)
expect(!is.null(p) && grepl("_sl-3_grp-ensowest\\.yml$", p), "groupedLagConfigPath finds the config")
expect(is.null(suppressWarnings(groupedLagConfigPath(cdir, "enso_west", 9L))),
       "groupedLagConfigPath returns NULL (warns) for a missing lag")
expect(identical(detect_model_from_yaml(p), "occurrence"), "detect_model_from_yaml: occurrence")
expect(identical(detect_model_from_yaml(groupedLagConfigPath(cdir, "enso_east", 3L)), "size"),
       "detect_model_from_yaml: size")
expect(inherits(try(detect_model_from_yaml(file.path(cdir, "manifest.yml")), silent = TRUE), "try-error"),
       "detect_model_from_yaml stops on a config with no model block")


# best_elpd rule + lag_selection_cfg()
cat("\n[best_elpd] best-ranked lag, leads allowed, downstream floor occ - offset\n")
rk_b <- data.frame(rank = 1:6, single_lag = c(-6L, 7L, 9L, -10L, 2L, 11L),
                   elpd_diff = c(0, -3, -5, -7, -9, -12), se_diff = c(0, 5, 5, 5, 5, 5))
expect(identical(select_best_lag_best_elpd(rk_b)$lag, -6L), "occurrence: rank-1 lead selected")
expect(identical(select_best_lag_best_elpd(rk_b, min_lag = 0)$lag, 7L), "min_lag = 0: best non-negative lag")
expect(identical(select_best_lag_best_elpd(rk_b, occ_lag = -6L, occ_floor_offset = 4)$lag, -6L),
       "downstream floor occ - 4 = -10 keeps -6")
expect(identical(select_best_lag_best_elpd(rk_b, occ_lag = 9L, occ_floor_offset = 4)$lag, 7L),
       "downstream floor 9 - 4 = 5 -> best lag >= 5 is 7")
expect(identical(select_best_lag_best_elpd(rk_b, occ_lag = 9L, occ_floor_offset = NULL)$lag, -6L),
       "no offset: no occurrence floor")
expect(identical(select_best_lag(rk_b, "duration", occ_lag = 9L, rule = "best_elpd", occ_floor_offset = 4,
                                 min_lag = -Inf)$lag, 7L), "select_best_lag(rule = best_elpd) dispatches")
expect(identical(select_best_lag(rk_b, "occurrence", rule = "longest_within_se")$lag,
                 select_best_lag(rk_b, "occurrence")$lag), "default rule unchanged")
c0 <- lag_selection_cfg(NULL)
expect(identical(c0$rule, "longest_within_se") && identical(c0$min_lag, 0L) && is.null(c0$occ_floor_offset) &&
         identical(c0$se_mult, 2), "lag_selection_cfg(NULL) = the current rule")
c1 <- lag_selection_cfg(list(rule = "best_elpd", se_mult = 2, min_lag = "rule_default", occ_floor_offset = 4))
expect(identical(c1$min_lag, -Inf) && identical(c1$occ_floor_offset, 4), "best_elpd defaults: min_lag -Inf, offset kept")
expect(is.null(lag_selection_cfg(list(rule = "longest_within_se", occ_floor_offset = "none"))$occ_floor_offset),
       "occ_floor_offset 'none' = no floor")
expect(inherits(try(lag_selection_cfg(list(rule = "rank1")), silent = TRUE), "try-error"), "unknown rule refused")

# Pareto-k columns: a ranking cache that carries them (as 05
# writes) resolves to the same per-group lag, and read_latest_ranking() keeps them.
cat("\n[pareto-k] rankings with the Pareto-k columns\n")
k_cols <- as.data.frame(setNames(replicate(length(PARETO_K_COLS), rep(0.5, nrow(rk_a)), simplify = FALSE),
                                 PARETO_K_COLS))
invisible(save_rk(dplyr::bind_cols(rk_a, k_cols), "enso_south", "20260301"))
invisible(save_rk(rk_a, "enso_nga", "20260301"))
for (rl in c("longest_within_se", "best_elpd"))
  expect(identical(unname(best_lags_by_group(batch, "enso_south", "occurrence", opt = opt, rule = rl, min_lag = -Inf)),
                   unname(best_lags_by_group(batch, "enso_nga", "occurrence", opt = opt, rule = rl, min_lag = -Inf))),
         sprintf("best_lags_by_group(rule = %s): same lag with and without the Pareto-k columns", rl))
rk_k <- read_latest_ranking(batch, "enso_south", opt = opt)
expect(all(PARETO_K_COLS %in% names(rk_k)) &&
         identical(attr(ranking_pareto(rk_k), "pareto_source"), "ranking"),
       "read_latest_ranking() keeps the k columns; ranking_pareto() passes them through")
expect(identical(attr(ranking_pareto(read_latest_ranking(batch, "enso_nga", opt = opt)), "pareto_source"), "none"),
       "an older ranking without them reads as NA columns (pareto_source 'none')")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
