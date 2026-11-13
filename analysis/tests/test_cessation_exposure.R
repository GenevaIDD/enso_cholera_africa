# test_cessation_exposure.R -- simulator.cessation_exposure:
# the cessation month's monthly case exposure in the v1 kernel
# (simulate_one_tu) and the spatPomp C kernel (build_spatpomp_tu_c).
#
# HARD gate (see analysis/tests/run_all.R). Synthetic, no MCMC; section R reads
# the monthly fit panel when it is present and is skipped otherwise.
#
#   H  helpers: cessationExposureOf() / checkCessationExposure(),
#      .cessation_log_offset_adj().
#   G  golden (rule 19): a fixture with short episodes (cessation p ~ 0.35 per
#      month, monthly mean on the log_offset path, latent and detected feeds),
#      per-month Y / active / mu / cases in v1 and Y / active / mso / mu in C,
#      captured BEFORE the reorder existed. Under the default `full` both
#      kernels must reproduce it byte for byte.
#   C  C kernel under `half`: the cessation draw has moved before the monthly
#      mean but draws nothing new, so Yout / active / mso are identical to
#      `full` and mu differs by exactly 1/2 on the ceased unit-months only.
#      The emitted C under `full` carries no trace of the key.
#   V  v1 kernel under `half`: on every active unit-month mu / (20 x days) is
#      1/2 at onset, 1/2 on a ceased month, 1/4 when both, 1 otherwise;
#      under `full` the ceased factor is 1.
#   P  cache identity: `_cessexp-half` filter-cache token and the smoothing
#      meta.yml field only under `half`; a legacy cache refused by `half`.
#   R  replay exposure on the monthly fit panels (replayExposureRatios(), the
#      .ONSET_EXPOSURE_FRAC measurement): base panel, half moves the
#      case-weighted mu ratio (1.048 -> 1.018) and person-days towards 1;
#      gid_season panel, person-days towards 1 (1.081 -> 0.997), case-weighted
#      already ~1 (0.990 -> 0.988).
#
# Run: Rscript analysis/tests/test_cessation_exposure.R            (compare)
#      Rscript analysis/tests/test_cessation_exposure.R --capture  (write the golden)

suppressPackageStartupMessages({ library(tibble); library(dplyr); library(lubridate); library(tidyr) })
source("analysis/utils/smoothing_core.R")   # base, simulator_core, spatpomp_core{,_c}
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
CAPTURE <- "--capture" %in% commandArgs(trailingOnly = TRUE)
GOLDEN_PATH <- "analysis/tests/golden/cessation_exposure_golden.rds"
HAS_KEY <- exists("cessationExposureOf", mode = "function")   # FALSE only when capturing pre-change
cat("# test_cessation_exposure (mode:", if (CAPTURE) "capture" else "compare", ")\n\n")

# ===========================================================================
# Fixture: 4 units (2 countries x 2 ADM1), 12 months from 2024-01, frame in
# C-locale gid order (the pomp unit order), warm start [6,Inf].
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2) %>%
  mutate(gid = sprintf("%s_%d_1", country, adm1n), adm1_unit = sprintf("%s_%d", country, adm1n),
         AFRO_region = "REG_W", tu = "hmmW", min_seas = 1L, pop = 1e5) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>% as_tibble()
n_m    <- 12L
months <- as.Date("2024-01-01") %m+% months(0:(n_m - 1L))
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
  mutate(min_seas = 1L, country = sub("_[0-9]+_[0-9]+$", "", gid),
         t_start_seas = as.integer(yday(t_start) - 1L),
         mean_precip = 0, dev_precip = 0, mean_temp_daymean = 0, dev_temp_daymean = 0,
         outbreak = 0L, ongoing_outbreak_lag1 = 0L, ongoing_outbreak_country_lag1 = 0L,
         ongoing_outbreak_tu_lag1 = 0L, time_to_last_outbreak_cat = "[6,Inf]") %>% as_tibble()
horizon <- list(start = "2024-01-01", end = "2024-12-31")
window  <- list(start = "2024-01-01", end = "2025-01-01")
dim_days <- as.numeric(lubridate::days_in_month(months))

# Onset p ~ 0.3 / month; cessation lambda = exp(-4.3) x days ~ 0.41 -> p ~ 0.34.
d_occ <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -4.5),
                                     active_terms = "(Intercept)")
d_dur <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = -4.3),
                                     active_terms = "(Intercept)")
# Monthly mean 20 cases / day x exposure: the log_offset path (where the onset
# and cessation corrections live), no population offset, no mso basis.
d_mon <- make_draws_bundle_synthetic(frame_tu, S = 2L, fe_overrides = c("(Intercept)" = log(20)),
                                     active_terms = "(Intercept)")
d_mon$phi <- rep(5, 2L); d_mon$offset_pop <- FALSE; d_mon$offset_logoffset <- TRUE

R_REPS <- 6L
ce_arg <- function(ce) if (is.null(ce)) list() else list(cessation_exposure = ce)

run_v1 <- function(feed, ce = NULL, seed = 7L) {
  args <- c(list(tu = "hmmW", frame_tu = frame_tu, panel_tu = panel_tu,
                 draws_occ = d_occ, draws_dur = d_dur, draws_mon = d_mon,
                 replicate_ids = seq_len(R_REPS), seed = seed,
                 output_part = tempfile("cex_"), horizon = horizon, dry_run = TRUE,
                 covariate_feed = feed, delta_fn = 0.3, delta_fp = 0.1), ce_arg(ce))
  suppressWarnings(do.call(simulate_one_tu, args))
}
run_c <- function(feed, ce = NULL, seed = 7L) {
  args <- c(list(frame_tu, panel_tu, d_occ, window, fixed_draw_idx = 1L,
                 draws_dur = d_dur, draws_mon = d_mon,
                 measurement = "occurrence_duration_mon_NB",
                 covariate_feed = feed, delta_fn = 0.3, delta_fp = 0.1), ce_arg(ce))
  o <- do.call(build_spatpomp_tu_c, args)
  out <- lapply(seq_len(R_REPS), function(r) {
    set.seed(seed + r)
    sm <- spatPomp::simulate(o, format = "data.frame")
    sm <- sm[sm$time >= 1 & sm$time <= n_m, c("time", "unitname", "Yout", "active", "mso", "mu")]
    sm[order(sm$time, sm$unitname), ]
  })
  list(obj = o, sims = out)
}
# Emitted C with pomp's per-build header (file hash, time, salt) normalised.
src_of <- function(o) {
  x <- paste(vapply(o@solibs, `[[`, "", "src"), collapse = "\n")
  x <- gsub("pomp_[0-9a-f]{32}", "pomp_<hash>", x)
  gsub("/\\* (Time|Salt): [^*]*\\*/", "/* \\1 */", x)
}

# ===========================================================================
cat("## G: golden (default = full)\n")
cur <- list()
for (feed in c("latent", "detected")) {
  v <- run_v1(feed)
  cur[[paste0("v1_", feed)]] <- list(monthly_summary = v$monthly_summary, event_table = v$event_table)
  cur[[paste0("C_", feed)]]  <- run_c(feed)$sims
}
if (CAPTURE) {
  saveRDS(cur, GOLDEN_PATH)
  cat("  golden written:", GOLDEN_PATH, "\n")
  quit(status = 0L)
}
gold <- readRDS(GOLDEN_PATH)
for (nm in names(gold)) expect(identical(cur[[nm]], gold[[nm]]), sprintf("G1 %s identical to golden", nm))
n_ceased_gold <- sum(gold$v1_latent$event_table$censored == 0L, na.rm = TRUE)
cat(sprintf("  golden v1 latent: %d ceased episodes over %d replicates\n", n_ceased_gold, R_REPS))
expect(n_ceased_gold >= 10L, "G2 the fixture has short, ceasing episodes")
if (HAS_KEY) {
  for (feed in c("latent", "detected")) {
    v <- run_v1(feed, "full")
    expect(identical(list(monthly_summary = v$monthly_summary, event_table = v$event_table),
                     gold[[paste0("v1_", feed)]]) &&
           identical(run_c(feed, "full")$sims, gold[[paste0("C_", feed)]]),
           sprintf("G3 explicit full == golden (%s)", feed))
  }
}

# ===========================================================================
cat("\n## H: helpers\n")
expect(identical(cessationExposureOf(NULL), "full") && identical(cessationExposureOf(list()), "full") &&
       identical(cessationExposureOf(list(cessation_exposure = "half")), "half"),
       "H1 cessationExposureOf: default full, half accepted")
expect(all(vapply(list("Half", "", NA_character_, c("full", "half"), TRUE, 0.5),
                  function(v) !is.null(err_of(cessationExposureOf(list(cessation_exposure = v)))),
                  logical(1))), "H2 cessationExposureOf refuses anything else")
expect(identical(.cessation_log_offset_adj(c(TRUE, FALSE, NA)), c(log(0.5), 0, NA_real_)) &&
       identical(.cessation_log_offset_adj(c(TRUE, FALSE), use_logoffset = FALSE), 0),
       "H3 .cessation_log_offset_adj: log(1/2) on a ceased month, 0 otherwise / without log_offset")

# ===========================================================================
cat("\n## C: C kernel under half\n")
# "Ceased this month": active during month t (onset this month, or active
# entering it) and recorded inactive at t (the cessation reset is written back
# in the same step). Initial state: no unit active.
c_ceased <- function(s) {
  s <- s[order(s$unitname, s$time), ]
  prev <- ave(s$active, s$unitname, FUN = function(a) c(0, head(a, -1)))
  during <- s$Yout > 0.5 | prev > 0.5
  s$during <- during; s$ceased <- during & s$active < 0.5
  s[order(s$time, s$unitname), ]
}
for (feed in c("latent", "detected")) {
  cf <- run_c(feed, "full"); ch <- run_c(feed, "half")
  same_traj <- all(vapply(seq_len(R_REPS), function(r)
    identical(cf$sims[[r]][, c("time", "unitname", "Yout", "active", "mso")],
              ch$sims[[r]][, c("time", "unitname", "Yout", "active", "mso")]), logical(1)))
  expect(same_traj, sprintf("C1 %s: Yout / active / mso identical under full and half (no new RNG draw)", feed))
  ratio_ok <- TRUE; n_ceased <- 0L; n_active <- 0L
  for (r in seq_len(R_REPS)) {
    s  <- c_ceased(cf$sims[[r]]); h <- ch$sims[[r]]
    want <- ifelse(s$ceased, 0.5, 1) * s$mu
    ratio_ok <- ratio_ok && isTRUE(all.equal(h$mu, want, tolerance = 1e-14)) &&
                all((h$mu > 0) == s$during)
    n_ceased <- n_ceased + sum(s$ceased); n_active <- n_active + sum(s$during)
  }
  cat(sprintf("  %s: %d ceased of %d active unit-months\n", feed, n_ceased, n_active))
  expect(ratio_ok && n_ceased > 0L,
         sprintf("C2 %s: mu(half) = mu(full) / 2 on ceased unit-months, unchanged elsewhere, > 0 iff active", feed))
  if (feed == "latent") {
    expect(!grepl("av_act", src_of(cf$obj), fixed = TRUE) && grepl("av_act", src_of(ch$obj), fixed = TRUE) &&
           identical(src_of(cf$obj), src_of(run_c(feed)$obj)),
           "C3 emitted C: unchanged under full (== default), snapshot block only under half")
    expect(is.null(attr(cf$obj, "cessation_exposure")) && identical(attr(ch$obj, "cessation_exposure"), "half"),
           "C4 object attribute only under half")
  }
}
expect(!is.null(err_of(build_spatpomp_tu_c(frame_tu, panel_tu, d_occ, window, fixed_draw_idx = 1L,
                                           draws_dur = d_dur, draws_mon = d_mon,
                                           measurement = "occurrence_duration_mon_NB",
                                           cessation_exposure = "quarter"))),
       "C5 an unknown value is refused by the C builder")

# ===========================================================================
cat("\n## V: v1 kernel under half\n")
v1_factor <- function(v) {
  ms <- v$monthly_summary; ev <- v$event_table
  ms <- ms[ms$active == 1L, ]
  ms$mso <- NA_integer_; ms$ceased <- FALSE
  for (i in seq_len(nrow(ev))) {
    st <- match(ev$start_month[i], as.character(months))
    en <- if (is.na(ev$end_month[i])) n_m else match(ev$end_month[i], as.character(months))
    hit <- ms$replicate_id == ev$replicate_id[i] & ms$gid == ev$gid[i] &
           match(ms$month, as.character(months)) %in% st:en
    ms$mso[hit] <- match(ms$month[hit], as.character(months)) - st + 1L
    if (identical(ev$censored[i], 0L))
      ms$ceased[hit & ms$month == ev$end_month[i]] <- TRUE
  }
  ms$fac <- ms$mu_mon / (20 * dim_days[match(ms$month, as.character(months))])
  ms
}
for (feed in c("latent", "detected")) {
  vf <- v1_factor(run_v1(feed, "full")); vh <- v1_factor(run_v1(feed, "half"))
  want_f <- ifelse(vf$mso == 1L, 0.5, 1)
  want_h <- ifelse(vh$mso == 1L, 0.5, 1) * ifelse(vh$ceased, 0.5, 1)
  cat(sprintf("  %s: half run has %d ceased of %d active unit-months (%d single-month episodes)\n",
              feed, sum(vh$ceased), nrow(vh), sum(vh$ceased & vh$mso == 1L)))
  expect(!anyNA(vf$mso) && isTRUE(all.equal(vf$fac, want_f, tolerance = 1e-12)),
         sprintf("V1 %s full: exposure factor 1/2 at onset only (ceased months keep a full month)", feed))
  expect(!anyNA(vh$mso) && sum(vh$ceased) > 0L && isTRUE(all.equal(vh$fac, want_h, tolerance = 1e-12)),
         sprintf("V2 %s half: exposure factor 1/2 at onset x 1/2 on the ceased month", feed))
}
expect(!is.null(err_of(run_v1("latent", "quarter"))), "V3 an unknown value is refused by simulate_one_tu")

# ===========================================================================
cat("\n## P: cache identity\n")
ww <- list(start = "2023-01-01", end = "2024-01-01")
p0 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15)
p1 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, cessation_exposure = "full")
p2 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, cessation_exposure = "half")
p3 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, "season", "half")
expect(identical(p0, p1) && identical(p2, sub("\\.qs$", "_cessexp-half.qs", p0)) &&
       identical(p3, sub("\\.qs$", "_onset-season_cessexp-half.qs", p0)),
       "P1 filter cache path: `_cessexp-half` only under half")
expect(is.null(smoothing_cache_identity(list(NP = 10L, N_REPS = 2L),
                                        optional = smoothing_cache_optional_fields(list()))$cessation_exposure) &&
       identical(smoothing_cache_identity(list(NP = 10L, N_REPS = 2L),
                                          optional = smoothing_cache_optional_fields(
                                            list(CESSATION_EXPOSURE = "half")))$cessation_exposure, "half"),
       "P2 smoothing meta.yml: cessation_exposure written only when not full")
cd <- tempfile("cex_cache_"); dir.create(cd)
yaml::write_yaml(smoothing_cache_meta("transmission_unit"), smoothing_cache_meta_path(cd))
expect(!is.null(err_of(check_smoothing_cache_meta(cd, c(smoothing_cache_meta("transmission_unit"),
                                                        list(cessation_exposure = "half"))))) &&
       is.null(err_of(check_smoothing_cache_meta(cd, c(smoothing_cache_meta("transmission_unit"),
                                                       list(cessation_exposure = "full"))))),
       "P3 a legacy (full) cache is refused by a half run, accepted by a full one")

# ===========================================================================
cat("\n## R: replay exposure on the monthly fit panels\n")
# Base panel (episodes keyed as in the .ONSET_EXPOSURE_FRAC measurement): half
# must move both ratios towards 1. gid_season panel (episodes.key gid_season,
# the deployed fits): the onset fix alone already leaves the case-weighted ratio
# ~1 (cessation months carry few cases), so there half must move person-days
# towards 1 and keep the case-weighted ratio within 0.02 of 1.
for (pp in c(base = "generated_data/outbreaks_monthly_data_base-1990_2019.rds",
             gid_season = "generated_data/outbreaks_monthly_data_base-1990_2019_ep-gs.rds")) {
  if (!file.exists(pp)) { cat("  (skipped: ", pp, " not present)\n", sep = ""); next }
  rr <- replayExposureRatios(readRDS(pp))
  cat("  ", pp, "\n", sep = ""); print(as.data.frame(rr), digits = 4)
  cw <- setNames(rr$case_weighted, rr$mode); pd <- setNames(rr$person_days, rr$mode)
  if (grepl("_ep-gs", pp)) {
    expect(abs(pd[["half"]] - 1) < abs(pd[["full"]] - 1) && abs(cw[["half"]] - 1) < 0.02,
           "R2 gid_season panel: half moves person-days towards 1, case-weighted stays within 0.02 of 1")
  } else {
    expect(abs(cw[["half"]] - 1) < abs(cw[["full"]] - 1) && abs(pd[["half"]] - 1) < abs(pd[["full"]] - 1),
           "R1 base panel: half moves the case-weighted mu ratio and person-days towards 1")
  }
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
