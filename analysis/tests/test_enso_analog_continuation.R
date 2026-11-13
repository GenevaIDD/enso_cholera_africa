# test_enso_analog_continuation.R -- k-NN analog continuation of ENSO past the
# CFSv2 span (HARD gate, quick, synthetic inputs only).
#
#   [1] loadEnsoContinuationLibrary() indexes a synthetic library by calendar
#       month; ensoAnalogNeighbours() puts the planted analog first, keeps the
#       calendar month and honours min_sep within a model, max_per_model and
#       exclude_models.
#   [2] continueEnsoAnalog(): path(j) = analog(i + j) + (seam - analog(i)) * rho^j
#       (rho = exp(-1 / blend_months)); blend_months = Inf is a constant offset;
#       seed-reproducible.
#   [3] sampleEnsoCFSv2Trajectories(): neutral and ar1 are identical with or
#       without the analog arguments, and identical to the sampler before the
#       analog branch (newest commit without continueEnsoAnalog()), every lag
#       and lead matrix.
#   [4] analog in the sampler: months up to the last span month identical to
#       neutral; continued months finite and seed-reproducible; lead matrices
#       past the span filled from the same path.
#   [5] refusals: missing library (names the build command), horizon too short,
#       legacy_seam with analog, non-finite seam.
#   [6] ensoContinuationProvenance(): method only for neutral / ar1; analog adds
#       the library path and digest.
#
# Run: Rscript analysis/tests/test_enso_analog_continuation.R

suppressPackageStartupMessages({ library(lubridate); library(dplyr) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/enso.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
if (!requireNamespace("qs2", quietly = TRUE)) { cat("[SKIP] qs2 not available\n"); quit(status = 0) }

# Previous sampler: the newest commit whose enso.R has no analog continuation.
old_env <- new.env(); have_old <- FALSE
for (rev in sprintf("HEAD~%d", 0:8)) {
  txt <- tryCatch(system2("git", c("show", paste0(rev, ":analysis/utils/enso.R")), stdout = TRUE, stderr = FALSE),
                  warning = function(w) NULL, error = function(e) NULL)
  if (length(txt) && !any(grepl("continueEnsoAnalog", txt, fixed = TRUE))) {
    tf <- tempfile(fileext = ".R"); writeLines(txt, tf); sys.source(tf, envir = old_env)
    have_old <- TRUE; break
  }
}

# --- synthetic library -----------------------------------------------------------
H <- 24L
set.seed(11)
mk_series <- function(n) as.numeric(stats::arima.sim(list(ar = 0.9), n = n, sd = 0.3))
xa <- mk_series(600); xb <- mk_series(480)
mon_a <- rep_len(1:12, 600); mon_b <- rep_len(1:12, 480)
# Planted analog in model B at index i0 (a December): a distinctive state.
i0 <- 300L; stopifnot(mon_b[i0] == 12L)
xb[i0 - 6L] <- 2.9; xb[i0 - 3L] <- 2.2; xb[i0] <- 2.6
lib_raw <- list(series = list(A = xa, B = xb),
                features = rbind(ensoContinuationFeatures(xa, mon_a, "A", H),
                                 ensoContinuationFeatures(xb, mon_b, "B", H)))
ok <- lib_raw$features$ok
lib_raw$feature_sd <- vapply(c("v", "tend", "phase"), function(f) sd(lib_raw$features[[f]][ok]), numeric(1))
lib_raw$meta <- list(horizon_months = H)
lib_raw$meta$digest <- ensoContinuationLibraryDigest(lib_raw)
lib_path <- tempfile(fileext = ".qs"); qs2::qs_save(lib_raw, lib_path)

cat("\n[1] library + neighbours\n")
lib <- loadEnsoContinuationLibrary(lib_path)
expect(identical(names(lib$by_month), as.character(1:12)) && all(lib$by_month[["12"]]$mon == 12L),
       "library indexed by calendar month")
st <- list(mon = 12L, v = 2.6, tend = 2.6 - 2.2, phase = 2.9)
nn <- ensoAnalogNeighbours(lib, st, k = 8L, min_sep = 24L)
expect(nn$model[1] == "B" && nn$i[1] == i0, "planted analog ranked first")
expect(nrow(nn) == 8L, "k neighbours returned")
expect(all(lib_raw$features$mon[match(paste(nn$model, nn$i), paste(lib_raw$features$model, lib_raw$features$i))] == 12L),
       "neighbours share the seam calendar month")
sep_ok <- all(vapply(split(nn$i, nn$model), function(v) length(v) < 2L || min(diff(sort(v))) >= 24L, logical(1)))
expect(sep_ok, "at most one analog per min_sep months within a model")
nn_cap <- ensoAnalogNeighbours(lib, st, k = 8L, min_sep = 24L, max_per_model = 2)
expect(all(table(nn_cap$model) <= 2L) && nn_cap$model[1] == "B" && nn_cap$i[1] == i0,
       "max_per_model caps each model; nearest still first")
nn_ex <- ensoAnalogNeighbours(lib, st, k = 8L, min_sep = 24L, exclude_models = "B")
expect(all(nn_ex$model == "A") && nrow(nn_ex) == 8L, "exclude_models drops a model")
expect(identical(ensoAnalogNeighbours(lib, st, k = 8L, min_sep = 24L, max_per_model = Inf,
                                      exclude_models = character(0)), nn),
       "max_per_model = Inf, no exclusion: unchanged")
expect(!is.null(err_of(ensoAnalogNeighbours(lib, st, exclude_models = "CM2.1-NMME"))),
       "excluding a model the library lacks is refused")
expect(!is.null(err_of(ensoAnalogNeighbours(lib, st, max_per_model = 0))), "max_per_model < 1 refused")

cat("\n[2] join\n")
# k = 1 makes the draw deterministic (the planted analog).
seam <- 2.4
p <- continueEnsoAnalog(lib, seam, st$tend, st$phase, 12L, n_months = 10L, k = 1L, blend_months = 3)
want <- xb[i0 + 1:10] + (seam - xb[i0]) * exp(-(1:10) / 3)
expect(isTRUE(all.equal(as.numeric(p), want)), "faded join: analog(i+j) + (seam - analog(i)) * rho^j")
pc <- continueEnsoAnalog(lib, seam, st$tend, st$phase, 12L, n_months = 10L, k = 1L, blend_months = Inf)
expect(isTRUE(all.equal(as.numeric(pc), xb[i0 + 1:10] + (seam - xb[i0]))), "blend_months = Inf: constant offset")
set.seed(5); r1 <- continueEnsoAnalog(lib, rep(seam, 6), rep(st$tend, 6), rep(st$phase, 6), 12L, 10L, k = 5L)
set.seed(5); r2 <- continueEnsoAnalog(lib, rep(seam, 6), rep(st$tend, 6), rep(st$phase, 6), 12L, 10L, k = 5L)
expect(identical(r1, r2) && identical(dim(r1), c(10L, 6L)), "seed-reproducible, [n_months x n_rep]")

# --- synthetic sampler inputs (as test_enso_ar1_continuation.R) ------------------
hist_d  <- seq(as.Date("1985-01-01"), as.Date("2026-06-01"), by = "month")
nino    <- 1.2 * sin(seq_along(hist_d) / 7)
hist_n  <- tibble::tibble(year = year(hist_d), month = month(hist_d), enso_index = nino)
oni     <- as.numeric(stats::filter(nino, rep(1/3, 3), sides = 2))
hist_o  <- tibble::tibble(year = year(hist_d), month = month(hist_d), enso_index = oni) %>%
  dplyr::filter(!is.na(enso_index), hist_d[seq_along(oni)] <= as.Date("2026-05-01"))
cfs_d   <- seq(as.Date("2026-05-01"), as.Date("2026-12-01"), by = "month")
cfs     <- tidyr::expand_grid(target_date = cfs_d, member = 1:3) %>%
  dplyr::mutate(nino34_anomaly = 1.5 + 0.1 * member + 0.05 * as.numeric(factor(target_date)))
sim_m   <- seq(as.Date("2026-07-01"), as.Date("2027-06-01"), by = "month")
run <- function(fun = sampleEnsoCFSv2Trajectories, ...)
  suppressMessages(fun(cfs, n_rep = 6L, historical_oni = hist_o, historical_nino34 = hist_n,
                       sim_months = sim_m, max_lag = 12L, rebase = "none", seed = 7L, ...))
analog_args <- list(analog_library = lib_path, analog_k = 3L, analog_blend_months = 3, analog_min_sep = 24L)

cat("\n[3] neutral / ar1 unchanged\n")
for (lead in c(0L, 3L)) for (m in c("neutral", "ar1")) {
  base <- run(beyond_cfsv2 = m, ar1_phi = 0.8, max_lead = lead)
  with_args <- do.call(run, c(list(beyond_cfsv2 = m, ar1_phi = 0.8, max_lead = lead), analog_args))
  expect(identical(base, with_args), sprintf("%s, max_lead %d: analog args ignored", m, lead))
  if (have_old) {
    old <- run(old_env$sampleEnsoCFSv2Trajectories, beyond_cfsv2 = m, ar1_phi = 0.8, max_lead = lead)
    expect(identical(base, old), sprintf("%s, max_lead %d: identical to the pre-analog sampler", m, lead))
  }
}
if (!have_old) cat("  SKIP  pre-analog enso.R not found in git history\n")

cat("\n[4] analog in the sampler\n")
neu <- run(max_lead = 3L)
ana <- do.call(run, c(list(beyond_cfsv2 = "analog", max_lead = 3L), analog_args))
last  <- which(sim_m == as.Date("2026-12-01"))
after <- which(sim_m > as.Date("2026-12-01"))
expect(identical(ana$enso_lag0[seq_len(last), ], neu$enso_lag0[seq_len(last), ]),
       "lag0 up to the last span month identical to neutral")
expect(all(is.finite(ana$enso_lag0[after, ])) && any(ana$enso_lag0[after, ] != 0),
       "continued months finite and not neutral")
expect(identical(ana$enso_lag3[last + 3L, ], ana$enso_lag0[last, ]), "lag3 at span end + 3 reads the last-span ONI")
expect(all(is.finite(ana$enso_lead3)), "leads past the span filled")
expect(identical(ana$enso_lead1[-length(sim_m), ], ana$enso_lag0[-1L, ]), "lead1(t) = lag0(t + 1) per replicate")
ana2 <- do.call(run, c(list(beyond_cfsv2 = "analog", max_lead = 3L), analog_args))
expect(identical(ana, ana2), "seed-reproducible")
lib_obj <- do.call(run, c(list(beyond_cfsv2 = "analog", max_lead = 3L),
                          modifyList(analog_args, list(analog_library = lib))))
expect(identical(ana, lib_obj), "library passed as path or as a loaded object: same output")
ex_a <- do.call(run, c(list(beyond_cfsv2 = "analog", max_lead = 3L, analog_exclude_models = "A"), analog_args))
ex_b <- do.call(run, c(list(beyond_cfsv2 = "analog", max_lead = 3L, analog_exclude_models = "B"), analog_args))
expect(identical(ex_a$enso_lag0[seq_len(last), ], neu$enso_lag0[seq_len(last), ]) &&
       !identical(ex_a$enso_lag0[after, ], ex_b$enso_lag0[after, ]),
       "sampler passes exclude_models through (only model A vs only model B differ past the span)")

cat("\n[5] refusals\n")
e <- err_of(loadEnsoContinuationLibrary(file.path(tempdir(), "no_such_library.qs")))
expect(!is.null(e) && grepl("03_prepare_enso_continuation_library.R", e, fixed = TRUE),
       "missing library names the build command")
e <- err_of(continueEnsoAnalog(lib, 1, 0, 0, 12L, n_months = H + 1L))
expect(!is.null(e) && grepl("horizon_months", e, fixed = TRUE), "continuation longer than the library horizon refused")
expect(!is.null(err_of(do.call(run, c(list(beyond_cfsv2 = "analog", legacy_seam = TRUE), analog_args)))),
       "legacy_seam with analog refused")
expect(!is.null(err_of(continueEnsoAnalog(lib, NA_real_, 0, 0, 12L, 5L))), "non-finite seam refused")

cat("\n[6] provenance\n")
expect(identical(ensoContinuationProvenance(list(method = "neutral", library = "/nonexistent")),
                 list(method = "neutral")), "neutral: method only, library not read")
expect(identical(ensoContinuationProvenance(NULL), list(method = "neutral")), "absent block = neutral")
pv <- ensoContinuationProvenance(list(method = "analog", library = lib_path, k = 3L))
expect(identical(pv$library_digest, lib_raw$meta$digest) && pv$library == lib_path && pv$k == 3L,
       "analog: library path + digest + knobs")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
