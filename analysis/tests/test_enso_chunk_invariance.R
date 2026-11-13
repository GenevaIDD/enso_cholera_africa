# test_enso_chunk_invariance.R -- replicate chunking of the ENSO scenarios
# (07b --replicate_offset / --n_replicates_total; HARD gate).
#
# slurm/submit_forecast_groups.sh --chunk_size runs a TU's
# replicates as concurrent chunks (offset + 1 .. offset + n). The ENSO samplers
# draw columns 1 .. n_rep from a fixed seed, so a chunk that sampled only its
# own n gave replicate offset + j the trajectory of replicate j: every chunk
# repeated the first chunk's trajectories. The loaders now sample the TU's
# total and keep the chunk's columns (utils/simulations.R::.enso_chunk_cols() /
# .enso_chunk_slice()).
#
#   [1] .enso_chunk_cols(): NULL (no slicing) when unchunked or the total is
#       unknown; the chunk's columns otherwise; refuses a chunk past the total.
#   [2] with a real sampler (synthetic CFSv2 inputs), the chunks' slices bound
#       together equal the unchunked sample for every lag matrix; sampling only
#       the chunk's n does not (the bug this fixes).
#   [3] both loaders route every ENSO sampler's n_rep through n_enso and slice
#       the scenarios they return; 07b reads --n_replicates_total and the
#       payload passes it.
#
# Synthetic inputs only. Run: Rscript analysis/tests/test_enso_chunk_invariance.R

suppressPackageStartupMessages({ library(lubridate); library(dplyr) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/enso.R")
source("analysis/utils/simulations.R")
err_of <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

cat("\n[1] .enso_chunk_cols()\n")
expect(is.null(.enso_chunk_cols(list(), 10L)), "no offset / total -> NULL")
expect(is.null(.enso_chunk_cols(list(replicate_offset = 20L), 10L)), "total unknown -> NULL (earlier behaviour)")
expect(is.null(.enso_chunk_cols(list(replicate_offset = 0L, replicate_total = 10L), 10L)),
       "single chunk (offset 0, total = n) -> NULL")
ch <- .enso_chunk_cols(list(replicate_offset = 20L, replicate_total = 50L), 10L)
expect(identical(ch$n_build, 50L) && identical(ch$cols, 21:30), "chunk 21..30 of 50")
expect(!is.null(err_of(.enso_chunk_cols(list(replicate_offset = 45L, replicate_total = 50L), 10L))),
       "chunk past the total refused")
expect(is.null(.enso_chunk_slice(NULL, ch)), "slice of NULL scenarios is NULL")

cat("\n[2] chunks reproduce the unchunked sample (synthetic CFSv2)\n")
hist_d <- seq(as.Date("1985-01-01"), as.Date("2026-06-01"), by = "month")
nino   <- 1.2 * sin(seq_along(hist_d) / 7)
hist_n <- tibble::tibble(year = year(hist_d), month = month(hist_d), enso_index = nino)
oni    <- as.numeric(stats::filter(nino, rep(1/3, 3), sides = 2))
hist_o <- tibble::tibble(year = year(hist_d), month = month(hist_d), enso_index = oni) %>%
  dplyr::filter(!is.na(enso_index), hist_d[seq_along(oni)] <= as.Date("2026-05-01"))
cfs_d  <- seq(as.Date("2026-05-01"), as.Date("2026-12-01"), by = "month")
cfs    <- tidyr::expand_grid(target_date = cfs_d, member = 1:5) %>%
  dplyr::mutate(nino34_anomaly = 1.5 + 0.1 * member + 0.05 * as.numeric(factor(target_date)))
sim_m  <- seq(as.Date("2026-07-01"), as.Date("2027-06-01"), by = "month")
smp <- function(n) suppressMessages(sampleEnsoCFSv2Trajectories(
  cfs, n_rep = n, historical_oni = hist_o, historical_nino34 = hist_n, sim_months = sim_m, max_lag = 12L,
  rebase = "none", seed = 11L))
N <- 12L; n <- 4L
full <- smp(N)
chunks <- lapply(seq(0L, N - n, by = n), function(off) {
  cc <- .enso_chunk_cols(list(replicate_offset = off, replicate_total = N), n)
  .enso_chunk_slice(smp(cc$n_build %||% n), cc)
})
bound <- lapply(stats::setNames(names(full), names(full)), function(k) do.call(cbind, lapply(chunks, `[[`, k)))
expect(identical(unname(lapply(bound, unname)), unname(lapply(full, unname))),
       sprintf("%d chunks of %d == unchunked sample of %d, every lag matrix", N / n, n, N))
expect(identical(dim(chunks[[2]][[1]]), c(nrow(full[[1]]), n)), "a chunk's matrices have n columns")
own <- smp(n)   # the old behaviour for the second chunk: its own n from the same seed
expect(!identical(unname(own[[1]]), unname(full[[1]][, n + seq_len(n)])),
       "sampling only the chunk's n differs from the unchunked columns (the fixed bug)")

cat("\n[3] loaders, 07b and the payload\n")
src <- readLines("analysis/utils/simulations.R")
body_of <- function(fn) {
  i <- grep(sprintf("^%s <- function", fn), src)
  j <- c(grep("^[.A-Za-z_][A-Za-z0-9_.]* <- function", src), length(src) + 1L)
  paste(src[i:(min(j[j > i]) - 1L)], collapse = "\n")
}
for (fn in c("load_artefacts", ".load_artefacts_grouped")) {
  b <- body_of(fn)
  # Each sampler call: the first n_rep argument after its start.
  starts <- gregexpr("enso_scenarios <- sampleEnso[A-Za-z0-9]+\\(", b)[[1]]
  n_rep_of <- vapply(starts, function(i) regmatches(substring(b, i), regexpr("n_rep\\s*=\\s*[A-Za-z_.]+",
                                                                            substring(b, i))), "")
  expect(length(starts) == 4L, sprintf("%s: four ENSO sampler calls found (%d)", fn, length(starts)))
  expect(all(grepl("n_rep\\s*=\\s*n_enso$", n_rep_of)), sprintf("%s: every ENSO sampler gets n_rep = n_enso", fn))
  expect(grepl("enso_chunk\\s*<-\\s*\\.enso_chunk_cols\\(opt, n_replicates\\)", b),
         sprintf("%s: chunk columns resolved from opt", fn))
  expect(grepl("enso_scenarios\\s*=\\s*\\.enso_chunk_slice\\(enso_scenarios, enso_chunk\\)", b),
         sprintf("%s: returned scenarios sliced to the chunk", fn))
}
b07 <- paste(readLines("analysis/07b_simulator_spatpomp.R"), collapse = "\n")
expect(grepl('"--n_replicates_total"', b07, fixed = TRUE) && grepl("opt$replicate_total", b07, fixed = TRUE),
       "07b reads --n_replicates_total into opt$replicate_total")
PAYLOAD <- "analysis/bash/slurm/submit_forecast_groups.sh"
if (!file.exists("analysis/RELEASE") || file.exists(PAYLOAD)) {
  pay <- paste(readLines(PAYLOAD), collapse = "\n")
  expect(grepl('--n_replicates_total "${N_REPLICATES}"', pay, fixed = TRUE),
         "the forecast payload passes the TU's total to every chunk")
}

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
