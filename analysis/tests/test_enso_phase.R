# test_enso_phase.R — one ENSO phase function.
#
# HARD gate, quick, no data. Frozen reference copies of five phase code
# paths (.classifyEnsoEpisodes body, classifyEnsoAnalogSets' label block,
# ensoEventLaunch's in-phase test, report_smoothing's and report_dataset's
# case_when) must equal ensoPhase() bit for bit on random series with exact
# +/-0.5 values, NA and runs of 4, 5 and 6; only the episode rule may be
# labelled CPC / official; no ONI phase literal may reappear outside the
# constants; report_smoothing's month expansion of the oni runs and its rule
# resolution.
#
# Run: Rscript analysis/tests/test_enso_phase.R

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

# --- frozen reference code -------------------------------------------------
old_episodes <- function(index, threshold = 0.5, min_run = 5L) {
  n <- length(index)
  phase <- rep("neutral", n)
  phase[is.na(index)] <- NA_character_
  mark_runs <- function(flag, label) {
    flag[is.na(flag)] <- FALSE
    rl     <- rle(flag)
    ends   <- cumsum(rl$lengths)
    starts <- ends - rl$lengths + 1L
    hits   <- which(rl$values & rl$lengths >= min_run)
    for (i in hits) phase[starts[i]:ends[i]] <<- label
  }
  mark_runs(index >= threshold, "elnino")
  mark_runs(index <= -threshold, "nina")
  phase
}
old_sets_labels <- function(stat, thresholds = c(nina = -0.5, elnino = 0.5)) {
  labels <- rep(NA_character_, length(stat))
  ok     <- !is.na(stat)
  labels[ok & stat <= thresholds[["nina"]]]   <- "nina"
  labels[ok & stat >= thresholds[["elnino"]]] <- "elnino"
  labels[ok & is.na(labels)]                  <- "neutral"
  labels
}
old_in_phase <- function(x, phase, threshold = 0.5)
  if (identical(phase, "elnino")) x >= threshold else x <= -threshold
old_smoothing <- function(x) dplyr::case_when(x >=  0.5 ~ "elnino", x <= -0.5 ~ "nina", TRUE ~ "neutral")
old_dataset <- function(x) dplyr::case_when(x >=  0.5 ~ "El Niño", x <= -0.5 ~ "La Niña", !is.na(x) ~ "Neutral")

# --- series ------------------------------------------------------------------
set.seed(16)
series <- list(
  edges = c(0.5, -0.5, 0.4999, -0.4999, 0.5001, -0.5001, 0, NA),
  runs  = c(rep(0.6, 4), 0, rep(0.5, 5), 0, rep(-0.5, 6), NA, rep(-0.7, 4), rep(0.9, 5)),
  rand  = round(rnorm(400, 0, 0.8), 1),
  rand_na = { x <- round(rnorm(300, 0, 0.9), 2); x[sample(300, 25)] <- NA; x },
  ar    = as.numeric(stats::filter(rnorm(500, 0, 0.3), 0.9, method = "recursive")))

cat("# [1] bit-identical to the frozen code paths\n")
for (nm in names(series)) {
  x <- series[[nm]]
  expect(identical(ensoPhase(x, "episode"), old_episodes(x)), sprintf("%s: episode == old .classifyEnsoEpisodes", nm))
  expect(identical(.classifyEnsoEpisodes(x), old_episodes(x)), sprintf("%s: .classifyEnsoEpisodes wrapper", nm))
  expect(identical(.classifyEnsoEpisodes(x, threshold = 1, min_run = 3L), old_episodes(x, 1, 3L)),
         sprintf("%s: wrapper with non-default threshold / min_run", nm))
  expect(identical(ensoPhase(x, "threshold"), old_sets_labels(x)), sprintf("%s: threshold == old analog-set labels", nm))
  for (ph in c("elnino", "lanina")) {
    code <- if (ph == "elnino") "elnino" else "nina"
    expect(identical(ensoPhase(x, "threshold") == code, old_in_phase(x, ph)),
           sprintf("%s: ensoEventLaunch in-phase (%s)", nm, ph))
  }
  expect(identical(tidyr::replace_na(ensoPhase(x, "threshold"), "neutral"), old_smoothing(x)),
         sprintf("%s: report_smoothing case_when", nm))
  expect(identical(unname(ENSO_PHASE_LABELS[ensoPhase(x, "threshold")]), old_dataset(x)),
         sprintf("%s: report_dataset case_when", nm))
}
# classifyEnsoAnalogSets end to end (labels via ensoPhase) vs the frozen block.
traj <- matrix(round(rnorm(60 * 8, 0, 0.8), 2), nrow = 60)
traj[1:3, 4:6] <- c(0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 0.49, 0.5, 0.51)
cs <- classifyEnsoAnalogSets(traj, "2026-09-01")
expect(identical(cs$labels, old_sets_labels(cs$stat)), "classifyEnsoAnalogSets labels unchanged")
expect(identical(cs$definition$thresholds, list(nina = -0.5, elnino = 0.5)),
       "classifyEnsoAnalogSets definition (cache meta) unchanged")
# ensoEventLaunch on a synthetic record.
oni <- data.frame(date = seq(as.Date("2014-01-01"), by = "1 month", length.out = 24),
                  enso_index = c(0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.5, 0.6, 0.9, 1.2, 1.5, 2.0,
                                 2.4, 2.2, 1.5, 0.9, 0.5, 0.2, -0.1, -0.4, -0.5, -0.8, -1.0, -0.9))
expect(identical(ensoEventLaunch(oni, "2015-01-01", "elnino"), as.Date("2014-06-01")),
       "ensoEventLaunch: launch before the crossing (inclusive 0.5)")
expect(identical(ensoEventLaunch(oni, "2015-11-01", "lanina"), as.Date("2015-08-01")),
       "ensoEventLaunch: La Niña crossing at -0.5 inclusive")

cat("# [2] constants and definitions\n")
expect(identical(ENSO_PHASE_THRESHOLDS, c(nina = -0.5, elnino = 0.5)), "ENSO_PHASE_THRESHOLDS = +/-0.5")
expect(identical(ENSO_EPISODE_MIN_RUN, 5L), "ENSO_EPISODE_MIN_RUN = 5")
expect(setequal(unique(na.omit(ensoPhase(series$rand))), ENSO_PHASE_CODES), "codes are ENSO_PHASE_CODES")
for (st in names(ENSO_PHASE_STATISTICS)) {
  ep <- ensoPhaseDefinition("episode", st)$label
  th <- ensoPhaseDefinition("threshold", st)$label
  expect(grepl("CPC", ep), sprintf("episode definition (%s) names CPC", st))
  expect(!grepl("CPC|official", th, ignore.case = TRUE),
         sprintf("threshold definition (%s) never says CPC / official", st))
}
expect(inherits(try(ensoPhase(1, thresholds = c(nina = 1, elnino = -1)), silent = TRUE), "try-error"),
       "nina >= elnino refused")

cat("# [3] no ONI phase literal outside the constants\n")
files <- c("analysis/utils/enso.R", "analysis/utils/results_derive.R", "analysis/utils/results_core.R",
           "analysis/10_manuscript_figures.R", "analysis/utils/figures_manuscript_fig1.R",
           "analysis/utils/figures_manuscript_fig2.R", "analysis/utils/figures_manuscript_fig34.R",
           "analysis/11_manuscript_numbers.qmd")
pat <- ">= *0\\.5\\b|<= *-0\\.5\\b|> *0\\.5\\b|< *-0\\.5\\b"
for (f in files) {
  hits <- grep(pat, readLines(f, warn = FALSE), value = TRUE)
  hits <- hits[!grepl("^\\s*#", hits)]
  expect(!length(hits), sprintf("%s: no +/-0.5 comparison literal%s", f,
                                if (length(hits)) paste0(" (", hits[1], ")") else ""))
}

cat("# [4] report_smoothing phases\n")
env <- new.env()
for (f in c("analysis/utils/results_derive.R")) sys.source(f, envir = env, keep.source = FALSE)
ph <- tibble::tibble(phase = c("neutral", "elnino", "neutral"),
                     TL = as.Date(c("2023-01-01", "2023-06-01", "2024-05-01")) - 15,
                     TR = as.Date(c("2023-05-01", "2024-04-01", "2024-08-01")) + 15)
mp <- env$smoothing_phases_from_runs(ph)
expect(nrow(mp) == 20L && !anyDuplicated(paste(mp$year, mp$month)), "runs expand to one row per month")
expect(identical(mp$phase[mp$year == 2023L & mp$month == 6L], "elnino") &&
         identical(mp$phase[mp$year == 2024L & mp$month == 4L], "elnino") &&
         identical(mp$phase[mp$year == 2024L & mp$month == 5L], "neutral"),
       "episode runs map to their months (ends inclusive)")
expect(identical(env$smoothing_phase_rule(list(ENSO_PHASE_SMOOTHING = "threshold"))$rule, "threshold") &&
         identical(env$smoothing_phase_rule(list(ENSO_PHASE_SMOOTHING = "episode"))$rule, "episode"),
       "smoothing_phase_rule follows results.derivation.enso_phase.smoothing")

cat(sprintf("\n[test_enso_phase] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
