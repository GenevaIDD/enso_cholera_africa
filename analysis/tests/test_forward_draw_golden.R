# test_forward_draw_golden.R -- golden characterization of the forward
# simulation's posterior-draw indices and handoff map (golden first; HARD,
# quick).
#
# Captured on a small synthetic fixture (13 smoothing replicates, 2 of them failed; 40 forward
# replicates in two chunks; S = 7 / 5 / 3 draws for occ / dur / mon):
#   * .forecast_param_indices() under random / cycle / inherited and the legacy
#     spread_param_draws fallback (blank mode);
#   * select_handoff_warmup()'s idx under uniform / shuffled / loglik (seeded);
#   * handoff_map_table() for each weighting;
#   * smoothing_draw_idx_for_tu() (the per-unit random smoothing draws).
# Under the defaults these must stay byte-identical.
#
#   Rscript analysis/tests/test_forward_draw_golden.R --capture   # once, pre-change
#   Rscript analysis/tests/test_forward_draw_golden.R             # compare

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/spatpomp_core.R")
source("analysis/utils/smoothing_core.R")

CAPTURE <- "--capture" %in% commandArgs(trailingOnly = TRUE)
GOLDEN_PATH <- "analysis/tests/golden/forward_draw_golden.rds"
cat("# forward draw / handoff golden (mode:", if (CAPTURE) "capture" else "compare", ")\n")

# --- fixture -----------------------------------------------------------------
M_SM  <- 13L                                  # smoothing replicates
recs  <- lapply(seq_len(M_SM), function(i) list(rep_id = i, draw_d = c(5L, 2L, 7L, 7L, 1L, 3L, 6L,
                                                                       2L, 4L, 4L, 1L, 5L, 3L)[i]))
ok    <- setdiff(seq_len(M_SM), c(4L, 9L))    # two failed extractions
wl    <- lapply(ok, function(i) data.frame(rep = i, posterior_draw_idx = 1L))
ll    <- c(-120, -101, -99, -5000, -98, -100, -97, -300, -110, -96, -103, -102, -104)[ok]
SEED  <- 20261002L
CHUNKS <- list(c(n = 25L, off = 0L), c(n = 15L, off = 25L))

out <- list()
for (w in c("uniform", "shuffled", "loglik")) {
  for (ch in CHUNKS) {
    key <- sprintf("%s_off%d", w, ch[["off"]])
    set.seed(11L)
    sel <- select_handoff_warmup(wl, ll, n = ch[["n"]], offset = ch[["off"]], weighting = w, seed = SEED)
    out[[paste0("idx_", key)]] <- attr(sel, "idx")
    out[[paste0("map_", key)]] <- handoff_map_table("TUX", recs[ok], sel, ch[["off"]], as.Date("2026-09-01"), w)
    for (mode in list("random", "cycle", "inherited", "")) {
      for (spread in c(FALSE, TRUE)) {
        if (nzchar(mode) && spread) next
        out[[sprintf("pidx_%s_%s_%s", key, if (nzchar(mode)) mode else "blank", spread)]] <-
          .forecast_param_indices(sel, seq_len(ch[["n"]]) + ch[["off"]], S_occ = 7L, S_dur = 5L, S_mon = 3L,
                                  spread_param_draws = spread, param_draw_mode = mode, seed = SEED)
      }
    }
  }
}
S_sm <- list(SEED_ROOT = 20260624L, N_REPS = 9L)
tus  <- c("TU_B", "TU_A", "TU_C")
for (tu in tus) out[[paste0("smdraw_", tu)]] <- smoothing_draw_idx_for_tu(S_sm, tu, tus, 11L)

if (CAPTURE) {
  saveRDS(out, GOLDEN_PATH)
  cat(sprintf("captured %d entries -> %s\n", length(out), GOLDEN_PATH))
  quit(status = 0)
}

expect(file.exists(GOLDEN_PATH), paste("golden file exists:", GOLDEN_PATH))
gold <- readRDS(GOLDEN_PATH)
expect(identical(sort(names(gold)), sort(names(out))), "same captured entries as the golden")
for (k in intersect(names(gold), names(out)))
  expect(identical(serialize(gold[[k]], NULL), serialize(out[[k]], NULL)), paste("byte-identical:", k))

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
