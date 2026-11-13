# test_handoff_weighting.R -- forecast handoff state selection among the
# smoothing replicates (simulator.forecast.handoff_weighting; HARD, quick).
#
#   [1] uniform (default): balanced, deterministic assignment by global
#       replicate id; every smoothing replicate used n / m times; consumes no RNG.
#   [2] chunk safety: chunks (offset, n) concatenate to the unchunked selection.
#   [3] loglik: identical to resample_by_loglik() under the same seed.
#   [3b] attr(, "idx"): the position each forecast replicate took (07b's
#        handoff_map.csv), under both weightings and with an offset.
#   [4] handoff_loglik_ess(): 1 for a dominant replicate, m for equal logLiks.
#   [5] defaults: config_defaults.yml and get_default_options() say "uniform".
#   [6] shuffled: a seeded permutation of 1..m per block of m global ids --
#       balance (once per full block, +-1 over a partial block), chunk
#       invariance, determinism given the seed, no global RNG use, the seed
#       required; S = 200 cycle draws x m = 100 states: each draw meets many
#       states (uniform: one); idx -> handoff_map.csv.
#   [7] uniform unchanged: idx byte-identical to the pre-shuffled formula.
#
# Synthetic inputs only. Run: Rscript analysis/tests/test_handoff_weighting.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/spatpomp_core.R")

m  <- 7L
wl <- lapply(seq_len(m), function(i) data.frame(rep = i))
ll <- c(-100, -90, -5000, -95, -91, -300, -92)
ids <- function(x) vapply(x, function(d) d$rep, integer(1))

cat("\n[1] uniform\n")
set.seed(1); before <- .Random.seed
u <- select_handoff_warmup(wl, ll, n = 21L)
expect(identical(.Random.seed, before), "uniform consumes no RNG")
expect(identical(ids(u), rep(seq_len(m), 3L)), "balanced: replicate r -> ((r - 1) mod m) + 1")
expect(all(table(ids(u)) == 3L), "each smoothing replicate used n / m times")
expect(identical(ids(select_handoff_warmup(wl, ll, n = 5L)), 1:5), "n < m: the first n replicates")

cat("\n[2] chunk safety\n")
full <- ids(select_handoff_warmup(wl, ll, n = 30L))
chunks <- unlist(lapply(c(0L, 10L, 20L), function(o) ids(select_handoff_warmup(wl, ll, n = 10L, offset = o))))
expect(identical(chunks, full), "three chunks of 10 == one run of 30")

cat("\n[3] loglik\n")
set.seed(7); a <- select_handoff_warmup(wl, ll, n = 50L, weighting = "loglik")
set.seed(7); b <- resample_by_loglik(wl, ll, n = 50L)
a0 <- a; attr(a0, "idx") <- NULL
expect(identical(a0, b), "loglik == resample_by_loglik() under the same seed")
expect(all(ids(a) %in% c(2L, 5L, 7L, 4L, 1L)), "loglik never picks the -5000 / -300 replicates")
expect(!is.null(tryCatch(select_handoff_warmup(wl, ll, 3L, weighting = "nope"),
                         error = function(e) conditionMessage(e))), "unknown weighting refused")

cat("\n[3b] idx attribute\n")
expect(identical(attr(a, "idx"), ids(a)), "loglik: idx names the replicate each forecast replicate took")
uo <- select_handoff_warmup(wl, ll, n = 10L, offset = 12L)
expect(identical(attr(uo, "idx"), ((12L + seq_len(10L) - 1L) %% m) + 1L) &&
         identical(attr(uo, "idx"), ids(uo)),
       "uniform + offset: idx = ((offset + r - 1) mod m) + 1, matching the panels")

cat("\n[4] ESS\n")
expect(abs(handoff_loglik_ess(c(0, -200, -500)) - 1) < 1e-9, "dominant replicate: ESS 1")
expect(abs(handoff_loglik_ess(rep(-3, 4)) - 4) < 1e-9, "equal logLiks: ESS m")

cat("\n[5] defaults\n")
expect(identical(get_default_options()$simulator$forecast$handoff_weighting, "uniform"),
       "get_default_options(): uniform")
expect(identical(yaml::read_yaml("analysis/config_defaults.yml")$simulator$forecast$handoff_weighting, "uniform"),
       "config_defaults.yml: uniform")

cat("\n[6] shuffled\n")
M <- 100L
wm <- lapply(seq_len(M), function(i) data.frame(rep = i))
ll_m <- rep(-1, M)
sh <- function(n, offset = 0L, seed = 20261002L)
  attr(select_handoff_warmup(wm, ll_m, n = n, offset = offset, weighting = "shuffled", seed = seed), "idx")
set.seed(3); before <- .Random.seed
full <- sh(250L)
expect(identical(.Random.seed, before), "shuffled leaves the global RNG stream untouched")
rm(.Random.seed, envir = .GlobalEnv); invisible(sh(5L))
expect(!exists(".Random.seed", envir = .GlobalEnv), "shuffled creates no .Random.seed when none existed")
expect(all(vapply(0:1, function(j) setequal(full[j * M + seq_len(M)], seq_len(M)), logical(1))),
       "balance: each full block of m ids uses every smoothing replicate exactly once")
expect(!anyDuplicated(full[201:250]), "partial block: no state used twice")
expect(all(tabulate(full, M) %in% 2:3), "balance over 2.5 blocks: counts differ by at most one")
expect(!identical(full[1:100], full[101:200]), "blocks draw different permutations")
expect(!identical(full[1:100], seq_len(M)), "a block is not the uniform order")
sizes <- c(37L, 143L, 21L, 49L)
offs  <- c(0L, cumsum(sizes)[-length(sizes)])
chunks <- unlist(Map(sh, sizes, offs))
expect(identical(chunks, full), "chunk invariance: uneven slices across block edges concatenate to the full run")
expect(identical(sh(1L, offset = 149L), full[150L]), "a single replicate mid-block matches the full run")
expect(identical(sh(250L), full), "deterministic given the seed")
expect(!identical(sh(250L, seed = 20261003L), full), "another seed, another pairing")
expect(identical(ids(select_handoff_warmup(wm, ll_m, n = 250L, weighting = "shuffled", seed = 20261002L)), full),
       "the panels follow idx")
expect(inherits(try(select_handoff_warmup(wm, ll_m, n = 3L, weighting = "shuffled"), silent = TRUE), "try-error"),
       "shuffled without a seed is refused")
S <- 200L; R <- S * 50L
draw <- ((seq_len(R) - 1L) %% S) + 1L
k_sh <- tapply(sh(R), draw, function(x) length(unique(x)))
k_un <- tapply(attr(select_handoff_warmup(wm, ll_m, n = R), "idx"), draw, function(x) length(unique(x)))
expect(all(k_un == 1L), "uniform, S = 200 vs m = 100: each draw meets one state (the confounding)")
expect(min(k_sh) >= 25L && abs(mean(k_sh) - M * (1 - (1 - 1 / M)^50)) < 3,
       sprintf("shuffled, S = 200 vs m = 100: each draw meets %d-%d distinct states (mean %.1f, ~39.5 expected)",
               min(k_sh), max(k_sh), mean(k_sh)))
d_per_state <- tapply(draw, sh(R), function(x) length(unique(x)))
expect(all(tabulate(sh(R), M) == 100L) && min(d_per_state) >= 50L,
       "shuffled: every state used 100 times, by many draws (no state tied to a few draws)")
recs <- lapply(seq_len(M), function(i) list(rep_id = 1000L + i))
hm <- handoff_map_table("tu1", recs, select_handoff_warmup(wm, ll_m, n = 30L, offset = 90L,
                                                           weighting = "shuffled", seed = 20261002L),
                        90L, "2026-09-01", "shuffled")
expect(identical(hm$replicate_id, 91:120) && identical(hm$smoothing_rep_id, 1000L + full[91:120]) &&
         all(hm$handoff_weighting == "shuffled"),
       "handoff_map_table: smoothing_rep_id = the shuffled state's rep_id, by global replicate id")

cat("\n[7] uniform unchanged\n")
old_idx <- function(n, offset, m) ((as.integer(offset) + seq_len(n) - 1L) %% m) + 1L
same <- vapply(list(c(21L, 0L, 7L), c(10L, 12L, 7L), c(250L, 37L, 100L), c(5L, 0L, 100L)), function(a) {
  w <- lapply(seq_len(a[3]), function(i) data.frame(rep = i))
  u <- select_handoff_warmup(w, rep(0, a[3]), n = a[1], offset = a[2])
  ref <- w[old_idx(a[1], a[2], a[3])]; attr(ref, "idx") <- old_idx(a[1], a[2], a[3])
  identical(serialize(u, NULL), serialize(ref, NULL))
}, logical(1))
expect(all(same), "uniform: panels and idx byte-identical to ((offset + r - 1) mod m) + 1")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
