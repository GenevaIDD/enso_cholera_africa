# test_joint_draws.R — HARD: joint state-parameter draw design.
#
# Usage:
#   Rscript analysis/tests/test_joint_draws.R
#
# Details:
#   Quick, synthetic inputs only. Pins the smoothing draw design
#   (per_tu_random vs shared_cycle), its meta.yml identity and the cache
#   refusals, the joint forward draws (forecast replicate r starts from
#   smoothing replicate ((r - 1) mod S_s) + 1 and runs its draw), the 09
#   design check, the event generator's joint mode and the defaults.

suppressPackageStartupMessages({
  library(dplyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/spatpomp_core.R")
source("analysis/utils/smoothing_core.R")
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")
ev <- new.env()
sys.source("analysis/utils/make_event_attribution_configs.R", envir = ev)

err <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
warns <- function(expr) {
  w <- NULL
  withCallingHandlers(force(expr), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  w
}

# --- [1] smoothing draw design ------------------------------------------------
# [1] smoothing draw design: per_tu_random reproduces the per-unit seeded
#     sample; shared_cycle gives replicate j draw j in every unit, consumes no
#     RNG, and refuses n_replicates != max_draws; the key is validated; 07b's
#     own filter cache name carries a non-default design.
cat("\n[1] smoothing draw design\n")
tus <- c("TU_C", "TU_A", "TU_B")
S_old <- list(SEED_ROOT = 20260624L, N_REPS = 9L)
old_fn <- function(S, tu, tus, m) {
  set.seed(S$SEED_ROOT + match(tu, sort(tus)))
  sample.int(m, S$N_REPS, replace = TRUE)
}
expect(
  all(vapply(
    tus,
    function(tu) {
      identical(
        smoothing_draw_idx_for_tu(S_old, tu, tus, 11L),
        old_fn(S_old, tu, tus, 11L)
      )
    },
    logical(1)
  )),
  "per_tu_random (and an S without DRAW_DESIGN): the old seeded per-unit sample"
)
S_sc <- c(S_old, list(DRAW_DESIGN = "shared_cycle"))
set.seed(3)
before <- .Random.seed
sc <- lapply(tus, function(tu) smoothing_draw_idx_for_tu(S_sc, tu, tus, 9L))
expect(identical(.Random.seed, before), "shared_cycle consumes no RNG")
expect(
  all(vapply(sc, identical, logical(1), 1:9)),
  "shared_cycle: replicate j runs draw j in every unit"
)
expect(
  !is.null(err(smoothing_draw_idx_for_tu(S_sc, "TU_A", tus, 200L))),
  "shared_cycle refuses n_replicates != max_draws"
)
expect(
  !is.null(err(smoothing_draw_design_of(list(draw_design = "cycle")))),
  "an unknown draw_design is refused"
)
expect(
  identical(smoothing_draw_design_of(NULL), "per_tu_random"),
  "unset draw_design = per_tu_random"
)

ww <- list(start = "2009-01-01", end = "2011-01-01")
p0 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, "inactive")
expect(
  identical(
    .make_filter_cache_path(
      "TUA",
      "c",
      "occ",
      100L,
      4L,
      ww,
      0.15,
      "inactive",
      draw_design = "per_tu_random"
    ),
    p0
  ) &&
    grepl(
      "_draws-shared_cycle.qs$",
      .make_filter_cache_path(
        "TUA",
        "c",
        "occ",
        100L,
        4L,
        ww,
        0.15,
        "inactive",
        draw_design = "shared_cycle"
      )
    ),
  "07b filter cache name: unchanged under per_tu_random, `_draws-shared_cycle` token otherwise"
)

# --- [2] meta.yml identity ----------------------------------------------------
# [2] meta.yml identity: a default smoothing writes neither draw_design nor
#     max_draws; shared_cycle writes both; a joint expectation refuses a
#     per_tu_random cache or another max_draws and accepts the matching one;
#     a forecast that does not ask (no expectation) accepts any cache.
cat("\n[2] meta.yml identity\n")
S_meta <- list(
  PARTITION_BY = "transmission_unit",
  NP = 10L,
  N_REPS = 9L,
  run_tag = "t",
  MAX_DRAWS_OCC = 200L,
  DRAW_DESIGN = "per_tu_random"
)
m_def <- smoothing_cache_identity(
  S_meta,
  tus = tus,
  optional = smoothing_cache_optional_fields(S_meta)
)
expect(
  is.null(m_def$draw_design) && is.null(m_def$max_draws),
  "a default smoothing writes neither field"
)
S_meta2 <- modifyList(
  S_meta,
  list(MAX_DRAWS_OCC = 9L, DRAW_DESIGN = "shared_cycle")
)
m_sc <- smoothing_cache_identity(
  S_meta2,
  tus = tus,
  optional = smoothing_cache_optional_fields(S_meta2)
)
expect(
  identical(m_sc$draw_design, "shared_cycle") && identical(m_sc$max_draws, 9L),
  "shared_cycle writes draw_design and max_draws"
)
wd <- function(meta) {
  d <- tempfile("smc_")
  dir.create(d)
  yaml::write_yaml(meta, file.path(d, "meta.yml"))
  d
}
d_def <- wd(m_def)
d_sc <- wd(m_sc)
exp_base <- smoothing_cache_meta("transmission_unit")
exp_joint <- function(S) {
  c(exp_base, list(draw_design = "shared_cycle", max_draws = S))
}
expect(
  !is.null(err(check_smoothing_cache_meta(d_def, exp_joint(200L)))),
  "joint refuses a per_tu_random cache"
)
expect(
  !is.null(err(check_smoothing_cache_meta(d_sc, exp_joint(100L)))),
  "joint refuses a cache thinned to another max_draws"
)
expect(
  is.null(err(check_smoothing_cache_meta(d_sc, exp_joint(9L)))),
  "joint accepts the matching shared_cycle cache"
)
expect(
  is.null(err(check_smoothing_cache_meta(d_sc, exp_base))) &&
    is.null(err(check_smoothing_cache_meta(d_def, exp_base))),
  "a forecast that does not ask accepts either cache"
)

# --- [3] joint forward draws --------------------------------------------------
# [3] joint forward draws: in every unit and arm (factual / counterfactual)
#     and every chunk, forecast replicate r starts from smoothing replicate
#     ((r - 1) mod S_s) + 1 and runs the draw that replicate ran; the same r
#     runs the same draw in every unit; handoff_map.csv carries the draw;
#     a missing smoothing replicate and a draw outside the bundle are refused.
cat("\n[3] joint forward draws\n")
S_s <- 9L
R <- 30L
S_occ <- 9L
units <- c("TU_A", "TU_B", "TU_C")
# shared_cycle smoothing records: rep_id j ran draw j; TU_B's records arrive out
# of order and TU_C's replicate 4 failed its extraction in a second scenario.
recs_of <- function(u) {
  ord <- if (u == "TU_B") rev(seq_len(S_s)) else seq_len(S_s)
  draw <- smoothing_draw_idx_for_tu(S_sc, u, units, S_s)
  lapply(ord, function(j) list(tu = u, rep_id = j, draw_d = draw[j]))
}
forward <- function(recs, n, off, arm_seed) {
  wl <- lapply(recs, function(r) {
    data.frame(gid = r$tu, posterior_draw_idx = 1L)
  })
  rid <- vapply(recs, function(r) as.integer(r$rep_id), integer(1))
  sel <- select_handoff_warmup(
    wl,
    rep(0, length(wl)),
    n = n,
    offset = off,
    weighting = "joint",
    seed = arm_seed,
    rep_ids = rid,
    n_smoothing = S_s
  )
  hd <- handoff_draw_d(recs, sel)
  p <- .forecast_param_indices(
    sel,
    seq_len(n) + off,
    S_occ = S_occ,
    S_dur = 7L,
    S_mon = 9L,
    param_draw_mode = "joint",
    seed = arm_seed,
    handoff_draw_d = hd
  )
  hm <- handoff_map_table(
    recs[[1]]$tu,
    recs,
    sel,
    off,
    as.Date("2026-09-01"),
    "joint",
    draws = hd
  )
  list(occ = p$occ, dur = p$dur, mon = p$mon, map = hm)
}
res <- list()
for (arm in c("factual", "counterfactual")) {
  for (u in units) {
    recs <- recs_of(u)
    chunks <- list(forward(recs, 12L, 0L, 1L), forward(recs, 18L, 12L, 2L))
    res[[paste(arm, u)]] <- list(
      occ = unlist(lapply(chunks, `[[`, "occ")),
      dur = unlist(lapply(chunks, `[[`, "dur")),
      map = do.call(rbind, lapply(chunks, `[[`, "map"))
    )
  }
}
j_of <- ((seq_len(R) - 1L) %% S_s) + 1L
ok_pair <- vapply(
  res,
  function(x) {
    identical(x$map$smoothing_rep_id, j_of) &&
      identical(x$occ, x$map$smoothing_draw)
  },
  logical(1)
)
expect(
  all(ok_pair),
  "every unit, arm and chunk: replicate r starts from smoothing replicate ((r - 1) mod S_s) + 1 and runs its draw"
)
ok_rec <- vapply(
  names(res),
  function(k) {
    u <- sub("^\\S+ ", "", k)
    recs <- recs_of(u)
    d_by_rep <- setNames(
      vapply(recs, `[[`, integer(1), "draw_d"),
      vapply(recs, `[[`, integer(1), "rep_id")
    )
    identical(
      unname(d_by_rep[as.character(res[[k]]$map$smoothing_rep_id)]),
      res[[k]]$occ
    )
  },
  logical(1)
)
expect(
  all(ok_rec),
  "the forward draw is the smoothing record's draw_d of the handoff replicate"
)
expect(
  length(unique(lapply(res, `[[`, "occ"))) == 1L,
  "the same replicate runs the same draw in every unit and arm"
)
expect(
  identical(res[[1]]$occ, j_of),
  "shared_cycle: the joint draw of replicate r is ((r - 1) mod S_s) + 1"
)
expect(
  identical(res[[1]]$dur, ((res[[1]]$occ - 1L) %% 7L) + 1L),
  "dur / mon stay paired to the occurrence draw"
)
recs_c <- recs_of("TU_C")[-4]
expect(
  !is.null(err(forward(recs_c, 12L, 0L, 1L))),
  "a missing smoothing replicate is refused under joint"
)
expect(
  !is.null(err(.forecast_param_indices(
    list(1, 2),
    1:2,
    S_occ = 5L,
    param_draw_mode = "joint",
    handoff_draw_d = c(3L, 7L)
  ))),
  "a draw outside the bundle is refused"
)
expect(
  !is.null(err(.forecast_param_indices(
    list(1, 2),
    1:2,
    S_occ = 5L,
    param_draw_mode = "joint"
  ))),
  "joint without the handoff draws is refused"
)
hm_def <- handoff_map_table(
  "TUX",
  recs_of("TU_A"),
  select_handoff_warmup(as.list(1:9), rep(0, 9), 5L),
  0L,
  "2026-09-01",
  "uniform"
)
expect(
  identical(
    names(hm_def),
    c(
      "tu",
      "replicate_id",
      "smoothing_rep_id",
      "handoff_month",
      "handoff_weighting"
    )
  ),
  "handoff_map.csv without draws keeps its columns (non-joint runs)"
)

# --- [4] 09 design check ------------------------------------------------------
# [4] 09: af_posterior_check_design() accepts joint arms (S = S_s, no shared
#     factor warning), refuses mixed cycle / joint arms and S != S_s.
cat("\n[4] 09 design check\n")
cfg <- function(mode = "joint", S = 100L, h = "uniform") {
  list(
    simulator = list(
      max_draws = S,
      forecast = list(param_draw_mode = mode, handoff_weighting = h)
    )
  )
}
d_j <- tibble::tibble(
  replicate_id = 1:10000,
  draw = ((1:10000 - 1L) %% 100L) + 1L
)
w <- warns(
  chk <- af_posterior_check_design(
    "u1",
    list(factual = cfg(), neutral = cfg()),
    list(factual = d_j, neutral = d_j),
    n_smoothing = 100L
  )
)
expect(
  identical(attr(chk, "max_draws"), 100L) && nrow(chk) == 10000L,
  "joint arms pass with S = S_s = 100"
)
expect(
  is.null(w),
  "joint: no shared-factor warning (draw and state paired by design)"
)
expect(
  !is.null(err(af_posterior_check_design(
    "u1",
    list(factual = cfg(), neutral = cfg("cycle")),
    list(factual = d_j, neutral = d_j)
  ))),
  "mixed cycle / joint arms are refused"
)
expect(
  !is.null(err(af_posterior_check_design(
    "u1",
    list(factual = cfg(S = 200L)),
    list(factual = d_j),
    n_smoothing = 100L
  ))),
  "joint with max_draws != S_s is refused"
)
expect(
  identical(.sim_param_draw_mode(cfg()$simulator), "joint"),
  ".sim_param_draw_mode() returns joint"
)
cfg_n <- cfg()
cfg_n$simulator$n_replicates <- 10000L
expect(
  identical(run_design_row("forecast", cfg_n)$reps_per_draw, 100),
  "run_design_row(): reps_per_draw = n / S under joint"
)

# --- [5] event generator ------------------------------------------------------
# [5] event generator: event_param_draw_mode("joint") needs uniform handoff
#     and max_draws = the smoothing's nrep.
cat("\n[5] event generator\n")
ecfg <- function(S = 100L, h = "uniform") {
  list(
    simulator = list(
      max_draws = S,
      forecast = list(param_draw_mode = "random", handoff_weighting = h)
    )
  )
}
expect(
  identical(
    ev$event_param_draw_mode(
      ecfg(),
      "joint",
      100L
    )$simulator$forecast$param_draw_mode,
    "joint"
  ),
  "joint with uniform handoff and max_draws = nrep is set"
)
expect(
  !is.null(err(ev$event_param_draw_mode(ecfg(h = "shuffled"), "joint", 100L))),
  "joint refuses a shuffled handoff"
)
expect(
  !is.null(err(ev$event_param_draw_mode(ecfg(S = 200L), "joint", 100L))),
  "joint refuses max_draws != nrep"
)
expect(
  !is.null(err(ev$event_param_draw_mode(ecfg(), "joint", NA))),
  "joint needs the smoothing's nrep"
)

# --- [6] defaults ---------------------------------------------------------------
# [6] defaults: filter.draw_design per_tu_random in get_default_options()
#     and config_defaults.yml; param_draw_mode default unchanged (random).
cat("\n[6] defaults\n")
dflt <- get_default_options()$simulator
yml <- yaml::read_yaml("analysis/config_defaults.yml")$simulator
expect(
  identical(dflt$filter$draw_design, "per_tu_random") &&
    identical(yml$filter$draw_design, "per_tu_random"),
  "filter.draw_design defaults to per_tu_random in both"
)
expect(
  identical(dflt$forecast$param_draw_mode, "random") &&
    identical(yml$forecast$param_draw_mode, "random"),
  "param_draw_mode default unchanged (random)"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
