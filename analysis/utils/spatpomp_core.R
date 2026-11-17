# spatpomp_core.R — shared library of the spatPomp block particle filter of
#   the outbreak simulator
#
# Provides:
#   build_obs_for_tu(), build_covar_for_tu(), make_block_list()   filter inputs
#   run_bpfilter_tu(), run_one_bpfilter_replicate()   bpfilter runners
#   filtered_states_to_warmup(), filter_traj_to_warmup()   filter -> warm-up
#   select_handoff_warmup(), forecast_from_filtered()   handoff and forecast
#   load_smoothing_distribution(), smoothing_cache_meta(),
#   check_smoothing_cache_meta()   the smoothing cache
#   load_partition_map(), assert_partition_column()   grouping subsets
#
# Details:
#   Observation/covariate tables, measurement Csnippets, grouping subsets,
#   bpfilter runners, handoff and the smoothing cache. Sourced by
#   spatpomp_core_c.R (whose build_spatpomp_tu_c() builds the model),
#   07b_simulator_spatpomp.R and the smoothing modules.
#
# See: the "Measurement Csnippets" section for the spatPomp probe findings

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(spatPomp)
  library(pomp)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Evaluate `expr` while saving/restoring the global RNG state, so a seeded draw
# does not perturb the caller's RNG stream. Same idiom as sim_io.R:427-433
# (kept local — modules never source() each other).
.with_preserved_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  force(expr)
}

# Offset added to the per-TU `seed` before the per-replicate draw-index seeding
# in .forecast_param_indices()'s "random" mode, so the draw-selection RNG stream
# never collides with the kernel-noise stream (set.seed(seed + rep_id) inside
# simulate_one_tu(), simulator_core.R).
.DRAW_SEED_OFFSET <- 100000000L

# Offset added to simulator.seed before the per-block permutation seeding of
# select_handoff_warmup()'s "shuffled" mode, so the handoff-pairing stream never
# collides with the kernel-noise or draw-selection streams above.
.HANDOFF_SEED_OFFSET <- 200000000L

# Smoothing draw design (simulator.filter.draw_design): how
# the smoothing replicates pick their posterior draw. "per_tu_random" (default)
# = per unit, at random with replacement, seeded
# by (seed, unit); "shared_cycle" = smoothing replicate j runs draw j in every
# unit, so a forward run can inherit the draw of its handoff state
# (param_draw_mode "joint").
.SM_DRAW_DESIGNS <- c("per_tu_random", "shared_cycle")

#' The smoothing draw design of a `simulator.filter` block, validated.
smoothing_draw_design_of <- function(flt) {
  v <- (flt %||% list())[["draw_design"]] %||% "per_tu_random"
  if (
    !is.character(v) || length(v) != 1L || is.na(v) || !v %in% .SM_DRAW_DESIGNS
  ) {
    stop(
      "simulator.filter.draw_design must be one of: ",
      paste(.SM_DRAW_DESIGNS, collapse = ", "),
      "; got: ",
      paste(format(v), collapse = ", "),
      call. = FALSE
    )
  }
  v
}

#' The posterior draw of each smoothing replicate of one unit.
#'
#' "per_tu_random": set.seed(seed_root + position of `tu` in sort(tus)), then
#' n_reps draws at random with replacement from 1..max_draw_s (the default
#' assignment, used by both 08 and 07b's native filter). "shared_cycle":
#' replicate j runs draw j in every unit; n_reps must equal max_draw_s, so the
#' smoothing uses every draw of the thinned bundle exactly once (a shorter
#' prefix would take the first chains only). Consumes no RNG.
#' @param design .SM_DRAW_DESIGNS
#' @param tus    the full unit universe (the seed is keyed on it)
#' @return integer(n_reps) in 1..max_draw_s
smoothing_draw_idx <- function(design, seed_root, tu, tus, n_reps, max_draw_s) {
  design <- match.arg(design, .SM_DRAW_DESIGNS)
  if (identical(design, "shared_cycle")) {
    if (as.integer(n_reps) != as.integer(max_draw_s)) {
      stop(
        sprintf(
          paste0(
            "simulator.filter.draw_design shared_cycle: %d smoothing replicates but %d ",
            "posterior draws; set simulator.max_draws to simulator.n_replicates so ",
            "replicate j runs draw j and every draw is used once."
          ),
          as.integer(n_reps),
          as.integer(max_draw_s)
        ),
        call. = FALSE
      )
    }
    return(seq_len(as.integer(n_reps)))
  }
  set.seed(seed_root + match(tu, sort(tus)))
  sample.int(max_draw_s, n_reps, replace = TRUE)
}

# ---------------------------------------------------------------------------
# State-vector constants
# ---------------------------------------------------------------------------
# 11 unit-statenames; each expands to name1..nameU (spatPomp convention).
# Order matters for the packed/unpacked state in the step function.
# `mu` is the monthly-case (`mon`) NegBin mean for this month, carried as state
# so the joint dunit Csnippet can read it (R-function measures are
# ignored by bpfilter). It is held at 0 under the occurrence_bernoulli
# measurements (no extra work; constant shape).
# `size_phi` is the mon NegBin dispersion, carried as STATE (not covariate)
# because spatPomp unit-measure Csnippets expose states bare per unit but NOT
# covariates (covariates only appear as covarname+unitindex, e.g. "size_phi3").
# Seeded once in rinit; re-emitted unchanged by rprocess. Held at 0 under the
# Bernoulli path. The slot names come from the episode-total size track, retired
# downstream (mon replaces it). Names and positions are kept because the smoothing
# cache schema token is
# length(.SP_UNIT_STATENAMES) and the C variant indexes ST(7,u)/ST(8,u).
# `active`/`mso` (duration track): pure unobserved latent state — the
# ongoing-episode indicator and months-since-onset clock, as in the v1
# kernel. Appended, never inserted, so the Csnippet C variant's
# 0-based positional index macros (spatpomp_core_c.R) keep indices 0..8 stable.
# The measurement model is intentionally unchanged: .SP_DUNIT* still scores
# obs_Y against latent `p`/`Yout` only, never `active`/`mso` — zero Csnippet
# change, zero logLik-scale change, no invalidation of resample_by_loglik()
# weights or existing smoothing comparisons. `active`/`mso` default to 0 (held there
# by rinit/rprocess) whenever no duration draws bundle is supplied, so every
# existing (non-duration) caller is unaffected by the two extra slots.
#
# Monthly-case (`mon`) track: adds NO new slot. `mu` and `size_phi` are exactly
# a mean and a dispersion, which is what the mon measurement needs, so the mon
# path WRITES INTO THOSE TWO: mu is the per-active-month expected case count and
# size_phi the mon fit's dispersion. Consequences, all deliberate:
# .SM_CACHE_SCHEMA_VERSION (below) is unchanged so every existing smoothing
# handoff cache stays valid with no forced rebuild; the C variant's positional
# ST(i,u) macros are untouched; and `active`/`mso` (already maintained by the
# duration promotion block) are read directly by the mon measurement.
.SP_UNIT_STATENAMES <- c(
  "u_days",
  "o",
  "o_ct",
  "o_tu",
  "Yout",
  "draw_idx",
  "p",
  "mu",
  "size_phi",
  "active",
  "mso"
)

# `covariate_feed: detected` appends ONE slot, `det` (the per-episode detection
# flag, simulator_core.R::update_detection()), at position 12 / C index 11.
# Appended -- never inserted -- so every positional ST(k,u) macro and every
# existing consumer is unchanged, and it exists ONLY under the detected feed:
# a latent run keeps the 11-slot layout, the sv11 cache token, and byte-identical
# trajectories. The layout of any saved trajectory is self-describing (its
# first-dim names), so readers ask .sp_has_det() instead of assuming.
.SP_DET_STATENAME <- "det"
# The WHO-dashboard measurement (`.SP_WHO_MEAS`) appends TWO more slots after
# everything else: `cm`, this month's SIMULATED case count of the unit (the C
# twin of simulator_core.R::draw_monthly_cases(), drawn after mu), and `ccm`, the
# country sum of `cm`, replicated on each unit of the country so the per-unit
# dunit can score the national WHO count against it. Both are filter-step
# scratch: Phase B reads `cm` from the in-memory trajectory, then 08b strips
# both before states.qs is written (.sp_strip_scratch_slots()), so smoothing
# caches keep the sv11/sv12 layouts and no cache consumer changes.
.SP_WHO_STATENAMES <- c("cm", "ccm")
.sp_unit_statenames <- function(covariate_feed = "latent", who = FALSE) {
  c(
    .SP_UNIT_STATENAMES,
    if (identical(covariate_feed, "detected")) .SP_DET_STATENAME,
    if (isTRUE(who)) .SP_WHO_STATENAMES
  )
}
.sp_has_det <- function(state_names) {
  any(startsWith(state_names, .SP_DET_STATENAME))
}
#' Drop the WHO scratch slots (`cm<u>`, `ccm<u>`) from a filter_traj array.
#' Identity for every array without them.
.sp_strip_scratch_slots <- function(ft) {
  if (is.null(ft)) {
    return(ft)
  }
  rn <- dimnames(ft)[[1L]]
  keep <- !grepl(
    sprintf("^(%s)[0-9]+$", paste(.SP_WHO_STATENAMES, collapse = "|")),
    rn
  )
  if (all(keep)) {
    return(ft)
  }
  ft[keep, , , drop = FALSE]
}
#' Number of units in a saved per-unit state array, from its own names.
.sp_n_units <- function(state_names) sum(grepl("^Yout[0-9]+$", state_names))

# ---------------------------------------------------------------------------
# smoothing handoff-cache schema versioning
# ---------------------------------------------------------------------------
# The smoothing cache (hist_filter_tu=<tu>_*.qs, written by
# 08c_smoothing_compile.R) has no statename manifest, so a cache written
# before the .SP_UNIT_STATENAMES 9->11 change would silently misalign
# named-vector state lookups in downstream consumers. `.SM_CACHE_SCHEMA_VERSION` is
# derived from `length(.SP_UNIT_STATENAMES)`, not a hand-incremented counter,
# so it moves in lockstep with the statenames vector automatically.
#
# The empirical baseline-rate anchor (`offset(log(base_rate_<track>))`,
# utils/anchor.R) deliberately does NOT move it. The anchor is a STATIC
# per-unit constant, exactly like `pop`: the C builder folds it into the `*_lp_static`
# covariate column that already exists (precompute_*_lp_static(),
# spatpomp_core_c.R). No new state slot, no new covariate column, no change to
# the positional ST(i,u) macros -- so every existing smoothing handoff cache stays
# valid (the precompute-at-R-time design of the C variant is what avoids a bump).
.SM_CACHE_SCHEMA_VERSION <- length(.SP_UNIT_STATENAMES)
.SM_CACHE_SCHEMA_VERSIONS <- c(
  .SM_CACHE_SCHEMA_VERSION,
  length(.sp_unit_statenames("detected"))
)

#' Build the `_sv<NN>` schema token embedded in smoothing handoff-cache filenames.
#'
#' Filename layout: `hist_filter_tu=<tu>_Np=<Np>_sv<NN>_nrep=<N>_<start8>_<end8>.qs`
#' — the token sits *after* `_Np=<Np>_`, not between `tu=` and `_Np=`, so it
#' cannot be folded into the non-greedy TU-name capture at
#' `07b_simulator_spatpomp.R` (`^hist_filter_tu=(.*?)_Np=.*$`). Both existing
#' discovery globs (`hist_filter_tu=<tu>_*.qs`, `hist_filter_tu=*_Np=*.qs`)
#' still match old (token-less) and new (token-bearing) filenames alike —
#' discovery stays backward-compatible while a schema mismatch becomes
#' visible on `ls` and fails loudly at load time
#' (`check_smoothing_cache_schema()`).
sm_cache_schema_token <- function(covariate_feed = "latent") {
  sprintf("sv%02d", length(.sp_unit_statenames(covariate_feed)))
}

#' Parse the `_sv<NN>` schema token out of a cache filename, if present.
#'
#' Returns `NA_integer_` for filenames written before schema versioning was
#' introduced — these predate the .SP_UNIT_STATENAMES 9->11 change
#' too and must be treated as stale, not assumed compatible.
sm_cache_schema_token_parse <- function(path) {
  bn <- basename(path)
  m <- regmatches(bn, regexpr("_sv([0-9]+)_", bn))
  if (length(m) == 0L || !nzchar(m)) {
    return(NA_integer_)
  }
  as.integer(sub("^_sv([0-9]+)_$", "\\1", m))
}

#' Validate a smoothing handoff-cache filename against the current state schema.
#'
#' Fails loudly (a) when the filename predates schema versioning (no
#' `_sv<NN>` token — treated as stale, since it also predates the 9->11
#' state-vector change) and (b) when the embedded version does not match
#' `.SM_CACHE_SCHEMA_VERSION`. Intended to be the *first* thing that fires on
#' a stale smoothing cache — called by `load_smoothing_distribution()`
#' (`07b_simulator_spatpomp.R`) right after a cache path is resolved and
#' before it is read, so the "rebuild the smoothing cache" instruction reaches the user
#' immediately rather than surfacing several calls later as a misleading
#' unit-count mismatch.
check_smoothing_cache_schema <- function(path, tu = NULL) {
  sv <- sm_cache_schema_token_parse(path)
  who <- if (is.null(tu)) "" else sprintf("TU '%s': ", tu)
  if (is.na(sv)) {
    stop(
      sprintf(
        "%ssmoothing cache %s predates schema versioning (no _sv<NN> token) -- rebuild the smoothing cache (08a/08b/08c_smoothing_*.R) before use.",
        who,
        basename(path)
      ),
      call. = FALSE
    )
  }
  # sv11 = latent layout, sv12 = latent + det (covariate_feed: detected). Both
  # are current; the reader picks the layout up from the trajectory names.
  if (!(sv %in% .SM_CACHE_SCHEMA_VERSIONS)) {
    stop(
      sprintf(
        "%ssmoothing cache %s has schema sv%02d but this build expects sv%02d (.SP_UNIT_STATENAMES has %d entries) -- rebuild the smoothing cache (08a/08b/08c_smoothing_*.R).",
        who,
        basename(path),
        sv,
        .SM_CACHE_SCHEMA_VERSION,
        .SM_CACHE_SCHEMA_VERSION
      ),
      call. = FALSE
    )
  }
  invisible(sv)
}

# ---------------------------------------------------------------------------
# Measurement-class predicates (used by the Csnippet builder and the smoothing)
# ---------------------------------------------------------------------------
# Orthogonal axes:
#   - State    : the occurrence likelihood conditions on the LATENT state `Yout`
#                (near-delta vs obs_Y), not the deterministic hazard `p`. The
#                non-state variants score `p`, which is (nearly) constant across
#                particles and therefore cannot pull the filtered occurrence
#                trajectory toward the observations — the cause of the systemic
#                occurrence over-estimation. The `_state` variants fix this by
#                making obs_Y a (noisy) measurement of the realized Yout.
#   - Monthly  : the case term is scored at EVERY ACTIVE month against the mon
#                fit's per-month mu. The case term is gated on `mu > 0.0` (a
#                cessation-order-safe proxy for "active this month" -- see the
#                .SP_DUNIT_MON_NB comment), so it hard-requires the duration
#                track (nothing else maintains `mu` on a per-month basis).
# The episode-total NegBin (size track) measurement family is not provided
# (the size track is retired downstream); its names fail .check_measurement()
# like any other unknown string.
# "act" (activity-scored): the detect sibling whose case term is scored in
# EVERY month carrying an observed monthly count -- including post-onset months
# whose obs_Y is NA (the survival panel censors them) -- and whose zero counts
# (build_obs_panel(zero_fill = TRUE)) penalise a latent episode that is still
# running. See .sp_case_actpen_zi_c() / .make_dunit_mon(score_case_if_na_y=TRUE).
# "actn" (normalized): the `_act` sibling whose INACTIVE-state case
# term is a proper mixture (total mass 1 for every phi) instead of the
# 0-or-floor-NB pair whose mass is 2 - NB(0) -- see .sp_case_actpen_zi_norm_c()
# (and .sp_case_who_norm_c() for the `_who` form). Part of .SP_ACT_MEAS, so
# every `_act` predicate (detect, state, mon, monthly_case_obs guard) applies.
.SP_ACTN_MEAS <- c(
  "occurrence_duration_mon_NB_state_detect_actn",
  "occurrence_duration_mon_NB_state_detect_actn_who"
)
.SP_ACT_MEAS <- c(
  "occurrence_duration_mon_NB_state_detect_act",
  "occurrence_duration_mon_NB_state_detect_act_who",
  .SP_ACTN_MEAS
)
# "who": the `_act` measurement plus a country-level term scoring the WHO
# dashboard national monthly count (obs_W) against the particle's SIMULATED
# country case count (`ccm`) -- see .sp_case_who_c(). Needs filter.who_obs.
.SP_WHO_MEAS <- c(
  "occurrence_duration_mon_NB_state_detect_act_who",
  "occurrence_duration_mon_NB_state_detect_actn_who"
)
.SP_MON_MEAS <- c(
  "occurrence_duration_mon_NB",
  "occurrence_duration_mon_NB_state",
  "occurrence_duration_mon_NB_state_actpen",
  "occurrence_duration_mon_NB_state_detect",
  .SP_ACT_MEAS
)
.SP_STATE_MEAS <- c(
  "occurrence_bernoulli_state",
  "occurrence_bernoulli_state_detect",
  "occurrence_duration_mon_NB_state",
  "occurrence_duration_mon_NB_state_actpen",
  "occurrence_duration_mon_NB_state_detect",
  .SP_ACT_MEAS
)
# "actpen" (activity-penalized): floors the NegBin mean when a particle has
# wrongly ceased so a real observed count is still graded (not skipped) --
# see .SP_DUNIT_MON_NB_STATE_ACTPEN.
.SP_ACTPEN_MEAS <- c("occurrence_duration_mon_NB_state_actpen")
# "detect" (detection-error): replaces the near-delta onset-match term with an
# explicit sensitivity/specificity imperfect-detection measurement model,
# decoupled from the fitted hazard `p` -- see .sp_onset_detect_c(). One detect
# sibling per state-conditioning family; the mon one is built on _ACTPEN's
# floor-based case term, which strictly dominates plain _MON_NB_STATE's.
.SP_DETECT_MEAS <- c(
  "occurrence_bernoulli_state_detect",
  "occurrence_duration_mon_NB_state_detect",
  .SP_ACT_MEAS
)
.SP_ALL_MEAS <- c(
  "occurrence_bernoulli",
  "occurrence_bernoulli_state",
  "occurrence_bernoulli_state_detect",
  .SP_MON_MEAS
)
.is_mon_meas <- function(m) m %in% .SP_MON_MEAS
.is_state_meas <- function(m) m %in% .SP_STATE_MEAS
.is_actpen_meas <- function(m) m %in% .SP_ACTPEN_MEAS
.is_detect_meas <- function(m) m %in% .SP_DETECT_MEAS
.is_act_meas <- function(m) m %in% .SP_ACT_MEAS
.is_who_meas <- function(m) m %in% .SP_WHO_MEAS
.is_actn_meas <- function(m) m %in% .SP_ACTN_MEAS

# Namespace prefix for the mon track's PANEL cross-basis covariate columns.
# The covariate table is one flat namespace shared by every track, but occurrence
# and monthly cases pick their ENSO single-lag independently (the deployed best-lag
# fits differ per track and per group), so their `enso_cb<k>` columns are two
# different fitted bases under one name. Prefixing the mon copy keeps both in the
# table; the C builder maps the prefixed column back to the bare model term when it
# assembles the mon linear predictor. Chosen to contain no digits, so it cannot
# alter spatPomp's paste0(colname, unit_index) collision behaviour.
.SP_MON_CB_PREFIX <- "monCB_"

# Same rationale, same constraint (no digits), for the duration track. Every
# track picks its ENSO single-lag independently -- for the Horn group the
# deployed best-lag fits are occ=11, mon=11, dur=8 -- so `enso_cb1` means
# different fitted bases. Before this prefix existed the duration linear
# predictor read the BARE (occurrence) columns, i.e. it evaluated its own
# coefficients against the occurrence model's basis at the occurrence model's
# lag. Only `mon` was namespaced.
.SP_DUR_CB_PREFIX <- "durCB_"

#' Validate a `simulator.filter.measurement` string against the supported enum.
#'
#' Called by every builder so an unknown/typo'd measurement fails loudly at build
#' time naming the valid set, rather than silently degrading to the Bernoulli
#' path (which is what an unguarded `%in%` predicate chain would do — every
#' predicate returns FALSE and the run looks superficially fine).
.check_measurement <- function(m, where = "build_spatpomp_tu_c") {
  if (!(length(m) == 1L && is.character(m) && m %in% .SP_ALL_MEAS)) {
    stop(sprintf(
      "%s: unknown measurement '%s'. Supported: %s.",
      where,
      paste(m, collapse = ","),
      paste(.SP_ALL_MEAS, collapse = ", ")
    ))
  }
  invisible(m)
}

#' Resolve the measurement-side `simulator.filter.*` keys, ONCE, for both
#' drivers (historical smoothing 08 via smoothing_make_config(), forecast 07b).
#' Previously each driver repeated its own `flt$key %||% default` lines; every
#' new measurement key is added here and nowhere else.
#'
#'   measurement       default "occurrence_duration_mon_NB" (both drivers' default)
#'   delta_fn/delta_fp detection-error rates (.SP_STATE_DELTA_FN/_FP defaults)
#'   covariate_feed    "latent" (default) | "detected": which activity field the
#'                     kernel's spillover covariates (o, o_ct, o_tu, nbr, ttlo)
#'                     are built from -- see update_detection()
#'   monthly_case_obs  "none" (default) | "observed" | "zero_filled": whether the
#'                     filter's obs table carries observed monthly case counts
#'                     (obs_CM) and whether empty non-onset panel months are 0
#'                     (sim_io.R::build_obs_panel())
#'   mon_case_col      derived: "cases_month" unless monthly_case_obs = "none"
#'   case_fp_rate/case_fp_mean  the inactive-state mixture of the `_actn`
#'                     measurements (eps, mu_0; .SP_CASE_FP_RATE / .SP_MON_MU_FLOOR
#'                     defaults); inert for every other measurement
#'
#' Guards the "config computed but not threaded" bug class: the `_act`
#' measurement is inert without observed monthly cases, so asking for it with
#' monthly_case_obs = "none" is an error, not a silent no-op.
resolve_filter_cfg <- function(flt) {
  flt <- flt %||% list()
  fc <- list(
    measurement = flt$measurement %||% "occurrence_duration_mon_NB",
    delta_fn = as.numeric(flt$delta_fn %||% .SP_STATE_DELTA_FN),
    delta_fp = as.numeric(flt$delta_fp %||% .SP_STATE_DELTA_FP),
    case_fp_rate = as.numeric(flt[["case_fp_rate"]] %||% .SP_CASE_FP_RATE),
    case_fp_mean = as.numeric(flt[["case_fp_mean"]] %||% .SP_MON_MU_FLOOR),
    covariate_feed = flt$covariate_feed %||% "latent",
    monthly_case_obs = flt$monthly_case_obs %||% "none",
    traj_storage = flt$traj_storage %||% "pruned",
    traj_prune_every = as.integer(flt$traj_prune_every %||% 1L)
  )
  fc$traj_storage <- .check_traj_storage(fc$traj_storage)
  if (
    length(fc$traj_prune_every) != 1L ||
      is.na(fc$traj_prune_every) ||
      fc$traj_prune_every < 1L
  ) {
    stop(
      "simulator.filter.traj_prune_every must be a positive integer.",
      call. = FALSE
    )
  }
  if (!fc$covariate_feed %in% c("latent", "detected")) {
    stop(
      "simulator.filter.covariate_feed must be 'latent' or 'detected', got '",
      fc$covariate_feed,
      "'.",
      call. = FALSE
    )
  }
  if (!fc$monthly_case_obs %in% c("none", "observed", "zero_filled")) {
    stop(
      "simulator.filter.monthly_case_obs must be 'none', 'observed' or ",
      "'zero_filled', got '",
      fc$monthly_case_obs,
      "'.",
      call. = FALSE
    )
  }
  if (.is_act_meas(fc$measurement) && fc$monthly_case_obs == "none") {
    stop(
      "simulator.filter.measurement '",
      fc$measurement,
      "' scores observed ",
      "monthly cases; set simulator.filter.monthly_case_obs to 'zero_filled' ",
      "(or 'observed').",
      call. = FALSE
    )
  }
  fc$mon_case_col <- if (fc$monthly_case_obs == "none") NULL else "cases_month"
  fc$who_obs <- .resolve_who_obs_cfg(flt$who_obs, fc$measurement)
  fc
}

#' Resolve `simulator.filter.who_obs` (WHO dashboard national monthly counts).
#'
#' `enable` and the `_who` measurement go together: either one without the
#' other is an error (the measurement reads obs_W / ccm, which only exist
#' when the observation is wired in; an enabled block under another
#' measurement would be silently ignored). `phi` is the negative-binomial size
#' linking the simulated national count to the WHO count (reporting
#' overdispersion). Contribution shares are refused -- see
#' utils/who_dashboard.R::who_block_share(). The block-granularity rule is
#' checked where the effective granularity is known
#' (.check_who_obs_blocks(), called by smoothing_filter_one()).
.resolve_who_obs_cfg <- function(w, measurement) {
  w <- w %||% list()
  # [[ ]] throughout: `w$scale` would partial-match `scale_mode` / `scale_est`
  # when `scale` is absent.
  se <- w[["scale_est"]] %||% list()
  out <- list(
    enable = isTRUE(w[["enable"]]),
    path = w[["path"]] %||% "generated_data/who_dashboard_monthly_africa.rds",
    start = as.Date(w[["start"]] %||% "2023-01-01"),
    end = as.Date(w[["end"]] %||% "2025-01-01"),
    phi = as.numeric(w[["phi"]] %||% 5),
    mu_floor = as.numeric(w[["mu_floor"]] %||% 1),
    # eps of the `_actn_who` inactive-country mixture; inert otherwise.
    fp_rate = as.numeric(w[["fp_rate"]] %||% .SP_CASE_FP_RATE),
    scale = as.numeric(w[["scale"]] %||% 1),
    block_share = as.character(w[["block_share"]] %||% "pop"),
    allow_adm1_blocks = isTRUE(w[["allow_adm1_blocks"]]),
    # Per-country scale and the DB-overlap rule; inert unless
    # scale_mode = "estimated" / db_overlap = "skip_who".
    scale_mode = as.character(w[["scale_mode"]] %||% "fixed"),
    scale_est = list(
      start = as.Date(se[["start"]] %||% "2023-01-01"),
      end = as.Date(se[["end"]] %||% "2025-01-01"),
      min_db_cases = as.numeric(se[["min_db_cases"]] %||% 100),
      min_months = as.integer(se[["min_months"]] %||% 6L),
      shrink_cases = as.numeric(se[["shrink_cases"]] %||% 0),
      fallback = as.character(se[["fallback"]] %||% "pooled"),
      clamp = as.numeric(unlist(se[["clamp"]] %||% c(0.1, 10)))
    ),
    db_overlap = as.character(w[["db_overlap"]] %||% "score_both")
  )
  if (!out$scale_mode %in% c("fixed", "estimated")) {
    stop(
      "simulator.filter.who_obs.scale_mode must be 'fixed' or 'estimated', got '",
      out$scale_mode,
      "'.",
      call. = FALSE
    )
  }
  if (!out$db_overlap %in% c("score_both", "skip_who")) {
    stop(
      "simulator.filter.who_obs.db_overlap must be 'score_both' or 'skip_who', got '",
      out$db_overlap,
      "'.",
      call. = FALSE
    )
  }
  if (identical(out$scale_mode, "estimated")) {
    s <- out$scale_est
    if (anyNA(c(s$start, s$end)) || s$end <= s$start) {
      stop(
        "simulator.filter.who_obs.scale_est: start/end must be dates with start < end.",
        call. = FALSE
      )
    }
    if (!s$fallback %in% c("pooled", "fixed")) {
      stop(
        "simulator.filter.who_obs.scale_est.fallback must be 'pooled' or 'fixed'.",
        call. = FALSE
      )
    }
    if (
      length(s$clamp) != 2L ||
        !all(is.finite(s$clamp)) ||
        s$clamp[1] <= 0 ||
        s$clamp[1] >= s$clamp[2]
    ) {
      stop(
        "simulator.filter.who_obs.scale_est.clamp must be [lo, hi] with 0 < lo < hi.",
        call. = FALSE
      )
    }
    if (
      !is.finite(s$min_db_cases) ||
        s$min_db_cases < 0 ||
        is.na(s$min_months) ||
        s$min_months < 1L ||
        !is.finite(s$shrink_cases) ||
        s$shrink_cases < 0
    ) {
      stop(
        "simulator.filter.who_obs.scale_est: min_db_cases / shrink_cases must be >= 0 and ",
        "min_months >= 1.",
        call. = FALSE
      )
    }
  }
  if (out$enable && !.is_who_meas(measurement)) {
    stop(
      "simulator.filter.who_obs.enable is true but measurement '",
      measurement,
      "' does not score it; use ",
      .SP_WHO_MEAS,
      ".",
      call. = FALSE
    )
  }
  if (!out$enable && .is_who_meas(measurement)) {
    stop(
      "simulator.filter.measurement '",
      measurement,
      "' scores the WHO ",
      "dashboard counts; set simulator.filter.who_obs.enable: true.",
      call. = FALSE
    )
  }
  if (identical(out$block_share, "contribution")) {
    stop(
      "simulator.filter.who_obs.block_share 'contribution' is refused: a ",
      "state-dependent share multiplying a log-pmf <= 0 favours particles that ",
      "simulate no cases over one that matches the WHO count. Use 'pop' or 'equal'.",
      call. = FALSE
    )
  }
  if (!out$block_share %in% c("pop", "equal")) {
    stop(
      "simulator.filter.who_obs.block_share must be 'pop' or 'equal', got '",
      out$block_share,
      "'.",
      call. = FALSE
    )
  }
  if (out$enable) {
    if (anyNA(c(out$start, out$end)) || out$end <= out$start) {
      stop(
        "simulator.filter.who_obs: start/end must be dates with start < end.",
        call. = FALSE
      )
    }
    if (
      !is.finite(out$phi) ||
        out$phi <= 0 ||
        !is.finite(out$mu_floor) ||
        out$mu_floor <= 0 ||
        !is.finite(out$scale) ||
        out$scale <= 0
    ) {
      stop(
        "simulator.filter.who_obs: phi, mu_floor and scale must be > 0.",
        call. = FALSE
      )
    }
  }
  out
}

#' The WHO term must enter each country's bpfilter block exactly once. With
#' `block_granularity: country` it does. With ADM1 blocks each block gets a
#' tempered share (NB^w_u): the recombined particle's country sum is not the
#' one scored, inactive blocks are resampled by noise from other blocks'
#' states, and the term's strength falls with the number of ADM1s. Allowed
#' only as a stated sensitivity choice (who_obs.allow_adm1_blocks).
.check_who_obs_blocks <- function(who_obs, block_granularity, frame_tu = NULL) {
  if (!isTRUE(who_obs$enable) || identical(block_granularity, "country")) {
    return(invisible(TRUE))
  }
  if (!isTRUE(who_obs$allow_adm1_blocks)) {
    stop(
      "simulator.filter.who_obs needs block_granularity: country (got '",
      block_granularity,
      "'); set who_obs.allow_adm1_blocks: true to split ",
      "the country term across ADM1 blocks as a sensitivity run.",
      call. = FALSE
    )
  }
  kc <- if (!is.null(frame_tu)) table(as.character(frame_tu$country)) else NULL
  warning(
    "who_obs with ADM1 blocks: the country term is tempered per block ",
    "(strength ~ 1/K_c) and resamples inactive blocks by noise.",
    if (!is.null(kc)) {
      paste0(" K_c: ", paste(names(kc), kc, sep = "=", collapse = ", "))
    },
    call. = FALSE
  )
  invisible(FALSE)
}

#' The filter's observation panel for a resolved filter config: the shared
#' sim_io.R::build_obs_panel() with the monthly frame read only when needed.
#'
#' `surv_path` (optional): the path `surv_dat` was read from. When given and
#' the monthly frame is read, the two must share one episode definition
#' (utils/provenance.R::assert_episode_coherence(), silent for two historical
#' lp_season artefacts) -- obs_Y and obs_CM are then scored on the same
#' outbreaks. NULL keeps the previous behaviour exactly.
obs_panel_for_filter <- function(
  surv_dat,
  fc,
  mon_data_path,
  surv_path = NULL
) {
  if (identical(fc$monthly_case_obs, "none")) {
    return(build_obs_panel(surv_dat))
  }
  if (!file.exists(mon_data_path)) {
    stop(
      "monthly_case_obs = '",
      fc$monthly_case_obs,
      "' needs the monthly case ",
      "frame at ",
      mon_data_path,
      " (03_prepare_mon.R).",
      call. = FALSE
    )
  }
  if (!is.null(surv_path)) {
    assert_episode_coherence(
      list(
        surv = episode_stamp_of(surv_path),
        mon = episode_stamp_of(mon_data_path)
      ),
      what = "the filter's observation panels"
    )
  }
  build_obs_panel(
    surv_dat,
    readRDS(mon_data_path),
    zero_fill = identical(fc$monthly_case_obs, "zero_filled")
  )
}

#' Validate the detection-error rates spliced into a `_detect` dunit Csnippet.
#'
#' Both rates enter the onset term through `log()`, so a value of 0, 1 or outside
#' [0,1] would emit `-Inf`/`NaN` into the particle weights and silently collapse
#' the filter rather than erroring. Called from `.sp_onset_detect_c()` — i.e.
#' exactly when the values are used, so it stays inert for every non-detect
#' measurement (whose builders never touch delta_fn/delta_fp).
.check_deltas <- function(delta_fn, delta_fp, where = "build_spatpomp_tu_c") {
  ok <- function(v) {
    length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0 && v < 1
  }
  if (!ok(delta_fn)) {
    stop(sprintf(
      "%s: delta_fn must be a single number in (0,1); got %s",
      where,
      paste(format(delta_fn), collapse = ",")
    ))
  }
  if (!ok(delta_fp)) {
    stop(sprintf(
      "%s: delta_fp must be a single number in (0,1); got %s",
      where,
      paste(format(delta_fp), collapse = ",")
    ))
  }
  invisible(TRUE)
}

#' Validate the inactive-state mixture of the `_actn` case terms:
#' eps enters through log(eps) and log((1 - eps) + eps NB(0)), so it must be in
#' (0, 1]; the floor mean must be finite and > 0. Called from the `_norm`
#' factories only, so it stays inert for every other measurement.
.check_case_fp <- function(eps, mu_0, where = "build_spatpomp_tu_c") {
  if (
    !(length(eps) == 1L &&
      is.numeric(eps) &&
      is.finite(eps) &&
      eps > 0 &&
      eps <= 1)
  ) {
    stop(sprintf(
      "%s: case_fp_rate must be a single number in (0,1]; got %s",
      where,
      paste(format(eps), collapse = ",")
    ))
  }
  if (
    !(length(mu_0) == 1L && is.numeric(mu_0) && is.finite(mu_0) && mu_0 > 0)
  ) {
    stop(sprintf(
      "%s: case_fp_mean must be a single number > 0; got %s",
      where,
      paste(format(mu_0), collapse = ",")
    ))
  }
  invisible(TRUE)
}

# Sentinel for NA in o_tu (spatPomp state must be finite double; NA not allowed).
.SP_OTU_SENTINEL <- -1.0

# ---------------------------------------------------------------------------
# Build observation table — long, one row per (month_index, gid)
# ---------------------------------------------------------------------------

#' Footprint-aligned count aggregation -- the ONE implementation of the rule
#' shared by the filter's observed-count channel (obs_CM monthly cases) and by
#' sim_io.R::build_obs_panel(): sum the counts
#' WITHIN each admin level, then take the MAX across per-level totals, i.e. the
#' finest-resolution footprint total without double-counting an ADM2 report
#' against its ADM3 children (mirrors 03_prepare_siz.R's collapse). A frame
#' without an `admin_level` column degrades to a plain sum per key.
#'
#' @param df        data frame carrying `keys`, `value_col` and optionally
#'   `admin_level`.
#' @param value_col name of the count column.
#' @param out_col   name of the aggregated output column.
#' @param keys      grouping keys (default the filter's long-table keys).
#' @param keep_na   FALSE (default, the pre-existing rule): an all-NA level sums
#'   to 0. TRUE: an all-NA level stays NA, so a reported-but-unquantified month
#'   (e.g. a `cases_gap` row) reads "magnitude unobserved" rather than "zero".
#' @return tibble(keys..., out_col); non-finite results are NA.
.aggregate_footprint_counts <- function(
  df,
  value_col,
  out_col,
  keys = c("time", "unit"),
  keep_na = FALSE
) {
  alv <- if ("admin_level" %in% names(df)) "admin_level" else NULL
  lvl_sum <- function(x) {
    x <- as.numeric(x)
    if (keep_na && all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  }
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(keys, alv)))) %>%
    dplyr::summarise(.v = lvl_sum(.data[[value_col]]), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::summarise(
      !!out_col := suppressWarnings(max(.v, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      !!out_col := ifelse(
        is.finite(.data[[out_col]]),
        .data[[out_col]],
        NA_real_
      )
    )
}

#' Build the long observation table for spatPomp from the warmup panel.
#'
#' Aggregates outbreak occurrence as max(outbreak) per (gid, month) within
#' `window`. Month indices are 1-based integers (t0 = 0 is before the data).
#'
#' @param panel_tu_warmup  panel tibble for this TU, including the warmup tail.
#' @param frame_tu         location frame for this TU.
#' @param window           list(start = "YYYY-MM-DD", end = "YYYY-MM-DD"),
#'                         the filter (warmup) window — obs come from [start, end).
#' @param months           ordered Date vector of months in the window (t_start values).
#' @param obs_data         optional data frame with columns (gid, t_start, obs_Y)
#'                         pre-aggregated from the survival panel. When supplied,
#'                         outbreak occurrence is read from here rather than from
#'                         panel_tu_warmup; the panel need not have 'outbreak'/'obs_Y'.
#' @param mon_case_col optional character column name holding the per-(gid,
#'                   month) observed case count consumed by the monthly-case
#'                   (`mon`) measurement, emitted as the long column `obs_CM`.
#'                   Missing cells coalesce to NA, which the mon
#'                   hurdle treats as "active but magnitude unobserved" and
#'                   contributes 0 to the log-likelihood — matching the
#'                   `ISNA(obs_Y)` early return at the top of every dunit snippet.
#' @return long tibble: time (integer 1..M), unit (character gid), obs_Y (0/1),
#'         and obs_CM (this month's cases or NA) when `mon_case_col` is supplied.
#'         Column order is (obs_Y, obs_CM) so spatPomp's __obsindex is
#'         deterministic.
build_obs_for_tu <- function(
  panel_tu_warmup,
  frame_tu,
  window,
  months,
  obs_data = NULL,
  na_fill = NA_integer_,
  mon_case_col = NULL,
  who_tbl = NULL
) {
  wstart <- as.Date(window$start)
  wend <- as.Date(window$end)
  want_CM <- !is.null(mon_case_col)

  if (!is.null(obs_data)) {
    # External observation source (survival panel pre-aggregated by load_artefacts).
    # obs_data already has unique (gid, t_start) by max(outbreak).
    have_CM <- want_CM && (mon_case_col %in% names(obs_data))
    obs_raw <- obs_data %>%
      dplyr::filter(
        as.character(gid) %in% as.character(frame_tu$gid),
        t_start >= wstart,
        t_start < wend
      ) %>%
      dplyr::mutate(
        unit = as.character(gid),
        time = match(as.Date(t_start), months)
      ) %>%
      dplyr::filter(!is.na(time))
    obs_raw$obs_CM <- if (have_CM) {
      as.numeric(obs_raw[[mon_case_col]])
    } else {
      NA_real_
    }
    obs_raw <- dplyr::select(
      obs_raw,
      dplyr::all_of(
        c("time", "unit", "obs_Y", if (want_CM) "obs_CM")
      )
    )
  } else {
    obs_col <- if ("outbreak" %in% names(panel_tu_warmup)) {
      "outbreak"
    } else if ("obs_Y" %in% names(panel_tu_warmup)) {
      "obs_Y"
    } else {
      stop("build_obs_for_tu: panel must have column 'outbreak' or 'obs_Y'.")
    }
    have_CM <- want_CM && (mon_case_col %in% names(panel_tu_warmup))

    obs_pre <- panel_tu_warmup %>%
      dplyr::filter(
        gid %in% frame_tu$gid,
        t_start >= wstart,
        t_start < wend
      ) %>%
      dplyr::mutate(
        unit = as.character(gid),
        time = match(t_start, months) # integer 1..M
      ) %>%
      dplyr::filter(!is.na(time)) %>%
      dplyr::rename(obs_Y_raw = dplyr::all_of(obs_col))
    obs_raw <- obs_pre %>%
      dplyr::group_by(time, unit) %>%
      dplyr::summarise(
        obs_Y = as.integer(max(obs_Y_raw, na.rm = TRUE)),
        .groups = "drop"
      )

    # Monthly observed cases (mon track), footprint-aligned: sum within admin
    # level, then max across levels -- see .aggregate_footprint_counts().
    if (have_CM) {
      obs_raw <- dplyr::left_join(
        obs_raw,
        .aggregate_footprint_counts(obs_pre, mon_case_col, "obs_CM"),
        by = c("time", "unit")
      )
    } else if (want_CM) {
      obs_raw$obs_CM <- NA_real_
    }
  }

  # Ensure every (time, unit) combination is present (fill missing obs_Y with
  # na_fill; obs_CM always coalesces to NA -> "magnitude unobserved").
  grid <- tidyr::expand_grid(
    time = seq_along(months),
    unit = as.character(frame_tu$gid)
  )
  # .locale = "C" pins the unit axis to byte-order (ASCII) collation
  # regardless of the running session's system locale. This table's row
  # order IS the pomp's internal unit_names order, and every consumer that
  # re-derives that order downstream (state-array extraction in
  # smoothing_core.R / spatpomp_core_c.R / 07b_simulator_spatpomp.R) must use
  # sort(..., method = "radix") / order(..., method = "radix") to match it —
  # radix is R's C-locale-equivalent collation for character vectors. Under
  # en_US.UTF-8 (ICU) collation, "BDI_1_1" < "BDI_11_1" < "BDI_2_1", but under
  # C/radix collation "BDI_11_1" < "BDI_1_1" < "BDI_2_1" ('_' sorts before
  # '1' as a byte but after 'gid' digits are compared lexically) — a mismatch
  # here silently swaps two units' entire state trajectories.
  obs <- grid %>%
    dplyr::left_join(obs_raw, by = c("time", "unit")) %>%
    dplyr::mutate(obs_Y = dplyr::coalesce(obs_Y, na_fill)) %>%
    dplyr::arrange(time, unit, .locale = "C")

  if (want_CM) {
    if (!("obs_CM" %in% names(obs))) {
      obs$obs_CM <- NA_real_
    }
    obs <- dplyr::select(
      obs,
      dplyr::all_of(c("time", "unit", "obs_Y", "obs_CM"))
    )
  }

  # WHO dashboard national monthly count (the `_who` measurement): the SAME
  # value on every unit of the country in a scored month, NA elsewhere, so the
  # per-unit shares of the country term sum to the full term. Joined on the
  # grid, never on obs_raw, so it does not depend on survival-panel rows.
  # Column order (obs_Y, obs_CM, obs_W) fixes spatPomp's __obsindex.
  if (!is.null(who_tbl)) {
    ctry <- as.character(frame_tu$country)[match(
      obs$unit,
      as.character(frame_tu$gid)
    )]
    key_obs <- paste(ctry, format(as.Date(months[obs$time])))
    key_who <- paste(
      as.character(who_tbl$country),
      format(as.Date(who_tbl$t_start))
    )
    obs$obs_W <- as.numeric(who_tbl$cases_who)[match(key_obs, key_who)]
  }

  obs
}

# ---------------------------------------------------------------------------
# Build covariate table — long, one row per (month_index, gid), including t0
# ---------------------------------------------------------------------------

#' Build the spatPomp covariate table from the TU climate panel.
#'
#' All covariates are per-unit rows (shared_covarnames is unsupported alongside
#' per-unit covars in spatPomp v1.1.0, so ENSO is replicated across units).
#' Includes a t0 row (time = 0) required by spatPomp.
#'
#' @param panel_tu   panel tibble for this TU (warmup + horizon months).
#' @param frame_tu   location frame for this TU.
#' @param months     ordered Date vector of ALL months covered (warmup + horizon),
#'                   indexed 1..M. Month index 0 is t0 (the row before month 1).
#' @param enso_cols  optional character vector of ENSO column names to include.
#'                   When NULL (default), all enso_lag* columns in the panel are
#'                   used — but NOTE: spatPomp v1.1.0 internally forms covariate
#'                   names via paste0(colname, unit_index), so multi-lag panels
#'                   cause naming collisions (e.g. "enso_lag1" + unit 11 =
#'                   "enso_lag11" + unit 1 = "enso_lag111"). Pass only the lags
#'                   actually used by the model to avoid this.
#' @param cb_meta   optional DLNM cross-basis metadata (occurrence model). When
#'                  non-NULL, cross-basis design columns (`<var>_cb<k>`) are
#'                  reconstructed from the panel's lag columns via
#'                  build_cb_design() and appended to the covar table. ENSO
#'                  cross-bases use the panel's observed enso_lag* values
#'                  (deterministic in the filter window); climate cross-bases
#'                  use the `<var>_lag*` columns (emitted by 03c/03d).
#'                  NOTE: like enso_cols, the cross-basis column names are
#'                  suffixed with the unit index inside spatPomp — keep cb basis
#'                  dimension small to avoid paste0(colname, unit) collisions.
#' @param cb_meta_by_group optional named list of per-group cb_meta (occurrence
#'                  model) for grouped fits. When non-NULL and frame_tu carries a
#'                  `group` column, each unit's cross-basis columns are built from
#'                  ITS group's cb_meta — groups may carry different single_lags /
#'                  fitted bases, so a single shared cb_meta would mis-parse all
#'                  but one group. Takes precedence over `cb_meta`.
#' @param mon_harmonics logical; when TRUE, append the four calendar-month
#'                  harmonic columns (`sin_month`, `cos_month`, `sin2_month`,
#'                  `cos2_month`) required by the monthly-case (`mon`) track's
#'                  default `seasonal_term: "harmonic2"`, and by a duration fit
#'                  under `duration.seasonal_term: "harmonic2"` (the C builder
#'                  sets it TRUE for either). Built here rather than
#'                  derived from the step function's `t` on purpose: the covariate
#'                  table already carries the one-month index shift that gives
#'                  month k's transition month k's covariates (see the
#'                  "Covariate time-alignment" comment below), so building the
#'                  harmonics as covariates inherits that alignment instead of
#'                  re-deriving it and risking the exact off-by-one the seam-gap
#'                  investigation found. Default FALSE leaves the covariate table
#'                  of every non-mon run bit-for-bit unchanged.
#' @param cb_meta_mon,cb_meta_mon_by_group optional cb_meta for the monthly-case
#'                  (`mon`) fit, ungrouped / per-group. Its PANEL cross-bases
#'                  (everything except the latent `mso`, which is handled by
#'                  `build_mso_basis()` instead) are built exactly like the
#'                  occurrence ones but emitted under a `.SP_MON_CB_PREFIX`
#'                  prefix. The prefix is not cosmetic: the two tracks select
#'                  their ENSO single-lag independently (occ 4 / mon 9 for East,
#'                  say), so an unprefixed `enso_cb1` would be one column serving
#'                  two different fitted bases, and whichever was written last
#'                  would silently drive both linear predictors.
#' @param cb_meta_dur,cb_meta_dur_by_group same treatment for the duration fit,
#'                  under `.SP_DUR_CB_PREFIX`. This exists for the same reason
#'                  as the mon pair and is equally not optional: the deployed
#'                  per-group best-lag fits select ENSO lags per track (Horn:
#'                  occ 11, mon 11, dur 8), so before it was added the duration
#'                  linear predictor silently read the occurrence columns -- its
#'                  own coefficients evaluated against another track's basis at
#'                  another track's lag. NULL (the default) emits nothing and
#'                  leaves the covariate table of every occurrence-only run
#'                  bit-for-bit unchanged.
#' @param season_k  logical; TRUE appends the month-of-season column `seas_k`
#'                  read by the season onset gate (onset_eligibility: season).
#'                  FALSE (default) leaves the table unchanged.
#' @return long tibble sorted (time, unit): time (0..M), unit (gid), covariates.
build_covar_for_tu <- function(
  panel_tu,
  frame_tu,
  months,
  enso_cols = NULL,
  cb_meta = NULL,
  cb_meta_by_group = NULL,
  mon_harmonics = FALSE,
  cb_meta_mon = NULL,
  cb_meta_mon_by_group = NULL,
  cb_meta_dur = NULL,
  cb_meta_dur_by_group = NULL,
  season_k = FALSE
) {
  gids <- as.character(frame_tu$gid)
  U <- length(gids)
  M <- length(months)

  # Covariate column names present in the panel (besides gid, t_start, etc.)
  clim_cols <- c(
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean"
  )
  clim_cols <- intersect(clim_cols, names(panel_tu))
  if (is.null(enso_cols)) {
    enso_cols <- grep("^enso_(lag|lead)[0-9]+$", names(panel_tu), value = TRUE)
  } else {
    enso_cols <- intersect(enso_cols, names(panel_tu))
  }

  # t_start_seas and log_offset must be recomputed from frame_tu$min_seas to
  # avoid off-by-one issues when the panel is sparse (same logic as 09a:1420-1433).
  # log_offset = log(days_in_month(t_start)); t_start_seas per frame row.

  # Index panel rows to (gid, month_idx)
  panel_indexed <- panel_tu %>%
    dplyr::filter(gid %in% gids, t_start %in% months) %>%
    dplyr::mutate(
      unit = as.character(gid),
      time = match(t_start, months)
    ) %>%
    dplyr::filter(!is.na(time))

  # DLNM cross-basis design columns (occurrence model). Reconstructed from the
  # panel's lag columns exactly as at fit time (build_cb_design → make_lag_matrix →
  # dlnm::crossbasis). Empty when cb_meta is NULL (zero behaviour change).
  cb_cols <- character(0)
  if (
    !is.null(cb_meta_by_group) &&
      length(cb_meta_by_group) > 0L &&
      "group" %in% names(frame_tu)
  ) {
    # Group-aware: each unit's cross-basis columns use ITS group's cb_meta.
    # Groups may carry different single_lags / bases, so building one shared
    # design from a single group's cb_meta (the prior behaviour) mis-parses all
    # but that group. Build per-group on the group's panel rows, then stitch
    # back into panel_indexed row order (union of cb_col_names across groups).
    gid2grp <- stats::setNames(
      as.character(frame_tu$group),
      as.character(frame_tu$gid)
    )
    pi_grp <- gid2grp[panel_indexed$unit]
    cb_parts <- list()
    all_cb_names <- character(0)
    for (g in unique(stats::na.omit(pi_grp))) {
      cbm_g <- cb_meta_by_group[[g]]
      if (is.null(cbm_g) || length(cbm_g) == 0L) {
        next
      }
      rows_g <- which(pi_grp == g)
      cb_g <- build_cb_design(panel_indexed[rows_g, , drop = FALSE], cbm_g)
      if (is.null(cb_g) || ncol(cb_g) == 0L) {
        next
      }
      cb_parts[[g]] <- list(rows = rows_g, cb = cb_g)
      all_cb_names <- union(all_cb_names, names(cb_g))
    }
    if (length(all_cb_names) > 0L) {
      cb_full <- as.data.frame(matrix(
        0.0,
        nrow = nrow(panel_indexed),
        ncol = length(all_cb_names)
      ))
      names(cb_full) <- all_cb_names
      for (e in cb_parts) {
        cb_full[e$rows, names(e$cb)] <- e$cb
      }
      cb_cols <- all_cb_names
      panel_indexed <- dplyr::bind_cols(panel_indexed, cb_full)
    }
  } else if (!is.null(cb_meta) && length(cb_meta) > 0L) {
    cb_design <- build_cb_design(panel_indexed, cb_meta)
    if (!is.null(cb_design) && ncol(cb_design) > 0L) {
      cb_cols <- names(cb_design)
      panel_indexed <- dplyr::bind_cols(panel_indexed, cb_design)
    }
  }

  # Non-occurrence PANEL cross-bases (mon, dur), each emitted under its own
  # namespace prefix. `mso` is dropped first: it is latent, has no panel column,
  # and build_cb_design() would fail looking for one (the same reason
  # .lags_required_from_draws() skips .LATENT_CB_VARS). The mon track's mso basis
  # is built once per TU by build_mso_basis() instead.
  #
  # Returns the model-term -> covariate-column map and appends the columns to
  # `panel_indexed` (assigned in the enclosing frame). character(0) when the
  # track is off or carries no cross-basis, which reproduces the pre-prefix
  # behaviour for occurrence-only runs exactly.
  .strip_latent <- function(m) m[setdiff(names(m), .LATENT_CB_VARS)]
  .prefixed_cb <- function(cbm, cbm_by_group, prefix) {
    if (
      !is.null(cbm_by_group) &&
        length(cbm_by_group) > 0L &&
        "group" %in% names(frame_tu)
    ) {
      gid2grp_x <- stats::setNames(
        as.character(frame_tu$group),
        as.character(frame_tu$gid)
      )
      pi_grp_x <- gid2grp_x[panel_indexed$unit]
      parts_x <- list()
      all_x <- character(0)
      for (g in unique(stats::na.omit(pi_grp_x))) {
        cbm_g <- .strip_latent(cbm_by_group[[g]] %||% list())
        if (length(cbm_g) == 0L) {
          next
        }
        rows_g <- which(pi_grp_x == g)
        cb_g <- build_cb_design(panel_indexed[rows_g, , drop = FALSE], cbm_g)
        if (is.null(cb_g) || ncol(cb_g) == 0L) {
          next
        }
        parts_x[[g]] <- list(rows = rows_g, cb = cb_g)
        all_x <- union(all_x, names(cb_g))
      }
      if (length(all_x) == 0L) {
        return(character(0))
      }
      full_x <- as.data.frame(matrix(
        0.0,
        nrow = nrow(panel_indexed),
        ncol = length(all_x)
      ))
      names(full_x) <- all_x
      for (e in parts_x) {
        full_x[e$rows, names(e$cb)] <- e$cb
      }
      names(full_x) <- paste0(prefix, all_x)
      panel_indexed <<- dplyr::bind_cols(panel_indexed, full_x)
      return(stats::setNames(names(full_x), all_x))
    }
    if (!is.null(cbm) && length(.strip_latent(cbm)) > 0L) {
      cb_x <- build_cb_design(panel_indexed, .strip_latent(cbm))
      if (is.null(cb_x) || ncol(cb_x) == 0L) {
        return(character(0))
      }
      map_x <- stats::setNames(paste0(prefix, names(cb_x)), names(cb_x))
      names(cb_x) <- unname(map_x)
      panel_indexed <<- dplyr::bind_cols(panel_indexed, cb_x)
      return(map_x)
    }
    character(0)
  }
  mon_cb_map <- .prefixed_cb(
    cb_meta_mon,
    cb_meta_mon_by_group,
    .SP_MON_CB_PREFIX
  )
  dur_cb_map <- .prefixed_cb(
    cb_meta_dur,
    cb_meta_dur_by_group,
    .SP_DUR_CB_PREFIX
  )
  cb_cols <- c(cb_cols, unname(mon_cb_map), unname(dur_cb_map))

  # Full grid: time 1..M × unit. Fill missing with 0.
  grid <- tidyr::expand_grid(time = seq_len(M), unit = gids)

  # Climate + ENSO + cross-basis columns (per-unit)
  covar <- grid %>%
    dplyr::left_join(
      panel_indexed %>%
        dplyr::select(
          time,
          unit,
          dplyr::all_of(c(clim_cols, enso_cols, cb_cols))
        ),
      by = c("time", "unit")
    ) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(c(clim_cols, enso_cols, cb_cols)),
      ~ dplyr::coalesce(.x, 0.0)
    ))

  # t_start_seas per (unit, month): yday(month shifted by -(min_seas-1) months) - 1.
  # min_seas varies per gid, so the shift MUST be applied per-gid (one scalar
  # `period()` call per min_seas value) -- lubridate::period() given a VECTOR
  # shift argument does NOT apply per-element month shifts; it silently
  # collapses to a single (effectively min_seas=1 / no-shift) result reused for
  # every row, dropping the per-gid cholera-season rotation entirely. Build a
  # per-min_seas lookup table via lapply/vapply instead (mirrors the correct
  # per-gid computation in simulate_one_tu, simulator_core.R:1663-1665) and
  # index into it, rather than a single vectorized `%m+% period(<vector>)` mutate.
  ms_map <- setNames(as.integer(frame_tu$min_seas), as.character(frame_tu$gid))
  ms_vals <- sort(unique(ms_map[is.finite(ms_map)]))
  tss_lut <- vapply(
    ms_vals,
    function(ms) {
      as.integer(
        lubridate::yday(months %m+% lubridate::period(-(ms - 1L), "months")) -
          1L
      )
    },
    integer(M)
  ) # M x length(ms_vals)
  covar <- covar %>%
    dplyr::mutate(
      t_start_raw = months[time],
      min_seas = ms_map[unit],
      t_start_seas = tss_lut[cbind(time, match(min_seas, ms_vals))],
      log_offset = log(as.integer(lubridate::days_in_month(t_start_raw)))
    ) %>%
    dplyr::select(-t_start_raw, -min_seas)

  # Month-of-season index (season onset gate, simulator.onset_eligibility:
  # season): (month - min_seas) mod 12, 0 on the season's first month -- the
  # rotation t_start_seas uses. Read by the C rprocess as sk[u]; it rides the
  # same one-month time shift below, so month k's transition sees month k's
  # index. Emitted only when asked, so every other covariate table is unchanged.
  if (isTRUE(season_k)) {
    covar$seas_k <- as.numeric(
      (as.integer(lubridate::month(months[covar$time])) -
        ms_map[covar$unit]) %%
        12L
    )
  }

  # Monthly-case harmonics. Same construction as load_monthly_data()
  # (utils/mon_data.R) uses at fit time, so the runtime design columns match the
  # fitted coefficients' basis exactly. Note these are CALENDAR-month harmonics
  # (unrotated), unlike t_start_seas above which is rotated by each unit's
  # min_seas — the mon fit's seasonal term is calendar-based.
  if (isTRUE(mon_harmonics)) {
    mnum <- as.numeric(lubridate::month(months[covar$time]))
    covar <- covar %>%
      dplyr::mutate(
        sin_month = sin(2 * pi * mnum / 12),
        cos_month = cos(2 * pi * mnum / 12),
        sin2_month = sin(4 * pi * mnum / 12),
        cos2_month = cos(4 * pi * mnum / 12)
      )
  }

  # Covariate time-alignment (one-month lag fix). The spatPomp object runs under
  # pomp::discrete_time(delta.t = 1) with t0 = 0, so the transition (k-1) -> k
  # evaluates __covars at its LEFT endpoint, time = k-1. To make month-k's hazard
  # use month-k's covariates -- the SAME-MONTH convention the brms model was fit
  # under (make_surv_data joins climate by t_start == month) and the one
  # simulate_one_tu() uses (simulator_core.R reads climate[,m] same-month) --
  # the row carrying calendar-month m's covariates must sit at pomp time m-1.
  # Shift every month one time index earlier (months 1..M -> times 0..M-1) and
  # duplicate the final month at time M so the table still spans t0..M. The
  # previous convention (leading time=0 row duplicating month 1, months at
  # times 1..M) only corrected step 1 and left every later step one month behind
  # (see seam-gap investigation Step V).
  covar <- covar %>% dplyr::mutate(time = time - 1L)

  # Trailing row (time = M): replicate the last month so the covariate table
  # covers observation times 1..M. time = M is not read by any state transition
  # (no step k -> k+1 exists past M); it only keeps the table well-formed.
  tail_rows <- covar %>%
    dplyr::filter(time == M - 1L) %>%
    dplyr::mutate(time = M)

  # .locale = "C" — see the matching comment in build_obs_table(): this
  # ordering is the pomp's unit axis and must agree with every downstream
  # radix-sort re-derivation of it.
  covar_shifted <- dplyr::bind_rows(covar, tail_rows) %>%
    dplyr::arrange(time, unit, .locale = "C") %>%
    dplyr::select(time, unit, dplyr::everything())

  # Model-term -> covariate-column maps for the prefixed (non-occurrence)
  # cross-bases. Carried as attributes rather than changing the return type, so
  # every existing caller is untouched. character(0) whenever the track is off.
  attr(covar_shifted, "mon_cb_map") <- mon_cb_map
  attr(covar_shifted, "dur_cb_map") <- dur_cb_map
  covar_shifted
}

# ---------------------------------------------------------------------------
# Measurement Csnippets — read latent state `p`
# ---------------------------------------------------------------------------
# Probe findings that shape this implementation:
#   (A) rprocess R-function IS honored by bpfilter (scalar per-particle per step).
#   (B) shared_covarnames forbidden alongside per-unit covars — ENSO replicated.
#   (C) R-function dunit/eunit/vunit/runit measures are SILENTLY IGNORED by
#       bpfilter; only Csnippets (mode=3) are dispatched. Fix: carry predictive
#       hazard `p` as a latent state; measures are trivial Csnippets reading `p`.
#   (D) simulate() does not auto-build rmeasure from runit_measure.
#   (E) saved.states extracted via bf@saved.states (exported accessor errors).
#
# dunit output var: `lik`  (+ give_log flag)
# eunit output var: `ey`
# vunit output var: `vc`
# runit output var: `Y`    (NOT the obsname; template-hardcoded)
# State `p` per unit u: __x[__stateindex["p"] + u]  (0-based u inside C)
# Obs  `obs_Y` per unit u: __y[__obsindex["obs_Y"] + u]  (but "obs_Y1" after hitch)

.SP_DUNIT <- pomp::Csnippet(
  "if (ISNA(obs_Y)) { lik = give_log ? 0.0 : 1.0; return; }
   double p_clamped = (p < 1e-9 ? 1e-9 : (p > 1.0-1e-9 ? 1.0-1e-9 : p));
   lik = obs_Y * log(p_clamped) + (1.0 - obs_Y) * log(1.0 - p_clamped);
   if (!give_log) lik = exp(lik);"
)
.SP_EUNIT <- pomp::Csnippet("ey = p;")
.SP_VUNIT <- pomp::Csnippet("vc = p * (1.0 - p);")
.SP_RUNIT <- pomp::Csnippet("Y = (unif_rand() < p) ? 1.0 : 0.0;")

# ---------------------------------------------------------------------------
# State-conditioning measurement Csnippets (siblings).
# ---------------------------------------------------------------------------
# The .SP_DUNIT measure above scores the latent HAZARD `p`, which
# is (nearly) deterministic across particles within a unit-month — so the
# bpfilter cannot discriminate particles by their realized occurrence draw, and
# the filtered `Yout` is a free-running prior simulation that over-fires wherever
# the covariate hazard is high.  The variants below instead treat obs_Y as a
# (near-perfect) measurement of the LATENT occurrence state `Yout`, with a small
# surveillance error EPS, so resampling pins the filtered occurrence to observed
# history.  EPS keeps the weight from collapsing to exactly zero when no particle
# matches (avoids degenerate ESS); it is the false-positive/false-negative rate.
#
# ONSET TERM x FAMILY TEMPLATE
# ----------------------------
# Every state-conditioning dunit is (NA guard) + (ONSET term on Yout vs obs_Y) +
# (that family's CASE term) + (lik tail). The onset term is the ONLY thing that
# differs between an eps sibling and its detection-error ("_detect") sibling, so
# it is written exactly once per scheme (.sp_onset_eps_c / .sp_onset_detect_c)
# and each family's case term exactly once (.SP_TPL_*). The eps objects below
# are built from those templates at module load; the detect siblings are built
# by the SAME factories at dispatch time inside build_spatpomp_tu_c(), because
# delta_fn/delta_fp are per-run config parameters and a statically built
# Csnippet cannot vary per run.
#
# Contract for template authors: the case term must NOT reference `match` (the
# detect onset block does not declare it) — gate on `Yout`/`obs_Y`/`mu` instead.
.SP_MEAS_EPS <- 1e-6
# Detection-error rates, used only by the `_detect` siblings. These are DEFAULTS
# now, not fixed constants: simulator.filter.delta_fn / .delta_fp override them
# per run (analysis/config_defaults.yml). Provisional and uncalibrated.
.SP_STATE_DELTA_FN <- 0.10 # P(obs_Y=0 | Yout=1): surveillance miss rate
.SP_STATE_DELTA_FP <- 0.05 # P(obs_Y=1 | Yout=0): surveillance false-alarm rate
# eps of the `_actn` inactive-state case mixture: P(a count is
# reported | unit inactive). 1.0 = the plain floor NB (the panel estimate is
# 0 by construction). Default of
# simulator.filter.case_fp_rate and simulator.filter.who_obs.fp_rate.
.SP_CASE_FP_RATE <- 1.0

#' Onset term, eps scheme: near-delta match between latent Yout and obs_Y.
.sp_onset_eps_c <- function(eps = .SP_MEAS_EPS) {
  sprintf(
    "    double eps = %.1e;
    int match = ((Yout > 0.5) == (obs_Y > 0.5));
    double ll = match ? log(1.0 - eps) : log(eps);",
    eps
  )
}

#' Onset term, detection-error scheme: an explicit two-parameter imperfect-
#' detection measurement model (a standard sensitivity/specificity confusion-
#' matrix likelihood), decoupled from the fitted occurrence hazard `p`:
#'   P(obs_Y=1 | Yout=1) = 1 - delta_fn   (sensitivity)
#'   P(obs_Y=1 | Yout=0) = delta_fp       (false-alarm rate)
#'
#' Why this exists: the eps scheme's single tiny epsilon (~1e-6) has to be that
#' extreme only because it is silently overriding the PROCESS-level mismatch
#' rate implied by the fitted hazard itself (Yout ~ Bernoulli(p); for rare
#' cholera onset, p is small, so a particle's own draw disagrees with the one
#' historical realization at a rate of roughly 1-p on onset months -- large, and
#' nothing to do with surveillance quality). A single eps conflates that purely
#' aleatory process disagreement with genuine surveillance detection error and
#' ends up doing almost none of the latter's nominal job -- which is why
#' recovery under eps looks artificially perfect. delta_fn/delta_fp are meant as
#' real (if provisional) surveillance error rates, NOT tuned for numerical
#' convenience, so recovery under this scheme is EXPECTED to degrade gracefully
#' as `p` shrinks rather than staying near-exact.
.sp_onset_detect_c <- function(
  delta_fn = .SP_STATE_DELTA_FN,
  delta_fp = .SP_STATE_DELTA_FP
) {
  .check_deltas(delta_fn, delta_fp, where = ".sp_onset_detect_c")
  sprintf(
    "    double delta_fn = %.4f;
    double delta_fp = %.4f;
    double ll;
    if (obs_Y > 0.5) {
      ll = (Yout > 0.5) ? log(1.0 - delta_fn) : log(delta_fp);
    } else {
      ll = (Yout > 0.5) ? log(delta_fn) : log(1.0 - delta_fp);
    }",
    delta_fn,
    delta_fp
  )
}

# Family template: occurrence only, no case term.
.SP_TPL_STATE <- "if (ISNA(obs_Y)) { lik = give_log ? 0.0 : 1.0; return; }
%s
    lik = give_log ? ll : exp(ll);"
.make_dunit_state <- function(onset) {
  pomp::Csnippet(sprintf(.SP_TPL_STATE, onset))
}
.SP_DUNIT_STATE <- .make_dunit_state(.sp_onset_eps_c())

# ---------------------------------------------------------------------------
# Monthly-case (`mon`) measurement Csnippets.
# ---------------------------------------------------------------------------
# Selected when measurement is in .SP_MON_MEAS. Occurrence term as .SP_DUNIT /
# .SP_DUNIT_STATE, plus a NegBin case term gated on `mu > 0.0` that scores
# `obs_CM` (this month's observed case count). `mu` carries the mon fit's
# per-month mean and `size_phi` the mon fit's dispersion.
#
# Gating on `mu > 0.0` rather than the RECORDED `active` state is deliberate,
# not cosmetic. `mon_mu_block` (spatpomp_core_c.R) computes mu from the PRE-cessation
# av_new/mv_new -- correctly nonzero for an episode's final active month --
# but the cessation block that runs immediately after can reset active/mso to
# 0 in that SAME step before writing them into ST(9,u)/ST(10,u). So the
# RECORDED `active` state for a month reflects "entering next month", not
# "was active this month", and gating on `active > 0.5` silently drops the
# case observation for the final month of every single episode. `mu` has no
# such ambiguity: it is set to exactly 0.0 when inactive and to exp(lp) > 0.0
# otherwise (mon_mu_block's `if (av_new[u] < 0.5) { ST(7,u) = 0.0; ... }`), so
# `mu > 0.0` is an exact, order-independent proxy for "active this month".
#
# Hurdle structure per (unit, month):
#   obs_Y NA                     -> loglik 0    (no observation)
#   obs_Y == 0                   -> log(1 - p)
#   obs_Y == 1                   -> log(p)
#   + mu > 0, obs_CM observed    -> + dnbinom(obs_CM; mu, size_phi)
#   + mu > 0, obs_CM NA          -> + 0         (active but magnitude unobserved)
# The occurrence half is untouched, so with every obs_CM NA this degenerates
# EXACTLY to .SP_DUNIT.
.SP_DUNIT_MON_NB <- pomp::Csnippet(
  "if (ISNA(obs_Y)) { lik = give_log ? 0.0 : 1.0; return; }
   double pc = (p < 1e-9 ? 1e-9 : (p > 1.0-1e-9 ? 1.0-1e-9 : p));
   double ll = (obs_Y > 0.5) ? log(pc) : log(1.0 - pc);
   if (mu > 0.0 && !ISNA(obs_CM)) {
     double r = (size_phi <= 0.0 ? 1e-6 : size_phi);
     double m = (mu < 1e-9 ? 1e-9 : mu);
     ll += lgamma(obs_CM + r) - lgamma(r) - lgamma(obs_CM + 1.0)
         + r * (log(r) - log(r + m)) + obs_CM * (log(m) - log(r + m));
   }
   lik = give_log ? ll : exp(ll);"
)
# State-conditioning sibling: occurrence half is the near-delta on latent `Yout`
# (.SP_DUNIT_STATE), monthly case half is unchanged from .SP_DUNIT_MON_NB. The
# case term is NOT nested inside the `match` branch: `mu` can legitimately be
# > 0 in a month where `Yout` is 0, so
# gating the monthly term on `match` would drop every non-onset month of every
# episode -- i.e. almost all of the monthly signal.
# ---------------------------------------------------------------------------
# mon-family dunit combinator: (NA guard) + (ONSET term) + (CASE term)
# ---------------------------------------------------------------------------
# Every mon-family dunit is the same three pieces. The onset terms already have
# factories (.sp_onset_eps_c / .sp_onset_detect_c); these are the case-term
# factories, and .make_dunit_mon() splices them. The pre-existing factories
# below are thin wrappers over it and emit byte-identical C (pinned by
# test_spatpomp_dunit_detect.R's golden).
#
# score_case_if_na_y = FALSE (every pre-existing sibling): the whole density is
#   skipped when obs_Y is NA -- the classic `if (ISNA(obs_Y)) return` guard.
# score_case_if_na_y = TRUE (the `_act` sibling): only the ONSET term is skipped
#   when obs_Y is NA; the CASE term is still scored whenever obs_CM is observed.
#   Needed because the survival panel censors every month after an onset, so
#   obs_Y is NA in exactly the post-onset months whose counts carry the
#   cessation information.
#' The negative-binomial log-pmf as a C expression, log NB(x; mean m, size r),
#' with `r` and `m` C variables in scope. `cont` is the indentation of the
#' continuation line, so every case-term factory emits the same text it always
#' did (pinned by test_spatpomp_dunit_detect.R's golden).
.sp_lnb_c <- function(x, cont) {
  sprintf(
    "lgamma(%s + r) - lgamma(r) - lgamma(%s + 1.0)\n%s+ r * (log(r) - log(r + m)) + %s * (log(m) - log(r + m))",
    x,
    x,
    strrep(" ", cont),
    x
  )
}
.sp_case_mon_nb_c <- function() {
  paste0(
    "    if (mu > 0.0 && !ISNA(obs_CM)) {\n      double r = (size_phi <= 0.0 ? 1e-6 : size_phi);\n      double m = (mu < 1e-9 ? 1e-9 : mu);\n      ll += ",
    .sp_lnb_c("obs_CM", 10L),
    ";\n    }"
  )
}
.sp_case_actpen_c <- function(mu_floor = .SP_MON_MU_FLOOR) {
  sprintf(
    "    if (!ISNA(obs_CM)) {\n      double r = (size_phi <= 0.0 ? 1e-6 : size_phi);\n      double mu_floor = %.2f;\n      double m_eff = (mu > 0.0) ? mu : mu_floor;\n      double m = (m_eff < 1e-9 ? 1e-9 : m_eff);\n      ll += %s;\n    }",
    mu_floor,
    .sp_lnb_c("obs_CM", 10L)
  )
}

#' Zero-inflated activity case term (the `_act` sibling). With delta_fn the
#' probability that a running episode goes unrecorded:
#'   obs_CM NA                -> 0
#'   mu > 0,  obs_CM == 0     -> log(delta_fn + (1 - delta_fn) * NB(0; mu, phi))
#'   mu > 0,  obs_CM >  0     -> log(1 - delta_fn) + log NB(obs_CM; mu, phi)
#'   mu == 0, obs_CM >  0     -> log NB(obs_CM; mu_floor, phi)   (ACTPEN floor)
#'   mu == 0, obs_CM == 0     -> 0
#' The zero branch is what charges a latent episode that never ceases: each
#' empty observed month costs log(delta_fn + (1-delta_fn) NB(0)), bounded below
#' by log(delta_fn), so a long false run accumulates a graded penalty instead of
#' the zero cost it had before.
.sp_case_actpen_zi_c <- function(
  delta_fn = .SP_STATE_DELTA_FN,
  mu_floor = .SP_MON_MU_FLOOR
) {
  .check_deltas(delta_fn, 0.5, where = ".sp_case_actpen_zi_c")
  sprintf(
    "    if (!ISNA(obs_CM)) {
      double r = (size_phi <= 0.0 ? 1e-6 : size_phi);
      double dfn = %.6f;
      if (mu > 0.0) {
        double m = (mu < 1e-9 ? 1e-9 : mu);
        double lnb = %s;
        ll += (obs_CM <= 0.0) ? log(dfn + (1.0 - dfn) * exp(lnb))
                              : log(1.0 - dfn) + lnb;
      } else if (obs_CM > 0.0) {
        double m = %.2f;
        ll += %s;
      }
    }",
    delta_fn,
    .sp_lnb_c("obs_CM", 19L),
    mu_floor,
    .sp_lnb_c("obs_CM", 12L)
  )
}

#' Normalized zero-inflated activity case term (the `_actn` sibling).
#' Active branches are .sp_case_actpen_zi_c()'s, verbatim. The inactive branch
#' is a mixture -- no count reported with probability 1 - eps, else
#' NB(obs_CM; mu_0, phi) -- so its total mass is 1 for every phi:
#'   mu == 0, obs_CM == 0     -> log((1 - eps) + eps * NB(0; mu_0, phi))
#'   mu == 0, obs_CM >  0     -> log(eps) + log NB(obs_CM; mu_0, phi)
#' (the `_act` pair 0 / log NB(obs_CM; mu_floor, phi) has mass 2 - NB(0)).
#' eps = 1 is the plain floor NB at every count.
.sp_case_actpen_zi_norm_c <- function(
  delta_fn = .SP_STATE_DELTA_FN,
  eps = .SP_CASE_FP_RATE,
  mu_0 = .SP_MON_MU_FLOOR
) {
  .check_deltas(delta_fn, 0.5, where = ".sp_case_actpen_zi_norm_c")
  .check_case_fp(eps, mu_0, where = ".sp_case_actpen_zi_norm_c")
  sprintf(
    "    if (!ISNA(obs_CM)) {
      double r = (size_phi <= 0.0 ? 1e-6 : size_phi);
      double dfn = %.6f;
      if (mu > 0.0) {
        double m = (mu < 1e-9 ? 1e-9 : mu);
        double lnb = %s;
        ll += (obs_CM <= 0.0) ? log(dfn + (1.0 - dfn) * exp(lnb))
                              : log(1.0 - dfn) + lnb;
      } else {
        double eps = %.17g;
        double m = %.17g;
        double lnb = %s;
        ll += (obs_CM <= 0.0) ? log((1.0 - eps) + eps * exp(lnb))
                              : log(eps) + lnb;
      }
    }",
    delta_fn,
    .sp_lnb_c("obs_CM", 19L),
    eps,
    mu_0,
    .sp_lnb_c("obs_CM", 19L)
  )
}

#' WHO dashboard country term (the `_who` measurement). Scores the national
#' monthly count `obs_W` against the particle's SIMULATED country case count
#' `ccm` (the country sum of the per-unit draws `cm`, built in rprocess),
#' scaled to the national total by `scale / who_f` (who_f = population share
#' of the country covered by the simulated gids):
#'   obs_W NA        -> 0
#'   ccm > 0         -> w_u * log NB(obs_W; ccm * scale / who_f, phi)
#'   ccm == 0        -> w_u * log NB(obs_W; mu_floor, phi)   (graded, not -Inf)
#' `phi` (config simulator.filter.who_obs.phi) is the NB size linking the
#' simulated national count to the WHO data: reporting overdispersion only,
#' process noise is already in ccm. `w_u` (who_w) is the unit's share of the
#' term; the shares sum to 1 over the country, so with country blocks the term
#' enters its block exactly once. Per-unit constants are embedded as static
#' arrays in spatPomp's unit order (unit u = 0-based index of the radix-sorted
#' gids), so no covariate offset is involved. obs_W is replicated on every
#' unit of the country; DB-episode cases also enter per gid through obs_CM
#' (documented double use).
#'
#' @param who_w,who_s per-unit share and per-unit multiplier scale / who_f,
#'   both in spatPomp unit order.
.sp_case_who_c <- function(who_w, who_s, phi, mu_floor) {
  stopifnot(
    length(who_w) == length(who_s),
    all(is.finite(who_w)),
    all(is.finite(who_s)),
    phi > 0,
    mu_floor > 0
  )
  num <- function(v) paste(sprintf("%.17g", v), collapse = ", ")
  sprintf(
    "    if (!ISNA(obs_W)) {
      static const double WHO_W[%d] = { %s };
      static const double WHO_S[%d] = { %s };
      double r = %.17g;
      double m = (ccm > 0.0) ? ccm * WHO_S[u] : %.17g;
      if (m < 1e-9) m = 1e-9;
      ll += WHO_W[u] * (%s);
    }",
    length(who_w),
    num(who_w),
    length(who_s),
    num(who_s),
    phi,
    mu_floor,
    .sp_lnb_c("obs_W", 10L)
  )
}

#' Normalized WHO country term (the `_actn_who` measurement).
#' .sp_case_who_c() with the country-inactive branch (ccm == 0) as the same
#' mixture as the unit term, eps = `fp_rate`, mu_0 = `mu_floor`:
#'   ccm == 0, obs_W == 0 -> w_u * log((1 - eps) + eps * NB(0; mu_floor, phi))
#'   ccm == 0, obs_W >  0 -> w_u * (log(eps) + log NB(obs_W; mu_floor, phi))
#' The `_who` floor branch is already a plain NB (mass 1); fp_rate = 1 scores
#' exactly what it scores.
.sp_case_who_norm_c <- function(
  who_w,
  who_s,
  phi,
  mu_floor,
  fp_rate = .SP_CASE_FP_RATE
) {
  stopifnot(
    length(who_w) == length(who_s),
    all(is.finite(who_w)),
    all(is.finite(who_s)),
    phi > 0,
    mu_floor > 0
  )
  .check_case_fp(fp_rate, mu_floor, where = ".sp_case_who_norm_c")
  num <- function(v) paste(sprintf("%.17g", v), collapse = ", ")
  sprintf(
    "    if (!ISNA(obs_W)) {
      static const double WHO_W[%d] = { %s };
      static const double WHO_S[%d] = { %s };
      double r = %.17g;
      if (ccm > 0.0) {
        double m = ccm * WHO_S[u];
        if (m < 1e-9) m = 1e-9;
        ll += WHO_W[u] * (%s);
      } else {
        double eps = %.17g;
        double m = %.17g;
        if (m < 1e-9) m = 1e-9;
        double lnb = %s;
        ll += WHO_W[u] * ((obs_W <= 0.0) ? log((1.0 - eps) + eps * exp(lnb))
                                         : log(eps) + lnb);
      }
    }",
    length(who_w),
    num(who_w),
    length(who_s),
    num(who_s),
    phi,
    .sp_lnb_c("obs_W", 12L),
    fp_rate,
    mu_floor,
    .sp_lnb_c("obs_W", 19L)
  )
}

.make_dunit_mon <- function(
  onset,
  case,
  score_case_if_na_y = FALSE,
  extra = NULL
) {
  # `extra`: an additional independent term (the WHO country term), scored like
  # the case term whenever its own observation is present. Only the `_act`
  # form carries it; NULL emits exactly the previous C.
  if (!is.null(extra) && !isTRUE(score_case_if_na_y)) {
    stop(".make_dunit_mon: `extra` needs score_case_if_na_y = TRUE.")
  }
  if (!isTRUE(score_case_if_na_y)) {
    return(pomp::Csnippet(paste0(
      "if (ISNA(obs_Y)) { lik = give_log ? 0.0 : 1.0; return; }\n",
      onset,
      "\n",
      case,
      "\n    lik = give_log ? ll : exp(ll);"
    )))
  }
  pomp::Csnippet(paste0(
    "double ll_t = 0.0;\n",
    "if (!ISNA(obs_Y)) {\n",
    onset,
    "\n    ll_t += ll;\n}\n",
    "{\n    double ll = 0.0;\n",
    case,
    "\n    ll_t += ll;\n}\n",
    if (!is.null(extra)) {
      paste0("{\n    double ll = 0.0;\n", extra, "\n    ll_t += ll;\n}\n")
    },
    "lik = give_log ? ll_t : exp(ll_t);"
  ))
}

.make_dunit_mon_nb_state <- function(onset) {
  .make_dunit_mon(onset, .sp_case_mon_nb_c())
}
# No `_detect` sibling is built from this template: _ACTPEN's floor-based case
# term (below) strictly dominates this one, so the mon track's detect sibling is
# built on _ACTPEN instead.
.SP_DUNIT_MON_NB_STATE <- .make_dunit_mon_nb_state(.sp_onset_eps_c())
# Activity-penalized sibling ("actpen" -- Idea B from the mon-recovery
# investigation, replacing an earlier Idea-A boolean-match attempt): the
# problem being fixed is the same as documented for .SP_DUNIT_MON_NB_STATE
# above -- the NegBin case-count term is silently SKIPPED whenever a
# particle's own mu == 0, so a particle that wrongly declared cessation pays
# nothing for missing a real observed case count. Root-caused empirically on
# the BBB_2_1 fixture unit: a 12-month true episode collapsed to 0% of
# filtered passes still "active" by month 5, because nothing in the
# likelihood disfavored early spurious cessation.
#
# Idea A (tried first, reverted here) added a second near-deterministic
# eps/1-eps state-conditioning term on activity match, on top of the onset
# match. That made the activity penalty ~13.8 nats on mismatch -- 5-8x
# larger than the NegBin term's entire dynamic range in this fixture -- so it
# behaved as a near-hard constraint rather than a graded probabilistic
# penalty (flagged as "too strong").
#
# Idea B (this version) instead keeps the NegBin term itself as the only
# source of penalty, but stops SKIPPING it when mu == 0: it substitutes a
# small floor mean (`.SP_MON_MU_FLOOR`) for likelihood evaluation ONLY when
# a case count was actually observed and the particle's own mu is 0. This
# grades the penalty by how implausible the observed count is under a
# near-zero mean, instead of an all-or-nothing eps switch, and does not
# touch the actual mu state written to ST(7,u) (downstream smoothing/
# forecast code relies on that being exactly 0 when inactive).
#
# Coverage note (unchanged from Idea A's caveat): this only fixes the
# "falsely ceased, real count observed" direction. It does nothing for the
# opposite direction (falsely still active, nothing observed) -- there is no
# observation to grade against in that case, so it remains unpenalized here,
# same as in .SP_DUNIT_MON_NB_STATE.
.SP_MON_MU_FLOOR <- 0.1
.make_dunit_mon_nb_state_actpen <- function(
  onset,
  mu_floor = .SP_MON_MU_FLOOR
) {
  .make_dunit_mon(onset, .sp_case_actpen_c(mu_floor))
}
.SP_DUNIT_MON_NB_STATE_ACTPEN <- .make_dunit_mon_nb_state_actpen(.sp_onset_eps_c())
# The mon track's detection-error sibling is this same _ACTPEN template with
# .sp_onset_detect_c() in the onset slot (term C, the mu-floor NegBin magnitude
# term, is reused verbatim -- the detect scheme only ever changes term A). It is
# NOT a static object: delta_fn/delta_fp are per-run config parameters, so it is
# built at dispatch time in build_spatpomp_tu_c() via
# .make_dunit_mon_nb_state_actpen(.sp_onset_detect_c(delta_fn, delta_fp)).
# See .sp_onset_detect_c() above for the eps-vs-fitted-`p` rationale.
#' The ONE measurement -> dunit Csnippet dispatch used by build_spatpomp_tu_c().
#' Detect siblings are built here, at dispatch time, because delta_fn/delta_fp
#' are per-run config values; the eps siblings stay module-level objects.
.select_dunit <- function(
  measurement,
  delta_fn = .SP_STATE_DELTA_FN,
  delta_fp = .SP_STATE_DELTA_FP,
  who_c = NULL,
  case_fp_rate = .SP_CASE_FP_RATE,
  case_fp_mean = .SP_MON_MU_FLOOR
) {
  if (.is_who_meas(measurement) && is.null(who_c)) {
    stop(
      ".select_dunit: measurement '",
      measurement,
      "' needs the WHO term ",
      "(who_c, from .sp_case_who_c())."
    )
  }
  onset_c <- if (.is_detect_meas(measurement)) {
    .sp_onset_detect_c(delta_fn, delta_fp)
  } else {
    NULL
  }
  if (.is_mon_meas(measurement)) {
    if (.is_actn_meas(measurement)) {
      .make_dunit_mon(
        onset_c,
        .sp_case_actpen_zi_norm_c(delta_fn, case_fp_rate, case_fp_mean),
        score_case_if_na_y = TRUE,
        extra = if (.is_who_meas(measurement)) who_c
      )
    } else if (.is_act_meas(measurement)) {
      .make_dunit_mon(
        onset_c,
        .sp_case_actpen_zi_c(delta_fn),
        score_case_if_na_y = TRUE,
        extra = if (.is_who_meas(measurement)) who_c
      )
    } else if (.is_detect_meas(measurement)) {
      .make_dunit_mon_nb_state_actpen(onset_c)
    } else if (.is_actpen_meas(measurement)) {
      .SP_DUNIT_MON_NB_STATE_ACTPEN
    } else if (.is_state_meas(measurement)) {
      .SP_DUNIT_MON_NB_STATE
    } else {
      .SP_DUNIT_MON_NB
    }
  } else {
    if (.is_detect_meas(measurement)) {
      .make_dunit_state(onset_c)
    } else if (.is_state_meas(measurement)) {
      .SP_DUNIT_STATE
    } else {
      .SP_DUNIT
    }
  }
}

# ---------------------------------------------------------------------------
# Shared per-group draw-subset helpers (used by build_spatpomp_tu_c())
# ---------------------------------------------------------------------------

#' Subset each group's occurrence draws bundle for grouped (multi-group-per-TU)
#' filtering.
#'
#' Mirrors the single-bundle subsetting in build_spatpomp_tu_c: under the
#' fixed_draw_idx (N-independent-runs) path each group is reduced to a single
#' draw (the index clamped to that group's S); otherwise each group is
#' subsampled to min(n_draws_in_filter, S). Guards that every group present in
#' frame_tu$group has a bundle in grouping$draws_occ_by_group.
#'
#' @param grouping          grouping list carrying $draws_occ_by_group (named by
#'                          group label).
#' @param fixed_draw_idx    integer or NULL — same semantics as build_spatpomp_tu_c.
#' @param n_draws_in_filter integer — number of draws to subsample per group when
#'                          fixed_draw_idx is NULL.
#' @param frame_tu          location frame (used only for the membership guard);
#'                          NULL skips the guard.
#' @return list(draws_sub_by_group = <named list of 1-or-n-draw bundles>).
.build_grouping_sub <- function(
  grouping,
  fixed_draw_idx,
  n_draws_in_filter,
  frame_tu = NULL
) {
  if (is.null(grouping)) {
    return(NULL)
  }

  if (!is.null(frame_tu) && !is.null(frame_tu$group)) {
    needed <- unique(as.character(frame_tu$group))
    have <- names(grouping$draws_occ_by_group)
    missing <- setdiff(needed, have)
    if (length(missing) > 0L) {
      stop(
        ".build_grouping_sub: no occurrence draws bundle for group(s) ",
        paste(missing, collapse = ", "),
        " present in frame_tu$group (have: ",
        paste(have, collapse = ", "),
        ")"
      )
    }
  }

  draws_sub_by_group <- lapply(grouping$draws_occ_by_group, function(d_g) {
    if (!is.null(fixed_draw_idx)) {
      d_idx_g <- as.integer(min(as.integer(fixed_draw_idx), d_g$S))
      .subset_draws(d_g, d_idx_g)
    } else {
      n_g <- min(as.integer(n_draws_in_filter), d_g$S)
      sub_g <- if (n_g >= d_g$S) {
        seq_len(d_g$S)
      } else {
        sort(sample.int(d_g$S, n_g))
      }
      .subset_draws(d_g, sub_g)
    }
  })
  list(draws_sub_by_group = draws_sub_by_group)
}

#' Monthly-case-draws sibling of .build_grouping_sub.
#'
#' Subsets each group's MONTHLY (`mon`) draws bundle
#' (`grouping$draws_mon_by_group`) using the same fixed_draw_idx / n-draw
#' subsample strategy as the occurrence/duration bundles, so the monthly LP and
#' compute_phi_vec use the matching per-group draw.
#'
#' @return list(draws_sub_by_group = <named list of 1-or-n-draw monthly bundles>).
.build_grouping_mon_sub <- function(
  grouping,
  fixed_draw_idx,
  n_draws_in_filter,
  frame_tu = NULL
) {
  if (is.null(grouping)) {
    return(NULL)
  }
  mng <- grouping$draws_mon_by_group
  if (is.null(mng)) {
    stop(
      ".build_grouping_mon_sub: grouped monthly-case measurement requires ",
      "grouping$draws_mon_by_group."
    )
  }

  if (!is.null(frame_tu) && !is.null(frame_tu$group)) {
    needed <- unique(as.character(frame_tu$group))
    missing <- setdiff(needed, names(mng))
    if (length(missing) > 0L) {
      stop(
        ".build_grouping_mon_sub: no monthly draws bundle for group(s) ",
        paste(missing, collapse = ", "),
        " present in frame_tu$group (have: ",
        paste(names(mng), collapse = ", "),
        ")"
      )
    }
  }

  draws_sub_by_group <- lapply(mng, function(d_g) {
    if (!is.null(fixed_draw_idx)) {
      d_idx_g <- as.integer(min(as.integer(fixed_draw_idx), d_g$S))
      .subset_draws(d_g, d_idx_g)
    } else {
      n_g <- min(as.integer(n_draws_in_filter), d_g$S)
      sub_g <- if (n_g >= d_g$S) {
        seq_len(d_g$S)
      } else {
        sort(sample.int(d_g$S, n_g))
      }
      .subset_draws(d_g, sub_g)
    }
  })
  list(draws_sub_by_group = draws_sub_by_group)
}

#' Duration-draws sibling of .build_grouping_sub.
#'
#' Subsets each group's DURATION (cessation-hazard) draws bundle
#' (`grouping$draws_dur_by_group`) using the same fixed_draw_idx / n-draw
#' subsample strategy as the occurrence bundles, so the
#' cessation-hazard evaluation uses the matching per-group draw. Mirrors
#' `simulate_one_tu()`'s own `grouping$draws_dur_by_group` lookup
#' (`simulator_core.R:2226`), just subsetted for the filter window instead of
#' used verbatim.
#'
#' @return list(draws_sub_by_group = <named list of 1-or-n-draw duration bundles>).
.build_grouping_dur_sub <- function(
  grouping,
  fixed_draw_idx,
  n_draws_in_filter,
  frame_tu = NULL
) {
  if (is.null(grouping)) {
    return(NULL)
  }
  drg <- grouping$draws_dur_by_group
  if (is.null(drg)) {
    stop(
      ".build_grouping_dur_sub: grouped duration-aware filtering requires ",
      "grouping$draws_dur_by_group."
    )
  }

  if (!is.null(frame_tu) && !is.null(frame_tu$group)) {
    needed <- unique(as.character(frame_tu$group))
    missing <- setdiff(needed, names(drg))
    if (length(missing) > 0L) {
      stop(
        ".build_grouping_dur_sub: no duration draws bundle for group(s) ",
        paste(missing, collapse = ", "),
        " present in frame_tu$group (have: ",
        paste(names(drg), collapse = ", "),
        ")"
      )
    }
  }

  draws_sub_by_group <- lapply(drg, function(d_g) {
    if (!is.null(fixed_draw_idx)) {
      d_idx_g <- as.integer(min(as.integer(fixed_draw_idx), d_g$S))
      .subset_draws(d_g, d_idx_g)
    } else {
      n_g <- min(as.integer(n_draws_in_filter), d_g$S)
      sub_g <- if (n_g >= d_g$S) {
        seq_len(d_g$S)
      } else {
        sort(sample.int(d_g$S, n_g))
      }
      .subset_draws(d_g, sub_g)
    }
  })
  list(draws_sub_by_group = draws_sub_by_group)
}

# ---------------------------------------------------------------------------
# Builder blocks shared by build_spatpomp_tu() and build_spatpomp_tu_c()
# ---------------------------------------------------------------------------

#' .sp_filter_window_months
#'
#' The panel months inside the filter window [window$start, window$end) -- the
#' observed data times 1..M. Stops (naming `caller` and the TU) when empty.
#' @noRd
.sp_filter_window_months <- function(
  panel_tu_warmup,
  window,
  frame_tu,
  caller
) {
  wstart <- as.Date(window$start)
  wend <- as.Date(window$end)
  all_months <- sort(unique(panel_tu_warmup$t_start))
  filt_months <- all_months[all_months >= wstart & all_months < wend]
  if (length(filt_months) == 0L) {
    stop(
      caller,
      ": no panel months within the filter window ",
      window$start,
      " – ",
      window$end,
      " for TU ",
      paste(unique(frame_tu$tu), collapse = "/")
    )
  }
  filt_months
}

#' .sp_enso_covar_names
#'
#' ENSO covariate names: the panel's enso_lag* columns that the occurrence
#' draws actually use. Under grouping, the union across groups, so every
#' group's required lag column is in the covariate table. Restricting to model
#' lags avoids paste0(colname, unit) collisions inside
#' spatPomp::construct_spatPomp. Ungrouped with no model lag: the first panel
#' lag column (if any).
#' @noRd
.sp_enso_covar_names <- function(panel_tu_warmup, draws_occ, grouping) {
  all_panel_enso <- grep(
    "^enso_(lag|lead)[0-9]+$",
    names(panel_tu_warmup),
    value = TRUE
  )
  if (!is.null(grouping)) {
    fe_enso <- unique(unlist(lapply(grouping$draws_occ_by_group, function(d) {
      intersect(colnames(d$fe), all_panel_enso)
    })))
  } else {
    fe_enso <- intersect(colnames(draws_occ$fe), all_panel_enso)
  }
  enso_names <- if (length(fe_enso) > 0L) {
    fe_enso
  } else if (is.null(grouping) && length(all_panel_enso) > 0L) {
    all_panel_enso[1L]
  } else {
    character(0)
  }
  if (length(enso_names) == 0L) {
    enso_names <- character(0)
  }
  enso_names
}

#' .sp_cb_covar_names
#'
#' DLNM cross-basis covariate column names of the occurrence model (union
#' across groups under grouping; empty without a cross-basis).
#' @noRd
.sp_cb_covar_names <- function(cb_meta, grouping) {
  cb_covar_names <- if (!is.null(grouping)) {
    unique(unlist(lapply(grouping$cb_meta_occ_by_group, function(m) {
      if (!is.null(m) && length(m) > 0L) {
        unlist(lapply(m, function(x) x$cb_col_names), use.names = FALSE)
      } else {
        character(0)
      }
    })))
  } else if (!is.null(cb_meta) && length(cb_meta) > 0L) {
    unlist(lapply(cb_meta, function(m) m$cb_col_names), use.names = FALSE)
  } else {
    character(0)
  }
  cb_covar_names
}

# ---------------------------------------------------------------------------
# Internal: subset a draws bundle to a selected index vector
# ---------------------------------------------------------------------------

.subset_draws <- function(draws, idx) {
  sub <- draws
  sub$S <- length(idx)
  sub$fe <- draws$fe[idx, , drop = FALSE]

  sub_mat <- function(m) if (is.null(m)) NULL else m[idx, , drop = FALSE]
  sub$re_intercept <- lapply(draws$re_intercept, sub_mat)

  sub$re_slope <- lapply(draws$re_slope, function(tier) {
    if (is.null(tier)) {
      return(NULL)
    }
    if (is.matrix(tier)) {
      return(sub_mat(tier))
    }
    lapply(tier, sub_mat)
  })

  if (!is.null(draws$phi)) {
    sub$phi <- if (is.list(draws$phi)) {
      lapply(draws$phi, function(x) if (is.numeric(x)) x[idx] else x)
    } else {
      draws$phi[idx]
    }
  }

  if (!is.null(draws$smooth_coef)) {
    sub$smooth_coef <- lapply(draws$smooth_coef, function(m) {
      m[idx, , drop = FALSE]
    })
  }

  sub
}

# ---------------------------------------------------------------------------
# Block-list builder (exposed API)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Partition -> frame-column mapping (shared by the forward simulator 07b and
# the historical-smoothing phase 08). Single source of truth for how a
# filter.partition_by switch selects the frame column that keys each partition.
# ---------------------------------------------------------------------------

#' Map a filter.partition_by switch to the frame column that keys each
#' partition. "transmission_unit" keys on the HMM cluster already carried as
#' `tu`; "country" on ISO3; "enso_region" on the ENSO `group` label. Errors on
#' an unknown value.
partition_frame_column <- function(partition_by) {
  switch(
    partition_by,
    country = "country",
    enso_region = "group",
    registry = "part",
    transmission_unit = "tu",
    stop(
      sprintf("unknown filter.partition_by: '%s' ", partition_by),
      "(expected one of transmission_unit, country, enso_region, registry)."
    )
  )
}

#' load_partition_map
#'
#' ISO3 -> partition (TU) label from a country-level registry
#' (analysis/groupings/<source>.csv, utils/config.R::load_grouping()), with the
#' labels asserted path-safe (assert_partition_labels()). Named character.
load_partition_map <- function(source, groupings_dir = "analysis/groupings") {
  if (is.null(source) || !nzchar(source)) {
    stop(
      "filter.partition_by: registry requires simulator.filter.partition_source ",
      "(a registry name under analysis/groupings/).",
      call. = FALSE
    )
  }
  reg <- load_grouping(source, groupings_dir)
  assert_partition_labels(reg$group, source)
  stats::setNames(reg$group, reg$iso3)
}

#' stamp_partition_column
#'
#' For partition_by = "registry": write `frame$part` = each row's partition
#' label from the `partition_source` registry. A country present on the frame
#' but absent from the registry is a hard error (the policy of the fit-grouping
#' loader, smoothing_core.R). Other partition_by values: frame returned as is.
stamp_partition_column <- function(
  frame,
  partition_by,
  partition_source = NULL,
  groupings_dir = "analysis/groupings"
) {
  if (!identical(partition_by, "registry")) {
    return(frame)
  }
  pmap <- load_partition_map(partition_source, groupings_dir)
  ctry <- as.character(frame$country)
  miss <- sort(unique(stats::na.omit(ctry[!ctry %in% names(pmap)])))
  if (length(miss)) {
    stop(
      sprintf(
        paste0(
          "partition registry '%s' does not cover %d country(ies) on the ",
          "frame: %s. Add them to analysis/groupings/%s.csv or narrow ",
          "scope.countries."
        ),
        partition_source,
        length(miss),
        paste(miss, collapse = ", "),
        partition_source
      ),
      call. = FALSE
    )
  }
  frame$part <- unname(pmap[ctry])
  frame
}

#' Resolve the partition column for `partition_by` and assert it is present and
#' not all-NA on `frame`. Returns the column name. The enso_region case appends
#' the "set simulator.grouping.enable: true" hint (that key is what populates
#' the `group` column).
assert_partition_column <- function(frame, partition_by) {
  col <- partition_frame_column(partition_by)
  if (is.null(frame[[col]]) || all(is.na(frame[[col]]))) {
    stop(
      sprintf(
        "filter.partition_by=%s requires a `%s` column on the frame",
        partition_by,
        col
      ),
      if (identical(partition_by, "enso_region")) {
        " — set simulator.grouping.enable: true."
      } else if (identical(partition_by, "registry")) {
        " — set simulator.filter.partition_source (stamp_partition_column())."
      } else {
        "."
      }
    )
  }
  col
}

#' check_partition_vs_grouping
#'
#' Pre-flight for a partition whose TUs may mix fit groups (partition_by =
#' registry; also run for enso_region / country). Collects every violation and
#' stops ONCE with all of them, before any per-TU precompute:
#'   - each country on the frame carries a fit group (grouped mode);
#'   - within a mixed TU, the groups' duration fits share one basehaz kind and
#'     (factor) breaks -- the C duration kernel holds one table per TU
#'     (precompute_dur_lp_components());
#'   - the neighbour-contagion partition REFINES the TU partition: every
#'     contagion block that has a unit in a TU lies entirely inside that TU.
#'     The sim-time kernel is cut to the TU's ADM1s while keeping the fit-time
#'     group-wide row-normalisation (build_partition_nbr_kernel()), so a block
#'     spanning TUs silently loses its out-of-TU neighbours. An error under
#'     `registry`; a warning under the older country / transmission_unit keys,
#'     where the truncation already existed and is left as a known limitation.
#' Informational: the per-TU draw ceiling is min S over the TU's groups
#' (07b samples draw indices below it).
#'
#' @param frame        location frame with `country`, `tu` (already re-keyed to
#'                     the partition) and, in grouped mode, `group`.
#' @param partition_by filter.partition_by.
#' @param grouping     grouped routing list (draws_*_by_group) or NULL.
#' @param nbr_spec     resolve_nbr_spec() result or NULL.
#' @return invisible list(tu_groups = named list TU -> fit groups,
#'         draw_ceiling = named integer TU -> min S).
check_partition_vs_grouping <- function(
  frame,
  partition_by,
  grouping = NULL,
  nbr_spec = NULL
) {
  errs <- character(0)
  warns <- character(0)
  tu <- as.character(frame$tu)
  tu[is.na(tu)] <- "_NA_TU_"
  ctry <- as.character(frame$country)
  grouped <- !is.null(grouping) && !is.null(frame$group)
  tu_groups <- list()
  draw_ceiling <- integer(0)

  if (grouped) {
    nog <- sort(unique(ctry[is.na(frame$group)]))
    if (length(nog)) {
      errs <- c(
        errs,
        sprintf(
          "country(ies) without a fit group: %s",
          paste(nog, collapse = ", ")
        )
      )
    }
    grp <- as.character(frame$group)
    tu_groups <- lapply(split(grp, tu), function(x) {
      sort(unique(stats::na.omit(x)))
    })
    S_of <- vapply(
      grouping$draws_occ_by_group %||% list(),
      function(d) as.integer(d$S %||% NA_integer_),
      integer(1)
    )
    draw_ceiling <- vapply(
      tu_groups,
      function(gs) {
        if (length(gs) && all(gs %in% names(S_of))) {
          min(S_of[gs])
        } else {
          NA_integer_
        }
      },
      integer(1)
    )
    # Duration basehaz coherence inside each mixed TU.
    dd <- grouping$draws_dur_by_group
    if (
      !is.null(dd) && exists("precompute_dur_basehaz_table", mode = "function")
    ) {
      bh <- lapply(dd, function(d) {
        tryCatch(precompute_dur_basehaz_table(d), error = function(e) NULL)
      })
      for (t in names(tu_groups)) {
        gs <- intersect(tu_groups[[t]], names(bh))
        if (length(gs) < 2L) {
          next
        }
        kinds <- vapply(
          bh[gs],
          function(b) as.character(b$kind %||% NA),
          character(1)
        )
        if (length(unique(kinds)) > 1L) {
          errs <- c(
            errs,
            sprintf(
              "TU '%s' mixes duration basehaz kinds (%s)",
              t,
              paste(sprintf("%s=%s", gs, kinds), collapse = ", ")
            )
          )
        } else if (identical(kinds[[1L]], "factor")) {
          b1 <- bh[[gs[1L]]]$breaks
          same <- vapply(
            bh[gs],
            function(b) isTRUE(all.equal(b$breaks, b1)),
            logical(1)
          )
          if (!all(same)) {
            errs <- c(
              errs,
              sprintf(
                "TU '%s' mixes duration basehaz breaks (groups %s)",
                t,
                paste(gs, collapse = ", ")
              )
            )
          }
        }
      }
    }
  }

  # Neighbour-contagion partition must refine the TU partition.
  if (isTRUE(nbr_spec$use) && !is.null(nbr_spec$registry)) {
    blk_of <- stats::setNames(nbr_spec$registry$group, nbr_spec$registry$iso3)
    blk <- unname(blk_of[ctry])
    span <- vapply(split(tu, blk), function(x) length(unique(x)), integer(1))
    bad <- names(span)[span > 1L]
    if (length(bad)) {
      msg <- sprintf(
        paste0(
          "neighbour-contagion block(s) %s span several TUs (%s): the ",
          "sim-time kernel would drop their out-of-TU neighbours. Use a ",
          "partition registry that the contagion registry refines."
        ),
        paste(bad, collapse = ", "),
        paste(
          vapply(
            bad,
            function(b) {
              sprintf(
                "%s -> %s",
                b,
                paste(sort(unique(tu[blk %in% b])), collapse = "/")
              )
            },
            character(1)
          ),
          collapse = "; "
        )
      )
      if (identical(partition_by, "registry")) {
        errs <- c(errs, msg)
      } else {
        warns <- c(warns, msg)
      }
    }
  }

  for (w in warns) {
    warning("[partition pre-flight] ", w, call. = FALSE)
  }
  if (length(errs)) {
    stop(
      sprintf(
        "[partition pre-flight] filter.partition_by=%s:\n  - %s",
        partition_by,
        paste(errs, collapse = "\n  - ")
      ),
      call. = FALSE
    )
  }
  mixed <- names(tu_groups)[lengths(tu_groups) > 1L]
  if (length(mixed)) {
    message(sprintf(
      "[partition pre-flight] %d TU(s) mix fit groups: %s",
      length(mixed),
      paste(
        sprintf(
          "%s{%s; S<=%s}",
          mixed,
          vapply(tu_groups[mixed], paste, character(1), collapse = ","),
          draw_ceiling[mixed]
        ),
        collapse = " "
      )
    ))
  }
  invisible(list(tu_groups = tu_groups, draw_ceiling = draw_ceiling))
}

#' Build block_list for bpfilter from a spatPomp object's unit ordering.
#'
#' @param obj         spatPomp object (to extract unit_names slot).
#' @param frame_tu    location frame for this TU.
#' @param granularity block grouping key: "country" (default; one block per
#'                    country) or "adm1" (one block per ADM1 unit). ADM1 blocks
#'                    are finer and reduce the number of simultaneous outbreaks
#'                    a single block's particles must jointly match, mitigating
#'                    particle-weight degeneracy under near-delta measurements.
#' @return list of integer vectors; one block per grouping key.
make_block_list <- function(obj, frame_tu, granularity = c("country", "adm1")) {
  granularity <- match.arg(granularity)
  unit_names <- obj@unit_names
  gids <- as.character(frame_tu$gid)
  key_col <- if (granularity == "adm1") "adm1_unit" else "country"
  keys <- frame_tu[[key_col]][match(unit_names, gids)]
  key_ord <- unique(keys)

  blks <- lapply(stats::setNames(key_ord, key_ord), function(k) {
    which(keys == k)
  })
  # Drop any empty blocks (can occur if unit_names has entries not in frame_tu).
  blks[vapply(blks, length, integer(1)) > 0L]
}

# ---------------------------------------------------------------------------
# Run bpfilter for one TU
# ---------------------------------------------------------------------------

#' Run bpfilter on a per-TU spatPomp object.
#'
#' @param obj        spatPomp object from build_spatpomp_tu_c().
#' @param Np         integer; number of particles per block.
#' @param block_list list of integer vectors; one block per country (from
#'                   make_block_list()). NULL triggers automatic per-country
#'                   derivation from obj's metadata.
#' @param seed       integer RNG seed for reproducibility.
#' @param ...        additional args forwarded to spatPomp::bpfilter().
#' @return bpfilterd_spatPomp object.
run_bpfilter_tu <- function(obj, Np, block_list = NULL, seed = NULL, ...) {
  frame_tu <- attr(obj, "frame_tu")
  if (is.null(block_list)) {
    block_list <- make_block_list(obj, frame_tu)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  spatPomp::bpfilter(
    obj,
    Np = Np,
    block_list = block_list,
    save_states = TRUE,
    ...
  )
}

# ---------------------------------------------------------------------------
# Extract filtered states → warmup handoff list
# ---------------------------------------------------------------------------

#' .state_vec_to_warmup_row
#'
#' One warm-start mini-panel (the build_warmup_state() row shape) from one
#' particle's named state vector: lag states, the o_tu sentinel back to NA, the
#' u_days -> time-to-last-outbreak bin (NA -> "[6,Inf]", the reference-level
#' maximum), the posterior draw index (same for all units), the duration-aware
#' active/mso slots and, on the detected feed, `det` (absent column == latent
#' run). Shared by filtered_states_to_warmup() (live bpf, per particle) and
#' filter_traj_to_warmup() (smoothing trajectory).
#'
#' @param vec named state vector (one particle / one time)
#' @param pos state position of each gid
#' @param gids frame gids, in output row order
#' @param warmup_date t_start of the rows (one month before the horizon)
#' @param require_active_mso stop when the active/mso slots are absent (an
#'   older-schema smoothing cache) instead of returning NA
#' @noRd
.state_vec_to_warmup_row <- function(
  vec,
  pos,
  gids,
  warmup_date,
  require_active_mso = FALSE
) {
  u_days_j <- vec[paste0("u_days", pos)]
  o_j <- vec[paste0("o", pos)]
  o_ct_j <- vec[paste0("o_ct", pos)]
  o_tu_raw_j <- vec[paste0("o_tu", pos)]
  draw_j <- vec[paste0("draw_idx", pos[1L])] # same for all units
  active_names <- paste0("active", pos)
  mso_names <- paste0("mso", pos)
  if (
    isTRUE(require_active_mso) &&
      (!all(active_names %in% names(vec)) || !all(mso_names %in% names(vec)))
  ) {
    stop(
      "filter_traj_to_warmup: state vector has no 'active'/'mso' entries -- ",
      "this trajectory predates the duration-aware state schema ",
      "(.SP_UNIT_STATENAMES 9->11); if it came from a smoothing cache, rebuild ",
      "it (08a/08b/08c_smoothing_*.R)."
    )
  }
  active_j <- vec[active_names]
  mso_j <- vec[mso_names]
  det_names <- paste0(.SP_DET_STATENAME, pos)
  det_j <- if (all(det_names %in% names(vec))) vec[det_names] else NULL

  o_tu_j <- ifelse(
    o_tu_raw_j <= .SP_OTU_SENTINEL + 0.5,
    NA_integer_,
    as.integer(round(o_tu_raw_j))
  )
  # has_none: a sentinel state (TTLO level "none") hands off as
  # "none"; it exists only when the occurrence fit carries the level.
  ttlo_j <- bin_ttlo(u_days_j, has_none = TRUE)
  ttlo_j[is.na(ttlo_j)] <- "[6,Inf]"

  tibble::tibble(
    gid = gids,
    t_start = warmup_date,
    ongoing_outbreak_lag1 = as.integer(round(o_j)),
    ongoing_outbreak_country_lag1 = as.integer(round(o_ct_j)),
    ongoing_outbreak_tu_lag1 = o_tu_j,
    time_to_last_outbreak_cat = ttlo_j,
    posterior_draw_idx = as.integer(round(draw_j)),
    active = as.integer(round(active_j)),
    mso = as.integer(round(mso_j)),
    det = if (is.null(det_j)) {
      NULL
    } else {
      as.integer(round(det_j))
    }
  )
}

#' Convert final-time bpfilter particle states to a list of mini-panel tibbles
#' compatible with simulate_one_tu(..., warmup_panels_list=).
#'
#' The last element of bf@saved.states is a (7U) × Np matrix. Each column j is
#' one particle's final state. We extract the lag fields and draw_idx to build
#' the per-replicate warmup mini-panel that build_warmup_state() expects.
#'
#' @param bpf            bpfilterd_spatPomp object.
#' @param frame_tu       location frame for this TU.
#' @param horizon_start  Date; the forecast horizon start (= filter window end).
#'                       Mini-panel rows get t_start = month before this date.
#' @return list of Np tibbles, each with columns (gid, t_start,
#'         ongoing_outbreak_lag1, ongoing_outbreak_country_lag1,
#'         ongoing_outbreak_tu_lag1, time_to_last_outbreak_cat,
#'         posterior_draw_idx, active, mso).
filtered_states_to_warmup <- function(bpf, frame_tu, horizon_start) {
  # Extract final saved-state slice via slot (exported accessor errors in v1.1.0)
  saved <- bpf@saved.states
  if (length(saved) == 0L) {
    stop(
      "filtered_states_to_warmup: bpf@saved.states is empty — ",
      "run with save_states=TRUE."
    )
  }
  final_mat <- saved[[length(saved)]] # (7U) × Np matrix
  Np <- ncol(final_mat)
  U <- nrow(frame_tu)
  gids <- as.character(frame_tu$gid)
  # State position of each frame gid, by NAME: bpf@unit_names is the pomp's
  # sorted unit order, which need not match frame_tu's row order.
  pos <- match(gids, bpf@unit_names)
  if (anyNA(pos)) {
    stop(
      "filtered_states_to_warmup: frame gid(s) absent from the filtered object: ",
      paste(gids[is.na(pos)], collapse = ", ")
    )
  }

  # t_start for mini-panel rows: one month before horizon_start.
  warmup_date <- as.Date(horizon_start) %m-% months(1L)

  # Duration-aware state (active/mso): always present in a live bpf's
  # state vector (.SP_UNIT_STATENAMES 9->11), held at 0/0 on the legacy path,
  # so no schema check here (unlike filter_traj_to_warmup()'s cache input).
  lapply(seq_len(Np), function(j) {
    .state_vec_to_warmup_row(final_mat[, j], pos, gids, warmup_date)
  })
}

# ---------------------------------------------------------------------------
# N-independent-runs BPF design — new functions (M4 redesign)
# ---------------------------------------------------------------------------

#' Run bpfilter once with filter_traj = TRUE for one replicate.
#'
#' This is the per-replicate workhorse for the N-independent-runs design.
#' Call build_spatpomp_tu_c with a fixed_draw_idx BEFORE calling this function;
#' the resulting object has all particles pinned to one occurrence posterior draw.
#'
#' @param obj        spatPomp object from build_spatpomp_tu_c(fixed_draw_idx = d).
#' @param Np         integer; number of particles.
#' @param block_list list of integer vectors (from make_block_list()).
#' @param seed       integer RNG seed (optional).
#' @param traj_storage `simulator.filter.traj_storage`: "full" (spatPomp's own
#'                   bpfilter, every particle of every step kept), "pruned" or
#'                   "pruned_float32" (.bpfilter_traj_compact()). See
#'                   .SP_TRAJ_STORAGES.
#' @param prune_every `simulator.filter.traj_prune_every`: ancestry-pruning
#'                   period in filter steps (pruned modes only).
#' @return list(filter_traj = array(state_vars × 1 × M), logLik = numeric)
run_one_bpfilter_replicate <- function(
  obj,
  Np,
  block_list,
  seed = NULL,
  traj_storage = "pruned",
  prune_every = 1L
) {
  traj_storage <- .check_traj_storage(traj_storage)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  bpf <- if (identical(traj_storage, "full")) {
    spatPomp::bpfilter(
      obj,
      Np = Np,
      block_list = block_list,
      filter_traj = TRUE
    )
  } else {
    .bpfilter_traj_compact(
      obj,
      Np = Np,
      block_list = lapply(block_list, as.integer),
      float32 = identical(traj_storage, "pruned_float32"),
      prune_every = prune_every,
      verbose = getOption("spatPomp_verbose", FALSE)
    )
  }
  list(
    filter_traj = bpf@filter.traj, # (7U × 1 × M) array
    logLik = pomp::logLik(bpf)
  )
}

# ---------------------------------------------------------------------------
# Compact filter_traj storage (simulator.filter.traj_storage)
# ---------------------------------------------------------------------------
# spatPomp 1.1.0's bpfilter(filter_traj = TRUE) keeps every resampled particle
# of every step and block (`xparticles`, doubles) plus every ancestry vector
# (`pedigree`) until the end, then draws ONE particle per block and walks its
# ancestry back. That store, ntimes x state vars x Np x 8 bytes, is the whole
# memory peak of a smoothing replicate (deployed hmmCE unit: ~90 GB at Np = 10k
# over ~213 months), and almost all of it is particles with no descendant in
# the final population, which the backward walk can never reach.
#
# .bpfilter_traj_compact() is spatPomp:::bpfilter.internal (1.1.0) with two
# storage changes and nothing else:
#   (a) ancestry pruning -- every `prune_every` steps, per block (ancestry is
#       per block), drop the stored columns that are no longer an ancestor of
#       the current population and remap the pedigree indices onto the kept
#       columns. The current top step is never pruned, so the final draw
#       sample.int(Np, 1) indexes the same population as the original;
#   (b) float32 (optional) -- kept states stored as float::fl(), converted
#       back to double one column at a time in the backward walk. Integer
#       states are exact below 2^24; p / mu round to ~6e-8 relative.
# Pruning consumes no random numbers and the backward walk reaches the same
# particles, so "pruned" returns a filter.traj and logLik identical() to
# "full" for the same seed (tests/test_bpfilter_traj_compact.R).
#
# Its environment is set to the spatPomp namespace (below), so the shared
# lines resolve the same internals (pStop_, vec_dmeasure, the
# bpfilter_computations C routine) as the original. Re-diff it against
# deparse(spatPomp:::bpfilter.internal) on any spatPomp upgrade: the test
# checks the installed version is the one this copy was taken from.
.SP_TRAJ_STORAGES <- c("full", "pruned", "pruned_float32")
.SP_TRAJ_COMPACT_SPATPOMP_VERSION <- "1.1.0"

.check_traj_storage <- function(traj_storage) {
  traj_storage <- as.character(traj_storage)
  if (length(traj_storage) != 1L || !traj_storage %in% .SP_TRAJ_STORAGES) {
    stop(
      "simulator.filter.traj_storage must be one of ",
      paste(.SP_TRAJ_STORAGES, collapse = " | "),
      ", got '",
      paste(traj_storage, collapse = ","),
      "'.",
      call. = FALSE
    )
  }
  if (
    !identical(traj_storage, "full") &&
      as.character(utils::packageVersion("spatPomp")) !=
        .SP_TRAJ_COMPACT_SPATPOMP_VERSION
  ) {
    stop(
      "simulator.filter.traj_storage '",
      traj_storage,
      "' runs a copy of ",
      "spatPomp ",
      .SP_TRAJ_COMPACT_SPATPOMP_VERSION,
      "'s bpfilter.internal; ",
      "spatPomp ",
      as.character(utils::packageVersion("spatPomp")),
      " is installed. Re-diff .bpfilter_traj_compact() first.",
      call. = FALSE
    )
  }
  if (
    identical(traj_storage, "pruned_float32") &&
      !requireNamespace("float", quietly = TRUE)
  ) {
    stop(
      "simulator.filter.traj_storage 'pruned_float32' needs the 'float' ",
      "package: install.packages(\"float\").",
      call. = FALSE
    )
  }
  traj_storage
}

#' bpfilter(filter_traj = TRUE, save_states = FALSE) with pruned (and
#' optionally float32) particle storage. See the section header above.
#'
#' @param object      spatPomp object.
#' @param Np          particles per block (length 1 or ntimes + 1).
#' @param block_list  list of integer unit-index vectors.
#' @param float32     store kept states as float32 (float package).
#' @param prune_every prune the stored ancestry every this many steps.
#' @param diag_env    optional environment; receives `kept_cols`, the stored
#'                    columns per level (rows, 1 = init) and block (columns)
#'                    at the end of the run (tests / benchmarks).
#' @return bpfilterd_spatPomp, as spatPomp::bpfilter() (saved.states empty).
#' @noRd
.bpfilter_traj_compact <- function(
  object,
  Np,
  block_list,
  float32 = FALSE,
  prune_every = 1L,
  diag_env = NULL,
  ...,
  verbose,
  .gnsi = TRUE
) {
  ep <- paste0("in ", sQuote("bpfilter"), ": ")
  verbose <- as.logical(verbose)
  prune_every <- as.integer(prune_every)
  if (length(prune_every) != 1L || is.na(prune_every) || prune_every < 1L) {
    pStop_(ep, sQuote("prune_every"), " must be a positive integer")
  }
  p_object <- pomp(object, ...)
  object <- new(
    "spatPomp",
    p_object,
    unit_covarnames = object@unit_covarnames,
    shared_covarnames = object@shared_covarnames,
    runit_measure = object@runit_measure,
    dunit_measure = object@dunit_measure,
    eunit_measure = object@eunit_measure,
    munit_measure = object@munit_measure,
    vunit_measure = object@vunit_measure,
    unit_names = object@unit_names,
    unitname = object@unitname,
    unit_statenames = object@unit_statenames,
    unit_obsnames = object@unit_obsnames,
    unit_accumvars = object@unit_accumvars
  )
  pompLoad(object, verbose = verbose)
  on.exit(pompUnload(object, verbose = verbose))
  gnsi <- as.logical(.gnsi)
  times <- time(object, t0 = TRUE)
  ntimes <- length(times) - 1
  nblocks <- length(block_list)
  if (length(Np) == 1) {
    Np <- rep(Np, times = ntimes + 1)
  } else if (length(Np) != (ntimes + 1)) {
    pStop_(ep, sQuote("Np"), " must have length 1 or length ", ntimes + 1)
  }
  if (any(Np <= 0)) {
    pStop_(
      ep,
      "number of particles, ",
      sQuote("Np"),
      ", must always be positive"
    )
  }
  if (!is.numeric(Np)) {
    pStop_(
      ep,
      sQuote("Np"),
      " must be a number, a vector of numbers, or a function"
    )
  }
  Np <- as.integer(Np)
  params <- as.matrix(coef(object))
  init.x <- rinit(object, params = params, nsim = Np[1L], .gnsi = gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x
  weights <- array(data = numeric(0), dim = c(nblocks, Np[1L]))
  loglik <- rep(NA, ntimes)
  block.loglik <- matrix(NA, nblocks, ntimes)
  us <- object@unit_statenames
  block_sn <- lapply(block_list, function(block) {
    paste0(rep(us, length(block)), rep(block, each = length(us)))
  })

  # --- compact store (replaces xparticles / pedigree) ------------------------
  # Per block i: st[[i]][[lev]] = kept states at level lev (lev 1 = init.x,
  # lev nt + 1 = after step nt's resampling), one column per stored particle;
  # anc[[i]][[lev]] (lev >= 2) = for each stored column of level lev, its
  # parent COLUMN in st[[i]][[lev - 1]]; alive[[i]][[lev]] = which stored
  # columns are still ancestors of the top population (NULL = all of them).
  # A level is physically compacted only once its live columns are at most
  # half of its stored ones, so each level is copied O(1) times in total
  # instead of at every prune (re-copying the recent, still-shrinking levels at
  # every step cost ~30 % run time in GC on the test fixture) and holds at
  # most twice its live particles. Compaction remaps the parent indices of the
  # level above; a dead column's parent may remap to 0, which is never read
  # (the backward walk only visits live columns; tabulate() ignores 0) and
  # stays 0 at later compactions.
  # Invariant after a prune at top level T: for every lev <= T the live set of
  # lev - 1 is exactly the parent set of the live set of lev, so a later prune
  # can stop at the first level at or below T whose live set does not shrink.
  # Under float32 the stored matrix is the float32 bit pattern
  # (float::fl(m)@Data, a plain integer matrix), so pruning subsets it with
  # base `[` (float's own `[` method fails on one-column extraction from a
  # matrix with rownames); only the backward walk decodes, a column at a time.
  .enc <- if (float32) {
    function(m) {
      r <- suppressWarnings(range(m, finite = TRUE))
      if (all(is.finite(r)) && max(abs(r)) > 3.4e38) {
        pStop_(
          ep,
          "a state exceeds the float32 range; use traj_storage 'pruned'"
        )
      }
      float::fl(m)@Data
    }
  } else {
    identity
  }
  .dec <- if (float32) function(v) float::dbl(float::float32(v)) else identity
  st <- lapply(block_sn, function(sn) {
    l <- vector("list", ntimes + 1L)
    l[[1L]] <- .enc(init.x[sn, , drop = FALSE])
    l
  })
  anc <- lapply(block_sn, function(sn) vector("list", ntimes + 1L))
  alive <- lapply(block_sn, function(sn) vector("list", ntimes + 1L))
  pruned_top <- integer(nblocks) # top level of the last prune, per block
  .prune <- function(i, top) {
    s <- st[[i]]
    a <- anc[[i]]
    al <- alive[[i]]
    lev <- top
    while (lev >= 2L) {
      par <- if (is.null(al[[lev]])) a[[lev]] else a[[lev]][al[[lev]]]
      n_below <- NCOL(s[[lev - 1L]])
      hit <- tabulate(par, nbins = n_below) > 0L
      n_hit <- sum(hit)
      n_old <- if (is.null(al[[lev - 1L]])) n_below else length(al[[lev - 1L]])
      # hit is a subset of the old live set, so equal counts = equal sets.
      if (n_hit == n_old && lev - 1L <= pruned_top[i]) {
        break
      } # nothing below can change
      if (2L * n_hit <= n_below) {
        keep <- which(hit)
        remap <- integer(n_below)
        remap[keep] <- seq_along(keep)
        # c(0L, remap)[. + 1L], not remap[.]: a dead column's parent may
        # already be 0 from an earlier compaction, and R drops a 0 index,
        # which would shorten a[[lev]] and shift every later parent.
        a[[lev]] <- c(0L, remap)[a[[lev]] + 1L]
        s[[lev - 1L]] <- s[[lev - 1L]][, keep, drop = FALSE]
        if (lev - 1L >= 2L) {
          a[[lev - 1L]] <- a[[lev - 1L]][keep]
        }
        al[lev - 1L] <- list(NULL)
      } else {
        al[lev - 1L] <- list(if (n_hit == n_below) NULL else which(hit))
      }
      lev <- lev - 1L
    }
    st[[i]] <<- s
    anc[[i]] <<- a
    alive[[i]] <<- al
    pruned_top[i] <<- top
    invisible(NULL)
  }

  filt.t <- array(
    data = numeric(1),
    dim = c(nvars, 1, ntimes + 1),
    dimnames = list(name = statenames, rep = 1, time = NULL)
  )
  for (nt in seq_len(ntimes)) {
    max_log_d <- vector(mode = "numeric", length = nblocks)
    X <- tryCatch(
      rprocess(
        object,
        x0 = x,
        t0 = times[nt],
        times = times[nt + 1],
        params = params,
        .gnsi = gnsi
      ),
      error = function(e) {
        pStop_(ep, "process simulation error: ", conditionMessage(e))
      }
    )
    for (i in seq(nblocks)) {
      block <- block_list[[i]]
      log_vd <- tryCatch(
        vec_dmeasure(
          object,
          y = object@data[, nt, drop = FALSE],
          x = X,
          units = block,
          times = times[nt + 1],
          params = params,
          log = TRUE,
          .gnsi = gnsi
        ),
        error = function(e) {
          pStop_(ep, "error in calculation of weights: ", conditionMessage(e))
        }
      )
      log_d <- apply(log_vd[,, 1, drop = FALSE], 2, function(x) sum(x))
      max_log_d[i] <- max(log_d)
      if (is.infinite(max_log_d[i])) {
        weights[i, ] <- rep(1, Np[nt + 1])
      } else {
        log_d <- log_d - max_log_d[i]
        weights[i, ] <- exp(log_d)
      }
    }
    gnsi <- FALSE
    for (i in seq_len(nblocks)) {
      statenames <- block_sn[[i]]
      tempX <- X[statenames, , , drop = FALSE]
      xx <- tryCatch(
        .Call(
          "bpfilter_computations",
          x = tempX,
          params = params,
          Np = Np[nt + 1],
          trackancestry = TRUE,
          doparRS = FALSE,
          weights = weights[i, ]
        ),
        error = function(e) {
          stop(ep, conditionMessage(e), call. = FALSE)
        }
      )
      x[statenames, ] <- xx$states
      params <- xx$params
      st[[i]][[nt + 1L]] <- .enc(xx$states)
      anc[[i]][[nt + 1L]] <- xx$ancestry
      if (nt %% prune_every == 0L) .prune(i, nt + 1L)
    }
    rm(X, xx, tempX)
    log_weights <- matrix(0, nrow = nblocks, ncol = Np[nt + 1])
    for (i in seq_len(nblocks)) {
      if (is.infinite(max_log_d[i])) {
        block.loglik[i, nt] <- -Inf
      } else {
        log_weights[i, ] <- max_log_d[i] + log(weights[i, ])
        block.loglik[i, nt] <- logmeanexp(log_weights[i, ])
      }
    }
    loglik[nt] <- sum(block.loglik[, nt])
  }
  if (!is.null(diag_env)) {
    assign(
      "kept_cols",
      vapply(st, function(s) vapply(s, NCOL, integer(1)), integer(ntimes + 1L)),
      envir = diag_env
    )
  }
  # Backward walk: the original's, over the compact store. Level ntimes + 1 is
  # never pruned (it is always the top), so b indexes the full population.
  for (i in seq(nblocks)) {
    block_statenames <- block_sn[[i]]
    s <- st[[i]]
    a <- anc[[i]]
    b <- sample.int(n = ncol(weights), size = 1L, replace = TRUE)
    filt.t[block_statenames, 1L, ntimes + 1] <- .dec(s[[ntimes + 1L]][, b])
    for (lev in seq.int(from = ntimes + 1L, to = 3L, by = -1L)) {
      b <- a[[lev]][b]
      filt.t[block_statenames, 1L, lev - 1L] <- .dec(s[[lev - 1L]][, b])
    }
    if (times[2L] > times[1L]) {
      b <- a[[2L]][b]
      filt.t[block_statenames, 1L, 1L] <- .dec(s[[1L]][, b])
    }
  }
  if (times[2L] <= times[1L]) {
    filt.t <- filt.t[,, -1L, drop = FALSE]
  }
  new(
    "bpfilterd_spatPomp",
    object,
    block_list = block_list,
    Np = as.integer(Np),
    cond.loglik = loglik,
    block.cond.loglik = block.loglik,
    loglik = sum(loglik),
    saved.states = list(),
    filter.traj = filt.t
  )
}
environment(.bpfilter_traj_compact) <- asNamespace("spatPomp")

#' Extract the final-step smoothing state from a filter.traj array.
#'
#' Converts the result of run_one_bpfilter_replicate into a single mini-panel
#' tibble compatible with simulate_one_tu(..., warmup_panels_list=).
#' This is the single-trajectory counterpart of filtered_states_to_warmup.
#'
#' @param ft_array      3-D array (state_vars × 1 × M) from bpf@filter.traj.
#' @param frame_tu      location frame for this TU.
#' @param horizon_start Date; forecast horizon start (= filter window end).
#' @return single tibble with columns gid, t_start, ongoing_outbreak_lag1,
#'         ongoing_outbreak_country_lag1, ongoing_outbreak_tu_lag1,
#'         time_to_last_outbreak_cat, posterior_draw_idx, active, mso.
filter_traj_to_warmup <- function(
  ft_array,
  frame_tu,
  horizon_start,
  smoothed_gids = NULL,
  filt_months = NULL
) {
  U <- nrow(frame_tu)
  gids <- as.character(frame_tu$gid)
  warmup_date <- as.Date(horizon_start) %m-% months(1L)

  # Default: the IC is the last smoothed step (end-of-coverage handoff). For an
  # interior-month (hindcast) handoff, filt_months lets us slice the smoothed
  # state at the requested launch month (= warmup_date) rather than the final
  # column. In the end-of-coverage case warmup_date == max(filt_months), so this
  # resolves to the same last column — backward-compatible.
  M_idx <- dim(ft_array)[3]
  if (!is.null(filt_months) && length(filt_months)) {
    m <- match(warmup_date, as.Date(filt_months))
    if (is.na(m)) {
      stop(sprintf(
        "filter_traj_to_warmup: IC month %s not found in filt_months [%s, %s] — cannot slice interior handoff.",
        format(warmup_date),
        format(min(as.Date(filt_months), na.rm = TRUE)),
        format(max(as.Date(filt_months), na.rm = TRUE))
      ))
    }
    # filt_months[m] lives at filter_traj column m + toff, where toff skips the
    # leading t0 column(s): a pomp filter.traj has ntimes+1 columns (col 1 = t0,
    # cols 2..ntimes+1 = the M filter months). Mirrors .smoothing_active_matrix()
    # (smoothing_core.R): toff <- dim(ft)[3] - M. For end-of-coverage (m == M)
    # this resolves to dim(ft_array)[3], matching the filt_months = NULL branch.
    M_idx <- m + (dim(ft_array)[3] - length(filt_months))
  }

  # Position map: forecast unit j → smoothed spatPomp position pos[j].
  # When smoothed_gids is provided the forecast frame may be a strict subset
  # of the smoothed units (e.g. unseen_adm1 restrict drops gids that were
  # present in the smoothing). A NULL smoothed_gids falls back to the original
  # 1:U identity map (exact-size case, backward-compatible).
  if (is.null(smoothed_gids)) {
    pos <- seq_len(U)
  } else {
    pos <- match(gids, smoothed_gids)
    missing_g <- gids[is.na(pos)]
    if (length(missing_g) > 0L) {
      stop(
        "filter_traj_to_warmup: forecast gid(s) absent from smoothed unit list: ",
        paste(missing_g, collapse = ", ")
      )
    }
  }

  # Unlike filtered_states_to_warmup() (a live bpf object built in the current
  # session), ft_array here may come straight from an on-disk smoothing cache
  # written by an older schema -- plain
  # named-vector indexing would then return silent NA instead of erroring, so
  # the presence of the active/mso names is checked explicitly.
  .state_vec_to_warmup_row(
    ft_array[, 1, M_idx],
    pos,
    gids,
    warmup_date,
    require_active_mso = TRUE
  )
}

#' Stack several per-partition filter_traj arrays into one region-wide array.
#'
#' Each member array is a pomp `filter.traj` of shape `[length(statenames)*U_k ×
#' 1 × M]` whose first-dim names are `<statename><position>` for local positions
#' `1..U_k` (see `filter_traj_to_warmup()`, which reads state purely by these
#' names). This concatenates the members along the UNIT axis, re-indexing each
#' member's position suffixes to a contiguous `1..U_region` block so the merged
#' array's names resolve exactly like a single-partition array. The month (3rd)
#' axis is kept intact and MUST agree across members (they share one smoothing
#' run's `filt_months`, so the interior-month handoff slice in
#' `filter_traj_to_warmup()` still lands on the right column).
#'
#' Assignment is name-based, so the merged array is correct irrespective of the
#' statename/position row ordering spatPomp happens to emit.
#'
#' @param ft_list    list of per-member 3-D arrays (bpf@filter.traj shape).
#' @param gids_list  list of per-member gid vectors; `gids_list[[k]]` has length
#'                   `U_k` and maps member position → gid (region unit order is
#'                   the concatenation across members, in list order).
#' @param statenames per-unit state names (default `.SP_UNIT_STATENAMES`).
#' @return list(filter_traj = region array `[length(statenames)*U_region × 1 ×
#'         M]`, gids = concatenated region gid vector).
merge_filter_traj_by_unit <- function(ft_list, gids_list, statenames = NULL) {
  stopifnot(length(ft_list) == length(gids_list), length(ft_list) >= 1L)
  # NULL: read the layout off the trajectories themselves (11 slots, or 12
  # under covariate_feed = "detected" -- see .sp_unit_statenames()).
  statenames <- statenames %||%
    .sp_unit_statenames(
      if (.sp_has_det(dimnames(ft_list[[1L]])[[1L]])) "detected" else "latent"
    )
  Ms <- vapply(ft_list, function(a) dim(a)[3], integer(1))
  if (length(unique(Ms)) != 1L) {
    stop(
      "merge_filter_traj_by_unit: members disagree on filter_traj month count: ",
      paste(Ms, collapse = ", ")
    )
  }
  M <- Ms[1L]
  if (!all(vapply(ft_list, function(a) dim(a)[2], integer(1)) == 1L)) {
    stop(
      "merge_filter_traj_by_unit: every member filter_traj must have a ",
      "singleton 2nd axis (1 trajectory per array)."
    )
  }

  gids_region <- unlist(gids_list, use.names = FALSE)
  U_region <- length(gids_region)
  # Statename-major block layout over the contiguous region positions; dim[1]
  # == length(statenames)*U_region matches dim(ft)[1]/length(statenames) in
  # every downstream consumer (e.g. filter_traj_to_warmup()'s n_smooth).
  row_names <- paste0(
    rep(statenames, each = U_region),
    rep(seq_len(U_region), times = length(statenames))
  )
  out <- array(
    NA_real_,
    dim = c(length(statenames) * U_region, 1L, M),
    dimnames = list(row_names, NULL, NULL)
  )

  pos0 <- 0L
  for (k in seq_along(ft_list)) {
    ft <- ft_list[[k]]
    U_k <- length(gids_list[[k]])
    if (dim(ft)[1] != length(statenames) * U_k) {
      stop(sprintf(
        paste0(
          "merge_filter_traj_by_unit: member %d has %d state ",
          "rows but %d units × %d statenames = %d expected."
        ),
        k,
        dim(ft)[1],
        U_k,
        length(statenames),
        length(statenames) * U_k
      ))
    }
    rn <- dimnames(ft)[[1]]
    for (lp in seq_len(U_k)) {
      rp <- pos0 + lp # region position (contiguous re-index)
      src <- paste0(statenames, lp) # member row names for this unit
      dst <- paste0(statenames, rp) # region row names
      if (!all(src %in% rn)) {
        stop(sprintf(
          paste0(
            "merge_filter_traj_by_unit: member %d is missing ",
            "state row(s) for unit position %d."
          ),
          k,
          lp
        ))
      }
      out[dst, 1L, ] <- ft[src, 1L, ]
    }
    pos0 <- pos0 + U_k
  }
  list(filter_traj = out, gids = as.character(gids_region))
}

#' Assemble region-level smoothing records from per-member (per-country) records.
#'
#' The ENSO-region forward/hindcast run (`filter.partition_by: enso_region`)
#' reuses the EXISTING country-level smoothing caches rather than re-running
#' smoothing: for a region TU, each member country's per-replicate smoothing
#' records are merged into one region record shaped exactly like a country
#' record, so the whole downstream handoff path
#' (`filter_traj_to_warmup → resample_by_loglik → forecast_from_filtered`) runs
#' unchanged. The ADM1 neighbour-contagion kernel — built over the forecast
#' `frame_tu`'s centroids — then spans the region's units across borders
#' automatically (there is no bpfilter in the forward phase).
#'
#' Records are paired by replicate across members (member record lists must have
#' equal length; records are ordered by `rep_id` when present so pairing is by
#' replicate identity, not list position — they come from one smoothing run).
#' For each replicate the region record carries:
#'   - `gids`        = concatenation of member gids (region unit order),
#'   - `filter_traj` = member arrays stacked along the unit axis with contiguous
#'                     `1..U_region` position re-indexing (`merge_filter_traj_by_unit`),
#'   - `filt_months` = the common vector (asserted identical across members),
#'   - `logLik`      = SUM of member logLiks (independent country partitions ⇒
#'                     additive joint log-likelihood; used by `resample_by_loglik`).
#'
#' @param member_recs named list (name = member ISO3) of per-replicate record
#'        lists, each record a `list(filter_traj, gids, filt_months, logLik,
#'        rep_id, draw_d, ...)` as loaded from a `hist_filter_tu=<ISO3>_*.qs`
#'        cache. Records with a NULL `filter_traj` must be dropped by the caller.
#' @param region_tu character region label (e.g. "enso_east").
#' @param statenames per-unit state names (default `.SP_UNIT_STATENAMES`).
#' @return list of region records (one per replicate).
assemble_region_smoothing_records <- function(
  member_recs,
  region_tu,
  statenames = NULL
) {
  if (length(member_recs) == 0L) {
    stop(sprintf(
      "assemble_region_smoothing_records: region '%s' has no member records.",
      region_tu
    ))
  }
  member_tus <- names(member_recs)
  if (is.null(member_tus) || any(!nzchar(member_tus))) {
    stop(
      "assemble_region_smoothing_records: `member_recs` must be a named list ",
      "(name = member ISO3)."
    )
  }

  n_reps <- vapply(member_recs, length, integer(1))
  if (length(unique(n_reps)) != 1L) {
    stop(sprintf(
      paste0(
        "assemble_region_smoothing_records: region '%s' members ",
        "disagree on replicate count (%s) — cannot pair by replicate."
      ),
      region_tu,
      paste(sprintf("%s=%d", member_tus, n_reps), collapse = ", ")
    ))
  }
  R <- n_reps[[1L]]
  if (R == 0L) {
    stop(sprintf(
      "assemble_region_smoothing_records: region '%s' has 0 replicates.",
      region_tu
    ))
  }

  # Order each member's records by rep_id (when present) so replicate r refers
  # to the same smoothing replicate across all members.
  member_recs <- lapply(member_recs, function(recs) {
    rid <- vapply(
      recs,
      function(r) as.integer(r$rep_id %||% NA_integer_),
      integer(1)
    )
    if (!any(is.na(rid))) recs[order(rid)] else recs
  })
  # Pairing by position is only valid when every member holds the SAME
  # replicate set: failed replicates are dropped per member by the caller, so
  # two members that each lost a different replicate have equal counts but
  # different rep_ids and would be silently mis-paired.
  rid_sets <- lapply(member_recs, function(recs) {
    vapply(recs, function(r) as.integer(r$rep_id %||% NA_integer_), integer(1))
  })
  if (all(vapply(rid_sets, function(x) !anyNA(x), logical(1)))) {
    same <- vapply(
      rid_sets,
      function(x) identical(x, rid_sets[[1L]]),
      logical(1)
    )
    if (!all(same)) {
      stop(
        sprintf(
          paste0(
            "assemble_region_smoothing_records: partition '%s' members ",
            "hold different replicate sets (member '%s' differs from '%s'; ",
            "rep_ids %s vs %s) -- a failed replicate in one member cannot ",
            "be paired. Re-run smoothing natively for this partition."
          ),
          region_tu,
          member_tus[which(!same)[1L]],
          member_tus[[1L]],
          paste(rid_sets[[which(!same)[1L]]], collapse = ","),
          paste(rid_sets[[1L]], collapse = ",")
        ),
        call. = FALSE
      )
    }
  }

  out <- vector("list", R)
  for (r in seq_len(R)) {
    recs_r <- lapply(member_recs, function(recs) recs[[r]])
    ft_list <- lapply(recs_r, function(x) x$filter_traj)
    gids_list <- lapply(recs_r, function(x) as.character(x$gids))

    # filt_months must be identical across members for this replicate (shared
    # smoothing coverage) — the interior-month handoff slice depends on it.
    fm0 <- as.Date(recs_r[[1L]]$filt_months)
    for (j in seq_along(recs_r)) {
      fmj <- as.Date(recs_r[[j]]$filt_months)
      if (!identical(as.character(fmj), as.character(fm0))) {
        stop(sprintf(
          paste0(
            "assemble_region_smoothing_records: region '%s' ",
            "replicate %d: member '%s' filt_months differ from '%s'."
          ),
          region_tu,
          r,
          member_tus[[j]],
          member_tus[[1L]]
        ))
      }
    }

    merged <- merge_filter_traj_by_unit(ft_list, gids_list, statenames)
    ll <- sum(vapply(
      recs_r,
      function(x) as.numeric(x$logLik %||% NA_real_),
      numeric(1)
    ))
    # The members' common smoothing rep_id (asserted identical above), so an
    # assembled record names the same smoothing replicate as a native one --
    # 07b's handoff_map.csv pairs forecast replicates with smoothing events by
    # it. Positional only when the members' records carry no rep_id.
    rid_r <- as.integer(recs_r[[1L]]$rep_id %||% NA_integer_)
    out[[r]] <- list(
      tu = region_tu,
      rep_id = if (is.na(rid_r)) r else rid_r,
      draw_d = as.integer(recs_r[[1L]]$draw_d %||% NA_integer_),
      filter_traj = merged$filter_traj,
      gids = merged$gids,
      filt_months = fm0,
      logLik = ll
    )
  }
  out
}

# Partition-neutral name for the same stitcher (used by partition_by: registry
# as well as enso_region): a TU of any registry is assembled from its member
# countries' caches exactly like an ENSO region.
assemble_partition_smoothing_records <- assemble_region_smoothing_records

#' Resample a list of warmup mini-panels by log-likelihood weights (SIR step).
#'
#' Weights are w_i = exp(loglik_i - max(loglik)) (log-sum-exp shift, avoids
#' overflow that would otherwise zero out and drop the best particle).
#' Resampling is with
#' replacement so the output list has exactly n elements.  This step bridges
#' the N-independent filter runs to the N forecast replicates.
#'
#' @param warmup_list  list of n mini-panel tibbles (from filter_traj_to_warmup).
#' @param logliks      numeric vector of length n; one logLik per filter run.
#' @param n            integer; number of replicates to return (default = length).
#' @return list of n tibbles sampled from warmup_list with probability ∝ exp(loglik).
resample_by_loglik <- function(warmup_list, logliks, n = length(warmup_list)) {
  stopifnot(length(warmup_list) == length(logliks))
  warmup_list[.loglik_resample_index(logliks, n)]
}

#' The SIR draw behind resample_by_loglik(): n indices into the m filter runs,
#' with probability proportional to exp(logLik). One sample() call, so the RNG
#' stream is the same whether the caller wants the panels or the indices.
.loglik_resample_index <- function(logliks, n) {
  wts <- exp(logliks - max(logliks, na.rm = TRUE))
  wts[!is.finite(wts)] <- 0
  if (sum(wts) == 0) {
    warning(
      "resample_by_loglik: all weights zero or non-finite — resampling uniformly."
    )
    wts <- rep(1, length(wts))
  }
  sample(length(logliks), n, replace = TRUE, prob = wts)
}

#' Pick the forecast replicates' handoff states from the smoothing replicates.
#'
#' `weighting = "uniform"` (default): every smoothing
#' replicate counts equally. Forecast replicate r (global id = offset + r) takes
#' smoothing replicate ((offset + r - 1) mod m) + 1, a balanced, deterministic
#' assignment: chunk-safe (the global id, not the chunk-local one, picks the
#' state), and the same replicate id starts from the same smoothing replicate
#' in every arm (factual / counterfactual / hindcast). Consumes no RNG.
#'
#' `weighting = "shuffled"`: balanced like "uniform", but the global ids fall
#' into blocks of m consecutive replicates and block j = floor((id - 1) / m)
#' takes a permutation of 1..m drawn under set.seed(seed + .HANDOFF_SEED_OFFSET
#' + j): each smoothing replicate is used exactly once per full block (counts
#' differ by at most one over a partial block). Under "uniform" the state is
#' (id - 1) mod m, so a cycle posterior draw (id - 1) mod S always meets the
#' same states when gcd(S, m) > 1 (S = 200, m = 100: one state per draw);
#' shuffling breaks that tie for any S. Still chunk-safe and shared across arms
#' (it depends only on the global id and `seed`); the global RNG stream is
#' preserved.
#'
#' `weighting = "loglik"`: SIR resampling by exp(logLik) (resample_by_loglik()),
#' The smoothing replicates are independent filter runs
#' whose whole-period log-likelihoods differ by thousands of log units, so the
#' effective sample size is ~1 and every forecast replicate starts from one
#' smoothing replicate.
#'
#' `weighting = "joint"` (param_draw_mode "joint"): forecast
#' replicate r (global id) takes the smoothing replicate whose rep_id is
#' j(r) = ((r - 1) mod n_smoothing) + 1, matched by rep_id rather than by
#' position, so the pairing is the same in every unit and arm even when a unit
#' lost a replicate. A missing rep_id is an error: its draw has no state in
#' that unit. Consumes no RNG.
#'
#' @param warmup_list list of m warmup mini-panels (one per smoothing replicate).
#' @param logliks     numeric(m); used by "loglik" and for the ESS message.
#' @param n           number of forecast replicates in this invocation.
#' @param offset      global replicate-id offset of this invocation (chunking).
#' @param weighting   "uniform" | "shuffled" | "loglik" | "joint".
#' @param seed        simulator.seed; read only by "shuffled" (required there).
#' @param rep_ids     integer(m): the smoothing rep_id of each element of
#'   `warmup_list`; read only by "joint" (required there).
#' @param n_smoothing the smoothing cache's replicate count (meta.yml nrep);
#'   read only by "joint" (required there).
#' @return list of n warmup mini-panels, with attribute `idx`: the position in
#'   `warmup_list` each forecast replicate took (07b writes it to the
#'   partition's handoff_map.csv, which the results readers use to pair a
#'   forecast replicate with its smoothing replicate).
select_handoff_warmup <- function(
  warmup_list,
  logliks,
  n,
  offset = 0L,
  weighting = c("uniform", "shuffled", "loglik", "joint"),
  seed = NULL,
  rep_ids = NULL,
  n_smoothing = NULL
) {
  weighting <- match.arg(weighting)
  m <- length(warmup_list)
  stopifnot(m >= 1L, length(logliks) == m)
  ids <- as.integer(offset) + seq_len(n)
  idx <- switch(
    weighting,
    loglik = .loglik_resample_index(logliks, n),
    uniform = ((ids - 1L) %% m) + 1L,
    shuffled = .shuffled_handoff_index(ids, m, seed),
    joint = .joint_handoff_index(ids, rep_ids, n_smoothing, m)
  )
  out <- warmup_list[idx]
  attr(out, "idx") <- idx
  out
}

#' The "shuffled" handoff index: global replicate id r -> perm_j[((r - 1) mod m) + 1]
#' with j = floor((r - 1) / m) and perm_j = a permutation of 1..m seeded by
#' (seed, j) alone, so any slice of ids gets the same pairing as the full run.
#' Consumes no global RNG (.with_preserved_rng()).
#' @param ids  integer global replicate ids (>= 1).
#' @param m    number of smoothing replicates.
#' @param seed simulator.seed.
#' @return integer(length(ids)) in 1..m.
.shuffled_handoff_index <- function(ids, m, seed) {
  if (
    is.null(seed) ||
      length(seed) != 1L ||
      is.na(suppressWarnings(as.integer(seed)))
  ) {
    stop(
      "select_handoff_warmup: handoff_weighting 'shuffled' needs the simulator seed.",
      call. = FALSE
    )
  }
  ids <- as.integer(ids)
  if (!length(ids)) {
    return(integer(0))
  }
  blk <- (ids - 1L) %/% m
  pos <- ((ids - 1L) %% m) + 1L
  blocks <- sort(unique(blk))
  perms <- .with_preserved_rng(lapply(blocks, function(j) {
    set.seed(as.integer(seed) + .HANDOFF_SEED_OFFSET + j)
    sample.int(m)
  }))
  vapply(
    seq_along(ids),
    function(i) perms[[match(blk[i], blocks)]][pos[i]],
    integer(1)
  )
}

#' The "joint" handoff index: global replicate id r -> the position in the
#' handoff list of smoothing rep_id ((r - 1) mod n_smoothing) + 1.
#' @param ids         integer global replicate ids (>= 1).
#' @param rep_ids     integer(m) smoothing rep_ids of the handoff list.
#' @param n_smoothing smoothing replicates per unit (meta.yml nrep).
#' @param m           length of the handoff list.
#' @return integer(length(ids)) in 1..m.
.joint_handoff_index <- function(ids, rep_ids, n_smoothing, m) {
  n_sm <- suppressWarnings(as.integer(n_smoothing %||% NA))
  if (length(n_sm) != 1L || is.na(n_sm) || n_sm < 1L) {
    stop(
      "select_handoff_warmup: the joint pairing needs the smoothing replicate count (meta.yml nrep).",
      call. = FALSE
    )
  }
  rep_ids <- as.integer(rep_ids %||% integer(0))
  if (length(rep_ids) != m || anyNA(rep_ids) || anyDuplicated(rep_ids)) {
    stop(
      "select_handoff_warmup: the joint pairing needs one distinct smoothing rep_id per handoff state.",
      call. = FALSE
    )
  }
  want <- ((as.integer(ids) - 1L) %% n_sm) + 1L
  idx <- match(want, rep_ids)
  if (anyNA(idx)) {
    stop(
      sprintf(
        paste0(
          "select_handoff_warmup: joint pairing: smoothing replicate(s) %s have no handoff ",
          "state in this unit (failed smoothing run or extraction), so their posterior ",
          "draw has no reconstructed state. Rerun those smoothing replicates."
        ),
        paste(sort(unique(want[is.na(idx)])), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  idx
}

#' The forecast-replicate -> smoothing-replicate pairing of one TU's handoff,
#' as 07b writes it to the partition's handoff_map.csv.
#'
#' @param tu               TU label.
#' @param ok_records       the smoothing records whose handoff extraction
#'                         succeeded, in the order passed to
#'                         select_handoff_warmup() (each carries `rep_id`, the
#'                         smoothing run's global replicate id).
#' @param warmup_list      select_handoff_warmup()'s result (attribute `idx`).
#' @param replicate_offset global replicate-id offset of this invocation.
#' @param handoff_month    the TU's handoff month.
#' @param weighting        simulator.forecast.handoff_weighting ("joint" under
#'                         param_draw_mode joint).
#' @param draws            NULL (default: the table without a draw column)
#'                         or integer(length(warmup_list)): each forecast
#'                         replicate's posterior draw, written as column
#'                         `smoothing_draw` (param_draw_mode joint: the draw its
#'                         handoff state was reconstructed under).
#' @return data.frame(tu, replicate_id, smoothing_rep_id, handoff_month,
#'   handoff_weighting[, smoothing_draw]), one row per forecast replicate.
handoff_map_table <- function(
  tu,
  ok_records,
  warmup_list,
  replicate_offset,
  handoff_month,
  weighting,
  draws = NULL
) {
  idx <- attr(warmup_list, "idx")
  if (is.null(idx) || length(idx) != length(warmup_list)) {
    stop(
      "handoff_map_table: `warmup_list` carries no `idx` attribute -- pass ",
      "select_handoff_warmup()'s result before anything re-subsets it.",
      call. = FALSE
    )
  }
  rid <- vapply(
    ok_records,
    function(r) as.integer(r$rep_id %||% NA_integer_),
    integer(1)
  )
  if (max(idx) > length(rid)) {
    stop(
      "handoff_map_table: `ok_records` is shorter than the list select_handoff_warmup() drew from.",
      call. = FALSE
    )
  }
  out <- data.frame(
    tu = tu,
    replicate_id = as.integer(replicate_offset %||% 0L) + seq_along(idx),
    smoothing_rep_id = rid[idx],
    handoff_month = format(as.Date(handoff_month)),
    handoff_weighting = weighting,
    stringsAsFactors = FALSE
  )
  if (!is.null(draws)) {
    if (length(draws) != length(idx)) {
      stop(
        "handoff_map_table: `draws` must hold one draw per forecast replicate.",
        call. = FALSE
      )
    }
    out$smoothing_draw <- as.integer(draws)
  }
  out
}

#' The posterior draw of each handoff state: the cache record's global
#' `draw_d` (the draw its smoothing replicate ran), not the particle state's
#' `posterior_draw_idx` (the local index into the one-draw subset, always 1).
#' @param ok_records  the records select_handoff_warmup() drew from, in order.
#' @param warmup_list select_handoff_warmup()'s result (attribute `idx`).
#' @return integer(length(warmup_list)).
handoff_draw_d <- function(ok_records, warmup_list) {
  idx <- attr(warmup_list, "idx")
  if (is.null(idx) || length(idx) != length(warmup_list)) {
    stop(
      "handoff_draw_d: `warmup_list` carries no `idx` attribute.",
      call. = FALSE
    )
  }
  d <- vapply(
    ok_records,
    function(r) as.integer(r$draw_d %||% NA_integer_),
    integer(1)
  )
  out <- d[idx]
  if (anyNA(out)) {
    stop(
      "handoff_draw_d: a handoff state's smoothing record carries no draw_d.",
      call. = FALSE
    )
  }
  out
}

#' Effective sample size of exp(logLik) weights (diagnostic for the handoff).
handoff_loglik_ess <- function(logliks) {
  w <- exp(logliks - max(logliks, na.rm = TRUE))
  w[!is.finite(w)] <- 0
  if (sum(w) == 0) {
    return(NA_real_)
  }
  w <- w / sum(w)
  1 / sum(w^2)
}

#' Anchor the forecast handoff active/mso to the OBSERVED active footprint (N4).
#'
#' Diagnostic-only (hindcast) IC anchor. Replaces every warmup mini-panel's
#' `active`/`mso` columns with the observed active set at the launch month, so
#' the forecast starts from the KNOWN observed active configuration instead of
#' the (cold-start under-recruited) smoothed one. This bypasses the PF joint-
#' coverage saturation at t0 (Deficit 1: ~4.6 smoothed vs ~16 observed active for
#' NGA) and, run OFF vs ON, partitions the "handoff recruitment deficit" from the
#' "forecast-ramp deficit". NOT valid for a true out-of-sample forecast — the
#' observed footprint is only knowable for a hindcast.
#'
#' The override is identical across replicates (it does not depend on the
#' particle), so applying it to the resampled `warmup_list` is equivalent to
#' applying it to the raw per-trajectory panels, and it composes cleanly with the
#' `spread_param_draws` posterior-draw spreading (which only touches draw_idx).
#'
#' @param warmup_list list of mini-panels (each carries `gid`, `active`, `mso`,
#'   ... as built by filter_traj_to_warmup()/filtered_states_to_warmup()).
#' @param obs_dur_panel observed duration panel; one row per observed
#'   active (gid, month), carrying columns `gid`, `month`, `mso`.
#' @param launch_month Date; the last observed month before the forecast
#'   (= handoff_date - 1 month, i.e. each mini-panel's t_start). A gid with an
#'   observed episode covering this month is active=1 at IC with its observed
#'   mso; every other gid is active=0/mso=0.
#' @return `warmup_list` with `active`/`mso` overwritten from the observed
#'   footprint. Other columns (lags, TTLO, posterior_draw_idx) are untouched.
anchor_warmup_active_to_observed <- function(
  warmup_list,
  obs_dur_panel,
  launch_month
) {
  if (length(warmup_list) == 0L) {
    return(warmup_list)
  }
  need <- c("gid", "month", "mso")
  miss <- setdiff(need, names(obs_dur_panel))
  if (length(miss)) {
    stop(
      "anchor_warmup_active_to_observed: obs_dur_panel is missing column(s): ",
      paste(miss, collapse = ", ")
    )
  }

  lm <- as.Date(launch_month)
  obs <- obs_dur_panel[
    as.Date(obs_dur_panel$month) == lm,
    c("gid", "mso"),
    drop = FALSE
  ]
  # One mso per gid: an episode covers a given month at most once, but guard a
  # duplicated (gid, month) row by keeping the largest mso (the longest-running
  # episode's month index).
  mso_by_gid <- if (nrow(obs) == 0L) {
    integer(0)
  } else {
    tapply(as.integer(obs$mso), as.character(obs$gid), max)
  }

  lapply(warmup_list, function(mp) {
    g <- as.character(mp$gid)
    hit <- g %in% names(mso_by_gid)
    mp$active <- as.integer(hit)
    mp$mso <- ifelse(hit, as.integer(mso_by_gid[g]), 0L)
    mp$mso <- as.integer(mp$mso)
    mp
  })
}

# ---------------------------------------------------------------------------
# Forecast from filtered particles (calls simulate_one_tu from 09a)
# ---------------------------------------------------------------------------

#' Derive the per-replicate posterior draw index for each track at the forecast
#' handoff. Pulled out of forecast_from_filtered() so the branch logic can be
#' unit-tested without running a forecast (mirrors .occ_bundle_from_draws_df()).
#'
#' Three modes select the per-replicate occurrence draw (`param_draw_mode`, with
#' the legacy `spread_param_draws` boolean kept as a fallback):
#'
#'   "inherited" (spread_param_draws = FALSE): each replicate inherits its
#'     handoff particle's posterior_draw_idx — the oldest behaviour. When the
#'     smoothing cache was built with n_draws_in_filter=1 every particle carries
#'     draw_idx=1, so the whole ensemble collapses to one posterior draw and the
#'     mon predictive tail flattens.
#'   "cycle" (spread_param_draws = TRUE): each replicate gets a distinct draw
#'     cycling 1..S over the full bundle (via replicate_ids), restoring
#'     posterior-parameter spread deterministically.
#'   "joint": each replicate runs the global posterior draw
#'     its handoff state was reconstructed under (`handoff_draw_d`, the
#'     smoothing record's draw_d), so state and parameters are one joint
#'     posterior draw. Needs a shared_cycle smoothing cache (07b checks
#'     meta.yml) and draws within 1..S_occ (no wrap: a wrapped index would be
#'     another draw).
#'   "random" (new global default): each replicate draws a FRESH random index in
#'     1..S_occ, seeded per replicate by seed + .DRAW_SEED_OFFSET + replicate_id
#'     so the mapping is chunk-safe (invariant to replicate_offset) and, in an
#'     enso_region partition, shared across all member countries' units (one
#'     per-replicate index broadcast by the kernel). The handoff draw_idx is
#'     ignored.
#'
#' In every mode dur/mon wrap off the clamped occurrence index, so the three
#' tracks stay paired (one posterior sample per replicate across
#' occurrence/duration/monthly).
#'
#' @param handoff_list list of Np warmup mini-panels (each carries
#'   posterior_draw_idx); read only in the "inherited" mode.
#' @param replicate_ids integer[Np] = seq_len(Np) + replicate_offset.
#' @param S_occ number of draws in the occurrence bundle.
#' @param S_dur,S_mon number of draws in the dur/mon bundles, or NULL when that
#'   track is inactive (then its index is NULL, matching the legacy path).
#' @param spread_param_draws logical; legacy alias honoured only when
#'   param_draw_mode is NULL/blank (TRUE -> "cycle", FALSE -> "inherited").
#' @param param_draw_mode character one of "random"/"cycle"/"inherited"/"joint",
#'   or NULL to fall back to spread_param_draws.
#' @param seed integer per-TU seed reused (with .DRAW_SEED_OFFSET) by "random".
#' @param handoff_draw_d integer[Np]: each handoff state's global draw
#'   (handoff_draw_d()); read only (and required) in the "joint" mode.
#' @return list(occ, dur, mon) of integer[Np] index vectors; dur/mon are
#'   NULL when their S argument is NULL.
.forecast_param_indices <- function(
  handoff_list,
  replicate_ids,
  S_occ,
  S_dur = NULL,
  S_mon = NULL,
  spread_param_draws = FALSE,
  param_draw_mode = NULL,
  seed = 20260101L,
  handoff_draw_d = NULL
) {
  mode <- param_draw_mode
  if (is.null(mode) || !nzchar(mode)) {
    mode <- if (isTRUE(spread_param_draws)) "cycle" else "inherited"
  }
  mode <- match.arg(mode, c("random", "cycle", "inherited", "joint"))
  if (identical(mode, "joint")) {
    hd <- as.integer(handoff_draw_d %||% integer(0))
    if (length(hd) != length(handoff_list) || anyNA(hd)) {
      stop(
        "param_draw_mode joint: one global smoothing draw per handoff state is required ",
        "(handoff_draw_d).",
        call. = FALSE
      )
    }
    if (any(hd < 1L | hd > S_occ)) {
      stop(
        sprintf(
          paste0(
            "param_draw_mode joint: handoff draws span %d..%d but the forecast bundle ",
            "holds %d draws; the forecast simulator.max_draws must equal the smoothing's."
          ),
          min(hd),
          max(hd),
          as.integer(S_occ)
        ),
        call. = FALSE
      )
    }
  }

  base_idx <- switch(
    mode,
    joint = as.integer(handoff_draw_d),
    cycle = as.integer(replicate_ids),
    inherited = vapply(
      handoff_list,
      function(mp) as.integer(mp$posterior_draw_idx[1L]),
      integer(1)
    ),
    random = .with_preserved_rng(vapply(
      replicate_ids,
      function(r) {
        set.seed(as.integer(seed) + .DRAW_SEED_OFFSET + as.integer(r))
        sample.int(S_occ, 1L)
      },
      integer(1)
    ))
  )
  idx_occ <- ((base_idx - 1L) %% S_occ) + 1L
  list(
    occ = idx_occ,
    dur = if (!is.null(S_dur)) ((idx_occ - 1L) %% S_dur) + 1L else NULL,
    mon = if (!is.null(S_mon)) ((idx_occ - 1L) %% S_mon) + 1L else NULL
  )
}

#' Run the forecast phase after bpfilter, using each particle as a replicate.
#'
#' Each of the Np particles becomes one forecast replicate. The posterior draw
#' index for each track is derived by .forecast_param_indices() via param_draw_mode:
#' "random" (default) = a fresh chunk-safe random draw per replicate (shared
#' across a region's member units); "cycle" = deterministic 1..S spread;
#' "inherited" = each particle's filtered draw_idx. The legacy spread_param_draws
#' boolean is honoured only when param_draw_mode is NULL/blank. ENSO scenarios are
#' passed through unchanged.
#'
#' This function is the M3 handoff: it bridges bpfilter → simulate_one_tu,
#' preserving the existing 09b/09c output schema.
#'
#' @param handoff_list     list of Np mini-panels from filtered_states_to_warmup().
#' @param tu               TU label (character).
#' @param frame_tu         location frame for this TU.
#' @param panel_tu_full    climate panel for this TU covering the forecast horizon.
#' @param draws_occ        full occurrence draws bundle (not subsetted).
#' @param draws_dur        Optional duration (cessation-hazard) draws bundle
#'                          (`extract_duration_draws()`). NULL
#'                          (default) leaves the forecast on
#'                          `simulate_one_tu()`'s legacy (non-duration-aware)
#'                          path. When supplied, `posterior_idx_dur` is
#'                          derived by wrapping `posterior_idx_occ`, so the
#'                          draws stay paired to the filtered particle's
#'                          occurrence draw.
#' @param draws_mon        Optional monthly-case (`mon`) draws bundle
#'                          (`extract_monthly_draws()`). NULL (default) leaves
#'                          the forecast without case counts. Supplied,
#'                          `posterior_idx_mon` is wrapped off
#'                          `posterior_idx_occ` exactly like `_dur`, keeping
#'                          all parameter sets on the filtered particle's draw.
#' @param cb_meta_mon      Optional cb_meta sidecar for the mon fit; forwarded to
#'                          `simulate_one_tu()`, which collapses its `mso` entry
#'                          into the months-since-onset basis table.
#' @param cb_meta_dur      Optional cb_meta sidecar for the `dur` fit; forwarded
#'                          to `simulate_one_tu()`. The deployed best-lag
#'                          duration fits are `_cb` fits with `enso_cb<k>` fixed
#'                          effects, and the duration best lag is selected
#'                          independently of occ/mon, so the duration basis
#'                          must come from its own sidecar. Omitting it leaves
#'                          the cessation hazard ENSO-blind (`evaluate_lp()`
#'                          silently skips a FE whose covariate is absent).
#' @param enso_scenarios   named list of ENSO scenario matrices (or NULL).
#' @param horizon          list(start=, end=) — the FORECAST horizon.
#' @param output_part      output directory for this TU.
#' @param out_fmt          output format ("qs" or "rds").
#' @param dry_run          logical; if TRUE, no files written.
#' @param seed             integer root seed.
#' @param force_o_tu       logical; forwarded to simulate_one_tu(). See that
#'                          function's docstring for the partition_by=country
#'                          rationale. Default FALSE.
#' @param incidence_cap    NULL or a number > 0; forwarded to simulate_one_tu()
#'                          (monthly cases <= incidence_cap x pop).
#' @param onset_eligibility "inactive" / "season"; forwarded to simulate_one_tu()
#'                          (simulator.onset_eligibility). Must match the
#'                          smoothing cache's (meta.yml identity field).
#' @param handoff_draw_d   NULL, or integer[Np]: each handoff state's global
#'                          smoothing draw (handoff_draw_d()); required by
#'                          param_draw_mode "joint", ignored otherwise.
#' @param cessation_exposure "full" / "half"; forwarded to simulate_one_tu()
#'                          (simulator.cessation_exposure). Must match the
#'                          smoothing cache's (meta.yml identity field).
#' @param replicate_offset Global replicate-ID offset for this call (default 0).
#'   Supports parallel replicate-chunking: when multiple concurrent processes
#'   each forecast a disjoint slice of one TU's replicates (e.g. one SLURM
#'   array task per country, fanned out further over replicate chunks), each
#'   chunk's `handoff_list` is only Np-long (its own chunk size) but must not
#'   reuse replicate_id 1..Np — that would collide with every other chunk of
#'   the same TU downstream (07b/07c key on replicate_id as a per-TU-unique
#'   ensemble-member label). `replicate_ids` becomes seq_len(Np) + offset.
#'   `warmup_panels_list`/`enso_scenarios` are looked up by simulate_one_tu()
#'   using the *value* of replicate_id (simulator_core.R:
#'   `warmup_panels_list[[rep_id]]`, `enso_scenarios[[ec]][, rep_id]`), not by
#'   loop position, so they must be left-padded by `offset` placeholder
#'   entries/columns (never read, since this call's replicate_ids never index
#'   into the padding) to keep the value-indexed lookups landing correctly.
forecast_from_filtered <- function(
  handoff_list,
  tu,
  frame_tu,
  panel_tu_full,
  draws_occ,
  draws_dur = NULL,
  draws_mon = NULL,
  enso_scenarios = NULL,
  horizon,
  output_part,
  out_fmt = "qs",
  dry_run = FALSE,
  seed = 20260101L,
  cb_meta_occ = NULL,
  cb_meta_mon = NULL,
  cb_meta_dur = NULL,
  grouping = NULL,
  force_o_tu = FALSE,
  replicate_offset = 0L,
  drop_country_tu_spillover = FALSE,
  spread_param_draws = FALSE,
  param_draw_mode = NULL,
  handoff_draw_d = NULL,
  use_nbr = FALSE,
  nbr_rho_km = 500L,
  nbr_gadm_path = "data/geodata/gadm_valid.gpkg",
  nbr_cache_dir = "generated_data/",
  nbr_spec = NULL,
  covariate_feed = "latent",
  delta_fn = .SP_STATE_DELTA_FN,
  delta_fp = .SP_STATE_DELTA_FP,
  clamp_enso_to_knots = FALSE,
  clamp_enso_margin = 0,
  incidence_cap = NULL,
  onset_eligibility = "inactive",
  cessation_exposure = "full"
) {
  Np <- length(handoff_list)
  if (Np == 0L) {
    stop("forecast_from_filtered: empty handoff_list.")
  }

  replicate_offset <- as.integer(replicate_offset %||% 0L)
  replicate_ids <- seq_len(Np) + replicate_offset
  # Per-replicate posterior draw index for each track, via three modes
  # (param_draw_mode, falling back to the legacy spread_param_draws boolean):
  # "random" (new default) draws a fresh chunk-safe random index per replicate;
  # "cycle" spreads draws deterministically 1..S; "inherited" reuses each handoff
  # particle's posterior_draw_idx (collapses to one draw when the smoothing cache
  # was built with n_draws_in_filter=1). dur/mon always wrap off the clamped
  # occurrence index, so cross-track pairing is preserved in every mode. See
  # .forecast_param_indices().
  .pidx <- .forecast_param_indices(
    handoff_list,
    replicate_ids,
    S_occ = draws_occ$S,
    S_dur = if (!is.null(draws_dur)) draws_dur$S else NULL,
    S_mon = if (!is.null(draws_mon)) draws_mon$S else NULL,
    spread_param_draws = spread_param_draws,
    param_draw_mode = param_draw_mode,
    seed = seed,
    handoff_draw_d = handoff_draw_d
  )
  posterior_idx_occ <- .pidx$occ
  posterior_idx_dur <- .pidx$dur
  posterior_idx_mon <- .pidx$mon

  warmup_panels_list <- handoff_list
  if (replicate_offset > 0L) {
    warmup_panels_list <- c(vector("list", replicate_offset), handoff_list)
    if (!is.null(enso_scenarios)) {
      enso_scenarios <- lapply(enso_scenarios, function(mat) {
        pad <- matrix(NA_real_, nrow = nrow(mat), ncol = replicate_offset)
        cbind(pad, mat)
      })
    }
  }

  simulate_one_tu(
    tu = tu,
    frame_tu = frame_tu,
    panel_tu = panel_tu_full,
    draws_occ = draws_occ,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    replicate_ids = replicate_ids,
    posterior_idx_occ = posterior_idx_occ,
    posterior_idx_dur = posterior_idx_dur,
    posterior_idx_mon = posterior_idx_mon,
    seed = seed,
    output_part = output_part,
    horizon = horizon,
    out_fmt = out_fmt,
    dry_run = dry_run,
    enso_scenarios = enso_scenarios,
    warmup_panels_list = warmup_panels_list,
    cb_meta_occ = cb_meta_occ,
    cb_meta_mon = cb_meta_mon,
    cb_meta_dur = cb_meta_dur,
    grouping = grouping,
    force_o_tu = force_o_tu,
    drop_country_tu_spillover = drop_country_tu_spillover,
    use_nbr = use_nbr,
    nbr_rho_km = nbr_rho_km,
    nbr_gadm_path = nbr_gadm_path,
    nbr_cache_dir = nbr_cache_dir,
    nbr_spec = nbr_spec,
    covariate_feed = covariate_feed,
    delta_fn = delta_fn,
    delta_fp = delta_fp,
    clamp_enso_to_knots = clamp_enso_to_knots,
    clamp_enso_margin = clamp_enso_margin,
    incidence_cap = incidence_cap,
    onset_eligibility = onset_eligibility,
    cessation_exposure = cessation_exposure
  )
}

# ---------------------------------------------------------------------------
# Smoothing cache-dir identity (meta.yml)
# ---------------------------------------------------------------------------
# 08c writes <cache_dir>/meta.yml next to the hist_filter_tu=*.qs files: the
# partition and fit-grouping identity of the smoothing run (registry names +
# content digests) plus Np / nrep / schema. The record format and
# .SM_CACHE_SCHEMA_VERSION are unchanged, so every existing cache stays
# readable; a dir without meta.yml loads with one warning.

.SM_META_IDENTITY_KEYS <- c(
  "partition_by",
  "partition_source",
  "partition_digest"
)

smoothing_cache_meta_path <- function(cache_dir) {
  file.path(cache_dir, "meta.yml")
}

#' smoothing_cache_meta
#'
#' The identity block for a smoothing run (written by 08c, expected by 07b).
#' `partition_source` / `partition_digest` are "" unless partition_by is
#' "registry"; `grouping_*` are "" when grouping is off.
smoothing_cache_meta <- function(
  partition_by,
  partition_source = "",
  grouping_source = "",
  fits_dir = "",
  Np = NA_integer_,
  nrep = NA_integer_,
  schema = sprintf("sv%02d", .SM_CACHE_SCHEMA_VERSION)
) {
  reg <- identical(partition_by, "registry")
  ps <- if (reg) as.character(partition_source %||% "") else ""
  gs <- as.character(grouping_source %||% "")
  list(
    partition_by = as.character(partition_by),
    partition_source = ps,
    partition_digest = if (nzchar(ps)) registry_digest(ps) else "",
    grouping_source = gs,
    grouping_digest = if (nzchar(gs)) registry_digest(gs) else "",
    fits_dir = as.character(fits_dir %||% ""),
    Np = as.integer(Np),
    nrep = as.integer(nrep),
    schema = schema
  )
}

# Optional identity fields of a smoothing cache dir:
# field name -> legacy value, the value every cache built before the field
# existed was built under. Each field is registered here once, and both
# sides follow one rule:
#   - written (smoothing_cache_identity()) only when the value is not the
#     legacy one, so a default run's meta.yml is unchanged;
#   - an absent field reads as the legacy value;
#   - a mismatch is an error naming the field and both values
#     (.cache_meta_optional_check(), called by check_smoothing_cache_meta()).
# Registered:
#   pop_year: simulator.location_frame.pop_year of the
#     smoothing; legacy 2020 (simulator_core.R::.SIM_POP_YEAR_DEFAULT). A
#     forecast reading the cache must run under the same population year.
#   who_scale_digest / db_overlap: the estimated WHO-scale
#     table's content digest (legacy "": fixed scale) and the DB-overlap rule
#     (legacy "score_both") of the smoothing. RECORD-ONLY (.SM_META_RECORD_ONLY):
#     a forecast does not know them, so the read side checks them only when
#     the expectation names them.
#   onset_eligibility: simulator.onset_eligibility of the
#     smoothing; legacy "inactive". A forecast handing off from the cache must
#     run under the same gate (checked always, like pop_year).
#   cessation_exposure: simulator.cessation_exposure of the
#     smoothing; legacy "full". Checked always, like onset_eligibility.
#   draw_design / max_draws: simulator.filter.draw_design
#     (legacy "per_tu_random") and simulator.max_draws (legacy 200, the
#     simulator default) of the smoothing. RECORD-ONLY: checked only when the
#     expectation names them, i.e. by a param_draw_mode joint forecast, which
#     needs a shared_cycle cache thinned to its own max_draws.
.SM_META_OPTIONAL_FIELDS <- list(
  pop_year = 2020L,
  who_scale_digest = "",
  db_overlap = "score_both",
  onset_eligibility = "inactive",
  cessation_exposure = "full"
)
.SM_META_OPTIONAL_FIELDS$draw_design <- "per_tu_random"
.SM_META_OPTIONAL_FIELDS$max_draws <- 200L
.SM_META_RECORD_ONLY <- c("who_scale_digest", "db_overlap")
.SM_META_RECORD_ONLY <- c(.SM_META_RECORD_ONLY, "draw_design", "max_draws")

#' smoothing_cache_optional_fields
#'
#' The registered optional identity fields' values for a smoothing config `S`
#' (smoothing_core.R::smoothing_make_config()), passed by 08c as
#' smoothing_cache_identity(optional =): each is written to meta.yml only when
#' not its legacy value.
smoothing_cache_optional_fields <- function(S) {
  who <- S[["WHO_OBS"]] %||% list()
  list(
    pop_year = S[["POP_YEAR"]] %||% .SM_META_OPTIONAL_FIELDS$pop_year,
    who_scale_digest = S[["WHO_SCALE_DIGEST"]] %||%
      .SM_META_OPTIONAL_FIELDS$who_scale_digest,
    db_overlap = if (isTRUE(who[["enable"]])) {
      who[["db_overlap"]] %||% .SM_META_OPTIONAL_FIELDS$db_overlap
    } else {
      .SM_META_OPTIONAL_FIELDS$db_overlap
    },
    onset_eligibility = S[["ONSET_ELIGIBILITY"]] %||%
      .SM_META_OPTIONAL_FIELDS$onset_eligibility,
    cessation_exposure = S[["CESSATION_EXPOSURE"]] %||%
      .SM_META_OPTIONAL_FIELDS$cessation_exposure,
    draw_design = S[["DRAW_DESIGN"]] %||% .SM_META_OPTIONAL_FIELDS$draw_design,
    max_draws = S[["MAX_DRAWS_OCC"]] %||% .SM_META_OPTIONAL_FIELDS$max_draws
  )
}

# Same value for the identity rule: compared as text, so an integer written to
# YAML and read back as a double still matches.
.cache_meta_same <- function(a, b) {
  identical(
    as.character(unlist(a, use.names = FALSE)),
    as.character(unlist(b, use.names = FALSE))
  )
}

#' .cache_meta_add_optional
#'
#' Add the optional identity fields in `optional` (named list, each name
#' registered in `registry`) to `meta`, each only when its value differs from
#' the field's legacy value. An empty `optional` returns `meta` unchanged.
.cache_meta_add_optional <- function(
  meta,
  optional = list(),
  registry = .SM_META_OPTIONAL_FIELDS
) {
  for (f in names(optional)) {
    if (!f %in% names(registry)) {
      stop(
        sprintf(
          "[smoothing cache] optional identity field '%s' is not registered in .SM_META_OPTIONAL_FIELDS",
          f
        ),
        call. = FALSE
      )
    }
    if (!.cache_meta_same(optional[[f]], registry[[f]])) {
      meta[[f]] <- optional[[f]]
    }
  }
  meta
}

#' .cache_meta_optional_check
#'
#' The read side of the optional-field rule: the cache's value of `field` in
#' `meta` (a read meta.yml list; an absent field reads as `legacy`) must equal
#' `expected` (the run's value, `legacy` when the run does not set it). A
#' mismatch is an error naming the field and both values.
#' @return invisible(TRUE)
.cache_meta_optional_check <- function(
  meta,
  field,
  expected,
  legacy,
  where = "smoothing cache"
) {
  got <- meta[[field]] %||% legacy
  if (!.cache_meta_same(got, expected)) {
    stop(
      sprintf(
        paste0(
          "[%s] identity field '%s' differs: cache '%s' vs run '%s'%s. Rebuild the ",
          "smoothing cache for this setting, or point filter.smoothing_input at ",
          "the matching one."
        ),
        where,
        field,
        paste(as.character(unlist(got)), collapse = ","),
        paste(as.character(unlist(expected)), collapse = ","),
        if (is.null(meta[[field]])) {
          " (absent in the cache's meta.yml: read as the legacy value)"
        } else {
          ""
        }
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' smoothing_cache_identity
#'
#' The identity list 08c writes to <cache_dir>/meta.yml, built in one place
#' from the smoothing config `S` (smoothing_core.R::smoothing_make_config()):
#' smoothing_cache_meta() (partition + fit-grouping registries and digests,
#' Np, nrep, schema token) + `run_tag` + `tus` (the TUs whose handoff caches
#' were written), then the optional fields in `optional` (named list; see
#' .SM_META_OPTIONAL_FIELDS: each written only when not its legacy value).
#' With `optional = list()` the list holds the identity keys only
#' (pinned in tests/test_t2_shared_infra.R).
#' `software`: the smoothing environment
#' (provenance.R::sim_software_versions(), passed by 08c), appended last as
#' `software`. Informational only: registered in neither
#' .SM_META_OPTIONAL_FIELDS nor .SM_META_RECORD_ONLY, so
#' check_smoothing_cache_meta() never reads it and a cache built under other
#' versions still loads. NULL (default) = not written.
smoothing_cache_identity <- function(
  S,
  tus = NULL,
  optional = list(),
  software = NULL
) {
  grouped <- isTRUE(S[["GROUPED"]])
  meta <- smoothing_cache_meta(
    S[["PARTITION_BY"]] %||% "transmission_unit",
    S[["PARTITION_SOURCE"]] %||% "",
    grouping_source = if (grouped) S[["GRP_SOURCE"]] else "",
    fits_dir = if (grouped) S[["GRP_FITS_DIR"]] else "",
    Np = S[["NP"]],
    nrep = S[["N_REPS"]],
    schema = sm_cache_schema_token(S[["COVARIATE_FEED"]] %||% "latent")
  )
  meta$run_tag <- S[["run_tag"]]
  meta$tus <- tus
  meta <- .cache_meta_add_optional(meta, optional)
  if (!is.null(software)) {
    meta$software <- software
  }
  meta
}

#' check_smoothing_cache_meta
#'
#' Compare <cache_dir>/meta.yml with `expect` (a smoothing_cache_meta() list)
#' on the partition identity keys, and on every registered optional field
#' (.SM_META_OPTIONAL_FIELDS, via .cache_meta_optional_check()). NULL `expect`
#' = no check. Missing meta.yml = one warning (older caches).
#' Mismatch = error naming the differing keys.
check_smoothing_cache_meta <- function(cache_dir, expect = NULL) {
  if (is.null(expect) || is.null(cache_dir) || !nzchar(cache_dir)) {
    return(invisible(NULL))
  }
  mp <- smoothing_cache_meta_path(cache_dir)
  if (!file.exists(mp)) {
    warning(
      sprintf(
        paste0(
          "[smoothing cache] %s has no meta.yml (cache built before ",
          "2026-10-16): partition identity not checked."
        ),
        cache_dir
      ),
      call. = FALSE
    )
    return(invisible(NULL))
  }
  got <- yaml::read_yaml(mp)
  # Optional identity fields (.SM_META_OPTIONAL_FIELDS; none registered yet):
  # they describe the smoothing content, not the partition, so they are checked
  # for every cache, before the per-country / enso_region stitch exemptions.
  for (f in names(.SM_META_OPTIONAL_FIELDS)) {
    leg <- .SM_META_OPTIONAL_FIELDS[[f]]
    # record-only fields: checked only when expected
    if (f %in% .SM_META_RECORD_ONLY && is.null(expect[[f]])) {
      next
    }
    .cache_meta_optional_check(got, f, expect[[f]] %||% leg, leg)
  }
  g_pb <- as.character(got$partition_by %||% "")
  e_pb <- as.character(expect$partition_by %||% "")
  # Per-country caches carry no registry: they are the legitimate stitching
  # source for an enso_region / registry forecast (members from the frame).
  if (identical(g_pb, "country")) {
    return(invisible(got))
  }
  # enso_region caches may be assembled into registry units (a unit that is a
  # union of whole fit groups); the per-TU gid check in
  # load_smoothing_distribution() guards their content.
  if (identical(g_pb, "enso_region") && identical(e_pb, "registry")) {
    return(invisible(got))
  }
  keys <- switch(
    e_pb,
    registry = .SM_META_IDENTITY_KEYS,
    enso_region = c("partition_by", "grouping_source", "grouping_digest"),
    "partition_by"
  )
  same_val <- function(k) {
    a <- as.character(got[[k]] %||% "")
    b <- as.character(expect[[k]] %||% "")
    if (grepl("_digest$", k) && (a %in% c("NA", NA) || b %in% c("NA", NA))) {
      return(TRUE)
    }
    identical(a, b)
  }
  diff <- keys[!vapply(keys, same_val, logical(1))]
  if (length(diff)) {
    stop(
      sprintf(
        paste0(
          "[smoothing cache] %s was built under a different partition identity ",
          "(%s). Rebuild the smoothing cache for this partition, or point ",
          "filter.smoothing_input at the matching one."
        ),
        cache_dir,
        paste(
          sprintf(
            "%s: cache '%s' vs run '%s'",
            diff,
            vapply(
              diff,
              function(k) as.character(got[[k]] %||% ""),
              character(1)
            ),
            vapply(
              diff,
              function(k) as.character(expect[[k]] %||% ""),
              character(1)
            )
          ),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }
  invisible(got)
}

# Smoothing caches of exactly TU `tu` in `dir`: the glob is anchored on the
# `_Np=` token that follows the TU label in every 08c filename, so TU "grp_a"
# never matches "grp_a_x"'s cache (a bare `hist_filter_tu=<tu>_*.qs` prefix
# glob would). A legacy name WITHOUT `_Np=` is still found by the prefix glob, as
# a warned fallback. Quiet: returns character(0) when nothing matches.
.smoothing_cache_hits <- function(tu, dir) {
  hits <- Sys.glob(file.path(dir, sprintf("hist_filter_tu=%s_Np=*.qs", tu)))
  if (length(hits)) {
    return(hits)
  }
  legacy <- Sys.glob(file.path(dir, sprintf("hist_filter_tu=%s_*.qs", tu)))
  legacy <- legacy[!grepl("_Np=", basename(legacy), fixed = TRUE)]
  if (length(legacy)) {
    warning(sprintf(
      paste0(
        "[07_simulate_spatpomp] smoothing_input: TU '%s' resolved by the ",
        "legacy prefix glob (no _Np= token): %s"
      ),
      tu,
      paste(basename(legacy), collapse = ", ")
    ))
  }
  legacy
}

# Resolve the per-TU smoothing cache file (hist_filter_tu=<tu>_Np=*.qs) inside `dir`.
# Returns NULL (warn) when no file matches; takes the first match (warn) on >1.
.resolve_smoothing_cache <- function(tu, dir) {
  hits <- .smoothing_cache_hits(tu, dir)
  if (length(hits) == 0L) {
    warning(sprintf(
      "[07_simulate_spatpomp] smoothing_input: no cache for TU '%s' in %s",
      tu,
      dir
    ))
    return(NULL)
  }
  if (length(hits) > 1L) {
    warning(sprintf(
      "[07_simulate_spatpomp] smoothing_input: %d caches for TU '%s' — using %s",
      length(hits),
      tu,
      basename(hits[1L])
    ))
  }
  hits[1L]
}

# Derive the forecast handoff date for a loaded smoothing cache. Primary source
# is the trajectory's filt_months (last smoothed month + 1). When that field is
# empty/NA in the cached records (some smoothing runs do not persist it), fall back to
# the FILTER_END token in the cache filename:
#   hist_filter_tu=<tu>_Np=<Np>_sv<NN>_nrep=<N>_<burninStart>_<filterEnd>.qs
# (the `_sv<NN>` schema token sits between Np and nrep so it
# doesn't get folded into the TU capture at line ~396 below.)
# filt_months covers [start, FILTER_END), so handoff == FILTER_END.
.handoff_from_cache <- function(recs, cache_path) {
  fm <- do.call(c, lapply(recs, function(r) r$filt_months))
  if (length(fm) && !all(is.na(fm))) {
    return(as.Date(max(fm, na.rm = TRUE)) %m+% months(1L))
  }
  toks <- regmatches(
    basename(cache_path),
    gregexpr("[0-9]{8}", basename(cache_path))
  )[[1]]
  if (length(toks) >= 1L) {
    d <- tryCatch(
      as.Date(toks[length(toks)], format = "%Y%m%d"),
      error = function(e) NA
    )
    if (!is.na(d)) return(d)
  }
  NA
}

# ---------------------------------------------------------------------------
# TU list
# ---------------------------------------------------------------------------

get_tu_list <- function(art, opt) {
  if (identical(opt$simulator$scope$transmission_units, "all")) {
    real_tus <- sort(unique(stats::na.omit(art$frame$tu)))
    has_na <- any(is.na(art$frame$tu))
    if (has_na) c(real_tus, "_NA_TU_") else real_tus
  } else {
    as.character(opt$simulator$scope$transmission_units)
  }
}

# ---------------------------------------------------------------------------
# Cache path helper
# ---------------------------------------------------------------------------

# A finite simulator.incidence_cap changes the mon mean inside the filter, hence
# the particle weights, so it is part of the key: `_cap-<value>`. NULL / Inf
# (uncapped) keeps the historical name, so a cache from before the cap existed
# is reused only by an uncapped run.
.make_filter_cache_path <- function(
  tu,
  cache_dir,
  occ_stem,
  Np,
  n_reps,
  warmup_window,
  incidence_cap = NULL,
  onset_eligibility = "inactive",
  cessation_exposure = "full",
  draw_design = "per_tu_random"
) {
  cap_tok <- if (!is.null(incidence_cap) && is.finite(incidence_cap)) {
    paste0("_cap-", format(incidence_cap))
  } else {
    ""
  }
  # simulator.onset_eligibility changes the filter's onset draws: `_onset-<v>`
  # for any non-default value, the historical name under "inactive".
  ons_tok <- if (
    !is.null(onset_eligibility) && !identical(onset_eligibility, "inactive")
  ) {
    paste0("_onset-", onset_eligibility)
  } else {
    ""
  }
  # simulator.cessation_exposure changes the filter's monthly mean (mu):
  # `_cessexp-<v>` for any non-default value, the historical name under "full".
  cex_tok <- if (
    !is.null(cessation_exposure) && !identical(cessation_exposure, "full")
  ) {
    paste0("_cessexp-", cessation_exposure)
  } else {
    ""
  }
  # simulator.filter.draw_design changes the filter's draws:
  # `_draws-<v>` for any non-default value, the historical name otherwise.
  drw_tok <- if (
    !is.null(draw_design) && !identical(draw_design, "per_tu_random")
  ) {
    paste0("_draws-", draw_design)
  } else {
    ""
  }
  file.path(
    cache_dir,
    sprintf(
      "filter_traj_tu=%s_occ=%s_Np=%d_nrep=%d_warmup=%s_%s%s%s%s%s.qs",
      tu,
      occ_stem,
      Np,
      n_reps,
      gsub("-", "", warmup_window$start),
      gsub("-", "", warmup_window$end),
      cap_tok,
      ons_tok,
      cex_tok,
      drw_tok
    )
  )
}

# ---------------------------------------------------------------------------
# Smoothing-input loader (Phase A replacement)
# ---------------------------------------------------------------------------
# Load an existing smoothing distribution (one hist_filter_tu=<tu>_*.qs per
# TU) as filter_results_by_tu, and derive the per-TU forecast handoff date from
# each trajectory's filt_months (last smoothed month + 1). Records already carry
# the fields Phase B consumes (tu, rep_id, draw_d, filt_months, filter_traj,
# logLik), so no reshaping is needed.

load_smoothing_distribution <- function(
  tus,
  smoothing_input,
  forecast_window,
  handoff_month = NULL,
  partition_by = "country",
  region_members = NULL,
  expect_meta = NULL,
  region_pieces = NULL,
  tu_gids = NULL
) {
  suppressPackageStartupMessages({
    library(qs2)
  })
  # Identity of the smoothing run behind this cache dir (08c's meta.yml) vs the
  # run that reads it: refuses a cache built under another partition key or
  # another version of the partition registry (same TU label, other members).
  check_smoothing_cache_meta(smoothing_input, expect_meta)
  filter_results_by_tu <- setNames(vector("list", length(tus)), tus)
  handoff_by_tu <- setNames(vector("list", length(tus)), tus)

  # Resolve + schema-check + read one per-country handoff cache. Returns the raw
  # record list (null filter_traj entries NOT dropped — the country path keeps
  # its existing empty/n_ok accounting; region assembly drops them itself) and
  # the resolved cache path, or NULL when no cache matches the key.
  .load_member_recs <- function(member_tu) {
    cp <- .resolve_smoothing_cache(member_tu, smoothing_input)
    if (is.null(cp)) {
      return(NULL)
    }
    # Fail loudly on a stale/unversioned smoothing cache BEFORE reading it, so a
    # "rebuild the smoothing cache" instruction is the first thing a schema mismatch produces
    # (see check_smoothing_cache_schema() in spatpomp_core.R).
    check_smoothing_cache_schema(cp, member_tu)
    list(recs = qs2::qs_read(cp), cp = cp)
  }

  # Hindcast (interior-month handoff): when handoff_month is set, launch the
  # forecast from that fixed month for EVERY TU instead of the end of the
  # smoothing coverage. The smoothed state at (handoff_month - 1 month) is later
  # sliced out of filter_traj by filter_traj_to_warmup(); here we only pin the
  # handoff date and verify that month is inside each cache's filt_months.
  ovr_handoff <- if (
    !is.null(handoff_month) && nzchar(as.character(handoff_month))
  ) {
    as.Date(handoff_month)
  } else {
    as.Date(NA)
  }
  if (!is.na(ovr_handoff)) {
    message(sprintf(
      "[07_simulate_spatpomp] hindcast handoff override: launching all TUs at %s (interior smoothing month).",
      format(ovr_handoff)
    ))
  }

  for (tu in tus) {
    ml_native_cp <- NULL # set when a native per-TU cache is loaded (log label)
    if (partition_by %in% c("enso_region", "registry")) {
      # Prefer a NATIVE region smoothing cache (hist_filter_tu=<region>_*.qs),
      # produced when the bpfilter smoothing itself ran as one enso_region
      # partition (member countries coupled inside that particle filter). When
      # present, load it DIRECTLY — there is nothing to assemble. A quiet glob
      # probe (not .resolve_smoothing_cache) avoids a spurious "no cache" warning
      # on the per-country fallback below.
      native_hits <- .smoothing_cache_hits(tu, smoothing_input)
      if (length(native_hits) > 0L) {
        ml <- .load_member_recs(tu) # schema-checked read, same shape as country path
        recs <- ml$recs
        cp <- ml$cp
        ml_native_cp <- cp
      } else {
        # Fallback: compile the region IC from the EXISTING per-country smoothing
        # caches (older country-partition style; smoothing NOT re-run). Load each member
        # country's cache, drop failed-filter records, then merge paired
        # replicates into region records shaped exactly like a country record
        # (assemble_region_smoothing_records).
        members <- region_members[[tu]]
        if (is.null(members) || length(members) == 0L) {
          warning(sprintf(
            "[07_simulate_spatpomp] region TU '%s': no member countries — skipping.",
            tu
          ))
          next
        }
        # Pieces to assemble from, in order of preference: one cache per member
        # country; else (registry mode) one cache per member FIT GROUP when the
        # unit is a union of whole fit groups (region_pieces[[tu]]$groups), e.g.
        # a unit merging two enso_region smoothing partitions. Each piece is an
        # independent smoothing partition, so the joint logLik is their sum.
        quiet_hit <- function(m) {
          length(.smoothing_cache_hits(m, smoothing_input)) > 0L
        }
        alt <- region_pieces[[tu]]$groups
        if (
          !all(vapply(members, quiet_hit, logical(1))) &&
            length(alt) &&
            all(vapply(alt, quiet_hit, logical(1)))
        ) {
          message(sprintf(
            "[07_simulate_spatpomp] TU '%s': assembling from fit-group caches [%s]",
            tu,
            paste(alt, collapse = ",")
          ))
          members <- alt
        }
        member_loaded <- setNames(lapply(members, .load_member_recs), members)
        miss <- members[vapply(member_loaded, is.null, logical(1))]
        if (length(miss)) {
          warning(sprintf(
            "[07_simulate_spatpomp] region TU '%s': missing smoothing cache for member(s) %s — skipping.",
            tu,
            paste(miss, collapse = ", ")
          ))
          next
        }
        member_recs <- lapply(member_loaded, function(x) {
          Filter(function(r) !is.null(r$filter_traj), x$recs)
        })
        recs <- assemble_partition_smoothing_records(member_recs, tu)
        cp <- member_loaded[[1L]]$cp # representative, for the filename handoff fallback
      }
    } else {
      ml <- .load_member_recs(tu)
      if (is.null(ml)) {
        next
      }
      recs <- ml$recs
      cp <- ml$cp
    }
    if (length(recs) == 0L) {
      warning(sprintf(
        "[07_simulate_spatpomp] TU '%s': empty smoothing cache.",
        tu
      ))
      next
    }
    # Content check: the loaded records must cover exactly this TU's units (a
    # cache of the same label from another registry version, or pieces that do
    # not tile the unit, would otherwise hand off a wrong initial state).
    if (!is.null(tu_gids[[tu]])) {
      rg <- unique(unlist(lapply(recs, function(r) as.character(r$gids))))
      want <- unique(as.character(tu_gids[[tu]]))
      if (length(rg) && !setequal(rg, want)) {
        stop(
          sprintf(
            paste0(
              "[07_simulate_spatpomp] TU '%s': smoothing records cover %d unit(s), ",
              "the TU frame %d (%d missing from the cache, %d extra). The cache ",
              "was built for another partition."
            ),
            tu,
            length(rg),
            length(want),
            length(setdiff(want, rg)),
            length(setdiff(rg, want))
          ),
          call. = FALSE
        )
      }
    }
    n_ok <- sum(vapply(recs, function(r) !is.null(r$filter_traj), logical(1)))
    handoff <- if (!is.na(ovr_handoff)) {
      ovr_handoff
    } else {
      .handoff_from_cache(recs, cp)
    }
    if (is.na(handoff)) {
      warning(sprintf(
        "[07_simulate_spatpomp] TU '%s': cannot derive handoff (no filt_months and no FILTER_END in filename) — skipping.",
        tu
      ))
      next
    }
    # Interior-month handoff must fall inside the smoothed coverage: the IC is
    # the state at (handoff - 1 month), which must be a column of filter_traj
    # (i.e. present in filt_months). Fail fast here — before per-TU precompute —
    # rather than let filter_traj_to_warmup silently pick the wrong/last column.
    if (!is.na(ovr_handoff)) {
      fm <- do.call(c, lapply(recs, function(r) r$filt_months))
      wd <- ovr_handoff %m-% months(1L)
      if (!length(fm) || all(is.na(fm))) {
        stop(sprintf(
          "[07_simulate_spatpomp] TU '%s': handoff_month override needs per-record filt_months, but the cache has none — rebuild the smoothing cache (08a/08b/08c_smoothing_*.R).",
          tu
        ))
      }
      if (!(wd %in% as.Date(fm))) {
        stop(sprintf(
          "[07_simulate_spatpomp] TU '%s': handoff_month %s implies IC month %s, outside smoothed coverage [%s, %s].",
          tu,
          format(ovr_handoff),
          format(wd),
          format(min(as.Date(fm), na.rm = TRUE)),
          format(max(as.Date(fm), na.rm = TRUE))
        ))
      }
    }
    fws <- tryCatch(as.Date(forecast_window$start), error = function(e) NA)
    if (!is.na(fws) && handoff != fws) {
      warning(sprintf(
        "[07_simulate_spatpomp] TU '%s': handoff %s != forecast_window.start %s — using %s.",
        tu,
        format(handoff),
        forecast_window$start,
        if (!is.na(ovr_handoff)) "override" else "auto-derived"
      ))
    }
    src_label <- if (
      partition_by %in%
        c("enso_region", "registry") &&
        !identical(cp, ml_native_cp)
    ) {
      sprintf(
        "%d member caches [%s]",
        length(members),
        paste(members, collapse = ",")
      )
    } else {
      basename(cp)
    }
    message(sprintf(
      "[07_simulate_spatpomp] TU '%s': loaded %d smoothing replicates (%d ok), handoff=%s ← %s",
      tu,
      length(recs),
      n_ok,
      format(handoff),
      src_label
    ))
    filter_results_by_tu[[tu]] <- recs
    handoff_by_tu[[tu]] <- handoff
  }
  list(
    filter_results_by_tu = filter_results_by_tu,
    handoff_by_tu = handoff_by_tu
  )
}

#' .cache_handoff_label
#'
#' Dry-run label of a smoothing cache's forecast handoff month: the hindcast
#' interior-month override (`handoff_month`) when set, else the cache's own
#' end-of-coverage handoff (.handoff_from_cache()), "?" when neither resolves.
#' The cache is read either way, as 07b's two dry-run branches did.
#' @noRd
.cache_handoff_label <- function(cp, handoff_month) {
  recs <- tryCatch(qs2::qs_read(cp), error = function(e) NULL)
  if (!is.null(handoff_month) && nzchar(as.character(handoff_month))) {
    format(as.Date(handoff_month))
  } else if (!is.null(recs)) {
    h <- .handoff_from_cache(recs, cp)
    if (is.na(h)) "?" else format(h)
  } else {
    "?"
  }
}
