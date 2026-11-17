# 07b_simulator_spatpomp.R — spatPomp bpfilter cholera outbreak simulator
#
# Usage:
#   Rscript analysis/07_simulator_run.R -C <spatpomp cfg> --dry-run
#   Rscript analysis/07b_simulator_spatpomp.R -C <cfg> --Np 50 \
#     --n_replicates 5 --workers 1
#
# Details:
#   Backend "spatpomp" of 07_simulator_run.R; run directly on the SLURM path.
#   Runs AFTER phase 08: forecast/hindcast configs point
#   simulator.filter.smoothing_input at the handoff cache 08c writes.
#   Phase A filters each TU x replicate (cached), phase B resamples by logLik
#   and forecasts. Writes outdir/sim_spatpomp/<run_tag>/ in the same schema
#   as 07b_simulator_v1.R, so 07c postprocessing works unchanged.
#
# See: analysis/PIPELINE.md (run order)

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R") # drop_overlapping_adm1_from_frame() -- required by simulations.R's load_artefacts()
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/sim_io.R") # sim_io_ext()/sim_io_save() -- required by simulator_core.R's simulate_one_tu() non-dry-run write path
source("analysis/utils/dlnm.R") # build_cb_design()/build_cb_columns_from_Q() -- required whenever the loaded fit is a DLNM crossbasis (brms_cb)
source("analysis/utils/simulator_core.R") # kernel to reuse (never modified)
source("analysis/utils/spatpomp_core.R") # spatPomp builder + N-runs functions
source("analysis/utils/spatial_contagion.R") # build_partition_nbr_kernel()/nbr_from_adm1_active()/resolve_use_nbr() -- required by simulator_core.R, spatpomp_core.R and spatpomp_core_c.R whenever simulator.nbr_contagion.enable is set (modules never source each other)
source("analysis/utils/enso.R") # fetch_enso_index() etc. -- required by simulations.R's ENSO scenario builder
source("analysis/utils/population.R") # load_gid_population()/resolve_gid_population() -- required by simulator_core.R's build_location_frame() gid_pop_path branch
source("analysis/utils/composite_lp.R") # mutual dependency with population.R (README: "source both together")
source("analysis/utils/anchor.R") # resolve_base_rate_paths() -- required by simulations.R's load_artefacts() base_rate_table_paths resolution
source("analysis/utils/simulations.R") # load_artefacts, check_frame_vs_draws

# 2. CLI -----------------------------------------------------------------------

option_list <- c(
  simulator_common_options(),
  list(
    make_option(
      c("--Np"),
      default = NULL,
      type = "integer",
      help = "Override simulator.filter.Np (particles per filter run)."
    ),
    make_option(
      c("--n_draws"),
      default = NULL,
      type = "integer",
      help = "Override simulator.filter.n_draws_in_filter (legacy; ignored when fixed_draw_idx active)."
    ),
    make_option(
      c("--only_country"),
      default = NULL,
      type = "character",
      help = "Restrict the TU list to the TU holding this ISO3 (one SLURM array task = one country under filter.partition_by=country, where TU==ISO3; under enso_region / registry the country's TU is run whole)."
    ),
    make_option(
      c("--only_tu"),
      default = NULL,
      type = "character",
      help = "Restrict the TU list to one partition (TU) label, e.g. a registry group under filter.partition_by=registry."
    ),
    make_option(
      c("--replicate_offset"),
      default = NULL,
      type = "integer",
      help = "Global replicate-ID offset for this invocation. Use together with --n_replicates to run one chunk of a country's replicates per process (e.g. --n_replicates 50 --replicate_offset 100 for chunk 3 of 6x50): keeps replicate_id unique and output partitions non-colliding across concurrently-running chunks of the same TU. Default 0 (unchunked, byte-identical to prior behaviour)."
    ),
    make_option(
      c("--n_replicates_total"),
      default = NULL,
      type = "integer",
      help = "The TU's total replicate count when this invocation is one chunk of it (with --replicate_offset). The ENSO scenarios are then sampled for the whole TU and sliced to this chunk's replicate ids, so replicate r gets the same ENSO trajectory however the TU is chunked (utils/simulations.R::.enso_chunk_cols()). Unset: the chunk samples its own n (earlier behaviour)."
    )
  )
)

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = c("overwrite"),
  model = NULL
)

# Late CLI overrides (shared common flags via utils/config.R helper, then spatpomp-specific).
opt <- apply_simulator_common_cli(opt)
if (!is.null(.flag_val("--Np"))) {
  opt$simulator$filter$Np <- as.integer(.flag_val("--Np"))
}
if (!is.null(.flag_val("--n_draws"))) {
  opt$simulator$filter$n_draws_in_filter <- as.integer(.flag_val("--n_draws"))
}
if (!is.null(.flag_val("--only_country"))) {
  opt$only_country <- .flag_val("--only_country")
}
if (!is.null(.flag_val("--only_tu"))) {
  opt$only_tu <- .flag_val("--only_tu")
}
# NULL (flag absent) means "unchunked" — preserve prior seed/output-path
# behaviour exactly; only a passed --replicate_offset flips chunked-mode on.
opt$replicate_offset <- if (!is.null(.flag_val("--replicate_offset"))) {
  as.integer(.flag_val("--replicate_offset"))
} else {
  0L
}
opt$.replicate_chunked <- !is.null(.flag_val("--replicate_offset"))
opt$replicate_total <- if (!is.null(.flag_val("--n_replicates_total"))) {
  as.integer(.flag_val("--n_replicates_total"))
} else {
  NULL
}

require_simulator_block(opt)

print_options(opt)

# 3. Validation ----------------------------------------------------------------

validate_spatpomp_options <- function(opt) {
  sim <- opt$simulator
  flt <- sim$filter %||% list()

  required_paths <- c(
    occurrence_fit = sim$occurrence_fit,
    climate_panel = sim$climate_panel
  )
  for (nm in names(required_paths)) {
    p <- required_paths[[nm]]
    if (is.null(p) || !nzchar(p)) {
      stop(sprintf("simulator.%s is required.", nm))
    }
    if (!file.exists(p)) {
      stop(sprintf("simulator.%s not found: %s", nm, p))
    }
  }

  # Smoothing-input mode: Phase A is skipped, so warmup_window is irrelevant.
  # Require the cache directory to exist and a valid forecast_window; do not
  # validate (or require) the warmup_window.
  smoothing_input <- flt$smoothing_input
  if (!is.null(smoothing_input) && nzchar(smoothing_input)) {
    if (!dir.exists(smoothing_input)) {
      stop(sprintf(
        "simulator.filter.smoothing_input not found: %s",
        smoothing_input
      ))
    }
    fw <- flt$forecast_window
    if (!is.null(fw)) {
      fs <- tryCatch(as.Date(fw$start), error = function(e) NA)
      fe <- tryCatch(as.Date(fw$end), error = function(e) NA)
      if (is.na(fs) || is.na(fe) || fe <= fs) {
        stop("simulator.filter.forecast_window: end must be after start.")
      }
    }
    return(invisible(opt))
  }

  # Filter windows
  ww <- flt$warmup_window
  fw <- flt$forecast_window
  if (!is.null(ww)) {
    ws <- tryCatch(as.Date(ww$start), error = function(e) NA)
    we <- tryCatch(as.Date(ww$end), error = function(e) NA)
    if (is.na(ws) || is.na(we) || we <= ws) {
      stop("simulator.filter.warmup_window: end must be after start.")
    }
  }
  if (!is.null(fw)) {
    fs <- tryCatch(as.Date(fw$start), error = function(e) NA)
    fe <- tryCatch(as.Date(fw$end), error = function(e) NA)
    if (is.na(fs) || is.na(fe) || fe <= fs) {
      stop("simulator.filter.forecast_window: end must be after start.")
    }
    if (!is.null(ww)) {
      we <- as.Date(ww$end)
      if (!is.na(fs) && abs(as.numeric(fs - we)) > 2L) {
        warning(
          sprintf(
            "filter.forecast_window.start (%s) != warmup_window.end (%s) — ",
            fw$start,
            ww$end
          ),
          "handoff may be misaligned."
        )
      }
    }
  }

  Np <- as.integer(flt$Np %||% 2000L)
  if (Np < 2L) {
    stop("simulator.filter.Np must be >= 2.")
  }

  meas <- flt$measurement %||% "occurrence_duration_mon_NB"
  # .SP_ALL_MEAS is the single source of truth (spatpomp_core.R); a local copy
  # here would drift the moment a measurement is added on the kernel side.
  if (!meas %in% .SP_ALL_MEAS) {
    stop(sprintf(
      "simulator.filter.measurement '%s' not recognised. Valid: %s",
      meas,
      paste(.SP_ALL_MEAS, collapse = ", ")
    ))
  }
  if (.is_mon_meas(meas)) {
    # Hard requirements, not warnings: the mon measurement scores cases at every
    # ACTIVE month, and `active`/`mso` exist only because the duration model
    # maintains them. Without either artefact there is nothing to degrade to.
    if (!isTRUE(sim$grouping$enable)) {
      if (!nzchar(sim$monthly_fit %||% "")) {
        stop(
          "simulator.filter.measurement = '",
          meas,
          "' requires ",
          "simulator.monthly_fit."
        )
      }
      if (!nzchar(sim$duration_fit %||% "")) {
        stop(
          "simulator.filter.measurement = '",
          meas,
          "' requires ",
          "simulator.duration_fit (`active`/`mso` are duration-model latent ",
          "states)."
        )
      }
    } else {
      if (!isTRUE(sim$grouping$monthly_enable)) {
        stop(
          "simulator.filter.measurement = '",
          meas,
          "' requires ",
          "simulator.grouping.monthly_enable in grouped mode."
        )
      }
      if (!isTRUE(sim$grouping$duration_enable)) {
        stop(
          "simulator.filter.measurement = '",
          meas,
          "' requires ",
          "simulator.grouping.duration_enable (`active`/`mso` are ",
          "duration-model latent states)."
        )
      }
    }
  }

  invisible(opt)
}

# 4. Csnippet builder ----------------------------------------------------------
# The only spatPomp builder.

.spdir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "analysis")
source(file.path(.spdir, "utils", "spatpomp_core_c.R"))

validate_spatpomp_options(opt)

# 5. Output paths --------------------------------------------------------------
# Under sim_spatpomp/<run_tag>/.

if (is.null(opt$simulator$run_tag)) {
  opt$simulator$run_tag <- derive_sim_run_tag("spatpomp")
  message("[07_simulate_spatpomp] derived run_tag = ", opt$simulator$run_tag)
}

# Delegates to make_run_file(opt, kind = "spatpomp"...)/utils.R::make_spatpomp_run_file(),
# which mirrors this exact (flat, non-run_id-scoped) layout — see that
# function's header comment for why run_id is deliberately NOT honored here.
run_dir <- make_run_file(opt, kind = "spatpomp")
parts_dir <- make_run_file(opt, kind = "spatpomp_parts_dir")
diag_dir <- make_run_file(opt, kind = "spatpomp_diag_dir")
# Filter cache:
#   generated_data/sim_spatpomp/<run_tag>/cache/
#     filter_traj_tu=<tu>_occ=<occ_stem>_Np=<Np>_nrep=<N>_warmup=<s>_<e>.qs
#   Controlled by filter.rerun_filter (default FALSE = use cache if present).
cache_dir <- make_run_file(opt, kind = "spatpomp_cache_dir")

if (!isTRUE(opt$dry_run)) {
  dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
}

message("[07_simulate_spatpomp] run_dir   = ", run_dir)

# 6. Load artefacts ------------------------------------------------------------

# In smoothing-input mode the effective per-TU forecast window starts at the
# auto-derived handoff (from the smoothing cache), which can be earlier than
# sim$horizon$start. The shared ENSO scenario matrices built inside
# load_artefacts() must span from that handoff; otherwise they end up shorter
# than each TU's climate panel and get silently recycled by the per-replicate
# ENSO override in simulate_one_tu() (simulator_core.R). Pre-scan the cache
# filenames (cheap — no payload read) for the earliest handoff and pass it down
# as enso_sim_months_start so the scenario window extends back to cover it.
.enso_start_override <- NULL
{
  .si <- opt$simulator$filter$smoothing_input %||% NULL
  if (!is.null(.si) && nzchar(.si) && dir.exists(.si)) {
    .hits <- Sys.glob(file.path(.si, "hist_filter_tu=*_Np=*.qs"))
    .hoffs <- as.Date(
      vapply(
        .hits,
        function(f) {
          toks <- regmatches(basename(f), gregexpr("[0-9]{8}", basename(f)))[[
            1
          ]]
          if (length(toks) >= 1L) toks[length(toks)] else NA_character_
        },
        character(1)
      ),
      format = "%Y%m%d"
    )
    .hoffs <- .hoffs[!is.na(.hoffs)]
    if (length(.hoffs)) {
      .enso_start_override <- min(.hoffs)
      message(sprintf(
        "[07_simulate_spatpomp] smoothing-input: earliest handoff across %d cache(s) = %s (ENSO scenario window extends back to cover it).",
        length(.hoffs),
        format(.enso_start_override)
      ))
      if (length(unique(.hoffs)) > 1L) {
        warning(
          "[07_simulate_spatpomp] smoothing caches have NON-UNIFORM handoff dates (",
          paste(sort(unique(format(.hoffs))), collapse = ", "),
          "). The shared ENSO scenario matrices are positionally aligned to a single ",
          "window, so per-TU forecast windows that differ in length will still be ",
          "mis-aligned. Calendar-date alignment in the kernel would be required to fix ",
          "this fully — results for TUs whose handoff != the earliest may be wrong."
        )
      }
    }
  }
}

art <- if (isTRUE(opt$dry_run)) {
  NULL
} else {
  message("[07_simulate_spatpomp] loading artefacts...")
  load_artefacts(opt, enso_sim_months_start = .enso_start_override)
}

if (!is.null(art) && all(is.na(art$frame$tu))) {
  message(
    "[07_simulate_spatpomp] every location has NA TU — collapsing to _ALL_."
  )
  art$frame$tu <- "_ALL_"
}

# Grouped mode: the Csnippet kernel reads
# only precomputed per-unit covariate columns (lp_static, ttlo/spill coefs); it
# performs no runtime LP evaluation, so multi-group routing happens entirely at
# R build time inside build_spatpomp_tu_c(grouping=) — no codegen change needed.

# Build the grouped routing list passed to build_spatpomp_tu_c / forecast_from_filtered.
grp_arg_08 <- if (!is.null(art) && isTRUE(art$grouped)) {
  list(
    draws_occ_by_group = art$draws_occ_by_group,
    draws_dur_by_group = art$draws_dur_by_group, # NULL unless grouping.duration_enable
    draws_mon_by_group = art$draws_mon_by_group, # NULL unless grouping.monthly_enable
    cb_meta_occ_by_group = art$cb_meta_occ_by_group,
    cb_meta_dur_by_group = art$cb_meta_dur_by_group, # NULL unless grouping.duration_enable
    cb_meta_mon_by_group = art$cb_meta_mon_by_group
  )
} else {
  NULL
}

# Neighbour-contagion build-time guard, resolved once here (art is loaded) and
# threaded into both the native Phase-A build_spatpomp_tu_c() filter build and
# the forecast_from_filtered() handoff. On only when the config flag is set AND
# the (per-group) occurrence fit carries an nbr_lag1_r* FE; FALSE ⇒ every nbr
# branch in the C kernel and the R forecast kernel is inert, byte-identical to
# pre-change. See resolve_use_nbr() (spatial_contagion.R).
nbr_cfg_08 <- nbr_cfg_from_sim(opt$simulator)
nbr_spec_08 <- if (!is.null(art)) {
  resolve_nbr_spec(nbr_cfg_08, art$draws_occ, grp_arg_08)
} else {
  NULL
}
use_nbr_08 <- isTRUE(nbr_spec_08$use)

# 7. Resolve filter parameters -------------------------------------------------

flt <- opt$simulator$filter %||% list()
# `%||% TRUE` disagrees with get_default_options()$simulator$filter$enable
# (FALSE, the effective default: simulations.R's isTRUE() is the only live
# reader). Left as is: filter_enabled is read nowhere in this script, and
# get_default_options() fills `enable` before this line for every config
# resolved by make_options_from_config().
filter_enabled <- isTRUE(flt$enable %||% TRUE)
Np <- as.integer(flt$Np %||% 2000L)
rerun_filter <- isTRUE(flt$rerun_filter %||% FALSE)
n_draws_in_filter <- as.integer(flt$n_draws_in_filter %||% 1L) # 1 with fixed-draw design
n_replicates <- as.integer(opt$simulator$n_replicates %||% 100L)
warmup_window <- flt$warmup_window %||%
  list(start = "2021-07-01", end = "2023-05-01")
forecast_window <- flt$forecast_window %||%
  list(
    start = warmup_window$end,
    end = opt$simulator$horizon$end %||%
      get_default_options()$simulator$horizon$end
  )
seed_root <- as.integer(opt$simulator$seed %||% 20260101L)
# Measurement-side keys: the shared resolve_filter_cfg() (spatpomp_core.R),
# the same function the smoothing driver uses.
fc_08 <- resolve_filter_cfg(flt)
measurement <- fc_08$measurement
delta_fn <- fc_08$delta_fn
delta_fp <- fc_08$delta_fp
case_fp_rate <- fc_08$case_fp_rate
case_fp_mean <- fc_08$case_fp_mean
# Monthly incidence cap: monthly cases per unit <= incidence_cap x pop, on the
# mean in the filter kernel and on the mean and the drawn count in
# the forecast kernel (simulate_one_tu()). Part of the filter-cache key.
incidence_cap <- incidence_cap_of(opt$simulator)
# Onset eligibility (simulator.onset_eligibility): the onset
# gate of the filter kernel and the forecast kernel alike. A non-default value
# is part of the filter-cache key and of the smoothing-cache identity.
onset_eligibility <- onset_eligibility_of(opt$simulator)
# Cessation-month exposure (simulator.cessation_exposure): the
# monthly mean of a ceasing month, in the filter and forecast kernels alike. A
# non-default value is part of the filter-cache key and of the smoothing-cache
# identity.
cessation_exposure <- cessation_exposure_of(opt$simulator)

# Top-level partition unit. Must match the smoothing run that produced the
# handoff caches: partition_by=country keys every artefact by ISO3 country;
# partition_by=enso_region keys the FORWARD/hindcast run by the ENSO fit group;
# partition_by=registry keys it by the group of an arbitrary country-level
# registry (filter.partition_source), independent of the fit grouping, so one TU
# may mix fit groups (each unit keeps its own group's draws, cb_meta and lag).
# For enso_region / registry a native per-TU smoothing cache is used when
# present, else the IC is assembled from per-country caches
# (load_smoothing_distribution()). transmission_unit (default) keeps the HMM-TU
# key already carried on the frame and does not re-key.
partition_by <- flt$partition_by %||% "transmission_unit"
partition_source <- as.character(flt$partition_source %||% "")
if (identical(partition_by, "registry")) {
  # Fails here, before any artefact work, on a missing / path-unsafe registry.
  .pmap_07b <- load_partition_map(partition_source)
  message(sprintf(
    "[07_simulate_spatpomp] partition_by=registry: source '%s' (%d countries, %d TUs, digest %s)",
    partition_source,
    length(.pmap_07b),
    length(unique(.pmap_07b)),
    registry_digest(partition_source)
  ))
}

# The partition_by -> frame-column map and its presence check live in
# spatpomp_core.R (partition_frame_column() / assert_partition_column()), shared
# with the historical-smoothing phase 08. "transmission_unit" is a no-op
# (frame$tu is already the HMM cluster) and is handled by the absence of a re-key
# below, not by that map.

if (identical(partition_by, "registry") && !is.null(art)) {
  art$frame <- stamp_partition_column(art$frame, "registry", partition_source)
}

if (
  partition_by %in% c("country", "enso_region", "registry") && !is.null(art)
) {
  # Re-key the forecast frame's `tu` to the partition column so get_tu_list()/
  # precompute_tu_inputs() enumerate and slice by that key:
  #   country     → each ISO3 is its own TU (mirrors the per-country smoothing
  #                 run; the leave-one-country-out o_tu term resolves to 0 within
  #                 a single-country partition).
  #   enso_region → each ENSO group is one TU spanning all member-country ADM1
  #                 units, so forecast_from_filtered()'s neighbour-contagion
  #                 kernel (built over frame_tu's ADM1 centroids) couples across
  #                 country borders within the region. The handoff IC is compiled
  #                 from the per-country smoothing caches in
  #                 load_smoothing_distribution() (smoothing is NOT re-run).
  #   registry    → each group of filter.partition_source is one TU; units keep
  #                 their fit `group`, so a TU may mix fit groups.
  col <- assert_partition_column(art$frame, partition_by)
  art$frame$tu <- as.character(art$frame[[col]])
  message(sprintf(
    "[07_simulate_spatpomp] partition_by=%s: %d %s TU(s).",
    partition_by,
    length(unique(stats::na.omit(art$frame$tu))),
    switch(
      partition_by,
      enso_region = "region",
      registry = "registry",
      "country"
    )
  ))
}

# Pre-flight on the (re-keyed) frame: fit-group coverage, mixed-TU duration
# basehaz coherence, contagion partition refining the TU partition.
if (!is.null(art)) {
  check_partition_vs_grouping(art$frame, partition_by, grp_arg_08, nbr_spec_08)
}

# Run identity record (partition + fit grouping, registry names + content
# digests), next to the run's parts: what a reader needs to know which registry
# version keyed these TUs. Same fields as the smoothing cache meta.yml.
if (!isTRUE(opt$dry_run)) {
  .pm <- smoothing_cache_meta(
    partition_by,
    partition_source,
    grouping_source = if (isTRUE(opt$simulator$grouping$enable)) {
      opt$simulator$grouping$source
    } else {
      ""
    },
    fits_dir = opt$simulator$grouping$fits_dir %||% ""
  )
  .pm$run_tag <- opt$simulator$run_tag %||% ""
  # The arm's launch (= its smoothing handoff month; sim_io.R::assert_arm_launch()).
  .pm$horizon_start <- format(as.Date(as.character(
    opt$simulator$horizon$start
  )))
  # ENSO past the CFSv2 span (enso.R::enso_continuation_provenance(): method, and
  # for analog the library path + content digest).
  .pm$enso_continuation <- enso_continuation_provenance(
    opt$simulator$enso$continuation
  )
  # Counterfactual analog cache identity: only for a library-
  # mode cache (lenssen_cache_provenance() is NULL for a legacy one), checked by
  # 09 (sim_io.R::assert_counterfactual_cache_identity()).
  if (
    identical(
      (opt$simulator$enso$counterfactual %||% list())[["scenario"]],
      "lenssen"
    )
  ) {
    .pm$enso_counterfactual <- lenssen_cache_provenance(
      lenssen_counterfactual_args(
        opt$simulator$enso$counterfactual,
        opt$simulator$horizon$start
      )$cache_path
    )
  }
  .pm$Np <- NULL
  .pm$nrep <- NULL
  .pm$schema <- NULL
  # Population year: recorded only when not the legacy 2020,
  # so a default run's partition_meta.yml is unchanged.
  .pm <- .cache_meta_add_optional(
    .pm,
    list(
      pop_year = sim_pop_year_spec(opt$simulator),
      onset_eligibility = onset_eligibility,
      cessation_exposure = cessation_exposure
    )
  )
  # Simulation environment (table S13; provenance.R::sim_software_versions()):
  # informational, read by 09's report_software, checked by nothing.
  .pm$software <- sim_software_versions()
  yaml::write_yaml(.pm, file.path(run_dir, "partition_meta.yml"))
}

# Smoothing-input mode: load an existing smoothing distribution (cache dir)
# and skip Phase A entirely.
smoothing_input <- flt$smoothing_input %||% NULL
use_smoothing <- !is.null(smoothing_input) && nzchar(smoothing_input)

# param_draw_mode "joint": forecast replicate r starts from
# smoothing replicate j(r) = ((r - 1) mod nrep) + 1 and runs the posterior draw
# that replicate was reconstructed under, so state and parameters are one joint
# posterior draw. Needs a smoothing cache built with
# simulator.filter.draw_design shared_cycle under this run's max_draws (checked
# against its meta.yml below), and the uniform handoff (the joint pairing
# replaces it; a shuffled / loglik pairing would contradict the design).
joint_draws <- identical(opt$simulator$forecast$param_draw_mode, "joint")
joint_n_smoothing <- NULL
if (joint_draws) {
  if (!use_smoothing) {
    stop(
      "simulator.forecast.param_draw_mode joint needs filter.smoothing_input (a shared_cycle ",
      "smoothing cache).",
      call. = FALSE
    )
  }
  .hw <- opt$simulator$forecast$handoff_weighting %||% "uniform"
  if (!identical(.hw, "uniform")) {
    stop(
      "simulator.forecast.param_draw_mode joint pairs each forecast replicate with its smoothing ",
      "replicate itself: set simulator.forecast.handoff_weighting to uniform, not '",
      .hw,
      "'.",
      call. = FALSE
    )
  }
  .mp <- smoothing_cache_meta_path(smoothing_input)
  if (file.exists(.mp)) {
    joint_n_smoothing <- as.integer(yaml::read_yaml(.mp)$nrep %||% NA)
  }
  # A dry run reports it and goes on (the cache may not be compiled yet); a run stops.
  .jm <- if (!file.exists(.mp)) {
    paste0(smoothing_input, " has no meta.yml, so its draw design is unknown.")
  } else if (is.na(joint_n_smoothing) || joint_n_smoothing < 1L) {
    paste0(
      "meta.yml of ",
      smoothing_input,
      " carries no replicate count (nrep)."
    )
  } else {
    NULL
  }
  if (!is.null(.jm)) {
    if (!isTRUE(opt$dry_run)) {
      stop("simulator.forecast.param_draw_mode joint: ", .jm, call. = FALSE)
    }
    message("[07_simulate_spatpomp] [DRY RUN] param_draw_mode joint: ", .jm)
  }
}

# Occurrence model stem for cache filenames
occ_stem <- tools::file_path_sans_ext(basename(
  opt$simulator$occurrence_fit %||% "occ"
))

message(sprintf(
  "[07_simulate_spatpomp] filter: Np=%d  n_replicates=%d  rerun_filter=%s",
  Np,
  n_replicates,
  rerun_filter
))
if (use_smoothing) {
  message(sprintf(
    "[07_simulate_spatpomp]         smoothing_input=%s (Phase A skipped)  forecast=%s – %s",
    smoothing_input,
    forecast_window$start,
    forecast_window$end
  ))
} else {
  message(sprintf(
    "[07_simulate_spatpomp]         warmup=%s – %s  forecast=%s – %s",
    warmup_window$start,
    warmup_window$end,
    forecast_window$start,
    forecast_window$end
  ))
}

if (!is.null(art)) {
  tus <- get_tu_list(art, opt)
  message(sprintf(
    "[07_simulate_spatpomp] %d TU(s): %s",
    length(tus),
    paste(tus, collapse = ", ")
  ))
} else {
  tu_scope <- opt$simulator$scope$transmission_units
  if (identical(partition_by, "registry")) {
    # Registry TUs of this config's scope.countries, from the partition registry
    # (artefacts are not loaded in dry-run).
    scope_iso <- opt$simulator$scope$countries %||% "all"
    pm <- .pmap_07b
    if (!identical(scope_iso, "all")) {
      pm <- pm[
        names(pm) %in%
          trimws(strsplit(paste(scope_iso, collapse = ","), ",")[[1L]])
      ]
    }
    tus <- sort(unique(unname(pm)))
  } else if (identical(partition_by, "enso_region")) {
    # Region TUs are the ENSO groups of this config's scope.countries, read from
    # the grouping registry — the per-country smoothing caches carry no region
    # key, so (unlike the country branch below) cache-filename globbing does not
    # recover them. Artefacts are not loaded in dry-run, so the registry is the
    # only source of the country→group map here.
    src <- opt$simulator$grouping$source
    reg <- tryCatch(load_grouping(src), error = function(e) NULL)
    scope_iso <- opt$simulator$scope$countries %||% "all"
    if (!is.null(reg)) {
      if (!identical(scope_iso, "all")) {
        reg <- reg[
          reg$iso3 %in%
            trimws(strsplit(paste(scope_iso, collapse = ","), ",")[[1L]]),
          ,
          drop = FALSE
        ]
      }
      tus <- sort(unique(reg$group))
    } else {
      tus <- character(0)
    }
  } else if (
    use_smoothing &&
      identical(partition_by, "country") &&
      dir.exists(smoothing_input)
  ) {
    # Artefacts are not loaded in dry-run, so derive the country TU list from the
    # handoff caches themselves (hist_filter_tu=<ISO3>_Np=...qs).
    hits <- Sys.glob(file.path(smoothing_input, "hist_filter_tu=*_Np=*.qs"))
    tus <- sort(unique(sub(
      "^hist_filter_tu=(.*?)_Np=.*$",
      "\\1",
      basename(hits)
    )))
    if (length(tus) == 0L) tus <- c("hmmW", "hmmCE", "hmmS", "_NA_TU_")
  } else if (!identical(tu_scope, "all") && !is.null(tu_scope)) {
    tus <- as.character(tu_scope)
  } else {
    tus <- c("hmmW", "hmmCE", "hmmS", "_NA_TU_")
  }
  message(sprintf(
    "[07_simulate_spatpomp] [DRY] TU scope: %s",
    paste(tus, collapse = ", ")
  ))
}

# `tus_all` = the full group's TU list BEFORE any --only_country restriction.
# Kept around so per-TU seed striding (run_forecast_phase) stays anchored to a
# TU's stable position within the whole group even when this process only
# handles one of them — otherwise every --only_country invocation would see
# match(tu, tus) == 1 and collide on the same seed range. See run_forecast_phase().
tus_all <- tus
if (!is.null(opt$only_country) && nzchar(opt$only_country)) {
  # Under enso_region / registry the country's whole TU is run (a TU couples
  # its member countries); under country TU == ISO3.
  want <- opt$only_country
  if (partition_by %in% c("enso_region", "registry")) {
    iso_map <- if (identical(partition_by, "registry")) {
      .pmap_07b
    } else {
      reg <- load_grouping(opt$simulator$grouping$source)
      stats::setNames(reg$group, reg$iso3)
    }
    want <- unname(iso_map[opt$only_country])
    if (is.na(want)) {
      stop(
        sprintf(
          "--only_country '%s' is not in the %s partition registry.",
          opt$only_country,
          partition_by
        ),
        call. = FALSE
      )
    }
  }
  tus <- intersect(tus, want)
  if (length(tus) == 0L) {
    stop(
      sprintf(
        "--only_country '%s' not found among this config's TU scope (%s). ",
        opt$only_country,
        paste(tus_all, collapse = ", ")
      ),
      "Its TU is not in scope."
    )
  }
  message(sprintf(
    "[07_simulate_spatpomp] --only_country %s restricts run to TU '%s'.",
    opt$only_country,
    tus
  ))
}
if (!is.null(opt$only_tu) && nzchar(opt$only_tu)) {
  if (!opt$only_tu %in% tus) {
    stop(
      sprintf(
        "--only_tu '%s' not found among this config's TU scope (%s).",
        opt$only_tu,
        paste(tus, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  tus <- opt$only_tu
  message(sprintf(
    "[07_simulate_spatpomp] --only_tu restricts run to TU '%s'.",
    tus
  ))
}

# 8. Pre-compute per-TU inputs -------------------------------------------------
# Sequential, done once per TU.

precompute_tu_inputs <- function(tu, art, warmup_window, forecast_window, opt) {
  frame_tu <- if (identical(tu, "_NA_TU_")) {
    art$frame %>% dplyr::filter(is.na(.data$tu))
  } else {
    art$frame %>% dplyr::filter(.data$tu == .env$tu)
  }
  if (nrow(frame_tu) == 0L) {
    return(NULL)
  }

  scenario <- opt$simulator$scenario$climate %||% "factual"

  # Full panel covering both warmup tail and forecast horizon
  panel_tu_full <- build_climate_panel_for_tu(
    art$panel,
    frame_tu,
    horizon = forecast_window,
    scenario = scenario,
    end_mode = "closed" # forward horizon [start, end] (utils/base.R)
  )
  panel_tu_warmup <- panel_tu_full %>%
    dplyr::filter(t_start < as.Date(warmup_window$end))

  ep_opt <- list(
    baseline = opt$baseline %||% "1990_2019",
    episodes = opt$episodes
  )
  obs_data_warmup <- if (!is.null(art$surv_panel)) {
    obs_panel_for_filter(
      art$surv_panel,
      fc_08,
      make_monthly_dat_file(ep_opt),
      surv_path = make_surv_dat_file(ep_opt)
    )
  } else {
    NULL
  }

  list(
    tu = tu,
    frame_tu = frame_tu,
    panel_tu_warmup = panel_tu_warmup,
    panel_tu_full = panel_tu_full,
    obs_data_warmup = obs_data_warmup
  )
}

# 9. Phase A: filter -----------------------------------------------------------
# N bpfilter runs per TU, flat furrr over TU × replicate.

# Pipeline (N-independent-runs design):
#   1. Load artefacts (load_artefacts from utils/simulations.R).
#   2. Phase A — Filter (parallel over TU × replicate):
#        For each replicate i of each TU:
#          build_spatpomp_tu_c(fixed_draw_idx = d_i) → bpfilter(filter_traj=TRUE)
#        → filter_traj array (9U × 1 × M) + logLik per replicate, saved to cache.
#   3. Phase B — Resample + Forecast (sequential per TU):
#        filter_traj_to_warmup → resample_by_loglik → forecast_from_filtered.
run_filter_phase <- function(
  tus,
  art,
  Np,
  n_replicates,
  n_draws_in_filter,
  warmup_window,
  tu_precomp,
  cache_dir,
  occ_stem,
  rerun_filter,
  seed_root,
  dry_run = FALSE
) {
  suppressPackageStartupMessages({
    library(qs2)
  })

  filter_results_by_tu <- vector("list", length(tus))
  names(filter_results_by_tu) <- tus

  # Check per-TU caches
  uncached_tus <- character(0)
  for (tu in tus) {
    cp <- .make_filter_cache_path(
      tu,
      cache_dir,
      occ_stem,
      Np,
      n_replicates,
      warmup_window,
      incidence_cap,
      onset_eligibility,
      cessation_exposure,
      smoothing_draw_design_of(flt)
    )
    if (file.exists(cp) && !rerun_filter) {
      message(sprintf(
        "[10] TU '%s': loading filter cache (%d replicates)...",
        tu,
        n_replicates
      ))
      filter_results_by_tu[[tu]] <- qs2::qs_read(cp)
      message(sprintf("[10] TU '%s': cache loaded.", tu))
    } else {
      # A cache built before simulator.incidence_cap / onset_eligibility /
      # cessation_exposure existed (no `_cap-` / `_onset-` / `_cessexp-` token)
      # ran an uncapped, inactive-gate, full-exposure filter: it is not reused
      # under a cap, the season gate or half cessation exposure.
      legacy <- .make_filter_cache_path(
        tu,
        cache_dir,
        occ_stem,
        Np,
        n_replicates,
        warmup_window
      )
      if (!rerun_filter && legacy != cp && file.exists(legacy)) {
        message(sprintf(
          paste0(
            "[10] TU '%s': ignoring filter cache %s -- it was built ",
            "without simulator.incidence_cap = %s / onset_eligibility = %s / ",
            "cessation_exposure = %s; re-running the filter."
          ),
          tu,
          basename(legacy),
          format(incidence_cap),
          onset_eligibility,
          cessation_exposure
        ))
      }
      uncached_tus <- c(uncached_tus, tu)
    }
  }

  if (dry_run) {
    for (tu in tus) {
      pc <- tu_precomp[[tu]]
      n_units <- if (!is.null(pc)) nrow(pc$frame_tu) else NA
      message(sprintf(
        "[10] [DRY] TU '%s': %d units  Np=%d  n_replicates=%d  warmup=%s – %s",
        tu,
        n_units,
        Np,
        n_replicates,
        warmup_window$start,
        warmup_window$end
      ))
      cached <- tu %in% setdiff(tus, uncached_tus)
      message(sprintf(
        "[10] [DRY] TU '%s': cache %s",
        tu,
        if (cached) "HIT (will load)" else "MISS (will run)"
      ))
    }
    return(invisible(NULL))
  }

  if (length(uncached_tus) == 0L) {
    message("[10] All TUs loaded from cache.")
    return(filter_results_by_tu)
  }

  message(sprintf(
    "[10] Running filter for %d TU(s): %s",
    length(uncached_tus),
    paste(uncached_tus, collapse = ", ")
  ))

  # Pre-assign draw indices per TU (seeded for reproducibility).
  # In grouped mode each group has its own S; use the smallest S so every
  # sampled index is valid for all groups (build_spatpomp_tu_c clamps per-group
  # defensively, but sampling within the true range is cleaner).
  .max_draw_s <- if (isTRUE(art$grouped)) {
    min(vapply(art$draws_occ_by_group, function(d) d$S, integer(1)))
  } else {
    art$draws_occ$S
  }
  # simulator.filter.draw_design: "per_tu_random" = a
  # seeded per-TU sample; "shared_cycle" = replicate j
  # runs draw j in every TU (spatpomp_core.R::smoothing_draw_idx()).
  draw_assignments_by_tu <- setNames(
    lapply(uncached_tus, function(tu) {
      smoothing_draw_idx(
        smoothing_draw_design_of(flt),
        seed_root,
        tu,
        tus,
        n_replicates,
        .max_draw_s
      )
    }),
    uncached_tus
  )

  # Build flat task table: (tu, rep_id, draw_d, task_seed) for all uncached TU × rep
  task_rows <- do.call(
    rbind,
    lapply(uncached_tus, function(tu) {
      data.frame(
        tu = tu,
        rep_id = seq_len(n_replicates),
        draw_d = draw_assignments_by_tu[[tu]],
        stringsAsFactors = FALSE
      )
    })
  )
  task_rows$task_seed <- seed_root + seq_len(nrow(task_rows))

  n_tasks <- nrow(task_rows)
  message(sprintf(
    "[10] Dispatching %d filter tasks (%d workers)...",
    n_tasks,
    getOption("mc.cores", 1L)
  ))

  # Parallel filter runs
  raw_results <- furrr::future_map(
    seq_len(n_tasks),
    function(i) {
      tu_i <- task_rows$tu[i]
      rep_i <- task_rows$rep_id[i]
      draw_i <- task_rows$draw_d[i]
      seed_i <- task_rows$task_seed[i]

      pc <- tu_precomp[[tu_i]]
      if (is.null(pc)) {
        return(NULL)
      }

      tryCatch(
        {
          build_args_i <- list(
            pc$frame_tu,
            pc$panel_tu_warmup,
            art$draws_occ,
            window = warmup_window,
            n_draws_in_filter = n_draws_in_filter,
            fixed_draw_idx = draw_i,
            obs_data = pc$obs_data_warmup,
            cb_meta = art$cb_meta_occ,
            grouping = grp_arg_08,
            measurement = measurement,
            # Native cold-start Phase A: same partition_by=country o_tu collapse
            # as the smoothing bpfilter's rprocess; force it to 1 for consistency
            # (build_spatpomp_tu_c skips .fold_observed_otu_into_lp automatically
            # when force_o_tu is set).
            force_o_tu = identical(partition_by, "country"),
            # Suppresses o_ct/o_tu state tracking AND LP contribution,
            # independent of force_o_tu. See config_defaults.yml's
            # simulator.drop_country_tu_spillover.
            drop_country_tu_spillover = isTRUE(
              opt$simulator$drop_country_tu_spillover
            ),
            # Neighbour-contagion: build the in-kernel WNBR matvec. Guard resolved
            # once above (use_nbr_08). Gated ONLY by use_nbr -- independent of
            # force_o_tu / drop_country_tu_spillover.
            use_nbr = use_nbr_08,
            nbr_rho_km = nbr_cfg_08$rho_km %||% 500L,
            nbr_gadm_path = nbr_cfg_08$gadm_path %||%
              "data/geodata/gadm_valid.gpkg",
            nbr_cache_dir = nbr_cfg_08$cache_dir %||% "generated_data/",
            nbr_spec = nbr_spec_08
          )
          # Duration draws: build_spatpomp_tu_c() wires them end to end (per-draw
          # subsetting, precompute_dur_* covariates, "draws_dur_sub" attr).
          build_args_i$draws_dur <- art$draws_dur
          # mon (monthly-case) draws/cb_meta, mirroring draws_dur above; ignored
          # unless `measurement` is one of .SP_MON_MEAS
          # (see .check_measurement()/use_mon inside build_spatpomp_tu_c()).
          build_args_i$draws_mon <- art$draws_mon
          build_args_i$cb_meta_mon <- art$cb_meta_mon
          # Detection-error rates: consumed only by `*_state_detect` measurements,
          # inert otherwise.
          build_args_i$delta_fn <- delta_fn
          build_args_i$delta_fp <- delta_fp
          # Inactive-state case mixture: consumed only by the `_actn` measurements.
          build_args_i$case_fp_rate <- case_fp_rate
          build_args_i$case_fp_mean <- case_fp_mean
          build_args_i$mon_case_col <- fc_08$mon_case_col
          build_args_i$covariate_feed <- fc_08$covariate_feed
          # Monthly incidence cap (the mon_mu_cap covariate).
          build_args_i$incidence_cap <- incidence_cap
          # Onset gate (simulator.onset_eligibility).
          build_args_i$onset_eligibility <- onset_eligibility
          # Cessation-month exposure (simulator.cessation_exposure).
          build_args_i$cessation_exposure <- cessation_exposure
          obj_i <- do.call(build_spatpomp_tu_c, build_args_i)
          bl_i <- make_block_list(obj_i, pc$frame_tu)
          res_i <- run_one_bpfilter_replicate(
            obj_i,
            Np = Np,
            block_list = bl_i,
            seed = seed_i,
            traj_storage = fc_08$traj_storage,
            prune_every = fc_08$traj_prune_every
          )
          list(
            tu = tu_i,
            rep_id = rep_i,
            draw_d = draw_i,
            filter_traj = res_i$filter_traj,
            logLik = res_i$logLik,
            # pomp unit order of this filter (C-locale gid order), so the
            # handoff below maps states by NAME (filter_traj_to_warmup's
            # smoothed_gids), never by frame_tu row position. Same field the
            # smoothing driver records (smoothing_filter_one()).
            gids = sort(as.character(pc$frame_tu$gid), method = "radix")
          )
        },
        error = function(e) {
          warning(sprintf(
            "[10] TU '%s' rep %d: bpfilter failed: %s",
            tu_i,
            rep_i,
            conditionMessage(e)
          ))
          list(
            tu = tu_i,
            rep_id = rep_i,
            draw_d = draw_i,
            filter_traj = NULL,
            logLik = NA_real_
          )
        }
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )

  # Group raw results by TU, save caches
  for (tu in uncached_tus) {
    tu_raw <- Filter(function(r) !is.null(r) && r$tu == tu, raw_results)
    if (length(tu_raw) == 0L) {
      warning(sprintf(
        "[10] TU '%s': no successful filter runs — skipping cache.",
        tu
      ))
      next
    }
    n_ok <- sum(vapply(tu_raw, function(r) !is.null(r$filter_traj), logical(1)))
    message(sprintf(
      "[10] TU '%s': %d/%d filter runs succeeded.",
      tu,
      n_ok,
      n_replicates
    ))

    cp <- .make_filter_cache_path(
      tu,
      cache_dir,
      occ_stem,
      Np,
      n_replicates,
      warmup_window,
      incidence_cap,
      onset_eligibility,
      cessation_exposure,
      smoothing_draw_design_of(flt)
    )
    qs2::qs_save(tu_raw, cp)
    message(sprintf("[10] TU '%s': filter cache saved → %s", tu, basename(cp)))
    filter_results_by_tu[[tu]] <- tu_raw
  }

  filter_results_by_tu
}

# 10. Phase B: resample + forecast ---------------------------------------------
# Sequential per TU.

run_forecast_phase <- function(
  tus,
  filter_results_by_tu,
  art,
  warmup_window,
  forecast_window,
  n_replicates,
  tu_precomp,
  parts_dir,
  diag_dir,
  out_fmt,
  seed_root,
  opt,
  handoff_by_tu = NULL,
  tus_all = tus,
  replicate_offset = 0L,
  chunked = FALSE
) {
  # Diagnostic (hindcast-only): observation-anchored handoff. When enabled,
  # each TU's handoff active/mso is overwritten with the OBSERVED active
  # footprint at the launch month, bypassing PF cold-start under-recruitment.
  # Load the observed duration panel once here (not per TU).
  anchor_obs_active <- isTRUE(opt$simulator$forecast$observed_active_handoff)
  obs_dur_panel <- NULL
  if (anchor_obs_active) {
    dur_path <- make_duration_dat_file(opt)
    if (!file.exists(dur_path)) {
      stop(
        "simulator.forecast.observed_active_handoff = TRUE requires the ",
        "observed duration panel (",
        dur_path,
        ") — run ",
        "analysis/03_prepare_dur.R first."
      )
    }
    obs_dur_panel <- readRDS(dur_path)
    miss <- setdiff(c("gid", "month", "mso"), names(obs_dur_panel))
    if (length(miss)) {
      stop(
        "observed_active_handoff: duration panel ",
        dur_path,
        " is missing column(s): ",
        paste(miss, collapse = ", ")
      )
    }
    obs_dur_panel <- obs_dur_panel[, c("gid", "month", "mso"), drop = FALSE]
    message(
      "[10] observed_active_handoff ON (N4 diagnostic): anchoring handoff ",
      "active/mso to the observed footprint from ",
      dur_path
    )
  }

  for (tu in tus) {
    tu_results <- filter_results_by_tu[[tu]]
    pc <- tu_precomp[[tu]]

    if (is.null(tu_results) || length(tu_results) == 0L) {
      warning(sprintf(
        "[10] TU '%s': no filter results — skipping forecast.",
        tu
      ))
      next
    }
    if (is.null(pc)) {
      warning(sprintf(
        "[10] TU '%s': no precomputed inputs — skipping forecast.",
        tu
      ))
      next
    }

    # Handoff date: in smoothing-input mode it is auto-derived per TU from the
    # loaded trajectory's filt_months (handoff_by_tu); otherwise warmup_window$end.
    # The forecast horizon starts at the handoff so the panel covers the forecast.
    handoff_date <- if (!is.null(handoff_by_tu)) {
      handoff_by_tu[[tu]]
    } else {
      as.Date(warmup_window$end)
    }
    forecast_window_tu <- if (!is.null(handoff_by_tu)) {
      list(start = format(handoff_date), end = forecast_window$end)
    } else {
      forecast_window
    }

    # spatPomp sorts unit_names alphabetically, so filter_traj row u_days{i} ↔
    # sort(gid)[i]. filter_traj_to_warmup maps row i → frame_tu$gid[i] unsorted,
    # so frame_tu must be gid-sorted before the handoff.
    # method = "radix" pins this to byte-order (C-locale-equivalent) collation,
    # matching the dplyr::arrange(..., .locale = "C") that built the pomp's
    # unit axis in spatpomp_core.R — a plain order() under a system ICU locale
    # (e.g. en_US.UTF-8) disagrees with it (e.g. "BDI_1_1" vs "BDI_11_1") and
    # silently swaps two units' state trajectories under their labels.
    frame_tu <- pc$frame_tu
    frame_tu <- frame_tu[
      order(as.character(frame_tu$gid), method = "radix"),
      ,
      drop = FALSE
    ]

    # Extract logLiks
    logliks <- vapply(
      tu_results,
      function(r) r$logLik %||% NA_real_,
      numeric(1)
    )
    n_ok <- sum(!is.na(logliks))
    message(sprintf(
      "[10] TU '%s': logLik range [%.2f, %.2f] (%d/%d runs ok).",
      tu,
      min(logliks, na.rm = TRUE),
      max(logliks, na.rm = TRUE),
      n_ok,
      length(logliks)
    ))

    # Save per-TU filter diagnostics
    diag_path <- file.path(diag_dir, paste0("bpfilter_diag_tu=", tu, ".rds"))
    saveRDS(
      list(
        tu = tu,
        logLiks = logliks,
        Np = Np,
        n_replicates = n_replicates,
        draw_assignments = vapply(tu_results, function(r) r$draw_d, integer(1)),
        warmup_window = warmup_window,
        generated_at = Sys.time()
      ),
      diag_path
    )

    # Build warmup mini-panels from final filter.traj states. The handoff_date is
    # warmup_window$end in native mode, or the auto-derived per-TU date in
    # smoothing-input mode. Guard the unit count against the gid-sorted frame.
    warmup_raw <- lapply(tu_results, function(r) {
      if (is.null(r$filter_traj)) {
        return(NULL)
      }
      n_smooth <- .sp_n_units(dimnames(r$filter_traj)[[1]])
      if (n_smooth < nrow(frame_tu)) {
        stop(sprintf(
          "TU '%s': smoothed trajectory has %d units but forecast frame needs %d — ",
          tu,
          as.integer(n_smooth),
          nrow(frame_tu),
          "forecast frame cannot be larger than the smoothed frame."
        ))
      }
      # When the smoothed frame is a superset of the forecast frame (e.g. the
      # smoothing did not apply unseen_adm1 restrict), pass the smoothed unit
      # ordering so filter_traj_to_warmup can extract by gid name rather than
      # by position (positional extraction would mis-map when a middle unit is
      # absent from the forecast frame).
      smoothed_gids <- r$gids %||% NULL
      tryCatch(
        filter_traj_to_warmup(
          r$filter_traj,
          frame_tu,
          format(handoff_date),
          smoothed_gids = smoothed_gids,
          filt_months = r$filt_months
        ),
        error = function(e) {
          warning(sprintf(
            "[10] TU '%s': filter_traj_to_warmup failed: %s",
            tu,
            conditionMessage(e)
          ))
          NULL
        }
      )
    })

    # Drop failed trajectories
    ok_idx <- which(!vapply(warmup_raw, is.null, logical(1)))
    warmup_raw <- warmup_raw[ok_idx]
    logliks_ok <- logliks[ok_idx]

    if (length(warmup_raw) == 0L) {
      warning(sprintf(
        "[10] TU '%s': all warmup extractions failed — skipping.",
        tu
      ))
      next
    }

    # Handoff states (simulator.forecast.handoff_weighting): "uniform" (default)
    # = balanced assignment over the smoothing replicates by global replicate
    # id; "shuffled" = balanced per block of m ids, each block a seeded
    # permutation (seed_root, block), so a cycle draw count need not be coprime
    # to m; "loglik" = SIR resampling by exp(logLik) (ESS ~1 in practice).
    # spatpomp_core.R::select_handoff_warmup().
    handoff_weighting <- opt$simulator$forecast$handoff_weighting %||% "uniform"
    if (!handoff_weighting %in% c("uniform", "shuffled", "loglik")) {
      stop(
        "simulator.forecast.handoff_weighting must be 'uniform', 'shuffled' or 'loglik', not '",
        handoff_weighting,
        "'",
        call. = FALSE
      )
    }
    # param_draw_mode joint: the joint pairing by smoothing rep_id
    # replaces the uniform one (spatpomp_core.R::.joint_handoff_index()).
    if (joint_draws) {
      handoff_weighting <- "joint"
    }
    warmup_list <- select_handoff_warmup(
      warmup_raw,
      logliks_ok,
      n = n_replicates,
      offset = replicate_offset,
      weighting = handoff_weighting,
      seed = seed_root,
      rep_ids = vapply(
        tu_results[ok_idx],
        function(r) {
          as.integer(r$rep_id %||% NA_integer_)
        },
        integer(1)
      ),
      n_smoothing = joint_n_smoothing
    )
    # Each handoff state's global smoothing draw (the record's draw_d); the
    # forecast draw under joint, recorded in handoff_map.csv.
    warmup_draw_d <- if (joint_draws) {
      handoff_draw_d(tu_results[ok_idx], warmup_list)
    } else {
      NULL
    }
    message(sprintf(
      "[10] TU '%s': %d → %d warmup states, handoff weighting '%s' (logLik ESS %.1f of %d).",
      tu,
      length(warmup_raw),
      length(warmup_list),
      handoff_weighting,
      handoff_loglik_ess(logliks_ok),
      length(warmup_raw)
    ))
    # Forecast replicate -> smoothing replicate (the smoothing run's global
    # rep_id, which its part.qs event tables carry). Written to the partition's
    # handoff_map.csv below: the results readers stitch the factual series from
    # the smoothing output before its launch by this pairing, which they cannot
    # reconstruct (failed extractions were dropped above).
    handoff_map <- handoff_map_table(
      tu,
      tu_results[ok_idx],
      warmup_list,
      replicate_offset,
      handoff_date,
      handoff_weighting,
      draws = warmup_draw_d
    )

    # Observed-active handoff: override the (under-recruited) smoothed handoff active/mso with the
    # observed active footprint at the launch month (= handoff_date - 1 month,
    # each mini-panel's t_start). Post-resample, so it applies to exactly the
    # replicates that will be forecast; composes with spread_param_draws below.
    if (anchor_obs_active) {
      launch_month <- as.Date(handoff_date) %m-% months(1L)
      warmup_list <- anchor_warmup_active_to_observed(
        warmup_list,
        obs_dur_panel,
        launch_month
      )
      message(sprintf(
        "[10] TU '%s': observation-anchored handoff — %d active unit(s) at %s.",
        tu,
        sum(warmup_list[[1]]$active),
        format(launch_month)
      ))
    }

    # Forecast. In chunked mode (--replicate_offset passed) this invocation only
    # simulates one replicate-chunk of `tu`; nest under a repoff= subdirectory so
    # concurrently-running chunks of the same tu never overwrite each other's
    # event_table/monthly_summary — aggregate_partitions() (07b) already globs
    # `_SUCCESS` sentinels recursively, so this needs no change on the read side.
    # Unchunked calls keep the exact prior flat "tu=<tu>" path.
    out_part <- file.path(parts_dir, paste0("tu=", tu))
    if (isTRUE(chunked)) {
      out_part <- file.path(
        out_part,
        sprintf("repoff=%06d", as.integer(replicate_offset))
      )
    }
    dir.create(out_part, recursive = TRUE, showWarnings = FALSE)
    # Before the forecast writes its _SUCCESS sentinel, so a successful
    # partition always carries its map (utils/results_core.R::read_handoff_map()).
    utils::write.csv(
      handoff_map,
      file.path(out_part, "handoff_map.csv"),
      row.names = FALSE
    )

    # Seed: the kernel uses set.seed(seed + rep_id), so every TU/chunk needs a
    # disjoint seed range. Unchunked: preserve the exact prior formula (strided
    # by n_replicates over `tus`, shifted past the filter phase's seed range).
    # Chunked: `tus` may contain a single TU (--only_country), so match(tu, tus)
    # would collapse to a constant across different countries' array tasks —
    # anchor on `tus_all` (the whole group's TU list) instead, with a fixed
    # wide stride per TU (independent of chunk n_replicates) plus the chunk's
    # own replicate_offset, so no two (country, chunk) pairs can collide.
    seed_val <- if (isTRUE(chunked)) {
      as.integer(seed_root) +
        (match(tu, tus_all) - 1L) * 1000000L +
        as.integer(replicate_offset)
    } else {
      # tus_all, not tus: --only_country / --only_tu leave one TU in `tus`, and
      # match(tu, tus) == 1 would give every country the same seeds. On a full
      # run tus == tus_all, so the seeds are unchanged.
      as.integer(seed_root) +
        length(tus_all) * as.integer(n_replicates) +
        (match(tu, tus_all) - 1L) * as.integer(n_replicates)
    }

    tryCatch(
      forecast_from_filtered(
        handoff_list = warmup_list,
        tu = tu,
        frame_tu = frame_tu,
        panel_tu_full = pc$panel_tu_full,
        draws_occ = art$draws_occ,
        # forecast_from_filtered() always runs the plain-R forecast
        # phase (simulate_one_tu()), so draws_dur is passed unconditionally here.
        draws_dur = art$draws_dur,
        # mon (monthly cases): forecast_from_filtered() derives posterior_idx_mon
        # from posterior_idx_occ internally (same wrap-from-occurrence-draw
        # pattern as duration) once draws_mon is non-NULL.
        draws_mon = art$draws_mon,
        enso_scenarios = art$enso_scenarios,
        horizon = forecast_window_tu,
        output_part = out_part,
        out_fmt = out_fmt,
        dry_run = FALSE,
        seed = seed_val,
        replicate_offset = replicate_offset,
        cb_meta_occ = art$cb_meta_occ,
        cb_meta_dur = art$cb_meta_dur,
        cb_meta_mon = art$cb_meta_mon,
        grouping = grp_arg_08,
        # partition_by=country collapses tu to a single country, so
        # update_lags() inside simulate_one_tu always yields o_tu=0 (no other
        # country in the partition to spill from). Force it to 1 instead —
        # the forecast-period counterpart of the same fix applied to the
        # smoothing bpfilter's rprocess (see force_o_tu in smoothing_core.R).
        force_o_tu = identical(partition_by, "country"),
        # Suppresses o_ct/o_tu state tracking AND LP contribution in the
        # forecast-phase kernel too, independent of force_o_tu (which takes
        # lower precedence -- see simulate_one_tu()'s gate).
        drop_country_tu_spillover = isTRUE(
          opt$simulator$drop_country_tu_spillover
        ),
        # Neighbour-contagion in the forecast-phase kernel (simulate_one_tu via
        # forecast_from_filtered). Same guard as the filter build above; gated
        # ONLY by use_nbr, independent of force_o_tu / drop_country_tu_spillover.
        use_nbr = use_nbr_08,
        nbr_rho_km = nbr_cfg_08$rho_km %||% 500L,
        nbr_gadm_path = nbr_cfg_08$gadm_path %||%
          "data/geodata/gadm_valid.gpkg",
        nbr_cache_dir = nbr_cfg_08$cache_dir %||% "generated_data/",
        nbr_spec = nbr_spec_08,
        # Detected-activity covariate feed in the forecast kernel too, with the
        # handoff's own det flag (filter_traj_to_warmup()).
        covariate_feed = fc_08$covariate_feed,
        delta_fn = delta_fn,
        delta_fp = delta_fp,
        # Hold ENSO at each fit's crossbasis boundary knots before the basis is
        # evaluated (simulator.enso.clamp_to_knots; off = earlier behaviour).
        clamp_enso_to_knots = isTRUE(opt$simulator$enso$clamp_to_knots),
        clamp_enso_margin = as.numeric(
          opt$simulator$enso$clamp_margin %||% 0.25
        ),
        # Monthly incidence cap on the forecast mean and drawn cases
        # (simulator.incidence_cap).
        incidence_cap = incidence_cap,
        # Onset gate in the forecast kernel (simulator.onset_eligibility); the
        # smoothing cache's own value is checked against it (expect_meta).
        onset_eligibility = onset_eligibility,
        # Cessation-month exposure in the forecast kernel
        # (simulator.cessation_exposure); checked against the cache likewise.
        cessation_exposure = cessation_exposure,
        # Per-replicate posterior-draw spreading: when the smoothing cache was
        # built with n_draws_in_filter=1 every handoff particle carries
        # draw_idx=1, collapsing the ensemble to one posterior draw and
        # flattening the mon tail. TRUE reassigns each replicate a distinct
        # draw cycling 1..S (see .forecast_param_indices()).
        spread_param_draws = isTRUE(opt$simulator$forecast$spread_param_draws),
        # Draw-selection mode (default "random"): fresh chunk-safe random
        # posterior draw per replicate, shared across a region's member
        # countries. Falls back to spread_param_draws when NULL/blank. See
        # .forecast_param_indices().
        param_draw_mode = opt$simulator$forecast$param_draw_mode,
        # joint: the handoff states' global smoothing draws.
        handoff_draw_d = warmup_draw_d
      ),
      error = function(e) {
        warning(sprintf(
          "[10] TU '%s': forecast failed: %s",
          tu,
          conditionMessage(e)
        ))
      }
    )

    message(sprintf("[07_simulate_spatpomp] TU '%s': done.", tu))
  }
}

# 11. Dispatch -----------------------------------------------------------------

if (isTRUE(opt$dry_run)) {
  message("[07_simulate_spatpomp] [DRY RUN] plan:")
  message(sprintf(
    "  incidence_cap = %s  (monthly cases per unit <= cap x pop)",
    format(incidence_cap)
  ))
  for (tu in tus) {
    if (
      use_smoothing &&
        partition_by %in% c("enso_region", "registry") &&
        length(.smoothing_cache_hits(tu, smoothing_input))
    ) {
      cp <- .resolve_smoothing_cache(tu, smoothing_input)
      message(sprintf(
        "  %s TU '%s': Phase A SKIPPED (native per-TU smoothing cache)  forecast end=%s",
        partition_by,
        tu,
        forecast_window$end
      ))
      message(sprintf(
        "  smoothing cache: %s  [handoff=%s]",
        basename(cp),
        .cache_handoff_label(cp, flt$handoff_month)
      ))
    } else if (
      use_smoothing && partition_by %in% c("enso_region", "registry")
    ) {
      # Region dry-run: report the per-country caches that will be merged into
      # this TU's handoff IC (members read from the partition registry -- the
      # fit grouping for enso_region -- since artefacts are not loaded in
      # dry-run), and the handoff from the first member cache's filt_months.
      members <- if (identical(partition_by, "registry")) {
        sort(names(.pmap_07b)[.pmap_07b == tu])
      } else {
        reg <- tryCatch(
          load_grouping(opt$simulator$grouping$source),
          error = function(e) NULL
        )
        if (!is.null(reg)) {
          sort(unique(reg$iso3[reg$group == tu]))
        } else {
          character(0)
        }
      }
      # Registry units that are unions of whole fit groups can be assembled from
      # the groups' caches when per-country caches are absent (as the real run,
      # load_smoothing_distribution()).
      has_c <- function(m) {
        length(.smoothing_cache_hits(m, smoothing_input)) > 0L
      }
      if (
        identical(partition_by, "registry") &&
          !all(vapply(members, has_c, logical(1)))
      ) {
        greg <- tryCatch(
          load_grouping(opt$simulator$grouping$source),
          error = function(e) NULL
        )
        if (!is.null(greg)) {
          gs <- sort(unique(greg$group[greg$iso3 %in% members]))
          scope_iso <- opt$simulator$scope$countries %||% "all"
          in_scope <- if (identical(scope_iso, "all")) {
            greg$iso3
          } else {
            trimws(strsplit(paste(scope_iso, collapse = ","), ",")[[1L]])
          }
          whole <- all(vapply(
            gs,
            function(g) {
              all(intersect(greg$iso3[greg$group == g], in_scope) %in% members)
            },
            logical(1)
          ))
          if (whole && length(gs) && all(vapply(gs, has_c, logical(1)))) {
            members <- gs
          }
        }
      }
      member_cps <- stats::setNames(
        lapply(members, function(m) {
          .resolve_smoothing_cache(m, smoothing_input)
        }),
        members
      )
      member_cps <- member_cps[!vapply(member_cps, is.null, logical(1))]
      message(sprintf(
        "  region TU '%s': Phase A SKIPPED (assemble from %d member cache(s): %s)  forecast end=%s",
        tu,
        length(member_cps),
        paste(names(member_cps), collapse = ","),
        forecast_window$end
      ))
      if (length(member_cps)) {
        cp0 <- member_cps[[1L]]
        message(sprintf(
          "  member cache (repr.): %s  [handoff=%s]",
          basename(cp0),
          .cache_handoff_label(cp0, flt$handoff_month)
        ))
      }
    } else if (use_smoothing) {
      cp <- .resolve_smoothing_cache(tu, smoothing_input)
      message(sprintf(
        "  TU '%s': Phase A SKIPPED (smoothing input)  forecast end=%s",
        tu,
        forecast_window$end
      ))
      if (!is.null(cp)) {
        message(sprintf(
          "  smoothing cache: %s  [handoff=%s]",
          basename(cp),
          .cache_handoff_label(cp, flt$handoff_month)
        ))
      }
    } else {
      message(sprintf(
        "  TU '%s': Np=%d  n_replicates=%d  warmup=%s – %s  forecast=%s – %s",
        tu,
        Np,
        n_replicates,
        warmup_window$start,
        warmup_window$end,
        forecast_window$start,
        forecast_window$end
      ))
      cp <- .make_filter_cache_path(
        tu,
        cache_dir,
        occ_stem,
        Np,
        n_replicates,
        warmup_window,
        incidence_cap,
        onset_eligibility,
        cessation_exposure,
        smoothing_draw_design_of(flt)
      )
      message(sprintf(
        "  cache: %s  [%s]",
        basename(cp),
        if (file.exists(cp) && !rerun_filter) "will LOAD" else "will RUN"
      ))
    }
  }
  message("[07_simulate_spatpomp] [DRY] done.")
} else if (length(tus) == 0L) {
  message("[07_simulate_spatpomp] no TUs to process — exiting.")
} else {
  # --- set up parallel backend ---
  backend <- opt$simulator$parallel$backend %||% "furrr"
  suppressPackageStartupMessages({
    library(furrr)
    library(future)
  })
  workers <- as.integer(opt$simulator$parallel$workers %||% 1L)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  message(sprintf("[07_simulate_spatpomp] furrr: %d worker(s).", workers))

  # --- Phase A source: native bpfilter OR loaded smoothing distribution ---
  # In smoothing mode the per-TU caches are loaded first so the auto-derived
  # handoff date can drive the per-TU forecast panel built in precompute below.
  handoff_by_tu <- NULL
  if (use_smoothing) {
    message(
      "[07_simulate_spatpomp] smoothing-input mode: skipping Phase A filter."
    )
    if (!is.null(art) && isTRUE(art$grouped)) {
      warning(
        "[07_simulate_spatpomp] grouped mode: smoothing_input must come from ",
        "a smoothing run with a MATCHING grouping config (source/fits_dir); ",
        "draws are routed by this run's grouping at forecast time."
      )
    }
    # Region mode: map each region TU → its member ISO3s from the (already
    # group-attached) frame, so load_smoothing_distribution() knows which
    # per-country caches to merge into each region's handoff IC.
    region_members <- NULL
    if (partition_by %in% c("enso_region", "registry") && !is.null(art)) {
      key_col <- if (identical(partition_by, "registry")) "part" else "group"
      region_members <- lapply(
        split(
          as.character(art$frame$country),
          as.character(art$frame[[key_col]])
        ),
        function(x) sort(unique(x))
      )
    }
    # Registry units that are unions of whole fit groups can be assembled from
    # those groups' (enso_region) smoothing caches; tu_gids lets the loader
    # check that the records tile each TU exactly.
    region_pieces <- NULL
    tu_gids <- NULL
    if (!is.null(art)) {
      tu_gids <- lapply(
        split(as.character(art$frame$gid), as.character(art$frame$tu)),
        unique
      )
      if (identical(partition_by, "registry") && !is.null(art$frame$group)) {
        ctry_of_grp <- split(
          as.character(art$frame$country),
          as.character(art$frame$group)
        )
        region_pieces <- lapply(
          stats::setNames(nm = names(region_members)),
          function(u) {
            gs <- sort(unique(as.character(art$frame$group[
              art$frame$part == u
            ])))
            whole <- all(vapply(
              gs,
              function(g) all(ctry_of_grp[[g]] %in% region_members[[u]]),
              logical(1)
            ))
            list(groups = if (whole) gs else character(0))
          }
        )
      }
    }
    sm <- load_smoothing_distribution(
      tus,
      smoothing_input,
      forecast_window,
      handoff_month = flt$handoff_month,
      partition_by = partition_by,
      region_members = region_members,
      region_pieces = region_pieces,
      tu_gids = tu_gids,
      expect_meta = c(
        smoothing_cache_meta(
          partition_by,
          partition_source,
          grouping_source = if (isTRUE(opt$simulator$grouping$enable)) {
            opt$simulator$grouping$source
          } else {
            ""
          },
          fits_dir = opt$simulator$grouping$fits_dir %||% ""
        ),
        # optional identity: the cache's population
        # year must be this run's
        list(
          pop_year = sim_pop_year_spec(opt$simulator),
          # same onset gate as this run
          onset_eligibility = onset_eligibility,
          # same cessation exposure
          cessation_exposure = cessation_exposure
        ),
        # joint draws: a shared_cycle cache
        # thinned to this run's max_draws
        if (joint_draws) {
          list(
            draw_design = "shared_cycle",
            max_draws = as.integer(opt$simulator$max_draws %||% 200L)
          )
        }
      )
    )
    filter_results_by_tu <- sm$filter_results_by_tu
    handoff_by_tu <- sm$handoff_by_tu

    # The ENSO scenario matrices are positional: their month 1 is the first
    # panel month >= min(horizon.start, cache-filename end token) (see the
    # pre-scan above and load_artefacts()), and the kernel reads row 1 as the
    # handoff month. If horizon.start is left nominal while the smoothing cache
    # ends later than it (e.g. a smoothing extended to the ENSO seam), ENSO
    # would be shifted by months with no error. Refuse instead.
    .es <- as.Date(opt$simulator$horizon$start)
    if (!is.null(.enso_start_override)) {
      .es <- min(.es, .enso_start_override)
    }
    .es_m <- as.Date(format(.es, "%Y-%m-01"))
    if (.es != .es_m) {
      .es_m <- seq(.es_m, by = "month", length.out = 2L)[2L]
    }
    .ho <- Filter(Negate(is.null), handoff_by_tu)
    .bad <- names(.ho)[vapply(.ho, function(h) as.Date(h) != .es_m, logical(1))]
    if (length(.bad)) {
      stop(
        sprintf(
          paste0(
            "[07_simulate_spatpomp] ENSO scenario starts %s but the smoothing handoff is %s ",
            "(TU %s): set simulator.horizon.start and filter.forecast_window.start to the ",
            "handoff month, or ENSO is misaligned by months."
          ),
          format(.es_m),
          format(as.Date(.ho[[.bad[1]]])),
          .bad[1]
        ),
        call. = FALSE
      )
    }
  }

  # Per-TU effective forecast window: starts at the handoff date in smoothing
  # mode (so panel_tu_full covers [handoff, forecast end]); otherwise constant.
  .fw_for_tu <- function(tu) {
    if (use_smoothing && !is.null(handoff_by_tu[[tu]])) {
      list(start = format(handoff_by_tu[[tu]]), end = forecast_window$end)
    } else {
      forecast_window
    }
  }

  # --- pre-compute per-TU inputs (sequential; lightweight) ---
  message("[07_simulate_spatpomp] pre-computing per-TU panels...")
  tu_precomp <- setNames(
    lapply(tus, function(tu) {
      precompute_tu_inputs(tu, art, warmup_window, .fw_for_tu(tu), opt)
    }),
    tus
  )

  # --- Phase A: filter (native bpfilter; skipped in smoothing mode) ---
  if (!use_smoothing) {
    filter_results_by_tu <- run_filter_phase(
      tus = tus,
      art = art,
      Np = Np,
      n_replicates = n_replicates,
      n_draws_in_filter = n_draws_in_filter,
      warmup_window = warmup_window,
      tu_precomp = tu_precomp,
      cache_dir = cache_dir,
      occ_stem = occ_stem,
      rerun_filter = rerun_filter,
      seed_root = seed_root,
      dry_run = FALSE
    )
  }

  # --- Phase B: resample + forecast ---
  run_forecast_phase(
    tus = tus,
    filter_results_by_tu = filter_results_by_tu,
    art = art,
    warmup_window = warmup_window,
    forecast_window = forecast_window,
    n_replicates = n_replicates,
    tu_precomp = tu_precomp,
    parts_dir = parts_dir,
    diag_dir = diag_dir,
    out_fmt = opt$simulator$output$format %||% "qs",
    seed_root = seed_root,
    opt = opt,
    handoff_by_tu = handoff_by_tu,
    tus_all = tus_all,
    replicate_offset = opt$replicate_offset %||% 0L,
    chunked = isTRUE(opt$.replicate_chunked)
  )

  message("[07_simulate_spatpomp] all TUs complete.")
}

message("[07_simulate_spatpomp] done.")
