# paths.R — the one authority for artefact and data paths: every file name
#   the pipeline reads or writes is built here, never by hand
#
# Provides:
#   make_run_file()          universal router to the per-track run-file builders
#   make_<track>_run_file()  occ / duration / size / monthly / sim / smooth /
#                            spatpomp / pop / clim / results / deploy builders
#   make_<artefact>_file()   prep-time data files (surv, clim, panels, WHO, ...)
#   make_loo_compare_run_file(), make_kl_sweep_run_file()   comparison outputs
#   occ_panel_spec_of(), episode_key_of(), national_zero_fill_of()   name keys
#   nbr_covariate_name(), parse_nbr_covariate()   contagion column names
#   read_deployed_lag()      the lag baked into a deployed best-fit file name
#   smoothing_cache_coverage()   coverage of a smoothing cache from its names
#   assert_basename_length() the 255-byte file-name guard
#
# Details:
#   SECURE_OUTBREAK_ROOT (env CHOLERA_OUTBREAK_DATA_DIR overrides it) anchors
#   every outbreak-identifying path; OUTBREAK_DATA_VINTAGE is the one default
#   outbreak-database vintage of the data leg. Bash never builds an artefact
#   name: it asks preprocess_paths.R or deploy_prod_run.R --stage resolve.

# Secure outbreak-data root ------------------------------------------------

#' Root directory for outbreak raw/derived data moved out of the repo.
#'
#' Outbreak-identifying data (location-period IDs, dates, case counts) lives
#' outside the repo tree in a restricted folder that pipeline code may
#' reference by path but that is otherwise off-limits to automated tooling.
#' All secure-outbreak-data path builders below are anchored on this root.
#'
#' The environment variable CHOLERA_OUTBREAK_DATA_DIR, when set, replaces the
#' default root `~/secure/cholera_outbreak_data` (another machine, or the
#' synthetic world of 00_ingest_synthetic.R). It is read when this module is
#' sourced. config.R restates the same rule for the two defaults it holds
#' (population.lp_shapefiles, simulator.location_frame.outbreak_metadata_path),
#' since modules never source each other.
SECURE_OUTBREAK_ROOT <- local({
  d <- Sys.getenv("CHOLERA_OUTBREAK_DATA_DIR")
  if (nzchar(d)) sub("/+$", "", d) else "~/secure/cholera_outbreak_data"
})
SECURE_OUTBREAK_DIR <- paste0(path.expand(SECURE_OUTBREAK_ROOT), "/")

#' Outbreak-database vintage the data leg reads by default.
#'
#' One value for every script that reads an OutbreakExtractR / legacy-RDS
#' snapshot (00 ingest, 01 outbreaks, 01b monthly, 01 national-zero months,
#' 01 data-overview figure, 02 climate, 02 crosswalk), so a run without
#' `--indir` cannot mix two vintages (LP ids of a newer refresh unmatched in an
#' older crosswalk are dropped silently by the 03 scripts). Change the vintage here,
#' and only here; `--indir` still overrides it per run, and 00/01/02 stamp the
#' directory they read into generated_data/prep_counts/<step>.yml
#' (provenance.R::write_step_counts(indir =)).
OUTBREAK_DATA_VINTAGE <- "outbreak_data_Nov_2026"

#' Default `--indir` of the data-leg scripts: the outbreak database snapshot
#' `OUTBREAK_DATA_VINTAGE` under the secure outbreak-data root.
#'
#' Returned unexpanded (`~/...`), as the per-script literals it replaces were,
#' so `--help` shows a portable path.
#' @return Character scalar, "<SECURE_OUTBREAK_ROOT>/<vintage>"
#'   ("~/secure/cholera_outbreak_data/<vintage>" unless CHOLERA_OUTBREAK_DATA_DIR is set).
default_outbreak_data_dir <- function() {
  file.path(SECURE_OUTBREAK_ROOT, OUTBREAK_DATA_VINTAGE)
}

# Filenames ---------------------------------------------------------------

#' Path to a per-(variable, gid) raw climate-stats file.
#'
#' @param var Climate variable name (e.g. "precip", "temp_daymean").
#' @param gid GADM admin unit identifier.
#' @return Character path under generated_data/climate_stats/.
make_clim_stats_file <- function(var, gid) {
  assert_basename_length(
    glue::glue("generated_data/climate_stats/clim_stats_{var}_{gid}.rds"),
    "make_clim_stats_file"
  )
}

#' Path to a per-(variable, gid, baseline) baseline climate-stats file.
#'
#' @param var Climate variable name.
#' @param gid GADM admin unit identifier.
#' @param baseline Baseline period tag (e.g. "1990_2019").
#' @return Character path under generated_data/baseline_climate_stats/.
make_clim_baseline_stats_file <- function(var, gid, baseline) {
  assert_basename_length(
    glue::glue(
      "generated_data/baseline_climate_stats/baseline_clim_stats_{var}_{gid}_{baseline}.rds"
    ),
    "make_clim_baseline_stats_file"
  )
}

#' Path to a per-(variable, gid, baseline) climate-extremes file.
#'
#' @param var Climate variable name.
#' @param gid GADM admin unit identifier.
#' @param baseline Baseline period tag.
#' @return Character path under generated_data/climate_stats/.
make_clim_extremes_file <- function(var, gid, baseline) {
  assert_basename_length(
    glue::glue(
      "generated_data/climate_stats/clim_extremes_{var}_{gid}_{baseline}.rds"
    ),
    "make_clim_extremes_file"
  )
}

#' Path to the ADM1 x month rainfall series (03_prepare_rain_adm1.R), read by
#' 09's rain_enso step (Figure 1d). A few MB, so it can live on the cluster
#' where climate_stats.sqlite (8.6 GB) does not.
#'
#' @param outdir Output root.
#' @return Character path; `sub("\\.qs$", "_meta.yml", .)` is its sidecar.
make_adm1_precip_file <- function(outdir = "generated_data/") {
  file.path(outdir, "adm1_monthly_precip.qs")
}

#' Path to the prepared occurrence-track (survival) model dataset.
#'
#' `make_options_from_config()` always fills `opt$baseline` (`"1990_2019"`,
#' the only baseline for which `outbreaks_survival_data_base-<baseline>.rds`
#' exists on disk); call sites that hand-build an option list default it to
#' the same literal. This builder has no default of its own.
#'
#' There is no defensible default here: unlike `size_mean_rate_estimator()`
#' below, `baseline` has no "no-block" case to license a fallback — it is a
#' shared top-level key that `make_options_from_config()` always populates
#' (`utils/config.R`), so a NULL here can only mean the caller built `opt` by
#' hand and forgot it. Hard error instead of a guess.
#'
#' @param opt Resolved options list; only `opt$baseline` is used.
#' @return Character path: generated_data/outbreaks_survival_data_base-<baseline>.rds.
# ---------------------------------------------------------------------------
# Episode-key token (cross-track episode definition)
# ---------------------------------------------------------------------------
# `episodes.key` selects what ONE outbreak episode is, for all four tracks at
# once (see utils/footprint.R::build_canonical_episodes()):
#   lp_season   today's behaviour (DEFAULT): occurrence censors per location
#               period and unions to gid afterwards; size/duration/monthly key
#               on (gid, year_seas, admin_level) with first()-taken dates.
#   gid_season  one canonical episode per (gid, year_seas): onset = earliest
#               member TL, end = latest member TR, one geometry-derived
#               admin_level per gid.
# The token is "" under the default, so EVERY existing artefact name -- prep
# frames, fits, base-rate tables -- is unchanged. Under gid_season it is
# "_ep-gs" everywhere, so the two definitions coexist on disk instead of
# overwriting each other (prep artefacts included: an episode key changes
# the ROWS of a prep frame, not a policy applied to the same rows -- see the
# note at prep_policy_stamp_path() in utils/provenance.R).
EPISODE_KEYS <- c("lp_season", "gid_season")

#' Resolve the episode key from an opt list (config-resolved or CLI).
#' A bare `list(baseline = ...)` resolves to the default, by design.
episode_key_of <- function(opt) {
  ep <- opt$episodes %||% list()
  k <- ep$key %||% "lp_season"
  if (!k %in% EPISODE_KEYS) {
    stop(
      "episodes.key = '",
      k,
      "' is not one of {",
      paste(EPISODE_KEYS, collapse = ", "),
      "}. Refusing to guess.",
      call. = FALSE
    )
  }
  k
}

.episode_key_path_token <- function(opt) {
  if (identical(episode_key_of(opt), "gid_season")) "_ep-gs" else ""
}

#' Resolve the national zero-fill flag (`episodes.national_zero_fill`) from an opt list.
#' Adds no-onset rows to the occurrence panel for the national-zero months of
#' each gid's country (03_prepare_occ.R --national_zero_fill); defined for
#' gid_season only. A bare `list(baseline = ...)` resolves to FALSE.
national_zero_fill_of <- function(opt) {
  nzf <- (opt$episodes %||% list())$national_zero_fill %||% FALSE
  if (!is.logical(nzf) || length(nzf) != 1L || is.na(nzf)) {
    stop("episodes.national_zero_fill must be true or false.", call. = FALSE)
  }
  if (nzf && !identical(episode_key_of(opt), "gid_season")) {
    stop(
      "episodes.national_zero_fill = true requires episodes.key = 'gid_season'.",
      call. = FALSE
    )
  }
  nzf
}

.national_zero_threshold_check <- function(threshold) {
  thr <- suppressWarnings(as.integer(threshold %||% 0L))
  if (
    length(thr) != 1L ||
      is.na(thr) ||
      thr < 0L ||
      !isTRUE(all.equal(as.numeric(threshold %||% 0L), as.numeric(thr)))
  ) {
    stop(
      "episodes.national_zero_threshold must be a whole number >= 0, got: ",
      paste(format(threshold), collapse = ", "),
      call. = FALSE
    )
  }
  thr
}

#' Resolve the national-zero case threshold (`episodes.national_zero_threshold`).
#' 0 (default) = exact WHO zeros only; N > 0 = WHO country-years with fewer
#' than N cases also count (01_preprocess_national_zero_months.R
#' --zero_threshold; an outbreak needs >= 50 cases). Only meaningful with
#' `episodes.national_zero_fill = true`.
# The threshold of a raw `episodes` block (deploy YAML / sim template), NULL-safe.
.episodes_nz_threshold <- function(ep) {
  .national_zero_threshold_check((ep %||% list())$national_zero_threshold)
}

national_zero_threshold_of <- function(opt) {
  thr <- .national_zero_threshold_check(
    (opt$episodes %||% list())$national_zero_threshold
  )
  if (thr > 0L && !national_zero_fill_of(opt)) {
    stop(
      "episodes.national_zero_threshold = ",
      thr,
      " requires episodes.national_zero_fill = true.",
      call. = FALSE
    )
  }
  thr
}

# "_nz" changes the occurrence panel's ROWS, so it is a filename token, placed
# AFTER "_ep-gs" so the "_ep-gs" substring filters keep matching. Only the
# survival panel and what is estimated on it (occurrence fit, occurrence base
# rate) carry it; size / duration / monthly frames and the sim climate panel
# (same gid set and min_seas) do not change. A threshold variant carries the
# threshold: "_nz50".
.national_zero_fill_path_token <- function(opt) {
  .occ_panel_nz_token(occ_panel_spec_of(opt))
}

# Fit-file filter for the "_nz" token, the companion of the "_ep-gs" filters
# in read_deployed_lag() / simulations.R::.discover_group_fit(): an occurrence
# fit on the national zero-fill panel is invisible to a lookup without it and vice
# versa, and "_nz" / "_nz50" never match each other (exact token). Only
# occurrence artefacts carry the token, so callers apply it to occurrence
# lookups only. A wrapper over the exact whole-suffix matcher
# .occ_panel_variant_file_match(), so a "_nz60" lookup never matches
# "..._nz60_ae_a1c" files nor reads "..._nz60-r2" files as "no zero-fill".
.national_zero_fill_file_match <- function(
  files,
  national_zero_fill = FALSE,
  national_zero_threshold = 0L
) {
  thr <- .national_zero_threshold_check(national_zero_threshold)
  .occ_panel_variant_file_match(
    files,
    .occ_panel_spec(
      key = "gid_season",
      national_zero_fill = isTRUE(national_zero_fill),
      national_zero_threshold = if (isTRUE(national_zero_fill)) thr else 0L
    )
  )
}

# --- Occurrence-panel spec ----------------------------------------------------
# ONE object for every episodes.* key that changes the occurrence panel's ROWS
# (and so its filename): the national zero-fill pair and its rule, the at-risk
# start and the ADM1-overlap treatment. It is threaded in place of the paired
# (national_zero_fill, national_zero_threshold) arguments; at the field defaults
# the token is "" / "_ep-gs" / "_ep-gs_nz" / "_ep-gs_nz<thr>".
#   national_zero_rule  ""  (= the `legacy` entry of analysis/national_zero_rules.yml)
#                       | a registry name; token "_nz<thr>-<rule>". Needs the fill.
#   at_risk_start       season (default) | after_episode; token "_ae". gid_season only.
#   adm1_overlap        drop (default) | censor_orphan; token "_a1c". gid_season only.
OCC_PANEL_AT_RISK_STARTS <- c("season", "after_episode")
OCC_PANEL_ADM1_OVERLAPS <- c("drop", "censor_orphan")
.OCC_PANEL_SPEC_DEFAULTS <- list(
  national_zero_fill = FALSE,
  national_zero_threshold = 0L,
  national_zero_rule = "",
  at_risk_start = "season",
  adm1_overlap = "drop"
)

# A rule name: "" (legacy) or a short lowercase registry key; "legacy" itself
# normalises to "" so one rule has one token (none).
.national_zero_rule_check <- function(rule) {
  r <- rule %||% .OCC_PANEL_SPEC_DEFAULTS$national_zero_rule
  if (
    !is.character(r) ||
      length(r) != 1L ||
      is.na(r) ||
      !(identical(r, "") || grepl("^[a-z][a-z0-9]{0,7}$", r))
  ) {
    stop(
      "episodes.national_zero_rule must be \"\" or a lowercase name of at most 8 ",
      "characters (analysis/national_zero_rules.yml), got: ",
      paste(format(rule), collapse = ", "),
      call. = FALSE
    )
  }
  if (identical(r, "legacy")) "" else r
}

.occ_panel_enum_check <- function(x, choices, key) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    stop(
      "episodes.",
      key,
      " = '",
      paste(format(x), collapse = ", "),
      "' is not one of {",
      paste(choices, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
  x
}

# Validated spec from explicit fields. Defaults come from
# .OCC_PANEL_SPEC_DEFAULTS (value-checked against get_default_options()$episodes
# by test_occ_panel_spec.R), never from literals. The key / fill / threshold checks reuse
# episode_key_of() / national_zero_fill_of() / national_zero_threshold_of() so their
# error messages are unchanged.
.occ_panel_spec <- function(
  key = "lp_season",
  national_zero_fill = FALSE,
  national_zero_threshold = 0L,
  national_zero_rule = "",
  at_risk_start = "season",
  adm1_overlap = "drop"
) {
  ep <- list(
    key = key,
    national_zero_fill = national_zero_fill,
    national_zero_threshold = national_zero_threshold
  )
  o <- list(episodes = ep)
  k <- episode_key_of(o)
  f <- national_zero_fill_of(o)
  t <- national_zero_threshold_of(o)
  r <- .national_zero_rule_check(national_zero_rule)
  if (nzchar(r) && !f) {
    stop(
      "episodes.national_zero_rule = '",
      r,
      "' requires episodes.national_zero_fill = true.",
      call. = FALSE
    )
  }
  a <- .occ_panel_enum_check(
    at_risk_start %||% .OCC_PANEL_SPEC_DEFAULTS$at_risk_start,
    OCC_PANEL_AT_RISK_STARTS,
    "at_risk_start"
  )
  d <- .occ_panel_enum_check(
    adm1_overlap %||% .OCC_PANEL_SPEC_DEFAULTS$adm1_overlap,
    OCC_PANEL_ADM1_OVERLAPS,
    "adm1_overlap"
  )
  for (kv in list(
    c("at_risk_start", a, "season"),
    c("adm1_overlap", d, "drop")
  )) {
    if (!identical(kv[[2]], kv[[3]]) && !identical(k, "gid_season")) {
      stop(
        "episodes.",
        kv[[1]],
        " = '",
        kv[[2]],
        "' requires episodes.key = 'gid_season'.",
        call. = FALSE
      )
    }
  }
  structure(
    list(
      key = k,
      national_zero_fill = f,
      national_zero_threshold = t,
      national_zero_rule = r,
      at_risk_start = a,
      adm1_overlap = d
    ),
    class = "occ_panel_spec"
  )
}

#' Occurrence-panel spec of a raw `episodes` block (config, deploy YAML, sim
#' template); NULL-safe, every absent key at its default (read with [[ ]]).
.occ_panel_spec_from_episodes <- function(ep) {
  ep <- ep %||% list()
  d <- .OCC_PANEL_SPEC_DEFAULTS
  .occ_panel_spec(
    key = ep[["key"]] %||% EPISODE_KEYS[[1L]],
    national_zero_fill = ep[["national_zero_fill"]] %||% d$national_zero_fill,
    national_zero_threshold = ep[["national_zero_threshold"]] %||%
      d$national_zero_threshold,
    national_zero_rule = ep[["national_zero_rule"]] %||% d$national_zero_rule,
    at_risk_start = ep[["at_risk_start"]] %||% d$at_risk_start,
    adm1_overlap = ep[["adm1_overlap"]] %||% d$adm1_overlap
  )
}

#' Occurrence-panel spec of an opt list (config-resolved or CLI). A bare
#' `list(baseline = ...)` resolves to every default.
occ_panel_spec_of <- function(opt) {
  .occ_panel_spec_from_episodes((opt %||% list())[["episodes"]])
}

.as_occ_panel_spec <- function(x) {
  if (inherits(x, "occ_panel_spec")) x else occ_panel_spec_of(x)
}

# "_nz[<thr>][-<rule>]" or "".
.occ_panel_nz_token <- function(spec) {
  if (!isTRUE(spec$national_zero_fill)) {
    return("")
  }
  paste0(
    "_nz",
    if (spec$national_zero_threshold > 0L) spec$national_zero_threshold else "",
    if (nzchar(spec$national_zero_rule)) {
      paste0("-", spec$national_zero_rule)
    } else {
      ""
    }
  )
}

#' Panel-variant filename suffix, in the fixed order
#' "_ep-gs" -> "_nz<thr>[-<rule>]" -> "_ae" -> "_a1c". `x` is an opt list or a spec.
.occ_panel_variant_token <- function(x) {
  s <- .as_occ_panel_spec(x)
  paste0(
    if (identical(s$key, "gid_season")) "_ep-gs" else "",
    .occ_panel_nz_token(s),
    if (identical(s$at_risk_start, "after_episode")) "_ae" else "",
    if (identical(s$adm1_overlap, "censor_orphan")) "_a1c" else ""
  )
}

# The variant tokens a filename carries after "_ep-gs" ("" when none, or when
# the name has no "_ep-gs"): the WHOLE run of variant tokens, in grammar order,
# ending at "_" or ".rds".
.OCC_PANEL_VARIANT_RE <- "_ep-gs((?:_nz[0-9]*(?:-[a-z][a-z0-9]*)?)?(?:_ae)?(?:_a1c)?)(?=_|\\.rds$)"
.occ_panel_file_variant <- function(files) {
  b <- basename(files)
  m <- regmatches(b, regexec(.OCC_PANEL_VARIANT_RE, b, perl = TRUE))
  vapply(m, function(x) if (length(x) >= 2L) x[[2L]] else "", character(1))
}

#' Exact whole-suffix fit-file filter for the occurrence-panel variant: TRUE
#' where the file's variant tokens equal the spec's exactly ("_nz60" never
#' matches "_nz60-r2" or "_nz60_ae_a1c" files, nor the reverse). The "_ep-gs"
#' key filter stays with the callers. `spec` is a spec or an opt list.
.occ_panel_variant_file_match <- function(files, spec) {
  s <- .as_occ_panel_spec(spec)
  want <- sub(
    "^_ep-gs",
    "",
    .occ_panel_variant_token(
      .occ_panel_spec(
        key = "gid_season",
        national_zero_fill = s$national_zero_fill,
        national_zero_threshold = s$national_zero_threshold,
        national_zero_rule = s$national_zero_rule,
        at_risk_start = s$at_risk_start,
        adm1_overlap = s$adm1_overlap
      )
    )
  )
  .occ_panel_file_variant(files) == want
}

#' A spec's `episodes` block, written the way every generated config / opt
#' always wrote it: `key`, then each other field only when non-default, except
#' the names in `always` (written at any value, for the sites that always did).
.episodes_block <- function(spec, always = character()) {
  s <- .as_occ_panel_spec(spec)
  aw <- function(f) f %in% always
  out <- list(key = s$key)
  if (aw("national_zero_fill")) {
    out$national_zero_fill <- s$national_zero_fill
  } else if (s$national_zero_fill) {
    out$national_zero_fill <- TRUE
  }
  if (aw("national_zero_threshold")) {
    out$national_zero_threshold <- s$national_zero_threshold
  } else if (s$national_zero_fill && s$national_zero_threshold > 0L) {
    out$national_zero_threshold <- s$national_zero_threshold
  }
  if (aw("national_zero_rule") || nzchar(s$national_zero_rule)) {
    out$national_zero_rule <- s$national_zero_rule
  }
  if (aw("at_risk_start") || !identical(s$at_risk_start, "season")) {
    out$at_risk_start <- s$at_risk_start
  }
  if (aw("adm1_overlap") || !identical(s$adm1_overlap, "drop")) {
    out$adm1_overlap <- s$adm1_overlap
  }
  out
}

#' Sidecar of an occurrence panel: its full LP -> gid member bridge
#' (`<panel>.lp_members.csv`; columns location_period_id, gid, country), written
#' by 03_prepare_occ.R under gid_season and read by coverage.R::coverage_lp_bridge().
occ_panel_lp_members_path <- function(surv_path) {
  sub("\\.rds$", ".lp_members.csv", as.character(surv_path))
}

#' Sidecars of an occurrence panel (03_prepare_occ.R, gid_season):
#' `<panel>.adm1_overlap.csv` (adm1_overlap_parents()), `.discarded_onsets.csv`
#' (discarded_onset_table()), `.adm0_dropped.csv` (adm0_episode_counts()) and
#' `.masked_rows.csv` (mask_occ_panel_rows(); only when a mask is on).
occ_panel_sidecar_paths <- function(surv_path) {
  p <- as.character(surv_path)
  list(
    adm1_overlap = sub("\\.rds$", ".adm1_overlap.csv", p),
    discarded_onsets = sub("\\.rds$", ".discarded_onsets.csv", p),
    adm0_dropped = sub("\\.rds$", ".adm0_dropped.csv", p),
    masked_rows = sub("\\.rds$", ".masked_rows.csv", p)
  )
}

make_surv_dat_file <- function(opt) {
  if (is.null(opt$baseline)) {
    stop(
      "make_surv_dat_file: opt$baseline is NULL. `baseline` is a shared ",
      "top-level key that make_options_from_config() always populates ",
      "(default \"1990_2019\"), so NULL means this `opt` was hand-built ",
      "and the caller forgot to set it -- there is no defensible default ",
      "to guess here. Pass opt$baseline explicitly.",
      call. = FALSE
    )
  }
  assert_basename_length(
    glue::glue(
      "generated_data/outbreaks_survival_data_base-{opt$baseline}",
      "{.occ_panel_variant_token(opt)}.rds"
    ),
    "make_surv_dat_file"
  )
}

#' Path to the prepared duration-track model dataset.
#'
#' @param opt Resolved options list; only `opt$baseline` is used.
#' @return Character path: generated_data/outbreaks_duration_data_base-<baseline>.rds.
#' @details Named on the baseline alone, exactly like `make_surv_dat_file()` above
#'   -- the two are the same kind of artefact (a prepared model frame whose
#'   only naming-relevant parameter is the climate baseline).
make_duration_dat_file <- function(opt) {
  assert_basename_length(
    glue::glue(
      "generated_data/outbreaks_duration_data_base-{opt$baseline}",
      "{.episode_key_path_token(opt)}.rds"
    ),
    "make_duration_dat_file"
  )
}

#' Path to the prepared size-track episode table.
#'
#' @param opt options list; only `episodes$key` is read (via the token), so
#'   `NULL` / a bare list resolves to the historical literal
#'   `generated_data/outbreak_size_analysis_data.rds` exactly.
#' @details Named WITHOUT a baseline, as it always has been (the size frame's
#'   climate is joined at the outbreak start month and its baseline was a
#'   hardcoded "1990_2019"); kept that way so every config that names the
#'   literal path (`size.data_path` / `simulator.size_data_path`) still
#'   resolves to the default artefact.
make_size_dat_file <- function(opt = NULL) {
  assert_basename_length(
    glue::glue(
      "generated_data/outbreak_size_analysis_data",
      "{.episode_key_path_token(opt %||% list())}.rds"
    ),
    "make_size_dat_file"
  )
}

#' Path to the per-(uid, location_period_id, month) outbreak monthly-cases file.
#'
#' @return Character path: <SECURE_OUTBREAK_DIR>outbreak_monthly_cases.rds.
#' @details Written by `01b_preprocess_outbreaks_monthly.R`. Deliberately takes
#'   NO `opt`: unlike the prepared model frames above, this artefact is a pure
#'   re-aggregation of the raw interval series and carries no climate
#'   baseline, no population policy and no composite policy -- nothing that
#'   would legitimately fork the filename. The parameters that DID shape it
#'   (`--indir` above all) are recorded in the `.policy.rds` sidecar instead,
#'   same as the duration track (`03_prepare_dur.R:183`). Lives under
#'   `SECURE_OUTBREAK_DIR` (outbreak-identifying data), not `generated_data/`.
make_monthly_cases_file <- function() {
  glue::glue("{SECURE_OUTBREAK_DIR}outbreak_monthly_cases.rds")
}

#' Path to the outbreak-level metadata table.
#'
#' @return Character path: <SECURE_OUTBREAK_DIR>outbreak_metadata.rds.
#' @details Written by `01_preprocess_outbreaks.R`; carries outbreak identity
#'   (location-period IDs, dates, case counts), so it lives under
#'   `SECURE_OUTBREAK_DIR`, not `generated_data/`.
make_outbreak_metadata_file <- function() {
  glue::glue("{SECURE_OUTBREAK_DIR}outbreak_metadata.rds")
}

#' Path to the aggregated pre-outbreak-extraction time-series file.
#'
#' @return Character path: <SECURE_OUTBREAK_DIR>aggregated_ts_preoutbreak_extraction.rds.
#' @details Written by `01_preprocess_outbreaks.R`; carries outbreak-level
#'   time series, so it lives under `SECURE_OUTBREAK_DIR`, not `generated_data/`.
make_agg_ts_file <- function() {
  glue::glue("{SECURE_OUTBREAK_DIR}aggregated_ts_preoutbreak_extraction.rds")
}

#' Path to the aggregated admin-level pre-outbreak-extraction time-series file.
#'
#' @return Character path: <SECURE_OUTBREAK_DIR>aggregated_admin_level_ts_preoutbreak_extraction.rds.
#' @details Written by `01_preprocess_outbreaks.R`; carries outbreak-level
#'   time series, so it lives under `SECURE_OUTBREAK_DIR`, not `generated_data/`.
make_agg_adm_ts_file <- function() {
  glue::glue(
    "{SECURE_OUTBREAK_DIR}aggregated_admin_level_ts_preoutbreak_extraction.rds"
  )
}

#' Path to the national-zero months (utils/national_zero.R).
#'
#' @return Character path: <SECURE_OUTBREAK_DIR>national_zero_months.rds.
#' @details Written by `01_preprocess_national_zero_months.R`; country-months
#'   derived from the outbreak time series, so it lives under
#'   `SECURE_OUTBREAK_DIR` beside the aggregated series. Read by
#'   `coverage.R::build_coverage_mask()`.
#'   A positive `threshold` (01_preprocess_national_zero_months.R
#'   --zero_threshold) names the threshold variant
#'   `national_zero_months_t<threshold>.rds`; 0 is the exact-zero file.
#'   A named rule (analysis/national_zero_rules.yml; 01 --rule; "" / "legacy" =
#'   no token) appends `_<rule>`: `national_zero_months_t60_r2.rds`,
#'   `national_zero_months_r2.rds`.
make_national_zero_months_file <- function(threshold = 0L, rule = "") {
  threshold <- .national_zero_threshold_check(threshold)
  rule <- .national_zero_rule_check(rule)
  rtok <- if (nzchar(rule)) paste0("_", rule) else ""
  if (threshold == 0L && !nzchar(rule)) {
    glue::glue("{SECURE_OUTBREAK_DIR}national_zero_months.rds")
  } else if (threshold == 0L) {
    glue::glue("{SECURE_OUTBREAK_DIR}national_zero_months{rtok}.rds")
  } else {
    glue::glue(
      "{SECURE_OUTBREAK_DIR}national_zero_months_t{threshold}{rtok}.rds"
    )
  }
}

#' The national-zero months file an occurrence-panel spec (or opt) reads: its
#' threshold and rule (the exact-zero legacy file when the fill is off).
national_zero_months_file_of <- function(x) {
  s <- .as_occ_panel_spec(x)
  make_national_zero_months_file(
    s$national_zero_threshold,
    s$national_zero_rule
  )
}

#' Sidecars of a national-zero months file (01_preprocess_national_zero_months.R),
#' in the secure dir beside the months: the audit (`<months>.audit.csv`,
#' rule-specific) and the persisted db_who_consistency() table. The latter is
#' named on the RULE-FREE months stem of the threshold --
#' `sub("\\.rds$", "_db_who_consistency.csv", make_national_zero_months_file(threshold))`,
#' the exact expression report_who_scale reads -- because the
#' table depends on the database, the WHO mask and db_tolerance only, not on the
#' rule (every registered rule keeps db_tolerance 1.2).
national_zero_audit_path <- function(months_path) {
  sub("\\.rds$", ".audit.csv", as.character(months_path))
}
national_zero_consistency_path <- function(threshold = 0L) {
  sub(
    "\\.rds$",
    "_db_who_consistency.csv",
    as.character(make_national_zero_months_file(threshold))
  )
}

#' Path to the WHO dashboard national monthly case counts.
#'
#' @param outdir output directory (default `generated_data`).
#' @return Character path: <outdir>/who_dashboard_monthly_africa.rds.
#' @details Written by `03_prepare_who_monthly.R` from the public WHO dashboard
#'   weekly series (national totals, no outbreak identity), so it lives under
#'   `generated_data/`. Read by `smoothing_core.R` when
#'   `simulator.filter.who_obs.enable` is set, and by `08c_smoothing_compile.R`.
make_who_monthly_file <- function(outdir = "generated_data") {
  file.path(outdir, "who_dashboard_monthly_africa.rds")
}

#' Paths to the WHO national annual and weekly tables.
#'
#' @param kind "report" (WHO annual report, 1970-2024), "dashboard_weekly"
#'   (WHO dashboard national weekly counts, 2023 onwards) or
#'   "dashboard_annual" (the dashboard's annual totals).
#' @param dir the folder (default `data/who_annual_report_africa`, the default
#'   of 01_preprocess_national_zero_months.R's --who_annual_dir).
#' @return Character path. Read by utils/who_annual.R's builders (01's national
#'   zero mask; 09's report_who_enso_ratio).
make_who_annual_file <- function(
  kind = c("report", "dashboard_weekly", "dashboard_annual"),
  dir = "data/who_annual_report_africa"
) {
  kind <- match.arg(kind)
  file.path(
    dir,
    c(
      report = "who_annual_report_africa.rds",
      dashboard_weekly = "who_dashboard_weekly_africa.rds",
      dashboard_annual = "who_dashboard_annual_africa.rds"
    )[[kind]]
  )
}

#' Path to the per-footprint-member membership file.
#'
#' @return Character path: generated_data/outbreak_footprint_members.rds.
#' @details Records which (uid, location_period_id) rows fed each (gid,
#'   year_seas, admin_level) footprint, written by `03_prepare_siz.R` alongside
#'   the size frame. Takes no `opt` for the same reason the size frame itself
#'   is a hardcoded literal at `03_prepare_siz.R:834` -- the two are produced
#'   by one invocation and must be read as a pair, so a second naming scheme
#'   would only let them drift apart. The pop/composite policy that shaped
#'   both lives in the size frame's own `.policy.rds` sidecar.
make_footprint_members_file <- function(opt = NULL) {
  assert_basename_length(
    glue::glue(
      "generated_data/outbreak_footprint_members",
      "{.episode_key_path_token(opt %||% list())}.rds"
    ),
    "make_footprint_members_file"
  )
}

#' Path to the prepared monthly-cases track model dataset.
#'
#' @param opt Resolved options list; only `opt$baseline` is used.
#' @return Character path: generated_data/outbreaks_monthly_data_base-<baseline>.rds.
#' @details Named on the baseline alone, exactly like `make_surv_dat_file()` /
#'   `make_duration_dat_file()` above -- the same kind of artefact (a prepared
#'   model frame whose only naming-relevant parameter is the climate
#'   baseline), and it is built by expanding the duration panel, which is
#'   named that way too.
make_monthly_dat_file <- function(opt) {
  assert_basename_length(
    glue::glue(
      "generated_data/outbreaks_monthly_data_base-{opt$baseline}",
      "{.episode_key_path_token(opt)}.rds"
    ),
    "make_monthly_dat_file"
  )
}

#' Path to the continuous per-(gid, month) simulator climate panel.
#'
#' @param opt Resolved options list; only `opt$baseline` is used.
#' @return Character path: generated_data/sim_climate_panel_base-<baseline>.rds.
make_sim_panel_file <- function(opt) {
  # Episode-key token: the panel's gid scope and per-gid season calendar come
  # from the survival panel, so a gid_season panel is a different artefact.
  assert_basename_length(
    glue::glue(
      "generated_data/sim_climate_panel_base-{opt$baseline}",
      "{.episode_key_path_token(opt)}.rds"
    ),
    "make_sim_panel_file"
  )
}

#' Path to the cached ENSO analog matchings from Lenssen et al. ensembles.
#'
#' @param cache_dir Directory to write/read the cache file (default "generated_data/").
#' @return Character path: <cache_dir>/lenssen_nino34_analogs.rds.
make_lenssen_analog_file <- function(cache_dir = "generated_data/") {
  file.path(cache_dir, "lenssen_nino34_analogs.rds")
}

#' Directory of one past ENSO event's Lenssen analog cache (event attribution,
#' deploy.event_attribution.lenssen_dir): 03_prepare_lenssen_analogs.R
#' --output_dir <this> writes make_lenssen_analog_file(<this>) + its meta YAML.
make_event_lenssen_dir <- function(lenssen_dir, event_id) {
  assert_basename_length(
    file.path(lenssen_dir, event_id),
    "make_event_lenssen_dir"
  )
}

#' smoothing_cache_coverage
#'
#' Months covered by every TU of a historical smoothing cache, parsed from the
#' cache filenames hist_filter_tu=<tu>_..._<start8>_<end8>.qs: filtered months
#' cover [start, end), so the last available handoff IC month is end - 1 month.
#' The intersection across TUs (max start, min end). NULL when the dir holds no
#' cache file. Shared by make_hindcast_configs.R and
#' make_event_attribution_configs.R.
smoothing_cache_coverage <- function(cache_dir) {
  hits <- Sys.glob(file.path(cache_dir, "hist_filter_tu=*.qs"))
  if (!length(hits)) {
    return(NULL)
  }
  starts <- ends <- as.Date(character(0))
  for (h in hits) {
    toks <- regmatches(basename(h), gregexpr("[0-9]{8}", basename(h)))[[1]]
    if (length(toks) < 2L) {
      next
    }
    d <- as.Date(toks[(length(toks) - 1L):length(toks)], format = "%Y%m%d")
    starts <- c(starts, d[1])
    ends <- c(ends, d[2])
  }
  if (!length(starts)) {
    return(NULL)
  }
  list(start = max(starts), end = min(ends)) # intersection
}

#' Path to the seam-agnostic ENSO continuation library (whole Lenssen piControl
#' series + per-month analog features), built by
#' 03_prepare_enso_continuation_library.R and read by
#' enso.R::load_enso_continuation_library().
#'
#' @param cache_dir Directory (default "generated_data/").
#' @return Character path: <cache_dir>/enso_continuation_library.qs.
make_enso_continuation_library_file <- function(cache_dir = "generated_data/") {
  file.path(cache_dir, "enso_continuation_library.qs")
}

#' Path to the skill-by-lead table of the ENSO analog continuation (fair CRPS,
#' CRPSS vs climatology / neutral on observed seams), written by
#' `03_prepare_enso_continuation_library.R --validate`
#' (utils/enso_continuation_skill.R); its provenance sits beside it as
#' enso_continuation_skill_meta.yml.
#'
#' @param cache_dir Directory (default "generated_data/").
#' @return Character path: <cache_dir>/enso_continuation_skill.csv.
make_enso_continuation_skill_file <- function(cache_dir = "generated_data/") {
  file.path(cache_dir, "enso_continuation_skill.csv")
}

#' Path to the cached GADM ADM1 centroid distance matrix (km).
#'
#' Rho-independent: the kernel `exp(-d / rho)` is cheap to recompute, so only
#' the distance matrix itself is cached (see utils/spatial_contagion.R).
#'
#' @param cache_dir Directory to write/read the cache file (default "generated_data/").
#' @return Character path: <cache_dir>/adm1_centroid_dist_km.rds.
make_adm1_dist_file <- function(cache_dir = "generated_data/") {
  file.path(cache_dir, "adm1_centroid_dist_km.rds")
}

# ---------------------------------------------------------------------------
# Neighbour-contagion covariate naming -- single source of truth
# ---------------------------------------------------------------------------
# The occurrence track's distance-weighted neighbour-contagion covariate is a
# DYNAMIC column name, nbr_lag1_r<rho>[_rn][_pg]: `_rn` is the row-normalised
# (share-of-neighbours, bounded in [0, 1]) sibling and `_pg` the
# partition-restricted sibling (neighbour set AND normalising denominator
# limited to the unit's grouping partition, e.g. its ENSO group). Prep
# (03_prepare_occ.R via spatial_contagion.R::build_contagion_table()) writes the
# columns, the formula builder (formula.R::build_formula_terms()) selects one,
# make_occ_run_file() below derives the artefact token, and the simulator kernels
# (simulator_core.R::evaluate_lp(), spatpomp_core_c.R::precompute_nbr_coefs(),
# spatial_contagion.R::resolve_nbr_spec()) recognise the fixed effect by
# pattern. Every one of them goes through these four helpers, so the convention
# lives in exactly one place. NULL `scale` / `partition` mean the defaults
# ("raw" / "none"), preserving byte-identical names for every pre-existing
# config. Call-time dependency only: callers source paths.R themselves
# (utils/README.md dependency table).
# partition = "registry": neighbour set and denominator limited to
# the unit's group in an ARBITRARY registry `partition_source`, named in the
# covariate itself as `_part<token>` (token = lower-case alphanumeric stem of
# the registry name, config.R::grouping_label_token()), so the simulator
# recovers the registry from the fit's fixed effect
# (config.R::grouping_source_from_token()). `_pg` stays the fit-grouping
# partition, byte-identical.
.nbr_partition_token <- function(partition_source) {
  if (is.null(partition_source) || !nzchar(partition_source)) {
    stop(
      "contagion.partition = 'registry' requires contagion.partition_source ",
      "(a registry under analysis/groupings/).",
      call. = FALSE
    )
  }
  tolower(gsub("[^A-Za-z0-9]", "", partition_source))
}

nbr_covariate_name <- function(
  rho_km = 500L,
  scale = NULL,
  partition = NULL,
  partition_source = NULL
) {
  sprintf(
    "nbr_lag1_r%d%s%s",
    as.integer(rho_km),
    if (identical(scale, "rownorm")) "_rn" else "",
    if (identical(partition, "group")) {
      "_pg"
    } else if (identical(partition, "registry")) {
      paste0("_part", .nbr_partition_token(partition_source))
    } else {
      ""
    }
  )
}

nbr_covariate_pattern <- function() {
  "^nbr_lag1_r[0-9]+(_rn)?(_pg|_part[a-z0-9]+)?$"
}

#' Name of the contagion-coverage SHARE column (a report column,
#' never a model term): the kernel-weighted share of a unit's neighbours with
#' surveillance coverage in the previous month, over the raw kernel of the same
#' partition as nbr_covariate_name(rho, "raw", partition, source):
#' "nbr_cov_lag1_r<rho>[_pg|_part<tok>]". Never matches nbr_covariate_pattern().
nbr_coverage_share_name <- function(
  rho_km = 500L,
  partition = NULL,
  partition_source = NULL
) {
  sub(
    "^nbr_lag1_",
    "nbr_cov_lag1_",
    nbr_covariate_name(rho_km, "raw", partition, partition_source)
  )
}

#' Sidecar of an occurrence panel built with 03 --contagion_coverage_report:
#' `<panel>.contagion_coverage.csv` (ADM1 x month surveillance coverage).
contagion_coverage_csv_path <- function(surv_path) {
  sub("\\.rds$", ".contagion_coverage.csv", as.character(surv_path))
}

#' Inverse of nbr_covariate_name() for ONE name; NULL when it is not an nbr term.
#' `partition_token` is the registry token of a `_part<token>` term (else NULL).
parse_nbr_covariate <- function(name) {
  m <- regmatches(
    name,
    regexec("^nbr_lag1_r([0-9]+)(_rn)?(_pg|_part([a-z0-9]+))?$", name)
  )[[1L]]
  if (!length(m)) {
    return(NULL)
  }
  list(
    rho_km = as.integer(m[2L]),
    scale = if (nzchar(m[3L])) "rownorm" else "raw",
    partition = if (identical(m[4L], "_pg")) {
      "group"
    } else if (nzchar(m[4L])) {
      "registry"
    } else {
      "none"
    },
    partition_token = if (nzchar(m[5L])) m[5L] else NULL
  )
}

nbr_artefact_token <- function(
  rho_km = 500L,
  scale = NULL,
  partition = NULL,
  partition_source = NULL
) {
  sprintf(
    "_nbr%d%s%s",
    as.integer(rho_km),
    if (identical(scale, "rownorm")) "n" else "",
    if (identical(partition, "group")) {
      "g"
    } else if (identical(partition, "registry")) {
      paste0("_part-", .nbr_partition_token(partition_source))
    } else {
      ""
    }
  )
}

# Config-driven filename builder -----------------------------------------

#' Deterministic filename builder for cholera pipeline outputs (occurrence track).
#'
#' @param opt resolved options list from make_options_from_config()
#' @param kind one of "fit", "re_stats", "coef_summary", "loo", "compare",
#'   "cb_meta", "cb_surfaces"
#' @param outdir output directory (defaults to opt$outdir or "generated_data/")
#' @return character path
#' @details Tokens are appended only when config values deviate from the
#'   canonical defaults (e.g. non-default basehaz, non-default covariate set).
#'
#'   Naming convention:
#'   `outbreaks_survival_{backend}_{kind}_base-{baseline}[_bh-{basehaz}][_ctry-{country}][_ext].rds`
#'
#'   where `{ext}` captures any non-default parameter deviations as short tokens.
make_occ_run_file <- function(
  opt,
  kind = c(
    "fit",
    "re_stats",
    "coef_summary",
    "loo",
    "compare",
    "cb_meta",
    "cb_surfaces",
    "occ_base_rate",
    "kl_meta",
    "kl_summary",
    "kl_resp",
    "kl_comembership",
    "kl_partition",
    "kl_pin",
    "kl_pin_plan"
  ),
  outdir = NULL
) {
  kind <- match.arg(kind)

  # Global baseline-rate anchor table: run-invariant by construction, so it
  # takes none of the per-run tokens assembled below. See `.base_rate_path()`.
  if (kind == "occ_base_rate") {
    return(.base_rate_path(
      opt,
      stem = "base_rate_table_occ",
      estimator = opt$occurrence$base_rate$estimator %||% "pooled",
      pre_tokens = paste0("base-", opt$baseline %||% "1990_2019"),
      pop_token = FALSE,
      suffix = paste0(
        .base_rate_flag_token(opt, "occurrence"),
        .occ_panel_variant_token(opt)
      )
    ))
  }

  if (is.null(outdir)) {
    outdir <- if (!is.null(opt$outdir)) opt$outdir else "generated_data/"
  }

  # Canonical defaults (must mirror get_default_options())
  DEFAULT_BACKEND <- "brms_cb"
  DEFAULT_BASEHAZ_POISSON <- "seasonal"
  DEFAULT_COUNTRY <- "all"

  # Read from the nested occurrence block (new schema)
  occ <- opt$occurrence %||% list()
  backend <- occ$backend %||% DEFAULT_BACKEND
  baseline <- if (!is.null(opt$baseline)) opt$baseline else "1990_2019"
  country <- if (!is.null(opt$country)) opt$country else DEFAULT_COUNTRY

  # Backend label used in filenames (matches existing convention)
  backend_label <- switch(
    backend,
    brms = "poisson",
    backend # fallback: use raw value
  )

  # Basehaz token — only append when non-default (preserves legacy names)
  basehaz_token <- ""
  if (backend == "brms") {
    bh <- occ$brms$basehaz_poisson %||% DEFAULT_BASEHAZ_POISSON
    if (bh != DEFAULT_BASEHAZ_POISSON) {
      basehaz_token <- paste0("_bh-", bh)
    } else {
      basehaz_token <- paste0("_bh-", bh)
    }
  } else if (backend %in% c("brms_cb", "brms_kl")) {
    bh <- occ[[backend]]$basehaz_poisson %||%
      occ$brms$basehaz_poisson %||%
      DEFAULT_BASEHAZ_POISSON
    basehaz_token <- paste0("_bh-", bh) # always appended (same convention as brms)
  }

  # Country token — always appended (omitting it caused country=all vs
  # country=ETH collisions)
  country_token <- paste0("_ctry-", country)

  # Grouping token — in country-grouping mode the group identifies the country
  # set, so replace the (long comma-joined) country token with a stable group
  # token. This token is the discovery key used by the grouped simulator.
  grp_token <- ""
  if (isTRUE(opt$grouping$enable)) {
    grp_token <- sprintf("_grp-%s-%s", opt$grouping$source, opt$grouping$group)
    country_token <- ""
  }

  # Admin-level token — always appended
  adm_token <- paste0("_adm-", tolower(occ$admin_level %||% "adm2"))

  res <- occ$random_effects %||% list()

  # brms-specific computational tokens (distinguish runs with different init strategies)
  uncorr_token <- ""
  init_token <- ""
  rinit_token <- ""
  if (backend %in% c("brms", "brms_cb", "brms_kl")) {
    brms_sub <- if (backend != "brms") {
      occ[[backend]] %||% occ$brms %||% list()
    } else {
      occ$brms %||% list()
    }
    if (isTRUE(brms_sub$re_uncorrelated)) {
      uncorr_token <- "_uncorr"
    }
    init_val <- brms_sub$init %||% "random"
    if (!identical(init_val, "random")) {
      init_token <- paste0("_init-", init_val)
    }
    if (!is.null(brms_sub$re_init_r)) {
      rinit_token <- paste0("_rinit-", brms_sub$re_init_r)
    }
  }

  # CB token — appended for brms_cb (and the brms_kl lag-kernel basis) when
  # dlnm is enabled
  cb_token <- ""
  if (backend %in% c("brms_cb", "brms_kl")) {
    dlnm_cfg <- occ$dlnm %||% list()
    cb_vars <- names(dlnm_cfg$variables %||% list())
    if (isTRUE(dlnm_cfg$enable) && length(cb_vars) > 0) {
      cb_token <- paste0("_cb-", paste(gsub("_", "", cb_vars), collapse = "-"))
      # Single-lag onebasis token — disambiguates per-lag fits in a sweep.
      sl <- .dlnm_single_lag(dlnm_cfg)
      if (!is.null(sl)) {
        cb_token <- paste0(cb_token, "_", sl_tag(sl))
      }
      # Multi-lag (`lags:`) basis token, e.g. "_lags-m11to11".
      lg <- .dlnm_lags_range(dlnm_cfg)
      if (!is.null(lg)) {
        cb_token <- paste0(
          cb_token,
          "_lags-",
          sl_token(lg[1L]),
          "to",
          sl_token(lg[2L])
        )
      }
      # Basis-type token — "" for bs (back-compat), "ns4" etc. for natural spline.
      bt <- .dlnm_basis_token(dlnm_cfg)
      if (nzchar(bt)) cb_token <- paste0(cb_token, "_", bt)
    } else {
      cb_token <- "_cb"
    }
  }

  # Lag-kernel token (brms_kl): "_kl" for the single kernel, plus "-K<k>" for
  # latent classes, "-pg<registry>[-kap<kappa>]" for a prior class partition
  # and "-sp" for the spatial (ICAR) membership term.
  kl_token <- if (backend == "brms_kl") {
    .kl_token(occ$kernel %||% list())
  } else {
    ""
  }

  # Non-ENSO lag-1 RE slopes (e.g. mean_precip_lag1_slope, dev_temp_daymean_lag1_slope, ...)
  # Appended when any lag-1 RE slope key (excluding ENSO lags) is active.
  lag1_keys <- grep("^re_adm1_.*_lag1_slope$", names(res), value = TRUE)
  lag1_token <- if (
    length(lag1_keys) > 0 &&
      any(sapply(lag1_keys, function(k) isTRUE(res[[k]])))
  ) {
    "_lag1"
  } else {
    ""
  }

  # ENSO token — appended whenever occurrence.covariates.enso is TRUE
  enso_token <- ""
  if (isTRUE(occ$covariates$enso)) {
    enso_src <- opt$enso$source %||% "oni"
    enso_lag <- as.integer(opt$enso$lag_months %||% 11L)
    enso_slope_key_rc <- sprintf("re_region_country_enso_lag%d_slope", enso_lag)
    enso_slope_key_adm1 <- sprintf("re_adm1_enso_lag%d_slope", enso_lag)
    slope_sfx <- if (
      isTRUE(res[[enso_slope_key_rc]]) ||
        isTRUE(res[[enso_slope_key_adm1]])
    ) {
      "-slope"
    } else {
      ""
    }
    enso_token <- sprintf("_enso-%slag%d%s", enso_src, enso_lag, slope_sfx)
  }

  # noloc token — ENSO is active but all 4 core local climate slope REs are off.
  # Distinguishes ENSO-only configs from ENSO + local climate configs.
  noloc_token <- ""
  if (isTRUE(occ$covariates$enso)) {
    local_slope_keys <- c(
      "re_adm1_mean_precip_slope",
      "re_adm1_dev_precip_slope",
      "re_adm1_mean_temp_daymean_slope",
      "re_adm1_dev_temp_daymean_slope"
    )
    if (!any(sapply(local_slope_keys, function(k) isTRUE(res[[k]])))) {
      noloc_token <- "_noloc"
    }
  }

  # Distance-weighted neighbour-contagion token. MUST be non-empty whenever the
  # covariate is on: two runs can differ by this covariate alone, so
  # without a token every fit/loo/coef_summary/re_stats path would be byte-identical
  # to its covariate-off sibling — overwriting archived fits, and (worse) tripping the
  # fit scripts' file.exists() early-exit so the new run would mostly not execute.
  # Empty for every pre-existing config (the covariate defaults FALSE), so no
  # existing artefact name changes.
  nbr_token <- ""
  if (isTRUE(occ$covariates$nbr_contagion_lag1)) {
    nbr_token <- nbr_artefact_token(
      opt$contagion$rho_km %||% 500L,
      opt$contagion$scale,
      opt$contagion$partition,
      opt$contagion$partition_source
    )
  }

  # Non-standard ADM1 climate slope token — fires when any re_adm1_mean_*_slope
  # beyond the two canonical defaults (mean_precip, mean_temp_daymean) is active.
  # Uses mean_ keys only (dev_ always pairs with mean_); vars sorted for stability.
  # Excludes lag-1 and ENSO keys (handled by lag1_token and enso_token).
  STD_ADM1_CLIM_SLOPES <- c(
    "re_adm1_mean_precip_slope",
    "re_adm1_mean_temp_daymean_slope"
  )
  nonstd_adm1_keys <- sort(grep("^re_adm1_mean_", names(res), value = TRUE))
  nonstd_adm1_keys <- nonstd_adm1_keys[
    !nonstd_adm1_keys %in% STD_ADM1_CLIM_SLOPES &
      !grepl("_lag1_slope$|_enso_", nonstd_adm1_keys) &
      vapply(nonstd_adm1_keys, function(k) isTRUE(res[[k]]), logical(1L))
  ]
  nonstd_clim_token <- if (length(nonstd_adm1_keys) > 0) {
    vars <- gsub("^re_adm1_mean_(.+)_slope$", "\\1", nonstd_adm1_keys)
    paste0("_clim-", paste(gsub("_", "", vars), collapse = "-"))
  } else {
    ""
  }

  # Non-default RE structural tokens
  adm2i_token <- if (isTRUE(res$re_adm2_intercept)) "_adm2i" else ""
  bhslp_token <- if (isTRUE(res$re_region_country_basehaz_slope)) {
    "_bhslp"
  } else {
    ""
  }

  # Kind suffix
  kind_label <- switch(
    kind,
    fit = "fit",
    re_stats = "re_stats",
    coef_summary = "coef_summary",
    loo = "loo",
    compare = "compare",
    cb_meta = "cb_meta",
    cb_surfaces = "cb_surfaces",
    kl_meta = "kl_meta",
    kl_summary = "kl_summary",
    kl_resp = "kl_resp",
    kl_comembership = "kl_comembership",
    kl_partition = "kl_partition",
    kl_pin = "kl_pin",
    kl_pin_plan = "kl_pin_plan"
  )

  fname <- paste0(
    "outbreaks_survival_",
    backend_label,
    "_",
    kind_label,
    "_base-",
    baseline,
    basehaz_token,
    country_token,
    grp_token,
    adm_token,
    uncorr_token,
    init_token,
    rinit_token,
    enso_token,
    noloc_token,
    nbr_token,
    nonstd_clim_token,
    lag1_token,
    adm2i_token,
    bhslp_token,
    cb_token,
    kl_token,
    # occurrence.ttlo_no_prior_level: the TTLO level "none" is
    # a fit-level choice (load-time recode), so a fit token, before the panel's.
    if (isTRUE(occ[["ttlo_no_prior_level"]])) "_ttn" else "",
    .occ_panel_variant_token(opt),
    ".rds"
  )

  assert_basename_length(
    file.path(outdir, .cap_file_name(fname)),
    "make_occ_run_file"
  )
}

# File-name length limit -----------------------------------------------------

#' NAME_MAX: the longest file or directory name, in bytes, the filesystems the
#' pipeline writes to accept (local ext4, cluster scratch). The one
#' place the limit is stated; `.cap_file_name()` and `assert_basename_length()` read it.
BASENAME_MAX_BYTES <- 255L

#' Stop when a built path has a file or directory name over BASENAME_MAX_BYTES.
#'
#' Every name builder in this file returns through it, so a name that cannot be
#' created fails where it is built instead of at saveRDS() after a fit.
#' The message names the builder, the length, the
#' limit and the longest `_`-separated tokens of the name, the ones to shorten.
#'
#' @param path character vector of paths (every element and component checked)
#' @param builder name of the calling builder, for the message
#' @param max_bytes limit [BASENAME_MAX_BYTES]
#' @return `path`, unchanged (class kept, e.g. glue)
assert_basename_length <- function(
  path,
  builder,
  max_bytes = BASENAME_MAX_BYTES
) {
  for (p in as.character(path)) {
    comp <- strsplit(p, "/", fixed = TRUE)[[1L]]
    n <- nchar(comp, type = "bytes")
    if (!any(n > max_bytes)) {
      next
    }
    nm <- comp[which.max(n)]
    tok <- strsplit(sub("\\.[A-Za-z0-9]+$", "", nm), "_", fixed = TRUE)[[1L]]
    tl <- nchar(tok, type = "bytes")
    top <- head(order(-tl), 5L)
    stop(
      sprintf(
        paste0(
          "%s: file name is %d bytes, over the %d-byte limit by %d ",
          "(BASENAME_MAX_BYTES, utils/paths.R):\n  %s\n  longest tokens: %s"
        ),
        builder,
        max(n),
        max_bytes,
        max(n) - max_bytes,
        nm,
        paste(sprintf("%s (%d)", tok[top], tl[top]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  path
}

#' Every artefact name a fit config's track can emit, built up front.
#'
#' The post-fit kinds (cb_meta, re_stats, coef_summary, cb_surfaces, pp_check)
#' are longer than the fit's own name, so checking the fit name alone lets a
#' config through that dies after its MCMC. The 04 fit workers call this before
#' fitting and `make_batch_configs.R::generate_batch()` on every config it
#' writes. The run-invariant anchor tables (`*_base_rate`, `size_mean_rate`)
#' are left out: short fixed stems, and their population token needs
#' provenance.R, which not every worker sources.
#'
#' @param opt resolved options
#' @param track occurrence | duration | monthly | size (or occ / dur / mon / siz)
#' @return invisibly, the names built (stops via assert_basename_length() first)
check_run_file_names <- function(opt, track) {
  builders <- c(
    occurrence = "make_occ_run_file",
    duration = "make_duration_run_file",
    monthly = "make_monthly_run_file",
    size = "make_size_run_file",
    occ = "make_occ_run_file",
    dur = "make_duration_run_file",
    mon = "make_monthly_run_file",
    siz = "make_size_run_file"
  )
  if (length(track) != 1L || !track %in% names(builders)) {
    stop(
      "check_run_file_names: unknown track '",
      paste(track, collapse = ","),
      "'.",
      call. = FALSE
    )
  }
  b <- get(builders[[track]])
  kinds <- eval(formals(b)$kind)
  kinds <- kinds[!grepl("(base|mean)_rate$", kinds)]
  invisible(vapply(kinds, function(k) as.character(b(opt, kind = k)), ""))
}

#' A file name the filesystem accepts (NAME_MAX, 255 bytes): returned unchanged
#' when it fits; otherwise its first characters plus `_h<xxhash64 of the full
#' name>` and the extension, so it is deterministic and distinct per full name.
#' No name that fits is altered, so no file already on disk is renamed.
#' Occurrence only (make_occ_run_file()); every other builder stops instead
#' (assert_basename_length()).
.cap_file_name <- function(fname, max_bytes = BASENAME_MAX_BYTES) {
  if (nchar(fname, type = "bytes") <= max_bytes) {
    return(fname)
  }
  ext <- regmatches(fname, regexpr("\\.[A-Za-z0-9]+$", fname))
  hash <- paste0(
    "_h",
    digest::digest(fname, algo = "xxhash64", serialize = FALSE)
  )
  paste0(substr(fname, 1L, max_bytes - nchar(hash) - nchar(ext)), hash, ext)
}


#' Build the file path for the simulation projection climate panel.
#'
#' Analogous to make_sim_panel_file(); used by 03d_prepare_sim_panel_projection.R
#' and by 09_simulator.R when the config specifies a projection panel.
#'
#' @param opt List with at least `opt$baseline` (e.g. "1990_2019").
#' @return Character string: path to the projection panel RDS.
make_sim_panel_projection_file <- function(opt) {
  assert_basename_length(
    glue::glue(
      "generated_data/sim_climate_panel_projection_base-{opt$baseline}",
      "{.episode_key_path_token(opt)}.rds"
    ),
    "make_sim_panel_projection_file"
  )
}

#' Path to the historical simulator panel extended past the climate DB to a
#' given month (03_prepare_sim_panel_extended.R): the historical panel's rows
#' through its last month, then the projection panel's rows up to `end`.
#'
#' @param opt Resolved options (`baseline`, `episodes`).
#' @param end Last month included (Date or "YYYY-MM-DD"); named in the file as
#'   `_to-YYYY-MM` so panels of different extents never overwrite each other.
#' @return Character path: generated_data/sim_climate_panel_base-<baseline>[_ep-gs]_to-YYYY-MM.rds.
make_sim_panel_extended_file <- function(opt, end) {
  assert_basename_length(
    glue::glue(
      "generated_data/sim_climate_panel_base-{opt$baseline}",
      "{.episode_key_path_token(opt)}_to-{format(as.Date(end), '%Y-%m')}.rds"
    ),
    "make_sim_panel_extended_file"
  )
}


# ============================================================================
# Outbreak-size scaffold
# ============================================================================

# Estimand-coherence path token.
# `build_mean_rate_table()`/`read_size_data()` behavior depends on
# opt$size$estimand$denominator, opt$size$pop_elasticity,
# opt$size$canonical_admin_level, opt$size$mean_rate$level_factor, and
# opt$size$mean_rate$cross_absent_levels, so they are encoded into every size
# output path: otherwise two configs differing only in these fields would
# collide on one filename, and persist_mean_rate_table()'s `if (!file.exists(path))`
# skip would silently reuse the first config's table. Reads raw opt$size fields directly (not
# resolve_size_estimand()) per the "modules never source() each other" rule --
# paths.R must not depend on provenance.R's validation/defaulting logic.
# Deployed-default config (denominator=footprint, pop_elasticity=fixed,
# canonical_admin_level unset/FALSE, level_factor=TRUE, cross_absent_levels=TRUE)
# resolves to "" -- bit-inert, zero path changes for every pre-existing config.
#
# Track-parameterised: the monthly track carries the same estimand keys
# (`monthly.estimand`, `monthly.pop_elasticity`, `monthly.canonical_admin_level`
# -- see `.ESTIMAND_SPEC` in utils/provenance.R). The `nolf`/`noxal` pair is
# emitted for SIZE ONLY: those two flags live under `size$mean_rate`, whereas
# the other three tracks carry their equivalents under `<block>$base_rate` and
# already get them from `.base_rate_flag_token(opt, block)`. Emitting both here
# would double-encode them into one filename.
.estimand_path_token <- function(opt, block = "size") {
  blk <- opt[[block]] %||% list()
  parts <- c()

  denom <- blk$estimand$denominator %||% "footprint"
  if (!identical(denom, "footprint")) {
    parts <- c(parts, denom)
  }

  urm <- blk$estimand$underreporting_multiplier %||% 1.0
  if (!isTRUE(all.equal(urm, 1.0))) {
    parts <- c(parts, paste0("urm", urm))
  }

  pop_elas <- blk$pop_elasticity %||% "fixed"
  if (!identical(pop_elas, "fixed")) {
    parts <- c(
      parts,
      switch(pop_elas, free = "free", free_by_level = "freelvl", pop_elas)
    )
  }

  if (isTRUE(blk$canonical_admin_level)) {
    parts <- c(parts, "canon")
  }

  if (identical(block, "size")) {
    if (!isTRUE(blk$mean_rate$level_factor %||% TRUE)) {
      parts <- c(parts, "nolf")
    }
    if (!isTRUE(blk$mean_rate$cross_absent_levels %||% TRUE)) {
      parts <- c(parts, "noxal")
    }
  }

  # Prefix "_estnd-" (not "_est-") -- the size_mean_rate path already uses
  # "_est-<estimator>" for opt$size$mean_rate$estimator (e.g. "geometric");
  # reusing that prefix here would make two structurally different tokens
  # look like variants of the same field.
  if (length(parts) == 0) {
    ""
  } else {
    paste0("_estnd-", paste(parts, collapse = "-"))
  }
}

# Original name kept per the side-by-side rule: it is the documented token in
# every archived size artefact name.
.size_estimand_path_token <- function(opt) .estimand_path_token(opt, "size")

# Shared path builder for the baseline-rate anchor sidecars (`utils/anchor.R`).
#
# Every track's anchor table is GLOBAL -- built once over the whole corpus, then
# shared by every country/group/backend config -- so its filename must carry the
# inputs that change its CONTENT and none of the run-scoped tokens the rest of
# this file assembles. Get that wrong in either direction and the failure is
# silent: run-scoped tokens fragment one table into many (each a different
# subset, none reproducible), while a missing content token makes two genuinely
# different tables collide on one name, whereupon `persist_base_rate_table()`'s
# `if (!file.exists(path))` skip hands the second run the first run's table.
#
# Content-determining inputs, in the order they appear in the name:
#   `pre_tokens`  per-track quality guards (size: min_pop / max_attack_rate)
#   `est-`        the estimator (geometric | arithmetic | pooled)
#   `pop-<hash>`  the population digest -- ONLY for tracks whose denominator
#                 contains `pop` (mon, siz). occ/dur divide by exposure days
#                 alone, so no population can change their table and a pop token
#                 would only invalidate their cache spuriously.
#   `suffix`      remaining per-track flags (size: `.size_estimand_path_token()`,
#                 which already encodes `nolf`/`noxal`; the new tracks get the
#                 same two flags via `.base_rate_flag_token()`).
#
# `make_size_run_file(opt, "size_mean_rate")` delegates here and its filename is
# byte-identical to the hand-rolled version it replaces.
.base_rate_path <- function(
  opt,
  stem,
  estimator,
  pre_tokens = character(0),
  pop_token = FALSE,
  suffix = ""
) {
  base <- paste(c(stem, pre_tokens, paste0("est-", estimator)), collapse = "_")
  if (isTRUE(pop_token)) {
    tok <- current_pop_digest_token(opt)
    if (!is.na(tok)) base <- paste0(base, "_", tok)
  }
  assert_basename_length(
    file.path(opt$outdir %||% "generated_data/", paste0(base, suffix, ".rds")),
    paste0(".base_rate_path(", stem, ")")
  )
}

# The occ/dur/mon analogue of `.size_estimand_path_token()`'s `nolf`/`noxal` pair:
# `level_factor` and `cross_absent_levels` both change every cell of the table,
# and both default TRUE, so this resolves to "" for every config that exists.
# `random_effects.re_adm2_leaf_scope` (size / monthly / duration): "all" (default, the
# historical `(1 | .../ADM1/<gid>)` leaf over every unit) | "adm2_only" (leaf
# only on canonical-ADM2 units; see `.intercept_re_term()` in utils/config.R).
# Lives here because the fit filename must change with it -- two fits that
# differ only in leaf structure would otherwise collide on the fit scripts'
# file.exists() early exit.
RE_LEAF_SCOPES <- c("all", "adm2_only")
.leaf_scope_tag <- function(blk) {
  re <- (blk %||% list())$random_effects %||% list()
  s <- re$re_adm2_leaf_scope %||% "all"
  if (!s %in% RE_LEAF_SCOPES) {
    stop(
      "random_effects.re_adm2_leaf_scope = '",
      s,
      "' is not one of {",
      paste(RE_LEAF_SCOPES, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
  s
}

# Fit-name token of a block's leaf scope ("leaf-adm2" | none), shared by the
# size / monthly / duration run-file builders.
.leaf_token <- function(blk) {
  if (identical(.leaf_scope_tag(blk), "adm2_only")) {
    "leaf-adm2"
  } else {
    character(0)
  }
}

# Loader side of a block's leaf scope, shared by read_size_data(),
# load_monthly_data() and load_duration_data(): a no-op under the default
# `all`, else utils/formula.R's .leaf_columns() (call-time only, so a loader
# needs formula.R sourced only when the leaf is on) with the canonical-level
# check satisfied by the canonical re-key (`admin_level_canonical`) or a gid_season panel.
# `ep_rec` is the loader's assert_episode_key() record.
.leaf_columns_for <- function(d, blk, ep_rec, tag = "leaf") {
  if (identical(.leaf_scope_tag(blk), "all")) {
    return(d)
  }
  .leaf_columns(
    d,
    canonical_ok = "admin_level_canonical" %in%
      names(d) ||
      identical(ep_rec$key, "gid_season"),
    tag = tag
  )
}

.base_rate_flag_token <- function(opt, block) {
  br <- ((opt[[block]] %||% list())$base_rate %||% list())
  parts <- c()
  if (!isTRUE(br$level_factor %||% TRUE)) {
    parts <- c(parts, "nolf")
  }
  if (!isTRUE(br$cross_absent_levels %||% TRUE)) {
    parts <- c(parts, "noxal")
  }
  if (length(parts) == 0) {
    ""
  } else {
    paste0("_brflg-", paste(parts, collapse = "-"))
  }
}

#' Resolve the single-lag index for a DLNM config block.
#'
#' The onebasis single-lag fit must carry a `sl<N>` token else every lag in a
#' sweep overwrites the same file. The lag is read from `defaults$single_lag`,
#' falling back to the first per-variable `single_lag`. Factored out of the
#' three runfile builders (occ/dur/mon), which each inlined this verbatim.
.dlnm_single_lag <- function(dlnm_cfg) {
  sl <- dlnm_cfg$defaults$single_lag
  if (is.null(sl)) {
    for (vv in names(dlnm_cfg$variables %||% list())) {
      if (!is.null(dlnm_cfg$variables[[vv]]$single_lag)) {
        sl <- dlnm_cfg$variables[[vv]]$single_lag
        break
      }
    }
  }
  sl
}

#' Exposure-variable token of a duration / monthly crossbasis fit name:
#' "" when the fit carries an ENSO basis (every deployed fit, including the
#' monthly variants with extra bases, so their names are unchanged), else
#' "_cb-<variables sans '_'>" -- the occurrence builder's `_cb-<vars>` token --
#' so a local-climate fit (which replaces ENSO) never shares a name with the
#' ENSO fit of the same lag. `mso` (the monthly months-since-onset basis, never
#' an exposure) is left out.
.dlnm_vars_token <- function(dlnm_cfg) {
  v <- setdiff(names(dlnm_cfg$variables %||% list()), "mso")
  if (!length(v) || "enso" %in% v) {
    return("")
  }
  paste0("_cb-", paste(gsub("_", "", v), collapse = "-"))
}

#' Lag range of a multi-lag (`lags:`) DLNM spec: c(min, max) from the shared
#' defaults or the first variable that sets it; NULL when the mode is off.
.dlnm_lags_range <- function(dlnm_cfg) {
  lg <- dlnm_cfg$defaults$lags
  if (is.null(lg)) {
    for (vv in names(dlnm_cfg$variables %||% list())) {
      if (!is.null(dlnm_cfg$variables[[vv]]$lags)) {
        lg <- dlnm_cfg$variables[[vv]]$lags
        break
      }
    }
  }
  if (is.null(lg)) {
    return(NULL)
  }
  range(as.integer(unlist(lg)))
}

#' Lag-kernel (brms_kl) filename token from the occurrence `kernel:` block.
.kl_token <- function(kcfg) {
  tok <- "_kl"
  if (identical(kcfg$engine, "brms")) {
    tok <- paste0(tok, "-nl")
  }
  K <- as.integer(kcfg$classes %||% 1L)
  if (K >= 2L) {
    tok <- paste0(tok, "-K", K)
  }
  pg <- kcfg$prior_grouping$source
  if (!is.null(pg) && nzchar(pg)) {
    tok <- paste0(tok, "-pg", gsub("_", "", pg))
    kap <- kcfg$prior_grouping[["kappa"]] # [[ ]]: `$` would partially match kappa_prior_sd
    if (!is.null(kap)) {
      tok <- paste0(tok, "-kap", kap)
    }
    if (isTRUE(kcfg$prior_grouping[["membership"]])) {
      tok <- paste0(tok, "-mem")
    }
    cg <- kcfg$prior_grouping[["contiguity"]]
    if (isTRUE(cg[["enable"]])) {
      tok <- paste0(
        tok,
        "-cmix",
        gsub(
          ".",
          "p",
          format(as.numeric(cg[["weight"]] %||% 0.5), drop0trailing = TRUE),
          fixed = TRUE
        )
      )
    }
  }
  if (isTRUE(kcfg$spatial$enable)) {
    tok <- paste0(
      tok,
      "-sp",
      if (identical(kcfg$spatial$prior, "dagar")) "dagar" else ""
    )
  }
  if (identical(kcfg$form, "dlnm")) {
    tok <- paste0(
      tok,
      "-dlnm",
      kcfg$lag_df %||% 10L,
      if (identical(kcfg$lag_prior, "normal")) "nm" else "hs"
    )
  }
  ksdf <- kcfg[["ksd_fixed"]]
  if (
    !is.null(ksdf) &&
      !identical(kcfg$form, "discrete") &&
      !identical(kcfg$form, "dlnm")
  ) {
    tok <- paste0(
      tok,
      "-ksd",
      gsub(
        ".",
        "p",
        format(as.numeric(ksdf), drop0trailing = TRUE),
        fixed = TRUE
      )
    )
  }
  if (isTRUE(kcfg$re_slope$enable) && !identical(kcfg$form, "discrete")) {
    # -rea = ADM1 deviations only, -rec = country deviations only (re_slope.adm1: false)
    tok <- paste0(
      tok,
      if (isFALSE(kcfg$re_slope[["country"]])) {
        "-rea"
      } else if (isFALSE(kcfg$re_slope[["adm1"]])) {
        "-rec"
      } else {
        "-re"
      }
    )
    sdf <- kcfg$re_slope[["sd_country_fixed"]]
    if (!isFALSE(kcfg$re_slope[["country"]]) && !is.null(sdf)) {
      tok <- paste0(
        tok,
        "-sdc",
        gsub(
          ".",
          "p",
          format(as.numeric(sdf), drop0trailing = TRUE),
          fixed = TRUE
        )
      )
    }
    # half-normal prior scale of the estimated country sd; the default 0.35 adds
    # no token
    sdp <- kcfg$re_slope[["sd_country"]]
    if (
      !isFALSE(kcfg$re_slope[["country"]]) &&
        is.null(sdf) &&
        !is.null(sdp) &&
        !isTRUE(all.equal(as.numeric(sdp), 0.35))
    ) {
      tok <- paste0(
        tok,
        "-sdcp",
        gsub(
          ".",
          "p",
          format(as.numeric(sdp), drop0trailing = TRUE),
          fixed = TRUE
        )
      )
    }
    sda <- kcfg$re_slope[["sd_adm1"]]
    if (
      !isFALSE(kcfg$re_slope[["adm1"]]) &&
        !is.null(sda) &&
        !isTRUE(all.equal(as.numeric(sda), 0.35))
    ) {
      tok <- paste0(
        tok,
        "-sda",
        gsub(
          ".",
          "p",
          format(as.numeric(sda), drop0trailing = TRUE),
          fixed = TRUE
        )
      )
    }
  }
  if (identical(kcfg$form, "discrete")) {
    dc <- kcfg$discrete %||% list()
    w <- gsub(
      ".",
      "p",
      format(
        as.numeric(unlist(dc$widths %||% 0)),
        trim = TRUE,
        drop0trailing = TRUE
      ),
      fixed = TRUE
    )
    tok <- paste0(tok, "-disc-w", paste(w, collapse = "-"))
    st <- as.integer(dc$centre_step %||% 1L)
    if (st != 1L) {
      tok <- paste0(tok, "-st", st)
    }
    if (identical(dc$combo_prior, "uniform")) {
      tok <- paste0(tok, "-cpu")
    }
    if (isTRUE(kcfg$re_slope$enable)) {
      tok <- paste0(
        tok,
        if (isFALSE(kcfg$re_slope[["country"]])) {
          "-rea"
        } else if (isFALSE(kcfg$re_slope[["adm1"]])) {
          "-rec"
        } else {
          "-re"
        }
      )
    }
    pn <- kcfg$pin %||% list()
    if (isTRUE(pn$enable)) {
      tok <- paste0(
        tok,
        "-pin",
        if (!is.null(pn[["seeds"]])) "x" else as.integer(pn$seeds_n %||% 1L),
        if (identical(pn$arrangements, "all")) "all" else "",
        if (identical(pn$order_source, "sweep")) "sw" else ""
      )
    }
  }
  pb <- kcfg[["prior_beta_sd"]]
  if (!is.null(pb) && !isTRUE(all.equal(as.numeric(pb), 0.5))) {
    tok <- paste0(
      tok,
      "-pb",
      gsub(".", "p", format(as.numeric(pb), drop0trailing = TRUE), fixed = TRUE)
    )
  }
  id <- kcfg$identify
  if (!is.null(id) && !identical(id, "auto")) {
    tok <- paste0(tok, "-id", id)
  }
  seed_rp <- identical(id, "seed") && isTRUE(kcfg$repel[["enable"]])
  if (seed_rp) {
    tok <- paste0(tok, "-rp")
  } # identify seed + the repel prior (kernel.repel.enable)
  if (identical(id, "repel") || seed_rp) {
    # identify repel (or seed + repulsion): the repulsion scales / grid and the
    # class-weight concentration, each only when not at its default
    num <- function(x) {
      gsub(
        ".",
        "p",
        gsub(
          "-",
          "m",
          format(as.numeric(x), drop0trailing = TRUE, trim = TRUE),
          fixed = TRUE
        ),
        fixed = TRUE
      )
    }
    rl <- kcfg$repel[["lag_scale"]]
    rs <- kcfg$repel[["erf_scale"]]
    rg <- kcfg$repel[["erf_grid"]]
    if (!is.null(rl) && !isTRUE(all.equal(as.numeric(rl), 3))) {
      tok <- paste0(tok, "-rl", num(rl))
    }
    if (!is.null(rs) && !isTRUE(all.equal(as.numeric(rs), 0.5))) {
      tok <- paste0(tok, "-rs", num(rs))
    }
    if (
      !is.null(rg) &&
        !isTRUE(all.equal(as.numeric(unlist(rg)), c(-1.5, -0.75, 0.75, 1.5)))
    ) {
      tok <- paste0(tok, "-rg", paste(num(unlist(rg)), collapse = "_"))
    }
    al <- kcfg[["alpha"]]
    if (!is.null(al) && !isTRUE(all.equal(as.numeric(al), 1))) {
      tok <- paste0(tok, "-a", num(al))
    }
  }
  # kernel.lag_re: a lag offset per country (half-normal sd of the offsets;
  # -lre = scale month, lags; -lrl = scale logit, logit units of the lag window)
  lr <- kcfg[["lag_re"]]
  if (isTRUE(lr[["enable"]])) {
    tok <- paste0(
      tok,
      if (identical(lr[["scale"]] %||% "month", "logit")) "-lrl" else "-lre",
      gsub(
        ".",
        "p",
        format(as.numeric(lr[["sd"]] %||% 2), drop0trailing = TRUE),
        fixed = TRUE
      )
    )
  }
  tok
}

#' Basis-type disambiguation token for a DLNM argvar.
#'
#' Returns "" for `bs` (or an absent argvar) so every existing B-spline
#' artefact keeps its byte-identical filename; returns e.g. "ns4" for a
#' natural-spline argvar so an `ns` fit does not overwrite the `bs` one.
.dlnm_basis_token <- function(dlnm_cfg) {
  av <- dlnm_cfg$defaults$argvar %||% list()
  fun <- av$fun %||% "bs"
  if (identical(fun, "bs")) {
    return("")
  }
  paste0(fun, av$df %||% "")
}

#' The size track's `mean_rate` estimator, resolved once for every layer.
#'
#' Three places need this answer, and their fallbacks for an unset key differ:
#'   * `make_size_run_file(opt, "size_mean_rate")` below  -> "geometric"
#'   * `build_mean_rate_table()` (utils/provenance.R)     -> "geometric"
#'   * `build_base_rate_table(opt, "siz")` (utils/anchor.R), which is what
#'     `read_size_data()` actually builds the table with -> `.ANCHOR_SPEC$siz$
#'     estimator`, i.e. "pooled".
#'
#' That disagreement is invisible while the key is set -- which it is for every
#' config in the repo, because `config_defaults.yml` ships
#' `size.mean_rate.estimator: pooled` and the merge injects the whole
#' `mean_rate` sub-block when a config does not write one. It stops being
#' invisible the moment a config writes a PARTIAL `mean_rate:` block (say
#' `level_factor: false` alone): the merge replaces that sub-block wholesale
#' rather than filling in the missing key, so `estimator` comes out genuinely
#' NULL, the table is built `pooled`, and the path it is looked up under says
#' `est-geometric` -- where `persist_base_rate_table()`'s `if (!file.exists(path))`
#' skip finds the deployed geometric table and loads THAT. A silently wrong
#' anchor on every row, not a crash.
#'
#' So the partial-block case is a hard error rather than a fourth guess: there
#' is no defensible default for a config that configured this block and then
#' left out the one field that changes every cell of the table.
#'
#' The no-`size`-block case keeps returning "geometric", and must. It is not a
#' modelling default -- it is the NAME of the already-deployed sidecar, which
#' the simulator/smoothing entry points (`07a_simulator_tasks.R:90`,
#' `07c_simulator_postprocess.R:1818`) resolve from their own simulator `opt`,
#' an object that has no `size:` block at all. 567 configs reach the size
#' mean-rate path that way. Nothing is being built on that route -- it is pure
#' path algebra against an artefact that already exists on disk -- so the
#' estimator disagreement cannot arise there.
#'
#' @param opt Resolved options list (or any list; `opt$size` may be absent).
#' @return "geometric" | "arithmetic" | "pooled"
size_mean_rate_estimator <- function(opt) {
  mr <- (opt$size %||% list())$mean_rate
  if (is.null(mr)) {
    return("geometric")
  } # no size block: legacy sidecar name
  if (is.null(mr$estimator)) {
    stop(
      "size.mean_rate is configured but size.mean_rate.estimator is not. ",
      "A partial mean_rate: block REPLACES the config_defaults.yml block ",
      "rather than merging into it, so the estimator does not fall through ",
      "to its default -- and the table builder and the path builder would ",
      "then disagree about which estimator this run uses. Set ",
      "size.mean_rate.estimator explicitly (pooled | geometric | ",
      "arithmetic), or drop the mean_rate: block entirely to inherit the ",
      "default.",
      call. = FALSE
    )
  }
  mr$estimator
}

#' Deterministic filename builder for the outbreak-size track.
#'
#' @param opt Resolved options list from make_options_from_config().
#' @param kind One of "size_fit", "size_re_stats", "size_coef_summary",
#'   "size_loo", "size_pp_check", "size_cb_meta", "size_cb_surfaces",
#'   "size_mean_rate".
#' @return Character path under `opt$outdir` (default "generated_data/").
#' @details Extends `make_run_file()` informally; encodes the most load-bearing
#'   fields (backend, kind, baseline, country, family, seasonal_term,
#'   admin_level, enso lag).
make_size_run_file <- function(
  opt,
  kind = c(
    "size_fit",
    "size_re_stats",
    "size_coef_summary",
    "size_loo",
    "size_pp_check",
    "size_cb_meta",
    "size_cb_surfaces",
    "size_mean_rate"
  )
) {
  kind <- match.arg(kind)
  # The mean_rate table is built once, globally (country x admin_level, all
  # groups/backends pooled) — it must NOT carry any of the run-scoped tokens
  # (country/group, family, seasonal term, ENSO lag, dlnm lag, phi_re) that the
  # rest of this function assembles below. Every fit config that shares the
  # same (min_pop, max_attack_rate, estimator) triple must resolve to the
  # identical path, so this is a special-cased early return.
  #
  # Population is ALSO content-addressed into this path (`_pop-<digest>`),
  # via `current_pop_digest_token()` above -- see that helper's doc comment for
  # the fragility this closes (a legacy- and a crosswalk-built table colliding
  # on one filename, so a never-rewritten sidecar pairs new populations with old
  # rates). When no digest sidecar is on disk (older artefact tree), the token is NA
  # and the path degrades to the
  # original, un-addressed filename -- unchanged behaviour for old runs.
  if (kind == "size_mean_rate") {
    return(.base_rate_path(
      opt,
      stem = "mean_rate_table",
      estimator = size_mean_rate_estimator(opt),
      pre_tokens = c(
        sprintf("minpop-%s", opt$size$min_pop %||% 10),
        sprintf("maxar-%s", opt$size$max_attack_rate %||% 0.10)
      ),
      pop_token = TRUE,
      suffix = paste0(
        .size_estimand_path_token(opt),
        .episode_key_path_token(opt)
      )
    ))
  }
  # Country-grouping mode: replace the (long comma-joined) country token with a
  # stable group token — the discovery key for the grouped simulator.
  if (isTRUE(opt$grouping$enable)) {
    loc_token <- sprintf("grp-%s-%s", opt$grouping$source, opt$grouping$group)
  } else {
    loc_token <- paste0("ctry-", opt$country %||% "all")
  }
  parts <- c(
    "outbreak_size",
    opt$size$backend %||% "brms",
    sub("^size_", "", kind),
    paste0("base-", opt$baseline %||% "1990_2019"),
    loc_token,
    paste0("fam-", opt$size$family %||% "negbinomial")
  )
  if (!identical(opt$size$seasonal_term, "none")) {
    parts <- c(parts, paste0("seas-", opt$size$seasonal_term))
  }
  parts <- c(parts, paste0("adm-", opt$size$admin_level %||% "ADM1"))
  if (isTRUE(opt$size$random_effects$re_adm2_intercept)) {
    parts <- c(parts, "adm2i")
  }
  # ENSO toggle now lives under size.covariates.enso (new schema)
  if (isTRUE(opt$size$covariates$enso)) {
    parts <- c(parts, paste0("enso-lag", opt$enso$lag_months))
  }
  # DLNM crossbasis single-lag token — disambiguates per-lag onebasis fits
  # (the distributed-lag cb fit carries no lag token; single-lag must, else
  # all lags in a sweep overwrite the same file).
  dlnm_cfg <- opt$size$dlnm %||% list()
  if (identical(opt$size$backend, "brms_cb") && isTRUE(dlnm_cfg$enable)) {
    sl <- dlnm_cfg$defaults$single_lag
    if (is.null(sl)) {
      for (vv in names(dlnm_cfg$variables %||% list())) {
        if (!is.null(dlnm_cfg$variables[[vv]]$single_lag)) {
          sl <- dlnm_cfg$variables[[vv]]$single_lag
          break
        }
      }
    }
    if (!is.null(sl)) parts <- c(parts, sl_tag(sl))
  }
  phi_re <- opt$size$phi_re %||% "constant"
  if (phi_re != "constant") {
    parts <- c(parts, paste0("phire-", gsub("_", "-", phi_re)))
  }
  parts <- c(parts, .leaf_token(opt$size))
  fname <- paste0(
    paste(parts, collapse = "_"),
    .size_estimand_path_token(opt),
    .episode_key_path_token(opt),
    ".rds"
  )
  assert_basename_length(
    file.path(opt$outdir %||% "generated_data/", fname),
    "make_size_run_file"
  )
}


# ============================================================================
# Outbreak-duration track (cessation hazard)
# ============================================================================

#' Deterministic filename builder for the outbreak-duration track.
#'
#' @param opt Resolved options list from make_options_from_config().
#' @param kind One of "dur_fit", "dur_re_stats", "dur_coef_summary", "dur_loo",
#'   "dur_pp_check", "dur_cb_meta", "dur_cb_surfaces".
#' @return Character path under `opt$outdir` (default "generated_data/").
#' @details Same shape as `make_size_run_file()`: a deterministic path under
#'   `opt$outdir` built from the load-bearing fields.
#'
#'   There is no estimand token here and no population digest: the duration
#'   response carries a TIME offset (log days of exposure), not a population
#'   one, so none of the pop/mean_rate provenance that `make_size_run_file()`
#'   must encode applies. What DOES have to be in the name is `mso_breaks` --
#'   the band set is the primary baseline hazard's entire parameterisation,
#'   `load_duration_data()` rebuilds `mso_band` from config at every load, and
#'   two configs differing only in the breaks would otherwise write the same
#'   file.
make_duration_run_file <- function(
  opt,
  kind = c(
    "dur_fit",
    "dur_re_stats",
    "dur_coef_summary",
    "dur_loo",
    "dur_pp_check",
    "dur_cb_meta",
    "dur_cb_surfaces",
    "dur_base_rate"
  )
) {
  kind <- match.arg(kind)
  dur <- opt$duration %||% list()

  # Global baseline-rate anchor table: run-invariant, so none of the per-run
  # tokens below apply. No pop digest -- the duration denominator is exposure
  # days alone. See `.base_rate_path()`.
  if (kind == "dur_base_rate") {
    return(.base_rate_path(
      opt,
      stem = "base_rate_table_dur",
      estimator = dur$base_rate$estimator %||% "pooled",
      pre_tokens = paste0("base-", opt$baseline %||% "1990_2019"),
      pop_token = FALSE,
      suffix = paste0(
        .base_rate_flag_token(opt, "duration"),
        .episode_key_path_token(opt)
      )
    ))
  }

  # Country-grouping mode: stable group token instead of the comma-joined
  # country list (mirrors make_size_run_file / make_run_file).
  if (isTRUE(opt$grouping$enable)) {
    loc_token <- sprintf("grp-%s-%s", opt$grouping$source, opt$grouping$group)
  } else {
    loc_token <- paste0("ctry-", opt$country %||% "all")
  }

  parts <- c(
    "outbreak_duration",
    dur$backend %||% "brms",
    sub("^dur_", "", kind),
    paste0("base-", opt$baseline %||% "1990_2019"),
    loc_token,
    paste0("bhd-", dur$basehaz_duration %||% "factor")
  )

  bhp <- dur$basehaz_poisson %||% "none"
  if (!identical(bhp, "none")) {
    parts <- c(parts, paste0("bhp-", bhp))
  }
  # Calendar-month seasonality: a token only when on (the monthly track's
  # "seas-" token), so every existing name is unchanged.
  dseas <- dur[["seasonal_term"]] %||% "none"
  if (!identical(dseas, "none")) {
    parts <- c(parts, paste0("seas-", dseas))
  }
  # Own prior on the band effects (brms block prior_mso_band): "mbp<sd>" for
  # normal(0, <sd>), e.g. "mbp1p5"; unset, no token.
  mbp <- (dur[[dur$backend %||% "brms"]] %||% dur[["brms"]] %||% list())[[
    "prior_mso_band"
  ]]
  if (!is.null(mbp)) {
    sd_mbp <- sub(
      "^\\s*normal\\(\\s*0\\s*,\\s*([0-9.]+)\\s*\\)\\s*$",
      "\\1",
      mbp
    )
    if (identical(sd_mbp, mbp)) {
      stop(
        "make_duration_run_file: prior_mso_band must be 'normal(0, <sd>)' to be named; got '",
        mbp,
        "'.",
        call. = FALSE
      )
    }
    parts <- c(parts, paste0("mbp", gsub(".", "p", sd_mbp, fixed = TRUE)))
  }

  # Only non-default break sets earn a token, so every pre-existing artefact
  # name stays stable if the default edges are ever revisited via config.
  brk <- as.numeric(unlist(dur$mso_breaks %||% c(1, 2, 3, 4, 6, 9, Inf)))
  if (!isTRUE(all.equal(brk, c(1, 2, 3, 4, 6, 9, Inf)))) {
    parts <- c(
      parts,
      paste0("mso-", paste(brk[is.finite(brk)], collapse = "-"))
    )
  }

  parts <- c(parts, paste0("adm-", dur$admin_level %||% "ADM1"))
  # ADM2 intercept / leaf tokens, as on the monthly track; no duration config
  # before ADM1-ADM2 sets either key, so every existing name is unchanged.
  if (isTRUE(dur$random_effects$re_adm2_intercept)) {
    parts <- c(parts, "adm2i")
  }
  parts <- c(parts, .leaf_token(dur))
  if (isTRUE(dur$covariates$enso)) {
    parts <- c(parts, paste0("enso-lag", opt$enso$lag_months))
  }

  # DLNM crossbasis single-lag token — mirrors make_size_run_file: disambiguates
  # per-lag onebasis fits in a sweep (the distributed-lag cb fit carries no
  # lag token; single-lag must, else all lags overwrite the same file).
  dlnm_cfg <- dur$dlnm %||% list()
  if (
    identical(dur$backend %||% "brms", "brms_cb") && isTRUE(dlnm_cfg$enable)
  ) {
    # Exposure-variable token ("" for ENSO: every existing name unchanged).
    vt <- .dlnm_vars_token(dlnm_cfg)
    if (nzchar(vt)) {
      parts <- c(parts, sub("^_", "", vt))
    }
    sl <- .dlnm_single_lag(dlnm_cfg)
    if (!is.null(sl)) {
      parts <- c(parts, sl_tag(sl))
    }
    bt <- .dlnm_basis_token(dlnm_cfg) # "" for bs (back-compat), "ns4" etc. for ns
    if (nzchar(bt)) parts <- c(parts, bt)
  }

  assert_basename_length(
    file.path(
      opt$outdir %||% "generated_data/",
      paste0(paste(parts, collapse = "_"), .episode_key_path_token(opt), ".rds")
    ),
    "make_duration_run_file"
  )
}

#' Monthly-size track artefact paths (fourth model track: per-calendar-month
#' size, response `cases_month`). Structural copy of make_duration_run_file(),
#' plus the estimand token: monthly carries the same denominator
#' keys as size (`monthly.estimand`, `monthly.pop_elasticity`,
#' `monthly.canonical_admin_level` -- see `.ESTIMAND_SPEC` in
#' utils/provenance.R), and those change both the fitted model and every cell
#' of the `base_rate_mon` anchor. The token goes on BOTH the per-run artefacts
#' and the anchor sidecar, because each fit script early-exits on
#' `file.exists()`: a footprint-denominated and a gid-denominated mon fit
#' sharing one filename would not overwrite and would not error -- the second
#' run would report success over the first run's fit.
#' The size-only `nolf`/`noxal` flags are NOT double-encoded here; monthly's
#' equivalents come from `.base_rate_flag_token(opt, "monthly")` (see
#' `.estimand_path_token()`).
make_monthly_run_file <- function(
  opt,
  kind = c(
    "mon_fit",
    "mon_re_stats",
    "mon_coef_summary",
    "mon_loo",
    "mon_kfold",
    "mon_pp_check",
    "mon_cb_meta",
    "mon_cb_surfaces",
    "mon_base_rate"
  )
) {
  kind <- match.arg(kind)
  mon <- opt$monthly %||% list()

  # Global baseline-rate anchor table: run-invariant, so none of the per-run
  # tokens below apply. UNLIKE occ/dur this one DOES carry the population
  # digest: the monthly anchor's denominator is `pop * exposure_days`, so a
  # population change silently rewrites every cell -- exactly the collision
  # `current_pop_digest_token()` was introduced to close for size.
  if (kind == "mon_base_rate") {
    return(.base_rate_path(
      opt,
      stem = "base_rate_table_mon",
      estimator = mon$base_rate$estimator %||% "pooled",
      pre_tokens = paste0("base-", opt$baseline %||% "1990_2019"),
      pop_token = TRUE,
      suffix = paste0(
        .base_rate_flag_token(opt, "monthly"),
        .estimand_path_token(opt, "monthly"),
        .episode_key_path_token(opt)
      )
    ))
  }

  # Country-grouping mode: stable group token instead of the comma-joined
  # country list (mirrors make_size_run_file / make_duration_run_file).
  if (isTRUE(opt$grouping$enable)) {
    loc_token <- sprintf("grp-%s-%s", opt$grouping$source, opt$grouping$group)
  } else {
    loc_token <- paste0("ctry-", opt$country %||% "all")
  }

  parts <- c(
    "outbreak_monthly",
    mon$backend %||% "brms_cb",
    sub("^mon_", "", kind),
    paste0("base-", opt$baseline %||% "1990_2019"),
    loc_token,
    paste0("fam-", mon$family %||% "negbinomial")
  )

  if (!identical(mon$seasonal_term %||% "none", "none")) {
    parts <- c(parts, paste0("seas-", mon$seasonal_term))
  }

  parts <- c(parts, paste0("adm-", mon$admin_level %||% "ADM1-ADM2"))
  if (isTRUE(mon$random_effects$re_adm2_intercept)) {
    parts <- c(parts, "adm2i")
  }

  # Two configs differing only in prev_cases would otherwise collide on the
  # same filename (same rationale as duration's mso- token).
  if (isTRUE(mon$prev_cases)) {
    parts <- c(parts, "prevcases")
  }

  if (isTRUE(mon$covariates$enso)) {
    parts <- c(parts, paste0("enso-lag", opt$enso$lag_months))
  }

  # DLNM crossbasis single-lag token — mirrors make_size_run_file/make_duration_run_file:
  # disambiguates per-lag onebasis fits in a sweep (the distributed-lag cb fit
  # carries no lag token; single-lag must, else all lags overwrite one file).
  dlnm_cfg <- mon$dlnm %||% list()
  if (
    identical(mon$backend %||% "brms_cb", "brms_cb") && isTRUE(dlnm_cfg$enable)
  ) {
    # Exposure-variable token ("" for ENSO: every existing name unchanged).
    vt <- .dlnm_vars_token(dlnm_cfg)
    if (nzchar(vt)) {
      parts <- c(parts, sub("^_", "", vt))
    }
    sl <- .dlnm_single_lag(dlnm_cfg)
    if (!is.null(sl)) {
      parts <- c(parts, sl_tag(sl))
    }
    bt <- .dlnm_basis_token(dlnm_cfg) # "" for bs (back-compat), "ns4" etc. for ns
    if (nzchar(bt)) parts <- c(parts, bt)
  }

  phi_re <- mon$phi_re %||% "constant"
  if (phi_re != "constant") {
    parts <- c(parts, paste0("phire-", gsub("_", "-", phi_re)))
  }
  parts <- c(parts, .leaf_token(mon))

  assert_basename_length(
    file.path(
      opt$outdir %||% "generated_data/",
      paste0(
        paste(parts, collapse = "_"),
        .estimand_path_token(opt, "monthly"),
        .episode_key_path_token(opt),
        ".rds"
      )
    ),
    "make_monthly_run_file"
  )
}

#' Derive a simulator run_tag from the config filename.
#'
#' @param ts_prefix Prefix used for the timestamped fallback tag (default "sim").
#' @return Character run_tag: the config's basename sans extension, or a
#'   timestamped tag when no `-C`/`--config` was passed.
#' @details Shared by the two simulator offspring (`07b_simulator_v1.R` uses
#'   `ts_prefix` "sim"; `07b_simulator_spatpomp.R` uses "spatpomp"). Depends on
#'   `.flag_val()` (utils/config.R).
derive_sim_run_tag <- function(ts_prefix = "sim") {
  cfg_path <- .flag_val("-C") %||% .flag_val("--config")
  if (!is.null(cfg_path)) {
    tools::file_path_sans_ext(basename(cfg_path))
  } else {
    format(Sys.time(), paste0(ts_prefix, "_%Y%m%d_%H%M%S"))
  }
}

#' Derive a config/batch run tag for the figure scripts.
#'
#' @param config_path Explicit config (or batch directory) path; falls back to
#'   `-C`/`--config` on the command line (`.flag_val()`, utils/config.R).
#' @param default Fallback tag when neither `config_path` nor a CLI flag is
#'   available (default "run").
#' @return Character run tag: the basename of `config_path`, sans extension.
#' @details Sibling of `derive_sim_run_tag()`. Unlike `derive_sim_run_tag()`
#'   there is no timestamped fallback: figure directories are named after the
#'   config that produced the fit being plotted, not the plotting invocation
#'   itself. `tools::file_path_sans_ext(basename(...))` is a no-op on an
#'   extension-less directory basename, so passing a batch *directory* (as
#'   `06_fit_figures_lag_compare.R`'s `--batch-dir` already does) reproduces
#'   that script's existing `basename(normalizePath(...))` behaviour
#'   unchanged.
derive_config_run_tag <- function(config_path = NULL, default = "run") {
  path <- config_path %||% .flag_val("-C") %||% .flag_val("--config")
  if (is.null(path)) {
    return(default)
  }
  path <- normalizePath(path, mustWork = FALSE)
  tools::file_path_sans_ext(basename(path))
}

#' Figure output directory builder.
#'
#' @param opt Resolved options list; only `opt$figdir` is used (default "figures").
#' @param track Track segment (e.g. "cb"), or NULL to omit it.
#' @param config_path Passed through to `derive_config_run_tag()`.
#' @param subdir Optional trailing subdirectory (e.g. "adm1").
#' @return Character directory path. Path algebra only -- does not create the
#'   directory (callers / `figures_io.R::write_figure_files()` do that once, at write
#'   time).
#' @details Resolves to `<opt$figdir %||% "figures/">/[<track>/]<batch>[/<subdir>]`,
#'   where `<batch> = derive_config_run_tag(config_path)` (the config's own
#'   basename, sans extension). `track = "cb"` + no subdir reproduces
#'   `figures/cb/<run_tag>/` exactly, matching `06_fit_figures_cb.R`'s current
#'   hand-rolled `fig_dir`; `track = "cb"` + `subdir = "adm1"` matches
#'   `06_fit_figures_cb_adm1.R`'s `fig_dir`. `track = NULL` (`06_fit_figures.R`'s
#'   own figures/ root, disambiguated by filename tag rather than a directory)
#'   is supported for callers that want the bare figdir with no track segment.
make_figure_dir <- function(
  opt = list(),
  track = NULL,
  config_path = NULL,
  subdir = NULL
) {
  figdir <- sub("/+$", "", opt$figdir %||% "figures")
  tag <- derive_config_run_tag(config_path)
  parts <- c(figdir, track, tag, subdir)
  assert_basename_length(
    do.call(file.path, as.list(parts[nzchar(parts)])),
    "make_figure_dir"
  )
}

#' Results-postprocess artefact paths (Phase 09 -> Phase 10).
#'
#' @param opt Resolved options list; uses `opt$outdir` only.
#' @param kind One of "results" (the run directory), "results_tasks_dir"
#'   (its 09a task-YAML dir), "results_enso_relation",
#'   "results_track_params", "results_onset_rate" (unkeyed),
#'   "results_attribution", "results_smoothing", "results_oni", or one of the
#'   phase-11 "key numbers" report bundles: "results_report_dataset" (unkeyed),
#'   "results_report_effects" (track-keyed), "results_report_smoothing"
#'   (unkeyed), "results_report_forecast" (sets-keyed),
#'   "results_report_enso_exposure" (unkeyed), "results_report_attrib_exposure"
#'   (sets-keyed), the hindcast and
#'   supplementary-table kinds listed below, and "results_rain_enso"
#'   (Figure 1d's rainfall tables, unkeyed).
#' @param track Modelling track for `results_enso_relation` /
#'   `results_report_effects` ("occurrence" | "size" | "duration" | "monthly").
#' @param sets Character vector of ENSO analog sets, for the set-keyed
#'   kinds. Order-insensitive: encoded via `sets_key()`'s sorted-and-hyphenated
#'   form (`utils/results_core.R`), reproduced here so `paths.R` stays free of
#'   cross-module dependencies.
#' @param schema Schema-version token for this artefact
#'   (`RESULTS_SCHEMA[[kind]]`, utils/results_core.R). Required for every file
#'   kind: it is the ONLY invalidation channel for a change of estimand or an
#'   added column, neither of which a modification-time comparison can see.
#' @param config_path Passed through to `derive_config_run_tag()`.
#' @return Character path. Path algebra only -- creates nothing.
#' @details Resolves to `<outdir>/results/<tag>/<file>`, where
#'   `<tag> = derive_config_run_tag(config_path)` -- the same tag
#'   `make_figure_dir()` uses, so `analysis/configs/results/manuscript.yml`
#'   produces `generated_data/results/manuscript/` alongside
#'   `figures/manuscript/`.
#'
#'   Deliberately NOT routed through `make_run_file()`: that dispatcher
#'   (`make_run_file()`, below) forwards only `kind` to every branch and drops
#'   `...`, so `track`/`sets`/`schema` could not survive the hop. Callers
#'   invoke this builder directly; nothing here needs generic dispatch.
make_results_run_file <- function(
  opt,
  kind = c(
    "results",
    "results_tasks_dir",
    "results_enso_relation",
    "results_track_params",
    "results_onset_rate",
    "results_attribution",
    "results_smoothing",
    "results_oni",
    "results_report_dataset",
    "results_report_effects",
    "results_report_smoothing",
    "results_report_forecast",
    "results_report_enso_exposure",
    "results_report_attrib_exposure",
    "results_report_attrib_metrics",
    "results_report_event_attrib",
    "results_report_event_af_model",
    "results_af_posterior",
    "results_event_af_posterior",
    "results_attrib_exposure_posterior",
    "results_report_attrib_information",
    "results_report_attrib_incidence",
    "results_hindcast_runs",
    "results_hindcast_occ",
    "results_hindcast_cases",
    "results_report_hindcast",
    "results_convergence",
    "results_report_convergence",
    "results_report_model_spec",
    "results_report_lags",
    "results_report_anchor",
    "results_report_software",
    "results_report_kl_sweep",
    "results_rain_enso",
    "results_group_labels",
    "results_lag_averaging",
    "results_report_membership",
    "results_report_who_scale",
    "results_report_analog_sets",
    "results_report_national_zero_audit",
    "results_report_occ_panel_exclusions",
    "results_report_contagion_coverage",
    "results_modelled_units",
    "results_occ_seasonality",
    "results_seasonality",
    "results_oni_composite",
    "results_report_who_enso_ratio",
    "results_report_confirmed_agreement"
  ),
  track = NULL,
  sets = NULL,
  schema = NULL,
  config_path = NULL
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_dir <- assert_basename_length(
    file.path(outdir, "results", derive_config_run_tag(config_path)),
    "make_results_run_file"
  )
  if (kind == "results") {
    return(run_dir)
  }
  # 09a_postprocess_results_tasks.R's per-step task YAMLs (one SLURM array
  # task each; analysis/bash/stages/09_postprocess_results.sh).
  if (kind == "results_tasks_dir") {
    return(file.path(run_dir, "tasks"))
  }

  if (is.null(schema) || !nzchar(schema)) {
    stop(
      "make_results_run_file: `schema` is required for kind '",
      kind,
      "'. Pass RESULTS_SCHEMA[[<stem>]] (utils/results_core.R) -- an ",
      "unversioned filename cannot be invalidated when the derived ",
      "table's estimand or columns change.",
      call. = FALSE
    )
  }

  key <- function() {
    if (is.null(sets) || !length(sets)) {
      stop(
        "make_results_run_file: `sets` is required for kind '",
        kind,
        "'.",
        call. = FALSE
      )
    }
    paste(sort(sets), collapse = "-")
  }

  assert_basename_length(
    file.path(
      run_dir,
      switch(
        kind,
        results_enso_relation = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind 'results_enso_relation'.",
              call. = FALSE
            )
          }
          sprintf("enso_relation_%s_track-%s.qs", schema, track)
        },
        results_track_params = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind 'results_track_params'.",
              call. = FALSE
            )
          }
          sprintf("track_params_%s_track-%s.qs", schema, track)
        },
        # Figure 1b inset: one unkeyed per-gid table.
        results_onset_rate = sprintf("onset_rate_%s.qs", schema),
        results_attribution = sprintf(
          "attribution_%s_metric-cases-mon_sets-%s.qs",
          schema,
          key()
        ),
        results_smoothing = sprintf("smoothing_%s_metric-cases-mon.qs", schema),
        results_oni = sprintf("oni_envelope_%s_sets-%s.qs", schema, key()),
        # --- Manuscript "key numbers" report artefacts (phase 11) -------------
        # dataset/smoothing are single (unkeyed) bundles; effects is track-keyed
        # (one artefact per track); forecast is sets-keyed like attribution/oni.
        results_report_dataset = sprintf("report_dataset_%s.qs", schema),
        results_report_effects = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind ",
              "'results_report_effects'.",
              call. = FALSE
            )
          }
          sprintf("report_effects_%s_track-%s.qs", schema, track)
        },
        results_report_smoothing = sprintf("report_smoothing_%s.qs", schema),
        results_report_forecast = sprintf(
          "report_forecast_%s_sets-%s.qs",
          schema,
          key()
        ),
        # Cross-track (the track set is in the stamped scope.fig1_tracks): unkeyed.
        results_report_enso_exposure = sprintf(
          "report_enso_exposure_%s.qs",
          schema
        ),
        # Forecast vs counterfactual (same runs and window as report_forecast): sets-keyed.
        results_report_attrib_exposure = sprintf(
          "report_attrib_exposure_%s_sets-%s.qs",
          schema,
          key()
        ),
        # Attribution metrics of the same paired runs: sets-keyed.
        results_report_attrib_metrics = sprintf(
          "report_attrib_metrics_%s_sets-%s.qs",
          schema,
          key()
        ),
        # Past-ENSO-event attribution: one bundle (events and set from the event manifest).
        results_report_event_attrib = sprintf(
          "report_event_attrib_%s.qs",
          schema
        ),
        results_report_event_af_model = sprintf(
          "report_event_af_model_%s.qs",
          schema
        ),
        # Posterior AF / AN: forecast sets-keyed (same runs as report_attrib_metrics),
        # events unkeyed (set from the event manifest).
        results_af_posterior = sprintf(
          "af_posterior_%s_sets-%s.qs",
          schema,
          key()
        ),
        results_event_af_posterior = sprintf(
          "event_af_posterior_%s.qs",
          schema
        ),
        # Posterior Childs population: sets-keyed (same runs as af_posterior).
        results_attrib_exposure_posterior = sprintf(
          "attrib_exposure_posterior_%s_sets-%s.qs",
          schema,
          key()
        ),
        # Information component of the forecast attribution: sets-keyed (same runs).
        results_report_attrib_information = sprintf(
          "report_attrib_information_%s_sets-%s.qs",
          schema,
          key()
        ),
        # Window totals / incidence by Africa, group and country of the same runs: sets-keyed.
        results_report_attrib_incidence = sprintf(
          "report_attrib_incidence_%s_sets-%s.qs",
          schema,
          key()
        ),
        # --- Hindcast skill (unkeyed: one bundle per kind; every scoring choice
        # is in the stamped scope/derivation.hindcast blocks) -----------------
        results_hindcast_runs = sprintf("hindcast_runs_%s.qs", schema),
        results_hindcast_occ = sprintf("hindcast_occ_%s.qs", schema),
        results_hindcast_cases = sprintf("hindcast_cases_%s.qs", schema),
        results_report_hindcast = sprintf("report_hindcast_%s.qs", schema),
        # --- Supplementary tables ----------------------------------------------
        # convergence is track-keyed (one artefact per track's deployed fits);
        # the reports are cross-track bundles.
        results_convergence = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind 'results_convergence'.",
              call. = FALSE
            )
          }
          sprintf("convergence_%s_track-%s.qs", schema, track)
        },
        results_report_convergence = sprintf(
          "report_convergence_%s.qs",
          schema
        ),
        # Fit-level model specification, track-keyed like convergence.
        results_report_model_spec = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind 'results_report_model_spec'.",
              call. = FALSE
            )
          }
          sprintf("report_model_spec_%s_track-%s.qs", schema, track)
        },
        results_report_lags = sprintf("report_lags_%s.qs", schema),
        results_report_anchor = sprintf("report_anchor_%s.qs", schema),
        results_report_software = sprintf("report_software_%s.qs", schema),
        results_report_kl_sweep = sprintf("report_kl_sweep_%s.qs", schema),
        # Figure 1d's rainfall tables: unkeyed (the grouping is the stamped
        # inputs.grouping, the method results.derivation.rain_enso).
        results_rain_enso = sprintf("rain_enso_%s.qs", schema),
        # Group display lags + per-country season ends (S30): unkeyed.
        results_group_labels = sprintf("group_labels_%s.qs", schema),
        # Lag averaging over the single-lag fits: track-keyed.
        results_lag_averaging = {
          if (is.null(track) || !nzchar(track)) {
            stop(
              "make_results_run_file: `track` is required for kind 'results_lag_averaging'.",
              call. = FALSE
            )
          }
          sprintf("lag_averaging_%s_track-%s.qs", schema, track)
        },
        # Membership uncertainty of the deployed classes: unkeyed.
        results_report_membership = sprintf("report_membership_%s.qs", schema),
        # Per-country WHO term scale of the smoothing: unkeyed.
        results_report_who_scale = sprintf("report_who_scale_%s.qs", schema),
        # Counterfactual analog sets: unkeyed.
        results_report_analog_sets = sprintf(
          "report_analog_sets_%s.qs",
          schema
        ),
        # Occurrence-panel tables: unkeyed.
        results_report_national_zero_audit = sprintf(
          "report_national_zero_audit_%s.qs",
          schema
        ),
        results_report_occ_panel_exclusions = sprintf(
          "report_occ_panel_exclusions_%s.qs",
          schema
        ),
        results_report_contagion_coverage = sprintf(
          "report_contagion_coverage_%s.qs",
          schema
        ),
        # Supplement maps MU / GM / KC: modelled units + fit membership, unkeyed.
        results_modelled_units = sprintf("modelled_units_%s.qs", schema),
        # Figure 1c's onset-season column: occurrence seasonal spline, unkeyed.
        results_occ_seasonality = sprintf("occ_seasonality_%s.qs", schema),
        # Supplement fig. SE: rainfall vs the three tracks by month, unkeyed.
        results_seasonality = sprintf("seasonality_%s.qs", schema),
        # Figure 1c: the composite ONI cycle, unkeyed.
        results_oni_composite = sprintf("oni_composite_%s.qs", schema),
        # Supplement fig. WR: observed vs modelled ENSO ratio of 2026, unkeyed.
        results_report_who_enso_ratio = sprintf(
          "report_who_enso_ratio_%s.qs",
          schema
        ),
        # Supplementary table: ENSO contrasts, all vs confirmed outbreaks, unkeyed.
        results_report_confirmed_agreement = sprintf(
          "report_confirmed_agreement_%s.qs",
          schema
        )
      )
    ),
    "make_results_run_file"
  )
}

#' Deploy-tool (utils/deploy_prod_run.R) output paths.
#'
#' Mirror of make_results_run_file above for the prod_run<N> deployment orchestrator
#' (analysis/utils/deploy_prod_run.R): the run tag is the deploy config's own
#' basename-sans-extension (e.g. `prod_run`), so artefacts land under
#' `<outdir>/deploy/<stem>/`. "deploy" returns the run-root directory;
#' "deploy_manifest" returns the provenance manifest path inside it.
#'
#' @param opt Resolved options list; only `opt$outdir` is used.
#' @param kind "deploy" (run root) | "deploy_manifest".
#' @param config_path Path to the deploy YAML (for the run-tag stem). Falls back
#'   to the active -C/--config flag via derive_config_run_tag().
#' @return Character path (path algebra only; does not create the directory).
make_deploy_run_file <- function(
  opt,
  kind = c("deploy", "deploy_manifest"),
  config_path = NULL
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_dir <- assert_basename_length(
    file.path(outdir, "deploy", derive_config_run_tag(config_path)),
    "make_deploy_run_file"
  )
  switch(
    kind,
    deploy = run_dir,
    deploy_manifest = file.path(run_dir, "deploy_manifest.qs")
  )
}

#' LOO lag-compare / batch coef-compare output paths (phase-05 derived,
#' per-batch artefacts written by `05_fit_postprocess_loo_lags.R` and the
#' batch-mode branch of `05_fit_postprocess.R`).
#'
#' Mirror of make_deploy_run_file above: the run tag is the batch name, so both
#' a batch config DIRECTORY (`analysis/configs/<batch>/`, what every `-d`
#' caller passes) and a bare batch name (`"<batch>"`, what `deploy_core.R`
#' holds in `dep$fits$*_dir`) resolve to `<outdir>/loo_compare/<batch>/`.
#' `utils/loo_lags.R::loo_compare_ranking_files()` also reads the legacy
#' location inside the batch config dir
#' (`analysis/configs/<batch>/loo_compare_lags_*.rds`) as a fallback; this
#' builder only ever names the new one.
#'
#' @param opt Resolved options list; only `opt$outdir` is used.
#' @param kind "loo_compare_dir" (batch root) | "loo_compare_ranking"
#'   (`loo_compare_lags[_<group>]_<YYYYMMDD>.rds`) | "loo_compare_png"
#'   (`loo_compare_lags[_<group>].png`) | "loo_compare_coef"
#'   (`coef_compare_<YYYYMMDD>.rds`) | "loo_compare_lag_sidecar"
#'   (`lag_loo_sidecar[_<group>]_<YYYYMMDD>.rds`: pointwise elpd, Pareto-k and
#'   per-lag crossbasis draws, `05_fit_postprocess_loo_lags.R --lag_sidecar`;
#'   named so it never matches the ranking regex, utils/loo_lags.R).
#' @param batch Batch config dir path or bare batch name (see above).
#' @param group Optional group token (e.g. "enso_east"); NULL/"" = ungrouped.
#' @param date Date stamped into the dated kinds [default today].
#' @return Character path (path algebra only; does not create the directory).
make_loo_compare_run_file <- function(
  opt = list(),
  kind = c(
    "loo_compare_dir",
    "loo_compare_ranking",
    "loo_compare_png",
    "loo_compare_coef",
    "loo_compare_lag_sidecar"
  ),
  batch,
  group = NULL,
  date = Sys.Date()
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_dir <- file.path(outdir, "loo_compare", derive_config_run_tag(batch))
  grp_sfx <- if (!is.null(group) && nzchar(group)) paste0("_", group) else ""
  d <- format(as.Date(date), "%Y%m%d")
  assert_basename_length(
    switch(
      kind,
      loo_compare_dir = run_dir,
      loo_compare_ranking = file.path(
        run_dir,
        sprintf("loo_compare_lags%s_%s.rds", grp_sfx, d)
      ),
      loo_compare_png = file.path(
        run_dir,
        sprintf("loo_compare_lags%s.png", grp_sfx)
      ),
      loo_compare_coef = file.path(run_dir, sprintf("coef_compare_%s.rds", d)),
      loo_compare_lag_sidecar = file.path(
        run_dir,
        sprintf("lag_loo_sidecar%s_%s.rds", grp_sfx, d)
      )
    ),
    "make_loo_compare_run_file"
  )
}

#' brms_kl class-sweep comparison paths (05_fit_postprocess_kl_sweep.R ->
#' 06_fit_figures_kl_sweep.R and the 09 report_kl_sweep step).
#'
#' Mirror of make_loo_compare_run_file(): one directory per sweep under
#' `<outdir>/kl_sweep_compare/<tag>/`, `<tag> = derive_config_run_tag(batch)`
#' (a batch config dir path, or a bare tag when several dirs are compared).
#'
#' @param opt Resolved options list; only `opt$outdir` is used.
#' @param kind "kl_sweep_dir" or one of the tables: kl_sweep_fits,
#'   kl_sweep_classes, kl_sweep_partition, kl_sweep_loo, kl_sweep_loo_country,
#'   kl_sweep_diagnostics, kl_sweep_diagnostics_chains,
#'   kl_sweep_responsibilities, kl_sweep_country_mixture, kl_sweep_stability
#'   (.csv), kl_sweep_rds (loo.rds: inventory, agreeing chains, loo / waic
#'   objects, settings).
#' @param batch Batch config dir path or bare tag.
#' @return Character path (path algebra only).
make_kl_sweep_run_file <- function(
  opt = list(),
  kind = c(
    "kl_sweep_dir",
    "kl_sweep_fits",
    "kl_sweep_classes",
    "kl_sweep_partition",
    "kl_sweep_loo",
    "kl_sweep_loo_country",
    "kl_sweep_diagnostics",
    "kl_sweep_diagnostics_chains",
    "kl_sweep_responsibilities",
    "kl_sweep_country_mixture",
    "kl_sweep_stability",
    "kl_sweep_rds"
  ),
  batch
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_dir <- assert_basename_length(
    file.path(outdir, "kl_sweep_compare", derive_config_run_tag(batch)),
    "make_kl_sweep_run_file"
  )
  if (kind == "kl_sweep_dir") {
    return(run_dir)
  }
  if (kind == "kl_sweep_rds") {
    return(file.path(run_dir, "loo.rds"))
  }
  file.path(run_dir, paste0(sub("^kl_sweep_", "", kind), ".csv"))
}

# Shared helper for the simulator path builders below (make_sim_run_file,
# make_smooth_run_file, make_spatpomp_run_file). Resolves a sub-directory under a run
# root, optionally nested under a run_id: run_id NULL/empty -> flat
# file.path(run_dir, sub); otherwise file.path(run_dir, sub, run_id). Only v1
# (make_sim_run_file) passes a run_id; the two spatPomp builders stay flat.
.sim_scoped <- function(run_dir, sub, run_id = NULL) {
  if (!is.null(run_id) && nzchar(run_id)) {
    file.path(run_dir, sub, run_id)
  } else {
    file.path(run_dir, sub)
  }
}

#' Simulator (v1) output paths.
#'
#' @param opt Resolved options list; uses `opt$outdir`, `opt$simulator$run_tag`,
#'   `opt$simulator$run_id`, `opt$simulator$output$format`.
#' @param kind One of "sim", "sim_tasks_dir", "sim_parts_dir", "sim_merged_dir",
#'   "sim_valid_dir", "sim_event_table", "sim_monthly_summary".
#' @return Character path. "sim" returns the run root directory; the other
#'   kinds return sub-paths for the materialised pieces.
#' @details v1 uses `opt$simulator$run_tag` (or a fallback "default") under
#'   `<outdir>/sim/`.
make_sim_run_file <- function(
  opt,
  kind = c(
    "sim",
    "sim_tasks_dir",
    "sim_parts_dir",
    "sim_merged_dir",
    "sim_valid_dir",
    "sim_event_table",
    "sim_monthly_summary"
  )
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_tag <- opt$simulator$run_tag %||% "default"
  run_dir <- file.path(outdir, "sim", run_tag)
  fmt_ext <- switch(
    opt$simulator$output$format %||% "qs",
    parquet = ".parquet",
    rds = ".rds",
    qs = ".qs",
    ".qs"
  )
  # Optional run_id (e.g. a timestamp set in the config) isolates partition
  # directories per invocation — prevents stale partitions from a previous
  # scope/parameter change from polluting aggregation. The merged and
  # validation outputs are scoped by the SAME run_id so a bumped run_id never
  # silently reads stale merged/validation results from a previous invocation.
  # When absent, falls back to the flat layout for backward compat.
  run_id <- opt$simulator$run_id
  scoped <- function(sub) .sim_scoped(run_dir, sub, run_id)
  assert_basename_length(
    switch(
      kind,
      sim = run_dir,
      sim_tasks_dir = file.path(run_dir, "tasks"),
      sim_parts_dir = scoped("parts"),
      sim_merged_dir = scoped("merged"),
      sim_valid_dir = scoped("validation"),
      sim_event_table = file.path(run_dir, paste0("event_table", fmt_ext)),
      sim_monthly_summary = file.path(
        run_dir,
        paste0("monthly_summary", fmt_ext)
      )
    ),
    "make_sim_run_file"
  )
}

#' Path builder for the spatPomp historical-smoothing chunked pipeline.
#'
#' @param opt Resolved options list; uses `opt$outdir`, `opt$simulator$run_tag`.
#' @param kind One of "smooth", "smooth_tasks_dir", "smooth_parts_dir",
#'   "smooth_merged_dir", "smooth_cache_dir".
#' @return Character path.
#'   \describe{
#'     \item{smooth}{run_dir}
#'     \item{smooth_tasks_dir}{per-(TU, rep-chunk) task YAMLs (tasks/)}
#'     \item{smooth_parts_dir}{per-task partition outputs (parts/)}
#'     \item{smooth_merged_dir}{compiled ribbons + diagnostics (merged/)}
#'     \item{smooth_cache_dir}{per-TU handoff trajectory caches (cache/) --
#'       the forecast runner's `smoothing_input` source}
#'   }
#' @details Pipeline: `08a_smoothing_tasks.R` -> `08b_smoothing_spatpomp.R` ->
#'   `08c_smoothing_compile.R`. Mirrors `make_sim_run_file()` but roots everything
#'   under `<outdir>/sim_spatpomp/<run_tag>/`.
make_smooth_run_file <- function(
  opt,
  kind = c(
    "smooth",
    "smooth_tasks_dir",
    "smooth_parts_dir",
    "smooth_merged_dir",
    "smooth_cache_dir"
  )
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_tag <- opt$simulator$run_tag %||% "s10_historical_smoothing"
  run_dir <- file.path(outdir, "sim_spatpomp", run_tag)
  # Flat layout (run_id = NULL).
  assert_basename_length(
    switch(
      kind,
      smooth = run_dir,
      smooth_tasks_dir = .sim_scoped(run_dir, "tasks"),
      smooth_parts_dir = .sim_scoped(run_dir, "parts"),
      smooth_merged_dir = .sim_scoped(run_dir, "merged"),
      smooth_cache_dir = .sim_scoped(run_dir, "cache")
    ),
    "make_smooth_run_file"
  )
}

#' Path builder for the spatPomp forecast/counterfactual runner.
#'
#' @param opt Resolved options list; uses `opt$outdir`, `opt$simulator$run_tag`.
#' @param kind One of "spatpomp", "spatpomp_parts_dir", "spatpomp_diag_dir",
#'   "spatpomp_cache_dir", "spatpomp_merged_dir", "spatpomp_valid_dir".
#' @return Character path.
#'   \describe{
#'     \item{spatpomp}{run_dir}
#'     \item{spatpomp_parts_dir}{per-(TU, [repoff-chunk]) partition outputs (parts/)}
#'     \item{spatpomp_diag_dir}{bpfilter ESS/diagnostics (filter_diagnostics/)}
#'     \item{spatpomp_cache_dir}{per-TU/replicate warmup+forecast caches (cache/)}
#'     \item{spatpomp_merged_dir}{07b --aggregate output (merged/)}
#'     \item{spatpomp_valid_dir}{07b --validate output (validation/)}
#'   }
#' @details For `07b_simulator_spatpomp.R`. Mirrors that script's own (manual,
#'   hardcoded) path construction EXACTLY -- same root
#'   `<outdir>/sim_spatpomp/<run_tag>/`, FLAT (no run_id-scoping of
#'   parts/merged/validation, unlike `make_sim_run_file()`'s v1 layout). This is
#'   deliberate, not an oversight: both pre-existing real spatPomp forecast
#'   runs (`sim_africa_brms_elag5_spatpomp_forecast`,
#'   `sim_africa_brms_elag9_spatpomp_forecast`) set `run_id` in their configs,
#'   but their actual on-disk `parts/` layout is flat -- retrofitting
#'   run_id-scoping here would silently orphan that existing data. See the
#'   cholera-spatpomp-forecast-runid-scoping-gap fleeting note for the open
#'   question of whether 08 should later be changed to honor run_id (in which
#'   case this helper's `scoped()` would need to match, analogous to
#'   `make_sim_run_file()` above).
make_spatpomp_run_file <- function(
  opt,
  kind = c(
    "spatpomp",
    "spatpomp_parts_dir",
    "spatpomp_diag_dir",
    "spatpomp_cache_dir",
    "spatpomp_merged_dir",
    "spatpomp_valid_dir"
  )
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  run_tag <- opt$simulator$run_tag %||% "default"
  run_dir <- file.path(outdir, "sim_spatpomp", run_tag)
  # Flat layout (run_id = NULL) — deliberate; see header note above.
  assert_basename_length(
    switch(
      kind,
      spatpomp = run_dir,
      spatpomp_parts_dir = .sim_scoped(run_dir, "parts"),
      spatpomp_diag_dir = .sim_scoped(run_dir, "filter_diagnostics"),
      spatpomp_cache_dir = .sim_scoped(run_dir, "cache"),
      spatpomp_merged_dir = .sim_scoped(run_dir, "merged"),
      spatpomp_valid_dir = .sim_scoped(run_dir, "validation")
    ),
    "make_spatpomp_run_file"
  )
}

#' make_sim_engine_dirs
#'
#' The four run directories of a forward-simulator run, dispatched on the
#' `--engine` flag of 07c_simulator_postprocess.R / 07d_simulator_figures.R
#' (v1 -> make_sim_run_file() kinds, spatpomp -> make_spatpomp_run_file() kinds).
#' One copy of the dispatch for both scripts, so figures always resolve against the
#' layout --aggregate/--validate wrote.
#' Pure path algebra -- creates nothing.
#'
#' @param opt resolved options list
#' @param engine "v1" (default) | "spatpomp" -- callers resolve the CLI's
#'   `opt$engine %||% "v1"` themselves (kept at the CLI, where the flag lives)
#' @param merged_dir optional override of the merged dir (07c's --merged-dir)
#' @return list(run_dir, parts_dir, merged_dir, valid_dir)
make_sim_engine_dirs <- function(opt, engine = "v1", merged_dir = NULL) {
  if (length(engine) != 1L || !engine %in% c("v1", "spatpomp")) {
    stop("--engine must be 'v1' or 'spatpomp', got: ", engine)
  }
  if (identical(engine, "spatpomp")) {
    list(
      run_dir = make_run_file(opt, kind = "spatpomp"),
      parts_dir = make_run_file(opt, kind = "spatpomp_parts_dir"),
      merged_dir = merged_dir %||%
        make_run_file(opt, kind = "spatpomp_merged_dir"),
      valid_dir = make_run_file(opt, kind = "spatpomp_valid_dir")
    )
  } else {
    list(
      run_dir = make_run_file(opt, kind = "sim"),
      parts_dir = make_run_file(opt, kind = "sim_parts_dir"),
      merged_dir = merged_dir %||% make_run_file(opt, kind = "sim_merged_dir"),
      valid_dir = make_run_file(opt, kind = "sim_valid_dir")
    )
  }
}

#' partition_sentinels
#'
#' The `_SUCCESS` sentinel files under a partitioned output tree, i.e. the
#' partitions whose writer finished (`dirname()` gives the partition dirs).
#' One scan for every reader of a partition tree: 07c's aggregate_partitions(),
#' results_core.R's load_sim_run_monthly(), 09's smoothing-events loader
#' Each caller keeps its own "none found" handling.
#'
#' @param parts_dir root of the partition tree
#' @param tus optional character vector of TU labels: scan only the existing
#'   `tu=<TU>` subdirectories (stale partitions of another scope are ignored);
#'   NULL scans everything.
#' @return character vector of sentinel paths (possibly empty)
partition_sentinels <- function(parts_dir, tus = NULL) {
  if (is.null(tus)) {
    return(list.files(
      parts_dir,
      pattern = "_SUCCESS$",
      recursive = TRUE,
      full.names = TRUE
    ))
  }
  tu_dirs <- file.path(parts_dir, paste0("tu=", tus))
  tu_dirs <- tu_dirs[dir.exists(tu_dirs)]
  if (length(tu_dirs) == 0L) {
    return(character(0))
  }
  unlist(lapply(tu_dirs, function(d) {
    list.files(d, pattern = "_SUCCESS$", recursive = TRUE, full.names = TRUE)
  }))
}

#' Path builder for the LP x gid population-crosswalk pipeline.
#'
#' @param opt Resolved options list; uses `opt$outdir`, `opt$population`.
#' @param kind One of "pop_crosswalk", "pop_gid", "pop_wpp", "pop_tasks_dir",
#'   "pop_parts_dir", "pop_gid_cache_dir", "pop_qc".
#' @return Character path.
#'   \describe{
#'     \item{pop_crosswalk}{long H[l,g,t] + P_l, Q_g, w, c, gid_cover, n_gid_eff}
#'     \item{pop_gid}{canonical per-(gid, year) population Q_g}
#'     \item{pop_wpp}{cached WPP2024 national totals (for out-of-range rescale)}
#'     \item{pop_tasks_dir}{per-(country, year) task YAMLs}
#'     \item{pop_parts_dir}{per-task outputs, iso3=<ISO3>/year=<YYYY>/}
#'     \item{pop_gid_cache_dir}{vintage-independent Q_g cache of the per-country
#'       runner, iso3=<ISO3>/year=<YYYY>.qs}
#'     \item{pop_qc}{QC/reconciliation reports}
#'   }
#' @details Pipeline: `02a_preprocess_population_tasks.R` ->
#'   `02b_preprocess_population_run.R` -> `02c_preprocess_population_compile.R`.
#'
#'   The compiled artefacts are keyed ONLY by the extraction year range, not by
#'   any model/run tag. This follows the deliberate precedent of the
#'   `size_mean_rate` table (see `make_size_run_file()`'s early return): the
#'   crosswalk is a property of the geography and the WorldPop rasters, not of
#'   any particular fit, so every fit and every simulator run must read the
#'   SAME artefact rather than each deriving its own. Overriding
#'   `population.crosswalk_path` / `population.gid_pop_path` in a config
#'   bypasses this and is the supported escape hatch for a side-by-side
#'   comparison.
make_pop_run_file <- function(
  opt,
  kind = c(
    "pop_crosswalk",
    "pop_gid",
    "pop_wpp",
    "pop_tasks_dir",
    "pop_parts_dir",
    "pop_gid_cache_dir",
    "pop_qc"
  )
) {
  kind <- match.arg(kind)
  outdir <- opt$outdir %||% "generated_data/"
  pop <- opt$population %||% list()

  yrs <- pop$years %||% c(2015L, 2030L)
  yrs <- as.integer(unlist(yrs))
  tag <- sprintf("%d-%d", min(yrs), max(yrs))

  run_dir <- file.path(outdir, "population")

  # Explicit config paths win — the escape hatch for side-by-side artefacts.
  assert_basename_length(
    switch(
      kind,
      pop_crosswalk = pop$crosswalk_path %||%
        file.path(outdir, sprintf("population_crosswalk_%s.qs", tag)),
      pop_gid = pop$gid_pop_path %||%
        file.path(outdir, sprintf("gid_population_%s.qs", tag)),
      pop_wpp = pop$wpp_table %||%
        file.path(outdir, "wpp2024_national_totals.csv"),
      pop_tasks_dir = file.path(run_dir, "tasks"),
      pop_parts_dir = file.path(run_dir, "parts"),
      pop_gid_cache_dir = file.path(run_dir, "gid_pop_cache"),
      pop_qc = file.path(outdir, "population_qc")
    ),
    "make_pop_run_file"
  )
}

#' Universal artefact-path router: dispatches to the per-track builders.
#'
#' @param opt Resolved options list, passed through to the matched builder.
#' @param kind Kind string; the `size_`/`dur_`/`mon_`/`sim`/`smooth`/`spatpomp`/`pop_`
#'   prefix selects the builder, everything else falls through to
#'   `make_occ_run_file()`.
#' @param ... Passed through to `make_occ_run_file()`.
#' @return Character path, from whichever builder `kind` routes to.
#' @details Intercepts `size_*`, `dur_*`, `mon_*`, `sim*`, `smooth*`, `spatpomp*` and
#'   `pop_*` kinds and forwards to their helpers. Occurrence kinds fall through
#'   to `make_occ_run_file()` (defined earlier in this file). This is what lets
#'   `05_fit_postprocess.R` build every track's paths through one
#'   `make_run_file(opt, paste0(kpfx, "fit"))` call.
make_run_file <- function(opt, kind, ...) {
  if (is.character(kind)) {
    if (grepl("^size_", kind)) {
      return(make_size_run_file(opt, kind = kind))
    }
    if (grepl("^dur_", kind)) {
      return(make_duration_run_file(opt, kind = kind))
    }
    if (grepl("^mon_", kind)) {
      return(make_monthly_run_file(opt, kind = kind))
    }
    if (kind == "sim" || grepl("^sim_", kind)) {
      return(make_sim_run_file(opt, kind = kind))
    }
    if (kind == "smooth" || grepl("^smooth_", kind)) {
      return(make_smooth_run_file(opt, kind = kind))
    }
    if (kind == "spatpomp" || grepl("^spatpomp_", kind)) {
      return(make_spatpomp_run_file(opt, kind = kind))
    }
    if (grepl("^pop_", kind)) {
      return(make_pop_run_file(opt, kind = kind))
    }
  }
  make_occ_run_file(opt, kind = kind, ...)
}

#' Config YAML for one (group, lag) of a grouped single-lag ENSO batch:
#' `<batch_dir>/<batch>_sl-<lag>_grp-<group sans non-alphanumerics>.yml`
#' (e.g. group "enso_east" -> `_grp-ensoeast`; unlike the ranking-cache
#' filenames, which keep the underscores).
#'
#' Matches on the `_sl-<lag>_grp-<group>.yml` suffix rather than building the
#' full stem, so a batch dir whose configs carry a different stem still
#' resolves. `manifest.yml` is never a candidate.
#'
#' @return Path, or NULL with a warning when no config matches.
grouped_lag_config_path <- function(batch_dir, group, lag) {
  grp_fname <- gsub("[^A-Za-z0-9]", "", group)
  pattern <- sprintf("_sl-%s_grp-%s\\.yml$", sl_token(lag), grp_fname)
  candidates <- list.files(batch_dir, pattern = "\\.yml$", full.names = TRUE)
  candidates <- candidates[basename(candidates) != "manifest.yml"]
  matches <- candidates[grepl(pattern, basename(candidates))]
  if (!length(matches)) {
    warning(
      sprintf(
        "No config found for group=%s lag=%d in %s",
        group,
        as.integer(lag),
        batch_dir
      ),
      call. = FALSE
    )
    return(NULL)
  }
  if (length(matches) > 1L) {
    warning("Multiple matches -- using first: ", matches[1L], call. = FALSE)
  }
  matches[1L]
}

# Deployed best-lag fits ----------------------------------------------------

# Filename prefix of each deployed track's fit (backend-agnostic: brms and
# brms_cb both match). Same prefixes as utils/simulations.R::.discover_group_fit();
# only `_fit_` files match (re_stats / loo / cb_meta are excluded).
# Fit-file stem of each simulator track -- the one table behind both
# .DEPLOYED_FIT_PREFIX (regex, read_deployed_lag()) and simulations.R's
# .discover_group_fit() (glob).
.FIT_FILE_STEM <- c(
  occ = "outbreaks_survival",
  dur = "outbreak_duration",
  mon = "outbreak_monthly"
)
.DEPLOYED_FIT_PREFIX <- setNames(
  paste0("^", .FIT_FILE_STEM, "_.*_fit_"),
  names(.FIT_FILE_STEM)
)

#' Deployed single lag of one (group, track): the `sl<N>` tag of the fit file
#' present in the deployed best-lag dir (deploy.fits.dest == the simulator's
#' grouping.fits_dir) -- the lag the grouped simulator loads at runtime. Used by
#' the forecast / hindcast generators for the horizon arithmetic.
#'
#' The episode-key filter mirrors simulations.R::.discover_group_fit(): a
#' gid_season fit carries "_ep-gs" (.episode_key_path_token()) and is invisible to
#' an lp_season lookup and vice versa, so both may share one dir.
#'
#' @param dir         Deployed best-lag fit dir.
#' @param group       Group label (e.g. "east").
#' @param grouping    Grouping source (e.g. "enso_grouping").
#' @param kind        Track: "occ", "dur" or "mon".
#' @param episode_key "lp_season" (default) or "gid_season".
#' @param national_zero_fill occurrence only: TRUE looks up the "_nz" fit (national
#'   zero-fill panel), FALSE (default) the fit without it.
#' @param national_zero_threshold occurrence only, with national_zero_fill: a
#'   positive threshold looks up the "_nz<threshold>" fit instead of "_nz".
#' @param occ_panel occurrence only: an occurrence-panel spec (occ_panel_spec_of() /
#'   .occ_panel_spec_from_episodes()); when given it replaces the two arguments above
#'   and carries the panel variant fields.
#' @return list(lag = <int>, file = <path>), or NULL when no fit matches; stops
#'   when more than one matches or the sl<N> tag cannot be parsed.
read_deployed_lag <- function(
  dir,
  group,
  grouping,
  kind = c("occ", "dur", "mon"),
  episode_key = "lp_season",
  national_zero_fill = FALSE,
  national_zero_threshold = 0L,
  occ_panel = NULL
) {
  kind <- match.arg(kind)
  pat <- sprintf(
    "%s.*_grp-%s-%s_.*\\.rds$",
    .DEPLOYED_FIT_PREFIX[[kind]],
    grouping,
    group
  )
  files <- list.files(dir, pattern = pat, full.names = TRUE)
  files <- files[
    grepl("_ep-gs", basename(files), fixed = TRUE) ==
      identical(episode_key, "gid_season")
  ]
  if (identical(kind, "occ")) {
    files <- files[
      if (is.null(occ_panel)) {
        .national_zero_fill_file_match(
          files,
          national_zero_fill,
          national_zero_threshold
        )
      } else {
        .occ_panel_variant_file_match(files, occ_panel)
      }
    ]
  }
  if (!length(files)) {
    return(NULL)
  }
  if (length(files) > 1L) {
    stop(
      "read_deployed_lag: multiple ",
      kind,
      " fit files matched for group '",
      group,
      "' in ",
      dir,
      " (episodes.key = ",
      episode_key,
      "): ",
      paste(basename(files), collapse = ", "),
      call. = FALSE
    )
  }
  f <- files[[1L]]
  # sl<N> = lag N; slm<k> = lead k, i.e. lag -k (base.R::sl_tag()). Leads are
  # deployable under deploy.lag_selection.rule best_elpd.
  tag <- sub(".*_sl(m?[0-9]+)(?:_|\\.rds$).*", "\\1", basename(f), perl = TRUE)
  if (!grepl("^m?[0-9]+$", tag)) {
    stop(
      "read_deployed_lag: could not parse sl<N> / slm<k> tag from ",
      basename(f),
      call. = FALSE
    )
  }
  lag <- if (startsWith(tag, "m")) {
    -as.integer(substring(tag, 2L))
  } else {
    as.integer(tag)
  }
  list(lag = lag, file = f)
}
