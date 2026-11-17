# provenance.R — provenance stamps and guards: population, prep policy,
#   episode definition, size/monthly estimand, generated configs, step counts
#
# Provides:
#   pop_provenance_stamp(), assert_pop_provenance(),
#   assert_fit_pop_provenance()     population stamp (mean_rate ordering guard)
#   prep_policy_stamp(), persist_prep_policy_stamp(),
#   read_prep_policy_stamp(), prep_policy_stamp_path()   prep policy stamp
#   episode_stamp_of(), assert_episode_key(), sim_episode_check(),
#   apply_episode_key_to_sim_config()   cross-track episode fingerprint
#   prep_cli_options(), validate_prep_cli()   shared 03_prepare_* CLI surface
#   estimand_spec(), build_size_estimand_stamp(), assert_size_estimand()
#   build_mean_rate_table(), resolve_mean_rate(), persist_mean_rate_table()
#   generated_yaml_header(), write_yaml_with_header()   generated configs
#   write_step_counts(), read_step_counts(), software_versions()
#
# Details:
#   Each stamp is written beside its artefact and checked by the reader, so a
#   fit or a simulation cannot silently mix inputs built under different
#   policies. Sources nothing.

#' fit_training_mean_rate
#'
#' RETIRED — not called by any production loader. Recovers the per-fit
#' `mean_rate` scalar baked into a size fit's training data. `mean_rate` is
#' resolved per unit from a country x admin-level table, which already encodes
#' group differences, so no loader needs this override.
#'
#' Kept side-by-side (no-deprecation rule). Hard-errors on a non-constant
#' `mean_rate` column rather than collapsing it to its median, which would
#' silently reintroduce a single scalar.
#'
#' @param fit a fitted size model carrying `$data$mean_rate`.
#' @return numeric scalar (the single training mean_rate) or NA_real_ when the
#'   column is absent (older fits).
fit_training_mean_rate <- function(fit) {
  d <- tryCatch(fit$data, error = function(e) NULL)
  if (is.null(d) || is.null(d$mean_rate)) {
    return(NA_real_)
  }
  mr_vals <- as.numeric(d$mean_rate)
  mr_uniq <- unique(mr_vals[is.finite(mr_vals)])
  if (length(mr_uniq) > 1L) {
    stop(
      "fit_training_mean_rate: fit$data$mean_rate is not constant (",
      length(mr_uniq),
      " distinct values) -- this fit was trained with a ",
      "per-unit country x admin-level mean_rate (D1), which this retired ",
      "function cannot collapse to one scalar without silently discarding ",
      "the per-unit resolution. Do not call fit_training_mean_rate() on fits ",
      "trained after the mean-rate-table change."
    )
  }
  mr <- if (length(mr_uniq) == 1L) mr_uniq else NA_real_
  if (!is.finite(mr)) {
    return(NA_real_)
  }
  mr
}

# ===========================================================================
# Population provenance stamp -- the mean_rate ordering guard
# ===========================================================================
#
# build_mean_rate_table() divides by `pop`, so a population change alters every
# cell, but the table is persisted once and read as a frozen sidecar by the
# fits and simulator_core.R. The stamp records which populations built each
# artefact, so a new population never pairs silently with an old rate table.
# It is a content digest, not a timestamp: rebuilding from unchanged inputs or
# copying between machines must not invalidate it.

#' pop_provenance_stamp
#'
#' Summarise the population column of a dataset into a comparable fingerprint.
#'
#' @param pop numeric population vector.
#' @param pop_source "legacy" (upstream LP `pop` carried onto the centroid gid)
#'   or "crosswalk" (apportioned via `H[l,g]`, footprints reconciled against
#'   `Q_g`). This is the field that catches the ordering error in the common
#'   case, before any digest arithmetic is needed.
#' @param extra optional named list merged into the stamp (e.g. crosswalk path).
#'
#' @return a list: `pop_source`, `n`, `pop_total`, `pop_digest`, `stamped_at`.
pop_provenance_stamp <- function(pop, pop_source, extra = list()) {
  pop <- as.numeric(pop)
  ok <- is.finite(pop)
  c(
    list(
      pop_source = as.character(pop_source),
      n = length(pop),
      n_finite = sum(ok),
      pop_total = if (any(ok)) sum(pop[ok]) else NA_real_,
      # Round before digesting: identical inputs must produce an identical digest
      # across platforms, and the last bits of a double are not portable.
      pop_digest = digest_or_null(round(sort(pop[ok]), 6)),
      stamped_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ),
    extra
  )
}

#' digest_or_null
#' `digest` is not a hard dependency of this project. Where it is unavailable
#' the stamp degrades to NA and comparison falls back to (pop_source, n,
#' pop_total) -- which still catches the legacy/crosswalk swap, the case this
#' guard exists for. A missing optional package must weaken a check, never
#' break the pipeline.
digest_or_null <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(NA_character_)
  }
  digest::digest(x, algo = "xxhash64")
}

#' pop_stamp_sidecar_path / persist_pop_stamp_sidecar / current_pop_digest_token
#'
#' Content-addresses the `size_mean_rate` sidecar path (`make_size_run_file()`)
#' by the population that built it, so a population change gives a new path and
#' the table rebuilds instead of `persist_mean_rate_table()`'s
#' `file.exists()` skip reusing a stale one. The digest must come from `opt`
#' alone (most readers never load the size table), so `03_prepare_siz.R`, the
#' only writer of that population, persists the stamp here right after
#' computing it.
#'
#' `current_pop_digest_token()` returns `NA_character_` (the pre-fix filename,
#' no pop token) when the sidecar is missing or `digest` was unavailable:
#' fail-open, like `digest_or_null()`. The digest covers `pop` only; a change
#' elsewhere in the size table with `pop` unchanged is caught only by
#' `persist_mean_rate_table()`'s reconciliation assert.
pop_stamp_sidecar_path <- function(opt) {
  # Episode-key token ("" by default): a gid_season size frame carries a
  # different population vector, so it gets its own stamp beside its own frame.
  file.path(
    opt$outdir %||% "generated_data/",
    paste0(
      "outbreak_size_analysis_data",
      .episode_key_path_token(opt),
      ".pop_stamp.rds"
    )
  )
}

persist_pop_stamp_sidecar <- function(opt, stamp) {
  path <- pop_stamp_sidecar_path(opt)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(stamp, path)
  invisible(path)
}

current_pop_digest_token <- function(opt) {
  path <- pop_stamp_sidecar_path(opt)
  if (!file.exists(path)) {
    return(NA_character_)
  }
  st <- tryCatch(readRDS(path), error = function(e) NULL)
  dg <- st$pop_digest
  if (is.null(dg) || is.na(dg) || !nzchar(dg)) {
    return(NA_character_)
  }
  paste0("pop-", dg)
}

#' attach_pop_provenance / get_pop_provenance
#' Carried as an attribute so the stamp travels with the artefact through
#' `saveRDS()`/`readRDS()` without changing any column schema -- no consumer
#' that does not know about the stamp can be broken by its presence.
attach_pop_provenance <- function(x, stamp) {
  attr(x, "pop_provenance") <- stamp
  x
}

get_pop_provenance <- function(x) attr(x, "pop_provenance", exact = TRUE)

#' attach_results_provenance / get_results_provenance
#'
#' The phase-09 (`09b_postprocess_results.R`) analogue of the pop stamp above,
#' and carried the same way and for the same reason: as an attribute, so the
#' stamp travels with the derived table through `qs2::qs_save()`/`qs_read()`
#' without changing a single column, and no consumer that predates the stamp
#' can be broken by its presence.
#'
#' What it records, and why it must exist: phase 10 plots tables it does not
#' derive, and the freshness rule it inherited (`read_cache()`, mtime-only) is
#' structurally blind to a changed *decision*. Edit `prob_threshold` and re-run,
#' and every input file is still older than the cache, so the stale table is
#' silently reused. That blindness is exactly what the hand-maintained
#' `v2`/`v3`/`v4` filename tokens were invented to paper over, one manual bump
#' at a time. Stamping the resolved `results` config block onto the artefact
#' lets `results_require()` (analysis/utils/results_core.R) compare decisions
#' rather than timestamps, and name the key that differs.
#'
#' `stamp` is built by `.results_stamp()`; see it for the recorded fields.
attach_results_provenance <- function(x, stamp) {
  attr(x, "results_provenance") <- stamp
  x
}

get_results_provenance <- function(x) {
  attr(x, "results_provenance", exact = TRUE)
}

#' attach_deploy_provenance / get_deploy_provenance
#'
#' One-line mirrors of the results stamp above, for the prod_run<N> deployment
#' orchestrator (analysis/utils/deploy_prod_run.R). A DISTINCT attribute name
#' ("deploy_provenance") so a deploy manifest and a results artefact can never
#' collide if the two ever share a container. The stamp records the resolved
#' `deploy` block + config path + schema token + derivation time + R version +
#' the input/generated tables; see deploy_core.R's manifest builder for fields.
attach_deploy_provenance <- function(x, stamp) {
  attr(x, "deploy_provenance") <- stamp
  x
}

get_deploy_provenance <- function(x) attr(x, "deploy_provenance", exact = TRUE)

#' pop_provenance_of
#'
#' Layer 2 of the mean_rate sidecar fragility fix (see `pop_stamp_sidecar_path()`):
#' the SINGLE canonical way to read a stamp back off a `build_mean_rate_table()`
#' result. Structural (`x$meta$pop_provenance`) is checked FIRST and is
#' authoritative -- `build_mean_rate_table()` always sets it as plain list data,
#' so it cannot be silently stripped by an attribute-dropping operation the
#' way `attr(x, "pop_provenance")` (set by `attach_pop_provenance()`) could be
#' if a future refactor ever routed `mrt` through a verb that rebuilds the
#' object. The attribute is checked only as a fallback, for artefacts old
#' enough to predate `$meta$pop_provenance` but that still carry the
#' attribute (round-tripping an attribute through saveRDS()/readRDS() is
#' lossless -- see analysis/tests/test_pop_provenance_guard.R section E --
#' so this only matters for pre-migration objects still in memory).
pop_provenance_of <- function(x) {
  if (!is.null(x[["meta"]][["pop_provenance"]])) {
    return(x$meta$pop_provenance)
  }
  get_pop_provenance(x)
}

#' compare_pop_provenance
#'
#' @return `character(0)` when the two stamps are compatible, otherwise one
#'   string per disagreeing field.
#'
#' An ABSENT stamp on either side is not a mismatch. Artefacts built before
#' this guard existed carry no stamp, and treating that as a failure would make
#' the guard's introduction a breaking change for every existing run -- the
#' caller decides whether unstamped is acceptable (see `assert_pop_provenance()`).
compare_pop_provenance <- function(a, b) {
  if (is.null(a) || is.null(b)) {
    return(character(0))
  }
  out <- character(0)
  if (!identical(a$pop_source, b$pop_source)) {
    out <- c(out, sprintf("pop_source: %s vs %s", a$pop_source, b$pop_source))
  }
  # Only compare digests when both are real; NA means `digest` was unavailable.
  if (
    !is.na(a$pop_digest %||% NA) &&
      !is.na(b$pop_digest %||% NA) &&
      !identical(a$pop_digest, b$pop_digest)
  ) {
    out <- c(out, sprintf("pop_digest: %s vs %s", a$pop_digest, b$pop_digest))
  } else if (
    (is.na(a$pop_digest %||% NA) || is.na(b$pop_digest %||% NA)) &&
      is.finite(a$pop_total %||% NA_real_) &&
      is.finite(b$pop_total %||% NA_real_) &&
      !isTRUE(all.equal(a$pop_total, b$pop_total, tolerance = 1e-6))
  ) {
    out <- c(out, sprintf("pop_total: %.0f vs %.0f", a$pop_total, b$pop_total))
  }
  out
}

#' assert_pop_provenance
#'
#' Hard-fail when a mean_rate sidecar was built from a different population
#' than the data now being paired with it.
#'
#' @param require_stamp when TRUE an artefact with no stamp is itself fatal.
#'   Default TRUE (Layer 3 of the mean_rate sidecar fragility fix -- see
#'   `pop_stamp_sidecar_path()`): an unverified pairing halts rather than
#'   warning. This is safe because Layer 1 (population folded into `make_size_run_file()`'s path,
#'   above) means every artefact built from here on always carries a stamp --
#'   `build_mean_rate_table()` stamps unconditionally, falling back to
#'   `pop_source = "unstamped"` rather than omitting the stamp entirely. An
#'   unstamped artefact reaching this assert therefore means a genuinely
#'   pre-Layer-1 artefact is still on disk; pass `require_stamp = FALSE`
#'   explicitly at any call site that must stay permissive for such artefacts.
assert_pop_provenance <- function(
  data_stamp,
  sidecar_stamp,
  what,
  require_stamp = TRUE
) {
  if (is.null(sidecar_stamp) || is.null(data_stamp)) {
    msg <- paste0(
      "population provenance is UNVERIFIED for ",
      what,
      " (",
      if (is.null(data_stamp)) "the data " else "",
      if (is.null(data_stamp) && is.null(sidecar_stamp)) "and " else "",
      if (is.null(sidecar_stamp)) "the sidecar " else "",
      "carries no stamp). Rebuild both after the population fix to restore ",
      "the check: pop fix -> mean_rate rebuild -> size refit."
    )
    if (isTRUE(require_stamp)) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  diffs <- compare_pop_provenance(data_stamp, sidecar_stamp)
  if (length(diffs)) {
    stop(
      "population provenance MISMATCH for ",
      what,
      ":\n  ",
      paste(diffs, collapse = "\n  "),
      "\n",
      "The mean_rate table is computed as suspected_cases / pop, so a ",
      "population change invalidates EVERY cell of it. This sidecar was ",
      "built from a different population than the data it is now being ",
      "paired with -- new populations against old rates.\n",
      "Required order: pop fix -> mean_rate rebuild ",
      "(delete the sidecar and re-run the size prep/fit) -> size refit -> ",
      "Phase-A smoothing rebuild -> forecast.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' fit_pop_stamp_path / persist_fit_pop_stamp / assert_fit_pop_provenance
#'
#' A population stamp next to each size FIT, whose path is not keyed on
#' population, so a fit trained on a stale population is not silently reused.
#' `persist_fit_pop_stamp()` is called only from the fresh-fit branch of each
#' `04_fit_siz_*.R` script, never from a skip-reuse branch (that would credit
#' the current population to a fit that never saw it).
#' `assert_fit_pop_provenance()` defaults to `require_stamp = FALSE`, narrower
#' than `assert_pop_provenance()`, because fits from before the stamp existed
#' have none; pass `TRUE` where every fit is known to carry one.
fit_pop_stamp_path <- function(fit_path) {
  paste0(tools::file_path_sans_ext(fit_path), ".pop_stamp.rds")
}

persist_fit_pop_stamp <- function(fit_path, stamp) {
  path <- fit_pop_stamp_path(fit_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(stamp, path)
  invisible(path)
}

assert_fit_pop_provenance <- function(
  opt,
  fit_path,
  what,
  require_stamp = FALSE
) {
  stamp_path <- fit_pop_stamp_path(fit_path)
  fit_stamp <- if (file.exists(stamp_path)) {
    tryCatch(readRDS(stamp_path), error = function(e) NULL)
  } else {
    NULL
  }
  current_sidecar_path <- pop_stamp_sidecar_path(opt)
  current_stamp <- if (file.exists(current_sidecar_path)) {
    tryCatch(readRDS(current_sidecar_path), error = function(e) NULL)
  } else {
    NULL
  }
  assert_pop_provenance(
    data_stamp = current_stamp,
    sidecar_stamp = fit_stamp,
    what = what,
    require_stamp = require_stamp
  )
}

#' prep_policy_stamp_path / persist_prep_policy_stamp / read_prep_policy_stamp
#'
#' Sidecar recording which `--composite_policy` built a `03_prepare_*.R` output.
#' The policy changes a prep file's content but not its name, and renaming is
#' not practical (`make_surv_dat_file()` has many bare `list(baseline = ...)`
#' callers), so the stamp makes a policy mismatch name its cause instead of
#' surfacing as a bare `pop_digest` mismatch. A file survives where the
#' `attach_pop_provenance()` attribute is stripped.
#'
#' Exception: `episodes.key` changes what an episode is, so `gid_season` prep
#' files carry an `_ep-gs` path token and coexist with the default ones. The
#' key and an episode-table digest (`episode_key`, `episode_digest`) are also
#' stamped here; `assert_episode_coherence()` compares them.
prep_policy_stamp_path <- function(path) {
  paste0(tools::file_path_sans_ext(path), ".policy.rds")
}

persist_prep_policy_stamp <- function(path, stamp) {
  out <- prep_policy_stamp_path(path)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(stamp, out)
  invisible(out)
}

read_prep_policy_stamp <- function(path) {
  out <- prep_policy_stamp_path(path)
  if (!file.exists(out)) {
    return(NULL)
  }
  tryCatch(readRDS(out), error = function(e) NULL)
}

#' assert_prep_stamp_overwrite
#'
#' Overwrite guard of a prep artefact (03_prepare_occ.R): when the
#' artefact at `path` exists, its prep stamp must equal `new_stamp` in every
#' field but `stamped_at`; otherwise the write is refused (naming the differing
#' fields) unless `force`. An existing artefact without a stamp counts as
#' different. A columns-only change (e.g. --contagion_coverage_report) keeps the
#' filename, so without this guard it would replace a deployed run's panel.
#' It is reused by the other 03 scripts
#' that write fixed-name artefacts (siz, dur, mon, the three sim panels, the
#' WHO monthly table): `ignore` names the volatile fields beside `stamped_at`
#' (03_prepare_dur.R's `restamped_at`), and `old_stamp` supplies the existing
#' artefact's record when it is not a `.policy.rds` sidecar (the WHO table's,
#' built from its `_meta.yml`). Both defaults give the occurrence behaviour.
#' @return invisibly, the differing field names (character(0) when none).
assert_prep_stamp_overwrite <- function(
  path,
  new_stamp,
  force = FALSE,
  rerun_hint = "",
  ignore = "stamped_at",
  old_stamp = read_prep_policy_stamp(path)
) {
  if (!file.exists(path)) {
    return(invisible(character(0)))
  }
  old <- old_stamp
  strip <- function(s) s[setdiff(names(s %||% list()), ignore)]
  a <- strip(old)
  b <- strip(new_stamp)
  diff <- if (is.null(old)) {
    "<no stamp>"
  } else {
    sort(union(
      setdiff(union(names(a), names(b)), intersect(names(a), names(b))),
      Filter(
        function(k) !identical(a[[k]], b[[k]]),
        intersect(names(a), names(b))
      )
    ))
  }
  if (!length(diff)) {
    return(invisible(character(0)))
  }
  msg <- sprintf(
    paste0(
      "%s exists and was prepped differently (field(s) %s). Refusing to ",
      "overwrite it%s."
    ),
    path,
    paste(utils::head(diff, 12), collapse = ", "),
    if (nzchar(rerun_hint)) {
      paste0("; to replace it on purpose: ", rerun_hint)
    } else {
      ""
    }
  )
  if (!isTRUE(force)) {
    stop(msg, call. = FALSE)
  }
  message("[overwrite --force] ", msg)
  invisible(diff)
}

#' prep_force_option
#'
#' The `--force` flag of the 03 scripts guarded by assert_prep_stamp_overwrite()
#' (03_prepare_occ.R keeps its own, identical option). `what`
#' names the artefact in the help text.
prep_force_option <- function(what = "output") {
  optparse::make_option(
    "--force",
    default = FALSE,
    action = "store_true",
    help = sprintf(
      "overwrite an existing %s built under a different prep stamp (or with none) [default: %%default]",
      what
    )
  )
}

#' who_monthly_guard_stamp
#'
#' The record assert_prep_stamp_overwrite() compares for the WHO monthly table
#' (03_prepare_who_monthly.R), which has a `_meta.yml` sidecar but no
#' `.policy.rds`: the meta without its volatile fields (`written`, the inputs'
#' `mtime`) plus a content digest of the table. NULL when `meta` is NULL (an
#' existing table without its meta is refused, as one without a stamp).
who_monthly_guard_stamp <- function(meta, data) {
  if (is.null(meta)) {
    return(NULL)
  }
  meta$written <- NULL
  for (k in names(meta$inputs %||% list())) {
    if (is.list(meta$inputs[[k]])) meta$inputs[[k]]$mtime <- NULL
  }
  # Through YAML and back, so a freshly built meta compares like one read from disk.
  meta <- yaml::yaml.load(yaml::as.yaml(meta))
  c(meta, list(content_digest = content_digest(data)))
}

#' prep_policy_stamp
#'
#' Build the stamp itself. Kept separate from the writer so both 03 scripts
#' assemble an IDENTICALLY-SHAPED record -- a reader must not have to know which
#' track wrote the file it is holding.
#'
#' @param pop_source        "legacy" or "crosswalk".
#' @param composite_policy  keep | collapse | drop | drop_collapse.
#' @param cover_threshold   the dominance threshold actually applied.
#' @param composite_tiling  "derived" | "corpus_label" — which GADM tiling the
#'                          policy was evaluated at. Decisive for reproducing a
#'                          run: the two answers differ for 938 of 4760 LPs.
#' @param counts            named list of run counts (n_lp_total, n_lp_dropped,
#'                          n_rows_in, n_rows_out, ...), merged in as-is.
#' @param pop_digest        digest of the population vector, so this stamp can be
#'                          matched to a `pop_provenance_stamp()` without opening
#'                          the (large) artefact it sits beside.
prep_policy_stamp <- function(
  pop_source,
  composite_policy,
  cover_threshold,
  composite_tiling = NA_character_,
  counts = list(),
  pop_digest = NA_character_,
  episode_key = "lp_season",
  episode_digest = NA_character_
) {
  c(
    list(
      pop_source = as.character(pop_source),
      composite_policy = as.character(composite_policy),
      cover_threshold = as.numeric(cover_threshold),
      composite_tiling = as.character(composite_tiling),
      pop_digest = as.character(pop_digest),
      episode_key = as.character(episode_key),
      episode_digest = as.character(episode_digest),
      stamped_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ),
    counts
  )
}

#' prep_stamp_inherit
#'
#' The population/composite/episode fields a DOWNSTREAM prep (duration from
#' size, monthly from duration) copies from its upstream stamp, so all four
#' tracks carry the same record shape and `assert_episode_coherence()` can name
#' the policy that produced each side. A missing upstream stamp yields NA
#' fields and the default episode key, which is exactly what those tracks
#' wrote before this helper existed.
prep_stamp_inherit <- function(upstream_path) {
  st <- read_prep_policy_stamp(upstream_path) %||% list()
  list(
    pop_source = st$pop_source %||% NA_character_,
    composite_policy = st$composite_policy %||% NA_character_,
    cover_threshold = st$cover_threshold %||% NA_real_,
    composite_tiling = st$composite_tiling %||% NA_character_,
    episode_key = st$episode_key %||% "lp_season",
    episode_digest = st$episode_digest %||% NA_character_
  )
}

# ===========================================================================
# Episode fingerprint (cross-track episode definition)
# ===========================================================================
#
# Four prep scripts, four fits and the simulator all have to agree on what one
# outbreak episode is. Under `episodes.key: gid_season` every prep output is
# stamped with a digest of the canonical (gid, year_seas, TL, TR) table built
# by utils/footprint.R::build_canonical_episodes(); every loader checks its own
# key against the config, and every consumer that pairs artefacts checks that
# their digests agree. Under the default `lp_season` nothing here can fail:
# the digest is NA and the guards return immediately unless some participant
# claims gid_season.

#' digest_episodes
#'
#' xxhash64 of the sorted `gid|year_seas|TL|TR` strings. A character vector,
#' not a data frame, so the result is independent of column order, classes and
#' attributes; gids are flattened so dotted (crosswalk) and underscore (legacy)
#' spellings of the same unit hash identically. Cases are deliberately NOT
#' hashed: the episode IDENTITY is the shared object; each track's response is
#' its own business.
digest_episodes <- function(ep) {
  need <- c("gid", "year_seas", "TL", "TR")
  miss <- setdiff(need, names(ep))
  if (length(miss)) {
    stop(
      "digest_episodes: missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  key <- paste(
    chartr(".", "_", as.character(ep$gid)),
    as.character(ep$year_seas),
    format(as.Date(ep$TL)),
    format(as.Date(ep$TR)),
    sep = "|"
  )
  digest_or_null(sort(unique(key)))
}

#' episode_stamp_of
#'
#' The episode record of a prep artefact, read from its `.policy.rds` sidecar.
#' An artefact with no sidecar, or one written before the episode fields
#' existed, reads as `lp_season` with no digest -- which is what it is.
episode_stamp_of <- function(path) {
  st <- read_prep_policy_stamp(path) %||% list()
  list(
    path = as.character(path),
    key = st$episode_key %||% "lp_season",
    digest = st$episode_digest %||% NA_character_,
    n = st$n_episodes_canonical %||% NA_integer_,
    pop_source = st$pop_source %||% NA_character_,
    composite_policy = st$composite_policy %||% NA_character_,
    cover_threshold = st$cover_threshold %||% NA_real_,
    composite_tiling = st$composite_tiling %||% NA_character_,
    stamped_at = st$stamped_at %||% NA_character_
  )
}

#' fit_episode_stamp_of
#'
#' The episode record a fit was trained on: `attr(fit, "episode_stamp")`, set
#' by the 04_fit_* scripts from the loader's record. A fit predating the field
#' reads as `lp_season`, with its path left NA.
fit_episode_stamp_of <- function(fit, label = "fit") {
  st <- tryCatch(attr(fit, "episode_stamp"), error = function(e) NULL)
  if (is.null(st)) {
    o <- tryCatch(attr(fit, "opt"), error = function(e) NULL)
    st <- list(
      key = (o$episodes %||% list())$key %||% "lp_season",
      digest = NA_character_
    )
  }
  utils::modifyList(
    list(
      path = label,
      key = "lp_season",
      digest = NA_character_,
      n = NA_integer_,
      pop_source = NA_character_,
      composite_policy = NA_character_,
      cover_threshold = NA_real_,
      composite_tiling = NA_character_,
      stamped_at = NA_character_
    ),
    st
  )
}

.episode_rerun_hint <- function(
  track,
  baseline = "1990_2019",
  key = "gid_season"
) {
  occ <- sprintf(
    "Rscript analysis/03_prepare_occ.R -b %s --episode_key %s",
    baseline,
    key
  )
  chain <- sprintf(
    paste0(
      "Rscript analysis/03_prepare_siz.R -b %s --episode_key %s\n",
      "  Rscript analysis/03_prepare_dur.R -b %s --episode_key %s\n",
      "  Rscript analysis/03_prepare_mon.R -b %s --episode_key %s"
    ),
    baseline,
    key,
    baseline,
    key,
    baseline,
    key
  )
  switch(
    track,
    occurrence = occ,
    size = ,
    duration = ,
    monthly = chain,
    paste(occ, chain, sep = "\n  ")
  )
}

#' assert_episode_key
#'
#' Loader-side guard (pattern: `assert_contagion_partition()`): the prep
#' artefact at `path` must have been built under the episode key the config
#' asks for. Returns the artefact's episode record, for the loader to attach
#' to the data (`attr(d, "episode_stamp")`) and the fit script to carry onto
#' the fit.
assert_episode_key <- function(opt, path, track) {
  want <- episode_key_of(opt)
  rec <- episode_stamp_of(path)
  if (!identical(rec$key, want)) {
    stop(
      sprintf(
        paste0(
          "%s: config asks for episodes.key = '%s' but %s was prepped under '%s'.\n",
          "Re-run:\n  %s"
        ),
        track,
        want,
        path,
        rec$key,
        .episode_rerun_hint(track, opt$baseline %||% "1990_2019", want)
      ),
      call. = FALSE
    )
  }
  if (identical(want, "gid_season") && is.na(rec$digest)) {
    stop(
      sprintf(
        paste0(
          "%s: %s claims episodes.key = 'gid_season' but carries no episode ",
          "digest. Re-run:\n  %s"
        ),
        track,
        path,
        .episode_rerun_hint(track, opt$baseline %||% "1990_2019")
      ),
      call. = FALSE
    )
  }
  rec
}

#' assert_national_zero_fill
#'
#' Loader-side guard for the occurrence panel's national zero-fill variant: the
#' artefact at `path` must have been prepped with --national_zero_fill exactly when
#' the config's episodes.national_zero_fill asks for it. The "_nz" filename token
#' already separates the two, so this catches a hand-renamed or overwritten
#' file. A stamp without the field reads as FALSE -- what it is.
#' A wrapper of assert_occ_panel_variant() (side-by-side rule).
assert_national_zero_fill <- function(opt, path, track = "occurrence") {
  assert_occ_panel_variant(opt, path, track = track)
}

#' assert_occ_panel_variant
#'
#' Generalises assert_national_zero_fill() to the whole occurrence-panel spec
#' (utils/paths.R::occ_panel_spec_of()): national zero-fill, its threshold and rule
#' the at-risk start and the ADM1-overlap treatment. Each
#' must equal the panel's prep stamp; a stamp WITHOUT a field reads as that
#' field's default (it was prepped before the field existed, i.e. at today's
#' behaviour). The fill / threshold messages are unchanged.
#' @return invisibly, the stamp's national_zero_fill (as assert_national_zero_fill did).
assert_occ_panel_variant <- function(opt, path, track = "occurrence") {
  spec <- occ_panel_spec_of(opt)
  want <- spec$national_zero_fill
  stamp <- read_prep_policy_stamp(path) %||% list()
  have <- isTRUE(stamp$national_zero_fill)
  if (!identical(want, have)) {
    stop(
      sprintf(
        paste0(
          "%s: config asks for episodes.national_zero_fill = %s but %s was prepped %s ",
          "--national_zero_fill. Re-run:\n  Rscript analysis/03_prepare_occ.R -b %s ",
          "--episode_key gid_season%s"
        ),
        track,
        tolower(want),
        path,
        if (have) "with" else "without",
        opt$baseline %||% "1990_2019",
        if (want) " --national_zero_fill" else ""
      ),
      call. = FALSE
    )
  }
  # Threshold variant (episodes.national_zero_threshold): a stamp without the
  # field was prepped before the threshold existed, i.e. at 0.
  if (want) {
    want_t <- national_zero_threshold_of(opt)
    have_t <- as.integer(stamp$national_zero_threshold %||% 0L)
    if (!identical(want_t, have_t)) {
      stop(
        sprintf(
          paste0(
            "%s: config asks for episodes.national_zero_threshold = %d but %s was prepped ",
            "with threshold %d. Re-run:\n  Rscript analysis/03_prepare_occ.R -b %s ",
            "--episode_key gid_season --national_zero_fill%s"
          ),
          track,
          want_t,
          path,
          have_t,
          opt$baseline %||% "1990_2019",
          if (want_t > 0L) paste0(" --national_zero_threshold ", want_t) else ""
        ),
        call. = FALSE
      )
    }
  }
  # Item 11/12 fields: absent from a stamp = the default.
  flag <- c(
    national_zero_rule = "--national_zero_rule",
    at_risk_start = "--at_risk_start",
    adm1_overlap = "--adm1_overlap"
  )
  for (f in names(flag)) {
    have_f <- as.character(stamp[[f]] %||% .OCC_PANEL_SPEC_DEFAULTS[[f]])
    if (identical(f, "national_zero_rule")) {
      have_f <- .national_zero_rule_check(have_f)
    }
    if (!identical(spec[[f]], have_f)) {
      stop(
        sprintf(
          paste0(
            "%s: config asks for episodes.%s = '%s' but %s was prepped with '%s'. Re-run:\n",
            "  Rscript analysis/03_prepare_occ.R -b %s --episode_key gid_season %s %s"
          ),
          track,
          f,
          spec[[f]],
          path,
          have_f,
          opt$baseline %||% "1990_2019",
          flag[[f]],
          spec[[f]]
        ),
        call. = FALSE
      )
    }
  }
  invisible(have)
}

#' assert_episode_coherence
#'
#' Consumer-side guard for anything that pairs artefacts (the simulator, the
#' smoothing filter, 07c/09/10). `records` is a named list of
#' `episode_stamp_of()` / `fit_episode_stamp_of()` outputs. Silent unless some
#' record claims `gid_season` (or `require_gid_season`); then every record must
#' be gid_season and every digest present and equal. The failure lists each
#' side's policy fields, because the usual cause is two preps run under
#' different --pop_source / --composite_policy, not a corrupt file.
assert_episode_coherence <- function(
  records,
  what = "artefacts",
  require_gid_season = FALSE
) {
  records <- Filter(Negate(is.null), records)
  if (!length(records)) {
    return(invisible(TRUE))
  }
  keys <- vapply(records, function(r) r$key %||% "lp_season", character(1))
  if (!require_gid_season && !any(keys == "gid_season")) {
    return(invisible(TRUE))
  }
  dg <- vapply(
    records,
    function(r) as.character(r$digest %||% NA_character_),
    character(1)
  )
  ok <- all(keys == "gid_season") && !anyNA(dg) && length(unique(dg)) == 1L
  if (ok) {
    return(invisible(TRUE))
  }
  # R truncates error messages at getOption("warning.length") (default 1000
  # chars), which would cut off the per-record table and the re-run commands.
  op <- options(warning.length = 8170L)
  on.exit(options(op), add = TRUE)
  fld <- function(r, f) as.character(r[[f]] %||% NA)
  tab <- vapply(
    names(records) %||% seq_along(records),
    function(nm) {
      r <- records[[nm]]
      sprintf(
        "  %-10s key=%-10s digest=%-16s pop_source=%-9s policy=%-13s thr=%s tiling=%s  (%s)",
        nm,
        fld(r, "key"),
        fld(r, "digest"),
        fld(r, "pop_source"),
        fld(r, "composite_policy"),
        fld(r, "cover_threshold"),
        fld(r, "composite_tiling"),
        fld(r, "path")
      )
    },
    character(1)
  )
  stop(
    paste0(
      "assert_episode_coherence: ",
      what,
      " do not share one episode definition.\n",
      paste(tab, collapse = "\n"),
      "\n",
      "Every participant must be episodes.key = 'gid_season' with the same ",
      "episode digest. Rebuild the four preps with the SAME --pop_source / ",
      "--composite_policy / --cover_threshold / --composite_tiling, then refit:\n  ",
      .episode_rerun_hint("all")
    ),
    call. = FALSE
  )
}

#' sim_episode_check
#'
#' One call for every multi-artefact consumer (simulator / smoothing loaders,
#' 07c, 09, 10): gathers the episode records of the fits it pairs (from the
#' fit objects) and of the prep artefacts it reads (from their stamps; paths
#' that do not exist are skipped), then runs assert_episode_coherence(). With
#' `episode_key = "gid_season"` every participant MUST be gid_season, so a
#' gid_season run cannot silently consume a historical fit or panel.
#'
#' @param fits named list of fit objects (NULL entries skipped).
#' @param preps named list of prep-artefact paths (NULL / missing skipped).
sim_episode_check <- function(
  fits = list(),
  preps = list(),
  what = "artefacts",
  episode_key = "lp_season"
) {
  recs <- list()
  for (n in names(fits)) {
    if (!is.null(fits[[n]])) {
      recs[[paste0(n, "_fit")]] <- fit_episode_stamp_of(
        fits[[n]],
        paste0(n, "_fit")
      )
    }
  }
  for (n in names(preps)) {
    p <- preps[[n]]
    if (!is.null(p) && length(p) == 1L && !is.na(p) && file.exists(p)) {
      recs[[n]] <- episode_stamp_of(p)
    }
  }
  # The run's own declared key is a participant too: a historical (lp_season)
  # config must not consume gid_season fits/panels even when those agree with
  # each other. (gid_season is already enforced by require_gid_season below.)
  if (!identical(episode_key, "gid_season") && length(recs)) {
    recs[["config"]] <- list(
      path = "episodes.key in this config",
      key = episode_key,
      digest = NA_character_
    )
  }
  assert_episode_coherence(
    recs,
    what = what,
    require_gid_season = identical(episode_key, "gid_season")
  )
}

#' apply_episode_key_to_sim_config
#'
#' Deploy-side authority for a generated simulator config's episode
#' definition (deploy_prod_run.R -> make_forecast_configs.R /
#' make_hindcast_configs.R). The deploy YAML's `deploy.episodes.key` is the one
#' declaration for a rollout; this writes it into every generated config
#' (top-level `episodes.key`) so no template can silently fall back to the
#' historical definition. Rules:
#'   * `episode_key` NULL/"" (standalone CLI without --episode_key): the
#'     template is inherited untouched -- previous behaviour.
#'   * a template that DECLARES a different key is a hard error (two
#'     authorities disagreeing is never resolved by guessing).
#'   * `lp_season` is written only if the template already carries an
#'     `episodes` block, so every default-generated config stays
#'     byte-identical to before.
#'   * `episodes.national_zero_fill` (national zero-fill occurrence panel, gid_season
#'     only): `national_zero_fill` NULL keeps the template's own value (the
#'     forecast / hindcast generators, whose templates gen_templates wrote);
#'     TRUE/FALSE is the deploy's declaration, and a template declaring the
#'     other value is a hard error. Written only when TRUE, so a config
#'     without zero-fill stays byte-identical to before.
apply_episode_key_to_sim_config <- function(
  cfg,
  episode_key,
  label = "config",
  national_zero_fill = NULL,
  national_zero_threshold = NULL,
  occ_panel = NULL
) {
  if (is.null(episode_key) || !nzchar(episode_key)) {
    return(cfg)
  }
  # occ_panel (an occurrence-panel spec, utils/paths.R) is the deploy's whole
  # declaration: it supplies the fill / threshold pair and the panel fields
  # (rule, at-risk start, ADM1 overlap), each under the same authority rule as
  # the pair. NULL keeps the template's own values.
  if (!is.null(occ_panel)) {
    national_zero_fill <- occ_panel$national_zero_fill
    national_zero_threshold <- occ_panel$national_zero_threshold
  }
  if (!episode_key %in% EPISODE_KEYS) {
    stop(
      label,
      ": episode key '",
      episode_key,
      "' is not one of {",
      paste(EPISODE_KEYS, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
  tpl_key <- (cfg$episodes %||% list())$key
  if (!is.null(tpl_key) && !identical(tpl_key, episode_key)) {
    stop(
      label,
      ": the template declares episodes.key = '",
      tpl_key,
      "' but the deploy declares '",
      episode_key,
      "'. Fix one of them -- ",
      "the deploy YAML (deploy.episodes.key) is meant to be the authority.",
      call. = FALSE
    )
  }
  tpl_nz <- (cfg$episodes %||% list())$national_zero_fill
  if (
    !is.null(national_zero_fill) &&
      !is.null(tpl_nz) &&
      !identical(isTRUE(tpl_nz), isTRUE(national_zero_fill))
  ) {
    stop(
      label,
      ": the template declares episodes.national_zero_fill = ",
      tolower(isTRUE(tpl_nz)),
      " but the deploy declares ",
      tolower(isTRUE(national_zero_fill)),
      ". Fix one of them -- ",
      "the deploy YAML (deploy.episodes.national_zero_fill) is meant to be the authority.",
      call. = FALSE
    )
  }
  nzf <- isTRUE(national_zero_fill %||% tpl_nz)
  if (nzf && !identical(episode_key, "gid_season")) {
    stop(
      label,
      ": episodes.national_zero_fill = true requires episodes.key = 'gid_season'.",
      call. = FALSE
    )
  }
  # Threshold (episodes.national_zero_threshold): same authority rule; written
  # only when > 0, so configs without it stay byte-identical.
  tpl_nzt <- (cfg$episodes %||% list())$national_zero_threshold
  if (
    !is.null(national_zero_threshold) &&
      !is.null(tpl_nzt) &&
      !identical(as.integer(tpl_nzt), as.integer(national_zero_threshold))
  ) {
    stop(
      label,
      ": the template declares episodes.national_zero_threshold = ",
      tpl_nzt,
      " but the deploy declares ",
      national_zero_threshold,
      ". Fix one of them -- ",
      "the deploy YAML (deploy.episodes.national_zero_threshold) is meant to be the authority.",
      call. = FALSE
    )
  }
  nzt <- as.integer(national_zero_threshold %||% tpl_nzt %||% 0L)
  if (nzt > 0L && !nzf) {
    stop(
      label,
      ": episodes.national_zero_threshold = ",
      nzt,
      " requires episodes.national_zero_fill = true.",
      call. = FALSE
    )
  }
  # Item 11/12 panel fields: same authority rule; written only when non-default
  # (.episodes_block()), so configs without them stay byte-identical.
  tpl_ep <- cfg$episodes %||% list()
  extra <- list()
  for (f in c("national_zero_rule", "at_risk_start", "adm1_overlap")) {
    dv <- if (!is.null(occ_panel)) occ_panel[[f]] else NULL
    tv <- tpl_ep[[f]]
    if (identical(f, "national_zero_rule") && !is.null(tv)) {
      tv <- .national_zero_rule_check(tv)
    }
    if (
      !is.null(dv) &&
        !is.null(tv) &&
        !identical(as.character(tv), as.character(dv))
    ) {
      stop(
        label,
        ": the template declares episodes.",
        f,
        " = '",
        tv,
        "' but the deploy ",
        "declares '",
        dv,
        "'. Fix one of them -- the deploy YAML (deploy.episodes.",
        f,
        ") is meant to be the authority.",
        call. = FALSE
      )
    }
    extra[[f]] <- dv %||% tv %||% .OCC_PANEL_SPEC_DEFAULTS[[f]]
  }
  spec <- .occ_panel_spec(
    key = episode_key,
    national_zero_fill = nzf,
    national_zero_threshold = nzt,
    national_zero_rule = extra$national_zero_rule,
    at_risk_start = extra$at_risk_start,
    adm1_overlap = extra$adm1_overlap
  )
  if (!identical(episode_key, "lp_season") || !is.null(cfg$episodes)) {
    cfg$episodes <- .episodes_block(spec)
  }
  cfg
}

#' preflight_sim_episode_inputs
#'
#' Config-generation-time version of the simulator's own episode check
#' (sim_episode_check() in load_artefacts / smoothing): the prep artefacts a
#' generated config will read -- its climate panel, the size frame it takes
#' population/admin level from, and the survival panel of its episode key --
#' must share one episode definition. Fails while generating configs locally,
#' not hours later in a SLURM array. Silent for a historical (lp_season)
#' config whose artefacts carry no gid_season stamp.
preflight_sim_episode_inputs <- function(
  cfg,
  label = "config",
  surv_path = NULL
) {
  ek <- episode_key_of(cfg)
  sim <- cfg$simulator %||% list()
  ep_opt <- list(
    baseline = cfg$baseline %||% "1990_2019",
    episodes = cfg$episodes
  )
  # national zero-fill panel: the survival panel the config resolves must have been
  # prepped with exactly the declared episodes.national_zero_fill.
  surv_nz <- surv_path %||% as.character(make_surv_dat_file(ep_opt))
  if (file.exists(surv_nz)) {
    assert_occ_panel_variant(
      ep_opt,
      surv_nz,
      track = paste0(label, " survival panel")
    )
  }
  tryCatch(
    sim_episode_check(
      preps = list(
        climate_panel = sim$climate_panel,
        size_data = (sim$location_frame %||% list())$size_data_path %||%
          as.character(make_size_dat_file(ep_opt)),
        surv = surv_path %||% as.character(make_surv_dat_file(ep_opt))
      ),
      what = paste0(label, " inputs"),
      episode_key = ek
    ),
    error = function(e) {
      op <- options(warning.length = 8170L)
      on.exit(options(op), add = TRUE)
      stop(
        conditionMessage(e),
        sprintf(
          paste0(
            "\nFor THIS config (episodes.key = '%s'), the usual fix is in the template ",
            "or deploy YAML, not the preps:\n",
            "  * simulator.climate_panel -> the %s panel: Rscript analysis/03_prepare_sim_panel.R ",
            "(historical) / 03_prepare_sim_panel_projection.R (forecast) --episode_key %s, ",
            "or deploy.hindcast.climate_panel for hindcasts;\n",
            "  * simulator.location_frame.size_data_path -> %s (or remove the key to use that default)."
          ),
          ek,
          ek,
          ek,
          as.character(make_size_dat_file(ep_opt))
        ),
        call. = FALSE
      )
    }
  )
}

# ===========================================================================
# Shared 03_prepare_* CLI surface
# ===========================================================================
#
# One definition of the prep scripts' population/composite options, so they
# cannot drift apart; each script appends its track-specific extras.
PREP_POP_SOURCES <- c("legacy", "crosswalk")
PREP_COMPOSITE_POLICIES <- c("keep", "collapse", "drop", "drop_collapse")
PREP_COMPOSITE_TILINGS <- c("derived", "corpus_label")

#' prep_cli_options
#'
#' @param track "occ" | "siz" (full population/composite surface) or
#'   "dur" | "mon" (baseline + episode key only; they inherit everything else
#'   from the upstream artefact's stamp).
#' @return list of optparse options.
prep_cli_options <- function(track = c("occ", "siz", "dur", "mon")) {
  track <- match.arg(track)
  mo <- optparse::make_option
  base <- list(
    mo(
      c("-b", "--baseline"),
      default = "1990_2019",
      type = "character",
      help = "climate baseline period [default: %default]"
    )
  )
  if (track %in% c("dur", "mon")) {
    return(c(
      base,
      list(
        mo(
          "--episode_key",
          default = "inherit",
          type = "character",
          help = paste0(
            "inherit | lp_season | gid_season. 'inherit' reads the ",
            "upstream artefact's stamp; an explicit value only ",
            "asserts it [default: %default]"
          )
        )
      )
    ))
  }
  c(
    base,
    list(
      mo(
        "--enso_source",
        default = "oni",
        type = "character",
        help = "ENSO index source: 'oni' (default) or 'nino34' (sensitivity) [default: %default]"
      ),
      mo(
        "--refresh_enso",
        default = FALSE,
        action = "store_true",
        help = "Re-download the ENSO index even if a cache exists [default: %default]"
      ),
      mo(
        "--min_pop",
        default = 10L,
        type = "integer",
        help = paste0(
          "Drop location periods whose WorldPop population is below ",
          "this many people (broken raster extraction guard). Size ",
          "applies it only under --episode_key gid_season, so both ",
          "tracks build episodes from the same LP set ",
          "[default: %default]"
        )
      ),
      # crosswalk is the default on every track: legacy emits no
      # denominator triple, so a legacy size prep cannot serve an
      # `estimand.denominator: gid` fit.
      mo(
        "--pop_source",
        default = "crosswalk",
        type = "character",
        help = "'legacy' (centroid + sum(pop)) or 'crosswalk' [default: %default]"
      ),
      mo(
        "--crosswalk_path",
        default = "",
        type = "character",
        help = "override for population_crosswalk_<Y0>-<Y1>.qs [default: auto]"
      ),
      mo(
        "--gid_pop_path",
        default = "",
        type = "character",
        help = "override for gid_population_<Y0>-<Y1>.qs [default: auto]"
      ),
      mo(
        "--wpp_path",
        default = "",
        type = "character",
        help = "override for wpp2024_national_totals.csv [default: auto]"
      ),
      mo(
        "--composite_policy",
        default = "drop_collapse",
        type = "character",
        help = "keep | collapse | drop | drop_collapse [default: %default]"
      ),
      mo(
        "--cover_threshold",
        default = 0.90,
        type = "double",
        help = "gid_cover strictly below this => true composite [default: %default]"
      ),
      mo(
        "--composite_tiling",
        default = "derived",
        type = "character",
        help = "derived | corpus_label [default: %default]"
      ),
      mo(
        "--episode_key",
        default = "lp_season",
        type = "character",
        help = paste0(
          "lp_season (per-location-period episodes, historical) | ",
          "gid_season (one canonical episode per gid-season, ",
          "shared by all four tracks) [default: %default]"
        )
      )
    )
  )
}

#' validate_prep_cli
#'
#' Validate the shared options and attach `opt$episodes$key` so the path
#' builders (which read `opt$episodes`) resolve the right artefact from the
#' CLI `opt` directly. For dur/mon, `--episode_key inherit` is resolved by the
#' caller from the upstream stamp (it needs the upstream path).
validate_prep_cli <- function(opt, track = c("occ", "siz", "dur", "mon")) {
  track <- match.arg(track)
  if (track %in% c("occ", "siz")) {
    if (!opt$pop_source %in% PREP_POP_SOURCES) {
      stop(
        "--pop_source must be 'legacy' or 'crosswalk', got: ",
        opt$pop_source
      )
    }
    if (!opt$composite_policy %in% PREP_COMPOSITE_POLICIES) {
      stop(
        "--composite_policy must be one of ",
        paste(PREP_COMPOSITE_POLICIES, collapse = " | "),
        ", got: ",
        opt$composite_policy
      )
    }
    if (
      !is.finite(opt$cover_threshold) ||
        opt$cover_threshold < 0 ||
        opt$cover_threshold > 1
    ) {
      stop(
        "--cover_threshold must be a number in [0, 1], got: ",
        opt$cover_threshold
      )
    }
    if (!opt$composite_tiling %in% PREP_COMPOSITE_TILINGS) {
      stop(
        "--composite_tiling must be 'derived' or 'corpus_label', got: ",
        opt$composite_tiling
      )
    }
    if (!opt$episode_key %in% EPISODE_KEYS) {
      stop(
        "--episode_key must be one of ",
        paste(EPISODE_KEYS, collapse = " | "),
        ", got: ",
        opt$episode_key
      )
    }
    if (
      identical(opt$episode_key, "gid_season") &&
        !identical(opt$pop_source, "crosswalk")
    ) {
      stop(
        "--episode_key gid_season requires --pop_source crosswalk: the ",
        "canonical admin level and the P1 denominator triple both come from ",
        "the population crosswalk."
      )
    }
    opt$episodes <- list(key = opt$episode_key)
  } else if (!opt$episode_key %in% c("inherit", EPISODE_KEYS)) {
    stop(
      "--episode_key must be one of inherit | ",
      paste(EPISODE_KEYS, collapse = " | "),
      ", got: ",
      opt$episode_key
    )
  }
  opt
}

#' resolve_inherited_episode_key
#'
#' dur/mon: the episode key comes from the upstream artefact. An explicit
#' --episode_key only asserts it.
resolve_inherited_episode_key <- function(opt, upstream_path, track) {
  up <- episode_stamp_of(upstream_path)$key
  if (
    !identical(opt$episode_key, "inherit") && !identical(opt$episode_key, up)
  ) {
    stop(
      sprintf(
        paste0(
          "%s: --episode_key %s, but the upstream artefact %s was prepped under '%s'. ",
          "Re-run the upstream prep, or drop the flag to inherit."
        ),
        track,
        opt$episode_key,
        upstream_path,
        up
      ),
      call. = FALSE
    )
  }
  opt$episodes <- list(key = up)
  opt
}

# ===========================================================================
# The size estimand stamp
# ===========================================================================
#
# The size mean mu_i = pop_i * mean_rate_{c(i)} * exp(eta_i) has both offsets
# fixed, so its absolute scale lives entirely in the two denominators. A fit
# and a mean_rate table built on different denominator UNITS (whole gid vs
# reporting footprint) scale every prediction by 1/coverage, and
# assert_pop_provenance() cannot see it (it compares digests, not meaning).
# This stamp is that semantic check; its design mirrors prep_policy_stamp().

SIZE_ESTIMAND_DENOMINATORS <- c("gid", "footprint")
SIZE_POP_ELASTICITY_MODES <- c("fixed", "free", "free_by_level")

#' Per-track estimand wiring.
#'
#' The size track got this machinery first, so every symbol here is named
#' `*Size*`. The MONTHLY track has the identical problem -- its `pop` is the
#' reporting footprint H[l,g] while the simulator replays at whole-gid Q_g, and
#' `offset(log(pop))` imposes an elasticity of exactly 1 on top of that
#' mismatch -- so it shares the same apparatus rather than growing a second
#' copy that could drift. The only per-track facts are in this table.
#'
#'   `block`     the `opt` sub-list the keys live under.
#'   `attr_name` the attribute/`$meta` slot the stamp travels in. Deliberately
#'               DISTINCT per track: a mon stamp and a size stamp describe
#'               different denominators on different panels, and letting
#'               `assert_size_estimand()` pair them across tracks would compare
#'               two things that were never supposed to match.
#'   `prep`      the script to re-run when a required population column is missing.
.ESTIMAND_SPEC <- list(
  size = list(
    block = "size",
    attr_name = "size_estimand",
    prep = "analysis/03_prepare_siz.R"
  ),
  monthly = list(
    block = "monthly",
    attr_name = "mon_estimand",
    prep = "analysis/03_prepare_dur.R then analysis/03_prepare_mon.R"
  )
)

estimand_spec <- function(track = "size") {
  spec <- .ESTIMAND_SPEC[[track]]
  if (is.null(spec)) {
    stop(
      "estimand_spec: unknown track '",
      paste(track, collapse = ","),
      "' (expected one of ",
      paste(names(.ESTIMAND_SPEC), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  c(spec, list(track = track))
}

#' resolve_size_estimand / resolve_monthly_estimand
#'
#' Resolve `opt$<block>$estimand` + `opt$<block>$pop_elasticity` into a
#' validated record. This is the ONLY sanctioned way to read those keys: every
#' consumer goes through here so the validation cannot be bypassed by reading
#' the raw option.
#'
#' Two things are deliberately NOT defaulted silently:
#'
#'  1. `rate_denominator` must EQUAL `denominator`. This is the whole point of
#'     the block, so a mismatch is a hard error rather than a warning -- the
#'     failure mode it guards is invisible in every diagnostic except total
#'     mass, which is exactly why it survived this long.
#'  2. An unrecognised value is an error, not a fallback to the default. A typo
#'     in `denominator:` that quietly resolved to "footprint" would reintroduce
#'     the bug while the config claims otherwise.
#'
#' Absent keys DO default (to today's deployed behaviour), so every existing
#' config keeps running unchanged -- on both tracks. `footprint` + `fixed` is
#' what mon has always done implicitly.
.resolve_track_estimand <- function(opt, track = "size") {
  spec <- estimand_spec(track)
  blk <- opt[[spec$block]] %||% list()
  es <- blk$estimand %||% list()

  den <- as.character(es$denominator %||% "footprint")
  rden <- as.character(es$rate_denominator %||% den)
  urm <- as.numeric(es$underreporting_multiplier %||% 1.0)
  ela <- as.character(blk$pop_elasticity %||% "fixed")
  b <- spec$block

  if (!den %in% SIZE_ESTIMAND_DENOMINATORS) {
    stop(
      sprintf(
        "%s.estimand.denominator = '%s' is not one of {%s}. Refusing to guess.",
        b,
        den,
        paste(SIZE_ESTIMAND_DENOMINATORS, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!rden %in% SIZE_ESTIMAND_DENOMINATORS) {
    stop(
      sprintf(
        "%s.estimand.rate_denominator = '%s' is not one of {%s}. Refusing to guess.",
        b,
        rden,
        paste(SIZE_ESTIMAND_DENOMINATORS, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!identical(den, rden)) {
    stop(
      sprintf(
        paste0(
          "%s.estimand: denominator = '%s' but rate_denominator = '%s'.\n",
          "  These MUST be equal: mu = pop x rate x exp(eta) with both offsets\n",
          "  fixed, so mixing units multiplies every prediction by a coverage ratio\n",
          "  that nothing in the fit can absorb. Set both to '%s' or both to '%s'."
        ),
        b,
        den,
        rden,
        den,
        rden
      ),
      call. = FALSE
    )
  }
  if (!ela %in% SIZE_POP_ELASTICITY_MODES) {
    stop(
      sprintf(
        "%s.pop_elasticity = '%s' is not one of {%s}. Refusing to guess.",
        b,
        ela,
        paste(SIZE_POP_ELASTICITY_MODES, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.finite(urm) || urm <= 0) {
    stop(
      sprintf(
        "%s.estimand.underreporting_multiplier = %s must be a positive finite number.",
        b,
        format(urm)
      ),
      call. = FALSE
    )
  }

  list(
    denominator = den,
    rate_denominator = rden,
    underreporting_multiplier = urm,
    pop_elasticity = ela,
    track = track
  )
}

#' canonical_admin_level_enabled / .apply_p3
#'
#' THE single gate + call for the canonical admin-level re-key, shared by every
#' loader and anchor view (read_size_data(), load_monthly_data(), .siz_anchor_view(),
#' .mon_anchor_view(), build_mean_rate_table(), load_surv_data(),
#' load_duration_data()), so all four tracks can key
#' `admin_level` on the same, geometry-derived definition.
#'
#' Gate: `opt[[block]]$canonical_admin_level`, defaulting -- for the two tracks
#' that carry an estimand (size, monthly) -- to "on for any non-`footprint`
#' denominator", exactly the expression the inline copies used. Occurrence and
#' duration have no estimand, so they default OFF and every existing config is
#' untouched.
#'
#' Lives in provenance.R (not population.R, where apply_canonical_admin_level()
#' is): every loader already sources provenance.R for its estimand stamp, and
#' the gate must be callable on the DEFAULT path without utils/population.R --
#' which is then needed only when the re-key actually runs (call-time).
#'
#' @param opt resolved options
#' @param block "occurrence" | "size" | "duration" | "monthly"
#' @param estimand resolved estimand (size/monthly) or NULL
canonical_admin_level_enabled <- function(opt, block, estimand = NULL) {
  dflt <- !is.null(estimand) && !identical(estimand$denominator, "footprint")
  isTRUE((opt[[block]] %||% list())$canonical_admin_level %||% dflt)
}

.apply_p3 <- function(
  d,
  opt,
  block,
  estimand = NULL,
  tag = paste0(block, "/P3")
) {
  if (!canonical_admin_level_enabled(opt, block, estimand)) {
    return(d)
  }
  gid_pop_path <- opt$population$gid_pop_path %||% opt$gid_pop_path
  if (is.null(gid_pop_path) || !nzchar(gid_pop_path)) {
    stop(
      block,
      ": canonical admin_level requested but no ",
      "population.gid_pop_path is configured. The canonical gadm_level ",
      "lives in gid_population_<Y0>-<Y1>.qs; there is no fallback -- ",
      "parsing the gid string does NOT recover the level (trailing GADM ",
      "version component).",
      call. = FALSE
    )
  }
  apply_canonical_admin_level(d, gid_pop_path, overwrite = TRUE, tag = tag)
}

resolve_size_estimand <- function(opt) .resolve_track_estimand(opt, "size")
resolve_monthly_estimand <- function(opt) {
  .resolve_track_estimand(opt, "monthly")
}

#' size_estimand_pop_column
#'
#' The single mapping from a declared denominator to the column prep emits.
#' Kept as a function, not an inline `switch`, so there is exactly one place
#' where "gid means pop_gid" is written down.
size_estimand_pop_column <- function(denominator) {
  switch(
    as.character(denominator),
    gid = "pop_gid",
    footprint = "pop_footprint",
    stop(
      sprintf("no population column mapped for denominator '%s'", denominator),
      call. = FALSE
    )
  )
}

#' resolve_size_pop_column
#'
#' Point `pop` at the population the declared estimand asks for.
#'
#' A single `pop` column whose meaning is decided upstream in prep hides the
#' `Q x rbar^H` mismatch: nothing at the fit seam can tell which population it
#' was handed. With prep emitting `pop_footprint` and `pop_gid` side by side,
#' the choice is explicit and checkable here.
#'
#' Asymmetric fallback, deliberately:
#'   - `footprint` with no `pop_footprint` column falls back to the existing
#'     `pop`, so an artefact without the split columns keeps running
#'     bit-for-bit.
#'   - `gid` with no `pop_gid` column is a HARD ERROR. Falling back would hand
#'     back a footprint population while the stamp claims whole-gid -- the exact
#'     silent-substitution failure this machinery exists to prevent.
resolve_size_pop_column <- function(d, estimand, tag = "read_size_data") {
  want <- size_estimand_pop_column(estimand$denominator)
  spec <- estimand_spec(estimand$track %||% "size")
  if (want %in% names(d)) {
    d$pop <- as.numeric(d[[want]])
    message(sprintf(
      "[%s] pop <- %s (estimand denominator '%s')",
      tag,
      want,
      estimand$denominator
    ))
    return(d)
  }
  if (identical(estimand$denominator, "gid")) {
    stop(
      "resolve_size_pop_column: ",
      spec$block,
      ".estimand.denominator = 'gid' ",
      "but the prepared data has no `pop_gid` column.\n  Re-run ",
      spec$prep,
      " (P1 emits pop_footprint/pop_gid/pop_cover). ",
      "Refusing to fall back to `pop`, whose denominator is unknown at this ",
      "seam -- that substitution is the Q x rbar^H bug.",
      call. = FALSE
    )
  }
  message(sprintf(
    paste0(
      "[%s] no `%s` column (pre-P1 artefact); using `pop` as-is, which ",
      "reproduces current behaviour."
    ),
    tag,
    want
  ))
  d
}

#' build_size_estimand_stamp
#'
#' Build the stamp. `extra` is merged in as-is so a producer can record what it
#' actually did (row counts, the resolved pop column, the fitted beta) beside
#' the declaration, without every consumer needing to know those fields exist.
build_size_estimand_stamp <- function(opt, extra = list(), track = "size") {
  es <- .resolve_track_estimand(opt, track)
  c(
    list(
      track = track,
      denominator = es$denominator,
      rate_denominator = es$rate_denominator,
      underreporting_multiplier = es$underreporting_multiplier,
      pop_elasticity = es$pop_elasticity,
      pop_column = size_estimand_pop_column(es$denominator),
      stamped_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ),
    extra
  )
}

monthly_estimand_stamp <- function(opt, extra = list()) {
  build_size_estimand_stamp(opt, extra, track = "monthly")
}

#' attach_size_estimand / get_size_estimand / size_estimand_of
#'
#' Carriage mirrors `attach_pop_provenance()`: an attribute, so no column schema
#' changes and no consumer that is unaware of the stamp can be broken by its
#' presence. `size_estimand_of()` prefers the STRUCTURAL slot for the same reason
#' `pop_provenance_of()` does -- an attribute can be dropped by a verb that
#' rebuilds the object, a `$meta` entry cannot.
attach_size_estimand <- function(x, stamp, track = "size") {
  attr(x, estimand_spec(track)$attr_name) <- stamp
  x
}

get_size_estimand <- function(x, track = "size") {
  attr(x, estimand_spec(track)$attr_name, exact = TRUE)
}

get_monthly_estimand <- function(x) get_size_estimand(x, track = "monthly")

size_estimand_of <- function(x, track = "size") {
  nm <- estimand_spec(track)$attr_name
  if (!is.null(x[["meta"]][[nm]])) {
    return(x$meta[[nm]])
  }
  get_size_estimand(x, track)
}

monthly_estimand_of <- function(x) size_estimand_of(x, track = "monthly")

size_estimand_stamp_path <- function(path) {
  paste0(tools::file_path_sans_ext(path), ".estimand.rds")
}

persist_size_estimand_stamp <- function(path, stamp) {
  out <- size_estimand_stamp_path(path)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(stamp, out)
  invisible(out)
}

read_size_estimand_stamp <- function(path) {
  out <- size_estimand_stamp_path(path)
  if (!file.exists(out)) {
    return(NULL)
  }
  tryCatch(readRDS(out), error = function(e) NULL)
}

#' compare_size_estimand
#'
#' @return `character(0)` when compatible, else one string per disagreeing
#'   field. As with `compare_pop_provenance()`, an ABSENT stamp on either side is
#'   not a mismatch -- artefacts predating this guard carry none, and making
#'   their absence fatal here would turn the guard's introduction into a
#'   breaking change for every artefact already on disk. The caller decides
#'   (see `assert_size_estimand(require_stamp = )`).
#'
#' `pop_elasticity` is compared but is NOT a dimensional field: a `free` fit
#' and a `fixed` simulator frame are still on the same units. It is reported so
#' a mismatch is legible, and the caller can choose to tolerate it.
compare_size_estimand <- function(a, b) {
  if (is.null(a) || is.null(b)) {
    return(character(0))
  }
  out <- character(0)
  for (f in c("denominator", "rate_denominator", "pop_elasticity")) {
    av <- a[[f]] %||% NA_character_
    bv <- b[[f]] %||% NA_character_
    if (!identical(as.character(av), as.character(bv))) {
      out <- c(out, sprintf("%s: %s vs %s", f, av, bv))
    }
  }
  au <- a$underreporting_multiplier %||% NA_real_
  bu <- b$underreporting_multiplier %||% NA_real_
  if (
    is.finite(au) &&
      is.finite(bu) &&
      !isTRUE(all.equal(au, bu, tolerance = 1e-9))
  ) {
    out <- c(
      out,
      sprintf("underreporting_multiplier: %s vs %s", format(au), format(bu))
    )
  }
  out
}

#' assert_size_estimand
#'
#' Hard-fail when two artefacts on the same seam declare different estimands.
#'
#' @param tolerate_elasticity when TRUE, a `pop_elasticity` disagreement alone
#'   downgrades to a warning. Defaults to TRUE because it is not a units
#'   mismatch -- see `compare_size_estimand()`. The denominator fields never
#'   downgrade.
assert_size_estimand <- function(
  a,
  b,
  what,
  require_stamp = FALSE,
  tolerate_elasticity = TRUE
) {
  if (is.null(a) || is.null(b)) {
    msg <- sprintf(
      "[%s] size-estimand stamp absent on %s side; pairing is UNVERIFIED.",
      what,
      if (is.null(a)) "the first" else "the second"
    )
    if (isTRUE(require_stamp)) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  diffs <- compare_size_estimand(a, b)
  if (!length(diffs)) {
    return(invisible(TRUE))
  }

  soft <- all(grepl("^pop_elasticity:", diffs))
  msg <- sprintf(
    "[%s] size-estimand mismatch:\n  - %s",
    what,
    paste(diffs, collapse = "\n  - ")
  )
  if (soft && isTRUE(tolerate_elasticity)) {
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  stop(
    paste0(
      msg,
      "\n",
      "  mu = pop x mean_rate x exp(eta) with both offsets fixed, so a denominator\n",
      "  disagreement rescales every prediction by a coverage ratio the fit cannot\n",
      "  absorb. Rebuild the stale side under one declared estimand."
    ),
    call. = FALSE
  )
}

#' build_mean_rate_table
#'
#' One global country x admin-level `mean_rate` table, persisted by
#' `persist_mean_rate_table()` so every fit and the simulator read the same
#' artefact. It reads the raw analysis data from `path` itself, not
#' `read_size_data()`'s output, because that admin-level filter would drop the
#' rows that key the level contrasts. The min_pop / max_attack_rate /
#' non-NA-TL-TR guards are applied as in `read_size_data()`; `opt$country` is
#' ignored.
#'
#' Centres are always geometric (log scale), since the offset enters as
#' `log(mean_rate)`; `opt$size$mean_rate$estimator` changes only the reported
#' `mr_raw` diagnostic. Fallback ladder (tagged in `mr_tier`):
#'   T1 own (country, level) cell, n >= 1: EB-shrunk toward the country's
#'      level-adjusted parent, `w = n/(n+k)`, `k = sigma2_within / tau2_between`
#'      estimated from the data.
#'   T2 country present, level absent: country parent + `L(level)`.
#'   T3/T4/T5: resolved at join time by `resolve_mean_rate()`.
#' `L(level)` is the country-demeaned admin-level contrast (ADM1 = 0); a level
#' whose 95% CI crosses 0 is set to 0.
#'
#' @param opt options list (`opt$size$min_pop`, `opt$size$max_attack_rate`,
#'   `opt$size$mean_rate$estimator`; `opt$country` ignored, see above).
#' @param path path to the raw outbreak-size analysis data (same file
#'   `read_size_data()` reads).
#' @param enso_grouping_path path to the iso3 -> ENSO-group registry, used
#'   for the T3 (group) fallback tier and baked into `$meta$country_group`.
#' @return list(cell, group, global, meta):
#'   `cell`   tibble(country, admin_level, mean_rate, mr_tier, n_cell,
#'            w_shrink, mr_raw, mr_parent, estimator, k, level_factor) --
#'            one row per (country, admin_level) for every country with
#'            >=1 training row, crossed with all 3 levels (T1 where the cell
#'            has data, T2 otherwise).
#'   `group`  tibble(enso_group, admin_level, mean_rate, mr_tier, n_group)
#'            for the T3 fallback (country absent entirely from training).
#'   `global` tibble(admin_level, mean_rate, mr_tier) for the T4 fallback.
#'   `meta`   list(estimator, k, level_factor, level_factor_ci, n_rows,
#'            n_cells, n_countries, country_group, build_time).
build_mean_rate_table <- function(
  opt,
  path = make_size_dat_file(opt),
  enso_grouping_path = "analysis/groupings/enso_grouping.csv"
) {
  if (!file.exists(path)) {
    stop("build_mean_rate_table: missing input ", path)
  }
  d <- readRDS(path)

  # Capture the population stamp IMMEDIATELY, before the quality guards below.
  # dplyr's attribute handling is verb-dependent -- filter/mutate/select/
  # arrange/left_join currently carry non-standard attributes through, but
  # summarise() constructs a fresh frame and drops them (verified against dplyr
  # 1.2.1 in analysis/tests/test_pop_provenance_guard.R). Reading the stamp
  # here rather than downstream means this function does not depend on which
  # side of that line any particular verb happens to fall on today.
  pop_prov_in <- get_pop_provenance(d)

  min_pop <- opt$size$min_pop %||% 10
  max_ar <- opt$size$max_attack_rate %||% 0.10
  estimator <- size_mean_rate_estimator(opt) # shared resolver, utils/paths.R
  if (!estimator %in% c("geometric", "arithmetic", "pooled")) {
    stop(
      "build_mean_rate_table: unknown estimator '",
      estimator,
      "' (must be geometric | arithmetic | pooled)"
    )
  }

  # --- build the rate on the DECLARED denominator ---------------------------
  # This is the single most consequential line in the whole coherence problem.
  # `mean_rate` enters the model as `offset(log(mean_rate))`, so it carries a
  # unit: cases per person of WHICHEVER population it was divided by. Pairing a
  # rate built on a reporting footprint (H) with an offset built on whole-gid
  # population (Q) multiplies every prediction by Q/H -- the coverage ratio --
  # and because both offsets are fixed, nothing in the likelihood can absorb it.
  # Making the denominator an explicit, stamped choice is what turns that from
  # an invisible error into a checkable one.
  estimand <- resolve_size_estimand(opt)
  d <- resolve_size_pop_column(d, estimand, tag = "build_mean_rate_table")

  # --- canonical admin level -------------------------------------------------
  # `admin_level` is this function's POOLING KEY, so a mis-keyed row is averaged
  # into the wrong cell and then handed back to every consumer of that cell.
  # This is also why the phantom-ADM3 problem matters here specifically: the
  # canonical table has no ADM3 units at all, so an "ADM3" cell is a mean over
  # rows that are canonically ADM1 or ADM2.
  use_canonical_admin_level <- canonical_admin_level_enabled(
    opt,
    "size",
    estimand
  )
  d <- .apply_p3(d, opt, "size", estimand, tag = "build_mean_rate_table/P3")

  # Same quality guards as read_size_data() -- deliberately NOT the
  # admin_level/country subsetting; see the deviation note above.
  d <- d %>%
    dplyr::filter(!is.na(TL), !is.na(TR)) %>%
    dplyr::filter(is.na(pop) | pop >= min_pop) %>%
    dplyr::filter(is.na(pop) | !(pop > 0) | suspected_cases <= max_ar * pop)

  n_before <- nrow(d)
  d <- d %>%
    dplyr::filter(
      is.finite(pop),
      pop > 0,
      is.finite(suspected_cases),
      suspected_cases > 0
    )
  n_dropped <- n_before - nrow(d)
  if (n_dropped > 0) {
    warning(
      "build_mean_rate_table: dropped ",
      n_dropped,
      " row(s) with ",
      "non-positive pop/suspected_cases -- log(rate) undefined."
    )
  }
  if (!nrow(d)) {
    stop("build_mean_rate_table: no rows survive the quality guards")
  }

  d <- d %>%
    dplyr::mutate(
      country = as.character(country),
      admin_level = as.character(admin_level),
      rate = suspected_cases / pop,
      log_rate = log(rate)
    )

  allowed_levels <- c("ADM1", "ADM2", "ADM3")
  if (!all(d$admin_level %in% allowed_levels)) {
    bad <- setdiff(unique(d$admin_level), allowed_levels)
    stop(
      "build_mean_rate_table: unexpected admin_level value(s): ",
      paste(bad, collapse = ", ")
    )
  }

  # Which levels the tables are CROSSED over. Crossing a level that no row
  # occupies manufactures T2/T3/T4 cells for it, and every one of those is a
  # pure extrapolation from `level_factor` -- a number with no data behind it.
  # That is how the phantom ADM3 stratum acquired `mean_rate` values.
  #
  # But suppressing it is NOT free, and NOT a no-op even under corpus labels:
  # the crossed rows are also the fallback that lets resolve_mean_rate() answer a
  # query at a level absent from training (level-absent -> T2/T3/T4 instead of
  # T5). test_mean_rate_table.R exercises exactly that path with an ADM1/ADM2
  # fixture and an ADM3 query, and suppression turns it into a hard T5 error.
  # So this is a switch, defaulting to the existing crossing behaviour; set it
  # FALSE alongside `canonical_admin_level` where the absent level is phantom
  # (canonical keys have no ADM3 at all) rather than merely unobserved.
  cross_absent <- opt$size$mean_rate$cross_absent_levels %||% TRUE
  valid_levels <- if (isTRUE(cross_absent)) {
    allowed_levels
  } else {
    intersect(allowed_levels, unique(d$admin_level))
  }
  if (!identical(valid_levels, allowed_levels)) {
    message(sprintf(
      paste0(
        "[build_mean_rate_table] levels present: %s (of %s) -- absent level(s) are ",
        "NOT crossed into the cell/group/global tables ",
        "(size.mean_rate.cross_absent_levels = FALSE)."
      ),
      paste(valid_levels, collapse = "/"),
      paste(allowed_levels, collapse = "/")
    ))
  }

  # --- k: EB shrinkage constant, sigma2_within / tau2_between, log scale ----
  cell_var <- d %>%
    dplyr::group_by(country, admin_level) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_log = mean(log_rate),
      var_log = if (dplyr::n() > 1) stats::var(log_rate) else NA_real_,
      .groups = "drop"
    )
  within_ok <- !is.na(cell_var$var_log)
  sigma2_within <- if (any(within_ok)) {
    stats::weighted.mean(
      cell_var$var_log[within_ok],
      w = cell_var$n[within_ok] - 1
    )
  } else {
    NA_real_
  }
  tau2_between <- stats::var(cell_var$mean_log)
  if (!is.finite(sigma2_within) || sigma2_within <= 0) {
    sigma2_within <- tau2_between
  }
  k <- sigma2_within / tau2_between
  if (!is.finite(k) || k <= 0) {
    k <- 1
  } # degenerate corpus: fall back to equal-weight shrink

  # --- L(level): country-demeaned admin-level contrast vs ADM1, CI-gated ----
  #
  # `level_factor` is a fallback and borrowing device, not a model term: it
  # scales a rate borrowed across admin levels (level-adjusting rows before
  # forming parents, re-inflating borrowed parents at T2-T4). On the deployed
  # frame every gid resolves at T1, so dropping it is near-free.
  # WARNING: that holds only while the frame equals the outbreak corpus. If the
  # frame grows to gids with no outbreak history, T2-T4 fire and level_factor is
  # what adjusts a borrowed ADM1 rate: restore it or replace it with an in-model
  # level term. Default stays TRUE.
  use_level_factor <- opt$size$mean_rate$level_factor %||% TRUE

  multi_level <- d %>%
    dplyr::group_by(country) %>%
    dplyr::filter(dplyr::n_distinct(admin_level) > 1) %>%
    dplyr::ungroup()
  level_factor <- stats::setNames(
    rep(0, length(allowed_levels)),
    allowed_levels
  )
  level_factor_ci <- NULL
  if (!isTRUE(use_level_factor)) {
    message(
      "[build_mean_rate_table] level_factor DISABLED (D2): all levels 0; ",
      "parents are unadjusted geometric means and borrowed rates are not ",
      "re-inflated across levels."
    )
  } else if (
    nrow(multi_level) > 0 && dplyr::n_distinct(multi_level$country) >= 2
  ) {
    dm <- multi_level %>%
      dplyr::group_by(country) %>%
      dplyr::mutate(log_rate_dm = log_rate - mean(log_rate)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(admin_level = factor(admin_level, levels = valid_levels))
    lvl_fit <- stats::lm(log_rate_dm ~ admin_level, data = dm)
    lvl_ci <- stats::confint(lvl_fit)
    for (lv in c("ADM2", "ADM3")) {
      coef_name <- paste0("admin_level", lv)
      if (coef_name %in% names(stats::coef(lvl_fit))) {
        ci <- lvl_ci[coef_name, ]
        if (ci[1] <= 0 && ci[2] >= 0) {
          warning(
            "build_mean_rate_table: level factor for ",
            lv,
            " CI crosses 0 (",
            round(ci[1], 3),
            ", ",
            round(ci[2], 3),
            ") -- gated to 0 per D3."
          )
        } else {
          level_factor[[lv]] <- stats::coef(lvl_fit)[[coef_name]]
        }
      }
    }
    level_factor_ci <- lvl_ci
  } else {
    warning(
      "build_mean_rate_table: fewer than 2 countries span multiple admin ",
      "levels -- level_factor left at 0 for all levels."
    )
  }

  # Level-adjust every row to its ADM1-equivalent, for country/group/global parents.
  d <- d %>% dplyr::mutate(log_rate_adj = log_rate - level_factor[admin_level])

  # --- country parent: level-adjusted geometric mean, one per country -------
  country_tbl <- d %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      log_r_country = mean(log_rate_adj),
      n_country = dplyr::n(),
      .groups = "drop"
    )

  # --- ENSO group registry + parents (T3) ------------------------------------
  if (!file.exists(enso_grouping_path)) {
    stop(
      "build_mean_rate_table: ENSO grouping registry not found: ",
      enso_grouping_path
    )
  }
  grouping <- utils::read.csv(enso_grouping_path, stringsAsFactors = FALSE) %>%
    dplyr::transmute(
      country = as.character(iso3),
      enso_group = as.character(group)
    )
  d_grp <- d %>% dplyr::left_join(grouping, by = "country")
  if (any(is.na(d_grp$enso_group))) {
    missing_ctry <- unique(d_grp$country[is.na(d_grp$enso_group)])
    warning(
      "build_mean_rate_table: country/countries not in ENSO grouping registry ",
      "(T3 unavailable for them if they later drop out of the cell table): ",
      paste(missing_ctry, collapse = ", ")
    )
  }
  group_tbl <- d_grp %>%
    dplyr::filter(!is.na(enso_group)) %>%
    dplyr::group_by(enso_group) %>%
    dplyr::summarise(
      log_r_group = mean(log_rate_adj),
      n_group = dplyr::n(),
      .groups = "drop"
    )
  group_level_tbl <- d_grp %>%
    dplyr::filter(!is.na(enso_group)) %>%
    dplyr::group_by(enso_group, admin_level) %>%
    dplyr::summarise(log_r = mean(log_rate), n = dplyr::n(), .groups = "drop")

  # --- global parent (T4) -----------------------------------------------------
  log_r_global <- mean(d$log_rate_adj)

  # --- T1/T2 cell table: every (country, level) for every country present ----
  own_cell <- d %>%
    dplyr::group_by(country, admin_level) %>%
    dplyr::summarise(
      n_cell = dplyr::n(),
      log_r_raw_geo = mean(log_rate),
      mr_raw = switch(
        estimator,
        geometric = exp(mean(log_rate)),
        arithmetic = mean(rate),
        pooled = sum(suspected_cases) / sum(pop)
      ),
      .groups = "drop"
    )

  cell <- tidyr::crossing(
    country = unique(d$country),
    admin_level = valid_levels
  ) %>%
    dplyr::left_join(own_cell, by = c("country", "admin_level")) %>%
    dplyr::left_join(country_tbl, by = "country") %>%
    dplyr::mutate(
      lvl_fac = level_factor[admin_level],
      log_r_parent = log_r_country + lvl_fac,
      w_shrink = ifelse(is.na(n_cell), NA_real_, n_cell / (n_cell + k)),
      mr_tier = ifelse(is.na(n_cell), "T2", "T1"),
      log_mean_rate = ifelse(
        is.na(n_cell),
        log_r_parent,
        w_shrink * log_r_raw_geo + (1 - w_shrink) * log_r_parent
      ),
      mean_rate = exp(log_mean_rate),
      mr_parent = exp(log_r_parent),
      estimator = estimator,
      k = k,
      level_factor = lvl_fac
    ) %>%
    dplyr::select(
      country,
      admin_level,
      mean_rate,
      mr_tier,
      n_cell,
      w_shrink,
      mr_raw,
      mr_parent,
      estimator,
      k,
      level_factor
    )

  group <- tidyr::crossing(
    enso_group = unique(group_tbl$enso_group),
    admin_level = valid_levels
  ) %>%
    dplyr::left_join(group_level_tbl, by = c("enso_group", "admin_level")) %>%
    dplyr::left_join(group_tbl, by = "enso_group") %>%
    dplyr::mutate(
      lvl_fac = level_factor[admin_level],
      log_mean_rate = ifelse(!is.na(log_r), log_r, log_r_group + lvl_fac),
      mean_rate = exp(log_mean_rate),
      mr_tier = ifelse(!is.na(log_r), "T3-own", "T3-parent")
    ) %>%
    dplyr::select(enso_group, admin_level, mean_rate, mr_tier, n_group)

  global <- tibble::tibble(
    admin_level = valid_levels,
    mean_rate = exp(log_r_global + level_factor[valid_levels]),
    mr_tier = "T4"
  )

  # Every cell above is a function of `pop`. Record WHICH populations produced
  # them, so a later consumer can refuse a table built from different ones.
  # `pop_provenance` is stamped onto the size data by 03_prepare_siz.R; taking
  # it from the input rather than recomputing here means the stamp describes
  # the FULL dataset, not the quality-guard-filtered subset used above.
  pop_prov <- pop_prov_in %||%
    pop_provenance_stamp(d$pop, pop_source = "unstamped")

  mrt <- list(
    cell = cell,
    group = group,
    global = global,
    meta = list(
      estimator = estimator,
      k = k,
      level_factor = level_factor,
      level_factor_ci = level_factor_ci,
      n_rows = nrow(d),
      n_cells = nrow(own_cell),
      n_countries = length(unique(d$country)),
      country_group = grouping,
      build_time = Sys.time(),
      pop_provenance = pop_prov,
      # structural, not an attribute -- `summarise()` drops
      # attributes, and this table is the side of the seam most
      # likely to be reached via a rebuilt object.
      size_estimand = build_size_estimand_stamp(
        opt,
        extra = list(
          use_level_factor = isTRUE(use_level_factor),
          levels_present = paste(valid_levels, collapse = "/"),
          # Which admin_level column this table's POOLING KEY is --
          # `build_location_frame()`'s replay frame carries both the
          # corpus label and the canonical one, and
          # `.resolve_frame_anchor()` (utils/simulator_core.R) reads
          # this field to join on the same convention the pooling
          # used. Without it the simulator has to re-derive the rule.
          canonical_admin_level = use_canonical_admin_level,
          n_rows = nrow(d)
        )
      )
    )
  )
  attach_pop_provenance(mrt, pop_prov)
}

#' resolve_mean_rate
#'
#' The ONLY lookup path a consumer should ever use to obtain a `mean_rate`
#' from a table built by `build_mean_rate_table()`. Fit time (`read_size_data()`)
#' and simulate time (`build_location_frame()`) both call this same function
#' on the same persisted table, so the two sides cannot drift. Resolves
#' T1/T2 by exact (country, admin_level) match against `mrt$cell`; falls
#' through to T3 (`mrt$group`, via `mrt$meta$country_group`) and T4
#' (`mrt$global`) for countries entirely absent from the training data; T5
#' (anything else -- unmapped country AND unmapped group, or an
#' unrecognised `admin_level`) is a hard `stop()`, never a silent default.
#'
#' @param country character vector of ISO3 codes.
#' @param admin_level character vector, same length as `country`.
#' @param mrt list returned by `build_mean_rate_table()` (or read back from its
#'   persisted sidecar via `make_size_run_file(opt, "size_mean_rate")`).
#' @return tibble(mean_rate, mr_tier), one row per input pair, in input order.
resolve_mean_rate <- function(country, admin_level, mrt) {
  stopifnot(length(country) == length(admin_level))
  key <- tibble::tibble(
    .row = seq_along(country),
    country = as.character(country),
    admin_level = as.character(admin_level)
  )
  out <- key %>% dplyr::left_join(mrt$cell, by = c("country", "admin_level"))

  need_fallback <- is.na(out$mean_rate)
  if (any(need_fallback)) {
    fb <- out[need_fallback, c(".row", "country", "admin_level")] %>%
      dplyr::left_join(mrt$meta$country_group, by = "country") %>%
      dplyr::left_join(mrt$group, by = c("enso_group", "admin_level")) %>%
      dplyr::left_join(
        mrt$global,
        by = "admin_level",
        suffix = c("", ".global")
      ) %>%
      dplyr::mutate(
        mean_rate = dplyr::coalesce(mean_rate, mean_rate.global),
        mr_tier = dplyr::coalesce(mr_tier, mr_tier.global)
      )
    still_na <- is.na(fb$mean_rate)
    if (any(still_na)) {
      stop(
        "resolve_mean_rate: T5 -- no cell/group/global mean_rate resolves for: ",
        paste(
          unique(paste(
            fb$country[still_na],
            fb$admin_level[still_na],
            sep = "/"
          )),
          collapse = ", "
        )
      )
    }
    out$mean_rate[need_fallback] <- fb$mean_rate
    out$mr_tier[need_fallback] <- fb$mr_tier
  }
  out[order(out$.row), c("mean_rate", "mr_tier")]
}

#' persist_mean_rate_table
#'
#' Write the `build_mean_rate_table()` result carried as an attribute on
#' `read_size_data()`'s return value to the shared, run-invariant sidecar path
#' (`make_size_run_file(opt, "size_mean_rate")` -- deliberately excludes the
#' per-run tokens country/admin_level/family/... so every 04_fit_siz_*.R
#' invocation resolves to the same file and none of them recomputes a second
#' copy), and assert that `data$mean_rate` reconciles exactly against a fresh
#' lookup in the ON-DISK sidecar -- catching any
#' drift between what this run computed in-memory and what is persisted for
#' every other consumer. Call once, immediately after `read_size_data()`, from
#' each of the three size fit scripts.
#'
#' The sidecar path is content-addressed by population (see
#' `current_pop_digest_token()`), so the `if (!file.exists(path))` skip below is
#' safe: a given path's content is immutable (same min_pop/max_attack_rate/
#' estimator/pop_digest quadruple -> same table), and a population change
#' resolves to a DIFFERENT path, so `file.exists()` misses and the table
#' rebuilds automatically instead of silently reusing a stale one.
#'
#' @param opt options list.
#' @param data the tibble returned by `read_size_data()` (must carry the
#'   `mean_rate_table` attribute).
#' @return invisibly, the sidecar path.
persist_mean_rate_table <- function(opt, data) {
  mrt <- attr(data, "mean_rate_table")
  if (is.null(mrt)) {
    stop(
      "persist_mean_rate_table: `data` carries no mean_rate_table attribute ",
      "-- was it produced by read_size_data()?"
    )
  }
  path <- make_size_run_file(opt, kind = "size_mean_rate")
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(mrt, path)
    message("size_mean_rate (global table) -> ", path)
  }
  on_disk <- readRDS(path)

  # Provenance BEFORE arithmetic. The numeric reconciliation below does detect a
  # stale sidecar -- different populations give different rates -- but it
  # reports the symptom ("mean_rate does not match") and leaves the operator to
  # infer the cause. Checking the stamp first names the cause directly, and
  # names the rebuild order that fixes it.
  assert_pop_provenance(
    data_stamp = pop_provenance_of(mrt),
    sidecar_stamp = pop_provenance_of(on_disk),
    what = paste0("the mean_rate sidecar at ", path)
  )

  check <- resolve_mean_rate(
    data$country,
    as.character(data$admin_level),
    on_disk
  )
  if (!isTRUE(all.equal(check$mean_rate, data$mean_rate, tolerance = 1e-8))) {
    stop(
      "persist_mean_rate_table: reconciliation FAILED -- data$mean_rate does not ",
      "match a fresh lookup against the persisted sidecar at ",
      path,
      ". The ",
      "in-memory table used to build this fit's data has drifted from what is ",
      "on disk; delete the stale sidecar and re-run.\n",
      "If the population definition changed (the crosswalk wiring), this is ",
      "the expected symptom: every cell of the table is suspected_cases / ",
      "pop. Required order: pop fix -> mean_rate rebuild -> size refit."
    )
  }
  invisible(path)
}

# ===========================================================================
# Generated-config files
# ===========================================================================
#
# One writer for every config a tool generates (forecast / hindcast configs,
# deploy sim templates): a "# GENERATED by ... DO NOT EDIT BY HAND." banner,
# caller-specific comment lines, a blank line, then the YAML body.

#' Header lines for a generated YAML file.
#'
#' @param generator What wrote the file (script path, optionally with a stage).
#' @param body      Further comment lines (each already starting with "#").
#' @param when      Timestamp to stamp.
#' @return character vector: banner, body, "".
generated_yaml_header <- function(
  generator,
  body = character(0),
  when = Sys.time()
) {
  c(
    sprintf("# GENERATED by %s on %s — DO NOT EDIT BY HAND.", generator, when),
    body,
    ""
  )
}

#' Write a config list as YAML under a header (see generated_yaml_header()).
write_yaml_with_header <- function(cfg, path, header_lines) {
  con <- file(path, open = "w")
  on.exit(close(con))
  writeLines(header_lines, con)
  writeLines(yaml::as.yaml(cfg), con)
  invisible(path)
}

#' manifest_add_provenance
#'
#' Add a provenance field to a generated manifest only when it is not the
#' default: `manifest[[name]] <- value` when `!identical(value, default)`,
#' else `manifest` unchanged. A NULL `value` never adds (or removes) anything.
#' The rule that lets a new setting be recorded without moving any existing
#' manifest: a run at the default regenerates byte-identically. Lives
#' here because all four manifest generators (utils/make_forecast_configs.R,
#' make_hindcast_configs.R, make_event_attribution_configs.R,
#' make_forecast_attrib_configs.R) source provenance.R. Users include
#' `enso_counterfactual`, `pop_year`, `firm_oni_end` and `firm_horizon_rule`.
#' `identical()` is strict: pass `value`
#' and `default` with the same type (e.g. both character).
#' @return the manifest list.
manifest_add_provenance <- function(manifest, name, value, default) {
  if (is.null(value) || identical(value, default)) {
    return(manifest)
  }
  manifest[[name]] <- value
  manifest
}

# ---------------------------------------------------------------------------
# Processing-step counts (supplementary table S1, the processing cascade)
# ---------------------------------------------------------------------------
# Steps 00-02 write record / location-period counts to
# generated_data/prep_counts/<step>.yml for 09's report_dataset cascade table
# (03 keeps its counts in the .policy.rds stamps). Each file also records the
# outbreak-database `indir` and `data_vintage`: metadata that
# read_step_counts() leaves out and read_step_provenance() returns (NA in
# older files).

#' Keys of a count file that are metadata, not counts (reserved in `counts`).
.STEP_COUNTS_META <- c("step", "written", "indir", "data_vintage")
.STEP_COUNTS_PROVENANCE <- c("indir", "data_vintage")

#' Path of one step's count file.
step_counts_path <- function(step, dir = "generated_data/prep_counts") {
  file.path(dir, paste0(step, ".yml"))
}

#' Write a named list of scalar counts for `step` (overwrites).
#' @param step the script's file stem, e.g. "02_preprocess_crosswalk".
#' @param counts named list of length-one numbers (integers preferred); the
#'   names in `.STEP_COUNTS_META` are reserved.
#' @param indir the outbreak-data directory the step read (its resolved
#'   `--indir`), recorded as `indir` (tilde-expanded, normalised, no trailing
#'   slash) and `data_vintage` (its basename). NULL records neither.
write_step_counts <- function(
  step,
  counts,
  dir = "generated_data/prep_counts",
  indir = NULL
) {
  ok <- vapply(
    counts,
    function(x) length(x) == 1L && (is.numeric(x) || is.character(x)),
    logical(1)
  )
  if (!length(counts) || is.null(names(counts)) || !all(ok)) {
    stop(
      "write_step_counts: `counts` must be a named list of scalars.",
      call. = FALSE
    )
  }
  clash <- intersect(names(counts), .STEP_COUNTS_META)
  if (length(clash)) {
    stop(
      "write_step_counts: reserved name(s) in `counts`: ",
      paste(clash, collapse = ", "),
      call. = FALSE
    )
  }
  meta <- list(step = step, written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  if (!is.null(indir)) {
    if (
      !is.character(indir) ||
        length(indir) != 1L ||
        is.na(indir) ||
        !nzchar(indir)
    ) {
      stop(
        "write_step_counts: `indir` must be one non-empty path.",
        call. = FALSE
      )
    }
    meta$indir <- normalizePath(
      sub("/+$", "", path.expand(indir)),
      mustWork = FALSE
    )
    meta$data_vintage <- basename(meta$indir)
  }
  path <- step_counts_path(step, dir)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_yaml_with_header(c(meta, counts), path, generated_yaml_header(step))
  message("[", step, "] counts written: ", path)
  invisible(path)
}

#' Read one step's counts; NULL when the step has not written them.
#' The provenance keys (`indir`, `data_vintage`) are left out, so the result is
#' `step`, `written` and the counts whether or not the file records them (see
#' read_step_provenance()).
read_step_counts <- function(step, dir = "generated_data/prep_counts") {
  path <- step_counts_path(step, dir)
  if (!file.exists(path)) {
    return(NULL)
  }
  rec <- yaml::read_yaml(path)
  rec[setdiff(names(rec), .STEP_COUNTS_PROVENANCE)]
}

#' Provenance of one step's counts: list(step, written, indir, data_vintage).
#' NULL when the step has not written them; `indir` / `data_vintage` are NA in
#' a file written before they were recorded.
read_step_provenance <- function(step, dir = "generated_data/prep_counts") {
  path <- step_counts_path(step, dir)
  if (!file.exists(path)) {
    return(NULL)
  }
  rec <- yaml::read_yaml(path)
  stats::setNames(
    lapply(.STEP_COUNTS_META, function(k) {
      if (is.null(rec[[k]])) NA_character_ else rec[[k]]
    }),
    .STEP_COUNTS_META
  )
}

#' software_versions
#'
#' Versions of R and of the packages the analysis relies on, in the session
#' that calls it (table S13). A package that is not installed is reported NA.
#' @param cmdstan FALSE drops the CmdStan row (sim_software_versions()).
#' @return tibble(software, version)
software_versions <- function(
  pkgs = c(
    "brms",
    "cmdstanr",
    "rstan",
    "StanHeaders",
    "posterior",
    "loo",
    "dlnm",
    "mgcv",
    "spatPomp",
    "pomp",
    "dplyr",
    "qs2",
    "ggplot2",
    "sf",
    "exactextractr"
  ),
  cmdstan = TRUE
) {
  pv <- vapply(
    pkgs,
    function(p) {
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) {
        NA_character_
      })
    },
    character(1)
  )
  r_rows <- c(
    R = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform
  )
  if (isTRUE(cmdstan)) {
    r_rows <- c(
      r_rows,
      CmdStan = tryCatch(
        as.character(cmdstanr::cmdstan_version()),
        error = function(e) NA_character_
      )
    )
  }
  tibble::tibble(
    software = c(names(r_rows), pkgs),
    version = c(unname(r_rows), unname(pv))
  )
}

#' Packages of the forward simulation / smoothing environment recorded by
#' sim_software_versions(): the spatPomp kernel (spatPomp, pomp), the posterior
#' draws it reads (posterior, brms), the ENSO crossbasis (dlnm), the frame and
#' output layer (dplyr, qs2) and the replicate parallelism (furrr, future).
SIM_SOFTWARE_PKGS <- c(
  "spatPomp",
  "pomp",
  "posterior",
  "brms",
  "dlnm",
  "dplyr",
  "qs2",
  "furrr",
  "future"
)

#' sim_software_versions
#'
#' The simulation environment (table S13): R, platform and
#' SIM_SOFTWARE_PKGS as a named list (software -> version; NA when not
#' installed), written by 07b spatPomp into each run's partition_meta.yml and
#' by 08c into the smoothing cache meta.yml (informational: no reader checks it).
sim_software_versions <- function(pkgs = SIM_SOFTWARE_PKGS) {
  sv <- software_versions(pkgs, cmdstan = FALSE)
  stats::setNames(as.list(sv$version), sv$software)
}
