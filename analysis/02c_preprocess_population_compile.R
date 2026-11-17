# 02c_preprocess_population_compile.R — bind the population-extraction
#   partitions and derive the crosswalk shares over the whole corpus.
#
# Usage:
#   Rscript analysis/02c_preprocess_population_compile.R -C <cfg.yml>
#   Rscript analysis/02c_preprocess_population_compile.R -C <cfg.yml> \
#     --allow_incomplete
#
# Details:
#   Scans the 02b partitions for _SUCCESS sentinels, binds them, then runs the
#   derivations that cannot be done per task (a cross-border LP's shares need
#   the bound corpus). Writes population_crosswalk_<Y0>-<Y1>.qs,
#   gid_population_<Y0>-<Y1>.qs, wpp2024_national_totals.csv (if absent) and
#   the population_qc/ CSVs. Mirrors 08c_smoothing_compile.R.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(qs2)
  library(dplyr)
  library(tibble)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R") # flatten_gid() — gid key normalisation
source("analysis/utils/paths.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/population_core.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    "--allow_incomplete",
    action = "store_true",
    default = FALSE,
    help = "compile even if some planned tasks have no _SUCCESS"
  )
)
opt <- make_options_from_config(option_list = option_list, model = NULL)
pop <- opt$population

parts_dir <- make_run_file(opt, kind = "pop_parts_dir")
tasks_dir <- make_run_file(opt, kind = "pop_tasks_dir")
qc_dir <- make_run_file(opt, kind = "pop_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Collect successful partitions ---------------------------------------------
if (!dir.exists(parts_dir)) {
  stop("[pop_compile] parts dir not found: ", parts_dir)
}
sentinels <- partition_sentinels(parts_dir) # utils/paths.R
if (!length(sentinels)) {
  stop(
    "[pop_compile] no _SUCCESS sentinels under ",
    parts_dir,
    " — run 02b_preprocess_population_run.R first."
  )
}

n_planned <- length(list.files(tasks_dir, pattern = "^task_\\d+\\.yml$"))
message(sprintf(
  "[pop_compile] %d successful partition(s) of %d planned",
  length(sentinels),
  n_planned
))
if (n_planned > 0L && length(sentinels) < n_planned) {
  # A partial corpus does not fail visibly — it produces a crosswalk that is
  # simply missing countries, and every downstream lookup for them falls back.
  msg <- sprintf(
    "[pop_compile] %d planned task(s) have no _SUCCESS",
    n_planned - length(sentinels)
  )
  if (isTRUE(opt$allow_incomplete)) {
    warning(msg, " — compiling anyway")
  } else {
    stop(msg, ". Re-run with --resume, or pass --allow_incomplete.")
  }
}

read_part <- function(sent, fname) {
  p <- file.path(dirname(sent), fname)
  if (!file.exists(p)) {
    warning("[pop_compile] missing ", fname, " next to ", sent)
    return(NULL)
  }
  qs2::qs_read(p)
}

crosswalk <- bind_rows(lapply(sentinels, read_part, "crosswalk.qs"))
lp_pop <- bind_rows(lapply(sentinels, read_part, "lp_pop.qs"))
gid_pop <- bind_rows(lapply(sentinels, read_part, "gid_pop.qs"))
manifests <- bind_rows(lapply(sentinels, function(s) {
  m <- yaml::read_yaml(file.path(dirname(s), "manifest.yml"))
  as_tibble(lapply(m, function(x) if (is.null(x)) NA else x[1]))
}))

if (!nrow(crosswalk)) {
  stop("[pop_compile] bound crosswalk is empty.")
}
message(sprintf(
  "[pop_compile] bound: %d pieces, %d gid-years, %d LP-years",
  nrow(crosswalk),
  nrow(gid_pop),
  nrow(lp_pop)
))

# 3a. gid key normalisation (single convention boundary) -----------------------
# 02b extracts against GADM polygons, so parts carry GADM-native DOTTED ids
# ("NAM.1_1", "AGO.1.1_1"). Every outbreak-side consumer keys on the FLAT
# underscore form ("AGO_10_1"): both 03_prepare_* tracks overwrite gid with
# flatten_gid() (03_prepare_occ.R:148,172), and the deployed brms fits'
# RE level names are underscore-joined (occ_post.R:478).
#
# Leaving the two conventions to meet at read time breaks
# resolve_gid_population(): almost no survival gid matches the dotted keys, and
# the tier ladder then refuses to invent a denominator -- a naming mismatch
# presenting as total absence of population data.
#
# Normalising here, at the single point where the canonical artefacts are
# written, means every downstream reader sees one convention instead of each
# read seam having to re-derive it. Flattening is injective over this corpus
# (verified: 0 collisions across all 4861 gids), and the (gid, level, year)
# duplicate guard immediately below doubles as the standing assertion that
# it stays injective -- if flattening ever merged two distinct GADM units, the
# compile hard-fails here rather than silently summing them.
#
# `lp_pop` is keyed by lctn_pr and carries no gid column, so it is untouched.
# Applied BEFORE derive_crosswalk_shares() so every derived artefact (xwalk,
# census, tiling) inherits the flat keys.
crosswalk$gid <- flatten_gid(crosswalk$gid)
gid_pop$gid <- flatten_gid(gid_pop$gid)

# A gid belongs to exactly one country, so a duplicate (gid, level, year) means
# two tasks extracted the same polygon — the sum over gids would double-count.
dup_gid <- gid_pop %>% count(gid, gadm_level, year) %>% filter(n > 1L)
if (nrow(dup_gid)) {
  stop(
    "[pop_compile] ",
    nrow(dup_gid),
    " duplicated (gid, level, year) rows — ",
    "a stale parts/ tree from an earlier grid? First: ",
    paste(head(dup_gid$gid, 3), collapse = ", ")
  )
}

# 4. Derive shares over the BOUND corpus ---------------------------------------
# Why the derivation lives here and not in the worker: a cross-border LP's
# pieces arrive from two different (country, year) tasks. P_l = sum_g H is
# therefore only whole once the parts are bound, and every share (w,
# gid_cover, n_gid_eff) is a function of P_l. Deriving them per task would
# give each half of such an LP w = 1 and gid_cover = 1 — i.e. it would report
# the cross-border composites as clean.
cover_threshold <- pop$composite$cover_threshold %||% 0.90
xwalk <- derive_crosswalk_shares(
  crosswalk,
  gid_pop,
  cover_threshold = cover_threshold
)

# Extraction only ever runs inside WorldPop's constrained coverage. Years
# outside it are NOT materialised here: resolve_gid_population() rescales at
# lookup time by the WPP national ratio with these shares frozen, which is why
# the artefact carries a tier column rather than 1950-2100 of rows.
xwalk$pop_tier <- "direct"
gid_pop$pop_tier <- "direct"

out_of_range <- setdiff(
  unique(gid_pop$year),
  POP_RASTER_YEAR_MIN:POP_RASTER_YEAR_MAX
)
if (length(out_of_range)) {
  warning(
    "[pop_compile] extracted years outside WorldPop coverage [",
    POP_RASTER_YEAR_MIN,
    ", ",
    POP_RASTER_YEAR_MAX,
    "]: ",
    paste(sort(out_of_range), collapse = ", ")
  )
}

# 5. QC ------------------------------------------------------------------------
tol <- pop$reconcile_tol %||% 0.01
rec <- reconcile_crosswalk(crosswalk, lp_pop, tol = tol)
frac_ok <- mean(rec$ok, na.rm = TRUE)
message(
  sprintf(
    "[pop_compile] reconciliation: %.2f%% of %d (LP, level, year) ",
    100 * frac_ok,
    nrow(rec)
  ),
  sprintf("within %.1f%%", 100 * tol)
)
if (any(!rec$ok, na.rm = TRUE)) {
  worst <- rec %>% filter(!ok) %>% arrange(desc(rel_diff)) %>% head(10)
  message("[pop_compile] worst offenders:")
  print(as.data.frame(worst[, c(
    "lctn_pr",
    "gadm_level",
    "year",
    "P_l",
    "P_direct",
    "rel_diff"
  )]))
}

# Checkpoint A fails corpus-wide (94.90% vs the required 99% on the full
# extraction), concentrated in cross-border LPs at complex/FCV borders (GADM's
# ADM1-union boundary vs. WorldPop's own country raster mask disagree there).
# The gap is real but small in
# population-mass terms (~0.045% of total P_direct); rather than relax the
# gate or force P_l/P_direct into artificial agreement, we propagate the
# per-(LP, level, year) reconciliation result onto the persisted crosswalk so
# every downstream consumer (size model, simulator) can see, per row, whether
# that LP's gid-split is independently verified.
xwalk <- xwalk %>%
  dplyr::left_join(
    rec %>%
      dplyr::select(
        lctn_pr,
        gadm_level,
        year,
        reconcile_ok = ok,
        reconcile_rel_diff = rel_diff,
        cross_border
      ),
    by = c("lctn_pr", "gadm_level", "year")
  )

# Composite census — the whole point of the exercise, and it needs no string
# parsing: an LP is composite iff its mass spreads over more than one gid.
census <- compute_gid_cover(crosswalk, cover_threshold = cover_threshold)
comp <- census %>%
  filter(gadm_level == "ADM1", is_composite) %>%
  arrange(desc(n_gid_eff))
message(sprintf(
  "[pop_compile] composites at ADM1: %d of %d LP-years (%.1f%%)",
  nrow(comp),
  sum(census$gadm_level == "ADM1"),
  100 * nrow(comp) / max(sum(census$gadm_level == "ADM1"), 1L)
))
if (nrow(comp)) {
  message("[pop_compile] top composites by n_gid_eff:")
  print(as.data.frame(head(
    comp[, c("lctn_pr", "year", "n_gid", "gid_cover", "n_gid_eff")],
    10
  )))
}
message("[pop_compile] gid_cover distribution (ADM1):")
print(round(
  quantile(
    census$gid_cover[census$gadm_level == "ADM1"],
    probs = c(0, .05, .1, .25, .5, .75, .9, .95, 1),
    na.rm = TRUE
  ),
  3
))

# Tiling depth — the level each LP is actually READ at downstream. Frozen here
# rather than derived at run time so every consumer (both 03 tracks, all four
# --composite_policy values, 02_preprocess_crosswalk.R) shares one decision.
# The outbreak corpus's own `admin_level` string is not usable for this: its
# admin depth is not GADM's, which manufactures false composites (256 of 1410
# corpus LPs). See derive_lp_tiling_level().
#
# The area tie-break used WITHIN derive_lp_tiling_level() (when both levels are
# population-admissible) is itself configurable: "tolerant" (default) vs
# "strict" (the bare acov_ADM2 >= acov_ADM1 comparison, structurally biased
# toward ADM1).
tiebreak <- pop$composite$tiebreak %||% "tolerant" # = get_default_options()
tiebreak_floor <- pop$composite$tiebreak_floor %||% 0.80
tiebreak_tol <- pop$composite$tiebreak_tol %||% 0.10
tiling <- derive_lp_tiling_level(
  xwalk,
  cover_threshold = cover_threshold,
  tiebreak = tiebreak,
  tiebreak_floor = tiebreak_floor,
  tiebreak_tol = tiebreak_tol
)
message(sprintf(
  "[pop_compile] tiling level (tiebreak=%s): %d LP(s) — %s",
  tiebreak,
  nrow(tiling),
  paste(
    sprintf(
      "%s %d",
      names(table(tiling$lp_level)),
      as.integer(table(tiling$lp_level))
    ),
    collapse = ", "
  )
))
print(table(tiling$lp_level_reason))

xwalk <- xwalk %>%
  dplyr::mutate(lctn_pr = as.character(.data$lctn_pr)) %>%
  dplyr::left_join(
    tiling %>% dplyr::select("lctn_pr", "lp_level"),
    by = "lctn_pr"
  )

# `is_composite` stays per (LP, level, year) — it is the per-level fact and the
# QC census is built on it. `is_composite_lp` is the ONE flag that answers "is
# this LP a genuine multi-unit aggregate", evaluated at the LP's own selected
# tiling, so it cannot be read at the wrong level. Constant within an LP.
tiling <- tiling %>%
  dplyr::mutate(
    pcov_lp = dplyr::if_else(
      .data$lp_level == "ADM1",
      .data$pcov_ADM1,
      .data$pcov_ADM2
    ),
    is_composite_lp = !is.na(.data$pcov_lp) & .data$pcov_lp < cover_threshold
  )
xwalk <- xwalk %>%
  dplyr::left_join(
    tiling %>% dplyr::select("lctn_pr", "is_composite_lp"),
    by = "lctn_pr"
  )
message(sprintf(
  "[pop_compile] composite LPs at the selected tiling: %d of %d",
  sum(tiling$is_composite_lp, na.rm = TRUE),
  nrow(tiling)
))

# Overlapping sibling LPs: sum_l H[l,g] > Q_g means two LPs report the same
# people, so BOTH their cases and their population are double-counted in a
# collapsed (gid, period) cell. Reported, never silently fixed.
overlap <- xwalk %>%
  group_by(gid, gadm_level, year) %>%
  summarise(
    sum_H = sum(H, na.rm = TRUE),
    Q = dplyr::first(Q),
    n_lp = n_distinct(lctn_pr),
    .groups = "drop"
  ) %>%
  filter(!is.na(Q), Q > 0, sum_H > 1.05 * Q) %>%
  mutate(excess = sum_H / Q) %>%
  arrange(desc(excess))
message(sprintf(
  "[pop_compile] gid-years with sum_l H > 1.05 x Q: %d",
  nrow(overlap)
))
if (nrow(overlap)) {
  print(as.data.frame(head(overlap, 10)))
}

# Duplicate/parent-inherited geometry (downstream-detectable symptom):
# distinct LPs whose crosswalk pieces are byte-identical at a given GADM
# level. The real fix is a patched OutbreakExtractR re-extraction; this is
# the detect-and-quarantine flag made from data already on disk —
# reported and propagated onto the persisted crosswalk, never auto-corrected.
geom_dup <- detect_geom_dup_lps(crosswalk)
dup_clusters <- geom_dup %>%
  filter(!is.na(geom_dup_group)) %>%
  distinct(gadm_level, geom_dup_group, geom_dup_n)
message(sprintf(
  "[pop_compile] duplicate-geometry LPs: %d LP(s) in %d cluster(s) (max cluster size %d)",
  sum(geom_dup$geom_dup_n[!is.na(geom_dup$geom_dup_group)]),
  nrow(dup_clusters),
  suppressWarnings(max(dup_clusters$geom_dup_n, 0L))
))
if (nrow(dup_clusters)) {
  worst_dup <- geom_dup %>%
    filter(!is.na(geom_dup_group)) %>%
    arrange(desc(geom_dup_n), geom_dup_group)
  print(as.data.frame(head(worst_dup, 15)))
}

xwalk <- xwalk %>%
  dplyr::left_join(geom_dup, by = c("lctn_pr", "gadm_level"))

# 6. Persist -------------------------------------------------------------------
# Writes:
#   generated_data/population_crosswalk_<Y0>-<Y1>.qs   (carries reconcile_ok /
#     reconcile_rel_diff / cross_border per (lctn_pr, gadm_level, year) — see
#     "Checkpoint A" above — plus geom_dup_n / geom_dup_group per
#     (lctn_pr, gadm_level), the downstream-detected duplicate/parent-inherited
#     geometry flag: detect_geom_dup_lps(), analysis/utils/population_core.R —
#     plus lp_level / is_composite_lp per lctn_pr, the geometry-derived GADM
#     tiling depth and the composite flag evaluated at it:
#     derive_lp_tiling_level())
#   generated_data/gid_population_<Y0>-<Y1>.qs
#   generated_data/wpp2024_national_totals.csv        (refreshed if absent)
#   generated_data/population_qc/qc_crosswalk_reconciliation.csv
#   generated_data/population_qc/qc_composite_census.csv
#   generated_data/population_qc/qc_lp_tiling_level.csv
#   generated_data/population_qc/qc_gid_overlap.csv
#   generated_data/population_qc/qc_geom_dup_lps.csv
#   generated_data/population_qc/qc_manifest.csv
xwalk_path <- make_run_file(opt, kind = "pop_crosswalk")
gid_path <- make_run_file(opt, kind = "pop_gid")
wpp_path <- make_run_file(opt, kind = "pop_wpp")

dir.create(dirname(xwalk_path), recursive = TRUE, showWarnings = FALSE)
qs2::qs_save(xwalk, xwalk_path)
qs2::qs_save(gid_pop, gid_path)
tryCatch(write_wpp_national_totals(wpp_path), error = function(e) {
  warning("[pop_compile] WPP cache not refreshed: ", conditionMessage(e))
})

write.csv(
  rec,
  file.path(qc_dir, "qc_crosswalk_reconciliation.csv"),
  row.names = FALSE
)
write.csv(
  census,
  file.path(qc_dir, "qc_composite_census.csv"),
  row.names = FALSE
)
write.csv(
  tiling,
  file.path(qc_dir, "qc_lp_tiling_level.csv"),
  row.names = FALSE
)
write.csv(overlap, file.path(qc_dir, "qc_gid_overlap.csv"), row.names = FALSE)
write.csv(
  geom_dup %>% filter(!is.na(geom_dup_group)),
  file.path(qc_dir, "qc_geom_dup_lps.csv"),
  row.names = FALSE
)
write.csv(manifests, file.path(qc_dir, "qc_manifest.csv"), row.names = FALSE)

message("[pop_compile] crosswalk  -> ", xwalk_path)
message("[pop_compile] gid pop    -> ", gid_path)
message("[pop_compile] QC         -> ", qc_dir)

if (frac_ok < 0.99) {
  warning(
    sprintf(
      "[pop_compile] only %.2f%% of LPs reconcile within %.1f%% ",
      100 * frac_ok,
      100 * tol
    ),
    "— Checkpoint A of the plan requires >= 99%."
  )
}
