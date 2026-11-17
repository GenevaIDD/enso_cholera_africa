# 02_preprocess_crosswalk.R — match location periods (LPs) to GADM units via
#   the population crosswalk.
#
# Usage:
#   Rscript analysis/02_preprocess_crosswalk.R -C <cfg.yml>
#   Rscript analysis/02_preprocess_crosswalk.R -C <cfg.yml> --strict
#
# Details:
#   Alternative to the single-centroid match in 02_preprocess_climate.R: an
#   LP's footprint is every gid holding WorldPop mass H[l,g] > 0, split by
#   w = H / P_l, so a composite LP is measured, never parsed from its name.
#   Reads <indir>/outbreak_shapefiles.gpkg (attributes only), the population
#   crosswalk and gid population (.qs), optionally the WPP national totals.
#   Writes matched_ob_gids_crosswalk.rds (many rows per LP),
#   u_gids_crosswalk.rds and the population_qc/ match / unmatched CSVs.
#   --strict makes any unmatched LP fatal.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(sf)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  # The data leg's one vintage (utils/paths.R::default_outbreak_data_dir()), shared
  # with 00/01.
  make_option(
    "--indir",
    default = default_outbreak_data_dir(),
    type = "character",
    help = "Folder containing outbreak_shapefiles.gpkg [default: %default]"
  ),
  make_option(
    "--outdir_gen",
    default = "generated_data",
    type = "character",
    help = "Folder for the matched artefacts [default: %default]"
  ),
  make_option(
    "--strict",
    action = "store_true",
    default = FALSE,
    help = "treat any unmatched location period as fatal"
  )
)
opt <- make_options_from_config(option_list = option_list, model = NULL)

xwalk_path <- make_run_file(opt, kind = "pop_crosswalk")
gid_path <- make_run_file(opt, kind = "pop_gid")
wpp_path <- make_run_file(opt, kind = "pop_wpp")
qc_dir <- make_run_file(opt, kind = "pop_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outdir_gen, recursive = TRUE, showWarnings = FALSE)

# 3. The location-period roster (attributes only — no geometry is used) --------
# Match location periods (LPs) to GADM via the population crosswalk, not a
# centroid: a composite LP (a union of several units) has one centroid, which
# puts all its cases and population on one arbitrary gid. The crosswalk holds
# the WorldPop mass H[l,g] of each LP in each gid: footprint = {g : H > 0},
# split w[l,g] = H / P_l, and composite status is measured (gid_cover), never
# parsed from the LP name. No geometry is evaluated here (02a/02b did the
# intersection); `sf` only reads the outbreak layer's attributes.
shp_path <- file.path(opt$indir, "outbreak_shapefiles.gpkg")
if (!file.exists(shp_path)) {
  stop("[lp_match] outbreak shapefiles not found: ", shp_path)
}

ob_lp <- sf::st_read(shp_path, quiet = TRUE) %>%
  sf::st_drop_geometry() %>%
  as_tibble() %>%
  filter(admin_level != "ADM0") %>%
  # One row per LP: admin_level and country are LP properties, so the duplicate
  # rows a multi-part LP contributes carry identical values.
  distinct(lctn_pr, .keep_all = TRUE) %>%
  mutate(
    admin_level_simple = if_else(admin_level == "ADM1", "ADM1", "ADM2")
  ) %>%
  select(lctn_pr, admin_level, admin_level_simple, country)

message(sprintf(
  "[lp_match] %d location period(s) in %s (ADM0 excluded)",
  nrow(ob_lp),
  basename(shp_path)
))

# 4. The crosswalk -------------------------------------------------------------
xwalk <- load_population_crosswalk(xwalk_path)
gid_pop <- load_gid_population(gid_path)
wpp <- load_wpp_table(wpp_path)
if (is.null(wpp)) {
  message(
    "[lp_match] no WPP cache at ",
    wpp_path,
    " — T2 year-rescaling is unavailable (in-range years are unaffected)"
  )
}

req <- c("lctn_pr", "gid", "gadm_level", "year", "H")
if (length(setdiff(req, names(xwalk)))) {
  stop(
    "[lp_match] crosswalk is missing column(s): ",
    paste(setdiff(req, names(xwalk)), collapse = ", ")
  )
}

# A zero- or NA-mass piece is not a match. These are the sliver-filter residue;
# keeping them would put gids with no people into the climate universe and give
# them a w-share of the cases.
xwalk <- xwalk %>% filter(is.finite(H), H > 0)

message(sprintf(
  "[lp_match] crosswalk: %d piece(s), %d LP(s), %d gid(s), years %s",
  nrow(xwalk),
  n_distinct(xwalk$lctn_pr),
  n_distinct(xwalk$gid),
  paste(range(xwalk$year), collapse = "-")
))

# 5. Level selection, then the match -------------------------------------------
# The crosswalk is extracted at both tilings; an LP must be read at exactly one
# or its pieces are counted twice.
#
# The tiling comes from the LP's OWN GEOMETRY (`lp_level`, frozen by
# derive_lp_tiling_level() at 02c), not from the corpus's `admin_level` string.
# The corpus's admin depth is not GADM's: Malawi's corpus "ADM2" is a district,
# which GADM calls ADM1 (GADM's ADM2 is the Traditional Authority). Reading
# districts at the ADM2 tiling shattered each into 8-25 units and made
# `gid_cover` measure within-district concentration instead of compositeness.
# `admin_level` / `admin_level_simple` are retained as reported QC columns.
lvl_avail <- xwalk %>% distinct(lctn_pr, gadm_level)

if (!"lp_level" %in% names(xwalk)) {
  stop(
    "[lp_match] crosswalk has no `lp_level` column — rebuild it with ",
    "02c_preprocess_population_compile.R (derive_lp_tiling_level())."
  )
}
lp_derived <- xwalk %>%
  distinct(lctn_pr, lp_level) %>%
  filter(!is.na(lp_level))
if (anyDuplicated(lp_derived$lctn_pr)) {
  stop(
    "[lp_match] the crosswalk carries more than one `lp_level` for some LP(s) ",
    "— it must be one decision per location period."
  )
}

lp_lvl <- ob_lp %>%
  left_join(lvl_avail, by = "lctn_pr", relationship = "many-to-many") %>%
  group_by(lctn_pr) %>%
  summarise(
    have = list(sort(unique(gadm_level[!is.na(gadm_level)]))),
    any_level = if (all(is.na(gadm_level))) {
      NA_character_
    } else {
      sort(unique(gadm_level[!is.na(gadm_level)]))[1]
    },
    .groups = "drop"
  ) %>%
  left_join(ob_lp, by = "lctn_pr") %>%
  left_join(lp_derived, by = "lctn_pr") %>%
  # An LP with no derived level (absent from the crosswalk, or level undecidable)
  # keeps the corpus label — there is nothing better available for it.
  mutate(
    want_level = coalesce(lp_level, admin_level_simple),
    has_own = mapply(function(w, h) isTRUE(w %in% h), want_level, have),
    use_level = if_else(has_own, want_level, any_level)
  )

n_derived_diff <- sum(
  !is.na(lp_lvl$lp_level) &
    lp_lvl$lp_level != lp_lvl$admin_level_simple
)
message(
  sprintf(
    "[lp_match] tiling: %d of %d LP(s) read at a GADM level that ",
    n_derived_diff,
    nrow(lp_lvl)
  ),
  "disagrees with their corpus admin_level (the corpus label is not GADM depth)"
)
n_no_derived <- sum(is.na(lp_lvl$lp_level))
if (n_no_derived) {
  message(sprintf(
    "[lp_match] %d LP(s) have no derived level; corpus label used",
    n_no_derived
  ))
}
n_fallback <- sum(
  !is.na(lp_lvl$use_level) & lp_lvl$use_level != lp_lvl$want_level
)
if (n_fallback) {
  message(
    sprintf(
      "[lp_match] %d LP(s) read at a level other than the one chosen ",
      n_fallback
    ),
    "(that tiling was not extracted for them)"
  )
}

matched <- xwalk %>%
  inner_join(
    lp_lvl %>%
      select(lctn_pr, admin_level, admin_level_simple, country, use_level),
    by = "lctn_pr",
    relationship = "many-to-many"
  ) %>%
  filter(gadm_level == use_level) %>%
  select(-use_level)

# w is re-derived from H rather than carried, for the same reason
# apportion_composite_cases() re-derives it: the sliver filter above can drop
# pieces AFTER 02b computed w, and every conservation claim downstream rests on
# sum(w) == 1 exactly.
matched <- matched %>%
  group_by(lctn_pr, gadm_level, year) %>%
  mutate(
    P_l_kept = sum(H, na.rm = TRUE),
    w = if_else(P_l_kept > 0, H / P_l_kept, NA_real_),
    n_member = n()
  ) %>%
  ungroup()

# 6. Gates ---------------------------------------------------------------------
# (iv) shares sum to one. Not a tautology after the sliver filter + the level
# selection; if it fails, apportionment does not conserve cases.
bad_w <- matched %>%
  group_by(lctn_pr, gadm_level, year) %>%
  summarise(sw = sum(w, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.finite(sw) | abs(sw - 1) > 1e-8)
if (nrow(bad_w)) {
  stop(
    "[lp_match] ",
    nrow(bad_w),
    " (LP, level, year) group(s) have sum(w) != 1, ",
    "e.g. ",
    paste(head(bad_w$lctn_pr, 3), collapse = ", "),
    ". Apportionment would not conserve cases."
  )
}

# (iii) no piece may claim more people than the gid holds.
if ("c" %in% names(matched)) {
  over <- matched %>% filter(is.finite(c), c > 1 + 1e-6)
  if (nrow(over)) {
    warning(
      "[lp_match] ",
      nrow(over),
      " piece(s) with H > Q (c > 1) — ",
      "an LP claiming more people than its gid contains. Worst c = ",
      signif(max(over$c), 4)
    )
  }
}

unmatched <- ob_lp %>% filter(!lctn_pr %in% matched$lctn_pr)
if (nrow(unmatched)) {
  msg <- sprintf(
    "[lp_match] %d of %d location period(s) (%.1f%%) have no crosswalk piece",
    nrow(unmatched),
    nrow(ob_lp),
    100 * nrow(unmatched) / nrow(ob_lp)
  )
  if (isTRUE(opt$strict)) {
    stop(msg, " — see qc_lp_gid_unmatched.csv")
  }
  warning(
    msg,
    " — they keep their upstream population; see qc_lp_gid_unmatched.csv"
  )
}

# 7. P_l audit column ----------------------------------------------------------
# resolve_lp_population() is the sanctioned LP denominator and never returns NA.
# Recording its tier here means the fallback count is a reported QC number at
# match time rather than a surprise inside 03_prepare_siz.R.
# resolve_lp_population() reads one tiling per call, so iterate the levels rather
# than passing gadm_level[1] — a mixed-level roster would otherwise be resolved
# entirely at whichever level happened to sort first.
lp_years <- matched %>% distinct(lctn_pr, year, gadm_level)
lp_audit <- bind_rows(lapply(split(lp_years, lp_years$gadm_level), function(d) {
  resolve_lp_population(
    d$lctn_pr,
    d$year,
    xwalk,
    pop_upstream = NA_real_,
    gadm_level = d$gadm_level[1],
    wpp = wpp
  ) %>%
    transmute(
      lctn_pr,
      year,
      gadm_level = d$gadm_level,
      P_l = pop,
      pop_tier_lp = pop_tier
    )
}))

matched <- matched %>%
  left_join(lp_audit, by = c("lctn_pr", "year", "gadm_level")) %>%
  mutate(pop_tier = "crosswalk_H")

keep <- intersect(
  c(
    "lctn_pr",
    "admin_level",
    "admin_level_simple",
    "country",
    "gid",
    "gid_iso3",
    "gadm_level",
    "year",
    "H",
    "P_l",
    "P_l_kept",
    "Q",
    "w",
    "c",
    "gid_cover",
    "n_gid_eff",
    "is_composite",
    "is_composite_lp",
    "lp_level",
    "n_member",
    "reconcile_ok",
    "reconcile_rel_diff",
    "cross_border",
    "geom_dup_n",
    "geom_dup_group",
    "pop_tier",
    "pop_tier_lp"
  ),
  names(matched)
)
matched <- matched %>% select(all_of(keep)) %>% arrange(lctn_pr, year, desc(H))

# 8. Report + write ------------------------------------------------------------
u_gids_crosswalk <- sort(unique(matched$gid[!is.na(matched$gid)]))

# Reported over every matched LP-year, not only the ADM1-tiled ones: `matched`
# carries exactly one level per LP, so an ADM1 filter would silently omit
# every LP whose geometry put it at ADM2.
comp <- matched %>%
  distinct(
    lctn_pr,
    gadm_level,
    year,
    n_member,
    gid_cover,
    n_gid_eff,
    is_composite
  )
message(sprintf(
  "[lp_match] matched %d LP(s) -> %d piece(s), %d unique gid(s)",
  n_distinct(matched$lctn_pr),
  nrow(matched),
  length(u_gids_crosswalk)
))
message(sprintf(
  "[lp_match] composite at the selected tiling: %d of %d LP-year(s); max members = %d",
  sum(comp$is_composite, na.rm = TRUE),
  nrow(comp),
  suppressWarnings(max(comp$n_member, na.rm = TRUE))
))
if ("reconcile_ok" %in% names(matched)) {
  message(sprintf(
    "[lp_match] reconcile_ok on matched pieces: %.2f%%",
    100 * mean(matched$reconcile_ok, na.rm = TRUE)
  ))
}
if ("geom_dup_group" %in% names(matched)) {
  dup_lp <- matched %>%
    filter(!is.na(geom_dup_group)) %>%
    distinct(lctn_pr, gadm_level, geom_dup_group)
  message(sprintf(
    "[lp_match] duplicate-geometry LPs (see qc_geom_dup_lps.csv): %d of %d matched LP-level pair(s)",
    nrow(dup_lp),
    n_distinct(matched$lctn_pr, matched$gadm_level)
  ))
}
if ("pop_tier_lp" %in% names(matched)) {
  print(table(matched$pop_tier_lp, useNA = "ifany"))
}

# Processing cascade (table S1): location periods in, matched, composite.
write_step_counts(
  "02_preprocess_crosswalk",
  list(
    n_lp_in = nrow(ob_lp),
    n_lp_crosswalk = n_distinct(xwalk$lctn_pr),
    n_lp_matched = n_distinct(matched$lctn_pr),
    n_lp_unmatched = nrow(unmatched),
    n_countries_matched = n_distinct(substr(matched$gid, 1, 3)),
    n_gids_matched = length(u_gids_crosswalk),
    n_lp_year_matched = nrow(comp),
    n_lp_year_composite = sum(comp$is_composite, na.rm = TRUE)
  ),
  indir = opt$indir
)

saveRDS(matched, file.path(opt$outdir_gen, "matched_ob_gids_crosswalk.rds"))
saveRDS(u_gids_crosswalk, file.path(opt$outdir_gen, "u_gids_crosswalk.rds"))
write.csv(comp, file.path(qc_dir, "qc_lp_gid_match.csv"), row.names = FALSE)
write.csv(
  unmatched,
  file.path(qc_dir, "qc_lp_gid_unmatched.csv"),
  row.names = FALSE
)

message(
  "[lp_match] wrote ",
  file.path(opt$outdir_gen, "matched_ob_gids_crosswalk.rds"),
  " and ",
  file.path(opt$outdir_gen, "u_gids_crosswalk.rds")
)
