# population_qc.R — data-quality gates on the WorldPop population
#   denominator of the extracted outbreak data (read-only, diagnostic)
#
# Usage:
#   Rscript analysis/utils/population_qc.R
#   Rscript analysis/utils/population_qc.R --indir <dir> --outdir <dir>
#   Rscript analysis/utils/population_qc.R --ar_warn 100 --ar_impossible 1000
#
# Details:
#   Reads the two legacy files of 00_ingest_outbreakextractr.R:
#   time_series_outbreak_extraction.rds (outbreak-period rows; has `pop`) and
#   outbreak_shapefiles.rds (one geometry per location period). Writes a
#   diagnostic report + CSV; modifies no pipeline data. Also sourced: the gate
#   functions (run_population_gates(), check_*()) are exercised on synthetic
#   fixtures by tests/test_population_gates.R; the CLI block runs only when
#   the file is executed directly.

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(sf)
})

rule <- function(title) {
  cat(
    "\n",
    strrep("=", 70),
    "\n  ",
    title,
    "\n",
    strrep("=", 70),
    "\n",
    sep = ""
  )
}

# Admin level from a HopkinsIDD location string ("AFR::AGO" -> 0, "...::X::Y" -> n)
admin_level_from_loc <- function(loc) {
  if (requireNamespace("taxdat", quietly = TRUE)) {
    return(vapply(loc, function(x) taxdat::get_admin_level(x), numeric(1)))
  }
  str_count(loc, "::") - 1L # fallback
}


# Motivation: attack rates can come out nonsensical (max ~1.7e11 per 1,000),
# caused by bad `pop` values inherited from the upstream raster extraction:
# (1) a hard pop == 0 location period (cases / 0 -> Inf), and (2) ADM2
# health-district polygons whose WorldPop population is far too small for the
# observed case counts (attack rate > 100%, i.e. more cases than people).
# These are *data* problems, not a bug in the rate formula.

# A. Population-value integrity (per location period) ---------------------
#
# Gates 1-2: `pop` must be finite and non-NA (gate 1), and > 0 for
# every LP (gate 2). The known upstream failure mode is a tiny-but-positive
# pop (e.g. ~0.005), not an exact 0, so `flag_low_pop` uses a person-count
# floor rather than `pop == 0` alone.

check_population_integrity <- function(ob, shp, pop_floor = 1) {
  shp_area <- st_drop_geometry(shp) %>%
    select(lctn_pr, area_per_1km2) %>%
    mutate(lctn_pr = as.character(lctn_pr))

  ob %>%
    mutate(admin_level = admin_level_from_loc(location)) %>%
    group_by(country, location, location_period_id, admin_level) %>%
    summarise(
      pop_min = min(pop, na.rm = TRUE),
      pop_max = max(pop, na.rm = TRUE),
      n_rows = n(),
      na_pop = sum(is.na(pop)),
      .groups = "drop"
    ) %>%
    left_join(shp_area, by = c("location_period_id" = "lctn_pr")) %>%
    mutate(
      flag_low_pop = is.finite(pop_min) & pop_min < pop_floor,
      flag_na_pop = na_pop > 0,
      flag_no_geom = is.na(area_per_1km2)
    )
}


# B. Implausible attack rates (per outbreak event) ------------------------

check_attack_rates <- function(ob, shp, ar_warn = 100, ar_impossible = 1000) {
  shp_area <- st_drop_geometry(shp) %>%
    select(lctn_pr, area_per_1km2) %>%
    mutate(lctn_pr = as.character(lctn_pr))

  ob %>%
    filter(outbreak_number > 0) %>%
    mutate(admin_level = admin_level_from_loc(location)) %>%
    group_by(
      country,
      location,
      location_period_id,
      outbreak_number,
      admin_level
    ) %>%
    summarise(
      TL = min(TL),
      TR = max(TR),
      cases = sum(sCh, na.rm = TRUE),
      pop = mean(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(shp_area, by = c("location_period_id" = "lctn_pr")) %>%
    mutate(
      attack_rate_per_1k = ifelse(
        is.finite(pop) & pop > 0,
        cases / pop * 1e3,
        NA_real_
      ),
      flag_warn = is.finite(attack_rate_per_1k) & attack_rate_per_1k > ar_warn,
      flag_impossible = is.finite(attack_rate_per_1k) &
        attack_rate_per_1k > ar_impossible,
      flag_inf_rate = is.finite(pop) & pop <= 0 & cases > 0 # cases / 0
    )
}


# D. Duplicate-geometry provenance (gates 4 / 5a / 5b) ---------------------
#
# Consumes `pop_geom_dup_n` / `pop_geom_dup_class` if a patched OutbreakExtractR
# extraction produced them.
# `pop_geom_dup_class` taxonomy: unique | alias | parent_inherited | cross_unit |
# unknown. `pop_geom_dup_n == 1` should coincide with class == "unique"; any
# other class implies `pop_geom_dup_n > 1`. Returns NULL (with a message,
# not a warning/error) when the corpus predates these columns, so this gate
# never fails on pre-patch data.

check_duplicate_geometry_provenance <- function(ob) {
  has_cols <- all(c("pop_geom_dup_n", "pop_geom_dup_class") %in% names(ob))
  if (!has_cols) {
    message(
      "[population_qc] pop_geom_dup_n / pop_geom_dup_class absent — ",
      "corpus predates the OutbreakExtractR population-fix patch; ",
      "skipping Section D."
    )
    return(NULL)
  }

  lp_dup <- ob %>%
    distinct(
      country,
      location,
      location_period_id,
      pop_geom_dup_n,
      pop_geom_dup_class
    ) %>%
    mutate(
      known_class = pop_geom_dup_class %in%
        c("unique", "alias", "parent_inherited", "cross_unit", "unknown"),
      flag_unknown_class = !is.na(pop_geom_dup_class) & !known_class,
      flag_n_class_mismatch = case_when(
        is.na(pop_geom_dup_n) | is.na(pop_geom_dup_class) ~ NA,
        pop_geom_dup_n <= 1L ~ pop_geom_dup_class != "unique",
        TRUE ~ pop_geom_dup_class == "unique"
      )
    )

  lp_dup
}


# E. adj_factor banding (gate 3) -------------------------------------------
#
# Bands: [0.67, 1.5] accept silently;
# [0.5, 0.67) u (1.5, 2.0] accept but flag; outside [0.5, 2.0] flag as a
# clamp-to-1.0 candidate (the package itself is expected to clamp — this gate
# just reports what reached downstream so a clamp failure is visible).

check_adj_factor_banding <- function(ob, adj_col = "adj_factors") {
  if (!adj_col %in% names(ob)) {
    message("[population_qc] '", adj_col, "' absent — skipping Section E.")
    return(NULL)
  }

  ob %>%
    distinct(
      country,
      location,
      location_period_id,
      adj_factor = .data[[adj_col]]
    ) %>%
    mutate(
      band = case_when(
        is.na(adj_factor) ~ "missing",
        adj_factor >= 0.67 & adj_factor <= 1.5 ~ "accept",
        (adj_factor >= 0.5 & adj_factor < 0.67) |
          (adj_factor > 1.5 & adj_factor <= 2.0) ~ "accept_flag",
        TRUE ~ "out_of_band"
      ),
      flag_out_of_band = band == "out_of_band"
    )
}


# F. Crosswalk reconciliation (LP x gid population crosswalk) --------------
#
# P_l == sum_g H[l,g] is not a tautology once slivers can be dropped and
# cross-border pieces summed across tasks -- it is a real check that the
# vector-intersection crosswalk (analysis/utils/population_core.R) reproduces
# an INDEPENDENT whole-geometry extraction (P_direct, the `lp_whole` role in
# the same extraction call). Checkpoint A's pass condition is exactly this
# gate: >=99% of LPs reconcile within `tol`.
#
# Graceful-skip idiom matches check_duplicate_geometry_provenance(): returns
# NULL (with a message, not a warning/error) when the crosswalk / lp_pop
# artefacts haven't been built yet, so this never fails on a pre-Stage-1
# corpus.

check_crosswalk_reconciliation <- function(xwalk, lp_pop, tol = 0.01) {
  need_xw <- c("lctn_pr", "year", "H")
  need_lp <- c("lctn_pr", "year", "P_direct")
  has_cols <- all(need_xw %in% names(xwalk)) && all(need_lp %in% names(lp_pop))
  if (!has_cols) {
    message(
      "[population_qc] crosswalk/lp_pop reconciliation columns absent -- ",
      "Stage 1 population artefacts not built yet; skipping Section F."
    )
    return(NULL)
  }

  group_cols <- c("lctn_pr", "year")
  if ("gadm_level" %in% names(xwalk)) {
    group_cols <- c(group_cols, "gadm_level")
  }

  p_l_sum <- xwalk %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(P_l_sum = sum(H, na.rm = TRUE), .groups = "drop")

  p_l_sum %>%
    dplyr::left_join(
      lp_pop %>%
        dplyr::select(lctn_pr, year, P_direct, dplyr::any_of("cross_border")),
      by = c("lctn_pr", "year")
    ) %>%
    dplyr::mutate(
      rel_diff = ifelse(
        is.finite(P_direct) & P_direct > 0,
        abs(P_l_sum - P_direct) / P_direct,
        NA_real_
      ),
      flag_no_p_direct = is.na(P_direct),
      flag_reconcile_fail = dplyr::case_when(
        is.na(rel_diff) ~ NA,
        TRUE ~ rel_diff > tol
      )
    )
}


# Orchestration -------------------------------------------------------------

#' Run all population QC gates against a legacy-format outbreak table + shapefile.
#'
#' @param ob legacy `time_series_outbreak_extraction.rds`-shaped tibble.
#' @param shp legacy `outbreak_shapefiles.rds`-shaped sf object.
#' @param opt list with optional `ar_warn`, `ar_impossible`, `pop_floor`,
#'   `reconcile_tol`.
#' @param xwalk optional crosswalk table (from `load_population_crosswalk()`) for
#'   Section F. NULL skips Section F (pre-Stage-1 corpus).
#' @param lp_pop_direct optional per-LP `P_direct` table (Stage 1 `lp_pop.qs`)
#'   for Section F. NULL skips Section F. Named distinctly from the `lp_pop`
#'   list element below (Section A's per-LP integrity table) -- same-looking
#'   name, unrelated content.
#' @return list(lp_pop, evt, geom_dup, adj_band, reconcile) — `geom_dup`,
#'   `adj_band` and `reconcile` are NULL when their inputs/columns are absent.
run_population_gates <- function(
  ob,
  shp,
  opt = list(),
  xwalk = NULL,
  lp_pop_direct = NULL
) {
  ar_warn <- opt$ar_warn %||% 100
  ar_impossible <- opt$ar_impossible %||% 1000
  pop_floor <- opt$pop_floor %||% 1
  reconcile_tol <- opt$reconcile_tol %||% 0.01

  list(
    lp_pop = check_population_integrity(ob, shp, pop_floor = pop_floor),
    evt = check_attack_rates(
      ob,
      shp,
      ar_warn = ar_warn,
      ar_impossible = ar_impossible
    ),
    geom_dup = check_duplicate_geometry_provenance(ob),
    adj_band = check_adj_factor_banding(ob),
    reconcile = if (is.null(xwalk) || is.null(lp_pop_direct)) {
      NULL
    } else {
      check_crosswalk_reconciliation(xwalk, lp_pop_direct, tol = reconcile_tol)
    }
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# CLI block — only runs when this file is executed directly (Rscript), not
# when sourced by a caller such as analysis/tests/test_population_gates.R.
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(optparse))

  option_list <- list(
    make_option(
      "--indir",
      default = "generated_data/new_ingest_check",
      type = "character",
      help = "Folder with the legacy RDS files [default: %default]"
    ),
    make_option(
      "--outdir",
      default = "generated_data/population_qc",
      type = "character",
      help = "Folder for CSV QC reports [default: %default]"
    ),
    make_option(
      "--ar_warn",
      default = 100,
      type = "double",
      help = "Attack-rate warn threshold, cases per 1,000 (>10%% pop) [default: %default]"
    ),
    make_option(
      "--ar_impossible",
      default = 1000,
      type = "double",
      help = "Attack-rate impossible threshold, cases per 1,000 (>100%% pop) [default: %default]"
    ),
    make_option(
      "--pop_floor",
      default = 1,
      type = "double",
      help = "Flag LPs whose population is below this many people (broken extraction) [default: %default]"
    ),
    make_option(
      "--crosswalk",
      default = NULL,
      type = "character",
      help = "Optional path to population_crosswalk_<Y0>-<Y1>.qs (Stage 1). Enables Section F."
    ),
    make_option(
      "--lp_pop_direct",
      default = NULL,
      type = "character",
      help = "Optional path to a Stage 1 lp_pop.qs (P_direct table). Enables Section F."
    ),
    make_option(
      "--reconcile_tol",
      default = 0.01,
      type = "double",
      help = "Section F: relative-difference tolerance for |P_l_sum - P_direct| [default: %default]"
    )
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  ob_path <- file.path(opt$indir, "time_series_outbreak_extraction.rds")
  shp_path <- file.path(opt$indir, "outbreak_shapefiles.rds")
  stopifnot(file.exists(ob_path), file.exists(shp_path))

  ob <- readRDS(ob_path)
  shp <- readRDS(shp_path)

  xwalk_in <- if (!is.null(opt$crosswalk)) qs2::qs_read(opt$crosswalk) else NULL
  lp_pop_direct_in <- if (!is.null(opt$lp_pop_direct)) {
    qs2::qs_read(opt$lp_pop_direct)
  } else {
    NULL
  }

  gates <- run_population_gates(
    ob,
    shp,
    opt,
    xwalk = xwalk_in,
    lp_pop_direct = lp_pop_direct_in
  )
  lp_pop <- gates$lp_pop
  evt <- gates$evt

  rule("A. Population (pop) value integrity, per location period")
  cat("Distinct location periods:", nrow(lp_pop), "\n")
  cat(
    "  pop < ",
    opt$pop_floor,
    " person(s): ",
    sum(lp_pop$flag_low_pop),
    "\n",
    sep = ""
  )
  cat("  any NA pop              :", sum(lp_pop$flag_na_pop), "\n")
  cat("  no matching geometry    :", sum(lp_pop$flag_no_geom), "\n")

  bad_pop <- lp_pop %>% filter(flag_low_pop | flag_na_pop)
  if (nrow(bad_pop)) {
    cat("\nLPs with near-zero / NA population:\n")
    print(
      bad_pop %>%
        mutate(area_per_1km2 = round(area_per_1km2, 2)) %>%
        select(
          country,
          admin_level,
          location,
          location_period_id,
          pop_min,
          pop_max,
          na_pop,
          area_per_1km2
        ) %>%
        as.data.frame()
    )
  }

  rule("B. Implausible attack rates, per outbreak event")
  cat("Outbreak events:", nrow(evt), "\n")
  cat(
    "  rate > ",
    opt$ar_warn,
    "/1k (>",
    round(opt$ar_warn / 10, 1),
    "% pop): ",
    sum(evt$flag_warn),
    "\n",
    sep = ""
  )
  cat(
    "  rate > ",
    opt$ar_impossible,
    "/1k (>",
    round(opt$ar_impossible / 10, 1),
    "% pop, impossible): ",
    sum(evt$flag_impossible),
    "\n",
    sep = ""
  )
  cat(
    "  cases with pop <= 0 (infinite rate): ",
    sum(evt$flag_inf_rate),
    "\n",
    sep = ""
  )

  cat("\nFlagged events by country x admin level:\n")
  print(
    evt %>%
      filter(flag_warn | flag_inf_rate) %>%
      count(country, admin_level, sort = TRUE) %>%
      as.data.frame()
  )

  cat("\nWorst 20 events by attack rate:\n")
  print(
    evt %>%
      filter(flag_warn | flag_inf_rate) %>%
      arrange(desc(replace_na(attack_rate_per_1k, Inf))) %>%
      mutate(
        pop = round(pop, 1),
        attack_rate_per_1k = round(attack_rate_per_1k, 1),
        area_per_1km2 = round(area_per_1km2, 2)
      ) %>%
      select(
        country,
        admin_level,
        location,
        outbreak_number,
        cases,
        pop,
        attack_rate_per_1k,
        area_per_1km2
      ) %>%
      head(20) %>%
      as.data.frame()
  )

  rule("D. Duplicate-geometry provenance (pop_geom_dup_n / pop_geom_dup_class)")
  if (is.null(gates$geom_dup)) {
    cat("  [SKIP] corpus predates the OutbreakExtractR population-fix patch\n")
  } else {
    print(
      gates$geom_dup %>%
        count(pop_geom_dup_class, sort = TRUE) %>%
        as.data.frame()
    )
    n_mismatch <- sum(gates$geom_dup$flag_n_class_mismatch, na.rm = TRUE)
    cat("  pop_geom_dup_n / pop_geom_dup_class mismatches:", n_mismatch, "\n")
  }

  rule("E. adj_factor banding")
  if (is.null(gates$adj_band)) {
    cat("  [SKIP] no adj_factors column present\n")
  } else {
    print(gates$adj_band %>% count(band, sort = TRUE) %>% as.data.frame())
  }

  rule("F. Crosswalk reconciliation (P_l = sum_g H vs. independent P_direct)")
  if (is.null(gates$reconcile)) {
    cat(
      "  [SKIP] --crosswalk / --lp_pop_direct not supplied (Stage 1 artefacts not built yet)\n"
    )
  } else {
    n_tot <- nrow(gates$reconcile)
    n_fail <- sum(gates$reconcile$flag_reconcile_fail, na.rm = TRUE)
    n_nopd <- sum(gates$reconcile$flag_no_p_direct)
    cat(
      "  (lctn_pr, year",
      if ("gadm_level" %in% names(gates$reconcile)) ", gadm_level" else "",
      ") rows: ",
      n_tot,
      "\n",
      sep = ""
    )
    cat(
      "  reconcile within tol (",
      opt$reconcile_tol * 100,
      "%): ",
      n_tot - n_nopd - n_fail,
      "/",
      n_tot - n_nopd,
      "\n",
      sep = ""
    )
    cat("  no matching P_direct   : ", n_nopd, "\n", sep = "")
    if (n_fail) {
      cat("\nWorst 20 reconciliation failures:\n")
      print(
        gates$reconcile %>%
          filter(flag_reconcile_fail) %>%
          arrange(desc(rel_diff)) %>%
          mutate(
            P_l_sum = round(P_l_sum, 1),
            P_direct = round(P_direct, 1),
            rel_diff = round(rel_diff, 4)
          ) %>%
          head(20) %>%
          as.data.frame()
      )
    }
  }

  rule("Reports written")

  bad_pop_out <- file.path(opt$outdir, "qc_population_zero_na_lps.csv")
  write.csv(
    bad_pop %>%
      select(
        country,
        admin_level,
        location,
        location_period_id,
        pop_min,
        pop_max,
        na_pop,
        area_per_1km2
      ),
    bad_pop_out,
    row.names = FALSE
  )

  bad_evt_out <- file.path(opt$outdir, "qc_implausible_attack_rates.csv")
  write.csv(
    evt %>%
      filter(flag_warn | flag_inf_rate) %>%
      arrange(desc(replace_na(attack_rate_per_1k, Inf))) %>%
      select(
        country,
        admin_level,
        location,
        location_period_id,
        outbreak_number,
        TL,
        TR,
        cases,
        pop,
        attack_rate_per_1k,
        area_per_1km2,
        flag_warn,
        flag_impossible,
        flag_inf_rate
      ),
    bad_evt_out,
    row.names = FALSE
  )

  cat(normalizePath(bad_pop_out), "\n")
  cat(normalizePath(bad_evt_out), "\n")

  if (!is.null(gates$geom_dup)) {
    dup_out <- file.path(opt$outdir, "qc_population_geom_dup.csv")
    write.csv(gates$geom_dup, dup_out, row.names = FALSE)
    cat(normalizePath(dup_out), "\n")
  }

  if (!is.null(gates$adj_band)) {
    adj_out <- file.path(opt$outdir, "qc_adj_factor_bands.csv")
    write.csv(gates$adj_band, adj_out, row.names = FALSE)
    cat(normalizePath(adj_out), "\n")
  }

  if (!is.null(gates$reconcile)) {
    reconcile_out <- file.path(opt$outdir, "qc_crosswalk_reconciliation.csv")
    write.csv(gates$reconcile, reconcile_out, row.names = FALSE)
    cat(normalizePath(reconcile_out), "\n")
  }
}
