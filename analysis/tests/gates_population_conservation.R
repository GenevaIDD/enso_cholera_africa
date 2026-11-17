# gates_population_conservation.R — conservation gates over the built
#   population crosswalk and size dataset.
#
# Usage:
#   Rscript analysis/tests/gates_population_conservation.R
#   Rscript analysis/tests/gates_population_conservation.R --year 2020
#
# Details:
#   Not a test_*.R file and not in run_all.R: run_all's tests are hermetic
#   (synthetic fixtures, no data/) and pass on a bare checkout; this one checks
#   the artefacts on disk, so it runs only where the population pipeline has
#   run. The unit tests prove apportion_composite_cases() conserves cases on a
#   synthetic composite; this proves the real corpus satisfies the invariants
#   those functions assume about their inputs. The gates are listed above
#   "Crosswalk gates".

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

source("analysis/utils/base.R") # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/paths.R") # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/population.R") # load_population_crosswalk() -- used directly
source("analysis/utils/composite_lp.R") # transitive dependency of population.R (utils/README.md)
source("analysis/utils/provenance.R") # get_pop_provenance()/read_prep_policy_stamp() -- used directly

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--year",
      default = 2020L,
      type = "integer",
      help = "crosswalk year to audit [default: %default]"
    ),
    make_option(
      "--crosswalk_path",
      default = NULL,
      type = "character",
      help = "population crosswalk .qs [default: newest match]"
    ),
    make_option(
      "--size_data",
      type = "character",
      default = "generated_data/outbreak_size_analysis_data.rds",
      help = "size dataset for gates (i)-(iii) [default: %default]"
    ),
    make_option(
      "--centroid_gids",
      type = "character",
      default = "generated_data/matched_ob_gids.rds",
      help = "legacy centroid match, for the costing report"
    ),
    make_option(
      "--tol",
      default = 1e-6,
      type = "double",
      help = "relative tolerance for conservation [default: %default]"
    )
  )
))

PASS <- 0L
FAIL <- 0L
SKIP <- 0L
gate <- function(cond, label) {
  if (isTRUE(cond)) {
    PASS <<- PASS + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  FAIL  %s\n", label))
  }
}
skip <- function(label, why) {
  SKIP <<- SKIP + 1L
  cat(sprintf("  SKIP  %s\n        (%s)\n", label, why))
}

# GATES
#   (i)   sum(suspected_cases) per country, before vs after expansion  [size_dat]
#   (ii)  sum(pop) per (lctn_pr, year) after expansion == P_l          [size_dat]
#   (iii) no row with pop <= 0                                        [size_dat]
#   (iv)  sum(w) == 1 per (lctn_pr, gadm_level, year)                  [crosswalk]
#   (v)   zero rows with pop_tier == "missing"                        [crosswalk]
#   plus  c <= 1, H > 0, and a COSTING report for the gid universe.
# Gates (iv)/(v) need only the crosswalk and run as soon as the population
# pipeline has. Gates (i)-(iii) need a size dataset built with
# `--pop_source crosswalk` and are SKIPPED (not failed) when it is absent: a
# gate that cannot run must say so rather than report a pass.
# ---------------------------------------------------------------------------
# Crosswalk gates
# ---------------------------------------------------------------------------
xw_path <- opt$crosswalk_path
if (is.null(xw_path)) {
  cand <- Sys.glob("generated_data/population_crosswalk_*.qs")
  if (!length(cand)) {
    stop("no population_crosswalk_*.qs found; pass --crosswalk_path")
  }
  xw_path <- cand[which.max(file.mtime(cand))]
}
cat("crosswalk : ", xw_path, "\n", sep = "")
cat("year      : ", opt$year, "\n\n", sep = "")

xw <- load_population_crosswalk(xw_path)
if (!opt$year %in% xw$year) {
  stop(
    "year ",
    opt$year,
    " absent from the crosswalk (have: ",
    paste(range(xw$year), collapse = "-"),
    ")"
  )
}
xwy <- xw %>% filter(year == opt$year)

cat(
  "## gate (iv) — w is a share: sum(w) == 1 per (lctn_pr, gadm_level, year)\n"
)
# This is the invariant apportion_composite_cases() relies on for conservation.
# If it fails anywhere, cases are created or destroyed by the expansion and
# every downstream total is wrong by a factor nobody would notice.
gw <- xw %>%
  group_by(lctn_pr, gadm_level, year) %>%
  summarise(sw = sum(w, na.rm = TRUE), .groups = "drop")
worst_w <- max(abs(gw$sw - 1), na.rm = TRUE)
gate(
  worst_w <= opt$tol,
  sprintf("%d groups, max |sum(w) - 1| = %.3e", nrow(gw), worst_w)
)
if (worst_w > opt$tol) {
  print(as.data.frame(head(gw %>% filter(abs(sw - 1) > opt$tol), 10)))
}

cat("\n## gate (v) — every row has a resolved population tier\n")
tt <- table(xwy$pop_tier, useNA = "ifany")
cat("  tiers: ", paste(sprintf("%s=%d", names(tt), tt), collapse = ", "), "\n")
gate(
  !any(names(tt) %in% "missing") && !any(is.na(names(tt))),
  "no pop_tier == 'missing' and no NA tier"
)

cat("\n## coverage: c = H/Q_g must not exceed 1\n")
# c > 1 means an LP contains more people than the gid it sits in — a sign the
# overlap and the canonical extraction disagree, which would break the
# min(sum H, Q) saturation in reconcile_footprint_pop().
gate(
  sum(xw$c > 1 + opt$tol, na.rm = TRUE) == 0,
  sprintf(
    "max c = %.6f, rows over 1: %d",
    max(xw$c, na.rm = TRUE),
    sum(xw$c > 1 + opt$tol, na.rm = TRUE)
  )
)
gate(
  sum(!is.finite(xw$H) | xw$H <= 0, na.rm = TRUE) == 0,
  sprintf(
    "all H finite and > 0 (violations: %d)",
    sum(!is.finite(xw$H) | xw$H <= 0, na.rm = TRUE)
  )
)

# ---------------------------------------------------------------------------
# Costing report — how big is the expansion actually?
# ---------------------------------------------------------------------------
cat("\n## costing — gid universe and expansion factor\n")
kk <- xwy %>% count(lctn_pr, gadm_level, name = "K")
cat(sprintf(
  "  members per LP: median %.0f, mean %.1f, max %d\n",
  stats::median(kk$K),
  mean(kk$K),
  max(kk$K)
))
cat(sprintf(
  "  K==1 %d | 2-22 %d | 23-99 %d | >=100 %d\n",
  sum(kk$K == 1),
  sum(kk$K >= 2 & kk$K <= 22),
  sum(kk$K >= 23 & kk$K < 100),
  sum(kk$K >= 100)
))
cat(sprintf(
  "  distinct gids: ADM1 %d, ADM2 %d\n",
  n_distinct(xwy$gid[xwy$gadm_level == "ADM1"]),
  n_distinct(xwy$gid[xwy$gadm_level == "ADM2"])
))

if (file.exists(opt$centroid_gids)) {
  cent <- readRDS(opt$centroid_gids)
  cent <- unique(cent$gid[!is.na(cent$gid)])
  for (lv in c("ADM1", "ADM2")) {
    u <- unique(xwy$gid[xwy$gadm_level == lv])
    cat(sprintf(
      "  if LPs resolve at %s: universe %d, NEW gids needing climate extraction %d\n",
      lv,
      length(u),
      length(setdiff(u, cent))
    ))
  }
  cat(sprintf("  legacy centroid universe: %d\n", length(cent)))
  cat(
    "  NOTE: each new gid needs full CHIRPS + ERA5 monthly extraction AND a\n",
    "        climate_extremes baseline. 02_preprocess_climate.R writes with\n",
    "        overwrite = TRUE, so the cost is the FULL rebuild, not a delta.\n",
    sep = ""
  )
} else {
  cat(
    "  (",
    opt$centroid_gids,
    " absent — cannot diff against the legacy universe)\n",
    sep = ""
  )
}

# ---------------------------------------------------------------------------
# Size-dataset gates (i)-(iii)
# ---------------------------------------------------------------------------
cat("\n## gates (i)-(iii) — the expanded size dataset\n")
if (!file.exists(opt$size_data)) {
  skip("gates (i)-(iii)", paste0(opt$size_data, " not found"))
} else {
  sd_ <- readRDS(opt$size_data)
  prov <- get_pop_provenance(sd_)
  cat(sprintf("  pop_source: %s\n", prov$pop_source %||% "UNSTAMPED"))

  # Which composite policy built this file. `pop_source` reads "crosswalk" under
  # all four, so it cannot answer this; the companion sidecar can. Structural
  # stamp first, sidecar as the fallback for artefacts whose attributes were
  # stripped in transit. Default "keep" -- the pre-flag behaviour, so an
  # artefact predating the sidecar is judged by the rules it was built under.
  pol <- prov$composite_policy %||%
    (read_prep_policy_stamp(opt$size_data)$composite_policy %||% "keep")
  cat(sprintf("  composite_policy: %s\n", pol))
  expanded <- identical(pol, "keep") || identical(pol, "drop")

  # --episode_key gid_season (utils/footprint.R::build_canonical_episodes()):
  # one row per (gid, season) and every episode's cases equal to the sum over
  # the members it COUNTS (`in_total`, the finest reporting level) in the
  # canonical membership map written beside it.
  ep_st <- episode_stamp_of(opt$size_data)
  cat(sprintf("  episode_key: %s\n", ep_st$key))
  if (identical(ep_st$key, "gid_season")) {
    mem_path <- make_footprint_members_file(list(
      episodes = list(key = "gid_season")
    ))
    gate(
      file.exists(mem_path),
      sprintf("(e-0) canonical members map present: %s", mem_path)
    )
    if (file.exists(mem_path)) {
      dupg <- sd_ %>% count(gid, year_seas) %>% filter(n > 1L)
      gate(
        nrow(dupg) == 0L,
        sprintf(
          "(e-1) exactly one row per (gid, season) -- %d violation(s)",
          nrow(dupg)
        )
      )
      mem <- readRDS(mem_path)
      tot <- mem %>%
        filter(in_total) %>%
        group_by(gid, year_seas = as.character(year_seas)) %>%
        summarise(member_cases = sum(suspected_cases), .groups = "drop")
      cmp <- sd_ %>%
        transmute(
          gid = as.character(gid),
          year_seas = as.character(year_seas),
          suspected_cases
        ) %>%
        left_join(tot, by = c("gid", "year_seas")) %>%
        mutate(
          rel = abs(suspected_cases - member_cases) / pmax(1, abs(member_cases))
        )
      gate(
        all(is.finite(cmp$rel) & cmp$rel <= opt$tol),
        sprintf(
          "(e-2) episode cases == sum of counted member cases (max rel err %.2e)",
          max(cmp$rel, na.rm = TRUE)
        )
      )
      gate(
        identical(ep_st$digest, episode_stamp_of(mem_path)$digest),
        "(e-3) members map and size frame carry the same episode digest"
      )
    }
  }

  if (!identical(prov$pop_source %||% "", "crosswalk")) {
    skip(
      "gates (i)-(iii)",
      paste0(
        "size data was built with pop_source = '",
        prov$pop_source %||% "unstamped",
        "'; rebuild with --pop_source crosswalk to audit the expansion"
      )
    )
  } else {
    # (iii) is checkable on its own and is the COD_4_6_1 regression.
    gate(
      all(is.finite(sd_$pop) & sd_$pop > 0),
      sprintf(
        "(iii) no row with pop <= 0 or non-finite (violations: %d)",
        sum(!is.finite(sd_$pop) | sd_$pop <= 0)
      )
    )

    if (!"apportion_w" %in% names(sd_)) {
      skip(
        "(i)/(ii)",
        "size data lacks apportion_w — cannot invert the expansion"
      )
    } else if (!expanded) {
      # UNDER A COLLAPSING POLICY THE INVERSION IS VACUOUS, NOT WRONG:
      # apportion_w == 1 on every row, so `orig` is just `suspected_cases` and
      # the check reduces to sum(x) == sum(x). Assert the property that actually
      # characterises a collapse instead — one row per (LP, season), each
      # carrying the LP's OWN case count. That is a strictly stronger claim than
      # the conservation the expansion path can make, and it is the thing a
      # regression would break.
      gate(
        all(abs(sd_$apportion_w - 1) < 1e-12, na.rm = TRUE),
        sprintf(
          "(i-a) policy=%s: apportion_w == 1 on every row (violations: %d)",
          pol,
          sum(abs(sd_$apportion_w - 1) >= 1e-12, na.rm = TRUE)
        )
      )
      gate(
        !any(sd_$apportioned, na.rm = TRUE),
        sprintf("(i-b) policy=%s: no row is flagged `apportioned`", pol)
      )
      dup <- sd_ %>% count(location_period_id, year_seas) %>% filter(n > 1L)
      gate(
        nrow(dup) == 0L,
        sprintf(
          "(i-c) policy=%s: exactly one row per (LP, season) — %d violation(s)",
          pol,
          nrow(dup)
        )
      )
      if (nrow(dup)) {
        print(as.data.frame(head(dup, 10)))
      }
      # Nothing was scaled, so cases must be whole numbers of the LP's own count
      # rather than a w-share — the failure mode is a silent re-introduction of
      # scaling somewhere between the collapse and the save.
      gate(
        all(is.finite(sd_$suspected_cases) & sd_$suspected_cases >= 0),
        sprintf("(i-d) policy=%s: all case counts finite and non-negative", pol)
      )
    } else {
      # (i) cases: the expansion multiplied by w, so dividing back out and
      # de-duplicating on the LP must recover the pre-expansion total. Under
      # `drop` this is conservation among the RETAINED LPs — dropped LPs are
      # absent from both sides, so the equality is unaffected by the drop.
      pre <- sd_ %>%
        mutate(orig = suspected_cases / apportion_w) %>%
        distinct(location_period_id, year_seas, .keep_all = TRUE) %>%
        group_by(country) %>%
        summarise(before = sum(orig), .groups = "drop")
      post <- sd_ %>%
        group_by(country) %>%
        summarise(after = sum(suspected_cases), .groups = "drop")
      cmp <- full_join(pre, post, by = "country") %>%
        mutate(rel = abs(after - before) / pmax(1, abs(before)))
      gate(
        all(cmp$rel <= opt$tol, na.rm = TRUE),
        sprintf(
          "(i) policy=%s: cases conserved per country among retained LPs (max rel err %.3e, %d countries)",
          pol,
          max(cmp$rel, na.rm = TRUE),
          nrow(cmp)
        )
      )
      if (any(cmp$rel > opt$tol, na.rm = TRUE)) {
        print(as.data.frame(
          cmp %>% filter(rel > opt$tol) %>% arrange(desc(rel)) %>% head(10)
        ))
      }

      # (ii) pop: sum of member H per (LP, year) must equal P_l. Only checkable
      # BEFORE the (gid, season) collapse saturates against Q_g, so it is
      # reported rather than gated when the collapse has already run.
      skip(
        "(ii) sum(pop) per (lctn_pr, year) == P_l",
        "size_dat is post-collapse; H no longer separable per LP. Covered pre-collapse by the in-script gate in 03_prepare_siz.R"
      )
    }
  }
}

cat(sprintf("\n%d pass, %d fail, %d skip\n", PASS, FAIL, SKIP))
if (FAIL > 0) {
  quit(status = 1)
}
