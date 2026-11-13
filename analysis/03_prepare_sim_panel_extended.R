# 03_prepare_sim_panel_extended.R — historical simulator panel extended to the ENSO seam
#
# Purpose: the historical smoothing (phase 08) can run past the climate DB
#   (which ends with the historical panel, 2024-12) up to the last month with
#   OBSERVED ENSO, so that later observations (the WHO dashboard national
#   counts, simulator.filter.who_obs) condition the smoothing up to the forecast
#   handoff. This step binds:
#     - the historical panel (03_prepare_sim_panel.R) through its last month, and
#     - the projection panel (03_prepare_sim_panel_projection.R) after it, up to
#       --end (default: the latest month of the ONI cache).
#   The projection rows carry zero-filled climate (no climate data past the DB)
#   and ENSO lags from observed ONI up to --end; every lag of every appended
#   month must be observed, so --end may not exceed the latest ONI month. The
#   deployed fits are ENSO-only, so the zero-filled climate enters no model term.
#   With --enso_max_lead, lead months past the observed ONI come from the
#   projection panel's forecast fill: CFSv2, then its --enso_continuation past
#   the last CFSv2 target (only reached when --end + max lead lies beyond it).
# Inputs:  makeSimPanelFile(opt), makeSimPanelProjectionFile(opt) (rebuild the
#          projection panel with current code first), the ONI cache.
# Outputs: makeSimPanelExtendedFile(opt, end) + its .policy.rds stamp (the
#          historical panel's episode record).
# Usage:   Rscript analysis/03_prepare_sim_panel_extended.R --episode_key gid_season
#          Rscript analysis/03_prepare_sim_panel_extended.R --episode_key gid_season --end 2026-08-01

suppressMessages({
  library(dplyr)
  library(optparse)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/enso.R")

option_list <- list(
  make_option(c("-b", "--baseline"), default = "1990_2019", type = "character",
              help = "Baseline period of both input panels [default: %default]"),
  make_option("--episode_key", default = "lp_season", type = "character",
              help = "lp_season | gid_season (selects the _ep-gs panels) [default: %default]"),
  make_option("--enso_max_lead", default = 0L, type = "integer",
              help = paste("Also write enso_lead1..K (deployed ENSO-lead fits) after stitching, from the",
                           "stitched panel's enso_lag0 (observed ONI) and, for months after --end, the",
                           "projection panel's filled series (CFSv2 point estimate, then continuation)",
                           "[default: %default]")),
  make_option("--end", default = "auto", type = "character",
              help = "Last month to include (YYYY-MM-01), or 'auto' = latest ONI month [default: %default]"),
  # Overwrite guard: a panel of another prep stamp is not replaced
  # without --force.
  prepForceOption("extended sim panel")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$episode_key %in% EPISODE_KEYS)
  stop("--episode_key must be one of ", paste(EPISODE_KEYS, collapse = " | "))
opt$episodes <- list(key = opt$episode_key)

oni_last <- latest_oni_month()
end <- if (identical(opt$end, "auto")) oni_last else as.Date(opt$end)
if (is.na(end) || format(end, "%d") != "01") stop("--end must be a first-of-month date, got: ", opt$end)
if (end > oni_last)
  stop(sprintf("--end %s is after the latest observed ONI month (%s): the appended months' ENSO would not be observed.",
               format(end), format(oni_last)))

hist_path <- as.character(makeSimPanelFile(opt))
proj_path <- as.character(makeSimPanelProjectionFile(opt))
for (p in c(hist_path, proj_path)) if (!file.exists(p)) stop("missing input panel: ", p)
hist <- readRDS(hist_path)
proj <- readRDS(proj_path)
hist_last <- max(hist$t_start)
message(sprintf("[sim_panel_extended] historical %s .. %s; projection %s .. %s; end %s",
                format(min(hist$t_start)), format(hist_last), format(min(proj$t_start)),
                format(max(proj$t_start)), format(end)))
if (end <= hist_last) stop("--end ", format(end), " is not after the historical panel's last month ", format(hist_last))

# Schema and scope must agree exactly: the smoothing reads both halves as one panel.
if (!setequal(names(hist), names(proj)))
  stop("column sets differ: only historical = ", paste(setdiff(names(hist), names(proj)), collapse = ", "),
       "; only projection = ", paste(setdiff(names(proj), names(hist)), collapse = ", "))
if (!setequal(unique(hist$gid), unique(proj$gid)))
  stop("gid sets differ between the historical and projection panels")
cls_h <- vapply(hist, function(x) class(x)[1], ""); cls_p <- vapply(proj[names(hist)], function(x) class(x)[1], "")
if (any(cls_h != cls_p))
  stop("column classes differ: ", paste(names(hist)[cls_h != cls_p], collapse = ", "))
for (st in c("episode_key", "episode_digest")) {
  a <- readPrepPolicyStamp(hist_path)[[st]]; b <- readPrepPolicyStamp(proj_path)[[st]]
  if (!identical(a, b)) stop("policy stamp '", st, "' differs: ", a, " vs ", b)
}
ms <- unique(hist[, c("gid", "min_seas")]); mp <- unique(proj[, c("gid", "min_seas")])
if (!isTRUE(all.equal(ms[order(ms$gid), "min_seas"], mp[order(mp$gid), "min_seas"], check.attributes = FALSE)))
  stop("per-gid min_seas differs between the panels")

ext <- proj[proj$t_start > hist_last & proj$t_start <= end, names(hist)]
need <- seq(as.Date(format(hist_last, "%Y-%m-01")), end, by = "month")[-1]
miss <- setdiff(format(need), format(unique(ext$t_start)))
if (length(miss)) stop("projection panel lacks month(s): ", paste(miss, collapse = ", "))
# enso_lag* and, when both input panels were prepped with --enso_max_lead,
# enso_lead* (the column-set check above makes them agree).
enso_cols <- grep("^enso_(lag|lead)[0-9]+$", names(ext), value = TRUE)
if (anyNA(ext[, enso_cols])) stop("NA ENSO lags / leads in the appended months")

# Appended ENSO must be the OBSERVED ONI at every lag (a stale projection panel
# built before the raw-Nino3.4 seam fix carries rolling means instead).
oni <- fetchEnsoIndex(source = "oni", cache_dir = "generated_data/", refresh = FALSE)
lag0 <- unique(ext[, c("t_start", "enso_lag0")])
dev <- abs(lag0$enso_lag0 - oni$enso_index[match(lag0$t_start, as.Date(oni$date))])
if (anyNA(dev) || max(dev) > 1e-8)
  stop(sprintf("appended enso_lag0 differs from the ONI cache (max |diff| %.3f): rebuild the projection panel (03_prepare_sim_panel_projection.R) with the current code and ENSO caches.",
               suppressWarnings(max(dev, na.rm = TRUE))))

panel <- dplyr::bind_rows(hist, ext) %>% dplyr::arrange(gid, t_start)
if (isTRUE(opt$enso_max_lead > 0L)) {
  panel <- addEnsoLeadsFromLag0(panel, opt$enso_max_lead, extra = list(proj))
  lead_cols <- grep("^enso_lead[0-9]+$", names(panel), value = TRUE)
  if (anyNA(panel[, lead_cols])) {
    bad <- panel$t_start[!stats::complete.cases(panel[, lead_cols])]
    stop("NA ENSO lead(s) in the extended panel (months ", format(min(bad)), " .. ", format(max(bad)),
         "): the projection panel does not reach --end + --enso_max_lead.")
  }
  message("Wrote ", length(lead_cols), " enso_lead columns (observed ONI to ", format(end),
          ", projection-panel ENSO after)")
}
out <- as.character(makeSimPanelExtendedFile(opt, end))
# Stamp built before the panel is written (same fields and values).
.up <- prepStampInherit(hist_path)
.stamp <- prepPolicyStamp(
  pop_source = .up$pop_source, composite_policy = .up$composite_policy,
  cover_threshold = .up$cover_threshold, composite_tiling = .up$composite_tiling,
  episode_key = .up$episode_key, episode_digest = .up$episode_digest,
  counts = list(track = "sim_panel_extended", baseline = opt$baseline,
                source_hist_panel = hist_path, source_projection_panel = proj_path,
                hist_last = format(hist_last), end = format(end),
                enso_max_lead = as.integer(opt$enso_max_lead %||% 0L)))
# Overwrite guard (the rule of 03_prepare_occ.R): a panel of
# another stamp (another vintage's historical panel, another lead count) or
# with none is not replaced without --force.
assertPrepStampOverwrite(out, .stamp, force = isTRUE(opt$force),
                         rerun_hint = "Rscript analysis/03_prepare_sim_panel_extended.R ... --force")
saveRDS(panel, out)
persistPrepPolicyStamp(out, .stamp)
message(sprintf("[sim_panel_extended] wrote %s: %d rows, %s .. %s (%d appended months, zero-filled climate, observed ONI)",
                out, nrow(panel), format(min(panel$t_start)), format(max(panel$t_start)), length(need)))
