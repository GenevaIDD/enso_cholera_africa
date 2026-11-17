# test_supp_hindcast_clim_skill.R — HARD: the supplementary hindcast
#   skill-vs-climatology figure, its frame helper and builder.
#
# Usage:
#   Rscript analysis/tests/test_supp_hindcast_clim_skill.R
#
# Details:
#   Synthetic phase-09 skill tables: no artefact, no data/. The figure is
#   10 --figures H --supplementary true (hindcast/figH_supp_skill_vs_
#   climatology); helper and builder in utils/figures_supp.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/figures_supp.R")
TH <- ggplot2::theme_bw(base_size = 7)
LEADS <- c("0-3", "4-8", "9+", "all")

# [1] hindcast_clim_skill_frame(): occurrence rows are the window grain, the
#     fair rule and the given reference only; case rows the country grain, the
#     given support and reference, crps_raw -> cases_raw and crps_log1p ->
#     cases_log; country / region_year scopes, other leads and twcrps out;
#     lead_class a factor in lead_levels order; a missing reference stops.
cat("# [1] hindcast_clim_skill_frame()\n")
occ <- expand_grid(
  scope_level = c("africa", "region", "region_year"),
  grain = c("window", "calendar_year"),
  rule = c("crps_emp", "crps_fair"),
  reference = c("loo_gid", "insample_gid"),
  lead_class = LEADS
) |>
  mutate(
    scope = ifelse(scope_level == "africa", "africa", "enso_west"),
    crpss = 0.2,
    crpss_lo = 0.1,
    crpss_hi = 0.3,
    n_scored = 10L
  )
occ$crpss[
  occ$scope_level == "region" &
    occ$grain == "window" &
    occ$rule == "crps_fair" &
    occ$reference == "loo_gid" &
    occ$lead_class == "0-3"
] <- 0.42
cas <- expand_grid(
  scope_level = c("africa", "region", "country"),
  grain = c("country", "adm1"),
  support = c("atrisk_gids", "exact_gid"),
  reference = c("clim_pooled", "clim_launch_state"),
  rule = c("crps_raw", "crps_log1p", "twcrps"),
  lead_class = LEADS
) |>
  mutate(
    scope = case_when(
      scope_level == "africa" ~ "africa",
      scope_level == "region" ~ "enso_west",
      TRUE ~ "NGA"
    ),
    crpss = ifelse(rule == "crps_raw", -0.5, 0.3),
    crpss_lo = ifelse(rule == "crps_raw", -2.5, 0.1),
    crpss_hi = 0.6,
    n_scored = 5L
  )
d <- hindcast_clim_skill_frame(
  occ,
  cas,
  "loo_gid",
  "clim_pooled",
  "atrisk_gids",
  lead_levels = LEADS
)
expect(
  identical(sort(unique(d$target)), c("cases_log", "cases_raw", "onsets")) &&
    nrow(d) == 3L * 2L * length(LEADS) &&
    setequal(d$scope, c("africa", "enso_west")),
  "frame: three targets x (Africa + region) x leads; country / region_year scopes out"
)
expect(
  identical(
    d$crpss[
      d$target == "onsets" & d$scope == "enso_west" & d$lead_class == "0-3"
    ],
    0.42
  ),
  "frame: occurrence row is window grain, crps_fair, the given reference"
)
expect(
  all(d$crpss[d$target == "cases_raw"] == -0.5) &&
    all(d$crpss[d$target == "cases_log"] == 0.3),
  "frame: crps_raw -> cases_raw, crps_log1p -> cases_log (twcrps out)"
)
expect(
  is.factor(d$lead_class) && identical(levels(d$lead_class), LEADS),
  "frame: lead_class a factor in lead_levels order"
)
d2 <- hindcast_clim_skill_frame(occ, cas, lead_levels = c("0-3", "all"))
expect(
  identical(sort(unique(as.character(d2$lead_class))), c("0-3", "all")),
  "frame: other leads dropped"
)
expect(
  inherits(
    try(hindcast_clim_skill_frame(occ, cas, "regional"), silent = TRUE),
    "try-error"
  ),
  "frame: a reference absent from the skill table stops"
)

# [2] build_figS_hindcast_clim_skill(): builds; scopes in display order with
#     their labels; an interval below ylim_lo is clipped there and marked,
#     aligned with its interval (same dodge).
cat("# [2] build_figS_hindcast_clim_skill()\n")
lab <- c(africa = "Africa", enso_west = "GW")
col <- c(africa = "grey25", enso_west = "#66a61e")
p <- build_figS_hindcast_clim_skill(d, lab, col, TH, ylim_lo = -1)
b <- try(suppressMessages(ggplot2::ggplot_build(p)), silent = TRUE)
expect(!inherits(b, "try-error"), "builder: builds")
expect(
  identical(levels(p$data$scope), c("Africa", "GW")) &&
    length(levels(p$data$target)) == 3L,
  "builder: scope labels in display order, three panels"
)
rng <- b$data[[2]]
tri <- b$data[[4]]
tri <- tri[!is.na(tri$y), ] # unclipped rows carry NA
expect(
  all(rng$ymin >= -1) &&
    sum(p$data$clipped) == 2L * length(LEADS) &&
    nrow(tri) == 2L * length(LEADS),
  "builder: intervals below ylim_lo clipped there and marked once each"
)
expect(
  all(tri$y == -1) && all(round(tri$x, 6) %in% round(rng$x[rng$ymin == -1], 6)),
  "builder: the clip markers sit on their intervals (same dodge)"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
