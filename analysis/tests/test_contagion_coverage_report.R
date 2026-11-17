# test_contagion_coverage_report.R — HARD: the contagion coverage report
#   (stamp only).
#
# Usage:
#   Rscript analysis/tests/test_contagion_coverage_report.R
#
# Details:
#   Checks the share-column names [1], contagion_coverage_shares() [2], the
#   coverage grids [3], the global activity grid [4] and table D33 [5]. Quick:
#   synthetic grids only.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lubridate)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c(
  "base",
  "gid",
  "paths",
  "occ_data",
  "coverage",
  "spatial_contagion"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
d <- as.Date

# Share-column names nbr_cov_lag1_r<rho>[_pg|_part<tok>] never match
# nbr_covariate_pattern() (they are report columns, not model terms).
cat("\n[1] names\n")
nm <- c(
  nbr_coverage_share_name(500L),
  nbr_coverage_share_name(500L, "group"),
  nbr_coverage_share_name(500L, "registry", "hmm_clusters")
)
expect(
  identical(
    nm,
    c(
      "nbr_cov_lag1_r500",
      "nbr_cov_lag1_r500_pg",
      "nbr_cov_lag1_r500_parthmmclusters"
    )
  ),
  "nbr_cov_lag1_r500[_pg|_parthmmclusters]"
)
expect(
  !any(grepl(nbr_covariate_pattern(), c(nm, "own_cov_lag1"))),
  "no share column matches nbr_covariate_pattern()"
)
expect(
  identical(
    sub("^nbr_cov_lag1_", "nbr_lag1_", nm[3]),
    nbr_covariate_name(500L, "raw", "registry", "hmm_clusters")
  ),
  "same tail as the raw covariate"
)
expect(
  identical(
    as.character(contagion_coverage_csv_path("g/x_ep-gs_nz60.rds")),
    "g/x_ep-gs_nz60.contagion_coverage.csv"
  ),
  "<panel>.contagion_coverage.csv"
)

# contagion_coverage_shares(): the raw-kernel coverage share over the same
# masked universe as the covariate (full coverage -> 1, isolated -> 0,
# partition masks respected).
cat("\n[2] coverage shares\n")
u <- c("AAA_1", "AAA_2", "BBB_1", "CCC_1")
D <- matrix(
  c(
    Inf,
    100,
    300,
    1e6,
    100,
    Inf,
    250,
    1e6,
    300,
    250,
    Inf,
    1e6,
    1e6,
    1e6,
    1e6,
    Inf
  ),
  4,
  dimnames = list(u, u)
)
months <- d(c("2016-02-01", "2016-03-01"))
full <- tidyr::expand_grid(ADM1 = u, t_start = months) %>%
  mutate(covered = TRUE)
sh <- contagion_coverage_shares(full, D, rhos = 500L)
expect(
  all(abs(sh$nbr_cov_lag1_r500[sh$ADM1 != "CCC_1"] - 1) < 1e-12),
  "full coverage -> share 1 (CCC is effectively isolated: exp(-2000) underflows to 0)"
)
half <- full %>% mutate(covered = ADM1 != "BBB_1")
sh2 <- contagion_coverage_shares(half, D, rhos = 500L)
K <- exp(-D / 500)
w <- K["AAA_1", ] / sum(K["AAA_1", ])
expect(
  abs(sh2$nbr_cov_lag1_r500[sh2$ADM1 == "AAA_1"][1] - (1 - w[["BBB_1"]])) <
    1e-12,
  "share = kernel weight of the covered neighbours (AAA_1: 1 - w_BBB)"
)
part <- data.frame(iso3 = c("AAA", "BBB"), group = c("g1", "g1"))
sh3 <- contagion_coverage_shares(
  half,
  D,
  rhos = 500L,
  partitions = list(hmm_clusters = part)
)
expect(
  "nbr_cov_lag1_r500_parthmmclusters" %in%
    names(sh3) &&
    all(sh3$nbr_cov_lag1_r500_parthmmclusters[sh3$ADM1 == "CCC_1"] == 0),
  "partition variant: a unit outside the registry gets 0 (as the covariate's NA -> 0)"
)

# adm1_coverage_grid(): an ADM1 is covered at t when a gid of it is covered at
# t-1; contagion_coverage_gid(): coverage_mask (member bridge, ADM1 windows)
# UNION recorded episode months, national zeros only when on.
cat("\n[3] coverage grids\n")
cg <- tibble(
  gid = c("AAA.1.1_1", "AAA.2_1"),
  month = d(c("2016-01-01", "2016-02-01"))
)
gr <- adm1_coverage_grid(
  cg,
  tidyr::expand_grid(ADM1 = c("AAA_1", "AAA_2"), t_start = months),
  lag = 1L
)
expect(
  identical(gr$covered, c(TRUE, FALSE, FALSE, TRUE)),
  "covered at t when a gid of the ADM1 is covered at t-1 (AAA_1 Feb <- Jan; AAA_2 Mar <- Feb)"
)
td <- tempfile("cov")
dir.create(td)
agg <- file.path(td, "agg.rds")
nzp <- file.path(td, "nz.rds")
saveRDS(
  tibble(
    location_period_id = c("lp1", "lp2"),
    admin_level = c("ADM2", "ADM1"),
    TL = d(c("2016-01-01", "2016-04-01")),
    TR = d(c("2016-01-31", "2016-04-30"))
  ),
  agg
)
saveRDS(
  tibble(country = "AAA", month = d("2016-06-01"), source = "national_and_who"),
  nzp
)
br <- tibble(
  location_period_id = c("lp1", "lp2"),
  gid = c("AAA.1.1_1", "AAA.2_1"),
  country = "AAA"
)
ep <- tibble(gid = "AAA.1.1_1", TL = d("2016-02-10"), TR = d("2016-03-05"))
cv <- suppressMessages(contagion_coverage_gid(
  br,
  agg,
  NULL,
  ep,
  nz_on = FALSE,
  nz_path = nzp
))
expect(
  setequal(
    paste(cv$gid, cv$month),
    paste(
      c("AAA.1.1_1", "AAA.1.1_1", "AAA.1.1_1", "AAA.2_1"),
      d(c("2016-01-01", "2016-02-01", "2016-03-01", "2016-04-01"))
    )
  ),
  "reporting windows (ADM1 window honoured) UNION episode months; no national zeros when off"
)
cv2 <- suppressMessages(contagion_coverage_gid(
  br,
  agg,
  NULL,
  ep,
  nz_on = TRUE,
  nz_path = nzp
))
expect(
  sum(cv2$month == d("2016-06-01")) == 2L,
  "national-zero month covers both gids when on (explicit path)"
)

# adm1_activity_grid(): the global-month activity grid.
cat("\n[4] global activity grid (M1)\n")
spans <- tibble(ADM1 = "BBB_1", TL = d("2016-02-05"), TR = d("2016-02-20"))
ag <- adm1_activity_grid(
  spans,
  c("AAA_1", "BBB_1"),
  d(c("2016-03-01", "2016-04-01"))
)
expect(
  identical(ag$ongoing_outbreak_lag1[ag$ADM1 == "BBB_1"], c(TRUE, FALSE)) &&
    !any(ag$ongoing_outbreak_lag1[ag$ADM1 == "AAA_1"]),
  "BBB_1 active in Feb reads as lag-1 activity in March on the global grid"
)

# Table D33 columns pinned; the step is off without the sidecar.
cat("\n[5] table D33\n")
for (m in c("config", "provenance", "results_core", "national_zero")) {
  source(sprintf("analysis/utils/%s.R", m))
}
source("analysis/utils/results_derive.R")
pan <- tibble(
  country = "AAA",
  t_start = months,
  ADM1 = "AAA_1",
  own_cov_lag1 = c(1, 0),
  nbr_cov_lag1_r500 = c(0.8, 0)
)
side <- tibble(
  ADM1 = "AAA_1",
  t_start = months,
  month = months %m-% months(1),
  covered = c(TRUE, FALSE),
  in_universe = TRUE
)
s39 <- contagion_coverage_table(
  pan,
  side,
  "nbr_cov_lag1_r500",
  data.frame(iso3 = "AAA", group = "enso_east")
)
expect(identical(names(s39), CONTAGION_COVERAGE_COLS), "D33 columns pinned")
expect(
  s39$panel_rows == 2L &&
    abs(s39$nbr_cov_share_mean - 0.4) < 1e-12 &&
    s39$rows_zero_cover == 1L &&
    s39$own_uncovered_rows == 1L &&
    s39$adm1_months == 2L &&
    s39$adm1_months_covered == 1L &&
    s39$group == "enso_east",
  "D33 counts"
)
ids <- function(p) {
  names(results_steps(list2env(
    list(
      FIG1_TRACKS = "occurrence",
      HC_SCOPE = list(enable = FALSE),
      LAG_AVG = list(enable = FALSE),
      CONTAGION_COV_PATH = p
    ),
    parent = emptyenv()
  )))
}
expect(
  "report_contagion_coverage" %in%
    ids("x.csv") &&
    !"report_contagion_coverage" %in% ids(""),
  "D33 step on only with the sidecar"
)

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
