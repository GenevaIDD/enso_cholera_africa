# test_results_rain_enso.R — HARD: Figure 1d's rainfall tables (09 step
#   rain_enso; helpers in utils/clim_compile.R).
#
# Usage:
#   Rscript analysis/tests/test_results_rain_enso.R
#
# Details:
#   Quick-safe (no fit, a few seconds). Input from 03_prepare_rain_adm1.R.
#   Checks adm1_monthly_precip(), precip_dev_series(), the ONI correlations,
#   the ADM1 population and rain seasons on synthetic panels, and [5]
#   (data-dependent, SKIPs without the inputs) that the helpers reproduce the
#   reference CSVs on the grouping those were built with.

suppressPackageStartupMessages({
  library(dplyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "gid", "paths", "clim_compile", "rain_enso")) {
  source(sprintf("analysis/utils/%s.R", m))
}

# [1] adm1_monthly_precip(): an ADM1's own row wins, an ADM1 without one takes
#     the mean of its ADM2 rows; Ghana-style gids get their separator.
cat("# [1] adm1_monthly_precip\n")
d0 <- as.numeric(as.Date("2000-01-01"))
# flattened GADM gids: ADM1 "AAA_1_1" (-> ADM1 id "AAA_1"), ADM2 "AAA_1_2_1"
pr <- data.frame(
  gid = c("AAA_1_1", "AAA_1_1_1", "AAA_2_1_1", "AAA_2_2_1", "GHA14_2"),
  month = d0,
  precip = c(10, 99, 4, 8, 7)
)
a <- adm1_monthly_precip(pr)
expect(
  identical(sort(a$adm1), c("AAA_1", "AAA_2", "GHA_14")),
  "one row per ADM1 (Ghana separator inserted)"
)
expect(
  a$precip[a$adm1 == "AAA_1"] == 10 && a$src[a$adm1 == "AAA_1"] == "adm1_row",
  "an ADM1's own row wins over its ADM2 rows"
)
expect(
  a$precip[a$adm1 == "AAA_2"] == 6 && a$src[a$adm1 == "AAA_2"] == "adm2_mean",
  "an ADM1 without its own row takes the mean of its ADM2 rows"
)
expect(
  all(a$month == as.Date("2000-01-01")) &&
    all(a$mo == 1L) &&
    all(a$year == 2000L),
  "month stored as days since 1970 becomes a Date with year / mo"
)

# --- synthetic panel ---------------------------------------------------------
set.seed(11)
months <- seq(as.Date("1980-01-01"), as.Date("2021-12-01"), by = "month")
oni_v <- as.numeric(stats::arima.sim(
  list(ar = 0.9),
  length(months) + 2,
  sd = 0.4
))
oni <- data.frame(
  year = as.integer(format(months, "%Y")),
  month = as.integer(format(months, "%m")),
  oni = oni_v[-(1:2)]
)
oni_lag2 <- oni_v[seq_along(months)] # ONI two months earlier
mo <- as.integer(format(months, "%m"))
clim_of <- function(peak) 70 + 50 * cos(2 * pi * (mo - peak) / 12) # driest month 20 mm, never "dry"
units <- list(
  AAA_1 = list(clim = clim_of(7), k = 0.4), # g1, answers ONI(t - 2)
  AAA_2 = list(clim = clim_of(2), k = 0.4), # g1
  BBB_1 = list(clim = clim_of(2), k = 0.4), # g1
  CCC_1 = list(clim = clim_of(7), k = 0.4), # not in the grouping
  DDD_1 = list(clim = clim_of(1), k = 0)
) # g2, no ENSO signal
units$DDD_1$clim[mo == 7] <- 5 # one dry calendar month
adm1_m <- bind_rows(lapply(names(units), function(u) {
  x <- units[[u]]
  tibble(
    adm1 = u,
    month = months,
    precip = pmax(
      0,
      x$clim * (1 + x$k * oni_lag2) + rnorm(length(months), 0, 3)
    ),
    src = "adm1_row",
    iso3 = substr(u, 1, 3),
    year = as.integer(format(months, "%Y")),
    mo = mo
  )
}))
G <- c(AAA = "g1", BBB = "g1", DDD = "g2")

# [2] precip_dev_series(): baseline sd units, dry calendar months dropped, the
#     centred 3-month running mean (>= 2 of 3 months) -- against a direct
#     computation on one unit.
cat("# [2] precip_dev_series\n")
dev <- precip_dev_series(adm1_m, c(1990L, 2019L), dry_mm_month = 10)
expect(
  identical(dev$units, names(units)) && nrow(dev$D) == length(months),
  "months x ADM1 matrix"
)
expect(
  isTRUE(all.equal(dev$dry_share, (1 / 12) / length(units))),
  "the one dry calendar month is dropped"
)
expect(
  all(
    is.na(dev$D[
      mo == 7 & seq_along(mo) > 1 & seq_along(mo) < length(mo),
      "DDD_1"
    ]) ==
      FALSE
  ),
  "a dropped month inside the window still has a 3-month mean (2 of 3 months)"
)
u <- adm1_m[adm1_m$adm1 == "AAA_1", ]
bl <- u$year >= 1990 & u$year <= 2019
mu <- tapply(u$precip[bl], u$mo[bl], mean)
sdv <- tapply(u$precip[bl], u$mo[bl], sd)
z <- (u$precip - mu[u$mo]) / sdv[u$mo]
D_ref <- vapply(
  seq_along(z),
  function(i) mean(z[max(1, i - 1):min(length(z), i + 1)]),
  0
)
expect(
  isTRUE(all.equal(unname(dev$D[, "AAA_1"]), D_ref)),
  "D = centred 3-month mean of the baseline-sd deviation"
)

# [3] precip_oni_corr() / precip_oni_group_corr(), on a synthetic panel whose
#     rainfall answers ONI two months earlier: r peaks at lag 2, p_r equals
#     cor.test(), min_n leaves r NA, countries outside the grouping are no
#     series, agree_r counts the countries significant with the group's sign.
cat("# [3] precip_oni_corr / precip_oni_group_corr\n")
st <- precip_oni_group_corr(dev, G, oni, lags = -6:6, min_n = 30)
expect(
  identical(sort(unique(st$series)), c("g1", "g2")),
  "one series per group; ungrouped countries are none"
)
expect(nrow(st) == 2 * 12 * 13, "series x centre month x lag rows")
g1 <- st[st$series == "g1" & st$centre == 6, ]
expect(
  g1$lag[which.max(g1$r)] == 2L && max(g1$r) > 0.8,
  sprintf("g1's r peaks at lag 2 (r = %.2f)", max(g1$r))
)
cell <- g1[g1$lag == 2, ]
rows <- which(dev$mo == 6)
yrs <- dev$yr[rows]
y <- rowMeans(dev$D[rows, c("AAA_1", "AAA_2", "BBB_1")], na.rm = TRUE)
x <- oni$oni[match(paste(yrs, 4), paste(oni$year, oni$month))]
ok <- is.finite(x) & is.finite(y)
ct <- cor.test(y[ok], x[ok])
expect(
  isTRUE(all.equal(cell$r, unname(ct$estimate))) &&
    isTRUE(all.equal(cell$p_r, ct$p.value)) &&
    cell$n == sum(ok),
  "r / p_r / n of one cell == cor.test() on the group mean"
)
expect(
  isTRUE(all.equal(
    c(cell$r_lo, cell$r_hi),
    unname(ct$conf.int[1:2]),
    tolerance = 1e-3
  )),
  "r_lo / r_hi = the Fisher-z 95% interval"
)
expect(
  cell$n_ctry == 2L && cell$agree_r == 1,
  "agree_r: both g1 countries significant with g1's sign"
)
expect(all(st$n_ctry[st$series == "g2"] == 1L), "g2 has one country")
st_n <- precip_oni_corr(
  dev,
  ifelse(substr(dev$units, 1, 3) %in% names(G), G[substr(dev$units, 1, 3)], NA),
  oni,
  lags = 0,
  min_n = 1000
)
expect(
  all(is.na(st_n$r)) && all(is.na(st_n$p_r)),
  "r and p_r are NA below min_n years"
)

# [4] adm1_population_from_gid() / adm1_rain_seasons(): the ADM1 population
#     row wins over its ADM2 sum; population vs ADM1-count weighting pick
#     different seasons; a season crossing the year end stays one.
cat("# [4] adm1_population_from_gid / adm1_rain_seasons\n")
pop <- data.frame(
  gid = c(
    "AAA_1_1",
    "AAA_1_1_1",
    "AAA_2_1_1",
    "AAA_2_2_1",
    "BBB_1_1",
    "DDD_1_1",
    "AAA_1_1"
  ),
  gadm_level = c("ADM1", "ADM2", "ADM2", "ADM2", "ADM1", "ADM1", "ADM1"),
  year = c(2020L, 2020L, 2020L, 2020L, 2020L, 2020L, 2015L),
  Q = c(1e6, 5, 10, 20, 50, 100, 1)
)
pa <- adm1_population_from_gid(pop, 2020L)
expect(
  pa$Q[pa$adm1 == "AAA_1"] == 1e6,
  "an ADM1's population row wins over its ADM2 sum (and the year is filtered)"
)
expect(
  pa$Q[pa$adm1 == "AAA_2"] == 30,
  "an ADM1 without a row takes the sum of its ADM2 rows"
)
rs <- adm1_rain_seasons(
  adm1_m,
  G,
  pa,
  c(1990L, 2019L),
  share_cut = 0.5,
  min_annual_mm = 100
)
sp <- rs$seasons[rs$seasons$weighting == "pop" & rs$seasons$group == "g1", ]
sa <- rs$seasons[rs$seasons$weighting == "adm1" & rs$seasons$group == "g1", ]
expect(
  nrow(sp) == 1L && sp$first <= 7 && sp$last >= 7,
  sprintf("pop-weighted g1 season around July (%s)", sp$months)
)
expect(
  nrow(sa) == 1L && (sa$first >= 11 || sa$first <= 2),
  sprintf("ADM1-count g1 season around February (%s)", sa$months)
)
sd_ <- rs$seasons[rs$seasons$weighting == "pop" & rs$seasons$group == "g2", ]
expect(
  nrow(sd_) == 1L && sd_$first > sd_$last,
  sprintf("g2's season crosses the year end as one season (%s)", sd_$months)
)
expect(
  rs$n_units == 4L && rs$n_pop_matched == 4L,
  "grouped ADM1s only (CCC left out), all with population"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
