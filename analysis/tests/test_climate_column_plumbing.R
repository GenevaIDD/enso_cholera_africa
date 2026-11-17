# test_climate_column_plumbing.R — HARD: the column plumbing the local-climate
#   crossbases rely on (utils/occ_data.R, dur_data.R, clim_db.R, clim_compile.R)
#
# Usage:
#   Rscript analysis/tests/test_climate_column_plumbing.R
#
# Details:
#   Checks scale_precip() [1], .clim_dev_wide() [2] and get_stats_gid() [3].
#   Quick-safe (synthetic inputs, a temp SQLite; no fit).
suppressPackageStartupMessages({
  library(dplyr)
  library(DBI)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "gid", "clim_compile", "clim_db", "occ_data")) {
  source(sprintf("analysis/utils/%s.R", m))
}

# Standardised anomalies (dev_sd_*) and SPEI are unitless and never rescaled;
# precipitation levels / anomalies / their lags are still divided by 100.
cat("# [1] scale_precip\n")
df <- tibble(
  uid = 1L,
  precip = 200,
  mean_precip = 100,
  dev_precip = 50,
  mean_precip_lag1 = 300,
  max_cum7 = 40,
  dev_sd_precip = 1.5,
  dev_sd_precip_lag1 = -0.5,
  spei1 = 0.7,
  hlq_precip = 1
)
s <- scale_precip(df)
expect(
  s$precip == 2 &&
    s$mean_precip == 1 &&
    s$dev_precip == 0.5 &&
    s$mean_precip_lag1 == 3 &&
    s$max_cum7 == 0.4,
  "levels, anomalies and lags still / 100"
)
expect(
  s$dev_sd_precip == 1.5 && s$dev_sd_precip_lag1 == -0.5,
  "dev_sd_* not rescaled"
)
expect(s$spei1 == 0.7 && s$hlq_precip == 1, "spei and hlq untouched")

# Identical to the pivot make_surv_data() ran inline (frozen copy below), and
# the multi-gid form keys on gid.
cat("# [2] .clim_dev_wide\n")
months <- seq(as.Date("2000-01-01"), by = "month", length.out = 3)
clim_stats <- tibble(
  gid = "AAA_1_1",
  month = months,
  precip = c(10, 20, 30),
  max_precip = 1:3,
  mean_cum7 = 4:6,
  max_cum7 = 7:9,
  temp_daymean = c(25, 26, 27),
  temp_daymax = 30:32,
  temp_daymin = 18:20
)
vars <- c(get_precip_vars(), get_temp_vars())
clim_extrms <- tidyr::expand_grid(tibble(month = months), name = vars) %>%
  mutate(gid = "AAA_1_1", mean = seq_len(n()), sd = 2, hilo_sd = 0, hilo_q = 0)
frozen <- clim_stats %>% # make_surv_data()'s inline pivot, verbatim
  tidyr::pivot_longer(cols = c(get_precip_vars(), get_temp_vars())) %>%
  inner_join(
    clim_extrms %>% select(month, name, mean, sd),
    by = c("month", "name")
  ) %>%
  mutate(dev = value - mean, dev_sd = (value - mean) / sd) %>%
  select(gid, month, name, mean, dev, dev_sd) %>%
  tidyr::pivot_wider(
    names_from = "name",
    values_from = c("mean", "dev", "dev_sd")
  )
expect(
  identical(
    .clim_dev_wide(clim_stats, clim_extrms, vars, by = c("month", "name")),
    frozen
  ),
  "single-gid form identical to the old inline pivot"
)
cs2 <- bind_rows(
  clim_stats,
  mutate(clim_stats, gid = "BBB_1_1", precip = precip + 1)
)
ce2 <- bind_rows(
  clim_extrms,
  mutate(clim_extrms, gid = "BBB_1_1", mean = mean + 100)
)
w <- .clim_dev_wide(cs2, ce2, vars)
expect(
  nrow(w) == 6 &&
    all(w$mean_precip[w$gid == "BBB_1_1"] > 100) &&
    isTRUE(all.equal(
      w$dev_sd_precip[w$gid == "AAA_1_1"],
      (clim_stats$precip - clim_extrms$mean[clim_extrms$name == "precip"]) / 2
    )),
  "multi-gid form keys on gid; dev_sd = (value - mean) / sd"
)

cat("# [3] get_stats_gid()\n")
db <- tempfile(fileext = ".sqlite")
con <- dbConnect(RSQLite::SQLite(), db, extended_types = TRUE)
dbWriteTable(con, "climate_stats", cs2)
dbDisconnect(con)
r1 <- get_stats_gid(c("AAA.1_1", "BBB_1_1"), path = db)
expect(
  identical(names(r1), c("gid", "month", vars)),
  "default vars: the stats columns"
)
expect(
  nrow(r1) == 6 && identical(sort(unique(r1$gid)), c("AAA.1_1", "BBB_1_1")),
  "caller's gid spelling returned"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
