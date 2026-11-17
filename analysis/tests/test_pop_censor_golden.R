# test_pop_censor_golden.R — HARD: golden characterization of
#   resolve_gid_population() (population.R) and censor_sim_by_country().
#
# Usage:
#   Rscript analysis/tests/test_pop_censor_golden.R
#   Rscript analysis/tests/test_pop_censor_golden.R --capture
#
# Details:
#   All four population tiers. The T2 WPP rescale is vectorised (one
#   rescale_population_wpp() call) and the censor derives ISO3 once per
#   unique gid; every output must stay bit-identical. One synthetic fixture
#   is run through both and compared with identical() against
#   tests/golden/pop_censor_golden.rds; --capture writes the golden. No
#   generated_data, no network.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
})
cli <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--capture",
      action = "store_true",
      default = FALSE,
      help = "write the golden file instead of comparing"
    )
  )
))
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c(
  "base",
  "config",
  "paths",
  "provenance",
  "gid",
  "population",
  "sim_io"
)) {
  source(sprintf("analysis/utils/%s.R", m))
}
err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

GOLDEN_PATH <- "analysis/tests/golden/pop_censor_golden.rds"
cat(
  "# population tiers + country censor golden (mode:",
  if (cli$capture) "capture" else "compare",
  ")\n"
)

# gid_pop rows deliberately out of year order: a nearest-year tie (2016 vs 2020
# for 2018) is decided by row order. ZZZ is absent from the WPP table (T2 NA).
# Fixture: nearest-year ties, unsorted gid_pop rows, an unknown gid, NA gid
# and year, a country missing from the WPP table, a T3 parent residual,
# wpp = NULL.
gp <- tibble::tribble(
  ~gid        , ~iso3 , ~year , ~Q   , ~gadm_level ,
  "KEN_1_1"   , "KEN" , 2020L , 1200 , "ADM1"      ,
  "KEN_1_1"   , "KEN" , 2016L , 1000 , "ADM1"      ,
  "KEN_2_1"   , "KEN" , 2016L ,  500 , "ADM1"      ,
  "KEN_2_1"   , "KEN" , 2017L ,  520 , "ADM1"      ,
  "KEN_2_1"   , "KEN" , 2030L ,  900 , "ADM1"      ,
  "ETH_3_1"   , "ETH" , 2015L , 3000 , "ADM1"      ,
  "ETH_3_1"   , "ETH" , 2018L , 3300 , "ADM1"      ,
  "ETH_3_1_1" , "ETH" , 2018L , 1000 , "ADM2"      ,
  "ETH_3_2_1" , "ETH" , 2018L ,  800 , "ADM2"      ,
  "ZZZ_1_1"   , "ZZZ" , 2016L ,   50 , "ADM1"
)
wpp <- tidyr::expand_grid(iso3 = c("KEN", "ETH"), year = 2005:2035) %>%
  dplyr::mutate(
    pop_total = 1e6 * (1 + 0.03 * (year - 2005)) * ifelse(iso3 == "KEN", 1, 2.5)
  )
q_gid <- c(
  "KEN_1_1",
  "KEN_1_1",
  "KEN_1_1",
  "KEN_2_1",
  "KEN_2_1",
  "ETH_3_1",
  "ETH_3_1",
  "ETH_3_3_1",
  "ETH_3_3_1",
  "ZZZ_1_1",
  "NOPE_1_1",
  NA,
  "KEN_2_1",
  "KEN_1_1"
)
q_year <- c(
  2016L,
  2018L,
  2010L,
  2023L,
  2035L,
  2016L,
  2012L,
  2018L,
  2019L,
  2014L,
  2018L,
  2018L,
  NA,
  2018L
)

pop <- list(
  na_wpp = resolve_gid_population(q_gid, q_year, gp, wpp, on_missing = "na"),
  na_nowpp = resolve_gid_population(q_gid, q_year, gp, NULL, on_missing = "na"),
  recycled = resolve_gid_population(
    q_gid[1:7],
    2019L,
    gp,
    wpp,
    on_missing = "na"
  ),
  factor = resolve_gid_population(factor(q_gid[1:7]), q_year[1:7], gp, wpp),
  err = err_of(resolve_gid_population(q_gid, q_year, gp, wpp)),
  empty = resolve_gid_population(character(), integer(), gp, wpp),
  rescale = rescale_population_wpp(
    c(10, 20, 30),
    c("ken", "ETH", "ZZZ"),
    2016L,
    c(2018L, 2010L, 2018L),
    wpp
  )
)

hz <- c(KEN = "2027-01-01", ETH = "2027-06-01")
months <- seq(as.Date("2026-11-01"), as.Date("2027-08-01"), by = "month")
ms <- tidyr::expand_grid(
  replicate_id = 1:3,
  gid = c("KEN_1_1", "ETH_2_1", "KEN_1_2_1", "KEN.3.1_1"),
  month = months
) %>%
  dplyr::mutate(Y = rep_len(0:1, dplyr::n()))
ev <- tibble::tibble(
  gid = c("KEN_1_1", "ETH_2_1", "ETH_2_1"),
  start_month = as.Date(c("2026-12-01", "2027-05-01", "2027-07-01")),
  end_month = as.Date(c("2027-03-01", "2027-05-01", "2027-09-01"))
)
cen <- list(
  drop = censor_sim_by_country(ms, hz),
  flag = censor_sim_by_country(ms, hz, keep_beyond = TRUE),
  chr_month = censor_sim_by_country(
    dplyr::mutate(ms, month = as.character(month)),
    hz
  ),
  factor_gid = censor_sim_by_country(dplyr::mutate(ms, gid = factor(gid)), hz),
  events = censor_sim_by_country(
    ev,
    hz,
    month_col = "start_month",
    kind = "events",
    keep_beyond = TRUE
  ),
  country = censor_sim_by_country(
    dplyr::mutate(ms, iso = substr(gid, 1, 3)),
    hz,
    country_col = "iso"
  ),
  empty = censor_sim_by_country(ms[0, ], hz),
  err_miss = err_of(censor_sim_by_country(
    dplyr::bind_rows(ms, tibble::tibble(gid = "UGA_1_1", month = months[1])),
    hz
  )),
  err_na = err_of(censor_sim_by_country(
    dplyr::bind_rows(
      ms,
      tibble::tibble(gid = NA_character_, month = months[1])
    ),
    hz
  ))
)
got <- list(pop = pop, censor = cen)

if (cli$capture) {
  saveRDS(got, GOLDEN_PATH)
  cat("captured", GOLDEN_PATH, "\n")
  quit(status = 0)
}
gold <- readRDS(GOLDEN_PATH)
for (k in names(gold$pop)) {
  expect(
    identical(got$pop[[k]], gold$pop[[k]]),
    paste0("resolve_gid_population / rescale: ", k, " identical to golden")
  )
}
for (k in names(gold$censor)) {
  expect(
    identical(got$censor[[k]], gold$censor[[k]]),
    paste0("censor_sim_by_country: ", k, " identical to golden")
  )
}
expect(
  all(
    c("direct", "wpp_rescale", "parent_share", "missing") %in%
      gold$pop$na_wpp$pop_tier
  ),
  "fixture reaches all four tiers"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
