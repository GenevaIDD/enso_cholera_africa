# test_enso_leads.R -- ENSO lead columns for the -11..11 single-lag sweep.
#
# Soft gate (see analysis/tests/run_all.R). No fit, no data, no MCMC.
#   * joinEnsoLags(max_lead = 0) writes exactly enso_lag0..max_lag (prior output).
#   * joinEnsoLags(max_lead = 2) adds enso_lead1/2 = the ONI of month m+1/m+2.
#   * a lead month past the end of the ENSO series is an error, not an NA (the
#     single-lag fit would read an NA as neutral ENSO).
#
# Run: Rscript analysis/tests/test_enso_leads.R

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R")
source("analysis/utils/enso.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# ENSO lead columns (joinEnsoLags max_lead)\n")

months_seq <- seq(as.Date("2020-01-01"), as.Date("2021-12-01"), by = "month")
enso_dat <- tibble::tibble(year = year(months_seq), month = month(months_seq),
                           enso_index = seq_along(months_seq) / 10)
oni_at <- function(d) enso_dat$enso_index[enso_dat$year == year(d) & enso_dat$month == month(d)]

df <- tibble::tibble(t_start = as.Date(c("2020-06-01", "2021-03-01")))

j0 <- joinEnsoLags(df, enso_dat, max_lag = 3L)
expect(identical(names(j0), c("t_start", sprintf("enso_lag%d", 0:3))),
       "max_lead = 0 writes exactly enso_lag0..max_lag")

j2 <- joinEnsoLags(df, enso_dat, max_lag = 3L, max_lead = 2L)
expect(all(c("enso_lead1", "enso_lead2") %in% names(j2)), "max_lead = 2 adds enso_lead1/2")
expect(isTRUE(all.equal(j2$enso_lead2, vapply(df$t_start %m+% months(2), oni_at, numeric(1)))),
       "enso_lead2 is the ONI two months after t_start")
expect(isTRUE(all.equal(j2$enso_lag3, vapply(df$t_start %m-% months(3), oni_at, numeric(1)))),
       "enso_lag3 unchanged by the leads")
expect(isTRUE(all.equal(j2[names(j0)], j0)), "lag columns identical with and without leads")

late <- tibble::tibble(t_start = as.Date("2021-11-01"))
err <- tryCatch({ joinEnsoLags(late, enso_dat, max_lag = 1L, max_lead = 2L); NULL },
                error = function(e) conditionMessage(e))
expect(!is.null(err) && grepl("past the end of the ENSO series", err),
       "a lead past the series end stops instead of writing NA")


# addEnsoLeadsFromLag0() (sim panels)
pn <- tibble::tibble(gid = rep(c("a", "b"), each = 4),
                     t_start = rep(seq(as.Date("2021-01-01"), by = "month", length.out = 4), 2),
                     enso_lag0 = rep(c(0.1, 0.2, 0.3, 0.4), 2))
pl <- addEnsoLeadsFromLag0(pn, 2L)
expect(isTRUE(all.equal(pl$enso_lead1[pl$gid == "a"], c(0.2, 0.3, 0.4, NA))), "enso_lead1(t) = enso_lag0(t + 1); NA past the panel")
expect(isTRUE(all.equal(pl$enso_lead2[pl$gid == "b"], c(0.3, 0.4, NA, NA))), "enso_lead2 per gid")
ex <- tibble::tibble(gid = "a", t_start = seq(as.Date("2021-04-01"), by = "month", length.out = 3), enso_lag0 = c(9, 0.5, 0.6))
pe <- addEnsoLeadsFromLag0(pn, 2L, extra = list(ex))
expect(isTRUE(all.equal(pe$enso_lead2[pe$gid == "a"], c(0.3, 0.4, 0.5, 0.6))),
       "extra panel extends the series after the panel's last month (its overlapping month ignored)")
expect(identical(addEnsoLeadsFromLag0(pn, 0L), pn), "max_lead = 0 is a no-op")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
