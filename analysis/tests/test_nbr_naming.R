# test_nbr_naming.R -- neighbour-contagion covariate naming helpers (paths.R).
#
# HARD gate for the single-source-of-truth naming convention used by prep,
# the formula builder, artefact tokens and the three simulator kernels:
#   nbrCovariateName() / nbrCovariatePattern() / parseNbrCovariate() /
#   nbrArtefactToken()
# Also asserts that simulator_core.R's source-time literal in
# .LP_CONSUMED_PATTERNS stays identical to nbrCovariatePattern().
#
# Run: Rscript analysis/tests/test_nbr_naming.R   (from the project root)

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
source("analysis/utils/base.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/simulator_core.R")
source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L

cat("\n## Names (defaults byte-identical to the pre-helper sprintf)\n")
expect(identical(nbrCovariateName(500L), "nbr_lag1_r500"), "raw default")
expect(identical(nbrCovariateName(500L, "raw", "none"), "nbr_lag1_r500"), "raw explicit")
expect(identical(nbrCovariateName(200, "rownorm"), "nbr_lag1_r200_rn"), "rownorm")
expect(identical(nbrCovariateName(500L, "rownorm", "group"), "nbr_lag1_r500_rn_pg"), "rownorm + group")
expect(identical(nbrCovariateName(500L, "raw", "group"), "nbr_lag1_r500_pg"), "raw + group")
expect(identical(nbrCovariateName(500L, NULL, NULL), "nbr_lag1_r500"), "NULL scale/partition")

cat("\n## Tokens\n")
expect(identical(nbrArtefactToken(500L), "_nbr500"), "token raw")
expect(identical(nbrArtefactToken(500L, "rownorm"), "_nbr500n"), "token rownorm")
expect(identical(nbrArtefactToken(500L, "rownorm", "group"), "_nbr500ng"), "token rownorm+group")
expect(identical(nbrArtefactToken(500L, "raw", "group"), "_nbr500g"), "token raw+group")

cat("\n## Pattern + parse round-trip\n")
for (sc in c("raw", "rownorm")) for (pt in c("none", "group")) for (rho in c(200L, 500L)) {
  nm <- nbrCovariateName(rho, sc, pt)
  pr <- parseNbrCovariate(nm)
  expect(grepl(nbrCovariatePattern(), nm), sprintf("pattern matches %s", nm))
  expect(identical(pr[c("rho_km", "scale", "partition")],
                   list(rho_km = rho, scale = sc, partition = pt)) && is.null(pr$partition_token),
         sprintf("parse round-trip %s", nm))
}
for (bad in c("nbr_lag1", "nbr_lag1_r500_x", "ongoing_outbreak_lag1", "nbr_lag1_r500_pg_rn")) {
  expect(!grepl(nbrCovariatePattern(), bad), sprintf("pattern rejects %s", bad))
  expect(is.null(parseNbrCovariate(bad)), sprintf("parse NULL for %s", bad))
}

cat("\n## simulator_core literal twin\n")
expect(nbrCovariatePattern() %in% .LP_CONSUMED_PATTERNS,
       ".LP_CONSUMED_PATTERNS carries nbrCovariatePattern() verbatim")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
