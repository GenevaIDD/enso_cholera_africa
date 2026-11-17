# test_nbr_naming.R — neighbour-contagion covariate naming helpers
#   (paths.R).
#
# Usage:
#   Rscript analysis/tests/test_nbr_naming.R
#
# Details:
#   Run from the project root. Pins the single-source-of-truth naming used by
#   prep, the formula builder, artefact tokens and the simulator kernels:
#   nbr_covariate_name() / nbr_covariate_pattern() / parse_nbr_covariate() /
#   nbr_artefact_token(). Also asserts that simulator_core.R's source-time
#   literal in .LP_CONSUMED_PATTERNS stays identical to
#   nbr_covariate_pattern().

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("analysis/utils/base.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/simulator_core.R")
source("analysis/tests/helpers.R")
PASS <- 0L
FAIL <- 0L

cat("\n## Names (defaults byte-identical to the pre-helper sprintf)\n")
expect(identical(nbr_covariate_name(500L), "nbr_lag1_r500"), "raw default")
expect(
  identical(nbr_covariate_name(500L, "raw", "none"), "nbr_lag1_r500"),
  "raw explicit"
)
expect(
  identical(nbr_covariate_name(200, "rownorm"), "nbr_lag1_r200_rn"),
  "rownorm"
)
expect(
  identical(
    nbr_covariate_name(500L, "rownorm", "group"),
    "nbr_lag1_r500_rn_pg"
  ),
  "rownorm + group"
)
expect(
  identical(nbr_covariate_name(500L, "raw", "group"), "nbr_lag1_r500_pg"),
  "raw + group"
)
expect(
  identical(nbr_covariate_name(500L, NULL, NULL), "nbr_lag1_r500"),
  "NULL scale/partition"
)

cat("\n## Tokens\n")
expect(identical(nbr_artefact_token(500L), "_nbr500"), "token raw")
expect(
  identical(nbr_artefact_token(500L, "rownorm"), "_nbr500n"),
  "token rownorm"
)
expect(
  identical(nbr_artefact_token(500L, "rownorm", "group"), "_nbr500ng"),
  "token rownorm+group"
)
expect(
  identical(nbr_artefact_token(500L, "raw", "group"), "_nbr500g"),
  "token raw+group"
)

cat("\n## Pattern + parse round-trip\n")
for (sc in c("raw", "rownorm")) {
  for (pt in c("none", "group")) {
    for (rho in c(200L, 500L)) {
      nm <- nbr_covariate_name(rho, sc, pt)
      pr <- parse_nbr_covariate(nm)
      expect(
        grepl(nbr_covariate_pattern(), nm),
        sprintf("pattern matches %s", nm)
      )
      expect(
        identical(
          pr[c("rho_km", "scale", "partition")],
          list(rho_km = rho, scale = sc, partition = pt)
        ) &&
          is.null(pr$partition_token),
        sprintf("parse round-trip %s", nm)
      )
    }
  }
}
for (bad in c(
  "nbr_lag1",
  "nbr_lag1_r500_x",
  "ongoing_outbreak_lag1",
  "nbr_lag1_r500_pg_rn"
)) {
  expect(
    !grepl(nbr_covariate_pattern(), bad),
    sprintf("pattern rejects %s", bad)
  )
  expect(is.null(parse_nbr_covariate(bad)), sprintf("parse NULL for %s", bad))
}

cat("\n## simulator_core literal twin\n")
expect(
  nbr_covariate_pattern() %in% .LP_CONSUMED_PATTERNS,
  ".LP_CONSUMED_PATTERNS carries nbr_covariate_pattern() verbatim"
)

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
