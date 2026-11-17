# test_nbr_partition_token.R — HARD: neighbour contagion over an arbitrary
#   registry.
#
# Usage:
#   Rscript analysis/tests/test_nbr_partition_token.R
#
# Details:
#   No data, no fit. Pins the _part<token> covariate names and their parse
#   round trip, token -> registry uniqueness, prep columns == the simulator
#   kernel for that registry, resolve_nbr_spec() and the formula / artefact
#   names under contagion.partition_source.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/paths.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatial_contagion.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

cat("# registry-partitioned neighbour contagion\n")

# [1] naming: nbr_lag1_r<rho>[_rn]_part<token>; _pg and unpartitioned names
#     byte-identical; artefact token _part-<token>; empty source refused.
cat("\n[1] naming\n")
expect(
  identical(
    nbr_covariate_name(500L, "rownorm", "registry", "smoke_grouping"),
    "nbr_lag1_r500_rn_partsmokegrouping"
  ),
  "rownorm + registry"
)
expect(
  identical(
    nbr_covariate_name(200L, "raw", "registry", "enso_grouping"),
    "nbr_lag1_r200_partensogrouping"
  ),
  "raw + registry"
)
expect(
  identical(
    nbr_covariate_name(500L, "rownorm", "group"),
    "nbr_lag1_r500_rn_pg"
  ),
  "_pg unchanged"
)
expect(identical(nbr_covariate_name(500L), "nbr_lag1_r500"), "plain unchanged")
expect(
  identical(
    nbr_artefact_token(500L, "rownorm", "registry", "smoke_grouping"),
    "_nbr500n_part-smokegrouping"
  ),
  "artefact token"
)
expect(
  identical(nbr_artefact_token(500L, "rownorm", "group"), "_nbr500ng"),
  "group artefact token unchanged"
)
expect(
  !is.null(err_of(nbr_covariate_name(500L, "raw", "registry", ""))),
  "registry without source refused"
)

# [2] parse round trip for all three partition kinds; pattern and the
#     simulator_core literal twin agree; malformed names rejected.
cat("\n[2] parse round trip\n")
for (sc in c("raw", "rownorm")) {
  for (pt in c("none", "group", "registry")) {
    nm <- nbr_covariate_name(
      500L,
      sc,
      pt,
      if (pt == "registry") "hmm_clusters" else NULL
    )
    p <- parse_nbr_covariate(nm)
    expect(
      grepl(nbr_covariate_pattern(), nm) &&
        identical(p$scale, sc) &&
        identical(p$partition, pt) &&
        identical(
          p$partition_token,
          if (pt == "registry") "hmmclusters" else NULL
        ),
      sprintf("round trip %s", nm)
    )
  }
}
for (bad in c(
  "nbr_lag1_r500_part",
  "nbr_lag1_r500_partX",
  "nbr_lag1_r500_pg_parta",
  "nbr_lag1_r500_part_a",
  "nbr_lag1_r500_parta_rn"
)) {
  expect(
    !grepl(nbr_covariate_pattern(), bad) && is.null(parse_nbr_covariate(bad)),
    sprintf("rejects %s", bad)
  )
}
expect(
  nbr_covariate_pattern() %in% .LP_CONSUMED_PATTERNS,
  "simulator_core literal twin updated"
)

# [3] token -> registry: unique for every shipped registry.
cat("\n[3] token -> registry\n")
srcs <- sub(
  "\\.csv$",
  "",
  list.files("analysis/groupings", pattern = "\\.csv$")
)
for (s in srcs) {
  tok <- parse_nbr_covariate(nbr_covariate_name(
    500L,
    "raw",
    "registry",
    s
  ))$partition_token
  expect(
    identical(grouping_source_from_token(tok), s),
    sprintf("%s -> %s -> %s", s, tok, s)
  )
}

# [4] build_contagion_table(partitions=): _part columns equal the kernel the
#     simulator builds for that registry; _pg / plain columns unchanged.
cat("\n[4] prep columns == simulator kernel\n")
ids <- c("ETH_1", "ETH_2", "KEN_1", "NGA_1", "CMR_1")
D <- matrix(
  c(
    0,
    100,
    300,
    3000,
    3100,
    100,
    0,
    250,
    2900,
    3000,
    300,
    250,
    0,
    3200,
    3300,
    3000,
    2900,
    3200,
    0,
    400,
    3100,
    3000,
    3300,
    400,
    0
  ),
  5,
  5,
  dimnames = list(ids, ids)
)
diag(D) <- Inf
months <- as.Date(c("2015-01-01", "2015-02-01"))
active <- tidyr::expand_grid(ADM1 = ids, t_start = months) %>%
  dplyr::mutate(
    ongoing_outbreak_lag1 = ADM1 %in% c("KEN_1", "NGA_1") & t_start == months[1]
  )
sg <- load_grouping("smoke_grouping")
tb <- build_contagion_table(
  active,
  D,
  rhos = 500L,
  partition = load_grouping("enso_grouping"),
  partitions = list(smoke_grouping = sg)
)
tb0 <- build_contagion_table(
  active,
  D,
  rhos = 500L,
  partition = load_grouping("enso_grouping")
)
expect(
  all(
    c(
      "nbr_lag1_r500_partsmokegrouping",
      "nbr_lag1_r500_rn_partsmokegrouping"
    ) %in%
      names(tb)
  ),
  "_part columns written for the listed registry"
)
expect(
  isTRUE(all.equal(tb[names(tb0)], tb0)),
  "plain and _pg columns unchanged by partitions="
)
K <- build_nbr_kernel(D, 500, scale = "rownorm", partition = sg)
# pressure at month 1 = K %*% the (already lag-1) activity indicator of month 1
a1 <- as.numeric(ids %in% c("KEN_1", "NGA_1"))
want <- stats::setNames(as.numeric(K %*% a1), ids)
got <- tb %>%
  dplyr::filter(t_start == months[1]) %>%
  dplyr::arrange(match(ADM1, ids))
expect(
  isTRUE(all.equal(
    unname(got$nbr_lag1_r500_rn_partsmokegrouping),
    unname(want)
  )),
  "_rn_part column = rownorm(registry kernel) %*% lag-1 activity (the simulator kernel)"
)
expect(
  got$nbr_lag1_r500_rn_partsmokegrouping[got$ADM1 == "ETH_1"] > 0 &&
    got$nbr_lag1_r500_rn_pg[got$ADM1 == "ETH_1"] == 0,
  "ETH sees KEN through smoke_grouping (grp_a) but not through enso_grouping (horn vs east)"
)

# [5] resolve_nbr_spec(): a _part term loads its registry from the token;
#     a configured partition_source that disagrees is refused; nbr_cfg_from_sim()
#     defaults the contagion registry to filter.partition_source.
cat("\n[5] resolve_nbr_spec\n")
fe_draws <- function(col) list(fe = matrix(0, 1, 1, dimnames = list(NULL, col)))
sp <- resolve_nbr_spec(
  list(enable = TRUE, rho_km = 500L),
  fe_draws("nbr_lag1_r500_rn_partsmokegrouping")
)
expect(
  isTRUE(sp$use) &&
    identical(sp$partition, "registry") &&
    identical(sp$partition_source, "smoke_grouping") &&
    identical(sp$registry$group, sg$group),
  "registry loaded from the fit's token"
)
expect(
  !is.null(err_of(resolve_nbr_spec(
    list(enable = TRUE, rho_km = 500L, partition_source = "hmm_clusters"),
    fe_draws("nbr_lag1_r500_rn_partsmokegrouping")
  ))),
  "configured partition_source disagreeing with the token refused"
)
cfg <- nbr_cfg_from_sim(list(
  nbr_contagion = list(enable = TRUE),
  filter = list(partition_source = "smoke_grouping"),
  grouping = list(source = "enso_grouping")
))
expect(
  identical(cfg$partition_source, "smoke_grouping") &&
    identical(cfg$grouping_source, "enso_grouping"),
  "nbr_cfg_from_sim(): contagion registry defaults to filter.partition_source"
)
sp_pg <- resolve_nbr_spec(cfg, fe_draws("nbr_lag1_r500_rn_pg"))
expect(
  identical(sp_pg$partition, "group") &&
    identical(
      unique(sp_pg$registry$group),
      unique(load_grouping("enso_grouping")$group)
    ),
  "a _pg term still reads the fit grouping"
)

# [6] formula + artefact naming follow contagion.partition_source.
cat("\n[6] formula + artefact names\n")
source("analysis/utils/formula.R")
opt <- get_default_options()
opt$occurrence$covariates$nbr_contagion_lag1 <- TRUE
opt$contagion$scale <- "rownorm"
opt$contagion$partition <- "registry"
opt$contagion$partition_source <- "smoke_grouping"
ft <- tryCatch(build_formula_terms(opt), error = function(e) NULL)
expect(
  !is.null(ft) &&
    any(grepl("nbr_lag1_r500_rn_partsmokegrouping", unlist(ft), fixed = TRUE)),
  "build_formula_terms() selects the _part column"
)
expect(
  grepl("_nbr500n_part-smokegrouping", make_run_file(opt, "fit"), fixed = TRUE),
  "occurrence artefact name carries _part-<token>"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
