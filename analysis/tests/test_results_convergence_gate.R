# test_results_convergence_gate.R — HARD: report_convergence v2 (table D27 /
#   S8), each deployed fit's deploy-gate verdict.
#
# Usage:
#   Rscript analysis/tests/test_results_convergence_gate.R
#
# Details:
#   Pure, synthetic: no fit, no artefact; quick. Pins the gate rows built
#   from a convergence artefact's stored tables and the deploy's
#   convergence block.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/base.R")
source("analysis/utils/deploy_core.R") # config, paths, diagnostics, deploy_convergence_cfg()
source("analysis/utils/results_core.R") # TRACK_CODE
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
env <- new.env()
sys.source("analysis/utils/results_derive.R", envir = env, keep.source = FALSE)

mk_d <- function(rhat_sd = 1.002, ess_tail_b = 1000, div = c(0L, 0L)) {
  list(
    summary = tibble::tibble(
      variable = c(
        "b_Intercept",
        "sd_country__Intercept",
        "r_country[KEN,Intercept]"
      ),
      rhat = c(1.001, rhat_sd, 1.004),
      ess_bulk = 1000,
      ess_tail = c(ess_tail_b, 1000, 1000)
    ),
    sampler = tibble::tibble(
      chain = seq_along(div),
      divergent = div,
      treedepth_hits = 0L,
      ebfmi = 0.9
    )
  )
}
# A convergence artefact of one track as derive_convergence() keys it.
ds <- list(
  enso_east = mk_d(),
  enso_west = mk_d(rhat_sd = 1.02),
  enso_horn = mk_d(ess_tail_b = 300)
)
lags <- c(enso_east = 0L, enso_west = -6L, enso_horn = 10L)
key <- function(x, g) {
  dplyr::mutate(
    tibble::as_tibble(x),
    track = "occurrence",
    enso_group = g,
    lag = lags[[g]],
    .before = 1
  )
}
conv <- list(
  scalars = bind_rows(lapply(names(ds), function(g) {
    key(tibble::tibble(n_chains = 2L), g)
  })),
  summary = bind_rows(lapply(names(ds), function(g) key(ds[[g]]$summary, g))),
  sampler = bind_rows(lapply(names(ds), function(g) key(ds[[g]]$sampler, g)))
)

# [1] .convergence_gate_rows(): from the keyed summary / sampler tables a
#     convergence artefact stores, the verdict equals
#     diagnostics.R::fit_diagnostics_gate() on that fit's fit_diagnostics()
#     result; pass / fail / allowed (deploy.convergence.allow ids)
cat("# [1] gate rows from the stored tables\n")
gc0 <- list(
  cfg = deploy_convergence_cfg(list(), validate_allow = FALSE),
  dep = NULL
)
r <- env$.convergence_gate_rows(conv, gc0)
expect(identical(r$enso_group, names(ds)), "one row per fit, in order")
same <- all(vapply(
  names(ds),
  function(g) {
    dg <- fit_diagnostics_gate(ds[[g]])
    rg <- r[r$enso_group == g, ]
    identical(rg$gate_reasons, dg$reasons) &&
      identical(rg$gate_max_rhat_scope, dg$max_rhat_scope) &&
      identical(rg$gate_max_rhat_variable, dg$max_rhat_variable)
  },
  logical(1)
))
expect(same, "verdict = fit_diagnostics_gate() on each fit's own diagnostics")
expect(
  identical(r$gate_verdict, c("pass", "fail", "fail")),
  "pass / fail (R-hat) / fail (tail ESS)"
)
expect(
  identical(r$gate_reasons[2], "rhat") &&
    identical(r$gate_reasons[3], "ess_tail"),
  "reasons"
)
gc1 <- gc0
gc1$cfg$allow <- "occ:west"
expect(
  identical(
    env$.convergence_gate_rows(conv, gc1)$gate_verdict,
    c("pass", "allowed", "fail")
  ),
  "an allow-listed failing fit is `allowed`"
)
gc2 <- gc0
gc2$cfg$scope <- "all"
expect(
  identical(env$.convergence_gate_rows(conv, gc2)$gate_scope, rep("all", 3L)),
  "scope from the deploy block"
)

# [2] .convergence_gate_cfg(): no deploy -> the in-function defaults (off,
#     population); the shipped prod_run.yml -> warn at rhat_max 1.05; a
#     fixture copy with refuse and default thresholds -> refuse
cat("# [2] the deploy's convergence block\n")
g0 <- env$.convergence_gate_cfg(list(RES = list(deploy_config = "")))
expect(
  is.null(g0$dep) &&
    identical(g0$cfg$action, "off") &&
    identical(g0$cfg$scope, "population"),
  "no deploy: action off, scope population"
)
# Fixture: prod_run.yml with convergence.action refuse and the default R-hat.
dep_refuse <- tempfile("conv_refuse_", fileext = ".yml")
on.exit(unlink(dep_refuse), add = TRUE)
yr <- yaml::read_yaml("analysis/configs/deploy/prod_run.yml")
yr$deploy$convergence$action <- "refuse"
yr$deploy$convergence$rhat_max <- "default"
yaml::write_yaml(yr, dep_refuse)
gr <- env$.convergence_gate_cfg(list(RES = list(deploy_config = dep_refuse)))
expect(
  identical(gr$cfg$action, "refuse") &&
    identical(gr$cfg$rhat_max, unname(DIAG_THRESHOLDS[["rhat"]])),
  "action refuse, rhat_max default: refuse at the default thresholds"
)
gw <- env$.convergence_gate_cfg(list(
  RES = list(deploy_config = "analysis/configs/deploy/prod_run.yml")
))
expect(
  identical(gw$cfg$action, "warn") && identical(gw$cfg$rhat_max, 1.05),
  "prod_run.yml: warn, rhat_max 1.05"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
