# test_deploy_convergence_gate.R — HARD: the convergence gate of
#   `deploy_prod_run.R --stage deploy_fits` and of `--stage gen_grouping`.
#
# Usage:
#   Rscript analysis/tests/test_deploy_convergence_gate.R
#
# Details:
#   Covers diagnostics.R::fit_diagnostics_gate(), the one source of the
#   thresholds, deploy_core.R deploy_convergence_cfg() /
#   deploy_check_convergence() / deploy_collect_fits() and
#   kl_grouping.R::kl_grouping_convergence(). No sampling, no data; quick-safe.

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/deploy_core.R")
for (m in c("dlnm", "dlnm_kernel", "kl_grouping")) {
  source(sprintf("analysis/utils/%s.R", m))
}
source("analysis/tests/helpers_synth.R")

errs <- function(expr) inherits(tryCatch(expr, error = identity), "error")
err_msg <- function(expr) {
  tryCatch(
    {
      expr
      ""
    },
    error = conditionMessage
  )
}
warns <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  w
}

# A synthetic fit_diagnostics() result: one row per parameter class, sampler
# table of `n_ch` chains.
mk_d <- function(
  rhat = c(b = 1.001, sd = 1.002, sds = 1.003, r = 1.004, z = 1.004),
  ess_bulk = 1000,
  ess_tail = 1000,
  div = c(0L, 0L, 0L, 0L),
  ebfmi = c(0.9, 0.9, 0.9, 0.9)
) {
  vars <- c(
    b = "b_Intercept",
    sd = "sd_country__Intercept",
    sds = "sds_st_start_seas_1",
    r = "r_country[KEN,Intercept]",
    z = "z_1[1,1]"
  )
  rh <- rep(1.001, length(vars))
  names(rh) <- names(vars)
  rh[names(rhat)] <- rhat
  eb <- rep(1000, length(vars))
  names(eb) <- names(vars)
  et <- eb
  if (!is.null(names(ess_bulk))) {
    eb[names(ess_bulk)] <- ess_bulk
  } else {
    eb[] <- ess_bulk
  }
  if (!is.null(names(ess_tail))) {
    et[names(ess_tail)] <- ess_tail
  } else {
    et[] <- ess_tail
  }
  list(
    summary = tibble::tibble(
      variable = unname(vars),
      rhat = unname(rh),
      ess_bulk = unname(eb),
      ess_tail = unname(et)
    ),
    sampler = tibble::tibble(
      chain = seq_along(div),
      divergent = as.integer(div),
      treedepth_hits = 0L,
      ebfmi = ebfmi
    )
  )
}
gate <- function(d, ...) fit_diagnostics_gate(d, ...)

# On synthetic fit_diagnostics() lists: scope population vs all, R-hat / bulk /
# tail ESS, the gross R-hat bound on every parameter, divergences, E-BFMI,
# all-NA scope.
cat("# [1] fit_diagnostics_gate()\n")
g <- gate(mk_d())
expect(isTRUE(g$pass) && identical(g$reasons, ""), "a clean fit passes")
expect(
  identical(
    names(g),
    c(
      "scope",
      "n_params_scope",
      "max_rhat_scope",
      "max_rhat_variable",
      "min_ess_bulk_scope",
      "min_ess_tail_scope",
      "max_rhat_all",
      "divergences",
      "treedepth_hits",
      "min_ebfmi",
      "pass",
      "reasons"
    )
  ) &&
    nrow(g) == 1L,
  "one row, columns in contract order"
)
expect(
  g$n_params_scope == 3L,
  "population scope drops r_ / z_ (3 of 5 parameters)"
)
g <- gate(mk_d(rhat = c(sd = 1.02)))
expect(
  !g$pass &&
    identical(g$reasons, "rhat") &&
    g$max_rhat_variable == "sd_country__Intercept",
  "sd_ at R-hat 1.02 fails with reason rhat, naming the parameter"
)
expect(!gate(mk_d(rhat = c(sds = 1.0206)))$pass, "sds_ at 1.0206 fails")
expect(
  isTRUE(gate(mk_d(rhat = c(r = 1.02)))$pass),
  "r_ at 1.02 passes under population"
)
g <- gate(mk_d(rhat = c(r = 1.02)), scope = "all")
expect(
  !g$pass && identical(g$reasons, "rhat") && g$n_params_scope == 5L,
  "r_ at 1.02 fails under all"
)
expect(
  isTRUE(gate(mk_d(rhat = c(z = 1.02)))$pass),
  "z_ at 1.02 passes under population"
)
g <- gate(mk_d(rhat = c(r = 1.06)))
expect(
  !g$pass && identical(g$reasons, "rhat_any") && g$max_rhat_all == 1.06,
  "r_ at 1.06 fails under population (rhat_any, every parameter)"
)
g <- gate(mk_d(ess_tail = c(b = 300)))
expect(
  !g$pass && identical(g$reasons, "ess_tail"),
  "tail ESS 300 fails (ess_tail)"
)
g <- gate(mk_d(ess_bulk = c(sd = 361), ess_tail = c(sd = 231)))
expect(
  !g$pass && identical(g$reasons, "ess_bulk;ess_tail"),
  "bulk 361 + tail 231 fail together"
)
expect(
  isTRUE(gate(mk_d(ess_bulk = c(r = 50)))$pass),
  "low ESS on r_ is out of the population scope"
)
g <- gate(mk_d(div = c(1L, 0L, 0L, 0L)))
expect(
  !g$pass && identical(g$reasons, "divergences") && g$divergences == 1L,
  "1 divergence fails at max_divergent 0"
)
expect(
  isTRUE(gate(mk_d(div = c(1L, 0L, 0L, 0L)), max_divergent = 1L)$pass),
  "1 divergence passes at max_divergent 1"
)
g <- gate(mk_d(ebfmi = c(0.9, 0.2, 0.9, 0.9)))
expect(!g$pass && identical(g$reasons, "ebfmi"), "E-BFMI 0.2 fails")
d_td <- mk_d()
d_td$sampler$treedepth_hits <- c(5L, 0L, 0L, 0L)
g <- gate(d_td)
expect(
  isTRUE(g$pass) && g$treedepth_hits == 5L,
  "treedepth hits reported, not gated"
)
d_na <- mk_d()
d_na$summary$rhat <- NA_real_
g <- gate(d_na)
expect(!g$pass && grepl("na", g$reasons), "an all-NA scope fails (na)")
d_e <- mk_d()
d_e$summary <- d_e$summary[grepl("^(r_|z_)", d_e$summary$variable), ]
g <- gate(d_e)
expect(
  !g$pass && g$n_params_scope == 0L && grepl("na", g$reasons),
  "an empty scope fails (na)"
)
g <- gate(mk_d(rhat = c(sd = 1.02)), rhat_max = 1.03)
expect(isTRUE(g$pass), "rhat_max is honoured")
expect(
  errs(diag_scope_summary(mk_d()$summary, "fixed")),
  "unknown scope is an error"
)
expect(
  identical(diag_scope_summary(mk_d()$summary, "all"), mk_d()$summary),
  "scope all keeps every row"
)

# The "default" sentinel resolves to DIAG_THRESHOLDS; no 1.01 / 400 literal in
# deploy_core.R.
cat("# [2] one source of the thresholds\n")
dep0 <- list(
  grouping = "enso_grouping",
  fits = list(groups = c("east", "west"))
)
c0 <- deploy_convergence_cfg(dep0)
expect(
  identical(c0$rhat_max, DIAG_THRESHOLDS[["rhat"]]) &&
    identical(c0$ess_min, DIAG_THRESHOLDS[["ess"]]),
  "missing thresholds resolve to DIAG_THRESHOLDS"
)
c1 <- deploy_convergence_cfg(c(
  dep0,
  list(convergence = list(rhat_max = "default", ess_min = "default"))
))
expect(
  identical(c1$rhat_max, DIAG_THRESHOLDS[["rhat"]]) &&
    identical(c1$ess_min, DIAG_THRESHOLDS[["ess"]]),
  "\"default\" sentinels resolve to DIAG_THRESHOLDS"
)
expect(
  identical(c0$rhat_any, DIAG_GATE[["rhat_any"]]) &&
    identical(c0$ebfmi_min, DIAG_GATE[["ebfmi"]]),
  "gross bounds come from DIAG_GATE"
)
src <- readLines("analysis/utils/deploy_core.R")
code <- sub("#.*$", "", src)
expect(
  !any(grepl(
    "(^|[^0-9.])1\\.01([^0-9]|$)|(^|[^0-9.])400([^0-9]|$)|1\\.05|0\\.3([^0-9]|$)",
    code
  )),
  "no 1.01 / 400 / 1.05 / 0.3 literal in deploy_core.R code"
)

# Defaults (a missing block = off), bad action / scope / threshold /
# max_divergent / allow id.
cat("# [3] deploy_convergence_cfg()\n")
expect(
  identical(c0$action, "off") &&
    identical(c0$scope, "population") &&
    identical(c0$max_divergent, 0L) &&
    identical(c0$allow, character(0)),
  "a missing block: off, population, 0, no allow"
)
expect(
  identical(
    c0,
    deploy_convergence_cfg(c(
      dep0,
      list(
        convergence = list(
          action = "off",
          scope = "population",
          rhat_max = "default",
          ess_min = "default",
          max_divergent = 0L,
          allow = list()
        )
      )
    ))
  ),
  "the restated default block = the missing block"
)
cfg_of <- function(...) {
  deploy_convergence_cfg(c(dep0, list(convergence = list(...))))
}
expect(errs(cfg_of(action = "stop")), "bad action is an error")
expect(errs(cfg_of(scope = "fixed")), "bad scope is an error")
expect(errs(cfg_of(rhat_max = "lax")), "non-numeric rhat_max is an error")
expect(errs(cfg_of(max_divergent = -1)), "negative max_divergent is an error")
expect(
  errs(cfg_of(max_divergent = 1.5)),
  "fractional max_divergent is an error"
)
expect(
  grepl(
    "occ:north",
    err_msg(cfg_of(action = "refuse", allow = list("occ:north")))
  ),
  "unknown allow id is an error naming it"
)
expect(
  errs(cfg_of(action = "refuse", allow = list("siz:east"))),
  "allow id of a retired track is an error"
)
expect(
  errs(cfg_of(action = "refuse", allow = list("occ:enso_east"))),
  "allow ids take the short group"
)
ca <- cfg_of(
  action = "refuse",
  rhat_max = 1.02,
  ess_min = 300,
  max_divergent = 2,
  allow = list("occ:west", "mon:east")
)
expect(
  identical(ca$allow, c("occ:west", "mon:east")) &&
    ca$rhat_max == 1.02 &&
    ca$ess_min == 300 &&
    identical(ca$max_divergent, 2L),
  "explicit values are taken"
)
expect(
  identical(
    deploy_convergence_cfg(
      c(
        list(grouping = "enso_grouping_absent", fits = list(groups = "x")),
        list(convergence = list(action = "warn", allow = list("occ:x")))
      ),
      validate_allow = FALSE
    )$allow,
    "occ:x"
  ),
  "validate_allow = FALSE skips the id check (gen_grouping)"
)

# With injected read / diagnostics functions: off reads nothing, warn warns,
# refuse stops naming track:group, allow-listed fits pass as `allowed`,
# --dry-run never stops.
cat("# [4] deploy_check_convergence()\n")
plan <- data.frame(
  track = c("occ", "occ", "dur", "dur", "occ"),
  group = c("ensoeast", "ensowest", "ensoeast", "ensowest", "ensoeast"),
  lag = c(0L, -6L, 1L, -4L, 0L),
  kind = c("fit", "fit", "dur_fit", "dur_fit", "loo"),
  src = c("s/occ_e", "s/occ_w", "s/dur_e", "s/dur_w", "s/occ_e_loo"),
  file = c("occ_e.rds", "occ_w.rds", "dur_e.rds", "dur_w.rds", "occ_e_loo.rds"),
  is_fit = c(TRUE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
diag_by_file <- list(
  occ_e.rds = mk_d(),
  occ_w.rds = mk_d(rhat = c(sds = 1.0206)),
  dur_e.rds = mk_d(),
  dur_w.rds = mk_d(rhat = c(sd = 1.0166))
)
dep_g <- function(...) {
  c(
    dep0,
    list(
      fits = list(groups = c("east", "west"), dest = "dest"),
      convergence = list(...)
    )
  )
}
stop_fn <- function(...) stop("must not be called")
expect(
  is.null(deploy_check_convergence(
    plan,
    dep_g(action = "off"),
    read_fn = stop_fn,
    diag_fn = stop_fn
  )),
  "off returns NULL without reading or diagnosing"
)
expect(
  is.null(deploy_check_convergence(
    plan,
    dep0,
    read_fn = stop_fn,
    diag_fn = stop_fn
  )),
  "a missing block is off (nothing read)"
)
w <- warns(
  r <- deploy_check_convergence(
    plan,
    dep_g(action = "warn"),
    diags = diag_by_file
  )
)
expect(
  length(w) == 1L && grepl("occ:west", w) && grepl("dur:west", w),
  "warn warns naming the failing fits"
)
expect(
  nrow(r$table) == 4L &&
    identical(r$failed, c("occ:west", "dur:west")) &&
    identical(r$table$verdict, c("pass", "fail", "pass", "fail")) &&
    all(r$table$action == "warn"),
  "warn returns the table: one row per fit, verdicts"
)
expect(
  identical(r$table$id, c("occ:east", "occ:west", "dur:east", "dur:west")),
  "ids are <track>:<short group>"
)
m <- err_msg(deploy_check_convergence(
  plan,
  dep_g(action = "refuse"),
  diags = diag_by_file
))
expect(
  grepl("refusing", m) && grepl("occ:west", m) && grepl("dur:west", m),
  "refuse stops naming track:group"
)
out <- capture.output(
  r <- deploy_check_convergence(
    plan,
    dep_g(action = "refuse"),
    dry_run = TRUE,
    say = function(...) cat(sprintf(...), "\n", sep = ""),
    diags = diag_by_file
  )
)
expect(
  any(grepl("WOULD refuse", out)) &&
    identical(r$failed, c("occ:west", "dur:west")),
  "dry-run prints WOULD refuse and never stops"
)
expect(
  any(grepl("FAIL +occ:west .*sds_st_start_seas_1", out)),
  "the verdict line names the worst parameter"
)
r <- deploy_check_convergence(
  plan,
  dep_g(action = "refuse", allow = list("occ:west", "dur:west")),
  diags = diag_by_file
)
expect(
  identical(r$allowed, c("occ:west", "dur:west")) &&
    !length(r$failed) &&
    identical(r$table$verdict, c("pass", "allowed", "pass", "allowed")) &&
    identical(r$table$allowed, c(FALSE, TRUE, FALSE, TRUE)),
  "allow-listed failing fits pass as allowed (no refusal)"
)
reads <- 0L
rd <- function(p) {
  reads <<- reads + 1L
  list(file = basename(p))
}
dfn <- function(fit) diag_by_file[[paste0(fit$file, ".rds")]]
td <- tempfile("srcfits_")
dir.create(td)
plan_t <- plan
plan_t$src <- file.path(td, sub("\\.rds$", "", plan$file))
for (p in plan_t$src[1:3]) {
  file.create(p)
} # dur_w not produced yet
w <- warns(
  r <- deploy_check_convergence(
    plan_t,
    dep_g(action = "warn"),
    read_fn = rd,
    diag_fn = dfn
  )
)
expect(
  reads == 3L,
  "one read per existing source fit (the loo artefact is not read)"
)
expect(
  identical(r$table$verdict, c("pass", "fail", "pass", "missing")) &&
    identical(r$failed, "occ:west"),
  "a source fit not produced yet is `missing`, not a failure"
)
sr <- .deploy_read_source_fits(
  plan_t,
  episode = FALSE,
  diagnostics = FALSE,
  read_fn = stop_fn
)
expect(
  !length(sr$episode) && !length(sr$diagnostics),
  ".deploy_read_source_fits: neither -> nothing read"
)

# fit_diagnostics(make_fake_kl_fit()) -> gate.
cat("# [5] Stan-free end to end (fake kl_fit)\n")
fk <- make_fake_kl_fit(rbind(c(0, 9), c(0, 9), c(-6, 4), c(0.1, 9)))
gm <- fit_diagnostics_gate(fit_diagnostics(fk), scope = "all")
expect(
  !gm$pass && grepl("rhat_any", gm$reasons),
  "multimodal kl fit fails (R-hat > 1.05)"
)
ga <- fit_diagnostics_gate(
  fit_diagnostics(fk, chains = c(1, 2)),
  scope = "all",
  rhat_max = 1.1,
  ess_min = 10
)
expect(isTRUE(ga$pass), "its agreeing chains (1, 2) pass at lax thresholds")
fd <- make_fake_kl_fit(rbind(c(0, 9), c(0, 9)), divergent = 1L)
expect(
  grepl(
    "divergences",
    fit_diagnostics_gate(
      fit_diagnostics(fd),
      scope = "all",
      rhat_max = 1.1,
      ess_min = 10
    )$reasons
  ),
  "kl fit divergences are gated"
)

# On synthetic fits.csv / diagnostics_chains.csv: off, pass, warn, refuse
# (< 2 agreeing chains, agreeing R-hat > 1.05), dry-run, warn action never
# stops.
cat("# [6] kl_grouping_convergence()\n")
cd <- tempfile("klcmp_")
dir.create(cd)
fits_row <- function(
  id,
  agreeing,
  rhat_agree,
  essb = 500,
  esst = 500,
  rhat_all = 1.6
) {
  data.frame(
    id = id,
    basis = "ns5",
    K = 5L,
    mix = FALSE,
    config = "c",
    chains = 6L,
    n_modes = 2L,
    mode_sizes = "3+3",
    agreeing = agreeing,
    kmu_spread_agree = 0.5,
    n_params = 50L,
    divergences = 1L,
    treedepth_hits = 0L,
    min_ebfmi = 0.8,
    max_rhat_all = rhat_all,
    n_rhat105_all = 10L,
    min_ess_bulk_all = 20,
    min_ess_tail_all = 30,
    max_rhat_agree = rhat_agree,
    n_rhat105_agree = 0L,
    min_ess_bulk_agree = essb,
    min_ess_tail_agree = esst,
    lp_agree = -8900,
    lp_spread = 10,
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  rbind(
    fits_row("ok", "3,5,6", 1.004),
    fits_row("lax", "3,5,6", 1.045, essb = 52), # a lax cell
    fits_row("one", "2", NA),
    fits_row("wide", "1,2,3", 1.08)
  ),
  file.path(cd, "fits.csv"),
  row.names = FALSE
)
dc <- do.call(
  rbind,
  lapply(c("ok", "lax", "one", "wide"), function(id) {
    data.frame(
      id = id,
      agreeing = TRUE,
      chain = 1:6,
      divergent = c(1L, 0L, 0L, 0L, 0L, 0L),
      treedepth_hits = 0L,
      ebfmi = 0.9
    )
  })
)
utils::write.csv(dc, file.path(cd, "diagnostics_chains.csv"), row.names = FALSE)
kcfg <- function(action) {
  list(
    action = action,
    rhat_max = DIAG_THRESHOLDS[["rhat"]],
    ess_min = DIAG_THRESHOLDS[["ess"]],
    max_divergent = 0L
  )
}
expect(
  is.null(kl_grouping_convergence(cd, "ok", NULL, kcfg("off"))),
  "off returns NULL"
)
expect(
  is.null(kl_grouping_convergence(cd, "ok", NULL, list())),
  "no action = off"
)
w <- warns(
  k <- kl_grouping_convergence(cd, "ok", c(3L, 5L, 6L), kcfg("refuse"))
)
expect(
  identical(k$verdict, "pass") && !length(w) && k$n_agree == 3L,
  "a converged cell passes"
)
expect(
  identical(k$table$chains, c("all", "agree")) &&
    k$table$divergences[[1]] == 1L &&
    k$table$divergences[[2]] == 0L,
  "all-chain and agreeing-chain sampler rows (agreeing-only divergences)"
)
w <- warns(k <- kl_grouping_convergence(cd, "lax", c(1L, 2L), kcfg("refuse")))
expect(
  identical(k$verdict, "warn") &&
    !length(k$refuse) &&
    length(w) == 1L &&
    grepl("1.0450", w) &&
    grepl("bulk ESS 52", w),
  "agreeing R-hat 1.045 / ESS 52 warns (beyond 1.01 / 400), no refusal"
)
expect(
  identical(k$table$chains, c("all", "agree", "deploy")),
  "a deploy chain set other than the agreeing one gets its own row"
)
m <- err_msg(kl_grouping_convergence(cd, "one", NULL, kcfg("refuse")))
expect(
  grepl("refusing", m) && grepl("1 chain\\(s\\) agree", m),
  "fewer than 2 agreeing chains: refuse"
)
m <- err_msg(suppressWarnings(kl_grouping_convergence(
  cd,
  "wide",
  NULL,
  kcfg("refuse")
)))
expect(
  grepl("refusing", m) && grepl("1.0800 > 1.05", m),
  "agreeing R-hat above 1.05: refuse"
)
out <- capture.output(
  w <- warns(
    k <- kl_grouping_convergence(
      cd,
      "wide",
      NULL,
      kcfg("refuse"),
      dry_run = TRUE,
      say = function(...) cat(sprintf(...), "\n", sep = "")
    )
  )
)
expect(
  identical(k$verdict, "refuse") && any(grepl("WOULD refuse", out)),
  "dry-run prints WOULD refuse, never stops"
)
w <- warns(k <- kl_grouping_convergence(cd, "one", NULL, kcfg("warn")))
expect(
  identical(k$verdict, "warn") && any(grepl("1 chain\\(s\\) agree", w)),
  "warn action: refusal conditions only warn"
)
expect(
  errs(kl_grouping_convergence(cd, "absent", NULL, kcfg("warn"))),
  "unknown cell is an error"
)
invisible(file.remove(file.path(cd, "diagnostics_chains.csv")))
w <- warns(k <- kl_grouping_convergence(cd, "ok", NULL, kcfg("warn")))
expect(
  k$table$divergences[[1]] == 1L && is.na(k$table$divergences[[2]]),
  "without diagnostics_chains.csv the all-chain row falls back to fits.csv"
)

# On mocked lag selection / configs / fits in a tempdir: refuse stops BEFORE
# the prune (a stale fit and no MANIFEST.csv), warn deploys and writes
# CONVERGENCE.csv, off reads no fit, calls no diagnostics and removes a stale
# CONVERGENCE.csv.
cat("# [7] deploy_collect_fits() end to end (mocked)\n")
root <- tempfile("deploy_e2e_")
dir.create(root)
srcd <- file.path(root, "src")
dest <- file.path(root, "dest")
cfgd <- file.path(root, "cfg")
for (d in c(srcd, dest, cfgd)) {
  dir.create(d)
}
LAGS <- list(
  east = c(occ = 0L, dur = 1L, mon = 2L),
  west = c(occ = -6L, dur = -4L, mon = 3L)
)
deploy_select_lags <- function(dep) LAGS
.deploy_group_n_countries <- function(dep, g) 2L
.deploy_cfg_path <- function(batch_dir, group, lag) {
  p <- file.path(
    cfgd,
    sprintf("%s__%s__%d.yml", basename(batch_dir), group, lag)
  )
  file.create(p)
  p
}
make_options_from_config <- function(option_list, cli_argv, model) {
  parts <- strsplit(sub("\\.yml$", "", basename(cli_argv[[2]])), "__")[[1]]
  list(track = parts[[1]], group = parts[[2]], lag = as.integer(parts[[3]]))
}
episode_key_of <- function(opt) "lp_season"
national_zero_fill_of <- function(opt) FALSE
national_zero_threshold_of <- function(opt) 0L
mk_spec <- function(track, model, prefix, kinds) {
  list(
    track = track,
    model = model,
    prefix = prefix,
    batch_key = paste0(track, "_dir"),
    make_fn = function(opt, kind) {
      file.path(
        srcd,
        sprintf(
          "%sa_%s_x_grp-enso_grouping-%s_sl%d.rds",
          prefix,
          kind,
          opt$group,
          opt$lag
        )
      )
    },
    kinds = kinds,
    anchor_of = function(o) FALSE
  )
}
.deploy_track_specs <- function() {
  list(
    occ = mk_spec("occ", "occurrence", "outbreaks_survival_", c("fit", "loo")),
    dur = mk_spec(
      "dur",
      "duration",
      "outbreak_duration_",
      c("dur_fit", "dur_loo")
    ),
    mon = mk_spec(
      "mon",
      "monthly",
      "outbreak_monthly_",
      c("mon_fit", "mon_loo")
    )
  )
}
for (sp in .deploy_track_specs()) {
  for (g in names(LAGS)) {
    for (k in sp$kinds) {
      f <- sp$make_fn(
        list(group = paste0("enso_", g), lag = LAGS[[g]][[sp$track]]),
        k
      )
      saveRDS(list(id = paste(sp$track, g, sep = ":")), f)
    }
  }
}
dep_e <- function(...) {
  list(
    grouping = "enso_grouping",
    episodes = list(key = "lp_season"),
    fits = list(
      groups = c("east", "west"),
      dest = dest,
      occ_dir = "occ",
      dur_dir = "dur",
      mon_dir = "mon"
    ),
    convergence = list(...)
  )
}
diag_e <- function(fit) {
  if (fit$id %in% c("occ:west", "mon:east")) {
    mk_d(rhat = c(sd = 1.02))
  } else {
    mk_d()
  }
}
stale <- file.path(
  dest,
  "outbreaks_survival_old_fit_x_grp-enso_grouping-enso_east_sl9.rds"
)
saveRDS(list(), stale)
q <- function(expr) invisible(capture.output(expr))

m <- err_msg(q(deploy_collect_fits(
  dep_e(action = "refuse"),
  verbose = FALSE,
  diag_fn = diag_e
)))
expect(
  grepl("refusing", m) && grepl("occ:west", m) && grepl("mon:east", m),
  "refuse stops naming the failing fits"
)
expect(
  file.exists(stale) &&
    !file.exists(file.path(dest, "MANIFEST.csv")) &&
    length(list.files(dest)) == 1L,
  "refusal before the prune: the stale fit is still there, nothing copied"
)

q(
  res <- deploy_collect_fits(
    dep_e(action = "refuse"),
    dry_run = TRUE,
    verbose = FALSE,
    diag_fn = diag_e
  )
)
expect(
  identical(res$convergence$failed, c("occ:west", "mon:east")) &&
    file.exists(stale) &&
    !file.exists(file.path(dest, "CONVERGENCE.csv")),
  "dry-run refuse: verdicts returned, nothing written"
)

q(
  res <- deploy_collect_fits(
    dep_e(action = "refuse", allow = list("occ:west", "mon:east")),
    verbose = FALSE,
    diag_fn = diag_e
  )
)
cv <- utils::read.csv(
  file.path(dest, "CONVERGENCE.csv"),
  stringsAsFactors = FALSE
)
expect(
  !file.exists(stale) && file.exists(file.path(dest, "MANIFEST.csv")),
  "allow-listed failures deploy (stale pruned, MANIFEST written)"
)
expect(
  nrow(cv) == 6L &&
    sum(cv$verdict == "allowed") == 2L &&
    all(cv$action == "refuse") &&
    identical(sort(cv$id[cv$allowed]), c("mon:east", "occ:west")),
  "CONVERGENCE.csv: one row per fit with verdict / action / allowed"
)

w <- warns(q(
  res <- deploy_collect_fits(
    dep_e(action = "warn"),
    verbose = FALSE,
    diag_fn = diag_e
  )
))
expect(
  length(w) == 1L &&
    identical(res$convergence$failed, c("occ:west", "mon:east")) &&
    all(utils::read.csv(file.path(dest, "CONVERGENCE.csv"))$action == "warn"),
  "warn deploys, warns and rewrites CONVERGENCE.csv"
)

reads <- 0L
q(
  res <- deploy_collect_fits(
    dep_e(action = "off"),
    verbose = FALSE,
    read_fn = function(p) {
      reads <<- reads + 1L
      readRDS(p)
    },
    diag_fn = stop_fn
  )
)
expect(
  is.null(res$convergence) && reads == 0L,
  "off (lp_season): no fit read, no diagnostics"
)
expect(
  !file.exists(file.path(dest, "CONVERGENCE.csv")) &&
    file.exists(file.path(dest, "MANIFEST.csv")),
  "off removes a stale CONVERGENCE.csv and deploys as before"
)
q(
  res <- deploy_collect_fits(
    list(
      grouping = "enso_grouping",
      episodes = list(key = "lp_season"),
      fits = list(
        groups = c("east", "west"),
        dest = dest,
        occ_dir = "occ",
        dur_dir = "dur",
        mon_dir = "mon"
      )
    ),
    verbose = FALSE,
    read_fn = stop_fn,
    diag_fn = stop_fn
  )
)
expect(is.null(res$convergence), "a deploy block without `convergence` = off")
expect(
  identical(
    names(res),
    c(
      "selected",
      "manifest",
      "dest",
      "copied",
      "pruned",
      "patched",
      "dry_run",
      "convergence"
    )
  ),
  "result fields"
)

unlink(c(td, cd, root), recursive = TRUE)

# T5 (integration): the config.R default block equals the in-function
# defaults, so a deploy YAML that restates it at back-compat values behaves as
# one without the block (off).
cat(
  "\n[T5] get_default_options()$deploy$convergence = the in-function default\n"
)
dflt <- get_default_options()$deploy
expect(
  identical(
    deploy_convergence_cfg(list(fits = dflt$fits)),
    deploy_convergence_cfg(list(
      convergence = dflt$convergence,
      fits = dflt$fits
    ))
  ),
  "default block resolves exactly as a missing block"
)
expect(identical(dflt$convergence$action, "off"), "default action off")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
