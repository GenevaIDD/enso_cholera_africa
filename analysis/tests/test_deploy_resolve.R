# test_deploy_resolve.R — HARD: the deploy YAML is the single source of run
#   identity.
#
# Usage:
#   Rscript analysis/tests/test_deploy_resolve.R
#
# Details:
#   Covers the deploy-tool half of "no hardcoded run numbers"
#   (utils/deploy_core.R + utils/deploy_prod_run.R): the defaults name no run
#   [1], shipped configs resolve [2], the CLI KEY=value contract [3] and
#   gen_batches planning [4], [5]. No fit, no data; [3] and [5] shell out to
#   Rscript (a few seconds each).

suppressPackageStartupMessages(library(yaml))

source("analysis/utils/deploy_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
# The group-id prefix of a fit registry (its sidecar's label_prefix; enso_ when absent).
grouping_prefix <- function(reg) {
  f <- file.path("analysis/groupings", paste0(reg, ".yml"))
  p <- if (file.exists(f)) yaml::read_yaml(f)$label_prefix else NULL
  if (is.null(p)) "enso_" else p
}

cat("# deploy YAML is the single source of run identity\n")

DEPLOY_DIR <- "analysis/configs/deploy"
CONFIG_PATHS <- sort(list.files(
  DEPLOY_DIR,
  pattern = "\\.yml$",
  full.names = TRUE
))
PROD_PATHS <- CONFIG_PATHS[grepl("^prod_run", basename(CONFIG_PATHS))]
REAL_CFG <- PROD_PATHS[length(PROD_PATHS)] # newest production run by name
expect(length(CONFIG_PATHS) >= 1L, "at least one deploy config exists")

run_cli <- function(args) run_rscript("analysis/utils/deploy_prod_run.R", args) # helpers.R
parse_keys <- function(lines) {
  kv <- grep("^[A-Z_]+=", lines, value = TRUE)
  stats::setNames(as.list(sub("^[A-Z_]+=", "", kv)), sub("=.*$", "", kv))
}

# --- [1] defaults name no run -------------------------------------------------
# Every run-specific leaf is the "" sentinel and deploy_require_keys() refuses
# it, naming each key.
cat("\n[1] get_default_options()$deploy names no run\n")
defaults <- get_default_options()$deploy
leaf_vals <- vapply(
  .DEPLOY_REQUIRED_KEYS,
  function(k) paste(.deploy_leaf(defaults, k), collapse = ""),
  character(1)
)
expect(
  all(leaf_vals == ""),
  sprintf(
    "every run-specific default is \"\" (non-empty: %s)",
    if (any(leaf_vals != "")) {
      paste(names(leaf_vals)[leaf_vals != ""], collapse = ", ")
    } else {
      "none"
    }
  )
)
err <- tryCatch(
  {
    deploy_require_keys(defaults, "defaults")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(err) &&
    all(vapply(
      .DEPLOY_REQUIRED_KEYS,
      function(k) grepl(paste0("deploy.", k), err, fixed = TRUE),
      logical(1)
    )),
  "deploy_require_keys() on the defaults fails, naming every run-specific key"
)

# --- [2] shipped configs resolve to their own values --------------------------
# Every shipped deploy/*.yml passes deploy_require_keys(), and
# deploy_resolved_keys() returns exactly the YAML's values (occ/dur/mon only:
# the size track is retired downstream, so no SIZ_* key and no fits.siz_dir in
# any shipped YAML).
cat("\n[2] shipped deploy configs pass and resolve to their YAML values\n")
for (cfg in CONFIG_PATHS) {
  raw <- yaml::read_yaml(cfg)$deploy
  dep <- load_deploy_config(cfg)
  ok <- tryCatch(
    {
      deploy_require_keys(dep, cfg)
      TRUE
    },
    error = function(e) FALSE
  )
  expect(ok, sprintf("%s sets every run-specific key", basename(cfg)))
  k <- deploy_resolved_keys(dep, cfg)
  expect(
    identical(k$OCC_BATCH, raw$fits$occ_dir) &&
      identical(k$DUR_BATCH, raw$fits$dur_dir) &&
      identical(k$MON_BATCH, raw$fits$mon_dir) &&
      identical(k$OCC_DIR, file.path("analysis/configs", raw$fits$occ_dir)) &&
      identical(k$DUR_DIR, file.path("analysis/configs", raw$fits$dur_dir)) &&
      identical(k$MON_DIR, file.path("analysis/configs", raw$fits$mon_dir)) &&
      identical(k$FORECAST_CONFIG_DIR, raw$outputs$forecast_config_outdir) &&
      identical(k$FORECAST_MANIFEST, raw$outputs$forecast_manifest) &&
      identical(k$HINDCAST_CONFIG_DIR, raw$outputs$hindcast_config_outdir) &&
      identical(k$HINDCAST_MANIFEST, raw$outputs$hindcast_manifest) &&
      identical(k$SMOOTHING_INPUT, raw$smoothing_input) &&
      identical(k$HISTORICAL_TEMPLATE, raw$templates$historical) &&
      identical(k$FITS_DEST, raw$fits$dest) &&
      # "" while a kl sweep has yet to derive the registry (before gen_grouping)
      identical(
        k$GROUPS,
        if (
          !file.exists(registry_path(raw$grouping)) &&
            nzchar(raw$kl_sweep$batch %||% "")
        ) {
          ""
        } else {
          paste0(grouping_prefix(raw$grouping), raw$fits$groups, collapse = " ")
        }
      ) &&
      identical(k$PARTITION_BY, raw$partition_by) &&
      identical(k$DEPLOY_TAG, tools::file_path_sans_ext(basename(cfg))) &&
      identical(
        k$OBS_PANEL,
        paste0(
          "generated_data/outbreaks_survival_data_base-1990_2019",
          if (identical(raw$episodes$key, "gid_season")) "_ep-gs" else "",
          if (isTRUE(raw$episodes$national_zero_fill)) {
            paste0(
              "_nz",
              if (isTRUE((raw$episodes$national_zero_threshold %||% 0) > 0)) {
                raw$episodes$national_zero_threshold
              } else {
                ""
              }
            )
          } else {
            ""
          },
          # national-zero rule and panel tokens: -<rule>, _ae, _a1c
          if (nzchar(raw$episodes$national_zero_rule %||% "")) {
            paste0("-", raw$episodes$national_zero_rule)
          } else {
            ""
          },
          if (identical(raw$episodes$at_risk_start, "after_episode")) {
            "_ae"
          } else {
            ""
          },
          if (identical(raw$episodes$adm1_overlap, "censor_orphan")) {
            "_a1c"
          } else {
            ""
          },
          ".rds"
        )
      ),
    sprintf("%s: deploy_resolved_keys() mirrors the YAML", basename(cfg))
  )
  # T-int 7(c): the counterfactual Lenssen cache = the deploy overlay's
  # cache_path, else the skeleton's shared legacy cache.
  cf_ov <- (raw$sim$counterfactual %||% list())[[
    "simulator.enso.counterfactual.cache_path"
  ]] %||%
    (raw$sim$shared %||% list())[["simulator.enso.counterfactual.cache_path"]]
  expect(
    identical(
      k$COUNTERFACTUAL_CACHE,
      cf_ov %||% "generated_data/lenssen_nino34_analogs.rds"
    ),
    sprintf(
      "%s: COUNTERFACTUAL_CACHE = %s",
      basename(cfg),
      k$COUNTERFACTUAL_CACHE
    )
  )
  expect(
    !any(grepl("^SIZ_", names(k))) && is.null(raw$fits$siz_dir),
    sprintf("%s: no size batch (no fits.siz_dir, no SIZ_* key)", basename(cfg))
  )
}

# --- [3] the CLI contract -----------------------------------------------------
# `deploy_prod_run.R --stage resolve` prints those values as KEY=value lines
# (the contract bash drivers parse), and a YAML missing a run key fails.
cat("\n[3] deploy_prod_run.R --stage resolve prints the KEY=value contract\n")
r <- run_cli(c("-C", REAL_CFG, "--stage", "resolve"))
got <- parse_keys(r$out)
want <- deploy_resolved_keys(load_deploy_config(REAL_CFG), REAL_CFG)
expect(
  r$status == 0L,
  sprintf("--stage resolve on %s exits 0", basename(REAL_CFG))
)
expect(
  identical(names(got), names(want)) &&
    identical(unname(unlist(got)), unname(unlist(want))),
  "stdout KEY=value lines equal deploy_resolved_keys(), in order"
)
expect(
  !any(grepl("^\\[deploy\\]", r$out)),
  "--stage resolve prints no [deploy] banner"
)

tmp_cfg <- tempfile("prod_runX_", fileext = ".yml")
on.exit(unlink(tmp_cfg), add = TRUE)
bad <- yaml::read_yaml(REAL_CFG)
bad$deploy$fits$occ_dir <- NULL
yaml::write_yaml(bad, tmp_cfg)
r_bad <- run_cli(c("-C", tmp_cfg, "--stage", "resolve"))
expect(
  r_bad$status != 0L &&
    any(grepl("deploy.fits.occ_dir", r_bad$out, fixed = TRUE)),
  "a deploy YAML without fits.occ_dir fails --stage resolve, naming the key"
)

# --- [4] gen_batches planning (tempdir) ----------------------------------------
# "generate" for a registered batch with no dir, "exists" once its
# manifest.yml is written, "unregistered" for an unknown name. Generation is
# redirected to a tempdir (never analysis/configs/).
cat("\n[4] gen_batches plan: generate / exists / unregistered\n")
env <- deploy_batch_env()
expect(
  exists("generate_batch", envir = env, inherits = FALSE) &&
    !exists("generate_batch", envir = globalenv(), inherits = FALSE),
  "make_batch_configs.R loads into its own env (nothing leaks into globalenv)"
)
tmp_dir <- tempfile("cfgdir_")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
real_config_dir <- env$get_config_dir
assign("get_config_dir", function() tmp_dir, envir = env)

dep <- load_deploy_config(REAL_CFG)
plan <- deploy_plan_batches(dep, env)
fit_rows <- plan$track %in% c("occ", "dur", "mon")
has_kl <- nzchar(dep$kl_sweep$batch %||% "")
expect(
  sum(fit_rows) ==
    length(unique(unlist(dep$fits[c("occ_dir", "dur_dir", "mon_dir")]))) &&
    all(plan$track[!fit_rows] == "kl") &&
    sum(!fit_rows) == as.integer(has_kl),
  "one plan row per distinct occ/dur/mon batch dir, plus the kl sweep row when the deploy has one"
)
expect(
  all(plan$status == "generate"),
  "every registered batch is 'generate' in an empty config dir"
)

invisible(capture.output(env$generate_batch(plan$batch[1])))
plan2 <- deploy_plan_batches(dep, env)
expect(
  file.exists(file.path(tmp_dir, plan$batch[1], "manifest.yml")) &&
    identical(plan2$status[1], "exists") &&
    all(plan2$status[-1] == "generate"),
  sprintf("after generate_batch('%s') its row turns 'exists'", plan$batch[1])
)

dep_unreg <- dep
dep_unreg$fits$occ_dir <- "no_such_batch_grouped_enso_singlelag"
plan3 <- deploy_plan_batches(dep_unreg, env)
expect(
  identical(plan3$status[plan3$track == "occ"], "unregistered"),
  "a batch name missing from BATCH_REGISTRY is 'unregistered'"
)

# --- [5] gen_batches --dry-run on the real tree -----------------------------------
# On a real deploy: skips every existing batch dir and writes nothing.
cat("\n[5] --stage gen_batches --dry-run writes nothing\n")
prov <- make_deploy_run_file(
  list(outdir = "generated_data/"),
  "deploy_manifest",
  config_path = REAL_CFG
)
mtime_before <- if (file.exists(prov)) file.mtime(prov) else NA
assign("get_config_dir", real_config_dir, envir = env)
n_exist <- sum(deploy_plan_batches(dep, env)$status == "exists")
r_gb <- run_cli(c("-C", REAL_CFG, "--stage", "gen_batches", "--dry-run"))
expect(
  r_gb$status == 0L,
  sprintf("gen_batches --dry-run on %s exits 0", basename(REAL_CFG))
)
# n_exist is 0 only in a code-release checkout (analysis/RELEASE: no generated batches).
expect(
  (n_exist > 0L || file.exists("analysis/RELEASE")) &&
    sum(grepl(
      "exists (manifest.yml present), skipped",
      r_gb$out,
      fixed = TRUE
    )) ==
      n_exist,
  sprintf("the %d existing batch dir(s) are reported skipped", n_exist)
)
mtime_after <- if (file.exists(prov)) file.mtime(prov) else NA
expect(
  identical(mtime_before, mtime_after),
  "the deploy provenance manifest is untouched"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
