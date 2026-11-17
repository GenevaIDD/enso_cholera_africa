# test_deploy_sim_templates.R — HARD: `deploy_prod_run.R --stage
#   gen_templates` generates a rollout's three spatPomp templates.
#
# Usage:
#   Rscript analysis/tests/test_deploy_sim_templates.R
#
# Details:
#   Templates are generated from analysis/configs/sim_templates/<arm>.yml +
#   deploy sim.shared + sim.<arm> + leaves derived from the deploy block
#   (deploy_core.R), instead of being copied and hand-patched per run. No fit,
#   no data; [6] shells out to Rscript (a few seconds per call).

suppressPackageStartupMessages(library(yaml))

source("analysis/utils/deploy_core.R")
source("analysis/utils/sim_io.R") # assert_arm_launch() (deploy_check_sim_template(), call-time)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# deploy sim templates are generated from skeleton + sim: overlay\n")

DEPLOY_DIR <- "analysis/configs/deploy"
CONFIG_PATHS <- sort(list.files(
  DEPLOY_DIR,
  pattern = "\\.yml$",
  full.names = TRUE
))
TODAY <- as.Date("2026-10-06")

errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
templates_exist <- function(dep) all(file.exists(unlist(dep$templates)))

DEPS <- stats::setNames(
  lapply(CONFIG_PATHS, load_deploy_config),
  basename(CONFIG_PATHS)
)
WITH_TPL <- names(DEPS)[vapply(DEPS, templates_exist, logical(1))]
# The synthetic run (synth_run*) is not a reference: its registry partition
# would change section [8]'s premise. Templates are deploy-generated, so a
# code-release checkout (analysis/RELEASE) has nothing to compare.
WITH_TPL <- WITH_TPL[!grepl("^synth_run", WITH_TPL)]
if (!length(WITH_TPL) && !file.exists("analysis/RELEASE")) {
  cat("FAIL: no deploy YAML has its templates on disk\n")
  quit(status = 1)
}
if (!length(WITH_TPL)) {
  cat(
    "SKIP: no deploy YAML has its templates on disk (run --stage gen_templates first)\n"
  )
  quit(status = 0)
}
REF_NAME <- WITH_TPL[length(WITH_TPL)]
REF_CFG <- file.path(DEPLOY_DIR, REF_NAME)
REF <- DEPS[[REF_NAME]]

IGNORE <- .DEPLOY_SIM_IGNORED_LEAVES

# --- [1] skeletons ---------------------------------------------------------------
# The skeletons exist, set no derived leaf and resolve as configs.
cat("\n[1] run-agnostic skeletons\n")
derived_names <- unique(unlist(lapply(.DEPLOY_SIM_ARMS, function(arm) {
  names(.deploy_sim_derived_leaves(REF, arm, TODAY))
})))
for (arm in .DEPLOY_SIM_ARMS) {
  sk <- deploy_sim_skeleton_path(arm)
  expect(file.exists(sk), sprintf("%s exists", sk))
  s <- yaml::read_yaml(sk)
  leaked <- intersect(
    config_leaf_paths(s),
    c(derived_names, .DEPLOY_SIM_REQUIRED_LEAVES[[arm]], "episodes.key")
  )
  expect(
    !length(leaked),
    sprintf(
      "%s sets no derived / required leaf (found: %s)",
      arm,
      if (length(leaked)) paste(leaked, collapse = ", ") else "none"
    )
  )
  e <- errmsg(invisible(capture.output(
    make_options_from_config(list(), cli_argv = c("--config", sk))
  )))
  expect(
    is.null(e),
    sprintf(
      "%s resolves through make_options_from_config() (%s)",
      arm,
      e %||% "ok"
    )
  )
}

# --- [2] backfill acceptance -----------------------------------------------------
# BACKFILL ACCEPTANCE: every deploy YAML whose templates exist reproduces them
# leaf for leaf (all but simulator.run_id): the skeleton-plus-overlay model
# expresses every real run, including country-level smoothing.
cat(
  "\n[2] every deploy with templates on disk reproduces them (all leaves but run_id)\n"
)
for (nm in names(DEPS)) {
  if (!nm %in% WITH_TPL) {
    cat(sprintf("  SKIP  %s: templates not written yet\n", nm))
    next
  }
  for (arm in .DEPLOY_SIM_ARMS) {
    built <- deploy_build_sim_template(DEPS[[nm]], arm, TODAY)
    disk <- yaml::read_yaml(DEPS[[nm]]$templates[[arm]])
    d <- deploy_leaf_diff(disk, built, ignore = IGNORE)
    bad <- d[d$status != "match", , drop = FALSE]
    expect(
      !nrow(bad),
      sprintf(
        "%s %-14s reproduced (%d leaves; differing: %s)",
        nm,
        arm,
        nrow(d),
        if (nrow(bad)) paste(bad$leaf, collapse = ", ") else "none"
      )
    )
  }
}

# --- [3] refusals ----------------------------------------------------------------
# The overlay cannot set a derived leaf, must set the required ones, and a
# deploy whose smoothing_input / templates.* disagree with what it derives is
# refused; no arm derives a location_frame leaf and no skeleton or built
# template carries a retired size leaf (simulator.size.*,
# location_frame.mean_rate_table_path).
cat(
  "\n[3] derived leaves are protected, required leaves enforced, deploy keys cross-checked\n"
)
d1 <- REF
d1$sim$historical[["simulator.run_tag"]] <- "x"
e <- errmsg(deploy_build_sim_template(d1, "historical", TODAY))
expect(
  !is.null(e) && grepl("simulator.run_tag", e, fixed = TRUE),
  "an overlay setting simulator.run_tag is refused, naming the leaf"
)
d1 <- REF
d1$sim$shared[["simulator.grouping"]] <- list(source = "x")
e <- errmsg(deploy_build_sim_template(d1, "forecast", TODAY))
expect(
  !is.null(e) && grepl("simulator.grouping", e, fixed = TRUE),
  "an overlay replacing a subtree that holds derived leaves is refused"
)
d1 <- REF
d1$sim$shared[["episodes.key"]] <- "gid_season"
expect(
  !is.null(errmsg(deploy_build_sim_template(d1, "forecast", TODAY))),
  "an overlay setting episodes.key is refused (deploy.episodes.key is the authority)"
)
d1 <- REF
d1$sim$shared <- list()
d1$sim$forecast <- list()
e <- errmsg(deploy_build_sim_template(d1, "forecast", TODAY))
expect(
  !is.null(e) && grepl("simulator.filter.measurement", e, fixed = TRUE),
  "a run whose overlay omits simulator.filter.measurement is refused"
)
d1 <- REF
d1$sim$historical[["simulator.filter.partition_by"]] <- NULL
e <- errmsg(deploy_build_sim_template(d1, "historical", TODAY))
expect(
  !is.null(e) && grepl("simulator.filter.partition_by", e, fixed = TRUE),
  "a historical overlay without filter.partition_by is refused"
)
d1 <- REF
d1$smoothing_input <- "generated_data/sim_spatpomp/other/cache"
e <- errmsg(deploy_build_sim_template(d1, "counterfactual", TODAY))
expect(
  !is.null(e) && grepl("smoothing_input", e, fixed = TRUE),
  "deploy.smoothing_input != the historical run's cache is refused"
)
d1 <- REF
d1$templates$forecast <- "analysis/configs/elsewhere.yml"
e <- errmsg(deploy_sim_template_path(d1, "forecast"))
expect(
  !is.null(e) && grepl("templates.forecast", e, fixed = TRUE),
  "templates.<arm> must name analysis/configs/<run_tag>/<run_tag>.yml"
)
RETIRED_SIZE <- function(leaves) {
  leaves[
    startsWith(leaves, "simulator.size.") |
      leaves == "simulator.size" |
      grepl("mean_rate_table_path", leaves, fixed = TRUE)
  ]
}
for (arm in .DEPLOY_SIM_ARMS) {
  dl <- names(.deploy_sim_derived_leaves(REF, arm, TODAY))
  expect(
    !any(startsWith(dl, "simulator.location_frame.")),
    sprintf(
      "%s derives no location_frame leaf (derived: %s)",
      arm,
      paste(dl, collapse = ", ")
    )
  )
  sk <- RETIRED_SIZE(config_leaf_paths(yaml::read_yaml(deploy_sim_skeleton_path(
    arm
  ))))
  expect(
    !length(sk),
    sprintf(
      "%s skeleton carries no retired size leaf (found: %s)",
      arm,
      if (length(sk)) paste(sk, collapse = ", ") else "none"
    )
  )
  for (nm in names(DEPS)) {
    bt <- RETIRED_SIZE(config_leaf_paths(deploy_build_sim_template(
      DEPS[[nm]],
      arm,
      TODAY
    )))
    expect(
      !length(bt),
      sprintf(
        "%s %-14s built template carries no retired size leaf (found: %s)",
        nm,
        arm,
        if (length(bt)) paste(bt, collapse = ", ") else "none"
      )
    )
  }
}
t_ok <- errmsg(deploy_check_sim_template(REF, "forecast"))
expect(
  is.null(t_ok),
  sprintf("%s forecast template passes deploy_check_sim_template()", REF_NAME)
)
d1 <- REF
d1$smoothing_input <- "generated_data/sim_spatpomp/other/cache"
e <- errmsg(deploy_check_sim_template(d1, "forecast"))
expect(
  !is.null(e) && grepl("smoothing_input", e, fixed = TRUE),
  "deploy_check_sim_template() catches a forward template on another smoothing cache"
)
# Per-arm launch: handoff_month / forecast_window.start must equal
# horizon.start (sim_io.R::assert_arm_launch()); a differing analog anchor only reports.
for (arm in c("forecast", "counterfactual")) {
  tpl <- yaml::read_yaml(REF$templates[[arm]])
  tp <- tempfile(paste0("launch_", arm, "_"), fileext = ".yml")
  bad <- tpl
  bad$simulator$filter$handoff_month <- "2031-01-01"
  yaml::write_yaml(bad, tp)
  e <- errmsg(deploy_check_sim_template(REF, arm, tp))
  expect(
    !is.null(e) && grepl("handoff_month", e, fixed = TRUE),
    sprintf("%s template with handoff_month != horizon.start refused", arm)
  )
  bad <- tpl
  bad$simulator$filter$forecast_window$start <- "2031-01-01"
  yaml::write_yaml(bad, tp)
  e <- errmsg(deploy_check_sim_template(REF, arm, tp))
  expect(
    !is.null(e) && grepl("forecast_window.start", e, fixed = TRUE),
    sprintf(
      "%s template with forecast_window.start != horizon.start refused",
      arm
    )
  )
}
tpl <- yaml::read_yaml(REF$templates$counterfactual)
tpl$simulator$enso$counterfactual$lenssen_horizon_start <- "2020-01-01"
tp <- tempfile("launch_anchor_", fileext = ".yml")
yaml::write_yaml(tpl, tp)
msg <- character(0)
e <- withCallingHandlers(
  errmsg(deploy_check_sim_template(REF, "counterfactual", tp)),
  message = function(m) {
    msg <<- c(msg, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
expect(
  is.null(e) && any(grepl("analog anchor", msg, fixed = TRUE)),
  "counterfactual analog anchor != launch is reported, not refused"
)

# --- [4] run_id ------------------------------------------------------------------
# run_id is date + occ/dur/mon fit lineage and never makes templates differ.
cat("\n[4] run_id = date + fit lineage; provenance only\n")
d1 <- REF
d1$fits$occ_dir <- "prod_run9_occ_grouped_enso_singlelag"
d1$fits$dur_dir <- "prod_run8_dur_grouped_enso_singlelag"
d1$fits$mon_dir <- "prod_run8_mon_grouped_enso_singlelag"
rid <- deploy_sim_run_id(d1, "historical", TODAY)
expect(
  identical(rid, "20261006_historical_prodrun9occ_prodrun8durmon"),
  sprintf("run_id groups tracks by the run that fit them (%s)", rid)
)
d1$fits$siz_dir <- "prod_run5_size_grouped_enso_singlelag"
expect(
  identical(deploy_sim_run_id(d1, "historical", TODAY), rid),
  "a stray fits.siz_dir does not enter the lineage"
)
a <- deploy_build_sim_template(REF, "forecast", TODAY)
b <- deploy_build_sim_template(REF, "forecast", TODAY + 30)
expect(
  !identical(a$simulator$run_id, b$simulator$run_id) &&
    deploy_sim_templates_equal(a, b),
  "templates differing only in run_id compare equal"
)
expect(
  identical(names(a$simulator)[1:3], c("backend", "run_tag", "run_id")),
  "run_tag / run_id are written at the top of the simulator block"
)

# --- [5] write / read-back round trip ---------------------------------------------
# A written template (header + as.yaml) reads back equal, resolves through
# make_options_from_config() (the retired-size refusal does not fire), and a
# changed leaf shows up in the diff.
cat("\n[5] a written template reads back equal; a changed leaf is reported\n")
tmp <- tempfile("tpl_", fileext = ".yml")
on.exit(unlink(tmp), add = TRUE)
for (arm in .DEPLOY_SIM_ARMS) {
  cfg <- deploy_build_sim_template(REF, arm, TODAY)
  write_yaml_with_header(cfg, tmp, deploy_sim_template_header(REF_CFG, arm))
  back <- yaml::read_yaml(tmp)
  expect(
    nrow(deploy_leaf_diff(cfg, back)) > 0L &&
      all(deploy_leaf_diff(cfg, back)$status == "match"),
    sprintf(
      "%s: written + re-read template equals the built one on every leaf",
      arm
    )
  )
  expect(
    grepl(
      "^# GENERATED by analysis/utils/deploy_prod_run.R --stage gen_templates",
      readLines(tmp, 1)
    ),
    sprintf("%s: file starts with the GENERATED banner", arm)
  )
  e <- errmsg(invisible(capture.output(
    make_options_from_config(list(), cli_argv = c("--config", tmp))
  )))
  expect(
    is.null(e),
    sprintf(
      "%s: the generated template resolves through make_options_from_config() (%s)",
      arm,
      e %||% "ok"
    )
  )
}
mut <- back
mut$simulator$filter$Np <- 1L
d <- deploy_leaf_diff(back, mut, ignore = .DEPLOY_SIM_IGNORED_LEAVES)
expect(
  identical(d$leaf[d$status != "match"], "simulator.filter.Np"),
  "changing one leaf yields exactly that leaf in the diff"
)

# --- [6] the CLI -----------------------------------------------------------------
# --dry-run reports the real templates unchanged and writes nothing; a
# differing overlay is refused before anything is written; an absent template
# is planned, not written, under --dry-run.
cat("\n[6] deploy_prod_run.R --stage gen_templates\n")
run_cli <- function(args) run_rscript("analysis/utils/deploy_prod_run.R", args) # helpers.R
tpl_paths <- unlist(REF$templates)
prov <- make_deploy_run_file(
  list(outdir = "generated_data/"),
  "deploy_manifest",
  config_path = REF_CFG
)
mt_before <- c(
  file.mtime(tpl_paths),
  prov = if (file.exists(prov)) file.mtime(prov) else NA
)

r <- run_cli(c("-C", REF_CFG, "--stage", "gen_templates", "--dry-run"))
expect(r$status == 0L, sprintf("--dry-run on %s exits 0", REF_NAME))
expect(
  sum(grepl("unchanged (left as is)", r$out, fixed = TRUE)) == 3L,
  "all three real templates are reported unchanged"
)

raw <- yaml::read_yaml(REF_CFG)
tmp_cfg <- tempfile("prod_runX_", fileext = ".yml")
on.exit(unlink(tmp_cfg), add = TRUE)
mut_cfg <- raw
mut_cfg$deploy$sim$shared[["simulator.filter.delta_fn"]] <- 0.5
yaml::write_yaml(mut_cfg, tmp_cfg)
r <- run_cli(c("-C", tmp_cfg, "--stage", "gen_templates"))
expect(
  r$status != 0L &&
    any(grepl("DIFFERS", r$out, fixed = TRUE)) &&
    any(grepl("simulator.filter.delta_fn", r$out, fixed = TRUE)),
  "a differing overlay is refused (non-dry), printing the differing leaf"
)

mt_after <- c(
  file.mtime(tpl_paths),
  prov = if (file.exists(prov)) file.mtime(prov) else NA
)
expect(
  identical(mt_before, mt_after),
  "no real template and no provenance manifest was touched"
)

abs_cfg <- raw
tag <- function(arm) paste0("zz_test_b5c_", arm, "_", Sys.getpid())
for (arm in .DEPLOY_SIM_ARMS) {
  abs_cfg$deploy$run_tags[[arm]] <- tag(arm)
  abs_cfg$deploy$templates[[arm]] <- file.path(
    "analysis/configs",
    tag(arm),
    paste0(tag(arm), ".yml")
  )
}
abs_cfg$deploy$smoothing_input <- file.path(
  .DEPLOY_SIM_ROOT,
  tag("historical"),
  "cache"
)
yaml::write_yaml(abs_cfg, tmp_cfg)
r <- run_cli(c("-C", tmp_cfg, "--stage", "gen_templates", "--dry-run"))
expect(
  r$status == 0L &&
    sum(grepl("WOULD write (absent)", r$out, fixed = TRUE)) == 3L,
  "absent templates are planned under --dry-run"
)
expect(
  !any(dir.exists(file.path(
    "analysis/configs",
    vapply(.DEPLOY_SIM_ARMS, tag, "")
  ))),
  "--dry-run creates no template dir"
)

# --- [7] shared helpers ----------------------------------------------------------
# The shared helpers it is built on (config.R leaf helpers, deploy_leaf_diff(),
# generated_yaml_header()).
cat("\n[7] shared helpers\n")
x <- list(a = list(b = 1L, c = list()), d = c("p", "q"), e = list("u", "v"))
expect(
  identical(config_leaf_paths(x), c("a.b", "a.c", "d", "e")),
  "config_leaf_paths(): {} / atomic vector / sequence are single leaves"
)
y <- config_set_leaf(x, "a.f.g", 2)
expect(
  identical(config_get_leaf(y, "a.f.g"), 2),
  "config_set_leaf() creates intermediate levels"
)
y <- config_set_leaf(y, "a.b", NULL)
expect(
  is.null(y$a$b) && !"b" %in% names(y$a),
  "config_set_leaf(..., NULL) deletes the leaf"
)
expect(
  is.null(config_get_leaf(x, "a.zz.q")),
  "config_get_leaf() is NULL on a missing path"
)
d <- deploy_leaf_diff(
  list(a = 1, b = list(c = 2), x = 9),
  list(a = 1, b = list(c = 3), y = 9, x = 9L),
  ignore = "x"
)
expect(
  identical(d$status, c("match", "differs", "only_b")) &&
    identical(d$leaf, c("a", "b.c", "y")),
  "deploy_leaf_diff(): match / differs / only_b, ignore drops a leaf, int == double"
)
h <- generated_yaml_header("tool.R", "# body")
expect(
  grepl("^# GENERATED by tool.R on .* DO NOT EDIT BY HAND\\.$", h[1]) &&
    identical(h[-1], c("# body", "")),
  "generated_yaml_header(): banner, body, trailing blank line"
)

# --- [8] registry-partitioned deploys -------------------------------------------
cat(
  "\n[8] partition_by: registry derives the partition + contagion registries\n"
)
dn <- REF
dn$partition_by <- "country"
dn$partition_source <- ""
dn$sim$historical[["simulator.filter.partition_by"]] <- "transmission_unit"
for (arm in .DEPLOY_SIM_ARMS) {
  b <- deploy_build_sim_template(dn, arm, TODAY)
  expect(
    is.null(config_get_leaf(b, "simulator.filter.partition_source")) &&
      is.null(config_get_leaf(b, "simulator.nbr_contagion.partition_source")),
    sprintf("%s: no registry leaves unless partition_by = registry", arm)
  )
}
dr <- REF
dr$partition_by <- "registry"
dr$partition_source <- "smoke_grouping"
dr$sim$historical[["simulator.filter.partition_by"]] <- "registry"
for (arm in .DEPLOY_SIM_ARMS) {
  b <- deploy_build_sim_template(dr, arm, TODAY)
  expect(
    identical(
      config_get_leaf(b, "simulator.filter.partition_source"),
      "smoke_grouping"
    ) &&
      identical(
        config_get_leaf(b, "simulator.nbr_contagion.partition_source"),
        "smoke_grouping"
      ),
    sprintf(
      "%s: partition_source + nbr_contagion.partition_source derived (default = partition)",
      arm
    )
  )
}
dc <- dr
dc$contagion$partition_source <- "enso_grouping"
expect(
  identical(
    config_get_leaf(
      deploy_build_sim_template(dc, "forecast", TODAY),
      "simulator.nbr_contagion.partition_source"
    ),
    "enso_grouping"
  ),
  "deploy.contagion.partition_source overrides the contagion registry"
)
de <- dr
de$partition_source <- ""
e <- errmsg(deploy_build_sim_template(de, "forecast", TODAY))
expect(
  !is.null(e) && grepl("partition_source", e),
  "registry deploy without partition_source refused"
)
do <- dr
do$sim$shared[["simulator.filter.partition_source"]] <- "x"
e <- errmsg(deploy_build_sim_template(do, "forecast", TODAY))
expect(
  !is.null(e) && grepl("partition_source", e),
  "overlay setting the derived partition_source refused"
)
# deploy_check_sim_template on a written registry template
tdir <- tempfile("regtpl_")
dir.create(tdir)
tp <- file.path(tdir, "hist.yml")
hb <- deploy_build_sim_template(dr, "historical", TODAY)
yaml::write_yaml(hb, tp)
expect(
  is.null(errmsg(deploy_check_sim_template(dr, "historical", tp))),
  "matching registry template accepted"
)
hb2 <- config_set_leaf(hb, "simulator.filter.partition_source", "hmm_clusters")
yaml::write_yaml(hb2, tp)
e <- errmsg(deploy_check_sim_template(dr, "historical", tp))
expect(
  !is.null(e) && grepl("partition_source", e),
  "template with another partition registry refused"
)
hb3 <- config_set_leaf(hb, "simulator.filter.partition_by", "enso_region")
yaml::write_yaml(hb3, tp)
e <- errmsg(deploy_check_sim_template(dr, "historical", tp))
expect(
  !is.null(e) && grepl("historical smoothing", e),
  "registry deploy with enso_region smoothing refused"
)
k <- deploy_resolved_keys(dr, REF_CFG)
expect(
  identical(k$PARTITION_BY, "registry") &&
    identical(k$PARTITION_SOURCE, "smoke_grouping") &&
    identical(k$GROUPING, REF$grouping),
  "resolve emits PARTITION_SOURCE and GROUPING"
)
unlink(tdir, recursive = TRUE)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
