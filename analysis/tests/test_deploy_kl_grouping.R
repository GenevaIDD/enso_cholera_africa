# test_deploy_kl_grouping.R — HARD: the fit registry derived from a pinned cell
#   of a brms_kl class sweep.
#
# Usage:
#   Rscript analysis/tests/test_deploy_kl_grouping.R
#
# Details:
#   Covers utils/kl_grouping.R and `deploy_prod_run.R --stage gen_grouping`
#   (relay phase `classes` of bash/tools/run_pipeline.sh) on Stan-free fake
#   kl_fit objects (tests/helpers_synth.R::make_fake_kl_fit()), everything in
#   a tempdir. No sampling, no data.

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/deploy_core.R")
source("analysis/utils/results_core.R") # group_display_meta() ([3b])
for (m in c("dlnm", "dlnm_kernel", "kl_grouping")) {
  source(sprintf("analysis/utils/%s.R", m))
}
source("analysis/tests/helpers_synth.R")

gdir <- tempfile("groupings_")
dir.create(gdir)
on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
labs <- c("enso_east", "enso_west", "enso_nga")
ctry <- c("TZA", "KEN", "NER", "MLI", "NGA", "SEN")
cls <- c(1L, 1L, 2L, 2L, 3L, 2L)
# chains 1-3 put every country in its class; chain 4 swaps east and west
cen <- rbind(c(0, 5, 9), c(0, 5, 9), c(0, 5, 9), c(5, 0, 9))
fit <- make_fake_kl_fit(cen, countries = ctry, class_of = cls, labels = labs)
# meta_from registry: enso_east / enso_west / enso_south with a sidecar
utils::write.csv(
  data.frame(
    iso3 = c("TZA", "KEN", "NER", "MLI", "SEN", "NGA", "GHA", "ZWE"),
    group = c(
      "enso_east",
      "enso_east",
      "enso_west",
      "enso_west",
      "enso_west",
      "enso_west",
      "enso_west",
      "enso_south"
    )
  ),
  file.path(gdir, "ref.csv"),
  row.names = FALSE,
  quote = FALSE
)
writeLines(
  c(
    'label_prefix: "enso_"',
    "order: [enso_east, enso_west, enso_south]",
    "groups:",
    '  enso_east:  {label: "East Africa", short: "East", short_fig: "East", colour: "#1b9e77"}',
    '  enso_west:  {label: "West Africa", short: "West", short_fig: "West", colour: "#d95f02"}',
    '  enso_south: {label: "Southern Africa", short: "Southern", short_fig: "South", colour: "#66a61e"}'
  ),
  file.path(gdir, "ref.yml")
)

# Argmax partition, chain subset, fill_from.
cat("# [1] kl_grouping_from_fit\n")
g <- kl_grouping_from_fit(fit, groupings_dir = gdir)
expect(
  identical(
    setNames(g$partition$group, g$partition$iso3)[ctry],
    setNames(labs[cls], ctry)
  ),
  "argmax partition = the class each country was built in"
)
expect(
  identical(g$partition, dplyr::arrange(g$partition, group, iso3)),
  "partition sorted by (group, iso3)"
)
expect(
  all(c("iso3", "group", "p_max", labs) %in% names(g$membership)) &&
    g$n_fit == 6L &&
    g$n_filled == 0L,
  "membership has p_max + one column per class"
)
g2 <- kl_grouping_from_fit(
  fit,
  chains = c(1L, 3L),
  fill_from = "ref",
  groupings_dir = gdir
)
expect(
  identical(g2$fit$source_chains, c(1L, 3L)) && g2$fit$chains == 2L,
  "chain subset kept (source ids recorded)"
)
expect(
  g2$n_filled == 2L &&
    setequal(setdiff(g2$partition$iso3, ctry), c("GHA", "ZWE")) &&
    g2$partition$group[g2$partition$iso3 == "ZWE"] == "enso_south",
  "fill_from carries the unfitted countries over unchanged"
)

# It reproduces make_kl_grouping.R's former inline body byte for byte.
cat("# [2] same bytes as make_kl_grouping.R's former inline body\n")
old_body <- function(fit, fill) {
  # make_kl_grouping.R before the extraction (:60-78)
  resp <- kl_responsibilities(fit)
  part <- kl_map_partition(fit, "argmax", resp)
  if (nzchar(fill)) {
    reg <- load_grouping(fill, gdir)
    keep <- reg[!reg$iso3 %in% part$iso3, c("iso3", "group")]
    part <- bind_rows(part, keep)
  }
  part %>% arrange(group, iso3)
}
f_old <- file.path(tempdir(), "old_part.csv")
f_new <- file.path(tempdir(), "new_part.csv")
write_kl_partition(old_body(kl_subset_chains(fit, c(1L, 3L)), "ref"), f_old)
write_kl_partition(g2$partition, f_new)
expect(
  identical(readBin(f_old, "raw", 1e5), readBin(f_new, "raw", 1e5)),
  "registry CSV bytes identical"
)

# Shared groups copy meta_from, a new group takes an unused fallback colour;
# the sidecar reads back (load_grouping_meta).
cat("# [3] kl_grouping_sidecar\n")
side <- kl_grouping_sidecar(g$partition, "ref", gdir)
expect(
  identical(unlist(side$order), c("enso_east", "enso_west", "enso_nga")),
  "order: meta_from's for shared groups, then the new ones"
)
expect(
  identical(side$groups$enso_west$label, "West Africa") &&
    identical(side$groups$enso_east$colour, "#1b9e77"),
  "shared groups copy label / colour"
)
expect(
  identical(side$groups$enso_nga$label, "nga") &&
    identical(
      side$groups$enso_nga$colour,
      setdiff(.GROUP_FALLBACK_COLOURS, c("#1b9e77", "#d95f02", "#66a61e"))[1]
    ),
  "a new group: name without prefix + first colour meta_from does not use"
)
expect(identical(side$label_prefix, "enso_"), "label_prefix kept")
out <- write_grouping_registry(
  g$partition,
  "derived",
  side,
  c("# test header"),
  g$membership,
  groupings_dir = gdir,
  membership_dir = file.path(gdir, "memb")
)
meta <- load_grouping_meta("derived", gdir)
expect(
  identical(unname(meta$colour["enso_nga"]), side$groups$enso_nga$colour) &&
    identical(readLines(out$yml, n = 1), "# test header") &&
    file.exists(out$membership_csv),
  "registry + sidecar + membership written and read back"
)

# Centroid names: KL_CENTROID_NAMES fallback, a previous sidecar's kept,
# unknown id -> centroid mode stops naming it.
cat(
  "# [3b] kl_grouping_sidecar: centroid names (results.style.group_labels centroid)\n"
)
expect(
  identical(
    unlist(side$centroid_order),
    c("enso_west", "enso_nga", "enso_east")
  ) &&
    identical(side$groups$enso_nga$centroid_short, "GG") &&
    identical(side$groups$enso_east$centroid_label, "East"),
  "no previous centroid names: KL_CENTROID_NAMES by id, its order"
)
cm <- group_display_meta(meta, mode = "centroid")
expect(
  identical(unname(cm$label), c("West [GW]", "Guinean Gulf [GG]", "East [GE]")),
  "written sidecar works under group_labels centroid"
)
prev <- yaml::read_yaml(out$yml)
prev$groups$enso_east$centroid_short <- "GX"
prev$groups$enso_east$centroid_label <- "Kept"
prev$centroid_order <- list("enso_east", "enso_nga", "enso_west")
writeLines(yaml::as.yaml(prev), out$yml)
side2 <- kl_grouping_sidecar(
  g$partition,
  "ref",
  gdir,
  centroid_from = "derived"
)
expect(
  identical(side2$groups$enso_east$centroid_short, "GX") &&
    identical(side2$groups$enso_east$centroid_label, "Kept") &&
    identical(
      unlist(side2$centroid_order),
      c("enso_east", "enso_nga", "enso_west")
    ),
  "a previous sidecar's centroid fields and order are preserved (ids match)"
)
odd <- g$partition
odd$group[odd$group == "enso_nga"] <- "enso_zzz"
side3 <- kl_grouping_sidecar(odd, "ref", gdir)
expect(
  is.null(side3$groups$enso_zzz$centroid_short) &&
    identical(
      unlist(side3$centroid_order),
      c("enso_west", "enso_east", "enso_zzz")
    ),
  "unknown id: no centroid names, still listed in centroid_order"
)
out3 <- write_grouping_registry(odd, "derived_odd", side3, groupings_dir = gdir)
e3 <- tryCatch(
  {
    group_display_meta(
      load_grouping_meta("derived_odd", gdir),
      mode = "centroid"
    )
    ""
  },
  error = function(e) conditionMessage(e)
)
expect(
  grepl("enso_zzz", e3) && grepl("KL_CENTROID_NAMES", e3),
  "centroid mode on an unknown id stops naming it and KL_CENTROID_NAMES"
)
writeLines(c("# test header", yaml::as.yaml(side)), out$yml) # restore for [4]-[5]

# Absent / identical (row order ignored) / differs.
cat("# [4] kl_grouping_status\n")
expect(
  identical(kl_grouping_status("nope", g$partition, gdir)$status, "absent"),
  "absent"
)
expect(
  identical(
    kl_grouping_status("derived", g$partition[6:1, ], gdir)$status,
    "identical"
  ),
  "identical (row order ignored)"
)
alt <- g$partition
alt$group[alt$iso3 == "SEN"] <- "enso_east"
st <- kl_grouping_status(
  "derived",
  rbind(alt, data.frame(iso3 = "GHA", group = "enso_west")),
  gdir
)
expect(
  identical(st$status, "differs") &&
    setequal(st$diff$iso3, c("SEN", "GHA")) &&
    identical(st$diff$change[st$diff$iso3 == "GHA"], "added"),
  "differs: changed + added rows in the diff"
)

# Missing comparison, empty / unknown cell, agreeing chains, comparison
# settings (n_agree / mode_tol) vs the deploy's, fits.groups check, write then
# re-plan identical.
cat("# [5] kl_grouping_plan_from_sweep\n")
cmp <- file.path(gdir, "cmp")
dir.create(cmp)
err <- function(expr) {
  tryCatch(
    {
      expr
      ""
    },
    error = function(e) conditionMessage(e)
  )
}
expect(
  grepl(
    "comparison incomplete",
    err(kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      name = "x",
      fits_groups = "east",
      groupings_dir = gdir
    ))
  ),
  "no comparison yet -> error naming the sweep command"
)
fit_rds <- file.path(gdir, "fit.rds")
saveRDS(fit, fit_rds)
saveRDS(
  list(
    inv = tibble(
      id = c("ns5_K3", "ns5_K2"),
      fit_file = c(fit_rds, "absent.rds")
    ),
    agree = list(ns5_K3 = c(1L, 2L, 3L), ns5_K2 = 1L),
    settings = list(mode_tol = 1.5, n_agree = 0L)
  ),
  file.path(cmp, "loo.rds")
)
invisible(file.create(file.path(cmp, c("fits.csv", "partition.csv"))))
expect(
  grepl(
    "cell is empty.*ns5_K3",
    err(kl_grouping_plan_from_sweep(
      cmp,
      "",
      name = "x",
      fits_groups = "east",
      groupings_dir = gdir
    ))
  ),
  "empty cell -> error listing the cells"
)
expect(
  grepl(
    "not in the sweep",
    err(kl_grouping_plan_from_sweep(
      cmp,
      "zz_K9",
      name = "x",
      fits_groups = "east",
      groupings_dir = gdir
    ))
  ),
  "unknown cell -> error"
)
pl <- kl_grouping_plan_from_sweep(
  cmp,
  "ns5_K3",
  "agree",
  "",
  "ref",
  "derived2",
  c("east", "west", "nga"),
  gdir
)
expect(
  identical(pl$chains, c(1L, 2L, 3L)) &&
    identical(pl$status, "absent") &&
    isTRUE(pl$groups_ok),
  "agree -> the recorded agreeing chains; short fits.groups resolve"
)
expect(
  identical(
    setNames(pl$grouping$partition$group, pl$grouping$partition$iso3)[ctry],
    setNames(labs[cls], ctry)
  ),
  "chain 4 (label swap) excluded by the agreeing chains"
)
expect(
  identical(
    kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      "agree",
      "",
      "ref",
      "derived2",
      c("east", "west", "nga"),
      gdir,
      n_agree = 0L,
      mode_tol = 1.5
    )$chains,
    c(1L, 2L, 3L)
  ),
  "comparison settings equal to the deploy's -> accepted"
)
e_n <- err(kl_grouping_plan_from_sweep(
  cmp,
  "ns5_K3",
  "agree",
  "",
  "ref",
  "derived2",
  c("east", "west", "nga"),
  gdir,
  n_agree = 2L,
  mode_tol = 1.5
))
expect(
  grepl("n_agree 0 \\(deploy 2\\)", e_n) &&
    grepl("--n_agree 2 --mode_tol 1.5", e_n),
  "comparison run with another n_agree -> refused, naming the re-run command"
)
expect(
  grepl(
    "mode_tol 1.5 \\(deploy 2\\)",
    err(kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      "agree",
      "",
      "ref",
      "derived2",
      c("east", "west", "nga"),
      gdir,
      n_agree = 0L,
      mode_tol = 2
    ))
  ),
  "comparison run with another mode_tol -> refused"
)
expect(
  identical(
    kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      "1,2",
      "",
      "ref",
      "derived2",
      c("east", "west", "nga"),
      gdir,
      n_agree = 2L
    )$chains,
    c(1L, 2L)
  ),
  "explicit chains: the comparison settings are not checked"
)
expect(
  is.null(
    kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      "all",
      "",
      "ref",
      "derived2",
      c("east", "west", "nga"),
      gdir
    )$chains
  ),
  "all -> no subset"
)
expect(
  !kl_grouping_plan_from_sweep(
    cmp,
    "ns5_K3",
    "agree",
    "",
    "ref",
    "derived2",
    c("east", "west"),
    gdir
  )$groups_ok,
  "fits.groups missing a group -> not ok"
)
invisible(write_grouping_registry(
  pl$grouping$partition,
  "derived2",
  pl$sidecar,
  groupings_dir = gdir,
  membership_dir = file.path(gdir, "memb")
))
expect(
  identical(
    kl_grouping_plan_from_sweep(
      cmp,
      "ns5_K3",
      "agree",
      "",
      "ref",
      "derived2",
      c("east", "west", "nga"),
      gdir
    )$status,
    "identical"
  ),
  "write then re-plan -> identical"
)

# Resolve keys (GROUPS "" before the registry exists, KL_SWEEP_*), the sweep
# row + "deferred" in deploy_plan_batches(), deploy_check_batch_groupings(),
# results.kl_sweep_dir alignment target.
cat("# [6] deploy_core: resolve, batch plan, guard, results target\n")
# Fixture: prod_run.yml without a kl sweep (kl_sweep.batch / cell "", the
# comparison settings at their defaults), its grouping an existing registry.
NOKL_DEP <- tempfile("nokl_deploy_", fileext = ".yml")
on.exit(unlink(NOKL_DEP), add = TRUE)
raw0 <- yaml::read_yaml("analysis/configs/deploy/prod_run.yml")
raw0$deploy$kl_sweep$batch <- ""
raw0$deploy$kl_sweep$cell <- ""
raw0$deploy$kl_sweep$n_agree <- 0L
raw0$deploy$kl_sweep$mode_tol <- 1.5
yaml::write_yaml(raw0, NOKL_DEP)
dep <- load_deploy_config(NOKL_DEP)
k0 <- deploy_resolved_keys(dep, NOKL_DEP)
kl_paths <- c(
  "KL_SWEEP_BATCH",
  "KL_SWEEP_DIR",
  "KL_SWEEP_COMPARE_DIR",
  "KL_SWEEP_CELL"
)
expect(
  all(kl_paths %in% names(k0)) &&
    all(unlist(k0[kl_paths]) == "") &&
    nzchar(k0$GROUPS),
  "no sweep: KL_SWEEP_* paths empty, GROUPS from the registry"
)
expect(
  identical(k0$KL_SWEEP_N_AGREE, "0") && identical(k0$KL_SWEEP_MODE_TOL, "1.5"),
  "KL_SWEEP_N_AGREE / KL_SWEEP_MODE_TOL resolve to the defaults"
)
dk <- dep
dk$grouping <- "zz_not_written_yet"
dk$kl_sweep$batch <- "prod_run_kl_sweep_classes"
dk$kl_sweep$cell <- "ns5_K4"
k1 <- deploy_resolved_keys(dk, "x.yml")
expect(
  identical(k1$GROUPS, "") &&
    identical(
      k1$KL_SWEEP_COMPARE_DIR,
      "generated_data/kl_sweep_compare/prod_run_kl_sweep_classes"
    ) &&
    identical(k1$KL_SWEEP_DIR, "analysis/configs/prod_run_kl_sweep_classes") &&
    identical(k1$KL_SWEEP_CELL, "ns5_K4") &&
    identical(
      deploy_resolved_keys(
        within(dk, kl_sweep$n_agree <- 3L),
        "x.yml"
      )$KL_SWEEP_N_AGREE,
      "3"
    ),
  "sweep before the registry exists: GROUPS \"\", KL_SWEEP_* set"
)
dn <- dep
dn$grouping <- "zz_not_written_yet"
expect(
  grepl("registry not found", err(deploy_resolved_keys(dn, "x.yml"))),
  "missing registry without a sweep is still an error"
)
env <- deploy_batch_env()
tmp_cfg <- tempfile("cfgdir_")
dir.create(tmp_cfg)
on.exit(unlink(tmp_cfg, recursive = TRUE), add = TRUE)
assign("get_config_dir", function() tmp_cfg, envir = env)
env$BATCH_REGISTRY$zz_fit_batch <- list(grouping = "zz_not_written_yet")
dk$fits$occ_dir <- dk$fits$dur_dir <- dk$fits$mon_dir <- "zz_fit_batch"
pb <- deploy_plan_batches(dk, env)
expect(
  identical(pb$status[pb$track == "kl"], "generate") &&
    identical(pb$status[pb$track == "occ"], "deferred"),
  "sweep row planned; a fit batch whose registry is absent is deferred"
)
expect(
  isTRUE(deploy_check_batch_groupings(dk, pb)),
  "batch grouping == deploy.grouping passes"
)
dk$grouping <- "enso_grouping"
expect(
  grepl("another registry", err(deploy_check_batch_groupings(dk, pb))),
  "batch grouped by another registry is refused"
)
tmp_dep <- tempfile("klsweep_deploy_", fileext = ".yml")
raw <- yaml::read_yaml(NOKL_DEP)
raw$deploy$kl_sweep$batch <- "prod_run_kl_sweep_classes"
yaml::write_yaml(raw, tmp_dep)
al <- deploy_plan_results_alignment(
  "analysis/configs/results/manuscript.yml",
  tmp_dep
)
expect(
  identical(
    al$targets$kl_sweep_dir,
    "generated_data/kl_sweep_compare/prod_run_kl_sweep_classes"
  ),
  "a sweep deploy aligns results.kl_sweep_dir"
)
expect(
  is.null(
    deploy_plan_results_alignment(
      "analysis/configs/results/manuscript.yml",
      NOKL_DEP
    )$targets$kl_sweep_dir
  ),
  "no sweep: results.kl_sweep_dir left alone"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
