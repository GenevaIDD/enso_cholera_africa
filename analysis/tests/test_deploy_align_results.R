# test_deploy_align_results.R — HARD: the results YAML follows the deploy
#   YAML, and says so.
#
# Usage:
#   Rscript analysis/tests/test_deploy_align_results.R
#
# Details:
#   Covers the deploy<->results coupling in utils/deploy_core.R
#   (align_results_to_deploy()): the line finder and setter [1], the write [2],
#   dry-run and re-align [3], the resolve guard [4] and the CLI dry-run [5].
#   No fit, no data; [5] shells out to Rscript (a few seconds).

suppressPackageStartupMessages({
  library(yaml)
  library(tools)
})

source("analysis/utils/deploy_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# deploy<->results alignment\n")

RESULTS <- "analysis/configs/results/manuscript.yml"
DEPLOYS <- sort(list.files(
  "analysis/configs/deploy",
  pattern = "\\.yml$",
  full.names = TRUE
))
expect(
  length(DEPLOYS) >= 2L,
  "at least two deploy configs exist (need a target != current)"
)
current <- yaml::read_yaml(RESULTS)$results$deploy_config
TARGET <- setdiff(DEPLOYS, current)[1] # some other run to align to
expect(
  !is.na(TARGET) && file.exists(TARGET),
  sprintf("target deploy for the test: %s", TARGET)
)

tmpdir <- tempfile("align_test_")
dir.create(tmpdir)
# One copy per fixture name -- [2]'s aligned copy must survive [3]/[4] making theirs.
tmp_copy <- function(tag) {
  p <- file.path(tmpdir, paste0(tag, "_", basename(RESULTS)))
  file.copy(RESULTS, p, overwrite = TRUE)
  p
}

# --- [1] line finder + scalar setter -----------------------------------------
# The indentation walk picks the right line (never a same-named grandchild)
# and the setter keeps padding + trailing comment; logicals and numbers are
# written bare (a quoted "60" is a string); .deploy_leaf_same() is typed.
cat("\n[1] .yaml_leaf_line() / .yaml_set_scalar_line()\n")
lines <- c(
  "# head comment",
  "results:",
  "  default: \"top\"", # same name as a grandchild
  "  deploy_config: \"a.yml\"   # trailing",
  "  inputs:",
  "    # a comment inside",
  "    manifest:   \"m.yml\"",
  "    lag_loo_batch:",
  "      default:    \"d\"",
  "      occurrence: \"o\"",
  "  scope:",
  "    default: \"scope\""
)
expect(
  identical(.yaml_leaf_line(lines, "deploy_config"), 4L),
  "finds results.deploy_config"
)
expect(
  identical(.yaml_leaf_line(lines, "inputs.manifest"), 7L),
  "finds results.inputs.manifest past a comment line"
)
expect(
  identical(.yaml_leaf_line(lines, "inputs.lag_loo_batch.default"), 9L),
  "finds inputs.lag_loo_batch.default, not results.default nor scope.default"
)
expect(
  identical(.yaml_leaf_line(lines, "default"), 3L),
  "finds results.default (the child, not a grandchild)"
)
expect(is.na(.yaml_leaf_line(lines, "inputs.nope")), "NA for an absent leaf")
expect(
  identical(
    .yaml_set_scalar_line(lines[4], "b.yml"),
    "  deploy_config: \"b.yml\"   # trailing"
  ),
  "setter keeps padding and the trailing comment"
)
expect(
  identical(.yaml_set_scalar_line(lines[7], "x"), "    manifest:   \"x\""),
  "setter keeps alignment padding"
)
err <- tryCatch(
  {
    .yaml_set_scalar_line(lines[8], "x")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(!is.null(err), "setter refuses a mapping (non-scalar) line")
nz_line <- "    episode_national_zero_threshold: \"0\"   # c"
expect(
  identical(
    .yaml_set_scalar_line(nz_line, 60L),
    "    episode_national_zero_threshold: 60   # c"
  ) &&
    identical(
      .yaml_set_scalar_line(nz_line, 2e5),
      "    episode_national_zero_threshold: 200000   # c"
    ) &&
    identical(
      .yaml_set_scalar_line(nz_line, TRUE),
      "    episode_national_zero_threshold: true   # c"
    ),
  "setter writes a number (60L, 2e5) and a logical bare, never quoted"
)
expect(
  identical(
    yaml::yaml.load(sub("#.*$", "", .yaml_set_scalar_line(nz_line, 60L)))[[1]],
    60L
  ),
  "a bare 60 re-reads as the integer 60"
)
expect(
  .deploy_leaf_same(60L, 60L) &&
    .deploy_leaf_same(60, 60L) &&
    !.deploy_leaf_same("60", 60L) &&
    !.deploy_leaf_same(50L, 60L) &&
    !.deploy_leaf_same(NULL, 60L) &&
    .deploy_leaf_same("a", "a") &&
    !.deploy_leaf_same("TRUE", TRUE) &&
    .deploy_leaf_same(TRUE, TRUE),
  ".deploy_leaf_same(): a number must be a number (int == double), a logical a logical"
)
expect(
  grepl("YAML number", .deploy_leaf_type_hint("60", 60L), fixed = TRUE) &&
    grepl("YAML boolean", .deploy_leaf_type_hint("TRUE", TRUE), fixed = TRUE) &&
    identical(.deploy_leaf_type_hint(50L, 60L), "") &&
    identical(.deploy_leaf_type_hint("a", "b"), ""),
  ".deploy_leaf_type_hint(): the unquoting hint for a string number / logical only"
)

# --- [2] write: owned leaves move, nothing else changes ----------------------
# On a copy of manuscript.yml: every owned leaf (deploy_config + the eight
# inputs.*) carries the target deploy's value, every other line is
# byte-identical (comments kept), the pre-edit copy lands under
# generated_data/deploy/<stem>/, and comment lines still naming another run
# are reported (a run id's letter suffix is part of the id).
cat("\n[2] align_results_to_deploy(write = TRUE) on a copy\n")
p <- tmp_copy("aligned")
# A comment naming another run (the shipped file names none), to be reported.
cat("# carried over from prod_run7\n", file = p, append = TRUE)
before <- readLines(p, warn = FALSE)
plan <- capture.output(res <- align_results_to_deploy(p, TARGET, write = TRUE))
after <- readLines(p, warn = FALSE)
dep_derived <- .deploy_derived_results_inputs(
  load_deploy_config(TARGET),
  TARGET
)
raw <- yaml::read_yaml(p)$results

expect(
  nrow(res$changes) >= 1L,
  sprintf("%d line(s) changed", nrow(res$changes))
)
expect(
  identical(raw$deploy_config, TARGET),
  "results.deploy_config now names the target deploy"
)
for (leaf in names(dep_derived)) {
  expect(
    identical(
      .deploy_leaf_str(raw$inputs[[leaf]]),
      .deploy_leaf_str(dep_derived[[leaf]])
    ) &&
      .deploy_leaf_same(raw$inputs[[leaf]], dep_derived[[leaf]]),
    sprintf("inputs.%s carries the deploy-derived value (typed)", leaf)
  )
}
expect(
  is.integer(raw$inputs$episode_national_zero_threshold),
  sprintf(
    "inputs.episode_national_zero_threshold re-reads as an integer (%s)",
    class(raw$inputs$episode_national_zero_threshold)[1]
  )
)
tdep <- load_deploy_config(TARGET)
want_cell <- if (
  nzchar(tdep$kl_sweep$batch %||% "") && nzchar(tdep$kl_sweep$cell %||% "")
) {
  tdep$kl_sweep$cell
} else {
  yaml::read_yaml(RESULTS)$results$kl_sweep_cell
}
expect(
  identical(raw$kl_sweep_cell, want_cell),
  sprintf(
    "results.kl_sweep_cell follows deploy.kl_sweep.cell when the deploy derives its grouping (%s)",
    want_cell
  )
)
# Fixtures: prod_run.yml with its pinned kl cell set to ns5_K5, and without a
# kl sweep (kl_sweep.batch "").
kl_fixture <- function(name, edit) {
  y <- yaml::read_yaml("analysis/configs/deploy/prod_run.yml")
  y$deploy$kl_sweep <- edit(y$deploy$kl_sweep)
  f <- file.path(tmpdir, paste0(name, ".yml"))
  yaml::write_yaml(y, f)
  f
}
p_kl <- deploy_plan_results_alignment(
  RESULTS,
  kl_fixture("deploy_kl_cell", function(k) {
    k$cell <- "ns5_K5"
    k
  })
)
p_nokl <- deploy_plan_results_alignment(
  RESULTS,
  kl_fixture("deploy_no_kl", function(k) {
    k$batch <- ""
    k$cell <- ""
    k
  })
)
expect(
  identical(p_kl$targets$kl_sweep_cell, "ns5_K5") &&
    is.null(p_nokl$targets$kl_sweep_cell),
  "kl_sweep_cell is a target only for a deploy with a kl sweep (pinned ns5_K5; no sweep: none)"
)
expect(length(after) == length(before), "same number of lines")
untouched <- setdiff(seq_along(before), res$changes$line)
expect(
  all(after[untouched] == before[untouched]),
  "every non-owned line is byte-identical"
)
expect(
  sum(grepl("^\\s*#", after)) == sum(grepl("^\\s*#", before)) &&
    sum(grepl("^\\s*#", before)) > 20L,
  sprintf("all %d comment lines kept", sum(grepl("^\\s*#", before)))
)
expect(
  all(grepl("^results\\.", paste0("results.", res$changes$leaf))) &&
    all(res$changes$leaf %in% names(res$targets)),
  "only owned leaves were changed"
)
bak_dir <- make_deploy_run_file(
  list(outdir = "generated_data/"),
  "deploy",
  config_path = TARGET
)
baks <- list.files(
  bak_dir,
  pattern = paste0(
    "^",
    gsub(".", "\\.", basename(p), fixed = TRUE),
    "\\..*\\.bak$"
  ),
  full.names = TRUE
)
expect(length(baks) >= 1L, sprintf("pre-edit copy written under %s", bak_dir))
if (length(baks)) {
  newest <- baks[which.max(file.mtime(baks))]
  expect(
    identical(readLines(newest, warn = FALSE), before),
    "pre-edit copy is the original text"
  )
  invisible(file.remove(newest))
}
tag_n <- sub("^.*run", "", derive_config_run_tag(TARGET))
expect(
  nrow(res$stale_comments) >= 1L &&
    all(grepl("#", after[res$stale_comments$line], fixed = TRUE)),
  sprintf(
    "comment lines still naming a run other than run%s are reported (%d)",
    tag_n,
    nrow(res$stale_comments)
  )
)
tab <- capture.output(cmp <- compare_deploy_results(TARGET, p))
expect(all(cmp$status == "match"), "post-align compare is all-match")

# A run id may carry one letter suffix (a run3b re-forecasting run3): under
# prod_run3b, "run3b" is current and "run3" stale; "run9occ" names run9;
# "prodrun9" (a letter before "run") names no run. Only the deploy path's stem
# is used (the files are hypothetical, never read).
lines_sfx <- c(
  "a: 1  # prod_run3b fits",
  "# run3 and run3b",
  "# lineage run9occ",
  "# prodrun9 only",
  "b: run3  # no run named"
)
st_b <- .results_stale_run_comments(
  lines_sfx,
  "analysis/configs/deploy/prod_run3b.yml"
)
st_3 <- .results_stale_run_comments(
  lines_sfx,
  "analysis/configs/deploy/prod_run3.yml"
)
expect(
  identical(st_b$line, c(2L, 3L)),
  sprintf(
    "run-letter suffix under prod_run3b: stale lines %s",
    paste(st_b$line, collapse = ",")
  )
)
expect(
  identical(st_3$line, c(1L, 2L, 3L)),
  sprintf(
    "run-letter suffix under prod_run3: stale lines %s",
    paste(st_3$line, collapse = ",")
  )
)

# --- [3] dry-run is inert; re-align is a no-op --------------------------------
# write = FALSE (dry-run) changes nothing; a second write is a no-op.
cat("\n[3] dry-run writes nothing; a second write changes nothing\n")
q <- tmp_copy("dry")
md5_0 <- md5sum(q)
invisible(capture.output(
  dry <- align_results_to_deploy(q, TARGET, write = FALSE)
))
expect(
  identical(md5sum(q), md5_0),
  "write = FALSE leaves the file byte-identical"
)
expect(
  nrow(dry$changes) == nrow(res$changes),
  "dry-run plans the same changes the write applied"
)
md5_1 <- md5sum(p)
invisible(capture.output(
  again <- align_results_to_deploy(p, TARGET, write = TRUE)
))
expect(
  nrow(again$changes) == 0L && identical(md5sum(p), md5_1),
  "already-aligned file: no change, no write"
)

# --- [4] the resolve guard -----------------------------------------------------
# An explicit inputs leaf that disagrees with the deploy is refused (naming the
# leaf and the align command); naming it under results.allow_overrides admits
# it; an unknown allow_overrides name is refused; a quoted number or logical
# is refused with the unquoting hint; the shipped results YAMLs resolve.
cat(
  "\n[4] resolve_results_inputs_from_deploy() refuses an undeclared mismatch\n"
)
resolve <- function(path) {
  tryCatch(
    list(
      opt = make_options_from_config(
        option_list = list(),
        cli_argv = c("--config", path)
      ),
      err = NULL
    ),
    error = function(e) list(opt = NULL, err = conditionMessage(e))
  )
}

for (rp in sort(list.files(
  "analysis/configs/results",
  pattern = "\\.yml$",
  full.names = TRUE
))) {
  r <- resolve(rp)
  expect(
    is.null(r$err) &&
      is.integer(r$opt$results$inputs$episode_national_zero_threshold),
    sprintf(
      "shipped %s resolves, zero-fill threshold an integer (%s)",
      rp,
      r$err %||% "no error"
    )
  )
}

g <- tmp_copy("guard") # aligned to `current`; now break one leaf
gl <- readLines(g, warn = FALSE)
i_bl <- .yaml_leaf_line(gl, "inputs.best_lag_dir")
gl[i_bl] <- .yaml_set_scalar_line(gl[i_bl], "generated_data/SOMEWHERE_ELSE")
writeLines(gl, g)
r <- resolve(g)
expect(
  !is.null(r$err) &&
    grepl("inputs.best_lag_dir", r$err, fixed = TRUE) &&
    grepl("--stage align_results", r$err, fixed = TRUE),
  "explicit best_lag_dir != deploy-derived is refused, naming the leaf and the align command"
)
i_ao <- .yaml_leaf_line(gl, "allow_overrides")
expect(!is.na(i_ao), "manuscript.yml carries results.allow_overrides")
gl[i_ao] <- sub("\\[\\]\\s*$", '["best_lag_dir"]', gl[i_ao])
writeLines(gl, g)
r <- resolve(g)
expect(
  is.null(r$err) &&
    identical(
      r$opt$results$inputs$best_lag_dir,
      "generated_data/SOMEWHERE_ELSE"
    ),
  sprintf(
    "the same mismatch is admitted once named under allow_overrides (%s)",
    r$err %||% "ok"
  )
)
gl[i_ao] <- sub('\\["best_lag_dir"\\]', '["not_a_leaf"]', gl[i_ao])
writeLines(gl, g)
r <- resolve(g)
expect(
  !is.null(r$err) && grepl("not_a_leaf", r$err, fixed = TRUE),
  "an unknown allow_overrides name is refused"
)


# A logical leaf written as a quoted string prints like the deploy's logical but
# reads as FALSE through isTRUE(): refused, with the fix in the message.
q <- tmp_copy("guard_logical")
ql <- readLines(q, warn = FALSE)
i_nz <- .yaml_leaf_line(ql, "inputs.episode_national_zero_fill")
nz_dep <- isTRUE(
  load_deploy_config(
    yaml::read_yaml(q)$results$deploy_config
  )$episodes$national_zero_fill
)
ql[i_nz] <- .yaml_set_scalar_line(ql[i_nz], toupper(as.character(nz_dep))) # e.g. "TRUE" (quoted)
writeLines(ql, q)
expect(
  grepl('"', ql[i_nz], fixed = TRUE),
  "test setup: the zero-fill leaf is now a quoted string"
)
r <- resolve(q)
expect(
  !is.null(r$err) &&
    grepl("inputs.episode_national_zero_fill", r$err, fixed = TRUE) &&
    grepl("YAML boolean", r$err, fixed = TRUE),
  "a quoted logical leaf equal in print to the deploy's is refused as a string"
)

# Same for a number: a quoted "60" against the deploy's 60 is refused,
# with the unquoting hint, and compare_deploy_results() calls it a MISMATCH.
q <- tmp_copy("guard_number")
ql <- readLines(q, warn = FALSE)
i_nt <- .yaml_leaf_line(ql, "inputs.episode_national_zero_threshold")
q_dep <- yaml::read_yaml(q)$results$deploy_config
nt_dep <- .deploy_derived_results_inputs(
  load_deploy_config(q_dep),
  q_dep
)$episode_national_zero_threshold
ql[i_nt] <- .yaml_set_scalar_line(ql[i_nt], as.character(nt_dep)) # e.g. "60" (quoted)
writeLines(ql, q)
expect(
  grepl('"', ql[i_nt], fixed = TRUE),
  sprintf(
    "test setup: the threshold leaf is now a quoted string (%s)",
    trimws(ql[i_nt])
  )
)
r <- resolve(q)
expect(
  !is.null(r$err) &&
    grepl("inputs.episode_national_zero_threshold", r$err, fixed = TRUE) &&
    grepl("YAML number", r$err, fixed = TRUE),
  "a quoted number equal in print to the deploy's is refused, with the unquoting hint"
)
invisible(capture.output(cq <- compare_deploy_results(q_dep, q)))
expect(
  identical(
    cq$status[cq$leaf == "episode_national_zero_threshold"],
    "MISMATCH"
  ) &&
    all(cq$status[cq$leaf != "episode_national_zero_threshold"] == "match"),
  "compare_deploy_results(): the quoted number is the only MISMATCH (typed status)"
)

# --- [4b] episode key: deploy-derived, copied to the top-level episodes.key ----
cat(
  "\n[4b] results.inputs.episode_key follows the deploy and drives episodes.key\n"
)
r <- resolve(RESULTS)
dep_key <- load_deploy_config(r$opt$results$deploy_config)$episodes$key %||%
  "lp_season"
expect(
  identical(r$opt$results$inputs$episode_key, dep_key) &&
    identical(r$opt$episodes$key, dep_key),
  sprintf(
    "shipped config: inputs.episode_key == episodes.key == deploy.episodes.key (%s)",
    dep_key
  )
)
e <- tmp_copy("epkey")
el <- readLines(e, warn = FALSE)
writeLines(c(el, "episodes:", "  key: \"lp_season\""), e)
r <- resolve(e)
expect(
  !is.null(r$err) && grepl("top-level episodes.key", r$err, fixed = TRUE),
  "a top-level episodes.key disagreeing with results.inputs.episode_key is refused"
)
i_ek <- .yaml_leaf_line(el, "inputs.episode_key")
expect(!is.na(i_ek), "manuscript.yml carries results.inputs.episode_key")
el[i_ek] <- .yaml_set_scalar_line(el[i_ek], "lp_season")
writeLines(el, e)
r <- resolve(e)
expect(
  !is.null(r$err) && grepl("inputs.episode_key", r$err, fixed = TRUE),
  "an explicit episode_key != deploy.episodes.key is refused"
)

nr <- file.path(tmpdir, "not_results.yml")
writeLines(
  c(
    "episodes:",
    "  key: \"gid_season\"",
    "results:",
    "  series: latent",
    "  inputs:",
    "    manifest: analysis/configs/forecast_horizons.yml"
  ),
  nr
)
r <- resolve(nr)
expect(
  is.null(r$err) && identical(r$opt$episodes$key, "gid_season"),
  "a config not stating results.inputs.episode_key (e.g. a dumped results: block) keeps its own episodes.key"
)

# --- [5] CLI dry-run on the shipped files --------------------------------------
# On the shipped files: exits 0 and leaves manuscript.yml byte-identical.
cat("\n[5] deploy_prod_run.R --stage align_results --dry-run\n")
md5_ship <- md5sum(RESULTS)
cli <- run_rscript(
  "analysis/utils/deploy_prod_run.R",
  c("-C", current, "--stage", "align_results", "--dry-run")
)
expect(cli$status == 0L, sprintf("exit status 0 (got %d)", cli$status))
expect(any(grepl("dry-run", cli$out, fixed = TRUE)), "reports the dry-run")
expect(
  identical(md5sum(RESULTS), md5_ship),
  "shipped manuscript.yml left byte-identical"
)

unlink(tmpdir, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
