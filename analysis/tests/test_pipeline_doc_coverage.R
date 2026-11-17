# test_pipeline_doc_coverage.R — HARD: analysis/PIPELINE.md stays the
#   complete run-order record.
#
# Usage:
#   Rscript analysis/tests/test_pipeline_doc_coverage.R
#
# Details:
#   There is no executed entry point (no Makefile, no targets, no umbrella
#   bash script, by design), so the run order lives only in
#   analysis/pipeline_steps.yml and the files generated from it (PIPELINE.md's
#   step table, notebooks/pipeline_run_graph.dot). This test stands in for
#   the ordering an entry point would have enforced. Pure text/YAML -- no
#   fit, no data, no network.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/make_pipeline_docs.R") # functions only (CLI guarded by sys.nframe())

cat("# pipeline_steps.yml / PIPELINE.md cover the whole pipeline\n")

steps <- pipeline_steps_load()
ids <- vapply(steps, `[[`, "", "id")
files <- vapply(steps, `[[`, "", "file")
kinds <- vapply(steps, `[[`, "", "kind")
expect(length(steps) > 40L, sprintf("loaded %d entries", length(steps)))

# 1. pipeline_steps.yml is structurally valid: unique ids, known enum values,
#    every file / driver exists, and the list order is a valid run order
#    (every `after` / `dispatched_by` entry appears earlier).
# 1. structure ------------------------------------------------------------------
probs <- pipeline_validate(steps)
expect(
  length(probs) == 0L,
  sprintf(
    "pipeline_steps.yml is valid (problems: %s)",
    if (length(probs)) paste(probs, collapse = "; ") else "none"
  )
)

# 2. every numbered script and the numbered .qmd report has an entry, and
#    every step/report entry is a numbered file.
# 2. numbered scripts <-> entries -----------------------------------------------
numbered <- list.files("analysis", pattern = "^[0-9][0-9].*\\.(R|qmd)$")
missing <- setdiff(numbered, files[kinds %in% c("step", "report")])
expect(
  length(missing) == 0L,
  sprintf(
    "every numbered script/report has an entry (missing: %s)",
    if (length(missing)) paste(missing, collapse = ", ") else "none"
  )
)
not_numbered <- files[
  kinds %in% c("step", "report") & !grepl("^[0-9][0-9]", files)
]
expect(
  length(not_numbered) == 0L,
  sprintf(
    "every step/report entry is a numbered file (offenders: %s)",
    if (length(not_numbered)) paste(not_numbered, collapse = ", ") else "none"
  )
)
step_ids_ok <- all(ids[kinds != "tool"] == files[kinds != "tool"])
expect(step_ids_ok, "a step/report entry's id is its filename")

# 3. every tier-1 stage driver (analysis/bash/stages/*.sh) drives some entry.
# 3. tier-1 stage drivers -------------------------------------------------------
stage_drivers <- list.files(
  "analysis/bash/stages",
  pattern = "\\.sh$",
  full.names = TRUE
)
stage_drivers <- stage_drivers[!nzchar(Sys.readlink(stage_drivers))]
stage_drivers <- file.path("stages", basename(stage_drivers))
expect(
  length(stage_drivers) > 5L,
  sprintf("found %d tier-1 stage drivers", length(stage_drivers))
)
used_drivers <- unique(unlist(lapply(steps, `[[`, "drivers")))
undriven <- setdiff(stage_drivers, used_drivers)
expect(
  length(undriven) == 0L,
  sprintf(
    "every stages/*.sh drives some entry (missing: %s)",
    if (length(undriven)) paste(undriven, collapse = ", ") else "none"
  )
)

# 4. every deploy_prod_run.R stage except the composites (all, resolve) has a
#    `deploy:<stage>` entry.
# 4. deploy tool stages ---------------------------------------------------------
pd <- parse("analysis/utils/deploy_prod_run.R", keep.source = FALSE)
vs <- Filter(
  function(e) {
    is.call(e) &&
      identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("VALID_STAGES"))
  },
  as.list(pd)
)
expect(length(vs) == 1L, "deploy_prod_run.R defines VALID_STAGES once")
if (length(vs) == 1L) {
  deploy_stages <- setdiff(eval(vs[[1]][[3]]), c("all", "resolve"))
  missing_stage <- setdiff(paste0("deploy:", deploy_stages), ids)
  expect(
    length(missing_stage) == 0L,
    sprintf(
      "every deploy stage has an entry (missing: %s)",
      if (length(missing_stage)) {
        paste(missing_stage, collapse = ", ")
      } else {
        "none"
      }
    )
  )
  unknown_stage <- setdiff(
    grep("^deploy:", ids, value = TRUE),
    paste0("deploy:", deploy_stages)
  )
  expect(
    length(unknown_stage) == 0L,
    sprintf(
      "no deploy:<stage> entry for a stage the tool lacks (stale: %s)",
      if (length(unknown_stage)) {
        paste(unknown_stage, collapse = ", ")
      } else {
        "none"
      }
    )
  )
}

# 5. PIPELINE.md names every numbered script and every stage driver.
# 5. PIPELINE.md names everything -----------------------------------------------
md <- paste(readLines(PIPELINE_MD_FILE, warn = FALSE), collapse = "\n")
unnamed <- c(numbered, stage_drivers)[
  !vapply(
    c(numbered, stage_drivers),
    function(x) grepl(x, md, fixed = TRUE),
    TRUE
  )
]
expect(
  length(unnamed) == 0L,
  sprintf(
    "PIPELINE.md names every numbered script and stage driver (missing: %s)",
    if (length(unnamed)) paste(unnamed, collapse = ", ") else "none"
  )
)

# 6. the generated block of PIPELINE.md and the .dot are up to date.
# 6. generated outputs are fresh ------------------------------------------------
want <- pipeline_render_all(steps)
expect(
  identical(want$md, readLines(PIPELINE_MD_FILE, warn = FALSE)),
  "PIPELINE.md's generated block is up to date (else: Rscript analysis/utils/make_pipeline_docs.R)"
)
expect(
  file.exists(PIPELINE_DOT_FILE) &&
    identical(want$dot, readLines(PIPELINE_DOT_FILE, warn = FALSE)),
  "notebooks/pipeline_run_graph.dot is up to date (else: Rscript analysis/utils/make_pipeline_docs.R)"
)

# The validator itself catches an out-of-order dependency -----------------------
swapped <- steps
i <- match("08a_smoothing_tasks.R", ids)
j <- match("08b_smoothing_spatpomp.R", ids)
swapped[c(i, j)] <- swapped[c(j, i)]
expect(
  any(grepl("comes later in the list", pipeline_validate(swapped))),
  "pipeline_validate() flags a step listed before its dependency"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
