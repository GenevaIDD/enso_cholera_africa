# make_pipeline_docs.R — render the pipeline run graph from
#   analysis/pipeline_steps.yml
#
# Usage:
#   Rscript analysis/utils/make_pipeline_docs.R          # rewrite both outputs
#   Rscript analysis/utils/make_pipeline_docs.R --check  # exit 1 if stale
#   dot -Tsvg notebooks/pipeline_run_graph.dot \
#     -o notebooks/pipeline_run_graph.svg
#
# Details:
#   One table feeds two outputs, so they cannot drift apart: the step table
#   between the GENERATED markers in analysis/PIPELINE.md (the hand-written
#   prose around it is left untouched) and notebooks/pipeline_run_graph.dot, a
#   file-level graphviz run graph. notebooks/pipeline_diagram.dot is
#   hand-maintained and not touched here. Can also be sourced
#   (tests/test_pipeline_doc_coverage.R does) for pipeline_steps_load(),
#   pipeline_validate() and pipeline_render_*().

suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

source("analysis/utils/base.R") # %||%

PIPELINE_STEPS_FILE <- "analysis/pipeline_steps.yml"
PIPELINE_MD_FILE <- "analysis/PIPELINE.md"
PIPELINE_DOT_FILE <- "notebooks/pipeline_run_graph.dot"
PIPELINE_MD_BEGIN <- "<!-- BEGIN GENERATED: pipeline_steps.yml via analysis/utils/make_pipeline_docs.R -- do not edit by hand -->"
PIPELINE_MD_END <- "<!-- END GENERATED -->"

.PIPELINE_ENUMS <- list(
  kind = c("step", "tool", "report"),
  track = c(
    "data",
    "pop",
    "occ",
    "siz",
    "dur",
    "mon",
    "clim",
    "all",
    "sim",
    "results"
  ),
  rollout = c("core", "optional"),
  runs = c("local", "slurm", "local+slurm")
)

# ---------------------------------------------------------------------------
# Load + validate
# ---------------------------------------------------------------------------

#' Read pipeline_steps.yml into a list of steps with every optional field filled.
pipeline_steps_load <- function(path = PIPELINE_STEPS_FILE) {
  raw <- yaml::read_yaml(path)$steps
  lapply(raw, function(s) {
    s$after <- as.character(unlist(s$after %||% character(0)))
    s$drivers <- as.character(unlist(s$drivers %||% character(0)))
    s$inputs <- as.character(unlist(s$inputs %||% character(0)))
    s$outputs <- as.character(unlist(s$outputs %||% character(0)))
    s$dispatched_by <- s$dispatched_by %||% NA_character_
    s$note <- s$note %||% ""
    s
  })
}

#' Structural problems with the step list (character(0) when clean).
#' Paths are relative to the repository root (the working directory).
pipeline_validate <- function(steps) {
  probs <- character(0)
  ids <- vapply(steps, `[[`, "", "id")
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) {
    probs <- c(probs, sprintf("duplicate id: %s", dup))
  }
  for (i in seq_along(steps)) {
    s <- steps[[i]]
    for (f in c("id", "file", "kind", "phase", "track", "rollout", "runs")) {
      if (is.null(s[[f]]) || !nzchar(s[[f]])) {
        probs <- c(probs, sprintf("%s: missing field '%s'", s$id, f))
      }
    }
    for (f in names(.PIPELINE_ENUMS)) {
      if (!is.null(s[[f]]) && !s[[f]] %in% .PIPELINE_ENUMS[[f]]) {
        probs <- c(
          probs,
          sprintf(
            "%s: %s '%s' not one of %s",
            s$id,
            f,
            s[[f]],
            paste(.PIPELINE_ENUMS[[f]], collapse = "|")
          )
        )
      }
    }
    if (!file.exists(file.path("analysis", s$file))) {
      probs <- c(
        probs,
        sprintf("%s: file analysis/%s does not exist", s$id, s$file)
      )
    }
    for (d in s$drivers) {
      if (!file.exists(file.path("analysis/bash", d))) {
        probs <- c(
          probs,
          sprintf("%s: driver analysis/bash/%s does not exist", s$id, d)
        )
      }
    }
    # List order is run order: everything a step needs must come earlier.
    earlier <- ids[seq_len(i - 1L)]
    for (a in c(s$after, stats::na.omit(s$dispatched_by))) {
      if (!a %in% ids) {
        probs <- c(probs, sprintf("%s: unknown dependency '%s'", s$id, a))
      } else if (!a %in% earlier) {
        probs <- c(
          probs,
          sprintf("%s: depends on '%s', which comes later in the list", s$id, a)
        )
      }
    }
  }
  probs
}

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

.md_cell <- function(x) gsub("|", "\\|", x, fixed = TRUE)

# Free text (outputs / inputs / notes): `<placeholder>` would be eaten as an HTML
# tag and a `*` glob would open an emphasis span.
.md_text <- function(x) {
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("*", "\\*", x, fixed = TRUE)
}

.step_label_md <- function(s) {
  if (s$kind == "tool" && startsWith(s$id, "deploy:")) {
    sprintf("`%s --stage %s`", s$file, sub("^deploy:", "", s$id))
  } else {
    sprintf("`%s`", s$file)
  }
}

#' Markdown lines for the generated block of PIPELINE.md.
pipeline_render_table <- function(steps) {
  hdr <- c(
    "| # | Step | Phase | Track | Rollout | Runs | After | Driver(s) | Writes |",
    "|---|---|---|---|---|---|---|---|---|"
  )
  rows <- vapply(
    seq_along(steps),
    function(i) {
      s <- steps[[i]]
      after <- c(
        s$after,
        if (!is.na(s$dispatched_by)) sprintf("(dispatch: %s)", s$dispatched_by)
      )
      cells <- c(
        as.character(i),
        .step_label_md(s),
        s$phase,
        s$track,
        s$rollout,
        s$runs,
        if (length(after)) paste(after, collapse = ", ") else "—",
        if (length(s$drivers)) {
          paste(sprintf("`%s`", s$drivers), collapse = ", ")
        } else {
          "—"
        },
        paste(
          c(
            .md_text(s$outputs),
            if (length(s$inputs)) sprintf("*input:* %s", .md_text(s$inputs)),
            if (nzchar(s$note)) sprintf("*note:* %s", .md_text(s$note))
          ),
          collapse = "<br>"
        )
      )
      paste0("| ", paste(.md_cell(cells), collapse = " | "), " |")
    },
    ""
  )

  drivers <- sort(unique(unlist(lapply(steps, `[[`, "drivers"))))
  drv_rows <- vapply(
    drivers,
    function(d) {
      who <- vapply(Filter(function(s) d %in% s$drivers, steps), `[[`, "", "id")
      sprintf("| `%s` | %s |", d, paste(who, collapse = ", "))
    },
    ""
  )

  c(
    PIPELINE_MD_BEGIN,
    "",
    sprintf("### Step table (%d entries, in run order)", length(steps)),
    "",
    hdr,
    rows,
    "",
    "### Bash driver → steps",
    "",
    "| Driver (`analysis/bash/`) | Runs |",
    "|---|---|",
    unname(drv_rows),
    "",
    PIPELINE_MD_END
  )
}

.dot_id <- function(id) paste0("n_", gsub("[^A-Za-z0-9]", "_", id))

.dot_label <- function(s) {
  lab <- if (startsWith(s$id, "deploy:")) {
    sprintf("deploy_prod_run.R\\n--stage %s", sub("^deploy:", "", s$id))
  } else {
    sub("\\.(R|qmd)$", "", basename(s$file))
  }
  if (s$runs != "local") {
    lab <- sprintf("%s\\n[%s]", lab, s$runs)
  }
  lab
}

.TRACK_FILL <- c(
  data = "#e8e8e8",
  pop = "#e3d7ee",
  occ = "#d6e6f5",
  siz = "#f8dfc8",
  dur = "#d9eed6",
  mon = "#f5d6dc",
  clim = "#e6f0c8",
  all = "#f2f0d0",
  sim = "#d3eeee",
  results = "#ffffff"
)

.PHASE_TITLE <- c(
  "00" = "00 ingest",
  "01" = "01 preprocess outbreaks",
  "02" = "02 preprocess population / climate",
  "03" = "03 prepare",
  "04" = "04 fit",
  "05" = "05 fit postprocess",
  "06" = "06 fit figures",
  "07" = "07 simulator",
  "08" = "08 smoothing (before 07b spatPomp)",
  "09" = "09 postprocess results",
  "10" = "10 manuscript figures",
  "11" = "11 manuscript numbers",
  "12" = "12 country reports"
)

#' Graphviz lines for notebooks/pipeline_run_graph.dot.
pipeline_render_dot <- function(steps) {
  node_line <- function(s, indent) {
    style <- if (s$rollout == "core") {
      "\"rounded,filled\""
    } else {
      "\"rounded,filled,dashed\""
    }
    shape <- switch(s$kind, tool = "hexagon", report = "note", "box")
    fill <- if (s$rollout == "core") .TRACK_FILL[[s$track]] else "#f7f7f7"
    sprintf(
      "%s%s [label=\"%s\", shape=%s, style=%s, fillcolor=\"%s\"];",
      indent,
      .dot_id(s$id),
      .dot_label(s),
      shape,
      style,
      fill
    )
  }
  out <- c(
    "// GENERATED by analysis/utils/make_pipeline_docs.R from analysis/pipeline_steps.yml.",
    "// Do not edit by hand. Render: dot -Tsvg notebooks/pipeline_run_graph.dot -o pipeline_run_graph.svg",
    "// Solid edge = consumes artefacts of; dashed edge = dispatcher routes to.",
    "// Dashed node = optional (not in a prod_run rollout); hexagon = non-step tool; fill = track.",
    "digraph pipeline_run_graph {",
    "  rankdir=TB;",
    "  fontname=\"Helvetica\";",
    "  node [fontname=\"Helvetica\", fontsize=11, color=\"#555555\"];",
    "  edge [color=\"#555555\", arrowsize=0.7];",
    ""
  )
  phases <- unique(vapply(steps, `[[`, "", "phase"))
  for (p in phases) {
    members <- Filter(function(s) s$phase == p, steps)
    if (p == "--") {
      out <- c(
        out,
        "  // non-step tools (utils/), placed by their edges",
        vapply(members, node_line, "", indent = "  "),
        ""
      )
    } else {
      out <- c(
        out,
        sprintf("  subgraph cluster_%s {", p),
        sprintf(
          "    label=\"%s\"; style=\"rounded\"; color=\"#999999\"; fontsize=12;",
          .PHASE_TITLE[[p]]
        ),
        vapply(members, node_line, "", indent = "    "),
        "  }",
        ""
      )
    }
  }
  for (s in steps) {
    for (a in s$after) {
      out <- c(out, sprintf("  %s -> %s;", .dot_id(a), .dot_id(s$id)))
    }
    if (!is.na(s$dispatched_by)) {
      out <- c(
        out,
        sprintf(
          "  %s -> %s [style=dashed];",
          .dot_id(s$dispatched_by),
          .dot_id(s$id)
        )
      )
    }
  }
  c(out, "}")
}

#' Replace the GENERATED block of PIPELINE.md's lines with `block`.
pipeline_splice_md <- function(md_lines, block) {
  b <- which(md_lines == PIPELINE_MD_BEGIN)
  e <- which(md_lines == PIPELINE_MD_END)
  if (length(b) != 1L || length(e) != 1L || e < b) {
    stop(
      "[make_pipeline_docs] ",
      PIPELINE_MD_FILE,
      " must contain exactly one BEGIN/END GENERATED marker pair"
    )
  }
  c(md_lines[seq_len(b - 1L)], block, md_lines[-seq_len(e)])
}

#' The two outputs as they should be on disk: list(md = <lines>, dot = <lines>).
pipeline_render_all <- function(steps = pipeline_steps_load()) {
  list(
    md = pipeline_splice_md(
      readLines(PIPELINE_MD_FILE, warn = FALSE),
      pipeline_render_table(steps)
    ),
    dot = pipeline_render_dot(steps)
  )
}

# ---------------------------------------------------------------------------
# CLI (Rscript top level only; sourcing defines the functions)
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  opt <- parse_args(OptionParser(
    option_list = list(
      make_option(
        "--check",
        action = "store_true",
        default = FALSE,
        help = "Do not write; exit 1 if PIPELINE.md or the .dot is stale"
      )
    )
  ))

  steps <- pipeline_steps_load()
  probs <- pipeline_validate(steps)
  if (length(probs)) {
    stop(
      "[make_pipeline_docs] ",
      PIPELINE_STEPS_FILE,
      ":\n  ",
      paste(probs, collapse = "\n  ")
    )
  }

  want <- pipeline_render_all(steps)
  have <- list(
    md = readLines(PIPELINE_MD_FILE, warn = FALSE),
    dot = if (file.exists(PIPELINE_DOT_FILE)) {
      readLines(PIPELINE_DOT_FILE, warn = FALSE)
    } else {
      character(0)
    }
  )
  files <- c(md = PIPELINE_MD_FILE, dot = PIPELINE_DOT_FILE)
  stale <- names(files)[
    !mapply(identical, want[names(files)], have[names(files)])
  ]

  if (opt$check) {
    if (length(stale)) {
      message(
        "[make_pipeline_docs] stale: ",
        paste(files[stale], collapse = ", "),
        " -- run: Rscript analysis/utils/make_pipeline_docs.R"
      )
      quit(status = 1)
    }
    message("[make_pipeline_docs] up to date (", length(steps), " entries)")
  } else {
    for (k in stale) {
      writeLines(want[[k]], files[[k]])
    }
    message(
      "[make_pipeline_docs] ",
      if (length(stale)) {
        paste("wrote", paste(files[stale], collapse = ", "))
      } else {
        "nothing to write"
      },
      " (",
      length(steps),
      " entries)"
    )
  }
}
