# test_figure_path_collisions.R — SOFT: no two fit configs write the same
#   figure files.
#
# Usage:
#   Rscript analysis/tests/test_figure_path_collisions.R
#
# Details:
#   make_figure_dir() (utils/paths.R) keys 06_fit_figures.R's output directory
#   on the config's enclosing directory as well as its basename
#   (figures/<track>/<batch>/), so two configs collide only if two config
#   directories share a basename AND hold configs sharing a basename; the
#   second run's figures would then silently clobber the first's. Ratchet, not
#   a hard invariant: an allowlist (empty) that can only shrink.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

# Ratchet: any new collision fails; any allowlist entry that stops colliding
# must be deleted (a stale entry also fails), so the count can only go down.
# (figure dir, stem) keys known to be shared by configs with genuinely
# different underlying fits. Empty at introduction.
KNOWN_COLLISIONS <- character(0)

cfgs <- list.files(
  "analysis/configs",
  pattern = "\\.ya?ml$",
  recursive = TRUE,
  full.names = TRUE
)
cfgs <- cfgs[basename(cfgs) != "manifest.yml"] # documentation, not a fit config

rows <- lapply(cfgs, function(cfg) {
  raw <- tryCatch(yaml::read_yaml(cfg), error = function(e) NULL)
  if (is.null(raw)) {
    return(NULL)
  }

  model <- if (!is.null(raw$occurrence)) {
    "occurrence"
  } else if (!is.null(raw$size)) {
    "size"
  } else if (!is.null(raw$duration)) {
    "duration"
  } else {
    NA_character_
  }
  if (is.na(model)) {
    return(NULL)
  } # e.g. a simulator-only or grouping config

  opt <- tryCatch(
    suppressMessages(make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", cfg),
      model = model
    )),
    error = function(e) NULL
  )
  if (is.null(opt)) {
    return(NULL)
  }

  # Exactly 06_fit_figures.R's own tag/fig_dir computation (config_path =
  # dirname(cfg), so <batch> is the config's ENCLOSING directory, not the
  # config's own basename).
  tag <- tools::file_path_sans_ext(basename(cfg))
  fig_dir <- make_figure_dir(
    opt = list(figdir = NULL),
    track = model,
    config_path = dirname(cfg)
  )
  fit_kind <- switch(
    model,
    occurrence = "fit",
    size = "size_fit",
    duration = "dur_fit"
  )
  fit_path <- tryCatch(
    make_run_file(opt, kind = fit_kind),
    error = function(e) NA_character_
  )

  data.frame(
    config = cfg,
    model = model,
    key = file.path(fig_dir, tag),
    fit_path = fit_path,
    stringsAsFactors = FALSE
  )
})
rows <- rows[!vapply(rows, is.null, logical(1))]
rows <- do.call(rbind, rows)

cat(sprintf(
  "# figure-path collision gate: %d fit configs resolved (of %d yaml files found)\n",
  nrow(rows),
  length(cfgs)
))

by_key <- split(rows, rows$key)
collisions <- Filter(function(d) length(unique(d$fit_path)) > 1L, by_key)
collision_keys <- sort(names(collisions))

cat(sprintf("\n  distinct (figure dir, stem) keys : %d\n", length(by_key)))
cat(sprintf(
  "  keys with >1 underlying fit      : %d (allowlist %d)\n",
  length(collision_keys),
  length(KNOWN_COLLISIONS)
))

# --- 1. no NEW collision -----------------------------------------------------
# Every fit config under analysis/configs/** gets its (make_figure_dir(), stem)
# key exactly as 06_fit_figures.R computes it; a key shared by configs
# resolving to DIFFERENT underlying fits (make_run_file(opt, "fit" |
# "size_fit" | "dur_fit")) is flagged. A key shared by configs that resolve to
# the SAME fit is not: that is a harmless re-run of the same result under two
# names, not a silent overwrite of two different ones.
cat("\n[1] no new collision\n")

new_collisions <- setdiff(collision_keys, KNOWN_COLLISIONS)
if (length(new_collisions)) {
  cat("        new collisions (different fits sharing one figure path):\n")
  for (k in new_collisions) {
    d <- collisions[[k]]
    cat(sprintf("          - %s\n", k))
    for (i in seq_len(nrow(d))) {
      cat(sprintf("              %s -> %s\n", d$config[i], d$fit_path[i]))
    }
  }
}
expect(
  length(new_collisions) == 0L,
  "no config pair collides on (make_figure_dir, stem) outside KNOWN_COLLISIONS"
)

# --- 2. allowlist is not stale ------------------------------------------------
cat("\n[2] allowlist is not stale\n")

stale <- setdiff(KNOWN_COLLISIONS, collision_keys)
if (length(stale)) {
  cat("        repaired -- delete from KNOWN_COLLISIONS:\n")
  cat(paste0("          - ", stale, "\n"), sep = "")
}
expect(length(stale) == 0L, "KNOWN_COLLISIONS carries no already-repaired key")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
