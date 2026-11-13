# test_figure_path_collisions.R — figure-path collision gate.
#
# SOFT gate (see analysis/tests/run_all.R).
#
# WHY THIS EXISTS
# ---------------
# Naming the output directory of 06_fit_figures.R from the
# CONFIG'S OWN BASENAME alone (figures/<track>/, filenames prefixed
# "<config-basename>_...") lets two configs sharing a basename in different batch
# directories (e.g. set_A/run_brms.yml and set_B/run_brms.yml) write the
# SAME figure files -- the second run's figures silently clobber the first
# run's, with no error. makeFigureDir() (utils/paths.R) prevents this by
# keying the output directory on the config's ENCLOSING directory as well
# (figures/<track>/<batch>/), so the only way two configs can still collide is
# if two DIFFERENT config directories anywhere under analysis/configs/ share a
# basename AND the configs inside them also share a basename.
#
# This test enumerates every fit config under analysis/configs/**, computes
# each one's (makeFigureDir(), stem) key exactly as 06_fit_figures.R does, and
# flags any key that is shared by configs resolving to DIFFERENT underlying
# fits (makeRunFile(opt, "fit" | "size_fit" | "dur_fit")). A key shared by
# configs that resolve to the SAME fit is not flagged -- that is a harmless
# re-run of the same result under two names, not a silent overwrite of two
# different ones.
#
# Ratchet, not a hard invariant: the check is ALLOWLISTED, and the allowlist is
# empty (0 collisions across all 594 configs at introduction). Any new
# collision fails; any allowlist entry that stops colliding must be deleted (a
# stale allowlist entry also fails), so the count can only go down.
#
# Run: Rscript analysis/tests/test_figure_path_collisions.R

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

# (figure dir, stem) keys known to be shared by configs with genuinely
# different underlying fits. Empty at introduction.
KNOWN_COLLISIONS <- character(0)

cfgs <- list.files("analysis/configs", pattern = "\\.ya?ml$",
                    recursive = TRUE, full.names = TRUE)
cfgs <- cfgs[basename(cfgs) != "manifest.yml"]  # documentation, not a fit config

rows <- lapply(cfgs, function(cfg) {
  raw <- tryCatch(yaml::read_yaml(cfg), error = function(e) NULL)
  if (is.null(raw)) return(NULL)

  model <- if (!is.null(raw$occurrence))    "occurrence"
           else if (!is.null(raw$size))     "size"
           else if (!is.null(raw$duration)) "duration"
           else NA_character_
  if (is.na(model)) return(NULL)  # e.g. a simulator-only or grouping config

  opt <- tryCatch(
    suppressMessages(make_options_from_config(option_list = list(),
                                               cli_argv    = c("--config", cfg),
                                               model       = model)),
    error = function(e) NULL)
  if (is.null(opt)) return(NULL)

  # Exactly 06_fit_figures.R's own tag/fig_dir computation (config_path =
  # dirname(cfg), so <batch> is the config's ENCLOSING directory, not the
  # config's own basename).
  tag      <- tools::file_path_sans_ext(basename(cfg))
  fig_dir  <- makeFigureDir(opt = list(figdir = NULL), track = model,
                            config_path = dirname(cfg))
  fit_kind <- switch(model, occurrence = "fit", size = "size_fit", duration = "dur_fit")
  fit_path <- tryCatch(makeRunFile(opt, kind = fit_kind), error = function(e) NA_character_)

  data.frame(config = cfg, model = model, key = file.path(fig_dir, tag),
             fit_path = fit_path, stringsAsFactors = FALSE)
})
rows <- rows[!vapply(rows, is.null, logical(1))]
rows <- do.call(rbind, rows)

cat(sprintf("# figure-path collision gate: %d fit configs resolved (of %d yaml files found)\n",
            nrow(rows), length(cfgs)))

by_key     <- split(rows, rows$key)
collisions <- Filter(function(d) length(unique(d$fit_path)) > 1L, by_key)
collision_keys <- sort(names(collisions))

cat(sprintf("\n  distinct (figure dir, stem) keys : %d\n", length(by_key)))
cat(sprintf("  keys with >1 underlying fit      : %d (allowlist %d)\n",
            length(collision_keys), length(KNOWN_COLLISIONS)))

# --- 1. no NEW collision -----------------------------------------------------
cat("\n[1] no new collision\n")

new_collisions <- setdiff(collision_keys, KNOWN_COLLISIONS)
if (length(new_collisions)) {
  cat("        new collisions (different fits sharing one figure path):\n")
  for (k in new_collisions) {
    d <- collisions[[k]]
    cat(sprintf("          - %s\n", k))
    for (i in seq_len(nrow(d)))
      cat(sprintf("              %s -> %s\n", d$config[i], d$fit_path[i]))
  }
}
expect(length(new_collisions) == 0L,
       "no config pair collides on (makeFigureDir, stem) outside KNOWN_COLLISIONS")

# --- 2. allowlist is not stale ------------------------------------------------
cat("\n[2] allowlist is not stale\n")

stale <- setdiff(KNOWN_COLLISIONS, collision_keys)
if (length(stale)) {
  cat("        repaired -- delete from KNOWN_COLLISIONS:\n")
  cat(paste0("          - ", stale, "\n"), sep = "")
}
expect(length(stale) == 0L, "KNOWN_COLLISIONS carries no already-repaired key")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
