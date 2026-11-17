# test_config_block_registry.R — HARD: the shared top-level config block merge
#   (.CONFIG_TOPLEVEL_BLOCKS).
#
# Usage:
#   Rscript analysis/tests/test_config_block_registry.R
#
# Details:
#   make_options_from_config() merges every top-level block through the
#   registry. A block added to get_default_options() without a registry entry
#   would resolve to NULL for every config that omits it, silently, with no
#   error. Checks registry <-> defaults completeness in both directions,
#   well-formed entries, a round trip that materialises nested keys,
#   fill_when_absent = FALSE, and YAML values winning over defaults. Pure
#   config resolution: no fit, no data, no network.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

resolve_yaml <- function(yaml_body) {
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(yaml_body, tmp)
  make_options_from_config(option_list = list(), cli_argv = c("--config", tmp))
}

defaults <- get_default_options()

# Model blocks are merged by .fill_model_block(), not by the registry.
MODEL_BLOCKS <- c("occurrence", "size", "duration", "monthly", "climate")

cat("# config top-level block registry\n")

# --- 1. Registry <-> defaults completeness -----------------------------------
cat("\n[1] registry covers exactly the shared top-level blocks\n")

default_blocks <- setdiff(
  names(defaults)[vapply(defaults, is.list, logical(1))],
  MODEL_BLOCKS
)
registry_blocks <- names(.CONFIG_TOPLEVEL_BLOCKS)

expect(
  setequal(default_blocks, registry_blocks),
  sprintf(
    "registry == list-valued top-level defaults minus model blocks (defaults: %s | registry: %s)",
    paste(sort(default_blocks), collapse = ","),
    paste(sort(registry_blocks), collapse = ",")
  )
)
expect(
  length(setdiff(default_blocks, registry_blocks)) == 0L,
  "no default block is missing from the registry (would never be filled)"
)
expect(
  length(setdiff(registry_blocks, default_blocks)) == 0L,
  "no registry block is missing from get_default_options() (would merge against NULL)"
)
expect(
  length(intersect(registry_blocks, MODEL_BLOCKS)) == 0L,
  "no model block is in the registry (they use .fill_model_block)"
)

# --- 2. Registry entries are well-formed -------------------------------------
cat("\n[2] registry entries are well-formed\n")

expect(
  all(vapply(
    .CONFIG_TOPLEVEL_BLOCKS,
    function(s) identical(s$merge, "flat") || identical(s$merge, "recursive"),
    logical(1)
  )),
  "every entry declares merge = flat | recursive"
)
expect(
  all(vapply(
    .CONFIG_TOPLEVEL_BLOCKS,
    function(s) {
      is.logical(s$fill_when_absent) && length(s$fill_when_absent) == 1L
    },
    logical(1)
  )),
  "every entry declares a scalar logical fill_when_absent"
)

# --- 3. Round-trip: population block is materialised, nested keys included ---
cat("\n[3] minimal config materialises the population block\n")

opt <- resolve_yaml("occurrence:\n  backend: brms\n")

expect(
  !is.null(opt$population),
  "opt$population present even though the YAML never mentions it"
)
expect(
  identical(
    opt$population$sliver_mass_frac,
    defaults$population$sliver_mass_frac
  ),
  "population.sliver_mass_frac filled (flat key)"
)
expect(
  identical(
    opt$population$slurm$array_concurrency,
    defaults$population$slurm$array_concurrency
  ),
  "population.slurm.array_concurrency filled (nested — .merge_recursive ran)"
)
expect(
  identical(
    opt$population$composite$cover_threshold,
    defaults$population$composite$cover_threshold
  ),
  "population.composite.cover_threshold filled (nested)"
)

# --- 4. fill_when_absent = FALSE is honoured ---------------------------------
cat("\n[4] fill_when_absent = FALSE leaves the block NULL\n")

expect(
  identical(.CONFIG_TOPLEVEL_BLOCKS$simulator$fill_when_absent, FALSE),
  "simulator is registered with fill_when_absent = FALSE"
)
expect(
  is.null(opt$simulator),
  "opt$simulator stays NULL for a non-simulator config (downstream tests is.null)"
)
expect(
  !is.null(opt$enso) && !is.null(opt$grouping),
  "enso/grouping are still always materialised (fill_when_absent = TRUE)"
)

opt_sim <- resolve_yaml(paste(
  "seed: 999",
  "simulator:",
  "  n_replicates: 3",
  sep = "\n"
))
expect(
  !is.null(opt_sim$simulator),
  "opt$simulator present when the YAML declares it"
)
expect(
  identical(opt_sim$simulator$n_replicates, 3L) ||
    identical(opt_sim$simulator$n_replicates, 3),
  "YAML value survives the merge"
)
expect(
  identical(opt_sim$simulator$seed, 999L) ||
    identical(opt_sim$simulator$seed, 999),
  "simulator.seed still inherits the top-level seed after the loop"
)
expect(
  !is.null(opt_sim$simulator$horizon),
  "simulator sibling sub-blocks still filled from defaults"
)

# --- 5. YAML wins over defaults; siblings still filled -----------------------
cat("\n[5] partial population block: YAML wins, siblings default\n")

opt_pop <- resolve_yaml(paste(
  "population:",
  "  reconcile_tol: 0.25",
  "  slurm:",
  "    array_concurrency: 4",
  sep = "\n"
))
expect(
  identical(opt_pop$population$reconcile_tol, 0.25),
  "explicit population.reconcile_tol overrides the default"
)
expect(
  identical(opt_pop$population$slurm$array_concurrency, 4L) ||
    identical(opt_pop$population$slurm$array_concurrency, 4),
  "explicit nested population.slurm.array_concurrency overrides the default"
)
expect(
  identical(
    opt_pop$population$slurm$partition,
    defaults$population$slurm$partition
  ),
  "unspecified nested sibling population.slurm.partition still defaults"
)
expect(
  identical(
    opt_pop$population$sliver_mass_min,
    defaults$population$sliver_mass_min
  ),
  "unspecified flat sibling population.sliver_mass_min still defaults"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
