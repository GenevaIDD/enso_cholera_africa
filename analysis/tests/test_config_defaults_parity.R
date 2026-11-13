# test_config_defaults_parity.R — config_defaults.yml vs get_default_options():
# same key paths AND the same values.
#
# HARD gate (see analysis/tests/run_all.R).
#
# WHY THIS EXISTS
# ---------------
# analysis/config_defaults.yml and get_default_options() (analysis/utils/config.R)
# express the same defaults twice, by hand. Only the R side fills missing keys
# at run time; the YAML's one code reader is utils/make_batch_configs.R, which
# copies it into every generated fit-batch config. So the YAML is the
# documentation of the parameter surface plus the batch generator's base, and
# get_default_options() is the effective default; the runtime R defaults win.
# The two must agree on every key path and every value.
#
# Values are compared after normalisation: integer and double are one numeric
# type (YAML reads `5` as an integer, R writes `5`), and a YAML sequence equals
# the same atomic vector (`[1, 2, Inf]` vs `c(1, 2, Inf)`). NULL (YAML `~`) and
# an empty list/vector are distinct: under the recursive merge a NULL default
# fills nothing, an empty one fills an empty value.
#
# Allowlists below are the only way a difference can pass, and the test fails
# on a stale entry too, so they cannot silently rot. All three are empty.
#
# THE `deploy` SUBTREE IS STRUCTURALLY OUT OF SCOPE (see drop_deploy() below).
# get_default_options()$deploy exists only to back the deploy orchestrator; it
# is DELIBERATELY absent from config_defaults.yml (each deploy run is
# self-contained in analysis/configs/deploy/<stem>.yml, never filled from the
# defaults file). Its parity is owned by the stronger, bidirectional
# test_deploy_config_coverage.R gate, so this test excises the subtree.
#
# Run: Rscript analysis/tests/test_config_defaults_parity.R

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
# The defaults are pinned with the outbreak-data root at its default
# (CHOLERA_OUTBREAK_DATA_DIR unset; test_data_vintage_default.R [1b] covers the override).
Sys.unsetenv("CHOLERA_OUTBREAK_DATA_DIR")

# Keys present in config_defaults.yml with NO default in get_default_options().
# A config omitting any of these would get NULL, not the documented value.
KNOWN_YAML_ONLY <- character(0)

# Keys with an R default but no YAML documentation.
KNOWN_R_ONLY <- character(0)

# Leaves present on both sides whose values may differ, each with the reason:
#   list(`<dotted.path>` = "<why the difference is intentional>")
KNOWN_VALUE_DIFF <- list()

#' Recursive key-path enumeration of a nested named list.
#' Unnamed lists (YAML sequences, e.g. population.years) are leaves.
key_paths <- function(l, prefix = "") {
  out <- character(0)
  for (nm in names(l)) {
    p <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    out <- c(out, p)
    if (is.list(l[[nm]]) && !is.null(names(l[[nm]]))) {
      out <- c(out, key_paths(l[[nm]], p))
    }
  }
  out
}

#' Leaf values of a nested named list, as a named list keyed by dotted path.
#' A leaf is anything that is not a non-empty named list (so NULL leaves, empty
#' lists and YAML sequences are leaves).
leaf_values <- function(l, prefix = "") {
  out <- list()
  for (nm in names(l)) {
    p <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    v <- l[[nm]]
    if (is.list(v) && !is.null(names(v)) && length(v)) {
      out <- c(out, leaf_values(v, p))
    } else {
      out[p] <- list(v)
    }
  }
  out
}

#' Normal form for comparing a YAML leaf with an R leaf (see header).
norm_value <- function(v) {
  if (is.null(v)) return(list(kind = "NULL"))
  u <- unlist(v, use.names = FALSE)
  if (!length(u)) return(list(kind = "empty"))
  if (is.numeric(u)) u <- as.numeric(u)
  list(kind = "value", value = u)
}

fmt_value <- function(v) paste(deparse(v, width.cutoff = 500L), collapse = "")

yml <- yaml::read_yaml("analysis/config_defaults.yml")
dfl <- get_default_options()

#' Drop the `deploy` subtree: its parity is gated by test_deploy_config_coverage.R
#' (config_defaults.yml deliberately carries no deploy block — see header).
drop_deploy <- function(paths) paths[!grepl("^deploy(\\.|$)", paths)]

yaml_only <- sort(setdiff(drop_deploy(key_paths(yml)), drop_deploy(key_paths(dfl))))
r_only    <- sort(setdiff(drop_deploy(key_paths(dfl)), drop_deploy(key_paths(yml))))

cat("# config defaults parity (key paths and values)\n")
cat(sprintf("\n  YAML paths with no R default : %d (allowlist %d)\n",
            length(yaml_only), length(KNOWN_YAML_ONLY)))
cat(sprintf("  R defaults with no YAML doc  : %d (allowlist %d)\n",
            length(r_only), length(KNOWN_R_ONLY)))

# --- 1. No key-path drift ----------------------------------------------------
cat("\n[1] no key-path drift\n")

new_yaml_only <- setdiff(yaml_only, KNOWN_YAML_ONLY)
new_r_only    <- setdiff(r_only,    KNOWN_R_ONLY)

if (length(new_yaml_only)) {
  cat("        YAML-only paths (documented but never defaulted):\n")
  cat(paste0("          - ", new_yaml_only, "\n"), sep = "")
}
if (length(new_r_only)) {
  cat("        R-only paths (defaulted but undocumented):\n")
  cat(paste0("          - ", new_r_only, "\n"), sep = "")
}

expect(length(new_yaml_only) == 0L,
       "no key in config_defaults.yml without a get_default_options() entry")
expect(length(new_r_only) == 0L,
       "no key in get_default_options() without a config_defaults.yml entry")

# --- 2. Values agree -----------------------------------------------------------
cat("\n[2] values agree on every common leaf\n")

ly <- leaf_values(yml); ly <- ly[drop_deploy(names(ly))]
lr <- leaf_values(dfl); lr <- lr[drop_deploy(names(lr))]
common <- intersect(names(ly), names(lr))
value_diff <- character(0)
for (k in common) {
  if (!identical(norm_value(ly[[k]]), norm_value(lr[[k]]))) value_diff <- c(value_diff, k)
}
cat(sprintf("        %d common leaves compared, %d differ (allowlist %d)\n",
            length(common), length(value_diff), length(KNOWN_VALUE_DIFF)))

new_value_diff <- setdiff(value_diff, names(KNOWN_VALUE_DIFF))
if (length(new_value_diff)) {
  cat("        differing values (YAML | get_default_options()):\n")
  for (k in new_value_diff)
    cat(sprintf("          - %s: %s | %s\n", k, fmt_value(ly[[k]]), fmt_value(lr[[k]])))
}
expect(length(new_value_diff) == 0L,
       "every config_defaults.yml value equals its get_default_options() value")
expect(length(common) > 500L,
       sprintf("the comparison covers the whole surface (%d leaves)", length(common)))

# --- 3. Allowlists have not gone stale --------------------------------------
cat("\n[3] allowlists are not stale\n")

stale_yaml  <- setdiff(KNOWN_YAML_ONLY, yaml_only)
stale_r     <- setdiff(KNOWN_R_ONLY,    r_only)
stale_value <- setdiff(names(KNOWN_VALUE_DIFF), value_diff)

for (s in list(list("KNOWN_YAML_ONLY", stale_yaml), list("KNOWN_R_ONLY", stale_r),
               list("KNOWN_VALUE_DIFF", stale_value))) {
  if (length(s[[2]])) {
    cat(sprintf("        repaired — delete from %s:\n", s[[1]]))
    cat(paste0("          - ", s[[2]], "\n"), sep = "")
  }
}
expect(length(stale_yaml) == 0L,  "KNOWN_YAML_ONLY carries no already-repaired path")
expect(length(stale_r) == 0L,     "KNOWN_R_ONLY carries no already-repaired path")
expect(length(stale_value) == 0L, "KNOWN_VALUE_DIFF carries no already-equal leaf")

# --- 4. The population and simulator.filter blocks are fully defaulted --------
# The blocks this gate protects most directly (population; simulator.filter).
# Explicit rather than implied by [1], so a regression names itself.
cat("\n[4] population and simulator.filter blocks are fully defaulted (both sides)\n")

for (blk in c("population", "simulator.filter")) {
  rx <- paste0("^", gsub(".", "\\.", blk, fixed = TRUE), "(\\.|$)")
  expect(length(grep(rx, yaml_only)) == 0L,
         sprintf("no %s path is documented-but-undefaulted", blk))
  expect(length(grep(rx, r_only)) == 0L,
         sprintf("no %s path is defaulted-but-undocumented", blk))
}

# --- 5. The normalisation is not hiding a real difference ----------------------
cat("\n[5] value normalisation\n")
expect(identical(norm_value(5L), norm_value(5)), "integer 5 == double 5")
expect(identical(norm_value(list(1L, 2L, Inf)), norm_value(c(1, 2, Inf))),
       "a YAML sequence == the same atomic vector")
expect(!identical(norm_value("rds"), norm_value("qs")), "\"rds\" != \"qs\"")
expect(!identical(norm_value(NULL), norm_value(list())), "NULL != empty list")
expect(!identical(norm_value(TRUE), norm_value(1)), "TRUE != 1")
expect(!identical(norm_value(c(1, 2)), norm_value(c(2, 1))), "order matters")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
