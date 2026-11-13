# test_config_split_symbols.R -- which module owns each former config.R symbol.
#
# HARD gate (see analysis/tests/run_all.R). utils/config.R's symbols live in
# config.R (options core) + formula.R + formula_legacy.R + priors.R. This test pins the
# ownership table:
# it fails when a symbol moves between the four modules, disappears, appears in
# two of them, or a new top-level definition is added without being listed here.
# A deliberate move or addition is a one-line edit to OWNERS in the same commit.
#
# Pure parse -- no fit, no data, no network.
# Run: Rscript analysis/tests/test_config_split_symbols.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# config.R split: symbol ownership across config / formula / formula_legacy / priors\n")

OWNERS <- list(
  config = c(
    ".flag_val", "simulator_common_options", "apply_simulator_common_cli",
    "require_simulator_block", ".validate_basehaz", ".merge1", ".merge_recursive",
    "config_get_leaf", "config_set_leaf", "config_leaf_paths", "config_flatten_leaves",
    ".CONFIG_TOPLEVEL_BLOCKS", ".fill_model_block", ".outbreak_data_root", "get_default_options",
    "load_grouping", ".read_registry_csv", "load_adjacency", "resolve_group_countries", "make_options_from_config",
    # grouping-registry helpers
    "registry_path", "registry_digest", "grouping_label_token", ".GROUP_FALLBACK_COLOURS",
    "load_grouping_meta", "resolve_group_label", "grouping_source_from_token",
    "assert_partition_labels",
    ".MCMC_ROUTE_KEYS", ".route_fragment", ".CLI_ROUTE",
    ".OCC_CLI_ROUTE", ".SIZE_CLI_ROUTE", ".DUR_CLI_ROUTE", ".MON_CLI_ROUTE",
    ".route_cli_to_nested", "resolve_brms_mcmc", "resolve_brms_init", "resolve_brms_threads", "print_options", "detect_model_from_yaml"),
  formula = c(
    ".occ_basehaz", ".occ_baseline_block",
    "integrate_cb_terms", "re_group_by_level", "nested_slope_re_terms",
    "filter_re_terms_for_data", "filter_fixed_factor_terms_for_data",
    ".single_country_keep", ".flag_single_country",
    ".resolve_fixed_covariates", ".slope_re_terms", ".re_terms_for_column",
    ".intercept_re_term", ".scoped_intercept_re_terms", ".leaf_scope_of", ".leafColumns",
    ".OFFSET_SPEC", ".offset_terms", ".offset_spec", ".anchor_gate_on",
    "build_formula_terms", "validate_size_options", "resolve_slopes",
    ".lvl_set", ".seasonal_terms", ".phi_re_term",
    ".pop_elasticity_terms", ".build_count_track_parts", ".monthly_cyclic_slope_check",
    ".build_nb_formula",
    "build_size_formula_parts", "build_size_formula",
    ".DUR_MS_PREFIX", ".DUR_SEAS_MS_PREFIX", ".dur_ms_cols",
    "build_duration_formula_parts", "build_duration_formula",
    "build_monthly_formula_parts", "build_monthly_formula"),
  formula_legacy = c(".occ_group_by_level", "build_formula_terms_legacy"),
  priors = c(
    ".PRIOR_DEFAULTS", ".prior_strings", ".intercept_prior_string",
    ".prior_formula_gates", ".re_term_group_names", ".re_group_tiers",
    ".tiered_sd_priors", ".has_smooth_terms", ".prior_config_slopes",
    ".shape_dpar_priors", ".build_hazard_track_priors", ".build_count_track_priors",
    ".SINGLE_COUNTRY_SD", ".single_country_strings",
    "build_occ_priors", "build_duration_priors",
    "build_size_priors", "build_monthly_priors"))

defined_in <- function(m) {
  ex <- parse(file.path("analysis/utils", paste0(m, ".R")), keep.source = FALSE)
  vapply(ex, function(e) as.character(e[[2]]), "")
}

all_listed <- unlist(OWNERS, use.names = FALSE)
expect(!anyDuplicated(all_listed), "OWNERS lists each symbol once")

actual <- lapply(setNames(names(OWNERS), names(OWNERS)), defined_in)
for (m in names(OWNERS)) {
  extra   <- setdiff(actual[[m]], OWNERS[[m]])
  missing <- setdiff(OWNERS[[m]], actual[[m]])
  expect(length(extra) == 0L && length(missing) == 0L,
         sprintf("%s.R defines exactly its listed symbols (unlisted: %s; missing: %s)", m,
                 if (length(extra)) paste(extra, collapse = ", ") else "none",
                 if (length(missing)) paste(missing, collapse = ", ") else "none"))
  expect(!anyDuplicated(actual[[m]]), sprintf("%s.R defines no symbol twice", m))
}
dup <- unique(unlist(actual)[duplicated(unlist(actual))])
expect(length(dup) == 0L, sprintf("no symbol defined in two of the four modules (dupes: %s)",
                                  if (length(dup)) paste(dup, collapse = ", ") else "none"))
expect(!"%||%" %in% unlist(actual), "`%||%` lives in base.R only, not in the config modules")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
