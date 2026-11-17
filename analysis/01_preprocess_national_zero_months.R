# 01_preprocess_national_zero_months.R — list the country-months with no
#   suspected cholera anywhere in the country (national + WHO annual zeros).
#
# Usage:
#   Rscript analysis/01_preprocess_national_zero_months.R
#   Rscript analysis/01_preprocess_national_zero_months.R --end 2024-12 \
#     --exclude ETH:2011-2014
#   Rscript analysis/01_preprocess_national_zero_months.R --zero_threshold 60
#   Rscript analysis/01_preprocess_national_zero_months.R --zero_threshold 60 \
#     --rule r2
#
# Details:
#   The surveillance-coverage mask (utils/sim_io.R::build_coverage_mask())
#   counts these months as observed zeros for every gid of the country.
#   Inputs: <indir>/time_series_preoutbreak_extraction.rds, the outbreak
#   metadata, the WHO annual report + dashboard. Output:
#   make_national_zero_months_file() (+ _t<N>[_<rule>] variants, .policy.rds
#   sidecar, audit CSV, database-vs-WHO consistency CSV).
#
# See: analysis/utils/national_zero.R (the rule)

# 1. Preamble ------------------------------------------------------------------

library(dplyr)
library(optparse)

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/who_annual.R")
source("analysis/utils/national_zero.R")

# 2. Options -------------------------------------------------------------------
# Rule (utils/national_zero.R): ((national 0 AND no WHO > 0) OR (WHO 0 AND
# no national > 0)) AND no sub-national cases AND no active database outbreak
# AND not excluded (default ETH 2011-2014).
# --zero_threshold N: WHO country-years under N cases also count (see
# national_zero.R::db_who_consistency() for the exclusions). --rule: a named
# parameter set from analysis/national_zero_rules.yml (`legacy` = defaults);
# its parameters win, and a passed option that differs is an error.
# Outputs: make_national_zero_months_file() (+ _t<N>[_<rule>] variants), the
# .policy.rds sidecar, the audit CSV and the database-vs-WHO table read by
# report_who_scale.
NZ_LEGACY <- national_zero_rule_params("legacy")

option_list <- list(
  make_option(
    "--indir",
    default = default_outbreak_data_dir(), # utils/paths.R, the one data vintage
    type = "character",
    help = "Folder with time_series_preoutbreak_extraction.rds (as 01_preprocess_outbreaks.R) [default: %default]"
  ),
  make_option(
    "--who_annual_dir",
    default = "data/who_annual_report_africa",
    type = "character",
    help = "Folder with who_annual_report_africa.rds and who_dashboard_annual_africa.rds [default: %default]"
  ),
  make_option(
    "--rule",
    default = "",
    type = "character",
    help = "Named rule of analysis/national_zero_rules.yml ('' = legacy); its parameters replace the options below [default: legacy]"
  ),
  make_option(
    "--exclude",
    default = NZ_LEGACY$exclude,
    type = "character",
    help = "Excluded country-years, ISO3:Y0-Y1[,...]; '' for none [default: %default]"
  ),
  make_option(
    "--end",
    default = NZ_LEGACY$end,
    type = "character",
    help = "Last month (YYYY-MM); 2024-12 = the WHO annual report's last full year [default: %default]"
  ),
  make_option(
    "--zero_day_frac",
    default = NZ_LEGACY$zero_day_frac,
    type = "double",
    help = "Share of a month's days zero national weeks must cover [default: %default]"
  ),
  make_option(
    "--zero_threshold",
    default = 0L,
    type = "integer",
    help = "WHO country-years with fewer cases than this (max over sources) also count; 0 = exact zeros only; e.g. 60 (the 50-case outbreak floor plus a margin) [default: %default]"
  ),
  make_option(
    "--window",
    default = NZ_LEGACY$window,
    type = "integer",
    help = "Last months of a below-threshold year claimed only if the next year is also below (assumed episode length) [default: %default]"
  ),
  make_option(
    "--head_window",
    default = NZ_LEGACY$head_window,
    type = "integer",
    help = "First months of a below-threshold year claimed only if the previous year is also below (head guard; 0 = off) [default: %default]"
  ),
  make_option(
    "--head_missing_quiet",
    default = NZ_LEGACY$head_missing_quiet,
    type = "logical",
    help = "Head guard: a previous year with no WHO figure counts as below the threshold (rule r3) [default: %default]"
  ),
  make_option(
    "--consistency_scope",
    default = NZ_LEGACY$consistency_scope,
    type = "character",
    help = "threshold | all: the database-vs-WHO check also removes months resting on a WHO exact zero (all) [default: %default]"
  ),
  make_option(
    "--db_tolerance",
    default = NZ_LEGACY$db_tolerance,
    type = "double",
    help = "Threshold branch: drop country-years whose database episode cases exceed this x the WHO total [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

# Named rule: its parameters win; an EXPLICITLY passed option that disagrees is
# refused (two authorities never resolved by guessing).
opt$rule <- .national_zero_rule_check(opt$rule)
RULE_PARAMS <- NULL
if (nzchar(opt$rule)) {
  RULE_PARAMS <- national_zero_rule_params(opt$rule)
  argv <- commandArgs(trailingOnly = TRUE)
  for (k in setdiff(names(RULE_PARAMS), "rule")) {
    explicit <- any(
      argv == paste0("--", k) | startsWith(argv, paste0("--", k, "="))
    )
    if (
      explicit &&
        !identical(as.character(opt[[k]]), as.character(RULE_PARAMS[[k]]))
    ) {
      stop(
        "--",
        k,
        " = ",
        opt[[k]],
        " disagrees with rule '",
        opt$rule,
        "' (",
        RULE_PARAMS[[k]],
        " in ",
        NATIONAL_ZERO_RULES_PATH,
        "). Drop the option or the rule."
      )
    }
    opt[[k]] <- RULE_PARAMS[[k]]
  }
}
print(opt)

end_month <- as.Date(paste0(opt$end, "-01"))
if (is.na(end_month)) {
  stop("--end must be YYYY-MM, got: ", opt$end)
}
if (
  !is.finite(opt$zero_day_frac) ||
    opt$zero_day_frac <= 0 ||
    opt$zero_day_frac > 1
) {
  stop("--zero_day_frac must be in (0, 1], got: ", opt$zero_day_frac)
}
exclude <- parse_national_zero_exclude(opt$exclude)
if (is.na(opt$zero_threshold) || opt$zero_threshold < 0L) {
  stop("--zero_threshold must be >= 0, got: ", opt$zero_threshold)
}
if (is.na(opt$window) || opt$window < 0L || opt$window > 12L) {
  stop("--window must be in 0..12, got: ", opt$window)
}
if (!is.finite(opt$db_tolerance) || opt$db_tolerance <= 0) {
  stop("--db_tolerance must be > 0, got: ", opt$db_tolerance)
}
if (is.na(opt$head_window) || opt$head_window < 0L || opt$head_window > 12L) {
  stop("--head_window must be in 0..12, got: ", opt$head_window)
}
if (!opt$consistency_scope %in% NATIONAL_ZERO_CONSISTENCY_SCOPES) {
  stop(
    "--consistency_scope must be one of {",
    paste(NATIONAL_ZERO_CONSISTENCY_SCOPES, collapse = ", "),
    "}, got: ",
    opt$consistency_scope
  )
}

# 3. Data ----------------------------------------------------------------------

ts_path <- file.path(
  path.expand(opt$indir),
  "time_series_preoutbreak_extraction.rds"
)
ts <- readRDS(ts_path) %>%
  ungroup() %>%
  mutate(TL = as.Date(TL), TR = as.Date(TR))
ob_meta <- readRDS(make_outbreak_metadata_file())
who_ann <- file.path(opt$who_annual_dir, "who_annual_report_africa.rds")
who_dash <- file.path(opt$who_annual_dir, "who_dashboard_annual_africa.rds")
who_mask <- build_who_annual_zero_mask(who_ann, who_dash)
who_dash_tab <- readRDS(who_dash) # read once for the guard and the coverage message below
who_ann_tab <- readRDS(who_ann)
# A named rule refuses an end inside the dashboard's partial (latest) year.
if (nzchar(opt$rule)) {
  .national_zero_end_guard(
    end_month,
    max(as.integer(who_dash_tab$report_year), na.rm = TRUE)
  )
}

# 4. Rule ----------------------------------------------------------------------

nat <- national_month_status(
  filter(ts, spatial_scale == "country"),
  opt$zero_day_frac
)
subpos <- subnational_positive_months(filter(ts, spatial_scale != "country"))
active <- active_episode_months(ob_meta)
start <- lubridate::floor_date(
  min(ts$TL[ts$spatial_scale == "country"], na.rm = TRUE),
  "month"
)

# The database-vs-WHO table: always computed and persisted (the
# report_who_scale cross-check reads it); it shapes the months only on the threshold branch, or
# on every WHO-zero month under --consistency_scope all.
consistency_all <- db_who_consistency(
  ob_meta,
  who_mask,
  tolerance = opt$db_tolerance
)
consistency <- if (
  opt$zero_threshold > 0L || identical(opt$consistency_scope, "all")
) {
  consistency_all
} else {
  NULL
}
zero <- build_national_zero_months(
  nat,
  subpos,
  active,
  who_mask,
  start = start,
  end = end_month,
  exclude = exclude,
  zero_threshold = opt$zero_threshold,
  window = opt$window,
  auto_exclude = consistency,
  head_window = opt$head_window,
  consistency_scope = opt$consistency_scope,
  head_missing_quiet = isTRUE(opt$head_missing_quiet),
  edge_detail = TRUE
)
audit <- national_zero_audit(
  zero,
  who_mask,
  consistency_all,
  years = seq(
    as.integer(lubridate::year(start)),
    as.integer(lubridate::year(end_month))
  )
)
attr(zero, "edge_detail") <- NULL # the saved months carry the pre-registry attributes
dropped <- attr(zero, "dropped")
auto_excl <- if (is.null(consistency)) {
  consistency
} else {
  consistency[consistency$excluded, , drop = FALSE]
}

message(sprintf(
  "[01/national-zero] %s..%s: %d national-zero country-months in %d countries",
  format(start),
  format(end_month),
  nrow(zero),
  n_distinct(zero$country)
))
src_tab <- table(zero$source)
message(
  "[01/national-zero] by source: ",
  paste(sprintf("%s=%d", names(src_tab), as.integer(src_tab)), collapse = ", ")
)
message(
  "[01/national-zero] rule-zero months dropped: ",
  paste(sprintf("%s=%d", names(dropped), dropped), collapse = ", ")
)
if (opt$zero_threshold > 0L) {
  message(sprintf(
    paste0(
      "[01/national-zero] threshold %d (window %d): %d below-threshold month(s); ",
      "%d tail month(s) left unclaimed (next year not below)"
    ),
    opt$zero_threshold,
    opt$window,
    sum(zero$source == "who_below_threshold"),
    attr(zero, "tail_unclaimed")
  ))
  if (opt$head_window > 0L) {
    message(sprintf(
      "[01/national-zero] head guard (head_window %d%s): %d head month(s) left unclaimed (previous year not below)",
      opt$head_window,
      if (isTRUE(opt$head_missing_quiet)) {
        ", a missing previous WHO year counts as below"
      } else {
        ""
      },
      attr(zero, "head_unclaimed")
    ))
  }
  message(sprintf(
    "[01/national-zero] database vs WHO: %d country-year(s) compared, %d over %.2fx WHO (excluded from the threshold branch)",
    nrow(consistency),
    nrow(auto_excl),
    opt$db_tolerance
  ))
  if (nrow(auto_excl)) {
    print(
      as.data.frame(
        auto_excl %>% mutate(across(c(db_deep, who_max, ratio), ~ round(.x, 1)))
      ),
      row.names = FALSE
    )
  }
  below_tab <- zero %>%
    filter(source == "who_below_threshold") %>%
    count(country, year = as.integer(lubridate::year(month)), name = "months")
  if (nrow(below_tab)) print(as.data.frame(below_tab), row.names = FALSE)
}
yrs <- as.integer(lubridate::year(zero$month))
used <- who_mask %>%
  filter(year %in% unique(yrs)) %>%
  left_join(
    who_dash_tab %>%
      transmute(
        country = iso_3_code,
        year = as.integer(report_year),
        dash = TRUE
      ),
    by = c("country", "year")
  )
message(sprintf(
  paste0(
    "[01/national-zero] WHO country-years in the window: %d (%d with a ",
    "dashboard figure, %d dashboard-only)"
  ),
  nrow(used),
  sum(used$dash %in% TRUE),
  sum(
    used$dash %in%
      TRUE &
      !paste(used$country, used$year) %in%
        with(who_ann_tab, paste(iso3, format(TL, "%Y")))
  )
))

# 5. Write ---------------------------------------------------------------------

out_path <- make_national_zero_months_file(opt$zero_threshold, opt$rule)
saveRDS(zero, out_path)
AUDIT_PATH <- national_zero_audit_path(out_path)
# Rule-free stem (paths.R::national_zero_consistency_path(); read by
# report_who_scale): the table does not depend on the rule, only on db_tolerance.
CONSISTENCY_PATH <- national_zero_consistency_path(opt$zero_threshold)
if (!isTRUE(all.equal(opt$db_tolerance, NZ_LEGACY$db_tolerance))) {
  warning(
    "--db_tolerance ",
    opt$db_tolerance,
    " differs from the legacy ",
    NZ_LEGACY$db_tolerance,
    ": ",
    CONSISTENCY_PATH,
    " now holds this tolerance's flags."
  )
}
utils::write.csv(audit, AUDIT_PATH, row.names = FALSE)
utils::write.csv(consistency_all, CONSISTENCY_PATH, row.names = FALSE)
if (identical(opt$consistency_scope, "all")) {
  message(sprintf(
    "[01/national-zero] consistency scope all: %d WHO-zero month(s) dropped (database cases > %.2fx WHO)",
    dropped[["auto_excluded_exact"]],
    opt$db_tolerance
  ))
}
message(sprintf(
  "[01/national-zero] audit: %d country-year(s) -> %s",
  nrow(audit),
  AUDIT_PATH
))
persist_prep_policy_stamp(
  out_path,
  list(
    ts_path = ts_path,
    ts_mtime = format(file.mtime(ts_path)),
    ob_meta_path = as.character(make_outbreak_metadata_file()),
    ob_meta_mtime = format(file.mtime(make_outbreak_metadata_file())),
    who_annual_dir = opt$who_annual_dir,
    who_mask_digest = digest::digest(who_mask, algo = "xxhash64"),
    exclude = opt$exclude,
    start = format(start),
    end = format(end_month),
    zero_day_frac = opt$zero_day_frac,
    n_country_months = nrow(zero),
    n_by_source = as.list(stats::setNames(as.integer(src_tab), names(src_tab))),
    n_dropped = as.list(dropped),
    zero_threshold = opt$zero_threshold,
    window = if (opt$zero_threshold > 0L) opt$window else NULL,
    db_tolerance = if (opt$zero_threshold > 0L) opt$db_tolerance else NULL,
    n_tail_unclaimed = attr(zero, "tail_unclaimed"),
    auto_excluded = if (!is.null(auto_excl) && nrow(auto_excl)) {
      as.data.frame(auto_excl[, c(
        "country",
        "year",
        "db_deep",
        "who_max",
        "ratio"
      )])
    } else {
      NULL
    },
    # Named rule: written only under one, so a legacy stamp keeps
    # its field set; window / db_tolerance are then always recorded.
    rule = if (nzchar(opt$rule)) opt$rule else NULL,
    rule_window = if (nzchar(opt$rule)) opt$window else NULL,
    rule_db_tolerance = if (nzchar(opt$rule)) opt$db_tolerance else NULL,
    head_window = if (nzchar(opt$rule)) opt$head_window else NULL,
    consistency_scope = if (nzchar(opt$rule)) opt$consistency_scope else NULL,
    head_missing_quiet = if (
      nzchar(opt$rule) && isTRUE(opt$head_missing_quiet)
    ) {
      TRUE
    } else {
      NULL
    },
    n_head_unclaimed = if (nzchar(opt$rule)) {
      (attr(zero, "head_unclaimed") %||% 0L)
    } else {
      NULL
    },
    rule_params_digest = if (nzchar(opt$rule)) {
      national_zero_rule_digest(RULE_PARAMS)
    } else {
      NULL
    },
    audit_path = AUDIT_PATH,
    consistency_path = CONSISTENCY_PATH,
    stamped_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
)
message("[01/national-zero] wrote ", out_path)
