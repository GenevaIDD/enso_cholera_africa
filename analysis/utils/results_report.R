# analysis/utils/results_report.R
#
# Display helpers for the phase-11 key-numbers report
# (analysis/11_manuscript_numbers.qmd): number formatters, label/level
# registries and kable table emission. Moved out of the notebook's init chunk so
# the .qmd is display code only and the formatting rules are one testable copy.
#
# Sourced-only; sources nothing (README §5). Call-time dependencies: the
# registries GROUP_SHORT / GROUP_ORDER and afro_subregion_of() (utils/results_core.R) and
# enso_conds_named() (utils/figures_cb.R, sourced by results_core.R); knitr,
# kableExtra and tidyr are used via `::`. Never computes a quantile: the phase-09
# artefacts carry every summary, this module only formats them.
#
# Sourced by: analysis/11_manuscript_numbers.qmd

# ---- formatters: "reasonable digits", enforced centrally --------------------
# Every formatter returns "—" for NA.

fmt_num   <- function(x, digits = 2) ifelse(is.na(x), "—",
                 formatC(as.numeric(x), format = "f", digits = digits, big.mark = ","))
fmt_int   <- function(x) ifelse(is.na(x), "—",
                 formatC(round(as.numeric(x)), format = "d", big.mark = ","))
# Case / onset counts: nearest 10 below 1,000, nearest 100 at/above 1,000.
fmt_cases <- function(x) {
  v <- as.numeric(x)
  r <- ifelse(abs(v) < 1000, round(v / 10) * 10, round(v / 100) * 100)
  ifelse(is.na(v), "—", formatC(r, format = "d", big.mark = ","))
}
# Populations: millions, 0 digits by default (descriptive tables pass digits = 2).
fmt_pop   <- function(x, digits = 0) ifelse(is.na(x), "—",
                 paste0(formatC(as.numeric(x) / 1e6, format = "f", digits = digits), "M"))
fmt_pop_ci <- function(m, lo, hi, digits = 0) ifelse(is.na(m), "—",
                 sprintf("%s (%s, %s)", fmt_pop(m, digits), fmt_pop(lo, digits), fmt_pop(hi, digits)))
fmt_pct   <- function(x, digits = 0) ifelse(is.na(x), "—",
                 paste0(formatC(100 * as.numeric(x), format = "f", digits = digits), "%"))
# Population share (a fraction) with its 95% interval, as percentages.
fmt_pct_ci <- function(m, lo, hi, digits = 0) ifelse(is.na(m), "\u2014",
                 sprintf("%s (%s, %s)", fmt_pct(m, digits), fmt_pct(lo, digits), fmt_pct(hi, digits)))
fmt_ci    <- function(m, lo, hi, digits = 2)
  sprintf("%s (%s, %s)", fmt_num(m, digits), fmt_num(lo, digits), fmt_num(hi, digits))
# Ratios (IRR, RMD ratio, cumulative/peak ratio): 1 digit.
fmt_ratio <- function(m, lo, hi, digits = 1) fmt_ci(m, lo, hi, digits)
# Bold (raw HTML) when the 95% interval excludes 1.
.bold_if_excludes_one <- function(txt, lo, hi) {
  sig <- !is.na(lo) & !is.na(hi) & (lo > 1 | hi < 1)
  ifelse(sig, paste0("<strong>", txt, "</strong>"), txt)
}
fmt_ratio_bold <- function(m, lo, hi, digits = 1)
  .bold_if_excludes_one(fmt_ratio(m, lo, hi, digits), lo, hi)
# Median + 95% of a case count (magnitude-rounded) / an onset count (exact).
fmt_ci_cases <- function(m, lo, hi) sprintf("%s (%s, %s)", fmt_cases(m), fmt_cases(lo), fmt_cases(hi))
fmt_ci_int   <- function(m, lo, hi) sprintf("%s (%s, %s)", fmt_int(m), fmt_int(lo), fmt_int(hi))
# Incidence per 100,000 person-years: 1 digit.
fmt_rate     <- function(m, lo, hi, digits = 1) fmt_ci(m, lo, hi, digits)
# Ratio with P(increase) bracketed in the same cell; bold when the CI excludes 1.
fmt_ratio_prob <- function(m, lo, hi, prob, digits = 1)
  .bold_if_excludes_one(sprintf("%s (%s, %s) [%s]", fmt_num(m, digits), fmt_num(lo, digits),
                                fmt_num(hi, digits), fmt_pct(prob)), lo, hi)
# Probability ratio (report_attrib_metrics) with its 95% interval: Inf as "∞"
# (the counterfactual never has the event), bold when the interval excludes 1.
.fmt_num_inf <- function(x, digits) ifelse(is.infinite(x) & x > 0, "\u221e", fmt_num(x, digits))
fmt_ratio_inf <- function(m, lo, hi, digits = 1)
  .bold_if_excludes_one(ifelse(is.na(m), "\u2014", sprintf("%s (%s, %s)", .fmt_num_inf(m, digits),
                                .fmt_num_inf(lo, digits), .fmt_num_inf(hi, digits))), lo, hi)
fmt_iqr   <- function(med, q25, q75, f = fmt_num, ...)
  sprintf("%s [%s–%s]", f(med, ...), f(q25, ...), f(q75, ...))

#' Median + 95% cell for a simulation quantity, formatted by its kind: onset
#' counts exact, incidence rates 1 digit, everything else as a case count.
#' @param quantity e.g. "cases", "onsets", "tot_cases", "mean_onsets", "inc_ann_100k"
fmt_by_quantity <- function(quantity, med, lo, hi)
  dplyr::case_when(grepl("onsets$", quantity) ~ fmt_ci_int(med, lo, hi),
                   grepl("^inc_", quantity)   ~ fmt_rate(med, lo, hi),
                   TRUE                       ~ fmt_ci_cases(med, lo, hi))

# ---- labels and level registries -------------------------------------------

#' Display-only shortening of the ENSO condition labels ("La Niña (−2.0)" ->
#' "La Niña (−2)"), so wide-table headers fit on one line.
short_enso <- function(x) sub("\\.0(?=\\))", "", as.character(x), perl = TRUE)

#' Scenario label (the forecast arm carries analog_set = NA).
scenario_lab <- function(run, analog_set)
  ifelse(run == "forecast", "forecast", paste0("counterfactual: ", analog_set))

#' Display orderings, all read from the results config / registries so the
#' report follows the configured ONI conditions and analog sets.
#' @param res the resolved `results:` block
#' @param group_meta load_grouping_meta() list of the fit registry; NULL = the
#'   registry named by the forecast manifest (res$inputs$manifest), else enso_grouping
#' @return list(region, cond, scope, scenario(sets))
report_levels <- function(res, group_meta = NULL) {
  if (is.null(group_meta)) {
    mp <- as.character(res$inputs$manifest %||% "")
    group_meta <- results_group_meta(
      if (length(mp) == 1L && nzchar(mp) && file.exists(mp)) list(MAN = yaml::read_yaml(mp)) else NULL)
  }
  region <- unname(group_meta$short[group_meta$order])
  list(region = region,
       cond   = short_enso(names(enso_conds_named(res$derivation$enso_conds))),
       scope  = c("Africa", region),
       scenario = function(sets) c("forecast", paste0("counterfactual: ", sets)))
}

QUANTITY_LABELS <- c(tot_cases = "Total cases", tot_onsets = "Total onsets",
                     mean_cases = "Mean monthly cases", mean_onsets = "Mean monthly onsets",
                     inc_ann_100k = "Incidence /100k/yr")

# ---- table shaping + emission ------------------------------------------------

#' Pivot a long table so `axis` becomes columns: factor `rows`/`axis` columns
#' with the given levels, sort by `rows`, keep rows + axis + value, pivot wide
#' (columns in level order), then rename (new = old).
#' @param levels named list column -> levels (columns not listed keep their order)
report_wide <- function(df, rows, axis, value, levels = list(), rename = character(0)) {
  for (nm in names(levels)) df[[nm]] <- factor(df[[nm]], levels = levels[[nm]])
  out <- df %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(rows))) %>%
    dplyr::select(dplyr::all_of(c(rows, axis, value))) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(axis), values_from = dplyr::all_of(value),
                       names_sort = TRUE)
  if (length(rename)) out <- dplyr::rename(out, dplyr::all_of(rename))
  out
}

#' Append an Africa-total frame to a per-region frame as a trailing "Africa"
#' group (region factor levels = regions then Africa). Extra named arguments set
#' columns on the Africa rows (e.g. country = "All Africa").
report_with_africa <- function(base, africa, region_levels, ...) {
  ex <- list(...)
  af <- africa %>% dplyr::mutate(region = "Africa")
  for (nm in names(ex)) af[[nm]] <- ex[[nm]]
  dplyr::bind_rows(base %>% dplyr::mutate(region = as.character(region)), af) %>%
    dplyr::mutate(region = factor(region, levels = c(region_levels, "Africa")))
}

#' Running table counter: every emitted table gets "Table N." in document order.
report_table_counter <- function() local({ n <- 0L; function() { n <<- n + 1L; n } })

#' File stem of a table's CSV: an explicit `slug` (e.g. "D27_convergence") or
#' the first words of its caption, prefixed by the table number.
report_csv_name <- function(n, cap, slug = NULL) {
  s <- slug %||% {
    w <- strsplit(tolower(gsub("[^A-Za-z0-9]+", " ", cap)), " +")[[1]]
    paste(utils::head(w[nzchar(w)], 6), collapse = "_")
  }
  sprintf("%02d_%s.csv", n, s)
}

#' Write one emitted table as CSV (machine-readable copy for the manuscript):
#' HTML markup and the "—" placeholder stripped from character cells.
report_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in names(df)) if (is.factor(df[[nm]]) || is.character(df[[nm]])) {
    v <- gsub("<[^>]+>", "", as.character(df[[nm]]))
    v[v == "\u2014"] <- NA_character_
    df[[nm]] <- v
  }
  utils::write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}

#' Emit a kable-styled HTML table (works inside results='asis' loops).
#' `pack_by` names a column (already sorted so groups are contiguous) whose
#' consecutive repeats become bold group-header rows; it is dropped from the body.
#' With `csv_dir`, the table is also written as `<NN>_<slug>.csv` there (the
#' grouping column kept); `raw` replaces the displayed frame in the CSV (e.g.
#' unformatted numbers).
report_emit_html <- function(df, cap, pack_by = NULL, counter, slug = NULL, csv_dir = NULL,
                             raw = NULL) {
  n   <- counter()
  if (!is.null(csv_dir) && nzchar(csv_dir) && !is.null(df) && nrow(df))
    report_write_csv(raw %||% df, file.path(csv_dir, report_csv_name(n, cap, slug)))
  cap <- sprintf("Table %d. %s", n, cap)
  if (is.null(df) || !nrow(df)) { cat("*", cap, ": no rows.*\n\n"); return(invisible()) }
  grp <- NULL
  if (!is.null(pack_by)) {
    grp <- as.character(df[[pack_by]])
    df  <- df[, setdiff(names(df), pack_by), drop = FALSE]
  }
  k <- knitr::kable(df, caption = cap, format = "html", align = "l", escape = FALSE) |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                              full_width = FALSE, position = "left")
  if (!is.null(grp)) {
    r <- rle(grp)
    starts <- cumsum(c(1L, utils::head(r$lengths, -1)))
    for (i in seq_along(r$values))
      k <- kableExtra::pack_rows(k, r$values[i], starts[i], starts[i] + r$lengths[i] - 1L)
  }
  cat(as.character(k), "\n\n")
}

#' Artefact loader for the report: R_(kind, ...) resolves any results_<kind>
#' path at its RESULTS_SCHEMA token and reads it through results_require()
#' (stops with the exact 09 re-run command when missing, stale or derived under
#' other decisions).
report_loader <- function(opt, config_path, res) function(kind, ...)
  results_require(makeResultsRunFile(opt, kind, ...,
                                     schema = RESULTS_SCHEMA[[sub("^results_", "", kind)]],
                                     config_path = config_path),
                  config_path, res)

# ---- country rows: WHO-AFRO subregion ---------------------------------------

# WHO-AFRO subregion beside the country (supplement convention; fixed map
# results_core.R::AFRO_SUBREGION, the pipeline's AFRO_region convention).
# Printed tables get an "AFRO subregion" column right after `col`
# ("‡" on the WHO EMRO members SDN, SOM, DJI; emit() then appends AFRO_NOTE to
# the caption); raw = TRUE (the CSV frames) adds afro_subregion and who_region
# instead. `col` NULL = the first of iso3 / country / Country / ISO3 present.
# Rows that are not a country (Africa, a group) get an empty cell.
AFRO_COL  <- "AFRO subregion"
AFRO_NOTE <- paste("‡ WHO EMRO member state, listed under Eastern Africa; subregions follow the",
                   "pipeline's AFRO_region convention (taxdat), not an official WHO subdivision.")
with_afro <- function(df, col = NULL, raw = FALSE) {
  if (is.null(df) || !nrow(df)) return(df)
  col <- col %||% intersect(c("iso3", "country", "Country", "ISO3"), names(df))[1]
  if (is.na(col) || !col %in% names(df) || any(c(AFRO_COL, "afro_subregion") %in% names(df))) return(df)
  a <- afro_subregion_of(df[[col]])
  new <- if (isTRUE(raw)) a else
    stats::setNames(data.frame(ifelse(is.na(a$afro_subregion), "",
                                      paste0(a$afro_subregion, ifelse(a$who_region %in% "EMRO", "‡", ""))),
                               stringsAsFactors = FALSE), AFRO_COL)
  at <- match(col, names(df))
  dplyr::bind_cols(df[seq_len(at)], new, df[-seq_len(at)])
}

# ---- data release split (results.style.data_release) -------------------------

#' Data id of a table CSV: "D12" from "57_D12_enso_exposure_by_country.csv" or
#' "D30_who_scale.csv", "T3" from "T3a_effects_occurrence.csv"; NA for a CSV
#' without a D / T slug (main-text tables).
report_data_id <- function(file) {
  b <- basename(file)
  ifelse(grepl("^([0-9]+_)?[DT][0-9]+[a-z]?_", b),
         sub("^([0-9]+_)?([DT][0-9]+)[a-z]?_.*$", "\\2", b), NA_character_)
}

#' Copy one table CSV into a release folder without its all-empty columns
#' (e.g. the bootstrap intervals when results.thresholds.attrib_boot_b is 0).
#' Cells are copied as text, so no value is re-formatted; a column is left
#' unquoted only when every value prints back as itself as a number or a
#' logical (what write.csv leaves unquoted), so "+6" or "0.10" stay quoted.
#' @noRd
.release_copy <- function(src, dst) {
  d <- utils::read.csv(src, colClasses = "character", na.strings = character(0), check.names = FALSE)
  keep <- vapply(d, function(x) any(nzchar(x)), logical(1))
  if (!ncol(d) || !nrow(d) || all(keep)) return(invisible(file.copy(src, dst, overwrite = TRUE)))
  d <- d[, keep, drop = FALSE]
  txt <- vapply(d, function(x) {
    v <- x[nzchar(x)]
    if (!length(v)) return(FALSE)
    if (all(v %in% c("TRUE", "FALSE"))) return(FALSE)
    n <- suppressWarnings(as.numeric(v))
    anyNA(n) || !identical(as.character(n), v)
  }, logical(1))
  utils::write.table(d, dst, sep = ",", quote = which(txt), row.names = FALSE, na = "", qmethod = "double")
  invisible(TRUE)
}

#' report_release_split
#'
#' Copy every D / T table CSV of `tables_dir` into the release folder its data
#' id maps to under results.style.data_release (supp -> `dirs[["supp"]]`,
#' github -> `dirs[["github"]]`), under its slug without the table-number
#' prefix, and write the manifest `data_release_manifest.csv` (file, source,
#' data_id, release, caption). Columns empty in every row are left out of the
#' released copy (.release_copy()); the numbered table keeps them. A D / T CSV whose id the map does not name is
#' listed with an empty release and not copied. A CSV whose slug is in
#' `exclude` (results.style.data_release_exclude) is listed with release
#' "withheld" and not copied. The release folders are cleared first.
#' @param release named list / character (data id -> "supp" | "github")
#' @param captions named character (slug -> caption), from emit()
#' @param exclude character, slugs withheld from the release
#' @return the manifest data frame (invisibly)
report_release_split <- function(tables_dir, release, captions = character(0),
                                 dirs = c(supp = "supp_data", github = "github"),
                                 exclude = character(0)) {
  release <- unlist(release %||% list())
  bad <- setdiff(unique(release), names(dirs))
  if (length(bad))
    stop("report_release_split: results.style.data_release values must be one of ",
         paste(names(dirs), collapse = " | "), "; got ", paste(bad, collapse = ", "), call. = FALSE)
  for (d in dirs) {
    p <- file.path(tables_dir, d)
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(p, "\\.csv$", full.names = TRUE))
  }
  f <- sort(list.files(tables_dir, "\\.csv$"))
  f <- f[!is.na(report_data_id(f))]
  id <- report_data_id(f)
  slug <- sub("\\.csv$", "", sub("^[0-9]+_", "", f))
  out_name <- paste0(make.unique(slug, sep = "_"), ".csv")
  rel <- if (length(f)) unname(release[id]) else character(0)
  rel[is.na(rel)] <- ""
  rel[nzchar(rel) & slug %in% exclude] <- "withheld"
  ok <- rel %in% names(dirs)
  for (i in which(ok))
    .release_copy(file.path(tables_dir, f[i]), file.path(tables_dir, dirs[[rel[i]]], out_name[i]))
  cap <- if (length(f)) unname(captions[slug]) else character(0)
  cap[is.na(cap)] <- ""
  man <- data.frame(file = ifelse(ok, file.path(unname(dirs[ifelse(ok, rel, names(dirs)[1])]),
                                                out_name), ""),
                    source = f, data_id = id, release = rel, caption = cap, stringsAsFactors = FALSE)
  man <- man[order(!nzchar(man$release), man$release, substr(man$data_id, 1, 1),
                   suppressWarnings(as.integer(substring(man$data_id, 2))), man$file), , drop = FALSE]
  utils::write.csv(man, file.path(tables_dir, "data_release_manifest.csv"), row.names = FALSE, na = "")
  invisible(man)
}

# ---- supplement tables -------------------------------------------------------

RUN_DESIGN_LABELS <- c(forecast = "Forecast", counterfactual = "Counterfactual",
                       attribution_factual = "Attribution factual", events = "Past events",
                       hindcast = "Hindcast", smoothing = "Smoothing")

#' D34: the run design (report_analog_sets$run_design) as printed: one row
#' per arm present (config read), labelled.
report_run_design_table <- function(rd) {
  if (is.null(rd) || !nrow(rd)) return(NULL)
  rd <- rd[!is.na(rd$config), , drop = FALSE]
  na_dash <- function(x) ifelse(is.na(x), "—", as.character(x))
  tibble::tibble(
    Arm = unname(ifelse(rd$arm %in% names(RUN_DESIGN_LABELS), RUN_DESIGN_LABELS[rd$arm], rd$arm)),
    Replicates = fmt_int(rd$n_replicates), `Posterior draws (max_draws)` = fmt_int(rd$max_draws),
    `Replicates per draw` = ifelse(is.na(rd$reps_per_draw), "—", fmt_num(rd$reps_per_draw, 0)),
    `Draw mode` = na_dash(rd$draw_mode), Handoff = na_dash(rd$handoff),
    `Analogs per model` = fmt_int(rd$n_analogs_per_model), Models = fmt_int(rd$n_models),
    `Clamp margin` = ifelse(is.na(rd$clamp_margin), "off", fmt_num(rd$clamp_margin, 2)))
}

#' D35: the analysis set by admin level (report_dataset$table1, Zheng et al.
#' 2022 layout), ADM0 rows dropped, as printed.
report_analysis_set_table <- function(t1) {
  if (is.null(t1) || !nrow(t1)) return(NULL)
  t1 <- t1[as.character(t1$admin_level) != "ADM0", , drop = FALSE]
  has_dur <- all(c("dur_mo_med", "n_episodes") %in% names(t1))
  out <- tibble::tibble(
    `Admin level` = as.character(t1$admin_level),
    Outbreaks = fmt_int(t1$n_outbreaks), Locations = fmt_int(t1$n_gids),
    `Season-years` = fmt_int(t1$n_years), `Total cases` = fmt_cases(t1$total_cases),
    `Cases, median [IQR]` = fmt_iqr(t1$cases_med, t1$cases_q25, t1$cases_q75, fmt_int),
    `Duration (weeks), median [IQR]` = fmt_iqr(t1$dur_wk_med, t1$dur_wk_q25, t1$dur_wk_q75, fmt_num, digits = 1),
    `Population (thousands), median [IQR]` = fmt_iqr(t1$pop_med / 1e3, t1$pop_q25 / 1e3, t1$pop_q75 / 1e3, fmt_int),
    `Attack rate /1,000, median [IQR]` = fmt_iqr(t1$ar_med, t1$ar_q25, t1$ar_q75, fmt_num, digits = 2))
  if (has_dur) {
    out$`Episodes (duration set)` <- fmt_int(t1$n_episodes)
    out$`Duration (months), median [IQR]` <- fmt_iqr(t1$dur_mo_med, t1$dur_mo_q25, t1$dur_mo_q75, fmt_num, digits = 1)
  }
  out
}

#' D35b: the analysis set by admin level (incl. All) for Africa and each
#' WHO-AFRO subregion (report_dataset v5 table1_region), as printed; `region`
#' kept for 11's pack_by.
report_analysis_set_region_table <- function(t1r) {
  if (is.null(t1r) || !nrow(t1r)) return(NULL)
  tibble::tibble(
    region = t1r$region, `Admin level` = t1r$admin_level,
    Outbreaks = fmt_int(t1r$n_outbreaks), Locations = fmt_int(t1r$n_gids),
    `Season-years` = fmt_int(t1r$n_years), `Total cases` = fmt_cases(t1r$total_cases),
    `Cases, median [IQR]` = fmt_iqr(t1r$cases_med, t1r$cases_q25, t1r$cases_q75, fmt_int),
    `Duration (weeks), median [IQR]` = fmt_iqr(t1r$dur_wk_med, t1r$dur_wk_q25, t1r$dur_wk_q75, fmt_num, digits = 1),
    `Population (thousands), median [IQR]` = fmt_iqr(t1r$pop_med / 1e3, t1r$pop_q25 / 1e3, t1r$pop_q75 / 1e3, fmt_int),
    `Attack rate /1,000, median [IQR]` = fmt_iqr(t1r$ar_med, t1r$ar_q25, t1r$ar_q75, fmt_num, digits = 2))
}

#' D35c: outbreaks by season-year ENSO state for Africa and each WHO-AFRO
#' subregion (report_dataset v5 enso_region + enso_region_p). `raw = TRUE`: the
#' CSV frame, one row per block x state with the block's Kruskal-Wallis p
#' (p_cases, p_duration_weeks, p_attack_rate_1k) joined on; otherwise as
#' printed, `region` kept for 11's pack_by.
report_enso_region_table <- function(er, erp, raw = FALSE) {
  if (is.null(er) || !nrow(er)) return(NULL)
  p <- if (is.null(erp) || !nrow(erp)) NULL else
    tidyr::pivot_wider(erp, names_from = "variable", values_from = "p_value", names_prefix = "p_")
  out <- if (is.null(p)) er else dplyr::left_join(er, p, by = "region")
  if (isTRUE(raw)) return(out)
  fp <- function(x) ifelse(is.na(x), "—", ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f")))
  tibble::tibble(
    region = out$region, `ENSO state` = out$enso_state,
    `Outbreaks (share)` = sprintf("%s (%s%%)", fmt_int(out$n_outbreaks), fmt_num(100 * out$share, 1)),
    `Outbreaks per season-year` = if ("per_year" %in% names(out)) fmt_num(out$per_year, 1) else "—",
    `Season-years` = fmt_int(out$n_years),
    `Cases, median [IQR]` = fmt_iqr(out$cases_med, out$cases_q25, out$cases_q75, fmt_int),
    `Duration (weeks), median [IQR]` = fmt_iqr(out$dur_wk_med, out$dur_wk_q25, out$dur_wk_q75, fmt_num, digits = 1),
    `Attack rate /1,000, median [IQR]` = fmt_iqr(out$ar_med, out$ar_q25, out$ar_q75, fmt_num, digits = 2),
    `p cases / duration / attack rate` = if (is.null(p)) "—" else
      sprintf("%s / %s / %s", fp(out$p_cases), fp(out$p_duration_weeks), fp(out$p_attack_rate_1k)))
}

# ---- supplement tables T1-T5 ------------------------------------------------

#' T1a: the curated ENSO-rainfall mechanisms (results.inputs.enso_rain_mechanisms,
#' copied by 09 into rain_enso$mechanisms) with the members of each group from
#' the fit registry, in `order`.
#' @param mech rain_enso$mechanisms: list(path, groups)
#' @param members named list / character: group id -> ISO3 codes
report_t1a_mechanisms <- function(mech, members, order) {
  y <- mech$groups
  miss <- setdiff(order, names(y))
  if (length(miss))
    stop("report_t1a_mechanisms: ", paste(mech$path, collapse = ""), " has no entry for group(s) ",
         paste(miss, collapse = ", "), call. = FALSE)
  f <- function(k) vapply(order, function(g) as.character(y[[g]][[k]] %||% ""), "", USE.NAMES = FALSE)
  tibble::tibble(enso_group = as.character(order),
                 members = vapply(order, function(g) paste(sort(unlist(members[[g]])), collapse = " "), "",
                                  USE.NAMES = FALSE),
                 members_note = f("members_note"), season = f("season"), stage = f("stage"),
                 main_lag = f("main_lag"), el_nino_signal = f("el_nino_signal"),
                 mechanism = f("mechanism"), references = f("references"))
}

#' T1b: per group (in `order`), the rainfall columns 09 derives (rain_enso
#' v2 `group_lags`, rain_enso_lag_summary(), `weighting`) and the deployed
#' cholera lags (report_lags$selected) looked up beside them.
report_t1b_rain_lags <- function(rain, selected, order, weighting = "pop") {
  gl <- rain$group_lags
  if (is.null(gl))
    stop("report_t1b_rain_lags: the rain_enso artefact has no group_lags (schema v2); re-run 09 step rain_enso.",
         call. = FALSE)
  gl <- gl[gl$weighting == weighting, , drop = FALSE]
  gl <- gl[match(order, gl$enso_group), setdiff(names(gl), "weighting"), drop = FALSE]
  gl$enso_group <- as.character(order)
  dl <- function(trk) vapply(order, function(g) {
    x <- selected$deployed_lag[selected$enso_group == g & selected$track == trk]
    if (length(x)) as.integer(x[1]) else NA_integer_
  }, integer(1), USE.NAMES = FALSE)
  gl$lag_occ <- dl("occurrence"); gl$lag_dur <- dl("duration"); gl$lag_mon <- dl("monthly")
  tibble::as_tibble(gl)
}

#' T3: a condition outside the group crossbasis's fitted ONI range
#' (report_effects v6 range_lo / range_hi; NA when the artefact predates them).
report_outside_range <- function(enso_val, range_lo, range_hi, tol = 1e-9)
  !is.na(range_lo) & !is.na(range_hi) & (enso_val < range_lo - tol | enso_val > range_hi + tol)

#' T3 cell mark of a condition outside the fitted range (the ERF is flat past
#' its boundary knots there); "§", since "†" marks a low p_max and "‡" an EMRO
#' member in the same tables.
RANGE_MARK <- "§"
