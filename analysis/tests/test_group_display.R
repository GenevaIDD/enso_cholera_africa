# test_group_display.R — HARD: group display labels, the group_labels step and
#   the per-country season-end and event-window tables.
#
# Usage:
#   Rscript analysis/tests/test_group_display.R
#
# Details:
#   Pins group_display_meta() [1], [1b], relabel_groups() [2], the labels from
#   the deployed fits [3], the season-end tables of table S30 [4] and the
#   per-country event windows of report_event_af_model [5]. Synthetic: no fit,
#   no data/, no simulation output.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("analysis/utils/results_core.R")
for (.m in c("season_calendar", "loo_lags", "kl_sweep")) {
  source(file.path("analysis/utils", paste0(.m, ".R")))
}
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

REG <- "enso_grouping_prod_run_prior"
meta <- load_grouping_meta(REG)
LAGS_FX <- c(
  enso_east = 0L,
  enso_horn = 10L,
  enso_west = -6L,
  enso_nga = 7L,
  enso_south = 4L
)

# registry = the sidecar's names; lag = letters by the deployed OCCURRENCE
# lag, longest lead first, ties in sidecar order; "Class A (ENSO lead 6)" /
# "A -6" (true minus) / "Class D (ENSO lag +7)"; ids and colours unchanged;
# bad mode / missing lag refused.
cat("# [1] group_display_meta()\n")
dr <- group_display_meta(meta, mode = "registry")
expect(
  identical(unname(dr$label), unname(meta$label[meta$order])) &&
    identical(unname(dr$short_fig), unname(meta$short_fig[meta$order])) &&
    identical(dr$order, meta$order),
  "registry: the sidecar's labels, short names and order"
)
dl <- group_display_meta(meta, LAGS_FX, "lag")
expect(
  identical(
    unname(dl$letter[c(
      "enso_west",
      "enso_east",
      "enso_south",
      "enso_nga",
      "enso_horn"
    )]),
    LETTERS[1:5]
  ),
  "lag: A West, B East, C South, D NGA, E Horn (ascending occurrence lag)"
)
expect(
  identical(
    dl$order,
    c("enso_west", "enso_east", "enso_south", "enso_nga", "enso_horn")
  ),
  "lag: display order = letter order"
)
expect(
  identical(dl$label[["enso_west"]], "Class A (ENSO lead 6)") &&
    identical(dl$label[["enso_nga"]], "Class D (ENSO lag +7)") &&
    identical(dl$label[["enso_east"]], "Class B (ENSO lag 0)"),
  "lag: 'Class A (ENSO lead 6)', 'Class D (ENSO lag +7)', 'Class B (ENSO lag 0)'"
)
expect(
  identical(dl$short[["enso_west"]], "A −6") &&
    identical(dl$short[["enso_horn"]], "E +10") &&
    identical(dl$short[["enso_east"]], "B 0") &&
    identical(dl$short_fig, dl$short),
  "lag: short 'A −6' (true minus), 'E +10', 'B 0'; short_fig = short"
)
expect(
  identical(names(dl$label), dl$order) &&
    identical(names(dl$short_fig), dl$order),
  "lag: label vectors named in display order (legends / facets list A, B, C, ...)"
)
expect(
  identical(dl$colour, dr$colour) && identical(names(dl$colour), meta$order),
  "colours keyed by id, unchanged by the mode"
)
tie <- group_display_meta(
  meta,
  c(
    enso_east = 2L,
    enso_horn = 2L,
    enso_west = 1L,
    enso_nga = 3L,
    enso_south = 2L
  ),
  "lag"
)
expect(
  identical(
    unname(tie$letter[c(
      "enso_west",
      "enso_east",
      "enso_horn",
      "enso_south",
      "enso_nga"
    )]),
    LETTERS[1:5]
  ),
  "ties keep the sidecar order"
)
expect(
  throws(
    group_display_meta(meta, LAGS_FX[-1], "lag"),
    "no deployed occurrence lag"
  ),
  "lag mode with a group missing its lag stops"
)
expect(
  throws(group_display_meta(meta, NULL, "lag"), "needs each group"),
  "lag mode without lags stops"
)
expect(
  throws(
    group_display_meta(meta, LAGS_FX, "letters"),
    "group_labels must be one of"
  ),
  "an unknown mode stops"
)
gl_tab <- tibble(enso_group = names(LAGS_FX), lag_occ = unname(LAGS_FX))
expect(
  identical(group_display_meta(meta, gl_tab, "lag")$label, dl$label),
  "lags given as the group_labels artefact's table"
)

# Centroid mode: GW / GG / GH / GE / GS from the sidecar, "West [GW]", display
# order centroid_order; missing sidecar fields refused.
cat("# [1b] group_display_meta(): centroid\n")
dc <- group_display_meta(meta, mode = "centroid")
expect(
  identical(
    dc$order,
    c("enso_west", "enso_nga", "enso_horn", "enso_east", "enso_south")
  ) &&
    identical(names(dc$label), dc$order) &&
    identical(names(dc$short), dc$order),
  "centroid: display order GW, GG, GH, GE, GS (sidecar centroid_order), vectors named in that order"
)
expect(
  identical(unname(dc$short), c("GW", "GG", "GH", "GE", "GS")) &&
    identical(dc$short_fig, dc$short),
  "centroid: short GW / GG / GH / GE / GS; short_fig = short"
)
expect(
  identical(dc$label[["enso_nga"]], "Guinean Gulf [GG]") &&
    identical(dc$label[["enso_west"]], "West [GW]") &&
    identical(dc$label[["enso_south"]], "South [GS]"),
  "centroid: labels 'West [GW]', 'Guinean Gulf [GG]', 'South [GS]'"
)
expect(
  identical(dc$colour, dr$colour) && is.null(dc$lag) && all(is.na(dc$letter)),
  "centroid: colours keyed by id, no lags or letters"
)
expect(
  identical(group_display_meta(meta, LAGS_FX, "centroid")$label, dc$label),
  "centroid ignores lags"
)
expect(
  identical(
    relabel_groups(
      c("enso_horn", "Nigeria+Ghana", "NGA+GHA", "South", "Africa"),
      dc
    ),
    c("Horn [GH]", "GG", "GG", "GS", "Africa")
  ),
  "centroid relabel: id -> label, stored short / short_fig -> GG / GS, Africa untouched"
)
meta_nc <- meta
meta_nc$centroid_short[["enso_east"]] <- NA_character_
expect(
  throws(
    group_display_meta(meta_nc, mode = "centroid"),
    "needs centroid_short"
  ),
  "centroid without a sidecar centroid_short stops"
)
meta_no <- meta
meta_no$centroid_order <- character(0)
expect(
  identical(group_display_meta(meta_no, mode = "centroid")$order, meta$order),
  "centroid without centroid_order: sidecar order"
)
m_prior <- load_grouping_meta("enso_grouping_prod_run_prior")
expect(
  identical(group_display_meta(m_prior, mode = "centroid")$label, dc$label),
  "enso_grouping_prod_run_prior sidecar carries the same centroid names"
)

# Identity under registry; ids -> label, stored short and short_fig -> display
# short, factor levels, column names; other strings untouched.
cat("# [2] relabel_groups()\n")
df <- tibble(
  enso_group = c("enso_west", "enso_south", NA),
  region = c("West", "Southern", "Africa"),
  scope = factor(
    c("South", "West", "Africa"),
    levels = c("Africa", "West", "South")
  ),
  West = 1:3,
  n = 1:3
)
expect(identical(relabel_groups(df, dr), df), "registry: identity")
rl <- relabel_groups(df, dl)
expect(
  identical(
    rl$enso_group,
    c("Class A (ENSO lead 6)", "Class C (ENSO lag +4)", NA)
  ),
  "ids -> label"
)
expect(
  identical(rl$region, c("A −6", "C +4", "Africa")),
  "stored short -> display short; Africa untouched"
)
expect(
  identical(levels(rl$scope), c("Africa", "A −6", "C +4")) &&
    identical(as.character(rl$scope), c("C +4", "A −6", "Africa")),
  "factor levels mapped in place (short_fig 'South' and short 'West')"
)
expect(
  identical(names(rl), c("enso_group", "region", "scope", "A −6", "n")),
  "column names mapped"
)
expect(
  identical(names(relabel_groups(df, dl, rename_cols = FALSE)), names(df)),
  "rename_cols = FALSE keeps names"
)
expect(
  identical(relabel_groups(df, dl, cols = "region")$enso_group, df$enso_group),
  "cols restricts the columns"
)
expect(
  identical(relabel_groups(c("Horn", "x"), dl), c("E +10", "x")),
  "a character vector is mapped"
)

# group_deployed_lags() / derive_group_labels() on a synthetic best-lag archive
# built from a forecast manifest's fit file names: the lags are read from the
# fit files, a disagreeing archive or manifest stops; the fixture gives
# A = West (lead 6), B = East, C = South, D = NGA, E = Horn.
cat("# [3] group_labels from the deployed fits\n")
man <- yaml::read_yaml("analysis/tests/fixtures/forecast_horizons_fixture.yml")
td <- tempfile("gl")
dir.create(td)
arch <- file.path(td, "arch")
dir.create(arch)
code <- c(occ = "", dur = "dur_", mon = "mon_")
rows <- list()
for (g in names(man$groups)) {
  for (k in names(code)) {
    f <- basename(man$groups[[g]][[paste0(k, "_fit_file")]])
    file.create(file.path(arch, f))
    lag <- as.integer(man$groups[[g]][[paste0("best_lag_", k)]])
    rows[[length(rows) + 1L]] <- data.frame(
      track = k,
      group = gsub("_", "", g),
      lag = lag,
      kind = c(paste0(code[[k]], "fit"), paste0(code[[k]], "cb_meta")),
      file = c(f, sub("_fit_", "_cb_meta_", f))
    )
  }
}
write.csv(
  do.call(rbind, rows),
  file.path(arch, "MANIFEST.csv"),
  row.names = FALSE
)
mk_ctx <- function(MAN = man, dir = arch) {
  list2env(
    list(
      MAN = MAN,
      ENSO_BEST_LAG_DIR = dir,
      ENSO_BEST_LAG_MANIFEST = file.path(dir, "MANIFEST.csv"),
      GROUPING_PATH = file.path("analysis/groupings", paste0(REG, ".csv")),
      MANIFEST_PATH = "analysis/tests/fixtures/forecast_horizons_fixture.yml",
      opt = list(
        episodes = list(
          key = "gid_season",
          national_zero_fill = TRUE,
          national_zero_threshold = 60L
        )
      ),
      HINDCAST_MANIFEST_PATH = "",
      EVENT_MANIFEST_PATH = "",
      SURV_PATH = file.path(td, "surv.rds"),
      RES = list(deploy_config = "")
    ),
    parent = emptyenv()
  )
}
ctx <- mk_ctx()
lg <- group_deployed_lags(ctx, meta$order)
expect(
  identical(lg$lag_occ, unname(LAGS_FX[meta$order])) &&
    identical(lg$lag_dur[lg$enso_group == "enso_west"], -4L) &&
    identical(lg$lag_mon[lg$enso_group == "enso_east"], 2L),
  "occ / dur / mon lags read from the fit files (slm<k> = lead)"
)
bad_man <- man
bad_man$groups$enso_horn$best_lag_occ <- 9L
expect(
  throws(group_deployed_lags(mk_ctx(bad_man), meta$order), "forecast manifest"),
  "a forecast manifest disagreeing with the fit files stops"
)
m2 <- read.csv(file.path(arch, "MANIFEST.csv"))
m2$lag[m2$track == "mon" & m2$group == "ensonga"] <- 3L
arch2 <- file.path(td, "arch2")
dir.create(arch2)
invisible(file.copy(list.files(arch, full.names = TRUE), arch2))
write.csv(m2, file.path(arch2, "MANIFEST.csv"), row.names = FALSE)
expect(
  throws(
    group_deployed_lags(mk_ctx(dir = arch2), meta$order),
    "archive MANIFEST"
  ),
  "an archive MANIFEST disagreeing with the fit files stops"
)
unlink(list.files(arch2, "outbreaks_survival.*enso_horn", full.names = TRUE))
expect(
  throws(
    group_deployed_lags(mk_ctx(dir = arch2), meta$order),
    "no deployed occurrence fit"
  ),
  "a missing occurrence fit stops"
)
gl <- suppressMessages(derive_group_labels(ctx))
expect(
  identical(gl$groups$enso_group, meta$order) &&
    all(gl$groups$in_manifest) &&
    sum(gl$groups$n_countries) == nrow(gl$members) &&
    nrow(gl$members) == 30L,
  "group_labels: groups in sidecar order, 30 members, counts"
)
expect(
  !any(vapply(
    gl$groups,
    function(x) is.character(x) && any(grepl("Class|ENSO lead", x)),
    logical(1)
  )),
  "group_labels stores no display string"
)
expect(
  all(is.na(gl$season_ends$season_end_month)) &&
    all(is.na(gl$season_ends$calendar_default)) &&
    !nrow(gl$season_windows),
  "without panel / manifests: season columns NA, no windows"
)
dl2 <- group_display_meta(meta, gl$groups, "lag")
expect(
  identical(dl2$label, dl$label),
  "labels built from the artefact equal the fixture's A-E"
)

# season_ends_country_table() / season_windows_table(): per-country season
# ends, the calendar-default flag NA without the 03 sidecar, per-country firm
# horizons preferred over the group's.
cat("# [4] season-end tables\n")
cal <- data.frame(
  gid = c("KEN.1", "KEN.2", "KEN.3", "TZA.1", "TZA.2", "MOZ.1"),
  country = c("KEN", "KEN", "KEN", "TZA", "TZA", "MOZ"),
  min_seas = c(8L, 8L, 3L, 1L, 6L, 10L)
)
se <- season_end_by_country(cal)
cgr <- c(
  KEN = "enso_east",
  TZA = "enso_east",
  MOZ = "enso_south",
  BDI = "enso_east"
)
st <- season_ends_country_table(se, cgr, c(enso_east = 7L, enso_south = 4L))
expect(
  identical(st$iso3, c("BDI", "KEN", "MOZ", "TZA")) &&
    identical(st$season_end_month, c(NA, 7L, 9L, 12L)) &&
    identical(st$group_season_end_month, c(7L, 7L, 4L, 7L)) &&
    isTRUE(all.equal(st$modal_share[2], 2 / 3)),
  "per country: end = modal min_seas - 1, ties -> the later end month (TZA 1/6 -> 12), group month kept"
)
expect(
  all(is.na(st$calendar_default)),
  "calendar_default NA without the 03 .min_seas_source.csv"
)
src <- data.frame(
  gid = cal$gid,
  min_seas_source = c(
    "estimate_gid",
    "default_january",
    "estimate_adm1",
    "default_january",
    "default_january",
    "estimate_gid"
  )
)
st2 <- season_ends_country_table(
  se,
  cgr,
  NULL,
  src,
  stats::setNames(cal$country, cal$gid)
)
expect(
  identical(st2$n_gid_default, c(NA, 1L, 0L, 2L)) &&
    identical(st2$calendar_default, c(NA, FALSE, FALSE, TRUE)),
  "with the sidecar: default counts and an all-default flag"
)
hc <- list(
  groups = list(
    enso_east = list(
      iso3 = list("KEN", "TZA"),
      years = list(
        `2015` = list(
          launch = "2015-05-01",
          horizon_end = "2016-07-01",
          season_end = "2016-07-01"
        )
      )
    )
  ),
  units = list(
    u1 = list(
      members = list("KEN", "TZA"),
      groups = list("enso_east"),
      years = list(`2015` = list(horizon_end = "2016-09-01"))
    )
  ),
  countries = list(
    KEN = list(
      group = "enso_east",
      unit = "u1",
      years = list(
        `2015` = list(
          season_end = "2016-07-01",
          firm_horizon_end = "2016-07-01"
        )
      )
    ),
    TZA = list(
      group = "enso_east",
      unit = "u1",
      years = list(
        `2015` = list(
          season_end = "2016-05-01",
          firm_horizon_end = "2016-05-01"
        )
      )
    )
  )
)
hw <- season_windows_table(hc, "hindcast")
expect(
  identical(hw$firm_horizon_end, as.Date(c("2016-07-01", "2016-05-01"))) &&
    identical(hw$unit_horizon_end, as.Date(c("2016-09-01", "2016-09-01"))) &&
    identical(hw$season_end, as.Date(c("2016-07-01", "2016-05-01"))) &&
    identical(hw$launch, as.Date(c("2015-05-01", "2015-05-01"))),
  "hindcast windows: each country's own firm horizon and season end, the unit's horizon"
)
hc_g <- hc
hc_g$countries <- lapply(hc$countries, function(x) x["group"])
hc_g$countries$KEN$unit <- hc_g$countries$TZA$unit <- "u1"
expect(
  identical(
    season_windows_table(hc_g, "hindcast")$firm_horizon_end,
    as.Date(c("2016-07-01", "2016-07-01"))
  ),
  "by_group manifest: every country at its group's horizon"
)

# manifest_country_horizons() on an event view.
cat("# [5] event windows per country (report_event_af_model, group_windows)\n")
em <- read_event_manifest(
  "analysis/tests/fixtures/event_attribution_horizons_fixture.yml"
)
ew <- season_windows_table(em, "event")
expect(
  nrow(ew) == length(em$events) * 30L &&
    all(
      ew$firm_horizon_end ==
        as.Date(vapply(
          seq_len(nrow(ew)),
          function(i) {
            em$events[[ew$window[i]]]$groups[[ew$enso_group[i]]]$horizon_end
          },
          ""
        ))
    ),
  "by_group event manifest: country windows = their group's horizon_end"
)
ev1 <- names(em$events)[1]
em_c <- em
em_c$events[[ev1]]$countries$KEN$firm_horizon_end <- "2016-07-01"
ce <- as.Date(manifest_country_horizons(event_manifest_view(em_c, ev1), "firm"))
expect(
  identical(ce[["KEN"]], as.Date("2016-07-01")) &&
    identical(
      ce[["TZA"]],
      as.Date(em$events[[ev1]]$groups$enso_east$horizon_end)
    ),
  "a country's own firm_horizon_end wins over the group's (the AF model's row filter)"
)
obs <- tibble(
  country = c("KEN", "KEN", "TZA"),
  t_start = as.Date(c("2016-07-01", "2016-08-01", "2016-08-01"))
)
kept <- dplyr::filter(obs, t_start <= ce[country])
expect(
  identical(kept$country, c("KEN", "TZA")),
  "rows past a country's own window end are dropped"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
