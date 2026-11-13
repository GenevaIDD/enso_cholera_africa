# test_supp_missing_figs.R -- the supplement figures
# (10 --figures DC,SV,FA,AA), their builders on synthetic phase-09 artefacts.
#
# HARD gate (see analysis/tests/run_all.R). Pins:
#   [1] DC: data_coverage_frames() counts per country x season year and per
#       admin level x season year; build_figS_data_coverage() builds, rows in
#       ENSO-group blocks with the display labels
#   [2] SV: 09's smoothing_validation_tables() (smoothing v4; v5 onsets per
#       unit-season) on a hand-checked
#       two-replicate example -- durations, active fraction, the three
#       contagion levels, PIT counts, seasonality -- pit_hist_rebin() and the
#       builder
#   [3] FA: forecast_analogue_frame() (attribution factual vs per-set
#       counterfactual, window) and the builder
#   [4] AA: attribution_analogue_map_sf() (headlined cell, magnitudes
#       collapsed, 0-vs-0 out) and the builder (infinite ratio at the cap)
#   [5] PC: processing_cascade_frames() (stamps, nested differences, missing
#       stamps as --) and the flow-diagram builder
#
# Synthetic: no artefact, no data/.
# Run: Rscript analysis/tests/test_supp_missing_figs.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/figures_supp.R")
ok <- function(p) !inherits(try(suppressMessages(ggplot2::ggplot_build(p)), silent = TRUE), "try-error")
ok_patch <- function(p) {
  f <- tempfile(fileext = ".png")
  r <- try(suppressMessages(ggplot2::ggsave(f, p, width = 180, height = 200, units = "mm", dpi = 40)),
           silent = TRUE)
  !inherits(r, "try-error") && file.exists(f)
}
TH <- ggplot2::theme_bw(base_size = 7)
GL <- c(enso_west = "West [GW]", enso_south = "South [GS]")

cat("# [1] DC: data coverage\n")
ob <- tibble(country = c("NGA", "NGA", "NGA", "GHA", "MOZ", "MOZ", "ZZZ"),
             enso_group = c("enso_west", "enso_west", "enso_west", "enso_west", "enso_south", "enso_south", NA),
             admin_level = c("ADM1", "ADM2", "ADM1", "ADM1", "ADM2", "ADM2", "ADM1"),
             year_seas = c(2010L, 2010L, 2012L, 2011L, 2015L, NA, 2016L),
             cases = c(100, 50, 10, 5, 1000, 3, 7))
fr <- data_coverage_frames(ob)
expect(identical(fr$country$n_episodes[fr$country$country == "NGA" & fr$country$year_seas == 2010L], 2L) &&
         identical(fr$country$cases[fr$country$country == "NGA" & fr$country$year_seas == 2010L], 150) &&
         !anyNA(fr$country$year_seas) && sum(fr$country$n_episodes) == 6L,
       "DC: episodes / cases per country x season year, no-season rows dropped")
expect(identical(sum(fr$level$n_episodes[fr$level$admin_level == "ADM2"]), 2L),
       "DC: episodes per admin level x season year")
p <- build_figS_data_coverage(fr, c("enso_west", "enso_south"), GL, TH,
                              country_lab = c(NGA = "Nigeria", GHA = "Ghana", MOZ = "Mozambique"))
pa <- p[[1]]; ba <- ggplot2::ggplot_build(pa)
expect(ok_patch(p) && identical(levels(pa$data$enso_group), c("enso_west", "enso_south", "(ungrouped)")) &&
         all(c("Nigeria", "Ghana", "ZZZ") %in% levels(pa$data$cname)),
       "DC: builds; blocks in display order (ungrouped last), country names (ISO3 fallback)")
expect(ok(build_figS_data_coverage(fr, c("enso_west", "enso_south"), GL, TH, fill = "cases")[[1]]),
       "DC: case fill builds")

cat("# [2] SV: smoothing validation (09 smoothing v4 + builder)\n")
.md <- new.env()
source("analysis/utils/base.R", local = .md)
ex <- parse("analysis/utils/results_derive.R", keep.source = FALSE)
for (e in ex) if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
                  as.character(e[[2]]) %in% c("smoothing_validation_tables", ".month_index", ".with_seed")) eval(e, .md)
months <- seq(as.Date("2015-01-01"), by = "month", length.out = 6)
frame <- tibble(gid = c("AAA_1_1", "AAA_1_2", "AAA_2_1", "BBB_1_1", "CCC_1_1"),
                adm1_unit = c("AAA_1", "AAA_1", "AAA_2", "BBB_1", "CCC_1"),
                country = c("AAA", "AAA", "AAA", "BBB", "CCC"), tu = c("t1", "t1", "t1", "t1", NA))
# Replicate 1: AAA_1_1 active Jan-Feb (onset Jan), AAA_1_2 onset + active Feb.
# Replicate 2: nothing.
ev_on <- tibble(rep_id = 1L, gid = c("AAA_1_1", "AAA_1_2"), t_start = months[c(1, 2)])
ev_mon <- tibble(rep_id = 1L, gid = c("AAA_1_1", "AAA_1_1", "AAA_1_2"), t_start = months[c(1, 2, 2)],
                 cases_month = c(5, 7, 3))
cells <- tidyr::expand_grid(gid = frame$gid, t_start = months) |>
  mutate(obs_onset = as.integer(gid == "AAA_1_1" & t_start == months[1]))
obs_mon <- tibble(gid = "AAA_1_1", t_start = months[1:3], obs_cases = c(10, 20, NA))
scopes <- list(Africa = frame$gid, GW = c("AAA_1_1", "AAA_1_2", "AAA_2_1"))
v <- .md$smoothing_validation_tables(ev_on, ev_mon, 1:2, frame, cells, obs_mon, scopes)
du <- filter(v$duration, scope == "Africa")
expect(isTRUE(all.equal(du$q50[du$duration %in% 1:2], c(0.5, 0.5))) && identical(du$obs[du$duration == 3], 1) &&
         identical(du$n_obs[du$duration == 3], 1),
       "SV: episode durations (runs of active months): replicate shares, observed run of 3")
ac <- filter(v$active, scope == "Africa")
expect(isTRUE(all.equal(c(ac$obs, ac$q50, ac$n_cells), c(3 / 30, 1.5 / 30, 30))),
       "SV: fraction of covered unit-months active, observed and across replicates (zeros counted)")
co <- filter(v$contagion, scope == "Africa")
g <- function(i, s, col) co[[col]][co$indicator == i & co$state == s]
# Replicate 1, same ADM1: only AAA_1_2 in Feb is free (AAA_1_1 was active
# itself) -> 1 onset / 1 cell; replicate 2 has no such cell (NaN, dropped).
# Observed: AAA_1_1 active Jan-Mar -> AAA_1_2 Feb-Apr, 3 cells, no onset.
expect(identical(g("adm1", "ongoing", "q50"), 1) && identical(g("adm1", "ongoing", "obs"), 0) &&
         identical(g("adm1", "ongoing", "n_obs"), 3),
       "SV: onset given an ongoing outbreak in the same ADM1 (unit itself active last month excluded)")
# TU "none": 25 cells after January, minus CCC (no TU, 5), AAA_1_1 Feb-Apr
# (active itself, 3) and BBB_1_1 Feb-Apr (ongoing, 3) = 14.
expect(identical(g("country", "ongoing", "q50"), 0) && identical(g("country", "ongoing", "n_obs"), 3) &&
         identical(g("tu", "ongoing", "n_obs"), 3) && identical(g("tu", "none", "n_obs"), 14) &&
         identical(g("adm1", "none", "n_obs"), 19),
       "SV: other ADM1 of the country / other country of the TU; no-TU units left out")
ph <- filter(v$pit_hist, scope == "Africa")
expect(identical(sum(ph$n[ph$metric == "onsets"]), v$n_at_risk) && v$n_at_risk == 30L &&
         identical(sum(ph$n[ph$metric == "cases"]), 29L) && max(ph$bin) == 100L,
       "SV: PIT counts every at-risk cell (onsets) and every cell with known cases")
on_jan <- filter(v$seasonal, scope == "Africa", month == 1)
expect(isTRUE(all.equal(c(on_jan$obs, on_jan$q50), c(1 / 5, 0.5 / 5))) && all(is.na(filter(v$seasonal, month > 6)$obs)),
       "SV: onset frequency by calendar month; months without cells NA")
expect(is.null(v$season_onsets) && is.na(v$max_season_onsets),
       "SV: no onsets-per-unit-season table without min_seas")
# Onsets per unit-season (smoothing v5). Seasons from min_seas: AAA_1_1
# starts in March (Jan-Feb 2015 belong to the 2014 season), the rest in January.
# Replicate 1: AAA_1_1 onsets in Jan and Feb (one season), AAA_1_2 in Feb and
# April (one season), BBB_1_1 in Jan; replicate 2: AAA_1_1 in Jan and March.
frame_ms <- mutate(frame, min_seas = c(3L, 1L, 1L, 1L, 1L))
ev_on2 <- tibble(rep_id = c(1L, 1L, 1L, 1L, 1L, 2L, 2L),
                 gid = c("AAA_1_1", "AAA_1_1", "AAA_1_2", "AAA_1_2", "BBB_1_1", "AAA_1_1", "AAA_1_1"),
                 t_start = months[c(1, 2, 2, 4, 1, 1, 3)])
v2 <- .md$smoothing_validation_tables(ev_on2, ev_mon, 1:2, frame_ms, cells, obs_mon, scopes)
so <- filter(v2$season_onsets, scope == "Africa")
# Unit-seasons: AAA_1_1 x {2014, 2015}, the others x 2015 -> 6.
# Rep 1: 2014/AAA_1_1 = 2, AAA_1_2 = 2, BBB_1_1 = 1, three with 0 -> 0: 3/6, 1: 1/6, 2: 2/6.
# Rep 2: 2014/AAA_1_1 = 1, 2015/AAA_1_1 = 1 -> 0: 4/6, 1: 2/6, 2: 0.
expect(identical(v2$n_unit_seasons, 6L) && identical(v2$max_season_onsets, 2L) &&
         isTRUE(all.equal(so$q50, c(3.5, 1.5, 1, 0) / 6)) &&
         isTRUE(all.equal(so$obs, c(5, 1, 0, 0) / 6)) && identical(so$n_obs, c(5, 1, 0, 0)) &&
         all(so$n_unit_seasons == 6),
       "SV: onsets per unit-season from min_seas, observed and across replicates")
expect(identical(filter(v2$season_onsets, scope == "GW")$n_unit_seasons[1], 4),
       "SV: unit-seasons counted per scope")
expect(identical(v2[setdiff(names(v), c("season_onsets", "max_season_onsets", "n_unit_seasons"))][c("pit_hist", "active", "contagion", "duration")],
                 .md$smoothing_validation_tables(ev_on2, ev_mon, 1:2, frame, cells, obs_mon, scopes)[c("pit_hist", "active", "contagion", "duration")]),
       "SV: min_seas changes only the unit-season elements")
rb <- pit_hist_rebin(v$pit_hist, 20L)
expect(all(abs(tapply(rb$density, paste(rb$scope, rb$metric), mean) - 1) < 1e-9) &&
         inherits(try(pit_hist_rebin(v$pit_hist, 30L), silent = TRUE), "try-error"),
       "SV: PIT re-binned to densities averaging 1; a non-dividing bin count refused")
sv <- build_figS_smoothing_validation(v, c("Africa", "GW"), TH)
expect(ok_patch(sv), "SV: figure builds")

cat("# [3] FA: forecast vs counterfactual by analogue set\n")
qcols <- function(d, base) {
  for (m in c("cases", "onsets")) for (k in c("q025", "q10", "q25", "q50", "q75", "q90", "q975"))
    d[[paste0(m, "_", k)]] <- base * (if (m == "cases") 100 else 1) *
      c(q025 = 0.5, q10 = 0.6, q25 = 0.8, q50 = 1, q75 = 1.2, q90 = 1.4, q975 = 1.6)[[k]]
  d
}
mo <- seq(as.Date("2025-10-01"), by = "month", length.out = 8)
SETS4 <- c("all", "nina", "neutral", "elnino")
reg <- bind_rows(
  tidyr::expand_grid(run = "forecast", analog_set = NA_character_, month = mo, enso_group = c("enso_west", "enso_south")),
  tidyr::expand_grid(run = "counterfactual", analog_set = SETS4, month = mo, enso_group = c("enso_west", "enso_south")))
reg <- qcols(reg, ifelse(reg$run == "forecast", 3, 2))
afr <- qcols(bind_rows(tidyr::expand_grid(run = "forecast", analog_set = NA_character_, month = mo),
                       tidyr::expand_grid(run = "counterfactual", analog_set = SETS4, month = mo)), 5)
f23 <- list(region_ribbon = reg, africa_ribbon = afr, attrib_start = mo[2], attrib_factual = "forecast")
d <- forecast_analogue_frame(f23, c("enso_west", "enso_south"), c("nina", "elnino"))
expect(setequal(unique(d$analog_set), c("nina", "elnino")) && min(d$month) == mo[2] &&
         nrow(d) == 2 * 2 * 3 * 7 * 2 && identical(unique(d$q50[d$arm == "factual" & d$scope == "enso_west" & d$metric == "cases"]), 300),
       "FA: factual repeated under each set, counterfactual per set, from the attribution start, both metrics")
f23a <- c(f23, list(attrib_factual = "attrib",
                    attrib_ribbon = qcols(dplyr::filter(reg, run == "forecast"), 7),
                    attrib_africa_ribbon = qcols(dplyr::filter(afr, run == "forecast"), 9)))
f23a$attrib_factual <- "attrib"
da <- forecast_analogue_frame(f23a, "enso_west", "all", end = mo[5])
expect(identical(unique(da$q50[da$arm == "factual" & da$scope == "Africa" & da$metric == "onsets"]), 9) &&
         max(da$month) == mo[5] && !"enso_south" %in% da$scope,
       "FA: the attribution factual's ribbons when the artefact has one; window end and groups honoured")
fa <- build_figS_forecast_analogues(d, c(Africa = "Africa", enso_west = "GW", enso_south = "GS"),
                                    c(nina = "La Niña analogues", elnino = "El Niño analogues"), TH)
expect(ok_patch(fa), "FA: figure builds")

cat("# [4] AA: attribution maps by analogue set\n")
sq <- function(x, y) sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1), c(x, y + 1), c(x, y))))
adm1_sf <- sf::st_sf(GID_1 = c("AAA.1_1", "AAA.2_1", "BBB.1_1"), geometry = sf::st_sfc(sq(0, 0), sq(1, 0), sq(0, 1)))
adm1 <- tidyr::expand_grid(window = c("attribution", "post_fork"), explosive_guard = c(TRUE, FALSE),
                           metric = c("cases", "onsets"), analog_set = SETS4, magnitude = c(0.1, 0.5),
                           adm1_unit = c("AAA.1_1", "AAA.2_1", "BBB.1_1")) |>
  mutate(ratio = ifelse(adm1_unit == "BBB.1_1", NA, ifelse(adm1_unit == "AAA.2_1", Inf, 2)),
         p_sign = ifelse(analog_set == "nina", 0.9, 0.4) + ifelse(window == "attribution", 0, 0.05))
mm <- attribution_analogue_map_sf(adm1, "cases", c("nina", "elnino"), adm1_sf, window = "attribution", guard = TRUE)
expect(nrow(mm) == 4L && setequal(mm$analog_set, c("nina", "elnino")) && !"BBB.1_1" %in% mm$GID_1 &&
         identical(sort(unique(mm$p_sign)), c(0.4, 0.9)),
       "AA: one row per ADM1 x set at the headlined cell (magnitudes collapsed), 0-vs-0 ADM1s out")
afro <- sf::st_sf(ISO3_CODE = c("AAA", "BBB"), geometry = sf::st_sfc(sq(0, 0), sq(0, 1)))
rf <- ggplot2::scale_fill_gradient2(limits = log10(c(1 / 3, 3)), oob = scales::squish)
pf <- ggplot2::scale_fill_gradient2(limits = c(0, 1), midpoint = 0.5)
aa <- build_figS_attribution_analogues(
  list(cases = mm, onsets = attribution_analogue_map_sf(adm1, "onsets", c("nina", "elnino"), adm1_sf)),
  c(nina = "La Niña analogues", elnino = "El Niño analogues"), sf::st_bbox(afro), rf, pf,
  overlay_layers = list(ggplot2::geom_sf(data = afro, fill = NA)))
v1 <- aa[[1]]$layers[[1]]$data$value
expect(ok_patch(aa) && length(aa$patches$plots) + 1L == 4L && isTRUE(all.equal(max(v1), log10(3))),
       "AA: four map columns (ratio + probability per metric); an infinite ratio drawn at the cap")

cat("# [5] PC: processing cascade\n")
cas <- tibble(step = c(rep("00_ingest_outbreakextractr", 2), rep("01_preprocess_outbreaks", 2), "03_prepare_siz",
                       rep("03_prepare_occ", 3)),
              metric = c("n_rows_read", "n_rows_afr_emr", "n_outbreaks_in", "n_outbreaks_africa",
                         "n_episodes_canonical", "n_rows_out", "n_nz_rows", "n_events_out"),
              value = c(1000, 990, 50, 45, 30, 1234567, 12, 7))
pc <- processing_cascade_frames(cas)
expect(nrow(pc$main) == 6L && nrow(pc$term) == 3L &&
         grepl("Outside the two WHO regions: 10 records", pc$excl$label[1], fixed = TRUE) &&
         grepl("Outside Africa: 5 location-period episodes", pc$excl$label[2], fixed = TRUE) &&
         grepl("1,234,567 unit-months", pc$term$label[1], fixed = TRUE) &&
         grepl("30 episodes from -- location-period episodes", pc$main$label[6], fixed = TRUE),
       "PC: spine, side and dataset boxes from the stamps; nested differences; a missing stamp prints --")
expect(ok(build_figS_processing_cascade(pc)), "PC: the flow diagram builds")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
