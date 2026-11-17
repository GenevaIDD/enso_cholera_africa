# test_national_zero_fill.R — HARD: national zero-fill of the occurrence
#   panel.
#
# Usage:
#   Rscript analysis/tests/test_national_zero_fill.R
#
# Details:
#   `episodes.national_zero_fill` (03_prepare_occ.R --national_zero_fill,
#   gid_season only) adds a no-onset row for every national-zero month of a
#   panel gid's country not yet a panel row, never at or after an onset in
#   its season-year (occ_data.R::national_zero_surv_rows()). The variant
#   carries an "_nz" token on occurrence panel / fit / base-rate artefacts,
#   the loader checks the prep stamp, and the flag rides through generated
#   sim configs. Off (default), nothing changes. make_surv_data() is stubbed
#   (no SQLite); no data/, no network, no MCMC.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lubridate)
  library(purrr)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_data.R")
source("analysis/tests/helpers.R")

PASS <- 0L
FAIL <- 0L
errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      ""
    },
    error = function(e) conditionMessage(e)
  )
}
d <- as.Date

# Stub of make_surv_data(): one row per calendar month from the season-year start
# to TL_seas, last row = event (as make_surv_data()).
stub_build <- function(x, baseline) {
  ms <- x$min_seas[1]
  start_seas <- floor_date(x$TL_seas[1], "year")
  m <- seq(start_seas, floor_date(x$TL_seas[1], "month"), by = "month") %m+%
    months(ms - 1)
  tibble(
    uid = x$uid[1],
    gid = x$gid[1],
    location_period_id = x$location_period_id[1],
    t_start = m,
    t_end = ceiling_date(m, "month") - 1,
    t_start_seas = as.numeric(m - (start_seas %m+% months(ms - 1))),
    outbreak = seq_along(m) == length(m),
    stub_col = "x"
  )
}

# g1 (AAA, min_seas 1): panel = season 2015 Jan-Dec (no onset), season 2016
#   Jan-Mar with an onset in March.
# g2 (AAA, min_seas 7): panel = season 2015 = Jul 2015 - Jun 2016, no onset,
#   coverage ending Sep 2015 (rows Jul-Sep).
# g3 (BBB): no national-zero months.
mk <- function(gid, months, onset = NA) {
  tibble(
    gid = gid,
    t_start = d(months),
    outbreak = !is.na(onset) & d(months) == d(onset)
  )
}
surv <- bind_rows(
  mk("g1", seq(d("2015-01-01"), d("2015-12-01"), by = "month")),
  mk(
    "g1",
    seq(d("2016-01-01"), d("2016-03-01"), by = "month"),
    onset = "2016-03-01"
  ),
  mk("g2", seq(d("2015-07-01"), d("2015-09-01"), by = "month")),
  mk("g3", d("2015-01-01"))
)
ids <- tibble(
  gid = c("g1", "g2", "g3"),
  uid = c("u1", "u2", "u3"),
  location_period_id = c("lp1", "lp2", "lp3"),
  min_seas = c(1L, 7L, 1L),
  country = c("AAA", "AAA", "BBB")
)
nz <- tibble(
  country = "AAA",
  month = d(c(
    "2015-05-01", # g1: already a panel row; g2: before its season (season 2014)
    "2015-11-01", # g1: panel row; g2: after coverage end, no onset -> add
    "2016-05-01",
    "2016-06-01", # g1: after the 2016 onset -> drop; g2: add
    "2017-02-01",
    "2017-03-01"
  )), # new season-years for both
  source = "national_and_who"
)

cat("\n[1] national_zero_surv_rows()\n")
r <- national_zero_surv_rows(
  surv,
  nz,
  ids,
  baseline = "1990_2019",
  build = stub_build
)
key <- function(x) paste(x$gid, format(x$t_start))
expect(
  setequal(
    key(r),
    c(
      "g2 2015-05-01",
      "g2 2015-11-01",
      "g2 2016-05-01",
      "g2 2016-06-01",
      "g1 2017-02-01",
      "g1 2017-03-01",
      "g2 2017-02-01",
      "g2 2017-03-01"
    )
  ),
  "adds new season-years, extends a no-onset season, keeps a pre-season month of the gid's own calendar"
)
expect(!any(key(r) %in% key(surv)), "no existing panel cell is duplicated")
expect(
  !any(key(r) %in% c("g1 2016-05-01", "g1 2016-06-01")) &&
    attr(r, "n_after_onset") == 2L,
  "months after an onset in the same season-year are dropped (2)"
)
expect(
  !any(r$gid == "g3"),
  "a country without national-zero months gets nothing"
)
expect(
  all(!r$outbreak) && all(r$national_zero_fill),
  "added rows: outbreak FALSE, national_zero_fill TRUE"
)
ref <- stub_build(
  tibble(
    gid = "g1",
    uid = "u1",
    location_period_id = "lp1",
    min_seas = 1L,
    TL_seas = d("2017-12-31")
  ),
  "1990_2019"
) %>%
  filter(t_start %in% d(c("2017-02-01", "2017-03-01")))
got <- r %>% filter(gid == "g1") %>% arrange(t_start)
g_cmp <- as.data.frame(select(got, -outbreak, -national_zero_fill))
r_cmp <- as.data.frame(select(ref, -outbreak))
expect(
  identical(names(g_cmp), names(r_cmp)) &&
    isTRUE(all.equal(g_cmp, r_cmp, check.attributes = FALSE)),
  "rows are the builder's rows for those months (columns, t_start_seas, ids)"
)
g2s <- r %>% filter(gid == "g2", t_start == d("2015-11-01"))
expect(
  identical(g2s$t_start_seas, as.numeric(d("2015-11-01") - d("2015-07-01"))),
  "season month follows the gid's own calendar (min_seas 7)"
)
expect(
  attr(r, "n_season_years") == 4L && attr(r, "n_new_season_years") == 3L,
  "season-years: 4 touched (g2 2014, 2015, 2016; g1 2017), 3 new (g2 2015 already had rows)"
)
r0 <- national_zero_surv_rows(
  surv,
  nz[0, ],
  ids,
  baseline = "1990_2019",
  build = stub_build
)
expect(
  nrow(r0) == 0L && "national_zero_fill" %in% names(r0),
  "no national-zero months -> no rows"
)

cat("\n[2] path tokens\n")
o_gs <- list(baseline = "1990_2019", episodes = list(key = "gid_season"))
o_nz <- list(
  baseline = "1990_2019",
  episodes = list(key = "gid_season", national_zero_fill = TRUE)
)
o_bad <- list(
  baseline = "1990_2019",
  episodes = list(key = "lp_season", national_zero_fill = TRUE)
)
expect(
  identical(
    as.character(make_surv_dat_file(list(baseline = "1990_2019"))),
    "generated_data/outbreaks_survival_data_base-1990_2019.rds"
  ) &&
    identical(
      as.character(make_surv_dat_file(o_gs)),
      "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs.rds"
    ),
  "default and gid_season survival paths unchanged"
)
expect(
  identical(
    as.character(make_surv_dat_file(o_nz)),
    "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs_nz.rds"
  ),
  "national zero-fill survival path carries _ep-gs_nz"
)
expect(
  grepl("_ep-gs_nz", make_occ_run_file(o_nz, "occ_base_rate"), fixed = TRUE) &&
    !grepl("_nz", make_occ_run_file(o_gs, "occ_base_rate"), fixed = TRUE),
  "occurrence base-rate table carries _nz only with the flag"
)
expect(
  identical(
    as.character(make_monthly_dat_file(o_nz)),
    as.character(make_monthly_dat_file(o_gs))
  ) &&
    identical(
      as.character(make_sim_panel_file(o_nz)),
      as.character(make_sim_panel_file(o_gs))
    ),
  "monthly frame and sim climate panel take no _nz token"
)
expect(
  grepl(
    "requires episodes.key = 'gid_season'",
    errmsg(make_surv_dat_file(o_bad))
  ),
  "national_zero_fill under lp_season is refused"
)
fits <- c(
  "outbreaks_survival_brms_cb_fit_x_ep-gs.rds",
  "outbreaks_survival_brms_cb_fit_x_ep-gs_nz.rds"
)
expect(
  identical(fits[.national_zero_fill_file_match(fits, TRUE)], fits[2]) &&
    identical(fits[.national_zero_fill_file_match(fits, FALSE)], fits[1]),
  ".national_zero_fill_file_match() separates the two occurrence fits"
)
o_t60 <- list(
  baseline = "1990_2019",
  episodes = list(
    key = "gid_season",
    national_zero_fill = TRUE,
    national_zero_threshold = 60L
  )
)
expect(
  identical(
    as.character(make_surv_dat_file(o_t60)),
    "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs_nz60.rds"
  ) &&
    grepl(
      "_ep-gs_nz60",
      make_occ_run_file(o_t60, "occ_base_rate"),
      fixed = TRUE
    ),
  "threshold variant: survival path and base rate carry _ep-gs_nz60"
)
expect(
  grepl("national_zero_months_t60.rds$", make_national_zero_months_file(60L)) &&
    grepl("national_zero_months.rds$", make_national_zero_months_file()),
  "make_national_zero_months_file(): _t<N> for a threshold, unchanged at 0"
)
fits3 <- c(
  "outbreaks_survival_brms_cb_fit_x_ep-gs.rds",
  "outbreaks_survival_brms_cb_fit_x_ep-gs_nz.rds",
  "outbreaks_survival_brms_cb_fit_x_ep-gs_nz60.rds",
  "outbreaks_survival_brms_cb_fit_x_sl0_ep-gs_nz60_extra.rds"
)
expect(
  identical(which(.national_zero_fill_file_match(fits3, TRUE, 60L)), 3:4) &&
    identical(which(.national_zero_fill_file_match(fits3, TRUE)), 2L) &&
    identical(which(.national_zero_fill_file_match(fits3, FALSE)), 1L),
  "exact token match: _nz and _nz60 never match each other"
)
expect(
  grepl(
    "requires episodes.national_zero_fill = true",
    errmsg(make_surv_dat_file(list(
      baseline = "1990_2019",
      episodes = list(key = "gid_season", national_zero_threshold = 60L)
    )))
  ) &&
    nzchar(errmsg(national_zero_threshold_of(list(
      episodes = list(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = -1
      )
    )))),
  "threshold without zero-fill, or negative, is refused"
)

cat("\n[3] assert_national_zero_fill()\n")
td <- tempfile("nzf")
dir.create(td)
p_nz <- file.path(td, "surv_ep-gs_nz.rds")
saveRDS(tibble(), p_nz)
persist_prep_policy_stamp(
  p_nz,
  prep_policy_stamp(
    "crosswalk",
    "drop_collapse",
    0.9,
    episode_key = "gid_season",
    counts = list(national_zero_fill = TRUE)
  )
)
p_gs <- file.path(td, "surv_ep-gs.rds")
saveRDS(tibble(), p_gs)
persist_prep_policy_stamp(
  p_gs,
  prep_policy_stamp(
    "crosswalk",
    "drop_collapse",
    0.9,
    episode_key = "gid_season"
  )
)
expect(
  isTRUE(assert_national_zero_fill(o_nz, p_nz)) &&
    identical(assert_national_zero_fill(o_gs, p_gs), FALSE),
  "matching stamp and config pass"
)
expect(
  grepl("prepped without", errmsg(assert_national_zero_fill(o_nz, p_gs))) &&
    grepl("prepped with", errmsg(assert_national_zero_fill(o_gs, p_nz))),
  "a stamp disagreeing with the config is refused"
)
p_t60 <- file.path(td, "surv_ep-gs_nz60.rds")
saveRDS(tibble(), p_t60)
persist_prep_policy_stamp(
  p_t60,
  prep_policy_stamp(
    "crosswalk",
    "drop_collapse",
    0.9,
    episode_key = "gid_season",
    counts = list(national_zero_fill = TRUE, national_zero_threshold = 60L)
  )
)
expect(
  isTRUE(assert_national_zero_fill(o_t60, p_t60)) &&
    isTRUE(assert_national_zero_fill(o_nz, p_nz)),
  "threshold stamp matches; a stamp without the threshold field reads as 0"
)
expect(
  grepl(
    "threshold = 60 but .* threshold 0",
    errmsg(assert_national_zero_fill(o_t60, p_nz))
  ) &&
    grepl(
      "threshold = 0 but .* threshold 60",
      errmsg(assert_national_zero_fill(o_nz, p_t60))
    ),
  "a threshold disagreeing with the stamp is refused"
)
unlink(td, recursive = TRUE)

cat("\n[4] apply_episode_key_to_sim_config()\n")
cfg <- list(simulator = list(run_tag = "x"))
expect(
  identical(
    apply_episode_key_to_sim_config(cfg, "gid_season")$episodes,
    list(key = "gid_season")
  ),
  "default: episodes block exactly as before"
)
c_nz <- apply_episode_key_to_sim_config(
  cfg,
  "gid_season",
  national_zero_fill = TRUE
)
expect(
  identical(c_nz$episodes, list(key = "gid_season", national_zero_fill = TRUE)),
  "national_zero_fill = TRUE is written"
)
expect(
  identical(
    apply_episode_key_to_sim_config(c_nz, "gid_season")$episodes,
    c_nz$episodes
  ),
  "national_zero_fill = NULL keeps the template's flag"
)
expect(
  grepl(
    "deploy declares false",
    errmsg(apply_episode_key_to_sim_config(
      c_nz,
      "gid_season",
      national_zero_fill = FALSE
    ))
  ),
  "a template disagreeing with the deploy is refused"
)
expect(
  grepl(
    "requires episodes.key",
    errmsg(apply_episode_key_to_sim_config(
      cfg,
      "lp_season",
      national_zero_fill = TRUE
    ))
  ),
  "national_zero_fill under lp_season is refused"
)
c_t60 <- apply_episode_key_to_sim_config(
  cfg,
  "gid_season",
  national_zero_fill = TRUE,
  national_zero_threshold = 60L
)
expect(
  identical(
    c_t60$episodes,
    list(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L
    )
  ) &&
    identical(
      apply_episode_key_to_sim_config(
        c_nz,
        "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 0L
      )$episodes,
      c_nz$episodes
    ),
  "threshold written only when > 0"
)
expect(
  grepl(
    "deploy declares 0",
    errmsg(apply_episode_key_to_sim_config(
      c_t60,
      "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 0L
    ))
  ) &&
    grepl(
      "requires episodes.national_zero_fill",
      errmsg(apply_episode_key_to_sim_config(
        cfg,
        "gid_season",
        national_zero_threshold = 60L
      ))
    ),
  "a template threshold disagreeing with the deploy, or a threshold without zero-fill, is refused"
)


cat(
  "\n[item 12] named rule: fill kept out of own episode spans; post-mask gate\n"
)
# g2 (min_seas 7) had an episode in season 2015 running Jun 2016 - Jun 2016's
# season boundary: TL 2016-05-15, TR 2016-08-31 -> season 2015 (Jul15-Jun16),
# running into season 2016 (Jul-Aug 2016). Candidate months 2016-05/06 of g2 lie
# inside it (the legacy fill adds them: no onset row of g2 in the panel).
ep12 <- tibble(
  gid = c("g1", "g2"),
  TL = d(c("2016-03-01", "2016-05-15")),
  TR = d(c("2016-04-30", "2016-08-31")),
  year_seas = c(2016L, 2015L),
  min_seas = c(1L, 7L)
)
r_leg <- national_zero_surv_rows(
  surv,
  nz,
  ids,
  baseline = "1990_2019",
  build = stub_build
)
r_new <- national_zero_surv_rows(
  surv,
  nz,
  ids,
  baseline = "1990_2019",
  build = stub_build,
  episodes = ep12
)
expect(
  identical(r_leg, r),
  "episodes = NULL (legacy) -> identical to the default call (attributes included)"
)
expect(
  is.null(attr(r_leg, "n_in_episode")),
  "legacy: no n_in_episode attribute"
)
expect(
  setequal(
    key(r_leg)[!key(r_leg) %in% key(r_new)],
    c("g2 2016-05-01", "g2 2016-06-01")
  ) &&
    identical(attr(r_new, "n_in_episode"), 2L),
  "named rule: g2's months inside its own running episode are not added (2)"
)
full_leg <- bind_rows(
  mutate(
    surv,
    t_end = ceiling_date(t_start, "month") - 1,
    national_zero_fill = FALSE
  ),
  r_leg
)
full_new <- bind_rows(
  mutate(
    surv,
    t_end = ceiling_date(t_start, "month") - 1,
    national_zero_fill = FALSE
  ),
  r_new
)
expect(
  grepl(
    "inside their gid's own episode span",
    errmsg(assert_nz_rows_outside_episodes(full_leg, ep12))
  ),
  "gate asserts under a named rule (legacy fill would trip it)"
)
n_rep <- suppressMessages(assert_nz_rows_outside_episodes(
  full_leg,
  ep12,
  assert = FALSE
))
expect(
  identical(n_rep, 2L),
  "gate is report-only under the legacy rule (count returned)"
)
expect(
  identical(errmsg(assert_nz_rows_outside_episodes(full_new, ep12)), ""),
  "named-rule fill passes the gate"
)

cat("\n[item 11/12] .own_episode_span_hits()\n")
rows <- tibble(
  gid = c("g2", "g2", "g2", "g2", "g1", "g9"),
  t_start = d(c(
    "2016-05-01",
    "2016-07-01",
    "2016-08-01",
    "2016-09-01",
    "2016-04-01",
    "2016-07-01"
  )),
  t_end = d(c(
    "2016-05-15",
    "2016-07-31",
    "2016-08-31",
    "2016-09-30",
    "2016-04-30",
    "2016-07-31"
  ))
)
expect(
  identical(
    .own_episode_span_hits(rows, ep12, "span"),
    c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
  ),
  "span: month-level overlap with an own episode; other gids never hit"
)
expect(
  identical(
    .own_episode_span_hits(rows, ep12, "running"),
    c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  ),
  "running: only rows of a LATER season inside (TL, TR] (the onset row and same-season rows are not hit)"
)
expect(
  identical(.own_episode_span_hits(rows[0, ], ep12, "span"), logical(0)) &&
    identical(
      .own_episode_span_hits(rows, ep12[0, ], "running"),
      rep(FALSE, 6L)
    ),
  "empty inputs"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
