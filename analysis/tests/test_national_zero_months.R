# test_national_zero_months.R — HARD: national-zero months in the coverage
#   mask.
#
# Usage:
#   Rscript analysis/tests/test_national_zero_months.R
#
# Details:
#   A country-month is a national zero (utils/national_zero.R) when
#     ((national 0 AND no WHO > 0) OR (WHO 0 AND no national > 0))
#       AND no sub-national cases AND no active outbreak AND not excluded,
#   and build_coverage_mask() (utils/sim_io.R) covers every panel gid of that
#   country in those months -- on by default. A wrong month here counts
#   simulated cholera against an "observed zero" that never was. Pure
#   synthetic tables; no data/, no network, no MCMC.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lubridate)
  library(purrr)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/coverage.R")
source("analysis/utils/national_zero.R")
source("analysis/utils/who_annual.R")
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

cat("\n[1] national_month_status()\n")
weeks <- tibble(
  country = "XXX",
  TL = d(c(
    "2020-01-06",
    "2020-01-13",
    "2020-01-20",
    "2020-01-27",
    "2020-03-02"
  )),
  sCh = c(0, 0, 0, 0, 3)
) %>%
  mutate(TR = TL + 6)
ns <- national_month_status(weeks)
st <- function(m) ns[ns$month == d(m), ]
expect(
  isTRUE(st("2020-01-01")$nat_zero),
  "Jan: zero weeks cover 26/31 days -> national zero"
)
expect(
  nrow(st("2020-02-01")) == 1L &&
    !st("2020-02-01")$nat_zero &&
    !st("2020-02-01")$nat_pos,
  "Feb: zero weeks cover 2/29 days -> neither zero nor positive"
)
expect(
  isTRUE(st("2020-03-01")$nat_pos) && !st("2020-03-01")$nat_zero,
  "Mar: a week with cases -> positive"
)
expect(
  !isTRUE(national_month_status(weeks, zero_day_frac = 0.9)[1, ]$nat_zero),
  "zero_day_frac 0.9: 26/31 days is not enough"
)

cat("\n[2] build_national_zero_months()\n")
nat <- tibble(
  country = c("AAA", "BBB", "CCC", "DDD", "DDD", "DDD"),
  month = d(c(
    "2020-01-01",
    "2020-01-01",
    "2020-02-01",
    "2020-01-01",
    "2020-02-01",
    "2020-03-01"
  )),
  nat_pos = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
  nat_zero = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
)
who <- tibble(
  country = c("BBB", "CCC", "DDD", "ETH", "ETH"),
  year = c(2020L, 2020L, 2020L, 2012L, 2015L),
  suspected_cases = c(5, 0, 0, 0, 0),
  is_zero = c(FALSE, TRUE, TRUE, TRUE, TRUE)
)
sub_pos <- tibble(country = "DDD", month = d("2020-02-01"))
active <- tibble(country = "DDD", month = d("2020-03-01"))
z <- build_national_zero_months(
  nat,
  sub_pos,
  active,
  who,
  start = d("2012-01-01"),
  end = d("2020-03-01"),
  exclude = parse_national_zero_exclude("ETH:2011-2014")
)
has <- function(cc, m) any(z$country == cc & z$month == d(m))
src <- function(cc, m) z$source[z$country == cc & z$month == d(m)]
expect(
  has("AAA", "2020-01-01") &&
    identical(src("AAA", "2020-01-01"), "national_who_missing"),
  "national 0, WHO missing -> zero (national_who_missing)"
)
expect(!has("BBB", "2020-01-01"), "national 0 but WHO > 0 -> not zero")
expect(
  has("CCC", "2020-01-01") &&
    identical(src("CCC", "2020-01-01"), "who_no_national"),
  "WHO 0, no national report -> zero (who_no_national)"
)
expect(!has("CCC", "2020-02-01"), "WHO 0 but national > 0 -> not zero")
expect(
  has("DDD", "2020-01-01") &&
    identical(src("DDD", "2020-01-01"), "national_and_who"),
  "national 0 and WHO 0 -> zero (national_and_who)"
)
expect(!has("DDD", "2020-02-01"), "sub-national cases that month -> not zero")
expect(
  !has("DDD", "2020-03-01"),
  "active database outbreak that month -> not zero"
)
expect(
  !any(z$country == "ETH" & year(z$month) == 2012L) &&
    sum(z$country == "ETH" & year(z$month) == 2015L) == 12L,
  "ETH 2012 excluded, ETH 2015 kept (all 12 months)"
)
expect(
  !has("CCC", "2020-04-01") && max(z$month) <= d("2020-03-01"),
  "nothing after `end`"
)
dr <- attr(z, "dropped")
expect(
  identical(unname(dr), c(1L, 1L, 12L)),
  "dropped counts: 1 active, 1 sub-national, 12 excluded"
)
expect(
  nrow(parse_national_zero_exclude("")) == 0L &&
    nrow(parse_national_zero_exclude("ETH:2011-2014,SOM:2020")) == 5L &&
    nzchar(errmsg(parse_national_zero_exclude("ETH2011"))),
  "parse_national_zero_exclude(): empty, ranges + single years, malformed refused"
)

cat("\n[3] build_coverage_mask()\n")
td <- tempfile("nz")
dir.create(td)
agg <- file.path(td, "agg.rds")
saveRDS(
  tibble(
    location_period_id = "lp1",
    admin_level = "ADM2",
    TL = d("2010-01-01"),
    TR = d("2010-03-31")
  ),
  agg
)
panel <- tibble(
  location_period_id = c("lp1", "lp2", "lp3"),
  gid = c("g1", "g2", "g3"),
  country = c("AAA", "AAA", "BBB"),
  t_start = d("2010-01-01")
)
nzp <- file.path(td, "nz.rds")
saveRDS(
  tibble(
    country = c("AAA", "ZZZ"),
    month = d(c("2011-01-01", "2011-01-01")),
    source = "national_and_who"
  ),
  nzp
)
m_off <- suppressMessages(build_coverage_mask(
  panel,
  agg,
  national_zero = FALSE
))
expect(
  nrow(m_off) == 3L && all(m_off$gid == "g1"),
  "off: window cells only (g1, Jan-Mar 2010)"
)
m_on <- suppressMessages(build_coverage_mask(
  panel,
  agg,
  national_zero_path = nzp
))
added <- anti_join(m_on, m_off, by = c("gid", "month"))
expect(
  identical(sort(added$gid), c("g1", "g2")) &&
    all(added$month == d("2011-01-01")),
  "on: the national-zero month covers every panel gid of the country (g1, g2), not BBB, not ZZZ"
)
expect(
  nrow(semi_join(m_off, m_on, by = c("gid", "month"))) == nrow(m_off),
  "on is a superset of off"
)
expect(
  grepl(
    "national-zero months not found",
    errmsg(suppressMessages(
      build_coverage_mask(
        panel,
        agg,
        national_zero_path = file.path(td, "missing.rds")
      )
    ))
  ),
  "on with a missing file is an error"
)
# national_zero_path is REQUIRED while national_zero is on (a default
# threshold-0 file would silently mismatch a thresholded run).
expect(
  grepl(
    "national_zero_path is required",
    errmsg(suppressMessages(
      build_coverage_mask(panel, agg)
    ))
  ),
  "on without national_zero_path is an error"
)
expect(
  identical(formals(build_coverage_mask)$national_zero_path, quote(expr = )),
  "build_coverage_mask() has no national_zero_path default"
)
# The LP member bridge. The panel carries one representative LP per
# gid; g2's own reporting window sits on its second member LP lp2b.
saveRDS(
  tibble(
    location_period_id = c("lp1", "lp2b"),
    admin_level = "ADM2",
    TL = d(c("2010-01-01", "2010-05-01")),
    TR = d(c("2010-03-31", "2010-06-30"))
  ),
  agg
)
sp <- file.path(td, "outbreaks_survival_data_base-1990_2019_ep-gs.rds")
expect(
  identical(coverage_lp_bridge(panel, sp), panel) &&
    identical(coverage_lp_bridge(panel, NULL), panel),
  "no sidecar (or no path): the panel itself is the bridge (pre-T3 behaviour)"
)
m_rep <- suppressMessages(build_coverage_mask(
  panel,
  agg,
  national_zero = FALSE
))
m_rep2 <- suppressMessages(build_coverage_mask(
  coverage_lp_bridge(panel, sp),
  agg,
  national_zero = FALSE
))
expect(
  identical(m_rep, m_rep2),
  "no sidecar: mask identical to the panel-bridged mask"
)
expect(
  !"g2" %in% m_rep$gid,
  "representative-LP bridge misses g2's member-LP window (M3)"
)
utils::write.csv(
  tibble(
    location_period_id = c("lp1", "lp2", "lp2b", "lp9"),
    gid = c("g1", "g2", "g2", "g9"),
    country = c("AAA", "AAA", "AAA", "CCC")
  ),
  occ_panel_lp_members_path(sp),
  row.names = FALSE
)
br <- suppressMessages(coverage_lp_bridge(panel, sp))
expect(
  setequal(
    paste(br$location_period_id, br$gid),
    c("lp1 g1", "lp2 g2", "lp2b g2")
  ),
  "sidecar: every member LP of the panel's gids (g9, not in the panel, dropped)"
)
m_mem <- suppressMessages(build_coverage_mask(br, agg, national_zero = FALSE))
expect(
  sum(m_mem$gid == "g2") == 2L && all(m_rep$gid %in% m_mem$gid),
  "sidecar: g2 covered May-Jun 2010 through lp2b; a superset of the representative mask"
)
expect(
  identical(
    as.character(occ_panel_lp_members_path("a/b_ep-gs_nz60.rds")),
    "a/b_ep-gs_nz60.lp_members.csv"
  ),
  "occ_panel_lp_members_path(): <panel>.lp_members.csv"
)
unlink(td, recursive = TRUE)

# Threshold variant ([4]): WHO country-years under zero_threshold cases (max
# over sources) add "who_below_threshold" months, the last `window` months only
# when the next year is also below; db_who_consistency() keeps under-reported
# WHO years out of it; zero_threshold = 0 is exactly the rule above ([5]).
cat("\n[4] threshold branch + db_who_consistency()\n")
# EEE: 2020 has 30 cases (below 60), 2021 has 500 -> Jan-Aug 2020 claimed,
#   Sep-Dec (window 4) unclaimed. FFF: 2020 = 40, 2021 = 0 (exact zero) -> all
#   12 months of 2020 claimed. GGG: 2020 = 40, 2021 missing -> tail unclaimed.
#   HHH: 2020 = 40 with an active outbreak in Mar 2020 -> March dropped.
#   III: 2020 = 40 but the database has 200 cases -> auto-excluded.
#   JJJ: 2020 = 0 (exact zero): no national rows -> who_no_national, not below.
who4 <- tibble(
  country = c("EEE", "EEE", "FFF", "FFF", "GGG", "HHH", "III", "JJJ"),
  year = c(2020L, 2021L, 2020L, 2021L, 2020L, 2020L, 2020L, 2020L),
  suspected_cases = c(30, 500, 40, 0, 40, 40, 40, 0),
  is_zero = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE),
  suspected_cases_max = c(30, 500, 45, 0, 40, 40, 40, 0)
)
nat0 <- tibble(
  country = character(),
  month = as.Date(character()),
  nat_pos = logical(),
  nat_zero = logical()
)
act4 <- tibble(country = "HHH", month = d("2020-03-01"))
ob4 <- tibble(
  country = c("III", "III"),
  TL = d(c("2020-02-01", "2019-12-01")),
  TR = d(c("2020-04-30", "2019-12-31")),
  suspected_cases = c(200, 999),
  admin_level = factor(c("ADM2", "ADM1"), levels = paste0("ADM", 0:3))
)
cons <- db_who_consistency(ob4, who4, tolerance = 1.2)
expect(
  nrow(cons) == 1L &&
    cons$country == "III" &&
    cons$year == 2020L &&
    abs(cons$db_deep - 200) < 1e-9 &&
    isTRUE(cons$excluded),
  "db_who_consistency: year-apportioned deepest-level cases vs WHO; III 2020 excluded (200 > 1.2 x 40)"
)
ob_span <- tibble(
  country = "KKK",
  TL = d("2020-12-01"),
  TR = d("2021-01-31"),
  suspected_cases = 62,
  admin_level = factor("ADM2", levels = paste0("ADM", 0:3))
)
cs <- db_who_consistency(
  ob_span,
  tibble(
    country = "KKK",
    year = c(2020L, 2021L),
    suspected_cases = 31,
    is_zero = FALSE,
    suspected_cases_max = 31
  )
)
expect(
  nrow(cs) == 2L && all(abs(cs$db_deep - 31) < 1e-9) && !any(cs$excluded),
  "db_who_consistency: an episode across the year end is split by days (31 + 31)"
)
z4 <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4,
  who4,
  start = d("2020-01-01"),
  end = d("2021-12-01"),
  zero_threshold = 60L,
  window = 4L,
  auto_exclude = cons
)
bm <- function(cc) {
  format(z4$month[z4$country == cc & z4$source == "who_below_threshold"], "%m")
}
expect(
  identical(bm("EEE"), sprintf("%02d", 1:8)),
  "EEE: below-threshold year, next year above -> Jan-Aug only"
)
expect(
  identical(bm("FFF"), sprintf("%02d", 1:12)) &&
    all(
      z4$source[z4$country == "FFF" & year(z4$month) == 2021L] ==
        "who_no_national"
    ),
  "FFF: next year an exact zero counts as below -> all 12 months; 2021 stays who_no_national"
)
expect(
  identical(bm("GGG"), sprintf("%02d", 1:8)),
  "GGG: next year missing -> tail unclaimed"
)
expect(
  identical(bm("HHH"), sprintf("%02d", c(1:2, 4:8))),
  "HHH: active outbreak month dropped"
)
expect(
  length(bm("III")) == 0L,
  "III: auto-excluded by the database consistency check"
)
expect(
  !any(z4$country == "JJJ" & z4$source == "who_below_threshold") &&
    sum(z4$country == "JJJ") == 12L,
  "JJJ: exact zero stays in its own branch (who_no_national)"
)
expect(
  identical(attr(z4, "tail_unclaimed"), 16L),
  "tail_unclaimed: 4 each for EEE (next above), GGG, HHH, III (next missing)"
)
dr4 <- attr(z4, "dropped")
expect(
  identical(
    names(dr4),
    c("active_outbreak", "subnational_cases", "excluded", "auto_excluded")
  ) &&
    dr4[["active_outbreak"]] == 1L &&
    dr4[["auto_excluded"]] == 8L,
  "dropped gains auto_excluded (III's 8 months); HHH's active month counted"
)
z4w <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4,
  who4,
  start = d("2020-01-01"),
  end = d("2021-12-01"),
  zero_threshold = 60L,
  window = 6L
)
expect(
  sum(z4w$country == "EEE") == 6L && sum(z4w$country == "III") == 6L,
  "window 6: six months claimed; no auto_exclude -> III kept"
)
expect(
  nzchar(errmsg(build_national_zero_months(
    nat0,
    sub_pos[0, ],
    act4,
    who4[, 1:4],
    start = d("2020-01-01"),
    end = d("2020-12-01"),
    zero_threshold = 60L
  ))),
  "threshold without suspected_cases_max is an error"
)
nat_pos <- tibble(
  country = "EEE",
  month = d("2020-02-01"),
  nat_pos = TRUE,
  nat_zero = FALSE
)
z4p <- build_national_zero_months(
  nat_pos,
  sub_pos[0, ],
  act4,
  who4[who4$country == "EEE", ],
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L
)
expect(
  !any(z4p$month == d("2020-02-01")),
  "a national report with cases blocks the below-threshold month"
)

cat("\n[5] zero_threshold = 0 reproduces the exact-zero rule\n")
z0 <- build_national_zero_months(
  nat,
  sub_pos,
  active,
  who,
  start = d("2012-01-01"),
  end = d("2020-03-01"),
  exclude = parse_national_zero_exclude("ETH:2011-2014"),
  zero_threshold = 0L,
  auto_exclude = cons
)
expect(
  identical(z0, z),
  "fixture: threshold 0 (with an auto_exclude table) is identical to the default call"
)
nz_disk <- as.character(make_national_zero_months_file())
ts_disk <- file.path(
  path.expand(default_outbreak_data_dir()),
  "time_series_preoutbreak_extraction.rds"
)
who_dir <- "data/who_annual_report_africa"
if (
  file.exists(nz_disk) &&
    file.exists(ts_disk) &&
    file.exists(make_outbreak_metadata_file()) &&
    file.exists(file.path(who_dir, "who_annual_report_africa.rds"))
) {
  st <- readRDS(sub("\\.rds$", ".policy.rds", nz_disk))
  ts <- readRDS(ts_disk) %>%
    ungroup() %>%
    mutate(TL = as.Date(TL), TR = as.Date(TR))
  wm <- build_who_annual_zero_mask(
    file.path(who_dir, "who_annual_report_africa.rds"),
    file.path(who_dir, "who_dashboard_annual_africa.rds")
  )
  zd <- build_national_zero_months(
    national_month_status(
      filter(ts, spatial_scale == "country"),
      st$zero_day_frac
    ),
    subnational_positive_months(filter(ts, spatial_scale != "country")),
    active_episode_months(readRDS(make_outbreak_metadata_file())),
    wm,
    start = as.Date(st$start),
    end = as.Date(st$end),
    exclude = parse_national_zero_exclude(st$exclude)
  )
  expect(
    identical(zd, readRDS(nz_disk)),
    "real inputs: threshold 0 reproduces the on-disk national_zero_months.rds"
  )
} else {
  cat("  (skip) real national-zero inputs not present\n")
}


# ---- rule registry, head guard, consistency scope ------
# Frozen pre-registry copy of build_national_zero_months(), for the identical()
# back-compat checks.
old_buildNationalZeroMonths <- function(
  nat,
  sub_pos,
  active,
  who_mask,
  start,
  end,
  exclude = tibble::tibble(country = character(), year = integer()),
  zero_threshold = 0L,
  window = 4L,
  auto_exclude = NULL
) {
  zero_threshold <- as.integer(zero_threshold)
  window <- as.integer(window)
  if (is.na(zero_threshold) || zero_threshold < 0L) {
    stop(
      "build_national_zero_months: zero_threshold must be >= 0.",
      call. = FALSE
    )
  }
  if (is.na(window) || window < 0L || window > 12L) {
    stop("build_national_zero_months: window must be in 0..12.", call. = FALSE)
  }
  months <- seq(
    lubridate::floor_date(as.Date(start), "month"),
    lubridate::floor_date(as.Date(end), "month"),
    by = "month"
  )
  countries <- sort(unique(c(nat$country, who_mask$country)))
  key <- function(d) paste(d$country, d$month)
  grid <- tidyr::expand_grid(country = countries, month = months) %>%
    dplyr::left_join(nat, by = c("country", "month")) %>%
    dplyr::mutate(
      year = as.integer(lubridate::year(month)),
      nat_pos = dplyr::coalesce(nat_pos, FALSE),
      nat_zero = dplyr::coalesce(nat_zero, FALSE)
    ) %>%
    dplyr::left_join(
      who_mask %>%
        dplyr::transmute(
          country,
          year = as.integer(year),
          who_zero = is_zero,
          who_pos = !is_zero
        ),
      by = c("country", "year")
    ) %>%
    dplyr::mutate(
      who_zero = dplyr::coalesce(who_zero, FALSE),
      who_pos = dplyr::coalesce(who_pos, FALSE),
      rule = (nat_zero & !who_pos) | (who_zero & !nat_pos),
      sub = key(.) %in% key(sub_pos),
      act = key(.) %in% key(active),
      excl = paste(country, year) %in% paste(exclude$country, exclude$year)
    )
  # Threshold branch (zero_threshold > 0): WHO below the threshold on every
  # source, not an exact zero (that is the who_zero branch), no national
  # report with cases; the last `window` months need the next year below too.
  tail_unclaimed <- 0L
  if (zero_threshold > 0L) {
    if (!"suspected_cases_max" %in% names(who_mask)) {
      stop(
        "build_national_zero_months: zero_threshold > 0 needs who_mask$suspected_cases_max ",
        "(build_who_annual_zero_mask()).",
        call. = FALSE
      )
    }
    below <- who_mask %>%
      dplyr::transmute(
        country,
        year = as.integer(year),
        below = suspected_cases_max < zero_threshold & !is_zero
      )
    # The next year counts as below when it is under the threshold, an exact
    # zero included (that year's own months come from the who_zero branch).
    nxt <- who_mask %>%
      dplyr::transmute(
        country,
        year = as.integer(year) - 1L,
        next_below = suspected_cases_max < zero_threshold
      )
    grid <- grid %>%
      dplyr::left_join(below, by = c("country", "year")) %>%
      dplyr::left_join(nxt, by = c("country", "year")) %>%
      dplyr::mutate(
        below = dplyr::coalesce(below, FALSE),
        next_below = dplyr::coalesce(next_below, FALSE),
        in_tail = as.integer(lubridate::month(month)) > 12L - window,
        below_rule = below & !nat_pos & !rule & (!in_tail | next_below)
      )
    tail_unclaimed <- sum(
      grid$below & !grid$nat_pos & !grid$rule & grid$in_tail & !grid$next_below
    )
    auto <- if (is.null(auto_exclude)) {
      character(0)
    } else {
      with(
        auto_exclude[auto_exclude$excluded, , drop = FALSE],
        paste(country, year)
      )
    }
    grid$auto <- grid$below_rule & paste(grid$country, grid$year) %in% auto
  } else {
    grid$below_rule <- FALSE
    grid$auto <- FALSE
  }
  r <- dplyr::filter(grid, rule | below_rule)
  dropped <- c(
    active_outbreak = sum(r$act),
    subnational_cases = sum(r$sub & !r$act),
    excluded = sum(r$excl & !r$sub & !r$act)
  )
  if (zero_threshold > 0L) {
    dropped <- c(
      dropped,
      auto_excluded = sum(r$auto & !r$excl & !r$sub & !r$act)
    )
  }
  out <- r %>%
    dplyr::filter(!sub, !act, !excl, !auto) %>%
    dplyr::transmute(
      country,
      month,
      source = dplyr::case_when(
        !rule ~ "who_below_threshold",
        nat_zero & who_zero ~ "national_and_who",
        nat_zero ~ "national_who_missing",
        TRUE ~ "who_no_national"
      )
    ) %>%
    dplyr::arrange(country, month)
  attr(out, "dropped") <- dropped
  if (zero_threshold > 0L) {
    attr(out, "tail_unclaimed") <- as.integer(tail_unclaimed)
  }
  out
}

cat("\n[6] head guard (head_window)\n")
# LLL: 2019 = 900 (high), 2020 = 30 -> head Jan-Apr unclaimed; MMM: 2019 = 20
#   (below) -> head claimed; NNN: 2019 = 0 (exact zero counts as below);
#   OOO: 2019 missing -> head unclaimed (as the tail). All 2021 = 0.
who6 <- tibble(
  country = rep(c("LLL", "MMM", "NNN", "OOO"), each = 3L),
  year = rep(2019:2021, 4L),
  suspected_cases = c(900, 30, 0, 20, 30, 0, 0, 30, 0, NA, 30, 0),
  is_zero = c(
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    NA,
    FALSE,
    TRUE
  ),
  suspected_cases_max = c(900, 30, 0, 20, 30, 0, 0, 30, 0, NA, 30, 0)
) %>%
  filter(!is.na(suspected_cases))
z6 <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  window = 4L,
  head_window = 4L
)
bm6 <- function(cc) {
  format(z6$month[z6$country == cc & z6$source == "who_below_threshold"], "%m")
}
expect(
  identical(bm6("LLL"), sprintf("%02d", 5:12)),
  "LLL: previous year high -> Jan-Apr unclaimed (head guard)"
)
expect(
  identical(bm6("MMM"), sprintf("%02d", 1:12)),
  "MMM: previous year below -> head claimed"
)
expect(
  identical(bm6("NNN"), sprintf("%02d", 1:12)),
  "NNN: previous year an exact zero counts as below"
)
expect(
  identical(bm6("OOO"), sprintf("%02d", 5:12)),
  "OOO: previous year missing -> head unclaimed"
)
expect(
  identical(attr(z6, "head_unclaimed"), 8L) &&
    identical(attr(z6, "tail_unclaimed"), 0L),
  "head_unclaimed = 8 (LLL, OOO); tail claimed (2021 an exact zero)"
)
z6t <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6[
    who6$country == "LLL" &
      who6$year != 2021L,
  ],
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  window = 4L,
  head_window = 4L
)
expect(
  identical(format(z6t$month, "%m"), sprintf("%02d", 5:8)) &&
    identical(attr(z6t, "head_unclaimed"), 4L) &&
    identical(attr(z6t, "tail_unclaimed"), 4L),
  "head and tail together: only May-Aug; a month is counted once (tail first)"
)
# head_missing_quiet (rule r3): OOO's missing 2019 now counts as below; LLL (a
# known high year) stays unclaimed; a missing NEXT year still blocks the tail.
z6q <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  window = 4L,
  head_window = 4L,
  head_missing_quiet = TRUE
)
bm6q <- function(cc) {
  format(
    z6q$month[z6q$country == cc & z6q$source == "who_below_threshold"],
    "%m"
  )
}
expect(
  identical(bm6q("OOO"), sprintf("%02d", 1:12)) &&
    identical(bm6q("LLL"), sprintf("%02d", 5:12)) &&
    identical(attr(z6q, "head_unclaimed"), 4L),
  "head_missing_quiet: a missing previous year counts as below (OOO), a high one does not (LLL)"
)
z6qt <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6[
    who6$country == "OOO" &
      who6$year != 2021L,
  ],
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  window = 4L,
  head_window = 4L,
  head_missing_quiet = TRUE
)
expect(
  identical(format(z6qt$month, "%m"), sprintf("%02d", 1:8)) &&
    identical(attr(z6qt, "tail_unclaimed"), 4L),
  "head_missing_quiet leaves the tail guard alone: a missing next year still blocks Sep-Dec"
)
expect(
  identical(
    build_national_zero_months(
      nat0,
      sub_pos[0, ],
      act4[0, ],
      who6,
      start = d("2020-01-01"),
      end = d("2020-12-01"),
      zero_threshold = 60L,
      window = 4L,
      head_window = 4L,
      head_missing_quiet = FALSE
    ),
    z6
  ),
  "head_missing_quiet = FALSE is the default (r2 output unchanged)"
)
el <- .year_edge_claim(
  d(sprintf("2020-%02d-01", 1:12)),
  4L,
  "tail",
  rep(FALSE, 12L)
)
eh <- .year_edge_claim(
  d(sprintf("2020-%02d-01", 1:12)),
  4L,
  "head",
  rep(FALSE, 12L)
)
expect(
  identical(which(el$in_edge), 9:12) &&
    identical(which(eh$in_edge), 1:4) &&
    identical(rev(el$in_edge), eh$in_edge),
  ".year_edge_claim(): head is the mirror of the tail"
)

cat("\n[7] consistency_scope = all\n")
# PPP 2020: WHO exact zero, no national report, database 300 cases -> its
#   who_no_national months are dropped under "all"; QQQ 2020: WHO zero AND a
#   national zero in Jan, database cases -> Jan (national_and_who) kept.
who7 <- tibble(
  country = c("PPP", "QQQ"),
  year = 2020L,
  suspected_cases = 0,
  is_zero = TRUE,
  suspected_cases_max = 0
)
ob7 <- tibble(
  country = c("PPP", "QQQ"),
  TL = d("2020-06-01"),
  TR = d("2020-06-30"),
  suspected_cases = 300,
  admin_level = factor("ADM2", levels = paste0("ADM", 0:3))
)
nat7 <- tibble(
  country = "QQQ",
  month = d("2020-01-01"),
  nat_pos = FALSE,
  nat_zero = TRUE
)
cons7 <- db_who_consistency(ob7, who7, tolerance = 1.2)
expect(
  all(cons7$excluded),
  "db_who_consistency flags a WHO exact zero with database cases"
)
z7t <- build_national_zero_months(
  nat7,
  sub_pos[0, ],
  act4[0, ],
  who7,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  auto_exclude = cons7
)
z7a <- build_national_zero_months(
  nat7,
  sub_pos[0, ],
  act4[0, ],
  who7,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  auto_exclude = cons7,
  consistency_scope = "all"
)
expect(
  sum(z7t$country == "PPP") == 12L,
  "scope threshold (legacy): the WHO exact zero is still claimed"
)
expect(
  !any(z7a$country == "PPP"),
  "scope all: WHO-zero months contradicted by the database are dropped"
)
expect(
  identical(z7a$month[z7a$country == "QQQ"], d("2020-01-01")) &&
    identical(z7a$source[z7a$country == "QQQ"], "national_and_who"),
  "scope all: a national-zero month resting on the national series is kept"
)
expect(
  identical(attr(z7a, "dropped")[["auto_excluded_exact"]], 23L) &&
    !"auto_excluded_exact" %in% names(attr(z7t, "dropped")),
  "dropped gains auto_excluded_exact (12 PPP + 11 QQQ) only under scope all"
)

cat(
  "\n[8] legacy values reproduce the pre-registry output (attributes included)\n"
)
calls <- list(
  list(
    nat,
    sub_pos,
    active,
    who,
    start = d("2012-01-01"),
    end = d("2020-03-01"),
    exclude = parse_national_zero_exclude("ETH:2011-2014")
  ),
  list(
    nat0,
    sub_pos[0, ],
    act4,
    who4,
    start = d("2020-01-01"),
    end = d("2021-12-01"),
    zero_threshold = 60L,
    window = 4L,
    auto_exclude = cons
  ),
  list(
    nat0,
    sub_pos[0, ],
    act4,
    who4,
    start = d("2020-01-01"),
    end = d("2021-12-01"),
    zero_threshold = 60L,
    window = 6L
  ),
  list(
    nat0,
    sub_pos[0, ],
    act4[0, ],
    who6,
    start = d("2020-01-01"),
    end = d("2020-12-01"),
    zero_threshold = 60L,
    window = 4L
  ),
  list(
    nat7,
    sub_pos[0, ],
    act4[0, ],
    who7,
    start = d("2020-01-01"),
    end = d("2020-12-01"),
    auto_exclude = cons7
  )
)
for (i in seq_along(calls)) {
  a <- calls[[i]]
  expect(
    identical(
      do.call(
        build_national_zero_months,
        c(a, list(head_window = 0L, consistency_scope = "threshold"))
      ),
      do.call(old_buildNationalZeroMonths, a)
    ),
    sprintf(
      "call %d: head_window 0 + scope threshold == frozen pre-registry function",
      i
    )
  )
}
zd6 <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  head_window = 4L,
  edge_detail = TRUE
)
ed <- attr(zd6, "edge_detail")
expect(
  identical(sort(ed$country[ed$head_unclaimed > 0L]), c("LLL", "OOO")),
  "edge_detail: per country-year head/tail unclaimed"
)
attr(zd6, "edge_detail") <- NULL
expect(
  identical(zd6, z6),
  "edge_detail stripped -> identical to the call without it"
)

cat("\n[9] rule registry\n")
lg <- national_zero_rule_params("")
expect(
  identical(lg, national_zero_rule_params("legacy")) &&
    identical(lg$rule, "legacy"),
  "\"\" reads the legacy entry"
)
src01 <- readLines("analysis/01_preprocess_national_zero_months.R")
expect(
  any(grepl("national_zero_rule_params\\(\"legacy\"\\)", src01)),
  "01 reads its option defaults from the registry's legacy entry"
)
legacy_cli <- list(
  zero_day_frac = 0.8,
  window = 4L,
  head_window = 0L,
  db_tolerance = 1.2,
  consistency_scope = "threshold",
  end = "2024-12",
  exclude = "ETH:2011-2014"
)
expect(
  identical(lg[names(legacy_cli)], legacy_cli),
  "legacy = the pre-registry 01 CLI defaults (0.8 / 4 / 0 / 1.2 / threshold / 2024-12 / ETH:2011-2014)"
)
r2 <- national_zero_rule_params("r2")
expect(
  identical(r2$head_window, 4L) &&
    identical(r2$consistency_scope, "all") &&
    identical(r2$exclude, "ETH:2010-2016") &&
    identical(r2$end, "2024-12"),
  "r2 = head 4, scope all, ETH:2010-2016 only (D9), end 2024-12"
)
expect(
  nrow(parse_national_zero_exclude(r2$exclude)) == 7L,
  "the r2 exclusion list parses (ETH 2010-2016)"
)
r3 <- national_zero_rule_params("r3")
expect(
  isTRUE(r3$head_missing_quiet) &&
    identical(lg$head_missing_quiet, FALSE) &&
    identical(r2$head_missing_quiet, FALSE) &&
    identical(
      r3[setdiff(names(r3), c("rule", "head_missing_quiet"))],
      r2[setdiff(names(r2), c("rule", "head_missing_quiet"))]
    ),
  "r3 = r2 + head_missing_quiet (legacy and r2: false, the key absent)"
)
old_digest <- function(params) {
  # national_zero_rule_digest() before the optional key
  p <- params[setdiff(names(params), c("rule", "head_missing_quiet"))]
  digest::digest(p[order(names(p))], algo = "xxhash64")
}
expect(
  identical(national_zero_rule_digest(lg), old_digest(lg)) &&
    identical(national_zero_rule_digest(r2), old_digest(r2)) &&
    !identical(national_zero_rule_digest(r3), national_zero_rule_digest(r2)),
  "legacy / r2 digests unchanged by the optional key (existing stamps stay valid); r3 differs"
)
expect(
  grepl("is not in", errmsg(national_zero_rule_params("zz9"))),
  "unknown rule refused"
)
expect(
  !identical(national_zero_rule_digest(lg), national_zero_rule_digest(r2)) &&
    identical(
      national_zero_rule_digest(lg),
      national_zero_rule_digest(national_zero_rule_params("legacy"))
    ),
  "rule digest: stable, distinguishes rules"
)
expect(
  identical(
    as.character(make_national_zero_months_file(60L, "r2")),
    paste0(SECURE_OUTBREAK_DIR, "national_zero_months_t60_r2.rds")
  ) &&
    identical(
      as.character(make_national_zero_months_file(0L, "r2")),
      paste0(SECURE_OUTBREAK_DIR, "national_zero_months_r2.rds")
    ) &&
    identical(
      as.character(make_national_zero_months_file(60L, "")),
      as.character(make_national_zero_months_file(60L))
    ) &&
    identical(
      as.character(make_national_zero_months_file(60L, "legacy")),
      as.character(make_national_zero_months_file(60L))
    ),
  "make_national_zero_months_file(threshold, rule): _t60_r2 / _r2; legacy unchanged"
)
expect(
  identical(
    as.character(national_zero_consistency_path(60L)),
    sub(
      "\\.rds$",
      "_db_who_consistency.csv",
      as.character(make_national_zero_months_file(60L))
    )
  ) &&
    identical(
      as.character(national_zero_consistency_path(0L)),
      sub(
        "\\.rds$",
        "_db_who_consistency.csv",
        as.character(make_national_zero_months_file(0L))
      )
    ),
  "db_who_consistency sidecar: the rule-free stem item 10's report_who_scale reads"
)
expect(
  grepl("partial", errmsg(.national_zero_end_guard(d("2026-03-01"), 2026L))) &&
    identical(errmsg(.national_zero_end_guard(d("2024-12-01"), 2026L)), ""),
  "end guard: an end inside the dashboard's partial year is refused"
)

cat("\n[10] national_zero_audit()\n")
zd6 <- build_national_zero_months(
  nat0,
  sub_pos[0, ],
  act4[0, ],
  who6,
  start = d("2020-01-01"),
  end = d("2020-12-01"),
  zero_threshold = 60L,
  head_window = 4L,
  edge_detail = TRUE
)
au <- national_zero_audit(zd6, who6, consistency = NULL, years = 2020L)
expect(
  identical(
    names(au),
    c(
      "country",
      "year",
      "who_cases_max",
      "who_exact_zero",
      "db_cases_deep",
      "db_who_ratio",
      "db_exceeds_who",
      "months_national_and_who",
      "months_national_who_missing",
      "months_who_no_national",
      "months_who_below_threshold",
      "months_claimed",
      "head_unclaimed",
      "tail_unclaimed",
      "auto_excluded_exact",
      "auto_excluded_threshold",
      "static_excluded"
    )
  ),
  "D32 audit columns"
)
auL <- au[au$country == "LLL", ]
expect(
  nrow(au) == 4L && auL$months_claimed == 8L && auL$head_unclaimed == 4L,
  "LLL 2020: 8 months claimed, 4 head months unclaimed"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
