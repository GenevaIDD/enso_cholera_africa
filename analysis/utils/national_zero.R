# national_zero.R — national-zero months: country-months in which a whole
#   country had no suspected cholera (national weekly series + WHO annual)
#
# Provides:
#   build_national_zero_months()   apply the rule on a country x month grid
#   national_zero_rule_params()    one named rule of the rule registry
#   national_zero_rule_digest()    digest of a rule's parameters
#   national_month_status(), subnational_positive_months(),
#   active_episode_months()        the rule's inputs
#   parse_national_zero_exclude()  "ISO3:Y0-Y1,..." -> excluded country-years
#   db_who_consistency()           database-vs-WHO consistency check
#   national_zero_audit()          per-country audit of the output
#
# Details:
#   Built by 01_preprocess_national_zero_months.R; read by
#   coverage.R::build_coverage_mask(), which adds the months, expanded to
#   every gid of the country, to the surveillance-coverage mask as observed
#   zeros. WHO zeros come from who_annual.R::build_who_annual_zero_mask().
#
# See: build_national_zero_months() for the rule

# Rule registry: analysis/national_zero_rules.yml (named parameter sets;
# `legacy` = the old 01 CLI defaults). Options beyond legacy, off there:
#   * head_window: the first `head_window` months of a below-threshold year
#     count only if the previous year is also below the threshold;
#   * consistency_scope = "all": db_who_consistency() also drops WHO exact
#     zeros that the database contradicts;
#   * head_missing_quiet (r3; absent = false): under the head guard, a
#     previous year with no WHO figure counts as quiet. The tail guard is
#     unchanged.
# At the legacy values every output and attribute is unchanged, so legacy and
# r2 keep their digests.
NATIONAL_ZERO_RULES_PATH <- "analysis/national_zero_rules.yml"
.NATIONAL_ZERO_RULE_KEYS <- c(
  "zero_day_frac",
  "window",
  "head_window",
  "db_tolerance",
  "consistency_scope",
  "end",
  "exclude"
)
.NATIONAL_ZERO_RULE_OPTIONAL_KEYS <- c("head_missing_quiet") # absent = false
NATIONAL_ZERO_CONSISTENCY_SCOPES <- c("threshold", "all")

#' national_zero_rule_params
#'
#' One named parameter set of the national-zero rule registry.
#'
#' @param rule registry name; "" (and NULL) = "legacy".
#' @param path the registry YAML.
#' @return list(rule, zero_day_frac, window, head_window, db_tolerance,
#'   consistency_scope, end, exclude, head_missing_quiet) -- `exclude` as the comma-joined
#'   "ISO3:Y0-Y1,..." string 01 and parse_national_zero_exclude() take, `rule` the
#'   name ("legacy" for "").
national_zero_rule_params <- function(
  rule = "",
  path = NATIONAL_ZERO_RULES_PATH
) {
  r <- rule %||% ""
  if (!nzchar(r)) {
    r <- "legacy"
  }
  if (!file.exists(path)) {
    stop(
      "national_zero_rule_params: rule registry not found: ",
      path,
      call. = FALSE
    )
  }
  reg <- yaml::read_yaml(path)
  if (!r %in% names(reg)) {
    stop(
      "national-zero rule '",
      r,
      "' is not in ",
      path,
      " (has: ",
      paste(names(reg), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  e <- reg[[r]]
  miss <- setdiff(.NATIONAL_ZERO_RULE_KEYS, names(e))
  extra <- setdiff(
    names(e),
    c(.NATIONAL_ZERO_RULE_KEYS, .NATIONAL_ZERO_RULE_OPTIONAL_KEYS)
  )
  if (length(miss) || length(extra)) {
    stop(
      "national-zero rule '",
      r,
      "' in ",
      path,
      ": ",
      if (length(miss)) {
        paste0("missing ", paste(miss, collapse = ", "), "; ")
      } else {
        ""
      },
      if (length(extra)) {
        paste0("unknown ", paste(extra, collapse = ", "))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  excl <- vapply(
    e$exclude %||% list(),
    function(x) as.character(x$entry %||% x),
    character(1)
  )
  parse_national_zero_exclude(paste(excl, collapse = ",")) # validates the entries
  if (!e$consistency_scope %in% NATIONAL_ZERO_CONSISTENCY_SCOPES) {
    stop(
      "national-zero rule '",
      r,
      "': consistency_scope must be one of {",
      paste(NATIONAL_ZERO_CONSISTENCY_SCOPES, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
  hmq <- e[["head_missing_quiet"]]
  if (is.null(hmq)) {
    hmq <- FALSE
  } # optional key, absent = false
  if (!is.logical(hmq) || length(hmq) != 1L || is.na(hmq)) {
    stop(
      "national-zero rule '",
      r,
      "': head_missing_quiet must be a YAML boolean.",
      call. = FALSE
    )
  }
  list(
    rule = r,
    zero_day_frac = as.numeric(e$zero_day_frac),
    window = as.integer(e$window),
    head_window = as.integer(e$head_window),
    db_tolerance = as.numeric(e$db_tolerance),
    consistency_scope = as.character(e$consistency_scope),
    end = as.character(e$end),
    exclude = paste(excl, collapse = ","),
    head_missing_quiet = hmq
  )
}

#' Digest of a rule's parameters (stamped by 01, compared by 03 and the deploy).
#' An optional key at its default (head_missing_quiet false) is left out, so the
#' rules written before it (legacy, r2) keep their digests.
national_zero_rule_digest <- function(params) {
  p <- params[setdiff(names(params), "rule")]
  if (!isTRUE(p[["head_missing_quiet"]])) {
    p[["head_missing_quiet"]] <- NULL
  }
  digest::digest(p[order(names(p))], algo = "xxhash64")
}

# Refuse an `end` inside a partial WHO year: the dashboard's latest year is
# partial (who_annual.R), so a below-threshold total there means nothing.
.national_zero_end_guard <- function(end, partial_year) {
  if (is.null(partial_year) || !length(partial_year) || is.na(partial_year)) {
    return(invisible(TRUE))
  }
  if (as.integer(lubridate::year(as.Date(end))) >= as.integer(partial_year)) {
    stop(
      "national-zero rule: end = ",
      format(as.Date(end), "%Y-%m"),
      " lies in ",
      partial_year,
      ", the WHO dashboard's partial (latest) year.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Edge guard of the threshold branch, one helper for both ends of the year so
# the head guard is the tail guard's mirror by construction. `side` "tail": the
# last `width` months, claimable only when `adjacent_below` (the NEXT year is
# below); "head": the first `width` months, claimable only when the PREVIOUS
# year is below. Returns list(in_edge, ok) over `month`.
.year_edge_claim <- function(
  month,
  width,
  side = c("tail", "head"),
  adjacent_below
) {
  side <- match.arg(side)
  m <- as.integer(lubridate::month(month))
  in_edge <- if (side == "tail") m > 12L - width else m <= width
  list(in_edge = in_edge, ok = !in_edge | adjacent_below)
}

#' .month_span
#'
#' Vectorised `map2(from, to, ~ seq(.x, .y, by = "month"))` + unnest: the row
#' index of each expanded month and the month itself, in the same order.
#' @param from,to first-of-month Dates, `to >= from`.
#' @return list(idx = integer row index into from/to, month = Date).
.month_span <- function(from, to) {
  m0 <- as.integer(lubridate::year(from)) *
    12L +
    as.integer(lubridate::month(from)) -
    1L
  n <- as.integer(lubridate::year(to)) *
    12L +
    as.integer(lubridate::month(to)) -
    1L -
    m0 +
    1L
  if (any(n < 1L, na.rm = TRUE) || anyNA(n)) {
    stop(
      ".month_span: every `to` must be a month on or after its `from`.",
      call. = FALSE
    )
  }
  idx <- rep(seq_along(from), n)
  mm <- m0[idx] + sequence(n) - 1L
  list(idx = idx, month = lubridate::make_date(mm %/% 12L, mm %% 12L + 1L, 1L))
}

#' national_month_status
#'
#' @param ts national rows of the pre-extraction weekly series
#'   (`spatial_scale == "country"`; columns country, TL, TR, sCh).
#' @param zero_day_frac a month is a national zero when no overlapping week
#'   has cases and zero weeks cover at least this share of its days.
#' @return tibble(country, month, nat_pos, nat_zero); one row per month any
#'   national week overlaps.
national_month_status <- function(ts, zero_day_frac = 0.8) {
  ts <- ts %>%
    dplyr::filter(!is.na(sCh)) %>%
    dplyr::mutate(TL = as.Date(TL), TR = as.Date(TR))
  if (!nrow(ts)) {
    return(tibble::tibble(
      country = character(),
      month = as.Date(character()),
      nat_pos = logical(),
      nat_zero = logical()
    ))
  }
  sp <- .month_span(
    lubridate::floor_date(ts$TL, "month"),
    lubridate::floor_date(ts$TR, "month")
  )
  ts[sp$idx, , drop = FALSE] %>%
    dplyr::mutate(month = sp$month) %>%
    dplyr::mutate(
      ov = pmax(
        0,
        as.numeric(
          pmin(TR, lubridate::ceiling_date(month, "month") - 1) -
            pmax(TL, month)
        ) +
          1
      )
    ) %>%
    dplyr::group_by(country, month) %>%
    dplyr::summarise(
      nat_pos = any(sCh > 0),
      zero_days = sum(ov[sCh == 0]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      nat_zero = !nat_pos &
        zero_days >= zero_day_frac * lubridate::days_in_month(month)
    ) %>%
    dplyr::select(country, month, nat_pos, nat_zero)
}

#' subnational_positive_months
#'
#' @param ts sub-national rows of the pre-extraction weekly series.
#' @return tibble(country, month) of months with any sub-national report of
#'   cases: every month the reporting week overlaps, as national_month_status()
#'   expands national weeks (a week straddling a month boundary is positive in
#'   both months, so the second one cannot be declared a national zero).
subnational_positive_months <- function(ts) {
  ts <- ts %>%
    dplyr::filter(!is.na(sCh), sCh > 0) %>%
    dplyr::mutate(TL = as.Date(TL), TR = as.Date(dplyr::coalesce(TR, TL)))
  sp <- .month_span(
    lubridate::floor_date(ts$TL, "month"),
    lubridate::floor_date(pmax(ts$TL, ts$TR), "month")
  )
  tibble::tibble(country = ts$country[sp$idx], month = sp$month) %>%
    dplyr::distinct()
}

#' active_episode_months
#'
#' @param ob_meta database outbreaks (country, TL, TR), e.g.
#'   readRDS(make_outbreak_metadata_file()) -- independent of the episode key.
#' @return tibble(country, month) of months any outbreak of the country spans.
active_episode_months <- function(ob_meta) {
  ob_meta <- ob_meta %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(TL), !is.na(TR)) %>%
    dplyr::distinct(country, TL = as.Date(TL), TR = as.Date(TR))
  if (!nrow(ob_meta)) {
    return(tibble::tibble(country = character(), month = as.Date(character())))
  }
  sp <- .month_span(
    lubridate::floor_date(ob_meta$TL, "month"),
    lubridate::floor_date(ob_meta$TR, "month")
  )
  tibble::tibble(country = ob_meta$country[sp$idx], month = sp$month) %>%
    dplyr::distinct(country, month)
}

#' parse_national_zero_exclude
#'
#' @param x "ISO3:Y0-Y1[,ISO3:Y0-Y1...]" (a single year "ISO3:Y" also works);
#'   "" for none.
#' @return tibble(country, year).
parse_national_zero_exclude <- function(x) {
  x <- trimws(strsplit(x %||% "", ",")[[1]])
  x <- x[nzchar(x)]
  if (!length(x)) {
    return(tibble::tibble(country = character(), year = integer()))
  }
  ok <- grepl("^[A-Z]{3}:[0-9]{4}(-[0-9]{4})?$", x)
  if (!all(ok)) {
    stop(
      "exclude: expected ISO3:Y0-Y1 entries, got: ",
      paste(x[!ok], collapse = ", "),
      call. = FALSE
    )
  }
  purrr::map_dfr(x, function(e) {
    p <- strsplit(e, "[:-]")[[1]]
    y0 <- as.integer(p[2])
    y1 <- as.integer(p[length(p)])
    tibble::tibble(country = p[1], year = seq(y0, y1))
  })
}

# Rule (one country-month):
#   ((national 0 AND no WHO > 0) OR (WHO 0 AND no national > 0))
#     AND no sub-national cases AND no active database outbreak that month
#     AND not excluded (default ETH 2011-2014, AWD-era WHO zeros)
# A missing WHO row counts as "no WHO > 0"; a month without national reports
# as "no national > 0".
#
# Threshold variant (zero_threshold > 0): a WHO country-year under the
# threshold (max over WHO sources) cannot hold a >= 50-case onset, so its
# months become national-zero months (source "who_below_threshold") with the
# same drops. Its last `window` months count only if the next year is also
# below the threshold. Country-years whose database cases exceed the WHO
# total (db_who_consistency()) are excluded. zero_threshold = 0 gives the
# exact-zero rule.
#' build_national_zero_months
#'
#' Apply the national-zero rule (see above) on a country x month grid.
#'
#' @param nat national_month_status() output.
#' @param sub_pos subnational_positive_months() output.
#' @param active active_episode_months() output.
#' @param who_mask build_who_annual_zero_mask() output (country, year, is_zero).
#' @param start,end first / last month (Dates; floored to the month).
#' @param exclude tibble(country, year) of excluded country-years.
#' @param zero_threshold WHO country-years with fewer cases than this (max over
#'   sources, `who_mask$suspected_cases_max`) also count; 0 = exact zeros only.
#' @param window months at the end of a below-threshold year claimed only when
#'   the next year is also below the threshold (episode length assumed).
#' @param auto_exclude db_who_consistency() output (or NULL): its `excluded`
#'   country-years are removed from the threshold branch.
#' @param head_window threshold branch: the first months of a below-threshold
#'   year claimed only when the previous year is also below (0 = off, legacy).
#' @param head_missing_quiet head guard: a previous year with no WHO figure
#'   counts as below the threshold (TRUE, rule r3) or not (FALSE, legacy / r2).
#' @param consistency_scope "threshold" (legacy) | "all": `auto_exclude` also
#'   removes rule-zero months resting on a WHO exact zero (who_zero and no
#'   national zero).
#' @param edge_detail TRUE adds attribute `edge_detail` (country, year,
#'   head_unclaimed, tail_unclaimed) for national_zero_audit(); strip it before
#'   saving (the saved file must equal the pre-registry output).
#' @return tibble(country, month, source) with source one of
#'   "national_and_who" | "national_who_missing" | "who_no_national" |
#'   "who_below_threshold". Attribute `dropped`: named integer, rule-zero
#'   country-months removed for an active outbreak, sub-national cases, the
#'   exclusion list, or (threshold branch) the database consistency check;
#'   attribute `tail_unclaimed`: below-threshold months left out because the
#'   next year is not below the threshold; `head_unclaimed` (only when
#'   head_window > 0) the same at the head; `dropped` gains
#'   `auto_excluded_exact` only under consistency_scope = "all".
build_national_zero_months <- function(
  nat,
  sub_pos,
  active,
  who_mask,
  start,
  end,
  exclude = tibble::tibble(country = character(), year = integer()),
  zero_threshold = 0L,
  window = 4L,
  auto_exclude = NULL,
  head_window = 0L,
  consistency_scope = "threshold",
  edge_detail = FALSE,
  head_missing_quiet = FALSE
) {
  if (
    !is.logical(head_missing_quiet) ||
      length(head_missing_quiet) != 1L ||
      is.na(head_missing_quiet)
  ) {
    stop(
      "build_national_zero_months: head_missing_quiet must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  zero_threshold <- as.integer(zero_threshold)
  window <- as.integer(window)
  head_window <- as.integer(head_window)
  if (is.na(head_window) || head_window < 0L || head_window > 12L) {
    stop(
      "build_national_zero_months: head_window must be in 0..12.",
      call. = FALSE
    )
  }
  if (!consistency_scope %in% NATIONAL_ZERO_CONSISTENCY_SCOPES) {
    stop(
      "build_national_zero_months: consistency_scope must be one of {",
      paste(NATIONAL_ZERO_CONSISTENCY_SCOPES, collapse = ", "),
      "}.",
      call. = FALSE
    )
  }
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
  auto_cy <- if (is.null(auto_exclude)) {
    character(0)
  } else {
    with(
      auto_exclude[auto_exclude$excluded, , drop = FALSE],
      paste(country, year)
    )
  }
  head_unclaimed <- 0L
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
    # ... and, mirrored, the previous year (head guard; head_window > 0). A
    # previous year with no WHO figure is quiet only under head_missing_quiet.
    prv <- who_mask %>%
      dplyr::transmute(
        country,
        year = as.integer(year) + 1L,
        prev_below = suspected_cases_max < zero_threshold
      )
    grid <- grid %>%
      dplyr::left_join(below, by = c("country", "year")) %>%
      dplyr::left_join(nxt, by = c("country", "year")) %>%
      dplyr::left_join(prv, by = c("country", "year")) %>%
      dplyr::mutate(
        below = dplyr::coalesce(below, FALSE),
        next_below = dplyr::coalesce(next_below, FALSE),
        prev_below = dplyr::coalesce(prev_below, head_missing_quiet)
      )
    tl <- .year_edge_claim(grid$month, window, "tail", grid$next_below)
    hd <- .year_edge_claim(grid$month, head_window, "head", grid$prev_below)
    base_below <- grid$below & !grid$nat_pos & !grid$rule
    grid$in_tail <- tl$in_edge
    grid$in_head <- hd$in_edge
    grid$below_rule <- base_below & tl$ok & hd$ok
    tail_unclaimed <- sum(base_below & !tl$ok)
    head_unclaimed <- sum(base_below & tl$ok & !hd$ok)
    grid$tail_unc <- base_below & !tl$ok
    grid$head_unc <- base_below & tl$ok & !hd$ok
    grid$auto <- grid$below_rule & paste(grid$country, grid$year) %in% auto_cy
  } else {
    grid$below_rule <- FALSE
    grid$auto <- FALSE
    grid$tail_unc <- FALSE
    grid$head_unc <- FALSE
  }
  # consistency_scope = "all": a WHO exact zero contradicted by the database
  # (db_who_consistency() excluded) no longer carries rule-zero months resting on
  # it alone; a national zero (nat_zero) still does.
  grid$auto_exact <- if (identical(consistency_scope, "all")) {
    grid$rule &
      grid$who_zero &
      !grid$nat_zero &
      paste(grid$country, grid$year) %in% auto_cy
  } else {
    FALSE
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
  if (identical(consistency_scope, "all")) {
    dropped <- c(
      dropped,
      auto_excluded_exact = sum(r$auto_exact & !r$excl & !r$sub & !r$act)
    )
  }
  out <- r %>%
    dplyr::filter(!sub, !act, !excl, !auto, !auto_exact) %>%
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
  if (zero_threshold > 0L && head_window > 0L) {
    attr(out, "head_unclaimed") <- as.integer(head_unclaimed)
  }
  if (isTRUE(edge_detail)) {
    attr(out, "edge_detail") <- grid %>%
      dplyr::group_by(country, year) %>%
      dplyr::summarise(
        head_unclaimed = sum(head_unc),
        tail_unclaimed = sum(tail_unc),
        auto_excluded_exact = sum(auto_exact & !sub & !act & !excl),
        auto_excluded_threshold = sum(auto & !sub & !act & !excl),
        static_excluded = any(excl),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        head_unclaimed > 0L |
          tail_unclaimed > 0L |
          auto_excluded_exact > 0L |
          auto_excluded_threshold > 0L |
          static_excluded
      )
  }
  out
}

#' national_zero_audit
#'
#' The national-zero audit: one row per WHO country-year in the rule's
#' window -- WHO total (max over sources) and exact-zero flag, the database's
#' deepest-level episode cases and their ratio to WHO, the consistency flag,
#' the national-zero months claimed by source, the head / tail months left
#' unclaimed and the exclusions.
#'
#' @param zero build_national_zero_months(..., edge_detail = TRUE) output.
#' @param who_mask build_who_annual_zero_mask() output.
#' @param consistency db_who_consistency() output (or NULL).
#' @param years the rule window's calendar years.
#' @return tibble(country, year, who_cases_max, who_exact_zero, db_cases_deep,
#'   db_who_ratio, db_exceeds_who, months_national_and_who,
#'   months_national_who_missing, months_who_no_national,
#'   months_who_below_threshold, months_claimed, head_unclaimed,
#'   tail_unclaimed, auto_excluded_exact, auto_excluded_threshold,
#'   static_excluded).
national_zero_audit <- function(zero, who_mask, consistency = NULL, years) {
  srcs <- c(
    "national_and_who",
    "national_who_missing",
    "who_no_national",
    "who_below_threshold"
  )
  cnt <- zero %>%
    dplyr::mutate(year = as.integer(lubridate::year(month))) %>%
    dplyr::count(country, year, source) %>%
    tidyr::pivot_wider(names_from = source, values_from = n, values_fill = 0L)
  for (s in setdiff(srcs, names(cnt))) {
    cnt[[s]] <- integer(nrow(cnt))
  }
  cnt <- cnt %>%
    dplyr::rename_with(~ paste0("months_", .x), dplyr::all_of(srcs))
  edge <- attr(zero, "edge_detail") %||%
    tibble::tibble(
      country = character(),
      year = integer(),
      head_unclaimed = integer(),
      tail_unclaimed = integer(),
      auto_excluded_exact = integer(),
      auto_excluded_threshold = integer(),
      static_excluded = logical()
    )
  cons <- (consistency %||%
    tibble::tibble(
      country = character(),
      year = integer(),
      db_deep = numeric(),
      ratio = numeric(),
      excluded = logical()
    )) %>%
    dplyr::transmute(
      country,
      year = as.integer(year),
      db_cases_deep = db_deep,
      db_who_ratio = ratio,
      db_exceeds_who = excluded
    )
  who_mask %>%
    dplyr::filter(as.integer(year) %in% as.integer(years)) %>%
    dplyr::transmute(
      country,
      year = as.integer(year),
      who_cases_max = suspected_cases_max,
      who_exact_zero = is_zero
    ) %>%
    dplyr::left_join(cons, by = c("country", "year")) %>%
    dplyr::left_join(cnt, by = c("country", "year")) %>%
    dplyr::left_join(edge, by = c("country", "year")) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          dplyr::starts_with("months_"),
          head_unclaimed,
          tail_unclaimed,
          auto_excluded_exact,
          auto_excluded_threshold
        ),
        ~ as.integer(dplyr::coalesce(.x, 0L))
      ),
      static_excluded = dplyr::coalesce(static_excluded, FALSE),
      db_exceeds_who = dplyr::coalesce(db_exceeds_who, FALSE),
      months_claimed = months_national_and_who +
        months_national_who_missing +
        months_who_no_national +
        months_who_below_threshold
    ) %>%
    dplyr::select(
      country,
      year,
      who_cases_max,
      who_exact_zero,
      db_cases_deep,
      db_who_ratio,
      db_exceeds_who,
      months_national_and_who,
      months_national_who_missing,
      months_who_no_national,
      months_who_below_threshold,
      months_claimed,
      head_unclaimed,
      tail_unclaimed,
      auto_excluded_exact,
      auto_excluded_threshold,
      static_excluded
    ) %>%
    dplyr::arrange(country, year)
}

#' db_who_consistency
#'
#' Is the WHO annual total an upper bound on the database's own episode cases?
#' Each episode's `suspected_cases` is apportioned to calendar years by day
#' overlap of TL..TR; per country-year only episodes at the deepest admin level
#' present are summed (nested levels would double count). A country-year whose
#' database cases exceed `tolerance` x the WHO total (`suspected_cases_max`) is
#' `excluded`: its WHO figure is an under-report, so a below-threshold reading
#' of it cannot be trusted (e.g. SSD 2023, database 11x WHO).
#'
#' @param ob_meta database outbreaks (country, TL, TR, suspected_cases,
#'   admin_level), e.g. readRDS(make_outbreak_metadata_file()).
#' @param who_mask build_who_annual_zero_mask() output.
#' @param tolerance excluded when db_deep > tolerance * who_max.
#' @return tibble(country, year, db_deep, who_max, ratio, excluded), one row per
#'   country-year present in both.
db_who_consistency <- function(ob_meta, who_mask, tolerance = 1.2) {
  if (!is.finite(tolerance) || tolerance <= 0) {
    stop("db_who_consistency: tolerance must be > 0.", call. = FALSE)
  }
  ep <- ob_meta %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(TL), !is.na(TR), !is.na(suspected_cases)) %>%
    dplyr::transmute(
      country = as.character(country),
      TL = as.Date(TL),
      TR = as.Date(TR),
      cases = as.numeric(suspected_cases),
      lvl = as.integer(sub("^ADM", "", as.character(admin_level)))
    )
  empty <- tibble::tibble(
    country = character(),
    year = integer(),
    db_deep = numeric(),
    who_max = numeric(),
    ratio = numeric(),
    excluded = logical()
  )
  if (!nrow(ep)) {
    return(empty)
  }
  ep$days <- as.numeric(ep$TR - ep$TL) + 1
  yrs <- ep %>%
    dplyr::mutate(
      year = purrr::map2(
        as.integer(lubridate::year(TL)),
        as.integer(lubridate::year(TR)),
        seq
      )
    ) %>%
    tidyr::unnest(year) %>%
    dplyr::mutate(
      y0 = pmax(TL, as.Date(sprintf("%d-01-01", year))),
      y1 = pmin(TR, as.Date(sprintf("%d-12-31", year))),
      cases_y = cases * (as.numeric(y1 - y0) + 1) / days
    )
  db <- yrs %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(db_deep = sum(cases_y[lvl == max(lvl)]), .groups = "drop")
  out <- db %>%
    dplyr::inner_join(
      who_mask %>%
        dplyr::transmute(
          country,
          year = as.integer(year),
          who_max = suspected_cases_max
        ),
      by = c("country", "year")
    ) %>%
    dplyr::mutate(
      ratio = db_deep / pmax(who_max, 1),
      excluded = db_deep > tolerance * who_max
    ) %>%
    dplyr::arrange(country, year)
  if (!nrow(out)) empty else out
}
