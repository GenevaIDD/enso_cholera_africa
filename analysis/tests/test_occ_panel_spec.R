# test_occ_panel_spec.R — HARD: one occurrence-panel spec, one token, one
#   matcher.
#
# Usage:
#   Rscript analysis/tests/test_occ_panel_spec.R
#
# Details:
#   paths.R::occ_panel_spec_of() holds every episodes.* key that changes the
#   panel's ROWS: national_zero_{fill,threshold,rule}, at_risk_start,
#   adm1_overlap. .occ_panel_variant_token() builds the suffix in the fixed
#   order "_ep-gs" -> "_nz<thr>[-<rule>]" -> "_ae" -> "_a1c";
#   .occ_panel_variant_file_match() is an EXACT whole-suffix filter. At the
#   defaults every name, block and generated config is byte-identical to the
#   frozen reference copies in this file. Pure synthetic; no data/.

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_parse.R")

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

# ---- frozen reference code (c032a12) ---------------------------------------
old_nz_token <- function(opt) {
  thr <- national_zero_threshold_of(opt)
  if (!national_zero_fill_of(opt)) {
    return("")
  }
  if (thr > 0L) paste0("_nz", thr) else "_nz"
}
old_file_match <- function(
  files,
  national_zero_fill = FALSE,
  national_zero_threshold = 0L
) {
  b <- basename(files)
  thr <- .national_zero_threshold_check(national_zero_threshold)
  tok <- regmatches(b, regexpr("_ep-gs_nz[0-9]*(?=_|\\.rds$)", b, perl = TRUE))
  have <- rep(NA_character_, length(b))
  have[grepl("_ep-gs_nz[0-9]*(_|\\.rds$)", b, perl = TRUE)] <- tok
  if (!isTRUE(national_zero_fill)) {
    return(is.na(have))
  }
  want <- if (thr > 0L) paste0("_ep-gs_nz", thr) else "_ep-gs_nz"
  !is.na(have) & have == want
}
old_apply <- function(
  cfg,
  episode_key,
  label = "config",
  national_zero_fill = NULL,
  national_zero_threshold = NULL
) {
  if (is.null(episode_key) || !nzchar(episode_key)) {
    return(cfg)
  }
  tpl_nz <- (cfg$episodes %||% list())$national_zero_fill
  nzf <- isTRUE(national_zero_fill %||% tpl_nz)
  tpl_nzt <- (cfg$episodes %||% list())$national_zero_threshold
  nzt <- as.integer(national_zero_threshold %||% tpl_nzt %||% 0L)
  if (!identical(episode_key, "lp_season") || !is.null(cfg$episodes)) {
    cfg$episodes <- c(
      list(key = episode_key),
      if (nzf) list(national_zero_fill = TRUE),
      if (nzt > 0L) list(national_zero_threshold = nzt)
    )
  }
  cfg
}

ep <- function(...) list(baseline = "1990_2019", episodes = list(...))

# ---- [1] token ---------------------------------------------------------------
cat("\n[1] variant token\n")
legacy <- list(
  ep(key = "lp_season"),
  ep(key = "gid_season"),
  ep(key = "gid_season", national_zero_fill = TRUE),
  ep(
    key = "gid_season",
    national_zero_fill = TRUE,
    national_zero_threshold = 60L
  ),
  ep(
    key = "gid_season",
    national_zero_fill = FALSE,
    national_zero_threshold = 0L
  ),
  list(baseline = "1990_2019")
)
for (o in legacy) {
  expect(
    identical(
      .occ_panel_variant_token(o),
      paste0(.episode_key_path_token(o), old_nz_token(o))
    ),
    sprintf(
      "legacy spec %s: token == pre-P1 '_ep-gs' + '_nz' composition",
      deparse(o$episodes, width.cutoff = 500L)
    )
  )
}
expect(
  identical(
    .occ_panel_variant_token(ep(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L,
      national_zero_rule = "r2",
      at_risk_start = "after_episode",
      adm1_overlap = "censor_orphan"
    )),
    "_ep-gs_nz60-r2_ae_a1c"
  ),
  "headline token _ep-gs_nz60-r2_ae_a1c (fixed order)"
)
expect(
  identical(
    .occ_panel_variant_token(ep(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_rule = "r2"
    )),
    "_ep-gs_nz-r2"
  ),
  "threshold 0 + rule -> _nz-r2"
)
expect(
  identical(
    .occ_panel_variant_token(ep(
      key = "gid_season",
      at_risk_start = "after_episode"
    )),
    "_ep-gs_ae"
  ),
  "at_risk_start alone -> _ep-gs_ae"
)
expect(
  identical(
    .occ_panel_variant_token(ep(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L,
      national_zero_rule = "legacy"
    )),
    "_ep-gs_nz60"
  ),
  "rule 'legacy' normalises to '' (no token)"
)
expect(
  identical(
    as.character(make_surv_dat_file(ep(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L
    ))),
    "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs_nz60.rds"
  ),
  "make_surv_dat_file: the _nz60 panel name (no rule tokens) is unchanged"
)
expect(
  identical(
    as.character(make_surv_dat_file(ep(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L,
      national_zero_rule = "r2",
      at_risk_start = "after_episode",
      adm1_overlap = "censor_orphan"
    ))),
    "generated_data/outbreaks_survival_data_base-1990_2019_ep-gs_nz60-r2_ae_a1c.rds"
  ),
  "make_surv_dat_file: the headline panel name"
)

# validation
expect(
  grepl(
    "requires episodes.national_zero_fill",
    errmsg(occ_panel_spec_of(
      ep(key = "gid_season", national_zero_rule = "r2")
    ))
  ),
  "rule without the fill is refused"
)
expect(
  grepl(
    "requires episodes.key = 'gid_season'",
    errmsg(occ_panel_spec_of(
      ep(key = "lp_season", at_risk_start = "after_episode")
    ))
  ),
  "at_risk_start needs gid_season"
)
expect(
  grepl(
    "requires episodes.key = 'gid_season'",
    errmsg(occ_panel_spec_of(
      ep(key = "lp_season", adm1_overlap = "censor_orphan")
    ))
  ),
  "adm1_overlap needs gid_season"
)
expect(
  grepl(
    "is not one of",
    errmsg(occ_panel_spec_of(ep(key = "gid_season", at_risk_start = "x")))
  ),
  "unknown at_risk_start refused"
)
expect(
  grepl(
    "lowercase name",
    errmsg(occ_panel_spec_of(
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_rule = "R-2"
      )
    ))
  ),
  "malformed rule name refused"
)
expect(
  identical(
    errmsg(occ_panel_spec_of(ep(key = "lp_season", national_zero_fill = TRUE))),
    errmsg(national_zero_fill_of(ep(
      key = "lp_season",
      national_zero_fill = TRUE
    )))
  ),
  "fill-without-gid_season message unchanged"
)

# every config in the repo: the survival-panel name is the frozen code's name
cfgs <- list.files(
  "analysis/configs",
  pattern = "\\.ya?ml$",
  recursive = TRUE,
  full.names = TRUE
)
n_ok <- 0L
n_bad <- character()
for (f in cfgs) {
  y <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  if (!is.list(y)) {
    next
  }
  blocks <- list(y$episodes, (y$deploy %||% list())$episodes)
  for (e in blocks) {
    if (is.null(e)) {
      next
    }
    # A block that sets a panel field off its default (the deployed spec) has
    # no frozen-code name to keep; its token is pinned by test_prod_run_batches.R.
    if (
      nzchar(e[["national_zero_rule"]] %||% "") ||
        !identical(e[["at_risk_start"]] %||% "season", "season") ||
        !identical(e[["adm1_overlap"]] %||% "drop", "drop")
    ) {
      next
    }
    o <- list(baseline = y$baseline %||% "1990_2019", episodes = e)
    old <- tryCatch(
      paste0(.episode_key_path_token(o), old_nz_token(o)),
      error = function(e) NA
    )
    if (is.na(old)) {
      next
    }
    if (identical(.occ_panel_variant_token(o), old)) {
      n_ok <- n_ok + 1L
    } else {
      n_bad <- c(n_bad, f)
    }
  }
}
expect(
  length(n_bad) == 0L && n_ok > 0L,
  sprintf(
    "every episodes / deploy.episodes block in analysis/configs (%d) keeps its token%s",
    n_ok,
    if (length(n_bad)) {
      paste0(": ", paste(head(n_bad, 3), collapse = ", "))
    } else {
      ""
    }
  )
)

# ---- [2] exact whole-suffix matcher -----------------------------------------
# A "_nz60" lookup must not match "_nz60_ae_a1c" files nor read "_nz60-r2"
# files as "no zero-fill" (as the frozen lookahead matcher above does).
cat("\n[2] exact matcher\n")
stem <- "outbreaks_survival_brms_cb_fit_base-1990_2019_grp-enso_kl-east_sl3"
toks <- c(
  none = "",
  gs = "_ep-gs",
  nz = "_ep-gs_nz",
  nz60 = "_ep-gs_nz60",
  nzr2 = "_ep-gs_nz-r2",
  nz60r2 = "_ep-gs_nz60-r2",
  nz60ae = "_ep-gs_nz60_ae",
  nz60a1c = "_ep-gs_nz60_a1c",
  nz60aea1c = "_ep-gs_nz60_ae_a1c",
  nz60r2aea1c = "_ep-gs_nz60-r2_ae_a1c",
  ae = "_ep-gs_ae",
  a1c = "_ep-gs_a1c"
)
files <- paste0(stem, toks, ".rds")
specs <- list(
  none = .occ_panel_spec(),
  gs = .occ_panel_spec("gid_season"),
  nz = .occ_panel_spec("gid_season", TRUE),
  nz60 = .occ_panel_spec("gid_season", TRUE, 60L),
  nzr2 = .occ_panel_spec("gid_season", TRUE, 0L, "r2"),
  nz60r2 = .occ_panel_spec("gid_season", TRUE, 60L, "r2"),
  nz60ae = .occ_panel_spec(
    "gid_season",
    TRUE,
    60L,
    at_risk_start = "after_episode"
  ),
  nz60a1c = .occ_panel_spec(
    "gid_season",
    TRUE,
    60L,
    adm1_overlap = "censor_orphan"
  ),
  nz60aea1c = .occ_panel_spec(
    "gid_season",
    TRUE,
    60L,
    at_risk_start = "after_episode",
    adm1_overlap = "censor_orphan"
  ),
  nz60r2aea1c = .occ_panel_spec(
    "gid_season",
    TRUE,
    60L,
    "r2",
    "after_episode",
    "censor_orphan"
  ),
  ae = .occ_panel_spec("gid_season", at_risk_start = "after_episode"),
  a1c = .occ_panel_spec("gid_season", adm1_overlap = "censor_orphan")
)
for (nm in names(specs)) {
  m <- .occ_panel_variant_file_match(files, specs[[nm]])
  # "none" and "gs" share the empty variant: the "_ep-gs" filter is the callers'
  want <- if (nm %in% c("none", "gs")) {
    names(toks) %in% c("none", "gs")
  } else {
    names(toks) == nm
  }
  expect(
    identical(m, want),
    sprintf("spec %-11s matches exactly its own files", nm)
  )
  expect(
    identical(.occ_panel_variant_token(specs[[nm]]), toks[[nm]]) ||
      nm == "none",
    sprintf("spec %-11s token == '%s'", nm, toks[[nm]])
  )
}
# the bug the old matcher had (documented, then refused by the wrapper)
expect(
  old_file_match(paste0(stem, "_ep-gs_nz60_ae_a1c.rds"), TRUE, 60L),
  "pre-P1 matcher: a _nz60 lookup matched _nz60_ae_a1c (the bug)"
)
expect(
  old_file_match(paste0(stem, "_ep-gs_nz60-r2.rds"), FALSE),
  "pre-P1 matcher: a no-fill lookup matched _nz60-r2 (the bug)"
)
expect(
  !.national_zero_fill_file_match(
    paste0(stem, "_ep-gs_nz60_ae_a1c.rds"),
    TRUE,
    60L
  ),
  "wrapper: _nz60 does NOT match _nz60_ae_a1c"
)
expect(
  !.national_zero_fill_file_match(
    paste0(stem, "_ep-gs_nz60-r2.rds"),
    TRUE,
    60L
  ),
  "wrapper: _nz60 does NOT match _nz60-r2"
)
expect(
  !.national_zero_fill_file_match(paste0(stem, "_ep-gs_nz60-r2.rds"), FALSE),
  "wrapper: no-fill does NOT match _nz60-r2"
)
# legacy filenames (no rule / _ae / _a1c token): wrapper == frozen old matcher
legacy_files <- c(
  paste0(
    stem,
    c("", "_ep-gs", "_ep-gs_nz", "_ep-gs_nz60", "_ep-gs_nz50"),
    ".rds"
  ),
  "outbreaks_survival_data_base-1990_2019_ep-gs_nz60.policy.rds",
  "x_ep-gs_nz60_sl3.rds",
  "x_ep-gs_nz.rds",
  "x.rds"
)
for (a in list(
  list(FALSE, 0L),
  list(TRUE, 0L),
  list(TRUE, 60L),
  list(TRUE, 50L)
)) {
  expect(
    identical(
      .national_zero_fill_file_match(legacy_files, a[[1]], a[[2]]),
      old_file_match(legacy_files, a[[1]], a[[2]])
    ),
    sprintf(
      "legacy filenames: wrapper == pre-P1 matcher (fill=%s, thr=%d)",
      a[[1]],
      a[[2]]
    )
  )
}
expect(
  identical(
    .occ_panel_variant_file_match(
      files,
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L
      )
    ),
    .occ_panel_variant_file_match(files, specs$nz60)
  ),
  "matcher accepts an opt list as well as a spec"
)

# ---- [3] .episodes_block() -----------------------------------------------------
cat("\n[3] .episodes_block()\n")
for (s in list(specs$none, specs$gs, specs$nz, specs$nz60)) {
  nzf <- s$national_zero_fill
  nzt <- s$national_zero_threshold
  ek <- s$key
  expect(
    identical(
      .episodes_block(s),
      c(
        list(key = ek),
        if (nzf) list(national_zero_fill = TRUE),
        if (nzt > 0L) list(national_zero_threshold = nzt)
      )
    ),
    sprintf(
      "%s: default block == provenance/hindcast hand-built list",
      .occ_panel_variant_token(s)
    )
  )
  expect(
    identical(
      .episodes_block(s, always = "national_zero_fill"),
      c(
        list(key = ek, national_zero_fill = nzf),
        if (nzt > 0L) list(national_zero_threshold = nzt)
      )
    ),
    sprintf(
      "%s: always=fill block == deploy_core hand-built list",
      .occ_panel_variant_token(s)
    )
  )
  expect(
    identical(
      .episodes_block(
        s,
        always = c("national_zero_fill", "national_zero_threshold")
      ),
      list(key = ek, national_zero_fill = nzf, national_zero_threshold = nzt)
    ),
    sprintf(
      "%s: always=fill,threshold block == the full triple",
      .occ_panel_variant_token(s)
    )
  )
}
expect(
  identical(
    .episodes_block(specs$nz60r2aea1c),
    list(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L,
      national_zero_rule = "r2",
      at_risk_start = "after_episode",
      adm1_overlap = "censor_orphan"
    )
  ),
  "non-default fields are written"
)
expect(
  identical(
    occ_panel_spec_of(list(episodes = .episodes_block(specs$nz60r2aea1c))),
    specs$nz60r2aea1c
  ),
  "block -> spec round trip"
)

# ---- [4] assert_occ_panel_variant() ------------------------------------------------
cat("\n[4] assert_occ_panel_variant()\n")
td <- tempfile("occpanel")
dir.create(td)
p <- file.path(td, "outbreaks_survival_data_base-1990_2019_ep-gs_nz60.rds")
saveRDS(1, p)
persist_prep_policy_stamp(
  p,
  list(national_zero_fill = TRUE, national_zero_threshold = 60L)
)
o60 <- ep(
  key = "gid_season",
  national_zero_fill = TRUE,
  national_zero_threshold = 60L
)
expect(
  identical(errmsg(assert_occ_panel_variant(o60, p)), ""),
  "a stamp without the new fields reads as the defaults"
)
expect(
  identical(errmsg(assert_national_zero_fill(o60, p)), ""),
  "assert_national_zero_fill wrapper passes"
)
expect(
  grepl(
    "national_zero_rule = 'r2'",
    errmsg(assert_occ_panel_variant(
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L,
        national_zero_rule = "r2"
      ),
      p
    ))
  ),
  "rule asked, stamp has none -> refused"
)
expect(
  grepl(
    "at_risk_start = 'after_episode'",
    errmsg(assert_occ_panel_variant(
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L,
        at_risk_start = "after_episode"
      ),
      p
    ))
  ),
  "at_risk_start asked, stamp default -> refused"
)
persist_prep_policy_stamp(
  p,
  list(
    national_zero_fill = TRUE,
    national_zero_threshold = 60L,
    national_zero_rule = "r2",
    at_risk_start = "after_episode",
    adm1_overlap = "censor_orphan"
  )
)
expect(
  identical(
    errmsg(assert_occ_panel_variant(
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L,
        national_zero_rule = "r2",
        at_risk_start = "after_episode",
        adm1_overlap = "censor_orphan"
      ),
      p
    )),
    ""
  ),
  "matching headline stamp passes"
)
expect(
  grepl(
    "adm1_overlap = 'drop'",
    errmsg(assert_occ_panel_variant(
      ep(
        key = "gid_season",
        national_zero_fill = TRUE,
        national_zero_threshold = 60L,
        national_zero_rule = "r2",
        at_risk_start = "after_episode"
      ),
      p
    ))
  ),
  "default asked, stamp censor_orphan -> refused"
)
expect(
  grepl(
    "national_zero_threshold = 0",
    errmsg(assert_occ_panel_variant(
      ep(key = "gid_season", national_zero_fill = TRUE),
      p
    ))
  ),
  "threshold message unchanged"
)

# ---- [5] apply_episode_key_to_sim_config(): byte-identical at the defaults ----------
cat("\n[5] apply_episode_key_to_sim_config()\n")
tpls <- list(
  list(simulator = list(run_tag = "a")),
  list(episodes = list(key = "lp_season"), simulator = list()),
  list(episodes = list(key = "gid_season", national_zero_fill = TRUE)),
  list(
    episodes = list(
      key = "gid_season",
      national_zero_fill = TRUE,
      national_zero_threshold = 60L
    )
  )
)
for (tp in tpls) {
  for (k in c("", "lp_season", "gid_season")) {
    if (isTRUE(tp$episodes$national_zero_fill) && !identical(k, "gid_season")) {
      next
    }
    if (
      !is.null(tp$episodes$key) && nzchar(k) && !identical(tp$episodes$key, k)
    ) {
      next
    }
    new <- apply_episode_key_to_sim_config(tp, k)
    old <- old_apply(tp, k)
    expect(
      identical(yaml::as.yaml(new), yaml::as.yaml(old)),
      sprintf(
        "template %s, key '%s': YAML identical to pre-P1",
        deparse(tp$episodes),
        k
      )
    )
    if (identical(k, "gid_season")) {
      sp <- .occ_panel_spec_from_episodes(utils::modifyList(
        tp$episodes %||% list(),
        list(key = k)
      ))
      expect(
        identical(
          yaml::as.yaml(apply_episode_key_to_sim_config(tp, k, occ_panel = sp)),
          yaml::as.yaml(old_apply(
            tp,
            k,
            national_zero_fill = sp$national_zero_fill,
            national_zero_threshold = sp$national_zero_threshold
          ))
        ),
        sprintf(
          "template %s: occ_panel at its own values == pre-P1 pair",
          deparse(tp$episodes)
        )
      )
    }
  }
}
sp19 <- specs$nz60r2aea1c
cfg19 <- apply_episode_key_to_sim_config(
  list(simulator = list()),
  "gid_season",
  occ_panel = sp19
)
expect(
  identical(cfg19$episodes, .episodes_block(sp19)),
  "headline spec written into a config"
)
expect(
  identical(
    apply_episode_key_to_sim_config(cfg19, "gid_season")$episodes,
    cfg19$episodes
  ),
  "template's own item 11/12 fields kept when no occ_panel is given"
)
expect(
  grepl(
    "national_zero_rule",
    errmsg(apply_episode_key_to_sim_config(
      cfg19,
      "gid_season",
      occ_panel = specs$nz60
    ))
  ),
  "template rule vs deploy default -> refused"
)

cat("\n[6] 03 overwrite guard (tier-2 P6)\n")
td7 <- tempfile("guard")
dir.create(td7)
pp <- file.path(td7, "outbreaks_survival_data_base-1990_2019_ep-gs_nz60.rds")
st1 <- list(
  pop_source = "crosswalk",
  episode_key = "gid_season",
  n_rows_out = 10L,
  stamped_at = "2026-10-01T00:00:00"
)
expect(
  identical(errmsg(assert_prep_stamp_overwrite(pp, st1)), ""),
  "no existing panel: nothing to guard"
)
saveRDS(1, pp)
expect(
  grepl("no stamp", errmsg(assert_prep_stamp_overwrite(pp, st1))),
  "an existing panel without a stamp is refused"
)
persist_prep_policy_stamp(pp, st1)
expect(
  identical(
    errmsg(assert_prep_stamp_overwrite(
      pp,
      modifyList(st1, list(stamped_at = "2026-10-27T10:00:00"))
    )),
    ""
  ),
  "same stamp but stamped_at: rebuilding is allowed"
)
st2 <- c(st1, list(contagion_coverage_report = TRUE))
expect(
  grepl(
    "contagion_coverage_report",
    errmsg(assert_prep_stamp_overwrite(pp, st2))
  ),
  "a columns-only change (new stamp field) is refused, naming the field"
)
expect(
  grepl(
    "n_rows_out",
    errmsg(assert_prep_stamp_overwrite(
      pp,
      modifyList(st1, list(n_rows_out = 11L))
    ))
  ),
  "a changed field is refused"
)
expect(
  identical(
    suppressMessages(assert_prep_stamp_overwrite(pp, st2, force = TRUE)),
    "contagion_coverage_report"
  ),
  "--force overrides (and reports the fields)"
)
src03 <- joined_source("analysis/03_prepare_occ.R")
i_guard <- grep("^assert_prep_stamp_overwrite\\(SURV_DAT_PATH", src03)
i_first_write <- min(grep(
  "^saveRDS\\(surv_dat_all|write\\.csv\\(|^persist_prep_policy_stamp\\(SURV_DAT_PATH",
  src03
))
expect(
  length(i_guard) == 1L && i_guard < i_first_write,
  "03 runs the guard before writing the panel, any sidecar or the stamp"
)
unlink(td7, recursive = TRUE)

cat("\n[7] spec defaults == get_default_options()$episodes\n")
source("analysis/utils/config.R")
gdo <- get_default_options()
for (f in names(.OCC_PANEL_SPEC_DEFAULTS)) {
  expect(
    identical(gdo$episodes[[f]], .OCC_PANEL_SPEC_DEFAULTS[[f]]),
    sprintf(
      "episodes.%s default: paths.R == config.R (%s)",
      f,
      format(.OCC_PANEL_SPEC_DEFAULTS[[f]])
    )
  )
}
expect(
  identical(occ_panel_spec_of(gdo), .occ_panel_spec(gdo$episodes$key)),
  "the resolved default opt gives the default spec"
)
expect(
  identical(
    gdo$results$inputs[c(
      "episode_national_zero_rule",
      "episode_at_risk_start",
      "episode_adm1_overlap"
    )],
    list(
      episode_national_zero_rule = "",
      episode_at_risk_start = "season",
      episode_adm1_overlap = "drop"
    )
  ),
  "results.inputs.episode_{national_zero_rule, at_risk_start, adm1_overlap} defaults"
)
expect(
  identical(gdo$occurrence$ttlo_no_prior_level, FALSE),
  "occurrence.ttlo_no_prior_level default FALSE"
)

cat(sprintf("\n[test_occ_panel_spec] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
