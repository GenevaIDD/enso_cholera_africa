# test_deployed_model_terms.R — HARD: the deployed model's term set is the one
#   written down in analysis/model_terms_deployed.yml.
#
# Usage:
#   Rscript analysis/tests/test_deployed_model_terms.R
#
# Details:
#   Follows analysis/configs/results/manuscript.yml -> results.deploy_config ->
#   deploy.fits.{occ,dur,mon}_dir to the three fit batches and resolves every
#   config in them; rebuilds each formula with the track's own builder on the
#   synthetic frame of test_formula_builder_golden.R; asserts one term set per
#   batch (the configs differ only in the ENSO lag and the group) and equality
#   with the yml: fixed effects, baseline, RE terms, RE groups, offsets, the
#   monthly shape terms, backend, family. Quick, no fit, no data/ access.

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/formula.R")
source("analysis/utils/paths.R") # nbr_covariate_name(): call-time dependency of build_formula_terms()

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

SPEC_PATH <- "analysis/model_terms_deployed.yml"
RESULTS_PATH <- "analysis/configs/results/manuscript.yml"
TRACK_OF <- c(occ_dir = "occurrence", dur_dir = "duration", mon_dir = "monthly")

cat("# deployed model terms vs analysis/model_terms_deployed.yml\n")

# --- Synthetic frame (test_formula_builder_golden.R's AMPLE_DAT) -------------
# A frame is needed: with dat = NULL the duration builder drops mso_band.
# Builders: utils/formula.R build_formula_terms() + integrate_cb_terms() for
# occurrence; build_duration_formula_parts() / build_monthly_formula_parts()
# for duration / monthly.
.make_dat <- function(n, region, country, adm1, adm2_gid, ttlo, mso_band) {
  d <- data.frame(
    AFRO_region = region,
    country = country,
    ADM1 = adm1,
    ADM2 = adm2_gid,
    gid = adm2_gid,
    time_to_last_outbreak = factor(ttlo),
    mso_band = mso_band,
    stringsAsFactors = FALSE
  )
  ms_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(ms_cols) <- paste0("msoms", 1:5)
  seas_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(seas_cols) <- paste0("seasms", 1:5)
  dplyr::bind_cols(d, ms_cols, seas_cols)
}
AMPLE_DAT <- .make_dat(
  n = 16,
  region = rep(c("East", "West"), each = 8),
  country = rep(c("KEN", "TZA", "NGA", "GHA"), each = 4),
  adm1 = rep(paste0("A", 1:8), each = 2),
  adm2_gid = paste0("G", 1:16),
  ttlo = rep(c("[0,3)", "[3,6)"), length.out = 16),
  mso_band = rep(c(1, 2, 3), length.out = 16)
)

# --- Crossbasis column names from the dlnm spec -------------------------------
# The crossbasis columns are named from the config's dlnm spec
# (dlnm::onebasis on a synthetic exposure, so the column count is dlnm's own),
# as build_cb_block() names them.
# build_cb_block()'s single-lag mode: one dlnm::onebasis(x, argvar) per
# variable, columns <var>_cb1..J. Only the single-lag mode is deployed; any
# other mode stops so the gate cannot pass on a spec it does not understand.
.spec_cb_meta <- function(dlnm_cfg) {
  if (!isTRUE(dlnm_cfg$enable) || length(dlnm_cfg$variables) == 0L) {
    return(list())
  }
  dflt <- dlnm_cfg$defaults %||% list()
  out <- list()
  x <- seq(-2.5, 2.5, length.out = 101)
  for (v in names(dlnm_cfg$variables)) {
    spec <- dflt
    vs <- dlnm_cfg$variables[[v]] %||% list()
    for (k in names(vs)) {
      spec[[k]] <- vs[[k]]
    }
    if (is.null(spec$single_lag) || !is.null(spec$lags)) {
      stop(
        "dlnm variable '",
        v,
        "' is not in single-lag mode; extend ",
        "test_deployed_model_terms.R before deploying it."
      )
    }
    if (!is.null(spec[["by"]])) {
      stop(
        "dlnm variable '",
        v,
        "' has a phase-split basis (by: ",
        spec[["by"]],
        "); extend ",
        "test_deployed_model_terms.R before deploying it."
      )
    }
    argvar <- spec$argvar %||% list(fun = "bs", degree = 3L, df = 4L)
    J <- ncol(do.call(dlnm::onebasis, c(list(x), argvar)))
    out[[v]] <- list(cb_col_names = sprintf("%s_cb%d", v, seq_len(J)))
  }
  out
}

# The duration / monthly crossbasis integration is inline in
# 04_fit_dur_brms_cb.R and 04_fit_mon_brms_cb.R, not a helper;
# .integrate_cb_count_track() below mirrors those lines. A change there must be mirrored here (occurrence uses
# the shared integrate_cb_terms()).
# Mirror of the inline crossbasis integration of 04_fit_dur_brms_cb.R section 6
# (honour_levels = FALSE) and 04_fit_mon_brms_cb.R section 5 (TRUE: the monthly
# re_slope_levels knob).
.integrate_cb_count_track <- function(
  parts,
  dlnm_cfg,
  cb_meta,
  re_sep,
  drop_region,
  dat,
  opt,
  honour_levels
) {
  cb_cols <- unlist(lapply(cb_meta, `[[`, "cb_col_names"), use.names = FALSE)
  if (length(cb_cols) == 0L) {
    return(parts)
  }
  dflt <- dlnm_cfg$defaults %||% list()
  for (v in names(dlnm_cfg$variables)) {
    vs <- dlnm_cfg$variables[[v]] %||% list()
    klt <- vs$keep_linear_term %||% dflt$keep_linear_term %||% FALSE
    if (!isTRUE(klt)) {
      parts$fixed_terms <- parts$fixed_terms[parts$fixed_terms != v]
      parts$re_terms <- parts$re_terms[
        !grepl(sprintf("^\\(%s[[:space:]]", v), parts$re_terms)
      ]
    }
  }
  parts$fixed_terms <- c(parts$fixed_terms, cb_cols)
  cb_groups <- re_group_by_level("gid", drop_region = drop_region)
  for (v in names(dlnm_cfg$variables)) {
    vs <- dlnm_cfg$variables[[v]] %||% list()
    if (!isTRUE(vs$re_adm1_slope %||% (dflt$re_adm1_slope %||% FALSE))) {
      next
    }
    lv <- if (honour_levels) {
      as.character(unlist(
        vs$re_slope_levels %||%
          dflt$re_slope_levels %||%
          c("region", "country", "adm1")
      ))
    } else {
      c("region", "country", "adm1")
    }
    parts$re_terms <- c(
      parts$re_terms,
      filter_re_terms_for_data(
        nested_slope_re_terms(
          cb_meta[[v]]$cb_col_names,
          cb_groups,
          re_sep,
          levels = lv
        ),
        dat,
        keep = .single_country_keep(opt, dat)
      )
    )
  }
  parts
}

# One config -> its normalised term record.
.term_record <- function(opt, track) {
  blk <- opt[[track]]
  brms_cfg <- blk$brms_cb %||% blk$brms %||% list()
  re_sep <- if (isTRUE(brms_cfg$re_uncorrelated)) "||" else "|"
  drop_reg <- isTRUE(opt$grouping$enable)
  dlnm_cfg <- blk$dlnm %||% list()
  cb_meta <- .spec_cb_meta(dlnm_cfg)
  shape <- character(0)

  if (track == "occurrence") {
    terms <- build_formula_terms(
      opt,
      n_country = 4L,
      n_region = 2L,
      dat = AMPLE_DAT
    )
    cb_cols <- unlist(lapply(cb_meta, `[[`, "cb_col_names"), use.names = FALSE)
    terms <- integrate_cb_terms(
      terms,
      dlnm_cfg,
      cb_cols,
      cb_meta,
      re_sep,
      drop_region = drop_reg
    )
    basehaz <- .occ_basehaz(brms_cfg)
    if (identical(basehaz, "ms")) {
      stop(
        "occurrence basehaz 'ms' needs the M-spline basis; extend this test."
      )
    }
    fixed <- terms$fixed_terms
    re <- terms$re_terms
    off <- terms$off_terms
    base <- .occ_baseline_block(basehaz)
    family <- "poisson"
  } else if (track == "duration") {
    parts <- build_duration_formula_parts(opt, AMPLE_DAT)
    parts <- .integrate_cb_count_track(
      parts,
      dlnm_cfg,
      cb_meta,
      re_sep,
      drop_reg,
      AMPLE_DAT,
      opt,
      honour_levels = FALSE
    )
    fixed <- parts$fixed_terms
    re <- parts$re_terms
    off <- parts$off_terms
    base <- c(parts$bh_terms, parts$seas_terms)
    family <- "poisson"
  } else {
    parts <- build_monthly_formula_parts(opt, AMPLE_DAT)
    parts <- .integrate_cb_count_track(
      parts,
      dlnm_cfg,
      cb_meta,
      re_sep,
      drop_reg,
      AMPLE_DAT,
      opt,
      honour_levels = TRUE
    )
    fixed <- parts$fixed_terms
    re <- parts$re_terms
    off <- parts$off_terms
    base <- parts$seas_terms
    shape <- .phi_re_term(
      parts$phi_re,
      drop_reg,
      AMPLE_DAT,
      keep = .single_country_keep(opt, AMPLE_DAT)
    )
    family <- blk$family %||% "negbinomial"
  }
  ctg <- opt$contagion %||% list()
  list(
    backend = blk$backend,
    nbr = nbr_covariate_name(
      ctg$rho_km %||% 500L,
      ctg$scale,
      ctg$partition,
      ctg$partition_source
    ),
    family = family,
    link = blk$link %||% "log",
    fixed = sort(setdiff(fixed, "1")),
    baseline = sort(base),
    re = sort(re),
    offsets = sort(off),
    shape = sort(shape),
    dup = anyDuplicated(c(fixed, re)) > 0L
  )
}

# yml entry list -> character vector ({k} expanded over 1..count; {nbr} = the
# config's nbr_covariate_name(), so a renamed kernel column needs no yml edit).
.expand <- function(entries, nbr = NULL) {
  out <- character(0)
  for (e in entries %||% list()) {
    if (is.list(e)) {
      if (is.null(e[["term"]]) || !all(names(e) %in% c("term", "count"))) {
        stop(
          "yml entry must be a string or {term: ..., count: ...}: ",
          paste(names(e), collapse = ", ")
        )
      }
      n <- as.integer(e[["count"]] %||% 1L)
      out <- c(
        out,
        vapply(
          seq_len(n),
          function(k) gsub("{k}", k, e$term, fixed = TRUE),
          character(1)
        )
      )
    } else {
      out <- c(out, as.character(e))
    }
  }
  if (!is.null(nbr)) {
    out <- gsub("{nbr}", nbr, out, fixed = TRUE)
  }
  if (any(grepl("{nbr}", out, fixed = TRUE))) {
    stop("unresolved {nbr} placeholder")
  }
  sort(out)
}

# RE groups brms creates from the RE terms: the part after `|` / `||`, with
# `a/b/c` nesting expanded to a, a:b, a:b:c.
.re_groups <- function(re) {
  g <- trimws(sub("^.*\\|\\s*", "", sub("\\)\\s*$", "", re)))
  sort(unique(unlist(lapply(strsplit(g, "/", fixed = TRUE), function(p) {
    vapply(
      seq_along(p),
      function(i) paste(p[seq_len(i)], collapse = ":"),
      character(1)
    )
  }))))
}

.same_set <- function(got, want, what) {
  extra <- setdiff(got, want)
  miss <- setdiff(want, got)
  expect(
    length(extra) == 0L && length(miss) == 0L,
    sprintf(
      "%s match the yml%s%s",
      what,
      if (length(extra)) {
        paste0(" | NOT in yml: ", paste(extra, collapse = ", "))
      } else {
        ""
      },
      if (length(miss)) {
        paste0(" | missing: ", paste(miss, collapse = ", "))
      } else {
        ""
      }
    )
  )
}

# --- [1] follow the results config to the fit batches -------------------------
cat("\n[1] results config -> deploy config -> fit batches\n")
spec <- yaml::read_yaml(SPEC_PATH)
res_cfg <- yaml::read_yaml(RESULTS_PATH)
deploy_path <- res_cfg$results$deploy_config %||% ""
expect(
  nzchar(deploy_path) && file.exists(deploy_path),
  sprintf("results.deploy_config resolves (%s)", deploy_path)
)
if (!nzchar(deploy_path) || !file.exists(deploy_path)) {
  quit(status = 1)
}
cat(sprintf(
  "  deploy: %s (yml written against %s)\n",
  deploy_path,
  spec$describes$deploy_config %||% "?"
))
fits <- yaml::read_yaml(deploy_path)$deploy$fits

# The deployed formulas also depend on code the builders do not see.
occ_src <- readLines("analysis/04_fit_occ_brms_cb.R", warn = FALSE)
dur_src <- readLines("analysis/04_fit_dur_brms_cb.R", warn = FALSE)
expect(
  any(grepl("family\\s*=\\s*poisson\\(\\)", occ_src)) &&
    any(grepl("family\\s*=\\s*poisson\\(\\)", dur_src)),
  "occurrence and duration brms_cb scripts fit family = poisson()"
)

# --- [2] every config of every batch -----------------------------------------
for (key in names(TRACK_OF)) {
  track <- TRACK_OF[[key]]
  want <- spec$tracks[[track]]
  cat(sprintf("\n[2] %s (%s)\n", track, key))
  expect(
    identical(want$batch_key, key),
    sprintf("yml %s.batch_key is %s", track, key)
  )

  batch <- fits[[key]] %||% ""
  bdir <- file.path("analysis/configs", batch)
  cfgs <- sort(setdiff(
    list.files(bdir, pattern = "\\.yml$", full.names = TRUE),
    file.path(bdir, "manifest.yml")
  ))
  expect(
    nzchar(batch) && length(cfgs) > 0L,
    sprintf("batch %s has %d config(s)", batch, length(cfgs))
  )
  if (length(cfgs) == 0L) {
    next
  }

  recs <- lapply(cfgs, function(f) {
    tryCatch(
      {
        invisible(capture.output(
          opt <- suppressWarnings(suppressMessages(
            make_options_from_config(
              option_list = list(),
              cli_argv = c("--config", f)
            )
          ))
        ))
        if (is.null(opt[[track]])) {
          stop("no ", track, " block")
        }
        suppressWarnings(suppressMessages(.term_record(opt, track)))
      },
      error = function(e) list(.error = conditionMessage(e))
    )
  })
  errs <- vapply(recs, function(r) !is.null(r$.error), logical(1))
  expect(
    !any(errs),
    sprintf(
      "all %d configs resolve and build (errors: %s)",
      length(cfgs),
      paste(
        head(unique(unlist(lapply(recs[errs], `[[`, ".error"))), 3),
        collapse = " / "
      )
    )
  )
  if (any(errs)) {
    next
  }

  u <- unique(recs)
  expect(
    length(u) == 1L,
    sprintf("one term set across the batch (%d distinct)", length(u))
  )
  r <- u[[1L]]

  expect(
    identical(r$backend, want$backend),
    sprintf("backend = %s", want$backend)
  )
  expect(
    identical(r$family, want$family) && identical(r$link, want$link %||% "log"),
    sprintf("family = %s(%s)", want$family, want$link %||% "log")
  )
  expect(!isTRUE(r$dup), "no duplicated term")
  .same_set(r$fixed, .expand(want$fixed_effects, r$nbr), "fixed effects")
  .same_set(r$baseline, .expand(want$baseline, r$nbr), "baseline terms")
  .same_set(r$re, .expand(want$re_terms, r$nbr), "RE terms")
  .same_set(.re_groups(r$re), sort(unlist(want$re_groups)), "RE groups")
  .same_set(r$offsets, .expand(want$offsets, r$nbr), "offsets")
  .same_set(
    r$shape,
    .expand(want$shape_terms, r$nbr),
    "shape (dispersion) terms"
  )

  absent <- as.character(unlist(want$not_in_model))
  hit <- absent[vapply(
    absent,
    function(a) {
      any(startsWith(r$fixed, a)) || any(startsWith(r$re, paste0("(", a)))
    },
    logical(1)
  )]
  expect(
    length(hit) == 0L,
    sprintf(
      "none of not_in_model in any term%s",
      if (length(hit)) {
        paste0(" (found: ", paste(hit, collapse = ", "), ")")
      } else {
        ""
      }
    )
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
