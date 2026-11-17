# test_basename_length.R — HARD: artefact and config names fit NAME_MAX (255
#   bytes).
#
# Usage:
#   Rscript analysis/tests/test_basename_length.R
#   Rscript analysis/tests/test_basename_length.R --dump <file.tsv>
#
# Details:
#   A file name over 255 bytes cannot be created, and the failure can come late
#   (writing cb_meta after the fit is saved, or inside 05_fit_postprocess.R's
#   tryCatch). paths.R::assert_basename_length() stops in every name builder,
#   and check_run_file_names() builds every kind of a config's track up front
#   (the fit workers and generate_batch() call it), so a long name fails before
#   a fit. --dump writes every enumerated name (the golden for a naming
#   refactor). Pure path algebra: no fit, no data.

suppressPackageStartupMessages({
  library(yaml)
})
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R") # current_pop_digest_token(), fit_pop_stamp_path(), prep_policy_stamp_path()
source("analysis/utils/enso.R") # ENSO_ANALOG_SETS
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
dump <- if ("--dump" %in% args) args[[which(args == "--dump") + 1L]] else NULL
LIMIT <- if (exists("BASENAME_MAX_BYTES")) BASENAME_MAX_BYTES else 255L

CONFIG_ROOT <- "analysis/configs"
TRACK_BUILDERS <- list(
  occurrence = "make_occ_run_file",
  duration = "make_duration_run_file",
  monthly = "make_monthly_run_file",
  size = "make_size_run_file"
)
DATA_BUILDERS <- c(
  "make_surv_dat_file",
  "make_duration_dat_file",
  "make_monthly_dat_file",
  "make_size_dat_file",
  "make_sim_panel_file",
  "make_sim_panel_projection_file",
  "make_footprint_members_file"
)
kinds_of <- function(fn) eval(formals(get(fn))$kind)

# --- Enumerate every name -----------------------------------------------------
rows <- list()
add <- function(src, builder, kind, path) {
  rows[[length(rows) + 1L]] <<- data.frame(
    source = src,
    builder = builder,
    kind = kind,
    path = as.character(path),
    stringsAsFactors = FALSE
  )
}
# A builder call; an error is recorded as a row so that it fails [2] below.
try_add <- function(src, builder, kind, expr) {
  p <- tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "err")
  })
  if (inherits(p, "err")) {
    add(src, builder, kind, paste0("<error> ", unclass(p)))
  } else {
    for (q in p) {
      add(src, builder, kind, q)
    }
  }
}

yml <- sort(list.files(
  CONFIG_ROOT,
  "\\.ya?ml$",
  recursive = TRUE,
  full.names = TRUE
))
yml <- yml[!basename(yml) %in% c("manifest.yml", "LAG_SELECTION.yml")]
is_results <- startsWith(yml, file.path(CONFIG_ROOT, "results/"))
is_deploy <- startsWith(yml, file.path(CONFIG_ROOT, "deploy/"))
opt_errors <- character()
batch_groups <- list()
fit_opts <- list()

for (f in yml[!is_results & !is_deploy]) {
  raw <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  if (!is.list(raw)) {
    next
  }
  o <- tryCatch(
    {
      out <- NULL
      invisible(capture.output(
        out <- suppressMessages(suppressWarnings(
          make_options_from_config(
            option_list = list(),
            cli_argv = c("--config", f)
          )
        ))
      ))
      out
    },
    error = function(e) e
  )
  # A file without a track or simulator block (the generators' *_horizons
  # parameter files) is not a pipeline config; one with either must resolve.
  if (inherits(o, "error")) {
    if (length(intersect(c(names(TRACK_BUILDERS), "simulator"), names(raw)))) {
      opt_errors[f] <- conditionMessage(o)
    }
    next
  }
  tracks <- intersect(names(TRACK_BUILDERS), names(raw))
  if (!length(tracks) && is.null(raw[["simulator"]])) {
    tracks <- "occurrence"
  }
  for (tr in tracks) {
    fn <- TRACK_BUILDERS[[tr]]
    for (k in kinds_of(fn)) {
      try_add(f, fn, k, get(fn)(opt = o, kind = k))
    }
    fit_opts[[length(fit_opts) + 1L]] <- list(file = f, track = tr, opt = o)
  }
  for (fn in DATA_BUILDERS) {
    try_add(f, fn, "", get(fn)(o))
  }
  if (!is.null(raw[["simulator"]])) {
    for (fn in c(
      "make_sim_run_file",
      "make_smooth_run_file",
      "make_spatpomp_run_file"
    )) {
      for (k in kinds_of(fn)) {
        try_add(f, fn, k, get(fn)(o, kind = k))
      }
    }
  }
  # Sidecars built outside paths.R by sub / paste on a builder's name.
  if ("occurrence" %in% tracks) {
    fit <- tryCatch(make_occ_run_file(o, "fit"), error = function(e) NULL)
    # 07b_simulator_v1.R / simulations.R / smoothing_core.R: cb_meta from the fit
    if (!is.null(fit)) {
      add(
        f,
        "derived",
        "fit_to_cb_meta",
        sub("_fit_base-", "_cb_meta_base-", fit, fixed = TRUE)
      )
    }
    if (identical(o$occurrence$backend, "brms_kl")) {
      # sidecars of the kl fits only
      pin <- tryCatch(make_occ_run_file(o, "kl_pin"), error = function(e) NULL)
      part <- tryCatch(
        make_occ_run_file(o, "kl_partition"),
        error = function(e) NULL
      )
      ch <- paste(
        seq_len(max(4L, as.integer(o$chains %||% 4L))),
        collapse = "-"
      )
      # make_kl_grouping.R --write_fit: <fit>_ch<chains>.rds
      if (!is.null(fit)) {
        add(
          f,
          "derived",
          "fit_ch",
          sub("\\.rds$", sprintf("_ch%s.rds", ch), fit)
        )
      }
      if (!is.null(pin)) {
        add(f, "derived", "kl_pin_arr", sub("\\.rds$", "_arr99.rds", pin))
      } # dlnm_kernel.R
      if (!is.null(part)) {
        add(f, "derived", "kl_partition_csv", sub("\\.rds$", ".csv", part))
      } # 04_fit_occ_brms_kl.R
    }
    surv <- tryCatch(as.character(make_surv_dat_file(o)), error = function(e) {
      NULL
    })
    if (!is.null(surv)) {
      add(
        f,
        "derived",
        "surv_sidecars",
        c(
          unlist(occ_panel_sidecar_paths(surv)),
          occ_panel_lp_members_path(surv),
          contagion_coverage_csv_path(surv),
          prep_policy_stamp_path(surv),
          sub("\\.rds$", ".min_seas_source.csv", surv)
        )
      )
    }
  }
  for (tr in intersect(c("duration", "monthly"), tracks)) {
    fit <- tryCatch(
      make_run_file(o, if (tr == "duration") "dur_fit" else "mon_fit"),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      add(
        f,
        "derived",
        "fit_to_cb_meta",
        sub("_fit_base-", "_cb_meta_base-", fit, fixed = TRUE)
      )
    }
  }
  if ("size" %in% tracks) {
    fit <- tryCatch(make_size_run_file(o, "size_fit"), error = function(e) NULL)
    mrt <- tryCatch(
      make_size_run_file(o, "size_mean_rate"),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      add(
        f,
        "derived",
        "size_fit_sidecars",
        c(
          fit_pop_stamp_path(fit),
          paste0(tools::file_path_sans_ext(fit), ".stan")
        )
      )
    }
    if (!is.null(mrt)) {
      add(
        f,
        "derived",
        "mean_rate_csv", # 05_fit_postprocess.R
        c(sub("\\.rds$", ".csv", mrt), sub("\\.rds$", "_tier_summary.csv", mrt))
      )
    }
  }
  bd <- dirname(f)
  if (!identical(bd, CONFIG_ROOT)) {
    batch_groups[[bd]] <- unique(c(
      batch_groups[[bd]],
      o$grouping$group %||% ""
    ))
  }
}

# Batch-keyed builders: one per batch dir, every group seen in it.
for (bd in names(batch_groups)) {
  for (k in kinds_of("make_loo_compare_run_file")) {
    for (g in batch_groups[[bd]]) {
      try_add(
        bd,
        "make_loo_compare_run_file",
        k,
        make_loo_compare_run_file(
          list(),
          k,
          batch = bd,
          group = g,
          date = as.Date("2026-11-08")
        )
      )
    }
  }
  for (k in kinds_of("make_kl_sweep_run_file")) {
    try_add(
      bd,
      "make_kl_sweep_run_file",
      k,
      make_kl_sweep_run_file(list(), k, batch = bd)
    )
  }
}

# Results configs: every kind, the longest schema token and every analog set.
schema_max <- "v99"
for (f in yml[is_results]) {
  for (k in kinds_of("make_results_run_file")) {
    try_add(
      f,
      "make_results_run_file",
      k,
      make_results_run_file(
        list(),
        k,
        track = "occurrence",
        sets = ENSO_ANALOG_SETS,
        schema = schema_max,
        config_path = f
      )
    )
  }
  try_add(
    f,
    "make_figure_dir",
    "",
    make_figure_dir(list(), track = "cb", config_path = f, subdir = "adm1")
  )
}

# Deploy YAMLs: the run root, and every run tag / batch name they state.
for (f in yml[is_deploy]) {
  for (k in kinds_of("make_deploy_run_file")) {
    try_add(
      f,
      "make_deploy_run_file",
      k,
      make_deploy_run_file(list(), k, config_path = f)
    )
  }
  d <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  dep <- d$deploy %||% d
  tags <- unlist(c(dep$run_tags, dep$fits), use.names = FALSE)
  tags <- tags[is.character(tags)]
  if (length(tags)) {
    add(f, "deploy_tags", "", file.path("generated_data", basename(tags)))
  }
}

names_df <- do.call(rbind, rows)
names_df$bytes <- nchar(basename(names_df$path), type = "bytes")
names_df$max_component <- vapply(
  strsplit(names_df$path, "/", fixed = TRUE),
  function(p) max(nchar(p, type = "bytes")),
  0L
)

if (!is.null(dump)) {
  write.table(
    names_df[, c("source", "builder", "kind", "path")],
    dump,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  cat(sprintf(
    "wrote %d names (%d configs did not resolve) to %s\n",
    nrow(names_df),
    length(opt_errors),
    dump
  ))
  quit(save = "no", status = 0L)
}

# --- [1] the guard ---------------------------------------------------------------
# assert_basename_length() passes at the limit, stops one byte over naming the
# builder, the length and the longest tokens, and checks every path component.
cat("\n[1] assert_basename_length()\n")
at <- paste0(strrep("a", LIMIT - 4L), ".rds")
over <- paste0("outbreaks_x_", strrep("b", LIMIT), "_ep-gs.rds")
expect(
  identical(
    assert_basename_length(file.path("generated_data", at), "makeTest"),
    file.path("generated_data", at)
  ),
  sprintf("a %d-byte name passes and is returned unchanged", LIMIT)
)
e <- tryCatch(
  {
    assert_basename_length(file.path("generated_data", over), "makeTest")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e) &&
    grepl("makeTest", e) &&
    grepl(sprintf("%d bytes", nchar(over, type = "bytes")), e) &&
    grepl(sprintf("%d-byte limit", LIMIT), e) &&
    grepl(strrep("b", 20L), e),
  "one over the limit stops, naming the builder, the length, the limit and the longest token"
)
e2 <- tryCatch(
  {
    assert_basename_length(
      file.path("generated_data", strrep("d", LIMIT + 1L), "x.rds"),
      "makeDir"
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e2) && grepl("makeDir", e2),
  "a directory component over the limit stops too"
)
expect(
  identical(
    assert_basename_length(c("a/b.rds", "c.qs"), "makeTest"),
    c("a/b.rds", "c.qs")
  ),
  "vectorised: every element checked, input returned"
)
multi <- paste0(strrep("é", 128L), ".rds") # 260 bytes, 132 characters
expect(
  inherits(
    tryCatch(assert_basename_length(multi, "makeTest"), error = identity),
    "error"
  ),
  "the limit is in bytes, not characters"
)
expect(
  identical(
    .cap_file_name(paste0(strrep("a", 300L), ".rds")),
    assert_basename_length(
      .cap_file_name(paste0(strrep("a", 300L), ".rds")),
      "makeTest"
    )
  ),
  ".cap_file_name() output (occurrence) always passes the guard"
)

# --- [2] every enumerated name ------------------------------------------------
# Every name every builder emits for every config in analysis/configs/ (each
# track builder x kind, the data-file builders, the sim / smooth / spatpomp,
# results, deploy, loo_compare and kl_sweep builders, and the sidecars derived
# from them by sub / paste) fits the limit.
cat("\n[2] every name the builders emit for analysis/configs/\n")
errs <- names_df[startsWith(names_df$path, "<error>"), ]
expect(
  nrow(errs) == 0L,
  sprintf(
    "no builder errors (%d)%s",
    nrow(errs),
    if (nrow(errs)) {
      paste0(
        ": ",
        paste(
          head(unique(paste(errs$builder, errs$kind, errs$path)), 5L),
          collapse = " | "
        )
      )
    } else {
      ""
    }
  )
)
ok_df <- names_df[!startsWith(names_df$path, "<error>"), ]
bad <- ok_df[ok_df$max_component > LIMIT, ]
expect(
  nrow(bad) == 0L,
  sprintf(
    "%d names, every path component <= %d bytes%s",
    nrow(ok_df),
    LIMIT,
    if (nrow(bad)) {
      paste0(": ", paste(head(basename(bad$path), 3L), collapse = " | "))
    } else {
      ""
    }
  )
)
expect(
  length(opt_errors) == 0L,
  sprintf(
    "every fit / simulator config resolves%s",
    if (length(opt_errors)) {
      paste0(": ", names(opt_errors)[[1L]], " -- ", opt_errors[[1L]])
    } else {
      ""
    }
  )
)
summ <- aggregate(bytes ~ builder + kind, ok_df, max)
summ <- summ[order(-summ$bytes), ]
cat("  longest per builder x kind (bytes, margin):\n")
for (i in seq_len(min(12L, nrow(summ)))) {
  cat(sprintf(
    "    %-24s %-26s %3d  %3d\n",
    summ$builder[i],
    summ$kind[i],
    summ$bytes[i],
    LIMIT - summ$bytes[i]
  ))
}

# --- [3] check_run_file_names() -----------------------------------------------------
# check_run_file_names() passes on every fit config and stops on a config
# whose longest kind is over the limit.
cat("\n[3] check_run_file_names()\n")
chk <- vapply(
  fit_opts,
  function(x) {
    tryCatch(
      {
        check_run_file_names(x$opt, x$track)
        ""
      },
      error = function(e) conditionMessage(e)
    )
  },
  ""
)
expect(
  all(!nzchar(chk)),
  sprintf(
    "passes on all %d fit configs%s",
    length(chk),
    if (any(nzchar(chk))) paste0(": ", chk[nzchar(chk)][[1L]]) else ""
  )
)
mon <- fit_opts[vapply(
  fit_opts,
  function(x) identical(x$track, "monthly"),
  TRUE
)][[1L]]$opt
mon$grouping$enable <- TRUE
mon$grouping$source <- strrep("g", 240L)
mon$grouping$group <- "east"
e3 <- tryCatch(
  {
    check_run_file_names(mon, "monthly")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e3) && grepl("make_monthly_run_file", e3),
  "stops on a monthly config whose names are over the limit, naming the builder"
)
expect(
  inherits(
    tryCatch(check_run_file_names(mon, "nope"), error = identity),
    "error"
  ),
  "an unknown track stops"
)

# --- [4] the config tree itself -------------------------------------------------
# Every file and directory name under analysis/configs/ fits the limit.
cat("\n[4] file and directory names under analysis/configs/\n")
tree <- list.files(CONFIG_ROOT, recursive = TRUE, include.dirs = TRUE)
comp <- unique(unlist(strsplit(tree, "/", fixed = TRUE)))
lc <- nchar(comp, type = "bytes")
expect(
  all(lc <= LIMIT),
  sprintf(
    "%d names, longest %d bytes (%s)",
    length(comp),
    max(lc),
    comp[which.max(lc)]
  )
)

cat(sprintf("\nPASS: %d   FAIL: %d\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(save = "no", status = 1L)
}
