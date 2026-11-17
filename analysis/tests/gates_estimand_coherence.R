# gates_estimand_coherence.R — checks on the built artefacts that the size and
#   monthly estimands mean the same thing at every seam they cross.
#
# Usage:
#   Rscript analysis/tests/gates_estimand_coherence.R
#   Rscript analysis/tests/gates_estimand_coherence.R --config <size cfg.yml>
#   Rscript analysis/tests/gates_estimand_coherence.R --frame_path <frame.rds>
#   Rscript analysis/tests/gates_estimand_coherence.R --eta_target 1.0
#
# Details:
#   Not a test_*.R file and not in run_all.R: it checks the artefacts on disk,
#   not the code, so it runs only where they exist. The size mean
#   mu_i = pop_i * mean_rate_{c(i)} * exp(eta_i) has BOTH offsets fixed, so the
#   absolute scale lives in the two denominators: a whole-gid population (Q)
#   offset paired with a reporting-footprint (H) mean_rate multiplies every
#   prediction by Q/H, and nothing in the model can absorb it. The gates are
#   listed above "Resolve the declared estimand"; (7)/(8) SKIP until
#   --eta_target / --frame_target are passed.

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

source("analysis/utils/base.R") # transitive dependency of provenance.R, paths.R, siz_data.R (utils/README.md)
source("analysis/utils/gid.R") # transitive dependency of siz_data.R (utils/README.md)
source("analysis/utils/paths.R") # make_pop_run_file() -- used directly
source("analysis/utils/clim_compile.R") # transitive dependency of siz_data.R (utils/README.md)
source("analysis/utils/clim_db.R") # transitive dependency of siz_data.R (utils/README.md)
source("analysis/utils/siz_data.R") # read_size_data() -- used directly
source("analysis/utils/population.R") # apply_canonical_admin_level()/canonical_admin_level_table() -- used directly
source("analysis/utils/composite_lp.R") # transitive dependency of population.R (utils/README.md)
source("analysis/utils/provenance.R") # compare_size_estimand()/resolve_size_estimand()/size_estimand_of() -- used directly
source("analysis/utils/config.R") # NOT in the facade's module list
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/anchor.R") # transitive dependency of siz_data.R (utils/README.md); anchor_enabled() -- used directly
source("analysis/utils/dur_data.R") # transitive dependency of mon_data.R (the monthly panel is the duration panel + monthly cases)
source("analysis/utils/mon_data.R") # load_monthly_data() -- used directly by the MONTHLY TRACK section below

# NOTE two prefix hazards, both hit while writing this file:
#   * `--frame_path` is not called `--frame` because optparse rejects a long flag
#     that is a strict prefix of another (`--frame_target`) as ambiguous.
#   * every read below uses `opt_cli[["k"]]`, never `opt_cli$k`: `$` partial-
#     matches on lists, so `opt_cli$frame` silently returned `frame_target`.
opt_cli <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("-C", "--config"),
      type = "character",
      default = NULL,
      help = "size config YAML; declares the estimand to audit against"
    ),
    make_option(
      "--size_data",
      type = "character",
      default = "generated_data/outbreak_size_analysis_data.rds",
      help = "prepared size dataset [default: %default]"
    ),
    make_option(
      "--mon_data",
      type = "character",
      default = NULL,
      help = "prepared monthly panel [default: make_monthly_dat_file(opt)]"
    ),
    make_option(
      "--gid_pop_path",
      type = "character",
      default = NULL,
      help = "canonical gid population .qs [default: from config/newest]"
    ),
    make_option(
      "--frame_path",
      type = "character",
      default = NULL,
      help = "simulator location frame .rds, for gates (1)/(3)"
    ),
    make_option(
      "--eta_target",
      type = "double",
      default = NA_real_,
      help = "gate (7) centre; UNSET => skip"
    ),
    make_option(
      "--frame_target",
      type = "double",
      default = NA_real_,
      help = "gate (8) centre; UNSET => skip"
    ),
    make_option(
      "--tol",
      type = "double",
      default = 0.25,
      help = "relative band for gates (7)/(8) [default: %default]"
    )
  )
))

PASS <- 0L
FAIL <- 0L
SKIP <- 0L
gate <- function(cond, label, detail = NULL) {
  if (isTRUE(cond)) {
    PASS <<- PASS + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  FAIL  %s\n", label))
  }
  if (!is.null(detail)) cat(sprintf("        %s\n", detail))
}
skip <- function(label, why) {
  SKIP <<- SKIP + 1L
  cat(sprintf("  SKIP  %s\n        (%s)\n", label, why))
}

# Why gates and not provenance: assert_pop_provenance() compares population
# DIGESTS, so it fires when the two sides used different population *numbers*;
# `Q x rbar^H` uses the right numbers on both sides and is still wrong, because
# it uses them as different *units*. Digests see bytes, not semantics.
#
# GATES
#   (1) estimand identity   fit stamp == mean_rate stamp == frame stamp
#   (2) dimensional         offset denominator units == rate denominator units
#   (3) level keying        zero admin_level disagreements vs canonical gadm_level
#   (4) no silent defaults  zero pop / mean_rate NA->constant substitutions
#   (5) phantom levels      no mean_rate cell at a level with no canonical units
#   (6) elasticity          log_pop_c present and centred iff free
#   (7) training calibration    eta_fixed in [lo, hi]        -- NEEDS A REFIT
#   (8) frame extrapolation     frame ratio in [lo, hi]      -- NEEDS A REFIT
# Gates (1)/(2)/(3)/(4)/(6) are repeated on the MONTHLY panel (the
# `=== MONTHLY TRACK ===` section at the bottom): mon has the same H-vs-Q
# problem in the same shape. Gate (5) does not transfer: it is about the
# mean_rate TABLE, and mon has none.
# ---------------------------------------------------------------------------
# Resolve the declared estimand
# ---------------------------------------------------------------------------
# `cli_argv` is passed EXPLICITLY, never left NULL: make_options_from_config()
# runs its own parse_args() over commandArgs() when cli_argv is NULL, and that
# parser knows nothing of this file's flags -- so `--frame_path` died there with
# "long flag is invalid" long after our own parse_args() had accepted it.
opt <- if (!is.null(opt_cli[["config"]])) {
  make_options_from_config(
    list(),
    cli_argv = c("--config", opt_cli[["config"]])
  )
} else {
  # No config: audit the DEFAULTS, which is still meaningful -- it proves the
  # default declaration is self-consistent and that the artefacts on disk match
  # what an unconfigured run would produce.
  cat("[note] no --config; auditing analysis/config_defaults.yml as declared\n")
  make_options_from_config(list(), cli_argv = character(0))
}

es <- resolve_size_estimand(opt) # itself gate (2) for the DECLARATION
cat(sprintf(
  "\nDeclared estimand: denominator=%s  rate_denominator=%s  urm=%s  elasticity=%s\n",
  es$denominator,
  es$rate_denominator,
  format(es$underreporting_multiplier),
  es$pop_elasticity
))

gid_pop_path <- opt_cli[["gid_pop_path"]] %||%
  opt$population$gid_pop_path %||%
  tryCatch(make_pop_run_file(list(), "pop_gid"), error = function(e) NULL)

# ---------------------------------------------------------------------------
# Gate (2) — dimensional. Declaration first, then the artefacts.
# ---------------------------------------------------------------------------
cat("\n== (2) dimensional: offset denominator == rate denominator ==\n")
# resolve_size_estimand() already hard-errors on a mismatch, so reaching here IS
# the pass. Asserted again rather than assumed: this file must fail loudly if
# that validation is ever relaxed, not silently stop covering it.
gate(
  identical(es$denominator, es$rate_denominator),
  sprintf(
    "declaration is dimensionally consistent (%s x %s)",
    es$denominator,
    es$rate_denominator
  )
)

# ---------------------------------------------------------------------------
# Gates (1), (4), (5), (6) — against the built artefacts
# ---------------------------------------------------------------------------
if (!file.exists(opt_cli[["size_data"]])) {
  skip(
    "(1) estimand identity",
    sprintf("size data absent: %s", opt_cli[["size_data"]])
  )
  skip("(4) no silent defaults", "size data absent")
  skip("(5) phantom levels", "size data absent")
  skip("(6) elasticity column", "size data absent")
  sd <- NULL
  mrt <- NULL
} else {
  sd <- suppressWarnings(read_size_data(opt, path = opt_cli[["size_data"]]))
  mrt <- attr(sd, "mean_rate_table")

  cat("\n== (1) estimand identity: fit data vs mean_rate table ==\n")
  s_fit <- size_estimand_of(sd)
  s_mrt <- size_estimand_of(mrt)
  diffs <- compare_size_estimand(s_fit, s_mrt)
  gate(
    length(diffs) == 0L,
    "fit-data stamp == mean_rate-table stamp",
    if (length(diffs)) {
      paste(diffs, collapse = " ; ")
    } else {
      sprintf("both: %s x %s", s_fit$denominator, s_fit$rate_denominator)
    }
  )

  cat("\n== (4) no silent defaults ==\n")
  # A constant substituted for a missing population or rate is the failure mode
  # that makes every downstream number look plausible while being wrong. Both
  # columns feed log() offsets, so NA/0/Inf are all fatal, not merely untidy.
  bad_pop <- sum(!is.finite(sd$pop) | sd$pop <= 0)
  bad_mr <- sum(!is.finite(sd$mean_rate) | sd$mean_rate <= 0)
  gate(
    bad_pop == 0L,
    "every fitted row has a finite positive pop",
    sprintf("%d offending row(s)", bad_pop)
  )
  gate(
    bad_mr == 0L,
    "every fitted row has a finite positive mean_rate",
    sprintf("%d offending row(s)", bad_mr)
  )
  if ("mr_tier" %in% names(sd)) {
    tt <- table(sd$mr_tier)
    cat(sprintf(
      "        mr_tier: %s\n",
      paste(names(tt), tt, sep = "=", collapse = " ")
    ))
    # T5 means "nothing resolved". T4 is the global mean -- not a default, but a
    # borrow so distant it should never carry real mass silently.
    gate(!any(names(tt) == "T5"), "no row falls through to T5 (unresolved)")
  }

  cat("\n== (5) phantom levels: no mean_rate cell at an unpopulated level ==\n")
  if (is.null(gid_pop_path) || !file.exists(gid_pop_path)) {
    skip("(5) phantom levels", "canonical gid population table not found")
  } else {
    canon_lv <- sort(unique(
      canonical_admin_level_table(gid_pop_path)$admin_level_canonical
    ))
    cell_lv <- sort(unique(as.character(mrt$cell$admin_level)))
    phantom <- setdiff(cell_lv, canon_lv)
    gate(
      length(phantom) == 0L,
      "mean_rate cells exist only at levels the canonical table has",
      sprintf(
        "canonical: %s | cells: %s%s",
        paste(canon_lv, collapse = "/"),
        paste(cell_lv, collapse = "/"),
        if (length(phantom)) {
          sprintf(" | PHANTOM: %s", paste(phantom, collapse = "/"))
        } else {
          ""
        }
      )
    )
  }

  cat("\n== (6) elasticity column present and centred iff free ==\n")
  has_c <- "log_pop_c" %in% names(sd)
  if (identical(es$pop_elasticity, "fixed")) {
    gate(
      !has_c,
      "pop_elasticity=fixed => no log_pop_c column (beta == 1 imposed)"
    )
  } else {
    gate(has_c, "pop_elasticity!=fixed => log_pop_c present")
    if (has_c) {
      # Centring is load-bearing, not cosmetic: uncentred, the deviation term
      # shifts the intercept as well as the slope and `beta = 1` no longer
      # reproduces the offset-only fit, so the two paths stop being nested.
      ctr <- if ("admin_level" %in% names(sd)) {
        tapply(sd$log_pop_c, as.character(sd$admin_level), mean)
      } else {
        mean(sd$log_pop_c)
      }
      gate(
        max(abs(ctr)) < 1e-8,
        "log_pop_c is centred within admin_level",
        sprintf("max |mean| = %.2e", max(abs(ctr)))
      )
    }
  }
}

# ---------------------------------------------------------------------------
# Gate (3) — level keying, on both the fit data and the frame
# ---------------------------------------------------------------------------
cat("\n== (3) level keying: admin_level == canonical gadm_level ==\n")
if (is.null(gid_pop_path) || !file.exists(gid_pop_path)) {
  skip("(3) level keying", "canonical gid population table not found")
} else {
  chk <- function(x, what) {
    y <- apply_canonical_admin_level(
      x,
      gid_pop_path,
      overwrite = FALSE,
      tag = "gate3"
    )
    r <- attr(y, "admin_level_recode")
    gate(
      r$n_changed == 0L && r$n_unmatched == 0L,
      sprintf("%s is canonically keyed", what),
      sprintf(
        "%d/%d row(s) mis-keyed, %d unmatched%s",
        r$n_changed,
        r$n_rows,
        r$n_unmatched,
        if (!is.null(r$crosstab)) {
          sprintf(
            " [%s]",
            paste(
              sprintf(
                "%s->%s:%d",
                r$crosstab$corpus,
                r$crosstab$canonical,
                r$crosstab$Freq
              ),
              collapse = " "
            )
          )
        } else {
          ""
        }
      )
    )
  }
  if (!is.null(sd)) {
    chk(sd, "fitted size data")
  } else {
    skip("(3) level keying [size data]", "size data absent")
  }

  if (
    !is.null(opt_cli[["frame_path"]]) && file.exists(opt_cli[["frame_path"]])
  ) {
    fr <- readRDS(opt_cli[["frame_path"]])
    chk(fr, "simulator frame")
    if (!is.null(sd)) {
      # The two sides can each be canonical and still disagree with each other
      # if one was built before the re-keying landed. Compare directly.
      j <- dplyr::inner_join(
        dplyr::distinct(dplyr::transmute(
          as.data.frame(sd),
          gid = as.character(gid),
          lv_fit = as.character(admin_level)
        )),
        dplyr::distinct(dplyr::transmute(
          as.data.frame(fr),
          gid = as.character(gid),
          lv_frame = as.character(admin_level)
        )),
        by = "gid"
      )
      nd <- sum(j$lv_fit != j$lv_frame)
      gate(
        nd == 0L,
        "fit data and frame agree on every shared gid's level",
        sprintf("%d of %d shared gid(s) disagree", nd, nrow(j))
      )
    }
  } else {
    skip(
      "(3) level keying [frame]",
      if (is.null(opt_cli[["frame_path"]])) {
        "--frame_path not given"
      } else {
        "frame file absent"
      }
    )
  }
}

# ---------------------------------------------------------------------------
# Gates (7)/(8) — calibration. Unset by construction; see the note below.
# ---------------------------------------------------------------------------
# GATES (7) AND (8) ARE DELIBERATELY UNSET. Their thresholds cannot be derived
# from the current artefacts: the composed post-fix prediction depends on which
# moment a refit is assumed to re-normalise, and that assumption is not
# testable without the refit. They SKIP (never PASS) until --eta_target /
# --frame_target are passed explicitly, because a gate that cannot run must
# say so rather than report a pass.
# GATE (8) IS PER-COUNTRY, NOT GLOBAL. The global aggregate can be a
# cancellation of ~5x opposite-signed country errors, and a region can get
# WORSE under a fix while the global number improves: a global-only
# calibration gate would pass on a model that is broken everywhere.
cat("\n== (7)/(8) calibration ==\n")
if (is.na(opt_cli[["eta_target"]])) {
  skip(
    "(7) training calibration (eta_fixed)",
    paste0(
      "no --eta_target. This threshold CANNOT be derived from the ",
      "current artefacts -- it must come from an actual refit."
    )
  )
} else if (is.null(sd)) {
  skip("(7) training calibration (eta_fixed)", "size data absent")
} else {
  # eta_fixed = case-weighted geometric mean of the offset multiplier relative
  # to the config the deployed fit used. A refit drives it to 1, so its
  # distance from 1 is how far eta will move.
  eta <- exp(stats::weighted.mean(
    log(sd$pop * sd$mean_rate),
    sd$suspected_cases,
    na.rm = TRUE
  ))
  lo <- opt_cli[["eta_target"]] * (1 - opt_cli[["tol"]])
  hi <- opt_cli[["eta_target"]] * (1 + opt_cli[["tol"]])
  gate(
    eta >= lo && eta <= hi,
    sprintf("eta_fixed in [%.3f, %.3f]", lo, hi),
    sprintf("eta_fixed = %.4f", eta)
  )
}
if (is.na(opt_cli[["frame_target"]])) {
  skip(
    "(8) frame extrapolation, PER COUNTRY",
    paste0(
      "no --frame_target. Needs a refit AND simulated output; and it ",
      "must be asserted per country, not on the aggregate -- the ",
      "global ratio is a cancellation of ~5x opposite-signed errors."
    )
  )
} else {
  skip(
    "(8) frame extrapolation, PER COUNTRY",
    "requires a completed simulation run; wire in 07c --validate output"
  )
}

# ---------------------------------------------------------------------------
# Gates (2)/(3)/(4)/(6) — repeated on the MONTHLY panel
# ---------------------------------------------------------------------------
# The mon track has size's H-vs-Q problem in exactly the same shape (a ~76x
# CMR over-prediction), so it needs the same audit. Bolted on ADDITIVELY
# rather than by making the whole file `--track`-aware, for two reasons:
#
#   * the size gates above are the reference run whose pass/fail counts are
#     compared across commits; rewriting their control flow to thread a track
#     through would make that comparison meaningless exactly when it matters;
#   * gates (1) and (5) do not transfer as written. They are about the
#     mean_rate TABLE, and mon has no mean_rate table at all -- its analogue is
#     the base-rate anchor (utils/anchor.R), which is OFF by default
#     (`monthly.offset_baserate`), so a mon run usually has no second
#     denominator artefact to be incoherent with. Gate (1) below is therefore
#     the weaker, honest version: fit-data stamp vs anchor stamp, skipped when
#     the anchor is not enabled.
#
# Skips by default: pass --mon_data to audit the monthly panel.
cat("\n\n=== MONTHLY TRACK ===\n")

mon_path <- opt_cli[["mon_data"]]
if (is.null(mon_path)) {
  mon_path <- tryCatch(make_monthly_dat_file(opt), error = function(e) NULL)
}

es_mon <- resolve_monthly_estimand(opt) # itself gate (2) for the DECLARATION
cat(sprintf(
  "Declared monthly estimand: denominator=%s  rate_denominator=%s  urm=%s  elasticity=%s\n",
  es_mon$denominator,
  es_mon$rate_denominator,
  format(es_mon$underreporting_multiplier),
  es_mon$pop_elasticity
))

cat("\n== (2-mon) dimensional: offset denominator == rate denominator ==\n")
gate(
  identical(es_mon$denominator, es_mon$rate_denominator),
  sprintf(
    "monthly declaration is dimensionally consistent (%s x %s)",
    es_mon$denominator,
    es_mon$rate_denominator
  )
)

if (is.null(mon_path) || !file.exists(mon_path)) {
  for (g in c(
    "(1-mon) estimand identity",
    "(3-mon) level keying",
    "(4-mon) no silent defaults",
    "(6-mon) elasticity column"
  )) {
    skip(g, sprintf("monthly panel absent: %s", mon_path %||% "<unresolved>"))
  }
} else {
  md <- suppressWarnings(suppressMessages(load_monthly_data(
    opt,
    path = mon_path
  )))

  cat("\n== (1-mon) estimand identity: fit data vs base-rate anchor ==\n")
  s_mon <- monthly_estimand_of(md)
  anchor_tbl <- attr(md, "base_rate_table") # attach_base_rate(), gated on monthly.offset_baserate
  if (is.null(anchor_tbl)) {
    skip(
      "(1-mon) estimand identity",
      paste0(
        "monthly.offset_baserate is off, so there is no second ",
        "denominator artefact to be incoherent with"
      )
    )
  } else {
    diffs <- compare_size_estimand(
      s_mon,
      size_estimand_of(anchor_tbl, track = "monthly")
    )
    gate(
      length(diffs) == 0L,
      "monthly fit-data stamp == base-rate-anchor stamp",
      if (length(diffs)) {
        paste(diffs, collapse = " ; ")
      } else {
        sprintf("both: %s x %s", s_mon$denominator, s_mon$rate_denominator)
      }
    )
  }

  cat("\n== (4-mon) no silent defaults ==\n")
  # Only `pop` here: unlike size, mon's mean is pop x exp(eta) with no rate
  # offset in the default configuration, so `pop` is the ENTIRE absolute scale
  # and a single non-finite row is a log() of a non-number, not a rounding
  # nuisance.
  bad_pop <- sum(!is.finite(md$pop) | md$pop <= 0)
  gate(
    bad_pop == 0L,
    "every fitted monthly row has a finite positive pop",
    sprintf("%d offending row(s) of %d", bad_pop, nrow(md))
  )
  # The response is rounded for family = negbinomial(); a negative or
  # non-finite count would be a hard brms error, so catch it here instead.
  bad_y <- sum(!is.finite(md$cases_month) | md$cases_month < 0)
  gate(
    bad_y == 0L,
    "every fitted monthly row has a finite non-negative cases_month",
    sprintf("%d offending row(s)", bad_y)
  )

  cat("\n== (6-mon) elasticity column present and centred iff free ==\n")
  has_c <- "log_pop_c" %in% names(md)
  if (identical(es_mon$pop_elasticity, "fixed")) {
    gate(
      !has_c,
      "monthly pop_elasticity=fixed => no log_pop_c column (beta == 1 imposed)"
    )
  } else {
    gate(has_c, "monthly pop_elasticity!=fixed => log_pop_c present")
    if (has_c) {
      ctr <- if ("admin_level" %in% names(md)) {
        tapply(md$log_pop_c, as.character(md$admin_level), mean)
      } else {
        mean(md$log_pop_c)
      }
      gate(
        max(abs(ctr)) < 1e-8,
        "monthly log_pop_c is centred within admin_level",
        sprintf("max |mean| = %.2e", max(abs(ctr)))
      )
      # The centres must also have TRAVELLED, or the simulator re-centres frame
      # units on their own mean and applies the elasticity around the wrong
      # pivot -- silently rescaling every prediction by exp((beta-1)*shift).
      gate(
        !is.null(monthly_estimand_of(md)[["log_pop_centres"]]),
        "monthly stamp carries log_pop_centres for the simulator to re-centre on"
      )
    }
  }

  cat("\n== (3-mon) level keying: admin_level == canonical gadm_level ==\n")
  # EXPECTED TO FAIL under the deployed default, and that is the point: P3
  # (`monthly.canonical_admin_level`) is off unless the estimand opts in, so the
  # panel still carries OUTBREAK-CORPUS labels. Measured on the current panel:
  # 2626/7875 rows mis-keyed, including a phantom ADM3 the canonical table does
  # not have; with P3 on, 0/6673 and the ADM3 label is gone. (The row count
  # drops because apply_admin_level_filter() then partitions on canonical labels
  # instead of corpus ones -- which is exactly why P3 must run BEFORE it.)
  if (is.null(gid_pop_path) || !file.exists(gid_pop_path)) {
    skip("(3-mon) level keying", "canonical gid population table not found")
  } else {
    y <- apply_canonical_admin_level(
      md,
      gid_pop_path,
      overwrite = FALSE,
      tag = "gate3-mon"
    )
    r <- attr(y, "admin_level_recode")
    gate(
      r$n_changed == 0L && r$n_unmatched == 0L,
      "fitted monthly data is canonically keyed",
      sprintf(
        "%d/%d row(s) mis-keyed, %d unmatched%s%s",
        r$n_changed,
        r$n_rows,
        r$n_unmatched,
        if (!is.null(r$crosstab)) {
          sprintf(
            " [%s]",
            paste(
              sprintf(
                "%s->%s:%d",
                r$crosstab$corpus,
                r$crosstab$canonical,
                r$crosstab$Freq
              ),
              collapse = " "
            )
          )
        } else {
          ""
        },
        if (r$n_changed > 0L && !isTRUE(opt$monthly$canonical_admin_level)) {
          " -- set monthly.canonical_admin_level: true"
        } else {
          ""
        }
      )
    )
  }
}

cat(sprintf("\n%d pass, %d fail, %d skip\n", PASS, FAIL, SKIP))
if (FAIL > 0) {
  quit(status = 1)
}
