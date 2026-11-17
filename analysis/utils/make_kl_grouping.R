# make_kl_grouping.R — country grouping registry from a brms_kl latent-class
#   fit: each fitted country goes to its class of maximum mean probability
#
# Usage:
#   Rscript analysis/utils/make_kl_grouping.R -C <kl config> --out <name> \
#     [--chains 2,6] [--write_fit]
#   Rscript analysis/utils/make_kl_grouping.R --fit <fit.rds> --out <name> \
#     [--fill_from <registry>] [--force]
#
# Details:
#   Reads a saved kl_fit (04_fit_occ_brms_kl.R); --chains keeps a subset
#   (dlnm_kernel.R::kl_subset_chains(); --write_fit saves it beside the fit as
#   <fit>_ch<a>-<b>.rds). Writes analysis/groupings/<out>.csv (iso3,group; the
#   fit's class labels), validated by load_grouping(), and the per-country
#   probabilities to generated_data/kl_groupings/<out>_membership.csv. Only
#   fitted countries are written unless --fill_from carries the others over.
#
# See: utils/kl_grouping.R::kl_grouping_from_fit(), shared with gen_grouping

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})
args <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("-C", "--config"),
      default = "",
      type = "character",
      help = "brms_kl config whose fit to read (make_run_file(opt, 'fit'))"
    ),
    make_option(
      "--fit",
      default = "",
      type = "character",
      help = "path to a saved kl_fit .rds (instead of -C)"
    ),
    make_option(
      "--chains",
      default = "",
      type = "character",
      help = "comma-separated chain ids to keep, e.g. 2,6 [default: all]"
    ),
    make_option(
      "--write_fit",
      default = FALSE,
      action = "store_true",
      help = "save the chain subset as <fit>_ch<ids>.rds beside the fit"
    ),
    make_option(
      "--out",
      default = "",
      type = "character",
      help = "registry name to write under analysis/groupings/ (required)"
    ),
    make_option(
      "--fill_from",
      default = "",
      type = "character",
      help = "registry whose countries absent from the fit are carried over [default: none]"
    ),
    make_option(
      "--force",
      default = FALSE,
      action = "store_true",
      help = "overwrite an existing analysis/groupings/<out>.csv"
    )
  )
))
for (m in c("base", "config", "paths", "dlnm", "dlnm_kernel", "kl_grouping")) {
  source(sprintf("analysis/utils/%s.R", m))
}

if (!nzchar(args$out)) {
  stop("--out <registry name> is required.")
}
out_csv <- file.path("analysis/groupings", paste0(args$out, ".csv"))
if (file.exists(out_csv) && !args$force) {
  stop(out_csv, " exists (use --force to overwrite).")
}

fit_path <- if (nzchar(args$fit)) {
  args$fit
} else if (nzchar(args$config)) {
  opt <- NULL
  invisible(capture.output(suppressMessages(
    opt <- make_options_from_config(
      option_list = list(),
      cli_argv = c("--config", args$config),
      model = "occurrence"
    )
  )))
  make_run_file(opt, "fit")
} else {
  stop("give -C <config> or --fit <path>.")
}
if (!file.exists(fit_path)) {
  stop("fit not found: ", fit_path)
}
fit <- readRDS(fit_path)
if (!inherits(fit, "kl_fit")) {
  stop(fit_path, " is not a brms_kl latent-class fit (class kl_fit).")
}
cat(
  "Fit:",
  fit_path,
  "\n  chains:",
  fit$chains,
  " classes:",
  paste(fit$class_labels, collapse = ", "),
  "\n"
)

ch <- if (nzchar(args$chains)) {
  as.integer(strsplit(args$chains, ",")[[1]])
} else {
  NULL
}
if (args$write_fit && is.null(ch)) {
  stop("--write_fit needs --chains.")
}
# Class of maximum mean posterior probability: kl_responsibilities() ->
# kl_map_partition(method = "argmax").
g <- kl_grouping_from_fit(fit, ch, args$fill_from)
if (!is.null(ch)) {
  cat("  kept chains:", paste(g$fit$source_chains, collapse = ", "), "\n")
  if (args$write_fit) {
    sub_path <- sub(
      "\\.rds$",
      sprintf("_ch%s.rds", paste(g$fit$source_chains, collapse = "-")),
      fit_path
    )
    saveRDS(g$fit, sub_path)
    cat("  wrote chain subset:", sub_path, "\n")
  }
}
if (nzchar(args$fill_from)) {
  cat(sprintf(
    "  carried over from %s (not in the fit): %d countries\n",
    args$fill_from,
    g$n_filled
  ))
}
part <- g$partition
memb <- g$membership
write_kl_partition(part, out_csv)

dir.create(
  "generated_data/kl_groupings",
  recursive = TRUE,
  showWarnings = FALSE
)
memb_csv <- file.path(
  "generated_data/kl_groupings",
  paste0(args$out, "_membership.csv")
)
utils::write.csv(memb, memb_csv, row.names = FALSE)

print(
  as.data.frame(memb %>% mutate(across(where(is.numeric), ~ round(.x, 2)))),
  row.names = FALSE
)
cat(sprintf(
  "\nWrote %s (%d countries: %s)\n      %s (per-country probabilities)\n",
  out_csv,
  nrow(part),
  paste(
    sprintf("%s %d", names(table(part$group)), table(part$group)),
    collapse = ", "
  ),
  memb_csv
))
