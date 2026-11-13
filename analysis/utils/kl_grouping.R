# analysis/utils/kl_grouping.R — country grouping registry from a brms_kl latent-class fit
#
# Library (sourced only). The body of utils/make_kl_grouping.R (kept as its
# CLI, same output) plus what the deploy tool needs to derive a fit registry
# from a pinned cell of a kl class sweep (deploy_prod_run.R --stage
# gen_grouping; relay phase `classes` -> `fits` in bash/tools/run_pipeline.sh):
#
#   kl_grouping_from_fit(fit, chains, fill_from)   argmax partition + membership
#   kl_grouping_sidecar(partition, meta_from)      display metadata (.yml) for it,
#                                                  incl. centroid names (KL_CENTROID_NAMES)
#   kl_grouping_status(name, partition)            absent | identical | differs + diff
#   kl_grouping_header(...)                        provenance comment lines
#   write_grouping_registry(...)                   CSV + .yml + membership CSV, re-read
#   kl_grouping_plan_from_sweep(compare_dir, ...)  everything gen_grouping prints/decides
#   kl_grouping_convergence(compare_dir, cell, ...) the pinned cell's convergence
#                                                  report / refusal
#
# Call-time dependencies (sources nothing; tests/test_utils_module_purity.R):
# base.R (%||%), config.R (load_grouping, load_grouping_meta, registry_path,
# registry_digest, resolve_group_label, .GROUP_FALLBACK_COLOURS), dlnm_kernel.R
# (kl_subset_chains, kl_responsibilities, kl_map_partition, write_kl_partition),
# which needs dlnm.R; diagnostics.R (DIAG_GATE, kl_grouping_convergence() only).
# Gates: tests/test_deploy_kl_grouping.R, tests/test_deploy_convergence_gate.R.

#' kl_grouping_from_fit
#'
#' Each fitted country in its class of maximum mean posterior responsibility
#' (kl_responsibilities() -> kl_map_partition(method = "argmax")); groups are
#' the fit's class labels (e.g. enso_east). Optionally on a chain subset
#' (kl_subset_chains()), and with the countries of registry `fill_from` that
#' the fit does not cover carried over unchanged.
#'
#' @param fit       kl_fit (04_fit_occ_brms_kl.R)
#' @param chains    integer chain ids to keep, NULL = all
#' @param fill_from registry name ("" = none)
#' @return list(fit (the chain subset or `fit`), partition = tibble(iso3, group)
#'   arranged by (group, iso3), membership = tibble(iso3, group, p_max,
#'   <class labels>), n_fit, n_filled)
kl_grouping_from_fit <- function(fit, chains = NULL, fill_from = "",
                                 groupings_dir = "analysis/groupings") {
  if (!inherits(fit, "kl_fit")) stop("kl_grouping_from_fit: not a brms_kl latent-class fit (class kl_fit).")
  if (length(chains)) fit <- kl_subset_chains(fit, chains)
  resp <- kl_responsibilities(fit)
  part <- kl_map_partition(fit, "argmax", resp)
  memb <- resp %>% dplyr::select(iso3 = country, label, mean) %>%
    tidyr::pivot_wider(names_from = label, values_from = mean) %>%
    dplyr::left_join(part, by = "iso3") %>%
    dplyr::mutate(p_max = purrr::pmap_dbl(dplyr::across(dplyr::all_of(fit$class_labels)), max)) %>%
    dplyr::relocate(iso3, group, p_max) %>% dplyr::arrange(group, dplyr::desc(p_max))
  n_fit <- nrow(part); n_filled <- 0L
  if (nzchar(fill_from %||% "")) {
    reg  <- load_grouping(fill_from, groupings_dir)
    keep <- reg[!reg$iso3 %in% part$iso3, c("iso3", "group")]
    n_filled <- nrow(keep)
    part <- dplyr::bind_rows(part, keep)
  }
  part <- part %>% dplyr::arrange(group, iso3)
  list(fit = fit, partition = part, membership = memb, n_fit = n_fit, n_filled = n_filled)
}

#' kl_grouping_sidecar
#'
#' Display metadata (utils/config.R::load_grouping_meta() format) for a
#' derived partition. Groups shared with registry `meta_from` copy its label /
#' short / short_fig / colour and keep its order; a new group is labelled by
#' its name without the prefix and takes the first fallback colour (config.R
#' .GROUP_FALLBACK_COLOURS) that `meta_from` does not use -- the rule
#' enso_grouping_kl.yml was written with by hand (enso_nga, #e6ab02).
#'
#' The display-only centroid names of results.style.group_labels "centroid"
#' (results_core.R::group_display_meta()) are carried too, per group, from the
#' first source that has them: the previous sidecar of the registry being
#' (re)written (`centroid_from`, ids matching), then `meta_from`, then
#' KL_CENTROID_NAMES (keyed by registry id). centroid_order keeps the sources'
#' order for the groups it lists, then KL_CENTROID_NAMES order, then the rest.
#' A group none of them names gets no centroid fields: centroid mode then stops
#' in group_display_meta() naming it (add it to KL_CENTROID_NAMES or the sidecar).
#' No group with centroid names -> no centroid fields at all.
#'
#' @param centroid_from registry whose existing sidecar's centroid fields are
#'   kept ("" = none; a missing sidecar is skipped)
#' @return list(label_prefix, order, [centroid_order,] groups = list(<g> =
#'   list(label, short, short_fig, colour[, centroid_short, centroid_label])))
kl_grouping_sidecar <- function(partition, meta_from = "enso_grouping",
                                groupings_dir = "analysis/groupings", centroid_from = "") {
  groups <- sort(unique(partition$group))
  meta   <- if (nzchar(meta_from %||% "")) load_grouping_meta(meta_from, groupings_dir) else NULL
  prefix <- meta$label_prefix %||% ""
  if (nzchar(prefix) && !all(startsWith(groups, prefix))) prefix <- ""
  shared <- intersect(meta$order %||% character(), groups)
  ord    <- c(shared, setdiff(groups, shared))
  free   <- setdiff(.GROUP_FALLBACK_COLOURS, unname(meta$colour %||% character()))
  side <- list(label_prefix = prefix, order = as.list(ord), groups = list())
  for (g in ord) {
    if (g %in% shared) {
      side$groups[[g]] <- list(label = unname(meta$label[g]), short = unname(meta$short[g]),
                               short_fig = unname(meta$short_fig[g]), colour = unname(meta$colour[g]))
    } else {
      nm <- if (nzchar(prefix)) sub(paste0("^", prefix), "", g) else g
      col <- if (length(free)) free[1] else .GROUP_FALLBACK_COLOURS[1]
      free <- free[-1]
      side$groups[[g]] <- list(label = nm, short = nm, short_fig = nm, colour = col)
    }
  }
  cen <- kl_centroid_fields(ord, centroid_from, meta, groupings_dir)
  if (length(cen$order)) {
    for (g in names(cen$short))
      side$groups[[g]] <- c(side$groups[[g]], list(centroid_short = cen$short[[g]],
                                                    centroid_label = cen$label[[g]]))
    side <- c(side[c("label_prefix", "order")], list(centroid_order = as.list(cen$order)),
              side["groups"])
  }
  side
}

#' Display-only centroid names by registry id (results.style.group_labels
#' "centroid"; supplement grouping section): the fallback of
#' kl_grouping_sidecar() for a group no previous sidecar names. Row order is
#' the centroid display order (west to east, then south). Ids never change.
KL_CENTROID_NAMES <- data.frame(
  id    = c("enso_west", "enso_nga",     "enso_horn", "enso_east", "enso_south"),
  short = c("GW",        "GG",           "GH",        "GE",        "GS"),
  label = c("West",      "Guinean Gulf", "Horn",      "East",      "South"),
  stringsAsFactors = FALSE)

#' kl_centroid_fields
#'
#' Centroid names of groups `ids` for kl_grouping_sidecar(): per group from the
#' sidecar of `centroid_from` (when it exists), then `meta` (load_grouping_meta()
#' of meta_from), then KL_CENTROID_NAMES.
#' @return list(short, label (named character, groups that have both), order
#'   (every id of `ids`, or character(0) when no group has names))
kl_centroid_fields <- function(ids, centroid_from = "", meta = NULL,
                               groupings_dir = "analysis/groupings") {
  srcs <- list()
  if (nzchar(centroid_from %||% "") && file.exists(registry_path(centroid_from, groupings_dir, "yml"))) {
    side <- yaml::read_yaml(registry_path(centroid_from, groupings_dir, "yml"))
    g <- names(side$groups)
    pick <- function(f) stats::setNames(vapply(g, function(x) as.character(side$groups[[x]][[f]] %||% NA_character_),
                                               character(1)), g)
    srcs[[length(srcs) + 1L]] <- list(short = pick("centroid_short"), label = pick("centroid_label"),
                                      order = as.character(unlist(side[["centroid_order"]])))
  }
  if (!is.null(meta))
    srcs[[length(srcs) + 1L]] <- list(short = meta[["centroid_short"]], label = meta[["centroid_label"]],
                                      order = as.character(meta[["centroid_order"]] %||% character(0)))
  srcs[[length(srcs) + 1L]] <- list(short = stats::setNames(KL_CENTROID_NAMES$short, KL_CENTROID_NAMES$id),
                                    label = stats::setNames(KL_CENTROID_NAMES$label, KL_CENTROID_NAMES$id),
                                    order = KL_CENTROID_NAMES$id)
  ok <- function(v) !is.null(v) && !is.na(v) && nzchar(v)
  short <- label <- character()
  for (g in ids) for (s in srcs) {
    cs <- unname(s$short[g]); cl <- unname(s$label[g])
    if (length(cs) && length(cl) && ok(cs) && ok(cl)) { short[[g]] <- cs; label[[g]] <- cl; break }
  }
  if (!length(short)) return(list(short = short, label = label, order = character(0)))
  ord <- character()
  for (s in srcs) ord <- c(ord, setdiff(intersect(s$order, ids), ord))
  ord <- c(ord, setdiff(ids, ord))
  list(short = short, label = label, order = ord)
}

#' kl_grouping_status
#'
#' Compare a derived partition with registry `name` on disk (iso3 -> group
#' map, row order ignored).
#' @return list(status = "absent" | "identical" | "differs",
#'   diff = data.frame(iso3, old, new, change in changed | added | removed))
kl_grouping_status <- function(name, partition, groupings_dir = "analysis/groupings") {
  empty <- data.frame(iso3 = character(), old = character(), new = character(),
                      change = character(), stringsAsFactors = FALSE)
  if (!file.exists(registry_path(name, groupings_dir)))
    return(list(status = "absent", diff = empty))
  old <- load_grouping(name, groupings_dir)
  m <- merge(data.frame(iso3 = old$iso3, old = old$group, stringsAsFactors = FALSE),
             data.frame(iso3 = partition$iso3, new = partition$group, stringsAsFactors = FALSE),
             by = "iso3", all = TRUE)
  m$change <- ifelse(is.na(m$old), "added", ifelse(is.na(m$new), "removed",
                     ifelse(m$old == m$new, "same", "changed")))
  d <- m[m$change != "same", c("iso3", "old", "new", "change"), drop = FALSE]
  rownames(d) <- NULL
  list(status = if (nrow(d)) "differs" else "identical", diff = if (nrow(d)) d else empty)
}

#' kl_grouping_header
#'
#' Provenance comment lines for the derived registry's sidecar: what the
#' hand-written enso_grouping_kl.yml records by hand.
kl_grouping_header <- function(name, deploy_config, batch, cell, chains, fit_path,
                               settings = list(), n_fit = NA, n_filled = 0L, fill_from = "",
                               meta_from = "") {
  sha <- if (file.exists(fit_path) && requireNamespace("digest", quietly = TRUE))
    digest::digest(file = fit_path, algo = "sha256") else NA_character_
  c(sprintf("# Display metadata for %s.csv (utils/config.R::load_grouping_meta()).", name),
    "# Generated by analysis/utils/deploy_prod_run.R --stage gen_grouping",
    "# (utils/kl_grouping.R): each fitted country in its class of maximum mean",
    "# posterior probability. Do not edit by hand; re-run the stage.",
    sprintf("#   deploy config : %s", deploy_config),
    sprintf("#   sweep batch   : %s (cell %s)", batch, cell),
    sprintf("#   chains        : %s", if (length(chains)) paste(chains, collapse = ",") else "all"),
    sprintf("#   fit           : %s", fit_path),
    sprintf("#   fit sha256    : %s", sha),
    sprintf("#   compare       : mode_tol %s, n_agree %s",
            settings$mode_tol %||% NA, settings$n_agree %||% NA),
    sprintf("#   countries     : %s fitted%s", n_fit,
            if (n_filled > 0L) sprintf(" + %d carried over from %s", n_filled, fill_from) else ""),
    sprintf("#   labels from   : %s", if (nzchar(meta_from %||% "")) meta_from else "(none)"),
    sprintf("#   written       : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

#' write_grouping_registry
#'
#' Write `<groupings_dir>/<name>.csv` (dlnm_kernel.R::write_kl_partition(),
#' the make_kl_grouping.R format) and its `.yml` sidecar, both read back
#' (load_grouping() / load_grouping_meta()), plus the per-country
#' probabilities to `<membership_dir>/<name>_membership.csv`.
#' @return list(csv, yml, membership_csv, digest)
write_grouping_registry <- function(partition, name, sidecar, header = character(),
                                    membership = NULL, groupings_dir = "analysis/groupings",
                                    membership_dir = "generated_data/kl_groupings") {
  csv <- registry_path(name, groupings_dir)
  yml <- registry_path(name, groupings_dir, "yml")
  write_kl_partition(partition, csv)
  writeLines(c(header, yaml::as.yaml(sidecar)), yml)
  load_grouping_meta(name, groupings_dir)
  memb_csv <- NA_character_
  if (!is.null(membership)) {
    dir.create(membership_dir, recursive = TRUE, showWarnings = FALSE)
    memb_csv <- file.path(membership_dir, paste0(name, "_membership.csv"))
    utils::write.csv(membership, memb_csv, row.names = FALSE)
  }
  list(csv = csv, yml = yml, membership_csv = memb_csv,
       digest = registry_digest(name, groupings_dir))
}

#' kl_grouping_plan_from_sweep
#'
#' Everything `deploy_prod_run.R --stage gen_grouping` prints and decides,
#' writing nothing: reads the compare outputs of 05_fit_postprocess_kl_sweep.R
#' (`loo.rds` holds the inventory `inv` and the agreeing chains `agree` per fit
#' id), loads the pinned cell's fit, derives the partition and its sidecar,
#' compares it with registry `name` on disk, and checks the deploy's
#' `fits.groups` against the derived groups.
#'
#' @param compare_dir generated_data/kl_sweep_compare/<batch>/
#' @param cell        pinned fit id (kl_sweep.R::kl_sweep_inventory() `id`, e.g. "ns5_K4")
#' @param chains      "agree" | "all" | "2,6"
#' @param n_agree,mode_tol the chain-selection settings the comparison must
#'   have been run with when chains = "agree" (NULL = not checked)
#' @param fits_groups deploy.fits.groups (short forms or full labels)
#' @return list(cell, fit_path, chains, settings, grouping, sidecar, status,
#'   diff, groups_ok, groups_short, fits_groups)
kl_grouping_plan_from_sweep <- function(compare_dir, cell, chains = "agree", fill_from = "",
                                        meta_from = "enso_grouping", name, fits_groups,
                                        groupings_dir = "analysis/groupings",
                                        fit_cmd = "bash analysis/bash/stages/04_fit_kl.sh --deploy <D> --postprocess",
                                        n_agree = NULL, mode_tol = NULL) {
  need <- file.path(compare_dir, c("fits.csv", "partition.csv", "loo.rds"))
  if (!all(file.exists(need)))
    stop("kl_grouping: sweep comparison incomplete in ", compare_dir, " (missing ",
         paste(basename(need[!file.exists(need)]), collapse = ", "), "). Fit + compare the sweep first: ",
         fit_cmd, call. = FALSE)
  loo <- readRDS(file.path(compare_dir, "loo.rds"))
  inv <- loo$inv
  ids <- inv$id
  if (!nzchar(cell %||% ""))
    stop("kl_grouping: deploy.kl_sweep.cell is empty. Pin one of: ", paste(ids, collapse = ", "),
         " (see ", file.path(compare_dir, "fits.csv"), " and loo.csv).", call. = FALSE)
  if (!cell %in% ids)
    stop("kl_grouping: cell '", cell, "' is not in the sweep (", paste(ids, collapse = ", "), ").",
         call. = FALSE)
  fit_path <- inv$fit_file[inv$id == cell]
  if (!file.exists(fit_path)) stop("kl_grouping: fit of cell '", cell, "' not found: ", fit_path, call. = FALSE)
  if (identical(chains %||% "agree", "agree") && (!is.null(n_agree) || !is.null(mode_tol))) {
    st  <- loo$settings %||% list()
    bad <- c(if (!is.null(n_agree) && !identical(as.integer(st$n_agree %||% NA), as.integer(n_agree)))
               sprintf("n_agree %s (deploy %s)", st$n_agree %||% "unset", n_agree),
             if (!is.null(mode_tol) && !isTRUE(all.equal(as.numeric(st$mode_tol %||% NA), as.numeric(mode_tol))))
               sprintf("mode_tol %s (deploy %s)", st$mode_tol %||% "unset", mode_tol))
    if (length(bad))
      stop("kl_grouping: the comparison in ", compare_dir, " picked its agreeing chains with ",
           paste(bad, collapse = ", "), ". Re-run it with the deploy's settings: Rscript ",
           "analysis/05_fit_postprocess_kl_sweep.R -d ", paste(st$batch_dirs %||% "<sweep dir>", collapse = ","),
           " --n_agree ", n_agree %||% 0L, " --mode_tol ", mode_tol %||% 1.5,
           " (or rerun ", fit_cmd, ": finished fits are skipped).", call. = FALSE)
  }
  ch <- switch(chains %||% "agree",
    agree = {
      a <- loo$agree[[cell]]
      if (is.null(a)) stop("kl_grouping: no agreeing chains recorded for '", cell, "' in loo.rds.", call. = FALSE)
      as.integer(a)
    },
    all = NULL,
    as.integer(strsplit(chains, ",")[[1]]))
  if (length(ch) && anyNA(ch)) stop("kl_grouping: kl_sweep.chains must be agree, all or ids like 2,6 (got '", chains, "').", call. = FALSE)
  fit <- readRDS(fit_path)
  if (length(ch) && setequal(ch, unique(fit$draws$.chain))) ch <- NULL   # every chain agrees: no subset
  grp  <- kl_grouping_from_fit(fit, ch, fill_from, groupings_dir)
  side <- kl_grouping_sidecar(grp$partition, meta_from, groupings_dir, centroid_from = name)
  st   <- kl_grouping_status(name, grp$partition, groupings_dir)
  labels <- unlist(side$order)
  full <- ifelse(fits_groups %in% labels, fits_groups, paste0(side$label_prefix, fits_groups))
  short <- if (nzchar(side$label_prefix)) sub(paste0("^", side$label_prefix), "", labels) else labels
  list(cell = cell, fit_path = fit_path, chains = ch, settings = loo$settings,
       grouping = grp, sidecar = side, status = st$status, diff = st$diff,
       groups_ok = setequal(full, labels) && length(full) == length(labels),
       groups_short = short, fits_groups = fits_groups)
}

#' Minimum number of agreeing chains gen_grouping derives a registry from
#' under a refusing convergence gate.
KL_GROUPING_MIN_AGREE <- 2L

#' kl_grouping_convergence
#'
#' Convergence report of the pinned sweep cell for `deploy_prod_run.R --stage
#' gen_grouping`, read from the 05
#' comparison only (`fits.csv`, columns of kl_sweep.R::kl_fit_row(), and
#' `diagnostics_chains.csv` when present): R-hat / ESS on all chains and on the
#' agreeing chains, divergences / max-treedepth hits / E-BFMI on all chains,
#' on the agreeing chains and on the chain subset the deploy derives the
#' registry from (`chains`, when it is another set). R-hat on chains chosen
#' BECAUSE they agree is optimistic and cannot certify a multimodal posterior;
#' membership uncertainty is report_membership's job. So the gate is coarse:
#' under cfg$action `refuse` it stops when fewer than KL_GROUPING_MIN_AGREE
#' chains agree or the agreeing-chain R-hat exceeds DIAG_GATE[["rhat_any"]]
#' (1.05); under `warn` those conditions only warn; beyond cfg$rhat_max /
#' cfg$ess_min / cfg$max_divergent / DIAG_GATE[["ebfmi"]] it warns under
#' either action. `off` (the default) returns NULL and prints nothing.
#' `dry_run` prints "WOULD refuse" and never stops.
#'
#' Call-time dependencies: diagnostics.R (DIAG_GATE); `cfg` is
#' deploy_core.R::deploy_convergence_cfg().
#'
#' @param compare_dir generated_data/kl_sweep_compare/<batch>/
#' @param cell   pinned fit id
#' @param chains the chain ids the registry is derived from (NULL = all)
#' @param cfg    list(action, rhat_max, ess_min, max_divergent)
#' @return NULL (off) or list(table (one row per chain set: all, agree[,
#'   deploy]), n_agree, agreeing, refuse (reasons), warn (reasons),
#'   verdict = "pass" | "warn" | "refuse")
kl_grouping_convergence <- function(compare_dir, cell, chains = NULL, cfg,
                                    dry_run = FALSE, say = function(...) invisible(NULL)) {
  if (identical(cfg$action %||% "off", "off")) return(NULL)
  fits <- utils::read.csv(file.path(compare_dir, "fits.csv"), stringsAsFactors = FALSE)
  r <- fits[fits$id == cell, , drop = FALSE]
  if (nrow(r) != 1L)
    stop("kl_grouping_convergence: cell '", cell, "' has ", nrow(r), " row(s) in ",
         file.path(compare_dir, "fits.csv"), " (expected 1).", call. = FALSE)
  agreeing <- suppressWarnings(as.integer(strsplit(as.character(r$agreeing %||% ""), ",")[[1]]))
  agreeing <- agreeing[!is.na(agreeing)]
  dc_path <- file.path(compare_dir, "diagnostics_chains.csv")
  dc <- if (file.exists(dc_path)) {
    x <- utils::read.csv(dc_path, stringsAsFactors = FALSE)
    x[x$id == cell, , drop = FALSE]
  } else NULL
  samp <- function(ch) {
    if (is.null(dc) || !nrow(dc))
      return(c(divergences = NA_real_, treedepth_hits = NA_real_, min_ebfmi = NA_real_))
    s <- if (is.null(ch)) dc else dc[dc$chain %in% ch, , drop = FALSE]
    c(divergences = sum(s$divergent), treedepth_hits = sum(s$treedepth_hits),
      min_ebfmi = if (nrow(s)) min(s$ebfmi) else NA_real_)
  }
  row_of <- function(set, ch, rhat, essb, esst, n_ch, smp) data.frame(
    chains = set, ids = if (is.null(ch)) "all" else paste(ch, collapse = ","), n_chains = n_ch,
    max_rhat = rhat, min_ess_bulk = essb, min_ess_tail = esst,
    divergences = smp[["divergences"]], treedepth_hits = smp[["treedepth_hits"]],
    min_ebfmi = smp[["min_ebfmi"]], stringsAsFactors = FALSE)
  s_all <- if (is.null(dc) || !nrow(dc))
    c(divergences = r$divergences, treedepth_hits = r$treedepth_hits, min_ebfmi = r$min_ebfmi)
  else samp(NULL)
  tab <- rbind(
    row_of("all", NULL, r$max_rhat_all, r$min_ess_bulk_all, r$min_ess_tail_all, r$chains, s_all),
    row_of("agree", agreeing, r$max_rhat_agree, r$min_ess_bulk_agree, r$min_ess_tail_agree,
           length(agreeing), samp(agreeing)))
  if (length(chains) && !setequal(as.integer(chains), agreeing))
    tab <- rbind(tab, row_of("deploy", as.integer(chains), NA_real_, NA_real_, NA_real_,
                             length(chains), samp(as.integer(chains))))
  ag <- tab[tab$chains == "agree", ]
  refuse <- c(
    if (length(agreeing) < KL_GROUPING_MIN_AGREE)
      sprintf("%d chain(s) agree (< %d)", length(agreeing), KL_GROUPING_MIN_AGREE),
    if (isTRUE(ag$max_rhat > DIAG_GATE[["rhat_any"]]))
      sprintf("agreeing-chain R-hat %.4f > %s", ag$max_rhat, format(DIAG_GATE[["rhat_any"]])),
    if (length(agreeing) >= KL_GROUPING_MIN_AGREE && is.na(ag$max_rhat))
      "agreeing-chain R-hat missing from fits.csv")
  warn <- c(
    if (isTRUE(ag$max_rhat > cfg$rhat_max))
      sprintf("agreeing-chain R-hat %.4f > %s", ag$max_rhat, format(cfg$rhat_max)),
    if (isTRUE(ag$min_ess_bulk < cfg$ess_min))
      sprintf("agreeing-chain bulk ESS %.0f < %s", ag$min_ess_bulk, format(cfg$ess_min)),
    if (isTRUE(ag$min_ess_tail < cfg$ess_min))
      sprintf("agreeing-chain tail ESS %.0f < %s", ag$min_ess_tail, format(cfg$ess_min)),
    if (isTRUE(ag$divergences > cfg$max_divergent))
      sprintf("%d divergence(s) on the agreeing chains", as.integer(ag$divergences)),
    if (isTRUE(ag$min_ebfmi < DIAG_GATE[["ebfmi"]]))
      sprintf("agreeing-chain E-BFMI %.2f < %s", ag$min_ebfmi, format(DIAG_GATE[["ebfmi"]])))
  say("[deploy]   convergence (action %s; cell %s, %s mode(s) %s):", cfg$action, cell,
      r$n_modes %||% NA, r$mode_sizes %||% "")
  for (i in seq_len(nrow(tab)))
    say("[deploy]     %-6s chains %-12s rhat=%s ess_bulk=%s ess_tail=%s div=%s td=%s ebfmi=%s",
        tab$chains[[i]], tab$ids[[i]], formatC(tab$max_rhat[[i]], format = "f", digits = 4),
        formatC(tab$min_ess_bulk[[i]], format = "f", digits = 0),
        formatC(tab$min_ess_tail[[i]], format = "f", digits = 0),
        tab$divergences[[i]], tab$treedepth_hits[[i]],
        formatC(tab$min_ebfmi[[i]], format = "f", digits = 2))
  if (any(tab$chains == "deploy"))
    say("[deploy]     (the registry is derived from chains %s; their R-hat / ESS are not in the comparison)",
        tab$ids[tab$chains == "deploy"])
  verdict <- if (length(refuse) && identical(cfg$action, "refuse")) "refuse"
             else if (length(refuse) || length(warn)) "warn" else "pass"
  if (length(warn))
    warning("[deploy] gen_grouping: cell ", cell, ": ", paste(warn, collapse = "; "), call. = FALSE)
  if (length(refuse)) {
    msg <- sprintf("cell %s fails the grouping convergence gate: %s.", cell,
                   paste(refuse, collapse = "; "))
    if (identical(cfg$action, "refuse")) {
      if (dry_run) say("[deploy] WOULD refuse: %s", msg)
      else stop("[deploy] gen_grouping: refusing to derive the registry: ", msg,
                " Pin another cell (deploy.kl_sweep.cell) or refit the sweep; nothing written.",
                call. = FALSE)
    } else {
      warning("[deploy] gen_grouping: ", msg, call. = FALSE)
    }
  }
  list(table = tab, n_agree = length(agreeing), agreeing = agreeing,
       refuse = refuse, warn = warn, verdict = verdict)
}
