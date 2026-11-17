# figures_io.R — shared figure-saving helpers: one saver for the figure
#   scripts instead of per-script copies.
#
# Provides:
#   write_figure_files()   write one plot to one or more formats under a dir
#   figure_saver()         a saver closure bound to a directory
#   true_minus(), ps_true_minus()   typographic minus in labels
#
# Details:
#   Filenames are left to the callers; directories come from
#   paths.R::make_figure_dir(), which callers call (this file only writes).
#   Callers source base.R, config.R and paths.R; ggplot2 via `::`.

# Deliberately NOT an extension of make_run_file(): its `kind` is a closed
# match.arg() set, and the figure scripts have dozens of distinct stems --
# adding figure kinds there would mean a huge match.arg() plus a paths.R edit
# per new figure. Path algebra (make_figure_dir(), paths.R) stays separate
# from writing (write_figure_files(), here) so paths.R never needs a ggplot2
# dependency.
#' write_figure_files
#'
#' Write one plot to one or more formats under `dir`, creating `dir` if
#' needed. Replaces the ~15 repeated `dir.create(recursive = TRUE)` +
#' `ggsave()` + `message("Saved: ...")` sequences across the figure scripts.
#'
#' @param plot A ggplot object (or anything `ggplot2::ggsave()` accepts).
#' @param stem File basename, without extension or leading `prefix`.
#' @param dir Output directory. Created (recursively, quietly) if absent.
#' @param width,height Passed through to `ggplot2::ggsave()`.
#' @param units Passed through to `ggplot2::ggsave()`. Default `"in"`; pass
#'   `"mm"` for 10_manuscript_figures.R's millimetre-based figures.
#' @param dpi Passed through to `ggplot2::ggsave()`. Ignored for vector
#'   formats (pdf) exactly as `ggsave()` itself ignores it.
#' @param formats Character vector of file extensions to write, e.g. `"png"`
#'   or `c("pdf", "png")`. One `ggsave()` call per format.
#' @param prefix Prepended to `stem` verbatim (no separator inserted) -- pass
#'   `"<tag>_"` to reproduce 06_fit_figures.R's `<tag>_<stem>.png` naming.
#' @param device Optional explicit device (e.g. `grDevices::cairo_pdf`,
#'   `ragg::agg_png`), passed through to every format's `ggsave()` call.
#'   NULL lets `ggsave()` infer the device from the extension, as it does
#'   for every existing bare `ggsave()` call site being migrated.
#' @param ... Forwarded verbatim to every `ggsave()` call -- e.g.
#'   `limitsize = FALSE` for the wide multi-country facets in
#'   `06_fit_figures_lag_compare.R`.
#' @return Invisibly, the character vector of paths written.
write_figure_files <- function(
  plot,
  stem,
  dir,
  width,
  height,
  units = "in",
  dpi = 150,
  formats = "png",
  prefix = "",
  device = NULL,
  ...
) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- vapply(
    formats,
    function(fmt) {
      path <- file.path(dir, paste0(prefix, stem, ".", fmt))
      ggplot2::ggsave(
        path,
        plot,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        device = device,
        ...
      )
      message("Saved: ", path)
      path
    },
    character(1)
  )
  invisible(unname(paths))
}


#' figure_saver
#'
#' Closure factory over `write_figure_files()`, signature-compatible with the four
#' pre-existing local `save_fig()` closures it replaces
#' (06_fit_figures_cb.R, _cb_adm1.R, _lag_compare.R; 10_manuscript_figures.R's
#' `save_figure()` is kept as its own thin wrapper -- see that file) so each
#' call site becomes `save <- figure_saver(dir, ...)` followed by unchanged
#' `save(plot, name)` calls.
#'
#' @param dir Output directory, fixed for every call the returned closure
#'   makes (matches the pre-existing `save_fig()` closures, which capture
#'   `fig_dir` from their enclosing scope).
#' @param width,height,dpi,formats As in `write_figure_files()`.
#' @param prefix As in `write_figure_files()` -- pass a batch tag to reproduce
#'   06_fit_figures.R's `<tag>_<stem>.png` naming.
#' @return A function `function(plot, stem, width = width, height = height,
#'   ...)` -- per-call `width`/`height` override the closure's defaults,
#'   matching the pre-existing `save_fig(p, name, width = 8, height = 5)`
#'   call shape.
figure_saver <- function(
  dir,
  width = 8,
  height = 5,
  dpi = 150,
  formats = "png",
  prefix = ""
) {
  # Captured under different names before being used as inner defaults --
  # `function(width = width, ...)` would make the inner formal's default
  # expression refer to itself (R evaluates defaults in the call's own
  # execution environment, which the formal shadows), raising "promise
  # already under evaluation" the first time a caller omits `width`.
  default_width <- width
  default_height <- height
  function(plot, stem, width = default_width, height = default_height, ...) {
    # modifyList(), not a bare `..., dpi = dpi, ...` splice: a caller passing
    # e.g. dpi = 300 in `...` would otherwise land as a SECOND `dpi=` actual
    # alongside the closure's own `dpi = dpi`, and write_figure_files() would raise
    # "formal argument 'dpi' matched by multiple actual arguments". This way
    # a `...` override cleanly replaces the closure's default instead.
    args <- utils::modifyList(
      list(dir = dir, dpi = dpi, formats = formats, prefix = prefix),
      list(...)
    )
    do.call(
      write_figure_files,
      c(list(plot = plot, stem = stem, width = width, height = height), args)
    )
  }
}

# ---------------------------------------------------------------------------
# True minus signs (Figures 1-4)
# ---------------------------------------------------------------------------
#' Every ASCII hyphen-minus right before a digit or a decimal point -> the true
#' minus sign (U+2212), e.g. "-0.5" -> "−0.5"; other hyphens are kept.
true_minus <- function(x) gsub("-(?=[0-9.])", "−", x, perl = TRUE)

#' A ggplot whose continuous scales (axes and colourbars) label negatives with
#' the true minus: explicit labels are rewritten, default ones are formatted by
#' the scale's transformation first. The scales are cloned, so the builder's
#' own scale objects are untouched.
ps_true_minus <- function(p) {
  p$scales$scales <- lapply(p$scales$scales, function(sc) {
    if (!inherits(sc, "ScaleContinuous")) {
      return(sc)
    }
    sc <- sc$clone()
    old <- sc$labels
    sc$labels <-
      if (is.function(old)) {
        function(x) true_minus(old(x))
      } else if (inherits(old, "waiver")) {
        function(x) true_minus(sc$get_transformation()$format(x))
      } else if (is.null(old)) {
        NULL
      } else {
        true_minus(as.character(old))
      }
    sc
  })
  p
}
