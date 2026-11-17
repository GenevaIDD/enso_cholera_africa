# test_size_backend_dispatch.R — 05_fit_postprocess.R's size-track backend
#   dispatch, checked by parsing the runner's source.
#
# Usage:
#   Rscript analysis/tests/test_size_backend_dispatch.R
#
# Details:
#   Static parse, no fit. The runner reads opt[[model]]$backend and for the
#   size track routes "brms" -> compute_nested_size_re_stats_brms(), else a
#   warning + skip. Why parse instead of run: a runtime check needs a valid
#   brms fit (very slow); the regression risk is a new backend not wiring all
#   4 output kinds, or a broken "brms" string match. Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(rlang)
})

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

path <- "analysis/05_fit_postprocess.R"
if (!file.exists(path)) {
  fail(paste("missing runner:", path))
}

# Parse the runner (keep.source = TRUE so deparse gives back readable code).
parsed <- tryCatch(parse(file = path, keep.source = TRUE), error = function(e) {
  fail(paste("parse error:", e$message))
})

src <- paste(as.character(parsed), collapse = "\n")

# --- 1. reads the backend from opt[[model]]$backend -----------------------
# 1. The runner resolves the backend from opt[[model]]$backend.
# 05_fit_postprocess.R uses opt[[model]]$backend (model-generic) rather than the
# hard-coded opt$size$backend form.
if (
  !grepl(
    "opt\\[\\[model\\]\\]\\$backend|opt\\$size\\$backend|opt\\$occurrence\\$backend",
    src,
    perl = TRUE
  )
) {
  fail("05_fit_postprocess.R does not read a model-scoped backend key")
}

# --- 2. brms branch exists + calls its adapter -----------------------------
# 2. The "brms" branch exists and calls compute_nested_size_re_stats_brms().
brms_branch_re <- '(identical\\(backend,\\s*"brms"\\)|backend\\s*==\\s*"brms")'

if (!grepl(brms_branch_re, src)) {
  fail("05_fit_postprocess.R: no brms branch (identical/==)")
}

if (!grepl("compute_nested_size_re_stats_brms\\(", src)) {
  fail(
    "05_fit_postprocess.R: brms branch must call compute_nested_size_re_stats_brms()"
  )
}

# --- 3. Fixed-effect extraction is present for all backends ---------------
# 3. The fixed-effect extractor (.extract_fixed_effects) is called for both
#    backends (backend-agnostic dispatch).
# 05_fit_postprocess.R dispatches via .extract_fixed_effects(fit, backend) which
# is backend-agnostic; no separate per-backend summariseSizeCoefs* calls.
if (!grepl("\\.extract_fixed_effects\\(", src)) {
  fail(
    "05_fit_postprocess.R: missing .extract_fixed_effects() call for coef summary"
  )
}

# --- 4. Fallback branch emits a warning -----------------------------------
# 4. The fallback `else` branch emits a warning referring to 'backend'.
if (!grepl("warning\\(.*backend", src, ignore.case = TRUE)) {
  fail(
    "05_fit_postprocess.R: unknown-backend fallback must emit a warning() referring to 'backend'"
  )
}

# --- 5. All four output kinds are referenced ------------------------------
# 5. All four output kinds are referenced (size_re_stats / size_coef_summary /
#    size_loo / size_pp_check): size_pp_check literally, the others via the
#    kpfx dynamic pattern.
# size_pp_check appears as a literal string; the others are built via
# paste0(kpfx, ...) where kpfx = "size_" for the size track.
if (!grepl('"size_pp_check"', src, fixed = TRUE)) {
  fail("05_fit_postprocess.R: missing literal kind reference: size_pp_check")
}

# kpfx-based kinds: verify the four suffix tokens appear alongside kpfx
for (suffix in c("re_stats", "coef_summary", "loo")) {
  pat <- sprintf(
    'kpfx.*"%s"|"%s.*kpfx|paste0.*kpfx.*%s|paste0.*%s.*kpfx',
    suffix,
    suffix,
    suffix,
    suffix
  )
  if (!grepl(suffix, src, fixed = TRUE)) {
    fail(paste("05_fit_postprocess.R: no reference to kind suffix", suffix))
  }
}

# --- 6. Sanity: runner has access to the size RE-stat adapters ------------
# 6. The runner sources the module where the adapters live.
# The runner sources analysis/utils/siz_post.R (where
# compute_nested_size_re_stats_brms lives); the analysis/utils.R facade form is
# also accepted.
if (
  !grepl(
    'source\\("analysis/utils\\.R"\\)|source\\("analysis/utils/siz_post\\.R"\\)',
    src
  )
) {
  fail(
    "05_fit_postprocess.R: does not source analysis/utils.R or analysis/utils/siz_post.R"
  )
}

message("PASS  test_size_backend_dispatch")
