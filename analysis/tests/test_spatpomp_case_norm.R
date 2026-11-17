# test_spatpomp_case_norm.R — HARD: the normalized inactive-state case term of
#   the `_actn` measurements (..._state_detect_actn and its `_who` form).
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_case_norm.R
#
# Details:
#   Pure C/R, run from the project root; needs a C compiler (R CMD SHLIB). No
#   fits, no MCMC. N1 the `_act` siblings emit byte-identical C against
#   golden/spatpomp_case_act_golden.rds (side-by-side rule); N2 the
#   inactive-state term has mass 1; N3 active branches; N4 enum / dispatch /
#   config.

source("analysis/utils/base.R")
source("analysis/utils/spatpomp_core.R")
source("analysis/tests/helpers.R")
.passes <- 0L
.fails <- 0L

cat("# test_spatpomp_case_norm -- inactive-state case term\n\n")

ACT <- "occurrence_duration_mon_NB_state_detect_act"
ACTW <- "occurrence_duration_mon_NB_state_detect_act_who"
ACTN <- "occurrence_duration_mon_NB_state_detect_actn"
ACTNW <- "occurrence_duration_mon_NB_state_detect_actn_who"
WHO_ARGS <- list(
  who_w = c(0.25, 0.75),
  who_s = c(2.5, 1 / 3),
  phi = 5,
  mu_floor = 1
)

# ---------------------------------------------------------------------------
# N1: the case terms and the full `_act` / `_act_who` dunit text equal the
#     golden (golden/spatpomp_case_act_golden.rds).
# REGENERATING THE GOLDEN: only if the `_act` C is meant to change, with the
# module version that defines the new reference:
#   g <- list(case_act = .sp_case_actpen_zi_c(0.1), case_act_dfn25 = .sp_case_actpen_zi_c(0.25),
#             case_who = <.sp_case_who_c(c(0.25, 0.75), c(2.5, 1/3), 5, 1)>,
#             dunit_act = as.character(.select_dunit("..._act", 0.1, 0.001)),
#             dunit_act_who = as.character(.select_dunit("..._act_who", 0.1, 0.001, who_c = case_who)))
cat("## N1: `_act` siblings unchanged\n")
g <- readRDS("analysis/tests/golden/spatpomp_case_act_golden.rds")
who <- do.call(.sp_case_who_c, WHO_ARGS)
expect_true(
  identical(.sp_case_actpen_zi_c(0.1), g$case_act) &&
    identical(.sp_case_actpen_zi_c(0.25), g$case_act_dfn25),
  "N1: .sp_case_actpen_zi_c() emits the golden C"
)
expect_true(
  identical(who, g$case_who),
  "N1: .sp_case_who_c() emits the golden C"
)
expect_true(
  identical(as.character(.select_dunit(ACT, 0.1, 0.001)), g$dunit_act),
  "N1: the `_act` dunit is byte-identical"
)
expect_true(
  identical(
    as.character(.select_dunit(ACTW, 0.1, 0.001, who_c = who)),
    g$dunit_act_who
  ),
  "N1: the `_act_who` dunit is byte-identical"
)
expect_true(
  identical(
    as.character(.select_dunit(
      ACT,
      0.1,
      0.001,
      case_fp_rate = 0.3,
      case_fp_mean = 2
    )),
    g$dunit_act
  ),
  "N1: case_fp_rate / case_fp_mean are inert for `_act`"
)

# ---------------------------------------------------------------------------
# Compile the emitted case terms as plain C functions: f(x, mu, phi, ccm) -> ll.
compile_terms <- function(terms) {
  d <- tempfile("casenorm_")
  dir.create(d)
  body <- vapply(
    names(terms),
    function(nm) {
      sprintf(
        "void %s(double *x, double *mu_in, double *phi_in, double *ccm_in, double *out) {
       double obs_CM = *x, obs_W = *x, mu = *mu_in, size_phi = *phi_in, ccm = *ccm_in;
       int u = 0; double ll = 0.0;
       (void)obs_CM; (void)obs_W; (void)mu; (void)size_phi; (void)ccm; (void)u;
       {
%s
       }
       *out = ll;
     }",
        nm,
        terms[[nm]]
      )
    },
    character(1)
  )
  src <- file.path(d, "casenorm.c")
  writeLines(
    c("#include <R.h>", "#include <Rmath.h>", "#include <math.h>", body),
    src
  )
  out <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", shQuote(src)),
    stdout = TRUE,
    stderr = TRUE
  )
  so <- sub("\\.c$", .Platform$dynlib.ext, src)
  if (!file.exists(so)) {
    stop("compile failed:\n", paste(out, collapse = "\n"))
  }
  dyn.load(so)
  function(nm, x, mu = 0, phi = 1, ccm = 0) {
    vapply(
      x,
      function(xi) {
        .C(
          nm,
          as.double(xi),
          as.double(mu),
          as.double(phi),
          as.double(ccm),
          out = double(1)
        )$out
      },
      numeric(1)
    )
  }
}
EPS <- c(1, 0.5, 0.01)
PHI <- c(0.2, 0.4, 2)
MU0 <- 0.1
terms <- c(
  list(
    act = .sp_case_actpen_zi_c(0.1),
    who = do.call(.sp_case_who_c, WHO_ARGS)
  ),
  setNames(
    lapply(EPS, function(e) .sp_case_actpen_zi_norm_c(0.1, e, MU0)),
    sprintf("actn_%d", seq_along(EPS))
  ),
  setNames(
    lapply(EPS, function(e) {
      do.call(.sp_case_who_norm_c, c(WHO_ARGS, fp_rate = e))
    }),
    sprintf("whon_%d", seq_along(EPS))
  )
)
f <- compile_terms(terms)
K <- 0:400

# ---------------------------------------------------------------------------
# N2: normalization: the emitted inactive-state C (mu == 0), compiled
#     standalone and summed over counts 0..400, has mass 1 within 1e-6 for phi
#     in {0.2, 0.4, 2} and eps in {1, 0.5, 0.01}; the old `_act` term has mass
#     2 - NB(0; mu_floor, phi) (documents the defect). Same for the WHO term
#     (ccm == 0).
cat("\n## N2: total mass at mu = 0 (unit) and ccm = 0 (WHO)\n")
for (phi in PHI) {
  nb0 <- dnbinom(0, size = phi, mu = MU0)
  m_old <- sum(exp(f("act", K, mu = 0, phi = phi)))
  expect_true(
    abs(m_old - (2 - nb0)) < 1e-6,
    sprintf(
      "N2: `_act` mass at phi %.1f is 2 - NB(0) = %.6f (got %.6f)",
      phi,
      2 - nb0,
      m_old
    )
  )
  for (i in seq_along(EPS)) {
    m <- sum(exp(f(sprintf("actn_%d", i), K, mu = 0, phi = phi)))
    expect_true(
      abs(m - 1) < 1e-6,
      sprintf(
        "N2: `_actn` mass at phi %.1f, eps %.2f is 1 (got %.9f)",
        phi,
        EPS[i],
        m
      )
    )
  }
}
# WHO term: tempered by w_u (unit 0, w = 0.25), so normalize the untempered
# density: exp(ll / w_u). Mass over 0..5000 at mu_floor 1, phi 5.
KW <- 0:5000
wu <- WHO_ARGS$who_w[1]
m_who <- sum(exp(f("who", KW, ccm = 0) / wu))
expect_true(
  abs(m_who - 1) < 1e-6,
  sprintf("N2: `_who` floor branch already has mass 1 (got %.9f)", m_who)
)
for (i in seq_along(EPS)) {
  m <- sum(exp(f(sprintf("whon_%d", i), KW, ccm = 0) / wu))
  expect_true(
    abs(m - 1) < 1e-6,
    sprintf(
      "N2: `_actn_who` country-inactive mass at fp_rate %.2f is 1 (got %.9f)",
      EPS[i],
      m
    )
  )
}

# ---------------------------------------------------------------------------
# N3: active branches of `_actn` equal `_act`'s; eps = 1 equals the plain
#     floor NB; the `_actn_who` term at fp_rate = 1 scores what `_who` scores.
cat("\n## N3: branch values\n")
kk <- c(0, 1, 3, 12, 50.5)
for (phi in PHI) {
  expect_true(
    max(abs(
      f("actn_1", kk, mu = 40, phi = phi) - f("act", kk, mu = 40, phi = phi)
    )) <
      1e-12,
    sprintf("N3: active branches of `_actn` equal `_act` (phi %.1f)", phi)
  )
  ki <- c(0, 1, 3, 12)
  expect_true(
    max(abs(
      f("actn_1", ki, mu = 0, phi = phi) -
        dnbinom(ki, size = phi, mu = MU0, log = TRUE)
    )) <
      1e-9,
    sprintf("N3: eps = 1 is the plain floor NB at every count (phi %.1f)", phi)
  )
  e <- EPS[3]
  want <- c(
    log((1 - e) + e * dnbinom(0, size = phi, mu = MU0)),
    log(e) + dnbinom(c(1, 3, 12), size = phi, mu = MU0, log = TRUE)
  )
  expect_true(
    max(abs(f("actn_3", c(0, 1, 3, 12), mu = 0, phi = phi) - want)) < 1e-9,
    sprintf("N3: eps = 0.01 mixture branches (phi %.1f)", phi)
  )
}
kw <- c(0, 1, 7, 300)
expect_true(
  max(abs(f("whon_1", kw, ccm = 0) - f("who", kw, ccm = 0))) < 1e-12 &&
    max(abs(f("whon_3", kw, ccm = 12) - f("who", kw, ccm = 12))) < 1e-12,
  "N3: `_actn_who` at fp_rate 1 = `_who`; ccm > 0 branch unchanged at any fp_rate"
)

# ---------------------------------------------------------------------------
# N4: `_actn` names registered, dispatched to the norm factories, keys
#     resolved by resolve_filter_cfg() and validated.
cat("\n## N4: enum, dispatch, config\n")
expect_true(
  all(c(ACTN, ACTNW) %in% .SP_ALL_MEAS) &&
    all(c(ACTN, ACTNW) %in% .SP_ACT_MEAS) &&
    all(c(ACTN, ACTNW) %in% .SP_DETECT_MEAS) &&
    all(c(ACTN, ACTNW) %in% .SP_MON_MEAS) &&
    ACTNW %in% .SP_WHO_MEAS &&
    !(ACTN %in% .SP_WHO_MEAS),
  "N4: `_actn` names registered in the act / detect / mon / who sets"
)
expect_true(
  !any(c(ACT, ACTW) %in% .SP_ACTN_MEAS),
  "N4: `_act` names are not `_actn`"
)
d_n <- as.character(.select_dunit(
  ACTN,
  0.1,
  0.001,
  case_fp_rate = 0.02,
  case_fp_mean = 0.3
))
expect_true(
  grepl(.sp_case_actpen_zi_norm_c(0.1, 0.02, 0.3), d_n, fixed = TRUE),
  "N4: `_actn` dispatches the norm case term with the configured eps / mu_0"
)
whon <- do.call(.sp_case_who_norm_c, c(WHO_ARGS, fp_rate = 0.5))
expect_true(
  grepl(
    whon,
    as.character(.select_dunit(ACTNW, 0.1, 0.001, who_c = whon)),
    fixed = TRUE
  ),
  "N4: `_actn_who` carries the WHO term it is given"
)
fc <- resolve_filter_cfg(list(
  measurement = ACTN,
  monthly_case_obs = "zero_filled",
  case_fp_rate = 0.05,
  case_fp_mean = 0.2,
  who_obs = list(fp_rate = 0.3)
))
expect_true(
  fc$case_fp_rate == 0.05 && fc$case_fp_mean == 0.2,
  "N4: resolve_filter_cfg() carries case_fp_rate / case_fp_mean"
)
fc0 <- resolve_filter_cfg(list())
expect_true(
  fc0$case_fp_rate == .SP_CASE_FP_RATE &&
    fc0$case_fp_mean == .SP_MON_MU_FLOOR &&
    fc0$who_obs$fp_rate == .SP_CASE_FP_RATE,
  "N4: defaults are eps = 1, mu_0 = the `_act` floor, who fp_rate = 1"
)
expect_true(
  inherits(
    try(resolve_filter_cfg(list(measurement = ACTN)), silent = TRUE),
    "try-error"
  ),
  "N4: `_actn` without monthly_case_obs is an error (the `_act` guard)"
)
bad <- list(
  c(0, 0.1),
  c(1.5, 0.1),
  c(NA, 0.1),
  c(0.5, 0),
  c(0.5, -1),
  c(0.5, Inf)
)
expect_true(
  all(vapply(
    bad,
    function(b) {
      inherits(
        try(.sp_case_actpen_zi_norm_c(0.1, b[1], b[2]), silent = TRUE),
        "try-error"
      )
    },
    logical(1)
  )),
  "N4: .check_case_fp() rejects eps outside (0,1] and mu_0 <= 0 / non-finite"
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
