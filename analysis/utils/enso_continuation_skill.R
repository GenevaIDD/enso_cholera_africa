# enso_continuation_skill.R — skill by lead of the ENSO continuation past the
#   CFSv2 span, scored on observed seams (perfect-seam hindcast).
#
# Provides:
#   compute_enso_continuation_skill()   scores by lead of every method
#   enso_continuation_hindcast()        continuation ensembles per seam
#   enso_skill_by_lead()                fair CRPS by lead
#   enso_skill_seams(), enso_skill_seam_state()   usable seams + seam state
#   enso_skill_method_neutral(), enso_skill_method_ar1(),
#   enso_skill_method_clim(), enso_skill_method_analog()   method closures
#   enso_continuation_skill_settings()  the scoring settings
#
# Details:
#   Behind `03_prepare_enso_continuation_library.R --validate`, which writes
#   generated_data/enso_continuation_skill.csv for 09. The seam state is the
#   OBSERVED ONI, so the scores are optimistic against a CFSv2-set seam.
#   Callers source base.R, enso.R (continue_enso_analog(),
#   enso_continuation_features()) and obs_sim_compare.R (crps_fair()).

# Design: each observed seam month s is continued to leads 1..H by each
# method; lead k is the k-th month after the seam. The seam state is the
# observed ONI, so scores are optimistic against a CFSv2-set seam. Scored with
# the fair CRPS (obs_sim_compare.R::crps_fair(); 1 member = |x - y|).
# Methods are closures function(ctx) -> [H x n] ensemble matrix, ctx = list(s,
# x0, state, H, obs, obs_month, obs_date), called in the given order once per
# seam (RNG stream fixed by seed + method list): enso_skill_method_neutral(),
# _ar1(phi), _clim() (leave-out climatology) and _analog(lib) (the deployed
# k-NN continuation, enso.R::continue_enso_analog()).

#' enso_skill_seams
#'
#' Indices of the observed months usable as seams: inside [from, to], at least
#' `history` months before (the phase feature reads ONI(s - 6)) and `horizon`
#' months of observed future.
enso_skill_seams <- function(
  obs_date,
  from = "1982-01-01",
  to = "2020-12-01",
  horizon = 24L,
  history = 6L
) {
  obs_date <- as.Date(obs_date)
  i <- seq_along(obs_date)
  which(
    obs_date >= as.Date(from) &
      obs_date <= as.Date(to) &
      i > history &
      i + horizon <= length(obs_date)
  )
}

#' enso_skill_seam_state
#'
#' The analog features of observed month s: list(mon, v, tend, phase).
enso_skill_seam_state <- function(obs, obs_month, s) {
  list(
    mon = as.integer(obs_month[s]),
    v = obs[s],
    tend = obs[s] - obs[s - 3L],
    phase = obs[s - 6L]
  )
}

#' enso_skill_method_neutral
#'
#' ONI = 0 at every lead (the simulator's `neutral` continuation).
enso_skill_method_neutral <- function() {
  function(ctx) matrix(0, ctx$H, 1L)
}

#' enso_skill_method_ar1
#'
#' Deterministic AR(1) decay of the seam value, phi^k * ONI(s) (the simulator's
#' `ar1`, here applied in ONI space).
enso_skill_method_ar1 <- function(phi = 0.969) {
  function(ctx) matrix(phi^seq_len(ctx$H) * ctx$x0, ctx$H, 1L)
}

#' enso_skill_method_clim
#'
#' Climatology reference: at lead k, n_ens draws (with replacement) from the
#' observed ONI of the target calendar month, over months from `pool_from` on,
#' leaving out the `exclude` months either side of the seam (so the verifying
#' observations of leads <= exclude are never in their own reference).
enso_skill_method_clim <- function(
  n_ens = 100L,
  pool_from = "1950-01-01",
  exclude = 36L
) {
  force(n_ens)
  force(pool_from)
  force(exclude)
  function(ctx) {
    s <- ctx$s
    pool_i <- setdiff(
      which(ctx$obs_date >= as.Date(pool_from)),
      seq(s - exclude, s + exclude)
    )
    t(vapply(
      seq_len(ctx$H),
      function(k) {
        mk <- ((ctx$obs_month[s] - 1L + k) %% 12L) + 1L
        sample(
          ctx$obs[pool_i][ctx$obs_month[pool_i] == mk],
          n_ens,
          replace = TRUE
        )
      },
      numeric(n_ens)
    ))
  }
}

#' enso_skill_method_analog
#'
#' The deployed analog continuation: n_ens replicates of the seam state passed
#' to enso.R::continue_enso_analog() (one neighbour set, one analog draw per
#' replicate, faded join). `lib` is a loaded library
#' (load_enso_continuation_library() or enso_continuation_library_rehorizon()).
enso_skill_method_analog <- function(
  lib,
  n_ens = 100L,
  k = 30L,
  blend_months = 3,
  min_sep = 24L,
  max_per_model = Inf,
  exclude_models = character(0)
) {
  force(lib)
  force(n_ens)
  force(k)
  force(blend_months)
  force(min_sep)
  force(max_per_model)
  force(exclude_models)
  function(ctx) {
    st <- ctx$state
    continue_enso_analog(
      lib,
      seam_v = rep(st$v, n_ens),
      seam_tend = rep(st$tend, n_ens),
      seam_phase = rep(st$phase, n_ens),
      seam_month = st$mon,
      n_months = ctx$H,
      k = k,
      blend_months = blend_months,
      min_sep = min_sep,
      max_per_model = max_per_model,
      exclude_models = exclude_models
    )
  }
}

#' enso_continuation_library_rehorizon
#'
#' The same library with the candidate rule `ok` recomputed for another
#' horizon (a start needs `horizon` months of future), the feature sds over the
#' new candidates and the by-month index rebuilt. The skill scoring uses a
#' 24-month candidate rule; the deployed library uses 48 (its --horizon_months).
enso_continuation_library_rehorizon <- function(lib, horizon) {
  horizon <- as.integer(horizon)
  feats <- do.call(
    rbind,
    lapply(names(lib$series), function(m) {
      f <- lib$features[lib$features$model == m, , drop = FALSE]
      f <- f[order(f$i), , drop = FALSE]
      enso_continuation_features(
        lib$series[[m]],
        f$mon,
        model = m,
        horizon = horizon
      )
    })
  )
  rownames(feats) <- NULL
  ok <- feats$ok
  out <- lib
  out$features <- feats
  out$feature_sd <- vapply(
    c("v", "tend", "phase"),
    function(f) stats::sd(feats[[f]][ok]),
    numeric(1)
  )
  cand <- feats[ok, c("model", "i", "mon", "v", "tend", "phase")]
  out$by_month <- split(cand, factor(cand$mon, levels = 1:12))
  out$meta$horizon_months <- horizon
  out
}

#' enso_continuation_hindcast
#'
#' Score every method on every seam. `obs` / `obs_date` / `obs_month` are the
#' observed ONI series in time order; `methods` a named list of method
#' closures (called in list order on each seam). Sets the seed once, when given.
#' Returns one row per seam x method x lead: seam, x0, method, lead, y (the
#' observed ONI), crps (fair CRPS; |mean - y| for a one-member ensemble), mean,
#' sd (0 for a one-member ensemble).
enso_continuation_hindcast <- function(
  obs,
  obs_date,
  obs_month,
  methods,
  seams,
  horizon = 24L,
  seed = NULL,
  progress = 100L
) {
  stopifnot(
    length(obs) == length(obs_date),
    length(obs) == length(obs_month),
    is.list(methods),
    length(names(methods)) == length(methods),
    all(nzchar(names(methods)))
  )
  H <- as.integer(horizon)
  if (any(seams + H > length(obs)) || any(seams <= 6L)) {
    stop(
      "enso_continuation_hindcast: seam(s) without 6 months of history or ",
      H,
      " months of observed future",
      call. = FALSE
    )
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  obs_date <- as.Date(obs_date)
  obs_month <- as.integer(obs_month)
  rows <- vector("list", length(seams))
  for (j in seq_along(seams)) {
    s <- seams[j]
    y <- obs[s + seq_len(H)]
    ctx <- list(
      s = s,
      x0 = obs[s],
      state = enso_skill_seam_state(obs, obs_month, s),
      H = H,
      obs = obs,
      obs_month = obs_month,
      obs_date = obs_date
    )
    rows[[j]] <- do.call(
      rbind,
      lapply(names(methods), function(m) {
        e <- methods[[m]](ctx)
        if (!is.matrix(e) || nrow(e) != H) {
          stop(
            "enso_continuation_hindcast: method '",
            m,
            "' did not return an [H x n] matrix",
            call. = FALSE
          )
        }
        data.frame(
          seam = obs_date[s],
          x0 = obs[s],
          method = m,
          lead = seq_len(H),
          y = y,
          crps = vapply(
            seq_len(H),
            function(k) crps_fair(e[k, ], y[k]),
            numeric(1)
          ),
          mean = rowMeans(e),
          sd = if (ncol(e) > 1L) apply(e, 1L, stats::sd) else rep(0, H),
          stringsAsFactors = FALSE
        )
      })
    )
    if (progress > 0L && j %% progress == 0L) {
      message("  seam ", j, "/", length(seams))
    }
  }
  do.call(rbind, rows)
}

#' enso_skill_by_lead
#'
#' Per method x lead: mean fair CRPS, RMSE of the ensemble mean, spread
#' (root mean ensemble variance), correlation of the ensemble mean with the
#' seam value (NA for a constant forecast), number of seams, and for every
#' reference method present the reference CRPS and the CRPSS 1 - crps / crps_ref
#' (columns crps_<ref>, crpss_vs_<ref>; a reference scored against itself is 0).
enso_skill_by_lead <- function(sc, references = c("clim", "neutral")) {
  key <- paste(sc$method, sc$lead, sep = "\r")
  out <- do.call(
    rbind,
    lapply(split(seq_len(nrow(sc)), factor(key, unique(key))), function(ix) {
      d <- sc[ix, , drop = FALSE]
      r <- if (isTRUE(stats::sd(d$mean) > 0) && isTRUE(stats::sd(d$x0) > 0)) {
        stats::cor(d$x0, d$mean)
      } else {
        NA_real_
      }
      data.frame(
        method = d$method[1L],
        lead = d$lead[1L],
        n_seams = length(unique(d$seam)),
        crps = mean(d$crps),
        rmse = sqrt(mean((d$mean - d$y)^2)),
        spread = sqrt(mean(d$sd^2)),
        r_seam = r,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  for (ref in intersect(references, unique(out$method))) {
    rr <- out[out$method == ref, c("lead", "crps")]
    cr <- rr$crps[match(out$lead, rr$lead)]
    out[[paste0("crps_", ref)]] <- cr
    out[[paste0("crpss_vs_", ref)]] <- 1 - out$crps / cr
  }
  out
}

#' compute_enso_continuation_skill
#'
#' The skill table over several seeds: enso_continuation_hindcast() once per seed
#' (same methods, same seams), enso_skill_by_lead() per seed, then per method x
#' lead the seed means of crps / rmse / spread / r_seam, the CRPSS of the seed-
#' mean CRPS against each reference, and the min / max of the per-seed CRPSS
#' (the Monte-Carlo range). Columns: method, lead, n_seams, n_seeds, crps,
#' rmse, spread, r_seam, then per reference crps_<ref>, crpss_vs_<ref>,
#' crpss_vs_<ref>_min, crpss_vs_<ref>_max.
compute_enso_continuation_skill <- function(
  obs,
  obs_date,
  obs_month,
  methods,
  seams,
  horizon = 24L,
  seeds = 20261018L,
  references = c("clim", "neutral"),
  progress = 0L
) {
  per_seed <- lapply(seeds, function(seed_i) {
    sc <- enso_continuation_hindcast(
      obs,
      obs_date,
      obs_month,
      methods,
      seams,
      horizon = horizon,
      seed = seed_i,
      progress = progress
    )
    cbind(seed = seed_i, enso_skill_by_lead(sc, references = references))
  })
  all <- do.call(rbind, per_seed)
  key <- paste(all$method, all$lead, sep = "\r")
  out <- do.call(
    rbind,
    lapply(split(seq_len(nrow(all)), factor(key, unique(key))), function(ix) {
      d <- all[ix, , drop = FALSE]
      data.frame(
        method = d$method[1L],
        lead = d$lead[1L],
        n_seams = d$n_seams[1L],
        n_seeds = nrow(d),
        crps = mean(d$crps),
        rmse = mean(d$rmse),
        spread = mean(d$spread),
        r_seam = mean(d$r_seam),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  for (ref in intersect(references, unique(out$method))) {
    rr <- out[out$method == ref, c("lead", "crps")]
    cr <- rr$crps[match(out$lead, rr$lead)]
    out[[paste0("crps_", ref)]] <- cr
    out[[paste0("crpss_vs_", ref)]] <- 1 - out$crps / cr
    rng <- tapply(
      all[[paste0("crpss_vs_", ref)]],
      factor(key, unique(key)),
      range
    )
    out[[paste0("crpss_vs_", ref, "_min")]] <- vapply(rng, `[`, numeric(1), 1L)
    out[[paste0("crpss_vs_", ref, "_max")]] <- vapply(rng, `[`, numeric(1), 2L)
  }
  out
}

#' enso_continuation_skill_settings
#'
#' The analog knobs to score, as the simulator resolves them: the defaults'
#' simulator.enso.continuation block, overlaid leaf by leaf by a deploy YAML's
#' deploy.enso_continuation when one is given. `max_per_model` .inf / NULL ->
#' Inf; `exclude_models` -> character. Returns list(method, library, k,
#' blend_months, min_sep, max_per_model, exclude_models, source).
enso_continuation_skill_settings <- function(
  defaults_yaml = "analysis/config_defaults.yml",
  deploy_yaml = NULL
) {
  d <- yaml::read_yaml(defaults_yaml)[["simulator"]][["enso"]][[
    "continuation"
  ]] %||%
    list()
  src <- defaults_yaml
  if (!is.null(deploy_yaml) && nzchar(deploy_yaml)) {
    dep <- yaml::read_yaml(deploy_yaml)[["deploy"]][["enso_continuation"]]
    if (is.null(dep)) {
      stop(
        "enso_continuation_skill_settings: no deploy.enso_continuation in ",
        deploy_yaml,
        call. = FALSE
      )
    }
    for (nm in names(dep)) {
      d[[nm]] <- dep[[nm]]
    }
    src <- deploy_yaml
  }
  mpm <- suppressWarnings(as.numeric(d[["max_per_model"]] %||% Inf))
  if (!length(mpm) || is.na(mpm)) {
    mpm <- Inf
  }
  list(
    method = d[["method"]] %||% "neutral",
    library = d[["library"]] %||% "generated_data/enso_continuation_library.qs",
    k = as.integer(d[["k"]] %||% 30L),
    blend_months = as.numeric(d[["blend_months"]] %||% 3),
    min_sep = as.integer(d[["min_sep"]] %||% 24L),
    max_per_model = mpm,
    exclude_models = as.character(unlist(d[["exclude_models"]])),
    source = src
  )
}
