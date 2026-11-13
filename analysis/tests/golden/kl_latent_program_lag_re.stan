// patched by utils/dlnm_kernel.R::patch_kl_latent_class_stancode (brms_kl latent classes; J = 3, 7 lags; lag_re)
// generated with brms 2.23.0
functions {

  // ---- brms_kl latent classes (utils/dlnm_kernel.R) ----
  // Row log-likelihood of every row under every class: eta = lin + U[tix, k],
  // U = KBu * Theta the kernel term of every distinct exposure row (theta_k =
  // the class per-lag exposure coefficients, lag-major, whatever the form --
  // Gaussian kernel or class lag basis -- that produced them), tix = the
  // distinct exposure row of each row.
  matrix kl_class_rowll(vector nlp_lin, matrix U, array[] int tix,
                        vector Yv, vector lgY1) {
    int N = rows(nlp_lin);
    int K = cols(U);
    matrix[N, K] ll;
    for (k in 1:K) {
      vector[N] eta = nlp_lin + col(U, k)[tix];
      ll[:, k] = Yv .* eta - exp(eta) - lgY1;
    }
    return ll;
  }
  // Non-centred DAGAR field (Datta et al. 2019), one field per column of z:
  // a deterministic map from standard normals, so no Jacobian / log-det term.
  // Ported verbatim from archive/cholera_regime_sim (stan/regime_mixture_nc.stan,
  // dagar_field()). A node with no parent (first node of a component, island)
  // gets sd tau.
  matrix kl_dagar_field(matrix z, real rho, real tau,
                        array[] int n_par, array[] int par_nodes,
                        array[] int par_start) {
    int Ii = rows(z);
    int Uc = cols(z);
    matrix[Ii, Uc] out;
    real rho2 = square(rho);
    for (i in 1:Ii) {
      real denom = 1 + (n_par[i] - 1) * rho2;  // == 1 - rho2 when n_par[i]==0
      real b     = rho / denom;
      real sd_i  = tau * sqrt((1 - rho2) / denom);
      row_vector[Uc] parent_sum = rep_row_vector(0.0, Uc);
      for (p in 1:n_par[i])
        parent_sum += out[par_nodes[par_start[i] + p - 1], ];
      out[i, ] = b * parent_sum + sd_i * z[i, ];
    }
    return out;
  }
  // Membership log-probabilities of one country: log_softmax of its logits, or
  // a point mass on class `fix` when the country is a seed (kernel.pin).
  row_vector kl_member_logp(row_vector logit, int fix) {
    row_vector[cols(logit)] out;
    if (fix == 0) return log_softmax(logit')';
    out = rep_row_vector(negative_infinity(), cols(logit));
    out[fix] = 0;
    return out;
  }
  // Sum row log-likelihoods within each country (rows sorted by country).
  matrix kl_country_ll(matrix ll, array[] int cs, array[] int ce) {
    int C = size(cs);
    int K = cols(ll);
    matrix[C, K] out;
    for (c in 1:C) for (k in 1:K) out[c, k] = sum(ll[cs[c]:ce[c], k]);
    return out;
  }
  // Form discrete: country log-likelihood of every class under every kernel
  // q (column q of KBq[j] = exposure basis j weighted over the lags by kernel
  // q). out[k][c, q] = sum_{i in c} log Poisson(y_i | lin_i + KBq[, q] kbeta[, k]
  // + E[i, q]); lgY1c = country sums of lgamma(y + 1); E (used when use_E) =
  // the country/ADM1 deviations of the exposure coefficients times the kernel
  // basis (kl_disc_re_term()), the same for every class.
  array[] matrix kl_disc_country_ll(vector nlp_lin, array[] matrix KBq, matrix kbeta,
                                    vector Yv, vector lgY1c, array[] int cs, array[] int ce,
                                    matrix E, int use_E) {
    int N = rows(nlp_lin);
    int Q = cols(KBq[1]);
    int J = size(KBq);
    int K = cols(kbeta);
    int C = size(cs);
    array[K] matrix[C, Q] out;
    matrix[N, Q] eta0 = rep_matrix(nlp_lin, Q);
    if (use_E) eta0 += E;   // random-slope deviations (class-independent)
    for (k in 1:K) {
      matrix[N, Q] eta = eta0;
      matrix[N, Q] ll;
      for (j in 1:J) eta += kbeta[j, k] * KBq[j];
      ll = diag_pre_multiply(Yv, eta) - exp(eta);
      for (c in 1:C)
        out[k][c] = rep_row_vector(1.0, ce[c] - cs[c] + 1) * ll[cs[c]:ce[c], ] - lgY1c[c];
    }
    return out;
  }
  // Forms gaussian / dlnm, country path (kernel threading or random slopes):
  // class log-probabilities of one country, lsm[k] + sum_i log Poisson(y_i |
  // lin_i + B_i theta_k), or with random slopes lin_i + (B_i T_k)(kbeta_k + D_i),
  // where B_i T_k = the exposure basis of row i weighted over the lags by the
  // class k kernel (T_k maps lag-major basis columns to J exposure columns).
  // The kernel term depends on the row only through its exposure basis B_i,
  // which takes few distinct values (ENSO is one series: one per calendar
  // month), so it is computed once per distinct row t by the caller and
  // indexed here by tix: U[t, k] = B_t theta_k, BW[k][t] = B_t T_k.
  vector kl_ctry_logp(vector lin, array[] int tix, vector y, real lgy, row_vector lsm,
                    matrix U, array[] matrix BW, matrix kbeta, int use_re, matrix D) {
    int K = cols(lsm);
    int n = rows(lin);
    vector[K] lp;
    for (k in 1:K) {
      vector[n] eta;
      if (use_re) {
        matrix[n, rows(kbeta)] Bw = BW[k][tix, ];
        eta = lin + Bw * col(kbeta, k) + rows_dot_product(Bw, D);
      } else {
        eta = lin + col(U, k)[tix];
      }
      lp[k] = lsm[k] + dot_product(y, eta) - sum(exp(eta)) - lgy;
    }
    return lp;
  }
  // reduce_sum partial sum over countries cidx[start..end]: sum of the country
  // mixture log-likelihoods. lin_pad = the linear predictor of each of those
  // countries, padded to the largest one (the sliced argument, so it is not
  // copied per task).
  real kl_partial_ctry(array[] vector lin_pad, int start, int end, array[] int cidx,
                       array[] int tix, vector Yv,
                       vector lgY1c, array[] int cs, array[] int ce, matrix lsm,
                       matrix U, array[] matrix BW, matrix kbeta, int use_re,
                       matrix kre_c, matrix kre_a, array[] int adm1) {
    real tot = 0;
    for (i in start:end) {
      int c = cidx[i];
      int n = ce[c] - cs[c] + 1;
      matrix[use_re ? n : 0, use_re ? rows(kbeta) : 0] D;
      if (use_re)
        D = (rep_matrix(col(kre_c, c), n) + kre_a[:, adm1[cs[c]:ce[c]]])';
      tot += log_sum_exp(kl_ctry_logp(lin_pad[i - start + 1][1:n], tix[cs[c]:ce[c]],
                                    Yv[cs[c]:ce[c]], lgY1c[c], lsm[c], U, BW, kbeta,
                                    use_re, D));
    }
    return tot;
  }
  // reduce_sum partial sum over rows r0 + start - 1 .. r0 + end - 1 of one large
  // country c, for class k: sum_i [y_i eta_ik - exp(eta_ik)] (lgamma term and
  // membership added by the caller). Large countries are split by rows so one
  // country does not bound the threading (kernel.split_big).
  real kl_partial_rows(array[] real lin_s, int start, int end, int r0, int k, int c,
                       array[] int tix, vector Yv, matrix U, array[] matrix BW, matrix kbeta,
                       int use_re, matrix kre_c, matrix kre_a, array[] int adm1) {
    int a = r0 + start - 1;
    int b = r0 + end - 1;
    int n = end - start + 1;
    vector[n] eta;
    if (use_re) {
      matrix[n, rows(kbeta)] Bw = BW[k][tix[a:b], ];
      matrix[n, rows(kbeta)] D = (rep_matrix(col(kre_c, c), n) + kre_a[:, adm1[a:b]])';
      eta = to_vector(lin_s) + Bw * col(kbeta, k) + rows_dot_product(Bw, D);
    } else {
      eta = to_vector(lin_s) + col(U, k)[tix[a:b]];
    }
    return dot_product(Yv[a:b], eta) - sum(exp(eta));
  }
  // Pair path (kernel.pair_collapse): the kernel term of a row depends on the
  // row only through its pair p = (distinct exposure row t, ADM1 unit a) --
  // (t, country) without random slopes -- so with g_pk the class k kernel term
  // of pair p, sum_{i in p} [y_i (lin_i + g_pk) - exp(lin_i + g_pk)]
  //   = sum_{i in p} y_i lin_i + g_pk Yp_p - exp(g_pk) S_p,
  // with Yp_p = sum y_i and S_p = sum exp(lin_i) over the pair rows. This
  // returns the class-dependent part over the pairs of country c:
  // sum_p [g_pk Yp_p - exp(g_pk) S_p] for every k (the row part
  // sum_i y_i lin_i - lgamma(y_i + 1) is the same for every class and is added
  // by the caller). g_pk = U[t_p, k], or with random slopes
  // BW[k][t_p] (kbeta_k + kre_c[, c] + kre_a[, a_p]).
  vector kl_pair_logp(vector S, vector Yp, array[] int pt, array[] int pa, int c,
                      matrix U, array[] matrix BW, matrix kbeta, int use_re,
                      matrix kre_c, matrix kre_a) {
    int K = cols(kbeta);
    int np = rows(S);
    vector[K] out;
    matrix[use_re ? np : 0, use_re ? rows(kbeta) : 0] D;
    if (use_re) D = (rep_matrix(col(kre_c, c), np) + kre_a[:, pa])';
    for (k in 1:K) {
      vector[np] g;
      if (use_re) g = rows_dot_product(BW[k][pt, ], D + rep_matrix(col(kbeta, k)', np));
      else g = col(U, k)[pt];
      out[k] = dot_product(Yp, g) - dot_product(S, exp(g));
    }
    return out;
  }
  // reduce_sum partial sum over countries cidx[start..end] on the pair path:
  // sum of log_sum_exp_k(lsm[c, k] + kl_pair_logp(...)[k]). S_pad = the pair
  // sums S_p of each of those countries, padded to the largest (sliced).
  real kl_partial_pairs(array[] vector S_pad, int start, int end, array[] int cidx,
                        vector Ypair, array[] int pcs, array[] int pce,
                        array[] int pt, array[] int pa, matrix lsm,
                        matrix U, array[] matrix BW, matrix kbeta, int use_re,
                        matrix kre_c, matrix kre_a) {
    real tot = 0;
    for (i in start:end) {
      int c = cidx[i];
      int np = pce[c] - pcs[c] + 1;
      tot += log_sum_exp(lsm[c]' + kl_pair_logp(S_pad[i - start + 1][1:np], Ypair[pcs[c]:pce[c]],
                                                pt[pcs[c]:pce[c]], pa[pcs[c]:pce[c]], c,
                                                U, BW, kbeta, use_re, kre_c, kre_a));
    }
    return tot;
  }
  // Form discrete, kernel.re_slope: E[i, q] = sum_j KBq[j][i, q] D[i, j], with D
  // the country + ADM1 deviation of the row for exposure coefficient j.
  matrix kl_disc_re_term(array[] matrix KBq, matrix D) {
    int N = rows(D);
    int Q = cols(KBq[1]);
    matrix[N, Q] E = rep_matrix(0, N, Q);
    for (j in 1:size(KBq)) E += diag_pre_multiply(col(D, j), KBq[j]);
    return E;
  }
  // kernel.re_slope: ADM1 deviations of the ENSO exposure coefficients (J x NA),
  // zero when that tier is off (kernel.re_slope.adm1: false; sda, za empty).
  matrix kl_adm1_dev(vector sda, matrix za, int use, int J, int NA) {
    if (use) return diag_pre_multiply(sda, za);
    return rep_matrix(0, J, NA);
  }
  // Form discrete: log-likelihood of every (width s, centre combination m),
  // classes marginalised per country: element (s - 1) * M + m is
  //   sum_c log_sum_exp_k( lsm[c, k] + LLq[k][c, qidx[s, k, m]] ).
  // Vectorised over combinations: the log_sum_exp over k is shifted by its
  // elementwise max (the shift cancels in the gradient).
  vector kl_disc_combo_ll(array[] matrix LLq, matrix lsm, array[,,] int qidx) {
    int S = dims(qidx)[1];
    int K = dims(qidx)[2];
    int M = dims(qidx)[3];
    int C = rows(lsm);
    vector[S * M] tot;
    for (s in 1:S) {
      array[K] matrix[C, M] G;
      matrix[C, M] mx;
      matrix[C, M] acc = rep_matrix(0, C, M);
      for (k in 1:K) G[k] = LLq[k][:, qidx[s, k]] + rep_matrix(col(lsm, k), M);
      mx = G[1];
      for (k in 2:K) mx = fmax(mx, G[k]);
      for (k in 1:K) acc += exp(G[k] - mx);
      tot[((s - 1) * M + 1):(s * M)] = (rep_row_vector(1.0, C) * (mx + log(acc)))';
    }
    return tot;
  }

  // ---- kernel.lag_re: country lag offsets (utils/dlnm_kernel.R) ----
  // Lag-weighted exposure basis of every class on the country-stacked distinct
  // exposure rows: BW[k][(c - 1) * T + t] = sum_l w_ckl KBl[l][t], with
  // w_ck = softmax(-0.5 ((lags - kmu_k - d_c) / ksd_k)^2) (KBl[l] = the J
  // exposure-basis columns of lag l).
  array[] matrix kl_lag_re_bw(array[] matrix KBl, vector lags, vector kmu, vector kksd, vector d) {
    int NL = size(KBl);
    int T = rows(KBl[1]);
    int J = cols(KBl[1]);
    int K = rows(kmu);
    int C = rows(d);
    array[K] matrix[C * T, J] BW;
    for (k in 1:K) {
      for (c in 1:C) {
        vector[NL] w = softmax(-0.5 * square((lags - kmu[k] - d[c]) / kksd[k]));
        matrix[T, J] acc = w[1] * KBl[1];
        for (l in 2:NL) acc += w[l] * KBl[l];
        BW[k][((c - 1) * T + 1):(c * T)] = acc;
      }
    }
    return BW;
  }
  // Kernel term without random slopes: U[(c - 1) * T + t, k] = BW[k][.] kbeta[, k].
  matrix kl_lag_re_u(array[] matrix BW, matrix kbeta) {
    int K = cols(kbeta);
    matrix[rows(BW[1]), K] U;
    for (k in 1:K) U[:, k] = BW[k] * col(kbeta, k);
    return U;
  }
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int<lower=1> K_lin;  // number of population-level effects
  matrix[N, K_lin] X_lin;  // population-level design matrix
  // data for spline 1
  int nb_lin_1;  // number of bases
  array[nb_lin_1] int knots_lin_1;  // number of knots
  // basis function matrices
  matrix[N, knots_lin_1[1]] Zs_lin_1_1;
  vector[N] offsets_lin;
  int<lower=1> K_kmu;  // number of population-level effects
  matrix[N, K_kmu] X_kmu;  // population-level design matrix
  int<lower=1> K_ksd;  // number of population-level effects
  matrix[N, K_ksd] X_ksd;  // population-level design matrix
  int<lower=1> K_kb1;  // number of population-level effects
  matrix[N, K_kb1] X_kb1;  // population-level design matrix
  int<lower=1> K_kb2;  // number of population-level effects
  matrix[N, K_kb2] X_kb2;  // population-level design matrix
  int<lower=1> K_kb3;  // number of population-level effects
  matrix[N, K_kb3] X_kb3;  // population-level design matrix
  // covariates for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  vector[N] C_6;
  vector[N] C_7;
  vector[N] C_8;
  vector[N] C_9;
  vector[N] C_10;
  vector[N] C_11;
  vector[N] C_12;
  vector[N] C_13;
  vector[N] C_14;
  vector[N] C_15;
  vector[N] C_16;
  vector[N] C_17;
  vector[N] C_18;
  vector[N] C_19;
  vector[N] C_20;
  vector[N] C_21;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_lin_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_lin_1;
  int prior_only;  // should the likelihood be ignored?

  // ---- brms_kl latent classes (utils/dlnm_kernel.R) ----
  int<lower=1> kl_K;                                // classes (1 = single kernel)
  int<lower=1> kl_C;                                // countries
  array[N] int<lower=1, upper=kl_C> kl_ctry;        // country of each row (rows sorted)
  array[kl_C] int<lower=1> kl_cs;                   // first row of each country
  array[kl_C] int<lower=1> kl_ce;                   // last row of each country
  int<lower=1> kl_J;                                // exposure basis functions
  int<lower=2> kl_NL;                               // lags
  int<lower=1> kl_T;                                // distinct rows of the kernel exposure basis (ENSO: calendar months)
  array[N] int<lower=1, upper=kl_T> kl_tix;         // distinct exposure-basis row of each row
  int<lower=0, upper=1> kl_pairs;                   // country path summed per pair (kernel.pair_collapse)
  int<lower=1> kl_P;                                // pairs: (exposure row, ADM1) with random slopes, else (exposure row, country)
  array[N] int<lower=1, upper=kl_P> kl_pix;         // pair of each row
  array[kl_P] int<lower=1, upper=kl_T> kl_pt;       // exposure row of each pair
  array[kl_P] int<lower=1> kl_pa;                   // ADM1 unit of each pair (pairs sorted by country)
  array[kl_C] int<lower=1> kl_pcs;                  // first pair of each country
  array[kl_C] int<lower=1> kl_pce;                  // last pair of each country
  int<lower=1> kl_maxp;                             // pairs of the country with the most
  vector[kl_NL] kl_lags;                            // lag grid (ascending)
  int<lower=0, upper=1> kl_dlnm;                    // form: 0 = Gaussian kernel, 1 = class lag basis
  int<lower=0, upper=1> kl_disc;                    // form discrete: grid centres, ordered (ties allowed), summed out
  int<lower=1> kl_S;                                // discrete: kernel widths (shared by the classes, summed out)
  int<lower=1> kl_NCt;                              // discrete: candidate centres
  vector[kl_NCt] kl_cent;                           // discrete: centre values (lags)
  vector[kl_S] kl_width;                            // discrete: width values (lags; 0 = single lag)
  matrix[kl_NL, kl_S * kl_NCt] kl_Wq;               // discrete: lag weights of kernel q = (s - 1) * NCt + centre
  int<lower=1> kl_M;                                // discrete: ordered centre combinations
  array[kl_M, kl_K] int<lower=1, upper=kl_NCt> kl_combo;   // discrete: centre index per class, non-decreasing
  vector[kl_M] kl_combo_lp;                         // discrete: log prior of each combination
  int<lower=0, upper=1> kl_re;                      // discrete: country + ADM1 deviations of the exposure coefficients
  int<lower=1> kl_NA;                               // ADM1 units (country:ADM1)
  array[N] int<lower=1, upper=kl_NA> kl_adm1;       // ADM1 unit of each row
  vector<lower=0>[2] kl_re_sd;                      // half-normal sds of the country / ADM1 deviation sds
  int<lower=0, upper=1> kl_re_ctry;                 // country-level deviations on (0 = ADM1 only)
  int<lower=0, upper=1> kl_re_adm1;                 // ADM1-level deviations on (0 = country only)
  int<lower=0, upper=1> kl_re_sdc_fix;              // country deviation sd fixed (kl_re_sdc_value), not estimated
  real<lower=0> kl_re_sdc_value;
  int<lower=0, upper=1> kl_rs;                      // gaussian/dlnm: country path kl_partial_ctry (reduce_sum)
  int<lower=1> kl_grain;                            // reduce_sum grainsize (countries)
  int<lower=1> kl_maxn;                             // rows of the largest country summed whole
  int<lower=0> kl_nsmall;                           // country path: countries summed whole (country reduce_sum)
  array[kl_nsmall] int<lower=1, upper=kl_C> kl_small;
  int<lower=0> kl_nbig;                             // country path: countries split by rows (one reduce_sum per class)
  array[kl_nbig] int<lower=1, upper=kl_C> kl_big;
  int<lower=1> kl_grain_rows;                       // reduce_sum grainsize (rows) for the split countries
  int<lower=0, upper=1> kl_ksd_fix;                 // gaussian: kernel width fixed (kl_ksd_value)
  real<lower=0> kl_ksd_value;
  int<lower=1> kl_BL;                               // lag-basis functions (form dlnm)
  matrix[kl_NL, kl_BL] kl_Lb;                       // B-spline lag basis, rows = lags (form dlnm)
  int<lower=0, upper=1> kl_hs;                      // dlnm prior: 1 = regularised group horseshoe, 0 = normal
  real<lower=0> kl_hs_scale_global;                 // half-Cauchy scale of the per-class global shrinkage
  real<lower=0> kl_hs_slab_scale;                   // slab scale (caps active coefficients)
  real<lower=0> kl_hs_slab_df;
  real<lower=0> kl_alpha;                           // Dirichlet concentration
  int<lower=0, upper=1> kl_mu_normal;               // normal (1) or uniform (0) kmu prior
  vector[kl_K] kl_mu_mean;
  real<lower=0> kl_mu_sd;
  real<lower=0> kl_sd_lb;                           // lower bound of ksd
  real<lower=0> kl_sd_sd;                           // half-normal sd of ksd
  real<lower=0> kl_beta_sd;
  int<lower=0, upper=1> kl_use_pg;                  // prior-partition anchoring (kappa * A) on?
  int<lower=0, upper=1> kl_mu_ordered;              // identify by ordered kmu (lag centre)
  int<lower=0, upper=1> kl_erf_ordered;             // identify by ordered ERF contrast
  vector[kl_J] kl_erf_v;                            // B(erf_at) - B(cen): ERF contrast direction
  matrix[kl_J, kl_J - 1] kl_erf_Q;                  // orthonormal complement of kl_erf_v
  matrix[kl_C, kl_K] kl_A;                          // one-hot prior membership
  int<lower=0, upper=1> kl_kappa_fixed;
  real<lower=0> kl_kappa_value;
  real<lower=0> kl_kappa_sd;
  int<lower=0, upper=1> kl_use_sp;                  // spatial membership term on?
  int<lower=0> kl_field_cols;                       // field columns: K - 1 (reference class) or K (symmetric)
  int<lower=0> kl_n_edges;
  array[kl_n_edges] int<lower=1, upper=kl_C> kl_node1;
  array[kl_n_edges] int<lower=1, upper=kl_C> kl_node2;
  int<lower=0> kl_n_comp;                           // connected components
  array[kl_C] int<lower=0, upper=kl_n_comp> kl_comp;
  vector[kl_n_comp] kl_comp_n;
  real<lower=0> kl_tau_sd;
  int<lower=0, upper=1> kl_sp_dagar;                // spatial prior: 0 = ICAR, 1 = DAGAR
  int<lower=0> kl_NP;                               // DAGAR parent links
  array[kl_C] int<lower=0> kl_n_par;
  array[kl_NP] int<lower=1, upper=kl_C> kl_par_nodes;
  array[kl_C] int<lower=1> kl_par_start;
  real<lower=0> kl_rho_a;                           // beta prior on the DAGAR rho
  real<lower=0> kl_rho_b;
  int<lower=0, upper=1> kl_gq_row_loglik;
  array[kl_C] int<lower=0, upper=kl_K> kl_fix;      // seed anchoring (kernel.pin): 0 = free, k = country fixed to class k

  // ---- kernel.lag_re (utils/dlnm_kernel.R::kl_lag_re_cfg()) ----
  real<lower=0> kl_lag_re_sd;                       // half-normal sd of the country lag-offset sd (lags)
}
transformed data {

  // ---- brms_kl latent classes ----
  real kl_lo = kl_lags[1];
  real kl_hi = kl_lags[kl_NL];
  vector[N] kl_Yv = to_vector(Y);
  vector[N] kl_lgY1;
  matrix[N, kl_NL * kl_J] kl_KB;   // exposure basis at every lag, lag-major (brms C_1..C_<NL*J>)
  matrix[kl_T, kl_NL * kl_J] kl_KBu;   // its distinct rows (kl_KB[n] = kl_KBu[kl_tix[n]])
  vector[kl_P] kl_Ypair = rep_vector(0, kl_P);   // pair sums of Y
  array[kl_C] int kl_call;         // 1..C (pair path: every country summed whole)
  vector[kl_C] kl_lgY1c;           // country sums of lgamma(Y + 1)
  array[kl_disc ? kl_J : 0] matrix[N, kl_S * kl_NCt] kl_KBq;   // discrete: kernel-weighted exposure basis
  array[kl_disc ? kl_S : 0, kl_K, kl_M] int kl_qidx;          // discrete: kernel of class k in (s, m)
  vector[kl_disc ? kl_S * kl_M : 0] kl_qlp;                   // discrete: log prior of (s, m)
  // ---- brms_kl latent classes ----
  for (n in 1:N) kl_lgY1[n] = lgamma(Y[n] + 1);
  kl_KB[:, 1] = C_1;
  kl_KB[:, 2] = C_2;
  kl_KB[:, 3] = C_3;
  kl_KB[:, 4] = C_4;
  kl_KB[:, 5] = C_5;
  kl_KB[:, 6] = C_6;
  kl_KB[:, 7] = C_7;
  kl_KB[:, 8] = C_8;
  kl_KB[:, 9] = C_9;
  kl_KB[:, 10] = C_10;
  kl_KB[:, 11] = C_11;
  kl_KB[:, 12] = C_12;
  kl_KB[:, 13] = C_13;
  kl_KB[:, 14] = C_14;
  kl_KB[:, 15] = C_15;
  kl_KB[:, 16] = C_16;
  kl_KB[:, 17] = C_17;
  kl_KB[:, 18] = C_18;
  kl_KB[:, 19] = C_19;
  kl_KB[:, 20] = C_20;
  kl_KB[:, 21] = C_21;
  for (c in 1:kl_C) kl_lgY1c[c] = sum(kl_lgY1[kl_cs[c]:kl_ce[c]]);
  for (n in 1:N) kl_KBu[kl_tix[n]] = kl_KB[n];
  for (n in 1:N) kl_Ypair[kl_pix[n]] += kl_Yv[n];
  for (c in 1:kl_C) kl_call[c] = c;
  if (kl_disc) {
    for (j in 1:kl_J) {
      matrix[N, kl_NL] bj;
      for (l in 1:kl_NL) bj[:, l] = col(kl_KB, (l - 1) * kl_J + j);
      kl_KBq[j] = bj * kl_Wq;
    }
    for (s in 1:kl_S) for (k in 1:kl_K) for (m in 1:kl_M)
      kl_qidx[s, k, m] = (s - 1) * kl_NCt + kl_combo[m, k];
    for (s in 1:kl_S) kl_qlp[((s - 1) * kl_M + 1):(s * kl_M)] = kl_combo_lp - log(kl_S);
  }

  // ---- kernel.lag_re ----
  array[kl_NL] matrix[kl_T, kl_J] kl_KBl;           // exposure basis of each lag (J columns of kl_KBu)
  array[N] int kl_tixc;                             // row -> country-stacked distinct exposure row
  array[kl_P] int kl_ptc;                           // pair -> country-stacked distinct exposure row

  // ---- kernel.lag_re ----
  for (l in 1:kl_NL) kl_KBl[l] = kl_KBu[:, ((l - 1) * kl_J + 1):(l * kl_J)];
  for (n in 1:N) kl_tixc[n] = kl_tix[n] + (kl_ctry[n] - 1) * kl_T;
  for (c in 1:kl_C) for (q in kl_pcs[c]:kl_pce[c]) kl_ptc[q] = kl_pt[q] + (c - 1) * kl_T;
}
parameters {

  // ---- brms_kl latent classes ----
  // Gaussian-kernel form (kl_dlnm = 0; form discrete keeps only kbeta_free)
  ordered[(kl_mu_ordered && !kl_dlnm) ? kl_K : 0] kmu_ord;   // identify = kmu: lag order fixes labels
  vector<lower=kl_lo, upper=kl_hi>[(kl_mu_ordered || kl_dlnm || kl_disc) ? 0 : kl_K] kmu_free;
  vector<lower=kl_sd_lb>[(kl_dlnm || kl_disc || kl_ksd_fix) ? 0 : kl_K] ksd;
  matrix[kl_J, (kl_erf_ordered || kl_dlnm) ? 0 : kl_K] kbeta_free;
  ordered[kl_erf_ordered ? kl_K : 0] kc;            // identify = erf: log RR at erf_at vs cen, increasing
  matrix[kl_J - 1, kl_erf_ordered ? kl_K : 0] keta; // free basis directions orthogonal to the contrast
  // Class lag-basis form (kl_dlnm = 1): non-centred coefficients, column-major
  // J x BL per class (exposure basis j within lag-basis function b)
  matrix[kl_dlnm ? kl_J * kl_BL : 0, kl_K] kz;
  vector<lower=0>[(kl_dlnm && kl_hs) ? kl_K : 0] khs_tau;              // per-class global shrinkage
  matrix<lower=0>[(kl_dlnm && kl_hs) ? kl_BL : 0, kl_K] khs_lambda;    // per lag window (group)
  vector<lower=0>[(kl_dlnm && kl_hs) ? kl_K : 0] khs_c2;               // slab
  // Random slopes on the exposure coefficients (form discrete, kernel.re_slope)
  matrix[(kl_re && kl_re_ctry) ? kl_J : 0, kl_C] kre_zc;   // country deviations, standardised
  matrix[(kl_re && kl_re_adm1) ? kl_J : 0, kl_NA] kre_za;   // ADM1 deviations, standardised
  vector<lower=0>[(kl_re && kl_re_ctry && !kl_re_sdc_fix) ? kl_J : 0] kre_sdc;   // country sd per exposure coefficient
  vector<lower=0>[(kl_re && kl_re_adm1) ? kl_J : 0] kre_sda;   // ADM1 sd per exposure coefficient
  // Membership
  simplex[kl_K] kpi;
  vector<lower=0>[(kl_use_pg && !kl_kappa_fixed) ? 1 : 0] kkappa_free;
  matrix[kl_use_sp ? kl_C : 0, kl_field_cols] kphi; // ICAR field / DAGAR innovations
  vector<lower=0>[kl_use_sp ? 1 : 0] ktau;          // spatial scale
  vector<lower=0, upper=1>[(kl_use_sp && kl_sp_dagar) ? 1 : 0] krho;   // DAGAR dependence

  vector[kl_C] klag_z;                              // kernel.lag_re: country lag offsets, standardised
  vector<lower=0>[1] klag_sd;                       // kernel.lag_re: their sd (lags)
  vector[K_lin] b_lin;  // regression coefficients
  // parameters for spline 1
  // standardized penalized spline coefficients
  vector[knots_lin_1[1]] zs_lin_1_1;
  vector<lower=0>[nb_lin_1] sds_lin_1;  // SDs of penalized spline coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
}
transformed parameters {

  // ---- brms_kl latent classes ----
  vector[(kl_dlnm || kl_disc) ? 0 : kl_K] kmu;      // Gaussian form: class lag centres
  matrix[(kl_dlnm || kl_disc) ? 0 : kl_NL, kl_K] kw;   // Gaussian form: class lag weights
  vector[(kl_dlnm || kl_disc) ? 0 : kl_K] kksd;     // Gaussian form: kernel widths (fixed or ksd)
  array[(kl_rs && kl_re && !kl_dlnm && !kl_disc) ? kl_K : 0] matrix[kl_NL * kl_J, kl_J] kTw;   // lag-weighting maps
  matrix[kl_J, kl_K] kbeta;                         // sustained-exposure coefficients (sum over lags)
  matrix[kl_disc ? 0 : kl_NL * kl_J, kl_K] ktheta;  // per-lag exposure coefficients, lag-major
  vector[kl_K] kcon;                                // log RR of sustained ENSO erf_at vs cen
  real kkappa;                                      // prior-partition strength
  matrix[kl_C, kl_K] klogit;                        // membership logits
  matrix[kl_re ? kl_J : 0, kl_C] kre_c;             // country deviations of the exposure coefficients
  matrix[kl_use_sp ? kl_C : 0, kl_field_cols] kfield;   // spatial membership field (scaled)

  vector[kl_C] kclag;                               // kernel.lag_re: country lag offsets (lags)
  // penalized spline coefficients
  vector[knots_lin_1[1]] s_lin_1_1;
  vector[N_1] r_1_lin_1;  // actual group-level effects
  vector[N_2] r_2_lin_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  // compute penalized spline coefficients
  s_lin_1_1 = sds_lin_1[1] * zs_lin_1_1;
  r_1_lin_1 = (sd_1[1] * (z_1[1]));
  r_2_lin_1 = (sd_2[1] * (z_2[1]));
  lprior += normal_lpdf(b_lin[1] | 0, 1.5);
  lprior += student_t_lpdf(sds_lin_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += normal_lpdf(sd_1 | 0, 2.5)
    - 1 * normal_lccdf(0 | 0, 2.5);
  lprior += normal_lpdf(sd_2 | 0, 2.5)
    - 1 * normal_lccdf(0 | 0, 2.5);

  // ---- brms_kl latent classes ----
  if (kl_disc) {
    // Discrete centres: the lag weights are data (kl_Wq); only the exposure
    // coefficients are parameters (centres and width are summed out).
    kbeta = kbeta_free;
  } else if (!kl_dlnm) {
    // Gaussian kernel: theta at lag l = w_l(kmu, ksd) * kbeta (weights sum to 1)
    if (kl_mu_ordered) {
      kmu = kl_lo + (kl_hi - kl_lo) * inv_logit(kmu_ord);
      lprior += sum(log(kl_hi - kl_lo) + log_inv_logit(kmu_ord) + log1m_inv_logit(kmu_ord));  // Jacobian
    } else {
      kmu = kmu_free;
    }
    // kbeta: free, or the linear (constant-Jacobian) map from the ordered ERF
    // contrast kc and the orthogonal directions keta, so dot(v, kbeta[, k]) = kc[k].
    if (kl_erf_ordered) {
      for (k in 1:kl_K)
        kbeta[:, k] = kc[k] * kl_erf_v / dot_self(kl_erf_v) + kl_erf_Q * keta[:, k];
    } else {
      kbeta = kbeta_free;
    }
    if (kl_ksd_fix) kksd = rep_vector(kl_ksd_value, kl_K);
    else kksd = ksd;
    for (k in 1:kl_K) {
      kw[:, k] = softmax(-0.5 * square((kl_lags - kmu[k]) / kksd[k]));
      for (l in 1:kl_NL) ktheta[((l - 1) * kl_J + 1):(l * kl_J), k] = kw[l, k] * kbeta[:, k];
    }
    if (kl_rs && kl_re)
      for (k in 1:kl_K) {
        kTw[k] = rep_matrix(0, kl_NL * kl_J, kl_J);
        for (l in 1:kl_NL) for (j in 1:kl_J) kTw[k][(l - 1) * kl_J + j, j] = kw[l, k];
      }
  } else {
    // Class lag basis: J x NL coefficient surface = (J x BL coefficients) * Lb'.
    // Regularised group horseshoe (Piironen & Vehtari 2017): one local scale per
    // lag-basis function b, shared by its J exposure coefficients, so a lag
    // window is either active or shrunk to zero as a whole.
    for (k in 1:kl_K) {
      matrix[kl_J, kl_BL] th;
      matrix[kl_J, kl_NL] bl;
      if (kl_hs) {
        vector[kl_BL] lt = sqrt(khs_c2[k] * square(khs_lambda[:, k])
                                ./ (khs_c2[k] + square(khs_tau[k]) * square(khs_lambda[:, k])));
        th = to_matrix(kz[:, k], kl_J, kl_BL) .* rep_matrix((khs_tau[k] * lt)', kl_J);
      } else {
        th = kl_beta_sd * to_matrix(kz[:, k], kl_J, kl_BL);
      }
      bl = th * kl_Lb';
      ktheta[:, k] = to_vector(bl);
      kbeta[:, k]  = bl * rep_vector(1.0, kl_NL);
    }
  }
  for (k in 1:kl_K) kcon[k] = dot_product(kl_erf_v, kbeta[:, k]);
  if (kl_re) {
    if (!kl_re_ctry) kre_c = rep_matrix(0, kl_J, kl_C);
    else if (kl_re_sdc_fix) kre_c = kl_re_sdc_value * kre_zc;
    else kre_c = diag_pre_multiply(kre_sdc, kre_zc);
  }
  kkappa = 0;
  if (kl_use_pg) kkappa = kl_kappa_fixed ? kl_kappa_value : kkappa_free[1];
  if (kl_use_sp) {
    if (kl_sp_dagar)
      kfield = kl_dagar_field(kphi, krho[1], ktau[1], kl_n_par, kl_par_nodes, kl_par_start);
    else
      kfield = ktau[1] * kphi;
  }
  for (c in 1:kl_C) {
    klogit[c] = log(kpi)';
    if (kl_use_pg) klogit[c] = klogit[c] + kkappa * kl_A[c];
    if (kl_use_sp && kl_field_cols > 0)
      klogit[c, 1:kl_field_cols] = klogit[c, 1:kl_field_cols] + kfield[c];
  }
  if (kl_re) {
    if (kl_re_ctry && !kl_re_sdc_fix)
      lprior += normal_lpdf(kre_sdc | 0, kl_re_sd[1]) - kl_J * normal_lccdf(0 | 0, kl_re_sd[1]);
    if (kl_re_adm1)
      lprior += normal_lpdf(kre_sda | 0, kl_re_sd[2]) - kl_J * normal_lccdf(0 | 0, kl_re_sd[2]);
  }
  if (kl_disc) {
    lprior += normal_lpdf(to_vector(kbeta) | 0, kl_beta_sd);   // centre priors: kl_combo_lp
  } else if (!kl_dlnm) {
    if (kl_mu_normal) lprior += normal_lpdf(kmu | kl_mu_mean, kl_mu_sd);
    else lprior += -kl_K * log(kl_hi - kl_lo);      // uniform on the lag range (as brms)
    if (!kl_ksd_fix)
      lprior += normal_lpdf(ksd | 0, kl_sd_sd) - kl_K * normal_lccdf(kl_sd_lb | 0, kl_sd_sd);
    lprior += normal_lpdf(to_vector(kbeta) | 0, kl_beta_sd);
  } else {
    lprior += std_normal_lpdf(to_vector(kz));
    if (kl_hs) {
      lprior += student_t_lpdf(khs_tau | 1, 0, kl_hs_scale_global)
                - kl_K * student_t_lccdf(0 | 1, 0, kl_hs_scale_global);           // half-Cauchy
      lprior += student_t_lpdf(to_vector(khs_lambda) | 1, 0, 1)
                - kl_K * kl_BL * student_t_lccdf(0 | 1, 0, 1);                    // half-Cauchy
      lprior += inv_gamma_lpdf(khs_c2 | 0.5 * kl_hs_slab_df,
                               0.5 * kl_hs_slab_df * square(kl_hs_slab_scale));
    }
  }
  lprior += dirichlet_lpdf(kpi | rep_vector(kl_alpha, kl_K));
  if (kl_use_pg && !kl_kappa_fixed)
    lprior += normal_lpdf(kkappa_free | 0, kl_kappa_sd) - normal_lccdf(0 | 0, kl_kappa_sd);
  if (kl_use_sp)
    lprior += normal_lpdf(ktau | 0, kl_tau_sd) - normal_lccdf(0 | 0, kl_tau_sd);
  if (kl_use_sp && kl_sp_dagar)
    lprior += beta_lpdf(krho | kl_rho_a, kl_rho_b);

  // ---- kernel.lag_re: country lag offsets (non-centred) ----
  kclag = klag_sd[1] * klag_z;
  lprior += normal_lpdf(klag_sd | 0, kl_lag_re_sd) - normal_lccdf(0 | 0, kl_lag_re_sd);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_lin = rep_vector(0.0, N);
    // initialize linear predictor term
    // initialize linear predictor term
    // initialize linear predictor term
    // initialize linear predictor term
    // initialize linear predictor term
    // initialize non-linear predictor term
    nlp_lin += X_lin * b_lin + Zs_lin_1_1 * s_lin_1_1 + offsets_lin;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_lin[n] += r_1_lin_1[J_1[n]] * Z_1_lin_1[n] + r_2_lin_1[J_2[n]] * Z_2_lin_1[n];
    }

    // ---- brms_kl latent classes: country-level mixture ----
    if (kl_disc) {
      // discrete centres: classes per country, then (width, combination) summed out
      matrix[kl_C, kl_K] kl_lsm;
      matrix[kl_re ? N : 0, kl_re ? kl_S * kl_NCt : 0] kl_E;
      if (kl_re)
        kl_E = kl_disc_re_term(kl_KBq, (kre_c[:, kl_ctry] + kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA)[:, kl_adm1])');
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      target += log_sum_exp(kl_qlp + kl_disc_combo_ll(
        kl_disc_country_ll(nlp_lin, kl_KBq, kbeta, kl_Yv, kl_lgY1c, kl_cs, kl_ce, kl_E, kl_re),
        kl_lsm, kl_qidx));
    } else if (kl_rs && kl_pairs) {
      // pair path (kernel.pair_collapse): the row part once, the class part per
      // (exposure row, ADM1) pair, reduce_sum over countries (no row split)
      matrix[kl_C, kl_K] kl_lsm;
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      matrix[kl_re ? 0 : kl_C * kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_C * kl_T, kl_J] kl_BW;
      vector[kl_P] kl_esum = rep_vector(0, kl_P);
      array[kl_C] vector[kl_maxp] kl_esum_pad;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        kl_BW = kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag);
      } else {
        kl_U = kl_lag_re_u(kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag), kbeta);
      }
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      {
        vector[N] kl_el = exp(nlp_lin);
        for (n in 1:N) kl_esum[kl_pix[n]] += kl_el[n];
      }
      for (c in 1:kl_C) {
        kl_esum_pad[c] = rep_vector(0, kl_maxp);
        kl_esum_pad[c][1:(kl_pce[c] - kl_pcs[c] + 1)] = kl_esum[kl_pcs[c]:kl_pce[c]];
      }
      target += dot_product(kl_Yv, nlp_lin) - sum(kl_lgY1c);
      target += reduce_sum(kl_partial_pairs, kl_esum_pad, kl_grain, kl_call, kl_Ypair, kl_pcs, kl_pce,
                           kl_ptc, kl_pa, kl_lsm, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a);
    } else if (kl_rs) {
      // country path: reduce_sum over countries (threaded when compiled with
      // STAN_THREADS), random slopes on the exposure coefficients when kl_re
      matrix[kl_C, kl_K] kl_lsm;
      array[kl_nsmall] vector[kl_maxn] kl_lin_pad;
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      // kernel term once per distinct exposure row, not per row
      matrix[kl_re ? 0 : kl_C * kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_C * kl_T, kl_J] kl_BW;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        kl_BW = kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag);
      } else {
        kl_U = kl_lag_re_u(kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag), kbeta);
      }
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      for (i in 1:kl_nsmall) {
        int c = kl_small[i];
        kl_lin_pad[i] = rep_vector(0, kl_maxn);
        kl_lin_pad[i][1:(kl_ce[c] - kl_cs[c] + 1)] = nlp_lin[kl_cs[c]:kl_ce[c]];
      }
      if (kl_nsmall > 0)
        target += reduce_sum(kl_partial_ctry, kl_lin_pad, kl_grain, kl_small, kl_tixc, kl_Yv, kl_lgY1c,
                             kl_cs, kl_ce, kl_lsm, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a, kl_adm1);
      // large countries: each class log-likelihood summed over row blocks in
      // parallel, then the country mixture
      for (b in 1:kl_nbig) {
        int c = kl_big[b];
        array[kl_ce[c] - kl_cs[c] + 1] real kl_lin_s = to_array_1d(nlp_lin[kl_cs[c]:kl_ce[c]]);
        vector[kl_K] kl_lpk;
        for (k in 1:kl_K)
          kl_lpk[k] = kl_lsm[c, k] - kl_lgY1c[c]
                      + reduce_sum(kl_partial_rows, kl_lin_s, kl_grain_rows, kl_cs[c], k, c,
                                   kl_tixc, kl_Yv, kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a, kl_adm1);
        target += log_sum_exp(kl_lpk);
      }
    } else {
      matrix[kl_C, kl_K] kl_llc = kl_country_ll(
        kl_class_rowll(nlp_lin, kl_lag_re_u(kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag), kbeta), kl_tixc, kl_Yv, kl_lgY1),
        kl_cs, kl_ce);
      for (c in 1:kl_C) target += log_sum_exp(kl_member_logp(klogit[c], kl_fix[c])' + kl_llc[c]');
    }
  }
  // priors including constants
  target += lprior;

  // ---- brms_kl latent classes: random-slope deviations (non-centred) ----
  if (kl_re) target += std_normal_lpdf(to_vector(kre_zc)) + std_normal_lpdf(to_vector(kre_za));
  // ---- brms_kl latent classes: spatial membership prior ----
  if (kl_use_sp && kl_sp_dagar) target += std_normal_lpdf(to_vector(kphi));   // DAGAR innovations
  if (kl_use_sp && !kl_sp_dagar) {                                             // ICAR
    for (k in 1:kl_field_cols) {
      vector[kl_n_comp] csum = rep_vector(0, kl_n_comp);
      target += -0.5 * dot_self(col(kphi, k)[kl_node1] - col(kphi, k)[kl_node2]);
      for (c in 1:kl_C) csum[kl_comp[c]] += kphi[c, k];
      // Soft sum-to-zero per connected component (island singletons -> ~0).
      // sd 0.01 * size, not the usual 0.001: with a mixture likelihood the field
      // is weakly informed and 0.001 made every transition hit max treedepth
      // (synthetic case D); 0.01 removed that and moves the field mean by ~0.01.
      target += normal_lpdf(csum | 0, 0.01 * kl_comp_n);
    }
  }

  // ---- kernel.lag_re: standardised country lag offsets ----
  target += std_normal_lpdf(klag_z);
  target += std_normal_lpdf(zs_lin_1_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {

  // ---- brms_kl latent classes ----
  matrix[kl_K, kl_C] kl_log_resp;               // log class responsibilities
  vector[kl_C] kl_log_lik_ctry;                 // country log-likelihood (LOCO unit)
  vector[kl_gq_row_loglik ? N : 0] log_lik;     // optional row-level log-lik
  vector[kl_disc ? kl_K : 0] klag;              // discrete: drawn class centres (non-decreasing)
  vector[kl_disc ? 1 : 0] kwidth;               // discrete: drawn kernel width
  matrix[kl_disc ? kl_K : 0, kl_NCt] kl_lag_prob;   // discrete: P(centre of class k = i | parameters)

  // ---- brms_kl latent classes ----
  {
    vector[N] nlp_lin = rep_vector(0.0, N);
    matrix[N, kl_K] kl_ll;
    matrix[kl_C, kl_K] kl_llc;
      nlp_lin += X_lin * b_lin + Zs_lin_1_1 * s_lin_1_1 + offsets_lin;
      for (n in 1:N) {
        nlp_lin[n] += r_1_lin_1[J_1[n]] * Z_1_lin_1[n] + r_2_lin_1[J_2[n]] * Z_2_lin_1[n];
      }
    if (kl_disc) {
      // draw (width, combination) given the parameters; class quantities below
      // are conditional on the draw
      matrix[kl_re ? N : 0, kl_re ? kl_J : 0] kl_D;
      matrix[kl_re ? N : 0, kl_re ? kl_S * kl_NCt : 0] kl_E;
      array[kl_K] matrix[kl_C, kl_S * kl_NCt] kl_LLq;
      if (kl_re) {
        kl_D = (kre_c[:, kl_ctry] + kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA)[:, kl_adm1])';
        kl_E = kl_disc_re_term(kl_KBq, kl_D);
      }
      kl_LLq = kl_disc_country_ll(nlp_lin, kl_KBq, kbeta, kl_Yv, kl_lgY1c, kl_cs, kl_ce, kl_E, kl_re);
      matrix[kl_C, kl_K] kl_lsm;
      vector[kl_S * kl_M] kl_w;
      int kl_d;
      int kl_s = 1;
      int kl_m;
      for (c in 1:kl_C) kl_lsm[c] = kl_member_logp(klogit[c], kl_fix[c]);
      kl_w = softmax(kl_qlp + kl_disc_combo_ll(kl_LLq, kl_lsm, kl_qidx));
      kl_d = categorical_rng(kl_w);
      while (kl_d > kl_s * kl_M) kl_s += 1;
      kl_m = kl_d - (kl_s - 1) * kl_M;
      kwidth[1] = kl_width[kl_s];
      kl_lag_prob = rep_matrix(0, kl_K, kl_NCt);
      for (s in 1:kl_S) for (m in 1:kl_M) for (k in 1:kl_K)
        kl_lag_prob[k, kl_combo[m, k]] += kl_w[(s - 1) * kl_M + m];
      for (k in 1:kl_K) {
        klag[k] = kl_cent[kl_combo[kl_m, k]];
        kl_llc[:, k] = col(kl_LLq[k], kl_qidx[kl_s, k, kl_m]);
        if (kl_gq_row_loglik) {
          vector[N] eta = nlp_lin;
          for (j in 1:kl_J) eta += kbeta[j, k] * col(kl_KBq[j], kl_qidx[kl_s, k, kl_m]);
          if (kl_re) eta += col(kl_E, kl_qidx[kl_s, k, kl_m]);
          kl_ll[:, k] = kl_Yv .* eta - exp(eta) - kl_lgY1;
        }
      }
    } else if (kl_rs) {
      // country path (row log-lik not available: kernel.gq_row_loglik refused)
      matrix[kl_re ? kl_J : 0, kl_NA] kl_kre_a;
      matrix[kl_re ? 0 : kl_C * kl_T, kl_re ? 0 : kl_K] kl_U;
      array[kl_re ? kl_K : 0] matrix[kl_C * kl_T, kl_J] kl_BW;
      if (kl_re) {
        kl_kre_a = kl_adm1_dev(kre_sda, kre_za, kl_re_adm1, kl_J, kl_NA);
        kl_BW = kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag);
      } else {
        kl_U = kl_lag_re_u(kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag), kbeta);
      }
      if (kl_pairs) {
        vector[kl_P] kl_esum = rep_vector(0, kl_P);
        vector[N] kl_el = exp(nlp_lin);
        for (n in 1:N) kl_esum[kl_pix[n]] += kl_el[n];
        for (c in 1:kl_C)
          kl_llc[c] = (dot_product(kl_Yv[kl_cs[c]:kl_ce[c]], nlp_lin[kl_cs[c]:kl_ce[c]]) - kl_lgY1c[c]
                       + kl_pair_logp(kl_esum[kl_pcs[c]:kl_pce[c]], kl_Ypair[kl_pcs[c]:kl_pce[c]],
                                      kl_ptc[kl_pcs[c]:kl_pce[c]], kl_pa[kl_pcs[c]:kl_pce[c]], c,
                                      kl_U, kl_BW, kbeta, kl_re, kre_c, kl_kre_a))';
      } else {
        for (c in 1:kl_C) {
          int n = kl_ce[c] - kl_cs[c] + 1;
          matrix[kl_re ? n : 0, kl_re ? kl_J : 0] D;
          if (kl_re) D = (rep_matrix(col(kre_c, c), n) + kl_kre_a[:, kl_adm1[kl_cs[c]:kl_ce[c]]])';
          kl_llc[c] = kl_ctry_logp(nlp_lin[kl_cs[c]:kl_ce[c]], kl_tixc[kl_cs[c]:kl_ce[c]],
                                 kl_Yv[kl_cs[c]:kl_ce[c]], kl_lgY1c[c], rep_row_vector(0, kl_K),
                                 kl_U, kl_BW, kbeta, kl_re, D)';
        }
      }
    } else {
      kl_ll  = kl_class_rowll(nlp_lin, kl_lag_re_u(kl_lag_re_bw(kl_KBl, kl_lags, kmu, kksd, kclag), kbeta), kl_tixc, kl_Yv, kl_lgY1);
      kl_llc = kl_country_ll(kl_ll, kl_cs, kl_ce);
    }
    for (c in 1:kl_C) {
      vector[kl_K] lp = kl_member_logp(klogit[c], kl_fix[c])' + kl_llc[c]';
      kl_log_lik_ctry[c] = log_sum_exp(lp);
      kl_log_resp[:, c] = lp - kl_log_lik_ctry[c];
    }
    if (kl_gq_row_loglik)
      for (n in 1:N) log_lik[n] = log_sum_exp(kl_log_resp[:, kl_ctry[n]] + kl_ll[n]');
  }
}
