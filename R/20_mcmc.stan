data {
  int<lower=1> T0;                  // pre periods
  int<lower=1> N;                   // donors
  int<lower=1> K;                   // covariates per unit per period

  matrix[T0, N] Yc_raw;             // outcomes (unscaled, row=t, col=i)
  array[T0] matrix[N, K] X;         // <-- 3次元はNG。配列(時点)×行列(N×K)にする

  vector[N] w_orig;                 // spatial loadings of treated on donors
  matrix[N, N] W;                   // spatial weights among donors
  real<lower=0> bnd;                // rho stability bound

  // horseshoe hyperparams
  real<lower=0> scale_global_alpha;
  real<lower=0> scale_global_beta;
  real<lower=0> nu_global;          // typically 1
  real<lower=0> nu_local;           // typically 1

  real<lower=0> sigma_prior_scale;  // prior scale for residual sd
}

transformed data {
  matrix[N, N] I_N = diag_matrix(rep_vector(1.0, N));
  vector[N] w = w_orig / sqrt(dot_self(w_orig));

  // scale Yc by unit-specific sd (列ごと)
  matrix[T0, N] Yc;
  vector[N] col_sd;
  for (i in 1:N) {
    real sd_i = sd(Yc_raw[, i]);
    if (is_nan(sd_i) || is_inf(sd_i) || (sd_i < 1e-8)) sd_i = 1.0;
    col_sd[i] = sd_i;
    Yc[, i] = Yc_raw[, i] / sd_i;
  }
}

parameters {
  real rho_raw;                       // unconstrained -> tanh -> (-bnd, bnd)
  vector[N] alpha;                    // synthetic weights
  vector[K] beta0;                    // coefficients on X
  real<lower=0> sigma;                // residual sd

  // horseshoe for alpha
  vector<lower=0>[N] lambda_alpha;    // local scales
  real<lower=0> tau_alpha;            // global scale

  // horseshoe for beta0
  vector<lower=0>[K] lambda_beta0;
  real<lower=0> tau_beta0;
}

transformed parameters {
  real rho = bnd * tanh(rho_raw);
  matrix[N, N] Mtilde = I_N - rho * W - rho * w * alpha';
  real log_det_Mtilde = log_determinant(Mtilde);  // LUベース、非対称でもOK
}

model {
  // priors
  rho_raw ~ std_normal();

  lambda_alpha ~ student_t(nu_local, 0, 1);
  tau_alpha    ~ student_t(nu_global, 0, scale_global_alpha);
  {
    vector[N] sigma_alpha_eff = lambda_alpha * tau_alpha;
    alpha ~ normal(0, sigma_alpha_eff);
  }

  lambda_beta0 ~ student_t(nu_local, 0, 1);
  tau_beta0    ~ student_t(nu_global, 0, scale_global_beta);
  {
    vector[K] sigma_beta0_eff = lambda_beta0 * tau_beta0;
    beta0 ~ normal(0, sigma_beta0_eff);
  }

  sigma ~ student_t(nu_global, 0, sigma_prior_scale);

  // likelihood
  target += T0 * log_det_Mtilde;  // SARのヤコビアン項

  for (t in 1:T0) {
    matrix[N, K] Xt = X[t];                 // 時点tの N×K 行列
    vector[N] u_t = Mtilde * to_vector(Yc[t]) - Xt * beta0;  // Yc[t]はrow_vector -> vectorへ
    target += normal_lpdf(u_t | 0, sigma);
  }
}

generated quantities {
  real rho_constrained = rho;
  vector[N] alpha_unscaled;
  for (i in 1:N)
    alpha_unscaled[i] = alpha[i] / col_sd[i];
}