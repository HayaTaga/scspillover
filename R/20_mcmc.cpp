// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// =========================== Utilities ===========================
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1 / scale); }

inline double logdet(const arma::mat& A) {
  double sign=0.0, val=0.0;
  arma::log_det(val, sign, A);
  return val;
}

struct SARDims { int T0, N, K, p; };

// N×K: Xt, K×1: beta → N×1
inline arma::vec Xbeta(const arma::mat& Xt, const arma::vec& beta) {
  if (Xt.n_cols == 0) return arma::zeros<arma::vec>(Xt.n_rows);
  return Xt * beta;
}

// N×p: Lambda, p×1: f_t → N×1
inline arma::vec Lf(const arma::mat& Lambda, const arma::vec& f) {
  if (Lambda.n_cols == 0) return arma::zeros<arma::vec>(Lambda.n_rows);
  return Lambda * f;
}

inline double clip1(double z, double lo = 1e-12, double hi = 1e12) {
  return std::min(std::max(z, lo), hi);
}

inline arma::vec clip_vec(const arma::vec& v, double lo = 1e-12, double hi = 1e12) {
  arma::vec out = v;
  for (arma::uword i = 0; i < out.n_elem; ++i) out(i) = clip1(out(i), lo, hi);
  return out;
}


inline void symmetrize_inplace(mat& A) { A = 0.5 * (A + A.t()); }


inline bool robust_chol(mat& R, mat& A, double base_eps=1e-12) {
  symmetrize_inplace(A);
  for (int k=0; k<6; ++k) {
    if (chol(R, A)) return true;
    A.diag() += base_eps * std::pow(10.0, k);
    symmetrize_inplace(A);
  }

  vec d; mat V;
  if (eig_sym(d, V, A)) {
    for (uword i=0;i<d.n_elem;++i) if (d(i) < base_eps) d(i) = base_eps;
    A = V * diagmat(d) * V.t();
    symmetrize_inplace(A);
  }
  return chol(R, A);
}


// ======================== 1) HS-Gibbs for alpha ===================
// Makalic & Schmidt (2015) augmentation
// y: T0, X: T0 x N
// returns: M x N draws (burn discarded)
// [[Rcpp::export]]
arma::mat hs_alpha_gibbs_cpp(const arma::vec& Y0_pre,
                                           const arma::mat& control_outcome_pre,
                                           const int iteration,
                                           const int burn,
                                           const bool verbose=false) {
  Rcpp::RNGScope scope;

  const int T0 = (int)Y0_pre.n_elem;
  const int N  = (int)control_outcome_pre.n_cols;
  const int iters = iteration;

  const double FLO = 1e-12, FHI = 1e12;

  
  vec alpha_curr = 1e-4 * randn<vec>(N);

  arma::vec sx(N);
  for (int j=0; j<N; ++j){
    double sdj = std::sqrt(arma::var(control_outcome_pre.col(j)));
    if(!arma::is_finite(sdj) || sdj < 1e-8) sdj = 1e-8;
    sx(j) = sdj;
  }

  double sy = std::sqrt(arma::var(Y0_pre));
  if(!arma::is_finite(sy) || sy < 1e-8) sy = 1e-8;

  arma::mat X = control_outcome_pre;
  for (int j=0; j<N; ++j) X.col(j) /= sx(j);
  arma::vec y = Y0_pre / sy;

  vec sigma2_i_curr(N, fill::ones);
  vec nu_sigma_i_curr(N, fill::ones);
  double tau2_curr   = 1.0;
  double nu_tau_curr = 1.0;
  double sigma2_curr = as_scalar(var(y));
  if (!arma::is_finite(sigma2_curr) || sigma2_curr <= 0.0) sigma2_curr = 1.0;
  double nu_sigma_curr = 1.0;

  mat XtX = X.t() * X;
  vec Xty = X.t() * y;

  const int M = std::max(0, iters - burn);
  mat draws(M, N, fill::none);

  for (int iter = 1; iter <= iters; ++iter) {
    if (verbose && (iter % 2000 == 0)) {
      Rcpp::Rcout << "Now, iteration = " << iter << "\n";
      Rcpp::checkUserInterrupt();
    }

    // ===== D = D_temp + sigma2 * Diagonal(1 ./ sigma2_i) =====
    vec inv_sigma2_i = 1.0 / clamp(sigma2_i_curr, FLO, FHI);
    double sig2 = std::min(std::max(sigma2_curr, FLO), FHI);
    double tau2 = std::min(std::max(tau2_curr, FLO), FHI);
    // Rcpp::Rcout << sig2 << '\n';
    mat D = XtX;
    D.diag() += sig2 * (inv_sigma2_i);

    // cholesky(D)
    mat R;                                  // 上三角：D = R^T R
    bool ok = robust_chol(R, D);
    if (!ok) Rcpp::stop("chol(D) failed even after stabilization.");

    // ===== D_inv = Symmetric( D_chol \ I ) =====
    mat I_N = eye<mat>(N, N);
    // mat D_inv = solve(trimatu(R), I_N, solve_opts::fast);
    // D_inv = solve(trimatl(R.t()), D_inv, solve_opts::fast);
    mat D_inv = solve(D, I_N, solve_opts::fast);
    symmetrize_inplace(D_inv);

    // ===== alpha_mean = D_inv * X' * Y0,  alpha_cov = sigma2 * D_inv =====
    vec alpha_mean = D_inv * Xty;

    // サンプリング：alpha ~ N(alpha_mean, sigma2 * D_inv)
    // L = sqrt(sigma2) * chol(D_inv)
    // mat Rinv;
    // bool ok2 = robust_chol(Rinv, D_inv);
    // if (!ok2) Rcpp::stop("chol(D_inv) failed even after stabilization.");

    // vec z = randn<vec>(N);
    // vec step = solve(trimatu(Rinv), z, solve_opts::fast);
    // vec alpha_new = alpha_mean + std::sqrt(sig2) * step;
    // alpha_curr = alpha_new;

    arma::mat Sigma = std::max(sig2, 1e-12) * D_inv;  // 共分散行列 Σ = σ² · A⁻¹
    // Rcpp::Rcout << alpha_mean << "\n";
    // Rcpp::Rcout << Sigma << "\n";
    arma::vec alpha_new = arma::mvnrnd(alpha_mean, Sigma, 1);  // 多変量正規 N(μ, Σ) から1サンプル
    alpha_curr = alpha_new;

    // sigma2_i[i] ~ IG(1, 0.5*alpha[i]^2 + 1/nu_sigma_i[i])
    vec sigma2_i_next(N);
    for (int i=0;i<N;++i) {
      double sc = 0.5 * alpha_curr[i]*alpha_curr[i] + 1.0 / std::max(nu_sigma_i_curr[i], FLO);
      double draw = rinvgamma(1.0, sc);
      sigma2_i_next[i] = std::min(std::max(draw, FLO), FHI);
    }
    sigma2_i_curr = sigma2_i_next;

    // nu_sigma_i[i] ~ IG(1, 1/sigma2_i[i] + 1/tau2)
    vec nu_sigma_i_next(N);
    for (int i=0;i<N;++i) {
      double sc = 1.0 / std::max(sigma2_i_curr[i], FLO) + 1.0 / std::max(tau2_curr, FLO);
      double draw = rinvgamma(1.0, sc);
      nu_sigma_i_next[i] = std::min(std::max(draw, FLO), FHI);
    }
    nu_sigma_i_curr = nu_sigma_i_next;

    // tau2 ~ IG((N+1)/2, sum(1./nu_sigma_i) + 1/nu_tau)
    double sc_tau = accu(1.0 / clamp(nu_sigma_i_curr, FLO, FHI)) + 1.0 / std::max(nu_tau_curr, FLO);
    tau2_curr = std::min(std::max(rinvgamma(0.5 * (N + 1.0), sc_tau), FLO), FHI);

    // nu_tau ~ IG(1, 1/tau2 + 1/sigma2)
    double sc_nutau = 1.0 / std::max(tau2_curr, FLO) + 1.0 / std::max(sigma2_curr, FLO);
    nu_tau_curr = std::min(std::max(rinvgamma(1.0, sc_nutau), FLO), FHI);

    // sigma2 ~ IG(1 + T0/2, 1/nu_tau + 1/nu_sigma + 0.5 * SSE)
    vec res = y - X * alpha_curr;
    double sse = dot(res, res);
    double shape_sig = 1.0 + 0.5 * T0;
    double sc_sig = 1.0 / std::max(nu_tau_curr, FLO) + 1.0 / std::max(nu_sigma_curr, FLO) + 0.5 * sse;
    sigma2_curr = std::min(std::max(rinvgamma(shape_sig, sc_sig), FLO), FHI);

    // nu_sigma ~ IG(1, 1/sigma2 + 1/10^2)
    double sc_nus = 1.0 / std::max(sigma2_curr, FLO) + 1.0 / (10.0 * 10.0);
    nu_sigma_curr = std::min(std::max(rinvgamma(1.0, sc_nus), FLO), FHI);

    // keep
    if (iter > burn){
      for (int j=0; j<N; ++j) draws(iter - burn - 1, j) = (sy / sx(j)) * alpha_curr[j];
    }
  }

  return draws;
}

// ===================== 2) Full SAR with factors ===================
//
// Pre-treatment model for controls:
//   A(rho) * y_c,t = rho * w * y0_t + X_t beta + Lambda f_t + eps_t
// where A(rho) = I - rho W
//
// Priors:
//   beta   ~ N(0, c_beta I_K)
//   Lambda ~ N(0, c_lambda I_{N*p})   (row-wise independent)
//   f_t    ~ N(0, I_p)
//   sigma2 ~ IG(a0, b0)
//   rho    ~ Unif((-bnd, bnd)) with spectral bound
//
// Gibbs steps:
//   (1) f_t | rest       : N(m_t, V); V=(I + (1/s2) Lambda'Lambda)^{-1}
//   (2) Lambda_i | rest  : N(m_i, V); V=( (1/s2) sum_t f_t f_t' + (1/c_lambda)I )^{-1}
//   (3) beta | rest      : N(m, V);   V=(X'X/s2 + I/c_beta)^{-1}
//   (4) sigma2 | rest    : IG( (T0*N)/2 + a0, 0.5*sum_t||res_t||^2 + b0 )
//   (5) rho | rest       : RW-MH against full likelihood
//
// Xc_pre_ のレイアウト：((t * N + i) * K + k) でベクトル化
//
// [[Rcpp::export]]
Rcpp::List sar_full_sampler_cpp(const arma::vec& Y0_pre,    // T0
                                const arma::mat& Yc_pre,    // T0 x N
                                Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K), idx=(t*N+i)*K+k
                                int T0, int N, int K, int p,
                                const arma::vec& w,         // N  (adj_vec)
                                const arma::mat& W,         // N x N (adj_mat)
                                int iteration, int burn,
                                double step_rho = 0.01,
                                double c_beta_ridge = 0.0,   // 未使用(HSに切替)。0固定可
                                double c_lambda_ridge = 0.0, // 未使用(階層化に切替)。0固定可
                                double a0 = 1.0, double b0 = 1.0,
                                bool verbose=false) {

  Rcpp::RNGScope scope;

  const int M = std::max(0, iteration - burn);

  // --- data
  arma::vec Y0 = Y0_pre;              // T0
  arma::mat Yc = Yc_pre;              // T0 x N
  arma::vec a  = w;                   // N
  arma::mat A  = W;                   // N x N

  // --- X accessor (idx=(t*N+i)*K+k)
  const bool useX = (K > 0) && Xc_pre_.isNotNull();
  Rcpp::NumericVector Xvec;
  if (useX) Xvec = Xc_pre_.get();
  auto X_get_row = [&](int t)->arma::mat {
    arma::mat Xt(N, K, arma::fill::zeros);
    if (!useX) return Xt;
    for (int i=0;i<N;++i) {
      for (int k=0;k<K;++k) {
        int idx = (t * N + i) * K + k;
        Xt(i,k) = Xvec[idx];
      }
    }
    return Xt;
  };

  // --- spectral bound for rho
  arma::cx_vec evals = arma::eig_gen(A);
  double maxabs = 0.0;
  for (arma::uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // --- storage
  arma::vec rho_draws(M, arma::fill::none);
  arma::vec s2_draws(M, arma::fill::none);
  arma::mat beta0_draws(M, K, arma::fill::zeros);
  arma::cube Lambda_draws(N, p, M, arma::fill::zeros); // = Eta
  arma::cube F_draws(p, T0, M, arma::fill::zeros);     // = Gamma

  // --- states
  double rho = 0.0;       // lam は Julia に倣い rho と同一視
  double s2  = 1.0;       // observation variance

  // beta0 ~ HS（Makalic-Schmidt）
  arma::vec beta0 = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());
  arma::vec sig2_b0 = (K>0 ? arma::ones<arma::vec>(K) : arma::vec());     // λ_j^2
  arma::vec nu_sig_b0 = (K>0 ? arma::ones<arma::vec>(K) : arma::vec());   // ν_j
  double tau2_b0 = 1.0, nu_tau_b0 = 1.0;                                  // τ^2, ν_τ

  // sigma^2 のハイパー（Julia の層）
  double nu_sigma2 = 1.0;

  // 因子: Eta (N x p), Gamma (p x T0)
  arma::mat Eta   = (p>0 ? arma::zeros<arma::mat>(N,p) : arma::mat());
  arma::mat Gamma = (p>0 ? arma::zeros<arma::mat>(p,T0) : arma::mat());
  double phi_g = 0.0, s2_g = 1.0, nu_s2_g = 1.0;        // gamma_t の AR(1) 事前
  arma::vec omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());     // shrink for Eta
  arma::vec nu_omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  double s2_eta = 1.0, nu_s2_eta = 1.0;

  arma::mat I_N = arma::eye(N,N);
  arma::mat I_K = (K>0 ? arma::eye(K,K) : arma::mat());
  arma::mat I_p = (p>0 ? arma::eye(p,p) : arma::mat());

  auto ll_rho = [&](double r)->double {
    if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
    arma::mat Mmat = I_N - r * A;
    double ldet = logdet(Mmat);
    if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    double ss = 0.0;
    for (int t=0;t<T0;++t) {
      arma::vec mu = r * a * Y0[t];                   // lam = rho
      if (useX) mu += X_get_row(t) * beta0;
      if (p>0) mu += Eta * Gamma.col(t);
      arma::vec u = Mmat * Yc.row(t).t() - mu;
      ss += arma::dot(u,u);
    }
    return T0 * ldet - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
  };

  int iters = M + burn;
  int acc = 0;

  for (int it=0; it<iters; ++it) {

    // ===== (1) gamma_t | rest  （Julia 準拠：AR(1) prior ベースの簡便更新） =====
    if (p>0) {
      // データ込みの厳密条件は重いので、Julia コードに倣い
      // V_g = (Eta'Eta/s2 + I/s2_g)^(-1),  m_t = V_g * (Eta' r_t)/s2 で近似
      arma::mat EtE = Eta.t() * Eta;                      // p x p
      arma::mat Vg  = arma::inv_sympd( EtE / s2 + I_p / s2_g );
      arma::mat Lg  = arma::chol(Vg, "lower");
      for (int t=0; t<T0; ++t) {
        arma::vec r = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
        if (useX) r -= X_get_row(t) * beta0;
        arma::vec mg = Vg * (Eta.t() * r / s2);
        arma::vec z  = arma::randn<arma::vec>(p);
        Gamma.col(t) = mg + Lg * z;
      }
      // phi_gamma ~ truncNorm, s2_g ~ IG, nu_s2_g ~ IG（Juliaの式）
      double den = 0.0, num = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::vec(p, arma::fill::zeros) : arma::vec(Gamma.col(t-1));
        den += arma::dot(gl, gl);
        num += arma::dot(gl, Gamma.col(t));
      }
      double mean_phi = (den>0? num/den : 0.0);
      double var_phi  = (den>0? s2_g/den : 1.0);
      double cand;
      do { cand = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand) > 1.0);
      phi_g = cand;

      double sc = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::vec(p, arma::fill::zeros) : arma::vec(Gamma.col(t-1));
        arma::vec diff = Gamma.col(t) - phi_g * gl;
        sc += 0.5 * arma::dot(diff, diff);
      }
      s2_g     = rinvgamma(0.5 + 0.5 * p * T0, sc + 1.0/clip1(nu_s2_g));
      nu_s2_g  = rinvgamma(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
    }

    // ===== (2) Eta | rest  （行ごと回帰＋階層縮退） =====
    if (p>0) {
      arma::mat GtG = Gamma * Gamma.t();                    // p x p
      arma::mat Domega = arma::diagmat(omega_k);            // p x p
      arma::mat Vrow = arma::inv_sympd( GtG / s2 + Domega / clip1(s2_eta) );
      arma::mat Lrow = arma::chol(Vrow, "lower");
      for (int i=0; i<N; ++i) {
        arma::vec rhs = arma::zeros<arma::vec>(p);
        for (int t=0; t<T0; ++t) {
          arma::vec r = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
          if (useX) r -= X_get_row(t) * beta0;
          rhs += Gamma.col(t) * r(i);
        }
        arma::vec m = Vrow * (rhs / s2);
        arma::vec z = arma::randn<arma::vec>(p);
        Eta.row(i) = (m + Lrow * z).t();
      }
      double sc = 0.0;
      for (int i=0;i<N;++i) {
        arma::vec ei = Eta.row(i).t();
        sc += arma::dot(ei, Domega * ei);
      }
      s2_eta    = rinvgamma(0.5 + 0.5 * p * N, 0.5*sc + 1.0/clip1(nu_s2_eta));
      nu_s2_eta = rinvgamma(1.0, 1.0/clip1(s2_eta) + 1.0/100.0);
      for (int k=0;k<p;++k) {
        double tmp=0.0; for (int i=0;i<N;++i) tmp += 0.5 * Eta(i,k)*Eta(i,k) / clip1(s2_eta);
        double rate_ok = 1.0/clip1(nu_omega_k(k)) + tmp;
        omega_k(k)     = rinvgamma(0.5*(N+1.0), rate_ok);
        nu_omega_k(k)  = rinvgamma(1.0, 1.0 + 1.0/clip1(omega_k(k)));
      }
    }

    // ===== (3) beta0 | rest  （Horseshoe / Makalic–Schmidt） =====
    if (useX) {
      arma::mat Ab = arma::zeros<arma::mat>(K,K);
      arma::vec Bb = arma::zeros<arma::vec>(K);
      for (int t=0; t<T0; ++t) {
        arma::mat Xt = X_get_row(t);          // N x K
        Ab += Xt.t() * Xt;
        arma::vec Btmp = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
        if (p>0) Btmp -= Eta * Gamma.col(t);
        Bb += Xt.t() * Btmp;
      }
      // A_beta0 = Σ X'X + σ² * Diag(1/σ²_{β0})
      Ab.diag() += clip1(s2) * (1.0 / clip_vec(sig2_b0));
      arma::mat Ainv = arma::inv_sympd(Ab);
      arma::vec m    = Ainv * Bb;
      arma::mat S    = clip1(s2) * Ainv;
      beta0 = arma::mvnrnd(m, 0.5*(S+S.t()), 1);

      // HS のスケール群を更新
      for (int j=0;j<K;++j) {
        double rate_l = 0.5 * beta0(j)*beta0(j) + 1.0/clip1(nu_sig_b0(j));
        sig2_b0(j)    = rinvgamma(1.0, rate_l);
        double rate_nu= 1.0/clip1(sig2_b0(j)) + 1.0/clip1(tau2_b0);
        nu_sig_b0(j)  = rinvgamma(1.0, rate_nu);
      }
      double sum_inv_nu = 0.0; for (int j=0;j<K;++j) sum_inv_nu += 1.0/clip1(nu_sig_b0(j));
      tau2_b0    = rinvgamma(1.0, 1.0/clip1(nu_tau_b0) + sum_inv_nu);
      nu_tau_b0  = rinvgamma(1.0, 1.0/clip1(tau2_b0) + 1.0/clip1(s2));
    }

    // ===== (4) sigma^2 | rest  （nu_sigma2 の層を含む） =====
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec mu = rho * a * Y0[t];
      if (useX) mu += X_get_row(t) * beta0;
      if (p>0)  mu += Eta * Gamma.col(t);
      arma::vec u = (I_N - rho*A) * Yc.row(t).t() - mu;
      ss += arma::dot(u,u);
    }
    double shape = 1.0 + 0.5 * (T0 * N);
    double rate  = 1.0/clip1(nu_sigma2) + 0.5 * ss;
    // Rcpp::Rcout << "[sigma2] rate=" << rate << "\n";
    s2 = rinvgamma(shape, rate);
    // Rcpp::Rcout << "[sigma2] draw=" << s2 << "\n";
    nu_sigma2 = rinvgamma(1.0, 1.0/clip1(s2) + 1.0/100.0);

    // ===== (5) rho | rest  （RW-MH, lam=rho） =====
    auto ll = [&](double r)->double {
      if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
      arma::mat Mmat = I_N - r * A;
      double ldet = logdet(Mmat);
      if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
      double ss2 = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec mu = r * a * Y0[t];
        if (useX) mu += X_get_row(t) * beta0;
        if (p>0)  mu += Eta * Gamma.col(t);
        arma::vec u = Mmat * Yc.row(t).t() - mu;
        ss2 += arma::dot(u,u);
      }
      return T0 * ldet - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss2 / s2;
    };
    double prop = R::rnorm(rho, step_rho);
    double lcur = ll(rho);
    double lprp = ll(prop);
    if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) { rho = prop; if (it>=burn) acc++; }

    // ===== store after burn =====
    if (it > burn) {
      int m = it - burn;
      rho_draws[m] = rho;
      s2_draws[m]  = s2;
      if (K>0) beta0_draws.row(m) = beta0.t();
      if (p>0) {
        Lambda_draws.slice(m) = Eta;     // 互換：Lambda := Eta
        F_draws.slice(m)      = Gamma;   // 互換：F      := Gamma
      }
    }

    if (verbose && (it % 2000 == 0)) Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    _["rho"]     = rho_draws,
    _["beta"]    = beta0_draws,     // HS(beta0). K==0なら空行列
    _["sigma2"]  = s2_draws,
    _["Lambda"]  = Lambda_draws,    // = Eta
    _["F"]       = F_draws,         // = Gamma
    _["acc_rate"]= acc / std::max(1, M)
  );
}


// [[Rcpp::export]]
Rcpp::List sar_full_one_step_cpp(const arma::vec& Y0_pre,    // T0
                                 const arma::mat& Yc_pre,    // T0 x N
                                 Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K)
                                 int T0, int N, int K, int p,
                                 const arma::vec& w,         // N
                                 const arma::mat& W,         // N x N
                                 Rcpp::List state,           // 現在の状態（下で仕様化）
                                 double step_rho = 0.01,
                                 double a0 = 1.0, double b0 = 1.0,
                                 bool verbose=false) {
  Rcpp::RNGScope scope;

  // ---- unpack state (必須ブロック) ----
  double rho    = Rcpp::as<double>(state["rho"]);
  double s2     = Rcpp::as<double>(state["sigma2"]);

  arma::vec beta = (K>0 && state.containsElementNamed("beta")) ? 
                    Rcpp::as<arma::vec>(state["beta"]) : arma::vec();

  arma::mat Eta   = (p>0 && state.containsElementNamed("Lambda")) ? 
                    Rcpp::as<arma::mat>(state["Lambda"]) : arma::mat();

  arma::mat Gamma = (p>0 && state.containsElementNamed("F")) ? 
                    Rcpp::as<arma::mat>(state["F"]) : arma::mat();

  // ---- unpack auxiliaries (存在しなければデフォルト) ----
  arma::vec sig2_b0   = (K>0 && state.containsElementNamed("sig2_b0"))  ? Rcpp::as<arma::vec>(state["sig2_b0"])  : arma::ones<arma::vec>(K);
  arma::vec nu_sig_b0 = (K>0 && state.containsElementNamed("nu_sig_b0"))? Rcpp::as<arma::vec>(state["nu_sig_b0"]): arma::ones<arma::vec>(K);
  double tau2_b0      = (K>0 && state.containsElementNamed("tau2_b0"))  ? Rcpp::as<double>(state["tau2_b0"])     : 1.0;
  double nu_tau_b0    = (K>0 && state.containsElementNamed("nu_tau_b0"))? Rcpp::as<double>(state["nu_tau_b0"])   : 1.0;

  double nu_sigma2    = state.containsElementNamed("nu_sigma2") ? Rcpp::as<double>(state["nu_sigma2"]) : 1.0;

  double phi_g        = (p>0 && state.containsElementNamed("phi_g"))     ? Rcpp::as<double>(state["phi_g"])     : 0.0;
  double s2_g         = (p>0 && state.containsElementNamed("s2_g"))      ? Rcpp::as<double>(state["s2_g"])      : 1.0;
  double nu_s2_g      = (p>0 && state.containsElementNamed("nu_s2_g"))   ? Rcpp::as<double>(state["nu_s2_g"])   : 1.0;

  arma::vec omega_k   = (p>0 && state.containsElementNamed("omega_k"))   ? Rcpp::as<arma::vec>(state["omega_k"]) : arma::ones<arma::vec>(p);
  arma::vec nu_omega_k= (p>0 && state.containsElementNamed("nu_omega_k"))? Rcpp::as<arma::vec>(state["nu_omega_k"]) : arma::ones<arma::vec>(p);
  double s2_eta       = (p>0 && state.containsElementNamed("s2_eta"))    ? Rcpp::as<double>(state["s2_eta"])    : 1.0;
  double nu_s2_eta    = (p>0 && state.containsElementNamed("nu_s2_eta")) ? Rcpp::as<double>(state["nu_s2_eta"]) : 1.0;

  // ---- X のアクセサ ----
  const bool useX = (K > 0) && Xc_pre_.isNotNull();
  Rcpp::NumericVector Xvec;
  if (useX) Xvec = Xc_pre_.get();
  auto X_get_row = [&](int t)->arma::mat {
    arma::mat Xt(N, K, arma::fill::zeros);
    if (!useX) return Xt;
    for (int i=0;i<N;++i) {
      for (int k=0;k<K;++k) {
        int idx = (t * N + i) * K + k;
        Xt(i,k) = Xvec[idx];
      }
    }
    return Xt;
  };

  // ---- 便利な定数など ----
  arma::vec Y0 = Y0_pre;
  arma::mat Yc = Yc_pre;
  arma::vec a  = w;
  arma::mat A  = W;

  arma::mat I_N = arma::eye(N,N);
  arma::mat I_K = (K>0 ? arma::eye(K,K) : arma::mat());
  arma::mat I_p = (p>0 ? arma::eye(p,p) : arma::mat());

  // ---- spectral bound for rho ----
  arma::cx_vec evals = arma::eig_gen(A);
  double maxabs = 0.0;
  for (arma::uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // ============================================================
  //  (1) gamma_t | rest  （p>0 の場合）
  // ============================================================
  if (p>0) {
    arma::mat EtE = Eta.t() * Eta;                      // p x p
    arma::mat Vg  = arma::inv_sympd( EtE / s2 + I_p / s2_g );
    arma::mat Lg  = arma::chol(Vg, "lower");
    for (int t=0; t<T0; ++t) {
      arma::vec r = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
      if (useX) r -= X_get_row(t) * beta;
      arma::vec mg = Vg * (Eta.t() * r / s2);
      arma::vec z  = arma::randn<arma::vec>(p);
      Gamma.col(t) = mg + Lg * z;
    }
    // AR(1) hyper: phi_g (truncNorm), s2_g (IG), nu_s2_g (IG)
    double den = 0.0, num = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec gl = (t==0) ? arma::vec(p, arma::fill::zeros) : arma::vec(Gamma.col(t-1));
      den += arma::dot(gl, gl);
      num += arma::dot(gl, Gamma.col(t));
    }
    double mean_phi = (den>0? num/den : 0.0);
    double var_phi  = (den>0? s2_g/den : 1.0);
    double cand;
    do { cand = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand) > 1.0);
    phi_g = cand;

    double sc = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec gl = (t==0) ? arma::vec(p, arma::fill::zeros) : arma::vec(Gamma.col(t-1));
      arma::vec diff = Gamma.col(t) - phi_g * gl;
      sc += 0.5 * arma::dot(diff, diff);
    }
    s2_g     = rinvgamma(0.5 + 0.5 * p * T0, sc + 1.0/std::max(nu_s2_g, 1e-12));
    nu_s2_g  = rinvgamma(1.0, 1.0/std::max(s2_g,1e-12) + 1.0/100.0);
  }

  // ============================================================
  //  (2) Eta | rest  （p>0 の場合）
  // ============================================================
  if (p>0) {
    arma::mat GtG = Gamma * Gamma.t();                    // p x p
    arma::mat Domega = arma::diagmat(omega_k);            // p x p
    arma::mat Vrow = arma::inv_sympd( GtG / s2 + Domega / std::max(s2_eta,1e-12) );
    arma::mat Lrow = arma::chol(Vrow, "lower");
    for (int i=0; i<N; ++i) {
      arma::vec rhs = arma::zeros<arma::vec>(p);
      for (int t=0; t<T0; ++t) {
        arma::vec r = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
        if (useX) r -= X_get_row(t) * beta;
        rhs += Gamma.col(t) * r(i);
      }
      arma::vec m = Vrow * (rhs / s2);
      arma::vec z = arma::randn<arma::vec>(p);
      Eta.row(i) = (m + Lrow * z).t();
    }
    double sc = 0.0;
    for (int i=0;i<N;++i) {
      arma::vec ei = Eta.row(i).t();
      sc += arma::dot(ei, Domega * ei);
    }
    s2_eta    = rinvgamma(0.5 + 0.5 * p * N, 0.5*sc + 1.0/std::max(nu_s2_eta,1e-12));
    nu_s2_eta = rinvgamma(1.0, 1.0/std::max(s2_eta,1e-12) + 1.0/100.0);
    for (int k=0;k<p;++k) {
      double tmp=0.0; for (int i=0;i<N;++i) tmp += 0.5 * Eta(i,k)*Eta(i,k) / std::max(s2_eta,1e-12);
      double rate_ok = 1.0/std::max(nu_omega_k(k),1e-12) + tmp;
      omega_k(k)     = rinvgamma(0.5*(N+1.0), rate_ok);
      nu_omega_k(k)  = rinvgamma(1.0, 1.0 + 1.0/std::max(omega_k(k),1e-12));
    }
  }

  // ============================================================
  //  (3) beta | rest  （K>0 の場合; Horseshoe 補助も更新）
  // ============================================================
  if (useX) {
    arma::mat Ab = arma::zeros<arma::mat>(K,K);
    arma::vec Bb = arma::zeros<arma::vec>(K);
    for (int t=0; t<T0; ++t) {
      arma::mat Xt = X_get_row(t);
      Ab += Xt.t() * Xt;
      arma::vec Btmp = (I_N - rho*A) * Yc.row(t).t() - rho * a * Y0[t];
      if (p>0) Btmp -= Eta * Gamma.col(t);
      Bb += Xt.t() * Btmp;
    }
    Ab.diag() += std::max(s2,1e-12) * (1.0 / arma::clamp(sig2_b0, 1e-12, 1e12));
    arma::mat Ainv = arma::inv_sympd(Ab);
    arma::vec m    = Ainv * Bb;
    arma::mat S    = std::max(s2,1e-12) * Ainv;
    beta           = arma::mvnrnd(m, 0.5*(S+S.t()), 1);

    // HS 補助を更新（Makalic–Schmidt）
    for (int j=0;j<K;++j) {
      double rate_l = 0.5 * beta(j)*beta(j) + 1.0/std::max(nu_sig_b0(j),1e-12);
      sig2_b0(j)    = rinvgamma(1.0, rate_l);
      double rate_nu= 1.0/std::max(sig2_b0(j),1e-12) + 1.0/std::max(tau2_b0,1e-12);
      nu_sig_b0(j)  = rinvgamma(1.0, rate_nu);
    }
    double sum_inv_nu = 0.0; for (int j=0;j<K;++j) sum_inv_nu += 1.0/std::max(nu_sig_b0(j),1e-12);
    tau2_b0    = rinvgamma(1.0, 1.0/std::max(nu_tau_b0,1e-12) + sum_inv_nu);
    nu_tau_b0  = rinvgamma(1.0, 1.0/std::max(tau2_b0,1e-12) + 1.0/std::max(s2,1e-12));
  }

  // ============================================================
  //  (4) sigma^2 | rest  （nu_sigma2 層も更新）
  // ============================================================
  double ss = 0.0;
  for (int t=0; t<T0; ++t) {
    arma::vec mu = rho * a * Y0[t];
    if (useX) mu += X_get_row(t) * beta;
    if (p>0)  mu += Eta * Gamma.col(t);
    arma::vec u = (I_N - rho*A) * Yc.row(t).t() - mu;
    ss += arma::dot(u,u);
  }
  double shape = a0 + 0.5 * (T0 * N);
  double rate  = 1.0/std::max(nu_sigma2,1e-12) + 0.5 * ss;
  s2 = rinvgamma(shape, rate);
  nu_sigma2 = rinvgamma(1.0, 1.0/std::max(s2,1e-12) + 1.0/100.0);

  // ============================================================
  //  (5) rho | rest  （RW-MH）
  // ============================================================
  auto ll = [&](double r)->double {
    if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
    arma::mat Mmat = I_N - r * A;
    double ldet = logdet(Mmat);
    if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    double ss2 = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec mu = r * a * Y0[t];
      if (useX) mu += X_get_row(t) * beta;
      if (p>0)  mu += Eta * Gamma.col(t);
      arma::vec u = Mmat * Yc.row(t).t() - mu;
      ss2 += arma::dot(u,u);
    }
    return T0 * ldet - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss2 / s2;
  };
  double prop = R::rnorm(rho, step_rho);
  double lcur = ll(rho);
  double lprp = ll(prop);
  bool accepted = false;
  if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) { rho = prop; accepted = true; }

  // ---- return updated state ----
  Rcpp::List out;

  // 必須フィールド
  out["rho"]      = rho;
  out["sigma2"]   = s2;
  out["acc_rho"]  = accepted;

  // beta ブロック（K 次元）
  if (K > 0) {
    out["beta"]      = Rcpp::wrap(beta);
    out["sig2_b0"]   = Rcpp::wrap(sig2_b0);
    out["nu_sig_b0"] = Rcpp::wrap(nu_sig_b0);
  } else {
    out["beta"]      = Rcpp::NumericVector(0);
    out["sig2_b0"]   = Rcpp::NumericVector(0);
    out["nu_sig_b0"] = Rcpp::NumericVector(0);
  }
  out["tau2_b0"]   = tau2_b0;
  out["nu_tau_b0"] = nu_tau_b0;

  // sigma^2 のハイパー
  out["nu_sigma2"] = nu_sigma2;

  // 因子・荷重ブロック（p 次元）
  if (p > 0) {
    out["Lambda"]     = Rcpp::wrap(Eta);     // N x p
    out["F"]          = Rcpp::wrap(Gamma);   // p x T0
    out["phi_g"]      = phi_g;
    out["s2_g"]       = s2_g;
    out["nu_s2_g"]    = nu_s2_g;
    out["omega_k"]    = Rcpp::wrap(omega_k);
    out["nu_omega_k"] = Rcpp::wrap(nu_omega_k);
    out["s2_eta"]     = s2_eta;
    out["nu_s2_eta"]  = nu_s2_eta;
  } else {
    out["Lambda"]     = Rcpp::NumericMatrix(N, 0);
    out["F"]          = Rcpp::NumericMatrix(0, T0);
    out["phi_g"]      = 0.0;
    out["s2_g"]       = 1.0;
    out["nu_s2_g"]    = 1.0;
    out["omega_k"]    = Rcpp::NumericVector(0);
    out["nu_omega_k"] = Rcpp::NumericVector(0);
    out["s2_eta"]     = 1.0;
    out["nu_s2_eta"]  = 1.0;
  }

  return out;
}