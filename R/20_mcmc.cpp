// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// =========================== Utilities ===========================
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1 / scale); }

inline double logdet_signed(const arma::mat& A) {
  double sign=0.0, val=0.0;
  arma::log_det(val, sign, A);
  if (sign <= 0.0 || !arma::is_finite(val)) return -std::numeric_limits<double>::infinity();
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
                                Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K)
                                int T0, int N, int K, int p,
                                const arma::vec& w,         // N
                                const arma::mat& W,         // N x N
                                int M, int burn,
                                double step_rho = 0.01,
                                double c_beta = 10.0,
                                double c_lambda = 10.0,
                                double a0 = 1.0, double b0 = 1.0,
                                bool verbose=false) {

  Rcpp::RNGScope scope;
  // --- data
  vec Y0 = conv_to<vec>::from(Y0_pre);
  mat Yc = conv_to<mat>::from(Yc_pre);
  vec ww = conv_to<vec>::from(w);
  mat WW = conv_to<mat>::from(W);

  // --- X accessor
  const bool useX = (K > 0) && Xc_pre_.isNotNull();
  Rcpp::NumericVector Xvec;
  if (useX) Xvec = Xc_pre_.get();

  auto X_get_row = [&](int t)->arma::mat{
    mat Xt(N, K, fill::zeros);
    if (!useX) return Xt;
    for (int i=0;i<N;++i){
      for (int k=0;k<K;++k){
        int idx = (t * N + i) * K + k;
        Xt(i,k) = Xvec[idx];
      }
    }
    return Xt;
  };

  // --- spectral bound for rho
  cx_vec evals = eig_gen(WW);
  double maxabs = 0.0; for (uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // --- storage
  vec rho_draws(M, fill::none);
  mat beta_draws(M, K, fill::zeros);
  vec s2_draws(M, fill::none);
  cube Lambda_draws(N, p, M, fill::zeros); // optional
  cube F_draws(p, T0, M, fill::zeros);     // optional

  // --- states
  double rho = 0.0;
  double s2  = 1.0;
  vec beta   = zeros<vec>(std::max(0, K));
  mat Lambda = (p>0 ? zeros<mat>(N, p) : mat());
  mat F      = (p>0 ? zeros<mat>(p, T0) : mat());
  int acc = 0;

  int iters = M + burn;

  // --- helpers (cached)
  mat I_N = eye(N,N);
  mat I_p = (p>0 ? eye(p,p) : mat());
  mat I_K = (K>0 ? eye(K,K) : mat());

  auto stack_design = [&](double rho, mat& A, mat& Xstack, mat& Fstack, vec& ystar){
    A = I_N - rho * WW; // N x N
    ystar.set_size(T0 * N);
    if (useX) Xstack.set_size(T0 * N, K); else Xstack.reset();
    if (p>0)  Fstack.set_size(T0 * N, p); else Fstack.reset();

    for (int t=0; t<T0; ++t) {
      vec yc = Yc.row(t).t();
      vec base = A * yc - rho * ww * Y0[t]; // N
      ystar.subvec(t*N, t*N+N-1) = base;

      if (useX) {
        mat Xt = X_get_row(t);               // N x K
        Xstack.rows(t*N, t*N+N-1) = Xt;
      }
      if (p>0) {
        mat Ft = repmat(F.col(t).t(), N, 1); // N x p (each row = f_t')
        Fstack.rows(t*N, t*N+N-1) = Ft;
      }
    }
  };

  auto cond_loglik_rho = [&](double r)->double{
    if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
    mat A = I_N - r * WW;
    double logdetA = logdet_signed(A);
    if (!is_finite(logdetA)) return -std::numeric_limits<double>::infinity();

    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      vec yc = Yc.row(t).t();
      vec rhs = r * ww * Y0[t];
      if (useX) {
        mat Xt = X_get_row(t);
        rhs += Xt * beta;
      }
      if (p>0) rhs += Lf(Lambda, F.col(t));
      vec u = A * yc - rhs;
      ss += dot(u,u);
    }
    // exact Gaussian likelihood with sigma2=s2
    double ll = T0 * logdetA - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
    return ll;
  };

  for (int it=0; it<iters; ++it) {
    // ---------- Stack residual base ----------
    mat A, Xstack, Fstack;
    vec ystar;
    stack_design(rho, A, Xstack, Fstack, ystar); // A Yc - rho w Y0

    // ---------- (1) sample F | rest ----------
    if (p>0) {
      // V_f = (I + (1/s2) Lambda'Lambda)^-1 ; m_t = V_f * (1/s2) Lambda' (ystar_block - X_t beta)
      mat LtL = Lambda.t() * Lambda;                      // p x p
      mat Vf  = inv_sympd(I_p + LtL / s2);               // p x p
      mat Lty = zeros<mat>(p, T0);
      for (int t=0; t<T0; ++t) {
        vec base = ystar.subvec(t*N, t*N+N-1);           // N
        if (useX) {
          mat Xt = Xstack.rows(t*N, t*N+N-1);
          base -= Xt * beta;
        }
        // base ≈ Lambda f_t + eps
        Lty.col(t) = Lambda.t() * base;
      }
      for (int t=0; t<T0; ++t) {
        vec mt = Vf * (Lty.col(t) / s2);
        vec z  = randn<vec>(p);
        F.col(t) = mt + chol(Vf, "lower") * z;
      }
    }

    // ---------- (2) sample Lambda | rest ----------
    if (p>0) {
      // V_L = ( (1/s2) sum_t f_t f_t' + (1/c_lambda) I )^-1
      mat FtF = F * F.t();                                 // p x p
      mat Vrow = inv_sympd( FtF / s2 + I_p / c_lambda );   // p x p
      mat cholV = chol(Vrow, "lower");
      // row-wise regression for each i
      for (int i=0; i<N; ++i) {
        // rhs_i = sum_t f_t * r_{ti}, where r_{t} = ystar_block - X_t beta
        vec rhs = zeros<vec>(p);
        for (int t=0; t<T0; ++t) {
          vec base = ystar.subvec(t*N, t*N+N-1);          // N
          if (useX) {
            mat Xt = Xstack.rows(t*N, t*N+N-1);
            base -= Xt * beta;
          }
          double rti = base[i];
          rhs += F.col(t) * rti;
        }
        vec mi = Vrow * (rhs / s2);
        vec zi = randn<vec>(p);
        Lambda.row(i) = (mi + cholV * zi).t();
      }
    }

    // ---------- (3) sample beta | rest ----------
    if (useX) {
      // ystar - Fstack * vec(F) = Xstack * beta + eps
      vec ytilde = ystar;
      if (p>0) {
        // subtract Lambda f_t per block
        for (int t=0; t<T0; ++t) {
          vec Lft = Lambda * F.col(t);                     // N
          ytilde.subvec(t*N, t*N+N-1) -= Lft;
        }
      }
      mat XtX = Xstack.t() * Xstack;
      vec Xty = Xstack.t() * ytilde;
      mat Vb  = inv_sympd(XtX / s2 + I_K / c_beta);
      vec mb  = Vb * (Xty / s2);
      vec zb  = randn<vec>(K);
      beta    = mb + chol(Vb, "lower") * zb;
    }

    // ---------- (4) sample sigma2 | rest ----------
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      vec yc = Yc.row(t).t();
      vec mu = rho * ww * Y0[t];
      if (useX) {
        mat Xt = X_get_row(t);
        mu += Xt * beta;
      }
      if (p>0) mu += Lambda * F.col(t);
      vec u = (I_N - rho * WW) * yc - mu;
      ss += dot(u,u);
    }
    double shape = 0.5 * (T0*N) + a0;
    double rate  = 0.5 * ss + b0;
    s2 = rinvgamma(shape, rate);

    // ---------- (5) sample rho | rest (RW-MH) ----------
    auto ll_rho = [&](double r)->double{
      if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
      mat A = I_N - r * WW;
      double logdetA = logdet_signed(A);
      if (!is_finite(logdetA)) return -std::numeric_limits<double>::infinity();
      double ss2 = 0.0;
      for (int t=0; t<T0; ++t) {
        vec yc = Yc.row(t).t();
        vec mu = r * ww * Y0[t];
        if (useX) {
          mat Xt = X_get_row(t);
          mu += Xt * beta;
        }
        if (p>0) mu += Lambda * F.col(t);
        vec u = A * yc - mu;
        ss2 += dot(u,u);
      }
      double ll = T0 * logdetA - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss2 / s2;
      return ll;
    };

    double prop = R::rnorm(rho, step_rho);
    double lcur = ll_rho(rho);
    double lprp = ll_rho(prop);
    if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) { rho = prop; if (it>=burn) acc++; }

    // ---------- store after burn ----------
    if (it >= burn) {
      int m = it - burn;
      rho_draws[m] = rho;
      s2_draws[m]  = s2;
      if (K>0) beta_draws.row(m) = beta.t();
      if (p>0) {
        Lambda_draws.slice(m) = Lambda;
        F_draws.slice(m)      = F;
      }
    }
    if (verbose && (it % 2000 == 0)) Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    _["rho"]     = rho_draws,
    _["beta"]    = beta_draws,      // M x K (K==0なら空行列)
    _["sigma2"]  = s2_draws,
    _["Lambda"]  = Lambda_draws,    // N x p x M  (p==0なら空)
    _["F"]       = F_draws,         // p x T0 x M (p==0なら空)
    _["acc_rate"]= acc / std::max(1, M)
  );
}