// 40_geweke.cpp
// Geweke (2004) Joint Distribution Test helpers:
//  - Prior sampler for theta
//  - Forward simulator y | theta  (pre-treatment panel Yc given Y0, X, W, w)
// 
// This file is self-contained and does NOT depend on your posterior samplers.
// It exports the C++ pieces used by 40_geweke.R, which orchestrates the JDT.
//
// Compile from R via: Rcpp::sourceCpp("40_geweke.cpp")
//
// Notes:
// - Uses the same data layout for X as in your MCMC: flattened as ((t*N + i)*K + k)
// - RNG comes from R (Rcpp::RNGScope), so set.seed() in R controls randomness.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// =========================== Utilities ===========================
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1.0 / scale); }

inline double rhalfcauchy() { return std::abs(R::rcauchy(0.0, 1.0)); }

inline void symmetrize_inplace(mat& A) { A = 0.5 * (A + A.t()); }

inline bool robust_chol(mat& R, mat A, double base_eps=1e-12) {
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

inline vec X_get_row_from_flattened(const Rcpp::Nullable<Rcpp::NumericVector>& Xc_pre_,
                                    int t, int N, int K) {
  // returns a vectorized Xt * beta result is computed elsewhere;
  // Here we return a flattened row-major view: actually we build a matrix row
  vec out; // not used; kept for compatibility
  return out;
}

// Build the N x K regressor matrix Xt for time t from flattened Xc (length T0*N*K)
inline mat X_row_from_flat(const Rcpp::Nullable<Rcpp::NumericVector>& Xc_pre_,
                           int t, int N, int K) {
  mat Xt(N, K, fill::zeros);
  if (!Xc_pre_.isNotNull() || K==0) return Xt;
  Rcpp::NumericVector Xvec = Xc_pre_.get();
  for (int i=0;i<N;++i) {
    for (int k=0;k<K;++k) {
      int idx = (t * N + i) * K + k;
      Xt(i,k) = Xvec[idx];
    }
  }
  return Xt;
}

// ======================== 1) Prior sampler =======================
// Prior used here is chosen to be compatible with common choices in the posterior:
//   rho    ~ Uniform( -0.95/|lambda_max(W)| , +0.95/|lambda_max(W)| )
//   sigma2 ~ Inv-Gamma(a0, b0)
//   beta_j | tau, lambda_j ~ N(0, sigma2 * tau^2 * lambda_j^2),
//     tau ~ half-Cauchy(0,1), lambda_j ~ half-Cauchy(0,1)
//   Lambda_i (row i) ~ N(0, s2_eta * diag(1 / omega)), with
//     s2_eta ~ Inv-Gamma(1,1), omega_k ~ Inv-Gamma(1,1)  (k = 1..p)
//   F (Gamma) follows AR(1): F_t = phi*F_{t-1} + eps_t, eps_t ~ N(0, s2_g I_p),
//     phi ~ Uniform(-0.95, 0.95), s2_g ~ Inv-Gamma(1,1), F_0 = 0.
//
// Returned object: list(rho, sigma2, beta, Lambda, F)
//
// [[Rcpp::export]]
Rcpp::List sample_prior_theta_cpp(int T0, int N, int K, int p,
                                  const arma::mat& W,
                                  double a0 = 1.0, double b0 = 1.0) {
  Rcpp::RNGScope scope;

  // --- spectral bound for rho
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs = 0.0;
  for (arma::uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);
  double rho = R::runif(-bnd, bnd);

  // --- sigma^2
  double nu_sigma2 = rinvgamma(1.0, 1.0 / 100.0);   // nu_sigma2 ~ IG(1, 1/100)
  double s2        = rinvgamma(1.0, 1.0 / nu_sigma2); // sigma2 | nu ~ IG(1, scale=1/nu)

  // --- beta with Horseshoe prior
  arma::vec beta(K, fill::zeros);
  if (K > 0) {
    double tau = rhalfcauchy();
    double tau2 = tau * tau;
    for (int j=0;j<K;++j) {
      double lj = rhalfcauchy();
      double var_j = s2 * tau2 * lj * lj;
      beta(j) = std::sqrt( (var_j > 1e-12) ? var_j : 1e-12 ) * R::rnorm(0.0, 1.0);
    }
  }

  // --- Lambda (N x p), hierarchical normal with diagonal precision
  arma::mat Lambda(N, p, fill::zeros);
  if (p > 0) {
    double s2_eta = rinvgamma(1.0, 1.0);
    arma::vec omega(p, fill::ones);
    for (int k=0;k<p;++k) omega(k) = rinvgamma(1.0, 1.0);
    for (int i=0;i<N;++i) {
      for (int k=0;k<p;++k) {
        double var_ik = (s2_eta / std::max(omega(k), 1e-12));
        Lambda(i,k) = std::sqrt( (var_ik > 1e-12) ? var_ik : 1e-12 ) * R::rnorm(0.0, 1.0);
      }
    }
  }

  // --- F (p x T0): AR(1)
  arma::mat F(p, T0, fill::zeros);
  if (p > 0) {
    double phi = R::runif(-0.95, 0.95);
    double s2_g = rinvgamma(1.0, 1.0);
    // initialize F_1 ~ N(0, s2_g I)
    for (int k=0;k<p;++k) F(k,0) = std::sqrt(s2_g) * R::rnorm(0.0, 1.0);
    for (int t=1; t<T0; ++t) {
      for (int k=0;k<p;++k) {
        double mean = phi * F(k, t-1);
        F(k,t) = mean + std::sqrt(s2_g) * R::rnorm(0.0, 1.0);
      }
    }
  }

  return Rcpp::List::create(
    _["rho"]    = rho,
    _["sigma2"] = s2,
    _["beta"]   = beta,
    _["Lambda"] = Lambda,
    _["F"]      = F
  );
}

// ===================== 2) Forward simulator y|theta ===============
// Simulate pre-treatment controls panel Yc given parameters.
// Model per time t:
//   (I - rho W) * Yc_t = rho * w * Y0_t + X_t beta + Lambda f_t + eps_t
// with eps_t ~ N(0, sigma^2 I_N).
//
// Arguments:
//  - Y0_pre: length T0 vector
//  - Xc_pre_: flattened regressors of length T0*N*K (idx=(t*N+i)*K+k), or NULL if K=0
//  - T0, N, K, p: dimensions
//  - w: length N vector
//  - W: N x N matrix
//  - rho, sigma2, beta (K), Lambda (N x p), F (p x T0)
//
// Returns: Yc (T0 x N)
//
// [[Rcpp::export]]
arma::mat simulate_Yc_given_theta_cpp(const arma::vec& Y0_pre,
                                      const Rcpp::Nullable<Rcpp::NumericVector>& Xc_pre_,
                                      int T0, int N, int K, int p,
                                      const arma::vec& w,
                                      const arma::mat& W,
                                      double rho,
                                      double sigma2,
                                      const arma::vec& beta,
                                      const arma::mat& Lambda,
                                      const arma::mat& F) {
  Rcpp::RNGScope scope;

  arma::mat Yc(T0, N, fill::zeros);
  arma::mat I_N = arma::eye(N, N);
  arma::mat A = I_N - rho * W;

  // Factor term per t: Lambda * f_t, where F is p x T0 (columns are f_t)
  for (int t=0; t<T0; ++t) {
    arma::vec mu = rho * w * Y0_pre[t];                // N
    if (K > 0) {
      arma::mat Xt = X_row_from_flat(Xc_pre_, t, N, K);
      mu += Xt * beta;                                 // N
    }
    if (p > 0) {
      mu += Lambda * F.col(t);                         // N
    }
    arma::vec eps = std::sqrt( (sigma2 > 1e-12) ? sigma2 : 1e-12 ) * arma::randn<vec>(N);

    // Solve A * Yc_t = mu + eps
    // Solve A * Yc_t = mu + eps
    arma::vec rhs = mu + eps;

    // Primary: direct solve for general (possibly non-symmetric) A
    arma::vec Yt;
    bool ok = false;
    try {
      Yt = arma::solve(A, rhs, arma::solve_opts::fast); // LU-based
      ok = Yt.is_finite();
    } catch(...) {
      ok = false;
    }

    // Fallback: Tikhonov regularization (A' A + eps I) y = A' rhs
    if (!ok) {
      arma::mat AtA = A.t() * A;
      AtA.diag() += 1e-10;
      arma::vec Atrhs = A.t() * rhs;
      try {
        Yt = arma::solve(AtA, Atrhs, arma::solve_opts::fast);
        ok = Yt.is_finite();
      } catch(...) {
        ok = false;
      }
    }

    if (!ok) Rcpp::stop("simulate_Yc_given_theta_cpp: failed to solve (I - rho W) * y = rhs");
        Yc.row(t) = Yt.t();
  }

  return Yc;
}
