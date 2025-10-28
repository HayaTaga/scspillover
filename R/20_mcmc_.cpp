// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// =========================== Utilities ===========================
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1.0/scale); }

// stable log|det|
inline double logdet_stable(const arma::mat& A) {
  double sign = 0.0, val = 0.0;
  arma::log_det(val, sign, A);
  return val;
}

inline double clip1(double z, double lo = 1e-12, double hi = 1e12) {
  return std::min(std::max(z, lo), hi);
}
inline arma::vec clip_vec(const arma::vec& v, double lo = 1e-12, double hi = 1e12) {
  arma::vec out = v;
  for (arma::uword i = 0; i < out.n_elem; ++i) out(i) = clip1(out(i), lo, hi);
  return out;
}
inline void symmetrize_inplace(arma::mat& A) { A = 0.5 * (A + A.t()); }

inline bool robust_chol(arma::mat& R, arma::mat A, double base_eps=1e-10) {
  symmetrize_inplace(A);
  for (int k=0; k<6; ++k) {
    if (chol(R, A)) return true;
    A.diag() += base_eps * std::pow(10.0, k);
    symmetrize_inplace(A);
  }
  arma::vec d; arma::mat V;
  if (eig_sym(d, V, A)) {
    for (uword i=0;i<d.n_elem;++i) if (d(i) < base_eps) d(i) = base_eps;
    A = V * diagmat(d) * V.t();
    symmetrize_inplace(A);
  }
  return chol(R, A);
}

// ================================================================
// Step 1: BSCM  (y0 = Yc * alpha + eps), Horseshoe prior on alpha
// Returns M x N matrix of alpha draws (post-burn)
// ================================================================
// [[Rcpp::export]]
arma::mat hs_alpha_gibbs_cpp(const arma::vec& Y0_pre,
                             const arma::mat& control_outcome_pre,
                             int iteration,
                             int burn,
                             double a0 = 1.0,     // IG prior for sigma2: shape a0
                             double b0 = 1.0,     // IG prior for sigma2: scale b0
                             bool verbose = false) {
  Rcpp::RNGScope scope;

  const int T0 = control_outcome_pre.n_rows;
  const int N  = control_outcome_pre.n_cols;
  const int M  = std::max(0, iteration - burn);

  arma::mat X = control_outcome_pre; // T0 x N
  arma::vec y = Y0_pre;              // T0

  // States
  arma::vec alpha = arma::zeros<arma::vec>(N);
  double sigma2   = 1.0;

  // Horseshoe auxiliaries (Makalic–Schmidt)
  arma::vec lambda2    = arma::ones<arma::vec>(N); // local scales
  arma::vec nu_lambda  = arma::ones<arma::vec>(N);
  double tau2          = 1.0;                      // global scale
  double nu_tau        = 1.0;

  arma::mat XtX = X.t() * X;                       // N x N
  arma::vec Xty = X.t() * y;                       // N

  arma::mat out(M, N, arma::fill::zeros);

  for (int it = 0; it < iteration; ++it) {
    // (1) alpha | sigma2, tau2, lambda2, y
    arma::vec inv_prior = 1.0 / clip_vec(tau2 * lambda2); // diag prior precision (since prior Var = tau2*lambda2)
    arma::mat Prec = XtX / clip1(sigma2) + diagmat(inv_prior);
    symmetrize_inplace(Prec);

    arma::mat L;
    if (!robust_chol(L, Prec)) Rcpp::stop("chol failed in alpha update");
    arma::vec m = solve(Prec, Xty / clip1(sigma2)); // mean
    // draw via solving L L' z = ...
    arma::vec z = arma::randn<arma::vec>(N);
    // Use precision factorization: alpha = m + L^{-T} z
    alpha = m + solve(trimatu(L.t()), z);

    // (2) sigma2 | rest
    arma::vec resid = y - X * alpha;
    double ss = arma::dot(resid, resid);
    sigma2 = rinvgamma(a0 + 0.5 * T0, b0 + 0.5 * ss);

    // (3) local scales lambda2_j
    for (int j=0; j<N; ++j) {
      double rate = 1.0/clip1(nu_lambda(j)) + 0.5 * (alpha(j)*alpha(j)) / clip1(sigma2 * tau2);
      lambda2(j) = rinvgamma(1.0, rate);
      double rate_nu = 1.0 + 1.0/clip1(lambda2(j));
      nu_lambda(j) = rinvgamma(1.0, rate_nu);
    }

    // (4) global scale tau2
    double sum_term = 0.0;
    for (int j=0; j<N; ++j) sum_term += (alpha(j)*alpha(j)) / clip1(sigma2 * lambda2(j));
    double rate_tau = 1.0/clip1(nu_tau) + 0.5 * sum_term;
    tau2 = rinvgamma(0.5*(N + 1.0), rate_tau);
    nu_tau = rinvgamma(1.0, 1.0 + 1.0/clip1(tau2));

    // store
    if (it >= burn) {
      out.row(it - burn) = alpha.t();
      if (verbose && (((it - burn + 1) % 2000) == 0)) Rcpp::checkUserInterrupt();
    }
  }

  return out; // M x N
}

// ================================================================
// Step 2: SAR with fixed alpha_hat
//  (I - rho W - rho w alpha^T) Yc_t = X_t beta + Lambda F_t + u_t
//  alpha is fixed at alpha_hat (from Step 1 posterior mean)
//  p: number of latent factors (if 0, no factors)
// ================================================================

// [[Rcpp::export]]
Rcpp::List sar_full_sampler_cpp_step2(const arma::mat& Yc_pre,           // T0 x N
                                      const arma::vec& alpha_hat_in,     // N (unscaled, matches Yc_pre)
                                      Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K)
                                      int T0, int N, int K, int p,
                                      const arma::vec& w_in,             // N
                                      const arma::mat& W,                // N x N
                                      int iteration, int burn,
                                      double step_rho = 0.01,
                                      double a0 = 1.0, double b0 = 1.0,
                                      bool verbose=false) {
  Rcpp::RNGScope scope;

  const int M = std::max(0, iteration - burn);

  // ---------- scale Yc columns by SD (improves numeric stability) ----------
  arma::mat Yc_orig = Yc_pre; // T0 x N
  arma::vec sds_Yc(N);
  arma::mat Yc = Yc_orig;
  for (int j = 0; j < N; ++j) {
    double sdj = arma::stddev(Yc_orig.col(j));
    if (!arma::is_finite(sdj) || sdj < 1e-8) sdj = 1.0;
    sds_Yc(j) = sdj;
    Yc.col(j) = Yc_orig.col(j) / sdj;
  }

  // alpha used internally must be consistent with scaled Yc:
  // alpha_int^T * (Yc_orig / sds) = (alpha_unscaled^T / sds) Yc_orig
  arma::vec alpha = alpha_hat_in / sds_Yc; // N

  // normalize w (reduce confounding scale with rho)
  arma::vec w = w_in;
  double wnorm = std::sqrt(arma::dot(w, w));
  if (wnorm > 0.0) w /= wnorm;

  // X accessor
  const bool useX = (K > 0) && Xc_pre_.isNotNull();
  Rcpp::NumericVector Xvec;
  if (useX) Xvec = Xc_pre_.get();
  auto X_get_row = [&](int t)->arma::mat {
    arma::mat Xt(N, K, arma::fill::zeros);
    if (!useX) return Xt;
    for (int i=0;i<N;++i)
      for (int k=0;k<K;++k)
        Xt(i,k) = Xvec[(t * N + i) * K + k];
    return Xt;
  };

  // spectral bound for rho (stability)
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs = 0.0;
  for (uword i=0; i<evals.n_elem; ++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // storage
  arma::vec rho_draws(M, arma::fill::none);
  arma::vec s2_draws(M, arma::fill::none);
  arma::mat beta_draws(M, K, arma::fill::zeros);
  arma::cube Lambda_draws(N, p, M, arma::fill::zeros); // Eta
  arma::cube F_draws(p, T0, M, arma::fill::zeros);     // Gamma

  // states
  double rho = 0.0;
  double s2  = 1.0;

  arma::vec beta = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());

  arma::mat Eta   = (p>0 ? arma::zeros<arma::mat>(N,p) : arma::mat()); // N x p
  arma::mat Gamma = (p>0 ? arma::zeros<arma::mat>(p,T0) : arma::mat()); // p x T0
  double phi_g = 0.0, s2_g = 1.0, nu_s2_g = 1.0;
  arma::vec omega_k    = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  arma::vec nu_omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  double s2_eta = 1.0, nu_s2_eta = 1.0;

  double nu_sigma2 = 1.0;

  arma::mat I_N = arma::eye(N,N);
  arma::mat I_K = (K>0 ? arma::eye(K,K) : arma::mat());
  arma::mat I_p = (p>0 ? arma::eye(p,p) : arma::mat());

  auto Mtilde = [&](double r)->arma::mat {
    // Ã(r) = I − r W − r w αᵀ  (alpha is fixed)
    return I_N - r * W - r * w * alpha.t();
  };

  auto loglik_core_pair = [&](double r)->std::pair<double,double> {
    if (std::abs(r) >= bnd) return {-std::numeric_limits<double>::infinity(), 0.0};
    arma::mat M = Mtilde(r);
    double ldetM = logdet_stable(M);
    if (!std::isfinite(ldetM)) return {-std::numeric_limits<double>::infinity(), 0.0};

    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = M * Yc.row(t).t();
      if (useX)  u -= X_get_row(t) * beta;
      if (p>0)   u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    double ll = T0 * ldetM - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
    return {ll, ss};
  };

  // RWMH adaptation for rho
  int acc_rho = 0;
  double log_step_rho = std::log(step_rho);
  const double target_accept_rho = 0.234;
  const double adapt_gamma = 0.6;

  for (int it=0; it<iteration; ++it) {

    // (1) Gamma | rest (FFBS) when p>0
    if (p > 0) {
      arma::mat Phi   = I_p * phi_g;
      arma::mat Q_mat = I_p * s2_g;
      arma::mat H_mat = Eta;

      arma::mat HtH_s2_inv = (H_mat.t() * H_mat) / clip1(s2);
      arma::mat Ht_s2_inv  = H_mat.t() / clip1(s2);

      std::vector<arma::vec> gamma_t_t(T0);
      std::vector<arma::mat> P_t_t(T0);
      std::vector<arma::vec> gamma_t_t1(T0);
      std::vector<arma::mat> P_t_t1(T0);

      arma::vec gamma_prev = arma::zeros<arma::vec>(p);
      arma::mat P_prev = (clip1(s2_g) / std::max(1e-6, 1.0 - phi_g*phi_g)) * I_p;

      for (int t = 0; t < T0; ++t) {
        arma::mat M = Mtilde(rho);
        arma::vec Y_star_t = M * Yc.row(t).t();
        if (useX) Y_star_t -= X_get_row(t) * beta;

        arma::vec pred = Phi * gamma_prev;
        arma::mat P_pred = Phi * P_prev * Phi.t() + Q_mat;

        arma::mat P_inv_pred = arma::inv_sympd(P_pred);
        arma::mat V_inv_t = P_inv_pred + HtH_s2_inv;
        arma::mat V_t = arma::inv_sympd(V_inv_t);

        arma::vec m_t = V_t * (P_inv_pred * pred + Ht_s2_inv * Y_star_t);

        gamma_t_t[t]  = m_t;
        P_t_t[t]      = V_t;
        gamma_t_t1[t] = pred;
        P_t_t1[t]     = P_pred;

        gamma_prev = m_t;
        P_prev     = V_t;
      }

      arma::vec m_T = gamma_t_t[T0 - 1];
      arma::mat V_T = P_t_t[T0 - 1];
      arma::mat L_T = arma::chol(0.5 * (V_T + V_T.t()), "lower");
      Gamma.col(T0 - 1) = m_T + L_T * arma::randn<arma::vec>(p);

      for (int t = T0 - 2; t >= 0; --t) {
        arma::vec g_next = Gamma.col(t + 1);
        arma::mat P_inv_next_pred = arma::inv_sympd(P_t_t1[t + 1]);
        arma::mat J_t = P_t_t[t] * Phi.t() * P_inv_next_pred;

        arma::vec m_s = gamma_t_t[t] + J_t * (g_next - gamma_t_t1[t + 1]);
        arma::mat V_s = P_t_t[t] - J_t * Phi * P_t_t[t];

        arma::mat Ls = arma::chol(0.5 * (V_s + V_s.t()), "lower");
        Gamma.col(t) = m_s + Ls * arma::randn<arma::vec>(p);
      }

      // phi_g
      double den = 0.0, num = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        den += arma::dot(gl, gl);
        num += arma::dot(gl, Gamma.col(t));
      }
      double mean_phi = (den > 0 ? num / den : 0.0);
      double var_phi  = (den > 0 ? s2_g / den : 1.0);
      double cand_phi;
      do { cand_phi = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand_phi) > 1.0);
      phi_g = cand_phi;

      // s2_g
      double sc_g = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        arma::vec diff = Gamma.col(t) - phi_g * gl;
        sc_g += 0.5 * arma::dot(diff, diff);
      }
      s2_g     = rinvgamma(0.5 + 0.5 * p * T0, sc_g + 1.0/clip1(nu_s2_g));
      nu_s2_g  = rinvgamma(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
    }

    // (2) Eta | rest
    if (p>0) {
      arma::mat GtG = Gamma * Gamma.t();     // p x p
      arma::mat Domega = diagmat(omega_k);   // p x p

      arma::mat Vrow = arma::inv_sympd( GtG / s2 + Domega / clip1(s2_eta) );
      arma::mat Lrow = chol(Vrow, "lower");

      for (int i=0; i<N; ++i) {
        arma::vec rhs = arma::zeros<arma::vec>(p);
        for (int t=0; t<T0; ++t) {
          arma::vec r = Mtilde(rho) * Yc.row(t).t();
          if (useX) r -= X_get_row(t) * beta;
          rhs += Gamma.col(t) * r(i);
        }
        arma::vec m = Vrow * (rhs / s2);
        arma::vec z = arma::randn<arma::vec>(p);
        Eta.row(i) = (m + Lrow * z).t();
      }

      double sc_eta = 0.0;
      for (int i=0;i<N;++i) {
        arma::vec ei = Eta.row(i).t();
        sc_eta += arma::dot(ei, Domega * ei);
      }
      s2_eta    = rinvgamma(0.5 + 0.5 * p * N, 0.5*sc_eta + 1.0/clip1(nu_s2_eta));
      nu_s2_eta = rinvgamma(1.0, 1.0/clip1(s2_eta) + 1.0/100.0);

      for (int k_idx=0;k_idx<p;++k_idx) {
        double tmp = 0.0;
        for (int i=0;i<N;++i) tmp += 0.5 * Eta(i,k_idx)*Eta(i,k_idx) / clip1(s2_eta);
        double rate_ok = 1.0/clip1(nu_omega_k(k_idx)) + tmp;
        omega_k(k_idx)     = rinvgamma(0.5*(N+1.0), rate_ok);
        nu_omega_k(k_idx)  = rinvgamma(1.0, 1.0 + 1.0/clip1(omega_k(k_idx)));
      }
    }

    // (3) beta | rest (Gaussian, with weak ridge for stability)
    if (useX) {
      arma::mat Ab = arma::zeros<arma::mat>(K,K);
      arma::vec Bb = arma::zeros<arma::vec>(K);

      for (int t=0; t<T0; ++t) {
        arma::mat Xt = X_get_row(t);          // N x K
        Ab += Xt.t() * Xt;
        arma::vec Btmp = Mtilde(rho) * Yc.row(t).t();
        if (p>0) Btmp -= Eta * Gamma.col(t);
        Bb += Xt.t() * Btmp;
      }
      // ridge with tiny prior variance ~ s2 * 1e6
      Ab += (s2 * 1e-6) * I_K;

      arma::mat Ainv = arma::inv_sympd(Ab);
      arma::vec m_b  = Ainv * Bb;
      arma::mat S_b  = clip1(s2) * Ainv;

      beta = arma::mvnrnd(m_b, 0.5*(S_b+S_b.t()), 1);
    }

    // (4) sigma^2 | rest
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = Mtilde(rho) * Yc.row(t).t();
      if (useX) u -= X_get_row(t) * beta;
      if (p>0)  u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    s2 = rinvgamma(a0 + 0.5 * (T0 * N), b0 + 0.5 * ss);
    nu_sigma2 = rinvgamma(1.0, 1.0/clip1(s2) + 1.0/100.0); // weakly inform

    // (5) rho | rest (Adaptive RWMH)
    {
      double current_step = std::exp(log_step_rho);
      double prop_rho = R::rnorm(rho, current_step);

      double lcur = loglik_core_pair(rho).first;
      double lprp = loglik_core_pair(prop_rho).first;

      bool accepted = false;
      double loga = lprp - lcur;
      if (std::log(R::runif(0.0,1.0)) < loga) {
        rho = prop_rho;
        accepted = true;
        if (it >= burn) acc_rho++;
      }
      double adapt_step = std::pow(it + 1.0, -adapt_gamma);
      log_step_rho += adapt_step * ((accepted ? 1.0 : 0.0) - target_accept_rho);
      log_step_rho = std::max(-10.0, std::min(log_step_rho, 3.0));
    }

    // store
    if (it >= burn) {
      int m = it - burn;
      rho_draws[m] = rho;
      s2_draws[m]  = s2;
      if (K>0) beta_draws.row(m) = beta.t();
      if (p>0) { Lambda_draws.slice(m) = Eta; F_draws.slice(m) = Gamma; }
      if (verbose && (((m+1) % 2000) == 0)) Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(
    _["rho"]        = rho_draws,
    _["sigma2"]     = s2_draws,
    _["beta"]       = beta_draws,        // M x K (if K=0: empty)
    _["Lambda"]     = Lambda_draws,      // N x p x M (if p=0: empty)
    _["F"]          = F_draws,           // p x T0 x M (if p=0: empty)
    _["acc_rho"]    = acc_rho / std::max(1, M),
    _["final_log_step_rho"] = log_step_rho,
    _["sds_Yc"]     = sds_Yc            // 参考：内部で使用したスケール
  );
}