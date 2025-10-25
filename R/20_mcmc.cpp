// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// =========================== Utilities ===========================
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1 / scale); }

inline double logdet(const arma::mat& A) {
  double sign=0.0, val=0.0;
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

inline bool robust_chol(arma::mat& R, arma::mat& A, double base_eps=1e-12) {
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
//  Full SAR sampler with latent factors and sparse alpha (HS prior)
//  alpha update: blocked Independent MH using Gaussian approximation
// ================================================================

// [[Rcpp::export]]
Rcpp::List sar_full_sampler_cpp(const arma::vec& Y0_pre,    // 未使用
                                const arma::mat& Yc_pre,    // T0 x N
                                Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K)
                                int T0, int N, int K, int p,
                                const arma::vec& w_in,      // N
                                const arma::mat& W,         // N x N
                                int iteration, int burn,
                                double step_rho = 0.01,     // rho用の初期RW幅
                                double step_alpha = 0.01,   // 互換性のため残すが未使用
                                double a0 = 1.0, double b0 = 1.0,
                                bool verbose=false) {

  Rcpp::RNGScope scope;

  const int M = std::max(0, iteration - burn);

  // ----------------------------
  // data preprocessing / scaling
  // ----------------------------
  arma::mat Yc_orig = Yc_pre; // store original scale
  arma::vec sds_Yc(N);
  arma::mat Yc = Yc_orig;
  for (int j = 0; j < N; ++j) {
    double sdj = arma::stddev(Yc_orig.col(j));
    if (!arma::is_finite(sdj) || sdj < 1e-8) sdj = 1.0;
    sds_Yc(j) = sdj;
    Yc.col(j) = Yc_orig.col(j) / sdj;
  }

  // X accessor
  const bool useX = (K > 0) && Xc_pre_.isNotNull();
  Rcpp::NumericVector Xvec;
  if (useX) Xvec = Xc_pre_.get();
  auto X_get_row = [&](int t)->arma::mat {
    arma::mat Xt(N, K, arma::fill::zeros);
    if (!useX) return Xt;
    for (int i=0;i<N;++i){
      for (int k=0;k<K;++k){
        Xt(i,k) = Xvec[(t * N + i) * K + k];
      }
    }
    return Xt;
  };

  // ----------------------------------------
  // spectral bound for rho (stability region)
  // ----------------------------------------
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs = 0.0;
  for (arma::uword i=0;i<evals.n_elem;++i)
    maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // -----------------
  // storage for draws
  // -----------------
  arma::vec rho_draws(M, arma::fill::none);
  arma::vec s2_draws(M, arma::fill::none);
  arma::mat beta0_draws(M, K, arma::fill::zeros);
  arma::cube Lambda_draws(N, p, M, arma::fill::zeros); // Eta
  arma::cube F_draws(p, T0, M, arma::fill::zeros);     // Gamma
  arma::mat alpha_draws(M, N, arma::fill::zeros);

  // ---------------
  // initial states
  // ---------------
  double rho = 0.0;
  double s2  = 1.0;

  arma::vec beta0      = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());
  arma::vec sig2_b0    = (K>0 ? arma::ones<arma::vec>(K)  : arma::vec()); // λ_j^2
  arma::vec nu_sig_b0  = (K>0 ? arma::ones<arma::vec>(K)  : arma::vec()); // ν_j
  double tau2_b0 = 1.0, nu_tau_b0 = 1.0;

  double nu_sigma2 = 1.0;

  arma::mat Eta   = (p>0 ? arma::zeros<arma::mat>(N,p) : arma::mat()); // N x p
  arma::mat Gamma = (p>0 ? arma::zeros<arma::mat>(p,T0) : arma::mat()); // p x T0
  double phi_g = 0.0, s2_g = 1.0, nu_s2_g = 1.0;
  arma::vec omega_k    = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  arma::vec nu_omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  double s2_eta = 1.0, nu_s2_eta = 1.0;

  arma::vec alpha = 1e-4 * arma::randn<arma::vec>(N); // sparse weights
  arma::vec sigma2_i   = arma::ones<arma::vec>(N);    // local HS scales λ_i^2
  arma::vec nu_sigma_i = arma::ones<arma::vec>(N);    // ν_{λ_i}
  double tau2 = 1.0, nu_tau = 1.0;                    // global HS scale τ^2

  // normalize w once (to reduce scaling confounding between rho and alpha)
  arma::vec w = w_in;
  double wnorm = std::sqrt(arma::dot(w,w));
  if (wnorm > 0.0) w /= wnorm;

  arma::mat I_N = arma::eye(N,N);
  arma::mat I_K = (K>0 ? arma::eye(K,K) : arma::mat());
  arma::mat I_p = (p>0 ? arma::eye(p,p) : arma::mat());

  // helper: Ã(r, α) = I − r W − r w αᵀ
  auto Mtilde = [&](double r, const arma::vec& a)->arma::mat {
    return I_N - r * W - r * w * a.t();
  };

  // log-likelihood core (includes log|det| term and Gaussian part)
  // returns pair( loglik , sumsq )
  auto loglik_core_pair = [&](double r, const arma::vec& a)->std::pair<double,double> {
    if (std::abs(r) >= bnd) return {-std::numeric_limits<double>::infinity(), 0.0};
    arma::mat M = Mtilde(r, a);
    double ldetM = logdet(M);
    if (!std::isfinite(ldetM)) return {-std::numeric_limits<double>::infinity(), 0.0};

    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = M * Yc.row(t).t();
      if (useX)  u -= X_get_row(t) * beta0;
      if (p>0)   u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    double ll = T0 * ldetM - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
    return {ll, ss};
  };

  // full log posterior for alpha up to additive constant:
  //   log p(alpha | rest) ∝ loglik_core(rho, alpha)  +  log N(alpha | 0, diag(sigma2_i))
  auto logpost_alpha = [&](const arma::vec& a)->double {
    auto pr = loglik_core_pair(rho, a);
    if (!std::isfinite(pr.first)) return -std::numeric_limits<double>::infinity();
    double ll = pr.first;
    double lp = 0.0;
    for (int i=0;i<N;++i) {
      double s2i = clip1(sigma2_i(i));
      lp += -0.5 * std::log(s2i) - 0.5 * (a(i)*a(i))/s2i;
    }
    return ll + lp;
  };

  int iters = iteration;
  int acc_rho = 0;
  int acc_alpha = 0; // acceptance count for alpha block MH

  // rho step size adaptation
  double log_step_rho = std::log(step_rho);
  double target_accept_rho = 0.234;
  double adapt_gamma = 0.6;

  // 係数（alpha提案の温度スケール）
  const double c_alpha_scale = 0.5; // <1で保守的に

  // 定数
  const double log2pi = std::log(2.0 * 3.14159265358979323846);

  arma::vec log_step_alpha(N);                     // N次元ベクトルとして宣言
  log_step_alpha.fill(std::log(step_alpha));       // 初期ステップ幅の対数を全要素に設定
  double target_accept_alpha = 0.44;              // 目標受理率
  // double adapt_gamma = 0.6;                     // rho と共通の適応減衰率を使う
  arma::ivec acc_alpha_counts = arma::zeros<arma::ivec>(N);

  for (int it=0; it<iters; ++it) {

    // ============================================================
    // (1) Gamma (factor scores) | rest  via FFBS
    // ============================================================
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
        arma::mat M = Mtilde(rho, alpha);
        arma::vec Y_star_t = M * Yc.row(t).t();
        if (useX) Y_star_t -= X_get_row(t) * beta0;

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

      // phi_g | rest (AR(1) coeff)
      double den = 0.0, num = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        den += arma::dot(gl, gl);
        num += arma::dot(gl, Gamma.col(t));
      }
      double mean_phi = (den > 0 ? num / den : 0.0);
      double var_phi  = (den > 0 ? s2_g / den : 1.0);
      double cand_phi;
      do {
        cand_phi = R::rnorm(mean_phi, std::sqrt(var_phi));
      } while (std::abs(cand_phi) > 1.0);
      phi_g = cand_phi;

      // s2_g | rest
      double sc_g = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        arma::vec diff = Gamma.col(t) - phi_g * gl;
        sc_g += 0.5 * arma::dot(diff, diff);
      }
      s2_g     = rinvgamma(0.5 + 0.5 * p * T0, sc_g + 1.0/clip1(nu_s2_g));
      nu_s2_g  = rinvgamma(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
    }

    // ============================================================
    // (2) Eta (factor loadings) | rest
    // ============================================================
    if (p>0) {
      arma::mat GtG = Gamma * Gamma.t();     // p x p
      arma::mat Domega = arma::diagmat(omega_k); // p x p

      arma::mat Vrow = arma::inv_sympd( GtG / s2 + Domega / clip1(s2_eta) );
      arma::mat Lrow = arma::chol(Vrow, "lower");

      for (int i=0; i<N; ++i) {
        arma::vec rhs = arma::zeros<arma::vec>(p);
        for (int t=0; t<T0; ++t) {
          arma::vec r = Mtilde(rho, alpha) * Yc.row(t).t();
          if (useX) r -= X_get_row(t) * beta0;
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
        for (int i=0;i<N;++i) {
          tmp += 0.5 * Eta(i,k_idx)*Eta(i,k_idx) / clip1(s2_eta);
        }
        double rate_ok = 1.0/clip1(nu_omega_k(k_idx)) + tmp;
        omega_k(k_idx)     = rinvgamma(0.5*(N+1.0), rate_ok);
        nu_omega_k(k_idx)  = rinvgamma(1.0, 1.0 + 1.0/clip1(omega_k(k_idx)));
      }
    }

    // ============================================================
    // (3) beta0 | rest (Horseshoe prior)
    // ============================================================
    if (useX) {
      arma::mat Ab = arma::zeros<arma::mat>(K,K);
      arma::vec Bb = arma::zeros<arma::vec>(K);

      for (int t=0; t<T0; ++t) {
        arma::mat Xt = X_get_row(t);          // N x K
        Ab += Xt.t() * Xt;

        arma::vec Btmp = Mtilde(rho, alpha) * Yc.row(t).t();
        if (p>0) Btmp -= Eta * Gamma.col(t);
        Bb += Xt.t() * Btmp;
      }

      Ab.diag() += clip1(s2) * (1.0 / clip_vec(sig2_b0));
      arma::mat Ainv = arma::inv_sympd(Ab);
      arma::vec m_b  = Ainv * Bb;
      arma::mat S_b  = clip1(s2) * Ainv;

      beta0 = arma::mvnrnd(m_b, 0.5*(S_b+S_b.t()), 1);

      // Horseshoe local scales
      for (int j=0;j<K;++j) {
        double rate_l = 0.5 * beta0(j)*beta0(j) + 1.0/clip1(nu_sig_b0(j));
        sig2_b0(j)    = rinvgamma(1.0, rate_l);

        double rate_nu= 1.0/clip1(sig2_b0(j)) + 1.0/clip1(tau2_b0);
        nu_sig_b0(j)  = rinvgamma(1.0, rate_nu);
      }

      // Horseshoe global scale
      double sum_inv_nu = 0.0;
      for (int j=0;j<K;++j) sum_inv_nu += 1.0/clip1(nu_sig_b0(j));
      tau2_b0    = rinvgamma(1.0, 1.0/clip1(nu_tau_b0) + sum_inv_nu);
      nu_tau_b0  = rinvgamma(1.0, 1.0/clip1(tau2_b0) + 1.0/clip1(s2));
    }

    // ============================================================
    // (4) sigma^2 | rest
    // ============================================================
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = Mtilde(rho, alpha) * Yc.row(t).t();
      if (useX) u -= X_get_row(t) * beta0;
      if (p>0)  u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    double shape = a0 + 0.5 * (T0 * N);
    double rate  = 1.0/clip1(nu_sigma2) + 0.5 * ss;
    s2 = rinvgamma(shape, rate);
    nu_sigma2 = rinvgamma(1.0, 1.0/clip1(s2) + 1.0/100.0);

    // ============================================================
    // (5) rho | rest  (Adaptive RW–MH)
    // ============================================================
    double current_step_rho = std::exp(log_step_rho);
    double prop_rho = R::rnorm(rho, current_step_rho);

    auto lpair_cur = loglik_core_pair(rho, alpha);
    double lcur_rho = lpair_cur.first;
    double lprp_rho = loglik_core_pair(prop_rho, alpha).first;

    bool accepted_rho = false;
    double log_alpha_accept_rho = lprp_rho - lcur_rho;

    if (std::log(R::runif(0.0, 1.0)) < log_alpha_accept_rho) {
      rho = prop_rho;
      accepted_rho = true;
      if (it >= burn) acc_rho++;
    }

    // adapt rho step size
    {
      double adapt_step = std::pow(it + 1.0, -adapt_gamma);
      log_step_rho += adapt_step * ( (accepted_rho ? 1.0 : 0.0) - target_accept_rho );
    }

    // ============================================================
    // (6) alpha | rest  — Blocked Independent MH
    // ============================================================
    int current_acc_alpha_iter = 0;
    for (int i = 0; i < N; ++i) {
      double lcur_a = logpost_alpha(alpha); // 現在の対数事後

      arma::vec prop_alpha = alpha;
      double current_step_i = std::exp(log_step_alpha(i)); // i番目のステップ幅
      prop_alpha(i) = alpha(i) + R::rnorm(0.0, current_step_i);

      double lprp_a = logpost_alpha(prop_alpha); // 提案点の対数事後

      bool accepted_i = false;
      double log_alpha_accept_i = lprp_a - lcur_a;
      if (std::log(R::runif(0.0, 1.0)) < log_alpha_accept_i) {
        alpha = prop_alpha; // ★★★ 受理されたら alpha ベクトル全体を更新 ★★★
        accepted_i = true;
        current_acc_alpha_iter++;
      }

      // --- i番目のステップ幅を適応 ---
      double adapt_step_i = std::pow(it + 1.0, -adapt_gamma);
      log_step_alpha(i) += adapt_step_i * ( (accepted_i ? 1.0 : 0.0) - target_accept_alpha );
      // log_step_alpha(i) の範囲制限 (オプションだが推奨)
      log_step_alpha(i) = std::max(-10.0, std::min(log_step_alpha(i), 3.0));
    }
    if (it >= burn) acc_alpha += current_acc_alpha_iter;
//     {
//       // 準備:
//       //   Prec_q = (rho^2 * ||w||^2 / s2) * Σ_t (y_t y_t') + Diag(1/σ_i^2)
//       //   b_q    = -(rho / s2) * Σ_t [ y_t * (w^T R_t) ]
//       // where R_t = (I - rho W) y_t - X_t β0 - Eta γ_t
//       arma::vec inv_sigma2_i = 1.0 / clip_vec(sigma2_i);
//       arma::mat Syy(N, N, arma::fill::zeros);
//       arma::vec bq(N, arma::fill::zeros);

//       arma::mat I_N_local = I_N; // just alias for clarity

//       for (int t=0; t<T0; ++t) {
//         arma::vec y_t = Yc.row(t).t();
//         arma::vec R_t = (I_N_local - rho * W) * y_t;
//         if (useX) R_t -= X_get_row(t) * beta0;
//         if (p>0)  R_t -= Eta * Gamma.col(t);

//         Syy += y_t * y_t.t();
//         double wtR = arma::dot(w, R_t);
//         bq   -= (rho / clip1(s2)) * (y_t * wtR);
//       }

//       // precision matrix of Gaussian approx
//       arma::mat Prec = ((rho*rho) / clip1(s2)) * Syy; // note w already normalized
//       // multiply by ||w||^2 = 1 if w normalized; if not, incorporate wnorm^2:
//       // but we already normalized w -> ||w||=1
//       Prec.diag() += inv_sigma2_i;
//       // ridge for numerical stability
//       Prec.diag() += 1e-2;

//       // symmetrize then invert
//       arma::mat Prec_sym = 0.5*(Prec + Prec.t());
//       arma::mat Rchol_prec;
//       if (!robust_chol(Rchol_prec, Prec_sym)) {
//         Rcpp::stop("chol(Prec_sym) failed in alpha block MH.");
//       }
//       arma::mat Sig_q = arma::inv_sympd(Prec_sym); // covariance of base proposal (before scaling)

//       arma::vec mu_q = Sig_q * bq;

//       // scaled proposal covariance: c_alpha_scale * Sig_q
//       arma::mat Sig_prop = c_alpha_scale * 0.5*(Sig_q + Sig_q.t());
//       arma::mat L_prop;
//       {
//         arma::mat Sig_prop_sym = 0.5*(Sig_prop + Sig_prop.t());
//         if (!robust_chol(L_prop, Sig_prop_sym)) {
//           Rcpp::stop("chol(Sig_prop) failed in alpha block MH.");
//         }
//       }

//       // draw proposal alpha_prop ~ N(mu_q, Sig_prop)
//       arma::vec z_draw = arma::randn<arma::vec>(N);
//       arma::vec alpha_prop = mu_q + L_prop * z_draw;

//       // 対数事後 (target) の差
//       double lp_cur = logpost_alpha(alpha);
//       double lp_prp = logpost_alpha(alpha_prop);

//       // 提案密度 q の対数
//       // q(x) = N(x | mu_q, c_alpha_scale * Sig_q)
//       // Prec_prop = (1/c_alpha_scale) * Prec_sym
//       arma::mat Prec_prop = (1.0 / c_alpha_scale) * Prec_sym;

//       auto log_q_eval = [&](const arma::vec& x)->double {
//         arma::vec diff = x - mu_q;
//         double quad = arma::as_scalar(diff.t() * Prec_prop * diff); // (x-mu)' Prec_prop (x-mu)
//         double ldPrec = logdet(Prec_sym);           // log|Prec_sym|
//         double ldSig_q = - ldPrec;                  // log|Sig_q|
//         double ldSig_prop = (double)N * std::log(c_alpha_scale) + ldSig_q; // log|c * Sig_q|
//         double log_q = -0.5 * ( quad
//                                 + (double)N * log2pi
//                                 + ldSig_prop );
//         return log_q;
//       };

//       double log_q_cur = log_q_eval(alpha);
//       double log_q_prp = log_q_eval(alpha_prop);

//       // MH acceptance
//       double log_acc_alpha = (lp_prp - lp_cur) + (log_q_cur - log_q_prp);

//         // Rcpp::Rcout << "--- Iteration " << it << " ---" << std::endl;
//         // Rcpp::Rcout << "lp_cur: " << lp_cur << ", lp_prp: " << lp_prp << std::endl;
//         // Rcpp::Rcout << "log_q_cur: " << log_q_cur << ", log_q_prp: " << log_q_prp << std::endl;
//         // Rcpp::Rcout << "log_acc_alpha: " << log_acc_alpha << std::endl;
//     double log_u = std::log(R::runif(0.0,1.0)); // Store the random draw
// // Rcpp::Rcout << "Iter " << it << ": log_u=" << log_u << ", log_acc=" << log_acc_alpha << std::endl;


//       bool accepted_alpha = false;
//       if (log_u < log_acc_alpha) {
//         Rcpp::Rcout << "  ACCEPTED! Updating alpha." << std::endl;
//            Rcpp::Rcout << "    alpha BEFORE update (first elem): " << alpha(0) << std::endl;
//            Rcpp::Rcout << "    alpha_prop      (first elem): " << alpha_prop(0) << std::endl;
//         alpha = alpha_prop;
//         accepted_alpha = true;
//       }
//       if (accepted_alpha && (it >= burn)) {
//         acc_alpha++;
//       }
//     } // end alpha block

    // ============================================================
    // (7) Horseshoe hyperparams for alpha
    // ============================================================
    // sigma2_i[i] | alpha, nu_sigma_i[i]
    for (int i=0;i<N;++i) {
      double sc = 0.5 * alpha(i)*alpha(i) + 1.0 / clip1(nu_sigma_i(i));
      sigma2_i(i) = rinvgamma(1.0, sc);
    }
    // nu_sigma_i[i] | sigma2_i[i], tau2
    for (int i=0;i<N;++i) {
      double sc = 1.0 / clip1(sigma2_i(i)) + 1.0 / clip1(tau2);
      nu_sigma_i(i) = rinvgamma(1.0, sc);
    }
    // tau2 | nu_sigma_i[*], nu_tau
    {
      double sc_tau = 0.0;
      for (int i=0;i<N;++i) sc_tau += 1.0 / clip1(nu_sigma_i(i));
      sc_tau += 1.0 / clip1(nu_tau);
      tau2 = rinvgamma(0.5*(N+1.0), sc_tau);
    }
    // nu_tau | tau2, s2
    nu_tau = rinvgamma(1.0, 1.0/clip1(tau2) + 1.0/clip1(s2));

    // ============================================================
    // Store draws after burn-in
    // ============================================================
    if (it >= burn) {
      int m = it - burn;
      rho_draws[m] = rho;
      s2_draws[m]  = s2;

      if (K>0) beta0_draws.row(m) = beta0.t();

      if (p>0) {
        Lambda_draws.slice(m) = Eta;
        F_draws.slice(m)      = Gamma;
      }

      // alpha back to original scale of Y
      arma::vec alpha_unscaled = alpha / sds_Yc;
      alpha_draws.row(m) = alpha_unscaled.t();

      if (verbose && ((m+1) % 2000 == 0)) {
        Rcpp::checkUserInterrupt();
      }
    }

  } // end main MCMC loop

  // acceptance rates
  double acc_rho_rate   = acc_rho   / std::max(1, M);
  double acc_alpha_rate = acc_alpha / std::max(1, M);

  return Rcpp::List::create(
    _["rho"]        = rho_draws,
    _["beta"]       = beta0_draws,
    _["sigma2"]     = s2_draws,
    _["Lambda"]     = Lambda_draws,
    _["F"]          = F_draws,
    _["alpha"]      = alpha_draws,
    _["acc_rho"]    = acc_rho_rate,
    _["acc_alpha"]  = acc_alpha_rate
  );
}

// ---------------------------------------------------------------
// Gibbs sampler for alpha under horseshoe prior (BSCM part only)
//   y = X alpha + e,  e ~ N(0, s2 I)
//   alpha_i ~ N(0, sigma2_i), horseshoe via Makalic–Schmidt
//   ※ X の各列を標準化してからサンプリングし、返却時に元スケールへ戻す
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::mat hs_alpha_gibbs_cpp(const arma::vec& Y0_pre,              // T0
                             const arma::mat& control_outcome_pre, // T0 x N (X)
                             int iteration,
                             int burn,
                             bool verbose = false) {

  Rcpp::RNGScope scope;

  const int T0 = Y0_pre.n_elem;
  const int N  = control_outcome_pre.n_cols;
  if (control_outcome_pre.n_rows != (unsigned)T0) {
    stop("Dimension mismatch: nrow(Yc_pre) must equal length(Y0_pre).");
  }
  const int M = std::max(0, iteration - burn);

  // ----- scale X (column-wise) -----
  arma::mat X = control_outcome_pre;
  arma::vec sds_X(N, fill::zeros);
  for (int j = 0; j < N; ++j) {
    double sdj = arma::stddev(X.col(j));
    if (!arma::is_finite(sdj) || sdj < 1e-8) sdj = 1.0;
    sds_X(j) = sdj;
    X.col(j) = X.col(j) / sdj;
  }

  // precompute XtX, Xty
  arma::mat XtX = X.t() * X;            // N x N
  arma::vec Xty = X.t() * Y0_pre;       // N

  // ----- storage -----
  arma::mat alpha_draws(M, N, fill::zeros);

  // ----- states -----
  arma::vec alpha = 1e-4 * arma::randn<arma::vec>(N);
  double s2 = 1.0;

  // horseshoe (Makalic–Schmidt)
  arma::vec sigma2_i(N, fill::ones);   // local λ_i^2
  arma::vec nu_sigma_i(N, fill::ones); // ν_{λ_i}
  double tau2 = 1.0, nu_tau = 1.0;     // global τ^2, ν_τ

  // s2 hyper
  double nu_sigma2 = 1.0;

  // small ridge for numeric stability
  const double ridge = 1e-10;

  for (int it = 0; it < iteration; ++it) {

    // ---- (1) alpha | rest  ~  N(m, V) ----
    // V^{-1} = XtX / s2 + Diag(1/sigma2_i)
    arma::vec inv_sig2 = 1.0 / clip_vec(sigma2_i);
    arma::mat Prec = XtX / clip1(s2);
    Prec.diag() += inv_sig2;
    Prec.diag() += ridge; // numeric guard
    symmetrize_inplace(Prec);

    arma::mat Prec_chol;
    if (!robust_chol(Prec_chol, Prec)) {
      stop("Cholesky(Precision) failed in alpha-step.");
    }
    // mean = V * (X'y / s2)  with V = Prec^{-1}
    // solve(Prec, Xty/s2)
    arma::vec rhs = Xty / clip1(s2);
    arma::vec m = solve(trimatu(Prec_chol.t()), solve(trimatl(Prec_chol), rhs));

    // sample from N(m, V): solve(Prec, z) + m, z~N(0,I)
    arma::vec z = arma::randn<arma::vec>(N);
    arma::vec v = solve(trimatu(Prec_chol.t()), z); // v ~ N(0, V)
    alpha = m + v;

    // ---- (2) s2 | rest  ~ IG ----
    arma::vec resid = Y0_pre - X * alpha;
    double ss = arma::dot(resid, resid);
    double shape = 1.0 + 0.5 * T0;              // a0=1 と同等（弱情報）
    double rate  = 1.0 / clip1(nu_sigma2) + 0.5 * ss;
    s2 = rinvgamma(shape, rate);
    // hyper for s2
    nu_sigma2 = rinvgamma(1.0, 1.0/clip1(s2) + 1.0/100.0);

    // ---- (3) horseshoe locals & global ----
    // sigma2_i | alpha, nu_sigma_i
    for (int j = 0; j < N; ++j) {
      double sc = 0.5 * alpha(j) * alpha(j) + 1.0 / clip1(nu_sigma_i(j));
      sigma2_i(j) = rinvgamma(1.0, sc);
    }
    // nu_sigma_i | sigma2_i, tau2
    for (int j = 0; j < N; ++j) {
      double sc = 1.0 / clip1(sigma2_i(j)) + 1.0 / clip1(tau2);
      nu_sigma_i(j) = rinvgamma(1.0, sc);
    }
    // tau2 | nu_sigma_i[*], nu_tau
    {
      double sc_tau = 0.0;
      for (int j = 0; j < N; ++j) sc_tau += 1.0 / clip1(nu_sigma_i(j));
      sc_tau += 1.0 / clip1(nu_tau);
      tau2 = rinvgamma(0.5 * (N + 1.0), sc_tau);
    }
    // nu_tau | tau2, s2（論文整合のため s2 を参照）
    nu_tau = rinvgamma(1.0, 1.0/clip1(tau2) + 1.0/clip1(s2));

    // ---- store ----
    if (it >= burn) {
      int m_ix = it - burn;
      // 元スケールへ戻す：alpha_original = alpha_scaled / sd(X_j)
      arma::vec alpha_unscaled = alpha / sds_X;
      alpha_draws.row(m_ix) = alpha_unscaled.t();

      if (verbose && ((m_ix + 1) % 2000 == 0)) {
        Rcpp::checkUserInterrupt();
      }
    }
  }

  return alpha_draws;
}