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

// ===================== 2) Full SAR with factors (+ alpha) ===================
//
//  Ã(ρ,α) * y_{c,t} = X_t β0 + Λ f_t + ε_t
//  Ã(ρ,α) = I_N − ρ W − ρ w αᵀ
//
// 事後カーネル（ρ,α 共通）:
//  T0 * log|det(Ã)| − (1/2σ²) ∑_t || Ã y_{c,t} − X_t β0 − Λ f_t ||²  + log p(α | HS)
//
// [[Rcpp::export]]
Rcpp::List sar_full_sampler_cpp(const arma::vec& Y0_pre,    // 未使用（旧モデルの ρ a y0 は捨象）
                                const arma::mat& Yc_pre,    // T0 x N
                                Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_, // (T0*N*K), idx=(t*N+i)*K+k
                                int T0, int N, int K, int p,
                                const arma::vec& w,         // N
                                const arma::mat& W,         // N x N
                                int iteration, int burn,
                                double step_rho = 0.01,
                                double step_alpha = 0.01,   // ← 追加: α の RW–MH ステップ幅
                                double a0 = 1.0, double b0 = 1.0,
                                bool verbose=false) {

  Rcpp::RNGScope scope;

  const int M = std::max(0, iteration - burn);

  // --- data
  arma::mat Yc_orig = Yc_pre; // store original data

  arma::vec sds_Yc(N); // 各対照ユニットの標準偏差を格納
  arma::mat Yc = Yc_orig; // スケーリング用のコピーを作成
  for (int j = 0; j < N; ++j) {
    double sdj = arma::stddev(Yc_orig.col(j));
    if (!arma::is_finite(sdj) || sdj < 1e-8) sdj = 1.0; // 0割やNaNを回避
    sds_Yc(j) = sdj;
    Yc.col(j) = Yc_orig.col(j) / sdj; // ★ Yc を標準偏差でスケーリング
  }

  // --- X accessor (idx=(t*N+i)*K+k)
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

  // --- spectral bound for rho (W の固有値に基づく保守的境界)
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs = 0.0;
  for (arma::uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // --- storage
  arma::vec rho_draws(M, arma::fill::none);
  arma::vec s2_draws(M, arma::fill::none);
  arma::mat beta0_draws(M, K, arma::fill::zeros);
  arma::cube Lambda_draws(N, p, M, arma::fill::zeros); // = Eta
  arma::cube F_draws(p, T0, M, arma::fill::zeros);     // = Gamma
  arma::mat alpha_draws(M, N, arma::fill::zeros);

  // --- states
  double rho = 0.0;
  double s2  = 1.0;

  // beta0 ~ HS（Makalic–Schmidt）
  arma::vec beta0 = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());
  arma::vec sig2_b0   = (K>0 ? arma::ones<arma::vec>(K) : arma::vec());   // λ_j^2
  arma::vec nu_sig_b0 = (K>0 ? arma::ones<arma::vec>(K) : arma::vec());   // ν_j
  double tau2_b0 = 1.0, nu_tau_b0 = 1.0;

  // sigma^2 のハイパー
  double nu_sigma2 = 1.0;

  // 因子: Eta (N x p), Gamma (p x T0)
  arma::mat Eta   = (p>0 ? arma::zeros<arma::mat>(N,p) : arma::mat());
  arma::mat Gamma = (p>0 ? arma::zeros<arma::mat>(p,T0) : arma::mat());
  double phi_g = 0.0, s2_g = 1.0, nu_s2_g = 1.0;
  arma::vec omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  arma::vec nu_omega_k = (p>0 ? arma::ones<arma::vec>(p) : arma::vec());
  double s2_eta = 1.0, nu_s2_eta = 1.0;

  // -------- α（提案手法の核） --------
  arma::vec alpha = 1e-4 * arma::randn<arma::vec>(N);     // 初期値は微小乱数
  // Horseshoe の階層パラメータ（Makalic–Schmidt）
  arma::vec sigma2_i = arma::ones<arma::vec>(N);          // 局所スケール λ_i^2
  arma::vec nu_sigma_i = arma::ones<arma::vec>(N);        // ν_{λ_i}
  double tau2 = 1.0, nu_tau = 1.0;                        // 大域スケール τ^2, ν_τ

  arma::mat I_N = arma::eye(N,N);
  arma::mat I_K = (K>0 ? arma::eye(K,K) : arma::mat());
  arma::mat I_p = (p>0 ? arma::eye(p,p) : arma::mat());

  auto Mtilde = [&](double r, const arma::vec& a)->arma::mat {
    // Ã(r,α) = I − r W − r w αᵀ
    return I_N - r * W - r * w * a.t();
  };

  auto loglik_core = [&](double r, const arma::vec& a)->std::pair<double,double> {
    if (std::abs(r) >= bnd) return {-std::numeric_limits<double>::infinity(), 0.0};
    arma::mat M = Mtilde(r, a);
    double ldet = logdet(M);
    if (!std::isfinite(ldet)) return {-std::numeric_limits<double>::infinity(), 0.0};
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = M * Yc.row(t).t();
      if (useX)  u -= X_get_row(t) * beta0;
      if (p>0)   u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    double ll = T0 * ldet - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
    return {ll, ss};
  };

  auto logpost_alpha = [&](const arma::vec& a)->double {
    // 共同尤度 + HS 事前（正規 N(0, diag(sigma2_i))）の対数（定数は省略）
    auto pr = loglik_core(rho, a);
    if (!std::isfinite(pr.first)) return -std::numeric_limits<double>::infinity();
    double lprior = 0.0;
    for (int i=0;i<N;++i) {
      double s2i = clip1(sigma2_i(i));
      lprior += -0.5 * std::log(s2i) - 0.5 * (a(i)*a(i))/s2i;
    }
    return pr.first + lprior;
  };

  int iters = iteration;
  int acc_rho = 0, acc_alpha = 0;

  double log_step_rho = std::log(step_rho); // 初期ステップ幅の対数
  double target_accept_rho = 0.234;       // 目標受理率 (1次元の場合)
  double adapt_gamma = 0.6;              // 適応の減衰率 (0.5 < gamma <= 1)

  arma::vec log_step_alpha(N);
  log_step_alpha.fill(std::log(step_alpha)); // Use the input 'step_alpha' as initial
  double target_accept_alpha = 0.44;       // Target for multi-dim RWM (or 0.234)
  // double adapt_gamma = 0.6; // Reuse from rho adaptation
  arma::ivec acc_alpha_counts = arma::zeros<arma::ivec>(N);

  for (int it=0; it<iters; ++it) {

    // ===== (1) gamma_t | rest (FFBS; 観測式は Ã を用いる) =====
    if (p > 0) {
      arma::mat Phi = I_p * phi_g;
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

      double den = 0.0, num = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        den += arma::dot(gl, gl);
        num += arma::dot(gl, Gamma.col(t));
      }
      double mean_phi = (den > 0 ? num / den : 0.0);
      double var_phi  = (den > 0 ? s2_g / den : 1.0);
      double cand; do { cand = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand) > 1.0);
      phi_g = cand;

      double sc = 0.0;
      for (int t=0; t<T0; ++t) {
        arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
        arma::vec diff = Gamma.col(t) - phi_g * gl;
        sc += 0.5 * arma::dot(diff, diff);
      }
      s2_g     = rinvgamma(0.5 + 0.5 * p * T0, sc + 1.0/clip1(nu_s2_g));
      nu_s2_g  = rinvgamma(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
    }

    // ===== (2) Eta | rest =====
    if (p>0) {
      arma::mat GtG = Gamma * Gamma.t();
      arma::mat Domega = arma::diagmat(omega_k);
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
      double sc = 0.0;
      for (int i=0;i<N;++i) { arma::vec ei = Eta.row(i).t(); sc += arma::dot(ei, Domega * ei); }
      s2_eta    = rinvgamma(0.5 + 0.5 * p * N, 0.5*sc + 1.0/clip1(nu_s2_eta));
      nu_s2_eta = rinvgamma(1.0, 1.0/clip1(s2_eta) + 1.0/100.0);
      for (int k=0;k<p;++k) {
        double tmp=0.0; for (int i=0;i<N;++i) tmp += 0.5 * Eta(i,k)*Eta(i,k) / clip1(s2_eta);
        double rate_ok = 1.0/clip1(nu_omega_k(k)) + tmp;
        omega_k(k)     = rinvgamma(0.5*(N+1.0), rate_ok);
        nu_omega_k(k)  = rinvgamma(1.0, 1.0 + 1.0/clip1(omega_k(k)));
      }
    }

    // ===== (3) beta0 | rest （HS） =====
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
      arma::vec m    = Ainv * Bb;
      arma::mat S    = clip1(s2) * Ainv;
      beta0 = arma::mvnrnd(m, 0.5*(S+S.t()), 1);

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

    // ===== (4) sigma^2 | rest =====
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

    // // ===== (5) rho | rest （RW–MH; Ã を用いた尤度） =====
    // auto lpair_cur = loglik_core(rho, alpha);
    // double lcur_rho = lpair_cur.first;
    // // double prop_rho = R::rnorm(rho, step_rho);
    // double prop_rho = R::runif(-0.99, 0.99);
    // Rcpp::Rcout << prop_rho << "\n";
    // double lprp_rho = loglik_core(prop_rho, alpha).first;
    // if (std::log(R::runif(0.0,1.0)) < (lprp_rho - lcur_rho)) {Rcpp::Rcout << "Accept!" << prop_rho << "\n";}
    // if (std::log(R::runif(0.0,1.0)) < (lprp_rho - lcur_rho)) { rho = prop_rho; if (it>=burn) acc_rho++; }

    // ===== (5) rho | rest （Adaptive RW–MH） =====
    double current_step_rho = std::exp(log_step_rho); // 現在のステップ幅
    double prop_rho = R::rnorm(rho, current_step_rho);

    auto lpair_cur = loglik_core(rho, alpha); // 現在の対数尤度+ヤコビアン
    double lcur_rho = lpair_cur.first;
    double lprp_rho = loglik_core(prop_rho, alpha).first; // 提案点の対数尤度+ヤコビアン

    bool accepted_rho = false;
    double log_alpha_accept = lprp_rho - lcur_rho;

    if (std::log(R::runif(0.0, 1.0)) < log_alpha_accept) {
      rho = prop_rho;
      accepted_rho = true;
      if (it >= burn) acc_rho++;
    }

    // --- ステップ幅の適応更新 (Roberts & Rosenthal style) ---
    // バーンイン中にのみ適応させるか、全体で適応させるかは選択肢あり
    // ここでは単純化のため、毎回更新する例を示す
    double adapt_step = std::pow(it + 1.0, -adapt_gamma); // 適応量を徐々に減らす
    log_step_rho += adapt_step * ( (accepted_rho ? 1.0 : 0.0) - target_accept_rho );
    // log_step_rho の値が極端にならないように制限を加える方が安全な場合がある
    // log_step_rho = std::max(-10.0, std::min(log_step_rho, 2.0)); // 例

      // ============================================================
  //  (6) alpha | rest
  // ============================================================
  // ===== (6) alpha | rest  — Independent MH =====
// 1) 提案分布 q(α) = N(μ_q, Σ_q) の計算
//    R_t = (I - ρW) y_{c,t} - X_t β0 - Λ γ_t
//    Z_t = ρ (w y_{c,t}^T)    （N×N のランク1行列）
//    ガウス部分（尤度×HS事前）による事後の精度：
//      Prec_q = (1/s²) * Σ_t Z_t'Z_t + Diag(1/σ_i²)
//             = (ρ² * ||w||² / s²) * Σ_t (y_t y_t') + Diag(1/σ_i²)
//    平均：
//      μ_q = Σ_q * b_q,   b_q = -(1/s²) * Σ_t Z_t' R_t
//          = -(ρ / s²) * Σ_t y_t * (w^T R_t)          （N×1）
// Rcpp::Rcout << "===Hello World!!===" << '\n';
// {
//   // 事前の対角（HS）
//   arma::vec inv_sigma2_i = 1.0 / clip_vec(sigma2_i);

//   // 事前計算：||w||², および Σ_t y_t y_t'（N×N, ランク ≤ T0）
//   double wnorm2 = arma::dot(w, w);
//   arma::mat Syy(N, N, arma::fill::zeros);
//   arma::vec bq(N, arma::fill::zeros);  // b_q

//   for (int t=0; t<T0; ++t) {
//     arma::vec y_t = Yc.row(t).t();
//     // R_t = (I - ρW) y_t - X_t β0 - Λ γ_t
//     arma::vec R_t = (I_N - rho * W) * y_t;
//     if (useX) R_t -= X_get_row(t) * beta0;
//     if (p>0)  R_t -= Eta * Gamma.col(t);

//     Syy += y_t * y_t.t();                              // Σ y_t y_t'
//     double wtR = arma::dot(w, R_t);                    // w^T R_t （スカラー）
//     bq   -= (rho / clip1(s2)) * (y_t * wtR);           // b_q に加算
//   }

//   // Prec_q = (ρ² * ||w||² / s²) * Syy + Diag(1/σ_i²)
//   arma::mat Prec = ((rho*rho) * wnorm2 / clip1(s2)) * Syy;
//   Prec.diag() += inv_sigma2_i;
//   Prec.diag() += 1e-1; // for stability

//   // 数値安定化して共分散 Σ_q を得る
//   arma::mat Sig_q;
//   {
//     // Prec は対称正定値のはずだが、念のため安定化
//     arma::mat Prec_sym = 0.5 * (Prec + Prec.t());
//     arma::mat Rchol;
//     if (!robust_chol(Rchol, Prec_sym)) Rcpp::stop("chol(Prec_q) failed.");
//     // Σ_q = Prec_q^{-1} だが、mvnrnd に直接 Prec は使えないので逆行列を作る
//     Sig_q = arma::inv_sympd(Prec_sym);
//   }

//   // μ_q = Σ_q * b_q
//   arma::vec mu_q = Sig_q * bq;

//   // 2) 提案：α_prop ~ N(μ_q, Σ_q)
//   arma::vec alpha_prop = arma::mvnrnd(mu_q, 0.5*(Sig_q + Sig_q.t()), 1);

//   // 3) 受理判定：ヤコビアンだけ（q はガウス項＋HS事前を包含）
//   //    r = min{1, J(α_prop,ρ)/J(α,ρ)},  J = |det(I - ρW - ρ w α^T)|^{T0}
//   auto logJ = [&](const arma::vec& a)->double {
//     arma::mat M = I_N - rho * W - rho * w * a.t();
//     double ldet = logdet(M);
//     if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
//     return T0 * ldet;
//   };

//   double lcur = logJ(alpha);
//   double lprp = logJ(alpha_prop);
//   if (it == 100){
//   // Rcpp::Rcout << "===Start===" << '\n';
//   // Rcpp::Rcout << "alpha:" << alpha << '\n';
//   // Rcpp::Rcout << "alpha_prop:" << alpha_prop << '\n';
//   // Rcpp::Rcout << "mu_q:" << mu_q << '\n';
//   // Rcpp::Rcout << "diagvec(Sig_q):" << diagvec(Sig_q) << '\n';
//   // Rcpp::Rcout << "lprp - lcur:" << lprp - lcur << '\n';
//   // Rcpp::Rcout << "Syy:" << Syy << '\n';
//   // Rcpp::Rcout << "===End===" << '\n';
//   // Rcpp::Rcout << '\n';
// }
//   if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) {
//     alpha = alpha_prop;
//   }
// }

// ===== (6) alpha | rest (Component-wise Adaptive RW–MH) =====
    int current_acc_alpha_iter = 0; // Track acceptance in this iteration
    for (int i = 0; i < N; ++i) {
      // (a) Get current log posterior
      double lcur_a = logpost_alpha(alpha); // log(Likelihood * Jacobian * Prior)

      // (b) Propose a change only for alpha[i]
      arma::vec prop_alpha = alpha;
      double current_step_i = std::exp(log_step_alpha(i));
      prop_alpha(i) = alpha(i) + R::rnorm(0.0, current_step_i);

      // (c) Get log posterior at proposal
      double lprp_a = logpost_alpha(prop_alpha);

      // (d) Accept / Reject
      bool accepted_i = false;
      double log_alpha_accept_i = lprp_a - lcur_a;
      if (std::log(R::runif(0.0, 1.0)) < log_alpha_accept_i) {
        alpha = prop_alpha; // Update the *whole* alpha vector if accepted
        accepted_i = true;
        if (it >= burn) acc_alpha_counts(i)++; // Count acceptance for component i
        current_acc_alpha_iter++; // Count for overall rate in this iteration
      }

      // (e) Adapt step size for alpha[i]
      double adapt_step_i = std::pow(it + 1.0, -adapt_gamma);
      log_step_alpha(i) += adapt_step_i * ( (accepted_i ? 1.0 : 0.0) - target_accept_alpha );
      // Optional: Clamp log_step_alpha(i) to prevent extreme values
      // log_step_alpha(i) = std::max(-10.0, std::min(log_step_alpha(i), 3.0));
    }
    // Store overall acceptance rate for this iteration if needed
    if (it >= burn) acc_alpha += current_acc_alpha_iter; // Total acceptances across all components
    // ===== (7) HS hyper for alpha （Gibbs） =====
    // sigma2_i[i] ~ IG(1, 0.5*alpha[i]^2 + 1/nu_sigma_i[i])
    for (int i=0;i<N;++i) {
      double sc = 0.5 * alpha(i)*alpha(i) + 1.0 / clip1(nu_sigma_i(i));
      sigma2_i(i) = rinvgamma(1.0, sc);
    }
    // nu_sigma_i[i] ~ IG(1, 1/sigma2_i[i] + 1/tau2)
    for (int i=0;i<N;++i) {
      double sc = 1.0 / clip1(sigma2_i(i)) + 1.0 / clip1(tau2);
      nu_sigma_i(i) = rinvgamma(1.0, sc);
    }
    // tau2 ~ IG((N+1)/2, sum(1./nu_sigma_i) + 1/nu_tau)
    {
      double sc_tau = 0.0; for (int i=0;i<N;++i) sc_tau += 1.0 / clip1(nu_sigma_i(i));
      sc_tau += 1.0 / clip1(nu_tau);
      tau2 = rinvgamma(0.5*(N+1.0), sc_tau);
    }
    // nu_tau ~ IG(1, 1/tau2 + 1/sigma2)   ← 論文側の一貫性のため s2 を参照
    nu_tau = rinvgamma(1.0, 1.0/clip1(tau2) + 1.0/clip1(s2));

    // ===== store after burn =====
    if (it >= burn) {
      int m = it - burn;
      rho_draws[m] = rho;
      s2_draws[m]  = s2;
      if (K>0) beta0_draws.row(m) = beta0.t();
      if (p>0) {
        Lambda_draws.slice(m) = Eta;
        F_draws.slice(m)      = Gamma;
      }
      arma::vec alpha_unscaled = alpha / sds_Yc;
      alpha_draws.row(m) = alpha_unscaled.t();
      if (verbose && ((m+1) % 2000 == 0)) Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(
    _["rho"]      = rho_draws,
    _["beta"]     = beta0_draws,
    _["sigma2"]   = s2_draws,
    _["Lambda"]   = Lambda_draws,
    _["F"]        = F_draws,
    _["alpha"]    = alpha_draws,
    _["acc_rho"]  = acc_rho / std::max(1, M),
    _["acc_alpha"]= static_cast<double>(acc_alpha) / std::max(1, M * N),
    _["final_log_step_alpha"] = Rcpp::wrap(log_step_alpha)
  );
}

// [[Rcpp::export]]
Rcpp::List sar_full_one_step_cpp(const arma::vec& Y0_pre,    // 未使用
                                 const arma::mat& Yc_pre,    // T0 x N
                                 Rcpp::Nullable<Rcpp::NumericVector> Xc_pre_,
                                 int T0, int N, int K, int p,
                                 const arma::vec& w,
                                 const arma::mat& W,         // N x N
                                 Rcpp::List state,
                                 double step_rho = 0.01,
                                 double step_alpha = 0.05,
                                 double a0 = 1.0, double b0 = 1.0,
                                 bool verbose=false) {
  Rcpp::RNGScope scope;

  // ---- unpack state (必須) ----
  double rho    = Rcpp::as<double>(state["rho"]);
  double s2     = Rcpp::as<double>(state["sigma2"]);
  arma::vec beta = (K>0 && state.containsElementNamed("beta")) ? Rcpp::as<arma::vec>(state["beta"]) : arma::vec();
  arma::mat Eta   = (p>0 && state.containsElementNamed("Lambda")) ? Rcpp::as<arma::mat>(state["Lambda"]) : arma::mat();
  arma::mat Gamma = (p>0 && state.containsElementNamed("F")) ? Rcpp::as<arma::mat>(state["F"]) : arma::mat();

  // ---- HS (beta) 補助 ----
  arma::vec sig2_b0   = (K>0 && state.containsElementNamed("sig2_b0"))  ? Rcpp::as<arma::vec>(state["sig2_b0"])  : arma::ones<arma::vec>(K);
  arma::vec nu_sig_b0 = (K>0 && state.containsElementNamed("nu_sig_b0"))? Rcpp::as<arma::vec>(state["nu_sig_b0"]): arma::ones<arma::vec>(K);
  double tau2_b0      = (K>0 && state.containsElementNamed("tau2_b0"))  ? Rcpp::as<double>(state["tau2_b0"])     : 1.0;
  double nu_tau_b0    = (K>0 && state.containsElementNamed("nu_tau_b0"))? Rcpp::as<double>(state["nu_tau_b0"])   : 1.0;

  double nu_sigma2    = state.containsElementNamed("nu_sigma2") ? Rcpp::as<double>(state["nu_sigma2"]) : 1.0;

  // ---- 因子の超パラ ----
  double phi_g        = (p>0 && state.containsElementNamed("phi_g"))     ? Rcpp::as<double>(state["phi_g"])     : 0.0;
  double s2_g         = (p>0 && state.containsElementNamed("s2_g"))      ? Rcpp::as<double>(state["s2_g"])      : 1.0;
  double nu_s2_g      = (p>0 && state.containsElementNamed("nu_s2_g"))   ? Rcpp::as<double>(state["nu_s2_g"])   : 1.0;
  arma::vec omega_k   = (p>0 && state.containsElementNamed("omega_k"))   ? Rcpp::as<arma::vec>(state["omega_k"]) : arma::ones<arma::vec>(p);
  arma::vec nu_omega_k= (p>0 && state.containsElementNamed("nu_omega_k"))? Rcpp::as<arma::vec>(state["nu_omega_k"]) : arma::ones<arma::vec>(p);
  double s2_eta       = (p>0 && state.containsElementNamed("s2_eta"))    ? Rcpp::as<double>(state["s2_eta"])    : 1.0;
  double nu_s2_eta    = (p>0 && state.containsElementNamed("nu_s2_eta")) ? Rcpp::as<double>(state["nu_s2_eta"]) : 1.0;

  // ---- α と HS(α) ----
  arma::vec alpha     = state.containsElementNamed("alpha") ? Rcpp::as<arma::vec>(state["alpha"]) : 1e-4*arma::randn<arma::vec>(N);
  arma::vec sigma2_i  = state.containsElementNamed("sigma2_i") ? Rcpp::as<arma::vec>(state["sigma2_i"]) : arma::ones<arma::vec>(N);
  arma::vec nu_sigma_i= state.containsElementNamed("nu_sigma_i") ? Rcpp::as<arma::vec>(state["nu_sigma_i"]) : arma::ones<arma::vec>(N);
  double tau2         = state.containsElementNamed("tau2") ? Rcpp::as<double>(state["tau2"]) : 1.0;
  double nu_tau       = state.containsElementNamed("nu_tau") ? Rcpp::as<double>(state["nu_tau"]) : 1.0;

  // ---- data, helpers ----
  arma::mat Yc = Yc_pre;
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

  arma::mat I_N = arma::eye(N,N);
  auto Mtilde = [&](double r, const arma::vec& a)->arma::mat { return I_N - r * W - r * w * a.t(); };

  // ---- spectral bound ----
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs = 0.0; for (arma::uword i=0;i<evals.n_elem;++i) maxabs = std::max(maxabs, std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  auto loglik_core = [&](double r, const arma::vec& a)->double {
    if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity();
    arma::mat M = Mtilde(r, a);
    double ldet = logdet(M);
    if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    double ss = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec u = M * Yc.row(t).t();
      if (useX)  u -= X_get_row(t) * beta;
      if (p>0)   u -= Eta * Gamma.col(t);
      ss += arma::dot(u,u);
    }
    return T0 * ldet - 0.5 * (N*T0) * std::log(s2) - 0.5 * ss / s2;
  };
  auto logpost_alpha = [&](const arma::vec& a)->double {
    double ll = loglik_core(rho, a);
    if (!std::isfinite(ll)) return ll;
    double lp=0.0; for (int i=0;i<N;++i){ double s2i=clip1(sigma2_i(i)); lp += -0.5*std::log(s2i) - 0.5*(a(i)*a(i))/s2i; }
    return ll + lp;
  };

  // ============================================================
  //  (1) gamma_t | rest  （p>0 の場合; FFBS）
  // ============================================================
  if (p > 0) {
    arma::mat I_p = arma::eye(p,p);
    arma::mat Phi = I_p * phi_g;
    arma::mat Q_mat = I_p * s2_g;
    arma::mat H_mat = Eta;

    arma::mat HtH_s2_inv = (H_mat.t() * H_mat) / clip1(s2);
    arma::mat Ht_s2_inv  = H_mat.t() / clip1(s2);

    std::vector<arma::vec> gamma_t_t(T0), gamma_t_t1(T0);
    std::vector<arma::mat> P_t_t(T0), P_t_t1(T0);

    arma::vec gprev = arma::zeros<arma::vec>(p);
    arma::mat Pprev = (clip1(s2_g) / std::max(1e-6, 1.0 - phi_g*phi_g)) * I_p;

    for (int t=0; t<T0; ++t) {
      arma::vec Y_star_t = Mtilde(rho, alpha) * Yc.row(t).t();
      if (useX) Y_star_t -= X_get_row(t) * beta;

      arma::vec pred = Phi * gprev;
      arma::mat Ppred = Phi * Pprev * Phi.t() + Q_mat;

      arma::mat Pinv = arma::inv_sympd(Ppred);
      arma::mat Vinv = Pinv + HtH_s2_inv;
      arma::mat V = arma::inv_sympd(Vinv);
      arma::vec m = V * (Pinv * pred + Ht_s2_inv * Y_star_t);

      gamma_t_t[t]=m; P_t_t[t]=V; gamma_t_t1[t]=pred; P_t_t1[t]=Ppred;
      gprev=m; Pprev=V;
    }
    arma::mat L = arma::chol(0.5*(P_t_t[T0-1]+P_t_t[T0-1].t()), "lower");
    Gamma.col(T0-1) = gamma_t_t[T0-1] + L * arma::randn<arma::vec>(p);
    for (int t=T0-2; t>=0; --t) {
      arma::mat P_inv_next = arma::inv_sympd(P_t_t1[t+1]);
      arma::mat J = P_t_t[t] * (Phi.t()) * P_inv_next;
      arma::vec m = gamma_t_t[t] + J * (Gamma.col(t+1) - gamma_t_t1[t+1]);
      arma::mat V = P_t_t[t] - J * Phi * P_t_t[t];
      arma::mat Ls = arma::chol(0.5*(V+V.t()), "lower");
      Gamma.col(t) = m + Ls * arma::randn<arma::vec>(p);
    }
    double den=0.0, num=0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
      den += arma::dot(gl, gl); num += arma::dot(gl, Gamma.col(t));
    }
    double mean_phi = (den>0 ? num/den : 0.0);
    double var_phi  = (den>0 ? s2_g/den : 1.0);
    double cand; do { cand = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand)>1.0);
    phi_g = cand;

    double sc=0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec gl = (t==0) ? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1));
      arma::vec d = Gamma.col(t) - phi_g * gl;
      sc += 0.5 * arma::dot(d,d);
    }
    s2_g     = rinvgamma(0.5 + 0.5*p*T0, sc + 1.0/clip1(nu_s2_g));
    nu_s2_g  = rinvgamma(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
  }

  // ============================================================
  //  (2) Eta | rest  （p>0）
  // ============================================================
  if (p>0) {
    arma::mat GtG = Gamma * Gamma.t();
    arma::mat Domega = arma::diagmat(omega_k);
    arma::mat Vrow = arma::inv_sympd( GtG / s2 + Domega / clip1(s2_eta) );
    arma::mat Lrow = arma::chol(Vrow, "lower");
    for (int i=0; i<N; ++i) {
      arma::vec rhs = arma::zeros<arma::vec>(p);
      for (int t=0; t<T0; ++t) {
        arma::vec r = Mtilde(rho, alpha) * Yc.row(t).t();
        if (useX) r -= X_get_row(t) * beta;
        rhs += Gamma.col(t) * r(i);
      }
      arma::vec m = Vrow * (rhs / s2);
      Eta.row(i) = (m + Lrow * arma::randn<arma::vec>(p)).t();
    }
    double sc=0.0; for (int i=0;i<N;++i){ arma::vec ei = Eta.row(i).t(); sc += arma::dot(ei, Domega*ei); }
    s2_eta    = rinvgamma(0.5 + 0.5*p*N, 0.5*sc + 1.0/clip1(nu_s2_eta));
    nu_s2_eta = rinvgamma(1.0, 1.0/clip1(s2_eta) + 1.0/100.0);
    for (int k=0;k<p;++k) {
      double tmp=0.0; for (int i=0;i<N;++i) tmp += 0.5 * Eta(i,k)*Eta(i,k) / clip1(s2_eta);
      double rate_ok = 1.0/clip1(nu_omega_k(k)) + tmp;
      omega_k(k)     = rinvgamma(0.5*(N+1.0), rate_ok);
      nu_omega_k(k)  = rinvgamma(1.0, 1.0 + 1.0/clip1(omega_k(k)));
    }
  }

  // ============================================================
  //  (3) beta | rest  （K>0; HS）
  // ============================================================
  if (useX) {
    arma::mat Ab = arma::zeros<arma::mat>(K,K);
    arma::vec Bb = arma::zeros<arma::vec>(K);
    for (int t=0; t<T0; ++t) {
      arma::mat Xt = X_get_row(t);
      Ab += Xt.t() * Xt;
      arma::vec Btmp = Mtilde(rho, alpha) * Yc.row(t).t();
      if (p>0) Btmp -= Eta * Gamma.col(t);
      Bb += Xt.t() * Btmp;
    }
    Ab.diag() += clip1(s2) * (1.0 / clip_vec(sig2_b0));
    arma::mat Ainv = arma::inv_sympd(Ab);
    arma::vec m    = Ainv * Bb;
    arma::mat S    = clip1(s2) * Ainv;
    beta           = arma::mvnrnd(m, 0.5*(S+S.t()), 1);

    for (int j=0;j<K;++j) {
      double rate_l = 0.5 * beta(j)*beta(j) + 1.0/clip1(nu_sig_b0(j));
      sig2_b0(j)    = rinvgamma(1.0, rate_l);
      double rate_nu= 1.0/clip1(sig2_b0(j)) + 1.0/clip1(tau2_b0);
      nu_sig_b0(j)  = rinvgamma(1.0, rate_nu);
    }
    double sum_inv_nu = 0.0; for (int j=0;j<K;++j) sum_inv_nu += 1.0/clip1(nu_sig_b0(j));
    tau2_b0    = rinvgamma(1.0, 1.0/clip1(nu_tau_b0) + sum_inv_nu);
    nu_tau_b0  = rinvgamma(1.0, 1.0/clip1(tau2_b0) + 1.0/clip1(s2));
  }

  // ============================================================
  //  (4) sigma^2 | rest
  // ============================================================
  double ss = 0.0;
  for (int t=0; t<T0; ++t) {
    arma::vec u = Mtilde(rho, alpha) * Yc.row(t).t();
    if (useX) u -= X_get_row(t) * beta;
    if (p>0)  u -= Eta * Gamma.col(t);
    ss += arma::dot(u,u);
  }
  double shape = a0 + 0.5 * (T0 * N);
  double rate  = 1.0/clip1(nu_sigma2) + 0.5 * ss;
  s2 = rinvgamma(shape, rate);
  nu_sigma2 = rinvgamma(1.0, 1.0/clip1(s2) + 1.0/100.0);

  // ============================================================
  //  (5) rho | rest
  // ============================================================
  double lcur = loglik_core(rho, alpha);
  double prop_r = R::rnorm(rho, step_rho);
  double lprp = loglik_core(prop_r, alpha);
  bool acc_r = false;
  if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) { rho = prop_r; acc_r = true; }

  // ============================================================
  //  (6) alpha | rest
  // ============================================================
  // ===== (6) alpha | rest  — Independent MH =====
// 1) 提案分布 q(α) = N(μ_q, Σ_q) の計算
//    R_t = (I - ρW) y_{c,t} - X_t β0 - Λ γ_t
//    Z_t = ρ (w y_{c,t}^T)    （N×N のランク1行列）
//    ガウス部分（尤度×HS事前）による事後の精度：
//      Prec_q = (1/s²) * Σ_t Z_t'Z_t + Diag(1/σ_i²)
//             = (ρ² * ||w||² / s²) * Σ_t (y_t y_t') + Diag(1/σ_i²)
//    平均：
//      μ_q = Σ_q * b_q,   b_q = -(1/s²) * Σ_t Z_t' R_t
//          = -(ρ / s²) * Σ_t y_t * (w^T R_t)          （N×1）
Rcpp::Rcout << "===Hello World!!===" << '\n';
{
  // 事前の対角（HS）
  arma::vec inv_sigma2_i = 1.0 / clip_vec(sigma2_i);

  // 事前計算：||w||², および Σ_t y_t y_t'（N×N, ランク ≤ T0）
  double wnorm2 = arma::dot(w, w);
  arma::mat Syy(N, N, arma::fill::zeros);
  arma::vec bq(N, arma::fill::zeros);  // b_q

  for (int t=0; t<T0; ++t) {
    arma::vec y_t = Yc.row(t).t();
    // R_t = (I - ρW) y_t - X_t β0 - Λ γ_t
    arma::vec R_t = (I_N - rho * W) * y_t;
    if (useX) R_t -= X_get_row(t) * beta;
    if (p>0)  R_t -= Eta * Gamma.col(t);

    Syy += y_t * y_t.t();                              // Σ y_t y_t'
    double wtR = arma::dot(w, R_t);                    // w^T R_t （スカラー）
    bq   -= (rho / clip1(s2)) * (y_t * wtR);           // b_q に加算
  }

  // Prec_q = (ρ² * ||w||² / s²) * Syy + Diag(1/σ_i²)
  arma::mat Prec = ((rho*rho) * wnorm2 / clip1(s2)) * Syy;
  Prec.diag() += inv_sigma2_i;

  // 数値安定化して共分散 Σ_q を得る
  arma::mat Sig_q;
  {
    // Prec は対称正定値のはずだが、念のため安定化
    arma::mat Prec_sym = 0.5 * (Prec + Prec.t());
    arma::mat Rchol;
    if (!robust_chol(Rchol, Prec_sym)) Rcpp::stop("chol(Prec_q) failed.");
    // Σ_q = Prec_q^{-1} だが、mvnrnd に直接 Prec は使えないので逆行列を作る
    Sig_q = arma::inv_sympd(Prec_sym);
  }

  // μ_q = Σ_q * b_q
  arma::vec mu_q = Sig_q * bq;

  // 2) 提案：α_prop ~ N(μ_q, Σ_q)
  arma::vec alpha_prop = arma::mvnrnd(mu_q, 0.5*(Sig_q + Sig_q.t()), 1);

  // 3) 受理判定：ヤコビアンだけ（q はガウス項＋HS事前を包含）
  //    r = min{1, J(α_prop,ρ)/J(α,ρ)},  J = |det(I - ρW - ρ w α^T)|^{T0}
  auto logJ = [&](const arma::vec& a)->double {
    arma::mat M = I_N - rho * W - rho * w * a.t();
    double ldet = logdet(M);
    if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    return T0 * ldet;
  };

  double lcur = logJ(alpha);
  double lprp = logJ(alpha_prop);
  Rcpp::Rcout << "===Start===" << '\n';
  Rcpp::Rcout << "alpha:" << alpha << '\n';
  Rcpp::Rcout << "alpha_prop:" << alpha_prop << '\n';
  Rcpp::Rcout << "mu_q:" << mu_q << '\n';
  Rcpp::Rcout << "diagvec(Sig_q):" << diagvec(Sig_q) << '\n';
  Rcpp::Rcout << "lprp - lcur:" << lprp - lcur << '\n';
  Rcpp::Rcout << "===End===" << '\n';
  Rcpp::Rcout << '\n';
  if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) {
    alpha = alpha_prop;
  }
}

  // ============================================================
  //  (7) HS hyper for alpha （Gibbs）
  // ============================================================
  for (int i=0;i<N;++i) {
    double sc = 0.5 * alpha(i)*alpha(i) + 1.0 / clip1(nu_sigma_i(i));
    sigma2_i(i) = rinvgamma(1.0, sc);
  }
  for (int i=0;i<N;++i) {
    double sc = 1.0 / clip1(sigma2_i(i)) + 1.0 / clip1(tau2);
    nu_sigma_i(i) = rinvgamma(1.0, sc);
  }
  {
    double sc_tau = 0.0; for (int i=0;i<N;++i) sc_tau += 1.0 / clip1(nu_sigma_i(i));
    sc_tau += 1.0 / clip1(nu_tau);
    tau2 = rinvgamma(0.5*(N+1.0), sc_tau);
  }
  nu_tau = rinvgamma(1.0, 1.0/clip1(tau2) + 1.0/clip1(s2));

  // ---- return updated state ----
  Rcpp::List out;
  out["rho"]      = rho;
  out["sigma2"]   = s2;
  out["acc_rho"]  = acc_r;

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

  out["nu_sigma2"] = nu_sigma2;

  if (p > 0) {
    out["Lambda"]     = Rcpp::wrap(Eta);
    out["F"]          = Rcpp::wrap(Gamma);
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

  // α と HS(α) の現在値を返す（外部でモニタ・継続可能）
  out["alpha"]      = Rcpp::wrap(alpha);
  out["sigma2_i"]   = Rcpp::wrap(sigma2_i);
  out["nu_sigma_i"] = Rcpp::wrap(nu_sigma_i);
  out["tau2"]       = tau2;
  out["nu_tau"]     = nu_tau;
  // out["acc_alpha"]  = acc_a;

  return out;
}