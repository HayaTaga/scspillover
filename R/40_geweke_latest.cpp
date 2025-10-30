// 40_geweke_full.cpp
// Geweke (2004) Joint Distribution Test – FULL version
// モデル： (I - rho * A) Yc_t = X_t beta + Eta * Gamma_t + eps_t
//          A = W + w * alpha_hat^T
//
// 提供：
//  (1) simulate_Yc_forward_cpp : 前進シミュレータ（X と因子を完全反映）
//  (2) scspill_one_step_cpp    : 1 ステップ posterior カーネル（実装一貫）
//
// 依存: RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//-------------------- utilities --------------------
inline double rinvgamma1(double shape, double rate) {
  // InvGamma(shape, rate): density ∝ rate^shape x^{-(shape+1)} exp(-rate/x)
  double g = R::rgamma(shape, 1.0 / rate); // Gamma(shape, scale=1/rate)
  return 1.0 / g;
}
inline double clip1(double z, double lo=1e-12, double hi=1e12) {
  if (z < lo) z = lo; if (z > hi) z = hi; return z;
}
inline vec clip_vec(const vec& v, double lo=1e-12, double hi=1e12) {
  vec out = v;
  for (uword i=0;i<out.n_elem;++i){ double x=out(i); if(x<lo)x=lo; if(x>hi)x=hi; out(i)=x; }
  return out;
}
inline bool solve_safe(vec& x, const mat& M, const vec& b) {
  try {
    x = solve(M, b, solve_opts::fast);
    return x.is_finite();
  } catch (...) { return false; }
}

//----------------------------------------------
// (1) 前進シミュレータ（MC 側に使用）
//----------------------------------------------
// [[Rcpp::export]]
arma::mat simulate_Yc_forward_cpp(
  int T0,
  const arma::mat& W_use,            // N x N
  const arma::vec& w_use,            // N
  const arma::vec& alpha_hat_scaled, // N
  double rho,
  double sigma2,
  const arma::cube& Xc_pre,          // T0 x N x K  (K==0 可)
  const arma::vec& beta,             // K ベクトル（K==0なら長さ0）
  const arma::mat& Eta,              // N x p（p==0 可）
  const arma::mat& Gamma             // p x T0
) {
  Rcpp::RNGScope scope;

  const int N = W_use.n_rows;
  const int K = static_cast<int>(Xc_pre.n_slices);
  const int p = static_cast<int>(Eta.n_cols);

  arma::mat Yc(T0, N, arma::fill::zeros);
  arma::mat A = W_use + w_use * alpha_hat_scaled.t();
  arma::mat I = arma::eye<arma::mat>(N, N);
  const double sd = std::sqrt(std::max(1e-12, sigma2));

  for (int t = 0; t < T0; ++t) {
    arma::vec mu(N, arma::fill::zeros);

    // ---- X_t beta ----
    if (K > 0 && beta.n_elem == static_cast<uword>(K)) {
      for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int k = 0; k < K; ++k) s += Xc_pre(t, i, k) * beta(k);
        mu(i) += s;
      }
    }

    // ---- factors ----
    if (p > 0) mu += Eta * Gamma.col(t);

    // eps ~ N(0, sigma2 I)
    arma::vec eps = sd * arma::randn<arma::vec>(N);

    // (I - rho A) Y_t = mu + eps
    arma::vec rhs = mu + eps;
    arma::vec Yt;
    if (!solve_safe(Yt, I - rho * A, rhs)) {
      arma::mat M = I - rho * A;
      arma::mat AtA = M.t() * M;  AtA.diag() += 1e-10;
      arma::vec Atr = M.t() * rhs;
      Yt = arma::solve(AtA, Atr, arma::solve_opts::fast);
    }
    Yc.row(t) = Yt.t();
  }
  return Yc;
}

//----------------------------------------------
// (2) 1 ステップ posterior カーネル（SC 側に使用）
//    ブロック：Gamma(FFBS)→Eta→beta→sigma2→rho
//----------------------------------------------
// [[Rcpp::export]]
Rcpp::List scspill_one_step_cpp(
  const arma::mat& Yc_data,          // T0 x N
  const arma::mat& W_use,            // N x N
  const arma::vec& w_use,            // N
  const arma::vec& alpha_hat_scaled, // N
  int T0,
  int N,
  const arma::cube& Xc_pre,          // T0 x N x K
  int K,
  int p,
  Rcpp::List state_in,               // list(rho,sigma2,beta,Eta,Gamma,phi_g,s2_g,nu_s2_g,omega_k,nu_omega_k,s2_eta,nu_s2_eta,nu_sigma2)
  double a0,
  double b0,
  double step_rho
) {
  Rcpp::RNGScope scope;

  // ---- unpack state ----
  double rho    = Rcpp::as<double>(state_in["rho"]);
  double sigma2 = Rcpp::as<double>(state_in["sigma2"]);
  arma::vec beta = (K > 0 ? Rcpp::as<arma::vec>(state_in["beta"]) : arma::vec());
  arma::mat Eta  = (p > 0 ? Rcpp::as<arma::mat>(state_in["Eta"])  : arma::mat(N, 0, arma::fill::zeros));
  arma::mat Gamma= (p > 0 ? Rcpp::as<arma::mat>(state_in["Gamma"]): arma::mat(0, T0, arma::fill::zeros));

  double phi_g   = state_in.containsElementNamed("phi_g")   ? Rcpp::as<double>(state_in["phi_g"])   : 0.0;
  double s2_g    = state_in.containsElementNamed("s2_g")    ? Rcpp::as<double>(state_in["s2_g"])    : 1.0;
  double nu_s2_g = state_in.containsElementNamed("nu_s2_g") ? Rcpp::as<double>(state_in["nu_s2_g"]) : 1.0;
  arma::vec omega_k    = (p>0 && state_in.containsElementNamed("omega_k"))    ? Rcpp::as<arma::vec>(state_in["omega_k"])    : arma::ones<arma::vec>(p);
  arma::vec nu_omega_k = (p>0 && state_in.containsElementNamed("nu_omega_k")) ? Rcpp::as<arma::vec>(state_in["nu_omega_k"]) : arma::ones<arma::vec>(p);
  double s2_eta  = state_in.containsElementNamed("s2_eta")  ? Rcpp::as<double>(state_in["s2_eta"])  : 1.0;
  double nu_s2_eta= state_in.containsElementNamed("nu_s2_eta") ? Rcpp::as<double>(state_in["nu_s2_eta"]) : 1.0;
  double nu_sigma2= state_in.containsElementNamed("nu_sigma2") ? Rcpp::as<double>(state_in["nu_sigma2"]) : 1.0;

  arma::mat I = arma::eye<arma::mat>(N, N);
  arma::mat A = W_use + w_use * alpha_hat_scaled.t();

  // ---- 共通の μ_t 構築子：mu_t = X_t beta + Eta * Gamma_t ----
  auto build_mu_t = [&](int t) -> arma::vec {
    arma::vec mu(N, arma::fill::zeros);
    if (K > 0 && beta.n_elem == static_cast<uword>(K)) {
      for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int k = 0; k < K; ++k) s += Xc_pre(t, i, k) * beta(k);
        mu(i) += s;
      }
    }
    if (p > 0) mu += Eta * Gamma.col(t);
    return mu;
  };

  //---------- (1) Gamma | rest（FFBS, p>0） ----------
  if (p > 0) {
    // 観測方程式: (I - rho*A) Yt = X_t beta + Eta * gamma_t + eps
    // ⇒ y*_t = (I - rho*A) Yt - X_t beta = Eta * gamma_t + eps
    arma::mat H = Eta;                     // N x p
    arma::mat Q = s2_g * arma::eye<arma::mat>(p, p);
    arma::mat Rm = sigma2 * arma::eye<arma::mat>(N, N);

    std::vector<arma::vec> a(T0), m(T0);
    std::vector<arma::mat> R(T0), C(T0);

    arma::vec m_prev = arma::zeros<arma::vec>(p);
    arma::mat C_prev = (s2_g / std::max(1e-6, 1.0 - std::min(0.9999, phi_g*phi_g)))
                       * arma::eye<arma::mat>(p, p);

    // forward pass
    for (int t = 0; t < T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec Xb(N, arma::fill::zeros);
      if (K > 0 && beta.n_elem == static_cast<uword>(K)) {
        for (int i=0;i<N;++i){ double s=0.0; for(int k=0;k<K;++k) s += Xc_pre(t,i,k) * beta(k); Xb(i)=s; }
      }
      arma::vec ystar = (I - rho*A) * Yt - Xb;

      a[t] = phi_g * m_prev;
      R[t] = phi_g * C_prev * phi_g + Q;

      arma::mat Rinv = arma::inv_sympd(R[t]);
      arma::mat S = H.t()*arma::inv_sympd(Rm)*H + Rinv;
      arma::mat Sinv = arma::inv_sympd(S);
      arma::vec b = H.t()*arma::inv_sympd(Rm)*ystar + Rinv*a[t];
      m[t] = Sinv * b;
      C[t] = Sinv;

      m_prev = m[t]; C_prev = C[t];
    }

    // backward sampling
    arma::mat Gamma_new(p, T0, arma::fill::zeros);
    arma::mat L_T = arma::chol(0.5*(C[T0-1]+C[T0-1].t()), "lower");
    Gamma_new.col(T0-1) = m[T0-1] + L_T * arma::randn<arma::vec>(p);
    for (int t=T0-2; t>=0; --t) {
      arma::mat Rinv_next = arma::inv_sympd(R[t+1]);
      arma::mat Jt = C[t] * phi_g * Rinv_next;
      arma::vec mean = m[t] + Jt * (Gamma_new.col(t+1) - a[t+1]);
      arma::mat V = C[t] - Jt * R[t+1] * Jt.t();
      arma::mat LV = arma::chol(0.5*(V+V.t()), "lower");
      Gamma_new.col(t) = mean + LV * arma::randn<arma::vec>(p);
    }
    Gamma = Gamma_new;

    // phi_g | Gamma
    double den = 0.0, num = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec prev = (t==0? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1)));
      den += arma::dot(prev, prev);  num += arma::dot(prev, Gamma.col(t));
    }
    double mean_phi = (den > 0 ? num/den : 0.0);
    double var_phi  = (den > 0 ? s2_g/den : 1.0);
    double cand_phi;
    do { cand_phi = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand_phi) > 0.999);
    phi_g = cand_phi;

    // s2_g | Gamma, phi_g
    double sc_g = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec prev = (t==0? arma::zeros<arma::vec>(p) : arma::vec(Gamma.col(t-1)));
      arma::vec diff = Gamma.col(t) - phi_g * prev;
      sc_g += 0.5 * arma::dot(diff, diff);
    }
    s2_g    = rinvgamma1(0.5 + 0.5 * p * T0, sc_g + 1.0/clip1(nu_s2_g));
    nu_s2_g = rinvgamma1(1.0, 1.0/clip1(s2_g) + 1.0/100.0);
  }

  //---------- (2) Eta | rest（p>0） ----------
  if (p > 0) {
    arma::mat GtG = Gamma * Gamma.t();                // p x p
    arma::mat Domega = arma::diagmat(omega_k);        // p x p
    arma::mat Vrow = arma::inv_sympd(GtG / clip1(sigma2) + Domega / clip1(s2_eta));
    arma::mat Lrow = arma::chol(0.5*(Vrow+Vrow.t()), "lower");

    for (int i = 0; i < N; ++i) {
      arma::vec rhs(p, arma::fill::zeros);
      for (int t = 0; t < T0; ++t) {
        arma::vec Yt = Yc_data.row(t).t();
        arma::vec Xb(N, arma::fill::zeros);
        if (K > 0 && beta.n_elem == static_cast<uword>(K)) {
          for (int ii=0; ii<N; ++ii){ double s=0.0; for(int k=0;k<K;++k) s += Xc_pre(t,ii,k)*beta(k); Xb(ii)=s; }
        }
        arma::vec ystar = (I - rho*A) * Yt - Xb;
        rhs += Gamma.col(t) * ystar(i);
      }
      arma::vec mean = Vrow * (rhs / clip1(sigma2));
      arma::vec z = arma::randn<arma::vec>(p);
      Eta.row(i) = (mean + Lrow * z).t();
    }
    // s2_eta, omega_k, nu_omega_k
    double sc_eta = 0.0;
    for (int i=0;i<N;++i) {
      arma::vec ei = Eta.row(i).t();
      sc_eta += 0.5 * arma::dot(ei, Domega * ei);
    }
    s2_eta = rinvgamma1(0.5 + 0.5*p*N, sc_eta + 1.0/clip1(nu_s2_eta));
    nu_s2_eta = rinvgamma1(1.0, 1.0/clip1(s2_eta) + 1.0/100.0);
    for (int k=0;k<p;++k) {
      double rate = 1.0/clip1(nu_omega_k(k));
      for (int i=0;i<N;++i) rate += 0.5 * Eta(i,k)*Eta(i,k) / clip1(s2_eta);
      omega_k(k) = rinvgamma1(0.5*(N+1.0), rate);
      nu_omega_k(k) = rinvgamma1(1.0, 1.0 + 1.0/clip1(omega_k(k)));
    }
  }

  //---------- (3) beta | rest（K>0） ----------
  if (K > 0) {
    arma::mat XtX(K,K,arma::fill::zeros);
    arma::vec XtY(K,arma::fill::zeros);
    for (int t=0; t<T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec fac = (p>0 ? Eta * Gamma.col(t) : arma::vec(N,arma::fill::zeros));
      arma::vec lhs = (I - rho*A) * Yt - fac; // = X_t beta + eps
      arma::mat Xt(N,K,arma::fill::zeros);
      for (int i=0;i<N;++i) for (int k=0;k<K;++k) Xt(i,k) = Xc_pre(t,i,k);
      XtX += Xt.t() * Xt;
      XtY += Xt.t() * lhs;
    }
    XtX += 1e-6 * clip1(sigma2) * arma::eye<arma::mat>(K,K);
    arma::mat Vb = arma::inv_sympd(XtX);
    arma::vec mb = Vb * XtY;
    arma::mat Sb = clip1(sigma2) * Vb;
    beta = arma::mvnrnd(mb, 0.5*(Sb+Sb.t()), 1);
  }

  //---------- (4) sigma2 | rest ----------
  double ss = 0.0;
  for (int t=0; t<T0; ++t) {
    arma::vec Yt = Yc_data.row(t).t();
    arma::vec mu = build_mu_t(t);                   // Xb + factor
    arma::vec eps = (I - rho*A) * Yt - mu;
    ss += arma::dot(eps, eps);
  }
  sigma2 = rinvgamma1(a0 + 0.5*(T0*N), b0 + 0.5*ss);
  double nu_sigma2_new = rinvgamma1(1.0, 1.0/clip1(sigma2) + 1.0/100.0);
  (void)nu_sigma2_new; // 必要なら state に戻す

  //---------- (5) rho | rest（RWMH；安定域は行列式で吸収） ----------
  auto logpost_rho = [&](double r) -> double {
    arma::mat M = I - r * A;
    double sign = 0.0, ldet = 0.0;
    bool ok = arma::log_det(ldet, sign, M);
    if (!ok || !std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    double ssum = 0.0;
    for (int t=0; t<T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec mu = build_mu_t(t);
      arma::vec res = M * Yt - mu;
      ssum += arma::dot(res, res);
    }
    return T0 * ldet - 0.5 * (T0 * N) * std::log(clip1(sigma2)) - 0.5 * ssum / clip1(sigma2);
  };
  double rho_prop = rho + step_rho * R::rnorm(0.0, 1.0);
  double lp_cur   = logpost_rho(rho);
  double lp_prop  = logpost_rho(rho_prop);
  bool moved_rho = false;
  if (std::log(R::runif(0.0,1.0)) < (lp_prop - lp_cur)) { rho = rho_prop; moved_rho = true; }

  //---------- return ----------
  return Rcpp::List::create(
    _["rho"]        = rho,
    _["sigma2"]     = sigma2,
    _["beta"]       = beta,
    _["Eta"]        = Eta,
    _["Gamma"]      = Gamma,
    _["phi_g"]      = phi_g,
    _["s2_g"]       = s2_g,
    _["nu_s2_g"]    = nu_s2_g,
    _["omega_k"]    = omega_k,
    _["nu_omega_k"] = nu_omega_k,
    _["s2_eta"]     = s2_eta,
    _["nu_s2_eta"]  = nu_s2_eta,
    _["nu_sigma2"]  = nu_sigma2,
    _["moved_rho"]  = moved_rho
  );
}