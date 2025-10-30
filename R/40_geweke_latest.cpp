// 40_geweke_full.cpp
// Geweke (2004) Joint Distribution Test – FULL version with unified priors
//
// モデル（MC/SC 共通）:
// (I - rho * A) Yc_t = X_t beta + Eta * Gamma_t + eps_t,   eps_t ~ N(0, sigma2 I_N)
// A = W + w * alpha^T,  alpha は固定（BSCM 推定に基づく外生入力）
// 事前：
//   rho ~ Unif(-b, b)   （b は A のスペクトル半径から決定）
//   sigma2 ~ InvGamma(a0, b0)
//   beta ~ N(0, I_K)
//   Gamma_t ~ N(0, I_p)  独立（t毎）
//   各行 Eta_i ~ N(0, I_p) 独立
//
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

// det(I - rho*A) の安定な log|det|（LU で十分）
// （固有値版でもよいが LU の方が実装が簡単）
inline bool logdet_I_minus_rA(const arma::mat& A, double rho, double& out) {
  arma::mat M = arma::eye<arma::mat>(A.n_rows, A.n_cols) - rho * A;
  double sign=0.0, ldet=0.0;
  bool ok = arma::log_det(ldet, sign, M);
  if (!ok || !std::isfinite(ldet)) return false;
  out = ldet;
  return true;
}

// A のスペクトル半径から一様事前サポート b を決める
inline double rho_bound_from_A(const arma::mat& A, double c_stab = 0.95) {
  arma::cx_vec ev;
  arma::eig_gen(ev, A);
  double m = 0.0;
  for (uword i=0;i<ev.n_elem;++i) m = std::max(m, std::abs(ev(i)));
  if (m < 1e-12) m = 1e-12;
  return c_stab / m;
}

//-------------------- Forward simulator --------------------
// [[Rcpp::export]]
arma::mat simulate_Yc_forward_cpp(
  int T0,
  const arma::mat& W_use,            // N x N
  const arma::vec& w_use,            // N
  const arma::vec& alpha_hat_scaled, // N
  double rho,
  double sigma2,
  const arma::cube& Xc_pre,          // T0 x N x K（K==0 可）
  const arma::vec& beta,             // 長さK（K==0なら長さ0）
  const arma::mat& Eta,              // N x p（p==0 可）
  const arma::mat& Gamma             // p x T0（p==0 可）
) {
  Rcpp::RNGScope scope;

  const int N = W_use.n_rows;
  const int K = Xc_pre.n_slices;
  const int p = Eta.n_cols;

  arma::mat Yc(T0, N, arma::fill::zeros);
  const arma::mat A = W_use + w_use * alpha_hat_scaled.t();
  const arma::mat I = arma::eye<arma::mat>(N,N);
  const double sd = std::sqrt(std::max(1e-12, sigma2));

  for (int t=0; t<T0; ++t) {
    arma::vec mu(N, arma::fill::zeros);

    // X_t beta
    if (K > 0 && (int)beta.n_elem == K) {
      for (int i=0;i<N;++i) {
        double s = 0.0;
        for (int k=0;k<K;++k) s += Xc_pre(t,i,k) * beta(k);
        mu(i) += s;
      }
    }
    // factor
    if (p > 0) mu += Eta * Gamma.col(t);

    arma::vec eps = sd * arma::randn<arma::vec>(N);
    arma::vec rhs = mu + eps;

    // (I - rho*A) Y_t = rhs
    arma::vec Yt;
    bool ok=false;
    try {
      Yt = arma::solve( I - rho*A, rhs, arma::solve_opts::fast);
      ok = Yt.is_finite();
    } catch(...) { ok=false; }
    if (!ok) {
      arma::mat M = I - rho*A;
      arma::mat AtA = M.t()*M;  AtA.diag() += 1e-10;
      arma::vec Atr = M.t()*rhs;
      Yt = arma::solve(AtA, Atr, arma::solve_opts::fast);
    }
    Yc.row(t) = Yt.t();
  }
  return Yc;
}

//-------------------- One-step posterior kernel (unified priors) --------------------
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
  Rcpp::List state_in,               // list(rho,sigma2,beta,Eta,Gamma)
  double a0,
  double b0,
  double step_rho,
  double rho_lo,
  double rho_hi
) {
  Rcpp::RNGScope scope;

  // unpack
  double rho    = as<double>(state_in["rho"]);
  double sigma2 = as<double>(state_in["sigma2"]);
  arma::vec beta = (K>0 ? as<arma::vec>(state_in["beta"]) : arma::vec());
  arma::mat Eta  = (p>0 ? as<arma::mat>(state_in["Eta"])  : arma::mat(N,0,arma::fill::zeros));
  arma::mat Gamma= (p>0 ? as<arma::mat>(state_in["Gamma"]): arma::mat(0,T0,arma::fill::zeros));

  const arma::mat A = W_use + w_use * alpha_hat_scaled.t();
  const arma::mat I = arma::eye<arma::mat>(N,N);
  const double bnd  = rho_bound_from_A(A);

  // ---- 1) Gamma | rest （iid ガウス事前） ----
  if (p>0) {
    arma::mat EtE = Eta.t() * Eta; // p x p
    arma::mat Vt  = arma::inv_sympd(EtE / clip1(sigma2) + arma::eye<arma::mat>(p,p));
    arma::mat Lt  = arma::chol(0.5*(Vt+Vt.t()), "lower");

    for (int t=0; t<T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec Xb(N, arma::fill::zeros);
      if (K>0) {
        for (int i=0;i<N;++i) {
          double s=0.0; for (int k=0;k<K;++k) s += Xc_pre(t,i,k) * beta(k);
          Xb(i)=s;
        }
      }
      arma::vec ystar = (I - rho*A) * Yt - Xb;                 // Eta * gamma_t + eps
      arma::vec mt = Vt * (Eta.t() * ystar / clip1(sigma2));
      Gamma.col(t) = mt + Lt * arma::randn<arma::vec>(p);
    }
  }

  // ---- 2) Eta | rest （iid ガウス事前） ----
  if (p>0) {
    arma::mat GtG = Gamma * Gamma.t(); // p x p
    arma::mat Vi  = arma::inv_sympd(GtG / clip1(sigma2) + arma::eye<arma::mat>(p,p));
    arma::mat Li  = arma::chol(0.5*(Vi+Vi.t()), "lower");

    for (int i=0; i<N; ++i) {
      arma::vec rhs(p, arma::fill::zeros);
      for (int t=0; t<T0; ++t) {
        arma::vec Yt = Yc_data.row(t).t();
        arma::vec Xb(N, arma::fill::zeros);
        if (K>0) {
          for (int ii=0; ii<N; ++ii) {
            double s=0.0; for (int k=0;k<K;++k) s += Xc_pre(t,ii,k) * beta(k);
            Xb(ii)=s;
          }
        }
        arma::vec ystar = (I - rho*A) * Yt - Xb; // = Eta * gamma_t + eps
        rhs += Gamma.col(t) * ystar(i);
      }
      arma::vec mi = Vi * (rhs / clip1(sigma2));
      Eta.row(i) = (mi + Li * arma::randn<arma::vec>(p)).t();
    }
  }

  // ---- 3) beta | rest （N(0,I) 事前） ----
  if (K > 0) {
    arma::mat XtX_sum(K, K, arma::fill::zeros);
    arma::vec XtY_sum(K, arma::fill::zeros);

    for (int t = 0; t < T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec fac = (p > 0 ? Eta * Gamma.col(t) : arma::vec(N, arma::fill::zeros));
      arma::vec lhs = (I - rho * A) * Yt - fac; // = X_t beta + eps

      arma::mat Xt(N, K, arma::fill::zeros);
      for (int i = 0; i < N; ++i) for (int k = 0; k < K; ++k) Xt(i, k) = Xc_pre(t, i, k);
      
      XtX_sum += Xt.t() * Xt;
      XtY_sum += Xt.t() * lhs;
    }

    // 正しい事後分布の計算
    // P_beta = (XtX_sum / sigma2) + I_K
    arma::mat P_beta = (XtX_sum / clip1(sigma2)) + arma::eye<arma::mat>(K, K);
    
    // V_beta = inv(P_beta)
    arma::mat V_beta = arma::inv_sympd(P_beta);
    
    // m_beta = V_beta * (XtY_sum / sigma2)
    arma::vec m_beta = V_beta * (XtY_sum / clip1(sigma2));

    // サンプリング
    beta = arma::mvnrnd(m_beta, 0.5 * (V_beta + V_beta.t()), 1);
  }

  // ---- 4) sigma2 | rest ----
  double ss = 0.0;
  for (int t=0; t<T0; ++t) {
    arma::vec Yt = Yc_data.row(t).t();
    arma::vec muX(N, arma::fill::zeros);
    if (K>0) { for (int i=0;i<N;++i){ double s=0.0; for(int k=0;k<K;++k) s+=Xc_pre(t,i,k)*beta(k); muX(i)=s; } }
    arma::vec fac = (p>0 ? Eta * Gamma.col(t) : arma::vec(N, arma::fill::zeros));
    arma::vec resid = (I - rho*A) * Yt - (muX + fac);
    ss += arma::dot(resid, resid);
  }
  sigma2 = rinvgamma1(a0 + 0.5*(T0*N), b0 + 0.5*ss);

  // ---- 5) rho | rest （Unif(-bnd,bnd) 事前） ----
  auto logpost_rho = [&](double r) {
    if (r < rho_lo || r > rho_hi) {
    return -std::numeric_limits<double>::infinity();
    }
    if (std::abs(r) >= bnd) return -std::numeric_limits<double>::infinity(); // support 外
    double ldet=0.0;
    if (!logdet_I_minus_rA(A, r, ldet)) return -std::numeric_limits<double>::infinity();

    double ssum = 0.0;
    arma::mat M = I - r*A;
    for (int t=0; t<T0; ++t) {
      arma::vec Yt = Yc_data.row(t).t();
      arma::vec mu(N, arma::fill::zeros);
      if (K>0) { for(int i=0;i<N;++i){ double s=0.0; for(int k=0;k<K;++k) s+=Xc_pre(t,i,k)*beta(k); mu(i)=s; } }
      if (p>0) mu += Eta * Gamma.col(t);
      arma::vec resid = M * Yt - mu;
      ssum += arma::dot(resid, resid);
    }
    double ll = T0 * ldet - 0.5*(T0*N)*std::log(clip1(sigma2)) - 0.5*ssum/clip1(sigma2);
    // 一様事前は support 内で定数 → 加算不要
    return ll;
  };

  double rho_prop = rho + step_rho * R::rnorm(0.0, 1.0);
  double lcur = logpost_rho(rho);
  double lprp = logpost_rho(rho_prop);
  bool moved = false;
  if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) { rho = rho_prop; moved = true; }

  // return
  return Rcpp::List::create(
    _["rho"]    = rho,
    _["sigma2"] = sigma2,
    _["beta"]   = beta,
    _["Eta"]    = Eta,
    _["Gamma"]  = Gamma,
    _["moved_rho"] = moved
  );
}