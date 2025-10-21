// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ========= 共通ユーティリティ（20_mcmc.cpp と整合） =========
inline double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1.0/scale); }
inline double clip1(double z, double lo=1e-12, double hi=1e12){ return std::min(std::max(z,lo),hi); }
inline double logdet_signed(const arma::mat& A){
  double sign=0.0,val=0.0; arma::log_det(val,sign,A);
  if(sign<=0.0 || !arma::is_finite(val)) return -std::numeric_limits<double>::infinity();
  return val;
}

// N×K: Xt, K×1: beta → N×1
inline arma::vec Xbeta(const arma::mat& Xt, const arma::vec& beta) {
  if (Xt.n_cols == 0) return arma::zeros<arma::vec>(Xt.n_rows);
  return Xt * beta;
}
// N×p: Eta, p×1: g → N×1
inline arma::vec Lf(const arma::mat& Eta, const arma::vec& g) {
  if (Eta.n_cols == 0) return arma::zeros<arma::vec>(Eta.n_rows);
  return Eta * g;
}

// 行和1の確認などは省略（テスト側で実施）

// ============== 1) 周辺生成器：prior 抽出（θ ~ prior） ==============
// [[Rcpp::export]]
Rcpp::List scspill_gir_prior_draw_cpp(const int N, const int K, const int p,
                                      const arma::mat& W, const arma::vec& w,
                                      double a0=1.0, double b0=1.0) {
  Rcpp::RNGScope scope;

  // ρ のサポート（20_mcmc と同じ）
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs=0.0; for(uword i=0;i<evals.n_elem;++i) maxabs=std::max(maxabs,std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  // --- 事前 ---
  // sigma2
  double sigma2 = rinvgamma(a0, b0);

  // beta0（horseshoe の階層に整合：variance = sigma2 * diag(sig2_b0)）
  arma::vec beta0 = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());
  arma::vec sig2_b0 = (K>0 ? arma::zeros<arma::vec>(K) : arma::vec());
  if (K>0) {
    for (int j=0;j<K;++j){
      // IG(1,1) 程度の素朴な層から draw（20_mcmc の更新式に合わせ“相対スケール”だけ合えば良い）
      double nu_j = rinvgamma(1.0, 1.0);                 // ν_j
      double tau2 = rinvgamma(1.0, 1.0);                 // τ^2
      double lam2 = rinvgamma(1.0, 1.0/clip1(nu_j));     // λ_j^2 | ν_j
      sig2_b0(j)  = clip1(lam2 * tau2);                  // 20_mcmc では τ を対角に入れていないが、prior 側は軽量化
      // β0 ~ N(0, σ² * sig2_b0)
      beta0(j)    = R::rnorm(0.0, std::sqrt(clip1(sigma2)*sig2_b0(j)));
    }
  }

  // 因子（Eta/Gamma）。20_mcmc は近似ガウス更新なので、ここでは素直に N(0,1)
  arma::mat Eta   = (p>0 ? arma::randn<arma::mat>(N,p) : arma::mat());
  arma::mat Gamma = (p>0 ? arma::randn<arma::mat>(p,1) : arma::mat()); // 初期カラムのみ（データ生成で展開可）
  double phi_g = (p>0 ? std::min(std::max(R::rnorm(0.0,0.3),-0.95),0.95) : 0.0);
  double s2_g  = 1.0;
  double s2_eta= 1.0;

  // rho ~ Unif((-bnd,bnd))
  double rho = R::runif(-bnd, bnd);

  return Rcpp::List::create(
    _["rho"]=rho, _["sigma2"]=sigma2, _["beta0"]=beta0,
    _["Eta"]=Eta, _["Gamma"]=Gamma,
    _["phi_g"]=phi_g, _["s2_g"]=s2_g,
    _["s2_eta"]=s2_eta
  );
}

// ============== 2) データ生成器：y | θ を生成 ==============
// Y0_t は外生系列として N(0,1) で生成（20_mcmc も観測として取り扱っている）
// [[Rcpp::export]]
Rcpp::List scspill_gir_data_sim_cpp(const Rcpp::List& theta, const int T0,
                                    const arma::mat& W, const arma::vec& w) {
  Rcpp::RNGScope scope;
  const int N = W.n_rows;

  double rho    = as<double>(theta["rho"]);
  double sigma2 = as<double>(theta["sigma2"]);
  arma::vec beta0 = (theta.containsElementNamed("beta0") ? as<arma::vec>(theta["beta0"]) : arma::vec());
  arma::mat Eta   = (theta.containsElementNamed("Eta")   ? as<arma::mat>(theta["Eta"])   : arma::mat());
  arma::mat Gamma0= (theta.containsElementNamed("Gamma") ? as<arma::mat>(theta["Gamma"]) : arma::mat());
  double phi_g    = theta.containsElementNamed("phi_g") ? as<double>(theta["phi_g"]) : 0.0;

  const int K = (beta0.n_elem>0 ? (int)beta0.n_elem : 0);
  const int p = Eta.n_cols;

  // X_t は GIR では外生に固定：N(0,1) iid
  auto X_get_row = [&](int /*t*/)->arma::mat{
    if (K==0) return arma::mat();
    return arma::randn<arma::mat>(N, K);
  };

  // Gamma の時系列生成（AR(1)）
  arma::mat Gamma = (p>0 ? arma::zeros<arma::mat>(p, T0) : arma::mat());
  if (p>0) {
    arma::vec gprev = (Gamma0.n_elem>0 ? Gamma0.col(0) : arma::randn<arma::vec>(p));
    for (int t=0;t<T0;++t){
      arma::vec noise = arma::randn<arma::vec>(p);
      arma::vec gt = (t==0 ? gprev : phi_g * Gamma.col(t-1)) + noise;
      Gamma.col(t) = gt;
    }
  }

  arma::vec Y0 = arma::randn<arma::vec>(T0);      // treated の事前系列 ~ N(0,1)
  arma::mat Yc(T0, N, arma::fill::none);

  arma::mat IN = arma::eye(N,N);
  arma::mat A  = IN - rho * W;

  for (int t=0;t<T0;++t){
    arma::vec rhs = rho * w * Y0[t];
    if (K>0) rhs += X_get_row(t) * beta0;
    if (p>0) rhs += Eta * Gamma.col(t);
    arma::vec eps = std::sqrt(sigma2) * arma::randn<arma::vec>(N);
    rhs += eps;
    // y_c,t = A^{-1} rhs
    Yc.row(t) = arma::solve(A, rhs).t();
  }
  return Rcpp::List::create(_["Y0_pre"]=Y0, _["Yc_pre"]=Yc);
}

// ============== 3) 逐次一歩：20_mcmc.cpp と同等の更新核 ==============
// 入力 theta, sim を更新して次の theta を返す
// [[Rcpp::export]]
Rcpp::List scspill_gir_sar_step_cpp(Rcpp::List theta, Rcpp::List sim,
                                    const arma::mat& W, const arma::vec& w,
                                    double a0=1.0, double b0=1.0,
                                    double step_rho=0.01,
                                    bool verbose=false) {
  Rcpp::RNGScope scope;

  // ----- データと寸法 -----
  arma::vec Y0 = as<arma::vec>(sim["Y0_pre"]);
  arma::mat Yc = as<arma::mat>(sim["Yc_pre"]);
  const int T0 = (int)Yc.n_rows, N=(int)Yc.n_cols;

  // θ の取り出し
  double rho    = as<double>(theta["rho"]);
  double s2     = as<double>(theta["sigma2"]);
  arma::vec beta0 = (theta.containsElementNamed("beta0") ? as<arma::vec>(theta["beta0"]) : arma::vec());
  arma::mat Eta   = (theta.containsElementNamed("Eta")   ? as<arma::mat>(theta["Eta"])   : arma::mat());
  arma::mat Gamma = (theta.containsElementNamed("Gamma") ? as<arma::mat>(theta["Gamma"]) : arma::mat());

  double phi_g = theta.containsElementNamed("phi_g") ? as<double>(theta["phi_g"]) : 0.0;
  double s2_g  = theta.containsElementNamed("s2_g")  ? as<double>(theta["s2_g"])  : 1.0;
  double s2_eta= theta.containsElementNamed("s2_eta")? as<double>(theta["s2_eta"]): 1.0;

  const int K = beta0.n_elem;
  const int p = Eta.n_cols;

  // X_t を逐次側でも外生生成（周辺生成器と対称）
  auto X_get_row = [&](int /*t*/)->arma::mat{
    if (K==0) return arma::mat();
    return arma::randn<arma::mat>(N, K);
  };

  arma::mat IN = arma::eye(N,N);

  // ----- (1) Gamma | rest（20_mcmc の近似ガウス更新に整合） -----
  if (p>0) {
    arma::mat EtE = Eta.t() * Eta;
    arma::mat Vg  = arma::inv_sympd( EtE / s2 + arma::eye(p,p) / s2_g );
    arma::mat Lg  = arma::chol(Vg, "lower");
    for (int t=0;t<T0;++t){
      arma::vec r = (IN - rho*W) * Yc.row(t).t() - rho * w * Y0[t];
      if (K>0) r -= X_get_row(t) * beta0;
      arma::vec mg = Vg * (Eta.t() * r / s2);
      arma::vec z  = arma::randn<arma::vec>(p);
      Gamma.col(t) = mg + Lg * z;
    }
    // phi_g, s2_g の更新（簡略：AR(1) の片側最尤に基づく事後）
    double den=0.0, num=0.0;
    for (int t=0;t<T0;++t){
      arma::vec gl;
      if (t == 0) {
        gl = arma::zeros<arma::vec>(p);
      } else {
        gl = arma::vec(Gamma.col(t-1));  // 明示的に vec にコピー
      }
      den += arma::dot(gl,gl);
      num += arma::dot(gl, Gamma.col(t));
    }
    double mean_phi = (den>0 ? num/den : 0.0);
    double var_phi  = (den>0 ? s2_g/den : 1.0);
    double cand; do{ cand = R::rnorm(mean_phi, std::sqrt(var_phi)); } while (std::abs(cand)>1.0);
    phi_g = cand;

    double sc=0.0;
    for (int t=0;t<T0;++t){
      arma::vec gl;
      if (t == 0) {
        gl = arma::zeros<arma::vec>(p);
      } else {
        gl = arma::vec(Gamma.col(t-1));  // 明示的に vec にコピー
      }
      arma::vec diff = Gamma.col(t) - phi_g * gl;
      sc += 0.5 * arma::dot(diff,diff);
    }
    s2_g = rinvgamma(0.5 + 0.5*p*T0, sc + 1.0/100.0);
  }

  // ----- (2) Eta | rest -----
  if (p>0){
    arma::mat GtG = Gamma * Gamma.t();    // p x p
    arma::mat Vrow= arma::inv_sympd( GtG / s2 + arma::eye(p,p) / s2_eta );
    arma::mat Lrow= arma::chol(Vrow,"lower");
    for (int i=0;i<N;++i){
      arma::vec rhs = arma::zeros<arma::vec>(p);
      for (int t=0;t<T0;++t){
        arma::vec r = (IN - rho*W) * Yc.row(t).t() - rho*w*Y0[t];
        if (K>0) r -= X_get_row(t) * beta0;
        rhs += Gamma.col(t) * r(i);
      }
      arma::vec m = Vrow * (rhs / s2);
      arma::vec z = arma::randn<arma::vec>(p);
      Eta.row(i) = (m + Lrow * z).t();
    }
    // s2_eta は簡略更新
    double sc=0.0;
    for (int i=0;i<N;++i){
      arma::vec ei = Eta.row(i).t();
      sc += arma::dot(ei,ei);
    }
    s2_eta = rinvgamma(0.5 + 0.5*p*N, 0.5*sc + 1.0/100.0);
  }

  // ----- (3) beta0 | rest（ridge + HS 近似：20_mcmc と同じ形） -----
  if (K>0){
    arma::mat Ab = arma::zeros<arma::mat>(K,K);
    arma::vec Bb = arma::zeros<arma::vec>(K);
    for (int t=0;t<T0;++t){
      arma::mat Xt = X_get_row(t);
      Ab += Xt.t()*Xt;
      arma::vec Btmp = (IN - rho*W) * Yc.row(t).t() - rho*w*Y0[t];
      if (p>0) Btmp -= Eta * Gamma.col(t);
      Bb += Xt.t() * Btmp;
    }
    // HS の局所分散を「既定の1」に固定（GIR 用の核一致が目的）
    Ab.diag() += clip1(s2);                 // 20_mcmc の Ab.diag() += s2 * (1/σ²_{β0}) と同型の簡便化
    arma::mat Ainv = arma::inv_sympd(Ab);
    arma::vec m    = Ainv * Bb;
    arma::mat S    = clip1(s2) * Ainv;
    beta0 = arma::mvnrnd(m, 0.5*(S+S.t()), 1);
  }

  // ----- (4) sigma2 | rest -----
  double ss=0.0;
  for (int t=0;t<T0;++t){
    arma::vec mu = rho*w*Y0[t];
    if (K>0) mu += X_get_row(t) * beta0;
    if (p>0) mu += Eta * Gamma.col(t);
    arma::vec u = (IN - rho*W) * Yc.row(t).t() - mu;
    ss += arma::dot(u,u);
  }
  double shape = 0.5*(T0*N) + a0;
  double rate  = 0.5*ss + b0;
  s2 = rinvgamma(shape, rate);

  // ----- (5) rho | rest（RW-MH） -----
  auto ll = [&](double r)->double{
    arma::mat Mmat = IN - r*W;
    double ldet = logdet_signed(Mmat);
    if (!std::isfinite(ldet)) return -std::numeric_limits<double>::infinity();
    double ss2=0.0;
    for (int t=0;t<T0;++t){
      arma::vec mu = r*w*Y0[t];
      if (K>0) mu += X_get_row(t) * beta0;
      if (p>0) mu += Eta * Gamma.col(t);
      arma::vec u = Mmat * Yc.row(t).t() - mu;
      ss2 += arma::dot(u,u);
    }
    return T0*ldet - 0.5*(N*T0)*std::log(s2) - 0.5*ss2/s2;
  };
  // ρ の境界
  arma::cx_vec evals = arma::eig_gen(W);
  double maxabs=0.0; for(uword i=0;i<evals.n_elem;++i) maxabs=std::max(maxabs,std::abs(evals[i]));
  double bnd = 0.95 / std::max(1.0, maxabs);

  double prop = R::rnorm(rho, step_rho);
  if (std::abs(prop) < bnd) {
    double lcur = ll(rho), lprp = ll(prop);
    if (std::log(R::runif(0.0,1.0)) < (lprp - lcur)) rho = prop;
  }

  // 返却
  return Rcpp::List::create(
    _["rho"]=rho, _["sigma2"]=s2, _["beta0"]=beta0,
    _["Eta"]=Eta, _["Gamma"]=Gamma,
    _["phi_g"]=phi_g, _["s2_g"]=s2_g, _["s2_eta"]=s2_eta
  );
}