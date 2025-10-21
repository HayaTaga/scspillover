#' Geweke (2004) joint distribution test for scspill (20_mcmc kernels aligned)
#' @param M1 周辺生成器サンプル数（prior→data）
#' @param M2 逐次1ステップサンプル数（posterior 1-step）
#' @param dims list(T0, N, K, p)
#' @param W, w 行和1を前提（必要なら正規化して渡す）
#' @param priors list(a0=1, b0=1)  # βのHS層は逐次側で近似
#' @param verbose ログ有無
#' @return list(z=Zベクトル, pval=両側p, mc=周辺行列, sc=逐次行列)
#' @export
scspill_geweke <- function(M1=2000L, M2=2000L,
                           dims=list(T0=12L, N=6L, K=2L, p=0L),
                           W, w,
                           priors=list(a0=1, b0=1),
                           verbose=FALSE) {
  T0 <- dims$T0; N <- dims$N; K <- dims$K; p <- dims$p

  # ---- (A) 周辺系列 ----
  mc <- matrix(NA_real_, nrow=M1, ncol=0)
  for (m in seq_len(M1)) {
    th  <- scspill_gir_prior_draw_cpp(N=N, K=K, p=p, W=W, w=w,
                                      a0=priors$a0, b0=priors$b0)
    sim <- scspill_gir_data_sim_cpp(th, T0=T0, W=W, w=w)
    g   <- .scspill_gir_stats(th, sim)   # 統計量の抽出（下で定義）
    if (ncol(mc) == 0) {
      mc <- matrix(NA_real_, nrow=M1, ncol=length(g), dimnames=list(NULL, names(g)))
    }
    mc[m, ] <- g
    if (verbose && (m %% 500 == 0)) message("marginal: ", m, "/", M1)
  }

  # ---- (B) 逐次系列 ----
  th  <- scspill_gir_prior_draw_cpp(N=N, K=K, p=p, W=W, w=w,
                                    a0=priors$a0, b0=priors$b0)
  sc <- matrix(NA_real_, nrow=M2, ncol=ncol(mc), dimnames=dimnames(mc))
  for (m in seq_len(M2)) {
    sim <- scspill_gir_data_sim_cpp(th, T0=T0, W=W, w=w)
    g   <- .scspill_gir_stats(th, sim)
    sc[m, ] <- g
    th  <- scspill_gir_sar_step_cpp(th, sim, W=W, w=w,
                                    a0=priors$a0, b0=priors$b0,
                                    step_rho=0.05, verbose=FALSE)
    if (verbose && (m %% 500 == 0)) message("successive: ", m, "/", M2)
  }

  # ---- (C) Z と p 値 ----
  z  <- .gir_zstats(mc, sc)
  pv <- 2 * pnorm(-abs(z))
  list(z=z, pval=pv, mc=mc, sc=sc)
}

# テスト関数 g(θ,y)：rho, sigma2, beta の一部、二次・交差モーメント
.scspill_gir_stats <- function(theta, sim) {
  rho <- theta$rho
  s2  <- theta$sigma2
  b   <- theta$beta0
  out <- c(rho=rho, sigma2=s2)
  if (length(b)) {
    nm <- paste0("beta_", seq_along(b))
    out <- c(out, setNames(as.numeric(b), nm))
  }
  out <- c(out,
           rho2.rho = rho^2,
           rhos2.rho = rho * s2)
  out
}

# Newey–West の簡易 HAC 分散で Z を計算
.gir_zstats <- function(mc, sc, maxlag=NULL) {
  stopifnot(ncol(mc)==ncol(sc))
  G1 <- scale(mc, center=TRUE, scale=FALSE)
  G2 <- scale(sc, center=TRUE, scale=FALSE)
  if (is.null(maxlag)) {
    maxlag <- floor(4 * nrow(G1)^(1/3))
  }
  v1 <- .hac_var_diag(G1, maxlag)
  v2 <- .hac_var_diag(G2, maxlag)
  m1 <- colMeans(mc); m2 <- colMeans(sc)
  denom <- sqrt(v1 / nrow(mc) + v2 / nrow(sc))
  z <- as.numeric((m1 - m2) / denom)
  names(z) <- colnames(mc)
  z
}

# 対角だけ使う簡易 HAC（本件のモーメント比較に十分）
.hac_var_diag <- function(X, L) {
  T <- nrow(X); K <- ncol(X)
  v <- numeric(K)
  w <- function(h, L) 1 - h/(L+1)              # Bartlett
  for (k in seq_len(K)) {
    x <- X[,k]
    s0 <- sum((x - mean(x))^2) / T
    sc <- 0
    for (h in 1:min(L, T-1)) {
      g <- sum( (x[1:(T-h)]-mean(x)) * (x[(1+h):T]-mean(x)) ) / T
      sc <- sc + 2 * w(h,L) * g
    }
    v[k] <- s0 + sc
  }
  v
}