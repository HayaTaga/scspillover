# 40_geweke_full.R

# --- 正規化（w と W を行同時に） ---
normalize_joint_wW <- function(W, w) {
  stopifnot(is.matrix(W))
  stopifnot(length(w) == nrow(W), ncol(W) == nrow(W))
  N <- nrow(W)
  joint <- cbind(w, W)
  rs <- rowSums(joint)
  rs[rs == 0] <- 1
  joint_norm <- sweep(joint, 1, rs, "/")
  list(
    W_new = joint_norm[, -1, drop = FALSE],
    w_new = joint_norm[, 1, drop = TRUE]
  )
}

# --- g の既定 ---
default_g_fn <- function(theta, Yc, Y0_pre, W_use, w_use) {
  yc <- as.numeric(Yc)
  N <- ncol(Yc)
  T0 <- nrow(Yc)
  wyc <- as.numeric(Yc %*% as.numeric(w_use))
  spatial_q <- sum(diag(Yc %*% W_use %*% t(Yc))) / (N * T0)
  corr <- if (sd(Y0_pre) > 0 && sd(wyc) > 0) cor(Y0_pre, wyc) else NA_real_
  c(
    rho = unname(theta$rho),
    log_sigma2 = log(pmax(theta$sigma2, 1e-12)),
    yc_mean = mean(yc),
    log_yc_var = log(pmax(stats::var(yc), 1e-12)),
    spatial_quadratic = spatial_q,
    corr_y0_wyc = corr
  )
}

# --- 事前に整合する初期状態 ---
rinvgamma1 <- function(a, b) 1 / rgamma(1, shape = a, rate = b)

draw_initial_state <- function(T0, N, K, p, a0, b0) {
  sigma2 <- rinvgamma1(a0, b0)
  beta <- if (K > 0) rnorm(K, 0, sqrt(sigma2)) else numeric()
  Eta <- if (p > 0) matrix(rnorm(N * p, 0, 1), N, p) else matrix(0, N, 0)
  Gamma <- if (p > 0) matrix(rnorm(p * T0, 0, 1), p, T0) else matrix(0, 0, T0)
  list(
    rho = runif(1, -0.5, 0.5), # 初期値（広く取りすぎない）
    sigma2 = sigma2,
    beta = beta,
    Eta = Eta,
    Gamma = Gamma,
    phi_g = 0,
    s2_g = 1,
    nu_s2_g = 1,
    omega_k = if (p > 0) rep(1, p) else numeric(),
    nu_omega_k = if (p > 0) rep(1, p) else numeric(),
    s2_eta = 1,
    nu_s2_eta = 1,
    nu_sigma2 = 1
  )
}

# --- バッチ平均 SE ---
var_mcmc_batchmeans <- function(x, b = NULL) {
  x <- x[is.finite(x)]
  M <- length(x)
  if (is.null(b)) {
    b <- max(2L, floor(sqrt(M)))
  }
  a <- floor(M / b)
  if (a < 2L) {
    return(list(var_mean = stats::var(x) / M))
  }
  x2 <- x[seq_len(a * b)]
  bm <- colMeans(matrix(x2, nrow = b, ncol = a))
  tau2 <- b * stats::var(bm)
  list(var_mean = as.numeric(tau2) / (a * b))
}

# --- 主関数：Geweke の JDT（完全版） ---
geweke_jdt_full <- function(
  Y0_pre, # 長さ T0
  Yc_dims, # c(T0, N)
  W_raw,
  w_raw, # 空間重み（正規化前）
  alpha_hat_scaled, # 長さ N（BSCM 推定の α を列SDで割ったスケール）
  Xc_pre = NULL, # T0 x N x K（なければ NULL）
  p = 0L, # 因子の次元
  M1 = 20000L, # MC 側サンプル
  M2 = 20000L, # SC 側サンプル（keep）
  burn_in = 5000L, # SC 側バーンイン
  a0 = 1.0,
  b0 = 1.0, # Inv-Gamma 事前
  step_rho = 0.05, # RW 提案幅
  g_fn = default_g_fn,
  batch_size = NULL,
  verbose = TRUE
) {
  stopifnot(is.numeric(Y0_pre), length(Yc_dims) == 2L)
  T0 <- Yc_dims[1]
  N <- Yc_dims[2]
  K <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # 正規化（行同時）
  normed <- normalize_joint_wW(W_raw, w_raw)
  W_use <- normed$W_new
  w_use <- normed$w_new

  # X を arma::cube と互換な配列に保証
  if (is.null(Xc_pre)) {
    Xc_pre <- array(0, c(T0, N, 0L))
  } else {
    stopifnot(dim(Xc_pre)[1] == T0, dim(Xc_pre)[2] == N)
  }

  # ---------- MC 側 ----------
  if (verbose) {
    message("[JDT] MC side (iid) ...")
  }
  g_iid <- NULL
  for (m in seq_len(M1)) {
    st <- draw_initial_state(T0, N, K, p, a0, b0)
    Yc_draw <- simulate_Yc_forward_cpp(
      T0 = T0,
      W_use = W_use,
      w_use = w_use,
      alpha_hat_scaled = alpha_hat_scaled,
      rho = st$rho,
      sigma2 = st$sigma2,
      Xc_pre = Xc_pre,
      beta = if (K > 0) st$beta else numeric(0),
      Eta = if (p > 0) st$Eta else matrix(0, N, 0),
      Gamma = if (p > 0) st$Gamma else matrix(0, 0, T0)
    )
    gval <- g_fn(
      list(rho = st$rho, sigma2 = st$sigma2),
      Yc_draw,
      Y0_pre,
      W_use,
      w_use
    )
    if (is.null(g_iid)) {
      g_iid <- matrix(NA_real_, nrow = M1, ncol = length(gval))
      colnames(g_iid) <- names(gval)
    }
    g_iid[m, ] <- gval
    if (verbose && (m %% max(1L, M1 %/% 5L) == 0L)) {
      message("  ... ", m, "/", M1)
    }
  }

  # ---------- SC 側 ----------
  if (verbose) {
    message("[JDT] SC side (successive-conditional) ...")
  }
  state <- draw_initial_state(T0, N, K, p, a0, b0)

  # burn-in
  for (m in seq_len(burn_in)) {
    # ここでは y|theta を前進シミュレータで再生成してもよいが、
    # 「JDTの定義」に忠実に、y|theta -> theta'|y とする
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_pre,
      if (K > 0) state$beta else numeric(0),
      if (p > 0) state$Eta else matrix(0, N, 0),
      if (p > 0) state$Gamma else matrix(0, 0, T0)
    )
    state <- scspill_one_step_cpp(
      Yc_data = Yc_draw,
      W_use = W_use,
      w_use = w_use,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0,
      N = N,
      Xc_pre = Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho
    )
    if (verbose && (m %% max(1L, burn_in %/% 5L) == 0L)) {
      message("  ... burn ", m, "/", burn_in)
    }
  }

  # keep
  g_mcmc <- matrix(NA_real_, nrow = M2, ncol = ncol(g_iid))
  colnames(g_mcmc) <- colnames(g_iid)
  moved <- logical(M2)

  for (m in seq_len(M2)) {
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_pre,
      if (K > 0) state$beta else numeric(0),
      if (p > 0) state$Eta else matrix(0, N, 0),
      if (p > 0) state$Gamma else matrix(0, 0, T0)
    )
    state2 <- scspill_one_step_cpp(
      Yc_data = Yc_draw,
      W_use = W_use,
      w_use = w_use,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0,
      N = N,
      Xc_pre = Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho
    )
    moved[m] <- isTRUE(state2$moved_rho)
    state <- state2

    g_mcmc[m, ] <- g_fn(
      list(rho = state$rho, sigma2 = state$sigma2),
      Yc_draw,
      Y0_pre,
      W_use,
      w_use
    )
    if (verbose && (m %% max(1L, M2 %/% 5L) == 0L)) {
      message("  ... keep ", m, "/", M2)
    }
  }

  # ---------- 統計量 ----------
  mean_iid <- colMeans(g_iid, na.rm = TRUE)
  mean_mcmc <- colMeans(g_mcmc, na.rm = TRUE)
  n_iid <- colSums(is.finite(g_iid))
  se_iid <- sqrt(apply(g_iid, 2, stats::var, na.rm = TRUE) / pmax(1L, n_iid))

  if (is.null(batch_size)) {
    batch_size <- max(2L, floor(sqrt(M2)))
  }
  se_mcmc <- vapply(
    seq_len(ncol(g_mcmc)),
    function(j) sqrt(var_mcmc_batchmeans(g_mcmc[, j], batch_size)$var_mean),
    numeric(1L)
  )

  Z <- (mean_iid - mean_mcmc) / sqrt(se_iid^2 + se_mcmc^2)
  pval <- 2 * pnorm(-abs(Z))

  list(
    summary = data.frame(
      g = colnames(g_iid),
      mean_iid = as.numeric(mean_iid),
      mean_mcmc = as.numeric(mean_mcmc),
      se_iid = as.numeric(se_iid),
      se_mcmc = as.numeric(se_mcmc),
      Z = as.numeric(Z),
      pval = as.numeric(pval),
      row.names = NULL,
      check.names = FALSE
    ),
    details = list(
      M1 = M1,
      M2 = M2,
      burn_in = burn_in,
      batch_size = batch_size,
      acc_rate_rho_proxy = mean(moved)
    )
  )
}
