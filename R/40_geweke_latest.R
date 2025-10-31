# 40_geweke_full.R

# --- 前処理ユーティリティ ---
#' @keywords internal
normalize_joint_wW <- function(W, w) {
  stopifnot(is.matrix(W))
  stopifnot(length(w) == nrow(W), ncol(W) == nrow(W))
  N <- nrow(W)
  joint <- cbind(w, W)
  rs <- rowSums(joint)
  rs[rs == 0] <- 1
  joint_norm <- joint / rs
  w_new <- joint_norm[, 1, drop = TRUE]
  W_new <- joint_norm[, -1, drop = FALSE]
  list(W_new = W_new, w_new = w_new)
}

#' @keywords internal
compute_bnd <- function(W_use, w_use, alpha_hat_scaled, c_stability = 0.95) {
  A <- W_use + w_use %*% t(alpha_hat_scaled)
  ev <- eigen(A, symmetric = FALSE, only.values = TRUE)$values
  maxabs <- max(Mod(ev))
  if (!is.finite(maxabs) || maxabs < 1e-12) {
    maxabs <- 1e-12
  }
  c_stability / maxabs
}

# --- g(θ, y) ---
#' @keywords internal
default_g_fn <- function(theta, Yc, Y0_pre, W_use, w_use) {
  yc_vec <- as.numeric(Yc)
  N <- ncol(Yc)
  T0 <- nrow(Yc)
  wyc <- as.numeric(Yc %*% as.numeric(w_use))
  spatial_q <- sum(diag(Yc %*% W_use %*% t(Yc))) / (N * T0)
  corr_y0_wyc <- if (stats::sd(Y0_pre) > 0 && stats::sd(wyc) > 0) {
    cor(Y0_pre, wyc)
  } else {
    NA_real_
  }

  # パラメータ側の統計量
  beta_mean <- if (length(theta$beta) > 0) mean(theta$beta) else NA_real_
  Eta_mean <- if (length(theta$Eta) > 0) mean(theta$Eta) else NA_real_
  Gamma_mean <- if (length(theta$Gamma) > 0) mean(theta$Gamma) else NA_real_

  c(
    # 既存
    rho = unname(theta$rho),
    log_sigma2 = log(pmax(theta$sigma2, 1e-12)),

    # データ側の統計量
    yc_mean = mean(yc_vec),
    log_yc_var = log(pmax(stats::var(yc_vec), 1e-12)),
    spatial_quadratic = spatial_q,
    corr_y0_wyc = corr_y0_wyc,

    # パラメータ側の統計量（追加）
    beta_mean = beta_mean,
    Eta_mean = Eta_mean,
    Gamma_mean = Gamma_mean
  )
}

# --- 事前からの初期値（MC/SC 共通の事前と整合） ---
#' @keywords internal
draw_initial_state <- function(
  T0,
  N,
  K,
  p,
  a0,
  b0,
  W_use,
  w_use,
  alpha_hat_scaled
) {
  bnd <- compute_bnd(W_use, w_use, alpha_hat_scaled)
  rho0 <- stats::runif(1, -bnd, bnd)
  sigma2_0 <- 1 / stats::rgamma(1, shape = a0, rate = b0)
  beta0 <- if (K > 0) stats::rnorm(K, 0, 1) else numeric()
  Eta0 <- if (p > 0) {
    matrix(stats::rnorm(N * p, 0, 1), N, p)
  } else {
    matrix(0, N, 0)
  }
  Gamma0 <- if (p > 0) {
    matrix(stats::rnorm(p * T0, 0, 1), p, T0)
  } else {
    matrix(0, 0, T0)
  }
  list(
    rho = as.numeric(rho0),
    sigma2 = as.numeric(sigma2_0),
    beta = beta0,
    Eta = Eta0,
    Gamma = Gamma0
  )
}

# --- バッチ平均分散（MCMC側 SE） ---
#' @keywords internal
var_mcmc_batchmeans <- function(x, b = NULL) {
  x <- x[is.finite(x)]
  M <- length(x)
  if (is.null(b)) {
    b <- max(2L, floor(sqrt(M)))
  }
  a <- floor(M / b)
  if (a < 2L) {
    vx <- stats::var(x)
    return(list(var_mean = vx / M))
  }
  x2 <- x[seq_len(a * b)]
  bm <- colMeans(matrix(x2, nrow = b, ncol = a))
  tau2 <- b * stats::var(bm)
  list(var_mean = as.numeric(tau2) / (a * b))
}

# --- Geweke JDT 本体（完全版） ---
#' @keywords internal
geweke_jdt_full <- function(
  Y0_pre,
  Yc_pre_like_dims, # c(T0, N)
  W_raw,
  w_raw,
  alpha_hat_scaled, # N
  Xc_pre = NULL, # T0 x N x K array or NULL
  p = 0L,
  M1 = 20000L,
  M2 = 20000L,
  burn_in = 5000L,
  a0 = 1.0,
  b0 = 1.0,
  step_rho = 0.05,
  g_fn = default_g_fn,
  batch_size = NULL,
  verbose = TRUE
) {
  stopifnot(is.numeric(Y0_pre))
  T0 <- length(Y0_pre)
  N <- Yc_pre_like_dims[2]
  K <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # 正規化（行ごとに w|W を同時に）
  normed <- normalize_joint_wW(W_raw, w_raw)
  W_use <- normed$W_new
  w_use <- normed$w_new

  # ---------- MC (iid) side ----------
  if (verbose) {
    message("[JDT] MC side (iid) ...")
  }
  g_iid_mat <- NULL
  for (m in seq_len(M1)) {
    st <- draw_initial_state(
      T0,
      N,
      K,
      p,
      a0,
      b0,
      W_use,
      w_use,
      alpha_hat_scaled
    )
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      st$rho,
      st$sigma2,
      Xc_used,
      if (K > 0) st$beta else numeric(),
      if (p > 0) st$Eta else matrix(0, N, 0),
      if (p > 0) st$Gamma else matrix(0, 0, T0)
    )
    g_val <- g_fn(st, Yc_draw, Y0_pre, W_use, w_use)
    if (is.null(g_iid_mat)) {
      g_iid_mat <- matrix(NA_real_, nrow = M1, ncol = length(g_val))
      colnames(g_iid_mat) <- names(g_val)
    }
    g_iid_mat[m, ] <- g_val
  }

  # ---------- SC (successive-conditional) side ----------
  if (verbose) {
    message("[JDT] SC side (successive-conditional) ...")
  }
  state <- draw_initial_state(
    T0,
    N,
    K,
    p,
    a0,
    b0,
    W_use,
    w_use,
    alpha_hat_scaled
  )

  # burn-in
  for (m in seq_len(burn_in)) {
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_used,
      if (K > 0) state$beta else numeric(),
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
      Xc_pre = if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho
    )
  }

  # keep
  g_mcmc_mat <- matrix(NA_real_, nrow = M2, ncol = ncol(g_iid_mat))
  colnames(g_mcmc_mat) <- colnames(g_iid_mat)

  for (m in seq_len(M2)) {
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_used,
      if (K > 0) state$beta else numeric(),
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
      Xc_pre = if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho
    )
    g_mcmc_mat[m, ] <- g_fn(state, Yc_draw, Y0_pre, W_use, w_use)
  }

  # ---------- 統計量 ----------
  mean_iid <- colMeans(g_iid_mat, na.rm = TRUE)
  mean_mcmc <- colMeans(g_mcmc_mat, na.rm = TRUE)

  n_iid <- colSums(is.finite(g_iid_mat))
  se_iid <- sqrt(
    apply(g_iid_mat, 2, stats::var, na.rm = TRUE) / pmax(n_iid, 1L)
  )

  if (is.null(batch_size)) {
    batch_size <- max(2L, floor(sqrt(M2)))
  }
  se_mcmc <- vapply(
    seq_len(ncol(g_mcmc_mat)),
    function(j) sqrt(var_mcmc_batchmeans(g_mcmc_mat[, j], batch_size)$var_mean),
    numeric(1L)
  )

  Z <- (mean_iid - mean_mcmc) / sqrt(se_iid^2 + se_mcmc^2)
  pval <- 2 * stats::pnorm(-abs(Z))

  summary <- data.frame(
    g = colnames(g_iid_mat),
    mean_iid = as.numeric(mean_iid),
    mean_mcmc = as.numeric(mean_mcmc),
    se_iid = as.numeric(se_iid),
    se_mcmc = as.numeric(se_mcmc),
    Z = as.numeric(Z),
    pval = as.numeric(pval),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    summary = summary,
    details = list(M1 = M1, M2 = M2, burn_in = burn_in, batch_size = batch_size)
  )
}

# --- 最小実行例（合格確認の段階1: K=0, p=0, W=0, alpha=0） ---
#' @keywords internal
example_jdt_minimal <- function() {
  set.seed(1)
  T0 <- 10
  N <- 5
  W <- matrix(0, N, N)
  w <- rep(0, N)
  alpha_hat_scaled <- rep(0, N)
  Y0_pre <- rnorm(T0)

  out <- geweke_jdt_full(
    Y0_pre = Y0_pre,
    Yc_pre_like_dims = c(T0, N),
    W_raw = W,
    w_raw = w,
    alpha_hat_scaled = alpha_hat_scaled,
    Xc_pre = NULL,
    p = 0L,
    M1 = 2000L,
    M2 = 2000L,
    burn_in = 500L,
    a0 = 3,
    b0 = 2,
    step_rho = 0.05,
    g_fn = default_g_fn,
    verbose = TRUE
  )
  out$summary
}
