# =========================================================
# Utilities
# =========================================================
`%||%` <- function(x, y) if (!is.null(x)) x else y

.q025q975 <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE, type = 8)

# rook 型の隣接行列（行標準化）
rook_W <- function(nrow, ncol, normalize = TRUE) {
  N <- nrow * ncol
  nb <- matrix(0, N, N)
  id <- function(r, c) (r - 1L) * ncol + c
  for (r in 1:nrow) for (c in 1:ncol) {
    i <- id(r, c)
    if (r > 1)    nb[i, id(r - 1, c)] <- 1
    if (r < nrow) nb[i, id(r + 1, c)] <- 1
    if (c > 1)    nb[i, id(r, c - 1)] <- 1
    if (c < ncol) nb[i, id(r, c + 1)] <- 1
  }
  if (normalize) {
    rs <- rowSums(nb); rs[rs == 0] <- 1
    nb <- nb / rs
  }
  nb
}
make_w <- function(N, treated = 1L) { w <- numeric(N); w[treated] <- 1; w }

# =========================================================
# SCM: convex QP solver
# =========================================================
if (!requireNamespace("quadprog", quietly = TRUE)) install.packages("quadprog")

.scM_qp <- function(y0_pre, Yc_pre, ridge = 1e-8) {
  N <- ncol(Yc_pre)
  Dmat <- crossprod(Yc_pre) + diag(ridge, N)
  dvec <- crossprod(Yc_pre, y0_pre)
  Amat <- cbind(rep(1, N), diag(N))  # sum(a)=1 (eq), a>=0 (ineq)
  bvec <- c(1, rep(0, N))
  sol  <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  as.numeric(sol$solution)
}

# =========================================================
# Counterfactual under SCSPILL structure (post)
# =========================================================
.scspill_cf_post <- function(Y0_post, Yc_post, W, w, alpha_hat, rho_hat) {
  N  <- length(alpha_hat)
  IN <- diag(N)
  Ainv <- solve(IN - rho_hat * (w %*% t(alpha_hat) + W))
  B    <- (IN - rho_hat * W)
  T1 <- nrow(Yc_post)
  ycf <- numeric(T1)
  for (t in seq_len(T1)) {
    tmp <- Ainv %*% (B %*% Yc_post[t, ] - rho_hat * w * Y0_post[t])
    ycf[t] <- as.numeric(crossprod(alpha_hat, tmp))
  }
  ycf
}

# =========================================================
# DGP (指示の4ステップで素直に生成)
#   1) Yc^0_t = (I - ρW - wαᵀ)^{-1}(X_tβ + ε_t)
#   2) Y0^0_t = αᵀ Yc^0_t     （perfect fit）
#   3) Y0^1_t = Y0^0_t + τ_t,  τ_t ~ N(μ_τ, σ_τ^2)
#   4) Yc^1_t = (I - ρW)^{-1}(w Y0^1_t + X_tβ + ε_t)
# =========================================================
scspill_sim_dgp <- function(
  T0, T1, N, W, w,
  rho, sigma2,
  alpha,                 # 長さN（合成重みの真値）
  K = 0, beta = NULL,    # K=0ならbetaはNULLでOK
  seed = NULL,
  mu_tau = 1.0, sd_tau = 1.0
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.matrix(W), nrow(W) == N, ncol(W) == N)
  stopifnot(length(w) == N, length(alpha) == N)
  if (K > 0 && is.null(beta)) stop("K>0 なら beta を与えてください。")
  if (is.null(beta)) beta <- numeric(0)

  IN <- diag(N)

  # (A) no-treatment world 全期間
  A_pre <- IN - rho * W - rho * (w %*% t(alpha))
  rcA <- tryCatch(rcond(A_pre), error = function(e) NA_real_)
  if (!is.finite(rcA) || rcA < 1e-10) stop("A_pre is near singular.")

  A_pre_inv <- solve(A_pre)

  X_pre  <- if (K > 0) array(rnorm(T0 * N * K), dim = c(T0, N, K)) else NULL
  X_post <- if (K > 0) array(rnorm(T1 * N * K), dim = c(T1, N, K)) else NULL

  TT <- T0 + T1
  Yc0_all <- matrix(NA_real_, TT, N)
  Y00_all <- numeric(TT)

  gen_yc0 <- function(Xt) {
    rhs <- if (K > 0) as.numeric(Xt %*% beta) else rep(0, N)
    rhs <- rhs + rnorm(N, 0, sqrt(max(sigma2, 1e-12)))
    as.numeric(A_pre_inv %*% rhs)
  }
  for (t in 1:TT) {
    Xt <- if (t <= T0) {
      if (K > 0) matrix(X_pre[t, , , drop = FALSE],  N, K) else NULL
    } else {
      if (K > 0) matrix(X_post[t - T0, , , drop = FALSE], N, K) else NULL
    }
    yc0 <- gen_yc0(Xt)
    Yc0_all[t, ] <- yc0
    Y00_all[t]   <- sum(alpha %*% yc0)
  }

  # (B) post: treatment shock & SAR
  A_post <- IN - rho * W
  rcB <- tryCatch(rcond(A_post), error = function(e) NA_real_)
  if (!is.finite(rcB) || rcB < 1e-10) stop("A_post is near singular.")
  A_post_inv <- solve(A_post)

  tau_post  <- rnorm(T1, mean = mu_tau, sd = sd_tau)
  Y01_post  <- Y00_all[(T0 + 1):TT] + tau_post

  Yc1_post <- matrix(NA_real_, T1, N)
  for (tt in 1:T1) {
    Xt <- if (K > 0) matrix(X_post[tt, , , drop = FALSE], N, K) else NULL
    rhs <- rho * w * Y01_post[tt]
    if (K > 0) rhs <- rhs + as.numeric(Xt %*% beta)
    rhs <- rhs + rnorm(N, 0, sqrt(max(sigma2, 1e-12)))
    Yc1_post[tt, ] <- as.numeric(A_post_inv %*% rhs)
  }

  # 出力
  Yc_pre  <- Yc0_all[1:T0, , drop = FALSE]
  Y0_pre  <- Y00_all[1:T0]
  Yc_post <- Yc1_post
  Y0_post <- Y01_post

  colnames(Yc_pre) <- colnames(Yc_post) <- paste0("u", seq_len(N))
  names(alpha) <- colnames(Yc_pre)
  if (K > 0 && is.null(names(beta))) names(beta) <- paste0("beta_", seq_len(K))

  vec_Xc_pre  <- if (K > 0) as.numeric(aperm(X_pre,  c(1, 2, 3))) else NULL
  vec_Xc_post <- if (K > 0) as.numeric(aperm(X_post, c(1, 2, 3))) else NULL

  list(
    data = list(
      Y0_pre  = Y0_pre,   Y0_post = Y0_post,
      Yc_pre  = Yc_pre,   Yc_post = Yc_post,
      Xc_pre  = vec_Xc_pre,
      Xc_post = vec_Xc_post
    ),
    truth = list(
      rho = rho, sigma2 = sigma2,
      alpha = alpha, beta = beta,
      # tau_post = tau_post,
      tau_post = Y0_post - Y00_all[(T0 + 1):TT],
      y0_cf_post = Y00_all[(T0 + 1):TT]
    ),
    W = W, w = w,
    dims = list(T0 = T0, T1 = T1, N = N, K = K)
  )
}

# =========================================================
# 1回分のシミュレーション（SCM / BSCM / SCSPILL）
#   - dgp が NULL の場合は dgp_args から自動生成：
#       * grid=c(nrow,ncol) で rook W
#       * treated で w
#       * seed は run_one_sim の引数から dgp に伝搬
# =========================================================
run_one_sim <- function(
  dgp = NULL, dgp_args = NULL,
  M = 2000, burn = 1000,
  step_rho = 0.02,
  seed = NULL
) {
  if (is.null(dgp)) {
    if (is.null(dgp_args)) stop("Provide either `dgp` or `dgp_args`.")
    args <- as.list(dgp_args)

    # W / w を補完
    if (is.null(args$W)) {
      grid <- args$grid %||% stop("dgp_args: specify `W` or `grid = c(nrow, ncol)`.")
      args$W <- rook_W(grid[1], grid[2], normalize = TRUE)
      args$N <- nrow(args$W)
    } else {
      args$N <- nrow(args$W)
    }
    args$w     <- args$w %||% make_w(args$N, treated = args$treated %||% 1L)
    args$K     <- args$K %||% 0L
    if (args$K <= 0) args$beta <- NULL
    args$seed  <- seed

    dgp <- do.call(scspill_sim_dgp, args)
  }

  T0 <- dgp$dims$T0; T1 <- dgp$dims$T1; N <- dgp$dims$N; K <- dgp$dims$K
  W <- dgp$W; w <- dgp$w
  Y0_pre  <- dgp$data$Y0_pre;  Yc_pre  <- dgp$data$Yc_pre
  Y0_post <- dgp$data$Y0_post; Yc_post <- dgp$data$Yc_post

  # 真の TE / ATE
  y0_cf_true <- dgp$truth$y0_cf_post
  te_true    <- Y0_post - y0_cf_true
  ate_true   <- mean(te_true)

  # --- SCM（点）---
  alpha_scm <- .scM_qp(Y0_pre, Yc_pre)
  ycf_scm   <- as.numeric(Yc_post %*% alpha_scm)
  te_scm    <- Y0_post - ycf_scm
  ate_scm   <- mean(te_scm)

  # --- BSCM（αの事後）---
  alpha_draws <- hs_alpha_gibbs_cpp(
    Y0_pre = Y0_pre,
    control_outcome_pre = Yc_pre,
    iteration = M, burn = burn, verbose = FALSE
  )
  colnames(alpha_draws) <- colnames(Yc_pre)
  te_bscm_mat     <- matrix(Y0_post, nrow = T1, ncol = nrow(alpha_draws)) - (Yc_post %*% t(alpha_draws))
  te_bscm_mean    <- rowMeans(te_bscm_mat)
  ate_bscm_draws  <- colMeans(te_bscm_mat)
  ate_bscm_mean   <- mean(te_bscm_mean)
  ci_ate_bscm     <- .q025q975(ate_bscm_draws)
  cover_ate_bscm  <- as.numeric(ci_ate_bscm[1] <= ate_true && ate_true <= ci_ate_bscm[2])
  cover_pt_bscm   <- mean(vapply(seq_len(T1), function(t) {
    ci <- .q025q975(te_bscm_mat[t, ]); as.numeric(ci[1] <= te_true[t] && te_true[t] <= ci[2])
  }, numeric(1)))

  # --- SCSPILL（α & ρ を同期）---
  sar <- sar_full_sampler_cpp(
    Y0_pre = Y0_pre, Yc_pre = Yc_pre,
    Xc_pre_ = if (K > 0) dgp$data$Xc_pre else NULL,
    T0 = T0, N = N, K = K, p = 0,
    w = as.numeric(w), W = W,
    iteration = M, burn = burn,
    step_rho = step_rho, a0 = 1.0, b0 = 1.0, verbose = FALSE
  )
  rho_draws <- as.numeric(sar$rho)

  M_pair <- min(nrow(alpha_draws), length(rho_draws))
  te_spill_mat <- matrix(NA_real_, nrow = T1, ncol = M_pair)
  for (m in seq_len(M_pair)) {
    ycf_m <- .scspill_cf_post(Y0_post, Yc_post, W, w, alpha_draws[m, ], rho_draws[m])
    te_spill_mat[, m] <- Y0_post - ycf_m
  }
  te_spill_mean    <- rowMeans(te_spill_mat)
  ate_spill_draws  <- colMeans(te_spill_mat)
  ate_spill_mean   <- mean(te_spill_mean)
  ci_ate_spill     <- .q025q975(ate_spill_draws)
  cover_ate_spill  <- as.numeric(ci_ate_spill[1] <= ate_true && ate_true <= ci_ate_spill[2])
  cover_pt_spill   <- mean(vapply(seq_len(T1), function(t) {
    ci <- .q025q975(te_spill_mat[t, ]); as.numeric(ci[1] <= te_true[t] && te_true[t] <= ci[2])
  }, numeric(1)))

  # --- metrics（SCM の coverage は NA）---
  effect_metrics <- function(te_hat, te_true) {
    c(
      bias_ate   = mean(te_hat) - mean(te_true),
      rmse_ate   = sqrt(mean((mean(te_hat) - mean(te_true))^2)),
      bias_point = mean(te_hat - te_true),
      rmse_point = sqrt(mean((te_hat - te_true)^2))
    )
  }
  metrics <- rbind(
    SCM     = c(effect_metrics(te_scm,        te_true), cover95_ate = NA_real_,        cover95_point = NA_real_),
    BSCM    = c(effect_metrics(te_bscm_mean,  te_true), cover95_ate = cover_ate_bscm,  cover95_point = cover_pt_bscm),
    SCSPILL = c(effect_metrics(te_spill_mean, te_true), cover95_ate = cover_ate_spill, cover95_point = cover_pt_spill)
  )
  metrics <- as.data.frame(metrics)
  metrics$method <- rownames(metrics); rownames(metrics) <- NULL

  list(
    truth = dgp$truth,
    draws = list(
      alpha  = alpha_draws,
      rho    = rho_draws,
      ate    = list(bscm = ate_bscm_draws, scspill = ate_spill_draws)
    ),
    effects = list(
      true    = te_true,
      scm     = te_scm,
      bscm    = te_bscm_mean,
      scspill = te_spill_mean
    ),
    metrics = metrics
  )
}

# =========================================================
# Monte Carlo
# =========================================================
run_many_sim <- function(
  n_sims,
  dgp_args,            # scspill_sim_dgp に渡す引数
  seeds = NULL,
  ...                  # run_one_sim の引数（M, burn, step_rho など）
) {
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, n_sims)
  } else {
    stopifnot(length(seeds) == n_sims)
  }
  res <- vector("list", n_sims)
  for (i in seq_len(n_sims)) {
    res[[i]] <- run_one_sim(dgp = NULL, dgp_args = dgp_args, seed = seeds[i], ...)
  }
  res
}

summarize_many <- function(results) {
  stopifnot(is.list(results), length(results) > 0)
  tab <- do.call(rbind, lapply(results, function(r) r$metrics))
  keep <- c("bias_ate","rmse_ate","bias_point","rmse_point","cover95_ate","cover95_point")
  agg_mean <- aggregate(. ~ method, data = tab, FUN = mean, na.rm = TRUE)
  agg_sd   <- aggregate(. ~ method, data = tab, FUN = sd,   na.rm = TRUE)
  out <- merge(
    agg_mean[, c("method", keep)],
    setNames(agg_sd[, c("method", keep)], c("method", paste0(keep, "_sd"))),
    by = "method", sort = FALSE
  )
  out[order(match(out$method, c("SCSPILL","BSCM","SCM"))), ]
}