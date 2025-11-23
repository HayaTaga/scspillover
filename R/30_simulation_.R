# =========================================================
# Utilities
# =========================================================
`%||%` <- function(x, y) if (!is.null(x)) x else y

.q025q975 <- function(x) {
  stats::quantile(x, c(0.025, 0.975), names = FALSE)
}

#' @keywords internal
# rook 型の隣接行列（行標準化）
rook_W <- function(nrow, ncol, normalize = FALSE) {
  N <- nrow * ncol
  nb <- matrix(0, N, N)
  id <- function(r, c) (r - 1L) * ncol + c
  for (r in 1:nrow) {
    for (c in 1:ncol) {
      i <- id(r, c)
      if (r > 1) {
        nb[i, id(r - 1, c)] <- 1
      }
      if (r < nrow) {
        nb[i, id(r + 1, c)] <- 1
      }
      if (c > 1) {
        nb[i, id(r, c - 1)] <- 1
      }
      if (c < ncol) nb[i, id(r, c + 1)] <- 1
    }
  }
  if (normalize) {
    rs <- rowSums(nb)
    rs[rs == 0] <- 1
    nb <- nb / rs
  }
  nb
}

#' @keywords internal
make_w <- function(N, treated = 1L) {
  w <- numeric(N)
  w[treated] <- 1
  w
}

.scspill_cf_post <- function(
    Y0_post,
    Yc_post,
    W,
    w,
    alpha_hat,
    rho_hat) {
  N <- length(alpha_hat)
  IN <- diag(N)

  Ainv <- solve(IN - rho_hat * (w %*% t(alpha_hat) + W))
  B <- (IN - rho_hat * W)
  T1 <- nrow(Yc_post)
  ycf <- numeric(T1)
  for (t in seq_len(T1)) {
    tmp <- Ainv %*% (B %*% Yc_post[t, ] - rho_hat * w * Y0_post[t])
    ycf[t] <- as.numeric(crossprod(alpha_hat, tmp))
  }
  ycf
}

# =========================================================
# SCM: convex QP solver
# =========================================================
if (!requireNamespace("quadprog", quietly = TRUE)) {
  install.packages("quadprog")
}

.scM_qp <- function(y0_pre, Yc_pre, ridge = 1e-8) {
  N <- ncol(Yc_pre)
  Dmat <- crossprod(Yc_pre) + diag(ridge, N)
  dvec <- crossprod(Yc_pre, y0_pre)
  Amat <- cbind(rep(1, N), diag(N)) # sum(a)=1 (eq), a>=0 (ineq)
  bvec <- c(1, rep(0, N))
  sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  as.numeric(sol$solution)
}

#' @keywords internal
scspill_sim_dgp <- function(
    T0,
    T1,
    N,
    W,
    w,
    rho,
    sigma2,
    alpha, # length N (true synthetic weights)
    K = 0,
    beta = NULL, # if K=0, beta can be NULL
    seed = NULL,
    mu_tau = 1.0,
    sd_tau = 1.0) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  stopifnot(is.matrix(W), nrow(W) == N, ncol(W) == N)
  stopifnot(length(w) == N, length(alpha) == N)
  if (K > 0 && is.null(beta)) {
    stop("If K>0, beta must be provided.")
  }
  if (is.null(beta)) {
    beta <- numeric(0)
  }

  IN <- diag(N)

  W_use <- row_normalize(W)
  w_use <- as.numeric(w)
  wsum <- sum(w_use)
  if (is.finite(wsum) && wsum > 1e-12) {
    w_use <- w_use / wsum
  }
  A_pre <- IN - rho * W_use - rho * (w_use %*% t(alpha))
  rcA <- tryCatch(rcond(A_pre), error = function(e) NA_real_)
  if (!is.finite(rcA) || rcA < 1e-10) {
    stop("A_pre is near singular.")
  }

  A_pre_inv <- solve(A_pre)

  X_pre <- if (K > 0) array(rnorm(T0 * N * K), dim = c(T0, N, K)) else NULL
  X_post <- if (K > 0) array(rnorm(T1 * N * K), dim = c(T1, N, K)) else NULL

  TT <- T0 + T1
  Yc0_all <- matrix(NA_real_, TT, N)
  Y00_all <- numeric(TT)

  error_all <- matrix(
    rnorm(TT * N, 0, sqrt(max(sigma2, 1e-12))),
    nrow = TT,
    ncol = N
  )

  gen_yc0 <- function(Xt, U_t) {
    rhs <- if (K > 0) as.numeric(Xt %*% beta) else rep(0, N)
    rhs <- rhs + U_t
    as.numeric(A_pre_inv %*% rhs)
  }
  for (t in 1:TT) {
    Xt <- if (t <= T0) {
      if (K > 0) matrix(X_pre[t, , , drop = FALSE], N, K) else NULL
    } else {
      if (K > 0) matrix(X_post[t - T0, , , drop = FALSE], N, K) else NULL
    }
    yc0 <- gen_yc0(Xt, error_all[t, ])
    Yc0_all[t, ] <- yc0
    Y00_all[t] <- as.numeric(crossprod(alpha, yc0))
  }

  A_post <- IN - rho * W_use
  rcB <- tryCatch(rcond(A_post), error = function(e) NA_real_)
  if (!is.finite(rcB) || rcB < 1e-10) {
    stop("A_post is near singular.")
  }
  A_post_inv <- solve(A_post)

  tau_post <- rnorm(T1, mean = mu_tau, sd = sd_tau)
  Y01_post <- Y00_all[(T0 + 1):TT] + tau_post

  Yc1_post <- matrix(NA_real_, T1, N)
  for (tt in 1:T1) {
    Xt <- if (K > 0) matrix(X_post[tt, , , drop = FALSE], N, K) else NULL
    rhs <- rho * w_use * Y01_post[tt]
    if (K > 0) {
      rhs <- rhs + as.numeric(Xt %*% beta)
    }
    rhs <- rhs + error_all[T0 + tt, ]
    Yc1_post[tt, ] <- as.numeric(A_post_inv %*% rhs)
  }

  # 出力
  Yc_pre <- Yc0_all[1:T0, , drop = FALSE]
  Y0_pre <- Y00_all[1:T0]
  Yc_post <- Yc1_post
  Y0_post <- Y01_post

  colnames(Yc_pre) <- colnames(Yc_post) <- paste0("u", seq_len(N))
  names(alpha) <- colnames(Yc_pre)
  if (K > 0 && is.null(names(beta))) {
    names(beta) <- paste0("beta_", seq_len(K))
  }

  vec_Xc_pre <- if (K > 0) as.numeric(aperm(X_pre, c(1, 2, 3))) else NULL
  vec_Xc_post <- if (K > 0) as.numeric(aperm(X_post, c(1, 2, 3))) else NULL

  # print(A_pre_inv)
  # print(A_post_inv)

  list(
    data = list(
      Y0_pre = Y0_pre,
      Y0_post = Y0_post,
      Yc_pre = Yc_pre,
      Yc_post = Yc_post,
      Xc_pre = vec_Xc_pre,
      Xc_post = vec_Xc_post
    ),
    truth = list(
      rho = rho,
      sigma2 = sigma2,
      alpha = alpha,
      beta = beta,
      # tau_post = tau_post,
      tau_post = Y0_post - Y00_all[(T0 + 1):TT],
      y0_cf_post = Y00_all[(T0 + 1):TT]
    ),
    W = W_use,
    w = w_use,
    dims = list(T0 = T0, T1 = T1, N = N, K = K)
  )
}

# =========================================================
# Single simulation run (SCM / BSCM / SCSPILL)
#   - if dgp is NULL, generate from dgp_args:
#       * grid=c(nrow,ncol) for rook W
#       * treated for w
#       * seed passed from run_one_sim arguments to dgp
# =========================================================
#' @keywords internal
run_one_sim <- function(
    dgp = NULL,
    dgp_args = NULL,
    M = 2000,
    burn = 1000,
    step_rho = 0.02, # C++ 側の引数に合わせて残置（Stan は使いません）
    seed = NULL) {
  if (is.null(dgp)) {
    if (is.null(dgp_args)) {
      stop("Provide either `dgp` or `dgp_args`.")
    }
    args <- as.list(dgp_args)

    if (is.null(args$W)) {
      grid <- args$grid %||%
        stop("dgp_args: specify `W` or `grid = c(nrow, ncol)`.")
      args$W <- rook_W(grid[1], grid[2], normalize = FALSE)
      args$N <- nrow(args$W)
    } else {
      args$N <- nrow(args$W)
    }
    args$w <- args$w %||% make_w(args$N, treated = args$treated %||% 1L)
    args$K <- args$K %||% 0L
    if (args$K <= 0) {
      args$beta <- NULL
    }
    args$seed <- seed

    dgp <- do.call(scspill_sim_dgp, args)
  }

  T0 <- dgp$dims$T0
  T1 <- dgp$dims$T1
  N <- dgp$dims$N
  K <- dgp$dims$K
  W <- dgp$W
  w <- dgp$w
  Y0_pre <- dgp$data$Y0_pre
  Yc_pre <- dgp$data$Yc_pre
  Y0_post <- dgp$data$Y0_post
  Yc_post <- dgp$data$Yc_post

  # 真の効果
  y0_cf_true <- dgp$truth$y0_cf_post
  te_true <- Y0_post - y0_cf_true
  ate_true <- mean(te_true)

  # --- SCM（点推定）---
  alpha_scm <- .scM_qp(Y0_pre, Yc_pre)
  ycf_scm <- as.numeric(Yc_post %*% alpha_scm)
  te_scm <- Y0_post - ycf_scm
  ate_scm <- mean(te_scm)

  alpha_draws_bscm <- hs_alpha_gibbs_cpp(
    Y0_pre = Y0_pre,
    control_outcome_pre = Yc_pre,
    iteration = M,
    burn = burn,
    verbose = FALSE
  )
  colnames(alpha_draws_bscm) <- colnames(Yc_pre)
  te_bscm_mat <- matrix(Y0_post, nrow = T1, ncol = nrow(alpha_draws_bscm)) -
    (Yc_post %*% t(alpha_draws_bscm))
  te_bscm_mean <- rowMeans(te_bscm_mat)
  ate_bscm_draws <- colMeans(te_bscm_mat)
  ate_bscm_mean <- mean(te_bscm_mean)
  ci_ate_bscm <- .q025q975(ate_bscm_draws)
  cover_ate_bscm <- as.numeric(
    ci_ate_bscm[1] <= ate_true && ate_true <= ci_ate_bscm[2]
  )
  cover_pt_bscm <- mean(vapply(
    seq_len(T1),
    function(t) {
      ci <- .q025q975(te_bscm_mat[t, ])
      as.numeric(ci[1] <= te_true[t] && te_true[t] <= ci[2])
    },
    numeric(1)
  ))

  alpha_hat_bscm <- colMeans(alpha_draws_bscm)

  # --- SCSPILL（C++ サンプラの α・ρ を共同で使用）---
  sar <- sar_full_sampler_cpp_step2(
    Yc_pre = Yc_pre,
    alpha_hat_in = alpha_hat_bscm, # ★ Step 1 の alpha_hat を渡す
    Xc_pre_ = if (K > 0) dgp$data$Xc_pre else NULL,
    T0 = T0,
    N = N,
    K = K,
    p = 0, # simulation assumes p=0 (no factors)
    w = as.numeric(w),
    W = W,
    iteration = M,
    burn = burn,
    step_rho = step_rho,
    a0 = 1.0,
    b0 = 1.0,
    verbose = FALSE
  )

  rho_draws_step2 <- as.numeric(sar$rho)
  M_rho <- length(rho_draws_step2)

  S <- min(nrow(alpha_draws_bscm), M_rho)
  idx_a <- sample(seq_len(nrow(alpha_draws_bscm)), S, replace = (S > nrow(alpha_draws_bscm)))
  idx_r <- sample(seq_len(M_rho), S, replace = (S > M_rho))

  te_spill_mat <- matrix(NA_real_, nrow = T1, ncol = S)

  for (s in seq_len(S)) {
    ah <- as.numeric(alpha_draws_bscm[idx_a[s], ])
    rh <- rho_draws_step2[idx_r[s]]
    ycf_m <- .scspill_cf_post(
      Y0_post = Y0_post,
      Yc_post = Yc_post,
      W = W,
      w = w,
      alpha_hat = ah,
      rho_hat = rh
    )
    te_spill_mat[, s] <- (Y0_post - ycf_m)
  }

  # 1) 各時点の事後平均（推定量）と 95% CI、被覆
  te_spill_mean <- rowMeans(te_spill_mat) # 事後平均（時点別推定量）
  te_spill_ci <- t(apply(
    te_spill_mat,
    1,
    stats::quantile,
    probs = c(0.025, 0.975)
  ))
  colnames(te_spill_ci) <- c("lower", "upper")
  cover_pt_spill_vec <- as.numeric(
    te_spill_ci[, "lower"] <= te_true &
      te_true <= te_spill_ci[, "upper"]
  )
  cover_pt_spill <- mean(cover_pt_spill_vec)
  # print(te_true)
  # print("--------")
  # print(te_spill_ci)
  # print("================")

  # 2) ATE posterior (average over draws, then CI/coverage from distribution)
  ate_spill_draws <- colMeans(te_spill_mat)
  ate_spill_mean <- mean(ate_spill_draws)
  ci_ate_spill <- stats::quantile(
    ate_spill_draws,
    c(0.025, 0.975),
    names = FALSE
  )
  cover_ate_spill <- as.numeric(
    ci_ate_spill[1] <= mean(te_true) &
      mean(te_true) <= ci_ate_spill[2]
  )

  # 3) MSE/Bias at each time point (Python: te_true - mean)
  mse_spill_time <- (te_true - te_spill_mean)^2
  bias_spill_time <- (te_true - te_spill_mean)

  effect_metrics <- function(te_hat, te_true) {
    c(
      bias_point = mean(te_true - te_hat),
      mse_point = mean((te_true - te_hat)^2)
    )
  }

  met_SCM <- effect_metrics(te_scm, te_true)
  met_BSCM <- effect_metrics(te_bscm_mean, te_true)
  met_SP <- effect_metrics(te_spill_mean, te_true)

  ate_true_scalar <- mean(te_true)
  ate_err_SCM <- ate_true_scalar - ate_scm
  ate_err_BSCM <- ate_true_scalar - ate_bscm_mean
  ate_err_SP <- ate_true_scalar - ate_spill_mean

  metrics <- rbind(
    SCM = c(
      bias_ate = ate_err_SCM,
      mse_ate = ate_err_SCM^2,
      met_SCM,
      cover95_ate = NA_real_,
      cover95_point = NA_real_
    ),
    BSCM = c(
      bias_ate = ate_err_BSCM,
      mse_ate = ate_err_BSCM^2,
      met_BSCM,
      cover95_ate = cover_ate_bscm,
      cover95_point = cover_pt_bscm
    ),
    SCSPILL = c(
      bias_ate = ate_err_SP,
      mse_ate = ate_err_SP^2,
      met_SP,
      cover95_ate = cover_ate_spill,
      cover95_point = cover_pt_spill
    )
  )
  metrics <- as.data.frame(metrics)
  metrics$method <- rownames(metrics)
  rownames(metrics) <- NULL

  per_time_mean <- list(
    true = te_true,
    scm = te_scm,
    bscm = te_bscm_mean,
    scspill = te_spill_mean
  )

  per_time_mse <- list(
    scm = (te_true - te_scm)^2,
    bscm = (te_true - te_bscm_mean)^2,
    scspill = (te_true - te_spill_mean)^2
  )

  per_time_bias <- list(
    scm = te_true - te_scm,
    bscm = te_true - te_bscm_mean,
    scspill = te_true - te_spill_mean
  )

  bscm_ci <- t(apply(te_bscm_mat, 1, stats::quantile, probs = c(0.025, 0.975)))
  spill_ci <- t(apply(
    te_spill_mat,
    1,
    stats::quantile,
    probs = c(0.025, 0.975)
  ))
  per_time_ci <- list(
    bscm = list(
      lower = bscm_ci[, 1],
      upper = bscm_ci[, 2],
      cover = as.numeric(bscm_ci[, 1] <= te_true & te_true <= bscm_ci[, 2])
    ),
    scspill = list(
      lower = spill_ci[, 1],
      upper = spill_ci[, 2],
      cover = as.numeric(spill_ci[, 1] <= te_true & te_true <= spill_ci[, 2])
    )
  )

  list(
    truth = dgp$truth,
    draws = list(
      alpha_bscm = alpha_draws_bscm,
      scspill_step2 = list(alpha_hat = alpha_hat_bscm, rho = rho_draws_step2),
      ate = list(bscm = ate_bscm_draws, scspill = ate_spill_draws),
      te_path = list(bscm = te_bscm_mat, scspill = te_spill_mat) # フル行列も保持
    ),
    per_time = list(
      mean = per_time_mean,
      mse = per_time_mse,
      bias = per_time_bias,
      ci = per_time_ci
    ),
    effects = list(
      true = te_true,
      scm = te_scm,
      bscm = te_bscm_mean,
      scspill = te_spill_mean
    ),
    metrics = metrics
  )
}

# =========================================================
# Monte Carlo
# =========================================================
#' @keywords internal
run_many_sim <- function(
    n_sims,
    dgp_args, # arguments passed to scspill_sim_dgp
    seeds = NULL,
    ... # arguments for run_one_sim (M, burn, step_rho, etc.)
    ) {
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, n_sims)
  } else {
    stopifnot(length(seeds) == n_sims)
  }
  res <- vector("list", n_sims)
  for (i in seq_len(n_sims)) {
    res[[i]] <- run_one_sim(
      dgp = NULL,
      dgp_args = dgp_args,
      seed = seeds[i],
      ...
    )
  }
  res
}

#' @keywords internal
summarize_many <- function(results) {
  stopifnot(is.list(results), length(results) > 0)
  tab <- do.call(rbind, lapply(results, function(r) r$metrics))

  keep <- c(
    "bias_ate",
    "mse_ate",
    "bias_point",
    "mse_point",
    "cover95_ate",
    "cover95_point"
  )
  lev <- c("SCM", "BSCM", "SCSPILL")
  tab$method <- factor(tab$method, levels = lev)

  safe_mean <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  safe_sd <- function(x) if (all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE)

  agg_mean <- aggregate(
    tab[, keep, drop = FALSE],
    list(method = tab$method),
    safe_mean
  )
  agg_sd <- aggregate(
    tab[, keep, drop = FALSE],
    list(method = tab$method),
    safe_sd
  )
  names(agg_sd)[-1] <- paste0(names(agg_sd)[-1], "_sd")

  out <- merge(agg_mean, agg_sd, by = "method", all = TRUE)

  out$rmse_ate <- sqrt(out$mse_ate)
  out$rmse_point <- sqrt(out$mse_point)

  out$method <- factor(out$method, levels = lev)
  out <- out[order(as.integer(out$method)), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Row-normalize a spatial weights matrix W
#'
#' Ensures that the sum of each row is 1.
#' Sets the diagonal to 0 and handles rows that sum to 0.
#'
#' @param W A numeric matrix.
#' @param tol Tolerance for checking if a row sum is zero.
#' @param zero_policy How to handle rows that sum to zero (or are close to it).
#'   "keep" (default): leaves the row as all zeros.
#'   "uniform": (not implemented here, but common) sets to 1/N.
#' @return A row-normalized matrix.
#'
row_normalize <- function(W, tol = 1e-12, zero_policy = c("keep")) {
  zero_policy <- match.arg(zero_policy)

  if (!is.matrix(W) || !is.numeric(W)) {
    stop("W must be a numeric matrix.")
  }

  # Ensure diagonal is zero (no self-loops)
  diag(W) <- 0

  # Calculate row sums
  rs <- rowSums(W, na.rm = TRUE)

  # Find rows that are not zero (or very close to it)
  nz <- rs > tol

  # Normalize non-zero rows
  if (any(nz)) {
    W[nz, ] <- W[nz, , drop = FALSE] / rs[nz]
  }

  # Handle zero-sum rows (if any)
  if (any(!nz)) {
    if (zero_policy == "keep") {
      # Do nothing, leave the row as all zeros
      W[!nz, ] <- 0 # Ensure it's clean
    }
  }

  W
}
