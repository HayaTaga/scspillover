rm(list = ls())
source("R/30_simulation.R")
source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover.R")
source("R/03_utils_plot.R")
source("R/40_geweke.R")
Rcpp::sourceCpp("./R/20_mcmc.cpp")
Rcpp::sourceCpp("./R/40_geweke.cpp")
# # --- グリッドと重み ---
# N <- 8
# T0 <- 20
# K <- 2
# p <- 0
# W <- matrix(0, N, N)
# for (i in 1:N) {
#   W[i, if (i < N) i + 1 else 1] <- 1
# } # 環状1近傍
# W <- W / rowSums(W)
# w <- numeric(N)
# w[2] <- 1 # 弱いスピル先

# dims <- list(T0 = T0, N = N, K = K, p = p)

# # --- 真値 ---
# theta_true <- make_theta_true(
#   rho = 0.30,
#   sigma2 = 1.0,
#   beta = c(0.5, -0.3),
#   N = N,
#   p = p
# )

# # --- 1 回走らせてチェーン診断 ---
# set.seed(1)
# one <- run_one_sim(
#   dims = dims,
#   W = W,
#   w = w,
#   theta_true = theta_true,
#   iteration = 3000,
#   burn = 1500,
#   step_rho = 0.03
# )

# one$diag$rho # mean, sd, q025, q975, ESS, MCSE, GewekeZ, ACF1, ACF5
# one$diag$sigma2
# one$diag$acc_rate # RW-MH 受容率（rho）

# # --- 繰り返し（例：R=50）で精度（バイアス/RMSE/被覆） ---
# set.seed(42)
# R <- 1000
# many <- lapply(seq_len(R), function(r) {
#   run_one_sim(
#     dims,
#     W,
#     w,
#     theta_true,
#     iteration = 5000,
#     burn = 1000,
#     step_rho = 0.03,
#     seed = 100 + r
#   )
# })

# summarize_many(many)


## 1) 真値の設定

make_rook_W <- function(nrow, ncol, normalize = TRUE) {
  stopifnot(nrow >= 1L, ncol >= 1L)
  idx <- matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE)
  N <- nrow * ncol
  W <- matrix(0, N, N)

  for (r in 1:nrow) for (c in 1:ncol) {
    i <- idx[r, c]
    if (r > 1)     W[i, idx[r - 1, c]] <- 1
    if (r < nrow)  W[i, idx[r + 1, c]] <- 1
    if (c > 1)     W[i, idx[r, c - 1]] <- 1
    if (c < ncol)  W[i, idx[r, c + 1]] <- 1
  }

  if (normalize) {
    rs <- rowSums(W)
    W <- W / pmax(rs, 1)  # 0割回避
  }
  W
}

N <- 16
W <- make_rook_W(sqrt(N), sqrt(N))
W
w <- rep(0, N); w[1:4] <- 1;
alpha_true <- rep(0, N); alpha_true[1] <- 0.5
alpha_true[2] <- -0.2
alpha_true[3:4] <- 0.4
alpha_true[5:10] <- 0.1/6

## 2) DGP で 1 セット生成
dgp_args <- list(
  T0 = 30, T1 = 20, N = N,
  W = W, w = w,
  rho = 0.3, sigma2 = 1.0,
  alpha = alpha_true,
  K=1, beta=c(1.0)
)


# 1回実行（SCMとBayesian SCMの効果推定）
out1 <- run_one_sim(dgp = NULL, dgp_args = dgp_args, M=5000, burn=2000)
out1$metrics

# 多回
many <- run_many_sim(100, dgp_args, M=5000, burn=2000,)
summarize_many(many)



dgp <- NULL
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
    args$seed  <- 1
  }


dgp <- do.call(scspill_sim_dgp, args)

dgp

dgp$data$Y0_post - dgp$data$Y0_pre

dgp$truth$tau_post


rho <- 0.3
alpha <- alpha_true
IN <- diag(N)
K <- 1
T0 <- 10
T1 <- 10
beta <- c(1.0)
sigma2 <- 1.0

  # (A) no-treatment world 全期間
  A_pre <- IN - rho * W - rho * w %*% t(alpha)
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
    Y00_all[t]   <- sum(alpha * yc0)
  }
