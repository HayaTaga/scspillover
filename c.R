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
# --- グリッドと重み ---
N <- 8
T0 <- 20
K <- 2
p <- 0
W <- matrix(0, N, N)
for (i in 1:N) {
  W[i, if (i < N) i + 1 else 1] <- 1
} # 環状1近傍
W <- W / rowSums(W)
w <- numeric(N)
w[2] <- 1 # 弱いスピル先

dims <- list(T0 = T0, N = N, K = K, p = p)

# --- 真値 ---
theta_true <- make_theta_true(
  rho = 0.30,
  sigma2 = 1.0,
  beta = c(0.5, -0.3),
  N = N,
  p = p
)

# --- 1 回走らせてチェーン診断 ---
set.seed(1)
one <- run_one_sim(
  dims = dims,
  W = W,
  w = w,
  theta_true = theta_true,
  iteration = 3000,
  burn = 1500,
  step_rho = 0.03
)

one$diag$rho # mean, sd, q025, q975, ESS, MCSE, GewekeZ, ACF1, ACF5
one$diag$sigma2
one$diag$acc_rate # RW-MH 受容率（rho）

# --- 繰り返し（例：R=50）で精度（バイアス/RMSE/被覆） ---
set.seed(42)
R <- 1000
many <- lapply(seq_len(R), function(r) {
  run_one_sim(
    dims,
    W,
    w,
    theta_true,
    iteration = 5000,
    burn = 1000,
    step_rho = 0.03,
    seed = 100 + r
  )
})

summarize_many(many)


## 1) 真値の設定
set.seed(123)
N <- 8
T0 <- 12
T1 <- 8
W <- matrix(0, N, N)
W[cbind(seq_len(N), (seq_len(N) %% N) + 1)] <- 1
W <- W / rowSums(W)
w <- rep(0, N)
w[2] <- 1

alpha_true <- c(0.40, 0.30, rep(0, N - 2))
names(alpha_true) <- paste0("u", seq_len(N))
rho_true <- 0.30
sigma2_true <- 1.00
K <- 2
beta_true <- c(0.7, -0.5)
names(beta_true) <- paste0("beta_", seq_len(K))

## 2) DGP で 1 セット生成
dgp <- scspill_sim_dgp(
  T0 = T0,
  T1 = T1,
  N = N,
  W = W,
  w = w,
  rho = rho_true,
  sigma2 = sigma2_true,
  alpha = alpha_true,
  K = K,
  beta = beta_true,
  s2y0 = 0.10,
  tau_post = rep(0, T1), # 介入なし
  seed = 1
)


# 1回実行（SCMとBayesian SCMの効果推定）
out1 <- run_one_sim(
  dgp,
  M_alpha = 2000,
  burn_alpha = 1000,
  M_sar = 4000,
  burn_sar = 2000
)
out1$scm$ATE
out1$bscm$ATE
out1$prop$ATE

# 多回
many <- replicate(
  200,
  run_one_sim(
    dgp,
    M_alpha = 2000,
    burn_alpha = 1000,
    M_sar = 4000,
    burn_sar = 2000
  ),
  simplify = FALSE
)
sum_te <- summarize_many(many)
sum_te$summary_methods
