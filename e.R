rm(list = ls())
library(dplyr)
library(tidyr)
library(Matrix)

load("./data/california_smoking.rda")
california_smoking

Rcpp::sourceCpp("./R/20_mcmc_.cpp")
Rcpp::sourceCpp("./R/40_geweke.cpp")
source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover_.R")
source("R/03_utils_plot.R")
source("R/40_geweke.R")
source("R/04_utils_diagnostics.R")

panel_df <- california_smoking$panel
w_vec <- california_smoking$w
W_mat <- california_smoking$W
panel_df <- panel_df %>%
  mutate(
    treatment = ifelse((state == "California") & (year >= 1988), 1, 0)
  )

# 隣接情報は数値だけ渡すようにする
w <- as.matrix(w_vec[, 2])
W <- as.matrix(W_mat[, -1])

fit <- sc_spillover(
  data = panel_df,
  treated_unit = "California",
  w = w,
  W = W,
  treatment_dummy = "treatment",
  y = "cigsale",
  X = c("retprice"),
  p_factors = 1,
  M = 10000,
  burn = 5000,
  seed = 20251022,
  step_rho = 0.03,
  unit_col = "state",
  time_col = "year"
)

# 推定済みオブジェクト: fit
plot(fit, time_col = "year") # treated effect + 95%CI
plot(fit, type = "effect", time_col = "year")
plot(fit, type = "spill_top", top_n = 8, time_col = "year")
plot(fit, type = "weights") # 合成ウェイトの棒グラフ
plot(fit, type = "rho") # rho の事後分布
plot(fit, type = "beta") # beta の事後要約（あれば）
plot(fit, type = "trace") # rho のトレース














source("R/30_simulation_.R")

## 1) 真値の設定

make_rook_W <- function(nrow, ncol, normalize = TRUE) {
  stopifnot(nrow >= 1L, ncol >= 1L)
  idx <- matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE)
  N <- nrow * ncol
  W <- matrix(0, N, N)

  for (r in 1:nrow) {
    for (c in 1:ncol) {
      i <- idx[r, c]
      if (r > 1) {
        W[i, idx[r - 1, c]] <- 1
      }
      if (r < nrow) {
        W[i, idx[r + 1, c]] <- 1
      }
      if (c > 1) {
        W[i, idx[r, c - 1]] <- 1
      }
      if (c < ncol) W[i, idx[r, c + 1]] <- 1
    }
  }

  if (normalize) {
    rs <- rowSums(W)
    W <- W / pmax(rs, 1) # 0割回避
  }
  W
}

N <- 16
W <- make_rook_W(sqrt(N), sqrt(N))
W
w <- rep(0, N)
w[1:4] <- 1
alpha_true <- rep(0, N)
alpha_true[1] <- 0.5
alpha_true[2] <- -0.2
alpha_true[3:4] <- 0.4
alpha_true[5:10] <- 0.1 / 6

## 2) DGP で 1 セット生成
dgp_args <- list(
  T0 = 30,
  T1 = 20,
  N = N,
  W = W,
  w = w,
  rho = 0.0,
  sigma2 = 1.0,
  alpha = alpha_true,
  K = 1,
  beta = c(1.0)
)


# 1回実行（SCMとBayesian SCMの効果推定）
out1 <- run_one_sim(dgp = NULL, dgp_args = dgp_args, M = 5000, burn = 2000)
out1$metrics

# 多回
many <- run_many_sim(100, dgp_args, M = 5000, burn = 2000, )
summarize_many(many)


run_scenario_many <- function(
  n_sims,
  dgp_args,
  seeds = NULL,
  M = 5000,
  burn = 2000,
  step_rho = 0.02
) {
  res <- run_many_sim(
    n_sims = n_sims,
    dgp_args = dgp_args,
    seeds = seeds,
    M = M,
    burn = burn,
    step_rho = step_rho
  )
  list(
    summary = summarize_many(res),
    raw = res
  )
}

# =========================================================
# 引数:
#   Ns, T0s, rhos: ベクトル
#   T1: ポスト長（共通）
#   sims_per: 各シナリオの繰り返し回数
#   K, beta, sigma2: DGPのその他ハイパラ
#   treated_idx: w で 1 にするインデックス（例: 1:4）
#   alpha_fn: 真の α を返す関数。既定は mc_default_alpha
#   M, burn, step_rho: MCMC設定
# 返り値:
#   list(
#     summary = シナリオ別×手法別の要約テーブル結合,
#     details = 各シナリオの raw 結果のリスト
#   )
# =========================================================
mc_grid_study <- function(
  Ns,
  T0s,
  rhos,
  T1,
  sims_per = 100,
  K = 1,
  beta = c(1.0),
  sigma2 = 1.0,
  treated_idx = 1:4,
  M = 5000,
  burn = 2000,
  step_rho = 0.02,
  seeds = NULL
) {
  # グリッド展開
  grid <- expand.grid(N = Ns, T0 = T0s, rho = rhos, stringsAsFactors = FALSE)

  # 乱数種
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, nrow(grid) * sims_per)
  }
  # シナリオごとの種をスライス
  get_seeds_i <- function(i) {
    idx <- ((i - 1L) * sims_per + 1L):(i * sims_per)
    seeds[idx]
  }

  # 実行
  details <- vector("list", nrow(grid))
  summaries <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    N <- grid$N[i]
    T0 <- grid$T0[i]
    rho <- grid$rho[i]

    # N は正方格子前提（論文の体裁を踏襲）。
    m <- round(sqrt(N))
    if (m * m != N) {
      stop(sprintf(
        "N=%d は完全平方数ではありません（%dx%dにできません）。",
        N,
        m,
        m
      ))
    }

    # W と w
    W <- make_rook_W(m, m, normalize = TRUE)
    w <- numeric(N)
    w[intersect(treated_idx, seq_len(N))] <- 1

    # 真の α
    alpha_true <- rep(0, N)
    alpha_true[1] <- 0.5
    alpha_true[2] <- -0.2
    alpha_true[3:4] <- 0.4
    alpha_true[5:10] <- 0.1 / 6

    # DGP 引数
    dgp_args <- list(
      T0 = T0,
      T1 = T1,
      N = N,
      W = W,
      w = w,
      rho = rho,
      sigma2 = sigma2,
      alpha = alpha_true,
      K = K,
      beta = beta
    )

    # 多回実行
    out <- run_scenario_many(
      n_sims = sims_per,
      dgp_args = dgp_args,
      seeds = get_seeds_i(i),
      M = M,
      burn = burn,
      step_rho = step_rho
    )

    # シナリオ情報を列として付与
    smry <- out$summary
    smry$N <- N
    smry$T0 <- T0
    smry$T1 <- T1
    smry$rho <- rho
    # 列順を調整（読みやすさ）
    smry <- smry[, c(
      "N",
      "T0",
      "T1",
      "rho",
      "method",
      "bias_point",
      "rmse_point",
      "cover95_point"
    )]

    details[[i]] <- out$raw
    summaries[[i]] <- smry
  }

  list(
    summary = do.call(rbind, summaries),
    details = details,
    grid = grid
  )
}

study <- mc_grid_study(
  # Ns = c(16, 36),
  Ns = c(16),
  # T0s = c(30, 60),
  T0s = c(30),
  rhos = c(0.3),
  T1 = 20,
  sims_per = 1000,
  K = 1,
  beta = c(1.0),
  sigma2 = 1.0,
  treated_idx = 1:4,
  M = 5000,
  burn = 2500,
  step_rho = 1
)


a <- study$details[[1]]
a[[1]]

rho_draws <- a[[1]]$draws$scspill_step2$rho
alpha_draws <- a[[1]]$draws$scspill_step2$alpha_hat
plot(rho_draws, type='l')
