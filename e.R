rm(list = ls())
library(dplyr)
library(tidyr)
library(Matrix)

load("./data/california_smoking.rda")
california_smoking

Rcpp::sourceCpp("./R/20_mcmc.cpp")
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
  step_rho = 0.5,
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

make_rook_W <- function(nrow, ncol, normalize = FALSE) {
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
  rho = 0.3,
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
    W <- make_rook_W(m, m)
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
    process_time <- system.time({
      out <- run_scenario_many(
        n_sims = sims_per,
        dgp_args = dgp_args,
        seeds = get_seeds_i(i),
        M = M,
        burn = burn,
        step_rho = step_rho
      )
    })
    print(process_time)

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

# source("R/30_simulation_.R")
study <- mc_grid_study(
  # Ns = c(16, 36),
  Ns = c(16, 36),
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
  burn = 1000,
  step_rho = 0.15
)
study$summary

a <- study$details[[1]]
a[[1]]

rho_draws <- a[[1]]$draws$scspill_step2$rho
alpha_draws <- a[[1]]$draws$scspill_step2$alpha_hat
plot(rho_draws, type = 'l')

mean(rho_draws)
quantile(rho_draws, probs = c(0.025, 0.5, 0.975))
length(unique(rho_draws)) / length(rho_draws)

mean(a[[1]]$draws$ate$scspill)
quantile(a[[1]]$draws$ate$scspill)

a[[1]]$draws$ate$scspill

a[[1]]$truth$tau_post


rm(list = ls())
Rcpp::sourceCpp("R/20_mcmc.cpp")
Rcpp::sourceCpp("R/40_geweke_latest.cpp")
source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover.R")
source("R/03_utils_plot.R")
source("R/40_geweke_latest.R")
source("R/04_utils_diagnostics.R")

load("./data/california_smoking.rda")
california_smoking

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
  M = 5000,
  burn = 1000,
  seed = 20251022,
  step_rho = 0.05,
  unit_col = "state",
  time_col = "year"
)
plot(fit, type = "trace")

alpha_hat <- colMeans(fit$alpha_draws)
sds_Yc <- apply(fit$inputs$Yc_pre, 2, stats::sd)
sds_Yc[!is.finite(sds_Yc) | sds_Yc < 1e-8] <- 1.0
alpha_hat_scaled <- alpha_hat / sds_Yc

# Xを入れるには N x T x K の行列で渡す必要あり

out_jdt <- geweke_jdt_full(
  Y0_pre = fit$inputs$Y0_pre,
  Yc_pre_like_dims = c(length(fit$inputs$Y0_pre), ncol(fit$inputs$Yc_pre)),
  W_raw = as.matrix(fit$inputs$W),
  w_raw = as.matrix(fit$inputs$w),
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = NULL,
  p = 1,
  M1 = 100000,
  M2 = 100000,
  burn_in = 50000,
  a0 = 1.0,
  b0 = 1.0,
  step_rho = 0.05,
  g_fn = default_g_fn,
  verbose = TRUE
)
out_jdt$summary


## =========================================================
##  Geweke JDT 実装の健全性検証（自己完結スクリプト）
##  前提：40_geweke_full.cpp / 40_geweke_full.R がロード済み
## =========================================================

set.seed(20251030) # 再現性

## -------------------------------
## 1) 小規模設定の生成ルーチン
## -------------------------------
make_small_setup <- function(T0 = 8L, N = 6L, K = 2L, p = 1L) {
  # 線形鎖の隣接行列
  W_raw <- matrix(0, N, N)
  for (i in 1:N) {
    if (i > 1) {
      W_raw[i, i - 1] <- 1
    }
    if (i < N) W_raw[i, i + 1] <- 1
  }
  # treated に接続する重み
  w_raw <- rep(0, N)
  w_raw[1] <- 1

  # 行同時正規化（JDT 内でも行われるが、ここでは安定性のため事前確認）
  normed <- normalize_joint_wW(W_raw, w_raw)
  W_use <- normed$W_new
  w_use <- normed$w_new

  # alpha_hat_scaled：任意のスパースでないベクトル
  alpha_hat_scaled <- drop(scale(runif(N))) # 長さ N

  # X（あってもなくても検証できる）。T0 x N x K
  if (K > 0) {
    Xc_pre <- array(0, c(T0, N, K))
    for (t in 1:T0) {
      for (i in 1:N) {
        for (k in 1:K) {
          Xc_pre[t, i, k] <- rnorm(1)
        }
      }
    }
  } else {
    Xc_pre <- NULL
  }

  # g() の一部で使う pre の treated 系列
  Y0_pre <- rnorm(T0, 0, 1)

  list(
    T0 = T0,
    N = N,
    K = K,
    p = p,
    W_raw = W_raw,
    w_raw = w_raw,
    alpha_hat_scaled = alpha_hat_scaled,
    Xc_pre = Xc_pre,
    Y0_pre = Y0_pre
  )
}

## ----------------------------------------
## 2) 1 回の JDT を実行（整合ケース）
## ----------------------------------------
run_one_jdt <- function(
  setup,
  M1 = 6000L,
  M2 = 6000L,
  burn_in = 2000L,
  a0 = 1.0,
  b0 = 1.0,
  step_rho = 0.05,
  verbose = FALSE
) {
  out <- geweke_jdt_full(
    Y0_pre = setup$Y0_pre,
    Yc_dims = c(setup$T0, setup$N),
    W_raw = setup$W_raw,
    w_raw = setup$w_raw,
    alpha_hat_scaled = setup$alpha_hat_scaled,
    Xc_pre = setup$Xc_pre,
    p = setup$p,
    M1 = M1,
    M2 = M2,
    burn_in = burn_in,
    a0 = a0,
    b0 = b0,
    step_rho = step_rho,
    g_fn = default_g_fn,
    verbose = verbose
  )
  out
}

## ------------------------------------------------------
## 3) 帰無（整合）・ストレス（混合悪化）の 2 条件を比較
##    - 帰無： step_rho を通常（例 0.05）
##    - ストレス： step_rho を極小（例 0.001）にして SC の混合を悪化
## ------------------------------------------------------
evaluate_once <- function(
  setup,
  M1 = 4000L,
  M2 = 4000L,
  burn_in = 1500L,
  step_ok = 0.05,
  step_bad = 0.001,
  verbose = FALSE
) {
  out_ok <- run_one_jdt(
    setup,
    M1,
    M2,
    burn_in,
    step_rho = step_ok,
    verbose = verbose
  )
  out_bad <- run_one_jdt(
    setup,
    M1,
    M2,
    burn_in,
    step_rho = step_bad,
    verbose = verbose
  )

  ok_sig <- with(out_ok$summary, setNames(pval < 0.05, g))
  bad_sig <- with(out_bad$summary, setNames(pval < 0.05, g))

  list(
    ok = out_ok,
    bad = out_bad,
    sig_flags = data.frame(
      g = out_ok$summary$g,
      reject_ok = ok_sig,
      reject_bad = bad_sig,
      stringsAsFactors = FALSE
    )
  )
}

## ------------------------------------------------------
## 4) 多回試行で経験的サイズ（帰無）と偽陽性増大（ストレス）を観察
## ------------------------------------------------------
replicate_assessment <- function(
  B = 20L,
  T0 = 8L,
  N = 6L,
  K = 2L,
  p = 1L,
  M1 = 3000L,
  M2 = 3000L,
  burn_in = 1000L,
  step_ok = 0.05,
  step_bad = 0.001,
  verbose_each = FALSE
) {
  rej_ok <- NULL
  rej_bad <- NULL

  for (b in 1:B) {
    set.seed(20251030 + b)
    setup <- make_small_setup(T0 = T0, N = N, K = K, p = p)
    res <- evaluate_once(
      setup,
      M1 = M1,
      M2 = M2,
      burn_in = burn_in,
      step_ok = step_ok,
      step_bad = step_bad,
      verbose = verbose_each
    )
    if (is.null(rej_ok)) {
      rej_ok <- as.data.frame(t(res$sig_flags$reject_ok))
      rej_bad <- as.data.frame(t(res$sig_flags$reject_bad))
      colnames(rej_ok) <- res$sig_flags$g
      colnames(rej_bad) <- res$sig_flags$g
    } else {
      rej_ok <- rbind(rej_ok, as.data.frame(t(res$sig_flags$reject_ok)))
      rej_bad <- rbind(rej_bad, as.data.frame(t(res$sig_flags$reject_bad)))
    }
    if (verbose_each) message("rep ", b, "/", B)
  }

  size_ok <- colMeans(rej_ok) # 帰無：名目 5% 付近が望ましい
  size_bad <- colMeans(rej_bad) # ストレス：偽陽性が増えるのが自然

  list(
    size_ok = size_ok,
    size_bad = size_bad,
    rej_ok = rej_ok,
    rej_bad = rej_bad
  )
}

## ================== 実行例 ==================

# まずは 1 回の結果（帰無 vs ストレス）を確認
setup <- make_small_setup(T0 = 8L, N = 6L, K = 2L, p = 1L)

cat("== Smoke test: 1 run (OK) ==\n")
out_ok <- run_one_jdt(
  setup,
  M1 = 4000L,
  M2 = 4000L,
  burn_in = 1500L,
  step_rho = 0.05,
  verbose = TRUE
)
print(out_ok$summary)

cat("\n== Smoke test: 1 run (BAD mixing) ==\n")
out_bad <- run_one_jdt(
  setup,
  M1 = 4000L,
  M2 = 4000L,
  burn_in = 1500L,
  step_rho = 0.001,
  verbose = TRUE
)
print(out_bad$summary)

# 簡易の多回（デフォルト 20 回）で経験的サイズを見る
cat("\n== Replicated assessment (B=20) ==\n")
assess <- replicate_assessment(
  B = 20L,
  T0 = 8L,
  N = 6L,
  K = 2L,
  p = 1L,
  M1 = 3000L,
  M2 = 3000L,
  burn_in = 1000L,
  step_ok = 0.05,
  step_bad = 0.001,
  verbose_each = FALSE
)

cat("\nEmpirical size at 5% (correct implementation, OK step_rho):\n")
print(round(assess$size_ok, 3))

cat("\nEmpirical rejection under stressed mixing (very small step_rho):\n")
print(round(assess$size_bad, 3))

## 期待される振る舞い：
## - size_ok は 0.05 付近（統計量ごとに多少のばらつきは許容）
## - size_bad は 0.05 より大きくなりやすい（SC 側の混合劣化で偽陽性増大）
