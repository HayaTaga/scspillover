rm(list = ls())
library(dplyr)
library(tidyr)
library(Matrix)

load("./data/california_smoking.rda")
california_smoking


source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover.R")
source("R/03_utils_plot.R")
source("R/40_geweke.R")
Rcpp::sourceCpp("./R/20_mcmc.cpp")
Rcpp::sourceCpp("./R/40_geweke.cpp")


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
  y = "cigsale", # 例: アウトカム列が "smoking_rate" の場合
  X = c("retprice"), # 共変量（列名ベクトル）
  p_factors = 1, # Appendixの因子レイヤを1つ使用
  M = 2000,
  burn = 500,
  seed = 20251022,
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

# トレース描画（従来通り）
p <- diagnostics.scspill(
  fit,
  which_alpha = NULL,
  top_n_alpha = 6,
  which_beta = NULL,
  top_n_beta = 6
)
print(p)

# 診断表（ESS, MCSE, ACT, split-Rhat, Geweke Z など）
diag_tab <- attr(p, "summary")
diag_tab[order(diag_tab$ess, decreasing = TRUE), ]


set.seed(1)
N <- 16L
K <- 0L
p <- 0L
T0 <- 12L

# Ring adjacency, row-stochastic W
W <- matrix(0, N, N)
diag(W) <- 0
for (i in 1:N) {
  j <- if (i < N) i + 1 else 1
  W[i, j] <- 1
}
W <- W / rowSums(W)

# w: unit vector selecting the 2nd control
w <- rep(0, N)
w[2] <- 1

# Pre-treatment treated outcome Y0_pre: AR(1) with phi=0.6
Y0_pre <- numeric(T0)
eps <- rnorm(T0, 0, 1)
phi <- 0.6
Y0_pre[1] <- eps[1]
for (t in 2:T0) {
  Y0_pre[t] <- phi * Y0_pre[t - 1] + eps[t]
}

# No regressors in this example
Xc <- NULL

out <- geweke_jdt(
  Y0_pre,
  W,
  w,
  Xc = Xc,
  p = p,
  M1 = 10000,
  M2 = 10000, # 推奨規模は適宜調整
  a0 = 1,
  b0 = 1, # sigma^2 の IG 事前
  step_rho = 0.01, # 既存 RW-MH の提案幅
  g_fn = robust_g_fn, # 必要に応じて差し替え可
  verbose = TRUE
)

out$summary # g ごとの平均・SE・Z・p 値

# baseline estimation
library(spdep)
library(splm)

lw <- mat2listw(W, style = "W")
lw

df_long <- panel_df %>% filter(state_id != 0)

sar_pool <- spml(
  cigsale ~ retprice,
  data = df_long,
  index = c("state_id", "year"),
  listw = lw,
  model = "pooling",
  lag = TRUE,
  spatial.error = "none",
  method = "eigne"
)
summary(sar_pool)
