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
source("R/04_utils_diagnostics.R")
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
  M = 50000,
  burn = 25000,
  seed = 20251022,
  step_rho = 0.1,
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


load("./data/sudan_secession_latest.rda")

panel_df <- sudan_secession_latest$panel
w_vec <- sudan_secession_latest$w
W_mat <- sudan_secession_latest$W
panel_df <- panel_df %>%
  mutate(
    treatment = ifelse((country == "Sudan") & (year >= 2011), 1, 0)
  )

w <- as.matrix(w_vec[, 2])
W <- as.matrix(W_mat[, -1])

fit <- sc_spillover(
  data = panel_df,
  treated_unit = "Sudan",
  w = w,
  W = W,
  treatment_dummy = "treatment",
  y = colnames(panel_df)[11],
  X = colnames(panel_df)[4:10],
  p_factors = 1,
  M = 5000,
  burn = 2500,
  seed = 20251022,
  step_rho = 0.1,
  unit_col = "country",
  time_col = "year"
)

plot(fit, time_col = "year")
plot(fit, type = "spill_top", top_n = 8, time_col = "year")
p <- diagnostics.scspill(
  fit,
  which_alpha = NULL,
  top_n_alpha = 6,
  which_beta = NULL,
  top_n_beta = 6
)
print(p)
