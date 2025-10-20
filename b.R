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
Rcpp::sourceCpp("./R/20_mcmc.cpp")


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
  M = 10000,
  burn = 2000,
  seed = 1,
  unit_col = "state",
  time_col = "year"
)

# 推定済みオブジェクト: fit
plot(fit) # treated effect + 95%CI
plot(fit, type = "spill_top") # スピルオーバー上位ユニットの推移
plot(fit, type = "weights") # 合成ウェイトの棒グラフ
plot(fit, type = "rho") # rho の事後分布
plot(fit, type = "beta") # beta の事後要約（あれば）
plot(fit, type = "trace") # rho のトレース

pp_check(fit, what = "treat_last_dist")
pp_check(fit, what = "rho_trace")

diagnostics(fit, what = "trace")

# data <- panel_df
# w <- as.matrix(w_vec[, 2])
# W <- as.matrix(W_mat[, -1])
# unit_col <- "state"
# time_col <- "year"
# y <- "cigsale"
# treated_unit <- "California"
# X <- c("retprice")
# treatment_dummy <- "treatment"
# p_factors <- 1
# M <- 1000
# burn <- 500
# seed <- 1
# verbose <- TRUE

# data <- data %>%
#   mutate(
#     treatment = ifelse((state == treated_unit) & (year >= 1988), 1, 0)
#   )

# treatment_dummy <- "treatment"

# times <- sort(unique(data[[time_col]]))
# treat_series <- data[
#   data[[unit_col]] == treated_unit,
#   c(time_col, treatment_dummy)
# ]
# t_start <- min(treat_series[[time_col]][
#   treat_series[[treatment_dummy]] == 1
# ])
# T0 <- sum(times < t_start)

# # 前処理（列名を指定）
# prep <- scspill_prep_X(
#   data,
#   treated_unit = treated_unit,
#   T0 = T0,
#   X = X,
#   y_col = y,
#   unit_col = unit_col,
#   time_col = time_col
# )

# Y0_pre <- prep$Y0_pre
# Yc_pre <- prep$Yc_pre
# Y0_post <- prep$Y0_post
# Yc_post <- prep$Yc_post
# Xc_pre <- prep$Xc_pre
# N <- ncol(Yc_pre)
# T1 <- nrow(Yc_post)

# # (i) α推定（horseshoe）
# hs <- hs_alpha_gibbs(
#   y = Y0_pre,
#   X = Yc_pre,
#   M = M,
#   burn = burn,
#   verbose = verbose
# )
# alpha_draws <- hs$alpha
# alpha_hat <- colMeans(alpha_draws)

# sar <- sar_gibbs_sampler(
#   Y0_pre = Y0_pre,
#   Yc_pre = Yc_pre,
#   Xc_pre = Xc_pre,
#   W = W,
#   w = w,
#   M = M,
#   burn = burn,
#   verbose = verbose,
#   p_factors = p_factors
# )
# rho_draws <- sar$rho
# rho_hat <- mean(rho_draws)

# cred <- 0.95

# T1 <- length(Y0_post)
# N <- ncol(Yc_post)
# M <- nrow(alpha_draws)
# Aeff <- array(NA_real_, dim = c(T1, M)) # treatment effects per draw
# Spill <- array(NA_real_, dim = c(T1, N, M)) # spillovers per draw

# IN <- diag(N)

# m <- 1
# a <- as.matrix(alpha_draws[m, ])
# r <- rho_draws[m]
# Ainv <- solve(IN - r * (w %*% t(a) + W))
# B <- (IN - r * W)
# for (t in 1:T1) {
#   yc <- Yc_post[t, ]
#   y0 <- Y0_post[t]
#   tmp <- Ainv %*% (B %*% yc - r * w * y0)
#   # (5) treatment
#   Aeff[t, m] <- y0 - as.numeric(crossprod(a, tmp))
#   # (6) spillovers
#   Spill[t, , m] <- yc - as.vector(tmp)
# }
# M <- 2000
# burn <- 1000
# verbose <- TRUE
# p_factors <- 1
# phi_gamma <- 0.7
# sig2_g_init <- 1.0
# sig2_e_init <- 1.0
# step_rho <- 0.05

# T0 <- nrow(Yc_pre)
# N <- ncol(Yc_pre)
# K <- if (is.null(Xc_pre)) 0 else dim(Xc_pre)[3]
# p <- max(0, p_factors)

# rho <- 0
# sig2_e <- sig2_e_init
# sig2_g <- sig2_g_init
# Eta <- if (p > 0) matrix(0, N, p) else matrix(0, N, 0)
# gamma <- if (p > 0) matrix(0, p, T0) else matrix(0, 0, T0)
# beta_state <- if (K > 0) hs_beta_init(K) else NULL

# Mtot <- M + burn
# acc <- 0L
# rho_draws <- numeric(M)
# beta_draws <- if (K > 0) matrix(NA_real_, M, K) else NULL
# eta_draws <- if (p > 0) array(NA_real_, c(N, p, M)) else NULL
# gam_draws <- if (p > 0) array(NA_real_, c(p, T0, M)) else NULL
# sig2e_draws <- numeric(M)

# pb <- if (verbose) {
#   progress::progress_bar$new(total = Mtot, clear = FALSE)
# } else {
#   NULL
# }

# it <- 1

# if (verbose) {
#   pb$tick()
# }
# A <- diag(N) - rho * W
# R_list <- vector("list", T0)

# library(Matrix)

# t <- 1
# yc <- Matrix(as.numeric(Yc_pre[t, ]))
# Xt <- Matrix(
#   matrix(as.numeric(Xc_pre[t, , , drop = FALSE]), ncol = K),
#   sparse = TRUE
# )
# R_list[[t]] <- as.vector(A %*% yc - rho * w * Y0_pre[t] - xterm)

# for (t in 1:T0) {
#   yc <- as.matrix(Yc_pre[t, ])
#   Xt <- Matrix(
#     matrix(as.numeric(Xc_pre[t, , , drop = FALSE]), ncol = K),
#     sparse = TRUE
#   )
#   b <- Matrix(matrix(beta_state$beta, ncol = 1), sparse = TRUE)
#   xterm <- if (K > 0) Xt %*% b else rep(0, N)
#   R_list[[t]] <- as.vector(A %*% yc - rho * w * Y0_pre[t] - xterm)
# }

# gamma <- ffbs_ar1(
#   R_list = R_list,
#   Eta = Eta,
#   phi = phi_gamma,
#   sig2_g = sig2_g,
#   sig2_e = sig2_e
# )

# R_list_eta <- vector("list", T0)
# for (t in 1:T0) {
#   R_list_eta[[t]] <- R_list[[t]] + as.vector(Eta %*% gamma[, t])
# }
# Sigma_eta <- diag(10, p)
# Eta <- update_eta(
#   R_list = R_list_eta,
#   gamma = gamma,
#   Sigma_eta = Sigma_eta,
#   sig2_e = sig2_e
# )

# y_stack <- numeric(T0 * N)
# X_stack <- matrix(0, T0 * N, K)
# for (t in 1:T0) {
#   idx <- ((t - 1) * N + 1):(t * N)
#   y_stack[idx] <- if (p > 0) {
#     R_list[[t]] + as.vector(Eta %*% gamma[, t])
#   } else {
#     R_list[[t]]
#   }
#   X_stack[idx, ] <- Xc_pre[t, , ]
# }
# beta_state <- hs_beta_update(y = y_stack, X = X_stack, state = beta_state)

# mh <- rho_mh_step(
#   rho_cur = rho,
#   Yc_pre = Yc_pre,
#   Y0_pre = Y0_pre,
#   Xc_pre = Xc_pre,
#   beta = if (K > 0) beta_state$beta else rep(0, K),
#   W = W,
#   w = w,
#   Eta = if (p > 0) Eta else matrix(0, N, 0),
#   gamma = if (p > 0) gamma else matrix(0, 0, T0),
#   sig2_e = sig2_e,
#   step = step_rho
# )
# rho <- mh$rho
# if (mh$accept && it > burn) {
#   acc <- acc + 1L
# }

# ss <- 0
# for (t in 1:T0) {
#   u <- if (p > 0) {
#     R_list[[t]] - as.vector(Eta %*% gamma[, t])
#   } else {
#     R_list[[t]]
#   }
#   ss <- ss + sum(u^2)
# }
# shape <- (N * T0) / 2 + 0.1
# rate <- ss / 2 + 0.1
# sig2_e <- 1 / rgamma(1, shape = shape, rate = rate)
# if (it > burn) {
#   k <- it - burn
#   rho_draws[k] <- rho
#   if (K > 0) {
#     beta_draws[k, ] <- beta_state$beta
#   }
#   if (p > 0) {
#     eta_draws[,, k] <- Eta
#     gam_draws[,, k] <- gamma
#   }
#   sig2e_draws[k] <- sig2_e
# }
# list(
#   rho = rho_draws,
#   beta = beta_draws,
#   Eta = eta_draws,
#   gamma = gam_draws,
#   sigma2_e = sig2e_draws,
#   acc_rate = acc / max(1, M)
# )
